!-----------------------------------------------------------
!! Module that will hopefully make tabulation easy
!!
!! Currently does not know anything about nf thresholds.
!!
!! But some form of support for these might be useful... Question is:
!! how can it be arranged in a fairly clean manner, without breaking
!! things that already use the system (e.g. caesar resum)?
!!
!! One possibility is to use an running_coupling as a source of information.
!! One can then either take a copy of the running_coupling and use its routines
!! or else just extract the relevant useful information
!!
!! Aim is to be able to write -- i.e. convolutions with full nf info
!!
!!   Ptab%tab = dh%P_LO(iloop,tab%nf_int(:)) .conv. tab%tab
!!
!! This probably requires some modifications also to the convolution 
!! routines (which do not currently have this overloading for .conv.).
!!
!! [NB at some point maybe include nf association in some simpler manner?]
module pdf_tabulate
  use types; use consts_dp
  use convolution; use dglap_objects
  use pdf_representation; use pdf_general
  use interpolation
  use warnings_and_errors
  !-- needed only for EvolvePdfTable [separate it out one day?]
  use qcd_coupling; use evolution; use dglap_holders
  implicit none
  private

  
  type pdfseginfo
     real(dp) :: lnlnQ_lo, lnlnQ_hi, dlnlnQ
     integer :: ilnlnQ_lo, ilnlnQ_hi
  end type pdfseginfo
  public :: pdfseginfo

  type pdf_table
     ! basic elements of a pdf_table, common regardless of whether we
     ! additionally have the nf segments...
     type(grid_def) :: grid
     real(dp) :: default_dlnlnQ
     real(dp) :: lnlnQ_min, lnlnQ_max, lambda_eff
     real(dp), pointer :: tab(:,:,:)
     real(dp), pointer :: lnlnQ_vals(:), Q_vals(:)
     integer :: nQ, lnlnQ_order
     logical :: freeze_at_Qmin
     ! this is useful only in absence of nf info.
     real(dp) :: dlnlnQ
     !
     ! Stuff to do with variable nf and alpha_s; not always available.
     ! In cases with variable nf, the table will be broken into multiple
     ! segments, each one of which potentially different spacings.
     !
     logical :: nf_info_associated
     integer :: nflo, nfhi
     type(pdfseginfo), pointer :: seginfo(:)
     integer, pointer :: nf_int(:)
     real(dp), pointer :: as2pi(:)
     !
     ! Elements needed in case we want to do precalculation of
     ! of evolution. Not always available.
     type(evln_operator), pointer :: evops(:)
     integer                      :: StartScale_iQlo
     real(dp)                     :: StartScale
  end type pdf_table
  public :: pdf_table

  !-- for calculating ln ln Q/lambda_eff
  real(dp), parameter :: default_lambda_eff = 0.1_dp
  real(dp), parameter :: default_dlnlnQ = 0.07_dp
  real(dp), parameter :: warn_tolerance = 1e-3_dp
  ! used in various contexts for deciding when an interval is
  ! sufficiently small that it can be ignored...
  real(dp), parameter :: min_dlnlnQ_singleQ = 1e-10_dp
  integer, parameter :: def_lnlnQ_order = 4
  !integer, parameter :: lnlnQ_order = 2

  interface AllocPdfTable
     module procedure pdftab_AllocTab_, pdftab_AllocTab_fromorig,&
          & pdftab_AllocTab_1d, pdftab_AllocTab_fromorig_1d
  end interface

  interface AddNfInfoToPdfTable
     module procedure AddNfInfoToPdfTable, pdftab_AssocNfInfo_1d
  end interface

  interface FillPdfTable_f90sub
     module procedure pdftab_InitTabSub_, pdftab_InitTabSub_iset
  end interface

  public :: FillPdfTable_LHAPDF

  interface EvolvePdfTable
     module procedure pdftab_InitTabEvolve_frompre, EvolvePdfTable
  end interface

  interface Delete
     module procedure pdftab_DelTab_0d, pdftab_DelTab_1d
  end interface

  public :: AllocPdfTable, FillPdfTable_f90sub
  public :: AddNfInfoToPdfTable
  public :: EvolvePdfTable, PreEvolvePdfTable
  public :: EvolvePdfTableGen
  public :: EvalPdfTable_yQ, EvalPdfTable_xQ, EvalPdfTable_Q
  public :: Delete

contains

  !---------------------------------------------------------
  !! Allocate a pdf_table, which covers a range Qmin to Qmax, using
  !! tabulation uniform in (ln ln Q/Lambda), where Lambda is currently
  !! a module parameter.
  !!
  !! If freeze_at_Qmin is present and .true., distributions are frozen
  !! at their Qmin value for Q < Qmin. Otherwise they are set to zero
  !! there.
  !!
  !! More sensible extrapolation beyond Q range offers scope for future 
  !! improvement here!
  !!
  subroutine pdftab_AllocTab_(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin)
    type(grid_def),    intent(in)  :: grid
    type(pdf_table),      intent(out) :: tab
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    integer,  intent(in), optional :: lnlnQ_order
    logical,  intent(in), optional :: freeze_at_Qmin
    !----------------------------------------------
    integer :: iQ
    tab%grid = grid
    tab%lambda_eff = min(half*Qmin, default_lambda_eff)
    tab%lnlnQ_min = lnln(tab,Qmin)
    tab%lnlnQ_max = lnln(tab,Qmax)

    tab%default_dlnlnQ = default_or_opt(default_dlnlnQ, dlnlnQ)
    tab%nQ = ceiling((tab%lnlnQ_max - tab%lnlnQ_min)/tab%default_dlnlnQ)
    tab%dlnlnQ = (tab%lnlnQ_max - tab%lnlnQ_min)/tab%nQ
    tab%freeze_at_Qmin = default_or_opt(.false.,freeze_at_Qmin)
    tab%lnlnQ_order = default_or_opt(def_lnlnQ_order,lnlnQ_order)

    !-- by default, no extra info is given
    tab%nf_info_associated = .false.
    nullify(tab%as2pi)
    nullify(tab%nf_int)
    nullify(tab%evops)

    !write(0,*) 'pdf_table info: Number of Q bins is ',tab%nQ
    call AllocPDF(grid,tab%tab,0,tab%nQ)
    allocate(tab%lnlnQ_vals(0:tab%nQ))
    allocate(tab%Q_vals(0:tab%nQ))
    do iQ = 0, tab%nQ
       tab%lnlnQ_vals(iQ) = tab%lnlnQ_min + iQ * tab%dlnlnQ
       tab%Q_vals(iQ) = invlnln(tab,tab%lnlnQ_vals(iQ))
    end do
    
  end subroutine pdftab_AllocTab_


  
  !---------------------------------------------------------
  !! 1d overloaded version of AllocPdfTable
  subroutine pdftab_AllocTab_1d(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin)
    type(grid_def), intent(in)     :: grid
    type(pdf_table),   intent(out)    :: tab(:)
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    integer,  intent(in), optional :: lnlnQ_order
    logical,  intent(in), optional :: freeze_at_Qmin
    integer :: i

    do i = 1, size(tab)
       call pdftab_AllocTab_(grid, tab(i), Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin)
    end do
  end subroutine pdftab_AllocTab_1d


  !---------------------------------------------------------------
  !! Associate a tab with the nf info from an alpha_s holder
  !! (i.e. read that info, and copy it into tab).
  !!
  !! This involves REALLOCATING the tab pointer so as to allow for
  !! extra information.
  !!
  !! The new tab separates stretches with different nf values. It also
  !! has extra info in the form of arrays nf_int and as2pi, which make
  !! it possible to write expressions such as
  !!
  !!   Dtab%tab = tab%as2pi(:) * (dh%P_LO(iloop,tab%nf_int(:)) .conv. tab%tab)
  !!
  !! This is one of the rare cases in which direct access to structure 
  !! components is still allowed... (Lack of finalisation makes it
  !! difficult to do otherwise).
  !!
  !! KNOWN LIMITATIONS: What happens with muM_mQ /= 1? All hell breaks
  !! loose?
  subroutine AddNfInfoToPdfTable(tab,coupling)
    type(pdf_table), intent(inout) :: tab
    type(running_coupling), intent(in) :: coupling
    !-----------------------------------
    integer  :: nflcl, iQ_prev, iQ
    real(dp) :: Qlo, Qhi, Qhi_test
    type(pdfseginfo), pointer :: seginfo

    
    !if (tab%nf_info_associated) call wae_error('AddNfInfoToPdfTable',&
    !     &'nf info already associated: delete it first')

    ! We will be reallocating everything here, so first clean up
    call Delete(tab)
    tab%dlnlnQ = zero ! will no longer be useful...

    tab%nflo = NfAtQ(coupling, invlnln(tab,tab%lnlnQ_min))
    tab%nfhi = NfAtQ(coupling, invlnln(tab,tab%lnlnQ_max))

    allocate(tab%seginfo(tab%nflo:tab%nfhi))

    ! figure out how we are going to bin things...
    iQ_prev = -1
    do nflcl = tab%nflo, tab%nfhi
       !write(0,*) 'nflcl',nflcl
       seginfo => tab%seginfo(nflcl)

       call QRangeAtNf(coupling, nflcl, Qlo, Qhi_test)
       call QRangeAtNf(coupling, nflcl, Qlo, Qhi, muM_mQ=one)
       ! if one weakens this restriction, then one should think about
       ! the consequences for the determination of alpha_s/2pi, here
       ! and elsewhere....
       if (Qhi_test /= Qhi) call wae_error('AddNfInfoToPdfTable',&
         &'it seems that coupling has muM_mQ /= one. Currently unsupported.',&
         &dbleval=Qhi_test/Qhi)

       ! Include min_dlnlnQ_singleQ to ensure we do not go EXACTLY to the
       ! mass threshold, where, in evolution one might run into problems.
       ! BUT, we will later have to worry about what to do when we
       ! are in between thresholds...
       seginfo%lnlnQ_lo = max(lnln(tab,Qlo)+min_dlnlnQ_singleQ, tab%lnlnQ_min)
       seginfo%lnlnQ_hi = min(lnln(tab,Qhi)-min_dlnlnQ_singleQ, tab%lnlnQ_max)

       seginfo%ilnlnQ_lo = iQ_prev + 1
       !write(0,*) 'ill_lo', seginfo%ilnlnQ_lo
       if ((seginfo%lnlnQ_hi - seginfo%lnlnQ_lo) < two*min_dlnlnQ_singleQ) then
          ! use just one point
          seginfo%ilnlnQ_hi = seginfo%ilnlnQ_lo
          seginfo%dlnlnQ = zero
          seginfo%lnlnQ_hi = seginfo%lnlnQ_lo
       else
          seginfo%ilnlnQ_hi = seginfo%ilnlnQ_lo + max(tab%lnlnQ_order,&
               & ceiling((seginfo%lnlnQ_hi - seginfo%lnlnQ_lo)/&
               & tab%default_dlnlnQ))
          seginfo%dlnlnQ = (seginfo%lnlnQ_hi - seginfo%lnlnQ_lo)/&
               & (seginfo%ilnlnQ_hi - seginfo%ilnlnQ_lo)
       end if
       !write(0,*) 'ill_hi', seginfo%ilnlnQ_hi, seginfo%dlnlnQ, &
       !     &invlnln(tab,seginfo%lnlnQ_lo),invlnln(tab,seginfo%lnlnQ_hi), tab%default_dlnlnQ
       iQ_prev = seginfo%ilnlnQ_hi
    end do
    
    ! this should not happen too often! But check it just in
    ! case...
    if (tab%seginfo(tab%nflo)%lnlnQ_lo /= tab%lnlnQ_min .or.&
         &tab%seginfo(tab%nfhi)%lnlnQ_hi /= tab%lnlnQ_max) &
         &call wae_error('AddNfInfoToPdfTable',&
         & 'mismatch in segment and global lnlnQ limits.',&
         & 'Could be due to coupling having more restricted range?')


    ! now reallocate things?
    tab%nQ = tab%seginfo(tab%nfhi)%ilnlnQ_hi
    call AllocPDF(tab%grid,tab%tab,0,tab%nQ)
    allocate(tab%lnlnQ_vals(0:tab%nQ))
    allocate(tab%Q_vals(0:tab%nQ))
    allocate(tab%nf_int(0:tab%nQ))
    allocate(tab%as2pi(0:tab%nQ))

    ! set up complementary info...
    do nflcl = tab%nflo, tab%nfhi
       !write(0,*) 'nflcl',nflcl
       seginfo => tab%seginfo(nflcl)
       do iQ = seginfo%ilnlnQ_lo, seginfo%ilnlnQ_hi
          tab%nf_int(iQ) = nflcl
          tab%lnlnQ_vals(iQ) = seginfo%lnlnQ_lo &
               & + (iQ-seginfo%ilnlnQ_lo)*seginfo%dlnlnQ
          tab%Q_vals(iQ)     = invlnln(tab,tab%lnlnQ_vals(iQ))
          tab%as2pi(iQ) = Value(coupling,invlnln(tab,tab%lnlnQ_vals(iQ)))/twopi
       end do
    end do
    
    ! REMEMBER TO COMPLETE FROM ORIG...
    tab%nf_info_associated = .true.
    !write(0,*) 'pdf_table info: Number of Q bins changed to ',tab%nQ
  end subroutine AddNfInfoToPdfTable


  !---------------------------------------------------------------
  !! 1d-overloaded verseion of AddNfInfoToPdfTable
  subroutine pdftab_AssocNfInfo_1d(tab,coupling)
    type(pdf_table), intent(inout) :: tab(:)
    type(running_coupling), intent(in) :: coupling
    !-----------------------------------
    integer :: i
    do i = 1, size(tab)
       call AddNfInfoToPdfTable(tab(i),coupling)
    end do
    
  end subroutine pdftab_AssocNfInfo_1d


  !-----------------------------------------------------------------------
  !! Allocate the memory for a new tab, using as a template an
  !! preexistent tab (origtab).
  !!
  !! Additionally, information concerning any varnf and alphas
  !! structure copied from origtab to tab. Actual PDF contents of the
  !! tab are not however copied.
  !! 
  subroutine pdftab_AllocTab_fromorig(tab, origtab)
    type(pdf_table), intent(out) :: tab
    type(pdf_table), intent(in)  :: origtab

    tab = origtab
    !-- this is the only thing that is not taken care of...
    call AllocPDF(tab%grid,tab%tab,0,tab%nQ)
    allocate(tab%lnlnQ_vals(0:tab%nQ))
    tab%lnlnQ_vals = origtab%lnlnQ_vals
    allocate(tab%Q_vals(0:tab%nQ))
    tab%Q_vals = origtab%Q_vals

    if (origtab%nf_info_associated) then
       allocate(tab%seginfo(tab%nflo:tab%nfhi))
       allocate(tab%nf_int(0:tab%nQ))
       allocate(tab%as2pi(0:tab%nQ))
       tab%seginfo = origtab%seginfo
       tab%nf_int = origtab%nf_int
       tab%as2pi = origtab%as2pi
    end if
  end subroutine pdftab_AllocTab_fromorig

  !---------------------------------------------------------
  !! 1d-overloaded version of pdftab_AllocTab_fromorig
  subroutine pdftab_AllocTab_fromorig_1d(tab, origtab)
    type(pdf_table), intent(out) :: tab(:)
    type(pdf_table), intent(in)  :: origtab
    integer :: i
    do i = 1, size(tab)
       call pdftab_AllocTab_fromorig(tab(i), origtab)
    end do
  end subroutine pdftab_AllocTab_fromorig_1d


  !------------------------------------------------------
  !! Initialise the tab with the results of a subroutine (see
  !! interface). Note that the subroutine takes y = ln1/x and Q as its
  !! arguments.
  subroutine pdftab_InitTabSub_(tab, sub)
    type(pdf_table), intent(inout) :: tab
    interface
       subroutine sub(y, Q, res)
         use types; implicit none
         real(dp), intent(in) :: y, Q
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-------------
    integer :: iQ
    real(dp) :: Q
    do iQ = 0, tab%nQ
       Q = invlnln(tab,tab%lnlnQ_min + iQ*tab%dlnlnQ)
       call InitPDFSub(tab%grid, tab%tab(:,:,iQ), sub, Q)
    end do
    
  end subroutine pdftab_InitTabSub_
  
  !------------------------------------------------------
  !! Initialise the tab with the results of a subroutine (see
  !! interface). In addition y = ln1/x and Q, the subroutine takes the
  !! argument iset, enabling one to initialise from a subroutine that
  !! provides several "PDF sets".
  !!
  subroutine pdftab_InitTabSub_iset(tab, sub, iset)
    type(pdf_table), intent(inout) :: tab
    integer,      intent(in)    :: iset
    interface
       subroutine sub(y, Q, iset, res)
         use types; implicit none
         real(dp), intent(in) :: y, Q
         integer, intent(in) :: iset
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-------------
    integer :: iQ
    real(dp) :: Q
    do iQ = 0, tab%nQ
       !Q = invlnln(tab,tab%lnlnQ_min + iQ*tab%dlnlnQ)
       Q = invlnln(tab,tab%lnlnQ_vals(iQ))
       call InitPDFSub(tab%grid, tab%tab(:,:,iQ), sub, Q, iset)
    end do
    
  end subroutine pdftab_InitTabSub_iset

  !------------------------------------------------------
  !! Initialise the tab with the results of a subroutine (see
  !! interface). Note that the subroutine takes y = ln1/x and Q as its
  !! arguments.
  subroutine FillPdfTable_LHAPDF(tab, LHAsub)
    type(pdf_table), intent(inout) :: tab
    interface
       subroutine LHAsub(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine LHAsub
    end interface
    !-------------
    integer :: iQ
    real(dp) :: Q
    do iQ = 0, tab%nQ
       Q = invlnln(tab,tab%lnlnQ_min + iQ*tab%dlnlnQ)
       call InitPDF_LHAPDF(tab%grid, tab%tab(:,:,iQ), LHAsub, Q)
    end do
  end subroutine FillPdfTable_LHAPDF

  !---------------------------------------------------------------------
  !! Given a starting distribution, StartDist, at StartScale and an
  !! "coupling" determining the behaviour of alphas, fill in the tab with
  !! the evolution of the starting distribution. Most of the arguments
  !! have meanings similar to those in the evolution routines
  !!
  subroutine EvolvePdfTable(tab, StartScale, StartDist, &
                                    & dh, coupling, muR_Q,nloop, untie_nf)
    type(pdf_table),    intent(inout) :: tab
    real(dp),           intent(in) :: StartScale
    real(dp),           intent(in) :: StartDist(0:,iflv_min:)
    type(dglap_holder), intent(in) :: dh
    type(running_coupling),  intent(in) :: coupling
    real(dp), intent(in), optional :: muR_Q
    integer,  intent(in), optional :: nloop
    logical,  intent(in), optional :: untie_nf
    !-----------------------------------------------------
    
    call EvolvePdfTableGen(tab, StartScale, dh, coupling, &
       &StartDist=StartDist, muR_Q=muR_Q, nloop=nloop, untie_nf=untie_nf)
  end subroutine EvolvePdfTable
  
  !---------------------------------------------------------------------
  !! Given a starting scale, precalculate the evolution operators that
  !! are needed to generate the full tab from a distribution at that 
  !! starting scale . Most of the arguments have meanings similar
  !! to those in the evolution routines
  !!
  subroutine PreEvolvePdfTable(tab,StartScale,dh,coupling,muR_Q,nloop,untie_nf)
    type(pdf_table),    intent(inout) :: tab
    real(dp),           intent(in) :: StartScale
    type(dglap_holder), intent(in) :: dh
    type(running_coupling),    intent(in) :: coupling
    real(dp), intent(in), optional :: muR_Q
    integer,  intent(in), optional :: nloop
    logical,  intent(in), optional :: untie_nf
    !-----------------------------------------------------
    
    call EvolvePdfTableGen(tab, StartScale, dh, coupling, &
       &precalc = .true., muR_Q=muR_Q, nloop = nloop, untie_nf=untie_nf)
  end subroutine PreEvolvePdfTable

  !---------------------------------------------------------------------
  !! Given a starting distribution, StartDist, and assuming that the
  !! pre-evolution has been carried out for the tab, then generate the
  !! contents of the table from the StartDist.
  !!
  !! NB: alphas will not be updated relative to last true evolution or
  !!     preevolution. So if at any point since doing the preevolution,
  !!     an evolution has been carried out with a different alphas, then 
  !!     as2pi will not be correct here. SHOULD THIS BE FIXED?
  !!
  subroutine pdftab_InitTabEvolve_frompre(tab, StartDist)
    type(pdf_table),    intent(inout) :: tab
    real(dp),           intent(in) :: StartDist(0:,iflv_min:)
    !-----------------------------------------------------
    real(dp) :: dist(0:ubound(StartDist,1),iflv_min:ubound(StartDist,2))
    integer :: iQ
    
    if (.not. associated(tab%evops)) call wae_error(&
         &'pdftab_InitTabEvolve_frompre',&
         &'No precalculated evolution is available')

    dist = StartDist
    do iQ = tab%StartScale_iQlo, 0, -1
       dist = tab%evops(iQ) .conv. dist
       tab%tab(:,:,iQ) = dist
    end do
    dist = StartDist
    do iQ = tab%StartScale_iQlo+1, tab%nQ
       dist = tab%evops(iQ) .conv. dist
       tab%tab(:,:,iQ) = dist
    end do
    
  end subroutine pdftab_InitTabEvolve_frompre

  !---------------------------------------------------------------------
  !! General internal routine, which serves both to carry out evolution
  !! of a parton distribution AND to determine the evops. That way
  !! everything to do with evolution is concentrated in a single location.
  !!
  !! For the user, this routine should probably be accessed via the more
  !! specific routines, EvolvePdfTable and pdftab_PreCalc.
  !!
  !! Given a starting distribution, StartDist, at StartScale and an
  !! "coupling" determining the behaviour of alphas, fill in the tab with
  !! the evolution of the starting distribution.
  !!
  subroutine EvolvePdfTableGen(tab, StartScale, dh, coupling, &
       &StartDist, precalc, muR_Q, nloop, untie_nf)
    type(pdf_table),           intent(inout) :: tab
    real(dp),               intent(in)    :: StartScale
    type(dglap_holder),     intent(in)    :: dh
    type(running_coupling), intent(in)    :: coupling
    !real(dp), optional, intent(in) :: StartDist(0:,iflv_min:)
    real(dp),     optional, intent(in)    :: StartDist(:,:)
    logical,      optional, intent(in)    :: precalc
    real(dp),     optional, intent(in)    :: muR_Q
    integer,      optional, intent(in)    :: nloop
    logical,      optional, intent(in)    :: untie_nf
    !-----------------------------------------------------
    real(dp), allocatable :: dist(:,:)
    real(dp) :: lnlnQ_norm, lnlnQ, Q_init, Q_end, last_Q
    integer :: i, iQ_lo, iQ_hi
    logical :: precalc_lcl

    precalc_lcl = default_or_opt(.false.,precalc)
    if (precalc_lcl) then
       if (associated(tab%evops)) call wae_error('EvolvePdfTableGen',&
            &'tab%evops has already been calculated. Delete the tab first,',&
            &'if you want to recalculated it.')
       allocate(tab%evops(0:tab%nQ))
    end if
    
    
    lnlnQ = lnln(tab,StartScale)
    call request_iQrange(tab, lnlnQ, 1, iQ_lo, iQ_hi, lnlnQ_norm)
    tab%StartScale = StartScale
    tab%StartScale_iQlo = iQ_lo
    ! force this...
    iQ_hi = iQ_lo + 1
    !write(0,*) iQ_lo, iQ_hi
    !lnlnQ_norm_start = lnln_norm(tab,StartScale)


    if (present(StartDist)) then
       allocate(dist(size(StartDist,1),size(StartDist,2)))
       dist = StartDist
    end if
    
    last_Q = StartScale
    !do i = floor(lnlnQ_norm_start), 0, -1
    do i = iQ_lo, 0, -1
       Q_init = last_Q
       Q_end = invlnln(tab,tab%lnlnQ_vals(i))
       !write(0,*) 'doing ev from ',Q_init,' to', Q_end

       if (present(StartDist)) then
          call EvolvePDF(dh, dist, &
               & coupling, Q_init, Q_end, muR_Q, nloop, untie_nf)
          tab%tab(:,:,i) = dist
       end if
       
       if (precalc_lcl) call EvolveGeneric(dh, &
            &coupling, Q_init, Q_end, evop=tab%evops(i), &
            &muR_Q = muR_Q, nloop = nloop, untie_nf = untie_nf)

       if (tab%nf_info_associated) tab%as2pi(i) = Value(coupling,Q_end)/twopi
       last_Q = Q_end
    end do
    

    if (present(StartDist)) dist = StartDist
    last_Q = StartScale
    !do i = ceiling(lnlnQ_norm_start), tab%nQ
    do i = iQ_hi, tab%nQ
       Q_init = last_Q
       Q_end = invlnln(tab,tab%lnlnQ_vals(i))
       !write(0,*) 'doing ev from ',Q_init,' to', Q_end

       !write(0,*) 'doing ev from ',Q_init,' to', Q_end
       if (present(StartDist)) then
          call EvolvePDF(dh, dist,&
               & coupling, Q_init, Q_end, muR_Q, nloop, untie_nf)
          tab%tab(:,:,i) = dist
       end if
       
       if (precalc_lcl) call EvolveGeneric(dh, &
            &coupling, Q_init, Q_end, evop=tab%evops(i), &
            &muR_Q = muR_Q, nloop = nloop, untie_nf = untie_nf)

       if (tab%nf_info_associated) tab%as2pi(i) = Value(coupling,Q_end)/twopi
       last_Q = Q_end
    end do
    
    if (present(StartDist)) deallocate(dist)
  end subroutine EvolvePdfTableGen


  !--------------------------------------------------------------------
  !! Sets the vector val(iflv_min:iflv_max) for the PDF at this
  !! y=ln1/x,Q.
  subroutine EvalPdfTable_yQ(tab,y,Q,val)
    type(pdf_table), intent(in) :: tab
    real(dp),     intent(in) :: y, Q
    real(dp),    intent(out) :: val(iflv_min:)
    !----------------------------------------
    real(dp) :: lnlnQ_wgts(0:tab%lnlnQ_order)
    real(dp), pointer :: y_wgts(:)
    real(dp), allocatable :: wgts(:,:)
    integer :: ilnlnQ_lo, ilnlnQ_hi, nQ,iylo, iQ, iflv
    integer, save :: warn_id = warn_id_INIT

    if (ubound(val,dim=1) < iflv_max) call wae_error('pdftab_ValTab',&
         &'upper bound of val is too low', intval=ubound(val,dim=1))

    !-- y weights taken care of elsewhere....
    call WgtGridQuant(tab%grid, y, iylo, y_wgts)

    !-- new routine for getting Q weights; ilnlnQ_lo > ilnlnQ_hi is the
    !   signal for Q being out of range
    call get_lnlnQ_wgts(tab, Q, lnlnQ_wgts, ilnlnQ_lo, ilnlnQ_hi)
    if (ilnlnQ_lo > ilnlnQ_lo) then
      val = zero
      return
    end if
    nQ = ilnlnQ_hi - ilnlnQ_lo

    allocate(wgts(lbound(y_wgts,dim=1):ubound(y_wgts,dim=1),0:nQ))
    do iQ = 0, nQ
       wgts(:,iQ) = y_wgts(:) * lnlnQ_wgts(iQ)
    end do

    !-- is this order more efficient, or should we not bother to
    !   calculate wgts?
    do iflv = iflv_min, iflv_max
       val(iflv) = sum(wgts*tab%tab(iylo:iylo+size(y_wgts)-1,&
            & iflv,ilnlnQ_lo:ilnlnQ_hi))
    end do
    !write(0,*) ilnlnQ_lo, ilnlnQ_hi, real(lnlnQ_wgts), val(1)
    
    deallocate(y_wgts, wgts)
  end subroutine EvalPdfTable_yQ


  !----------------------------------------------------------------
  !! sets the vector val(iflv_min:iflv_max) for the PDF at this x,Q.
  subroutine EvalPdfTable_xQ(tab,x,Q,val)
    type(pdf_table), intent(in) :: tab
    real(dp),     intent(in) :: x, Q
    real(dp),    intent(out) :: val(iflv_min:)
    call EvalPdfTable_yQ(tab,-log(x),Q,val)
  end subroutine EvalPdfTable_xQ


  !--------------------------------------------------------------------
  !! Sets the pdf(0:iflv_min:) for the PDF at this value of Q.
  !!
  subroutine EvalPdfTable_Q(tab,Q,pdf)
    type(pdf_table), intent(in) :: tab
    real(dp),    intent(in)  :: Q
    real(dp),    intent(out) :: pdf(0:,iflv_min:)
    !-----------------------------------------------
    real(dp) :: lnlnQ_wgts(0:tab%lnlnQ_order)
    integer  :: ilnlnQ_lo, ilnlnQ_hi
    integer  :: iQ

    ! safety checks
    if (tab%grid%ny /= ubound(pdf,dim=1) .or.&
         & ubound(pdf,dim=2) /= ncompmax) call wae_error("EvalPdfTable_Q",&
         & "pdf argument was not of size consistent with the table")

    !-- new routine for getting Q weights; ilnlnQ_lo > ilnlnQ_hi is the
    !   signal for Q being out of range
    call get_lnlnQ_wgts(tab, Q, lnlnQ_wgts, ilnlnQ_lo, ilnlnQ_hi)
    if (ilnlnQ_lo > ilnlnQ_lo) then
      pdf = zero
      return
    end if

    ! by setting the pdf to zero we will then get whatever representation
    ! is implicity in the table once we do the additions
    pdf(:,:) = zero
    do iQ = ilnlnQ_lo, ilnlnQ_hi
       pdf(:,:) = pdf(:,:) + lnlnQ_wgts(iQ-ilnlnQ_lo) * tab%tab(:, : ,iQ)
    end do
    
  end subroutine EvalPdfTable_Q


  
  !-----------------------------------------------------------
  !! Deletes all allocated info associated with the tabulation
  subroutine pdftab_DelTab_0d(tab)
    type(pdf_table), intent(inout) :: tab
    integer :: i
    deallocate(tab%tab)
    deallocate(tab%lnlnQ_vals)
    deallocate(tab%Q_vals)
    if (tab%nf_info_associated) then
       deallocate(tab%seginfo)
       deallocate(tab%nf_int)
       deallocate(tab%as2pi)
    end if
    if (associated(tab%evops)) then
       do i = 0, tab%nQ
          call Delete(tab%evops(i))
       end do
       !write(0,*) "CURRENTLY UNABLE TO DELETE EVOP CONTENTS"
       deallocate(tab%evops)
    end if
  end subroutine pdftab_DelTab_0d
  subroutine pdftab_DelTab_1d(tab)
    type(pdf_table), intent(inout) :: tab(:)
    integer :: i
    do i = 1, size(tab)
       call pdftab_DelTab_0d(tab(i))
    end do
  end subroutine pdftab_DelTab_1d
  

  !--------------------------------------------------------------
  !! Service routine used in evaluating the PDF table so as to obtain,
  !! given Q, the relevant range of ilnlnQ values and weights for
  !! interpolating the table
  !!
  subroutine get_lnlnQ_wgts(tab, Q, lnlnQ_wgts, ilnlnQ_lo, ilnlnQ_hi)
    type(pdf_table), intent(in) :: tab
    real(dp), intent(in)    :: Q
    real(dp), intent(out)   :: lnlnQ_wgts(0:)
    integer,  intent(out)   :: ilnlnQ_lo, ilnlnQ_hi
    !------------------------------------------------
    real(dp) :: lnlnQ, lnlnQ_norm
    integer  :: nQ
    integer, save :: warn_id = warn_id_INIT

    !-- Q weights need some help in finding location etc.
    lnlnQ = lnln(tab,Q)

    if (tab%freeze_at_Qmin .and. lnlnQ < tab%lnlnQ_min) then
       lnlnQ = tab%lnlnQ_min
    else if (lnlnQ < (one-warn_tolerance)*tab%lnlnQ_min &
         &.or. lnlnQ > (one+warn_tolerance)*tab%lnlnQ_max) then
      call wae_warn(default_max_warn, warn_id, &
           &'get_lnlnQ_wgts: Q out of range; result will be set to zero; Q was:',&
           &dbleval=Q)
      ilnlnQ_lo  = 0
      ilnlnQ_hi  = 0
      lnlnQ_wgts = zero ! set this to avoid warning
      return
    end if
    
    call request_iQrange(tab,lnlnQ, tab%lnlnQ_order,&
         &               ilnlnQ_lo, ilnlnQ_hi, lnlnQ_norm)

    nQ = ilnlnQ_hi - ilnlnQ_lo
    if (nQ < ubound(lnlnQ_wgts,1)) call wae_error('get_lnlnQ_wgts',&
         & 'lnlnQ_wgts to small for requested Q interpolation')
    call uniform_interpolation_weights(lnlnQ_norm, lnlnQ_wgts(0:nQ))

  end subroutine get_lnlnQ_wgts


  !------------------------------------------------------------------
  !! Given a tab and lnlnQ value, determine the range of iQ entries,
  !! iQ_lo:iQ_hi to be used for interpolating the grid.
  !!
  !! Where possible (i.e. if there are sufficient Q values in the
  !! grid) ensure that iQ_hi-iQ_lo = nrequest; for grids with varnf,
  !! all iQ_lo:iQ_hi values should correspond to the same nf value.
  !! 
  !! return value of lnlnQ_norm = (lnlnQ - lnlnQ(iQ_lo))/dlnlnQ
  !!
  subroutine request_iQrange(tab, lnlnQ, nrequest, iQ_lo, iQ_hi, lnlnQ_norm)
    type(pdf_table), intent(in) :: tab
    real(dp),     intent(in) :: lnlnQ
    integer,      intent(in) :: nrequest
    integer,     intent(out) :: iQ_lo, iQ_hi
    real(dp),    intent(out) :: lnlnQ_norm
    !-------------------------------------
    real(dp) :: dist, distclosest
    integer  :: nfclosest, nflcl
    type(pdfseginfo), pointer :: seginfo
    if (nrequest < 1) call wae_error('request_iQrange',&
         &'nrequest should be >= 1',intval=nrequest)

    if (.not. tab%nf_info_associated) then
       lnlnQ_norm = (lnlnQ-tab%lnlnQ_min)/tab%dlnlnQ
       iQ_lo = floor(lnlnQ_norm) - nrequest/2
       iQ_lo = max(0, min(tab%nQ-nrequest,iQ_lo))
       iQ_hi = min(tab%nQ,iQ_lo + nrequest)
       lnlnQ_norm = lnlnQ_norm - iQ_lo
    else
       ! need to find range to which we are closest. There is certainly
       ! a better way of doing it...
       distclosest = 1e10_dp
       do nflcl = tab%nflo, tab%nfhi
          dist = max(zero,(lnlnQ-tab%seginfo(nflcl)%lnlnQ_hi),&
               & (tab%seginfo(nflcl)%lnlnQ_lo-lnlnQ))
          if (dist < distclosest) then
             nfclosest = nflcl
             distclosest = dist
          end if
       end do
       seginfo => tab%seginfo(nfclosest)
       if (seginfo%ilnlnQ_lo == seginfo%ilnlnQ_hi) then
          iQ_lo = seginfo%ilnlnQ_lo
          iQ_hi = iQ_lo
          lnlnQ_norm = zero
       else
          lnlnQ_norm = (lnlnQ-seginfo%lnlnQ_lo)/seginfo%dlnlnQ &
               &+ seginfo%ilnlnQ_lo
          !write(0,*) lnlnQ_norm, seginfo%dlnlnQ
          iQ_lo = floor(lnlnQ_norm) - nrequest/2
          iQ_lo = max(seginfo%ilnlnQ_lo, min(seginfo%ilnlnQ_hi-nrequest,iQ_lo))
          iQ_hi = min(seginfo%ilnlnQ_hi,iQ_lo + nrequest)
          lnlnQ_norm = lnlnQ_norm - iQ_lo
          !write(0,*) iQ_lo, iQ_hi, invlnln(tab,lnlnQ), lnlnQ_norm
       end if
    end if
  end subroutine request_iQrange
  

  !-------------------------------------
  !! conversion from Q to lnlnQ
  function lnln(tab,Q)
    type(pdf_table), intent(in) :: tab
    real(dp),     intent(in) :: Q
    real(dp)                 :: lnln

    lnln = log(log(Q/tab%lambda_eff))
  end function lnln

  !! conversion from lnlnQ to Q
  function invlnln(tab,lnlnQ)
    type(pdf_table), intent(in) :: tab
    real(dp),     intent(in) :: lnlnQ
    real(dp)                 :: invlnln

    invlnln = exp(exp(lnlnQ))*tab%lambda_eff
  end function invlnln

!OLD! function lnln_norm(tab,Q)
!OLD! type(pdf_table), intent(in) :: tab
!OLD! real(dp), intent(in) :: Q
!OLD! real(dp) :: lnln_norm
!OLD!
!OLD! if (tab%nf_info_associated) call wae_error('lnln_norm')
!OLD! lnln_norm = (lnln(tab,Q) - tab%lnlnQ_min)/tab%dlnlnQ
!OLD! end function lnln_norm
!OLD!
!OLD! function invlnln_norm(tab,lnlnQ_norm)
!OLD! type(pdf_table), intent(in) :: tab
!OLD! real(dp), intent(in) :: lnlnQ_norm
!OLD! real(dp) :: invlnln_norm
!OLD!
!OLD! if (tab%nf_info_associated) call wae_error('invlnln_norm')
!OLD! invlnln_norm = invlnln(tab, lnlnQ_norm*tab%dlnlnQ + tab%lnlnQ_min)
!OLD! end function invlnln_norm
end module pdf_tabulate
