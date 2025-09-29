!
! This include sets __HOPPET_tab_min_subgrid_ny__ and
! __HOPPET_tab_min_nQ__, which are the minimum values for the subgrid ny
! and nQ respectively for various interpolation assumptions to be
! correct, avoiding the need for runtime checks
#include "inc/pdf_tabulate.inc"

!-----------------------------------------------------------
!! Module to make tabulation easy
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
  !use pdf_tabulate_types
  !-- needed only for EvolvePdfTable [separate it out one day?]
  use qcd_coupling; use evolution; use dglap_holders
  implicit none
  private

  type pdfseginfo
  real(dp) :: lnlnQ_lo, lnlnQ_hi, dlnlnQ, inv_dlnlnQ
    integer :: ilnlnQ_lo, ilnlnQ_hi
  end type pdfseginfo
  public :: pdfseginfo

  type pdf_table
    ! basic elements of a pdf_table, common regardless of whether we
    ! additionally have the nf segments...
    logical :: allocated = .false.
    type(grid_def) :: grid
    real(dp) :: default_dlnlnQ
    real(dp) :: lnlnQ_min, lnlnQ_max, lambda_eff
    real(dp), pointer :: tab(:,:,:)
    real(dp), pointer :: lnlnQ_vals(:), Q_vals(:)
    integer :: nQ, lnlnQ_order
    logical :: freeze_at_Qmin
    ! this is useful only in absence of nf info.
    real(dp) :: dlnlnQ    

    ! a seginfo-like object when there is no nf info,
    ! makes it easier to have common code across nf and no-nf 
    ! cases
    type(pdfseginfo) :: seginfo_no_nf

    ! Stuff to do with variable nf and alpha_s; not always available.
    ! In cases with variable nf, the table will be broken into multiple
    ! segments, each one of which potentially different spacings.
    !
    logical :: nf_info_associated
    integer :: nflo, nfhi
    type(pdfseginfo), pointer :: seginfo(:)
    integer,          pointer :: nf_int(:)
    real(dp),         pointer :: as2pi(:)
    !
    ! Elements needed in case we want to do precalculation of
    ! of evolution. Not always available.
    type(evln_operator), pointer :: evops(:)
    integer                      :: StartScale_iQlo
    real(dp)                     :: StartScale

    ! for some purposes (e.g. structure function work), it is useful
    ! for the table to have a user-tuned maximum iflv value (this
    ! affects the range used for interpolation rather than for
    ! allocation of the table itself)
    integer :: tab_iflv_max = iflv_max

    ! pointers to evaluation routines; the abstract interfaces
    ! are defined below
    procedure(EvalPdfTable_yQ_interface ), pointer, nopass :: EvalPdfTable_yQ_ptr  => null()
    procedure(EvalPdfTable_yQf_interface), pointer, nopass :: EvalPdfTable_yQf_ptr => null()
    procedure(EvalPdfTable1D_yQ_interface), pointer, nopass :: EvalPdfTable1D_yQ_ptr => null()
  end type pdf_table
  public :: pdf_table


  integer :: override_npnt_y = -1
  integer :: override_order_y = -1
  integer :: override_order_Q = -1
  logical :: override_on = .false.
  

  !-- for calculating ln ln Q/lambda_eff
  real(dp), parameter :: default_lambda_eff = 0.1_dp
  real(dp), parameter :: default_dlnlnQ = 0.07_dp
  real(dp), parameter :: warn_tolerance = 1e-3_dp
  !! used in various contexts for deciding when an interval is
  !! sufficiently small that it can be ignored...
  real(dp), parameter :: min_dlnlnQ_singleQ = 1e-10_dp
  integer, public :: def_lnlnQ_order = 4
  !! the max_lnlnQ_order is used in certain hard-coded array sizes
  !! so that they can be created on the stack rather than the heap,
  !! which has a speed advantage
  integer, parameter :: max_lnlnQ_order = 6

  interface AllocPdfTable
    module procedure pdftab_AllocTab_, pdftab_AllocTab_fromorig,&
          & pdftab_AllocTab_1d, pdftab_AllocTab_fromorig_1d
  end interface

  interface AllocPdfTableWithPhoton
    module procedure pdftab_AllocTabWithPhoton_, pdftab_AllocTabWithPhoton_1d
  end interface
  public :: AllocPdfTableWithPhoton

  interface AllocPdfTableWithLeptons
    module procedure pdftab_AllocTabWithLeptons_, pdftab_AllocTabWithLeptons_1d
  end interface
  public :: AllocPdfTableWithLeptons  


  interface AddNfInfoToPdfTable
    module procedure AddNfInfoToPdfTable, pdftab_AssocNfInfo_1d
  end interface

  public :: PdfTableSetInterpPointers

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
  public :: EvolvePdfTable, PreEvolvePdfTable, EvolvePdfTableQED
  public :: EvolvePdfTableGen

  interface EvalPdfTable_yQ
    module procedure EvalPdfTable0D_yQ, EvalPdfTable1D_yQ
  end interface
  interface EvalPdfTable_xQ
    module procedure EvalPdfTable0D_xQ, EvalPdfTable1D_xQ
  end interface
  public :: EvalPdfTable0D_yQ, EvalPdfTable0D_xQ, EvalPdfTable1D_yQ, EvalPdfTable1D_xQ
  public :: EvalPdfTable_yQ, EvalPdfTable_xQ

  
  public :: EvalPdfTable_Q
  public :: EvalPdfTable_yQf, EvalPdfTable_xQf
  public :: WriteLHAPDFFromPdfTable
  public :: Delete

  public :: PdfTableOverrideInterpOrders


  public :: EvalPdfTable_yQf_order22
  public :: EvalPdfTable_yQf_order33
  public :: EvalPdfTable_yQf_order44

  public :: EvalPdfTable_yQ_order22
  public :: EvalPdfTable_yQ_order33
  public :: EvalPdfTable_yQ_order44
  public :: EvalPdfTable_yQ_any_order

  abstract interface
    subroutine EvalPdfTable_yQ_interface(tab,y,Q,val)
      use types
      import :: pdf_table ! ensures pdf_table is visible
      implicit none
      type(pdf_table), intent(in), target  :: tab
      real(dp),        intent(in)          :: y, Q
      real(dp),        intent(out)         :: val(:)
    end subroutine EvalPdfTable_yQ_interface
  end interface

  abstract interface
    function EvalPdfTable_yQf_interface(tab,y,Q,iflv) result(res)
      use types
      import :: pdf_table ! ensures pdf_table is visible
      implicit none
      type(pdf_table), intent(in), target  :: tab
      real(dp),        intent(in)          :: y, Q
      integer,         intent(in)          :: iflv
      real(dp)                             :: res
    end function EvalPdfTable_yQf_interface
  end interface

  abstract interface
    subroutine EvalPdfTable1D_yQ_interface(tab,y,Q,val)
      use types
      import :: pdf_table ! ensures pdf_table is visible
      implicit none
      type(pdf_table), intent(in), target  :: tab(:)
      real(dp),        intent(in)          :: y, Q
      real(dp),        intent(out)         :: val(:,:)
    end subroutine EvalPdfTable1D_yQ_interface
  end interface

  procedure(EvalPdfTable_yQ_interface ), pointer :: EvalPdfTable_yQ_ptr  => null()
  procedure(EvalPdfTable_yQf_interface), pointer :: EvalPdfTable_yQf_ptr => null()
  procedure(EvalPdfTable1D_yQ_interface), pointer :: EvalPdfTable1D_yQ_ptr  => null()

  public :: EvalPdfTable_yQ_interface, EvalPdfTable_yQf_interface, EvalPdfTable1D_yQ_interface

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
  !! If the pdf_table had previously been allocated, then delete its
  !! contents before the new allocations
  !!
  subroutine pdftab_AllocTab_(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table)
    type(grid_def),    intent(in)    :: grid
    type(pdf_table),   intent(inout) :: tab
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    integer,  intent(in), optional :: lnlnQ_order
    logical,  intent(in), optional :: freeze_at_Qmin
    integer,  intent(in), optional :: iflv_max_table
    !----------------------------------------------
    integer :: iQ, igd
    logical :: grids_bad

    ! clear memory if the table was already allocated
    if (tab%allocated) call Delete(tab)

    tab%grid = grid
    tab%lambda_eff = min(half*Qmin, default_lambda_eff)
    tab%lnlnQ_min = lnln(tab,Qmin)
    tab%lnlnQ_max = lnln(tab,Qmax)

    ! some of our interpolation needs a minimum number of ny points in
    ! the grid or each sub-grid, so check that
    if (associated(tab%grid%subgd)) then
      grids_bad = any(tab%grid%subgd(:)%ny < __HOPPET_tab_min_subgrid_ny__)
    else
      grids_bad = tab%grid%ny < __HOPPET_tab_min_subgrid_ny__
    end if
    if (grids_bad) call wae_error("pdftab_AllocTab_",&
                  "grid or subgrid has ny < __HOPPET_tab_min_subgrid_ny__=",intval=__HOPPET_tab_min_subgrid_ny__)

    tab%default_dlnlnQ = default_or_opt(default_dlnlnQ, dlnlnQ)
    tab%nQ = ceiling((tab%lnlnQ_max - tab%lnlnQ_min)/tab%default_dlnlnQ)
    tab%nQ = max(tab%nQ, __HOPPET_tab_min_nQ__)
    tab%dlnlnQ = (tab%lnlnQ_max - tab%lnlnQ_min)/tab%nQ

    tab%freeze_at_Qmin = default_or_opt(.false.,freeze_at_Qmin)
    tab%lnlnQ_order = default_or_opt(def_lnlnQ_order,lnlnQ_order)
    if (tab%lnlnQ_order > max_lnlnQ_order) call wae_error('pdftab_AllocTab_',&
         &'too large a value of lnlnQ_order =',intval=tab%lnlnQ_order)

    !-- by default, no extra info is given
    tab%nf_info_associated = .false.
    nullify(tab%as2pi)
    nullify(tab%nf_int)
    nullify(tab%evops)

    !write(0,*) 'pdf_table info: Number of Q bins is ',tab%nQ
    ! start off with this default
    tab%tab_iflv_max = iflv_max
    if (present(iflv_max_table)) then
      if (iflv_max_table < iflv_max) call wae_error('pdftab_AllocTab_','iflv_max_table too small (<iflv_max)',intval=iflv_max_table)
      call AllocGridQuant(grid, tab%tab, iflv_min, iflv_max_table, 0, tab%nQ)
      if (iflv_max_table /= ncompmax) tab%tab_iflv_max = iflv_max_table
    else
       call AllocPDF(grid,tab%tab,0,tab%nQ)
    end if
    allocate(tab%lnlnQ_vals(0:tab%nQ))
    allocate(tab%Q_vals(0:tab%nQ))
    do iQ = 0, tab%nQ
      tab%lnlnQ_vals(iQ) = tab%lnlnQ_min + iQ * tab%dlnlnQ
      tab%Q_vals(iQ) = invlnln(tab,tab%lnlnQ_vals(iQ))
    end do


    tab%seginfo_no_nf%lnlnQ_lo = tab%lnlnQ_min
    tab%seginfo_no_nf%lnlnQ_hi = tab%lnlnQ_max
    tab%seginfo_no_nf%dlnlnQ = tab%dlnlnQ
    tab%seginfo_no_nf%inv_dlnlnQ = 1.0_dp / tab%dlnlnQ
    tab%seginfo_no_nf%ilnlnQ_lo = 0
    tab%seginfo_no_nf%ilnlnQ_hi = tab%nQ

    call pdftab_SetupEvalPtrs(tab)

    ! label this tab as allocated
    tab%allocated = .true.
  end subroutine pdftab_AllocTab_


  
  !---------------------------------------------------------
  !! 1d overloaded version of AllocPdfTable
  subroutine pdftab_AllocTab_1d(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table)
    type(grid_def),    intent(in)     :: grid
    type(pdf_table),   intent(inout)  :: tab(:)
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    integer,  intent(in), optional :: lnlnQ_order
    logical,  intent(in), optional :: freeze_at_Qmin
    integer,  intent(in), optional :: iflv_max_table
    integer :: i

    do i = 1, size(tab)
       call pdftab_AllocTab_(grid, tab(i), Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table)
    end do
  end subroutine pdftab_AllocTab_1d


  subroutine pdftab_AllocTabWithPhoton_(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin)
    use qed_objects
    type(grid_def),    intent(in)    :: grid
    type(pdf_table),   intent(inout) :: tab
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    integer,  intent(in), optional :: lnlnQ_order
    logical,  intent(in), optional :: freeze_at_Qmin
    call AllocPdfTable(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table=ncompmaxPhoton)
  end subroutine pdftab_AllocTabWithPhoton_
  
  subroutine pdftab_AllocTabWithPhoton_1d(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin)
    use qed_objects
    type(grid_def),    intent(in)    :: grid
    type(pdf_table),   intent(inout) :: tab(:)
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    integer,  intent(in), optional :: lnlnQ_order
    logical,  intent(in), optional :: freeze_at_Qmin
    call AllocPdfTable(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table=ncompmaxPhoton)
  end subroutine pdftab_AllocTabWithPhoton_1d
  

  subroutine pdftab_AllocTabWithLeptons_(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin)
   use qed_objects
   type(grid_def),    intent(in)    :: grid
    type(pdf_table),   intent(inout) :: tab
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    integer,  intent(in), optional :: lnlnQ_order
    logical,  intent(in), optional :: freeze_at_Qmin
    call AllocPdfTable(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table=ncompmaxLeptons)
  end subroutine pdftab_AllocTabWithLeptons_
  
  subroutine pdftab_AllocTabWithLeptons_1d(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin)
    use qed_objects
    type(grid_def),    intent(in)    :: grid
    type(pdf_table),   intent(inout) :: tab(:)
    real(dp), intent(in)           :: Qmin, Qmax
    real(dp), intent(in), optional :: dlnlnQ
    integer,  intent(in), optional :: lnlnQ_order
    logical,  intent(in), optional :: freeze_at_Qmin
    call AllocPdfTable(grid, tab, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table=ncompmaxLeptons)
  end subroutine pdftab_AllocTabWithLeptons_1d


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
    type(pdf_table), intent(inout), target :: tab
    type(running_coupling), intent(in) :: coupling
    !-----------------------------------
    integer  :: nflcl, iQ_prev, iQ, iflv_max_table
    real(dp) :: Qlo, Qhi, Qhi_test
    type(pdfseginfo), pointer :: seginfo

    
    !if (tab%nf_info_associated) call wae_error('AddNfInfoToPdfTable',&
    !     &'nf info already associated: delete it first')

    ! We will be reallocating everything here, so first clean up
    iflv_max_table = ubound(tab%tab,dim=2)
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
         seginfo%ilnlnQ_hi = seginfo%ilnlnQ_lo + &
            & max(__HOPPET_tab_min_nQ__,&
                & tab%lnlnQ_order,&
                & ceiling((seginfo%lnlnQ_hi - seginfo%lnlnQ_lo)/tab%default_dlnlnQ))
         seginfo%dlnlnQ = (seginfo%lnlnQ_hi - seginfo%lnlnQ_lo)/&
              & (seginfo%ilnlnQ_hi - seginfo%ilnlnQ_lo)
      end if
      !write(0,*) 'ill_hi', seginfo%ilnlnQ_hi, seginfo%dlnlnQ, &
      !     &invlnln(tab,seginfo%lnlnQ_lo),invlnln(tab,seginfo%lnlnQ_hi), tab%default_dlnlnQ
      iQ_prev = seginfo%ilnlnQ_hi
      seginfo%inv_dlnlnQ = one / seginfo%dlnlnQ
    end do
    

    ! this should not happen too often! But check it just in
    ! case...
    if (tab%seginfo(tab%nflo)%lnlnQ_lo /= tab%lnlnQ_min .or.&
         &tab%seginfo(tab%nfhi)%lnlnQ_hi /= tab%lnlnQ_max) &
         &call wae_error('AddNfInfoToPdfTable',&
         & 'mismatch in segment and global lnlnQ limits.',&
         & 'Could be due to coupling having more restricted range?')


    ! now reallocate things. 
    tab%nQ = tab%seginfo(tab%nfhi)%ilnlnQ_hi
    ! copy across the iflv_max used in the original version of the table
    call AllocGridQuant(tab%grid, tab%tab, iflv_min, iflv_max_table, 0, tab%nQ)
    allocate(tab%lnlnQ_vals(0:tab%nQ))
    allocate(tab%Q_vals(0:tab%nQ))
    allocate(tab%nf_int(0:tab%nQ))
    allocate(tab%as2pi(0:tab%nQ))

    ! carry out a check that it's always safe to use up to __HOPPET_tab_min_nQ__+1 points
    ! in the Q interpolation, even for segments with a single Q (for those segments, the
    ! interpolation weights in some high-efficiency routines will be just 1 for that iQ
    ! and 0 at higher iQ, but table must actually contain valid higher iQ points to avoid
    ! accessing invalid memory)
    do nflcl = tab%nflo, tab%nfhi      
      seginfo => tab%seginfo(nflcl)
      if (seginfo%ilnlnQ_lo == seginfo%ilnlnQ_hi .and. tab%nQ - seginfo%ilnlnQ_hi < __HOPPET_tab_min_nQ__) then
        call wae_error("AddNfInfoToPdfTable",&
              "too few ilnlnQ points for safe interpolation, above a segment with a single ilnlnQ.",&
              "Try raising the table's Qmax.",&
              "The problem occurred for nf=",intval=nflcl)
      end if
    end do


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

    call pdftab_SetupEvalPtrs(tab)

    ! REMEMBER TO COMPLETE FROM ORIG...
    tab%nf_info_associated = .true.
  end subroutine AddNfInfoToPdfTable

  !---------------------------------------------------------------
  subroutine pdftab_SetupEvalPtrs(tab)
    type(pdf_table), intent(inout) :: tab
    integer :: grid_order, npnt_y, order_y

    if (associated(tab%grid%subgd)) then
      grid_order = tab%grid%subgd(1)%order
      if (any(tab%grid%subgd(:)%order /= grid_order)) then
        ! different per-grid orders not yet supported, so set the
        ! orders to -1, which sets the pointers to null
        call PdfTableSetInterpPointers(-1,-1, &
                                   tab%EvalPdfTable_yQ_ptr, tab%EvalPdfTable_yQf_ptr, tab%EvalPdfTable1D_yQ_ptr)
      end if
    else
      grid_order = tab%grid%order
    end if

    ! convolution.f90's default is to set the interpolation order to abs(grid%order)-1
    ! with limits on the range; so table-specific interpolation, we adopt the 
    ! same defaults here
    npnt_y = min(grid_interp_npnt_max, max(grid_interp_npnt_min, abs(grid_order)))
    order_y = npnt_y - 1
    call PdfTableSetInterpPointers(order_y, tab%lnlnQ_order, &
                                   tab%EvalPdfTable_yQ_ptr, tab%EvalPdfTable_yQf_ptr, tab%EvalPdfTable1D_yQ_ptr)
  end subroutine pdftab_SetupEvalPtrs

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
    integer                      :: iflv_max_table
    
    tab = origtab
    !-- the next things are those not taken care of automatically by the assignment;
    !   we allow for arbitrary numbers of flavours
    iflv_max_table = ubound(tab%tab,dim=2)
    call AllocGridQuant(tab%grid, tab%tab, iflv_min, iflv_max_table, 0, tab%nQ)
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
    type(pdf_table), intent(inout) :: tab(:)
    type(pdf_table), intent(in)    :: origtab
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
       Q = tab%Q_vals(iQ)
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


  !---------------------------------------------------------------------
  !! Routine for creating a QED evolution table. It assumes the tabulation
  !! has been set up with the correct size, etc.
  !!
  subroutine EvolvePdfTableQED(tab, StartScale, StartDist, dh, qed_sm, &
       & coupling_qcd, coupling_qed, &
       & nloop_qcd, nqcdloop_qed_arg, with_Plq_nnloqed)
    use qed_evolution; use qed_coupling_module; use qed_objects
    type(pdf_table),           intent(inout) :: tab
    real(dp),               intent(in)    :: StartScale
    type(dglap_holder),     intent(in)    :: dh
    type(qed_split_mat),    intent(in)    :: qed_sm
    type(running_coupling), intent(in)    :: coupling_qcd
    type(qed_coupling),     intent(in)    :: coupling_qed
    real(dp),               intent(in)    :: StartDist(:,:)
    integer,                intent(in)    :: nloop_qcd, nqcdloop_qed_arg
    logical,  optional,     intent(in)    :: with_Plq_nnloqed
    !logical,      optional, intent(in)    :: untie_nf
    !-----------------------------------------------------
    real(dp), allocatable :: dist(:,:)
    real(dp) :: lnlnQ_norm, lnlnQ, Q_init, Q_end, last_Q
    integer :: i, iQ_lo, iQ_hi

    lnlnQ = lnln(tab,StartScale)
    call request_iQrange(tab, lnlnQ, 1, iQ_lo, iQ_hi, lnlnQ_norm)
    tab%StartScale = StartScale
    tab%StartScale_iQlo = iQ_lo
    ! force this...
    iQ_hi = iQ_lo + 1

    allocate(dist(size(StartDist,1),size(StartDist,2)))
    dist = StartDist
    
    last_Q = StartScale
    !do i = floor(lnlnQ_norm_start), 0, -1
    do i = iQ_lo, 0, -1
       Q_init = last_Q
       Q_end = invlnln(tab,tab%lnlnQ_vals(i))

       call QEDQCDEvolvePDF(dh, qed_sm, dist, coupling_qcd, coupling_qed,&
                   &        Q_init, Q_end, nloop_qcd, nqcdloop_qed_arg, with_Plq_nnloqed)
       tab%tab(:,:,i) = dist
       
       last_Q = Q_end
    end do

    dist = StartDist
    last_Q = StartScale
    do i = iQ_hi, tab%nQ
       Q_init = last_Q
       Q_end = invlnln(tab,tab%lnlnQ_vals(i))
       call QEDQCDEvolvePDF(dh, qed_sm, dist, coupling_qcd, coupling_qed,&
                   &        Q_init, Q_end, nloop_qcd, nqcdloop_qed_arg, with_Plq_nnloqed)
       tab%tab(:,:,i) = dist
       
       last_Q = Q_end
    end do
    
    deallocate(dist)
  end subroutine EvolvePdfTableQED

  !--------------------------------------------------------------------
  !! Override all tables' default interpolation orders (which are
  !! normally set by the grid and the olnlnQ in the table constructor)
  !!
  !! Either both orders should be >0 or both negative (negative means
  !! use the table's defaults)
  !!
  !! Specific settings give significant speed improvements. Currently:
  !! (2,2), (3,3) and (4,4), because these correspond to cases where
  !! there are specific routines with hard-coded orders
  subroutine PdfTableOverrideInterpOrders(yorder, lnlnQorder)
    integer, intent(in) :: yorder, lnlnQorder

    if (yorder * lnlnQorder <= 0) call wae_error("PdfTableOverrideInterpOrders: both orders must be >0 or <0")

    if (yorder <= 0) then
      override_order_y = -1
      override_npnt_y  = -1
    else
      override_order_y = yorder
      override_npnt_y  = yorder + 1
    end if

    override_order_Q = lnlnQorder
    if (override_order_Q <= 0) override_order_Q = -1

    ! it's enough for one to be >0 for override to be turned on
    override_on = (override_order_y > 0 .or. override_order_Q > 0)

    ! this will only set pointers if there is a hard-coded routine for
    ! the override; if either of the overrides <=0, then no pointers
    ! are set, consistent with override_on being .false.
    call PdfTableSetInterpPointers(override_order_y, override_order_Q, &
                                   EvalPdfTable_yQ_ptr, EvalPdfTable_yQf_ptr, EvalPdfTable1D_yQ_ptr)
  end subroutine PdfTableOverrideInterpOrders

  !--------------------------------------------------------------------
  !! Set up the pointers if we have a hard-coded routine for this specific
  !! order_y and order_Q, otherwise set them to null.
  subroutine PdfTableSetInterpPointers(order_y, order_Q, yQ_ptr, yQf_ptr, t1D_yQ_ptr)
    integer, intent(in) :: order_y, order_Q
    procedure(EvalPdfTable_yQ_interface ), pointer, intent(out) :: yQ_ptr
    procedure(EvalPdfTable_yQf_interface), pointer, intent(out) :: yQf_ptr
    procedure(EvalPdfTable1D_yQ_interface), pointer, intent(out) :: t1D_yQ_ptr

    if      (order_y == 2 .and. order_Q == 2) then
      yQ_ptr  => EvalPdfTable_yQ_order22
      yQf_ptr => EvalPdfTable_yQf_order22
      t1D_yQ_ptr => EvalPdfTable1D_yQ_order22
    else if (order_y == 3 .and. order_Q == 3) then
      yQ_ptr  => EvalPdfTable_yQ_order33
      yQf_ptr => EvalPdfTable_yQf_order33
      t1D_yQ_ptr => EvalPdfTable1D_yQ_order33
    else if (order_y == 4 .and. order_Q == 4) then
      yQ_ptr  => EvalPdfTable_yQ_order44
      yQf_ptr => EvalPdfTable_yQf_order44
      t1D_yQ_ptr => EvalPdfTable1D_yQ_order44
    else if (order_y == 5 .and. order_Q == 4) then
      yQ_ptr  => EvalPdfTable_yQ_order54
      yQf_ptr => EvalPdfTable_yQf_order54
      t1D_yQ_ptr => EvalPdfTable1D_yQ_order54
    else if (order_y == 6 .and. order_Q == 4) then
      yQ_ptr  => EvalPdfTable_yQ_order64
      yQf_ptr => EvalPdfTable_yQf_order64
      t1D_yQ_ptr => EvalPdfTable1D_yQ_order64
    else
      yQ_ptr  => null()
      yQf_ptr => null()
      t1D_yQ_ptr => null()
    end if
  end subroutine PdfTableSetInterpPointers

  
  !--------------------------------------------------------------------
  !! Sets the vector val(iflv_min:) for the PDF at this
  !! y=ln1/x,Q.
  !!
  !! If this is a standard QCD table then entries val(iflv_min:iflv_max)
  !! get set.
  !!
  !! If this is a generalised table, then val(:) should be as big as
  !! the number of flavours in the table and all will get set
  !! (including the "representation" flavour).
  subroutine EvalPdfTable0D_yQ(tab,y,Q,val)
    type(pdf_table), intent(in), target :: tab
    real(dp),     intent(in) :: y, Q
    real(dp),    intent(out) :: val(iflv_min:)
    !----------------------------------------
    real(dp) :: lnlnQ_wgts(0:max_lnlnQ_order)
    real(dp) :: y_wgts(0:WeightGridQuand_npnt_max-1)
    real(dp) :: wgts(0:WeightGridQuand_npnt_max-1,0:max_lnlnQ_order)
    integer :: ilnlnQ_lo, ilnlnQ_hi, nQ,iylo, iQ, iflv, iflv_max_table, npnt_y
    integer, save :: warn_id = 4!warn_id_INIT    

    if (ubound(val,dim=1) < tab%tab_iflv_max) call wae_error('pdftab_ValTab',&
         &'upper bound of val is too low', intval=ubound(val,dim=1))

    if (associated(EvalPdfTable_yQ_ptr)) then
      call EvalPdfTable_yQ_ptr(tab,y,Q,val)
    else if (associated(tab%EvalPdfTable_yQ_ptr) .and. (.not. override_on)) then
      call tab%EvalPdfTable_yQ_ptr(tab,y,Q,val)
    else
      call wae_warn(warn_id, "EvalPdfTable_yQ: y & Q interpolation orders not available hard-coded, reverting to slower routine")
      call EvalPdfTable_yQ_any_order(tab,y,Q,val)      
    end if
    
  end subroutine EvalPdfTable0D_yQ

  subroutine EvalPdfTable1D_xQ(tab,x,Q,val)
    type(pdf_table), intent(in), target :: tab(:)
    real(dp),     intent(in) :: x, Q
    real(dp),    intent(out) :: val(iflv_min:,:)
    call EvalPdfTable1D_yQ(tab, -log(x), Q, val)
  end subroutine EvalPdfTable1D_xQ
    
  !! Sets the vector val(iflv_min:,:) at this
  !! y=ln1/x,Q, for a set of tables (e.g. error members of
  !! a PDF set)
  subroutine EvalPdfTable1D_yQ(tab,y,Q,val)
    type(pdf_table), intent(in), target :: tab(:)
    real(dp),     intent(in) :: y, Q
    real(dp),    intent(out) :: val(iflv_min:,:)
    !----------------------------------------
    !real(dp) :: lnlnQ_wgts(0:max_lnlnQ_order)
    !real(dp) :: y_wgts(0:WeightGridQuand_npnt_max-1)
    !real(dp) :: wgts(0:WeightGridQuand_npnt_max-1,0:max_lnlnQ_order)
    !integer :: ilnlnQ_lo, ilnlnQ_hi, nQ,iylo, iQ, iflv, iflv_max_table, npnt_y
    integer, save :: warn_id = 4!warn_id_INIT    
    integer :: itab

    if (ubound(val,dim=1) < tab(1)%tab_iflv_max) call wae_error('pdftab_ValTab',&
         &'upper bound of val is too low', intval=ubound(val,dim=1))
    if (size(val,2) < size(tab)) then
      call wae_error("EvalPdfTable1D_yQ","size of res second dim is smaller than tab array, size(res,2)=",intval=size(val,2))
    end if


    if (associated(EvalPdfTable1D_yQ_ptr)) then
      call EvalPdfTable1D_yQ_ptr(tab,y,Q,val)
    else if (associated(tab(1)%EvalPdfTable_yQ_ptr) .and. (.not. override_on)) then
      do itab = 2, size(tab)
         if (.not. associated(tab(itab)%EvalPdfTable_yQ_ptr, tab(1)%EvalPdfTable_yQ_ptr)) then
            call wae_error("EvalPdfTable1D_yQ","not all tables have EvalPdfTable_yQ_ptr associated")
         end if
      end do
      call tab(1)%EvalPdfTable1D_yQ_ptr(tab,y,Q,val)
    else
      call wae_warn(warn_id, "EvalPdfTable1D_yQ: y & Q interpolation orders not available hard-coded, reverting to slower routine")

      do itab = 1, size(tab)
         call EvalPdfTable_yQ_any_order(tab(itab),y,Q,val(:,itab))      
      end do
    end if
    
  end subroutine EvalPdfTable1D_yQ


  !! subsidiary routine that handles arbitrary interpolation orders.
  !! Note that this is quite a bit slower than the specialized versions
  !! (order2, etc.)
  subroutine EvalPdfTable_yQ_any_order(tab,y,Q,val)
    type(pdf_table), intent(in), target :: tab
    real(dp),     intent(in) :: y, Q
    real(dp),    intent(out) :: val(iflv_min:)
    !----------------------------------------
    real(dp) :: lnlnQ_wgts(0:max_lnlnQ_order)
    real(dp) :: y_wgts(0:WeightGridQuand_npnt_max-1)
    real(dp) :: wgts(0:WeightGridQuand_npnt_max-1,0:max_lnlnQ_order)
    integer :: ilnlnQ_lo, ilnlnQ_hi, nQ,iylo, iQ, iflv, iflv_max_table, npnt_y
    integer, save :: warn_id = warn_id_INIT    

    !-- y weights taken care of elsewhere....
    call WgtGridQuant_noalloc(tab%grid, y, iylo, y_wgts, npnt_y, npnt_in = override_npnt_y)

    !-- new routine for getting Q weights; ilnlnQ_lo > ilnlnQ_hi is the
    !   signal for Q being out of range
    call get_lnlnQ_wgts(tab, Q, lnlnQ_wgts(0:tab%lnlnQ_order), ilnlnQ_lo, ilnlnQ_hi)
    nQ = ilnlnQ_hi - ilnlnQ_lo

    ! diagnostics
    !print *, "any_order: Q weights:", lnlnQ_wgts(0:nQ), ilnlnQ_lo
    !print *, "any_order: y weights:", y_wgts(0:npnt_y-1), iylo

    !-- is this order more efficient, or should we not bother to
    !   calculate wgts? Not calculating it would imply significantly
    !   more operations, so probably good to stay as is
    !
    ! we specialise certain specific cases, because this appears to
    ! help the compiler with its optimisation.
    ! 
    ! timings below done 2024-11-29, on MacM2 gfortran-14.2 -O3 with 
    ! prec_and_timing -nrep 1 -nxQ 5000 -outputgrid -dy 0.05 -order -6 -nopreev -4grids -dlnlnQ 0.01 -du 0.01 -olnlnQ 2 -yinterp 2 -nrep-eval 10000000
    !
    ! (NB: nQ == olnlnQ ; npnt_y == y interpolation order + 1; so cubic -olnlnQ 3 -yinterp 3  gives npnt_y == 4, nQ==3, which means 4 actual points for each...)
    if (npnt_y == 4 .and. nQ == 3) then
      ! ! specialising this saves about 4ns (M2Pro -O3)
      ! do iQ = 0, 3
      !    wgts(:npnt_y-1,iQ) = y_wgts(0:npnt_y-1) * lnlnQ_wgts(iQ)
      ! end do
      ! !! 110ns/call with this specialisation (M2Pro)   -> 186ns/call with -O2 -g (RelWithDebInfo)
      ! !! 165ns/call without this specialisation
      ! do iflv = iflv_min, tab%tab_iflv_max
      !    val(iflv) = sum(wgts(0:3,0:3) &
      !                    * tab%tab(iylo:iylo+3, iflv,ilnlnQ_lo:ilnlnQ_lo+3))
      ! end do

      !! 2025-09-03 (GPS+AK): this variant seems to be the fastest, saving a further 10ns (M2Pro,O3)
      !!                      getting it to about 99ns
      !! Note that further unrolling the loop manually didn't seem to bring any advantages
      do iflv = iflv_min, tab%tab_iflv_max
        val(iflv) = sum(y_wgts(0:3) * tab%tab(iylo:iylo+3, iflv,ilnlnQ_lo  )) * lnlnQ_wgts(0) &
                +   sum(y_wgts(0:3) * tab%tab(iylo:iylo+3, iflv,ilnlnQ_lo+1)) * lnlnQ_wgts(1) &
                +   sum(y_wgts(0:3) * tab%tab(iylo:iylo+3, iflv,ilnlnQ_lo+2)) * lnlnQ_wgts(2) &
                +   sum(y_wgts(0:3) * tab%tab(iylo:iylo+3, iflv,ilnlnQ_lo+3)) * lnlnQ_wgts(3) 
      end do
    else if (npnt_y == 3 .and. nQ == 2) then
      !do iQ = 0, 2
      !   wgts(0:2,iQ) = y_wgts(0:2) * lnlnQ_wgts(iQ)
      !end do
      !  85ns/call with this specialisation    -> 150ns/call with -O2 -g (RelWithDebInfo)
      ! 123ns/call without this specialisation -> 160ns/call with -O2 -g (RelWithDebInfo)
      !do iflv = iflv_min, tab%tab_iflv_max
      !      val(iflv) = sum(wgts(0:2,0:2) &
      !                      * tab%tab(iylo:iylo+2, iflv,ilnlnQ_lo:ilnlnQ_lo+2))
      !end do
      !! 2025-09-03 (GPS+AK): the following variant shaves off a further couple of percent
      do iflv = iflv_min, tab%tab_iflv_max
        val(iflv) = sum(y_wgts(0:2) * tab%tab(iylo:iylo+2, iflv,ilnlnQ_lo  )) * lnlnQ_wgts(0) &
                +   sum(y_wgts(0:2) * tab%tab(iylo:iylo+2, iflv,ilnlnQ_lo+1)) * lnlnQ_wgts(1) &
                +   sum(y_wgts(0:2) * tab%tab(iylo:iylo+2, iflv,ilnlnQ_lo+2)) * lnlnQ_wgts(2) 
      end do
    else
      do iQ = 0, nQ
         wgts(:npnt_y-1,iQ) = y_wgts(0:npnt_y-1) * lnlnQ_wgts(iQ)
      end do
      ! 250ns/call with olnlnQ==4 and yinterp==6 -> 290ns/call with -O2 -g (RelWithDebInfo)
      do iflv = iflv_min, tab%tab_iflv_max
          val(iflv) = sum(wgts(0:npnt_y-1,0:nQ) &
                          * tab%tab(iylo:iylo+npnt_y-1, iflv,ilnlnQ_lo:ilnlnQ_hi))
      end do
       !write(0,*) ilnlnQ_lo, ilnlnQ_hi, real(lnlnQ_wgts), val(1)
   end if
    
  end subroutine EvalPdfTable_yQ_any_order

  !--------------------------------------------------------------------
  !! Sets the vector val(iflv_min:) for the PDF at this
  !! y=ln1/x,Q.
  !!
  !! If this is a standard QCD table then entries val(iflv_min:iflv_max)
  !! get set.
  !!
  !! If this is a generalised table, then val(:) should be as big as
  !! the number of flavours in the table and all will get set
  !! (including the "representation" flavour).
  function EvalPdfTable_yQf(tab,y,Q,iflv) result(val)
   type(pdf_table), intent(in), target :: tab
   real(dp),     intent(in) :: y, Q
   integer,      intent(in) :: iflv
   real(dp)                 :: val
   !----------------------------------------
   !real(dp), allocatable :: wgts(:,:)
   !real(dp) :: lnlnQ_wgts(0:tab%lnlnQ_order)
   real(dp) :: lnlnQ_wgts(0:max_lnlnQ_order)
   real(dp) :: y_wgts(0:WeightGridQuand_npnt_max-1)
   integer :: ilnlnQ_lo, ilnlnQ_hi, nQ,iylo, iQ, npnt_y
   integer, save :: warn_id = 5
   real(dp) :: wgts(0:WeightGridQuand_npnt_max-1,0:max_lnlnQ_order)    

   if (iflv > tab%tab_iflv_max) call wae_error('pdftab_ValTab',&
        &'iflv is too high', intval=iflv)
   if (iflv < iflv_min) call wae_error('pdftab_ValTab',&
        &'iflv is too low', intval=iflv)

   if (associated(EvalPdfTable_yQf_ptr)) then
     val = EvalPdfTable_yQf_ptr(tab,y,Q,iflv)
     return
   else if (associated(tab%EvalPdfTable_yQf_ptr) .and. (.not. override_on)) then
     val = tab%EvalPdfTable_yQf_ptr(tab,y,Q,iflv)
     return
   end if

   call wae_warn(warn_id, "EvalPdfTable_yQf: y & Q interpolation orders not available hard-coded, reverting to slower routine")

   !-- y weights taken care of elsewhere....
   call WgtGridQuant_noalloc(tab%grid, y, iylo, y_wgts, npnt_y, override_npnt_y)

   !-- new routine for getting Q weights; ilnlnQ_lo > ilnlnQ_hi is the
   !   signal for Q being out of range
   call get_lnlnQ_wgts(tab, Q, lnlnQ_wgts(0:tab%lnlnQ_order), ilnlnQ_lo, ilnlnQ_hi)
   nQ = ilnlnQ_hi - ilnlnQ_lo

   !write(6,*) "general: iylo = ", iylo, "y_wgts = ", y_wgts(0:npnt_y-1)
   !write(6,*) "general: ilnlnQ = ", ilnlnQ_lo, "lnlnQ_wgts = ", lnlnQ_wgts(0:nQ)
   !write(6,*) "general: Qvals = ", tab%Q_vals(ilnlnQ_lo:ilnlnQ_hi)


   ! 2024-11-29 (Mac M2, gfortran 14.2) explored various options here,
   ! including going via a 2D wgts vector, and having an if statement to
   ! choose between versions with hard-coded limits for specific npnt_y
   ! and nQ values (as in EvalPdfTable_yQ). The conclusion is that
   ! within measurement error, the simple code below remains about as
   ! efficient as anything else (if not more so, by O(5-10%)).
   val = zero
   do iQ = 0, nQ
         val = val + lnlnQ_wgts(iQ) * sum(y_wgts(0:npnt_y-1) * &
             & tab%tab(iylo:iylo+npnt_y-1, iflv,ilnlnQ_lo+iQ))
   end do
 
  end function EvalPdfTable_yQf

  !----------------------------------------------------------------
  !! sets the vector val(iflv_min:iflv_max) for the PDF at this x,Q.
  subroutine EvalPdfTable0D_xQ(tab,x,Q,val)
    type(pdf_table), intent(in), target :: tab
    real(dp),     intent(in) :: x, Q
    real(dp),    intent(out) :: val(iflv_min:)
    call EvalPdfTable_yQ(tab,-log(x),Q,val)
  end subroutine EvalPdfTable0D_xQ

  !----------------------------------------------------------------
  !! sets the vector val(iflv_min:iflv_max) for the PDF at this x,Q.
  function EvalPdfTable_xQf(tab,x,Q,iflv) result(val)
    type(pdf_table), intent(in), target :: tab
    real(dp),     intent(in) :: x, Q
    integer,      intent(in) :: iflv
    real(dp) :: val
    val = EvalPdfTable_yQf(tab,-log(x),Q,iflv)
  end function EvalPdfTable_xQf

  ! 2025-09-08 attempt to explore potential for speedup in PDF
  ! evaluation when we hard-code and inline as much as possible
  function EvalPdfTable_yQf_order2_hc(tab,y,Q,iflv) result(val)
    use interpolation_coeffs
    type(pdf_table), intent(in), target :: tab
    real(dp),        intent(in)         :: y, Q
    integer,         intent(in)         :: iflv
    real(dp)                            :: val
    !----------------------------------------
    integer, parameter :: NN = 2, halfNN=0
    real(dp) :: y_wgts(0:NN), lnlnQ_wgts(0:NN), lnlnQ
    integer  :: i_nf, iylo, ilnlnQ, igd
    real(dp) :: ynorm, lnlnQ_norm
    type(grid_def),   pointer :: subgd
    type(pdfseginfo), pointer :: seginfo

    !------ First deal with the y interpolation ------
    ! For now, to get a sense of the maximal possible speed, we 
    ! assume we have a standard 4-part locked grid -- we will revisit this later
    ! once we have a better picture of speed
    ! With that assumption, we know tab%grid%subgd(4) has largest ymax
    if (y > tab%grid%subgd(4)%ymax .or. y < 0) then
      call wae_error("EvalPdfTable_yQf_order2","y did not satisfy 0 <= y <= ymax, with y=",dbleval=y)
    endif
    do igd = 4, 2, -1
      if (y > tab%grid%subgd(igd-1)%ymax) exit
    end do

    subgd => tab%grid%subgd(igd)
    ynorm = y / subgd%dy
    iylo = int(ynorm) - halfNN
    if (iylo < 0) iylo = 0
    if (iylo + NN > subgd%ny) iylo = subgd%ny - NN
    call fill_interp_weights2(ynorm - iylo, y_wgts)
    iylo = iylo + tab%grid%subiy(igd)

    !----- next deal with the Q interpolation
    lnlnQ = lnln(tab,Q)
    if (lnlnQ < tab%lnlnQ_min) then
      lnlnQ = tab%lnlnQ_min
    else if (lnlnQ > tab%lnlnQ_max) then
      call wae_error("EvalPdfTable_yQf_order2","Q was too large",dbleval=Q)
    endif
    if (.not. tab%nf_info_associated) call wae_error("EvalPdfTable_yQf_order2",&
         &"tab%nf_info_associated was not set")
    do i_nf = tab%nflo, tab%nfhi-1
      if (lnlnQ < tab%seginfo(i_nf)%lnlnQ_lo) exit
    end do

    seginfo => tab%seginfo(i_nf)
    lnlnQ_norm = (lnlnQ - seginfo%lnlnQ_lo) / seginfo%dlnlnQ
    if (seginfo%ilnlnQ_hi - seginfo%ilnlnQ_lo < NN) then
      call wae_error("EvalPdfTable_yQf_order2","not enough points in Q segment")
    end if
    ilnlnQ = int(lnlnQ_norm) + seginfo%ilnlnQ_lo - halfNN
    if (ilnlnQ    < seginfo%ilnlnQ_lo) ilnlnQ = seginfo%ilnlnQ_lo
    if (ilnlnQ+NN > seginfo%ilnlnQ_hi) ilnlnQ = seginfo%ilnlnQ_hi - NN
    call fill_interp_weights2(lnlnQ_norm - (ilnlnQ-seginfo%ilnlnQ_lo), lnlnQ_wgts)

    val = sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0) &
        + sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ+1) * y_wgts) * lnlnQ_wgts(1) &
        + sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ+2) * y_wgts) * lnlnQ_wgts(2)
  end function EvalPdfTable_yQf_order2_hc

  !! given a table and a y value, sets up a pointer (grid_ptr) to the relevant
  !! grid or sub-grid, as well as the iy_offset to the start of that grid
  !! in the y dimension of the table.
  !!
  !! Note that this function overlaps substantially with conv_BestISub,
  !! Because of its location it can be inlined in other pdf_tabulate.f90 
  !! routines, which gives a small but relevant speed advantage
  subroutine tab_get_grid_ptr(tab, y, grid_ptr, iy_offset)
    type(pdf_table), intent(in), target :: tab
    real(dp), intent(in) :: y
    type(grid_def), pointer :: grid_ptr
    integer, intent(out) :: iy_offset
    !----------------------------------------
    type(grid_def), pointer :: grid
    integer :: igd

    grid => tab%grid
    ! note that relative to conv_BestIsub, we bail out if we are
    ! beyond grid%ymax
    if (y > grid%ymax .or. y < 0) then
      call wae_error("tab_get_grid_ptr","y did not satisfy 0 <= y <= ymax, with y=",dbleval=y)
    endif

    if (associated(grid%subgd)) then
      if (grid%locked) then
        do igd = size(grid%subgd), 2, -1
          if (y > grid%subgd(igd-1)%ymax) exit
        end do
      else
        igd = minloc(grid%subgd%ymax,dim=1,mask=(grid%subgd%ymax>=y))
      end if
      grid_ptr => grid%subgd(igd)
      iy_offset = grid%subiy(igd)
    else
      grid_ptr => grid
      iy_offset = 0
    endif
  end subroutine tab_get_grid_ptr

  !subroutine tab_get_seginfo_ptr(tab, lnlnQ, seginfo_ptr)
  !  type(pdf_table), intent(in), target :: tab
  !  real(dp),        intent(inout)      :: lnlnQ
  !  type(pdfseginfo), pointer           :: seginfo_ptr
  !  !----------------------------------------
  !  integer :: i_nf
  !
  !  if (lnlnQ < tab%lnlnQ_min) then
  !    lnlnQ = tab%lnlnQ_min
  !  else if (lnlnQ > tab%lnlnQ_max) then
  !    call wae_error("EvalPdfTable_yQ_orderNNNN","Q was too large",dbleval=invlnln(tab,lnlnQ))
  !  endif
  !  if (.not. tab%nf_info_associated) then
  !    seginfo_ptr => tab%seginfo_no_nf
  !  else
  !    do i_nf = tab%nflo, tab%nfhi-1
  !      if (lnlnQ < tab%seginfo(i_nf)%lnlnQ_lo) exit
  !    end do
  !    seginfo_ptr => tab%seginfo(i_nf)
  !  end if
  !end subroutine tab_get_seginfo_ptr

  ! the following code loads the appropriate PDF evaluation subroutine
  ! with the order hard-coded to NNNN
#define __HOPPET_InterpOrder__ 2
#include "inc/pdf_tabulate_OrderNNN.F90"

#define __HOPPET_InterpOrder__ 3
#include "inc/pdf_tabulate_OrderNNN.F90"

#define __HOPPET_InterpOrder__ 4
#include "inc/pdf_tabulate_OrderNNN.F90"

#define __HOPPET_InterpOrderY__ 5
#define __HOPPET_InterpOrderQ__ 4
#include "inc/pdf_tabulate_OrderNNN.F90"

#define __HOPPET_InterpOrderY__ 6
#define __HOPPET_InterpOrderQ__ 4
#include "inc/pdf_tabulate_OrderNNN.F90"

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
         & ubound(pdf,dim=2) /= ubound(tab%tab,dim=2)) call wae_error("EvalPdfTable_Q",&
         & "pdf argument was not of size consistent with the table")

    !-- new routine for getting Q weights; ilnlnQ_lo > ilnlnQ_hi is the
    !   signal for Q being out of range
    call get_lnlnQ_wgts(tab, Q, lnlnQ_wgts, ilnlnQ_lo, ilnlnQ_hi)
    if (ilnlnQ_lo > ilnlnQ_hi) then
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

  !----------------------------------------------------------------------
  !! Write a PDF table to an LHAPDF-format file, i.e. the .dat files
  !! used by LHAPDF
  !!
  !! - table: the pdf_table object that is to be written
  !!
  !! - coupling: the coupling object associated with this table
  !! 
  !! - basename: the basename of the LHAPDF file; the full name
  !!   will be basename_nnnn.dat, where nnnn is given by pdf_index.
  !!   
  !! - pdf_index is the index that will be used for the LHAPDF .dat file
  !!   If it is zero, the routing will also output the .info file and
  !!   and the .dat will be declared a central member; otherwise 
  !!   the .dat file will be declared an error member.
  !!
  !! - iy_increment (optional, default 1):
  !!   * anything positive: chooses the x-points based on table's own
  !!     grid; the value indicates the increment to use in going up in y.
  !!     The essence of the rule is that every time y goes up and iy has
  !!     gone up by at least iy_increment, it outputs an x point.
  !!   * -1: uses a fixed distribution of x points with a coarse ln x spacing
  !!     at small x, and a finer spacing at larger x. LHAPDF's resulting interpolation
  !!     accuracy should be about $10^{-4}$ for $x<0.5$, $10^{-3}$ for $x<0.9$
  !!
  !! - flav_indices (optional): indicates which indices to use from the PDF table
  !!
  !! - flav_pdg_ids (optional): indicates the corresponding PDG IDs
  !!
  !! - flav_rescale (optional): specifies by how much each flavour should be
  !!   rescaled (usually by one, unless you want to use an entry that
  !!   has a sum of flavour and anti-flavour and so rescale by 50% to get
  !!   just the flavour)
  !!
  !! If flav_indices and flav_pdg_ids are not supplied, these are chosen
  !! by establishing whether the table has photons and leptons based on the
  !! the upper bound of the table's flavour index dimension.
  !!
  subroutine WriteLHAPDFFromPdfTable(table, coupling, basename, pdf_index, &
                                 & iy_increment, &
                                 & flav_indices, flav_pdg_ids, flav_rescale)                                
    use qed_objects
    type(pdf_table), target, intent(in) :: table
    type(running_coupling), intent(in) :: coupling
    character(len=*),       intent(in) :: basename
    integer,                intent(in) :: pdf_index
    integer,  optional,     intent(in) :: iy_increment
    integer,  optional,     intent(in) :: flav_indices(:), flav_pdg_ids(:)
    real(dp), optional,     intent(in) :: flav_rescale(:)
    !------------------------------------
    integer          :: dat_unit, info_unit
    character(len=4) :: pdf_index_str
    character(len=:), allocatable :: dat_file, info_file
    real(dp) :: flavs(lbound(table%tab,2):ubound(table%tab,2))
    real(dp) :: xVals_orig(0:table%grid%ny), val, yVal
    real(dp) :: zeta, zetamax, dzeta
    ! these parameters appear to give OK accuracy, <=10^{-4} for x<0.5
    ! and <=10^{-3} for x<0.9, at least for Q=100
    real(dp), parameter :: dzeta_def = 0.30_dp, zeta_a = 10.0_dp, zeta_b = 5.0_dp
    integer  :: iyVals(0:table%grid%ny), nx, nx_max, iy, ix, iQ, iseg, iQlo, iQhi, ipdg, iflv
    type(pdfseginfo), pointer :: seginfos(:)
    integer  :: iy_inc, i

    integer,  allocatable :: lcl_flav_indices(:), lcl_flav_pdg_ids(:)
    real(dp), allocatable :: lcl_flav_rescale(:)
    ! for choosing the xvalues output
    real(dp), allocatable :: xVals_tmp(:), xVals(:)

    iy_inc = default_or_opt(1, iy_increment)

    if (present(flav_indices)) then
      if (.not. present(flav_pdg_ids)) call wae_error("WritePdfTableLHAPDF",&
                                    "flav_pdg_ids must be present if flav_indices is present")
      lcl_flav_indices = flav_indices
      lcl_flav_pdg_ids = flav_pdg_ids
    else
      if (present(flav_pdg_ids) .or. present(flav_rescale)) then
        call wae_error("WritePdfTableLHAPDF","flav_pdg_ids and flav_rescale must not be present if flav_indices is not present")
     end if
     !print*, "Entered here", ubound(table%tab,2), ncompmax, ncompmaxPhoton, ncompmaxLeptons
      if (ubound(table%tab,2) == ncompmax) then
        lcl_flav_indices = (/-6,-5,-4,-3,-2,-1, 0,1,2,3,4,5,6/)
        lcl_flav_pdg_ids = (/-6,-5,-4,-3,-2,-1,21,1,2,3,4,5,6/)
      else if (ubound(table%tab,2) == ncompmaxPhoton) then
        lcl_flav_indices = (/-6,-5,-4,-3,-2,-1, 0,1,2,3,4,5,6, 8/)
        lcl_flav_pdg_ids = (/-6,-5,-4,-3,-2,-1,21,1,2,3,4,5,6,22/)
      else if (ubound(table%tab,2) == ncompmaxLeptons) then
        ! recall that our internal representation has 
        lcl_flav_indices = (/-6,-5,-4,-3,-2,-1, 0,1,2,3,4,5,6, 8,&
             & iflv_tau, iflv_muon, iflv_electron, iflv_electron,&
             & iflv_muon, iflv_tau/)
        lcl_flav_pdg_ids = (/-6,-5,-4,-3,-2,-1,21,1,2,3,4,5,6,22,    &
             &  -15,       -13,           -11,            11,       &
             & 13,       15/)
        lcl_flav_rescale = (/ (one, i=1,14), half, half, half, half,&
             & half, half/)
      end if

    end if
    if (present(flav_rescale)) lcl_flav_rescale = flav_rescale

    write(pdf_index_str,'(i4.4)') pdf_index
    dat_unit = 0
    dat_file = trim(basename)//"_"//pdf_index_str//'.dat'
    open(newunit=dat_unit, file=dat_file, status='replace')
    if (dat_unit == 0) call wae_error("WritePdfTableLHAPDF","Could not open "//dat_file)
    write(6,'(a)') "# Writing LHAPDF dat file "//trim(dat_file)

    ! First work out what x values to use.
    ! In this simple version, we simply take the grid points (eliminating duplications)
    xVals_orig = xValues(table%grid)
    write(dat_unit,'(a)') "# data file generated by HOPPET"
    write(dat_unit,'(a,i4,a,f10.5)') '# iy_increment = ', iy_inc, ", dy(base) = ", table%grid%dy
    if (iy_inc > 0) then
      xVals_tmp = MonotonicUniqueXValues(table%grid)
      nx = size(xVals_tmp) / iy_inc
      nx_max = nx
      ! make sure we go all the way to the lowest point in x
      if (nx * iy_inc < size(xVals_tmp)) nx_max = nx + 1
      allocate(xVals(1:nx_max))
      xVals(1:nx) = xVals_tmp(::iy_inc)
      if (nx < nx_max) xVals(nx_max) = xVals_tmp(ubound(xVals_tmp,1))
    !else if (iy_inc == 0) then
    !  ! use the PDF4LHC15 points; note that we need to reverse their order
    !  do iy = size(pdf4lhc_yvalues), 1, -1
    !     yVal = pdf4lhc_yvalues(iy)
    !     if (yVal < table%grid%ymax) then
    !        nx_max = nx_max + 1
    !        xVals(nx_max) = exp(-yVal)
    !     end if
    !  end do
    !  ! finish off with the last grid point
    !  nx_max = nx_max + 1
    !  xVals(nx_max) = exp(-table%grid%ymax)
    else if (iy_inc == -1) then
      write(dat_unit,'(a,f10.5,a,f10.5,a,f10.5)') &
          & '# dzeta_def = ', dzeta_def,&
          & ', zeta_a = ', zeta_a,&
          & ', zeta_b = ', zeta_b

      ! use our own smoothly changing spacing
      zetamax = zetaext_of_y(table%grid%ymax, zeta_a, zeta_b)
      nx_max = ceiling(zetamax/dzeta_def)
      dzeta = zetamax / nx_max
      allocate(xVals(0:nx_max))
      do ix = 0, nx_max
        xVals(ix) = exp(-y_of_zetaext(dzeta*ix, zeta_a, zeta_b))
      end do
    else
      call wae_error('WritePdfTableLHAPDF','value of iy_increment could not be handled; it was',&
          &intval = iy_inc)
    end if
    
    !write(0,*) 'nx = ', nx_max+1, ', nQ = ', size(table%Q_vals)

    ! the official header
    write(dat_unit,'(a,a)') 'PdfType: ',trim(merge("central","error  ",pdf_index==0))!trim(pdf_type)
    write(dat_unit,'(a)'  ) 'Format: lhagrid1'
    write(dat_unit,'(a)'  ) '---'

    ! handle different seginfo scenarios
    if (table%nf_info_associated) then
       ! if the table has seginfo, then just point to it
       seginfos => table%seginfo
    else
       ! otherwise create a fictitious seginfo with the info we need...
       allocate(seginfos(1:1))
       seginfos(1)%ilnlnQ_lo = lbound(table%Q_vals,1)
       seginfos(1)%ilnlnQ_hi = ubound(table%Q_vals,1)
    end if

    ! now loop over the actual or "invented" seginfos
    do iseg = lbound(seginfos,1), ubound(seginfos,1)
       ! first we output info about x structure, Q structure and flavours
       write(dat_unit,'(4000es14.7)') xVals(ubound(xVals,1):lbound(xVals,1):-1)
       iQlo = seginfos(iseg)%ilnlnQ_lo
       iQhi = seginfos(iseg)%ilnlnQ_hi
       write(dat_unit,'(4000es14.7)') table%Q_vals(iQlo:iQhi)
       write(dat_unit,'(100i4)')      lcl_flav_pdg_ids

       ! then write out the PDF itself
       ! Note that LHAPDF wants _increasing_ x values
       do ix = ubound(xVals,1), lbound(xVals,1), -1
          do iQ = iQlo, iQhi
             flavs = table%tab(:,:,iQ) .atx. (xVals(ix).with.table%grid)
             do ipdg = 1, ubound(lcl_flav_pdg_ids,1)
                val = flavs(lcl_flav_indices(ipdg))
                !print*, val, lcl_flav_rescale(ipdg), lcl_flav_indices(ipdg), ipdg
                if (allocated(lcl_flav_rescale)) val = val * lcl_flav_rescale(lcl_flav_indices(ipdg))
                if (val == zero) then
                   write(dat_unit,'(a)',advance='no') '  0'
                else
                   write(dat_unit,'(es15.7)',advance='no') val
                end if
             end do
             !stop
             write(dat_unit,'(a)') ''
          end do
       end do

       ! and finish with a yaml subdocument separator
       write(dat_unit,'(a)') '---'
    end do
    close(dat_unit)

    if (pdf_index == 0) then
      info_file = trim(basename)//'.info'
      info_unit = 0
      open(newunit=info_unit, file=info_file, status='replace')
      if (info_unit == 0) call wae_error("WritePdfTableLHAPDF: could not open info file  "//info_file)
      write(6,'(a)') "# Writing LHAPDF info file "//trim(info_file)

      write(info_unit,'(a)') 'SetDesc: <SetDesc>'
      write(info_unit,'(a)') 'SetIndex: 00000'
      write(info_unit,'(a)') 'Authors: <Authors>'
      write(info_unit,'(a)') 'Reference: <Reference>'
      write(info_unit,'(a)') 'Format: lhagrid1'
      write(info_unit,'(a)') 'DataVersion: 1'
      write(info_unit,'(a)') 'NumMembers: 1'
      write(info_unit,'(a)') 'Particle: 2212'
      write(info_unit,'(a)',advance="no") 'Flavors: '
      call write_yaml_array(info_unit, '(i3)', lcl_flav_pdg_ids(:), 15)
      write(info_unit,'(a,i2)') 'OrderQCD: ', NumberOfLoops(coupling)-1
      if (size(seginfos) == 1) then
        write(info_unit,'(a)') 'FlavorScheme: fixed'
      else
        write(info_unit,'(a)') 'FlavorScheme: variable'
      endif 
      write(info_unit,'(a,i2)') 'NumFlavors: ', table%nfhi
      write(info_unit,'(a)') 'ErrorType: <ErrorType>'
      write(info_unit,'(a,es22.14)') 'XMin: ', minval(xVals)
      write(info_unit,'(a,es22.14)') 'XMax: ', maxval(xVals)
      write(info_unit,'(a,es22.14)') 'QMin: ', minval(table%Q_vals)
      write(info_unit,'(a,es22.14)') 'QMax: ', maxval(table%Q_vals)
      write(info_unit,'(a)') 'MZ: 0.911870E+02'
      write(info_unit,'(a)') 'MUp: 0'
      write(info_unit,'(a)') 'MDown: 0'
      write(info_unit,'(a)') 'MStrange: 0'
      write(info_unit,'(a,es22.14)') 'MCharm: ',  QuarkMass(coupling,4)
      write(info_unit,'(a,es22.14)') 'MBottom: ', QuarkMass(coupling,5)
      write(info_unit,'(a,es22.14)') 'MTop: ',    QuarkMass(coupling,6)
      write(info_unit,'(a,es22.14)') 'AlphaS_MZ:', value(coupling, 91.11870_dp)
      write(info_unit,'(a)') 'AlphaS_OrderQCD: 2'
      write(info_unit,'(a)') 'AlphaS_Type: ipol'
      ! allocate, fill and output Q_vals and alphas_vals, 
      block
        real(dp), allocatable :: Q_vals(:), alphas_vals(:)
        Q_vals = table%Q_vals
        allocate(alphas_vals(lbound(Q_vals,1):ubound(Q_vals,1)))
        if (table%nf_info_associated) then
          ! near thresholds, by default the pdf_table chooses Q_vals that are very slightly
          ! below/above the thresholds; but LHAPDF wants two identical Q values
          ! as a signal that it should do separate interpolation above/below
          do iseg = lbound(seginfos,1), ubound(seginfos,1)-1
            Q_vals(seginfos(iseg)%ilnlnQ_hi) = 0.5_dp*(&
              table%Q_vals(seginfos(iseg)%ilnlnQ_hi) + table%Q_vals(seginfos(iseg+1)%ilnlnQ_lo))
            Q_vals(seginfos(iseg+1)%ilnlnQ_lo) = Q_vals(seginfos(iseg)%ilnlnQ_hi)
          end do
        end if
        do iQ = lbound(Q_vals,1), ubound(Q_vals,1)
          if (table%nf_info_associated) then
            alphas_vals(iQ) = value(coupling, Q_vals(iQ), fixnf=table%nf_int(iQ))
          else
            alphas_vals(iQ) = value(coupling, table%Q_vals(iQ))
          end if
        end do
      write(info_unit,'(a)', advance="no") 'AlphaS_Qs: '
      call write_yaml_array(info_unit, '(es22.14)', Q_vals, 4)
        write(info_unit,'(a)', advance="no") 'AlphaS_Vals: '
        call write_yaml_array(info_unit, '(es22.14)', alphas_vals,4)
      end block
      close(info_unit)
    end if

    ! cleaning
    if (.not. table%nf_info_associated) deallocate(seginfos)
  end subroutine WriteLHAPDFFromPdfTable

  !! write a yaml array (integer or real(dp)), limiting the number
  !! of elements to max_per_line on each line
  subroutine write_yaml_array(iunit, format, array, max_per_line)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: format
    class(*),         intent(in) :: array(:)
    integer,          intent(in) :: max_per_line
    !------------
    integer :: i
    write(iunit,'(a)',advance="no") "["
    do i = 1, size(array)
      if (i > 1) write(iunit,'(a)',advance="no") ","
      if (mod(i,max_per_line) == 0) write(iunit,'(a)') ""
      ! now write the value
      select type (array)
      type is (integer)
        write(iunit,format,advance="no") array(i)
      type is (real(dp))
        write(iunit,format,advance="no") array(i)
      class default
        call wae_error("write_yaml_array: unsupported array type")
      end select
    end do
    write(iunit,'(a)') "]"
  end subroutine write_yaml_array

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
    ! reset the information about allocation so that we don't
    ! try to delete something that has already been deleted
    tab%allocated = .false.
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
    use hoppet_to_string
    type(pdf_table), intent(in) :: tab
    real(dp), intent(in)    :: Q
    real(dp), intent(out)   :: lnlnQ_wgts(0:)
    integer,  intent(out)   :: ilnlnQ_lo, ilnlnQ_hi
    !------------------------------------------------
    real(dp) :: lnlnQ, lnlnQ_norm
    integer  :: nQ, nQ_request
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

    
    nQ_request = tab%lnlnQ_order    
    if (override_order_Q > 0) then
      nQ_request = override_order_Q
    else
      nQ_request = tab%lnlnQ_order
    end if
    if (nQ_request > ubound(lnlnQ_wgts,1)) then
      call wae_error('get_lnlnQ_wgts',&
         & 'lnlnQ_wgts too small (ubound='//trim(to_string(ubound(lnlnQ_wgts,1)))&
         //') for requested Q interpolation (nQ_request='//trim(to_string(nQ_request))//')')
    end if
    call request_iQrange(tab,lnlnQ, nQ_request,&
         &               ilnlnQ_lo, ilnlnQ_hi, lnlnQ_norm)

    nQ = ilnlnQ_hi - ilnlnQ_lo
    if (nQ == 0) then
       if (abs(lnlnQ - tab%lnlnQ_vals(ilnlnQ_lo)) < two * min_dlnlnQ_singleQ) then
          lnlnQ_wgts(0) = 1
       else
         call wae_error('get_lnlnQ_wgts',&
            & 'nQ=0 but lnlnQ='//trim(to_string(lnlnQ))//&
            & ' not close enough to tab%lnlnQ_vals(ilnlnQ_lo)(='//trim(to_string(tab%lnlnQ_vals(ilnlnQ_lo)))//')')
       end if

    !else if (nQ < ubound(lnlnQ_wgts,1)) then
    !  ! nQ can be zero if the table had a very narrow range of Q values (< O(min_dlnlnQ_singleQ))
    !  if (nQ == 0 .and. abs(lnlnQ - tab%lnlnQ_vals(ilnlnQ_lo)) < two * min_dlnlnQ_singleQ) then
    !     lnlnQ_wgts(0) = 1
    !  else
    !     !write(6,*) "Q=", Q, 'nQ = ',nQ, 'lnlnQ_wgts = ',ubound(lnlnQ_wgts,1), tab%nf_int
    !     call wae_error('get_lnlnQ_wgts',&
    !        & 'lnlnQ_wgts too small (ubound='//trim(to_string(ubound(lnlnQ_wgts,1)))&
    !        //') for requested Q interpolation (nQ='//trim(to_string(nQ))//')')
    !  end if
    end if
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
    type(pdf_table), intent(in), target :: tab
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

    if (Q < tab%lambda_eff) then
       lnln = -1e300_dp
    else
       lnln = log(log(Q/tab%lambda_eff))
    end if
  end function lnln

  !! conversion from lnlnQ to Q
  function invlnln(tab,lnlnQ)
    type(pdf_table), intent(in) :: tab
    real(dp),     intent(in) :: lnlnQ
    real(dp)                 :: invlnln

    invlnln = exp(exp(lnlnQ))*tab%lambda_eff
  end function invlnln

  !-----------------------------------------------------------------
  !! return zeta = ln 1/x + a*(1-x) + b * (1-x)^4 (x = exp(-y))
  function zetaext_of_y(y, a, b) result(zeta)
    real(dp), intent(in) :: y, a, b
    real(dp)             :: zeta, x
    x = exp(-y)
    zeta = y + a*(one - x) + b * (one - x**4)
  end function zetaext_of_y

  
  !-----------------------------------------------------------------
  !! return inversion of zetaext = ln 1/x + a*(1-x) +b*(1-x^4) (x = exp(-y))
  function y_of_zetaext(zetaext, a, b) result(y)
    real(dp), intent(in) :: zetaext, a, b
    real(dp)             :: y, x, diff_from_zero, deriv
    integer              :: iter
    real(dp), parameter  :: eps = 1e-12_dp
    integer,  parameter  :: maxiter = 100

    ! starting condition (and soln if a = 0 and b=0)
    y = zetaext 
    if (a /= zero .or. b /= zero) then
      do iter = 0, maxiter
        x = exp(-y);
        diff_from_zero = zetaext - (y + a*(one-x) + b*(one-x**4));
        ! we have found good solution
        if (abs(diff_from_zero) < eps) exit
        deriv = -one  - a*x - four*b*x**4;
        y = y - diff_from_zero / deriv;
      end do
    end if
    
    if (iter > maxiter) write(0,*) "y_of_zeta reached maxiter"

  end function y_of_zetaext


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
