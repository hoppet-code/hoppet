!======================================================================
!!
!! Module for providing access to variable-flavour number DGLAP
!! evolution, both directly of a pdf and via evolution operators (with
!! facilities both for setting them and using them).
!!
module evolution
  use types; use consts_dp
  use dglap_objects; use qcd_coupling
  use assertions; use qcd
  use warnings_and_errors
  implicit none
  private

  public :: EvolvePDF, EvolveGeneric
  public :: SetDefaultEvolutionDt, SetDefaultEvolutionDu
  public :: DefaultEvolutionDu

  public :: ev_MSBar2DIS, ev_evolve  ! should these be public

  !!
  !! A type that allows one to store the result of the evolution as an
  !! operator that can be "applied" to any parton distribution, eliminating 
  !! the need to repeat the whole Runge-Kutta evolution for each new PDF
  !!
  type evln_operator
     type(split_mat)              :: P
     type(mass_threshold_mat)  :: MTM ! assume we have just one of these...
     real(dp)                :: MTM_coeff, Q_init, Q_end
     logical                 :: cross_mass_threshold
     type(evln_operator), pointer :: next
  end type evln_operator
  public :: evln_operator
  public :: InitEvlnOperator

  interface Delete
     module procedure Delete_evln_operator
  end interface
  public :: Delete

  
  !!
  !! related operation and operator
  !!
  public :: ev_conv_evop
  interface operator(*)
     module procedure ev_conv_evop
  end interface
  public :: operator(*)
  interface operator(.conv.)
     module procedure ev_conv_evop
  end interface
  public  :: operator(.conv.)

  !!
  !! fact that alpha_s/2pi is small means that we can quite comfortably
  !! take large steps. Have tested with 0.4 down to 1.5 GeV
  !!   
  !! Originally steps were uniform in dlnQ^2 (dt). As of September 2004,
  !! steps are now uniform in d(alpha_s * lnQ^2) = du. They are taken
  !! "matched" for alpha_s = 0.25
  !!
  real(dp), parameter :: du_dt = 0.25_dp, dt_ballpark_def = 0.4_dp
  real(dp) :: dt_ballpark = dt_ballpark_def
  real(dp) :: du_ballpark = du_dt * dt_ballpark_def
  real(dp) :: ev_tmp_du_ballpark
  !real(dp), parameter :: dt_ballpark = 0.2

  type(split_mat),      pointer :: ev_PLO, ev_PNLO, ev_PNNLO
  type(running_coupling), pointer :: ev_ash
  integer                  :: ev_nloop
  real(dp) :: ev_muR_Q, fourpibeta0_lnmuR_Q
  logical  :: ev_untie_nf = .false.

  ! this was introduced for testing --- should now never be
  ! necessary to set it to .true.
  logical,  parameter :: ev_force_old_dt = .false.
  real(dp), parameter :: ev_minalphadiff     = 0.02_dp
  integer,  parameter :: ev_du_is_dt         = 1
  integer,  parameter :: ev_du_is_dtas_fixed = 2
  integer,  parameter :: ev_du_is_dtas_run   = 3
  integer             :: ev_du_type
  real(dp)            :: ev_u_asref, ev_u_tref, ev_u_bval, ev_du_dt
  integer,  parameter :: ev_nloop_interp = -1

contains


  !======================================================================
  !! Set the spacing in t = ln Q^2 for the evolution. Actually this is
  !! now out of date, since t is no longer the evolution variable; instead
  !! it is a variable u, where du \simeq d(alphas*t). 
  !!
  !! By default du is set to be dt * du_dt, where an effective value
  !! of du_dt = 0.25 is used.
  subroutine SetDefaultEvolutionDt(dt)
    real(dp), intent(in) :: dt
    dt_ballpark = dt
    du_ballpark = dt * du_dt
  end subroutine SetDefaultEvolutionDt

  
  !======================================================================
  !! Set the spacing in u \simeq \int dlnQ^2 alphas(Q^2) for evolution
  !! steps.
  !!
  subroutine SetDefaultEvolutionDu(du)
    real(dp), intent(in) :: du
    dt_ballpark = du / du_dt
    du_ballpark = du
  end subroutine SetDefaultEvolutionDu

  !======================================================================
  !! Returns du_ballpark
  real(dp) function DefaultEvolutionDu() result(res)
    res = du_ballpark
  end function DefaultEvolutionDu
  

  !======================================================================
  !! This serves as a simplified entry point to the EvolveGeneric
  !! routine when you just want to evolve a PDF (but don't care about
  !! evln_operators). See the EvolveGeneric routine for an explanation
  !! of the various options.
  !! 
  subroutine EvolvePDF(dh, pdf, coupling, &
       &                    Q_init, Q_end, muR_Q, nloop, untie_nf, du)
    use dglap_holders; use pdf_representation
    type(dglap_holder),     intent(in), target   :: dh
    real(dp),               intent(inout)        :: pdf(:,:)
    type(running_coupling), intent(in), target   :: coupling
    real(dp),               intent(in)           :: Q_init, Q_end
    real(dp),               intent(in), optional :: muR_Q
    integer,                intent(in), optional :: nloop
    logical,                intent(in), optional :: untie_nf
    real(dp),               intent(in), optional :: du

    call EvolveGeneric(dh, coupling, Q_init, Q_end, &
         & pdf=pdf, muR_Q = muR_Q, nloop = nloop, untie_nf = untie_nf, du=du)
  end subroutine EvolvePDF


  !======================================================================
  !! This serves as a simplified entry point to the EvolveGeneric
  !! routine when you just want to evolve a PDF (but don't care about
  !! evln_operators). See the EvolveGeneric routine for an explanation
  !! of the various options.
  !! 
  subroutine InitEvlnOperator(dh, evop, coupling, &
       &                       Q_init, Q_end, muR_Q, nloop, untie_nf, du)
    use dglap_holders; use pdf_representation
    type(dglap_holder),     intent(in), target   :: dh
    type(evln_operator),    intent(out)          :: evop
    type(running_coupling), intent(in), target   :: coupling
    real(dp),               intent(in)           :: Q_init, Q_end
    real(dp),               intent(in), optional :: muR_Q
    integer,                intent(in), optional :: nloop
    logical,                intent(in), optional :: untie_nf
    real(dp),               intent(in), optional :: du

    call EvolveGeneric(dh, coupling, Q_init, Q_end, &
         & evop=evop, muR_Q = muR_Q, nloop = nloop, untie_nf = untie_nf, du=du)
  end subroutine InitEvlnOperator
  

  !======================================================================
  !! A routine that either evolves a PDF from scale Q_init to scale
  !! Q_end, or that returns an evln_operator that corresponds to the
  !! equivalent evolution.
  !!
  !! The decision on whether to evolve in a fixed flavour-number
  !! scheme (FFNS) are a variable flavour-number scheme (VFNS) is
  !! determined from the running coupling object (i.e. the DGLAP
  !! evolution does the same thing that was done for the running
  !! coupling).
  !!
  !! If the dglap_holder has been initialised with a smaller nf range
  !! than would be needed according to the coupling flavour
  !! information, then it is that smaller nf range that is used both
  !! for the DGLAP evolution and the coupling; e.g. if dglap_holder is
  !! initialised to have 3..5 flavours, then we will have a VFNS with
  !! flavours ranging from 3..5, but in the region above the top mass
  !! we will resort to a 5-flavour FFNS.
  !!
  !! Note that this is a potentially inefficient way of getting a
  !! mixed VFNS-FFNS (because one has to override the flavour number
  !! info in the coupling, involving significant recomputations) --
  !! instead you're better off setting the coupling to have a
  !! ficititious large top mass, so that the coupling is naturally in
  !! a mixed VFNS-FFNS.
  !!
  !! dh         a holder for dglap, initialised at least up to nloop
  !!
  !! coupling   a running_coupling object
  !!
  !! Q_init, Q_end: the initial and final scales for the evolution
  !!
  !! pdf     the pdf we evolve
  !!
  !! evop       or the evolution operator we set
  !!
  !! muR_Q      the ratio of the renormalisation scale to the PDF scale
  !!
  !! nloop      specifies the number of loops to be used; by default it
  !!            is taken equal to the number of loops used in the coupling.
  !!
  !! untie_nf   if (untie_nf) [default false] then alphas is allowed to 
  !!            have its natural value for nf, even if the dglap_holder
  !!            cannot reach this nf value.
  !!
  subroutine EvolveGeneric(dh, coupling, Q_init, Q_end, &
       &                     pdf, evop, muR_Q, nloop, untie_nf, du)
    use dglap_holders; use pdf_representation
    type(dglap_holder),     intent(in), target      :: dh
    type(running_coupling), intent(in),    target   :: coupling
    real(dp),               intent(in)              :: Q_init, Q_end
    real(dp),               intent(inout), optional :: pdf(:,:)
    type(evln_operator),    intent(inout), optional, target :: evop
    real(dp),               intent(in),    optional :: muR_Q
    integer,                intent(in),    optional :: nloop
    logical,                intent(in),    optional :: untie_nf
    real(dp),               intent(in),    optional :: du
    !-----------------------------------------------------
    !real(dp)           :: pdf_ev(size(pdf,dim=1), size(pdf,dim=2))
    type(dglap_holder) :: dhcopy
    integer            :: nf_init, nf_end, nflcl, nfuse
    integer            :: shnf_init, shnf_end, direction, i
    integer            :: nfstore_dhcopy, nfstore_qcd
    type(evln_operator), pointer :: this_evop
    real(dp),          pointer :: probes(:,:,:)
    !-- ranges
    real(dp) :: QRlo, QRhi, lcl_Q_init, lcl_Q_end
    
    !-- this will just copy pointers, so that we don't
    !   alter info held in original
    dhcopy = dh
    !   We want to be able to restore things at the end... (not that
    !   it is not so much for the dhcopy [which will be deleted] as for
    !   the global qcd nf. [maybe we could drop the dhcopy call?]
    nfstore_dhcopy = dh%nf
    nfstore_qcd    = nf_int


    call ev_SetModuleConsts(coupling, muR_Q, nloop, untie_nf)

    nf_init = NfAtQ(coupling, Q_init, muM_mQ=one)
    nf_end  = NfAtQ(coupling, Q_end,  muM_mQ=one)

    !write(0,*) real(Q_init), real(Q_end), nf_init, nf_end
    
    !-- could use sign, but cannot remember subtleties of signed
    !   zeroes in f95, and do not have manual at hand...
    if (Q_end >= Q_init) then
       direction = 1
    else
       direction = -1
    end if

    call ev_limit_nf(dhcopy, nf_init, shnf_init,'initial nf')
    call ev_limit_nf(dhcopy, nf_end , shnf_end ,'final nf')


    ! in loop we may want to refer to subsequent elements in chain, 
    ! so make use of a pointer which is kept updated
    if (present(evop)) then
       this_evop => evop
    else
       nullify(this_evop)
    end if

    do nflcl = shnf_init, shnf_end, direction
       call QRangeAtNf(coupling,nflcl,QRlo,QRhi, muM_mQ=one)
       if (direction == 1) then
          lcl_Q_init = max(Q_init, QRlo)
          lcl_Q_end  = min(Q_end,  QRhi)
       else
          lcl_Q_init = min(Q_init, QRhi)
          lcl_Q_end  = max(Q_end,  QRlo)
       end if
       ! make sure that end points are correct...
       if (nflcl == shnf_init) lcl_Q_init = Q_init
       if (nflcl == shnf_end)  lcl_Q_end  = Q_end
       !-- this will also set the global (qcd) nf_int which will
       !   be used elsewhere in this module
       call SetNfDglapHolder(dhcopy, nflcl)

       ! convention: cross mass thresholds before a step in the evolution
       if (nflcl /= shnf_init) then
          ! act directly
          if (present(pdf))&
               & call ev_CrossMassThreshold(dhcopy,coupling,direction,pdf)
          ! or store the information that will be needed to act subsequently.
          if (present(evop)) call ev_CrossMassThreshold(dhcopy,&
               &                  coupling,direction,evop=this_evop)
       else if (present(evop)) then
          this_evop%cross_mass_threshold = .false.
       end if
       
       ! do evolution
       if (present(pdf)) call ev_evolve(dhcopy, &
            &pdf, coupling, lcl_Q_init, lcl_Q_end, muR_Q, nloop, untie_nf, du)

       ! do the fake evolutions that allow us to do an accelerated 
       ! evolution later on.
       if (present(evop)) then
          ! recall: memory management of probes is done automatically
          !call GetDerivedSplitMatProbes(dh%grid,probes)
          call GetDerivedSplitMatProbes(dh%grid,nflcl,probes)
          do i = 1, size(probes,3)
             call ev_evolve(dhcopy, probes(:,:,i), &
                  & coupling, lcl_Q_init,lcl_Q_end, muR_Q, nloop, untie_nf, du)
          end do
          call AllocSplitMat(dh%grid, this_evop%P, nflcl)
          call SetDerivedSplitMat(this_evop%P,probes)
          if (nflcl == shnf_end) then
             nullify(this_evop%next)
          else
             allocate(this_evop%next)
             this_evop => this_evop%next
          end if
       end if
    end do

    !-- clean up
    call SetNfDglapHolder(dhcopy, nfstore_dhcopy)
    if (nfstore_qcd /= nfstore_dhcopy) call qcd_SetNf(nfstore_qcd)
  end subroutine EvolveGeneric


  !======================================================================
  !! Cross the mass threshold in the specified direction (+-1). It is
  !! assumed that the current (module qcd) nf_int value corrresponds
  !! to the number of flavours _after_ crossing the threshold.
  !!
  !! Currently only supports mass thresholds at muF = m_H.
  !!
  !! Currently direction = -1 is not supported at NNLO -- this needs
  !! to be sorted out, with special care about temporarily resetting
  !! the number of active flavours so that it corresponds to the
  !! number above threshold, rather than the number after crossing the
  !! threshold.
  !!
  subroutine ev_CrossMassThreshold(dh,coupling,direction,pdf,evop)
    use dglap_holders; use pdf_representation; use dglap_choices
    type(dglap_holder),       intent(inout) :: dh
    type(running_coupling),   intent(in)    :: coupling
    integer,                  intent(in)    :: direction
    real(dp),                 intent(inout), optional :: pdf(:,:)
    type(evln_operator),      intent(inout), optional :: evop
    !----------------------------------------------
    integer, save      :: warn_DIS = 2, warn_Direction = 2
    real(dp) :: as2pi, muR
    integer  :: nfstore
    
    !-- CHANGE THIS IF HAVE MATCHING AT MUF/=MH
    if (ev_nloop < 3) return
    if (.not. mass_steps_on) return
    if (dh%factscheme /= factscheme_MSBar) then
       call wae_Warn(warn_DIS,&
            &'ev_CrossMassThreshold',&
            &'Factscheme is not MSBar;&
            & mass thresholds requested but not implemented')
       return
    end if
    
    nfstore = nf_int ! keep record in case we change it
    select case(direction)
    case(1)
    case(-1)
       ! nf value is that after crossing threshold; but for MTM
       ! (and other things) we need nf value above threshold, i.e. before
       ! crossing the threshold; so put in the correct value temporarily
       call SetNfDglapHolder(dh, nfstore + 1)
    case default
       call wae_error('ev_CrossMassThreshold',&
            &  'direction had unsupported value of',intval=direction)
    end select
    
    !if (direction /= 1) then
    !   call wae_Warn(max_warn,warn_Direction,&
    !        &'ev_CrossMassThreshold',&
    !        &'Direction/=1; mass thresholds requested but not implemented')
    !   return
    !end if
    
    !-- now actually do something!
    !muR   = quark_masses(nf_int) * ev_MuR_Q
    muR   = QuarkMass(coupling,nf_int) * ev_MuR_Q
    !write(0,*) 'evolution crossing threshold ', nf_int, muR

    !-- fix nf so as to be sure of getting alpha value corresponding
    !   to the desired nf value, despite proximity to threshold.
    as2pi = Value(coupling, muR, fixnf=nf_int) / twopi
    if (present(pdf)) pdf = pdf + &
         &                   (direction*as2pi**2) * (dh%MTM2 .conv. pdf)
    if (present(evop)) then
       evop%cross_mass_threshold = .true.
       evop%MTM = dh%MTM2  ! stores current nf value
       evop%MTM_coeff = (direction*as2pi**2)
    end if
    
    if (nf_int /= nfstore) call SetNfDglapHolder(dh, nfstore)

  end subroutine ev_CrossMassThreshold
  

  !======================================================================
  !! Return the action of the evln_operator on the pdfdist 
  function ev_conv_evop(evop,pdfdist) result(res)
    type(evln_operator), intent(in), target :: evop
    real(dp),          intent(in)         :: pdfdist(:,:)
    real(dp) :: res(size(pdfdist,dim=1),size(pdfdist,dim=2))
    !------------
    type(evln_operator), pointer :: this_evop
    
    this_evop => evop
    res = pdfdist
    do
       if (this_evop%cross_mass_threshold) then
          ! NB: this never eccurs on first pass
          res = res + this_evop%MTM_coeff * (this_evop%MTM .conv. res)
       end if
       
       res = this_evop%P .conv. res

       if (associated(this_evop%next)) then
          this_evop => this_evop%next
       else
          return
       end if
    end do
  end function ev_conv_evop
  


  !======================================================================
  !! Return dhnflcl, which is nflcl modified to as to be
  !! within the supported limits of dh.
  !!
  !! Thus we can evolve with 5 flavours even into the 
  !! 6 flavour region, without the whole house crashing down
  !!
  !! Beware that this sort of thing is dangerous, because
  !! alphas might have one value for nf, beta0 a different value, 
  !! and overall ca sera la pagaille (if the untie_nf option is .true.)
  subroutine ev_limit_nf(dh, nflcl, dhnflcl, nfname)
    use dglap_holders
    type(dglap_holder), intent(in)  :: dh
    integer,            intent(in)  :: nflcl
    integer,            intent(out) :: dhnflcl
    character(len=*),   intent(in)  :: nfname
    !----------------------------------------
    integer, parameter :: max_warn = 4
    integer            :: warn_id = warn_id_INIT
    integer            :: nflo, nfhi
    character(len=80)  :: nfstring
    
    nflo = lbound(dh%allP,dim=2)
    nfhi = ubound(dh%allP,dim=2)

    dhnflcl = max(min(nflcl,nfhi),nflo)
    
    if (nflcl /= dhnflcl) then
       write(nfstring,'(a,i1,a,i1,a)') " changed from ",nflcl," to ",&
            &dhnflcl,"."
       call wae_warn(max_warn, warn_id, 'ev_limit_nf: '//&
            &nfname//trim(nfstring))
    end if
  end subroutine ev_limit_nf
  


  !======================================================================
  !! Internal entry routine for evolution with a fixed number of
  !! flavours. Takes a human distribution and deals with the necessary
  !! conversion to and from "evolution" format. The evolution order
  !! is set by nloop.
  subroutine ev_evolve(dh, pdf, coupling, Q_init, Q_end, muR_Q,&
       & nloop, untie_nf, du)
    use dglap_holders; use pdf_representation
    type(dglap_holder), intent(in), target :: dh
    real(dp),        intent(inout)         :: pdf(iflv_min:,0:)
    type(running_coupling), intent(in), target    :: coupling
    real(dp),        intent(in)            :: Q_init, Q_end
    real(dp),        intent(in), optional  :: muR_Q
    integer,         intent(in), optional  :: nloop
    logical,         intent(in), optional  :: untie_nf
    real(dp),        intent(in), optional  :: du
    !-----------------------------------------------------
    real(dp) :: pdf_ev(size(pdf,dim=1), size(pdf,dim=2))
    integer  :: pdfrep

    !-- make sure that number of flavours is correct?
    !if (nf_int /= dh%prep%nf) call wae_error('ev_evolve: &
    if (nf_int /= dh%nf) call wae_error('ev_evolve: &
         &global nf and representation nf are not equal.')

    !-- put things into the right format just once (it would
    !   be done automatically by the conv_object routines, but
    !   that would be wasteful)
    pdfrep = GetPdfRep(pdf)
    if (pdfrep == pdfr_Human) then
       !call CopyHumanPdfToEvln(dh%prep, pdf, pdf_ev)
       call CopyHumanPdfToEvln(dh%nf, pdf, pdf_ev)
    else
       pdf_ev = pdf
    end if
    
    ev_ash => coupling
    call ev_SetModuleConsts(coupling, muR_Q, nloop, untie_nf)

    if (ev_nloop > dh%nloop) &
         &call wae_error('ev_evolve: dh%nloop must be >= nloop')

    ev_PLO  => dh%P_LO
    if (ev_nloop >= 2) ev_PNLO => dh%P_NLO
    if (ev_nloop >= 3) ev_PNNLO => dh%P_NNLO


    call ev_evolveLocal(pdf_ev, Q_init, Q_end)

    !-- put things back into a "human" format
    if (pdfrep == pdfr_Human) then
       !call CopyEvlnPdfToHuman(dh%prep, pdf_ev, pdf)
       call CopyEvlnPdfToHuman(dh%nf, pdf_ev, pdf)
    else
       pdf = pdf_ev
    end if
  end subroutine ev_evolve


  !---------------------------------------------------------
  !! A shortcut for setting up copies of information
  !! which may be useful module-wide
  !!
  !! Is it really needed?
  subroutine ev_SetModuleConsts(coupling, muR_Q, nloop, untie_nf, du)
    type(running_coupling),    intent(in)  :: coupling
    real(dp),        intent(in), optional  :: muR_Q
    integer,         intent(in), optional  :: nloop
    logical,         intent(in), optional  :: untie_nf
    real(dp),        intent(in), optional  :: du

    ev_nloop = default_or_opt(NumberOfLoops(coupling),nloop)
    ev_muR_Q = default_or_opt(one,muR_Q)
    ev_untie_nf = default_or_opt(.false., untie_nf)
    ev_tmp_du_ballpark = default_or_opt(du_ballpark, du)

    !* !**** Additional for interpolated splitting-function sets
    !* if (dhcopy%using_interp) then
    !*    nloop = ev_nloop_interp
    !*    if (ev_muR_Q /= one) call wae_error('ev_SetModuleConsts',&
    !*         & "muR_Q must be 1 for interpolated evln; it was ",&
    !*         & dble_val=ev_muR_Q)
  end subroutine ev_SetModuleConsts
  


  !======================================================================
  !! Takes pdf in the MSBar scheme and converts it into the DIS scheme,
  !! hopefully correctly!
  subroutine ev_MSBar2DIS(dh,pdf,coupling,Q,nloop)
    use dglap_holders; use pdf_representation; use convolution
    type(dglap_holder),     intent(in)            :: dh
    real(dp),               intent(inout)         :: pdf(:,:)
    type(running_coupling), intent(in), target    :: coupling
    real(dp),               intent(in)            :: Q
    integer,                intent(in), optional  :: nloop
    !-----------------------------------------------------
    real(dp) :: pdf_ev(size(pdf,dim=1), -ncomponents:ncomponents )
    real(dp) :: Cq_x_q(size(pdf,dim=1)), Cg_x_g(size(pdf,dim=1))
    real(dp) :: as2pi
    integer  :: id

    !-- put things into the right format
    !call CopyHumanPdfToEvln(dh%prep, pdf, pdf_ev)
    call CopyHumanPdfToEvln(dh%nf, pdf, pdf_ev)

    as2pi = Value(coupling,Q) / twopi
    ev_nloop = default_or_opt(NumberOfLoops(coupling),nloop)
    
    if (ev_nloop /= 2) call wae_error('ev_MSBar2DIS',&
         &  'number of loops was not=2 [currently only case supported]')
    
    Cq_x_q = as2pi * (dh%C2_1%q .conv. pdf_ev(:,iflv_sigma))
    !Cg_x_g = (2*dh%prep%nf * as2pi) * (dh%C2_1%g .conv. pdf_ev(:,iflv_g))
    Cg_x_g = (2*dh%nf * as2pi) * (dh%C2_1%g .conv. pdf_ev(:,iflv_g))

    pdf_ev(:,iflv_sigma) = pdf_ev(:,iflv_sigma) + Cq_x_q + Cg_x_g
    pdf_ev(:,iflv_g)     = pdf_ev(:,iflv_g)     - Cq_x_q - Cg_x_g

    !do id = -dh%prep%nf, dh%prep%nf
    do id = -dh%nf, dh%nf
       if (id == iflv_sigma .or. id == iflv_g) cycle
       pdf_ev(:,id) = pdf_ev(:,id) + as2pi * (dh%C2_1%q .conv. pdf_ev(:,id))
    end do

    !-- put things back into a "human" format
    !call CopyEvlnPdfToHuman(dh%prep, pdf_ev, pdf)
    call CopyEvlnPdfToHuman(dh%nf, pdf_ev, pdf)
    !write(0,*) 'hello', lbound(pdf), lbound

  end subroutine ev_MSBar2DIS
  

  

  !======================================================================
  !! this bit here does the actual evolution of a PDF. Could make it
  !! more sophisticated later, if necessary (e.g. variable step size).
  subroutine ev_evolveLocal(pdf, Q_init, Q_end)
    use runge_kutta
    real(dp),        intent(inout)         :: pdf(:,:)
    real(dp),        intent(in)            :: Q_init, Q_end
    !------------------------------------------------------------
    real(dp) :: t1, t2, dt, t, u1, u2, du, u, as1, as2
    integer  :: n, i
    integer :: ntot=0

    fourpibeta0_lnmuR_Q = four*pi*beta0*log(ev_muR_Q)
    !write(0,*) fourpibeta0_lnmuR_Q
    t1 = two*log(Q_init)
    t2 = two*log(Q_end)

    if (t2 == t1) return

    !-- now work out jacobians...
    as1 = ev_asval(Q_init)
    as2 = ev_asval(Q_end)

    ! allow for both signs of coupling (and sign changes)
    ! See CCN25-95 for formulae (NB: 0->1 and 1->2)
    if (ev_force_old_dt) then
       ev_du_type = ev_du_is_dt
       u1 = t1; u2 = t2
       n  = ceiling(abs(t2-t1)/(ev_tmp_du_ballpark/du_dt))
    else if (abs(as1 - as2)/max(as1,as2) < ev_minalphadiff &
         &.or. as1*as2 <= zero) then
    !else if (.true.) then
       ev_du_type = ev_du_is_dtas_fixed
       ev_du_dt = max(abs(as1),abs(as2))
       u1 = t1 * ev_du_dt; u2 = t2 * ev_du_dt
       n  = ceiling(abs(u2-u1)/ev_tmp_du_ballpark)
    else
       ev_du_type = ev_du_is_dtas_run
       ev_u_asref = as1
       ev_u_tref  = t1
       ev_u_bval  = (as1/as2-1)/(as1*(t2-t1))
       u1 = zero
       u2 = log(1+ev_u_bval*ev_u_asref*(t2-t1)) / ev_u_bval
       n  = ceiling(abs(u2-u1)/ev_tmp_du_ballpark)
       !write(0,*) ev_u_asref, ev_u_bval
    end if
    
    du = (u2 - u1)/n

    ntot = ntot + n
    !write(0,*) 'Qinit,end, nsteps', Q_init, Q_end, n, ntot
    u = u1
    do i = 1, n
       call rkstp(du, u, pdf, ev_conv)
    end do
  end subroutine ev_evolveLocal


  !======================================================================
  !! Assuming the relevant module-wide variables have been set, this
  !! returns the derivative dpdf (wrt "u" -- defined precisely
  !! elsewhere) of pdf, as specified by the DGLAP equations.
  !! 
  subroutine ev_conv(u, pdf, dpdf)
    use pdf_representation
    ! The following fails with absoft compiler: see ABSOFT_BUG.txt
    real(dp), intent(in)  :: u, pdf(:,:)
    real(dp), intent(out) :: dpdf(:,:)
    ! The following fails with the intel compiler! See INTEL_BUG.txt
    !real(dp), intent(in)  :: t, pdf(0:,-ncomponents:)
    !real(dp), intent(out) :: dpdf(0:,-ncomponents:)
    !--------------------------------------
    real(dp) :: as2pi, Q, t, jacobian
    type(split_mat) :: Pfull

    ! for analysing Intel bug
    !write(0,*) 'X',lbound(pdf),lbound(pdf,dim=1),size(pdf,dim=1)

    select case(ev_du_type)
    case(ev_du_is_dt)
       t = u
       jacobian = one
    case(ev_du_is_dtas_fixed)
       t = u / ev_du_dt
       jacobian = one / ev_du_dt
    case(ev_du_is_dtas_run)
       t = ev_u_tref + (exp(ev_u_bval*u)-one)/(ev_u_bval*ev_u_asref)
       jacobian = (1+ev_u_bval*ev_u_asref*(t-ev_u_tref))/ev_u_asref
    case default
       call wae_error('evconv: unknown ev_du_type',intval=ev_du_type)
    end select
    
    Q     = exp(half*t)
    as2pi = ev_asval(Q)/twopi
    
    select case (ev_nloop)
    case(1)
       dpdf = (jacobian * as2pi) * (ev_PLO .conv. pdf)
    case(2)
       if (fourpibeta0_lnmuR_Q /= zero) then
          call InitSplitMat(Pfull, ev_PLO, one + as2pi*fourpibeta0_lnmuR_Q)
       else
          call InitSplitMat(Pfull, ev_PLO)
       end if
       call AddWithCoeff(Pfull, ev_PNLO, as2pi)
       dpdf = (jacobian * as2pi) * (Pfull .conv. pdf)
       call Delete(Pfull)
    case(3)
       if (fourpibeta0_lnmuR_Q /= zero) then
          call InitSplitMat(Pfull, ev_PLO, one + as2pi*fourpibeta0_lnmuR_Q&
               & + (as2pi*fourpibeta0_lnmuR_Q)**2&
               & + as2pi**2*(twopi**2*beta1)*two*log(ev_muR_Q))
          call AddWithCoeff(Pfull, ev_PNLO, &
               &as2pi*(one + two*as2pi*fourpibeta0_lnmuR_Q))
          call AddWithCoeff(Pfull, ev_PNNLO, as2pi**2)
          !call wae_error('ev_conv: NNL evolution not supported with muR_Q/=1')
       else
          call InitSplitMat(Pfull, ev_PLO)
          call AddWithCoeff(Pfull, ev_PNLO, as2pi)
          call AddWithCoeff(Pfull, ev_PNNLO, as2pi**2)
       end if
       dpdf = (jacobian * as2pi) * (Pfull .conv. pdf)
       call Delete(Pfull)
    case(ev_nloop_interp)
       ! *** SetCurrentP
       !
       ! dpdf = jacobian * (dhcopy%currentP .conv. pdf)
       stop ! not yet programmed !!!
    case default
       call wae_error('ev_conv','unrecognised value for ev_nloop',&
            &         intval=ev_nloop)
    end select
  end subroutine ev_conv
  
  
  !======================================================================
  !! Given module-wide variables (including a pointer to the running
  !! coupling object), this returns alphas for the appropriate number 
  !! of flavours and at a scale Q*muR_Q.
  !!
  !! (put there because nf handling and xmuR handling involve various
  !! options, and it is convenient to have a common call that handles
  !! all options correctly).
  function ev_asval(Q) result(res)
    real(dp), intent(in) :: Q
    real(dp)             :: res
    if (.not. ev_untie_nf) then
       !-- fixnf option here will be quite restrictive. It means that
       !   if we have ev_muR_Q/=1 and variable numbers of flavours,
       !   then either ev_ash supports "extrapolation" with the same nf
       !   beyond the strictly legal region, or else it has been defined 
       !   so as to have flavour thresholds at ev_muR_Q*masses
       res   = Value(ev_ash, Q*ev_muR_Q, fixnf=nf_int)
    else
       ! sometimes (e.g. comparisons with others) it is useful to
       ! allow nf in alpha_s to have a value different from the nf
       ! being used in the splitting functions...
       res   = Value(ev_ash, Q*ev_muR_Q)
    end if

  end function ev_asval
  

  !======================================================================
  !! Delete the objects allocated in evop, including any
  !! subsiduary objects.
  recursive subroutine Delete_evln_operator(evop)
    type(evln_operator), intent(inout) :: evop

    if (associated(evop%next)) then
       call Delete_evln_operator(evop%next)
       deallocate(evop%next)
    end if
    
    ! do nothing here since MTM has not actually been allocated
    ! but rather just set "equal" to dh%MTM2 -- i.e. it just points
    ! to the contents of dh%MTM2 (except for the nf info, which is
    ! local, but does not need dealocating)
    !if (evop%cross_mass_threshold) then ! do nothing

    call Delete(evop%P)
    return
  end subroutine Delete_evln_operator
  

end module evolution

