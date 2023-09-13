!======================================================================
!! Module that provides tools for constructing a general (vector of)
!! running couplings, including coupled evolution and mass thresholds. 
!! 
!! Aims here are:
!!   - reasonable flexbility
!!   - decent performance on evaluating a coupling at a given scale
!!   - guaranteed accuracy
!!
!! Initialisation speed will however not be good since the routine
!! explicitly solves the differential equation for evolution to build
!! a grid from which it then interpolates the coupling. While better
!! methods exist for simple cases (e.g. inverting an equation for
!! lambda in terms of alpha up to NNLO in QCD), the method adopted
!! here has the advantage of being easily generalised.
!! 
!! To "create" a coupling it is necessary first to create the
!! CouplingTerm objects that go to making a CouplingDefinition
!!
module couplings
  use types; use consts_dp
  use warnings_and_errors
  use runge_kutta
  use assertions
  implicit none
  private

  !! this provides the definition for a term in coupling's evolution
  type CouplingTerm
     integer               :: itarget
     integer, allocatable  :: alpha_powers(:)
     real(dp)              :: coeff
  end type CouplingTerm
  public :: CouplingTerm
  

  !! this provides the definition for a coupling's evolution
  type CouplingDefinition
     integer            :: n_couplings, nf_lo, nf_hi
     logical, allocatable :: use_reciprocal(:)
     type(CouplingTerm), allocatable :: beta_terms(:,:)
     type(CouplingTerm), allocatable :: threshold_terms(:,:)
     type(CouplingTerm), allocatable :: reverse_threshold_terms(:,:)
     real(dp)  :: muMatch_M !! matching scale / mass scale
  end type CouplingDefinition
  public :: CouplingDefinition

  
  type(CouplingTerm),       pointer :: current_operation_glb(:)
  type(CouplingDefinition), pointer :: cpl_def_glb
  integer                           :: nf_glb

  
  !! Sub-Segment (fixed nf) in the evolution of a coupling
  type CouplingSubSegment
     real(dp) :: dt, t_lo, t_hi
     real(dp), allocatable :: ra(:,:)
  end type CouplingSubSegment
  public :: CouplingSubSegment

  !! Segment (fixed nf, multiple sub-segments) in the evolution of a coupling
  type CouplingSegment
     real(dp) :: dt, t_lo, t_hi, min_sub_dt
     type(CouplingSubSegment), allocatable :: subseg(:)
  end type CouplingSegment
  public :: CouplingSegment


  !--------- go from 0.5 GeV to over 10^20 GeV
  !real(dp), parameter :: default_t_lo = -1.3862944 !-3.2188758_dp
  !real(dp), parameter :: default_t_hi = 93.0_dp 
  ! for speed and stability, go from 0.71 GeV to about 10^6 GeV
  ! (as of 18/08/05 c. 3 ms on 1.7GHz P4 with Lahey)
  real(dp), parameter :: default_t_lo = -0.71_dp
  real(dp), parameter :: default_t_hi = 30.0_dp
  real(dp), parameter :: initial_dt = 0.4_dp, default_eps = 1e-8_dp
  real(dp), parameter :: subsegment_size = 2.5_dp
  ! something that we can safely divide by but that should
  ! not cause problems when adding things on
  real(dp), parameter :: fake_zero_dt = 1e-200_dp
  ! fraction of dt by which we will consider it to be safe to extrapolate
  real(dp), parameter :: safe_extrapolation = 0.02
  ! NOTE: this is not free, but should be >= than order of Runge-Kutta 
  ! routine used...
  integer,  parameter :: interp_ord = 4
  integer :: iii
  real(dp), parameter :: interp_resc_t(0:interp_ord)=(/(iii,iii=0,interp_ord)/)

  !! an (array of) couplings
  type Coupling
     type(CouplingDefinition)            :: def
     real(dp)                            :: t_lo, t_hi, eps
     !! segment_boundaries(nf) is tofQ(mass of object(nf)) -- or else
     !! t_lo, t_hi for outer reach.
     real(dp),              allocatable  :: segment_boundaries(:)
     type(CouplingSegment), allocatable  :: segments(:)
     integer                             :: n_points
  end type Coupling
  public :: Coupling

  public :: new_CouplingTerm, new_CouplingDefinition
  public :: new_Coupling
  public :: coupling_array, coupling_value
  public :: nf_at_Q, nf_at_t, nf_range, Q_Range_at_nf, t_Range_at_nf

  interface apply_operation
     module procedure apply_operation_alpha
  end interface
  public :: apply_operation
  public :: tofQ, Qoft

  ! some private interfaces
  interface segment_alpha
     module procedure segment_alpha_all, segment_alpha_sctn
  end interface

  interface in_range
     module procedure in_range_int, in_range_dp
  end interface

contains

  !======================================================================
  !! Returns a term for renorm group evolution of the form
  !!
  !! d alpha(itarget) / dln Q^2 = sum(alpha(:) ** alpha_powers(:)) * coeff
  !!
  !! Analogously, for mass thresholds we will have:
  !!
  !! delta alpha(itarget) = sum(alpha(:) ** alpha_powers(:)) * coeff
  !!
  function new_CouplingTerm(itarget,alpha_powers,coeff) result(term)
    integer,  intent(in) :: itarget, alpha_powers(:)
    real(dp), intent(in) :: coeff
    type(CouplingTerm) :: term
    term%itarget = itarget
    allocate(term%alpha_powers(size(alpha_powers)))
    term%alpha_powers = alpha_powers
    term%coeff        = coeff
  end function new_CouplingTerm


  !======================================================================
  !! Returns a definition object for a coupling which should have a 
  !! range of flavors going from nf_lo to nf_hi, with 
  !! beta_terms(:,nf_lo:nf_hi) for describing the renormalisation group
  !! evolution. 
  !!
  !! - If beta_terms is not supplied then when creating the
  !!   coupling object itself, the user will need to supply a
  !!   subprogram for calculating the derivative.
  !!
  !! - If the optional threshold_terms(:,nf_lo+1:nf_hi) are supplied then
  !!   these are used for transitions going both up and down in flavor
  !!   (using the lower-nf coupling as reference when going up in nf,
  !!   and the higher-nf coupling as reference when going down).
  !!
  !! - If at the desired accuracy that behaviour is not good enough, one
  !!   can supply reverse_threshold_terms(:,nflo+1:nf_hi) for going down
  !!   in nf, using the higher-nf coupling to compute the change.
  !!
  !! - by default evolution is carried out in 1/alpha; in some cases
  !!   this may not be viable (e.g. in certain coupled equations) and
  !!   one can then provide the array use_reciprocal which specifies
  !!   which couplings evolve in their reciprocal.
  function new_CouplingDefinition(n_couplings, &
       & nf_lo, nf_hi, beta_terms, threshold_terms,&
       & reverse_threshold_terms, muMatch_M, use_reciprocal)&
       &     result(def)
    integer,            intent(in) :: n_couplings, nf_lo, nf_hi
    type(CouplingTerm), intent(in), optional :: beta_terms(:,:)
    type(CouplingTerm), intent(in), optional :: threshold_terms(:,:)
    type(CouplingTerm), intent(in), optional :: reverse_threshold_terms(:,:)
    real(dp),           intent(in), optional :: muMatch_M
    logical,            intent(in), optional :: use_reciprocal(:)
    type(CouplingDefinition) :: def
    !---------------------------
    integer :: dummy
    
    
    def%n_couplings = n_couplings
    def%nf_lo = nf_lo
    def%nf_hi = nf_hi

    if (present(beta_terms)) then
       call check_CouplingTerms(def, beta_terms, is_threshold=.false., &
            &                   extra_info = '(while checking beta_terms)')
       allocate(def%beta_terms(size(beta_terms,1), nf_lo:nf_hi))
       def%beta_terms = beta_terms
    end if
    
    if (present(threshold_terms)) then
       call check_CouplingTerms(def, threshold_terms, is_threshold=.true., &
            &                   extra_info ='(while checking threshold_terms)')
       allocate(def%threshold_terms(size(threshold_terms,1), nf_lo+1:nf_hi))
       def%threshold_terms = threshold_terms
    end if

    if (present(reverse_threshold_terms)) then
       if (.not.present(threshold_terms)) call wae_error(&
            &'new_CouplingDefinition',&
            &'reverse_threshold_terms may only be supplied if &
            &threshold_terms also present')
       call check_CouplingTerms(def, reverse_threshold_terms, &
            & is_threshold=.true., &
            & extra_info = '(while checking reverse_threshold_terms)')
       allocate(def%reverse_threshold_terms(&
            &       size(reverse_threshold_terms,1), nf_lo+1:nf_hi))
       def%reverse_threshold_terms = reverse_threshold_terms
    end if

    allocate(def%use_reciprocal(1:def%n_couplings))
    if (present(use_reciprocal)) then
       if (size(use_reciprocal) /= n_couplings) call wae_error(&
            &'new_CouplingDefinition',&
            &'Size of use_reciprocal did not correspond to n_couplings')
       def%use_reciprocal = use_reciprocal
    else
       def%use_reciprocal = .true.
    end if
    

    def%muMatch_M = default_or_opt(one, muMatch_M)
  end function new_CouplingDefinition
  

  !======================================================================
  !! Checks a bunch of coupling terms to ensure that they make sense
  !! in the context of the (partially initialised) CouplingDefinition;
  !! tests differ (marginally) according to whether is_threshold =
  !! .true. or .false.
  !!
  subroutine check_CouplingTerms(def, terms, is_threshold, extra_info)
    type(CouplingDefinition), intent(in) :: def
    type(CouplingTerm),       intent(in) :: terms(:,:)
    logical,                  intent(in) :: is_threshold
    character(len=*),         intent(in) :: extra_info
    !------------
    character(len=100) :: our_info
    integer :: right_size_2, nf_lbound, nf, iterm
    
    if (is_threshold) then
       nf_lbound = def%nf_lo + 1
    else
       nf_lbound = def%nf_lo
    end if
    right_size_2 = def%nf_hi - nf_lbound + 1

    ! size!
    if (size(terms,2) /= right_size_2) then
       write(our_info,'(a,i0,a,i0)') &
            & 'Size (dim 2) of terms was ', size(terms,2), &
            & ' but expected ', right_size_2 
       call wae_error('check_CouplingTerms', trim(our_info), extra_info)
    end if
    
    ! sanity of range of alpha
    if (any(terms%itarget < 1 .or. terms%itarget > def%n_couplings)) then
       write(our_info,'(a,i0)') &
         &'There were target terms with indices outside expected range 1..',&
         & def%n_couplings
       call wae_error('check_CouplingTerms', trim(our_info), extra_info)
    end if

    ! sanity of size of alpha_powers
    do nf = nf_lbound, def%nf_hi
       do iterm = 1, size(terms,1)
          if (size(terms(iterm,nf-nf_lbound+1)%alpha_powers) &
               &/= def%n_couplings) then
             write(our_info,'(a,"(iterm=",i0,",nf=",i0,")",a,i0,a,i0)') &
             & "size of terms", iterm,nf,"%alpha_powers(:) should have been ",&
                 & def%n_couplings, " but was ",&
                 & size(terms(iterm,nf-nf_lbound+1)%alpha_powers)
             call wae_error('check_CouplingTerms', trim(our_info), extra_info)
          end if
       end do
    end do
  end subroutine check_CouplingTerms
  

  !======================================================================
  !! Returns a coupling type based on the definition def, with values
  !! alpha(:) at initial scale Q with masses(nf_lo+1:nf_hi).
  !!
  function new_Coupling(def,alpha,Q,masses, deriv_or_threshold, &
       &                    Q_lo, Q_hi, eps) result(cpl)
    type(CouplingDefinition), intent(in) :: def
    real(dp),              intent(in) :: alpha(:), Q, masses(def%nf_lo+1:)
    real(dp), optional,    intent(in) :: Q_lo !!default= exp(.5*default_t_lo)
    real(dp), optional,    intent(in) :: Q_hi !!default= exp(.5*default_t_hi)
    real(dp), optional,    intent(in) :: eps  !!deafult= default_eps
    type(Coupling) :: cpl
    interface
       !! function that returns either dalpha/dt (for nf_in=nf_out)
       !! or delta alpha when crossing a threshold from nf_in to nf_out
       function deriv_or_threshold(t, alpha, nf_in, nf_out)
         use types; implicit none
         real(dp), intent(in) :: t
         real(dp), intent(in) :: alpha(:)
         integer,  intent(in) :: nf_in, nf_out
         real(dp)             :: deriv_or_threshold(size(alpha))
       end function deriv_or_threshold
    end interface
    optional :: deriv_or_threshold
    !--------------------
    integer  :: dummy, nf_Q, nf
    real(dp) :: alpha_here(def%n_couplings), this_t

    cpl%def   = def
    
    dummy =assert_eq(size(alpha),def%n_couplings,'new_Coupling: size of alpha')
    dummy =assert_eq(size(masses),def%nf_hi - def%nf_lo,&
         &                                      'new_Coupling: size of masses')
    
    cpl%t_lo = tofQ(default_or_opt(Qoft(default_t_lo), Q_lo))
    cpl%t_hi = tofQ(default_or_opt(Qoft(default_t_hi), Q_hi))
    cpl%eps  = default_or_opt(default_eps, eps)

    allocate(cpl%segment_boundaries(def%nf_lo:def%nf_hi+1))
    cpl%segment_boundaries(def%nf_lo) = cpl%t_lo
    cpl%segment_boundaries(def%nf_hi+1) = cpl%t_hi
    cpl%segment_boundaries(def%nf_lo+1:def%nf_hi) = &
         & in_range(tofQ(masses(:)) + tofQ(def%muMatch_M), cpl%t_lo, cpl%t_hi)
    
    ! the following error message should now be redundant since all
    ! masses are put at edges (this means that results at edges
    ! may not be predictable).
    if (any(cpl%segment_boundaries > cpl%t_hi .or.&
         &  cpl%segment_boundaries < cpl%t_lo)) call wae_error('new_Coupling',&
         & 'A segment boundary has fallen outside Q_lo:Q_hi range')

    ! keep track of number of points we are storing
    cpl%n_points = 0
    
    allocate(cpl%segments(def%nf_lo:def%nf_hi))
    nf_Q = nf_at_Q(cpl, Q)

    ! fill "central" segment
    call fill_segment(cpl,nf_Q, alpha, tofQ(Q), deriv_or_threshold)

    ! fill segments above this one
    do nf = nf_Q + 1, def%nf_hi
       this_t = cpl%segment_boundaries(nf)
       alpha_here = segment_alpha(cpl, nf-1, this_t)

       ! go through various options
       if (present(deriv_or_threshold)) then
          alpha_here = alpha_here &
               & + deriv_or_threshold(this_t, alpha_here, nf-1, nf)
       else if (allocated(def%threshold_terms)) then
          alpha_here = alpha_here + &
               &apply_operation(def%threshold_terms(:,nf),alpha_here)
       end if
       call fill_segment(cpl,nf, alpha_here, this_t, deriv_or_threshold)
    end do
    
    ! fill segments below this one
    do nf = nf_Q - 1, def%nf_lo, -1
       this_t = cpl%segment_boundaries(nf+1)
       alpha_here = segment_alpha(cpl, nf+1, this_t)
       if (present(deriv_or_threshold)) then
          alpha_here = alpha_here &
               & + deriv_or_threshold(this_t, alpha_here, nf+1, nf)
       else if (allocated(def%reverse_threshold_terms)) then
          alpha_here = alpha_here + &
               &apply_operation(def%reverse_threshold_terms(:,nf+1),alpha_here)
       else if (allocated(def%threshold_terms)) then
          alpha_here = alpha_here - &
               &apply_operation(def%threshold_terms(:,nf+1),alpha_here)
       end if
       call fill_segment(cpl,nf, alpha_here, this_t, deriv_or_threshold)
    end do

  end function new_Coupling
  
  
  !======================================================================
  !! Fill a segment to the required precision -- algorithm used is
  !! to subdivide the segment into "bite-sized" pieces (subsegment_size)
  !! and for each subsegment to use an adaptive uniform grid.
  !! 
  !! Probably not the most efficient general approach in the world, but 
  !! should adequate for likely real-life uses.
  subroutine fill_segment(cpl, nf, alpha, t, deriv)
    type(Coupling), intent(inout), target :: cpl
    integer,        intent(in)    :: nf
    real(dp),       intent(in)    :: alpha(:), t
    interface
       !! function that returns either dalpha/dt (for nf_in=nf_out)
       !! or delta alpha when crossing a threshold from nf_in to nf_out
       function deriv(t, alpha, nf_in, nf_out)
         use types; implicit none
         real(dp), intent(in) :: t
         real(dp), intent(in) :: alpha(:)
         integer,  intent(in) :: nf_in, nf_out
         real(dp)             :: deriv(size(alpha))
       end function deriv
    end interface
    optional :: deriv
    !-------------------------------------------
    type(CouplingSegment), pointer :: seg
    integer  :: n
    integer  :: iseg_init, iseg
    real(dp) :: this_alpha(size(alpha))

    

    seg => cpl%segments(nf)
    ! maybe extend these limits later?
    seg%t_lo = cpl%segment_boundaries(nf)
    seg%t_hi = cpl%segment_boundaries(nf+1)

    ! having n >= 1 makes life simpler elswhere
    n = max(interp_ord+1,ceiling((seg%t_hi - seg%t_lo)/subsegment_size))
    seg%dt   = max(fake_zero_dt,(seg%t_hi - seg%t_lo)/n)
    allocate(seg%subseg(0:n-1))

    ! now fill up the subsegments.
    ! first their bounds
    do iseg = 0, n-1
       seg%subseg(iseg)%t_lo = seg%t_lo + iseg * seg%dt; 
       seg%subseg(iseg)%t_hi = seg%t_lo + (iseg+1) * seg%dt
       !write(0,*) nf, iseg, seg%subseg(iseg)%t_lo, seg%subseg(iseg)%t_hi
    end do

    ! then fill their contents:
    ! first do the one that contains the specific t value
    iseg_init = in_range(floor((t - seg%t_lo)/seg%dt), 0, n-1)
    
    call fill_subsegment(cpl, nf, seg%subseg(iseg_init), alpha, t, deriv)

    ! then do the ones above it
    do iseg =  iseg_init+1, n-1
       this_alpha= reciprocal_or_not(&
            &       seg%subseg(iseg-1)%ra(:,ubound(seg%subseg(iseg-1)%ra,2)))
       call fill_subsegment(cpl, nf, seg%subseg(iseg), &
            & this_alpha, seg%subseg(iseg-1)%t_hi, deriv)
    end do

    ! then below
    do iseg =  iseg_init-1, 0, -1
       this_alpha= reciprocal_or_not(seg%subseg(iseg+1)%ra(:,0))
       call fill_subsegment(cpl, nf, seg%subseg(iseg), &
            & this_alpha, seg%subseg(iseg+1)%t_lo, deriv)
    end do
    
    seg%min_sub_dt = minval(seg%subseg%dt)
  end subroutine fill_segment
  

  !======================================================================
  !! Fill this sub-segment with alpha values. The algorithm used
  !! is to do the fill based on a coarse grid and then a fine grid
  !! and if they differ by too much then halve the grid size and
  !! try again. In the end it is the finer grid that is taken.
  !!
  !! No doubt there exist better approaches -- but hopefully this
  !! will be reasonably robust.
  subroutine fill_subsegment(cpl, nf, subseg, alpha, t, deriv)
    type(Coupling), intent(inout), target   :: cpl
    type(CouplingSubSegment), intent(inout) :: subseg
    integer,        intent(in)    :: nf
    real(dp),       intent(in)    :: alpha(:), t
    interface
       !! function that returns either dalpha/dt (for nf_in=nf_out)
       !! or delta alpha when crossing a threshold from nf_in to nf_out
       function deriv(t, alpha, nf_in, nf_out)
         use types; implicit none
         real(dp), intent(in) :: t
         real(dp), intent(in) :: alpha(:)
         integer,  intent(in) :: nf_in, nf_out
         real(dp)             :: deriv(size(alpha))
       end function deriv
    end interface
    optional :: deriv
    !--------------------------
    integer  :: idiv, n
    real(dp), pointer :: coarse(:,:), fine(:,:)

    subseg%dt = initial_dt
    ! having n >= 1 makes life simpler elswhere
    n  = max(1,ceiling((subseg%t_hi - subseg%t_lo)/subseg%dt))
    subseg%dt = max(fake_zero_dt, (subseg%t_hi - subseg%t_lo)/n)
    idiv = 0
    do
       allocate(fine(size(alpha), 0:n))
       call fill_tentative_subsegment(cpl, nf, subseg, fine, alpha, t, deriv)
       if (idiv > 0) then
          !write(0,*) 'dt,diff', subseg%dt, maxval(abs(coarse-fine(:,::2)))
          !-- the funny rel + abs error
          !if (all (abs(coarse-fine(:,::2)) &
          !     & <= cpl%eps*(one+abs(fine(:,::2))))) exit
          !-- here instead just have a plain relative error
          if (all (abs(coarse-fine(:,::2)) <= cpl%eps*abs(fine(:,::2)))) exit
          deallocate(coarse)
       end if
       coarse => fine
       subseg%dt = half * subseg%dt; n = 2*n
       idiv = idiv + 1
       ! these are quite arbitrary conditions...
       if (n > 10000 .or. &
            &(subseg%dt > fake_zero_dt .and.subseg%dt < 1e-4_dp)) then
          call wae_error(&
            &'fill_subsegment',&
            &'Need too fine a grid for storing alpha',&
            &'Try reducing the required precision &
            &or restricing the evolution range')
       end if
       
    end do
    allocate(subseg%ra(size(alpha), 0:n))
    subseg%ra = fine
    !write(0,*) 'final dt is', nf, subseg%dt
    deallocate(coarse)
    deallocate(fine)
    cpl%n_points = cpl%n_points + n+1
  end subroutine fill_subsegment
  
  !======================================================================
  !! Fill the supplied array candidate for a segment given
  !! the supplied alpha and t initial conditions.
  !!
  !! It assumes that set%dt and t_lo, t_hi correspond to 
  !! the n chosen.
  subroutine fill_tentative_subsegment(cpl, nf, subseg, fine, alpha, t, deriv)
    type(Coupling), intent(inout), target :: cpl
    integer,        intent(in)    :: nf
    type(CouplingSubSegment), intent(inout) :: subseg
    real(dp),       intent(out)   :: fine(:,0:)
    real(dp),       intent(in)    :: alpha(:), t
    interface
       !! function that returns either dalpha/dt (for nf_in=nf_out)
       !! or delta alpha when crossing a threshold from nf_in to nf_out
       function deriv(t, alpha, nf_in, nf_out)
         use types; implicit none
         real(dp), intent(in) :: t
         real(dp), intent(in) :: alpha(:)
         integer,  intent(in) :: nf_in, nf_out
         real(dp)             :: deriv(size(alpha))
       end function deriv
    end interface
    optional :: deriv
    !--------
    real(dp) :: ra(size(alpha)), this_t, new_t
    integer :: dir(2), init_it(2), last_it(2), idir, it, n

    ! used within the D.E.
    if (allocated(cpl%def%beta_terms))&
         & current_operation_glb => cpl%def%beta_terms(:,nf)
    cpl_def_glb => cpl%def
    nf_glb = nf
    n = ubound(fine,2)

    init_it(1) = max(0,floor((t-subseg%t_lo)/subseg%dt))
    init_it(2) = min(n, init_it(1) + 1)
    init_it(1) = init_it(2) - 1
    dir(1:2)     = (/-1, 1/)
    last_it(1:2) = (/0, n/)

    do idir = 1, 2
       this_t = t
       ra = reciprocal_or_not(alpha)
       !write(0,*) 'doing',init_it(idir), last_it(idir), dir(idir)
       do it = init_it(idir), last_it(idir), dir(idir)
          new_t = it*subseg%dt + subseg%t_lo
          if (present(deriv)) then
             call rkstp_1d_special(new_t - this_t, this_t, ra, deriv)
          else
             call rkstp(new_t - this_t, this_t, ra, calc_dra)
          end if
          fine(:,it) = ra
       end do
    end do
    
  end subroutine fill_tentative_subsegment
  
  
  !======================================================================
  !! Calculates dra, the derivative of 1/alpha (ra)
  subroutine calc_dra(t, ra, dra)
    real(dp), intent(in)  :: t, ra(:)
    real(dp), intent(out) :: dra(:)
    !------------------
    real(dp) :: alpha(size(dra))

    alpha = reciprocal_or_not(ra)
    dra = apply_operation(current_operation_glb, alpha)
    where(cpl_def_glb%use_reciprocal)
       dra = -ra**2 * dra 
    end where
  end subroutine calc_dra

  !======================================================================
  !! Calculates dra, the derivative of 1/alpha (ra)
  subroutine calc_dra_special(t, ra, dra, deriv)
    real(dp), intent(in)  :: t, ra(:)
    real(dp), intent(out) :: dra(:)
    interface
       !! function that returns either dalpha/dt (for nf_in=nf_out)
       !! or delta alpha when crossing a threshold from nf_in to nf_out
       function deriv(t, alpha, nf_in, nf_out)
         use types; implicit none
         real(dp), intent(in) :: t
         real(dp), intent(in) :: alpha(:)
         integer,  intent(in) :: nf_in, nf_out
         real(dp)             :: deriv(size(alpha))
       end function deriv
    end interface
    !------------------
    real(dp) :: alpha(size(dra))

    alpha = reciprocal_or_not(ra)
    dra = deriv(t,alpha,nf_glb, nf_glb)
    where(cpl_def_glb%use_reciprocal)
       dra = -ra**2 * dra 
    end where

    !write(0,'(i3,f10.3,3es13.5," ",3es13.5)') nf_glb, t, ra, dra
  end subroutine calc_dra_special


  !======================================================================
  !! Returns the value of coupling icpl (default = 1)
  function coupling_array(cpl, Q, fix_nf) result(alpha)
    type(Coupling), intent(in) :: cpl
    real(dp),       intent(in) :: Q
    integer, optional, intent(in) :: fix_nf
    real(dp)                   :: alpha(cpl%def%n_couplings)
    !-----------------------
    integer  :: nf
    real(dp) :: t

    t = tofQ(Q)
    ! avoid using default_or_opt since nf_at_t(cpl,t) might add extra
    ! time delay? (Maybe this is exaggerated?)
    if (present(fix_nf)) then
       nf = fix_nf
    else
       nf = nf_at_t(cpl,t)
    end if
    alpha = segment_alpha(cpl, nf, t)
  end function coupling_array
  

  !======================================================================
  !! Returns the value of the first of the couplings (since this
  !! is often enough) -- NB not efficient if large array of couplings.
  !!
  !! On 1.7GHz P4 takes about 1.5ns with lf95.
  function coupling_value(cpl, Q, fix_nf) result(alpha)
    type(Coupling), intent(in) :: cpl
    real(dp),       intent(in) :: Q
    integer, optional, intent(in) :: fix_nf
    real(dp)                   :: alpha
    real(dp)                   :: alpha_arr(1:1)
    !-----------------------
    integer  :: nf
    real(dp) :: t

    t = tofQ(Q)
    ! avoid using default_or_opt since nf_at_t(cpl,t) might add extra
    ! time dealy? (Maybe this is exaggerated?)
    if (present(fix_nf)) then
       nf = fix_nf
    else
       nf = nf_at_t(cpl,t)
    end if
    alpha_arr = segment_alpha(cpl, nf, t, icpl_lo = 1, icpl_hi = 1)
    alpha = alpha_arr(1)
    !write(0,'(f15.10,i,f15.12)') t, nf, alpha ! TESTING!!!
  end function coupling_value


  !======================================================================
  !! Returns the value of the coupling array from the sepcific segment , 
  !! (seg_nf) at the given t value for the coupling cpl.
  function segment_alpha_all(cpl, seg_nf, t) result(alpha)
    use interpolation
    type(Coupling),    intent(in), target :: cpl
    integer,           intent(in) :: seg_nf
    real(dp),          intent(in) :: t
    real(dp)                      :: alpha(1:cpl%def%n_couplings)
    alpha = segment_alpha(cpl, seg_nf, t, 1, cpl%def%n_couplings)
  end function segment_alpha_all
  

  !======================================================================
  !! Returns the value the coupling array section (icpl_lo:icpl_hi) from
  !! the sepcific segment (seg_nf), at the given t value for the
  !! coupling cpl.
  !!
  !! Note that we have two routines segment_alpha_sctn and segment_alpha_all
  !! rather than optional arguments, so as to be able to define the
  !! size of the result without difficulty.
  function segment_alpha_sctn(cpl, seg_nf, t, icpl_lo, icpl_hi) result(alpha)
    use interpolation
    type(Coupling),    intent(in), target :: cpl
    integer,           intent(in) :: seg_nf
    real(dp),          intent(in) :: t
    integer,           intent(in) :: icpl_lo, icpl_hi
    real(dp)                      :: alpha(icpl_lo:icpl_hi)
    !-------------------------------
    real(dp) :: dalpha, t_resc
    type(CouplingSegment), pointer :: seg
    type(CouplingSubSegment), pointer :: subseg
    integer :: it_lo, it_hi, n, i, iseg
    real(dp) :: weights(0:interp_ord)
    character(len=100) :: error_msg
    
    seg => cpl%segments(seg_nf)
    ! sometimes rounding errors take us beyond the allowed limits, 
    ! especially with fix_nf option -- so allow some margin for going
    ! outside range.
    if (t > seg%t_hi + safe_extrapolation*seg%min_sub_dt&
         & .or. t < seg%t_lo - safe_extrapolation*seg%min_sub_dt) then
       write(error_msg,*) 't was', t, ', range was', seg%t_lo, seg%t_hi
       call wae_error(&
         &'segment_alpha: requested t out of segment range',&
         &trim(error_msg))
    end if
    
    iseg = in_range(floor((t-seg%t_lo)/seg%dt), 0, ubound(seg%subseg,1))
    subseg => seg%subseg(iseg)

    it_lo = max(0,floor((t - subseg%t_lo)/subseg%dt) - (interp_ord)/2)
    it_hi = min(ubound(subseg%ra,2),it_lo+interp_ord)
    it_lo = max(0,it_hi - interp_ord)

    n = it_hi - it_lo

    t_resc = (t - (subseg%t_lo+it_lo*subseg%dt))/subseg%dt

    ! efficient, but how accurate?
    ! alternatively: go via polint?
    call uniform_interpolation_weights(t_resc,weights(0:n))
    do i = icpl_lo, icpl_hi  !1, cpl%def%n_couplings
       alpha(i) = dot_product(weights(0:n),subseg%ra(i,it_lo:it_hi))
       !call polint(interp_resc_t(0:n),subseg%ra(i,it_lo:it_hi),t_resc,
       !alpha(i),dalpha)
    end do
    
    ! remember that not all couplings are stored as reciprocal
    cpl_def_glb => cpl%def
    alpha = reciprocal_or_not(alpha)
  end function segment_alpha_sctn
  


  !******************* extracing numbers of flav ************************
  !-------------------------------------------------------------
  !! Returns the number of flavours relevant for scale Q
  !! Optionally, returns also the range of Qvalues (Qlo, Qhi) for which 
  !! this value of nf remains the same.
  !!
  !! Use convention that AT the quark mass that flavor is considered active
  integer function nf_at_Q(cpl, Q, Q_lo, Q_hi, muMatch_M)
    type(Coupling),  target, intent(in)  :: cpl
    real(dp),        intent(in)  :: Q
    real(dp),        intent(out), optional :: Q_lo, Q_hi
    real(dp),        intent(in),  optional :: muMatch_M
    nf_at_Q = nf_at_t(cpl, tofQ(Q), Q_lo, Q_hi, muMatch_M)
  end function nf_at_Q
  

  !-------------------------------------------------------------
  !! Returns the number of flavours relevant for scale Q
  !! Optionally, returns also the range of Qvalues (Qlo, Qhi) for which 
  !! this value of nf remains the same.
  !!
  !! Use convention that AT the quark mass that flavor is considered active
  integer function nf_at_t(cpl, t, Q_lo, Q_hi, muMatch_M)
    type(Coupling),  target, intent(in)  :: cpl
    real(dp),        intent(in)  :: t
    real(dp),        intent(out), optional :: Q_lo, Q_hi
    real(dp),        intent(in),  optional :: muMatch_M
    !-----------------------------
    type(CouplingDefinition), pointer :: def
    real(dp) :: deltat_match
    integer  :: nf
    real(dp), pointer :: mod_bound(:)
    real(dp), pointer :: inr_bound(:)

    def => cpl%def

    if (present(muMatch_M)) then
       allocate(mod_bound(cpl%def%nf_lo:cpl%def%nf_hi+1))
       deltat_match = tofQ(default_or_opt(&
            &def%muMatch_M, muMatch_M)/def%muMatch_M)

       ! redefine limits so as to account for requested matching
       ! thresholds muMatch_M. Outer limits should not be modified
       ! since they are rigidly fixed?
       mod_bound = cpl%segment_boundaries
       inr_bound => mod_bound(def%nf_lo+1:def%nf_hi)
       inr_bound = inr_bound + deltat_match
       inr_bound = min(cpl%segment_boundaries(def%nf_hi+1),&
            &          max(cpl%segment_boundaries(def%nf_lo),inr_bound))
    else
       mod_bound => cpl%segment_boundaries
    end if
    
    do nf = def%nf_lo, def%nf_hi
       if (t >= mod_bound(nf) .and. t < mod_bound(nf+1)) exit
    end do
    
    if (nf > def%nf_hi) then
       call wae_error('nf_at_t: Specified t is not in supported range',dbleval=t)
    end if

    if (present(Q_lo)) Q_lo = Qoft(mod_bound(nf))
    if (present(Q_hi)) Q_hi = Qoft(mod_bound(nf+1))

    if (present(muMatch_M)) deallocate(mod_bound)

    nf_at_t = nf
  end function nf_at_t
  

  !======================================================================
  !! returns values for nf_lo and nf_hi for this coupling
  subroutine nf_range(cpl, nf_lo, nf_hi)
    type(Coupling), intent(in)  :: cpl
    integer,        intent(out) :: nf_lo, nf_hi
    nf_lo = cpl%def%nf_lo
    nf_hi = cpl%def%nf_hi
  end subroutine nf_range
  
  !-------------------------------------------------------------
  ! returns the Q range for a given value of nf. If supported.
  subroutine t_Range_at_nf(cpl, nf, t_lo, t_hi, muMatch_M)
    type(Coupling),  intent(in)  :: cpl
    integer,         intent(in)  :: nf
    real(dp),        intent(out) :: t_lo, t_hi
    real(dp), optional, intent(in) :: muMatch_M
    !---------------------------------------
    real(dp) :: deltat_match
    character(len=60) :: string

    if (nf < cpl%def%nf_lo .or. nf > cpl%def%nf_hi) then
       write(string,'(a,i0,a)') 'nf value ',nf,' not supported'
       call wae_Error('t_Range_at_nf', trim(string))
    end if

    deltat_match = tofQ(default_or_opt(&
         &cpl%def%muMatch_M, muMatch_M)/cpl%def%muMatch_M)

    ! remember that only internal boundaries get modified 
    ! by deltat_match
    t_lo = cpl%segment_boundaries(nf)
    if (nf > cpl%def%nf_lo) t_lo = t_lo + deltat_match
    t_lo = max(cpl%t_lo,min(cpl%t_hi,t_lo))

    t_hi = cpl%segment_boundaries(nf+1)
    if (nf < cpl%def%nf_hi) t_hi = t_hi + deltat_match
    t_hi = max(cpl%t_lo,min(cpl%t_hi,t_hi))
 
  end subroutine t_Range_at_nf


  !-------------------------------------------------------------
  ! returns the Q range for a given value of nf. If supported.
  subroutine Q_Range_at_nf(cpl, nf, Q_lo, Q_hi, muMatch_M)
    type(Coupling),  intent(in)  :: cpl
    integer,         intent(in)  :: nf
    real(dp),        intent(out) :: Q_lo, Q_hi
    real(dp), optional, intent(in) :: muMatch_M
    !---------------------------------------
    real(dp) :: t_lo, t_hi

    call t_Range_at_nf(cpl, nf, t_lo, t_hi, muMatch_M)
    Q_lo = Qoft(t_lo)
    Q_hi = Qoft(t_hi)
  end subroutine Q_Range_at_nf

  !======================================================================
  !! Return the result of the application to alpha of the operation
  !! implicitly contained in the supplied terms.
  function apply_operation_alpha(terms,alpha) result(dalpha)
    type(CouplingTerm), intent(in) :: terms(:)
    real(dp),           intent(in) :: alpha(:)
    real(dp)                       :: dalpha(size(alpha))
    integer  :: i, ia
    real(dp) :: term
    dalpha = zero
    do i = 1, size(terms)
       term = terms(i)%coeff
       if (term == zero) cycle
       do ia = 1, size(terms(i)%alpha_powers)
          if (terms(i)%alpha_powers(ia) /= 0) then
             term = term * alpha(ia)**terms(i)%alpha_powers(ia)
          end if
       end do
       dalpha(terms(i)%itarget) = dalpha(terms(i)%itarget) + term
    end do
  end function apply_operation_alpha
  

  !************************* ancillary stuff ******************************
  !! conversion from t to Q (helps avoid risk of forgetting factors of two)
  elemental function tOfQ(Q)
    real(dp), intent(in) :: Q
    real(dp)             :: tOfQ
    if (Q <= zero) then
       !write(0,*) "resetting Q", Q
       tofQ = -1e200_dp
    else
       tOfQ = two * log(Q)
    end if
  end function tOfQ
  !! conversion from Q to t (helps avoid risk of forgetting factors of two)
  elemental function QOft(t)
    real(dp), intent(in) :: t
    real(dp)             :: QOft
    QOft = exp(half*t)
  end function QOft

  !----------------------------------------------------------------------
  subroutine rkstp_1d_special(h,x,y, deriv)
    real(dp), intent(in)    :: h
    real(dp), intent(inout) :: x, y(:)
    interface
       !! function that returns either dalpha/dt (for nf_in=nf_out)
       !! or delta alpha when crossing a threshold from nf_in to nf_out
       function deriv(t, alpha, nf_in, nf_out) result(res)
         use types; implicit none
         real(dp), intent(in) :: t
         real(dp), intent(in) :: alpha(:)
         integer,  intent(in) :: nf_in, nf_out
         real(dp)             :: res(size(alpha))
       end function deriv
    end interface
    real(dp), parameter :: third = one/three
    !------------------------------------------------------------
    real(dp) :: w1(size(y)), w2(size(y)), w3(size(y))
    real(dp) :: hh
    hh = half * h
    call calc_dra_special(x,y,w1,deriv); w1 = w1 * hh            ! w1 = k1/2
    call calc_dra_special(x+hh, y + w1, w2,deriv); w2 = w2 * hh  ! w2 = k2/2
    call calc_dra_special(x+hh, y + w2, w3,deriv); w3 = w3 * h   ! w3 = k3
    w2 = w1 + two*w2                                       ! w2 = half*k1 + k2
    call calc_dra_special(x+h , y + w3, w1,deriv); w1 = w1 * hh  ! w1 = k4/2
    
    x = x + h
    !                k1/2 + k2  + k3 + k4/2
    y = y + third * (w2         + w3 + w1)
  end subroutine rkstp_1d_special


  
  !======================================================================
  !! returns the reciprocal of arg where cpl_def_glb%use_reciprocal = .true.
  !! otherwise returns arg itself
  function reciprocal_or_not(arg)
    real(dp), intent(in) :: arg(:)
    real(dp)             :: reciprocal_or_not(size(arg))
    where (cpl_def_glb%use_reciprocal)
       reciprocal_or_not = one/arg
    elsewhere
       reciprocal_or_not = arg
    end where
  end function reciprocal_or_not
  
  !======================================================================
  !! returns i unless it is outside the bounds imin imax, in which case
  !! it returns the closer bound, assuming imax > imin
  elemental integer function in_range_int(i, imin, imax) result(in_range)
    integer, intent(in) :: i, imin, imax
    in_range = max(imin, min(imax, i))
  end function in_range_int
  elemental real(dp) function in_range_dp(i, imin, imax)  result(in_range)
    real(dp), intent(in) :: i, imin, imax
    in_range = max(imin, min(imax, i))
  end function in_range_dp

end module couplings
