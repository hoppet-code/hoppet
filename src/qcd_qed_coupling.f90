module qcd_qed_coupling
  use types
  use qcd
  use couplings
  use qed_coupling_module ! gives us lepton masses
  implicit none
  private

  integer, parameter, public :: qcd_mass_scheme_pole  = 1
  integer, parameter, public :: qcd_mass_scheme_MSbar = 2
  
  type QCDQEDThreshold
     real(dp) :: mu
     integer  :: pdg_id
     integer  :: nflav(3) = 0 ! (nleptons, ndown, nup)
     real(dp) :: qcd_mult = one
  end type QCDQEDThreshold
  public :: QCDQEDThreshold

  integer, parameter    :: n_thresholds = 8
  type QCDQEDCoupling
     type(Coupling)        :: Coupling
     type(QCDQEDThreshold) :: thresholds(0:n_thresholds)
     integer               :: quark_schemes(4:6)
  end type QCDQEDCoupling
  public :: QCDQEDCoupling

  public :: new_QCDQEDCoupling
  public :: coupling_value, coupling_array ! from couplings
  public :: coupling
contains

  ! what thresholds will we have?
  ! - 3 lepton masses
  ! - 3 heavy-quark masses (c,b,t)
  ! - effective uds mass for QED part
  ! - freezing scale for QCD coupling
  !
  ! loops is a string variable of the form "01 02 03 10"
  ! to indicate QCD running up to three loops and QED running up to one loop
  ! add "11" to get the mixed QCD-QED running terms
  function new_QCDQEDCoupling(alphas, alpha, mu, loops, quark_masses, quark_schemes) result(cpl)
    use sort; use assertions; use warnings_and_errors
    real(dp), intent(in) :: alphas, alpha, mu
    character(len=*), intent(in) :: loops
    real(dp), intent(in) :: quark_masses(4:6)
    integer,  intent(in) :: quark_schemes(4:6)
    type(QCDQEDCoupling), target :: cpl
    !----------------------------------------------------------------------
    type(QCDQEDThreshold), target  :: thresholds(0:n_thresholds)
    type(QCDQEDThreshold), pointer :: threshold, last
    integer  :: i, indices(0:n_thresholds)
    integer  :: nloops, nloops_tot, loops_int(20)=-1
    logical :: qcd_target, qed_target
    type(CouplingTerm), allocatable :: beta_terms(:,:), thr_terms(:,:), rthr_terms(:,:)
    type(CouplingDefinition) :: coupling_def
    integer :: ithr, iloop, nf_lcl, nf_old, nf_diff, nlept
    
    ! first we collect the thresholds
    thresholds(0) = QCDQEDThreshold(m_electron/10   , -1 ) ! lower edge
    thresholds(1) = QCDQEDThreshold(m_electron      , 11 )
    thresholds(2) = QCDQEDThreshold(m_muon          , 13 )
    thresholds(3) = QCDQEDThreshold(m_tau           , 15 )
    thresholds(4) = QCDQEDThreshold(0.5_dp          , 123) ! turn off quarks in QED coupling
    thresholds(5) = QCDQEDThreshold(1.0_dp          , 21 ) ! unfreeze QCD coupling
    thresholds(6) = QCDQEDThreshold(quark_masses(4) , 4  )
    thresholds(7) = QCDQEDThreshold(quark_masses(5) , 5  )
    thresholds(8) = QCDQEDThreshold(quark_masses(6) , 6  )

    ! then we sort them (and adjust output, which is based on numbers from 1 upwards)
    call indexx(thresholds%mu, indices); indices = indices - 1

    ! just check that lower bound is genuinely zero
    i = assert_eq(indices(0),0, 'new_QCDQEDCoupling internal inconsistent lower bound')  

    ! now we establish the number of each kind of flavour
    do i = 1, n_thresholds
       last      => thresholds(indices(i-1))
       threshold => thresholds(indices(i))
       if (123 == threshold%pdg_id) then
          ! handle light quarks
          threshold%nflav(:) = last%nflav(:) + (/ 0, 2, 1 /)
       else  if (1 <= threshold%pdg_id .and. threshold%pdg_id <= 6) then
          ! handle heavy quarks
          if (mod(threshold%pdg_id,2) == 1) then
             threshold%nflav(:) = last%nflav(:) + (/ 0, 1, 0 /)
          else
             threshold%nflav(:) = last%nflav(:) + (/ 0, 0, 1 /)
          end if
       else if (11 <= threshold%pdg_id .and. threshold%pdg_id <= 15) then
          ! handle leptons
          threshold%nflav(:) = last%nflav(:) + (/1, 0, 0/)
       else if (21 == threshold%pdg_id) then
          ! no changes in flavour, but everything below this
          ! has frozen QCD coupling, obtained by setting qcd_mult -> 0
          threshold%nflav(:) = last%nflav(:)
          thresholds(indices(:i-1))%qcd_mult = zero
       else
          call wae_error('new_QCDQEDCoupling','inconsistent pdg_id', intval=threshold%pdg_id)
       end if
    end do
    cpl%thresholds = thresholds(indices(:))

    ! next we figure out which loop terms we will need
    nloops = Count_Items(loops)
    read(loops,*) loops_int(1:nloops)
    nloops_tot = nloops
    do i = 1, nloops
       qcd_target = (digit(1,loops_int(i)) /= 0)
       qed_target = (digit(2,loops_int(i)) /= 0)
       if (qcd_target) then
          ! this multiplies QCD coupling (index 1)
          loops_int(i) = loops_int(i) + 100             ! target is QCD
          if (qed_target) then
             ! and also QED coupling, so we add an extra entry
             nloops_tot = nloops_tot + 1
             loops_int(nloops_tot) = loops_int(i) + 100 ! target for this entry is QED
          end if
       else if (qed_target) then
          loops_int(i) = loops_int(i) + 200             ! target is QED
       else
          call wae_error('new_QCDQEDCoupling','inconsistent loops entry', intval=loops_int(i))
       end if
       
    end do
    nloops = nloops_tot

    ! allocate the terms we need
    allocate(beta_terms(1:nloops, 0:n_thresholds))
    allocate(thr_terms (1:nloops, 0:n_thresholds))
    allocate(rthr_terms(1:nloops, 0:n_thresholds))
    ! set defaults
    beta_terms = new_CouplingTerm(1, (/0,0/), zero)
    thr_terms  = new_CouplingTerm(1, (/0,0/), zero)
    rthr_terms = new_CouplingTerm(1, (/0,0/), zero)
    
    write(6,*) loops_int(1:nloops)

    ! now fill up 
    do ithr = 0, n_thresholds
       threshold => cpl%thresholds(ithr)
       nf_old = sum(cpl%thresholds(max(0,ithr-1))%nflav(2:3))
       nf_lcl = sum(threshold%nflav(2:3))
       nf_diff = nf_lcl - nf_old
       call qcd_SetNf(nf_lcl)
       nlept  = threshold%nflav(1)
       do iloop = 1, nloops
          select case(loops_int(iloop))
          case(101)
             beta_terms(iloop,ithr) = new_CouplingTerm(1, (/1,0/), -beta0)
          case(102)
             beta_terms(iloop,ithr) = new_CouplingTerm(1, (/2,0/), -beta1)
          case(103)
             beta_terms(iloop,ithr) = new_CouplingTerm(1, (/3,0/), -beta2)
             if (nf_diff == 1) thr_terms (iloop,ithr) = new_CouplingTerm(1, (/2,0/),  alphastep20_pole)
             if (nf_diff == 1) rthr_terms(iloop,ithr) = new_CouplingTerm(1, (/2,0/), -alphastep20_pole)
          case(104)
             beta_terms(iloop,ithr) = new_CouplingTerm(1, (/4,0/), -beta3)
             ! note unconventional sign here, as used in new_as.f90
             if (nf_diff == 1) thr_terms (iloop,ithr) = new_CouplingTerm(1, (/3,0/), &
                  &                       -alphastep30_pole-(nf_lcl-1)*alphastep30_pole_nl)
             if (nf_diff == 1) rthr_terms(iloop,ithr) = new_CouplingTerm(1, (/3,0/),  &
                  &                        alphastep30_pole+(nf_lcl-1)*alphastep30_pole_nl)
          case(210)
             beta_terms(iloop,ithr) = new_CouplingTerm(2, (/0,1/), &
                  &two/(6*pi)*qed_squared_charge(threshold%nflav(1),&
            &                                    threshold%nflav(2),&
            &                                    threshold%nflav(3)))
          case(211)
             ! got this from MNSZ backup notes, Eq. 96
             beta_terms(iloop,ithr) = new_CouplingTerm(2, (/1,1/), &
                  & three*CF/(4*pi)*one/(three*pi)*&
                  &           qed_squared_charge(                 0,&
            &                                    threshold%nflav(2),&
            &                                    threshold%nflav(3)))
          case(111)
             ! deduced this from Eq.(30) in 1512.00612
             beta_terms(iloop,ithr) = new_CouplingTerm(1, (/1,1/), &
                  & one/(four*pisq) * TR/CA*&
                  &           qed_squared_charge(                 0,&
            &                                    threshold%nflav(2),&
            &                                    threshold%nflav(3)))
             write(0,*) 'CHECK 111 case'
          case default
             call wae_error('new_QCDQEDCoupling','loops entry cannot be handled',&
                  &intval=loops_int(iloop))
          end select

          ! now post-process to make sure it makes sense
          call check_and_amend(beta_terms(iloop,ithr), loops_int(iloop), .false., cpl%thresholds(ithr         )%qcd_mult, 'beta')
          call check_and_amend(thr_terms (iloop,ithr), loops_int(iloop), .true. , cpl%thresholds(max(0,ithr-1))%qcd_mult, 'thr ')
          call check_and_amend(rthr_terms(iloop,ithr), loops_int(iloop), .true. , cpl%thresholds(max(0,ithr-1))%qcd_mult, 'rthr')

          ! write(6,*) ithr, iloop, loops_int(iloop)
          ! write(6,*) 'beta = ', beta_terms(iloop,ithr)%itarget, beta_terms(iloop,ithr)%alpha_powers, beta_terms(iloop,ithr)%coeff
          ! write(6,*) 'thr  = ', thr_terms (iloop,ithr)%itarget, thr_terms (iloop,ithr)%alpha_powers, thr_terms (iloop,ithr)%coeff
          ! write(6,*) 'rthr = ', rthr_terms(iloop,ithr)%itarget, rthr_terms(iloop,ithr)%alpha_powers, rthr_terms(iloop,ithr)%coeff
       end do
    end do
    
    coupling_def = new_CouplingDefinition(&
         & n_couplings = 2,&
         & nf_lo       = 0,&
         & nf_hi       = n_thresholds,&
         & beta_terms  = beta_terms,&
         & threshold_terms         = thr_terms (:,1:),&
         & reverse_threshold_terms = rthr_terms(:,1:),&
         & muMatch_M   = one)

    write(0,*) one/alpha
    cpl%Coupling = new_Coupling(coupling_def, (/alphas, alpha /), mu, cpl%thresholds(1:)%mu,&
         &                      Q_lo = m_electron/10.0_dp, Q_hi = 1e20_dp, eps=1e-10_dp)
  end function new_QCDQEDCoupling

  !----------------------------------------------------------------------
  ! check that the coupling term corresponds to the loops_int
  ! expectation (unless the coeff is zero) and adjust the power of the
  ! coupling for the target coupling
  subroutine check_and_amend(coupling_term, loops_int, is_threshold, qcd_mult, term)
    use warnings_and_errors
    type(CouplingTerm), intent(inout) :: coupling_term
    integer,            intent(in)    :: loops_int
    logical,            intent(in)    :: is_threshold
    real(dp),           intent(in)    :: qcd_mult
    character(len=*),   intent(in)    :: term
    !-------------------------------
    logical :: issue
    integer :: itarget
    integer :: adj_powers(lbound(coupling_term%alpha_powers,1):ubound(coupling_term%alpha_powers,1))

    itarget = coupling_term%itarget
    adj_powers = coupling_term%alpha_powers
    if (is_threshold) adj_powers(itarget) = adj_powers(itarget) + 1
    
    issue = .false.
    issue = issue .or. .not.( adj_powers(1)         == digit(1,loops_int) .or. coupling_term%coeff == zero)
    issue = issue .or. .not.( adj_powers(2)         == digit(2,loops_int) .or. coupling_term%coeff == zero)
    issue = issue .or. .not.( coupling_term%itarget == digit(3,loops_int) .or. coupling_term%coeff == zero)
    if (issue) call wae_error('new_QCDQEDCoupling->check_and_amend','issue with loop term '//term,intval=loops_int)

    ! add 1 to alpha_power corresponding to target, because above we've assumed that
    !
    !   dalpha(itarget) = alpha(itarget) * (...)
    !
    ! whereas the couplings module doesn't include the leading power
    coupling_term%alpha_powers(itarget) = coupling_term%alpha_powers(itarget) + 1

    ! then apply qcd_mult (typically either one or zero) to anything that involves
    ! a QCD coupling
    if (coupling_term%alpha_powers(1) /= 0 .or. coupling_term%itarget == 1)&
         &coupling_term%coeff = coupling_term%coeff * qcd_mult
  end subroutine check_and_amend
  
  
  ! --------------------
  ! from http://fortranwiki.org/fortran/show/String_Functions
  INTEGER FUNCTION Count_Items(s1)  ! in string or C string that are blank or comma separated
    CHARACTER(*) :: s1
    CHARACTER(len(s1)) :: s
    INTEGER :: i, k
    
    s = s1                            ! remove possible last char null
    k = 0  ; IF (s /= ' ') k = 1      ! string has at least 1 item
    DO i = 1,LEN_TRIM(s)-1
       IF (s(i:i) /= ' '.AND.s(i:i) /= ',' &
            .AND.s(i+1:i+1) == ' '.OR.s(i+1:i+1) == ',') k = k+1
    END DO
    Count_Items = k
  END FUNCTION Count_Items

  ! returns the i^th digit of number n
  integer function digit(i,n)
    integer, intent(in) :: i, n
    digit = mod(n/10**(i-1),10)
  end function digit
  
  
end module qcd_qed_coupling
