!---------------------------------------------------------------
!! provides facilities related to alphas calculations
!!
!! there are two ways of accessing it: either with a handle which 
!! contains the details of the particular alfas which is being dealt
!! with (e.g. one loop etc...), or without a handle, in which case the
!! global handle is used.
!!
!! Original core based on BRW (analytical NLO running coupling)
!! Tested 25/04/00
!!
!! NB: Subsequently updated to makes use of the new_as module which
!! implements running coupling via the solution of the QCD
!! differential equation. Therefore the code present in this module is
!! mostly redundant and the real work is being done in new_as.f90
!!  
!-------------------------------------------------------------------
module qcd_coupling
  use new_as; use types; use consts_dp
  use warnings_and_errors
  implicit none
  private


  type running_coupling 
     private
     type(na_handle) :: nah      ! simply "redirect" to new coupling
     logical         :: use_nah
     !-- following components are related to old analytical approximations
     real(dp) :: QCDL5
     real(dp) :: BN(3:6), CN(3:6), CN5(3:6)
     integer  :: nloop, nf
  end type running_coupling
  public :: running_coupling

  interface InitRunningCoupling
     module procedure  as_Init_ash, as_Init_noash
  end interface
  public :: InitRunningCoupling
  interface Value
     module procedure  as_Value_ash, as_Value_noash
  end interface
  public :: Value
  interface Delete
     !-- for digital f90 the second one causes problems
     module procedure  as_Del_ash!, as_Del_noash
  end interface
  public :: Delete
  public :: NfRange, NfAtQ, QRangeAtNf, QuarkMass
  interface NumberOfLoops
     module procedure NumLoops_in_coupling
  end interface
  public :: NumberOfLoops
  public :: SetDefaultCouplingDt, DefaultCouplingDt

  ! A public instance of the coupling (as far as an outside user is
  ! concerned). This is set if the user calls the initialisation
  ! routine without providing a "running_coupling" object as an 
  ! argument. HIGHLY discouraged and may break in the future...
  type(running_coupling), target, save :: ash_global
  public :: ash_global


  !--- all constants below relate to the now somewhat obsolete
  !    NLO analytical approximation for the running coupling
  !-----------------------------------------------------------------
  !-- inclusion of CN is a rather tricky way of introducing 
  !   support for 1 and 2-loops in a fairly similary way
  !   and B(3:6) is a way of introducing support for fixed or variable
  !   flavour numbers
  real(dp), parameter :: rmass(6) = (/0D0,0D0,.15D0,1.7D0,5.3D0,175D0/)
  !real(dp), parameter :: rmass(6) = (/0D0,0D0,1D0,1.1D0,1.2D0,2005D0/)
  real(dp), parameter :: CAFAC = 3.0_dp, CFFAC = 4.0_dp/3.0_dp
  real(dp), parameter :: PIFAC = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: B3=((11.0_dp*CAFAC)- 6.0_dp)/(12.0_dp*PIFAC) 
  real(dp), parameter :: B4=((11.0_dp*CAFAC)- 8.0_dp)/(12.0_dp*PIFAC) 
  real(dp), parameter :: B5=((11.0_dp*CAFAC)-10.0_dp)/(12.0_dp*PIFAC) 
  real(dp), parameter :: B6=((11.0_dp*CAFAC)-12.0_dp)/(12.0_dp*PIFAC) 
  !-- with absoft there are some difficulties in accessing this quant
  !   from outside 
  real(dp), parameter, public :: as_BN(3:6) = (/ B3, B4, B5, B6/)
  real(dp), parameter :: trueC3=((17.0_dp*CAFAC**2)-(5.0_dp*CAFAC+3.0_dp&
       &*CFFAC)*3.0_dp) /(24.0_dp*PIFAC**2)/B3**2 
  real(dp), parameter :: trueC4=((17.0_dp*CAFAC**2)-(5.0_dp*CAFAC+3.0_dp&
       &*CFFAC)*4.0_dp) /(24.0_dp*PIFAC**2)/B4**2 
  real(dp), parameter :: trueC5=((17.0_dp*CAFAC**2)-(5.0_dp*CAFAC+3.0_dp&
       &*CFFAC)*5.0_dp) /(24.0_dp*PIFAC**2)/B5**2 
  real(dp), parameter :: trueC6=((17.0_dp*CAFAC**2)-(5.0_dp*CAFAC+3.0_dp&
       &*CFFAC)*6.0_dp) /(24.0_dp*PIFAC**2)/B6**2 
  real(dp), parameter :: as_CN(3:6) = (/ trueC3, trueC4, trueC5, trueC6 /)


contains

  !-- alternative versions to allow access to a global alpha
  subroutine as_Init_noash(alfas, Q, qcdl5, nloop, fixnf,  &
       & quark_masses, muMatch_mQuark, use_nah)
    real(dp), intent(in), optional :: alfas, Q, qcdl5
    integer,  intent(in), optional :: nloop, fixnf
    real(dp), intent(in), optional :: quark_masses(4:6)
    real(dp), intent(in), optional :: muMatch_mQuark
    logical,  intent(in), optional :: use_nah
    type(running_coupling), pointer :: coupling
    coupling => ash_global
    call as_Init_ash(coupling, alfas, Q, nloop, fixnf,  &
       & quark_masses, muMatch_mQuark, use_nah, qcdl5)
  end subroutine as_Init_noash
  
  function as_Value_noash(Q,fixnf) result(Value)
    real(dp), intent(in) :: Q
    real(dp)             :: Value
    type(running_coupling), pointer  :: coupling
    integer, intent(in), optional :: fixnf
    coupling => ash_global
    Value = as_Value_ash(coupling, Q, fixnf)
  end function as_Value_noash

  subroutine as_Del_noash()
    type(running_coupling), pointer :: coupling
    coupling => ash_global
    call as_Del_ash(coupling)
  end subroutine as_Del_noash
  

  !======================================================================
  ! The following are essentially interface routines to bryan's 
  ! code for alfas
  ! if present alfas then set coupling such alfas at M_Z or Q is equal to alpha_s
  ! Alternatively set qcdl5
  ! Can also optionally modify the number of loops and fix the number
  ! of flavours
  !
  ! As of 30/06/2001, by sepcifying use_nah=.true. it should give
  ! fairly transparent access to the differential equation solution
  ! for alfas, using new_as. The only proviso is that initialisation 
  ! may well be quite a bit slower. Also if we want to reinitialise
  ! then it is necessary to call Delete
  !
  ! As od 16/01/2001 use_nah=.true. will become the default
  subroutine as_Init_ash(coupling, alfas, Q, nloop, fixnf, &
       & quark_masses, muMatch_mQuark, use_nah, qcdl5)
    use assertions
    type(running_coupling) :: coupling
    real(dp), intent(in), optional :: alfas, Q, qcdl5
    integer,  intent(in), optional :: nloop, fixnf
    real(dp), intent(in), optional :: quark_masses(4:6)
    real(dp), intent(in), optional :: muMatch_mQuark
    logical,  intent(in), optional :: use_nah
    real(dp) :: dummy, lower, upper, middle
    real(dp) :: Qloc
    real(dp), parameter :: mz = 91.2_dp, eps = 1e-8_dp

    coupling%use_nah = default_or_opt(.true., use_nah)
    if (coupling%use_nah) then
       call na_Init(coupling%nah, alfas, Q, nloop, fixnf, &
            &           quark_masses, muMatch_mQuark)
       return
    else
       call wae_error('as_Init_ash','old alphas (analytic approx to solution) is obsolete and will not run...')
    end if

    if (present(nloop)) then
       coupling%nloop = nloop
    else
       coupling%nloop = 2
    end if

    !-- a fudge to fix the number of flavours
    if (present(fixnf)) then
       coupling%BN = as_BN(fixnf)
       coupling%CN = as_CN(fixnf)
    else
       coupling%BN = as_BN
       coupling%CN = as_CN
    end if
    
    !-- a fudge to set the number of loops
    select case (coupling%nloop)
    case(1)
       coupling%CN   = zero
    case(2)
       ! keep coupling%CN = as_CN as it was set before
    case default
       stop 'InitRunningCoupling: nloop must be either 1 or 2'
    end select

    if (present(alfas)) then
       if (present(Q)) then
          Qloc = Q
       else
          Qloc = mz
       end if
       !-- set up search limits ----------------------------
       lower = 0.01_dp
       upper = 1.0_dp
       !-- we may need to revise lower and upper...
       do
          call priv_as_init_2l(coupling,lower); dummy = Value(coupling,Qloc)
          if (dummy > alfas) then
             !call hwwarn(100)
             upper = lower
             lower = lower**2
             cycle
          end if
          call priv_as_init_2l(coupling,upper); dummy = Value(coupling,Qloc)
          if (dummy < alfas) call hwwarn(101)
          exit
       end do
       do 
          !-- could be more efficient -- but too lazy for now
          !write(0,'(4f20.15)') upper, lower, dummy, alfas
          middle = sqrt(lower*upper)
          call priv_as_init_2l(coupling,middle)
          dummy = Value(coupling,Qloc)
          if (dummy >= alfas) upper = middle
          if (dummy < alfas) lower = middle
          if (upper/lower-1.0_dp < eps ) exit
       end do
    elseif (present(qcdl5)) then
       call priv_as_init_2l(coupling,qcdl5)
    else
       !call priv_as_init_2l(coupling,0.280_dp)  ! alfas(91) = 0.124
       call priv_as_init_2l(coupling,0.214_dp)  ! alfas(91) = 0.117
       !call priv_as_init_2l(coupling,0.158_dp)  ! alfas(91) = 0.112
    end if
  end subroutine as_Init_ash


  !======================================================================
  function as_Value_ash(coupling, Q, fixnf) result(Value)
    real(dp), intent(in) :: Q
    real(dp)             :: Value
    type(running_coupling)      :: coupling
    integer, intent(in), optional :: fixnf
    if (coupling%use_nah) then
       Value = na_Value(coupling%nah, Q, fixnf)
    else
       if (present(fixnf)) call wae_error('as_Value_ash: &
            &fixnf not support for old alpha_s')
       Value = hwualf(coupling, 1, Q)
    end if
  end function as_Value_ash


  !======================================================================
  !! Returns the number of loops used for the evolution of this coupling.
  integer function NumLoops_in_coupling(coupling)
    type(running_coupling) :: coupling
    if (coupling%use_nah) then
       NumLoops_in_coupling = na_NumberOfLoops(coupling%nah)
    else
       NumLoops_in_coupling = coupling%nloop
    end if
    
  end function NumLoops_in_coupling
  

  !-----------------------------------------------------------------
  ! do any required cleaning up 
  subroutine as_Del_ash(coupling)
    type(running_coupling)      :: coupling
    if (coupling%use_nah) call na_Del(coupling%nah)
  end subroutine as_Del_ash

  !-----------------------------------------------------------------
  ! Indicate the range of nf values supported by this handle.
  subroutine NfRange(coupling, nflo, nfhi)
    type(running_coupling), intent(in)  :: coupling
    integer,         intent(out) :: nflo, nfhi
    if (coupling%use_nah) then
       call na_nfRange(coupling%nah,nflo,nfhi)
    else
       call wae_Error('NfRange: this routine&
            & is only supported with new alpha_s')
    end if
  end subroutine NfRange
  
  !-------------------------------------------------------------
  ! Returns the number of flavours relevant for scale Q
  ! Also returns the range of Qvalues (Qlo, Qhi) for which 
  ! this value of nf remains the same.
  function NfAtQ(coupling, Q, Qlo, Qhi, muM_mQ)
    integer                             :: nfAtQ
    type(running_coupling), intent(in)  :: coupling
    real(dp),        intent(in)  :: Q
    real(dp),        intent(out), optional :: Qlo, Qhi
    real(dp),        intent(in),  optional :: muM_mQ
    !-----------------------------
    if (coupling%use_nah) then
       call na_nfAtQ(coupling%nah, Q, nfAtQ, Qlo, Qhi, muM_mQ)
    else
       call wae_Error('NfAtQ: this routine&
            & is only supported with new alpha_s')
    end if
  end function NfAtQ
  
  
  !======================================================================
  real(dp) function QuarkMass(coupling, iflv) result(mass)
    type(running_coupling), intent(in) :: coupling
    integer,                intent(in) :: iflv
    if (coupling%use_nah) then
       mass = na_QuarkMass(coupling%nah, iflv)
    else
       call wae_Error('QuarkMass: this routine&
            & is only supported with new alpha_s')
    end if
  end function QuarkMass



  !-------------------------------------------------------------
  ! returns the Q range for a given value of nf. If supported.
  subroutine QRangeAtNf(coupling, nflcl, Qlo, Qhi, muM_mQ)
    type(running_coupling), intent(in)  :: coupling
    integer,         intent(in)  :: nflcl
    real(dp),        intent(out) :: Qlo, Qhi
    real(dp),        intent(in),  optional :: muM_mQ
    if (coupling%use_nah) then
       call na_QrangeAtNf(coupling%nah, nflcl, Qlo, Qhi, muM_mQ)
    else
       call wae_Error('QRangeAtNf: this routine&
            & is only supported with new alpha_s')
    end if
  end subroutine QRangeAtNf
  
  !----------------------------------------------------------------------
  ! the routine that does the real initialisation of alfas
  subroutine priv_as_init_2l(coupling, qcdl5)
    type(running_coupling)      :: coupling
    real(dp), intent(in) :: qcdl5
    real(dp) :: dummy
    coupling%qcdl5 = qcdl5
    dummy = HWUALF(coupling, 0, 1e2_dp)
  end subroutine priv_as_init_2l
  !-----------------------------------------------------------------------
  !     STRONG COUPLING CONSTANT
  !     IOPT == 0  INITIALIZES
  !          == 1  TWO-LOOP, FLAVOUR THRESHOLDS
  !
  ! BRW routine, with various bits and pieces shuffled around or
  ! removed by GPS on 24/04/00
  !---------------------------------------------------------------------
  FUNCTION HWUALF(coupling, IOPT, SCALE) 
    REAL(DP) :: HWUALF 
    type(running_coupling), target :: coupling
    real(dp), intent(in) :: scale
    integer,  intent(in) :: iopt
    !--------------------------------------------------------------------
    REAL(DP)          :: RHO,RAT,RLF
    real(dp), pointer :: CN5(:), CN(:), QCDL5, BN(:)
    integer           :: nf

    CN5   => coupling%CN5
    QCDL5 => coupling%QCDL5
    CN    => coupling%CN
    BN    => coupling%BN

    !-- the rest is init related to this value of lambda5
    IF (IOPT == 0) THEN 
       !---QCDL5 IS 5-FLAVOUR LAMBDA-MS-BAR                                    
       !---COMPUTE THRESHOLD MATCHING                                          
       RHO=2.0_dp*LOG(RMASS(6)/QCDL5) 
       RAT=LOG(RHO)/RHO 
       CN5(6)=(coupling%BN(5)/(1.0_dp-CN(5)*RAT)-coupling%BN(6)/(1.0_dp-CN(6)*RAT))*RHO 
       RHO=2.0_dp*LOG(RMASS(5)/QCDL5) 
       RAT=LOG(RHO)/RHO 
       CN5(4)=(coupling%BN(5)/(1.0_dp-CN(5)*RAT)-coupling%BN(4)/(1.0_dp-CN(4)*RAT))*RHO 
       RHO=2.0_dp*LOG(RMASS(4)/QCDL5) 
       RAT=LOG(RHO)/RHO 
       CN5(3)=(coupling%BN(4)/(1.0_dp-CN(4)*RAT)-&
            & coupling%BN(3)/(1.0_dp-CN(3)*RAT))*RHO+CN5(4) 
       CN5(5) = zero
    ENDIF
    IF (SCALE <= QCDL5) then; CALL HWWARN(51); hwualf=0.0_dp; return
    end IF
    RHO=2.0_dp*LOG(SCALE/QCDL5) 
    RAT=LOG(RHO)/RHO 
    !-- this will allow us to fiddle nf later on, according to some
    !   flag in coupling
    nf = priv_nf(coupling, scale)
    RLF = coupling%BN(nf)*RHO/(one - CN(nf)*RAT) + CN5(nf)
!!$    select case(nf)
!!$    case(6)
!!$       RLF=B6*RHO/(1.0_dp-CN(6)*RAT)+CN5(6) 
!!$    case(5)
!!$       RLF=B5*RHO/(1.0_dp-CN(5)*RAT) 
!!$    case(4)
!!$       RLF=B4*RHO/(1.0_dp-CN(4)*RAT)+CN5(4) 
!!$    case default
!!$       RLF=B3*RHO/(1.0_dp-CN(3)*RAT)+CN5(3) 
!!$    end select
    IF (RLF <= ZERO) then; CALL HWWARN(53); hwualf=0.0_dp; return
    end IF

    IF (IOPT == 1) THEN 
       HWUALF=1.0_dp/RLF 
    ENDIF
  end function hwualf

  !---------------------------------------------------------------
  ! returns the appropriate for the given scale
  function priv_nf(coupling, scale) result(nf)
    type(running_coupling), intent(in) :: coupling
    real(dp),        intent(in) :: scale
    integer                     :: nf
    IF (SCALE > RMASS(6)) THEN 
       nf = 6
    ELSEIF (SCALE > RMASS(5)) THEN 
       nf = 5
    ELSEIF (SCALE > RMASS(4)) THEN 
       nf = 4
    ELSE 
       nf = 3
    ENDIF
  end function priv_nf
  

  !-- a little primitive? --------------------------------------------
  SUBROUTINE HWWARN(ICODE) 
    integer, intent(in) :: icode
    real(dp) :: a, b
    write(0,*) ' ALFAS WARNING CODE',ICODE
    a = sqrt(0.0_dp)
    b = a**2
    !write(0,*) 1.0_dp/(a-b)
    stop
  END SUBROUTINE HWWARN


end module qcd_coupling

!!$program astest
!!$  use types; use consts_dp
!!$  use qcd_coupling
!!$  implicit none
!!$  integer  :: i
!!$  real(dp) :: Q
!!$
!!$  call InitRunningCoupling(alfas=0.118_dp,nloop=2)
!!$  do i = 0,100
!!$     Q = 91.2_dp**(i/100.0_dp)
!!$     write(6,*) Q, Value(Q)
!!$  end do
!!$  
!!$end program astest
