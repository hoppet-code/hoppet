! This file copyright the Hoppet authors, released under GPLv3

!======================================================================
! What follows is mainly interfaces for VFNS at order(as^3). The
! actual code is contained in VFNSmod.f, momnsmod.f, and AQg3mod.f90,
! which is based on code from Blumlein and collaborators.
module blumlein_interfaces
  use types; use consts_dp; use assertions; use vfnsmod; use aqg3mod;&
       & use momnsmod
  implicit none

  ! generic interface corresponding to functions of z in code from
  ! Johannes Blumlein, allowing us to easily have pointers to them
  abstract interface
    double precision function JBfunc(z,nf,a4pi,LL)
      double precision :: z,nf,a4pi,LL
    end function JBfunc
  end interface

  !abstract interface
  !  double precision function JBone(z)
  !    double precision, intent(in) :: z
  !  end function JBone
  !end interface

  ! generic interface corresponding to functions of N in code from
  ! Johannes Blumlein, allowing us to easily have pointers to them
  abstract interface
    double precision function JBmom(N,z,nf,a4pi,LL)
      integer :: N
      double precision :: z,nf,a4pi,LL
    end function JBmom
  end interface

  abstract interface
    double precision function JBMomInt(a1,a2,a3,a4,a5)
      integer :: a1,a2,a3,a4,a5
    end function JBMomInt
  end interface

contains
  !! a_Qq^(3,PS)
  double precision function APS(z,nf,as,LL)
    double precision :: z,nf,as,LL

    APS = APS1(z,nf,as,LL) + APS2(z,nf,as,LL)
  end function 

  double precision function AGGREG(z,nf,as,LL)
    double precision :: z,nf,as,LL

    IF(z.LT.0.5D0) AGGREG=AGGREG0(z,nf,as,LL)
    IF(z.GE.0.5D0) AGGREG=AGGREG1(z,nf,as,LL)
  end function 

  ! This is the plus part of the NS plus piece
  double precision function ANSPLU(z,nf,as,LL)
    double precision :: z,nf,as,LL

    ANSPLU = ANSPLU1(z,nf,as,LL) + ANSPLU2(z,nf,as,LL)
  end function

  ! This is the plus part of the light-flavour NS minus piece
  ! with arguments in our standard order
  double precision function OPLUS_znfasLL(z,nf,as,LL)  
    double precision :: z,nf,as,LL  
    OPLUS_znfasLL = OPLUS1(as, LL, NF, z) + OPLUS2(as, LL, NF, z)
  end function

  ! This is the regular part of the light-flavour NS minus piece
  ! with arguments in our standard order
  double precision function OREG_znfasLL(z,nf,as,LL)  
    double precision :: z,nf,as,LL  
    OREG_znfasLL = OREG(as, LL, NF, z)
  end function

  ! This is the delta-function part of the light-flavour NS minus piece
  ! with arguments in our standard order
  double precision function ODEL_znfasLL(z,nf,as,LL)  
    double precision :: z,nf,as,LL  
    ODEL_znfasLL = ODEL(as, LL, NF)
  end function

  !! 2024-03-22: an aproximation to AQG that gets the sum rule to within
  !! O(1) * (as/2pi)^3
  !!
  !! 2024-05-30: this is now the full AQG.
  double precision function AQGtotal(z,nf,as,LL)
    double precision :: z,nf,as,LL
    double precision :: RED

    IF(z.LT.0.5D0) RED=RED0(z)
    IF(z.GE.0.5D0) RED=RED1(z)
    
    AQGtotal = (aQg3(z)/2+RED+nf*AQG3NF(z))*as**3
    !AQGtotal = (RED+nf*AQG3NF(z))*as**3
  end function AQGtotal

end module blumlein_interfaces

module mass_thresholds_n3lo
  use types; use consts_dp; use assertions
  use qcd
  use convolution
  use convolution_communicator
  use blumlein_interfaces
  implicit none

  abstract interface
  real(dp) function HoppetPij(y)
    use types
    real(dp), intent(in) :: y
  end function HoppetPij
  end interface

  type comparator
    character(len=100) :: name = ""
    procedure(JBfunc),    nopass, pointer :: blum_fn   => null()
    procedure(HoppetPij), nopass, pointer :: hoppet_fn => null()
    integer            :: piece = -1
  contains
    procedure :: compare
  end type comparator



  type JB
    procedure(JBfunc), nopass, pointer :: reg    => null()
    procedure(JBfunc), nopass, pointer :: plus   => null()
    procedure(JBfunc), nopass, pointer :: delta  => null()
    !! - if order is 10n (n=2,3) then  JB_fn will return the coefficient
    !!   of (as/4pi)^n by evaluating the function at multiple alphas values
    !! - if order is 20n (n=2,3) then  JB_fn will return the coefficient
    !!   of (as/4pi)^n by evaluating the function at a single very small 
    !!   or large alphas value, such that the other order is relatively negligible
    !! - if order is n (n=2,3) then JB_fn will return the coefficient of
    !!   of (as/4pi)^n by evaluating the function at a single alphas value
    !! - care will be taken also with a 1*delta(1-x) component
    integer                    :: order = 3
  contains 
    !!
    !! a hoppet-style function that evalues the integrand needed for
    !! a grid_conv object (coefficient of as/2pi^n) taking y=ln(1/x)
    !! and returning xP(x). Call it as
    !!  
    !!    this%integrand(y)
    !!
    !! (it is normally accessed via the JB_ptr_integrand, which uses
    !! the module wide InitGridConv_JB_ptr)
    procedure :: integrand => JB_integrand
    !!
    !! underlying order separator function (deals with a4pi^n).
    !! To be called as 
    !! 
    !!       this%fn(this%reg, z)
    !!
    !! (or similarly with this%plus and this%delta) and it returns the
    !! coefficient of (as/4pi)^n, with n defined from order as above
    procedure :: fn => JB_fn
  end type JB

  interface InitGridConv
    module procedure InitGridConv_JB
  end interface
  interface AddWithCoeff
    module procedure AddWithCoeff_JB
  end interface
  class(JB), pointer :: InitGridConv_JB_ptr => null()

!  private
!  public init_MTM_N3LO
contains

  !! returns the integrand as needed by hoppet, assuming
  !! that InitGridConv_JB_ptr has been set up
  real(dp) function JB_ptr_integrand(y) result(res)
    real(dp), intent(in) :: y
    !-------------------------------
    res = InitGridConv_JB_ptr%integrand(y)
  end function JB_ptr_integrand

  !! performs InitGridConv with a JB object
  subroutine InitGridConv_JB(grid, gc, JB_obj) 
    type(grid_def),  intent(in)    :: grid
    type(grid_conv), intent(inout) :: gc
    class(JB), intent(in), target :: JB_obj
    !-------------------------------
    InitGridConv_JB_ptr => JB_obj
    call InitGridConv(grid, gc, JB_ptr_integrand)
  end subroutine InitGridConv_JB

  !! performs InitGridConv with a JB object
  subroutine AddWithCoeff_JB(gc, JB_obj) 
    type(grid_conv), intent(inout) :: gc
    class(JB), intent(in), target :: JB_obj
    !-------------------------------
    InitGridConv_JB_ptr => JB_obj
    call AddWithCoeff(gc, JB_ptr_integrand)
  end subroutine AddWithCoeff_JB

  !! returns the integrand as needed by hoppet
  !! Implements JB%integrand(...)
  real(dp) function JB_integrand(this, y) result(res)
    class(JB), intent(in) :: this
    real(dp),  intent(in) :: y
    !-------------------------------
    real(dp) :: z
    z = exp(-y)

    select case(cc_piece)
    case(cc_REAL)
      res = this%fn(this%reg,z)
      res = res + this%fn(this%plus,z)
    case(cc_REALVIRT)
      res = this%fn(this%reg,z)
    case(cc_VIRT)
      res = -this%fn(this%plus,z)
    case(cc_DELTA)
      res = this%fn(this%delta,z)
    end select
    if (cc_piece /= cc_DELTA) res = res * z
    ! convert from (as/4pi)**order normalisation to (as/2pi)**order
    res = res * half**mod(this%order,100)
  end function JB_integrand

  !! return the evaluation of the specific function at the point z
  !! working out the coefficient of as4pi**this%order
  !! Implements JB%fn(...)
  real(dp) function JB_fn(this, fn, z) result(res)
    use warnings_and_errors
    class(JB), intent(in) :: this
    procedure(JBfunc), pointer :: fn
    real(dp)                   :: z
    !-------------------------------
    integer  :: nf_light_int
    real(dp) :: orders(2:3), nf_light_dp
    real(dp) :: a4pi(0:2) =  (/ 0.0_dp, 0.10_dp, 1.0_dp /), this_a4pi
    real(dp) :: res_a4pi(0:2)
    logical  :: delta_one
    integer  :: i, imin, imax
    integer  :: mod_order

    ! diagnostics
    ! if (associated(fn,this%reg  )) write(6,*) 'JB_fn: reg'
    ! if (associated(fn,this%plus )) write(6,*) 'JB_fn: plus'
    ! if (associated(fn,this%delta)) write(6,*) 'JB_fn: delta'
    ! if (.not. associated(fn)     ) write(6,*) 'JB_fn: null'
    ! write(6,*) 'JB_fn: z = ', z

    res = zero

    ! if the function is a null ptr then return zero
    if (.not. associated(fn)) then
      return
    end if

    ! our convention is that nf_int = nf_light + 1
    nf_light_int = nf_int - 1
    nf_light_dp = real(nf_light_int, dp)

    ! otherwise, we need to extract the coefficient of 
    ! (as/4pi)**order from the function, with care for the
    ! the case where we are using the delta-function 
    ! component, which needs special treatment
    delta_one = associated(fn, this%delta)
    mod_order = mod(this%order, 100)

    ! if order > 200 and we're not dealing with the delta function piece
    ! then we make use of the fact that if we choose a very small or very 
    ! large value of alphas, we will be dominated by either a4pi^2 or a4pi^3
    ! term
    if (this%order > 200 .and. (.not. delta_one)) then
      if (mod_order == 2) then
        ! very small a4pi, so the relative contribution of a4pi^3 is
        ! negligible
        this_a4pi = 1e-80_dp
      else if (mod_order == 3) then
        ! very large a4pi, so the relative contribution of a4pi^2 is
        ! negligible
        this_a4pi = 1e80_dp
      else
        call wae_error('JB_fn', 'order > 200 and not mod_order 2 or 3; order = ', intval=this%order)
      end if
      res = fn(z, nf_light_dp, this_a4pi, zero)
      res = res / this_a4pi**mod_order
      return
    end if


    ! otherwise we will use up to two (or even three) finite a4pi values 
    ! solve the linear relation to extract the coefficients of a4pi^2 or a4pi^3

    ! fortran ternary, equivalent to C's delta_one ? 0 : 1
    imin = merge(0, 1, delta_one)
    ! special treatment that order > 100 means mixed orders
    imax = merge(2, 1, this%order > 100)

    ! for some reason the separate_orders function did not
    ! correctly receive the function pointer, so do things
    ! manually here    
    !print *, this%order, imin, imax, delta_one, mod_order
    do i = imin, imax
      ! last argument is LL = log(mQ^2/mu^2), taken to be zero
      res_a4pi(i) = fn(z, nf_light_dp, a4pi(i), zero)
    end do
    ! remove any constant piece
    if (delta_one) res_a4pi(1:imax) = res_a4pi(1:imax) - res_a4pi(0)

    if (imax == 2) then
      orders(3) = (res_a4pi(2)/a4pi(2)**2 - res_a4pi(1)/a4pi(1)**2)/(a4pi(2) - a4pi(1))
      orders(2) = (res_a4pi(1) - orders(3)*a4pi(1)**3)/a4pi(1)**2
      res = orders(mod_order)
    else
      res = res_a4pi(1) / a4pi(1)**mod_order
    end if
    !print*,res
  end function JB_fn
  
  ! Given a comparator object, compare the results of the Blumlein function
  ! and the Hoppet function at the point x
  subroutine compare(this, x)
    class(comparator), intent(in) :: this
    real(dp), intent(in) :: x
    real(dp) :: res_blum(2:3), res_hoppet
    !-------------------------------
    logical  :: delta_one
    real(dp) :: lcl_x, diff
    character(len=1) :: warning

    ! set the global cc_piece
    cc_piece = this%piece    
    ! nothing to compare
    if (cc_piece < 0) return

    if (associated(this%hoppet_fn)) then
      res_hoppet = this%hoppet_fn(log(one/x)) 
    else
      res_hoppet = zero
    end if
    ! for powers of x, take the convention from Blumlein
    if (this%piece /= cc_DELTA) res_hoppet = res_hoppet / x 
    if (this%piece == cc_VIRT)  res_hoppet = -res_hoppet

    ! establish if we will need to subtract 1 from the result 
    ! (because the terms with a delta function are of the form 1 + alpha^2...)
    delta_one = (this%piece == cc_DELTA)

    ! include a factor 0.5^n, so that the result multiplies 
    ! (as/2pi)^n
    res_blum = separate_orders(x, nf_int, this%blum_fn, delta_one) * (/ 0.25_dp, 0.125_dp /)

    lcl_x = x
    if (this%piece == cc_DELTA) lcl_x = one
    write(6,'(a8, f8.3)', advance='no') trim(this%name), lcl_x
    diff = res_blum(2) - res_hoppet
    if (abs(diff) > 1e-10) then
      warning = '*'
    else
      warning = ' '
    end if
    write(6,*) res_blum(2), res_hoppet, res_blum(2) - res_hoppet, warning
  end subroutine compare


  !! given a JB function fn that contains the sum of 2nd and 3rd order
  !! results, return the separate values of the 2nd and 3rd orders 
  !! (coefficients of (as/4pi)^n) by evaluating fn at two alphas values.
  !!
  !! if delta_one is .true., then we subtract one from the result
  function separate_orders(x, nf_light, fn, delta_one) result(res)
    real(dp), intent(in) :: x
    integer,  intent(in) :: nf_light
    procedure(JBfunc) :: fn
    !! if present, then we must subtract one from the result
    !! which has a delta(x-1)*(one)
    logical, intent(in), optional :: delta_one
    !interface 
    !  function fn(xx, nnf, aa4pi, LL) result(res)
    !    double precision, intent(in) :: xx
    !    integer,  intent(in) :: nnf
    !    double precision, intent(in) :: aa4pi, LL
    !    double precision :: res
    !  end function fn
    !end interface
    real(dp) :: res(2:3)
    !-------------------------------
    real(dp) :: a4pi(2) =  (/ 0.10_dp, 0.3_dp /)
    real(dp) :: res_a4pi(2), nf_light_dp
    integer  :: i

    !ffn => fn
    nf_light_dp = nf_light
    !write(0,*) fn(x, nf_light_dp, one, zero)
    do i = 1, size(a4pi)
      ! last argument is LL = log(mQ^2/mu^2), taken to be zero
      res_a4pi(i) = fn(x, nf_light_dp, a4pi(i), zero)
      if (default_or_opt(.false., delta_one)) res_a4pi(i) = res_a4pi(i) - one
    end do
    res(3) = (res_a4pi(2)/a4pi(2)**2 - res_a4pi(1)/a4pi(1)**2)/(a4pi(2) - a4pi(1))
    res(2) = (res_a4pi(1) - res(3)*a4pi(1)**3)/a4pi(1)**2
    !write(6,*) fn(x, nf, 0.0001_dp, zero), fn(x, nf, 0.00001_dp, zero)
  end function separate_orders

  !! Given a4pi(1:2) and res_a4pi(1:2) that contains a result at
  !! those two a4pi values, return the coefficients of the 2nd and 3rd
  !! powers of a4pi.
  function separate23(a4pi, res_a4pi, delta_one) result(res23)
    real(dp), intent(in) :: a4pi(2), res_a4pi(2)
    logical, intent(in), optional :: delta_one
    real(dp) :: res23(2:3), res_a4pi_lcl(2)
    !-------------------------------
    res_a4pi_lcl = res_a4pi
    if (default_or_opt(.false., delta_one)) res_a4pi_lcl = res_a4pi_lcl - one

    res23(3) = (res_a4pi_lcl(2)/a4pi(2)**2 - res_a4pi_lcl(1)/a4pi(1)**2)/(a4pi(2) - a4pi(1))
    res23(2) = (res_a4pi_lcl(1) - res23(3)*a4pi(1)**3)/a4pi(1)**2
  end function separate23

  ! a super inefficient (O N^2) way of getting the moment that also
  ! suffers from a truncation error...
  real(dp) function gc_moment(gc, momN) result(res)
    type(grid_conv), intent(in) :: gc
    real(dp), intent(in) :: momN
    !-------------------------------
    real(dp) :: powx(0:gc%grid%ny), conv_res(0:gc%grid%ny)
    powx = xValues(gc%grid)**(-(momN-1))
    conv_res = gc * powx
    res = conv_res(gc%grid%ny) / powx(gc%grid%ny)
  end function gc_moment

end module mass_thresholds_n3lo
