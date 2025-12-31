!! module with 
module hoppet_splitting_function_interfaces
  use types, only: dp
  use convolution, only: conv_ignd
  implicit none
  private
  

  !! the interface for MVV functions like X2NSPA etc in xpns2e etc.
  !! 
  !! Note that the arguments are const in practice, but not
  !! declared as such, so to ensure we match the exact form
  !! of the functions in the splitting function modules
  abstract interface
    !! returns P(x) for given nf
    real(dp) function P_x_nfint(x, nf) 
      import dp
      implicit none
      real(dp) :: x   !< x value
      integer  :: nf  !< number of light flavours
    end function P_x_nfint
  end interface

  !! the interface for MVV functions like X2NSPA etc in xpns2e etc.
  !! 
  !! Note that the arguments are const in practice, but not ! declared
  !as such, so to ensure we match the exact form ! of the functions in
  !the splitting function modules. This version is needed when the
  !splitting functions in question have an extra "mode" argument, as
  !is the case for the approximate splitting functions.
  abstract interface
    !! returns P(x) for given nf and imod
    real(dp) function P_x_nfint_imod(x, nf, imod) 
      import dp
      implicit none
      real(dp) :: x   !< x value
      integer  :: nf, imod  !< number of light flavours and imod
    end function P_x_nfint_imod
  end interface

  !! the interfaces for MVV n3lo coefficient functions functions.
  !! 
  abstract interface
    !! returns C(x) for given nf and imod
    real(dp) function C_x_y_nfint_imod(x, y, nf, imod) 
      import dp
      implicit none
      real(dp) :: x, y   !< x and y values
      integer  :: nf, imod  !< number of light flavours
    end function C_x_y_nfint_imod
  end interface
  abstract interface
    !! returns C(x) for given nf and imod
    real(dp) function C_x_y_nfint(x, y, nf) 
      import dp
      implicit none
      real(dp) :: x, y   !< x and y values
      integer  :: nf  !< number of light flavours
    end function C_x_y_nfint
  end interface

  !! a type that holds pointers to the regular, plus and delta parts
  !! of an MVV splitting function and implements the conv_ignd%f(y,piece) 
  !! interface
  !!
  !! Users would typically create an mvv_splitting_function via the
  !! as 
  !!
  !!     mvv_splitting_function(P2GGA, P2GGB , P2GGC , 0.5_dp**3)
  !!
  !! where the last argument is the multiplier, e.g. here to convert from
  !! (alpha_s/4pi)**order to (alpha_s/2pi)**order
  !!
  !! Where a part is not present, the corresponding argument should be null()
  !! e.g. 
  !!
  !!     mvv_splitting_function(P2QGA, null(), null(), 0.5_dp**3)
  !!
  type, extends(conv_ignd) :: mvv_splitting_function
    procedure(P_x_nfint), pointer, nopass :: reg   => null()  !< A, regular part 
    procedure(P_x_nfint), pointer, nopass :: plus  => null()  !< B, plus part
    procedure(P_x_nfint), pointer, nopass :: delta => null()  !< C, delta-function part (call with x=0)    
    real(dp) :: multiplier !< multiplier for the function
  contains
    procedure :: f => mvv_splitting_function__f  !< f(y=ln1/x, piece)
  end type mvv_splitting_function

  !! a type that holds pointers to the regular, plus and delta parts
  !! of an MVV splitting function with an imod parameter and implements 
  !! the conv_ignd%f(y,piece) interface
  !!
  !! This version is for splitting functions with an extra "mode" argument
  !!
  !! Users would typically create an mvv_splitting_function_imod via:
  !!
  !!     mvv_splitting_function_imod(P2GGA_imod, P2GGB_imod , P2GGC_imod , imod_val, 0.5_dp**3)
  !!
  !! where the imod_val is the mode argument passed to the splitting functions
  !!
  type, extends(conv_ignd) :: mvv_splitting_function_imod
    procedure(P_x_nfint_imod), pointer, nopass :: reg   => null()  !< A, regular part 
    procedure(P_x_nfint_imod), pointer, nopass :: plus  => null()  !< B, plus part
    procedure(P_x_nfint_imod), pointer, nopass :: delta => null()  !< C, delta-function part (call with x=0)    
    integer  :: imod !< mode argument for the splitting functions
    real(dp) :: multiplier !< multiplier for the function
  contains
    procedure :: f => mvv_splitting_function_imod__f  !< f(y=ln1/x, piece)
  end type mvv_splitting_function_imod
  
  ! Alias types for coefficient functions — identical to the splitting
  ! function types but with different names so callers can use a
  ! semantically clearer type without duplicating code.
  
  type, extends(mvv_splitting_function) :: mvv_coefficient_function
  end type mvv_coefficient_function
  type, extends(mvv_splitting_function_imod) :: mvv_coefficient_function_imod
  end type mvv_coefficient_function_imod

  type, extends(conv_ignd) :: mvv_coefficient_function_n3lo
    procedure(C_x_y_nfint_imod), pointer, nopass :: reg   => null()  !< A, regular part 
    procedure(C_x_y_nfint), pointer, nopass :: plus  => null()  !< B, plus part
    procedure(P_x_nfint), pointer, nopass :: delta => null()  !< C, delta-function part (call with x=0)    
    integer  :: imod !< mode argument for the splitting functions
    real(dp) :: multiplier !< multiplier for the function
  contains
    procedure :: f => mvv_coefficient_function_n3lo__f  !< f(y=ln1/x, piece)
  end type mvv_coefficient_function_n3lo

  type, extends(conv_ignd) :: mvv_coefficient_function_n3lo_exact
    procedure(C_x_y_nfint_imod), pointer, nopass :: reg   => null()  !< A, regular part 
    procedure(C_x_y_nfint), pointer, nopass :: plus  => null()  !< B, plus part
    procedure(C_x_y_nfint), pointer, nopass :: delta => null()  !< C, delta-function part (call with x=0)    
    integer  :: imod !< mode argument for the splitting functions
    real(dp) :: multiplier !< multiplier for the function
  contains
    procedure :: f => mvv_coefficient_function_n3lo_exact__f  !< f(y=ln1/x, piece)
  end type mvv_coefficient_function_n3lo_exact

  public :: mvv_splitting_function, mvv_splitting_function_imod, &
            mvv_coefficient_function, mvv_coefficient_function_n3lo, &
            mvv_coefficient_function_n3lo_exact, mvv_coefficient_function_imod

contains


  !! implementation of the mvv_splitting_function%f
  function mvv_splitting_function__f(this, y, piece) result(res)
    use qcd, only: nf_int
    use convolution_pieces
    class(mvv_splitting_function), intent(in) :: this
    real(dp),                      intent(in) :: y
    integer,                       intent(in) :: piece
    real(dp) :: res
    real(dp) :: x
    integer  :: nf_lcl

    x = exp(-y)
    nf_lcl = nf_int ! get the global nf_int
    res = 0.0_dp

    select case(piece)
    case(cc_REAL)
      if (associated(this%reg  )) res =       this%reg (x, nf_lcl)
      if (associated(this%plus )) res = res + this%plus(x, nf_lcl)
    case(cc_REALVIRT)
      if (associated(this%reg  )) res = res + this%reg (x, nf_lcl)
    case(cc_VIRT)
      if (associated(this%plus )) res =     - this%plus(x, nf_lcl)
    case(cc_DELTA)
      if (associated(this%delta)) res =       this%delta(0.0_dp, nf_lcl)
    end select
    if (piece /= cc_DELTA) res = res * x

    res = res * this%multiplier

    !write(6,*) y, piece, nf_lcl, res
  end function mvv_splitting_function__f

  !! implementation of the mvv_splitting_function_imod%f
  function mvv_splitting_function_imod__f(this, y, piece) result(res)
    use qcd, only: nf_int
    use convolution_pieces
    class(mvv_splitting_function_imod), intent(in) :: this
    real(dp),                          intent(in) :: y
    integer,                           intent(in) :: piece
    real(dp) :: res
    real(dp) :: x
    integer  :: nf_lcl

    x = exp(-y)
    nf_lcl = nf_int ! get the global nf_int
    res = 0.0_dp

    select case(piece)
    case(cc_REAL)
      if (associated(this%reg  )) res =       this%reg (x, nf_lcl, this%imod)
      if (associated(this%plus )) res = res + this%plus(x, nf_lcl, this%imod)
    case(cc_REALVIRT)
      if (associated(this%reg  )) res = res + this%reg (x, nf_lcl, this%imod)
    case(cc_VIRT)
      if (associated(this%plus )) res =     - this%plus(x, nf_lcl, this%imod)
    case(cc_DELTA)
      if (associated(this%delta)) res =       this%delta(0.0_dp, nf_lcl, this%imod)
    end select
    if (piece /= cc_DELTA) res = res * x

    res = res * this%multiplier

    !write(6,*) y, piece, nf_lcl, res
  end function mvv_splitting_function_imod__f

  !! implementation of the mvv_coefficient_function_n3lo%f
  function mvv_coefficient_function_n3lo__f(this, y, piece) result(res)
    use qcd, only: nf_int
    use convolution_pieces
    class(mvv_coefficient_function_n3lo), intent(in) :: this
    real(dp),                             intent(in) :: y
    integer,                              intent(in) :: piece
    real(dp) :: res
    real(dp) :: x
    integer  :: nf_lcl

    x = exp(-y)
    nf_lcl = nf_int ! get the global nf_int
    res = 0.0_dp

    select case(piece)
    case(cc_REAL)
      if (associated(this%reg  )) res =       this%reg (x, -y, nf_lcl, this%imod)
      if (associated(this%plus )) res = res + this%plus(x, -y, nf_lcl)
    case(cc_REALVIRT)
      if (associated(this%reg  )) res = res + this%reg (x, -y, nf_lcl, this%imod)
    case(cc_VIRT)
      if (associated(this%plus )) res =     - this%plus(x, -y, nf_lcl)
    case(cc_DELTA)
      if (associated(this%delta)) res =       this%delta(0.0_dp, nf_lcl)
    end select
    if (piece /= cc_DELTA) res = res * x

    res = res * this%multiplier

    !write(6,*) y, piece, nf_lcl, res
  end function mvv_coefficient_function_n3lo__f

  !! implementation of the mvv_coefficient_function_n3lo_exact%f
  function mvv_coefficient_function_n3lo_exact__f(this, y, piece) result(res)
    use qcd, only: nf_int
    use convolution_pieces
    class(mvv_coefficient_function_n3lo_exact), intent(in) :: this
    real(dp),                             intent(in) :: y
    integer,                              intent(in) :: piece
    real(dp) :: res
    real(dp) :: x
    integer  :: nf_lcl

    x = exp(-y)
    nf_lcl = nf_int ! get the global nf_int
    res = 0.0_dp

    select case(piece)
    case(cc_REAL)
      if (associated(this%reg  )) res =       this%reg (x, -y, nf_lcl, this%imod)
      if (associated(this%plus )) res = res + this%plus(x, -y, nf_lcl)
    case(cc_REALVIRT)
      if (associated(this%reg  )) res = res + this%reg (x, -y, nf_lcl, this%imod)
    case(cc_VIRT)
      if (associated(this%plus )) res =     - this%plus(x, -y, nf_lcl)
    case(cc_DELTA)
      if (associated(this%delta)) res =       this%delta(0.0_dp, -y, nf_lcl)
    end select
    if (piece /= cc_DELTA) res = res * x

    res = res * this%multiplier

    !write(6,*) y, piece, nf_lcl, res
  end function mvv_coefficient_function_n3lo_exact__f

end module hoppet_splitting_function_interfaces