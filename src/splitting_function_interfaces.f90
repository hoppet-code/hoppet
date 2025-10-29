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
  
  public :: mvv_splitting_function, P_x_nfint

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
end module hoppet_splitting_function_interfaces