module hoppet_libome_fortran
  use types
  use convolution, only: conv_ignd
  use, intrinsic :: iso_c_binding
  implicit none

  ! Pointers to the various A[...] objects in libome, initialised in 
  ! the hoppet_libome_interface.cc file
  type(c_ptr), bind(C, name="ome_AggQ"      ) :: AggQ_ptr ;  
  type(c_ptr), bind(C, name="ome_AgqQ"      ) :: AgqQ_ptr ; 
  type(c_ptr), bind(C, name="ome_AQg"       ) :: AQg_ptr ; 
  type(c_ptr), bind(C, name="ome_AqgQ"      ) :: AqgQ_ptr ; 
  type(c_ptr), bind(C, name="ome_AQqPS"     ) :: AQqPS_ptr ;  
  type(c_ptr), bind(C, name="ome_AqqQNSEven") :: AqqQNSEven_ptr ; 
  type(c_ptr), bind(C, name="ome_AqqQNSOdd ") :: AqqQNSOdd_ptr ;  
  type(c_ptr), bind(C, name="ome_AqqQPS"    ) :: AqqQPS_ptr ; 


  !! interface to the ome_piece_hoppet function in libome
  !! which provides the piecewise components of the splitting functions
  !! in the form required by hoppet
  interface 
    function ome_piece_hoppet(ptr,y,piece,order,LM,NF) bind(C, name="ome_piece_hoppet") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double)                    :: res     !< return value: coefficient of (as/2pi)**order for xP(x), for given hoppet piece
      type(c_ptr)   , value, intent(in) :: ptr     !< pointer to one of the A[...] objects, i.e. an rpd_distribution object
      real(c_double),        intent(in) :: y       !< y = ln(1/x)
      integer(c_int),        intent(in) :: piece   !< which piece (cc_REAL, cc_REALVIRT, cc_PLUS, cc_DELTA)
      integer(c_int),        intent(in) :: order   !< order in as (0,1,2,3); returns a coefficient of (as/2pi)**order
      real(c_double),        intent(in) :: LM      !< LM = log(m^2/mu^2)
      real(c_double),        intent(in) :: NF      !< number of light flavours
    end function
  end interface

  !! Derived type implementing the conv_ignd abstract type
  !! using the libome library.
  !!
  !! Initialise an object as conv_OME(A..._ptr, order, ...) with
  !! parameters as follows
  !!
  !! \param ptr      pointer to one of the A..._ptr objects above
  !! \param order    the order in alphas to use (0,1,2,3)
  !! \param nf_light if >0, use this number of light flavours, otherwise the global nf-1 (defaults to global)
  !! \param LM       the LM value to use (default 0.0)  
  !!
  !! Example usage:
  !!
  !!    type(grid_conv) :: A3ggQ
  !!    call InitGridConv(grid, A3ggQ, conv_OME(AggQ_ptr,order=3))
  !!
  type, extends(conv_ignd) :: conv_OME
    type(c_ptr)    :: ptr                !< pointer to one of the A[...]_ptr objects
    integer(c_int) :: order              !< the order in alphas    
    real(c_double) :: nf_light = -1.0_c_double !< if >0, use this nf instead of the one from qcd module
    real(c_double) :: LM = 0.0_c_double  !< the LM value to use (mainly for testing)
  contains
    procedure :: f => conv_OME__f !< f(y,piece) function
  end type

contains
  
  !! Implementation of the conv_ignd%f(y,piece) function for the
  !! conv_OME derived type
  function conv_OME__f(this, y, piece) result(res)
    use qcd, only: nf
    implicit none
    class(conv_OME), intent(in) :: this
    real(dp), intent(in) :: y
    integer , intent(in) :: piece
    real(c_double) :: res, nf_c, lm_c

    if (this%nf_light > 0.0_c_double) then
      nf_c = this%nf_light
    else
      nf_c = real(nf, c_double) - 1.0_c_double ! subtract 1 to get number of light flavours, as per hoppet convention
    end if

    ! Call the ome_piece_hoppet function with the stored pointer and order
    res = ome_piece_hoppet(this%ptr, y, piece, this%order, this%LM, nf_c)
  end function conv_OME__f
end module hoppet_libome_fortran