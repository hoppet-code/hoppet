module hoppet_libome_interfaces
  use types
  use, intrinsic :: iso_c_binding
  implicit none

  ! C getters for min/max powers of the generated OMEs
  interface
    integer(c_int) function ome_AqqQNSEven_reg_min_power() bind(C, name="ome_AqqQNSEven_reg_min_power")
      import c_int
    end function

    integer(c_int) function ome_AqqQNSEven_reg_max_power() bind(C, name="ome_AqqQNSEven_reg_max_power")
      import c_int
    end function

    integer(c_int) function ome_AqqQNSEven_plus_min_power() bind(C, name="ome_AqqQNSEven_plus_min_power")
      import c_int
    end function

    integer(c_int) function ome_AqqQNSEven_plus_max_power() bind(C, name="ome_AqqQNSEven_plus_max_power")
      import c_int
    end function

    integer(c_int) function ome_AqqQNSEven_delta_min_power() bind(C, name="ome_AqqQNSEven_delta_min_power")
      import c_int
    end function

    integer(c_int) function ome_AqqQNSEven_delta_max_power() bind(C, name="ome_AqqQNSEven_delta_max_power")
      import c_int
    end function
  end interface

  abstract interface
    !! a generic interface for functions in the libome library
    !!
    !! \param as4pi = (alpha_s/4pi), simply called "as" in libome
    !! \param LM    = log(m^2/mu^2)
    !! \param NF    = number of light flavours
    !! \param x     = Bjorken x
    function ome_as_LM_NF_x(as4pi, LM, NF, x) bind(C) result(res)
      import c_double
      real(c_double), value, intent(in) :: as4pi, LM, NF, x
      real(c_double) :: res
    end function

    function ome_as_LM_NF_delta(as, LM, NF) bind(C) result(res)
      import c_double
      real(c_double), value, intent(in) :: as, LM, NF
      real(c_double) :: res
    end function

    function ome_orderas_LM_NF_x(order_as, LM, NF, x) bind(C) result(res)
      import c_double, c_int
      integer(c_int), value, intent(in) :: order_as
      real(c_double), value, intent(in) :: LM, NF, x
      real(c_double) :: res
    end function

    function ome_orderas_LM_NF_delta(order_as, LM, NF) bind(C) result(res)
      import c_double, c_int
      integer(c_int), value, intent(in) :: order_as
      real(c_double), value, intent(in) :: LM, NF
      real(c_double) :: res
    end function
  end interface 

  procedure(ome_as_LM_NF_x    ), bind(C, name="ome_AqqQNSEven_reg")   :: ome_AqqQNSEven_reg
  procedure(ome_as_LM_NF_x    ), bind(C, name="ome_AqqQNSEven_plus")  :: ome_AqqQNSEven_plus
  procedure(ome_as_LM_NF_delta), bind(C, name="ome_AqqQNSEven_delta") :: ome_AqqQNSEven_delta

  procedure(ome_orderas_LM_NF_x),     bind(C, name="ome_AqqQNSEven_reg_coeff_as"  ) :: ome_AqqQNSEven_reg_coeff_as
  procedure(ome_orderas_LM_NF_x),     bind(C, name="ome_AqqQNSEven_plus_coeff_as" ) :: ome_AqqQNSEven_plus_coeff_as
  procedure(ome_orderas_LM_NF_delta), bind(C, name="ome_AqqQNSEven_delta_coeff_as") :: ome_AqqQNSEven_delta_coeff_as

  ! Declare a pointer to the AqqQNSEven object
  ! There are declared in hoppet_libome_interface.cc
  !type(c_ptr), bind(C, name="ome_AqqQNSEven") :: AqqQNSEven_ptr
  !type(c_ptr), bind(C, name="ome_AqqQNSOdd" ) :: AqqQNSOdd_ptr

  type(c_ptr), bind(C, name="ome_AggQ"      ) :: AggQ_ptr ;  
  type(c_ptr), bind(C, name="ome_AgqQ"      ) :: AgqQ_ptr ; 
  type(c_ptr), bind(C, name="ome_AQg"       ) :: AQg_ptr ; 
  type(c_ptr), bind(C, name="ome_AqgQ"      ) :: AqgQ_ptr ; 
  type(c_ptr), bind(C, name="ome_AQqPS"     ) :: AQqPS_ptr ;  
  type(c_ptr), bind(C, name="ome_AqqQNSEven") :: AqqQNSEven_ptr ; 
  type(c_ptr), bind(C, name="ome_AqqQNSOdd ") :: AqqQNSOdd_ptr ;  
  type(c_ptr), bind(C, name="ome_AqqQPS"    ) :: AqqQPS_ptr ; 


  interface 
    function ome_piece_hoppet(ptr,y,piece,order,LM,NF) bind(C, name="ome_piece_hoppet") result(res)
      use, intrinsic :: iso_c_binding
      real(c_double)                    :: res     !< return value: coefficient of (as/2pi)**order for xP(x), for given hoppet piece
      type(c_ptr)   ,        intent(in) :: ptr     !< pointer to one of the A[...] objects, i.e. an rpd_distribution object
      real(c_double),        intent(in) :: y       !< y = ln(1/x)
      integer(c_int),        intent(in) :: piece   !< which piece (cc_REAL, cc_REALVIRT, cc_PLUS, cc_DELTA)
      integer(c_int),        intent(in) :: order   !< order in as (0,1,2,3); returns a coefficient of (as/2pi)**order
      real(c_double),        intent(in) :: LM      !< LM = log(m^2/mu^2)
      real(c_double),        intent(in) :: NF      !< number of light flavours
    end function
  end interface
end module hoppet_libome_interfaces