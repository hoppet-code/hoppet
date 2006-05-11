!======================================================================
! $Id: coefficient_functions.f90,v 1.1 2001/06/27 14:15:16 gsalam Exp $
module coefficient_functions
  use types; use consts_dp; use convolution_communicator
  use qcd
  implicit none
  private


  public :: cf_CqF2MSbar, cf_CgF2MSbar
  public :: cf_CqFL, cf_CgFL

contains

  function cf_CqF2MSbar(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x); ln1mx = log(one - x)
       res = CF*(two*(ln1mx/(1-x)) - 1.5_dp/(1-x) - (1+x)*ln1mx&
            & - (1+x**2)/(1-x)*lnx + 3 + 2*x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ln1mx = log(one - x) ! spiacente di farlo due volte
       res = res - CF*(two*(ln1mx/(1-x)) - 1.5_dp/(1-x))
    case(cc_DELTA)
       res =  -CF*(pi**2/3.0_dp + 9.0_dp*half)
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function cf_CqF2MSbar


  function cf_CgF2MSbar(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x); ln1mx = log(one - x)
       res = TR*(((1-x)**2+x**2)*(ln1mx-lnx)-8*x**2+8*x-1)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function cf_CgF2MSbar

  function cf_CqFL(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = two*CF*x
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function cf_CqFL

  function cf_CgFL(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = four*TR * x*(1-x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function cf_CgFL
  
!.
!.
!.

end module coefficient_functions
