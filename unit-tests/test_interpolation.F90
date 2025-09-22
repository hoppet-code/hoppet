#include "timing.inc"

module test_interpolation
  use types
  use unit_tests
  implicit none

contains

  !! checks equivalence between interpolation weights from
  !! uniform_interpolation_weights and fill_interp_weightsN
  subroutine test_interpolation_coeffs()
    use interpolation
    use interpolation_coeffs
    real(dp) :: x
    real(dp) :: weights(0:4)
    integer   :: i, ix
    integer, parameter :: nlo = 1, nhi = 6
    real(dp) :: weights1(0:nhi), weights2(0:nhi)
    real(dp), parameter :: xvals(*) = [0.2_dp,1.5_dp,3.4_dp, 1.0_dp]

    do ix = 1, size(xvals)
      x = xvals(ix)
      do i = nlo, nhi
        weights1 = 0.0_dp
        call uniform_interpolation_weights(x, weights1(0:i))
        select case (i)
          case (1)
            call fill_interp_weights1(x, weights2)
          case (2)
            call fill_interp_weights2(x, weights2)
          case (3)
            call fill_interp_weights3(x, weights2)
          case (4)
            call fill_interp_weights4(x, weights2)
          case (5)
            call fill_interp_weights5(x, weights2)
          case (6)
            call fill_interp_weights6(x, weights2)
          case Default
            call fail("Unsupported order "//trim(to_string(i))//" in test_interpolation_coeffs")
            continue
        end select

        call check_approx_eq("fill_interp_weightsN, order = "//trim(to_string(i))//", x="//trim(to_string(x)), &
                              weights2(0:i), weights1(0:i), tol_abs = 1e-10_dp)
      end do
    end do

  ! some basic timing code
  ! NB: it only runs if timing_name=="fill_interp_weights6:x4"
  TIME_START("fill_interp_weights6:x4",10**5)
  real(dp) :: weights_sum(0:6) = 0.0_dp
  TIME_LOOP
    ! cover a range of x values to include special cases & not
    do ix = 1, size(xvals)
      call fill_interp_weights6(xvals(ix), weights2)
      weights_sum = weights_sum + weights2
    end do
  TIME_END(weights_sum(0))

  end subroutine test_interpolation_coeffs


end module test_interpolation  