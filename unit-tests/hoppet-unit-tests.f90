program hoppet_unit_tests
  use types
  use unit_tests
  use interpolation
  use interpolation_coeffs
  implicit none

  print '(a)', "Running HOPPET unit tests"

  call test_interpolation_coeffs()

  if (unit_test_failures > 0) then
    ! print a message in red
    print '(a)', red//trim(to_string(unit_test_failures))//' unit tests failed.'//reset
    stop 1
  else
    ! print a message in green
    print '(a)', green//'All '//trim(to_string(unit_test_successes))//' unit tests passed.'//reset
  end if

contains

  !! checks equivalence between interpolation weights from
  !!  uniform_interpolation_weights and fill_interp_weightsN
  subroutine test_interpolation_coeffs()
    real(dp) :: x
    real(dp) :: weights(0:4)
    integer   :: i, ix
    integer, parameter :: nlo = 1, nhi = 4
    real(dp) :: weights1(0:nhi), weights2(0:nhi)
    real(dp), parameter :: xvals(3) = [0.2_dp,1.5_dp,3.4_dp]

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
          case Default
            call fail("Unsupported order "//trim(to_string(i))//" in test_interpolation_coeffs")
            continue
        end select

        call check_approx_eq("fill_interp_weightsN, order = "//trim(to_string(i))//", x="//trim(to_string(x)), &
                              weights2(0:i), weights1(0:i), tol = 1e-10_dp)

      end do
    end do

  end subroutine test_interpolation_coeffs


end program hoppet_unit_tests
