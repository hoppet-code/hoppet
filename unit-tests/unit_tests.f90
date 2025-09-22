module unit_tests
  use types
  use hoppet_to_string
  use io_utils
  implicit none

  ! anything matching this will be timed
  character(len=max_arg_len):: timing_name = "fill_interp_weights6:x4"
  logical :: list_timing = .false.

  interface check_approx_eq
     module procedure check_approx_eq_0d, check_approx_eq_1d
  end interface

  character(len=*), parameter :: red = achar(27)//'[31m'
  character(len=*), parameter :: green = achar(27)//'[32m'
  character(len=*), parameter :: blue = achar(27)//'[34m'
  character(len=*), parameter :: reset = achar(27)//'[0m'

  integer :: unit_test_failures = 0
  integer :: unit_test_successes = 0

contains

  subroutine check_approx_eq_0d(testname, answer, expected, tol_abs, tol_rel)
    implicit none
    character(len=*), intent(in) :: testname
    real(dp), intent(in) :: answer
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tol_abs
    real(dp), intent(in), optional :: tol_rel
    logical :: is_equal

    is_equal = (abs(answer - expected) < tol_abs)
    if (present(tol_rel)) then
      is_equal = is_equal .and. (abs(answer - expected) < tol_rel*max(abs(expected), abs(answer)))
    end if

    if (.not. is_equal) then
      print *, red//"Failed: ", testname,reset
      print *, "  Expected: ", expected
      print *, "  Actual:   ", answer
      unit_test_failures = unit_test_failures + 1
    else
      unit_test_successes = unit_test_successes + 1
    end if
  end subroutine check_approx_eq_0d

  subroutine check_approx_eq_1d(testname, answer, expected, tol_abs, tol_rel)
    implicit none
    character(len=*), intent(in) :: testname
    real(dp), intent(in) :: answer(:)
    real(dp), intent(in) :: expected(:)
    real(dp), intent(in) :: tol_abs
    real(dp), intent(in), optional :: tol_rel
    logical :: is_equal

    is_equal = all(abs(answer - expected) < tol_abs)
    if (present(tol_rel)) then
      is_equal = is_equal .and. all(abs(answer - expected) < tol_rel*max(abs(expected), abs(answer)))
    end if

    if (.not. is_equal) then
      print *, red//"Failed: ", testname,reset
      print *, "  Expected: ", expected
      print *, "  Actual:   ", answer
      unit_test_failures = unit_test_failures + 1
    else
      unit_test_successes = unit_test_successes + 1
    end if

  end subroutine check_approx_eq_1d

  subroutine fail(testname)
    implicit none
    character(len=*), intent(in) :: testname
    print *, red//"Failed: ", testname, reset
    unit_test_failures = unit_test_failures + 1
  end subroutine fail



end module unit_tests

