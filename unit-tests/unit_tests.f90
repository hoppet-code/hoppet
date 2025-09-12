module unit_tests
  use types
  implicit none

  interface to_string
     module procedure to_string_dp, to_string_int
  end interface
  interface check_approx_eq
     module procedure check_approx_eq_1d
  end interface

  character(len=*), parameter :: red = achar(27)//'[31m'
  character(len=*), parameter :: green = achar(27)//'[32m'
  character(len=*), parameter :: reset = achar(27)//'[0m'

  integer :: unit_test_failures = 0
  integer :: unit_test_successes = 0

contains

  subroutine check_approx_eq_1d(testname, answer, expected, tol)

    implicit none
    character(len=*), intent(in) :: testname
    real(dp), intent(in) :: answer(:)
    real(dp), intent(in) :: expected(:)
    real(dp), intent(in) :: tol
    logical :: is_equal

    is_equal = all(abs(answer - expected) < tol)
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

  function to_string_dp(x) result(str)
    implicit none
    real(dp), intent(in) :: x
    character(len=100) :: str

    write(str, *) x
    str = adjustl(str)
  end function to_string_dp

  function to_string_int(x) result(str)
    implicit none
    integer, intent(in) :: x
    character(len=100) :: str

    write(str, *) x
    str = adjustl(str)
  end function to_string_int


end module unit_tests

