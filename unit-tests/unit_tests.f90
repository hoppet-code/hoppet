module unit_tests
  use types
  use hoppet_to_string
  use hoppet_term
  use io_utils
  implicit none

  ! anything matching this will be timed
  character(len=max_arg_len):: timing_name = "fill_interp_weights6:x4"
  logical :: list_timing = .false.

  interface check_approx_eq
     module procedure check_approx_eq_0d, check_approx_eq_1d
  end interface

  integer :: unit_test_failures = 0
  integer :: unit_test_successes = 0

  character(len=:), allocatable :: current_test_group
  character(len=:), allocatable :: test_only
  logical :: verbose = .false.

contains

  logical function do_test(test_group)
    character(len=*), intent(in) :: test_group
    current_test_group = test_group
    if (allocated(test_only)) then
      if (len_trim(test_only) == 0) then
        do_test = .true.
      else
        do_test = (index(test_group, trim(test_only)) > 0)
      end if
    else
      do_test = .true.
    end if
    if (do_test) then
      print '(a)', blue//bold//"Running test group: "//trim(test_group)//reset
    else
      print '(a)', yellow//bold//"Skipping test group: "//trim(test_group)//reset//" (not in -only)"
    end if

  end function do_test

  subroutine check_approx_eq_0d(testname, answer, expected, tol_abs, tol_rel, tol_choice_or)
    use assertions
    implicit none
    character(len=*), intent(in) :: testname
    real(dp), intent(in) :: answer
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tol_abs
    real(dp), intent(in), optional :: tol_rel
    logical, intent(in),  optional :: tol_choice_or ! if true, use OR for tol checks, else (default) AND
    logical :: is_equal, is_equal_abs, is_equal_rel

    is_equal_abs = (abs(answer - expected) < tol_abs)
    if (present(tol_rel)) then
      is_equal_rel = (abs(answer - expected) < tol_rel*max(abs(expected), abs(answer)))
      if (default_or_opt(.false., tol_choice_or)) then
        is_equal = is_equal_abs .or. is_equal_rel
      else
        is_equal = is_equal_abs .and. is_equal_rel
      end if
    else
      is_equal = is_equal_abs
    end if

    if (.not. is_equal) then
      print *, ""
      print '(a)', red//bold//"Failed: "//trim(testname)//reset
      print *, "  Expected: ", expected
      print *, "  Actual:   ", answer
      print *, "  abs tol:  ", tol_abs
      if (present(tol_rel)) print *, "  rel tol:  ", tol_rel
      print *, "  difference: ", answer - expected,  ", rel diff: ", (answer - expected)/max(abs(expected), abs(answer))
      print *, ""
      unit_test_failures = unit_test_failures + 1
    else
      unit_test_successes = unit_test_successes + 1
      if (verbose) print '(a)', green//"Passed: "//testname//reset//" ans = "//&
                          trim(to_string(answer))//" matches expected = "//trim(to_string(expected))//reset
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

