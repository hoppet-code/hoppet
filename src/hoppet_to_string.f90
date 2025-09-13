!! Module that provides a simple interface for converting various numerical
!! types to strings
module hoppet_to_string
  use types
  implicit none

  interface to_string
     module procedure to_string_dp, to_string_int
  end interface

contains

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

end module hoppet_to_string