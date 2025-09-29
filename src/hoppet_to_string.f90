!! Module that provides a simple interface for converting various numerical
!! types to strings
module hoppet_to_string
  use types
  implicit none

  interface to_string
     module procedure to_string_dp, to_string_int
  end interface

contains

  function to_string_dp(x,fmt) result(str)
    implicit none
    real(dp), intent(in) :: x
    character(len=*), intent(in), optional :: fmt
    character(len=:), allocatable :: str 
    character(len=100) :: tmp

    if (present(fmt)) then
       write(tmp, fmt) x
    else
       write(tmp, *) x
    end if
    tmp = adjustl(tmp)
    str = tmp(1:len_trim(tmp))
  end function to_string_dp

  function to_string_int(x) result(str)
    implicit none
    integer, intent(in) :: x
    character(len=:), allocatable :: str 
    character(len=100) :: tmp

    write(tmp, *) x
    tmp = adjustl(tmp)
    str = tmp(1:len_trim(tmp))
  end function to_string_int

end module hoppet_to_string