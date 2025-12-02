!=====================================================================
!! utilities for converting between C-style strings and Fortran strings  
module hoppet_c_f_string_utils
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    function hoppet_cstr_len(cstr) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: cstr
      integer(c_int) :: hoppet_cstr_len
    end function hoppet_cstr_len
  end interface

  interface
    function hoppet_allocate_cstr(size) bind(C) result(cstr)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: size
      type(c_ptr) :: cstr
    end function hoppet_allocate_cstr
  end interface


contains

  !! Take a c_ptr to a C-style array of characters, null terminated, and
  !! convert to a Fortran string, which is returned
  !!
  !! If the optional errcode is provided, it is set to 0 on success,
  !! or a non-zero value if an error occurred (e.g. string too long)
  function fortran_string_from_cstr(cstr_ptr) result(fstr)
    implicit none
    type(c_ptr), intent(in), value :: cstr_ptr
    character(len=:), allocatable :: fstr
    character(kind=c_char,len=1), pointer :: cchars(:)
    integer :: len, i

    len = hoppet_cstr_len(cstr_ptr)

    ! convert c_ptr to Fortran pointer to array of characters
    call c_f_pointer(cstr_ptr, cchars, shape=[len])

    ! allocate Fortran string and copy characters
    allocate(character(len=len) :: fstr)
    do i = 1, len
      fstr(i:i) = cchars(i)
    end do

  end function fortran_string_from_cstr

  !! Take a Fortran string and convert to a C-style array of characters,
  !! null terminated, returning a cstr_ptr 
  !!
  !! The allocation is carried out inside a dedicated C++ function
  !! and whatever C++ function received the c_chars takes responsibility
  !! for deallocating the memory with `delete[] cstr_ptr
  function cstr_from_fortran_string(fstr) result(cstr_ptr)
    implicit none
    character(len=*), intent(in) :: fstr
    type(c_ptr) :: cstr_ptr
    character(kind=c_char,len=1), pointer :: cchars(:)
    integer :: len, i

    len = len_trim(fstr)

    cstr_ptr = hoppet_allocate_cstr(len + 1)  ! +1 for null terminator
    call c_f_pointer(cstr_ptr, cchars, shape=[len+1]) 

    ! allocate C-style array of characters (with null terminator)
    ! copy characters from Fortran string to C-style array
    do i = 1, len
      cchars(i) = fstr(i:i)
    end do
    cchars(len+1) = c_null_char  ! null terminator

  end function cstr_from_fortran_string

end module hoppet_c_f_string_utils
