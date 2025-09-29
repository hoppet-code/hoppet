module hoppet_pdf_table_array
  use types
  use pdf_tabulate
  use iso_c_binding
  implicit none
  
  type pdf_table_array
    type(pdf_table), allocatable :: tables(:) ! array of tables
  end type pdf_table_array

contains 

  !! creates a pdf_table_array object with space for ntables tables, indexed 0..ntables-1
  !! and returns a c_ptr to it
  function hoppet_pta_create(ntables) bind(c,name="hoppet_pta_create") result(ptr)
    integer(c_int), intent(in) :: ntables
    type(c_ptr) :: ptr

    type(pdf_table_array), pointer :: obj
    allocate(obj)
    allocate(obj%tables(0:ntables-1))    
    ptr = c_loc(obj)
  end function hoppet_pta_create

  !! returns the size of the array of tables in the pdf_table_array object pointed to by ptr
  function hoppet_pta_size(ptr) bind(c,name="hoppet_pta_size") result(sz)
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int) :: sz

    type(pdf_table_array), pointer :: obj
    call c_f_pointer(ptr, obj)
    sz = size(obj%tables)
  end function hoppet_pta_size

end module hoppet_pdf_table_array
