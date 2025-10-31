module hoppet_cxx_oo
  use types, only: dp
  use, intrinsic :: iso_c_binding
  use convolution
  implicit none
  
  type grid_quant
    real(dp), pointer :: data(:) => null()
    type(grid_def)    :: grid
  end type grid_quant

contains


  !! return a c pointer to a new grid_def
  function hoppet_cxx__grid_def__new(dy,ymax,order,eps) bind(C) result(ptr)
    use convolution
    implicit none
    real(c_double), intent(in), value :: dy, ymax
    integer(c_int), intent(in), value :: order
    real(c_double), intent(in), value :: eps
    type(c_ptr) :: ptr
    !--
    real(dp)  :: f_dy, f_ymax
    integer   :: f_order
    real(dp)  :: f_eps
    type(grid_def), pointer :: f_ptr

    allocate(f_ptr)
    f_dy   = dy
    f_ymax = ymax
    f_order= order
    f_eps  = eps
    call InitGridDef(f_ptr, f_dy, f_ymax, f_order, f_eps)
    ptr = c_loc(f_ptr)
  end function hoppet_cxx__grid_def__new

  !! return a c pointer to a copy of the (grid_c_ptr) grid_def
  function hoppet_cxx__grid_def__copy(grid_c_ptr) bind(C) result(ptr)
    use convolution
    implicit none
    type(c_ptr),  intent(in), value :: grid_c_ptr
    type(c_ptr) :: ptr
    !----------------------
    type(grid_def), pointer :: grid_ptr
    type(grid_def), pointer :: f_ptr

    ! prepare fortran objects
    call c_f_pointer(grid_c_ptr, grid_ptr)
    allocate(f_ptr)

    if (associated(grid_ptr%subgd)) then
      call InitGridDef(f_ptr, grid_ptr%subgd, grid_ptr%locked)
    else
      f_ptr = grid_ptr
    end if
    ptr = c_loc(f_ptr)
  end function hoppet_cxx__grid_def__copy

  function hoppet_cxx__grid_def__new_from_grids(grid_c_ptrs, ngrids, locked) bind(C) result(ptr)
    use convolution
    implicit none
    type(c_ptr),     intent(in)        :: grid_c_ptrs(*)
    integer(c_int),  intent(in), value :: ngrids
    logical(c_bool), intent(in), value :: locked
    type(c_ptr) :: ptr
    !----------------------
    type(grid_def) :: grids(ngrids)
    type(grid_def), pointer :: grid_ptr
    type(grid_def), pointer :: f_ptr
    logical :: f_locked
    integer :: i

    ! Create a Fortran grid_def object for each grid
    do i = 1, ngrids
      call c_f_pointer(grid_c_ptrs(i), grid_ptr)
      grids(i) = grid_ptr
    end do
    f_locked = merge(.true., .false., locked)

    allocate(f_ptr)
    call InitGridDef(f_ptr, grids, f_locked)
    ptr = c_loc(f_ptr)
  end function hoppet_cxx__grid_def__new_from_grids

  function hoppet_cxx__grid_def__new_default(dy, ymax, order) bind(C) result(ptr)
    use convolution
    implicit none
    real(c_double), intent(in), value :: dy, ymax
    integer(c_int), intent(in), value :: order
    type(c_ptr) :: ptr
    !--
    real(dp)  :: f_dy, f_ymax
    integer   :: f_order
    type(grid_def), pointer :: f_ptr

    allocate(f_ptr)
    f_dy   = dy
    f_ymax = ymax
    f_order= order
    call InitGridDefDefault(f_ptr, f_dy, f_ymax, f_order)
    ptr = c_loc(f_ptr)
  end function hoppet_cxx__grid_def__new_default

  subroutine hoppet_cxx__grid_def__delete(ptr) bind(C)
    use convolution
    implicit none
    type(c_ptr), intent(inout) :: ptr
    !--
    type(grid_def), pointer :: f_ptr

    call c_f_pointer(ptr, f_ptr)
    !print * , associated(f_ptr)
    call delete(f_ptr)
    deallocate(f_ptr)
    ptr = c_null_ptr
  end subroutine hoppet_cxx__grid_def__delete

  function hoppet_cxx__grid_def__ny(ptr) bind(C) result(ny)
    use convolution
    implicit none
    type(c_ptr), intent(in), value :: ptr
    integer(c_int) :: ny
    type(grid_def), pointer :: f_ptr

    call c_f_pointer(ptr, f_ptr)
    ny = f_ptr%ny
  end function hoppet_cxx__grid_def__ny

  subroutine hoppet_cxx__grid_def__y_values(ptr, yvals) bind(C)
    use convolution
    implicit none
    type(c_ptr), intent(in), value :: ptr
    real(c_double), intent(out) :: yvals(*)
    type(grid_def), pointer :: grid

    call c_f_pointer(ptr, grid)
    yvals(1:grid%ny+1) = yValues(grid)
  end subroutine

  subroutine hoppet_cxx__grid_def__x_values(ptr, xvals) bind(C)
    use convolution
    implicit none
    type(c_ptr), intent(in), value :: ptr
    real(c_double), intent(out) :: xvals(*)
    type(grid_def), pointer :: grid

    call c_f_pointer(ptr, grid)
    xvals(1:grid%ny+1) = xValues(grid)
  end subroutine


  function hoppet_cxx__grid_quant__new(grid_def_ptr) bind(C) result(ptr)
    implicit none
    type(c_ptr), intent(in), value :: grid_def_ptr
    type(c_ptr) :: ptr
    !--
    type(grid_quant), pointer :: f_ptr
    type(grid_def), pointer :: grid

    allocate(f_ptr)
    call c_f_pointer(grid_def_ptr, grid)
    call AllocGridQuant(grid, f_ptr%data)
    f_ptr%grid = grid
    ptr = c_loc(f_ptr)
  end function hoppet_cxx__grid_quant__new

  subroutine hoppet_cxx__grid_quant__delete(ptr) bind(C)
    implicit none
    type(c_ptr), intent(inout) :: ptr
    !--
    type(grid_quant), pointer :: f_ptr

    call c_f_pointer(ptr, f_ptr)
    if (associated(f_ptr%data)) then
      deallocate(f_ptr%data)
    end if
    !print * , associated(f_ptr)
    deallocate(f_ptr)
    ptr = c_null_ptr
  end subroutine hoppet_cxx__grid_quant__delete

  function hoppet_cxx__grid_quant__at_y(ptr, y) bind(C) result(res)
    implicit none
    type(c_ptr), intent(in), value :: ptr
    real(c_double), intent(in), value :: y
    real(c_double) :: res
    !--
    type(grid_quant), pointer :: f_ptr

    call c_f_pointer(ptr, f_ptr)
    res = EvalGridQuant(f_ptr%grid, f_ptr%data, y)
  end function hoppet_cxx__grid_quant__at_y

  function hoppet_cxx__grid_quant__data_ptr(ptr) bind(C) result(data_ptr)
    implicit none
    type(c_ptr), intent(in), value :: ptr
    type(c_ptr) :: data_ptr
    !--
    type(grid_quant), pointer :: f_ptr

    call c_f_pointer(ptr, f_ptr)
    data_ptr = c_loc(f_ptr%data(0))
  end function hoppet_cxx__grid_quant__data_ptr


end module hoppet_cxx_oo