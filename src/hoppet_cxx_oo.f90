module hoppet_cxx_oo
  use types, only: dp
  use, intrinsic :: iso_c_binding
  implicit none
  
contains


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

  subroutine hoppet_cxx__grid_def__delete(ptr) bind(C)
    use convolution
    implicit none
    type(c_ptr), intent(inout) :: ptr
    !--
    type(grid_def), pointer :: f_ptr

    call c_f_pointer(ptr, f_ptr)
    print * , associated(f_ptr)
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
end module hoppet_cxx_oo