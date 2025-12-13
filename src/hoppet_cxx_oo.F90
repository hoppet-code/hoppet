
#include "inc/ftlMacros.inc"
#include "inc/extraMacros.inc"


#define BINDC(OBJ,NAME) bind(C,name="hoppet_cxx__"//STRINGIFY(OBJ)//"__"//STRINGIFY(NAME))

! define a macro to generate the functions that returns a generic object's integer member
#define DEFINE_RETURN_INT_MEMBER(OBJ,NAME) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    integer(c_int) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = obj_f%NAME;\
  end function

  ! define a macro to generate the functions that returns a generic object's real(dp) member
#define DEFINE_RETURN_DBL_MEMBER(OBJ,NAME) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    real(c_double) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = obj_f%NAME;\
  end function

  ! define a macro to generate the functions that returns a generic object's real(dp) member
#define DEFINE_RETURN_DBL_ARR_PTR_MEMBER(OBJ,NAME) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    type(c_ptr) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = c_loc(obj_f%NAME(lbound(obj_f%NAME)));\
  end function

! define a macro to generate the functions that returns a generic object's logical member
#define DEFINE_RETURN_LOG_MEMBER(OBJ,NAME) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    logical(c_bool) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = obj_f%NAME;\
  end function

! define a macro to generate the functions that returns a generic object's integer member
#define DEFINE_RETURN_INT_MEMBER_I(OBJ,NAME) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj,i) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    integer(c_int), intent(in), value :: i;\
    integer(c_int) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = obj_f%NAME(i);\
  end function

  ! define a macro to generate the functions that returns a generic object's real(dp) member
#define DEFINE_RETURN_DBL_MEMBER_I(OBJ,NAME) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj,i) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    integer(c_int), intent(in), value :: i;\
    real(c_double) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = obj_f%NAME(i);\
  end function

#define DEFINE_RETURN_OBJ_MEMBER(OBJ,NAME,TYPE) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    type(c_ptr) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = c_loc(obj_f%NAME);\
  end function

#define DEFINE_RETURN_OBJ_MEMBER_I(OBJ,NAME,TYPE) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj,i) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    integer(c_int), intent(in), value :: i;\
    type(c_ptr) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = c_loc(obj_f%NAME(i));\
  end function

#define DEFINE_RETURN_OBJ_MEMBER_IJ(OBJ,NAME,TYPE) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj,i,j) BINDC(OBJ,NAME) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    integer(c_int), intent(in), value :: i,j;\
    type(c_ptr) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = c_loc(obj_f%NAME(i,j));\
  end function


#define DEFINE_DELETE(OBJ) \
  subroutine CAT3(hoppet_cxx__,OBJ,__delete)(ptr) bind(C);\
    implicit none;\
    type(c_ptr), intent(inout) :: ptr;\
    type(OBJ), pointer :: f_ptr;\
    call c_f_pointer(ptr, f_ptr);\
    call delete(f_ptr);\
    deallocate(f_ptr);\
    ptr = c_null_ptr;\
  end subroutine



!=====================================================================
!! misc routines for the C++ OO interface, for example to help with debugging
module hoppet_cxx_oo_various
  use, intrinsic :: iso_c_binding
  implicit none
contains

  !! routine to help testing: print a C-style string passed from C++
  subroutine hoppet_cxx__print_c_str(cstr_ptr) bind(C)
    use hoppet_c_f_string_utils
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(in), value :: cstr_ptr
    character(len=:), allocatable :: fstr
    fstr = fortran_string_from_cstr(cstr_ptr)
    print '(a)', "START:"//fstr//":END"
  end subroutine hoppet_cxx__print_c_str

end module hoppet_cxx_oo_various

!=====================================================================
!! wrappers for C++ OO interface to grid_def
module hoppet_cxx_oo_grid_def
  use types, only: dp
  use, intrinsic :: iso_c_binding
  use convolution
  implicit none  

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

  function hoppet_cxx__grid_def__equiv(grid1, grid2) bind(C) result(res)
    use convolution
    implicit none
    type(c_ptr), intent(in), value :: grid1, grid2
    logical(c_bool) :: res
    type(grid_def), pointer :: grid1_f, grid2_f

    call c_f_pointer(grid1, grid1_f)
    call c_f_pointer(grid2, grid2_f)

    res = (grid1_f == grid2_f)
  end function hoppet_cxx__grid_def__equiv

  subroutine hoppet_cxx__grid_def__monotonic_indices(grid, indices, nindices) bind(C)
    use convolution
    implicit none
    type(c_ptr), intent(in), value :: grid
    integer(c_int), intent(out) :: indices(*)
    integer(c_int), intent(out) :: nindices
    !------
    type(grid_def), pointer :: grid_f
    real(dp), allocatable :: indices_f(:)
    integer :: i

    call c_f_pointer(grid, grid_f)
    indices_f = MonotonicUniqueIndices(grid_f)
    nindices = size(indices_f)
    indices(1:nindices) = indices_f(:)
  end subroutine hoppet_cxx__grid_def__monotonic_indices  

  DEFINE_RETURN_DBL_MEMBER(grid_def,dy)
  DEFINE_RETURN_DBL_MEMBER(grid_def,ymax)
  DEFINE_RETURN_DBL_MEMBER(grid_def,eps)
  DEFINE_RETURN_INT_MEMBER(grid_def,nsub)
  DEFINE_RETURN_INT_MEMBER(grid_def,order)
  DEFINE_RETURN_LOG_MEMBER(grid_def,locked)
  DEFINE_RETURN_OBJ_MEMBER_I(grid_def,subgd,grid_def)
  DEFINE_RETURN_INT_MEMBER_I(grid_def,subiy)

end module hoppet_cxx_oo_grid_def

!=====================================================================
!! wrappers for C++ OO interface to grid quantities, with an intermediate
!! grid_quant type holds both the data and the grid_def
module hoppet_cxx_oo_grid_quant
  use types, only: dp
  use, intrinsic :: iso_c_binding
  use hoppet_cxx_oo_grid_def
  implicit none

  type grid_quant
    real(dp), pointer :: data(:) => null()
    type(grid_def)    :: grid
  end type grid_quant

contains

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

  function hoppet_cxx__grid_quant__at_y(grid_ptr, data_ptr, y) bind(C) result(res)
    implicit none
    type(c_ptr), intent(in), value :: grid_ptr, data_ptr
    real(c_double), intent(in), value :: y
    real(c_double) :: res
    !--
    type(grid_def), pointer :: grid
    real(dp), pointer :: gq(:)

    call c_f_pointer(grid_ptr, grid)
    call c_f_pointer(data_ptr, gq, [grid%ny+1])
    res = EvalGridQuant(grid, gq, y)
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

  function hoppet_cxx__grid_quant__trunc_mom(grid, data, n, ymax) bind(C) result(res)
    implicit none
    type(c_ptr),    intent(in), value :: grid, data  
    real(c_double), intent(in), value :: n
    real(c_double), intent(in)        :: ymax
    real(c_double) :: res
    !--
    type(grid_def),   pointer :: grid_ptr
    real(dp),         pointer :: data_ptr(:)

    call c_f_pointer(grid, grid_ptr)
    call c_f_pointer(data, data_ptr, shape=[grid_ptr%ny+1])
    res = TruncatedMoment(grid_ptr, data_ptr, n, ymax)
  end function hoppet_cxx__grid_quant__trunc_mom

  !! result_data = lumi(gq1,gq2)
  subroutine hoppet_cxx__grid_quant__luminosity(grid, gq1, gq2, result_data) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: grid
    type(c_ptr), intent(in), value :: gq1
    type(c_ptr), intent(in), value :: gq2
    type(c_ptr), intent(in), value :: result_data
    !--
    type(grid_def), pointer :: grid_f
    real(dp), pointer :: gq1_f(:), gq2_f(:)
    real(dp), pointer :: result_data_f(:)

    call c_f_pointer(grid, grid_f)
    call c_f_pointer(gq1, gq1_f,                 shape=[grid_f%ny+1])
    call c_f_pointer(gq2, gq2_f,                 shape=[grid_f%ny+1])
    call c_f_pointer(result_data, result_data_f, shape=[grid_f%ny+1])

    result_data_f = PartonLuminosity(grid_f, gq1_f, gq2_f)
  end subroutine 

end module hoppet_cxx_oo_grid_quant

!=====================================================================
!! wrappers for C++ OO interface to grid quantities, with an intermediate
!! grid_quant type holds both the data and the grid_def
module hoppet_cxx_oo_grid_quant_2d
  use types, only: dp
  use, intrinsic :: iso_c_binding
  use hoppet_cxx_oo_grid_def
  implicit none

  type grid_quant_2d
    real(dp), pointer :: data(:,:) => null()
    type(grid_def)    :: grid
  end type grid_quant_2d

contains

  function hoppet_cxx__grid_quant_2d__new(grid_def_ptr, size) bind(C) result(ptr)
    implicit none
    type(c_ptr),    intent(in), value :: grid_def_ptr
    integer(c_int), intent(in), value :: size
    type(c_ptr) :: ptr
    !--
    type(grid_quant_2d), pointer :: f_ptr
    type(grid_def), pointer :: grid

    allocate(f_ptr)
    call c_f_pointer(grid_def_ptr, grid)
    call AllocGridQuant(grid, f_ptr%data, 0, size-1)
    f_ptr%grid = grid
    ptr = c_loc(f_ptr)
  end function hoppet_cxx__grid_quant_2d__new

  subroutine hoppet_cxx__grid_quant_2d__delete(ptr) bind(C)
    implicit none
    type(c_ptr), intent(inout) :: ptr
    !--
    type(grid_quant_2d), pointer :: f_ptr

    call c_f_pointer(ptr, f_ptr)
    if (associated(f_ptr%data)) then
      deallocate(f_ptr%data)
    end if
    !print * , associated(f_ptr)
    deallocate(f_ptr)
    ptr = c_null_ptr
  end subroutine hoppet_cxx__grid_quant_2d__delete

!  function hoppet_cxx__grid_quant__at_y(grid_ptr, data_ptr, y) bind(C) result(res)
!    implicit none
!    type(c_ptr), intent(in), value :: grid_ptr, data_ptr
!    real(c_double), intent(in), value :: y
!    real(c_double) :: res
!    !--
!    type(grid_def), pointer :: grid
!    real(dp), pointer :: gq(:)
!
!    call c_f_pointer(grid_ptr, grid)
!    call c_f_pointer(data_ptr, gq, [grid%ny+1])
!    res = EvalGridQuant(grid, gq, y)
!  end function hoppet_cxx__grid_quant__at_y

  function hoppet_cxx__grid_quant_2d__data_ptr(ptr) bind(C) result(data_ptr)
    implicit none
    type(c_ptr), intent(in), value :: ptr
    type(c_ptr) :: data_ptr
    !--
    type(grid_quant_2d), pointer :: f_ptr

    call c_f_pointer(ptr, f_ptr)
    data_ptr = c_loc(f_ptr%data(0, 0))
  end function hoppet_cxx__grid_quant_2d__data_ptr

!  function hoppet_cxx__grid_quant__trunc_mom(grid, data, n, ymax) bind(C) result(res)
!    implicit none
!    type(c_ptr),    intent(in), value :: grid, data  
!    real(c_double), intent(in), value :: n
!    real(c_double), intent(in)        :: ymax
!    real(c_double) :: res
!    !--
!    type(grid_def),   pointer :: grid_ptr
!    real(dp),         pointer :: data_ptr(:)
!
!    call c_f_pointer(grid, grid_ptr)
!    call c_f_pointer(data, data_ptr, shape=[grid_ptr%ny+1])
!    res = TruncatedMoment(grid_ptr, data_ptr, n, ymax)
!  end function hoppet_cxx__grid_quant__trunc_mom
end module hoppet_cxx_oo_grid_quant_2d

!=====================================================================
!! wrappers for C++ OO interface to grid_conv objects
module hoppet_cxx_oo_grid_conv
  use types, only: dp
  use, intrinsic :: iso_c_binding
  use hoppet_cxx_oo_grid_def
  use hoppet_cxx_oo_grid_quant
  implicit none


  abstract interface
    function conv_ignd_c_interface(y, piece) bind(C) result(res)
      import :: c_double, c_int
      real(c_double), intent(in), value :: y
      integer(c_int), intent(in), value :: piece
      real(c_double) :: res
    end function conv_ignd_c_interface
  end interface


  type, extends(conv_ignd) :: conv_ignd_from_c
    !type(c_funptr) :: fptr  = c_null_funptr
    type(c_ptr)    :: ctx   = c_null_ptr  
  contains
    procedure :: f => conv_ignd_from_c__f! f(this, y, piece)
  end type conv_ignd_from_c
  interface
    function hoppet_grid_conv_f__wrapper(y, piece, ctx) bind(C)
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in), value :: y
      integer(c_int), intent(in), value :: piece
      type(c_ptr),    intent(in), value :: ctx
      real(c_double) :: hoppet_grid_conv_f__wrapper
    end function hoppet_grid_conv_f__wrapper
  end interface

contains

  !! return a c pointer to a new grid_conv constructed from a conv_ignd_c_interface object
  function hoppet_cxx__grid_conv__new_from_fn(grid_ptr, conv_ignd_c_fn_obj, split_array_ptr, nsplit) bind(C) result(ptr)
    use warnings_and_errors
    implicit none
    type(c_ptr), intent(in), value :: grid_ptr
    type(c_ptr), intent(in), value :: conv_ignd_c_fn_obj
    type(c_ptr), intent(in), value :: split_array_ptr
    integer(c_int), intent(in), optional :: nsplit
    type(c_ptr) :: ptr
    !--
    type(grid_conv), pointer :: gc
    type(grid_def), pointer :: grid
    procedure(conv_ignd_c_interface), pointer :: conv_ignd_f_fn
    type(conv_ignd_from_c) :: lcl_conv_ignd
    real(dp), pointer :: split_array_f(:)

    call c_f_pointer(grid_ptr, grid)
    lcl_conv_ignd%ctx = conv_ignd_c_fn_obj

    allocate(gc)
    if (c_associated(split_array_ptr)) then
      if (.not. present(nsplit)) then
        call wae_error("hoppet_cxx__grid_conv__new_from_fn: nsplit must be provided if split_array_ptr is associated")
      end if
      call c_f_pointer(split_array_ptr, split_array_f, shape=[nsplit])
      call InitGridConv(grid, gc, lcl_conv_ignd, alloc=.true., split=split_array_f)
    else
      call InitGridConv(grid, gc, lcl_conv_ignd, alloc=.true.)
    end if
    call InitGridConv(grid, gc, lcl_conv_ignd, alloc=.true.)
    ptr = c_loc(gc)
  end function hoppet_cxx__grid_conv__new_from_fn

  !! return a c pointer to a new grid_conv that is a copy of gc_other
  function hoppet_cxx__grid_conv__new_from_gc(gc_other) bind(C) result(ptr)
    implicit none
    type(c_ptr), intent(in), value :: gc_other
    type(c_ptr) :: ptr
    !--
    type(grid_conv), pointer :: gc
    type(grid_conv), pointer :: gc_other_f

    call c_f_pointer(gc_other, gc_other_f)
    allocate(gc)
    call InitGridConv(gc, gc_other_f, alloc=.true.)
    ptr = c_loc(gc)
  end function hoppet_cxx__grid_conv__new_from_gc

  !! copy the contents of gc_src into gc_dest
  subroutine hoppet_cxx__grid_conv__copy_contents(dest, src) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: dest, src
    !--
    type(grid_conv), pointer :: dest_f, src_f

    call c_f_pointer(dest, dest_f)
    call c_f_pointer(src, src_f)

    call InitGridConv(dest_f, src_f)
  end subroutine hoppet_cxx__grid_conv__copy_contents


  !! implementation of conv_ignd_from_c%f
  function conv_ignd_from_c__f(this, y, piece) result(res)
    class(conv_ignd_from_c), intent(in) :: this
    real(dp), intent(in) :: y
    integer,  intent(in) :: piece
    real(dp) :: res

    res = hoppet_grid_conv_f__wrapper(y, piece, this%ctx)
  end function conv_ignd_from_c__f

  !! delete a grid_conv object (and the the associated fortran object)
  subroutine hoppet_cxx__grid_conv__delete(ptr) bind(C)
    implicit none
    type(c_ptr), intent(inout) :: ptr
    !--
    type(grid_conv), pointer :: gc

    call c_f_pointer(ptr, gc)
    call delete(gc)
    deallocate(gc)
    ptr = c_null_ptr
  end subroutine hoppet_cxx__grid_conv__delete

  !! result_data = conv * q_data
  subroutine hoppet_cxx__grid_conv__times_grid_quant(conv_ptr, q_data, result_data) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: conv_ptr
    type(c_ptr), intent(in), value :: q_data
    type(c_ptr), intent(in), value :: result_data
    !--
    type(grid_conv), pointer :: gc
    real(dp), pointer :: q(:)
    real(dp), pointer :: result(:)

    call c_f_pointer(conv_ptr, gc)
    call c_f_pointer(q_data, q, shape=[gc%grid%ny+1])
    call c_f_pointer(result_data, result, shape=[gc%grid%ny+1])

    result = gc * q
  end subroutine hoppet_cxx__grid_conv__times_grid_quant


  !! conv1 += conv2
  subroutine hoppet_cxx__grid_conv__add(conv1, conv2, factor) bind(C)
    implicit none
    type(c_ptr),    intent(in), value    :: conv1, conv2
    real(c_double), intent(in), optional :: factor
    type(grid_conv), pointer :: conv1_f, conv2_f
    call c_f_pointer(conv1, conv1_f)
    call c_f_pointer(conv2, conv2_f)
    call AddWithCoeff(conv1_f,conv2_f,fact=factor)
  end subroutine hoppet_cxx__grid_conv__add

  !! conv1 *=factor
  subroutine hoppet_cxx__grid_conv__multiply(conv1, factor) bind(C)
    implicit none
    type(c_ptr),    intent(in), value :: conv1
    real(c_double), intent(in), value :: factor
    type(grid_conv), pointer :: conv1_f
    call c_f_pointer(conv1, conv1_f)
    call Multiply(conv1_f,factor)
  end subroutine hoppet_cxx__grid_conv__multiply

  !! allocate a new grid_conv that is the convolution of conv1 and conv2
  function hoppet_cxx__grid_conv__alloc_and_conv(conv1, conv2) bind(C) result(res)
    implicit none
    type(c_ptr),    intent(in), value :: conv1, conv2
    type(c_ptr) :: res
    type(grid_conv), pointer :: conv1_f, conv2_f, res_f

    call c_f_pointer(conv1, conv1_f)
    call c_f_pointer(conv2, conv2_f)
    allocate(res_f)
    call AllocGridConv(conv1_f%grid, res_f)
    call SetToConvolution(res_f, conv1_f, conv2_f)
    res = c_loc(res_f)
  end function hoppet_cxx__grid_conv__alloc_and_conv

  function hoppet_cxx__grid_conv__moment(conv, momN) bind(C) result(res)
    implicit none
    type(c_ptr),    intent(in), value :: conv
    real(c_double), intent(in), value :: momN
    real(c_double) :: res
    type(grid_conv), pointer :: conv_f

    call c_f_pointer(conv, conv_f)
    res = gc_moment(conv_f, momN)
  end function hoppet_cxx__grid_conv__moment

  !! return a pointer to the convolution object inside the grid_conv
  !! or a null pointer if it not allocated
  function hoppet_cxx__grid_conv__conv_ptr(gc, sz_y, sz_order) bind(C) result(res)
    implicit none
    type(c_ptr),    intent(in), value :: gc
    integer(c_int), intent(out)       :: sz_y
    integer(c_int), intent(out)       :: sz_order
    type(c_ptr)                       :: res
    type(grid_conv), pointer :: gc_f

    call c_f_pointer(gc, gc_f)
    if (associated(gc_f%conv)) then
      res = c_loc(gc_f%conv(lbound(gc_f%conv,1), lbound(gc_f%conv,2)))
      sz_y = size(gc_f%conv,1)
      sz_order = size(gc_f%conv,2)
    else
      res = c_null_ptr
    end if
  end function hoppet_cxx__grid_conv__conv_ptr

  !! finally the member accessors
  DEFINE_RETURN_OBJ_MEMBER(grid_conv,grid,grid_def)
  DEFINE_RETURN_OBJ_MEMBER_I(grid_conv,subgc,grid_conv)

end module hoppet_cxx_oo_grid_conv


!=====================================================================
!! wrappers for C++ OO interface to split_mat objects
module hoppet_cxx_oo_split_mat
  use types, only: dp
  use, intrinsic :: iso_c_binding
  use hoppet_cxx_oo_grid_def
  use hoppet_cxx_oo_grid_quant
  use dglap_objects
  implicit none

contains
  !! return a c pointer to a new split_mat with the given lcl_nf
  function hoppet_cxx__split_mat__new(lcl_nf) bind(C) result(split_mat_ptr)
    implicit none
    integer(c_int), intent(in), value :: lcl_nf
    type(c_ptr) :: split_mat_ptr
    !--
    type(split_mat), pointer :: f_ptr

    allocate(f_ptr)
    call cobj_InitSplitLinks(f_ptr)
    f_ptr%nf_int = lcl_nf
    split_mat_ptr = c_loc(f_ptr)
  end function hoppet_cxx__split_mat__new

  !! return a c pointer to a copy of the (split_mat_c_ptr) split_mat
  function hoppet_cxx__split_mat__copy(other) bind(C) result(split_mat_ptr)
    implicit none
    type(c_ptr), intent(in), value :: other
    type(c_ptr) :: split_mat_ptr
    !--
    type(split_mat), pointer :: f_ptr, other_f_ptr

    call c_f_pointer(other, other_f_ptr)
    allocate(f_ptr)
    call InitSplitMat(f_ptr, other_f_ptr)
    split_mat_ptr = c_loc(f_ptr)
  end function hoppet_cxx__split_mat__copy

  !! copy the contents of the split mat 'src' into 'dest'
  subroutine hoppet_cxx__split_mat__copy_contents(dest, src) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: dest
    type(c_ptr), intent(in), value :: src
    !--
    type(split_mat), pointer :: dest_f, src_f

    call c_f_pointer(dest, dest_f)
    call c_f_pointer(src, src_f)
    call InitSplitMat(dest_f, src_f)
  end subroutine hoppet_cxx__split_mat__copy_contents

  !! sm1 += sm2
  subroutine hoppet_cxx__split_mat__add(sm1, sm2, factor) bind(C)
    implicit none
    type(c_ptr),    intent(in), value    :: sm1, sm2
    real(c_double), intent(in), optional :: factor
    type(split_mat), pointer :: sm1_f, sm2_f
    call c_f_pointer(sm1, sm1_f)
    call c_f_pointer(sm2, sm2_f)
    call AddWithCoeff(sm1_f,sm2_f,factor)
  end subroutine hoppet_cxx__split_mat__add

  !! sm1 *=factor
  subroutine hoppet_cxx__split_mat__multiply(sm1, factor) bind(C)
    implicit none
    type(c_ptr),    intent(in), value :: sm1
    real(c_double), intent(in), value :: factor
    type(split_mat), pointer :: sm1_f
    call c_f_pointer(sm1, sm1_f)
    call Multiply(sm1_f,factor)
  end subroutine hoppet_cxx__split_mat__multiply

  !! allocate a new split_mat that is the convolution of sm1 and sm2 and return a pointer to it
  function hoppet_cxx__split_mat__alloc_and_conv(sm1, sm2) bind(C) result(res)
    implicit none
    type(c_ptr),    intent(in), value :: sm1, sm2
    type(c_ptr) :: res
    type(split_mat), pointer :: sm1_f, sm2_f, res_f

    call c_f_pointer(sm1, sm1_f)
    call c_f_pointer(sm2, sm2_f)
    allocate(res_f)
    call AllocSplitMat(sm1_f%qq%grid, res_f, sm1_f%nf_int)
    call SetToConvolution(res_f, sm1_f, sm2_f)
    res = c_loc(res_f)
  end function hoppet_cxx__split_mat__alloc_and_conv

  !! allocate a new split_mat that is the convolution of sm1 and sm2 and return a pointer to it
  function hoppet_cxx__split_mat__alloc_and_commutate(sm1, sm2) bind(C) result(res)
    implicit none
    type(c_ptr),    intent(in), value :: sm1, sm2
    type(c_ptr) :: res
    type(split_mat), pointer :: sm1_f, sm2_f, res_f

    call c_f_pointer(sm1, sm1_f)
    call c_f_pointer(sm2, sm2_f)
    allocate(res_f)
    call AllocSplitMat(sm1_f%qq%grid, res_f, sm1_f%nf_int)
    call SetToCommutator(res_f, sm1_f, sm2_f)
    res = c_loc(res_f)
  end function hoppet_cxx__split_mat__alloc_and_commutate


  !! result_data = split_mat_ptr * q_data
  subroutine hoppet_cxx__split_mat__times_grid_quant_2d(split_mat_ptr, q_data, result_data) bind(C)
    use pdf_representation
    implicit none
    type(c_ptr), intent(in), value :: split_mat_ptr
    type(c_ptr), intent(in), value :: q_data
    type(c_ptr), intent(in), value :: result_data
    !--
    type(split_mat), pointer :: split_mat_f
    real(dp), pointer :: q(:,:)
    real(dp), pointer :: result(:,:)

    call c_f_pointer(split_mat_ptr, split_mat_f)
    call c_f_pointer(q_data, q, shape=[split_mat_f%qq%grid%ny+1, ncompmax-ncompmin+1])
    call c_f_pointer(result_data, result, shape=[split_mat_f%qq%grid%ny+1, ncompmax-ncompmin+1])

    result = split_mat_f * q
  end subroutine hoppet_cxx__split_mat__times_grid_quant_2d

  !! return the nf value of the split_mat
  function hoppet_cxx__split_mat__nf(other) bind(C) result(lcl_nf)
    implicit none
    type(c_ptr), intent(in), value :: other
    integer(c_int) :: lcl_nf
    !--
    type(split_mat), pointer :: other_f_ptr

    call c_f_pointer(other, other_f_ptr)
    lcl_nf = other_f_ptr%nf_int
  end function hoppet_cxx__split_mat__nf

  !! delete a split_mat object (and its associated fortran storage)
  subroutine hoppet_cxx__split_mat__delete(split_mat_c_ptr) bind(C)
    implicit none
    type(c_ptr), intent(inout) :: split_mat_c_ptr
    !--
    type(split_mat), pointer :: f_ptr

    call c_f_pointer(split_mat_c_ptr, f_ptr)
    call Delete(f_ptr)
    deallocate(f_ptr)
    split_mat_c_ptr = c_null_ptr
  end subroutine hoppet_cxx__split_mat__delete


  ! define a macro to generate the functions that return references to the grid_conv members
#define SPLIT_MAT__REF(NAME)  DEFINE_RETURN_OBJ_MEMBER(split_mat,NAME,grid_conv)

  SPLIT_MAT__REF(qq      )
  SPLIT_MAT__REF(qg      )
  SPLIT_MAT__REF(gq      )
  SPLIT_MAT__REF(gg      )
  SPLIT_MAT__REF(ns_plus )
  SPLIT_MAT__REF(ns_minus)
  SPLIT_MAT__REF(ns_v    )

#undef SPLIT_MAT__REF

end module hoppet_cxx_oo_split_mat


!=====================================================================
!! wrappers for C++ OO interface to mass_threshold_mat objects
module hoppet_cxx_oo_mass_threshold_mat
  use, intrinsic :: iso_c_binding
  use dglap_objects
  implicit none
contains

  !! return a C pointer to a new mass_threshold_mat with the given lcl_nf
  function hoppet_cxx__mass_threshold_mat__new(lcl_nf) bind(C) result(mass_threshold_mat_ptr)
    implicit none
    integer(c_int), intent(in), value :: lcl_nf
    type(c_ptr) :: mass_threshold_mat_ptr
    !--
    type(mass_threshold_mat), pointer :: f_ptr

    allocate(f_ptr)    
    f_ptr%nf_int = lcl_nf
    f_ptr%loops  = 0
    mass_threshold_mat_ptr = c_loc(f_ptr)
  end function hoppet_cxx__mass_threshold_mat__new

  !! return a C pointer to a copy of the other mass_threshold_mat
  function hoppet_cxx__mass_threshold_mat__copy(other) bind(C) result(mass_threshold_mat_ptr)
    implicit none
    type(c_ptr), intent(in), value :: other
    type(c_ptr) :: mass_threshold_mat_ptr
    !--
    type(mass_threshold_mat), pointer :: f_ptr, other_f_ptr

    call c_f_pointer(other, other_f_ptr)
    allocate(f_ptr)
    call InitMTM(f_ptr, other_f_ptr)
    mass_threshold_mat_ptr = c_loc(f_ptr)
  end function hoppet_cxx__mass_threshold_mat__copy

  !! return a C pointer to a copy of the other mass_threshold_mat
  subroutine hoppet_cxx__mass_threshold_mat__copy_contents(dest,other) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: dest,other
    type(c_ptr) :: mass_threshold_mat_ptr
    !--
    type(mass_threshold_mat), pointer :: f_ptr, other_f_ptr

    call c_f_pointer(other, other_f_ptr)
    call InitMTM(other_f_ptr, f_ptr)
  end subroutine hoppet_cxx__mass_threshold_mat__copy_contents

  DEFINE_DELETE(mass_threshold_mat)

  !! set the nf(heavy) value of the mass_threshold_mat
  subroutine hoppet_cxx__mass_threshold_mat__set_nf(mtm, nf_lcl) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: mtm
    integer(c_int), intent(in), value :: nf_lcl
    !--
    type(mass_threshold_mat), pointer :: mtm_f

    call c_f_pointer(mtm, mtm_f)
    call SetNfMTM(mtm_f, nf_lcl)
  end subroutine hoppet_cxx__mass_threshold_mat__set_nf

  !! result_data = mass_threshold_mat_ptr * q_data
  subroutine hoppet_cxx__mass_threshold_mat__times_grid_quant_2d(mass_threshold_mat_ptr, q_data, result_data) bind(C)
    use pdf_representation
    implicit none
    type(c_ptr), intent(in), value :: mass_threshold_mat_ptr
    type(c_ptr), intent(in), value :: q_data
    type(c_ptr), intent(in), value :: result_data
    !--
    type(mass_threshold_mat), pointer :: mass_threshold_mat_f
    real(dp), pointer :: q(:,:)
    real(dp), pointer :: result(:,:)

    call c_f_pointer(mass_threshold_mat_ptr, mass_threshold_mat_f)
    call c_f_pointer(q_data, q, shape=[mass_threshold_mat_f%PShq%grid%ny+1, ncompmax-ncompmin+1])
    call c_f_pointer(result_data, result, shape=[mass_threshold_mat_f%PShq%grid%ny+1, ncompmax-ncompmin+1])

    result = mass_threshold_mat_f * q
  end subroutine hoppet_cxx__mass_threshold_mat__times_grid_quant_2d

  !! mtm1 += mtm2
  subroutine hoppet_cxx__mass_threshold_mat__add(mtm1, mtm2, factor) bind(C)
    implicit none
    type(c_ptr),    intent(in), value    :: mtm1, mtm2
    real(c_double), intent(in), optional :: factor
    type(mass_threshold_mat), pointer :: mtm1_f, mtm2_f
    call c_f_pointer(mtm1, mtm1_f)
    call c_f_pointer(mtm2, mtm2_f)
    call AddWithCoeff(mtm1_f,mtm2_f,factor)
  end subroutine hoppet_cxx__mass_threshold_mat__add

  !! mtm1 *=factor
  subroutine hoppet_cxx__mass_threshold_mat__multiply(mtm1, factor) bind(C)
    implicit none
    type(c_ptr),    intent(in), value :: mtm1
    real(c_double), intent(in), value :: factor
    type(mass_threshold_mat), pointer :: mtm1_f
    call c_f_pointer(mtm1, mtm1_f)
    call Multiply(mtm1_f,factor)
  end subroutine hoppet_cxx__mass_threshold_mat__multiply


#define MTM_REF(NAME)  DEFINE_RETURN_OBJ_MEMBER(mass_threshold_mat,NAME,grid_conv)
  MTM_REF(pshq      ) !!< A^PS_Qq    Q+Qbar from singlet(nflight)
  MTM_REF(pshg      ) !!< A^PS_Qg    Q+Qbar from gluon  (nflight)
  MTM_REF(nsqq_h    ) !!< A^NS_qq,Q  ΔNS(nfheavy) from NS(nflight)
  MTM_REF(sgg_h     ) !!< A^S_gg,Q   Δg(nfheavy) from g(nflight)
  MTM_REF(sgq_H     ) !!< A^S_gq,Q   Δg(nfheavy) from singlet(nflight)
  MTM_REF(psqq_h    ) !!< A^PS_qq,Q  Δsinglet(nfheavy) from singlet(nflight)
  MTM_REF(sqg_h     ) !!< A^S_qg,Q   Δsinglet(nfheavy) from gluon(nflight)
  MTM_REF(nsmqq_h   ) !!< A^{NSm}_qq,Q ΔNSminus(1:nflight) from NSminus(1:nflight)
  MTM_REF(pshg_msbar) !!< replaces PShg when masses are MSbar

#define MTM_INT(NAME)  DEFINE_RETURN_INT_MEMBER(mass_threshold_mat,NAME)
  MTM_INT(nf_int)
  MTM_INT(loops)

end module hoppet_cxx_oo_mass_threshold_mat

!=====================================================================
!! module for the wrapping a running_coupling object
module hoppet_cxx_oo_running_coupling
  use, intrinsic :: iso_c_binding
  use qcd_coupling
  implicit none
contains

  function hoppet_cxx__running_coupling__new_fixnf(alphas_at_Q, Q, nloop, fixnf, Qmax) bind(C) result(res)
    real(c_double), intent(in), value :: alphas_at_Q, Q
    integer(c_int), intent(in), value :: nloop, fixnf
    real(c_double), intent(in), optional :: Qmax
    type(c_ptr) :: res
    !--
    type(running_coupling), pointer :: rc_f
    allocate(rc_f)
    call InitRunningCoupling(rc_f, alphas_at_Q, Q, nloop, fixnf=fixnf, Qmax=Qmax)
    res = c_loc(rc_f)
  end function hoppet_cxx__running_coupling__new_fixnf

  function hoppet_cxx__running_coupling__new_varnf(alphas_at_Q, Q, nloop, mc, mb, mt, masses_are_MSbar, muMatch_mQuark, Qmax) bind(C) result(res)
    real(c_double), intent(in), value :: alphas_at_Q, Q
    integer(c_int), intent(in), value :: nloop
    real(c_double), intent(in), value :: mc, mb, mt
    logical(c_bool), intent(in), value :: masses_are_MSbar
    real(c_double), intent(in), value :: muMatch_mQuark
    real(c_double), intent(in), optional :: Qmax
    type(c_ptr) :: res
    !--
    type(running_coupling), pointer :: rc_f
    logical :: masses_are_MSbar_f
    masses_are_MSbar_f = (masses_are_MSbar .neqv. .false._c_bool) ! safely convert to fortran logical
    allocate(rc_f)
    call InitRunningCoupling(rc_f, alphas_at_Q, Q, nloop, quark_masses = [mc, mb, mt], &
                             masses_are_MSbar = masses_are_MSbar_f, muMatch_mQuark = muMatch_mQuark, Qmax=Qmax)
    res = c_loc(rc_f)
  end function hoppet_cxx__running_coupling__new_varnf

  subroutine hoppet_cxx__running_coupling__delete(rc) bind(C)
    implicit none
    type(c_ptr), intent(inout) :: rc
    !--
    type(running_coupling), pointer :: f_ptr
    call c_f_pointer(rc, f_ptr)
    call Delete(f_ptr)
    deallocate(f_ptr)
    rc = c_null_ptr
  end subroutine hoppet_cxx__running_coupling__delete


  function hoppet_cxx__running_coupling__value(rc, Q, fixnf) bind(C) result(res)
    implicit none
    type(c_ptr),    intent(in), value    :: rc
    real(c_double), intent(in), value    :: Q
    integer(c_int), intent(in), optional :: fixnf
    real(c_double) :: res
    !--
    type(running_coupling), pointer :: rc_f

    call c_f_pointer(rc, rc_f)
    res = Value(rc_f, Q, fixnf=fixnf)
  end function hoppet_cxx__running_coupling__value

  !-- number of loops in the coupling
  function hoppet_cxx__running_coupling__num_loops(rc) bind(C) result(nloops)
    implicit none
    type(c_ptr), intent(in), value :: rc
    integer(c_int) :: nloops
    type(running_coupling), pointer :: rc_f

    call c_f_pointer(rc, rc_f)
    nloops = NumberOfLoops(rc_f)
  end function hoppet_cxx__running_coupling__num_loops

  !-- nf range: returns nflo,nfhi via out arguments
  subroutine hoppet_cxx__running_coupling__nf_range(rc, nflo, nfhi) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: rc
    integer(c_int), intent(out) :: nflo, nfhi
    type(running_coupling), pointer :: rc_f

    call c_f_pointer(rc, rc_f)
    call NfRange(rc_f, nflo, nfhi)
  end subroutine hoppet_cxx__running_coupling__nf_range

  !-- nf at Q: returns nf and optional Qlo,Qhi via out args
  function hoppet_cxx__running_coupling__nf_at_q(rc, Q, Qlo, Qhi, muM_mQ) bind(C) result(nf)
    implicit none
    type(c_ptr), intent(in), value :: rc
    real(c_double), intent(in), value :: Q
    real(c_double), intent(out), optional :: Qlo, Qhi
    real(c_double), intent(in),  optional :: muM_mQ
    integer(c_int) :: nf
    type(running_coupling), pointer :: rc_f

    call c_f_pointer(rc, rc_f)
    nf = NfAtQ(rc_f, Q, Qlo, Qhi, muM_mQ=muM_mQ)
  end function hoppet_cxx__running_coupling__nf_at_q

  !-- quark mass for a flavour
  function hoppet_cxx__running_coupling__quark_mass(rc, iflv) bind(C) result(mass)
    implicit none
    type(c_ptr), intent(in), value :: rc
    integer(c_int), intent(in), value :: iflv
    real(c_double) :: mass
    type(running_coupling), pointer :: rc_f

    call c_f_pointer(rc, rc_f)
    mass = QuarkMass(rc_f, iflv)
  end function hoppet_cxx__running_coupling__quark_mass

  !-- quark masses are MSbar?
  function hoppet_cxx__running_coupling__quark_masses_are_msbar(rc) bind(C) result(res)
    implicit none
    type(c_ptr), intent(in), value :: rc
    logical(c_bool) :: res
    type(running_coupling), pointer :: rc_f

    call c_f_pointer(rc, rc_f)
    res = QuarkMassesAreMSbar(rc_f)
  end function hoppet_cxx__running_coupling__quark_masses_are_msbar

  !-- Q range at nf
  subroutine hoppet_cxx__running_coupling__q_range_at_nf(rc, nflcl, Qlo, Qhi) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: rc
    integer(c_int), intent(in), value :: nflcl
    real(c_double), intent(out) :: Qlo, Qhi
    type(running_coupling), pointer :: rc_f

    call c_f_pointer(rc, rc_f)
    call QRangeAtNf(rc_f, nflcl, Qlo, Qhi)
  end subroutine hoppet_cxx__running_coupling__q_range_at_nf

end module hoppet_cxx_oo_running_coupling

!=====================================================================
! wrappers for C++ OO interface to dglap_holder
module hoppet_cxx_oo_dglap_holder
  use, intrinsic :: iso_c_binding
  use dglap_holders
  use convolution
  implicit none
contains

  function hoppet_cxx__dglap_holder__new(grid, factscheme, nloop, nflo, nfhi) bind(C) result(dh_ptr)
    implicit none
    type(c_ptr) :: dh_ptr
    type(c_ptr), intent(in), value :: grid
    integer(c_int), intent(in), optional :: factscheme
    integer(c_int), intent(in), optional :: nloop
    integer(c_int), intent(in), optional :: nflo
    integer(c_int), intent(in), optional :: nfhi
    !-- local Fortran pointers
    type(grid_def), pointer :: grid_f
    type(dglap_holder), pointer :: dh_f

    call c_f_pointer(grid, grid_f)
    allocate(dh_f)

    call InitDglapHolder(grid_f, dh_f, factscheme, nloop, nflo, nfhi)

    dh_ptr = c_loc(dh_f)
  end function hoppet_cxx__dglap_holder__new

  subroutine hoppet_cxx__dglap_holder__set_nf(dh, nf_lcl) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: dh
    integer(c_int), intent(in), value :: nf_lcl
    !--
    type(dglap_holder), pointer :: dh_f

    call c_f_pointer(dh, dh_f)
    call SetNfDglapHolder(dh_f, nf_lcl)
  end subroutine hoppet_cxx__dglap_holder__set_nf


  DEFINE_DELETE(dglap_holder)
  DEFINE_RETURN_INT_MEMBER(dglap_holder,nloop)
  DEFINE_RETURN_INT_MEMBER(dglap_holder,nf)
  DEFINE_RETURN_INT_MEMBER(dglap_holder,factscheme)
  DEFINE_RETURN_OBJ_MEMBER_IJ(dglap_holder,allp,split_mat)

  DEFINE_RETURN_OBJ_MEMBER(dglap_holder,grid,grid_def)
  DEFINE_RETURN_OBJ_MEMBER(dglap_holder,p_lo,split_mat)
  DEFINE_RETURN_OBJ_MEMBER(dglap_holder,p_nlo,split_mat)
  DEFINE_RETURN_OBJ_MEMBER(dglap_holder,p_nnlo,split_mat)
  DEFINE_RETURN_OBJ_MEMBER(dglap_holder,p_n3lo,split_mat)

  DEFINE_RETURN_OBJ_MEMBER_IJ(dglap_holder,allmtm,mass_threshold_mat)
  DEFINE_RETURN_OBJ_MEMBER(dglap_holder,mtm_nnlo,mass_threshold_mat)
  DEFINE_RETURN_OBJ_MEMBER(dglap_holder,mtm_n3lo,mass_threshold_mat)

end module hoppet_cxx_oo_dglap_holder


!=====================================================================
!! wrappers for C++ OO interface to pdf_table, as well as pdfseginfo
module hoppet_cxx_oo_pdf_table
  use, intrinsic :: iso_c_binding
  use pdf_tabulate
  use convolution
  use dglap_objects
  use types, only: dp
  implicit none

contains

  ! things related to the pdfseginfo type, will be 
  ! accessible only as a "view" so not need for allocation
  ! or deletion functions
  DEFINE_RETURN_INT_MEMBER(pdfseginfo,ilnlnQ_lo)
  DEFINE_RETURN_INT_MEMBER(pdfseginfo,ilnlnQ_hi)
  DEFINE_RETURN_DBL_MEMBER(pdfseginfo,lnlnQ_lo)
  DEFINE_RETURN_DBL_MEMBER(pdfseginfo,lnlnQ_hi)
  DEFINE_RETURN_DBL_MEMBER(pdfseginfo,dlnlnQ)
  DEFINE_RETURN_DBL_MEMBER(pdfseginfo,inv_dlnlnQ)

  !! allocate a new pdf_table object and return a pointer to it
  function hoppet_cxx__pdf_table__new(grid, Qmin, Qmax, dlnlnQ, lnlnQ_order, freeze_at_Qmin, iflv_max_table)  bind(C) result(tab_ptr)
    use pdf_representation
    implicit none
    type(c_ptr) :: tab_ptr
    type(c_ptr), intent(in), value :: grid
    real(c_double),  intent(in), value :: Qmin, Qmax 
    real(c_double),  intent(in), value :: dlnlnQ
    integer(c_int),  intent(in), value :: lnlnQ_order
    logical(c_bool), intent(in), value :: freeze_at_Qmin
    integer(c_int),  intent(in), value :: iflv_max_table
    !-- local Fortran pointers
    type(grid_def), pointer :: grid_f
    type(pdf_table), pointer :: tab_f
    type(logical) :: freeze_at_Qmin_f

    call c_f_pointer(grid, grid_f)
    allocate(tab_f)
    freeze_at_Qmin_f = (freeze_at_Qmin .neqv. .false._c_bool) ! safely convert to fortran logical
    call AllocPDFTable(grid_f, tab_f, Qmin, Qmax, dlnlnQ, lnlnQ_order, &
                      freeze_at_Qmin=freeze_at_Qmin_f, &
                      iflv_max_table=iflv_max_table + iflv_min)

    tab_ptr = c_loc(tab_f)
  end function hoppet_cxx__pdf_table__new

  function hoppet_cxx__pdf_table__copy(other) bind(C) result(tab_ptr)
    implicit none
    type(c_ptr), intent(in), value :: other
    type(c_ptr) :: tab_ptr
    !-- local Fortran pointers
    type(pdf_table), pointer :: other_f, tab_f

    call c_f_pointer(other, other_f)
    allocate(tab_f)
    call AllocPDFTable(tab_f, other_f)
    tab_f%tab = other_f%tab
    tab_ptr = c_loc(tab_f)
  end function hoppet_cxx__pdf_table__copy

  !! Note that as of 2025-11-12 this does not copy any evolution
  !! operators
  subroutine hoppet_cxx__pdf_table__copy_contents(dest, src) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: dest, src
    !-- local Fortran pointers
    type(pdf_table), pointer :: src_f, dest_f

    call c_f_pointer(src, src_f)
    call c_f_pointer(dest, dest_f)
    call AllocPDFTable(dest_f, src_f)
    dest_f%tab = src_f%tab
  end subroutine hoppet_cxx__pdf_table__copy_contents

  subroutine hoppet_cxx__pdf_table__add_nf_info(tab, coupling) bind(C)
    use qcd_coupling
    implicit none
    type(c_ptr), intent(in), value :: tab
    type(c_ptr), intent(in), value :: coupling
    !--
    type(pdf_table), pointer :: tab_f
    type(running_coupling), pointer :: coupling_f

    call c_f_pointer(tab, tab_f)
    call c_f_pointer(coupling, coupling_f)

    call AddNfInfoToPDFTable(tab_f, coupling_f)
  end subroutine hoppet_cxx__pdf_table__add_nf_info

  !! return a c pointer to the first element of the tab array
  function hoppet_cxx__pdf_table__tab_ptr(obj) result(res) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: obj
    type(c_ptr) :: res
    type(pdf_table), pointer :: obj_f

    call c_f_pointer(obj, obj_f)
    res = c_loc(obj_f%tab(lbound(obj_f%tab,1),lbound(obj_f%tab,2),lbound(obj_f%tab,3)))
  end function

  !! return the size of the flavour dimension of the tab array (including the "info" index)
  function hoppet_cxx__pdf_table__size_flv(obj) result(res) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: obj
    integer(c_int) :: res
    type(pdf_table), pointer :: obj_f

    call c_f_pointer(obj, obj_f)
    res = size(obj_f%tab,2)
  end function

  subroutine hoppet_cxx__pdf_table__at_Q_into(tab, Q, pdf) bind(C)
    implicit none
    type(c_ptr),    intent(in), value :: tab
    real(c_double), intent(in), value :: Q
    type(c_ptr),    intent(in), value :: pdf
    !--------
    type(pdf_table), pointer :: tab_f
    real(c_double),  pointer :: pdf_f(:,:)

    call c_f_pointer(tab, tab_f)
    !                                           y             flv
    call c_f_pointer(pdf, pdf_f, shape=[size(tab_f%tab,1), size(tab_f%tab,2)])
    call EvalPdfTable_Q(tab_f, Q, pdf_f)
  end subroutine


  !! return the interpolation of the pdf table at (y,Q,iflv)
  !! watch out iflv_cxx starts from zero, while iflv in fortran pdf_table starts from ncompmin
  function hoppet_cxx__pdf_table__at_yQf(obj, y, Q, iflv_cxx) result(res) bind(C)
    use pdf_representation
    implicit none
    type(c_ptr), intent(in), value :: obj
    real(c_double), intent(in), value :: y, Q
    integer(c_int), intent(in), value :: iflv_cxx
    real(c_double) :: res
    type(pdf_table), pointer :: obj_f

    call c_f_pointer(obj, obj_f)
    res = EvalPdfTable_yQf(obj_f, y, Q, iflv_cxx + ncompmin)
  end function

  !! fill res with the interpolation of the pdf table at (y,Q) for all iflv up to and including iflv_max
  subroutine hoppet_cxx__pdf_table__at_yQ_into(obj, y, Q, res) bind(C)
    use pdf_representation
    implicit none
    type(c_ptr), intent(in), value :: obj
    real(c_double), intent(in), value :: y, Q
    real(c_double), intent(out) :: res(*)
    type(pdf_table), pointer :: obj_f
    integer, parameter :: iflv_max_offset = - ncompmin + 1

    call c_f_pointer(obj, obj_f)
    call EvalPdfTable_yQ(obj_f, y, Q, res( 1:obj_f%tab_iflv_max+iflv_max_offset ) )
  end subroutine

  !! carry out the evolution
  subroutine hoppet_cxx__pdf_table__evolve(tab, &
              Q0, pdf_at_Q0, &
              dh, coupling, muR_over_Q, nloop, untie_nf) bind(C)
    use, intrinsic :: iso_c_binding
    use qcd_coupling,    only: running_coupling
    use dglap_holders,   only: dglap_holder
    implicit none
    type(c_ptr), intent(in), value :: tab
    real(c_double), intent(in), value :: Q0
    type(c_ptr), intent(in), value :: pdf_at_Q0
    type(c_ptr), intent(in), value :: dh, coupling
    real(c_double), intent(in), value :: muR_over_Q
    integer(c_int), intent(in), value :: nloop
    logical(c_bool), intent(in), value :: untie_nf
    !-----------
    real(c_double), pointer :: pdf_at_Q0_f(:,:)
    type(pdf_table), pointer        :: tab_f
    type(running_coupling), pointer :: coupling_f
    type(dglap_holder), pointer    :: dh_f
    logical :: untie_nf_f

    call c_f_pointer(tab, tab_f)
    call c_f_pointer(pdf_at_Q0, pdf_at_Q0_f, shape=[size(tab_f%tab,1), size(tab_f%tab,2)])
    call c_f_pointer(coupling, coupling_f)
    call c_f_pointer(dh, dh_f)
    untie_nf_f = (untie_nf .neqv. .false._c_bool) ! safely convert to fortran logical

    call EvolvePdfTable(tab_f, &
                        Q0, pdf_at_Q0_f, &
                        dh_f, coupling_f, muR_over_Q, nloop, untie_nf=untie_nf_f)
  end subroutine hoppet_cxx__pdf_table__evolve

  !! pre-evolve: set up evolution operators but do not apply them
  subroutine hoppet_cxx__pdf_table__pre_evolve(tab, &
              Q0, &
              dh, coupling, muR_over_Q, nloop, untie_nf) bind(C)
    use, intrinsic :: iso_c_binding
    use qcd_coupling,    only: running_coupling
    use dglap_holders,   only: dglap_holder
    implicit none
    type(c_ptr), intent(in), value :: tab
    real(c_double), intent(in), value :: Q0
    type(c_ptr), intent(in), value :: dh, coupling
    real(c_double), intent(in), value :: muR_over_Q
    integer(c_int), intent(in), value :: nloop
    logical(c_bool), intent(in), value :: untie_nf
    !-----------
    type(pdf_table), pointer        :: tab_f
    type(running_coupling), pointer :: coupling_f
    type(dglap_holder), pointer    :: dh_f
    logical :: untie_nf_f

    call c_f_pointer(tab, tab_f)
    call c_f_pointer(coupling, coupling_f)
    call c_f_pointer(dh, dh_f)
    untie_nf_f = (untie_nf .neqv. .false._c_bool) ! safely convert to fortran logical

    call PreEvolvePdfTable(tab_f, &
                        Q0, &
                        dh_f, coupling_f, muR_over_Q, nloop, untie_nf=untie_nf_f)
  end subroutine hoppet_cxx__pdf_table__pre_evolve


  !! carry out the evolution
  subroutine hoppet_cxx__pdf_table__evolve_frompre(tab,pdf_at_Q0) bind(C)
    use, intrinsic :: iso_c_binding
    use qcd_coupling,    only: running_coupling
    use dglap_holders,   only: dglap_holder
    implicit none
    type(c_ptr), intent(in), value :: tab
    type(c_ptr), intent(in), value :: pdf_at_Q0
    !-----------
    real(c_double), pointer :: pdf_at_Q0_f(:,:)
    type(pdf_table), pointer        :: tab_f

    call c_f_pointer(tab, tab_f)
    call c_f_pointer(pdf_at_Q0, pdf_at_Q0_f, shape=[size(tab_f%tab,1), size(tab_f%tab,2)])

    call EvolvePdfTable(tab_f, pdf_at_Q0_f)
  end subroutine 

  !! interface to write LHAPDF format files
  subroutine hoppet_cxx__pdf_table__write_lhapdf(table, coupling, basename, pdf_index, &
                                 & iy_increment, &
                                 & n_flav, flav_indices, flav_pdg_ids, flav_rescale) bind(C)                                
    use qcd_coupling,    only: running_coupling
    use hoppet_c_f_string_utils
    use warnings_and_errors
    implicit none
    type(c_ptr),    intent(in), value :: table
    type(c_ptr),    intent(in), value :: coupling
    type(c_ptr),    intent(in), value :: basename
    integer(c_int), intent(in), value :: pdf_index
    integer(c_int), intent(in), optional :: iy_increment
    integer(c_int), intent(in), optional :: n_flav !! the size of the flav_indices, etc. arrays
    type(c_ptr),    intent(in)           :: flav_indices !! integer array
    type(c_ptr),    intent(in)           :: flav_pdg_ids !! integer array
    type(c_ptr),    intent(in)           :: flav_rescale !! real(c_double) array
    !--
    type(pdf_table), pointer        :: table_f
    type(running_coupling), pointer :: coupling_f
    character(len=:), allocatable :: basename_f
    integer, pointer :: flav_indices_f(:)
    integer, pointer :: flav_pdg_ids_f(:)
    real(dp), pointer :: flav_rescale_f(:)

    call c_f_pointer(table, table_f)
    call c_f_pointer(coupling, coupling_f)
    basename_f = fortran_string_from_cstr(basename)
    if (present(n_flav)) then
      if (.not. c_associated(flav_indices) .or. &
          .not. c_associated(flav_pdg_ids) .or. &
          .not. c_associated(flav_rescale)) then
        call wae_error("hoppet_cxx__pdf_table__write_lhapdf: "//&
                      "if n_flav is present, flav_indices, flav_pdg_ids and flav_rescale must also be present")
      end if
      call c_f_pointer(flav_indices, flav_indices_f, shape=[n_flav])
      call c_f_pointer(flav_pdg_ids, flav_pdg_ids_f, shape=[n_flav])
      call c_f_pointer(flav_rescale, flav_rescale_f, shape=[n_flav])
      call WriteLHAPDFFromPdfTable(table_f, coupling_f, basename_f, pdf_index, &
                                iy_increment=iy_increment, &
                                flav_indices = flav_indices_f, &
                                flav_pdg_ids = flav_pdg_ids_f, flav_rescale = flav_rescale_f)
    else
      call WriteLHAPDFFromPdfTable(table_f, coupling_f, basename_f, pdf_index, &
                                iy_increment=iy_increment)
    end if
  end subroutine hoppet_cxx__pdf_table__write_lhapdf


  DEFINE_DELETE(pdf_table)

  ! think carefully which of these interfaces should be public, which perhaps renamed
  DEFINE_RETURN_OBJ_MEMBER(pdf_table,grid,grid_def)
  DEFINE_RETURN_INT_MEMBER(pdf_table,nQ)
  DEFINE_RETURN_INT_MEMBER(pdf_table,tab_iflv_max)
  DEFINE_RETURN_INT_MEMBER(pdf_table,lnlnQ_order)
  DEFINE_RETURN_LOG_MEMBER(pdf_table,freeze_at_Qmin)
  DEFINE_RETURN_LOG_MEMBER(pdf_table,nf_info_associated)
  DEFINE_RETURN_INT_MEMBER(pdf_table,nflo)
  DEFINE_RETURN_INT_MEMBER(pdf_table,nfhi)
  DEFINE_RETURN_DBL_MEMBER(pdf_table,lnlnQ_min)
  DEFINE_RETURN_DBL_MEMBER(pdf_table,lnlnQ_max)
  DEFINE_RETURN_DBL_MEMBER(pdf_table,lambda_eff)
  DEFINE_RETURN_DBL_MEMBER(pdf_table,dlnlnQ)
  DEFINE_RETURN_OBJ_MEMBER(pdf_table,seginfo_no_nf,pdfseginfo)


  DEFINE_RETURN_OBJ_MEMBER_I(pdf_table,seginfo,pdfseginfo)
  DEFINE_RETURN_DBL_MEMBER_I(pdf_table,as2pi)
  DEFINE_RETURN_DBL_MEMBER_I(pdf_table,lnlnQ_vals)
  DEFINE_RETURN_DBL_MEMBER_I(pdf_table,Q_vals)
  DEFINE_RETURN_INT_MEMBER_I(pdf_table,nf_int)
end module hoppet_cxx_oo_pdf_table

!=====================================================================
!! wrappers for C++ OO interface to grid quantities, with an intermediate
!! grid_quant type holds both the data and the grid_def
module hoppet_cxx_oo_pdf_table_array
  use types, only: dp
  use, intrinsic :: iso_c_binding
  use pdf_tabulate
  implicit none

  type pdf_table_array
    type(pdf_table), pointer :: tables(:) => null()
  end type pdf_table_array

contains

  !! create a pdf_table_array of the specified size and set up pointers
  !! to the pdf_table_array and to its tables(0) element
  subroutine hoppet_cxx__pdf_table_array__new(sz, ptr, table0_ptr) bind(C)
    implicit none
    integer(c_int), intent(in), value :: sz !< size of the array
    type(c_ptr), intent(out) :: ptr         !< pointer to the pdf_table_array object
    type(c_ptr), intent(out) :: table0_ptr  !< pointer to zeroth element of pdf_table_array
    type(c_ptr) :: res
    !--
    type(pdf_table_array), pointer :: ptr_f

    allocate(ptr_f)
    allocate(ptr_f%tables(0:sz-1))
    ptr = c_loc(ptr_f)
    table0_ptr = c_loc(ptr_f%tables(0))
  end subroutine hoppet_cxx__pdf_table_array__new

  !! delete a pdf_table_array created with hoppet_cxx__pdf_table_array__new
  subroutine hoppet_cxx__pdf_table_array__delete(ptr) bind(C)
    implicit none
    type(c_ptr), intent(out) :: ptr         !< pointer to the pdf_table_array object
    !--
    type(pdf_table_array), pointer :: ptr_f

    call c_f_pointer(ptr, ptr_f)
    call Delete(ptr_f%tables) ! clean up the storage in each pdf_table
    deallocate(ptr_f%tables)  ! deallocate the array of pdf_table pointers 
    deallocate(ptr_f)         ! deallocate the pdf_table_array itself
    ptr = c_null_ptr
  end subroutine hoppet_cxx__pdf_table_array__delete

  function hoppet_cxx__pdf_tables__table_i(ptr, sz, i) bind(C) result(res)
    implicit none
    type(c_ptr),    intent(in), value :: ptr !< pointer to zeroth element of pdf_table_array
    integer(c_int), intent(in), value :: sz  !< size of the array
    integer(c_int), intent(in), value :: i   !< index of the desired pdf_table (zero-indexed)
    type(c_ptr) :: res
    !--
    type(pdf_table), pointer :: ptr_f(:)

    call c_f_pointer(ptr, ptr_f, shape=[sz])
    ! i is zero indexed in C++, but Fortran arrays from c_f_pointer are one-indexed
    ! (and this can be changed only from Fortran 2023, which we are not yet using)
    res = c_loc(ptr_f(i+1))
  end function 

end module hoppet_cxx_oo_pdf_table_array
