
#include "inc/ftlMacros.inc"

! define a macro to generate the functions that returns a generic object's integer member
#define DEFINE_RETURN_INT_MEMBER(OBJ,NAME) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj) bind(C) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    integer(c_int) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = obj_f%NAME;\
  end function

  ! define a macro to generate the functions that returns a generic object's real(dp) member
#define DEFINE_RETURN_DOUBLE_MEMBER(OBJ,NAME) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj) bind(C) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    real(c_double) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = obj_f%NAME;\
  end function

#define DEFINE_RETURN_OBJ_MEMBER(OBJ,NAME,TYPE) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj) bind(C) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    type(c_ptr) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = c_loc(obj_f%NAME);\
  end function

#define DEFINE_RETURN_OBJ_MEMBER_I(OBJ,NAME,TYPE) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj,i) bind(C) result(res);\
    implicit none;\
    type(c_ptr), intent(in), value :: obj;\
    integer(c_int), intent(in), value :: i;\
    type(c_ptr) :: res;\
    type(OBJ), pointer :: obj_f;\
    call c_f_pointer(obj, obj_f);\
    res = c_loc(obj_f%NAME(i));\
  end function

#define DEFINE_RETURN_OBJ_MEMBER_IJ(OBJ,NAME,TYPE) \
  function CAT4(hoppet_cxx__,OBJ,__,NAME)(obj,i,j) bind(C) result(res);\
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
  function hoppet_cxx__grid_conv__new_from_fn(grid_ptr, conv_ignd_c_fn_obj) bind(C) result(ptr)
    implicit none
    type(c_ptr), intent(in), value :: grid_ptr
    type(c_ptr), intent(in), value :: conv_ignd_c_fn_obj
    !type(c_ptr), intent(in) :: split_array_ptr
    !integer(c_int), intent(in), value :: nsplit
    type(c_ptr) :: ptr
    !--
    type(grid_conv), pointer :: gc
    type(grid_def), pointer :: grid
    procedure(conv_ignd_c_interface), pointer :: conv_ignd_f_fn
    type(conv_ignd_from_c) :: lcl_conv_ignd

    call c_f_pointer(grid_ptr, grid)
    lcl_conv_ignd%ctx = conv_ignd_c_fn_obj

    allocate(gc)
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


  !! result_data = conv * q_data
  subroutine hoppet_cxx__split_mat__times_grid_quant_2d(split_max_ptr, q_data, result_data) bind(C)
    use pdf_representation
    implicit none    
    type(c_ptr), intent(in), value :: split_max_ptr
    type(c_ptr), intent(in), value :: q_data
    type(c_ptr), intent(in), value :: result_data
    !--
    type(split_mat), pointer :: split_mat_f
    real(dp), pointer :: q(:,:)
    real(dp), pointer :: result(:,:)

    call c_f_pointer(split_max_ptr, split_mat_f)
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
#define HOPPET_CXX__SPLIT_MAT__REF(NAME) \
  function CAT(hoppet_cxx__split_mat__,NAME)(split_mat_c_ptr) bind(C) result(grid_conv_c_ptr);\
    implicit none;\
    type(c_ptr), intent(in), value :: split_mat_c_ptr;\
    type(c_ptr) :: grid_conv_c_ptr;\
    type(split_mat), pointer :: f_ptr;\
    type(grid_conv), pointer :: gg_ptr;\
    call c_f_pointer(split_mat_c_ptr, f_ptr);\
    grid_conv_c_ptr = c_loc(CAT(f_ptr%,NAME));\
  end function CAT(hoppet_cxx__split_mat__,NAME)

  HOPPET_CXX__SPLIT_MAT__REF(qq      )
  HOPPET_CXX__SPLIT_MAT__REF(qg      )
  HOPPET_CXX__SPLIT_MAT__REF(gq      )
  HOPPET_CXX__SPLIT_MAT__REF(gg      )
  HOPPET_CXX__SPLIT_MAT__REF(ns_plus )
  HOPPET_CXX__SPLIT_MAT__REF(ns_minus)
  HOPPET_CXX__SPLIT_MAT__REF(ns_v    )

#undef HOPPET_CXX__SPLIT_MAT__REF

end module hoppet_cxx_oo_split_mat

!! module for the wrapping a running_coupling object
module hoppet_cxx_oo_running_coupling
  use, intrinsic :: iso_c_binding
  use qcd_coupling
  implicit none
contains

  function hoppet_cxx__running_coupling__new_fixnf(alphas_at_Q, Q, nloop, fixnf) bind(C) result(res)
    real(c_double), intent(in), value :: alphas_at_Q, Q
    integer(c_int), intent(in), value :: nloop, fixnf
    type(c_ptr) :: res
    !--
    type(running_coupling), pointer :: rc_f
    allocate(rc_f)
    call InitRunningCoupling(rc_f, alphas_at_Q, Q, nloop, fixnf=fixnf)
    res = c_loc(rc_f)
  end function hoppet_cxx__running_coupling__new_fixnf

  function hoppet_cxx__running_coupling__new_varnf(alphas_at_Q, Q, nloop, mc, mb, mt, masses_are_MSbar, muMatch_mQuark) bind(C) result(res)
    real(c_double), intent(in), value :: alphas_at_Q, Q
    integer(c_int), intent(in), value :: nloop
    real(c_double), intent(in), value :: mc, mb, mt
    logical(c_bool), intent(in), value :: masses_are_MSbar
    real(c_double), intent(in), value :: muMatch_mQuark
    type(c_ptr) :: res
    !--
    type(running_coupling), pointer :: rc_f
    logical :: masses_are_MSbar_f
    masses_are_MSbar_f = (masses_are_MSbar .neqv. .false._c_bool) ! safely convert to fortran logical
    allocate(rc_f)
    call InitRunningCoupling(rc_f, alphas_at_Q, Q, nloop, quark_masses = [mc, mb, mt], &
                             masses_are_MSbar = masses_are_MSbar_f, muMatch_mQuark = muMatch_mQuark)
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

  subroutine hoppet_cxx__dglap_holder__set_nf(dh, nf) bind(C)
    implicit none
    type(c_ptr), intent(in), value :: dh
    integer(c_int), intent(in), value :: nf
    !--
    type(dglap_holder), pointer :: dh_f

    call c_f_pointer(dh, dh_f)
    call SetNfDglapHolder(dh_f, nf)
  end subroutine hoppet_cxx__dglap_holder__set_nf


  DEFINE_DELETE(dglap_holder)
  DEFINE_RETURN_INT_MEMBER(dglap_holder,nloop)
  DEFINE_RETURN_INT_MEMBER(dglap_holder,nf)
  DEFINE_RETURN_INT_MEMBER(dglap_holder,factscheme)
  DEFINE_RETURN_OBJ_MEMBER_IJ(dglap_holder,allp,split_mat)


end module hoppet_cxx_oo_dglap_holder