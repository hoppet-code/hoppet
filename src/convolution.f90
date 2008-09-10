! $Id: convolution.f90,v 1.26 2005/07/08 21:19:17 salam Exp $


!======================================================================
! Module exlcusively for communication between convolution and 
! any routine which might supply convolutions
!======================================================================
module convolution_communicator
  use types
  implicit none
  private

  integer, public , parameter :: cc_REAL=1,cc_VIRT=2,&
       & cc_REALVIRT=3,cc_DELTA =4 
  integer, public, save :: cc_piece
end module convolution_communicator



!======================================================================
! All the base types and routines for happy convolution
!
! Decide on the following policy: grid details will be copied rather
! than pointed to (avoids danger of object that is pointed to being 
! changed). Not clear if this is really the best policy, but we will
! give it a try...
!
! GPS 04/04/2000
!======================================================================
module convolution
  use types; use consts_dp; use assertions; use warnings_and_errors
  implicit none
  private

  integer, parameter :: conv_UndefinedInt = 1000000004

  !-------------------------------
  ! definition of a grid
  ! Includes possibility of one layer of subsiduary definitions
  type grid_def
     real(dp) :: dy, ymax, eps
     integer  :: ny, order, nsub
     logical  :: locked
     integer,        pointer :: subiy(:) ! starting points of subsiduary grid
     type(grid_def), pointer :: subgd(:) ! subsiduary grid defs
  end type grid_def

  !--------------------------------------------------
  ! used for abbreviated access to EvalGridQuant
  
  
  type grid_conv
     type(grid_def)           :: grid
     real(dp),        pointer :: conv(:,:) => null()
     !-- support in construction...
     type(grid_conv), pointer :: subgc(:) => null()
  end type grid_conv
  integer, parameter :: FULL=1, UPPR=2
  !-- for standard linear approach with proper end-point treatment.
  ! This must remain zero otherwise inconsistencies will arise
  integer, parameter :: LIN_ORDER=0
  

  public :: grid_def, grid_conv


  !-- public routines --
  interface InitGridDef
     module procedure conv_InitGridDef_single, conv_InitGridDef_multi
  end interface
  public :: conv_InitGridDef_single, conv_InitGridDef_multi
  public :: InitGridDef, ValidateGD, GetGridInfoString
  public :: SetDefaultConvolutionEps, DefaultConvolutionEps
  interface operator(==)
     module procedure conv_CmpGridDef
  end interface
  public :: operator(==)
  interface Delete ! made public later...
     module procedure Delete_grid_def_0d, Delete_grid_def_1d
  end interface
  public :: xValues, yValues

  !-- quant routines -----------------------------------------------
  interface AllocGridQuant
     module procedure conv_AllocGridQuant_1d,conv_AllocGridQuant_2d,&
          & conv_AllocGridQuant_3d
  end interface
  public :: AllocGridQuant
  interface InitGridQuant
     module procedure conv_InitGridQuant_func, conv_InitGridQuant_func2d,&
          & conv_InitGridQuant_func_a, conv_InitGridQuant_func2d_a,&
          & conv_InitGridQuant_func_ai, conv_InitGridQuant_func2d_ai
  end interface
  ! standard version gives problems with pgf90 5.1-2 (wrong answers)...
  interface InitGridQuantSub
     module procedure conv_InitGridQuantSub_2d, &
          &conv_InitGridQuantSub_2d_a, conv_InitGridQuantSub_2d_ai
  end interface
  ! This is the workaround suggested by Giulia, across all variants
  ! of the subroutine... Testing with pgf90 however shows that we still 
  ! have the wrong answer (though whether because of this routine
  ! or others has not been established).
  !interface InitGridQuantSub
  !   module procedure conv_InitGridQuantSub_2d_pgf, &
  !        &conv_InitGridQuantSub_2d_a_pgf, conv_InitGridQuantSub_2d_ai_pgf
  !end interface
  
  public :: InitGridQuant, InitGridQuantSub, PrintGridQuant
  public :: InitGridQuantLHAPDF
  interface PrintGridQuant
     module procedure conv_PrintGridQuant_1, conv_PrintGridQuant_2,&
          & conv_PrintGridQuant_3, conv_PrintGridQuant_4
  end interface
  interface EvalGridQuant
     module procedure conv_EvalGridQuant_0d, conv_EvalGridQuant_1d
  end interface
  public :: MomGridQuant, EvalGridQuant, WgtGridQuant
  interface Delete
     module procedure conv_DelGridQuant_1d, conv_DelGridQuant_2d
     module procedure conv_DelGridQuant_3d
  end interface
  
  !public :: Delete

  !-- for calculating moments ---
  real(dp)       :: conv_moment_index
  type(grid_def) :: conv_moment_gd
  real(dp), allocatable :: conv_moment_gq(:)
  interface TruncatedMoment
    module procedure conv_TruncatedMoment_1d
  end interface
  public :: TruncatedMoment

  !-- convolution routines -------------------------------------------
  interface AllocGridConv
     module procedure conv_AllocGridConv_0d, conv_AllocGridConv_1d, &
          & conv_AllocGridConv_2d
  end interface
  interface InitGridConv
     module procedure conv_InitGridConv_zero, conv_InitGridConv_zero_1d,&
          & conv_InitGridConv_zero_2d, conv_InitGridConv_func,&
          & conv_InitGridConv_gc, conv_InitGridConv_gc_1d, &
          & conv_InitGridConv_gc_2d, conv_InitGridConv_conv
  end interface
  !-- keep following for historical reasons? Not as of 28/12/01
  !interface conv_GridConvAdd
  !   module procedure conv_GridConvAdd_func, conv_GridConvAdd_gc
  !end interface
  interface SetToZero
     module procedure conv_ZeroGridConv_0d, conv_ZeroGridConv_1d,&
          & conv_ZeroGridConv_2d
  end interface
  interface AddWithCoeff
     module procedure conv_AddGridConv_func, conv_AddGridConv_gc,&
          & conv_AddGridConv_gc_1d, conv_AddGridConv_gc_2d
  end interface
  interface Multiply
     module procedure conv_MultGridConv_0d, conv_MultGridConv_1d, &
          & conv_MultGridConv_2d
  end interface
  interface Delete
     module procedure conv_DelGridConv_0d, conv_DelGridConv_1d, &
          & conv_DelGridConv_2d
  end interface
  interface SetToConvolution
     module procedure conv_ConvGridConv_0d, conv_ConvGridConv_2dx2d
  end interface
  interface SetToCommutator
     module procedure SetToCommutator_gc
  end interface
  

  public :: GridConvAllocated
  public :: AllocGridConv, InitGridConv, Delete
  public :: AddWithCoeff, Multiply, SetToConvolution
  public :: SetToZero, SetToCommutator

  !-- REMEMBER: .conv. seems to have a precedence similar to .or., .and.
  !             and so is lower precednece than any arithmetic operation
  interface operator(.conv.)
     module procedure conv_ConvGridQuant_mat, conv_ConvGridQuant_scalar
  end interface
  public :: operator(.conv.)
  !-- put in this version as well because to ease the writing of 
  !   expressions: * should have the right precedence properties
  !   where .conv. seems to have a very low precedence
  interface operator(*)
     module procedure conv_ConvGridQuant_mat, conv_ConvGridQuant_scalar
  end interface
  public :: operator(*)

  !-- keep this for historical reasons (hopefully not too long)
  interface conv_ConvConv
     module procedure conv_InitGridConv_conv
  end interface
  public :: conv_ConvConv

  !--------------------------
  type gdval
     type(grid_def) :: grid
     real(dp)       :: val
  end type gdval
  public :: gdval
  interface operator(.with.)
     module procedure conv_gdval_gdv, conv_gdval_vgd
  end interface
  public :: operator(.with.)  
  interface operator(.atx.)
     module procedure conv_EvalGridQuant_atx, conv_EvalGridQuant_atx_1d
  end interface
  public :: operator(.atx.)
  interface operator(.aty.)
     module procedure conv_EvalGridQuant_aty, conv_EvalGridQuant_aty_1d
  end interface
  public :: operator(.aty.)

  !-- precision used for integration
  real(dp) :: default_conv_eps=1e-7_dp
  !real(dp), parameter :: default_conv_eps=1e-8_dp
  real(dp), parameter :: warn_tolerance = 1e-3_dp
  
  !-- used for generation of derived convoluters (e.g. exponentiations
  !   of existsing convoluters)
  logical :: override_grid_locking = .false.
  integer :: nconv_with_override_off = 0 ! 
  public :: GetDerivedProbes, SetDerivedConv
  public :: SetDerivedConv_nodealloc

  ! actual subroutine is to be found in welcome_message.f90
  interface
     subroutine HoppetWelcomeMessage
     end subroutine HoppetWelcomeMessage
  end interface
  public :: HoppetWelcomeMessage

contains

  !======================================================================
  ! Things just for Grid defs.
  ! updated for multi.
  !
  ! order = 0 -> standard linear
  ! order > 0 -> full (order) treatment
  ! order < 0 -> treatment correct to |order| except at end points 
  !              (this is similar to Ratcliffes proposal, hep-ph/0012376)
  subroutine conv_InitGridDef_single(grid,dy,ymax,order,eps)
    type(grid_def), intent(out) :: grid
    real(dp),       intent(in)  :: dy, ymax
    integer,        intent(in), optional :: order
    real(dp),       intent(in), optional :: eps
    logical, save :: first_call = .true.
    character(len=80) :: string1, string2
    integer, save     :: warn_dy = 4

    if (first_call) then
       first_call = .false.
       call HoppetWelcomeMessage
    end if
    
    !-- this is a plain grid def ---------------------------
    grid%nsub = 0
    nullify(grid%subiy,grid%subgd)
    grid%locked = .false.
    !-------------------------------------------------------

    grid%ymax = ymax
    grid%ny   = nint(ymax / dy)
    if (grid%ny < 2) then
       
       write(string1,*) 'InitGridDef: requested too small a number of bins'
       write(string2,*) '                 dy and ymax were',dy,ymax
       call wae_error(trim(string1),trim(string2))
    end if
    
    grid%dy   = ymax / grid%ny
    if (abs(grid%dy/dy - one) > 0.001_dp) then
       write(string1,*) 'InitGridDef: requested dy of ', dy
       write(string2,*) '                         provided  dy of ', grid%dy
       call wae_warn(warn_dy,trim(string1),trim(string2))
    end if

    if (present(order)) then
       !if (order < LIN_ORDER) then
       !   write(0,*) 'ERROR in InitGridDef: order < 0&
       !        & is not valid in InitGridDef'
       !   stop
       if (abs(order)+1 > grid%ny) then
          call wae_error('InitGridDef: |order|+1 > number of grid points')
       end if
       grid%order = order
    else
       grid%order = LIN_ORDER
    end if
    
    grid%eps = default_or_opt(default_conv_eps, eps)
  end subroutine conv_InitGridDef_single


  !--------------------------------------------------------------
  ! Create a multi grid def
  subroutine conv_InitGridDef_multi(grid,gdarray,locked)
    use sort
    type(grid_def), intent(out) :: grid
    type(grid_def), intent(in)  :: gdarray(:)
    logical,        intent(in), optional :: locked
    !----------------------------------
    integer :: i,indx(size(gdarray))
    real(dp) :: dyratio, approx_ny, new_ymax
    type(grid_def), pointer :: subgd(:) ! shorthand
    ! temp needed for workaround on ifort 8.0.039
    real(dp) :: gdarraydy(size(gdarray))
    character(len=80) :: string1, string2
    integer, save     :: warn_locking = 4

    !-- enforce one layer only
    if (any(gdarray(:)%nsub /= 0)) then
       call wae_error(&
        'ERROR in conv_InitGridDef_multi:',&
        'One of grid defs in array was a compound grid def.',&
        'Only one layer of compounding is currently allowed.')
    end if
    
    grid%nsub = size(gdarray)
    allocate(grid%subiy(grid%nsub+1))
    allocate(grid%subgd(grid%nsub));   subgd => grid%subgd

    grid%locked = default_or_opt(.false.,locked)
    if (grid%locked) then
       ! this calls gives wrong results with ifort-8.0.039
       ! (dummy array gdarray%dy is corrupted in indexx). Issue
       ! submitted 27/01/2004: 225712 
       ! "Corrupted answer on call with array of derived type subcomponents"
       !call indexx(gdarray%dy, indx)

       ! workaround for ifort-8.0.039 
       gdarraydy = gdarray%dy
       call indexx(gdarraydy, indx)
       do i = 1, grid%nsub
          subgd(i) = gdarray(indx(i))
          if (i > 1) then 
             if (subgd(i)%ymax <  subgd(i-1)%ymax) then
                write(string1,*) &
                     &'ERROR in conv_InitGridDef_multi: for locking,'
                write(string2,*) 'gdarray with smaller dy should&
                     & also have smaller gdarray%ymax'
                call wae_error(trim(string1), trim(string2))
                ! for testing ifort_8_0_039...
                !write(0,*) 'dy   values (i-1,i)', subgd(i-1:i)%dy
                !write(0,*) 'ymax values (i-1,i)', subgd(i-1:i)%ymax
                !write(0,*) indx
                !write(0,*) gdarray%dy
                !write(0,*) gdarray(indx(:))%dy
                !write(0,*) i, subgd(i)%ymax, subgd(i-1)%ymax
                !stop
             end if
          end if
       end do
       
       !-- now ensure that there is proper matching between locked grids
       do i = grid%nsub-1, 1, -1
          ! dyratio must be an integer
          dyratio = subgd(i+1)%dy / subgd(i)%dy
          subgd(i)%dy = subgd(i+1)%dy / nint(dyratio)
          if (abs(dyratio-nint(dyratio)) > warn_tolerance*dyratio) then
             write(string1,'(a,i2,a,f18.14)') ' InitGridDef (locking):&
                  & redefined dy(', i,') to be ', subgd(i)%dy
             call wae_warn(warn_locking, trim(string1))
          end if
          ! after fixing dy one must still have an integer number of bins
          approx_ny = subgd(i)%ymax/subgd(i)%dy
          subgd(i)%ny = ceiling(approx_ny - warn_tolerance)
          new_ymax = subgd(i)%ny * subgd(i)%dy
          if (abs(new_ymax-subgd(i)%ymax) > warn_tolerance*new_ymax) then
             write(string1,'(a,i2,a,f18.14)') ' InitGridDef (locking):&
                  & redefined ymax(', i,') to be ', new_ymax
             call wae_warn(warn_locking, trim(string1))
          end if
          subgd(i)%ymax = new_ymax
          subgd(i)%ny   = nint(subgd(i)%ymax / subgd(i)%dy)
          ! condition on order must still hold
          if (abs(subgd(i)%order)+1 > subgd(i)%ny) then
             write(string1,'(a)') 'Error in InitGridDef (locking):'
             write(string2,'(a,i2,a)') '       For grid def ',i,' |order|+1 > ny'
             call wae_error(trim(string1),trim(string2))
          end if
       end do
    else
       !-- no questions asked!
       subgd(:) = gdarray(:)
    end if
    
    grid%dy       = zero
    grid%ymax     = maxval(subgd(:)%ymax)
    grid%order    = conv_UndefinedInt
    !-- arrays will stretch from 0:grid%ny; so must get this right!
    grid%ny       = sum(subgd(:)%ny) + grid%nsub-1
    
    !-- indicate starting points of arrays
    grid%subiy(1) = 0
    do i = 2, grid%nsub+1
       grid%subiy(i) = grid%subiy(i-1) + subgd(i-1)%ny + 1
    end do
    
  end subroutine conv_InitGridDef_multi
  

  !======================================================================
  !! Set a string with a compact summary about the grid
  !!
  subroutine GetGridInfoString(grid,string)
    use sort
    type(grid_def),   intent(in)  :: grid
    character(len=*), intent(out) :: string
    !----------------------------------------
    character(len=200) :: my_string ! nice and long!
    character(len=200) :: frmt
    integer, save :: nwarn_len = 4, nwarn_eps = 4
    integer, allocatable  :: indx(:)

    if (grid%nsub == 0) then
       write(my_string,&
            &'("dy=",f5.3,"|ymax=",f0.2,"|order=",i0,"|eps=",es7.1)') &
            & grid%dy,grid%ymax,grid%order,grid%eps
    else 
       allocate(indx(grid%nsub))
       call indexx(-grid%subgd%dy,indx)  ! grid%subgd(indx)%dy == decreasing dy
       if (any(grid%subgd(2:)%eps /= grid%subgd(1)%eps)) &
            &call wae_warn(nwarn_eps,'GetGridInfoString: &
            &subgrids have different eps; only largest printed')
       if (grid%locked) then
          write(frmt,'(a,i0,a,i0,a,i0,a)') &
               &'("dy=",f5.3,',grid%nsub-1 ,&
               &'("/",i0),"|ymax=",f0.2,',grid%nsub-1,&
               &'(":",f0.2),"|order=",sp,i0,',grid%nsub-1,&
               &'(":",i0),ss,"|eps=",es7.1)'
          write(my_string,trim(frmt)) grid%subgd(indx(1))%dy,&
               & nint(grid%subgd(indx(1:grid%nsub-1))%dy/grid%subgd(indx(2:))%dy),&
               & grid%subgd(indx)%ymax, grid%subgd(indx)%order,&
               & maxval(grid%subgd(:)%eps)
       else
          write(frmt,'(a,i0,a,i0,a,i0,a)') &
               &'("dy=",f5.3,',grid%nsub-1 ,&
               &'(":",f5.4),"|ymax=",f0.2,',grid%nsub-1,&
               &'(":",f0.2),"|order=",sp,i0,',grid%nsub-1,&
               &'ss,(":",i0),"|eps=",es7.1)'
          write(my_string,trim(frmt)) grid%subgd(indx(:))%dy,&
               & grid%subgd(indx)%ymax, grid%subgd(indx)%order,&
               & maxval(grid%subgd(:)%eps)
       end if
    end if
    
    if (len_trim(my_string) > len(string)) then
       call wae_warn(nwarn_len,&
            &"GetGridInfoString: too short a string was passed")
       string = my_string(1:len(string))
    else
       string = my_string(1:len_trim(my_string))
    end if
    
  end subroutine GetGridInfoString
  

  !======================================================================
  !! Delete the memory associated with a grid definition
  recursive subroutine Delete_grid_def_0d(grid)
    type(grid_def), intent(inout) :: grid

    if (grid%nsub /= 0) then
       call Delete(grid%subgd(:)) 
       !do isub = 1, grid%nsub
       !   ! following line not strictly needed since (a) we have no more
       !   ! than one level of depth and (b) we do not have anything to
       !   ! delete in the lowest level.
       !   call Delete(grid%subgd(isub)) 
       !end do
       deallocate(grid%subiy)
       deallocate(grid%subgd)
    end if
    
  end subroutine Delete_grid_def_0d
  

  !======================================================================
  recursive subroutine Delete_grid_def_1d(grid_array)
    type(grid_def), intent(inout) :: grid_array(:)
    integer :: igrid
    do igrid = 1, size(grid_array)
       call Delete(grid_array(igrid))
    end do
  end subroutine Delete_grid_def_1d
  
  !======================================================================
  !! Returns .true. iff the two grid definitions are equivalent (i.e. they
  !! may be different objects, but correspond to grids with identical 
  !! parameters). 
  !!
  !! Note that we do not check equality of the value of eps
  !! (integration precision) -- whether this is the right choice is
  !! arguable, but since it in no way affects interoperability of two
  !! grids (e.g. for convolutions between grid_conv objects), it's not
  !! too unreasonable.
  recursive function conv_CmpGridDef(gd1,gd2) result(equal)
    type(grid_def), intent(in) :: gd1, gd2
    logical :: equal
    integer :: i
    logical, parameter :: verbose = .false.

    if (gd1%nsub /= gd2%nsub) then
       equal = .false.
       return
    else if (gd1%nsub == 0) then
       if (verbose) write(0,*) gd1%dy,gd2%dy, gd1%ny, gd2%ny, gd1%ymax,gd2%ymax
       equal = (gd1%dy   == gd2%dy) .and. (gd1%ny == gd2%ny) .and.&
            &  (gd1%ymax == gd2%ymax) .and. (gd1%order == gd2%order)
    else
       equal = .true.
       if (.not.associated(gd1%subgd,gd2%subgd)) then
          do i = 1, gd1%nsub
             equal = equal .and. conv_CmpGridDef(gd1%subgd(i),gd2%subgd(i))
          end do
       end if
       equal = equal .and. (gd1%locked .eqv. gd2%locked)
    end if
  end function conv_CmpGridDef

  
  !-- useful to be able to check easily -------------------------
  ! no problem with multi
  subroutine ValidateGD(gd1,gd2,source)
    type(grid_def), intent(in) :: gd1, gd2
    character(len=*), intent(in) :: source
    if (.not. (gd1 == gd2)) then
       write(0,*) 'Problem validating two grid defs in ',source
       stop
    end if
  end subroutine ValidateGD

  !======================================================================
  !! Return an array containing the y-value of each point on the grid
  recursive function yValues(grid) result(res)
    type(grid_def), intent(in) :: grid
    real(dp)                   :: res(0:grid%ny)
    integer :: iy, isub

    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          res(grid%subiy(isub):grid%subiy(isub+1)-1) = yValues(grid%subgd(isub))
       end do
    else
       forall(iy=0:grid%ny) res(iy) = iy*grid%dy
    end if
  end function yValues

  
  !======================================================================
  !! Return an array containing the x-value of each point on the grid
  function xValues(grid) result(res)
    type(grid_def), intent(in) :: grid
    real(dp)                   :: res(0:grid%ny)

    res = exp(-yValues(grid))
  end function xValues


  !======================================================================
  ! Things just for grid quants
  ! multi makes no difference
  subroutine conv_AllocGridQuant_1d(grid,gq)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: gq(:)
    !integer :: istat
    ! this form of deallocate ought to be OK (?), but on lahey,
    ! compaq+condor (and perhaps elsewhere) it causes problems
    !deallocate(gq,stat=istat)
    allocate(gq(0:grid%ny))
  end subroutine conv_AllocGridQuant_1d

  !----------------------------------------------------------------------
  ! multi makes no difference
  subroutine conv_AllocGridQuant_2d(grid,gq,nl,nh)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: gq(:,:)
    integer,        intent(in) :: nl,nh
    !integer :: istat
    ! this form of deallocate ought to be OK (?), but on lahey,
    ! compaq+condor (and perhaps elsewhere) it causes problems
    !deallocate(gq,stat=istat)
    allocate(gq(0:grid%ny,nl:nh))
  end subroutine conv_AllocGridQuant_2d

  !----------------------------------------------------------------------
  ! multi makes no difference
  subroutine conv_AllocGridQuant_3d(grid,gq, nl2,nh2, nl3,nh3)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: gq(:,:,:)
    integer,        intent(in) :: nl2,nh2,nl3,nh3
    !integer :: istat
    ! this form of deallocate ought to be OK (?), but on lahey,
    ! compaq+condor (and perhaps elsewhere) it causes problems
    !deallocate(gq,stat=istat)
    allocate(gq(0:grid%ny,nl2:nh2,nl3:nh3))
  end subroutine conv_AllocGridQuant_3d
  
  !----------------------------------------------------------
  ! multi makes no difference
  subroutine conv_DelGridQuant_1d(gq)
    real(dp),       pointer    :: gq(:)
    integer :: istat
    deallocate(gq,stat=istat)
  end subroutine conv_DelGridQuant_1d
  !----------------------------------------------------------
  ! multi makes no difference
  subroutine conv_DelGridQuant_2d(gq)
    real(dp),       pointer    :: gq(:,:)
    integer :: istat
    deallocate(gq,stat=istat)
  end subroutine conv_DelGridQuant_2d
  !----------------------------------------------------------
  ! multi makes no difference
  subroutine conv_DelGridQuant_3d(gq)
    real(dp),       pointer    :: gq(:,:,:)
    integer :: istat
    deallocate(gq,stat=istat)
  end subroutine conv_DelGridQuant_3d
  

  
  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func(grid, gq, func)
    real(dp),         intent(inout) :: gq(0:)
    type(grid_def),   intent(in)    :: grid
    interface
       function func(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: func
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuant_func(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1), func)
       end do
    else
       do iy = 0, ny
          gq(iy) = func(iy*grid%dy)
       end do
    end if
    
  end subroutine conv_InitGridQuant_func


  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func_a(grid, gq, func, axtra)
    real(dp),         intent(inout) :: gq(0:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    interface
       function func(x,axtra)
         use types; implicit none
         real(dp), intent(in) :: x, axtra
         real(dp)             :: func
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuant_func_a(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1), func, axtra)
       end do
    else
       do iy = 0, ny
          gq(iy) = func(iy*grid%dy, axtra)
       end do
    end if
    
  end subroutine conv_InitGridQuant_func_a


  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func_ai(grid, gq, func, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       function func(x,axtra,ixtra)
         use types; implicit none
         real(dp), intent(in) :: x, axtra
         integer,  intent(in) :: ixtra
         real(dp)             :: func
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuant_func_ai(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1), func, axtra, ixtra)
       end do
    else
       do iy = 0, ny
          gq(iy) = func(iy*grid%dy, axtra, ixtra)
       end do
    end if
    
  end subroutine conv_InitGridQuant_func_ai


  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func2d(grid, gq, func)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    interface
       function func(x,n)
         use types; implicit none
         real(dp), intent(in) :: x
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny, n

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuant_func2d(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1,:), func)
       end do
    else
       n  = ubound(gq, dim=2)
       do iy = 0, ny
          gq(iy,:) = func(iy*grid%dy, n)
       end do
    end if
  end subroutine conv_InitGridQuant_func2d

  !----------------------------------------------------------------------
  ! version intended for use when there is an extra argument whose
  ! value is fixed and needs to be passed to func
  !
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func2d_a(grid, gq, func, axtra)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    interface
       function func(x,axtra,n)
         use types; implicit none
         real(dp), intent(in) :: x,axtra
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny, n

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func2d_a")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuant_func2d_a(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1,:), func, axtra)
       end do
    else
       n  = ubound(gq, dim=2)
       do iy = 0, ny
          gq(iy,:) = func(iy*grid%dy, axtra, n)
       end do
    end if
  end subroutine conv_InitGridQuant_func2d_a

  !----------------------------------------------------------------------
  ! version intended for use when there is an extra argument whose
  ! value is fixed and needs to be passed to func
  !
  ! updated for multi
  recursive subroutine conv_InitGridQuant_func2d_ai(grid, gq, func, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       function func(x,axtra,ixtra,n)
         use types; implicit none
         real(dp), intent(in) :: x,axtra
         integer , intent(in) :: ixtra,n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny, n

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func2d_a")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuant_func2d_ai(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1,:), func, axtra, ixtra)
       end do
    else
       n  = ubound(gq, dim=2)
       do iy = 0, ny
          gq(iy,:) = func(iy*grid%dy, axtra, ixtra, n)
       end do
    end if
  end subroutine conv_InitGridQuant_func2d_ai


  !----------------------------------------------------------------------
  !! Version added specially for initialising a PDF from a LHAPDF style 
  !! interface.
  recursive subroutine InitGridQuantLHAPDF(grid, gq, LHAsub, Q)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: Q
    interface
       subroutine LHAsub(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine LHAsub
    end interface
    !-----------------------------------------
    integer  :: iy, isub, ny
    real(dp) :: f(size(gq,dim=2))

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call InitGridQuantLHAPDF(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1,:), LHAsub, Q)
       end do
    else
       do iy = 0, ny
          call LHAsub(exp(-iy*grid%dy), Q, f)
          gq(iy,:) = f
       end do
    end if
  end subroutine InitGridQuantLHAPDF

  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuantSub_2d(grid, gq, sub)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    interface
       subroutine sub(y,res)
         use types; implicit none
         real(dp), intent(in)  :: y
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuantSub_2d(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1,:), sub)
       end do
    else
       do iy = 0, ny
          call sub(iy*grid%dy, gq(iy,:))
       end do
    end if
  end subroutine conv_InitGridQuantSub_2d

  recursive subroutine conv_InitGridQuantSub_2d_a(grid, gq, sub, axtra)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    interface
       subroutine sub(y, axtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuantSub_2d_a(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1,:), sub, axtra)
       end do
    else
       do iy = 0, ny
          call sub(iy*grid%dy, axtra, gq(iy,:))
       end do
    end if
  end subroutine conv_InitGridQuantSub_2d_a

  recursive subroutine conv_InitGridQuantSub_2d_ai(grid, gq, sub, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       subroutine sub(y, axtra, ixtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         integer,  intent(in)  :: ixtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuantSub_2d_ai(grid%subgd(isub), &
               &gq(grid%subiy(isub):grid%subiy(isub+1)-1,:), sub, axtra, ixtra)
       end do
    else
       do iy = 0, ny
          call sub(iy*grid%dy, axtra, ixtra, gq(iy,:))
       end do
    end if
  end subroutine conv_InitGridQuantSub_2d_ai


  !======================================================================
  ! Variants for buggy pgf90-5.1.2, using workaround suggested by
  ! Giulia
  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine conv_InitGridQuantSub_2d_pgf(grid, gq, sub)
    real(dp),         intent(out)   :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    interface
       subroutine sub(y,res)
         use types; implicit none
         real(dp), intent(in)  :: y
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny
    ! -- GZ:make a copy shifting indices to make it compatible with pgf90 
    real(dp) :: gq_cpy(1:size(gq,dim=1),1:size(gq,dim=2))

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuantSub_2d(grid%subgd(isub), &
               &gq_cpy(grid%subiy(isub)+1:grid%subiy(isub+1),:), sub)
       end do
    else
       do iy = 0, ny
          call sub(iy*grid%dy, gq_cpy(iy+1,:))
       end do
    end if
    do iy=0,ny
       gq(iy,:) = gq_cpy(iy+1,:)
    end do
  end subroutine conv_InitGridQuantSub_2d_pgf

  recursive subroutine conv_InitGridQuantSub_2d_a_pgf(grid, gq, sub, axtra)
    real(dp),         intent(out)   :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    interface
       subroutine sub(y, axtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny
    real(dp) :: gq_cpy(1:size(gq,dim=1),1:size(gq,dim=2))

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuantSub_2d_a(grid%subgd(isub), &
               &gq_cpy(grid%subiy(isub)+1:grid%subiy(isub+1),:), sub, axtra)
       end do
    else
       do iy = 0, ny
          call sub(iy*grid%dy, axtra, gq_cpy(iy+1,:))
       end do
    end if
    do iy=0,ny
       gq(iy,:) = gq_cpy(iy+1,:)
    end do
  end subroutine conv_InitGridQuantSub_2d_a_pgf

  recursive subroutine conv_InitGridQuantSub_2d_ai_pgf(grid, gq, sub, axtra, ixtra)
    real(dp),         intent(out)   :: gq(0:,:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       subroutine sub(y, axtra, ixtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         integer,  intent(in)  :: ixtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    integer :: iy, isub, ny
    real(dp) :: gq_cpy(1:size(gq,dim=1),1:size(gq,dim=2))

    ny = assert_eq(grid%ny,ubound(gq,dim=1),"conv_InitGridQuant_func")
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_InitGridQuantSub_2d_ai(grid%subgd(isub), &
               &gq_cpy(grid%subiy(isub)+1:grid%subiy(isub+1),:), sub, axtra, ixtra)
       end do
    else
       do iy = 0, ny
          call sub(iy*grid%dy, axtra, ixtra, gq_cpy(iy+1,:))
       end do
    end if
    do iy=0,ny
       gq(iy,:) = gq_cpy(iy+1,:)
    end do
  end subroutine conv_InitGridQuantSub_2d_ai_pgf


  !--------------------------------------------------------------------
  ! multi-grid done
  recursive function conv_EvalGridQuant_0d(grid, gq, y) result(f)
    use interpolation
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq(0:)
    real(dp), intent(in) :: y
    real(dp) :: f
    !-----------------------------------------
    integer, parameter :: npnt_min = 4, npnt_max = 10
    integer :: i, ny, npnt, isub
    real(dp), parameter :: resc_yvals(npnt_max) = (/ (i,i=0,npnt_max-1) /)
    real(dp) :: wgts(npnt_max)

    !write(0,*) y,grid%ny, ubound(gq,dim=1)
    ny = assert_eq(grid%ny,ubound(gq,dim=1),"EvalGridQuant")
    if (y > grid%ymax*(one+warn_tolerance)) then
       write(0,*) 'EvalGridQuant: &
            &requested function value beyond maximum'
       write(0,*) 'y = ', y, 'ymax=',grid%ymax
       stop
    end if
    if (grid%nsub /= 0) then
       isub = conv_BestIsub(grid,y)
       f = EvalGridQuant(grid%subgd(isub), &
            & gq(grid%subiy(isub):grid%subiy(isub+1)-1), y)
    else
       npnt = min(npnt_max, max(npnt_min, abs(grid%order)))
       
       i = min(max(floor(y / grid%dy)-(npnt-1)/2,0),ny-npnt+1)
       call uniform_interpolation_weights(y/grid%dy - i, wgts(1:npnt))
       f = sum(wgts(1:npnt)*gq(i:i+npnt-1))
       !-- this was less efficient...
       !call polint(resc_yvals(1:npnt),gq(i:i+npnt-1),y/grid%dy-i,f,df)
!!$    i = min(grid%ny - 1, floor(y / grid%dy))
!!$    ey = y/grid%dy - i
!!$    f  = (one-ey)*gq(i) + ey*gq(i+1)
    end if
    
  end function conv_EvalGridQuant_0d

  !--------------------------------------------------------------------
  !! 1-D version of grid evaluation.
  !! 
  !! NB: we rewrite everything above so as to avoid unnecessary looping.
  !!     It would be better to have a common routine that gives the
  !!     required set of points and weights?
  recursive function conv_EvalGridQuant_1d(grid, gq, y) result(f)
    use interpolation
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq(0:,1:)
    real(dp), intent(in) :: y
    real(dp) :: f(size(gq,dim=2))
    !-----------------------------------------
    integer, parameter :: npnt_min = 4, npnt_max = 10
    integer :: i, j, ny, npnt, isub
    !real(dp) :: ey, df
    real(dp), parameter :: resc_yvals(npnt_max) = (/ (i,i=0,npnt_max-1) /)
    real(dp) :: wgts(npnt_max)

    !write(0,*) y,grid%ny, ubound(gq,dim=1)
    ny = assert_eq(grid%ny,ubound(gq,dim=1),"EvalGridQuant")
    if (y > grid%ymax*(one+warn_tolerance)) then
       write(0,*) 'EvalGridQuant: &
            &requested function value beyond maximum'
       write(0,*) 'y = ', y, 'ymax=',grid%ymax
       stop
    end if
    if (grid%nsub /= 0) then
       isub = conv_BestIsub(grid,y)
       f = conv_EvalGridQuant_1d(grid%subgd(isub), &
            & gq(grid%subiy(isub):grid%subiy(isub+1)-1,:), y)
    else
       npnt = min(npnt_max, max(npnt_min, abs(grid%order)))
       
       i = min(max(floor(y / grid%dy)-(npnt-1)/2,0),ny-npnt+1)
       call uniform_interpolation_weights(y/grid%dy - i, wgts(1:npnt))
       do j = 1, size(f)
          f(j) = sum(wgts(1:npnt)*gq(i:i+npnt-1,j))
       end do
       
       !-- this was less efficient...
       !call polint(resc_yvals(1:npnt),gq(i:i+npnt-1),y/grid%dy-i,f,df)
!!$    i = min(grid%ny - 1, floor(y / grid%dy))
!!$    ey = y/grid%dy - i
!!$    f  = (one-ey)*gq(i) + ey*gq(i+1)
    end if
  end function conv_EvalGridQuant_1d

!O  !--------------------------------------------------------------------
!O  !! 1-D version of grid evaluation.
!O  !! 
!O  !! NB: this is currently VERY slow because interpolation is recalculated
!O  !!     from scratch for each component of the grid quantity
!O  function conv_EvalfGridQuant_1d(grid, gq, y) result(f)
!O    type(grid_def), intent(in) :: grid
!O    real(dp), intent(in) :: gq(0:,1:)
!O    real(dp), intent(in) :: y
!O    real(dp) :: f(size(gq,dim=2))
!O    integer :: i
!O    do i = 1, size(gq,dim=2)
!O       f(i) = conv_EvalGridQuant_0d(grid,gq(:,i),y)
!O    end do
!O  end function conv_EvalfGridQuant_1d
  

  !--------------------------------------------------------------------
  ! Returns starting iymin point and a set of weights in order to calculate
  ! the value of the function at y -- one day we might introduce some
  ! option of setting the number of points; but not for now...
  !
  ! Qu: is the relation between number of points and order correct? It seems
  !     like we ought to have abs(grid%order)+1...
  recursive subroutine WgtGridQuant(grid, y, iylo, wgts)
    use interpolation
    type(grid_def), intent(in) :: grid
    real(dp), intent(in)  :: y
    integer,  intent(out) :: iylo
    real(dp), pointer     :: wgts(:)
    !-----------------------------------------
    integer, parameter :: npnt_min = 4, npnt_max = 10
    integer :: ny, npnt, isub

    ny = grid%ny
    if (grid%nsub /= 0) then
       isub = conv_BestIsub(grid,y)
       call WgtGridQuant(grid%subgd(isub), y, iylo, wgts)
       iylo = iylo + grid%subiy(isub)
    else
       if (y > grid%ymax*(one+warn_tolerance) .or. y < -warn_tolerance) then
          write(0,*) 'WgtGridQuant: &
               &requested function value outside y range'
          write(0,*) 'y = ', y, ' but should be 0 < y < ymax=',grid%ymax
          stop
       end if
       
       npnt = min(npnt_max, max(npnt_min, abs(grid%order)))
       allocate(wgts(0:npnt-1))
       
       iylo = min(max(floor(y / grid%dy)-(npnt-1)/2,0),ny-npnt+1)
       call uniform_interpolation_weights(y/grid%dy-iylo, wgts)
    end if
  end subroutine WgtGridQuant
  

  

  !-- for internal use only
  function conv_BestIsub(grid,y) result(isub)
    type(grid_def), intent(in) :: grid
    real(dp),       intent(in) :: y
    integer                    :: isub
       !-- this will probably slow things down, but
       !   for the time being accept this penalty
       ! find the grid with the smallest ymax > y
       if (y>grid%ymax) then
          isub = sum(maxloc(grid%subgd%ymax))
       else
          isub = sum(minloc(grid%subgd%ymax,mask=(grid%subgd%ymax>=y)))
       end if
  end function conv_BestIsub
  
  !----------------------------------------------------------------------
  ! not yet updated for multi-grid purposes, because would require
  ! a little bit of thought as to how to treat different grid regions
  function MomGridQuant(grid,gq,omega) result(res)
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq(0:)
    real(dp), intent(in) :: omega
    real(dp) :: res
    !-----------------------------------------
    real(dp) :: weight, weightprod, dy
    integer :: i, ny

    if (grid%nsub /= 0) then
       write(0,*) 'ERROR in MomGridQuant:&
            & multiple grids not yet supported'
    end if
    
    ny = assert_eq(grid%ny,ubound(gq,dim=1),"MomGridQuant")
    dy = grid%dy
    if (omega == zero) then
       weight = grid%dy
       weightprod = one
       res = half*weight*gq(0)
    else
       weightprod = exp(-dy*omega)
       weight = (exp(-dy*omega) - one + dy*omega)/(dy*omega**2)
       res = weight*gq(0)
       weight = weight + &
            & (exp(+dy*omega) - one - dy*omega)/(dy*omega**2)
    end if
    do i = 1, ny
       weight = weight * weightprod
       if (i == ny) weight = weight * half
       res = res + weight*gq(i)
    end do
  end function MomGridQuant
  

  !----------------------------------------------------------------------
  ! Use a non-sophisticated fix for the multi-grid option, just take the
  ! largest of the dy values if none is specified (so as to avoid information
  ! overload).
  subroutine conv_PrintGridQuant_1(grid,gq,dy,dev)
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq(0:)
    real(dp), intent(in), optional :: dy
    integer,  intent(in), optional :: dev
    real(dp) :: dy_local, y, q
    integer :: ny, i, dev_local
    
    ny = assert_eq(grid%ny,ubound(gq,dim=1),'PrintGridQuant')
    if (grid%nsub /= 0) then
       dy_local = default_or_opt(maxval(grid%subgd%dy),dy)
    else
       dy_local = default_or_opt(grid%dy,dy)
    end if
    dev_local = default_or_opt(6,dev)

    ny = floor(grid%ymax / dy_local)
    do i = 0, ny
       y = i*dy_local
       q = EvalGridQuant(grid, gq, y)
       write(dev_local,*) y, exp(-y),q
    end do
    
  end subroutine conv_PrintGridQuant_1

  !----------------------------------------------------------------------
  ! See conv_PrintGridQuant_1 re multigrid
  subroutine conv_PrintGridQuant_2(grid,gq,gq2,dy,dev)
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq(0:),gq2(0:)
    real(dp), intent(in), optional :: dy
    integer,  intent(in), optional :: dev
    real(dp) :: dy_local, y,  q,q2
    integer :: ny, i, dev_local
    
    ny = assert_eq(grid%ny,ubound(gq,dim=1),&
         & ubound(gq2,dim=1),'PrintGridQuant')
    if (grid%nsub /= 0) then
       dy_local = default_or_opt(maxval(grid%subgd%dy),dy)
    else
       dy_local = default_or_opt(grid%dy,dy)
    end if
    dev_local = default_or_opt(6,dev)

    ny = floor(grid%ymax / dy_local)
    do i = 0, ny
       y = i*dy_local
       q  = EvalGridQuant(grid, gq, y)
       q2 = EvalGridQuant(grid, gq2, y)
       write(dev_local,'(25es25.16)') y, exp(-y),q, q2
    end do
    
  end subroutine conv_PrintGridQuant_2

  !----------------------------------------------------------------------
  ! See conv_PrintGridQuant_1 re multigrid
  subroutine conv_PrintGridQuant_3(grid,gq,gq2,gq3,dy,dev)
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq(0:),gq2(0:),gq3(0:)
    real(dp), intent(in), optional :: dy
    integer,  intent(in), optional :: dev
    real(dp) :: dy_local, y, x, q,q2,q3
    integer :: ny, i, dev_local
    
    ny = assert_eq(grid%ny,ubound(gq,dim=1),&
         & ubound(gq2,dim=1),ubound(gq3,dim=1),&
         & 'PrintGridQuant: distributions must be same size')
    if (grid%nsub /= 0) then
       dy_local = default_or_opt(maxval(grid%subgd%dy),dy)
    else
       dy_local = default_or_opt(grid%dy,dy)
    end if
    dev_local = default_or_opt(6,dev)

    ny = floor(grid%ymax / dy_local)
    do i = 0, ny
       y = i*dy_local
       q  = EvalGridQuant(grid, gq, y)
       q2 = EvalGridQuant(grid, gq2, y)
       q3 = EvalGridQuant(grid, gq3, y)
       write(dev_local,'(25es25.16)') y, exp(-y),q, q2, q3
    end do
    
  end subroutine conv_PrintGridQuant_3
  
  !----------------------------------------------------------------------
  ! See conv_PrintGridQuant_1 re multigrid
  subroutine conv_PrintGridQuant_4(grid,gq,gq2,gq3,gq4,dy,dev)
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq(0:),gq2(0:),gq3(0:),gq4(0:)
    real(dp), intent(in), optional :: dy
    integer,  intent(in), optional :: dev
    real(dp) :: dy_local, y, q,q2,q3, q4
    integer :: ny, i, dev_local
    
    ny = assert_eq(grid%ny,ubound(gq,dim=1),&
         & ubound(gq2,dim=1),ubound(gq3,dim=1),ubound(gq4,dim=1),&
         & 'PrintGridQuant: distributions must be same size')
    if (grid%nsub /= 0) then
       dy_local = default_or_opt(maxval(grid%subgd%dy),dy)
    else
       dy_local = default_or_opt(grid%dy,dy)
    end if
    dev_local = default_or_opt(6,dev)

    ny = floor(grid%ymax / dy_local)
    do i = 0, ny
       y = i*dy_local
       q  = EvalGridQuant(grid, gq, y)
       q2 = EvalGridQuant(grid, gq2, y)
       q3 = EvalGridQuant(grid, gq3, y)
       q4 = EvalGridQuant(grid, gq4, y)
       write(dev_local,'(25es25.16)') y, exp(-y),q, q2, q3, q4
    end do
    
  end subroutine conv_PrintGridQuant_4


  !======================================================================
  ! Routines related to calculations of moments
  !======================================================================

  !---------------------------------------------------------------------
  !! Computes the truncated moment between 0 and min(y,ymax)
  !! of the grid quantity gq with the moment index defined as
  !!
  !! Define the function which is going to be integrated
  !! Multiply the pdf by exp( - y * moment_index)
  !! Moment index are defined as
  !! 
  !!  M(moment_index) = \int_0^ymax dy exp( - y * moment_index) * xpdf(y)
  !! 
  !! or in terms of x, where y = ln(1/x)
  !!
  !!   M(moment_index) = \int_xmin^1 dx x^( moment_index - 1) * xpdf(x)
  !! 
  !! where xpdf(y) is gq .aty. (y .with. gd)
  !-------------------------------------------------------------
  function conv_TruncatedMoment_1d(gd, gq, moment_index, y) result(res)
    use integrator
    type(grid_def), intent(in) :: gd
    real(dp),       intent(in) :: gq(0:), moment_index
    real(dp), intent(in), optional :: y
    real(dp)                   :: res, ymax,yl,yh,eps
    integer :: ny
    !-----------------------------------------------------
    
    res = zero
    
    ! global variables to pass to auxiliary conv_TruncatedMoment_helper
    ! inefficient to make copies, but works around impossibility
    ! of having pointers to intent(in) objects???
    conv_moment_index = moment_index
    conv_moment_gd    = gd
    allocate( conv_moment_gq(0:gd%ny) )
    conv_moment_gq    = gq
    
    ! Check dimensions are correct
    ! write(6,*) gd%ny, ubound(gq,dim=1)
    ny = assert_eq(gd%ny,ubound(gq,dim=1),"EvalGridQuant")
    if( gd%ny .ne. ubound(gq,dim=1) ) &
         & call wae_error("Different dimensions for ",&
         & "grid quantity and grid in conv_TruncatedMoment_1d")

  
    ! Get ymax, upper integration limit
    ymax = default_or_opt(gd%ymax, y)
    ! Warning: the upper integration limit cannot be larger than
    ! ymax from the grid definition
    !write(6,*) ymax
    if (ymax > gd%ymax) call wae_error("Input value for y",&
         & "for computation of truncated moment",&
         & "larger than ymax from grid")
    
    ! Set convolution for integration
    eps = default_conv_eps

    ! Integrate between 0 and ymax
    ! Do the integration in two steps to improve accuracy
    ! Check ymax > two, the indermediate limit
    yl=zero
    yh=ymax
    res = ig_LinWeight(conv_TruncatedMoment_helper,&
         & yl,min(two,ymax),one,one,eps)
    res = res +  ig_LinWeight(conv_TruncatedMoment_helper,&
         & min(two,ymax),ymax, one,one,eps)

    ! Dealocate auxiliary grid quantity
    deallocate( conv_moment_gq )

  end function conv_TruncatedMoment_1d

  
  !------------------------------------------------
  ! Auxiliary function for conv_TruncatedMoment_1d
  !------------------------------------------------
  function conv_TruncatedMoment_helper(y) result(res)
    real(dp), intent(in) :: y 
    real(dp) :: res

    res = exp( -y * conv_moment_index ) &
         & * EvalGridQuant(conv_moment_gd,conv_moment_gq,y) 
   
  end function conv_TruncatedMoment_helper
  

  !======================================================================
  ! Routines related to convolutions
  !======================================================================

  !----------------------------------------------------------------------
  !! set the default integration precision
  subroutine SetDefaultConvolutionEps(eps)
    real(dp), intent(in) :: eps
    default_conv_eps = eps
  end subroutine SetDefaultConvolutionEps

  !----------------------------------------------------------------------
  !! set the default integration precision
  real(dp) function DefaultConvolutionEps() result(res)
    res = default_conv_eps
  end function DefaultConvolutionEps
  

  !------------------------------------------------------------------
  ! Just does the memory allocation
  recursive subroutine conv_AllocGridConv_0d(grid,gc)
    type(grid_def),  intent(in)    :: grid
    type(grid_conv), intent(inout) :: gc
    !-------------------------------------------
    integer :: isub
    
    gc%grid = grid
    if (gc%grid%nsub /= 0) then
       nullify(gc%conv) ! avoid doing anything too stupid...
       allocate(gc%subgc(gc%grid%nsub))
       do isub = 1, gc%grid%nsub
          call conv_AllocGridConv_0d(grid%subgd(isub),gc%subgc(isub))
       end do
    else
       nullify(gc%subgc)
       ! this form of deallocate ought to be OK (?), but on lahey,
       ! compaq+condor (and perhaps elsewhere) it causes problems
       !deallocate(gc%conv,stat=istat)
       if (grid%order == LIN_ORDER) then
          allocate(gc%conv(0:grid%ny,2))
       else if (grid%order > LIN_ORDER) then
          allocate(gc%conv(0:grid%ny,0:grid%order+1))
       else
          allocate(gc%conv(0:grid%ny,1))
       end if
    end if
  end subroutine conv_AllocGridConv_0d

  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_AllocGridConv_1d(grid,gc)
    type(grid_def),  intent(in)    :: grid
    type(grid_conv), intent(inout) :: gc(:)
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_AllocGridConv_0d(grid,gc(i))
    end do
  end subroutine conv_AllocGridConv_1d
  subroutine conv_AllocGridConv_2d(grid,gc)
    type(grid_def),  intent(in)    :: grid
    type(grid_conv), intent(inout) :: gc(:,:)
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_AllocGridConv_0d(grid,gc(j,i))
       end do
    end do
  end subroutine conv_AllocGridConv_2d


  !======================================================================
  !! Return .true. if the gc is currently allocated, .false. otherwise
  function GridConvAllocated(gc)
    logical :: GridConvAllocated
    type(grid_conv), intent(in) :: gc
    GridConvAllocated = associated(gc%conv) .or. associated(gc%subgc)
    !GridConvAllocated = .false.
  end function GridConvAllocated
  
  !--------------------------------------------------------------------
  ! initialise a grid convolution with zero
  ! Default for alloc is .true.; writing this this way allows
  ! in principle the possibility of a more levels of recursion
  ! in the grid def, should one ever want to have the option...
  recursive subroutine conv_InitGridConv_zero(grid,gc,alloc)
    type(grid_def),  intent(in)    :: grid
    type(grid_conv), intent(inout) :: gc
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: isub

    if (default_or_opt(.not.GridConvAllocated(gc),alloc)) then
       call AllocGridConv(grid,gc)
    else
       call ValidateGD(grid, gc%grid, 'conv_InitGridConv_zero')
    end if

    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          ! remember not to do double allocate
          call conv_InitGridConv_zero(grid%subgd(isub),&
               &gc%subgc(isub),alloc=.false.)
       end do
    else
       gc%conv = zero
    end if
    
  end subroutine conv_InitGridConv_zero

  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_InitGridConv_zero_1d(grid,gc,alloc)
    type(grid_def),  intent(in)    :: grid
    type(grid_conv), intent(inout) :: gc(:)
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_InitGridConv_zero(grid,gc(i),alloc)
    end do
  end subroutine conv_InitGridConv_zero_1d
  subroutine conv_InitGridConv_zero_2d(grid,gc,alloc)
    type(grid_def),  intent(in)    :: grid
    type(grid_conv), intent(inout) :: gc(:,:)
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_InitGridConv_zero(grid,gc(j,i),alloc)
       end do
    end do
  end subroutine conv_InitGridConv_zero_2d
  

  !--------------------------------------------------------------------
  ! initialise a grid convolution with another gc, potentially multiplied
  ! by a factor. Qu: is alloc supposed to be used from outside?
  recursive subroutine conv_InitGridConv_gc(gc,gcin,fact,alloc)
    type(grid_conv), intent(inout) :: gc
    type(grid_conv), intent(in)    :: gcin
    real(dp),        intent(in), optional :: fact
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: isub

    if (default_or_opt(.not.GridConvAllocated(gc),alloc)) then
       call AllocGridConv(gcin%grid,gc)
    else
       call ValidateGD(gcin%grid,gc%grid,'conv_InitGridConv_gc')
    end if

    if (gcin%grid%nsub /= 0) then
       do isub = 1, gcin%grid%nsub
          ! remember not to do double allocate
          call conv_InitGridConv_gc(gc%subgc(isub),&
               &gcin%subgc(isub),fact=fact,alloc=.false.)
       end do
    else
       if (present(fact)) then
          gc%conv = gcin%conv * fact
       else
          gc%conv = gcin%conv
       end if
    end if
  end subroutine conv_InitGridConv_gc

  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_InitGridConv_gc_1d(gc,gcin,fact,alloc)
    type(grid_conv), intent(inout) :: gc(:)
    type(grid_conv), intent(in)    :: gcin(:)
    real(dp),        intent(in), optional :: fact
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: i, ni
    ni = assert_eq(size(gc),size(gcin),'conv_InitGridConv_gc_1d')
    do i = 1, ni
       call conv_InitGridConv_gc(gc(i),gcin(i),fact,alloc)
    end do
  end subroutine conv_InitGridConv_gc_1d
  subroutine conv_InitGridConv_gc_2d(gc,gcin,fact,alloc)
    type(grid_conv), intent(inout) :: gc(:,:)
    type(grid_conv), intent(in)    :: gcin(:,:)
    real(dp),        intent(in), optional :: fact
    logical,         intent(in), optional :: alloc
    !-------------------------------------------
    integer :: i,j,ni,nj
    ni = assert_eq(size(gc,dim=2),size(gcin,dim=2),'conv_InitGridConv_gc_1d')
    nj = assert_eq(size(gc,dim=1),size(gcin,dim=1),'conv_InitGridConv_gc_1d')
    do i = 1, ni
       do j = 1, nj
          call conv_InitGridConv_gc(gc(j,i),gcin(j,i),fact,alloc)
       end do
    end do
  end subroutine conv_InitGridConv_gc_2d


  !----------------------------------------------------------------------
  ! Initialise a convoluter with the function to use in the 
  ! convolution. 
  subroutine conv_InitGridConv_func(grid,gc,func,alloc)
    type(grid_def),  intent(in)    :: grid
    type(grid_conv), intent(inout) :: gc
    interface
       function func(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: func
       end function func
    end interface
    logical,         intent(in), optional :: alloc
    !-------------------------------------------

    !if (default_or_opt(.not.GridConvAllocated(gc),alloc)) call conv_InitGridConv_zero(grid,gc)
    call conv_InitGridConv_zero(grid,gc)
    
    call AddWithCoeff(gc,func)
  end subroutine conv_InitGridConv_func


  !----------------------------------------------------------------------
  ! Initialise a convoluter with the function to use in the 
  ! convolution. 
  subroutine conv_InitGridConv_conv(gc,gca,gcb,alloc)
    type(grid_conv),  intent(inout) :: gc
    type(grid_conv),  intent(in)    :: gca, gcb
    logical,          intent(in), optional :: alloc
    !-------------------------------------------

    if (default_or_opt(.not.GridConvAllocated(gc),alloc)) then
       call AllocGridConv(gca%grid,gc)
    else
       call ValidateGD(gc%grid, gca%grid,&
            & 'conv_InitGridConv_conv: gc and gca')
    end if
    
    call SetToConvolution(gc,gca, gcb)
  end subroutine conv_InitGridConv_conv

  !------------------------------------------------------------------
  ! Zero the contents of the grid convolution
  !
  recursive subroutine conv_ZeroGridConv_0d(gc)
    type(grid_conv), intent(inout) :: gc
    integer :: isub
    if (gc%grid%nsub /= 0) then
       do isub = 1, gc%grid%nsub
          call SetToZero(gc%subgc(isub))
       end do
    else
       gc%conv = zero
    end if
  end subroutine conv_ZeroGridConv_0d
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_ZeroGridConv_1d(gc)
    type(grid_conv), intent(inout) :: gc(:)
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_ZeroGridConv_0d(gc(i))
    end do
  end subroutine conv_ZeroGridConv_1d
  subroutine conv_ZeroGridConv_2d(gc)
    type(grid_conv), intent(inout) :: gc(:,:)
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_ZeroGridConv_0d(gc(j,i))
       end do
    end do
  end subroutine conv_ZeroGridConv_2d
  
  !--------------------------------------------------------
  ! Remove memory associated with a given gc.
  recursive subroutine conv_DelGridConv_0d(gc)
    type(grid_conv), intent(inout) :: gc
    !-------------------------------------------
    integer :: istat, isub
    
    if (.not. GridConvAllocated(gc)) return

    if (gc%grid%nsub /= 0) then
       do isub = 1, gc%grid%nsub
          call conv_DelGridConv_0d(gc%subgc(isub))
       end do
       deallocate(gc%subgc,stat=istat)
       nullify(gc%subgc)  ! to make sure it isn't left dangling (necessary?)
    else
       deallocate(gc%conv,stat=istat)
       nullify(gc%conv) ! to make sure it isn't left dangling (necessary?)
    end if

    if (istat /= 0) then
       write(0,*) 'ERROR: problems with deallocation in conv_DelGridConv_0d'
       !write(0,*) one/sqrt(sin(zero))
       stop
    end if
    
  end subroutine conv_DelGridConv_0d
  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_DelGridConv_1d(gc)
    type(grid_conv), intent(inout) :: gc(:)
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_DelGridConv_0d(gc(i))
    end do
  end subroutine conv_DelGridConv_1d
  subroutine conv_DelGridConv_2d(gc)
    type(grid_conv), intent(inout) :: gc(:,:)
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_DelGridConv_0d(gc(j,i))
       end do
    end do
  end subroutine conv_DelGridConv_2d


  !--------------------------------------------------------------------
  ! Multiply the given GridConv by the relevant factor
  recursive subroutine conv_MultGridConv_0d(gc,fact)
    type(grid_conv), intent(inout) :: gc
    real(dp),        intent(in)    :: fact
    !-------------------------------------------
    integer :: isub

    if (gc%grid%nsub /= 0) then
       do isub = 1, gc%grid%nsub
          call Multiply(gc%subgc(isub),fact)
       end do
    else
       gc%conv = gc%conv * fact
    end if
    
  end subroutine conv_MultGridConv_0d
  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_MultGridConv_1d(gc,fact)
    type(grid_conv), intent(inout) :: gc(:)
    real(dp),        intent(in)    :: fact
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_MultGridConv_0d(gc(i),fact)
    end do
  end subroutine conv_MultGridConv_1d
  subroutine conv_MultGridConv_2d(gc,fact)
    type(grid_conv), intent(inout) :: gc(:,:)
    real(dp),        intent(in)    :: fact
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_MultGridConv_0d(gc(j,i),fact)
       end do
    end do
  end subroutine conv_MultGridConv_2d

  

  !-------------------------------------------------------------
  ! To gc add a function for convolution
  recursive subroutine conv_AddGridConv_func(gc,func)
    use integrator; use convolution_communicator
    type(grid_conv), intent(inout), target :: gc
    interface
       function func(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: func
       end function func
    end interface
    !----------------------------------------------------
    real(dp), pointer :: dy
    integer,  pointer :: ny
    integer  :: order
    real(dp) :: upper,lower, yl,ym,yh
    integer  :: i,j, k
    !------------------------------------------------------
    real(dp) :: res,eps!, nodes(gc%grid%order+1)
    integer  :: il, ih, iy, jy
    real(dp), allocatable :: nodes(:)
    !------------------------------------------------------
    integer :: isub

    if (gc%grid%nsub /= 0) then
       do isub = 1, gc%grid%nsub
          call conv_AddGridConv_func(gc%subgc(isub),func)
       end do
       return
    end if
    
    dy => gc%grid%dy
    ny => gc%grid%ny
    order = gc%grid%order
    eps   = gc%grid%eps

    !-- this used to be an automatic, but ran into problems 
    !   with multi-grid splitting functions because compound holder
    !   had an undefined value for gc%grid%order
    allocate(nodes(abs(gc%grid%order)+1))

    cc_piece = cc_REAL
    if (gc%grid%order == LIN_ORDER) then
       do i = 1, ny
          ym = i*dy
          yl = ym - dy
          yh = ym + dy
          lower = ig_LinWeight(func,yl,ym,zero,one,eps) ! specify 0 at yl
          upper = ig_LinWeight(func,ym,yh,one,zero,eps) ! specify 0 at yh
          gc%conv(i,FULL) = gc%conv(i,FULL) + lower + upper
          gc%conv(i,UPPR) = gc%conv(i,UPPR) + upper
       end do
       !-- see CCN17-62 ----
       ym = zero
       yh = dy
       cc_piece = cc_REAL
       upper    = ig_LinWeight(func,ym,yh,zero,-one,eps)
       cc_piece = cc_REALVIRT
       upper    = upper + ig_LinWeight(func,ym,yh,one,one,eps) 
       cc_piece = cc_VIRT
       upper    = upper + ig_LinWeight(func,yh,-two*log(eps),one,one,eps)
       cc_piece = cc_DELTA
       upper    = upper + func(zero)
       gc%conv(0,FULL) = gc%conv(0,FULL) + upper
       gc%conv(0,UPPR) = gc%conv(0,UPPR) + upper
    else if (gc%grid%order < LIN_ORDER) then
       ! should be similar to Ratcliffes proposal (and the sort of thing
       ! regularly done in BFKL). NB Not documented in any CCN -- hopefully
       ! straightforward enough that it can be done in one's head?
       order = -gc%grid%order
       do i = 1, ny
          !-- this is the range of interest
          yl = (i-1) * dy
          yh =  i    * dy
          !-- this is range of nodes for that interval
          il = i-1
          ih = il + order
          nodes = (/ (j,j=il,ih) /) * dy
          do iy = il, min(ih,ny)
             res = conv_GCAf_Helper(func,yl,yh,il,iy,nodes,eps)
             gc%conv(iy,1) = gc%conv(iy,1) + res
          end do
       end do
       !write(0,*) 'conv_AddGridConv_func: Negative orders not yet supported'
       !stop
    else
       !-- CCN19 p.4 -------------
       ! index 0:          central points
       ! index 1..order+1: for dealing with last order+1 points
       !-- first do central points. Do an interval at a time
       !   i = 1, means interval from 0 to 1
       !-- the first loop is missing something:: it should also be filling up
       !   some pieces in the second loop, at least in some cases
       do i = 1, ny
          !-- this is the range of interest
          yl = (i-1) * dy
          yh =  i    * dy
          !-- these are the interpolation points to be used
          il = max(0,i - (order+2)/2)
          ih = min(ny, il+order)
          il = ih - order
          !-- fill things up?
          nodes = (/ (j,j=il,ih) /) * dy
          do iy = il, ih
             res = conv_GCAf_Helper(func,yl,yh,il,iy,nodes,eps)
             gc%conv(iy,0) = gc%conv(iy,0) + res
             !-- deal properly with special end array --
             !   try to be clever with limits, but there are no guarantees.
             !   It does however seem to work!
             do jy = max(i+order,iy), iy+order
                if (jy <= ny) then
                   gc%conv(jy,order+1-(jy-iy)) &
                        &= gc%conv(jy,order+1-(jy-iy)) + res
                end if
             end do
          end do
       end do
       !-- now deal with integrations close to end points
       ! k represents the end points
       do k = 1, ny
          if (k <= order) then
             yl = max(zero, (k-order)*dy)
             yh = k*dy
             !-- the interpolation points & yes, it is normal for them 
             !   to be negative
             !   it is related to the way the convolution is done for the case 
             !   i < order, i.e. effectively conv(0:order)*gq(order:0:-1)
             !   (actually conv(i, 1:order+1))
             ih = k
             il = k - order
             !-- fill things up
             nodes = (/ (j,j=il,ih) /) * dy
             do iy = il, ih
                res = conv_GCAf_Helper(func,yl,yh,il,iy,nodes,eps)
                gc%conv(k,iy-il+1) = gc%conv(k,iy-il+1) + res
             end do
             cycle
          end if
          ! now do things region by region
          !-- i represents the region that we will study
          do i = max(1,k-order+1), k
             !-- this is the range of interest
             yl = (i-1) * dy
             yh =  i    * dy
             !-- these are the interpolation points to be used
             !   note extra min(k,...) compared to above formula
             il = max(0,i - (order+2)/2)
             ih = min(k,min(ny, il+order)) !-- could simplify expression
             il = ih - order
             nodes = (/ (j,j=il,ih) /) * dy
             do iy = max(il,k-order), ih
                res = conv_GCAf_Helper(func,yl,yh,il,iy,nodes,eps)
                gc%conv(k,order+1-(k-iy)) = &
                     &gc%conv(k,order+1-(k-iy)) + res
             end do
          end do
       end do
    end if
    
    deallocate(nodes)

  end subroutine conv_AddGridConv_func


  !---------------------------------------------------------------------
  ! look after some repetitive work...
  !
  ! Guess that this function may do the following:
  !   Work out the weight corresponding to the integral, between yl and yh
  !   of the function func(y) mutiplied by a polynomial which is zero at 
  !   all the nodes (which start from il) except that indicated by iy.
  !
  function conv_GCAf_Helper(func,yl,yh,il,iy,nodes,eps) result(res)
    use integrator; use convolution_communicator
    interface
       function func(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: func
       end function func
    end interface
    real(dp), intent(in) :: yl, yh, nodes(:), eps
    integer,  intent(in) :: il, iy
    real(dp)             :: res
    integer :: inode
    inode = iy - il + 1
    res = zero
    if (yl == zero .and. iy == 0) then
       cc_piece = cc_VIRT
       res = res + ig_LinWeight(func, yh,-two*log(eps), one,one, eps)
       cc_piece = cc_DELTA
       res = res + func(zero)
       cc_piece = cc_REALVIRT
       res = res + ig_LinWeight(func, yl, yh, one, one, eps)
       cc_piece = cc_REAL
       res = res + ig_PolyWeight(func, yl, yh, nodes, inode, eps,wgtadd=-one)
    else
       cc_piece = cc_REAL
       res = ig_PolyWeight(func, yl, yh, nodes, inode, eps)
    end if
  end function conv_GCAf_Helper
  


  !-------------------------------------------------------------
  ! To gc add another grid convolution.
  !
  ! Perhaps there are some inefficiences here (validation may sometimes
  ! be carried out twice), but for the time being, leave it as is.
  recursive subroutine conv_AddGridConv_gc(gc,gcadd,fact)
    type(grid_conv), intent(inout) :: gc
    type(grid_conv), intent(in)  :: gcadd
    real(dp),        intent(in), optional  :: fact
    !----------------------------------------------
    integer :: isub

    call ValidateGD(gc%grid, gcadd%grid, 'conv_AddGridConv_gc')
    if (gc%grid%nsub /= 0) then
       do isub = 1, gc%grid%nsub
          call conv_AddGridConv_gc(gc%subgc(isub),gcadd%subgc(isub),fact)
       end do
    else
       if (present(fact)) then
          gc%conv = gc%conv + gcadd%conv * fact
       else
          gc%conv = gc%conv + gcadd%conv
       end if
    end if
    
  end subroutine conv_AddGridConv_gc
  !-----------------------------------------------------------------
  ! Handy to have multi-dimensional versions of the above
  subroutine conv_AddGridConv_gc_1d(gc,gcadd,fact)
    type(grid_conv), intent(inout) :: gc(:)
    type(grid_conv), intent(in)    :: gcadd(:)
    real(dp),        intent(in), optional :: fact
    !-------------------------------------------
    integer :: i
    do i = 1, size(gc)
       call conv_AddGridConv_gc(gc(i),gcadd(i),fact)
    end do
  end subroutine conv_AddGridConv_gc_1d
  subroutine conv_AddGridConv_gc_2d(gc,gcadd,fact)
    type(grid_conv), intent(inout) :: gc(:,:)
    type(grid_conv), intent(in)    :: gcadd(:,:)
    real(dp),        intent(in), optional    :: fact
    !-------------------------------------------
    integer :: i,j
    do i = 1, size(gc,dim=2)
       do j = 1, size(gc,dim=1)
          call conv_AddGridConv_gc(gc(j,i),gcadd(j,i),fact)
       end do
    end do
  end subroutine conv_AddGridConv_gc_2d


  !--------------------------------------------------------------
  ! Carry out the convolution of gc on gq
  recursive function conv_ConvGridQuant_scalar(gc,gq) result(gqout)
    type(grid_conv),  intent(in)    :: gc
    real(dp),         intent(in)    :: gq(0:)
    real(dp)                        :: gqout(0:ubound(gq,dim=1))
    !---------------------------------------------
    integer :: i, ny !, j
    integer :: order
    integer :: isub, iy, dy_ratio

    if (gc%grid%nsub /= 0) then
       do isub = 1, gc%grid%nsub
          gqout(gc%grid%subiy(isub):gc%grid%subiy(isub+1)-1) = &
               & conv_ConvGridQuant_scalar(gc%subgc(isub),&
               &        gq(gc%grid%subiy(isub):gc%grid%subiy(isub+1)-1))
       end do
        
       if (gc%grid%locked .and. .not.override_grid_locking) then
          !-- decant information from finer grids into coarser grids
          ! (remember: finest grid has lowest isub)
          do isub = 2, gc%grid%nsub
             ! the ratio should be an exact integer, but use
             ! nint() to avoid the dangers of rounding errors
             dy_ratio = nint(gc%grid%subgd(isub)%dy / gc%grid%subgd(isub-1)%dy)
             do iy = 0, gc%grid%subgd(isub-1)%ny / dy_ratio
                gqout(gc%grid%subiy(isub)+iy) = &
                     &gqout(gc%grid%subiy(isub-1)+iy*dy_ratio)
             end do
          end do
          nconv_with_override_off = nconv_with_override_off + 1
       end if
       return
    end if

    ny = assert_eq(gc%grid%ny,ubound(gq,dim=1),"conv_ConvGridQuant")
    order = gc%grid%order

    !-- Hopefully this will avoid some wasted convolutions?
    if (all(gq == zero)) then
       gqout = zero
       return
    end if
    

    if (order == LIN_ORDER) then
       !-- a test to avoid N^2 operations...
       if (all(gc%conv(:,FULL) == zero)) then
          gqout = zero
          return
       end if
       
       gqout(0) = zero
       do i = 1, ny
          !-- following is legal fortran but it is very slow sometimes
          gqout(i) =  sum(gq(0:i)*gc%conv(i:0:-1,FULL))
          !-- maybe the version that follows will be faster on some compilers
          !   on absoft it definitely is
!!$       gqout(i) = zero
!!$       do j = 0, i
!!$          gqout(i) = gqout(i) + gq(j)*gc%conv(i-j,FULL)
!!$       end do
          !write(0,*) i*gc%grid%dy,gqout(i), gc%conva(1)
          gqout(i) = (gqout(i)-gq(0)*gc%conv(i,UPPR))
       end do
    else if (order < 0) then
       ! Ratcliffes proposal
       !gqout(0) = zero
       ! NB: explicitly do also the i=0 case. Under normal circumstances
       !     it should give zero; however we still let the user enjoy
       !     their folly if they really want to have gq(0)/=0 --- this
       !     is something that comes in handy when generating "derived"
       !     convolutors...
       do i = 0, ny
          gqout(i) = sum(gc%conv(0:i,1)*gq(i:0:-1))
       end do
    else
       gqout = zero
       !-- a test to avoid N^2 operations...
       if (all(gc%conv(:,0) == zero)) return
       !-- current writing is designed to reduce thrashing of cache. 
       !   It is not clear that it actually helps in any way
       !   Commented version theoretically has higher cache thrashing.
       do i = 1, ny
          if (i > order) then
             gqout(i) = sum(gc%conv(0:i-order-1,0)*gq(i:order+1:-1))
          end if
          !-- HIGH CACHE THRASH
          !gqout(i) = gqout(i) + sum(gc%conv(i,1:order+1)*gq(order:0:-1))
       end do
       !-- LOW CACHE THRASH
       do i = 1, order+1
          gqout(:) = gqout(:) + gc%conv(:,i)*gq(order+1-i)
       end do
    end if
    
  end function conv_ConvGridQuant_scalar



  !--------------------------------------------------------------
  ! Carry out the convolution of two convolution functions
  ! It is recommended that if one of them is singular at x=1, then this be 
  ! gca: this helps minimise the discretisation error (in some cases).
  !
  ! NB: subroutine name is not quite consistent with things like
  !     Add or Mult since the operands are completely distinct from the
  !     result.
  recursive subroutine conv_ConvGridConv_0d(gc,gca, gcb, reorder)
    type(grid_conv),  intent(inout) :: gc
    type(grid_conv),  intent(in)    :: gca, gcb
    logical,          intent(in), optional :: reorder
    !---------------------------------------------
    integer :: i, j, order, ix
    !-- these might be better if on the fly, given recursive nature?
    real(dp) :: deltafn(0:gca%grid%ny), res(0:gca%grid%ny)
    integer :: isub

    call ValidateGD(gca%grid, gcb%grid, 'conv_ConvGridConv_0d: gca and gcb')
    if (.not.GridConvAllocated(gc)) then
       call AllocGridConv(gca%grid,gc)
    else
       call ValidateGD(gc%grid, gca%grid, 'conv_ConvGridConv_0d: gc and gca')
    end if
    
    if (gc%grid%nsub /= 0) then
       do isub = 1, gc%grid%nsub
          call conv_ConvGridConv_0d(gc%subgc(isub),&
               &  gca%subgc(isub), gcb%subgc(isub),reorder)
       end do
       return
    end if

    if (default_or_opt(.true.,reorder)) then
       !-- decide which part of conv matrix to use
       select case (gca%grid%order)
       case(LIN_ORDER)
          ix = FULL
       case(LIN_ORDER+1:)
          ix = 0
       case(:LIN_ORDER-1)
          ix = 1
       end select
       
       !-- use as only criterion of singularity that there should be
       !   opposite signs for the first two grid points? This is not
       !   by any means foolproof, but in basic cases it seems to work.
       !   Put the less singular convolutor to the right.
       if (product(gca%conv(0:1,ix)) >= zero .and. &
            & product(gcb%conv(0:1,ix)) < zero) then
          call conv_ConvGridConv_0d(gc,gcb, gca, reorder=.false.)
          return
       end if
    end if


    if (gca%grid%order == LIN_ORDER) then
       !-- this is just a guess. Now we will see if it is right
       gc%conv(:,FULL) = (gca .conv. gcb%conv(:,FULL)) +&
            &gca%conv(:,UPPR) * gcb%conv(0,FULL)
       
       !-- next 3 lines actually calculate LOWER. Then correct it to upper.
       !   Do the convolution of a deltafn by hand in order to
       !   be more efficient: tested that it gave identical answers
       deltafn = gcb%conv(:,FULL) - gcb%conv(:,UPPR)
       deltafn(0) = zero
       gc%conv(:,UPPR) = gca .conv. deltafn
       !gc%conv(:,UPPR) = gca .conv. (gcb .conv. deltafn)
       gc%conv(:,UPPR) = gc%conv(:,FULL) - gc%conv(:,UPPR)
    else if (gca%grid%order < LIN_ORDER) then
       ! the ratcliffe approach
       do i = 0, gc%grid%ny
          gc%conv(i,1) = sum(gca%conv(0:i,1)*gcb%conv(i:0:-1,1))
       end do
       !write(0,*) 'SetToConvolution: Negative orders not yet supported'
       !stop
    else
       order = gca%grid%order
       !-- let us hope that the following works?
       !   in this first incarnation it will be roughly half the optimum
       !   speed because the inner convolution will be done "stupidly",
       !   using O(N^2) steps rather than the possible O(N)
       do j = 0, order+1
          !-- set up a delta function and its convolution
          deltafn = zero; deltafn(j) = one
          !-- try to avoid too much rubbish in the result?
          res = zero
          res = gca .conv. (gcb .conv. deltafn)
          !-- redistribute the result
          if (j /= order+1) then
             gc%conv(:,order+1-j) = res(:)
          else
             gc%conv(0:gca%grid%ny-(order+1),0) = res(order+1:gca%grid%ny)
             !-- The region conv(gca%grid%ny-order:) remains unfilled. This 
             !   is OK since in practice (current incarnation) it never 
             !   gets used.
          end if
       end do
       !-- now remember to test it!
    end if


  end subroutine conv_ConvGridConv_0d

  !--------------------------------------------------------------------
  ! version for automatically doing convolution of two 2x2 arrays
  ! with array multiplication in the usual fortran sense
  subroutine conv_ConvGridConv_2dx2d(gc,gca, gcb)
    type(grid_conv),  intent(inout) :: gc(:,:)
    type(grid_conv),  intent(in)    :: gca(:,:), gcb(:,:)
    !----------------------------------------------------
    integer :: ir, im, ic, nr, nm, nc
    type(grid_conv) :: conv
    
    nr = assert_eq(size(gc,dim=1),size(gca,dim=1) ,'conv_ConvGridConv_2dx2d')
    nc = assert_eq(size(gc,dim=2),size(gcb,dim=2) ,'conv_ConvGridConv_2dx2d')
    nm = assert_eq(size(gca,dim=2),size(gcb,dim=1),'conv_ConvGridConv_2dx2d')
    
    ! allocate gc if need be, with default initialization of zero
    call InitGridConv(gca(1,1)%grid, gc) 
    !call SetToZero(gc)

    call AllocGridConv(gc(1,1)%grid,conv)
    do ic = 1, nc
       do ir = 1, nr
          do im = 1, nm
             call conv_ConvGridConv_0d(conv,gca(ir,im),gcb(im,ic))
             call AddWithCoeff(gc(ir,ic),conv)
          end do
       end do
    end do
    call Delete(conv)
    
  end subroutine conv_ConvGridConv_2dx2d
  
  !--------------------------------------------------------------------------
  ! Returns in gc (assumed preallocated) the commutator of
  ! gca and gcb. Does nothing clever about ordering of
  ! the arrays...
  !
  ! This is a very inefficient version because it does many more convolutions
  ! than are needed. Only for 2x2 arrays does it do a reasonably sensible job.
  !
  subroutine SetToCommutator_gc(gc,gca, gcb)
    type(grid_conv),  intent(inout) :: gc(:,:)
    type(grid_conv),  intent(in)    :: gca(:,:), gcb(:,:)
    !----------------------------------------------------
    type(grid_conv) :: cc(size(gc,dim=1),size(gc,dim=2))
    type(grid_conv) :: prod

    ! allocate gc if need be, with default initialization of zero
    call InitGridConv(gca(1,1)%grid, gc) 
    !call SetToZero(gc)

    ! allocate some more memory
    call AllocGridConv(gc(1,1)%grid,cc)

    !-- use a shortcut for 2-dim arrays; do things long-hand
    !   for others
    if (size(gc,dim=1) == 2 .and. size(gc,dim=2) == 2) then
    !-- for 2x2 arrays we will want to do something along the following lines
    ! (with a similar definition for B).
!> A:=matrix(2,2,[a11,a12,a21,a22]);
!                                                     [a11    a12]
!                                                A := [          ]
!                                                     [a21    a22]
!----------------------------------------------------------------------- 
!> matadd(multiply(A,B),-multiply(B,A));
!  [         a12 b21 - a21 b12           a11 b12 + a12 b22 - b11 a12 - b12 a22]
!  [                                                                          ]
!  [a21 b11 + a22 b21 - b21 a11 - b22 a21          a21 b12 - a12 b21          ]
       call AllocGridConv(gc(1,1)%grid,prod)
       
       ! res11 = a12 b21 - a21 b12
       call SetToConvolution(prod,gca(1,2),gcb(2,1))
       call AddWithCoeff(gc(1,1),prod)
       call SetToConvolution(prod,gca(2,1),gcb(1,2))
       call AddWithCoeff(gc(1,1),prod,-one)
       ! res22 = -res 11
       call AddWithCoeff(gc(2,2),gc(1,1),-one)
       ! res12 = a11 b12 + a12 b22 - b11 a12 - b12 a22
       call SetToConvolution(prod,gca(1,1),gcb(1,2))
       call AddWithCoeff(gc(1,2),prod)
       call SetToConvolution(prod,gca(1,2),gcb(2,2))
       call AddWithCoeff(gc(1,2),prod)
       call SetToConvolution(prod,gca(1,2),gcb(1,1))
       call AddWithCoeff(gc(1,2),prod,-one)
       call SetToConvolution(prod,gca(2,2),gcb(1,2))
       call AddWithCoeff(gc(1,2),prod,-one)
       ! res21 = a21 b11 + a22 b21 - b21 a11 - b22 a21
       call SetToConvolution(prod,gca(2,1),gcb(1,1))
       call AddWithCoeff(gc(2,1),prod)
       call SetToConvolution(prod,gca(2,2),gcb(2,1))
       call AddWithCoeff(gc(2,1),prod)
       call SetToConvolution(prod,gca(1,1),gcb(2,1))
       call AddWithCoeff(gc(2,1),prod,-one)
       call SetToConvolution(prod,gca(2,1),gcb(2,2))
       call AddWithCoeff(gc(2,1),prod,-one)
       
       call Delete(prod)
    else
       call SetToConvolution(cc,gca,gcb)
       call AddWithCoeff(gc,cc)
       call SetToConvolution(cc,gcb,gca)
       call AddWithCoeff(gc,cc,-one)
    end if
    
    call Delete(cc)
  end subroutine SetToCommutator_gc
  


  !--------------------------------------------------------------
  ! Carry out the convolution of gc on gq: but not always the
  ! most efficient way of dealing with the problem...
  function conv_ConvGridQuant_mat(gc,gq) result(gqout)
    type(grid_conv),  intent(in)    :: gc(:,:)
    real(dp),         intent(in)    :: gq(0:,:)
    real(dp)                        :: gqout(0:ubound(gq,dim=1),size(gc,dim=1))
    !---------------------------------------------
    integer :: ny, ic, ir, ncol, nrow

    ny = assert_eq(gc(1,1)%grid%ny,ubound(gq,dim=1),"conv_ConvGridQuant")
    ncol = assert_eq(size(gc,dim=2),size(gq,dim=2),"conv_ConvGridQuant")
    nrow = size(gc,dim=1)

    gqout = zero
    do ir = 1, nrow
       do ic = 1, ncol
          !gqout(:,ir) = gqout(:,ir) + gc(ir,ic) .conv. gq(:,ic)
          gqout(:,ir) = gqout(:,ir) +&
               & conv_ConvGridQuant_scalar(gc(ir,ic), gq(:,ic))
       end do
    end do
  end function conv_ConvGridQuant_mat

  

  !======================================================================
  ! things for easier access to a given grid point
  ! (but may not be as fast as just calling the routine)
  function conv_gdval_gdv(grid,val) result(gdv)
    type(grid_def), intent(in) :: grid
    real(dp),       intent(in) :: val
    type(gdval) :: gdv
    gdv%grid  = grid
    gdv%val = val
  end function conv_gdval_gdv
  function conv_gdval_vgd(val,grid) result(gdv)
    type(grid_def), intent(in) :: grid
    real(dp),       intent(in) :: val
    type(gdval) :: gdv
    gdv%grid  = grid
    gdv%val = val
  end function conv_gdval_vgd
  function conv_EvalGridQuant_atx(gq, gdv) result(res)
    real(dp),    intent(in) :: gq(:)
    type(gdval), intent(in) :: gdv
    real(dp) :: res
    res = EvalGridQuant(gdv%grid, gq, -log(gdv%val))
  end function conv_EvalGridQuant_atx
  function conv_EvalGridQuant_atx_1d(gq, gdv) result(res)
    real(dp),    intent(in) :: gq(:,:)
    type(gdval), intent(in) :: gdv
    real(dp) :: res(size(gq,dim=2))
    res = EvalGridQuant(gdv%grid, gq, -log(gdv%val))
  end function conv_EvalGridQuant_atx_1d
  function conv_EvalGridQuant_aty(gq, gdv) result(res)
    real(dp),    intent(in) :: gq(:)
    type(gdval), intent(in) :: gdv
    real(dp) :: res
    res = EvalGridQuant(gdv%grid, gq, gdv%val)
  end function conv_EvalGridQuant_aty
  function conv_EvalGridQuant_aty_1d(gq, gdv) result(res)
    real(dp),    intent(in) :: gq(:,:)
    type(gdval), intent(in) :: gdv
    real(dp) :: res(size(gq,dim=2))
    res = EvalGridQuant(gdv%grid, gq, gdv%val)
  end function conv_EvalGridQuant_aty_1d
  

  !======================================================================
  ! Routines for getting derived convolution operators
  !
  ! See also notes on CCN25-96
  !======================================================================
  ! Allocates and returns an array of probes, each of which must 
  ! be operated on with the convolver, after which the user should
  ! call SetDerivedConv, which formalises the results
  subroutine GetDerivedProbes(grid,probes)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: probes(:,:)
    !---------------------
    integer :: nprobes

    nprobes = 0
    call conv_GetNProbes(grid,nprobes)
    call AllocGridQuant(grid,probes,1,nprobes)
    call conv_SetProbes(grid,probes)

    ! safety measures -- this is a NASTY business, but cannot
    ! thing of another way to ensure & check that grid locking is
    ! off...
    override_grid_locking = .true.
    nconv_with_override_off = 0
  end subroutine GetDerivedProbes

  !----------------------------------------------------------------------
  ! Sets up the "derived" convolution operator.
  !
  ! This routine actually just does the housekeeping -- the actual "hard"
  ! work is done by conv_SetDerivedConv_rec
  !
  ! NB: gc must previously have been allocated
  subroutine SetDerivedConv(gc,probes)
    real(dp),        pointer       :: probes(:,:)
    type(grid_conv), intent(inout) :: gc

    call SetDerivedConv_nodealloc(gc,probes)
    deallocate(probes)
  end subroutine SetDerivedConv

  !----------------------------------------------------------------------
  ! runs through each of the elements of grid and establishes
  ! the number of probes needed. Currently only supports cases
  ! in which there is just a single probe...
  recursive subroutine conv_GetNProbes(grid,nprobes)
    use warnings_and_errors
    type(grid_def), intent(in)  :: grid
    integer,        intent(out) :: nprobes
    !---------
    integer :: nprobes_tmp, isub
    if (grid%nsub /= 0) then
       nprobes = 0
       do isub = 1, grid%nsub
          call conv_GetNProbes(grid%subgd(isub),nprobes_tmp)
          nprobes = max(nprobes, nprobes_tmp)
       end do
    else
       select case(grid%order)
       case(:LIN_ORDER-1)
          nprobes = 1
       case(LIN_ORDER)
          nprobes = 2
       case(LIN_ORDER+1:)
          nprobes = grid%order+2
       end select
    end if
  end subroutine conv_GetNProbes

  !----------------------------------------------------------------------
  ! Sets the exact form of the probes... Only works for order<0 for now.
  recursive subroutine conv_SetProbes(grid,probes)
    type(grid_def), intent(in)  :: grid
    real(dp),       intent(out) :: probes(0:,:)
    integer :: isub, iprobe
    if (grid%nsub /= 0) then
       do isub = 1, grid%nsub
          call conv_SetProbes(grid%subgd(isub),&
               &probes(grid%subiy(isub):grid%subiy(isub+1)-1,:))
       end do
    else
       select case(grid%order)
       case(:LIN_ORDER-1)
          probes(0,1) = one
          probes(1:,:) = zero
       case(LIN_ORDER)
          probes(:,:) = zero
          probes(0,1) = one
          probes(1,2) = one
       case(LIN_ORDER+1:)
          probes = zero
          do iprobe = 1, grid%order+2
             probes(iprobe-1,iprobe) = one
          end do
       end select
    end if
  end subroutine conv_SetProbes
  
  
  !----------------------------------------------------------------------
  ! Does the work in the setup of the "derived" convolution operator.
  recursive subroutine SetDerivedConv_nodealloc(gc,probes)
    use warnings_and_errors
    real(dp),        intent(in)    :: probes(0:,:)
    type(grid_conv), intent(inout) :: gc
    integer :: isub, order, iprobe

    override_grid_locking = .false.
    if (nconv_with_override_off /= 0) call wae_error(&
         &'SetDerivedConv_nodealloc',&
         &'Detected convolutions while lock overried off')

    if (gc%grid%nsub /= 0) then
       do isub = 1, gc%grid%nsub
          call SetDerivedConv_nodealloc(gc%subgc(isub),&
               &probes(gc%grid%subiy(isub):gc%grid%subiy(isub+1)-1,:))
       end do
    else
       select case(gc%grid%order)
       case(:LIN_ORDER-1)
          gc%conv = probes
       case(LIN_ORDER)
          gc%conv(0:gc%grid%ny-1,FULL) = probes(1:,2)
          ! fake value here -- but since we will actually only be 
          ! interested in lower for this point, it does not
          ! matter!
          gc%conv(gc%grid%ny,FULL) = zero
          ! actually it is lower that has been calculated here
          gc%conv(0:,UPPR) = probes(0:,1)
          ! now convert it to upper
          gc%conv(0:,UPPR) = gc%conv(0:,FULL) - gc%conv(0:,UPPR)
       case(LIN_ORDER+1:)
          order = gc%grid%order
          ! a safety check -- if I forget anything it should
          ! jump out...
          gc%conv = 1e90_dp ! dummy value to detect uninitialized pieces...
          do iprobe = 1, order+1
             gc%conv(:,iprobe) = probes(:,order+2-iprobe)
          end do
          gc%conv(0:gc%grid%ny-order-1,0) = probes(order+1:,order+2)
          ! originally known, but unknowable with this method (and
          ! irrelevant).
          gc%conv(gc%grid%ny-order:,0) = zero
       end select
    end if
  end subroutine SetDerivedConv_nodealloc
  
end module convolution
