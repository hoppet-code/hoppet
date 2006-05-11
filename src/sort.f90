!! A module that provides access to a quicksort form of sorting
!! routine. 
!!
!! Not clear that it is very efficient (e.g. compared to C++ stdlib
!! sort). It is based on the ideas behind NR indexx, but has been
!! rewritten (without looking at the NR f90 code) to avoid copyright
!! issues. Considerable tests have been performed [but no guarantees
!! are provided...!]
!!
module sort
  use types
  use warnings_and_errors
  use assertions
  implicit none
  private

  

  interface indexx
     module procedure indexx_dp, indexx_int
  end interface
  public :: indexx

  interface swap
     module procedure swap_int, swap_dp
  end interface
  

contains

  !-------------------------------------------------------------------------
  !! A non-recursive version of the quick sort algorithm. It returns
  !! an array idx such that array(idx) is sorted.
  !!
  !! Following NR outline, it uses insertion for for sorting small
  !! sub-segments, sentinels to reduce the number of checks when
  !! pivoting and median of first, middle and last elements in order
  !! to select the pivot, to reduce the chances of the worst case
  !! (N^2) occurring.
  !!
  !! Note it uses a fixed-size stack, which is not so great an idea...
  !! (recursion would have had the advantage that this would not have
  !! been a problem).
  subroutine indexx_dp(array,idx)
    real(dp), intent(in)  :: array(:)
    integer,  intent(out) :: idx(size(array))
    !----------------------------------
    integer, parameter :: insertion_threshold = 7, stack_size = 50
    type lrpair
       integer :: l, r
    end type lrpair
    type(lrpair) :: stack(stack_size)
    integer  :: nmid, i, j, l, r, ixstack, medidx, n
    real(dp) :: medval
    
    !n = size(array)
    !if (n/=size(idx)) stop 'array and idx are different sizes'
    n = assert_eq(size(idx),size(array),'indexx_dp')

    forall(i=1:n) idx(i) = i
    
    ! push the whole array section onto the stack
    ixstack = 1
    stack(1) = lrpair(1,n); ixstack = 1
    do 
       ! pull an array section off the stack
       l = stack(ixstack)%l
       r = stack(ixstack)%r
       ixstack = ixstack - 1
       ! get to work on it
       if (r-l < insertion_threshold) then
          ! run an inlined insertion sort.
          do i = l+1, r
             medidx = idx(i)
             medval = array(medidx)
             ! look for the place where we should insert the element "medval", 
             ! and as we go about it, shift the other elements up by one
             do j = i-1, l, -1
                if (medval < array(idx(j))) then
                   idx(j+1) = idx(j)
                else
                   exit
                end if
             end do
             ! insert medval in its correct location
             idx(j+1) = medidx
          end do
       else
          ! NR idea is to take the median value of first, middle and
          ! last elements, so do a little insertion sort to get these
          ! three elements into order
          nmid = (l+r)/2
          ! first place the middle element in position l+1
          call swap(idx(l+1),idx(nmid))
          if (array(idx(l+1)) < array(idx(l))  ) call swap(idx(l+1),idx(l))
          if (array(idx(r))   < array(idx(l+1))) then
             call swap(idx(r),idx(l+1))
             if (array(idx(l+1)) < array(idx(l))) call swap(idx(l+1),idx(l))
          end if
          medidx = idx(l+1)
          medval = array(medidx)

          ! Then run up from first (non-trivial) element, and down from last
          ! element, swapping things as need be until i and j cross
          i = l+1; j = r
          do
             do 
                i = i + 1
                if (array(idx(i)) >= medval) exit
             end do
             do 
                j = j - 1
                if (array(idx(j)) <= medval) exit
             end do
             if (i > j) exit
             call swap(idx(i), idx(j))
          end do
          ! now insert medval (the partitioning element) into the
          ! correct position
          idx(l+1) = idx(j)
          idx(j)   = medidx

          ! now put things onto the stack
          if (ixstack+2 > stack_size) stop 'stack too small'
          stack(ixstack+1) = lrpair(l,j-1)
          stack(ixstack+2) = lrpair(j,r)
          ixstack = ixstack + 2
       end if
       
       if (ixstack == 0) exit
    end do
    
  end subroutine indexx_dp



  !-------------------------------------------------------------------------
  !! A non-recursive version of the quick sort algorithm. It returns
  !! an array idx such that array(idx) is sorted.
  !!
  !! Following NR outline, it uses insertion for for sorting small
  !! sub-segments, sentinels to reduce the number of checks when
  !! pivoting and median of first, middle and last elements in order
  !! to select the pivot, to reduce the chances of the worst case
  !! (N^2) occurring.
  !!
  !! Note it uses a fixed-size stack, which is not so great an idea...
  !! (recursion would have had the advantage that this would not have
  !! been a problem).
  subroutine indexx_int(array,idx)
    integer,  intent(in)  :: array(:)
    integer,  intent(out) :: idx(size(array))
    !----------------------------------
    integer, parameter :: insertion_threshold = 7, stack_size = 50
    type lrpair
       integer :: l, r
    end type lrpair
    type(lrpair) :: stack(stack_size)
    integer :: nmid, i, j, l, r, ixstack, medidx, n
    real    :: medval
    
    !n = size(array)
    !if (n/=size(idx)) stop 'array and idx are different sizes'
    n = assert_eq(size(idx),size(array),'indexx_dp')

    forall(i=1:n) idx(i) = i
    
    ! push the whole array section onto the stack
    ixstack = 1
    stack(1) = lrpair(1,n); ixstack = 1
    do 
       ! pull an array section off the stack
       l = stack(ixstack)%l
       r = stack(ixstack)%r
       ixstack = ixstack - 1
       ! get to work on it
       if (r-l < insertion_threshold) then
          ! run an inlined insertion sort.
          do i = l+1, r
             medidx = idx(i)
             medval = array(medidx)
             ! look for the place where we should insert the element "medval", 
             ! and as we go about it, shift the other elements up by one
             do j = i-1, l, -1
                if (medval < array(idx(j))) then
                   idx(j+1) = idx(j)
                else
                   exit
                end if
             end do
             ! insert medval in its correct location
             idx(j+1) = medidx
          end do
       else
          ! NR idea is to take the median value of first, middle and
          ! last elements, so do a little insertion sort to get these
          ! three elements into order
          nmid = (l+r)/2
          ! first place the middle element in position l+1
          call swap(idx(l+1),idx(nmid))
          if (array(idx(l+1)) < array(idx(l))  ) call swap(idx(l+1),idx(l))
          if (array(idx(r))   < array(idx(l+1))) then
             call swap(idx(r),idx(l+1))
             if (array(idx(l+1)) < array(idx(l))) call swap(idx(l+1),idx(l))
          end if
          medidx = idx(l+1)
          medval = array(medidx)

          ! Then run up from first (non-trivial) element, and down from last
          ! element, swapping things as need be until i and j cross
          i = l+1; j = r
          do
             do 
                i = i + 1
                if (array(idx(i)) >= medval) exit
             end do
             do 
                j = j - 1
                if (array(idx(j)) <= medval) exit
             end do
             if (i > j) exit
             call swap(idx(i), idx(j))
          end do
          ! now insert medval (the partitioning element) into the
          ! correct position
          idx(l+1) = idx(j)
          idx(j)   = medidx

          ! now put things onto the stack
          if (ixstack+2 > stack_size) stop 'stack too small'
          stack(ixstack+1) = lrpair(l,j-1)
          stack(ixstack+2) = lrpair(j,r)
          ixstack = ixstack + 2
       end if
       
       if (ixstack == 0) exit
    end do
    
  end subroutine indexx_int


  subroutine swap_dp(a,b)
    real, intent(inout) :: a,b
    real                :: dummy
    dummy=a
    a=b
    b=dummy
  end subroutine swap_dp

  subroutine swap_int(a,b)
    integer, intent(inout) :: a,b
    integer                :: dummy
    dummy=a
    a=b
    b=dummy
  end subroutine swap_int

end module sort
