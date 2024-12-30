!! Test of array bounds bug that is present in a family of compilers that, as of 2024-12
!! includes AOCC (5.0) nvfortran (24.11). See the bug report at:
!! https://forums.developer.nvidia.com/t/bugreport-for-nvfortran-treatment-of-array-bounds-in-the-specification-part/278324
program array_bound_test
  implicit none
  real :: arr_in1(0:10, 2:5), arr1(0:10,2:5)
  integer  i, j
  do i=0,10 
  do j=2,5
    arr_in1(i,j) = 1000*i+j
  end do
  end do

  call evaluate(arr_in1,arr1)
contains

  function justcopy(arr_in)
   implicit none
   real, intent(in) :: arr_in(0:,:)
   real :: justcopy(0:ubound(arr_in,dim=1), 2:5)  !This works for "A" GNU,Intel/InteLLVM,lfortran,NAG,IBM. Does not work for "B" AOCC, NVidia, ARM.
   !real :: justcopy(0:size(arr_in,dim=1)-1, 2:5) !This works for all compilers
   write(*,*) "arr_in, dim 1 (lb,ub,sz):", lbound(arr_in,dim=1),ubound(arr_in,dim=1),size(arr_in,1)
   write(*,*) "arr_in, dim 2 (lb,ub,sz):", lbound(arr_in,dim=2),ubound(arr_in,dim=2),size(arr_in,2)
   write(*,*) "justcp, dim 1 (lb,ub,sz):", lbound(justcopy,dim=1),ubound(justcopy,dim=1),size(justcopy,1)
   write(*,*) "justcp, dim 2 (lb,ub,sz):", lbound(justcopy,dim=2),ubound(justcopy,dim=2),size(justcopy,2)
   justcopy=arr_in
  end function justcopy
  
  subroutine evaluate(arr_in1,arr_out1)
    implicit none
    real,    intent(in)  :: arr_in1(:,:)
    real,    intent(out) :: arr_out1(:,:)
    real, allocatable  :: X(:,:)

    !arr_out1 = justcopy(arr_in1)
    X = justcopy(arr_in1)       
    write(*,*) "X     , dim 1 (lb,ub,sz):", lbound(X,dim=1),ubound(X,dim=1),size(X,1)
    write(*,*) "X     , dim 2 (lb,ub,sz):", lbound(X,dim=2),ubound(X,dim=2),size(X,2)

    if (size(X,1) /= size(arr_in1,1)) then
       print "(A,I2,A,I2)","Error: size of X(dim1)=",size(X,1),&
                           " was not equal to size of arr_in1(dim1)=",size(arr_in1,1)
      error stop 1
    end if
    
    arr_out1 = X
  end subroutine evaluate  
end program array_bound_test
