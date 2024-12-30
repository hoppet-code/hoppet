program aocc_run_test
  implicit none
  real :: arr_in1(0:10, 2:5), arr1(0:10,2:5)
  real :: arr_in2(0:10, 2:5), arr2(0:10,2:5)
  real :: arr_in3(0:10, 2:5), arr3(0:10,2:5)
  integer  i, j
  do i=0,10 
  do j=2,5
    arr_in1(i,j) = 1000*i+j
    arr_in2(i,j) = 1000*i+j
    arr_in3(i,j) = 1000*i+j
  end do
  end do

  call evaluate(arr_in1,arr_in2,arr_in3,arr1,arr2,arr3)
contains

  function justcopy(arr_in) 
   real, intent(in) :: arr_in(0:,:)
   real :: justcopy(0:ubound(arr_in,dim=1), 2:5)  !This works for "A" GNU,Intel/InteLLVM,lfortran,NAG,IBM. Does not work for "B" AOCC, NVidia, ARM.
   !real :: justcopy(0:size(arr_in,dim=1)-1, 2:5) !This works for all compilers
   write(*,*) "arr_in, dim 1 (lb,ub,sz):", lbound(arr_in,dim=1),ubound(arr_in,dim=1),size(arr_in,1)
   write(*,*) "arr_in, dim 2 (lb,ub,sz):", lbound(arr_in,dim=2),ubound(arr_in,dim=2),size(arr_in,2)
   write(*,*) "justcp, dim 1 (lb,ub,sz):", lbound(justcopy,dim=1),ubound(justcopy,dim=1),size(justcopy,1)
   write(*,*) "justcp, dim 2 (lb,ub,sz):", lbound(justcopy,dim=2),ubound(justcopy,dim=2),size(justcopy,2)
   justcopy=arr_in
   !write(*,*)lbound(justcopy,dim=1),ubound(justcopy,dim=1),size(justcopy,1)
   !write(*,*) "justcp, dim 1 (lb,ub,sz):", lbound(justcopy,dim=1),ubound(justcopy,dim=1),size(justcopy,1)
  end function justcopy
  
  subroutine evaluate(arr_in1,arr_in2,arr_in3,arr_out1,arr_out2,arr_out3)
    real,    intent(in)  :: arr_in1(:,:)
    real,    intent(in)  :: arr_in2(:,:)
    real,    intent(in)  :: arr_in3(:,:)
    real,    intent(out) :: arr_out1(:,:)
    real,    intent(out) :: arr_out2(:,:)
    real,    intent(out) :: arr_out3(:,:)
    real, allocatable  :: X(:,:)

    !arr_out1 = justcopy(arr_in1)
    X = justcopy(arr_in2)       
    write(*,*) "X     , dim 1 (lb,ub,sz):", lbound(X,dim=1),ubound(X,dim=1),size(X,1)
    write(*,*) "X     , dim 2 (lb,ub,sz):", lbound(X,dim=2),ubound(X,dim=2),size(X,2)

    if (size(X,1) /= 11) error stop 1
    
    !arr_out2 = X
    !arr_out3 = 1.0*justcopy(arr_in3)

    !write(*,*)"X--       ", lbound(X), ubound(X),size(X,1) !This line is different for "A" and "B" compilers with "ubound"
    !write(*,*)X
    ! write(*,*)"arr_out1--", lbound(arr_out1), ubound(arr_out1),size(arr_out1,1)
    ! !write(*,*)arr_out1
    ! write(*,*)"arr_out2--", lbound(arr_out2), ubound(arr_out2),size(arr_out2,1)
    ! !write(*,*)arr_out2
    ! write(*,*)"arr_out3--", lbound(arr_out3), ubound(arr_out3),size(arr_out3,1)
    ! !write(*,*)arr_out3
    ! !With ubound arr_out3=arr_out2 != arr_out1 != X
  end subroutine evaluate  
end program aocc_run_test
