! Check of the various functionalities of the hoppet toolkit
! i) Fast parton evolution
! ii) PDF manipulation
! iii) Convolutions

!
! ToDo
! Understand DGLAP evolution
! Solve problem with F77 program
!


program test_hoppet

! Use the main module
use hoppet_v1
! Use our own modules
use jrc_modules
! Communicator of terms in convolutions
use convolution_communicator

implicit none 

integer i,sizevec,k,ntimes
! Define a grid object
type(grid_def) :: grid
real(dp),pointer :: gluon(:),yvals(:),Pgg_x_gluon(:),pdf_ns(:)
real(dp),pointer :: evolved_gluon(:),evolved_pdf_ns(:)
real(dp),pointer :: unit_grid(:)
real(dp) :: x,g,y,gluon_at_y
! Define a grid convolution operator
type(grid_conv) :: Pgg,GammaNS
real(dp) :: eps,Q,factor
type(running_coupling) :: coupling
real(dp) gammax,gammax_nsp,q2i,q2f,t2,t1,pdf_in,pdf_out

type(split_mat) ::PA,PB,PC

write(6,*) " "
write(6,*) "Checks of the HOPPET parton toolkit"
write(6,*) " "

! Initialize a grid object
! order -> Interpolation order associated to the grid object
call InitGridDef(grid,dy=0.2_dp,ymax=8.0_dp,order=3)

factor=2_dp
! Allocate space
call AllocSplitMat(grid,PB,5)
call SetToZero(PB)
!call InitSplitMat(PA,PB)




call exit(-10)

! Initialize the gluon grid with our subroutine
! Allocate array
call AllocGridQuant(grid,gluon)
! call InitGridQuant(grid,gluon,example_gluon)
call InitGridQuant(grid,gluon,toy_gluon)
y = 4_dp
gluon_at_y = EvalGridQuant(grid,gluon,y)
print *,y,exp(-y),gluon_at_y

call AllocGridQuant(grid,yvals)
yvals = yValues(grid)
!print *,xvals(1),xvals(2)
!print *,size(xvals)

sizevec = size(yvals)-1
do i=1,sizevec
   y=yvals(i)
   call example_gluon(y,g)
   gluon(i)=g
!   print *,i,x,g
enddo
! deallocate array
deallocate(yvals)

! Access the gluon grid
y=4.0_dp
gluon_at_y = EvalGridQuant(grid,gluon,y)
print *,y,exp(-y),gluon_at_y

print *,"Checking running coupling"
! Initialize running coupling
call InitRunningCoupling(coupling)
! Evaluate running coupling
Q=91.2
print *,Q,Value(coupling,Q)
Q=2
print *,Q,Value(coupling,Q)

! Check convolutions of functions
! which are distributions
! Checking convolutions
print *,"Check convolutions"
! Function to convolute 
! Must select through convolution communicator which part of
! it must be returned
cc_piece = cc_REAL
print *,y,exp(-y),Pgg_func(y)
y=0.1_dp
print *,y,exp(-y),Pgg_func(y)
! Initialize unit grid for checks of convolution
call AllocGridQuant(grid,unit_grid)
! call InitGridQuant(grid,gluon,example_gluon)
call InitGridQuant(grid,unit_grid,unit)
! Grid convolutions
call AllocGridConv(grid,Pgg)
call InitGridConv(grid,Pgg,Pgg_func)
! Use the convolution
call AllocGridQuant(grid,Pgg_x_gluon)
print *,"Convolution initialized"
! Standard notation
! Pgg_x_gluon = Pgg.conv.gluon
! Alternative notation for convolution operator
Pgg_x_gluon = Pgg*unit_grid
y=0.6931_dp
print *,y,exp(-y),EvalGridQuant(grid,Pgg_x_gluon,y)
y=4.0_dp
print *,y,exp(-y),EvalGridQuant(grid,Pgg_x_gluon,y)
! Compare with expected analytical result
! Check dgauss
! print *,y,exp(-y),Pgg_func_int(y)

! Check convolutions with a pdf
call AllocGridQuant(grid,gluon)
call InitGridQuant(grid,gluon,toy_gluon)
!call InitGridQuantSub(grid,gluon,example_gluon)
! Plot toy gluon used for convolutions
! Save the NS evolution factor
do k=1,100
   x=0.001_dp*(1000_dp)**(dfloat(k-1)/100)
   y=log(1_dp/x)
   !write(6,*) x,x*toy_gluon(y),x*EvalGridQuant(grid,gluon,y)
enddo
call AllocGridQuant(grid,Pgg_x_gluon)
Pgg_x_gluon = Pgg.conv.gluon
y=0.693_dp
print *,"Check convolution with a gluon pdf"
! Set number of calls to zero
number_of_calls = 0
print *,y,exp(-y),EvalGridQuant(grid,Pgg_x_gluon,y)
print *,"Number of calls = ",number_of_calls
y=4.0_dp
print *,"Check convolution with a gluon pdf"
print *,y,exp(-y),EvalGridQuant(grid,Pgg_x_gluon,y)


! Check multiflavour grids
! If we want hoppet to deal with automatic allocation of the
! pdfs, they should be pointer arrays


end program test_hoppet

!-----------------------------------------------------------


