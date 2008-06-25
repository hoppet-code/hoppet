!-------------------------------------------------------
! 
!    The hoppet tools for performing and precomputing
!    fast integrations and interpolations are used
!    with the NNPDF x-space evolution kernels
!
!    JRC, 11/07
!
!-------------------------------------------------

! ToDo:
! Use correct evolution kernel for the u_v distribution
! Compare with a version of checkkin.f

program nnpdf_hoppet

! Use the main module
use hoppet_v1
! Use our own modules
use jrc_modules
! Communicator of terms in convolutions
use convolution_communicator

implicit none 

integer i,k,ntimes
! Define a grid object
type(grid_def) :: grid,subgrids(2)
real(dp),pointer :: pdf_ns(:),evolved_pdf_ns(:)
real(dp) :: x,g,y
! Define a grid convolution operator
type(grid_conv) :: Pgg,GammaNS
real(dp) :: eps,Q
type(running_coupling) :: coupling
real(dp) :: gammax,gammax_nsp,q2i,q2f,t2,t1,pdf_in,pdf_out,pdf_out_exact
real(dp) :: rel_diff


! Initialize a grid object
! order -> Interpolation order associated to the grid object
call InitGridDef(subgrids(1),dy=0.20_dp,ymax=7.0_dp,order=3)
call InitGridDef(subgrids(2),dy=0.05_dp,ymax=2.0_dp,order=3)
! Lock all subgrids together on a single grid
call InitGridDef(grid,subgrids,locked=.true.)


! Check the NNPDF exact evolution with the Hoppet tools
print *,"                     "
print *,"Check NNPDF evolution"

! Check our NS evolution kernel
! First all constants and variables must be initializated
q2i = 2_dp
q2f = 10000_dp

! Save the NS evolution factor
do k=1,100
   x=0.001_dp*(1000_dp)**(dfloat(k-1)/100)
   gammax = gammax_nsp(x,q2i,q2f)
   !write(6,*) x,gammax
enddo

! Check dgauss
y=4.0_dp
print *,y,exp(-y),Pgg_func_int(y)

! Check the convolution of a precomputed evolution
! kernel
! (Accuracy and number of calls to input pdf for various grids)

! Grid convolutions
call AllocGridQuant(grid,pdf_ns)
call InitGridQuant(grid,pdf_ns,input_pdf_ns)
call AllocGridConv(grid,GammaNS)
print *,"Initialize the grid for the convolutions"
call InitGridConv(grid,GammaNS,GammaNS_func)
print *,"Grid for convolution initialized!"
! Use the convolution
call AllocGridQuant(grid,evolved_pdf_ns)
evolved_pdf_ns = GammaNS*pdf_ns

! Save evolved pdf
! Compare with exact result for the parton evolution
do k=1,50
   x=0.001_dp*(1000_dp)**(dfloat(k-1)/50)
   y=log(1_dp/x)
   pdf_in = x*EvalGridQuant(grid,pdf_ns,y)
   pdf_out = x*EvalGridQuant(grid,evolved_pdf_ns,y)
   pdf_out_exact = x*evolved_uv(y)
   rel_diff=abs( (pdf_out -pdf_out_exact)/pdf_out_exact )
   write(6,'(5f14.9)') x,pdf_in,pdf_out,pdf_out_exact,rel_diff
enddo

!call exit(-10)

call CPU_time(t1)
ntimes=100000
! Plot evolved NS pdf
do i=1,ntimes
   ! Each time different initial condition
   call InitGridQuant(grid,pdf_ns,input_pdf_ns)
   ! Do the convolution
   evolved_pdf_ns = GammaNS*pdf_ns
   ! Compute the output of evolution
   do k=1,100
      x=0.001_dp*(1000_dp)**(dfloat(k-1)/100)
      y=log(1_dp/x)
      pdf_out = x*EvalGridQuant(grid,evolved_pdf_ns,y)
   enddo
enddo
call cpu_time(t2)
print *,"CPU time = ",(t2-t1)/ntimes/100

! Cross check with exact evolution in our case

end program nnpdf_hoppet

!---------------------------------------------------


