! Program which checks different ways of initializing pdfs in a grid



program test_hoppet

! Use the main module
use hoppet_v1
! Use our own modules
use jrc_modules

implicit none

type(grid_def) :: grid
real(dp),pointer :: gluon(:)
real(dp),pointer :: gq(:,:)
real(dp) gluon_at_y,y,Q
real(dp) pdf_at_y(-6:6)

! Initialize a grid object
! order -> Interpolation order associated to the grid object
call InitGridDef(grid,dy=0.2_dp,ymax=8.0_dp,order=3)


! Initialize the gluon grid with our subroutine
! Allocate array
call AllocGridQuant(grid,gluon)

! Initialize gluon with a function
call InitGridQuant(grid,gluon,toy_gluon)

y = 4_dp
gluon_at_y = EvalGridQuant(grid,gluon,y)
print *,y,exp(-y),exp(-y)*gluon_at_y

! Initialize gluon with a subroutine
call InitGridQuantSub(grid,gq,LHAsub2)

!y = 4_dp
!gluon_at_y = EvalGridQuant(grid,gluon,y)
!print *,y,exp(-y),exp(-y)*gluon_at_y

! Initialize with LHAPDF
Q=2
call AllocPDF(grid,gq)
call InitGridQuantLHAPDF(grid, gq, LHAsub, Q)
pdf_at_y = EvalGridQuant(grid,gq(:,-6:6),y)
print *,y,exp(-y),exp(-y)*pdf_at_y(1),exp(-y)*pdf_at_y(0)








end program test_hoppet

!-----------------------------------------------------------


