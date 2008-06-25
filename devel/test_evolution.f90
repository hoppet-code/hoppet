program test_hoppet
! Use the main module
use hoppet_v1
use jrc_modules
implicit none 

type(grid_def) :: grid
type(running_coupling) :: coupling
real(dp) :: pole_mass
type(dglap_holder) :: dh
real(dp) Q_init,Q_end,pdf_at_y(-6:6),y
real(dp),pointer :: PDF(:,:)
integer nloop

write(6,*) " "
write(6,*) "Checks of the HOPPET parton toolkit"
write(6,*) " "

! Initialize a grid object
! order -> Interpolation order associated to the grid object
call InitGridDef(grid,dy=0.2_dp,ymax=8.0_dp,order=3)

! Set running coupling object
call InitRunningCoupling(coupling)

! Check pole masses
pole_mass = QuarkMass(coupling,5)
print *,5,pole_mass

! Init DGLAP holder
nloop=2
call InitDglapHolder(grid,dh,1,nloop)

! Allocate memory to a multi-flavoured grid pdf
call AllocPDF(grid,PDF)   ! Allocated PDF(0:-6,7)
Q_init=1.414_dp
! Initalize the pdf grid using a LH-like routine
call InitPDF_LHAPDF(grid,PDF,LHASUB,Q_init)

y=6_dp
pdf_at_y = EvalGridQuant(grid,PDF(:,-6:6),y)
print *,y,pdf_at_y(0)
! Evolve pdfs
Q_end=100
call EvolvePDF(dh,pdf,coupling,Q_init,Q_end)
pdf_at_y = EvalGridQuant(grid,PDF(:,-6:6),y)
print *,y,pdf_at_y(0)



end program test_hoppet

!-----------------------------------------------------------


