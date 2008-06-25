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

implicit none 

real(dp),pointer :: PDF(:,:),PDF_evln(:,:)
type(grid_def) :: grid
real(dp) :: Q,y,pdf_at_y(-6:6),x,pdf_evln_at_y(-6:6)
real(dp) :: pdf_nnpdf(1:13),as2pi,dt,gluon_at_y
integer :: nf_lcl
! Splitting functions
type(split_mat) :: P_sum,P_LO,P_NLO
real(dp), pointer :: dq(:,:)
type(running_coupling) :: coupling

write(6,*) " "
write(6,*) "Checks of the HOPPET parton toolkit"
write(6,*) " "

! Initialize a grid object
call InitGridDef(grid,dy=0.5_dp,ymax=10.0_dp,order=2)

! Allocate memory to a multi-flavoured grid pdf
call AllocPDF(grid,PDF)   ! Allocated PDF(0:-6,7)

Q=1.414_dp
! Initalize the pdf grid using a LH-like routine
call InitPDF_LHAPDF(grid,PDF,LHASUB,Q)

! Extract all flavours simulmatneously
y=3_dp
! The option below (doc) does not work!
! pdf_at_y = 2 * EvalGridQuant(grid,PDF(:,-6:7),y)
pdf_at_y = EvalGridQuant(grid,PDF(:,-6:6),y)
print *,y,pdf_at_y(0)
pdf_at_y = PDF(:,-6:6).aty.(y.with.grid)
print *,y,pdf_at_y(0)
print *,y,pdf_at_y(iflv_g)

gluon_at_y = EvalGridQuant(grid,PDF(:,iflv_g),y)
print *,y,gluon_at_y
gluon_at_y = PDF(:,iflv_g).aty.(y.with.grid)
print *,y,gluon_at_y

x=exp(-y)
write(6,'(5f15.12)'),x,x*pdf_at_y(0),x*pdf_at_y(1),&
     x*pdf_at_y(2),x*pdf_at_y(-1)
! Compare with exact result
call LHASUB(x,Q,pdf_at_y)
write(6,'(5f15.12)'),x,x*pdf_at_y(0),x*pdf_at_y(1),&
     x*pdf_at_y(2),x*pdf_at_y(-1)

! Copy the pdf array from the 'human' to the 'evln' representations
call AllocPDF(grid,pdf_evln)
! Set number of active flavours
nf_lcl=4
call CopyHumanPdfToEvln(nf_lcl,pdf,pdf_evln)
! Use now the evln hoppet pdf basis convention
! ipdf = 1 -> Singlet
! ipdf = 0 -> gluon
! ipdf = -1 -> uV+dV+sV+cV
pdf_evln_at_y = EvalGridQuant(grid,PDF_evln(:,:),y)
x=exp(-y)
write(6,'(5(f15.12,1x))'),x,x*pdf_evln_at_y(0),x*pdf_evln_at_y(1),&
     x*pdf_evln_at_y(-1)
! Compare with exact result
! Use now the NNPdf evolution pdf basis convention
! ipdf = 1 -> Singlet
! ipdf = 2 -> gluon
! ipdf = 3 -> uV
! ipdf = 4 -> dV
call PDFINLH(x,pdf_nnpdf)
write(6,'(5(f15.12,1x))'),x,x*pdf_nnpdf(2),x*pdf_nnpdf(1),&
     x*pdf_nnpdf(3)+x*pdf_nnpdf(4)

! Set value of the number of flavours
nf_lcl=4
call qcd_SetNf(nf_lcl)

! Play a bit with splitting functions
call InitSplitMatLO(grid,P_LO)
call InitSplitMatNLO(grid,P_NLO)
call InitSplitMat(P_sum,P_LO)

! Set value of the coupling constant
call InitRunningCoupling(coupling)
Q=10_dp
as2pi=Value(coupling,Q)

call AddWithCoeff(P_sum,P_NLO,as2pi)

! One step for the DGLAP evolution
dt=0.1_dp
call AllocPDF(grid,dq)
dq = ( as2pi * dt ) * ( P_sum.conv.pdf_evln )
print *,dq(1,0),dq(10,1)

! Delete the dynamically allocated object
call Delete(P_sum)


end program test_hoppet

!-----------------------------------------------------------


