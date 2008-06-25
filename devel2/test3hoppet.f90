

module pdf_initial_condition
  use hoppet_v1
  implicit none
contains
  function unpolarized_dummy_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),-6:7)

    ! clean method for labelling as PDF as being in the human representation
    ! call LabelPdfAsHuman(pdf)
    ! by setting everything to zero (notably pdf(:,7)) representation
    ! is automatically set to be human
      pdf(:,:) = 0
   

    ! iflv_g is pre-defined integer parameter (=0) for symbolic ref. to gluon
    pdf(:,iflv_g) = 1.7_dp * xvals**(-0.1_dp) * (1-xvals)**5
    pdf(:,iflv_u) = 1.7_dp * xvals**(-0.1_dp) * (1-xvals)**5


  end function unpolarized_dummy_pdf
end module pdf_initial_condition


program test_hoppet

! Use the main module
use hoppet_v1
! Use our own modules
use jrc_modules
use pdf_initial_condition
implicit none 

type(grid_def) :: grid
type(pdf_rep) :: prep
real(dp), pointer :: pdf(:,:),pdf_evln(:,:)
real(dp) y, gluon_at_y
integer nf_lcl

! Initialize a grid object
! order -> Interpolation order associated to the grid object
call InitGridDef(grid,dy=0.2_dp,ymax=8.0_dp,order=3)


call AllocPDF(grid,pdf)
pdf = unpolarized_dummy_pdf(xValues(grid))

y=6_dp
gluon_at_y = EvalGridQuant(grid,pdf(:,iflv_g),y)
print *,y,gluon_at_y

call AllocPDF(grid,pdf_evln)

nf_lcl=5
prep%nf = nf_lcl


print *,"Converting pdf representations"
print *,GetPdfRep(pdf)
call CopyHumanPdfToEvln(prep,pdf,pdf_evln)
print *,GetPdfRep(pdf_evln)
call CopyEvlnPdfToHuman(prep,pdf_evln,pdf)
print *,GetPdfRep(pdf)




end program test_hoppet

!-----------------------------------------------------------


