program test_qed_lhapdf
  use hoppet; use convolution
  use io_utils
  implicit none
  real(dp) :: ymax = 20.0_dp, dy = 0.10_dp
  type(grid_def) :: grid
  character(len=200) :: pdfname
  real(dp), pointer :: pdf(:,:)
  real(dp)          :: moments(-6:8), Q
  integer           :: imem
  
  call InitGridDefDefault(grid, dy, ymax)
  call AllocGridQuant(grid, pdf, -6, 8)


  
  pdfname = string_val_opt("-pdf")
  imem    = int_val_opt("-imem",0)
  Q       = dble_val_opt("-Q", 2.0_dp)
  call InitPDFsetByName(trim(pdfname))
  call InitPDFm(1,imem)
  call fill_from_lhapdf(pdf, Q)

  moments = TruncatedMoment(grid, pdf(:,:), one)
  write(6,*) moments
  write(6,*) sum(moments) - one
contains

  !----------------------------------------------------------------------
  subroutine fill_from_lhapdf(this_pdf, Q)
    real(dp), intent(out) :: this_pdf(0:,-6:)
    real(dp), intent(in)  :: Q
    !---------
    real(dp) :: xv(0:grid%ny)
    integer  :: i

    xv = xValues(grid)
    this_pdf(:,7) = zero
    do i = 0, grid%ny
       call EvolvePDFPhoton(xv(i), Q, this_pdf(i,-6:6), this_pdf(i,8))
       !write(6,*) xv(i), this_pdf(i,8)
    end do
    
  end subroutine fill_from_lhapdf
  
end program test_qed_lhapdf
