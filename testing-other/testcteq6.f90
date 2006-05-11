program testcteq6
  use types; use consts_dp
  implicit none
  !--- user info
  external :: evolvePDF
  real(dp) :: alphasPDF
  real(dp) :: Q(3), x, lha(-6:6,3), pdf(-6:6), pdf_deriv(-6:6), dummy(-6:6)
  real(dp) :: pdf_deriv_lo(-6:6), pdf_deriv_nlo(-6:6)
  real(dp) :: lha_deriv(-6:6)
  real(dp) :: alphas
  integer  :: pdf_handle, P_LO_handle, P_NLO_handle
  external :: pdfconv_new_pdf, pdfconv_P_LO, pdfconv_P_NLO
  integer  :: pdfconv_new_pdf, pdfconv_P_LO, pdfconv_P_NLO  ! return pdf handles
  integer  :: i

  ! select cteq6m1
  call InitPDFsetByName("cteq61.LHgrid")
  !call InitPDFsetByName("cteq6l.LHpdf")
  call InitPDF(0)

  ! start our 
  call pdfconv_start(0.1_dp, 5)  ! dy and nf

  x = 0.1_dp
  Q(:) = 30.0_dp*(/ 1.01_dp,1.0_dp,1.0_dp/1.01_dp/)
  do i = 1, 3
     call evolvePDF(x,Q(i),lha(:,i))
  end do
  write(0,*) real(lha(:,2))

  pdf_handle = pdfconv_new_pdf(evolvePDF,Q(2))
  call pdfconv_eval_pdf(pdf_handle, x, pdf)
  write(0,*) real(pdf)


  write(0,*) "deriv:"
  lha_deriv = (lha(:,3)-lha(:,1))/(two*log(Q(3)/Q(1)))
  write(0,*) real(lha_deriv)


  P_LO_handle = pdfconv_P_LO(pdf_handle)
  P_NLO_handle = pdfconv_P_NLO(pdf_handle)
  call pdfconv_eval_pdf(P_LO_handle, x, pdf_deriv_lo)
  call pdfconv_eval_pdf(P_NLO_handle, x, pdf_deriv_nlo)
  pdf_deriv = pdf_deriv_lo * alphasPDF(Q(2))/twopi + &
       &      pdf_deriv_nlo * (alphasPDF(Q(2))/twopi)**2
  write(0,*) alphasPDF(Q(2))
  write(0,*) real(pdf_deriv)


  call pdfconv_release_pdf(pdf_handle)
  call pdfconv_release_pdf(P_LO_handle)
end program testcteq6
