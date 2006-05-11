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
  integer  :: i, nf

  ! select cteq6m1
  call InitPDFsetByName("cteq61.LHgrid")
  !call InitPDFsetByName("cteq6l.LHpdf")
  call InitPDF(0)

  ! start our 
  call dglapStart(0.1_dp, 2)  ! dy and nloop
  call dglapAssign(evolvePDF)

  x = 0.1_dp
  Q(:) = 30.0_dp*(/ 1.01_dp,1.0_dp,1.0_dp/1.01_dp/)
  do i = 1, 3
     call evolvePDF(x,Q(i),lha(:,i))
  end do
  write(6,*) real(lha(:,2))

  call dglapEval(x,Q(2),pdf)
  write(6,*) real(pdf)


  write(6,*) "deriv:"
  write(6,*) "deriv:"
  lha_deriv = (lha(:,3)-lha(:,1))/(two*log(Q(3)/Q(1)))
  write(6,*) real(lha_deriv)


  nf = 5
  call dglapEvalSplit(x,Q(2), 1, nf, pdf_deriv_lo)
  call dglapEvalSplit(x,Q(2), 2, nf, pdf_deriv_nlo)
  pdf_deriv = pdf_deriv_lo * alphasPDF(Q(2))/twopi + &
       &      pdf_deriv_nlo * (alphasPDF(Q(2))/twopi)**2
  write(6,*) alphasPDF(Q(2))
  write(6,*) real(pdf_deriv)

end program testcteq6
