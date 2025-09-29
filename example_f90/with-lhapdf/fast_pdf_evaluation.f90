! An example in Fortran90 showing how to replace LHAPDF interpolation with
! (faster) hoppet interpolation using the streamlined interface
!
program fast_pdf_evaluation
  use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF ! Avoid namespace clashes with LHAPDF
  use streamlined_interface
  use hoppet_term
  use hoppet_lhapdf
  use hoppet_to_string
  use new_as
  use pdf_tabulate
  implicit none
  character(len=200) :: pdfname
  integer :: imem
  real(dp) :: x, Q, lhapdf(-6:6), hoppetpdf(-6:6)
  integer, parameter :: npoints = 3000
  real(dp), parameter :: xmin = 1d-5, xmax = 0.9_dp, Qmin = 1.5_dp, Qmax = 100000.0_dp
  real(dp), allocatable :: xvals(:), qvals(:)
  real(dp), external :: hoppetAlphaS, alphasPDF, hoppetEvalPID
  real(dp) :: aslhapdf, ashoppet
  integer :: i, j
  real(dp) :: t1, t2
  character(len=*), parameter :: tfmt = '("'//blue//'",a,f8.2,a,"'//reset//'")'
  character(len=*), parameter :: dfmt = '(f10.6)'

  ! Interface to LHAPDF as needed by hoppetAssign

  pdfname = "PDF4LHC21_40"
  imem = 0
  call LoadLHAPDF(pdfname,imem)
  !call load_lhapdf_assign_hoppet(trim(pdfname), imem)
  x = 0.01_dp
  Q = 13.0_dp
 
  ! Standard LHAPDF call
  call EvolvePDF(x, Q, lhapdf) ! Get the PDF from LHAPDF

  ! Equivalent hoppet call
  call hoppetEval(x, Q, hoppetpdf)

  ! Print the two PDFs
  write(*,'(a)') bold // "PDFs at x = " // to_string(x,dfmt) // ", Q = " // to_string(Q,dfmt) // reset
  write(*,'(a,11f12.8)') "LHAPDF(-5:5): ", lhapdf(-5:5)
  write(*,'(a,11f12.8)') "Hoppet(-5:5): ", hoppetpdf(-5:5)
  write(*,'(a,11f12.8)') "Difference  : ", hoppetpdf(-5:5)-lhapdf(-5:5)
  write(*,'(a)') 'LHAPDF AlphaS(Q) = ' // to_string(alphasPDF(Q),dfmt)
  write(*,'(a)') 'Hoppet AlphaS(Q) = ' // to_string(hoppetAlphaS(Q),dfmt)
  write(*,*) ! a blank line for clarity

  write(6,'(a)') bold // "Timing tests for PDF and alphaS evaluations" // reset
  allocate(xvals(npoints), qvals(npoints))
  do i = 1, npoints
    xvals(i) = exp(log(xmin) + (i-1) * (log(xmax) - log(xmin)) / real(npoints-1, dp))
    qvals(i) = exp(log(Qmin) + (i-1) * (log(Qmax) - log(Qmin)) / real(npoints-1, dp))
  end do

  ! Benchmark hoppetEval
  ! All flavours at once
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      call hoppetEval(xvals(i), qvals(j), hoppetpdf)
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) blue//"hoppetEval    time (all flav): ", (t2-t1)/npoints/npoints*1d9, " ns"//reset
  ! One flavour at a time
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      hoppetpdf(0) = hoppetEvalPID(xvals(i), qvals(j), mod(j,6))
      ! if pdf_tabulate is "use"d, this is slightly faster
      !hoppetpdf(0) = EvalPdfTable_xQf(tables(0),xvals(i),qvals(j),mod(j,6))
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) "hoppetEvalPID time (one flav): ", (t2-t1)/npoints/npoints*1d9, " ns"

  ! Benchmark EvolvePDF (LHAPDF)
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      call EvolvePDF(xvals(i), qvals(j), lhapdf)
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) "LHAPDF EvolvePDF time (all flav): ", (t2-t1)/npoints/npoints*1d9, " ns"

  ! Now let's check the timing of hoppetAlphaS vs LHAPDF alphaS
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      ashoppet = hoppetAlphaS(qvals(j))
      !! GPS M2Pro 2025-09-06: this variant saves only 0.4ns per call
      !! out of about 31.7
      !ashoppet = na_Value(coupling%nah, qvals(j))
      !ashoppet = na_Value_faster(coupling%nah, qvals(j))
      !if (i == 1) write(6,*) qvals(j), hoppetAlphaS(qvals(j)), na_Value_faster(coupling%nah, qvals(j))
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) "hoppet alphaS time: ", (t2-t1)/npoints/npoints*1d9, " ns"
  
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      aslhapdf = alphasPDF(qvals(j))
      !alphasPDF(qvals(j))
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) "LHAPDF alphaS time: ", (t2-t1)/npoints/npoints*1d9, " ns" 
  write(*,*) ! a blank line for clarity

  call hoppetDeleteAll()

end program fast_pdf_evaluation
