! An example in Fortran90 showing how to replace LHAPDF interpolation with
! (faster) hoppet interpolation using the streamlined interface
!

module time_format
  use hoppet_term
  implicit none

  character(len=*), parameter :: tfmt = '("'//blue//'",a,f8.2,a,"'//reset//'")'
end module time_format


program fast_manypdf_evaluation
  use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF ! Avoid namespace clashes with LHAPDF
  use streamlined_interface, InitPDF_hoppet => InitPDF
  use hoppet_lhapdf
  use hoppet_to_string
  use time_format
  use new_as
  use pdf_tabulate

  implicit none
  character(len=200) :: pdfname
  integer :: imem
  real(dp) :: x, Q, lhapdf(-6:6), hoppetpdf(-6:6)
  real(dp), allocatable :: hoppetpdf1D(:,:), values(:)
  integer, parameter :: npoints = 500
  real(dp), parameter :: xmin = 1d-5, xmax = 0.9_dp, Qmin = 1.5_dp, Qmax = 100000.0_dp
  real(dp) :: central, err_up, err_down, err_symm
  real(dp), allocatable :: xvals(:), qvals(:)
  real(dp), external :: hoppetAlphaS, alphasPDF, hoppetEvalPID
  real(dp) :: aslhapdf, ashoppet
  integer :: i, j
  character(len=*), parameter :: dfmt = '(f10.6)'
  real(dp) :: t1, t2
  external :: InitPDF

  ! Interface to LHAPDF as needed by hoppetAssign

  pdfname = "PDF4LHC21_40"
  call LoadLHAPDFSet(pdfname)
  !call load_lhapdf_assign_many_hoppet(trim(pdfname))
  x = 0.01_dp
  Q = 13.0_dp
 
  ! Standard LHAPDF call
  call initPDF(0) ! load the central member
  call EvolvePDF(x, Q, lhapdf) ! Get the PDF from LHAPDF

  ! Equivalent hoppet call
  call EvalPdfTable_xQ(lhapdf_set%tables(0), x, Q, hoppetpdf)  
  !call hoppetEval(x, Q, hoppetpdf)

  ! Print the two PDFs
  write(*,'(a)') bold // "PDFs at x = " // to_string(x,dfmt) // ", Q = " // to_string(Q,dfmt) // reset
  write(*,'(a,11f12.8)') "LHAPDF(-5:5): ", lhapdf(-5:5)
  write(*,'(a,11f12.8)') "Hoppet(-5:5): ", hoppetpdf(-5:5)
  write(*,'(a,11f12.8)') "Difference  : ", hoppetpdf(-5:5)-lhapdf(-5:5)


  !!------- check that we get the same uncertainty as direct calls to LHAPDF
  allocate(hoppetpdf1D(-6:6,0:ubound(lhapdf_set%tables,1) ))
  allocate(values(0:ubound(lhapdf_set%tables,1)))
  call EvalPdfTable_xQ(lhapdf_set%tables(:), x, Q, hoppetpdf1D(:,:))
  call GetPDFUncertainty(hoppetpdf1D(iflv_g,:), central, err_up, err_down, err_symm)
  do imem = 0, ubound(values,1)
    call initPDF(imem)
    call EvolvePDF(x, Q, lhapdf)
    values(imem) = lhapdf(iflv_g)
  end do
  call GetPDFUncertainty(values(0:), central, err_up, err_down, err_symm)
  write(*,'(a)') 'LHAPDF gluon: central = '//to_string(central,dfmt)// &
                ', err_up = '//to_string(err_up,dfmt)// &
                ', err_down = '//to_string(err_down,dfmt)// &
                ', err_symm = '//to_string(err_symm,dfmt)

  write(*,'(a)') 'HOPPET gluon: central = '//to_string(central,dfmt)// &
                ', err_up = '//to_string(err_up,dfmt)// &
                ', err_down = '//to_string(err_down,dfmt)// &
                ', err_symm = '//to_string(err_symm,dfmt)

  write(*,*) ! a blank line for clarity

  !!------- carry out timing tests
  allocate(xvals(npoints), qvals(npoints))
  do i = 1, npoints
    xvals(i) = exp(log(xmin) + (i-1) * (log(xmax) - log(xmin)) / real(npoints-1, dp))
    qvals(i) = exp(log(Qmin) + (i-1) * (log(Qmax) - log(Qmin)) / real(npoints-1, dp))
  end do

  ! Benchmark hoppetEval
  ! All flavours at once
  call cpu_time(t1)
  write(*,'(a,i4,a)') bold//"Timing evaluation of all", size(lhapdf_set%tables), " PDF members with HOPPET"//reset
  do i = 1, npoints
    do j = 1, npoints
      do imem = 0, ubound(lhapdf_set%tables,1)
        call EvalPdfTable_xQ(lhapdf_set%tables(imem), xvals(i), qvals(j), hoppetpdf)
      end do
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) blue//"EvalPdfTable_xQ(0D) time per mem (all flav): ", (t2-t1)/npoints/npoints*1d9/size(lhapdf_set%tables), " ns"//reset


  call cpu_time(t1)
  !write(*,*) "Benchmarking all", size(many_tables), " PDF members"
  do i = 1, npoints
    do j = 1, npoints
      call EvalPdfTable_xQ(lhapdf_set%tables(:), xvals(i), qvals(j), hoppetpdf1D(:,:))
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) blue//"EvalPdfTable_xQ(1D) time per mem (all flav): ", (t2-t1)/npoints/npoints*1d9/size(lhapdf_set%tables), " ns"//reset


  !! ! One flavour at a time
  !! call cpu_time(t1)
  !! do i = 1, npoints
  !!   do j = 1, npoints
  !!     hoppetpdf(0) = hoppetEvalPID(xvals(i), qvals(j), mod(j,6))
  !!     ! if pdf_tabulate is "use"d, this is slightly faster
  !!     !hoppetpdf(0) = EvalPdfTable_xQf(tables(0),xvals(i),qvals(j),mod(j,6))
  !!   end do
  !! end do
  !! call cpu_time(t2)
  !! write(*,tfmt) "hoppetEvalPID time (one flav): ", (t2-t1)/npoints/npoints*1d9, " ns"

  ! Benchmark EvolvePDF (LHAPDF)
  write(*,'(a,i4,a)') bold//"Timing evaluation of all", size(lhapdf_set%tables), " PDF members with LHAPDF"//reset
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      do imem = 0, ubound(lhapdf_set%tables,1)
        call InitPDF(imem) ! load the table into hoppet
        call EvolvePDF(xvals(i), qvals(j), lhapdf)
      end do
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) "LHAPDF EvolvePDF time per mem (all flav): ", (t2-t1)/npoints/npoints*1d9/size(lhapdf_set%tables), " ns"

  ! Now let's check the timing of hoppetAlphaS vs LHAPDF alphaS
  write(*,'(a)') bold//"Timing alpha_s"//reset

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
contains
  ! Routine that loads an LHAPDF set, extracts some information from
  ! it and transfers the PDF to hoppet. 

end program fast_manypdf_evaluation
