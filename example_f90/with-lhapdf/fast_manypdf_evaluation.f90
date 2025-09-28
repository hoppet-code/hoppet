! An example in Fortran90 showing how to replace LHAPDF interpolation with
! (faster) hoppet interpolation using the streamlined interface
!

module time_format
  implicit none
  character(len=*), parameter :: blue = char(27)//'[34m', bold = char(27)//'[1m', reset = char(27)//'[0m'
  character(len=*), parameter :: tfmt = '("'//blue//'",a,f8.2,a,"'//reset//'")'
end module time_format

module many_pdf_tables
  use types
  use pdf_tabulate
  use time_format
  implicit none
  type(pdf_table), allocatable :: many_tables(:)

contains

  subroutine load_lhapdf_assign_many_hoppet(pdfname)
    use streamlined_interface, InitPDF_hoppet => InitPDF 
    character(len=*), intent(in) :: pdfname
    real(dp), external :: alphasPDF
    interface
       subroutine EvolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine EvolvePDF
    end interface
    real(dp) :: mc, mb, mt, Q2minPDF, Q2maxPDF, xmin, xmax, Qmin, QMax, Q0
    real(dp) :: ymax, dy, dlnlnQ
    integer :: orderPDF, nloop, order, yorder, lnlnQorder,nfmax
    real(dp) :: ta1, ta2
    integer :: imem, max_mem
    integer, external :: numberpdf
    external :: initPDFSetByName, initPDF

    ! Load LHAPDF set
    call cpu_time(ta1)
    call initPDFSetByName(pdfname)
    max_mem = numberpdf() - 1
    do imem = 1, max_mem
      call initPDF(imem)
    end do
    call cpu_time(ta2)
    write(*,tfmt) "Time to load LHAPDF sets: ", (ta2-ta1)*1e3_dp, " ms"
    write(*,'(a,i4,a)') "Number of members in set: 1 + ", max_mem

    call getQ2min(0,Q2minPDF)
    call getQ2max(0,Q2maxPDF)
    Qmin = sqrt(Q2minPDF)
    Qmax = sqrt(Q2maxPDF)
    call getxmin(0,xmin)
    call getxmax(0,xmax)
    call getorderas(orderPDF) ! NB: LHAPDF returns 0 for 1-loop running, 1 for 2-loop etc.
    nloop = 1 + orderPDF
    call getthreshold(4,mc)
    call getthreshold(5,mb)
    call getthreshold(6,mt)
    call getnf(nfmax)
    if(nfmax .lt. 6) mt = 2.0d0*Qmax ! If no top in PDF set threshold beyond table max

    write(*,*) "LHAPDF set: ", pdfname, " loaded successfully"

    ! Now let us define some hoppet specific parameters. These are
    ! typical values, and should guarantee similar accuracy as can be
    ! expected from LHAPDF
    ymax = real(ceiling(log(1.0d0/xmin)), kind=dp) ! To get a nice value of ymax that can contain the full LHAPDF grid
    dy = 0.05_dp
    dlnlnQ = dy/4.0_dp
    if(ymax > 15.0) dlnlnQ = dy/8.0_dp ! for large ymax we need a finer grid in Q
    order = -6 ! Default
    yorder = 2 ! Quadratic interpolation in y
    lnlnQorder = 2 ! Quadratic interpolation in lnlnQ

    write(*,*) "Hoppet starting with:"
    write(*,*) " ymax:       ", ymax
    write(*,*) " dy:         ", dy
    write(*,*) " Qmin:       ", Qmin
    write(*,*) " Qmax:       ", Qmax
    write(*,*) " dlnlnQ:     ", dlnlnQ
    write(*,*) " nloop:      ", nloop
    write(*,*) " order:      ", order
    write(*,*) " yorder:     ", yorder
    write(*,*) " lnlnQorder: ", lnlnQorder

    call hoppetSetPoleMassVFN(mc,mb,mt) ! set the pole masses
    call hoppetSetYLnlnQInterpOrders(yorder, lnlnQorder) ! Set the interpolation orders
    call cpu_time(ta1)
    call hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme_MSbar) ! Start hoppet
    call cpu_time(ta2)
    write(*,tfmt) "Time to start HOPPET: ", (ta2-ta1)*1e3_dp, " ms"

    ! Now we fill the hoppet grid using the LHAPDF grid directly,
    ! rather than evolving ourselves
    Q0 = Qmin
    call hoppetSetCoupling(alphasPDF(Q0), Q0, nloop)
    call cpu_time(ta1)
    call hoppetAssign(EvolvePDF)
    call cpu_time(ta2)
    write(*,tfmt) "Time to fill HOPPET grid from LHAPDF: ", (ta2-ta1)*1e3_dp, " ms"
    write(*,*) ! a blank line for clarity

    allocate(many_tables(0:max_mem))
    call AllocPdfTable(many_tables, tables(0)) ! allocate the tables based on tables(0) structure
       
    
    call cpu_time(ta1)
    do imem = 0, max_mem
      write(*,*) "Filling table for member ", imem, " of ", max_mem
      call initPDF(imem)
      call FillPdfTable_LHAPDF(many_tables(imem), EvolvePDF)
      !call tabulatePDFSet(many_tables(imem), pdfname, imem, xmin, xmax, Qmin, Qmax, dy, dlnlnQ)
    end do
    call cpu_time(ta2)
    write(*,tfmt) "Time to fill all sets: ", (ta2-ta1)*1e3_dp, " ms"


    ! If instead we want to evolve the PDF with hoppet starting from
    ! some low scale Q0 (>= Qmin) make a call to hoppetEvolve instead
    ! of hoppetAssign
    ! call hoppetEvolve(alphasPDF(Q0), Q0, nloop, 1.0, EvolvePDF, Q0)

  end subroutine load_lhapdf_assign_many_hoppet
end module many_pdf_tables

program fast_manypdf_evaluation
  use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF ! Avoid namespace clashes with LHAPDF
  use streamlined_interface, InitPDF_hoppet => InitPDF
  use many_pdf_tables
  use time_format
  use new_as
  use pdf_tabulate

  implicit none
  character(len=200) :: pdfname
  integer :: imem
  real(dp) :: x, Q, lhapdf(-6:6), hoppetpdf(-6:6)
  real(dp), allocatable :: hoppetpdf1D(:,:)
  integer, parameter :: npoints = 500
  real(dp), parameter :: xmin = 1d-5, xmax = 0.9_dp, Qmin = 1.5_dp, Qmax = 100000.0_dp
  real(dp), allocatable :: xvals(:), qvals(:)
  real(dp), external :: hoppetAlphaS, alphasPDF, hoppetEvalPID
  real(dp) :: aslhapdf, ashoppet
  integer :: i, j
  real(dp) :: t1, t2
  external :: InitPDF

  ! Interface to LHAPDF as needed by hoppetAssign

  pdfname = "PDF4LHC21_40"
  call load_lhapdf_assign_many_hoppet(trim(pdfname))
  x = 0.01
  Q = 13.0
 
  ! Standard LHAPDF call
  call initPDF(0) ! load the central member
  call EvolvePDF(x, Q, lhapdf) ! Get the PDF from LHAPDF

  ! Equivalent hoppet call
  call EvalPdfTable_xQ(many_tables(0), x, Q, hoppetpdf)  
  !call hoppetEval(x, Q, hoppetpdf)

  ! Print the two PDFs
  write(*,*) "PDFs at x = ", x, ", Q = ", Q
  write(*,'(a,11f12.8)') "LHAPDF(-5:5): ", lhapdf(-5:5)
  write(*,'(a,11f12.8)') "Hoppet(-5:5): ", hoppetpdf(-5:5)
  write(*,'(a,11f12.8)') "Difference  : ", hoppetpdf(-5:5)-lhapdf(-5:5)
  write(*,*) ! a blank line for clarity

  allocate(xvals(npoints), qvals(npoints))
  do i = 1, npoints
    xvals(i) = exp(log(xmin) + (i-1) * (log(xmax) - log(xmin)) / real(npoints-1, dp))
    qvals(i) = exp(log(Qmin) + (i-1) * (log(Qmax) - log(Qmin)) / real(npoints-1, dp))
  end do

  ! Benchmark hoppetEval
  ! All flavours at once
  call cpu_time(t1)
  write(*,'(a,i4,a)') bold//"Benchmarking evaluation of all", size(many_tables), " PDF members with HOPPET"//reset
  do i = 1, npoints
    do j = 1, npoints
      do imem = 0, ubound(many_tables,1)
        call EvalPdfTable_xQ(many_tables(imem), xvals(i), qvals(j), hoppetpdf)
      end do
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) blue//"EvalPdfTable_xQ   time per mem (all flav): ", (t2-t1)/npoints/npoints*1d9/size(many_tables), " ns"//reset


  allocate(hoppetpdf1D(-6:6,0:ubound(many_tables,1) ))
  call cpu_time(t1)
  !write(*,*) "Benchmarking all", size(many_tables), " PDF members"
  do i = 1, npoints
    do j = 1, npoints
      call EvalPdfTable1D_xQ(many_tables(:), xvals(i), qvals(j), hoppetpdf1D(:,:))
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) blue//"EvalPdfTable1D_xQ time per mem (all flav): ", (t2-t1)/npoints/npoints*1d9/size(many_tables), " ns"//reset


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
  write(*,'(a,i4,a)') bold//"Benchmarking evaluation of all", size(many_tables), " PDF members with LHAPDF"//reset
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      do imem = 0, ubound(many_tables,1)
        call InitPDF(imem) ! load the table into hoppet
        call EvolvePDF(xvals(i), qvals(j), lhapdf)
      end do
    end do
  end do
  call cpu_time(t2)
  write(*,tfmt) "LHAPDF EvolvePDF time per mem (all flav): ", (t2-t1)/npoints/npoints*1d9/size(many_tables), " ns"

  ! Now let's check the timing of hoppetAlphaS vs LHAPDF alphaS
  write(*,'(a)') bold//"Benchmarking alpha_s"//reset

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
