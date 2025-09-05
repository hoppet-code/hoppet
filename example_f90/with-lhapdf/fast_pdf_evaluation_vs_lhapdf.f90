! An example in Fortran90 showing how to replace LHAPDF interpolation with
! (faster) hoppet interpolation using the streamlined interface
!
program fast_pdf_evaluation_vs_lhapdf
  use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF ! Avoid namespace clashes with LHAPDF
  use streamlined_interface
  implicit none
  character(len=200) :: pdfname
  integer :: imem
  real(dp) :: x, Q, lhapdf(-6:6), hoppetpdf(-6:6)
  integer, parameter :: npoints = 3000
  real(dp), parameter :: xmin = 1d-5, xmax = 0.9_dp, Qmin = 1.5_dp, Qmax = 100000.0_dp
  real(dp), allocatable :: xvals(:), qvals(:)
  integer :: i, j
  real(dp) :: t1, t2

  ! Interface to LHAPDF as needed by hoppetAssign

  pdfname = "PDF4LHC21_40"
  imem = 0
  call load_lhapdf_assign_hoppet(trim(pdfname), imem)

  x = 0.01
  Q = 13.0
 
  ! Standard LHAPDF call
  call EvolvePDF(x, Q, lhapdf) ! Get the PDF from LHAPDF

  ! Equivalent hoppet call
  call hoppetEval(x, Q, hoppetpdf)

  ! Print the two PDFs
  write(*,*) "PDFs at x = ", x, ", Q = ", Q
  write(*,*) "LHAPDF(-5:5): ", lhapdf(-5:5)
  write(*,*) "Hoppet(-5:5): ", hoppetpdf(-5:5)

  allocate(xvals(npoints), qvals(npoints))
  do i = 1, npoints
    xvals(i) = exp(log(xmin) + (i-1) * (log(xmax) - log(xmin)) / real(npoints-1, dp))
    qvals(i) = exp(log(Qmin) + (i-1) * (log(Qmax) - log(Qmin)) / real(npoints-1, dp))
  end do

  ! Benchmark hoppetEval
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      call hoppetEval(xvals(i), qvals(j), hoppetpdf)
    end do
  end do
  call cpu_time(t2)
  write(*,*) "hoppetEval time: ", (t2-t1)/npoints/npoints*1d9, " ns"

  ! Benchmark EvolvePDF (LHAPDF)
  call cpu_time(t1)
  do i = 1, npoints
    do j = 1, npoints
      call EvolvePDF(xvals(i), qvals(j), lhapdf)
    end do
  end do
  call cpu_time(t2)
  write(*,*) "LHAPDF time: ", (t2-t1)/npoints/npoints*1d9, " ns"
  
contains
  ! Routine that loads an LHAPDF set, extracts some information from
  ! it and transfers the PDF to hoppet. 
  subroutine load_lhapdf_assign_hoppet(pdfname, imem)
    use streamlined_interface
    character(len=*), intent(in) :: pdfname
    integer, intent(in) :: imem
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
    integer :: orderPDF, nloop, order, yorder, lnlnQorder

    ! Load LHAPDF set
    call initPDFSetByName(pdfname)

    call initPDFSetByName(pdfname)
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
    call hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme_MSbar) ! Start hoppet

    ! Now we fill the hoppet grid using the LHAPDF grid directly,
    ! rather than evolving ourselves
    Q0 = Qmin
    call hoppetAssignWithCoupling(EvolvePDF, alphasPDF(Q0), Q0, nloop)

    ! If instead we want to evolve the PDF with hoppet starting from
    ! some low scale Q0 (>= Qmin) make a call to hoppetEvolve instead
    ! of hoppetAssign
    ! call hoppetEvolve(alphasPDF(Q0), Q0, nloop, 1.0, EvolvePDF, Q0)

  end subroutine load_lhapdf_assign_hoppet
end program fast_pdf_evaluation_vs_lhapdf
