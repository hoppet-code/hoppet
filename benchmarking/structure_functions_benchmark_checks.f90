!! An example program which prints out the DIS structure functions up
!! to N3LO. This program produces the checks that were carried out
!! with Valerio Bertone in (INSERT ARXIV IDENTIFIER HERE)
!!

program structure_functions_benchmark_checks
  use hoppet
  use pdfs_for_benchmarks
  use streamlined_interface
  use structure_functions
  use io_utils !< for command-line access
  implicit none
  real(dp) :: Qmax, xmur, xmuf, Qmin, ymax, Q, mc, mb, mt, asQ, Q0alphas, muR_Q, Q0pdf
  real(dp) :: xpdf_at_xQ(-6:6), dy, dlnlnQ, minQval, maxQval
  integer  :: order_max, sc_choice,nloop,ix,nf_lcl,order
  real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  ! Couplings needed for the F3 structure functions
  real(dp), parameter :: mz = 91.1876_dp, mw = 80.398_dp, sin_thw_sq = 1.0_dp - (mw/mz)**2
  real(dp), parameter :: Ae = - 0.5_dp, Ve = - 0.5_dp + 2.0_dp * sin_thw_sq, Ae2 = Ae**2, &
       &  Ve2 = Ve**2, Ve2_Ae2 = Ve2 + Ae2, two_Ve_Ae = 2.0_dp * Ve * Ae, &
       &  sin_2thw_sq = 4.0_dp * (1.0_dp - sin_thw_sq) * sin_thw_sq
  logical :: param_coefs = .true.
  character(len=1024) :: outdir

  Qmax = 13000.0_dp 
  order_max = 4
  xmur = one
  xmuf = one
  sc_choice = 1 ! Uses Q as the central scale choice
  Qmin = one

  ! Set heavy flavour scheme
  mc = 1.414213563_dp   ! sqrt(2.0_dp) + epsilon
  mb = 4.5_dp
  mt = 175.0_dp
  call hoppetSetPoleMassVFN(mc, mb, mt)

  ! Streamlined initialization
  ! including  parameters for x-grid
  order = -6 
  ymax  = 16.0_dp
  dy    = dble_val_opt("-dy",0.05_dp)
  param_coefs = .not. log_val_opt("-exact-coefs")
  outdir = string_val_opt("-outdir","structure-functions-benchmarks")
  dlnlnQ = dy/4.0_dp
  nloop = 3 
  minQval = min(xmuF*Qmin, Qmin)
  maxQval = max(xmuF*Qmax, Qmax)
  
  ! initialise the grid and dglap holder
  write(6,'(a)') "Initialising splitting functions"
  call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
       &         order,factscheme_MSbar)

  ! Setup all constants and parameters needed by the structure functions
  write(6,'(a)') "Initialising coefficient functions for structure functions"
  if (.not. param_coefs) write(6,('(a)')) "NB: this is being done with exact coefficient functions and will take about a minute "
  call StartStrFct(order_max, scale_choice = sc_choice, param_coefs =&
       & param_coefs, wmass = mw, zmass = mz)

  ! Evolve the PDF
  nloop = 3 ! NNLO evolution
  asQ = 0.35_dp
  Q0alphas = sqrt(2.0_dp)
  muR_Q = 1.0_dp
  Q0pdf = sqrt(2.0_dp) ! The initial evolution scale
  call hoppetEvolve(asQ, Q0alphas, nloop,muR_Q, benchmark_pdf_unpolarized_lha, Q0pdf)

  write(6,'(a)') "Quick output of alphas and PDFs, e.g. to enable consistency checks with other codes"
  Q = 100.0_dp
  print*, 'Q = ', Q, ' GeV'
  !AV commented out because it is an error in non-GNU compilers
  !print*, 'αS(Q) = ', Value(coupling, Q)
  
  write(6,'(a)')
  write(6,'(a,f8.3,a)') "Evaluating PDFs at Q = ",Q," GeV"
  write(6,'(a5,2a12,a14,a10,a12)') "x",&
       & "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
  do ix = 1, size(heralhc_xvals)
     call EvalPdfTable_xQ(tables(0),heralhc_xvals(ix),Q,xpdf_at_xQ)
     write(6,'(es7.1,5es12.4)') heralhc_xvals(ix), &
          &  xpdf_at_xQ(2)-xpdf_at_xQ(-2), &
          &  xpdf_at_xQ(1)-xpdf_at_xQ(-1), &
          &  2*(xpdf_at_xQ(-1)+xpdf_at_xQ(-2)), &
          &  (xpdf_at_xQ(-4)+xpdf_at_xQ(4)), &
          &  xpdf_at_xQ(0)
  end do
  print*, ''
  
  ! Initialise the structure functions using separate order
  write(6,'(a)') "Initialising tabulated structure functions"
  call InitStrFct(order_max, .true., xR = xmur, xF = xmuf)

  write(6,'(a)') "Writing structure function benchmarks to ** " // trim(outdir) // "/ ** directory"
  open(unit = 99, file = trim(outdir) // '/structure-functions-Q-2.0-GeV.dat')
  Q = 2.0_dp
  write(6,'(a,f5.1,a)') "Printing structure functions at Q = ",Q," GeV"
  call write_f1_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f2_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f3_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  close(99)
  
  open(unit = 99, file = trim(outdir) // '/structure-functions-Q-5.0-GeV.dat')
  Q = 5.0_dp
  write(6,'(a,f5.1,a)') "Printing structure functions at Q = ",Q," GeV"
  call write_f1_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f2_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f3_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  close(99)
  
  open(unit = 99, file = trim(outdir) // '/structure-functions-Q-10.0-GeV.dat')
  Q = 10.0_dp
  write(6,'(a,f5.1,a)') "Printing structure functions at Q = ",Q," GeV"
  call write_f1_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f2_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f3_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  close(99)
  
  open(unit = 99, file = trim(outdir) // '/structure-functions-Q-50.0-GeV.dat')
  Q = 50.0_dp
  write(6,'(a,f5.1,a)') "Printing structure functions at Q = ",Q," GeV"
  call write_f1_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f2_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f3_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  close(99)
  
  open(unit = 99, file = trim(outdir) // '/structure-functions-Q-100.0-GeV.dat')
  Q = 100.0_dp
  write(6,'(a,f5.1,a)') "Printing structure functions at Q = ",Q," GeV"
  call write_f1_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f2_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f3_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  close(99)
  
  open(unit = 99, file = trim(outdir) // '/structure-functions-Q-500.0-GeV.dat')
  Q = 500.0_dp
  write(6,'(a,f5.1,a)') "Printing structure functions at Q = ",Q," GeV"
  call write_f1_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f2_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  call write_f3_benchmark(99, Q, -3.0_dp, 10.0_dp, 200)
  close(99)
  
contains 
 !----------------------------------------------------------------------
  ! write the F1 structure function to idev
  subroutine write_f1_benchmark (idev, Qtest, logxomx_min, logxomx_max, nx)
  implicit none
    real(dp), intent(in) :: Qtest, logxomx_min, logxomx_max
    integer, intent(in)  :: idev, nx
    real(dp) :: ytest, xval, mR, mF, F1Z_LO, F1Z_NLO, F1Z_NNLO, F1Z_N3LO, res(-6:7)
    real(dp) :: F1g_LO, F1g_NLO, F1g_NNLO, F1g_N3LO
    integer  :: iy
    !F1 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    write(idev,'(a,a)') '# x  F1Wp(LO) F1Wm(LO) F1Wp(NLO) F1Wm(NLO) F1Wp(NNLO) F1Wm(NNLO)', &
          & ' F1Wp(N3LO) F1Wm(N3LO) F1Z(LO) F1Z(NLO) F1Z(NNLO) F1Z(N3LO) F1γ(LO) F1γ(NLO) F1γ(NNLO) F1γ(N3LO)'
    mF = sf_muF(Qtest)
    mR = sf_muR(Qtest)
    do iy = nx, 0, -1
       ytest = logxomx_max + iy * (logxomx_min - logxomx_max) / nx
       xval = exp(-ytest)/(one+exp(-ytest))
       ytest = - log(xval)
       res = F_LO(xval, Qtest, mR, mF)
       write(idev,'(3es22.12)',advance='no') xval, res(iF1Wp),res(iF1Wm)
       F1Z_LO = res(iF1Z)
       F1g_LO = res(iF1EM)
       res = F_NLO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF1Wp), res(iF1Wm)
       F1Z_NLO = res(iF1Z)
       F1g_NLO = res(iF1EM)
       res = F_NNLO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF1Wp), res(iF1Wm)
       F1Z_NNLO = res(iF1Z)
       F1g_NNLO = res(iF1EM)
       res = F_N3LO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF1Wp), res(iF1Wm)
       F1Z_N3LO = res(iF1Z)
       F1g_N3LO = res(iF1EM)
       write(idev,'(4es22.12)',advance='no') F1Z_LO, F1Z_NLO, F1Z_NNLO, F1Z_N3LO
       write(idev,'(4es22.12)',advance='no') F1g_LO, F1g_NLO, F1g_NNLO, F1g_N3LO
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f1_benchmark

  !----------------------------------------------------------------------
  ! write the F2 structure function to idev
  subroutine write_f2_benchmark (idev, Qtest, logxomx_min, logxomx_max, nx)
  implicit none
    real(dp), intent(in) :: Qtest, logxomx_min, logxomx_max
    integer, intent(in)  :: idev, nx
    real(dp) :: ytest, xval, mR, mF, F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO, res(-6:7)
    real(dp) :: F2g_LO, F2g_NLO, F2g_NNLO, F2g_N3LO
    integer  :: iy
    !F2 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    write(idev,'(a,a)') '# x  F2Wp(LO) F2Wm(LO) F2Wp(NLO) F2Wm(NLO) F2Wp(NNLO) F2Wm(NNLO)', &
          & ' F2Wp(N3LO) F2Wm(N3LO) F2Z(LO) F2Z(NLO) F2Z(NNLO) F2Z(N3LO) F2γ(LO) F2γ(NLO) F2γ(NNLO) F2γ(N3LO)'
    mF = sf_muF(Qtest)
    mR = sf_muR(Qtest)
    do iy = nx, 0, -1
       ytest = logxomx_max + iy * (logxomx_min - logxomx_max) / nx
       xval = exp(-ytest)/(one+exp(-ytest))
       ytest = - log(xval)
       res = F_LO(xval, Qtest, mR, mF)
       write(idev,'(3es22.12)',advance='no') xval, res(iF2Wp),res(iF2Wm)
       F2Z_LO = res(iF2Z)
       F2g_LO = res(iF2EM)
       res = F_NLO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF2Wp), res(iF2Wm)
       F2Z_NLO = res(iF2Z)
       F2g_NLO = res(iF2EM)
       res = F_NNLO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF2Wp), res(iF2Wm)
       F2Z_NNLO = res(iF2Z)
       F2g_NNLO = res(iF2EM)
       res = F_N3LO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF2Wp), res(iF2Wm)
       F2Z_N3LO = res(iF2Z)
       F2g_N3LO = res(iF2EM)
       write(idev,'(4es22.12)',advance='no') F2Z_LO, F2Z_NLO, F2Z_NNLO, F2Z_N3LO
       write(idev,'(4es22.12)',advance='no') F2g_LO, F2g_NLO, F2g_NNLO, F2g_N3LO
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f2_benchmark

  !----------------------------------------------------------------------
  ! write the F3 structure function to idev
  subroutine write_f3_benchmark (idev, Qtest, logxomx_min, logxomx_max, nx)
  implicit none
    real(dp), intent(in) :: Qtest, logxomx_min, logxomx_max
    integer, intent(in)  :: idev, nx
    real(dp) :: ytest, xval, mR, mF, F3Z_LO, F3Z_NLO, F3Z_NNLO, F3Z_N3LO, res(-6:7)
    real(dp) :: F3gZ_LO, F3gZ_NLO, F3gZ_NNLO, F3gZ_N3LO, F3LO, F3NLO, F3NNLO, F3N3LO, &
            & prop
    integer  :: iy
    !F3 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    write(idev,'(a,a)') '# x  F3Wp(LO) F3Wm(LO) F3Wp(NLO) F3Wm(NLO) F3Wp(NNLO) F3Wm(NNLO)', &
          & ' F3Wp(N3LO) F3Wm(N3LO) F3Z(LO) F3Z(NLO) F3Z(NNLO) F3Z(N3LO)'
    mF = sf_muF(Qtest)
    mR = sf_muR(Qtest)
    do iy = nx, 0, -1
       ytest = logxomx_max + iy * (logxomx_min - logxomx_max) / nx
       xval = exp(-ytest)/(one+exp(-ytest))
       ytest = - log(xval)
       res = F_LO(xval, Qtest, mR, mF)
       write(idev,'(3es22.12)',advance='no') xval, res(iF3Wp),res(iF3Wm)
       F3Z_LO = res(iF3Z)
       F3gZ_LO = res(iF3gZ)
       res = F_NLO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF3Wp), res(iF3Wm)
       F3Z_NLO = res(iF3Z)
       F3gZ_NLO = res(iF3gZ)
       res = F_NNLO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF3Wp), res(iF3Wm)
       F3Z_NNLO = res(iF3Z)
       F3gZ_NNLO = res(iF3gZ)
       res = F_N3LO(xval, Qtest, mR, mF)
       write(idev,'(2es22.12)',advance='no') res(iF3Wp), res(iF3Wm)
       F3Z_N3LO = res(iF3Z)
       F3gZ_N3LO = res(iF3gZ)
       ! Before printing we dress the F3 structure functions with the appropriate couplings 
       ! as defined in for instance the pink book (notice that there are typos in these 
       ! equations in the pink book) eqs. 4.20-4.21
       prop   = Qtest**2 / (Qtest**2 + MZ**2) / sin_2thw_sq
       F3LO   = - (Ae * prop * F3gZ_LO   - two_Ve_Ae * prop**2 * F3Z_LO)   * 0.5_dp
       F3NLO  = - (Ae * prop * F3gZ_NLO  - two_Ve_Ae * prop**2 * F3Z_NLO)  * 0.5_dp
       F3NNLO = - (Ae * prop * F3gZ_NNLO - two_Ve_Ae * prop**2 * F3Z_NNLO) * 0.5_dp
       F3N3LO = - (Ae * prop * F3gZ_N3LO - two_Ve_Ae * prop**2 * F3Z_N3LO) * 0.5_dp
       write(idev,'(4es22.12)',advance='no') F3LO, F3NLO, F3NNLO, F3N3LO
       write(idev,*)
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_f3_benchmark

end program structure_functions_benchmark_checks
