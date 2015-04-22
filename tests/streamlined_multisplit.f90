!----------------------------------------------------------------------
! program for testing the streamlined interface for convolution
! with multiple splitting functions
program streamlined_multisplit
  use streamlined_interface
  use dummy_pdfs
  implicit none
  real(dp) :: dy
  integer  :: nloop
  real(dp) :: asQ, Q0alphas, muR_Q = 1.0_dp, Q0pdf
  real(dp) :: mc, mb, mt
  !----------------------------------------------------------------------
  real(dp) :: x = 0.3_dp, Q = 100.0_dp
  real(dp), pointer :: pdf(:,:), pdfconv(:,:)
  real(dp) :: res(-6:6)
  integer  :: nf_lcl
  
  nloop = 3

  call startup()


  ! very basic test ------------
  call hoppetEval(x, Q, res)
  if (nloop == 3) call compare('pdf v. known ref', res(0), 7.8898E-02_dp)

  Q = 100.0_dp
  nf_lcl = -1 ! dynamic nf
  !call hoppetEvalSplit(x, Q, 7, nf_lcl, res) ! for testing errors
  call SetNfDglapHolder(dh, 5)
  call do_full_set

  Q = 3.3_dp
  call SetNfDglapHolder(dh, 4)
  call do_full_set

  nf_lcl = 5 ! fixed nf
  call SetNfDglapHolder(dh, nf_lcl)
  call do_full_set

contains

  ! run a full suite of checks for a given x, Q value
  subroutine do_full_set()

    write(6,*)
    write(6,*) 'Running full set with Q = ',Q, ', x = ', x, 'nf = ', nf_lcl

    call hoppetEval(x, Q, res)
    call EvalPdfTable_Q(tables(0), Q, pdf)
    call compare('internal v. table', res(0), pdf(:,0).atx.(x.with.grid))
    
    call hoppetEvalSplit(x, Q, 1, nf_lcl, res)
    pdfconv = dh%P_LO * pdf
    call compare('PLO * pdf', res(0), pdfconv(:,0).atx.(x.with.grid))
    
    call hoppetEvalSplit(x, Q, 2, nf_lcl, res)
    pdfconv = dh%P_NLO * pdf
    call compare('PNLO * pdf', res(0), pdfconv(:,0).atx.(x.with.grid))
    
    call hoppetEvalSplit(x, Q, 3, nf_lcl, res)
    pdfconv = dh%P_NNLO * pdf
    call compare('PNNLO * pdf', res(0), pdfconv(:,0).atx.(x.with.grid))
    
    call hoppetEvalSplit(x, Q, 11, nf_lcl, res)
    pdfconv = dh%P_LO * (dh%P_LO * pdf)
    call compare('PLO * PLO * pdf', res(0), pdfconv(:,0).atx.(x.with.grid))
    
    call hoppetEvalSplit(x, Q, 12, nf_lcl, res)
    pdfconv = dh%P_LO * (dh%P_NLO * pdf)
    call compare('PLO * PNLO * pdf', res(0), pdfconv(:,0).atx.(x.with.grid))
    
    call hoppetEvalSplit(x, Q, 21, nf_lcl, res)
    pdfconv = dh%P_NLO * (dh%P_LO * pdf)
    call compare('PNLO * PLO * pdf', res(0), pdfconv(:,0).atx.(x.with.grid))
    
    call hoppetEvalSplit(x, Q, 111, nf_lcl, res)
    pdfconv = dh%P_LO * (dh%P_LO * (dh%P_LO * pdf))
    call compare('PLO^3 * pdf', res(0), pdfconv(:,0).atx.(x.with.grid))

  end subroutine do_full_set
  
  
  !----------------------------------------------------------------------
  ! check two numbers are equal
  subroutine compare(testname, v1, v2)
    character(len=*), intent(in) :: testname
    real(dp), intent(in) :: v1, v2
    character(len=4)     :: status
    real(dp), parameter  :: tolerance = 3e-5

    if (abs(v1-v2) > tolerance*max(abs(v1),abs(v2))) then
       status = 'FAIL'
    else
       status = 'OK'
    end if

    write(6,'(a20,2es15.5,a6)') trim(testname), v1, v2, status
  end subroutine compare
  
  
  subroutine startup()
    dy = 0.1_dp
    call hoppetStart(dy, nloop)
    
    ! Set heavy flavour scheme
    mc = 1.414213563_dp   ! sqrt(2.0_dp) + epsilon
    mb = 4.5_dp
    mt = 175.0_dp
    call hoppetSetVFN(mc, mb, mt)
    
    ! Set parameters of running coupling
    asQ = 0.35_dp
    Q0alphas = sqrt(2.0_dp)
    muR_Q = 1.0_dp
    Q0pdf = Q0alphas
    call hoppetEvolve(asQ, Q0alphas, nloop, muR_Q, lha_unpolarized_dummy_pdf, Q0pdf)

    call AllocPDF(grid, pdf)
    call AllocPDF(grid, pdfconv)
  end subroutine startup
  


end program streamlined_multisplit
