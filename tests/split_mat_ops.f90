!======================================================================
!! Program for checking that operations on splitting matrices behave
!! as expected. 
program split_mat_ops
  use hoppet
  implicit none
  type(grid_def)     :: grid
  type(split_mat)    :: PA, PB, PC
  type(dglap_holder) :: dh
  real(dp), pointer  :: pdf(:,:), pdfres(:,:)
  integer i

  ! choose grid in which multiplications of grid_conv are exactly
  ! commutative, and in which PA*PB*pdf is exactly associative (both
  ! to within machine precision) -- negative order ensures the former,
  ! the fact that we have a single grid (rather than subgrids) ensures
  ! the latter
  call InitGridDef(grid, dy=0.1_dp, ymax = 10.0_dp, order=-5)

  ! get some convenient ready-made splitting matrices
  do i = 1, 3
     call InitDglapHolder(grid, dh, nloop=2)
     call delete(dh)
  end do
  
  call InitDglapHolder(grid, dh, nloop=2)
  call AllocPDF(grid, pdf)
  pdf = unpolarized_dummy_pdf(xValues(grid))

  ! allocate space for result
  call AllocPDF(grid, pdfres)

  write(6,*) "The following tests are passed if the results are all"
  write(6,*) "zero to within machine precision"
  
  ! multiplication test
  call InitSplitMat(PA, dh%P_LO); call Multiply(PA,two)
  call InitSplitMat(PB, dh%P_LO, two)
  pdfres = PA*pdf - PB*pdf
  write(6,*) 'Multiplication test: ', maxval(pdfres(:,-6:6))

  ! convolution test
  call SetToConvolution(PC, dh%P_LO, dh%P_NLO)
  pdfres = PC*pdf - dh%P_LO*(dh%P_NLO*pdf)
  write(6,*) '(PA*PB)*pdf - PA*(PB*pdf): ', maxval(pdfres(:,-6:6))
  
  !! commutation test
  call SetToCommutator(PC, dh%P_LO, dh%P_NLO)
  pdfres = PC*pdf - (dh%P_LO*(dh%P_NLO*pdf)  - dh%P_NLO*(dh%P_LO*pdf))
  write(6,*) '[PA,PB]*pdf - (PA*(PB*pdf)-PB*(PA*pdf)): ', maxval(pdfres(:,-6:6))

  ! cleaning up
  call delete(PA)
  call delete(PB)
  call delete(PC)
  call delete(dh)
  call delete(grid)
  deallocate(pdf, pdfres)

contains

  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolazrized evolution
  function unpolarized_dummy_pdf(xvals) result(dummy_pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: dummy_pdf(size(xvals),-6:7)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    dummy_pdf = zero ! automatically in human rep

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)
    dummy_pdf(:,iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
        
    dummy_pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    dummy_pdf(:, iflv_s) = dummy_pdf(:,-iflv_s)
    dummy_pdf(:, iflv_u) = uv + ubar
    dummy_pdf(:,-iflv_u) = ubar
    dummy_pdf(:, iflv_d) = dv + dbar
    dummy_pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf

end program split_mat_ops
