!! Module containing dummy PDFs that are widely used for benchmarking
!! and other tests. Used notably as part of hep-ph/0511119 and 
!! follow-up work.
module pdfs_for_benchmarks
  use types; use consts_dp; use pdf_representation
  implicit none
  private

  public :: benchmark_pdf_unpolarized, benchmark_pdf_unpolarized_lha
  
contains
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution, subsequently adopted in a wide range of
  !! of benchmarks
  pure subroutine benchmark_pdf_unpolarized_lha(x, Q, xpdf)
    real(dp), intent(in)  :: x,Q
    real(dp), intent(out) :: xpdf(-6:6)
    real(dp) :: uv, dv
    real(dp) :: ubar, dbar
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    ! Set to zero the xpdf array
    xpdf = zero
    
  
    !-- remember that these are all x*q(x)
    uv = N_uv * x**0.8_dp * (1-x)**3
    dv = N_dv * x**0.8_dp * (1-x)**4
    dbar = N_db * x**(-0.1_dp) * (1-x)**6
    ubar = dbar * (1-x)
  
    ! labels iflv_g, etc., come from the hoppet module, inherited
    ! from the main program
    xpdf(iflv_g) = N_g * x**(-0.1_dp) * (1-x)**5
    xpdf(-iflv_s) = 0.2_dp*(dbar + ubar)
    xpdf(iflv_s) = xpdf(-iflv_s)
    xpdf(iflv_u) = uv + ubar
    xpdf(-iflv_u) = ubar
    xpdf(iflv_d) = dv + dbar
    xpdf(-iflv_d) = dbar
  end subroutine benchmark_pdf_unpolarized_lha
  

  !----------------------------------------------------------------------
  ! the same dummy PDF, in a slightly different interface
  function benchmark_pdf_unpolarized(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmax)
    !---------------------
    integer  :: iy
    real(dp) :: dummy_Q

    pdf = zero
    call LabelPdfAsHuman(pdf)
    dummy_Q = 10.0_dp

    do iy = 1, size(xvals)
       call benchmark_pdf_unpolarized_lha(xvals(iy), dummy_Q, pdf(iy, iflv_min:iflv_max))
    end do
  end function benchmark_pdf_unpolarized

end module pdfs_for_benchmarks
