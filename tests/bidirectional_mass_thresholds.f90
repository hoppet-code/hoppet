!! Program to test correct implementation of bidirectional mass thresholds.
!!
!! Run it as follows
!!
!! ./bidirectional_mass_thresholds -as 0.005 [-evop] >! bi-c0.005
!! ./bidirectional_mass_thresholds -as 0.01  [-evop] >! bi-c0.01
!!
!! It evolves from 1.4GeV to 1.5GeV (crossing charm threshold) and
!! back. Col  1 = y; 
!!       cols 2-4 = initial (/g,s,c/);
!!       cols 5-7 = final   (/g,s,c/);
!!
!! Given that upwards threshold were working before, test that 
!! downwards ones also work is that the following combinations should
!! all be O(as^4) [since threshold = O(as^2), so non-cancelling piece
!! between up and down directions should be the square of that]
!!
!!       ($5-$2) [initial glue    - final glue]
!!       ($6-$3) [initial strange - final strange]
!!       ($7-$4) [initial charm   - final charm]
!!
!! Note that in the evop formulation, because the result of a
!! convolution has frozen flavours to zero, one loses the left-over
!! "intrinsic charm" component with the -evop option.
!!
!! Tests carried out 7 May 2006 (will be revision 44)
program bidirectional_mass_thresholds
  use hoppet; use io_utils
  implicit none
  type(grid_def) :: grid
  type(dglap_holder) :: dh
  real(dp), pointer  :: pdf_init(:,:), pdf_final(:,:)
  real(dp)           :: Q_init, Q_mid, asQ
  type(running_coupling) :: coupling
  type(evln_operator)    :: evop_up, evop_down
  integer, parameter     :: nloop = 3
  integer  :: iy
  logical  :: use_evop
  real(dp) :: y

  Q_init = 1.4_dp
  Q_mid  = 1.5_dp
  asQ    = dble_val_opt('-as')
  use_evop = log_val_opt('-evop')

  call InitGridDef(grid,ymax=10.0_dp, dy=0.1_dp,order=-5)
  call InitDglapHolder(grid,dh,nloop=nloop,nflo=3,nfhi=4)
  
  call AllocPDF(grid, pdf_init)
  call AllocPDF(grid, pdf_final)
  pdf_init = unpolarized_dummy_pdf(xValues(grid))

  call InitRunningCoupling(coupling, asQ, Q_init, nloop)

  ! go up and down
  if (use_evop) then
     call InitEvlnOperator(dh,evop_up,   coupling, Q_init, Q_mid)
     call InitEvlnOperator(dh,evop_down, coupling, Q_mid, Q_init)
     pdf_final = evop_down * (evop_up * pdf_init)
     call Delete(evop_up)
     call Delete(evop_down)
  else
     pdf_final = pdf_init
     call EvolvePDF(dh, pdf_final, coupling, Q_init, Q_mid)
     call EvolvePDF(dh, pdf_final, coupling, Q_mid, Q_init)
  end if
  
  ! and see what's left
  do iy = 1, 100
     y = iy*0.1_dp
     write(6,'(f10.5,6es25.15)') y, &
          &        pdf_init (:,(/0,3,4/)).aty.(y.with.grid),&
          &        pdf_final(:,(/0,3,4/)).aty.(y.with.grid)
  end do
     
  ! cleanup
  !call Delete(grid)
  call Delete(dh)
  call Delete(coupling)
  call Delete(pdf_init)
  call Delete(pdf_final)

contains
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolazrized evolution
  function unpolarized_dummy_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmax)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    pdf = zero
    ! clean method for labelling as PDF as being in the human representation
    call LabelPdfAsHuman(pdf)

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)
    pdf(:,iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
        
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf

end program bidirectional_mass_thresholds
