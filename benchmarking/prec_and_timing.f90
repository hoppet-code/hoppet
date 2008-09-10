!======================================================================
!
! Program for preparing info for accuracy tests and for carrying out
! timing tests.
!
! A high-accuracy reference grid as used in the documentation can be
! prepared with the instructions
!
!   ./prec_and_timing -nrep 1 -nxQ 5000 -outputgrid -dy 0.025 -order 6 \
!                     -nopreev -4grids -dlnlnQ 0.005 -du 0.005 -olnlnQ 4
!
! Fig. 1 and fig.2(left) of the doc have been obtained by comparing
! this with the output from 
!
!   ./prec_and_timing -nrep 1 -nxQ 5000 -outputgrid -dy 0.2 -order -6 \
!                     -nopreev -4grids -dlnlnQ 0.05 -du 0.4 -olnlnQ 4
!
! using the command
!
!   test_acc/compare2files_v2 resA resB -channel 11 [-protect]
!
! Other options of interest in prec_and_timing include
!
!   -eps [1e-7]       precision of adaptive integration for preparing split.fn.
!   -asdt [0.2]       step for evolution of alpha_s
!   -exactsp          use exact 3-loop split-fns
!   -exactth          use exact 3-loop mass thresholds
!
! BEWARE: older copies of executables have ymax=5 and Qmax=100(?). Explicitly
!         specify -ymax 11.5 and -Qmax 1e4 to get sensible timings.
! 
!======================================================================
module pdf_initial_condition
  use hoppet_v1
  implicit none
  
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


  ! an alternative way of providing the same thing
  subroutine VogtInitSub(y,res)
    real(dp), intent(in)  :: y
    real(dp), intent(out) :: res(ncompmin:)
    real(dp) :: x
    real(dp) :: uv, dv, ubar, dbar
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
    
    res = zero
    x = exp(-y)
    !-- remember that my definitions that these are all x*q(x)
    uv = N_uv * x**0.8_dp * (1-x)**3
    dv = N_dv * x**0.8_dp * (1-x)**4
    dbar = N_db * x**(-0.1_dp) * (1-x)**6
    ubar = dbar * (1-x)
    res(iflv_g) = N_g * x**(-0.1_dp) * (1-x)**5

    res(-iflv_s) = 0.2_dp*(dbar + ubar)
    res( iflv_s) = res(-iflv_s)
    res( iflv_u) = uv + ubar
    res(-iflv_u) = ubar
    res( iflv_d) = dv + dbar
    res(-iflv_d) = dbar
  end subroutine VogtInitSub
end module pdf_initial_condition


!======================================================================
program prec_and_timing
  use hoppet_v1
  use pdf_initial_condition
  use sub_defs_io
  implicit none
  !-------------------------------
  type(grid_def)     :: grid, gridarray(4)
  type(dglap_holder) :: dh
  type(running_coupling)    :: coupling
  integer            :: order, order1, order2, nloop, i, nrep, nxQ,olnlnQ
  integer            :: hires
  real(dp)           :: dy, Qinit, Qmax, du, dlnlnQ
  real(dp)           :: ymax
  real(dp)           :: time_start, time_init_done, time_ev_done, time_end
  real(dp), pointer  :: vogt_init(:,:)
  type(pdf_table)       :: table
  logical            :: output, outputgrid, preev

  ! set the details of the y=ln1/x grid
  dy    = dble_val_opt('-dy',0.25_dp)
  ymax  = dble_val_opt('-ymax',11.5_dp)
  order = int_val_opt('-order',-6)
  order2 = int_val_opt('-order2',order)
  order1 = int_val_opt('-order1',order)

  preev = log_val_opt('-preev',.true.)

  if (log_val_opt('-eps')) call SetDefaultConvolutionEps(dble_val_opt('-eps'))
  if (log_val_opt('-asdt')) call SetDefaultCouplingDt(dble_val_opt('-asdt'))

  if (log_val_opt('-alt7')) then
     call InitGridDef(gridarray(3),dy/7.0_dp, 1.0_dp, order=order2)
     call InitGridDef(gridarray(2),dy/2.0_dp, 3.0_dp, order=order1)
     call InitGridDef(gridarray(1),dy,        ymax, order=order)
     call InitGridDef(grid,gridarray(1:3),locked=.true.)
  else if (log_val_opt('-4grids')) then
     call InitGridDef(gridarray(4),dy/27.0_dp, 0.2_dp, order=order2)
     call InitGridDef(gridarray(3),dy/9.0_dp, 0.5_dp, order=order1)
     call InitGridDef(gridarray(2),dy/3.0_dp, 2.0_dp, order=order )
     call InitGridDef(gridarray(1),dy,        ymax,   order=order )
     call InitGridDef(grid,gridarray(1:4),locked=.true.)
  else
     hires = int_val_opt('-hires',9)
     call InitGridDef(gridarray(3),dy/hires, 0.5_dp, order=order2)
     call InitGridDef(gridarray(2),dy/3.0_dp, 2.0_dp, order=order1)
     call InitGridDef(gridarray(1),dy,        ymax,   order=order )
     call InitGridDef(grid,gridarray(1:3),locked=.true.)
  end if
  ! set parameter for evolution step in Q
  du = dble_val_opt('-du',0.1_dp)
  call SetDefaultEvolutionDu(du)

  ! set up the splitting functions
  nloop = int_val_opt('-nloop',3)
  if (log_val_opt('-exactth')) &
       &   call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact)
  if (log_val_opt('-exactsp')) &
       &   call dglap_Set_nnlo_splitting(nnlo_splitting_exact)

  call cpu_time(time_start)
  call InitDglapHolder(grid, dh, factscheme=factscheme_MSbar, &
       &                              nloop=nloop, nflo=3, nfhi=6)
  call cpu_time(time_init_done)


  ! first way to get the initial distribution
  !call AllocInitPDFSub(grid,vogt_init,VogtInitSub)
  ! alternative way to get the initial distribution
  call AllocPDF(grid,vogt_init)
  vogt_init = unpolarized_dummy_pdf(xValues(grid))

  ! set up the coupling
  Qinit = sqrt(two); Qmax = dble_val_opt('-Qmax',1e4_dp)
  do i = 1, int_val_opt("-nas",1)
     if (i /= 1) call Delete(coupling)
     call InitRunningCoupling(coupling, alfas=0.35_dp, Q=Qinit, &
          &nloop=nloop)
  end do
  write(0,*) Value(coupling,91.2d0)

  ! set up the table
  dlnlnQ = dble_val_opt('-dlnlnQ',0.07_dp)
  olnlnQ = int_val_opt('-olnlnQ',4)
  call AllocPdfTable(grid,table,Qinit,Qmax,dlnlnQ,lnlnQ_order=olnlnQ)
  call AddNfInfoToPdfTable(table,coupling)
  if (preev) call PreEvolvePdfTable(table,Qinit,dh,coupling)
  call cpu_time(time_ev_done)

  ! decide 
  nrep  = int_val_opt('-nrep',1)
  nxQ = int_val_opt('-nxQ',0)
  output = log_val_opt('-output') .or. log_val_opt('-outputgrid')
  outputgrid  = log_val_opt('-outputgrid')
  if (output) call output_info

  !-- security ----------------------
  if (.not. CheckAllArgsUsed(0)) stop
  !----------------------------------

  ! evolution & output
  do i = 1, nrep
     if (preev) then
        call EvolvePdfTable(table,vogt_init)
     else
        call EvolvePdfTable(table,Qinit,vogt_init,dh,coupling)
     end if

     ! one form of output
     if (outputgrid) then
        call eval_output_grid()
     else
        call eval_output_lines()
     end if
     
  end do
  call cpu_time(time_end)
  write(0,'(a,4f10.5)') "Timings (init, preevln, evln) = ", &
       &   time_init_done-time_start, &
       &   time_ev_done-time_init_done, &
       &   (time_end-time_ev_done)/nrep
  write(6,'(a)',advance='no') "# "
  call time_stamp(6)
  call system("echo '# host:' `hostname`")
  ! record info about the cpu
  call system("grep -e name -e cache -e MHz /proc/cpuinfo | sed 's/^/# /'")
  if (output) write(6,'(a,4f10.5)') "# Timings (init, preevln, evln) = ", &
       &   time_init_done-time_start, &
       &   time_ev_done-time_init_done, &
       &   (time_end-time_ev_done)/nrep

  ! clean up
  call Delete(table)
  call Delete(vogt_init)  ! here can also use deallocate(vogt_init)
  call Delete(dh)
  call Delete(coupling)
  call Delete(grid)
  call Delete(gridarray)

contains

  !---------------------------------------------------------------
  !! some basic commentary
  subroutine output_info
    write(6,'(a)') '# '//trim(command_line())
    write(6,'(a)') '# Grid properties:'
    write(6,'(a,4f10.6)')  '# dymax = ', maxval(grid%subgd(:)%dy)
    write(6,'(a,4i10)')    '# dnsty = ', nint(maxval(grid%subgd(:)%dy)/grid%subgd(:)%dy)
    write(6,'(a,4f10.6)')  '# ymax  = ', grid%subgd(:)%ymax
    write(6,'(a,4i10)'  )  '# order = ', grid%subgd(:)%order
    write(6,'(a,4es10.2)') '# eps   = ', grid%subgd(:)%eps
    
    write(6,'(a,f10.5,a,f10.5)') '# Evolution: du = ',du, ",   coupling dt = ", DefaultCouplingDt()
    write(6,'(a,f10.5,a,i5)') '# Tabulation: dlnlnQ = ',dlnlnQ, ', olnlnQ = ',olnlnQ
  end subroutine output_info


  !-------------------------------------------------------------------
  !! output lines in the y,Q plane
  subroutine eval_output_lines()
    integer nn, j
    real(dp) :: y, Q, pdfval(-6:6)
    nn = nxQ/4
    do j = 1, nn
       y = j*ymax/nn
       Q = Qmax - j*(Qmax-Qinit)/nn
       call EvalPdfTable_yQ(table,y,Q,pdfval)
       if (output .and. i==1) write(6,'(20es20.10)') y,Q,pdfval(0:4)
       !if (output .and. i==1) write(6,'(20es20.8)') y,Q,vogt_init(:,0:3).atx.(grid.with.exp(-y))
    end do

    if (output .and. i==1) write(6,*)
    if (output .and. i==1) write(6,*)
    do j = nn,1,-1
       y = j*ymax/nn
       Q = 4.0_dp + j*5.0_dp/nn
       call EvalPdfTable_yQ(table,y,Q,pdfval) 
       if (output .and. i==1) write(6,'(20es20.10)') y,Q,pdfval(0:4)
    end do

    if (output .and. i==1) write(6,*)
    if (output .and. i==1) write(6,*)
    do j = nn,1,-1
       y = j*ymax/nn
       Q = Qmax*(1-j*0.2_dp/nn)
       call EvalPdfTable_yQ(table,y,Q,pdfval) 
       if (output .and. i==1) write(6,'(20es20.10)') y,Q,pdfval(0:4)
    end do

    if (output .and. i==1) write(6,*)
    if (output .and. i==1) write(6,*)
    do j = nn,1,-1
       y = j*ymax/nn
       Q = sqrt(Qinit*Qmax)*(1+j*0.2_dp/nn)
       call EvalPdfTable_yQ(table,y,Q,pdfval) 
       if (output .and. i==1) write(6,'(20es20.10)') y,Q,pdfval(0:4)
    end do
  end subroutine eval_output_lines


  !----------------------------------------------------------------
  !! output the results on a grid uniform in 
  !! zeta = ln 1/x + grid_a*(1-x)
  subroutine eval_output_grid()
    integer  :: nz, nQ, iz, iQ
    real(dp) :: zmax, zeta, y, Q, zQ, zQmax, pdfval(-6:6)
    real(dp), parameter :: grid_a = 9, gridQ_a = 3

    nz = nint(sqrt(four*nxQ))
    nQ = nint(sqrt(0.25_dp*nxQ))-1

    zmax = zeta_of_y(ymax, grid_a)
    zQmax = zeta_of_y(log(Qmax/Qinit), gridQ_a)
    do iQ = 0, nQ
       do iz = 1, nz
          zeta = iz * zmax/nz
          y    = y_of_zeta(zeta, grid_a)
          zQ = (iQ+zeta/zmax) * zQmax / nQ
          Q = max(Qinit,min(Qmax,Qinit * exp(y_of_zeta(zQ, gridQ_a))))
          call EvalPdfTable_yQ(table,y,Q,pdfval) 
          if (output .and. i == 1) write(6,'(20es20.10)') y,Q,pdfval(-5:5)
       end do
          if (output .and. i == 1) write(6,'(a)') 
    end do
    
  end subroutine eval_output_grid
  

  
  !-----------------------------------------------------------------
  !! return zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  function zeta_of_y(y, a) result(zeta)
    real(dp), intent(in) :: y, a
    real(dp)             :: zeta
    zeta = y + a*(one - exp(-y))
  end function zeta_of_y

  
  !-----------------------------------------------------------------
  !! return inversion of zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  function y_of_zeta(zeta, a) result(y)
    real(dp), intent(in) :: zeta, a
    real(dp)             :: y, x, diff_from_zero, deriv
    integer             :: iter
    real(dp), parameter :: eps = 1e-12_dp
    integer,  parameter :: maxiter = 100

    ! starting condition (and soln if a = 0)
    y = zeta 
    if (a /= zero) then
       do iter = 0, maxiter
          x = exp(-y);
          diff_from_zero = zeta - y - a*(one-x);
          ! we have found good solution
          if (abs(diff_from_zero) < eps) exit
          deriv = -one  - a*x;
          y = y - diff_from_zero / deriv;
       end do
    end if
    
    if (iter > maxiter) write(0,*) "y_of_zeta reached maxiter"

  end function y_of_zeta
  
end program prec_and_timing
