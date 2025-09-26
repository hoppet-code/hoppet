!======================================================================
!
! Program for preparing info for accuracy tests and for carrying out
! timing tests.
!
! A high-accuracy reference grid as used in the documentation can be
! prepared with the instructions
!
!  ./prec_and_timing -outputgrid -nxQ 5000 -nopreev -nloop 3 -dy 0.025 -o refgrid.dat
!
! Fig. 2 and fig.3(left) of the doc can be obtained by comparing
! this with the output from 
!
!   ./prec_and_timing -nrep 1 -nxQ 5000 -outputgrid -dy 0.2 -nopreev -o testgrid.dat
!
! To compare the outputs, use the command
!
!   ./compare2files_v2 testgrid.dat refgrid.dat -channel 11 [-protect]
!
! Other options of interest in prec_and_timing include
!   -eps [1e-7]       precision of adaptive integration for preparing split.fn.
!   -asdt [0.2]       step for evolution of alpha_s
!   -exact-nnlo-sp    use exact 3-loop split-fns
!   -exact-nnlo-th    use exact 3-loop mass thresholds
!
! BEWARE: older copies of executables have ymax=5 and Qmax=100(?). Explicitly
!         specify -ymax 11.5 and -Qmax 1e4 to get sensible timings.
!
! NB: hoppet v1 used the following command for the test grid
!
!   ./prec_and_timing -nrep 1 -nxQ 5000 -outputgrid -dy 0.2 -order -6 \
!                     -nopreev -4grids -dlnlnQ 0.05 -du 0.4 -olnlnQ 4 -o testgrid.dat
!
! LHAPDF output
! -------------
! This code can also produce an LHAPDF output file, for example, assuming
! your LHAPDF sets are kept in ~/LHAPDF, then the following will produce a dummy
! set corresponding to the standard Les Houches PDF benchmark
!
! ./benchmarking/prec_and_timing  -nloop 3 -dy 0.05 -lhapdf-out ~/LHAPDF/dummy/dummy -o ~/LHAPDF/dummy/hoppet.log
!
program prec_and_timing
  use hoppet
  use hoppet_git_state
  use pdfs_for_benchmarks
  use io_utils
  use NameSelect

  implicit none
  !-------------------------------
  type(grid_def)     :: grid, gridarray(4)
  type(dglap_holder) :: dh
  type(running_coupling) :: coupling
  real(dp), allocatable  :: yvals(:), Qvals(:)
  integer            :: order, order1, order2, nloop, i, j, nrep, nrep_orig, nxQ,olnlnQ
  integer            :: nrep_eval, n_alphas, ny, nQ, nrep_interp
  integer            :: hires
  real(dp)           :: dy, Qinit, Qmin, Qmax, du, dlnlnQ
  real(dp)           :: ymax, pdfval(-6:6)
  real(dp)           :: time_start, time_init_done, time_prev_done, time_ev_done, time_interp_start, time_end
  real(dp), pointer  :: initial_condition(:,:)
  type(pdf_table)    :: table
  logical            :: output, outputgrid, output_benchmark, preev, auto_nrep
  logical            :: write_lhapdf
  character(len=999) :: lhapdf_out = ""
  integer            :: idev, y_interp_order
  integer, parameter :: lhapdf_unit = 99
  integer            :: lhapdf_iyinc
  character(len=300) :: hostname

  idev = idev_open_opt("-o","/dev/stdout")    ! output file (required)

  ! set the details of the y=ln1/x grid
  dy    = dble_val_opt('-dy',0.10_dp)         ! grid spacing in y = ln(1/x)
  ymax  = dble_val_opt('-ymax',12.0_dp)       ! maximum y = ln(1/x)
  order = int_val_opt('-order',-6)            ! interpolation order for splitting-function rep
  dlnlnQ = dble_val_opt('-dlnlnQ',dy/4.0_dp)  ! table spacing in ln(ln(Q))
  olnlnQ = int_val_opt('-olnlnQ',4)           ! table interpolation order for ln(ln(Q))
  y_interp_order = int_val_opt('-yinterp-order',-1) ! table interpolation order for y
  Qmax = dble_val_opt('-Qmax',1e4_dp)

  du = dble_val_opt('-du',0.1_dp) ! evolution step, overriden by dlnlnQ

  order2 = int_val_opt('-order2',order) ! override interp order for finest grid
  order1 = int_val_opt('-order1',order) ! override interp order for finest second grid

  preev = log_val_opt('-preev',.true.)        ! use pre-evolution (cached evolution)

  lhapdf_out = string_val_opt("-lhapdf-out","") ! write out an LHAPDF file
  lhapdf_iyinc = int_val_opt("-lhapdf-iyinc",1)

  if (log_val_opt('-eps')) call SetDefaultConvolutionEps(dble_val_opt('-eps')) ! precision for split.fn. adaptive integration 
  if (log_val_opt('-asdt')) call SetDefaultCouplingDt(dble_val_opt('-asdt')) ! spacing for coupling ev.

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
  else if (log_val_opt('-4grids-wider')) then
    call InitGridDef(gridarray(4),dy/27.0_dp, 0.24_dp, order=order2)
    call InitGridDef(gridarray(3),dy/9.0_dp, 0.6_dp, order=order1)
    call InitGridDef(gridarray(2),dy/3.0_dp, 2.4_dp, order=order )
    call InitGridDef(gridarray(1),dy,        ymax,   order=order )
    call InitGridDef(grid,gridarray(1:4),locked=.true.)
  else if (log_val_opt('-3grids')) then
    hires = int_val_opt('-hires',9)
    call InitGridDef(gridarray(3),dy/hires, 0.5_dp, order=order2)
    call InitGridDef(gridarray(2),dy/3.0_dp, 2.0_dp, order=order1)
    call InitGridDef(gridarray(1),dy,        ymax,   order=order )
    call InitGridDef(grid,gridarray(1:3),locked=.true.)
  else 
    if (order1 /= order .or. order2 /= order) then
       call wae_error("order1 and order2 != order not suppoed with default grid")
    end if
    ! this is like 4 grids and is the standard, without control over order1 and order2
    call InitGridDefDefault(grid, dy, ymax, order=order)
  end if
  ! set parameter for evolution step in Q
  call SetDefaultEvolutionDu(du)


  ! set up the splitting functions
  nloop = int_val_opt('-nloop',3)
  if (log_val_opt('-exact-nnlo-th')) &
       &   call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact)
  if (log_val_opt('-exact-nnlo-sp')) &
       &   call dglap_Set_nnlo_splitting(nnlo_splitting_exact)
 

  ! decide how many repetitions, and what form of output
  n_alphas = int_val_opt("-n-alphas",1)
  nrep  = int_val_opt('-nrep',1)
  auto_nrep = log_val_opt('-auto-nrep',.false.)
  nrep_eval = int_val_opt('-nrep-eval',10)
  nrep_interp = int_val_opt('-nrep-interp',1000)
  nxQ = int_val_opt('-nxQ',0)
  output = log_val_opt('-output') .or. log_val_opt('-outputgrid')
  outputgrid  = log_val_opt('-outputgrid')
  output_benchmark = log_val_opt('-output-benchmark')
  !-- security ----------------------
  if (.not. CheckAllArgsUsed(0)) error stop
  !----------------------------------

  if (output .or. trim(lhapdf_out) /= "") call output_info
  write(idev,'(a)',advance='no') "# "
  call time_stamp(idev)
  call get_hostname(hostname)  
  write(idev,'(a)') "# host: "//trim(hostname)
  call hoppet_print_git_state(idev,prefix="#")

  call cpu_time(time_start)
  call InitDglapHolder(grid, dh, factscheme=factscheme_MSbar, &
       &                              nloop=nloop, nflo=3, nfhi=6)
  call cpu_time(time_init_done)


  ! first way to get the initial distribution
  !call AllocInitPDFSub(grid,vogt_init,VogtInitSub)
  ! alternative way to get the initial distribution
  call AllocPDF(grid,initial_condition)
  initial_condition = benchmark_pdf_unpolarized(xValues(grid))

  ! set up the coupling
  Qinit = sqrt(two)
  Qmin = one
  do i = 1, n_alphas 
     if (i /= 1) call Delete(coupling)
     call InitRunningCoupling(coupling, alfas=0.35_dp, Q=Qinit, &
          &nloop=nloop)
  end do
  !AV commented out because it is an error in non-GNU compilers 
  !write(0,*) Value(coupling,91.2d0)

  ! set up the table
  call AllocPdfTable(grid,table,Qmin,Qmax,dlnlnQ,lnlnQ_order=olnlnQ)
  call AddNfInfoToPdfTable(table,coupling)
  if (preev) call PreEvolvePdfTable(table,Qinit,dh,coupling)
  call cpu_time(time_prev_done)

  !call PDFTableSetYInterpOrder(y_interp_order)
  if (y_interp_order > 0 .and. olnlnQ >0) then
    call PdfTableOverrideInterpOrders(y_interp_order,olnlnQ)
  end if

  ! record info about the cpu
  !call system("grep -e name -e cache -e MHz /proc/cpuinfo | sed 's/^/# /'")

  ! evolution & output
  do i = 1, nrep
     if (preev) then
        call EvolvePdfTable(table,initial_condition)
     else
        call EvolvePdfTable(table,Qinit,initial_condition,dh,coupling)
     end if     
  end do
  call cpu_time(time_ev_done)

  ! auto-timing refinement of nrep
  if (auto_nrep .and. time_ev_done-time_prev_done < 0.05_dp) then
     nrep_orig = nrep
     nrep = 0.1_dp / ((time_ev_done-time_prev_done)/nrep_orig)
     do i = nrep_orig+1, nrep
        if (preev) then
           call EvolvePdfTable(table,initial_condition)
        else
           call EvolvePdfTable(table,Qinit,initial_condition,dh,coupling)
        end if     
     end do
     call cpu_time(time_ev_done)
  end if

  ! Here we compute a list of y and Q values to use when interpolating
  ! the grid. We use ceiling(sqrt(nrep)) points in each direction.
  ny = nint(sqrt(four*nrep_interp))
  nQ = nint(sqrt(0.25_dp*nrep_interp))-1
  if (ny > 0) then
    ALLOCATE(yvals(ny))
    do i = 1, ny
       yvals(i) = i*ymax/ny
    end do
  end if 
  if (nQ > 0) then
    ALLOCATE(Qvals(nQ))
    do i = 1, nQ
       Qvals(i) = Qinit + i*(Qmax-Qinit)/nQ
    end do
  end if

  call cpu_time(time_interp_start)
  do i = 1, ny
     do j = 1, nQ
        call EvalPdfTable_yQ(table,yvals(i),Qvals(j),pdfval)
     end do
  end do
  call cpu_time(time_end)

  if (output_benchmark) then
    call write_benchmark_output()
  end if

  ! one form of output
  if (outputgrid) then
     call eval_output_grid()
  else
     call eval_output_lines()
  end if

  write(6,'(a,4f16.11," s, nrep=",i7, ", nrep_interp=", i7)') "Timings (init, preevln, evln, interp) = ", &
       &   time_init_done-time_start, &
       &   time_prev_done-time_init_done, &
       &   (time_ev_done-time_prev_done)/nrep, &
       &   (time_end - time_interp_start)/ny/nQ, nrep , nrep_interp
  if (output) write(idev,'(a,4f16.11," s, nrep=",i7, ", nrep_interp=", i7)') "# Init Timings (init, preevln, evln, interp) = ", &
       &   time_init_done-time_start, &
       &   time_prev_done-time_init_done, &
       &   (time_ev_done-time_prev_done)/nrep,  &
       &   (time_end - time_interp_start)/ny/nQ, nrep, nrep_interp
  write(idev,'(a,f10.5," s")') "# Initialisation time = ", time_init_done-time_start
  write(idev,'(a,f10.5," s")') "# Pre-evolution time = ", time_prev_done-time_init_done
  write(idev,'(a,f10.5," s, nrep = ",i7)') "# Evolution time = ", (time_ev_done-time_prev_done)/nrep, nrep
  write(idev,'(a,f10.5," ns, nrep_interp =", i7)') "# Interpolation time = ", (time_end - time_interp_start)/ny/nQ*1d9, nrep_interp

  !call get_evaluation_times()
  call get_evaluation_times_new()

  if (trim(lhapdf_out) /= "") then
    call WriteLHAPDFFromPdfTable(table, coupling, lhapdf_out, &
           pdf_index = 0, iy_increment = lhapdf_iyinc)
  end if


  ! clean up
  call Delete(table)
  call Delete(initial_condition)  ! here can also use deallocate(vogt_init)
  call Delete(dh)
  call Delete(coupling)
  call Delete(grid)
  call Delete(gridarray)

contains

  !---------------------------------------------------------------
  !! some basic commentary
  subroutine output_info
    write(idev,'(a)') '# '//trim(command_line())
    write(idev,'(a)') '# Grid properties:'
    write(idev,'(a,4f10.6)')  '# dymax = ', maxval(grid%subgd(:)%dy)
    write(idev,'(a,4i10)')    '# dnsty = ', nint(maxval(grid%subgd(:)%dy)/grid%subgd(:)%dy)
    write(idev,'(a,4f10.6)')  '# ymax  = ', grid%subgd(:)%ymax
    write(idev,'(a,4i10)'  )  '# order = ', grid%subgd(:)%order
    write(idev,'(a,4es10.2)') '# eps   = ', grid%subgd(:)%eps
    
    write(idev,'(a,f10.5,a,f10.5)') '# Evolution: du = ',du, ",   coupling dt = ", DefaultCouplingDt()
    write(idev,'(a,f10.5,a,i5)') '# Tabulation: dlnlnQ = ',dlnlnQ, ', olnlnQ = ',olnlnQ

    if (nloop >= 3) then
       write(idev,'(a)') '# nnlo_spltting = '//trim(NameOfCode(nnlo_splitting_variant,'nnlo_splitting'))
    else   
       write(idev,'(a)') '# nnlo_spltting = N/A'
    end if
    if (nloop >= 4) then
       !write(idev,'(a)') '# n3lo_spltting = '//trim(NameOfCode(n3lo_splitting_variant,'n3lo_splitting'))
      write(idev,'(a)') '# n3lo_spltting = '// &
      &       trim(NameOfCode(n3lo_splitting_variant,'n3lo_splitting'))
      write(idev,'(a)') '# n3lo_spltting_approx = '// &
      &       trim(NameOfCode(n3lo_splitting_approximation,'n3lo_splitting_approx'))
    else 
      write(idev,'(a)') '# n3lo_spltting = N/A'
      write(idev,'(a)') '# n3lo_spltting_approx = N/A'
    end if
  end subroutine output_info


  !-------------------------------------------------------------------
  !! output 4 lines in the y,Q plane, each with Q relating to y in a different
  !! way
  subroutine eval_output_lines()
    integer nn, j
    real(dp) :: y, Q, pdfval(-6:6)
    character(len=*), parameter :: header = "# y=ln1/x Q xf(x,0:4)"
    if (.not. output) return

    nn = nxQ/4
    write(idev,'(a)') "# Q = Qmax - j*(Qmax-Qinit)/nn"
    write(idev,'(a)') header
    do j = 1, nn
       y = j*ymax/nn
       Q = Qmax - j*(Qmax-Qinit)/nn
       call EvalPdfTable_yQ(table,y,Q,pdfval)
       write(idev,'(20es20.10)') y,Q,pdfval(0:4)
       !if (output .and. i==1) write(idev,'(20es20.8)') y,Q,vogt_init(:,0:3).atx.(grid.with.exp(-y))
    end do

    write(idev,*)
    write(idev,*)
    write(idev,'(a)') "# Q = 4.0_dp + j*5.0_dp/nn"
    write(idev,'(a)') header
    do j = nn,1,-1
       y = j*ymax/nn
       Q = 4.0_dp + j*5.0_dp/nn
       call EvalPdfTable_yQ(table,y,Q,pdfval) 
       write(idev,'(20es20.10)') y,Q,pdfval(0:4)
    end do

    write(idev,*)
    write(idev,*)
    write(idev,'(a)') "# Q = Qmax*(1-j*0.2_dp/nn)"
    write(idev,'(a)') header
    do j = nn,1,-1
       y = j*ymax/nn
       Q = Qmax*(1-j*0.2_dp/nn)
       call EvalPdfTable_yQ(table,y,Q,pdfval) 
       write(idev,'(20es20.10)') y,Q,pdfval(0:4)
    end do

    write(idev,*)
    write(idev,*)
    write(idev,'(a)') "# Q = sqrt(Qinit*Qmax)*(1+j*0.2_dp/nn)"
    write(idev,'(a)') header
    do j = nn,1,-1
       y = j*ymax/nn
       Q = sqrt(Qinit*Qmax)*(1+j*0.2_dp/nn)
       call EvalPdfTable_yQ(table,y,Q,pdfval) 
       write(idev,'(20es20.10)') y,Q,pdfval(0:4)
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
    write(6,*) 'nQ = ', nQ, ' nz = ', nz
    write(idev,'(a)') "# y=ln1/x Q pdf(-5:5)"
    do iQ = 0, nQ
       do iz = 1, nz
          zeta = iz * zmax/nz
          y    = y_of_zeta(zeta, grid_a)
          zQ = (iQ+zeta/zmax) * zQmax / nQ
          Q = max(Qinit,min(Qmax,Qinit * exp(y_of_zeta(zQ, gridQ_a))))
          call EvalPdfTable_yQ(table,y,Q,pdfval) 
          write(idev,'(20es20.10)') y,Q,pdfval(-5:5)
       end do
       write(idev,'(a)') 
    end do
    
  end subroutine eval_output_grid
  

  
  !-----------------------------------------------------------------
  !! return zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  pure function zeta_of_y(y, a) result(zeta)
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
  

  !-----------------------------------------------------------------
  !! new (2025) routine for getting evaluation times, which uses a table
  !! of y,Q values, to reduce the likelihood of being biased by some
  !! special case.
  subroutine get_evaluation_times_new()
    real(dp) :: pdf_g, pdf_g_sum, pdf_all(-6:6), pdf_sum(-6:6)
    integer, parameter :: n = 100
    real(dp) :: yvals(n), Qvals(n), y, Q
    integer :: i, irep, iQ, iy, nvals
    character(len=*), parameter :: fmt = '(a,f6.1," ",a)'

    ! set up the y and Q values
    do i = 1, n
      yvals(i) = i*ymax/n
      Qvals(i) = Qinit + i*(Qmax-Qinit)/n
    end do

    ! then loop over 
    pdf_sum = zero
    call cpu_time(time_start)
    do irep = 1, nrep_eval
      do iQ = 1, n
        Q = Qvals(iQ)
        do iy = 1, n
          y = yvals(iy)
          call EvalPdfTable_yQ(table,y,Q,pdf_all)
          pdf_sum = pdf_sum + pdf_all ! prevent compiler from optimising this out
        end do
      end do
    end do
    call cpu_time(time_end)
    write(idev,fmt) "# Evaluation (all flav)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"
    write(0   ,fmt) "# Evaluation (all flav)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"

    ! then loop over 
    pdf_sum = zero
    call cpu_time(time_start)
    do irep = 1, nrep_eval
      do iQ = 1, n
        Q = Qvals(iQ)
        do iy = 1, n
          y = yvals(iy)
          call EvalPdfTable_yQ_order22(table,y,Q,pdf_all)
          pdf_sum = pdf_sum + pdf_all ! prevent compiler from optimising this out
        end do
      end do
    end do
    call cpu_time(time_end)
    write(idev,fmt) "# Evaluation (all flav,o2)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"
    write(0   ,fmt) "# Evaluation (all flav,o2)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"
    
    ! then loop over 
    pdf_sum = zero
    call cpu_time(time_start)
    do irep = 1, nrep_eval
      do iQ = 1, n
        Q = Qvals(iQ)
        do iy = 1, n
          y = yvals(iy)
          call EvalPdfTable_yQ_order33(table,y,Q,pdf_all)
          pdf_sum = pdf_sum + pdf_all ! prevent compiler from optimising this out
        end do
      end do
    end do
    call cpu_time(time_end)
    write(idev,fmt) "# Evaluation (all flav,o3)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"
    write(0   ,fmt) "# Evaluation (all flav,o3)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"
    
    ! then loop over 
    pdf_sum = zero
    call cpu_time(time_start)
    do irep = 1, nrep_eval
      do iQ = 1, n
        Q = Qvals(iQ)
        do iy = 1, n
          y = yvals(iy)
          call EvalPdfTable_yQ_order44(table,y,Q,pdf_all)
          pdf_sum = pdf_sum + pdf_all ! prevent compiler from optimising this out
        end do
      end do
    end do
    call cpu_time(time_end)
    write(idev,fmt) "# Evaluation (all flav,o4)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"
    write(0   ,fmt) "# Evaluation (all flav,o4)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"


    call cpu_time(time_start)
    pdf_g_sum = zero
    do irep = 1, nrep_eval
      do iQ = 1, n
        Q = Qvals(iQ)
        do iy = 1, n
          y = yvals(iy)
          pdf_g_sum = pdf_g_sum + EvalPdfTable_yQf(table,y,Q,0)
          !call EvalPdfTable_yQ(table,y,Q,pdf_all)
          !pdf_sum = pdf_sum + pdf_all ! prevent compiler from optimising this out
        end do
      end do
    end do
    !do i = 1, nrep_eval
    !  pdf_g_sum = pdf_g_sum + EvalPdfTable_yQf(table,y,Q,0)
    !end do  
    call cpu_time(time_end)
    write(idev,fmt) "# Evaluation (one flav)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"
    write(0   ,fmt) "# Evaluation (one flav)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"


    call cpu_time(time_start)
    pdf_g_sum = zero
    do irep = 1, nrep_eval
      do iQ = 1, n
        Q = Qvals(iQ)
        do iy = 1, n
          y = yvals(iy)
          pdf_g_sum = pdf_g_sum + EvalPdfTable_yQf_order22(table,y,Q,0)
          !call EvalPdfTable_yQ(table,y,Q,pdf_all)
          !pdf_sum = pdf_sum + pdf_all ! prevent compiler from optimising this out
        end do
      end do
    end do
    call cpu_time(time_end)
    write(idev,fmt) "# Evaluation (1fl,ord2)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"
    write(0   ,fmt) "# Evaluation (1fl,ord2)", (time_end-time_start)/(nrep_eval*n**2)*1e9_dp, "ns"

  end subroutine get_evaluation_times_new


  ! old routine for getting evaluation times, which uses just a single y,Q combination
  subroutine get_evaluation_times()
    real(dp) :: pdf_g, pdf_g_sum, pdf_all(-6:6), pdf_sum(-6:6)
    real(dp) :: y = 0.5d0, Q = 10.0d0
    character(len=*), parameter :: fmt = '(a,f6.1," ",a)'
    call cpu_time(time_start)
    pdf_sum = zero
    do i = 1, nrep_eval
      call EvalPdfTable_yQ(table,y,Q,pdf_all)
      pdf_sum = pdf_sum + pdf_all
    end do  
    call cpu_time(time_end)
    write(idev,fmt) "# Evaluation (all flav)", (time_end-time_start)/nrep_eval*1e9_dp, "ns"
    write(0   ,fmt) "# Evaluation (all flav)", (time_end-time_start)/nrep_eval*1e9_dp, "ns"

    call cpu_time(time_start)
    pdf_g_sum = zero
    do i = 1, nrep_eval
      pdf_g_sum = pdf_g_sum + EvalPdfTable_yQf(table,y,Q,0)
    end do  
    call cpu_time(time_end)
    write(idev,fmt) "# Evaluation (one flav)", (time_end-time_start)/nrep_eval*1e9_dp, "ns"
    write(0   ,fmt) "# Evaluation (one flav)", (time_end-time_start)/nrep_eval*1e9_dp, "ns"

  end subroutine get_evaluation_times

  subroutine write_benchmark_output

    real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
    integer  :: ix    ! get the value of the tabulation at some point
    real(dp) :: Q, pdf_at_xQ(-6:6)
    Q = 100.0_dp
    write(6,'(a)')
    write(6,'(a,f8.3,a)') "           Evaluating PDFs at Q = ",Q," GeV"
    write(6,'(a5,2a12,a14,a10,a12)') "x",&
         & "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
    do ix = 1, size(heralhc_xvals)
      call EvalPdfTable_xQ(table,heralhc_xvals(ix),Q,pdf_at_xQ)
      write(6,'(es7.1,5es12.4)') heralhc_xvals(ix), &
            &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
            &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
            &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
            &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
            &  pdf_at_xQ(0)
    end do
  end subroutine write_benchmark_output
end program prec_and_timing
