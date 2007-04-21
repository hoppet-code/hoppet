
! on tycho the following combination runs in 4ms per pdf and has
! a relative accuracy (on g,d,u,s) that is typically \lesssim 1e-4, though
! one might want to check the accuracy more systematically.
! 
! ./small_fast_tab -nrep 1 -nxQ 500 -output -dy 0.25 -dlnlnQ 0.07 -dt 0.4 -olnlnQ 4 -order 6 -order2 -6
!
! roughly speaking: evolution takes about 2.5ms and evaluations 1.5ms
!
! Funny thing is that earlier it seemed that runs were much faster
! (factor of two for evolution) and I can't figure out what has
! changed...

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
program small_fast_tab
  use hoppet_v1
  use pdf_initial_condition
  use sub_defs_io
  implicit none
  !-------------------------------
  type(grid_def)     :: grid, gridarray(3)
  type(dglap_holder) :: dh
  type(running_coupling)    :: coupling
  integer            :: order, order2, nloop, i, j, nrep, nxQ,nn, olnlnQ
  real(dp)           :: dy, Qinit, Qmax, y, Q, pdfval(-6:6), dt, dlnlnQ
  real(dp)           :: ymax
  real(dp), pointer  :: vogt_init(:,:)
  type(pdf_table)       :: table
  logical            :: output

  ! set the details of the y=ln1/x grid
  dy    = dble_val_opt('-dy',0.25_dp)
  ymax  = dble_val_opt('-ymax',5.0_dp)
  order = int_val_opt('-order',5)
  order2 = int_val_opt('-order2',order)
  if (log_val_opt('-alt')) then
     call InitGridDef(gridarray(3),dy/7.0_dp, 1.0_dp, order=order2)
     call InitGridDef(gridarray(2),dy/2.0_dp, 3.0_dp, order=order2)
     call InitGridDef(gridarray(1),dy,        ymax, order=order)
  else
     call InitGridDef(gridarray(3),dy/9.0_dp, 0.5_dp, order=order2)
     call InitGridDef(gridarray(2),dy/3.0_dp, 2.0_dp, order=order2)
     call InitGridDef(gridarray(1),dy,        ymax,   order=order )
  end if
  call InitGridDef(grid,gridarray(1:3),locked=.true.)
  ! set parameter for evolution step in Q
  dt = dble_val_opt('-dt',0.4_dp)
  call SetDefaultEvolutionDt(dt)

  ! set up the splitting functions
  nloop = int_val_opt('-nloop',3)
  if (log_val_opt('-exactth')) &
       &   call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact)
  if (log_val_opt('-exactsp')) &
       &   call dglap_Set_nnlo_splitting(nnlo_splitting_exact)
  call InitDglapHolder(grid, dh, factscheme=factscheme_MSbar, &
       &                              nloop=nloop, nflo=3, nfhi=6)
  
  ! first way to get the initial distribution
  !call AllocInitPDFSub(grid,vogt_init,VogtInitSub)
  ! alternative way to get the initial distribution
  call AllocPDF(grid,vogt_init)
  vogt_init = unpolarized_dummy_pdf(xValues(grid))

  ! set up the coupling
  Qinit = sqrt(two); Qmax = dble_val_opt('-Qmax',50.0_dp)
  call InitRunningCoupling(coupling, alfas=0.35_dp, Q=Qinit, &
       &nloop=nloop, use_nah=.true.)

  ! set up the table
  dlnlnQ = dble_val_opt('-dlnlnQ',0.07_dp)
  olnlnQ = int_val_opt('-olnlnQ',4)
  call AllocPdfTable(grid,table,Qinit,Qmax,dlnlnQ,lnlnQ_order=olnlnQ)
  call AddNfInfoToPdfTable(table,coupling)
  call PreEvolvePdfTable(table,Qinit,dh,coupling)

  nrep  = int_val_opt('-nrep',1)
  nxQ = int_val_opt('-nxQ',0); y = 3.14_dp; Q = 13.354_dp
  output = log_val_opt('-output')
  if (output) call output_info
 
  ! output
  do i = 1, nrep
     call EvolvePdfTable(table,vogt_init)
     
     nn = nxQ/3
     do j = 1, nn
        y = j*5.0_dp/nn
        Q = Qmax - j*(Qmax-Qinit)/nn
        call EvalPdfTable_yQ(table,y,Q,pdfval)
        if (output .and. i==1) write(6,'(20es20.8)') y,Q,pdfval(0:3)
        !if (output .and. i==1) write(6,'(20es20.8)') y,Q,vogt_init(:,0:3).atx.(grid.with.exp(-y))
     end do

     if (output .and. i==1) write(6,*)
     if (output .and. i==1) write(6,*)
     do j = nn,1,-1
        y = j*5.0_dp/nn
        Q = 4.0_dp + j*5.0_dp/nn
        call EvalPdfTable_yQ(table,y,Q,pdfval) 
        if (output .and. i==1) write(6,'(20es20.8)') y,Q,pdfval(0:3)
     end do
     
     if (output .and. i==1) write(6,*)
     if (output .and. i==1) write(6,*)
     do j = nn,1,-1
        y = j*5.0_dp/nn
        Q = Qmax*(1-j*0.2_dp/nn)
        call EvalPdfTable_yQ(table,y,Q,pdfval) 
        if (output .and. i==1) write(6,'(20es20.8)') y,Q,pdfval(0:3)
     end do
  end do
  
  ! clean up
  call Delete(table)
  call Delete(vogt_init)  ! here can also use deallocate(vogt_init)
  call Delete(dh)
  call Delete(coupling)
  call Delete(grid)
  call Delete(gridarray)

contains
  subroutine output_info
    write(6,'(a)') '# '//trim(command_line())
    write(6,'(a,f10.5,a,f10.5)') '# dy = ',dy,      ';    ymax = ',ymax
    write(6,'(a,f10.5,a,f10.5)') '# dt = ',dt,      ';    Qmax = ',Qmax
    write(6,'(a,i5,a,i5)')    '# order  = ',order,  ';    order2 = ',order2
    write(6,'(a,f10.5,a,i5)') '# dlnlnQ = ',dlnlnQ, ';    olnlnQ = ',olnlnQ
  end subroutine output_info
  
end program small_fast_tab
