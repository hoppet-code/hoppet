!! An example program using structure functions up to N3LO. !!
!!
program structure_functions_example
  use hoppet
  use pdfs_for_benchmarks
  use streamlined_interface
  use structure_functions
  use io_utils, int_value => value
  implicit none
  character (len=400) :: filename,dirname
  real(dp) :: Qmax, Qmin, ymax, Q, mc, mb, mt, asQ0, Q0,&
  & muR_Q, dy, dlnlnQ, minQval, maxQval
  integer  :: order, nloop_max = 4
  real(dp) :: xmur_evolv_vals(3) = (/1.0_dp, 0.5_dp, 2.0_dp/)
  integer  :: ixmur_evolv
  integer  :: idev
  logical  :: vfn = .false.

  idev = idev_open_opt("-o")
  write(idev,'(a,a)') "# ", trim(command_line())
  nloop_max = int_val_opt("-nloop", 4); write(idev,'(a,i4)') "# nloop = ", nloop_max
  Q  = dble_val_opt("-Q", 100.0_dp); write(idev,'(a,f10.4)') "# Q = ", Q
  Q0 = dble_val_opt("-Q0", sqrt(two)) ; write(idev,'(a,f22.16)') "# Q0 = ", Q0
  Qmax = Q
  Qmin = Q0
  asQ0 =dble_val_opt("-asQ0", 0.30_dp); write(idev,'(a,f10.4)') "# asQ0 = ", asQ0

  vfn = log_val_opt("-vfn", .true.); write(idev,'(a,l1)') "# vfn = ", vfn

  if (.not. CheckAllArgsUsed(0)) stop "Unused arguments"


  ! Set heavy flavour scheme
  if (vfn) then
    mc = 1.414213563_dp   ! sqrt(2.0_dp) + epsilon
    mb = 4.5_dp
    mt = 175.0_dp
    call hoppetSetPoleMassVFN(mc, mb, mt)
  else
    call hoppetSetFFN(5)
  end if
  
  ! Streamlined initialization
  ! including  parameters for x-grid
  order = -6 ! interpolation order, not perturbative order in alphas!
  ymax  = 12.0_dp
  dy    = dble_val_opt("-dy",0.05_dp)
  dlnlnQ = dy/4.0_dp
  minQval = Qmin
  maxQval = Qmax 

  ! initialise the grid and dglap holder, using the streamlined
  ! interface for simplicity
  call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop_max,&
    &         order,factscheme_MSbar)

  do ixmur_evolv = 1, size(xmur_evolv_vals)
     write(idev,"(a,f10.4)") "# xmuR = ", xmur_evolv_vals(ixmur_evolv)
     call hoppetEvolve(asQ0, Q0, nloop_max, xmur_evolv_vals(ixmur_evolv), benchmark_pdf_unpolarized_lha, Q0)
     call write_pdf(Q, 100)
  end do



contains
  
  ! write the time since the last call to this routine (or since start of the program)
  subroutine write_time()
    implicit none
    real(dp), save :: time, last_time = 0.0_dp
    call cpu_time(time)
    write(6,'(a,f10.4)') "Time since last call: ", time - last_time
    last_time = time
  end subroutine write_time

  !> @brief write the F1 structure function to idev
  !!
  !! Writes the W and Z F1 structure functions to a file with equal spacing in log(x)
  !!
  !! @param[in]     idev     file/device to write to
  !! @param[in]     Qtest    Value of Q to use
  !! @param[in]     ymax     Largest value of y = - log(x) to print
  !! @param[in]     ny       number of bins in y to print
  !!
  subroutine write_pdf(Qtest, ny)
    real(dp), intent(in) :: Qtest
    integer, intent(in)  :: ny
    real(dp) :: ytest, xval, xpdf_at_xQ(-6:6)
    integer  :: iy
    !F1 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    write(idev,'(a)') '# x pdf(-6:6)'
    
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
       call hoppetEval(xval,Qtest,xpdf_at_xQ)
       write(idev,'(14es16.8)') xval, xpdf_at_xQ
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_pdf
  

  
end program structure_functions_example


