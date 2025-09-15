!! An example program using structure functions up to N3LO. !!
!!
program structure_functions_example
  use hoppet
  use pdfs_for_benchmarks
  use streamlined_interface
  use structure_functions
  implicit none
  character (len=400) :: filename,dirname
  real(dp) :: Qmax, Qmin, ymax, Q, mc, mb, mt, asQ, Q0,&
  & muR_Q, dy, dlnlnQ, minQval, maxQval
  integer  :: order, nloop_max = 4

  Qmax = 13000.0_dp 
  Qmin = one
  
  ! Set heavy flavour scheme
  !mc = 1.414213563_dp   ! sqrt(2.0_dp) + epsilon
  !mb = 4.5_dp
  !mt = 175.0_dp
  !call hoppetSetPoleMassVFN(mc, mb, mt)
  call hoppetSetFFN(5)
!  call hoppetSetApproximateDGLAPN3LO(n3lo_splitting_approximation_up_to_2310_05744)
!  call hoppetSetApproximateDGLAPN3LO(n3lo_splitting_approximation_up_to_2404_09701)
  call hoppetSetApproximateDGLAPN3LO(n3lo_splitting_approximation_up_to_2410_08089)
  
  ! Streamlined initialization
  ! including  parameters for x-grid
  order = -6 ! interpolation order, not perturbative order in alphas!
  ymax  = 12.0_dp
  dy    = 0.05_dp  ! dble_val_opt("-dy",0.1_dp)
  dlnlnQ = dy/4.0_dp
  minQval = Qmin
  maxQval = Qmax 

  ! initialise the grid and dglap holder, using the streamlined
  ! interface for simplicity
  call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop_max,&
    &         order,factscheme_MSbar)

  ! Setup all constants and parameters needed by the structure functions
  write(6,'(a)') "Filling structure function coefficient tables"
   call StartStrFct(nloop_max, nflav = 5, scale_choice =&
        & scale_choice_Q, param_coefs = .true.)


  dirname = 'structure_functions_example_n3lo_tests_results/'

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-1.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,4,1.0_dp,1.0_dp,1.0_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-0.5.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,4,1.0_dp,1.0_dp,0.5_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-n3lo-evol-xmur-2.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,4,1.0_dp,1.0_dp,2.0_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,3,1.0_dp,1.0_dp,1.0_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-0.5.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,3,1.0_dp,1.0_dp,0.5_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-2.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,3,1.0_dp,1.0_dp,2.0_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-nlo-evol-xmur-1.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,2,1.0_dp,1.0_dp,1.0_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-nlo-evol-xmur-2.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,2,1.0_dp,1.0_dp,2.0_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-1.0-xmuf-1.0-nlo-evol-xmur-0.5.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,2,1.0_dp,1.0_dp,0.5_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-2.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,3,2.0_dp,1.0_dp,1.0_dp,99)

  filename = trim(dirname)//'n3lo-coef-xmur-0.5-xmuf-1.0-nnlo-evol-xmur-1.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(4,3,0.5_dp,1.0_dp,1.0_dp,99)

  filename = trim(dirname)//'nnlo-coef-xmur-1.0-xmuf-1.0-nnlo-evol-xmur-1.0.dat'
  open(unit = 99, file = trim(filename))
  call start_evolve_init_write(3,3,1.0_dp,1.0_dp,1.0_dp,99)

contains

  subroutine start_evolve_init_write(nloop_coefs,nloop_evolv,xmur,xmuf,xmur_evolv,idev)
    implicit none
    integer,  intent(in) :: nloop_coefs, nloop_evolv, idev
    real(dp), intent(in) :: xmur, xmuf, xmur_evolv
    
    ! Set up all constants and parameters needed by the structure functions
    ! call StartStrFct(nloop_coefs, nflav = 5, xR = xmur, xF = xmuf, scale_choice = scale_choice_Q, &
    !   param_coefs = .true.)

    ! Evolve the PDF
    ! asQ = 0.35_dp
    ! Q0 = sqrt(2.0_dp)
    asQ = 0.30_dp
    Q0 = 2.0_dp
    
    call write_time()
    write(6,'(a)') "Doing evolution to fill PDF table"
    call hoppetEvolve(asQ, Q0, nloop_evolv, xmur_evolv, benchmark_pdf_unpolarized_lha, Q0)
    ! NB: this uses the PDFs that were set up in the streamlined interface
    ! with the hoppetEvolve routine
    call write_time()
    write(6,'(a)') "Filling StrFct tables"
    call InitStrFct(nloop_coefs, .true., xmur, xmuf)
    
    ymax = log(1e5) !ymax=20
    Q = 100.0_dp
    !Q = Q0

    call write_time()
    call write_f(idev, Q, ymax, 100, xmur, xmuf)
    
  end subroutine start_evolve_init_write
  
  ! write the time since the last call to this routine (or since start of the program)
  subroutine write_time()
    implicit none
    real(dp), save :: time, last_time = 0.0_dp
    call cpu_time(time)
    write(6,'(a,f10.4)') "Time since last call: ", time - last_time
    last_time = time
  end subroutine write_time

  subroutine write_f(idev, Qtest, ymax, ny, xmur, xmuf)
    implicit none
    real(dp), intent(in) :: Qtest, ymax
    integer,  intent(in) :: idev, ny
    real(dp), intent(in) :: xmur, xmuf
    
    call write_f1 (idev, Qtest, ymax, ny, muF = xmuf*Qtest, muR = xmur*Qtest)
    call write_f2 (idev, Qtest, ymax, ny, muF = xmuf*Qtest, muR = xmur*Qtest)
    call write_f3 (idev, Qtest, ymax, ny, muF = xmuf*Qtest, muR = xmur*Qtest)
    call write_pdf(idev, Qtest, ymax, ny)
    
  end subroutine write_f

  !> @brief write the F1 structure function to idev
  !!
  !! Writes the W and Z F1 structure functions to a file with equal spacing in log(x)
  !!
  !! @param[in]     idev     file/device to write to
  !! @param[in]     Qtest    Value of Q to use
  !! @param[in]     ymax     Largest value of y = - log(x) to print
  !! @param[in]     ny       number of bins in y to print
  !!
  subroutine write_pdf(idev, Qtest, ymax, ny)
    real(dp), intent(in) :: Qtest, ymax
    integer, intent(in)  :: idev, ny
    real(dp) :: ytest, xval, xpdf_at_xQ(-6:6)
    integer  :: iy
    !F1 Wp Wm Z
    write(idev,'(a,f10.4,a,f10.4)') '# Q = ', Qtest
    write(idev,'(a)') '# x pdf(-6:6)'
    
    do iy = ny, 1, -1
       ytest = iy * ymax / ny
       xval = exp(-ytest)
       call hoppetEval(xval,Qtest,xpdf_at_xQ)
       write(idev,'(14es12.4)') xval, xpdf_at_xQ
    end do
    write(idev,*)
    write(idev,*)
  end subroutine write_pdf
  

  
end program structure_functions_example


