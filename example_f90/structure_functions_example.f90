!! An example program using structure functions up to N3LO. !!
!!
program structure_functions_example
  use hoppet_v1
  use dummy_pdfs
  use streamlined_interface
  use structure_functions
  implicit none
  real(dp) :: Qmax, xmur, xmuf, Qmin, ymax, Q, mc, mb, mt, asQ, Q0,&
       & muR_Q
    integer  :: order_max, sc_choice,nloop
  
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

  ! Setup all constants and parameters needed by the structure functions
  call StartStrFct(Qmax, order_max, xR = xmur, xF = xmuf, sc_choice = sc_choice, &
       param_coefs = .true., Qmin_PDF = Qmin)

  ! Evolve the PDF
  nloop = 3 ! NNLO evolution
  asQ = 0.35_dp
  Q0 = sqrt(2.0_dp)
  muR_Q = 1.0_dp

  call hoppetEvolve(asQ, Q0, nloop,muR_Q, lha_unpolarized_dummy_pdf, Q0)
  
  ! Initialise the structure functions using separate order
  call InitStrFct(order_max, .true.)

  ! write out the structure functions
  write(6,*) "Writing structure functions to structure_functions_example.dat"
  open(unit = 99, file = 'structure_functions_example.dat')

  ymax = log(1e5) !ymax=20
  Q = 100.0_dp
  call write_f1(99, Q, ymax, 10)
  call write_f2(99, Q, ymax, 10)
  call write_f3(99, Q, ymax, 10)

end program structure_functions_example


