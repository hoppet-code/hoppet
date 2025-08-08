!! An example program using structure functions up to N3LO. !!
!!
program structure_functions_example
  use hoppet_v1
  use dummy_pdfs
  use streamlined_interface
  use structure_functions
  implicit none
  real(dp) :: Qmax, xmur, xmuf, Qmin, ymax, Q, mc, mb, mt, asQ, Q0,&
       & muR_Q, dy, dlnlnQ, minQval, maxQval
    integer  :: nloop_coefs, sc_choice,nloop, order
  
  Qmax = 13000.0_dp 
  nloop_coefs = 4
  xmur = one
  xmuf = one
  sc_choice = scale_choice_Q ! Uses Q as the central scale choice
  Qmin = one

  ! Set heavy flavour scheme
  mc = 1.414213563_dp   ! sqrt(2.0_dp) + epsilon
  mb = 4.5_dp
  mt = 175.0_dp
  call hoppetSetPoleMassVFN(mc, mb, mt)

  ! Streamlined initialization
  ! including  parameters for x-grid
  order = -6 ! interpolation order, not perturbative order in alphas!
  ymax  = 16.0_dp
  dy    = 0.05_dp  ! dble_val_opt("-dy",0.1_dp)
  !dlnlnQ = dy/4.0_dp
  dlnlnQ = dy/8.0_dp
  nloop = 3 
  minQval = min(xmuF*Qmin, Qmin)
  maxQval = max(xmuF*Qmax, Qmax)
  
  ! initialise the grid and dglap holder, using the streamlined
  ! interface for simplicity
  call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
         &         order,factscheme_MSbar)

  ! Setup all constants and parameters needed by the structure functions
  call StartStrFct(nloop_coefs, scale_choice = sc_choice, param_coefs&
       & = .true.)
!   For exact coeffs. If one wants to run with the parametrized
!   coefficients for the NS piece, one should set "tiny" in
!   coefficient_functions_holder.f90 to one.
!  call StartStrFct(nloop_coefs, scale_choice = sc_choice, param_coefs&
!       & = .false.)

  ! Evolve the PDF
  asQ = 0.35_dp
  Q0 = sqrt(2.0_dp)
  muR_Q = 1.0_dp

  call hoppetEvolve(asQ, Q0, nloop,muR_Q, lha_unpolarized_dummy_pdf, Q0)
  
  ! Initialise the structure functions using separate orders 
  ! NB: this uses the PDFs that were set up in the streamlined interface
  ! with the hoppetEvolve routine
  call InitStrFct(nloop_coefs, .true., xR = xmur, xF = xmuf)

  ! write out the structure functions
  ymax = log(1e5) !ymax=20
  Q = 100.0_dp
  call write_f1(12, Q, ymax, 100)
  call write_f2(12, Q, ymax, 100)
  call write_f3(12, Q, ymax, 100)

end program structure_functions_example


