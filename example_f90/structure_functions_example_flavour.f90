!! An example program using structure functions up to N3LO. !!
!!
program structure_functions_example
  use hoppet_v1
  use dummy_pdfs
  use streamlined_interface
  use structure_functions
  implicit none
  real(dp) :: Qmax, xmur, xmuf, Qmin, ymax, Q, mc, mb, mt, asQ, Q0,&
       & muR_Q, dy, dlnlnQ, minQval, maxQval, muR, muF, x
  integer  :: nloop_coefs, sc_choice,nloop, order, i
  real(dp) :: F(-6:7), FLO(-6:7), FNLO(-6:7)
  real(dp) :: Fflav(-6:6,3), FLOflav(-6:6,3), FNLOflav(-6:6,3)
  
  Qmax = 13000.0_dp 
  nloop_coefs = 2
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
  dlnlnQ = dy/4.0_dp
  nloop = 3 
  minQval = min(xmuF*Qmin, Qmin)
  maxQval = max(xmuF*Qmax, Qmax)
  
  ! initialise the grid and dglap holder, using the streamlined
  ! interface for simplicity
  call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
         &         order,factscheme_MSbar)

  ! Setup all constants and parameters needed by the structure functions
  call StartStrFct(nloop_coefs, scale_choice = sc_choice, param_coefs&
       & = .true., constant_mu = 91.1876_dp, nflav = 5)

  ! Evolve the PDF
  asQ = 0.35_dp
  Q0 = sqrt(2.0_dp)
  muR_Q = 1.0_dp

  call hoppetEvolve(asQ, Q0, nloop,muR_Q, lha_unpolarized_dummy_pdf, Q0)
  
  ! Initialise the structure functions using separate orders 
  ! NB: this uses the PDFs that were set up in the streamlined interface
  ! with the hoppetEvolve routine
  call InitStrFct(nloop_coefs, .true., xR = xmur, xF = xmuf, flavour_decomposition = .true.)

  x = 0.01_dp
  Q = 50.0_dp
  muR = 91.1876_dp
  muF = 91.1876_dp
  do i =-6,6
     FLOflav(i,:) = F_LO_flav(x,Q,muR,muF,i)
     FNLOflav(i,:) = F_NLO_flav(x,Q,muR,muF,i)
     Fflav(i,:)   = StrFct_flav(x,Q,muR,muF,i)
  enddo

  FLO  = F_LO(x,Q,muR,muF)
  FNLO = F_NLO(x,Q,muR,muF)
  F    = StrFct(x,Q,muR,muF)
  print*, 'OLD'
  print*, 'FLOW  = ', FLO(iF2Wm) - 2.0_dp * x *FLO(iF1Wm), FLO(iF2Wm), FLO(iF3Wm)
  print*, 'FNLOW = ', FNLO(iF2Wm) - 2.0_dp * x *FNLO(iF1Wm), FNLO(iF2Wm), FNLO(iF3Wm)
  print*, 'FOW   = ', F(iF2Wm) - 2.0_dp * x *F(iF1Wm), F(iF2Wm), F(iF3Wm)
  Print*, 'NEW'
  print*, 'FLOW  = ', FLOflav(-3,1) + FLOflav(-1,1) + FLOflav(2,1) +&
       & FLOflav(4,1), FLOflav(-3,2) + FLOflav(-1,2) + FLOflav(2,2) +&
       & FLOflav(4,2), (-FLOflav(-3,3) - FLOflav(-1,3) + FLOflav(2,3)&
       & + FLOflav(4,3))/x
  print*, 'FNLOW = ', FNLOflav(-3,1) + FNLOflav(-1,1) + FNLOflav(2,1) +&
       & FNLOflav(4,1), FNLOflav(-3,2) + FNLOflav(-1,2) + FNLOflav(2&
       &,2) + FNLOflav(4,2), (-FNLOflav(-3,3) - FNLOflav(-1,3) +&
       & FNLOflav(2,3) + FNLOflav(4,3))/x
  print*, 'FW    = ', Fflav(-3,1) + Fflav(-1,1) + Fflav(2,1) + Fflav(4&
       &,1), Fflav(-3,2) + Fflav(-1,2) + Fflav(2,2) + Fflav(4,2), (&
       &-Fflav(-3,3) - Fflav(-1,3) + Fflav(2,3) + Fflav(4,3))/x
  !  ! write out the structure functions
!  ymax = log(1e5) !ymax=20
!  Q = 100.0_dp
!  call write_f1(6, Q, ymax, 10)
!  call write_f2(6, Q, ymax, 10)
!  call write_f3(6, Q, ymax, 10)

end program structure_functions_example

subroutine fill_structure_functions(FL,F2,F3,SF)
  use hoppet_v1
  use streamlined_interface
  use structure_functions
  implicit none
  real(dp) :: FL(-6:6),F2(-6:6),F3(-6:6),SF(6:-7)
  
end subroutine fill_structure_functions

!----------------------------------------------------------------------
function ulike(f) result(res)
  use hoppet_v1
  real(dp), intent(in) :: f(0:,ncompmin:)
  real(dp)             :: res(0:ubound(f,dim=1))
  res = sum(f(:, 2: nf_lcl: 2), dim=2)
end function ulike
!----------------------------------------------------------------------
function dlike(f) result(res)
  use hoppet_v1
  real(dp), intent(in) :: f(0:,ncompmin:)
  real(dp)             :: res(0:ubound(f,dim=1))
  res = sum(f(:, 1: nf_lcl: 2), dim=2)
end function dlike
!----------------------------------------------------------------------
function ubarlike(f) result(res)
  use hoppet_v1
  real(dp), intent(in) :: f(0:,ncompmin:)
  real(dp)             :: res(0:ubound(f,dim=1))
  res = sum(f(:,-2:-nf_lcl:-2), dim=2)
end function ubarlike
!----------------------------------------------------------------------
function dbarlike(f) result(res)
  use hoppet_v1
  real(dp), intent(in) :: f(0:,ncompmin:)
  real(dp)             :: res(0:ubound(f,dim=1))
  res = sum(f(:,-1:-nf_lcl:-2), dim=2)
end function dbarlike

  
