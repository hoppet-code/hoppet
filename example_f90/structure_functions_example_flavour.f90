!! An example program using structure functions up to NLO with flavour decomposition. !!
!!
program structure_functions_example
  use hoppet
  use pdfs_for_benchmarks
  use streamlined_interface
  use structure_functions
  implicit none
  real(dp) :: Qmax, xmur, xmuf, Qmin, ymax, Q, mc, mb, mt, asQ, Q0,&
       & muR_Q, dy, dlnlnQ, minQval, maxQval, muR, muF, x, ytest
  integer  :: nloop_coefs, sc_choice,nloop, order, i, ny, iy 
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

  call hoppetEvolve(asQ, Q0, nloop,muR_Q, benchmark_pdf_unpolarized_lha, Q0)

  ! Initialise the structure functions using separate orders
  ! NB: this uses the PDFs that were set up in the streamlined interface
  ! with the hoppetEvolve routine
  call InitStrFct(nloop_coefs, .true., xR = xmur, xF = xmuf, flavour_decomposition = .true.) ! Ask
                                                                                             ! for
                                                                                             ! flavour
                                                                                             ! decomposed
                                                                                             ! structure
                                                                                             ! functions
                                                                                             ! here

  Q = 100.0_dp
  muR = Q
  muF = Q

  ymax = log(1e5)
  ny = 10
  print*, 'Standard structure functions at NLO for Q=100 GeV and x=0.01'
  do iy = ny, 1, -1
     ytest = iy * ymax / ny
     x = exp(-ytest)
     FLO  = F_LO(x,Q,muR,muF)
     FNLO = F_NLO(x,Q,muR,muF)
     F    = StrFct(x,Q,muR,muF)
     write(*,'(15es12.4)') x, F
  enddo
  print*, ''
  print*, 'Struture functions obtained from the flavour decomposed ones'
  do iy = ny, 1, -1
     ytest = iy * ymax / ny
     x = exp(-ytest)
     do i =-6,6
        FLOflav(i,:) = F_LO_flav(x,Q,muR,muF,i)
        FNLOflav(i,:) = F_NLO_flav(x,Q,muR,muF,i)
        Fflav(i,:)   = StrFct_flav(x,Q,muR,muF,i)
     enddo
     call fill_structure_functions(x,Fflav(:,1),Fflav(:,2),Fflav(:,3),F)
     write(*,'(15es12.4)') x, F
  enddo

end program structure_functions_example

! Below a routine that fills all the NC and CC structure functions
! from the flavour decomposed ones. These routines are essentially
! adapted from those that can be found in structure_functions.f90
subroutine fill_structure_functions(x,FL,F2,F3,SF)
     use hoppet
  use streamlined_interface
  use structure_functions
  use coefficient_functions_holder
  implicit none
  real(dp) :: FL(-6:6),F2(-6:6),F3(-6:6),SF(-6:7),x
  real(dp), parameter :: viW = 1/sqrt(two), aiW = viW
  real(dp), save      :: vi2_ai2_avg_W, two_vi_ai_avg_W
  real(dp), save      :: vi2_ai2_Z_down, vi2_ai2_Z_up
  real(dp), save      :: two_vi_ai_Z_down, two_vi_ai_Z_up
  real(dp), save      :: two_vi_Z_down, two_vi_Z_up
  real(dp), save      :: two_ai_Z_down, two_ai_Z_up
  real(dp), parameter :: e_up = 2.0_dp/3.0_dp, e_down = - 1.0_dp/3.0_dp
  real(dp), parameter :: e2_up = 4.0_dp/9.0_dp, e2_down = 1.0_dp/9.0_dp
  real(dp) :: sin_thw_sq, mw, mz , two_xvals
  real(dp) :: ulike, dlike, ubarlike, dbarlike
  integer :: nf_lcl
  
  mw                = 80.377_dp
  mz                = 91.1876_dp

  two_xvals = 2.0_dp * x

  ! compute sin(\theta_w) from W/Z mass
  sin_thw_sq = 1.0_dp - (mw/mz)**2

  ! evaluate parameters needed for the structure functions
  ! cf. Eqs. (3.10+11) of 1109.3717
  vi2_ai2_Z_up     = one/four + (half - (four/three) * sin_thw_sq)**2
  vi2_ai2_Z_down   = one/four + (half - (two/three)  * sin_thw_sq)**2
  two_vi_ai_Z_up   = half - (four/three) * sin_thw_sq
  two_vi_ai_Z_down = half - (two/three)  * sin_thw_sq

  ! cf. Eq.3.20 + 3.17 & 3.18
  ! 1/nf \sum_j=1^nf (vj^2 + aj^2)
  vi2_ai2_avg_W = viW**2 + aiW**2
  ! 1/nf \sum_j=1^nf 2*vj*aj
  two_vi_ai_avg_W = two * viW * aiW

  ! these parameters are needed explicitly for the gamma-Z structure function
  two_vi_Z_up = one - (8.0_dp/three) * sin_thw_sq
  two_vi_Z_down   = -one + (four/three) * sin_thw_sq
  two_ai_Z_up = one
  two_ai_Z_down   = -one

  SF = zero
  nf_lcl = 5
  !--- deal with Z case -----------------------------------------
  SF(iF2Z) = (dlike(F2,nf_lcl) + dbarlike(F2,nf_lcl))*vi2_ai2_Z_down + &
       &       (ulike(F2,nf_lcl) + ubarlike(F2,nf_lcl))*vi2_ai2_Z_up

  ! temporarily put FL into F1;
  SF(iF1Z) = (dlike(FL,nf_lcl) + dbarlike(FL,nf_lcl))*vi2_ai2_Z_down +&
       & (ulike(FL,nf_lcl) + ubarlike(FL,nf_lcl))*vi2_ai2_Z_up
  ! then convert to F1
  SF(iF1Z) = (SF(iF2Z) - SF(iF1Z)) / two_xvals

  SF(iF3Z) = (dlike(F3,nf_lcl) - dbarlike(F3,nf_lcl))*two_vi_ai_Z_down + (ulike(F3,nf_lcl)&
       & - ubarlike(F3,nf_lcl))*two_vi_ai_Z_up
  SF(iF3Z) = SF(iF3Z)/x

  !--- deal with EM case -----------------------------------------
  SF(iF2EM) = (dlike(F2,nf_lcl) + dbarlike(F2,nf_lcl))*e2_down + (ulike(F2,nf_lcl) +&
       & ubarlike(F2,nf_lcl))*e2_up

  ! temporarily put FL into F1;
  SF(iF1EM) = (dlike(FL,nf_lcl) + dbarlike(FL,nf_lcl))*e2_down +&
       & (ulike(FL,nf_lcl) + ubarlike(FL,nf_lcl))*e2_up
  ! then convert to F1
  SF(iF1EM) = (SF(iF2EM) - SF(iF1EM)) / two_xvals

  !--- deal with gamma-Z case
  !-----------------------------------------
  ! equations taken from section 19 of PDG (19.18)
  SF(iF2gZ) = (dlike(F2,nf_lcl) + dbarlike(F2,nf_lcl))*e_down*two_vi_Z_down +&
       & (ulike(F2,nf_lcl) + ubarlike(F2,nf_lcl))*e_up * two_vi_Z_up

  ! temporarily put FL into F1;
  SF(iF1gZ) = (dlike(FL,nf_lcl) + dbarlike(FL,nf_lcl))*e_down*two_vi_Z_down +&
       & (ulike(FL,nf_lcl) + ubarlike(FL,nf_lcl))*e_up * two_vi_Z_up
  ! then convert to F1
  SF(iF1gZ) = (SF(iF2gZ) - SF(iF1gZ)) / two_xvals

  SF(iF3gZ) = (dlike(F3,nf_lcl) - dbarlike(F3,nf_lcl))*e_down*two_ai_Z_down + &
       &        (ulike(F3,nf_lcl) - ubarlike(F3,nf_lcl))*e_up * two_ai_Z_up
  SF(iF3gZ) = SF(iF3gZ)/x

  ! for the W cases, it only makes sense to sum over an even number
  ! of light flavours; so save the actual number of flavours, switch
  ! the module-local nf_lcl variable to the nearest (lower) event
  ! number
  ! for our W calculations; switch back later
  nf_lcl = 4

  ! --- deal with Wm case
  !-----------------------------------------
  SF(iF2Wm) = (ulike(F2,nf_lcl) + dbarlike(F2,nf_lcl))*vi2_ai2_avg_W
  ! temporarily put FL into F1;
  SF(iF1Wm) = (ulike(FL,nf_lcl) + dbarlike(FL,nf_lcl))*vi2_ai2_avg_W
  ! then convert to F1
  SF(iF1Wm) = (SF(iF2Wm) - SF(iF1Wm)) / two_xvals
  SF(iF3Wm) = (ulike(F3,nf_lcl) - dbarlike(F3,nf_lcl))*two_vi_ai_avg_W
  SF(iF3Wm) = SF(iF3Wm)/x

  !--- deal with Wp case -----------------------------------------
  SF(iF2Wp) = (dlike(F2,nf_lcl) + ubarlike(F2,nf_lcl))*vi2_ai2_avg_W
  ! temporarily put FL into F1;
  SF(iF1Wp) = (dlike(FL,nf_lcl) + ubarlike(FL,nf_lcl))*vi2_ai2_avg_W
  ! then convert to F1
  SF(iF1Wp) = (SF(iF2Wp) - SF(iF1Wp)) / two_xvals
  SF(iF3Wp) = (dlike(F3,nf_lcl) - ubarlike(F3,nf_lcl))*two_vi_ai_avg_W
  SF(iF3Wp) = SF(iF3Wp)/x

end subroutine fill_structure_functions

!----------------------------------------------------------------------
function ulike(f,nf) result(res)
  double precision, intent(in) :: f(-6:6)
  integer, intent(in)  :: nf
  double precision             :: res
  res = sum(f(2: nf: 2))
end function ulike
!----------------------------------------------------------------------
function dlike(f,nf) result(res)
  double precision, intent(in) :: f(-6:6)
  integer, intent(in)  :: nf
  double precision             :: res
  res = sum(f(1: nf: 2))
end function dlike
!----------------------------------------------------------------------
function ubarlike(f,nf) result(res)
  double precision, intent(in) :: f(-6:6)
  integer, intent(in)  :: nf
  double precision             :: res
  res = sum(f(-2:-nf:-2))
end function ubarlike
!----------------------------------------------------------------------
function dbarlike(f,nf) result(res)
  double precision, intent(in) :: f(-6:6)
  integer, intent(in)  :: nf
  double precision             :: res
  res = sum(f(-1:-nf:-2))
end function dbarlike

  
