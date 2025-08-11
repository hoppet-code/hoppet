

!! An example program using HOPPET's streamlined interface. 
!! It outputs a subset of
!! table 15 of hep-ph/0511119 and this output should be identical
!! to the contents of the file tabulation_example_streamlined.default_output
!!
!! NB: for the full functionality used in generating the HeraLHC and
!!     Les Houches comparison tables, see ../benchmarks/benchmarks.f90
!!     and carefully read the comments at the beginning. Subtleties
!!     exist in particular wrt the treatment of scales muF/=muR.
!!
!! NB: commented code shows usage with LHAPDF (e.g. for cross-checking 
!!     some public PDF set's evolution) -- to use this part of the code
!!     you will also need to link with LHAPDF
!!
!! This program uses the streamlined HOPPET interface rather
!! than the more general ones, used in the example programs
!! tabulation_example.f90 and tabulation_example_2.f90
!!
!!



module lhasub4streamlined
  implicit none

contains
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
  subroutine LHAsub(x,Q,xpdf)
  use types; use consts_dp; use hoppet; implicit none
    real(dp), intent(in)  :: x,Q
    real(dp), intent(out) :: xpdf(-6:6)
    real(dp) :: uv, dv
    real(dp) :: ubar, dbar, glu
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    ! Set to zero the xpdf array
    xpdf = zero
    
  
    !-- remember that these are all x*q(x)
    uv   = 1* N_uv * x**0.8_dp * (1-x)**3
    dv   = 1* N_dv * x**0.8_dp * (1-x)**4
    dbar = 1* N_db * x**(-0.1_dp) * (1-x)**6
    ubar = 1* dbar * (1-x)
    glu  = 1* N_g * x**(-0.1_dp) * (1-x)**5
    
    ! labels iflv_g, etc., come from the hoppet module, inherited
    ! from the main program
    xpdf(iflv_g) = glu
    xpdf(-iflv_s) = 0.2_dp*(dbar + ubar)
    xpdf(iflv_s) = xpdf(-iflv_s)
    xpdf(iflv_u) = uv + ubar
    xpdf(-iflv_u) = ubar
    xpdf(iflv_d) = dv + dbar
    xpdf(-iflv_d) = dbar
  end subroutine LHAsub
end module lhasub4streamlined




program tabulation_example_streamlined
  use hoppet; use lhasub4streamlined
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  !use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  implicit none
  real(dp) :: dy, ymax, dlnlnQ, Qmin, Qmax, muR_Q
  real(dp) :: asQ, Q0alphas, Q0pdf
  real(dp) :: mc,mb,mt
  integer  :: order, nloop, scheme
  !! hold results at some x, Q
  real(dp) :: Q, xpdf_at_xQ(-6:6)
  real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: ix
  !! define the interfaces for LHA pdf (by default not used)
  !! (NB: unfortunately this conflicts with an internal hoppet name,
  !! so make sure that you "redefine" the internal hoppet name,
  !! as illustrated in the commented "use" line above:
  !! use hoppet, EvolvePDF_hoppet => EvolvePDF, ...)
  ! interface
  !    subroutine EvolvePDF(x,Q,res)
  !      use types; implicit none
  !      real(dp), intent(in)  :: x,Q
  !      real(dp), intent(out) :: res(*)
  !    end subroutine EvolvePDF
  ! end interface


  ! Streamlined initialization
  ! including  parameters for x-grid
  order = -6
  ymax  = 12.0_dp
  dy    = 0.025_dp
  ! and parameters for Q tabulation
  Qmin=1.0_dp
  Qmax=28000.0_dp
  dlnlnQ = dy/4.0_dp
  ! and number of loops to initialise! (1=LO, 2=NLO)
  nloop = 2
  scheme = factscheme_FragMSbar
  !scheme = factscheme_MSbar
  call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,&
       &         order,scheme)
  write(6,'(a)') "Streamlined initialization completed!"
  
  ! Set heavy flavour scheme
  mc = 1.414213563_dp   ! sqrt(2.0_dp) + epsilon
  mb = 4.5_dp
  mt = 175.0_dp
  call hoppetSetVFN(mc, mb, mt)

  ! Streamlined evolution
 
  ! Set parameters of running coupling
  asQ = 0.35_dp
  Q0alphas = sqrt(2.0_dp)
  muR_Q = 1.0_dp
  Q0pdf = 5.0_dp ! The initial evolution scale
  call hoppetEvolve(asQ, Q0alphas, nloop,muR_Q, LHAsub, Q0pdf)
  write(6,'(a)') "Evolution done!"
  
  
  ! get the value of the tabulation at some point
  Q = 100.0_dp
  !Q = 5.01_dp
  write(6,'(a,i1)') "# nloop = ",nloop
  write(6,'(a,f8.3,a)') "           Evaluating PDFs at Q = ",Q," GeV"
  write(6,'(a5,2a12,a14,a10,a12)') "x",&
       & "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
  do ix = 1, size(heralhc_xvals)
     call hoppetEval(heralhc_xvals(ix),Q,xpdf_at_xQ)
     write(6,'(es7.1,5es12.4)') heralhc_xvals(ix), &
          &  xpdf_at_xQ(2)-xpdf_at_xQ(-2), &
          &  xpdf_at_xQ(1)-xpdf_at_xQ(-1), &
          &  2*(xpdf_at_xQ(-1)+xpdf_at_xQ(-2)), &
          &  (xpdf_at_xQ(-4)+xpdf_at_xQ(4)), &
          &  xpdf_at_xQ(0)
  end do
  
 
! ! Now perform cached evolution
!  call hoppetPreEvolve(asQ, Q0alphas, nloop,muR_Q,Q0pdf)
!  write(6,'(a)') "Pre-evolution prepared!"
!
!  call hoppetCachedEvolve(LHAsub)
!  write(6,'(a)') "Cached evolution performed!"
!
!  ! Same output as before
!  Q = 100.0_dp
!  write(6,'(a)')
!  write(6,'(a,f8.3,a)') "           Evaluating PDFs at Q = ",Q," GeV"
!  write(6,'(a5,2a12,a14,a10,a12)') "x",&
!       & "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
!  do ix = 1, size(heralhc_xvals)
!     call hoppetEval(heralhc_xvals(ix),Q,xpdf_at_xQ)
!     write(6,'(es7.1,5es12.4)') heralhc_xvals(ix), &
!          &  xpdf_at_xQ(2)-xpdf_at_xQ(-2), &
!          &  xpdf_at_xQ(1)-xpdf_at_xQ(-1), &
!          &  2*(xpdf_at_xQ(-1)+xpdf_at_xQ(-2)), &
!          &  (xpdf_at_xQ(-4)+xpdf_at_xQ(4)), &
!          &  xpdf_at_xQ(0)
!  end do
  
end program tabulation_example_streamlined

