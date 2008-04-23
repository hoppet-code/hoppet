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
program tabulation_example_streamlined
  use hoppet_v1
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  !use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  implicit none
  real(dp) :: dy, ymax, dlnlnQ, Qmin, Qmax, muR_Q
  real(dp) :: asQ, Q0alphas, Q0pdf
  real(dp) :: mc,mb,mt
  integer  :: order, nloop
  !! holds information about the grid
  type(grid_def) :: grid, gdarray(4)
  !! hold results at some x, Q
  real(dp) :: Q, xpdf_at_xQ(-6:6)
  real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: ix

  !! define the interfaces for LHA pdf (by default not used)
  !! (NB: unfortunately this conflicts with an internal hoppet name,
  !! so make sure that you "redefine" the internal hoppet name, 
  !! as illustrated in the commented "use" line above:
  !! use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, ...)
  ! interface
  !    subroutine EvolvePDF(x,Q,res)
  !      use types; implicit none
  !      real(dp), intent(in)  :: x,Q
  !      real(dp), intent(out) :: res(*)
  !    end subroutine EvolvePDF
  ! end interface

  ! set up parameters for grid
  order = -6
  ymax  = 12.0_dp
  dy    = 0.1_dp

  ! set up the grid itself -- we use 4 nested subgrids
  call InitGridDef(gdarray(4),dy/27.0_dp,0.2_dp, order=order)
  call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
  call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  call InitGridDef(grid,gdarray(1:4),locked=.true.)

  ! Streamlined initialization
  Qmin=1_dp
  Qmax=28000_dp
  dlnlnQ = dy/4.0_dp
  nloop = 3
  call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,&
       &         order,factscheme_MSbar)
  write(6,'(a)') "Streamlined initialization completed!"
  
  ! Set heavy flavour scheme
  mc = 1.414213563_dp
  mb = 4.5_dp
  mt = 175.0_dp
  call hoppetSetVFN(mc, mb, mt)
  write(6,'(a)') "Heavy flavour scheme set!"

  ! Streamlined evolution
 
  ! Set parameters of running coupling
  asQ = 0.35_dp
  Q0alphas = sqrt(2.0_dp)
  muR_Q = 1.0_dp
  Q0pdf = sqrt(2.0_dp) ! The initial evolution scale
  call hoppetEvolve(asQ, Q0alphas, nloop,muR_Q,&
       &       LHAsub, Q0pdf)
  write(6,'(a)') "Evolution prepared!"
  
  
  ! get the value of the tabulation at some point
  Q = 100.0_dp
  write(6,'(a)')
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
  
 
 ! Now perform cached evolution
  call hoppetPreEvolve(asQ, Q0alphas, nloop,muR_Q,Q0pdf)
  write(6,'(a)') "Pre-evolution prepared!"

  call hoppetCachedEvolve(LHAsub)
  write(6,'(a)') "Cached evolution ready!"

  ! Same output as before
  Q = 100.0_dp
  write(6,'(a)')
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
  

contains 
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
  subroutine LHAsub(x,Q,xpdf)
    use types; implicit none
    real(dp), intent(in)  :: x,Q
    real(dp), intent(out) :: xpdf(-6:6)
    real(dp) :: uv, dv
    real(dp) :: ubar, dbar
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    ! Set to zero the xpdf array
    xpdf = zero
    

    !-- remember that these are all x*q(x)
    uv = N_uv * x**0.8_dp * (1-x)**3
    dv = N_dv * x**0.8_dp * (1-x)**4
    dbar = N_db * x**(-0.1_dp) * (1-x)**6
    ubar = dbar * (1-x)

    ! labels iflv_g, etc., come from the hoppet_v1 module, inherited
    ! from the main program
    xpdf(iflv_g) = N_g * x**(-0.1_dp) * (1-x)**5
    xpdf(-iflv_s) = 0.2_dp*(dbar + ubar)
    xpdf(iflv_s) = xpdf(-iflv_s)
    xpdf(iflv_u) = uv + ubar
    xpdf(-iflv_u) = ubar
    xpdf(iflv_d) = dv + dbar
    xpdf(-iflv_d) = dbar
  end subroutine LHAsub

end program tabulation_example_streamlined


