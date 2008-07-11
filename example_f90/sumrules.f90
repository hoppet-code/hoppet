!! ------------------------------------------------------------
!!
!! An example program which shows how to compute
!! sum rules (and general moments) for a given set of PDFS
!!
!! -----------------------------------------------------

program sumrules

  use hoppet_v1
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  !use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF

  implicit none
  real(dp) :: dy, ymax
  integer  :: order, nloop
  !! holds information about the grid
  type(grid_def) :: grid, gdarray(4)
  !! holds the splitting functions
  type(dglap_holder) :: dh
  !! hold the PDF tabulation
  type(pdf_table)       :: table
  !! hold the coupling
  real(dp)               :: quark_masses(4:6)
  type(running_coupling) :: coupling
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)
  real(dp), pointer :: pdf_flav(:)
  real(dp) :: Q0
  !! hold results at some x, Q
  real(dp) :: Q, pdf_at_xQ(-6:6)
  real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: ix

  real(dp) :: sum_rule, moment_index
  integer  :: nf_rep,nf_lcl

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

  ! initialise the splitting-function holder
  nloop = 3
  call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop,nflo=3,nfhi=6)
  write(6,'(a)') "Splitting functions initialised!"

  ! initialise a PDF from the function below (must be contained,
  ! in a "used" module, or with an explicitly defined interface)
  call AllocPDF(grid, pdf0)
  pdf0 = unpolarized_dummy_pdf(xValues(grid))
  Q0 = sqrt(2.0_dp)  ! the initial scale


  ! Compute sum rules of the initial PDF set
  write(6,'(a)') "  "
  write(6,'(a,1x,f3.1,1x,a)') " Computation of  sum rules at Q02 =  ",&
       & Q0**2," GeV"
  write(6,'(a)') "  "

  ! Allocate a grid quantity for the pdf combination
  call AllocGridQuant(grid,pdf_flav)
  
  !
  ! (Truncated) moments are defined in terms of y as
  ! 
  !  M(moment_index) = \int_0^ymax dy exp( - y * moment_index) * xpdf(y)
  ! 
  ! or in terms of x, where y = ln(1/x)
  !
  !   M(moment_index) = \int_xmin^1 dx x^( moment_index - 1) * xpdf(x)
  !

  ! Define the PDF combinations whose moments we are to compute
  ! Check PDF representtation of input PDF set
  
  ! Momentum sum rule
  nf_rep = GetPdfRep(pdf0)
  ! write(6,*) "nf_rep = ",nf_rep
  if(nf_rep.lt.0) then
     ! PDF set in human representation -> Sum all pdfs
     pdf_flav(:)=sum(pdf0(:,:),dim=2)
  elseif(nf_rep.ge.0) then
     ! PDF in evln representation -> Singlet + Gluon
     pdf_flav(:)=pdf0(:,0) + pdf0(:,1)
  endif
  moment_index=1
  sum_rule = GetTruncMoment(grid,pdf_flav,moment_index)
  
  write(6,*) "Momentum sum rule = ",sum_rule

  ! Valence sum rules
  ! uv + dv sum rule
  nf_rep = GetPdfRep(pdf0)
  if(nf_rep.lt.0) then
     ! PDF set in human representation 
     pdf_flav(:)= pdf0(:,2) - pdf0(:,-2) + pdf0(:,1) - pdf0(:,-1)
  endif
  moment_index=0
  sum_rule = GetTruncMoment(grid,pdf_flav,moment_index)
  
  write(6,*) "uv + dv sum rule = ",sum_rule
  

  !
  ! Sum rules at arbitrary values of Q
  !

  Q=100_dp
  write(6,'(a)') "  "
  write(6,'(a,1x,f7.1,1x,a)') " Computation of  sum rules at Q2 =  ",&
       & Q**2," GeV"
  write(6,'(a)') "  "
  

  ! allocate and initialise the running coupling with a given
  ! set of quark masses (NB: charm mass just above Q0).
  quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
  call InitRunningCoupling(coupling,alfas=0.35_dp,Q=Q0,nloop=nloop,&
       &                   quark_masses = quark_masses)

  ! Evolve the initial PDF up to scale Q to test sum rules
  call EvolvePDF(dh,pdf0,coupling,Q0,Q)

  ! Momentum sum rule
  nf_rep = GetPdfRep(pdf0)
  ! write(6,*) "nf_rep = ",nf_rep
  if(nf_rep.lt.0) then
     ! PDF set in human representation -> Sum all pdfs
     pdf_flav(:)=sum(pdf0(:,:),dim=2)
  elseif(nf_rep.ge.0) then
     ! PDF in evln representation -> Singlet + Gluon
     pdf_flav(:)=pdf0(:,0) + pdf0(:,1)
  endif
  moment_index=1
  sum_rule = GetTruncMoment(grid,pdf_flav,moment_index)
  
  write(6,*) "Momentum sum rule = ",sum_rule


  ! Valence sum rules
  ! uv + dv sum rule
  nf_rep = GetPdfRep(pdf0)
  if(nf_rep.lt.0) then
     ! PDF set in human representation 
     pdf_flav(:)= pdf0(:,2) - pdf0(:,-2) + pdf0(:,1) - pdf0(:,-1)
  endif
  moment_index=0
  sum_rule = GetTruncMoment(grid,pdf_flav,moment_index)
  
  write(6,*) "uv + dv sum rule = ",sum_rule
  write(6,*) " "
 
  ! Some cleaning
  call Delete(pdf0)
  call Delete(dh)
  call Delete(coupling)
  call Delete(grid)

contains 
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
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
    ! (not actually needed after setting pdf=0
    call LabelPdfAsHuman(pdf)

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)

    ! labels iflv_g, etc., come from the hoppet_v1 module, inherited
    ! from the main program
    pdf(:, iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf

end program sumrules


