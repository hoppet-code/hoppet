!! An example program using a tabulation
!!
!! NB: commented code shows usage with LHAPDF (e.g. for cross-checking 
!!     some public PDF set's evolution) -- to use this code
!!     you will also need to link with LHAPDF
!!
program tabulation_example
  !! use hoppet_v1, and rename a couple of internal functions which
  !! would otherwise conflict with LHAPDF (if one uses it)
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
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
  type(running_coupling) :: coupling
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)
  real(dp) :: Q0
  !! hold results at some x, Q
  real(dp) :: x, Q, pdf_at_xQ(-6:6)

  !! define the interfaces for LHA pdf
  interface
     subroutine EvolvePDF(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine EvolvePDF
  end interface

  ! set up parameters for grid
  order = -6
  ymax = 12
  dy = 0.1_dp
  
  ! set up the grid itself
  call InitGridDef(gdarray(4),dy/27.0_dp,0.2_dp, order=order)
  call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
  call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  call InitGridDef(grid,gdarray(1:4),locked=.true.)

  ! initialise splitting-function holder
  nloop = 2
  call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop,nflo=3,nfhi=6)
  write(0,'(a)') "Splitting function initialised!"

  !! set up LHAPDF with cteq
  !call InitPDFsetByName("cteq61.LHgrid")
  !call InitPDF(0)
  ! allocate and set up the initial pdf from LHAPDF ...
  !Q0 = 5.0_dp 
  !call AllocPDF(grid, pdf0)
  !call InitPDF_LHAPDF(grid, pdf0, EvolvePDF, Q0)

  ! initialise a PDF from the function below (must be contained,
  ! in a "used" module, or with an explicitly defined interface)
  call AllocPDF(grid, pdf0)
  pdf0 = unpolarized_dummy_pdf(xValues(grid))
  Q0 = 1.4_dp  ! the initial scale

  ! allocate and initialise the running coupling
  call InitRunningCoupling(coupling,alfas=0.118_dp)

  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.
  call AllocPdfTable(grid, table, Qmin=1.0_dp, Qmax=28000.0_dp, & 
       & dlnlnQ = dy/4.0_dp, freeze_at_Qmin=.true.)
  ! add information about the nf transitions to the table (improves
  ! interpolation quality)
  call AddNfInfoToPdfTable(table,coupling)

  ! create the tabulation
  call EvolvePdfTable(table, Q0, pdf0, dh, coupling, nloop=nloop)
  ! alternatively "pre-evolve" so that subsequent evolutions are faster
  !call PreEvolvePdfTable(table, Q0, dh, coupling)
  !call EvolvePdfTable(table,pdf0)

  ! get the value of the tabulation at some point
  x = 0.1_dp
  Q = 100.0_dp
  call EvalPdfTable_xQ(table,x,Q,pdf_at_xQ)
  write(6,*) "Gluon at x,Q ", x, Q, "is ", pdf_at_xQ(0)

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
    ! (not actually needed after setting pdf=0
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

end program tabulation_example


