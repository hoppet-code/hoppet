!! An example program using a tabulation
!!
!! NB: to compile it you will also need to link with LHAPDF, because
!! we use that for the initial condition.
program tabulation_example
  !! use hoppet_v1, and rename a couple of internal functions which
  !! would otherwise conflict with LHAPDF
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  implicit none
  real(dp) :: dy, ymax
  integer  :: order, nloop
  !! holds information about the grid
  type(grid_def) :: grid, gdarray(3)
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

  !! define the interface for LHA pdf
  interface
     subroutine EvolvePDF(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine EvolvePDF
  end interface

  ! set up parameters for grid
  order = -5
  ymax = 12
  dy = 0.1_dp
  
  ! set up the grid itself
  call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
  call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  call InitGridDef(grid,gdarray(1:3),locked=.true.)

  ! initialise splitting-function holder
  nloop = 2
  call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop,nflo=3,nfhi=6)

  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.
  call AllocPdfTable(grid, table, Qmin=1.0_dp, Qmax=28000.0_dp, & 
       & dlnlnQ = 0.05_dp, freeze_at_Qmin=.true.)

  ! set up LHAPDF with cteq
  call InitPDFsetByName("cteq61.LHgrid")
  call InitPDF(0)

  ! allocate and set up the initial pdf from LHAPDF ...
  Q0 = 5.0_dp
  call AllocPDF(grid, pdf0)
  call InitPDF_LHAPDF(grid, pdf0, EvolvePDF, Q0)

  ! allocate the running coupling
  call InitRunningCoupling(coupling,alfas=0.118_dp)

  ! create the tabulation
  call EvolvePdfTable(table, Q0, pdf0, dh, coupling, nloop=nloop)

  ! get the value of the tabulation at some point
  x = 0.1_dp
  Q = 100.0_dp
  call EvalPdfTable_xQ(table,x,Q,pdf_at_xQ)
  write(6,*) "Gluon at x,Q ", x, Q, "is ", pdf_at_xQ(0)
end program tabulation_example
