module f77_pdftab
  use types; use consts_dp
  use pdf_tabulate
  use convolution; use pdf_general; use dglap_objects
  use dglap_holders; use pdf_general; use dglap_choices
  implicit none

  !! holds information about the grid
  type(grid_def),     save :: grid, gdarray(3)

  !! holds the splitting functions
  type(dglap_holder), save :: dh

  !! 0 is main pdf table, while i=1:3 contain convolutions with the
  !! i-loop splitting function
  type(pdf_table), save :: tables(0:3)
  logical,      save :: setup_done(0:3) = .false.
  integer,      save :: setup_nf(3)     = 0
end module f77_pdftab


!======================================================================
!! initialise the underlying grid, splitting functions and pdf-table
!! objects, using the dy and nloop parameters as explained below.
subroutine dglapStart(dy,nloop)
  use f77_pdftab
  implicit none
  !--------------------------------------
  real(dp), intent(in) :: dy     !! internal grid spacing: 0.1 is a sensible value
  integer,  intent(in) :: nloop  !! the maximum number of loops we'll want (<=3)
  !--------------------------------------
  real(dp) :: ymax, Qmin, Qmax, dlnlnQ
  integer  :: order
  ymax = 12.0d0
  Qmin = 1.0d0
  Qmax = 28000d0 ! twice LHC c.o.m.
  dlnlnQ = min(dy,0.05_dp)
  order = -5
  call dglapStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order)
end subroutine dglapStart


!======================================================================
!! initialise the underlying grid, splitting functions and pdf-table
!! objects, using an extended set of parameters, as described below
subroutine dglapStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order)
  use f77_pdftab
  implicit none
  real(dp), intent(in) :: ymax   !! highest value of ln1/x user wants to access
  real(dp), intent(in) :: dy     !! internal grid spacing: 0.1 is a sensible value
  real(dp), intent(in) :: Qmin, Qmax !! range in Q
  real(dp), intent(in) :: dlnlnQ !! internal table spacing in lnlnQ
  integer,  intent(in) :: nloop  !! the maximum number of loops we'll want (<=3)
  integer,  intent(in) :: order  !! order of numerical interpolation (+ve v. -ve: see below)
  !-------------------------------------

  ! initialise our grids

  ! the internal interpolation order (with a minus sign allows
  ! interpolation to take fake zero points beyond x=1 -- convolution
  ! times are unchanged, initialisation time is much reduced and
  ! accuracy is slightly reduced)
  !order = -5 
  ! Now create a nested grid
  call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
  call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  call InitGridDef(grid,gdarray(1:3),locked=.true.)

  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.
  call AllocPdfTable(grid, tables(:), Qmin, Qmax, & 
       & dlnlnQ = dlnlnQ, freeze_at_Qmin=.true.)

  ! initialise splitting-function holder
  call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop,nflo=3,nfhi=6)
  ! choose a sensible default number of flavours.
  call SetNfDglapHolder(dh,nflcl=5)

  ! indicate the pdfs and convolutions have not been initialised...
  setup_done = .false.
end subroutine dglapStartExtended


!======================================================================
!! Given a pdf_subroutine with the interface shown below, initialise
!! our internal pdf table.
subroutine dglapAssign(pdf_subroutine)
  use f77_pdftab ! this module which provides access to the array of tables
  implicit none
  interface ! indicate what "interface" pdf_subroutine is expected to have
     subroutine pdf_subroutine(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine pdf_subroutine
  end interface
  !-----------------------------------

  ! set up table(0) by copying the values returned by pdf_subroutine onto 
  ! the x-Q grid in table(0)
  call FillPdfTable_LHAPDF(tables(0), pdf_subroutine)
  ! indicate that table(0) has been set up
  setup_done(0)  = .true.
  ! indicate that table(1), table(2), etc... (which will contain the
  ! convolutions with splitting matrices) have yet to be set up [they
  ! get set up "on demand" later].
  setup_done(1:) = .false.
end subroutine dglapAssign


!======================================================================
!! Given a pdf_subroutine with the interface shown below, initialise
!! our internal pdf table.
subroutine dglapEvolve(asmz, nloop, pdf_subroutine, Q0)
  use f77_pdftab ! this module which provides access to the array of tables
  use qcd_coupling
  implicit none
  real(dp), intent(in) :: asmz, Q0
  integer,  intent(in) :: nloop
  interface ! indicate what "interface" pdf_subroutine is expected to have
     subroutine pdf_subroutine(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine pdf_subroutine
  end interface
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)
  !! hold the running coupling
  type(running_coupling) :: coupling

  ! create our internal pdf object for the initial condition
  call AllocPDF(grid, pdf0)
  call InitPDF_LHAPDF(grid, pdf0, pdf_subroutine, Q0)

  ! get a running coupling with the desired scale
  call InitRunningCoupling(coupling, alfas=0.118_dp, nloop=nloop)

  ! create the tabulation
  call EvolvePdfTable(tables(0), Q0, pdf0, dh, coupling, nloop=nloop)

  ! indicate that table(0) has been set up
  setup_done(0)  = .true.
  ! indicate that table(1), table(2), etc... (which will contain the
  ! convolutions with splitting matrices) have yet to be set up [they
  ! get set up "on demand" later].
  setup_done(1:) = .false.

  ! clean up
  call Delete(pdf0)
  call Delete(coupling)
end subroutine dglapEvolve


!======================================================================
!! Return in f(-6:6) the value of the internally stored pdf at the
!! given x,Q, with the usual LHApdf meanings for the indices -6:6.
subroutine dglapEval(x,Q,f)
  use f77_pdftab
  implicit none
  real(dp), intent(in)  :: x, Q
  real(dp), intent(out) :: f(-6:6)
  
  call EvalPdfTable_xQ(tables(0),x,Q,f)
end subroutine dglapEval



!======================================================================
!! Return in f(-6:6) the value of 
!!
!!    [P(iloop,nf) \otimes pdf] (x,Q)
!!
!! where P(iloop,nf) is the iloop-splitting function for the given
!! value of nf, and pdf is our internally stored pdf.
!!
!! The normalisation is such that the nloop dglap evolution equation is
!!
!!     dpdf/dlnQ^2 = sum_{iloop=1}^nloop 
!!                        (alphas/(2*pi))^iloop * P(iloop,nf) \otimes pdf
!!
!! Note that each time nf changes relative to a previous call for the
!! same iloop, the convolution has to be repeated for the whole
!! table. So for efficient results when requiring multiple nf values,
!! calls with the same nf value should be grouped together.
!!
!! In particular, for repeated calls with the same value of nf, the
!! convolutions are carried out only on the first call (i.e. once for
!! each value of iloop). Multiple calls with different values for
!! iloop can be carried out without problems.
!!
subroutine dglapEvalSplit(x,Q,iloop,nf,f)
  use f77_pdftab; use warnings_and_errors
  implicit none
  real(dp), intent(in)  :: x, Q
  integer,  intent(in)  :: iloop, nf
  real(dp), intent(out) :: f(-6:6)
  
  if (.not. setup_done(iloop) .or. setup_nf(iloop) /= nf) then
     if (iloop > size(dh%allP,dim=1) .or. iloop < 1) &
          &call wae_error('dglapeval_split','illegal value for iloop:',&
          &intval=iloop)

     if (nf < lbound(dh%allP,dim=2) .or. nf > ubound(dh%allP,dim=2)) &
          &call wae_error('dglapeval_split','illegal value for nf:',&
          &intval=nf)

     tables(iloop)%tab = dh%allP(iloop, nf) .conv. tables(0)%tab
     
     setup_done(iloop) = .true.
     setup_nf(iloop)   = nf
  end if
  
  call EvalPdfTable_xQ(tables(iloop),x,Q,f)
end subroutine dglapEvalSplit


