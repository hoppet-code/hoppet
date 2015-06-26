module streamlined_interface
  use types; use consts_dp
  use pdf_tabulate
  use convolution; use pdf_general; use dglap_objects
  use dglap_holders; use pdf_general; use dglap_choices
  use qcd_coupling
  use qcd, only: quark_masses_def
  implicit none

  !! holds information about the grid
  type(grid_def),     save :: grid, gdarray(4)

  !! holds the splitting functions
  type(dglap_holder), save :: dh

  !! table_index_from_iloop tells us which table to loop up for a give
  !! iloop value; most of the iloop values are illegal (index -1), but
  !! for those that are allowed (e.g. 1, 2, 3, 12, 111) this gives
  !! quick access to the relevant table (without having to store 111
  !! tables, most of which would be redundant)
  !!
  !! The array gets initialised in hoppetStartExtended
  integer, save :: table_index_from_iloop(1:111) = -1
  integer, parameter :: max_table_index = 7

  !!
  !! NB (2012-12-24): 
  !!     in future evolution of the code, we should aim to guarantee
  !!     that tables(0) remains the default PDF, even if the other
  !!     entries change
  type(pdf_table), save :: tables(0:max_table_index)
  logical,         save :: setup_done(0:max_table_index) = .false.
  integer,         save :: setup_nf(max_table_index)     = 0
  logical,         save :: alloc_already_done = .false.
  
  !! coupling
  logical,                save :: coupling_initialised = .false.
  type(running_coupling), save :: coupling
  integer,  save :: ffn_nf = -1
  logical,  save :: quark_masses_are_MSbar = .false.
  real(dp), save :: masses(4:6) = quark_masses_def(4:6)

contains

  !----------------------------------------------------------------------
  ! returns the value of the table index for the given iloop and
  ! ensures that the table is actually set up (working down recursively
  ! when iloop is composite)
  recursive integer function tableIndexValue(iloop, nf) result(tabindex)
    use warnings_and_errors
    integer, intent(in) :: iloop, nf
    !------------------------------
    integer :: table_index_source, iloop_P, iQ

    ! first establish if iloop is "legal"
    if (iloop < 1 .or. iloop > ubound(table_index_from_iloop,dim=1)) then
       call wae_error('tableIndexValue','illegal value for iloop (out of bounds):',&
            &intval=iloop)
    else if (table_index_from_iloop(iloop) < 0) then
       call wae_error('tableIndexValue','illegal value for iloop (lookup < 0):',&
            &intval=iloop)
    end if
    tabindex = table_index_from_iloop(iloop)

    !write(0,*) 'iloop, tabindex ', iloop, tabindex
    if (.not. setup_done(tabindex) .or. setup_nf(tabindex) /= nf) then

       if (iloop < 10) then
          table_index_source = 0
          iloop_P = iloop
       else if (iloop < 100) then
          table_index_source = tableIndexValue(mod(iloop,10), nf)
          iloop_P = iloop / 10
       else if (iloop < 1000) then
          table_index_source = tableIndexValue(mod(iloop,100), nf)
          iloop_P = iloop / 100
       else 
          call wae_error('tableIndexValue','unsupported value for iloop (source & iloop_P determination):',&
            &intval=iloop)
       end if
       !write(0,*) 'setting up with source and iloop_P = ', table_index_source, iloop_P
       
       if (nf < 0) then
          ! use whatever nf is relevant 
          if (.not. tables(0)%nf_info_associated) call wae_error(&
               & 'hoppetEvalSplit','automatic choice of nf (nf<0) but the tabulation has no information on nf')
          do iQ = 0, tables(0)%nQ
             tables(tabindex)%tab(:,:,iQ) = dh%allP(iloop_P, tables(0)%nf_int(iQ)) &
                  &             .conv. tables(table_index_source)%tab(:,:,iQ)
          end do
       else
          ! use a fixed nf
          if (nf < lbound(dh%allP,dim=2) .or. nf > ubound(dh%allP,dim=2)) &
               &call wae_error('hoppetEvalSplit','illegal value for nf:',&
               &intval=nf)
          
          tables(tabindex)%tab = dh%allP(iloop_P, nf) .conv. tables(table_index_source)%tab
       end if
       
       setup_done(tabindex) = .true.
       setup_nf(tabindex)   = nf
  end if

  end function tableIndexValue

  
  
end module streamlined_interface


!======================================================================
!! initialise the underlying grid, splitting functions and pdf-table
!! objects, using the dy and nloop parameters as explained below.
subroutine hoppetStart(dy,nloop)
  use streamlined_interface
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
  dlnlnQ = dy/4.0_dp  ! min(0.5_dp*dy,0.07_dp)
  order = -6
  call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,factscheme_MSbar)
end subroutine hoppetStart


!======================================================================
!! initialise the underlying grid, splitting functions and pdf-table
!! objects, using an extended set of parameters, as described below
subroutine hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,factscheme)
  use streamlined_interface
  implicit none
  real(dp), intent(in) :: ymax   !! highest value of ln1/x user wants to access
  real(dp), intent(in) :: dy     !! internal grid spacing: 0.1 is a sensible value
  real(dp), intent(in) :: Qmin, Qmax !! range in Q
  real(dp), intent(in) :: dlnlnQ !! internal table spacing in lnlnQ
  integer,  intent(in) :: nloop  !! the maximum number of loops we'll want (<=3)
  integer,  intent(in) :: order  !! order of numerical interpolation (+ve v. -ve: see below)
  integer,  intent(in) :: factscheme !! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar
  !-------------------------------------
  integer :: iloop
  ! initialise our grids

  ! the internal interpolation order (with a minus sign allows
  ! interpolation to take fake zero points beyond x=1 -- convolution
  ! times are unchanged, initialisation time is much reduced and
  ! accuracy is slightly reduced)
  !order = -5 
  ! Now create a nested grid
  call InitGridDef(gdarray(4),dy/27.0_dp,0.2_dp, order=order)
  call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
  call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  call InitGridDef(grid,gdarray(1:4),locked=.true.)

  ! Fill the array that will be used for table index lookup (e.g. 21 is PNLO*PLO).
  ! For now do it by hand; one day we might automate this;
  ! entries that aren't filled are automatically -1
  do iloop = 1, nloop
     table_index_from_iloop(iloop) = iloop
  end do
  table_index_from_iloop(11)  = 4
  table_index_from_iloop(111) = 5
  if (nloop >= 2) table_index_from_iloop(12)  = 6
  if (nloop >= 2) table_index_from_iloop(21)  = 7

  ! if the allocation has already been done previously, delete
  ! the existing tables and dglap holder to avoid a memory leak
  if (alloc_already_done) then
     call Delete(tables)
     call Delete(dh)
  end if
  

  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.
  call AllocPdfTable(grid, tables(:), Qmin, Qmax, & 
       & dlnlnQ = dlnlnQ, freeze_at_Qmin=.true.)

  ! initialise splitting-function holder
  call InitDglapHolder(grid,dh,factscheme=factscheme,&
       &                      nloop=nloop,nflo=3,nfhi=6)
  ! choose a sensible default number of flavours.
  call SetNfDglapHolder(dh,nflcl=5)

  ! indicate that allocations have already been done,
  ! to allow for cleanup if hoppetStartExtended is
  ! called a second time
  alloc_already_done = .true.

  ! indicate the pdfs and convolutions have not been initialised...
  setup_done = .false.
end subroutine hoppetStartExtended


!======================================================================
!! Given a pdf_subroutine with the interface shown below, initialise
!! our internal pdf table.
subroutine hoppetAssign(pdf_subroutine)
  use streamlined_interface ! this module which provides access to the array of tables
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
end subroutine hoppetAssign


!======================================================================
!! Given a pdf_subroutine with the interface shown below, fill the 
!! table by evolving the PDF from scale Q0pdf, with alphas provided 
!! at scale Q0alphas
subroutine hoppetEvolve(asQ0, Q0alphas, nloop, muR_Q, pdf_subroutine, Q0pdf)
  use streamlined_interface ! this module which provides access to the array of tables implicit none
  implicit none
  real(dp), intent(in) :: asQ0, Q0alphas, muR_Q, Q0pdf
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

  ! create our internal pdf object for the initial condition
  call AllocPDF(grid, pdf0)
  call InitPDF_LHAPDF(grid, pdf0, pdf_subroutine, Q0pdf)

  ! get a running coupling with the desired scale
  if (coupling_initialised) call Delete(coupling) 
  if (ffn_nf > 0) then
     call InitRunningCoupling(coupling, alfas=asQ0, Q=Q0alphas, nloop=nloop, &
          &                   fixnf=ffn_nf)
  else 
     call InitRunningCoupling(coupling, alfas=asQ0, Q=Q0alphas, nloop=nloop, &
          &                   quark_masses=masses, &
          &                   masses_are_MSbar = quark_masses_are_MSbar)
  end if
  call AddNfInfoToPdfTable(tables,coupling)
  coupling_initialised = .true.

  ! create the tabulation
  call EvolvePdfTable(tables(0), Q0pdf, pdf0, dh, muR_Q=muR_Q, &
       &              coupling=coupling, nloop=nloop)

  ! indicate that table(0) has been set up
  setup_done(0)  = .true.
  ! indicate that table(1), table(2), etc... (which will contain the
  ! convolutions with splitting matrices) have yet to be set up [they
  ! get set up "on demand" later].
  setup_done(1:) = .false.

  ! clean up
  call Delete(pdf0)
end subroutine hoppetEvolve


!======================================================================
!! Prepare a cached evolution
subroutine hoppetPreEvolve(asQ0, Q0alphas, nloop, muR_Q, Q0pdf)
  use streamlined_interface
  implicit none
  real(dp), intent(in) :: asQ0, Q0alphas, muR_Q, Q0pdf
  integer,  intent(in) :: nloop

  ! get a running coupling with the desired scale
  if (coupling_initialised) call Delete(coupling) 
  if (ffn_nf > 0) then
     call InitRunningCoupling(coupling, alfas=asQ0, Q=Q0alphas, nloop=nloop, &
          &                   fixnf=ffn_nf)
  else 
     call InitRunningCoupling(coupling, alfas=asQ0, Q=Q0alphas, nloop=nloop, &
          &                   quark_masses=masses)
  end if
  call AddNfInfoToPdfTable(tables,coupling)
  coupling_initialised = .true.

  ! create the tabulation
  call PreEvolvePdfTable(tables(0), Q0pdf, dh, muR_Q=muR_Q, &
       &                 coupling=coupling, nloop=nloop)
end subroutine hoppetPreEvolve


!======================================================================
!! Carry out a cached evolution based on the initial condition
!! that can be obtained from pdf_subroutine at the scale Q0pdf set in
!! PreEvolve
subroutine hoppetCachedEvolve(pdf_subroutine)
  use streamlined_interface
  implicit none
  interface ! indicate what "interface" pdf_subroutine is expected to have
     subroutine pdf_subroutine(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine pdf_subroutine
  end interface
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)

  ! create our internal pdf object for the initial condition
  call AllocPDF(grid, pdf0)
  call InitPDF_LHAPDF(grid, pdf0, pdf_subroutine, tables(0)%StartScale)

  ! create the tabulation
  call EvolvePdfTable(tables(0), pdf0)

  ! indicate that table(0) has been set up
  setup_done(0)  = .true.
  ! indicate that table(1), table(2), etc... (which will contain the
  ! convolutions with splitting matrices) have yet to be set up [they
  ! get set up "on demand" later].
  setup_done(1:) = .false.

  ! clean up
  call Delete(pdf0)
end subroutine hoppetCachedEvolve


!======================================================================
!! Return the coupling at scale Q
function hoppetAlphaS(Q)
  use streamlined_interface ! this module which provides access to the array of tables
  use warnings_and_errors
  implicit none
  real(dp)             :: hoppetAlphaS
  real(dp), intent(in) :: Q

  if (.not. coupling_initialised) call wae_error('hoppetAlphaS',&
       &'coupling is not yet initialised (and will remain so until',&
       &'first call to an evolution routine).')
  hoppetAlphaS = Value(coupling, Q)
end function hoppetAlphaS


!======================================================================
!! Set up things to be a fixed-flavour number scheme with the given
!! fixed_nf number of flavours
subroutine hoppetSetFFN(fixed_nf)
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  integer,  intent(in) :: fixed_nf

  ffn_nf = fixed_nf
end subroutine hoppetSetFFN


!======================================================================
!! Set up things to be a variable-flavour number scheme with the given
!! quark masses in the pole mass scheme.
!!
!! This interface retained for legacy purposes. Preferred interface is
!! via hoppetSetPoleMassVFN(mc,mb,mt)
subroutine hoppetSetVFN(mc,mb,mt)
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  real(dp) :: mc, mb, mt
  call hoppetSetPoleMassVFN(mc,mb,mt)
end subroutine hoppetSetVFN


!======================================================================
!! Set up things to be a variable-flavour number scheme with the given
!! quark masses in the pole mass scheme. Thresholds are crossed at the
!! pole masses, both for the coupling and the PDF evolution.
subroutine hoppetSetPoleMassVFN(mc,mb,mt)
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  real(dp) :: mc, mb, mt

  ffn_nf = -1
  masses = (/mc,mb,mt/)
  quark_masses_are_MSbar = .false.
end subroutine hoppetSetPoleMassVFN


!======================================================================
!! Set up things to be a variable-flavour number scheme with the given
!! quark masses in the MSbar mass scheme, i.e. m(m). Thresholds are
!! crossed at the MSbar masses, both for the coupling and the PDF
!! evolution.
subroutine hoppetSetMSbarMassVFN(mc,mb,mt)
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  real(dp) :: mc, mb, mt

  ffn_nf = -1
  masses = (/mc,mb,mt/)
  quark_masses_are_MSbar = .true.
end subroutine hoppetSetMSbarMassVFN


!======================================================================
!! Return in f(-6:6) the value of the internally stored pdf at the
!! given x,Q, with the usual LHApdf meanings for the indices -6:6.
subroutine hoppetEval(x,Q,f)
  use streamlined_interface
  implicit none
  real(dp), intent(in)  :: x, Q
  real(dp), intent(out) :: f(-6:6)
  
  call EvalPdfTable_xQ(tables(0),x,Q,f)
end subroutine hoppetEval



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
!! Note that iloop can also be of the form ij or ijk, which means
!! P(i)*P(j)*pdf or P(i)*P(j)*P(k)*pdf. The sum of i+j+k is currently
!! bounded to be <= 3.
!!
!! The number of loops must be consistent with iloop
subroutine hoppetEvalSplit(x,Q,iloop,nf,f)
  use streamlined_interface; use warnings_and_errors
  implicit none
  real(dp), intent(in)  :: x, Q
  integer,  intent(in)  :: iloop, nf
  real(dp), intent(out) :: f(-6:6)
  integer :: iQ, tabindex

  tabindex = tableIndexValue(iloop, nf)
  call EvalPdfTable_xQ(tables(tabindex),x,Q,f)

end subroutine hoppetEvalSplit


