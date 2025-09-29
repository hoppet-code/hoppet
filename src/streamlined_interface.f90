module streamlined_interface
  use types; use consts_dp
  use pdf_tabulate
  use convolution; use pdf_general; use dglap_objects
  use dglap_holders; use pdf_general; use dglap_choices
  use qcd_coupling
  use qed_evolution
  use qed_objects
  use qed_coupling_module
  use qcd, only: quark_masses_def
  use pdfs_for_benchmarks
  use warnings_and_errors
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
  integer, save :: table_index_from_iloop(1:1111) = -1
  integer, parameter :: max_table_index = 15

  !!
  !! NB (2012-12-24): 
  !!     in future evolution of the code, we should aim to guarantee
  !!     that tables(0) remains the default PDF, even if the other
  !!     entries change
  type(pdf_table), save, target :: tables(0:max_table_index)
  logical,         save :: setup_done(0:max_table_index) = .false.
  integer,         save :: setup_nf(max_table_index)     = 0
  logical,         save :: alloc_already_done = .false.
  
  !! coupling
  logical,                save :: coupling_initialised = .false.
  type(running_coupling), save :: coupling
  integer,  save :: ffn_nf = -1
  logical,  save :: quark_masses_are_MSbar = .false.
  real(dp), save :: masses(4:6) = quark_masses_def(4:6)
  !! Additional things for QED
  !integer,  save            :: nqcdloop_qed
  type(qed_coupling),  save :: coupling_qed
  type(qed_split_mat), save :: qed_split
  logical,  save :: with_qed = .false.
  real(dp), save :: effective_light_quark_masses = 0.109_dp
  integer,  save :: with_nqcdloop_qed = 0
  logical,  save :: with_Plq_nnloqed = .false.

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
       else if (iloop < 10000) then
          table_index_source = tableIndexValue(mod(iloop,1000), nf)
          iloop_P = iloop / 1000
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
!! This routine should be called before hoppetStart if one wants
!! to include QED evolution. Note that it affects expected upper
!! bounds for all routines that set and access the PDFs
!!
!! - use_qed: if true, QED evolution will be included
!! - use_qcd_qed: if true, the mixed alpha_s*alpha splitting functions will be included
!! - use_Plq_nnlo: if true, the NNLO Plq splitting functions will be included
subroutine hoppetSetQED(use_qed, use_qcd_qed, use_Plq_nnlo)
  use streamlined_interface
  use qed_evolution
  implicit none
  logical, intent(in)  :: use_qed
  logical, intent(in)  :: use_qcd_qed, use_Plq_nnlo
  write(*,*) 'hoppetSetQED: use_qed, use_qcd_qed, use_Plq_nnlo=', &
       & use_qed, use_qcd_qed, use_Plq_nnlo
  with_qed = use_qed
  if (use_qcd_qed) then 
    with_nqcdloop_qed = 1
  else
    with_nqcdloop_qed = 0
  end if
  with_Plq_nnloqed = use_Plq_nnlo
end subroutine hoppetSetQED

!======================================================================
!! initialise the underlying grid, splitting functions and pdf-table
!! objects, using the dy and nloop parameters as explained below.
subroutine hoppetStart(dy,nloop)
  use streamlined_interface
  implicit none
  !--------------------------------------
  real(dp), intent(in) :: dy     !! internal grid spacing: 0.1 is a sensible value
  integer,  intent(in) :: nloop  !! the maximum number of loops we'll want (<=4)
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
  integer,  intent(in) :: nloop  !! the maximum number of loops we'll want (<=4)
  integer,  intent(in) :: order  !! order of numerical interpolation (+ve v. -ve: see below)
  integer,  intent(in) :: factscheme !! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar
  !-------------------------------------
  integer :: iloop

  ! if the allocation has already been done previously, delete
  ! the existing tables and dglap holder, etc. to avoid a memory leak
  if (alloc_already_done) call hoppetDeleteAll()  
  
  ! initialise our grids

  ! the internal interpolation order (with a minus sign allows
  ! interpolation to take fake zero points beyond x=1 -- convolution
  ! times are unchanged, initialisation time is much reduced and
  ! accuracy is slightly reduced)
  !order = -5 
  ! Now create a nested grid
  !call InitGridDef(gdarray(4),dy/27.0_dp,0.2_dp, order=order)
  !call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
  !call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  !call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  !call InitGridDef(grid,gdarray(1:4),locked=.true.)

  ! as of 2016-05-26, better to use the "central" facility for
  ! producing a nested grid (which has the same nested grid structure
  ! as above)
  call InitGridDefDefault(grid, dy, ymax, order)
  
  ! Fill the array that will be used for table index lookup (e.g. 21 is PNLO*PLO).
  ! For now do it by hand; one day we might automate this;
  ! entries that aren't filled are automatically -1

  ! First fill the bare splitting functions
  do iloop = 1, nloop
     table_index_from_iloop(iloop) = iloop
  end do
  ! Now fill entries where the sum is 2
  table_index_from_iloop(11)                   = nloop + 1
  ! Then the entries where the sum is 3
  table_index_from_iloop(111)                  = nloop + 2
  if (nloop >= 2) table_index_from_iloop(12)   = nloop + 3
  if (nloop >= 2) table_index_from_iloop(21)   = nloop + 4
  ! Finally the entries where the sum is 4
  table_index_from_iloop(1111)                 = nloop + 5
  if (nloop >= 2) table_index_from_iloop(112)  = nloop + 6
  if (nloop >= 2) table_index_from_iloop(121)  = nloop + 7
  if (nloop >= 2) table_index_from_iloop(211)  = nloop + 8
  if (nloop >= 2) table_index_from_iloop(22)   = nloop + 9
  if (nloop >= 3) table_index_from_iloop(13)   = nloop + 10
  if (nloop >= 3) table_index_from_iloop(31)   = nloop + 11 ! which is 15 if nloop = 14
  
  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.

  if(with_qed) then
     call AllocPdfTableWithLeptons(grid, tables(:), Qmin, Qmax, & 
          & dlnlnQ = dlnlnQ, freeze_at_Qmin=.true.)
  else
     call AllocPdfTable(grid, tables(:), Qmin, Qmax, & 
          & dlnlnQ = dlnlnQ, freeze_at_Qmin=.true.)
  endif

  ! initialise splitting-function holder
  call InitDglapHolder(grid,dh,factscheme=factscheme,&
       &                      nloop=nloop,nflo=3,nfhi=6)

  if(with_qed) then
     call InitQEDSplitMat(grid, qed_split)
  endif
  
  ! choose a sensible default number of flavours.
  call SetNfDglapHolder(dh,nflcl=5)

  ! indicate that allocations have already been done,
  ! to allow for cleanup if hoppetStartExtended is
  ! called a second time
  alloc_already_done = .true.

  ! indicate the pdfs and convolutions have not been initialised...
  setup_done = .false.
end subroutine hoppetStartExtended


!! Delete all hoppet objects associated with the streamlined interface,
!! including grid definitions, PDF tables, couplings, etc.
!!
!! NB: this does not delete anything associated with the structure function
!! part of the interface
subroutine hoppetDeleteAll()
  use streamlined_interface
  implicit none
  if (alloc_already_done) then
    call Delete(tables)
    call Delete(dh)
    if (with_qed) call Delete(qed_split)
    call Delete(grid)
    if (coupling_initialised) then
      call Delete(coupling) 
      if (with_qed) call Delete(coupling_qed)
      coupling_initialised = .false.
    end if
    alloc_already_done = .false.
  end if
end subroutine hoppetDeleteAll

!======================================================================
!! Given a pdf_subroutine with the interface shown below, initialise
!! our internal pdf table.
subroutine hoppetAssign(pdf_subroutine)
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  interface ! indicate what "interface" pdf_subroutine is expected to have
    !! It expects pdf_subroutine to set components -6:6 with plain QCD evolution,
	  !! but the full size of the flavour dimension if the upper limit
	  !! is anything other than ncompax=7 (e.g. for QED evolution 
	  !! it should go from -6 to ncompmaxLeptons=11).
    subroutine pdf_subroutine(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine pdf_subroutine
  end interface
  !-----------------------------------

  if (ffn_nf < 0 .and. (.not. coupling_initialised)) then
    call wae_error('hoppetAssign','Since v2 of hoppet, if you assign a PDF with hoppet''s VFN option turned on (default) '//&
    'you must first set the quark masses (hoppetSetPoleMassVFN(...)) and then the coupling (hoppetSetCoupling(...)). ')
  end if

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
!! Set up the strong coupling such that alphas(Q)=alphas_Q, with the
!! given number of loops (nloop).
!!
!! The user should have set the quark masses or requested a FFN scheme
!! prior to calling this function.
!!
!! This function is provided mainly for use in conjunction with hoppetAssign.
!! In particular, it has the side effect of modifying the structure of the
!! PDF tables to make sure they know about the mass thresholds.
!!
!! If QED has been requested, a QED coupling will also be set up
!! (its value is not currently configurable from this interface).
!!
!! If you call hoppetEvolve (below), there is no need to separately call
!! hoppetSetCoupling.
subroutine hoppetSetCoupling(alphas_Q,Q,nloop)
  use streamlined_interface
  implicit none
  real(dp), intent(in) :: alphas_Q, Q
  integer, intent(in) :: nloop

  ! get a running coupling with the desired scale
  if (coupling_initialised) then
    call Delete(coupling) 
    if (with_qed) call Delete(coupling_qed)
  end if
  if (ffn_nf > 0) then
     call InitRunningCoupling(coupling, alfas=alphas_Q, Q=Q, nloop=nloop, &
          &                   fixnf=ffn_nf)
  else 
     call InitRunningCoupling(coupling, alfas=alphas_Q, Q=Q, nloop=nloop, &
          &                   quark_masses=masses, &
          &                   masses_are_MSbar = quark_masses_are_MSbar)
  end if
  if(with_qed) then
     call InitQEDCoupling(coupling_qed, effective_light_quark_masses, masses(4:6))
  endif
  call AddNfInfoToPdfTable(tables,coupling)
  coupling_initialised = .true.

end subroutine hoppetSetCoupling

!======================================================================
!! Given a pdf_subroutine with the interface shown below, fill the 
!! table by evolving the PDF from scale Q0pdf, with alphas provided 
!! at scale Q0alphas
subroutine hoppetEvolve(asQ0, Q0alphas, nloop,  muR_Q, pdf_subroutine, Q0pdf)
  use streamlined_interface ! this module which provides access to the array of tables
  use pdf_representation
  implicit none
  real(dp), intent(in) :: asQ0, Q0alphas, muR_Q, Q0pdf
  integer,  intent(in) :: nloop
  interface ! indicate what "interface" pdf_subroutine is expected to have
    !! It expects pdf_subroutine to set components -6:6 with plain QCD evolution,
	  !! but the full size of the flavour dimension if the upper limit
	  !! is anything other than ncompax=7 (e.g. for QED evolution 
	  !! it should go from -6 to ncompmaxLeptons=11).
    subroutine pdf_subroutine(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine pdf_subroutine
  end interface
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)
  ! create our internal pdf object for the initial condition
  if(with_qed) then
     !write(*,*)'***********HoppetEvolve, before AllocPDFWithLeptons'
     call AllocPDFWithLeptons(grid, pdf0)
  else
     call AllocPDF(grid, pdf0)
  endif

  !write(*,*)'***********HoppetEvolve, before InitPDF_LHAPDF'
  call InitPDF_LHAPDF(grid, pdf0, pdf_subroutine, Q0pdf)

  call hoppetSetCoupling(asQ0,Q0alphas,nloop)

  ! create the tabulation
  if(with_qed) then
     if(muR_Q /= 1) then
        call wae_error("hoppetEvolve", "muR_Q /= 1 not allowed if qed is included")
     endif
     !write(*,*)'***********HoppetEvolve, before EvolvePdfTableQED'
     call EvolvePdfTableQED(tables(0), Q0pdf, pdf0, dh, qed_split, &
          coupling, coupling_qed, nloop, with_nqcdloop_qed, with_Plq_nnloqed)
  else
     call EvolvePdfTable(tables(0), Q0pdf, pdf0, dh, muR_Q=muR_Q, &
          &              coupling=coupling, nloop=nloop)
  endif

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

  if (with_qed) call wae_error("hoppetPreEvolve","Does not support QED evolution")

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

  if (with_qed) call wae_error("hoppetCachedEvolve","Does not support QED evolution")

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
!! Arrange for the use of exact NNLO splitting and mass-threshold
!! functions.
subroutine hoppetSetExactDGLAP(exact_nfthreshold, exact_splitting)
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  logical :: exact_nfthreshold, exact_splitting

  if(exact_nfthreshold) call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact)
  if(exact_splitting) call dglap_Set_nnlo_splitting(nnlo_splitting_exact)

end subroutine hoppetSetExactDGLAP

!======================================================================
!! Arrange for the use of various approximate N3LO splitting functions.
subroutine hoppetSetApproximateDGLAPN3LO(splitting_variant)
  ! splitting_variant can be one of
  ! n3lo_splitting_approximation_up_to_2310_05744
  ! n3lo_splitting_approximation_up_to_2404_09701 
  ! n3lo_splitting_approximation_up_to_2410_08089 (default)
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  integer :: splitting_variant

  call dglap_Set_n3lo_splitting_approximation(splitting_variant)

end subroutine hoppetSetApproximateDGLAPN3LO
!======================================================================

!! Arrange for the use of N3LO splitting functions, either exact (not currently available), parametrised not currently available) or based on upper or lower values of the aproximation
subroutine hoppetSetSplittingN3LO(splitting_variant)
  ! splitting_variant can be one of
  !n3lo_splitting_exact = -2;
  !n3lo_splitting_param = -1;
  !n3lo_splitting_Nfitav   =  0;
  !n3lo_splitting_Nfiterr1 =  1;
  !n3lo_splitting_Nfiterr2 =  2;
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  integer :: splitting_variant

  call dglap_Set_n3lo_splitting(splitting_variant)

end subroutine hoppetSetSplittingN3LO

!! Arrange for the use of N3LO splitting functions, either exact (not currently available), parametrised not currently available) or based on upper or lower values of the aproximation
subroutine hoppetSetSplittingNNLO(splitting_variant)
  ! splitting_variant can be one of
  !nnlo_splitting_exact = -2;
  !nnlo_splitting_param = -1;
  !nnlo_splitting_Nfitav   =  0;
  !nnlo_splitting_Nfiterr1 =  1;
  !nnlo_splitting_Nfiterr2 =  2;
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  integer :: splitting_variant

  call dglap_Set_nnlo_splitting(splitting_variant)

end subroutine hoppetSetSplittingNNLO


!======================================================================
!! Arrange for the use of N3LO mass thresholds or not.
subroutine hoppetSetN3LOnfthresholds(variant)
  ! variant can be one of
  ! integer, parameter, public :: n3lo_nfthreshold_on  = 1
  ! integer, parameter, public :: n3lo_nfthreshold_off = 0
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  integer :: variant

  call dglap_Set_n3lo_nfthreshold(variant)

end subroutine hoppetSetN3LOnfthresholds

!======================================================================
!! Override the default interpolation order in y and lnlnQ. Should be
!! called before hoppetStart
subroutine hoppetSetYLnlnQInterpOrders(yorder, lnlnQorder)
  ! Defaults are
  ! yorder = -1 (hoppet makes a sensible choice in the ranges 4 to 10)
  ! lnlnQorder = 4
  ! For fast cubic interpolation with similar accuracy to LHAPDF one can set
  ! yorder = 3
  ! lnlnQorder = 3
  use streamlined_interface ! this module which provides access to the array of tables
  implicit none
  integer :: yorder, lnlnQorder

  call PdfTableOverrideInterpOrders(yorder, lnlnQorder)
  !call PDFTableSetYInterpOrder(yorder)
  def_lnlnQ_order = lnlnQorder

end subroutine hoppetSetYLnlnQInterpOrders




!======================================================================
!! Assuming this is called with an array starting at index -6, return in f(-6:6) 
!! the value of the internally stored pdf at the
!! given x,Q, with the usual LHApdf meanings for the indices -6:6.
!!
!! If QED has been enabled, the indices that will be set are 
!! f(-6:ncompmaxLeptons), where ncompmaxLeptons=11.
!! NB: f(7) is a dummy entry, f(8) is the photon and f(9:11) are 
!! of (e+ + e-), (mu+ + mu-) and (tau+ + tau-)
subroutine hoppetEval(x,Q,f)
  use streamlined_interface; use pdf_representation
  implicit none
  real(dp), intent(in)  :: x, Q
  ! the interface does pass the size of the array, but the functions we
  ! call have interfaces that do need the size; so here give it a dummy
  ! value that is the largest that could ever be used. The functions
  ! that we call will fill -6:6 with normal QCD evolution set and
  ! -6:ncompmaxLeptons with QED evolution set.
  !real(dp), intent(out) :: f(-6:ncompmaxLeptons) ! QED-TBD [check it still works without QED]
  !call EvalPdfTable_xQ(tables(0),x,Q,f)

  real(dp), intent(out) :: f(*) ! QED-TBD [check it still works without QED]
  call EvalPdfTable_xQ(tables(0),x,Q,f(1:tables(0)%tab_iflv_max-ncompmin+1))
end subroutine hoppetEval

!======================================================================
!! Assuming this is called with an array starting at index -6, return in f(-6:6) 
!! the value of the internally stored pdf at the
!! given x,Q, and for the hoppet PID
function hoppetEvalPID(x,Q,pid) 
  use streamlined_interface; use pdf_representation
  implicit none
  real(dp), intent(in)  :: x, Q
  integer, intent(in) :: pid
  real(dp) :: hoppetEvalPID

  hoppetEvalPID = EvalPdfTable_xQf(tables(0),x,Q,pid)
end function hoppetEvalPID



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
!! bounded to be <= 4.
!!
!! The number of loops must be consistent with iloop
subroutine hoppetEvalSplit(x,Q,iloop,nf,f)
  use streamlined_interface; use warnings_and_errors
  implicit none
  real(dp), intent(in)  :: x, Q
  integer,  intent(in)  :: iloop, nf
  real(dp), intent(out) :: f(-6:6)
  integer :: iQ, tabindex

  if (with_qed) call wae_error("hoppetEvalSplit","Does not support QED evolution")

  tabindex = tableIndexValue(iloop, nf)
  call EvalPdfTable_xQ(tables(tabindex),x,Q,f)

end subroutine hoppetEvalSplit

!======================================================================
!! Write out the contents of tables(0) (assumed to be the PDF) in the
!! LHAPDF format
!! 
!! - basename: the basename of the LHAPDF file; the full name will be
!!   basename_nnnn.dat, where nnnn is given by pdf_index.
!!   
!! - pdf_index is the index that will be used for the LHAPDF .dat file
!!   If it is zero, the routing will also output the .info file and
!!   and the .dat will be declared a central member; otherwise the
!!   .dat file will be declared an error member.
subroutine hoppetWriteLHAPDFGrid(basename, pdf_index)
  use streamlined_interface; use warnings_and_errors
  implicit none
  character(len=*),       intent(in) :: basename
  integer,                intent(in) :: pdf_index

  if(.not.setup_done(0)) call wae_error("hoppetWriteLHAPDFgrid:"&
       &,"Trying to write out contents of tables(0) from the&
       & streamlined interface, but the table has not been filled.&
       & Did you forget to run hoppetAssign or hoppetEvolve?")

  call WriteLHAPDFFromPdfTable(tables(0), coupling, basename, pdf_index)
  
end subroutine hoppetWriteLHAPDFGrid


subroutine hoppetWriteLHAPDFWithLen(basename_len, basename, pdf_index)
  use streamlined_interface; use warnings_and_errors
  implicit none
  integer,                     intent(in) :: basename_len
  character(len=basename_len), intent(in) :: basename
  integer,                     intent(in) :: pdf_index

  call hoppetWriteLHAPDFGrid(basename, pdf_index)

end subroutine hoppetWriteLHAPDFWithLen

!======================================================================
!! The dummy PDF suggested by Vogt as the initial condition for the 
!! unpolarized evolution, subsequently adopted in a wide range of
!! of benchmarks
subroutine hoppetBenchmarkPDFunpol(x, Q, xpdf)
  use streamlined_interface; use warnings_and_errors
  implicit none
  real(dp) :: x, Q, xpdf(-6:6)
  
  call benchmark_pdf_unpolarized_lha(x, Q, xpdf)

end subroutine hoppetBenchmarkPDFunpol

  
