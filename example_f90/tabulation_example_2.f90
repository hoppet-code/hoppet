!! A second example program using a tabulation. It differs from
!! tabulation_example.f90 in that it shows how to consistently carry
!! out initialisation and evolution from separate sub-programs, as
!! might for example be needed in a PDF fit.
!!
!! It outputs a subset of table 15 of hep-ph/0511119 and this output
!! should be identical to the contents of the file
!! tabulation_example.default_output
!!
!! It also outputs (truncated) sum rules, in order to illustrate
!! how to obtain them from a tabulation
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

!----------------------------------------------------------------------
!! This module holds the bits and pieces needed to store the grid,
!! the dglap_holder, the coupling and a table
!! 
module table_module
  use hoppet_v1
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  !use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  implicit none
  !! holds information about the grid (not strictly needed here, but handy)
  type(grid_def),  save :: grid, gdarray(4)
  !! holds the splitting functions
  type(dglap_holder), save :: dh
  !! hold the PDF tabulation in the module table_module
  type(pdf_table)       :: table
  !! hold the coupling
  type(running_coupling), save :: coupling
end module table_module


!----------------------------------------------------------------------
!! Either place the intial condition subroutine in a module and "use" that
!! module, otherwise define an explicit interface for the subroutine.
!!
!! Here we go for the first option
module initial_condition
  use hoppet_v1 ! 
  implicit none
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
end module initial_condition


!----------------------------------------------------------------------
!! the main program that does the work once the evolution and table 
!! have been initialised
program main
  use table_module
  implicit none
  !! hold results at some x, Q
  real(dp) :: Q, pdf_at_xQ(-6:6)
  real(dp), parameter :: heralhc_xvals(9) = &
       & (/1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: ix, iloop

  ! do the common initialisation (including allocating a table)
  ! (call this just once)
  call tabulation_example_startup

  ! the following might be replaced by the "evolve and calculate
  ! chi^2" loop in a fit program; here we just loop over it once.
  ! NB: if you do plan to run a fit program, then you'd be advised
  !     to use the precalculated evolution (see below)
  do iloop = 1, 1
    ! evolve an initial condition to fill the table
    ! (do this for each initial condition)
    call tabulation_example_evolve_table
  
    ! get the value of the tabulation at some point
    Q = 100.0_dp
    write(6,'(a)')
    write(6,'(a,f8.3,a)') "           Evaluating PDFs at Q = ",Q," GeV"
    write(6,'(a5,2a12,a14,a10,a12)') "x",&
         & "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
    do ix = 1, size(heralhc_xvals)
       call EvalPdfTable_xQ(table,heralhc_xvals(ix),Q,pdf_at_xQ)
       write(6,'(es7.1,5es12.4)') heralhc_xvals(ix), &
            &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
            &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
            &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
            &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
            &  pdf_at_xQ(0)
    end do
  end do

  !! The following routine shows how to print out the sum rules for
  !! the table at the given scale
  call tabulation_print_sum_rules(Q)

  !! deallocate things where relevant (right at end of program,
  !! not strictly needed)
  call tabulation_example_cleanup

end program main


!----------------------------------------------------------------------
!! do the initialisation
subroutine tabulation_example_startup
  use table_module
  implicit none
  real(dp) :: dy, ymax
  integer  :: order, nloop
  real(dp) :: quark_masses(4:6)

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

  !--- if doing many evolutions with different couplings, the
  !    following lines should be placed inside
  !    tabulation_example_evolve_table, and the coupling should be
  !    deleted each time around once you're done with it.
  !
  ! allocate and initialise the running coupling with a given
  ! set of quark masses (NB: charm mass just above Q0).
  quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
  call InitRunningCoupling(coupling,alfas=0.35_dp,Q=sqrt(2.0_dp),nloop=nloop,&
       &                   quark_masses = quark_masses)

  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.
  call AllocPdfTable(grid, table, Qmin=1.0_dp, Qmax=10000.0_dp, & 
       & dlnlnQ = dy/4.0_dp, freeze_at_Qmin=.true.)
  ! add information about the nf transitions to the table (improves
  ! interpolation quality)
  call AddNfInfoToPdfTable(table,coupling)

  ! "pre-evolve" so that subsequent evolutions are faster,
  ! specifying the initial scale for PDFs (here sqrt(2.0))
  !call PreEvolvePdfTable(table, sqrt(2.0_dp), dh, coupling)

end subroutine tabulation_example_startup

!----------------------------------------------------------------------
!!
subroutine tabulation_example_evolve_table
  use table_module
  use initial_condition
  implicit none
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)
  real(dp) :: Q0_pdf
  integer  :: nloop_ev
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


  ! initialise a PDF from the function below (must be contained,
  ! in a "used" module, or with an explicitly defined interface)
  call AllocPDF(grid, pdf0)

  ! get the PDF from the routine contained in the initial_condition
  ! module
  pdf0 = unpolarized_dummy_pdf(xValues(grid))
  Q0_pdf   = sqrt(2.0_dp)

  !! or get it from LHAPDF
  ! set up LHAPDF with cteq
  !call InitPDFsetByName("cteq61.LHgrid")
  !call InitPDF(0)
  ! allocate and set up the initial pdf from LHAPDF ...
  !Q0_pdf = 5.0_dp 
  !call InitPDF_LHAPDF(grid, pdf0, EvolvePDF, Q0_pdf)

  ! choose the number of loops (must not be higher than what's
  ! been initialised)
  nloop_ev = 3

  ! create the tabulation based on the evolution of pdf0 from scale Q0
  call EvolvePdfTable(table, Q0_pdf, pdf0, dh, coupling, nloop=nloop_ev)
  ! use the following if instead the table has been pre-evolved
  !call EvolvePdfTable(table,pdf0)
  write(6,'(a)') "Evolution done!"

  call Delete(pdf0)
end subroutine tabulation_example_evolve_table


!----------------------------------------------------------------------
!! print out the truncated sum rules
subroutine tabulation_print_sum_rules(Q)
  use table_module
  implicit none
  real(dp), intent(in) :: Q
  !-------------------------
  !! hold the pdf at scale Q
  real(dp), pointer :: pdf(:,:)

  ! allocate place to store the pdf at scale Q
  call AllocPDF(grid, pdf)
  ! evaluate the pdf at scale Q
  call EvalPDFTable_Q(table, Q, pdf)
  
  ! now output the results
  write(6,*)
  write(6,'(a,f8.3,a,es11.4)') "Evaluating sum rules at scale Q = ",Q, &
       &         " GeV, truncated at x = ", exp(-grid%ymax)
  write(6,'(a)') &
     &"[NB: O(0.5%) differences from expectations are due to truncation]"
  write(6,'(a,f7.4)') '  momentum sum: ', &
       &   TruncatedMoment(grid, sum(pdf(:,-6:6),dim=2), one )
  write(6,'(a,f7.4)') '  down valence: ', &
       &   TruncatedMoment(grid, pdf(:,1)-pdf(:,-1),     zero)
  write(6,'(a,f7.4)') '  up   valence: ', &
       &   TruncatedMoment(grid, pdf(:,2)-pdf(:,-2),     zero)

  ! cleaning up
  call Delete(pdf)
  
end subroutine tabulation_print_sum_rules

!----------------------------------------------------------------------
!! deallocate things 
subroutine tabulation_example_cleanup
  use table_module
  implicit none
  ! some cleaning up (not strictly speaking needed, but illustrates
  ! how it's done)
  call Delete(table)
  call Delete(dh)
  call Delete(coupling)
  call Delete(grid)

end subroutine tabulation_example_cleanup

