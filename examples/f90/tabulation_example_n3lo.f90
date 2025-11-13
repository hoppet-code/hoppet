!! An example program using a tabulation. It outputs a subset of
!! table 15 of hep-ph/0511119 and this output should be identical
!! to the contents of the file tabulation_example.default_output
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
!! Pipe the output through
!!
!!   pcresed 's/E(.).(.)/\$^\{$1$2\}\$/g' | pcresed 's/  ([0-9])/ &  $1/g' | pcresed 's/ -/ & -/g'
!!
!! to get the number in LaTeX table format
program tabulation_example_n3lo
  use hoppet
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  !use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  implicit none
  real(dp) :: dy, ymax, dlnlnQ
  integer  :: order, nloop
  !! holds information about the grid
  type(grid_def) :: grid
  !! holds the splitting functions
  type(dglap_holder) :: dh
  !! hold the PDF tabulation
  type(pdf_table)       :: table
  !! hold the coupling
  real(dp)               :: quark_masses(4:6)
  type(running_coupling) :: coupling
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)
  real(dp) :: Q0
  !! hold results at some x, Q
  real(dp) :: Q, pdf_at_xQ(-6:6)
  real(dp), parameter :: heralhc_xvals(11) = &
       & (/1e-7_dp, 1e-6_dp, 1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: i, ix
  real(dp) :: x
  logical :: vfn = .true. ! Change to false to run nf = 4 FFN for comparison with 2406.16188

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

  ! set up parameters for grid
  order = -6
  ymax  = 17.0_dp
  dy    = 0.05_dp      ! a large value of 0.20 gives results that are good to within relative 10^{-4}, 0.1 agrees with 0.05 to within relative 10^{-5}
  dlnlnQ = dy / 4.0_dp ! a good default as long as ymax is not too large

  ! set up the grid itself (this call sets up a nested grid composed of 4 subgrids)
  call InitGridDefDefault(grid, dy, ymax, order=order)

  ! initialise the splitting-function holder
  nloop = 4 ! 1; LO, 2; NLO, 3; NNLO; 4; N3LO
  if(vfn) then 
     !n3lo_splitting_approximation = n3lo_splitting_approximation_up_to_2310_05744 ! Controls the approximation that we use the in n3lo splitting functions
     call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact)
     call dglap_Set_n3lo_nfthreshold(n3lo_nfthreshold_on) ! To turn on/off the n3lo mass thresholds
     call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
          &                      nloop=nloop,nflo=3,nfhi=5)
     write(6,'(a)') "Splitting functions initialised!"
     write(6,'(a,i1)') "Using nloop = ",nloop
     write(6,'(a,f6.4,a,f8.6)') "Using dy = ",dy," dlnlnQ = ",dlnlnQ

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
     Q0 = sqrt(2.0_dp)  ! the initial scale
     
     ! allocate and initialise the running coupling with a given
     ! set of quark masses (NB: charm mass just above Q0).
     quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 1e100_dp/)
     call InitRunningCoupling(coupling,alfas=0.35_dp,Q=Q0,nloop=nloop,&
          &                   quark_masses = quark_masses)
     
  else
     ! To perform the nf = 4 FFN run that compares against 2406.16188
     ! tables 1 & 2 uncomment the below and comment the above
     n3lo_splitting_approximation = n3lo_splitting_approximation_up_to_2310_05744
     call dglap_Set_n3lo_nfthreshold(n3lo_nfthreshold_off) ! To turn on/off the n3lo mass thresholds
     call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
          &                      nloop=nloop,nflo=4,nfhi=4)
     write(6,'(a)') "Doing FFN nf = 4 tabulation for comparison with 2406.16188"
     write(6,'(a)') "Splitting functions initialised, using 2310_05744 approximation!"
     write(6,'(a,i1)') "Using nloop = ",nloop
     write(6,'(a,f6.4,a,f8.6)') "Using dy = ",dy," dlnlnQ = ",dlnlnQ
     
     ! initialise a PDF from the function below (must be contained,
     ! in a "used" module, or with an explicitly defined interface)
     call AllocPDF(grid, pdf0)
     pdf0 = unpolarized_dummy_pdf(xValues(grid))
     Q0 = sqrt(2.0_dp)  ! the initial scale
     
     ! allocate and initialise the running coupling with a given
     ! set of quark masses (NB: charm mass just above Q0).
     quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 1e100_dp/)
     call InitRunningCoupling(coupling,alfas=0.35_dp,Q=Q0,nloop=nloop,&
          &                   quark_masses = quark_masses, fixnf = 4)
  endif
  ! create the tables that will contain our copy of the user's pdf
  ! as well as the convolutions with the pdf.
  call AllocPdfTable(grid, table, Qmin=1.0_dp, Qmax=10000.0_dp, & 
       & dlnlnQ = dlnlnQ, freeze_at_Qmin=.true.)
  ! add information about the nf transitions to the table (improves
  ! interpolation quality)
  call AddNfInfoToPdfTable(table,coupling)

  ! create the tabulation based on the evolution of pdf0 from scale Q0
  call EvolvePdfTable(table, Q0, pdf0, dh, coupling, nloop=nloop)
  ! alternatively "pre-evolve" so that subsequent evolutions are faster
  !call PreEvolvePdfTable(table, Q0, dh, coupling)
  !call EvolvePdfTable(table,pdf0)
  write(6,'(a)') "Evolution done!"

  ! get the value of the tabulation at some point
  Q = 100.0_dp
  write(6,'(a)')
  write(6,'(a,f8.3,a)') "                                   Evaluating PDFs at Q = ",Q," GeV"
  write(6,'(a5,3a12,a14,4a11,a12)') "x",&
       & "u-ubar","d-dbar","dbr-ubr","2(ubr+dbr)","s-sbar","s+sbar","c+cbar","b+bbar","gluon"
  do ix = 1, size(heralhc_xvals)
     call EvalPdfTable_xQ(table,heralhc_xvals(ix),Q,pdf_at_xQ)
     write(6,'(es7.1,9es12.4)') heralhc_xvals(ix), &
          &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
          &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
          &  (pdf_at_xQ(-1)-pdf_at_xQ(-2)), &
          &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
          &  (pdf_at_xQ(3)-pdf_at_xQ(-3)), &
          &  (pdf_at_xQ(-3)+pdf_at_xQ(3)), &
          &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
          &  (pdf_at_xQ(-5)+pdf_at_xQ(5)), &
          &  pdf_at_xQ(0)
  end do
  ! Write to a file
  
  write(12,'(a5,3a12,a14,4a11,a12)') "x",&
       & "u-ubar","d-dbar","dbr-ubr","2(ubr+dbr)","s-sbar","s+sbar","c+cbar","b+bbar","gluon"
  do i = 1, 200
     x     = 0.0001_dp ** (i/200.0_dp)
     call EvalPdfTable_xQ(table,x,Q,pdf_at_xQ)
     write(12,'(es10.4,9es14.7)') x, &
          &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
          &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
          &  (pdf_at_xQ(-1)-pdf_at_xQ(-2)), &
          &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
          &  (pdf_at_xQ(3)-pdf_at_xQ(-3)), &
          &  (pdf_at_xQ(-3)+pdf_at_xQ(3)), &
          &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
          &  (pdf_at_xQ(-5)+pdf_at_xQ(5)), &
          &  pdf_at_xQ(0)
   end do
   write(12,*) ''
   write(12,*) ''
  
  ! some cleaning up (not strictly speaking needed, but illustrates
  ! how it's done)
  call Delete(table)
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

    ! labels iflv_g, etc., come from the hoppet module, inherited
    ! from the main program
    pdf(:, iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf

end program tabulation_example_n3lo


