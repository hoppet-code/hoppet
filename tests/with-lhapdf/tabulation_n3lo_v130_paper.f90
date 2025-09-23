!! A program that can be used to produce the results presented for
!! n3lo evolution in the v130 manual.
program tabulation_n3lo
  use io_utils, int_value => value
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
     use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use streamlined_interface, HoppetInitPDF => initPDF
  use hoppet_git_state
  implicit none
  real(dp) :: dy, ymax, dlnlnQ
  integer  :: order, nloop
  real(dp) :: quark_masses(4:6)
  !! hold the initial pdf
  real(dp), pointer :: pdf0(:,:)
  real(dp) :: Q0
  !! hold results at some x, Q
  real(dp) :: Q, pdf_at_xQ(-6:6)
  real(dp), parameter :: heralhc_xvals(11) = &
       & (/1e-7_dp, 1e-6_dp, 1e-5_dp,1e-4_dp,1e-3_dp,1e-2_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp/)
  integer  :: i, ix
  real(dp) :: x
  integer :: idev, n3lo_split_approx, n3lo_vfns_on, scheme
  character(200)   pdfname
  real(dp) :: mc, mb, mt, asQ0, QMin, Qmax, xmu, lhapdf_qmin, Qevolve
  double precision alphasPDF, hoppetAlphas
  external         evolvePDF, alphasPDF
  
  ! Get output and pdf name here
  idev = idev_open_opt("-o","")
  pdfname = string_val_opt("-pdf","toyHERA")
  write(idev,'(a,a)') "# ", trim(command_line())
  call hoppet_print_git_state(idev,prefix="# ")
    
  ! set up parameters for grid
  Qmin  = dble_val_opt("-Qmin",1.0_dp)  ! smallest Q value in tabulation
  Qmax  = 1d4       ! largest Q value in tabulation
  order  = int_val_opt("-order",-6) ! numerical interpolation order (-6 is a good choice) 
  ymax   = dble_val_opt("-ymax",17.0_dp) ! max value of ln 1/x
  dy     = dble_val_opt("-dy",0.05_dp) ! the internal grid spacing (smaller->higher accuarcy)
  dlnlnQ = dble_val_opt("-dlnlnQ",dy/8.0_dp) ! tabulation spacing in dlnlnQ (dy/4 recommended) 
  scheme = 1 ! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar
  Qevolve = dble_val_opt("-Qevolve",100.0_dp)
  ! Determine evolution order
  nloop = int_val_opt("-nloop",4) ! 1; LO, 2; NLO, 3; NNLO; 4; N3LO
  n3lo_split_approx = int_val_opt("-n3lo-split-approx",2) ! 0: n3lo_splitting_approximation_up_to_2310_05744
                                                          ! 1: n3lo_splitting_approximation_up_to_2404_09701
                                                          ! 2: n3lo_splitting_approximation_up_to_2410_08089
                                                          ! see dglap_choices.f90 for details
  n3lo_vfns_on = int_val_opt("-n3lo-vfns-on",1) ! 1: on, 0: off
  call dglap_Set_n3lo_splitting_approximation(n3lo_split_approx)
  call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact) ! Exact nnlo thresholds
  call dglap_Set_n3lo_nfthreshold(n3lo_vfns_on) ! To turn on/off the n3lo mass thresholds
  xmu = 1.0_dp ! Renormlization scale
  call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,scheme)
  write(6,'(a)') "Hoppet grid setup done!"

  call AllocPDF(grid, pdf0)
  if(trim(pdfname).eq."toyHERA") then
     mc = 1.414213563_dp
     mb = 4.5_dp
     mt = 175.0_dp
     ! initialise a PDF from the function below (must be contained,
     ! in a "used" module, or with an explicitly defined interface)
     pdf0 = unpolarized_dummy_pdf(xValues(grid))
     Q0 = dble_val_opt("-Q0",sqrt(2.0_dp))  ! the initial scale
     asQ0 = dble_val_opt("-asQ0",0.35_dp) ! alphas at Q0
  else
     call InitPDFsetByName(pdfname)
     call InitPDF(0)

     ! Retrieve starting scale and coupling from LHAPDF, and also the
     ! quark masses
     call getQ2min(0,lhapdf_qmin)
     lhapdf_qmin = sqrt(lhapdf_qmin)
     Q0 = dble_val_opt("-Q0",lhapdf_qmin)
     asQ0 = alphasPDF(Q0)
     call getthreshold(4,mc)
     call getthreshold(5,mb)
     call getthreshold(6,mt)
     call InitPDF_LHAPDF(grid, pdf0, evolvePDF, Q0)
  endif

  if (.not.CheckAllArgsUsed(0)) then
     call exit()
  endif

    
  call hoppetSetPoleMassVFN(mc,mb,mt)
  quark_masses(4:6) = (/mc,mb,mt/)
  call InitRunningCoupling(coupling,alfas=asQ0,Q=Q0,nloop=nloop,&
       &                   quark_masses = quark_masses)
  call AddNfInfoToPdfTable(tables,coupling)
  coupling_initialised = .true.
  write(6,'(a)') "Coupling initialised!"
  write(6,'(a,f7.4,a,f7.4,a,f7.3,f7.3,f9.3)') "Q0 = ", Q0, ", αS(Q0) = ", asQ0, ", (mc,mb,mt) = ", quark_masses
  write(idev,'(a,f8.4,a,f7.4,a,f7.3,f7.3,f9.3)') "# Q0 = ", Q0, ", αS(Q0) = ", asQ0, ", (mc,mb,mt) = ", quark_masses
  write(6,'(a,f7.4)') "Hoppet αS(Q0) = ", hoppetAlphas(Q0)
! create the tabulation based on the evolution of pdf0 from scale Q0
  call EvolvePdfTable(tables(0), Q0, pdf0, dh, coupling, nloop=nloop)
  
    
  write(6,'(a)') "Evolution done!"

  ! get the value of the tabulation at some point
  Q = Q0
  write(6,'(a)')
  write(6,'(a,f10.5,a)') "                                   Evaluating PDFs at Q = ",Q," GeV"
  write(6,'(a5,3a12,a14,3a11,a12)') "x",&
       & "u-ubar","d-dbar","dbr-ubr","2(ubr+dbr)","s+sbar","c+cbar","b+bbar","gluon"
  do ix = 1, size(heralhc_xvals)
     call EvalPdfTable_xQ(tables(0),heralhc_xvals(ix),Q,pdf_at_xQ)
     write(6,'(es7.1,8es12.4)') heralhc_xvals(ix), &
          &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
          &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
          &  (pdf_at_xQ(-1)-pdf_at_xQ(-2)), &
          &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
          &  (pdf_at_xQ(-3)+pdf_at_xQ(3)), &
          &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
          &  (pdf_at_xQ(-5)+pdf_at_xQ(5)), &
          &  pdf_at_xQ(0)
  end do

  if(trim(pdfname).ne."toyHERA") then
     ! get the value of the tabulation at some point
     Q = Q0
     write(6,'(a)')
     write(6,'(a,f10.5,a)') "                                   LHAPDF interpolation at Q = ",Q," GeV"
     write(6,'(a5,3a12,a14,3a11,a12)') "x",&
          & "u-ubar","d-dbar","dbr-ubr","2(ubr+dbr)","s+sbar","c+cbar","b+bbar","gluon"
     do ix = 1, size(heralhc_xvals)
        call evolvePDF(heralhc_xvals(ix),Q,pdf_at_xQ)
        write(6,'(es7.1,8es12.4)') heralhc_xvals(ix), &
             &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
             &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
             &  (pdf_at_xQ(-1)-pdf_at_xQ(-2)), &
             &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
             &  (pdf_at_xQ(-3)+pdf_at_xQ(3)), &
             &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
             &  (pdf_at_xQ(-5)+pdf_at_xQ(5)), &
             &  pdf_at_xQ(0)
     end do
  end if

  ! get the value of the tabulation at some point
  Q = Qevolve
  write(6,'(a)')
  write(6,'(a,f10.5,a)') "                                   Evaluating PDFs at Q = ",Q," GeV"
  write(6,'(a5,3a12,a14,3a11,a12)') "x",&
       & "u-ubar","d-dbar","dbr-ubr","2(ubr+dbr)","s+sbar","c+cbar","b+bbar","gluon"
  do ix = 1, size(heralhc_xvals)
     call EvalPdfTable_xQ(tables(0),heralhc_xvals(ix),Q,pdf_at_xQ)
     write(6,'(es7.1,8es12.4)') heralhc_xvals(ix), &
          &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
          &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
          &  (pdf_at_xQ(-1)-pdf_at_xQ(-2)), &
          &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
          &  (pdf_at_xQ(-3)+pdf_at_xQ(3)), &
          &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
          &  (pdf_at_xQ(-5)+pdf_at_xQ(5)), &
          &  pdf_at_xQ(0)
  end do

  if(trim(pdfname).ne."toyHERA") then
     ! get the value of the tabulation at some point
     Q = Qevolve
     write(6,'(a)')
     write(6,'(a,f10.5,a)') "                                   LHAPDF interpolation at Q = ",Q," GeV"
     write(6,'(a5,3a12,a14,3a11,a12)') "x",&
          & "u-ubar","d-dbar","dbr-ubr","2(ubr+dbr)","s+sbar","c+cbar","b+bbar","gluon"
     do ix = 1, size(heralhc_xvals)
        call evolvePDF(heralhc_xvals(ix),Q,pdf_at_xQ)
        write(6,'(es7.1,8es12.4)') heralhc_xvals(ix), &
             &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
             &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
             &  (pdf_at_xQ(-1)-pdf_at_xQ(-2)), &
             &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
             &  (pdf_at_xQ(-3)+pdf_at_xQ(3)), &
             &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
             &  (pdf_at_xQ(-5)+pdf_at_xQ(5)), &
             &  pdf_at_xQ(0)
     end do
  endif

  ! Write to a file
  Q = Qevolve
  write(idev,'(a, f10.5)') "# hoppet at Q  = ", Q
  write(idev,'(a8,3a15,a17,5a14,a15)') "x",&
       & "u-ubar","d-dbar","dbr-ubr","2(ubr+dbr)","s+sbar","s-sbar","c+cbar","c-cbar","b+bbar","gluon"
  do i = 1, 200
     x     = 0.00001_dp ** (i/200.0_dp)
     call EvalPdfTable_xQ(tables(0),x,Q,pdf_at_xQ)
     write(idev,'(es10.4,10es15.7)') x, &
          &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
          &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
          &  (pdf_at_xQ(-1)-pdf_at_xQ(-2)), &
          &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
          &  (pdf_at_xQ(-3)+pdf_at_xQ(3)), &
          &  (pdf_at_xQ(-3)-pdf_at_xQ(3)), &
          &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
          &  (pdf_at_xQ(-4)-pdf_at_xQ(4)), &
          &  (pdf_at_xQ(-5)+pdf_at_xQ(5)), &
          &  pdf_at_xQ(0)
   end do
   write(idev,*) ''
   write(idev,*) ''

  write(idev,'(a, f10.5)') "# LHAPDF at Q  = ", Q
  write(idev,'(a8,3a15,a17,5a14,a15)') "x",&
       & "u-ubar","d-dbar","dbr-ubr","2(ubr+dbr)","s+sbar","s-sbar","c+cbar","c-cbar","b+bbar","gluon"
  do i = 1, 200
     x     = 0.00001_dp ** (i/200.0_dp)
     call evolvePDF(x,Q,pdf_at_xQ)
     write(idev,'(es10.4,10es15.7)') x, &
          &  pdf_at_xQ(2)-pdf_at_xQ(-2), &
          &  pdf_at_xQ(1)-pdf_at_xQ(-1), &
          &  (pdf_at_xQ(-1)-pdf_at_xQ(-2)), &
          &  2*(pdf_at_xQ(-1)+pdf_at_xQ(-2)), &
          &  (pdf_at_xQ(-3)+pdf_at_xQ(3)), &
          &  (pdf_at_xQ(-3)-pdf_at_xQ(3)), &
          &  (pdf_at_xQ(-4)+pdf_at_xQ(4)), &
          &  (pdf_at_xQ(-4)-pdf_at_xQ(4)), &
          &  (pdf_at_xQ(-5)+pdf_at_xQ(5)), &
          &  pdf_at_xQ(0)
   end do
   write(idev,*) ''
   write(idev,*) ''
   

contains 


  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
  function unpolarized_dummy_pdf(xvals) result(pdf)
    use types; use consts_dp; use pdf_representation
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

end program tabulation_n3lo


