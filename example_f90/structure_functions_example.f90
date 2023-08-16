!! An example program using structure functions up to N3LO. The
!! program uses LHAPDF. 
!!
!!

program structure_functions_example
  use hoppet_v1
  use dummy_pdfs
  use streamlined_interface
  use structure_functions
  real(dp) :: Qmax, xmur, xmuf, Qmin, ymax
  integer  :: ipdf, order_max, sc_choice
  
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  ! use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  !implicit none

  ipdf = 91200
  Qmax = 13000.0_dp 
  order_max = 4
  xmur = one
  xmuf = one
  sc_choice = 1
  
  ! initialize PDF set
  call PDFSET('DEFAULT', dble(ipdf))
  call getQ2min(0,Qmin)
  Qmin = sqrt(Qmin)

  ! initialise hoppet
  call StartStrFct(Qmax, order_max, xR = xmur, xF = xmuf, sc_choice = sc_choice, &
       param_coefs = .true., Qmin_PDF = Qmin)
  call read_PDF()
  call InitStrFct(order_max, .true.)

  open(unit = 99, file = 'structure-functions.dat')

  ymax = log(1e5) !ymax=20
  call write_f1(99, 100.0_dp, ymax, 100)
  call write_f2(99, 100.0_dp, ymax, 100)
  call write_f3(99, 100.0_dp, ymax, 100)

contains 
  !----------------------------------------------------------------------
  ! fill the streamlined interface PDF table (possibly using hoppet's
  ! evolution)
  subroutine read_PDF(toyQ0, dglapQ0, xR_PDF)
    real(dp), external :: alphasPDF
    real(dp), optional :: toyQ0, dglapQ0, xR_PDF
    interface
       subroutine EvolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine EvolvePDF
    end interface
    !----------------
    real(dp) :: toy_Q0, Q0pdf, xmuR_PDF
    real(dp) :: res_lhapdf(-6:6), x, Q
    real(dp) :: res_hoppet(-6:6)
    real(dp) :: toy_pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)
    real(dp), parameter :: mz = 91.2_dp
    real(dp) :: pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)

    toy_Q0       = -one
    Q0pdf        = -one
    xmuR_PDF     = one  
    if(present(toyQ0)) toy_Q0=toyQ0
    if(present(dglapQ0)) Q0pdf=dglapQ0
    if(present(xR_PDF)) xmuR_PDF=xR_PDF

    if (toy_Q0 > zero) then
       write(6,*) "WARNING: Using toy PDF"
       toy_pdf_at_Q0 = unpolarized_dummy_pdf(xValues(grid))
       call InitRunningCoupling(toy_coupling, alfas=toy_alphas_Q0, &
            &                   nloop = 3, Q = toy_Q0, fixnf=nf_int)
       call EvolvePdfTable(tables(0), toy_Q0, toy_pdf_at_Q0, dh, toy_coupling, nloop=3)
    elseif (Q0pdf > zero) then

       write(6,*) "WARNING: Using internal HOPPET DGLAP evolution"
       call InitPDF_LHAPDF(grid, pdf_at_Q0, EvolvePDF, Q0pdf)
       call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , 4,&
            & -1000000045, quark_masses_sf(4:6), .true.)
       call EvolvePdfTable(tables(0), Q0pdf, pdf_at_Q0, dh, coupling, &
            &  muR_Q=xmuR_PDF, nloop=3)

    else
       ! InitRunningCoupling has to be called for the HOPPET coupling to be initialised 
       ! Default is to ask for 4 loop running and threshold corrections at quark masses.  
       call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , 4,&
            & -1000000045, quark_masses_sf(4:6), masses_are_MSbar = .true.)
       ! fixnf can be set to a positive number for
       ! fixed nf. -1000000045 gives variable nf
       ! and threshold corrections at quarkmasses.
       call hoppetAssign(EvolvePDF)
    end if

    ! quickly test that we have read in the PDFs correctly
    write(6,*) "Quick test that PDFs have been read in correctly"
    x = 0.08_dp
    Q = 17.0_dp
    call EvolvePDF(x, Q, res_lhapdf)
    call EvalPdfTable_xQ(tables(0), x, Q, res_hoppet)
    write(6,*) 'lhapdf: ', res_lhapdf
    write(6,*) 'hoppet: ', res_hoppet
  end subroutine read_PDF

end program structure_functions_example


