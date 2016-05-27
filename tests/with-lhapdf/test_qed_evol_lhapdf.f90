program test_qed_evol
  use hoppet_v1
  use qed_evolution; use qed_objects; use qed_coupling_module
  use sub_defs_io
  implicit none
  character(len=200) :: pdfname
  type(grid_def)      :: grid
  type(qed_split_mat) :: qed_split
  type(dglap_holder)  :: dh
  real(dp)               :: quark_masses(4:6)
  type(running_coupling) :: coupling
  type(qed_coupling)     :: coupling_qed
  real(dp)               :: ymax, dy
  real(dp), pointer :: pdf_in(:,:), pdf_out(:,:), pdf_lhapdf_out(:,:)
  real(dp)         :: Q0, Qfinal, moments(ncompmin:ncompmaxLeptons)
  real(dp) :: alphaspdf ! from LHAPDF
  integer :: nloop_qcd, nqcdloop_qed, imem
  logical :: no_qed

  pdfname = string_val_opt("-pdf")
  imem    = int_val_opt("-imem",0)
  call InitPDFsetByName(trim(pdfname))
  call InitPDFm(1,imem)
  
  ymax  = 15.0_dp
  dy    = 0.1_dp
  nloop_qcd = int_val_opt("-nloop-qcd",3)
  nqcdloop_qed = int_val_opt('-nqcdloop-qed',0)
  Q0 = dble_val_opt("-Qlo",sqrt(2.0_dp))  ! the initial scale
  Qfinal = dble_val_opt("-Qhi",40.0_dp)
  no_qed = log_val_opt("-no-qed",.false.)

  if (.not.CheckAllArgsUsed(0)) stop
  
  call InitGridDefDefault(grid, dy, ymax)
  call SetDefaultEvolutionDu(dy/3.0_dp)  ! generally a good choice
  call InitQEDSplitMat(grid, qed_split)
  call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop_qcd,nflo=3,nfhi=6)

  !quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
  quark_masses(4:6) = (/lhapdf_qmass(4), lhapdf_qmass(6), lhapdf_qmass(6)/)

  call InitRunningCoupling(coupling,alfas=alphaspdf(91.2_dp),Q=91.2_dp,nloop=nloop_qcd,&
       &                   quark_masses = quark_masses)
  call InitQEDCoupling(coupling_qed, quark_masses(4), quark_masses(5), quark_masses(6))
  write(6,*) "# QED coupling at 1 GeV = ", Value(coupling_qed,one)
  write(6,*) "# QCD coupling at 91.2 GeV = ", Value(coupling,91.2_dp)
  
  call AllocPDFWithLeptons(grid, pdf_in)
  call AllocPDFWithLeptons(grid, pdf_out)
  call AllocPDFWithLeptons(grid, pdf_lhapdf_out)
  
  call fill_from_lhapdf(pdf_in, Q0)
  call fill_from_lhapdf(pdf_lhapdf_out, Qfinal)

  pdf_out = pdf_in
  if (no_qed) then
     call EvolvePDF(dh, pdf_out(:,ncompmin:ncompmax), coupling, Q0, Qfinal, nloop=nloop_qcd)
  else
     call QEDQCDEvolvePDF(dh, qed_split, pdf_out, coupling, coupling_qed,&
          &               Q0, Qfinal, nloop_qcd, nqcdloop_qed)
  end if
  

  call write_moments(pdf_in, ' in')
  call write_moments(pdf_lhapdf_out,'lha')
  call write_pdf(pdf_in,'in')
  call write_pdf(pdf_out,'out')
  call write_pdf(pdf_lhapdf_out,'lhapdf out')
  
  !! TESTS
  !! * Test against plain qcd evolution, with QED alpha turned to zero
  !!   -> works OK, but watch out for the fact that the tau introduces
  !!      an extra threshold, so that numerically there are differences
  !!      at the evolution precision; reducing du by a factor of 10
  !!      significantly reduces these differences
  !!
  !! * test momentum sum rule, with both nqcdloop_qed values
  !!   -> this seems to work fine
  !!
  !! * compare to photon evolution from another source?
contains

  function lhapdf_qmass(iflv) result(res)
    integer, intent(in) :: iflv
    real(dp)            :: res
    call getqmass(iflv,res)
  end function lhapdf_qmass
  
  
  !----------------------------------------------------------------------
  subroutine fill_from_lhapdf(this_pdf, Q)
    real(dp), intent(out) :: this_pdf(0:,-6:)
    real(dp), intent(in)  :: Q
    !---------
    real(dp) :: xv(0:grid%ny)
    integer  :: i

    xv = xValues(grid)
    this_pdf = zero
    do i = 0, grid%ny
       call EvolvePDFPhoton(xv(i), Q, this_pdf(i,-6:6), this_pdf(i,8))
       !write(6,*) xv(i), this_pdf(i,8)
    end do
    
  end subroutine fill_from_lhapdf

  !----------------------------------------------------------------------
  subroutine write_moments(pdf,label)
    real(dp),         intent(in) :: pdf(0:,ncompmin:)
    character(len=*), intent(in) :: label
    !-----------------
    real(dp) :: moments(ncompmin:ubound(pdf,dim=2))
    
    write(6,"(a)", advance="no"), "# total momentum "//trim(label)//" (& components)"
    moments = TruncatedMoment(grid, pdf, one)
    write(6,"(40f10.7)") sum(moments), moments

    write(6,"(a)", advance="no"), "# total number "//trim(label)//" (& components)"
    moments(1:6) = TruncatedMoment(grid, pdf(:,1:6)-pdf(:,-1:-6:-1), zero)
    write(6,"(40f11.7)") sum(moments(1:6)), moments(1:6)


  end subroutine write_moments
  
  
  !----------------------------------------------------------------------
  subroutine write_pdf(pdf,label)
    real(dp), intent(in) :: pdf(0:,ncompmin:)
    character(len=*), intent(in) :: label
    !----------------------
    real(dp) :: y, yVals(0:grid%ny), last_y
    integer  :: i

    write(6,"(a)"), "# pdf "//trim(label)
    !do i = 0, 100
    !   y = 0.1_dp * i
    !   write(6,*) y, pdf.aty.(y.with.grid)
    !end do
    last_y = -one
    yVals = yValues(grid)
    do i = 0, grid%ny
       y = yVals(i)
       if (y < 1.000000001_dp * last_y) cycle
       write(6,*) y, pdf(i,:)
       last_y = y
    end do
    write(6,*)
    write(6,*)
    
  end subroutine write_pdf
  
end program test_qed_evol
