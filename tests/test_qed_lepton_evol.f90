!
! run this with
!
! ./test_qed_evol_lhapdf -Qlo 3.0 -Qhi 100  -pdf CT14qed_proton -nqcdloop-qed 0  -x 1e-3 >! a
!
!
!! TESTS
!! * Test against plain qcd evolution, with QED alpha turned to zero
!!   -> works OK, but watch out for the fact that the tau introduces
!!      an extra threshold, so that numerically there are differences
!!      at the evolution precision; reducing du by a factor of 10
!!      significantly reduces these differences. See basic test output with
!!   ./test_qed_evol_lhapdf -no-qed -x 0.5 | grep moment
!!   ./test_qed_evol_lhapdf -alpha-qed 0 -x 0.5 | grep moment
!!
!! * test momentum sum rule, with both nqcdloop_qed values
!!   -> this seems to work fine, e.g. test with
!!      ./test_qed_evol_lhapdf -ymax 25 -nqcdloop-qed 0 |& grep moment
!!      ./test_qed_evol_lhapdf -ymax 25 -nqcdloop-qed 1 |& grep moment
!!
!! * compare to photon evolution from another source?
!!   -> done in with-lhapdf/test_qed_evol_lhapdf.f90
!!      generally, answers are to within O(1%) of the MRST2004qed and
!!      CT14qed_proton results for evolution from 3 to 100 GeV. E.g.
!!    ./test_qed_evol_lhapdf -Qlo 3.0 -Qhi 100  -pdf CT14qed_proton -nqcdloop-qed 0  -x 1e-3 >! a
!!
program test_qed_evol
   use hoppet
  use qed_evolution; use qed_objects; use qed_coupling_module
  use io_utils
  implicit none
  character(len=200) :: pdfname
  type(grid_def)      :: grid
  type(qed_split_mat) :: qed_split
  type(dglap_holder)  :: dh
  real(dp)               :: quark_masses(4:6)
  type(running_coupling) :: coupling
  type(qed_coupling)     :: coupling_qed
  real(dp)               :: ymax, dy, x, alphas_mz, lha_qmin
  real(dp), parameter    :: mz = 91.1876_dp
  real(dp), pointer :: pdf_in(:,:), pdf_out(:,:)
  real(dp)         :: Qlo, Qhi, moments(ncompmin:ncompmaxLeptons)
  integer :: compact_output_n = -1
  integer :: nloop_qcd = 3, nqcdloop_qed, imem, iflv, iunit
  logical :: no_qed

  if (log_val_opt("-out")) then
     iunit = idev_open_opt("-out")
  else
     iunit = 6
  end if
  compact_output_n = int_val_opt("-compact-output", compact_output_n)
  
  write(iunit,'(a)') "# "//trim(command_line())

  ymax  = dble_val_opt("-ymax",20.0_dp)
  dy    = dble_val_opt("-dy",0.1_dp)
  call InitGridDefDefault(grid, dy, ymax)
  write(iunit,'(a,2f10.6)') "# ymax, dy ", ymax, dy

  ! allow option to just write out the x values of the grid
  ! and exit
  if (log_val_opt("-xvalues")) then
     write(iunit,"(a)") "# next lines contain number of x values, and actual x values"
     write(iunit,*) grid%ny+1
     write(iunit,*) xValues(grid)
     stop
  end if

  Qlo = m_electron * 1.00001_dp ! the initial scale
  Qhi = 40.0_dp
  Qlo = sqrt(dble_val_opt("-Q2lo", Qlo**2))
  Qhi = sqrt(dble_val_opt("-Q2hi", Qhi**2))
  Qlo = dble_val_opt("-Qlo",Qlo)
  Qhi = dble_val_opt("-Qhi",Qhi)
  write(iunit,'(a,es14.6)') "# Qlo = ", Qlo
  write(iunit,'(a,es14.6)') "# Qhi = ", Qhi


  quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
  alphas_mz = 0.118_dp
  write(iunit,*) '# quark masses (c,b,t) = ', quark_masses
  

  nloop_qcd = int_val_opt("-nloop-qcd",nloop_qcd)
  write(iunit,*) "# nloop_qcd = ", nloop_qcd

  nqcdloop_qed = int_val_opt('-nqcdloop-qed',0)
  write(iunit,*) "# nqcdloop_qed = ", nqcdloop_qed
  
  no_qed = log_val_opt("-no-qed",.false.)

  x = dble_val_opt("-x",1e-3_dp)

  
  call SetDefaultEvolutionDu(dy/3.0_dp)  ! generally a good choice
  call InitQEDSplitMat(grid, qed_split)
  call InitDglapHolder(grid,dh,factscheme=factscheme_MSbar,&
       &                      nloop=nloop_qcd,nflo=3,nfhi=6)


  call InitRunningCoupling(coupling,alfas=alphas_mz,Q=mz,nloop=nloop_qcd,&
       &                   quark_masses = quark_masses)
  if (log_val_opt("-alpha-qed")) then
     call InitQEDCoupling(coupling_qed, quark_masses(4), quark_masses(5), quark_masses(6), &
          &               dble_val_opt("-alpha-qed"))
  else
     call InitQEDCoupling(coupling_qed, quark_masses(4), quark_masses(5), quark_masses(6))
  end if

  
  write(iunit,*) "# 1/QED coupling at Me = ", one/Value(coupling_qed,m_electron)
  write(iunit,*) "# 1/QED coupling at MZ = ", one/Value(coupling_qed,MZ)
  write(iunit,*) "# QCD coupling at MZ = ", Value(coupling,MZ)
  write(iunit,*) "# Qlo, Qhi ", Qlo, Qhi
  
  call AllocPDFWithLeptons(grid, pdf_in)
  call AllocPDFWithLeptons(grid, pdf_out)
  
  pdf_in = zero
  ! then add in a lepton of some kind
  pdf_in(:,iflv_electron) = electron_dummy_pdf(xValues(grid))
  

  if (.not.CheckAllArgsUsed(0)) stop

  pdf_out = pdf_in
  if (no_qed) then
     call EvolvePDF(dh, pdf_out(:,ncompmin:ncompmax), coupling, Qlo, Qhi, nloop=nloop_qcd)
  else
     call QEDQCDEvolvePDF(dh, qed_split, pdf_out, coupling, coupling_qed,&
          &               Qlo, Qhi, nloop_qcd, nqcdloop_qed)
  end if
     

  call write_moments(pdf_in, ' in')
  call write_moments(pdf_out,'out')
  call write_pdf(pdf_in,'in')
  call write_pdf(pdf_out,'out')

  write(0,*) 'x = ', x, 'Q = ', Qhi
  write(0,'(a5,4a17)') 'iflv','hoppet'
  do iflv = -6, 11
     write(0,'(i5,4f17.7)') iflv, pdf_out(:,iflv).atx.(x.with.grid)
  end do
  
contains
  
  !----------------------------------------------------------------------
  subroutine write_moments(pdf,label)
    real(dp),         intent(in) :: pdf(0:,ncompmin:)
    character(len=*), intent(in) :: label
    !-----------------
    real(dp) :: moments(ncompmin:ubound(pdf,dim=2))
    
    write(iunit,"(a)", advance="no") "# total momentum "//trim(label)//" (& components)"
    moments = TruncatedMoment(grid, pdf, one)
    write(iunit,"(40f10.7)") sum(moments), moments

    write(iunit,"(a)", advance="no") "# total number "//trim(label)//" (& components)"
    moments(1:6) = TruncatedMoment(grid, pdf(:,1:6)-pdf(:,-1:-6:-1), zero)
    write(iunit,"(40f11.7)") sum(moments(1:6)), moments(1:6)


  end subroutine write_moments
  
  
  !----------------------------------------------------------------------
  subroutine write_pdf(pdf,label)
    real(dp), intent(in) :: pdf(0:,ncompmin:)
    character(len=*), intent(in) :: label
    !----------------------
    real(dp) :: y, yVals(0:grid%ny), last_y
    real(dp) :: zeta, zeta_max, a_stretch=7.0_dp, pdf_val
    integer  :: i,j

    write(iunit,"(a)") "# pdf "//trim(label)
    !do i = 0, 100
    !   y = 0.1_dp * i
    !   write(iunit,*) y, pdf.aty.(y.with.grid)
    !end do
    if (compact_output_n > 0) then
       ! value of compact_output_n sets number of points being written
       zeta_max = zeta_of_y(grid%ymax, a_stretch)
       do i = 0, compact_output_n
          zeta = (i*zeta_max)/compact_output_n
          y = y_of_zeta(zeta, a_stretch)
          write(iunit,'(f10.6)',advance="no") y
          do j = lbound(pdf,dim=2), ubound(pdf,dim=2)
             pdf_val = pdf(:,j).aty.(y.with.grid)
             if (abs(pdf_val) < 1e-100_dp) then
                write(iunit,'(a)',advance='no') ' 0 '
             else
                write(iunit,'(es14.6)',advance='no') pdf_val
             end if
          end do
          write(iunit,*)
       end do
    else
       last_y = -one
       yVals = yValues(grid)
       do i = 0, grid%ny
          y = yVals(i)
          if (y < 1.000000001_dp * last_y) cycle
          write(iunit,*) y, pdf(i,:)
          last_y = y
       end do
    end if
    
    write(iunit,*)
    write(iunit,*)
    
  end subroutine write_pdf
  

  function electron_dummy_pdf(xvals) result(pdf)
    implicit none
    real(dp), intent(in) :: xvals(0:)
    real(dp)             :: pdf(0:ubound(xvals,1))
    pdf = xvals**2 * (one-xvals)**4 * 100.0_dp
  end function electron_dummy_pdf
  
  
  !-----------------------------------------------------------------
  !! return zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  function zeta_of_y(y, a) result(zeta)
    real(dp), intent(in) :: y, a
    real(dp)             :: zeta
    zeta = y + a*(one - exp(-y))
  end function zeta_of_y

  
  !-----------------------------------------------------------------
  !! return inversion of zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  function y_of_zeta(zeta, a) result(y)
    real(dp), intent(in) :: zeta, a
    real(dp)             :: y, x, diff_from_zero, deriv
    integer             :: iter
    real(dp), parameter :: eps = 1e-12_dp
    integer,  parameter :: maxiter = 100

    ! starting condition (and soln if a = 0)
    y = zeta 
    if (a /= zero) then
       do iter = 0, maxiter
          x = exp(-y);
          diff_from_zero = zeta - y - a*(one-x);
          ! we have found good solution
          if (abs(diff_from_zero) < eps) exit
          deriv = -one  - a*x;
          y = y - diff_from_zero / deriv;
       end do
    end if
    
    if (iter > maxiter) write(0,*) "y_of_zeta reached maxiter"

  end function y_of_zeta

end program test_qed_evol
