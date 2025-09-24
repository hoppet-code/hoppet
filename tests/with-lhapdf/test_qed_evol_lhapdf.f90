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
  real(dp)               :: ymax, dy, x, alphas_mz, lha_qmin, scale_qg
  real(dp), parameter    :: mz = 91.1876_dp
  real(dp), pointer :: pdf_in(:,:), pdf_out(:,:), pdf_lhapdf_out(:,:)
  real(dp), pointer :: photon_lumi(:,:)
  real(dp)         :: Qlo, Qhi, moments(ncompmin:ncompmaxLeptons)
  real(dp) :: alphaspdf ! from LHAPDF
  integer :: compact_output_n = -1
  integer :: nloop_qcd = 3, nqcdloop_qed, imem, iflv, iunit
  logical :: no_qed, use_lhapdf

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

  Qlo = sqrt(two) ! the initial scale
  Qhi = 40.0_dp
  Qlo = sqrt(dble_val_opt("-Q2lo", Qlo**2))
  Qhi = sqrt(dble_val_opt("-Q2hi", Qhi**2))
  Qlo = dble_val_opt("-Qlo",Qlo)
  Qhi = dble_val_opt("-Qhi",Qhi)
  write(iunit,'(a,es14.6)') "# Qlo = ", Qlo
  write(iunit,'(a,es14.6)') "# Qhi = ", Qhi


  scale_qg = dble_val_opt("-scale-qg",one)
  
  use_lhapdf = log_val_opt("-pdf")
  if (use_lhapdf) then
     pdfname = string_val_opt("-pdf")
     imem    = int_val_opt("-imem",0)
     call InitPDFsetByName(trim(pdfname))
     call InitPDFm(1,imem)
     write(iunit,*) "# using LHAPDF ", trim(pdfname), ' imem = ', imem
     quark_masses(4:6) = (/lhapdf_qmass(4), lhapdf_qmass(5), lhapdf_qmass(6)/)
     alphas_mz = alphaspdf(mz)
     call getorderas(nloop_qcd); nloop_qcd = nloop_qcd + 1 ! LHAPDF convention is that LO=0
     ! the following lines are problematic for some reason
     ! if (Qlo < lhapdf_qmin()) &
    !      &write(0,*) "WARNING: Qlo = ", Qlo, " < lhapdf Qmin = ", lhapdf_qmin() 
  else
     write(iunit,*) "# using toy pdf"
     quark_masses(4:6) = (/1.414213563_dp, 4.5_dp, 175.0_dp/)
     alphas_mz = 0.118_dp
  end if
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

  
  write(iunit,*) "# 1/QED coupling at MZ = ", one/Value(coupling_qed,MZ)
  write(iunit,*) "# QCD coupling at MZ = ", Value(coupling,MZ)
  write(iunit,*) "# Qlo, Qhi ", Qlo, Qhi
  
  call AllocPDFWithLeptons(grid, pdf_in)
  call AllocPDFWithLeptons(grid, pdf_out)
  call AllocPDFWithLeptons(grid, pdf_lhapdf_out)
  call AllocGridQuant(grid, photon_lumi, 1, 1)
  
  if (use_lhapdf) then
     call fill_from_lhapdf(pdf_in, Qlo)
     call fill_from_lhapdf(pdf_lhapdf_out, Qhi)
  else
     pdf_in = zero
     ! fill in the QCD part
     pdf_in(:,:ncompmax) = unpolarized_dummy_pdf(xValues(grid))
     ! then add in a photon of some kind
     pdf_in(:,iflv_photon) = alpha_qed_scale_0 * (one - xValues(grid))**4
  end if
  if (scale_qg /= one) then
     pdf_in(:,:6) = pdf_in(:,:6) * scale_qg
  end if
  
  ! allow user to override the photon from a suitably formatted file
  if (log_val_opt("-read-photon")) pdf_in(:,iflv_photon) = read_photon(grid)

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
  if (use_lhapdf) call write_moments(pdf_lhapdf_out,'lha')
  call write_pdf(pdf_in,'in')
  call write_pdf(pdf_out,'out')
  if (use_lhapdf) call write_pdf(pdf_lhapdf_out,'lhapdf out')

  photon_lumi(:,1) = PartonLuminosity(grid, pdf_in(:,8),pdf_in(:,8))
  call write_pdf(photon_lumi,'photon-lumi-in')
  photon_lumi(:,1) = PartonLuminosity(grid, pdf_out(:,8),pdf_out(:,8))
  call write_pdf(photon_lumi,'photon-lumi-out')
  if (use_lhapdf)  then
     photon_lumi(:,1) = PartonLuminosity(grid, pdf_lhapdf_out(:,8),pdf_lhapdf_out(:,8))
     call write_pdf(photon_lumi,'photon-lumi-lhapdf-out')
  end if
  
  write(0,*) 'x = ', x, 'Q = ', Qhi
  write(0,'(a5,4a17)') 'iflv','lhapdf','hoppet'
  do iflv = -6, 8
     write(0,'(i5,4f17.7)') iflv, pdf_lhapdf_out(:,iflv).atx.(x.with.grid),&
          &                       pdf_out       (:,iflv).atx.(x.with.grid)
  end do
  
contains

  function lhapdf_qmin() result(qmin)
    real(dp) :: qmin, q2min
    call getq2min(q2min)
    write(0,*) q2min
    qmin = sqrt(q2min)
  end function lhapdf_qmin
  
  
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
       !write(iunit,*) xv(i), this_pdf(i,8)
    end do
  end subroutine fill_from_lhapdf

  !----------------------------------------------------------------------
  ! this assumes a 6-column format where column 1 is x and column 6 is
  ! xf_{gamma/p}(x)
  function read_photon(grid) result(photon)
    type(grid_def), intent(in) :: grid
    !---------------------------------
    integer             :: iy, iread
    character(len=1000) :: line
    real(dp)            :: photon(0:grid%ny), xVals(0:grid%ny), numbers(6)
    real(dp), parameter :: tolerance = 1e-7_dp
    
    iread = idev_open_opt("-read-photon",status='OLD')
    xVals = xValues(grid)
    write(0,*) 'reading photon from file',trim(string_val_opt("-read-photon"))
    do iy = 0, grid%ny
       ! skip comment lines
       do while (.true.)
          read(iread,'(a)') line
          if (line(1:1) /= '#') exit
       end do
       ! then get & check the numbers
       read(line,*) numbers
       if (abs(log(numbers(1)/xVals(iy))) > tolerance) then
          write(0,*) 'grid  x value = ', xVals(iy)
          write(0,*) 'input x value = ', numbers(1)
          call wae_error('read_photon', 'x value is inconsistent with grid')
       end if
       photon(iy) = numbers(6)
    end do
    
  end function read_photon
  
  
  !----------------------------------------------------------------------
  subroutine write_moments(pdf,label)
    real(dp),         intent(in) :: pdf(0:,ncompmin:)
    character(len=*), intent(in) :: label
    !-----------------
    real(dp) :: moments(ncompmin:ubound(pdf,dim=2))
    
    write(iunit,"(a)", advance="no"), "# total momentum "//trim(label)//" (& components)"
    moments = TruncatedMoment(grid, pdf, one)
    write(iunit,"(40f10.7)") sum(moments), moments

    write(iunit,"(a)", advance="no"), "# total number "//trim(label)//" (& components)"
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

    write(iunit,"(a)"), "# pdf "//trim(label)
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
