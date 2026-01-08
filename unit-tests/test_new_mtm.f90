module test_new_mtm
  use hoppet
  use hoppet_mass_threshold_mat
  use unit_tests
  implicit none

  private

  public :: test_mass_threshold_mat

contains

  subroutine test_mass_threshold_mat()
    use streamlined_interface
    type(mass_threshold_mat) :: mtm_from_P, mtm_from_mtm_p, mtm_from_p_mtm, mtm_tmp
    type(split_mat) :: tmp_split_mat, Pnf3, Pnf4a, Pnf4b
    type(grid_conv) :: test_conv
    integer :: nf_light, nloop, iflv, ix, iprobe
    real(dp), pointer :: P_q(:,:), mtm_q(:,:), q(:,:), tmp(:,:), probes(:,:,:)
    real(dp) :: xvals(3) = [0.01_dp, 0.1_dp, 0.4_dp], x
    !type(grid_conv) :: delta !! for diagnostics

    if (.not. do_test("test_mtm")) return

    nf_light = 3
    nloop    = dh%nloop
    if (nloop > dh%nloop .or. nloop < 2) call wae_error(&
          "test_mass_threshold_mat: nloop(="//trim(to_string(nloop))//&
          ") exceeds dh%nloop(="//trim(to_string(dh%nloop))//") or is below 2")

    !call InitGridConv(dh%grid, delta, delta_fn)

    ! set up some "dummy" splitting matrices that have all the structure we need
    ! Pnf4a (4-flavour) = P_4 + 0.36 * P_3
    call InitSplitMat(Pnf4a, dh%allP(nloop,nf_light+1))
    call AddWithCoeff(Pnf4a, dh%allP(nloop-1,nf_light+1), 0.36_dp)
    ! Pnf4b (4-flavour) = P_4 + pi * P_3
    call InitSplitMat(Pnf4b, dh%allP(nloop,nf_light+1))
    call AddWithCoeff(Pnf4b, dh%allP(nloop-1,nf_light+1), pi)
    ! Pnf3a (3-flavour) = P_4 + P_3 (3-flav)
    call InitSplitMat(Pnf3 , dh%allP(nloop,nf_light))
    call AddWithCoeff(Pnf3 , dh%allP(nloop-1,nf_light))


    call AllocPDF(grid, q)
    call AllocPDF(grid, P_q)
    call AllocPDF(grid, mtm_q)
    q = tables(0)%tab(:,:,0)

    ! first check that mtm_from_P works the same as just multiplying by the corresponding P
    call InitMTM(mtm_from_P, Pnf4a)
    P_q   = Pnf4a      * q
    mtm_q = mtm_from_P * q
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (mtm from sm) iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    ! check multiplication
    call InitMTM(mtm_tmp, mtm_from_P)
    call Multiply(mtm_tmp, two)
    P_q   = two * (mtm_from_p * q)
    mtm_q = mtm_tmp * q
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (multiply) iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    ! check initialisation to zero & addition of MTM
    call InitMTM(grid, mtm_tmp, nf_light)
    call AddWithCoeff(mtm_tmp, mtm_from_P, 1.5_dp)
    P_q   = 1.5_dp * (mtm_from_p * q)
    mtm_q = mtm_tmp * q
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (zero & add mtm) iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    ! check addition of nf=3 split mat
    call InitMTM(mtm_tmp, mtm_from_P)
    call AddWithCoeff(mtm_tmp, Pnf3, 0.5_dp)
    P_q   = mtm_from_p * q + 0.5_dp * (Pnf3 * q)
    mtm_q = mtm_tmp * q
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (add Pnf3) iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    ! check addition of nf=4 split mat
    call InitMTM(mtm_tmp, mtm_from_P)
    call AddWithCoeff(mtm_tmp, Pnf4b, 0.5_dp)
    P_q   = mtm_from_p * q + 0.5_dp * (Pnf4b * q)
    mtm_q = mtm_tmp * q
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (add Pnf4b) iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    ! grid locking breaks the exact numerical correspondence (within machine prec) between 
    ! (P_A * P_B) * q and P_A * (P_B * q), so we temporarily disable it for the next few tests
    call DisableGridLocking()

    ! now check that we can multiply an mtm and a splitting matrix correctly
    call SetToConvolution(mtm_from_mtm_p, mtm_from_P, Pnf3)
    P_q = mtm_from_P * (Pnf3 * q)
    mtm_q = mtm_from_mtm_p * q
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (mtm_from_P * Pnf3) element iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    !! for diagnostics
    !call Multiply(mtm_from_P, zero)
    !call Multiply(Pnf4b, zero)
    !call InitGridConv(mtm_from_P%PShg, delta)
    !call InitGridConv(Pnf4b%qq,delta)
    !q(:,:)=0
    !q(0,0)=1
    
    call SetToConvolution(mtm_from_p_mtm, Pnf4b, mtm_from_P)
    P_q =  Pnf4b * (mtm_from_P * q)
    call AllocPDF(grid,tmp)
    tmp = mtm_from_P * q
    P_q = Pnf4b * tmp
    mtm_q = mtm_from_p_mtm * q
    !! useful printouts if doing the kinds of diagnostics
    !! set up above
    !write(6,*) "q    :", real(q  (0,-4:4))
    !write(6,*) "tmp  :", real(tmp(0,-4:4))
    !write(6,*) "P_q  :", real(P_q  (0,-4:4))
    !write(6,*) "mtm_q:", real(mtm_q(0,-4:4))
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (Pnf4b * mtm_from_P) iflv="//trim(to_string(iflv)), &
           answer=mtm_q(:,iflv), expected=P_q(:,iflv), tol_abs=1.0e-8_dp, tol_rel=1.0e-10_dp, tol_choice_or=.true.)
    end do
    ! restore grid locking
    call RestoreGridLocking()

    ! we'll want to test derived probes with MTMs, but first do it with
    ! splitting matrices to make sure the probes are correct
    call GetDerivedSplitMatProbes(dh%grid, nf_light, probes)
    do iprobe = 1, size(probes,dim=3)
      call CopyEvlnPdfToHuman(nf_light, probes(:,:,iprobe), tmp)
      tmp = Pnf3 * (Pnf3 * tmp)
      call CopyHumanPdfToEvln(nf_light, tmp, probes(:,:,iprobe))
    end do
    call InitSplitMat(grid, tmp_split_mat, nf_light)
    call SetDerivedSplitMat(tmp_split_mat, probes)
    ! and then test it
    call DisableGridLocking()
    P_q = Pnf3 * (Pnf3 * q)
    tmp = tmp_split_mat * q
    call RestoreGridLocking()
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (SetDerivedSplitMat) iflv="//trim(to_string(iflv)), &
           answer=tmp(:,iflv), expected=P_q(:,iflv), tol_abs=1.0e-8_dp, tol_rel=1.0e-10_dp, tol_choice_or=.true.)
    end do

    ! now test SetDerivedMTM
    call GetDerivedMTMProbes(dh%grid, nf_light, probes)
    do iprobe = 1, size(probes,dim=3)
      probes(:,:,iprobe) = Pnf4b * (mtm_from_P * probes(:,:,iprobe))
      !probes(:,:,iprobe) = Pnf3 * probes(:,:,iprobe)
    end do
    call InitMTM(grid, mtm_tmp, nf_light)
    call SetDerivedMTM(mtm_tmp, probes)
    call DisableGridLocking()
    !P_q = Pnf3 * q
    P_q = Pnf4b * (mtm_from_P * q)
    mtm_q = mtm_tmp * q
    call RestoreGridLocking()
    do iflv = iflv_tbar, iflv_t
      call check_approx_eq_1d("mass_threshold_mat (SetDerivedMTM) iflv="//trim(to_string(iflv)), &
           answer=mtm_q(:,iflv), expected=P_q(:,iflv), tol_abs=1.0e-8_dp, tol_rel=1.0e-10_dp, tol_choice_or=.true.)
    end do

    ! next, check that we can get a new mtm from an old one and that they act in the same way
    if (nloop >= 3) then
      call InitMTM(mtm_tmp, dh%allMTM(nloop,nf_light+1))
      P_q = dh%allMTM(nloop,nf_light+1) * q
      mtm_Q = (mtm_tmp * q)
      do iflv = iflv_tbar, iflv_t
        call check_approx_eq_1d("mass_threshold_mat (new from old) iflv="//trim(to_string(iflv)), &
             answer=mtm_q(:,iflv), expected=P_q(:,iflv), tol_abs=1.0e-8_dp, tol_rel=1.0e-10_dp, tol_choice_or=.true.)
      end do
    end if


    call Delete(q)
    call Delete(P_q)
    call Delete(mtm_q)
    call Delete(tmp)

    call Delete(mtm_from_P)
    call Delete(mtm_from_mtm_p)
    call Delete(mtm_from_p_mtm)
    call Delete(mtm_tmp)

    call Delete(Pnf3)
    call Delete(Pnf4a)
    call Delete(Pnf4b)
    call Delete(tmp_split_mat)
  end subroutine test_mass_threshold_mat


  function delta_fn(y) result(res)
    use convolution_communicator
    real(dp), intent(in) :: y
    real(dp)             :: res
    if (cc_piece == cc_delta) then
      res = one
    else
      res = zero
    end if
  end function delta_fn
end module test_new_mtm