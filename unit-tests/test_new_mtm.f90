module test_new_mtm
  use hoppet
  use dglap_objects_new_mtm
  use unit_tests
  implicit none

  private

  public :: test_new_mass_threshold_mat

contains

  subroutine test_new_mass_threshold_mat()
    use streamlined_interface
    type(new_mass_threshold_mat) :: mtm_from_P, mtm_from_mtm_p
    type(split_mat) :: tmp_split_mat
    type(grid_conv) :: test_conv
    integer :: nf_light, nloop, iflv, ix
    real(dp), pointer :: P_q(:,:), mtm_q(:,:), q(:,:)
    real(dp) :: xvals(3) = [0.01_dp, 0.1_dp, 0.4_dp], x

    if (.not. do_test("new_mtm")) return
    
    nf_light = 3
    nloop   = 4
    if (nloop > dh%nloop) call wae_error(&
          "test_new_mass_threshold_mat: nloop(="//trim(to_string(nloop))//&
          ") exceeds dh%nloop(="//trim(to_string(dh%nloop))//")")

    call AllocPDF(grid, q)
    call AllocPDF(grid, P_q)
    call AllocPDF(grid, mtm_q)

    ! first check that mtm_from_P works the same as just multiplying by the corresponding P
    call InitMTMFromSplitMat(mtm_from_P, dh%allP(nloop,nf_light+1))
    q = tables(0)%tab(:,:,0)
    P_q   = dh%allP(nloop,nf_light+1) * q
    mtm_q = mtm_from_P                * q

    do iflv = -nf_light-1, nf_light+1
      call check_approx_eq_1d("new_mass_threshold_mat check matrix element iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do


    ! grid locking breaks the exact numerical correspondence (within machine prec) between 
    ! (P_A * P_B) * q and P_A * (P_B * q), so we temporarily disable for the next few tests
    call DisableGridLocking()

    ! now check that we can multiply an mtm and a splitting matrix correctly
    !q(:,2:6) = 0.0_dp
    !q(:,-6:0) = 0.0_dp
    !call SetToConvolution(tmp_split_mat, dh%allP(nloop-1,nf_light), dh%allP(nloop,nf_light))
    !P_q = dh%allP(nloop-1,nf_light) * (dh%allP(nloop,nf_light) * q)
    !mtm_q = tmp_split_mat * q

    call SetToConvolution(mtm_from_mtm_p, mtm_from_P, dh%allP(3,nf_light))
    P_q = mtm_from_P * (dh%allP(3,nf_light) * q)
    mtm_q = mtm_from_mtm_p * q
    do iflv = -nf_light-1, nf_light+1
      call check_approx_eq_1d("new_mass_threshold_mat check matrix element iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    ! restore grid locking
    call RestoreGridLocking()

    call Delete(q)
    call Delete(P_q)
    call Delete(mtm_q)

    call Delete(mtm_from_P)
    !call Delete(mtm_from_mtm_p)
  end subroutine test_new_mass_threshold_mat
end module test_new_mtm