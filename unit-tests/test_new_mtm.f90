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
    type(new_mass_threshold_mat) :: mtm_from_P
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
    call InitMTMFromSplitMat(mtm_from_P, dh%allP(nloop,nf_light+1))

    call AllocPDF(grid, q)
    call AllocPDF(grid, P_q)
    call AllocPDF(grid, mtm_q)
    q = tables(0)%tab(:,:,0)
    P_q   = dh%allP(nloop,nf_light+1) * q
    mtm_q = mtm_from_P                * q

    do iflv = -nf_light-1, nf_light+1
      call check_approx_eq_1d("new_mass_threshold_mat check matrix element iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    call Delete(q)
    call Delete(P_q)
    call Delete(mtm_q)

    call Delete(mtm_from_P)

  end subroutine test_new_mass_threshold_mat
end module test_new_mtm