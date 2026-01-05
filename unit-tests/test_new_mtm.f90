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
    type(new_mass_threshold_mat) :: mtm_from_P, mtm_from_mtm_p, mtm_from_p_mtm
    type(split_mat) :: tmp_split_mat, Pnf3, Pnf4a, Pnf4b
    type(grid_conv) :: test_conv
    integer :: nf_light, nloop, iflv, ix
    real(dp), pointer :: P_q(:,:), mtm_q(:,:), q(:,:), tmp(:,:)
    real(dp) :: xvals(3) = [0.01_dp, 0.1_dp, 0.4_dp], x

    type(grid_conv) :: delta

    if (.not. do_test("new_mtm")) return
    


    nf_light = 3
    nloop   = 4
    if (nloop > dh%nloop) call wae_error(&
          "test_new_mass_threshold_mat: nloop(="//trim(to_string(nloop))//&
          ") exceeds dh%nloop(="//trim(to_string(dh%nloop))//")")

    call InitGridConv(dh%grid, delta, delta_fn)

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

    ! first check that mtm_from_P works the same as just multiplying by the corresponding P
    call InitMTMFromSplitMat(mtm_from_P, Pnf4a)

    q = tables(0)%tab(:,:,0)
    P_q   = Pnf4a      * q
    mtm_q = mtm_from_P * q

    do iflv = -nf_light-1, nf_light+1
      call check_approx_eq_1d("new_mass_threshold_mat check matrix element iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do


    ! grid locking breaks the exact numerical correspondence (within machine prec) between 
    ! (P_A * P_B) * q and P_A * (P_B * q), so we temporarily disable for the next few tests
    call DisableGridLocking()

    ! now check that we can multiply an mtm and a splitting matrix correctly
    !call SetToConvolution(tmp_split_mat, dh%allP(nloop-1,nf_light), dh%allP(nloop,nf_light))
    !P_q = dh%allP(nloop-1,nf_light) * (dh%allP(nloop,nf_light) * q)
    !mtm_q = tmp_split_mat * q

    call SetToConvolution(mtm_from_mtm_p, mtm_from_P, Pnf3)
    P_q = mtm_from_P * (Pnf3 * q)
    mtm_q = mtm_from_mtm_p * q
    do iflv = -nf_light-1, nf_light+1
      call check_approx_eq_1d("new_mass_threshold_mat (mtm_from_P * Pnf3) element iflv="//trim(to_string(iflv)), &
           mtm_q(:,iflv), P_q(:,iflv), 1.0e-10_dp, 1.0e-10_dp, tol_choice_or=.true.)
    end do

    !call Multiply(mtm_from_P, zero)
    !call Multiply(Pnf4b, zero)
    !call InitGridConv(mtm_from_P%PShg, delta)
    !call InitGridConv(Pnf4b%qq,delta)
    !q(:,:)=0
    !q(0,0)=1
    
    ! and next, check that we convolute a splitting matrix and an mtm correctly
    !q(:,0) = 0.0_dp
    !q(:,+2:+6:+1) = 0.0_dp
    !q(:,-2:-6:-1) = 0.0_dp
    !q(:, 1) =  q(:,1)
    !q(:, 2) =  q(:,1)
    !q(:,-3:-1) = -q(:,1:3)

    
    !q(:,1) = q(:,-1)
    !q(:,2) = 0.0*q(:,-1)
    !q(:,3) = 0.0*q(:,-1)
    !q(:,-3:-1) = q(:,3:1:-1)
    !q(:,-6:0) = 0.0_dp
    call SetToConvolution(mtm_from_p_mtm, Pnf4b, mtm_from_P)
    !write(0,*) mtm_from_p_mtm%PShq%subgc(1)%conv
    !P_q =  Pnf4b * (mtm_from_P * q)
    call AllocPDF(grid,tmp)
    tmp = mtm_from_P * q
    !tmp(:,+4)=0.0_dp
    !tmp(:,-4)=0.0_dp
    P_q = Pnf4b * tmp
    mtm_q = mtm_from_p_mtm * q
    write(6,*) "q    :", real(q  (0,-4:4))
    write(6,*) "tmp  :", real(tmp(0,-4:4))
    write(6,*) "P_q  :", real(P_q  (0,-4:4))
    write(6,*) "mtm_q:", real(mtm_q(0,-4:4))
    !do iflv = -4,4!-4,4,8!-nf_light, nf_light
    do iflv = -nf_light-1, nf_light+1
      call check_approx_eq_1d("new_mass_threshold_mat check matrix element iflv="//trim(to_string(iflv)), &
           answer=mtm_q(:,iflv), expected=P_q(:,iflv), tol_abs=1.0e-8_dp, tol_rel=1.0e-10_dp, tol_choice_or=.true.)
           !answer=mtm_q(:,iflv)+mtm_q(:,-iflv), expected=P_q(:,iflv)+P_q(:,-iflv), tol_abs=1.0e-8_dp, tol_rel=1.0e-10_dp, tol_choice_or=.true.)
    end do

    !write(6,*) mtm_q(:,iflv)+mtm_q(:,-iflv)

    ! restore grid locking
    call RestoreGridLocking()

    call Delete(q)
    call Delete(P_q)
    call Delete(mtm_q)

    call Delete(mtm_from_P)
    call Delete(mtm_from_mtm_p)
    call Delete(mtm_from_p_mtm)
  end subroutine test_new_mass_threshold_mat


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