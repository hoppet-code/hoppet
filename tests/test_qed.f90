program test_qed
  use qed_coupling_module
  use types; use consts_dp
  implicit none
  type(qed_coupling) :: qed_alpha

  
  call InitQEDCoupling(qed_alpha, 1.5_dp, 4.5_dp, 173.0_dp)
  write(6,*) IQEDThresholdAtQ(qed_alpha, 0.2_dp)

  write(6,*) one/Value(qed_alpha, 91.2_dp)
  
  
end program test_qed
