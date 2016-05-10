program test_qed
  use qed_coupling_module
  use types; use consts_dp
  implicit none
  type(qed_coupling) :: qed_alpha

  
  call InitQEDCoupling(qed_alpha, 1.5_dp, 4.5_dp, 173.0_dp)
  write(6,*) IQEDThresholdAtQ(qed_alpha, 0.2_dp)

  
  write(6,*) 'alpha QED at 0       = ', one/Value(qed_alpha, 0.0_dp)
  write(6,*) 'alpha QED at MZ      = ', one/Value(qed_alpha, 91.2_dp)
  write(6,*) 'alpha QED at 0.1 GeV = ', one/Value(qed_alpha, 0.1_dp)
  write(6,*) 'alpha QED at 1 GeV   = ', one/Value(qed_alpha, 1.0_dp)
  write(6,*) 'alpha QED at 10 GeV  = ', one/Value(qed_alpha, 10.0_dp)
  write(6,*) 'alpha QED at 100 GeV = ', one/Value(qed_alpha, 100.0_dp)
  write(6,*) 'alpha QED at 1 TeV   = ', one/Value(qed_alpha, 1000.0_dp)
  
  
end program test_qed
