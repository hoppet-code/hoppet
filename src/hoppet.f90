module hoppet
  use types; use consts_dp
  use convolution; use dglap_objects
  use pdf_representation; use pdf_general
  use dglap_choices
  use warnings_and_errors
  use dglap_holders 
  use evolution
  use qcd_coupling
  use qcd
  use pdf_tabulate

  ! extra QED things
  use qed_evolution
  use qed_objects
  use qed_coupling_module

  use structure_functions
  use coefficient_functions_holder_exact
  implicit none
end module hoppet
