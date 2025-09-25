# CMake generated Testfile for 
# Source directory: /home/runner/work/hoppet/hoppet
# Build directory: /home/runner/work/hoppet/hoppet/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(tabulation_example "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check example_f90/tabulation_example")
set_tests_properties(tabulation_example PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;435;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(tabulation_example_n3lo_takes_up_to_30s "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check example_f90/tabulation_example_n3lo")
set_tests_properties(tabulation_example_n3lo_takes_up_to_30s PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;435;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(tabulation_example_qed "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check example_f90/tabulation_example_qed")
set_tests_properties(tabulation_example_qed PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;435;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(structure_functions_example "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check example_f90/structure_functions_example")
set_tests_properties(structure_functions_example PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;435;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(structure_functions_example_flavour "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check example_f90/structure_functions_example_flavour")
set_tests_properties(structure_functions_example_flavour PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;435;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(tabulation_example_qed_streamlined "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check example_f90/tabulation_example_qed_streamlined")
set_tests_properties(tabulation_example_qed_streamlined PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;435;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(tabulation_example_cpp "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check example_cpp/tabulation_example")
set_tests_properties(tabulation_example_cpp PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;441;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:2:-yinterp-order:2 "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:2:-yinterp-order:2")
set_tests_properties(benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:2:-yinterp-order:2 PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;479;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:3:-yinterp-order:3 "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:3:-yinterp-order:3")
set_tests_properties(benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:3:-yinterp-order:3 PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;479;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:4:-yinterp-order:4 "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:4:-yinterp-order:4")
set_tests_properties(benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:4:-yinterp-order:4 PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;479;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
add_test(benchmarking/prec_and_timing:-dy:0.1:-exact-nnlo-th:-exact-nnlo-sp "sh" "-c" "/home/runner/work/hoppet/hoppet/scripts/check benchmarking/prec_and_timing:-dy:0.1:-exact-nnlo-th:-exact-nnlo-sp")
set_tests_properties(benchmarking/prec_and_timing:-dy:0.1:-exact-nnlo-th:-exact-nnlo-sp PROPERTIES  _BACKTRACE_TRIPLES "/home/runner/work/hoppet/hoppet/CMakeLists.txt;479;add_test;/home/runner/work/hoppet/hoppet/CMakeLists.txt;0;")
subdirs("unit-tests")
subdirs("pyinterface")
