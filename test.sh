#!/bin/bash
#
# script for testing that things are OK before a release
#
# Usage: ./test.sh compiler-name

# generate makefiles
./configure $*

# just in case
make realclean

#
make

#
make check

# # build lib
# cd src/
# make realclean
# make -j2
# cd ..
# 
# # build and run test prog
# cd example_f90/
# if [[ -e tabulation_example.test_output ]] ; then
#     rm tabulation_example.test_output
# fi
#   
# make realclean
# make -j2
# ./tabulation_example > tabulation_example.test_output
# diffcmd="diff tabulation_example.test_output tabulation_example.default_output"
# res=`$diffcmd | wc -l`
# if [[ $res -eq 0 ]]
# then
#   echo 
#   echo '************ TEST PASSED OK ************'
# else
#   echo 
#   echo '************ TEST FAILED -- SEE DIFFERENCES BELOW ************'
#   echo 
#   echo '% ' $diffcmd
#   $diffcmd
#   echo 
#   echo '************ TEST FAILED -- SEE DIFFERENCES ABOVE ************'
# fi


