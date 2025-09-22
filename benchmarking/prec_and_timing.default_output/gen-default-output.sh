#!/usr/bin/env bash
# 
# Run this from the CMake build directory
#
# Usage PATH/gen-default-output.sh [-s]
#
# Options:
#
# -s: generate actual output (default is dry-run)
#
# If you change the list of tests, remember to make sure:
#
# - that the list is the same in the top-level CMakeLists.txt
# - that no unrequired files are left in this directory


runs=(\
  benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:2:-yinterp-order:2 \
  benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:3:-yinterp-order:3 \
  benchmarking/prec_and_timing:-dy:0.2:-olnlnQ:4:-yinterp-order:4 \
  benchmarking/prec_and_timing:-dy:0.1:-exact-nnlo-th:-exact-nnlo-sp \
)

dummy=yes
# if the user provides a -s argument, then set dummy to no 
# and generate actual output
if [[ "$1" == "-s" ]]; then
    dummy=no
fi

# Colors for output
blue="\033[0;34m"
yellow="\033[0;33m"
red="\033[0;31m"
green="\033[0;32m"
bold="\033[1m"
reset="\033[0m"

scriptdir=$(dirname "$(realpath "$0")")
echo "Output files will go to $scriptdir"
extraargs="-o /dev/null -output-benchmark"

for run in "${runs[@]}"
do 
    echo "============================================================="
    echo -e "${blue}${bold}Preparing $run${reset}"
    testcmd=$(echo "$run" | sed -e 's/:/ /g')
    testcmd="${testcmd} ${extraargs}"
    outfile=$(echo "$run" | sed -e 's/:/_/g' -e s'|benchmarking/||').default_output
    echo -e "Command is: ${bold}$testcmd${reset}"

    if [[ "$dummy" == "no" ]]; then
        echo -e "Output saved to: $scriptdir/$outfile"
        echo "=============================================================\n"
        $testcmd | tee $scriptdir/$outfile 
    else
        echo -e "${bold}Dry-run${reset}, output would be saved to: $scriptdir/$outfile"
        echo "=============================================================\n"
        #$testcmd
    fi
done

if [[ "$dummy" == "no" ]]; then
    echo -e "${green}${bold}Output saved...${reset}"
else
    echo -e "${red}${bold}Dry run mode: no output was saved (re-run with -s to save output).${reset}"
fi