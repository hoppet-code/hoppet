#!/usr/bin/env bash
##
## Script to produce a a series of runs of varying accuracy, for easy
## comparisons of accuracy and timing.
##
## This script takes about 5 minutes on an M2Pro and produces
## around 50-90MB of output
##
## Reference results are stored in the 
## https://github.com/hoppet-code/2025-prec-and-timing
## repository

usage="$0 buildDir outputDir"
if [ "$#" -ne 2 ]; then
  echo "Usage: $usage"
  exit 1
fi
build_dir=$1
outdir=$2
nloop=4
checkinterp=no
#nloop=3
#checkinterp=yes


# Set color variables only if output is a terminal
if [ -t 1 ]; then
  term_red='\033[0;31m'
  term_blue='\033[0;34m'
  term_reset='\033[0m'
else
  term_red=''
  term_blue=''
  term_reset=''
fi



dyvals="0.25 0.2 0.15 0.1 0.08 0.07 0.05 0.04 0.03 0.025"
echo -e $term_blue"Running with nloop = $nloop"$term_reset
echo -e $term_blue"Output directory: $outdir"$term_reset
echo "======================"
accdir=$(dirname "$0")
#build_dir=$(dirname "$0")/../../build
compare="$build_dir/benchmarking/compare2files_v2"
common="$build_dir/benchmarking/prec_and_timing -outputgrid -nxQ 5000 -auto-nrep -nrep-eval 100 -nloop $nloop"

# a shortcut for running comamnd and checking output
# Usage: `runthis outname command args`
runthis() {
  outname=$1
  shift 
  $* -o $outname > $outname.log 2>&1 || \
    (echo -e $term_red ERROR when running$term_reset; \
    echo "> $*" ; \
    echo -e $term_red"tail of output in $outname.log:"$term_reset; \
    echo "======================" ; \
    tail -n 10 $outname.log;
    exit 1)
}

if [ ! -d $outdir ]; then
  mkdir $outdir
fi


dyref=0.02
fname_ref=$outdir/nloop$nloop-ref-dy$dyref.dat
echo "Getting a reference run with dy=$dyref ($fname_ref) and default eps"
runthis $fname_ref $common -dy $dyref -nopreev 


fname_ref_smalleps=$outdir/nloop$nloop-ref-eps1e-9-dy$dyref.dat
echo "Getting a reference run with dy=$dyref ($fname_ref_smalleps) and smaller eps"
runthis $fname_ref_smalleps $common -dy $dyref -nopreev -eps 1e-9

# tests with exact splitting functions
dyexact=0.1
echo "Running checks with exact NNLO splitting or threshold functions (dy=$dyexact)"
runthis $outdir/nloop$nloop-exactsp-nopreev-dy$dyexact.dat $common -dy $dyexact -nopreev -exact-nnlo-sp
runthis $outdir/nloop$nloop-exactspth-nopreev-dy$dyexact.dat $common -dy $dyexact -nopreev -exact-nnlo-sp -exact-nnlo-th

# doing runs to get accuracy with default interpolation settings
for dy in $dyvals
do 
  # get a basic "accuracy" run without pre-evolution
  echo -n Running with dy=$dy
  echo -n " nopreev"
  fname=$outdir/nloop$nloop-nopreev-dy$dy.dat
  $common -dy $dy -nopreev -o $fname  > $fname.log 2>&1
  $compare $fname $fname_ref -summary -protect > $fname.acc.smry 2>&1 
  #$compare $fname $fname_ref -channel 11 -protect > $fname.acc.ch11 2>&1

  # get a basic "accuracy" run with pre-evolution
  echo -n " preev"
  fname=$outdir/nloop$nloop-preev-dy$dy.dat
  $common -dy $dy -preev -o $fname  > $fname.log 2>&1
  $compare $fname $fname_ref -summary -protect > $fname.acc.smry 2>&1
  #$compare $fname $fname_ref -channel 11 -protect > $fname.acc.ch11 2>&1

  if [ x"$checkinterp" = x"yes" ]; then
    #for olnlnQ in 2 3 4
    for olnlnQ in 2 3
    do
      for oy in 2 3
      #for oy in $olnlnQ
      do
        echo -n " oQ$olnlnQ-oY$oy"
        fname=$outdir/nloop$nloop-preev-oQ$olnlnQ-oY$oy-dy$dy.dat
        runthis $fname $common -dy $dy -preev -olnlnQ $olnlnQ -yinterp-order $oy
        $compare $fname $fname_ref -summary -protect > $fname.acc.smry 2>&1
        #$compare $fname $fname_ref -channel 11 -protect > $fname.acc.ch11 2>&1
      done
    done
  fi
  echo    

done
