#!/usr/bin/env bash

usage="$0 buildDir ref outputDir"
if [ "$#" -ne 3 ]; then
  echo "Usage: $usage"
  exit 1
fi
builddir=$1
ref=$2
outputdir=$3

for file in nloop3-preev-dy0.025 nloop3-preev-dy0.04 nloop3-preev-dy0.07 nloop3-preev-dy0.1 nloop3-preev-dy0.2 nloop3-ref-dy0.02 nloop3-preev-dy0.03 nloop3-preev-dy0.05 nloop3-preev-dy0.08 nloop3-preev-dy0.15 nloop3-preev-dy0.25 nloop3-ref-eps1e-9-dy0.02
do
    $builddir/example_cpp/lhapdf_timings $file $outputdir
    $buildir/benchmarking/compare2files_v2 ${file}.dat $ref -summary -protect > $outputdir/$file.acc.smry 2>&1 
done

