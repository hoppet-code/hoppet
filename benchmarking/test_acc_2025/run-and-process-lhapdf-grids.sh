#!/usr/bin/env bash

usage="$0 buildDir refFile outputDir"
if [ "$#" -ne 3 ]; then
    echo "Usage: $usage"
    exit 1
fi
builddir=$1
ref=$2
outputdir=$3

mkdir -p $outputdir
logfile=$outputdir/LHAPDF-run.log
if [ -f $logfile ]; then
    echo "Warning: $logfile already exists, replacing it"
    mv $logfile $logfile.old
fi
touch $logfile
echo $* >> $logfile
date >> $logfile
echo "Running from: $(pwd)" >> $logfile
echo "Using build dir: $builddir" >> $logfile
echo "Using reference file: $ref" >> $logfile
echo "Using output dir: $outputdir" >> $logfile



for file in nloop3-preev-dy0.025 nloop3-preev-dy0.04 nloop3-preev-dy0.07 nloop3-preev-dy0.1 nloop3-preev-dy0.2 nloop3-ref-dy0.02 nloop3-preev-dy0.03 nloop3-preev-dy0.05 nloop3-preev-dy0.08 nloop3-preev-dy0.15 nloop3-preev-dy0.25 nloop3-ref-eps1e-9-dy0.02
do
    echo "------- Running PDFset $file through LHAPDF" | tee -a $logfile 
    $builddir/benchmarking/lhapdf_timings $file $outputdir >> $logfile 
    $builddir/benchmarking/compare2files_v2 $outputdir/${file}.lhapdf.dat $ref -summary -protect > $outputdir/$file.lhapdf.dat.acc.smry 2>&1 
done

