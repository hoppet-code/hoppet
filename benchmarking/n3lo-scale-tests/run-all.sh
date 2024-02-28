#!/usr/bin/env bash
#
# Run the evolution_n3lo_tests program for different scales and alpha_s values
# so as to test that we do not break the N3LO accuracy when changing x_muR.

# location of the executable relative to the directory of this script
execdir=`dirname $0`/../

$execdir/evolution_n3lo_tests -nloop 4 -Q 100000000 -Q0 1 -asQ0 0.0125 -o nloop4-asQ0-0.0125-range100000000.dat
$execdir/evolution_n3lo_tests -nloop 4 -Q 10000     -Q0 1 -asQ0 0.025  -o nloop4-asQ0-0.025-range10000.dat
$execdir/evolution_n3lo_tests -nloop 4 -Q 100       -Q0 1 -asQ0 0.05   -o nloop4-asQ0-0.05-range100.dat
$execdir/evolution_n3lo_tests -nloop 4 -Q 10        -Q0 1 -asQ0 0.1    -o nloop4-asQ0-0.1-range10.dat