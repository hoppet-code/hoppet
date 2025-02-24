#!/bin/bash

dirname=../with-lhapdf/
prog=$dirname/tabulation_n3lo_v130_paper


$prog -o MSHT20nnlo_as118-n3lo-evolu-n3lo-mass -nloop 4 -n3lo-vfns-on 1 -pdf MSHT20nnlo_as118 > MSHT20nnlo_as118-n3lo-evolu-n3lo-mass.log &

exit 0

$prog -o MSHT20nnlo_as118-n3lo-evolu-nnlo-mass -nloop 4 -n3lo-vfns-on 0 -pdf MSHT20nnlo_as118 > MSHT20nnlo_as118-n3lo-evolu-nnlo-mass.log &

$prog -o MSHT20nnlo_as118-nnlo-evolu-nnlo-mass -nloop 3 -n3lo-vfns-on 0 -pdf MSHT20nnlo_as118 > MSHT20nnlo_as118-nnlo-evolu-nnlo-mass.log &

$prog -o NNPDF40_nnlo_pch_as_01180-n3lo-evolu-n3lo-mass -nloop 4 -n3lo-vfns-on 1 -pdf NNPDF40_nnlo_pch_as_01180 > NNPDF40_nnlo_pch_as_01180-n3lo-evolu-n3lo-mass.log &

$prog -o NNPDF40_nnlo_pch_as_01180-n3lo-evolu-nnlo-mass -nloop 4 -n3lo-vfns-on 0 -pdf NNPDF40_nnlo_pch_as_01180 > NNPDF40_nnlo_pch_as_01180-n3lo-evolu-nnlo-mass.log &

$prog -o NNPDF40_nnlo_pch_as_01180-nnlo-evolu-nnlo-mass -nloop 3 -n3lo-vfns-on 0 -pdf NNPDF40_nnlo_pch_as_01180 > NNPDF40_nnlo_pch_as_01180-nnlo-evolu-nnlo-mass.log &

$prog -o toyHERA-n3lo-evolu-n3lo-mass -nloop 4 -n3lo-vfns-on 1 > toyHERA-n3lo-evolu-n3lo-mass.log&

$prog -o toyHERA-n3lo-evolu-nnlo-mass -nloop 4 -n3lo-vfns-on 0 > toyHERA-n3lo-evolu-nnlo-mass.log &

$prog -o toyHERA-nnlo-evolu-nnlo-mass -nloop 3 -n3lo-vfns-on 0 > toyHERA-nnlo-evolu-nnlo-mass.log
