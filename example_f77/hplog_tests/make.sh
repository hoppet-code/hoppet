#!/bin/bash

gfortran -c *.f; gfortran *.o -o xc2ns3e_vs_xc2ns3p
