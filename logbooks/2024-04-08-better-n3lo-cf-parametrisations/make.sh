#!/bin/bash

gfortran -c hplog5.f xc2ns3e.f  xc2ns3p.f  xc3ns3e.f  xc3ns3p.f  xclns3e.f  xclns3p.f N3LOCNS.f c2nsreg.f print_coefficient_functions_for_fit.f; gfortran *.o -o print_coefficient_functions_for_fit
