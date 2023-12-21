#!/bin/bash

gfortran -c hplog5.f xc2ns3e.f xc2ns3p.f xc2ns3e_vs_xc2ns3p.f; gfortran *.o -o xc2ns3e_vs_xc2ns3p
