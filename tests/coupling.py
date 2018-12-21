#!/usr/bin/env python

from __future__ import print_function, division
import sys

sys.path.append('/Users/gsalam/work/alphas')
import alphas as al

mu_vals = [1.1, 2.0, 6.0, al.pdg_mZ, 1000.0 ]
coupling = al.alphas(0.118, nloop = 4)

print (coupling.thresholds())
for mu in mu_vals:
    print(mu, coupling(mu))

msbar_thresholds = al.thresholds(
    *[al.mass_threshold(al.pdg_mc_MSbarM, al.MSbar),
      al.mass_threshold(al.pdg_mb_MSbarM, al.MSbar),
      al.mass_threshold(161.0,            al.MSbar)])

print()
msbar_coupling = al.alphas(0.118, thresholds = msbar_thresholds, nloop = 4)
print (msbar_coupling.thresholds())
for mu in mu_vals:
    print(mu, msbar_coupling(mu))
