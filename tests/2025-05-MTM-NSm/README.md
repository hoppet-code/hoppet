Tests of mass thresholds at N3LO, in particular the effect of the non-singlet minus
piece (post v. pre-NSm) around the b-threshold.

## Sum rule checks, across charm threshold, toy PDF

See `sumrule-Qmax1.42-postNSm.dat`, `sumrule-Qmax1.42-preNSm.dat`

prior to inclusion of NSm part in the threshold:
```
 Momentum sum rule =    1.0000000789329637      , expected 1, difference =    7.8932963720745875E-008
 uv + dv sum rule =    2.9997982067332973      , expected 3, difference =   -2.0179326670266562E-004
 uv - dv sum rule =   0.99993273630074975      , expected 1, difference =   -6.7263699250252884E-005
```

after inclusion of NSm part in the threshold:
```
 Momentum sum rule =    1.0000000789329637      , expected 1, difference =    7.8932963720745875E-008
 uv + dv sum rule =    2.9999999580394734      , expected 3, difference =   -4.1960526608875170E-008
 uv - dv sum rule =   0.99999998693845016      , expected 1, difference =   -1.3061549841708597E-008
```

## Full PDF checks, comparing to NNPDF40_an3lo_as_01180

All of these checks take the NNPDF LHAPDF initial condition just below
(4.919) the b threshold (4.920), evolve just above (4.921) and then
study various flavours.

See [MTM-NSm-test-plots.pdf](MTM-NSm-test-plots.pdf)

- p.5 shows the effect of the NS-minus piece on the down, up and strange
  valence contributions, showing the "post-NSm" versus "pre-NSm"
  difference just above the b threshold. The effect is at the few times
  $10^{-6}$ level

- p.2 shows that there seems to be a small-x issue in NNPDF's treatment
  of the b threshold (at NNLO things were fine). At $x=10^{-4}$
  NNPDF=$-0.135$ whereas we find $-0.08$ for $b+\bar b$. For comparison
  at 50GeV, $b+\bar b\simeq3.5$, so we're talking a likely ballpark of
  1.5% there.

- any other NNPDF "issues" look like they are small. To within accuracy
  we cannot reliably tell whether they have the NSminus piece in the
  mass thresholds.
