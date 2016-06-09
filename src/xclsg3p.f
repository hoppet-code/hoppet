      MODULE XCLSG3P
      USE XC2NS3P
      CONTAINS
*
* ..File: xclsg3p.f    F_L^PS  and  F_L^G
*
*
* ..Parametrizations of the 3-loop MS(bar) pure-singlet and gluon coef-
*    ficient functions for the electromagnetic structure function F_L
*    at  mu_r = mu_f = Q. The expansion parameter is  alpha_s/(4 pi).
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is one part in thousand or better.
*
* ..Reference: S. Moch, J. Vermaseren and A. Vogt,
*              hep-ph/0411112 = Phys. Lett. B606 (2005) 123
* 
* =====================================================================
*
*
* ..The pure-singlet coefficient function
*
       FUNCTION CLS3A (Y, DL, NF, CC)
       IMPLICIT REAL*8 (A-Z)
       DIMENSION FLS(6)
       INTEGER CC ! charged current
       INTEGER NF
       DATA FLS / 1.d0, 0.1d0, 0.d0, 0.1d0, 0.01818181818d0, 0.1d0 /
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
       D27 = 1./27.D0
       D81 = 1./81.D0
       FLS11 = FLS(NF)
       CLS3A = 0.D0
*
       IF (CC.EQ.1) THEN
       CLS31 = (1568.*D27 * DL1**3 - 11904.*D27 * DL1**2 + 5124.* DL1)
     ,           * Y1**2  +  DL*DL1 * (2184.* DL + 6059.* Y1)
     ,         - (795.6 + 1036.* Y) * Y1**2  - 143.6 * DL*Y1
     ,         + 8544.*D27 * DL**2 - 1600.*D27 * DL**3 
     ,         - 885.53 /Y *Y1**2 - 182.00 * DL/Y * Y1
       CLS32 = ( - 96.*D27 * DL1**2 + 29.52 * DL1) * Y1**2  + 
     ,         +  DL*DL1 * (35.18 * DL + 73.06 * Y1) 
     ,         - (14.16 - 69.84 * Y) * Y1**2 - 35.24 * Y*DL**2 
     ,         - 69.41 * DL*Y1 - 384.*D27 * DL**2 + 40.239 /Y *Y1**2 
       CLS3A = NF * ( CLS31 + NF * CLS32 )
       ELSE
       CLS3F = ( (107.0 + 321.05 * Y - 54.62 * Y**2) *(1.-Y)
     ,         - 26.717 - 320.*D81 * DL**3 - 640.*D81 * DL**2
     ,         + 9.773 * DL + Y*DL * (363.8 + 68.32 * DL) ) * Y
       CLS3A = NF * FLS11 * CLS3F
       ENDIF   
*     
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The gluon coefficient function
*
       FUNCTION CLG3A (Y, DL, NF, CC)
       IMPLICIT REAL*8 (A-Z)
       DIMENSION FLG(6)
       INTEGER CC ! charged current
       INTEGER NF
       DATA FLG / 1.d0, 0.1d0, 0.d0, 0.1d0, 0.01818181818d0, 0.1d0 /
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
       D27 = 1./27.D0
       FLG11 = FLG(NF)
*
       CLG3A = 0.D0
       IF (CC.EQ.1) THEN
       CLG31 = (144.D0 * DL1**4 - 47024.*D27 * DL1**3 + 6319.* DL1**2 
     ,         + 53160.* DL1) * Y1 + DL*DL1 * (72549. + 88238.* DL)
     ,         + (3709. - 33514.* Y - 9533.* Y**2) * Y1 
     ,         + 66773.* Y*DL**2 - 1117.* DL + 45.37 * DL**2 
     ,         - 5360.*D27 * DL**3 - 2044.70 /Y*Y1 - 409.506 * DL/Y
       CLG32 = (288.*D27 * DL1**3 - 3648.*D27 * DL1**2 - 592.3 * DL1
     ,         + 1511.* Y*DL1) * Y1 + DL*DL1 * (311.3 + 14.24 * DL)
     ,         + (577.3 - 729.0 * Y) * Y1 + 30.78 * Y*DL**3
     ,         + 366.0 * DL + 3000.*D27 * DL**2 + 480.*D27 * DL**3
     ,         + 88.5037 /Y*Y1
       CLG3A = NF * ( CLG31 + NF * CLG32  )
       ELSE
       CLG3F = (-0.0105 * DL1**3 + 1.550 * DL1**2 + 19.72 *Y*DL1 
     ,         - 66.745 * Y + 0.615 * Y**2) * Y1 + 20.*D27 * Y*DL**4
     ,         + (280./81.D0 + 2.260* Y) * Y*DL**3 - (15.40 - 2.201* Y)
     ,           * Y*DL**2 - (71.66 - 0.121 * Y) * Y*DL   
       CLG3A = NF * NF * FLG11 * CLG3F
       ENDIF
*     
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XCLSG3P
