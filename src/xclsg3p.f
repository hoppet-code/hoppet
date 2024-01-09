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
       DIMENSION FL(6), FLS(6)
       INTEGER CC ! charged current
       INTEGER NF
       DATA FL  / -1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
       DATA FLS / 1.d0, 0.1d0, 0.d0, 0.1d0, 0.01818181818d0, 0.1d0 /
       DOUBLE PRECISION D27, D81
       PARAMETER (D27 = 1.D0/27.D0, D81 = 1.D0/81.D0)
*
       OY = 1.0D0/Y
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
       FLS11 = FLS(NF) - FL(NF)  
       CLS3A = 0.D0
*
       IF (CC.EQ.1) THEN
          CLS31 = (1568.D0*D27 * DL1**3 - 11904.D0*D27 * DL1**2 +
     $         5124.D0* DL1)* Y1**2  +  DL*DL1 * (2184.D0* DL + 6059.D0*
     $         Y1)- (795.6D0 + 1036.D0* Y) * Y1**2  - 143.6D0 * DL*Y1+
     $         8544.D0*D27 * DL**2 - 1600.D0*D27 * DL**3- 885.53D0 *OY
     $         *Y1**2 - 182.00D0 * DL*OY * Y1
       CLS32 = ( - 96.D0*D27 * DL1**2 + 29.52D0 * DL1) * Y1**2  +  DL
     $      *DL1 * (35.18D0 * DL + 73.06D0 * Y1) - (14.16D0 - 69.84D0 *
     $      Y) * Y1**2 - 35.24D0 * Y*DL**2 - 69.41D0 * DL*Y1 - 384.D0
     $      *D27 * DL**2 + 40.239D0 *OY *Y1**2 
       CLS3A = NF * ( CLS31 + NF * CLS32 )
       ELSE
          CLS3F = ( (107.0D0 + 321.05D0 * Y - 54.62D0 * Y**2) *(1.D0-Y)
     $         - 26.717D0 - 320.D0*D81 * DL**3 - 640.D0*D81 * DL**2 +
     $         9.773D0 * DL + Y*DL * (363.8D0 + 68.32D0 * DL) ) * Y
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
       DOUBLE PRECISION D27, D81
       PARAMETER (D27 = 1.D0/27.D0, D81 = 1.D0/81.D0)
*
       OY = 1.0D0/Y
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
       FLG11 = FLG(NF)
*
       CLG3A = 0.D0
       IF (CC.EQ.1) THEN
          CLG31 = (144.D0 * DL1**4 - 47024.D0*D27 * DL1**3 + 6319.D0*
     $         DL1**2+ 53160.D0* DL1) * Y1 + DL*DL1 * (72549.D0 +
     $         88238.D0* DL)+ (3709.D0 - 33514.D0* Y - 9533.D0* Y**2) *
     $         Y1+ 66773.D0* Y*DL**2 - 1117.D0* DL + 45.37D0 * DL**2-
     $         5360.D0*D27 * DL**3 - 2044.70D0 *OY*Y1 - 409.506D0 * DL
     $         *OY
          CLG32 = (288.D0*D27 * DL1**3 - 3648.D0*D27 * DL1**2 - 592.3D0
     $         *DL1+ 1511.D0* Y*DL1) * Y1 + DL*DL1 * (311.3D0 + 14.24D0
     $         *DL)+ (577.3D0 - 729.0D0 * Y) * Y1 + 30.78D0 * Y*DL**3
     $         +366.0D0 * DL + 3000.D0*D27 * DL**2 + 480.D0*D27 * DL**3
     $         +88.5037D0 *OY*Y1
       CLG3A = NF * ( CLG31 + NF * CLG32  )
       ELSE
          CLG3F = (-0.0105D0 * DL1**3 + 1.550D0 * DL1**2 + 19.72D0 *Y
     $         *DL1- 66.745D0 * Y + 0.615D0 * Y**2) * Y1 + 20.D0*D27 * Y
     $         *DL**4+ (280.D0*D81 + 2.260D0* Y) * Y*DL**3 - (15.40D0
     $         - 2.201D0* Y)* Y*DL**2 - (71.66D0 - 0.121D0 * Y) * Y*DL
       CLG3A = NF * NF * FLG11 * CLG3F
       ENDIF
*     
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XCLSG3P
