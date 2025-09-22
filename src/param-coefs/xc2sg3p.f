      MODULE XC2SG3P
      USE XC2NS3P
      CONTAINS
*     
* ..File: xc2sg3p.f    F_2^PS  and  F_2^G
*
*
* ..Parametrizations of the 3-loop MS(bar) pure-singlet and gluon coef-
*    ficient functions for the electromagnetic structure function F_2
*    at  mu_r = mu_f = Q. The expansion parameter is  alpha_s/(4 pi).
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is one part in thousand or better.
*
* ..Reference: J. Vermaseren, A. Vogt and S. Moch 
*              hep-ph/0504242 = Nucl. Phys. B724 (2005) 3
* 
* =====================================================================
*
*
* ..The pure-singlet coefficient function, regular piece
*
       FUNCTION C2S3A (Y, DL, NF, CC)
       IMPLICIT REAL*8 (A-Z)
       DIMENSION FL(6), FLS(6)
       INTEGER NF
       INTEGER CC ! charged current
       DATA FL  / -1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
       DATA FLS /  1.d0, 0.1d0, 0.d0, 0.1d0, 0.018181818d0, 0.1d0 /
       DOUBLE PRECISION D9, D81, D243
       PARAMETER (D9=1.D0/9.D0, D81=1.0D0/81.0D0, D243=1.0D0/243.0D0)
*     
       OY = 1.0D0/Y
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
       
       FL11 = FL(NF)
       FLS11 = FLS(NF)
*
       C2S3A = 0.D0
       
       IF (CC.EQ.1) THEN
          C2S31 = ( 856.D0*D81 * DL1**4 - 6032.D0*D81 * DL1**3 +
     $         130.57D0* DL1**2- 542.0D0 * DL1 + 8501.D0 - 4714.D0* Y +
     $         61.50D0 * Y**2 ) * Y1+ DL*DL1 * (8831.D0* DL + 4162.D0*
     $         Y1) - 15.44D0 * Y*DL**5+ 3333.D0* Y*DL**2 + 1615.D0* DL +
     $         1208.D0* DL**2- 333.73D0 * DL**3 + 4244.D0*D81 * DL**4 -
     $         40.D0*D9 * DL**5- 2731.82D0 * Y1*OY - 414.262D0 * DL*OY
          C2S32 = ( - 64.D0*D81 * DL1**3 + 208.D0*D81 * DL1**2 + 23.09D0
     $         * DL1- 220.27D0 + 59.80D0 * Y - 177.6D0 * Y**2) * Y1 -
     $         DL*DL1 * (160.3D0 * DL + 135.4D0 * Y1) - 24.14D0 * Y*DL
     $         **3- 215.4D0 * Y*DL**2 - 209.8D0 * DL - 90.38D0 * DL**2-
     $         3568.D0*D243* DL**3 - 184.D0*D81 * DL**4 + 40.2426D0 *
     $         Y1*OY
       C2S3A = NF * ( C2S31 +  NF * C2S32 )
       ELSE
          C2S3F = ( ( 126.42D0 - 50.29D0 * Y - 50.15D0 * Y**2) * Y1 -
     $         26.717D0- 320.D0*D81 * DL**2 * (DL+5.D0) + 59.59D0 * DL-
     $         Y*DL**2 * (101.8D0 + 34.79D0 * DL + 3.070D0 * DL**2)-
     $         9.075D0 * Y*Y1*DL1 ) * Y
       C2S3A = NF * (FLS11-FL11) * C2S3F
       ENDIF
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The (truncated) 'local' piece due to the FL11 contribution
*
       FUNCTION C2S3C (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       DIMENSION FL(6), FLS(6)
       DATA FL  / -1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
       DATA FLS /  1.d0, 0.1d0, 0.d0, 0.1d0, 0.018181818d0, 0.1d0 /
*
       FL11 = FL(NF)
       FLS11 = FLS(NF)
       C2S3C = - (FLS11-FL11) * NF * 11.8880D0
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The gluon coefficient function
*
       FUNCTION C2G3A (Y, DL, NF, CC)
       IMPLICIT REAL*8 (A-Z)
       DIMENSION FLG(6)
       INTEGER NF
       INTEGER CC ! charged current
       DATA FLG / 1.d0, 0.1d0, 0.d0, 0.1d0, 0.018181818d0, 0.1d0 /
       DOUBLE PRECISION D9, D81
       PARAMETER (D9=1.D0/9.D0, D81=1.0D0/81.0D0)
*
       YI  = 1.D0/Y
       DL1 = DL1VAL(Y, DL)
       FLG11 = FLG(NF)
*
       C2G3A = 0.D0
       IF (CC.EQ.1) THEN
          C2G31 = 966.D0*D81 * DL1**5 - 935.5D0*D9 * DL1**4 + 89.31D0 *
     $         DL1**3 + 979.2D0 * DL1**2 - 2405.D0 * DL1 + 1372.D0*
     $         (1.D0-Y)* DL1**4 - 15729.D0 - 310510.D0* Y + 331570.D0* Y
     $         **2 - 244150.D0* Y*DL**2 - 253.3D0* Y*DL**5 + DL*DL1 *
     $         (138230.D0 - 237010.D0* DL) - 11860.D0* DL - 700.8D0 * DL
     $         **2 - 1440.D0* DL**3 + 2480.5D0*D81 * DL**4 - 134.D0*D9 *
     $         DL**5 - 6362.54D0 * YI - 932.089D0 * DL*YI
          C2G32 = 131.D0*D81 * DL1**4 - 14.72D0 * DL1**3 + 3.607D0 * DL1
     $         **2 - 226.1D0 * DL1 + 4.762D0 - 190.0D0 * Y - 818.4D0 * Y
     $         **2 - 4019.D0* Y*DL**2 - DL*DL1 * (791.5D0 + 4646D0 * DL)
     $         + 739.0D0 * DL + 418.0D0 * DL**2 + 104.3D0 * DL**3 +
     $         809.D0*D81 * DL**4 + 12.D0*D9 * DL**5 + 84.423D0 * YI
       C2G3A = NF * ( C2G31 + NF * C2G32 )
       ELSE
          C2G3F =   3.211D0 * DL1**2 + 19.04D0 * Y*DL1 + 0.623D0 * (1.D0
     $         -Y)*DL1**3- 64.47D0 * Y + 121.6D0 * Y**2 - 45.82D0 * Y**3
     $         - Y*DL*DL1* ( 31.68D0 + 37.24D0 * DL) - Y*DL * (82.40D0 +
     $         16.08D0 * DL)+ Y*DL**3 * (520.D0*D81 + 11.27D0 * Y) +
     $         60.D0*D81 * Y*DL**4
       C2G3A = NF * NF * FLG11 * C2G3F
       ENDIF
*     
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The artificial 'local' piece, introduced to fine-tune the accuracy. 
*
       FUNCTION C2G3C (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       C2G3C = 0.625D0 * NF
*
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XC2SG3P
