      MODULE XCDIFF3PNEW
      USE XC2NS3P
      CONTAINS

*
* ..File: xcdiff3p.f    F_{i,NS},  i = 2,3,L  (even-N - odd-N)
*
*
* ..Parametrizations of the differences between the even-N and odd-N 
*    based three-loop coefficient functions for the structure functions
*    F_L, F_2 and F_3 in charged-current deep-inelastic scattering.
*    MS(bar) scheme, the expansion parameter is  a_s = alpha_s/(4 pi).
*
* ..Except where the values are very small, the relative accuracy of 
*     these parametrizations, as well as of the convolution results, 
*     is one part in thousand or better.
*
* ..These functions should be used together with the results for the 
*    even-N based F_L, F_2 and the odd-N F_3 respectively presented in 
*              S. Moch, J. Vermaseren and A. Vogt,
*              hep-ph/0411112, hep-ph/0504242 and arXiv:0812.4168
*    For F_2 and F_L the (NC) flavour class fl11 has to be deactivated 
*    there, and for F_3 the fl02 contributions have to be left out.
*
* ..Reference: J. Davies, S. Moch, J. Vermaseren and A. Vogt, 
*              arXiv:1606.08907
*
* ===================================================================== 
*
*
* ..F_2: third-order c_2^{nu+nubar} - c_2^{nu-nubar}
*
       FUNCTION c2q3dfP (Y, DL, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, NF
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
*
       C2Q30 =   273.59D0 - 44.95D0* Y - 73.56D0* Y**2 + 40.68D0* Y**3 +
     $      0.1356D0* DL**5 + 8.483D0* DL**4 + 55.90D0* DL**3 + 120.67D0
     $      * DL**2 + 388.0D0* DL - 329.8D0* DL*DL1 - Y*DL* (316.2D0 +
     $      71.63D0* DL) + 46.30D0*DL1 + 5.447D0* DL1**2
*
       C2Q31 = - 19.093D0 + 12.97D0* Y + 36.44D0* Y**2 - 29.256D0* Y**3
     $      - 0.76D0* DL**4 - 5.317D0* DL**3 - 19.82D0* DL**2 - 38.958D0
     $      * DL - 13.395D0* DL*DL1 + Y*DL* (14.44D0 + 17.74D0*DL) +
     $      1.395D0* DL1
*
       c2q3dfP = (C2Q30 + NF * C2Q31)* Y1
*
       RETURN
       END
*
* ..The `local' piece for F2, artificial but useful for maximal accuracy 
*    of moments and convolutions
*
       FUNCTION c2q3dfPC (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       c2q3dfPC = - 0.0008D0 + 0.0001D0* NF  
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
* ..F_L: third-order c_L^{nu+nubar} - c_L^{nu-nubar}
*
       FUNCTION cLq3dfP (Y, DL, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
*
       CLQ30 = - 620.53D0 - 394.5D0* Y + 1609.D0* Y**2 - 596.2D0* Y**3 +
     $      0.217D0* DL**3 + 62.18D0* DL**2 + 208.47D0* DL - 482.5D0* DL
     $      *DL1 - Y*DL* (1751.D0 - 197.5D0* DL) + 105.5D0* DL1 +
     $      0.442D0* DL1**2
*
       CLQ31 = - 6.500D0 - 12.435D0* Y + 23.66D0* Y**2 + 0.914D0* Y**3 +
     $      0.015D0* DL**3 - 6.627D0* DL**2 - 31.91D0* DL - Y*DL*
     $      (5.711D0 + 28.635D0* DL)
*
       cLq3dfP = (CLQ30 + NF * CLQ31) * Y1**2
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
* ..F_3: third-order c_3^{nu-nubar} - c_3^{nu+nubar}    Note the sign!
*
       FUNCTION c3q3dfP (Y, DL, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
*
       C3Q30 = - 553.5D0 + 1412.5D0* Y - 990.3D0* Y**2 + 361.1D0* Y**3 +
     $      0.1458D0* DL**5 + 9.688D0* DL**4 + 90.62D0* DL**3 + 83.684D0
     $      * DL**2 - 602.32D0* DL - 382.5D0* DL*DL1 - Y*DL* (2.805D0 +
     $      325.92D0* DL) + 133.5D0* DL1 + 10.135D0* DL1**2 
*
       C3Q31 = - 16.777D0 + 77.78D0* Y - 24.81D0* Y**2 - 28.89D0* Y**3 -
     $      0.7714D0* DL**4 - 7.701D0* DL**3 - 21.522D0* DL**2 - 7.897D0
     $      * DL - 16.17D0* DL*DL1 + Y*DL* (43.21D0 + 67.04D0*DL) +
     $      1.519D0* DL1
*
       c3q3dfP = (C3Q30 + NF * C3Q31)* Y1
*
       RETURN
       END
*
* ..The `local' piece for F3 - present for the same reason as for F2
*
       FUNCTION c3q3dfPC (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       c3q3dfPC = - 0.0029D0 + 0.00006D0* NF
*
       RETURN
       END
*
* =================================================================av==
       END MODULE XCDIFF3PNEW
