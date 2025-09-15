      MODULE XC3NS3P
      USE XC2NS3P
      CONTAINS
*     
* ..File: xc3ns3p.f    F3^NU+NU(BAR)
*
*
* ..Parametrization of the 3-loop MS(bar) non-singlet coefficient
*    functions for the structure function F_2 in electromagnetic DIS.
*    at  mu_r = mu_f = Q.  The expansion parameter is  alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is one part in thousand or better.
*
* ..Reference: A. Vogt, J. Vermaseren and S. Moch, 
*              hep-ph/0608307 = NP B (Proc. Suppl.) 160 (2006) 44
*
*
* =====================================================================
*
*
* ..The regular piece. The rational end-point coefficients are exact, 
*    the rest has been fitted for x between 10^-6 and 1 - 10^-6. 
*
       FUNCTION C3NM3A (Y, DL, NF, V) 
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       INTEGER V
       DOUBLE PRECISION D27, D81, D243, D405, D729
       PARAMETER (       D27  = 1.D0/27.D0, D81  = 1.D0/81.D0, D243 =
     $      1.D0/243.D0, D405 = 1.0D0/405.0D0, D729 = 1.0D0/729.0D0)
*     
       FL02 = 1.D0

*      if V==0 then return the NS minus coefficient instead of the NS valence
       IF (V.eq.0) FL02 = 0.D0
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
*
       C3NM3A = - 1853.D0 - 5709.D0* Y + Y*Y1* (5600.D0 - 1432.D0* Y) -
     $      536.D0*D405* DL**5 - 4036.D0*D81* DL**4 - 496.95D0 * DL **3
     $      - 1488.D0 * DL**2 - 293.3D0 * DL - 512.D0*D27 * DL1**5 +
     $      8896.D0*D27 * DL1**4 - 1396.D0* DL1**3 + 3990.D0* DL1**2 +
     $      14363.D0* DL1 - 0.463D0 * Y*DL**6 - DL*DL1 * (4007.D0 +
     $      1312.D0*DL) + NF * ( 516.1D0 - 465.2D0 * Y + Y*Y1* (635.3D0
     $      + 310.4D0 * Y) + 304.D0*D81 * DL**4 + 48512.D0*D729 * DL **3
     $      + 305.32D0 * DL**2 + 366.9D0 * DL - 1.200D0 * Y*DL**4 -
     $      640.D0*D81 * DL1**4 + 32576.D0*D243* DL1**3 - 660.7D0 * DL1
     $      **2 + 959.1D0 * DL1 + 31.95D0 * (1.D0-Y)*DL1**4 + DL*DL1 *
     $      (1496.D0 + 270.1D0 * DL - 1191.D0* DL1) ) + NF**2 * (11.32D0
     $      + 51.94D0 * Y - Y*Y1* (44.52D0 + 11.05D0 * Y) - 368.D0* D243
     $      * DL**3 - 2848.D0*D243* DL**2 - 16.00D0 * DL - 64.D0*D81*
     $      DL1**3 + 992.D0*D81* DL1**2 - 49.65D0 * DL1 - DL*DL1 * (
     $      39.99D0 + 5.103D0 * DL - 16.30D0 * DL1) + 0.0647D0 * Y*DL**4
     $      ) + FL02*NF * ( 48.79D0 - (242.4D0 - 150.7D0 * Y ) * Y1 -
     $      16.D0*D27 * DL**5 + 17.26D0* DL**3 - 113.4D0 * DL**2 -
     $      477.0D0 * DL + 2.147D0 * DL1**2 - 24.57D0 * DL1 + Y*DL *
     $      (218.1D0 + 82.27D0 * DL**2) - DL*DL1 * (81.70D0 + 9.412D0 *
     $      DL1) ) * Y1
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The exact singular piece (irrational coefficients truncated)
*
       FUNCTION C3NS3B (Y, DL, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
       DOUBLE PRECISION D81
       PARAMETER (D81 = 1.0D0/81.0D0)
*
       DL1 = DL1VAL(Y, DL)
       DM  = DMVAL(Y, DL)
*
       C3NS3B = + 1536.D0*D81 * DL1**5 - 16320.D0* D81 * DL1**4 +
     $      5.01099D+2 * DL1**3 + 1.17154D+3 * DL1**2 - 7.32845D+3 * DL1
     $      + 4.44276D+3 + NF * ( 640.D0* D81 * DL1**4 - 6592.D0* D81 *
     $      DL1**3 + 220.573D0 * DL1**2 + 294.906D0 * DL1 - 729.359D0 )
     $      + NF**2 * ( 64.D0* D81 * DL1**3 - 464.D0* D81 * DL1**2 +
     $      7.67505D0 * DL1 + 1.00830D0 )
*
       C3NS3B = DM * C3NS3B
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The 'local' piece.  The coefficients of delta(1-x) have been 
*    slightly shifted with respect to their (truncated) exact values.  
*
       FUNCTION C3NM3C (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       DOUBLE PRECISION D3, D81
       PARAMETER ( D81 = 1.D0/81.D0, D3  = 1.D0/3.D0)
*
       DL1 = LOG (1.D0-Y)
*
       C3NM3C = + 256.D0*D81 * DL1**6 - 3264.D0*D81 * DL1**5 + 1.252745D
     $      +2 * DL1**4 + 3.905133D+2 * DL1**3 - 3.664225D+3 * DL1**2 +
     $      4.44276D+3  * DL1 - 9195.48D0 + 22.80D0 + NF * ( 128.D0* D81
     $      * DL1**5 - 1648.D0* D81 * DL1**4 + 220.573D0 * D3 * DL1**3 +
     $      147.453D0 * DL1**2 - 729.359D0 * DL1 + 2575.074D0 + 0.386D0)
     $      + NF**2 * ( 16.D0* D81 * DL1**4 - 464.D0* D81*D3 * DL1**3 +
     $      7.67505D0 * 5.D-1 * DL1**2 + 1.0083D0 * DL1 - 103.2521D0 -
     $      0.0081D0 )
*
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XC3NS3P
