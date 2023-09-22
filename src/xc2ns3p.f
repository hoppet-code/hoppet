      MODULE XC2NS3P
      CONTAINS
*     
* ..File: xc2ns3p.f    F2_NS
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
* ..References: S. Moch, J. Vermaseren and A. Vogt, hep-ph/0209100
*               J. Vermaseren, A. Vogt and S. Moch, hep-ph/0504242
*
*
* =====================================================================
*
*
* ..The regular piece. The rational end-point coefficients are exact, 
*    the rest has been fitted for x between 10^-6 and 1 - 10^-6. 
*
       FUNCTION C2NP3A (Y, DL, NF, CC)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       INTEGER CC ! charged current
       DIMENSION FL(6)
       DATA FL / -1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
*
       FL11 = FL(NF)
*
       Y1 = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
       D27  = 1.D0/27.D0
       D243 = 1.D0/243.D0
*
       C2NP3A = 0.D0
       IF (CC.EQ.1) THEN
          C2NP3A = - 4926.D0 + 7725.D0* Y + 57256.D0* Y**2 + 12898.D0* Y
     $         **3 - 32.D0*D27 * DL**5 - 8796.D0*D243 * DL**4 - 309.1D0
     $         * DL**3 - 899.6D0 * DL**2 - 775.8D0 * DL + 4.719D0 * Y*DL
     $         **5 - 512.D0*D27 * DL1**5 + 6336.D0*D27 * DL1**4 -
     $         3368.D0* DL1**3 - 2978.D0* DL1**2 + 18832.D0* DL1 -
     $         56000.D0* (1.D0-Y)*DL1**2 - DL*DL1 * (6158.D0 + 1836.D0
     $         *DL) + NF * ( 831.6D0 - 6752.D0* Y - 2778.D0* Y**2 +
     $         728.D0* D243 * DL**4 + 12224.D0* D243 * DL**3 + 187.3D0 *
     $         DL**2 + 275.6D0 * DL + 4.102D0 * Y*DL**4 - 1920.D0* D243
     $         * DL1**4 + 153.5D0 * DL1**3 - 828.7D0 * DL1**2 - 501.1D0
     $         * DL1 + 171.0D0 * (1.D0-Y)*DL1**4 + DL*DL1 * (4365.D0 +
     $         716.2D0 * DL - 5983.D0* DL1) ) + NF**2 * ( 129.2D0 * Y +
     $         102.5D0 * Y**2 - 368.D0* D243 * DL**3 - 1984.D0* D243 *
     $         DL**2 - 8.042D0 * DL - 192.D0* D243 * DL1**3 + 18.21D0 *
     $         DL1**2 - 19.09D0 * DL1 + DL*DL1 * ( - 96.07D0 - 12.46D0 *
     $         DL + 85.88D0 * DL1) )
       ELSE
          C2NP3A = FL11*NF * ( ( 126.42D0 - 50.29D0 * Y - 50.15D0 * Y
     $         **2) * Y1- 26.717D0 - 960.D0*D243 * DL**2 * (DL+5.D0) +
     $         59.59D0 * DL- Y*DL**2 * (101.8D0 + 34.79D0 * DL + 3.070D0
     $         * DL**2)- 9.075D0 * Y*Y1*DL1 ) * Y
*
       ENDIF
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The exact singular piece (irrational coefficients truncated)
*
       FUNCTION C2NS3B (Y, DL, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = DL1VAL(Y, DL)
       DM  = DMVAL(Y, DL)
       D81 = 1.D0/81.D0
*
       C2NS3B = + 1536.D0*D81 * DL1**5 - 16320.D0* D81 * DL1**4 +
     $      5.01099D+2 * DL1**3 + 1.17154D+3 * DL1**2 - 7.32845D+3 * DL1
     $      + 4.44276D+3 + NF * ( 640.D0* D81 * DL1**4 - 6592.D0* D81 *
     $      DL1**3 + 220.573D0 * DL1**2 + 294.906D0 * DL1 - 729.359D0 )
     $      + NF**2 * ( 64.D0* D81 * DL1**3 - 464.D0* D81 * DL1**2 +
     $      7.67505D0 * DL1 + 1.00830D0 )
*
       C2NS3B = DM * C2NS3B
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
       FUNCTION C2NP3C (Y, NF, CC)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       INTEGER CC ! charged current
       DIMENSION FL(6)
       DATA FL / -1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
*
       FL11 = FL(NF)

       DL1 = LOG (1.D0-Y)
       D81 = 1.D0/81.D0
       D3  = 1.D0/3.D0
*
       C2NP3C = 0.D0

       IF (CC.EQ.1) THEN
          C2NP3C = + 256.D0*D81 * DL1**6 - 3264.D0*D81 * DL1**5 +
     $         1.252745D+2 * DL1**4 + 3.905133D+2 * DL1**3 - 3.664225D+3
     $         * DL1**2 + 4.44276D+3  * DL1 - 9195.48D0 + 25.10D0 + NF *
     $         ( 128.D0* D81 * DL1**5 - 1648.D0* D81 * DL1**4 +
     $         220.573D0 * D3 * DL1**3 + 147.453D0 * DL1**2 - 729.359D0
     $         * DL1 + 2575.074D0 - 0.387D0 ) + NF**2 * ( 16.D0* D81 *
     $         DL1**4 - 464.D0* D81*D3 * DL1**3 + 7.67505D0 * 5.D-1 *
     $         DL1**2 + 1.0083D0 * DL1 - 103.2521D0 + 0.0155D0 )
       ELSE
       C2NP3C = - FL11*NF * 11.8880D0
       ENDIF
      
*     
       RETURN
       END FUNCTION

* ..For Y values close to 1, use the series expansion close 
*   to Y1=0 instead of full value, for numerical convergence
*     
       FUNCTION Y1VAL (Y, DL)
       IMPLICIT REAL*8 (A - Z)

       IF (ABS(DL).LT.1D-4) THEN
       Y1VAL = - DL - DL**2/2.0D0 - DL**3/6.0D0 - DL**4/24.0D0
     ,         - DL**5/120.0D0
       ELSE
       Y1VAL = 1.0D0 - Y
       ENDIF

       RETURN
       END FUNCTION


* ..For Y values close to 1, use the series expansion close 
*   to Y1=0 instead of full value, for numerical convergence
*     
       FUNCTION DL1VAL (Y, DL)
       IMPLICIT REAL*8 (A - Z)

       IF (ABS(DL).LT.1D-4) THEN
       DL1VAL = LOG(-DL) +  DL/2.0D0 + DL**2/24.0D0 - DL**4/2880.0D0
       ELSE
       DL1VAL = LOG(1.0D0 - Y)
       ENDIF

       RETURN
       END FUNCTION

* ..For Y values close to 1, use the series expansion close 
*   to Y1=0 instead of full value, for numerical convergence
*     
       FUNCTION DMVAL (Y, DL)
       IMPLICIT REAL*8 (A - Z)

       IF (ABS(DL).LT.1D-3) THEN
       DMVAL  = 0.5D0 - 1.0D0/DL - DL/12.0D0 + DL**3/720.0D0
     ,         - DL**5/30240.0D0
       ELSE
       DMVAL = 1.0D0/(1.0D0-Y)
       ENDIF

       RETURN
       END FUNCTION

* =================================================================av==
      END MODULE XC2NS3P
