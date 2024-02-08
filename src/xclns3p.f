      MODULE XCLNS3P
      USE XC2NS3P
      CONTAINS
*     
* ..File: xclns3p.f    FL_NS
*
*
* ..Parametrization of the third-order MS(bar) non-singlet coefficient 
*    functions for the structure function F_L in electromagnetic DIS.
*    mu_r = mu_f = Q.  The expansion parameter is  alpha_s/(4 pi).
* 
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
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
* ..This is the regular piece. The rational end-point coefficients are
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact and requires the dilogarithm Li2(x).
*
       FUNCTION CLNP3A (Y, DL, NF, CC)
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 WGPLG
       INTEGER NF
       INTEGER CC ! charged current
       DIMENSION FL(6)
       DATA FL / -1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
       DOUBLE PRECISION D6, D27, D81
       PARAMETER ( Z2 = 1.6449 34066 84822 64365D0, D6 = 1.0D0/6.0D0,
     $      D27 = 1.0D0/27.0D0, D81 = 1.0D0/81.0D0 )
*     
       DL1 = DL1VAL(Y, DL)
*
       FL11 = FL(NF)
*
       CLNP3A = 0.D0
       IF (CC.EQ.1) THEN
          CLNP3A =  - 2220.5D0 - 7884.D0* Y + 4168.D0* Y**2 - 1280.D0
     $         *D81 *DL**3 - 7456.D0*D27 * DL**2 - 1355.7D0 * DL +
     $         512.D0*D27 * DL1**4 - 177.40D0 * DL1**3 + 650.6D0 *DL1
     $         **2 - 2729.D0* DL1 + 208.3D0 * Y*DL**3 - DL1**3*(1.D0-Y)*
     $         (125.3D0 - 195.6D0 *DL1) - DL*DL1 * (844.7D0 * DL +
     $         517.3D0 * DL1) + NF * ( 408.4D0 - 9.345D0 * Y - 919.3D0 *
     $         Y**2 + 1728.D0*D81 * DL**2 + 200.73D0 * DL - 1792.D0*D81*
     $         Y*DL**3 + 1024.D0*D81 * DL1**3 - 112.35D0 * DL1**2 +
     $         344.1D0 * DL1 + (1.D0-Y)*DL1**2 * (239.7D0 + 20.63D0 *
     $         DL1) + DL*DL1 * (887.3D0 + 294.5D0 * DL - 59.14D0 * DL1)
     $         ) + NF**2 * ( - 19.D0 + (317.D0*D6 - 12.D0*Z2) * Y +
     $         9.D0* Y*DL**2 + DL * (-6.D0 + 50.D0* Y) + 3.D0* Y*DL1**2
     $         + DL1 * (6.D0 - 25.D0* Y) - 6.D0* Y*DL*DL1 + 6.D0* Y*
     $         DDILOG(Y) ) * 64.D0* D81
       ELSE
          CLNP3A = FL11*NF * ( (107.0D0 + 321.05D0 * Y - 54.62D0 * Y**2)
     $         *(1.D0-Y)- 26.717D0 - 320.0D0*D81 * DL**3 - 640.D0*D81 *
     $         DL**2+ 9.773D0 * DL + Y*DL * (363.8D0 + 68.32D0 * DL) ) *
     $         Y
       ENDIF
*     
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece, introduced to fine-tune the accuracy.
*
       FUNCTION CLNP3C (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       CLNP3C = 0.113 + NF * 0.006
*
       RETURN
       END FUNCTION

      DOUBLE PRECISION FUNCTION DDILOG(X)
      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2
      INTEGER I
      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/
      DATA C( 0) / 0.42996 69356 08137 0D0/
      DATA C( 1) / 0.40975 98753 30771 1D0/
      DATA C( 2) /-0.01858 84366 50146 0D0/
      DATA C( 3) / 0.00145 75108 40622 7D0/
      DATA C( 4) /-0.00014 30418 44423 4D0/
      DATA C( 5) / 0.00001 58841 55418 8D0/
      DATA C( 6) /-0.00000 19078 49593 9D0/
      DATA C( 7) / 0.00000 02419 51808 5D0/
      DATA C( 8) /-0.00000 00319 33412 7D0/
      DATA C( 9) / 0.00000 00043 45450 6D0/
      DATA C(10) /-0.00000 00006 05784 8D0/
      DATA C(11) / 0.00000 00000 86121 0D0/
      DATA C(12) /-0.00000 00000 12443 3D0/
      DATA C(13) / 0.00000 00000 01822 6D0/
      DATA C(14) /-0.00000 00000 00270 1D0/
      DATA C(15) / 0.00000 00000 00040 4D0/
      DATA C(16) /-0.00000 00000 00006 1D0/
      DATA C(17) / 0.00000 00000 00000 9D0/
      DATA C(18) /-0.00000 00000 00000 1D0/
      IF(X .EQ. ONE) THEN
         DDILOG=PI6
         RETURN
      ELSE IF(X .EQ. MONE) THEN
         DDILOG=MALF*PI6
         RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
         Y=MONE/(ONE+T)
         S=ONE
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
         Y=MONE-T
         S=MONE
         A=LOG(-T)
         A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
         Y=(MONE-T)/T
         S=ONE
         A=LOG(-T)
         A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
         Y=-T/(ONE+T)
         S=MONE
         A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
         Y=T
         S=ONE
         A=ZERO
      ELSE
         Y=ONE/T
         S=MONE
         A=PI6+HALF*LOG(T)**2
      END IF
      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO I = 18,0,-1
         B0=C(I)+ALFA*B1-B2
         B2=B1
         B1=B0
      END DO
      DDILOG=-(S*(B0-H*B2)+A)
      RETURN
      END

* FD: This defines the dilogarithm function needed in CLNP3A.
*     Taken from cernroutines.
      DOUBLE PRECISION FUNCTION LI2 (X)
C      IMPLICIT REAL* 8 (A - Z)
C*
C      LI2 = DDILOG(Y)
C
C      FUNCTION DDILOG(X)
C
C      END FUNCTION
C
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I

      DIMENSION C(0:19)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)

      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/

      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO I = 19,0,-1
         B0=C(I)+ALFA*B1-B2
         B2=B1
         B1=B0
       END DO
       H=-(S*(B0-H*B2)+A)
      ENDIF

      LI2=H

      RETURN
      END FUNCTION
*
* =================================================================av==
      END MODULE XCLNS3P
