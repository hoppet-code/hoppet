      MODULE XC2NS2E
      CONTAINS
*
* ..File: xc2ns2e.f    F2_NS  (even-N)
*
*
* ..The exact 2-loop MS(bar) non-singlet coefficient functions for the 
*    structure function F_2 in electromagnetic DIS at mu_r = mu_f = Q. 
*    Expansion parameter: alpha_s/(4 pi).
* 
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
*
* =====================================================================
*     
      SUBROUTINE SET_C2SOFT(nf)
      IMPLICIT REAL*8 (A - Z)
      INTEGER NF
      COMMON / C2SOFT / A0, A1, A2, A3
      PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,     Z3 = 1.2020 56903 15959 42854 D0 )

*     
*     ...Colour factors
*     
      CF  = 4./3.D0
      CA  = 3.D0
*     
*     ...The soft (`+'-distribution) part of the coefficient function
*     
      A3 = 
     &     + 8.D0*cf**2
      A2 =
     &     - 22.D0/3.D0*ca*cf
     &     - 18.D0*cf**2
     &     + 4.D0/3.D0*cf*nf
      A1 =
     &     - 8.D0*z2*ca*cf
     &     - 32.D0*z2*cf**2
     &     + 367.D0/9.D0*ca*cf
     &     - 27.D0*cf**2
     &     - 58.D0/9.D0*cf*nf
      A0 =
     &     + 44.D0/3.D0*z2*ca*cf
     &     + 36.D0*z2*cf**2
     &     + 40.D0*z3*ca*cf
     &     - 8.D0*z3*cf**2
     &     - 3155.D0/54.D0*ca*cf
     &     + 51.D0/2.D0*cf**2
     &     + 247.D0/27.D0*cf*nf
     &     - 8.D0/3.D0*z2*cf*nf

      END SUBROUTINE
*     
* ..This is the regular piece. 
*
       FUNCTION X2NP2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4 
       INTEGER NF, NF2, N1, N2, NW, I1, I2, I3, N
       PARAMETER ( N1 = -1, N2 = 1, NW = 3 ) 
       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
*
* ..The soft coefficients for use in X2NP2B and X2NP2C
*
       COMMON / C2SOFT / A0, A1, A2, A3
*
* ...Colour factors
*
       CF  = 4./3.D0
       CA  = 3.D0
*
* ...Some abbreviations
*
       DX = 1.D0/X
       DM = 1.D0/(1.D0-X)
       DP = 1.D0/(1.D0+X)
       DL1 = LOG (1.D0-X)
*
* ...The harmonic polylogs up to weight NW=3 by Gehrmann and Remiddi
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2) 
*
* ...The coefficient function in terms of the harmonic polylogs
*    (without the delta(1-x) part, but with the soft contribution)
*
      c2qq2 =
     &  + nf*cf * (  - 158.D0/27.D0 - 16.D0/3.D0*z2*dm + 8.D0/3.D0*z2*x
     &     + 8.D0/3.D0*z2 + 247.D0/27.D0*dm - 488.D0/27.D0*x - 26.D0/3.D
     &    0*Hr1(0) + 38.D0/3.D0*Hr1(0)*dm - 38.D0/3.D0*Hr1(0)*x - 32.D0/
     &    9.D0*Hr1(1) + 58.D0/9.D0*Hr1(1)*dm - 68.D0/9.D0*Hr1(1)*x - 10.
     &    D0/3.D0*Hr2(0,0) + 20.D0/3.D0*Hr2(0,0)*dm - 10.D0/3.D0*Hr2(0,
     &    0)*x - 8.D0/3.D0*Hr2(0,1) + 16.D0/3.D0*Hr2(0,1)*dm - 8.D0/3.D0
     &    *Hr2(0,1)*x - 4.D0/3.D0*Hr2(1,0) + 8.D0/3.D0*Hr2(1,0)*dm - 4.D
     &    0/3.D0*Hr2(1,0)*x - 4.D0/3.D0*Hr2(1,1) + 8.D0/3.D0*Hr2(1,1)*
     &    dm - 4.D0/3.D0*Hr2(1,1)*x )
      c2qq2 = c2qq2 + cf*ca * ( 3709.D0/135.D0 + 88.D0/3.D0*z2*dm - 104.
     &    D0/3.D0*z2*x + 72.D0/5.D0*z2*x**3 - 44.D0/3.D0*z2 + 4.D0*z3*
     &    dm - 28.D0*z3*dp - 56.D0*z3*x + 12.D0*z3 - 3155.D0/54.D0*dm
     &     - 8.D0/5.D0*dx + 17626.D0/135.D0*x - 72.D0/5.D0*x**2 + 32.D0
     &    *Hr1(-1)*z2*dp + 36.D0*Hr1(-1)*z2*x - 12.D0*Hr1(-1)*z2 + 583.D
     &    0/15.D0*Hr1(0) + 8.D0*Hr1(0)*z2*dm - 8.D0*Hr1(0)*z2*dp - 8.D0
     &    *Hr1(0)*z2*x - 239.D0/3.D0*Hr1(0)*dm + 8.D0*Hr1(0)*dp + 8.D0/
     &    5.D0*Hr1(0)*dx + 1693.D0/15.D0*Hr1(0)*x - 72.D0/5.D0*Hr1(0)*
     &    x**2 + 56.D0/9.D0*Hr1(1) + 24.D0*Hr1(1)*z2*dm - 32.D0*Hr1(1)*
     &    z2*x - 8.D0*Hr1(1)*z2 - 367.D0/9.D0*Hr1(1)*dm + 668.D0/9.D0*
     &    Hr1(1)*x - 36.D0*Hr2(-1,0) - 8.D0/5.D0*Hr2(-1,0)*dx**2 - 20.D0
     &    *Hr2(-1,0)*x + 72.D0/5.D0*Hr2(-1,0)*x**3 + 55.D0/3.D0*Hr2(0,0
     &    ) - 110.D0/3.D0*Hr2(0,0)*dm + 115.D0/3.D0*Hr2(0,0)*x - 72.D0/
     &    5.D0*Hr2(0,0)*x**3 + 44.D0/3.D0*Hr2(0,1) - 88.D0/3.D0*Hr2(0,1
     &    )*dm + 44.D0/3.D0*Hr2(0,1)*x + 22.D0/3.D0*Hr2(1,0) - 44.D0/3.D
     &    0*Hr2(1,0)*dm )
      c2qq2 = c2qq2 + cf*ca * ( 22.D0/3.D0*Hr2(1,0)*x + 22.D0/3.D0*Hr2(
     &    1,1) - 44.D0/3.D0*Hr2(1,1)*dm + 22.D0/3.D0*Hr2(1,1)*x - 8.D0*
     &    Hr3(-1,-1,0) + 32.D0*Hr3(-1,-1,0)*dp + 56.D0*Hr3(-1,-1,0)*x
     &     + 16.D0*Hr3(-1,0,0) - 40.D0*Hr3(-1,0,0)*dp - 40.D0*Hr3(-1,0,
     &    0)*x + 8.D0*Hr3(-1,0,1) - 16.D0*Hr3(-1,0,1)*dp - 8.D0*Hr3(-1,
     &    0,1)*x + 16.D0*Hr3(0,-1,0) - 24.D0*Hr3(0,-1,0)*dm - 24.D0*
     &    Hr3(0,-1,0)*dp - 12.D0*Hr3(0,0,0)*dm + 12.D0*Hr3(0,0,0)*dp + 
     &    12.D0*Hr3(0,0,0)*x - 8.D0*Hr3(0,0,1)*dm + 8.D0*Hr3(0,0,1)*dp
     &     + 8.D0*Hr3(0,0,1)*x + 4.D0*Hr3(1,0,0) - 16.D0*Hr3(1,0,0)*dm
     &     + 28.D0*Hr3(1,0,0)*x + 4.D0*Hr3(1,0,1) - 8.D0*Hr3(1,0,1)*dm
     &     + 4.D0*Hr3(1,0,1)*x - 4.D0*Hr3(1,1,0) + 8.D0*Hr3(1,1,0)*dm
     &     - 4.D0*Hr3(1,1,0)*x )
      c2qq2 = c2qq2 + cf**2 * (  - 124.D0/5.D0 + 24.D0*z2*dm - 8.D0*z2*
     &    x - 144.D0/5.D0*z2*x**3 - 32.D0*z2 + 64.D0*z3*dm + 56.D0*z3*
     &    dp + 72.D0*z3*x - 64.D0*z3 + 51.D0/2.D0*dm + 16.D0/5.D0*dx - 
     &    461.D0/5.D0*x + 144.D0/5.D0*x**2 - 64.D0*Hr1(-1)*z2*dp - 72.D0
     &    *Hr1(-1)*z2*x + 24.D0*Hr1(-1)*z2 - 132.D0/5.D0*Hr1(0) + 48.D0
     &    *Hr1(0)*z2*dm + 16.D0*Hr1(0)*z2*dp - 24.D0*Hr1(0)*z2*x - 40.D0
     &    *Hr1(0)*z2 + 61.D0*Hr1(0)*dm - 16.D0*Hr1(0)*dp - 16.D0/5.D0*
     &    Hr1(0)*dx - 502.D0/5.D0*Hr1(0)*x + 144.D0/5.D0*Hr1(0)*x**2 + 
     &    16.D0*Hr1(1) + 16.D0*Hr1(1)*z2*dm + 32.D0*Hr1(1)*z2*x - 16.D0
     &    *Hr1(1)*z2 + 27.D0*Hr1(1)*dm - 68.D0*Hr1(1)*x + 72.D0*Hr2(-1,
     &    0) + 16.D0/5.D0*Hr2(-1,0)*dx**2 + 40.D0*Hr2(-1,0)*x - 144.D0/
     &    5.D0*Hr2(-1,0)*x**3 + 24.D0*Hr2(0,0) - 6.D0*Hr2(0,0)*dm - 4.D0
     &    *Hr2(0,0)*x + 144.D0/5.D0*Hr2(0,0)*x**3 + 32.D0*Hr2(0,1) - 24.
     &    D0*Hr2(0,1)*dm + 48.D0*Hr2(0,1)*x + 32.D0*Hr2(1,0) - 36.D0*
     &    Hr2(1,0)*dm + 32.D0*Hr2(1,0)*x + 28.D0*Hr2(1,1) - 36.D0*Hr2(1
     &    ,1)*dm )
      c2qq2 = c2qq2 + cf**2 * ( 36.D0*Hr2(1,1)*x + 16.D0*Hr3(-1,-1,0)
     &     - 64.D0*Hr3(-1,-1,0)*dp - 112.D0*Hr3(-1,-1,0)*x - 32.D0*Hr3(
     &    -1,0,0) + 80.D0*Hr3(-1,0,0)*dp + 80.D0*Hr3(-1,0,0)*x - 16.D0*
     &    Hr3(-1,0,1) + 32.D0*Hr3(-1,0,1)*dp + 16.D0*Hr3(-1,0,1)*x - 32.
     &    D0*Hr3(0,-1,0) + 48.D0*Hr3(0,-1,0)*dm + 48.D0*Hr3(0,-1,0)*dp
     &     + 30.D0*Hr3(0,0,0) - 16.D0*Hr3(0,0,0)*dm - 24.D0*Hr3(0,0,0)*
     &    dp + 6.D0*Hr3(0,0,0)*x + 40.D0*Hr3(0,0,1) - 48.D0*Hr3(0,0,1)*
     &    dm - 16.D0*Hr3(0,0,1)*dp + 24.D0*Hr3(0,0,1)*x + 28.D0*Hr3(0,1
     &    ,0) - 48.D0*Hr3(0,1,0)*dm + 28.D0*Hr3(0,1,0)*x + 32.D0*Hr3(0,
     &    1,1) - 56.D0*Hr3(0,1,1)*dm + 32.D0*Hr3(0,1,1)*x + 20.D0*Hr3(1
     &    ,0,0) - 24.D0*Hr3(1,0,0)*dm - 28.D0*Hr3(1,0,0)*x + 24.D0*Hr3(
     &    1,0,1) - 48.D0*Hr3(1,0,1)*dm + 24.D0*Hr3(1,0,1)*x + 32.D0*
     &    Hr3(1,1,0) - 64.D0*Hr3(1,1,0)*dm + 32.D0*Hr3(1,1,0)*x + 24.D0
     &    *Hr3(1,1,1) - 48.D0*Hr3(1,1,1)*dm + 24.D0*Hr3(1,1,1)*x )
*
       C2QQ2L = DM * ( DL1**3 * A3 + DL1**2 * A2 + DL1 * A1 + A0) 
*
* ...The regular piece of the coefficient function
*
       X2NP2A = C2QQ2 - C2QQ2L
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The singular (soft) piece.
*
       FUNCTION X2NS2B (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       COMMON / C2SOFT / A0, A1, A2, A3
*
       DL1 = LOG (1.D0-Y)
       DM  = 1.D0/(1.D0-Y)
*
       X2NS2B = DM * ( DL1**3 * A3 + DL1**2 * A2 + DL1 * A1 + A0)
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The 'local' piece.
*
       FUNCTION X2NS2C (Y, NF)
*
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF, NF2
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
*
       COMMON / C2SOFT / A0, A1, A2, A3
*
* ...Colour factors
*
       CF  = 4./3.D0
       CA  = 3.D0
*
* ...The coefficient of delta(1-x)
*
       C2DELT = 
     &     - 251.D0/3.D0*z2*ca*cf
     &     + 38.D0/3.D0*z2*cf*nf
     &     + 69.D0*z2*cf**2
     &     + 71.D0/5.D0*z2**2*ca*cf
     &     + 6.D0*z2**2*cf**2
     &     + 140.D0/3.D0*z3*ca*cf
     &     + 4.D0/3.D0*z3*cf*nf
     &     - 78.D0*z3*cf**2
     &     - 5465.D0/72.D0*ca*cf
     &     + 457.D0/36.D0*cf*nf
     &     + 331.D0/8.D0*cf**2
*
       DL1 = LOG (1.D0-Y)
*
       X2NS2C =   DL1**4 * A3/4.D0 + DL1**3 * A2/3.D0 
     ,          + DL1**2 * A1/2.D0 + DL1 * A0 + C2DELT
*
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XC2NS2E
