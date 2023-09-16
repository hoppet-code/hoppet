      MODULE XCDIFF2E
      CONTAINS
*     
* ..File: xcdiff2e.f    Fi_NS,  i = 2,3,L  (even-N - odd-N)
*
*
* ..The exact 2-loop MS(bar) even-odd differences of the NS coefficient 
*    functions F2, F3 and FL in charged-current DIS at mu_r = mu_f = Q. 
*    The expansion parameter is normalized as  a_s = alpha_s/(4 pi).
*
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* ..Equivalent expressions for F2 and F3 were first published in
*     E.B. Zijlstra and W.L van Neerven, Phys. Lett. B272 (1991) 127,
*                                        Phys. Lett. B297 (1992) 377.
*   The FL results were also calculated by them, but not published.
*
* ..Reference: M. Rogal, S. Moch and A. Vogt, 
*              arXiv:0807.3731 = Nucl.Phys. B790 (2008) 317-335 
*
* =====================================================================
*
* ..F2
*
       FUNCTION XC2DFF2 (X)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4 
       INTEGER N1, N2, NW, I1, I2, I3, N
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
* ...Colour factors and abbreviations
*
       CF = 4./3.D0
       CA = 3.D0
       CAM2CF = CA - 2.* CF
*
       DX = 1.D0/X
       DP = 1.D0/(1.D0+X)
*
* ...The harmonic polylogs up to weight 3 by Gehrmann and Remiddi
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2) 
*
* ...The coefficient function in terms of the harmonic polylogs
*
       c2qq2d =
     &  + cf*cam2cf * ( 162.D0/5.D0 - 82.D0/5.D0*x - 72.D0/5.D0*x**2
     &     - 8.D0/5.D0*dx + 20.D0*z3 - 68.D0*z3*x - 56.D0*z3*dp - 4.D0*
     &    z2 - 28.D0*z2*x - 48.D0*z2*x**2 + 72.D0/5.D0*z2*x**3 - 24.D0*
     &    Hr1(-1)*z2 + 72.D0*Hr1(-1)*z2*x + 64.D0*Hr1(-1)*z2*dp - 14.D0/
     &    5.D0*Hr1(0) + 146.D0/5.D0*Hr1(0)*x - 72.D0/5.D0*Hr1(0)*x**2
     &     + 8.D0/5.D0*Hr1(0)*dx + 16.D0*Hr1(0)*dp + 8.D0*Hr1(0)*z2 - 8.
     &    D0*Hr1(0)*z2*x - 16.D0*Hr1(0)*z2*dp + 16.D0*Hr1(1) - 16.D0*
     &    Hr1(1)*x - 32.D0*Hr2(-1,0) - 32.D0*Hr2(-1,0)*x - 48.D0*Hr2(-1
     &    ,0)*x**2 + 72.D0/5.D0*Hr2(-1,0)*x**3 - 8.D0/5.D0*Hr2(-1,0)*
     &    dx**2 + 8.D0*Hr2(0,0) + 32.D0*Hr2(0,0)*x + 48.D0*Hr2(0,0)*
     &    x**2 - 72.D0/5.D0*Hr2(0,0)*x**3 + 8.D0*Hr2(0,1) + 8.D0*Hr2(0,
     &    1)*x - 16.D0*Hr3(-1,-1,0) + 112.D0*Hr3(-1,-1,0)*x + 64.D0*
     &    Hr3(-1,-1,0)*dp + 32.D0*Hr3(-1,0,0) - 80.D0*Hr3(-1,0,0)*x -
     &    80.D0*Hr3(-1,0,0)*dp + 16.D0*Hr3(-1,0,1) - 16.D0*Hr3(-1,0,1)*
     &    x - 32.D0*Hr3(-1,0,1)*dp + 16.D0*Hr3(0,-1,0) - 64.D0*Hr3(0,-1
     &    ,0)*x )
       c2qq2d = c2qq2d + cf*cam2cf * (  - 48.D0*Hr3(0,-1,0)*dp - 12.D0
     &    *Hr3(0,0,0) + 12.D0*Hr3(0,0,0)*x + 24.D0*Hr3(0,0,0)*dp - 8.D0
     &    *Hr3(0,0,1) + 8.D0*Hr3(0,0,1)*x + 16.D0*Hr3(0,0,1)*dp )
*
       XC2DFF2 = C2QQ2D 
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
* ..FL
*
       FUNCTION XCLDFF2 (X)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4
       INTEGER N1, N2, NW, I1, I2, I3, N
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
* ...Colour factors and abbreviations
*
       CF = 4./3.D0
       CA = 3.D0
       CAM2CF = CA - 2.* CF
*
       DX = 1.D0/X
       DP = 1.D0/(1.D0+X)
*
* ...The harmonic polylogs up to weight 3 by Gehrmann and Remiddi
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the harmonic polylogs
*
       cLqq2d =
     &  + cf*cam2cf * ( 208.D0/5.D0 - 128.D0/5.D0*x - 48.D0/5.D0*x**2
     &     - 32.D0/5.D0*dx - 32.D0*z3*x - 16.D0*z2*x - 32.D0*z2*x**2 +
     &    48.D0/5.D0*z2*x**3 + 32.D0*Hr1(-1)*z2*x - 16.D0/5.D0*Hr1(0)
     &     + 224.D0/5.D0*Hr1(0)*x - 48.D0/5.D0*Hr1(0)*x**2 + 32.D0/5.D0
     &    *Hr1(0)*dx - 32.D0*Hr2(-1,0) - 32.D0*Hr2(-1,0)*x - 32.D0*Hr2(
     &    -1,0)*x**2 + 48.D0/5.D0*Hr2(-1,0)*x**3 + 16.D0*Hr2(-1,0)*dx
     &     - 32.D0/5.D0*Hr2(-1,0)*dx**2 + 16.D0*Hr2(0,0)*x + 32.D0*Hr2(
     &    0,0)*x**2 - 48.D0/5.D0*Hr2(0,0)*x**3 + 64.D0*Hr3(-1,-1,0)*x
     &     - 32.D0*Hr3(-1,0,0)*x - 32.D0*Hr3(0,-1,0)*x )
*
       XCLDFF2 = CLQQ2D
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
* ..F3
*
       FUNCTION XC3DFF2 (X)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4 
       INTEGER N1, N2, NW, I1, I2, I3, N
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
* ...Colour factors and abbreviations
*
       CF = 4./3.D0
       CA = 3.D0
       CAM2CF = CA - 2.* CF
*
       DX = 1.D0/X
       DP = 1.D0/(1.D0+X)
*
* ...The harmonic polylogs up to weight 3 by Gehrmann and Remiddi
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2) 
*
* ...The coefficient function in terms of the harmonic polylogs
*
       c3qq2d =
     &  + cf*cam2cf * (  - 30.D0 + 30.D0*x + 36.D0*z3 - 20.D0*z3*x -
     &    56.D0*z3*dp - 12.D0*z2 - 4.D0*z2*x - 8.D0*z2*x**2 - 40.D0*
     &    Hr1(-1)*z2 + 24.D0*Hr1(-1)*z2*x + 64.D0*Hr1(-1)*z2*dp - 14.D0
     &    *Hr1(0) - 30.D0*Hr1(0)*x + 16.D0*Hr1(0)*dp + 8.D0*Hr1(0)*z2
     &     - 8.D0*Hr1(0)*z2*x - 16.D0*Hr1(0)*z2*dp + 16.D0*Hr1(1) - 16.D
     &    0*Hr1(1)*x - 8.D0*Hr2(-1,0)*x**2 - 8.D0*Hr2(-1,0)*dx + 16.D0*
     &    Hr2(0,0) + 8.D0*Hr2(0,0)*x + 8.D0*Hr2(0,0)*x**2 + 8.D0*Hr2(0,
     &    1) + 8.D0*Hr2(0,1)*x - 48.D0*Hr3(-1,-1,0) + 16.D0*Hr3(-1,-1,0
     &    )*x + 64.D0*Hr3(-1,-1,0)*dp + 48.D0*Hr3(-1,0,0) - 32.D0*Hr3(
     &    -1,0,0)*x - 80.D0*Hr3(-1,0,0)*dp + 16.D0*Hr3(-1,0,1) - 16.D0*
     &    Hr3(-1,0,1)*x - 32.D0*Hr3(-1,0,1)*dp + 32.D0*Hr3(0,-1,0) - 16.
     &    D0*Hr3(0,-1,0)*x - 48.D0*Hr3(0,-1,0)*dp - 12.D0*Hr3(0,0,0) +
     &    12.D0*Hr3(0,0,0)*x + 24.D0*Hr3(0,0,0)*dp - 8.D0*Hr3(0,0,1) +
     &    8.D0*Hr3(0,0,1)*x + 16.D0*Hr3(0,0,1)*dp )
*
       XC3DFF2 = C3QQ2D 
*
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XCDIFF2E
