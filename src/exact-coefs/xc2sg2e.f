      MODULE XC2SG2E
      CONTAINS
*     
* ..File: xc2sg2e.f    F2_G,PS
*
*
* ..The one-loop and two-loop MS(bar) singlet coefficient functions 
*    for the structure function F_2 in e.m. DIS at mu_r = mu_f = Q. 
*    The expansion parameter is alpha_s/(4 pi).
* 
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* =====================================================================
*
*
* ..The one-loop gluonic coefficient function 
*
       FUNCTION X2G1A (X, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       lx  = log (x)
       l1x = log (1.d0 - x)
*
       X2G1A = 
     &   + nf * ( -2.D0 + 16.D0*x - 16.D0*x**2 + 2.D0*l1x - 4.D0*l1x*x
     &   + 4.D0*l1x*x**2 - 2.D0*lx + 4.D0*lnx*x - 4.D0*lx*x**2 )
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The two-loop pure-singlet coefficient function 
*
       FUNCTION X2S2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4, HC5
       INTEGER NF, NF2, N1, N2, NW
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
       CF  = 4./3.D0
       CA  = 3.D0
*
       DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=3 by Gehrmann and Remiddi 
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,             HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
      c2qps2 =
     &  + nf*cf * ( 158.D0/9.D0 - 422.D0/9.D0*x + 448.D0/27.D0*x**2 +
     &    344.D0/27.D0*dx - 8.D0*z3 - 8.D0*z3*x - 16.D0*z2*x + 16.D0*z2
     &    *x**2 - 16.D0/3.D0*z2*dx + 56.D0*Hr1(0) - 88.D0/3.D0*Hr1(0)*x
     &     - 128.D0/9.D0*Hr1(0)*x**2 - 16.D0*Hr1(0)*z2 - 16.D0*Hr1(0)*
     &    z2*x + 104.D0/3.D0*Hr1(1) - 80.D0/3.D0*Hr1(1)*x + 32.D0/9.D0*
     &    Hr1(1)*x**2 - 104.D0/9.D0*Hr1(1)*dx - 16.D0*Hr2(-1,0) - 16.D0
     &    *Hr2(-1,0)*x - 16.D0/3.D0*Hr2(-1,0)*x**2 - 16.D0/3.D0*Hr2(-1,
     &    0)*dx - 2.D0*Hr2(0,0) + 30.D0*Hr2(0,0)*x - 64.D0/3.D0*Hr2(0,0
     &    )*x**2 - 16.D0*Hr2(0,1)*x**2 + 4.D0*Hr2(1,0) - 4.D0*Hr2(1,0)*
     &    x - 16.D0/3.D0*Hr2(1,0)*x**2 + 16.D0/3.D0*Hr2(1,0)*dx + 4.D0*
     &    Hr2(1,1) - 4.D0*Hr2(1,1)*x - 16.D0/3.D0*Hr2(1,1)*x**2 + 16.D0/
     &    3.D0*Hr2(1,1)*dx + 20.D0*Hr3(0,0,0) + 20.D0*Hr3(0,0,0)*x + 16.
     &    D0*Hr3(0,0,1) + 16.D0*Hr3(0,0,1)*x + 8.D0*Hr3(0,1,0) + 8.D0*
     &    Hr3(0,1,0)*x + 8.D0*Hr3(0,1,1) + 8.D0*Hr3(0,1,1)*x )
* 
       X2S2A = C2QPS2 
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The two-loop gluonic coefficient function 
*
       FUNCTION X2G2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4, HC5
       INTEGER NF, NF2, N1, N2, NW
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
       CF  = 4./3.D0
       CA  = 3.D0
*
       DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=3 by Gehrmann and Remiddi 
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,             HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
      c2gg2 =
     &  + nf*ca * ( 239.D0/9.D0 + 1072.D0/9.D0*x - 4493.D0/27.D0*x**2
     &     + 344.D0/27.D0*dx + 4.D0*z3 - 48.D0*z3*x + 24.D0*z3*x**2 + 8.
     &    D0*z2 - 144.D0*z2*x + 148.D0*z2*x**2 - 16.D0/3.D0*z2*dx - 4.D0
     &    *Hr1(-1)*z2 - 8.D0*Hr1(-1)*z2*x - 16.D0*Hr1(-1)*z2*x**2 + 58.D
     &    0*Hr1(0) + 584.D0/3.D0*Hr1(0)*x - 2090.D0/9.D0*Hr1(0)*x**2 -
     &    8.D0*Hr1(0)*z2 - 64.D0*Hr1(0)*z2*x + 16.D0*Hr1(0)*z2*x**2 +
     &    62.D0/3.D0*Hr1(1) + 454.D0/3.D0*Hr1(1)*x - 1570.D0/9.D0*Hr1(1
     &    )*x**2 - 104.D0/9.D0*Hr1(1)*dx + 8.D0*Hr1(1)*z2 - 16.D0*Hr1(1
     &    )*z2*x + 8.D0*Hr1(1)*z2*x**2 - 24.D0*Hr2(-1,0) + 80.D0/3.D0*
     &    Hr2(-1,0)*x**2 - 16.D0/3.D0*Hr2(-1,0)*dx - 2.D0*Hr2(0,0) +
     &    176.D0*Hr2(0,0)*x - 388.D0/3.D0*Hr2(0,0)*x**2 - 8.D0*Hr2(0,1)
     &     + 144.D0*Hr2(0,1)*x - 148.D0*Hr2(0,1)*x**2 - 4.D0*Hr2(1,0)
     &     + 80.D0*Hr2(1,0)*x - 268.D0/3.D0*Hr2(1,0)*x**2 + 16.D0/3.D0*
     &    Hr2(1,0)*dx - 4.D0*Hr2(1,1) + 72.D0*Hr2(1,1)*x - 244.D0/3.D0*
     &    Hr2(1,1)*x**2 + 16.D0/3.D0*Hr2(1,1)*dx + 8.D0*Hr3(-1,-1,0) +
     &    16.D0*Hr3(-1,-1,0)*x )
      c2gg2 = c2gg2 + nf*ca * ( 8.D0*Hr3(-1,0,0) + 16.D0*Hr3(-1,0,0)*x
     &     + 24.D0*Hr3(-1,0,0)*x**2 + 8.D0*Hr3(-1,0,1) + 16.D0*Hr3(-1,0
     &    ,1)*x + 16.D0*Hr3(-1,0,1)*x**2 + 16.D0*Hr3(0,-1,0)*x**2 + 20.D
     &    0*Hr3(0,0,0) + 56.D0*Hr3(0,0,0)*x + 8.D0*Hr3(0,0,1) + 64.D0*
     &    Hr3(0,0,1)*x - 16.D0*Hr3(0,0,1)*x**2 + 48.D0*Hr3(0,1,0)*x -
     &    16.D0*Hr3(0,1,0)*x**2 + 48.D0*Hr3(0,1,1)*x - 16.D0*Hr3(0,1,1)
     &    *x**2 - 12.D0*Hr3(1,0,0) + 24.D0*Hr3(1,0,0)*x - 16.D0*Hr3(1,0
     &    ,0)*x**2 - 4.D0*Hr3(1,0,1) + 8.D0*Hr3(1,0,1)*x - 8.D0*Hr3(1,0
     &    ,1)*x**2 - 12.D0*Hr3(1,1,0) + 24.D0*Hr3(1,1,0)*x - 24.D0*Hr3(
     &    1,1,0)*x**2 - 4.D0*Hr3(1,1,1) + 8.D0*Hr3(1,1,1)*x - 8.D0*Hr3(
     &    1,1,1)*x**2 )
      c2gg2 = c2gg2 + nf*cf * (  - 647.D0/15.D0 + 239.D0/5.D0*x - 36.D0/
     &    5.D0*x**2 + 8.D0/15.D0*dx + 32.D0*z3 + 72.D0*z3*x**2 + 16.D0*
     &    z2 - 104.D0/3.D0*z2*x + 72.D0*z2*x**2 + 96.D0/5.D0*z2*x**3 -
     &    16.D0*Hr1(-1)*z2 - 32.D0*Hr1(-1)*z2*x - 16.D0*Hr1(-1)*z2*x**2
     &     - 236.D0/15.D0*Hr1(0) + 113.D0/5.D0*Hr1(0)*x - 216.D0/5.D0*
     &    Hr1(0)*x**2 - 8.D0/15.D0*Hr1(0)*dx + 16.D0*Hr1(0)*z2 - 32.D0*
     &    Hr1(0)*z2*x + 48.D0*Hr1(0)*z2*x**2 - 14.D0*Hr1(1) + 40.D0*
     &    Hr1(1)*x - 24.D0*Hr1(1)*x**2 + 8.D0*Hr1(1)*z2 - 16.D0*Hr1(1)*
     &    z2*x + 32.D0*Hr1(1)*z2*x**2 + 48.D0*Hr2(-1,0) + 64.D0/3.D0*
     &    Hr2(-1,0)*x + 96.D0/5.D0*Hr2(-1,0)*x**3 + 8.D0/15.D0*Hr2(-1,0
     &    )*dx**2 - 3.D0*Hr2(0,0) + 44.D0/3.D0*Hr2(0,0)*x - 72.D0*Hr2(0
     &    ,0)*x**2 - 96.D0/5.D0*Hr2(0,0)*x**3 - 16.D0*Hr2(0,1) + 56.D0*
     &    Hr2(0,1)*x - 72.D0*Hr2(0,1)*x**2 - 26.D0*Hr2(1,0) + 80.D0*
     &    Hr2(1,0)*x - 72.D0*Hr2(1,0)*x**2 - 26.D0*Hr2(1,1) + 80.D0*
     &    Hr2(1,1)*x - 72.D0*Hr2(1,1)*x**2 - 32.D0*Hr3(-1,-1,0) - 64.D0
     &    *Hr3(-1,-1,0)*x )
      c2gg2 = c2gg2 + nf*cf * (  - 32.D0*Hr3(-1,-1,0)*x**2 + 16.D0*Hr3(
     &    -1,0,0) + 32.D0*Hr3(-1,0,0)*x + 16.D0*Hr3(-1,0,0)*x**2 + 32.D0
     &    *Hr3(0,-1,0) + 32.D0*Hr3(0,-1,0)*x**2 - 10.D0*Hr3(0,0,0) + 20.
     &    D0*Hr3(0,0,0)*x - 40.D0*Hr3(0,0,0)*x**2 - 16.D0*Hr3(0,0,1) +
     &    32.D0*Hr3(0,0,1)*x - 48.D0*Hr3(0,0,1)*x**2 - 12.D0*Hr3(0,1,0)
     &     + 24.D0*Hr3(0,1,0)*x - 32.D0*Hr3(0,1,0)*x**2 - 16.D0*Hr3(0,1
     &    ,1) + 32.D0*Hr3(0,1,1)*x - 40.D0*Hr3(0,1,1)*x**2 - 4.D0*Hr3(1
     &    ,0,0) + 8.D0*Hr3(1,0,0)*x - 24.D0*Hr3(1,0,0)*x**2 - 24.D0*
     &    Hr3(1,0,1) + 48.D0*Hr3(1,0,1)*x - 48.D0*Hr3(1,0,1)*x**2 - 16.D
     &    0*Hr3(1,1,0) + 32.D0*Hr3(1,1,0)*x - 32.D0*Hr3(1,1,0)*x**2 -
     &    20.D0*Hr3(1,1,1) + 40.D0*Hr3(1,1,1)*x - 40.D0*Hr3(1,1,1)*x**2
     &     )
* 
       X2G2A = C2GG2 
*
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XC2SG2E
