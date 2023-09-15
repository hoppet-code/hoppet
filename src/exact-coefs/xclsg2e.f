      MODULE XCLSG2E
      CONTAINS
*     
* ..File: xclg2e.f    FL_G,PS
*
*
* ..The one-loop and two-loop MS(bar) singlet coefficient functions 
*    for the structure function F_L in e.m. DIS at mu_r = mu_f = Q. 
*    The expansion parameter is alpha_s/(4 pi).
* 
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* ..The two-loop results were first derived in
*     J. Sanchez Guillen et al, Nucl. Phys. B353 (1991) 337  and
*     E.B. Zijlstra and W.L. van Neerven, Phys. Lett. B273 (1991) 476
*
* =====================================================================
*
*
* ..The one-loop gluonic coefficient function 
*
       FUNCTION XLG1A (X, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       XLG1A = 8.* NF * X*(1.-X)
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The two-loop pure-singlet coefficient function 
*
       FUNCTION XLS2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4, HC5
       INTEGER NF, N1, N2, NW
       PARAMETER ( N1 = -1, N2 = 1, NW = 2 )
       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2),
     ,           HC5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2)
       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2),
     ,           HR5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2)
       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2),
     ,           HI5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2)
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0 )
*
* ...Colour factors and abbreviations
*
       CF  = 4./3.D0
       CA  = 3.D0
*
       DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=2 by Gehrmann and Remiddi 
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,             HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
       clqps2 =
     &  + nf*cf * ( 16.D0/3.D0 - 64.D0/3.D0*x + 160.D0/9.D0*x**2 - 16.D0
     &    /9.D0*dx - 16.D0*z2*x + 16.D0*Hr1(0) - 16.D0*Hr1(0)*x - 32.D0
     &    *Hr1(0)*x**2 + 16.D0*Hr1(1) - 32.D0/3.D0*Hr1(1)*x**2 - 16.D0/
     &    3.D0*Hr1(1)*dx + 32.D0*Hr2(0,0)*x + 16.D0*Hr2(0,1)*x )
* 
       XLS2A = CLQPS2 
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The two-loop gluonic coefficient function 
*
       FUNCTION XLG2A (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4, HC5
       INTEGER NF, N1, N2, NW
       PARAMETER ( N1 = -1, N2 = 1, NW = 2 )
       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2),
     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2),
     ,           HC5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2)
       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2),
     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2),
     ,           HR5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2)
       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2),
     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2),
     ,           HI5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2)
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0 )
*
* ...Colour factors and abbreviations
*
       CF  = 4./3.D0
       CA  = 3.D0
*
       DX = 1.D0/X
*
* ...Harmonic polylogs (HPLs) up to weight NW=2 by Gehrmann and Remiddi 
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,             HI1,HI2,HI3,HI4, N1, N2)
*
* ...The coefficient function in terms of the HPLs
*
       clgg2 =
     &  + nf*ca * ( 16.D0/3.D0 + 272.D0/3.D0*x - 848.D0/9.D0*x**2 - 16.D
     &    0/9.D0*dx - 64.D0*z2*x + 32.D0*z2*x**2 + 16.D0*Hr1(0) + 128.D0
     &    *Hr1(0)*x - 208.D0*Hr1(0)*x**2 + 16.D0*Hr1(1) + 144.D0*Hr1(1)
     &    *x - 464.D0/3.D0*Hr1(1)*x**2 - 16.D0/3.D0*Hr1(1)*dx + 32.D0*
     &    Hr2(-1,0)*x + 32.D0*Hr2(-1,0)*x**2 + 96.D0*Hr2(0,0)*x + 96.D0
     &    *Hr2(0,1)*x - 32.D0*Hr2(0,1)*x**2 + 32.D0*Hr2(1,0)*x - 32.D0*
     &    Hr2(1,0)*x**2 + 32.D0*Hr2(1,1)*x - 32.D0*Hr2(1,1)*x**2 )
       clgg2 = CLgg2 + nf*cf * ( - 128.D0/15.D0 - 304.D0/5.D0*x + 336.D0
     &    /5.D0*x**2 + 32.D0/15.D0*dx + 16.D0/3.D0*z2*x + 64.D0/5.D0*z2
     &    *x**3 - 104.D0/15.D0*Hr1(0) - 208.D0/5.D0*Hr1(0)*x + 96.D0/5.D
     &    0*Hr1(0)*x**2 - 32.D0/15.D0*Hr1(0)*dx - 8.D0*Hr1(1) - 24.D0*
     &    Hr1(1)*x + 32.D0*Hr1(1)*x**2 - 32.D0/3.D0*Hr2(-1,0)*x + 64.D0/
     &    5.D0*Hr2(-1,0)*x**3 + 32.D0/15.D0*Hr2(-1,0)*dx**2 - 64.D0/3.D0
     &    *Hr2(0,0)*x - 64.D0/5.D0*Hr2(0,0)*x**3 - 16.D0*Hr2(0,1)*x )
* 
       XLG2A = CLGG2 
*
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XCLSG2E
