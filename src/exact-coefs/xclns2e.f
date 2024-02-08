      MODULE XCLNS2E
      CONTAINS
*     
* ..File: xclns2e.f    FL_NS  (even-N)
*
*
* ..The 1-loop and 2-loop MS(bar) non-singlet coefficient functions 
*    for the structure function F_L in e.m. DIS at mu_r = mu_f = Q. 
*    The expansion parameter is alpha_s/(4 pi).
* 
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* ..The two-loop results were first derived in
*     J. Sanchez Guillen et al, Nucl. Phys. B353 (1991) 337 
*
* =====================================================================
*
*
* ..The one-loop coefficient function 
*
       FUNCTION XLNS1A (X, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       CF = 4./3.D0
       XLNS1A = 4.* CF * X
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..The two-loop coefficient function 
*
       FUNCTION XLNP2A (X, NF)
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
       clqq2 =
     &  + cf*ca * (  - 196.D0/15.D0 + 3458.D0/45.D0*x - 48.D0/5.D0*x**2
     &     - 32.D0/5.D0*dx - 32.D0*z3*x - 16.D0*z2*x + 48.D0/5.D0*z2*
     &    x**3 + 16.D0*Hr1(-1)*z2*x - 16.D0/5.D0*Hr1(0) + 752.D0/15.D0*
     &    Hr1(0)*x - 48.D0/5.D0*Hr1(0)*x**2 + 32.D0/5.D0*Hr1(0)*dx + 92.
     &    D0/3.D0*Hr1(1)*x - 16.D0*Hr1(1)*z2*x - 32.D0*Hr2(-1,0) - 16.D0
     &    *Hr2(-1,0)*x + 48.D0/5.D0*Hr2(-1,0)*x**3 - 32.D0/5.D0*Hr2(-1,
     &    0)*dx**2 + 16.D0*Hr2(0,0)*x - 48.D0/5.D0*Hr2(0,0)*x**3 + 32.D0
     &    *Hr3(-1,-1,0)*x - 16.D0*Hr3(-1,0,0)*x + 16.D0*Hr3(1,0,0)*x )
       clqq2 = cLqq2 + cf**2 * ( 44.D0/5.D0 - 374.D0/5.D0*x + 96.D0/5.D0
     &    *x**2 + 64.D0/5.D0*dx + 64.D0*z3*x + 8.D0*z2*x - 96.D0/5.D0*
     &    z2*x**3 - 32.D0*Hr1(-1)*z2*x - 8.D0/5.D0*Hr1(0) - 208.D0/5.D0
     &    *Hr1(0)*x + 96.D0/5.D0*Hr1(0)*x**2 - 64.D0/5.D0*Hr1(0)*dx - 8.
     &    D0*Hr1(1) - 28.D0*Hr1(1)*x + 32.D0*Hr1(1)*z2*x + 64.D0*Hr2(-1
     &    ,0) + 32.D0*Hr2(-1,0)*x - 96.D0/5.D0*Hr2(-1,0)*x**3 + 64.D0/5.
     &    D0*Hr2(-1,0)*dx**2 - 16.D0*Hr2(0,0)*x + 96.D0/5.D0*Hr2(0,0)*
     &    x**3 + 24.D0*Hr2(0,1)*x + 16.D0*Hr2(1,0)*x + 16.D0*Hr2(1,1)*x
     &     - 64.D0*Hr3(-1,-1,0)*x + 32.D0*Hr3(-1,0,0)*x - 32.D0*Hr3(1,0
     &    ,0)*x )
       clqq2 = cLqq2 + nf*cf * ( 8.D0/3.D0 - 100.D0/9.D0*x - 16.D0/3.D0*
     &    Hr1(0)*x - 8.D0/3.D0*Hr1(1)*x )
* 
       XLNP2A = CLQQ2 
*
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XCLNS2E
