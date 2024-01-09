module xpns3m
  character(len=*), parameter :: name_xpns3m = "xpns3m"
contains
!
! ..File: xpns3m.f   
!
! ..The 4-loop MSbar splitting function P_ns^(3)- for the evolution  
!    of non-singlet_- combinations of quark and anti-quark densities.
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!
! ..The nf^{0,1} leading large-nc contributions and the nf^2 part are
!    high-accuracy (0.1% or better) parametrizations of the exact 
!    results. The nf^3 expression is exact up to numerical truncations.
!
! ..The remaining nf^{0,1} terms are approximations based on the first
!    eight odd moments together with small-x and large-x constraints.
!    The two sets spanning the error estimate are called via  IMOD = 1
!    and  IMOD = 2.  Any other value of IMOD invokes their average.
!
! ..The distributions (in the mathematical sense) are given as in eq.
!    (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
!    The name-endings A, B, and C of the functions below correspond to 
!    the kernel superscripts [2], [3], and [1] in that equation.
!
!  ..Reference: 
!        S. Moch, B. Ruijl, T. Ueda, J.Vermaseren and A. Vogt,
!         DESY 17-106, Nikhef 2017-034, LTH 1139
!
!  ..Additional reference for the nf^2 and nf^3 parts
!        J. Davies, B. Ruijl, T. Ueda, J.Vermaseren and A. Vogt
!         arXiv:1610.07477, Nucl. Phys. B915 (2017) 335
!
! =====================================================================
!
!
! ..The regular piece of P_ns^(3)-. 
!
       FUNCTION P3NSMA (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf
!
       Y1  = 1.D0-Y
       DM  = 1.D0/Y1
       DL  = LOG (Y)
       DL1 = LOG (1.D0-Y)
!
! Leading large-n_c, nf^0 and nf^1, parametrized
!
       P3NSA0  = 2.5D+4* ( Y1* ( 3.5254D0 + 8.6935D0* Y - 1.5051D0* Y**2 &
               + 1.8300D0* Y**3 ) + 11.883D0* Y*DL - 0.09066D0* Y*DL**2&
               + 11.410D0* Y1*DL1 + 13.376D0 * DL*DL1 )&
               + 5.167133D+4*DL + 1.712095D+4*DL**2 + 2.863226D+3*DL**3&
               + 2.978255D+2*DL**4 + 1.6D+1*DL**5 + 5.D-1*DL**6&
               - 2.973385D+4 + 1.906980D+4*DL1 
!
       P3NSA1  = 2.5D+4* ( Y1* ( - 0.74077D0 + 1.4860D0* Y - 0.23631D0* Y**2&
               + 0.31584D0* Y**3 ) + 2.5251D0* Y1*DL1 + 2.5203D0* DL*DL1&
               + 2.2242D0* Y*DL - 0.02460D0* Y*DL**2 + 0.00310D0* Y*DL**3 )&
               - 9.239374D+3*DL - 2.917312D+3*DL**2 &
               - 4.305308D+2*DL**3 - 3.6D+1*DL**4 - 4.D0/3.D+0*DL**5&
               + 8.115605D+3 - 3.079761D+3*DL1
!
! Nonleading large-n_c, nf^0 and nf^1: two approximations
!
       P3NMA01 =  (5992.88D0* (1.D0+2.*Y) + 31321.44D0* Y*Y)*Y1 + 511.228D0&
                 - 1618.07D0* DL + 2.25480D0* DL**3 + 31897.82D0* DL1*Y1 &
                 + 4653.76D0* DL1**2*Y1 + 4.964335D-1* (DL**6 + 6.D0*DL**5)&
                 - 2.601749D+3 - 2.118867D+3*DL1
       P3NMA02 = ( 4043.59D0 - 15386.6D0* Y)* Y*Y1 + 502.481D0&
                 + 1532.96D0 * DL**2 + 31.6023D0* DL**3 - 3997.39D0 * DL1*Y1&
                 + 511.567D0* DL1**3*Y1 + 4.964335D-1* (DL**6 + 18.D0*DL**5)&
                 - 2.601749D+3 - 2.118867D+3*DL1
!
       P3NMA11 =  (114.457D0* (1.D0+2.D0*Y) + 2570.73D0* Y*Y)*Y1 - 7.08645D0&
                 - 127.012D0* DL**2 + 2.69618D0* DL**4 + 1856.63D0* DL1*Y1 &
                 + 440.17D0* DL1**2*Y1 + 3.121643D+2 + 3.379310D+2*DL1
       P3NMA12 = (-335.995D0* (2.D0+Y) -1605.91D0* Y*Y)*Y1 - 7.82077D0&
                 - 9.76627D0* DL**2 + 0.14218D0* DL**5 - 1360.04D0* DL1*Y1 &
                 + 38.7337D0* DL1**3*Y1 + 3.121643D+2 + 3.379310D+2*DL1
!
! nf^2 (parametrized) and nf^3 (exact)
!
       P3NSMA2 = 2.5D+2*  ( Y1* ( 3.2206D0 + 1.7507D0* Y + 0.13281D0* Y**2 &
               + 0.45969D0* Y**3 ) + 1.5641D0* Y*DL - 0.37902D0* Y*DL**2&
               - 0.03248D0* Y*DL**3 + 2.7511D0* Y1*DL1 + 3.2709D0 * DL*DL1 )&
               + 4.378810D+2*DL + 1.282948D+2*DL**2 + 1.959945D+1*DL**3&
               + 9.876543D-1*DL**4 - 3.760092D+2 + 2.668861D+1*DL1
!
       P3NSA3  = - 2.426296D+0 - 8.460488D-1* Y &
               + ( 5.267490D-1* DM - 3.687243D+0 + 3.160494D+0* Y )* DL&
               - ( 1.316872D+0* (DM+1.D-1) - 1.448560D+0*Y )* DL**2&
               - ( 2.633744D-1*DM - 1.31687D-1* (1.D0+Y) )* DL**3
!
! Assembly
!
       P3NSMAI = P3NSA0 + nf*P3NSA1 + nf**2*P3NSMA2 + nf**3*P3NSA3
       IF (IMOD .EQ. 1) THEN
         P3NSMA = P3NSMAI + P3NMA01 + nf* P3NMA11
       ELSE IF (IMOD .EQ. 2) THEN
         P3NSMA = P3NSMAI + P3NMA02 + nf* P3NMA12
       ELSE
         P3NSMA = P3NSMAI &
                + 0.5D0* ((P3NMA01+P3NMA02) + nf* (P3NMA11+P3NMA12))
       END IF
!*
       RETURN
     END FUNCTION P3NSMA
!
! ---------------------------------------------------------------------
!
! ..The singular piece.
!
       FUNCTION P3NSMB (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf
       D1 = 1.D0/(1.D0-Y)
!
       A4qI  = 2.120902D+4 &
             - 5.179372D+3* nf &
             + 1.955772D+2* nf**2&
             + 3.272344D+0* nf**3 
       A4ap1 = -511.228D0 + 7.08645D0*nf
       A4ap2 = -502.481D0 + 7.82077D0*nf
!
       IF (IMOD .EQ. 1) THEN
         P3NSMB = (A4qI + A4ap1)* D1
       ELSE IF (IMOD .EQ. 2) THEN
         P3NSMB = (A4qI + A4ap2)* D1
       ELSE
         P3NSMB = (A4qI + 0.5D0* (A4ap1+A4ap2) )* D1
       ENDIF
!
       RETURN
     END FUNCTION P3NSMB
!
! ---------------------------------------------------------------------
!
! ..The 'local' piece.
!
       FUNCTION P3NSMC (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf
       DL1 = LOG (1.D0-Y)
!
       A4qI  = 2.120902D+4 &
             - 5.179372D+3* nf &
             + 1.955772D+2* nf**2 &
             + 3.272344D+0* nf**3 
       A4ap1 = -511.228D0 + 7.08645D0*nf
       A4ap2 = -502.481D0 + 7.82077D0*nf
!
       B4qI =    2.579609D+4 + 0.08D0 &
             - ( 5.818637D+3 + 0.97D0)* nf&
             + ( 1.938554D+2 + 0.0037D0)* nf**2 &
             +   3.014982D+0* nf**3
       B4ap1 = -2426.05D0  + 266.674D0*nf - 0.05D0*nf
       B4ap2 = -2380.255D0 + 270.518D0*nf - 0.05D0*nf
!
       IF (IMOD .EQ. 1) THEN
         P3NSMC = (A4qI+A4ap1)* DL1 + B4qI+B4ap1
       ELSE IF (IMOD .EQ. 2) THEN
         P3NSMC = (A4qI+A4ap1)* DL1 + B4qI+B4ap2
       ELSE
         P3NSMC = (A4qI + 0.5D0*(A4ap1+A4ap2))* DL1 &
                 + B4qI + 0.5D0*(B4ap1+B4ap2)
       ENDIF  
!
       RETURN
     END FUNCTION P3NSMC
!
! =================================================================av==
   end module xpns3m
