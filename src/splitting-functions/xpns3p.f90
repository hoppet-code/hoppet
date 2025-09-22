module xpns3p
  character(len=*), parameter :: name_xpns3p = "xpns3p"
contains
!
! ..File: xpns3p.f   
!
! ..The 4-loop MSbar splitting function P_ns^(3)+ for the evolution  
!    of non-singlet_+ combinations of quark and anti-quark densities.
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!
! ..The nf^{0,1} leading large-nc contributions and the nf^2 part
!    are high-accuracy (0.1% or better) parametrizations of the exact
!    results. The nf^3 expression is exact up to numerical truncations.
!
! ..The remaining nf^{0,1} terms are approximations based on the first
!    eight even moments together with small-x and large-x constraints.
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
!  ..Additional references for the nf^2 and nf^3 part
!        J. Davies, B. Ruijl, T. Ueda, J.Vermaseren and A. Vogt
!         arXiv:1610.07477, Nucl. Phys. B915 (2017) 335
!
! =====================================================================
!
!
! ..The regular piece of P_ns^(3)+. 
!
       FUNCTION P3NSPA (Y, NF, IMOD)
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
               - 4.305308D+2*DL**3 - 3.6D+1*DL**4 - 4./3.D+0*DL**5&
               + 8.115605D+3 - 3.079761D+3*DL1
!
! Nonleading large-n_c, nf^0 and nf^1: two approximations
!
       P3NPA01 =   3948.16D0* Y1 - 2464.61D0* (2*Y-Y*Y)*Y1&
                 - 1839.44D0* DL**2 - 402.156D0* DL**3&
                 - 1777.27D0* DL1**2*Y1 - 204.183D0 * DL1**3*Y1 + 507.152D0&
                 - 5.587553D+1*DL**4 - 2.831276D+0*DL**5&
                 - 1.488340D-1*DL**6 - 2.601749D+3 - 2.118867D+3*DL1
       P3NPA02 =  (8698.39 - 10490.47*Y)* Y*Y1&
                 + 1389.73* DL + 189.576* DL**2&
                 - 173.936* DL1**2*Y1 + 223.078* DL1**3*Y1 + 505.209&
                 - 5.587553D+1*DL**4 - 2.831276D+0*DL**5&
                 - 1.488340D-1*DL**6 - 2.601749D+3 - 2.118867D+3*DL1
!
       P3NPA11 = (-1116.34D0 + 1071.24D0*Y)* Y*Y1&
                 - 59.3041D0* DL**2 - 8.4620D0* DL**3&
                 - 143.813D0* DL1*Y1 - 18.8803D0* DL1**3*Y1 - 7.33927D0&
                 + 4.658436D+0*DL**4 + 2.798354D-1*DL**5&
                 + 3.121643D+2 + 3.379310D+2*DL1
       P3NPA12 = (-690.151D0 - 656.386D0* Y*Y)* Y1&
                 + 133.702D0* DL**2 + 34.0569D0* DL**3&
                 - 745.573D0* DL1*Y1 + 8.61438D0* DL1**3*Y1 - 7.53662D0&
                 + 4.658437D+0*DL**4 + 2.798354D-1*DL**5&
                 + 3.121643D+2 + 3.379310D+2*DL1
!
! nf^2 (parametrized) and nf^3 (exact)
!
       P3NSPA2 = 2.5D+2*  ( Y1* ( 3.0008D0 + 0.8619D0* Y - 0.12411D0* Y**2 &
               + 0.31595D0* Y**3 ) - 0.37529D0* Y*DL - 0.21684D0* Y*DL**2&
               - 0.02295D0* Y*DL**3 + 0.03394D0* Y1*DL1 + 0.40431D0 * DL*DL1 )&
               + 3.930056D+2*DL + 1.125705D+2*DL**2 + 1.652675D+1*DL**3&
               + 7.901235D-1*DL**4 - 3.760092D+2 + 2.668861D+1*DL1
!
       P3NSA3  = - 2.426296D+0 - 8.460488D-1* Y &
               + ( 5.267490D-1* DM - 3.687243D+0 + 3.160494D+0* Y )* DL&
               - ( 1.316872D+0* (DM+1.D-1) - 1.448560D+0*Y )* DL**2&
               - ( 2.633745D-1*DM - 1.31687D-1* (1.D0+Y) )* DL**3
!
! Assembly
!
       P3NSPAI = P3NSA0 + nf*P3NSA1 + nf**2*P3NSPA2 + nf**3*P3NSA3
       IF (IMOD .EQ. 1) THEN
         P3NSPA = P3NSPAI + P3NPA01 + nf* P3NPA11
       ELSE IF (IMOD .EQ. 2) THEN
         P3NSPA = P3NSPAI + P3NPA02 + nf* P3NPA12
       ELSE
         P3NSPA = P3NSPAI &
                + 0.5D0* ((P3NPA01+P3NPA02) + nf* (P3NPA11+P3NPA12))
       END IF
!
       RETURN
     END FUNCTION P3NSPA
!
! ---------------------------------------------------------------------
!
! ..The singular piece.
!
       FUNCTION P3NSPB (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf
       D1 = 1./(1.-Y)
!
       A4qI  = 2.120902D+4 &
             - 5.179372D+3* nf &
             + 1.955772D+2* nf**2 &
             + 3.272344D+0* nf**3 
       A4ap1 = -507.152D0 + 7.33927D0*nf
       A4ap2 = -505.209D0 + 7.53662D0*nf
!
       IF (IMOD .EQ. 1) THEN
         P3NSPB = (A4qI + A4ap1)* D1
       ELSE IF (IMOD .EQ. 2) THEN
         P3NSPB = (A4qI + A4ap2)* D1
       ELSE
         P3NSPB = (A4qI + 0.5D0* (A4ap1+A4ap2) )* D1
       ENDIF
!
       RETURN
     END FUNCTION P3NSPB
!
! ---------------------------------------------------------------------
!
! ..The 'local' piece.
!
       FUNCTION P3NSPC (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf
       DL1 = LOG (1.-Y)
!
       A4qI  = 2.120902D+4 &
             - 5.179372D+3* nf &
             + 1.955772D+2* nf**2 &
             + 3.272344D+0* nf**3 
       A4ap1 = -507.152D0 + 7.33927D0*nf
       A4ap2 = -505.209D0 + 7.53662D0*nf
!
       B4qI =    2.579609D+4 + 0.08D0 &
             - ( 5.818637D+3 + 0.97D0)* nf&
             + ( 1.938554D+2 + 0.0037D0)* nf**2 &
             +   3.014982D+0* nf**3 
       B4ap1 = -2405.03D0 + 267.965D0*nf
       B4ap2 = -2394.47D0 + 269.028D0*nf
!
       IF (IMOD .EQ. 1) THEN
         P3NSPC = (A4qI+A4ap1)* DL1 + B4qI+B4ap1
       ELSE IF (IMOD .EQ. 2) THEN
         P3NSPC = (A4qI+A4ap1)* DL1 + B4qI+B4ap2
       ELSE
         P3NSPC = (A4qI + 0.5D0*(A4ap1+A4ap2))* DL1 &
                 + B4qI + 0.5D0*(B4ap1+B4ap2)
       ENDIF  
!
       RETURN
     END FUNCTION P3NSPC
!
! =================================================================av==
   end module xpns3p
