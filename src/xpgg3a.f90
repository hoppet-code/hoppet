module xpgg3a
  character(len=*), parameter :: name_xpgg3 = "xpgg3a"
contains
! ..File: xPgg3a.f   
!
! ..The 4-loop MSbar gluon-gluon splitting function P_{gg}^(3) 
!    for the  evolution of the flavour-singlet parton distributions.
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!
! ..These are approximations for fixed nf = 3, 4 and 5 based on the 
!    first five even moments together with small-x/large-x constraints.
!    The two sets providing the error estimate are called via IMOD = 1 
!    and IMOD = 2.  Any other value of IMOD invokes their average.
!
! ..The distributions (in the mathematical sense) are given as in eq.
!    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
!    The name-endings A, B, and C of the functions below correspond to 
!    the kernel superscripts [2], [3], and [1] in that equation.
!   
!    Reference: Additional moments and x-space approximations
!                 of four-loop splitting functions in QCD
!               S. Moch, B. Ruijl, T. Ueda, J. Vermaseren and A. Vogt
!               DESY 23-150, Nikhef 2023-016, LTH 1354
!
! =====================================================================
!
! ..The regular part
!
       FUNCTION P3GGA (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf
       COMMON / P3GSOFT / A4gluon
!
       YM  = 1.D0/Y
       Y1  = 1.D0-Y
       DL  = DLOG(Y)
       DL1 = DLOG(1.D0-Y)
!
       nf2 = nf*nf
       nf3 = nf*nf2
!
!  The known large-x coefficients [except delta(1-x)]
!
       A4gluon =  40880.330D0     - 11714.246D0*nf &
                + 440.04876D0*nf2 + 7.3627750D0*nf3
       Ccoeff = 8.5814120D+4 - 1.3880515D+4*nf + 1.3511111D+2*nf2
       Dcoeff = 5.4482808D+4 - 4.3411337D+3*nf - 2.1333333D+1*nf2 
!    
!  The known coefficients of 1/x*ln^a x terms, a = 3,2
!
       bfkl0  = - 8.308617314D+3
       bfkl1  = - 1.069119905D+5 - 9.963830436D+2*nf
!
!  The resulting part of the function
!
       P3gg01 =               &
            + bfkl0* DL**3*YM &
            + bfkl1* DL**2*YM &
            + Ccoeff* DL1     &
            + Dcoeff - A4gluon
!
!  The selected approximations for nf = 3, 4, 5
!
       IF ( NF .EQ. 3 ) THEN
         P3ggApp1 = P3gg01        &
            + 3.4D0 *bfkl1* DL*YM    &
            - 345063.D0 * Y1*YM      &
            + 86650.D0 * (1.D0 +Y*Y)*Y1 &
            + 158160.D0 * DL         &
            - 15741.D0 * Y1*DL1**2   &
            - 9417.D0 * Y1*DL1**3  
         P3ggApp2 = P3gg01        &
            + 5.4D0 *bfkl1* DL*YM    &
            - 1265632.D0 * Y1*YM     &
            - 656644.D0 * (1.D0 +Y*Y)*Y1&
            - 1352233.D0 * DL        &
            + 203298.D0 * Y1*DL1**2  &
            + 39112.D0 * Y1*DL1**3  
       ELSE IF ( NF .EQ. 4 ) THEN
         P3ggApp1 = P3gg01        &
            + 3.4D0 *bfkl1* DL*YM    &
            - 342625.D0 * Y1*YM      &
            + 100372.D0 * (1.D0 +Y*Y)*Y1&
            + 189167.D0 * DL         &
            - 29762.D0 * Y1*DL1**2   &
            - 12102.D0 * Y1*DL1**3  
         P3ggApp2 = P3gg01        &
            + 5.4D0 *bfkl1* DL*YM    &
            - 1271540.D0 * Y1*YM     &
            - 649661.D0 * (1.D0 +Y*Y)*Y1&
            - 1334919.D0 * DL        &
            + 191263.D0 * Y1*DL1**2  &
            + 36867.D0 * Y1*DL1**3  
       ELSE IF ( NF .EQ. 5 ) THEN
         P3ggApp1 = P3gg01        &
            + 3.4D0 *bfkl1* DL*YM    &
            - 337540.D0 * Y1*YM      &
            + 119366.D0 * (1.D0 +Y*Y)*Y1&
            + 223769.D0 * DL         &
            - 45129.D0 * Y1*DL1**2   &
            - 15046.D0 * Y1*DL1**3
         P3ggApp2 = P3gg01        &
            + 5.4D0 *bfkl1* DL*YM    &
            - 1274800.D0 * Y1*YM     &
            - 637406.D0 * (1.D0 +Y*Y)*Y1&
            - 1314010.D0 * DL        &
            + 177882.D0 * Y1*DL1**2  &
            + 34362.D0 * Y1*DL1**3  
       ELSE
         WRITE(6,*) '  Error in function P3GGA: choice of nf   '
         CALL ABORT
       END IF
!
!  We return (for now) one of the two error-band representatives
!  or the present best estimate, their average
       IF ( IMOD .EQ. 1 ) THEN
         P3GGA = P3ggApp1
       ELSE IF ( IMOD .EQ. 2 ) THEN
         P3GGA = P3ggApp2
       ELSE
         P3GGA = 0.5D0 * ( P3ggApp1 + P3ggApp2 )
       END IF
!
       RETURN
     END FUNCTION P3GGA
!
! ---------------------------------------------------------------------
!
! ..The singular (soft) piece of P_gg^(3).
!   Note: A4gluon is provided by a common block set in P3GGA
!
       FUNCTION P3GGB (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER nf
       COMMON / P3GSOFT / A4gluon
       !
       nf2 = nf*nf
       nf3 = nf*nf2
       !
       ! AK Added A4gluon here in case the A piece is not beign
       ! called. Should check later if this is ever the case.
       A4gluon =  40880.330D0     - 11714.246D0*nf &
            + 440.04876D0*nf2 + 7.3627750D0*nf3
       P3GGB  = A4gluon/(1.D0-Y)
       !
       RETURN
     END FUNCTION P3GGB
!
! ---------------------------------------------------------------------
!
!
! ..The 'local' piece of P_gg^(3).
!   Note: A4gluon is provided by a common block set in P3GGA
!
       FUNCTION P3GGC (Y, NF)
!
       IMPLICIT REAL*8 (A - Z)
       INTEGER nf
       COMMON / P3GSOFT / A4gluon
!
       nf2 = nf*nf
       nf3 = nf*nf2
!
!  The coefficient of delta(1-x), also called the virtual anomalous
!  dimension. nf^0 and nf^1 are still approximate, but the errors 
!  (+- 0.3 in nf^0) are too small to be relevant in this context.
!
       B4gluon =  68587.64D0        - 18143.983D0*nf &
                + 423.81135D0*nf2 + 9.0672154D-1*nf3
!
       P3GGC = LOG (1.D0-Y)* A4gluon + B4gluon
!
       RETURN
     END FUNCTION P3GGC
!
! =================================================================av==
   end module xpgg3a
