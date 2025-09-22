module xpgg3a_2410
!  character(len=*), parameter :: name_xpgg3 = "xpgg3a"
contains
!
! ..File: xPgg3a.f   
!
! ..The 4-loop MSbar gluon-gluon splitting function P_{gg}^(3) 
!    for the  evolution of the flavour-singlet parton distributions.
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!
! ..These are approximations for fixed nf = 3, 4 and 5 based on the 
!    first 10 even moments together with small-x/large-x constraints.
!    The two sets providing the error estimate are called via IMOD = 1 
!    and IMOD = 2.  Any other value of IMOD invokes their average.
!
! ..The distributions (in the mathematical sense) are given as in eq.
!    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
!    The name-endings A, B, and C of the functions below correspond to 
!    the kernel superscripts [2], [3], and [1] in that equation.
!   
! ..Reference: Four-loop splitting functions in QCD 
!                - The gluon-gluon case -
!
!              G. Falcioni, F. Herzog, S. Moch, A. Pelloni and A. Vogt
!              ZU-TH 47/24, DESY-24-144, LTH 1384
!
! =====================================================================
!
! ..The regular part
!
       FUNCTION P3GGA_2410 (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf,nf2,nf3
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
!  Known large-x coefficients [except delta(1-x)]
!
       A4gluon =  40880.330D0     - 11714.246D0*nf &
     &          + 440.04876D0*nf2 + 7.3627750D0*nf3
       Ccoeff  = 8.5814120D+4 - 1.3880515D+4*nf + 1.3511111D+2*nf2
       Dcoeff  = 5.4482808D+4 - 4.3411337D+3*nf - 2.1333333D+1*nf2 
!
       x1L4cff = 5.6460905D+1*nf - 3.6213992D+0*nf2
       x1L3cff = 2.4755054D+2*nf - 4.0559671D+1*nf2 + 1.5802469D+0*nf3
!    
!  Known small-x coefficients 
!
       bfkl0  = - 8.3086173D+3
       bfkl1  = - 1.0691199D+5 - 9.9638304D+2*nf
!
       x0L6cff =  1.44D+2 - 2.7786008D+1*nf + 7.9012346D-1*nf2
       x0L5cff = -1.44D+2 - 1.6208066D+2*nf + 1.4380247D+1*nf2
       x0L4cff =  2.6165784D+4     - 3.3447551D+3*nf &
     &          + 9.1522635D+1*nf2 - 1.9753086D-1*nf3
!
!  The resulting part of the function
!
       P3gg01 =                 &
     &      + bfkl0* DL**3*YM   &
     &      + bfkl1* DL**2*YM   &
     &      + x0L6cff* DL**6    &
     &      + x0L5cff* DL**5    &
     &      + x0L4cff* DL**4    &
     &      + Ccoeff* DL1       &
     &      + Dcoeff - A4gluon  &
     &      + x1L4cff* Y1*DL1**4&
     &      + x1L3cff* Y1*DL1**3
!
!  The selected approximations for nf = 3, 4, 5
!
       IF ( NF .EQ. 3 ) THEN
         P3ggApp1 = P3gg01        &
     &      - 421311. * Y1*DL*YM  &
     &      - 325557. * Y1*YM     &
     &      + 1679790.* Y1        &
     &      - 1456863.* Y1*Y      &
     &      + 3246307.* Y1*DL     &
     &      + 2026324.* DL*DL     &
     &      + 549188. * DL**3     &
     &      +   8337. * Y1*DL1    &
     &      +  26718. * Y1*DL1*DL1&
     &      -  27049. * Y1*Y1*DL1**3
         P3ggApp2 = P3gg01            &   
     &      - 700113. * Y1*DL*YM      &
     &      - 2300581.* Y1*YM         &
     &      + 896407. * Y1*(1.+2.*Y)  &
     &      - 162733. * Y1*Y*Y        &
     &      - 2661862.* Y1*DL         &
     &      + 196759. * DL*DL         &
     &      - 260607. * DL**3         &
     &      +  84068. * Y1*DL1        &
     &      + 346318. * Y1*DL1*DL1    &
     &      + 315725. * DL*DL1*DL1
       ELSE IF ( NF .EQ. 4 ) THEN
         P3ggApp1 = P3gg01        &
     &      - 437084. * Y1*DL*YM  &
     &      - 361570. * Y1*YM     &
     &      + 1696070.* Y1        &
     &      - 1457385.* Y1*Y      &
     &      + 3195104.* Y1*DL     &
     &      + 2009021.* DL*DL     &
     &      + 544380. * DL**3     &
     &      +  9938.  * Y1*DL1    &
     &      +  24376. * Y1*DL1*DL1&
     &      -  22143. * Y1*Y1*DL1**3
         P3ggApp2 = P3gg01           &
     &      - 706649. * Y1*DL*YM     &
     &      - 2274637.* Y1*YM        &
     &      + 836544. * Y1*(1.+2.*Y) &
     &      - 199929. * Y1*Y*Y       &
     &      - 2683760.* Y1*DL        &
     &      + 168802. * DL*DL        &
     &      - 250799. * DL**3        &
     &      +  36967. * Y1*DL1       &
     &      +  24530. * Y1*DL1*DL1   &
     &      -  71470. * Y1*Y1*DL1*DL1
       ELSE IF ( NF .EQ. 5 ) THEN
         P3ggApp1 = P3gg01        &
     &      - 439426. * Y1*DL*YM  &
     &      - 293679. * Y1*YM     &
     &      + 1916281.* Y1        &
     &      - 1615883.* Y1*Y      &
     &      + 3648786.* Y1*DL     &
     &      + 2166231.* DL*DL     &
     &      + 594588. * DL**3     &
     &      +  50406. * Y1*DL1    &
     &      +  24692. * Y1*DL1*DL1&
     &      + 174067. * Y1*Y1*DL1
         P3ggApp2 = P3gg01           &
     &      - 705978. *  Y1*DL*YM    &
     &      - 2192234.* Y1*YM        &
     &      + 1730508.* Y1*Y         &
     &      + 353143. * Y1*(2.-Y*Y)  &
     &      - 2602682.* Y1*DL        &
     &      + 178960. * DL*DL        &
     &      - 218133. * DL**3        &
     &      +   2285. * Y1*DL1       &
     &      +  19295. * Y1*DL1*DL1   &
     &      -  13719. * Y1*Y1*DL1*DL1
       ELSE
         WRITE(6,*) '  Error in function P3ggA: choice of nf   '
         CALL ABORT
       END IF
!
!  We return one of the two error-band representatives or the present 
!  best estimate, their average
       IF ( IMOD .EQ. 1 ) THEN
         P3GGA_2410 = P3ggApp1
       ELSE IF ( IMOD .EQ. 2 ) THEN
         P3GGA_2410 = P3ggApp2
       ELSE
         P3GGA_2410 = 0.5* ( P3ggApp1 + P3ggApp2 )
       END IF
!
       RETURN
       END
!
! ---------------------------------------------------------------------
!
! ..The singular (soft) piece of P_gg^(3).
!   Note: A4gluon is provided by a common block set in P3GGA
!
       FUNCTION P3GGB_2410 (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER nf
       COMMON / P3GSOFT / A4gluon
       !
       nf2 = nf*nf
       nf3 = nf*nf2
       !!
       A4gluon =  40880.330D0     - 11714.246D0*nf &
            &          + 440.04876D0*nf2 + 7.3627750D0*nf3
       P3GGB_2410  = A4gluon/(1.D0-Y)
!
       RETURN
       END
!
! ---------------------------------------------------------------------
!
!
! ..The 'local' piece of P_gg^(3).
!   Note: A4gluon is provided by a common block set in P3GGA
!
       FUNCTION P3GGC_2410 (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A - Z)
       INTEGER IMOD, nf,nf2,nf3
       COMMON / P3GSOFT / A4gluon
!
       nf2 = nf*nf
       nf3 = nf*nf2
!
!  The coefficient of delta(1-x), also called the virtual anomalous
!  dimension. nf^0 and nf^1 are still approximate, but the error at 
!  nf^1 is far too small to be relevant in this context.
!
       B4gluon =  68587.64        - 18143.983D0*nf &
     &          + 423.81135D0*nf2 + 9.0672154D-1*nf3
       IF ( IMOD .EQ. 1 ) B4gluon = B4gluon - 0.2
       IF ( IMOD .EQ. 2 ) B4gluon = B4gluon + 0.2
!
       P3GGC_2410 = DLOG (1.D0-Y)* A4gluon + B4gluon
!
       RETURN
       END
!
! =================================================================av==
end module xpgg3a_2410
