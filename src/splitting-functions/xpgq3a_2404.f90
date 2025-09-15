module xpgq3a_2404
  !character(len=*), parameter :: name_xpgq3 = "xpgq3a_2404"
contains
!
! ..File: xPgq3a.f   
!
! ..The 4-loop MSbar quark-to-gluon splitting function P_{gq}^(3) 
!    for the  evolution of the flavour-singlet parton distributions.
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!
! ..These are approximations for fixed nf = 3, 4 and 5 based on the 
!    first ten even moments together with small-x/large-x constraints.
!    The two sets indicating the error estimate are called via IMOD = 1 
!    and IMOD = 2.  Any other value of IMOD invokes their average.
!   
! ..Reference: Four-loop splitting functions in QCD 
!              - The quark-to-gluon case -
!
!              G. Falcioni, F. Herzog, S. Moch, A. Pelloni and A. Vogt
!              ZU-TH 20/24, DESY-24-053, LTH 1367
!
! =====================================================================
!
!
       FUNCTION P3GQA_2404 (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf,nf2,nf3
!
       YM  = 1.D0/Y
       Y1  = 1.D0-Y
       DL  = DLOG(Y)
       DL1 = DLOG(Y1)
!
       nf2 = nf*nf
       nf3 = nf*nf2
!
! Known large-x coefficients
       x1L5cff = 1.3443073E+1 - 5.4869684E-1*nf
       x1L4cff = 3.7539831E+2 - 3.4494742E+1*nf + 8.7791495E-1*nf2
       y1L5cff = 2.2222222E+1 - 5.4869684E-1*nf
       y1L4cff = 6.6242163E+2 - 4.7992684E+1*nf + 8.7791495E-1*nf2
!    
! x^-1 small-x coeff's, Casimir scaled from P_gg (approx. for bfkl1)
       bfkl0 =   - 8.3086173E+3 / 2.25
       bfkl1 = ( - 1.0691199E+5 - nf* 9.9638304E+2 ) / 2.25
!
! Small-x double-logs with x^0
       x0L6cff =   5.2235940E+1 - 7.3744856E+0*nf
       x0L5cff = - 2.9221399E+2 + 1.8436214E+0*nf
       x0L4cff =   7.3106077E+3 - 3.7887135E+2*nf - 3.2438957E+1*nf2
!*
!* The resulting part of the function
       P3gq01 =                  & 
     &      + bfkl0*   YM*DL**3  &
     &      + bfkl1*   YM*DL**2  &
     &      + x0L6cff* DL**6     &
     &      + x0L5cff* DL**5     &
     &      + x0L4cff* DL**4     &
     &      + x1L4cff* DL1**4    &
     &      + x1L5cff* DL1**5    &
     &      + y1L4cff* Y1*DL1**4 &
     &      + y1L5cff* Y1*DL1**5 
!
! The selected approximations for nf = 3, 4, 5
       IF ( NF .EQ. 3 ) THEN
         P3gqApp1 = P3GQ01       &
     &      + 6.*bfkl1*  YM*DL   &   
     &      - 744384.* YM*Y1     & 
     &      + 2453640.           &
     &      - 1540404.* Y*(2.+Y) &
     &      + 1933026.* DL       &
     &      + 1142069.* DL**2    &
     &      + 162196.* DL**3     &
     &      - 2172.1 * DL1**3    &
     &      - 93264.1* DL1**2    &
     &      - 786973.* DL1       &
     &      + 875383.* Y1*DL1**2 
         P3gqApp2  = P3GQ01      &
     &      + 3.*bfkl1*  YM*DL   &
     &      + 142414.* YM*Y1     &
     &      - 326525.            &
     &      + 2159787.* Y*(2.-Y) &
     &      - 289064.* DL        &
     &      - 176358.* DL**2     &
     &      + 156541.* DL**3     &
     &      + 9016.5*  DL1**3    &
     &      + 136063.* DL1**2    &
     &      + 829482.* DL1       &
     &      - 2359050.* DL*DL1
       ELSE IF ( NF .EQ. 4 ) THEN
         P3gqApp1 = P3gq01       &
     &      + 6.*bfkl1*  YM*DL   &    
     &      - 743535.* YM*Y1     &
     &      + 2125286.           &
     &      - 1332472.* Y*(2.+Y) &
     &      + 1631173.* DL       &
     &      + 1015255.* DL**2    &
     &      + 142612.* DL**3     &
     &      - 1910.4 * DL1**3    &
     &      - 80851.*  DL1**2    &
     &      - 680219.* DL1       &
     &      + 752733.* Y1*DL1**2
         P3gqApp2 = P3gq01       &
     &      + 3.*bfkl1*  YM*DL   &
     &      + 160568.* YM*Y1     &
     &      - 361207.            &
     &      + 2048948.* Y*(2.-Y) &
     &      - 245963.* DL        &
     &      - 171312.* DL**2     &
     &      + 163099.* DL**3     &
     &      + 8132.2*  DL1**3    &
     &      + 124425.* DL1**2    &
     &      + 762435.* DL1       &
     &      - 2193335.* DL*DL1
       ELSE IF ( NF .EQ. 5 ) THEN
         P3gqApp1 = P3GQ01       &
     &      + 6.*bfkl1* YM*DL    &   
     &      - 785864.* YM*Y1     &
     &      + 285034.            &
     &      - 131648.* Y*(2.+Y)  &
     &      - 162840.* DL        &
     &      + 321220.* DL**2     &
     &      + 12688.*  DL**3     &
     &      + 1423.4 * DL1**3    &
     &      + 1278.9*  DL1**2    &
     &      - 30919.9* DL1       &
     &      + 47588.*  Y1*DL1**2
         P3gqApp2  = P3GQ01      &
     &      + 3.*bfkl1*  YM*DL   &
     &      + 177094.* YM*Y1     &
     &      - 470694.            &
     &      + 1348823.* Y*(2.-Y) &
     &      - 52985.*  DL        &
     &      - 87354.*  DL**2     &
     &      + 176885.* DL**3     &
     &      + 4748.8*  DL1**3    &
     &      + 65811.9* DL1**2    &
     &      + 396390.* DL1       &
     &      - 1190212.* DL*DL1
       ELSE
         WRITE(6,*) '  Error in function P3GQA: choice of nf   '
         CALL ABORT
       END IF
!
! We return (for now) one of the two error-band representatives
! or the present best estimate, their average
       IF ( IMOD .EQ. 1 ) THEN
         P3GQA_2404 = P3gqApp1
       ELSE IF ( IMOD .EQ. 2 ) THEN
         P3GQA_2404 = P3gqApp2
       ELSE
         P3GQA_2404 = 0.5* ( P3gqApp1 + P3gqApp2 )
       END IF
!
       RETURN
       END
!
! =================================================================av==
     end module xpgq3a_2404
