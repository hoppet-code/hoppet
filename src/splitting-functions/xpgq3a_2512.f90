module xpgq3a_2512
  !character(len=*), parameter :: name_xpgq3 = "xpgq3a_2512"
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
       FUNCTION P3GQA_2512 (Y, NF, IMOD)
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
       x1L5cff = 1.3443073D+1 - 5.4869684D-1*nf
       x1L4cff = 3.7539831D+2 - 3.4494742D+1*nf + 8.7791495D-1*nf2
       y1L5cff = 2.2222222D+1 - 5.4869684D-1*nf
       y1L4cff = 6.6242163D+2 - 4.7992684D+1*nf + 8.7791495D-1*nf2
!    
! x^-1 small-x coeff's, Casimir scaled from P_gg (approx. for bfkl1)
       bfkl0 =   - 8.3086173D+3 / 2.25
       bfkl1 = ( - 1.0691199D+5 - nf* 9.9638304D+2 ) / 2.25
!
! Small-x double-logs with x^0
       x0L6cff =   5.2235940D+1 - 7.3744856D+0*nf
       x0L5cff = - 2.9221399D+2 + 1.8436214D+0*nf
       x0L4cff =   7.3106077D+3 - 3.7887135D+2*nf - 3.2438957D+1*nf2
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
! The selected approximations for nf = 3,...,6
       IF ( NF .EQ. 3 ) THEN
         P3gqApp1 = P3GQ01             &  
     &      + 3.5*bfkl1*  YM*DL        &
     &      - 27891.* YM*Y1            &
     &      - 309124.                  &
     &      + 1056866.* Y*(2.-Y)       &
     &      - 124735.* DL              &
     &      - 16246.* DL**2            &
     &      + 131175.* DL**3           &
     &      + 4970.1*  DL1**3          &
     &      + 60041.* DL1**2           &
     &      + 343181.* DL1             &
     &      - 958330.* DL*DL1          
         P3gqApp2 = P3GQ01             &
     &      + 7.*bfkl1* YM*DL          &
     &      - 1139334.* YM*Y1          &
     &      + 143008.                  &
     &      - 290390.* Y*(2.-Y)        &
     &      - 659492.* DL              &
     &      + 303685.* DL**2           &
     &      - 81867.*  DL**3           &
     &      + 1811.8* DL1**3           &
     &      -  465.9*  DL1**2          &
     &      - 51206.* DL1              &
     &      + 274249.* DL*DL1          
       ELSE IF ( NF .EQ. 4 ) THEN      
         P3gqApp1 = P3gq01             &
     &      + 3.5*bfkl1*  YM*DL        &
     &      - 8302.8* YM*Y1            &
     &      - 347706.                  &
     &      + 1105306.* Y*(2.-Y)       &
     &      - 127650.* DL              &
     &      - 29728.* DL**2            &
     &      + 137537.* DL**3           &
     &      + 4658.1*  DL1**3          &
     &      + 59205.* DL1**2           &
     &      + 345513.* DL1             &
     &      - 995120.* DL*DL1          
         P3gqApp2 = P3gq01             &
     &      + 7.*bfkl1* YM*DL          &
     &      - 1129822.* YM*Y1          &
     &      + 108527.                  &
     &      - 254166.* Y*(2.-Y)        &
     &      - 667254.* DL              &
     &      + 293099.* DL**2           &
     &      - 77437.*  DL**3           &
     &      + 1471.3* DL1**3           &
     &      - 1850.3*  DL1**2          &
     &      - 52451.* DL1              &
     &      + 248634.*  DL*DL1         
       ELSE IF ( NF .EQ. 5 ) THEN      
         P3gqApp1 = P3GQ01             &
     &      + 3.5*bfkl1*  YM*DL        &
     &      + 14035.* YM*Y1            &
     &      - 384003.                  &
     &      + 1152711.* Y*(2.-Y)       &
     &      - 126346.* DL              &
     &      - 42967.* DL**2            &
     &      + 144270.* DL**3           &
     &      + 4385.5*  DL1**3          &
     &      + 58688.* DL1**2           &
     &      + 348988.* DL1             &
     &      - 1031165.* DL*DL1         
         P3gqApp2  = P3GQ01            &
     &      + 7.*bfkl1* YM*DL          &
     &      - 1117561.* YM*Y1          &
     &      + 76329.                   &
     &      - 218973.* Y*(2.-Y)        &
     &      - 670799.* DL              &
     &      + 282763.* DL**2           &
     &      - 72633.*  DL**3           &
     &      + 1170.0* DL1**3           &
     &      - 2915.5*  DL1**2          &
     &      - 52548.* DL1              &
     &      + 223771.*  DL*DL1         
       ELSE IF ( NF .EQ. 6 ) THEN      
         P3gqApp1 = P3GQ01             &
     &      + 3.5*bfkl1*  YM*DL        &
     &      + 39203.* YM*Y1            &
     &      - 417914.                  &
     &      + 1199042.* Y*(2.-Y)       &
     &      - 120750.* DL              &
     &      - 55941.* DL**2            &
     &      + 151383.* DL**3           &
     &      + 4149.2*  DL1**3          &
     &      + 58466.* DL1**2           &
     &      + 353589.* DL1             &
     &      - 1066510.* DL*DL1         
         P3gqApp2  = P3GQ01            &
     &      + 7.*bfkl1* YM*DL          &
     &      - 1102470.* YM*Y1          &
     &      + 46517.                   &
     &      - 184858.* Y*(2.-Y)        &
     &      - 670056.* DL              &
     &      + 272689.* DL**2           &
     &      - 67453.*  DL**3           &
     &      + 905.0 * DL1**3           &
     &      - 3686.2*  DL1**2          &
     &      - 51523.* DL1              &
     &      + 199594.* DL*DL1
       ELSE
         WRITE(6,*) '  Error in function P3GQA: choice of nf   '
         CALL ABORT
      END IF

!
! We return (for now) one of the two error-band representatives
! or the present best estimate, their average
       IF ( IMOD .EQ. 1 ) THEN
         P3GQA_2512 = P3gqApp1
       ELSE IF ( IMOD .EQ. 2 ) THEN
         P3GQA_2512 = P3gqApp2
       ELSE
         P3GQA_2512 = 0.5* ( P3gqApp1 + P3gqApp2 )
       END IF
!
       RETURN
       END
!
! =================================================================av==
     end module xpgq3a_2512
