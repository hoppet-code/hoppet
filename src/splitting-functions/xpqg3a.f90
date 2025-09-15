module xpqg3a
  character(len=*), parameter :: name_xpqg3 = "xpqg3a"
contains!
! ..File: xpqg3p.f   
!
! ..The 4-loop MSbar gluon-to-quark splitting function P_{qg}^(3) 
!    for the  evolution of the flavour-singlet parton distributions.
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!
! ..These are approximations for fixed nf = 3, 4 and 5 based on the 
!    first ten even moments together with small-x/large-x constraints.
!    The two sets indicating the error estimate are called via IMOD = 1 
!    and IMOD = 2.  Any other value of IMOD invokes their average.
!   
! ..Reference: Four-loop splitting functions in QCD 
!              - The gluon-to-quark case -
!
!              G. Falcioni, F. Herzog, S. Moch and A. Vogt
!              DESY 23-096, LTH 1345 
!
! =====================================================================
!
!
       FUNCTION P3QGA (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf,nf2,nf3
!
       YM  = 1.D0/Y
       Y1  = 1.D0-Y
       DL  = DLOG(Y)
       DL1 = DLOG(1.D0-Y)
!
       nf2 = nf*nf
       nf3 = nf*nf2
!
! Known large-x coefficients
       x1L5cff =   1.8518519D+0*nf - 4.1152263D-1*nf2
       x1L4cff =   3.5687794D+1*nf - 3.5116598D+0*nf2 &
                 - 8.2304527D-2*nf3
       y1L5cff =   2.8806584D+0*nf + 8.2304527D-1*nf2
       y1L4cff = - 4.0511391D+1*nf + 5.5418381D+0*nf2 &
                 + 1.6460905D-1*nf3  
!    
! Known small-x coefficients
       bfkl1   =   3.9357613D+3*nf
       x0L6cff = - 1.9588477D+1*nf + 2.7654321D+0*nf2
       x0L5cff =   2.1573663D+1*nf + 1.7244444D+1*nf2
       x0L4cff = - 2.8667643D+3*nf + 3.0122403D+2*nf2 &
                 + 4.1316872D+0*nf3 
!
! The resulting part of the function
       P3QG01 = &
            + bfkl1*   YM*DL**2 &
            + x0L6cff* DL**6    &
            + x0L5cff* DL**5    &
            + x0L4cff* DL**4    &
            + x1L4cff* DL1**4   &
            + x1L5cff* DL1**5   &
            + y1L4cff* Y1*DL1**4&
            + y1L5cff* Y1*DL1**5
!
! The selected approximations for nf = 3, 4, 5
       IF ( NF .EQ. 3 ) THEN
         P3qgApp1 = P3QG01       &
            + 187500.D0 * YM*DL     &        
            + 826060.D0 * YM*Y1     &  
            - 150474.D0             &
            + 226254.D0 * Y*(2.D0-Y)  &
            + 577733.D0 * DL        &
            - 180747.D0 * DL**2     &
            +  95411.D0 * DL**3      &
            +   119.8D0 * DL1**3     &
            +  7156.3D0 * DL1**2     &
            +  45790.D0 * DL1        & 
            -  95682.D0 * DL*DL1
         P3qgApp2  = P3QG01     &
            + 135000.D0 * YM*DL    &
            + 484742.D0 * YM*Y1    &
            -  11627.D0             &
            - 187478.D0 * Y*(2.D0-Y) &
            + 413512.D0 * DL       &
            -  82500.D0 * DL**2     &
            +  29987.D0 * DL**3     &
            -   850.1D0 * DL1**3    &
            -  11425.D0 * DL1**2    &
            -  75323.D0 * DL1       &
            + 282836.D0 * DL*DL1
       ELSE IF ( NF .EQ. 4 ) THEN
         P3qgApp1 = P3QG01     &
            +  250000.D0 * YM*DL   &    
            + 1089180.D0 * YM*Y1  &     
            -  241088.D0           &
            +  342902.D0 * Y*(2.D0-Y)&
            +  720081.D0 * DL      &
            -  247071.D0 * DL**2   &
            +  126405.D0 * DL**3   &
            +    272.4D0 * DL1**3   & 
            +   10911.D0 * DL1**2   &
            +   60563.D0 * DL1      &
            -  161448.D0 * DL*DL1
         P3qgApp2  = P3QG01     &
            + 180000.D0 * YM*DL    &
            + 634090.D0 * YM*Y1    &
            -  55958.D0             &
            - 208744.D0 * Y*(2.D0-Y) &
            + 501120.D0 * DL       &
            - 116073.D0 * DL**2    &
            +  39173.D0 * DL**3     &
            -  1020.8D0 * DL1**3    &
            -  13864.D0 * DL1**2    &
            - 100922.D0 * DL1      &
            + 343243.D0 * DL*DL1
       ELSE IF ( NF .EQ. 5 ) THEN
         P3qgApp1 = P3QG01      & 
            +  312500.D0 * YM*DL    &   
            + 1345700.D0 * YM*Y1   &    
            -  350466.D0            &
            +  480028.D0 * Y*(2.D0-Y) &
            +  837903.D0 * DL       &
            -  315928.D0 * DL**2    &
            +  157086.D0 * DL**3    &
            +    472.7D0 * DL1**3    &
            +   15415.D0 * DL1**2    &
            +   75644.D0 * DL1       &
            -  244869.D0 * DL*DL1
         P3qgApp2  = P3QG01     &
            + 225000.D0 * YM*DL    &
            + 776837.D0 * YM*Y1    &
            - 119054.D0            &
            - 209530.D0 * Y*(2.D0-Y) &
            + 564202.D0 * DL       &
            - 152181.D0 * DL**2    &
            +  48046.D0 * DL**3     &
            -  1143.8D0 * DL1**3    &
            -  15553.D0 * DL1**2    &
            - 126212.D0 * DL1      &
            + 385995.D0 * DL*DL1
       ELSE
         WRITE(6,*) '  Error in function P3QGA: choice of nf   '
         CALL ABORT
       END IF
!
! We return (for now) one of the two error-band representatives
! or the present best estimate, their average
       IF ( IMOD .EQ. 1 ) THEN
         P3QGA = P3qgApp1
       ELSE IF ( IMOD .EQ. 2 ) THEN
         P3QGA = P3qgApp2
       ELSE
         P3QGA = 0.5D0* ( P3qgApp1 + P3qgApp2 )
       END IF
!
       RETURN
     END FUNCTION P3QGA
!
! =================================================================av==
   end module xpqg3a
