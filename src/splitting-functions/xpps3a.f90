module xpps3a
  character(len=*), parameter :: name_xpps3 = "xpps3a"
contains
!
! ..File: xPps3a.f   
!
! ..The 4-loop MSbar pure-singlet qq splitting function P_{qq,ps}^(3) 
!    for the  evolution of the flavour-singlet parton distributions.
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!    P_{qq}^(3) is obtained by adding P_{qq,ns}^(3)+
!
! ..These are approximations for fixed nf = 3, 4 and 5 based on the 
!    first ten even moments together with small-x/large-x constraints.
!    The two sets spanning the error estimate are called via IMOD = 1 
!    and IMOD = 2.  Any other value of IMOD invokes their average.
!   
! ..Reference: Four-loop splitting functions in QCD - The quark-quark case - 
!
!              G. Falcioni, F. Herzog, S. Moch and A. Vogt
!              DESY 23-022, LTH 1333
!
! =====================================================================
!
!
       FUNCTION P3PSA (Y, NF, IMOD)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf
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
       x1L4cff = - 5.6460905D+1*nf + 3.6213992D+0*nf2
       x1L3cff = - 2.4755054D+2*nf + 4.0559671D+1*nf2 &
     &           - 1.5802469D+0*nf3
       y1L4cff = - 1.3168724D+1*nf
       y1L3cff = - 1.9911111D+2*nf + 1.3695473D+1*nf2
!    
! Known small-x coefficients
       bfkl1   =   1.7492273D+3*nf
       x0L6cff = - 7.5061728D+0*nf + 7.9012346D-1*nf2
       x0L5cff =   2.8549794D+1*nf + 3.7925926D+0*nf2
       x0L4cff = - 8.5480010D+2*nf + 7.7366255D+1*nf2 & 
     &           - 1.9753086D-1*nf3
!
! The resulting part of the function
       P3ps01 =&
     &      + bfkl1* DL**2*YM &
     &      + x0L6cff* DL**6 &
     &      + x0L5cff* DL**5 &
     &      + x0L4cff* DL**4 &
     &      + x1L3cff* Y1*DL1**3 &
     &      + x1L4cff* Y1*DL1**4 & 
     &      + y1L3cff* Y1*Y1*DL1**3 & 
     &      + y1L4cff* Y1*Y1*DL1**4 
!
! The selected approximations for nf = 3, 4, 5
       IF ( NF .EQ. 3 ) THEN
         P3psApp1 = P3ps01&
     &      +  67731.D0 * Y1*DL*YM &
     &      + 274100.D0 * Y1*YM &
     &      - 104493.D0 * Y1*(1.D0+2.D0*Y) &
     &      +  34403.D0 * Y1*Y*Y &
     &      + 353656.D0 * Y1*DL &
     &      +  10620.D0 * DL**2 &
     &      +  40006.D0 * DL**3 &
     &      -  7412.1D0 * Y1*DL1&
     &      -  2365.1D0 * Y1*DL1**2&
     &      +  1533.0D0 * Y1*Y1*DL1**2
         P3psApp2  = P3ps01&
     &      +  54593.D0 * Y1*DL*YM&
     &      + 179748.D0 * Y1*YM&
     &      - 195263.D0 * Y1&
     &      +  12789.D0 * Y1*Y*(1.D0+Y)&
     &      +  4700.0D0 * Y1*DL&
     &      - 103604.D0 * DL**2&
     &      -  2758.3D0 * DL**3&
     &      -  2801.2D0 * Y1*DL1&
     &      -  1986.9D0 * Y1*DL1**2&
     &      -  6005.9D0 * Y1*Y1*DL1**2
       ELSE IF ( NF .EQ. 4 ) THEN
         P3psApp1 = P3ps01&
     &      +  90154.D0 * Y1*DL*YM&       
     &      + 359084.D0 * Y1*YM    &   
     &      - 136319.D0 * Y1*(1.D0+2.D0*Y)&
     &      +  45379.D0 * Y1*Y*Y   &
     &      + 461167.D0 * Y1*DL    &
     &      +  13869.D0 * DL**2&
     &      +  52525.D0 * DL**3 &  
     &      -  7498.2D0 * Y1*DL1 &   
     &      -  2491.5D0 * Y1*DL1**2& 
     &      +  1727.2D0 * Y1*Y1*DL1**2
         P3psApp2  = P3ps01&
     &      +  72987.D0 * Y1*DL*YM&       
     &      + 235802.D0 * Y1*YM    &   
     &      - 254921.D0 * Y1   &
     &      +  17138.D0 * Y1*Y*(1.D0+Y)&
     &      +  5212.9D0 * Y1*DL    &
     &      - 135378.D0 * DL**2   &
     &      -  3350.9D0 * DL**3&
     &      -  1472.7D0 * Y1*DL1&    
     &      -  1997.2D0 * Y1*DL1**2& 
     &      -  8123.3D0 * Y1*Y1*DL1**2
       ELSE IF ( NF .EQ. 5 ) THEN
         P3psApp1 = P3ps01&
     &      + 112481.D0 * Y1*DL*YM&
     &      + 440555.D0 * Y1*YM&
     &      - 166581.D0 * Y1*(1.D0+2.D0*Y)&
     &      +  56087.D0 * Y1*Y*Y&
     &      + 562992.D0 * Y1*DL&
     &      +  16882.D0 * DL**2&
     &      +  64577.D0 * DL**3&
     &      -  6570.1D0 * Y1*DL1&
     &      -  2365.7D0 * Y1*DL1**2&
     &      +  1761.7D0 * Y1*Y1*DL1**2
         P3psApp2  = P3ps01&
     &      +  91468.D0 * Y1*DL*YM&
     &      + 289658.D0 * Y1*YM&
     &      - 311749.D0 * Y1&
     &      +  21521.D0 * Y1*Y*(1.D0+Y)&
     &      +  4908.9D0 * Y1*DL&
     &      - 165795.D0 * DL**2&
     &      -  3814.9D0 * DL**3&
     &      +   804.5D0 * Y1*DL1&
     &      -  1760.8D0 * Y1*DL1**2&
     &      -  10295.D0 * Y1*Y1*DL1**2
       ELSE
         WRITE(6,*) '  Error in function P3PSA: choice of nf   '
       END IF
!
! We return (for now) one of the two error-band boundaries
! or the present best estimate, their average
       IF ( IMOD .EQ. 1 ) THEN
         P3psA = P3psApp1
       ELSE IF ( IMOD .EQ. 2 ) THEN
         P3psA = P3psApp2
       ELSE
         P3psA = 0.5D0* ( P3psApp1 + P3psApp2 )
       END IF
!
       RETURN
     END FUNCTION P3psA
!
! =================================================================av==
   end module xpps3a
