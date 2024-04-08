      MODULE c2nsreg
      CONTAINS

       FUNCTION C2NP3A_new (Y, DL, NF, CC)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       INTEGER CC ! charged current
       DIMENSION FL(6)
       DATA FL / -1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
       DOUBLE PRECISION largex, smallx, polyfit
*     
       FL11 = FL(NF)
*
       Y1 = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
       Y1 = 1d0 - Y
       DL1 = log(Y1)
*
       C2NP3A_new = 0.D0
       IF (CC.EQ.1) THEN
! The small x expansion. Coefficients are exact but truncated.
          smallx =  -787.1751967077007d0*DL - 899.5533538033412d0*DL**2
     $         -309.0785726886318d0*DL**3 - 36.19753086419753d0*DL**4
     $         -1.1851851851851851d0*DL**5 + (278.85556475324074d0*DL
     $         +187.42856459809343d0*DL**2 +50.304526748971185d0*DL**3
     $         +2.995884773662551d0*DL**4)*NF +(-8.121621121342653d0*DL
     $         -8.16460905349794d0*DL**2 - 1.5144032921810697d0*DL**3)
     $         *NF**2
          
! The large x expansion. Coefficients are exact but truncated.
          largex = 4988.494928013526d0*DL1 + 2319.6558207170438d0*DL1**2
     $         -1787.041827321787d0*DL1**3 +348.44444444444446d0*DL1**4
     $         -18.962962962962962d0*DL1**5 +(1199.6906563815378d0*DL1
     $         -787.5420539087113d0*DL1**2 + 146.69958847736626d0*DL1**3
     $         -7.901234567901234d0*DL1**4)*NF + (-65.15652656374833d0
     $         *DL1 +14.617283950617283d0*DL1**2 -0.7901234567901234d0
     $         *DL1**3)*NF**2

!     The polynomial fit
          polyfit = -5029.59226644672d0  -2410.28652896786d0*y 
     $         -144568.737096439d0*y**2 + 751258.100128793d0*y**3 
     $         -2086451.59207511d0*y**4 + 3145469.93778892d0*y**5 
     $         -2426120.54237392d0*y**6 + 753341.689545131d0*y**7 + nf
     $         *(845.207816391045d0 + 2367.4533133017d0*y
     $         +11823.2183636921d0*y**2 -88108.2103000565d0*y**3
     $         +264008.049094827d0*y**4 -417245.388478629d0*y**5
     $         +332058.821281488d0*y**6 -105879.821765296d0*y**7)+ nf**2
     $         *(-0.421306270142914d0 -109.688202074944d0*y
     $         +140.499130847744d0*y**2 + 669.634880823259d0*y**3
     $         -2834.96889483387d0*y**4 + 5102.08622034399d0*y**5
     $         -4341.02063218853d0*y**6 + 1448.15389768751d0*y**7)

!          polyfit = -800.378502722104d0  -3614.72704727467d0*y +
!     $         6577.16919408954d0*y**2 -15399.4977106532d0*y**3 +
!     $         17194.339313936d0*y1*DL1 + DL*DL1*(11183.5198653907d0 +
!     $         1223.96715526525d0*DL)

!          print*, y,dl,smallx,largex
          
          C2NP3A_new = largex + smallx + polyfit
       ELSE
          C2NP3A_new = FL11*NF * ( ( 126.42D0 - 50.29D0 * Y - 50.15D0 *
     $         Y**2) * Y1- 26.717D0 - 960.D0*D243 * DL**2 * (DL+5.D0)
     $         +59.59D0 * DL- Y*DL**2 * (101.8D0 + 34.79D0 * DL +
     $         3.070D0* DL**2)- 9.075D0 * Y*Y1*DL1 ) * Y
*
       ENDIF
       RETURN
       END FUNCTION

       FUNCTION C2NP3A_nopoly (Y, DL, NF, CC)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
       INTEGER CC ! charged current
       DIMENSION FL(6)
       DATA FL / -1.d0, 0.5d0, 0.d0, 0.5d0, 0.2d0, 0.5d0 /
       DOUBLE PRECISION largex, smallx, polyfit
*     
       FL11 = FL(NF)
*
       Y1 = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
       Y1 = 1d0 - Y
       DL1 = log(Y1)
*
       C2NP3A_nopoly = 0.D0
       IF (CC.EQ.1) THEN
! The small x expansion. Coefficients are exact but truncated.
          smallx =  -787.1751967077007d0*DL - 899.5533538033412d0*DL**2
     $         -309.0785726886318d0*DL**3 - 36.19753086419753d0*DL**4
     $         -1.1851851851851851d0*DL**5 + (278.85556475324074d0*DL
     $         +187.42856459809343d0*DL**2 +50.304526748971185d0*DL**3
     $         +2.995884773662551d0*DL**4)*NF +(-8.121621121342653d0*DL
     $         -8.16460905349794d0*DL**2 - 1.5144032921810697d0*DL**3)
     $         *NF**2
          
! The large x expansion. Coefficients are exact but truncated.
          largex = 4988.494928013526d0*DL1 + 2319.6558207170438d0*DL1**2
     $         -1787.041827321787d0*DL1**3 +348.44444444444446d0*DL1**4
     $         -18.962962962962962d0*DL1**5 +(1199.6906563815378d0*DL1
     $         -787.5420539087113d0*DL1**2 + 146.69958847736626d0*DL1**3
     $         -7.901234567901234d0*DL1**4)*NF + (-65.15652656374833d0
     $         *DL1 +14.617283950617283d0*DL1**2 -0.7901234567901234d0
     $         *DL1**3)*NF**2

!     The polynomial fit
          polyfit = 0d0

!          print*, y,dl,smallx,largex
          
          C2NP3A_nopoly = largex + smallx + polyfit
       ELSE
          C2NP3A_nopoly = FL11*NF * ( ( 126.42D0 - 50.29D0 * Y - 50.15D0
     $         *Y**2) * Y1- 26.717D0 - 960.D0*D243 * DL**2 * (DL+5.D0)
     $         +59.59D0 * DL- Y*DL**2 * (101.8D0 + 34.79D0 * DL +3.070D0
     $         * DL**2)- 9.075D0 * Y*Y1*DL1 ) * Y
*
       ENDIF
       RETURN
       END FUNCTION

* ..For Y values close to 1, use the series expansion close 
*   to Y1=0 instead of full value, for numerical convergence
*     
       FUNCTION Y1VAL (Y, DL)
       IMPLICIT REAL*8 (A - Z)
       DOUBLE PRECISION D2, D6, D24, D120
       PARAMETER (D2 = 0.5D0, D6 = 1.0D0/6.0D0, D24 = 1.0D0/24.0D0, D120
     $      = 1.0D0/120.0D0)
       IF (ABS(DL).LT.1D-4) THEN
          Y1VAL = - DL - DL**2*D2 - DL**3*D6 - DL**4*D24 - DL**5*D120
       ELSE
          Y1VAL = 1.0D0 - Y
       ENDIF

       RETURN
       END FUNCTION


* ..For Y values close to 1, use the series expansion close 
*   to Y1=0 instead of full value, for numerical convergence
*     
       FUNCTION DL1VAL (Y, DL)
       IMPLICIT REAL*8 (A - Z)
       DOUBLE PRECISION D2, D24, D2880
       PARAMETER (D2=0.5D0, D24=1.0D0/24.0D0, D2880=1.0D0/2880.0D0)

       IF (ABS(DL).LT.1D-4) THEN
          DL1VAL = LOG(-DL) +  DL*D2 + DL**2*D24 - DL**4*D2880
       ELSE
          DL1VAL = LOG(1.0D0 - Y)
       ENDIF

       RETURN
       END FUNCTION

* ..For Y values close to 1, use the series expansion close 
*   to Y1=0 instead of full value, for numerical convergence
*     
       FUNCTION DMVAL (Y, DL)
       IMPLICIT REAL*8 (A - Z)
       DOUBLE PRECISION D2, D12, D720, D30240
       PARAMETER (D2 = 0.5D0, D12 = 1.0D0/12.0D0, D720 = 1.0D0/720.0D0,
     $      D30240=1.0D0/30240.0D0)

       IF (ABS(DL).LT.1D-3) THEN
          DMVAL  = D2 - 1.0D0/DL - DL*D12 + DL**3*D720 - DL**5*D30240
       ELSE
          DMVAL = 1.0D0/(1.0D0-Y)
       ENDIF

       RETURN
       END FUNCTION

* =================================================================av==
      END MODULE c2nsreg
