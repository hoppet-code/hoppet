!======================================================================
! A collection of special functions -- currently taken from CERNLIB (GPL
! with restriction that military use is forbidden). Collection curently
! very limited.
!
! Does not follow the coding conventions set out for the disresum package.
!
! $Id: special_functions.f90,v 1.2 2004/09/21 18:50:24 salam Exp $
module special_functions
  private
  
  public :: ddilog, dpsipg, dgamma
  public :: wgplg !NEW: TO BE TESTED


  !--------------------------------------------
  ! want them accessible with simpler names?
  interface gamma
     module procedure dgamma
  end interface
  interface dilog
     module procedure ddilog
  end interface
  public :: dilog, gamma
  public :: psi
  

contains


  !----------------------------------------------------------------------
  ! Will need it often
  function psi(x,k)
    use types
    implicit none
    real(dp) :: psi
    real(dp), intent(in) :: x
    integer, intent(in), optional :: k
    !----------------------------------------------------------------------
    integer :: k_local
    if (present(k)) then
       k_local = k
    else
       k_local = 0
    end if
    psi = dpsipg(x,k_local)
  end function psi




!======================================================================
! BELOW: the cernlib codes converted to f90
!======================================================================
  
  !-- DDILOG aka C332 -------------------------------------
  FUNCTION DDILOG(X) 
! 1 "gen/imp64.inc" 1                                                   
!                                                                       
! imp64.inc                                                             
!                                                                       
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
! 12 "dilog64.F" 2                                                      
                                                                        
      DIMENSION C(0:19) 
                                                                        
      PARAMETER (Z1 = 1, HF = Z1/2) 
      PARAMETER (PI = 3.14159265358979324D0) 
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12) 
                                                                        
      DATA C( 0) / 0.4299669356813697D0/ 
      DATA C( 1) / 0.4097598753077105D0/ 
      DATA C( 2) /-0.0185884366014592D0/ 
      DATA C( 3) / 0.0014575108062268D0/ 
      DATA C( 4) /-0.0001430418442340D0/ 
      DATA C( 5) / 0.0000158841541880D0/ 
      DATA C( 6) /-0.0000019078959387D0/ 
      DATA C( 7) / 0.0000002419180854D0/ 
      DATA C( 8) /-0.0000000319341274D0/ 
      DATA C( 9) / 0.0000000043545063D0/ 
      DATA C(10) /-0.0000000006578480D0/ 
      DATA C(11) / 0.0000000000612098D0/ 
      DATA C(12) /-0.0000000000244332D0/ 
      DATA C(13) / 0.0000000000182256D0/ 
      DATA C(14) /-0.0000000000027007D0/ 
      DATA C(15) / 0.0000000000004042D0/ 
      DATA C(16) /-0.0000000000000610D0/ 
      DATA C(17) / 0.0000000000000093D0/ 
      DATA C(18) /-0.0000000000000014D0/ 
      DATA C(19) /+0.0000000000000002D0/ 
                                                                        
      IF(X .EQ. 1) THEN 
       H=PI6 
      ELSEIF(X .EQ. -1) THEN 
       H=-PI12 
      ELSE 
       T=-X 
       IF(T .LE. -2) THEN 
        Y=-1/(1+T) 
        S=1 
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2) 
       ELSEIF(T .LT. -1) THEN 
        Y=-1-T 
        S=-1 
        A=LOG(-T) 
        A=-PI6+A*(A+LOG(1+1/T)) 
       ELSE IF(T .LE. -HF) THEN 
        Y=-(1+T)/T 
        S=1 
        A=LOG(-T) 
        A=-PI6+A*(-HF*A+LOG(1+T)) 
       ELSE IF(T .LT. 0) THEN 
        Y=-T/(1+T) 
        S=-1 
        A=HF*LOG(1+T)**2 
       ELSE IF(T .LE. 1) THEN 
        Y=T 
        S=1 
        A=0 
       ELSE 
        Y=1/T 
        S=-1 
        A=PI6+HF*LOG(T)**2 
       ENDIF 
       H=Y+Y-1 
       ALFA=H+H 
       B1=0 
       B2=0 
       DO I = 19,0,-1 
          B0=C(I)+ALFA*B1-B2 
          B2=B1 
          B1=B0 
       END DO
       H=-(S*(B0-H*B2)+A) 
      ENDIF 
                                                                        
      DDILOG=H 
      RETURN 
      END  FUNCTION DDILOG


!======================================================================
! $Id: special_functions.f90,v 1.2 2004/09/21 18:50:24 salam Exp $           
!                                                                       
! $Log: special_functions.f90,v $
! Revision 1.2  2004/09/21 18:50:24  salam
! Various speed improvements in evaluation of grid quantities; added WGPLG to special functions -- no longer need CERNLIB linkage
!
! Revision 1.1  2001/06/27 13:40:17  gsalam
! Imported files from release-H1-1-0-7 (soon to become 1-0-8) of the disresum package
!
! Revision 1.4  2001/04/20 14:39:03  salam
! removed Id and Log entries from special functions
!
! Revision 1.3  2001/04/20 14:07:29  salam
! added new documentation figure
!
! Revision 1.2  2001/04/20 09:48:56  salam
! Added some Id keywords to files
!
! Revision 1.1  2001/04/19 15:09:16  salam
! imported all the basic files I hope!
!                                                   
! Revision 1.1.1.1  1996/04/01 15:01:54  mclareni                       
! Mathlib gen                                                           
!                                                                       
!                                                                       
      FUNCTION DGAMMA(X) 
!                                                                       
!                                                                       
!                                                                       
! imp64.inc                                                             
!                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
!                                                                       
      CHARACTER*(*) NAME 
      PARAMETER(NAME='GAMMA/DGAMMA') 
!                                                                       
      CHARACTER*80 ERRTXT 
                                                                        
      DIMENSION C(0:15) 
                                                                        
      DATA C( 0) /3.65738772508338244D0/ 
      DATA C( 1) /1.95754345666126827D0/ 
      DATA C( 2) /0.33829711382616039D0/ 
      DATA C( 3) /0.04208951276557549D0/ 
      DATA C( 4) /0.00428765048212909D0/ 
      DATA C( 5) /0.00036521216929462D0/ 
      DATA C( 6) /0.00002740064222642D0/ 
      DATA C( 7) /0.00000181240233365D0/ 
      DATA C( 8) /0.00000010965775866D0/ 
      DATA C( 9) /0.00000000598718405D0/ 
      DATA C(10) /0.00000000030769081D0/ 
      DATA C(11) /0.00000000001431793D0/ 
      DATA C(12) /0.00000000000065109D0/ 
      DATA C(13) /0.00000000000002596D0/ 
      DATA C(14) /0.00000000000000111D0/ 
      DATA C(15) /0.00000000000000004D0/ 
                                                                        
      U=X 
      IF(U .LE. 0) THEN 
       WRITE(ERRTXT,101) U 
       CALL MTLPRT(NAME,'C302.1',ERRTXT) 
       H=0 
       GO TO 9 
      ENDIF 
    8 F=1 
      IF(U .LT. 3) THEN 
       DO 1 I = 1,INT(4-U) 
       F=F/U 
    1  U=U+1 
      ELSE 
       DO 2 I = 1,INT(U-3) 
       U=U-1 
    2  F=F*U 
      END IF 
      H=U+U-7 
      ALFA=H+H 
      B1=0 
      B2=0 
      DO 3 I = 15,0,-1 
      B0=C(I)+ALFA*B1-B2 
      B2=B1 
    3 B1=B0 
                                                                        
    9 DGAMMA=F*(B0-H*B2) 
                                                                        
      RETURN 
  101 FORMAT('ARGUMENT IS NEGATIVE = ',1P,E15.1) 
      END FUNCTION DGAMMA                                


!======================================================================
!                                                                       
      FUNCTION dpsipg(X,K) 
!                                                                       
!                                                                       
! $Id: special_functions.f90,v 1.2 2004/09/21 18:50:24 salam Exp $           
!                                                                       
! $Log: special_functions.f90,v $
! Revision 1.2  2004/09/21 18:50:24  salam
! Various speed improvements in evaluation of grid quantities; added WGPLG to special functions -- no longer need CERNLIB linkage
!
! Revision 1.1  2001/06/27 13:40:17  gsalam
! Imported files from release-H1-1-0-7 (soon to become 1-0-8) of the disresum package
!
! Revision 1.4  2001/04/20 14:39:03  salam
! removed Id and Log entries from special functions
!
! Revision 1.3  2001/04/20 14:07:29  salam
! added new documentation figure
!
! Revision 1.2  2001/04/20 09:48:56  salam
! Added some Id keywords to files
!
! Revision 1.1  2001/04/19 15:09:16  salam
! imported all the basic files I hope!
!                                                   
! Revision 1.1.1.1  1996/04/01 15:02:59  mclareni                       
! Mathlib gen                                                           
!                                                                       
!                                                                       
! imp64.inc                                                             
!                                                                       
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
!                                                                       
      CHARACTER*(*) NAME 
      PARAMETER(NAME='RPSIPG/dpsipg') 
!                                                                       
      DIMENSION B(0:20,6),C(7,6),NB(6),P1(0:7),Q1(0:7),P2(0:4),Q2(0:4) 
      DIMENSION SGN(6),SGF(0:6),SGH(6) 
                                                                        
      PARAMETER (DELTA = 1D-13) 
      PARAMETER (Z1 = 1, HF = Z1/2) 
      PARAMETER (PI = 3.14159265358979324D0) 
      PARAMETER (C1 = -PI**2, C2 = 2*PI**3, C3 = 2*PI**4) 
      PARAMETER (C4 = -8*PI**5, C5 = -8*PI**6, C6 = 16*PI**7) 
                                                                        
      CHARACTER*80 ERRTXT 
                                                                        
      DATA NB /16,17,17,18,19,20/ 
      DATA SGN /-1,1,-1,1,-1,1/, SGF /1,-1,2,-6,24,-120,720/ 
      DATA SGH /-0.5D0,1,-3,12,-60,360/ 
      DATA X0 /1.46163214496836234D0/ 
                                                                        
      DATA (P1(J),Q1(J),J=0,7)                                          &
     &/ 1.35249996677263464D+4, 6.93891117537634444D-7,           &
     &  4.52856016995472897D+4, 1.97685742630467364D+4,           &
     &  4.51351684697366626D+4, 4.12551608353538323D+4,           &
     &  1.85290118185826102D+4, 2.93902871199326819D+4,           &
     &  3.32915251494069355D+3, 9.08196660748551703D+3,           &
     &  2.40680324743572018D+2, 1.24474777856708560D+3,           &
     &  5.15778920001390847D+0, 6.74291295163785938D+1,           &
     &  6.22835069189847458D-3, 1/                                   
                                                                        
      DATA (P2(J),Q2(J),J=0,4)                                          &
     &/-2.72817575131529678D-15,7.77788548522961604D+0,           &
     & -6.48157123766196510D-1, 5.46117738103215070D+1,           &
     & -4.48616543918019358D+0, 8.92920700481861370D+1,           &
     & -7.01677227766758664D+0, 3.22703493791143361D+1,           &
     & -2.12940445131010517D+0, 1/                                   
                                                                        
      DATA B( 0,1) / 0.334838697910949386D0/ 
      DATA B( 1,1) /-0.055187482048730095D0/ 
      DATA B( 2,1) / 0.004510190736011502D0/ 
      DATA B( 3,1) /-0.000365705888303721D0/ 
      DATA B( 4,1) / 0.000029434627468223D0/ 
      DATA B( 5,1) /-0.000002352776815151D0/ 
      DATA B( 6,1) / 0.000000186853176633D0/ 
      DATA B( 7,1) /-0.000000014750720184D0/ 
      DATA B( 8,1) / 0.000000001157993337D0/ 
      DATA B( 9,1) /-0.000000000090439179D0/ 
      DATA B(10,1) / 0.000000000007029627D0/ 
      DATA B(11,1) /-0.000000000000543989D0/ 
      DATA B(12,1) / 0.000000000000041925D0/ 
      DATA B(13,1) /-0.000000000000003219D0/ 
      DATA B(14,1) / 0.000000000000000246D0/ 
      DATA B(15,1) /-0.000000000000000019D0/ 
      DATA B(16,1) / 0.000000000000000001D0/ 
                                                                     
      DATA B( 0,2) /-0.112592935345473830D0/ 
      DATA B( 1,2) / 0.036557001742820941D0/ 
      DATA B( 2,2) /-0.004435942496027282D0/ 
      DATA B( 3,2) / 0.000475475854728926D0/ 
      DATA B( 4,2) /-0.000047471836382632D0/ 
      DATA B( 5,2) / 0.000004521815237353D0/ 
      DATA B( 6,2) /-0.000000416300079620D0/ 
      DATA B( 7,2) / 0.000000037338998165D0/ 
      DATA B( 8,2) /-0.000000003279914474D0/ 
      DATA B( 9,2) / 0.000000000283211377D0/ 
      DATA B(10,2) /-0.000000000024104028D0/ 
      DATA B(11,2) / 0.000000000002026297D0/ 
      DATA B(12,2) /-0.000000000000168524D0/ 
      DATA B(13,2) / 0.000000000000013885D0/ 
      DATA B(14,2) /-0.000000000000001135D0/ 
      DATA B(15,2) / 0.000000000000000092D0/ 
      DATA B(16,2) /-0.000000000000000007D0/ 
      DATA B(17,2) / 0.000000000000000001D0/ 
                                                                     
      DATA B( 0,3) / 0.076012604655110384D0/ 
      DATA B( 1,3) /-0.036257186481828739D0/ 
      DATA B( 2,3) / 0.005797202338937002D0/ 
      DATA B( 3,3) /-0.000769646513610481D0/ 
      DATA B( 4,3) / 0.000091492082189884D0/ 
      DATA B( 5,3) /-0.000010097131488364D0/ 
      DATA B( 6,3) / 0.000001055777442831D0/ 
      DATA B( 7,3) /-0.000000105929577481D0/ 
      DATA B( 8,3) / 0.000000010285494201D0/ 
      DATA B( 9,3) /-0.000000000972314310D0/ 
      DATA B(10,3) / 0.000000000089884635D0/ 
      DATA B(11,3) /-0.000000000008153171D0/ 
      DATA B(12,3) / 0.000000000000727572D0/ 
      DATA B(13,3) /-0.000000000000064010D0/ 
      DATA B(14,3) / 0.000000000000005562D0/ 
      DATA B(15,3) /-0.000000000000000478D0/ 
      DATA B(16,3) / 0.000000000000000041D0/ 
      DATA B(17,3) /-0.000000000000000003D0/ 
                                                                     
      DATA B( 0,4) /-0.077234724056994793D0/ 
      DATA B( 1,4) / 0.047867163451599467D0/ 
      DATA B( 2,4) /-0.009440702186674632D0/ 
      DATA B( 3,4) / 0.001489544740103448D0/ 
      DATA B( 4,4) /-0.000204944023348860D0/ 
      DATA B( 5,4) / 0.000025671425065297D0/ 
      DATA B( 6,4) /-0.000003001393581584D0/ 
      DATA B( 7,4) / 0.000000332766437356D0/ 
      DATA B( 8,4) /-0.000000035365412111D0/ 
      DATA B( 9,4) / 0.000000003630622927D0/ 
      DATA B(10,4) /-0.000000000362096951D0/ 
      DATA B(11,4) / 0.000000000035237509D0/ 
      DATA B(12,4) /-0.000000000003357440D0/ 
      DATA B(13,4) / 0.000000000000314068D0/ 
      DATA B(14,4) /-0.000000000000028908D0/ 
      DATA B(15,4) / 0.000000000000002623D0/ 
      DATA B(16,4) /-0.000000000000000235D0/ 
      DATA B(17,4) / 0.000000000000000021D0/ 
      DATA B(18,4) /-0.000000000000000002D0/ 
                                                                     
      DATA B( 0,5) / 0.104933034459278632D0/ 
      DATA B( 1,5) /-0.078877901652793557D0/ 
      DATA B( 2,5) / 0.018397415112159397D0/ 
      DATA B( 3,5) /-0.003352284159396504D0/ 
      DATA B( 4,5) / 0.000522878230918016D0/ 
      DATA B( 5,5) /-0.000073179785814740D0/ 
      DATA B( 6,5) / 0.000009449729612085D0/ 
      DATA B( 7,5) /-0.000001146339856723D0/ 
      DATA B( 8,5) / 0.000000132269366108D0/ 
      DATA B( 9,5) /-0.000000014646669180D0/ 
      DATA B(10,5) / 0.000000001566940742D0/ 
      DATA B(11,5) /-0.000000000162791157D0/ 
      DATA B(12,5) / 0.000000000016490345D0/ 
      DATA B(13,5) /-0.000000000001634028D0/ 
      DATA B(14,5) / 0.000000000000158807D0/ 
      DATA B(15,5) /-0.000000000000015171D0/ 
      DATA B(16,5) / 0.000000000000001427D0/ 
      DATA B(17,5) /-0.000000000000000132D0/ 
      DATA B(18,5) / 0.000000000000000012D0/ 
      DATA B(19,5) /-0.000000000000000001D0/ 
                                                                     
      DATA B( 0,6) /-0.178617622142502753D0/ 
      DATA B( 1,6) / 0.155776462200520579D0/ 
      DATA B( 2,6) /-0.041723637673831277D0/ 
      DATA B( 3,6) / 0.008597141303245400D0/ 
      DATA B( 4,6) /-0.001496227761073229D0/ 
      DATA B( 5,6) / 0.000231089608557137D0/ 
      DATA B( 6,6) /-0.000032632044778436D0/ 
      DATA B( 7,6) / 0.000004296097867090D0/ 
      DATA B( 8,6) /-0.000000534528790204D0/ 
      DATA B( 9,6) / 0.000000063478151644D0/ 
      DATA B(10,6) /-0.000000007248699714D0/ 
      DATA B(11,6) / 0.000000000800521979D0/ 
      DATA B(12,6) /-0.000000000085888793D0/ 
      DATA B(13,6) / 0.000000000008985442D0/ 
      DATA B(14,6) /-0.000000000000919356D0/ 
      DATA B(15,6) / 0.000000000000092225D0/ 
      DATA B(16,6) /-0.000000000000009090D0/ 
      DATA B(17,6) / 0.000000000000000882D0/ 
      DATA B(18,6) /-0.000000000000000084D0/ 
      DATA B(19,6) / 0.000000000000000008D0/ 
      DATA B(20,6) /-0.000000000000000001D0/ 
                                                                      
      DATA C(1,1) / 1.66666666666666667D-1/ 
      DATA C(2,1) /-3.33333333333333333D-2/ 
      DATA C(3,1) / 2.38095238095238095D-2/ 
      DATA C(4,1) /-3.33333333333333333D-2/ 
      DATA C(5,1) / 7.57575757575757576D-2/ 
      DATA C(6,1) /-2.53113553113553114D-1/ 
      DATA C(7,1) / 1.16666666666666667D0/ 
                                                                     
      DATA C(1,2) / 5.00000000000000000D-1/ 
      DATA C(2,2) /-1.66666666666666667D-1/ 
      DATA C(3,2) / 1.66666666666666667D-1/ 
      DATA C(4,2) /-3.00000000000000000D-1/ 
      DATA C(5,2) / 8.33333333333333333D-1/ 
      DATA C(6,2) /-3.29047619047619048D0/ 
      DATA C(7,2) / 1.75000000000000000D1/ 
                                                                     
      DATA C(1,3) / 2.00000000000000000D0/ 
      DATA C(2,3) /-1.00000000000000000D0/ 
      DATA C(3,3) / 1.33333333333333333D0/ 
      DATA C(4,3) /-3.00000000000000000D0/ 
      DATA C(5,3) / 1.00000000000000000D+1/ 
      DATA C(6,3) /-4.60666666666666667D+1/ 
      DATA C(7,3) / 2.80000000000000000D+2/ 
                                                                        
      DATA (C(J,4),J=1,7) /10,-7,12,-33,130,-691,4760/ 
      DATA (C(J,5),J=1,7) /60,-56,120,-396,1820,-11056,85680/ 
      DATA (C(J,6),J=1,7) /420,-504,1320,-5148,27300,-187952,1627920/ 
                                                                        
      A=ABS(X) 
      V=A 
      IX=X-DELTA 
      IF(K .LT. 0 .OR. K .GT. 6) THEN 
       H=0 
       WRITE(ERRTXT,101) K 
       CALL MTLPRT(NAME,'C316.1',ERRTXT) 
      ELSEIF(ABS(IX-X) .LE. DELTA) THEN 
       H=0 
       WRITE(ERRTXT,102) X 
       CALL MTLPRT(NAME,'C316.2',ERRTXT) 
      ELSEIF(K .EQ. 0) THEN 
       IF(A .LE. 3) THEN 
        S=0 
        IF(A .LT. HF) THEN 
         S=1/V 
         V=V+1 
        ENDIF 
        AP=P1(7) 
        AQ=Q1(7) 
        DO 11 I = 6,0,-1 
        AP=P1(I)+V*AP 
   11   AQ=Q1(I)+V*AQ 
        H=(V-X0)*AP/AQ-S 
       ELSE 
        R=1/V**2 
        AP=P2(4) 
        AQ=Q2(4) 
        DO 12 I = 3,0,-1 
        AP=P2(I)+R*AP 
   12   AQ=Q2(I)+R*AQ 
        H=LOG(V)-HF/V+AP/AQ 
       ENDIF 
       IF(X .LT. 0) H=H+1/A+PI/TAN(PI*A) 
      ELSE 
       K1=K+1 
       IF(A .LE. 10) THEN 
        IF(A .LT. 3) THEN 
         S=-1/V**K1 
         DO 1 J = 1,2-INT(A) 
         V=V+1 
    1    S=S-1/V**K1 
         V=V+1 
        ELSEIF(A .LE. 4) THEN 
         S=0 
        ELSE 
         V=V-1 
         S=1/V**K1 
         DO 5 J = 1,INT(A)-4 
         V=V-1 
    5    S=S+1/V**K1 
        ENDIF 
        H=2*V-7 
        ALFA=H+H 
        B1=0 
        B2=0 
        DO 2 J = NB(K),0,-1 
        B0=B(J,K)+ALFA*B1-B2 
        B2=B1 
    2   B1=B0 
        H=B0-H*B2+SGF(K)*S 
       ELSE 
        S=0 
        IF(A .LT. 15) THEN 
         S=1/V**K1 
         DO 3 J = 1,14-INT(A) 
         V=V+1 
    3    S=S+1/V**K1 
         V=V+1 
        ENDIF 
        R=1/V**2 
        P=R*C(7,K) 
        DO 4 J = 6,1,-1 
    4   P=R*(C(J,K)+P) 
        H=((SGF(K-1)-SGN(K)*P)*V-SGH(K))/V**K1-SGF(K)*S 
       ENDIF 
       IF(X .LT. 0) THEN 
        P=PI*A 
        IF(K .EQ. 1) THEN 
         V=C1/SIN(P)**2 
        ELSEIF(K .EQ. 2) THEN 
         V=C2*COS(P)/SIN(P)**3 
        ELSEIF(K .EQ. 3) THEN 
         S=SIN(P)**2 
         V=C3*(2*S-3)/S**2 
        ELSEIF(K .EQ. 4) THEN 
         S=SIN(P) 
         V=C4*COS(P)*(S**2-3)/S**5 
        ELSEIF(K .EQ. 5) THEN 
         S=SIN(P)**2 
         V=C5*(15-15*S+2*S**2)/S**3 
        ELSEIF(K .EQ. 6) THEN 
         S=SIN(P) 
         V=C6*COS(P)*(45-30*S**2+2*S**4)/S**7 
        ENDIF 
        H=SGN(K)*(H+V+SGF(K)/A**K1) 
       ENDIF 
      ENDIF 
                                                                        
      dpsipg=H 
                                                                        
      RETURN 
  101 FORMAT('K = ',I5,'  (< 0  OR  > 6)') 
  102 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER =',1P,E15.6) 
      END function dpsipg

!======================================================================
!# 1 "cgplg64.F"                                                        
!# 1 "<built-in>"                                                       
!# 1 "<command line>"                                                   
!# 1 "cgplg64.F"                                                        
!                                                                       
! $Id: special_functions.f90,v 1.2 2004/09/21 18:50:24 salam Exp $           
!                                                                       
! $Log: special_functions.f90,v $
! Revision 1.2  2004/09/21 18:50:24  salam
! Various speed improvements in evaluation of grid quantities; added WGPLG to special functions -- no longer need CERNLIB linkage
!                                                   
! Revision 1.1.1.1 1996/04/01 15:02:01 mclareni                         
! Mathlib gen                                                           
!                                                                       
!                                                                       
!# 1 "gen/pilot.h" 1                                                    
!# 10 "cgplg64.F" 2                                                     
                                                                        
      FUNCTION WGPLG(N,M,X)
!# 1 "gen/imp64.inc" 1                                                  
!                                                                       
! $Id: special_functions.f90,v 1.2 2004/09/21 18:50:24 salam Exp $           
!                                                                       
! $Log: special_functions.f90,v $
! Revision 1.2  2004/09/21 18:50:24  salam
! Various speed improvements in evaluation of grid quantities; added WGPLG to special functions -- no longer need CERNLIB linkage
!                                                   
! Revision 1.1.1.1 1996/04/01 15:02:59 mclareni                         
! Mathlib gen                                                           
!                                                                       
!                                                                       
! imp64.inc                                                             
!                                                                       
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
!# 13 "cgplg64.F" 2                                                     
!# 1 "gen/defc64.inc" 1                                                 
!                                                                       
! $Id: special_functions.f90,v 1.2 2004/09/21 18:50:24 salam Exp $          
!                                                                       
! $Log: special_functions.f90,v $
! Revision 1.2  2004/09/21 18:50:24  salam
! Various speed improvements in evaluation of grid quantities; added WGPLG to special functions -- no longer need CERNLIB linkage
!                                                  
! Revision 1.1.1.1 1996/04/01 15:02:59 mclareni                         
! Mathlib gen                                                           
!                                                                       
!                                                                       
! defc64.inc                                                            
!                                                                       
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      COMPLEX*16                                                        &
     & WGPLG                                                            
!# 14 "cgplg64.F" 2                                                     
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
!# 1 "gen/defc64.inc" 1                                                 
!                                                                       
! $Id: special_functions.f90,v 1.2 2004/09/21 18:50:24 salam Exp $          
!                                                                       
! $Log: special_functions.f90,v $
! Revision 1.2  2004/09/21 18:50:24  salam
! Various speed improvements in evaluation of grid quantities; added WGPLG to special functions -- no longer need CERNLIB linkage
!                                                  
! Revision 1.1.1.1 1996/04/01 15:02:59 mclareni                         
! Mathlib gen                                                           
!                                                                       
!                                                                       
! defc64.inc                                                            
!                                                                       
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      COMPLEX*16                                                        &
     & Z,I,V(0:5),SK,SJ                                                 
!# 22 "cgplg64.F" 2                                                     
!                                                                       
      CHARACTER NAME*(*) 
      CHARACTER*80 ERRTXT 
                                                                        
      PARAMETER (NAME = 'CGPLG/WGPLG') 
                                                                        
                                                                        
                                                                        
                                                                        
      DIMENSION FCT(0:4),SGN(0:4),U(0:4),S1(4,4),C(4,4),A(0:30,10) 
      DIMENSION NC(10),INDEX(31) 
                                                                        
      PARAMETER (I = (0,1)) 
      PARAMETER (Z0 = 0, Z1 = 1, HF = Z1/2, C1 = 4*Z1/3, C2 = Z1/3) 
                                                                        
      DATA FCT /1,1,2,6,24/, SGN /1,-1,1,-1,1/ 
                                                                        
      DATA S1(1,1) /1.6449340668482D0/ 
      DATA S1(1,2) /1.2020569031596D0/ 
      DATA S1(1,3) /1.0823232337111D0/ 
      DATA S1(1,4) /1.0369277551434D0/ 
      DATA S1(2,1) /1.2020569031596D0/ 
      DATA S1(2,2) /2.7058080842778D-1/ 
      DATA S1(2,3) /9.6551159989444D-2/ 
      DATA S1(3,1) /1.0823232337111D0/ 
      DATA S1(3,2) /9.6551159989444D-2/ 
      DATA S1(4,1) /1.0369277551434D0/ 
                                                                      
      DATA C(1,1) / 1.6449340668482D0/ 
      DATA C(1,2) / 1.2020569031596D0/ 
      DATA C(1,3) / 1.0823232337111D0/ 
      DATA C(1,4) / 1.0369277551434D0/ 
      DATA C(2,1) / 0.0000000000000D0/ 
      DATA C(2,2) /-1.8940656589945D0/ 
      DATA C(2,3) /-3.0142321054407D0/ 
      DATA C(3,1) / 1.8940656589945D0/ 
      DATA C(3,2) / 3.0142321054407D0/ 
      DATA C(4,1) / 0.0000000000000D0/ 
                                                                        
      DATA INDEX /1,2,3,4,6*0,5,6,7,7*0,8,9,8*0,10/ 
                                                                        
      DATA NC /24,26,28,30,22,24,26,19,22,17/ 
                                                                        
      DATA A( 0,1) / .96753215043498D0/ 
      DATA A( 1,1) / .16607303292785D0/ 
      DATA A( 2,1) / .02487932292423D0/ 
      DATA A( 3,1) / .00468636195945D0/ 
      DATA A( 4,1) / .00100162749616D0/ 
      DATA A( 5,1) / .00023200219609D0/ 
      DATA A( 6,1) / .00005681782272D0/ 
      DATA A( 7,1) / .00001449630056D0/ 
      DATA A( 8,1) / .00000381632946D0/ 
      DATA A( 9,1) / .00000102990426D0/ 
      DATA A(10,1) / .00000028357538D0/ 
      DATA A(11,1) / .00000007938705D0/ 
      DATA A(12,1) / .00000002253670D0/ 
      DATA A(13,1) / .00000000647434D0/ 
      DATA A(14,1) / .00000000187912D0/ 
      DATA A(15,1) / .00000000055029D0/ 
      DATA A(16,1) / .00000000016242D0/ 
      DATA A(17,1) / .00000000004827D0/ 
      DATA A(18,1) / .00000000001444D0/ 
      DATA A(19,1) / .00000000000434D0/ 
      DATA A(20,1) / .00000000000131D0/ 
      DATA A(21,1) / .00000000000040D0/ 
      DATA A(22,1) / .00000000000012D0/ 
      DATA A(23,1) / .00000000000004D0/ 
      DATA A(24,1) / .00000000000001D0/ 
                                                                        
      DATA A( 0,2) / .95180889127832D0/ 
      DATA A( 1,2) / .43131131846532D0/ 
      DATA A( 2,2) / .10002250714905D0/ 
      DATA A( 3,2) / .02442415595220D0/ 
      DATA A( 4,2) / .00622512463724D0/ 
      DATA A( 5,2) / .00164078831235D0/ 
      DATA A( 6,2) / .00044407920265D0/ 
      DATA A( 7,2) / .00012277494168D0/ 
      DATA A( 8,2) / .00003453981284D0/ 
      DATA A( 9,2) / .00000985869565D0/ 
      DATA A(10,2) / .00000284856995D0/ 
      DATA A(11,2) / .00000083170847D0/ 
      DATA A(12,2) / .00000024503950D0/ 
      DATA A(13,2) / .00000007276496D0/ 
      DATA A(14,2) / .00000002175802D0/ 
      DATA A(15,2) / .00000000654616D0/ 
      DATA A(16,2) / .00000000198033D0/ 
      DATA A(17,2) / .00000000060204D0/ 
      DATA A(18,2) / .00000000018385D0/ 
      DATA A(19,2) / .00000000005637D0/ 
      DATA A(20,2) / .00000000001735D0/ 
      DATA A(21,2) / .00000000000536D0/ 
      DATA A(22,2) / .00000000000166D0/ 
      DATA A(23,2) / .00000000000052D0/ 
      DATA A(24,2) / .00000000000016D0/ 
      DATA A(25,2) / .00000000000005D0/ 
      DATA A(26,2) / .00000000000002D0/ 
                                                                        
      DATA A( 0,3) / .98161027991365D0/ 
      DATA A( 1,3) / .72926806320726D0/ 
      DATA A( 2,3) / .22774714909321D0/ 
      DATA A( 3,3) / .06809083296197D0/ 
      DATA A( 4,3) / .02013701183064D0/ 
      DATA A( 5,3) / .00595478480197D0/ 
      DATA A( 6,3) / .00176769013959D0/ 
      DATA A( 7,3) / .00052748218502D0/ 
      DATA A( 8,3) / .00015827461460D0/ 
      DATA A( 9,3) / .00004774922076D0/ 
      DATA A(10,3) / .00001447920408D0/ 
      DATA A(11,3) / .00000441154886D0/ 
      DATA A(12,3) / .00000135003870D0/ 
      DATA A(13,3) / .00000041481779D0/ 
      DATA A(14,3) / .00000012793307D0/ 
      DATA A(15,3) / .00000003959070D0/ 
      DATA A(16,3) / .00000001229055D0/ 
      DATA A(17,3) / .00000000382658D0/ 
      DATA A(18,3) / .00000000119459D0/ 
      DATA A(19,3) / .00000000037386D0/ 
      DATA A(20,3) / .00000000011727D0/ 
      DATA A(21,3) / .00000000003687D0/ 
      DATA A(22,3) / .00000000001161D0/ 
      DATA A(23,3) / .00000000000366D0/ 
      DATA A(24,3) / .00000000000116D0/ 
      DATA A(25,3) / .00000000000037D0/ 
      DATA A(26,3) / .00000000000012D0/ 
      DATA A(27,3) / .00000000000004D0/ 
      DATA A(28,3) / .00000000000001D0/ 
                                                                        
      DATA A( 0,4) /1.0640521184614D0/ 
      DATA A( 1,4) /1.0691720744981D0/ 
      DATA A( 2,4) / .41527193251768D0/ 
      DATA A( 3,4) / .14610332936222D0/ 
      DATA A( 4,4) / .04904732648784D0/ 
      DATA A( 5,4) / .01606340860396D0/ 
      DATA A( 6,4) / .00518889350790D0/ 
      DATA A( 7,4) / .00166298717324D0/ 
      DATA A( 8,4) / .00053058279969D0/ 
      DATA A( 9,4) / .00016887029251D0/ 
      DATA A(10,4) / .00005368328059D0/ 
      DATA A(11,4) / .00001705923313D0/ 
      DATA A(12,4) / .00000542174374D0/ 
      DATA A(13,4) / .00000172394082D0/ 
      DATA A(14,4) / .00000054853275D0/ 
      DATA A(15,4) / .00000017467795D0/ 
      DATA A(16,4) / .00000005567550D0/ 
      DATA A(17,4) / .00000001776234D0/ 
      DATA A(18,4) / .00000000567224D0/ 
      DATA A(19,4) / .00000000181313D0/ 
      DATA A(20,4) / .00000000058012D0/ 
      DATA A(21,4) / .00000000018579D0/ 
      DATA A(22,4) / .00000000005955D0/ 
      DATA A(23,4) / .00000000001911D0/ 
      DATA A(24,4) / .00000000000614D0/ 
      DATA A(25,4) / .00000000000197D0/ 
      DATA A(26,4) / .00000000000063D0/ 
      DATA A(27,4) / .00000000000020D0/ 
      DATA A(28,4) / .00000000000007D0/ 
      DATA A(29,4) / .00000000000002D0/ 
      DATA A(30,4) / .00000000000001D0/ 
                                                                        
      DATA A( 0,5) / .97920860669175D0/ 
      DATA A( 1,5) / .08518813148683D0/ 
      DATA A( 2,5) / .00855985222013D0/ 
      DATA A( 3,5) / .00121177214413D0/ 
      DATA A( 4,5) / .00020722768531D0/ 
      DATA A( 5,5) / .00003996958691D0/ 
      DATA A( 6,5) / .00000838064065D0/ 
      DATA A( 7,5) / .00000186848945D0/ 
      DATA A( 8,5) / .00000043666087D0/ 
      DATA A( 9,5) / .00000010591733D0/ 
      DATA A(10,5) / .00000002647892D0/ 
      DATA A(11,5) / .00000000678700D0/ 
      DATA A(12,5) / .00000000177654D0/ 
      DATA A(13,5) / .00000000047342D0/ 
      DATA A(14,5) / .00000000012812D0/ 
      DATA A(15,5) / .00000000003514D0/ 
      DATA A(16,5) / .00000000000975D0/ 
      DATA A(17,5) / .00000000000274D0/ 
      DATA A(18,5) / .00000000000077D0/ 
      DATA A(19,5) / .00000000000022D0/ 
      DATA A(20,5) / .00000000000006D0/ 
      DATA A(21,5) / .00000000000002D0/ 
      DATA A(22,5) / .00000000000001D0/ 
                                                                      
      DATA A( 0,6) / .95021851963952D0/ 
      DATA A( 1,6) / .29052529161433D0/ 
      DATA A( 2,6) / .05081774061716D0/ 
      DATA A( 3,6) / .00995543767280D0/ 
      DATA A( 4,6) / .00211733895031D0/ 
      DATA A( 5,6) / .00047859470550D0/ 
      DATA A( 6,6) / .00011334321308D0/ 
      DATA A( 7,6) / .00002784733104D0/ 
      DATA A( 8,6) / .00000704788108D0/ 
      DATA A( 9,6) / .00000182788740D0/ 
      DATA A(10,6) / .00000048387492D0/ 
      DATA A(11,6) / .00000013033842D0/ 
      DATA A(12,6) / .00000003563769D0/ 
      DATA A(13,6) / .00000000987174D0/ 
      DATA A(14,6) / .00000000276586D0/ 
      DATA A(15,6) / .00000000078279D0/ 
      DATA A(16,6) / .00000000022354D0/ 
      DATA A(17,6) / .00000000006435D0/ 
      DATA A(18,6) / .00000000001866D0/ 
      DATA A(19,6) / .00000000000545D0/ 
      DATA A(20,6) / .00000000000160D0/ 
      DATA A(21,6) / .00000000000047D0/ 
      DATA A(22,6) / .00000000000014D0/ 
      DATA A(23,6) / .00000000000004D0/ 
      DATA A(24,6) / .00000000000001D0/ 
                                                                      
      DATA A( 0,7) / .95064032186777D0/ 
      DATA A( 1,7) / .54138285465171D0/ 
      DATA A( 2,7) / .13649979590321D0/ 
      DATA A( 3,7) / .03417942328207D0/ 
      DATA A( 4,7) / .00869027883583D0/ 
      DATA A( 5,7) / .00225284084155D0/ 
      DATA A( 6,7) / .00059516089806D0/ 
      DATA A( 7,7) / .00015995617766D0/ 
      DATA A( 8,7) / .00004365213096D0/ 
      DATA A( 9,7) / .00001207474688D0/ 
      DATA A(10,7) / .00000338018176D0/ 
      DATA A(11,7) / .00000095632476D0/ 
      DATA A(12,7) / .00000027313129D0/ 
      DATA A(13,7) / .00000007866968D0/ 
      DATA A(14,7) / .00000002283195D0/ 
      DATA A(15,7) / .00000000667205D0/ 
      DATA A(16,7) / .00000000196191D0/ 
      DATA A(17,7) / .00000000058018D0/ 
      DATA A(18,7) / .00000000017246D0/ 
      DATA A(19,7) / .00000000005151D0/ 
      DATA A(20,7) / .00000000001545D0/ 
      DATA A(21,7) / .00000000000465D0/ 
      DATA A(22,7) / .00000000000141D0/ 
      DATA A(23,7) / .00000000000043D0/ 
      DATA A(24,7) / .00000000000013D0/ 
      DATA A(25,7) / .00000000000004D0/ 
      DATA A(26,7) / .00000000000001D0/ 
                                                                      
      DATA A( 0,8) / .98800011672229D0/ 
      DATA A( 1,8) / .04364067609601D0/ 
      DATA A( 2,8) / .00295091178278D0/ 
      DATA A( 3,8) / .00031477809720D0/ 
      DATA A( 4,8) / .00004314846029D0/ 
      DATA A( 5,8) / .00000693818230D0/ 
      DATA A( 6,8) / .00000124640350D0/ 
      DATA A( 7,8) / .00000024293628D0/ 
      DATA A( 8,8) / .00000005040827D0/ 
      DATA A( 9,8) / .00000001099075D0/ 
      DATA A(10,8) / .00000000249467D0/ 
      DATA A(11,8) / .00000000058540D0/ 
      DATA A(12,8) / .00000000014127D0/ 
      DATA A(13,8) / .00000000003492D0/ 
      DATA A(14,8) / .00000000000881D0/ 
      DATA A(15,8) / .00000000000226D0/ 
      DATA A(16,8) / .00000000000059D0/ 
      DATA A(17,8) / .00000000000016D0/ 
      DATA A(18,8) / .00000000000004D0/ 
      DATA A(19,8) / .00000000000001D0/ 
                                                                      
      DATA A( 0,9) / .95768506546350D0/ 
      DATA A( 1,9) / .19725249679534D0/ 
      DATA A( 2,9) / .02603370313918D0/ 
      DATA A( 3,9) / .00409382168261D0/ 
      DATA A( 4,9) / .00072681707110D0/ 
      DATA A( 5,9) / .00014091879261D0/ 
      DATA A( 6,9) / .00002920458914D0/ 
      DATA A( 7,9) / .00000637631144D0/ 
      DATA A( 8,9) / .00000145167850D0/ 
      DATA A( 9,9) / .00000034205281D0/ 
      DATA A(10,9) / .00000008294302D0/ 
      DATA A(11,9) / .00000002060784D0/ 
      DATA A(12,9) / .00000000522823D0/ 
      DATA A(13,9) / .00000000135066D0/ 
      DATA A(14,9) / .00000000035451D0/ 
      DATA A(15,9) / .00000000009436D0/ 
      DATA A(16,9) / .00000000002543D0/ 
      DATA A(17,9) / .00000000000693D0/ 
      DATA A(18,9) / .00000000000191D0/ 
      DATA A(19,9) / .00000000000053D0/ 
      DATA A(20,9) / .00000000000015D0/ 
      DATA A(21,9) / .00000000000004D0/ 
      DATA A(22,9) / .00000000000001D0/ 
                                                                        
      DATA A( 0,10) / .99343651671347D0/ 
      DATA A( 1,10) / .02225770126826D0/ 
      DATA A( 2,10) / .00101475574703D0/ 
      DATA A( 3,10) / .00008175156250D0/ 
      DATA A( 4,10) / .00000899973547D0/ 
      DATA A( 5,10) / .00000120823987D0/ 
      DATA A( 6,10) / .00000018616913D0/ 
      DATA A( 7,10) / .00000003174723D0/ 
      DATA A( 8,10) / .00000000585215D0/ 
      DATA A( 9,10) / .00000000114739D0/ 
      DATA A(10,10) / .00000000023652D0/ 
      DATA A(11,10) / .00000000005082D0/ 
      DATA A(12,10) / .00000000001131D0/ 
      DATA A(13,10) / .00000000000259D0/ 
      DATA A(14,10) / .00000000000061D0/ 
      DATA A(15,10) / .00000000000015D0/ 
      DATA A(16,10) / .00000000000004D0/ 
      DATA A(17,10) / .00000000000001D0/ 
                                                                        
      IF(N .LT. 1 .OR. N .GT. 4 .OR. M .LT. 1 .OR. M .GT. 4 .OR.        &
     & N+M .GT. 5) THEN                                                 
       Z=0 
       WRITE(ERRTXT,101) N,M 
       CALL MTLPRT(NAME,'C321.1',ERRTXT) 
      ELSEIF(X .EQ. 1) THEN 
       Z=S1(N,M) 
      ELSEIF(X .GT. 2 .OR. X .LT. -1) THEN 
       X1=1/X 
       H=C1*X1+C2 
       ALFA=H+H 
       V(0)=1 
       V(1)=LOG(-X+I*Z0) 
       DO 33 L = 2,N+M 
   33 V(L)=V(1)*V(L-1)/L 
       SK=0 
       DO 34 K = 0,M-1 
       M1=M-K 
       R=X1**M1/(FCT(M1)*FCT(N-1)) 
       SJ=0 
       DO 35 J = 0,K 
       N1=N+K-J 
       L=INDEX(10*N1+M1-10) 
       B1=0 
       B2=0 
       DO 31 IT = NC(L),0,-1 
       B0=A(IT,L)+ALFA*B1-B2 
       B2=B1 
   31 B1=B0 
       Q=(FCT(N1-1)/FCT(K-J))*(B0-H*B2)*R/M1**N1 
   35 SJ=SJ+V(J)*Q 
   34 SK=SK+SGN(K)*SJ 
       SJ=0 
       DO 36 J = 0,N-1 
   36 SJ=SJ+V(J)*C(N-J,M) 
       Z=SGN(N)*SK+SGN(M)*(SJ+V(N+M)) 
      ELSEIF(X .GT. HF) THEN 
       X1=1-X 
       H=C1*X1+C2 
       ALFA=H+H 
       V(0)=1 
       U(0)=1 
       V(1)=LOG(X1+I*Z0) 
       U(1)=LOG(X) 
       DO 23 L = 2,M 
   23 V(L)=V(1)*V(L-1)/L 
       DO 26 L = 2,N 
   26 U(L)=U(1)*U(L-1)/L 
       SK=0 
       DO 24 K = 0,N-1 
       M1=N-K 
       R=X1**M1/FCT(M1) 
       SJ=0 
       DO 25 J = 0,M-1 
       N1=M-J 
       L=INDEX(10*N1+M1-10) 
       B1=0 
       B2=0 
       DO 12 IT = NC(L),0,-1 
       B0=A(IT,L)+ALFA*B1-B2 
       B2=B1 
   12 B1=B0 
       Q=SGN(J)*(B0-H*B2)*R/M1**N1 
   25 SJ=SJ+V(J)*Q 
   24 SK=SK+U(K)*(S1(M1,M)-SJ) 
       Z=SK+SGN(M)*U(N)*V(M) 
      ELSE 
       L=INDEX(10*N+M-10) 
       H=C1*X+C2 
       ALFA=H+H 
       B1=0 
       B2=0 
       DO 11 IT = NC(L),0,-1 
       B0=A(IT,L)+ALFA*B1-B2 
       B2=B1 
   11 B1=B0 
       Z=(B0-H*B2)*X**M/(FCT(M)*M**N) 
      ENDIF 
                                                                        
      WGPLG=Z 
                                                                        
                                                                        
                                                                        
                                                                        
      RETURN 
  101 FORMAT('ILLEGAL VALUES   N = ',I3,'   M = ',I3) 
      END function wgplg                                           


      !--------------------------------------------
      !-- a very abbreviated version of mtlprt
      SUBROUTINE MTLPRT(NAME,ERC,TEXT)
        CHARACTER(len=*) NAME,ERC,TEXT
        WRITE(  *,100) ERC(1:4),NAME,ERC,trim(TEXT)
        stop
100     FORMAT(7X,'***** CERN ',A,1X,A,' ERROR ',A,': ',A)
      END SUBROUTINE MTLPRT
      
end module special_functions

