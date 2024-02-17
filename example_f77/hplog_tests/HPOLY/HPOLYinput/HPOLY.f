*--------------------------------------------------------------    
      PROGRAM HPOLY
*--------------------------------------------------------------    
*          
*--------------------------------------------------------------    
*   
*     numerical evaluation of HPLs up to weight w = 8  
*   
*     J. Bluemlein, 17.09.2018   
*
*     J. Ablinger, J. Bl\"umlein, M. Round and C. Schneider
*     Numerical Implementation of Harmonic Polylogarithms to Weight w = 8
*     DESY 13-064, DO-TH 17/12, arXive:1809.07084 [hep-ph].
*   
*     Any use of this code and the files attached requires the citation of
*     the associated paper.    
*   
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER IW
*
      COMMON /WEIGHT/  IW
*
*---  default setting of IW
*     
      IW = 8
*
      WRITE(6,*) 
     &'**********************************************' 
      WRITE(6,*) 
     &'*  HPOLY.f: J. Bluemlein 17.09.2018 **********' 
      WRITE(6,*) 
     &'**********************************************' 
      WRITE(6,*) 
     &'**** HPOLY: calculation of HPLS up to w=8 ****' 
*
      CALL UHPOLYIN
* 
      CALL HPOLYIN
*
      CALL HBASIN
*
      WRITE(6,*) '**** HPOLY: initialization finished ****' 
*
      CALL UHPOLY
*  
      WRITE(6,*) 
     &'**********************************************' 
      WRITE(6,*) '**** HPOLY: calculation    finished ****' 
      WRITE(6,*) 
     &'**********************************************' 
*      
      STOP
      END
*--------------------------------------------------------------    
      SUBROUTINE UHPOLYIN
*--------------------------------------------------------------    
* 
*     USER INPUT, IW setting mandatory 
* 
*--------------------------------------------------------------    
* 
      IMPLICIT NONE
      INTEGER   IW 
      COMMON/WEIGHT/ IW
*
      IW = 8
      WRITE(6,*) '**** HPLs up to weight w =', IW
* 
      RETURN
      END  
*--------------------------------------------------------------    
      SUBROUTINE XHPOLY
*--------------------------------------------------------------    
* 
*     TEST HPL calculation
* 
*--------------------------------------------------------------    
* 
      IMPLICIT NONE
*
      REAL*8   H1,H2,H3,H4,H5,H6,H7,H8
      REAL*8   X,TT,START,END
      INTEGER  K
      INTEGER  IW,JJ 
      COMMON/WEIGHT/ IW
      COMMON /JJ/ JJ
*  
      EXTERNAL H1,H2,H3,H4,H5,H6,H7,H8 
*
      CALL CPU_TIME(START) 
      X   =  0.3D0
* 
       OPEN(UNIT=14,FILE="o1",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="WRITE")
      CALL CPU_TIME(START) 
      TT  =  H1(1,X)
      WRITE(14,*) 'T(1):=',TT,';'
      TT  =  H1(0,X) 
      WRITE(14,*) 'T(2):=',TT,';'
      TT  =  H1(-1,X) 
      WRITE(14,*) 'T(3):=',TT,';'
      CLOSE(14)
*---------------------------------------------
       OPEN(UNIT=14,FILE="o2",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="WRITE")
      DO 101 JJ = 1,3
      TT  =  H2(0,1,X)
      WRITE(14,*) 'T(',JJ,'):=',TT,';'
101   CONTINUE
      CLOSE(14)
*---------------------------------------------
       OPEN(UNIT=14,FILE="o3",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="WRITE")
      DO 102 JJ = 1,8
      TT  =  H3(0,0,1,X)
*      WRITE(14,*) X,TT
      WRITE(14,*) 'T(',JJ,'):=',TT,';'
102   CONTINUE
      CLOSE(14)
*---------------------------------------------
       OPEN(UNIT=14,FILE="o4",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="WRITE")
      DO 103 JJ = 1,18
      TT  =  H4(0,0,0,1,X)
*      WRITE(14,*) X,TT
      WRITE(14,*) 'T(',JJ,'):=',TT,';'
103   CONTINUE
      CLOSE(14)
*---------------------------------------------
       OPEN(UNIT=14,FILE="o5",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="WRITE")
      DO 104 JJ = 1,48
      TT  =  H5(0,0,0,0,1,X)
*      WRITE(14,*) X,TT
      WRITE(14,*) 'T(',JJ,'):=',TT,';'
104   CONTINUE
      CLOSE(14)
*---------------------------------------------
       OPEN(UNIT=14,FILE="o6",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="WRITE")
      DO 105 JJ = 1,116
      TT  =  H6(0,0,0,0,0,1,X)
*      WRITE(14,*) X,TT
      WRITE(14,*) 'T(',JJ,'):=',TT,';'
105   CONTINUE
      CLOSE(14)
*---------------------------------------------
       OPEN(UNIT=14,FILE="o7",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="WRITE")
      DO 106 JJ = 1,312
      TT  =  H7(0,0,0,0,0,0,1,X)
*      WRITE(14,*) X,TT
      WRITE(14,*) 'T(',JJ,'):=',TT,';'
106   CONTINUE
      CLOSE(14)
*---------------------------------------------
       OPEN(UNIT=14,FILE="o8",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="WRITE")
      DO 107 JJ = 1,810
      TT  =  H8(0,0,0,0,0,0,0,1,X)
*      WRITE(14,*) X,TT
      WRITE(14,*) 'T(',JJ,'):=',TT,';'
107   CONTINUE
      CLOSE(14)
*---------------------------------------------
      CALL CPU_TIME(END) 
      WRITE(6,*) 'CPU TIME=', END-START, 'SEC'
      STOP 
      END  
*--------------------------------------------------------------    
      SUBROUTINE UHPOLY
*--------------------------------------------------------------    
* 
*     USER HPL calculation
* 
*--------------------------------------------------------------    
* 
      IMPLICIT NONE
*
      REAL*8   H1,H2,H3,H4,H5,H6,H7,H8
      REAL*8   HH3
      REAL*8   T(8),X,TT,START,END
      INTEGER  K
      INTEGER  IW 
      COMMON/WEIGHT/ IW
*  
      EXTERNAL H1,H2,H3,H4,H5,H6,H7,H8 
*
      CALL CPU_TIME(START) 
      X   =  0.3D0
      IF(X.EQ.0) GOTO 100
* 
      T(1)  =  H1(-1,X) 
      T(2)  =  H2(0,1,X) 
      T(3)  =  H3(-1,0,0,X) 
      T(4)  =  H4(-1,-1,1,0,X) 
      T(5)  =  H5(-1,-1,1,0,1,X) 
      T(6)  =  H6(-1,0,-1,1,1,1,X)
      T(7)  =  H7(-1,-1,1,1,0,1,0,X)
      T(8)  =  H8(-1,0,-1,0,-1,0,1,1,X)
*
      DO 1 K = 1,8
      WRITE(6,*) 'K,X,T=', K,X,T(K)
1     CONTINUE
*
      CALL CPU_TIME(END)
      WRITE(6,*) 'CPU TIME=', END-START, '   SEC' 
* 
      RETURN
100   TT=0.0D0
*
      RETURN
      END  
*--------------------------------------------------------------    
      SUBROUTINE HPOLYIN
*--------------------------------------------------------------    
*
*     READ IN THE CONSTANT FIELDS 
*
*--------------------------------------------------------------    
      IMPLICIT NONE
*
      INTEGER   IW 
      INTEGER   LEN(7),K
      INTEGER   I1,I2,I3,I4,I5,I6,I7,I8,J
      REAL*8 DA2,DA3,DA4,DA5,DA6,DA7,DA8
      REAL*8 START,FINISH
*
      COMMON /WEIGHT/  IW
      COMMON /DLIST2/  DA2(360)
      COMMON /DLIST3/  DA3(1600)
      COMMON /DLIST4/  DA4(5040)
      COMMON /DLIST5/  DA5(17280)
      COMMON /DLIST6/  DA6(51040)
      COMMON /DLIST7/  DA7(162240)
      COMMON /DLIST8/  DA8(486000)

*
      DATA LEN/360,1600,5040,17280,51040,162240,486000/

*
      CALL CPU_TIME(START) 
      IF(IW.EQ.1) GOTO 10      
      IF(IW.EQ.2) GOTO 11     
      IF(IW.EQ.3) GOTO 12      
      IF(IW.EQ.4) GOTO 13      
      IF(IW.EQ.5) GOTO 14      
      IF(IW.EQ.6) GOTO 15      
      IF(IW.EQ.7) GOTO 16      
      IF(IW.EQ.8) GOTO 17      
*
      WRITE(6,*) '**** IW = ',IW,'  NOT ALLOWED, STOP: HPOLYIN ***' 
*     
10    GOTO 100
*
11    OPEN(UNIT=14,FILE="TTT2.m",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="READ")
      DO 21  K=1,LEN(1)
      READ(UNIT=14,FMT=*) DA2(K)
*     WRITE(6,*), K,DA2(K)
      DA2(K)=DA2(K)*1.0D0 
21    CONTINUE
      CLOSE(14)
      GOTO 100
*
12    OPEN(UNIT=14,FILE="TTT3.m",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="READ")
      DO 22  K=1,LEN(2)
      READ(UNIT=14,FMT=*) DA3(K)
      DA3(K)=DA3(K)*1.0D0 
*     WRITE(6,*), K,DA3(K)
22    CONTINUE
      CLOSE(14)
      GOTO 11
*
13    OPEN(UNIT=14,FILE="TTT4.m",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="READ")
      DO 23  K=1,LEN(3)
      READ(UNIT=14,FMT=*) DA4(K)
      DA4(K)=DA4(K)*1.0D0 
*     WRITE(6,*), K,DA4(K)
23    CONTINUE
      CLOSE(14)
      GOTO 12
14    OPEN(UNIT=14,FILE="TTT5.m",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="READ")
      DO 24  K=1,LEN(4)
      READ(UNIT=14,FMT=*) DA5(K)
      DA5(K)=DA5(K)*1.0D0 
*     WRITE(6,*), K,DA5(K)
24    CONTINUE
      CLOSE(14)
      GOTO 13
*
15    OPEN(UNIT=14,FILE="TTT6.m",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="READ")
      DO 25  K=1,LEN(5)
      READ(UNIT=14,FMT=*) DA6(K)
      DA6(K)=DA6(K)*1.0D0 
*     WRITE(6,*), K,DA6(K)
25    CONTINUE
      CLOSE(14)
      GOTO 14
*
16    OPEN(UNIT=14,FILE="TTT7.m",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="READ")
      DO 26  K=1,LEN(6)
      READ(UNIT=14,FMT=*) DA7(K)
      DA7(K)=DA7(K)*1.0D0 
*     WRITE(6,*), K,DA7(K)
26    CONTINUE
      CLOSE(14)
      GOTO 15
*
17    OPEN(UNIT=14,FILE="TTT8.m",FORM="FORMATTED",STATUS="OLD",
     &     ACTION="READ")
      DO 27  K=1,LEN(7)
      READ(UNIT=14,FMT=*) DA8(K)
      DA8(K)=DA8(K)*1.0D0 
*     WRITE(6,*), K,DA8(K)
27    CONTINUE
      CLOSE(14)
      GOTO 16
*
100   CONTINUE
      CALL CPU_TIME(FINISH)
      WRITE(6,*) 'READ-IN TIME=',-START+FINISH, ' SEC'  
      RETURN
      END
*--------------------------------------------------------------    
      SUBROUTINE HBASIN
*--------------------------------------------------------------    
*
*     READ IN THE CHOSEN BASIS FOR w = 5...8 
*
*--------------------------------------------------------------    
      IMPLICIT NONE
*
      INTEGER   IW,K 
      INTEGER   II15,II25,II35,II45,II55
      INTEGER   II16,II26,II36,II46,II56,II66
      INTEGER   II17,II27,II37,II47,II57,II67,II77
      INTEGER   II18,II28,II38,II48,II58,II68,II78,II88
*
      COMMON /WEIGHT/  IW
      COMMON/BAS5/ II15(48),II25(48),II35(48),II45(48),II55(48)
      COMMON/BAS6/ II16(116),II26(116),II36(116),II46(116),II56(116),
     &             II66(116)
      COMMON/BAS7/ II17(312),II27(312),II37(312),II47(312),II57(312),
     &             II67(312),II77(312)    
      COMMON/BAS8/ II18(810),II28(810),II38(810),II48(810),II58(810),
     &             II68(810),II78(810),II88(810)    
*
      IF(IW.GE.1.AND.IW.LE.4) RETURN
*
      IF(IW.EQ.5) GOTO 10      
      IF(IW.EQ.6) GOTO 11     
      IF(IW.EQ.7) GOTO 12      
      IF(IW.EQ.8) GOTO 13      
*
      WRITE(6,*) '**** IW = ',IW,'  NOT ALLOWED, STOP: HBASIN ***' 
*     
10    CONTINUE
      OPEN(UNIT=14,FILE="B51.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=15,FILE="B52.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=16,FILE="B53.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=17,FILE="B54.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=18,FILE="B55.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      DO 20  K=1,48
      READ(UNIT=14,FMT=*) II15(K)
      READ(UNIT=15,FMT=*) II25(K)
      READ(UNIT=16,FMT=*) II35(K)
      READ(UNIT=17,FMT=*) II45(K)
      READ(UNIT=18,FMT=*) II55(K)
20    CONTINUE
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      RETURN
*
11    CONTINUE
      OPEN(UNIT=14,FILE="B61.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=15,FILE="B62.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=16,FILE="B63.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=17,FILE="B64.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=18,FILE="B65.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=19,FILE="B66.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      DO 21  K=1,116
      READ(UNIT=14,FMT=*) II16(K)
      READ(UNIT=15,FMT=*) II26(K)
      READ(UNIT=16,FMT=*) II36(K)
      READ(UNIT=17,FMT=*) II46(K)
      READ(UNIT=18,FMT=*) II56(k)
      READ(UNIT=19,FMT=*) II66(k)
21    CONTINUE
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)
      GOTO  10
12    CONTINUE
      OPEN(UNIT=14,FILE="B71.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=15,FILE="B72.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=16,FILE="B73.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=17,FILE="B74.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=18,FILE="B75.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=19,FILE="B76.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=30,FILE="B77.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      DO 22  K=1,312
      READ(UNIT=14,FMT=*) II17(K)
      READ(UNIT=15,FMT=*) II27(K)
      READ(UNIT=16,FMT=*) II37(K)
      READ(UNIT=17,FMT=*) II47(K)
      READ(UNIT=18,FMT=*) II57(k)
      READ(UNIT=19,FMT=*) II67(k)
      READ(UNIT=30,FMT=*) II77(k)
22    CONTINUE
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)
      CLOSE(30)
      GOTO 11
13    CONTINUE
      OPEN(UNIT=14,FILE="B81.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=15,FILE="B82.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=16,FILE="B83.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=17,FILE="B84.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=18,FILE="B85.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=19,FILE="B86.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=30,FILE="B87.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      OPEN(UNIT=31,FILE="B88.m",FORM="FORMATTED",STATUS="OLD",
     &ACTION="READ")
      DO 23  K=1,810
      READ(UNIT=14,FMT=*) II18(K)
      READ(UNIT=15,FMT=*) II28(K)
      READ(UNIT=16,FMT=*) II38(K)
      READ(UNIT=17,FMT=*) II48(K)
      READ(UNIT=18,FMT=*) II58(k)
      READ(UNIT=19,FMT=*) II68(k)
      READ(UNIT=30,FMT=*) II78(k)
      READ(UNIT=31,FMT=*) II88(k)
23    CONTINUE
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)
      CLOSE(30)
      CLOSE(31)
      GOTO 12
*
      END
*--------------------------------------------------------------    
      SUBROUTINE HPOLYC(X)
*--------------------------------------------------------------    
*
*     CHECK whether x is in range or not  
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      REAL*8 X
*  
      IF(X.GE.0.0D0.AND.X.LE.DSQRT(2.0D0)-1.0D0) GOTO 1
*
      WRITE(6,*) 'X=',X,' *** NOT IN RANGE, STOP: HPOLYC ***'
      STOP
1     RETURN
      END
*--------------------------------------------------------------    
      REAL*8 FUNCTION H1(I1,X)
*--------------------------------------------------------------    
*
*     Weight w = 1 HPLs     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER I1
      REAL*8 X,T
*  
      IF(I1.EQ.-1) GOTO 1
      IF(I1.EQ. 0) GOTO 2
      IF(I1.EQ. 1) GOTO 3
*
      WRITE(6,*) 'I1=',I1,' *** W=1, NOT IN RANGE, STOP: H1 ***'
      STOP
*
1     T = DLOG(1.0D0+X)
      GOTO 10
2     T = DLOG(X)
      GOTO 10
3     T = -DLOG(1.0D0-X)
*
10    H1 = T
*
      RETURN
      END
*--------------------------------------------------------------    
      REAL*8 FUNCTION H2(I1,I2,X)
*--------------------------------------------------------------    
*
*     Weight w = 2 HPLs     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER I1,I2,J,K,JJ
      REAL*8 X,T1,T2,T3,U,V
      REAL*8 DA2
*
      COMMON /JJ/ JJ
      COMMON /DLIST2/  DA2(360)
*
      CALL HPOLYC(X)
*
      CALL CHECKB2(I1,I2,J)
*     J = JJ

      U= DLOG(1.0D0+X)
      V=-DLOG(1.0D0-X)
*
      T1=0.0D0
      DO 1 K=1,40
      T1=T1+U**K*DA2((J-1)*40+K)
1     CONTINUE
      T2=0.0D0
      DO 2 K=1,40
      T2=T2+U**K*DA2(120+(J-1)*40+K)
2     CONTINUE
      T2=T2*DLOG(U)
      T3=0.0D0
      DO 3 K=1,40
      T3=T3+V**K*DA2(240+(J-1)*40+K)
3     CONTINUE
*
      H2=T1+T2+T3
*
      RETURN
      END
*--------------------------------------------------------------    
      REAL*8 FUNCTION H3(I1,I2,I3,X)
*--------------------------------------------------------------    
*
*     Weight w = 3 HPLs     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER I1,I2,I3,J,K,JJ
      REAL*8 X,T1,T2,T3,T4,T5,U,V,DLU,DLV
      REAL*8 DA3
*
      COMMON /DLIST3/  DA3(1600)
      COMMON /JJ/ JJ
*
      H3=0.0D0
      CALL HPOLYC(X)
*
      CALL CHECKB3(I1,I2,I3,J)
*     J=JJ
*
      U  =  DLOG(1.0D0+X)
      V  = -DLOG(1.0D0-X)
      DLU=  DLOG(U)
      DLV=  DLOG(V)
*
      T1=0.0D0
      DO 1 K=1,40
      T1=T1+U**K*DA3((J-1)*40+K)
1     CONTINUE
      T2=0.0D0
      DO 2 K=1,40
      T2=T2+U**K*DA3(320+(J-1)*40+K)
2     CONTINUE
      T2=T2*DLU
      T3=0.0D0
      DO 3 K=1,40
      T3=T3+U**K*DA3(640+(J-1)*40+K)
3     CONTINUE
      T3=T3*DLU*DLU
      T4=0.0D0
      DO 4 K=1,40
      T4=T4+V**K*DA3(960+(J-1)*40+K)
4     CONTINUE
      T5=0.0D0
      DO 5 K=1,40
      T5=T5+V**K*DA3(1280+(J-1)*40+K)
5     CONTINUE
      T5=T5*DLV
*
      H3=T1+T2+T3+T4+T5
*      WRITE(6,*) 'H3=',H3
*
      RETURN
      END
*--------------------------------------------------------------    
      REAL*8 FUNCTION H4(I1,I2,I3,I4,X)
*--------------------------------------------------------------    
*
*     Weight w = 4 HPLs     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER I1,I2,I3,I4,J,K,JJ
      REAL*8 X,T1,T2,T3,T4,T5,T6,T7,U,V,DLU,DLV
      REAL*8 DA4
*
      COMMON /JJ/ JJ
      COMMON /DLIST4/  DA4(5040)
*
      CALL HPOLYC(X)
*
      CALL CHECKB4(I1,I2,I3,I4,J)
*     J = JJ
*
      U  =  DLOG(1.0D0+X)
      V  = -DLOG(1.0D0-X)
      DLU=  DLOG(U)
      DLV=  DLOG(V)
*
      T1=0.0D0
      DO 1 K=1,40
      T1=T1+U**K*DA4((J-1)*40+K)
1     CONTINUE
      T2=0.0D0
      DO 2 K=1,40
      T2=T2+U**K*DA4(720+(J-1)*40+K)
2     CONTINUE
      T2=T2*DLU
      T3=0.0D0
      DO 3 K=1,40
      T3=T3+U**K*DA4(1440+(J-1)*40+K)
3     CONTINUE
      T3=T3*DLU*DLU
      T4=0.0D0
      DO 4 K=1,40
      T4=T4+U**K*DA4(2160+(J-1)*40+K)
4     CONTINUE
      T4=T4*DLU*DLU*DLU
      T5=0.0D0
      DO 5 K=1,40
      T5=T5+V**K*DA4(2880+(J-1)*40+K)
5     CONTINUE
      T6=0.0D0
      DO 6 K=1,40
      T6=T6+V**K*DA4(3600+(J-1)*40+K)
6     CONTINUE
      T6=T6*DLV
      T7=0.0D0
      DO 7 K=1,40
      T7=T7+V**K*DA4(4320+(J-1)*40+K)
7     CONTINUE
      T7=T7*DLV*DLV
*
      H4=T1+T2+T3+T4+T5+T6+T7
*
      RETURN
      END
*--------------------------------------------------------------    
      REAL*8 FUNCTION H5(I1,I2,I3,I4,I5,X)
*--------------------------------------------------------------    
*
*     Weight w = 5 HPLs     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER I1,I2,I3,I4,I5,J,K,JJ
      REAL*8 X,T1,T2,T3,T4,T5,T6,T7,T8,T9,U,V,DLU,DLV
      REAL*8 DA5
*
      COMMON /DLIST5/  DA5(17280)
      COMMON /JJ/ JJ
*
      CALL HPOLYC(X)
*
      CALL CHECKB5(I1,I2,I3,I4,I5,J)
*     J = JJ
*
      U  =  DLOG(1.0D0+X)
      V  = -DLOG(1.0D0-X)
      DLU=  DLOG(U)
      DLV=  DLOG(V)
*
      T1=0.0D0
      DO 1 K=1,40
      T1=T1+U**K*DA5((J-1)*40+K)
1     CONTINUE
      T2=0.0D0
      DO 2 K=1,40
      T2=T2+U**K*DA5(1920+(J-1)*40+K)
2     CONTINUE
      T2=T2*DLU
      T3=0.0D0
      DO 3 K=1,40
      T3=T3+U**K*DA5(3840+(J-1)*40+K)
3     CONTINUE
      T3=T3*DLU*DLU
      T4=0.0D0
      DO 4 K=1,40
      T4=T4+U**K*DA5(5760+(J-1)*40+K)
4     CONTINUE
      T4=T4*DLU*DLU*DLU
      T5=0.0D0
      DO 5 K=1,40
      T5=T5+U**K*DA5(7680+(J-1)*40+K)
5     CONTINUE
      T5=T5*DLU*DLU*DLU*DLU
      T6=0.0D0
      DO 6 K=1,40
      T6=T6+V**K*DA5(9600+(J-1)*40+K)
6     CONTINUE
      T7=0.0D0
      DO 7 K=1,40
      T7=T7+V**K*DA5(11520+(J-1)*40+K)
7     CONTINUE
      T7=T7*DLV
      T8=0.0D0
      DO 8 K=1,40
      T8=T8+V**K*DA5(13440+(J-1)*40+K)
8     CONTINUE
      T8=T8*DLV*DLV
      T9=0.0D0
      DO 9 K=1,40
      T9=T9+V**K*DA5(15360+(J-1)*40+K)
9     CONTINUE
      T9=T9*DLV*DLV*DLV
*
      H5=T1+T2+T3+T4+T5+T6+T7+T8+T9
*
      RETURN
      END
*--------------------------------------------------------------    
      REAL*8 FUNCTION H6(I1,I2,I3,I4,I5,I6,X)
*--------------------------------------------------------------    
*
*     Weight w = 6 HPLs     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER I1,I2,I3,I4,I5,I6,J,K,JJ
      REAL*8 X,T1,T2,T3,T4,T5,T6,T7,T8,T9,U,V,DLU,DLV
      REAL*8 T10,T11  
      REAL*8 DA6
*
      COMMON /DLIST6/  DA6(51040)
      COMMON /JJ/ JJ
*
      CALL HPOLYC(X)
*
      CALL CHECKB6(I1,I2,I3,I4,I5,I6,J)
*     J = JJ
*
      U  =  DLOG(1.0D0+X)
      V  = -DLOG(1.0D0-X)
      DLU=  DLOG(U)
      DLV=  DLOG(V)
*
      T1=0.0D0
      DO 1 K=1,40
      T1=T1+U**K*DA6((J-1)*40+K)
1     CONTINUE
      T2=0.0D0
      DO 2 K=1,40
      T2=T2+U**K*DA6(4640+(J-1)*40+K)
2     CONTINUE
      T2=T2*DLU
      T3=0.0D0
      DO 3 K=1,40
      T3=T3+U**K*DA6(9280+(J-1)*40+K)
3     CONTINUE
      T3=T3*DLU*DLU
      T4=0.0D0
      DO 4 K=1,40
      T4=T4+U**K*DA6(13920+(J-1)*40+K)
4     CONTINUE
      T4=T4*DLU*DLU*DLU
      T5=0.0D0
      DO 5 K=1,40
      T5=T5+U**K*DA6(18560+(J-1)*40+K)
5     CONTINUE
      T5=T5*DLU*DLU*DLU*DLU
      T6=0.0D0
      DO 6 K=1,40
      T6=T6+U**K*DA6(23200+(J-1)*40+K)
6     CONTINUE
      T6=T6*DLU*DLU*DLU*DLU*DLU
      T7=0.0D0
      DO 7 K=1,40
      T7=T7+V**K*DA6(27840+(J-1)*40+K)
7     CONTINUE
      T7=T7
      T8=0.0D0
      DO 8 K=1,40
      T8=T8+V**K*DA6(32480+(J-1)*40+K)
8     CONTINUE
      T8=T8*DLV
      T9=0.0D0
      DO 9 K=1,40
      T9=T9+V**K*DA6(37120+(J-1)*40+K)
9     CONTINUE
      T9=T9*DLV*DLV
      T10=0.0D0
      DO 10 K=1,40
      T10=T10+V**K*DA6(41760+(J-1)*40+K)
10    CONTINUE
      T10=T10*DLV*DLV*DLV
      T11=0.0D0
      DO 11 K=1,40
      T11=T11+V**K*DA6(46400+(J-1)*40+K)
11    CONTINUE
      T11=T11*DLV*DLV*DLV*DLV
*
      H6=T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11
*
      RETURN
      END
*--------------------------------------------------------------    
      REAL*8 FUNCTION H7(I1,I2,I3,I4,I5,I6,I7,X)
*--------------------------------------------------------------    
*
*     Weight w = 7 HPLs     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER I1,I2,I3,I4,I5,I6,I7,J,K,JJ
      REAL*8 X,T1,T2,T3,T4,T5,T6,T7,T8,T9,U,V,DLU,DLV
      REAL*8 T10,T11,T12,T13  
      REAL*8 DA7
*
      COMMON /DLIST7/  DA7(162240)
      COMMON /JJ/ JJ
*
      CALL HPOLYC(X)
*
      CALL CHECKB7(I1,I2,I3,I4,I5,I6,I7,J)
*     J = JJ
*
      U  =  DLOG(1.0D0+X)
      V  = -DLOG(1.0D0-X)
      DLU=  DLOG(U)
      DLV=  DLOG(V)
*
      T1=0.0D0
      DO 1 K=1,40
      T1=T1+U**K*DA7((J-1)*40+K)
1     CONTINUE
      T2=0.0D0
      DO 2 K=1,40
      T2=T2+U**K*DA7(12480+(J-1)*40+K)
2     CONTINUE
      T2=T2*DLU
      T3=0.0D0
      DO 3 K=1,40
      T3=T3+U**K*DA7(24960+(J-1)*40+K)
3     CONTINUE
      T3=T3*DLU*DLU
      T4=0.0D0
      DO 4 K=1,40
      T4=T4+U**K*DA7(37440+(J-1)*40+K)
4     CONTINUE
      T4=T4*DLU*DLU*DLU
      T5=0.0D0
      DO 5 K=1,40
      T5=T5+U**K*DA7(49920+(J-1)*40+K)
5     CONTINUE
      T5=T5*DLU*DLU*DLU*DLU
      T6=0.0D0
      DO 6 K=1,40
      T6=T6+U**K*DA7(62400+(J-1)*40+K)
6     CONTINUE
      T6=T6*DLU*DLU*DLU*DLU*DLU
      T7=0.0D0
      DO 7 K=1,40
      T7=T7+U**K*DA7(74880+(J-1)*40+K)
7     CONTINUE
      T7=T7*DLU*DLU*DLU*DLU*DLU*DLU
      T8=0.0D0
      DO 8 K=1,40
      T8=T8+V**K*DA7(87360+(J-1)*40+K)
8     CONTINUE
      T8=T8
      T9=0.0D0
      DO 9 K=1,40
      T9=T9+V**K*DA7(99840+(J-1)*40+K)
9     CONTINUE
      T9=T9*DLV
      T10=0.0D0
      DO 10 K=1,40
      T10=T10+V**K*DA7(112320+(J-1)*40+K)
10    CONTINUE
      T10=T10*DLV*DLV
      T11=0.0D0
      DO 11 K=1,40
      T11=T11+V**K*DA7(124800+(J-1)*40+K)
11    CONTINUE
      T11=T11*DLV*DLV*DLV
      T12=0.0D0
      DO 12 K=1,40
      T12=T12+V**K*DA7(137280+(J-1)*40+K)
12    CONTINUE
      T12=T12*DLV*DLV*DLV*DLV
      T13=0.0D0
      DO 13 K=1,40
      T13=T13+V**K*DA7(149760+(J-1)*40+K)
13    CONTINUE
      T13=T13*DLV*DLV*DLV*DLV*DLV
*
      H7=T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11+T12+T13
*
      RETURN
      END
*--------------------------------------------------------------    
      REAL*8 FUNCTION H8(I1,I2,I3,I4,I5,I6,I7,I8,X)
*--------------------------------------------------------------    
*
*     Weight w = 8 HPLs     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,J,K,JJ
      REAL*8 X,T1,T2,T3,T4,T5,T6,T7,T8,T9,U,V,DLU,DLV
      REAL*8 T10,T11,T12,T13,T14,T15  
      REAL*8 DA8
*
      COMMON /DLIST8/  DA8(486000)
      COMMON /JJ/ JJ
*
      CALL HPOLYC(X)
*
      CALL CHECKB8(I1,I2,I3,I4,I5,I6,I7,I8,J)
*     J=JJ
*
      U  =  DLOG(1.0D0+X)
      V  = -DLOG(1.0D0-X)
      DLU=  DLOG(U)
      DLV=  DLOG(V)
*
      T1=0.0D0
      DO 1 K=1,40
      T1=T1+U**K*DA8((J-1)*40+K)
1     CONTINUE
      T2=0.0D0
      DO 2 K=1,40
      T2=T2+U**K*DA8(32400+(J-1)*40+K)
2     CONTINUE
      T2=T2*DLU
      T3=0.0D0
      DO 3 K=1,40
      T3=T3+U**K*DA8(64800+(J-1)*40+K)
3     CONTINUE
      T3=T3*DLU*DLU
      T4=0.0D0
      DO 4 K=1,40
      T4=T4+U**K*DA8(97200+(J-1)*40+K)
4     CONTINUE
      T4=T4*DLU*DLU*DLU
      T5=0.0D0
      DO 5 K=1,40
      T5=T5+U**K*DA8(129600+(J-1)*40+K)
5     CONTINUE
      T5=T5*DLU*DLU*DLU*DLU
      T6=0.0D0
      DO 6 K=1,40
      T6=T6+U**K*DA8(162000+(J-1)*40+K)
6     CONTINUE
      T6=T6*DLU*DLU*DLU*DLU*DLU
      T7=0.0D0
      DO 7 K=1,40
      T7=T7+U**K*DA8(194400+(J-1)*40+K)
7     CONTINUE
      T7=T7*DLU*DLU*DLU*DLU*DLU*DLU
      T8=0.0D0
      DO 8 K=1,40
      T8=T8+U**K*DA8(226800+(J-1)*40+K)
8     CONTINUE
      T8=T8*DLU*DLU*DLU*DLU*DLU*DLU*DLU
      T9=0.0D0
      DO 9 K=1,40
      T9=T9+V**K*DA8(259200+(J-1)*40+K)
9     CONTINUE
      T9=T9
      T10=0.0D0
      DO 10 K=1,40
      T10=T10+V**K*DA8(291600+(J-1)*40+K)
10    CONTINUE
      T10=T10*DLV
      T11=0.0D0
      DO 11 K=1,40
      T11=T11+V**K*DA8(324000+(J-1)*40+K)
11    CONTINUE
      T11=T11*DLV*DLV
      T12=0.0D0
      DO 12 K=1,40
      T12=T12+V**K*DA8(356400+(J-1)*40+K)
12    CONTINUE
      T12=T12*DLV*DLV*DLV
      T13=0.0D0
      DO 13 K=1,40
      T13=T13+V**K*DA8(388800+(J-1)*40+K)
13    CONTINUE
      T13=T13*DLV*DLV*DLV*DLV
      T14=0.0D0
      DO 14 K=1,40
      T14=T14+V**K*DA8(421200+(J-1)*40+K)
14    CONTINUE
      T14=T14*DLV*DLV*DLV*DLV*DLV
      T15=0.0D0
      DO 15 K=1,40
      T15=T15+V**K*DA8(453600+(J-1)*40+K)
15    CONTINUE
      T15=T15*DLV*DLV*DLV*DLV*DLV*DLV
*
      H8=T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11+T12+T13
     &  +T14+T15
*
      RETURN
      END
*--------------------------------------------------------------    
      SUBROUTINE CHECKB2(I1,I2,J)
*--------------------------------------------------------------    
*
*     Weight w = 2 basis check     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER II1,II2,I1,I2,J,K,L
*        
      DIMENSION II1(3),II2(3)
      DATA II1/0,-1,-1/
      DATA II2/1, 1, 0/
*
*     columns denote the HPL index pattern
*
*
      J = 0
      DO 1 K=1,3
      IF(I1.EQ.II1(K).AND.I2.EQ.II2(K)) J=K
1     CONTINUE
      IF(J.EQ.0) GOTO 10
*
      RETURN
*
10    WRITE(6,*) 'I1,I2 =',I1,I2,'not in basis STOP: H2 ****'
      STOP
      END
*--------------------------------------------------------------    
      SUBROUTINE CHECKB3(I1,I2,I3,J)
*--------------------------------------------------------------    
*
*     Weight w = 3 basis check     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER II1,II2,II3,I1,I2,I3,J,K
*        
      DIMENSION II1(8),II2(8),II3(8)
      DATA II1/0, 0,-1,-1,-1,-1,-1,-1/
      DATA II2/1, 0, 1, 1, 0, 0,-1,-1/
      DATA II3/1, 1, 1, 0, 1, 0, 1, 0/
*
*     columns denote the HPL index pattern
*
*
      J = 0
      DO 1 K=1,8
      IF(I1.EQ.II1(K).AND.I2.EQ.II2(K).AND.
     &I3.EQ.II3(K)) J=K
1     CONTINUE
      IF(J.EQ.0) GOTO 10
*
      RETURN
*
10    WRITE(6,*) 'I1,I2,I3 =',I1,I2,I3,' not in basis STOP:H3 ****'
      STOP
      END
*--------------------------------------------------------------    
      SUBROUTINE CHECKB4(I1,I2,I3,I4,J)
*--------------------------------------------------------------    
*
*     Weight w = 4 basis check     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER II1,II2,II3,II4,I1,I2,I3,I4,J,K
*        
      DIMENSION II1(18),II2(18),II3(18),II4(18)
      DATA II1/0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      DATA II2/1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1/
      DATA II3/1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0,-1, 1, 1, 0, 0,-1,-1/
      DATA II4/1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0/
*
*     columns denote the HPL index pattern
*
*
      J = 0
      DO 1 K=1,18
      IF(I1.EQ.II1(K).AND.I2.EQ.II2(K).AND.
     &I3.EQ.II3(K).AND.II4(K).EQ.I4) J=K
1     CONTINUE
      IF(J.EQ.0) GOTO 10
*
      RETURN
*
10    WRITE(6,*) 'I1,I2,I3,I4 =',I1,I2,I3,I4,'not in basis STOP:H4 ****'
      STOP
      END
*--------------------------------------------------------------    
      SUBROUTINE CHECKB5(I1,I2,I3,I4,I5,J)
*--------------------------------------------------------------    
*
*     Weight w = 5 basis check     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER II1,II2,II3,II4,II5,I1,I2,I3,I4,I5,J,K
*        
*
      COMMON/BAS5/ II1(48),II2(48),II3(48),II4(48),II5(48)
*
*     columns denote the HPL index pattern
*
*
      J = 0
      DO 1 K=1,48
      IF(I1.EQ.II1(K).AND.I2.EQ.II2(K).AND.
     &I3.EQ.II3(K).AND.II4(K).EQ.I4.AND.II5(K).EQ.I5) J=K
*      WRITE(6,*) 'J,I1,I2,I3,I4,I5 =',J,I1,I2,I3,I4,I5 
1     CONTINUE
      IF(J.EQ.0) GOTO 10
*
      RETURN
*
10    WRITE(6,*) 'I1,I2,I3,I4,I5 =',I1,I2,I3,I4,I5,
     &'not in basis STOP: CHECKB5 ****'
      STOP
      END
*--------------------------------------------------------------    
      SUBROUTINE CHECKB6(I1,I2,I3,I4,I5,I6,J)
*--------------------------------------------------------------    
*
*     Weight w = 6 basis check     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER II1,II2,II3,II4,II5,II6,I1,I2,I3,I4,I5,I6,J,K
*        
*
      COMMON/BAS6/ II1(116),II2(116),II3(116),II4(116),II5(116),
     &             II6(116)
*
*     columns denote the HPL index pattern
*
*
      J = 0
      DO 1 K=1,116
      IF(I1.EQ.II1(K).AND.I2.EQ.II2(K).AND.
     &I3.EQ.II3(K).AND.II4(K).EQ.I4.AND.II5(K).EQ.I5.
     &AND.II6(K).EQ.I6) J=K
1     CONTINUE
      IF(J.EQ.0) GOTO 10
*
      RETURN
*
10    WRITE(6,*) 'I1,I2,I3,I4,I5,I6 =',I1,I2,I3,I4,I5,I6,
     &'not in basis STOP: CHECKB6 ****'
      STOP
      END
*--------------------------------------------------------------    
      SUBROUTINE CHECKB7(I1,I2,I3,I4,I5,I6,I7,J)
*--------------------------------------------------------------    
*
*     Weight w = 7 basis check     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER II1,II2,II3,II4,II5,II6,II7,I1,I2,I3,I4,I5,I6,I7,J,K
*        
*
      COMMON/BAS7/ II1(312),II2(312),II3(312),II4(312),II5(312),
     &             II6(312),II7(312)
*
*     columns denote the HPL index pattern
*
*
      J = 0
      DO 1 K=1,312
      IF(I1.EQ.II1(K).AND.I2.EQ.II2(K).AND.
     &I3.EQ.II3(K).AND.II4(K).EQ.I4.AND.II5(K).EQ.I5.
     &AND.II6(K).EQ.I6.AND.II7(K).EQ.I7) J=K
1     CONTINUE
      IF(J.EQ.0) GOTO 10
*
      RETURN
*
10    WRITE(6,*) 'I1,I2,I3,I4,I5,I6,I7 =',
     &I1,I2,I3,I4,I5,I6,I7,
     &'not in basis STOP: CHECKB7 ****'
      STOP
      END
*--------------------------------------------------------------    
      SUBROUTINE CHECKB8(I1,I2,I3,I4,I5,I6,I7,I8,J)
*--------------------------------------------------------------    
*
*     Weight w = 8 basis check     
*
*--------------------------------------------------------------    
*
      IMPLICIT NONE
*
      INTEGER II1,II2,II3,II4,II5,II6,II7,II8,I1,I2,I3,I4,I5,I6
      INTEGER I7,I8,J,K
*        
*
      COMMON/BAS8/ II1(810),II2(810),II3(810),II4(810),II5(810),
     &             II6(810),II7(810),II8(810)
*
*     columns denote the HPL index pattern
*
*
      J = 0
      DO 1 K=1,810
      IF(I1.EQ.II1(K).AND.I2.EQ.II2(K).AND.
     &I3.EQ.II3(K).AND.II4(K).EQ.I4.AND.II5(K).EQ.I5.
     &AND.II6(K).EQ.I6.AND.II7(K).EQ.I7.AND.
     &II8(K).EQ.I8) J=K
1     CONTINUE
      IF(J.EQ.0) GOTO 10
*
      RETURN
*
10    WRITE(6,*) 'I1,I2,I3,I4,I5,I6,I7,I8 =',
     &I1,I2,I3,I4,I5,I6,I7,I8,
     &'not in basis STOP: CHECKB8 ****'
      STOP
      END
