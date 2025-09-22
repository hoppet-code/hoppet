!=======================================================================
!
! This source file has been provided to the Hoppet authors under GPLv3
! by Johannes Blumlein and corresponds to the following article, which
! should be cited whenever any routine in this source file is called:
!
!     DESY 24--037   
!     J. Ablinger, A. Behring, J. Bl\"umlein, A. De Freitas,
!     A. von Manteuffel, C.~Schneider, and K.~Sch\"onwald,
!     "The Variable Flavor Number Scheme at  Next-to-Next-to-Leading Order"
!
!     >>>>> Version of:   March 19, 2024 + May 30, 2024 update
!     >>>>> Author:       J. Bluemlein
!     >>>>> compilation:  gfortran VFNS.f
!     >>>>> f90 part of code: K. Sch\"onwald, March 15, 2024
!
! Accompanying files include AQg3mod.f90 and momnsmod.f
! 
! The code as distributed here with HOPPET has been adapted as follows
!
! - the main program has been turned into a subroutine VFNS3sub.
! - so has the program ODDNS (ODDNSsub)
! - the code has been placed in a module VFNSmod, to avoid filling the
!   the global namespace with short function names that may conflict with
!   user functions
! - the original code was distributed with the DAIND integration routine
!   by R. Piessens, which has was used for numerical integrations of the
!   the z-space functions, to be compared to direct moment calculations
!   for the threshold terms. DAIND has been replaced with a call to a
!   HOPPET CERNLib DGAUSS variant (wrapped through the daind_fake
!   module); the original calls to DAIND used a precision of 1e-18, that
!   has been reduced to 1e-10 or 1e-8.
!
!=========================================================================
!     Some comments (2025-02-11 AK):
!
!     *) I had to change the JB interface in blumlein_helpers. Not sure
!     that is what we want.
      
      module VFNSmod
c     For accessing the fortran90 routines that contain the regular part
c     of aQg3
      use aqg3mod
      use daind_fake
      
      private
c     See below in VFNS3sub for some guidance towards what is in each
c     piece
      public VFNS3sub

c     non-singlet pieces
      PUBLIC NSPLU, ANSPLU1, ANSPLU2, NSREG, ANSREG, NSDEL, ANSDEL
c     gluon pieces
      PUBLIC GGPLU, AGGPLU, GGREG, AGGREG0, AGGREG1, AGGDEL, GGDEL
c     pure singlet pieces
      PUBLIC APS1, APS2, PSL, PS
c     qg pieces
      PUBLIC QG, AQG3NF, RED0, RED1, QGL, AQG
c     gq pieces
      PUBLIC GQ, AGQ

c     For accessing moments 
      PUBLIC MOMAGQ, MOMAPS, MOMGG, MOMGQ, MOMNS, MOMPS, MOMPSL, MOMQG,
     $     MOMQGL

      PUBLIC MPS

!      interface
!!     The integration routine DAIND which is needed by the above moment
!!     routines.
!      DOUBLE PRECISION FUNCTION DAIND(A,B,FUN,EPS,KEY,MAX,KOUNT,EST)
!        IMPLICIT REAL*8 (A-H,O-Z)
!        interface 
!          DOUBLE PRECISION FUNCTION FUN(X)
!          IMPLICIT REAL*8 (A-H,O-Z)
!          END FUNCTION FUN
!        end interface
!      end function DAIND
!      end interface


      contains
      
*     ------------------------------------------------------------------------      
      subroutine VFNS3sub(z_in, nflightdp_in, a4pi_in, LL_in)
      implicit none
      real*8, intent(in) :: z_in, nflightdp_in, a4pi_in, LL_in
*     
*     ------------------------------------------------------------------------
*     The single mass variable flavor number scheme at 3-loop order
*     [unpolarized case].
*
*     Any use of this code or parts of it shall quote:
*     ------------------------------------------------     
*     DESY 24--037   
*     J. Ablinger, A. Behring, J. Bl\"umlein, A. De Freitas,
*     A. von Manteuffel, C.~Schneider, and K.~Sch\"onwald,
*     "The Variable Flavor Number Scheme at  Next-to-Next-to-Leading Order"
*
*     >>>>> Version of:   March 19, 2024:
*     >>>>> Author:       J. Bluemlein
*     >>>>> compilation:  gfortran VFNS.f
*     
*     The used integration routine shall also be quoted:
*     CF. R. PIESSENS, ANGEW. INFORMATIK, VOL. 9 (1973) 399.
*    
*     For evaluating the HPLs we use (in part):     
*     T. Gehrmann and E. Remiddi, "Numerical evaluation of harmonic polylogarithms",
*     Comput. Phys. Commun. {\bf 141} (2001) 296--312 [hep-ph/0107173].
*
*     Code: J. Bluemlein, March 15, 2024
*
*     aQg3 from: J. Ablinger, A. Behring, J. Bl\"mlein, A. De Freitasd,
*     A. von Manteuffel, C. Schneider, and K. Sch\"onwald
*     The non-first-order-factorizable contributions to the three-loop
*     single-mass operator matrix elements AQg and ∆AQg
*     Phys.Lett.B 854 (2024) 138713,  arXiv:2403.00513v1 [hep-ph] 
*
*     f90 part of Code: K. Sch\"onwald, March 15, 2024
*     [quote:
*      M. Fael, F. Lange, K. Sch\"onwald and M. Steinhauser,
*      "Singlet and nonsinglet three-loop massive form factors",
*      Phys. Rev. D \textbf{106} (2022) no.3, 034029 [arXiv:2207.00027 [hep-ph]].
*
*
*     Compilation:
*     ------------
*     gfortran VFNS.f                      [as for the present version]
*     gfortran -oVFNS AQg3-test.f90 VFNS.f [for a later, extended version]
*
*-----------------------------------------------------------------------
*     This code is available from https://www-zeuthen.desy.de/~blumlein
*     The licensing conditions are such, that it shall not be 
*     distributed with other software, but it shall be loaded as 
*     separate file in all applications.
*-----------------------------------------------------------------------
*      
*     The main VFNS routines in the single mass case:      
*
*      nf                  - number of light flavors
*      as = gs^2/(16 Pi^2) - strong coupling constant
*      LL = LOG(mQ^2/Q^2)  - phase space logarithm due to the heavy quark
*
*      QGL(z,nf,as,LL)    : A_qg,Q
*      PSL(z,nf,as,LL)    : A_qq,Q^PS
*      GQ(z,nf,as,LL)     : A_gq,Q     without a_gq,Q^(3)
*      QG(z,nf,as,LL)     : A_Qg       without a_Qg^(3)
*      PS(z,nf,as,LL)     : A_Qq^PS    without a_Qq^(3,PS)
*      GGDEL(z,nf,as,LL)  : A_gg,Q     without a_gg,Q^(3):    delta-part
*      GGPLU(z,nf,as,LL)  : A_gg,Q     without a_gg,Q^(3):    plus-part
*      GGREG(z,nf,as,LL)  : A_gg,Q     without a_gg,Q^(3):    regular part
*      NSDEL(z,nf,as,LL)  : A_qq,Q^NS  without a_qq,Q^(3,NS): delta-part
*      NSPLU(z,nf,as,LL)  : A_qq,Q^NS  without a_qq,Q^(3,NS): plus-part
*      NSREG(z,nf,as,LL)  : A_qq,Q^NS  without a_qq,Q^(3,NS): regular part
*      
*
*      AGQ(z,nf,as,LL)     : a_gq,Q^(3) 
*      AQG(z,nf,as,LL)     : a_Qg^(3)
*      APS1(z,nf,as,LL)    : a_Qq^(3,PS)  :  first piece  [formerly read "z \in [0,1/2]", GPS+AK verified with JB that was incorrect]
*      APS2(z,nf,as,LL)    : a_Qq^(3,PS)  :  second piece [formerly read "z \in [1/2,1]", GPS+AK verified with JB that was incorrect]
*      AGGDEL(z,nf,as,LL)  : a_gg,Q^(3,NS):  delta-part 
*      AGGPLU(z,nf,as,LL)  : a_gg,Q^(3,NS):  plus-part
*      AGGREG0(z,nf,as,LL) : a_gg,Q^(3,NS):  regular part z \in [0,1/2]
*      AGGREG1(z,nf,as,LL) : a_gg,Q^(3,NS):  regular part z \in [1/2,1]
*      ANSDEL(z,nf,as,LL)  : a_qq,Q^(3,NS):  delta-part
*      ANSPLU1(z,nf,as,LL) : a_qq,Q^(3,NS):  plus-part: one
*      ANSPLU2(z,nf,as,LL) : a_qq,Q^(3,NS):  plus-part: two
*      ANSREG(z,nf,as,LL)  : a_qq,Q^(3,NS):  regular part
*      
*-----------------------------------------------------------------------
*      
      REAL*8 T40,CONVAqqQNS
      REAL*8 T41,CONVAggQ
      REAL*8 T30,CONVAqqQPS
      REAL*8 T31,CONVAqgQ
      REAL*8 T32,CONVAQqPS
      REAL*8 T33,CONVAgq
      REAL*8 T34,CONVAAQQPS
      REAL*8 T35,CONVAAGQQ
      REAL*8 T36,CONVAQg
      REAL*8 T37,CONVAAQG
      REAL*8 T17,T18!RED0,RED1
!      REAL*8 aQg3,T16,AAQG3!,AQG3NF
      REAL*8 T16,AAQG3!,AQG3NF
      REAL*8 nf,as,z,T1,T2,T3!,QGL,PSL,GQ
      REAL*8 T4,T5a,T5b,T5c!,QG,GGDEL,GGPLU,GGREG
      REAL*8 T6,T7a,T7b,T7c!,PS,NSDEL,NSPLU,NSREG
      REAL*8 LL
      REAL*8 T8a,T8b,T8c,T8d!,AGGDEL,AGGPLU,AGGREG0,AGGREG1
      REAL*8 T9a,T9b!,APS1,APS2
      REAL*8 T10a,T10b1,T10b2,T10c!,ANSDEL,ANSPLU1,ANSPLU2
      ! REAL*8 ANSREG
      REAL*8 T11,T12,T13!,NSDIS,SIGDIS,GLUDIS
      REAL*8 PI,alphas,XMC2,Q2,XMC
      REAL*8 CSIGAqqQ,x,T20
      REAL*8 T14,T15!,AQG,AGQ
      REAL*8 U1,U2,U3,U4,U5,U7
      !REAL*8 MOMPSL,MOMQGL,MOMGQ,MOMAGQ,MOMAPS
      REAL*8 U6!,MOMQG,MOMPS
      REAL*8 REF!,MOMGQF,MOMPSLF,MOMQGLF,MOMPSF,MOMGGF 
      !REAL*8 MOMGG,MOMQGF,MOMQQNS,MOMNS
      INTEGER INLO,IQ2,N,K
      COMMON/PILOT/INLO,IQ2
      COMMON/VAR1/ nf,as,LL
!      EXTERNAL aQg3
*
*
      WRITE(6,*)
     &'***********************************************************'
      WRITE(6,*) '                 '
      WRITE(6,*)
     & '>>> Code VFNS <<<: single-mass VFNS to 3 loop order'
      WRITE(6,*)
     & 'Please quote the references given in the code at any use'
      WRITE(6,*)
     & 'preprint: DESY 24--037'
      WRITE(6,*)
     & 'based on work by: J. Ablinger, A. Behring, J. Bl\"umlein,'
      WRITE(6,*) 
     & 'A. De Freitas, A. von Manteuffel, C.~Schneider,' 
      WRITE(6,*) 
     & 'and K.~Sch\"onwald'
      WRITE(6,*) '                 '
      WRITE(6,*)
     &'***********************************************************'
*
*
*     ------------------------------------------------------------------------
*     Setting parameters and options:
*     ------------------------------- 
*     General:
*     alphas, as:  strong coupling constant at Q2
*     Q2        :  matching scale
*     XMC       :  heavy quark mass in the OMS [not MSbar scheme]
*                  default case: charm
*     INLO      :  QCD order, currently only INLO = 2:  NNLO
*     PDFs      :  NSDIS    
*                  SIGDIS
*                  GLUDIS
*     - effective fits in x for:
*     Q2 =    30 GeV^2  IQ2=1      
*     Q2 =   100 GeV^2  IQ2=2      
*     Q2 = 10000 GeV^2  IQ2=3      
*
*     z momentum fraction: operator matrix elements (OMEs)
*     x Bjorken-x in the case that convolutions are calculated      
*     nf - number of massless flavors [default nf=3]
*     LL = Log[XMC^2/Q2]  massive log in the OMEs
*     ------------------------------------------------------------------------
*
*      WRITE(6,*) '*** Version yet without aQG3(x) ***'
*
      INLO = 2
      IQ2  = 2 
      CALL PARAM(INLO,IQ2)
      IF(IQ2.EQ.1)  Q2 = 30.0D0
      IF(IQ2.EQ.2)  Q2 = 100.0D0
      IF(IQ2.EQ.3)  Q2 = 10000.0D0
*
      IF(IQ2.EQ.1.AND.INLO.EQ.2)  alphas = 0.1928
      IF(IQ2.EQ.2.AND.INLO.EQ.2)  alphas = 0.1673
      IF(IQ2.EQ.3.AND.INLO.EQ.2)  alphas = 0.1118
      XMC=1.59D0
      XMC2=XMC**2
      PI = 3.1415926535897932385D0
      nf = 3.0D0
*      as = alphas/(4.0D0*PI)
*     LL = LOG(XMC2/Q2)
*
      CALL INIT 
      ! GPS: overwrite as, LL and nf
      as = a4pi_in
      LL = LL_in
      nf=nflightdp_in
      WRITE(6,*) '>>> nf,as,LL=',nf,as,LL
*
*
      z  = z_in
      x  = 0.1D0
*
      T1 = QGL(z,nf,as,LL) 
      T2 = PSL(z,nf,as,LL)
      T3 = GQ(z,nf,as,LL) 
      T4 = QG(z,nf,as,LL) 
      T5a= GGDEL(z,nf,as,LL) 
      T5b= GGPLU(z,nf,as,LL) 
      T5c= GGREG(z,nf,as,LL) 
      T6 = PS(z,nf,as,LL) 
      T7a= NSDEL(z,nf,as,LL) 
      T7b= NSPLU(z,nf,as,LL) 
      T7c= NSREG(z,nf,as,LL) 
      T8a= AGGDEL(z,nf,as,LL) 
      T8b= AGGPLU(z,nf,as,LL) 
      T8c= AGGREG0(z,nf,as,LL) 
      T8d= AGGREG1(z,nf,as,LL) 
      T9a= APS1(z,nf,as,LL) 
      T9b= APS2(z,nf,as,LL) 
      T10a= ANSDEL(z,nf,as,LL) 
      T10b1= ANSPLU1(z,nf,as,LL) 
      T10b2= ANSPLU2(z,nf,as,LL) 
      T10c= ANSREG(z,nf,as,LL) 
      T15= AGQ(z,nf,as,LL) 
      T14= AQG(z,nf,as,LL) 
      z=1.0D-1
      T16=RED0(z)
      WRITE(6,*) 'z,RED0=',z,T16
      z=9.0D-1
      T16=RED1(z)
      WRITE(6,*) 'z,RED1=',z,T16
      z=1.0D-1
      T16=aQg3(z)
      WRITE(6,*) 'z,aQg3=',z,T16
      WRITE(6,*) '******************************'
      z=1.0D-1
      T17=T16/2.0d0 + RED0(z)  
      WRITE(6,*)'z,AQG3NF0=', z,T17
      T18=AQG3NF(z)
      WRITE(6,*)'z,AQG3NF=', z,T18
      WRITE(6,*)'z,AQG3NF3=', z,T17+3.0D0*T18
      z=1.0D-5
      T16=aQg3(z)
      WRITE(6,*) 'z,T16=',z,T16 
      z=1.0D-4
      T16=aQg3(z)
      WRITE(6,*) 'z,T16=',z,T16 
      z=1.0D-3
      T16=aQg3(z)
      WRITE(6,*) 'z,T16=',z,T16 
      z=1.0D-2
      T16=aQg3(z)
      WRITE(6,*) 'z,T16=',z,T16 
      z=1.0D-1
      T16=aQg3(z)
      WRITE(6,*) 'z,T16=',z,T16 
      z=5.0D-1
      T16=aQg3(z)
      WRITE(6,*) 'z,T16=',z,T16 
      z=0.9D0
      T16=aQg3(z)
      z=9.90D-1
      T16=aQg3(z)
      z=9.990D-1
      T16=aQg3(z)
      WRITE(6,*) 'z,T16=',z,T16 
      z=0.1D0
      WRITE(6,*) '----------------------'    
      T16=AQG3NF(z)
      WRITE(6,*) 
     &'***********************************************************'
      WRITE(6,*) 'Current implementation: NNLO only: INNLO = 2.'
      WRITE(6,*) 'INLO,IQ2=',INLO,IQ2
      WRITE(6,*) 
     &'***********************************************************'
      WRITE(6,*) 'Parameters:'
      WRITE(6,*) 'XMC,Q2, alphas=',XMC,Q2, alphas
      WRITE(6,*) 'z,nf,as,LL =',z,nf,as,LL
      WRITE(6,*) 
     &'***********************************************************'
*      
      T11=NSDIS(z) 
      T12=SIGDIS(z) 
      T13=GLUDIS(z)
*
      WRITE(6,*) 'z, NSDIS=', z,T11        
      WRITE(6,*) 'z, NSDIS=', z,T11        
      WRITE(6,*) 'z, SIGDIS=',z,T12        
      WRITE(6,*) 'z, GLUDIS=',z,T13        
      WRITE(6,*) 
     &'***********************************************************'
*
      WRITE(6,*) 'QGL: complete'
      WRITE(6,*) 'QGL=',T1
      WRITE(6,*) 'PSL: complete'
      WRITE(6,*) 'PSL=',T2
      WRITE(6,*) 'GQ: without agQ3'
      WRITE(6,*) 'GQ=', T3
      WRITE(6,*) 'QG: without aQg3'
      WRITE(6,*) 'QG=', T4
      WRITE(6,*) 'PS: without aQqPS3'
      WRITE(6,*) 'PS=', T6
      WRITE(6,*) 'GG: without agg3'
      WRITE(6,*) 'GGdel=', T5a
      WRITE(6,*) 'GGplu=', T5b
      WRITE(6,*) 'GGreg=', T5c
      WRITE(6,*) 'NS: without aqqQNS3'
      WRITE(6,*) 'NSdel=', T7a
      WRITE(6,*) 'NSplu=', T7b
      WRITE(6,*) 'NSreg=', T7c
      WRITE(6,*) 
     &'***********************************************************'
      WRITE(6,*) 'Known a_ij^(3) terms:'
      WRITE(6,*) 
     &'***********************************************************'
      WRITE(6,*) 'aGGdel=', T8a
      WRITE(6,*) 'aGGplu=', T8b
      WRITE(6,*) 'If z < 1/2'
      WRITE(6,*) 'aGGreg0=', T8c
      WRITE(6,*) 'If z > 1/2'
      WRITE(6,*) 'aGGreg1=', T8d
      WRITE(6,*) 'aPS1=', T9a
      WRITE(6,*) 'aPS2=', T9b
      WRITE(6,*) 'aNSdel=', T10a
      WRITE(6,*) 'aNSplu1=', z,T10b1
      WRITE(6,*) 'aNSplu2=', z,T10b2
      WRITE(6,*) 'aNSreg=', T10c
      WRITE(6,*) 'agQ=', T15
      WRITE(6,*) 'full AQg(nf=3) =', T14
      N=2
      T11=MOMQG(N,x,nf,as,LL)
      WRITE(6,*) 'N,MOMQG=',N,T11
      N=4
      T11=MOMQG(N,x,nf,as,LL)
      WRITE(6,*) 'N,MOMQG=',N,T11
      N=6
      T11=MOMQG(N,x,nf,as,LL)
      WRITE(6,*) 'N,MOMQG=',N,T11
      N=8
      T11=MOMQG(N,x,nf,as,LL)
      WRITE(6,*) 'N,MOMQG=',N,T11
      N=10
      T11=MOMQG(N,x,nf,as,LL)
      WRITE(6,*) 'N,MOMQG=',N,T11
      WRITE(6,*) 
     &'***********************************************************'
*      LL=10.0D0
*      WRITE(6,*) '>>> LL=',LL
      WRITE(6,*) 'Moment calculation:'
      WRITE(6,*) '-------------------'
      WRITE(6,*) 'Moments of AqqQPS'
      WRITE(6,*) '-----------------'
      DO 100 K=1,5
      N=2*K
      U1=MOMPSL(N,z,nf,as,LL)
      IF(K.EQ.1) REF=MOMPSLF(1,0,0,0,0)
      IF(K.EQ.2) REF=MOMPSLF(0,1,0,0,0)
      IF(K.EQ.3) REF=MOMPSLF(0,0,1,0,0)
      IF(K.EQ.4) REF=MOMPSLF(0,0,0,1,0)
      IF(K.EQ.5) REF=MOMPSLF(0,0,0,0,1)
      WRITE(6,*) 'N,MOMPSL,RAT=',N,U1,U1/REF-1.0D0
100   CONTINUE
      WRITE(6,*) 
     &'-----------------------------------------------------------'
      WRITE(6,*) 'Moments of AqgQ'
      WRITE(6,*) '---------------'
      DO 101 K=1,5
      N=2*K
      U2=MOMQGL(N,z,nf,as,LL)
      IF(K.EQ.1) REF=MOMQGLF(1,0,0,0,0)
      IF(K.EQ.2) REF=MOMQGLF(0,1,0,0,0)
      IF(K.EQ.3) REF=MOMQGLF(0,0,1,0,0)
      IF(K.EQ.4) REF=MOMQGLF(0,0,0,1,0)
      IF(K.EQ.5) REF=MOMQGLF(0,0,0,0,1)
      WRITE(6,*) 'N,MOMQGLF,RAT=',N,U2,U2/REF-1.0D0
101   CONTINUE
      WRITE(6,*) 
     &'-----------------------------------------------------------'
      WRITE(6,*) 'Moments of AgqQ'
      WRITE(6,*) '---------------'
      DO 103 K=1,5
      N=2*K
      U3=MOMGQ(N,z,nf,as,LL)
      U4=MOMAGQ(N,z,nf,as,LL)
      IF(K.EQ.1) REF=MOMGQF(1,0,0,0,0)
      IF(K.EQ.2) REF=MOMGQF(0,1,0,0,0)
      IF(K.EQ.3) REF=MOMGQF(0,0,1,0,0)
      IF(K.EQ.4) REF=MOMGQF(0,0,0,1,0)
      IF(K.EQ.5) REF=MOMGQF(0,0,0,0,1)
      WRITE(6,*) 'N,MOMAGQ,SUM=',N,U3+U4,(U3+U4)/REF-1.0D0
103   CONTINUE
      WRITE(6,*) 
     &'-----------------------------------------------------------'
      WRITE(6,*) 'Moments of AQqPS'
      WRITE(6,*) '----------------'
      DO 111 K=1,5
      N=2*K
      U7=MOMPS(N,z,nf,as,LL)
      U5=MOMAPS(N,z,nf,as,LL)
      IF(K.EQ.1) REF=MOMPSF(1,0,0,0,0)
      IF(K.EQ.2) REF=MOMPSF(0,1,0,0,0)
      IF(K.EQ.3) REF=MOMPSF(0,0,1,0,0)
      IF(K.EQ.4) REF=MOMPSF(0,0,0,1,0)
      IF(K.EQ.5) REF=MOMPSF(0,0,0,0,1)
      WRITE(6,*) 'N,MOMPS,RAT=',N,U7+U5,(U7+U5)/REF-1.0D0
111   CONTINUE
      WRITE(6,*) 
     &'-----------------------------------------------------------'
      WRITE(6,*) 'Moments of AggQ'
      WRITE(6,*) '---------------'
      DO 112 K=1,5
      N=2*K
      U7=MOMGG(N,z,nf,as,LL)
      IF(K.EQ.1) REF=MOMGGF(1,0,0,0,0)
      IF(K.EQ.2) REF=MOMGGF(0,1,0,0,0)
      IF(K.EQ.3) REF=MOMGGF(0,0,1,0,0)
      IF(K.EQ.4) REF=MOMGGF(0,0,0,1,0)
      IF(K.EQ.5) REF=MOMGGF(0,0,0,0,1)
      WRITE(6,*) 'N,MOMGG,RAT=',N,U7,U7/REF-1.0D0
      WRITE(6,*) 'N,MOMGGF=',N,REF
112   CONTINUE
      WRITE(6,*) 
     &'-----------------------------------------------------------'
      WRITE(6,*) 'Moments of AQg'
      WRITE(6,*) '--------------'
*      WRITE(6,*) '>>> not complete in this version <<<'
      DO 104 K=1,5
      N=2*K
      U7=MOMQG(N,z,nf,as,LL)
      IF(K.EQ.1) REF=MOMQGF(1,0,0,0,0)
      IF(K.EQ.2) REF=MOMQGF(0,1,0,0,0)
      IF(K.EQ.3) REF=MOMQGF(0,0,1,0,0)
      IF(K.EQ.4) REF=MOMQGF(0,0,0,1,0)
      IF(K.EQ.5) REF=MOMQGF(0,0,0,0,1)
      WRITE(6,*) 'N,MOMQG,RAT=',N,U7,U7/REF-1.0D0
104   CONTINUE
      WRITE(6,*) 
     &'-----------------------------------------------------------'
      WRITE(6,*) 'Moments of AqqQ^NS'
      WRITE(6,*) '------------------'
      DO 105 K=1,5
      N=2*K
      U7=MOMNS(N,z,nf,as,LL) 
      IF(K.EQ.1) REF=MOMQQNS(1,0,0,0,0)
      IF(K.EQ.2) REF=MOMQQNS(0,1,0,0,0)
      IF(K.EQ.3) REF=MOMQQNS(0,0,1,0,0)
      IF(K.EQ.4) REF=MOMQQNS(0,0,0,1,0)
      IF(K.EQ.5) REF=MOMQQNS(0,0,0,0,1)
      WRITE(6,*) 'N,MOMNS,RAT=',N,U7,U7/REF-1.0D0
105   CONTINUE
      WRITE(6,*) 
     &'-----------------------------------------------------------'
      END
      SUBROUTINE INIT
*     ---------------
*     ------------------------------------------------------------------------
*     TEST case: initialization
*     Code: J. Bluemlein, March 19, 2024
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 nf,as,LL
      COMMON/VAR1/ nf,as,LL
*
      nf = 3.0D0
      as = 0.1D0
      LL = 1.0D1
*
      RETURN
      END
      REAL*8 FUNCTION DIST(z)
*     --------------------------
*     ------------------------------------------------------------------------
*     TESt distributoon for Mellin convolutions
*     Code: J. Bluemlein, March 15, 2024
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 z,x,nf,as,LL
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT
*
      DIST=3.0D0*(1.0D0-z)**5
*
      RETURN
      END
      REAL*8 FUNCTION AQG3NF(z)
*     -------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 11, 2024
*     regular part of aqqQNS3: NF contribution
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
*
      REAL*8 w(46)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4,x
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
*
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      x=z
      w(1)=x**(-1)
      w(2)=Hr1(-1)
      w(3)=Hr1(0)
      w(4)=Hr2(0,-1)
      w(5)=Hr2(0,1)
      w(6)=Hr3(0,-1,-1)
      w(7)=Hr3(0,0,-1)
      w(8)=Hr3(0,0,1)
      w(9)=Hr1(1)
      w(10)=Hr3(0,-1,1)
      w(11)=Hr4(0,0,0,1)
      w(12)=Hr3(0,1,-1)
      w(13)=Hr3(0,1,1)
      w(14)=Hr4(0,-1,-1,-1)
      w(15)=Hr4(0,-1,0,1)
      w(16)=Hr4(0,0,-1,-1)
      w(17)=Hr4(0,0,-1,1)
      w(18)=Hr4(0,0,0,-1)
      w(19)=Hr5(0,0,0,0,1)
      w(20)=Hr4(0,0,1,-1)
      w(21)=Hr4(0,0,1,1)
      w(22)=Hr4(0,1,1,1)
      w(23)=2*x
      w(24)=w(23) - 1
      w(25)=4*w(3)
      w(25)=w(24)*w(25)
      w(26)=28*x
      w(27)=17 + w(26)
      w(27)=w(27)*w(23)
      w(27)= - 35 + w(27)
      w(27)=1.D0/3.D0*w(27) + w(25)
      w(28)=1.D0/3.D0*w(3)
      w(27)=w(27)*w(28)
      w(29)=x + 1
      w(29)=w(29)*w(23)
      w(30)=w(29) + 1
      w(31)=8*w(2)
      w(31)=w(30)*w(31)
      w(32)= - 260.D0/3.D0 - 307*x
      w(32)=w(32)*w(23)
      w(32)= - 709.D0/3.D0 + w(32)
      w(33)=x - 1
      w(33)=w(33)*w(23)
      w(34)=w(33) + 1
      w(35)=4.D0/9.D0*w(9)
      w(36)=w(34)*w(35)
      w(27)=w(27) + w(36) + 1.D0/9.D0*w(32) + w(31)
      w(32)=2*w(3)
      w(27)=w(27)*w(32)
      w(36)=3*w(2)
      w(37)=w(30)*w(36)
      w(29)=w(29) - 5
      w(38)= - w(37) + 2.D0/3.D0*w(29)
      w(38)=w(38)*w(2)
      w(39)=w(34)*w(9)
      w(33)=w(33) - 5
      w(40)=11.D0/3.D0*w(33) + 4*w(39)
      w(35)=w(40)*w(35)
      w(40)=w(30)*w(4)
      w(41)=628.D0/81.D0*w(1)
      w(42)= - 14738 + 41597*x
      w(42)=x*w(42)
      w(42)= - 4796 + w(42)
      w(43)=1427 + 148*x
      w(43)=x*w(43)
      w(43)=182 + w(43)
      w(43)=z2*w(43)
      w(27)=16.D0/135.D0*w(43) - 16*w(40) + w(27) + w(35) + 1.D0/243.D0
     & *w(42) - 4*w(38) - w(41)
      w(27)=z2*w(27)
      w(35)= - 439 - 65*x
      w(35)=w(35)*w(23)
      w(35)=277 + w(35)
      w(35)=w(3)*w(35)
      w(42)=136.D0/9.D0*w(1)
      w(43)=15104 + 193*x
      w(43)=x*w(43)
      w(43)= - 211 + w(43)
      w(31)=1.D0/27.D0*w(35) + w(42) - 110.D0/9.D0*w(39) + 1.D0/81.D0*
     & w(43) - w(31)
      w(31)=w(3)*w(31)
      w(35)= - 4*w(33) - 5*w(39)
      w(35)=w(9)*w(35)
      w(43)= - 15305 + 31733*x
      w(43)=x*w(43)
      w(43)= - 12203 + w(43)
      w(35)= - 628.D0/9.D0*w(1) + 1.D0/27.D0*w(43) + w(35)
      w(43)= - 5.D0/3.D0*w(5) - 14.D0/3.D0*z2
      w(43)=w(34)*w(43)
      w(31)=1.D0/9.D0*w(35) + w(31) + w(43)
      w(31)=w(5)*w(31)
      w(35)=w(30)*w(2)
      w(43)=17 + x
      w(43)=x*w(43)
      w(43)= - 265 + 82*w(43)
      w(43)=1.D0/27.D0*w(43) - w(25)
      w(32)=w(43)*w(32)
      w(43)= - 17996 - 8233*x
      w(43)=x*w(43)
      w(43)= - 4283 + w(43)
      w(32)=w(32) - w(42) + 172.D0/9.D0*w(39) + 1.D0/81.D0*w(43) + 16*
     & w(35)
      w(32)=w(8)*w(32)
      w(43)= - 23.D0/9.D0 - 7*x
      w(43)=x*w(43)
      w(44)=w(3)*w(24)
      w(43)= - 58.D0/3.D0*w(44) + 371.D0/9.D0 + 8*w(43)
      w(43)=w(3)*w(43)
      w(44)=10493 + 10970*x
      w(44)=x*w(44)
      w(44)=6575 + w(44)
      w(43)=2.D0/3.D0*w(43) + 476.D0/81.D0*w(1) - 262.D0/27.D0*w(39) + 
     & 1.D0/81.D0*w(44) - 10*w(35)
      w(43)=z3*w(43)
      w(31)=w(32) + w(43) + w(31)
      w(32)= - 2*w(33) - 7*w(39)
      w(32)=w(9)*w(32)
      w(33)=20731 - 33541*x
      w(33)=x*w(33)
      w(33)=10585 + w(33)
      w(32)=1.D0/27.D0*w(33) + w(32)
      w(32)=1.D0/9.D0*w(32) + w(41)
      w(32)=w(9)*w(32)
      w(33)=1.D0/3.D0*w(29)
      w(41)= - w(33) + w(35)
      w(41)=w(2)*w(41)
      w(43)=397 + 388*x
      w(43)=w(43)*x
      w(43)=w(43) + 56
      w(43)=1.D0/27.D0*w(43)
      w(41)=w(43) + w(41)
      w(41)=w(2)*w(41)
      w(44)=8*w(10) - 20.D0/3.D0*w(7)
      w(44)=w(30)*w(44)
      w(45)= - 76969.D0/3.D0 - 133924*x
      w(45)=x*w(45)
      w(45)= - 110314.D0/3.D0 + w(45)
      w(46)=w(11)*w(24)
      w(32)=32*w(46) + 1.D0/243.D0*w(45) + 2*w(41) + w(44) + w(32)
      w(41)=125.D0/81.D0 - 77*x
      w(41)=x*w(41)
      w(41)= - 158.D0/27.D0*w(39) - 20.D0/3.D0*w(35) - 3745.D0/81.D0 + 
     & w(41)
      w(44)=77 + 124*x
      w(44)=w(44)*w(23)
      w(44)= - 203 + w(44)
      w(25)=1.D0/9.D0*w(44) + w(25)
      w(25)=w(25)*w(28)
      w(25)=2*w(41) + w(25)
      w(25)=w(25)*w(28)
      w(28)=2.D0/9.D0*w(29) - w(35)
      w(28)=w(2)*w(28)
      w(29)= - 11996 + 8333*x
      w(29)=x*w(29)
      w(29)=4705 + w(29)
      w(29)=1.D0/27.D0*w(29) + 10*w(39)
      w(29)= - w(42) + 1.D0/3.D0*w(29)
      w(29)=w(9)*w(29)
      w(41)= - 14917 + 60043*x
      w(41)=x*w(41)
      w(41)= - 13342 + w(41)
      w(28)=1.D0/243.D0*w(41) + w(28) + w(29)
      w(25)=2*w(28) + w(25)
      w(25)=w(3)*w(25)
      w(28)=w(12)*w(30)
      w(25)=32*w(28) + 4*w(32) + w(25)
      w(25)=w(3)*w(25)
      w(28)=5.D0/3.D0*w(3) + w(36)
      w(28)=w(30)*w(28)
      w(26)=5 - w(26)
      w(26)=w(26)*w(23)
      w(26)=5 + w(26)
      w(26)=1.D0/9.D0*w(26) + w(28)
      w(26)=w(3)*w(26)
      w(26)=w(40) + w(26) - w(43) + w(38)
      w(26)=w(4)*w(26)
      w(28)= - w(33) + w(37)
      w(29)=5*w(30)
      w(29)= - w(3)*w(29)
      w(28)=2*w(28) + w(29)
      w(28)=w(6)*w(28)
      w(29)= - 1 + 5*x
      w(29)=x*w(29)
      w(29)= - 5 + 22*w(29)
      w(29)=1.D0/9.D0*w(29) - 5*w(35)
      w(29)=w(7)*w(29)
      w(26)=w(28) + w(29) + w(26)
      w(28)=8906 - 10157*x
      w(28)=w(28)*w(23)
      w(28)=5677 + w(28)
      w(29)= - 80 + 107*x
      w(29)=x*w(29)
      w(29)= - 227 + w(29)
      w(29)=1.D0/3.D0*w(29) + 20*w(39)
      w(29)=w(9)*w(29)
      w(32)=3221 - 3248*x
      w(32)=x*w(32)
      w(32)=65 + w(32)
      w(29)=1.D0/3.D0*w(32) + w(29)
      w(29)=w(9)*w(29)
      w(28)=2.D0/9.D0*w(28) + w(29)
      w(28)=w(9)*w(28)
      w(29)=104.D0/27.D0*w(22) - 896.D0/9.D0*w(21)
      w(29)=w(34)*w(29)
      w(32)=w(20) + w(17)
      w(32)= - 32*w(15) + 40*w(16) + 80.D0/3.D0*w(18) - 48*w(14) - 64*
     & w(32)
      w(30)=w(30)*w(32)
      w(32)=z5 - w(19)
      w(32)=128*w(32)
      w(24)=w(24)*w(32)
      w(32)= - 658337 + 2587252.D0/3.D0*x
      w(32)=w(32)*w(23)
      w(32)= - 381527 + w(32)
      w(33)= - 853 + 151*x
      w(23)=w(33)*w(23)
      w(23)=43 + w(23)
      w(23)=w(11)*w(23)
      w(33)=160*w(3) + 14*w(9)
      w(33)=w(34)*w(33)
      w(34)=739 + 665*x
      w(34)=x*w(34)
      w(34)= - 164 + w(34)
      w(33)=1.D0/9.D0*w(34) + w(33)
      w(33)=w(13)*w(33)
      w(34)= - 28351.D0/243.D0 + 2*w(9)
      w(34)=w(1)*w(34)
*
      AQG3NF = 8.D0/27.D0*w(23) + w(24) + w(25) + 8*w(26) + w(27) 
     &+ 2.D0/81.D 0*w(28) + w(29) + w(30) + 4*w(31) + 1.D0/729.D0*w(32) 
     &+ 4.D0/9.D0*w(33) + 8.D0/9.D0*w(34)
*
      RETURN
      END      
      REAL*8 FUNCTION RED0(z)
*     -----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 15, 2024
*     NF = 0 one particle reducible contribution to a_Qg^(3) for z \in [0,1/2]
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(23)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
*
      w(1)=z**(-1)
      w(2)=Log(z)
      w(3)=1.D0/49.D0*z
      w(4)=9129761.D0/752.D0 - 990319763.D0/31875.D0*z
      w(4)=w(4)*w(3)
      w(4)= - 6451613683.D0/9754944.D0 + w(4)
      w(4)=w(4)*w(3)
      w(4)=55495901.D0/10517049.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 289313603.D0/20511975.D0 + w(4)
      w(5)=1.D0/2.D0*z
      w(4)=w(4)*w(5)
      w(4)=222817321.D0/80776575.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 4163321497.D0/563797080.D0 + w(4)
      w(4)=z*w(4)
      w(4)=118081219.D0/40861051.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 823399957.D0/106256010.D0 + w(4)
      w(4)=z*w(4)
      w(4)=39783683.D0/13111800.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 858540853.D0/105320800.D0 + w(4)
      w(4)=w(4)*w(5)
      w(4)=831793.D0/520923.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 994309529.D0/231289812.D0 + w(4)
      w(4)=z*w(4)
      w(4)=920477.D0/546231.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 1513805509.D0/332866800.D0 + w(4)
      w(4)=z*w(4)
      w(4)=125032067.D0/70096950.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 23630939.D0/4895660.D0 + w(4)
      w(4)=z*w(4)
      w(4)=8433517.D0/4452096.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 834948043.D0/162370560.D0 + w(4)
      w(4)=z*w(4)
      w(4)=11254441.D0/5573800.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 301186007.D0/54749100.D0 + w(4)
      w(4)=z*w(4)
      w(4)=4809431.D0/2225286.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 47191337.D0/7980336.D0 + w(4)
      w(4)=z*w(4)
      w(4)=954511.D0/410670.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 72912691.D0/11407500.D0 + w(4)
      w(4)=z*w(4)
      w(4)=326693.D0/130000.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 194198647.D0/27931200.D0 + w(4)
      w(4)=z*w(4)
      w(4)=3674798.D0/1344189.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 20782343.D0/2727340.D0 + w(4)
      w(4)=z*w(4)
      w(4)=2300534.D0/768075.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 76653469.D0/9097200.D0 + w(4)
      w(4)=z*w(4)
      w(4)=40609.D0/12274.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 5585201.D0/593028.D0 + w(4)
      w(4)=z*w(4)
      w(4)=614167.D0/166464.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 8123873.D0/761600.D0 + w(4)
      w(4)=z*w(4)
      w(4)=794821.D0/191100.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 6102611.D0/496860.D0 + w(4)
      w(4)=z*w(4)
      w(4)=1357.D0/286.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 5452873.D0/377520.D0 + w(4)
      w(4)=z*w(4)
      w(4)=27091.D0/4950.D0 + w(4)
      w(6)=1.D0/3.D0*z
      w(4)=w(4)*w(6)
      w(4)= - 6401.D0/1100.D0 + w(4)
      w(4)=z*w(4)
      w(4)=19129.D0/9072.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 616783.D0/84672.D0 + w(4)
      w(4)=z*w(4)
      w(4)=341.D0/147.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 59537.D0/6300.D0 + w(4)
      w(4)=z*w(4)
      w(4)=323.D0/225.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 2599.D0/240.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 373.D0/18.D0 + w(4)
      w(4)=z*w(4)
      w(4)=509.D0/9.D0 + w(4)
      w(4)=w(4)*w(6)
      w(4)= - 281.D0/2.D0 + w(4)
      w(4)=z*w(4)
      w(7)=5.D0/2.D0 - 227.D0/9.D0*z
      w(7)=z*w(7)
      w(8)=3*z
      w(9)= - 323.D0/12.D0 - w(8)
      w(9)=w(2)*w(9)
      w(7)=w(9) + 3*w(1) + 511.D0/6.D0 + w(7)
      w(7)=w(2)*w(7)
      w(4)=1.D0/2.D0*w(7) - 47.D0/4.D0*w(1) - 1181.D0/72.D0 + w(4)
      w(4)=w(2)*w(4)
      w(7)=1.D0/43.D0*z
      w(9)=88460832466300759836991338472039100104920733.D0/6251.D0 + 
     & 1386343784310605469304173389032388616066469.D0/50.D0*z
      w(9)=z*w(9)
      w(9)=1275493636881995615401290024906767818617496933.D0/30890656.D0
     &  + 1.D0/175.D0*w(9)
      w(9)=z*w(9)
      w(9)=24468684672292036464800170528467507104359.D0/111013295.D0 + 
     & 1.D0/95697.D0*w(9)
      w(9)=z*w(9)
      w(9)=437450874070998231292559675049338411170939.D0/108733979475.D0
     &  + 1.D0/28.D0*w(9)
      w(9)=z*w(9)
      w(9)=6727782258957630142503886076054805871429.D0/1998250912350.D0
     &  + 1.D0/611.D0*w(9)
      w(9)=z*w(9)
      w(9)=156150794814724273909288666340236422387319.D0/90656879072760.
     & D0 + w(9)
      w(9)=z*w(9)
      w(9)=2172283540828910385075304011576160063189.D0/19711003252941.D0
     &  + 1.D0/8.D0*w(9)
      w(9)=z*w(9)
      w(9)=5143233374362244918630209616135228070683.D0/91123454079840.D0
     &  + w(9)
      w(9)=w(9)*w(7)
      w(9)=405194913785782181227348319599061153.D0/602382315600.D0 + 
     & w(9)
      w(10)=1.D0/11.D0*z
      w(9)=w(9)*w(10)
      w(9)=12744574392379059423478367308096300961.D0/406446448262400.D0
     &  + w(9)
      w(11)=1.D0/41.D0*z
      w(9)=w(9)*w(11)
      w(9)=376028841808034803600250758218759937.D0/958247068418640.D0
     &  + w(9)
      w(9)=w(9)*w(5)
      w(9)=312226606155393195044154963863287.D0/3099225658347.D0 + w(9)
      w(9)=z*w(9)
      w(9)=24006102816012623325533760779412311.D0/463754576244960.D0 + 
     & w(9)
      w(9)=z*w(9)
      w(9)=125375449592529062111798486608140701.D0/4710110489884800.D0
     &  + w(9)
      w(9)=z*w(9)
      w(9)=1380562171521135220749716146987.D0/3728877352200.D0 + 1.D0/
     & 37.D0*w(9)
      w(9)=z*w(9)
      w(9)=317051712051778818055174077983777.D0/1662582115434240.D0 + 
     & w(9)
      w(9)=z*w(9)
      w(9)=7569809929363781836388962717853.D0/76997266970624.D0 + w(9)
      w(9)=z*w(9)
      w(9)=29140602371572895812032704575387.D0/574391376599040.D0 + 
     & w(9)
      w(9)=z*w(9)
      w(9)=6200835575446750458724151183287.D0/236610084110400.D0 + w(9)
      w(9)=z*w(9)
      w(9)=2624800164411491998248264603311.D0/193676802719400.D0 + w(9)
      w(9)=z*w(9)
      w(9)=439944889746674625886555063.D0/62700317680.D0 + w(9)
      w(9)=z*w(9)
      w(9)=245363363613158326258218527.D0/67456893504.D0 + w(9)
      w(9)=z*w(9)
      w(9)=53394368936039294550260761.D0/820105936650.D0 + 1.D0/29.D0*
     & w(9)
      w(9)=z*w(9)
      w(9)=176232632370950299169096347.D0/5207021820000.D0 + w(9)
      w(9)=z*w(9)
      w(9)=2004186805569269340888553.D0/113733620000.D0 + w(9)
      w(9)=w(9)*w(6)
      w(9)=2169887827125789901456847.D0/708297990400.D0 + w(9)
      w(12)=1.D0/5.D0*z
      w(9)=w(9)*w(12)
      w(9)=2297410398434641314683.D0/7176177008.D0 + w(9)
      w(9)=w(9)*w(6)
      w(9)=2898208959213060971671.D0/51871279460.D0 + w(9)
      w(13)=1.D0/23.D0*z
      w(9)=w(9)*w(13)
      w(9)=893792831736206946301.D0/701184884400.D0 + w(9)
      w(14)=1.D0/17.D0*z
      w(9)=w(9)*w(14)
      w(9)=81892894404101737153.D0/2076235761600.D0 + w(9)
      w(9)=z*w(9)
      w(9)=7429809751364447813.D0/357161985180.D0 + w(9)
      w(9)=z*w(9)
      w(9)=27627156824025367117.D0/2510045568576.D0 + w(9)
      w(9)=z*w(9)
      w(9)=90250968973507169.D0/15451188480.D0 + w(9)
      w(9)=z*w(9)
      w(9)=1089680408637648643.D0/349923974400.D0 + w(9)
      w(9)=z*w(9)
      w(9)=115814433566393.D0/69560400.D0 + w(9)
      w(9)=z*w(9)
      w(9)=4016471773006753.D0/4476211740.D0 + w(9)
      w(9)=z*w(9)
      w(9)=1361067575001263.D0/2813618808.D0 + w(9)
      w(9)=z*w(9)
      w(9)=1805857773150037.D0/6802155360.D0 + w(9)
      w(15)=1.D0/16.D0*z
      w(9)=w(9)*w(15)
      w(9)=6290662261.D0/698775.D0 + w(9)
      w(9)=z*w(9)
      w(9)=2253903889057.D0/439084800.D0 + w(9)
      w(9)=z*w(9)
      w(9)=84411318161.D0/30481920.D0 + w(9)
      w(9)=z*w(9)
      w(9)=3536428759.D0/2032128.D0 + w(9)
      w(9)=z*w(9)
      w(9)=297682339.D0/329280.D0 + w(9)
      w(9)=z*w(9)
      w(9)=63778837.D0/88200.D0 + w(9)
      w(9)=z*w(9)
      w(9)=887407.D0/1800.D0 + w(9)
      w(9)=z*w(9)
      w(9)=163019.D0/2880.D0 + w(9)
      w(9)=w(9)*w(12)
      w(9)=3215.D0/24.D0 + w(9)
      w(9)=z*w(9)
      w(9)=75367.D0/144.D0 + w(9)
      w(9)=w(9)*w(6)
      w(9)=3701.D0/8.D0 + w(9)
      w(9)=z*w(9)
      w(9)=98681.D0/144.D0 + w(9)
      w(16)= - 2099 - 6233.D0/15.D0*z
      w(16)=z*w(16)
      w(17)= - 1777.D0/3.D0 + 66*z
      w(17)=w(2)*w(17)
      w(16)=1.D0/5.D0*w(17) + 2648.D0/15.D0*w(1) + 356.D0/5.D0 + 1.D0/
     & 12.D0*w(16)
      w(16)=z2*w(16)
      w(17)=2159.D0/12.D0 - 30*z
      w(17)=z3*w(17)
      w(4)=w(16) + w(4) + 1337.D0/24.D0*w(1) + 1.D0/3.D0*w(9) + w(17)
      w(4)=z2*w(4)
      w(9)=1.D0/1225.D0*z
      w(16)=58973888310883717493969307958141065703365379.D0/752.D0 + 
     & 195907706020375230740567263186195211581296457.D0/1275.D0*z
      w(16)=w(16)*w(9)
      w(16)=159436704610193865217700678269803785533430611.D0/4877472.D0
     &  + w(16)
      w(16)=z*w(16)
      w(16)=160002258365087305104371441259218825662831.D0/52585245.D0
     &  + 1.D0/5488.D0*w(16)
      w(17)=1.D0/7.D0*z
      w(16)=w(16)*w(17)
      w(16)=36454239505866075836440771876408721538061.D0/164095800.D0
     &  + w(16)
      w(16)=z*w(16)
      w(16)=874611693668333445977381213973479581393.D0/361879056.D0 + 1.
     & D0/47.D0*w(16)
      w(16)=w(16)*w(12)
      w(16)=390376987034667324392274412535570344729.D0/1578631824.D0 + 
     & w(16)
      w(16)=z*w(16)
      w(16)=724094513622209486200780428929679930047.D0/5720547140.D0 + 
     & w(16)
      w(16)=z*w(16)
      w(16)=1928712515343836123475687855457698055627.D0/29751682800.D0
     &  + w(16)
      w(7)=w(16)*w(7)
      w(7)=2836364396695023750629815645885685017.D0/3671304000.D0 + 
     & w(7)
      w(7)=z*w(7)
      w(7)=312367019394077515915598865259177.D0/147449120.D0 + 1.D0/187.
     & D0*w(7)
      w(7)=w(7)*w(11)
      w(7)=47003605238730031368273378773171341.D0/1772909338200.D0 + 
     & w(7)
      w(7)=z*w(7)
      w(7)=535245610368843630790290891668167.D0/39319268040.D0 + w(7)
      w(7)=z*w(7)
      w(7)=2000508570136953669471989025054023.D0/286006551600.D0 + w(7)
      w(7)=z*w(7)
      w(7)=211782853754065227053477666744279.D0/58881438000.D0 + w(7)
      w(11)=1.D0/703.D0*z
      w(7)=w(7)*w(11)
      w(7)=483196762064799969157207461437.D0/183513815100.D0 + w(7)
      w(7)=z*w(7)
      w(7)=5284195172455441831566650604143.D0/3896318715520.D0 + w(7)
      w(16)=1.D0/13.D0*z
      w(7)=w(7)*w(16)
      w(7)=93582816609690621780316535191.D0/1740012679680.D0 + w(7)
      w(7)=z*w(7)
      w(7)=219102268805147780760110188061.D0/7894456627200.D0 + w(7)
      w(7)=z*w(7)
      w(7)=1033472663122449564307750459609.D0/72085509496000.D0 + w(7)
      w(7)=z*w(7)
      w(7)=131239997155126873104994286639.D0/17701643259300.D0 + w(7)
      w(18)=1.D0/31.D0*z
      w(7)=w(7)*w(18)
      w(7)=1272698041515446595869870749.D0/10278373505400.D0 + w(7)
      w(7)=z*w(7)
      w(7)=144975678532463582458009.D0/2258601345.D0 + w(7)
      w(7)=z*w(7)
      w(7)=386916103544672987854157.D0/7745444957250.D0 + 1.D0/667.D0*
     & w(7)
      w(7)=z*w(7)
      w(7)=11014524969380487058729451.D0/424155319087500.D0 + w(7)
      w(7)=w(7)*w(6)
      w(7)=1002097315145038720769197.D0/222349227100000.D0 + w(7)
      w(19)=1.D0/9.D0*z
      w(7)=w(7)*w(19)
      w(7)=90411519884255557567627.D0/346180642808000.D0 + w(7)
      w(20)=1.D0/4.D0*z
      w(7)=w(7)*w(20)
      w(7)=1023081131550096217397.D0/29987898183243.D0 + w(7)
      w(7)=z*w(7)
      w(7)=143315057060673514697.D0/8023517908560.D0 + w(7)
      w(7)=z*w(7)
      w(7)=223461407626826227939.D0/23840286069600.D0 + w(7)
      w(7)=z*w(7)
      w(7)=70620178102211027.D0/14324678550.D0 + w(7)
      w(7)=z*w(7)
      w(7)=16720889803944896897.D0/6428915733240.D0 + w(7)
      w(7)=z*w(7)
      w(7)=10849991956691576989.D0/7888714644096.D0 + w(7)
      w(7)=z*w(7)
      w(7)=372612046864045169.D0/509889219840.D0 + w(7)
      w(7)=z*w(7)
      w(7)=120918848608956401.D0/311043532800.D0 + w(7)
      w(7)=z*w(7)
      w(7)=5752003296804451.D0/27545918400.D0 + w(7)
      w(7)=z*w(7)
      w(7)=99893106085397.D0/895242348.D0 + w(7)
      w(7)=w(7)*w(12)
      w(7)=22970678044885.D0/1875745872.D0 + w(7)
      w(7)=z*w(7)
      w(7)=44211891677651.D0/6802155360.D0 + w(7)
      w(7)=z*w(7)
      w(7)=711137974409.D0/188669250.D0 + w(7)
      w(7)=z*w(7)
      w(7)=46968958369.D0/24948000.D0 + w(7)
      w(7)=z*w(7)
      w(7)=1368747181.D0/1058400.D0 + w(7)
      w(7)=z*w(7)
      w(7)=16643806763.D0/35562240.D0 + w(7)
      w(7)=z*w(7)
      w(7)=689312651.D0/1234800.D0 + w(7)
      w(7)=z*w(7)
      w(7)= - 8348041.D0/147000.D0 + w(7)
      w(7)=w(7)*w(5)
      w(7)=186587.D0/1125.D0 + w(7)
      w(7)=z*w(7)
      w(7)= - 4379879.D0/14400.D0 + w(7)
      w(7)=z*w(7)
      w(7)= - 51577.D0/72.D0 + w(7)
      w(7)=z*w(7)
      w(7)= - 134167.D0/72.D0 + w(7)
      w(7)=z*w(7)
      w(7)= - 5147 + w(7)
      w(7)=z*w(7)
      w(7)=26459.D0/24.D0 + w(7)
      w(21)=1.D0/25.D0*z
      w(22)= - 5293193.D0/2961.D0 + 135473.D0/200.D0*z
      w(22)=w(22)*w(21)
      w(22)=11731439.D0/415104.D0 + w(22)
      w(22)=w(22)*w(3)
      w(22)= - 960874693.D0/631022940.D0 + w(22)
      w(22)=z*w(22)
      w(22)=22253021.D0/36921555.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 154301449.D0/96931890.D0 + w(22)
      w(22)=z*w(22)
      w(22)=2256701.D0/3575880.D0 + w(22)
      w(22)=w(22)*w(5)
      w(22)= - 306606745.D0/367749459.D0 + w(22)
      w(22)=z*w(22)
      w(22)=10297883.D0/31099320.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 80334461.D0/91782600.D0 + w(22)
      w(22)=z*w(22)
      w(22)=3907343.D0/11217600.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 124752169.D0/135439980.D0 + w(22)
      w(22)=w(22)*w(5)
      w(22)=23276.D0/126711.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 95557201.D0/196643160.D0 + w(22)
      w(22)=z*w(22)
      w(22)=11088563.D0/57062880.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 140444.D0/273105.D0 + w(22)
      w(22)=z*w(22)
      w(22)=777527.D0/3769920.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 40107715.D0/73459584.D0 + w(22)
      w(22)=z*w(22)
      w(22)=209237.D0/952320.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 116740321.D0/200656800.D0 + w(22)
      w(22)=z*w(22)
      w(22)=147913.D0/629300.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 83112457.D0/133517160.D0 + w(22)
      w(22)=z*w(22)
      w(22)=12097301.D0/47882016.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 28851919.D0/43120350.D0 + w(22)
      w(22)=z*w(22)
      w(22)=1495199.D0/5475600.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 778213.D0/1076400.D0 + w(22)
      w(22)=z*w(22)
      w(22)=2164201.D0/7286400.D0 + w(22)
      w(22)=w(22)*w(6)
      w(22)= - 938965.D0/3584504.D0 + w(22)
      w(22)=z*w(22)
      w(22)=51881.D0/478170.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 754931.D0/2633400.D0 + w(22)
      w(22)=z*w(22)
      w(22)=1033883.D0/8618400.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 785861.D0/2485485.D0 + w(22)
      w(22)=z*w(22)
      w(22)=336721.D0/2511648.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 658831.D0/1872720.D0 + w(22)
      w(22)=z*w(22)
      w(22)=73447.D0/483840.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 543833.D0/1375920.D0 + w(22)
      w(22)=z*w(22)
      w(22)=1541.D0/8820.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 114575.D0/255528.D0 + w(22)
      w(22)=z*w(22)
      w(22)=42277.D0/205920.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 250853.D0/490050.D0 + w(22)
      w(22)=z*w(22)
      w(22)=8027.D0/32400.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 78691.D0/136080.D0 + w(22)
      w(22)=z*w(22)
      w(22)=22417.D0/72576.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 4469.D0/7560.D0 + w(22)
      w(22)=z*w(22)
      w(22)=373.D0/945.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 37.D0/360.D0 + w(22)
      w(22)=z*w(22)
      w(22)=1673.D0/4320.D0 + w(22)
      w(22)=z*w(22)
      w(22)=10 + w(22)
      w(22)=z*w(22)
      w(22)=19871.D0/216.D0 + w(22)
      w(22)=z*w(22)
      w(22)=1163.D0/9.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 827.D0/8.D0 + w(22)
      w(23)= - 1001.D0/3.D0 + 163.D0/2.D0*z
      w(23)=z*w(23)
      w(23)=1295.D0/6.D0 + w(23)
      w(23)= - 161.D0/20.D0*w(2) + 1.D0/3.D0*w(23) + w(1)
      w(23)=w(2)*w(23)
      w(22)=1.D0/12.D0*w(23) + 1.D0/3.D0*w(22) + 1.D0/4.D0*w(1)
      w(22)=w(2)*w(22)
      w(8)=131.D0/4.D0 - w(8)
      w(8)=z3*w(8)
      w(7)=w(22) + 61.D0/8.D0*w(1) + 1.D0/18.D0*w(7) + w(8)
      w(7)=w(2)*w(7)
      w(8)=7373063994303022886557371476455187471130619052899966750211945
     & 64913.D0/6768.D0 + 9077128447003180654551044701099893838037562238
     & 1452878730411641883.D0/425.D0*z
      w(8)=z*w(8)
      w(8)=1495505922749518957747124586998049787146775559811038229906356
     & 119.D0/1219368.D0 + 1.D0/45325.D0*w(8)
      w(8)=z*w(8)
      w(8)=1170821911564083724648453467519104943081584649162971783672915
     & 3.D0/1167392439.D0 + 1.D0/62426.D0*w(8)
      w(8)=z*w(8)
      w(8)=4060948137695477682626293543635038512725757517507992271640821
     & .D0/88796339775.D0 + 1.D0/112.D0*w(8)
      w(8)=z*w(8)
      w(8)=2907286001932394948646859542831793497844283046864200013079.D0
     & /5221914778080.D0 + 1.D0/41971.D0*w(8)
      w(8)=z*w(8)
      w(8)=123180189112206299147751356135507875963981661105533289948579.
     & D0/432813487186080.D0 + w(8)
      w(8)=w(8)*w(13)
      w(8)=1191103705236810094692665250550202964376142172417784626677.D0
     & /188208289124856.D0 + w(8)
      w(8)=w(8)*w(10)
      w(8)=58140509727377222075363183324272043629440064408606570873.D0/
     & 197485720089840.D0 + w(8)
      w(8)=z*w(8)
      w(8)=387521522135070807744373640903517244254755917488549571.D0/
     & 4753654055230080.D0 + 1.D0/1849.D0*w(8)
      w(8)=z*w(8)
      w(8)=197137498276922832040508545210475082649181464721601697.D0/
     & 4720527546134400.D0 + w(8)
      w(8)=z*w(8)
      w(8)=75857476223404788882210515952962018770216900889.D0/
     & 23827380061485.D0 + 1.D0/6724.D0*w(8)
      w(8)=z*w(8)
      w(8)=9169211398077031752749544101024596290440842889651.D0/
     & 5615197042797342.D0 + w(8)
      w(8)=w(8)*w(12)
      w(8)=1707166509696792189220865726485604212607490199669.D0/
     & 10184648887340928.D0 + w(8)
      w(8)=z*w(8)
      w(8)=4281456149370951735889563189452525361496639105637.D0/
     & 49730819192208000.D0 + w(8)
      w(8)=w(8)*w(11)
      w(8)=4876055771181730609692768671818242179432243.D0/
     & 77464419832800.D0 + w(8)
      w(8)=z*w(8)
      w(8)=399400959244877352619145587323930933981529691.D0/
     & 12335286662899200.D0 + w(8)
      w(8)=z*w(8)
      w(8)=6513403887824109525272984976535421493795075659.D0/
     & 390748710961889280.D0 + w(8)
      w(8)=w(8)*w(14)
      w(8)=54545219435655765549464940551993408671081283.D0/
     & 107960873795174400.D0 + w(8)
      w(8)=z*w(8)
      w(8)=59099315456286294879305746147142263852659529.D0/
     & 226809847078214400.D0 + w(8)
      w(8)=z*w(8)
      w(8)=3468572773047951780695924293342755765788717.D0/
     & 25785393681047000.D0 + w(8)
      w(11)=1.D0/961.D0*z
      w(8)=w(8)*w(11)
      w(8)=6823231117040921380185076054801678718807.D0/
     & 94324633659055800.D0 + w(8)
      w(8)=z*w(8)
      w(8)=316490153070636033496319630614144610317.D0/8456691293570520.D
     & 0 + w(8)
      w(8)=z*w(8)
      w(8)=50155797340865397804643340438090563.D0/2175916786918875.D0
     &  + 1.D0/841.D0*w(8)
      w(8)=z*w(8)
      w(8)=1156124148342833749991143214056689787.D0/96707412751950000.D0
     &  + w(8)
      w(8)=z*w(8)
      w(8)=104921480124916946624785463381472733.D0/16898541259600000.D0
     &  + w(8)
      w(8)=z*w(8)
      w(8)=21240597230076172758754434776623433.D0/6577432213352000.D0
     &  + w(8)
      w(8)=w(8)*w(12)
      w(8)=4484595392027554276468823447773.D0/13327954748108.D0 + w(8)
      w(8)=z*w(8)
      w(8)=67658126361451956670396605594209.D0/385351735108340.D0 + 
     & w(8)
      w(8)=w(8)*w(13)
      w(8)=2114791572981325949668330837.D0/529784134880.D0 + w(8)
      w(8)=z*w(8)
      w(8)=2225225340233297292993873997.D0/1064482779360.D0 + w(8)
      w(8)=w(8)*w(6)
      w(8)=42466357901371580110508081.D0/116001328520.D0 + w(8)
      w(8)=z*w(8)
      w(8)=219873211981732581415024643.D0/1142918352576.D0 + w(8)
      w(8)=z*w(8)
      w(8)=2561781321020125831990763.D0/25179714560.D0 + w(8)
      w(8)=z*w(8)
      w(8)=1189231965390480743799751.D0/22217395200.D0 + w(8)
      w(8)=w(8)*w(14)
      w(8)=33473951609616979897.D0/19874400.D0 + w(8)
      w(8)=z*w(8)
      w(8)=72756357287243256371.D0/82892810.D0 + w(8)
      w(8)=z*w(8)
      w(8)=1955254059133301231.D0/4060056.D0 + w(8)
      w(8)=z*w(8)
      w(8)=7507240007803649843.D0/31491460.D0 + w(8)
      w(8)=w(8)*w(16)
      w(8)=45160804713054451.D0/4192650.D0 + w(8)
      w(8)=z*w(8)
      w(8)=932773102647509.D0/217800.D0 + w(8)
      w(8)=w(8)*w(10)
      w(8)=118766476660751.D0/423360.D0 + w(8)
      w(8)=z*w(8)
      w(8)=4577524649581.D0/197568.D0 + w(8)
      w(8)=w(8)*w(19)
      w(8)=40958403011.D0/6860.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 7167421903.D0/1225.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 5765051.D0/100.D0 + 1.D0/189.D0*w(8)
      w(8)=z*w(8)
      w(8)=51564763.D0/1620.D0 + w(8)
      w(8)=z*w(8)
      w(8)=128501.D0/27.D0 + 1.D0/50.D0*w(8)
      w(8)=w(8)*w(5)
      w(8)=518273.D0/81.D0 + w(8)
      w(8)=w(8)*w(15)
      w(8)=13630.D0/9.D0 + w(8)
      w(8)=z*w(8)
      w(10)=47 - 17.D0/3.D0*z
      w(10)=z*w(10)
      w(10)= - 1303 + 67.D0/2.D0*w(10)
      w(10)=z3*w(10)
      w(8)=w(10) - 598391.D0/432.D0 + w(8)
      w(10)= - 31817.D0/216.D0 + 18*z3
      w(10)=w(1)*w(10)
      w(7)=w(7) + 1.D0/3.D0*w(8) + w(10)
      w(7)=w(2)*w(7)
      w(8)=ln2**2
      w(10)=2*z2 - 1.D0/3.D0*w(8)
      w(8)=w(8)*w(10)
      w(8)=8*w(8) - 64*li4half
      w(10)= - 2*z - 1 + 2.D0/3.D0*w(1)
      w(8)=w(10)*w(8)
      w(10)= - 370857314159219731307454770205761641153623392045793419931
     & 387715038544131555516989102731.D0/2256.D0 - 
     & 17791607980075384601932264040579631453759475936089970289050395914
     & 79184632909137961107521.D0/5525.D0*z
      w(10)=w(10)*w(9)
      w(10)= - 144851805882074389081635630846665104546456362604941309380
     & 2252943361906943652213426885987.D0/21135712.D0 + w(10)
      w(10)=z*w(10)
      w(10)= - 712144618610261320587764527284065146580192454991095151118
     & 865110099349017933701935481.D0/683608185.D0 + 1.D0/33614.D0*
     & w(10)
      w(10)=w(10)*w(14)
      w(10)= - 227810123397456124310667258200858020604753744207499791973
     & 7471992120997873421407151.D0/72724275.D0 + w(10)
      w(10)=z*w(10)
      w(10)= - 550942374727067475168721828326010837642026064724357840323
     & 557943333731875415451.D0/3570324615.D0 + 1.D0/103823.D0*w(10)
      w(10)=w(10)*w(17)
      w(10)= - 100395150235041973842290836871753547057959662685919417104
     & 748091589401413006737.D0/8899939620.D0 + w(10)
      w(5)=w(10)*w(5)
      w(5)= - 5476202797295289828877925619324335976584951752946728509821
     & 71574551215557516267.D0/189636137691.D0 + w(5)
      w(5)=w(5)*w(18)
      w(5)= - 2733193342239836275820402328047600366995468959745671289817
     & 239375265554889.D0/57274581.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 5038011379469468741605250177629199465878088543999078520214
     & 41766113023.D0/26200108025.D0 + 1.D0/1272112.D0*w(5)
      w(5)=w(5)*w(12)
      w(5)= - 1912938188591979212637310853429413399176210841094143872174
     & 9136845161321.D0/9697669642752.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 4543576664696423491811777157831373289367149753831633224991
     & 4274961.D0/3092997676600.D0 + 1.D0/68921.D0*w(5)
      w(5)=z*w(5)
      w(5)= - 3347551451082308902673320647435611946869275448405862699469
     & 418146681.D0/443678620563360.D0 + w(5)
      w(5)=w(5)*w(20)
      w(5)= - 1904801679750545742143591818260825229726807710351661091644
     & 746103.D0/1964670005025.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 1819721578343005205269206133703624336939246709684068529553
     & 34096101.D0/364875271488000.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 5452748424894402799654322425808052760749486702412450899195
     & 17.D0/107572461796800.D0 + 1.D0/50653.D0*w(5)
      w(5)=z*w(5)
      w(5)= - 2354766817984087445596679616241057540364399753695971524969
     & 759.D0/901559679820800.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 9688105982641693420770540991659604444626852998778429240157
     & .D0/136645867438080.D0 + 1.D0/19.D0*w(5)
      w(5)=w(5)*w(12)
      w(5)= - 2310551277810740903280796997843016512367658763840203388269
     & .D0/315625160474624.D0 + w(5)
      w(5)=w(5)*w(13)
      w(5)= - 211052107355659824196010106961411981816282758714926902159.
     & D0/1283053087680000.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 130200088243077392765433936090099374533182475091753062123.
     & D0/1529851108349200.D0 + w(5)
      w(5)=w(5)*w(11)
      w(5)= - 138027638577736154983952233080788061017318687267075777.D0/
     & 3008760244308000.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 35958622059166496656874224841522922491894038882331.D0/
     & 1511209627200.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 27355712491615931315823857713800881814626396833.D0/
     & 53983404247500.D0 + 1.D0/24389.D0*w(5)
      w(5)=z*w(5)
      w(5)= - 101641697856363529958622455973316644022659974389.D0/
     & 385595744625000.D0 + w(5)
      w(5)=w(5)*w(19)
      w(5)= - 6169189534226915530429399487994430634588756857.D0/
     & 404271322000000.D0 + w(5)
      w(5)=w(5)*w(19)
      w(5)= - 557015771553833990254371366365708924779428307.D0/
     & 629419350560000.D0 + w(5)
      w(5)=w(5)*w(12)
      w(5)= - 1009317016508155789307682741538625624997409.D0/
     & 10904690248452.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 292668253610651746051415878345381574329043.D0/
     & 6034215947760.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 52187749111970057069829723964647361637.D0/1083649366800.D0
     &  + 1.D0/529.D0*w(5)
      w(5)=z*w(5)
      w(5)= - 1660356874987276237619136532045559689.D0/65484244800.D0
     &  + w(5)
      w(5)=w(5)*w(12)
      w(5)= - 14908001010863053991862145479994667.D0/5566160808.D0 + 
     & w(5)
      w(5)=z*w(5)
      w(5)= - 847339260397293579365798682850128911.D0/597629897280.D0
     &  + w(5)
      w(5)=z*w(5)
      w(5)= - 645692926726525441889261399541877.D0/309023769600.D0 + 1.D
     & 0/361.D0*w(5)
      w(5)=z*w(5)
      w(5)= - 235921621010622499014987289511779.D0/212075136000.D0 + 
     & w(5)
      w(5)=z*w(5)
      w(5)= - 4319684780438365743832852507.D0/2086812000.D0 + 1.D0/289.D
     & 0*w(5)
      w(5)=z*w(5)
      w(5)= - 41653774555827480955067329.D0/37678550.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 11747664276958461433785107.D0/19377540.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 27493238511238016419032719.D0/85885800.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 3107171950556824096847.D0/2858625.D0 + 1.D0/169.D0*w(5)
      w(5)=z*w(5)
      w(5)= - 555229350617841075383.D0/1039500.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 276045521555764133.D0/1058400.D0 + 1.D0/1331.D0*w(5)
      w(5)=z*w(5)
      w(5)= - 136614789497097049.D0/1481760.D0 + w(5)
      w(5)=w(5)*w(6)
      w(5)= - 496230248490439.D0/17150.D0 + w(5)
      w(5)=z*w(5)
      w(5)=26561002271.D0/6125.D0 + 1.D0/72.D0*w(5)
      w(3)=w(5)*w(3)
      w(3)=5317780441.D0/3000.D0 + w(3)
      w(3)=z*w(3)
      w(3)= - 2874732259.D0/1800.D0 + w(3)
      w(3)=w(3)*w(21)
      w(3)=9579313.D0/9.D0 + w(3)
      w(3)=w(3)*w(15)
      w(3)= - 2461039.D0/9.D0 + w(3)
      w(3)=z*w(3)
      w(3)=1634009.D0/2.D0 + w(3)
      w(3)=w(3)*w(19)
      w(3)=130505.D0/4.D0 + w(3)
      w(5)=607 + 261*z
      w(5)=z5*w(5)
      w(3)=1.D0/36.D0*w(3) + w(5)
      w(5)= - 2960722279.D0/376.D0 + 12206090957.D0/1275.D0*z
      w(5)=w(5)*w(9)
      w(5)=19882099099.D0/2438736.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 3603529601.D0/52585245.D0 + 1.D0/98.D0*w(5)
      w(5)=w(5)*w(20)
      w(5)=89170751.D0/4102395.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 5794854721.D0/323106300.D0 + w(5)
      w(5)=z*w(5)
      w(5)=12833983429.D0/563797080.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 69898595.D0/3714641.D0 + w(5)
      w(5)=z*w(5)
      w(5)=725340857.D0/30358860.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 3632923237.D0/183565200.D0 + w(5)
      w(5)=z*w(5)
      w(5)=2647585573.D0/105320800.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 1413419321.D0/67719990.D0 + w(5)
      w(5)=w(5)*w(20)
      w(5)=383375869.D0/57822453.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 241188473.D0/43698480.D0 + w(5)
      w(5)=z*w(5)
      w(5)=934146629.D0/133146720.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 821261363.D0/140193900.D0 + w(5)
      w(5)=z*w(5)
      w(5)=583485019.D0/78330560.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 1222605685.D0/195892224.D0 + w(5)
      w(5)=z*w(5)
      w(5)=2578039363.D0/324741120.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 148895147.D0/22295200.D0 + w(5)
      w(5)=z*w(5)
      w(5)=465201781.D0/54749100.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 639259297.D0/89011440.D0 + w(5)
      w(5)=z*w(5)
      w(5)=1602911.D0/175392.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 17177087.D0/2211300.D0 + w(5)
      w(5)=z*w(5)
      w(5)=90211831.D0/9126000.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 101169527.D0/11960000.D0 + w(5)
      w(5)=z*w(5)
      w(5)=601227799.D0/55862400.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 4341775.D0/467544.D0 + w(5)
      w(5)=z*w(5)
      w(5)=16103897.D0/1363670.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 11501507.D0/1117200.D0 + w(5)
      w(5)=z*w(5)
      w(5)=237952129.D0/18194400.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 74627.D0/6460.D0 + w(5)
      w(5)=z*w(5)
      w(5)=69494357.D0/4744224.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 43825717.D0/3329280.D0 + w(5)
      w(5)=z*w(5)
      w(5)=1013729.D0/60928.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 11699393.D0/764400.D0 + w(5)
      w(5)=z*w(5)
      w(5)=195079.D0/10140.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 3808555.D0/208208.D0 + w(5)
      w(5)=z*w(5)
      w(5)=17193493.D0/755040.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 2477411.D0/108900.D0 + w(5)
      w(5)=w(5)*w(6)
      w(5)=40807.D0/4400.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 1821157.D0/181440.D0 + w(5)
      w(5)=z*w(5)
      w(5)=2006719.D0/169344.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 87187.D0/5880.D0 + w(5)
      w(5)=z*w(5)
      w(5)=5078.D0/315.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 99841.D0/3600.D0 + w(5)
      w(5)=z*w(5)
      w(5)=10723.D0/480.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 4115.D0/36.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 40405.D0/72.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 10841.D0/6.D0 + w(5)
      w(5)=z*w(5)
      w(5)= - 11207.D0/12.D0 + w(5)
      w(5)=z3*w(5)
      w(6)= - 8409043.D0/648.D0 + 2077*z3
      w(6)=w(1)*w(6)
*
      RED0 = 1.D0/2.D0*w(3) + w(4) + 1.D0/3.D0*w(5) 
     &+ 1.D0/6.D0*w(6) + w(7) + w(8)
*
      RETURN
      END      
      REAL*8 FUNCTION RED1(z)
*     -----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 15, 2024
*     NF = 0 one particle reducible contribution to a_Qg^(3) for z \in [1/2,1]
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(29)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4,y
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      y = 1.0D0-z
*
      w(1)=Log(y)
      w(2)=1.D0/23.D0*y
      w(3)=297050441029112222257999.D0/12691.D0 + 470001863017695052897.
     & D0/20.D0*y
      w(3)=y*w(3)
      w(3)=6899861243595132606899.D0/296.D0 + w(3)
      w(3)=y*w(3)
      w(3)=80730810684744225784057.D0/3478.D0 + w(3)
      w(3)=y*w(3)
      w(3)=386264055228199206013.D0/10212.D0 + 1.D0/611.D0*w(3)
      w(3)=y*w(3)
      w(3)=326014372134704495951.D0/8658.D0 + w(3)
      w(3)=w(3)*w(2)
      w(3)=482836485529899468263.D0/296296.D0 + w(3)
      w(3)=y*w(3)
      w(3)=36366690826388890529.D0/6659926.D0 + 1.D0/297.D0*w(3)
      w(3)=y*w(3)
      w(3)=435166877393164019.D0/58545396.D0 + 1.D0/731.D0*w(3)
      w(3)=y*w(3)
      w(3)=21554194798479898787.D0/2914724358.D0 + w(3)
      w(3)=y*w(3)
      w(3)=185179993271737.D0/568726704.D0 + 1.D0/22591.D0*w(3)
      w(3)=y*w(3)
      w(3)=164894332224487951.D0/509223672594.D0 + w(3)
      w(3)=y*w(3)
      w(3)=159749969991537307.D0/496166655348.D0 + w(3)
      w(3)=y*w(3)
      w(3)=73243137390806753.D0/228841407522.D0 + w(3)
      w(3)=y*w(3)
      w(3)=12123096024532001.D0/38112374664.D0 + w(3)
      w(4)=1.D0/25.D0*y
      w(3)=w(3)*w(4)
      w(3)=1030541210317831.D0/81518134698.D0 + w(3)
      w(3)=y*w(3)
      w(3)=10358299710161387.D0/824885886825.D0 + w(3)
      w(3)=y*w(3)
      w(3)=558510443754469.D0/44790183900.D0 + w(3)
      w(3)=y*w(3)
      w(3)=230641941165243199.D0/18632716502400.D0 + w(3)
      w(3)=y*w(3)
      w(3)=55430860123276453.D0/4512611027925.D0 + w(3)
      w(5)=1.D0/31.D0*y
      w(3)=w(3)*w(5)
      w(3)=114433839435989.D0/291136195350.D0 + w(3)
      w(3)=y*w(3)
      w(3)=4935822956430911.D0/12664424497725.D0 + w(3)
      w(3)=y*w(3)
      w(3)=488608789040569.D0/1264936572900.D0 + w(3)
      w(3)=y*w(3)
      w(3)=155577690909673.D0/406586755575.D0 + w(3)
      w(3)=y*w(3)
      w(3)=49442242847929.D0/130509328950.D0 + w(3)
      w(6)=1.D0/2.D0*y
      w(3)=w(3)*w(6)
      w(3)=276702953861.D0/1476349875.D0 + w(3)
      w(3)=y*w(3)
      w(3)=148838872639.D0/803134332.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1691177501161.D0/9236044818.D0 + w(3)
      w(3)=y*w(3)
      w(3)=115720627817.D0/640179540.D0 + w(3)
      w(3)=y*w(3)
      w(3)=108942228487.D0/611080470.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1572564559.D0/8953560.D0 + w(3)
      w(3)=y*w(3)
      w(3)=50583789829.D0/292702410.D0 + w(3)
      w(7)=1.D0/4.D0*y
      w(3)=w(3)*w(7)
      w(3)=877669043.D0/20675655.D0 + w(3)
      w(3)=y*w(3)
      w(3)=26018080163.D0/624864240.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1410131293.D0/34594560.D0 + w(3)
      w(3)=y*w(3)
      w(3)=21519131.D0/540540.D0 + w(3)
      w(3)=y*w(3)
      w(3)=293409439.D0/7567560.D0 + w(3)
      w(3)=y*w(3)
      w(3)=44072233.D0/1171170.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1512169.D0/41580.D0 + w(3)
      w(3)=y*w(3)
      w(3)=15988631.D0/457380.D0 + w(3)
      w(3)=y*w(3)
      w(3)=75673.D0/2268.D0 + w(3)
      w(3)=y*w(3)
      w(3)=2146549.D0/68040.D0 + w(3)
      w(3)=y*w(3)
      w(3)=254389.D0/8640.D0 + w(3)
      w(3)=w(3)*w(6)
      w(3)=12739.D0/945.D0 + w(3)
      w(3)=y*w(3)
      w(3)=12943.D0/1080.D0 + w(3)
      w(3)=y*w(3)
      w(3)=10939.D0/1080.D0 + w(3)
      w(3)=y*w(3)
      w(3)=2227.D0/288.D0 + w(3)
      w(3)=y*w(3)
      w(3)=467.D0/108.D0 + w(3)
      w(3)=y*w(3)
      w(3)=161.D0/24.D0 + w(3)
      w(3)=w(3)*w(6)
      w(3)=28.D0/9.D0 + w(3)
      w(8)=1.D0/3.D0*y
      w(3)=w(3)*w(8)
      w(9)=1.D0/27.D0*y
      w(10)= - 475.D0/7.D0 - 679.D0/10.D0*y
      w(10)=y*w(10)
      w(10)= - 1085.D0/16.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 3185.D0/47.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 3115.D0/46.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 203.D0/3.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2975.D0/44.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2905.D0/43.D0 + w(10)
      w(10)=w(10)*w(9)
      w(10)= - 5.D0/2.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2765.D0/1107.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 539.D0/216.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 875.D0/351.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2555.D0/1026.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2485.D0/999.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 805.D0/324.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 67.D0/27.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2275.D0/918.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 245.D0/99.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2135.D0/864.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2065.D0/837.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 133.D0/54.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1925.D0/783.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 265.D0/108.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 595.D0/243.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1715.D0/702.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 329.D0/135.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 175.D0/72.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1505.D0/621.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1435.D0/594.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 65.D0/27.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 259.D0/108.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1225.D0/513.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 385.D0/162.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1085.D0/459.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1015.D0/432.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 7.D0/3.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 125.D0/54.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 805.D0/351.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 245.D0/108.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 665.D0/297.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 119.D0/54.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 175.D0/81.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 455.D0/216.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 55.D0/27.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 35.D0/18.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 49.D0/27.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 175.D0/108.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 35.D0/27.D0 + w(10)
      w(10)=y*w(10)
      w(10)=215.D0/54.D0 + w(10)
      w(11)=1.D0/16.D0*y
      w(10)=w(10)*w(11)
      w(10)= - 10.D0/27.D0 + w(10)
      w(10)=y*w(10)
      w(10)=5.D0/27.D0 + w(10)
      w(10)=w(1)*w(10)
      w(3)=w(3) + w(10)
      w(3)=w(1)*w(3)
      w(10)=1.D0/1849.D0*y
      w(12)= - 40384392388890914685437850260712760575601430532410885783.
     & D0/504.D0 - 11745722722451933296339012247806596511004302358838072
     & 61319.D0/7475.D0*y
      w(12)=y*w(12)
      w(12)= - 438017144907118607619308597384129742991913380255163063.D0
     & /21528.D0 + 1.D0/2009.D0*w(12)
      w(12)=y*w(12)
      w(12)= - 371529194803929597365624156302375855211023657484083136966
     & 3.D0/357803433.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 1127792629036274536426777518616603899247306499051221.D0/
     & 1409785.D0 + 1.D0/6627.D0*w(12)
      w(12)=y*w(12)
      w(12)= - 133895747092043301904818538314092664176676376745499089.D0
     & /327683070.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 1520513630400239397897225796617309479012326321918289.D0/
     & 7281846.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 1519428986413594715073307378128016885513663866354577.D0/
     & 14232699.D0 + w(12)
      w(12)=w(12)*w(10)
      w(12)= - 205284472686101272764427086109687281906455349511.D0/
     & 6950853.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 4102415673558519644420634415646255547679875742477.D0/
     & 271414260.D0 + w(12)
      w(13)=1.D0/41.D0*y
      w(12)=w(12)*w(13)
      w(12)= - 17227818501023876883668831922209671982079019.D0/91260.D0
     &  + w(12)
      w(12)=y*w(12)
      w(12)= - 578600070349005500008503353013614191987584233.D0/5982093.
     & D0 + w(12)
      w(14)=1.D0/7.D0*y
      w(12)=w(12)*w(14)
      w(12)= - 3259667624834307797524086956825554076852311.D0/460161.D0
     &  + w(12)
      w(12)=y*w(12)
      w(12)= - 2170900343704216407331112214658297167462307.D0/597402.D0
     &  + w(12)
      w(15)=1.D0/1369.D0*y
      w(12)=w(12)*w(15)
      w(12)= - 5390109983896293138462408909924493605947.D0/3955770.D0
     &  + w(12)
      w(12)=y*w(12)
      w(12)= - 23535050374602251255161717937977981159969.D0/33624045.D0
     &  + w(12)
      w(12)=y*w(12)
      w(12)= - 1803786741327213534448013510884852783.D0/621621.D0 + 1.D0
     & /124.D0*w(12)
      w(12)=y*w(12)
      w(12)= - 460274777531595134381495711612865826241.D0/308324016.D0
     &  + w(12)
      w(12)=y*w(12)
      w(12)= - 19150258711484439007930394373683409787.D0/24915072.D0 + 
     & w(12)
      w(12)=y*w(12)
      w(12)= - 4676118373095091251557522333457885683.D0/11805885.D0 + 
     & w(12)
      w(12)=w(12)*w(5)
      w(12)= - 1862066922071465676452443599361253.D0/282555.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 234182515843729512719658241271077289.D0/68830398.D0 + 
     & w(12)
      w(12)=y*w(12)
      w(12)= - 875695936247900979003782950693.D0/7120386.D0 + 1.D0/
     & 14297.D0*w(12)
      w(12)=y*w(12)
      w(12)= - 9194886993212254289538482107567.D0/144514773.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 137186943912650076780595956781.D0/4162977.D0 + w(12)
      w(12)=w(12)*w(8)
      w(12)= - 238439861732190965384631703.D0/41860.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 27831940620336746074397159.D0/1646892.D0 + 1.D0/175.D0*
     & w(12)
      w(12)=y*w(12)
      w(12)= - 14948089194255898887977159119.D0/1701376677.D0 + w(12)
      w(12)=w(12)*w(2)
      w(12)= - 213195654010044925636163.D0/1072071.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 1273075732951659779995283.D0/12280086.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 474694280771085070569437.D0/8771490.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 157135173438020021314723.D0/5555277.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 26967266705634112961.D0/2631447.D0 + 1.D0/1444.D0*w(12)
      w(12)=y*w(12)
      w(12)= - 213269146571280263147.D0/39764088.D0 + w(12)
      w(16)=1.D0/17.D0*y
      w(12)=w(12)*w(16)
      w(12)= - 907970638955160733.D0/5503680.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 5751980161526563.D0/66885.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 5356510205009461.D0/120393.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 10127558749265033.D0/447174.D0 + w(12)
      w(17)=1.D0/13.D0*y
      w(12)=w(12)*w(17)
      w(12)= - 2267681851481.D0/2646.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 9480198037193.D0/24255.D0 + w(12)
      w(12)=y*w(12)
      w(12)= - 5986973621.D0/19845.D0 + 1.D0/484.D0*w(12)
      w(12)=y*w(12)
      w(12)= - 457671871.D0/11907.D0 + w(12)
      w(12)=y*w(12)
      w(12)=1338986573.D0/14112.D0 + w(12)
      w(12)=y*w(12)
      w(12)=7684423.D0/49.D0 + w(12)
      w(12)=w(12)*w(6)
      w(12)=794123.D0/9.D0 + w(12)
      w(12)=y*w(12)
      w(12)=1636681.D0/18.D0 + w(12)
      w(12)=w(12)*w(4)
      w(12)=35719.D0/24.D0 + w(12)
      w(12)=w(12)*w(6)
      w(12)=2207.D0/9.D0 + w(12)
      w(12)=y*w(12)
      w(12)=17.D0/2.D0 + w(12)
      w(12)=w(12)*w(11)
      w(12)= - 82 + w(12)
      w(12)=y*w(12)
      w(12)= - 10 + w(12)
      w(18)=1.D0/49.D0*y
      w(19)= - 92003.D0/2.D0 - 1150657.D0/25.D0*y
      w(19)=w(19)*w(18)
      w(19)= - 352793.D0/376.D0 + w(19)
      w(19)=w(19)*w(6)
      w(19)= - 506843.D0/1081.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 107773.D0/230.D0 + w(19)
      w(19)=w(19)*w(6)
      w(19)= - 38633.D0/165.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 885391.D0/3784.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 140759.D0/602.D0 + w(19)
      w(20)=1.D0/9.D0*y
      w(19)=w(19)*w(20)
      w(19)= - 29803.D0/1148.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 191443.D0/7380.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 242609.D0/9360.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 8857.D0/342.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 654829.D0/25308.D0 + w(19)
      w(19)=w(19)*w(7)
      w(19)= - 2152.D0/333.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 195229.D0/30240.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 276281.D0/42840.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 173467.D0/26928.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 13589.D0/2112.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 458971.D0/71424.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 71617.D0/11160.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 133799.D0/20880.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 46757.D0/7308.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 12877.D0/2016.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 53711.D0/8424.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 22909.D0/3600.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 22861.D0/3600.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 27979.D0/4416.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 115127.D0/18216.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 69887.D0/11088.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 3959.D0/630.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 171367.D0/27360.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 8537.D0/1368.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 45643.D0/7344.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 30289.D0/4896.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 35449.D0/5760.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 1713.D0/280.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 79621.D0/13104.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 217.D0/36.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 18917.D0/3168.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 23381.D0/3960.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 4193.D0/720.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 2473.D0/432.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 22579.D0/4032.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 2741.D0/504.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 1253.D0/240.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 221.D0/45.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 1277.D0/288.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 259.D0/72.D0 + w(19)
      w(19)=y*w(19)
      w(19)=79.D0/48.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 20.D0/3.D0 + w(19)
      w(19)=y*w(19)
      w(19)=10.D0/3.D0 + w(19)
      w(19)=z2*w(19)
      w(3)=w(3) + 1.D0/9.D0*w(12) + w(19)
      w(3)=w(1)*w(3)
      w(12)=1.D0/343.D0*y
      w(19)= - 166729223244115286170556405759524919769237414773521404232
     & 2411725106461495020469.D0/182.D0 - 763945161121583226509531671284
     & 4449107962308407574947701528054340703748234227647.D0/425.D0*y
      w(19)=w(19)*w(12)
      w(19)= - 240744008363184549190874799130415025761325281436765761978
     & 43312731194341866231.D0/1768.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 165914491468954581343602022183317243772354232895840273483
     & 5279462971945437419867.D0/238901.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 166811254525766351262651076290716473754530105619053118917
     & 17797298813643.D0/1955.D0 + 1.D0/415292.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 127054545222597685374872684049858475170352990019368798032
     & 5173070636017459.D0/291720.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 218016652588550851804681675348487301374968773875147652441
     & 34916810329.D0/2431.D0 + 1.D0/248.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 302344672327358033691031019013926779385845928318752746116
     & 682399599683399.D0/65988832.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 376684282154100712345642149406635743420016672323055605555
     & 5544797.D0/383656.D0 + 1.D0/238521.D0*w(19)
      w(21)=1.D0/19.D0*y
      w(19)=w(19)*w(21)
      w(19)= - 142515827612175778873885420263550943664970823524769734298
     & 918821183.D0/539310720.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 159341191503025941440188870553057627806179834859066800668
     & 4853.D0/812254560.D0 + 1.D0/68921.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 182329712718345321016066290130428846736266346822752459849
     & 689.D0/181562784.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 196973587600432094579550650929059281718700372600907953381
     & 9.D0/3829488.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 243763019208260746673125397649850837831984414566498616696
     & 603.D0/924720576.D0 + w(19)
      w(22)=1.D0/50653.D0*y
      w(19)=w(19)*w(22)
      w(19)= - 500239214022001089115888920910311580116779796814264297.D0
     & /187443360.D0 + w(19)
      w(23)=1.D0/11.D0*y
      w(19)=w(19)*w(23)
      w(19)= - 925645860921182595517820241601049671795426218162813757.D0
     & /7435253280.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 298762534601372112402835314877166402157253976887735217.D0
     & /4673587776.D0 + w(19)
      w(19)=w(19)*w(14)
      w(19)= - 18921028237137228407240593730779239409590075714955667.D0/
     & 4032114944.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 148600083603314307839377846618731775584776799148957721.D0
     & /61581391872.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 161898697356512348290477393804933246682598663553507.D0/
     & 130369260.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 993224372356389822336629807414378550424411009.D0/1492260.
     & D0 + 1.D0/961.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 31793295157828177282366968279449305425632194913.D0/
     & 92660568.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 43811595466818565126274627094062639591707.D0/18106088.D0
     &  + 1.D0/73167.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 1701976803915396724427336483163010988931851.D0/
     & 1361836476.D0 + w(19)
      w(19)=w(19)*w(8)
      w(19)= - 3621285981586280492347100718911618638393.D0/16812796.D0
     &  + w(19)
      w(19)=y*w(19)
      w(19)= - 720917101773587961745612760887691860583.D0/6466460.D0 + 
     & w(19)
      w(19)=y*w(19)
      w(19)= - 4950258149956577384362498075895729.D0/739024.D0 + 1.D0/
     & 8625.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 4615831497858431446700348310400918747.D0/1328058732.D0
     &  + w(19)
      w(19)=y*w(19)
      w(19)= - 25497221244723461737405310543231.D0/74687613.D0 + 1.D0/
     & 5290.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 482502112623758554513255804635497.D0/2715913200.D0 + 
     & w(19)
      w(19)=y*w(19)
      w(19)= - 691046762027548602871041741997.D0/7461300.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 7133921134962346962192436112669.D0/147435288.D0 + w(19)
      w(24)=1.D0/361.D0*y
      w(19)=w(19)*w(24)
      w(19)= - 22058747989974543867870509.D0/314160.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 466040119618531924120651757.D0/12623520.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 1297849841282916373759837.D0/19219200.D0 + 1.D0/289.D0*
     & w(19)
      w(19)=w(19)*w(8)
      w(19)= - 420321936752374084283.D0/35035.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 185673995979093808129.D0/28665.D0 + w(19)
      w(19)=w(19)*w(7)
      w(19)= - 522436341166192800251.D0/585585.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 83122468529435399.D0/27720.D0 + 1.D0/169.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 1355053810021747037.D0/762300.D0 + w(19)
      w(25)=1.D0/121.D0*y
      w(19)=w(19)*w(25)
      w(19)= - 38546345660941.D0/4200.D0 + w(19)
      w(19)=y*w(19)
      w(19)= - 5709957317011.D0/945.D0 + w(19)
      w(19)=w(19)*w(8)
      w(19)= - 3073236534197.D0/2240.D0 + w(19)
      w(19)=w(19)*w(7)
      w(19)= - 1152421439.D0/5.D0 + w(19)
      w(19)=w(19)*w(18)
      w(19)= - 54757001.D0/20.D0 + w(19)
      w(19)=w(19)*w(8)
      w(19)= - 2530193.D0/4.D0 + w(19)
      w(19)=w(19)*w(4)
      w(19)=489233.D0/48.D0 + w(19)
      w(19)=y*w(19)
      w(19)=1153 + 1.D0/18.D0*w(19)
      w(19)=y*w(19)
      w(19)= - 4507.D0/6.D0 + w(19)
      w(19)=y*w(19)
      w(26)=1.D0/47.D0*y
      w(27)=10308579262870816745413723.D0 + 
     & 8280439412289559679263263019.D0/800.D0*y
      w(27)=y*w(27)
      w(27)=15433063687282555502901599.D0/18048.D0 + 1.D0/12005.D0*
     & w(27)
      w(27)=y*w(27)
      w(27)=73633189764306744541033487.D0/86480.D0 + w(27)
      w(27)=w(27)*w(26)
      w(27)=331862725481933679279841.D0/18400.D0 + w(27)
      w(27)=y*w(27)
      w(27)=355491460796685293101187.D0/19800.D0 + w(27)
      w(27)=y*w(27)
      w(27)=18933723532831274049809141.D0/1059520.D0 + w(27)
      w(27)=w(27)*w(23)
      w(27)=54502659858020329606037.D0/33712.D0 + w(27)
      w(27)=y*w(27)
      w(27)=2775069676163967586211.D0/964320.D0 + 1.D0/559.D0*w(27)
      w(27)=y*w(27)
      w(27)=8545305275036290295929.D0/2984800.D0 + w(27)
      w(13)=w(27)*w(13)
      w(13)=262933963922930872973.D0/3785600.D0 + w(13)
      w(13)=w(13)*w(2)
      w(13)=16200790378264677817.D0/5394480.D0 + w(13)
      w(13)=y*w(13)
      w(13)=140597419401423968047.D0/47084128.D0 + w(13)
      w(13)=y*w(13)
      w(13)=4597565041154820887.D0/1548820.D0 + w(13)
      w(27)=1.D0/37.D0*y
      w(13)=w(13)*w(27)
      w(13)=16822115871801276013.D0/210974400.D0 + w(13)
      w(13)=y*w(13)
      w(13)=3774957277975563697.D0/47647600.D0 + w(13)
      w(13)=y*w(13)
      w(13)=27106372047468396721.D0/344424080.D0 + w(13)
      w(28)=1.D0/5.D0*y
      w(13)=w(13)*w(28)
      w(13)=3040245698347848967.D0/194498304.D0 + w(13)
      w(13)=y*w(13)
      w(13)=283537121799818249873.D0/18271052800.D0 + w(13)
      w(13)=y*w(13)
      w(13)=10992114634971017089.D0/713713000.D0 + w(13)
      w(13)=w(13)*w(5)
      w(13)=1974532066209618241.D0/4006002000.D0 + w(13)
      w(13)=y*w(13)
      w(13)=17573448686216867.D0/35951300.D0 + w(13)
      w(29)=1.D0/29.D0*y
      w(13)=w(13)*w(29)
      w(13)=1292563797965317.D0/77357280.D0 + w(13)
      w(13)=y*w(13)
      w(13)=2058318478580207.D0/124324200.D0 + w(13)
      w(13)=w(13)*w(8)
      w(13)=16357877879742457.D0/2992990000.D0 + w(13)
      w(13)=y*w(13)
      w(13)=77840742660391.D0/14389375.D0 + w(13)
      w(13)=y*w(13)
      w(13)=90673269966215.D0/16944928.D0 + w(13)
      w(13)=y*w(13)
      w(13)=246483203655713.D0/46598552.D0 + w(13)
      w(13)=y*w(13)
      w(13)=4602273712919.D0/880880.D0 + w(13)
      w(9)=w(13)*w(9)
      w(9)=1338100013.D0/7007.D0 + w(9)
      w(9)=y*w(9)
      w(9)=8593976959343.D0/45645600.D0 + w(9)
      w(9)=y*w(9)
      w(9)=1269519340307.D0/6846840.D0 + w(9)
      w(9)=w(9)*w(21)
      w(9)=317520787621.D0/33081048.D0 + w(9)
      w(9)=y*w(9)
      w(9)=1385935407403.D0/147026880.D0 + w(9)
      w(9)=w(9)*w(16)
      w(9)=2411079773.D0/4435200.D0 + w(9)
      w(9)=y*w(9)
      w(9)=30187585381.D0/56756700.D0 + w(9)
      w(9)=y*w(9)
      w(9)=1042235543.D0/2007720.D0 + w(9)
      w(9)=y*w(9)
      w(9)=4550977.D0/9009.D0 + w(9)
      w(9)=y*w(9)
      w(9)=1344118921.D0/2744280.D0 + w(9)
      w(9)=y*w(9)
      w(9)=360352667.D0/762300.D0 + w(9)
      w(9)=w(9)*w(6)
      w(9)=39687.D0/175.D0 + w(9)
      w(9)=y*w(9)
      w(9)=8394649.D0/38880.D0 + w(9)
      w(9)=y*w(9)
      w(9)=34443323.D0/169344.D0 + w(9)
      w(9)=y*w(9)
      w(9)=2496799.D0/13230.D0 + w(9)
      w(9)=y*w(9)
      w(9)=308093.D0/1800.D0 + w(9)
      w(9)=y*w(9)
      w(9)=100987.D0/675.D0 + w(9)
      w(9)=y*w(9)
      w(9)=105875.D0/864.D0 + w(9)
      w(9)=y*w(9)
      w(9)=9731.D0/108.D0 + w(9)
      w(9)=y*w(9)
      w(9)=1675.D0/24.D0 + w(9)
      w(9)=w(9)*w(6)
      w(9)=8 + w(9)
      w(9)=z2*w(9)
      w(9)=w(9) + 4 + 1.D0/144.D0*w(19)
      w(9)=y*w(9)
      w(13)=45431.D0/4.D0 + 397783.D0/35.D0*y
      w(13)=w(13)*w(14)
      w(13)=609655.D0/376.D0 + w(13)
      w(13)=w(13)*w(7)
      w(13)=437876.D0/1081.D0 + w(13)
      w(13)=y*w(13)
      w(13)=186191.D0/460.D0 + w(13)
      w(13)=y*w(13)
      w(13)=53387.D0/132.D0 + w(13)
      w(13)=y*w(13)
      w(13)=1529177.D0/3784.D0 + w(13)
      w(13)=y*w(13)
      w(13)=121535.D0/301.D0 + w(13)
      w(13)=w(13)*w(20)
      w(13)=7351.D0/164.D0 + w(13)
      w(13)=y*w(13)
      w(13)=660967.D0/14760.D0 + w(13)
      w(13)=y*w(13)
      w(13)=83747.D0/1872.D0 + w(13)
      w(13)=y*w(13)
      w(13)=7642.D0/171.D0 + w(13)
      w(13)=y*w(13)
      w(13)=1129775.D0/25308.D0 + w(13)
      w(13)=y*w(13)
      w(13)=59393.D0/1332.D0 + w(13)
      w(13)=y*w(13)
      w(13)=336683.D0/7560.D0 + w(13)
      w(13)=y*w(13)
      w(13)=6805.D0/153.D0 + w(13)
      w(13)=y*w(13)
      w(13)=299009.D0/6732.D0 + w(13)
      w(13)=y*w(13)
      w(13)=46835.D0/1056.D0 + w(13)
      w(13)=y*w(13)
      w(13)=790709.D0/17856.D0 + w(13)
      w(13)=y*w(13)
      w(13)=61672.D0/1395.D0 + w(13)
      w(13)=y*w(13)
      w(13)=46073.D0/1044.D0 + w(13)
      w(13)=y*w(13)
      w(13)=321901.D0/7308.D0 + w(13)
      w(13)=y*w(13)
      w(13)=1055.D0/24.D0 + w(13)
      w(13)=y*w(13)
      w(13)=46187.D0/1053.D0 + w(13)
      w(13)=y*w(13)
      w(13)=39383.D0/900.D0 + w(13)
      w(13)=y*w(13)
      w(13)=15713.D0/360.D0 + w(13)
      w(13)=y*w(13)
      w(13)=48053.D0/1104.D0 + w(13)
      w(13)=y*w(13)
      w(13)=98810.D0/2277.D0 + w(13)
      w(13)=y*w(13)
      w(13)=119893.D0/2772.D0 + w(13)
      w(13)=y*w(13)
      w(13)=7757.D0/180.D0 + w(13)
      w(13)=y*w(13)
      w(13)=58717.D0/1368.D0 + w(13)
      w(13)=y*w(13)
      w(13)=7307.D0/171.D0 + w(13)
      w(13)=y*w(13)
      w(13)=78065.D0/1836.D0 + w(13)
      w(13)=y*w(13)
      w(13)=103507.D0/2448.D0 + w(13)
      w(13)=y*w(13)
      w(13)=60503.D0/1440.D0 + w(13)
      w(13)=y*w(13)
      w(13)=292.D0/7.D0 + w(13)
      w(13)=y*w(13)
      w(13)=19361.D0/468.D0 + w(13)
      w(13)=y*w(13)
      w(13)=1475.D0/36.D0 + w(13)
      w(13)=y*w(13)
      w(13)=32083.D0/792.D0 + w(13)
      w(13)=y*w(13)
      w(13)=19781.D0/495.D0 + w(13)
      w(13)=y*w(13)
      w(13)=1415.D0/36.D0 + w(13)
      w(13)=y*w(13)
      w(13)=8317.D0/216.D0 + w(13)
      w(13)=y*w(13)
      w(13)=37805.D0/1008.D0 + w(13)
      w(13)=y*w(13)
      w(13)=326.D0/9.D0 + w(13)
      w(13)=y*w(13)
      w(13)=2071.D0/60.D0 + w(13)
      w(13)=y*w(13)
      w(13)=1157.D0/36.D0 + w(13)
      w(13)=y*w(13)
      w(13)=2059.D0/72.D0 + w(13)
      w(13)=y*w(13)
      w(13)=205.D0/9.D0 + w(13)
      w(13)=y*w(13)
      w(13)= - 343.D0/36.D0 + w(13)
      w(13)=y*w(13)
      w(13)=400.D0/9.D0 + w(13)
      w(13)=y*w(13)
      w(13)= - 200.D0/9.D0 + w(13)
      w(13)=z3*w(13)
      w(3)=w(3) + w(13) - 40.D0/9.D0 + w(9)
      w(3)=w(1)*w(3)
      w(9)= - 1533711289160135489600915096708201168534422500043517371061
     & .D0 - 93237546970432758667735291103941532841449629180102910305513
     & 79.D0/3100.D0*y
      w(9)=y*w(9)
      w(9)= - 690593758546119254025622935049119909195302674436404205369.
     & D0/23312.D0 + 1.D0/26411.D0*w(9)
      w(9)=y*w(9)
      w(9)= - 5572937921895021333806071816715874745325293783818846229103
     & .D0/368621.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 11915374297409594793759471836724084436721480765889863.D0/
     & 3410.D0 + 1.D0/2209.D0*w(9)
      w(9)=y*w(9)
      w(9)= - 33473936769425270684557522604569375930066776165628259.D0/
     & 18755.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 588438774844883472937284614201276946099507303844587373.D0/
     & 645172.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 6837430436210497900324689581823240836782212747075719.D0/
     & 14663.D0 + w(9)
      w(9)=w(9)*w(10)
      w(9)= - 25249990121816018753626574117781594296726156293877.D0/
     & 195734.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 9230435252646298569930525618181639732642307714131.D0/
     & 139810.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 1783079210160203854943502581825632541402349451.D0/88660.D0
     &  + 1.D0/1681.D0*w(9)
      w(9)=w(9)*w(6)
      w(9)= - 433950050606279979000387513749079580543940521.D0/84227.D0
     &  + w(9)
      w(9)=y*w(9)
      w(9)= - 24353478084838598291848168322801266445367803.D0/239723.D0
     &  + 1.D0/26.D0*w(9)
      w(9)=y*w(9)
      w(9)= - 34191679813488825585558894080049019175303201.D0/656084.D0
     &  + w(9)
      w(9)=w(9)*w(15)
      w(9)= - 24255494130981753506980938673129073733557.D0/1241240.D0
     &  + w(9)
      w(9)=w(9)*w(6)
      w(9)= - 1203496820392717694272156623858074541809.D0/239785.D0 + 
     & w(9)
      w(9)=y*w(9)
      w(9)= - 68442876309159132058333680652790274142637.D0/26527072.D0
     &  + w(9)
      w(9)=y*w(9)
      w(9)= - 16569888489016379949657667838626993047491.D0/12483328.D0
     &  + w(9)
      w(9)=w(9)*w(14)
      w(9)= - 13738937388241637118001554585874121901313.D0/140721152.D0
     &  + w(9)
      w(9)=w(9)*w(29)
      w(9)= - 33377786509609823434525612209115312987.D0/19239220.D0 + 
     & w(9)
      w(5)=w(9)*w(5)
      w(5)= - 971997657499696496214185643140273999.D0/33673640.D0 + 
     & w(5)
      w(5)=y*w(5)
      w(5)= - 175636467399299873709796743449640523.D0/11785774.D0 + 
     & w(5)
      w(5)=y*w(5)
      w(5)= - 902228279828064739311760025519.D0/10192.D0 + 1.D0/87.D0*
     & w(5)
      w(5)=y*w(5)
      w(5)= - 64363709720732918313763216660303.D0/1405404.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 170064109654389815949827660851.D0/236600.D0 + 1.D0/33.D0*
     & w(5)
      w(5)=w(5)*w(2)
      w(5)= - 17843678509646019588358437283.D0/1101100.D0 + w(5)
      w(9)=1.D0/75.D0*y
      w(5)=w(5)*w(9)
      w(5)= - 10455087522119911238056822481.D0/93197104.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 22420491275859921078792518611.D0/384438054.D0 + w(5)
      w(5)=w(5)*w(2)
      w(5)= - 1220878169225605681960667.D0/924924.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 795508598661430292576999.D0/1156155.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 1502699294788191960489199.D0/4184180.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 117790533880990237462907.D0/627627.D0 + w(5)
      w(5)=w(5)*w(24)
      w(5)= - 3664764025659674228293.D0/13477464.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 2556625213558399639477.D0/17969952.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 111054648563553421.D0/3020160.D0 + 1.D0/2023.D0*w(5)
      w(5)=y*w(5)
      w(5)= - 22172282938113617.D0/1156155.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 59785678020229097.D0/6012006.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 90661976571391.D0/17787.D0 + w(5)
      w(5)=w(5)*w(17)
      w(5)= - 76886401285663.D0/391314.D0 + w(5)
      w(5)=w(5)*w(6)
      w(5)= - 46038880277441.D0/978285.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 59965738691.D0/2940.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 23204232361.D0/3528.D0 + w(5)
      w(5)=w(5)*w(20)
      w(5)=729023779.D0/10976.D0 + w(5)
      w(5)=w(5)*w(6)
      w(5)=254823103.D0/1029.D0 + w(5)
      w(5)=y*w(5)
      w(5)=3775969.D0/10.D0 + w(5)
      w(5)=y*w(5)
      w(5)=7752929.D0/15.D0 + w(5)
      w(5)=w(5)*w(4)
      w(5)=103349.D0/8.D0 + w(5)
      w(5)=w(5)*w(6)
      w(5)= - 21593.D0/3.D0 + w(5)
      w(5)=w(5)*w(8)
      w(5)=6791.D0/2.D0 + w(5)
      w(5)=y*w(5)
      w(5)= - 275 + 1.D0/8.D0*w(5)
      w(5)=y*w(5)
      w(5)=211.D0/2.D0 + w(5)
      w(10)= - 89995.D0/4.D0 - 2814209.D0/125.D0*y
      w(10)=w(10)*w(18)
      w(10)= - 862621.D0/1880.D0 + w(10)
      w(10)=w(10)*w(6)
      w(10)= - 1239124.D0/5405.D0 + w(10)
      w(10)=w(10)*w(8)
      w(10)= - 17563.D0/230.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 377689.D0/4950.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 432727.D0/5676.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 343918.D0/4515.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 218417.D0/2870.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 187037.D0/2460.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 592453.D0/7800.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 56224.D0/741.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1598453.D0/21090.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 84031.D0/1110.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 95269.D0/1260.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 673942.D0/8925.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 84607.D0/1122.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 66261.D0/880.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1118663.D0/14880.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 6980.D0/93.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 325903.D0/4350.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 91079.D0/1218.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 94027.D0/1260.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 130678.D0/1755.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 144853.D0/1950.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 111139.D0/1500.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 13595.D0/184.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 279544.D0/3795.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 169591.D0/2310.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 15361.D0/210.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 415259.D0/5700.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1378.D0/19.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 110411.D0/1530.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 146389.D0/2040.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 17113.D0/240.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 12388.D0/175.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 38329.D0/546.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 27113.D0/390.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 45361.D0/660.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 11186.D0/165.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 10001.D0/150.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 2351.D0/36.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 53423.D0/840.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 6448.D0/105.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 117.D0/2.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 8167.D0/150.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 581.D0/12.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 578.D0/15.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1931.D0/30.D0 + w(10)
      w(10)=y*w(10)
      w(10)=76.D0/15.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 38.D0/15.D0 + w(10)
      w(10)=z2*w(10)
      w(5)=1.D0/18.D0*w(5) + w(10)
      w(5)=z2*w(5)
      w(10)=284630875939348464831641017177815015571333522236372200655927
     & 39068557064438932960468639036488953574011.D0/60088.D0 + 
     & 30184659660218766039037308533225962822807066087067468008744360257
     & 4115598783980753901134001572622347.D0/325.D0*y
      w(10)=y*w(10)
      w(10)=320945486039599487862813241968276968704675112877031121664609
     & 6582122851270615851684010287732724667.D0/223184.D0 + 1.D0/16807.D
     & 0*w(10)
      w(10)=y*w(10)
      w(10)=911507053191460825969051457304905298431443106375927875221256
     & 3480515648676496194900474174108002187.D0/105552083.D0 + 1.D0/85.D
     & 0*w(10)
      w(10)=y*w(10)
      w(10)=172512580500546037176329014788272915464127660052838190927667
     & 70569310913571686534351057575091.D0/1908920650.D0 + 1.D0/4879681.
     & D0*w(10)
      w(10)=y*w(10)
      w(10)=972723993993716623588081877051328351007871669855015705196042
     & 125210492237935141056229230787.D0/210683550.D0 + w(10)
      w(10)=w(10)*w(21)
      w(10)=453573374733266972285349452296979064757188765867820139491472
     & 20892877068784328045620274763.D0/365184820.D0 + w(10)
      w(10)=y*w(10)
      w(10)=119161531374785911883741785020635917348781302203419161514562
     & 59328548846707809060742958263.D0/47465726945.D0 + 1.D0/253.D0*
     & w(10)
      w(10)=y*w(10)
      w(10)=399712980602654508233023008973408211770306917424101679767653
     & 86202342908568078849.D0/15958576634.D0 + 1.D0/51282015.D0*w(10)
      w(10)=y*w(10)
      w(10)=314599685689733294463204082860932953466864244138663764912323
     & 228573281663176875629641.D0/245363115747750.D0 + w(10)
      w(10)=y*w(10)
      w(10)=425432555981648045046547654108115320672973233243977602393063
     & 5895887581199691.D0/18305426139000.D0 + 1.D0/2825761.D0*w(10)
      w(10)=w(10)*w(11)
      w(10)=110038436363762226882257374087742672334938006456266557344442
     & 5811514714414931.D0/147816316072425.D0 + w(10)
      w(10)=y*w(10)
      w(10)=138901925639726608291810416203782877446809928647490286791457
     & 4986252074684457.D0/363855547255200.D0 + w(10)
      w(10)=y*w(10)
      w(10)=138801070243021609143842471002603262261437058053650988025340
     & 176787327463601.D0/70856080254960.D0 + w(10)
      w(10)=w(10)*w(22)
      w(10)=196218732046247824360707801173253429498644647654638955590064
     & 080886933.D0/9880986024000.D0 + w(10)
      w(10)=y*w(10)
      w(10)=785613217277155812356682925597498685901038926899692407869753
     & 493832759.D0/76989349437000.D0 + w(10)
      w(10)=y*w(10)
      w(10)=101574760366720056293647780285087765267443939737246614231861
     & 4426189209.D0/193573221441600.D0 + w(10)
      w(10)=y*w(10)
      w(10)=409996179535578476387757652348819166687533328379789727296280
     & 75435553.D0/15182213446400.D0 + w(10)
      w(10)=y*w(10)
      w(10)=921494978405921087074946155511935838598887162189292250702966
     & 3043113.D0/6624965867520.D0 + w(10)
      w(10)=y*w(10)
      w(10)=359533835399538971844339874904073388748214025926185314396621
     & 7794993.D0/5014012253250.D0 + w(10)
      w(10)=y*w(10)
      w(10)=2314303176309453968531405159740037524040308601569350572561.D
     & 0/11553023625.D0 + 1.D0/1847042.D0*w(10)
      w(10)=y*w(10)
      w(10)=69685384121205009168421245497719619098448783512477816925879.
     & D0/673511654200.D0 + w(10)
      w(10)=y*w(10)
      w(10)=34253745058662436892158232010288249743721978713159431.D0/
     & 15616500900.D0 + 1.D0/24389.D0*w(10)
      w(10)=w(10)*w(23)
      w(10)=2449111080526528998970152115592624644956910184481883.D0/
     & 23728968900.D0 + w(10)
      w(10)=y*w(10)
      w(10)=5015561225502933012860100637958336424407185295647.D0/
     & 7595781050.D0 + 1.D0/81.D0*w(10)
      w(10)=w(10)*w(8)
      w(10)=9341577228077414745533038877114440882447087073209.D0/
     & 81800719000.D0 + w(10)
      w(4)=w(10)*w(4)
      w(4)=27954193776867913017339035182003798463014571821.D0/
     & 11779303536.D0 + w(4)
      w(4)=y*w(4)
      w(4)=4645670566124071228678966116238070281879459871.D0/3762833074.
     & D0 + w(4)
      w(4)=y*w(4)
      w(4)=24819861697527890192043543891573898385417.D0/469464996.D0 + 
     & 1.D0/12167.D0*w(4)
      w(4)=y*w(4)
      w(4)=61815357681413681445211561759414812970691.D0/2240628390.D0
     &  + w(4)
      w(4)=y*w(4)
      w(4)=1678555127412704382591340069566615504169.D0/116396280.D0 + 
     & w(4)
      w(4)=w(4)*w(28)
      w(4)=918627770495124737509920548928345308243.D0/608170563.D0 + 
     & w(4)
      w(4)=y*w(4)
      w(4)=4673193630368451209379626135360653.D0/40432392.D0 + 1.D0/
     & 6859.D0*w(4)
      w(4)=y*w(4)
      w(4)=139360330690584984086420502984266371.D0/2291168880.D0 + w(4)
      w(4)=y*w(4)
      w(4)=296239687203990750215600977847.D0/45302400.D0 + 1.D0/4913.D0
     & *w(4)
      w(4)=w(4)*w(7)
      w(4)=1673089582581267640765421779.D0/1926925.D0 + w(4)
      w(4)=w(4)*w(6)
      w(4)=807109909486818476482400599.D0/3468465.D0 + w(4)
      w(4)=y*w(4)
      w(4)=1632695223787025936456488919.D0/12882870.D0 + w(4)
      w(4)=y*w(4)
      w(4)=985117100091854650609.D0/30492.D0 + 1.D0/2197.D0*w(4)
      w(4)=w(4)*w(6)
      w(4)=2208945006690565547233.D0/232925.D0 + w(4)
      w(4)=w(4)*w(25)
      w(4)=623319207738846559.D0/12600.D0 + w(4)
      w(4)=w(4)*w(11)
      w(4)=2030622472831957.D0/945.D0 + w(4)
      w(4)=w(4)*w(20)
      w(4)=6612599956725053.D0/35840.D0 + w(4)
      w(4)=y*w(4)
      w(4)=1570806046931.D0/10.D0 + w(4)
      w(4)=w(4)*w(12)
      w(4)=16800269251.D0/40.D0 + w(4)
      w(4)=y*w(4)
      w(4)=19938943633.D0/40.D0 + w(4)
      w(4)=y*w(4)
      w(4)=14113669.D0/32.D0 + 1.D0/1125.D0*w(4)
      w(4)=w(4)*w(6)
      w(4)=3637967.D0/9.D0 + w(4)
      w(4)=y*w(4)
      w(4)=497.D0/2.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 215 + 1.D0/216.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 587.D0/6.D0 + w(4)
      w(10)= - 143145370434312981628659151.D0/3116.D0 - 
     & 7982065782750398919244867.D0/175.D0*y
      w(10)=w(10)*w(12)
      w(10)= - 118537707625838640342382171.D0/878712.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 230694761661631613940211.D0/842099.D0 + 1.D0/496.D0*
     & w(10)
      w(10)=w(10)*w(26)
      w(10)= - 260893034708209430098621.D0/44434160.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 1131257709883699171150153.D0/191260080.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 15246163724630110712193073.D0/2558634848.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 8544343599612361346311.D0/203527772.D0 + 1.D0/143.D0*
     & w(10)
      w(10)=y*w(10)
      w(10)= - 1026869316766943469683.D0/30273572784.D0 + 1.D0/1247.D0*
     & w(10)
      w(10)=y*w(10)
      w(10)= - 14288324349436470084253.D0/418063624160.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 222689211933741420943.D0/6466183360.D0 + w(10)
      w(10)=y*w(10)
      w(10)= - 79945018683816636037.D0/2303577822.D0 + w(10)
      w(8)=w(10)*w(8)
      w(8)= - 117886744191443221.D0/10111952.D0 + w(8)
      w(8)=y*w(8)
      w(8)= - 10811715314500785613.D0/920187632.D0 + w(8)
      w(8)=w(8)*w(27)
      w(8)= - 15044120913871575473.D0/47004179040.D0 + w(8)
      w(2)=w(8)*w(2)
      w(2)= - 1712330390453570221.D0/122080298340.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 49907406175309545601.D0/3529864626288.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 3738211744386764597.D0/262280962944.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 134527222424307023261.D0/9362635586304.D0 + w(2)
      w(2)=w(2)*w(6)
      w(2)= - 9958583227330081.D0/1374917115.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 148793426010348137.D0/20375142480.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 2730479479866488731.D0/370827593136.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 6547517425177843.D0/881872992.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 13794179774580661.D0/1842484644.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 17173027975934143.D0/2274672400.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 40988963917547.D0/5383840.D0 + w(2)
      w(2)=w(2)*w(9)
      w(2)= - 9416540945275.D0/91986752.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 3264270603883.D0/31620446.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 3326277868439.D0/31951920.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 16770765752173.D0/159759600.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 275392211273903.D0/2601799200.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 2313729879169.D0/21681660.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 439323479363.D0/4084080.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 6487092431.D0/59840.D0 + w(2)
      w(2)=w(2)*w(16)
      w(2)= - 20578850309.D0/3203200.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 1457094629.D0/225225.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 71170027871.D0/10930920.D0 + w(2)
      w(2)=w(2)*w(20)
      w(2)= - 29127107.D0/40040.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 4010172877.D0/5488560.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 2513259433.D0/3430350.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 1319333.D0/1800.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 99443329.D0/136080.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 153428651.D0/211680.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 673769.D0/945.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 1867111.D0/2700.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 39187.D0/60.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 125615.D0/216.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 3899.D0/9.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 6895.D0/18.D0 + w(2)
      w(2)=w(2)*w(7)
      w(2)= - 454.D0/9.D0 + w(2)
      w(2)=y*w(2)
      w(2)=59.D0/9.D0 + w(2)
      w(2)=z3*w(2)
*
      RED1 = w(2) + w(3) + 1.D0/12.D0*w(4) + w(5)
*
      RETURN
      END      
      REAL*8 FUNCTION ANSREG1(z,nf,as,LL)
*     -----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 19, 2024
*     regular part of aqqQNS3 z \in [0,1/2] 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(44),G
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
*      WRITE(6,*) 'ANSREG1::LL=', LL !!
      w(1)=Log(z)
      w(2)=z + 1
      w(3)=w(2)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(3)=w(3)*z
      w(3)=w(3) + 1
      w(4)=w(3)*z
      w(5)=w(4) + 1
      w(5)=w(5)*z
      w(5)=w(5) + 1
      w(6)=832.D0/3.D0*LL
      w(6)=w(5)*w(6)
      w(7)=1.D0/49.D0*z
      w(8)= - 68239769087361229294201.D0 - 220350976406867345205044.D0*
     & z
      w(8)=w(8)*w(7)
      w(8)= - 4486574670402187298201.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 1403362286978859425524.D0 + w(8)
      w(9)=1.D0/47.D0*z
      w(8)=w(8)*w(9)
      w(8)= - 95228636650661002948.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 30096814915289926592.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 94987928479155706648.D0 + w(8)
      w(10)=1.D0/11.D0*z
      w(8)=w(8)*w(10)
      w(8)= - 19311013931662851724.D0/7.D0 + w(8)
      w(11)=1.D0/43.D0*z
      w(8)=w(8)*w(11)
      w(8)= - 1402011864831473132.D0/7.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 452958238450897618.D0/7.D0 + w(8)
      w(12)=1.D0/41.D0*z
      w(8)=w(8)*w(12)
      w(8)= - 34099986116334682.D0/7.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 45128898908444.D0/7.D0 + 1.D0/247.D0*w(8)
      w(8)=z*w(8)
      w(8)= - 33999627994450732.D0/1729.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 11251242425094968.D0/1729.D0 + w(8)
      w(13)=1.D0/37.D0*z
      w(8)=w(8)*w(13)
      w(8)= - 916048612991296.D0/1729.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 43867109413052.D0/247.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 10043257636065956.D0/19019.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 6825096168297713.D0/38038.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 20015887758873337.D0/38038.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 3449559875771344.D0/19019.D0 + w(8)
      w(14)=1.D0/31.D0*z
      w(8)=w(8)*w(14)
      w(8)= - 321623470229696.D0/19019.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 112551905761804.D0/19019.D0 + w(8)
      w(15)=1.D0/29.D0*z
      w(8)=w(8)*w(15)
      w(8)= - 11045703049324.D0/19019.D0 + w(8)
      w(16)=1.D0/3.D0*z
      w(8)=w(8)*w(16)
      w(8)= - 1309443748192.D0/19019.D0 + w(8)
      w(17)=1.D0/15.D0*z
      w(8)=w(8)*w(17)
      w(8)= - 244390754672.D0/19019.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 442146282874.D0/95095.D0 + w(8)
      w(18)=1.D0/5.D0*z
      w(8)=w(8)*w(18)
      w(8)= - 48646940170.D0/19019.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 17932025920.D0/19019.D0 + w(8)
      w(19)=1.D0/23.D0*z
      w(8)=w(8)*w(19)
      w(8)= - 2104143952.D0/19019.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 791368876.D0/19019.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 10460716492.D0/95095.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 12064633904.D0/285285.D0 + w(8)
      w(8)=w(8)*w(16)
      w(8)= - 182355616.D0/5005.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 646389761.D0/45045.D0 + w(8)
      w(20)=1.D0/17.D0*z
      w(8)=w(8)*w(20)
      w(8)= - 95849753.D0/45045.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 38782888.D0/45045.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 13580504.D0/6435.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 39658972.D0/45045.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 7242476.D0/3465.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 3130264.D0/3465.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 650824.D0/315.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 97802.D0/105.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 213818.D0/105.D0 + w(8)
      w(8)=w(8)*w(16)
      w(8)= - 11288.D0/35.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 5992.D0/9.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 5092.D0/15.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 1940.D0/3.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 3320.D0/9.D0 + w(8)
      w(8)=z*w(8)
      w(8)= - 5464.D0/9.D0 + w(8)
      w(8)=w(8)*w(16)
      w(8)= - 392 + w(8)
      w(8)=z*w(8)
      w(21)=z - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(21)=w(21)*z
      w(21)=w(21) - 1
      w(21)=w(21)*z
      w(21)=w(21) + 1
      w(22)=2.D0/3.D0*z
      w(23)= - w(21)*w(22)
      w(23)= - 7 + w(23)
      w(23)=z*w(23)
      w(23)= - 23.D0/3.D0 + w(23)
      w(23)=w(1)*w(23)
      w(8)=16.D0/3.D0*w(23) - w(6) - 6328.D0/27.D0 + w(8)
      w(8)=w(1)*w(8)
      w(23)=1.D0/343.D0*z
      w(24)= - 59891519246204288545559035271957.D0/8.D0 - 10359886093864
     & 74319156648744903741.D0/25.D0*z
      w(24)=w(24)*w(23)
      w(24)= - 406707854253070895237694394870253.D0/3384.D0 + w(24)
      w(25)=1.D0/4.D0*z
      w(24)=w(24)*w(25)
      w(24)= - 2376620254463830350823839478118.D0/423.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 4187098187940934591461460717.D0/270.D0 + 1.D0/1927.D0*
     & w(24)
      w(24)=w(24)*w(19)
      w(24)= - 642338739722554652054431477.D0/4920.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1276545223762845224103948531259.D0/1904040.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 448138178249654352957885616973.D0/3332070.D0 + w(24)
      w(24)=w(24)*w(11)
      w(24)= - 3517406932921564226091948562.D0/226935.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 81917752998851964880729365053.D0/25416720.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 232907333589662448959258173.D0/15120.D0 + w(24)
      w(26)=1.D0/19.D0*z
      w(24)=w(24)*w(26)
      w(24)= - 330483066015873260426977.D0/1890.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 535111645942838493137257372.D0/664335.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 50439707547051416481827161.D0/279720.D0 + w(24)
      w(24)=w(24)*w(13)
      w(24)= - 54486833023720500443639.D0/2520.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 527966913540316847006.D0/105.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1932011021596172014271.D0/90.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 26647667749600862593.D0/960.D0 + 1.D0/187.D0*w(24)
      w(27)=1.D0/7.D0*z
      w(24)=w(24)*w(27)
      w(24)= - 15983540033443961174873.D0/982080.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 879990797437847527199.D0/214830.D0 + w(24)
      w(24)=w(24)*w(14)
      w(24)= - 34887870536519064221.D0/66990.D0 + w(24)
      w(24)=w(24)*w(20)
      w(24)= - 1507165143381492697.D0/187572.D0 + w(24)
      w(24)=w(24)*w(15)
      w(24)= - 1036466574723820729.D0/989604.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 141849629092665629.D0/494802.D0 + w(24)
      w(24)=w(24)*w(26)
      w(24)= - 225222997120774841.D0/4123350.D0 + w(24)
      w(28)=1.D0/13.D0*z
      w(24)=w(24)*w(28)
      w(24)= - 376715489917739407.D0/313374600.D0 + w(24)
      w(24)=w(24)*w(18)
      w(24)= - 239901154511071013.D0/288304632.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 17970118011233503.D0/72076158.D0 + w(24)
      w(24)=w(24)*w(19)
      w(24)= - 112179204770903.D0/3133746.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 32071184884859.D0/2848860.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1915476723590809.D0/54128340.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 316969212127081.D0/27064170.D0 + w(24)
      w(24)=w(24)*w(28)
      w(24)= - 1141995672887.D0/424830.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 41500139919097.D0/44182320.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 135102939971.D0/50960.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1876605737.D0/1911.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 108093888911.D0/41405.D0 + w(24)
      w(24)=w(24)*w(16)
      w(24)= - 4062876053.D0/11830.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 5927127287.D0/6930.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1254614561.D0/3465.D0 + w(24)
      w(24)=w(24)*w(10)
      w(24)= - 24014957.D0/315.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 8806757.D0/252.D0 + w(24)
      w(24)=w(24)*w(16)
      w(24)= - 4869593.D0/196.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 5515906.D0/441.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1089602.D0/45.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 205859.D0/15.D0 + w(24)
      w(24)=w(24)*w(18)
      w(24)= - 4721 + w(24)
      w(24)=z*w(24)
      w(24)= - 28876.D0/9.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 44588.D0/9.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 14240.D0/3.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 36344.D0/9.D0 + w(24)
      w(29)= - 35264663766436512300791.D0 - 35326644656521432234919.D0*
     & z
      w(29)=w(29)*w(7)
      w(29)= - 718396284902600210759.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 717078663939910586309.D0 + w(29)
      w(29)=w(29)*w(9)
      w(29)= - 15228361882737299947.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 15199108503362137147.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 15169205048889748507.D0 + w(29)
      w(29)=w(29)*w(10)
      w(29)= - 1376238360950188937.D0 + w(29)
      w(29)=w(29)*w(11)
      w(29)= - 31939381940025659.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 223101517325332013.D0/7.D0 + w(29)
      w(29)=w(29)*w(12)
      w(29)= - 5429653567897093.D0/7.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 5417510541858313.D0/7.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 5405056156177513.D0/7.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 5392274023505113.D0/7.D0 + w(29)
      w(29)=w(29)*w(13)
      w(29)= - 145382335886149.D0/7.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 145017680449549.D0/7.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 144642606286189.D0/7.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 144256500529789.D0/7.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 1582445640588479.D0/77.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 1577933029560554.D0/77.D0 + w(29)
      w(29)=w(29)*w(14)
      w(29)= - 50750801626934.D0/77.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 50595528989414.D0/77.D0 + w(29)
      w(29)=w(29)*w(15)
      w(29)= - 1739134555966.D0/77.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 1733397882166.D0/77.D0 + w(29)
      w(30)=1.D0/9.D0*z
      w(29)=w(29)*w(30)
      w(29)= - 191938748774.D0/77.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 191252309174.D0/77.D0 + w(29)
      w(29)=w(29)*w(18)
      w(29)= - 38107682398.D0/77.D0 + w(29)
      w(29)=w(29)*w(16)
      w(29)= - 12652984606.D0/77.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 60875618.D0/77.D0 + 1.D0/207.D0*w(29)
      w(29)=z*w(29)
      w(29)= - 545529122.D0/693.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 1629197126.D0/2079.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 1621437374.D0/2079.D0 + w(29)
      w(29)=w(29)*w(26)
      w(29)= - 84908906.D0/2079.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 253365358.D0/6237.D0 + w(29)
      w(29)=w(29)*w(20)
      w(29)= - 14819054.D0/6237.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 14728964.D0/6237.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 14632868.D0/6237.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 14529908.D0/6237.D0 + w(29)
      w(29)=w(29)*w(28)
      w(29)= - 1109156.D0/6237.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 1099916.D0/6237.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 99076.D0/567.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 98068.D0/567.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 10772.D0/63.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 3544.D0/21.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 1496.D0/9.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 4408.D0/27.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 4312.D0/27.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 4192.D0/27.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 448.D0/3.D0 + w(29)
      w(29)=z*w(29)
      w(29)= - 1264.D0/9.D0 + w(29)
      w(29)=z*w(29)
      w(31)=w(5)*LL
      w(29)= - 32*w(31) - 368.D0/3.D0 + w(29)
      w(29)=LL*w(29)
      w(8)=1.D0/9.D0*w(8) + 1.D0/27.D0*w(24) + w(29)
      w(8)=w(1)*w(8)
      w(24)= - 164419643964810970413430887850696563449732161.D0/5.D0 - 
     & 32977827561187372442452103182078351864688909.D0*z
      w(24)=w(24)*w(23)
      w(24)= - 477960608210155291356012998631664167820327.D0/5.D0 + 
     & w(24)
      w(24)=z*w(24)
      w(24)= - 476534302606798836603964302175793656802827.D0/5.D0 + 
     & w(24)
      w(29)=1.D0/2209.D0*z
      w(24)=w(24)*w(29)
      w(24)= - 215064457221536822777748331392549461803.D0/5.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 214390455297757931784958357351297781803.D0/5.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 213701339477013106306652711912816194603.D0/5.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 212996416438118266509961452122236714603.D0/5.D0 + w(24)
      w(32)=1.D0/1849.D0*z
      w(24)=w(24)*w(32)
      w(24)= - 114805269895641228219395951702897347.D0/5.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 800839851969436475084242794240121429.D0/35.D0 + w(24)
      w(33)=1.D0/1681.D0*z
      w(24)=w(24)*w(33)
      w(24)= - 474701912289674277242026715602309.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 472953954311986541310544694545909.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 471160707686012353315576527505909.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 469319763453385551106212309745909.D0/35.D0 + w(24)
      w(34)=1.D0/1369.D0*z
      w(24)=w(24)*w(34)
      w(24)= - 341437922004547824132562449661.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 340017631117080499099042089661.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 338556287646749293696470288061.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 337051448103363653564239728061.D0/35.D0 + w(24)
      w(24)=w(24)*w(10)
      w(24)= - 30500040379949173166699262551.D0/35.D0 + w(24)
      w(24)=w(24)*w(25)
      w(24)= - 7588644415962303351549756419.D0/35.D0 + w(24)
      w(35)=1.D0/961.D0*z
      w(24)=w(24)*w(35)
      w(24)= - 7857534068216837683760579.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 7817135557717828510558979.D0/35.D0 + w(24)
      w(36)=1.D0/841.D0*z
      w(24)=w(24)*w(36)
      w(24)= - 9245332281721691383019.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1313402139621031724717.D0/5.D0 + w(24)
      w(24)=w(24)*w(30)
      w(24)= - 145085093605435649413.D0/5.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 144203472437290169413.D0/5.D0 + w(24)
      w(37)=1.D0/25.D0*z
      w(24)=w(24)*w(37)
      w(24)= - 1146288109118907497.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1138637448186725957.D0 + w(24)
      w(38)=1.D0/529.D0*z
      w(24)=w(24)*w(38)
      w(24)= - 2137331351317493.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 23336828310250183.D0/11.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 162082026302402161.D0/77.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 803705877724595621.D0/385.D0 + w(24)
      w(39)=1.D0/361.D0*z
      w(24)=w(24)*w(39)
      w(24)= - 2206762500292061.D0/385.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 2186081345397661.D0/385.D0 + w(24)
      w(40)=1.D0/289.D0*z
      w(24)=w(24)*w(40)
      w(24)= - 7488424893949.D0/385.D0 + w(24)
      w(24)=w(24)*w(25)
      w(24)= - 1851923304481.D0/385.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 366071757101.D0/77.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 361441954781.D0/77.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 2109142229.D0/77.D0 + 1.D0/169.D0*w(24)
      w(24)=z*w(24)
      w(24)= - 2077028609.D0/77.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 185626339.D0/7.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 910498751.D0/35.D0 + w(24)
      w(24)=w(24)*w(30)
      w(24)= - 98980039.D0/35.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 96506764.D0/35.D0 + w(24)
      w(24)=w(24)*w(27)
      w(24)= - 1911436.D0/5.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1843036.D0/5.D0 + w(24)
      w(24)=w(24)*w(18)
      w(24)= - 70396 + w(24)
      w(24)=z*w(24)
      w(24)= - 66160 + w(24)
      w(24)=w(24)*w(30)
      w(24)= - 6704 + w(24)
      w(24)=z*w(24)
      w(24)= - 17024.D0/3.D0 + w(24)
      w(24)=w(24)*w(16)
      w(41)= - 11944353871985387085131.D0 - 12006334762070307019259.D0*
     & z
      w(41)=w(41)*w(7)
      w(41)= - 242471593179107859419.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 241153972216418234969.D0 + w(41)
      w(41)=w(41)*w(9)
      w(41)= - 5102304612024696727.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 5073051232649533927.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 5043147778177145287.D0 + w(41)
      w(41)=w(41)*w(10)
      w(41)= - 455687699976315917.D0 + w(41)
      w(41)=w(41)*w(11)
      w(41)= - 10531227033656519.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 73244432980748033.D0/7.D0 + w(41)
      w(41)=w(41)*w(12)
      w(41)= - 1774602730224313.D0/7.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 1762459704185533.D0/7.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 1750005318504733.D0/7.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 1737223185832333.D0/7.D0 + w(41)
      w(41)=w(41)*w(13)
      w(41)= - 46597178111209.D0/7.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 46232522674609.D0/7.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 45857448511249.D0/7.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 45471342754849.D0/7.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 495808905064139.D0/77.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 491296294036214.D0/77.D0 + w(41)
      w(41)=w(41)*w(14)
      w(41)= - 15698003706794.D0/77.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 15542731069274.D0/77.D0 + w(41)
      w(41)=w(41)*w(15)
      w(41)= - 530417386306.D0/77.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 524680712506.D0/77.D0 + w(41)
      w(41)=w(41)*w(30)
      w(41)= - 57636841034.D0/77.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 56950401434.D0/77.D0 + w(41)
      w(41)=w(41)*w(37)
      w(41)= - 2249460170.D0/77.D0 + w(41)
      w(41)=w(41)*w(16)
      w(41)= - 739904818.D0/77.D0 + w(41)
      w(41)=w(41)*w(19)
      w(41)= - 31719934.D0/77.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 31249646.D0/77.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 92270890.D0/231.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 453594698.D0/1155.D0 + w(41)
      w(41)=w(41)*w(26)
      w(41)= - 23443502.D0/1155.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 68969146.D0/3465.D0 + w(41)
      w(41)=w(41)*w(20)
      w(41)= - 3972218.D0/3465.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 3882128.D0/3465.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 3786032.D0/3465.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 3683072.D0/3465.D0 + w(41)
      w(41)=w(41)*w(28)
      w(41)= - 274784.D0/3465.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 265544.D0/3465.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 23224.D0/315.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 22216.D0/315.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 2344.D0/35.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 2204.D0/35.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 292.D0/5.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 796.D0/15.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 140.D0/3.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 116.D0/3.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 28 + w(41)
      w(41)=w(41)*w(16)
      w(41)= - 4 + w(41)
      w(41)=z*w(41)
      w(41)=20.D0/3.D0 + w(41)
      w(41)=LL*w(41)
      w(24)=16*w(41) - 1088 + w(24)
      w(24)=LL*w(24)
      w(41)=139279196412485297168279595698026195969035739843727.D0/64.D0
     &  - 92256703182777045165066183972531293558362054440503.D0/5.D0*z
      w(41)=z*w(41)
      w(41)= - 2562516427912179908093498003039596318120330249.D0/3008.D0
     &  + 1.D0/21428925.D0*w(41)
      w(41)=z*w(41)
      w(41)=136906148404255449121402459235232464488129151.D0/324300.D0
     &  + 1.D0/217.D0*w(41)
      w(41)=z*w(41)
      w(41)= - 35320163675684697722489232155360188727931139.D0/
     & 381811500.D0 + 1.D0/41971.D0*w(41)
      w(41)=z*w(41)
      w(41)=5635640803919195310105146215356193663119551.D0/630819000.D0
     &  + w(41)
      w(41)=w(41)*w(15)
      w(41)= - 5702587899397025052639868346283944663913893.D0/
     & 1808347800.D0 + w(41)
      w(41)=z*w(41)
      w(41)=344448518624510423692255446200358926473663747.D0/
     & 1284831111900.D0 + w(41)
      w(41)=w(41)*w(32)
      w(41)= - 137639919895050514626474310401467863482593.D0/
     & 81671435020.D0 + w(41)
      w(41)=z*w(41)
      w(41)=1502794921589772892585768559872754715145171.D0/
     & 12250715253000.D0 + w(41)
      w(41)=w(41)*w(33)
      w(41)= - 147933604109266298931750144308128785317.D0/149398966500.D
     & 0 + w(41)
      w(41)=w(41)*w(16)
      w(41)=3737024917646235897990323176614554807.D0/189238690900.D0 + 
     & w(41)
      w(41)=w(41)*w(28)
      w(41)= - 1579414818958530707051590228466983924973.D0/
     & 63016484069700.D0 + w(41)
      w(41)=z*w(41)
      w(41)=22970538294187376109992636360738709101.D0/19899942337800.D0
     &  + w(41)
      w(41)=z*w(41)
      w(41)= - 3736443233368709539741048779509243.D0/2689181397000.D0
     &  + 1.D0/17797.D0*w(41)
      w(41)=w(41)*w(19)
      w(41)=187699530324643938581837822266129.D0/99051514789500.D0 + 
     & w(41)
      w(41)=z*w(41)
      w(41)= - 9046601297375258256066588414251507.D0/151878989343900.D0
     &  + w(41)
      w(41)=z*w(41)
      w(41)=319937220776707490262597728649857.D0/343067834753280.D0 + 
     & w(41)
      w(42)=1.D0/121.D0*z
      w(41)=w(41)*w(42)
      w(41)= - 25792977158171029144058284596502523.D0/53175514386758400.
     & D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 204015535238978061507094840769.D0/346194755122125.D0 + 
     & w(41)
      w(41)=w(41)*w(14)
      w(41)= - 5311378841779896693612099607523.D0/344753777993625.D0 + 
     & w(41)
      w(43)=1.D0/2.D0*z
      w(41)=w(41)*w(43)
      w(41)= - 77507108718671228798518223.D0/518426733825.D0 + w(41)
      w(41)=w(41)*w(15)
      w(41)= - 77132564665103276919114157.D0/295151867010.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 17830443369364622209694791.D0/1744079214150.D0 + w(41)
      w(44)=1.D0/27.D0*z
      w(41)=w(41)*w(44)
      w(41)= - 3313844918578533192240119.D0/348309211250.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 756248390946885989123507.D0/1311786075600.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 33826358404031943085320181.D0/3620529568656.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 3913233075569964398437223.D0/4978228156902.D0 + w(41)
      w(41)=w(41)*w(19)
      w(41)= - 59522709713151259307.D0/149375226.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 37543671640565343233.D0/855512658.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 58115851886960977643.D0/148852935.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 1891428917816880935323.D0/34831586790.D0 + w(41)
      w(41)=w(41)*w(26)
      w(41)= - 32981085326289966221.D0/1640268630.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 1830145705834128707.D0/530187840.D0 + w(41)
      w(41)=w(41)*w(20)
      w(41)= - 69967849885252469.D0/60540480.D0 + w(41)
      w(41)=w(41)*w(43)
      w(41)= - 799068609502102.D0/6621615.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 154072556678918.D0/273273.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 1397195162782219.D0/9837828.D0 + w(41)
      w(41)=w(41)*w(28)
      w(41)= - 27056565970523.D0/640332.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 162622008712.D0/12705.D0 + w(41)
      w(41)=w(41)*w(10)
      w(41)= - 14838428288.D0/3969.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 86726467517.D0/63504.D0 + w(41)
      w(41)=w(41)*w(37)
      w(41)= - 799070941.D0/5488.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 1668266578.D0/25725.D0 + w(41)
      w(41)=w(41)*w(16)
      w(41)= - 5938814.D0/125.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 11940899.D0/450.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 859969.D0/18.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 110356.D0/3.D0 + w(41)
      w(41)=z*w(41)
      w(41)= - 60404 + w(41)
      w(41)=z*w(41)
      w(41)= - 50264 + w(41)
      w(41)=z*w(41)
      w(41)= - 4664.D0/3.D0 + w(41)
      w(24)=1.D0/27.D0*w(41) + w(24)
      w(8)=1.D0/3.D0*w(24) + w(8)
      w(8)=w(1)*w(8)
      w(24)=8898324234708784680131.D0/7.D0 - 1020930991203502332913.D0/
     & 17.D0*z
      w(24)=w(24)*w(27)
      w(24)= - 162983405644019657017.D0/17.D0 + w(24)
      w(24)=z*w(24)
      w(24)=179003303547944292269.D0 + w(24)
      w(24)=w(24)*w(9)
      w(24)= - 3846950576668825931.D0/17.D0 + w(24)
      w(24)=z*w(24)
      w(24)=3751002719715501667.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 4242004511731263491.D0/17.D0 + w(24)
      w(24)=w(24)*w(10)
      w(24)=2348723459709475319.D0/7.D0 + w(24)
      w(24)=w(24)*w(11)
      w(24)= - 9839885894590267.D0/17.D0 + w(24)
      w(24)=z*w(24)
      w(24)=53689637185866413.D0/7.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 20171006767139.D0/17.D0 + 1.D0/533.D0*w(24)
      w(24)=z*w(24)
      w(24)=1285660128538033.D0/91.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 285508183021967.D0/221.D0 + w(24)
      w(24)=z*w(24)
      w(24)=1260586354716433.D0/91.D0 + w(24)
      w(24)=w(24)*w(13)
      w(24)= - 8377233637931.D0/221.D0 + w(24)
      w(24)=z*w(24)
      w(24)=567041033883173.D0/1547.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 99800592656761.D0/2431.D0 + w(24)
      w(24)=z*w(24)
      w(24)=358596571883489.D0/1001.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 107875859622436.D0/2431.D0 + w(24)
      w(24)=z*w(24)
      w(24)=349776702042434.D0/1001.D0 + w(24)
      w(24)=w(24)*w(14)
      w(24)= - 3755630816386.D0/2431.D0 + w(24)
      w(24)=z*w(24)
      w(24)=10980106931774.D0/1001.D0 + w(24)
      w(24)=w(24)*w(15)
      w(24)= - 139604866594.D0/2431.D0 + w(24)
      w(24)=z*w(24)
      w(24)=367449030206.D0/1001.D0 + w(24)
      w(24)=w(24)*w(30)
      w(24)= - 16708163282.D0/2431.D0 + w(24)
      w(24)=w(24)*w(16)
      w(24)=13164388546.D0/1001.D0 + w(24)
      w(24)=w(24)*w(37)
      w(24)= - 239856250.D0/2431.D0 + w(24)
      w(24)=z*w(24)
      w(24)=1522035290.D0/3003.D0 + w(24)
      w(24)=w(24)*w(19)
      w(24)= - 11227438.D0/2431.D0 + w(24)
      w(24)=z*w(24)
      w(24)=105744634.D0/5005.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 181372778.D0/36465.D0 + w(24)
      w(24)=w(24)*w(26)
      w(24)=47727926.D0/45045.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 30864874.D0/109395.D0 + w(24)
      w(24)=z*w(24)
      w(24)=766950362.D0/765765.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1959136.D0/6435.D0 + w(24)
      w(24)=z*w(24)
      w(24)=42190108.D0/45045.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 2117668.D0/6435.D0 + w(24)
      w(24)=z*w(24)
      w(24)=38870044.D0/45045.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 176476.D0/495.D0 + w(24)
      w(24)=z*w(24)
      w(24)=2694724.D0/3465.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 17428.D0/45.D0 + w(24)
      w(24)=z*w(24)
      w(24)=213172.D0/315.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 6328.D0/15.D0 + w(24)
      w(24)=z*w(24)
      w(24)=58064.D0/105.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 6896.D0/15.D0 + w(24)
      w(24)=w(24)*w(16)
      w(24)=656.D0/5.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1480.D0/9.D0 + w(24)
      w(24)=z*w(24)
      w(24)=520.D0/9.D0 + w(24)
      w(24)=z*w(24)
      w(24)= - 1256.D0/9.D0 + w(24)
      w(24)=w(24)*w(16)
      w(24)=140 + w(24)
      w(24)=z*w(24)
      w(30)=10*z
      w(41)=17 + w(30)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=17 + w(41)
      w(41)=z*w(41)
      w(41)=10 + w(41)
      w(41)=z*w(41)
      w(41)=14 + w(41)
      w(22)=w(41)*w(22)
      w(22)=25 + w(22)
      w(41)=8.D0/3.D0*w(1)
      w(22)=w(22)*w(41)
      w(6)=w(22) + w(6) - 2660.D0/27.D0 + w(24)
      w(6)=w(1)*w(6)
      w(22)=1182054044412164445230796028012776687049693.D0 - 37493634743
     & 5651745886440105510076973874329.D0*z
      w(22)=w(22)*w(23)
      w(22)= - 1134694831146471852325556716889758904251.D0 + w(22)
      w(22)=z*w(22)
      w(22)=81815262774565066895687290409711366461373.D0/23.D0 + w(22)
      w(22)=w(22)*w(29)
      w(22)= - 12226414128331734392525110598901959837.D0/23.D0 + w(22)
      w(22)=z*w(22)
      w(22)=38221015295999835350774406037478934637.D0/23.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 12613679370683995824799566025110576877.D0/23.D0 + w(22)
      w(22)=z*w(22)
      w(22)=276044111483699166228334929467916511339.D0/161.D0 + w(22)
      w(22)=w(22)*w(32)
      w(22)= - 49108552143251993956784007752310211.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)=154009841988407950954234582254060571.D0/161.D0 + w(22)
      w(22)=w(22)*w(33)
      w(22)= - 29942090482860542886501661136491.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)=94499633033000241079341029867251.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 30577668226743555691673883544051.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)=97462020056009825516044173189251.D0/161.D0 + w(22)
      w(22)=w(22)*w(34)
      w(22)= - 22719502133643384502967742139.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)=73418611176770408315586227339.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 23007010275584557662595567499.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)=75711442610595471133743514949.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 23179007035418903411859278249.D0/161.D0 + w(22)
      w(22)=w(22)*w(25)
      w(22)=19518615637342700760318676511.D0/161.D0 + w(22)
      w(22)=w(22)*w(35)
      w(22)= - 6038389931832579809645951.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)=20944696746207483659104631.D0/161.D0 + w(22)
      w(22)=w(22)*w(36)
      w(22)= - 21413047715548656624973.D0/483.D0 + w(22)
      w(22)=z*w(22)
      w(22)=77046689546020494863573.D0/483.D0 + w(22)
      w(22)=w(22)*w(16)
      w(22)= - 2343734161316719584397.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)=8828183215052736250777.D0/161.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 18254658283083851141.D0/161.D0 + 1.D0/125.D0*w(22)
      w(22)=z*w(22)
      w(22)=72830729700232812821.D0/161.D0 + w(22)
      w(22)=w(22)*w(19)
      w(22)= - 33073496028696869.D0/7.D0 + w(22)
      w(22)=z*w(22)
      w(22)=709814735020144129.D0/35.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 155051704588259569.D0/35.D0 + w(22)
      w(22)=z*w(22)
      w(22)=731742436345631801.D0/35.D0 + w(22)
      w(22)=w(22)*w(39)
      w(22)= - 389794556356721.D0/35.D0 + w(22)
      w(22)=w(22)*w(44)
      w(22)=77351121617233.D0/35.D0 + w(22)
      w(22)=w(22)*w(40)
      w(22)= - 25808265919.D0/21.D0 + w(22)
      w(22)=w(22)*w(25)
      w(22)=123955153363.D0/63.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 75678966359.D0/315.D0 + w(22)
      w(22)=w(22)*w(16)
      w(22)=23567568637.D0/35.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 25423201.D0/945.D0 + 1.D0/1859.D0*w(22)
      w(22)=z*w(22)
      w(22)=1281932849.D0/3465.D0 + w(22)
      w(22)=w(22)*w(10)
      w(22)= - 1333081.D0/2835.D0 + w(22)
      w(22)=z*w(22)
      w(22)=96212537.D0/2835.D0 + w(22)
      w(22)=z*w(22)
      w(22)=225059.D0/105.D0 + w(22)
      w(22)=z*w(22)
      w(22)=3476756.D0/105.D0 + w(22)
      w(22)=w(22)*w(27)
      w(22)=3772.D0/5.D0 + w(22)
      w(22)=z*w(22)
      w(22)=60704.D0/15.D0 + w(22)
      w(22)=w(22)*w(18)
      w(22)=1564.D0/9.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 56.D0/3.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 18608.D0/9.D0 + w(22)
      w(22)=z*w(22)
      w(22)=65776.D0/9.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 10960.D0/9.D0 + w(22)
      w(23)=2*z
      w(24)=z*w(21)
      w(24)= - 173.D0/9.D0 + w(24)
      w(24)=w(24)*w(23)
      w(24)= - 139.D0/9.D0 + w(24)
      w(24)=z2*w(24)
      w(6)=64.D0/5.D0*w(24) + 1.D0/9.D0*w(22) + 2*w(6)
      w(6)=z2*w(6)
      w(22)= - 81018177705009814936423.D0 + 78614441782775133052288.D0*
     & z
      w(22)=w(22)*w(7)
      w(22)=1582262929971172751977.D0 + w(22)
      w(22)=z*w(22)
      w(22)= - 1633406978800012064752.D0 + w(22)
      w(9)=w(22)*w(9)
      w(9)=33173228715273082496.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 34309781601734268716.D0 + w(9)
      w(9)=z*w(9)
      w(9)=32657805404636735396.D0 + w(9)
      w(9)=w(9)*w(10)
      w(9)= - 21539154829865208352.D0/7.D0 + w(9)
      w(9)=w(9)*w(11)
      w(9)=475297802272587664.D0/7.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 493759538646547114.D0/7.D0 + w(9)
      w(9)=w(9)*w(12)
      w(9)=11386956752385014.D0/7.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 11860346313377264.D0/7.D0 + w(9)
      w(9)=z*w(9)
      w(9)=11169836550236864.D0/7.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 11668829130441164.D0/7.D0 + w(9)
      w(9)=w(9)*w(13)
      w(9)=42239091071756.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 309931052207072.D0/7.D0 + w(9)
      w(9)=w(9)*w(20)
      w(9)=187045927167536.D0/77.D0 + w(9)
      w(9)=w(9)*w(23)
      w(9)= - 6692247939924299.D0/1309.D0 + w(9)
      w(9)=z*w(9)
      w(9)=6204612457330949.D0/1309.D0 + w(9)
      w(9)=w(9)*w(14)
      w(9)= - 211573536741104.D0/1309.D0 + w(9)
      w(9)=z*w(9)
      w(9)=194793383431184.D0/1309.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 207006761951384.D0/1309.D0 + w(9)
      w(9)=z*w(9)
      w(9)=2172716396632.D0/1309.D0 + 1.D0/87.D0*w(9)
      w(9)=z*w(9)
      w(9)= - 6970529230096.D0/3927.D0 + w(9)
      w(9)=w(9)*w(17)
      w(9)=140060944048.D0/1309.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 754595228204.D0/6545.D0 + w(9)
      w(10)=1.D0/45.D0*z
      w(9)=w(9)*w(10)
      w(9)=2996451604.D0/1309.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 3258754736.D0/1309.D0 + w(9)
      w(9)=w(9)*w(19)
      w(9)=1122610480.D0/11781.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 411680968.D0/3927.D0 + w(9)
      w(9)=z*w(9)
      w(9)=5335955528.D0/58905.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 661935152.D0/6545.D0 + w(9)
      w(9)=w(9)*w(26)
      w(9)=793356496.D0/176715.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 903056086.D0/176715.D0 + w(9)
      w(9)=z*w(9)
      w(9)=43369322.D0/10395.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 7240784.D0/1485.D0 + w(9)
      w(9)=z*w(9)
      w(9)=39533072.D0/10395.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 47979752.D0/10395.D0 + w(9)
      w(9)=w(9)*w(28)
      w(9)=2688824.D0/10395.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 3457424.D0/10395.D0 + w(9)
      w(9)=z*w(9)
      w(9)=68272.D0/315.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 290356.D0/945.D0 + w(9)
      w(9)=z*w(9)
      w(9)=5668.D0/35.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 29264.D0/105.D0 + w(9)
      w(9)=z*w(9)
      w(9)=752.D0/9.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 11176.D0/45.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 152.D0/3.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 2096.D0/9.D0 + w(9)
      w(9)=z*w(9)
      w(9)= - 1360.D0/3.D0 + w(9)
      w(9)=z*w(9)
      w(9)=2132.D0/9.D0 + w(9)
      w(9)=z*w(9)
      w(11)=w(21)*w(16)
      w(11)=22 + w(11)
      w(11)=z*w(11)
      w(11)= - 14.D0/3.D0 + w(11)
      w(11)=w(1)*w(11)
      w(12)=NF*w(2)
      w(9)=224.D0/3.D0*w(12) + 64*w(11) - 128*w(31) + 12068.D0/9.D0 + 
     & w(9)
      w(9)=z3*w(9)
      w(11)=w(3)*w(30)
      w(11)=19.D0/3.D0 + w(11)
      w(11)=z*w(11)
      w(3)=w(3)*w(23)
      w(12)=5.D0/3.D0 + w(3)
      w(12)=z*w(12)
      w(12)=5.D0/3.D0 + w(12)
      w(12)=w(1)*w(12)
      w(11)=w(12) + 31.D0/3.D0 + w(11)
      w(12)=2*LL
      w(12)=w(5)*w(12)
      w(11)=w(12) + 1.D0/3.D0*w(11)
      w(11)=w(1)*w(11)
      w(4)=13 + 56*w(4)
      w(4)=w(4)*w(16)
      w(4)=10*w(31) + 13 + w(4)
      w(4)=2.D0/3.D0*w(4) + w(11)
      w(4)=w(1)*w(4)
      w(3)=1 + w(3)
      w(3)=z*w(3)
      w(3)=1 + w(3)
      w(3)=w(1)*w(3)
      w(11)=1 - 11*z
      w(3)=1.D0/3.D0*w(11) + w(3)
      w(3)=z2*w(3)
      w(11)= - 283 - 1453*z
      w(3)=w(3) + 2.D0/81.D0*w(11) + w(4)
      w(3)=NF*w(3)
      w(4)= - 3147610505505388526231345011900288710042786012314622438132
     & 62177946761981119.D0/22.D0 + 335268293009688723381151527766036623
     & 232555567741931860881168533886773585429.D0/25.D0*z
      w(4)=z*w(4)
      w(4)=5392938046976681066594295281108588699905345499331579145460547
     & 2221972477.D0/1034.D0 + 1.D0/252105.D0*w(4)
      w(4)=w(4)*w(14)
      w(4)= - 1269032935208259623550699756558198862937900564984826588307
     & 013392623853.D0/705.D0 + w(4)
      w(11)=1.D0/103823.D0*z
      w(4)=w(4)*w(11)
      w(4)=1217809123723202306964486543559275916338668324060417133149588
     & 201551.D0/76725.D0 + w(4)
      w(4)=w(4)*w(19)
      w(4)= - 6250632198213434326548825350879345010205041668545892986605
     & 57009383.D0/843975.D0 + w(4)
      w(4)=w(4)*w(15)
      w(4)=5632889188506029210664303427541903703623351964338530530769652
     & 7479.D0/2419395.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 1231300905293248577929326184800218759498213928165805781523
     & 4564463201.D0/491137185.D0 + w(4)
      w(13)=1.D0/79507.D0*z
      w(4)=w(4)*w(13)
      w(4)=1276206914888036702855784372790592280969249871428186128926369
     & .D0/4459939.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 7241152432434618122008077451449725070805338419861699540455
     & 39903.D0/2341467975.D0 + w(4)
      w(14)=1.D0/68921.D0*z
      w(4)=w(4)*w(14)
      w(4)=231379795733597961590183593243402229465056271998714352927.D0/
     & 57108975.D0 + w(4)
      w(4)=w(4)*w(16)
      w(4)= - 5579329122725076991221310104151642844255330276272363693.D0
     & /3807265.D0 + w(4)
      w(4)=z*w(4)
      w(4)=1668313839921554117478465202636794593323013473575807109241.D0
     & /1267819245.D0 + w(4)
      w(4)=w(4)*w(26)
      w(4)= - 287097298197554730118312707097693014274654498611282829141.
     & D0/3803457735.D0 + w(4)
      w(15)=1.D0/50653.D0*z
      w(4)=w(4)*w(15)
      w(4)=97654565514946439770632759594163603530594869414071.D0/
     & 73425825.D0 + w(4)
      w(4)=w(4)*w(19)
      w(4)= - 1550029943257145713320082135239169763443395061271.D0/
     & 24475275.D0 + w(4)
      w(4)=z*w(4)
      w(4)=2106836697304499303300857659239627698047804025423.D0/
     & 37528755.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 5504455258127157945255768515689944177732425503.D0/
     & 180138024.D0 + 1.D0/2023.D0*w(4)
      w(4)=z*w(4)
      w(4)=5254306610500847055455972037495229090848892028797.D0/
     & 195449756040.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 636871167285987625434699395311409862639737687.D0/
     & 5816957025.D0 + 1.D0/272.D0*w(4)
      w(4)=w(4)*w(35)
      w(4)=214765628002810213924974409515176076660371.D0/2160929925.D0
     &  + w(4)
      w(4)=z*w(4)
      w(4)= - 255897426736529608953193770360705758653679.D0/87733754955.
     & D0 + 1.D0/38.D0*w(4)
      w(4)=z*w(4)
      w(4)=2138329505643903168889620183281333857.D0/3567769821.D0 + 1.D0
     & /4205.D0*w(4)
      w(4)=z*w(4)
      w(4)= - 120137265633023157237878238345072450877.D0/178388491050.D0
     &  + w(4)
      w(4)=z*w(4)
      w(4)=1864346925293615325969339140176283317.D0/3238716250.D0 + 
     & w(4)
      w(4)=w(4)*w(28)
      w(4)= - 5525757389478152391768277750734436297.D0/110116352500.D0
     &  + w(4)
      w(4)=w(4)*w(18)
      w(4)=34335983973945861736450780579959289.D0/4052281772.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 19644468977348507124266893078510619.D0/2026140886.D0 + 
     & w(4)
      w(4)=w(4)*w(38)
      w(4)=58432406396797640451269784719.D0/3830134.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 30727057263792636429177146909.D0/1740970.D0 + w(4)
      w(17)=1.D0/2197.D0*z
      w(4)=w(4)*w(17)
      w(4)=217603206874301573789323001.D0/33078430.D0 + w(4)
      w(4)=w(4)*w(16)
      w(4)= - 186680127379527804688053127681.D0/72673310710.D0 + w(4)
      w(4)=w(4)*w(39)
      w(4)=997907743273112855035441727.D0/174536732370.D0 + w(4)
      w(4)=w(4)*w(16)
      w(4)= - 95816584321981515303306469.D0/42311935120.D0 + w(4)
      w(4)=w(4)*w(20)
      w(4)=29757955932910301970343.D0/284203920.D0 + w(4)
      w(16)=1.D0/8.D0*z
      w(4)=w(4)*w(16)
      w(4)= - 56134531659613270517.D0/3552549.D0 + w(4)
      w(4)=w(4)*w(43)
      w(4)=466101048904968984427.D0/76971895.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 490647099706986845561.D0/65975910.D0 + w(4)
      w(4)=z*w(4)
      w(4)=281690776947721387.D0/50820.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 116269543968917.D0/605.D0 + 1.D0/36.D0*w(4)
      w(4)=w(4)*w(42)
      w(4)=26124317059159.D0/22680.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 52601455081901.D0/36288.D0 + w(4)
      w(4)=w(4)*w(10)
      w(4)=72276711203.D0/3136.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 124713125209.D0/4410.D0 + w(4)
      w(4)=w(4)*w(7)
      w(4)=97997074.D0/225.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 202217671.D0/450.D0 + w(4)
      w(4)=w(4)*w(37)
      w(4)=703687.D0/36.D0 + w(4)
      w(4)=z*w(4)
      w(4)= - 20530.D0/9.D0 + w(4)
      w(4)=z*w(4)
      w(4)=60884 + w(4)
      w(4)=z*w(4)
      w(4)=493676.D0/3.D0 + w(4)
      w(4)=z*w(4)
      w(4)=416348.D0/3.D0 + w(4)
      w(7)=3251936498016197890488580753947498071979865853113300948741325
     & 281.D0 + 32519581441700503854762800401372481950728247504475648761
     & 32422113.D0*z
      w(7)=z*w(7)
      w(7)=193485660698352322717069844610628641704331298565835334172983.
     & D0 + 1.D0/16807.D0*w(7)
      w(7)=z*w(7)
      w(7)=193484204981727857854516358562481963499884430557108252547983.
     & D0 + w(7)
      w(7)=w(7)*w(11)
      w(7)=1863581810928400887495614920004714515751957589167249921.D0
     &  + w(7)
      w(7)=z*w(7)
      w(7)=1863565880278810884700563441385986724071911791343249921.D0
     &  + w(7)
      w(7)=z*w(7)
      w(7)=1863548863810130582784962555538716618930963228632721921.D0
     &  + w(7)
      w(7)=z*w(7)
      w(7)=20498837266153624877569101189072513963848945492103941131.D0/
     & 11.D0 + w(7)
      w(7)=w(7)*w(13)
      w(7)=257821609833832348126575080735950781504161602488233.D0/11.D0
     &  + w(7)
      w(7)=z*w(7)
      w(7)=12633116994890115497971030645876882487620870265923417.D0/539.
     & D0 + w(7)
      w(7)=w(7)*w(14)
      w(7)=183296302588059550724553475570667194157737036577.D0/539.D0
     &  + w(7)
      w(7)=z*w(7)
      w(7)=183293919397146063185739716882058931672182724577.D0/539.D0
     &  + w(7)
      w(7)=z*w(7)
      w(7)=183291348143104976814019348095145315477110724577.D0/539.D0
     &  + w(7)
      w(7)=z*w(7)
      w(7)=183288568506524148659695092444341102051766724577.D0/539.D0
     &  + w(7)
      w(7)=w(7)*w(15)
      w(7)=3618454135942327005753424024374331114877909.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=3618389596362477493637792479927569538877909.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=3618319365257230640290459410826505793341909.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=3618242753242233214143835521713232129341909.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=3618158963314374493681167370846211841341909.D0/539.D0 + w(7)
      w(7)=w(7)*w(16)
      w(7)=452258383755676134382925365715715643032973.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=15180616603702475646118887931489459603.D0/539.D0 + 1.D0/
     & 29791.D0*w(7)
      w(7)=z*w(7)
      w(7)=15180148658710783705030088844500083603.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=622396597531053462802244872775927.D0/539.D0 + 1.D0/24389.D0*
     & w(7)
      w(7)=z*w(7)
      w(7)=622372998700608819911112463775927.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=2561097446355634460847287110189.D0/539.D0 + 1.D0/243.D0*w(7)
      w(7)=z*w(7)
      w(7)=2560976152652211536408111110189.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=4097343541970887059090818941.D0/539.D0 + 1.D0/625.D0*w(7)
      w(7)=z*w(7)
      w(7)=4097096799131609128541865541.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=336715414881752830592723.D0/539.D0 + 1.D0/12167.D0*w(7)
      w(7)=z*w(7)
      w(7)=336689086344438551975123.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=336658814638864376538323.D0/539.D0 + w(7)
      w(7)=z*w(7)
      w(7)=1683118856778495358491487.D0/2695.D0 + w(7)
      w(7)=z*w(7)
      w(7)=245358578943756786493.D0/2695.D0 + 1.D0/6859.D0*w(7)
      w(7)=z*w(7)
      w(7)=2207911834254200646437.D0/24255.D0 + w(7)
      w(7)=z*w(7)
      w(7)=449325761325072949.D0/24255.D0 + 1.D0/4913.D0*w(7)
      w(7)=w(7)*w(16)
      w(7)=56154295334575853.D0/24255.D0 + w(7)
      w(7)=z*w(7)
      w(7)=56140429821090029.D0/24255.D0 + w(7)
      w(7)=z*w(7)
      w(7)=56123375845866029.D0/24255.D0 + w(7)
      w(7)=w(7)*w(17)
      w(7)=25535765062457.D0/24255.D0 + w(7)
      w(7)=z*w(7)
      w(7)=25523438671457.D0/24255.D0 + w(7)
      w(7)=w(7)*w(42)
      w(7)=19164113947.D0/2205.D0 + w(7)
      w(7)=z*w(7)
      w(7)=19148110939.D0/2205.D0 + w(7)
      w(7)=w(7)*w(44)
      w(7)=78708473.D0/245.D0 + w(7)
      w(7)=w(7)*w(25)
      w(7)=19644962.D0/245.D0 + w(7)
      w(7)=w(7)*w(27)
      w(7)=57134.D0/5.D0 + w(7)
      w(7)=z*w(7)
      w(7)=512206.D0/45.D0 + w(7)
      w(7)=w(7)*w(37)
      w(7)=4070.D0/9.D0 + w(7)
      w(7)=z*w(7)
      w(7)=4016.D0/9.D0 + w(7)
      w(7)=w(7)*w(44)
      w(7)=16 + w(7)
      w(7)=z*w(7)
      w(7)=128.D0/9.D0 + w(7)
      w(7)=LL*z*w(7)
      w(2)=B4*w(2)
      w(2)=1.D0/9.D0*w(9) + 16.D0/9.D0*w(3) + 1.D0/3.D0*w(6) + 32.D0/9.D
     & 0*w(2) + w(8) + 1.D0/243.D0*w(4) + w(7)
      w(2)=as*w(2)
      w(3)= - 5.D0/3.D0*w(5) - w(12)
      w(4)= - w(1)*w(5)
      w(3)=2*w(3) + w(4)
      w(3)=w(3)*w(41)
      w(2)=w(3) + w(2)
      w(2)=1.D0/3.D0*w(2)*as**2
*
      ANSREG1 = w(2)
*
      RETURN
      END            
      REAL*8 FUNCTION ANSREG2(z,nf,as,LL)
*     ----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of aqqQNS3 z \in [1/2,1] 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(26),G,y
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
*      WRITE(6,*) 'ANSREG2::LL=', LL !!
*      LL=10.0D0
      y=1.0D0-z
      w(1)=Log(y)
      w(2)=4*LL
      w(3)=16909.D0/45.D0 - 2164*LL
      w(3)=w(3)*w(2)
      w(3)= - 208480967.D0/30375.D0 + w(3)
      w(4)=16.D0/3.D0*NF
      w(5)= - 11.D0/27.D0 + 47.D0/7.D0*LL
      w(5)=w(5)*w(4)
      w(3)=1.D0/7.D0*w(3) + w(5)
      w(5)=1.D0/9.D0*NF
      w(6)= - 67493625479047241221391709256343.D0/51400597248.D0 + 
     & 8716182513395807182631.D0/5.D0*LL
      w(6)=w(6)*w(5)
      w(7)=127792123416549923443817266116785243095613599.D0/
     & 3099044504245996706400.D0 - 15602948078386910974631.D0*LL
      w(7)=LL*w(7)
      w(7)= - 1556487034515711483687906938163898640602510991230917232807
     & 4471167009.D0/467517090705930589354473771262813060172390400.D0
     &  + w(7)
      w(6)=1.D0/5.D0*w(7) + w(6)
      w(7)= - 14925331762253820388271202721679.D0/1115464350.D0 + 
     & 17556326806961454233518.D0*LL
      w(7)=w(7)*w(5)
      w(8)=43188869880248758264551532425216874653530069.D0/
     & 516507417374332784400.D0 - 31329857936943661817518.D0*LL
      w(8)=LL*w(8)
      w(7)=w(7) - 461828309625864080381880893718768888073162042806097263
     & 19395839091167.D0/682861687718848116290947924567349787216000000.D
     & 0 + w(8)
      w(8)=1.D0/51.D0*y
      w(7)=w(7)*w(8)
      w(6)=1.D0/5.D0*w(6) + w(7)
      w(6)=y*w(6)
      w(7)=367325487875734133044215997244311781732473.D0/
     & 221360321731856907600.D0 - 634273562129706489838.D0*LL
      w(7)=LL*w(7)
      w(6)=w(6) - 878166547229221223270124658233624106178425016435006792
     & 38316941358579.D0/65663661354980192834261731576052534034416640000
     & .D0 + w(7)
      w(7)= - 776852787687916056491097467.D0/25700298624.D0 + 
     & 353181090089253273838.D0/8695.D0*LL
      w(7)=w(7)*w(5)
      w(6)=w(7) + 1.D0/8695.D0*w(6)
      w(7)=1.D0/49.D0*y
      w(6)=w(6)*w(7)
      w(9)= - 540915439090495467200256259.D0/20996976.D0 + 
     & 175272924081937012469.D0/5.D0*LL
      w(9)=w(9)*w(5)
      w(10)=40222662982877262218551298595550914211897.D0/
     & 49191182607079312800.D0 - 315819160102163620469.D0*LL
      w(10)=LL*w(10)
      w(10)= - 253554949901256858277045877121136714824332094947729177309
     & 255992272181.D0/385432012224284777730275893157350030504101840000.
     & D0 + w(10)
      w(9)=1.D0/5.D0*w(10) + w(9)
      w(6)=1.D0/41736.D0*w(9) + w(6)
      w(6)=y*w(6)
      w(9)=161432664501808434407643057895437719897.D0/
     & 4709794079401210800.D0 - 13381851261740686454.D0*LL
      w(9)=LL*w(9)
      w(10)= - 485102046561601472787654110083.D0/90352612350.D0 + 
     & 7401160367262958454.D0*LL
      w(10)=NF*w(10)
      w(9)=1.D0/423.D0*w(10) - 15655718650354990131944215518861949206089
     & 6086915527411267857819.D0/267292458127279157764463214387266811749
     & 760000.D0 + 1.D0/47.D0*w(9)
      w(6)=1.D0/185.D0*w(9) + w(6)
      w(6)=y*w(6)
      w(9)= - 3633571100543601342569172626399.D0/1382787806400.D0 + 
     & 3671326804256316427.D0*LL
      w(9)=w(9)*w(5)
      w(10)=158952130801768359380961786769111783097.D0/
     & 9419588158802421600.D0 - 6661672251495180427.D0*LL
      w(10)=LL*w(10)
      w(9)=w(9) - 227881424462071656493949742767364025900252408591483596
     & 40213037.D0/1684262496076113155415647223612267244800000.D0 + 
     & w(10)
      w(6)=1.D0/4255.D0*w(9) + w(6)
      w(6)=y*w(6)
      w(9)= - 1385611798595424713269058929.D0/269549280.D0 + 
     & 7282846699567855574.D0*LL
      w(9)=w(9)*w(5)
      w(10)=469302038176879399482011202043736022571.D0/
     & 14129382238203632400.D0 - 13263537594045583574.D0*LL
      w(10)=LL*w(10)
      w(9)=w(9) - 295255335736805674223968792775168707488775603082935497
     & 9773536109.D0/111049040574618394047071673610168820340480000.D0
     &  + w(10)
      w(6)=1.D0/8325.D0*w(9) + w(6)
      w(6)=y*w(6)
      w(9)= - 128945357280002551661490582089.D0/564638354280.D0 + 
     & 328258206486023417.D0*LL
      w(9)=w(9)*w(5)
      w(10)=41966455851575670634261779637585175761.D0/
     & 28258764476407264800.D0 - 600107792598647417.D0*LL
      w(10)=LL*w(10)
      w(9)=w(9) - 148279530059448909440606666441716951110926134972472162
     & 2680595169.D0/1249301706464456933029556328114399228830400000.D0
     &  + w(10)
      w(6)=1.D0/370.D0*w(9) + w(6)
      w(6)=y*w(6)
      w(9)= - 5856782799632006342708582453.D0/564638354280.D0 + 
     & 15135500881718038.D0*LL
      w(9)=w(9)*w(5)
      w(10)=22313741517747181321596094057949689.D0/328590284609386800.D0
     &  - 27779667677654038.D0*LL
      w(10)=LL*w(10)
      w(9)=w(9) - 768046812188300528666524408817154822469530166755624565
     & 389119.D0/14173219203729736420600196309245577174400000.D0 + 
     & w(10)
      w(6)=1.D0/185.D0*w(9) + 1.D0/11.D0*w(6)
      w(9)=1.D0/43.D0*y
      w(6)=w(6)*w(9)
      w(10)= - 3628135372543762865661161213.D0/102428726400.D0 + 
     & 52500096831165533.D0*LL
      w(10)=w(10)*w(5)
      w(11)=153471428205406992292227736196514223.D0/657180569218773600.D
     & 0 - 96754680616941533.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 7176654668621564463652748793322543594271636003070344
     & 8763353.D0/385665828672917997839461124061104140800000.D0 + w(11)
      w(6)=1.D0/27195.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 10091991963845938165912906879.D0/5992080494400.D0 + 
     & 2537286623883626.D0*LL
      w(10)=w(10)*w(5)
      w(11)=12807172573839974392030917535369.D0/1144913883656400.D0 - 
     & 4696046808555626.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 1359789792722997466743151218318448714258870466720838
     & 4363.D0/1527644277182667306193654863810470400000.D0 + w(11)
      w(6)=1.D0/1295.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)=1795411187722241465235575549047.D0/2289827767312800.D0 - 
     & 2335880378239033.D0/7.D0*LL
      w(10)=LL*w(10)
      w(11)= - 14451622959883679884549909.D0/17623766160.D0 + 
     & 1256500285903033.D0*LL
      w(11)=NF*w(11)
      w(10)=1.D0/63.D0*w(11) - 83291626869374132653338438540993488810881
     & 705078059121.D0/133668874253483389291944800583416160000.D0 + 
     & w(10)
      w(6)=1.D0/3700.D0*w(10) + 1.D0/41.D0*w(6)
      w(6)=y*w(6)
      w(10)= - 1021229271353694818549489.D0/640179540.D0 + 
     & 2488091800444466.D0*LL
      w(10)=w(10)*w(5)
      w(11)=12324493957416544323372383179729.D0/1144913883656400.D0 - 
     & 4646851985116466.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 4524833571764109065521462038433927373919324392846033
     & 537.D0/529553218372492660874822174190710400000.D0 + w(11)
      w(6)=1.D0/50505.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 3759801109930430781713921.D0/4851887040.D0 + 
     & 1231263767549833.D0*LL
      w(10)=w(10)*w(5)
      w(11)=635624153884224476594517589291.D0/120517250911200.D0 - 
     & 2310643859885833.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 4785029869274667737277074480719401013454317642339757
     & 57.D0/114418700922588441961872935143699200000.D0 + w(11)
      w(6)=1.D0/24605.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 70428305394080152125461.D0/1732816800.D0 + 
     & 65845198477418.D0*LL
      w(10)=w(10)*w(5)
      w(11)=8637547593039269316670252441.D0/30943618477200.D0 - 
     & 124190068333418.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 412642752993023943958870193075541018472803719693.D0/
     & 1866025633512569663907271699200000.D0 + w(11)
      w(6)=1.D0/1295.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 552902754619936428977.D0/28078050.D0 + 32557943802109.D0
     & *LL
      w(10)=w(10)*w(5)
      w(11)=8450068494235272039766368041.D0/61887236954400.D0 - 
     & 61730378730109.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 87623268441080531140762858687429018699271975747.D0/
     & 810684469670460820653048038208000.D0 + w(11)
      w(6)=1.D0/630.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 4035128413057280960461.D0/105894360.D0 + 64365739277498.D
     & 0*LL
      w(10)=w(10)*w(5)
      w(11)=2753028761751790627852350467.D0/10314539492400.D0 - 
     & 122710609133498.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 200090021397920701860515016715287690265136995823.D0/
     & 947553276238200959204861343360000.D0 + w(11)
      w(6)=1.D0/1225.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 917151256026532177199.D0/49832640.D0 + 31796763882349.D0
     & *LL
      w(10)=w(10)*w(5)
      w(11)=18817051258873156640872262069.D0/144403552893600.D0 - 
     & 60969198810349.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 4162912878910592057990300715382692696393566485661.D0
     & /40393966654860331548110352076800000.D0 + w(11)
      w(6)=1.D0/595.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 1276195637136675109.D0/22880.D0 + 690777074933278.D0/7.D0
     & *LL
      w(10)=w(10)*w(5)
      w(11)=18354001509868391567206274869.D0/6563797858800.D0 - 
     & 1332570643349278.D0*LL
      w(11)=LL*w(11)
      w(11)= - 5894340853639570206409589362378128179217346217491.D0/
     & 2665617383810592878268739737600000.D0 + w(11)
      w(10)=1.D0/7.D0*w(11) + w(10)
      w(6)=1.D0/1815.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 7254448834607777009.D0/77220.D0 + 170437963219357.D0*LL
      w(10)=w(10)*w(5)
      w(11)=5463837926714342168438197469.D0/8022419605200.D0 - 
     & 330886355323357.D0*LL
      w(11)=LL*w(11)
      w(10)=w(10) - 23137505302738412128117579685565559729761457859.D0/
     & 42951842610229279776791216475000.D0 + w(11)
      w(6)=1.D0/3080.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)=4.D0/9.D0*NF
      w(11)= - 56187908486266552.D0/19305.D0 + 5422866892147.D0*LL
      w(11)=w(11)*w(10)
      w(12)=5532295373421710399878279.D0/64696932300.D0 - 
     & 42394485904588.D0*LL
      w(12)=LL*w(12)
      w(11)=w(11) - 45269427371054318427329587369079222339383753.D0/
     & 668982409826672008876208400000.D0 + w(12)
      w(6)=1.D0/385.D0*w(11) + w(6)
      w(11)=1.D0/31.D0*y
      w(6)=w(6)*w(11)
      w(12)=48381980246088280050926831.D0/1164544781400.D0 - 
     & 21041970314774.D0*LL
      w(12)=LL*w(12)
      w(12)= - 91726676452856753088347062544981653820969.D0/
     & 2787426707611133370317535000.D0 + w(12)
      w(13)= - 2344977186032647.D0/21021.D0 + 5345230573387.D0/25.D0*LL
      w(13)=NF*w(13)
      w(12)=1.D0/25.D0*w(12) + 2.D0/9.D0*w(13)
      w(6)=1.D0/231.D0*w(12) + w(6)
      w(6)=y*w(6)
      w(12)= - 31299198369588007.D0/1702701.D0 + 181548866903.D0/5.D0*
     & LL
      w(10)=w(12)*w(10)
      w(12)=55817752799970056262391.D0/20078358300.D0 - 1440092651612.D0
     & *LL
      w(12)=LL*w(12)
      w(12)= - 120963365749066178866313701154395867559.D0/
     & 54859355432914183351200000.D0 + w(12)
      w(10)=1.D0/5.D0*w(12) + w(10)
      w(6)=1.D0/77.D0*w(10) + w(6)
      w(6)=y*w(6)
      w(10)= - 8522452202221537.D0/486486.D0 + 178680530003.D0/5.D0*LL
      w(5)=w(10)*w(5)
      w(10)=54068236923548084518591.D0/80313433200.D0 - 357154826003.D0
     & *LL
      w(10)=LL*w(10)
      w(10)= - 13626216002989625143108545626310121729.D0/
     & 25470415022424442270200000.D0 + w(10)
      w(5)=1.D0/5.D0*w(10) + w(5)
      w(5)=1.D0/539.D0*w(5) + 1.D0/29.D0*w(6)
      w(5)=y*w(5)
      w(6)=17426204357357618897797.D0/20078358300.D0 - 472240339204.D0*
     & LL
      w(6)=LL*w(6)
      w(6)= - 10006412895582619353053080832098501.D0/
     & 14471709252114114000000.D0 + w(6)
      w(10)= - 2498830686386303.D0/1216215.D0 + 4338418726.D0*LL
      w(10)=NF*w(10)
      w(6)=1.D0/3.D0*w(6) + 2*w(10)
      w(5)=1.D0/1155.D0*w(6) + w(5)
      w(6)=1.D0/3.D0*y
      w(5)=w(5)*w(6)
      w(10)=50446779877570819104991.D0/120470149800.D0 - 234060850802.D0
     & *LL
      w(10)=LL*w(10)
      w(12)= - 52402508634452567.D0/8981280.D0 + 12786442978.D0*LL
      w(12)=NF*w(12)
      w(10)=w(12) - 22744375821785033271219377138312527.D0/
     & 67954113009927144000000.D0 + w(10)
      w(5)=1.D0/5005.D0*w(10) + w(5)
      w(10)=1.D0/9.D0*y
      w(5)=w(5)*w(10)
      w(12)=388563735481220774987.D0/60235074900.D0 
     & - 3710706548.0D0*LL
      w(12)=LL*w(12)
      w(12)= - 19356337834672775086257269442047077.D0/
     & 3738563481354151754304000.D0 + w(12)
      w(13)= - 1976182397174329.D0/202078800.D0 + 22308404*LL
      w(13)=NF*w(13)
      w(12)=1.D0/9.D0*w(12) + w(13)
      w(5)=1.D0/77.D0*w(12) + w(5)
      w(5)=y*w(5)
      w(12)=1.D0/3.D0*NF
      w(13)= - 25758441486989.D0/83160.D0 + 738035777*LL
      w(13)=w(13)*w(12)
      w(14)=373177742000759371187.D0/48188059920.D0 - 4593764611.D0*LL
      w(14)=LL*w(14)
      w(13)=w(13) - 1982574953489661086734352189506613.D0/
     & 317316344868022139640000.D0 + w(14)
      w(5)=1.D0/4158.D0*w(13) + 1.D0/5.D0*w(5)
      w(13)=1.D0/13.D0*y
      w(5)=w(5)*w(13)
      w(14)= - 2533414824527.D0/654885.D0 + 125655004.D0/13.D0*LL
      w(14)=w(14)*w(12)
      w(15)=675604855461522443.D0/523783260.D0 - 790818452*LL
      w(15)=LL*w(15)
      w(15)= - 4425119375760251576381053250081.D0/
     & 4224973113225905040000.D0 + w(15)
      w(14)=1.D0/13.D0*w(15) + w(14)
      w(5)=1.D0/693.D0*w(14) + w(5)
      w(14)=1.D0/23.D0*y
      w(5)=w(5)*w(14)
      w(15)= - 555655940079883.D0/23814000.D0 + 61416638*LL
      w(12)=w(15)*w(12)
      w(15)=644976040219704923.D0/1047566520.D0 - 391176634*LL
      w(15)=LL*w(15)
      w(12)=w(12) - 35192080083278878115790305687.D0/
     & 69834266334312480000.D0 + w(15)
      w(5)=1.D0/99099.D0*w(12) + w(5)
      w(5)=y*w(5)
      w(12)=613506144953719643.D0/523783260.D0 - 773484980*LL
      w(12)=LL*w(12)
      w(15)= - 8126873820479.D0/567000.D0 + 39959060*LL
      w(15)=NF*w(15)
      w(12)=w(12) + w(15)
      w(12)= - 78468135995698084307592157.D0/567360997623763200.D0 + 1.D
     & 0/7.D0*w(12)
      w(5)=1.D0/27027.D0*w(12) + w(5)
      w(5)=y*w(5)
      w(12)=11860076155861331.D0/59860944.D0 - 955216597.D0/7.D0*LL
      w(12)=LL*w(12)
      w(12)= - 145143733495501790809276962659.D0/875721699832278499200.D
     & 0 + w(12)
      w(15)= - 14404708156.D0/76545.D0 + 559259*LL
      w(15)=NF*w(15)
      w(12)=1.D0/3.D0*w(12) + 29.D0/7.D0*w(15)
      w(5)=1.D0/10725.D0*w(12) + w(5)
      w(5)=y*w(5)
      w(12)=216791836487021.D0/787644.D0 - 198518812*LL
      w(12)=LL*w(12)
      w(15)= - 159150064811.D0/51030.D0 + 9956668*LL
      w(15)=NF*w(15)
      w(12)=w(15) - 362286686475108244772539039.D0/1550337293769264000.D
     & 0 + w(12)
      w(5)=1.D0/45045.D0*w(12) + w(5)
      w(12)=1.D0/19.D0*y
      w(5)=w(5)*w(12)
      w(15)=1.D0/27.D0*NF
      w(16)= - 10540589023.D0/840.D0 + 43443646*LL
      w(16)=w(16)*w(15)
      w(17)=203210665287997.D0/4725864.D0 - 32632682*LL
      w(17)=LL*w(17)
      w(16)=w(16) - 34157245625301312369488233.D0/918718396307712000.D0
     &  + w(17)
      w(5)=1.D0/135135.D0*w(16) + w(5)
      w(5)=y*w(5)
      w(16)=1426332151.D0/5148.D0 - 3782612.D0/17.D0*LL
      w(16)=LL*w(16)
      w(17)= - 24585163657.D0/18900.D0 + 4941436*LL
      w(17)=NF*w(17)
      w(16)=1.D0/459.D0*w(17) - 1502261819252573567479.D0/
     & 6130844463744000.D0 + w(16)
      w(5)=1.D0/15015.D0*w(16) + 1.D0/17.D0*w(5)
      w(5)=y*w(5)
      w(16)= - 120858193.D0/9450.D0 + 595157.D0/11.D0*LL
      w(16)=w(16)*w(15)
      w(17)=5598654097.D0/10296.D0 - 465319*LL
      w(17)=LL*w(17)
      w(17)= - 11380973800387978691.D0/23009449963500.D0 + w(17)
      w(16)=1.D0/11.D0*w(17) + w(16)
      w(5)=1.D0/2730.D0*w(16) + w(5)
      w(5)=y*w(5)
      w(16)=230397184207.D0/173745.D0 - 1219496*LL
      w(16)=LL*w(16)
      w(17)= - 1783688.D0/15.D0 + 571133*LL
      w(17)=NF*w(17)
      w(16)=8.D0/81.D0*w(17) - 179409825158187373493.D0/143578967772240.
     & D0 + w(16)
      w(5)=1.D0/75075.D0*w(16) + w(5)
      w(5)=w(5)*w(6)
      w(16)=1456802335889.D0/7297290.D0 - 199436*LL
      w(16)=LL*w(16)
      w(17)= - 210571399.D0/540.D0 + 2181572*LL
      w(17)=NF*w(17)
      w(16)=1.D0/243.D0*w(17) - 51881335762931863541.D0/263716471418400.
     & D0 + w(16)
      w(5)=1.D0/35035.D0*w(16) + w(5)
      w(5)=y*w(5)
      w(16)=7664643401.D0/93555.D0 - 90152*LL
      w(16)=LL*w(16)
      w(17)= - 12741569.D0/270.D0 + 318568*LL
      w(17)=NF*w(17)
      w(16)=1.D0/81.D0*w(17) - 263468253082233053.D0/3080893384800.D0
     &  + w(16)
      w(5)=1.D0/15015.D0*w(16) + w(5)
      w(5)=y*w(5)
      w(16)=2.D0/27.D0*NF
      w(17)= - 65683.D0/15.D0 + 37511*LL
      w(17)=w(17)*w(16)
      w(18)=6674811481.D0/124740.D0 - 66074*LL
      w(18)=LL*w(18)
      w(17)=w(17) - 7362680170055611.D0/122257674000.D0 + w(18)
      w(5)=1.D0/10395.D0*w(17) + w(5)
      w(5)=y*w(5)
      w(17)= - 20953.D0/135.D0 + 12724.D0/7.D0*LL
      w(16)=w(17)*w(16)
      w(17)=15562307.D0/945.D0 - 23416*LL
      w(17)=LL*w(17)
      w(17)= - 39548581810757.D0/1928934000.D0 + w(17)
      w(16)=1.D0/7.D0*w(17) + w(16)
      w(5)=1.D0/495.D0*w(16) + w(5)
      w(5)=y*w(5)
      w(16)= - 3539.D0/27.D0 + 11716.D0/5.D0*LL
      w(15)=w(16)*w(15)
      w(16)=12635579.D0/1890.D0 - 11372*LL
      w(16)=LL*w(16)
      w(15)=w(15) - 236436849433.D0/123451776.D0 + 1.D0/5.D0*w(16)
      w(5)=1.D0/315.D0*w(15) + w(5)
      w(5)=y*w(5)
      w(15)=2.D0/3.D0*NF
      w(16)= - 23 + 3532.D0/5.D0*LL
      w(16)=w(16)*w(15)
      w(17)=3203113.D0/105.D0 - 65992*LL
      w(17)=LL*w(17)
      w(17)= - 1802688043.D0/33600.D0 + w(17)
      w(16)=1.D0/5.D0*w(17) + w(16)
      w(5)=1.D0/1701.D0*w(16) + w(5)
      w(5)=y*w(5)
      w(16)= - 53.D0/27.D0 + 389.D0/5.D0*LL
      w(16)=w(16)*w(15)
      w(17)=1623331.D0/630.D0 - 7934*LL
      w(17)=LL*w(17)
      w(17)= - 1872782896.D0/297675.D0 + w(17)
      w(16)=1.D0/5.D0*w(17) + w(16)
      w(5)=1.D0/189.D0*w(16) + w(5)
      w(5)=y*w(5)
      w(3)=1.D0/135.D0*w(3) + w(5)
      w(3)=y*w(3)
      w(5)=2*LL
      w(16)=689.D0/45.D0 - 2044*LL
      w(16)=w(16)*w(5)
      w(16)= - 30984001.D0/9000.D0 + w(16)
      w(17)= - 5.D0/3.D0 + 37.D0/5.D0*LL
      w(17)=NF*w(17)
      w(16)=1.D0/5.D0*w(16) + 8.D0/3.D0*w(17)
      w(3)=1.D0/81.D0*w(16) + w(3)
      w(3)=y*w(3)
      w(16)= - 197.D0/15.D0 - 76*LL
      w(16)=w(16)*w(2)
      w(17)= - 41.D0/45.D0 + LL
      w(17)=w(17)*w(4)
      w(16)=w(17) - 102215.D0/324.D0 + w(16)
      w(3)=1.D0/27.D0*w(16) + w(3)
      w(3)=y*w(3)
      w(16)= - 46.D0/3.D0 - 43*LL
      w(16)=w(16)*w(2)
      w(17)=4.D0/3.D0*NF
      w(18)= - 19.D0/3.D0 + LL
      w(18)=w(18)*w(17)
      w(16)=w(18) - 26515.D0/81.D0 + w(16)
      w(3)=2.D0/27.D0*w(16) + w(3)
      w(3)=y*w(3)
      w(16)=8*LL
      w(18)= - 2 - 37.D0/9.D0*LL
      w(18)=w(18)*w(16)
      w(19)= - 74.D0/3.D0 - LL
      w(19)=NF*w(19)
      w(18)=8.D0/27.D0*w(19) - 33419.D0/81.D0 + w(18)
      w(3)=4.D0/9.D0*w(18) + w(3)
      w(3)=y*w(3)
      w(18)= - 5 - 28*LL
      w(18)=w(18)*w(2)
      w(19)=938.D0/27.D0 - LL
      w(4)=w(19)*w(4)
      w(4)=w(4) - 5111.D0/81.D0 + w(18)
      w(3)=4.D0/27.D0*w(4) + w(3)
      w(3)=y*w(3)
      w(4)=256*LL
      w(18)=214538028770087821079.D0/79057257761377467.D0 - w(4)
      w(18)=w(18)*w(5)
      w(18)= - 205668991374578074580834424540209637567351469777.D0/
     & 28228425847048363321607225875948510001568000.D0 + w(18)
      w(19)=10637590106028325140599.D0/3873805630307495883.D0 - w(4)
      w(19)=w(19)*w(5)
      w(19)= - 68883469663160217066460687951858543192638413919263.D0/
     & 9277908972374488468911386211222360129778500000.D0 + w(19)
      w(19)=w(19)*w(8)
      w(20)=10575609215943405206471.D0/3873805630307495883.D0 - w(4)
      w(20)=LL*w(20)
      w(20)= - 293097763031312768289199519376799749850443795929.D0/
     & 79692354962711395690859906824815219641044800.D0 + w(20)
      w(19)=1.D0/25.D0*w(20) + w(19)
      w(19)=y*w(19)
      w(18)=1.D0/49.D0*w(18) + w(19)
      w(18)=y*w(18)
      w(19)=32*LL
      w(20)=213220407807398196629.D0/632458062091019736.D0 - w(19)
      w(20)=LL*w(20)
      w(20)= - 208511922403713075957563019225645090868350118417.D0/
     & 462404592000341649178071854101656551915220000.D0 + w(20)
      w(18)=1.D0/3.D0*w(20) + w(18)
      w(18)=y*w(18)
      w(20)=4507973454385972507.D0/1682069314071861.D0 - w(4)
      w(20)=LL*w(20)
      w(18)=w(18) - 1309699764757090823636498416335797030657343901.D0/
     & 8618060405889456971572474670224376370900000.D0 + 2.D0/47.D0*
     & w(20)
      w(18)=y*w(18)
      w(20)=4478720075010809707.D0/1682069314071861.D0 - w(4)
      w(20)=LL*w(20)
      w(20)= - 115165955927957921869428717185508593903245769.D0/
     & 32582459001472427113695556409165884200000.D0 + w(20)
      w(18)=1.D0/23.D0*w(20) + w(18)
      w(18)=y*w(18)
      w(20)=4448816620538421067.D0/1682069314071861.D0 - w(4)
      w(20)=w(20)*w(5)
      w(20)= - 263613076004593899783840237664397913089791993.D0/
     & 37688949651995591129175327823000069560000.D0 + w(20)
      w(18)=1.D0/45.D0*w(20) + w(18)
      w(18)=y*w(18)
      w(20)=128*LL
      w(21)=401657594736431897.D0/305830784376702.D0 - w(20)
      w(21)=LL*w(21)
      w(21)= - 7843251488443750207221055625026685280939089.D0/4534766669
     & 357758290943555486724607300000.D0 + w(21)
      w(18)=1.D0/11.D0*w(21) + w(18)
      w(18)=y*w(18)
      w(21)=9274712958310379.D0/3556171911357.D0 - w(4)
      w(21)=LL*w(21)
      w(18)=w(18) - 1321806365648256675368393774729741439009457.D0/
     & 8308615153701595123908845239005794700000.D0 + 2.D0/43.D0*w(21)
      w(18)=y*w(18)
      w(21)=64448834453325053.D0/24893203379499.D0 - w(4)
      w(21)=LL*w(21)
      w(21)= - 764373432244623119448091450544586854176529.D0/22608476608
     & 7118234664186265006960400000.D0 + w(21)
      w(18)=1.D0/21.D0*w(21) + w(18)
      w(18)=y*w(18)
      w(21)=1560075936872533.D0/607151301939.D0 - w(4)
      w(21)=LL*w(21)
      w(18)=w(18) - 11225805844267985935913743346354027009.D0/6888727801
     & 1483915322585446600400000.D0 + 2.D0/41.D0*w(21)
      w(18)=y*w(18)
      w(21)=64*LL
      w(22)=1547932910833753.D0/2428605207756.D0 - w(21)
      w(22)=LL*w(22)
      w(22)= - 53984000190507718656527887017192981841.D0/654429141109097
     & 19556456174270380000.D0 + w(22)
      w(18)=1.D0/5.D0*w(22) + w(18)
      w(18)=y*w(18)
      w(22)=1535478525152953.D0/607151301939.D0 - w(4)
      w(22)=w(22)*w(5)
      w(22)= - 11968154594855919954016768072353715991.D0/183688956311914
     & 7036556575973300000.D0 + w(22)
      w(18)=1.D0/39.D0*w(22) + w(18)
      w(18)=y*w(18)
      w(22)=1522696392480553.D0/607151301939.D0 - w(4)
      w(22)=LL*w(22)
      w(22)= - 5529282918630845962471309804067711501.D0/1719857818100475
     & 618715021271400000.D0 + w(22)
      w(18)=1.D0/19.D0*w(22) + w(18)
      w(18)=y*w(18)
      w(22)=40799156669269.D0/16409494647.D0 - w(4)
      w(22)=LL*w(22)
      w(18)=w(18) - 2322335619382188354148147082942183.D0/13547534586738
     & 984302357086200000.D0 + 2.D0/37.D0*w(22)
      w(18)=y*w(18)
      w(22)=40434501232669.D0/32818989294.D0 - w(20)
      w(22)=LL*w(22)
      w(22)= - 20002663112109650328359877723597119.D0/127948937763645962
     & 85559470300000.D0 + w(22)
      w(18)=1.D0/9.D0*w(22) + w(18)
      w(18)=y*w(18)
      w(22)=40059427069309.D0/16409494647.D0 - w(4)
      w(22)=w(22)*w(5)
      w(22)= - 20273865701860204646786695374873493.D0/
     & 3290115542493753330572435220000.D0 + w(22)
      w(18)=1.D0/35.D0*w(22) + w(18)
      w(18)=y*w(18)
      w(22)=39673321312909.D0/16409494647.D0 - w(4)
      w(22)=LL*w(22)
      w(22)= - 9672398470923714729490650860215079.D0/
     & 3187655196879761012319314400000.D0 + w(22)
      w(18)=1.D0/17.D0*w(22) + w(18)
      w(18)=y*w(18)
      w(4)=432030669202799.D0/180504441117.D0 - w(4)
      w(4)=w(4)*w(5)
      w(4)= - 6910516149472855506447161379944297.D0/
     & 1156951989501125381193029400000.D0 + w(4)
      w(4)=1.D0/33.D0*w(4) + w(18)
      w(4)=y*w(4)
      w(18)=16*LL
      w(22)=213759029087437.D0/1444035528936.D0 - w(18)
      w(22)=LL*w(22)
      w(4)=w(4) - 21904905936361037617972261472318.D0/
     & 119310673917303554935531156875.D0 + w(22)
      w(4)=y*w(4)
      w(22)=6820320629827.D0/5822723907.D0 - w(20)
      w(22)=LL*w(22)
      w(4)=w(4) - 173077848501555499361489204767.D0/
     & 929142235870377790105845000.D0 + 4.D0/31.D0*w(22)
      w(4)=y*w(4)
      w(22)=6742684311067.D0/5822723907.D0 - w(20)
      w(22)=w(22)*w(5)
      w(22)= - 878268845220646765365668977423.D0/
     & 309714078623459263368615000.D0 + w(22)
      w(4)=1.D0/15.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=229736926823.D0/200783583.D0 - w(20)
      w(22)=LL*w(22)
      w(4)=w(4) - 73113976945966942896442583.D0/
     & 380967746061904051050000.D0 + 4.D0/29.D0*w(22)
      w(4)=y*w(4)
      w(22)=226868589923.D0/200783583.D0 - w(20)
      w(22)=LL*w(22)
      w(22)= - 241292213642979555979311647.D0/176877882100169737987500.D
     & 0 + w(22)
      w(4)=1.D0/7.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=921374561.D0/7436429.D0 - 128.D0/9.D0*LL
      w(22)=w(22)*w(2)
      w(22)= - 90774253136569800856471867.D0/152857428975455329125000.D0
     &  + w(22)
      w(4)=1.D0/3.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=8177964449.D0/7436429.D0 - w(20)
      w(22)=w(22)*w(5)
      w(22)= - 1182428460642212537765597.D0/452240914128566062500.D0 + 
     & w(22)
      w(4)=1.D0/13.D0*w(22) + w(4)
      w(4)=w(4)*w(6)
      w(22)=128.D0/15.D0*LL
      w(23)=537265439.D0/7436429.D0 - w(22)
      w(23)=w(23)*w(2)
      w(23)= - 1768935041877020063248819.D0/5192449279658544103200.D0
     &  + w(23)
      w(4)=1.D0/5.D0*w(23) + w(4)
      w(4)=y*w(4)
      w(23)=23805123305.D0/44618574.D0 - w(21)
      w(23)=LL*w(23)
      w(23)= - 1648435320054286592079079.D0/2644302873900184497000.D0
     &  + w(23)
      w(4)=1.D0/9.D0*w(23) + w(4)
      w(4)=y*w(4)
      w(23)=1018136335.D0/969969.D0 - w(20)
      w(23)=LL*w(23)
      w(23)= - 59043776657004209701.D0/279429438705417000.D0 + 4.D0/23.D
     & 0*w(23)
      w(4)=1.D0/3.D0*w(23) + w(4)
      w(4)=y*w(4)
      w(23)=1000500535.D0/969969.D0 - w(20)
      w(23)=w(23)*w(5)
      w(23)= - 2519747087061651241.D0/1065846555774000.D0 + w(23)
      w(4)=1.D0/33.D0*w(23) + w(4)
      w(4)=y*w(4)
      w(23)=327341645.D0/323323.D0 - w(20)
      w(23)=w(23)*w(2)
      w(23)= - 90451358385156453457.D0/19700034639714000.D0 + w(23)
      w(4)=1.D0/63.D0*w(23) + w(4)
      w(4)=y*w(4)
      w(22)=21391679.D0/323323.D0 - w(22)
      w(22)=LL*w(22)
      w(4)=w(4) - 56844223092223389481.D0/766563112598283000.D0 + w(22)
      w(4)=y*w(4)
      w(22)=16529915.D0/17017.D0 - w(20)
      w(22)=LL*w(22)
      w(22)= - 2437231154717898757.D0/10766231206731000.D0 + 4.D0/19.D0
     & *w(22)
      w(4)=1.D0/3.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=145365835.D0/153153.D0 - w(20)
      w(22)=w(22)*w(5)
      w(22)= - 1102269453970229029.D0/531665738604000.D0 + w(22)
      w(4)=1.D0/27.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(20)=8338955.D0/9009.D0 - w(20)
      w(20)=LL*w(20)
      w(20)= - 342408701086817.D0/1460917458000.D0 + 4.D0/17.D0*w(20)
      w(4)=1.D0/3.D0*w(20) + w(4)
      w(4)=y*w(4)
      w(20)=4056865.D0/18018.D0 - w(19)
      w(20)=LL*w(20)
      w(20)= - 30472983788123.D0/127830277575.D0 + w(20)
      w(4)=1.D0/3.D0*w(20) + w(4)
      w(4)=y*w(4)
      w(20)=64.D0/5.D0*LL
      w(22)=787349.D0/9009.D0 - w(20)
      w(22)=w(22)*w(16)
      w(22)= - 32212792942643.D0/44314496226.D0 + w(22)
      w(4)=1.D0/9.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=3808045.D0/9009.D0 - w(21)
      w(22)=w(22)*w(2)
      w(22)= - 4326208642747.D0/2512159650.D0 + w(22)
      w(4)=1.D0/21.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=282265.D0/693.D0 - w(21)
      w(22)=LL*w(22)
      w(22)= - 1185453064309.D0/4754465100.D0 + 8.D0/13.D0*w(22)
      w(4)=1.D0/3.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=270715.D0/693.D0 - w(21)
      w(22)=w(22)*w(5)
      w(22)= - 2994966351931.D0/3962054250.D0 + w(22)
      w(4)=1.D0/9.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=23465.D0/63.D0 - w(21)
      w(22)=LL*w(22)
      w(22)= - 2263398049.D0/8930250.D0 + 8.D0/11.D0*w(22)
      w(4)=1.D0/3.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(20)=4441.D0/63.D0 - w(20)
      w(20)=w(20)*w(2)
      w(20)= - 451868423.D0/1786050.D0 + w(20)
      w(4)=1.D0/3.D0*w(20) + w(4)
      w(4)=y*w(4)
      w(20)=6935.D0/21.D0 - w(21)
      w(20)=w(20)*w(16)
      w(20)= - 76945487.D0/34300.D0 + w(20)
      w(4)=1.D0/27.D0*w(20) + w(4)
      w(4)=y*w(4)
      w(20)=3205.D0/21.D0 - w(19)
      w(20)=LL*w(20)
      w(20)= - 16652719.D0/138915.D0 + w(20)
      w(4)=2.D0/3.D0*w(20) + w(4)
      w(4)=y*w(4)
      w(20)=415.D0/3.D0 - w(19)
      w(20)=LL*w(20)
      w(20)= - 555779.D0/10125.D0 + 4.D0/7.D0*w(20)
      w(4)=4.D0/3.D0*w(20) + w(4)
      w(4)=y*w(4)
      w(20)=8.D0/3.D0*LL
      w(22)=365.D0/3.D0 - w(19)
      w(22)=w(22)*w(20)
      w(22)= - 22183.D0/125.D0 + w(22)
      w(4)=1.D0/3.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=61.D0/3.D0 - 32.D0/5.D0*LL
      w(22)=w(22)*w(18)
      w(22)= - 6593.D0/81.D0 + w(22)
      w(4)=1.D0/3.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=115.D0/3.D0 - w(18)
      w(22)=w(22)*w(5)
      w(22)=4141.D0/81.D0 + w(22)
      w(4)=4.D0/3.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)=65.D0/3.D0 - w(18)
      w(22)=w(22)*w(16)
      w(22)=4567.D0/3.D0 + w(22)
      w(4)=4.D0/9.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)= - 5.D0/3.D0 - w(16)
      w(22)=LL*w(22)
      w(22)= - 3637.D0/27.D0 + w(22)
      w(4)=32.D0/3.D0*w(22) + w(4)
      w(4)=y*w(4)
      w(22)= - 1670060319654470380055477939.D0/464.D0 - 
     & 88562257476676995473072752.D0/25.D0*y
      w(22)=y*w(22)
      w(22)= - 90448162065138546568396313.D0/10904.D0 + 1.D0/441.D0*
     & w(22)
      w(22)=y*w(22)
      w(22)= - 9517360018230344573997907153.D0/1128564.D0 + w(22)
      w(22)=y*w(22)
      w(22)= - 10954205332882946180494529.D0/60030.D0 + 1.D0/47.D0*
     & w(22)
      w(22)=y*w(22)
      w(22)= - 726798496214482651847222.D0/3915.D0 + w(22)
      w(22)=y*w(22)
      w(22)= - 31805934126373566280113683.D0/168345.D0 + w(22)
      w(22)=y*w(22)
      w(22)= - 32381467964455104894052411.D0/168345.D0 + w(22)
      w(9)=w(22)*w(9)
      w(9)= - 731310385763903160172348.D0/160515.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 2980557847629019628498081.D0/642060.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 12212082186665827342.D0/3915.D0 + 1.D0/1517.D0*w(9)
      w(9)=y*w(9)
      w(9)= - 17510203547109574980031.D0/5504490.D0 + w(9)
      w(9)=w(9)*w(13)
      w(9)= - 25423756777343638319509.D0/101833065.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 1365904760258785862737.D0/5359635.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 339668339807063017.D0/1305.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 54467579001758413.D0/13311.D0 + 1.D0/65.D0*w(9)
      w(9)=y*w(9)
      w(9)= - 18099525303537844184.D0/4326075.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 69711966858561875321.D0/16286400.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 9141996864628144241.D0/252439200.D0 + 1.D0/121.D0*w(9)
      w(9)=w(9)*w(14)
      w(9)= - 25449264072856367.D0/15777450.D0 + w(9)
      w(9)=w(9)*w(11)
      w(9)= - 99594339247692787.D0/1867083075.D0 + w(9)
      w(11)=1.D0/7.D0*y
      w(9)=w(9)*w(11)
      w(9)= - 14596966172631413.D0/1867083075.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 374252596709956.D0/46621575.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 4999923546826426.D0/606080475.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 952416175990192.D0/112237125.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 301738637797453.D0/34534500.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 13006049925053.D0/1444170.D0 + w(9)
      w(9)=w(9)*w(6)
      w(9)= - 90220327106647.D0/29124095.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 226555711502.D0/70785.D0 + w(9)
      w(9)=w(9)*w(11)
      w(9)= - 202933637.D0/429.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 419455404416.D0/855855.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 1305804857198.D0/2567565.D0 + w(9)
      w(9)=w(9)*w(12)
      w(9)= - 63919699192.D0/2297295.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 39447325927.D0/1361360.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 2848817.D0/8008.D0 + 1.D0/85.D0*w(9)
      w(9)=y*w(9)
      w(9)= - 2133094.D0/5733.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 319815116.D0/819819.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 21832157.D0/53235.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 499426.D0/1155.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 7923808.D0/17325.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 2293672.D0/4725.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 1954871.D0/3780.D0 + w(9)
      w(9)=w(9)*w(6)
      w(9)= - 90367.D0/490.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 437212.D0/2205.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 48296.D0/225.D0 + w(9)
      w(9)=w(9)*w(10)
      w(9)= - 653.D0/25.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 2374.D0/81.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 2836.D0/81.D0 + w(9)
      w(9)=y*w(9)
      w(9)= - 488.D0/9.D0 + w(9)
      w(9)=y*w(9)
      w(11)= - 47.D0/9.D0 - 128.D0/25.D0*y
      w(11)=w(11)*w(7)
      w(11)= - 46.D0/423.D0 + w(11)
      w(11)=w(11)*w(6)
      w(11)= - 40.D0/1081.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 352.D0/9315.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 172.D0/4455.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 56.D0/1419.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 328.D0/8127.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 320.D0/7749.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 26.D0/615.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 76.D0/1755.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 296.D0/6669.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 32.D0/703.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 140.D0/2997.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 136.D0/2835.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 88.D0/1785.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 256.D0/5049.D0 + w(11)
      w(12)=2*y
      w(11)=w(11)*w(12)
      w(11)= - 31.D0/297.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 10.D0/93.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 464.D0/4185.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 448.D0/3915.D0 + w(11)
      w(11)=w(11)*w(6)
      w(11)= - 8.D0/203.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 208.D0/5103.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 400.D0/9477.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 128.D0/2925.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 92.D0/2025.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 88.D0/1863.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 112.D0/2277.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 320.D0/6237.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 152.D0/2835.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 16.D0/285.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 272.D0/4617.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 256.D0/4131.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 10.D0/153.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 28.D0/405.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 208.D0/2835.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 64.D0/819.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 88.D0/1053.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 80.D0/891.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 16.D0/165.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 128.D0/1215.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 28.D0/243.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 8.D0/63.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 80.D0/567.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 64.D0/405.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 8.D0/45.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 16.D0/81.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 16.D0/81.D0 + w(11)
      w(11)=w(11)*y**2
      w(11)= - 64.D0/27.D0 + w(11)
      w(11)=w(1)*w(11)
      w(9)=4*w(11) - 200.D0/9.D0 + w(9)
      w(9)=y*w(9)
      w(9)=88.D0/9.D0 + w(9)
      w(9)=w(1)*w(9)
      w(11)= - 10.D0/3.D0 - LL
      w(11)=LL*w(11)
      w(11)=1 + w(11)
      w(4)=w(9) + 512.D0/3.D0*w(11) + w(4)
      w(4)=w(1)*w(4)
      w(9)=68*LL
      w(11)= - w(17) + 2735.D0/27.D0 - w(9)
      w(13)=544*LL
      w(14)= - 1177.D0/147.D0*NF + 80800374295841332941899018891.D0/
     & 97717521785632645147851600.D0 - w(13)
      w(22)= - 2452.D0/20825.D0*NF + 3525956539430814929590131271.D0/
     & 288402408047874126304423125.D0 - w(16)
      w(22)=y*w(22)
      w(14)=1.D0/25.D0*w(14) + 8.D0/3.D0*w(22)
      w(14)=y*w(14)
      w(22)= - 2258.D0/987.D0*NF + 449410515185869958092052137.D0/
     & 1912837785891184640545200.D0 - 1088.D0/7.D0*LL
      w(14)=1.D0/7.D0*w(22) + w(14)
      w(14)=y*w(14)
      w(22)= - 1082.D0/1081.D0*NF + 720031077737617770641977697.D0/
     & 7044994067644983953903400.D0 - w(9)
      w(14)=1.D0/3.D0*w(22) + w(14)
      w(14)=y*w(14)
      w(22)=1088*LL
      w(23)= - 16576.D0/1035.D0*NF + 2379839720061693448378406341.D0/
     & 1463853452715730029353400.D0 - w(22)
      w(14)=1.D0/47.D0*w(23) + w(14)
      w(14)=y*w(14)
      w(23)= - 3964.D0/495.D0*NF + 1890641092547502122106283.D0/
     & 2339919201911333167125.D0 - w(13)
      w(14)=1.D0/23.D0*w(23) + w(14)
      w(14)=y*w(14)
      w(23)= - 7576.D0/473.D0*NF + 2625632990976866266174789.D0/
     & 1634868278706986423100.D0 - w(22)
      w(14)=1.D0/45.D0*w(23) + w(14)
      w(14)=y*w(14)
      w(23)=272*LL
      w(24)= - 3616.D0/903.D0*NF + 6892999697868398379038687.D0/
     & 17277585218153379244125.D0 - w(23)
      w(14)=1.D0/11.D0*w(24) + w(14)
      w(14)=y*w(14)
      w(24)= - 13792.D0/861.D0*NF + 1179046166383987911793787.D0/
     & 743695250219881084125.D0 - w(22)
      w(14)=1.D0/43.D0*w(24) + w(14)
      w(14)=y*w(14)
      w(24)= - 1642.D0/205.D0*NF + 3792128470652277611479061.D0/
     & 4816312096662087021000.D0 - w(13)
      w(14)=1.D0/21.D0*w(24) + w(14)
      w(14)=y*w(14)
      w(14)=w(14) - 3124.D0/7995.D0*NF + 28331226571971557603.D0/
     & 742815887294481000.D0 - 1088.D0/41.D0*LL
      w(14)=y*w(14)
      w(24)=136*LL
      w(25)= - 1484.D0/741.D0*NF + 22185494355252989875823.D0/
     & 114319365054620625900.D0 - w(24)
      w(14)=1.D0/5.D0*w(25) + w(14)
      w(14)=y*w(14)
      w(25)= - 11264.D0/703.D0*NF + 835666906657850453176943.D0/
     & 542284167566790148500.D0 - w(22)
      w(14)=1.D0/39.D0*w(25) + w(14)
      w(14)=y*w(14)
      w(25)= - 2668.D0/333.D0*NF + 196403158053532776189271.D0/
     & 256871447794795333500.D0 - w(13)
      w(14)=1.D0/19.D0*w(25) + w(14)
      w(14)=y*w(14)
      w(25)= - 5048.D0/315.D0*NF + 40664528991315158011.D0/
     & 26804909505874500.D0 - w(22)
      w(14)=1.D0/37.D0*w(25) + w(14)
      w(14)=y*w(14)
      w(25)= - 2384.D0/595.D0*NF + 584707924763740259.D0/
     & 1554475746198375.D0 - w(23)
      w(14)=1.D0/9.D0*w(25) + w(14)
      w(14)=y*w(14)
      w(25)= - 4496.D0/561.D0*NF + 5664704103786313441.D0/
     & 7594724359997775.D0 - w(13)
      w(14)=2.D0/35.D0*w(25) + w(14)
      w(14)=y*w(14)
      w(14)=w(14) - 529.D0/1122.D0*NF + 2802270677744838413.D0/
     & 64440085478769000.D0 - w(19)
      w(14)=y*w(14)
      w(25)= - 497.D0/31.D0*NF + 26483393387625457571.D0/
     & 18078214025718000.D0 - w(22)
      w(14)=1.D0/33.D0*w(25) + w(14)
      w(14)=y*w(14)
      w(25)=34*LL
      w(14)=w(14) - 233.D0/465.D0*NF + 31967600967246198187.D0/
     & 705050347003002000.D0 - w(25)
      w(14)=y*w(14)
      w(26)= - 6976.D0/435.D0*NF + 149813422077275609.D0/
     & 104295260569500.D0 - w(22)
      w(14)=1.D0/31.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 1628.D0/203.D0*NF + 1285881438117967.D0/1809203499675.D0
     &  - w(13)
      w(14)=1.D0/15.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 3032.D0/189.D0*NF + 40017355805241229.D0/28461072890250.D
     & 0 - w(22)
      w(14)=1.D0/29.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 704.D0/351.D0*NF + 4592055394870169.D0/26428139112375.D0
     &  - w(24)
      w(14)=2.D0/7.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 2608.D0/975.D0*NF + 7447604148739.D0/32534376875.D0 - 
     & 544.D0/3.D0*LL
      w(14)=2.D0/9.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 602.D0/75.D0*NF + 216623041666141.D0/319428427500.D0 - 
     & w(13)
      w(14)=1.D0/13.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 1108.D0/69.D0*NF + 10180753970819.D0/7606154556.D0 - 
     & w(22)
      w(14)=1.D0/25.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 508.D0/253.D0*NF + 9877137417599.D0/59863253450.D0 - 
     & w(24)
      w(14)=1.D0/3.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 3712.D0/231.D0*NF + 4857230439221.D0/3734380650.D0 - 
     & w(22)
      w(14)=1.D0/23.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 844.D0/105.D0*NF + 150489260983.D0/235030950.D0 - w(13)
      w(14)=1.D0/11.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 1528.D0/95.D0*NF + 3482154875279.D0/2764411650.D0 - 
     & w(22)
      w(14)=1.D0/21.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 344.D0/171.D0*NF + 18330323581.D0/118474785.D0 - w(24)
      w(14)=2.D0/5.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)= - 308.D0/153.D0*NF + 29653596157.D0/195270075.D0 - w(24)
      w(14)=8.D0/19.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(26)=1076638709.D0/7232225.D0 - w(24)
      w(26)=4*w(26) - 137.D0/17.D0*NF
      w(14)=1.D0/9.D0*w(26) + w(14)
      w(14)=y*w(14)
      w(14)=w(14) - 242.D0/255.D0*NF + 122766013.D0/1790100.D0 - w(21)
      w(14)=y*w(14)
      w(14)=w(14) - 106.D0/105.D0*NF + 4042578481.D0/56756700.D0 - w(9)
      w(14)=y*w(14)
      w(21)= - 1472.D0/91.D0*NF + 4558868983.D0/4099095.D0 - w(22)
      w(14)=1.D0/15.D0*w(21) + w(14)
      w(14)=y*w(14)
      w(21)= - 316.D0/39.D0*NF + 2855173081.D0/5270265.D0 - w(13)
      w(14)=1.D0/7.D0*w(21) + w(14)
      w(14)=y*w(14)
      w(21)= - 536.D0/33.D0*NF + 361464179.D0/343035.D0 - w(22)
      w(14)=1.D0/13.D0*w(21) + w(14)
      w(14)=y*w(14)
      w(21)= - 56.D0/55.D0*NF + 1107779.D0/17325.D0 - w(9)
      w(14)=4.D0/3.D0*w(21) + w(14)
      w(14)=y*w(14)
      w(21)= - 184.D0/45.D0*NF + 10546097.D0/42525.D0 - w(23)
      w(14)=4.D0/11.D0*w(21) + w(14)
      w(14)=y*w(14)
      w(13)= - 74.D0/9.D0*NF + 4088779.D0/8505.D0 - w(13)
      w(13)=1.D0/5.D0*w(13) + w(14)
      w(13)=y*w(13)
      w(14)= - 116.D0/7.D0*NF + 45743.D0/49.D0 - w(22)
      w(13)=1.D0/9.D0*w(14) + w(13)
      w(13)=y*w(13)
      w(13)=w(13) - 44.D0/21.D0*NF + 150839.D0/1323.D0 - w(24)
      w(13)=y*w(13)
      w(14)= - 16.D0/15.D0*NF + 38231.D0/675.D0 - w(9)
      w(13)=16.D0/7.D0*w(14) + w(13)
      w(13)=y*w(13)
      w(14)= - 22.D0/5.D0*NF + 5867.D0/25.D0 - w(23)
      w(13)=2.D0/3.D0*w(14) + w(13)
      w(13)=y*w(13)
      w(14)=911.D0/27.D0 - w(25)
      w(14)=4*w(14) - 7.D0/3.D0*NF
      w(13)=8.D0/5.D0*w(14) + w(13)
      w(13)=y*w(13)
      w(11)=4*w(11) + w(13)
      w(11)=y*w(11)
      w(13)= - w(15) + 101 - 68.D0/3.D0*LL
      w(11)=16*w(13) + w(11)
      w(11)=y*w(11)
      w(9)=28.D0/3.D0*NF - 5227.D0/9.D0 - w(9)
      w(9)=8*w(9) + w(11)
      w(9)=y*w(9)
      w(11)=11087.D0/3.D0 + 61672.D0/17.D0*y
      w(11)=y*w(11)
      w(11)=21241.D0/141.D0 + 1.D0/25.D0*w(11)
      w(11)=w(11)*w(7)
      w(11)=3388.D0/1081.D0 + w(11)
      w(11)=y*w(11)
      w(11)=1352.D0/423.D0 + w(11)
      w(11)=y*w(11)
      w(11)=37124.D0/11385.D0 + w(11)
      w(11)=y*w(11)
      w(11)=23612.D0/7095.D0 + w(11)
      w(11)=y*w(11)
      w(11)=33752.D0/9933.D0 + w(11)
      w(11)=y*w(11)
      w(11)=128504.D0/37023.D0 + w(11)
      w(11)=y*w(11)
      w(11)=1018.D0/287.D0 + w(11)
      w(11)=y*w(11)
      w(11)=28994.D0/7995.D0 + w(11)
      w(11)=y*w(11)
      w(11)=13744.D0/3705.D0 + w(11)
      w(11)=y*w(11)
      w(11)=34696.D0/9139.D0 + w(11)
      w(11)=y*w(11)
      w(11)=24596.D0/6327.D0 + w(11)
      w(11)=y*w(11)
      w(11)=9284.D0/2331.D0 + w(11)
      w(11)=y*w(11)
      w(11)=7288.D0/1785.D0 + w(11)
      w(11)=y*w(11)
      w(11)=82232.D0/19635.D0 + w(11)
      w(11)=w(11)*w(12)
      w(11)=4823.D0/561.D0 + w(11)
      w(11)=y*w(11)
      w(11)=3011.D0/341.D0 + w(11)
      w(11)=y*w(11)
      w(11)=844.D0/93.D0 + w(11)
      w(11)=y*w(11)
      w(11)=125872.D0/13485.D0 + w(11)
      w(11)=y*w(11)
      w(11)=9752.D0/1015.D0 + w(11)
      w(11)=y*w(11)
      w(11)=54248.D0/5481.D0 + w(11)
      w(11)=y*w(11)
      w(11)=25072.D0/2457.D0 + w(11)
      w(11)=y*w(11)
      w(11)=1232.D0/117.D0 + w(11)
      w(11)=y*w(11)
      w(11)=10604.D0/975.D0 + w(11)
      w(11)=y*w(11)
      w(11)=19396.D0/1725.D0 + w(11)
      w(11)=y*w(11)
      w(11)=128.D0/11.D0 + w(11)
      w(11)=y*w(11)
      w(11)=64048.D0/5313.D0 + w(11)
      w(11)=y*w(11)
      w(11)=2888.D0/231.D0 + w(11)
      w(11)=y*w(11)
      w(11)=8632.D0/665.D0 + w(11)
      w(11)=y*w(11)
      w(11)=11536.D0/855.D0 + w(11)
      w(11)=y*w(11)
      w(11)=40816.D0/2907.D0 + w(11)
      w(11)=y*w(11)
      w(11)=746.D0/51.D0 + w(11)
      w(11)=y*w(11)
      w(11)=778.D0/51.D0 + w(11)
      w(11)=y*w(11)
      w(11)=1672.D0/105.D0 + w(11)
      w(11)=y*w(11)
      w(11)=7568.D0/455.D0 + w(11)
      w(11)=y*w(11)
      w(11)=4744.D0/273.D0 + w(11)
      w(11)=y*w(11)
      w(11)=7784.D0/429.D0 + w(11)
      w(11)=y*w(11)
      w(11)=208.D0/11.D0 + w(11)
      w(11)=y*w(11)
      w(11)=9712.D0/495.D0 + w(11)
      w(11)=y*w(11)
      w(11)=908.D0/45.D0 + w(11)
      w(11)=y*w(11)
      w(11)=428.D0/21.D0 + w(11)
      w(11)=y*w(11)
      w(11)=416.D0/21.D0 + w(11)
      w(11)=y*w(11)
      w(11)=368.D0/21.D0 + w(11)
      w(11)=y*w(11)
      w(11)=56.D0/5.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 88.D0/15.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 176.D0/3.D0 + w(11)
      w(11)=y*w(11)
      w(11)= - 304 + w(11)
      w(11)=y*w(11)
      w(11)=1504 + w(11)
      w(11)=w(11)*w(6)
      w(11)= - 208 + w(11)
      w(11)=w(1)*w(11)
      w(13)= - w(17) + 154.D0/3.D0 - 17*LL
      w(9)=2*w(11) + 64*w(13) + w(9)
      w(11)=y - 2
      w(13)=z2*w(11)
      w(9)=1.D0/9.D0*w(9) + 368.D0/5.D0*w(13)
      w(9)=z2*w(9)
      w(13)=17 - 5.D0/3.D0*LL
      w(13)=w(13)*w(16)
      w(14)= - 224.D0/27.D0 - LL
      w(14)=NF*w(14)
      w(13)=40.D0/9.D0*w(14) + 105371.D0/243.D0 + w(13)
      w(11)=B4*w(11)
      w(3)= - 32.D0/27.D0*w(11) + 1.D0/3.D0*w(9) + 1.D0/9.D0*w(4) + 8.D0
     & /9.D0*w(13) + w(3)
      w(3)=as*w(3)
      w(4)= - 8778163403480727116759.D0/193690281515374794150.D0 + 
     & w(19)
      w(4)=w(4)*w(8)
      w(8)= - 8716182513395807182631.D0/387380563030749588300.D0 + 
     & w(18)
      w(4)=1.D0/25.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 176590545044626636919.D0/3952862888068873350.D0 + w(19)
      w(4)=1.D0/49.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 175272924081937012469.D0/63245806209101973600.D0 + w(5)
      w(4)=1.D0/3.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 3700580183631479227.D0/84103465703593050.D0 + w(19)
      w(4)=1.D0/47.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 3671326804256316427.D0/168206931407186100.D0 + w(18)
      w(4)=1.D0/23.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 3641423349783927787.D0/84103465703593050.D0 + w(19)
      w(4)=1.D0/45.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 328258206486023417.D0/30583078437670200.D0 + w(16)
      w(4)=1.D0/11.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 7567750440859019.D0/177808595567850.D0 + w(19)
      w(4)=1.D0/43.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 52500096831165533.D0/2489320337949900.D0 + w(18)
      w(4)=1.D0/21.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 1268643311941813.D0/30357565096950.D0 + w(19)
      w(4)=1.D0/41.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 1256500285903033.D0/242860520775600.D0 + w(2)
      w(4)=1.D0/5.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 1244045900222233.D0/30357565096950.D0 + w(19)
      w(4)=1.D0/39.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 1231263767549833.D0/60715130193900.D0 + w(18)
      w(4)=1.D0/19.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 32922599238709.D0/820474732350.D0 + w(19)
      w(4)=1.D0/37.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 32557943802109.D0/3281898929400.D0 + w(16)
      w(4)=1.D0/9.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 32182869638749.D0/820474732350.D0 + w(19)
      w(4)=1.D0/35.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 31796763882349.D0/1640949464700.D0 + w(18)
      w(4)=1.D0/17.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 345388537466639.D0/9025222055850.D0 + w(19)
      w(4)=1.D0/33.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(4)=w(4) - 170437963219357.D0/144403552893600.D0 + LL
      w(4)=y*w(4)
      w(8)= - 5422866892147.D0/145568097675.D0 + w(19)
      w(4)=1.D0/31.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 5345230573387.D0/291136195350.D0 + w(18)
      w(4)=1.D0/15.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 181548866903.D0/5019589575.D0 + w(19)
      w(4)=1.D0/29.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 178680530003.D0/20078358300.D0 + w(16)
      w(4)=1.D0/7.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 2169209363.D0/185910725.D0 + 32.D0/3.D0*LL
      w(4)=1.D0/9.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 6393221489.D0/371821450.D0 + w(18)
      w(4)=1.D0/13.D0*w(8) + w(4)
      w(4)=w(4)*w(10)
      w(8)= - 5577101.D0/7436429.D0 + 32.D0/45.D0*LL
      w(4)=1.D0/5.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 738035777.D0/178474296.D0 + w(2)
      w(4)=1.D0/27.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 31413751.D0/969969.D0 + w(19)
      w(4)=1.D0/207.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 30708319.D0/1939938.D0 + w(18)
      w(4)=1.D0/99.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 9989765.D0/323323.D0 + w(19)
      w(4)=1.D0/189.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 16218511.D0/6466460.D0 + w(20)
      w(4)=1.D0/15.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 2489167.D0/85085.D0 + w(19)
      w(4)=1.D0/171.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 21721823.D0/1531530.D0 + w(18)
      w(4)=1.D0/81.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 1235359.D0/45045.D0 + w(19)
      w(4)=1.D0/153.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 595157.D0/360360.D0 + w(5)
      w(4)=1.D0/9.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 571133.D0/45045.D0 + w(18)
      w(4)=2.D0/135.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 545393.D0/45045.D0 + w(18)
      w(4)=1.D0/63.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 39821.D0/3465.D0 + w(18)
      w(4)=2.D0/117.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 37511.D0/6930.D0 + w(16)
      w(4)=1.D0/27.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 3181.D0/315.D0 + w(18)
      w(4)=2.D0/99.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 2929.D0/315.D0 + w(18)
      w(4)=1.D0/45.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 883.D0/105.D0 + w(18)
      w(4)=2.D0/81.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 389.D0/210.D0 + w(2)
      w(4)=1.D0/9.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 47.D0/15.D0 + w(16)
      w(4)=4.D0/63.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 37.D0/15.D0 + w(16)
      w(4)=2.D0/27.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 1.D0/3.D0 + 8.D0/5.D0*LL
      w(4)=4.D0/9.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(8)= - 1.D0/3.D0 + w(2)
      w(4)=2.D0/9.D0*w(8) + w(4)
      w(4)=y*w(4)
      w(2)=1.D0/3.D0 + w(2)
      w(2)=8.D0/27.D0*w(2) + w(4)
      w(2)=y*w(2)
      w(4)=1.D0/3.D0 + LL
      w(2)=16.D0/9.D0*w(4) + w(2)
      w(2)=y*w(2)
      w(4)=5.D0/3.D0 + w(5)
      w(2)=w(3) + 16.D0/9.D0*w(4) + w(2)
      w(2)=w(2)*as**2
      w(3)=551.D0/9.D0 + 1496.D0/25.D0*y
      w(3)=w(3)*w(7)
      w(3)=541.D0/423.D0 + w(3)
      w(3)=w(3)*w(6)
      w(3)=472.D0/1081.D0 + w(3)
      w(3)=y*w(3)
      w(3)=4168.D0/9315.D0 + w(3)
      w(3)=y*w(3)
      w(3)=2044.D0/4455.D0 + w(3)
      w(3)=y*w(3)
      w(3)=668.D0/1419.D0 + w(3)
      w(3)=y*w(3)
      w(3)=3928.D0/8127.D0 + w(3)
      w(3)=y*w(3)
      w(3)=3848.D0/7749.D0 + w(3)
      w(3)=y*w(3)
      w(3)=314.D0/615.D0 + w(3)
      w(3)=y*w(3)
      w(3)=922.D0/1755.D0 + w(3)
      w(3)=y*w(3)
      w(3)=3608.D0/6669.D0 + w(3)
      w(3)=y*w(3)
      w(3)=392.D0/703.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1724.D0/2997.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1684.D0/2835.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1096.D0/1785.D0 + w(3)
      w(3)=y*w(3)
      w(3)=3208.D0/5049.D0 + w(3)
      w(3)=w(3)*w(12)
      w(3)=391.D0/297.D0 + w(3)
      w(3)=y*w(3)
      w(3)=127.D0/93.D0 + w(3)
      w(3)=y*w(3)
      w(3)=5936.D0/4185.D0 + w(3)
      w(3)=y*w(3)
      w(3)=5776.D0/3915.D0 + w(3)
      w(3)=w(3)*w(6)
      w(3)=104.D0/203.D0 + w(3)
      w(3)=y*w(3)
      w(3)=2728.D0/5103.D0 + w(3)
      w(3)=y*w(3)
      w(3)=5296.D0/9477.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1712.D0/2925.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1244.D0/2025.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1204.D0/1863.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1552.D0/2277.D0 + w(3)
      w(3)=y*w(3)
      w(3)=4496.D0/6237.D0 + w(3)
      w(3)=y*w(3)
      w(3)=2168.D0/2835.D0 + w(3)
      w(3)=y*w(3)
      w(3)=232.D0/285.D0 + w(3)
      w(3)=y*w(3)
      w(3)=4016.D0/4617.D0 + w(3)
      w(3)=y*w(3)
      w(3)=3856.D0/4131.D0 + w(3)
      w(3)=y*w(3)
      w(3)=154.D0/153.D0 + w(3)
      w(3)=y*w(3)
      w(3)=442.D0/405.D0 + w(3)
      w(3)=y*w(3)
      w(3)=3376.D0/2835.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1072.D0/819.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1528.D0/1053.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1448.D0/891.D0 + w(3)
      w(3)=y*w(3)
      w(3)=304.D0/165.D0 + w(3)
      w(3)=y*w(3)
      w(3)=2576.D0/1215.D0 + w(3)
      w(3)=y*w(3)
      w(3)=604.D0/243.D0 + w(3)
      w(3)=y*w(3)
      w(3)=188.D0/63.D0 + w(3)
      w(3)=y*w(3)
      w(3)=2096.D0/567.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1936.D0/405.D0 + w(3)
      w(3)=y*w(3)
      w(3)=296.D0/45.D0 + w(3)
      w(3)=y*w(3)
      w(3)=808.D0/81.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1456.D0/81.D0 + w(3)
      w(3)=w(3)*w(6)
      w(3)=16 + w(3)
      w(3)=y*w(3)
      w(4)= - 481 - 56.D0/3.D0*NF
      w(3)=2.D0/27.D0*w(4) + w(3)
      w(3)=y*w(3)
      w(4)=767 + 56*NF
      w(3)=4.D0/81.D0*w(4) + w(3)
      w(3)=z3*w(3)*as**3
      w(2)=w(2) + 2*w(3)
*
      ANSREG2 = w(2)
*
      RETURN
      END            
      REAL*8 FUNCTION BNSREG(z,nf,as,LL)
*     ----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of aqqQNS3: old file, not active
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(116),G
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
*      LL=10.0D0 !!
      w(1)=Hr1(-1)
      w(2)=1.0D0/(9.0D0 + 9.0D0*z)
      w(3)=1.0D0/(27.0D0 + 27.0D0*z)
      w(4)=Hr2(-1,-1)
      w(5)=Hr4(-1,-1,-1,0)
      w(6)=Hr3(-1,-1,0)
      w(7)=Hr4(-1,-1,0,0)
      w(8)=Hr2(-1,0)
      w(9)=1.0D0/(81.00 + 81.0D0*z)
      w(10)=Hr4(-1,0,-1,0)
      w(11)=Hr3(-1,0,0)
      w(12)=Hr4(-1,0,0,0)
      w(13)=Hr4(-1,0,0,1)
      w(14)=Hr3(-1,0,1)
      w(15)=Hr4(-1,0,1,1)
      w(16)=Hr1(0)
      w(17)=1.0D0/( - 81.0D0 + 81.0D0*z)
      w(18)=1.0D0/( - 27.0D0 + 27.0D0*z**2)
      w(19)=1.0D0/( - 9.0D0 + 9.0D0*z)
      w(20)=Hr2(0,-1)
      w(21)=Hr4(0,-1,-1,0)
      w(22)=Hr3(0,-1,0)
      w(23)=Hr4(0,-1,0,0)
      w(24)=Hr4(0,-1,0,1)
      w(25)=Hr2(0,0)
      w(26)=1.0D0/( - 81.0D0 + 81.0D0*z**2)
      w(27)=1.0D0/( - 9.0D0 + 9.0D0*z**2)
      w(28)=1.0D0/( - 3.0D0 + 3.0D0*z)
      w(29)=Hr4(0,0,-1,0)
      w(30)=Hr3(0,0,0)
      w(31)=1.0D0/( - 27.0D0 + 27.0D0*z)
      w(32)=Hr4(0,0,0,0)
      w(33)=Hr4(0,0,0,1)
      w(34)=Hr3(0,0,1)
      w(35)=Hr4(0,0,1,0)
      w(36)=Hr4(0,0,1,1)
      w(37)=Hr2(0,1)
      w(38)=Hr3(0,1,0)
      w(39)=Hr4(0,1,0,0)
      w(40)=Hr4(0,1,0,1)
      w(41)=Hr3(0,1,1)
      w(42)=Hr4(0,1,1,0)
      w(43)=Hr4(0,1,1,1)
      w(44)=Hr1(1)
      w(45)=Hr2(1,0)
      w(46)=Hr3(1,0,0)
      w(47)=Hr4(1,0,0,0)
      w(48)=Hr4(1,0,0,1)
      w(49)=Hr3(1,0,1)
      w(50)=Hr4(1,0,1,0)
      w(51)=Hr4(1,0,1,1)
      w(52)=Hr2(1,1)
      w(53)=Hr3(1,1,0)
      w(54)=Hr4(1,1,0,0)
      w(55)=Hr4(1,1,0,1)
      w(56)=Hr4(1,1,1,0)
      w(57)=1.0D0/(15.0D0 + 15.0D0*z)
      w(58)=1.0D0/(45.0D0 + 45.0D0*z)
      w(59)=CF**2
      w(60)=2*w(59)
      w(61)=CA*CF
      w(62)=w(60) - w(61)
      w(63)=18 + 29*z
      w(63)=w(63)*z
      w(63)=w(63) + 29
      w(64)=as**3
      w(65)=w(64)*TF
      w(66)=z2*w(65)
      w(63)= - w(3)*w(66)*w(63)*w(62)
      w(67)=z**2
      w(68)=w(67) + 1
      w(69)= - w(68)*w(62)
      w(70)=w(2)*w(65)
      w(71)=w(69)*w(70)
      w(72)=z3*w(71)
      w(63)=w(63) - 4*w(72)
      w(63)=32*w(63)
      w(72)=64*w(70)
      w(72)=w(69)*w(72)
      w(73)=z2*w(72)
      w(71)=128*w(71)
      w(74)=z + 1
      w(75)=w(60)*w(74)
      w(76)=w(74)*w(61)
      w(75)=w(75) - w(76)
      w(77)=64.D0/9.D0*w(65)
      w(78)=w(75)*w(77)
      w(79)=174 + 199*z
      w(79)=w(79)*z
      w(79)=w(79) + 199
      w(79)= - w(65)*w(9)*w(79)*w(62)
      w(80)=w(70)*z2
      w(69)=w(69)*w(80)
      w(79)=2*w(79) + 5*w(69)
      w(79)=16*w(79)
      w(81)=18 + 19*z
      w(81)=w(81)*z
      w(81)=w(81) + 19
      w(82)=w(65)*w(3)
      w(83)=32*w(82)
      w(81)=w(83)*w(62)*w(81)
      w(83)=4*z
      w(84)=w(83) + 3
      w(84)=w(84)*z
      w(84)=w(84) + 4
      w(62)=256*w(82)*w(84)*w(62)
      w(82)= - 473 + 1449*z
      w(82)=z*w(82)
      w(84)=5 + 6*LL
      w(85)=LL*w(84)
      w(82)= - 27*w(85) - 868 + w(82)
      w(85)=w(60)*w(17)
      w(82)=w(82)*w(85)
      w(86)=2086 - 1761*z
      w(86)=z*w(86)
      w(86)=1344*LL + 907 + w(86)
      w(86)=w(17)*w(86)*w(61)
      w(82)=w(82) + w(86)
      w(82)=as*w(82)
      w(84)=CF*w(19)*w(84)
      w(82)=2*w(84) + w(82)
      w(84)=as**2
      w(82)=w(82)*w(84)
      w(86)=52 - 55*z
      w(86)=z*w(86)
      w(87)=LL**2
      w(88)=26 - 43*z
      w(88)=z*w(88)
      w(88)= - 30*LL - 39 + w(88)
      w(88)=NF*w(88)
      w(86)=w(88) + 18*w(87) - 53 + w(86)
      w(86)=w(17)*w(86)
      w(87)=NF + 2
      w(87)=w(87)*z2
      w(68)= - w(19)*w(68)*w(87)
      w(68)=2*w(86) + w(68)
      w(86)=4*w(65)
      w(68)=w(86)*CF*w(68)
      w(68)=w(82) + w(68)
      w(68)=TF*w(68)
      w(82)=4 + z
      w(82)=z*w(82)
      w(82)=1 + w(82)
      w(82)=w(82)*w(59)
      w(88)=7*z
      w(89)=4 + w(88)
      w(89)=z*w(89)
      w(89)= - 2 + w(89)
      w(89)=w(89)*w(61)
      w(82)=w(82) + w(89)
      w(82)=z3*w(82)*w(70)
      w(89)=2*z
      w(90)=27 + 70*z
      w(90)=w(90)*w(89)
      w(90)= - 33 + w(90)
      w(90)=z*w(90)
      w(91)=w(74)*LL
      w(92)=72*w(91)
      w(90)= - w(92) - 107 + w(90)
      w(90)=w(90)*w(60)
      w(93)= - 27 + w(88)
      w(93)=z*w(93)
      w(93)= - 21 + w(93)
      w(93)=z*w(93)
      w(93)= - 144*w(91) + 173 + w(93)
      w(93)=w(93)*w(61)
      w(90)=w(90) + w(93)
      w(90)=w(18)*w(90)*w(66)
      w(68)=w(90) + w(68) + 8*w(82)
      w(68)=4*w(68)
      w(69)=32*w(69)
      w(82)=87 + 629*z
      w(82)=z*w(82)
      w(82)=591 + w(82)
      w(82)=z*w(82)
      w(82)= - 324*w(91) - 659 + w(82)
      w(82)=w(82)*w(59)
      w(90)=31*w(74) + 6*w(91)
      w(93)=18*LL
      w(90)=w(90)*w(93)
      w(94)= - 41 - 426*z
      w(94)=z*w(94)
      w(94)=308 + w(94)
      w(94)=z*w(94)
      w(90)=w(90) + 819 + w(94)
      w(90)=w(90)*w(61)
      w(82)=w(82) + w(90)
      w(82)=w(64)*w(26)*w(82)
      w(90)=11*z
      w(94)=w(90) - 12
      w(94)=w(94)*z
      w(95)= - w(93) - 31 - w(94)
      w(95)=NF*w(95)
      w(94)= - 16 - w(94)
      w(94)=2*w(94) + w(95)
      w(86)=w(17)*w(94)*w(86)
      w(84)=w(28)*w(84)
      w(84)=w(86) + w(84)
      w(84)=CF*w(84)
      w(82)=w(82) + w(84)
      w(82)=TF*w(82)
      w(84)=2*w(61)
      w(83)=9 - w(83)
      w(83)=z*w(83)
      w(83)= - 6 + w(83)
      w(83)=z*w(83)
      w(83)= - 5 + w(83)
      w(66)=w(27)*w(83)*w(66)*w(84)
      w(83)=2 + 15*z
      w(83)=z*w(83)
      w(83)=15 + w(83)
      w(80)=w(83)*w(59)*w(80)
      w(66)=w(80) + w(82) + w(66)
      w(66)=8*w(66)
      w(80)= - 15 - w(90)
      w(80)=z*w(80)
      w(80)=51 + w(80)
      w(80)=z*w(80)
      w(80)=18*w(91) - 25 + w(80)
      w(82)=4*w(59)
      w(80)=w(80)*w(82)
      w(83)=11 - 75*z
      w(83)=z*w(83)
      w(83)=61 + w(83)
      w(83)=z*w(83)
      w(83)=w(92) + 135 + w(83)
      w(83)=w(83)*w(61)
      w(80)=w(80) + w(83)
      w(83)=w(65)*w(18)
      w(80)=w(80)*w(83)
      w(86)= - 2 - w(67)
      w(92)= - 5 - w(67)
      w(92)=NF*w(92)
      w(86)=2*w(86) + w(92)
      w(92)=CF*w(64)*TF**2
      w(86)=w(31)*w(86)*w(92)
      w(80)=4*w(86) + w(80)
      w(80)=8*w(80)
      w(86)=w(74)*z
      w(86)=w(86) + 1
      w(86)=w(86)*w(61)
      w(94)=5*z
      w(95)= - 26 - w(94)
      w(95)=z*w(95)
      w(95)= - 5 + w(95)
      w(95)=w(95)*w(59)
      w(95)=w(95) - 8*w(86)
      w(96)=16*w(70)
      w(95)=w(95)*w(96)
      w(97)= - 6 - w(88)
      w(97)=z*w(97)
      w(97)= - 7 + w(97)
      w(97)=w(97)*w(60)
      w(98)= - 4 + w(94)
      w(98)=z*w(98)
      w(98)= - 1 + w(98)
      w(98)=w(98)*w(61)
      w(97)=w(97) + w(98)
      w(96)=w(97)*w(96)
      w(97)=36*w(91)
      w(98)= - 9 - 73*z
      w(98)=z*w(98)
      w(98)=33 + w(98)
      w(98)=z*w(98)
      w(98)=w(97) + 49 + w(98)
      w(98)=w(98)*w(82)
      w(99)=30 + 89*z
      w(99)=z*w(99)
      w(99)= - 9 + w(99)
      w(99)=z*w(99)
      w(97)=w(97) - 110 + w(99)
      w(97)=w(97)*w(61)
      w(97)=w(98) + w(97)
      w(97)=8*w(97)*w(83)
      w(98)=w(74)*w(59)
      w(99)=4.D0/3.D0*w(98)
      w(100)= - z*w(61)
      w(100)= - w(99) + w(100)
      w(101)=32.D0/3.D0*w(65)
      w(100)=w(100)*w(101)
      w(94)=2 + w(94)
      w(94)=z*w(94)
      w(94)=5 + w(94)
      w(94)=w(94)*w(59)
      w(86)=w(94) - 4*w(86)
      w(70)=32*w(86)*w(70)
      w(86)= - 1 - 5*w(67)
      w(86)=w(19)*w(86)*w(59)
      w(94)=2 - z
      w(94)=w(94)*w(61)
      w(86)=w(86) + 2.D0/9.D0*w(94)
      w(86)=z2*w(86)
      w(94)=10 + 3*LL
      w(102)=w(94)*w(93)
      w(103)=39 - 259*z
      w(103)=z*w(103)
      w(102)=w(102) + 220 + w(103)
      w(102)=w(102)*w(85)
      w(103)=1.D0/9.D0*w(61)
      w(104)=29 + 8*z
      w(104)=w(104)*w(103)
      w(86)=w(86) + w(102) + w(104)
      w(102)=16*w(65)
      w(86)=w(86)*w(102)
      w(104)=9 - w(90)
      w(104)=z*w(104)
      w(104)=9*LL + 2 + w(104)
      w(104)=w(31)*w(104)*w(59)
      w(105)=1.D0/3.D0*w(61)
      w(106)= - 1 + 1.D0/3.D0*z
      w(106)=w(106)*w(105)
      w(104)=8*w(104) + w(106)
      w(104)=w(104)*w(102)
      w(106)= - 5*w(98) + 2*w(76)
      w(107)=32.D0/9.D0*w(65)
      w(106)=w(106)*w(107)
      w(108)=w(89) - 1
      w(109)=w(108)*w(61)
      w(109)= - w(98) + w(109)
      w(109)=w(109)*w(107)
      w(110)=w(59)*w(77)
      w(108)= - w(108)*w(105)
      w(108)= - w(98) + w(108)
      w(101)=w(108)*w(101)
      w(107)= - w(98)*w(107)
      w(108)= - 328 + 337*z
      w(108)=w(108)*w(59)
      w(111)=91 - 103*z
      w(111)=w(111)*w(61)
      w(108)=1.D0/3.D0*w(108) + w(111)
      w(111)=z - 1
      w(112)=w(111)*w(59)
      w(113)= - 14.D0/3.D0 - z
      w(113)=w(113)*w(61)
      w(112)= - 8.D0/3.D0*w(112) + w(113)
      w(112)=z2*w(112)
      w(113)= - 8*w(98) + w(76)
      w(113)=z3*w(113)
      w(108)=2.D0/3.D0*w(113) + 1.D0/3.D0*w(108) + w(112)
      w(112)=8.D0/3.D0*w(65)
      w(108)=w(108)*w(112)
      w(67)=1 - 3*w(67)
      w(67)=w(28)*w(67)*w(60)
      w(67)=w(67) + 10.D0/9.D0*w(76)
      w(67)=z2*w(67)
      w(113)=69 - 263*z
      w(113)=z*w(113)
      w(94)=LL*w(94)
      w(94)=36*w(94) + 194 + w(113)
      w(85)=w(94)*w(85)
      w(94)=4.D0/9.D0 + z
      w(94)=w(94)*w(61)
      w(67)=w(85) + w(94) + w(67)
      w(67)=8*w(67)*w(65)
      w(85)=42 - 61*z
      w(85)=z*w(85)
      w(85)=72*LL + 19 + w(85)
      w(85)=w(31)*w(85)*w(60)
      w(94)= - 11 + 14*z
      w(94)=z*w(94)
      w(93)= - w(93) - 3 + w(94)
      w(93)=w(19)*w(93)*w(61)
      w(85)=w(85) + w(93)
      w(85)=w(85)*w(102)
      w(93)= - 16*w(98) + 7*w(76)
      w(94)=16.D0/9.D0*w(65)
      w(93)=w(93)*w(94)
      w(113)= - w(76) + 2.D0/3.D0*w(98)
      w(114)=16.D0/3.D0*w(65)
      w(113)=w(113)*w(114)
      w(88)=10 + w(88)
      w(88)=8.D0/9.D0*w(88)*w(61)*w(65)
      w(99)= - w(99) + w(76)
      w(99)=w(99)*w(114)
      w(114)= - w(98) + w(76)
      w(114)=z2*w(114)
      w(115)= - 4.D0/3.D0 + z
      w(115)=w(115)*w(82)
      w(116)=14 - w(90)
      w(105)=w(116)*w(105)
      w(105)=16.D0/3.D0*w(114) + w(115) + w(105)
      w(105)=w(105)*w(112)
      w(82)= - w(111)*w(82)
      w(89)= - 5 - w(89)
      w(89)=w(89)*w(61)
      w(82)=w(82) + w(89)
      w(82)=w(82)*w(94)
      w(89)=4*w(98) - w(76)
      w(89)=w(89)*w(94)
      w(76)= - w(76)*w(77)
      w(77)= - 16 + 7*NF
      w(74)=w(92)*w(74)*w(77)
      w(77)=288*w(91)
      w(91)= - 677 + 187*z
      w(91)=z*w(91)
      w(91)=293 + w(91)
      w(91)=z*w(91)
      w(91)= - w(77) + 197 + w(91)
      w(91)=w(91)*w(60)
      w(92)=1033 - 35*z
      w(92)=z*w(92)
      w(92)= - 445 + w(92)
      w(92)=z*w(92)
      w(77)=w(77) - 553 + w(92)
      w(77)=w(77)*w(61)
      w(77)=w(91) + w(77)
      w(77)=w(77)*w(83)
      w(74)=w(77) + 16.D0/27.D0*w(74)
      w(74)=z3*w(74)
      w(77)=1925 + 2689*z
      w(60)=w(77)*w(60)
      w(77)=18473 - 6079*z
      w(77)=w(77)*w(103)
      w(60)=w(60) + w(77)
      w(77)=1458 - 653*z
      w(77)=z*w(77)
      w(77)=319 + w(77)
      w(77)=w(9)*w(77)
      w(83)=34 - z
      w(83)=z*w(83)
      w(83)= - 1 + w(83)
      w(83)=z2*w(58)*w(83)
      w(77)=8*w(83) + w(77)
      w(59)=w(59)*w(77)
      w(77)=124 + 721*z
      w(77)=z*w(77)
      w(77)= - 149 + w(77)
      w(77)=w(9)*w(77)*w(84)
      w(59)=w(77) + w(59)
      w(59)=z2*w(59)
      w(59)=1.D0/81.D0*w(60) + 2*w(59)
      w(59)=w(59)*w(64)
      w(60)=w(90) - 1
      w(60)= - w(60)*w(87)
      w(64)= - 283 - 1453*z
      w(64)=NF*w(64)
      w(64)=2*w(64) + 617 - 1831*z
      w(60)=1.D0/27.D0*w(64) + w(60)
      w(60)=w(65)*CF*w(60)
      w(59)=w(59) + 8.D0/27.D0*w(60)
      w(59)=TF*w(59)
      w(60)= - B4*w(75)*w(112)
      w(64)= - 23 - 12*z
      w(64)=z*w(64)
      w(64)= - 5 + w(64)
      w(61)=w(57)*w(64)*z2**2*w(61)*w(102)
      w(59)=w(61) + w(60) + w(59) + w(74)
      w(59)=2*w(59)

      G = w(59) + w(62)*w(14) + w(63)*w(1) + w(66)*w(25) + w(67)*w(45)
     &  + w(68)*w(16) + w(69)*w(20) + w(70)*w(36) + w(71)*w(5) + w(71)*
     & w(15) - w(72)*w(7) - w(72)*w(10) + w(72)*w(12) - w(72)*w(13) - 
     & w(72)*w(21) + w(72)*w(23) - w(72)*w(24) + w(72)*w(29) + w(73)*
     & w(4) + w(76)*w(55) + w(78)*w(6) - w(78)*w(56) + w(79)*w(8) + 
     & w(80)*w(30) - w(81)*w(11) - w(81)*w(22) + w(82)*w(53) + w(85)*
     & w(46) + w(86)*w(37) + w(88)*w(49) + w(89)*w(54) + w(93)*w(47) + 
     & w(95)*w(32) + w(96)*w(33) + w(97)*w(34) + w(99)*w(50) + w(100)*
     & w(35) + w(101)*w(42) + w(104)*w(38) + w(105)*w(52) + w(106)*
     & w(39) + w(107)*w(43) + w(108)*w(44) + w(109)*w(40) + w(110)*
     & w(41) + w(113)*w(48) + w(113)*w(51)
*
      BNSREG = G
*
      RETURN
      END      
      REAL*8 FUNCTION ANSREG(z,nf,as,LL)
*     ----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of aqqQNS3: combination of cases
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(116),G
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4,ANSR!,ANSREG1,ANSREG2
!      EXTERNAL ANSREG1,ANSREG2
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
*      LL=10.0D0
      IF(z.LE.0.5D0) ANSR=ANSREG1(z,nf,as,LL)
      IF(z.GT.0.5D0) ANSR=ANSREG2(z,nf,as,LL)
*
      ANSREG = ANSR
*
      RETURN
      END      
      REAL*8 FUNCTION ANSPLU2(z,nf,as,LL)
*     -----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     plus part of aqqQNS3, second contribution
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(4),G
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      ANSPLU2=-4.0D0*as**3*CA*CF*TF*(8.0D0*z2 + Hr1(0)**2 
     &- 8.0D0*Hr2(0,1))/(3.0D0*(1.0D0 - z)**2)
*
      RETURN
      END      
      REAL*8 FUNCTION ANSPLU1(z,nf,as,LL)
*     -----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     plus part of aqqQNS3; first contribution
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(27)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*      LL=10.0D0
*
      w(1)=1.0D0/(1.0D0 - z)
      w(2)=Hr1(0)
      w(3)=Hr1(1)
      w(4)=Hr2(0,1)
      w(5)=Hr3(0,0,1)
      w(6)=Hr3(0,1,1)
      w(7)=Hr4(0,0,0,1)
      w(8)=Hr4(0,0,1,1)
      w(9)=Hr4(0,1,1,1)
      w(10)=1.D0/3.D0*w(2)
      w(11)=2*LL
      w(12)=w(10) + w(11)
      w(13)= - 8.D0/9.D0*w(3) + 23.D0/9.D0 + w(12)
      w(14)=2*w(2)
      w(13)=w(13)*w(14)
      w(15)= - 1 + 16.D0/9.D0*w(3)
      w(15)=w(15)*w(3)
      w(15)=w(15) - 4.D0/9.D0*w(4)
      w(13)=w(13) - 4.D0/9.D0 + w(15)
      w(13)=w(4)*w(13)
      w(16)=2.D0/3.D0*w(2)
      w(17)=4*LL
      w(18)=w(17) + 31.D0/9.D0
      w(19)= - 7.D0/3.D0*w(2) - 10.D0/3.D0*w(3) - w(18)
      w(19)=w(19)*w(16)
      w(15)=64.D0/15.D0*z2 + w(19) - 370.D0/81.D0 - w(15)
      w(15)=z2*w(15)
      w(19)=4*w(3)
      w(20)=20.D0/9.D0*w(2) + 3 - w(19)
      w(20)=w(6)*w(20)
      w(13)=w(20) + w(13) + w(15)
      w(15)=16*LL
      w(20)=8.D0/3.D0*w(3)
      w(21)=40.D0/3.D0*w(2) - w(20) - 343.D0/9.D0 + w(15)
      w(21)=z3*w(21)
      w(22)=8*B4
      w(23)= - 6197.D0/81.D0 + 56*w(9)
      w(21)=w(21) + 8*w(3) + 1.D0/3.D0*w(23) - w(22)
      w(12)= - 7.D0/3.D0*w(3) + 31.D0/9.D0 + w(12)
      w(23)=2.D0/9.D0*w(2)
      w(12)=w(12)*w(23)
      w(24)=31.D0/3.D0 + w(11)
      w(24)=LL*w(24)
      w(24)=400.D0/27.D0 + w(24)
      w(25)=LL + 10.D0/9.D0
      w(26)= - 2*w(25) + 1.D0/9.D0*w(3)
      w(26)=w(3)*w(26)
      w(12)=w(12) + 1.D0/3.D0*w(24) + w(26)
      w(12)=w(2)*w(12)
      w(24)=w(3)**2
      w(26)= - 1 - 2.D0/3.D0*w(24)
      w(26)=w(26)*w(19)
      w(27)=407.D0/3.D0 + 224*LL
      w(26)=1.D0/3.D0*w(27) + w(26)
      w(12)=1.D0/9.D0*w(26) + w(12)
      w(12)=w(12)*w(14)
      w(26)=w(16) + 20.D0/3.D0*w(3) - 67.D0/9.D0 - w(17)
      w(27)=4.D0/3.D0*w(5)
      w(26)=w(26)*w(27)
      w(12)= - 56.D0/9.D0*w(8) + w(26) + w(12) + 1.D0/3.D0*w(21) + 2*
     & w(13)
      w(12)=CA*w(12)
      w(13)= - w(16) + w(19) - w(18)
      w(13)=w(13)*w(14)
      w(16)= - 59 + 16*w(24)
      w(13)= - 136.D0/15.D0*z2 - 16.D0/3.D0*w(4) + 1.D0/3.D0*w(16) + 
     & w(13)
      w(13)=z2*w(13)
      w(15)= - 16.D0/3.D0*w(2) + 32.D0/3.D0*w(3) + 293.D0/9.D0 - w(15)
      w(15)=z3*w(15)
      w(13)=w(15) + w(13) - 769.D0/9.D0 + w(22)
      w(15)=4*w(25) - 1.D0/3.D0*w(3)
      w(15)=w(15)*w(20)
      w(16)=w(20) + 1.D0/9.D0 + LL
      w(16)=4*w(16) + w(2)
      w(16)=w(16)*w(23)
      w(15)=w(16) + w(15) - 1 - w(17)
      w(15)=w(2)*w(15)
      w(16)=10.D0/3.D0 + LL
      w(16)=LL*w(16)
      w(16)=2.D0/9.D0*w(24) + 112.D0/27.D0 + w(16)
      w(16)=w(16)*w(20)
      w(11)=w(11) + 5.D0/3.D0
      w(18)= - LL*w(11)
      w(16)=w(16) - 302.D0/27.D0 + w(18)
      w(15)=2*w(16) + w(15)
      w(15)=w(2)*w(15)
      w(16)= - w(2) + w(20) - 49.D0/9.D0 - w(17)
      w(14)=w(16)*w(14)
      w(16)=7 - 8*w(24)
      w(14)=2.D0/3.D0*w(4) + 1.D0/3.D0*w(16) + w(14)
      w(14)=w(4)*w(14)
      w(16)=8.D0/3.D0*w(2) - 28.D0/3.D0*w(3) + 107.D0/9.D0 + 8*LL
      w(16)=w(16)*w(27)
      w(17)=w(2) - w(3)
      w(17)= - 1 - 4.D0/3.D0*w(17)
      w(17)=w(6)*w(17)
      w(13)=16.D0/3.D0*w(8) + 4*w(17) + w(16) + 4.D0/3.D0*w(14) - 32.D0/
     & 9.D0*w(7) + w(15) + 2.D0/3.D0*w(13)
      w(13)=CF*w(13)
      w(14)=LL**2
      w(14)=1.D0/9.D0 + w(14)
      w(15)= - 5 - w(2)
      w(15)=w(2)*w(15)
      w(14)=2*w(14) + 1.D0/9.D0*w(15)
      w(14)=w(2)*w(14)
      w(15)= - w(23) - w(25)
      w(15)=w(2)*w(15)
      w(16)=2.D0/3.D0 - 5*LL
      w(15)=2.D0/3.D0*w(16) + w(15)
      w(15)=w(2)*w(15)
      w(15)= - 14.D0/3.D0*z3 + 5.D0/3.D0*z2 + 1736.D0/81.D0 + w(15)
      w(15)=NF*w(15)
      w(14)=w(15) + 32.D0/3.D0*z3 + 10.D0/3.D0*z2 + 607.D0/81.D0 + 
     & w(14)
      w(14)=w(14)*TF
      w(12)=8.D0/9.D0*w(14) + w(12) + w(13)
      w(12)=as*CF*w(12)
      w(11)=2*w(11) + w(2)
      w(10)=CF*w(11)*w(10)
      w(10)=w(10) + w(12)
      ANSPLU1 = 4.0D0*TF*as**2*w(10)*w(1)
*
      RETURN
      END      
      REAL*8 FUNCTION ANSDEL(z,nf,as,LL)
*     ----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Delta[1-z] part of aqqQNS3; This is a constant.
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(9)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*      LL=10.0D0
*
      w(1)= - 4163.D0/3.D0 - 1792*LL
      w(2)=4*LL
      w(3)=359.D0/45.D0 + w(2)
      w(3)=z2*w(3)
      w(1)=1.D0/9.D0*w(1) + 8*w(3)
      w(1)=z2*w(1)
      w(3)=16*LL
      w(4)=2*LL
      w(5)=31.D0/3.D0 + w(4)
      w(5)=w(5)*w(3)
      w(6)=88.D0/3.D0*z2
      w(5)= - w(6) - 7801.D0/27.D0 + w(5)
      w(5)=z3*w(5)
      w(1)=w(5) - 7663.D0/162.D0 + w(1)
      w(1)= - 176.D0/9.D0*z5 + 1.D0/3.D0*w(1) - 8*B4
      w(1)=CA*w(1)
      w(5)=w(4) + 5.D0/3.D0
      w(7)=8*LL
      w(8)=w(5)*w(7)
      w(9)=8.D0/5.D0*z2
      w(7)= - 961.D0/27.D0 + w(7)
      w(7)=w(7)*w(9)
      w(7)=w(7) - 119.D0/3.D0 + w(8)
      w(7)=z2*w(7)
      w(4)= - 29.D0/3.D0 - w(4)
      w(3)=w(4)*w(3)
      w(3)=w(6) + 8137.D0/27.D0 + w(3)
      w(3)=z3*w(3)
      w(3)=352.D0/9.D0*z5 + 16*B4 + 2.D0/3.D0*w(3) - 2477.D0/18.D0 + 
     & w(7)
      w(3)=CF*w(3)
      w(4)=2*z2
      w(6)=LL**2
      w(6)=56.D0/5.D0*z2 + 233.D0/9.D0 - 8*w(6)
      w(6)=w(6)*w(4)
      w(6)=112.D0/3.D0*z3 - 139.D0/27.D0 + w(6)
      w(7)=457.D0/3.D0 + 80*LL
      w(7)=1.D0/27.D0*w(7) + w(9)
      w(4)=w(7)*w(4)
      w(2)= - 41.D0/3.D0 - w(2)
      w(2)=z3*w(2)
      w(2)=8.D0/9.D0*w(2) + 3917.D0/243.D0 + w(4)
      w(2)=NF*w(2)
      w(2)=2.D0/9.D0*w(6) + w(2)
      w(2)=TF*w(2)
      w(1)=2*w(2) + w(3) + w(1)
      w(1)=as*CF*w(1)
      w(2)= - z2*w(5)
      w(2)=w(2) + z3
      w(2)=CF*w(2)
      w(1)=8.D0/3.D0*w(2) + w(1)
*
      ANSDEL= TF*as**2*w(1)
*      WRITE(6,*) 'LL,ANSDEL=',LL,ANSDEL
*
      RETURN
      END      
      REAL*8 FUNCTION NSREG(z,nf,as,LL)
*     ---------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of OME ANS without aqqQNS3
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(44)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
      w(1)=Hr1(0)
      w(2)=Hr2(0,1)
      w(3)=Hr1(1)
      w(4)=1.0D0/( - 1.0D0 + z)
      w(5)=Hr3(0,0,1)
      w(6)=1.0D0/(1.0D0 + z)
      w(7)=Hr1(-1)
      w(8)=Hr3(-1,0,1)
      w(9)=Hr2(0,-1)
      w(10)=Hr3(0,0,-1)
      w(11)=1.0D0/(3.0D0 + 3.0D0*z)
      w(12)=z**2
      w(12)=w(12) + 1
      w(13)=w(1)*w(6)
      w(14)=w(12)*w(13)
      w(15)=1 + 4.D0/3.D0*z
      w(15)=w(15)*z
      w(15)=w(15) + 4.D0/3.D0
      w(15)=w(15)*w(6)
      w(16)=w(14) + 8*w(15)
      w(16)=w(16)*w(1)
      w(17)=w(12)*w(6)
      w(18)=LL*w(17)
      w(19)=2*w(1)
      w(20)=w(18)*w(19)
      w(16)=w(16) + w(20)
      w(16)=w(16)*LL
      w(14)= - w(14) + 4*w(18)
      w(14)=w(14)*z2
      w(14)=w(16) + w(14)
      w(14)=w(14)*w(7)
      w(16)=4*w(4)
      w(20)= - 10 - z
      w(20)=z*w(20)
      w(20)= - 1 + w(20)
      w(20)=w(6)*w(20)
      w(20)= - w(16) + w(20)
      w(20)=w(1)*w(20)
      w(21)=5*z
      w(22)=2 + w(21)
      w(22)=z*w(22)
      w(22)=17 + w(22)
      w(22)=w(6)*w(22)
      w(22)=w(4) + 1.D0/9.D0*w(22)
      w(20)=2*w(22) + 1.D0/9.D0*w(20)
      w(20)=w(1)*w(20)
      w(20)=w(20) + 5.D0/3.D0*w(4) + 37.D0/9.D0 + w(21)
      w(20)=w(20)*w(19)
      w(22)=z - 2
      w(22)=w(22)*z
      w(22)=w(22) + 1
      w(13)=1.D0/3.D0*w(13)
      w(13)=w(22)*w(13)
      w(22)=z + 1
      w(23)=2.D0/3.D0*w(22)
      w(24)=w(6)*w(23)
      w(24)=w(4) + w(24)
      w(24)=2*w(24) + w(13)
      w(24)=w(24)*w(19)
      w(25)= - 2 + 1.D0/3.D0*z
      w(25)=z*w(25)
      w(25)= - 7.D0/3.D0 + w(25)
      w(25)=w(6)*w(25)
      w(24)=w(25) + w(24)
      w(24)=LL*w(24)
      w(25)=19*z
      w(26)=13 - w(25)
      w(20)=w(24) + 5.D0/9.D0*w(26) + w(20)
      w(24)=4*LL
      w(20)=w(20)*w(24)
      w(26)= - 2.D0/3.D0 - w(4)
      w(13)=2*w(26) - w(13)
      w(13)=w(1)*w(13)
      w(26)=w(22)*z
      w(26)=w(26) + 1
      w(27)=w(26)*w(6)
      w(27)=w(27) + w(4)
      w(28)=w(27)*w(1)
      w(29)= - 5.D0/3.D0*w(17) + w(28)
      w(29)=2*w(29) - w(18)
      w(29)=LL*w(29)
      w(30)= - 13 - 17*z
      w(13)=8.D0/3.D0*w(29) + 1.D0/3.D0*w(30) + w(13)
      w(29)=4*z2
      w(13)=w(13)*w(29)
      w(30)=z + 2
      w(30)=w(30)*z
      w(30)=w(30) + 1
      w(30)=w(30)*w(6)
      w(31)=2*w(4)
      w(30)=w(30) + w(31)
      w(32)=w(30)*w(1)
      w(33)=4*z
      w(34)= - 5*w(4) - 1 - w(33)
      w(34)=2.D0/3.D0*w(34) - w(32)
      w(34)=w(1)*w(34)
      w(35)=z - 1
      w(34)=2*w(35) + w(34)
      w(36)= - LL*w(32)
      w(34)=2*w(34) + w(36)
      w(36)=2*LL
      w(34)=w(34)*w(36)
      w(37)=w(31) + w(22)
      w(38)=w(37)*w(1)
      w(39)=z2*w(38)
      w(34)=w(34) + w(39)
      w(34)=w(3)*w(34)
      w(33)=5 + w(33)
      w(33)=z*w(33)
      w(33)=4 + w(33)
      w(33)=w(6)*w(33)
      w(33)=w(31) + w(33)
      w(33)=z3*LL*w(33)
      w(39)=w(18)*w(10)
      w(33)=w(39) + w(33)
      w(31)= - w(31) + 5 - 7*z
      w(31)=2.D0/3.D0*w(31) - w(32)
      w(31)=w(1)*w(31)
      w(32)=7 - 473*z
      w(12)=w(11)*w(12)*z2**2
      w(12)=16.D0/3.D0*w(34) - 32.D0/3.D0*w(14) + 16*w(12) + w(13) + 
     & w(20) + 1.D0/9.D0*w(32) + 8*w(31) + 64.D0/3.D0*w(33)
      w(13)=CF**2
      w(12)=w(12)*w(13)
      w(20)= - 1 + 11*z
      w(31)=w(20) + 10*w(4)
      w(32)=w(31) + w(38)
      w(32)=w(32)*w(1)
      w(33)= - 11 + 67*z
      w(34)=56*w(4)
      w(40)=w(33) + w(34)
      w(32)=w(32) + 2.D0/3.D0*w(40)
      w(32)=w(32)*w(1)
      w(40)=463 - 1880*z
      w(40)=4.D0/9.D0*w(40) - 11*w(32)
      w(41)=w(1)**2
      w(42)=w(27)*w(41)
      w(43)=w(22)*LL
      w(44)= - 11 - 35*z
      w(42)=22.D0/9.D0*w(43) + 1.D0/3.D0*w(44) - 4*w(42)
      w(42)=LL*w(42)
      w(25)= - 68.D0/3.D0 - w(25)
      w(25)=z*w(25)
      w(25)= - 17 + w(25)
      w(25)=w(6)*w(25)
      w(25)= - 8.D0/3.D0*w(28) - 62.D0/3.D0*w(4) + w(25)
      w(25)=w(1)*w(25)
      w(28)= - w(34) - 7 - 76*z
      w(25)=8.D0/9.D0*w(28) + w(25)
      w(25)=w(1)*w(25)
      w(28)=2071 - 2381*z
      w(25)=w(42) + 1.D0/27.D0*w(28) + w(25)
      w(25)=LL*w(25)
      w(19)=w(27)*w(19)
      w(19)= - 11.D0/3.D0*w(37) + w(19)
      w(19)=w(1)*w(19)
      w(28)=4 + z
      w(28)=z*w(28)
      w(28)=1 + w(28)
      w(28)=w(6)*w(28)
      w(16)=w(16) + w(28)
      w(16)=w(1)*w(16)
      w(28)=23*z
      w(34)=14 + w(28)
      w(34)=z*w(34)
      w(34)=11 + w(34)
      w(34)=w(6)*w(34)
      w(16)=1.D0/3.D0*w(34) + w(16)
      w(16)=w(16)*w(24)
      w(28)=25 - w(28)
      w(16)=w(16) + 2.D0/9.D0*w(28) + w(19)
      w(16)=z2*w(16)
      w(16)=w(16) + 1.D0/9.D0*w(40) + w(25)
      w(19)=LL**2
      w(25)=2*w(19) - z2
      w(25)=w(11)*w(29)*w(26)*w(25)
      w(26)= - 3 - 2*z
      w(26)=z*w(26)
      w(26)= - 2 + w(26)
      w(26)=w(6)*w(26)
      w(26)= - 2.D0/3.D0*w(4) + w(26)
      w(24)=w(26)*w(24)
      w(24)= - 11.D0/27.D0*w(22) + w(24)
      w(24)=z3*w(24)
      w(26)=w(30)*w(41)
      w(26)= - 4*w(35) + w(26)
      w(26)=w(3)*w(26)*w(36)
      w(14)= - 8.D0/3.D0*w(39) + w(26) + 4.D0/3.D0*w(14) + 2*w(24) + 1.D
     & 0/3.D0*w(16) + w(25)
      w(16)=CA*CF
      w(14)=w(14)*w(16)
      w(24)= - 49 + 377*z
      w(24)=w(32) + 2.D0/9.D0*w(24)
      w(24)=1.D0/3.D0*w(24)
      w(25)=w(38) + 2.D0/3.D0*w(31)
      w(25)=w(25)*w(1)
      w(19)= - w(19)*w(23)
      w(21)= - 73 + w(21)
      w(19)=w(19) + 2.D0/9.D0*w(21) + w(25)
      w(19)=LL*w(19)
      w(20)=w(38) + 1.D0/3.D0*w(20)
      w(21)=w(20)*z2
      w(26)=z3*w(23)
      w(19)=w(26) + w(21) + w(24) + w(19)
      w(19)=nf*w(19)
      w(23)= - LL*w(23)
      w(23)=w(23) - w(20)
      w(23)=LL*w(23)
      w(23)= - 31.D0/9.D0*w(22) + w(23)
      w(23)=w(23)*w(36)
      w(26)=z3*w(22)
      w(19)=w(19) + 4.D0/3.D0*w(26) + 2*w(21) + w(24) + w(23)
      w(21)=TF*CF
      w(19)=w(19)*w(21)
      w(23)= - w(16) + 2*w(13)
      w(24)=w(8)*w(18)*w(23)
      w(26)= - w(13)*w(36)
      w(28)=LL*w(16)
      w(26)=w(26) + w(28)
      w(26)=w(5)*w(27)*w(26)
      w(24)=w(24) + w(26)
      w(15)=w(18) + 4*w(15)
      w(15)=w(15)*w(36)
      w(17)=w(17)*z2
      w(15)=w(15) - w(17)
      w(15)=w(9)*w(15)*w(23)
      w(13)=4.D0/3.D0*w(13) - w(16)
      w(16)=w(38) - w(22)
      w(13)=w(2)*LL*w(16)*w(13)
      w(12)=16*w(13) + 16.D0/9.D0*w(19) + 16.D0/3.D0*w(15) + w(12) + 4*
     & w(14) + 64.D0/3.D0*w(24)
      w(12)=as*TF*w(12)
      w(13)= - 2*w(20) - w(43)
      w(13)=w(13)*w(36)
      w(13)=w(13) - 2.D0/9.D0*w(33) - w(25)
      w(13)=w(13)*w(21)
      w(12)=2.D0/3.D0*w(13) + w(12)
*
      NSREG = as**2*w(12)
*
      RETURN
      END      
      REAL*8 FUNCTION NSPLU(z,nf,as,LL)
*     ---------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     plus part of OME ANS without aqqQNS3
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(7),G
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5 
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5 
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5 
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half,x
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      w(1)=1.0D0/(243 - 243*z)
      w(2)=3*LL
      w(3)=69 - 22*LL
      w(3)=w(3)*w(2)
      w(3)=155 + w(3)
      w(3)=w(3)*w(2)
      w(4)=18*z3
      w(5)=11 + 126*LL
      w(5)=w(5)*w(4)
      w(6)=w(2) + 10
      w(7)=LL*w(6)
      w(7)=9*z2 - 1 - 6*w(7)
      w(7)=z2*w(7)
      w(3)=w(5) + 18*w(7) + 2834 + w(3)
      w(3)=CA*w(3)
      w(5)=5 + 2*LL
      w(5)=w(5)*w(2)
      w(5)=31 + w(5)
      w(5)=w(5)*w(2)
      w(4)=w(4) + 45*z2
      w(5)= - 164 + w(5) - w(4)
      w(7)=LL**2
      w(7)=34 + 3*w(7)
      w(7)=w(7)*w(2)
      w(7)= - 164 + w(7)
      w(4)=2*w(7) - w(4)
      w(4)=nf*w(4)
      w(4)=2*w(5) + w(4)
      w(4)=TF*w(4)
      w(5)=5 + w(2)
      w(5)= - 288*z3 + 12*w(5)
      w(5)=LL*w(5)
      w(5)=180*z2 + 233 + w(5)
      w(5)=CF*w(5)
      w(3)=16*w(4) + 4*w(3) + 27*w(5)
      w(3)=as*CF*w(3)
      w(2)=w(6)*w(2)
      w(2)=28 + w(2)
      w(2)=CF*w(2)
      w(2)=36*w(2) + w(3)

      NSPLU = 2.0D0*TF*as**2*w(2)*w(1)
*
      RETURN
      END      
      REAL*8 FUNCTION NSDEL(z,nf,as,LL)
*     ---------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Delta[1-z] part of OME ANS without aqqQNS3
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(9)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      w(1)=8*z2
      w(2)=224.D0/27.D0 - 1.D0/5.D0*z2
      w(2)=w(2)*w(1)
      w(3)=2*LL
      w(4)=8*z3
      w(5)= - 22.D0/9.D0*LL + 17.D0/3.D0 - w(4)
      w(5)=w(5)*w(3)
      w(2)=w(5) - 224.D0/9.D0*z3 - 1595.D0/27.D0 + w(2)
      w(2)=LL*w(2)
      w(5)=4*z3
      w(6)=2*z2
      w(7)= - 187.D0/27.D0 + w(6)
      w(7)=w(7)*w(5)
      w(8)=4163.D0/81.D0 + 88.D0/5.D0*z2
      w(8)=z2*w(8)
      w(2)=w(2) + w(7) + 2591.D0/162.D0 + w(8)
      w(2)=CA*w(2)
      w(7)=16*z3
      w(8)= - 2 - z2
      w(8)=w(8)*w(7)
      w(5)=w(5) + 1 - w(6)
      w(5)=LL*w(5)
      w(9)= - 5 - 58.D0/5.D0*z2
      w(9)=z2*w(9)
      w(5)=8*w(5) + 272.D0/3.D0*z3 + 1 + 8.D0/3.D0*w(9)
      w(5)=LL*w(5)
      w(9)=119.D0/3.D0 + w(1)
      w(9)=z2*w(9)
      w(5)=w(5) + w(8) + 365.D0/6.D0 + w(9)
      w(5)=CF*w(5)
      w(2)=w(2) + w(5)
      w(2)=CF*w(2)
      w(5)=35 - 16*z2
      w(8)=LL**2
      w(5)=4*w(8) + 5.D0/3.D0*w(5) + w(7)
      w(5)=LL*w(5)
      w(7)= - 457.D0/81.D0 - 8.D0/5.D0*z2
      w(6)=w(7)*w(6)
      w(5)=2.D0/9.D0*w(5) + 136.D0/27.D0*z3 - 517.D0/81.D0 + w(6)
      w(5)=nf*w(5)
      w(1)=w(1) + 1
      w(6)=4*LL
      w(7)=w(6) + w(1)
      w(7)=LL*w(7)
      w(7)=62.D0/3.D0 + w(7)
      w(6)=w(7)*w(6)
      w(7)= - 233.D0/9.D0 - 56.D0/5.D0*z2
      w(7)=z2*w(7)
      w(6)=w(6) + 112.D0/3.D0*z3 - 517.D0/9.D0 + 4*w(7)
      w(5)=1.D0/9.D0*w(6) + w(5)
      w(6)=TF*CF
      w(5)=w(5)*w(6)
      w(2)=w(2) + 2*w(5)
      w(2)=as*TF*w(2)
      w(5)=73.D0/2.D0 + 40*z2
      w(4)=1.D0/3.D0*w(5) - w(4)
      w(1)=1.D0/3.D0*w(1) + LL
      w(1)=w(1)*w(3)
      w(1)=1.D0/3.D0*w(4) + w(1)
      w(1)=w(1)*w(6)
      w(1)=w(1) + w(2)
      w(1)=w(1)*as**2
*
      NSDEL = 1.0D0 + w(1)
*
      RETURN
      END      
      REAL*8 FUNCTION AGGREG0(z,nf,as,LLQ)
*     ------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of aggQ3 in the range: representation for  z \in [0,1/2]
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(13),x,LLQ
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half,G,Lx,z4
      REAL*8 B4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
      x=z
      LL=Log(x)
*      Lx = LOG(z)
      G =
     &  + LL * (  - 4889048.D0/729.D0 - 185168.D0/243.D0*x**(-1) - 
     &    2116042543.D0/182250.D0*x - 13975819327.D0/5358150.D0*x**2 + 
     &    251062814.D0/637875.D0*x**3 + 27774049549.D0/10692000.D0*x**4
     &     + 2714121599828047.D0/2145403260000.D0*x**5 + 
     &    3044357932851601.D0/536350815000.D0*x**6 + 160067472875227.D0/
     &    78153975900.D0*x**7 + 44472696349419288454121.D0/
     &    4610159235025344000.D0*x**8 + 21667742908178434974299.D0/
     &    7244535940754112000.D0*x**9 + 287626693061430177489329.D0/
     &    19620618172875720000.D0*x**10 + 29634937290080233922298149.D0/
     &    7170245907843138120000.D0*x**11 + 535881993438848032482998981.
     &    D0/25812885268235297232000.D0*x**12 + 
     &    4855649803732904087214420991.D0/884677976920427914224000.D0*
     &    x**13 + 2852616615373467379859039218871.D0/
     &    101986783026858080486635500.D0*x**14 + 75370910879296031533015
     &    3134237529.D0/106693865320405376509095600000.D0*x**15 + 
     &    96118352570202312074054308145582011.D0/
     &    2647749800604345670103270400000.D0*x**16 )
      G = G + LL * ( 1177972822683428433567098496836199812791.D0/
     &    132901683041494647434297575111680000.D0*x**17 + 22806941803386
     &    0902256568324047877651327.D0/498381311405604927878615906668800
     &    0.D0*x**18 + 4144477967700582437630577263361322826377.D0/
     &    380626511357810038056315481661760000.D0*x**19 + 53036452030011
     &    488579026575922508693886403247.D0/9410990493321853190942400284
     &    08701600000.D0*x**20 + 998024490994404586563619301151641517997
     &    8431.D0/759483443320710959269035812400004800000.D0*x**21 + 
     &    20793540086494180039542812450086374463177481803.D0/30539673126
     &    4183662399404511675068596800000.D0*x**22 + 7053769527289522397
     &    39870819538992057165272203413.D0/45155088122632870083340524226
     &    242285384000000.D0*x**23 + 13293331972512262417068870302896137
     &    02323944880313.D0/16420032044593770939396554264088103776000000
     &    .D0*x**24 + 105513192986287572723025013786171922446823170243.D
     &    0/5756206665249161256204365746516466400000000.D0*x**25 + 
     &    29818901564392951462867419789923301619486888918187311.D0/
     &    313978987274637040505041193416945847417400000000.D0*x**26 )
      G = G + LL * ( 220088487695527447962430743620664183174414893896513
     &    069.D0/10348008646484331826715557638096848199471816000000.D0*
     &    x**27 + 635987886231206342175866176243090034890034260626367359
     &    229.D0/5775011358758875144657609834588248993495487754400000.D0
     &    *x**28 + 25796551305938089763610982383154847756514688612579046
     &    487463607.D0/1055652826343259846859595552396616621847663509878
     &    472000000.D0*x**29 + 39721743244618836060443245688872523366702
     &    81633705782093953899.D0/31418238879263685918440343821327875650
     &    228080651145000000.D0*x**30 + 21367215629696388504705242548948
     &    43930819174726891994557907941.D0/76765746227166444411558672588
     &    072050849818364546640000000.D0*x**31 + 61109598996620551856322
     &    8683330829304864149435359727540843967883.D0/424737444931728662
     &    7571143071711651019922853511948288000000.D0*x**32 + 
     &    14699421183442334755865731653490795561617273850687546509007134
     &    19.D0/46721118942490152903282573788828161219151388631431168000
     &    000.D0*x**33 )
      G = G + LL * ( 741555884081506471076364775286220050261906247128065
     &    88018053214994397.D0/45642591112402106779855860124521681135008
     &    0407944924741312000000.D0*x**34 + 1006379078731130530309400715
     &    581267858351668049918113089532096778695801.D0/2849321385083160
     &    7953557167533821200480506175235575694752320000000.D0*x**35 + 
     &    54055736471962830042063447760794093362481820805395987126551059
     &    677462583.D0/2966646383292467416340952149109619108852701774527
     &    58704185920000000.D0*x**36 + 162619257770679423063037146648689
     &    427611300485144753405746164657356019.D0/4126538202324108812203
     &    580057032176955922931039756418065744000000.D0*x**37 + 
     &    26132745579824181565204945594242545204288693128714680323987252
     &    6339105558699.D0/128668652338008843408655730280174049035856236
     &    8862758602532020953600000.D0*x**38 + 1350809470174788399310467
     &    267866982086624416162651403176672620292581429374519.D0/
     &    30892068331603024115681758666366111773023141558731997078809331
     &    904000000.D0*x**39 )
      G = G + LL * ( 130644487051146582509932654864001805752092646581690
     &    69852010980851028405534709381.D0/58028311513416417404667471936
     &    979291009431364327954998723300266087040000000.D0*x**40 + 
     &    16831361802294251970403330549223053882672221342167923028410516
     &    886944372593375553.D0/3486501423500837425467386909319675902372
     &    89419265478402747930177108160000000.D0*x**41 + 
     &    95756407134086764104129052082336828941011707430250775782247401
     &    58201036365511107141.D0/38560705743919261925669299217075615480
     &    244209770761911343921077588162496000000.D0*x**42 + 
     &    17360261346808802228897928431667694776886489764478335365221146
     &    016396382326776738871.D0/3272102458580465973160452507843644799
     &    17726365394026241049192824678265792000000.D0*x**43 + 
     &    21745007160199503915596715838019818964261094762994560838741125
     &    6498573615503791739796311.D0/797504857797733285258435431947434
     &    627410904217141068814077182683119413239616000000.D0*x**44 + 
     &    54230596107342065301649174229622539919154461847382911630746610
     &    822846300866731561178001.D0/9339550244973454586498288696394375
     &    62001141981867862979940886697008615172640000000.D0*x**45 )
      G = G + LL * ( 247757737322086646363185238473088846731883901104351
     &    33078800741009060828066096727607389.D0/83098512817945630750274
     &    909097512819829888841713581619104991853118754928320000000.D0*
     &    x**46 + 548419329646584377795714537369543004344975614758698089
     &    473873603152210483648616104832587161.D0/8662947942559723119502
     &    84236849449362282205917831139868762550969469083168975545600000
     &    0.D0*x**47 + 1533210164137149303093268105813011664909375388324
     &    4802327117525697453734333739183566903858033.D0/
     &    47206788788441389752653170008027965249001365957175158065901328
     &    191358735005044224000000.D0*x**48 + 48839819020964411870151193
     &    96493564744158400387755290741257914515789189337295869030905146
     &    7527.D0/710110631349533245853740238205867477256254590036656101
     &    118983809176396290395026944000000.D0*x**49 + 
     &    29897579965357668669500734135889870714615521934301379939977600
     &    834339718206097067469601131917050507.D0/
     &    84798304861749106011903049133002543463728927028346061225186650
     &    937680098265203639409600000000.D0*x**50 )
      G = G + LL * (  - 16250.D0/27.D0*z4 + 46942.D0/27.D0*z4*x - 1808.D
     &    0/27.D0*z2*x**(-1) - 67112.D0/81.D0*z2 - 2366.D0/5.D0*z2*x - 
     &    3410.D0/27.D0*z2*x**2 - 182696.D0/1215.D0*z2*x**3 + 40469.D0/
     &    1215.D0*z2*x**4 - 15615848.D0/42525.D0*z2*x**5 + 7550812.D0/
     &    42525.D0*z2*x**6 - 27276008.D0/59535.D0*z2*x**7 + 40141589.D0/
     &    119070.D0*z2*x**8 - 870429953.D0/1683990.D0*z2*x**9 + 
     &    38700071.D0/76545.D0*z2*x**10 - 112476338638.D0/200675475.D0*
     &    z2*x**11 + 9087314903.D0/13378365.D0*z2*x**12 - 9411184282.D0/
     &    15810795.D0*z2*x**13 + 31589882906.D0/36891855.D0*z2*x**14 - 
     &    41060793164.D0/65786175.D0*z2*x**15 + 428252997547.D0/
     &    413513100.D0*z2*x**16 - 9455415443329.D0/14570697960.D0*z2*
     &    x**17 + 146259523692257.D0/120208258170.D0*z2*x**18 - 
     &    4095354836749.D0/6106836645.D0*z2*x**19 + 10442723839456.D0/
     &    7463911455.D0*z2*x**20 - 218170126683737.D0/316234143225.D0*
     &    z2*x**21 + 1100970034722707.D0/695715115095.D0*z2*x**22 - 
     &    8083741879500787.D0/11429605462275.D0*z2*x**23 )
      G = G + LL * ( 3671475506298011.D0/2078110084050.D0*z2*x**24 - 
     &    148477294534757.D0/205346846250.D0*z2*x**25 + 
     &    85964523382553386.D0/44046898520625.D0*z2*x**26 - 
     &    92497552418660531.D0/125413532878725.D0*z2*x**27 + 
     &    317510160149090756.D0/148566800487105.D0*z2*x**28 - 
     &    5065338666200927071.D0/6745533012015525.D0*z2*x**29 + 
     &    3530544037967039.D0/1519730939727.D0*z2*x**30 - 
     &    961065647012005342.D0/1259018476791075.D0*z2*x**31 + 
     &    24780851520702691.D0/9874654719930.D0*z2*x**32 - 
     &    11078643816049283191.D0/14295951736466400.D0*z2*x**33 + 
     &    16382184601557398147.D0/6075779487998220.D0*z2*x**34 - 
     &    843142621711194450197.D0/1072927423221503850.D0*z2*x**35 + 
     &    49631203570403052763.D0/17212739409970650.D0*z2*x**36 - 
     &    1957535605396151563.D0/2458962772852950.D0*z2*x**37 + 
     &    25135751564280701.D0/8185481115210.D0*z2*x**38 - 
     &    2229452739491770961309.D0/2766879555631313850.D0*z2*x**39 )
      G = G + LL * ( 47450507151076368824.D0/14562523977006915.D0*z2*
     &    x**40 - 1149591262813325402662.D0/1410644492937537975.D0*z2*
     &    x**41 + 1048125394132600833486751.D0/304134952677333187410.D0
     &    *z2*x**42 - 1313618674252864443454531.D0/
     &    1594854020137235007150.D0*z2*x**43 + 
     &    1821659827415071708948153.D0/501239834900273859390.D0*z2*
     &    x**44 - 2930193566624079959042221.D0/3522000833269366154850.D0
     &    *z2*x**45 + 28149957078528902217594769.D0/
     &    7364183560472311051050.D0*z2*x**46 - 
     &    431650885109012970462459467.D0/513930749689931283350550.D0*z2
     &    *x**47 + 28575819063480613186689907.D0/7124433833130742428300.
     &    D0*z2*x**48 - 26323799803975131751442089.D0/
     &    31060877872101921534600.D0*z2*x**49 + 
     &    298927097606555805513802804.D0/71181178456900236850125.D0*z2*
     &    x**50 + 1280.D0/9.D0*z3*x**(-1) - 9760.D0/27.D0*z3 + 92192.D0/
     &    27.D0*z3*x )
      G = G + LL * ( 3584.D0/81.D0*z3*x**2 - 72*z3*x**3 - 1322*z3*x**4
     &     - 1998.D0/5.D0*z3*x**5 - 3464*z3*x**6 - 6480.D0/7.D0*z3*x**7
     &     - 46049.D0/7.D0*z3*x**8 - 1645*z3*x**9 - 10656*z3*x**10 - 
     &    140616.D0/55.D0*z3*x**11 - 172642.D0/11.D0*z3*x**12 - 47586.D0
     &    /13.D0*z3*x**13 - 1974152.D0/91.D0*z3*x**14 - 173472.D0/35.D0
     &    *z3*x**15 - 57307.D0/2.D0*z3*x**16 - 219105.D0/34.D0*z3*x**17
     &     - 621744.D0/17.D0*z3*x**18 - 154360.D0/19.D0*z3*x**19 - 
     &    863606.D0/19.D0*z3*x**20 - 349866.D0/35.D0*z3*x**21 - 4257544.
     &    D0/77.D0*z3*x**22 - 3051216.D0/253.D0*z3*x**23 - 1520131.D0/
     &    23.D0*z3*x**24 - 357903.D0/25.D0*z3*x**25 - 5060416.D0/65.D0*
     &    z3*x**26 - 653800.D0/39.D0*z3*x**27 - 1902022.D0/21.D0*z3*
     &    x**28 - 3939030.D0/203.D0*z3*x**29 - 3023320.D0/29.D0*z3*
     &    x**30 - 3446592.D0/155.D0*z3*x**31 - 14742653.D0/124.D0*z3*
     &    x**32 - 1111443.D0/44.D0*z3*x**33 - 25150064.D0/187.D0*z3*
     &    x**34 - 16943256.D0/595.D0*z3*x**35 - 1057366.D0/7.D0*z3*
     &    x**36 )
      G = G + LL * (  - 1179710.D0/37.D0*z3*x**37 - 118506296.D0/703.D0
     &    *z3*x**38 - 8764560.D0/247.D0*z3*x**39 - 2431679.D0/13.D0*z3*
     &    x**40 - 8051589.D0/205.D0*z3*x**41 - 59263264.D0/287.D0*z3*
     &    x**42 - 13021272.D0/301.D0*z3*x**43 - 107320006.D0/473.D0*z3*
     &    x**44 - 2608982.D0/55.D0*z3*x**45 - 5709800.D0/23.D0*z3*x**46
     &     - 56000160.D0/1081.D0*z3*x**47 - 25433783.D0/94.D0*z3*x**48
     &     - 5523675.D0/98.D0*z3*x**49 - 71993776.D0/245.D0*z3*x**50 )
      G = G + LL**2 * (  - 343736.D0/243.D0 - 9302411.D0/12150.D0*x + 
     &    12019933.D0/51030.D0*x**2 - 986516.D0/18225.D0*x**3 - 
     &    54974833.D0/145800.D0*x**4 - 5588464697.D0/17860500.D0*x**5
     &     - 1744385252.D0/1488375.D0*x**6 - 2414816837.D0/5000940.D0*
     &    x**7 - 90720296601889.D0/39607444800.D0*x**8 - 14706151995529.
     &    D0/23340101400.D0*x**9 - 11298372652645051.D0/3034213182000.D0
     &    *x**10 - 27628910417725271.D0/36157707085500.D0*x**11 - 
     &    3196696327171849.D0/584366983200.D0*x**12 - 20262535172954341.
     &    D0/22790312344800.D0*x**13 - 1276465394193926987.D0/
     &    169502948064450.D0*x**14 - 4465376465797767287.D0/
     &    4433154026301000.D0*x**15 - 101702370616426925831.D0/
     &    10268029135203840.D0*x**16 - 135694802985565738183.D0/
     &    121141074253399200.D0*x**17 - 78310457456991162298757.D0/
     &    6218575145007825600.D0*x**18 - 6404789003915655695999.D0/
     &    5212629165668324400.D0*x**19 - 74790357584699352428719907.D0/
     &    4795618832414858448000.D0*x**20 )
      G = G + LL**2 * (  - 2007299361434569889126039.D0/
     &    1505057762219867856000.D0*x**21 - 21348484576503777292088327.D
     &    0/1128793321664900892000.D0*x**22 - 
     &    97620930228238414274520137.D0/67996359614576172780000.D0*
     &    x**23 - 3344598038287969274126620999.D0/
     &    148355693704529831520000.D0*x**24 - 
     &    82512706196227677443368639.D0/53752062936423852000000.D0*
     &    x**25 - 1482380244035590460213446634851.D0/
     &    55957729973508244381500000.D0*x**26 - 
     &    5244371136103119928918169718479.D0/
     &    3213092855078843392385730000.D0*x**27 - 1319507471775090267769
     &    2567708005281.D0/429071476647451702244740560000.D0*x**28 - 
     &    74032534058010583993741599479883671.D0/42859473056228786702002
     &    418160000.D0*x**29 - 56900527583579963202410475035999.D0/
     &    1610580246521343896630081250.D0*x**30 - 4137734010478985925770
     &    30736637841.D0/227258426509149628448630775000.D0*x**31 )
      G = G + LL**2 * (  - 26742605440057545671436486715583773.D0/
     &    664893224986769198661136896000.D0*x**32 - 98705091613285700269
     &    4449225036363.D0/516096555685294640089793760000.D0*x**33 - 
     &    983131874047781696097622823710664947.D0/2164164890173668857443
     &    2018336000.D0*x**34 - 57409787999851209637386376690605316813.D
     &    0/28662888403379671060796042466600000.D0*x**35 - 7635924601001
     &    9302190352235208820381311.D0/149871311913096319272136169760000
     &    0.D0*x**36 - 89492778888814112717513422312868915659.D0/
     &    42777554457480920843675438168640000.D0*x**37 - 157696884069412
     &    709854394813144549530330513.D0/2776976243531469778101930527780
     &    880000.D0*x**38 - 29069229000852454297209522003767558662249.D0
     &    /13334489529750210736291251993758640000.D0*x**39 - 19041667224
     &    0564827848726697130210749228033.D0/302536894752280206653752188
     &    7516800000.D0*x**40 - 3710997748407323332279976490642979407117
     &    86521.D0/163704342795276716713765752614550403200000.D0*x**41
     &     - 33140341481660090171202185916775133364213761.D0/47747099981
     &    9557090415150111792438676000.D0*x**42 )
      G = G + LL**2 * (  - 17672886922295413023206253411580249656412381.
     &    D0/7511434021551568861409068831856657220000.D0*x**43 - 
     &    405782653473377438441745344005789699214453125777.D0/5325821333
     &    680678081849355775181565873472000.D0*x**44 - 27694754018685203
     &    90113373027464676250954877.D0/11360539700076501483109935611382
     &    28800000.D0*x**45 - 476403487738951897162240913262024790278374
     &    33639.D0/571978260434597295943262648525178502800000.D0*x**46
     &     - 88368052426761920857402883983977902753692944730079.D0/
     &    35039273838571343430025081196122730045827440000.D0*x**47 - 
     &    5713035299966877974705279558143744448170308276143.D0/
     &    62984876015856813265970261671493504183040000.D0*x**48 - 
     &    711919690149995043613017226263079417675943483431.D0/
     &    273251476524523621588330833236451914853120000.D0*x**49 - 
     &    77324492105497458147996004008988301784142779083257.D0/
     &    785540621377716108394505107445917021458000000.D0*x**50 - 2974.
     &    D0/27.D0*z2 )
      G = G + LL**2 * (  - 4240.D0/9.D0*z2*x - 328.D0/27.D0*z2*x**2 - 4
     &    *z2*x**3 + 222*z2*x**4 - 4*z2*x**5 + 3252.D0/5.D0*z2*x**6 - 4
     &    *z2*x**7 + 8913.D0/7.D0*z2*x**8 - 4*z2*x**9 + 10444.D0/5.D0*
     &    z2*x**10 - 4*z2*x**11 + 34062.D0/11.D0*z2*x**12 - 4*z2*x**13
     &     + 390972.D0/91.D0*z2*x**14 - 4*z2*x**15 + 56883.D0/10.D0*z2*
     &    x**16 - 4*z2*x**17 + 123628.D0/17.D0*z2*x**18 - 4*z2*x**19 + 
     &    859578.D0/95.D0*z2*x**20 - 4*z2*x**21 + 848244.D0/77.D0*z2*
     &    x**22 - 4*z2*x**23 + 303051.D0/23.D0*z2*x**24 - 4*z2*x**25 + 
     &    5046636.D0/325.D0*z2*x**26 - 4*z2*x**27 + 379514.D0/21.D0*z2*
     &    x**28 - 4*z2*x**29 + 3017172.D0/145.D0*z2*x**30 - 4*z2*x**31
     &     + 2943273.D0/124.D0*z2*x**32 - 4*z2*x**33 + 5022084.D0/187.D0
     &    *z2*x**34 - 4*z2*x**35 + 1055882.D0/35.D0*z2*x**36 - 4*z2*
     &    x**37 + 23671452.D0/703.D0*z2*x**38 - 4*z2*x**39 + 2428923.D0/
     &    65.D0*z2*x**40 - 4*z2*x**41 + 11840484.D0/287.D0*z2*x**42 - 4
     &    *z2*x**43 + 21443946.D0/473.D0*z2*x**44 - 4*z2*x**45 + 
     &    5704924.D0/115.D0*z2*x**46 )
      G = G + LL**2 * (  - 4*z2*x**47 + 5082771.D0/94.D0*z2*x**48 - 4*
     &    z2*x**49 + 71941836.D0/1225.D0*z2*x**50 + 160.D0/9.D0*z3 + 
     &    3880.D0/27.D0*z3*x )
      G = G + LL**3 * (  - 83680.D0/243.D0 + 535418.D0/1215.D0*x - 8560.
     &    D0/729.D0*x**2 + 13402.D0/3645.D0*x**3 + 409487.D0/14580.D0*
     &    x**4 + 2567543.D0/170100.D0*x**5 + 740044.D0/18225.D0*x**6 + 
     &    3378034.D0/178605.D0*x**7 + 406203451.D0/8573040.D0*x**8 + 
     &    863029313.D0/40415760.D0*x**9 + 94164239.D0/1804275.D0*x**10
     &     + 4651691827.D0/200675475.D0*x**11 + 116663401.D0/2084940.D0
     &    *x**12 + 4676705219.D0/189729540.D0*x**13 + 19610586293.D0/
     &    332026695.D0*x**14 + 56188102913.D0/2170943775.D0*x**15 + 
     &    174992564719.D0/2835518400.D0*x**16 + 6909927797111.D0/
     &    256444284096.D0*x**17 + 3298608158777.D0/51517824930.D0*x**18
     &     + 11237409815939.D0/403051218570.D0*x**19 + 12683193047453.D0
     &    /191929151700.D0*x**20 + 9906794148281.D0/344982701700.D0*
     &    x**21 + 25778245779233.D0/379480971870.D0*x**22 + 
     &    673716963903461.D0/22859210924550.D0*x**23 + 45087622665443.D0
     &    /647722623600.D0*x**24 + 125779500498769.D0/4170120570000.D0*
     &    x**25 )
      G = G + LL**3 * ( 8058543040844147.D0/113263453338750.D0*x**26 + 
     &    254910745101198061.D0/8277293169995850.D0*x**27 + 
     &    646885997702353879.D0/8914008029226300.D0*x**28 + 
     &    9314782703312140007.D0/296803452528683100.D0*x**29 + 
     &    74593853333590333.D0/1009535552818650.D0*x**30 + 
     &    48241265180570243.D0/1510822172149290.D0*x**31 + 
     &    10376672584574933639.D0/138132312882220800.D0*x**32 + 
     &    11130956876292766417.D0/343102841675193600.D0*x**33 + 
     &    8207307852232027.D0/107599400613900.D0*x**34 + 
     &    70647280376222382583.D0/2145854846443007700.D0*x**35 + 
     &    3994945013765874943.D0/51638218229911950.D0*x**36 + 
     &    9109785041814932759.D0/272944867786677450.D0*x**37 + 
     &    38717609672069406773.D0/493900236947321100.D0*x**38 + 
     &    561201791048146069037.D0/16601277333787883100.D0*x**39 + 
     &    39625968290440816267.D0/499286536354522800.D0*x**40 + 
     &    5405171830194832430533.D0/157992183209004253200.D0*x**41 )
      G = G + LL**3 * ( 66597565714590853030199.D0/
     &    829458961847272329300.D0*x**42 + 331079937243473748217603.D0/
     &    9569124120823410042900.D0*x**43 + 23777823551940963442627.D0/
     &    292932371045614593150.D0*x**44 + 73894324536312536540293.D0/
     &    2113200499961619692910.D0*x**45 + 3623761931865295616618381.D0
     &    /44185101362833866306300.D0*x**46 + 
     &    399354023821745288065024727.D0/11306476493178488233712100.D0*
     &    x**47 + 93056713187366638220075633.D0/
     &    1123624993110905662977600.D0*x**48 + 
     &    292404077423712364497028949.D0/8200071758234907285134400.D0*
     &    x**49 + 1071006519764146686926117069.D0/
     &    12812612122242042633022500.D0*x**50 - 2956.D0/81.D0*z2 - 364.D
     &    0/81.D0*z2*x )
      G = G + LL**4 * (  - 752.D0/27.D0 + 1117.D0/27.D0*x - 128.D0/243.D
     &    0*x**2 + 4.D0/3.D0*x**4 + 4.D0/3.D0*x**6 + 4.D0/3.D0*x**8 + 4.
     &    D0/3.D0*x**10 + 4.D0/3.D0*x**12 + 4.D0/3.D0*x**14 + 4.D0/3.D0
     &    *x**16 + 4.D0/3.D0*x**18 + 4.D0/3.D0*x**20 + 4.D0/3.D0*x**22
     &     + 4.D0/3.D0*x**24 + 4.D0/3.D0*x**26 + 4.D0/3.D0*x**28 + 4.D0/
     &    3.D0*x**30 + 4.D0/3.D0*x**32 + 4.D0/3.D0*x**34 + 4.D0/3.D0*
     &    x**36 + 4.D0/3.D0*x**38 + 4.D0/3.D0*x**40 + 4.D0/3.D0*x**42
     &     + 4.D0/3.D0*x**44 + 4.D0/3.D0*x**46 + 4.D0/3.D0*x**48 + 4.D0/
     &    3.D0*x**50 )
      G = G + LL**5 * (  - 344.D0/135.D0 + 502.D0/135.D0*x )
      G = G + NF * ( 81922.D0/729.D0 + 360376.D0/2187.D0*x**(-1) - 
     &    80146.D0/729.D0*x - 1437539.D0/4374.D0*x**2 + 3977.D0/486.D0*
     &    x**3 + 1431817.D0/729000.D0*x**4 + 1758719.D0/1822500.D0*x**5
     &     + 1105053911.D0/1875352500.D0*x**6 + 2049235103.D0/
     &    5250987000.D0*x**7 + 302460143.D0/1134213192.D0*x**8 + 
     &    44460453497.D0/243045684000.D0*x**9 + 1015227229589.D0/
     &    8252392995000.D0*x**10 + 2376342760193.D0/30258774315000.D0*
     &    x**11 + 131722658273489.D0/2954601207558000.D0*x**12 + 
     &    620208831636743.D0/34219653985717200.D0*x**13 - 
     &    14035653085661.D0/4949123510331000.D0*x**14 - 13609984606079.D
     &    0/692185106340000.D0*x**15 - 16855994928234239.D0/
     &    505937578386240000.D0*x**16 - 937204623056521.D0/
     &    21063523671590400.D0*x**17 - 45468172105502033.D0/
     &    846531497246710275.D0*x**18 - 51033052557002243.D0/
     &    831756868658694000.D0*x**19 - 160951971236962309.D0/
     &    2376524335094910000.D0*x**20 )
      G = G + NF * (  - 17413715326415003747.D0/238371644821427617500.D0
     &    *x**21 - 5480520227357927255671.D0/70699048926239131398000.D0
     &    *x**22 - 188786449235765885311.D0/2322968750433571460220.D0*
     &    x**23 - 36383469260609373278503.D0/430977504718658898000000.D0
     &    *x**24 - 626638999289739765699691.D0/
     &    7197153981961793850000000.D0*x**25 - 1427627551258930519772.D0
     &    /15989153944027158984375.D0*x**26 - 261742912938538704160913.D
     &    0/2871768333103234314750000.D0*x**27 - 
     &    24609317226964328223384581.D0/265503800933229203253548640.D0*
     &    x**28 - 70271947980637631631482249.D0/
     &    747825391128947833561860000.D0*x**29 - 
     &    4889031040442766334182518969.D0/51453583710915288067114950000.
     &    D0*x**30 - 3013174660543466219896597067.D0/
     &    31429775271692097144444600000.D0*x**31 - 
     &    442739925009520658158333.D0/4585590199709596965888000.D0*
     &    x**32 )
      G = G + NF * (  - 62862872720594200521906585361.D0/
     &    647527198199850505000848384000.D0*x**33 - 
     &    4274626745919619622696291671.D0/43850006186628143752089840000.
     &    D0*x**34 - 30692968647886649161758475039.D0/
     &    313926180654269665497915900000.D0*x**35 - 99959684579169597288
     &    459788261147.D0/1020404459594484793570006465200000.D0*x**36
     &     - 2246330516497346630627506567937429.D0/229067938937874510717
     &    89203959120000.D0*x**37 - 1560123805762281510879857117381.D0/
     &    15904812761470457116804220436000.D0*x**38 - 110499903578048040
     &    684235733343851.D0/1126949120441127434447434177740000.D0*
     &    x**39 - 193388117480085971958592035753413.D0/19742559096560976
     &    62674859707200000.D0*x**40 - 387831700498280983440442547955545
     &    9.D0/39652969910137249006194246993600000.D0*x**41 - 4677245710
     &    33420648759360772549696609.D0/47916174086110059112915264823290
     &    20000.D0*x**42 - 3034096739028902775190582469162500043.D0/
     &    31157200027699857950007828297290652000.D0*x**43 )
      G = G + NF * (  - 395995402166990633021644820940056807.D0/
     &    4077682335570475162999239872340240000.D0*x**44 - 2123933261530
     &    46467470379507435058501.D0/21938036024646905158586323402956000
     &    00.D0*x**45 - 1666956106205893665617533325383388359.D0/
     &    17275809081474184987415043285525900000.D0*x**46 - 676670875078
     &    840389793526751984386461.D0/7038157315010730332082196281439410
     &    000.D0*x**47 - 16869369810691383109471370072317303474477.D0/
     &    176134843030923693150400667052139848099840.D0*x**48 - 
     &    216358923527226576138784104771617440633243.D0/2268156407292085
     &    654163241286885788768000000.D0*x**49 - 38510578255497444850885
     &    8830899498293877.D0/4054218921706387120721263525821525000000.D
     &    0*x**50 + 208.D0/27.D0*z4 + 208.D0/27.D0*z4*x + 1052.D0/81.D0
     &    *z2*x**(-1) + 8644.D0/243.D0*z2 + 6226.D0/243.D0*z2*x + 260.D0
     &    /243.D0*z2*x**2 - 64.D0/243.D0*z2*x**3 - 98.D0/1215.D0*z2*
     &    x**4 - 242.D0/6075.D0*z2*x**5 - 1024.D0/42525.D0*z2*x**6 - 
     &    968.D0/59535.D0*z2*x**7 )
      G = G + NF * (  - 841.D0/71442.D0*z2*x**8 - 1369.D0/153090.D0*z2*
     &    x**9 - 4232.D0/601425.D0*z2*x**10 - 12544.D0/2205225.D0*z2*
     &    x**11 - 8978.D0/1911195.D0*z2*x**12 - 12482.D0/3162159.D0*z2*
     &    x**13 - 33856.D0/10061415.D0*z2*x**14 - 11236.D0/3869775.D0*
     &    z2*x**15 - 14641.D0/5783400.D0*z2*x**16 - 18769.D0/8427240.D0
     &    *z2*x**17 - 23716.D0/12008817.D0*z2*x**18 - 118336.D0/
     &    67108095.D0*z2*x**19 - 72962.D0/46054575.D0*z2*x**20 - 89042.D
     &    0/62214075.D0*z2*x**21 - 215296.D0/165685905.D0*z2*x**22 - 
     &    129032.D0/108879309.D0*z2*x**23 - 76729.D0/70700850.D0*z2*
     &    x**24 - 90601.D0/90821250.D0*z2*x**25 - 212552.D0/231001875.D0
     &    *z2*x**26 - 495616.D0/582124725.D0*z2*x**27 - 287282.D0/
     &    363604059.D0*z2*x**28 - 331298.D0/450620415.D0*z2*x**29 - 
     &    760384.D0/1108669275.D0*z2*x**30 - 108578.D0/169304175.D0*z2*
     &    x**31 - 247009.D0/411000480.D0*z2*x**32 - 279841.D0/495852192.
     &    D0*z2*x**33 - 157922.D0/297411345.D0*z2*x**34 - 1420864.D0/
     &    2838926475.D0*z2*x**35 )
      G = G + NF * (  - 796322.D0/1685138175.D0*z2*x**36 - 889778.D0/
     &    1991011995.D0*z2*x**37 - 1982464.D0/4683618693.D0*z2*x**38 - 
     &    1101128.D0/2742659595.D0*z2*x**39 - 609961.D0/1599559650.D0*
     &    z2*x**40 - 674041.D0/1858597650.D0*z2*x**41 - 1486088.D0/
     &    4303368405.D0*z2*x**42 - 3268864.D0/9929235393.D0*z2*x**43 - 
     &    1793618.D0/5708445435.D0*z2*x**44 - 1964162.D0/6542902575.D0*
     &    z2*x**45 - 4293184.D0/14953229775.D0*z2*x**46 - 1170724.D0/
     &    4259404845.D0*z2*x**47 - 1274641.D0/4839671592.D0*z2*x**48 - 
     &    1385329.D0/5484364200.D0*z2*x**49 - 1503076.D0/6199081875.D0*
     &    z2*x**50 - 2800.D0/81.D0*z3*x**(-1) - 160.D0/81.D0*z3 - 688.D0
     &    /81.D0*z3*x + 80.D0/3.D0*z3*x**2 )
      G = G + NF*LL * ( 195812.D0/729.D0 + 250016.D0/729.D0*x + 3383.D0/
     &    81.D0*x**2 + 251.D0/729.D0*x**3 + 109177.D0/36450.D0*x**4 + 
     &    669703.D0/182250.D0*x**5 + 17642171.D0/4465125.D0*x**6 + 
     &    25586674.D0/6251175.D0*x**7 + 93996677.D0/22504230.D0*x**8 + 
     &    102015766.D0/24111675.D0*x**9 + 8894873816.D0/2083937625.D0*
     &    x**10 + 4688728181.D0/1091586375.D0*x**11 + 743068926307.D0/
     &    172179557550.D0*x**12 + 1233929854321.D0/284878904310.D0*
     &    x**13 + 1968743054221.D0/453216438675.D0*x**14 + 
     &    1517956279157.D0/348628029750.D0*x**15 + 7025504707693.D0/
     &    1610445564000.D0*x**16 + 37596612209917.D0/8604380584800.D0*
     &    x**17 + 2206307721077.D0/504250225830.D0*x**18 + 
     &    19009369599379.D0/4339518119937.D0*x**19 + 587646519245893.D0/
     &    134014530174525.D0*x**20 + 264845901853952.D0/60345724113675.D
     &    0*x**21 + 81175430098739677.D0/18481672032498675.D0*x**22 + 
     &    10676101484101342.D0/2429019752842683.D0*x**23 + 
     &    231226067991452731.D0/52576185126465000.D0*x**24 )
      G = G + NF*LL * ( 14152039790023259.D0/3216122749125000.D0*x**25
     &     + 5104878493970742364.D0/1159534603555453125.D0*x**26 + 
     &    1170006668961716059.D0/265638836450885625.D0*x**27 + 
     &    233221853698108330271.D0/52929151175538462825.D0*x**28 + 
     &    26285784613350701024.D0/5963268780460821375.D0*x**29 + 
     &    3151499558899778897071.D0/714713313812373500625.D0*x**30 + 
     &    3369907010827990760024.D0/764003887178744086875.D0*x**31 + 
     &    850194533862564514003.D0/192694576470698016000.D0*x**32 + 
     &    312259412396985560483.D0/70753772959367097600.D0*x**33 + 
     &    2154405401266516342999.D0/488036987373457021500.D0*x**34 + 
     &    56566328453862906759517.D0/12810970918553246814375.D0*x**35
     &     + 451858950281236776156763.D0/102313383689836042132500.D0*
     &    x**36 + 800578277625197132701.D0/181236210635744491500.D0*
     &    x**37 + 314094367587642614990489.D0/71091629681019688643175.D0
     &    *x**38 + 26280298892691664088269.D0/5947176228141930195375.D0
     &    *x**39 )
      G = G + NF*LL * ( 35196861325261557451509113.D0/
     &    7963632737018047523070000.D0*x**40 + 449486352899980564504883.
     &    D0/101684503960678900170000.D0*x**41 + 
     &    72732187147815744821076257.D0/16451281684538337561253875.D0*
     &    x**42 + 1846240392777640185849981671.D0/
     &    417541554169038879664409325.D0*x**43 + 
     &    29723962582544591580341750851.D0/6721400628086967331183174500.
     &    D0*x**44 + 4867619674949061174676805737.D0/
     &    1100561564636926307194207500.D0*x**45 + 
     &    3533213863721674534576367551.D0/798757662626540124701675625.D0
     &    *x**46 + 37242221582988905788334351422.D0/
     &    8418421668530625920340690375.D0*x**47 + 1658961255132274076456
     &    2609651271.D0/3749589412385768409922411147200.D0*x**48 + 
     &    1709225460060903184276161749701.D0/
     &    386279289393033916588592520000.D0*x**49 + 
     &    2656724160212188321459136679673.D0/
     &    600350956752803710123560515625.D0*x**50 )
      G = G + NF*LL * (  - 472.D0/81.D0*z2 - 664.D0/81.D0*z2*x - 64.D0/
     &    81.D0*z2*x**2 - 640.D0/27.D0*z3 - 640.D0/27.D0*z3*x )
      G = G + NF*LL**2 * ( 5996.D0/243.D0 + 7019.D0/243.D0*x - 1940.D0/
     &    243.D0*x**2 - 128.D0/243.D0*x**3 - 196.D0/1215.D0*x**4 - 484.D
     &    0/6075.D0*x**5 - 2048.D0/42525.D0*x**6 - 1936.D0/59535.D0*
     &    x**7 - 841.D0/35721.D0*x**8 - 1369.D0/76545.D0*x**9 - 8464.D0/
     &    601425.D0*x**10 - 25088.D0/2205225.D0*x**11 - 17956.D0/
     &    1911195.D0*x**12 - 24964.D0/3162159.D0*x**13 - 67712.D0/
     &    10061415.D0*x**14 - 22472.D0/3869775.D0*x**15 - 14641.D0/
     &    2891700.D0*x**16 - 18769.D0/4213620.D0*x**17 - 47432.D0/
     &    12008817.D0*x**18 - 236672.D0/67108095.D0*x**19 - 145924.D0/
     &    46054575.D0*x**20 - 178084.D0/62214075.D0*x**21 - 430592.D0/
     &    165685905.D0*x**22 - 258064.D0/108879309.D0*x**23 - 76729.D0/
     &    35350425.D0*x**24 - 90601.D0/45410625.D0*x**25 - 425104.D0/
     &    231001875.D0*x**26 - 991232.D0/582124725.D0*x**27 - 574564.D0/
     &    363604059.D0*x**28 - 662596.D0/450620415.D0*x**29 - 1520768.D0
     &    /1108669275.D0*x**30 - 217156.D0/169304175.D0*x**31 - 247009.D
     &    0/205500240.D0*x**32 )
      G = G + NF*LL**2 * (  - 279841.D0/247926096.D0*x**33 - 315844.D0/
     &    297411345.D0*x**34 - 2841728.D0/2838926475.D0*x**35 - 1592644.
     &    D0/1685138175.D0*x**36 - 1779556.D0/1991011995.D0*x**37 - 
     &    3964928.D0/4683618693.D0*x**38 - 2202256.D0/2742659595.D0*
     &    x**39 - 609961.D0/799779825.D0*x**40 - 674041.D0/929298825.D0
     &    *x**41 - 2972176.D0/4303368405.D0*x**42 - 6537728.D0/
     &    9929235393.D0*x**43 - 3587236.D0/5708445435.D0*x**44 - 
     &    3928324.D0/6542902575.D0*x**45 - 8586368.D0/14953229775.D0*
     &    x**46 - 2341448.D0/4259404845.D0*x**47 - 1274641.D0/
     &    2419835796.D0*x**48 - 1385329.D0/2742182100.D0*x**49 - 
     &    3006152.D0/6199081875.D0*x**50 - 16.D0/27.D0*z2 - 16.D0/27.D0
     &    *z2*x )
      G = G + NF*LL**3 * ( 1000.D0/243.D0 + 1480.D0/243.D0*x - 128.D0/
     &    243.D0*x**2 )
      G = G + NF*LL**4 * ( 56.D0/81.D0 + 56.D0/81.D0*x )
      G = G - 617302.D0/1215.D0 - 15688534.D0/2187.D0*x**(-1) - 
     &    46388051359.D0/5467500.D0*x + 4939872908609.D0/562605750.D0*
     &    x**2 + 286488983.D0/2480625.D0*x**3 + 13194688587203.D0/
     &    88914672000.D0*x**4 - 72672035925126988267.D0/
     &    128852919795600000.D0*x**5 - 171062678786066236159.D0/
     &    193279379693400000.D0*x**6 - 1173892072693239587522581.D0/
     &    469205022143697840000.D0*x**7 - 
     &    1362664747703225965816196695997.D0/
     &    536605385164595747320320000.D0*x**8 - 
     &    2759744790663629500887054592243.D0/
     &    505942220298047418902016000.D0*x**9 - 970630931725201112194562
     &    00085167.D0/19697490087124502376523800000.D0*x**10 - 406722037
     &    2150108346395901237955352911.D0/431900299310349922109245188000
     &    000.D0*x**11 - 18113214801181912638007990158791207.D0/
     &    2243637918495324270697377600000.D0*x**12 - 2469780100763435783
     &    46746947275761792381.D0/17170785354036566176074100510560000.D0
     &    *x**13
      G = G - 156882989509321293950728472189003492439905123.D0/
     & 13090892282015117808864387068580873600000.D0*x**14 - 870943306259
     &    8643625930532949939223798051199627.D0/427971478450494236059028
     &    038780528560000000.D0*x**15 - 22300760816252417344334936067079
     &    664204916882540643.D0/1338205674332059691265715101832601303040
     &    000000.D0*x**16 - 10775498563832382588552424911294068060692296
     &    16271937667.D0/39449143501058031913447517582270164825823232000
     &    000.D0*x**17 - 14723437273381447751514046825168385770345456029
     &    178363.D0/665704296580354288539426859200809031435767040000.D0
     &    *x**18 - 11610086834036549333566638799029948172035595130309833
     &    7881.D0/3291320360445810467749342795048705828951542336000000.D
     &    0*x**19 - 4187545627736535624143597738368971307559127560487294
     &    7223849051.D0/147746091022494480748752762214204996831590143196
     &    5760000000.D0*x**20 - 1107450055400207825201116600248255080143
     &    818596124315296542189993.D0/2503907437328590673742020496472316
     &    2620932645320682880000000.D0*x**21
      G = G - 1535968307168818446154891196523690069944258642238443437660
     & 214291419.D0/4345889058546244627232668288841486632328873918909095
     & 0080000000.D0*x**22 - 7019664160762071133069812594827297290765996
     &    2924321314874917666804941.D0/129571877486286182404529554537681
     &    3607046201298045119067200000000.D0*x**23 - 7836204835778761864
     &    1421722747142166919471706445248109019874922486861.D0/
     &    18173717881193386622972976480609853189738927297256215488000000
     &    00.D0*x**24 - 107406821557232074424542027062948440882309137558
     &    628412967719581220901.D0/1649497671259377521064826130061038799
     &    888094103099076840000000000.D0*x**25 - 15321646646989357326883
     &    576930078439179059644460080333157778966156556091302687.D0/
     &    29653237204972668808077407790673849556605972931024303120915196
     &    0000000000.D0*x**26 - 1669616688186013523210054018512126815605
     &    488810126511826842108046259739003271353.D0/2167058574939402636
     &    4942969613424449255967645017992560720764825236800000000.D0*
     &    x**27
      G = G - 1600131451441508753298189564862501978009689645647188220319
     & 974651067341002214396051.D0/2623311226221155681646077813730724391
     & 9270948382622061762130082360413900000000.D0*x**28 - 
     &    22456795218798034447909115416982456810921194824995680185614009
     &    4147251563329151894384091.D0/249612675539036710171680567035958
     &    9231121717243438401022389491108418190494950400000.D0*x**29 - 
     &    10564028024144188480676029519308493655786748207272003708793619
     &    797543289153637446196748557.D0/1485789735351408989117146232356
     &    89835185816502585619108475564946929654196128000000000.D0*
     &    x**30 - 240503534311760090772265425081682837502282303852703098
     &    33691069045864329722137904958729.D0/23152448876994912366860125
     &    2654387498278893715700924518082965384277635472000000000.D0*
     &    x**31 - 475021916832211896807471260891503363449219486197732066
     &    09856447686679481679736291888866581157.D0/
     &    57942878520122913213177401042616011103486567081860048124640311
     &    6271991194243891200000000.D0*x**32
      G = G - 2197943619196809039815261316315579317184314504238036513185
     & 34212370357408026243857611658481053.D0/18504338624168285187434073
     & 23619027451369409722936820891722384145513778330004684800000000.D0
     & *x**33 - 30045877778704805033380859672312428453140219155452569285
     &    653034942933786671844336367381159167846763.D0/
     &    32086990119727775761582209648245755495703996095209246600863036
     &    6535831357347816965516800000000.D0*x**34 - 
     &    11934770545105293179475200594426322797919035031086897032903142
     &    952328828135992977484271553545172468729.D0/
     &    88621489044363142320191189563994401853963356521651390671410106
     &    823041709234264062952576000000000.D0*x**35 - 
     &    16149376568827003576251063982207170974480727399037457554888165
     &    79254726756188213111794536891937110510967.D0/
     &    15224650514944856302712845242449744153794116630675463909167836
     &    292746665401098129167823424000000000.D0*x**36 - 
     &    19513403366206458933638782332254480334030718971774340403718585
     &    99350328983756690284452823657231259670887.D0/
     &    12875704435496221330294291976471783627208738636228392334610512
     &    979008608453500132096216381440000000.D0*x**37
      G = G - 9877154581858128006552810657035646143916872509065701215638
     & 2100441426948574784916296455355130908559224057429.D0/
     & 82802003183234078828205874046825334470849452074812682434493360793
     & 7200172673429612663054683673600000000.D0*x**38 - 
     &    12125169475555902693611188536658349323481092716255165463614368
     &    76050365994536790357029855570712134224475849561.D0/
     &    71567785454049347052054698703174978280480148036554313087975612
     &    92386935546512507841341753725373440000000.D0*x**39 - 
     &    16265245594866454309251287292238945696351318344291827799660647
     &    6351152900113750915884976874122277974941587319479963.D0/
     &    12204126042095125946079198242534366779411393888731136638304462
     &    43644682014223464195906250474148304776000000000.D0*x**40 - 
     &    58905027451638170602906287633473010468456006618681403262452557
     &    784481120945339369190643320152764834930514704400805071.D0/
     &    31285659354338658731697684754013189572157821378686080052619875
     &    8176052361961062442410729754785365359104000000000.D0*x**41
      G = G - 3841972135838434406504464189224604223747934253494126474305
     & 230041238950875876845120821773621493147427202492785491099283.D0/
     & 25951454434423917417943229503453940750104912833620103403648186990
     & 703543424670129597970033159446056537676800000000.D0*x**42 - 
     &    44870329686370877094366823297574745409693418721681542456604620
     &    647481096938938301024545244480128237160375498722409814141629.D
     &    0/215588859383611940104968726794266404422145942094378245843726
     &    425209624365879827320767253346444068419387861734400000000.D0*
     &    x**43 - 723559480656449064130985822610981098033533599669403473
     &    09518019687221219254644313829885609516642307529526606045225079
     &    310709.D0/4422993662751086231518604434628402027232121113602918
     &    77385740324735618242697899701732817381395172114379065542400000
     &    000.D0*x**44 - 73192279385033482784083915211678608871606676696
     &    17107659377845566369533217597812933516005674674935788291329979
     &    140595245736317.D0/3196646772490438592005947239227399433824904
     &    61879189592255816500478120294247917259515298675077630800083740
     &    39296000000000.D0*x**45
      G = G - 4208850754608925940470287155227205541618398410670827996693
     & 62927310885006266120112793698834357163246487247178131328526813227
     & 77.D0/23393642289589118786952613887073241311173165619340692887812
     & 0257168078942608703085372559484943175267334010014848000000000.D0
     & *x**46 - 43946090238715131020235243235033463415845006356905177005
     &    47878686171254868865728406552860640727042522786286291600121821
     &    652431438653.D0/1752299204573104357985991703430732922161170010
     &    44485078331255596289255828954944380500194148033278262283768066
     &    18800788480000000.D0*x**47 - 540362209060569675406071655947081
     &    87179395620500470778213668361936805839044812763634641590395221
     &    930940857112390638633850740535116381.D0/
     &    27427291897665981255432913618915819651218313206962881825761745
     &    50614439061903477260003038838781746714006804514247079936000000
     &    00.D0*x**48 - 199520471291013584798863156913098596946203496179
     &    01396791996270804847917461136988900141576258443858767319988765
     &    728146761251791137776173.D0/7292289589232150566248343850798645
     &    49832011842456765399878999081356423803539037874965563276875647
     &    50704765278172440886113280000000.D0*x**49
      G = G - 5912081997326061927248578503258334631532005008610786586484
     & 82026607420414522645123404824397524325799438459088087376656695019
     & 169877859744967.D0/2751131718575339455138499755364673909594029272
     & 28302924548955962471716255355957863138762763440558543558932543712
     & 8659696800640000000000.D0*x**50 + 75616.D0/81.D0*z4*x**(-1) - 
     &    44002.D0/27.D0*z4 + 2468*z4*x - 520.D0/9.D0*z4*x**2 - 10*z4*
     &    x**3 + 1455.D0/2.D0*z4*x**4 - 919.D0/10.D0*z4*x**5 + 2334*z4*
     &    x**6 - 1564.D0/7.D0*z4*x**7 + 130755.D0/28.D0*z4*x**8 - 1613.D
     &    0/4.D0*z4*x**9 + 7728*z4*x**10 - 34714.D0/55.D0*z4*x**11 + 
     &    253155.D0/22.D0*z4*x**12 - 23585.D0/26.D0*z4*x**13 + 1456590.D
     &    0/91.D0*z4*x**14 - 43088.D0/35.D0*z4*x**15 + 169809.D0/8.D0*
     &    z4*x**16 - 218017.D0/136.D0*z4*x**17 + 461820.D0/17.D0*z4*
     &    x**18 - 38438.D0/19.D0*z4*x**19 + 1285377.D0/38.D0*z4*x**20
     &     - 174373.D0/70.D0*z4*x**21 + 3172830.D0/77.D0*z4*x**22 - 
     &    760780.D0/253.D0*z4*x**23 + 4536105.D0/92.D0*z4*x**24 - 
     &    357103.D0/100.D0*z4*x**25
      G = G + 3778152.D0/65.D0*z4*x**26 - 163138.D0/39.D0*z4*x**27 + 
     &    947315.D0/14.D0*z4*x**28 - 1966267.D0/406.D0*z4*x**29 + 
     &    2259834.D0/29.D0*z4*x**30 - 860408.D0/155.D0*z4*x**31 + 
     &    44097015.D0/496.D0*z4*x**32 - 1110035.D0/176.D0*z4*x**33 + 
     &    18813180.D0/187.D0*z4*x**34 - 4231054.D0/595.D0*z4*x**35 + 
     &    1582353.D0/14.D0*z4*x**36 - 589263.D0/74.D0*z4*x**37 + 
     &    88694130.D0/703.D0*z4*x**38 - 2189164.D0/247.D0*z4*x**39 + 
     &    7281309.D0/52.D0*z4*x**40 - 8045029.D0/820.D0*z4*x**41 + 
     &    44371680.D0/287.D0*z4*x**42 - 3252910.D0/301.D0*z4*x**43 + 
     &    160730265.D0/946.D0*z4*x**44 - 1303611.D0/110.D0*z4*x**45 + 
     &    4276278.D0/23.D0*z4*x**46 - 13991392.D0/1081.D0*z4*x**47 + 
     &    76202085.D0/376.D0*z4*x**48 - 5520539.D0/392.D0*z4*x**49 + 
     &    53930652.D0/245.D0*z4*x**50 + 224.D0/27.D0*ln2**4*x**(-1) + 
     &    928.D0/81.D0*ln2**4 - 800.D0/81.D0*ln2**4*x + 32.D0/9.D0*
     &    ln2**4*x**2 + 1792.D0/9.D0*li4half*x**(-1) + 7424.D0/27.D0*
     &    li4half
      G = G - 6400.D0/27.D0*li4half*x + 256.D0/3.D0*li4half*x**2 - 
     &    10174.D0/81.D0*z2*x**(-1) + 52994.D0/243.D0*z2 - 503738.D0/
     &    2025.D0*z2*x + 5433964.D0/8505.D0*z2*x**2 + 5162518.D0/18225.D
     &    0*z2*x**3 - 102376249.D0/145800.D0*z2*x**4 + 17525223367.D0/
     &    11907000.D0*z2*x**5 - 305231107.D0/297675.D0*z2*x**6 + 
     &    5372311859.D0/1786050.D0*z2*x**7 - 513856337717.D0/377213760.D
     &    0*z2*x**8 + 27432958463861.D0/5334880320.D0*z2*x**9 - 
     &    642652853729819.D0/379276647750.D0*z2*x**10 + 
     &    10382122248509177.D0/1314825712200.D0*z2*x**11 - 
     &    2784020541620449.D0/1377436460400.D0*z2*x**12 + 
     &    128486120627320249.D0/11395156172400.D0*z2*x**13 - 
     &    48146542386918319.D0/20545811886600.D0*z2*x**14 + 
     &    67748839012221044779.D0/4433154026301000.D0*z2*x**15 - 
     &    186361562062965493903.D0/70009289558208000.D0*z2*x**16 + 
     &    424662773154018801687899.D0/21320829068598259200.D0*z2*x**17
     &     - 2922419709173026382.D0/981467036775225.D0*z2*x**18
      G = G + 262533699333637822085551.D0/10425258331336648800.D0*z2*
     & x**19 - 1174232604097040050169.D0/356816877411819825.D0*z2*x**20
     &     + 1063012748613929390828477.D0/34205858232269724000.D0*z2*
     &    x**21 - 14908827490179401640363347.D0/4138908846104636604000.D
     &    0*z2*x**22 + 365251421509169249689677373.D0/
     &    9713765659225167540000.D0*z2*x**23 - 
     &    96720455139104347826093807.D0/24725948950754971920000.D0*z2*
     &    x**24 + 10825823275493038371129089537.D0/
     &    241884283213907334000000.D0*z2*x**25 - 
     &    24737765107314317233582778401.D0/5862238378177054173300000.D0
     &    *z2*x**26 + 67528227365219970087624086377681.D0/
     &    1285237142031537356954292000.D0*z2*x**27 - 3468403345405319218
     &    5549243313667.D0/7661990654418780397227510000.D0*z2*x**28 + 
     &    93307069188354827818719408492330619.D0/
     &    1530695466293885239357229220000.D0*z2*x**29 - 6322430871791344
     &    2019150323811799.D0/13082867233281070421856660000.D0*z2*x**30
      G = G + 63635884167999080655288268179503761.D0/
     & 909033706036598513794523100000.D0*z2*x**31 - 30655376074789749127
     &    26478562701073.D0/596699048065049280849738240000.D0*z2*x**32
     &     + 4537735613619457002975063422043071.D0/
     &    56948585454929063734046208000.D0*z2*x**33 - 173318257611865330
     &    917871859939279.D0/31850935393058313247299870000.D0*z2*x**34
     &     + 101151404173394693000906574385666961039.D0/1124034839348222
     &    394541021273200000.D0*z2*x**35 - 17882042564251269103683893300
     &    5716713.D0/31127118628104620156520589104000.D0*z2*x**36 + 
     &    269843022258001921023962614401376612223.D0/2673597153592557552
     &    729714885540000.D0*z2*x**37 - 42089420032546867406945948013099
     &    023659.D0/6959840209352054581709099067120000.D0*z2*x**38 + 
     &    12606036263806518209581563483899172027817.D0/11205453386344714
     &    9044464302468560000.D0*z2*x**39 - 6206158502981039500951817935
     &    8720630287.D0/9774268907381360522659686098131200.D0*z2*x**40
     &     + 35440992999217700421181407148683673777725319.D0/28420892846
     &    4022077628065542733594450000.D0*z2*x**41
      G = G - 190538145399282429124108066803815432413756487.D0/
     & 28648259989173425424909006707546320560000.D0*z2*x**42 + 
     &    1060919106282670558275107243423585326093463.D0/771392454074615
     &    5441755141290738544000.D0*z2*x**43 - 1542692816918387400204231
     &    4448662726747372753241.D0/221909222236694920077056490632565244
     &    7280000.D0*z2*x**44 + 1177217082175136810277076484402767659531
     &    228374581.D0/7796312375923750099052732187506237418600000.D0*
     &    z2*x**45 - 47289696394649725018661679560669530804360540673.D0/
     &    6520552168954409173753194193187034931920000.D0*z2*x**46 + 
     &    113425233787021131601059073127879936559343423642333.D0/
     &    687044585070026341765197670512210393055440000.D0*z2*x**47 - 
     &    130562259858228590324166901875963220146735749287.D0/
     &    17287338310735167896404603735367366041728000.D0*z2*x**48 + 
     &    9890897980537277575306383117195925660249377624359.D0/
     &    55005167352339170579469193703441619223680000.D0*z2*x**49 - 
     &    15429113769150995135460570644862358457903278545841.D0/
     &    1964951750392998556852361515263876471042000000.D0*z2*x**50
      G = G - 448.D0/9.D0*z2*ln2**2*x**(-1) - 1856.D0/27.D0*z2*ln2**2
     &     + 1600.D0/27.D0*z2*ln2**2*x - 64.D0/3.D0*z2*ln2**2*x**2 + 
     &    90950.D0/81.D0*z3*x**(-1) - 107233.D0/405.D0*z3 + 225994.D0/
     &    135.D0*z3*x - 76408.D0/1215.D0*z3*x**2 + 45562.D0/1215.D0*z3*
     &    x**3 - 2063797.D0/2430.D0*z3*x**4 - 4135897.D0/85050.D0*z3*
     &    x**5 - 60840952.D0/42525.D0*z3*x**6 - 8145077.D0/59535.D0*z3*
     &    x**7 - 2799580727.D0/1428840.D0*z3*x**8 - 1544673509.D0/
     &    6735960.D0*z3*x**9 - 10412406392.D0/4209975.D0*z3*x**10 - 
     &    64966207612.D0/200675475.D0*z3*x**11 - 6133347103.D0/2058210.D
     &    0*z3*x**12 - 13252525703.D0/31621590.D0*z3*x**13 - 
     &    77077079143.D0/22135113.D0*z3*x**14 - 745299866351.D0/
     &    1447295850.D0*z3*x**15 - 13170144371561.D0/3308104800.D0*z3*
     &    x**16 - 391788494697293.D0/641110710240.D0*z3*x**17 - 
     &    269141361039073.D0/60104129085.D0*z3*x**18 - 47521485529811.D0
     &    /67175203095.D0*z3*x**19 - 159072986807791.D0/31988191950.D0*
     &    z3*x**20
      G = G - 46218712822721.D0/57497116950.D0*z3*x**21 - 
     &    3803040311687081.D0/695715115095.D0*z3*x**22 - 
     &    10290329929901356.D0/11429605462275.D0*z3*x**23 - 
     &    49531387817001611.D0/8312440336200.D0*z3*x**24 - 
     &    9006673283366881.D0/9035261235000.D0*z3*x**25 - 
     &    121759138777261939.D0/18877242223125.D0*z3*x**26 - 
     &    116027327862898873.D0/106119143205075.D0*z3*x**27 - 
     &    10311358499908928873.D0/1485668004871050.D0*z3*x**28 - 
     &    176584622348341948837.D0/148401726264341550.D0*z3*x**29 - 
     &    17502774172219450189.D0/2355582956576850.D0*z3*x**30 - 
     &    588969754390827851.D0/457824900651300.D0*z3*x**31 - 
     &    1276260634149970475401.D0/161154365029257600.D0*z3*x**32 - 
     &    79084712015459573989.D0/57183806945865600.D0*z3*x**33 - 
     &    127714328960999524193.D0/15189448719995550.D0*z3*x**34 - 
     &    317484502493417086999.D0/214585484644300770.D0*z3*x**35 - 
     &    76564074421510895831.D0/8606369704985325.D0*z3*x**36
      G = G - 71695776483581632688.D0/45490811297779575.D0*z3*x**37 - 
     &    5407156363517403423757.D0/576216943105207950.D0*z3*x**38 - 
     &    243567709198423159871.D0/145625239770069150.D0*z3*x**39 - 
     &    5749971059930435537753.D0/582500959080276600.D0*z3*x**40 - 
     &    139749657535153519776131.D0/78996091604502126600.D0*z3*x**41
     &     - 15751323564919814341550543.D0/1520674763386665937050.D0*z3
     &    *x**42 - 2975304102324479576732701.D0/1594854020137235007150.D
     &    0*z3*x**43 - 40768582109619927919339693.D0/
     &    3759298761752053945425.D0*z3*x**44 - 
     &    3455168968225119947482124.D0/1761000416634683077425.D0*z3*
     &    x**45 - 83443966785602481528338377.D0/7364183560472311051050.D
     &    0*z3*x**46 - 5818657024729579405924885286.D0/
     &    2826619123294622058428025.D0*z3*x**47 - 
     &    15491002631421621635403605287.D0/1310895825296056606807200.D0
     &    *z3*x**48 - 2945170344579008019499961657.D0/
     &    1366678626372484547522400.D0*z3*x**49
      G = G - 26272098385583833496156387083.D0/
     & 2135435353707007105503750.D0*z3*x**50 - 688.D0/9.D0*z3*z2 - 968.D
     &    0/9.D0*z3*z2*x - 17264.D0/27.D0*z5 + 75376.D0/27.D0*z5*x
*
      AGGREG0=G*as**3
      RETURN
      END
      REAL*8 FUNCTION AGGREG1(z,nf,as,LL)
*     -----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of agg3 in the range: representation for  z \in [1/2,1]
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(33)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half,G,Lx,z4
      REAL*8 B4,y,Ly
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
      y=1.0D0-z
      Ly = LOG(1.0D0-z)
      w(1)=1.D0/49.D0*y
      w(2)= - 1294897166440581626339803701292983236251182720763.D0/144.D
     & 0 - 1139082380824112146858502531679383940015709617589.D0/125.D0*
     & y
      w(2)=w(2)*w(1)
      w(2)= - 612654786754690178779177558843461399324656946031.D0/3384.D
     & 0 + w(2)
      w(2)=y*w(2)
      w(2)= - 1875772444502617639095150511727337613430282291.D0/10810.D0
     &  + 1.D0/1029.D0*w(2)
      w(3)=1.D0/34.D0*y
      w(2)=w(2)*w(3)
      w(2)= - 234359430040253239973876758285401432831336457.D0/46575.D0
     &  + w(2)
      w(4)=1.D0/2209.D0*y
      w(2)=w(2)*w(4)
      w(2)= - 6801548850126460889178591071763579942390733.D0/3029400.D0
     &  + w(2)
      w(2)=y*w(2)
      w(2)= - 2490365678266864674489918929883207872679793.D0/1125740.D0
     &  + w(2)
      w(2)=y*w(2)
      w(2)= - 8429176761069408915129094645931023825452481.D0/3868452.D0
     &  + w(2)
      w(2)=y*w(2)
      w(2)= - 9891784090223710940166955282096857286152463.D0/4610655.D0
     &  + w(2)
      w(5)=1.D0/1849.D0*y
      w(2)=w(2)*w(5)
      w(2)= - 6686266862489513607865316621815637977529.D0/5854800.D0 + 
     & w(2)
      w(2)=y*w(2)
      w(2)= - 9385888713863091167041887913328318358427.D0/8353800.D0 + 
     & w(2)
      w(6)=1.D0/1681.D0*y
      w(2)=w(2)*w(6)
      w(2)= - 10432632816153009687276816741649605061.D0/15872220.D0 + 
     & w(2)
      w(7)=1.D0/5.D0*y
      w(2)=w(2)*w(7)
      w(2)= - 1544276174491475060922617411458019.D0/11951.D0 + w(2)
      w(8)=1.D0/11.D0*y
      w(2)=w(2)*w(8)
      w(2)= - 823220039941400081598524162981294293.D0/71328600.D0 + 
     & w(2)
      w(2)=y*w(2)
      w(2)= - 4205770359050612438936307635679111521.D0/371101500.D0 + 
     & w(2)
      w(2)=y*w(2)
      w(2)= - 27513167559871926418078433051663.D0/77885500.D0 + 1.D0/
     & 31487.D0*w(2)
      w(2)=y*w(2)
      w(2)= - 5267019860354915449454210765123999.D0/15201024300.D0 + 
     & w(2)
      w(2)=y*w(2)
      w(2)= - 12148203127196443834252308345223.D0/35767116.D0 + w(2)
      w(9)=1.D0/7.D0*y
      w(2)=w(2)*w(9)
      w(2)= - 9724781973152851190984331517.D0/204600.D0 + w(2)
      w(10)=1.D0/29.D0*y
      w(2)=w(2)*w(10)
      w(2)= - 3076165588153015989698519844787.D0/1917357750.D0 + w(2)
      w(11)=1.D0/4.D0*y
      w(2)=w(2)*w(11)
      w(2)= - 234678733093744733257912001922091.D0/598184692875.D0 + 
     & w(2)
      w(12)=1.D0/961.D0*y
      w(2)=w(2)*w(12)
      w(2)= - 2672970821932024075762939733243.D0/6699668560200.D0 + 
     & w(2)
      w(2)=y*w(2)
      w(2)= - 80996723223625928955202553729.D0/207920748420.D0 + w(2)
      w(2)=w(2)*w(10)
      w(2)= - 15860630502803576705639689.D0/1210465620.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 51138110452244775787787.D0/684934250.D0 + 1.D0/171.D0*
     & w(2)
      w(2)=y*w(2)
      w(2)= - 748160797060998007279783.D0/10296594000.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 238923635279210230045259.D0/3383166600.D0 + w(2)
      w(2)=w(2)*w(7)
      w(2)= - 166643565630675968618851.D0/12156845316.D0 + w(2)
      w(13)=1.D0/13.D0*y
      w(2)=w(2)*w(13)
      w(2)= - 14175706081506697204811.D0/13874660415.D0 + w(2)
      w(14)=1.D0/253.D0*y
      w(2)=w(2)*w(14)
      w(2)= - 10125684853140147517.D0/2592462600.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 4814809353988258331893.D0/1277158182300.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 2781928716581648790137.D0/766294909380.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 478118936643733874087.D0/137126457468.D0 + w(2)
      w(15)=1.D0/19.D0*y
      w(2)=w(2)*w(15)
      w(2)= - 33548590873985611.D0/190930740.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 169457065390705487.D0/1010809800.D0 + w(2)
      w(16)=1.D0/17.D0*y
      w(2)=w(2)*w(16)
      w(2)= - 974927063827007.D0/104053950.D0 + w(2)
      w(17)=1.D0/16.D0*y
      w(2)=w(2)*w(17)
      w(2)= - 1664125232426.D0/3006003.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 2251790462534603.D0/4328644320.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 111073455063829.D0/228918690.D0 + w(2)
      w(2)=w(2)*w(13)
      w(2)= - 24293466248171.D0/704365200.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 6049588496249.D0/192099600.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 1797342233.D0/63504.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 819520883.D0/32928.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 174520021.D0/8232.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 12603559.D0/735.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 186629.D0/15.D0 + w(2)
      w(2)=y*w(2)
      w(2)= - 221173.D0/40.D0 + w(2)
      w(18)=1.D0/9.D0*y
      w(2)=w(2)*w(18)
      w(2)=2018 + w(2)
      w(2)=y*w(2)
      w(2)= - 1268 + w(2)
      w(2)=w(2)*w(7)
      w(2)= - 10603.D0/6.D0 + w(2)
      w(19)=1.D0/3.D0*y
      w(2)=w(2)*w(19)
      w(20)=1656184863679869711274837.D0/36.D0 + 
     & 1151588487643419515283296.D0/25.D0*y
      w(20)=y*w(20)
      w(20)=113320841445612293182819.D0/846.D0 + 1.D0/343.D0*w(20)
      w(20)=w(20)*w(7)
      w(20)=86760192847481975047271.D0/3243.D0 + w(20)
      w(21)=1.D0/47.D0*y
      w(20)=w(20)*w(21)
      w(20)=8824434870808178453632.D0/15525.D0 + w(20)
      w(20)=y*w(20)
      w(20)=8428207586821300351979.D0/14850.D0 + w(20)
      w(20)=y*w(20)
      w(20)=1340187525930365961491.D0/2365.D0 + w(20)
      w(20)=w(20)*w(8)
      w(20)=4876573416966328652621.D0/94815.D0 + w(20)
      w(22)=1.D0/43.D0*y
      w(20)=w(20)*w(22)
      w(20)=21590165695504134788.D0/18081.D0 + w(20)
      w(20)=y*w(20)
      w(20)=34209424178109984101.D0/28700.D0 + w(20)
      w(23)=1.D0/41.D0*y
      w(20)=w(20)*w(23)
      w(20)=1188304912837559243.D0/40950.D0 + w(20)
      w(20)=y*w(20)
      w(20)=2253390372129062279.D0/77805.D0 + w(20)
      w(20)=y*w(20)
      w(20)=2133469872363516712.D0/73815.D0 + w(20)
      w(20)=y*w(20)
      w(20)=57624226845085091.D0/1998.D0 + w(20)
      w(24)=1.D0/37.D0*y
      w(20)=w(20)*w(24)
      w(20)=25723299508207387.D0/33075.D0 + w(20)
      w(20)=y*w(20)
      w(20)=177734046809297831.D0/229075.D0 + w(20)
      w(20)=y*w(20)
      w(20)=501471957347376602.D0/647955.D0 + w(20)
      w(20)=y*w(20)
      w(20)=941443785714565457.D0/1219680.D0 + w(20)
      w(20)=y*w(20)
      w(20)=29396785981438543.D0/38192.D0 + w(20)
      w(20)=y*w(20)
      w(20)=412161479445793754.D0/537075.D0 + w(20)
      w(25)=1.D0/31.D0*y
      w(20)=w(20)*w(25)
      w(20)=1771201755784684.D0/71775.D0 + w(20)
      w(20)=y*w(20)
      w(20)=17299330614802903.D0/703395.D0 + w(20)
      w(20)=w(20)*w(10)
      w(20)=184464640602122.D0/218295.D0 + w(20)
      w(20)=y*w(20)
      w(20)=341261763985622.D0/405405.D0 + w(20)
      w(20)=w(20)*w(18)
      w(20)=11654923448452.D0/125125.D0 + w(20)
      w(20)=y*w(20)
      w(20)=2677708209151.D0/28875.D0 + w(20)
      w(20)=w(20)*w(7)
      w(20)=98070404890.D0/5313.D0 + w(20)
      w(20)=y*w(20)
      w(20)=459953119922.D0/25047.D0 + w(20)
      w(26)=1.D0/23.D0*y
      w(20)=w(20)*w(26)
      w(20)=211833784124.D0/266805.D0 + w(20)
      w(20)=y*w(20)
      w(20)=95704047653.D0/121275.D0 + w(20)
      w(20)=y*w(20)
      w(20)=774143996554.D0/987525.D0 + w(20)
      w(20)=y*w(20)
      w(20)=153704944934.D0/197505.D0 + w(20)
      w(20)=w(20)*w(15)
      w(20)=7180331848.D0/176715.D0 + w(20)
      w(20)=y*w(20)
      w(20)=6326071523.D0/157080.D0 + w(20)
      w(20)=w(20)*w(16)
      w(20)=69666919.D0/29700.D0 + w(20)
      w(20)=y*w(20)
      w(20)=844022152.D0/363825.D0 + w(20)
      w(20)=y*w(20)
      w(20)=722410516.D0/315315.D0 + w(20)
      w(20)=y*w(20)
      w(20)=305244364.D0/135135.D0 + w(20)
      w(20)=w(20)*w(13)
      w(20)=19547708.D0/114345.D0 + w(20)
      w(20)=y*w(20)
      w(20)=95918024.D0/571725.D0 + w(20)
      w(20)=y*w(20)
      w(20)=2326396.D0/14175.D0 + w(20)
      w(20)=y*w(20)
      w(20)=129517.D0/810.D0 + w(20)
      w(20)=y*w(20)
      w(20)=16271.D0/105.D0 + w(20)
      w(20)=y*w(20)
      w(20)=46972.D0/315.D0 + w(20)
      w(20)=y*w(20)
      w(20)=31964.D0/225.D0 + w(20)
      w(20)=y*w(20)
      w(20)=29962.D0/225.D0 + w(20)
      w(20)=w(20)*w(19)
      w(20)=40 + w(20)
      w(20)=y*w(20)
      w(20)=244.D0/9.D0 + w(20)
      w(20)=y*w(20)
      w(20)=32.D0/3.D0 + w(20)
      w(20)=y*w(20)
      w(20)=431.D0/9.D0 + w(20)
      w(20)=y*w(20)
      w(20)= - 16 + w(20)
      w(20)=NF*w(20)
      w(2)=w(20) + 457 + w(2)
      w(20)=2*y
      w(27)=169 + 38048.D0/225.D0*y
      w(27)=w(27)*y
      w(28)=42.D0/47.D0 + 1.D0/189.D0*w(27)
      w(28)=y*w(28)
      w(28)=234584.D0/262683.D0 + w(28)
      w(28)=y*w(28)
      w(28)=74816.D0/83835.D0 + w(28)
      w(28)=y*w(28)
      w(28)=35756.D0/40095.D0 + w(28)
      w(28)=y*w(28)
      w(28)=102424.D0/114939.D0 + w(28)
      w(28)=y*w(28)
      w(28)=9304.D0/10449.D0 + w(28)
      w(28)=y*w(28)
      w(28)=8864.D0/9963.D0 + w(28)
      w(28)=y*w(28)
      w(28)=44282.D0/49815.D0 + w(28)
      w(28)=y*w(28)
      w(28)=4676.D0/5265.D0 + w(28)
      w(28)=y*w(28)
      w(28)=17752.D0/20007.D0 + w(28)
      w(28)=y*w(28)
      w(28)=151424.D0/170829.D0 + w(28)
      w(28)=y*w(28)
      w(28)=23884.D0/26973.D0 + w(28)
      w(28)=y*w(28)
      w(28)=3224.D0/3645.D0 + w(28)
      w(28)=y*w(28)
      w(28)=18248.D0/20655.D0 + w(28)
      w(28)=y*w(28)
      w(28)=40096.D0/45441.D0 + w(28)
      w(28)=w(28)*w(20)
      w(28)=4711.D0/2673.D0 + w(28)
      w(28)=y*w(28)
      w(28)=13258.D0/7533.D0 + w(28)
      w(28)=y*w(28)
      w(28)=22064.D0/12555.D0 + w(28)
      w(28)=y*w(28)
      w(28)=20608.D0/11745.D0 + w(28)
      w(28)=y*w(28)
      w(28)=12344.D0/7047.D0 + w(28)
      w(28)=y*w(28)
      w(28)=3824.D0/2187.D0 + w(28)
      w(28)=y*w(28)
      w(28)=49616.D0/28431.D0 + w(28)
      w(28)=y*w(28)
      w(28)=137536.D0/78975.D0 + w(28)
      w(28)=y*w(28)
      w(28)=10556.D0/6075.D0 + w(28)
      w(28)=y*w(28)
      w(28)=9688.D0/5589.D0 + w(28)
      w(28)=y*w(28)
      w(28)=106288.D0/61479.D0 + w(28)
      w(28)=y*w(28)
      w(28)=512.D0/297.D0 + w(28)
      w(28)=y*w(28)
      w(28)=232.D0/135.D0 + w(28)
      w(28)=y*w(28)
      w(28)=39536.D0/23085.D0 + w(28)
      w(28)=y*w(28)
      w(28)=23632.D0/13851.D0 + w(28)
      w(28)=y*w(28)
      w(28)=21056.D0/12393.D0 + w(28)
      w(28)=y*w(28)
      w(28)=6986.D0/4131.D0 + w(28)
      w(28)=y*w(28)
      w(28)=2044.D0/1215.D0 + w(28)
      w(28)=y*w(28)
      w(28)=2032.D0/1215.D0 + w(28)
      w(28)=y*w(28)
      w(28)=5248.D0/3159.D0 + w(28)
      w(28)=y*w(28)
      w(28)=1736.D0/1053.D0 + w(28)
      w(28)=y*w(28)
      w(28)=1456.D0/891.D0 + w(28)
      w(28)=y*w(28)
      w(28)=21616.D0/13365.D0 + w(28)
      w(28)=y*w(28)
      w(28)=5824.D0/3645.D0 + w(28)
      w(28)=y*w(28)
      w(28)=1148.D0/729.D0 + w(28)
      w(28)=y*w(28)
      w(28)=376.D0/243.D0 + w(28)
      w(28)=y*w(28)
      w(28)=368.D0/243.D0 + w(28)
      w(28)=y*w(28)
      w(28)=1792.D0/1215.D0 + w(28)
      w(28)=y*w(28)
      w(28)=1736.D0/1215.D0 + w(28)
      w(28)=y*w(28)
      w(28)=112.D0/81.D0 + w(28)
      w(28)=y*w(28)
      w(28)=112.D0/81.D0 + w(28)
      w(29)=y**2
      w(28)=w(28)*w(29)
      w(28)=112.D0/81.D0 + w(28)
      w(28)=NF*w(28)
      w(30)= - 845.D0/2.D0 - 19024.D0/45.D0*y
      w(30)=y*w(30)
      w(30)= - 15.D0/47.D0 + 1.D0/1323.D0*w(30)
      w(30)=y*w(30)
      w(30)= - 83780.D0/262683.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 5344.D0/16767.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 2554.D0/8019.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 36580.D0/114939.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 23260.D0/73143.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 22160.D0/69741.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 3163.D0/9963.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 334.D0/1053.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 6340.D0/20007.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 54080.D0/170829.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 8530.D0/26973.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 1612.D0/5103.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 9124.D0/28917.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 14320.D0/45441.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 3365.D0/10692.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 4735.D0/15066.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 788.D0/2511.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 736.D0/2349.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 15430.D0/49329.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 4780.D0/15309.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 8860.D0/28431.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 4912.D0/15795.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 377.D0/1215.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 1730.D0/5589.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 18980.D0/61479.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 640.D0/2079.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 58.D0/189.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 1412.D0/4617.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 4220.D0/13851.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 3760.D0/12393.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 2495.D0/8262.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 73.D0/243.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 508.D0/1701.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 6560.D0/22113.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 310.D0/1053.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 260.D0/891.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 772.D0/2673.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 208.D0/729.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 205.D0/729.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 470.D0/1701.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 460.D0/1701.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 64.D0/243.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 62.D0/243.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 20.D0/81.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 20.D0/81.D0 + w(30)
      w(30)=w(30)*w(29)
      w(30)= - 20.D0/81.D0 + w(30)
      w(30)=Ly*w(30)
      w(28)=w(30) + w(28)
      w(28)=y*w(28)
      w(30)=684300777839777498081819671.D0/792.D0 + 
     & 21697496590733107121389723.D0/25.D0*y
      w(30)=y*w(30)
      w(30)=2252951472543715527906139.D0/18612.D0 + 1.D0/7105.D0*w(30)
      w(30)=y*w(30)
      w(30)=6358982578625283554857528.D0/5172585.D0 + 1.D0/98.D0*w(30)
      w(30)=y*w(30)
      w(30)=12117781179564447210079213.D0/9904950.D0 + w(30)
      w(30)=y*w(30)
      w(30)=10668999380937281961851.D0/9474300.D0 + 1.D0/1081.D0*w(30)
      w(30)=y*w(30)
      w(30)=3887989797336440185889.D0/3470401.D0 + w(30)
      w(30)=y*w(30)
      w(30)=155050792597498532332088.D0/139131531.D0 + w(30)
      w(30)=y*w(30)
      w(30)=133668808823686417300241.D0/120600270.D0 + w(30)
      w(30)=w(30)*w(22)
      w(30)=140188642252933388381.D0/5469400.D0 + w(30)
      w(30)=y*w(30)
      w(30)=1392081947438224210811.D0/54627300.D0 + w(30)
      w(30)=w(30)*w(23)
      w(30)=6412596405631753957.D0/10379187.D0 + w(30)
      w(30)=y*w(30)
      w(30)=12092405978245563793.D0/19693842.D0 + w(30)
      w(30)=y*w(30)
      w(30)=56914129242928022987.D0/93286620.D0 + w(30)
      w(30)=y*w(30)
      w(30)=13370435889789628774.D0/22061025.D0 + w(30)
      w(30)=w(30)*w(24)
      w(30)=2485620434464062886.D0/152793025.D0 + w(30)
      w(30)=y*w(30)
      w(30)=1396148807406980191.D0/86437197.D0 + w(30)
      w(30)=y*w(30)
      w(30)=2608621248503411047.D0/162705312.D0 + w(30)
      w(30)=y*w(30)
      w(30)=2026400661157707001.D0/127370320.D0 + w(30)
      w(30)=y*w(30)
      w(30)=11307202602196765723.D0/716458050.D0 + w(30)
      w(30)=y*w(30)
      w(30)=476787165527213003.D0/30465225.D0 + w(30)
      w(30)=w(30)*w(25)
      w(30)=93915080501627227.D0/187665786.D0 + w(30)
      w(30)=w(30)*w(7)
      w(30)=262543318640294.D0/2647323.D0 + w(30)
      w(30)=y*w(30)
      w(30)=352255477287068.D0/3586275.D0 + w(30)
      w(30)=y*w(30)
      w(30)=21525407172337.D0/221375.D0 + w(30)
      w(30)=y*w(30)
      w(30)=85189925578139.D0/885500.D0 + w(30)
      w(30)=y*w(30)
      w(30)=193748342991817.D0/2036650.D0 + w(30)
      w(30)=y*w(30)
      w(30)=42121223039932.D0/448063.D0 + w(30)
      w(30)=y*w(30)
      w(30)=63294936346597.D0/681835.D0 + w(30)
      w(30)=y*w(30)
      w(30)=2468469298033.D0/26950.D0 + w(30)
      w(30)=y*w(30)
      w(30)=9907476886072.D0/109725.D0 + w(30)
      w(30)=y*w(30)
      w(30)=59132581904.D0/665.D0 + w(30)
      w(30)=y*w(30)
      w(30)=114496441858.D0/1309.D0 + w(30)
      w(30)=w(30)*w(15)
      w(30)=21526833001.D0/4760.D0 + w(30)
      w(30)=y*w(30)
      w(30)=102470505047.D0/23100.D0 + w(30)
      w(30)=w(30)*w(16)
      w(30)=10328378299.D0/40425.D0 + w(30)
      w(30)=y*w(30)
      w(30)=8746257326.D0/35035.D0 + w(30)
      w(30)=y*w(30)
      w(30)=730690007.D0/3003.D0 + w(30)
      w(30)=y*w(30)
      w(30)=3004022636.D0/12705.D0 + w(30)
      w(30)=w(30)*w(13)
      w(30)=1118591468.D0/63525.D0 + w(30)
      w(30)=y*w(30)
      w(30)=293980532.D0/17325.D0 + w(30)
      w(30)=y*w(30)
      w(30)=10239829.D0/630.D0 + w(30)
      w(30)=y*w(30)
      w(30)=756775.D0/49.D0 + w(30)
      w(30)=y*w(30)
      w(30)=3555278.D0/245.D0 + w(30)
      w(30)=y*w(30)
      w(30)=469232.D0/35.D0 + w(30)
      w(30)=y*w(30)
      w(30)=12038 + w(30)
      w(30)=w(30)*w(19)
      w(30)=16936.D0/5.D0 + w(30)
      w(30)=w(30)*w(19)
      w(30)=724 + w(30)
      w(30)=y*w(30)
      w(30)=824 + w(30)
      w(30)=w(30)*w(19)
      w(30)=32 + w(30)
      w(30)=w(30)*w(18)
      w(30)=16 + w(30)
      w(28)=1.D0/3.D0*w(30) + w(28)
      w(28)=Ly*w(28)
      w(30)= - 155971.D0/2.D0 - 1950224.D0/25.D0*y
      w(30)=w(30)*w(1)
      w(30)= - 74779.D0/47.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 1719356.D0/1081.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 182848.D0/115.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 87418.D0/55.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 751516.D0/473.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 478052.D0/301.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 455632.D0/287.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 325313.D0/205.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 103102.D0/65.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 391604.D0/247.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 1114016.D0/703.D0 + w(30)
      w(30)=w(30)*w(19)
      w(30)= - 19534.D0/37.D0 + w(30)
      w(31)=1.D0/27.D0*y
      w(30)=w(30)*w(31)
      w(30)= - 684.D0/35.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 55372.D0/2835.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 17392.D0/891.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 69523.D0/3564.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 97897.D0/5022.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 244564.D0/12555.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 228608.D0/11745.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 319786.D0/16443.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 33052.D0/1701.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 61324.D0/3159.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 510512.D0/26325.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 39227.D0/2025.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 36046.D0/1863.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 395996.D0/20493.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 120352.D0/6237.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 54622.D0/2835.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 148012.D0/7695.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 29548.D0/1539.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 1552.D0/81.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 3097.D0/162.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 7723.D0/405.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 53908.D0/2835.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 139712.D0/7371.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 19886.D0/1053.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 16756.D0/891.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 83372.D0/4455.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 2512.D0/135.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 499.D0/27.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 10394.D0/567.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 10292.D0/567.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 7264.D0/405.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 7162.D0/405.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 1412.D0/81.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 1412.D0/81.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 4 + w(30)
      w(30)=y*w(30)
      w(30)= - 244.D0/27.D0 + w(30)
      w(30)=y*w(30)
      w(30)= - 6 + w(30)
      w(30)=z2*w(30)
      w(2)=w(28) + w(30) + 1.D0/9.D0*w(2)
      w(2)=Ly*w(2)
      w(28)=1.D0/2.D0*y
      w(30)=288654698772131449369761442949248417745425444221133218512136
     & 92270927584683304285173732282313.D0/55872.D0 + 
     & 17478901909035654492101583737438200209461251024124578761369503178
     & 790467799000018278887115829.D0/33475.D0*y
      w(30)=y*w(30)
      w(30)=352596330701021355184910340138870978650958515226382360157339
     & 609314031106513123928385729557.D0/17068896.D0 + 1.D0/24745.D0*
     & w(30)
      w(30)=y*w(30)
      w(30)=758078576519826354513928196292019773534014188051272873719961
     & 7454952242683462886799272839.D0/54525640.D0 + 1.D0/147.D0*w(30)
      w(30)=w(30)*w(28)
      w(30)=416408373645649180179546537042179534724368302811909651238757
     & 34307685704446276197325831.D0/605475.D0 + w(30)
      w(4)=w(30)*w(4)
      w(4)=2854025479925843599350603195784582936675540128105268943822416
     & 79045086932426525223033.D0/9266400.D0 + w(4)
      w(4)=w(4)*w(25)
      w(4)=1450432667315256689392141839997585219140800426062050045493321
     & 326115607553765261871.D0/1475760.D0 + w(4)
      w(4)=w(4)*w(21)
      w(4)=1083954822269933349106430613735737422719261543939276054590723
     & 985991198805552412197.D0/52402896.D0 + w(4)
      w(4)=y*w(4)
      w(4)=1403831297963906646774762203243571530604580016770555887888224
     & 45472961114237571.D0/2935476180.D0 + 1.D0/427823.D0*w(4)
      w(4)=y*w(4)
      w(4)=1065976487625547617181330148717877511206419529029693070757971
     & 0581201972615921.D0/17918519361600.D0 + 1.D0/79507.D0*w(4)
      w(4)=y*w(4)
      w(4)=6467991797741183703438078840510834758724413282504923328985875
     & 61497015484522989.D0/1099366718392800.D0 + w(4)
      w(4)=y*w(4)
      w(4)=7918253747002369340607201207304481133891399640310669408760113
     & 60950497767.D0/189890614995120.D0 + 1.D0/139523.D0*w(4)
      w(4)=w(4)*w(7)
      w(4)=4539325071805632250902087267246135639652344056186210505530953
     & 096626493.D0/5504663839246.D0 + w(4)
      w(4)=y*w(4)
      w(4)=9687590039026827228127589739312995453051732298287454290625825
     & 7427925281.D0/9386900441661600.D0 + 1.D0/79.D0*w(4)
      w(4)=y*w(4)
      w(4)=4529634889131596973815384837791353388497784028375641354695494
     & 1481641551.D0/4439750208894000.D0 + w(4)
      w(30)=1.D0/1369.D0*y
      w(4)=w(4)*w(30)
      w(4)=2059215899167741181143064682643179781193296443696765342777109
     & 5764957.D0/2795398279674000.D0 + w(4)
      w(4)=y*w(4)
      w(4)=9236501649300286447705135274823968286477712512475242581839519
     & .D0/21486368764575.D0 + 1.D0/16936.D0*w(4)
      w(4)=y*w(4)
      w(4)=1603988060212181640827600971824752849781677908318247220591439
     & .D0/1876644725978880.D0 + 1.D0/497.D0*w(4)
      w(4)=y*w(4)
      w(4)=955663395838782847594996078451517769311231632767366028067139.
     & D0/1131691777507200.D0 + w(4)
      w(4)=y*w(4)
      w(4)=1656439296317794453606982633169210207309132977841272887290499
     & 1.D0/1330445145931902000.D0 + 1.D0/67.D0*w(4)
      w(4)=w(4)*w(11)
      w(4)=3105472324120745623911469289289935928271604895979485672239.D0
     & /1010235369502125.D0 + w(4)
      w(4)=w(4)*w(25)
      w(4)=6289577932118140343338280338197143223598714939580786973151.D0
     & /64237934205244800.D0 + w(4)
      w(4)=y*w(4)
      w(4)=3158722596258269977821778556287659142871626326213462391.D0/
     & 1993591061542080.D0 + 1.D0/61.D0*w(4)
      w(4)=y*w(4)
      w(4)=11905128824323963269824032457303882841220337187013.D0/
     & 638341965518400.D0 + 1.D0/83839.D0*w(4)
      w(4)=y*w(4)
      w(4)=3627912612182723108342800819652431966080975329.D0/
     & 10056184513375.D0 + 1.D0/51.D0*w(4)
      w(4)=y*w(4)
      w(4)=1839838493893804302799882049337103366537331663.D0/
     & 480928736288000.D0 + 1.D0/93.D0*w(4)
      w(4)=y*w(4)
      w(4)=7319804689326959169984134472696887098624649311.D0/
     & 102870656692003200.D0 + 1.D0/53.D0*w(4)
      w(4)=w(4)*w(10)
      w(4)=2969083418888719696584775806749883513019301.D0/
     & 1228068304695232.D0 + w(4)
      w(4)=w(4)*w(11)
      w(4)=1000384751421247306370358581319761797118113.D0/
     & 1680571942676319.D0 + w(4)
      w(4)=w(4)*w(26)
      w(4)=119818968190884783786586210406888363881423.D0/
     & 4703509734443520.D0 + w(4)
      w(4)=y*w(4)
      w(4)=15479244843430298324697452697371855765087.D0/617742061896960.
     & D0 + w(4)
      w(4)=y*w(4)
      w(4)=940185001451024846399649876042205666631.D0/38172765286656.D0
     &  + w(4)
      w(4)=w(4)*w(23)
      w(4)=339060830167831504869831721991010661.D0/574716481416.D0 + 
     & w(4)
      w(4)=w(4)*w(15)
      w(4)=1774732378372008931081397108952037.D0/58255861664.D0 + w(4)
      w(4)=w(4)*w(24)
      w(4)=99543795534584390576918115430819.D0/123365354112.D0 + w(4)
      w(4)=w(4)*w(16)
      w(4)=589938647749761430623678997567.D0/12699374688.D0 + w(4)
      w(4)=y*w(4)
      w(4)=2775604871251890749575081987.D0/4585885304.D0 + 1.D0/75.D0*
     & w(4)
      w(4)=y*w(4)
      w(4)=2694692699010663335980064821.D0/4564752192.D0 + w(4)
      w(4)=w(4)*w(7)
      w(4)=1093036257835777437723649.D0/9513504.D0 + w(4)
      w(4)=w(4)*w(13)
      w(4)=5491780001671811734980869.D0/640332000.D0 + w(4)
      w(4)=y*w(4)
      w(4)=45260480327056771745893.D0/5457375.D0 + w(4)
      w(4)=w(4)*w(14)
      w(4)=20038532292125308729.D0/635040.D0 + w(4)
      w(4)=y*w(4)
      w(4)=2071092336527802971.D0/68600.D0 + w(4)
      w(4)=y*w(4)
      w(4)=60265417202191.D0/17150.D0 + 1.D0/8151.D0*w(4)
      w(3)=w(4)*w(3)
      w(3)=16963368134177.D0/175175.D0 + w(3)
      w(3)=y*w(3)
      w(3)=7624596947887.D0/85800.D0 + w(3)
      w(3)=w(3)*w(18)
      w(3)=9506782273.D0/1100.D0 + w(3)
      w(3)=w(3)*w(28)
      w(3)=42780814.D0/15.D0 + w(3)
      w(3)=y*w(3)
      w(3)=23005457.D0/5.D0 + w(3)
      w(3)=y*w(3)
      w(3)= - 250255.D0/6.D0 + 1.D0/35.D0*w(3)
      w(3)=w(3)*w(18)
      w(3)=5491 + w(3)
      w(4)=1.D0/2401.D0*y
      w(14)= - 702568809154139649781950971.D0/9792.D0 - 
     & 1870492897534859699234258.D0/25.D0*y
      w(14)=w(14)*w(4)
      w(14)= - 85669797096570471207805667.D0/2991456.D0 + w(14)
      w(14)=w(14)*w(7)
      w(14)= - 220159914706694147031679831.D0/40135368.D0 + w(14)
      w(14)=w(14)*w(21)
      w(14)= - 2682782421707701124094323.D0/24017175.D0 + w(14)
      w(14)=w(14)*w(15)
      w(14)= - 29526564515745257441347.D0/5250960.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 2988945210321605976305407.D0/556115560.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 1257817030489719999048481.D0/245001960.D0 + w(14)
      w(14)=w(14)*w(22)
      w(14)= - 12361311522099001267151.D0/108460170.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 14974558589407009987261.D0/137727200.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 112189307020197976183.D0/1375592400.D0 + 1.D0/1271.D0*
     & w(14)
      w(14)=y*w(14)
      w(14)= - 6296263722753165209599.D0/81022392360.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 37423914517654220921.D0/505706565.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 5129300150036318266349.D0/72821745360.D0 + w(14)
      w(14)=w(14)*w(24)
      w(14)= - 62378359759061094803.D0/34442717400.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 5335598617641573139.D0/3098022200.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 7173561056769624523.D0/4381488540.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 44916789044118231161.D0/28866277440.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 6682599115728890099.D0/4519467680.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 8929697473883264059.D0/6355501425.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 512178098422781059.D0/383578650.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 2726539489217356643.D0/2148040440.D0 + w(14)
      w(14)=w(14)*w(10)
      w(14)= - 4625046581169533.D0/111105540.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 89831842345811.D0/2267460.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 257518666128536.D0/6823375.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 2116073560683947.D0/58786000.D0 + w(14)
      w(14)=w(14)*w(7)
      w(14)= - 37192663731767.D0/5408312.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 978824582943479.D0/148728580.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 640225192.D0/20995.D0 + 1.D0/207.D0*w(14)
      w(14)=y*w(14)
      w(14)= - 1206679472737.D0/41150200.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 6092317871473.D0/215408700.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 3537533050523.D0/129245220.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 6663339581.D0/250614.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 5777347585.D0/222768.D0 + w(14)
      w(14)=w(14)*w(8)
      w(14)= - 162155147.D0/70200.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 2150754106.D0/945945.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 9217491677.D0/4099095.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 2243050751.D0/1003860.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 509984633.D0/228690.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 2554296529.D0/1143450.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 392549.D0/175.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 2839229.D0/1260.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 4977863.D0/2205.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 990488.D0/441.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 496826.D0/225.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 481538.D0/225.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 18653.D0/9.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 19616.D0/9.D0 + w(14)
      w(14)=y*w(14)
      w(14)= - 7402 + w(14)
      w(14)=y*w(14)
      w(14)=9914 + w(14)
      w(14)=w(14)*w(18)
      w(14)= - 220 + w(14)
      w(14)=z2*w(14)
      w(3)=1.D0/3.D0*w(3) + w(14)
      w(14)=12907102468180463043077061912810659136255190657.D0/208.D0
     &  + 1552725207771221790851536004571329435319664429.D0/25.D0*y
      w(14)=w(14)*w(4)
      w(14)=126211796466762296453573000661050568216627877.D0/4888.D0 + 
     & w(14)
      w(14)=y*w(14)
      w(14)=55770387553221939531745803269851567786649213.D0/2162.D0 + 
     & w(14)
      w(14)=y*w(14)
      w(14)=498295949353204049050193963650064269397.D0/299.D0 + 1.D0/
     & 15463.D0*w(14)
      w(14)=y*w(14)
      w(14)=952260545326157849445026575396823848211.D0/572.D0 + w(14)
      w(14)=y*w(14)
      w(14)=143158238015787909953175515065027311726373.D0/86086.D0 + 
     & w(14)
      w(14)=y*w(14)
      w(14)=1857072453275497386980670684309975548999.D0/1118.D0 + w(14)
      w(5)=w(14)*w(5)
      w(5)=23434431626836946882571747526928875349.D0/26117.D0 + w(5)
      w(5)=y*w(5)
      w(5)=187241304292500035170402759540390200211.D0/208936.D0 + w(5)
      w(5)=w(5)*w(6)
      w(5)=17635803833733218038019998328466781.D0/33124.D0 + w(5)
      w(5)=y*w(5)
      w(5)=167310583299870868492413389571741251.D0/314678.D0 + w(5)
      w(5)=y*w(5)
      w(5)=237753573712582977125037618479345437.D0/447811.D0 + w(5)
      w(5)=y*w(5)
      w(5)=49977807660804353528095767361044677.D0/94276.D0 + w(5)
      w(5)=w(5)*w(30)
      w(5)=3447861960271801951667072444977.D0/8918.D0 + w(5)
      w(5)=y*w(5)
      w(5)=58515633379353974459984614995601.D0/151606.D0 + w(5)
      w(5)=y*w(5)
      w(5)=8344634059396783973023772629129.D0/21658.D0 + w(5)
      w(5)=y*w(5)
      w(5)=31356582601784180986456750553393.D0/81536.D0 + w(5)
      w(5)=y*w(5)
      w(5)=485069537916740290725642460913329.D0/1263808.D0 + w(5)
      w(5)=w(5)*w(11)
      w(5)=1890848009290310335547008832458.D0/19747.D0 + w(5)
      w(5)=w(5)*w(12)
      w(5)=1836565444550902893887621234.D0/18473.D0 + w(5)
      w(5)=y*w(5)
      w(5)=3664479051100838573715978977.D0/36946.D0 + w(5)
      w(6)=1.D0/841.D0*y
      w(5)=w(5)*w(6)
      w(5)=1348868279425261324360187.D0/11466.D0 + w(5)
      w(5)=y*w(5)
      w(5)=8744153595815723270282623.D0/74529.D0 + w(5)
      w(5)=y*w(5)
      w(5)=244643251713075923741317.D0/207025.D0 + 1.D0/99.D0*w(5)
      w(5)=w(5)*w(19)
      w(5)=50028915796528837661827.D0/127400.D0 + w(5)
      w(14)=1.D0/125.D0*y
      w(5)=w(5)*w(14)
      w(5)=45874396567508589781.D0/14651.D0 + w(5)
      w(5)=y*w(5)
      w(5)=33186039540721033724591.D0/10636626.D0 + w(5)
      w(5)=y*w(5)
      w(5)=47547910416673705328.D0/8093085.D0 + 1.D0/529.D0*w(5)
      w(5)=y*w(5)
      w(5)=86088840008782231997.D0/14714700.D0 + w(5)
      w(5)=y*w(5)
      w(5)=58151084176569625823.D0/9984975.D0 + w(5)
      w(5)=y*w(5)
      w(5)=7714826729187284327.D0/1331330.D0 + w(5)
      w(5)=y*w(5)
      w(5)=9508247676340058.D0/595595.D0 + 1.D0/361.D0*w(5)
      w(5)=y*w(5)
      w(5)=1814558317550527457.D0/114354240.D0 + w(5)
      w(32)=1.D0/289.D0*y
      w(5)=w(5)*w(32)
      w(5)=917160511227463.D0/16816800.D0 + w(5)
      w(5)=y*w(5)
      w(5)=199129509537089.D0/3678675.D0 + w(5)
      w(5)=y*w(5)
      w(5)=513390970939807.D0/9564555.D0 + w(5)
      w(5)=y*w(5)
      w(5)=145289654855327.D0/2732730.D0 + w(5)
      w(5)=w(5)*w(13)
      w(5)=719549192371.D0/177870.D0 + w(5)
      w(5)=y*w(5)
      w(5)=1776355687853.D0/444675.D0 + w(5)
      w(5)=w(5)*w(8)
      w(5)=3945231242.D0/11025.D0 + w(5)
      w(5)=y*w(5)
      w(5)=886096523.D0/2520.D0 + w(5)
      w(5)=w(5)*w(19)
      w(5)=56220827.D0/490.D0 + w(5)
      w(5)=w(5)*w(18)
      w(5)=3040889.D0/245.D0 + w(5)
      w(5)=y*w(5)
      w(5)=299509.D0/25.D0 + w(5)
      w(5)=y*w(5)
      w(5)=569059.D0/50.D0 + w(5)
      w(5)=w(5)*w(7)
      w(5)=12493.D0/6.D0 + w(5)
      w(5)=y*w(5)
      w(5)=1723 + w(5)
      w(5)=y*w(5)
      w(5)=1361 + w(5)
      w(5)=y*w(5)
      w(5)=5432.D0/3.D0 + w(5)
      w(5)=w(5)*w(31)
      w(5)= - 32 + w(5)
      w(27)= - 18.D0/47.D0 - 1.D0/441.D0*w(27)
      w(27)=y*w(27)
      w(27)= - 33512.D0/87561.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 10688.D0/27945.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 5108.D0/13365.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 14632.D0/38313.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 9304.D0/24381.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 8864.D0/23247.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 6326.D0/16605.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 668.D0/1755.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 2536.D0/6669.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 21632.D0/56943.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 3412.D0/8991.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 3224.D0/8505.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 18248.D0/48195.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 5728.D0/15147.D0 + w(27)
      w(27)=w(27)*w(20)
      w(27)= - 673.D0/891.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 1894.D0/2511.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 3152.D0/4185.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 2944.D0/3915.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 12344.D0/16443.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 3824.D0/5103.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 7088.D0/9477.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 19648.D0/26325.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 1508.D0/2025.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 1384.D0/1863.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 15184.D0/20493.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 512.D0/693.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 232.D0/315.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 5648.D0/7695.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 3376.D0/4617.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 3008.D0/4131.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 998.D0/1377.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 292.D0/405.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 2032.D0/2835.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 5248.D0/7371.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 248.D0/351.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 208.D0/297.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 3088.D0/4455.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 832.D0/1215.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 164.D0/243.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 376.D0/567.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 368.D0/567.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 256.D0/405.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 248.D0/405.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 16.D0/27.D0 + w(27)
      w(27)=y*w(27)
      w(27)= - 16.D0/27.D0 + w(27)
      w(27)=w(27)*w(29)
      w(27)= - 16.D0/27.D0 + w(27)
      w(33)=z2*y
      w(27)=w(27)*w(33)
      w(5)=1.D0/3.D0*w(5) + w(27)
      w(5)=NF*w(5)
      w(27)=767309 + 3837488.D0/5.D0*y
      w(27)=w(27)*w(1)
      w(27)=735802.D0/47.D0 + w(27)
      w(27)=y*w(27)
      w(27)=16918936.D0/1081.D0 + w(27)
      w(27)=y*w(27)
      w(27)=1079632.D0/69.D0 + w(27)
      w(27)=y*w(27)
      w(27)=516196.D0/33.D0 + w(27)
      w(27)=y*w(27)
      w(27)=7396568.D0/473.D0 + w(27)
      w(27)=w(27)*w(18)
      w(27)=522824.D0/301.D0 + w(27)
      w(27)=y*w(27)
      w(27)=36464.D0/21.D0 + w(27)
      w(27)=y*w(27)
      w(27)=15622.D0/9.D0 + w(27)
      w(27)=y*w(27)
      w(27)=203012.D0/117.D0 + w(27)
      w(27)=y*w(27)
      w(27)=3855752.D0/2223.D0 + w(27)
      w(27)=y*w(27)
      w(27)=10969648.D0/6327.D0 + w(27)
      w(27)=y*w(27)
      w(27)=1731316.D0/999.D0 + w(27)
      w(27)=y*w(27)
      w(27)=327400.D0/189.D0 + w(27)
      w(27)=y*w(27)
      w(27)=1854392.D0/1071.D0 + w(27)
      w(27)=y*w(27)
      w(27)=970864.D0/561.D0 + w(27)
      w(27)=y*w(27)
      w(27)=38053.D0/22.D0 + w(27)
      w(27)=y*w(27)
      w(27)=482311.D0/279.D0 + w(27)
      w(27)=y*w(27)
      w(27)=482024.D0/279.D0 + w(27)
      w(27)=y*w(27)
      w(27)=450640.D0/261.D0 + w(27)
      w(27)=y*w(27)
      w(27)=3152348.D0/1827.D0 + w(27)
      w(27)=y*w(27)
      w(27)=977608.D0/567.D0 + w(27)
      w(27)=y*w(27)
      w(27)=1814152.D0/1053.D0 + w(27)
      w(27)=y*w(27)
      w(27)=1007024.D0/585.D0 + w(27)
      w(27)=y*w(27)
      w(27)=25798.D0/15.D0 + w(27)
      w(27)=y*w(27)
      w(27)=118556.D0/69.D0 + w(27)
      w(27)=y*w(27)
      w(27)=3908248.D0/2277.D0 + w(27)
      w(27)=y*w(27)
      w(27)=1188112.D0/693.D0 + w(27)
      w(27)=y*w(27)
      w(27)=107876.D0/63.D0 + w(27)
      w(27)=y*w(27)
      w(27)=292408.D0/171.D0 + w(27)
      w(27)=y*w(27)
      w(27)=875912.D0/513.D0 + w(27)
      w(27)=y*w(27)
      w(27)=782416.D0/459.D0 + w(27)
      w(27)=y*w(27)
      w(27)=260327.D0/153.D0 + w(27)
      w(27)=w(27)*w(19)
      w(27)=566 + w(27)
      w(27)=y*w(27)
      w(27)=35576.D0/63.D0 + w(27)
      w(27)=y*w(27)
      w(27)=1383856.D0/2457.D0 + w(27)
      w(27)=y*w(27)
      w(27)=197108.D0/351.D0 + w(27)
      w(27)=y*w(27)
      w(27)=166216.D0/297.D0 + w(27)
      w(27)=y*w(27)
      w(27)=165560.D0/297.D0 + w(27)
      w(27)=y*w(27)
      w(27)=44944.D0/81.D0 + w(27)
      w(27)=y*w(27)
      w(27)=44698.D0/81.D0 + w(27)
      w(27)=y*w(27)
      w(27)=103612.D0/189.D0 + w(27)
      w(27)=y*w(27)
      w(27)=34264.D0/63.D0 + w(27)
      w(27)=y*w(27)
      w(27)=1616.D0/3.D0 + w(27)
      w(27)=y*w(27)
      w(27)=14380.D0/27.D0 + w(27)
      w(27)=y*w(27)
      w(27)=14216.D0/27.D0 + w(27)
      w(27)=y*w(27)
      w(27)=14216.D0/27.D0 + w(27)
      w(27)=w(27)*w(19)
      w(27)=48 + w(27)
      w(27)=y*w(27)
      w(27)=680.D0/9.D0 + w(27)
      w(27)=y*w(27)
      w(27)=72 + w(27)
      w(27)=z3*w(27)
      w(2)=w(2) + w(5) + 1.D0/3.D0*w(3) + w(27)
      w(2)=Ly*w(2)
      w(3)=2215691131739904442293512073447278280148734322179140350182439
     & 6048001.D0/4.D0 + 13853398047536771761444646963087235642133800054
     & 4478341844159434795247.D0/25.D0*y
      w(3)=w(3)*w(4)
      w(3)=4335510021608511863112895894150590335110627394939677911396462
     & 78247.D0/188.D0 + w(3)
      w(3)=y*w(3)
      w(3)=4983737595802723571524484824609050706243421733126619761331736
     & 534953.D0/2162.D0 + w(3)
      w(3)=y*w(3)
      w(3)=35488966378276435227167997454577309418919143772923301251133.D
     & 0/23.D0 + 1.D0/1493284.D0*w(3)
      w(3)=w(3)*w(21)
      w(3)=976036827284465864089030032111090748939785659315874837823487.
     & D0/29744.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1613431654708476083140653026552444758148508856410814839594751
     & .D0/49192.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1905811252681104011211915729188682056778977951890260582929933
     & .D0/58136.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1096270258905545699853216537596327794496093946372913499.D0/
     & 6929.D0 + 1.D0/207088.D0*w(3)
      w(3)=w(3)*w(22)
      w(3)=11415166052753951646565526552220351824903403234122224049.D0/
     & 3104192.D0 + w(3)
      w(3)=y*w(3)
      w(3)=226512569844029634229461170785944156079993545424781.D0/
     & 1968512.D0 + 1.D0/31939.D0*w(3)
      w(3)=w(3)*w(23)
      w(3)=996582807641834629821877447264999684327128502951853.D0/
     & 355316416.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1417275288850692354770330972028985529187159509267161.D0/
     & 505642592.D0 + w(3)
      w(3)=y*w(3)
      w(3)=74541226467114657125909012263437289950345388297069.D0/
     & 26612768.D0 + w(3)
      w(3)=w(3)*w(30)
      w(3)=20587206087303959412191140961580479001962930057.D0/10069696.D
     & 0 + w(3)
      w(3)=y*w(3)
      w(3)=4882000649187605693513566267606354506859087.D0/10699052.D0
     &  + 1.D0/4477.D0*w(3)
      w(3)=y*w(3)
      w(3)=5396380533687899130759710032768384983806184919.D0/
     & 11836208384.D0 + w(3)
      w(3)=y*w(3)
      w(3)=20297754330866304650378082453305640641415756923.D0/
     & 44559843328.D0 + w(3)
      w(3)=y*w(3)
      w(3)=28574337851752547227110255726869797981192898129.D0/
     & 62788870144.D0 + w(3)
      w(3)=y*w(3)
      w(3)=55752939491052661631845218981855971438498827.D0/122634512.D0
     &  + w(3)
      w(3)=w(3)*w(12)
      w(3)=108428112201399151295248855102576602020161.D0/229445216.D0
     &  + w(3)
      w(3)=w(3)*w(25)
      w(3)=998182477680717935968366321553128897201.D0/65555776.D0 + 
     & w(3)
      w(3)=w(3)*w(6)
      w(3)=1287618885144161364381598877983649171.D0/71207136.D0 + w(3)
      w(3)=w(3)*w(10)
      w(3)=144109486206952884660288440470785373.D0/231423192.D0 + w(3)
      w(3)=w(3)*w(31)
      w(3)=59218933985769450217964152124827891.D0/2571368800.D0 + w(3)
      w(3)=y*w(3)
      w(3)=4548211187293697959837094557329811.D0/197797600.D0 + w(3)
      w(3)=w(3)*w(14)
      w(3)=735203287065966391769841811780453.D0/4003423424.D0 + w(3)
      w(3)=y*w(3)
      w(3)=4036195561125789791130255148530679.D0/22018828832.D0 + w(3)
      w(3)=y*w(3)
      w(3)=92699228686439670655685823011.D0/6701382688.D0 + 1.D0/13225.D
     & 0*w(3)
      w(3)=w(3)*w(26)
      w(3)=9139920085183384297444276049.D0/15230415200.D0 + w(3)
      w(3)=y*w(3)
      w(3)=12374291888200930009105874707.D0/20669849200.D0 + w(3)
      w(3)=y*w(3)
      w(3)=24682789879500986964844837957.D0/41339698400.D0 + w(3)
      w(3)=w(3)*w(15)
      w(3)=30497958982907137389610933.D0/973372400.D0 + w(3)
      w(3)=w(3)*w(15)
      w(3)=102392801832094119166108307.D0/62295833600.D0 + w(3)
      w(3)=w(3)*w(32)
      w(3)=51911492960211043283317.D0/9161152000.D0 + w(3)
      w(3)=w(3)*w(16)
      w(3)=41575816670314021153.D0/125250125.D0 + w(3)
      w(3)=y*w(3)
      w(3)=430358501083046962103.D0/1302601300.D0 + w(3)
      w(3)=y*w(3)
      w(3)=244600516163831612051.D0/744343600.D0 + w(3)
      w(3)=y*w(3)
      w(3)=1217112949223714371.D0/3726800.D0 + w(3)
      w(3)=w(3)*w(13)
      w(3)=58095629764812401.D0/2329250.D0 + w(3)
      w(3)=y*w(3)
      w(3)=259688927808703.D0/10500.D0 + w(3)
      w(3)=w(3)*w(8)
      w(3)=10685065576273.D0/4800.D0 + w(3)
      w(3)=w(3)*w(31)
      w(3)=57029452841.D0/700.D0 + w(3)
      w(3)=y*w(3)
      w(3)=56223949651.D0/700.D0 + w(3)
      w(3)=w(3)*w(9)
      w(3)=5642285771.D0/500.D0 + w(3)
      w(3)=w(3)*w(9)
      w(3)=197991541.D0/125.D0 + w(3)
      w(5)=1.D0/15.D0*y
      w(3)=w(3)*w(5)
      w(3)=837895.D0/8.D0 + w(3)
      w(3)=y*w(3)
      w(3)=112075 + w(3)
      w(3)=w(3)*w(19)
      w(3)= - 11629.D0/2.D0 + w(3)
      w(3)=y*w(3)
      w(3)=213802.D0/3.D0 + w(3)
      w(3)=y*w(3)
      w(3)= - 14968 + w(3)
      w(6)=33563444462114597612535743.D0/144.D0 + 
     & 5826338148943620434839528.D0/25.D0*y
      w(4)=w(6)*w(4)
      w(4)=328543758439677196449863.D0/3384.D0 + w(4)
      w(4)=w(4)*w(7)
      w(4)=125957383962699316123427.D0/6486.D0 + w(4)
      w(4)=w(4)*w(21)
      w(4)=6415600350281629016306.D0/15525.D0 + w(4)
      w(4)=y*w(4)
      w(4)=12275021801198158377839.D0/29700.D0 + w(4)
      w(4)=y*w(4)
      w(4)=13686381680619761446037.D0/33110.D0 + w(4)
      w(4)=w(4)*w(8)
      w(4)=1018153404961231817903.D0/27090.D0 + w(4)
      w(4)=w(4)*w(22)
      w(4)=2258040275068920856.D0/2583.D0 + w(4)
      w(4)=y*w(4)
      w(4)=50187447779032339871.D0/57400.D0 + w(4)
      w(4)=w(4)*w(23)
      w(4)=1746883692744024953.D0/81900.D0 + w(4)
      w(4)=y*w(4)
      w(4)=3319735717846210979.D0/155610.D0 + w(4)
      w(4)=y*w(4)
      w(4)=1575078466910651306.D0/73815.D0 + w(4)
      w(4)=y*w(4)
      w(4)=597006135017010041.D0/27972.D0 + w(4)
      w(4)=w(4)*w(24)
      w(4)=38166980583699547.D0/66150.D0 + w(4)
      w(4)=y*w(4)
      w(4)=264409695283011611.D0/458150.D0 + w(4)
      w(4)=y*w(4)
      w(4)=34004907499941566.D0/58905.D0 + w(4)
      w(4)=y*w(4)
      w(4)=64028188772524181.D0/110880.D0 + w(4)
      w(4)=y*w(4)
      w(4)=22061175937059593.D0/38192.D0 + w(4)
      w(4)=y*w(4)
      w(4)=310341709172904247.D0/537075.D0 + w(4)
      w(4)=w(4)*w(25)
      w(4)=9368618314788584.D0/502425.D0 + w(4)
      w(4)=y*w(4)
      w(4)=26242638227406553.D0/1406790.D0 + w(4)
      w(4)=w(4)*w(10)
      w(4)=140479364995261.D0/218295.D0 + w(4)
      w(4)=y*w(4)
      w(4)=261013323575629.D0/405405.D0 + w(4)
      w(4)=w(4)*w(18)
      w(4)=8955673379564.D0/125125.D0 + w(4)
      w(4)=y*w(4)
      w(4)=8271413853653.D0/115500.D0 + w(4)
      w(4)=w(4)*w(7)
      w(4)=152287527577.D0/10626.D0 + w(4)
      w(4)=y*w(4)
      w(4)=2514439918129.D0/175329.D0 + w(4)
      w(4)=w(4)*w(26)
      w(4)=166485510328.D0/266805.D0 + w(4)
      w(4)=y*w(4)
      w(4)=30294977527.D0/48510.D0 + w(4)
      w(4)=y*w(4)
      w(4)=123456357217.D0/197505.D0 + w(4)
      w(4)=y*w(4)
      w(4)=123581825137.D0/197505.D0 + w(4)
      w(4)=w(4)*w(15)
      w(4)=5826259784.D0/176715.D0 + w(4)
      w(4)=y*w(4)
      w(4)=5185515721.D0/157080.D0 + w(4)
      w(4)=w(4)*w(16)
      w(4)=404299307.D0/207900.D0 + w(4)
      w(4)=y*w(4)
      w(4)=708683198.D0/363825.D0 + w(4)
      w(4)=y*w(4)
      w(4)=615341228.D0/315315.D0 + w(4)
      w(4)=y*w(4)
      w(4)=264284327.D0/135135.D0 + w(4)
      w(4)=w(4)*w(13)
      w(4)=17244646.D0/114345.D0 + w(4)
      w(4)=y*w(4)
      w(4)=86472046.D0/571725.D0 + w(4)
      w(4)=y*w(4)
      w(4)=2151164.D0/14175.D0 + w(4)
      w(4)=y*w(4)
      w(4)=863867.D0/5670.D0 + w(4)
      w(4)=y*w(4)
      w(4)=112501.D0/735.D0 + w(4)
      w(4)=y*w(4)
      w(4)=339308.D0/2205.D0 + w(4)
      w(4)=y*w(4)
      w(4)=34828.D0/225.D0 + w(4)
      w(4)=y*w(4)
      w(4)=35024.D0/225.D0 + w(4)
      w(4)=w(4)*w(19)
      w(4)=52 + w(4)
      w(4)=y*w(4)
      w(4)=464.D0/9.D0 + w(4)
      w(4)=y*w(4)
      w(4)=16.D0/3.D0 + w(4)
      w(4)=y*w(4)
      w(4)=958.D0/9.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 32 + w(4)
      w(4)=z2*w(4)
      w(3)=1.D0/27.D0*w(3) + w(4)
      w(4)= - 1531 - 68918.D0/45.D0*y
      w(4)=w(4)*w(1)
      w(4)= - 1468.D0/47.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 303766.D0/9729.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 19382.D0/621.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 9266.D0/297.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 132758.D0/4257.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 84446.D0/2709.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 80482.D0/2583.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 11492.D0/369.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 1214.D0/39.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 23054.D0/741.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 196738.D0/6327.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 31046.D0/999.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 5870.D0/189.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 33242.D0/1071.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 52202.D0/1683.D0 + w(4)
      w(4)=w(4)*w(20)
      w(4)= - 6137.D0/99.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 17282.D0/279.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 5756.D0/93.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 5380.D0/87.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 112876.D0/1827.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 34996.D0/567.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 64924.D0/1053.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 36028.D0/585.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 2768.D0/45.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 12716.D0/207.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 139676.D0/2277.D0 + w(4)
      w(4)=w(4)*w(18)
      w(4)= - 524.D0/77.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 428.D0/63.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 10436.D0/1539.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 31244.D0/4617.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 27892.D0/4131.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 9274.D0/1377.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 544.D0/81.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 3796.D0/567.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 49172.D0/7371.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 2332.D0/351.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 1964.D0/297.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 5860.D0/891.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 1588.D0/243.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 1576.D0/243.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 3644.D0/567.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 3604.D0/567.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 508.D0/81.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 500.D0/81.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 164.D0/27.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 164.D0/27.D0 + w(4)
      w(4)=w(4)*w(29)
      w(4)= - 248.D0/27.D0 + w(4)
      w(4)=y*w(4)
      w(4)=28.D0/9.D0 + w(4)
      w(4)=z3*w(4)
      w(3)=1.D0/9.D0*w(3) + 4*w(4)
      w(3)=NF*w(3)
      w(4)= - 1575299615668941392345610152553424542010442311491772284851
     & 15093.D0/11.D0 - 145957148741011690636520568353535856750106190297
     & 2038217611371091.D0/100.D0*y
      w(4)=w(4)*w(1)
      w(4)= - 1481861389873930570536247360992689556756002171257582943168
     & 10799.D0/517.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 4261386743959471180741277943749344886406434278378455245681
     & 267.D0/11891.D0 + 1.D0/784.D0*w(4)
      w(4)=w(4)*w(9)
      w(4)= - 6345070869122540279894303770349046979185535948058884941448
     & 3.D0/1265.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 5885253443787478666368962717254003552014456356191990363.D0
     & /4235.D0 + 1.D0/35344.D0*w(4)
      w(4)=w(4)*w(25)
      w(4)= - 1597876100648280260518632516434364048079144877283919603.D0
     & /36421.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 1989152847496956784123778851472831228581593447858645073.D0
     & /46354.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 94998373895354283145876509684182438998077132323358073.D0/
     & 124558.D0 + 1.D0/55.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 1084492801164571580855351934236208021651014030643707.D0/
     & 5382685.D0 + 1.D0/3698.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 1678480255559551262223129193688890603573933699451207.D0/
     & 8533525.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 27217512068050274302743504149951224765586430439.D0/
     & 32427395.D0 + 1.D0/228616.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 2121947743160469006907578830898799660993007683.D0/2593367.
     & D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 65816170184643208807940051925161301304734308569.D0/
     & 82578265.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 60595184983308786600398535641123068182341066279.D0/
     & 78114575.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 1146897897851858453694548401781600433377551.D0/120722525.D
     & 0 + 1.D0/79402.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 111712662886239169191750167320413871799682481.D0/
     & 12103295435.D0 + w(4)
      w(4)=w(4)*w(11)
      w(4)= - 1594527784960752329276723910414752740814309.D0/711958555.D
     & 0 + w(4)
      w(4)=y*w(4)
      w(4)= - 4357019503090431177490399635002652070597933.D0/2006428655.
     & D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 659354085925838165282496538954494981886573.D0/10032143275.
     & D0 + 1.D0/32.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 17048452315865812631004753574764764624449.D0/268140235.D0
     &  + w(4)
      w(4)=y*w(4)
      w(4)= - 104983063039530512642750394313465225781.D0/423834565.D0
     &  + 1.D0/248.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 10478729327359949061848920963541647067.D0/43844955.D0 + 
     & w(4)
      w(4)=w(4)*w(10)
      w(4)= - 22294372299364807044973568537612603.D0/2807805.D0 + w(4)
      w(4)=w(4)*w(18)
      w(4)= - 2836831726569044216823489297727343.D0/3342625.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 181284693781234228283709760763.D0/15125.D0 + 1.D0/68.D0*
     & w(4)
      w(4)=y*w(4)
      w(4)= - 2717905123934904125559208069817.D0/236555.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 17018413284993958890991054818781.D0/61930099.D0 + 1.D0/40.
     & D0*w(4)
      w(4)=y*w(4)
      w(4)= - 4944526741969892732766650295377.D0/18848291.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 2023124162797743555707553937.D0/8567405.D0 + 1.D0/1058.D0
     & *w(4)
      w(4)=w(4)*w(5)
      w(4)= - 584444830699646630231809.D0/39083.D0 + w(4)
      w(4)=w(4)*w(28)
      w(4)= - 211736830287217422256186037.D0/29898495.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 14729910398325761959549093.D0/2203047.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 175756787758656042509.D0/524535.D0 + 1.D0/18772.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 14985484216659198111059.D0/47732685.D0 + w(4)
      w(4)=w(4)*w(17)
      w(4)= - 105790278285229967977.D0/5780775.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 255188042049943461317.D0/15030015.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 131316432084055973.D0/2147145.D0 + 1.D0/256.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 203468326957823897.D0/3633630.D0 + w(4)
      w(4)=w(4)*w(13)
      w(4)= - 12491769169842367.D0/3194400.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 2381279539183289.D0/677600.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 873418016237.D0/280.D0 + w(4)
      w(4)=w(4)*w(18)
      w(4)= - 593907954551.D0/1960.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 191136301547.D0/735.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 114599958421.D0/525.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 478627393.D0/75.D0 + 1.D0/28.D0*w(4)
      w(4)=y*w(4)
      w(4)= - 76246367.D0/15.D0 + w(4)
      w(4)=y*w(4)
      w(4)= - 1425737.D0/3.D0 + 1.D0/8.D0*w(4)
      w(4)=w(4)*w(18)
      w(4)=74771 + w(4)
      w(4)=w(4)*w(7)
      w(4)= - 60652.D0/3.D0 + w(4)
      w(4)=w(4)*w(18)
      w(4)=784 + w(4)
      w(4)=z2*w(4)
      w(4)=w(4) - 2714*z4 + 2299927.D0/810.D0 + 32*B4
      w(5)=16*B4
      w(6)=12020998.D0/11025.D0*z4 + 74072100129305179617786661155874130
     & 51277907501158257643992511909279372257263592153539789685446476280
     & 8943019415328032747625307987849917.D0/
     & 26143989793959931700128194585304827774114430409973452236172707179
     & 46941955703103422442616611671972236005610298007394116240000000000
     & .D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 5768941.D0/5292.D0*z4 + 
     & 97119271466090377912792184087497647469775374667562214539309836956
     & 74745158092964445641291127766704483268867479814257856001295222443
     & 20167.D0/35002990028314322717992050483833498391936568437924739194
     & 19195590510834256987381799834703729003108033828733352277162533437
     & 4400000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 2766169.D0/2538.D0*z4 + 
     & 11427718408176752390227447046819120713749264406275208910326080909
     & 2982400807338230065648230792139212802737970175342795533270848757.
     & D0/42087915955497669445677105297569032712866976788178335282498842
     & 20891722345631934414326913819102424625074380328767897600000000.D0
     &  - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 10601380.D0/9729.D0*z4 + 
     & 10254342224770114823660587485677272228790364437100585049123991701
     & 05963774213681583652941487152918507111516663607858113559555284673
     & .D0/3862210491712148380867083754500390930477680839347834379546041
     & 7141305366381906108355144832627742963931932553363887452160000000.
     & D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 10147994.D0/9315.D0*z4 + 
     & 21240151467178731407848396137440669688275756224933053861414075349
     & 273428028322956918225111628781951435159677488120848589576861.D0/
     & 81877748013561915754334148604756344589106079667692425107342090008
     & 8276299130460798803958197301113435669035051968000000000.D0 - 
     & w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 4852259.D0/4455.D0*z4 + 718285812165952533809265780998
     & 54489480178997520002764131875753136731328992067178690026666352618
     & 211663206798489778108805809.D0/2836318548701339810018857668123903
     & 52661693664682304934731281277019970926539573485129011047149006957
     & 0684161024000000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 4635476.D0/4257.D0*z4 + 100416399034874751522981771723
     & 58262438150122785411175958700343979631564783068279390536474307652
     & 7501643058130424302990012757.D0/406538991980525372769369599097759
     & 50548176091937797040644816983039529166137338866201824916758024330
     & 51313964134400000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 8847296.D0/8127.D0*z4 + 441608642044321915046878411231
     & 49944924058698717843559356640593880942771061973426101887072329630
     & 1197936708828364279287691.D0/183479880326478246897845724931290556
     & 95501782305904531561168206400819094968495942192957731612261142075
     & 562700800000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 8433550.D0/7749.D0*z4 + 611383862259187363793929609521
     & 31074756379535438174267173744362818585529650929100715421151432191
     & 948726202732258467.D0/2609497680686165652885191503615278104585712
     & 70322977409790328677634022558317447255887079267566073972224000000
     & 0.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 4014857.D0/3690.D0*z4 + 353650940802178803947623776688
     & 30030471119130850998129139004571006446022911953867482185548700323
     & 0306295991040555929.D0/155230229424342967310923455469239713085044
     & 82632457511762181986221893140449366904276350170495116437504000000
     & 000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1908947.D0/1755.D0*z4 + 168814770996543696452758405248
     & 99215592489159936214653977003313021984914157426047865778629396845
     & 94008134219573875473.D0/76289978607329531495490429944028785821064
     & 71342537045414828463960271872684262027345572095987234054041600000
     & 0000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 7251772.D0/6669.D0*z4 + 206685705872548147649141322138
     & 54030963623650307155368824409592930236394008378776654790221846933
     & 940412272251.D0/9628393778195445949436727872023723797401346965060
     & 61778605446858537390312851634697915378023731200000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 6877666.D0/6327.D0*z4 + 926120771480266093417251072089
     & 18648507004431707588825153872724489976387949932578341477773531258
     & 432403391601.D0/4453132122415393751614486640810972256298122971340
     & 535726050191720735430196938810477858623359756800000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 3256735.D0/2997.D0*z4 + 725381855508827965734925248881
     & 59095415249682603214380235325773792290362302994497945969731102545
     & 42346017.D0/36051972419389419724824017534120994156184468181439498
     & 5369094363412241036698003698694058680320000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 3079592.D0/2835.D0*z4 + 150646699088273712028077834246
     & 28132178570957676783550461874230325014328017807320216098809689145
     & 76067281.D0/77507311712446541177447212143380515692042775574347816
     & 264854439308528478405590475763464704000000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 5814808.D0/5355.D0*z4 + 389133321283538580411327855148
     & 13096480403792701279350970669538677622189181583902176969256787876
     & 13209.D0/20759202065666227244994383578699432009069163168196395882
     & 9656399414125459571827424737792000000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 5480342.D0/5049.D0*z4 + 541527830983592426247184837962
     & 167375461584633268537008043503433750037404895460782443134798773.D
     & 0/300074722900287812228394366859120504027906070281579038631469528
     & 22952525703527257600000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 2577893.D0/2376.D0*z4 + 133697588415358133757449187331
     & 9351542237276939278326216850937891707829063640052233215133279.D0/
     & 77101410934034521614308638484126143807058738455700870488432672729
     & 740763750195200000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1210285.D0/1116.D0*z4 + 630501887941282791749400154780
     & 91063748452526170396097960844391863625890515340062320398931.D0/
     & 37919681579997335863194351790354661006090585798031679693756990642
     & 45371735183360000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 4536404.D0/4185.D0*z4 + 303171513045082889201854015295
     & 70337042071695469720404841712499262371952588033851917319917.D0/
     & 19059095915542211860399255118509178858318530676500106328589710433
     & 73495205504000000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 4241578.D0/3915.D0*z4 + 241643516072128869122418462101
     & 521199615328172296443998717167895781226554185392187499717.D0/
     & 15919175735907953454826566775252482341337482419887761622381958599
     & 605806728000000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1978331.D0/1827.D0*z4 + 267994889590109501821097871417
     & 8168739580605940742480346183968295894453960602818426163.D0/
     & 18553045602723109125292148583020582957646181384260450590080950709
     & 2180057600000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1840828.D0/1701.D0*z4 + 106020760599949648426548491989
     & 403077935307544705244566531591541127463535030416941.D0/
     & 77369547685785236693709207409569751651306022879622670634853699127
     & 48800000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 3416560.D0/3159.D0*z4 + 105095349355708714836566280272
     & 950035149727294853625910147352300657308287.D0/8112598380385797109
     & 077888708762295689467578801306164433194240000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 3161374.D0/2925.D0*z4 + 366664539113455206909158582842
     & 2437040205942227272519907326958752977991.D0/300565967331312907296
     & 695734665955619986275547153036784862000000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1458049.D0/1350.D0*z4 + 163773888232484333328707922967
     & 8466718053491796640234392099661472471487.D0/143197681706233444667
     & 408389716782826395446414524523728768000000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 670183.D0/621.D0*z4 + 77550927536109193242988282403030
     & 29071033697286214201489653600398333.D0/72694871524773546491891905
     & 9224394127589557091890248619520000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 2455276.D0/2277.D0*z4 + 211454373275178940622966750748
     & 6087846101677746305399640611187331.D0/213751292908983964557821696
     & 109819176783884736989915392000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 2239730.D0/2079.D0*z4 + 277157243249848964264948166047
     & 143082744256602734185281097000699643.D0/3042122340982371239062867
     & 8021890406426302117432363665056000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1017047.D0/945.D0*z4 + 1575709011305449630916951508167
     & 755037500251173822120528697139.D0/1893313752233338883736877502058
     & 46220196088055354880000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 919184.D0/855.D0*z4 + 20799756876149389535047256752599
     & 257486231072882312214927119.D0/2762900252875072103763492514524637
     & 621909119087360000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1652552.D0/1539.D0*z4 + 890456479728954059044572125076
     & 5410579241530791637162903.D0/132349390155492907697857170594551133
     & 8626017024000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1476646.D0/1377.D0*z4 + 182521416351630941183511654726
     & 3844461700669341411763.D0/308196433602015874323808731111485662701
     & 744000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 655325.D0/612.D0*z4 + 17453647372418847982181481986761
     & 213159186955684969.D0/3414918931878292236274888987384882689216000
     & 000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 288641.D0/270.D0*z4 + 28853715331525303151636856261621
     & 7157284043.D0/67185475423939440511621356166488000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 1008388.D0/945.D0*z4 + 8116422230054863987070661427320
     & 156920393.D0/2336678765520740012470550227776000000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 872122.D0/819.D0*z4 + 38232104685944821463935961953319
     & 398453313.D0/14438484134575497583306309266817140000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 372883.D0/351.D0*z4 + 12489412793499460979134398323574
     & 51528563.D0/686831414161462647042964020422400000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 314660.D0/297.D0*z4 + 34049846966652175534042174022161
     & 7507.D0/345520239448279937687396150400000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 522784.D0/495.D0*z4 + 30630573554717343941174375864189
     & 3.D0/2044498458273845785132521600000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 426158.D0/405.D0*z4 - 95829240569637444540144798450517
     & .D0/139397622155034939895399200000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 169721.D0/162.D0*z4 - 2656138027044827776640476181.D0/
     & 1742225276508427751040000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 65659.D0/63.D0*z4 - 11836079780316914712935197.D0/
     & 5018755940559256896000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 195740.D0/189.D0*z4 - 773653809631594649521.D0/
     & 243048444518880000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 138754.D0/135.D0*z4 - 2141099165834337671.D0/
     & 536887165815000.D0 - w(5)
      w(6)=y*w(6)
      w(6)=w(6) + 45839.D0/45.D0*z4 - 83230572524920841.D0/
     & 17531009496000.D0 - w(5)
      w(6)=y*w(6)
      w(8)= - w(5) + 27256.D0/27.D0*z4
      w(6)=w(6) - 23113031631653.D0/4322241000.D0 + w(8)
      w(6)=y*w(6)
      w(6)=w(6) - 490947223.D0/110250.D0 + w(8)
      w(6)=w(6)*w(19)
      w(6)=w(6) - 8774642209.D0/3969000.D0 - 4*z4
      w(6)=y*w(6)
      w(5)=w(6) + 33325.D0/27.D0*z4 - 27479633.D0/14580.D0 - w(5)
      w(5)=y*w(5)
      w(6)=22824310896752140481695112861.D0/24.D0 + 
     & 124064960806850657677039785091.D0/125.D0*y
      w(1)=w(6)*w(1)
      w(1)=104789301643186290852946623733.D0/5640.D0 + w(1)
      w(1)=w(1)*w(9)
      w(1)=13720373908360806185178962963.D0/5405.D0 + w(1)
      w(1)=y*w(1)
      w(1)=12550923310981232663773049383.D0/5175.D0 + w(1)
      w(1)=w(1)*w(21)
      w(1)=243788289108634879008860297.D0/4950.D0 + w(1)
      w(1)=y*w(1)
      w(1)=222092169198001900791583517.D0/4730.D0 + w(1)
      w(1)=w(1)*w(15)
      w(1)=1517905564388643916478009.D0/645.D0 + w(1)
      w(1)=y*w(1)
      w(1)=893105565776442149561453.D0/399.D0 + w(1)
      w(1)=w(1)*w(22)
      w(1)=93953472620113825971799.D0/1900.D0 + w(1)
      w(1)=w(1)*w(25)
      w(1)=8625731147810959384457.D0/5700.D0 + w(1)
      w(1)=w(1)*w(23)
      w(1)=76306370807474508181027.D0/2182245.D0 + w(1)
      w(1)=y*w(1)
      w(1)=1255441482479175211.D0/108965.D0 + 1.D0/2871.D0*w(1)
      w(1)=w(1)*w(9)
      w(1)=8754293869355915615617.D0/5631093270.D0 + w(1)
      w(1)=y*w(1)
      w(1)=7809080498721086408533.D0/5326709850.D0 + w(1)
      w(1)=w(1)*w(24)
      w(1)=52541397323695150543.D0/1408618827.D0 + w(1)
      w(1)=y*w(1)
      w(1)=1163780601372975368539.D0/33203158065.D0 + w(1)
      w(1)=y*w(1)
      w(1)=2367039042297413933.D0/72004680.D0 + w(1)
      w(1)=y*w(1)
      w(1)=451633879721114641553.D0/14678044920.D0 + w(1)
      w(1)=w(1)*w(28)
      w(1)=79091431071083291197.D0/5504266845.D0 + w(1)
      w(1)=y*w(1)
      w(1)=344737665063950778617.D0/25745764275.D0 + w(1)
      w(1)=y*w(1)
      w(1)=526290632606261491.D0/42280434.D0 + w(1)
      w(1)=y*w(1)
      w(1)=8329533067847045051.D0/721683270.D0 + w(1)
      w(1)=y*w(1)
      w(1)=246626568395707967.D0/23108085.D0 + w(1)
      w(1)=y*w(1)
      w(1)=255212998602083.D0/25935.D0 + w(1)
      w(1)=y*w(1)
      w(1)=6786184250701.D0/29260.D0 + 1.D0/39.D0*w(1)
      w(1)=y*w(1)
      w(1)=159339711907373.D0/749892.D0 + w(1)
      w(1)=y*w(1)
      w(1)=8401120811856523.D0/43306263.D0 + w(1)
      w(1)=y*w(1)
      w(1)=2325771987695837.D0/13180167.D0 + w(1)
      w(1)=y*w(1)
      w(1)=1633230607781.D0/3993990.D0 + 1.D0/391.D0*w(1)
      w(1)=y*w(1)
      w(1)=47081792107303.D0/127588230.D0 + w(1)
      w(1)=y*w(1)
      w(1)=359477006579.D0/1084083.D0 + w(1)
      w(1)=y*w(1)
      w(1)=132089864538197.D0/445215771.D0 + w(1)
      w(1)=y*w(1)
      w(1)=35748796607.D0/135252.D0 + w(1)
      w(1)=y*w(1)
      w(1)=165757998671.D0/706860.D0 + w(1)
      w(1)=w(1)*w(28)
      w(1)=98037789511.D0/945945.D0 + w(1)
      w(1)=y*w(1)
      w(1)=24962780833.D0/273273.D0 + w(1)
      w(1)=w(1)*w(7)
      w(1)=11298647309.D0/702702.D0 + w(1)
      w(1)=y*w(1)
      w(1)=42101102713.D0/2972970.D0 + w(1)
      w(1)=y*w(1)
      w(1)=1431883063.D0/114345.D0 + w(1)
      w(1)=y*w(1)
      w(1)=1741482289.D0/155925.D0 + w(1)
      w(1)=y*w(1)
      w(1)=28655251.D0/2835.D0 + w(1)
      w(1)=w(1)*w(28)
      w(1)=10310353.D0/2205.D0 + w(1)
      w(1)=w(1)*w(28)
      w(1)=4922669.D0/2205.D0 + w(1)
      w(1)=y*w(1)
      w(1)=702796.D0/315.D0 + w(1)
      w(1)=y*w(1)
      w(1)=86279.D0/36.D0 + w(1)
      w(1)=y*w(1)
      w(1)=529847.D0/180.D0 + w(1)
      w(1)=y*w(1)
      w(1)=446221.D0/90.D0 + w(1)
      w(1)=y*w(1)
      w(1)=115229.D0/10.D0 + w(1)
      w(1)=y*w(1)
      w(1)= - 104746.D0/5.D0 + w(1)
      w(1)=w(1)*w(19)
      w(1)=10736.D0/5.D0 + w(1)
      w(1)=z3*w(1)
      w(6)= - 913184725341446113223257.D0/72.D0 - 
     & 14224770851718786978751927.D0/1075.D0*y
      w(6)=y*w(6)
      w(6)= - 69762167489019084210823.D0/3096.D0 + 1.D0/539.D0*w(6)
      w(6)=y*w(6)
      w(6)= - 234511058158026120002891.D0/10879.D0 + w(6)
      w(6)=w(6)*w(21)
      w(6)= - 643736831974262111229113.D0/1468665.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 53438232225327044992183.D0/127710.D0 + w(6)
      w(6)=w(6)*w(16)
      w(6)= - 22203451578358738151.D0/946.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 7206449056580634956287.D0/322371.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 9377742633724880911.D0/441.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 6868513184834556227.D0/340.D0 + w(6)
      w(6)=w(6)*w(23)
      w(6)= - 130162432177975953167.D0/278460.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 117137882299081092107.D0/264537.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 946953019177353917.D0/2261.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 363058552609919813.D0/918.D0 + w(6)
      w(6)=w(6)*w(24)
      w(6)= - 2266173342700809251.D0/224910.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 671272878319208117.D0/70805.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 2801938644795077083.D0/314721.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 8649388716259389023.D0/1036728.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 81689283450932051.D0/10472.D0 + w(6)
      w(6)=w(6)*w(28)
      w(6)= - 214303974642825413.D0/58905.D0 + w(6)
      w(6)=w(6)*w(25)
      w(6)= - 83527814736461.D0/765.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 1856460081647779.D0/18326.D0 + w(6)
      w(6)=w(6)*w(10)
      w(6)= - 177646321592771.D0/54978.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 3834288355507.D0/1287.D0 + w(6)
      w(6)=w(6)*w(31)
      w(6)= - 156834723471.D0/1547.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 771630468691.D0/8316.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 478232533565.D0/141372.D0 + 1.D0/25.D0*w(6)
      w(6)=y*w(6)
      w(6)= - 56816363759.D0/18513.D0 + w(6)
      w(6)=w(6)*w(26)
      w(6)= - 1638890276861.D0/13607055.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 1337317912621.D0/12370050.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 8114878817.D0/84150.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 12350489027.D0/144585.D0 + w(6)
      w(6)=w(6)*w(15)
      w(6)= - 6284201681.D0/1590435.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 811008067.D0/235620.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 264482731.D0/89100.D0 + w(6)
      w(6)=w(6)*w(28)
      w(6)= - 125551247.D0/99225.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 25820903.D0/24255.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 7859497.D0/8910.D0 + w(6)
      w(6)=w(6)*w(13)
      w(6)= - 689077.D0/12474.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 254489.D0/5775.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 486887.D0/14175.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 3523.D0/135.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 676.D0/35.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 1468.D0/105.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 1352.D0/135.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 323.D0/45.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 43.D0/9.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 14.D0/27.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 146.D0/9.D0 + w(6)
      w(6)=y*w(6)
      w(6)= - 40.D0/3.D0 + w(6)
      w(6)=ln2*w(6)*w(33)
*
      G = 1.D0/9.D0*w(1) + w(2) + w(3) + 1.D0/3.D0*w(4) + w(5) + w(6)
*
      AGGREG1 = G*as**3
*
      RETURN
      END
      REAL*8 FUNCTION AGGPLU(z,nf,as,LL)
*     ---------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     plus part of agg3 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(5)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
      NF=nf
      w(1)=1.0D0/(729.0D0 - 729.0D0*z)
      w(2)=135*z2
      w(3)=270*z3 + 314 + w(2)
      w(4)= - 378*z3 + 1736 + w(2)
      w(4)=NF*w(4)
      w(3)=7*w(3) + 2*w(4)
      w(4)=4*TF
      w(3)=w(3)*w(4)
      w(2)= - 324*z4 + 453*z3 - w(2) - 769 + 72*B4
      w(2)=CF*w(2)
      w(5)=18954*z4 - 12123*z3 - 7308*z2 - 8141 - 1944*B4
      w(5)=CA*w(5)
      w(2)=w(3) + 54*w(2) + w(5)
*
      AGGPLU = CA*w(4)*w(2)*w(1)*as**3
*
      RETURN
      END      
      REAL*8 FUNCTION AGGDEL(z,nf,as,LL)
*     ---------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Delta[1-z] part of agg3 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(3)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      w(1)=44*z2
      w(2)= - 971.D0/27.D0 - 10*z2
      w(2)=NF*w(2)
      w(2)=w(2) - 739.D0/27.D0 - w(1)
      w(2)= - 7*z3 + 2.D0/3.D0*w(2)
      w(2)=CF*w(2)
      w(1)=199.D0/5.D0 + w(1)
      w(3)= - 89 + 98.D0/3.D0*z2
      w(3)=NF*w(3)
      w(1)=2.D0/9.D0*w(3) + 13.D0/27.D0*w(1) - 291.D0/10.D0*z3
      w(1)=CA*w(1)
      w(3)=TF*z3
      w(1)=64.D0/27.D0*w(3) + w(2) + w(1)
      w(1)=TF*w(1)
      w(2)=34315.D0/12.D0 + 992*z2
      w(3)=20435.D0/216.D0 + 24*z2
      w(3)=z3*w(3)
      w(2)=32.D0/3.D0*B4 - 3778.D0/27.D0*z4 - 304.D0/9.D0*z5 + 1.D0/27.D
     & 0*w(2) + w(3)
      w(2)=CA*w(2)
      w(3)= - 64.D0/3.D0*B4 + 128.D0/3.D0*z4 - 2617.D0/12.D0*z3 + 16541.
     & D0/162.D0 + 52*z2
      w(3)=CF*w(3)
      w(2)=w(3) + w(2)
      w(2)=CA*w(2)
      w(3)=274.D0/3.D0 + 95*z3
      w(3)=w(3)*CF**2
      w(1)=w(1) + 1.D0/3.D0*w(3) + w(2)
*
      AGGDEL = TF*w(1)*as**3
*
      RETURN
      END      
      REAL*8 FUNCTION APS1(z,nf,as,LL)
*     --------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of aQqPS3: contribution of only H[a,z]
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(124)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4,li5half
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      li5half = 0.50840057924226870746D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
      w(1)=z**(-1)
      w(2)=Hr1(-1)
      w(3)=Hr1(0)
      w(4)=Hr2(0,-1)
      w(5)=Hr2(0,1)
      w(6)=Hr3(0,-1,-1)
      w(7)=Hr3(0,-1,1)
      w(8)=Hr3(0,0,-1)
      w(9)=Hr3(0,0,1)
      w(10)=Hr3(0,1,-1)
      w(11)=Hr3(0,1,1)
      w(12)=Hr1(1)
      w(13)=Hr4(0,-1,-1,-1)
      w(14)=Hr4(0,-1,0,1)
      w(15)=Hr4(0,0,-1,-1)
      w(16)=Hr4(0,0,-1,1)
      w(17)=Hr4(0,0,0,-1)
      w(18)=Hr4(0,0,0,1)
      w(19)=Hr4(0,0,1,-1)
      w(20)=Hr4(0,0,1,1)
      w(21)=Hr4(0,1,1,1)
      w(22)=Hr5(0,-1,-1,0,1)
      w(23)=Hr4(0,-1,-1,1)
      w(24)=Hr5(0,-1,0,-1,-1)
      w(25)=Hr4(0,-1,1,-1)
      w(26)=Hr4(0,-1,1,1)
      w(27)=Hr5(0,0,-1,-1,-1)
      w(28)=Hr5(0,0,-1,0,-1)
      w(29)=Hr5(0,0,-1,0,1)
      w(30)=Hr5(0,0,0,-1,-1)
      w(31)=Hr5(0,0,0,-1,1)
      w(32)=Hr5(0,0,0,0,-1)
      w(33)=Hr5(0,0,0,0,1)
      w(34)=Hr5(0,0,0,1,-1)
      w(35)=Hr5(0,0,0,1,1)
      w(36)=Hr5(0,0,1,0,-1)
      w(37)=Hr5(0,0,1,0,1)
      w(38)=Hr5(0,0,1,1,1)
      w(39)=Hr4(0,1,-1,-1)
      w(40)=Hr4(0,1,-1,1)
      w(41)=Hr5(0,1,0,1,1)
      w(42)=Hr4(0,1,1,-1)
      w(43)=Hr5(0,1,1,1,1)
      w(44)=1/( - 81 + 81*z)
      w(45)=1/( - 81 + 81*z**2)
      w(46)=1/(81*z + 81*z**2)
      w(47)=4*w(3)
      w(48)=5*z
      w(49)=17.D0/3.D0 - w(48)
      w(49)=w(49)*w(47)
      w(50)=56.D0/3.D0*w(1)
      w(51)=8.D0/3.D0*z
      w(52)=3 - w(51)
      w(52)=z*w(52)
      w(49)=w(49) + w(50) - 9 + w(52)
      w(49)=w(3)*w(49)
      w(52)=z - 1
      w(53)=w(52)*w(4)
      w(54)=z + 1
      w(55)=w(54)*w(5)
      w(56)=721 - 401*z
      w(56)=z*w(56)
      w(49)=38*w(53) + 16.D0/3.D0*w(55) + w(49) + 256.D0/27.D0*w(1) - 
     & 57 + 2.D0/9.D0*w(56)
      w(49)=w(8)*w(49)
      w(56)=w(52)*w(3)
      w(57)= - 11 + 80.D0/3.D0*z
      w(57)=z*w(57)
      w(57)=146.D0/3.D0*w(56) - w(50) + 15 + w(57)
      w(57)=w(17)*w(57)
      w(58)=140.D0/9.D0*w(1)
      w(59)= - 21 - 148.D0/9.D0*z
      w(59)=z*w(59)
      w(59)=w(58) + 5 + w(59)
      w(60)=w(54)*w(3)
      w(59)=2*w(59) + 13*w(60)
      w(59)=w(3)*w(59)
      w(61)= - 3725 - 356*z
      w(61)=z*w(61)
      w(61)= - 272*w(1) + 3469 + w(61)
      w(59)=1.D0/27.D0*w(61) + w(59)
      w(59)=w(11)*w(59)
      w(61)=8*w(3)
      w(62)=3*z
      w(63)= - 7.D0/3.D0 + w(62)
      w(63)=w(63)*w(61)
      w(64)=136.D0/3.D0*z
      w(65)= - 61 + w(64)
      w(65)=z*w(65)
      w(65)=104.D0/3.D0*w(1) - 25 + w(65)
      w(63)=1.D0/3.D0*w(65) + w(63)
      w(63)=w(14)*w(63)
      w(65)= - 71 + 200.D0/3.D0*z
      w(65)=w(65)*z
      w(66)=200.D0/3.D0*w(1)
      w(65)=w(65) + w(66) - 71
      w(67)= - 1.D0/3.D0*w(65) + 40*w(56)
      w(67)=w(15)*w(67)
      w(68)=369179 - 1024228.D0/3.D0*z
      w(68)=z*w(68)
      w(68)=274870.D0/3.D0*w(1) - 119393 + w(68)
      w(69)= - 79 + 63*z
      w(70)=w(34)*w(69)
      w(49)=w(59) + 1.D0/243.D0*w(68) + 2*w(70) + w(57) + w(63) + w(67)
     &  + w(49)
      w(57)=w(47)*w(54)
      w(59)=4.D0/3.D0*z
      w(63)=w(59) - 5
      w(67)=2*z
      w(68)=w(63)*w(67)
      w(70)=w(54)*ln2
      w(71)=2.D0/3.D0*w(1)
      w(57)= - 8*w(70) - w(57) - 5 + w(68) - w(71)
      w(57)=w(57)*B4
      w(68)=68*z
      w(72)=1577 - w(68)
      w(72)=z*w(72)
      w(72)=364*w(1) + 260 + w(72)
      w(73)= - 1 - w(62)
      w(73)=w(3)*w(73)
      w(72)=1.D0/9.D0*w(72) + 31*w(73)
      w(72)=w(18)*w(72)
      w(73)=197.D0/3.D0 + 137*z
      w(73)=w(33)*w(73)
      w(74)= - 79.D0/3.D0 + 21*z
      w(74)=w(36)*w(74)
      w(69)=w(31)*w(69)
      w(69)=w(72) + w(57) + w(73) + w(74) + w(69)
      w(72)=3 - 112.D0/9.D0*z
      w(72)=z*w(72)
      w(72)=70.D0/9.D0*w(56) - w(50) + 19.D0/3.D0 + w(72)
      w(72)=w(3)*w(72)
      w(73)=2.D0/3.D0*z
      w(74)= - 359 + 551.D0/3.D0*z
      w(74)=w(74)*w(73)
      w(75)=256.D0/9.D0*w(1)
      w(74)= - w(75) + 85 + w(74)
      w(72)=2.D0/3.D0*w(74) + w(72)
      w(72)=w(3)*w(72)
      w(63)=w(63)*z
      w(74)=4.D0/3.D0*w(1)
      w(63)= - w(74) + w(63) + 5
      w(76)=2*w(3)
      w(77)=w(76)*w(54)
      w(78)= - w(77) - w(63)
      w(78)=w(5)*w(78)
      w(79)= - 9 - 8.D0/9.D0*z
      w(79)=z*w(79)
      w(79)= - 19.D0/3.D0*w(56) - 80.D0/9.D0*w(1) + 29.D0/3.D0 + w(79)
      w(80)=2*w(4)
      w(79)=w(79)*w(80)
      w(81)=2623.D0/3.D0 + 1400*z
      w(81)=w(81)*z
      w(81)=w(81) - 341.D0/3.D0 + 412*w(1)
      w(82)=w(52)*w(6)
      w(72)=w(79) + 112.D0/3.D0*w(82) + 8.D0/3.D0*w(78) + 2.D0/9.D0*
     & w(81) + w(72)
      w(72)=w(72)*w(80)
      w(78)=104*z
      w(79)=101 + w(78)
      w(79)=z*w(79)
      w(79)= - 104*w(1) - 64 + w(79)
      w(79)=1.D0/3.D0*w(79) - 16*w(60)
      w(79)=w(21)*w(79)
      w(83)=40*z
      w(84)=617 + w(83)
      w(84)=z*w(84)
      w(84)= - 164*w(1) + 119 + w(84)
      w(85)=37 + 35*z
      w(85)=w(3)*w(85)
      w(84)=1.D0/3.D0*w(84) + w(85)
      w(84)=w(20)*w(84)
      w(85)= - 347 - 341*z
      w(85)=w(35)*w(85)
      w(79)=w(85) + w(79) + w(84)
      w(84)=32.D0/3.D0*z
      w(85)=w(84) - 1
      w(85)=w(85)*z
      w(86)=40.D0/3.D0*w(1)
      w(85)=w(85) + w(86) - 11
      w(85)=w(85)*w(76)
      w(87)= - 13 + 41.D0/3.D0*z
      w(87)=w(87)*w(67)
      w(88)=80.D0/3.D0*w(1)
      w(87)=w(87) - w(88) - 89
      w(85)=w(85) - 1.D0/3.D0*w(87)
      w(89)=w(10) + w(7)
      w(89)= - 8.D0/3.D0*w(89)
      w(85)=w(85)*w(89)
      w(89)= - 178 + 77.D0/3.D0*z
      w(89)=w(89)*w(67)
      w(89)= - 844.D0/3.D0*w(1) - 217.D0/3.D0 + w(89)
      w(90)=185.D0/3.D0 - 16*z
      w(90)=z*w(90)
      w(90)= - 62.D0/9.D0*w(60) + 344.D0/9.D0*w(1) + 157.D0/3.D0 + 
     & w(90)
      w(90)=w(3)*w(90)
      w(89)=2.D0/3.D0*w(89) + w(90)
      w(89)=w(3)*w(89)
      w(90)= - 67 + 224.D0/3.D0*z
      w(90)=z*w(90)
      w(90)= - 98*w(60) - 112.D0/3.D0*w(1) - 35 + w(90)
      w(91)=1.D0/3.D0*w(5)
      w(90)=w(90)*w(91)
      w(89)=w(89) + w(90)
      w(90)=2*w(5)
      w(89)=w(89)*w(90)
      w(92)=17 + 56.D0/3.D0*z
      w(92)=z*w(92)
      w(66)= - 19*w(56) + w(66) - 95 + w(92)
      w(66)=w(3)*w(66)
      w(92)=35 + 188.D0/3.D0*z
      w(92)=w(92)*z
      w(92)=w(92) + w(88) - 1
      w(93)=8.D0/3.D0*w(92)
      w(66)=w(93) + w(66)
      w(66)=w(6)*w(66)
      w(94)=416*w(70)
      w(95)= - 3775 - 584*z
      w(95)=z*w(95)
      w(95)= - 1828.D0/3.D0*w(1) - 467 + w(95)
      w(96)=5 - 17.D0/3.D0*z
      w(96)=w(3)*w(96)
      w(95)= - w(94) + 1.D0/3.D0*w(95) + 38*w(96)
      w(95)=z4*w(95)
      w(96)=w(59) - 1
      w(96)=w(96)*z
      w(96)=w(96) + w(74) - 1
      w(97)= - w(52)*w(76)
      w(97)=w(97) + w(96)
      w(97)=w(13)*w(97)
      w(98)=w(53) - w(55)
      w(99)=25.D0/3.D0 + 23*z
      w(99)=w(99)*w(76)
      w(100)= - 535 + w(64)
      w(100)=z*w(100)
      w(100)= - 536.D0/3.D0*w(1) - 197 + w(100)
      w(99)=1.D0/3.D0*w(100) + w(99)
      w(99)=w(3)*w(99)
      w(100)=3913 + 1936*z
      w(100)=z*w(100)
      w(100)=3019 + w(100)
      w(100)=1.D0/3.D0*w(100) + 1004*w(1)
      w(98)=1.D0/9.D0*w(100) + w(99) - 142.D0/3.D0*w(98)
      w(99)=4*w(9)
      w(98)=w(98)*w(99)
      w(100)=16.D0/3.D0*w(43)
      w(101)=ln2**5
      w(102)=80*w(38) + 128.D0/15.D0*w(101) + 176.D0/3.D0*w(41) - 1024*
     & li5half + w(100)
      w(102)=w(54)*w(102)
      w(103)= - w(25) - w(23) - w(39)
      w(84)=w(84) + 31
      w(84)=w(84)*z
      w(84)=w(84) + 31 + 32.D0/3.D0*w(1)
      w(104)=8.D0/3.D0*w(84)
      w(103)=w(104)*w(103)
      w(104)=w(73) + 1
      w(104)=w(104)*z
      w(104)=w(104) + w(71) + 1
      w(105)=32.D0/3.D0*w(104)
      w(106)= - w(40) - w(42) - w(26)
      w(105)=w(105)*w(106)
      w(106)=88.D0/3.D0*z
      w(107)=w(106) - 27
      w(107)=w(107)*z
      w(108)=88.D0/3.D0*w(1)
      w(107)=w(107) + w(108) - 27
      w(109)=w(107) + 16.D0/3.D0*w(60)
      w(110)=w(16) + w(19)
      w(110)=4*w(110)
      w(109)=w(109)*w(110)
      w(110)= - 936*w(30) - 448.D0/3.D0*w(24) + 80.D0/3.D0*w(22) - 312*
     & w(28) - 896.D0/3.D0*w(27)
      w(110)=w(52)*w(110)
      w(111)=1 - 4.D0/5.D0*z
      w(111)=w(111)*w(76)
      w(112)= - 172 + 95*z
      w(111)=1.D0/9.D0*w(112) + w(111)
      w(111)=w(3)*w(111)
      w(112)=1213.D0/3.D0 + 430*z
      w(112)=z*w(112)
      w(112)=1045.D0/3.D0 + w(112)
      w(111)=2.D0/9.D0*w(112) + w(111)
      w(112)=w(3)**2
      w(111)=w(111)*w(112)
      w(113)=110993 + 503464*z
      w(113)=z*w(113)
      w(113)=20992*w(1) + 171290 + w(113)
      w(111)=2.D0/81.D0*w(113) + w(111)
      w(113)=1.D0/3.D0*w(3)
      w(111)=w(111)*w(113)
      w(114)=2957 - 11356*z
      w(114)=w(114)*w(67)
      w(114)=9627 + w(114)
      w(114)=z*w(114)
      w(114)= - 5671 + w(114)
      w(114)=z*w(114)
      w(114)=13652 + w(114)
      w(115)=2*w(112)
      w(114)=w(45)*w(114)*w(115)
      w(68)= - 67 - w(68)
      w(68)=w(37)*w(68)
      w(116)= - 41.D0/3.D0 + 11*z
      w(116)=w(29)*w(116)
      w(117)=145 - 47*z
      w(117)=z5*w(117)
      w(62)=1 - w(62)
      w(62)=w(32)*w(62)
      w(49)=w(98) + 6*w(117) + 112.D0/3.D0*w(97) + 16*w(116) + w(72) + 
     & w(95) + 4.D0/3.D0*w(66) + w(89) + 16.D0/3.D0*w(68) + w(111) + 
     & w(110) + w(109) + w(85) + w(105) + w(103) + w(102) + 8.D0/3.D0*
     & w(79) + 8*w(69) + 4*w(49) + 64*w(62) + w(114)
      w(49)=CA*w(49)
      w(62)=28*w(60)
      w(66)= - 109 + 196.D0/3.D0*z
      w(66)=z*w(66)
      w(66)= - w(62) + w(50) - 71 + w(66)
      w(68)=4.D0/3.D0*w(5)
      w(66)=w(66)*w(68)
      w(68)=ln2**3
      w(69)= - 304.D0/3.D0*w(9) - 248.D0/3.D0*w(11) + 512.D0/3.D0*w(68)
      w(69)=w(54)*w(69)
      w(72)=20*w(1)
      w(79)= - 7001 - 8896.D0/3.D0*z
      w(79)=z*w(79)
      w(79)= - 43 + w(79)
      w(79)=1.D0/27.D0*w(79) + w(72)
      w(85)= - 11 - 68.D0/9.D0*z
      w(85)=z*w(85)
      w(85)=19.D0/3.D0 + w(85)
      w(85)=2*w(85) + 47.D0/9.D0*w(60)
      w(85)=w(85)*w(76)
      w(89)=1879.D0/3.D0 + 512*z
      w(89)=z*w(89)
      w(89)=1943.D0/3.D0 + w(89)
      w(85)=1.D0/3.D0*w(89) + w(85)
      w(85)=w(3)*w(85)
      w(66)=w(66) + 2*w(79) + w(85) + w(69)
      w(66)=z2*w(66)
      w(69)=w(59) + 1
      w(69)=w(69)*z
      w(69)=w(69) - w(74) - 1
      w(79)=w(69)*w(3)
      w(85)=392.D0/3.D0*w(1)
      w(89)=593 + 904.D0/3.D0*z
      w(89)=z*w(89)
      w(89)= - 20*w(79) + w(85) - 1025 + w(89)
      w(89)=w(3)*w(89)
      w(95)= - 1015 - 1312.D0/3.D0*z
      w(95)=z*w(95)
      w(95)= - 722.D0/3.D0*w(1) + 1693 + w(95)
      w(89)=4.D0/3.D0*w(95) + w(89)
      w(89)=w(3)*w(89)
      w(95)= - 1387 + 67*z
      w(95)=z*w(95)
      w(95)=1490*w(1) - 170 + w(95)
      w(89)=4.D0/9.D0*w(95) + w(89)
      w(95)= - 148*z3 + 224*w(9) + 152*w(11)
      w(95)=w(69)*w(95)
      w(97)=1.D0/3.D0*z
      w(98)= - 17 + w(97)
      w(98)=z*w(98)
      w(98)=17 + w(98)
      w(98)= - 7*w(79) + 2*w(98) - 11.D0/3.D0*w(1)
      w(98)=w(5)*w(98)
      w(102)=1399 + 752.D0/3.D0*z
      w(102)=z*w(102)
      w(102)=652.D0/3.D0*w(1) - 1723 + w(102)
      w(102)=1.D0/3.D0*w(102) + 28*w(79)
      w(102)=z2*w(102)
      w(89)=w(102) + 16*w(98) + 1.D0/3.D0*w(89) + w(95)
      w(95)=w(69)*w(76)
      w(98)=59 - 122.D0/3.D0*z
      w(98)=z*w(98)
      w(98)= - w(86) - 5 + w(98)
      w(98)=1.D0/3.D0*w(98) - w(95)
      w(98)=w(3)*w(98)
      w(102)=335 + 913.D0/9.D0*z
      w(102)=z*w(102)
      w(102)=320.D0/9.D0*w(1) - 472 + w(102)
      w(98)=1.D0/3.D0*w(102) + w(98)
      w(102)=4*w(5)
      w(103)= - w(69)*w(102)
      w(98)=1.D0/3.D0*w(98) + w(103)
      w(103)=w(69)*z2
      w(98)=2*w(98) + 31.D0/3.D0*w(103)
      w(105)=4*w(1)
      w(83)=1 + w(83)
      w(83)=z*w(83)
      w(83)= - w(105) - 37 + w(83)
      w(109)= - w(69)*w(47)
      w(83)=1.D0/3.D0*w(83) + w(109)
      w(109)=w(69)*w(12)
      w(83)=2*w(83) + w(109)
      w(110)=1.D0/9.D0*w(12)
      w(83)=w(83)*w(110)
      w(83)=2*w(98) + w(83)
      w(83)=w(12)*w(83)
      w(83)=2.D0/3.D0*w(89) + w(83)
      w(83)=w(12)*w(83)
      w(89)=20.D0/9.D0*w(1)
      w(98)=4.D0/9.D0*z
      w(111)= - 25 - w(98)
      w(111)=z*w(111)
      w(111)=10.D0/9.D0*w(60) - w(89) - 6 + w(111)
      w(111)=w(3)*w(111)
      w(114)= - 509.D0/3.D0 - 98*z
      w(114)=z*w(114)
      w(114)= - 196.D0/9.D0*w(1) + 788.D0/3.D0 + w(114)
      w(111)=1.D0/3.D0*w(114) + w(111)
      w(111)=w(3)*w(111)
      w(51)=w(51) - 1
      w(114)=13*z
      w(116)= - w(51)*w(114)
      w(50)=w(62) + w(50) + 35 + w(116)
      w(50)=w(50)*w(91)
      w(62)=308 + 3293.D0/27.D0*z
      w(62)=z*w(62)
      w(62)=722.D0/27.D0*w(1) - 581.D0/3.D0 + w(62)
      w(91)=w(11)*w(54)
      w(50)=w(50) + 8*w(91) + 1.D0/3.D0*w(62) + w(111)
      w(50)=w(5)*w(50)
      w(62)=4.D0/3.D0*w(60)
      w(91)= - 27 - 100.D0/3.D0*z
      w(91)=z*w(91)
      w(91)=w(62) + 304.D0/9.D0*w(1) + 73.D0/3.D0 + w(91)
      w(91)=w(21)*w(91)
      w(50)=w(91) + w(50)
      w(91)= - 31 + w(73)
      w(91)=w(91)*w(59)
      w(91)= - 56.D0/3.D0*w(60) + w(86) - 9 + w(91)
      w(91)=w(20)*w(91)
      w(57)=w(57) - w(91)
      w(91)=w(59) - 3
      w(91)=w(91)*z
      w(70)=w(71) - w(91) + 2*w(70) + w(77) + 2
      w(70)=w(70)*ln2
      w(91)=275 - 76*z
      w(91)=z*w(91)
      w(91)=112 + w(91)
      w(91)=2*w(91) + 115*w(60)
      w(91)=w(91)*w(76)
      w(111)=1375 - 2824.D0/3.D0*z
      w(111)=z*w(111)
      w(111)= - 646*w(1) + 5777 + w(111)
      w(91)=1.D0/3.D0*w(111) + w(91)
      w(111)=z2*w(54)
      w(91)=236.D0/3.D0*w(111) + 296.D0/3.D0*w(55) + 1.D0/3.D0*w(91) + 
     & 224*w(70)
      w(111)=2*z3
      w(91)=w(91)*w(111)
      w(116)=23.D0/3.D0 + 12*z
      w(116)=z*w(116)
      w(116)=2.D0/3.D0*w(60) - w(86) - 32.D0/3.D0 + w(116)
      w(116)=w(3)*w(116)
      w(117)=152.D0/3.D0*w(1)
      w(118)=581 + 209.D0/3.D0*z
      w(118)=z*w(118)
      w(118)=w(117) - 524 + w(118)
      w(116)=1.D0/9.D0*w(118) + w(116)
      w(118)=16*w(11)
      w(116)=w(116)*w(118)
      w(119)=1223 - 628.D0/3.D0*z
      w(119)=w(119)*w(97)
      w(94)=w(94) + 232.D0/3.D0*w(60) - 60*w(1) + 39 + w(119)
      w(119)=2*z4
      w(94)=w(94)*w(119)
      w(120)=4*z
      w(121)= - 73 + w(114)
      w(121)=w(121)*w(120)
      w(122)=20.D0/3.D0*w(1)
      w(121)=20*w(60) - w(122) - 85 + w(121)
      w(123)=16.D0/3.D0*w(18)
      w(121)=w(121)*w(123)
      w(124)=1091 + 1652.D0/3.D0*z
      w(124)=z*w(124)
      w(85)=w(85) - 2555 + w(124)
      w(124)=179 - 28.D0/3.D0*z
      w(124)=z*w(124)
      w(124)=w(122) + 41 + w(124)
      w(124)=2.D0/3.D0*w(124) - 3*w(60)
      w(124)=w(124)*w(76)
      w(85)= - 224.D0/3.D0*w(55) + 1.D0/9.D0*w(85) + w(124)
      w(85)=w(85)*w(99)
      w(100)=2048*li5half + 1696*w(35) - w(100) - 1772*z5 - 416*w(38)
     &  + 2144.D0/3.D0*w(37) - 256.D0/15.D0*w(101) - 1280.D0/3.D0*w(33)
     &  - 800.D0/3.D0*w(41)
      w(100)=w(54)*w(100)
      w(101)=5 + 76.D0/3.D0*z
      w(101)=z*w(101)
      w(101)=11 + w(101)
      w(101)=2.D0/9.D0*w(101) - 3.D0/5.D0*w(60)
      w(101)=w(3)*w(101)
      w(124)=313 - 704*z
      w(124)=z*w(124)
      w(124)= - 151 + w(124)
      w(101)=1.D0/27.D0*w(124) + w(101)
      w(101)=w(101)*w(112)
      w(112)= - 1187 + 11500*z
      w(112)=z*w(112)
      w(112)=3989 + w(112)
      w(101)=2.D0/81.D0*w(112) + w(101)
      w(101)=w(3)*w(101)
      w(112)= - 5311 - 4376*z
      w(112)=z*w(112)
      w(112)=9879 + w(112)
      w(112)=z*w(112)
      w(112)= - 840 + w(112)
      w(112)=w(44)*w(112)*w(115)
      w(115)=1589.D0/27.D0 - 200*z
      w(115)=z*w(115)
      w(115)=1538.D0/9.D0*w(1) - 803.D0/27.D0 + w(115)
      w(50)=w(83) + w(91) + w(66) + w(85) + w(121) + w(94) + w(112) + 
     & w(116) + 4.D0/3.D0*w(115) + w(101) + w(100) - 16*w(57) + 8*w(50)
      w(50)=CF*w(50)
      w(57)=20*w(11) - 128*w(68)
      w(57)=w(54)*w(57)
      w(48)=4 - w(48)
      w(48)=w(48)*w(47)
      w(66)=8*z
      w(68)=23 + w(66)
      w(68)=z*w(68)
      w(48)=w(48) - 211 + w(68)
      w(48)=w(3)*w(48)
      w(68)= - 2183 + 1318*z
      w(68)=z*w(68)
      w(68)= - 275 + w(68)
      w(48)=w(48) + 1.D0/3.D0*w(68) + 80*w(1)
      w(48)=w(48)*w(113)
      w(48)=w(48) + w(57)
      w(57)=31 - w(66)
      w(57)=w(57)*w(120)
      w(57)=55 + w(57)
      w(57)=w(62) + 1.D0/3.D0*w(57) - 8*w(1)
      w(57)=w(57)*w(90)
      w(62)=15 + 76.D0/9.D0*z
      w(62)=z*w(62)
      w(56)= - 38.D0/3.D0*w(56) + 148.D0/9.D0*w(1) - 7.D0/3.D0 + w(62)
      w(56)=w(56)*w(80)
      w(62)=1 + 7.D0/3.D0*z
      w(62)=w(62)*w(99)
      w(52)=w(8)*w(52)
      w(48)=136.D0/3.D0*w(52) + w(62) + w(56) - 32*w(82) + 1.D0/3.D0*
     & w(48) + w(57)
      w(52)=z2*CA
      w(48)=w(48)*w(52)
      w(56)=80.D0/9.D0 - w(114)
      w(56)=w(56)*w(66)
      w(57)=17 + 152.D0/3.D0*z
      w(57)=z*w(57)
      w(57)= - w(117) - 17 + w(57)
      w(57)=w(3)*w(57)
      w(56)=1.D0/9.D0*w(57) + 76*w(1) - 379.D0/9.D0 + w(56)
      w(56)=w(3)*w(56)
      w(57)=w(4)*w(63)
      w(62)= - 7523 + 21512.D0/3.D0*z
      w(62)=z*w(62)
      w(62)= - 7436.D0/3.D0*w(1) + 2831 + w(62)
      w(56)= - 8.D0/3.D0*w(57) + 2.D0/27.D0*w(62) + w(56)
      w(56)=w(3)*w(56)
      w(57)=79 - 167.D0/9.D0*z
      w(57)=z*w(57)
      w(57)=w(58) - 73 + w(57)
      w(58)=184.D0/3.D0*w(1)
      w(62)=55 + 184.D0/3.D0*z
      w(62)=z*w(62)
      w(62)= - w(58) - 55 + w(62)
      w(62)=w(3)*w(62)
      w(57)=2*w(57) + w(62)
      w(57)=w(5)*w(57)
      w(62)= - 961 - 6776.D0/9.D0*z
      w(62)=w(62)*w(97)
      w(62)= - 11926.D0/27.D0*w(1) + 1013 + w(62)
      w(66)=w(11)*w(69)
      w(63)=w(8)*w(63)
      w(56)=16.D0/3.D0*w(63) + 2.D0/3.D0*w(57) - 44.D0/3.D0*w(66) + 1.D0
     & /9.D0*w(62) + w(56)
      w(56)=CA*w(56)
      w(57)=49 + w(106)
      w(57)=z*w(57)
      w(57)= - w(108) - 49 + w(57)
      w(57)=z3*w(57)
      w(62)= - 65 - 296.D0/3.D0*z
      w(62)=z*w(62)
      w(62)=296.D0/3.D0*w(1) + 65 + w(62)
      w(62)=w(9)*w(62)
      w(57)=w(57) + w(62)
      w(62)=2.D0/3.D0*CA
      w(57)=w(62)*w(57)
      w(63)=z**2
      w(66)= - w(63) + w(1)
      w(61)=w(66)*w(61)
      w(66)= - 1049 - 214.D0/3.D0*z
      w(66)=z*w(66)
      w(66)= - 326.D0/3.D0*w(1) + 1193 + w(66)
      w(61)=1.D0/3.D0*w(66) + w(61)
      w(61)=w(61)*w(52)
      w(56)=1.D0/3.D0*w(61) + w(57) + w(56)
      w(57)= - 1 - 154.D0/9.D0*z
      w(57)=z*w(57)
      w(57)=16*w(79) + 46.D0/9.D0*w(1) + 13 + w(57)
      w(57)=2*w(57) - w(109)
      w(61)=1.D0/3.D0*w(12)
      w(57)=w(61)*w(57)
      w(66)= - 47 + 481.D0/3.D0*z
      w(66)=z*w(66)
      w(58)= - w(58) - 52 + w(66)
      w(64)= - 49 - w(64)
      w(64)=z*w(64)
      w(64)=136.D0/3.D0*w(1) + 49 + w(64)
      w(64)=w(3)*w(64)
      w(58)=4.D0/3.D0*w(58) + w(64)
      w(58)=w(3)*w(58)
      w(64)= - 805 - 1754.D0/3.D0*z
      w(64)=z*w(64)
      w(64)= - 658.D0/3.D0*w(1) + 1609 + w(64)
      w(57)=w(57) + 2.D0/9.D0*w(64) + w(58)
      w(57)=CA*w(57)
      w(58)=w(69)*w(52)
      w(57)= - 10*w(58) + w(57)
      w(57)=w(57)*w(61)
      w(56)=2*w(56) + w(57)
      w(56)=w(12)*w(56)
      w(57)=1.D0/5.D0*z
      w(58)=w(57) - 2
      w(58)=w(58)*z
      w(58)=w(58) + 1
      w(64)=w(58)*w(97)
      w(64)=7 + w(64)
      w(64)=w(64)*w(67)
      w(64)= - 197.D0/9.D0 + w(64)
      w(64)=z*w(64)
      w(66)=z + 2.D0/3.D0
      w(66)=w(66)*w(67)
      w(66)=w(66) - 5.D0/3.D0
      w(68)= - w(60) + 2*w(66)
      w(80)= - w(68)*w(113)
      w(64)=w(80) + 64.D0/9.D0 + w(64)
      w(64)=w(3)*w(64)
      w(80)=w(67) - 19
      w(82)=w(80)*z
      w(83)= - 334 - w(82)
      w(83)=z*w(83)
      w(83)= - 3502.D0/9.D0 + w(83)
      w(83)=w(83)*w(57)
      w(83)= - 155.D0/9.D0 + w(83)
      w(64)=2.D0/3.D0*w(83) + w(64)
      w(64)=w(3)*w(64)
      w(80)= - w(80)*w(67)
      w(80)=12971.D0/27.D0 + w(80)
      w(80)=z*w(80)
      w(80)= - 4678.D0/9.D0 + w(80)
      w(80)=w(80)*w(57)
      w(80)= - 607.D0/27.D0*w(1) + 212.D0/9.D0 + w(80)
      w(83)=w(54)*w(21)
      w(64)=2*w(83) + 2.D0/3.D0*w(80) + w(64)
      w(80)=2*w(1)
      w(85)= - 53.D0/3.D0 - w(120)
      w(85)=z*w(85)
      w(85)= - 35.D0/3.D0 + w(85)
      w(85)=w(60) + 1.D0/3.D0*w(85) + w(80)
      w(85)=w(11)*w(85)
      w(91)= - w(119) + 2.D0/3.D0*w(20)
      w(91)=w(54)*w(91)
      w(63)= - w(58)*w(63)
      w(63)=103.D0/3.D0 + w(63)
      w(63)=z*w(63)
      w(63)= - 61.D0/15.D0 + w(63)
      w(63)=2*w(63) - 5*w(1)
      w(94)=11 + w(67)
      w(94)=z*w(94)
      w(94)=8 + w(94)
      w(94)=w(3)*w(94)
      w(63)=2.D0/3.D0*w(63) + w(94)
      w(63)=1.D0/3.D0*w(63) - w(55)
      w(63)=w(63)*w(90)
      w(63)=w(63) + 1.D0/3.D0*w(64) + 2*w(85) + w(91)
      w(64)= - w(3)*w(68)
      w(68)=37*z
      w(85)=64.D0/3.D0 + w(68)
      w(85)=z*w(85)
      w(64)=w(64) + 100.D0/3.D0 + w(85)
      w(64)=w(64)*w(113)
      w(85)=50*z
      w(91)=569.D0/27.D0 - w(85)
      w(91)=z*w(91)
      w(64)=w(64) + 218.D0/27.D0 + w(91)
      w(64)=w(3)*w(64)
      w(91)= - 2599 + 5165.D0/3.D0*z
      w(91)=z*w(91)
      w(91)= - 1736.D0/3.D0*w(1) + 1456 + w(91)
      w(64)= - w(83) + 1.D0/27.D0*w(91) + w(64)
      w(83)=5.D0/3.D0 + w(67)
      w(83)=z*w(83)
      w(83)= - w(105) - 4.D0/3.D0 + w(83)
      w(83)=1.D0/3.D0*w(83) + w(77)
      w(83)=w(11)*w(83)
      w(68)= - 415.D0/3.D0 + w(68)
      w(68)=z*w(68)
      w(68)= - w(72) + 89.D0/3.D0 + w(68)
      w(91)=1 + w(67)
      w(91)=z*w(91)
      w(91)=w(105) - 2 + w(91)
      w(91)=1.D0/3.D0*w(91) - w(60)
      w(76)=w(91)*w(76)
      w(68)=1.D0/9.D0*w(68) + w(76)
      w(68)=w(5)*w(68)
      w(76)=w(20)*w(54)
      w(64)= - 10.D0/3.D0*w(76) + w(68) + 1.D0/3.D0*w(64) + w(83)
      w(68)=z4*w(54)
      w(64)=4*w(64) + 35*w(68)
      w(64)=NF*w(64)
      w(68)= - NF - 1
      w(54)=w(123)*w(54)*w(68)
      w(68)= - 11.D0/9.D0 - w(67)
      w(68)=z*w(68)
      w(68)=w(77) - w(74) + 16.D0/9.D0 + w(68)
      w(68)=NF*w(68)
      w(74)= - 31.D0/3.D0 - z
      w(74)=w(74)*w(67)
      w(74)= - 53.D0/3.D0 + w(74)
      w(68)=1.D0/3.D0*w(74) + w(68)
      w(68)=w(9)*w(68)
      w(54)=8*w(68) + w(54) + 4*w(63) + w(64)
      w(63)=w(58)*z
      w(64)= - 11.D0/3.D0 - w(63)
      w(64)=w(64)*w(120)
      w(64)=157 + w(64)
      w(64)=z*w(64)
      w(64)= - w(86) - 661.D0/5.D0 + w(64)
      w(64)=1.D0/9.D0*w(64) - w(79)
      w(68)= - 13 + 74.D0/9.D0*z
      w(68)=w(68)*z
      w(68)=w(68) - w(89) + 7
      w(74)=1.D0/3.D0*w(68) - w(95)
      w(74)=NF*w(74)
      w(76)= - 2 + NF
      w(76)=w(110)*w(69)*w(76)
      w(64)=w(76) + 2*w(64) + w(74)
      w(61)=w(64)*w(61)
      w(64)= - 298.D0/3.D0 + w(82)
      w(64)=w(64)*w(67)
      w(64)=443.D0/3.D0 + w(64)
      w(57)=w(64)*w(57)
      w(63)= - 7 + w(63)
      w(63)=w(63)*w(120)
      w(63)= - 55 + w(63)
      w(63)=z*w(63)
      w(63)=w(72) + 331.D0/5.D0 + w(63)
      w(63)=w(3)*w(63)
      w(57)=w(63) + 58.D0/3.D0*w(1) - 7.D0/3.D0 + w(57)
      w(63)=w(79) - w(68)
      w(63)=w(3)*w(63)
      w(64)= - 47.D0/3.D0 + w(85)
      w(64)=z*w(64)
      w(64)= - w(72) - 43.D0/3.D0 + w(64)
      w(63)=1.D0/9.D0*w(64) + w(63)
      w(63)=NF*w(63)
      w(64)=w(69)*w(90)
      w(57)=w(63) + 1.D0/9.D0*w(57) + w(64)
      w(63)= - 2 + 5.D0/3.D0*NF
      w(63)=w(63)*w(103)
      w(57)=w(61) + 4.D0/3.D0*w(57) + w(63)
      w(57)=w(12)*w(57)
      w(58)=w(58)*w(59)
      w(58)=25 + w(58)
      w(58)=w(58)*w(67)
      w(58)= - 1453.D0/9.D0 + w(58)
      w(58)=z*w(58)
      w(59)= - 127.D0/3.D0 - 20*z
      w(59)=z*w(59)
      w(59)=7*w(60) - 37.D0/3.D0 + w(59)
      w(59)=w(3)*w(59)
      w(58)=w(59) - w(86) - 229.D0/9.D0 + w(58)
      w(58)=1.D0/9.D0*w(58) + 2*w(55)
      w(59)=w(60) - w(66)
      w(59)=w(3)*w(59)
      w(61)= - 95.D0/9.D0 + 74*z
      w(61)=z*w(61)
      w(59)=14*w(59) - w(122) + 589.D0/9.D0 + w(61)
      w(59)=1.D0/3.D0*w(59) - 10*w(55)
      w(61)=1.D0/3.D0*NF
      w(59)=w(59)*w(61)
      w(58)=2*w(58) + w(59)
      w(58)=z2*w(58)
      w(59)=251 + 64*z
      w(59)=z*w(59)
      w(59)= - 64*w(1) + 155 + w(59)
      w(59)=1.D0/27.D0*w(59) - w(77)
      w(63)=67 + w(78)
      w(63)=z*w(63)
      w(63)=28*w(1) - 89 + w(63)
      w(60)=1.D0/9.D0*w(63) - 10*w(60)
      w(60)=w(60)*w(61)
      w(59)=2*w(59) + w(60)
      w(59)=w(59)*w(111)
      w(54)=w(57) + w(59) + 1.D0/3.D0*w(54) + w(58)
      w(54)=TF*w(54)
      w(57)= - 7 - w(98)
      w(57)=w(57)*w(114)
      w(58)= - 2 + 5.D0/3.D0*z
      w(58)=w(58)*w(47)
      w(57)=w(58) - 28.D0/9.D0*w(1) + 16 + w(57)
      w(57)=w(57)*w(47)
      w(58)= - 7631 - 2306*z
      w(58)=z*w(58)
      w(58)= - 4262 + w(58)
      w(58)=2.D0/3.D0*w(58) - 853*w(1)
      w(53)=152.D0/3.D0*w(53) - 124.D0/3.D0*w(55) - 112*w(70) + 1.D0/9.D
     & 0*w(58) + w(57)
      w(53)=CA*w(53)
      w(55)=2.D0/3.D0*w(52)
      w(57)= - 127 + 31*z
      w(57)=w(57)*w(55)
      w(53)=w(53) + w(57)
      w(53)=w(53)*w(111)
      w(57)= - 3096 - 5729*z
      w(57)=z*w(57)
      w(57)= - 1145 + w(57)
      w(57)=z*w(57)
      w(57)= - 4805 + w(57)
      w(57)=z*w(57)
      w(57)= - 703 + w(57)
      w(57)=w(57)*w(52)
      w(58)= - 2819 - 20970*z
      w(58)=z*w(58)
      w(58)=10430 + w(58)
      w(58)=z*w(58)
      w(58)= - 857 + w(58)
      w(58)=z*w(58)
      w(58)=6540 + w(58)
      w(58)=CA*w(5)*w(58)
      w(57)=w(58) + w(57)
      w(57)=w(46)*w(57)
      w(48)=w(50) + 8*w(54) + w(56) + w(53) + 4*w(57) + 2*w(48) + w(49)
      w(48)=TF*w(48)
      w(49)= - 1 + 101.D0/9.D0*z
      w(49)=w(49)*w(67)
      w(50)=z*w(51)
      w(50)=8.D0/3.D0*w(1) - 1 + w(50)
      w(50)=w(3)*w(50)
      w(49)=7*w(50) + w(75) + 1 + w(49)
      w(49)=w(3)*w(49)
      w(49)= - 2.D0/3.D0*w(81) + w(49)
      w(49)=w(3)*w(49)
      w(50)=w(104)*w(118)
      w(51)=w(7)*w(84)
      w(49)=4*w(51) + w(49) + w(50)
      w(50)= - 1 + w(67)
      w(50)=z*w(50)
      w(50)=w(80) - 1 + w(50)
      w(47)=w(50)*w(47)
      w(47)= - 1.D0/9.D0*w(87) + w(47)
      w(47)=w(47)*w(102)
      w(50)=40.D0/3.D0*z
      w(51)=37 - w(50)
      w(51)=z*w(51)
      w(51)= - w(86) + 37 + w(51)
      w(51)=w(3)*w(51)
      w(51)= - w(93) + w(51)
      w(51)=w(4)*w(51)
      w(53)=56.D0/3.D0*w(96)
      w(54)= - w(6)*w(53)
      w(56)=w(9)*w(107)
      w(57)= - 13 + w(120)
      w(57)=z*w(57)
      w(57)=w(105) - 13 + w(57)
      w(57)=z3*w(57)
      w(58)=4.D0/3.D0*w(84)
      w(59)=w(10)*w(58)
      w(47)=4*w(57) - 2*w(56) + 2.D0/3.D0*w(51) + w(54) + w(59) + 1.D0/
     & 3.D0*w(49) + w(47)
      w(47)=CA*w(47)
      w(49)=w(8)*w(65)*w(62)
      w(51)=22 + 229.D0/3.D0*z
      w(51)=w(51)*w(67)
      w(51)=w(88) - 91 + w(51)
      w(54)= - 17 - 148.D0/3.D0*z
      w(54)=z*w(54)
      w(54)= - 148.D0/3.D0*w(1) - 17 + w(54)
      w(54)=w(3)*w(54)
      w(51)=2.D0/3.D0*w(51) + w(54)
      w(51)=w(51)*w(55)
      w(47)=w(51) + w(49) + w(47)
      w(49)= - 1 - w(50)
      w(49)=z*w(49)
      w(49)= - w(86) - 1 + w(49)
      w(49)=w(3)*w(49)
      w(49)=16.D0/9.D0*w(92) + w(49)
      w(49)=w(3)*w(49)
      w(50)=w(4)*w(53)
      w(51)=w(2)*w(3)*w(96)
      w(53)= - w(5)*w(58)
      w(49)= - 56.D0/9.D0*w(51) + w(50) + w(49) + w(53)
      w(49)=CA*w(49)
      w(50)=19 + w(73)
      w(50)=z*w(50)
      w(50)=w(71) + 19 + w(50)
      w(50)=w(50)*w(52)
      w(49)=8.D0/3.D0*w(50) + w(49)
      w(49)=w(2)*w(49)
      w(47)=2*w(47) + w(49)
      w(47)=w(2)*TF*w(47)
      w(47)=w(47) + w(48)
*
      APS1 = 2*CF*w(47)*as**3
*
      RETURN
      END      
      REAL*8 FUNCTION APS2(z,nf,as,LL)
*     --------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of aQqPS3: contribution of H[a,1-2*z] and H[a,z]
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(53)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      complex*16 HHc1,HHc2,HHc3,HHc4,HHc5
      real*8 HHr1,HHr2,HHr3,HHr4,HHr5
      real*8 HHi1,HHi2,HHi3,HHi4,HHi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension HHc1(-1:1),HHc2(-1:1,-1:1),HHc3(-1:1,-1:1,-1:1),
     $          HHc4(-1:1,-1:1,-1:1,-1:1),
     $          HHc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension HHr1(-1:1),HHr2(-1:1,-1:1),HHr3(-1:1,-1:1,-1:1),
     $          HHr4(-1:1,-1:1,-1:1,-1:1),
     $          HHr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension HHi1(-1:1),HHi2(-1:1,-1:1),HHi3(-1:1,-1:1,-1:1),
     $          HHi4(-1:1,-1:1,-1:1,-1:1),
     $          HHi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4,t
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      z5 = 1.0369277551433699263D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4 
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
      t=1.0D0-2.0D0*z
      nw = 5
      call hplog5(t,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
      call hplog5(z,nw,HHc1,HHc2,HHc3,HHc4,HHc5,
     $                       HHr1,HHr2,HHr3,HHr4,HHr5,
     $                       HHi1,HHi2,HHi3,HHi4,HHi5,-1,1)
*
      w(1)=Hr2(0,-1)
      w(2)=z**(-1)
      w(3)=HHr1(0)
      w(4)=HHr2(0,1)
      w(5)=HHr1(1)
      w(6)=Hr3(0,-1,-1)
      w(7)=Hr4(0,-1,-1,-1)
      w(8)=Hr5(0,-1,-1,-1,1)
      w(9)=Hr4(0,-1,-1,1)
      w(10)=Hr5(0,-1,-1,1,-1)
      w(11)=Hr5(0,-1,-1,1,1)
      w(12)=Hr3(0,-1,1)
      w(13)=Hr4(0,-1,1,-1)
      w(14)=Hr5(0,-1,1,-1,-1)
      w(15)=Hr4(0,-1,1,1)
      w(16)=Hr5(0,-1,1,1,-1)
      w(17)=Hr5(0,-1,1,1,1)
      w(18)=Hr2(0,1)
      w(19)=Hr3(0,1,-1)
      w(20)=Hr4(0,1,-1,-1)
      w(21)=Hr5(0,1,-1,-1,1)
      w(22)=Hr4(0,1,-1,1)
      w(23)=Hr5(0,1,-1,1,-1)
      w(24)=Hr5(0,1,-1,1,1)
      w(25)=Hr3(0,1,1)
      w(26)=Hr4(0,1,1,-1)
      w(27)=Hr5(0,1,1,-1,-1)
      w(28)=Hr4(0,1,1,1)
      w(29)=Hr5(0,1,1,1,-1)
      w(30)=Hr5(0,1,1,1,1)
      w(31)=4.D0/3.D0*z
      w(32)=w(31) - 15
      w(32)=w(32)*z
      w(33)=z + 1
      w(34)=w(33)*w(3)
      w(35)=2*w(34)
      w(36)=4.D0/3.D0*w(2)
      w(32)= - w(32) - w(36) + w(35) + 5
      w(32)=w(32)*w(26)
      w(37)=4*z
      w(38)=w(37) - 5
      w(38)=w(38)*z
      w(39)=8.D0/3.D0*w(2)
      w(40)=6*w(34)
      w(38)= - w(40) + w(38) - w(39) - 5
      w(38)=w(38)*w(20)
      w(41)=w(37) + 3
      w(41)=w(41)*z
      w(41)= - w(40) + w(41) - 3 - 4*w(2)
      w(41)=w(41)*w(7)
      w(42)=w(29) + w(16)
      w(43)=6*w(33)
      w(42)=w(43)*w(42)
      w(43)=2*w(33)
      w(44)=w(43)*w(14)
      w(45)=4*w(33)
      w(46)=w(45)*w(10)
      w(32)= - w(41) + w(32) - w(38) + w(42) + w(44) + w(46)
      w(38)=w(31) + 1
      w(38)=w(38)*z
      w(38)=w(38) - w(36) - 1
      w(41)=4*w(5)
      w(41)=w(38)*w(41)
      w(42)=w(31) - 7
      w(42)=w(42)*z
      w(42)=w(42) - 3
      w(34)=w(34) - w(42)
      w(44)=4*w(3)
      w(34)=w(34)*w(44)
      w(46)= - 8*w(4) + 4*z2
      w(46)=w(33)*w(46)
      w(47)=w(37) + 25
      w(47)=w(47)*z
      w(47)=w(47) - 13
      w(47)=w(2) + 8*w(47)
      w(34)=w(46) + w(41) - w(34) + 1.D0/3.D0*w(47)
      w(41)=w(34)*w(19)
      w(46)=24*w(24)
      w(46)=w(46)*w(33)
      w(47)=w(33)*w(17)
      w(32)=w(46) + w(41) + 48*w(47) - 4*w(32)
      w(41)=CF**2
      w(46)= - w(41)*w(32)
      w(47)=w(31) - 5
      w(48)=2*z
      w(47)=w(47)*w(48)
      w(44)=w(44)*w(33)
      w(49)=2.D0/3.D0*w(2)
      w(47)= - w(44) + w(47) - w(49) - 5
      w(47)=w(12)*w(47)
      w(50)=w(31) - 9
      w(50)=w(50)*w(48)
      w(50)=w(50) - w(44) + w(49) - 7
      w(50)=w(25)*w(50)
      w(47)=w(47) + w(50)
      w(50)=w(1) + w(18)
      w(51)=w(50)*w(34)
      w(52)=w(31) - 1
      w(52)=w(52)*w(48)
      w(52)= - w(44) + w(52) - 3 - 2*w(2)
      w(53)=4*w(19)
      w(52)=w(52)*w(53)
      w(53)=w(31) + 3
      w(48)=w(53)*w(48)
      w(44)= - w(44) + w(48) - 1 - 10.D0/3.D0*w(2)
      w(44)=w(6)*w(44)
      w(48)=w(20) + w(7)
      w(53)=12*w(33)
      w(48)=w(53)*w(48)
      w(45)=w(45)*w(26)
      w(44)=4*w(47) + 4*w(44) + w(52) + w(48) + w(51) - w(45)
      w(45)=CA*CF
      w(41)= - w(45) + 2*w(41)
      w(44)=w(41)*w(44)
      w(47)=w(28) + w(15)
      w(48)=w(9) + w(22)
      w(51)=4*w(13)
      w(47)= - w(51) - 20*w(48) - 36*w(47)
      w(33)=w(33)*w(41)
      w(47)=w(33)*w(47)
      w(48)= - w(37) + w(49) - 1
      w(48)=w(50)*w(48)
      w(43)=w(43)*w(19)
      w(43)=w(48) - w(43)
      w(43)=w(43)*w(41)
      w(48)=2*w(33)
      w(49)=w(25) + w(12)
      w(50)= - w(6) - w(49)
      w(48)=w(48)*w(50)
      w(43)=w(48) + w(43)
      w(43)=ln2*w(43)
      w(43)=8*w(43) + w(47) + w(44)
      w(43)=ln2*w(43)
      w(43)=w(46) + w(43)
      w(32)=w(32)*w(45)
      w(44)=w(37) - 13
      w(44)=w(44)*z
      w(36)= - w(36) + w(44) - w(40) - 7
      w(36)=w(15)*w(36)
      w(31)=w(31) + 9
      w(31)=w(31)*z
      w(31)= - w(31) + w(35) + w(39) - 1
      w(31)= - w(9)*w(31)
      w(37)=w(37) - 21
      w(37)=w(37)*z
      w(37)=w(37) - w(40) - 9
      w(37)=w(28)*w(37)
      w(38)=w(38) - w(35)
      w(38)=w(22)*w(38)
      w(31)=w(36) + w(37) + w(38) + w(31)
      w(36)= - w(6) + w(49)
      w(34)=w(36)*w(34)
      w(35)=w(35) - w(42)
      w(35)=w(51)*w(35)
      w(31)=4*w(31) + w(35) + w(34)
      w(31)=w(41)*w(31)
      w(34)=w(21) + w(8) - w(11)
      w(34)=8*w(27) + 16*w(23) - 48*w(30) + 24*w(34)
      w(33)=w(33)*w(34)
      w(31)=w(32) + w(33) + 2*w(43) + w(31)
*
      APS2 = 8*TF*w(31)*as**3
*
      RETURN
      END      
      REAL*8 FUNCTION GGREG(z,nf,as,LL)
*     ---------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of the OME AggQ without aggQ3 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(152)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half,G
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
      w(1)=z**(-1)
      w(2)=Hr1(-1)
      w(3)=Hr1(0)
      w(4)=Hr2(0,-1)
      w(5)=Hr2(0,1)
      w(6)=Hr3(0,-1,-1)
      w(7)=Hr3(0,0,-1)
      w(8)=Hr3(0,0,1)
      w(9)=Hr1(1)
      w(10)=Hr3(0,1,1)
      w(11)=Hr4(0,-1,-1,-1)
      w(12)=Hr4(0,-1,0,1)
      w(13)=Hr3(0,-1,1)
      w(14)=Hr4(0,0,0,-1)
      w(15)=Hr4(0,0,0,1)
      w(16)=Hr4(0,0,1,1)
      w(17)=Hr3(0,1,-1)
      w(18)=Hr4(0,1,1,1)
      w(19)=Hr5(0,-1,0,-1,-1)
      w(20)=Hr4(0,0,-1,-1)
      w(21)=Hr5(0,0,-1,-1,-1)
      w(22)=Hr5(0,0,-1,0,-1)
      w(23)=Hr5(0,0,-1,0,1)
      w(24)=Hr4(0,0,-1,1)
      w(25)=Hr5(0,0,0,-1,-1)
      w(26)=Hr5(0,0,0,-1,1)
      w(27)=Hr5(0,0,0,0,-1)
      w(28)=Hr5(0,0,0,0,1)
      w(29)=Hr5(0,0,0,1,-1)
      w(30)=Hr5(0,0,0,1,1)
      w(31)=Hr4(0,0,1,-1)
      w(32)=Hr5(0,0,1,0,-1)
      w(33)=Hr5(0,0,1,0,1)
      w(34)=Hr5(0,0,1,1,1)
      w(35)=Hr5(0,1,0,1,1)
      w(36)=Hr5(0,1,1,1,1)
      w(37)=1/( - 81 + 81*z)
      w(38)=1/( - 27 + 27*z)
      w(39)=1/( - 9 + 9*z**2)
      w(40)=1/( - 3 + 3*z)
      w(41)=1/( - 3 + 3*z**2)
      w(42)=1/(3 + 3*z)
      w(43)=1/( - 27*z + 27*z**2)
      w(44)=1/( - 9*z + 9*z**2)
      w(45)=1/( - 3*z + 3*z**2)
      w(46)=1/(3*z + 3*z**2)
      w(47)=1/(9*z + 9*z**2)
      w(48)=4.D0/3.D0*z
      w(49)=w(48) + 1
      w(50)=w(49)*z
      w(51)=4.D0/3.D0*w(1)
      w(52)=w(50) - w(51) - 1
      w(53)=2*w(3)
      w(54)=w(52)*w(53)
      w(55)=w(52)*w(9)
      w(56)=8.D0/3.D0*w(1)
      w(57)=w(56) + 9
      w(58)=8.D0/3.D0*z
      w(59)=9 + w(58)
      w(59)=z*w(59)
      w(59)=w(59) - w(57)
      w(60)=8*LL
      w(59)=w(59)*w(60)
      w(61)= - 1 - 154.D0/9.D0*z
      w(61)=z*w(61)
      w(59)= - w(55) + w(59) + w(54) + 46.D0/9.D0*w(1) + 13 + w(61)
      w(59)=w(9)*w(59)
      w(61)=z - 1
      w(62)=w(61)*w(53)
      w(63)=4*w(1)
      w(64)=w(63) + 11
      w(65)=4*z
      w(66)=17 + w(65)
      w(66)=z*w(66)
      w(66)=w(66) + w(64)
      w(67)=z + 1
      w(68)=w(67)*LL
      w(69)=2.D0/3.D0*w(68)
      w(66)=w(69) + 1.D0/3.D0*w(66) + w(62)
      w(66)=LL*w(66)
      w(70)=2.D0/3.D0*w(1)
      w(71)=16.D0/3.D0 + z
      w(71)=z*w(71)
      w(71)= - w(70) - 2.D0/3.D0 + w(71)
      w(72)=w(67)*w(3)
      w(73)=3*w(72)
      w(71)=2*w(71) - w(73)
      w(71)=w(71)*w(53)
      w(74)=10*z
      w(75)=539.D0/3.D0 + w(74)
      w(75)=z*w(75)
      w(75)= - 56.D0/3.D0*w(1) - 334.D0/3.D0 + w(75)
      w(66)=w(66) + 1.D0/3.D0*w(75) + w(71)
      w(71)=4*LL
      w(66)=w(66)*w(71)
      w(75)=w(65) - 3
      w(75)=w(75)*z
      w(75)=w(75) + w(63) - 3
      w(76)=w(61)*LL
      w(77)=w(61)*w(3)
      w(78)=6*w(77)
      w(79)=28*w(76) - w(78) + w(75)
      w(79)=w(4)*w(79)
      w(80)=w(53)*w(67)
      w(81)=23*z
      w(82)=14 + w(81)
      w(82)= - 32*w(68) + 1.D0/3.D0*w(82) - w(80)
      w(83)=2*w(5)
      w(82)=w(82)*w(83)
      w(84)= - 17 - z
      w(84)=w(84)*w(65)
      w(84)= - 152.D0/5.D0*w(1) + 137.D0/5.D0 + w(84)
      w(85)=1 - 38.D0/5.D0*z
      w(85)=w(85)*w(53)
      w(86)=73 + 89*z
      w(86)=LL*w(86)
      w(84)=2.D0/5.D0*w(86) + 1.D0/3.D0*w(84) + w(85)
      w(85)=2*z2
      w(84)=w(84)*w(85)
      w(86)=16*w(6)
      w(87)= - w(86) + 24*w(7)
      w(87)=w(61)*w(87)
      w(88)=32*z
      w(89)=5 + w(88)
      w(89)=1.D0/3.D0*w(89) - w(77)
      w(89)=w(89)*w(53)
      w(90)=121 + 298*z
      w(90)=z*w(90)
      w(90)=80*w(1) + 427 + w(90)
      w(89)=1.D0/9.D0*w(90) + w(89)
      w(89)=w(3)*w(89)
      w(90)=4*w(10)
      w(91)=w(67)*w(90)
      w(92)=16 - 845.D0/27.D0*z
      w(92)=z*w(92)
      w(92)=182.D0/27.D0*w(1) - 85.D0/9.D0 + w(92)
      w(93)= - 37 - 43*z
      w(93)=z3*w(93)
      w(59)=w(84) + w(91) + w(82) + 2.D0/3.D0*w(93) + w(59) + 2*w(79)
     &  + w(66) + 4*w(92) + w(89) + w(87)
      w(59)=w(59)*w(85)
      w(66)=4*w(4)
      w(79)= - w(85) + w(66)
      w(82)=w(48) - 1
      w(82)=w(82)*z
      w(84)=w(51) - 1
      w(87)=w(82) + w(84)
      w(79)=w(87)*w(79)
      w(89)=19.D0/3.D0*z
      w(91)=w(89) + 1
      w(91)=w(91)*z
      w(92)=10.D0/3.D0*w(1)
      w(93)=w(92) - 2
      w(91)=w(91) + w(93)
      w(94)=4.D0/3.D0*w(91)
      w(95)=w(87)*w(3)
      w(96)=w(94) + w(95)
      w(96)=w(3)*w(96)
      w(97)=w(65) - 17
      w(97)=w(97)*z
      w(97)=w(97) + w(63) - 17
      w(98)=w(53)*LL
      w(99)=w(97)*w(98)
      w(100)=w(2)*w(95)
      w(79)= - 4.D0/3.D0*w(100) + w(96) + w(99) + w(79)
      w(79)=w(2)*w(79)
      w(96)=2.D0/3.D0*w(91)
      w(99)= - LL*w(97)
      w(99)=w(99) - w(96) - w(95)
      w(99)=w(99)*w(66)
      w(100)= - w(3)*w(75)
      w(101)= - 25 + w(48)
      w(101)=z*w(101)
      w(101)=w(51) - 25 + w(101)
      w(102)=2*LL
      w(101)=w(101)*w(102)
      w(94)=w(101) + w(94) + w(100)
      w(94)=z2*w(94)
      w(91)= - 2*w(91) + w(95)
      w(91)=w(3)*w(91)
      w(100)=10 + 91.D0/3.D0*z
      w(101)=2*z
      w(100)=w(100)*w(101)
      w(100)=w(100) - 25 + 47.D0/3.D0*w(1)
      w(91)= - 4.D0/3.D0*w(100) + w(91)
      w(103)=1.D0/3.D0*w(3)
      w(91)=w(91)*w(103)
      w(104)=17 + w(58)
      w(104)=z*w(104)
      w(104)=w(56) + 17 + w(104)
      w(104)=w(3)*w(104)
      w(105)=32 + 29.D0/9.D0*z
      w(105)=w(105)*z
      w(105)=w(105) + 25 - 34.D0/9.D0*w(1)
      w(104)=4*w(105) + w(104)
      w(104)=w(3)*w(104)
      w(106)=w(71)*w(95)
      w(104)=w(104) + w(106)
      w(104)=LL*w(104)
      w(106)=1.D0/3.D0*z
      w(107)=w(106) + 1
      w(107)=w(107)*z
      w(108)=1.D0/3.D0*w(1)
      w(109)=w(108) + 1
      w(110)=w(107) + w(109)
      w(111)=w(110)*w(71)
      w(95)=w(111) + w(95)
      w(111)=4*w(5)
      w(112)=w(95)*w(111)
      w(113)= - 8*w(6) + 4*w(7)
      w(113)=w(87)*w(113)
      w(114)=2*z3
      w(75)=w(75)*w(114)
      w(75)=w(79) + w(94) + w(112) + w(75) + w(99) + w(91) + w(104) + 
     & w(113)
      w(79)=4*w(2)
      w(75)=w(75)*w(79)
      w(91)=w(101) - 1
      w(94)=w(91)*w(3)
      w(99)=w(94) - 3*w(61) - w(56)
      w(99)=w(99)*w(53)
      w(104)=w(102)*w(67)
      w(112)= - 49 + 16*w(1)
      w(113)=20*z
      w(115)= - 53 - w(113)
      w(115)=z*w(115)
      w(115)=w(115) - w(112)
      w(116)=z + 2
      w(117)=w(3)*w(116)
      w(115)= - w(104) + 1.D0/3.D0*w(115) - 12*w(117)
      w(115)=LL*w(115)
      w(117)=5*w(55)
      w(118)=18 + 49.D0/9.D0*z
      w(118)=w(118)*w(101)
      w(79)= - w(87)*w(79)
      w(119)=w(67)*w(5)
      w(120)=w(61)*w(4)
      w(121)=z2*w(67)
      w(79)=w(79) + w(121) + 10*w(119) - w(117) - 8*w(120) + w(115) + 
     & w(99) + 149.D0/9.D0*w(1) + 23.D0/3.D0 + w(118)
      w(79)=w(8)*w(79)
      w(99)=10.D0/3.D0*z
      w(115)=3 + w(99)
      w(115)=w(115)*w(101)
      w(115)= - w(104) + 8*w(72) - 20.D0/3.D0*w(1) - 5 + w(115)
      w(115)=w(16)*w(115)
      w(118)=6*LL
      w(121)=3*z
      w(122)= - 5 - w(121)
      w(122)=w(122)*w(118)
      w(122)=w(122) + w(78) - w(87)
      w(122)=w(14)*w(122)
      w(79)=w(79) + w(115) + w(122)
      w(115)=1.D0/3.D0*w(68)
      w(74)= - 7 - w(74)
      w(74)= - w(115) + 1.D0/3.D0*w(74) + w(72)
      w(74)=w(74)*w(102)
      w(122)=5*z
      w(123)=17*z
      w(124)= - 19 + w(123)
      w(124)=w(124)*w(122)
      w(124)= - 64*w(1) + 403 + w(124)
      w(125)=7*w(72)
      w(126)=35.D0/3.D0 + w(65)
      w(126)=z*w(126)
      w(126)=w(125) - 40.D0/3.D0 + w(126)
      w(126)=w(3)*w(126)
      w(74)=w(74) + 1.D0/9.D0*w(124) + w(126)
      w(74)=w(74)*w(102)
      w(124)=w(65) + 3
      w(124)=w(124)*z
      w(124)=w(124) - w(63) - 3
      w(126)=w(3)*w(124)
      w(127)=w(52)*w(102)
      w(126)=w(126) + w(127)
      w(127)=2*w(9)
      w(126)=w(126)*w(127)
      w(128)= - w(104) - w(73) + w(52)
      w(83)=w(128)*w(83)
      w(128)=245 - 1732.D0/9.D0*z
      w(128)=w(128)*w(106)
      w(129)=w(56) - w(67)
      w(129)=w(3)*w(129)
      w(99)= - 83 - w(99)
      w(99)=z*w(99)
      w(99)= - 109.D0/3.D0*w(1) + 16 + w(99)
      w(99)=2.D0/3.D0*w(99) + w(129)
      w(99)=w(3)*w(99)
      w(129)=w(67)*z3
      w(74)=w(83) - 20.D0/3.D0*w(129) + w(126) + w(74) + w(99) + 908.D0/
     & 27.D0*w(1) - 36 + w(128)
      w(74)=w(74)*w(111)
      w(83)=16.D0/3.D0*w(1)
      w(99)=32.D0/3.D0*z
      w(126)= - 3 - w(99)
      w(126)=z*w(126)
      w(126)= - 3*w(77) - w(83) + 9 + w(126)
      w(126)=w(3)*w(126)
      w(128)=w(62) - w(87)
      w(130)=w(128)*w(102)
      w(105)=w(130) - 2*w(105) + w(126)
      w(105)=w(105)*w(102)
      w(126)= - 8 + w(89)
      w(126)=z*w(126)
      w(126)=w(126) + w(93)
      w(130)=2.D0/3.D0*w(77) - w(87)
      w(130)=w(3)*w(130)
      w(126)=4.D0/3.D0*w(126) + w(130)
      w(126)=w(3)*w(126)
      w(76)= - w(77) - 7*w(76)
      w(76)=w(76)*w(66)
      w(76)=w(76) + w(105) + 4.D0/9.D0*w(100) + w(126)
      w(76)=w(76)*w(66)
      w(100)=w(52)*w(3)
      w(105)= - 1 + w(89)
      w(105)=z*w(105)
      w(92)= - w(92) - 2 + w(105)
      w(92)=4.D0/3.D0*w(92) - w(100)
      w(92)=w(3)*w(92)
      w(105)=2*w(1)
      w(126)=1 + w(113)
      w(126)=z*w(126)
      w(126)= - w(105) - 19 + w(126)
      w(124)= - w(124)*w(53)
      w(124)=1.D0/3.D0*w(126) + w(124)
      w(124)=LL*w(124)
      w(126)=5 - 328.D0/3.D0*z
      w(126)=z*w(126)
      w(126)=67.D0/3.D0*w(1) + 82 + w(126)
      w(92)=w(124) + 1.D0/9.D0*w(126) + w(92)
      w(124)=w(52)*LL
      w(126)=1 - 7.D0/3.D0*z
      w(126)=w(126)*w(65)
      w(126)=w(108) + 5 + w(126)
      w(126)= - w(124) + 1.D0/3.D0*w(126) + w(54)
      w(126)=4*w(126) - w(55)
      w(130)=1.D0/3.D0*w(9)
      w(126)=w(126)*w(130)
      w(92)=2*w(92) + w(126)
      w(92)=w(9)*w(92)
      w(126)=52.D0/3.D0*w(1)
      w(131)=44 - 115.D0/3.D0*z
      w(131)=z*w(131)
      w(131)=w(126) - 23 + w(131)
      w(132)=16.D0/3.D0*z
      w(133)= - 13 - w(132)
      w(133)=z*w(133)
      w(83)=w(83) + 13 + w(133)
      w(83)=w(3)*w(83)
      w(83)=2.D0/3.D0*w(131) + w(83)
      w(83)=w(3)*w(83)
      w(131)=2.D0/3.D0*w(124)
      w(133)=5 + 104.D0/3.D0*z
      w(133)=z*w(133)
      w(133)= - 32.D0/3.D0*w(1) - 29 + w(133)
      w(133)=w(131) + 1.D0/3.D0*w(133) - w(54)
      w(133)=LL*w(133)
      w(134)= - 449 + 340.D0/3.D0*z
      w(134)=z*w(134)
      w(134)=254.D0/3.D0*w(1) + 251 + w(134)
      w(83)=w(133) + 2.D0/9.D0*w(134) + w(83)
      w(83)=w(83)*w(102)
      w(133)=13*z
      w(134)=29 - w(133)
      w(134)=w(134)*w(101)
      w(135)=23*w(1)
      w(134)=w(135) - 55 + w(134)
      w(134)=w(3)*w(134)
      w(136)= - 178 + 461.D0/3.D0*z
      w(136)=z*w(136)
      w(136)= - 227.D0/3.D0*w(1) + 100 + w(136)
      w(134)=4.D0/3.D0*w(136) + w(134)
      w(134)=w(134)*w(53)
      w(136)=23 - 2500.D0/3.D0*z
      w(136)=z*w(136)
      w(136)=1771.D0/3.D0*w(1) + 328 + w(136)
      w(134)=1.D0/9.D0*w(136) + w(134)
      w(83)=1.D0/3.D0*w(134) + w(83)
      w(83)=2*w(83) + w(92)
      w(83)=w(9)*w(83)
      w(92)=2*w(55)
      w(134)=7*z
      w(58)=5 - w(58)
      w(58)=w(58)*w(134)
      w(58)=40.D0/3.D0*w(1) + 22 + w(58)
      w(136)= - 5*w(52) + w(72)
      w(136)=w(136)*w(53)
      w(137)=8*w(1)
      w(138)=8*z
      w(139)= - 23 - w(138)
      w(139)=z*w(139)
      w(139)=w(137) - 5 + w(139)
      w(139)=1.D0/3.D0*w(139) + 6*w(72)
      w(139)=w(139)*w(102)
      w(58)= - w(92) + w(139) + 1.D0/3.D0*w(58) + w(136)
      w(58)=w(58)*w(90)
      w(136)=2.D0/3.D0*z
      w(139)=w(136) + 1
      w(139)=w(139)*z
      w(140)=w(68) - w(80) - w(70) - 1 + w(139)
      w(140)=w(18)*w(140)
      w(62)=w(62) + w(87)
      w(62)=w(12)*w(62)
      w(62)=w(140) + w(62)
      w(140)=40*w(1)
      w(141)= - 103 - 184*z
      w(141)=z*w(141)
      w(141)= - w(140) - 205 + w(141)
      w(142)=4*w(3)
      w(143)=w(61)*w(142)
      w(144)=38*z
      w(143)=w(143) - 5 - w(144)
      w(143)=w(3)*w(143)
      w(141)=1.D0/3.D0*w(141) + w(143)
      w(141)=w(141)*w(103)
      w(143)= - 148 + 517.D0/3.D0*z
      w(143)=w(143)*w(106)
      w(64)= - w(138) - w(64)
      w(64)=1.D0/3.D0*w(64) - w(94)
      w(64)=w(3)*w(64)
      w(94)= - 5 + 44.D0/9.D0*z
      w(94)=z*w(94)
      w(64)=2.D0/3.D0*w(64) - 44.D0/9.D0*w(1) + 5 + w(94)
      w(64)=LL*w(64)
      w(64)=w(64) + w(141) - 127.D0/9.D0*w(1) + 5 + w(143)
      w(64)=w(64)*w(102)
      w(94)= - 4795 - 2384*z
      w(94)=z*w(94)
      w(94)= - 2230 + w(94)
      w(94)=1.D0/3.D0*w(94) + w(140)
      w(140)= - 23 + 37*z
      w(141)=11*z
      w(143)= - 7 + w(141)
      w(143)=w(3)*w(143)
      w(140)=4.D0/3.D0*w(140) + w(143)
      w(140)=w(140)*w(103)
      w(143)=26*z
      w(145)= - 1487.D0/9.D0 - w(143)
      w(145)=z*w(145)
      w(140)=w(140) - 1031.D0/9.D0 + w(145)
      w(140)=w(3)*w(140)
      w(94)=2.D0/9.D0*w(94) + w(140)
      w(94)=w(3)*w(94)
      w(140)= - 388 + 2021.D0/3.D0*z
      w(140)=z*w(140)
      w(140)= - 1313.D0/3.D0 + w(140)
      w(140)=1.D0/3.D0*w(140) + 49*w(1)
      w(64)=w(64) + 2*w(140) + w(94)
      w(64)=w(64)*w(102)
      w(94)=2 - w(122)
      w(94)=w(94)*w(142)
      w(140)=7 + z
      w(140)=w(140)*w(118)
      w(145)= - 5 + w(134)
      w(94)=w(140) + w(94) + 3*w(145) + w(137)
      w(140)=8*w(15)
      w(94)=w(94)*w(140)
      w(145)=77 - w(101)
      w(145)=w(145)*w(101)
      w(112)=w(145) - w(112)
      w(145)= - 3 - w(101)
      w(145)=w(145)*w(53)
      w(91)=w(91)*w(71)
      w(91)=w(91) + 1.D0/3.D0*w(112) + w(145)
      w(91)=w(91)*w(102)
      w(112)= - 83 - 208.D0/3.D0*z
      w(112)=w(112)*w(101)
      w(145)=w(63) + 83 - 145*z
      w(146)= - 13 + 14*z
      w(146)=w(3)*w(146)
      w(145)=1.D0/3.D0*w(145) + w(146)
      w(145)=w(145)*w(53)
      w(112)=w(145) - 76.D0/3.D0*w(1) - 185 + w(112)
      w(91)=10.D0/3.D0*w(55) + 12*w(120) + 1.D0/3.D0*w(112) + w(91)
      w(112)=4*z3
      w(91)=w(91)*w(112)
      w(145)=1 + w(121)
      w(145)=w(145)*w(142)
      w(146)= - w(61)*w(60)
      w(147)= - 11 + 56.D0/3.D0*z
      w(147)=z*w(147)
      w(145)=w(146) + w(145) + w(137) - 35 + w(147)
      w(145)=LL*w(145)
      w(146)= - w(3)*w(128)
      w(147)=17 - w(89)
      w(147)=z*w(147)
      w(93)=w(147) - w(93)
      w(147)=4*w(120)
      w(93)=w(147) + w(145) + 2.D0/3.D0*w(93) + w(146)
      w(145)=8*w(7)
      w(93)=w(93)*w(145)
      w(146)=w(77) + w(87)
      w(146)=w(3)*w(146)
      w(97)=14*w(77) + w(97)
      w(97)=LL*w(97)
      w(96)=w(147) + w(97) + w(96) + w(146)
      w(86)=w(96)*w(86)
      w(96)=w(13) + w(17)
      w(96)=16*w(96)
      w(95)= - w(95)*w(96)
      w(97)=w(11)*w(128)
      w(128)=w(31) + w(24)
      w(128)= - 16*w(20) + 32*w(128)
      w(87)=w(87)*w(128)
      w(128)=16*w(36)
      w(146)= - 144*w(33) + w(128) + 64*w(34) + 16*w(35) - 432*w(30)
      w(146)=w(67)*w(146)
      w(147)=w(27) - w(32)
      w(147)= - 96*w(25) - 32*w(22) + 192*w(29) - 256*w(21) - 128*w(19)
     &  - 64*w(147)
      w(147)=w(61)*w(147)
      w(143)=5 + w(143)
      w(143)=1.D0/3.D0*w(143) - 2.D0/5.D0*w(77)
      w(143)=w(3)*w(143)
      w(148)=20.D0/3.D0*z
      w(149)=47 + w(148)
      w(149)=z*w(149)
      w(149)=58 + w(149)
      w(143)=2.D0/3.D0*w(149) + w(143)
      w(143)=w(3)*w(143)
      w(149)=901.D0/3.D0 - 100*z
      w(149)=z*w(149)
      w(149)=169.D0/3.D0 + w(149)
      w(143)=2*w(149) + w(143)
      w(143)=w(3)*w(143)
      w(149)=3193 + 36604.D0/9.D0*z
      w(149)=z*w(149)
      w(149)=2624.D0/9.D0*w(1) + 9020.D0/3.D0 + w(149)
      w(143)=2.D0/3.D0*w(149) + w(143)
      w(143)=w(143)*w(103)
      w(149)=96173 - 417275.D0/3.D0*z
      w(149)=z*w(149)
      w(149)=102833.D0/3.D0*w(1) + 6544 + w(149)
      w(150)=w(65) - 1
      w(151)=w(28)*w(150)
      w(152)=3 - w(134)
      w(152)=z5*w(152)
      w(58)=w(75) + 40*w(152) + w(59) + w(86) + w(58) + w(93) - 32*
     & w(97) + w(74) + w(91) + w(83) + w(76) + w(94) + w(64) + w(143)
     &  + 1.D0/81.D0*w(149) + 64*w(151) + w(147) + w(146) + w(87) + 
     & w(95) + 16*w(62) + 8*w(79)
      w(59)=CA*TF
      w(58)=w(58)*w(59)
      w(62)=w(52)*w(142)
      w(64)=40*z
      w(74)= - 7 + w(64)
      w(74)=z*w(74)
      w(74)= - w(63) - 29 + w(74)
      w(75)=w(107) - w(109)
      w(76)=LL*w(75)
      w(74)=w(117) - 32*w(76) + 1.D0/3.D0*w(74) + w(62)
      w(74)=w(74)*w(127)
      w(76)=w(121) + 2
      w(79)= - w(115) - w(72) - w(76)
      w(79)=w(79)*w(102)
      w(83)= - 64 + w(134)
      w(83)=w(83)*w(136)
      w(86)= - w(72) - 3 - w(138)
      w(86)=w(86)*w(53)
      w(79)=w(79) + w(86) + 7 + w(83)
      w(79)=w(79)*w(60)
      w(83)=4*w(72)
      w(86)=w(83) - w(56)
      w(87)=16*w(68)
      w(91)=5 + w(65)
      w(91)=z*w(91)
      w(91)=w(87) + 3 + w(91) - w(86)
      w(91)=w(91)*w(111)
      w(93)=9 - 188.D0/3.D0*z
      w(93)=z*w(93)
      w(94)=w(133) + 23
      w(95)=LL*w(94)
      w(93)= - 4.D0/5.D0*w(95) + 56.D0/5.D0*w(72) + 7 + 1.D0/5.D0*w(93)
      w(95)=4*z2
      w(93)=w(93)*w(95)
      w(97)= - 1 + w(82)
      w(97)=2*w(97) - w(72)
      w(97)=w(97)*w(53)
      w(107)=469 + 80*z
      w(107)=z*w(107)
      w(107)=221 + w(107)
      w(97)=1.D0/3.D0*w(107) + w(97)
      w(97)=w(3)*w(97)
      w(107)=28*z
      w(109)= - 79.D0/3.D0 - w(107)
      w(109)=z*w(109)
      w(109)=20*w(1) - 25 + w(109)
      w(117)=LL*w(120)
      w(143)=w(10)*w(67)
      w(74)=w(93) - 40*w(143) + w(91) + 40.D0/3.D0*w(129) + w(74) - 64*
     & w(117) + w(79) + 2*w(109) + w(97)
      w(74)=z2*w(74)
      w(79)=7 + w(48)
      w(79)=z*w(79)
      w(69)=w(69) + w(80) - w(51) + 3 + w(79)
      w(69)=w(69)*w(102)
      w(79)=2 + w(108)
      w(79)=2*w(79) - w(72)
      w(91)=8*w(3)
      w(79)=w(79)*w(91)
      w(93)= - 23 + w(138)
      w(93)=z*w(93)
      w(93)= - 43 + w(93)
      w(69)=w(69) + w(79) + 1.D0/3.D0*w(93) - w(137)
      w(69)=LL*w(69)
      w(79)=2.D0/3.D0*w(72)
      w(93)= - 9 - w(48)
      w(93)=z*w(93)
      w(84)= - w(79) + w(93) + w(84)
      w(84)=w(3)*w(84)
      w(93)=w(101) - 7
      w(97)= - w(93)*w(122)
      w(97)= - w(105) + 49 + w(97)
      w(84)=2.D0/3.D0*w(97) + w(84)
      w(84)=w(3)*w(84)
      w(97)= - w(124) + 1.D0/3.D0*w(61) - w(100)
      w(97)=2*w(97) - w(55)
      w(97)=w(97)*w(127)
      w(109)= - w(139) + w(108)
      w(109)=w(68) + 2*w(109) + w(72)
      w(109)=w(109)*w(111)
      w(117)=17 + w(134)
      w(117)=w(117)*w(138)
      w(117)= - 121 + w(117)
      w(143)=6*w(1)
      w(69)=w(109) + 44.D0/3.D0*w(129) + w(97) + w(69) + w(84) + 1.D0/3.
     & D0*w(117) + w(143)
      w(69)=w(69)*w(111)
      w(84)=w(106) - 1
      w(97)= - w(84)*w(65)
      w(97)= - w(68) + 9*w(72) - w(51) - 9 + w(97)
      w(97)=w(97)*w(71)
      w(109)=11 - w(48)
      w(109)=z*w(109)
      w(73)=w(73) - w(51) + 5 + w(109)
      w(73)=w(73)*w(53)
      w(109)=8*w(55)
      w(113)= - 239 + w(113)
      w(113)=z*w(113)
      w(113)=w(63) - 187 + w(113)
      w(117)= - w(67)*w(95)
      w(73)=w(117) - 16*w(119) + w(109) + w(97) + 1.D0/3.D0*w(113) + 
     & w(73)
      w(97)=4*w(8)
      w(73)=w(73)*w(97)
      w(113)= - 49 - w(88)
      w(113)=w(113)*w(106)
      w(117)=3 - w(48)
      w(117)=z*w(117)
      w(117)=w(79) + 1 + w(117)
      w(117)=w(3)*w(117)
      w(129)= - 1 - w(88)
      w(129)=z*w(129)
      w(129)=31 + w(129)
      w(117)=1.D0/3.D0*w(129) + w(117)
      w(117)=w(117)*w(53)
      w(129)= - w(49)*w(101)
      w(129)=w(129) + w(72)
      w(129)=w(129)*w(53)
      w(129)= - 23.D0/3.D0*w(61) + w(129)
      w(146)=1.D0/3.D0*LL
      w(129)=w(129)*w(146)
      w(113)=w(129) + w(117) + w(63) + 23 + w(113)
      w(113)=w(113)*w(102)
      w(117)=128.D0/3.D0 - w(133)
      w(117)=w(117)*w(101)
      w(117)=45 + w(117)
      w(129)=2 + w(122)
      w(129)=8.D0/3.D0*w(129) + w(72)
      w(129)=w(3)*w(129)
      w(147)=599.D0/3.D0 - 24*z
      w(147)=z*w(147)
      w(129)=w(129) + 69 + w(147)
      w(129)=w(3)*w(129)
      w(117)=4*w(117) + w(129)
      w(117)=w(3)*w(117)
      w(129)=w(67)*z
      w(147)= - z*w(53)
      w(147)=w(147) + w(108) + 3 - w(129)
      w(147)=w(7)*w(147)
      w(84)=w(84)*z
      w(84)=w(108) + w(84) - 1
      w(149)= - w(77) - w(84)
      w(149)=w(6)*w(149)
      w(151)= - 109 - 268*z
      w(151)=z*w(151)
      w(151)= - 68*w(1) + 445 + w(151)
      w(113)=128*w(149) + 64*w(147) + w(113) + 2.D0/3.D0*w(151) + 
     & w(117)
      w(113)=LL*w(113)
      w(99)= - 1 - w(99)
      w(99)=z*w(99)
      w(57)= - w(131) - w(54) + w(99) + w(57)
      w(57)=LL*w(57)
      w(99)=w(75)*w(91)
      w(117)=12*z
      w(131)=11.D0/3.D0 - w(117)
      w(131)=z*w(131)
      w(99)=w(99) + w(137) + 1.D0/3.D0 + w(131)
      w(99)=w(3)*w(99)
      w(131)=28 - w(133)
      w(131)=z*w(131)
      w(131)= - 13.D0/3.D0*w(1) - 32.D0/3.D0 + w(131)
      w(57)=w(57) + 2*w(131) + w(99)
      w(57)=w(57)*w(102)
      w(99)=7 + w(148)
      w(99)=z*w(99)
      w(99)=2.D0/3.D0*w(100) + w(51) - 15 + w(99)
      w(99)=w(3)*w(99)
      w(88)= - 95 - w(88)
      w(88)=z*w(88)
      w(88)=145 + w(88)
      w(88)=1.D0/3.D0*w(88) - w(143)
      w(88)=2*w(88) + w(99)
      w(88)=w(3)*w(88)
      w(99)=80.D0/3.D0*z
      w(100)= - 87 + w(99)
      w(100)=z*w(100)
      w(57)=w(57) + w(88) + 3*w(1) + 172.D0/3.D0 + w(100)
      w(88)= - w(150)*w(134)
      w(88)=w(63) + 17 + w(88)
      w(88)=1.D0/3.D0*w(88) + w(62)
      w(88)=LL*w(88)
      w(100)= - 37.D0/3.D0 + w(138)
      w(100)=z*w(100)
      w(88)=w(88) - 17.D0/3.D0*w(1) + 10 + w(100)
      w(100)= - 7.D0/3.D0 + w(101)
      w(100)=z*w(100)
      w(100)=w(124) + w(105) - 5.D0/3.D0 + w(100)
      w(100)=4*w(100) + w(55)
      w(100)=w(100)*w(130)
      w(88)=2*w(88) + w(100)
      w(88)=w(9)*w(88)
      w(57)=2*w(57) + w(88)
      w(57)=w(9)*w(57)
      w(56)=w(56) - 1
      w(88)=9 + w(65)
      w(88)=z*w(88)
      w(88)= - w(83) + w(88) - w(56)
      w(88)=w(88)*w(102)
      w(100)=w(52)*w(91)
      w(131)= - 23.D0/3.D0 + w(65)
      w(131)=z*w(131)
      w(88)=8*w(119) + w(109) + w(88) + w(100) - 22.D0/3.D0 + w(131)
      w(88)=w(88)*w(90)
      w(90)=10*w(72)
      w(100)= - 7 + w(148)
      w(100)=z*w(100)
      w(100)= - 26*w(68) - w(90) + w(51) - 11 + w(100)
      w(100)=w(100)*w(140)
      w(109)= - 3 + z
      w(109)=w(109)*w(142)
      w(119)= - 17 + w(117)
      w(119)=z*w(119)
      w(56)=w(104) + w(109) + w(119) - w(56)
      w(56)=w(56)*w(71)
      w(109)= - 17 - 44.D0/3.D0*z
      w(109)=w(109)*w(106)
      w(109)=2 + w(109)
      w(109)=2*w(109) + 11.D0/3.D0*w(72)
      w(109)=w(109)*w(53)
      w(119)=1673.D0/9.D0 - w(138)
      w(119)=z*w(119)
      w(56)= - 44.D0/3.D0*w(55) + w(56) + w(109) + 1009.D0/9.D0 + 
     & w(119)
      w(56)=w(56)*w(114)
      w(109)= - w(2)*w(3)
      w(109)=w(109) - z2
      w(109)=w(71)*w(109)
      w(114)=w(4)*w(60)
      w(109)=w(114) + w(109)
      w(109)=w(84)*w(109)
      w(84)=w(84)*w(53)
      w(114)=67.D0/3.D0*w(67)
      w(84)= - w(114) + w(84)
      w(119)=LL*w(3)
      w(84)=w(84)*w(119)
      w(84)=w(84) + w(109)
      w(84)=w(2)*w(84)
      w(109)=z + 5
      w(131)= - w(86) + w(109)
      w(131)=w(16)*w(131)
      w(84)=w(131) + w(84)
      w(107)= - 181 - w(107)
      w(107)=w(107)*w(136)
      w(131)=5.D0/3.D0*w(1)
      w(140)=w(28)*w(67)
      w(107)=48*w(140) - w(131) + 141 + w(107)
      w(140)= - 23 - 64*z
      w(140)=w(140)*w(106)
      w(142)= - 1.D0/5.D0*w(72) + 1 + w(82)
      w(142)=w(3)*w(142)
      w(64)=149 + w(64)
      w(64)=z*w(64)
      w(64)=115 + w(64)
      w(64)=1.D0/3.D0*w(64) + w(142)
      w(64)=w(64)*w(103)
      w(64)=w(64) + 23 + w(140)
      w(64)=w(3)*w(64)
      w(99)=11 + w(99)
      w(99)=z*w(99)
      w(99)=82 + w(99)
      w(64)=2*w(99) + w(64)
      w(64)=w(3)*w(64)
      w(75)=2*w(75) + w(77)
      w(75)=w(75)*w(53)
      w(75)=w(114) + w(75)
      w(75)=LL*w(75)
      w(77)=w(71)*w(120)
      w(75)=w(75) + w(77)
      w(77)=16*w(4)
      w(75)=w(75)*w(77)
      w(50)= - w(104) + w(137) + 9 - 5*w(50)
      w(50)=w(18)*w(50)
      w(99)=160*w(33) - w(128) - 80*z5 - 208*w(34) - 96*w(35) + 384*
     & w(30)
      w(99)=w(67)*w(99)
      w(114)=w(14)*w(68)
      w(50)=w(73) + 192*w(114) + w(74) + w(88) + 8*w(50) + w(69) + 
     & w(56) + w(57) + w(75) + w(100) + 2*w(107) + w(64) + w(99) + 
     & w(113) + 16*w(84)
      w(56)=CF*TF
      w(50)=w(50)*w(56)
      w(57)=14*w(1)
      w(64)= - 223 - w(138)
      w(64)=z*w(64)
      w(64)= - w(57) - 271 + w(64)
      w(69)= - 37 - w(65)
      w(69)=z*w(69)
      w(69)= - 25 + w(69)
      w(69)=1.D0/3.D0*w(69) - w(72)
      w(69)=w(3)*w(69)
      w(73)=1 + w(101)
      w(73)=z*w(73)
      w(73)= - w(83) + 2 + w(73)
      w(73)=w(73)*w(71)
      w(64)= - w(55) + w(73) + 1.D0/9.D0*w(64) + w(69)
      w(69)=2*nf
      w(64)=w(64)*w(69)
      w(73)=z - 2
      w(74)=w(73)*w(101)
      w(74)=w(74) - 1
      w(74)=1.D0/3.D0*w(74)
      w(75)=w(80) + w(68)
      w(83)=w(74) - w(75)
      w(83)=w(83)*w(60)
      w(84)= - 61 + w(138)
      w(84)=z*w(84)
      w(84)= - 31 + w(84)
      w(84)=1.D0/3.D0*w(84) - w(125)
      w(84)=w(3)*w(84)
      w(88)=w(67)*nf
      w(99)=w(88) - w(67)
      w(99)=w(99)*w(111)
      w(100)=13 + 67.D0/9.D0*z
      w(100)=z*w(100)
      w(100)=41.D0/9.D0*w(1) - 25 + w(100)
      w(64)=w(99) + w(64) + w(92) + w(83) + 2*w(100) + w(84)
      w(67)=1.D0/3.D0*w(67) - 3.D0/5.D0*w(88)
      w(67)=w(67)*w(95)
      w(64)=1.D0/3.D0*w(64) + w(67)
      w(64)=w(64)*w(85)
      w(67)=1 + 38.D0/9.D0*z
      w(67)=w(67)*w(106)
      w(67)= - 1.D0/9.D0*w(55) - w(124) - 20.D0/27.D0*w(1) - 1 + w(67)
      w(67)=w(9)*w(67)
      w(83)= - 13 - w(136)
      w(83)=z*w(83)
      w(62)= - w(62) - 28.D0/3.D0*w(1) + 23 + w(83)
      w(62)=LL*w(62)
      w(83)= - 1 - 55.D0/9.D0*z
      w(83)=z*w(83)
      w(83)=28.D0/9.D0*w(1) + 4 + w(83)
      w(62)=4.D0/3.D0*w(83) + w(62)
      w(62)=2.D0/3.D0*w(62) + w(67)
      w(62)=w(62)*w(127)
      w(67)=22*z
      w(83)= - 151 - w(67)
      w(83)=w(83)*w(122)
      w(83)= - 518 + w(83)
      w(84)=5.D0/3.D0*z
      w(92)=w(84) + 1
      w(92)=1.D0/3.D0*w(72) + 2*w(92)
      w(99)= - w(92)*w(103)
      w(100)= - 2.D0/3.D0 - z
      w(99)=8*w(100) + w(99)
      w(99)=w(3)*w(99)
      w(83)=8.D0/81.D0*w(83) + w(99)
      w(83)=w(3)*w(83)
      w(52)=w(52) - w(80)
      w(99)=LL**2
      w(100)=w(52)*w(99)
      w(103)= - 4 - w(106)
      w(103)=w(103)*w(134)
      w(100)= - w(100) + 61.D0/3.D0*w(1) + 10 + w(103)
      w(103)=16*z
      w(106)= - 43 - w(103)
      w(106)=z*w(106)
      w(106)= - 37 + w(106)
      w(106)=1.D0/3.D0*w(106) + w(80)
      w(106)=w(3)*w(106)
      w(107)= - z*w(109)
      w(107)=41 + w(107)
      w(106)=4.D0/3.D0*w(107) + w(106)
      w(106)=w(3)*w(106)
      w(100)=w(106) + 4.D0/3.D0*w(100)
      w(106)=2.D0/3.D0*LL
      w(100)=w(100)*w(106)
      w(107)=74 + 569.D0/9.D0*z
      w(107)=w(107)*w(101)
      w(107)= - 193.D0/9.D0*w(1) - 253 + w(107)
      w(62)=w(62) + w(100) + 8.D0/27.D0*w(107) + w(83)
      w(62)=w(62)*w(69)
      w(83)=2 + 5.D0/9.D0*z
      w(83)=w(83)*w(101)
      w(83)=w(68) - w(79) - 4.D0/9.D0*w(1) + 7.D0/3.D0 + w(83)
      w(83)=nf*w(83)
      w(83)=w(68) + w(83)
      w(83)=5.D0/3.D0*w(52) + 4*w(83)
      w(83)=z3*w(83)
      w(100)=7 + 22.D0/9.D0*z
      w(100)=w(100)*z
      w(100)=w(100) - 13 + 32.D0/9.D0*w(1)
      w(54)= - w(124) - w(54) - w(100)
      w(54)=w(54)*w(102)
      w(55)= - LL*w(55)
      w(54)=w(54) + w(55)
      w(54)=w(9)*w(54)
      w(54)=w(54) + w(83)
      w(55)=w(101) + 7
      w(55)=w(55)*z
      w(63)= - w(63) - 2 + w(55)
      w(63)=1.D0/3.D0*w(63) + w(75)
      w(63)=LL*w(63)
      w(75)= - 4 + w(139) + w(86)
      w(75)=LL*w(75)
      w(83)=67 + 19*z
      w(83)=z*w(83)
      w(83)=43 + w(83)
      w(75)=1.D0/9.D0*w(83) + w(75)
      w(75)=nf*w(75)
      w(63)=w(63) + w(75)
      w(63)=w(5)*w(63)
      w(75)=w(18)*w(88)
      w(63)=w(63) + w(75)
      w(52)=w(52)*LL
      w(74)= - w(74) + w(72)
      w(74)=w(74)*w(53)
      w(74)=w(74) - w(100)
      w(74)=2*w(74) - 5.D0/3.D0*w(52)
      w(74)=w(74)*w(102)
      w(75)= - 13 - w(103)
      w(75)=z*w(75)
      w(75)=w(90) - 19 + w(75)
      w(75)=w(3)*w(75)
      w(83)=13 - w(141)
      w(83)=z*w(83)
      w(83)=76 + w(83)
      w(75)=8.D0/3.D0*w(83) + w(75)
      w(75)=w(3)*w(75)
      w(83)= - 62 - 73.D0/9.D0*z
      w(83)=z*w(83)
      w(83)=73.D0/9.D0*w(1) + 62 + w(83)
      w(75)=4*w(83) + w(75)
      w(74)=1.D0/3.D0*w(75) + w(74)
      w(74)=w(74)*w(106)
      w(75)=15 + w(89)
      w(75)=z*w(75)
      w(75)=w(131) - 23 + w(75)
      w(83)= - w(3)*w(92)
      w(76)=8*w(76)
      w(83)= - w(76) + w(83)
      w(83)=w(3)*w(83)
      w(86)= - 5 - w(134)
      w(83)=16*w(86) + w(83)
      w(83)=w(3)*w(83)
      w(86)= - 11 - w(101)
      w(86)=z*w(86)
      w(86)= - 8 + w(86)
      w(86)=1.D0/9.D0*w(86) + w(68)
      w(86)=nf*w(86)
      w(86)=w(115) + w(86)
      w(86)=w(10)*w(86)
      w(88)=w(104)*nf
      w(88)=w(88) + w(68)
      w(92)=w(8)*w(88)
      w(54)= - 32.D0/3.D0*w(92) + w(64) + 16*w(86) + w(62) + w(74) + 8*
     & w(75) + w(83) + 4.D0/3.D0*w(54) + 16.D0/3.D0*w(63)
      w(62)=TF**2
      w(62)=2*w(62)
      w(54)=w(54)*w(62)
      w(50)=w(50) + w(54) + w(58)
      w(50)=CF*w(50)
      w(54)= - 43 - w(138)
      w(54)=z*w(54)
      w(54)= - 34 + w(54)
      w(54)=z*w(54)
      w(54)=5 + w(54)
      w(54)=z*w(54)
      w(54)=8 + w(54)
      w(54)=w(54)*w(119)
      w(58)=w(116)*z
      w(58)=w(58) + 3
      w(58)=w(58)*z
      w(58)=w(58) + 2
      w(58)=w(58)*z
      w(58)=w(58) + 1
      w(63)=w(58)*w(66)
      w(54)=w(54) + w(63)
      w(54)=w(46)*w(54)
      w(63)=z2*w(94)
      w(63)= - 6.D0/5.D0*w(63) - 24*w(120)
      w(63)=LL*w(63)
      w(64)=w(122) + 3
      w(66)=w(64)*w(102)
      w(73)=w(73)*z
      w(74)=w(73) - 1
      w(74)=w(74)*z
      w(74)=w(74) + 2
      w(74)=w(74)*z
      w(74)=w(74) + 1
      w(74)=w(74)*w(41)
      w(66)=4*w(74) + w(66)
      w(75)=w(3)**2
      w(66)=w(75)*w(66)
      w(83)=1 + w(65)
      w(83)=w(83)*w(141)
      w(83)=47 + w(83)
      w(83)=w(83)*w(53)
      w(86)=437 - 416*z
      w(86)=z*w(86)
      w(86)=416*w(1) - 433 + w(86)
      w(83)=1.D0/3.D0*w(86) + w(83)
      w(86)= - 137 + 132*z
      w(86)=z*w(86)
      w(86)=156 + w(86)
      w(86)=z*w(86)
      w(86)=481 + w(86)
      w(86)=z*w(86)
      w(86)=96 + w(86)
      w(86)=w(47)*w(86)*w(102)
      w(55)=w(105) + 7 - w(55)
      w(55)=w(9)*w(55)*w(71)
      w(92)=w(5)*w(68)
      w(54)=2*w(54) + 24*w(92) + w(55) + w(86) + 1.D0/9.D0*w(83) + 
     & w(66) + w(63)
      w(54)=z2*w(54)
      w(55)=w(4)*w(71)
      w(63)= - w(2)*w(98)
      w(55)=w(63) + w(55)
      w(63)=w(93)*z
      w(63)=w(105) + w(63) - 7
      w(55)=w(63)*w(55)
      w(66)=11 - w(136)
      w(66)=z*w(66)
      w(66)= - w(70) + 11 + w(66)
      w(66)=LL*w(66)
      w(58)=w(58)*w(46)
      w(70)= - w(53)*w(58)
      w(66)=w(66) + w(70)
      w(66)=w(66)*w(85)
      w(48)= - 9 + w(48)
      w(48)=z*w(48)
      w(48)=w(51) - 9 + w(48)
      w(51)=w(75)*LL
      w(48)=w(48)*w(51)
      w(70)= - 77 + 114*z
      w(70)=w(70)*z
      w(70)=w(70) - 342
      w(70)=w(70)*z
      w(70)=w(70) - 77
      w(70)=w(70)*z
      w(70)=w(70) + 114
      w(70)=w(70)*w(47)
      w(83)=w(98)*w(70)
      w(85)= - w(5)*w(110)*w(60)
      w(58)=w(58)*w(99)
      w(86)=w(91)*w(58)
      w(48)=w(66) + w(86) + w(85) + w(83) + w(48) + w(55)
      w(48)=w(2)*w(48)
      w(55)= - 7 + w(132)
      w(55)=w(55)*w(133)
      w(55)= - 208.D0/3.D0*w(1) + 68 + w(55)
      w(66)=z**2
      w(66)= - 1 - 11.D0/9.D0*w(66)
      w(66)=w(66)*w(91)
      w(83)=w(61)*z
      w(83)= - w(1) + w(83) + 2
      w(85)=w(83)*LL
      w(55)=22.D0/27.D0*w(85) + 1.D0/9.D0*w(55) + w(66)
      w(55)=w(55)*w(102)
      w(66)=w(110)*w(96)
      w(86)= - 4 + w(89)
      w(89)= - w(3)*z
      w(86)=2*w(86) + w(89)
      w(86)=w(86)*w(3)**3
      w(89)= - 13763.D0/9.D0 + 1466*z
      w(89)=z*w(89)
      w(89)= - 1378*w(1) + 12661.D0/9.D0 + w(89)
      w(86)=1.D0/3.D0*w(89) + w(86)
      w(55)=w(66) + 1.D0/3.D0*w(86) + w(55)
      w(55)=LL*w(55)
      w(66)=w(61)*w(51)
      w(86)=w(102)*w(120)
      w(66)=w(66) + w(86)
      w(70)= - w(71)*w(70)
      w(58)= - 16*w(58) + w(70) + 12*w(66)
      w(58)=w(4)*w(58)
      w(66)=8 + z
      w(66)=z*w(66)
      w(66)= - 16 + w(66)
      w(66)=z*w(66)
      w(66)=6 + w(66)
      w(66)=w(66)*w(119)*w(77)
      w(70)= - 54 + 25*z
      w(70)=z*w(70)
      w(70)=25 + w(70)
      w(70)=w(9)*w(70)*w(51)
      w(77)=54 - w(123)
      w(77)=z*w(77)
      w(77)= - 33 + w(77)
      w(77)=w(5)*w(77)*w(98)
      w(66)=w(77) + w(66) + w(70)
      w(66)=w(40)*w(66)
      w(70)=w(73) + 3
      w(70)=w(70)*z
      w(70)=w(70) - 2
      w(70)=w(70)*z
      w(70)=w(70) + 1
      w(70)=w(70)*w(53)*w(9)
      w(73)= - w(99)*w(70)
      w(77)=1 - w(138)
      w(77)=z*w(77)
      w(77)= - 24 + w(77)
      w(77)=z*w(77)
      w(77)=41 + w(77)
      w(77)=z*w(77)
      w(77)= - 11 + w(77)
      w(77)=z3*LL*w(77)
      w(73)=w(73) + w(77)
      w(77)= - 11 - w(65)
      w(77)=w(77)*w(121)
      w(77)=128 + w(77)
      w(77)=z*w(77)
      w(77)= - 79 + w(77)
      w(77)=z*w(77)
      w(77)=4 + w(77)
      w(77)=w(7)*LL*w(77)
      w(70)=z2*w(70)
      w(70)=w(70) + 2*w(73) + w(77)
      w(70)=w(45)*w(70)
      w(73)=1 - w(122)
      w(73)=w(73)*w(145)
      w(77)=62 - 33*z
      w(77)=z*w(77)
      w(77)= - 78 + w(77)
      w(77)=z*w(77)
      w(77)=62 + w(77)
      w(77)=z*w(77)
      w(77)= - 33 + w(77)
      w(77)=w(44)*w(9)*w(77)
      w(86)=5413 - 3918*z
      w(86)=z*w(86)
      w(86)= - 5517 + w(86)
      w(86)=z*w(86)
      w(86)=3418 + w(86)
      w(86)=z*w(86)
      w(86)=568 + w(86)
      w(86)=w(43)*w(86)
      w(73)=w(73) + w(86) + 8*w(77)
      w(73)=w(119)*w(73)
      w(77)=131 + 82.D0/3.D0*z
      w(77)=w(77)*w(122)
      w(77)= - 410.D0/3.D0*w(1) - 677 + w(77)
      w(77)=LL*w(77)
      w(49)=w(49) - w(1)
      w(86)=w(9)*z
      w(77)= - 11*w(86) + 44*w(49) + w(77)
      w(89)=1.D0/9.D0*w(9)
      w(77)=w(77)*w(89)
      w(74)=w(75)*w(99)*w(74)
      w(91)=w(129) + 2
      w(91)=w(91)*w(101)
      w(91)=w(91) + 3
      w(92)=2*w(99) - z2
      w(91)=w(42)*w(95)*w(91)*w(92)
      w(92)=w(14) - w(15)
      w(92)=48*w(92)
      w(68)=w(68)*w(92)
      w(67)=w(67) + 13
      w(92)=w(67) + w(80)
      w(92)=w(92)*w(75)
      w(93)=9665 - 12283*z
      w(93)=z*w(93)
      w(93)=7927*w(1) - 8143 + w(93)
      w(93)=2.D0/9.D0*w(93) + 11*w(92)
      w(94)=w(119)*w(117)
      w(94)= - 11.D0/27.D0*w(83) + w(94)
      w(94)=w(94)*w(112)
      w(75)= - 229.D0/9.D0 - 3*w(75)
      w(75)=w(5)*w(75)*w(104)
      w(63)= - w(78) - w(63)
      w(60)=w(6)*w(63)*w(60)
      w(63)=583 - 124*z
      w(63)=z*w(63)
      w(63)=154 + w(63)
      w(63)=z*w(63)
      w(63)= - 583 + w(63)
      w(63)=z*w(63)
      w(63)= - 110 + w(63)
      w(51)=w(39)*w(63)*w(51)
      w(63)=w(72)*w(118)
      w(65)=27 + w(65)
      w(65)=z*w(65)
      w(65)= - 54 + w(65)
      w(65)=z*w(65)
      w(65)=19 + w(65)
      w(65)=w(40)*LL*w(65)
      w(63)=w(63) + w(65)
      w(63)=w(63)*w(97)
      w(65)= - 99 + 161*z
      w(65)=w(65)*z
      w(65)=w(65) - 71
      w(78)=w(3)*w(65)*w(37)
      w(48)=w(63) + w(91) + 4*w(70) + w(51) + 2*w(48) + w(54) + w(60)
     &  - 8*w(74) + w(66) + w(75) + w(94) + w(77) + 22*w(78) + 1.D0/27.D
     & 0*w(93) + w(73) + w(58) + w(55) + w(68)
      w(48)=w(48)*w(59)
      w(51)=w(81) - 19
      w(51)=w(51)*z
      w(51)= - w(51) + w(135) - 29
      w(51)=w(80) + 1.D0/3.D0*w(51)
      w(54)= - w(87) - w(51)
      w(54)=w(54)*w(69)
      w(54)=w(54) - 7*w(51) - w(87)
      w(54)=z2*w(54)
      w(55)=2*w(51) - 7.D0/3.D0*w(85)
      w(55)=w(55)*w(71)
      w(58)=52*z
      w(60)=89 - w(58)
      w(60)=z*w(60)
      w(60)=2 + w(60)
      w(60)=1.D0/3.D0*w(60) + w(90)
      w(60)=w(3)*w(60)
      w(63)=211 - 814.D0/3.D0*z
      w(63)=z*w(63)
      w(63)=598.D0/3.D0*w(1) - 259 + w(63)
      w(55)=w(55) + 1.D0/3.D0*w(63) + w(60)
      w(55)=LL*w(55)
      w(60)= - 949 + 1187*z
      w(60)=w(60)*z
      w(60)=w(60) + 881 - 791*w(1)
      w(60)= - w(92) + 2.D0/9.D0*w(60)
      w(54)=w(54) + w(55) + w(60)
      w(55)=85.D0/3.D0 - w(144)
      w(55)=z*w(55)
      w(55)=w(57) - 149.D0/3.D0 + w(55)
      w(57)=1 - w(58)
      w(57)=z*w(57)
      w(57)= - 50 + w(57)
      w(57)=1.D0/9.D0*w(57) + w(80)
      w(57)=w(3)*w(57)
      w(58)=w(83)*w(99)
      w(55)= - 4.D0/9.D0*w(58) + 1.D0/9.D0*w(55) + w(57)
      w(55)=LL*w(55)
      w(55)=1.D0/9.D0*w(60) + w(55)
      w(49)= - w(86) + 4*w(49)
      w(57)=w(126) - 11
      w(58)= - 13*w(82) + w(57)
      w(58)=LL*w(58)
      w(58)=w(58) - w(49)
      w(58)=w(58)*w(89)
      w(53)=w(65)*w(53)
      w(60)= - w(37)*w(53)
      w(55)=w(58) + 1.D0/3.D0*w(55) + w(60)
      w(55)=w(55)*w(69)
      w(58)=17 - 52.D0/3.D0*z
      w(58)=z*w(58)
      w(57)=w(58) + w(57)
      w(57)=w(57)*w(146)
      w(49)=w(57) - w(49)
      w(49)=w(49)*w(130)
      w(53)= - w(38)*w(53)
      w(57)=7 + w(69)
      w(57)=z3*w(83)*w(57)
      w(58)=w(5)*w(88)
      w(49)=16.D0/9.D0*w(58) + 4.D0/27.D0*w(57) + w(55) + w(49) + w(53)
     &  + 1.D0/9.D0*w(54)
      w(49)=w(49)*w(62)
      w(48)=w(49) + w(48)
      w(49)=64*w(23) + 192*w(26)
      w(49)=w(49)*w(56)*w(61)
      w(48)=2*w(48) + w(49)
      w(48)=CA*w(48)
      w(48)=w(50) + w(48)
      w(48)=as*w(48)
      w(49)= - w(85) + w(51)
      w(49)=w(49)*w(102)
      w(50)=1.D0/3.D0*w(67) + w(72)
      w(50)=w(3)*w(50)
      w(51)=137 - 175*z
      w(51)=z*w(51)
      w(51)=139*w(1) - 157 + w(51)
      w(49)=w(86) + w(49) + 1.D0/9.D0*w(51) + w(50)
      w(49)=w(49)*w(59)
      w(50)= - 2 - w(84)
      w(50)=z*w(50)
      w(50)= - w(108) + 4 + w(50)
      w(51)=w(72) + w(64)
      w(51)=w(3)*w(51)
      w(50)=2*w(50) + w(51)
      w(50)=2*w(50) - w(52)
      w(50)=w(50)*w(102)
      w(51)= - w(116)*w(121)
      w(51)= - w(1) + 10 + w(51)
      w(52)=w(79) + w(64)
      w(52)=w(3)*w(52)
      w(52)=w(76) + w(52)
      w(52)=w(3)*w(52)
      w(50)=w(50) + 4*w(51) + w(52)
      w(50)=w(50)*w(56)
      w(48)=w(48) + 2.D0/3.D0*w(49) + w(50)

      GGREG = 2*as**2*w(48)
*
      G =
     &  + as**2 * ( 278.D0/9.D0 + 230.D0/9.D0*z**(-1) + 244.D0/9.D0*
     &    z**(-1)*LL + 68.D0/9.D0*z**(-1)*LL**2 + 52.D0/3.D0*LL - 4.D0/
     &    3.D0*LL**2 - 14.D0/9.D0*z + 4*z*LL + 4.D0/3.D0*z*LL**2 - 494.D
     &    0/9.D0*z**2 - 436.D0/9.D0*z**2*LL - 68.D0/9.D0*z**2*LL**2 + 
     &    30*Hr1(0) + 24*Hr1(0)*LL + 16.D0/3.D0*Hr1(0)*LL**2 + 140.D0/3.
     &    D0*Hr1(0)*z + 104.D0/3.D0*Hr1(0)*z*LL + 16.D0/3.D0*Hr1(0)*z*
     &    LL**2 + 6*Hr1(0)**2 + 16.D0/3.D0*Hr1(0)**2*LL + 26.D0/3.D0*
     &    Hr1(0)**2*z + 16.D0/3.D0*Hr1(0)**2*z*LL + 8.D0/9.D0*Hr1(0)**3
     &     + 8.D0/9.D0*Hr1(0)**3*z + 2*Hr1(1)*z - 336/( - 27 + 27*z)*z
     &     - 120/( - 9 + 9*z)*z*LL - 12/( - 3 + 3*z)*z*LL**2 - 12/(3 - 
     &    3*z)*LL**2 - 120/(9 - 9*z)*LL - 336/(27 - 27*z) )
      G = G + as**3 * (  - 1772.D0/81.D0 + 671108.D0/243.D0*z**(-1) - 
     &    485680.D0/243.D0*z**(-1)*LL - 37312.D0/81.D0*z**(-1)*LL**2 - 
     &    7888.D0/81.D0*z**(-1)*LL**3 + 43454.D0/81.D0*z**(-1)*z2 + 
     &    1664.D0/9.D0*z**(-1)*z2*LL + 128.D0/3.D0*z**(-1)*z2*LL**2 - 
     &    2432.D0/15.D0*z**(-1)*z2**2 - 69304.D0/729.D0*z**(-1)*nf + 
     &    4408.D0/81.D0*z**(-1)*nf*LL + 400.D0/81.D0*z**(-1)*nf*LL**3
     &     - 1052.D0/81.D0*z**(-1)*nf*z2 - 9392.D0/81.D0*z**(-1)*z3 + 
     &    8624.D0/27.D0*z**(-1)*z3*LL - 400.D0/81.D0*z**(-1)*z3*nf + 
     &    10304.D0/9.D0*LL + 4988.D0/9.D0*LL**2 + 9328.D0/81.D0*LL**3
     &     - 1852.D0/3.D0*z2 + 7532.D0/27.D0*z2*LL + 448.D0/3.D0*z2*
     &    LL**2 + 320.D0/27.D0*z2*LL**3 + 2272.D0/15.D0*z2**2 - 1444.D0/
     &    9.D0*z2**2*LL - 11768.D0/81.D0*nf + 196.D0/9.D0*nf*LL + 16.D0/
     &    27.D0*nf*LL**3 - 3652.D0/81.D0*nf*z2 - 64.D0/9.D0*nf*z2*LL - 
     &    32.D0/5.D0*nf*z2**2 - 49744.D0/81.D0*z3 - 688*z3*LL - 896.D0/
     &    9.D0*z3*LL**2 - 4688.D0/27.D0*z3*z2 + 496.D0/27.D0*z3*nf + 64.
     &    D0/9.D0*z3*nf*LL )
      G = G + as**3 * ( 3040.D0/9.D0*z5 + 467492.D0/81.D0*z - 47552.D0/
     &    9.D0*z*LL - 32468.D0/27.D0*z*LL**2 - 9328.D0/81.D0*z*LL**3 + 
     &    18988.D0/27.D0*z*z2 + 5764.D0/27.D0*z*z2*LL + 800.D0/9.D0*z*
     &    z2*LL**2 + 320.D0/27.D0*z*z2*LL**3 - 15584.D0/45.D0*z*z2**2
     &     + 9668.D0/45.D0*z*z2**2*LL + 1880.D0/81.D0*z*nf - 484.D0/9.D0
     &    *z*nf*LL - 16.D0/27.D0*z*nf*LL**3 - 4252.D0/81.D0*z*nf*z2 - 
     &    128.D0/9.D0*z*nf*z2*LL - 32.D0/5.D0*z*nf*z2**2 - 16064.D0/81.D
     &    0*z*z3 + 1072*z*z3*LL + 2560.D0/9.D0*z*z3*LL**2 - 5552.D0/27.D
     &    0*z*z3*z2 + 80.D0/3.D0*z*z3*nf + 64.D0/9.D0*z*z3*nf*LL - 
     &    11360.D0/9.D0*z*z5 - 2068268.D0/243.D0*z**2 + 1491376.D0/243.D
     &    0*z**2*LL + 89824.D0/81.D0*z**2*LL**2 + 7888.D0/81.D0*z**2*
     &    LL**3 - 107678.D0/81.D0*z**2*z2 + 6352.D0/9.D0*z**2*z2*LL + 
     &    416.D0/3.D0*z**2*z2*LL**2 - 21392.D0/135.D0*z**2*z2**2 + 
     &    158296.D0/729.D0*z**2*nf - 1816.D0/81.D0*z**2*nf*LL - 400.D0/
     &    81.D0*z**2*nf*LL**3 + 700.D0/81.D0*z**2*nf*z2 + 128.D0/9.D0*
     &    z**2*nf*z2*LL )
      G = G + as**3 * (  - 63760.D0/81.D0*z**2*z3 - 256*z**2*z3*LL + 
     &    784.D0/81.D0*z**2*z3*nf + 640.D0/9.D0*Hr1(-1)*z**(-1)*z2 - 
     &    1168.D0/27.D0*Hr1(-1)*z**(-1)*z2*LL + 128*Hr1(-1)*z**(-1)*z3
     &     - 128.D0/3.D0*Hr1(-1)*z2 + 952.D0/9.D0*Hr1(-1)*z2*LL - 96*
     &    Hr1(-1)*z3 + 64.D0/3.D0*Hr1(-1)*z*z2 + 952.D0/9.D0*Hr1(-1)*z*
     &    z2*LL - 96*Hr1(-1)*z*z3 + 1216.D0/9.D0*Hr1(-1)*z**2*z2 - 1168.
     &    D0/27.D0*Hr1(-1)*z**2*z2*LL + 128*Hr1(-1)*z**2*z3 - 128.D0/3.D
     &    0*Hr1(-1)**2*z**(-1)*z2 + 32*Hr1(-1)**2*z2 + 32*Hr1(-1)**2*z*
     &    z2 - 128.D0/3.D0*Hr1(-1)**2*z**2*z2 + 64.D0/3.D0*Hr1(-1)**3*
     &    Hr1(0) - 256.D0/9.D0*Hr1(-1)**3*Hr1(0)*z**(-1) + 64.D0/3.D0*
     &    Hr1(-1)**3*Hr1(0)*z - 256.D0/9.D0*Hr1(-1)**3*Hr1(0)*z**2 - 
     &    128.D0/3.D0*Hr1(-1)**2*Hr1(0) + 640.D0/9.D0*Hr1(-1)**2*Hr1(0)
     &    *z**(-1) - 1456.D0/27.D0*Hr1(-1)**2*Hr1(0)*z**(-1)*LL + 664.D0
     &    /9.D0*Hr1(-1)**2*Hr1(0)*LL + 64.D0/3.D0*Hr1(-1)**2*Hr1(0)*z
     &     + 664.D0/9.D0*Hr1(-1)**2*Hr1(0)*z*LL + 1216.D0/9.D0*Hr1(-1)
     &    **2*Hr1(0)*z**2 )
      G = G + as**3 * (  - 1456.D0/27.D0*Hr1(-1)**2*Hr1(0)*z**2*LL - 16
     &    *Hr1(-1)**2*Hr1(0)**2 + 64.D0/3.D0*Hr1(-1)**2*Hr1(0)**2*
     &    z**(-1) - 16*Hr1(-1)**2*Hr1(0)**2*z + 64.D0/3.D0*Hr1(-1)**2*
     &    Hr1(0)**2*z**2 - 64*Hr1(-1)**2*Hr2(0,-1) + 256.D0/3.D0*Hr1(-1
     &    )**2*Hr2(0,-1)*z**(-1) - 64*Hr1(-1)**2*Hr2(0,-1)*z + 256.D0/3.
     &    D0*Hr1(-1)**2*Hr2(0,-1)*z**2 + 1600.D0/9.D0*Hr1(-1)*Hr1(0) - 
     &    3008.D0/27.D0*Hr1(-1)*Hr1(0)*z**(-1) + 6032.D0/9.D0*Hr1(-1)*
     &    Hr1(0)*z**(-1)*LL + 544.D0/3.D0*Hr1(-1)*Hr1(0)*z**(-1)*LL**2
     &     - 112*Hr1(-1)*Hr1(0)*z**(-1)*z2 - 15208.D0/27.D0*Hr1(-1)*
     &    Hr1(0)*LL + 32*Hr1(-1)*Hr1(0)*LL**2 - 1280.D0/9.D0*Hr1(-1)*
     &    Hr1(0)*z - 3112.D0/27.D0*Hr1(-1)*Hr1(0)*z*LL + 32*Hr1(-1)*
     &    Hr1(0)*z*LL**2 - 11648.D0/27.D0*Hr1(-1)*Hr1(0)*z**2 + 10064.D0
     &    /9.D0*Hr1(-1)*Hr1(0)*z**2*LL + 544.D0/3.D0*Hr1(-1)*Hr1(0)*
     &    z**2*LL**2 - 112*Hr1(-1)*Hr1(0)*z**2*z2 + 64.D0/3.D0*Hr1(-1)*
     &    Hr1(0)**2 - 320.D0/9.D0*Hr1(-1)*Hr1(0)**2*z**(-1) + 2960.D0/
     &    27.D0*Hr1(-1)*Hr1(0)**2*z**(-1)*LL )
      G = G + as**3 * (  - 980.D0/9.D0*Hr1(-1)*Hr1(0)**2*LL - 32.D0/3.D0
     &    *Hr1(-1)*Hr1(0)**2*z - 980.D0/9.D0*Hr1(-1)*Hr1(0)**2*z*LL - 
     &    608.D0/9.D0*Hr1(-1)*Hr1(0)**2*z**2 + 2960.D0/27.D0*Hr1(-1)*
     &    Hr1(0)**2*z**2*LL - 16.D0/3.D0*Hr1(-1)*Hr1(0)**3 + 64.D0/9.D0
     &    *Hr1(-1)*Hr1(0)**3*z**(-1) - 16.D0/3.D0*Hr1(-1)*Hr1(0)**3*z
     &     + 64.D0/9.D0*Hr1(-1)*Hr1(0)**3*z**2 + 256.D0/3.D0*Hr1(-1)*
     &    Hr2(0,-1) - 1280.D0/9.D0*Hr1(-1)*Hr2(0,-1)*z**(-1) + 2624.D0/
     &    27.D0*Hr1(-1)*Hr2(0,-1)*z**(-1)*LL - 1616.D0/9.D0*Hr1(-1)*
     &    Hr2(0,-1)*LL - 128.D0/3.D0*Hr1(-1)*Hr2(0,-1)*z - 1616.D0/9.D0
     &    *Hr1(-1)*Hr2(0,-1)*z*LL - 2432.D0/9.D0*Hr1(-1)*Hr2(0,-1)*z**2
     &     + 2624.D0/27.D0*Hr1(-1)*Hr2(0,-1)*z**2*LL + 128*Hr1(-1)*Hr3(
     &    0,-1,-1) - 512.D0/3.D0*Hr1(-1)*Hr3(0,-1,-1)*z**(-1) + 128*
     &    Hr1(-1)*Hr3(0,-1,-1)*z - 512.D0/3.D0*Hr1(-1)*Hr3(0,-1,-1)*
     &    z**2 - 64*Hr1(-1)*Hr3(0,0,-1) + 256.D0/3.D0*Hr1(-1)*Hr3(0,0,
     &    -1)*z**(-1) - 64*Hr1(-1)*Hr3(0,0,-1)*z + 256.D0/3.D0*Hr1(-1)*
     &    Hr3(0,0,-1)*z**2 )
      G = G + as**3 * ( 128*Hr1(-1)*Hr3(0,0,1) - 512.D0/3.D0*Hr1(-1)*
     &    Hr3(0,0,1)*z**(-1) + 128*Hr1(-1)*Hr3(0,0,1)*z - 512.D0/3.D0*
     &    Hr1(-1)*Hr3(0,0,1)*z**2 + 85672.D0/27.D0*Hr1(0) + 20992.D0/81.
     &    D0*Hr1(0)*z**(-1) - 2768.D0/9.D0*Hr1(0)*z**(-1)*LL - 640.D0/9.
     &    D0*Hr1(0)*z**(-1)*LL**2 - 128.D0/9.D0*Hr1(0)*z**(-1)*LL**3 + 
     &    640.D0/9.D0*Hr1(0)*z**(-1)*z2 + 32.D0/3.D0*Hr1(0)*z**(-1)*z2*
     &    LL + 128.D0/9.D0*Hr1(0)*z**(-1)*z3 - 291464.D0/81.D0*Hr1(0)*
     &    LL - 15280.D0/27.D0*Hr1(0)*LL**2 - 896.D0/27.D0*Hr1(0)*LL**3
     &     + 680*Hr1(0)*z2 - 1988.D0/9.D0*Hr1(0)*z2*LL - 832.D0/9.D0*
     &    Hr1(0)*z2*LL**2 + 5024.D0/45.D0*Hr1(0)*z2**2 - 38264.D0/243.D0
     &    *Hr1(0)*nf + 2024.D0/27.D0*Hr1(0)*nf*LL + 128.D0/27.D0*Hr1(0)
     &    *nf*LL**3 - 472.D0/27.D0*Hr1(0)*nf*z2 - 256.D0/9.D0*Hr1(0)*nf
     &    *z2*LL + 8576.D0/27.D0*Hr1(0)*z3 - 1088.D0/3.D0*Hr1(0)*z3*LL
     &     - 128.D0/27.D0*Hr1(0)*z3*nf + 30992.D0/9.D0*Hr1(0)*z - 97844.
     &    D0/81.D0*Hr1(0)*z*LL - 496.D0/3.D0*Hr1(0)*z*LL**2 - 736.D0/27.
     &    D0*Hr1(0)*z*LL**3 )
      G = G + as**3 * ( 10856.D0/27.D0*Hr1(0)*z*z2 + 188.D0/9.D0*Hr1(0)
     &    *z*z2*LL + 320.D0/9.D0*Hr1(0)*z*z2*LL**2 - 1472.D0/9.D0*Hr1(0
     &    )*z*z2**2 - 59912.D0/243.D0*Hr1(0)*z*nf - 308.D0/27.D0*Hr1(0)
     &    *z*nf*LL + 128.D0/27.D0*Hr1(0)*z*nf*LL**3 - 664.D0/27.D0*Hr1(
     &    0)*z*nf*z2 - 256.D0/9.D0*Hr1(0)*z*nf*z2*LL - 16256.D0/27.D0*
     &    Hr1(0)*z*z3 + 7136.D0/9.D0*Hr1(0)*z*z3*LL - 128.D0/27.D0*Hr1(
     &    0)*z*z3*nf + 300512.D0/81.D0*Hr1(0)*z**2 - 342620.D0/81.D0*
     &    Hr1(0)*z**2*LL - 20512.D0/27.D0*Hr1(0)*z**2*LL**2 - 512.D0/81.
     &    D0*Hr1(0)*z**2*LL**3 + 1472.D0/3.D0*Hr1(0)*z**2*z2 + 32*Hr1(0
     &    )*z**2*z2*LL - 7040.D0/243.D0*Hr1(0)*z**2*nf - 688.D0/27.D0*
     &    Hr1(0)*z**2*nf*LL - 64.D0/27.D0*Hr1(0)*z**2*nf*z2 - 5632.D0/
     &    81.D0*Hr1(0)*z**2*z3 + 2308.D0/9.D0*Hr1(0)**2 - 15464.D0/27.D0
     &    *Hr1(0)**2*LL + 320.D0/9.D0*Hr1(0)**2*LL**2 + 352.D0/27.D0*
     &    Hr1(0)**2*LL**3 - 32.D0/3.D0*Hr1(0)**2*z2 - 1012.D0/9.D0*Hr1(
     &    0)**2*z2*LL - 20*Hr1(0)**2*nf - 376.D0/27.D0*Hr1(0)**2*nf*LL
     &     - 16.D0/9.D0*Hr1(0)**2*nf*z2 )
      G = G + as**3 * (  - 3040.D0/27.D0*Hr1(0)**2*z3 + 24352.D0/27.D0*
     &    Hr1(0)**2*z + 202*Hr1(0)**2*z*LL - 704.D0/9.D0*Hr1(0)**2*z*
     &    LL**2 - 512.D0/27.D0*Hr1(0)**2*z*LL**3 + 328.D0/3.D0*Hr1(0)**
     &    2*z*z2 - 364.D0/9.D0*Hr1(0)**2*z*z2*LL - 280.D0/9.D0*Hr1(0)**
     &    2*z*nf - 472.D0/27.D0*Hr1(0)**2*z*nf*LL - 16.D0/9.D0*Hr1(0)**
     &    2*z*nf*z2 + 4736.D0/27.D0*Hr1(0)**2*z*z3 - 8224.D0/27.D0*Hr1(
     &    0)**2*z**2 - 13592.D0/27.D0*Hr1(0)**2*z**2*LL - 1552.D0/27.D0
     &    *Hr1(0)**2*z**2*LL**2 + 904.D0/27.D0*Hr1(0)**2*z**2*z2 - 256.D
     &    0/27.D0*Hr1(0)**2*z**2*nf*LL + 6880.D0/81.D0*Hr1(0)**3 - 352.D
     &    0/3.D0*Hr1(0)**3*LL - 448.D0/27.D0*Hr1(0)**3*LL**2 + 112.D0/9.
     &    D0*Hr1(0)**3*z2 - 8.D0/3.D0*Hr1(0)**3*nf + 32.D0/9.D0*Hr1(0)
     &    **3*nf*LL + 6488.D0/81.D0*Hr1(0)**3*z + 2108.D0/9.D0*Hr1(0)**
     &    3*z*LL + 704.D0/27.D0*Hr1(0)**3*z*LL**2 - 176.D0/9.D0*Hr1(0)
     &    **3*z*z2 - 104.D0/27.D0*Hr1(0)**3*z*nf + 32.D0/9.D0*Hr1(0)**3
     &    *z*nf*LL + 1120.D0/81.D0*Hr1(0)**3*z**2 + 64.D0/27.D0*Hr1(0)
     &    **4 )
      G = G + as**3 * (  - 152.D0/9.D0*Hr1(0)**4*LL - 8.D0/27.D0*Hr1(0)
     &    **4*nf + 284.D0/27.D0*Hr1(0)**4*z + 226.D0/9.D0*Hr1(0)**4*z*
     &    LL - 8.D0/27.D0*Hr1(0)**4*z*nf + 64.D0/81.D0*Hr1(0)**4*z**2
     &     + 56.D0/135.D0*Hr1(0)**5 - 88.D0/135.D0*Hr1(0)**5*z - 64.D0/
     &    27.D0*Hr1(0)**3*Hr1(1) - 256.D0/81.D0*Hr1(0)**3*Hr1(1)*
     &    z**(-1) + 64.D0/27.D0*Hr1(0)**3*Hr1(1)*z + 256.D0/81.D0*Hr1(0
     &    )**3*Hr1(1)*z**2 - 416.D0/27.D0*Hr1(0)**3*Hr2(0,-1) + 160.D0/
     &    27.D0*Hr1(0)**3*Hr2(0,-1)*z - 1040.D0/3.D0*Hr1(0)**2*Hr1(1)
     &     + 3440.D0/27.D0*Hr1(0)**2*Hr1(1)*z**(-1) + 1792.D0/27.D0*
     &    Hr1(0)**2*Hr1(1)*z**(-1)*LL + 10.D0/9.D0*Hr1(0)**2*Hr1(1)*LL
     &     + 3008.D0/9.D0*Hr1(0)**2*Hr1(1)*z - 10.D0/9.D0*Hr1(0)**2*
     &    Hr1(1)*z*LL - 3104.D0/27.D0*Hr1(0)**2*Hr1(1)*z**2 - 1792.D0/
     &    27.D0*Hr1(0)**2*Hr1(1)*z**2*LL + 8*Hr1(0)**2*Hr1(1)**2 + 32.D0
     &    /3.D0*Hr1(0)**2*Hr1(1)**2*z**(-1) - 8*Hr1(0)**2*Hr1(1)**2*z
     &     - 32.D0/3.D0*Hr1(0)**2*Hr1(1)**2*z**2 - 64.D0/9.D0*Hr1(0)**2
     &    *Hr2(0,-1) )
      G = G + as**3 * ( 832.D0/27.D0*Hr1(0)**2*Hr2(0,-1)*z**(-1) - 1060.
     &    D0/9.D0*Hr1(0)**2*Hr2(0,-1)*LL - 64*Hr1(0)**2*Hr2(0,-1)*z + 
     &    236*Hr1(0)**2*Hr2(0,-1)*z*LL - 832.D0/27.D0*Hr1(0)**2*Hr2(0,
     &    -1)*z**2 - 64*Hr1(0)**2*Hr3(0,-1,-1) + 64*Hr1(0)**2*Hr3(0,-1,
     &    -1)*z + 64*Hr1(0)**2*Hr3(0,0,-1) - 64*Hr1(0)**2*Hr3(0,0,-1)*z
     &     - 64.D0/3.D0*Hr1(0)**2*Hr3(0,0,1) + 512.D0/3.D0*Hr1(0)**2*
     &    Hr3(0,0,1)*z + 32*Hr1(0)**2*Hr3(0,1,1) + 32*Hr1(0)**2*Hr3(0,1
     &    ,1)*z + 28480.D0/27.D0*Hr1(0)*Hr1(1) - 15680.D0/27.D0*Hr1(0)*
     &    Hr1(1)*z**(-1) + 21040.D0/27.D0*Hr1(0)*Hr1(1)*z**(-1)*LL + 
     &    4256.D0/27.D0*Hr1(0)*Hr1(1)*z**(-1)*LL**2 - 2384.D0/27.D0*
     &    Hr1(0)*Hr1(1)*z**(-1)*z2 + 512.D0/27.D0*Hr1(0)*Hr1(1)*z**(-1)
     &    *nf*LL - 18896.D0/27.D0*Hr1(0)*Hr1(1)*LL - 448.D0/9.D0*Hr1(0)
     &    *Hr1(1)*LL**2 + 160.D0/9.D0*Hr1(0)*Hr1(1)*z2 + 128.D0/9.D0*
     &    Hr1(0)*Hr1(1)*nf*LL - 40256.D0/27.D0*Hr1(0)*Hr1(1)*z + 25712.D
     &    0/27.D0*Hr1(0)*Hr1(1)*z*LL + 448.D0/9.D0*Hr1(0)*Hr1(1)*z*
     &    LL**2 )
      G = G + as**3 * (  - 160.D0/9.D0*Hr1(0)*Hr1(1)*z*z2 - 128.D0/9.D0
     &    *Hr1(0)*Hr1(1)*z*nf*LL + 9152.D0/9.D0*Hr1(0)*Hr1(1)*z**2 - 
     &    27856.D0/27.D0*Hr1(0)*Hr1(1)*z**2*LL - 4256.D0/27.D0*Hr1(0)*
     &    Hr1(1)*z**2*LL**2 + 2384.D0/27.D0*Hr1(0)*Hr1(1)*z**2*z2 - 512.
     &    D0/27.D0*Hr1(0)*Hr1(1)*z**2*nf*LL - 64.D0/3.D0*Hr1(0)*Hr1(1)
     &    **2 - 320.D0/9.D0*Hr1(0)*Hr1(1)**2*z**(-1) + 1216.D0/27.D0*
     &    Hr1(0)*Hr1(1)**2*z**(-1)*LL + 304.D0/9.D0*Hr1(0)*Hr1(1)**2*LL
     &     - 32.D0/3.D0*Hr1(0)*Hr1(1)**2*z - 304.D0/9.D0*Hr1(0)*Hr1(1)
     &    **2*z*LL + 608.D0/9.D0*Hr1(0)*Hr1(1)**2*z**2 - 1216.D0/27.D0*
     &    Hr1(0)*Hr1(1)**2*z**2*LL - 32.D0/3.D0*Hr1(0)*Hr1(1)**3 - 128.D
     &    0/9.D0*Hr1(0)*Hr1(1)**3*z**(-1) + 32.D0/3.D0*Hr1(0)*Hr1(1)**3
     &    *z + 128.D0/9.D0*Hr1(0)*Hr1(1)**3*z**2 - 608.D0/9.D0*Hr1(0)*
     &    Hr1(1)*Hr2(0,-1) - 2432.D0/27.D0*Hr1(0)*Hr1(1)*Hr2(0,-1)*
     &    z**(-1) + 608.D0/9.D0*Hr1(0)*Hr1(1)*Hr2(0,-1)*z + 2432.D0/27.D
     &    0*Hr1(0)*Hr1(1)*Hr2(0,-1)*z**2 + 9728.D0/27.D0*Hr1(0)*Hr2(0,
     &    -1) )
      G = G + as**3 * (  - 8800.D0/27.D0*Hr1(0)*Hr2(0,-1)*z**(-1) - 512.
     &    D0/3.D0*Hr1(0)*Hr2(0,-1)*z**(-1)*LL - 572.D0/3.D0*Hr1(0)*Hr2(
     &    0,-1)*LL - 320.D0/9.D0*Hr1(0)*Hr2(0,-1)*LL**2 + 320.D0/9.D0*
     &    Hr1(0)*Hr2(0,-1)*z2 + 256.D0/9.D0*Hr1(0)*Hr2(0,-1)*nf*LL - 
     &    24032.D0/27.D0*Hr1(0)*Hr2(0,-1)*z + 3196.D0/3.D0*Hr1(0)*Hr2(0
     &    ,-1)*z*LL + 1984.D0/9.D0*Hr1(0)*Hr2(0,-1)*z*LL**2 - 1408.D0/9.
     &    D0*Hr1(0)*Hr2(0,-1)*z*z2 + 256.D0/9.D0*Hr1(0)*Hr2(0,-1)*z*nf*
     &    LL + 1408.D0/27.D0*Hr1(0)*Hr2(0,-1)*z**2 - 2144.D0/27.D0*Hr1(
     &    0)*Hr2(0,-1)*z**2*LL - 32.D0/9.D0*Hr1(0)*Hr2(0,-1)**2 - 1184.D
     &    0/9.D0*Hr1(0)*Hr2(0,-1)**2*z - 64*Hr1(0)*Hr3(0,-1,-1) + 256.D0
     &    /3.D0*Hr1(0)*Hr3(0,-1,-1)*z**(-1) + 1760.D0/9.D0*Hr1(0)*Hr3(0
     &    ,-1,-1)*LL - 64*Hr1(0)*Hr3(0,-1,-1)*z - 1760.D0/9.D0*Hr1(0)*
     &    Hr3(0,-1,-1)*z*LL + 256.D0/3.D0*Hr1(0)*Hr3(0,-1,-1)*z**2 + 64
     &    *Hr1(0)*Hr3(0,-1,1) - 256.D0/3.D0*Hr1(0)*Hr3(0,-1,1)*z**(-1)
     &     + 64*Hr1(0)*Hr3(0,-1,1)*z - 256.D0/3.D0*Hr1(0)*Hr3(0,-1,1)*
     &    z**2 )
      G = G + as**3 * (  - 32*Hr1(0)*Hr3(0,0,-1) + 128.D0/3.D0*Hr1(0)*
     &    Hr3(0,0,-1)*z**(-1) + 272*Hr1(0)*Hr3(0,0,-1)*LL - 32*Hr1(0)*
     &    Hr3(0,0,-1)*z - 5072.D0/9.D0*Hr1(0)*Hr3(0,0,-1)*z*LL + 128.D0/
     &    3.D0*Hr1(0)*Hr3(0,0,-1)*z**2 + 2368.D0/9.D0*Hr1(0)*Hr3(0,0,1)
     &     - 5120.D0/27.D0*Hr1(0)*Hr3(0,0,1)*z**(-1) - 80*Hr1(0)*Hr3(0,
     &    0,1)*LL - 320.D0/9.D0*Hr1(0)*Hr3(0,0,1)*z + 304*Hr1(0)*Hr3(0,
     &    0,1)*z*LL - 512.D0/27.D0*Hr1(0)*Hr3(0,0,1)*z**2 + 64*Hr1(0)*
     &    Hr3(0,1,-1) - 256.D0/3.D0*Hr1(0)*Hr3(0,1,-1)*z**(-1) + 64*
     &    Hr1(0)*Hr3(0,1,-1)*z - 256.D0/3.D0*Hr1(0)*Hr3(0,1,-1)*z**2 + 
     &    928.D0/9.D0*Hr1(0)*Hr3(0,1,1) + 3712.D0/27.D0*Hr1(0)*Hr3(0,1,
     &    1)*z**(-1) + 1216.D0/9.D0*Hr1(0)*Hr3(0,1,1)*LL - 928.D0/9.D0*
     &    Hr1(0)*Hr3(0,1,1)*z + 1216.D0/9.D0*Hr1(0)*Hr3(0,1,1)*z*LL - 
     &    3712.D0/27.D0*Hr1(0)*Hr3(0,1,1)*z**2 + 256*Hr1(0)*Hr4(0,-1,-1
     &    ,-1) - 256*Hr1(0)*Hr4(0,-1,-1,-1)*z - 128*Hr1(0)*Hr4(0,-1,0,1
     &    ) + 128*Hr1(0)*Hr4(0,-1,0,1)*z - 192*Hr1(0)*Hr4(0,0,0,-1) + 
     &    192*Hr1(0)*Hr4(0,0,0,-1)*z )
      G = G + as**3 * ( 1024.D0/9.D0*Hr1(0)*Hr4(0,0,0,1) - 7040.D0/9.D0
     &    *Hr1(0)*Hr4(0,0,0,1)*z + 1280.D0/9.D0*Hr1(0)*Hr4(0,0,1,1) + 
     &    1280.D0/9.D0*Hr1(0)*Hr4(0,0,1,1)*z - 128*Hr1(0)*Hr4(0,1,1,1)
     &     - 128*Hr1(0)*Hr4(0,1,1,1)*z + 10288.D0/27.D0*Hr1(1) + 8552.D0
     &    /81.D0*Hr1(1)*z**(-1) - 2836.D0/81.D0*Hr1(1)*z**(-1)*LL - 896.
     &    D0/27.D0*Hr1(1)*z**(-1)*LL**2 - 640.D0/81.D0*Hr1(1)*z**(-1)*
     &    LL**3 + 304.D0/9.D0*Hr1(1)*z**(-1)*z2 + 304.D0/27.D0*Hr1(1)*
     &    z**(-1)*z2*LL + 4880.D0/243.D0*Hr1(1)*z**(-1)*nf - 272.D0/27.D
     &    0*Hr1(1)*z**(-1)*nf*LL + 64.D0/27.D0*Hr1(1)*z**(-1)*nf*z2 - 
     &    128.D0/81.D0*Hr1(1)*z**(-1)*z3 - 15508.D0/27.D0*Hr1(1)*LL - 
     &    784.D0/9.D0*Hr1(1)*LL**2 - 160.D0/27.D0*Hr1(1)*LL**3 + 1832.D0
     &    /27.D0*Hr1(1)*z2 + 376.D0/9.D0*Hr1(1)*z2*LL + 368.D0/27.D0*
     &    Hr1(1)*nf + 604.D0/9.D0*Hr1(1)*nf*LL + 16.D0/9.D0*Hr1(1)*nf*
     &    z2 - 32.D0/27.D0*Hr1(1)*z3 - 5288.D0/27.D0*Hr1(1)*z + 884.D0/
     &    9.D0*Hr1(1)*z*LL + 16*Hr1(1)*z*LL**2 + 160.D0/27.D0*Hr1(1)*z*
     &    LL**3 )
      G = G + as**3 * (  - 392.D0/27.D0*Hr1(1)*z*z2 - 376.D0/9.D0*Hr1(1
     &    )*z*z2*LL - 320.D0/27.D0*Hr1(1)*z*nf - 260.D0/9.D0*Hr1(1)*z*
     &    nf*LL - 16.D0/9.D0*Hr1(1)*z*nf*z2 + 32.D0/27.D0*Hr1(1)*z*z3
     &     - 12320.D0/81.D0*Hr1(1)*z**2 + 38164.D0/81.D0*Hr1(1)*z**2*LL
     &     + 2816.D0/27.D0*Hr1(1)*z**2*LL**2 + 640.D0/81.D0*Hr1(1)*z**2
     &    *LL**3 - 784.D0/9.D0*Hr1(1)*z**2*z2 - 304.D0/27.D0*Hr1(1)*
     &    z**2*z2*LL - 7040.D0/243.D0*Hr1(1)*z**2*nf - 688.D0/27.D0*
     &    Hr1(1)*z**2*nf*LL - 64.D0/27.D0*Hr1(1)*z**2*nf*z2 + 128.D0/81.
     &    D0*Hr1(1)*z**2*z3 + 976.D0/9.D0*Hr1(1)**2 - 8.D0/27.D0*Hr1(1)
     &    **2*z**(-1) + 16.D0/9.D0*Hr1(1)**2*z**(-1)*LL - 352.D0/27.D0*
     &    Hr1(1)**2*z**(-1)*z2 - 320.D0/81.D0*Hr1(1)**2*z**(-1)*nf + 64.
     &    D0/9.D0*Hr1(1)**2*z**(-1)*nf*LL - 776.D0/27.D0*Hr1(1)**2*LL
     &     - 88.D0/9.D0*Hr1(1)**2*z2 - 16.D0/3.D0*Hr1(1)**2*nf + 16.D0/
     &    3.D0*Hr1(1)**2*nf*LL - 1604.D0/27.D0*Hr1(1)**2*z + 248.D0/27.D
     &    0*Hr1(1)**2*z*LL + 88.D0/9.D0*Hr1(1)**2*z*z2 + 28.D0/9.D0*
     &    Hr1(1)**2*z*nf )
      G = G + as**3 * (  - 16.D0/3.D0*Hr1(1)**2*z*nf*LL - 1856.D0/27.D0
     &    *Hr1(1)**2*z**2 + 160.D0/9.D0*Hr1(1)**2*z**2*LL + 352.D0/27.D0
     &    *Hr1(1)**2*z**2*z2 + 608.D0/81.D0*Hr1(1)**2*z**2*nf - 64.D0/9.
     &    D0*Hr1(1)**2*z**2*nf*LL + 400.D0/81.D0*Hr1(1)**3 + 16.D0/3.D0
     &    *Hr1(1)**3*z**(-1) + 320.D0/81.D0*Hr1(1)**3*z**(-1)*LL + 64.D0
     &    /81.D0*Hr1(1)**3*z**(-1)*nf + 80.D0/27.D0*Hr1(1)**3*LL + 16.D0
     &    /27.D0*Hr1(1)**3*nf + 128.D0/81.D0*Hr1(1)**3*z - 80.D0/27.D0*
     &    Hr1(1)**3*z*LL - 16.D0/27.D0*Hr1(1)**3*z*nf - 320.D0/27.D0*
     &    Hr1(1)**3*z**2 - 320.D0/81.D0*Hr1(1)**3*z**2*LL - 64.D0/81.D0
     &    *Hr1(1)**3*z**2*nf + 20.D0/27.D0*Hr1(1)**4 + 80.D0/81.D0*Hr1(
     &    1)**4*z**(-1) - 20.D0/27.D0*Hr1(1)**4*z - 80.D0/81.D0*Hr1(1)
     &    **4*z**2 + 128.D0/9.D0*Hr1(1)**2*Hr2(0,-1) + 512.D0/27.D0*
     &    Hr1(1)**2*Hr2(0,-1)*z**(-1) - 128.D0/9.D0*Hr1(1)**2*Hr2(0,-1)
     &    *z - 512.D0/27.D0*Hr1(1)**2*Hr2(0,-1)*z**2 - 256.D0/27.D0*
     &    Hr1(1)*Hr2(0,-1) - 1280.D0/27.D0*Hr1(1)*Hr2(0,-1)*z**(-1)*LL
     &     - 320.D0/9.D0*Hr1(1)*Hr2(0,-1)*LL )
      G = G + as**3 * ( 256.D0/27.D0*Hr1(1)*Hr2(0,-1)*z + 320.D0/9.D0*
     &    Hr1(1)*Hr2(0,-1)*z*LL + 1280.D0/27.D0*Hr1(1)*Hr2(0,-1)*z**2*
     &    LL + 928.D0/9.D0*Hr1(1)*Hr3(0,0,1) + 3712.D0/27.D0*Hr1(1)*
     &    Hr3(0,0,1)*z**(-1) - 928.D0/9.D0*Hr1(1)*Hr3(0,0,1)*z - 3712.D0
     &    /27.D0*Hr1(1)*Hr3(0,0,1)*z**2 - 224.D0/9.D0*Hr1(1)*Hr3(0,1,1)
     &     - 896.D0/27.D0*Hr1(1)*Hr3(0,1,1)*z**(-1) + 224.D0/9.D0*Hr1(1
     &    )*Hr3(0,1,1)*z + 896.D0/27.D0*Hr1(1)*Hr3(0,1,1)*z**2 - 28096.D
     &    0/27.D0*Hr2(0,-1) + 18688.D0/27.D0*Hr2(0,-1)*z**(-1) - 26032.D
     &    0/27.D0*Hr2(0,-1)*z**(-1)*LL - 5408.D0/27.D0*Hr2(0,-1)*
     &    z**(-1)*LL**2 + 3536.D0/27.D0*Hr2(0,-1)*z**(-1)*z2 - 512.D0/
     &    27.D0*Hr2(0,-1)*z**(-1)*nf*LL + 26572.D0/27.D0*Hr2(0,-1)*LL
     &     - 1184.D0/9.D0*Hr2(0,-1)*LL**2 - 320.D0/27.D0*Hr2(0,-1)*
     &    LL**3 + 832.D0/9.D0*Hr2(0,-1)*z2 + 1184.D0/9.D0*Hr2(0,-1)*z2*
     &    LL + 2752.D0/81.D0*Hr2(0,-1)*nf - 64.D0/9.D0*Hr2(0,-1)*nf*LL
     &     + 32.D0/9.D0*Hr2(0,-1)*nf*z2 - 5248.D0/27.D0*Hr2(0,-1)*z3 + 
     &    47824.D0/27.D0*Hr2(0,-1)*z )
      G = G + as**3 * (  - 10492.D0/9.D0*Hr2(0,-1)*z*LL - 416.D0/3.D0*
     &    Hr2(0,-1)*z*LL**2 - 320.D0/27.D0*Hr2(0,-1)*z*LL**3 + 464.D0/3.
     &    D0*Hr2(0,-1)*z*z2 - 64*Hr2(0,-1)*z*z2*LL + 4288.D0/81.D0*Hr2(
     &    0,-1)*z*nf + 256.D0/9.D0*Hr2(0,-1)*z*nf*LL + 32.D0/9.D0*Hr2(0
     &    ,-1)*z*nf*z2 + 5120.D0/27.D0*Hr2(0,-1)*z*z3 - 4160.D0/9.D0*
     &    Hr2(0,-1)*z**2 - 21392.D0/27.D0*Hr2(0,-1)*z**2*LL - 4384.D0/
     &    27.D0*Hr2(0,-1)*z**2*LL**2 + 1264.D0/9.D0*Hr2(0,-1)*z**2*z2
     &     + 1216.D0/81.D0*Hr2(0,-1)*z**2*nf + 128.D0/27.D0*Hr2(0,-1)*
     &    z**2*nf*LL - 32*Hr2(0,-1)**2 - 640.D0/27.D0*Hr2(0,-1)**2*
     &    z**(-1) - 400.D0/3.D0*Hr2(0,-1)**2*LL - 224.D0/9.D0*Hr2(0,-1)
     &    **2*z + 560.D0/9.D0*Hr2(0,-1)**2*z*LL + 128.D0/27.D0*Hr2(0,-1
     &    )**2*z**2 - 256*Hr2(0,-1)*Hr3(0,-1,-1) + 256*Hr2(0,-1)*Hr3(0,
     &    -1,-1)*z - 128*Hr2(0,-1)*Hr3(0,0,-1) + 128*Hr2(0,-1)*Hr3(0,0,
     &    -1)*z + 4160.D0/9.D0*Hr2(0,-1)*Hr3(0,0,1) - 448.D0/9.D0*Hr2(0
     &    ,-1)*Hr3(0,0,1)*z + 512.D0/9.D0*Hr2(0,-1)*Hr3(0,1,1) + 512.D0/
     &    9.D0*Hr2(0,-1)*Hr3(0,1,1)*z )
      G = G + as**3 * (  - 256.D0/3.D0*Hr3(0,-1,-1) + 1280.D0/9.D0*Hr3(
     &    0,-1,-1)*z**(-1) - 2912.D0/27.D0*Hr3(0,-1,-1)*z**(-1)*LL + 
     &    1328.D0/9.D0*Hr3(0,-1,-1)*LL + 128*Hr3(0,-1,-1)*z2 + 128.D0/3.
     &    D0*Hr3(0,-1,-1)*z + 1328.D0/9.D0*Hr3(0,-1,-1)*z*LL - 128*Hr3(
     &    0,-1,-1)*z*z2 + 2432.D0/9.D0*Hr3(0,-1,-1)*z**2 - 2912.D0/27.D0
     &    *Hr3(0,-1,-1)*z**2*LL + 32.D0/3.D0*Hr3(0,-1,1)*z**(-1)*LL + 
     &    32*Hr3(0,-1,1)*LL + 32*Hr3(0,-1,1)*z*LL + 32.D0/3.D0*Hr3(0,-1
     &    ,1)*z**2*LL + 128.D0/3.D0*Hr3(0,0,-1) - 640.D0/9.D0*Hr3(0,0,
     &    -1)*z**(-1) + 5344.D0/27.D0*Hr3(0,0,-1)*z**(-1)*LL + 3064.D0/
     &    3.D0*Hr3(0,0,-1)*LL + 256*Hr3(0,0,-1)*LL**2 - 192*Hr3(0,0,-1)
     &    *z2 + 1088.D0/3.D0*Hr3(0,0,-1)*z - 13912.D0/9.D0*Hr3(0,0,-1)*
     &    z*LL - 256*Hr3(0,0,-1)*z*LL**2 + 192*Hr3(0,0,-1)*z*z2 - 1216.D
     &    0/9.D0*Hr3(0,0,-1)*z**2 + 1760.D0/9.D0*Hr3(0,0,-1)*z**2*LL - 
     &    5344.D0/27.D0*Hr3(0,0,1) + 14560.D0/27.D0*Hr3(0,0,1)*z**(-1)
     &     - 5632.D0/27.D0*Hr3(0,0,1)*z**(-1)*LL - 1832.D0/9.D0*Hr3(0,0
     &    ,1)*LL )
      G = G + as**3 * (  - 832.D0/9.D0*Hr3(0,0,1)*LL**2 + 32.D0/9.D0*
     &    Hr3(0,0,1)*z2 - 256.D0/9.D0*Hr3(0,0,1)*nf*LL + 15808.D0/27.D0
     &    *Hr3(0,0,1)*z + 2504.D0/9.D0*Hr3(0,0,1)*z*LL - 832.D0/9.D0*
     &    Hr3(0,0,1)*z*LL**2 + 32.D0/9.D0*Hr3(0,0,1)*z*z2 - 256.D0/9.D0
     &    *Hr3(0,0,1)*z*nf*LL + 10688.D0/27.D0*Hr3(0,0,1)*z**2 - 4192.D0
     &    /27.D0*Hr3(0,0,1)*z**2*LL + 32.D0/3.D0*Hr3(0,1,-1)*z**(-1)*LL
     &     + 32*Hr3(0,1,-1)*LL + 32*Hr3(0,1,-1)*z*LL + 32.D0/3.D0*Hr3(0
     &    ,1,-1)*z**2*LL + 1760.D0/27.D0*Hr3(0,1,1) + 640.D0/9.D0*Hr3(0
     &    ,1,1)*z**(-1) + 1280.D0/27.D0*Hr3(0,1,1)*z**(-1)*LL - 32*Hr3(
     &    0,1,1)*LL - 352.D0/9.D0*Hr3(0,1,1)*z2 - 512.D0/27.D0*Hr3(0,1,
     &    1)*nf + 64.D0/3.D0*Hr3(0,1,1)*nf*LL + 3568.D0/27.D0*Hr3(0,1,1
     &    )*z - 992.D0/9.D0*Hr3(0,1,1)*z*LL - 352.D0/9.D0*Hr3(0,1,1)*z*
     &    z2 - 704.D0/27.D0*Hr3(0,1,1)*z*nf + 64.D0/3.D0*Hr3(0,1,1)*z*
     &    nf*LL - 640.D0/9.D0*Hr3(0,1,1)*z**2 - 256.D0/9.D0*Hr3(0,1,1)*
     &    z**2*LL - 128.D0/27.D0*Hr3(0,1,1)*z**2*nf - 128*Hr4(0,-1,-1,
     &    -1) )
      G = G + as**3 * ( 512.D0/3.D0*Hr4(0,-1,-1,-1)*z**(-1) - 128*Hr4(0
     &    ,-1,-1,-1)*z + 512.D0/3.D0*Hr4(0,-1,-1,-1)*z**2 - 64*Hr4(0,-1
     &    ,0,1) + 256.D0/3.D0*Hr4(0,-1,0,1)*z**(-1) - 64*Hr4(0,-1,0,1)*
     &    z + 256.D0/3.D0*Hr4(0,-1,0,1)*z**2 + 64*Hr4(0,0,-1,-1) - 256.D
     &    0/3.D0*Hr4(0,0,-1,-1)*z**(-1) + 64*Hr4(0,0,-1,-1)*z - 256.D0/
     &    3.D0*Hr4(0,0,-1,-1)*z**2 - 128*Hr4(0,0,-1,1) + 512.D0/3.D0*
     &    Hr4(0,0,-1,1)*z**(-1) - 128*Hr4(0,0,-1,1)*z + 512.D0/3.D0*
     &    Hr4(0,0,-1,1)*z**2 + 32*Hr4(0,0,0,-1) - 128.D0/3.D0*Hr4(0,0,0
     &    ,-1)*z**(-1) + 736.D0/3.D0*Hr4(0,0,0,-1)*LL + 32*Hr4(0,0,0,-1
     &    )*z + 1888.D0/3.D0*Hr4(0,0,0,-1)*z*LL - 128.D0/3.D0*Hr4(0,0,0
     &    ,-1)*z**2 - 5728.D0/9.D0*Hr4(0,0,0,1) + 7424.D0/27.D0*Hr4(0,0
     &    ,0,1)*z**(-1) + 992.D0/9.D0*Hr4(0,0,0,1)*LL + 5152.D0/9.D0*
     &    Hr4(0,0,0,1)*z - 9376.D0/9.D0*Hr4(0,0,0,1)*z*LL + 2560.D0/27.D
     &    0*Hr4(0,0,0,1)*z**2 - 128*Hr4(0,0,1,-1) + 512.D0/3.D0*Hr4(0,0
     &    ,1,-1)*z**(-1) - 128*Hr4(0,0,1,-1)*z + 512.D0/3.D0*Hr4(0,0,1,
     &    -1)*z**2 )
      G = G + as**3 * (  - 160.D0/9.D0*Hr4(0,0,1,1) - 3712.D0/27.D0*
     &    Hr4(0,0,1,1)*z**(-1) - 64*Hr4(0,0,1,1)*LL + 1984.D0/9.D0*Hr4(
     &    0,0,1,1)*z - 64*Hr4(0,0,1,1)*z*LL + 640.D0/3.D0*Hr4(0,0,1,1)*
     &    z**2 + 64*Hr4(0,1,1,1) + 640.D0/9.D0*Hr4(0,1,1,1)*z**(-1) + 
     &    320.D0/9.D0*Hr4(0,1,1,1)*LL + 64.D0/9.D0*Hr4(0,1,1,1)*nf - 64.
     &    D0/9.D0*Hr4(0,1,1,1)*z + 320.D0/9.D0*Hr4(0,1,1,1)*z*LL + 64.D0
     &    /9.D0*Hr4(0,1,1,1)*z*nf - 1408.D0/27.D0*Hr4(0,1,1,1)*z**2 + 
     &    512*Hr5(0,-1,0,-1,-1) - 512*Hr5(0,-1,0,-1,-1)*z + 1024*Hr5(0,
     &    0,-1,-1,-1) - 1024*Hr5(0,0,-1,-1,-1)*z + 128*Hr5(0,0,-1,0,-1)
     &     - 128*Hr5(0,0,-1,0,-1)*z - 256*Hr5(0,0,-1,0,1) + 256*Hr5(0,0
     &    ,-1,0,1)*z + 384*Hr5(0,0,0,-1,-1) - 384*Hr5(0,0,0,-1,-1)*z - 
     &    768*Hr5(0,0,0,-1,1) + 768*Hr5(0,0,0,-1,1)*z + 256*Hr5(0,0,0,0
     &    ,-1) - 256*Hr5(0,0,0,0,-1)*z - 256.D0/3.D0*Hr5(0,0,0,0,1) + 
     &    3584.D0/3.D0*Hr5(0,0,0,0,1)*z - 768*Hr5(0,0,0,1,-1) + 768*
     &    Hr5(0,0,0,1,-1)*z - 3136.D0/3.D0*Hr5(0,0,0,1,1) - 3136.D0/3.D0
     &    *Hr5(0,0,0,1,1)*z )
      G = G + as**3 * (  - 256*Hr5(0,0,1,0,-1) + 256*Hr5(0,0,1,0,-1)*z
     &     - 2624.D0/9.D0*Hr5(0,0,1,0,1) - 2624.D0/9.D0*Hr5(0,0,1,0,1)*
     &    z - 1024.D0/9.D0*Hr5(0,0,1,1,1) - 1024.D0/9.D0*Hr5(0,0,1,1,1)
     &    *z - 320.D0/3.D0*Hr5(0,1,0,1,1) - 320.D0/3.D0*Hr5(0,1,0,1,1)*
     &    z + 320.D0/9.D0*Hr5(0,1,1,1,1) + 320.D0/9.D0*Hr5(0,1,1,1,1)*z
     &     - 102024/( - 243 + 243*z)*z + 7872/( - 243 + 243*z)*z*nf + 
     &    3936/( - 81 + 81*z)*z - 5580/( - 81 + 81*z)*z*LL - 1632/( - 
     &    81 + 81*z)*z*nf*LL - 480/( - 27 + 27*z)*z*LL**2 + 624/( - 27
     &     + 27*z)*z*LL**3 + 492/( - 27 + 27*z)*z*z2 - 48/( - 27 + 27*z
     &    )*z*nf*LL**3 + 120/( - 27 + 27*z)*z*nf*z2 - 624/( - 27 + 27*z
     &    )*z*z3 + 48/( - 27 + 27*z)*z*z3*nf - 932/( - 9 + 9*z)*z - 240
     &    /( - 9 + 9*z)*z*LL - 828/( - 9 + 9*z)*z*LL**2 + 1440/( - 9 + 
     &    9*z)*z*z2*LL - 396/( - 9 + 9*z)*Hr1(0)*z + 24/( - 9 + 9*z)*
     &    Hr1(0)*z*nf - 720/( - 9 + 9*z)*Hr1(0)**2*z*LL - 2880/( - 9 + 
     &    9*z)*Hr1(0)*Hr1(1)*z*LL - 80/( - 3 + 3*z)*z*LL + 144/( - 3 + 
     &    3*z)*z*z2*LL**2 )
      G = G + as**3 * (  - 72/( - 3 + 3*z)*z*z2**2 - 1152/( - 3 + 3*z)*
     &    z*z3*LL + 12/( - 3 + 3*z)*Hr1(0)*z - 72/( - 3 + 3*z)*Hr1(0)*z
     &    *LL - 72/( - 3 + 3*z)*Hr1(0)**2*z*LL**2 + 36/( - 3 + 3*z)*
     &    Hr1(0)**2*z*z2 - 72/( - 3 + 3*z)*Hr1(0)**2*Hr1(1)*z*LL - 288
     &    /( - 3 + 3*z)*Hr1(0)*Hr1(1)*z*LL**2 + 144/( - 3 + 3*z)*Hr1(0)
     &    *Hr1(1)*z*z2 - 144/( - 3 + 3*z)*Hr1(0)*Hr2(0,-1)*z*LL + 576/(
     &     - 3 + 3*z)*Hr3(0,0,-1)*z*LL - 288/( - 3 + 3*z)*Hr3(0,0,1)*z*
     &    LL - 16/( - 1 + z)*z*LL**2 - 80/( - 1 + z)*z*z2 + 128/( - 1
     &     + z)*z*z3*LL - 16/(1 - z)*LL**2 - 80/(1 - z)*z2 + 128/(1 - z
     &    )*z3*LL - 80/(3 - 3*z)*LL + 144/(3 - 3*z)*z2*LL**2 - 72/(3 - 
     &    3*z)*z2**2 - 1008/(3 - 3*z)*z3*LL + 144/(3 + 3*z)*z*z2*LL**2
     &     - 72/(3 + 3*z)*z*z2**2 + 288/(3 + 3*z)*Hr1(-1)*Hr1(0)*z*
     &    LL**2 - 144/(3 + 3*z)*Hr1(-1)*Hr1(0)*z*z2 + 144/(3 + 3*z)*
     &    Hr1(0)*z*z2*LL - 72/(3 + 3*z)*Hr1(0)**2*z*LL**2 + 36/(3 + 3*z
     &    )*Hr1(0)**2*z*z2 - 288/(3 + 3*z)*Hr2(0,-1)*z*LL**2 + 144/(3
     &     + 3*z)*Hr2(0,-1)*z*z2 )
      G = G + as**3 * (  - 932/(9 - 9*z) - 240/(9 - 9*z)*LL - 828/(9 - 
     &    9*z)*LL**2 + 1440/(9 - 9*z)*z2*LL + 1440/(9 + 9*z)*z*z2*LL + 
     &    2880/(9 + 9*z)*Hr1(-1)*Hr1(0)*z*LL - 720/(9 + 9*z)*Hr1(0)**2*
     &    z*LL - 2880/(9 + 9*z)*Hr2(0,-1)*z*LL - 480/(27 - 27*z)*LL**2
     &     + 624/(27 - 27*z)*LL**3 + 492/(27 - 27*z)*z2 - 48/(27 - 27*z
     &    )*nf*LL**3 + 120/(27 - 27*z)*nf*z2 - 624/(27 - 27*z)*z3 + 48
     &    /(27 - 27*z)*z3*nf + 3936/(81 - 81*z) - 5580/(81 - 81*z)*LL
     &     - 1632/(81 - 81*z)*nf*LL - 102024/(243 - 243*z) + 7872/(243
     &     - 243*z)*nf )
*
      RETURN
      END      
      REAL*8 FUNCTION GGPLU(z,nf,as,LL)
*     ---------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     plus part of the OME AggQ without aggQ3 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(6)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      w(1)=1/(1 - z)
      w(2)=11.D0/9.D0 + 14*LL
      w(2)=z3*w(2)
      w(3)=LL + 10.D0/3.D0
      w(3)=w(3)*LL
      w(4)=z2 - 1.D0/9.D0 - 2*w(3)
      w(4)=z2*w(4)
      w(2)=w(2) + w(4)
      w(4)=23 - 22.D0/3.D0*LL
      w(4)=LL*w(4)
      w(4)=155.D0/9.D0 + w(4)
      w(4)=LL*w(4)
      w(4)=2834.D0/27.D0 + w(4)
      w(2)=1.D0/3.D0*w(4) + 2*w(2)
      w(4)=4.D0/3.D0*CA
      w(2)=w(2)*w(4)
      w(5)=20 + 7*LL
      w(5)=LL*w(5)
      w(5)=10 + 1.D0/3.D0*w(5)
      w(5)=LL*w(5)
      w(5)= - 7.D0/3.D0*z3 - 164.D0/9.D0 + w(5)
      w(6)=LL**2
      w(6)=34.D0/3.D0 + w(6)
      w(6)=LL*w(6)
      w(6)= - z3 - 164.D0/9.D0 + w(6)
      w(6)=2*w(6) - 5*z2
      w(6)=nf*w(6)
      w(5)=2.D0/3.D0*w(6) + 2*w(5) - 35.D0/3.D0*z2
      w(5)=w(5)*TF
      w(6)=5.D0/3.D0 + LL
      w(6)= - 32*z3 + 4*w(6)
      w(6)=LL*w(6)
      w(6)=20*z2 + 233.D0/9.D0 + w(6)
      w(6)=CF*w(6)
      w(2)=8.D0/9.D0*w(5) + w(6) + w(2)
      w(2)=as*CA*w(2)
      w(3)=28.D0/9.D0 + w(3)
      w(3)=w(3)*w(4)
      w(2)=w(3) + w(2)
*
      GGPLU = 2*TF*as**2*w(2)*w(1)
*
      RETURN
      END      
      REAL*8 FUNCTION GGDEL(z,nf,as,LL)
*     ---------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Delta[1-z] part of the OME AggQ without aggQ3 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(7)
*
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5 
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5 
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5 
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      w(1)=4*LL
      w(2)=2*LL
      w(3)=z2*w(2)
      w(3)=w(3) + 31.D0/3.D0 + w(1)
      w(3)=z2*w(3)
      w(4)=277.D0/3.D0 - w(2)
      w(4)=LL*w(4)
      w(3)=4*w(3) - 616.D0/9.D0 + w(4)
      w(4)=10.D0/3.D0 + LL
      w(4)=w(4)*w(2)
      w(4)=w(4) - z2
      w(4)=z3*w(4)
      w(3)=1.D0/3.D0*w(3) + 8*w(4)
      w(3)=w(3)*CA**2
      w(4)=11*LL
      w(5)=368.D0/3.D0 - w(4)
      w(5)=w(5)*w(2)
      w(6)=20*z2
      w(5)=w(6) - 1045.D0/2.D0 + w(5)
      w(7)=ln2*z2
      w(5)= - 64*w(7) + 1.D0/3.D0*w(5) + 16*z3
      w(5)=CA*w(5)
      w(2)=128*w(7) - 32*z3 - 80*z2 - 39 - w(2)
      w(2)=CF*w(2)
      w(2)=w(5) + w(2)
      w(2)=CF*w(2)
      w(5)= - 1 + 28.D0/3.D0*LL
      w(5)=LL*w(5)
      w(4)= - 11.D0/3.D0*z2 + 56.D0/9.D0 - w(4)
      w(4)=nf*w(4)
      w(4)=2.D0/3.D0*w(4) - 10*z2 - 4.D0/27.D0 + w(5)
      w(4)=CA*w(4)
      w(5)= - 73.D0/3.D0 + 5*LL
      w(5)=w(5)*w(1)
      w(5)= - w(6) + 391.D0/3.D0 + w(5)
      w(6)=59 - 134.D0/3.D0*LL
      w(6)=1.D0/3.D0*w(6) + 14*z2
      w(6)=nf*w(6)
      w(5)=1.D0/3.D0*w(5) + w(6)
      w(5)=CF*w(5)
      w(6)=LL**3
      w(6)=w(6) - z3
      w(6)=TF*w(6)
      w(4)=32.D0/27.D0*w(6) + w(4) + w(5)
      w(4)=TF*w(4)
      w(2)=2*w(4) + w(3) + w(2)
      w(2)=as*w(2)
      w(3)=5.D0/3.D0 + 8*LL
      w(3)=CA*w(3)
      w(1)= - 15 + w(1)
      w(1)=CF*w(1)
      w(4)=TF*LL**2
      w(1)=w(2) + 16.D0/9.D0*w(4) + 2.D0/3.D0*w(3) + w(1)
      w(1)=as*w(1)
      w(1)=4.D0/3.D0*LL + w(1)
      w(1)=as*TF*w(1)
*    
      GGDEL = 1.0D0 + w(1) 
*
      RETURN
      END      
      REAL*8 FUNCTION PS(z,nf,as,LL)
*     ------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of the OME AQq without aQqPS3 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(138)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
      w(1)=z**(-1)
      w(2)=Hr1(-1)
      w(3)=Hr1(0)
      w(4)=Hr2(-1,-1)
      w(5)=Hr3(-1,-1,-1)
      w(6)=Hr4(-1,-1,0,-1)
      w(7)=Hr3(-1,0,-1)
      w(8)=Hr4(-1,0,-1,-1)
      w(9)=Hr4(-1,0,0,-1)
      w(10)=Hr4(-1,0,0,1)
      w(11)=Hr3(-1,0,1)
      w(12)=Hr2(0,-1)
      w(13)=Hr2(0,1)
      w(14)=Hr1(1)
      w(15)=Hr3(0,-1,-1)
      w(16)=Hr3(0,0,-1)
      w(17)=Hr3(0,0,1)
      w(18)=Hr3(0,1,1)
      w(19)=Hr4(0,-1,-1,-1)
      w(20)=Hr4(0,-1,0,-1)
      w(21)=Hr4(0,-1,0,1)
      w(22)=Hr4(0,0,-1,-1)
      w(23)=Hr4(0,0,0,-1)
      w(24)=Hr4(0,0,0,1)
      w(25)=Hr4(0,0,1,1)
      w(26)=Hr4(0,1,0,1)
      w(27)=Hr4(0,1,1,1)
      w(28)=Hr5(0,-1,-1,0,-1)
      w(29)=Hr5(0,-1,0,-1,-1)
      w(30)=Hr5(0,-1,0,0,-1)
      w(31)=Hr5(0,-1,0,0,1)
      w(32)=Hr5(0,0,-1,-1,-1)
      w(33)=Hr5(0,0,-1,0,-1)
      w(34)=Hr5(0,0,-1,0,1)
      w(35)=Hr5(0,0,0,-1,-1)
      w(36)=Hr5(0,0,0,0,-1)
      w(37)=Hr5(0,0,0,0,1)
      w(38)=Hr5(0,0,0,1,1)
      w(39)=Hr5(0,0,1,0,1)
      w(40)=Hr5(0,0,1,1,1)
      w(41)=Hr5(0,1,0,0,1)
      w(42)=Hr5(0,1,0,1,1)
      w(43)=Hr5(0,1,1,0,1)
      w(44)=Hr5(0,1,1,1,1)
      w(45)=z - 1
      w(46)=2*w(3)
      w(47)=w(45)*w(46)
      w(48)=4*w(1)
      w(49)=w(48) + 11
      w(50)=4*z
      w(51)= - 17 - w(50)
      w(51)=z*w(51)
      w(51)=w(51) - w(49)
      w(52)=z + 1
      w(53)=w(52)*LL
      w(54)=2.D0/3.D0*w(53)
      w(51)= - w(54) + 1.D0/3.D0*w(51) - w(47)
      w(55)=2*LL
      w(51)=w(51)*w(55)
      w(56)=w(45)*w(3)
      w(57)=8*z
      w(58)= - 13 + w(57)
      w(58)=z*w(58)
      w(58)=59 + w(58)
      w(58)=1.D0/3.D0*w(58) + w(56)
      w(58)=w(3)*w(58)
      w(59)=4*w(13)
      w(60)= - w(52)*w(59)
      w(61)=w(45)*w(12)
      w(62)= - 313.D0/3.D0 - 44*z
      w(62)=z*w(62)
      w(62)= - 56*w(1) + 473.D0/3.D0 + w(62)
      w(51)=w(51) + 8*w(61) + w(60) + 1.D0/3.D0*w(62) + w(58)
      w(51)=w(51)*w(55)
      w(58)=4.D0/3.D0*z
      w(60)=w(58) + 1
      w(62)=w(60)*z
      w(63)=4.D0/3.D0*w(1)
      w(64)=w(62) - w(63) - 1
      w(65)=w(64)*w(46)
      w(66)=w(64)*w(14)
      w(67)=4*LL
      w(68)=w(64)*w(67)
      w(69)=1 + 154.D0/9.D0*z
      w(69)=z*w(69)
      w(68)=w(66) + w(68) - w(65) - 46.D0/9.D0*w(1) - 13 + w(69)
      w(68)=w(14)*w(68)
      w(69)=19.D0/3.D0*z
      w(70)=w(69) + 1
      w(70)=w(70)*z
      w(71)=10.D0/3.D0*w(1)
      w(72)=w(71) - 2
      w(70)=w(70) + w(72)
      w(73)=4.D0/3.D0*w(70)
      w(74)=w(50) - 3
      w(74)=w(74)*z
      w(75)=w(48) - 3
      w(74)=w(74) + w(75)
      w(76)=w(3)*w(74)
      w(76)= - 12*w(53) - w(73) + w(76)
      w(77)=2*w(2)
      w(76)=w(76)*w(77)
      w(78)=17 + z
      w(78)=w(78)*w(50)
      w(78)=152.D0/5.D0*w(1) - 137.D0/5.D0 + w(78)
      w(79)= - 1 + 38.D0/5.D0*z
      w(79)=w(79)*w(46)
      w(80)= - 87 - 209*z
      w(80)=LL*w(80)
      w(78)=1.D0/5.D0*w(80) + 1.D0/3.D0*w(78) + w(79)
      w(79)=2*z2
      w(78)=w(78)*w(79)
      w(80)=16*w(15)
      w(81)= - 24*w(16) + w(80)
      w(81)=w(45)*w(81)
      w(82)=10*z
      w(83)=17 - w(82)
      w(83)=1.D0/3.D0*w(83) + w(56)
      w(83)=w(83)*w(46)
      w(84)= - 121.D0/3.D0 - 158*z
      w(84)=z*w(84)
      w(84)= - 80.D0/3.D0*w(1) - 295.D0/3.D0 + w(84)
      w(83)=1.D0/3.D0*w(84) + w(83)
      w(83)=w(3)*w(83)
      w(84)=w(46)*w(52)
      w(85)=23*z
      w(86)= - 14 - w(85)
      w(86)=1.D0/3.D0*w(86) + w(84)
      w(87)=2*w(13)
      w(86)=w(86)*w(87)
      w(88)=w(58) - 1
      w(88)=w(88)*z
      w(89)=w(63) - 1
      w(90)=w(88) + w(89)
      w(91)=8*w(4)
      w(92)=w(90)*w(91)
      w(93)=6*w(56)
      w(94)=w(93) - w(74)
      w(94)=w(12)*w(94)
      w(95)=w(52)*w(17)
      w(96)= - 191 + 1841.D0/9.D0*z
      w(96)=z*w(96)
      w(96)= - 515.D0/9.D0*w(1) + 410.D0/3.D0 + w(96)
      w(51)=w(78) - 4*w(95) + w(76) + w(68) + w(51) + 2*w(94) + w(92)
     &  + w(86) + 2.D0/3.D0*w(96) + w(83) + w(81)
      w(51)=w(51)*w(79)
      w(68)=40*w(1)
      w(76)=103 + 272*z
      w(76)=z*w(76)
      w(76)=w(68) + 139 + w(76)
      w(78)=4*w(3)
      w(81)=w(78)*w(45)
      w(83)=16*z
      w(86)= - w(81) - 17 + w(83)
      w(86)=w(3)*w(86)
      w(76)=1.D0/3.D0*w(76) + w(86)
      w(86)=1.D0/3.D0*w(3)
      w(76)=w(76)*w(86)
      w(92)=w(90) - w(47)
      w(94)=4*w(12)
      w(96)=w(92)*w(94)
      w(49)=w(57) + w(49)
      w(97)=2*z
      w(98)=w(97) - 1
      w(99)=w(98)*w(3)
      w(49)=1.D0/3.D0*w(49) + w(99)
      w(49)=w(3)*w(49)
      w(100)=w(52)*w(13)
      w(101)=5 - 44.D0/9.D0*z
      w(101)=z*w(101)
      w(49)=4.D0/3.D0*w(100) + 2.D0/3.D0*w(49) + 44.D0/9.D0*w(1) - 5 + 
     & w(101)
      w(49)=LL*w(49)
      w(101)=w(52)*w(3)
      w(102)=7 + w(82)
      w(102)=1.D0/3.D0*w(102) - w(101)
      w(102)=w(102)*w(59)
      w(103)=87 - 1864.D0/27.D0*z
      w(103)=z*w(103)
      w(49)=w(49) + w(96) + w(102) + w(76) + 694.D0/27.D0*w(1) - 131.D0/
     & 3.D0 + w(103)
      w(49)=w(49)*w(55)
      w(76)=173.D0/3.D0 - 101*z
      w(76)=w(76)*w(97)
      w(96)=46 - 35*z
      w(102)= - 5 + w(50)
      w(102)=w(3)*w(102)
      w(96)=2.D0/3.D0*w(96) + w(102)
      w(96)=w(3)*w(96)
      w(76)=w(96) - 377.D0/3.D0 + w(76)
      w(76)=w(3)*w(76)
      w(96)=2093.D0/3.D0 + 1243*z
      w(96)=z*w(96)
      w(96)=224.D0/3.D0*w(1) + 2330.D0/3.D0 + w(96)
      w(76)=4.D0/3.D0*w(96) + w(76)
      w(76)=w(3)*w(76)
      w(96)=11143 - 11542*z
      w(96)=z*w(96)
      w(96)=3892*w(1) - 3493 + w(96)
      w(76)=2.D0/9.D0*w(96) + w(76)
      w(96)= - 13 + w(50)
      w(96)=z*w(96)
      w(47)=w(47) + w(96) - w(75)
      w(47)=w(3)*w(47)
      w(75)=25.D0/3.D0*z
      w(96)=w(75) + 7
      w(102)=5*z
      w(96)=w(96)*w(102)
      w(96)=w(96) + 14 + 62.D0/3.D0*w(1)
      w(47)=2.D0/3.D0*w(96) + w(47)
      w(47)=w(47)*w(94)
      w(103)= - w(92)*w(80)
      w(104)=2.D0/3.D0*z
      w(105)=w(104) + 1
      w(105)=w(105)*z
      w(106)=1 + 2.D0/3.D0*w(1)
      w(107)=w(105) + w(106)
      w(108)=16*w(11)
      w(107)=w(107)*w(108)
      w(108)=8*w(24)
      w(109)=5 + w(85)
      w(109)=w(109)*w(108)
      w(110)= - 41 + w(83)
      w(110)=z*w(110)
      w(110)= - w(68) - 77 + w(110)
      w(110)=1.D0/3.D0*w(110) + 7*w(101)
      w(110)=w(3)*w(110)
      w(111)=205 - 274*z
      w(111)=z*w(111)
      w(111)=200*w(1) - 659 + w(111)
      w(110)=1.D0/9.D0*w(111) + w(110)
      w(110)=w(110)*w(87)
      w(111)=w(90)*w(3)
      w(112)=16*w(111)
      w(113)=w(4)*w(112)
      w(114)=16*w(90)
      w(115)= - w(7)*w(114)
      w(116)=64*w(45)
      w(117)= - w(22)*w(116)
      w(47)=w(49) + w(103) + w(115) + w(47) + w(117) + w(113) + w(110)
     &  + w(109) + 1.D0/3.D0*w(76) + w(107)
      w(47)=w(47)*w(55)
      w(49)=w(81) - w(90)
      w(49)=w(22)*w(49)
      w(76)=w(55)*w(52)
      w(81)=3*w(101)
      w(103)=w(76) + w(81) - w(64)
      w(103)=w(26)*w(103)
      w(107)=w(67)*w(45)
      w(109)= - w(107) - w(92)
      w(109)=w(20)*w(109)
      w(110)=w(84) - w(53)
      w(113)= - w(105) + w(106) + w(110)
      w(113)=w(27)*w(113)
      w(49)=w(49) + w(103) + w(109) + w(113)
      w(79)=w(79)*w(52)
      w(103)=7*z
      w(109)=8.D0/3.D0*z
      w(113)= - 5 + w(109)
      w(113)=w(113)*w(103)
      w(113)= - 40.D0/3.D0*w(1) - 22 + w(113)
      w(115)=5*w(64) - w(101)
      w(115)=w(115)*w(46)
      w(117)=8*w(1)
      w(118)=23 + w(57)
      w(118)=z*w(118)
      w(118)= - w(117) + 5 + w(118)
      w(118)=1.D0/3.D0*w(118) - 6*w(101)
      w(118)=w(118)*w(55)
      w(113)= - w(79) + 2*w(66) + w(118) + 1.D0/3.D0*w(113) + w(115)
      w(115)=4*w(18)
      w(113)=w(113)*w(115)
      w(118)=2*w(1)
      w(119)=6*z
      w(120)=13 - w(119)
      w(120)=z*w(120)
      w(121)=3*z
      w(122)= - 1 - w(121)
      w(122)=w(3)*w(122)
      w(107)=w(107) + w(122) + w(118) - 3 + w(120)
      w(107)=w(107)*w(55)
      w(120)= - w(3)*w(92)
      w(122)= - 17 + w(69)
      w(122)=z*w(122)
      w(122)=w(122) + w(72)
      w(107)=w(107) + 2.D0/3.D0*w(122) + w(120)
      w(107)=w(16)*w(107)
      w(120)=3 + z
      w(120)=LL*w(120)
      w(93)=6*w(120) - w(93) + w(90)
      w(93)=w(23)*w(93)
      w(120)= - 5 - w(119)
      w(120)=w(120)*w(97)
      w(122)=4*w(101)
      w(120)=10*w(53) + w(122) + 12*w(1) + 9 + w(120)
      w(120)=w(25)*w(120)
      w(93)=w(120) + w(107) + w(93)
      w(107)=w(64)*LL
      w(120)=2.D0/3.D0*w(107)
      w(123)= - 5 - 104.D0/3.D0*z
      w(123)=z*w(123)
      w(123)=32.D0/3.D0*w(1) + 29 + w(123)
      w(123)= - w(120) + 1.D0/3.D0*w(123) + w(65)
      w(123)=w(123)*w(55)
      w(124)=w(117) + 9
      w(125)= - 9 - w(57)
      w(125)=z*w(125)
      w(125)=w(125) + w(124)
      w(125)=w(3)*w(125)
      w(126)= - 52 + 203.D0/3.D0*z
      w(126)=z*w(126)
      w(126)= - 140.D0/3.D0*w(1) + 31 + w(126)
      w(125)=4.D0/3.D0*w(126) + w(125)
      w(125)=w(3)*w(125)
      w(126)=8*w(13)
      w(127)= - w(64)*w(126)
      w(128)=769.D0/3.D0 - 70*z
      w(128)=z*w(128)
      w(128)= - 62*w(1) - 373.D0/3.D0 + w(128)
      w(123)=w(123) + w(127) + 2.D0/3.D0*w(128) + w(125)
      w(123)=LL*w(123)
      w(125)=17 + 500.D0/3.D0*z
      w(125)=w(125)*w(102)
      w(125)= - 1771.D0/3.D0*w(1) - 328 + w(125)
      w(127)=277 - 769.D0/3.D0*z
      w(127)=z*w(127)
      w(127)=337.D0/3.D0*w(1) - 133 + w(127)
      w(128)= - 47 + 122.D0/3.D0*z
      w(128)=z*w(128)
      w(128)= - 113.D0/3.D0*w(1) + 44 + w(128)
      w(128)=w(3)*w(128)
      w(127)=4.D0/3.D0*w(127) + w(128)
      w(127)=w(127)*w(46)
      w(125)=1.D0/9.D0*w(125) + w(127)
      w(127)=w(50) + 3
      w(127)=w(127)*z
      w(128)= - w(48) + w(127) - 3
      w(129)= - w(3)*w(128)*w(59)
      w(123)=w(123) + 1.D0/3.D0*w(125) + w(129)
      w(125)=w(64)*w(3)
      w(129)=1 - w(69)
      w(129)=z*w(129)
      w(71)=w(71) + 2 + w(129)
      w(71)=4.D0/3.D0*w(71) + w(125)
      w(71)=w(3)*w(71)
      w(129)=20*z
      w(130)= - 1 - w(129)
      w(130)=z*w(130)
      w(130)=w(118) + 19 + w(130)
      w(128)=w(128)*w(46)
      w(128)=1.D0/3.D0*w(130) + w(128)
      w(128)=LL*w(128)
      w(130)= - 5 + 328.D0/3.D0*z
      w(130)=z*w(130)
      w(130)= - 67.D0/3.D0*w(1) - 82 + w(130)
      w(71)=w(128) + 1.D0/9.D0*w(130) + w(71)
      w(128)=w(65) - w(107)
      w(130)= - 1 + 7.D0/3.D0*z
      w(130)=w(130)*w(50)
      w(131)=1.D0/3.D0*w(1)
      w(130)= - w(131) - 5 + w(130)
      w(130)=1.D0/3.D0*w(130) - w(128)
      w(130)=4*w(130) + w(66)
      w(132)=1.D0/3.D0*w(14)
      w(130)=w(130)*w(132)
      w(71)=2*w(71) + w(130)
      w(71)=w(14)*w(71)
      w(71)=2*w(123) + w(71)
      w(71)=w(14)*w(71)
      w(123)=11 - w(50)
      w(123)=z*w(123)
      w(123)=14*w(1) + 17 + w(123)
      w(130)= - 5 - 11*z
      w(130)=w(3)*w(130)
      w(123)=w(53) + 1.D0/3.D0*w(123) + w(130)
      w(123)=w(123)*w(55)
      w(130)=5*w(66)
      w(133)= - 3 - 10.D0/9.D0*z
      w(85)=w(133)*w(85)
      w(82)=w(82) + 1
      w(133)=w(48) + w(82)
      w(99)=2.D0/3.D0*w(133) - w(99)
      w(99)=w(99)*w(46)
      w(85)=w(130) + w(123) + w(99) - 193.D0/9.D0*w(1) - 56.D0/3.D0 + 
     & w(85)
      w(99)=8*w(17)
      w(85)=w(85)*w(99)
      w(123)=1.D0/3.D0*z
      w(133)= - 227 + 2132.D0/3.D0*z
      w(133)=w(133)*w(123)
      w(133)=95 + w(133)
      w(134)=17 - w(50)
      w(134)=1.D0/3.D0*w(134) + 2.D0/5.D0*w(56)
      w(134)=w(134)*w(86)
      w(135)= - 47.D0/9.D0 - w(50)
      w(135)=z*w(135)
      w(135)= - 4 + w(135)
      w(134)=2*w(135) + w(134)
      w(134)=w(3)*w(134)
      w(133)=2.D0/3.D0*w(133) + w(134)
      w(133)=w(3)*w(133)
      w(73)= - w(73) - w(111)
      w(73)=w(73)*w(91)
      w(91)= - 3113 - 18068*z
      w(91)=w(91)*w(123)
      w(91)= - 2624.D0/9.D0*w(1) - 1898 + w(91)
      w(73)=w(73) + 2.D0/9.D0*w(91) + w(133)
      w(73)=w(3)*w(73)
      w(91)=2*w(70) - w(111)
      w(91)=w(3)*w(91)
      w(133)=10 + 91.D0/3.D0*z
      w(133)=w(133)*w(97)
      w(133)=w(133) - 25 + 47.D0/3.D0*w(1)
      w(91)=4.D0/3.D0*w(133) + w(91)
      w(91)=w(91)*w(86)
      w(134)=z**2
      w(135)=w(134) + w(1)
      w(135)=w(3)*w(135)
      w(96)= - 1.D0/3.D0*w(96) + w(135)
      w(96)=w(3)*w(96)
      w(135)= - LL*w(111)
      w(96)=w(96) + w(135)
      w(96)=w(96)*w(67)
      w(91)=w(91) + w(96)
      w(91)=w(2)*w(91)
      w(96)=16*w(21) - 32*w(19)
      w(92)=w(92)*w(96)
      w(96)= - 2 + w(102)
      w(96)=w(96)*w(78)
      w(135)=5 - w(103)
      w(96)=w(96) + 3*w(135) - w(117)
      w(96)=w(96)*w(108)
      w(135)=13*z
      w(136)= - 29 + 76.D0/3.D0*z
      w(136)=w(136)*w(135)
      w(136)= - 1348.D0/9.D0*w(1) + 152 + w(136)
      w(137)= - 19*w(52) - w(117)
      w(137)=w(137)*w(86)
      w(119)=127.D0/3.D0 + w(119)
      w(119)=z*w(119)
      w(119)=17*w(1) + 2 + w(119)
      w(119)=2*w(119) + w(137)
      w(119)=w(3)*w(119)
      w(119)=1.D0/3.D0*w(136) + w(119)
      w(119)=w(119)*w(59)
      w(69)=8 - w(69)
      w(69)=z*w(69)
      w(69)=w(69) - w(72)
      w(72)= - 2.D0/3.D0*w(56) + w(90)
      w(72)=w(3)*w(72)
      w(69)=4.D0/3.D0*w(69) + w(72)
      w(69)=w(3)*w(69)
      w(69)= - 4.D0/9.D0*w(133) + w(69)
      w(69)=w(69)*w(94)
      w(70)=2.D0/3.D0*w(70)
      w(72)=w(70) + w(111)
      w(94)=16*w(7)
      w(72)=w(72)*w(94)
      w(56)= - w(56) + w(90)
      w(56)=w(3)*w(56)
      w(56)=w(70) + w(56)
      w(56)=w(56)*w(80)
      w(70)=32*w(90)
      w(80)= - w(8) - w(6) + w(10)
      w(70)=w(70)*w(80)
      w(80)= - 16*w(42) - 80*w(41)
      w(80)=w(52)*w(80)
      w(90)= - 128*w(32) - 32*w(30) - 96*w(35)
      w(90)=w(45)*w(90)
      w(94)=w(34) - w(29) + w(31) + w(36) - w(28) - w(33)
      w(94)=w(116)*w(94)
      w(116)= - 41585.D0/3.D0 + 12586*z
      w(116)=w(116)*w(97)
      w(116)= - 7255*w(1) + 29419.D0/3.D0 + w(116)
      w(133)=w(52)*w(40)
      w(116)=1.D0/27.D0*w(116) - 32*w(133)
      w(112)= - w(11)*w(112)
      w(114)= - w(9)*w(114)
      w(111)=w(5)*w(111)
      w(136)=w(50) - 1
      w(137)=w(37)*w(136)
      w(138)= - 3 + w(103)
      w(138)=z5*w(138)
      w(47)=40*w(138) + 32*w(111) + w(113) + w(51) + w(85) + 4*w(91) + 
     & w(114) - 64*w(137) + w(71) + w(47) + w(56) + w(72) + w(69) + 
     & w(119) + w(96) + w(112) + 2*w(116) + w(94) + w(90) + w(92) + 
     & w(80) + w(73) + w(70) + 8*w(93) + 16*w(49)
      w(49)=as**3
      w(51)=w(49)*CA
      w(47)=w(47)*w(51)
      w(56)=1.D0/3.D0*w(53) + w(101) + 2 + w(121)
      w(56)=w(56)*w(55)
      w(69)= - 1 - 37.D0/9.D0*z
      w(69)=w(69)*w(97)
      w(70)=3 + w(109)
      w(70)=w(70)*w(97)
      w(70)= - w(81) - 3 + w(70)
      w(70)=w(70)*w(46)
      w(56)=w(56) + w(70) - 33 + w(69)
      w(56)=LL*w(56)
      w(69)=w(64)*w(78)
      w(70)=40*z
      w(71)=7 - w(70)
      w(71)=z*w(71)
      w(71)=w(48) + 29 + w(71)
      w(71)= - w(130) + 1.D0/3.D0*w(71) - w(69)
      w(72)=2*w(14)
      w(71)=w(71)*w(72)
      w(73)=8.D0/3.D0*w(1)
      w(80)= - 5 - w(50)
      w(80)=z*w(80)
      w(80)=w(122) - w(73) - 3 + w(80)
      w(80)=w(80)*w(59)
      w(85)= - 9 + 188.D0/3.D0*z
      w(85)=z*w(85)
      w(85)=36.D0/5.D0*w(53) - 56.D0/5.D0*w(101) - 7 + 1.D0/5.D0*w(85)
      w(85)=z2*w(85)
      w(90)=20*w(1)
      w(91)=79.D0/3.D0 + 28*z
      w(91)=z*w(91)
      w(91)= - w(90) + 25 + w(91)
      w(92)=1 - w(88)
      w(92)=2*w(92) + w(101)
      w(92)=w(92)*w(46)
      w(93)= - 469 - 80*z
      w(93)=z*w(93)
      w(93)= - 221 + w(93)
      w(92)=1.D0/3.D0*w(93) + w(92)
      w(92)=w(3)*w(92)
      w(93)=16*w(95)
      w(56)=4*w(85) + w(93) + w(71) + 8*w(56) + w(80) + 2*w(91) + w(92)
      w(56)=z2*w(56)
      w(71)=1 - 12*z
      w(71)=z*w(71)
      w(71)=12*w(101) - 1 + w(71) - w(76) + w(117)
      w(71)=w(71)*w(67)
      w(80)=17 + 44.D0/3.D0*z
      w(80)=w(80)*w(123)
      w(80)= - 2 + w(80)
      w(80)=2*w(80) - 11.D0/3.D0*w(101)
      w(80)=w(80)*w(46)
      w(85)=w(52)*z2
      w(91)= - 809.D0/9.D0 + w(57)
      w(91)=z*w(91)
      w(71)= - 20.D0/3.D0*w(85) + 44.D0/3.D0*w(66) + w(71) - 88.D0/3.D0
     & *w(100) + w(80) - 145.D0/9.D0 + w(91)
      w(71)=z3*w(71)
      w(80)=19 + 352.D0/3.D0*z
      w(80)=w(80)*w(104)
      w(71)=w(71) + 8*w(133) - 233.D0/9.D0*w(1) - 65 + w(80)
      w(80)= - 7 - w(58)
      w(80)=z*w(80)
      w(80)= - w(84) + w(63) - 3 + w(80)
      w(80)=w(80)*w(59)
      w(60)=w(60)*w(97)
      w(60)=w(60) - w(101)
      w(60)=w(60)*w(46)
      w(60)= - 8*w(100) + 23.D0/3.D0*w(45) + w(60)
      w(91)=1.D0/3.D0*LL
      w(60)=w(60)*w(91)
      w(92)=32*z
      w(94)=49 + w(92)
      w(94)=w(94)*w(123)
      w(95)=2.D0/3.D0*w(101)
      w(96)= - 3 + w(58)
      w(96)=z*w(96)
      w(96)= - w(95) - 1 + w(96)
      w(96)=w(3)*w(96)
      w(92)=1 + w(92)
      w(92)=z*w(92)
      w(92)= - 31 + w(92)
      w(92)=1.D0/3.D0*w(92) + w(96)
      w(92)=w(92)*w(46)
      w(60)=w(60) + w(80) + w(92) - w(48) - 23 + w(94)
      w(60)=w(60)*w(55)
      w(80)= - 664 + 121.D0/3.D0*z
      w(80)=w(80)*w(97)
      w(80)= - 211 + w(80)
      w(92)= - 1 - w(50)
      w(92)=z*w(92)
      w(92)= - 1 + w(92)
      w(92)=16.D0/3.D0*w(92) + w(81)
      w(92)=w(3)*w(92)
      w(94)= - 37 + 136.D0/3.D0*z
      w(94)=z*w(94)
      w(92)=w(92) + 71 + w(94)
      w(92)=w(3)*w(92)
      w(80)=4.D0/9.D0*w(80) + w(92)
      w(80)=w(3)*w(80)
      w(92)=104.D0/3.D0*w(1)
      w(94)=101 + w(109)
      w(94)=z*w(94)
      w(94)=w(92) + 121 + w(94)
      w(96)= - 1 - w(127)
      w(96)=w(96)*w(78)
      w(94)=1.D0/3.D0*w(94) + w(96)
      w(94)=w(94)*w(59)
      w(96)= - 1579.D0/9.D0 + 356*z
      w(96)=z*w(96)
      w(90)= - w(90) - 1445.D0/9.D0 + w(96)
      w(96)=w(24)*w(52)
      w(60)=w(60) + w(94) - 48*w(96) + 2.D0/3.D0*w(90) + w(80)
      w(60)=LL*w(60)
      w(80)=1 + 32.D0/3.D0*z
      w(80)=z*w(80)
      w(65)=w(120) + w(65) - w(73) - 9 + w(80)
      w(65)=LL*w(65)
      w(80)=w(64)*w(59)
      w(90)= - 89 + 140.D0/3.D0*z
      w(90)=z*w(90)
      w(90)= - w(92) + 77 + w(90)
      w(90)=w(90)*w(86)
      w(92)= - 52 + 229.D0/27.D0*z
      w(92)=z*w(92)
      w(92)=239.D0/27.D0*w(1) + 104.D0/3.D0 + w(92)
      w(65)=w(65) + w(80) + 2*w(92) + w(90)
      w(65)=w(65)*w(55)
      w(90)=20.D0/3.D0*z
      w(92)= - 7 - w(90)
      w(92)=z*w(92)
      w(92)= - 2.D0/3.D0*w(125) - w(63) + 15 + w(92)
      w(92)=w(3)*w(92)
      w(94)=22*w(1)
      w(96)=59 - w(83)
      w(96)=z*w(96)
      w(96)= - 109 + w(96)
      w(96)=1.D0/3.D0*w(96) + w(94)
      w(92)=2*w(96) + w(92)
      w(92)=w(3)*w(92)
      w(45)= - 1.D0/3.D0*w(45) + w(125)
      w(45)=w(45)*w(126)
      w(96)=80.D0/3.D0*z
      w(111)=87 - w(96)
      w(111)=z*w(111)
      w(45)=w(65) + w(45) + w(92) - 3*w(1) - 172.D0/3.D0 + w(111)
      w(65)=w(136)*w(103)
      w(65)= - w(48) - 17 + w(65)
      w(65)=1.D0/3.D0*w(65) - w(69)
      w(65)=LL*w(65)
      w(57)=37.D0/3.D0 - w(57)
      w(57)=z*w(57)
      w(57)=w(65) + w(80) + 17.D0/3.D0*w(1) - 10 + w(57)
      w(65)=7.D0/3.D0 - w(97)
      w(65)=z*w(65)
      w(65)= - w(107) - w(118) + 5.D0/3.D0 + w(65)
      w(65)=4*w(65) - w(66)
      w(65)=w(65)*w(132)
      w(57)=2*w(57) + w(65)
      w(57)=w(14)*w(57)
      w(45)=2*w(45) + w(57)
      w(45)=w(14)*w(45)
      w(57)=w(81) + 5
      w(65)= - 11 + w(58)
      w(65)=z*w(65)
      w(65)=w(63) + w(65) - w(57)
      w(65)=w(65)*w(46)
      w(57)=w(53) + 8.D0/3.D0*w(134) + w(57)
      w(57)=w(57)*w(67)
      w(69)=8*w(66)
      w(80)=19 - w(50)
      w(80)=w(80)*w(102)
      w(80)= - w(48) + 43 + w(80)
      w(57)= - w(69) + w(57) + 1.D0/3.D0*w(80) + w(65)
      w(57)=w(17)*w(57)
      w(65)= - 9 - w(50)
      w(65)=z*w(65)
      w(65)=w(122) + w(73) - 1 + w(65)
      w(65)=w(65)*w(55)
      w(73)=23.D0/3.D0 - w(50)
      w(73)=z*w(73)
      w(65)=10*w(85) - w(69) + w(65) - 8*w(125) + 22.D0/3.D0 + w(73)
      w(65)=w(65)*w(115)
      w(69)=191 + 160*z
      w(69)=w(69)*w(123)
      w(73)=1.D0/5.D0*w(101) - 1 - w(88)
      w(73)=w(3)*w(73)
      w(70)= - 149 - w(70)
      w(70)=z*w(70)
      w(70)= - 115 + w(70)
      w(70)=1.D0/3.D0*w(70) + w(73)
      w(70)=w(70)*w(86)
      w(69)=w(70) - 39 + w(69)
      w(69)=w(3)*w(69)
      w(70)=9 - w(96)
      w(70)=w(70)*w(102)
      w(70)= - 74 + w(70)
      w(69)=2*w(70) + w(69)
      w(69)=w(3)*w(69)
      w(63)=w(63) - 11
      w(70)=7 - w(90)
      w(70)=z*w(70)
      w(70)=10*w(101) + w(70) - w(63)
      w(70)=w(70)*w(108)
      w(73)=9 + w(58)
      w(73)=z*w(73)
      w(73)=w(95) + w(73) - w(89)
      w(73)=w(3)*w(73)
      w(80)=z*w(82)
      w(80)=w(118) - 13 + w(80)
      w(73)=2.D0/3.D0*w(80) + w(73)
      w(73)=w(3)*w(73)
      w(80)= - 25 - w(97)
      w(80)=w(80)*w(50)
      w(80)=85 + w(80)
      w(73)=w(73) + 1.D0/3.D0*w(80) - w(94)
      w(73)=w(73)*w(59)
      w(80)=w(105) - w(131)
      w(80)= - w(53) + 2*w(80) - w(101)
      w(80)=w(26)*w(80)
      w(62)=w(76) + 5*w(62) - w(124)
      w(62)=w(27)*w(62)
      w(81)= - 96*w(37) + 64*w(41) - 32*w(43)
      w(81)=w(52)*w(81)
      w(82)=4*w(53)
      w(85)=7 + 16.D0/3.D0*z
      w(85)=z*w(85)
      w(85)= - w(82) - 16.D0/3.D0*w(1) - 5 + w(85)
      w(85)=w(25)*w(85)
      w(45)=w(65) + w(56) + 16*w(85) + 8*w(62) + 4*w(57) + w(45) + 32*
     & w(80) + w(60) + w(73) + w(70) + w(69) + w(81) + 2*w(71)
      w(45)=w(49)*w(45)
      w(56)=w(49)*w(52)
      w(57)=32*w(39) + 80*z5
      w(57)=w(56)*w(57)
      w(45)=w(57) + w(45)
      w(45)=CF*w(45)
      w(57)=w(64) - w(84)
      w(60)=w(57)*w(91)
      w(62)=w(50) - 7
      w(62)=w(62)*z
      w(62)=w(62) - 13
      w(65)=1.D0/3.D0*w(62)
      w(69)= - w(65) + w(101)
      w(69)=w(3)*w(69)
      w(70)=w(97) - 7
      w(70)=w(70)*z
      w(70)=w(70) + 5
      w(71)=2*w(100)
      w(60)=w(60) - w(71) + w(69) + w(70)
      w(60)=LL*w(60)
      w(69)= - 1 - w(97)
      w(69)=z*w(69)
      w(69)= - w(48) + 2 + w(69)
      w(69)=1.D0/3.D0*w(69) + w(84)
      w(69)=w(69)*w(87)
      w(62)=w(62) - w(84)
      w(73)=w(3)*w(62)
      w(80)=2 + z
      w(80)=z*w(80)
      w(80)= - 58 + w(80)
      w(73)=2.D0/3.D0*w(80) + w(73)
      w(73)=w(3)*w(73)
      w(75)=23 + w(75)
      w(75)=z*w(75)
      w(75)= - 34.D0/3.D0*w(1) - 20 + w(75)
      w(73)=4.D0/3.D0*w(75) + w(73)
      w(60)=w(60) + 1.D0/3.D0*w(73) + w(69)
      w(60)=LL*w(60)
      w(69)= - 13 + 74.D0/9.D0*z
      w(69)=w(69)*z
      w(69)=w(69) + 7 - 20.D0/9.D0*w(1)
      w(73)= - w(128) + w(69)
      w(73)=w(73)*w(91)
      w(75)=w(64)*w(86)
      w(80)= - 1 + 28.D0/27.D0*z
      w(80)=w(80)*z
      w(80)=w(80) + 1.D0/3.D0 - 10.D0/27.D0*w(1)
      w(75)= - w(75) + 4*w(80)
      w(75)=w(75)*w(3)
      w(73)=w(75) + w(73)
      w(81)=w(66)*w(91)
      w(73)=2*w(73) + w(81)
      w(73)=w(73)*w(72)
      w(81)=w(109) + 5
      w(81)=w(81)*z
      w(81)=w(81) + 1
      w(85)= - w(101) + 2*w(81)
      w(85)=w(85)*w(3)
      w(87)=11 + 56.D0/3.D0*z
      w(87)=w(87)*z
      w(87)=w(87) + 7
      w(85)=w(85) - 4*w(87)
      w(85)=w(85)*w(3)
      w(88)= - 19 + 800.D0/3.D0*z
      w(88)=z*w(88)
      w(88)=37 + w(88)
      w(88)=4.D0/3.D0*w(88) + w(85)
      w(88)=w(3)*w(88)
      w(89)=151 - 1156.D0/9.D0*z
      w(89)=z*w(89)
      w(89)=328.D0/9.D0*w(1) - 59 + w(89)
      w(88)=8.D0/3.D0*w(89) + w(88)
      w(89)=w(123) + 1
      w(89)=w(89)*w(97)
      w(89)=w(89) + w(106)
      w(89)= - w(101) + 2*w(89)
      w(90)= - w(89)*w(86)
      w(80)= - 2*w(80) + w(90)
      w(80)=w(80)*w(126)
      w(60)=w(73) + 4.D0/3.D0*w(60) + 1.D0/9.D0*w(88) + w(80)
      w(60)=nf*w(60)
      w(65)= - w(53) + w(65) - w(84)
      w(65)=w(65)*w(55)
      w(50)=w(50) + 37.D0/3.D0
      w(50)=w(50)*z
      w(50)=w(50) - w(101) + 19.D0/3.D0
      w(50)=w(50)*w(3)
      w(73)=w(71) - w(66)
      w(80)=179 - 56*z
      w(80)=z*w(80)
      w(68)=w(68) + 95 + w(80)
      w(65)=w(65) + 1.D0/9.D0*w(68) + w(50) - w(73)
      w(68)=z + 2.D0/3.D0
      w(68)=w(68)*w(97)
      w(68)=w(68) - 5.D0/3.D0
      w(80)= - w(110) + w(68)
      w(80)=w(80)*w(67)
      w(68)= - w(101) + w(68)
      w(68)=w(68)*w(46)
      w(68)=w(80) + w(68) + w(73) - w(69)
      w(68)=nf*w(68)
      w(69)=14.D0/5.D0 - nf
      w(69)=w(69)*w(79)
      w(65)=w(69) + 2*w(65) + w(68)
      w(65)=z2*w(65)
      w(57)=w(57)*LL
      w(68)= - 5 + 94.D0/9.D0*z
      w(68)=z*w(68)
      w(50)=4.D0/3.D0*w(57) + w(71) - w(50) - 40.D0/9.D0*w(1) - 1 + 
     & w(68)
      w(50)=LL*w(50)
      w(62)=w(62)*w(86)
      w(68)= - 65.D0/9.D0 - w(97)
      w(68)=z*w(68)
      w(68)= - 101.D0/9.D0 + w(68)
      w(62)=2*w(68) + w(62)
      w(62)=w(3)*w(62)
      w(68)=5*w(52) - w(118)
      w(68)=1.D0/3.D0*w(68) + w(101)
      w(68)=w(68)*w(59)
      w(69)=13 + 40.D0/3.D0*z
      w(69)=z*w(69)
      w(69)= - 31.D0/3.D0*w(1) - 16 + w(69)
      w(50)=w(50) + w(68) + 8.D0/9.D0*w(69) + w(62)
      w(50)=w(50)*w(67)
      w(62)=41.D0/3.D0 + w(129)
      w(62)=z*w(62)
      w(62)=607.D0/3.D0 + 17*w(62)
      w(62)=4.D0/3.D0*w(62) + w(85)
      w(62)=w(3)*w(62)
      w(68)=46 - 1781.D0/27.D0*z
      w(68)=z*w(68)
      w(68)=656.D0/27.D0*w(1) - 13.D0/3.D0 + w(68)
      w(62)=8*w(68) + w(62)
      w(68)= - w(3)*w(89)
      w(69)= - 13.D0/3.D0 - 25*z
      w(69)=z*w(69)
      w(69)=20.D0/3.D0*w(1) - 61.D0/3.D0 + w(69)
      w(68)=1.D0/3.D0*w(69) + w(68)
      w(68)=w(68)*w(126)
      w(50)=w(50) + 1.D0/3.D0*w(62) + w(68)
      w(62)= - 1 - 38.D0/9.D0*z
      w(62)=w(62)*w(123)
      w(62)=1.D0/9.D0*w(66) + w(120) + 20.D0/27.D0*w(1) + 1 + w(62)
      w(62)=w(14)*w(62)
      w(68)= - w(125) + w(70)
      w(68)=2*w(68) - w(107)
      w(68)=w(68)*w(91)
      w(69)=1 + 55.D0/9.D0*z
      w(69)=z*w(69)
      w(69)= - 28.D0/9.D0*w(1) - 4 + w(69)
      w(68)=w(68) + 4.D0/9.D0*w(69) + w(75)
      w(62)=2*w(68) + w(62)
      w(62)=w(62)*w(72)
      w(68)= - 43 - 62.D0/3.D0*z
      w(68)=w(68)*w(123)
      w(68)=w(82) + 16.D0/3.D0*w(101) + 8.D0/9.D0*w(1) - 4 + w(68)
      w(69)= - 61 - 100.D0/3.D0*z
      w(69)=z*w(69)
      w(63)=26*w(101) + w(69) + w(63)
      w(63)=1.D0/9.D0*w(63) + w(76)
      w(63)=nf*w(63)
      w(63)=2.D0/3.D0*w(68) + w(63)
      w(68)=4*z3
      w(63)=w(63)*w(68)
      w(69)=11 + w(97)
      w(69)=z*w(69)
      w(69)=8 + w(69)
      w(53)= - nf*w(53)
      w(53)=w(53) + 1.D0/3.D0*w(69) - w(76)
      w(53)=w(18)*w(53)
      w(69)=w(27)*w(52)
      w(53)=w(69) - w(53)
      w(58)=w(58) + 3
      w(58)=w(58)*z
      w(54)=w(95) - w(58) - 4.D0/9.D0*w(1) + w(54) - 1
      w(58)= - nf - 1
      w(54)=w(99)*w(54)*w(58)
      w(50)=w(63) + 2.D0/3.D0*w(65) + w(54) + w(60) + 1.D0/3.D0*w(50)
     &  + w(62) - 8.D0/3.D0*w(53)
      w(49)=TF*w(49)*w(50)
      w(50)=w(84) + w(64)
      w(50)=w(50)*w(59)
      w(53)=w(101) - w(81)
      w(53)=w(3)*w(53)
      w(54)= - 3 + 28.D0/9.D0*z
      w(54)=z*w(54)
      w(54)= - 10.D0/9.D0*w(1) + 1 + w(54)
      w(53)=2*w(54) + w(53)
      w(53)=2*w(53) + w(57)
      w(53)=w(53)*w(55)
      w(54)= - w(95) + w(81)
      w(54)=w(3)*w(54)
      w(54)= - 4.D0/3.D0*w(87) + w(54)
      w(54)=w(3)*w(54)
      w(57)= - w(78)*w(66)
      w(58)= - 31 + 400.D0/9.D0*z
      w(58)=z*w(58)
      w(58)= - 112.D0/9.D0*w(1) - 1 + w(58)
      w(50)= - w(93) + w(57) + w(53) + w(50) + 2.D0/3.D0*w(58) + w(54)
      w(53)=as**2
      w(50)=w(50)*w(53)
      w(54)=193 + 128*z
      w(54)=w(54)*w(97)
      w(48)= - w(48) - 215 + w(135)
      w(57)=13 - 14*z
      w(57)=w(3)*w(57)
      w(48)=1.D0/3.D0*w(48) + w(57)
      w(48)=w(48)*w(46)
      w(57)=37 + 43*z
      w(57)=z2*w(57)
      w(48)=w(57) + 20*w(100) + w(48) + 76.D0/3.D0*w(1) + 229 + w(54)
      w(54)= - 9 - w(83)
      w(46)=w(54)*w(46)
      w(54)= - w(98)*w(67)
      w(57)= - 1 + w(104)
      w(57)=z*w(57)
      w(46)=w(54) + w(46) - 16*w(1) - 11 + 26*w(57)
      w(46)=w(46)*w(55)
      w(54)= - w(74)*w(77)
      w(46)=w(54) - 10.D0/3.D0*w(66) + w(46) - 12*w(61) + 1.D0/3.D0*
     & w(48)
      w(46)=w(46)*w(51)
      w(48)=w(52)*w(53)
      w(46)=4*w(48) + w(46)
      w(46)=w(46)*w(68)
      w(48)=w(56)*CA
      w(51)=CF*w(56)
      w(51)= - w(48) + w(51)
      w(51)=w(44)*w(51)
      w(52)= - 48*w(38) - 96*w(39)
      w(48)=w(48)*w(52)
      w(45)=16*w(51) + 4*w(49) + w(45) + w(46) + w(50) + w(47) + w(48)
*
      PS = 2*TF*CF*w(45)
*
      RETURN
      END      
      REAL*8 FUNCTION QG(z,nf,as,LL)
*     ------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of the OME AQg without aQg3 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(346)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
      w(1)=z**(-1)
      w(2)=z**(-1)
      w(3)=Hr1(-1)
      w(4)=z**(-1)
      w(5)=Hr1(-1)
      w(6)=z**(-1)
      w(7)=Hr1(-1)
      w(8)=Hr1(0)
      w(9)=Hr2(-1,-1)
      w(10)=Hr2(-1,-1)
      w(11)=Hr3(-1,-1,-1)
      w(12)=Hr4(-1,-1,-1,-1)
      w(13)=Hr5(-1,-1,-1,0,-1)
      w(14)=Hr5(-1,-1,-1,0,1)
      w(15)=Hr5(-1,-1,-1,0,1)
      w(16)=Hr4(-1,-1,0,-1)
      w(17)=Hr5(-1,-1,0,-1,-1)
      w(18)=Hr5(-1,-1,0,-1,-1)
      w(19)=Hr5(-1,-1,0,0,1)
      w(20)=Hr4(-1,-1,0,1)
      w(21)=Hr3(-1,0,-1)
      w(22)=Hr4(-1,0,-1,-1)
      w(23)=Hr5(-1,0,-1,-1,-1)
      w(24)=Hr5(-1,0,-1,0,1)
      w(25)=Hr4(-1,0,0,-1)
      w(26)=Hr5(-1,0,0,0,-1)
      w(27)=Hr5(-1,0,0,0,1)
      w(28)=Hr4(-1,0,0,1)
      w(29)=Hr5(-1,0,0,1,1)
      w(30)=Hr3(-1,0,1)
      w(31)=Hr4(-1,0,1,1)
      w(32)=Hr1(1)
      w(33)=Hr2(1,1)
      w(34)=Hr3(1,0,1)
      w(35)=Hr3(1,0,-1)
      w(36)=Hr4(1,0,0,1)
      w(37)=Hr4(1,0,1,1)
      w(38)=Hr4(1,1,0,1)
      w(39)=Hr3(1,1,1)
      w(40)=Hr4(1,1,1,1)
      w(41)=Hr2(0,-1)
      w(42)=Hr3(0,-1,-1)
      w(43)=Hr4(0,-1,-1,-1)
      w(44)=Hr5(0,-1,-1,-1,-1)
      w(45)=Hr5(0,-1,-1,0,-1)
      w(46)=Hr5(0,-1,-1,0,1)
      w(47)=Hr4(0,-1,0,-1)
      w(48)=Hr5(0,-1,0,-1,-1)
      w(49)=Hr5(0,-1,0,0,-1)
      w(50)=Hr5(0,-1,0,0,-1)
      w(51)=Hr5(0,-1,0,0,1)
      w(52)=Hr4(0,-1,0,1)
      w(53)=Hr5(0,-1,0,1,1)
      w(54)=Hr2(0,0)
      w(55)=Hr3(0,0,-1)
      w(56)=Hr4(0,0,-1,-1)
      w(57)=Hr5(0,0,-1,-1,-1)
      w(58)=Hr5(0,0,-1,0,-1)
      w(59)=Hr5(0,0,-1,0,-1)
      w(60)=Hr5(0,0,-1,0,1)
      w(61)=Hr5(0,0,-1,0,1)
      w(62)=Hr3(0,0,0)
      w(63)=Hr4(0,0,0,-1)
      w(64)=Hr5(0,0,0,-1,-1)
      w(65)=Hr4(0,0,0,0)
      w(66)=Hr5(0,0,0,0,-1)
      w(67)=Hr5(0,0,0,0,0)
      w(68)=Hr5(0,0,0,0,0)
      w(69)=Hr5(0,0,0,0,1)
      w(70)=Hr5(0,0,0,0,1)
      w(71)=Hr4(0,0,0,1)
      w(72)=Hr5(0,0,0,1,1)
      w(73)=Hr3(0,0,1)
      w(74)=Hr5(0,0,1,0,1)
      w(75)=Hr4(0,0,1,1)
      w(76)=Hr5(0,0,1,1,1)
      w(77)=Hr2(0,1)
      w(78)=Hr4(0,1,0,-1)
      w(79)=Hr5(0,1,0,0,1)
      w(80)=Hr5(0,1,0,0,1)
      w(81)=Hr4(0,1,0,1)
      w(82)=Hr5(0,1,0,1,1)
      w(83)=Hr3(0,1,1)
      w(84)=Hr5(0,1,1,0,1)
      w(85)=Hr4(0,1,1,1)
      w(86)=Hr5(0,1,1,1,1)
      w(87)=Hr4(1,0,0,-1)
      w(88)=Hr5(1,0,0,0,1)
      w(89)=Hr5(1,0,0,1,1)
      w(90)=Hr5(1,0,1,0,1)
      w(91)=Hr5(1,0,1,1,1)
      w(92)=Hr5(1,1,0,0,1)
      w(93)=Hr5(1,1,0,1,1)
      w(94)=Hr5(1,1,1,0,1)
      w(95)=Hr5(1,1,1,1,1)
      w(96)=z + 1.D0
      w(97)=2.D0*z
      w(98)=w(96)*w(97)
      w(99)=w(98) + 1.D0
      w(100)=3.D0*w(99)
      w(101)=w(99)*w(8)
      w(102)= - w(100) + w(101)
      w(103)=8.D0*w(8)
      w(104)=w(103)*w(7)
      w(102)=w(102)*w(104)
      w(105)=20.D0*z
      w(106)=w(97) - 1.D0
      w(107)=w(106)*w(105)
      w(108)= - 47.D0 + 146.D0*z
      w(108)=w(108)*w(97)
      w(108)=11.D0 + w(108)
      w(109)=2.D0/3.D0*w(8)
      w(108)=w(108)*w(109)
      w(107)=w(108) - 31.D0 + w(107)
      w(108)=2.D0/3.D0*LL
      w(107)=w(107)*w(108)
      w(110)=4.D0*z
      w(111)= - 397.D0 + 586.D0*z
      w(111)=w(111)*w(110)
      w(111)= - 545.D0 + w(111)
      w(112)=26.D0*z
      w(113)=88.D0*z
      w(114)= - 109.D0 + w(113)
      w(114)=w(114)*w(112)
      w(114)= - 587.D0 + w(114)
      w(115)=2.D0/9.D0*w(8)
      w(114)=w(114)*w(115)
      w(116)=32.D0*z3
      w(117)=w(110) - 1.D0
      w(118)= - z*w(117)
      w(118)= - 2.D0 + w(118)
      w(118)=w(118)*w(116)
      w(119)=32.D0*w(99)
      w(120)=w(30)*w(119)
      w(121)=z - 8.D0
      w(121)=w(121)*z
      w(122)=1.D0 - w(121)
      w(122)=w(73)*w(122)
      w(123)=w(97) + 1.D0
      w(124)=w(83)*w(123)
      w(102)=w(107) - 48.D0*w(124) + w(120) + w(102) + w(118) + w(114)
     &  + 1.D0/3.D0*w(111) + 32.D0*w(122)
      w(102)=LL*w(102)
      w(107)=7.D0*z
      w(111)=1.D0 + w(107)
      w(111)=w(111)*w(97)
      w(111)= - 3.D0 + w(111)
      w(111)=w(56)*w(111)
      w(114)=w(110) + 1.D0
      w(118)=w(114)*w(97)
      w(120)= - 1.D0 + w(118)
      w(120)=w(47)*w(120)
      w(122)=2.D0*w(8)
      w(124)=w(99)*w(122)
      w(125)=z + 2.D0
      w(126)=w(125)*z
      w(127)=w(126) + 1.D0
      w(128)=7.D0*w(127) + w(124)
      w(128)=w(30)*w(128)
      w(129)=5.D0*z
      w(130)=w(129) + 2.D0
      w(131)=w(65)*w(130)
      w(111)=w(111) + w(120) + w(128) + w(131)
      w(120)=z - 1.D0
      w(128)=w(120)*w(97)
      w(131)=w(128) + 1.D0
      w(132)=w(131)*w(122)
      w(133)=10.D0*z
      w(134)= - 3.D0 + 8.D0/3.D0*z
      w(134)=w(134)*w(133)
      w(135)=8.D0/3.D0*w(6)
      w(134)=w(132) - w(135) + 3.D0 + w(134)
      w(134)=w(34)*w(134)
      w(136)=w(110)*w(96)
      w(137)=26.D0/3.D0*w(101) - 9.D0 - w(136)
      w(137)=w(8)*w(137)
      w(138)=23.D0*z
      w(139)=w(138) - 12.D0
      w(139)=w(139)*z
      w(139)=w(139) - 35.D0
      w(137)=4.D0*w(139) + w(137)
      w(137)=w(8)*w(137)
      w(140)=96.D0*w(99)
      w(141)= - z3*w(140)
      w(137)=w(137) + w(141)
      w(141)=4.D0*w(7)
      w(137)=w(137)*w(141)
      w(142)=124.D0*z
      w(143)= - 277.D0 + w(142)
      w(143)=w(143)*w(97)
      w(144)=32.D0*w(6)
      w(143)= - w(144) + 25.D0 + w(143)
      w(145)=22.D0*z
      w(146)=w(145) + 13.D0
      w(147)=w(146)*w(122)
      w(143)=1.D0/3.D0*w(143) + w(147)
      w(147)=8.D0*w(83)
      w(143)=w(143)*w(147)
      w(148)=w(131)*w(37)
      w(149)=192.D0*w(131)
      w(150)=w(38) + w(87)
      w(151)=w(149)*w(150)
      w(152)= - 256.D0*w(20) - 320.D0*w(22) + 64.D0*w(31)
      w(152)=w(99)*w(152)
      w(153)=w(101)*w(11)
      w(154)=17.D0*z
      w(155)=w(154) + 4.D0
      w(156)=w(155)*w(97)
      w(156)= - 25.D0 + w(156)
      w(156)=w(73)*w(156)
      w(157)=3.D0*z
      w(158)=w(157) - 5.D0
      w(159)=w(158)*w(97)
      w(160)= - 1.D0 + w(159)
      w(161)=32.D0*w(81)
      w(160)=w(160)*w(161)
      w(162)=13.D0*z
      w(163)= - 7.D0 + w(162)
      w(163)=w(163)*w(97)
      w(163)=15.D0 + w(163)
      w(164)=4.D0*w(73)
      w(163)=w(163)*w(164)
      w(165)= - 10579.D0 + 4322.D0*z
      w(165)=z*w(165)
      w(165)=200.D0 + w(165)
      w(163)=1.D0/27.D0*w(165) + w(163)
      w(165)=4.D0*w(8)
      w(163)=w(163)*w(165)
      w(166)=w(28)*w(140)
      w(167)= - 1.D0 - w(154)
      w(167)=w(167)*w(97)
      w(167)= - 25.D0 + w(167)
      w(168)=16.D0*w(71)
      w(167)=w(167)*w(168)
      w(169)=8.D0 + 45.D0*z
      w(169)=w(169)*w(122)
      w(170)=153.D0 - 445.D0/3.D0*z
      w(170)=z*w(170)
      w(169)=w(169) + 5.D0 + w(170)
      w(170)=16.D0*z3
      w(169)=w(169)*w(170)
      w(171)=w(6)*w(73)
      w(172)= - 48955.D0/3.D0 + 20458.D0*z
      w(172)=z*w(172)
      w(172)= - 30697.D0/6.D0 + w(172)
      w(173)=z**2
      w(174)=1.D0 + 2.D0*w(173)
      w(174)=w(52)*w(174)
      w(102)=w(102) + w(143) + 288.D0*w(148) + w(137) + 320.D0*w(153)
     &  + 16.D0*w(134) + w(169) + w(167) + w(166) + w(163) + 128.D0*
     & w(174) + w(160) + 128.D0/3.D0*w(171) + 32.D0/3.D0*w(156) + 1.D0/
     & 9.D0*w(172) + w(152) + w(151) + 32.D0*w(111)
      w(102)=LL*w(102)
      w(111)=w(131)*LL
      w(134)=24.D0*w(111)
      w(137)=w(131)*w(103)
      w(143)=4.D0/3.D0*w(6)
      w(151)=58.D0*z
      w(152)=137.D0/3.D0 - w(151)
      w(152)=z*w(152)
      w(152)=w(134) + w(137) + w(143) + 32.D0/3.D0 + w(152)
      w(152)=w(40)*w(152)
      w(156)=6.D0*w(101)
      w(160)=2.D0*LL
      w(163)=w(99)*w(160)
      w(166)=4.D0*w(173)
      w(167)=w(166) - 1.D0
      w(169)=w(163) - w(156) - w(167)
      w(169)=w(25)*w(169)
      w(172)=w(173)*w(165)
      w(174)=37.D0*z
      w(175)=164.D0/3.D0 - w(174)
      w(175)=z*w(175)
      w(175)=18.D0*w(111) + w(172) - 7.D0/3.D0 + w(175)
      w(175)=w(85)*w(175)
      w(176)=w(157) + 1.D0
      w(176)=w(176)*w(97)
      w(177)=w(176) + 1.D0
      w(178)= - w(177)*w(122)
      w(178)=w(178) - w(167)
      w(178)=w(56)*w(178)
      w(179)= - 2.D0 - w(157)
      w(179)=z*w(179)
      w(179)= - 1.D0 + w(179)
      w(179)=w(179)*w(165)
      w(167)=w(179) - w(167)
      w(167)=w(47)*w(167)
      w(179)=6.D0*z
      w(180)=29.D0 - w(179)
      w(180)=w(180)*w(97)
      w(180)=19.D0 + w(180)
      w(180)=w(80)*w(180)
      w(152)=w(178) + w(167) + w(152) + w(169) + w(175) + w(180)
      w(167)=w(131)*w(8)
      w(169)=6.D0*w(167)
      w(175)=12.D0*w(111)
      w(178)= - 166.D0 + 193.D0*z
      w(178)=z*w(178)
      w(178)= - 16.D0 + w(178)
      w(178)=w(175) + 1.D0/3.D0*w(178) - w(169)
      w(178)=w(178)*w(160)
      w(180)=w(157) + 4.D0
      w(181)= - z*w(180)
      w(181)=4.D0 + w(181)
      w(181)=w(181)*w(122)
      w(182)=4.D0*w(6)
      w(183)=166.D0 - 125.D0*z
      w(183)=z*w(183)
      w(178)=w(178) + w(181) - w(182) - 13.D0 + w(183)
      w(178)=w(39)*w(178)
      w(181)=w(110) + 5.D0
      w(183)=w(181)*w(110)
      w(183)=w(183) + 1.D0
      w(184)=w(118) + 1.D0
      w(185)=w(184)*w(122)
      w(186)=4.D0*LL
      w(187)=2.D0 + 29.D0*z
      w(187)=z*w(187)
      w(187)=16.D0 + w(187)
      w(187)=w(187)*w(186)
      w(185)=w(187) + w(185) - w(183)
      w(185)=w(63)*w(185)
      w(187)=191.D0*z
      w(188)=181.D0 - w(187)
      w(188)=w(188)*w(110)
      w(189)=64.D0*w(6)
      w(188)=w(189) - 11.D0 + w(188)
      w(188)=20.D0*w(111) + 1.D0/3.D0*w(188) - 16.D0*w(167)
      w(188)=w(36)*w(188)
      w(190)=10.D0*w(167)
      w(187)= - 184.D0 + w(187)
      w(187)=w(187)*w(97)
      w(144)= - w(144) + 7.D0 + w(187)
      w(144)=1.D0/3.D0*w(144) + w(190)
      w(144)=w(8)*w(144)
      w(187)= - 14.D0/3.D0 + w(129)
      w(187)=z*w(187)
      w(187)= - 1.D0/3.D0 + w(187)
      w(144)=2.D0*w(187) + w(144)
      w(144)=w(34)*w(144)
      w(187)=w(133) - 7.D0
      w(187)=w(187)*w(97)
      w(191)=11.D0 - w(187)
      w(191)=z5*w(191)
      w(192)=28.D0*w(167) + w(117)
      w(192)=w(37)*w(192)
      w(144)=w(191) + w(144) + w(192) + w(178) + w(185) + w(188)
      w(178)=8.D0*w(73)
      w(185)=28.D0/3.D0 + 55.D0*z
      w(185)=z*w(185)
      w(185)= - 35.D0/3.D0 + w(185)
      w(185)=w(185)*w(178)
      w(188)=8.D0*w(81)
      w(146)= - w(146)*w(188)
      w(191)=16.D0*w(52)
      w(192)= - w(123)*w(191)
      w(193)=16.D0*w(8)
      w(194)=w(193)*w(73)
      w(195)= - 3.D0 + w(121)
      w(195)=w(195)*w(194)
      w(196)=2497.D0 - 131320.D0/27.D0*z
      w(196)=z*w(196)
      w(196)= - 1673.D0/3.D0 + w(196)
      w(146)=w(195) + w(192) + w(146) + 32.D0*w(171) + 1.D0/3.D0*w(196)
     &  + w(185)
      w(146)=w(146)*w(122)
      w(185)=w(143) + 1.D0
      w(192)= - 16.D0 + 55.D0/3.D0*z
      w(192)=z*w(192)
      w(192)=w(192) - w(185)
      w(192)=2.D0*w(192) + w(167)
      w(192)=w(8)*w(192)
      w(195)=73.D0*z
      w(196)= - 71.D0 + w(195)
      w(196)=w(196)*w(97)
      w(197)=16.D0*w(6)
      w(196)= - w(197) + 5.D0 + w(196)
      w(196)=1.D0/3.D0*w(196) + w(137)
      w(198)=1.D0/3.D0*LL
      w(196)=w(196)*w(198)
      w(199)=28.D0*w(6)
      w(200)= - 937.D0 + 842.D0*z
      w(200)=z*w(200)
      w(200)= - w(199) + 83.D0 + w(200)
      w(192)=w(196) + 1.D0/9.D0*w(200) + w(192)
      w(192)=w(192)*w(160)
      w(196)=w(131)*w(109)
      w(200)=137.D0/3.D0*z
      w(201)=60.D0 - w(200)
      w(201)=w(201)*w(97)
      w(202)=16.D0/3.D0*w(6)
      w(201)=w(196) + w(202) - 37.D0 + w(201)
      w(201)=w(8)*w(201)
      w(203)=9.D0*z
      w(204)= - 262.D0/9.D0 + w(203)
      w(204)=w(204)*w(97)
      w(201)=w(201) - 208.D0/9.D0*w(6) + 415.D0/9.D0 + w(204)
      w(201)=w(8)*w(201)
      w(204)=w(131)*z3
      w(205)= - 1802.D0 + 1613.D0*z
      w(205)=w(205)*w(129)
      w(205)=1094.D0 + w(205)
      w(205)=1.D0/3.D0*w(205) + 34.D0*w(6)
      w(192)=w(192) - 168.D0*w(204) + 2.D0/9.D0*w(205) + w(201)
      w(192)=w(192)*w(160)
      w(201)=119.D0 - 137.D0*z
      w(201)=w(201)*w(97)
      w(201)=7.D0 + w(201)
      w(205)=8.D0*w(6)
      w(201)= - w(167) + 1.D0/3.D0*w(201) + w(205)
      w(201)=w(201)*w(122)
      w(206)= - 169.D0 + 310.D0/3.D0*z
      w(206)=w(206)*w(110)
      w(201)=w(201) - 184.D0/3.D0*w(6) + 323.D0 + w(206)
      w(206)=1.D0/3.D0*w(8)
      w(201)=w(201)*w(206)
      w(207)=44.D0*w(6)
      w(113)= - 457.D0/3.D0 + w(113)
      w(113)=z*w(113)
      w(113)= - 380.D0/3.D0 + w(113)
      w(113)=1.D0/3.D0*w(113) + w(207)
      w(113)=2.D0*w(113) + w(201)
      w(113)=w(8)*w(113)
      w(201)= - 727.D0 + 833.D0*z
      w(201)=w(201)*w(97)
      w(201)= - 128.D0*w(6) - 17.D0 + w(201)
      w(201)=1.D0/3.D0*w(201) + 76.D0*w(167)
      w(208)=4.D0/3.D0*z3
      w(201)=w(201)*w(208)
      w(209)=9977.D0 - 81500.D0/9.D0*z
      w(209)=z*w(209)
      w(209)=2381.D0/9.D0*w(6) - 1018.D0 + w(209)
      w(113)=w(192) + w(201) + 2.D0/9.D0*w(209) + w(113)
      w(192)=2.D0*w(32)
      w(113)=w(113)*w(192)
      w(201)=z - 5.D0
      w(209)=w(201)*w(97)
      w(210)=w(209) - 1.D0
      w(211)=w(8)*w(210)
      w(212)= - w(123)*w(160)
      w(213)=11.D0*z
      w(214)= - 26.D0/3.D0 - w(213)
      w(214)=z*w(214)
      w(211)=w(212) + w(211) + w(143) + 13.D0/3.D0 + w(214)
      w(211)=w(211)*w(186)
      w(212)=18.D0*z
      w(214)= - 173.D0/3.D0 - w(212)
      w(214)=w(214)*w(97)
      w(215)=10.D0 - w(107)
      w(215)=z*w(215)
      w(215)= - 2.D0 + w(215)
      w(215)=w(215)*w(122)
      w(214)=w(215) - w(202) + 101.D0/3.D0 + w(214)
      w(214)=w(214)*w(122)
      w(215)= - 973.D0 + 676.D0*z
      w(215)=w(215)*w(97)
      w(215)=208.D0*w(6) - 215.D0 + w(215)
      w(211)=w(211) + 1.D0/9.D0*w(215) + w(214)
      w(211)=LL*w(211)
      w(214)=16.D0*z
      w(215)= - 85.D0 - w(214)
      w(215)=w(215)*w(97)
      w(215)=25.D0 + w(215)
      w(216)=z + 4.D0
      w(217)=w(216)*w(97)
      w(218)=5.D0 + w(217)
      w(219)=4.D0/3.D0*w(8)
      w(218)=w(218)*w(219)
      w(215)=w(218) + 1.D0/3.D0*w(215) - w(205)
      w(215)=w(8)*w(215)
      w(218)=65.D0/3.D0 + z
      w(218)=z*w(218)
      w(218)=23.D0/9.D0*w(6) - 11.D0 + w(218)
      w(215)=8.D0*w(218) + w(215)
      w(215)=w(8)*w(215)
      w(218)= - w(123)*w(179)
      w(218)= - 13.D0 + w(218)
      w(220)=8.D0*z3
      w(218)=w(218)*w(220)
      w(221)=1181.D0 + 5540.D0/9.D0*z
      w(221)=z*w(221)
      w(221)=140.D0 + w(221)
      w(207)=w(211) + w(218) + w(215) + 1.D0/3.D0*w(221) - w(207)
      w(211)=4.D0*w(77)
      w(207)=w(207)*w(211)
      w(215)=w(101) - w(98)
      w(218)=w(215)*w(20)
      w(221)=w(101)*w(31)
      w(218)=w(221) - w(218)
      w(221)=w(106)*w(110)
      w(222)= - 5.D0 - w(221)
      w(222)=w(84)*w(222)
      w(130)=z*w(130)
      w(130)=1.D0 + w(130)
      w(130)=w(64)*w(130)
      w(223)=w(129) + 4.D0
      w(223)=w(223)*z
      w(224)=2.D0 + w(223)
      w(224)=w(50)*w(224)
      w(130)=w(222) + w(224) + w(218) + w(130)
      w(222)=z + 3.D0
      w(224)=w(222)*w(97)
      w(225)= - w(124) - 1.D0 - w(224)
      w(225)=w(28)*w(225)
      w(226)=w(97) + 3.D0
      w(227)=w(226)*w(110)
      w(227)=w(227) + 1.D0
      w(228)=w(99)*w(165)
      w(229)=w(227) + w(228)
      w(230)=w(99)*LL
      w(231)=10.D0*w(230)
      w(232)= - w(231) + w(229)
      w(232)=w(16)*w(232)
      w(233)=w(66)*w(177)
      w(225)=w(232) + w(225) - w(233)
      w(232)= - w(123)*w(110)
      w(232)=w(101) + 1.D0 + w(232)
      w(232)=w(232)*w(206)
      w(233)= - 7.D0 - w(129)
      w(233)=z*w(233)
      w(233)=4.D0 + w(233)
      w(232)=2.D0*w(233) + w(232)
      w(232)=w(8)*w(232)
      w(233)=w(212) + 23.D0
      w(233)=w(233)*z
      w(233)=w(233) + 13.D0
      w(233)=4.D0*w(233)
      w(232)= - w(233) + w(232)
      w(232)=w(8)*w(232)
      w(234)=w(99)*w(103)
      w(235)=10.D0 + w(107)
      w(235)=w(235)*w(110)
      w(235)=w(234) + 3.D0 + w(235)
      w(236)=2.D0*z3
      w(235)=w(235)*w(236)
      w(232)=w(232) + w(235)
      w(232)=w(232)*w(141)
      w(235)=w(129) - 13.D0
      w(237)=w(235)*w(110)
      w(238)=15.D0*z
      w(239)=1.D0 - w(238)
      w(239)=w(239)*w(97)
      w(239)= - 11.D0 + w(239)
      w(239)=w(239)*w(122)
      w(237)=w(163) + w(239) + 95.D0 + w(237)
      w(237)=LL*w(237)
      w(172)= - w(172) + w(227)
      w(172)=w(8)*w(172)
      w(239)=19.D0*z
      w(240)=5.D0 + w(239)
      w(240)=z*w(240)
      w(240)=4.D0 + w(240)
      w(172)=w(237) + 2.D0*w(240) + w(172)
      w(237)=8.D0*w(55)
      w(172)=w(172)*w(237)
      w(240)=w(97) - 7.D0
      w(241)= - w(240)*w(110)
      w(242)= - 6.D0 + z
      w(242)=z*w(242)
      w(242)= - 2.D0 + w(242)
      w(242)=w(242)*w(122)
      w(241)=w(242) - 43.D0 + w(241)
      w(241)=w(8)*w(241)
      w(100)=w(100) - w(124)
      w(100)=LL*w(100)
      w(100)=w(100) - 2.D0*w(139) + w(241)
      w(100)=w(100)*w(160)
      w(139)= - w(123)*w(109)
      w(139)=w(139) - w(114)
      w(139)=w(8)*w(139)
      w(241)=w(107) - 1.D0
      w(242)= - z*w(241)
      w(242)= - 4.D0 + w(242)
      w(139)=4.D0*w(242) + w(139)
      w(139)=w(8)*w(139)
      w(242)=w(106)*w(97)
      w(243)= - 1.D0 + w(242)
      w(244)=4.D0*z3
      w(243)=w(243)*w(244)
      w(100)=w(100) + w(243) + w(233) + w(139)
      w(139)=4.D0*w(41)
      w(100)=w(100)*w(139)
      w(233)=32.D0*w(8)
      w(243)= - w(11)*w(233)
      w(243)=w(243) + 32.D0*w(22)
      w(229)=w(229)*w(243)
      w(243)=w(53) - w(24)
      w(245)=w(243) - w(46)
      w(246)=w(23) + w(44)
      w(247)=w(246) + w(13)
      w(248)=128.D0*w(15)
      w(245)= - w(248) - 128.D0*w(29) - 384.D0*w(247) - 64.D0*w(245)
      w(245)=w(99)*w(245)
      w(247)=82.D0 - 317.D0*z
      w(247)=w(247)*w(97)
      w(247)=137.D0 + w(247)
      w(249)=z - 2.D0
      w(250)= - w(249)*w(129)
      w(250)=2.D0 + w(250)
      w(250)=w(250)*w(103)
      w(247)=w(250) + 1.D0/3.D0*w(247) - w(205)
      w(250)=8.D0*w(71)
      w(247)=w(247)*w(250)
      w(251)= - 995.D0 + 793.D0*z
      w(251)=w(251)*w(110)
      w(252)= - 751.D0 - 134.D0*z
      w(252)=w(252)*w(97)
      w(252)= - 149.D0 + w(252)
      w(252)=w(252)*w(109)
      w(251)=w(252) + 203.D0 + w(251)
      w(252)=2.D0/3.D0*z3
      w(251)=w(251)*w(252)
      w(253)=12.D0*z
      w(254)=w(253)*w(96)
      w(255)=w(124) + w(114)
      w(255)=w(8)*w(255)
      w(255)= - w(254) + w(255)
      w(256)=16.D0*w(30)
      w(255)=w(255)*w(256)
      w(257)=w(97) - 5.D0/3.D0
      w(258)=14.D0*z
      w(259)=w(257)*w(258)
      w(260)=w(135) - 5.D0/3.D0
      w(259)=w(259) - w(260)
      w(261)=20.D0*w(167)
      w(262)=w(261) + w(259)
      w(263)=16.D0*w(38)
      w(262)=w(262)*w(263)
      w(264)=7.D0/3.D0*z
      w(265)= - 15.D0 - w(264)
      w(265)=w(265)*w(110)
      w(226)=w(226)*w(97)
      w(266)=9.D0 + w(226)
      w(266)=w(8)*w(266)
      w(265)=w(266) - w(202) + 3.D0 + w(265)
      w(265)=w(8)*w(265)
      w(266)=80.D0/9.D0*w(6)
      w(267)=307.D0 - 2069.D0/9.D0*z
      w(267)=z*w(267)
      w(265)=w(265) + w(266) - 76.D0/3.D0 + w(267)
      w(265)=w(265)*w(147)
      w(267)=1.D0/3.D0*z
      w(268)=w(267) + 1.D0
      w(268)=w(268)*w(97)
      w(268)=w(268) + 1.D0
      w(269)= - w(114)*w(198)
      w(269)=w(269) - w(268)
      w(270)=8.D0*LL
      w(269)=w(269)*w(270)
      w(271)=622.D0*z
      w(272)=1883.D0/3.D0 - w(271)
      w(272)=z*w(272)
      w(272)=353.D0/3.D0 + w(272)
      w(269)=1.D0/3.D0*w(272) + w(269)
      w(269)=w(269)*w(160)
      w(257)=w(257)*w(97)
      w(257)=11.D0/3.D0 + w(257)
      w(257)=w(257)*w(220)
      w(272)=1153.D0 + 43264.D0/3.D0*z
      w(272)=z*w(272)
      w(272)= - 16.D0 + w(272)
      w(257)=w(269) + 1.D0/9.D0*w(272) + w(257)
      w(269)=2.D0*w(54)
      w(257)=w(257)*w(269)
      w(272)=8.D0*z
      w(273)= - 23.D0 - 634.D0/9.D0*z
      w(273)=w(273)*w(272)
      w(274)=1.D0 - 118.D0*z
      w(274)=w(274)*w(97)
      w(274)= - 11.D0 + w(274)
      w(274)=w(274)*w(160)
      w(273)=w(274) - 331.D0/3.D0 + w(273)
      w(274)=2.D0*w(62)
      w(273)=w(273)*w(274)
      w(275)=w(48) + w(45)
      w(276)= - 128.D0*w(57) - 64.D0*w(275)
      w(276)=w(123)*w(276)
      w(277)=w(111)*w(78)
      w(278)=w(93) - w(91)
      w(278)= - 544.D0*w(92) - 320.D0*w(95) - 416.D0*w(89) + 96.D0*
     & w(278)
      w(278)=w(131)*w(278)
      w(279)= - 8069.D0/27.D0 + 16.D0*w(27)
      w(279)=w(279)*w(97)
      w(280)=32.D0*w(27)
      w(279)=w(279) + 66073.D0/81.D0 + w(280)
      w(279)=w(279)*w(97)
      w(281)=35.D0 - 1042.D0/3.D0*z
      w(281)=w(281)*w(97)
      w(281)=269.D0 + w(281)
      w(282)=4.D0/3.D0*w(73)
      w(281)=w(281)*w(282)
      w(283)=54.D0 - 107.D0/3.D0*z
      w(283)=w(283)*w(97)
      w(284)=32.D0/3.D0*w(6)
      w(283)=w(284) + 9.D0 + w(283)
      w(283)=w(283)*w(188)
      w(285)= - w(114)*w(191)
      w(286)=w(258) - 11.D0
      w(286)=w(286)*w(97)
      w(287)=1.D0 - w(286)
      w(288)=16.D0*w(72)
      w(287)=w(287)*w(288)
      w(289)=w(179) - 7.D0
      w(290)= - w(289)*w(97)
      w(290)= - 5.D0 + w(290)
      w(290)=w(86)*w(290)
      w(291)=89.D0 - 38.D0*z
      w(291)=w(291)*w(97)
      w(291)= - 11.D0 + w(291)
      w(291)=w(65)*w(291)
      w(292)= - LL - 5.D0/3.D0 + w(220)
      w(292)=w(292)*w(186)
      w(292)= - 233.D0/9.D0 + w(292)
      w(292)=w(2)*w(292)
      w(293)=64.D0*w(74)
      w(294)=12.D0 - w(129)
      w(294)=z*w(294)
      w(294)=3.D0 + w(294)
      w(294)=w(294)*w(293)
      w(295)=128.D0*w(131)
      w(296)=w(94)*w(295)
      w(297)=192.D0*w(99)
      w(298)=w(19)*w(297)
      w(299)=w(61)*w(173)
      w(300)=32.D0*z
      w(301)=w(68)*w(300)
      w(100)=w(301) + w(100) + w(298) + w(273) + w(296) + w(294) + 
     & w(207) - 128.D0*w(299) + w(257) + w(172) + 96.D0*w(277) + 4.D0*
     & w(292) + w(113) + w(102) + w(265) + w(262) + 4.D0/3.D0*w(291) + 
     & w(255) + w(232) + 48.D0*w(290) + w(287) + w(251) + w(247) + 
     & w(146) + w(285) + w(283) - 736.D0/9.D0*w(171) + w(281) + w(279)
     &  + 16891.D0/81.D0 + w(280) + w(245) + w(278) + w(276) + w(229)
     &  + 32.D0*w(225) + 64.D0*w(130) + 8.D0*w(144) + 16.D0*w(152)
      w(100)=CF*w(100)
      w(102)=w(195) + 59.D0
      w(102)=w(102)*w(97)
      w(102)=w(102) + w(197) - 1.D0
      w(102)=1.D0/3.D0*w(102)
      w(113)=3.D0*w(101)
      w(130)= - w(102) - w(113)
      w(144)=w(7)*w(8)
      w(130)=w(130)*w(144)
      w(146)=2.D0*w(99)
      w(152)=w(73)*w(146)
      w(172)=142.D0/3.D0 + 89.D0*z
      w(172)=w(172)*w(129)
      w(207)=40.D0/3.D0*w(6)
      w(172)=w(207) + 266.D0/3.D0 + w(172)
      w(172)=w(172)*w(206)
      w(225)=2.D0 - w(154)
      w(225)=w(225)*w(244)
      w(229)=4.D0*w(30)
      w(232)= - w(99)*w(229)
      w(245)=887.D0 - 2695.D0/3.D0*z
      w(245)=z*w(245)
      w(245)= - 64.D0 + w(245)
      w(247)=w(83)*w(114)
      w(130)=4.D0*w(247) + w(232) + w(130) + w(225) + w(172) + 1.D0/3.D0
     & *w(245) + w(152)
      w(152)=76.D0/9.D0 - z
      w(152)=w(152)*w(97)
      w(172)=8.D0/9.D0*w(6)
      w(152)=w(172) + 11.D0/9.D0 + w(152)
      w(152)=w(152)*w(122)
      w(225)=220.D0 - 1883.D0/9.D0*z
      w(225)=z*w(225)
      w(225)= - 41.D0 + w(225)
      w(152)=1.D0/3.D0*w(225) + w(152)
      w(152)=LL*w(152)
      w(130)=2.D0*w(130) + w(152)
      w(130)=LL*w(130)
      w(152)=z - 17.D0
      w(152)=w(152)*w(97)
      w(225)=11.D0 + w(152)
      w(225)=w(56)*w(225)
      w(232)= - 11.D0 - w(97)
      w(232)=w(232)*w(97)
      w(232)=3.D0 + w(232)
      w(232)=w(47)*w(232)
      w(225)=w(225) + w(232)
      w(232)=8.D0*w(28)
      w(191)=80.D0*w(22) + w(232) + w(191) - 16.D0*w(31)
      w(191)=w(99)*w(191)
      w(245)=w(131)*w(165)
      w(247)=13.D0 - w(105)
      w(247)=w(247)*w(97)
      w(247)=w(205) + 17.D0 + w(247)
      w(247)=1.D0/3.D0*w(247) + w(245)
      w(247)=w(34)*w(247)
      w(251)=w(99)*w(109)
      w(255)=20.D0 + w(154)
      w(255)=w(255)*w(272)
      w(255)=59.D0 + w(255)
      w(255)=w(251) + 1.D0/3.D0*w(255) + w(205)
      w(255)=w(8)*w(255)
      w(257)=476.D0 + 541.D0*z
      w(257)=w(257)*z
      w(257)=w(257) - 50.D0 + 62.D0*w(6)
      w(255)= - 8.D0/9.D0*w(257) + w(255)
      w(255)=w(8)*w(255)
      w(262)=w(99)*z3
      w(255)=w(255) + 48.D0*w(262)
      w(255)=w(7)*w(255)
      w(265)=w(216)*w(110)
      w(265)=w(265) + 7.D0
      w(273)= - w(8)*w(265)
      w(276)=25.D0 - 32.D0/3.D0*z
      w(276)=z*w(276)
      w(273)=w(273) + w(135) + 2.D0 + w(276)
      w(273)=w(273)*w(147)
      w(276)= - w(263) - 32.D0*w(37) - 48.D0*w(87)
      w(276)=w(131)*w(276)
      w(278)=19.D0/3.D0*z
      w(279)=26.D0 - w(278)
      w(279)=w(279)*w(272)
      w(279)=25.D0 + w(279)
      w(279)=w(73)*w(279)
      w(280)=z - 3.D0
      w(281)= - w(81)*w(280)*w(300)
      w(283)= - 20.D0 - w(162)
      w(283)=z*w(283)
      w(283)= - 7.D0 + w(283)
      w(283)=w(283)*w(164)
      w(285)=13409.D0/3.D0 + 10324.D0*z
      w(285)=z*w(285)
      w(285)=3023.D0/3.D0 + w(285)
      w(283)=448.D0/27.D0*w(6) + 1.D0/9.D0*w(285) + w(283)
      w(283)=w(283)*w(122)
      w(285)=63.D0 + w(213)
      w(285)=w(285)*w(97)
      w(285)=13.D0 + w(285)
      w(287)=4.D0*w(71)
      w(285)=w(285)*w(287)
      w(290)=3.D0*w(8)
      w(291)= - 3.D0 - w(145)
      w(291)=w(291)*w(290)
      w(292)=60.D0*z
      w(294)= - 343.D0/3.D0 + w(292)
      w(294)=z*w(294)
      w(291)=w(291) - 14.D0/3.D0 + w(294)
      w(291)=w(291)*w(220)
      w(294)=5.D0/3.D0*z
      w(296)=3.D0 + w(294)
      w(296)=z*w(296)
      w(296)=2.D0/3.D0*w(6) + 2.D0 + w(296)
      w(256)=w(296)*w(256)
      w(296)=8.D0*w(65)
      w(298)= - 5.D0 + w(239)
      w(298)=w(298)*w(296)
      w(299)=98509.D0 - 290560.D0/3.D0*z
      w(299)=z*w(299)
      w(299)= - 10648.D0 + w(299)
      w(130)=w(130) + w(273) + w(298) + w(256) + w(255) - 80.D0*w(153)
     &  + 4.D0*w(247) + w(291) + w(285) + w(283) + w(281) + 80.D0/3.D0*
     & w(171) + 2.D0*w(279) + 1.D0/27.D0*w(299) + w(191) + w(276) + 8.D0
     & *w(225)
      w(130)=LL*w(130)
      w(191)= - 5569.D0 - 112981.D0/9.D0*z
      w(191)=w(191)*w(97)
      w(191)= - 5459.D0/3.D0 + w(191)
      w(225)=70.D0 - 43.D0*z
      w(225)=z*w(225)
      w(225)= - 28.D0/3.D0 + w(225)
      w(225)=w(225)*w(164)
      w(247)= - 82.D0/27.D0 + w(73)
      w(247)=w(247)*w(284)
      w(255)=w(265)*w(188)
      w(194)= - w(120)*w(194)
      w(256)=w(120)*z
      w(265)=w(256) + 1.D0
      w(273)=w(52)*w(265)
      w(191)=w(194) + 32.D0*w(273) + w(255) + w(247) + 1.D0/9.D0*w(191)
     &  + w(225)
      w(191)=w(8)*w(191)
      w(194)= - w(261) + w(259)
      w(194)=w(37)*w(194)
      w(225)=203.D0/3.D0*z
      w(247)=w(225) - 55.D0
      w(247)=w(247)*w(97)
      w(207)=w(247) - w(207) - 1.D0
      w(247)=w(81)*w(207)
      w(255)= - 13.D0/3.D0 + w(157)
      w(255)=w(255)*w(97)
      w(255)=23.D0/3.D0 + w(255)
      w(255)=w(65)*w(255)
      w(194)=w(194) + w(247) + w(255)
      w(247)= - 29.D0/3.D0 - w(272)
      w(247)=w(247)*w(97)
      w(255)=z - 7.D0
      w(255)=w(255)*w(97)
      w(255)=w(255) + 5.D0
      w(259)= - w(255)*w(122)
      w(273)=w(135) + 5.D0/3.D0
      w(247)=w(259) + w(247) - w(273)
      w(247)=w(56)*w(247)
      w(259)= - 29.D0/3.D0 - w(133)
      w(259)=w(259)*w(97)
      w(276)=z + 5.D0
      w(279)=w(276)*w(97)
      w(279)= - 1.D0 + w(279)
      w(279)=w(279)*w(122)
      w(259)=w(279) + w(259) - w(273)
      w(259)=w(47)*w(259)
      w(279)=w(97) + 5.D0/3.D0
      w(279)=w(279)*w(258)
      w(281)=w(279) + w(273)
      w(281)=w(52)*w(281)
      w(247)=w(259) + w(281) + w(247)
      w(259)= - 23.D0/3.D0 - w(133)
      w(259)=w(259)*w(97)
      w(259)= - w(101) + w(259) - w(273)
      w(259)=w(8)*w(259)
      w(281)=20.D0/3.D0*w(6)
      w(283)=w(281) - 4.D0
      w(284)=109.D0 + 359.D0/3.D0*z
      w(284)=z*w(284)
      w(284)=w(284) + w(283)
      w(259)=2.D0*w(284) + w(259)
      w(259)=w(8)*w(259)
      w(284)=733.D0 + 2443.D0/3.D0*z
      w(284)=w(284)*z
      w(284)=w(284) - 50.D0 + 94.D0/3.D0*w(6)
      w(259)=4.D0/3.D0*w(284) + w(259)
      w(259)=w(259)*w(206)
      w(285)= - 7.D0 - w(272)
      w(285)=w(285)*w(258)
      w(285)= - w(234) - 5.D0 + w(285)
      w(285)=w(285)*w(236)
      w(259)=w(259) + w(285)
      w(259)=w(7)*w(259)
      w(285)= - 11.D0 + z
      w(285)=w(285)*w(97)
      w(285)=7.D0 + w(285)
      w(285)=w(64)*w(285)
      w(255)=w(66)*w(255)
      w(255)= - w(285) + w(255) + w(67) + w(49)
      w(285)=w(46)*w(99)
      w(218)=w(218) + w(285)
      w(285)=w(212) + 47.D0/3.D0
      w(285)=w(285)*w(97)
      w(285)=w(285) + w(273)
      w(291)=w(285) + w(228)
      w(298)=w(11)*w(103)
      w(298)=w(298) - 8.D0*w(22)
      w(298)=w(291)*w(298)
      w(299)=w(124) + w(273)
      w(301)=38.D0/3.D0 + w(238)
      w(301)=w(301)*w(97)
      w(301)=w(301) + w(299)
      w(232)=w(301)*w(232)
      w(241)=w(241)*w(103)
      w(225)= - 136.D0 + w(225)
      w(225)=z*w(225)
      w(225)=w(241) - w(202) + 23.D0 + w(225)
      w(225)=w(225)*w(287)
      w(241)=4.D0/3.D0*z
      w(287)=698.D0 + 4091.D0/9.D0*z
      w(287)=w(287)*w(241)
      w(301)=34.D0/3.D0 + w(107)
      w(301)=w(301)*w(133)
      w(301)= - w(135) - 269.D0/3.D0 + w(301)
      w(301)=w(301)*w(109)
      w(302)= - w(5)*w(197)
      w(287)=w(301) + w(302) + 131.D0 + w(287)
      w(287)=z3*w(287)
      w(207)= - w(245) - w(207)
      w(207)=w(34)*w(207)*w(122)
      w(279)= - w(279) - w(299)
      w(279)=w(8)*w(279)
      w(254)=w(254) + w(279)
      w(229)=w(254)*w(229)
      w(254)=131.D0 - w(145)
      w(254)=w(254)*w(97)
      w(254)=25.D0 + w(254)
      w(279)= - z*w(276)
      w(279)= - 2.D0 + w(279)
      w(279)=w(279)*w(165)
      w(254)=w(279) + 1.D0/3.D0*w(254) + w(205)
      w(254)=w(254)*w(122)
      w(279)= - 368.D0 + 3425.D0/9.D0*z
      w(279)=z*w(279)
      w(254)=w(254) - 80.D0/3.D0*w(6) + 16.D0/3.D0 + w(279)
      w(254)=w(83)*w(254)
      w(246)=w(140)*w(246)
      w(279)=88.D0*w(92) + 40.D0*w(95) + 72.D0*w(89)
      w(279)=w(131)*w(279)
      w(299)=8.D0*w(27)
      w(301)=579401.D0/243.D0 - w(299)
      w(301)=z*w(301)
      w(301)=w(301) - 65491.D0/27.D0 - w(299)
      w(301)=w(301)*w(97)
      w(302)= - 656.D0 - 235.D0*z
      w(302)=w(302)*w(267)
      w(302)= - 11.D0 + w(302)
      w(302)=w(73)*w(302)
      w(303)= - 11.D0 - z
      w(303)=z*w(303)
      w(303)= - 2.D0 + w(303)
      w(288)=w(303)*w(288)
      w(303)=w(129) - 7.D0
      w(304)=w(303)*z
      w(304)=w(304) + 2.D0
      w(305)=16.D0*w(86)
      w(306)=w(304)*w(305)
      w(307)=w(29)*w(119)
      w(308)=w(38)*w(167)
      w(309)=755.D0/81.D0 - w(236)
      w(310)=61.D0 + 34.D0/3.D0*LL
      w(310)=LL*w(310)
      w(309)=4.D0*w(309) + 1.D0/9.D0*w(310)
      w(309)=w(309)*w(270)
      w(310)= - 66607.D0/9.D0 + 448.D0*z3
      w(309)=1.D0/27.D0*w(310) + w(309)
      w(309)=w(2)*w(309)
      w(310)=w(272) - 1.D0
      w(311)=z5*w(310)
      w(130)=w(309) + w(130) + w(254) - 48.D0*w(308) + w(229) + w(259)
     &  + w(307) + w(306) - 80.D0*w(79) + w(207) + w(288) + 24.D0*w(58)
     &  + 32.D0*w(14) + w(287) + w(225) + w(232) + w(191) - 628.D0/9.D0
     & *w(171) + 5.D0*w(302) + 68.D0*w(311) + w(301) + 3077.D0/9.D0 - 
     & w(299) + w(279) + w(298) + w(246) - 16.D0*w(218) - 8.D0*w(255)
     &  + 4.D0*w(247) + 2.D0*w(194)
      w(191)=19.D0 - w(258)
      w(191)=w(191)*w(110)
      w(191)= - 7.D0 + w(191)
      w(194)= - 16.D0 - w(157)
      w(194)=z*w(194)
      w(194)= - 3.D0 + w(194)
      w(194)=w(194)*w(165)
      w(207)=13.D0 - w(157)
      w(207)=w(207)*w(97)
      w(207)= - 11.D0 + w(207)
      w(207)=w(207)*w(160)
      w(191)=w(207) + w(194) + 7.D0/3.D0*w(191) + w(205)
      w(191)=LL*w(191)
      w(194)= - 47.D0/3.D0 - w(145)
      w(194)=w(194)*w(97)
      w(207)=w(280)*w(97)
      w(207)=w(207) + 3.D0
      w(218)= - w(207)*w(122)
      w(194)=w(218) + w(194) - w(273)
      w(194)=w(8)*w(194)
      w(218)= - 5.D0 + 41.D0/3.D0*z
      w(218)=w(218)*w(107)
      w(218)=w(218) + w(283)
      w(191)=w(191) + 2.D0/3.D0*w(218) + w(194)
      w(191)=w(55)*w(191)
      w(194)=w(261) + w(135)
      w(218)=34.D0*z
      w(225)= - 83.D0/3.D0 + w(218)
      w(225)=w(225)*w(97)
      w(225)= - w(134) - 19.D0/3.D0 + w(225) - w(194)
      w(225)=w(40)*w(225)
      w(229)=28.D0*z
      w(232)=59.D0/3.D0 + w(229)
      w(232)=w(232)*w(97)
      w(246)=z*w(201)
      w(246)=2.D0 + w(246)
      w(246)=w(246)*w(103)
      w(247)=w(129) + 23.D0
      w(254)=w(247)*w(97)
      w(254)=23.D0 + w(254)
      w(254)=w(254)*w(160)
      w(232)=w(254) + w(246) + w(232) + w(273)
      w(232)=w(63)*w(232)
      w(246)= - 295.D0 + 364.D0*z
      w(246)=w(246)*w(97)
      w(246)= - 11.D0 + w(246)
      w(254)=24.D0*w(6)
      w(246)= - 32.D0*w(111) + w(245) + 1.D0/3.D0*w(246) - w(254)
      w(246)=w(36)*w(246)
      w(191)=w(191) + w(225) + w(232) + w(246)
      w(225)= - 4.D0 + w(278)
      w(225)=z*w(225)
      w(185)=w(225) - w(185)
      w(185)=2.D0*w(185) + w(167)
      w(185)=w(8)*w(185)
      w(195)=59.D0 - w(195)
      w(195)=w(195)*w(97)
      w(195)=w(197) + 1.D0 + w(195)
      w(195)=1.D0/3.D0*w(195) - w(245)
      w(195)=w(195)*w(198)
      w(225)=604.D0 - 617.D0*z
      w(225)=z*w(225)
      w(225)=w(199) + 70.D0 + w(225)
      w(185)=w(195) + 1.D0/9.D0*w(225) + w(185)
      w(185)=w(185)*w(160)
      w(195)= - 59.D0/3.D0 - w(179)
      w(195)=z*w(195)
      w(195)= - 10.D0/3.D0*w(167) + w(135) + 58.D0/3.D0 + w(195)
      w(195)=w(8)*w(195)
      w(225)=1583.D0*z
      w(232)= - 1588.D0 + w(225)
      w(232)=z*w(232)
      w(232)=197.D0 + w(232)
      w(189)=w(195) + 1.D0/3.D0*w(232) - w(189)
      w(189)=w(8)*w(189)
      w(195)=13409.D0 - 9884.D0*z
      w(195)=z*w(195)
      w(195)= - 3625.D0 + w(195)
      w(195)=1.D0/3.D0*w(195) - 164.D0*w(6)
      w(185)=w(185) + 136.D0*w(204) + 1.D0/9.D0*w(195) + w(189)
      w(185)=w(185)*w(160)
      w(189)= - 596.D0 + 1619.D0/3.D0*z
      w(189)=w(189)*w(267)
      w(195)= - 16.D0 + 65.D0/3.D0*z
      w(232)=w(195)*w(267)
      w(246)=w(172) + 1.D0
      w(232)=w(232) - w(246)
      w(255)=w(232)*w(165)
      w(189)=w(255) - 308.D0/9.D0*w(6) + 53.D0 + w(189)
      w(189)=w(8)*w(189)
      w(255)=121.D0 - 3208.D0/27.D0*z
      w(255)=z*w(255)
      w(255)= - 35.D0/3.D0 + w(255)
      w(255)=5.D0*w(255) + 1388.D0/27.D0*w(6)
      w(189)=4.D0*w(255) + w(189)
      w(189)=w(8)*w(189)
      w(255)=349.D0 - 437.D0*z
      w(255)=w(255)*w(97)
      w(255)=80.D0*w(6) + 41.D0 + w(255)
      w(255)=1.D0/3.D0*w(255) - w(261)
      w(255)=w(255)*w(208)
      w(259)= - 18101.D0 + 18046.D0*z
      w(259)=z*w(259)
      w(259)= - 1465.D0*w(6) + 1286.D0 + w(259)
      w(185)=w(185) + w(255) + 2.D0/27.D0*w(259) + w(189)
      w(185)=w(32)*w(185)
      w(189)=40.D0*z
      w(255)= - 329.D0 + w(189)
      w(255)=z*w(255)
      w(255)= - 86.D0 + w(255)
      w(255)=1.D0/3.D0*w(255) - w(205)
      w(259)= - 1.D0 + w(212)
      w(259)=w(259)*w(97)
      w(259)=23.D0 + w(259)
      w(259)=w(8)*w(259)
      w(255)=2.D0*w(255) + w(259)
      w(255)=w(8)*w(255)
      w(224)= - 3.D0 - w(224)
      w(224)=w(8)*w(224)
      w(259)=w(268)*w(160)
      w(261)=12.D0 + 25.D0/3.D0*z
      w(261)=z*w(261)
      w(224)=w(259) + w(261) + w(224)
      w(224)=w(224)*w(186)
      w(259)=53.D0 - 427.D0/9.D0*z
      w(259)=w(259)*w(272)
      w(224)=w(224) + w(255) + 416.D0/9.D0*w(6) - 451.D0/3.D0 + w(259)
      w(224)=LL*w(224)
      w(255)=w(114)*w(165)
      w(259)= - 31.D0 + w(214)
      w(259)=w(259)*w(97)
      w(259)= - w(255) + 7.D0 + w(259)
      w(259)=w(259)*w(109)
      w(261)=1964.D0 + 1501.D0/3.D0*z
      w(261)=z*w(261)
      w(261)=31.D0 + w(261)
      w(259)=w(259) + 1.D0/3.D0*w(261) + 52.D0*w(6)
      w(259)=w(8)*w(259)
      w(261)=z3*w(268)
      w(268)= - 47.D0 + w(112)
      w(268)=z*w(268)
      w(268)= - 1388.D0/9.D0*w(6) + 125.D0 + 41.D0*w(268)
      w(224)=w(224) + 40.D0*w(261) + 2.D0/3.D0*w(268) + w(259)
      w(259)=2.D0*w(77)
      w(224)=w(224)*w(259)
      w(261)= - 29.D0/3.D0 - w(253)
      w(261)=w(261)*w(97)
      w(156)= - w(231) + w(156) + w(261) - w(273)
      w(156)=w(25)*w(156)
      w(261)=w(135) - 11.D0/3.D0
      w(174)= - 100.D0/3.D0 + w(174)
      w(174)=z*w(174)
      w(174)= - 10.D0*w(111) + w(255) + w(174) - w(261)
      w(268)=8.D0*w(85)
      w(174)=w(174)*w(268)
      w(279)=64.D0*z
      w(287)= - 173.D0 + w(279)
      w(287)=w(287)*w(97)
      w(287)= - 5.D0 + w(287)
      w(288)=z + 8.D0
      w(298)=z*w(288)
      w(298)= - 1.D0 + w(298)
      w(298)=w(298)*w(122)
      w(287)=w(298) + 1.D0/3.D0*w(287) - w(205)
      w(287)=w(8)*w(287)
      w(159)=w(159) + 7.D0
      w(298)=w(159)*w(122)
      w(102)=w(102) + w(298)
      w(102)=LL*w(102)
      w(102)=w(102) + 4.D0/9.D0*w(257) + w(287)
      w(102)=w(102)*w(160)
      w(257)=35.D0/3.D0 + w(214)
      w(257)=w(257)*w(97)
      w(287)=w(8)*w(265)
      w(257)=8.D0/3.D0*w(287) + w(257) + w(273)
      w(257)=w(8)*w(257)
      w(287)= - 37.D0 - 323.D0/3.D0*z
      w(287)=z*w(287)
      w(287)=w(287) - w(283)
      w(257)=4.D0/3.D0*w(287) + w(257)
      w(257)=w(8)*w(257)
      w(287)=w(249)*w(110)
      w(298)=5.D0 + w(287)
      w(298)=w(298)*w(220)
      w(102)=w(102) + w(298) - 4.D0/9.D0*w(284) + w(257)
      w(257)=2.D0*w(41)
      w(102)=w(102)*w(257)
      w(284)=103.D0*z
      w(298)=94.D0 - w(284)
      w(298)=z*w(298)
      w(298)= - 2.D0 + w(298)
      w(190)= - 6.D0*w(111) + 1.D0/3.D0*w(298) + w(190)
      w(190)=w(190)*w(186)
      w(298)=31.D0*z
      w(299)=82.D0/3.D0 - w(298)
      w(299)=z*w(299)
      w(299)=w(135) - 8.D0/3.D0 + w(299)
      w(299)=w(299)*w(165)
      w(301)= - 160.D0 + 1321.D0/9.D0*z
      w(301)=z*w(301)
      w(190)=w(190) + w(299) - 4.D0/9.D0*w(6) + 5.D0/3.D0 + w(301)
      w(299)=2.D0*w(39)
      w(190)=w(190)*w(299)
      w(291)=w(231) - w(291)
      w(291)=w(16)*w(291)
      w(301)= - 2191.D0 + 5248.D0/3.D0*z
      w(301)=w(301)*w(267)
      w(302)=5.D0 - w(298)
      w(302)=w(302)*w(170)
      w(301)=w(302) + 82.D0 + w(301)
      w(302)=w(310)*w(198)
      w(306)=1.D0/3.D0 + w(157)
      w(306)=z*w(306)
      w(302)=w(302) - 1.D0 + w(306)
      w(302)=w(302)*w(270)
      w(306)= - 24.D0 - 4303.D0/9.D0*z
      w(306)=z*w(306)
      w(302)=w(302) - 1303.D0/9.D0 + w(306)
      w(302)=LL*w(302)
      w(301)=1.D0/3.D0*w(301) + w(302)
      w(301)=w(301)*w(269)
      w(302)= - 20.D0/3.D0 - w(203)
      w(302)=w(302)*w(97)
      w(306)= - w(117)*w(270)
      w(302)=w(306) + 37.D0/3.D0 + w(302)
      w(302)=w(302)*w(160)
      w(306)= - 4.D0 + 75.D0*z
      w(306)=z*w(306)
      w(302)=w(302) - 35.D0/3.D0 + w(306)
      w(302)=w(302)*w(274)
      w(306)= - 7.D0 - w(157)
      w(306)=w(50)*w(306)
      w(307)= - 17.D0 - z
      w(307)=w(80)*w(307)
      w(306)=w(307) + w(306)
      w(306)=w(300)*w(306)
      w(119)=w(119)*w(243)
      w(243)= - w(61)*w(280)
      w(243)=w(68) + w(243)
      w(243)=w(279)*w(243)
      w(280)=96.D0*w(84) + 16.D0*w(93) + 48.D0*w(91)
      w(280)=w(131)*w(280)
      w(275)=256.D0*w(57) + 128.D0*w(275)
      w(265)=w(265)*w(275)
      w(275)=w(96)*z
      w(248)=w(275)*w(248)
      w(297)=w(13)*w(297)
      w(307)=z + 9.D0
      w(308)= - z*w(307)
      w(308)= - 3.D0 + w(308)
      w(293)=w(308)*w(293)
      w(308)= - w(19)*w(140)
      w(102)=w(102) + w(308) - 96.D0*w(60) + w(302) + w(174) + 8.D0*
     & w(156) + w(293) + w(224) + w(301) + 192.D0*w(17) + w(297) + 16.D0
     & *w(291) - 48.D0*w(277) + w(190) + w(248) + w(185) + w(265) + 
     & w(280) + w(243) + w(119) + w(306) + 2.D0*w(130) + 4.D0*w(191)
      w(119)=2.D0*CA
      w(102)=w(102)*w(119)
      w(130)=3.D0 + w(129)
      w(130)=w(130)*w(97)
      w(130)=3.D0 + w(130)
      w(130)=CF*w(130)
      w(156)=w(307)*w(97)
      w(174)= - CA*w(156)
      w(130)=w(130) + w(174)
      w(130)=w(59)*w(130)
      w(174)=CA - CF
      w(185)= - w(12)*w(101)*w(174)
      w(100)=32.D0*w(130) + 384.D0*w(185) + w(102) + w(100)
      w(100)=CA*w(100)
      w(102)= - 11.D0 + w(112)
      w(102)=w(102)*w(110)
      w(130)=w(110) - 3.D0
      w(185)=w(130)*w(97)
      w(185)=3.D0 + w(185)
      w(190)=16.D0*w(73)
      w(185)=w(185)*w(190)
      w(191)=13.D0 - 50.D0*z
      w(191)=w(191)*w(97)
      w(191)= - 27.D0 + w(191)
      w(191)=w(191)*w(122)
      w(224)=w(179) + 1.D0
      w(224)=w(224)*w(97)
      w(243)=3.D0 + w(224)
      w(243)=w(243)*w(170)
      w(248)=w(157) + 5.D0
      w(265)=w(248)*z
      w(265)=w(265) + 2.D0
      w(280)=w(101) + 2.D0*w(265)
      w(291)=w(193)*w(7)
      w(293)=w(280)*w(291)
      w(297)=16.D0*w(83)
      w(301)= - w(106)*w(297)
      w(302)=w(97) - 11.D0
      w(306)=w(110)*w(8)
      w(307)=1.D0/3.D0*w(302) + w(306)
      w(307)=w(307)*w(160)
      w(102)=w(307) + w(301) + w(293) + w(243) + w(191) + w(185) - 27.D0
     &  + w(102)
      w(102)=LL*w(102)
      w(185)= - LL*w(177)
      w(191)=7.D0 - w(157)
      w(191)=z*w(191)
      w(191)= - 7.D0 + w(191)
      w(113)=w(185) + 2.D0*w(191) - w(113)
      w(113)=w(55)*w(113)
      w(185)=3.D0*w(167)
      w(191)= - z*w(158)
      w(191)= - 2.D0 + w(191)
      w(191)=2.D0*w(191) - w(185)
      w(191)=w(34)*w(191)
      w(209)= - 1.D0 - w(209)
      w(209)=w(63)*w(209)
      w(113)=w(209) + w(191) + w(113)
      w(150)= - w(295)*w(150)
      w(191)=w(98) - 1.D0
      w(161)= - 64.D0*w(47) - 128.D0*w(56) - w(161)
      w(161)=w(191)*w(161)
      w(209)= - w(222)*w(110)
      w(209)=5.D0 + w(209)
      w(209)=w(73)*w(209)
      w(118)= - 7.D0 - w(118)
      w(118)=w(118)*w(164)
      w(222)= - 79.D0 - 94.D0*z
      w(222)=z*w(222)
      w(118)=w(118) - 67.D0 + w(222)
      w(118)=w(118)*w(165)
      w(222)=5.D0 + w(286)
      w(222)=w(8)*w(222)
      w(243)= - 11.D0 + w(105)
      w(243)=z*w(243)
      w(243)=1.D0 + w(243)
      w(222)=3.D0*w(243) + w(222)
      w(222)=w(222)*w(170)
      w(243)=w(99)*w(206)
      w(286)=9.D0 + w(129)
      w(286)=z*w(286)
      w(286)=4.D0 + w(286)
      w(286)=2.D0*w(286) + w(243)
      w(286)=w(8)*w(286)
      w(247)=w(247)*z
      w(247)=w(247) + 20.D0
      w(247)=2.D0*w(247)
      w(286)=w(247) + w(286)
      w(286)=w(286)*w(291)
      w(291)= - 1.D0 - w(253)
      w(291)=w(291)*w(97)
      w(291)= - 3.D0 + w(291)
      w(291)=w(291)*w(296)
      w(293)=w(191)*w(122)
      w(301)=41.D0 - 24.D0*z
      w(301)=z*w(301)
      w(293)=w(293) - 6.D0 + w(301)
      w(293)=w(293)*w(297)
      w(301)=z*w(240)
      w(307)=w(129) - 2.D0
      w(308)=z*w(307)
      w(308)=4.D0 + w(308)
      w(308)=w(71)*w(308)
      w(102)=w(102) + w(293) - 160.D0*w(148) + w(291) + w(286) + w(222)
     &  + 64.D0*w(308) + w(118) + 24.D0*w(209) - 959.D0/2.D0 - 81.D0*
     & w(301) + w(161) + w(150) + 32.D0*w(113)
      w(102)=LL*w(102)
      w(113)= - 21.D0 - w(300)
      w(113)=w(113)*w(129)
      w(118)=9.D0 - w(187)
      w(118)=w(118)*w(164)
      w(148)=16.D0*w(81)
      w(150)= - w(191)*w(148)
      w(161)=1.D0 - w(224)
      w(161)=w(161)*w(165)*w(73)
      w(113)=w(161) + w(150) + w(118) - 59.D0 + w(113)
      w(113)=w(113)*w(122)
      w(118)=107.D0 - 92.D0*z
      w(118)=w(118)*w(97)
      w(150)=8.D0/3.D0*w(167)
      w(161)= - w(150) + 11.D0 - w(214)
      w(161)=w(8)*w(161)
      w(118)=w(161) - 63.D0 + w(118)
      w(118)=w(8)*w(118)
      w(161)=w(167) + 1.D0
      w(187)= - w(287) - w(161)
      w(187)=w(187)*w(122)
      w(209)= - w(132) + w(117)
      w(209)=w(209)*w(108)
      w(222)= - 37.D0 + 25.D0*z
      w(222)=w(222)*z
      w(187)=w(209) + w(187) - 17.D0 - w(222)
      w(187)=w(187)*w(160)
      w(209)=47.D0*z
      w(224)=20.D0 - w(209)
      w(224)=z*w(224)
      w(224)=29.D0 + w(224)
      w(118)=w(187) + 56.D0*w(204) + 2.D0*w(224) + w(118)
      w(118)=w(118)*w(160)
      w(187)=w(131)*w(206)
      w(224)=w(303)*w(241)
      w(224)=w(187) + 1.D0 + w(224)
      w(224)=w(224)*w(122)
      w(286)=23.D0 - w(272)
      w(286)=w(286)*w(97)
      w(224)=w(224) + 1.D0 + w(286)
      w(224)=w(8)*w(224)
      w(286)= - 233.D0 + 264.D0*z
      w(286)=z*w(286)
      w(286)= - 6.D0 + w(286)
      w(224)=2.D0*w(286) + w(224)
      w(224)=w(8)*w(224)
      w(286)=w(110) - 17.D0/3.D0
      w(287)= - w(286)*w(97)
      w(287)= - 14.D0/3.D0*w(167) - 1.D0/3.D0 + w(287)
      w(287)=w(287)*w(170)
      w(291)=80.D0*z
      w(293)=111.D0 - w(291)
      w(293)=z*w(293)
      w(293)= - 41.D0 + w(293)
      w(118)=w(118) + w(287) + 2.D0*w(293) + w(224)
      w(118)=w(118)*w(192)
      w(130)=w(130)*w(110)
      w(130)= - 5.D0 + w(130)
      w(130)= - w(134) + 3.D0*w(130) + w(245)
      w(130)=w(40)*w(130)
      w(224)=w(117)*w(97)
      w(224)=w(224) + 1.D0
      w(287)=w(68)*w(224)
      w(130)=w(130) + w(287)
      w(287)=3.D0*w(111)
      w(293)=w(129) - 4.D0
      w(293)=w(293)*z
      w(301)=w(293) - 1.D0
      w(303)= - w(287) - 3.D0*w(301) - w(132)
      w(303)=w(303)*w(160)
      w(158)=w(158)*w(272)
      w(308)=w(157) - 2.D0
      w(308)=w(308)*z
      w(309)= - 1.D0 + w(308)
      w(309)=w(309)*w(122)
      w(158)=w(303) + w(309) + 7.D0 + w(158)
      w(303)=16.D0*w(39)
      w(158)=w(158)*w(303)
      w(309)= - 7.D0 - w(133)
      w(309)=w(309)*w(97)
      w(311)= - w(224)*w(198)
      w(309)=w(311) - 1.D0 + w(309)
      w(309)=LL*w(309)
      w(311)=w(107) + 2.D0
      w(312)= - w(311)*w(214)
      w(309)=w(309) - 25.D0 + w(312)
      w(309)=w(309)*w(186)
      w(312)= - 67.D0 + 176.D0*z
      w(312)=w(312)*w(157)
      w(313)=56.D0*z
      w(314)=17.D0 - w(313)
      w(314)=w(314)*w(97)
      w(314)= - 17.D0 + w(314)
      w(314)=w(314)*w(208)
      w(309)=w(309) + w(312) + w(314)
      w(309)=w(309)*w(269)
      w(312)=23.D0 - w(179)
      w(312)=w(312)*w(97)
      w(314)=w(97) - 3.D0
      w(315)=w(314)*w(97)
      w(316)=1.D0 - w(315)
      w(316)=w(8)*w(316)
      w(312)=w(316) - 17.D0 + w(312)
      w(312)=w(312)*w(122)
      w(242)=w(242) + 1.D0
      w(316)=w(242)*w(122)
      w(317)=4.D0/3.D0*LL
      w(318)= - w(173)*w(317)
      w(319)= - z*w(114)
      w(318)=w(318) - w(316) - 1.D0 + w(319)
      w(318)=w(318)*w(186)
      w(319)=w(154) + 5.D0
      w(320)= - w(319)*w(97)
      w(312)=w(318) + w(312) + 11.D0 + w(320)
      w(312)=LL*w(312)
      w(318)=w(106)*w(122)
      w(289)=w(318) + w(289)
      w(289)=w(8)*w(289)
      w(320)= - 11.D0 - w(258)
      w(289)=2.D0*w(320) + w(289)
      w(289)=w(8)*w(289)
      w(320)= - 5.D0 + w(264)
      w(320)=w(320)*w(97)
      w(320)=5.D0 + w(320)
      w(320)=w(320)*w(220)
      w(321)=271.D0 - 112.D0*z
      w(321)=z*w(321)
      w(289)=w(312) + w(320) + w(289) + 24.D0 + w(321)
      w(289)=w(289)*w(211)
      w(312)=16.D0*w(111)
      w(320)= - 17.D0 + w(203)
      w(320)=w(320)*w(110)
      w(320)= - w(312) + w(245) - 11.D0 + w(320)
      w(268)=w(320)*w(268)
      w(320)=w(129) - 8.D0
      w(321)=w(320)*z
      w(322)=w(321) + 1.D0
      w(322)=w(185) + 2.D0*w(322)
      w(287)=w(287) + w(322)
      w(287)=w(36)*w(287)
      w(323)=w(129) + 1.D0
      w(324)=w(323)*w(97)
      w(324)=w(324) - 1.D0
      w(325)=w(72)*w(324)
      w(287)=w(325) + w(287)
      w(325)=w(37)*w(167)
      w(277)=w(277) + w(325)
      w(238)= - 11.D0 + w(238)
      w(238)=w(238)*w(97)
      w(224)=w(224)*w(290)
      w(224)=w(224) - 5.D0 + w(238)
      w(168)=w(224)*w(168)
      w(224)=431.D0/3.D0 - w(279)
      w(224)=w(224)*w(97)
      w(238)=19.D0 - w(300)
      w(238)=z*w(238)
      w(238)= - 1.D0 + w(238)
      w(238)=w(238)*w(165)
      w(224)=w(238) - 277.D0/3.D0 + w(224)
      w(224)=w(224)*w(236)
      w(236)= - 1.D0 - w(214)
      w(236)=w(236)*w(97)
      w(184)= - w(184)*w(160)
      w(238)= - 23.D0 - w(189)
      w(238)=z*w(238)
      w(184)=w(184) - 4.D0 + w(238)
      w(184)=w(184)*w(186)
      w(184)=w(184) + 13.D0 + w(236)
      w(184)=w(184)*w(274)
      w(236)= - w(288)*w(97)
      w(236)=3.D0 + w(236)
      w(238)=3.D0 + w(176)
      w(238)=w(8)*w(238)
      w(236)=2.D0*w(236) + w(238)
      w(236)=w(8)*w(236)
      w(238)=w(173)*w(122)
      w(265)=w(238) - w(265)
      w(279)=w(265)*w(160)
      w(236)=w(279) - w(247) + w(236)
      w(247)=16.D0*LL
      w(236)=w(41)*w(236)*w(247)
      w(279)=w(176) - 1.D0
      w(288)=32.D0*w(80) - 64.D0*z5
      w(279)=w(279)*w(288)
      w(288)= - w(93) - w(94) + w(89)
      w(288)=w(295)*w(288)
      w(295)=281.D0*z
      w(325)= - 359.D0 + w(295)
      w(325)=w(325)*w(97)
      w(326)=39.D0 + w(110)
      w(326)=w(326)*w(164)
      w(320)=w(320)*w(97)
      w(320)=3.D0 + w(320)
      w(148)=w(320)*w(148)
      w(320)= - w(34)*w(322)*w(193)
      w(322)=w(272) - 7.D0
      w(322)=w(322)*w(97)
      w(322)=w(322) + 7.D0
      w(305)=w(322)*w(305)
      w(327)=160.D0*w(131)
      w(327)=w(95)*w(327)
      w(328)= - 7.D0 + w(105)
      w(328)=z*w(328)
      w(328)= - 1.D0 + w(328)
      w(328)=w(328)*w(296)
      w(329)= - w(137) + w(117)
      w(263)=w(329)*w(263)
      w(329)=w(249)*z
      w(330)=w(329) + 1.D0
      w(331)=w(330)*w(8)
      w(332)= - w(331) - 3.D0 - w(98)
      w(332)=w(332)*w(165)
      w(333)= - 95.D0 + w(313)
      w(333)=z*w(333)
      w(332)=w(332) + 29.D0 + w(333)
      w(332)=w(332)*w(147)
      w(333)=w(92)*w(149)
      w(334)=w(74)*w(173)
      w(335)=32.D0*w(106)
      w(336)= - w(84)*w(335)
      w(102)=w(236) + w(336) + w(184) + w(268) + 192.D0*w(334) + w(289)
     &  + w(309) + w(158) + w(118) + w(333) + w(332) + w(263) + w(328)
     &  + w(327) + w(305) + w(320) + w(224) + w(168) + w(113) + w(148)
     &  + w(326) + 133.D0 + w(325) + w(288) + w(102) + w(279) - 64.D0*
     & w(277) + 32.D0*w(287) + 8.D0*w(130)
      w(113)=CF**2
      w(102)=w(113)*w(102)
      w(118)= - w(52) + w(47) + w(56)
      w(130)=w(22) - w(28)
      w(148)=8.D0/3.D0*w(130) + 4.D0/3.D0*w(118)
      w(148)=w(99)*w(148)
      w(158)= - 38.D0 - 275.D0/3.D0*z
      w(158)=w(158)*w(97)
      w(168)=13.D0 + w(229)
      w(168)=w(168)*w(122)
      w(158)=w(168) - 11.D0 + w(158)
      w(158)=z3*w(158)
      w(168)=16.D0 + 37.D0/3.D0*z
      w(168)=w(168)*z
      w(168)=w(168) + 7.D0/3.D0
      w(172)=w(172) + w(168)
      w(172)=w(73)*w(172)
      w(184)=2.D0*w(131)
      w(224)= - w(81)*w(184)
      w(236)=w(97) + 1.D0/3.D0
      w(263)= - w(236)*w(164)
      w(268)=215.D0 + 3392.D0/3.D0*z
      w(268)=z*w(268)
      w(268)=37.D0 + w(268)
      w(263)=1.D0/27.D0*w(268) + w(263)
      w(263)=w(8)*w(263)
      w(268)=w(34)*w(132)
      w(277)=691.D0 - 5854.D0/9.D0*z
      w(277)=z*w(277)
      w(277)= - 59.D0 + w(277)
      w(279)=w(71)*w(272)
      w(148)= - 8.D0/3.D0*w(153) + w(268) + 1.D0/9.D0*w(158) + w(279)
     &  + w(263) + w(224) + 2.D0/27.D0*w(277) + w(148) + w(172)
      w(158)=2.D0*nf
      w(148)=w(148)*w(158)
      w(172)=w(131)*w(34)
      w(152)= - 7.D0 + w(152)
      w(152)=w(152)*w(164)
      w(224)= - 452.D0 - w(295)
      w(224)=w(224)*w(97)
      w(224)= - 79.D0 + w(224)
      w(224)=w(224)*w(115)
      w(155)=w(155)*w(220)
      w(263)= - 2894.D0/9.D0 + 365.D0*z
      w(263)=z*w(263)
      w(152)= - 8.D0*w(172) + w(155) + w(224) + w(152) + 95.D0/9.D0 + 
     & w(263)
      w(155)= - 67.D0 + 148.D0*z
      w(155)=w(155)*w(110)
      w(155)= - 119.D0 + w(155)
      w(224)=w(267) - 3.D0
      w(224)=w(224)*w(97)
      w(224)=w(224) - 1.D0
      w(263)=w(224)*w(164)
      w(268)=16.D0/3.D0*z3
      w(277)=w(323)*w(268)
      w(151)= - 269.D0 + w(151)
      w(151)=z*w(151)
      w(151)= - 2.D0 + w(151)
      w(151)=w(8)*w(151)
      w(151)= - 8.D0/3.D0*w(172) + w(277) + 2.D0/27.D0*w(151) + 1.D0/27.
     & D0*w(155) + w(263)
      w(151)=nf*w(151)
      w(155)=w(99)*nf
      w(263)=w(206)*w(155)
      w(263)=w(101) + w(263)
      w(263)=w(263)*w(141)
      w(277)=w(114)*w(122)
      w(279)= - 8.D0 + 31.D0/3.D0*z
      w(279)=w(279)*z
      w(279)=w(279) - 1.D0
      w(287)=w(279) - w(277)
      w(288)=nf + 7.D0
      w(289)=2.D0/9.D0*LL
      w(295)=w(289)*w(287)*w(288)
      w(305)= - 68.D0 - 79.D0*z
      w(305)=z*w(305)
      w(305)= - 11.D0 + w(305)
      w(305)=w(305)*w(165)
      w(309)= - 596.D0 + 1769.D0/3.D0*z
      w(309)=z*w(309)
      w(305)=w(305) + 32.D0 + w(309)
      w(309)=20.D0/9.D0 - z
      w(309)=z*w(309)
      w(309)=5.D0/9.D0 + w(309)
      w(309)=w(309)*w(165)
      w(320)=56.D0 - 205.D0/3.D0*z
      w(320)=z*w(320)
      w(320)= - 14.D0 + w(320)
      w(309)=1.D0/9.D0*w(320) + w(309)
      w(309)=nf*w(309)
      w(263)=w(295) + w(263) + 1.D0/9.D0*w(305) + w(309)
      w(263)=LL*w(263)
      w(295)=w(136) + 5.D0
      w(305)= - 4.D0/3.D0*w(295) - w(101)
      w(309)=nf*w(8)
      w(305)=w(305)*w(309)
      w(320)=w(275) + 5.D0
      w(323)= - w(320)*w(219)
      w(305)=w(323) + w(305)
      w(305)=w(7)*w(305)
      w(323)=w(83)*w(131)
      w(151)=w(263) + 4.D0/3.D0*w(323) + 2.D0/3.D0*w(305) + 1.D0/3.D0*
     & w(152) + w(151)
      w(151)=LL*w(151)
      w(152)= - 2.D0*w(232) + w(187)
      w(152)=w(8)*w(152)
      w(187)=w(187) - w(232)
      w(187)=w(8)*w(187)
      w(263)=29.D0 + 64.D0*w(256)
      w(187)=1.D0/27.D0*w(263) + w(187)
      w(187)=nf*w(187)
      w(263)=w(131)*nf
      w(305)=7.D0*w(131) + w(263)
      w(305)=LL*w(305)
      w(323)=w(110)*w(120)
      w(325)=w(323) + 5.D0
      w(326)= - nf*w(325)
      w(305)=w(305) + w(326) + 5.D0 + 28.D0*w(256)
      w(289)=w(305)*w(289)
      w(305)= - 221.D0 + 230.D0*z
      w(305)=z*w(305)
      w(305)=67.D0 + w(305)
      w(152)=w(289) + w(187) + 1.D0/27.D0*w(305) + w(152)
      w(152)=w(152)*w(160)
      w(187)=w(232)*w(8)
      w(232)= - 19.D0 + 499.D0/27.D0*z
      w(232)=w(232)*w(97)
      w(289)=80.D0/27.D0*w(6)
      w(232)= - w(187) - w(289) + 11.D0/3.D0 + w(232)
      w(232)=w(8)*w(232)
      w(305)=593.D0 - 602.D0*z
      w(305)=z*w(305)
      w(305)= - 391.D0 + w(305)
      w(305)=1.D0/9.D0*w(305) - w(182)
      w(232)=10.D0/9.D0*w(204) + 1.D0/9.D0*w(305) + w(232)
      w(305)= - 29.D0 + 254.D0/9.D0*z
      w(305)=w(305)*w(267)
      w(326)= - 1.D0 + 20.D0/27.D0*w(6)
      w(305)=w(305) - w(326)
      w(187)=4.D0*w(305) - w(187)
      w(187)=w(8)*w(187)
      w(305)= - 1.D0 + w(128)
      w(187)=20.D0/9.D0*w(204) + 4.D0/3.D0*w(305) + w(187)
      w(187)=nf*w(187)
      w(152)=w(152) + 2.D0*w(232) + w(187)
      w(152)=w(32)*w(152)
      w(187)= - 3.D0 + 34.D0/9.D0*z
      w(187)=w(187)*w(97)
      w(187)=w(187) - w(246)
      w(224)= - w(8)*w(224)
      w(224)=w(224) + w(187)
      w(232)= - w(210)*w(109)
      w(187)=w(232) + w(187)
      w(187)=nf*w(187)
      w(187)=2.D0*w(224) + w(187)
      w(187)=LL*w(187)
      w(224)=23.D0/3.D0*z
      w(232)=w(224) + 32.D0
      w(232)=w(232)*z
      w(232)= - w(232) - w(135) + w(277) - 5.D0
      w(246)=w(232)*w(122)
      w(277)=2.D0/3.D0*z
      w(305)=179.D0 - 562.D0/3.D0*z
      w(305)=w(305)*w(277)
      w(246)=w(246) + w(266) - 13.D0 + w(305)
      w(232)=w(232)*w(206)
      w(305)=29.D0/3.D0*z
      w(327)=1.D0 - 10.D0/9.D0*z
      w(327)=w(327)*w(305)
      w(326)=w(327) + w(326)
      w(232)=2.D0*w(326) + w(232)
      w(232)=nf*w(232)
      w(187)=w(187) + 1.D0/3.D0*w(246) + w(232)
      w(187)=w(187)*w(259)
      w(232)= - 52.D0/3.D0 + w(129)
      w(232)=w(232)*w(277)
      w(246)=1.D0/3.D0*nf
      w(264)= - 22.D0 + w(264)
      w(264)=z*w(264)
      w(264)= - 13.D0/3.D0 + w(264)
      w(264)=w(264)*w(246)
      w(326)=nf*z
      w(327)=w(326) + w(96)
      w(327)=LL*w(327)
      w(232)=8.D0/3.D0*w(327) + w(264) + 1.D0 + w(232)
      w(232)=LL*w(232)
      w(264)=191.D0 - 301.D0*z
      w(264)=w(264)*w(97)
      w(264)= - 221.D0 + w(264)
      w(327)= - 13.D0 - 292.D0/3.D0*z
      w(327)=z*w(327)
      w(327)= - 14.D0 + w(327)
      w(327)=nf*w(327)
      w(264)=1.D0/3.D0*w(264) + w(327)
      w(232)=1.D0/9.D0*w(264) + w(232)
      w(232)=w(232)*w(269)
      w(264)=w(101) - w(136)
      w(264)=w(264)*w(206)
      w(264)=w(264) - w(136)
      w(327)=nf + 2.D0
      w(264)= - w(264)*w(327)
      w(295)=2.D0/3.D0*w(295) + w(101)
      w(295)=nf*w(295)
      w(295)=2.D0/3.D0*w(320) + w(295)
      w(320)=w(246) + 1.D0
      w(328)= - w(320)*w(230)
      w(295)=1.D0/3.D0*w(295) + w(328)
      w(295)=w(295)*w(160)
      w(264)=w(295) + w(264)
      w(264)=w(264)*w(257)
      w(155)=w(155) + w(146)
      w(295)=8.D0/3.D0*w(25) - 4.D0/3.D0*w(63) + 16.D0/3.D0*w(16)
      w(295)=w(155)*w(295)
      w(328)=w(184) + w(263)
      w(332)=w(40) + w(37)
      w(332)= - 20.D0/3.D0*w(36) - 8.D0/3.D0*w(85) - 4.D0/3.D0*w(332)
      w(332)=w(328)*w(332)
      w(118)=32.D0/3.D0*w(130) + 16.D0/3.D0*w(118)
      w(118)=w(99)*w(118)
      w(130)=w(243) - w(98)
      w(130)=w(130)*w(206)
      w(130)=w(130) - w(136)
      w(130)=w(130)*w(8)
      w(243)=w(146)*z3
      w(130)=w(130) + w(243)
      w(243)=2.D0*w(7)
      w(130)=w(243)*w(130)*w(327)
      w(333)=w(101)*nf
      w(333)=w(333) + w(124)
      w(334)=w(30)*w(333)
      w(336)=1.D0 - w(294)
      w(336)=z*w(336)
      w(336)=w(336) + w(196)
      w(337)= - w(293) + w(132)
      w(337)=w(337)*w(246)
      w(336)=2.D0*w(336) + w(337)
      w(336)=w(83)*w(336)
      w(337)=w(131) + 2.D0/3.D0*w(263)
      w(338)=w(337)*w(160)
      w(339)=w(301) - w(245)
      w(339)= - w(339)*w(327)
      w(338)=1.D0/3.D0*w(339) + w(338)
      w(299)=w(338)*w(299)
      w(338)=z3 + 164.D0/9.D0
      w(339)=w(338)*w(158)
      w(339)=w(339) + 1447.D0/9.D0 + 14.D0*z3
      w(288)= - w(288)*w(160)
      w(288)=w(288) - 53.D0 + 13.D0*nf
      w(288)=w(288)*w(198)
      w(288)=w(288) - 22.D0 - w(246)
      w(288)=LL*w(288)
      w(288)=1.D0/3.D0*w(339) + w(288)
      w(339)=4.D0/9.D0*w(2)
      w(288)=w(288)*w(339)
      w(340)=w(215)*w(327)
      w(341)= - nf*w(230)
      w(340)=w(341) + w(340)
      w(340)=w(55)*w(340)
      w(341)=37.D0 + w(138)
      w(341)=w(341)*w(97)
      w(341)= - 13.D0 + w(341)
      w(342)=1.D0 - w(133)
      w(342)=nf*w(342)
      w(342)=w(342) + 1.D0 - w(258)
      w(342)=w(342)*w(160)
      w(224)=w(224) + 4.D0
      w(224)=w(224)*z
      w(224)=w(224) + 1.D0
      w(343)=nf*w(224)
      w(341)=w(342) + 1.D0/3.D0*w(341) + w(343)
      w(341)=w(62)*w(341)
      w(342)=14503.D0 - 41333.D0/3.D0*z
      w(342)=w(342)*w(97)
      w(342)= - 2945.D0 + w(342)
      w(168)=w(168)*w(164)
      w(188)= - w(131)*w(188)
      w(343)= - w(236)*w(178)
      w(344)=5083.D0/3.D0 + 2780.D0*z
      w(344)=z*w(344)
      w(344)=346.D0/3.D0 + w(344)
      w(343)=1.D0/27.D0*w(344) + w(343)
      w(343)=w(343)*w(122)
      w(344)=31.D0 + 76.D0*z
      w(344)=w(344)*w(122)
      w(345)= - 106.D0 - 1255.D0/3.D0*z
      w(345)=z*w(345)
      w(344)=w(344) - 17.D0 + w(345)
      w(344)=z3*w(344)
      w(345)=w(103)*w(172)
      w(346)= - w(123)*w(246)
      w(346)=w(346) - 1.D0 - w(277)
      w(346)=w(65)*w(346)
      w(300)=w(71)*w(300)
      w(118)=w(264) + 2.D0/3.D0*w(341) + w(187) + w(232) + 4.D0/3.D0*
     & w(340) + w(288) + w(299) + w(152) + w(151) + 2.D0*w(336) + 4.D0*
     & w(346) + 8.D0/3.D0*w(334) + w(130) + w(148) - 32.D0/3.D0*w(153)
     &  + w(345) + 2.D0/9.D0*w(344) + w(300) + w(343) + w(188) + 32.D0/
     & 9.D0*w(171) + 1.D0/81.D0*w(342) + w(168) + w(118) + w(332) + 
     & w(295)
      w(130)=4.D0*CA
      w(118)=w(118)*w(130)
      w(148)= - 14.D0 - z
      w(148)=w(148)*w(97)
      w(148)=17.D0 + w(148)
      w(148)=w(148)*w(164)
      w(151)=1643.D0 + 27824.D0/3.D0*z
      w(151)=z*w(151)
      w(151)=877.D0 + w(151)
      w(148)=1.D0/9.D0*w(151) + w(148)
      w(151)=w(106)*w(165)
      w(152)=w(73)*w(151)
      w(148)=1.D0/3.D0*w(148) + w(152)
      w(148)=w(148)*w(122)
      w(152)=37.D0 - w(107)
      w(152)=w(152)*w(97)
      w(152)= - 19.D0 + w(152)
      w(151)=1.D0/3.D0*w(152) - w(151)
      w(151)=w(151)*w(250)
      w(152)= - 31.D0 - 172.D0/3.D0*z
      w(152)=w(152)*w(145)
      w(152)= - 365.D0 + w(152)
      w(153)= - 133.D0/3.D0 + w(112)
      w(153)=z*w(153)
      w(153)= - 28.D0/3.D0 + w(153)
      w(153)=w(153)*w(165)
      w(152)=1.D0/3.D0*w(152) + w(153)
      w(152)=w(152)*w(252)
      w(153)=w(167)*w(34)
      w(168)=12557.D0 - 110579.D0/9.D0*z
      w(168)=w(168)*w(97)
      w(168)= - 1907.D0 + w(168)
      w(187)=109.D0 + 83.D0*z
      w(187)=w(187)*w(97)
      w(187)=1.D0 + w(187)
      w(187)=w(187)*w(282)
      w(188)= - z5*w(335)
      w(148)= - 16.D0/3.D0*w(153) + w(152) + w(151) + w(148) + 64.D0/9.D
     & 0*w(171) + w(187) + 1.D0/27.D0*w(168) + w(188)
      w(148)=w(148)*w(158)
      w(151)=113.D0 + 250.D0/3.D0*z
      w(151)=w(151)*w(110)
      w(152)= - 73.D0 - w(291)
      w(152)=w(152)*w(97)
      w(152)= - 305.D0 + w(152)
      w(152)=w(152)*w(122)
      w(151)=w(152) - 749.D0 + w(151)
      w(152)=1.D0/9.D0*nf
      w(151)=w(151)*w(152)
      w(168)=62.D0/3.D0*z
      w(171)= - 43.D0 + w(168)
      w(171)=w(171)*w(97)
      w(171)=53.D0 + w(171)
      w(187)=1.D0/3.D0 - w(97)
      w(187)=z*w(187)
      w(187)=4.D0/3.D0 + w(187)
      w(187)=w(187)*w(165)
      w(171)=1.D0/3.D0*w(171) + w(187)
      w(171)=nf*w(171)
      w(187)=1.D0/3.D0 - z
      w(187)=w(187)*w(272)
      w(187)=5.D0/3.D0 + w(187)
      w(187)=w(8)*w(187)
      w(188)= - 37.D0 + w(168)
      w(188)=z*w(188)
      w(188)=25.D0 + w(188)
      w(187)=1.D0/3.D0*w(188) + w(187)
      w(171)=2.D0*w(187) + w(171)
      w(171)=w(171)*w(108)
      w(187)=1.D0 - 287.D0/3.D0*z
      w(187)=w(187)*w(272)
      w(187)=559.D0 + w(187)
      w(188)=83.D0/9.D0 + w(272)
      w(188)=w(188)*w(97)
      w(188)=325.D0/9.D0 + w(188)
      w(188)=w(188)*w(122)
      w(151)=w(171) + w(151) + 1.D0/9.D0*w(187) + w(188)
      w(151)=w(151)*w(160)
      w(171)= - 12913.D0 + 10480.D0/3.D0*z
      w(171)=w(171)*w(97)
      w(171)=20789.D0 + w(171)
      w(187)= - w(96)*w(133)
      w(187)=23.D0 + w(187)
      w(187)=w(187)*w(190)
      w(171)=1.D0/9.D0*w(171) + w(187)
      w(187)=613.D0 - w(271)
      w(187)=w(187)*w(97)
      w(187)=1319.D0 + w(187)
      w(188)=w(106)*w(164)
      w(187)=1.D0/27.D0*w(187) + w(188)
      w(103)=w(187)*w(103)
      w(187)=w(318) - 3.D0 + 10.D0/3.D0*w(173)
      w(116)=w(187)*w(116)
      w(187)=w(71)*w(106)
      w(103)=w(116) - 96.D0*w(187) + 1.D0/3.D0*w(171) + w(103)
      w(103)=nf*w(103)
      w(116)= - 14185.D0 + 3446.D0/3.D0*z
      w(116)=w(116)*w(97)
      w(116)=28379.D0 + w(116)
      w(171)=2203.D0/3.D0 - w(291)
      w(171)=z*w(171)
      w(171)=4228.D0/3.D0 + w(171)
      w(171)=w(171)*w(219)
      w(187)=w(73)*w(330)
      w(116)=w(171) + 1.D0/9.D0*w(116) + 128.D0*w(187)
      w(171)= - w(337)*w(297)
      w(187)=z + 1.D0/3.D0
      w(188)=w(187)*w(97)
      w(188)= - 1.D0/3.D0 + w(188)
      w(170)=w(188)*w(170)
      w(188)=w(106)*nf
      w(190)=w(188) + w(106)
      w(219)=w(65)*w(190)
      w(103)=w(151) + w(171) - 80.D0*w(219) + w(103) + 1.D0/3.D0*w(116)
     &  + w(170)
      w(103)=LL*w(103)
      w(116)=4.D0*w(131)
      w(151)= - w(116) - w(263)
      w(151)=w(151)*w(198)
      w(170)= - w(120)*w(272)
      w(170)=5.D0 + w(170)
      w(171)=w(325)*w(246)
      w(151)=w(151) + w(171) + 1.D0/3.D0*w(170) - w(132)
      w(151)=w(151)*w(108)
      w(170)=101.D0 - w(298)
      w(170)=w(170)*w(97)
      w(170)=w(197) - 149.D0 + w(170)
      w(170)=1.D0/9.D0*w(170) + w(167)
      w(170)=w(8)*w(170)
      w(171)=77.D0 - 68.D0*z
      w(171)=z*w(171)
      w(171)= - 28.D0 + w(171)
      w(170)=2.D0/27.D0*w(171) + w(170)
      w(170)=nf*w(170)
      w(171)= - 26.D0 + w(138)
      w(171)=z*w(171)
      w(171)=4.D0 + w(171)
      w(171)=1.D0/3.D0*w(171) + w(167)
      w(171)=w(171)*w(165)
      w(219)=271.D0 - 244.D0*z
      w(219)=z*w(219)
      w(219)= - 56.D0 + w(219)
      w(171)=1.D0/9.D0*w(219) + w(171)
      w(151)=w(151) + 1.D0/3.D0*w(171) + w(170)
      w(151)=w(151)*w(160)
      w(170)=2.D0/9.D0*w(167)
      w(171)= - 2.D0 + w(294)
      w(171)=w(171)*w(97)
      w(171)=w(170) + 1.D0/3.D0 + w(171)
      w(171)=w(8)*w(171)
      w(219)=16.D0/3.D0*z
      w(232)=9.D0 - w(219)
      w(232)=z*w(232)
      w(232)= - 10.D0/3.D0 + w(232)
      w(171)=2.D0*w(232) + w(171)
      w(171)=w(8)*w(171)
      w(232)=w(105) - 17.D0
      w(232)=w(232)*z
      w(232)=w(232) - 1.D0
      w(232)=2.D0/3.D0*w(232)
      w(171)= - 40.D0/9.D0*w(204) + w(232) + w(171)
      w(250)=13.D0 - 47.D0/9.D0*z
      w(250)=w(250)*w(97)
      w(264)=16.D0/9.D0*w(6)
      w(170)=w(170) + w(264) - 53.D0/3.D0 + w(250)
      w(170)=w(8)*w(170)
      w(250)= - 57.D0 + 1556.D0/27.D0*z
      w(250)=z*w(250)
      w(250)= - w(289) + 8.D0/3.D0 + w(250)
      w(170)=2.D0*w(250) + w(170)
      w(170)=w(8)*w(170)
      w(170)= - 44.D0/9.D0*w(204) + w(232) + w(170)
      w(170)=nf*w(170)
      w(151)=w(151) + 2.D0*w(171) + w(170)
      w(170)=4.D0*w(32)
      w(151)=w(151)*w(170)
      w(171)= - 97.D0/9.D0 + w(110)
      w(171)=w(171)*w(97)
      w(217)= - 13.D0 + w(217)
      w(217)=w(217)*w(109)
      w(171)=w(217) - w(264) + 175.D0/9.D0 + w(171)
      w(171)=nf*w(171)
      w(217)=w(263) + w(131)
      w(232)= - w(217)*w(108)
      w(250)=w(277) - 1.D0
      w(264)=w(250)*w(97)
      w(271)= - 1.D0 - w(264)
      w(271)=w(8)*w(271)
      w(222)=4.D0 - w(222)
      w(222)=1.D0/9.D0*w(222) + w(271)
      w(171)=w(232) + 4.D0*w(222) + w(171)
      w(171)=LL*w(171)
      w(222)= - 7.D0 + w(110)
      w(107)=w(222)*w(107)
      w(222)= - w(331) - 2.D0 - w(293)
      w(222)=w(222)*w(122)
      w(107)=w(222) + 16.D0 + w(107)
      w(222)= - 74.D0 - 101.D0/3.D0*z
      w(222)=z*w(222)
      w(222)= - w(135) + 13.D0 + w(222)
      w(232)=2.D0/3.D0 + z
      w(232)=w(232)*w(97)
      w(232)= - 11.D0/3.D0 + w(232)
      w(232)=w(8)*w(232)
      w(222)=2.D0/3.D0*w(222) + w(232)
      w(222)=w(8)*w(222)
      w(232)=149.D0 - 1448.D0/9.D0*z
      w(232)=z*w(232)
      w(232)=w(266) - 2.D0 + w(232)
      w(222)=1.D0/3.D0*w(232) + w(222)
      w(222)=nf*w(222)
      w(107)=w(171) + 2.D0/3.D0*w(107) + w(222)
      w(171)=8.D0*w(77)
      w(107)=w(107)*w(171)
      w(222)= - 65.D0 + 74.D0*z
      w(222)=w(222)*w(110)
      w(232)=97.D0 + 145.D0*z
      w(232)=w(232)*w(97)
      w(232)=377.D0 + w(232)
      w(232)=w(232)*w(158)
      w(222)=w(232) + 931.D0 + w(222)
      w(232)=1.D0 + w(241)
      w(232)=w(232)*w(97)
      w(232)= - 3.D0 + w(232)
      w(232)=nf*w(232)
      w(190)= - w(190)*w(198)
      w(187)= - z*w(187)
      w(187)=2.D0/3.D0 + w(187)
      w(187)=w(190) + 4.D0*w(187) + w(232)
      w(187)=w(187)*w(186)
      w(187)=1.D0/9.D0*w(222) + w(187)
      w(187)=LL*w(187)
      w(190)=w(106)*z3
      w(219)=45.D0 - w(219)
      w(219)=z*w(219)
      w(219)=1.D0/3.D0*w(190) - 58.D0/3.D0 + w(219)
      w(222)=443.D0 - 3112.D0/3.D0*z
      w(222)=z*w(222)
      w(222)=70.D0 + w(222)
      w(190)=1.D0/3.D0*w(222) + 52.D0*w(190)
      w(190)=w(190)*w(246)
      w(187)=w(187) + 4.D0*w(219) + w(190)
      w(190)=4.D0*w(54)
      w(187)=w(187)*w(190)
      w(219)= - w(276)*w(214)
      w(219)=55.D0 + w(219)
      w(222)= - w(123)*w(157)
      w(222)=7.D0 + w(222)
      w(222)=w(222)*w(158)
      w(106)=w(188) - w(106)
      w(232)=w(106)*w(186)
      w(219)=w(232) + 1.D0/3.D0*w(219) + w(222)
      w(219)=w(219)*w(160)
      w(222)= - 3.D0 + 202.D0/9.D0*z
      w(222)=w(222)*w(97)
      w(222)=23.D0/3.D0 + w(222)
      w(222)=nf*w(222)
      w(232)=31.D0 + 20.D0/3.D0*z
      w(232)=z*w(232)
      w(232)= - 43.D0/3.D0 + w(232)
      w(219)=w(219) + 2.D0*w(232) + w(222)
      w(222)=4.D0*w(62)
      w(219)=w(219)*w(222)
      w(232)=w(85) + w(40)
      w(232)=64.D0/3.D0*w(36) - 32.D0/3.D0*w(38) + 16.D0/3.D0*w(232)
      w(232)=w(328)*w(232)
      w(241)=w(216)*w(272)
      w(241)= - 13.D0 + w(241)
      w(266)=w(97) - 11.D0/3.D0
      w(266)=w(266)*z
      w(271)=2.D0/3.D0 - w(266)
      w(271)=w(271)*w(158)
      w(241)=1.D0/3.D0*w(241) + w(271)
      w(241)=w(241)*w(296)
      w(271)=w(329) - 1.D0/3.D0
      w(271)=w(297)*w(271)*w(327)
      w(276)= - LL*w(337)
      w(282)=z - 2.D0/3.D0
      w(282)=w(282)*z
      w(282)=w(282) - 1.D0/3.D0
      w(282)=w(282)*w(327)
      w(276)=w(276) + w(282)
      w(276)=w(276)*w(303)
      w(282)= - 1.D0 - nf
      w(282)=w(282)*w(160)
      w(288)=1.D0 - nf
      w(282)=13.D0*w(288) + w(282)
      w(282)=LL*w(282)
      w(288)= - 50.D0 - 133.D0*nf
      w(282)=1.D0/3.D0*w(288) + w(282)
      w(198)=w(282)*w(198)
      w(282)=nf*w(338)
      w(198)=w(198) + 2.D0/3.D0*w(282) - 5.D0 + w(252)
      w(198)=w(2)*w(198)
      w(252)=3469.D0 - 4897.D0/3.D0*z
      w(252)=w(252)*w(277)
      w(277)=w(307)*w(97)
      w(277)=7.D0 + w(277)
      w(277)=w(73)*w(277)
      w(282)=3.D0 + w(305)
      w(282)=z*w(282)
      w(288)=w(126) - 1.D0
      w(289)=w(73)*w(288)
      w(282)= - 2.D0/3.D0*w(289) - 23.D0 + w(282)
      w(233)=w(282)*w(233)
      w(282)=w(71)*w(324)
      w(289)=5.D0/9.D0 - z
      w(289)=w(289)*w(272)
      w(289)= - 29.D0/9.D0 + w(289)
      w(289)=w(8)*w(289)
      w(291)=61.D0 - 170.D0/3.D0*z
      w(291)=z*w(291)
      w(291)= - 49.D0 + w(291)
      w(289)=1.D0/9.D0*w(291) + w(289)
      w(220)=w(289)*w(220)
      w(289)=w(68)*w(106)
      w(291)=w(70)*w(188)
      w(103)=64.D0*w(291) - 16.D0*w(289) + w(219) + w(107) + w(187) + 
     & 32.D0/9.D0*w(198) + w(276) + w(151) + w(103) + w(271) + w(241)
     &  + w(148) - 64.D0/3.D0*w(153) + w(220) + 32.D0/3.D0*w(282) + 
     & w(233) + 16.D0/3.D0*w(277) - 1335.D0 + w(252) + w(232)
      w(103)=CF*w(103)
      w(107)=w(155)*LL
      w(148)=w(101) + w(98)
      w(151)=w(327)*w(148)
      w(107)=w(107) - w(151)
      w(151)=w(21)*CA*w(107)
      w(153)= - w(131)*LL**3
      w(153)=w(204) + w(153)
      w(153)=TF*w(153)
      w(103)=32.D0/9.D0*w(153) + 32.D0/3.D0*w(151) + w(118) + w(103)
      w(118)=2.D0*TF
      w(103)=w(103)*w(118)
      w(151)= - 83.D0 - 44.D0*z
      w(151)=w(151)*w(267)
      w(153)= - w(120)*w(165)
      w(187)=w(99)*w(141)
      w(198)= - w(114)*w(317)
      w(151)=w(198) + w(187) + w(153) + 1.D0 + w(151)
      w(151)=w(151)*w(186)
      w(153)=w(319)*w(272)
      w(153)= - 85.D0 + w(153)
      w(153)=1.D0/3.D0*w(153) + w(234)
      w(153)=w(153)*w(243)
      w(187)= - 2860.D0 - 1333.D0*z
      w(187)=w(187)*w(267)
      w(187)=254.D0 + w(187)
      w(198)=30.D0*z
      w(219)= - 31.D0 + w(198)
      w(219)=z*w(219)
      w(219)=13.D0 + w(219)
      w(219)=w(219)*w(122)
      w(151)=w(151) + w(153) + 1.D0/3.D0*w(187) + w(219)
      w(151)=w(151)*w(160)
      w(153)=w(131)*w(186)
      w(187)=38.D0/3.D0 - w(203)
      w(187)=w(187)*w(97)
      w(187)=w(153) - 41.D0/3.D0 + w(187) + w(194)
      w(187)=w(187)*w(160)
      w(194)=32.D0 - w(200)
      w(194)=z*w(194)
      w(194)=w(281) + 7.D0 + w(194)
      w(194)=2.D0*w(194) - w(167)
      w(194)=w(8)*w(194)
      w(200)= - 3052.D0/3.D0 + 1037.D0*z
      w(200)=z*w(200)
      w(200)= - 242.D0/3.D0*w(6) + 116.D0/3.D0 + w(200)
      w(187)=w(187) + 1.D0/3.D0*w(200) + w(194)
      w(187)=w(187)*w(192)
      w(168)= - 43.D0 - w(168)
      w(168)=w(168)*w(97)
      w(194)=19.D0 + z
      w(194)=w(194)*w(97)
      w(194)=11.D0 + w(194)
      w(194)=w(8)*w(194)
      w(200)=13.D0 - w(133)
      w(200)=w(200)*w(97)
      w(200)= - 9.D0 + w(200)
      w(200)=w(200)*w(160)
      w(168)=w(200) + w(194) - w(202) - 11.D0 + w(168)
      w(168)=w(168)*w(211)
      w(194)=1.D0 - w(157)
      w(194)=w(194)*w(97)
      w(194)= - 5.D0 + w(194)
      w(194)=w(194)*w(290)
      w(200)= - w(310)*w(97)
      w(200)= - 13.D0 + w(200)
      w(200)=w(200)*w(160)
      w(219)= - 70.D0 - w(284)
      w(219)=z*w(219)
      w(219)= - 2.D0 + w(219)
      w(194)=w(200) + w(194) + 1.D0/3.D0*w(219) - w(182)
      w(139)=w(194)*w(139)
      w(194)=9.D0*w(101)
      w(200)=46.D0 + 67.D0*z
      w(200)=z*w(200)
      w(200)=2.D0 + w(200)
      w(200)=1.D0/3.D0*w(200) + w(182)
      w(200)=2.D0*w(200) + w(194)
      w(200)=w(8)*w(200)
      w(219)= - 127.D0 - 413.D0/3.D0*z
      w(219)=z*w(219)
      w(219)=4.D0 + w(219)
      w(200)=4.D0/3.D0*w(219) + w(200)
      w(200)=w(200)*w(243)
      w(219)= - 1.D0 + w(278)
      w(219)=w(219)*w(97)
      w(220)= - 11.D0 + w(253)
      w(220)=w(220)*w(97)
      w(220)= - 1.D0 + w(220)
      w(220)=w(220)*w(160)
      w(219)=w(220) + 47.D0/3.D0 + w(219)
      w(219)=w(219)*w(269)
      w(220)=64.D0*w(172)
      w(232)=w(131)*w(39)
      w(233)=48.D0*w(30) - 160.D0*w(11)
      w(233)=w(99)*w(233)
      w(241)= - 23.D0 - z
      w(241)=w(241)*w(97)
      w(241)= - 7.D0 + w(241)
      w(241)=w(241)*w(164)
      w(252)= - 932.D0/3.D0 - 1079.D0*z
      w(252)=w(252)*w(97)
      w(252)= - 160.D0/3.D0*w(6) - 1165.D0/3.D0 + w(252)
      w(252)=w(252)*w(206)
      w(271)=83.D0 + 278.D0*z
      w(271)=w(271)*w(208)
      w(147)=w(210)*w(147)
      w(159)=w(55)*w(159)
      w(210)= - 3403.D0 + 57548.D0/3.D0*z
      w(210)=z*w(210)
      w(276)=w(6)*w(9)
      w(277)= - 7.D0 - LL
      w(277)=LL*w(277)
      w(277)= - 68.D0/9.D0 + w(277)
      w(277)=w(2)*w(277)
      w(278)=w(4)*w(3)
      w(281)= - 5.D0 + w(272)
      w(281)=w(62)*w(281)
      w(139)=w(139) + 8.D0*w(281) + w(168) + w(219) - 160.D0/9.D0*
     & w(278) + 12.D0*w(159) + 32.D0/3.D0*w(277) + 32.D0*w(232) + 
     & w(187) + w(151) + w(147) + w(200) + w(220) + w(271) + w(252) + 
     & 64.D0/3.D0*w(276) + w(241) + 170.D0 + 1.D0/9.D0*w(210) + w(233)
      w(139)=w(139)*w(119)
      w(147)=w(146)*w(7)
      w(151)= - 1.D0 - w(315)
      w(151)=w(8)*w(151)
      w(159)=w(97) - 5.D0
      w(168)=w(159)*w(97)
      w(187)=1.D0 - w(168)
      w(187)=w(187)*w(108)
      w(200)=65.D0 - w(145)
      w(200)=z*w(200)
      w(200)= - 10.D0 + w(200)
      w(151)=w(187) - w(147) + 1.D0/3.D0*w(200) + w(151)
      w(151)=w(151)*w(160)
      w(187)= - 18.D0 - w(213)
      w(187)=w(187)*w(97)
      w(187)= - 44.D0*w(101) - 11.D0 + w(187)
      w(187)=w(7)*w(187)
      w(200)=62.D0 + 59.D0*z
      w(200)=w(200)*w(97)
      w(200)=3.D0 + w(200)
      w(200)=w(8)*w(200)
      w(210)=1019.D0 - 757.D0*z
      w(210)=z*w(210)
      w(210)= - 100.D0 + w(210)
      w(151)=w(151) + w(187) + 1.D0/9.D0*w(210) + w(200)
      w(151)=w(151)*w(270)
      w(112)=w(250)*w(112)
      w(112)= - w(175) - 32.D0*w(167) - w(202) + 15.D0 + w(112)
      w(112)=w(112)*w(160)
      w(138)= - 61.D0/3.D0 + w(138)
      w(138)=w(138)*w(97)
      w(138)=w(167) - w(202) - 23.D0/3.D0 + w(138)
      w(138)=w(8)*w(138)
      w(175)=1772.D0 - w(225)
      w(175)=w(175)*w(97)
      w(175)= - 293.D0 + w(175)
      w(112)=w(112) + w(138) + 1.D0/9.D0*w(175) + w(254)
      w(112)=w(112)*w(170)
      w(138)=253.D0/3.D0 - 21.D0*z
      w(138)=w(138)*w(97)
      w(175)=9.D0 - w(162)
      w(175)=w(175)*w(97)
      w(175)= - 21.D0 + w(175)
      w(175)=w(175)*w(122)
      w(187)= - 8.D0 + w(157)
      w(187)=w(187)*w(97)
      w(187)=3.D0 + w(187)
      w(187)=w(187)*w(270)
      w(138)=w(187) + w(175) + w(202) - 121.D0/3.D0 + w(138)
      w(138)=w(138)*w(211)
      w(175)= - 31.D0 + w(189)
      w(175)=w(175)*w(97)
      w(175)=1.D0 + w(175)
      w(187)= - 5.D0 - w(97)
      w(187)=w(187)*w(97)
      w(187)= - 9.D0 + w(187)
      w(187)=w(187)*w(270)
      w(175)=1.D0/3.D0*w(175) + w(187)
      w(175)=w(175)*w(269)
      w(187)=6.D0 + w(129)
      w(187)=z*w(187)
      w(187)=3.D0 + w(187)
      w(187)=w(187)*w(122)
      w(200)=13.D0 + w(97)
      w(200)=w(200)*w(97)
      w(200)=7.D0 + w(200)
      w(200)=w(200)*w(160)
      w(210)= - 9.D0 - w(213)
      w(210)=z*w(210)
      w(187)=w(200) + w(187) - 1.D0 + w(210)
      w(187)=w(41)*w(187)
      w(200)=1.D0 + w(162)
      w(200)=w(200)*w(97)
      w(200)=35.D0 + w(200)
      w(200)=w(200)*w(178)
      w(210)= - 2621.D0 - 1156.D0*z
      w(210)=z*w(210)
      w(210)= - 629.D0 + w(210)
      w(210)=w(210)*w(115)
      w(219)= - 151.D0 + 46.D0*z
      w(219)=w(219)*w(97)
      w(219)= - 89.D0 + w(219)
      w(208)=w(219)*w(208)
      w(219)=w(11)*w(99)
      w(225)=ln2*w(149)
      w(233)=7.D0*w(101)
      w(241)=17.D0 + w(239)
      w(241)=z*w(241)
      w(241)=1.D0 + w(241)
      w(241)=2.D0*w(241) - w(233)
      w(241)=w(8)*w(241)
      w(250)=13.D0 + w(213)
      w(250)=z*w(250)
      w(250)= - 4.D0 + w(250)
      w(241)=4.D0*w(250) + w(241)
      w(141)=w(241)*w(141)
      w(140)= - w(30)*w(140)
      w(241)=7.D0 + w(145)
      w(241)=w(241)*w(297)
      w(248)=w(248)*w(97)
      w(248)=w(248) + 5.D0
      w(237)= - w(248)*w(237)
      w(250)= - 14911.D0 - 27920.D0/3.D0*z
      w(250)=z*w(250)
      w(250)=467.D0/2.D0 + w(250)
      w(252)=3.D0 + w(110)
      w(252)=z*w(252)
      w(252)=3.D0 + w(252)
      w(252)=w(62)*w(252)
      w(112)=8.D0*w(187) + 16.D0*w(252) + w(138) + w(175) + w(237) - 80.
     & D0*w(2) - 256.D0*w(232) + w(112) + w(151) + w(241) + w(140) + 
     & w(141) + w(225) + 320.D0*w(219) - 192.D0*w(172) + w(208) + 
     & w(210) + 1.D0/9.D0*w(250) + w(200)
      w(112)=CF*w(112)
      w(112)=w(139) + w(112)
      w(112)=CA*w(112)
      w(138)= - 1.D0 + w(179)
      w(138)=w(138)*w(97)
      w(138)=3.D0 + w(138)
      w(138)=w(8)*w(138)
      w(108)=w(242)*w(108)
      w(139)=w(272) + 3.D0
      w(140)=z*w(139)
      w(108)=w(108) - w(147) + w(138) + 2.D0 + w(140)
      w(108)=w(108)*w(160)
      w(138)= - 8.D0 + w(203)
      w(138)=z*w(138)
      w(138)=2.D0 + w(138)
      w(138)=w(138)*w(165)
      w(140)=109.D0*z
      w(141)= - 10.D0 + w(140)
      w(141)=z*w(141)
      w(151)=w(7)*w(127)
      w(108)=w(108) - 16.D0*w(151) + w(138) + 26.D0 + w(141)
      w(108)=w(108)*w(270)
      w(138)=12.D0*w(167)
      w(141)= - 1.D0 + w(129)
      w(141)=w(141)*w(97)
      w(141)=w(153) + w(138) - 7.D0 + w(141)
      w(141)=w(141)*w(186)
      w(151)=5.D0 - w(179)
      w(151)=w(151)*w(110)
      w(151)= - 1.D0 + w(151)
      w(151)=w(8)*w(151)
      w(172)= - 36.D0 - w(213)
      w(172)=z*w(172)
      w(141)=w(141) + w(151) + 21.D0 + w(172)
      w(141)=w(141)*w(170)
      w(151)=w(55)*w(177)
      w(170)= - w(186) - w(265)
      w(170)=w(41)*w(170)
      w(151)=w(151) + w(170)
      w(170)=11.D0 + w(214)
      w(170)=w(170)*w(97)
      w(172)= - 19.D0 + w(212)
      w(172)=w(172)*w(97)
      w(172)=19.D0 + w(172)
      w(172)=w(172)*w(122)
      w(170)= - w(247) + w(172) + 23.D0 + w(170)
      w(170)=w(170)*w(211)
      w(172)= - 169.D0 - w(292)
      w(172)=w(172)*w(157)
      w(168)= - 5.D0 - w(168)
      w(168)=w(73)*w(168)
      w(175)=20.D0 - w(213)
      w(175)=w(175)*w(306)
      w(177)=w(240)*w(97)
      w(177)=1.D0 + w(177)
      w(177)=w(177)*w(268)
      w(187)=ln2*w(131)
      w(104)= - w(280)*w(104)
      w(149)=w(39)*w(149)
      w(200)=24.D0*LL
      w(203)=w(242)*w(200)
      w(208)=35.D0 - w(253)
      w(208)=z*w(208)
      w(203)=w(203) + 4.D0 + w(208)
      w(203)=w(203)*w(190)
      w(208)=w(181)*w(97)
      w(208)= - 1.D0 + w(208)
      w(208)=w(208)*w(222)
      w(210)=w(83)*w(173)
      w(104)=w(208) + w(170) + w(203) + w(149) + w(141) + w(108) + 96.D0
     & *w(210) + w(104) - 384.D0*w(187) + w(220) + w(177) + w(175) + 40.
     & D0*w(168) + 343.D0/2.D0 + w(172) + 16.D0*w(151)
      w(104)=w(113)*w(104)
      w(108)= - w(215)*w(158)
      w(141)=5.D0*w(101)
      w(96)=w(96)*w(272)
      w(96)=w(108) + w(96) - w(141)
      w(96)=w(7)*w(96)
      w(108)= - w(120)*w(129)
      w(108)= - 1.D0 + w(108)
      w(120)=w(128) - 5.D0
      w(149)= - nf*w(120)
      w(108)=10.D0*w(108) + w(149)
      w(149)=w(217)*w(160)
      w(108)=1.D0/3.D0*w(108) + w(149)
      w(108)=w(32)*w(108)
      w(149)=w(158) + 5.D0
      w(151)=w(41)*w(99)*w(149)
      w(96)=w(151) + w(96) + w(108)
      w(108)= - w(155)*w(243)
      w(151)=4.D0/3.D0 - z
      w(151)=w(151)*w(97)
      w(168)= - 10.D0/3.D0 - z
      w(168)=w(168)*w(326)
      w(108)=w(108) + w(151) + w(168)
      w(151)=nf*w(267)
      w(151)=z + w(151)
      w(151)=w(151)*w(186)
      w(108)=1.D0/3.D0*w(108) + w(151)
      w(108)=w(108)*w(186)
      w(151)=5.D0 - 56.D0/9.D0*z
      w(151)=w(151)*w(214)
      w(168)=31.D0 + w(218)
      w(168)=w(168)*w(272)
      w(168)=41.D0 + w(168)
      w(115)=w(168)*w(115)
      w(168)= - 8.D0 + w(298)
      w(168)=w(168)*w(97)
      w(168)= - 7.D0 + w(168)
      w(168)=w(168)*w(122)
      w(170)= - 50.D0 + w(239)
      w(170)=z*w(170)
      w(168)=w(168) + 2.D0 + w(170)
      w(152)=w(168)*w(152)
      w(168)=32.D0 - nf
      w(168)=w(168)*w(339)
      w(170)= - nf*w(236)
      w(170)=w(170) - 7.D0/3.D0 - w(97)
      w(170)=w(170)*w(190)
      w(96)=w(170) + w(168) + w(108) + w(152) + w(115) - 43.D0/9.D0 + 
     & w(151) + 4.D0/3.D0*w(96)
      w(96)=w(96)*w(119)
      w(108)=941.D0 - 3640.D0/3.D0*z
      w(108)=w(108)*w(97)
      w(108)=185.D0 + w(108)
      w(115)=83.D0/3.D0 + w(189)
      w(115)=w(115)*w(97)
      w(115)=193.D0/3.D0 + w(115)
      w(115)=w(115)*w(165)
      w(108)=1.D0/3.D0*w(108) + w(115)
      w(108)=nf*w(108)
      w(115)=553.D0 - 152.D0/3.D0*z
      w(115)=w(115)*w(110)
      w(115)= - 1489.D0 + w(115)
      w(151)= - 17.D0/3.D0 + z
      w(151)=w(151)*w(272)
      w(151)= - 559.D0/3.D0 + w(151)
      w(151)=w(151)*w(122)
      w(108)=w(108) + 1.D0/3.D0*w(115) + w(151)
      w(115)=1.D0/3.D0*w(120) + w(167)
      w(115)=nf*w(115)
      w(120)=w(323) - 1.D0
      w(115)=w(115) + 5.D0/3.D0*w(120) + w(245)
      w(115)=w(32)*w(115)
      w(151)=w(131)*w(320)
      w(152)=w(151)*LL
      w(168)=z*w(302)
      w(168)= - 8.D0 + w(168)
      w(170)= - 13.D0 - w(223)
      w(170)=nf*w(170)
      w(168)=2.D0*w(168) + w(170)
      w(168)=1.D0/9.D0*w(168) + w(152)
      w(168)=w(168)*w(247)
      w(170)=w(181)*w(272)
      w(172)=5.D0 - w(258)
      w(172)=z*w(172)
      w(172)=11.D0 + w(172)
      w(172)=w(172)*w(158)
      w(170)=w(172) - 29.D0 + w(170)
      w(170)=w(54)*w(170)
      w(172)=3.D0*w(131)
      w(175)=w(172) + 5.D0/3.D0*w(263)
      w(171)=w(175)*w(171)
      w(106)=w(62)*w(106)
      w(175)= - 16.D0 + 23.D0*nf
      w(175)=w(2)*w(175)
      w(106)= - 24.D0*w(106) + w(171) + 4.D0/3.D0*w(170) + 16.D0/27.D0*
     & w(175) + 8.D0/3.D0*w(115) + 1.D0/3.D0*w(108) + w(168)
      w(106)=CF*w(106)
      w(96)=w(106) + w(96)
      w(96)=w(96)*w(118)
      w(106)=w(131)*w(32)
      w(108)=37.D0 - 142.D0*z
      w(108)=w(108)*w(97)
      w(108)= - 47.D0 + w(108)
      w(115)=1.D0/5.D0*w(8)
      w(108)=w(108)*w(115)
      w(168)=47.D0 - 66.D0*z
      w(168)=w(168)*w(97)
      w(168)= - 63.D0 + w(168)
      w(168)=LL*w(168)
      w(170)=23.D0 - 54.D0*z
      w(170)=z*w(170)
      w(108)= - 132.D0/5.D0*w(106) + 1.D0/5.D0*w(168) + w(147) + w(108)
     &  - 23.D0/5.D0 + w(170)
      w(168)=2.D0*w(113)
      w(108)=w(108)*w(168)
      w(170)=w(99)*w(7)
      w(171)= - 2312.D0 + 1919.D0*z
      w(171)=z*w(171)
      w(171)=5.D0 + 1.D0/5.D0*w(171)
      w(175)= - 131.D0 + w(229)
      w(175)=z*w(175)
      w(175)= - 7.D0 + w(175)
      w(175)=w(8)*w(175)
      w(133)=19.D0/5.D0 + w(133)
      w(133)=LL*w(133)
      w(133)=247.D0/5.D0*w(106) + 3.D0*w(133) + 37.D0/5.D0*w(170) + 1.D0
     & /3.D0*w(171) + 2.D0/5.D0*w(175)
      w(171)=2.D0*CF
      w(133)=w(133)*w(171)
      w(142)= - 19.D0 + w(142)
      w(115)=w(142)*w(115)
      w(142)=1673.D0 + 86.D0*z
      w(142)=z*w(142)
      w(115)=w(115) - 6.D0 + 1.D0/15.D0*w(142)
      w(115)=2.D0*w(115) - 47.D0/5.D0*w(170)
      w(142)= - 879.D0 - w(313)
      w(142)=z*w(142)
      w(142)= - 49.D0 + 2.D0/5.D0*w(142)
      w(142)=LL*w(142)
      w(106)=304.D0/15.D0*w(2) - 46.D0*w(106) + 2.D0*w(115) + w(142)
      w(106)=CA*w(106)
      w(106)=w(133) + w(106)
      w(106)=CA*w(106)
      w(115)=61.D0 - w(140)
      w(115)=w(115)*w(97)
      w(133)= - 29.D0 - z
      w(133)=w(133)*w(97)
      w(133)=11.D0 + w(133)
      w(133)=nf*w(133)
      w(115)=4.D0*w(133) - 61.D0 + w(115)
      w(133)=w(188)*w(200)
      w(115)=1.D0/3.D0*w(115) + w(133)
      w(115)=CF*w(115)
      w(133)= - 27.D0*z - 38.D0/3.D0*w(326)
      w(133)=w(133)*w(119)
      w(115)=w(115) + w(133)
      w(115)=TF*w(115)
      w(106)=4.D0/5.D0*w(115) + w(108) + w(106)
      w(106)=z2*w(106)
      w(108)=w(99)*CF
      w(115)=w(99)*CA
      w(133)=3.D0*w(108) - w(115)
      w(133)=CA*w(133)
      w(142)=w(146)*w(113)
      w(133)= - w(142) + w(133)
      w(170)=16.D0*w(21)
      w(133)=w(133)*w(170)
      w(96)=4.D0*w(106) + w(96) + w(133) + w(112) + w(104)
      w(96)=z2*w(96)
      w(104)=w(125)*w(110)
      w(104)=w(104) + 1.D0
      w(106)=w(238) - w(104)
      w(106)=w(8)*w(106)
      w(112)=w(157) + 10.D0
      w(112)=w(112)*w(97)
      w(112)=w(112) + 17.D0
      w(125)=z*w(216)
      w(125)=3.D0 + w(125)
      w(125)=w(125)*w(165)
      w(125)=w(125) - w(112)
      w(125)=LL*w(125)
      w(133)=w(256) + 4.D0
      w(175)=2.D0*w(133)
      w(106)=w(125) + w(175) + w(106)
      w(106)=CF*w(106)
      w(125)=41.D0/3.D0 + w(258)
      w(125)=w(125)*w(97)
      w(177)=w(8)*w(207)
      w(125)=w(177) + w(125) + w(273)
      w(125)=w(8)*w(125)
      w(177)=w(99)*w(186)
      w(181)=w(213) + 7.D0
      w(181)=w(181)*w(214)
      w(181)=w(181) + w(197) - 37.D0
      w(181)=1.D0/3.D0*w(181)
      w(177)=w(177) + w(181)
      w(121)= - 3.D0 - w(121)
      w(121)=w(121)*w(165)
      w(121)=w(121) - w(177)
      w(121)=LL*w(121)
      w(187)=91.D0 + 305.D0/3.D0*z
      w(187)=w(187)*z
      w(187)=w(187) + w(283)
      w(188)=2.D0/3.D0*w(187)
      w(121)=w(121) + w(188) + w(125)
      w(121)=CA*w(121)
      w(106)=w(106) + w(121)
      w(106)=CA*w(106)
      w(121)=TF*CA
      w(125)=4.D0/3.D0*w(121)
      w(107)=w(107)*w(125)
      w(190)=w(8)*w(275)
      w(190)=w(190) + w(127)
      w(190)=4.D0*w(190) + w(230)
      w(197)=w(113)*w(186)
      w(190)=w(190)*w(197)
      w(200)=CF*w(248)
      w(203)=w(129) - 3.D0
      w(208)= - w(203)*w(97)
      w(208)= - 9.D0 + w(208)
      w(208)=CA*w(208)
      w(200)=w(200) + w(208)
      w(200)=CA*w(200)
      w(142)= - w(142) + w(200)
      w(142)=z2*w(142)
      w(106)=w(142) + w(107) + w(190) + w(106)
      w(106)=w(42)*w(106)
      w(107)=7.D0 - w(253)
      w(107)=w(107)*w(97)
      w(107)= - 7.D0 + w(107)
      w(107)=w(107)*w(113)
      w(142)= - 19.D0 + w(105)
      w(142)=w(142)*w(97)
      w(142)= - 1.D0 + w(142)
      w(142)=CF*w(142)
      w(190)=CA*z
      w(142)=w(142) - 108.D0*w(190)
      w(142)=CA*w(142)
      w(107)=w(107) + w(142)
      w(107)=w(70)*w(107)
      w(142)= - CF*w(322)
      w(190)=w(304)*w(130)
      w(142)=w(142) + w(190)
      w(142)=CA*w(142)
      w(190)= - w(113)*w(116)
      w(142)=w(190) + w(142)
      w(142)=w(82)*w(142)
      w(190)= - CF*w(120)
      w(130)= - w(114)*w(130)
      w(130)=w(190) + w(130)
      w(130)=w(130)*w(119)
      w(190)= - 3.D0 - w(315)
      w(190)=w(190)*w(113)
      w(130)=w(190) + w(130)
      w(130)=w(76)*w(130)
      w(107)=w(130) + w(107) + w(142)
      w(130)=3.D0 + w(226)
      w(130)=w(130)*w(122)
      w(130)= - w(231) + w(130) + w(227)
      w(130)=CF*w(130)
      w(142)= - 3.D0 - w(166)
      w(142)=w(142)*w(165)
      w(142)=w(231) + w(142) - w(285)
      w(142)=CA*w(142)
      w(130)=w(130) + w(142)
      w(130)=CA*w(130)
      w(142)=w(125)*w(155)
      w(130)=w(130) + w(142)
      w(130)=w(43)*w(130)
      w(155)=w(113)*w(122)*w(111)
      w(166)=CA*w(167)*LL*w(174)
      w(155)=w(155) + 3.D0*w(166)
      w(155)=w(35)*w(155)
      w(166)=w(69)*CA**2
      w(130)=w(130) + w(166) + w(155)
      w(104)=w(101) - w(104)
      w(104)=w(8)*w(104)
      w(155)=w(234) - w(112)
      w(155)=LL*w(155)
      w(104)=w(155) + w(175) + w(104)
      w(104)=CF*w(104)
      w(155)=41.D0/3.D0 + w(214)
      w(155)=w(155)*w(97)
      w(155)= - w(101) + w(155) + w(273)
      w(155)=w(8)*w(155)
      w(166)=w(228) - w(177)
      w(166)=LL*w(166)
      w(155)=w(166) + w(188) + w(155)
      w(155)=CA*w(155)
      w(104)=w(104) + w(155)
      w(104)=CA*w(104)
      w(127)=w(101) + 4.D0*w(127)
      w(155)=w(230) + w(127)
      w(155)=w(155)*w(197)
      w(104)=w(155) + w(104)
      w(104)=w(104)*w(170)
      w(155)=19.D0 + w(253)
      w(155)=w(155)*w(179)
      w(166)= - 1.D0 + w(329)
      w(166)=w(166)*w(193)
      w(170)=w(235)*w(97)
      w(170)= - 7.D0 + w(170)
      w(170)=w(170)*w(186)
      w(155)=w(170) + w(166) + w(202) - 13.D0 + w(155)
      w(155)=CF*w(155)
      w(166)= - 197.D0/3.D0 - w(105)
      w(166)=w(166)*w(97)
      w(156)=5.D0 + w(156)
      w(156)=w(156)*w(165)
      w(170)=7.D0 + z
      w(170)=w(170)*w(110)
      w(170)=9.D0 + w(170)
      w(170)=w(170)*w(186)
      w(156)=w(170) + w(156) - w(205) - 37.D0/3.D0 + w(166)
      w(156)=CA*w(156)
      w(155)=w(155) + w(156)
      w(155)=CA*w(155)
      w(156)= - w(288)*w(122)
      w(166)=7.D0 + w(97)
      w(166)=z*w(166)
      w(156)=w(156) + 5.D0 + w(166)
      w(166)=1.D0 - w(253)
      w(166)=w(166)*w(97)
      w(166)= - 1.D0 + w(166)
      w(166)=LL*w(166)
      w(156)=2.D0*w(156) + w(166)
      w(156)=w(156)*w(168)
      w(166)= - w(328)*w(125)
      w(155)=w(166) + w(156) + w(155)
      w(155)=w(75)*w(155)
      w(156)= - CF*w(116)
      w(166)=CA*w(172)
      w(156)=w(156) + w(166)
      w(156)=CA*w(156)
      w(166)=w(131)*w(113)
      w(156)=w(166) + w(156)
      w(156)=w(90)*w(156)
      w(96)=16.D0*w(106) + w(96) + 8.D0*w(155) + 64.D0*w(156) + w(103)
     &  + w(104) + w(102) + w(100) + 32.D0*w(130) + 16.D0*w(107)
      w(96)=TF*w(96)
      w(100)=w(180)*w(110)
      w(100)=w(251) + 1.D0 + w(100)
      w(100)=w(8)*w(100)
      w(100)= - 4.D0*w(133) + w(100)
      w(100)=w(8)*w(100)
      w(102)= - w(194) + w(112)
      w(102)=LL*w(102)*w(122)
      w(103)=40.D0*w(262)
      w(100)=w(102) + w(100) - w(103)
      w(100)=CF*w(100)
      w(102)= - 53.D0/3.D0 - w(105)
      w(102)=w(102)*w(97)
      w(102)= - w(251) + w(102) - w(273)
      w(102)=w(8)*w(102)
      w(102)= - 4.D0/3.D0*w(187) + w(102)
      w(102)=w(8)*w(102)
      w(104)=w(181) + w(101)
      w(104)=w(8)*w(104)
      w(106)=w(165)*w(230)
      w(104)=w(104) + w(106)
      w(104)=w(104)*w(160)
      w(102)=w(104) + w(102) + w(103)
      w(102)=CA*w(102)
      w(100)=w(100) + w(102)
      w(100)=CA*w(100)
      w(102)=w(101) + w(136)
      w(103)=w(122) + w(309)
      w(103)=w(102)*w(103)
      w(104)= - w(333)*w(160)
      w(103)=w(104) + w(103)
      w(103)=w(103)*w(125)
      w(104)= - w(8)*w(127)
      w(101)= - LL*w(101)
      w(101)=w(104) + w(101)
      w(101)=w(101)*w(113)*w(270)
      w(100)=w(103) + w(101) + w(100)
      w(100)=TF*w(100)
      w(101)=26.D0*w(230) - w(233) - w(183)
      w(101)=CF*w(101)
      w(103)=59.D0/3.D0 + w(145)
      w(103)=w(103)*w(97)
      w(103)= - w(231) + w(141) + 5.D0/3.D0 + w(103)
      w(103)=CA*w(103)
      w(101)=w(101) + w(103)
      w(101)=CA*w(101)
      w(103)=w(113)*w(124)
      w(101)= - w(142) + w(103) + w(101)
      w(103)=w(118)*z2
      w(101)=w(101)*w(103)
      w(100)=w(100) + w(101)
      w(100)=w(10)*w(100)
      w(101)= - w(143) + w(279)
      w(101)=5.D0*w(101) + w(132)
      w(101)=w(101)*w(165)
      w(104)=8.D0/3.D0*w(111)
      w(106)=59.D0/3.D0 - w(145)
      w(106)=w(106)*w(97)
      w(106)= - w(104) + w(106) + w(260)
      w(106)=w(106)*w(160)
      w(107)=2032.D0 - 2023.D0*z
      w(107)=z*w(107)
      w(107)=247.D0 + w(107)
      w(107)=1.D0/3.D0*w(107) + w(199)
      w(101)=w(106) + 1.D0/3.D0*w(107) + w(101)
      w(101)=LL*w(101)
      w(106)= - 71.D0 + 95.D0*z
      w(106)=w(106)*w(97)
      w(106)= - 13.D0 + w(106)
      w(106)=4.D0/3.D0*w(167) + 1.D0/3.D0*w(106) - w(205)
      w(106)=w(8)*w(106)
      w(107)=229.D0 - 241.D0*z
      w(107)=z*w(107)
      w(107)=20.D0*w(6) - 11.D0 + w(107)
      w(106)=4.D0/3.D0*w(107) + w(106)
      w(106)=w(8)*w(106)
      w(107)= - 259.D0 + 884.D0/3.D0*z
      w(107)=w(107)*w(110)
      w(107)= - 134.D0/3.D0*w(6) + 13.D0 + w(107)
      w(101)=w(101) - 104.D0/3.D0*w(204) + 1.D0/9.D0*w(107) + w(106)
      w(101)=CA*w(101)
      w(106)= - 89.D0/3.D0 + w(198)
      w(106)=z*w(106)
      w(106)=w(104) + w(169) - w(143) + 4.D0/3.D0 + w(106)
      w(106)=w(106)*w(160)
      w(107)=33.D0 - 97.D0/3.D0*z
      w(107)=z*w(107)
      w(107)= - w(185) + w(202) - 2.D0 + w(107)
      w(107)=w(107)*w(122)
      w(112)= - 2222.D0 + 2077.D0*z
      w(112)=z*w(112)
      w(112)= - 56.D0*w(6) + 121.D0 + w(112)
      w(106)=w(106) + 1.D0/9.D0*w(112) + w(107)
      w(106)=w(106)*w(160)
      w(107)=26.D0 - 83.D0/3.D0*z
      w(107)=w(107)*w(97)
      w(107)= - w(150) + w(202) + 3.D0 + w(107)
      w(107)=w(8)*w(107)
      w(112)= - 57.D0 + 445.D0/9.D0*z
      w(112)=z*w(112)
      w(112)= - 40.D0/9.D0*w(6) + 7.D0 + w(112)
      w(107)=4.D0*w(112) + w(107)
      w(107)=w(8)*w(107)
      w(112)=215.D0 - 202.D0*z
      w(112)=z*w(112)
      w(112)=17.D0/3.D0*w(6) - 151.D0/3.D0 + w(112)
      w(106)=w(106) + 160.D0/3.D0*w(204) + 2.D0*w(112) + w(107)
      w(106)=CF*w(106)
      w(101)=w(106) + w(101)
      w(101)=CA*w(101)
      w(106)=w(323) - w(167)
      w(107)=w(106)*w(109)
      w(106)=w(106)*w(206)
      w(106)=w(106) - 1.D0/3.D0 - w(264)
      w(106)=nf*w(106)
      w(109)= - 13.D0 + w(126)
      w(109)=1.D0/3.D0*w(109) - w(167)
      w(112)=16.D0 - w(162)
      w(112)=z*w(112)
      w(112)= - 23.D0 + w(112)
      w(112)=w(112)*w(246)
      w(109)=2.D0*w(109) + w(112)
      w(112)=w(151)*w(160)
      w(109)=1.D0/3.D0*w(109) + w(112)
      w(109)=LL*w(109)
      w(106)=w(109) + w(106) + w(107) - 1.D0 - w(266)
      w(106)=CA*w(106)
      w(107)= - z*w(123)
      w(107)=8.D0 + w(107)
      w(109)=13.D0 + w(321)
      w(109)=nf*w(109)
      w(107)=2.D0*w(107) + w(109)
      w(107)=1.D0/9.D0*w(107) - w(152)
      w(107)=w(107)*w(160)
      w(109)=w(286)*z
      w(109)=w(109) + 4.D0/3.D0
      w(109)=w(109)*w(327)
      w(107)=w(107) + w(109)
      w(107)=CF*w(107)
      w(106)=w(107) + w(106)
      w(106)=TF*w(106)
      w(107)=5.D0 - w(110)
      w(107)=w(107)*w(110)
      w(104)= - w(104) - w(138) - 1.D0 + w(107)
      w(104)=LL*w(104)
      w(107)= - w(185) + 7.D0 - w(176)
      w(107)=w(107)*w(122)
      w(109)=128.D0 - w(140)
      w(109)=z*w(109)
      w(104)=w(104) + w(107) - 29.D0 + w(109)
      w(104)=w(104)*w(160)
      w(107)=w(129) - 6.D0
      w(107)=w(107)*w(97)
      w(109)=w(196) + 1.D0 + w(107)
      w(109)=w(8)*w(109)
      w(112)= - z*w(159)
      w(112)= - 3.D0 + w(112)
      w(109)=4.D0*w(112) + w(109)
      w(109)=w(109)*w(122)
      w(112)= - 137.D0 + 152.D0*z
      w(112)=z*w(112)
      w(104)=w(104) - 56.D0/3.D0*w(204) + w(109) + 36.D0 + w(112)
      w(104)=w(104)*w(113)
      w(101)=4.D0*w(106) + w(104) + w(101)
      w(101)=TF*w(101)
      w(104)=45.D0 - w(209)
      w(104)=z*w(104)
      w(104)= - w(312) - w(137) + w(182) - 6.D0 + w(104)
      w(104)=w(104)*w(171)
      w(106)= - 47.D0/3.D0 + w(154)
      w(97)=w(106)*w(97)
      w(106)=w(131)*w(270)
      w(97)=w(106) - w(132) + w(97) - w(261)
      w(97)=CA*w(97)
      w(97)=w(104) + w(97)
      w(97)=CA*w(97)
      w(104)=w(201)*w(110)
      w(104)=w(134) + 18.D0*w(167) + 7.D0 + w(104)
      w(104)=w(104)*w(113)
      w(106)=13.D0/3.D0 + w(158)
      w(106)=CF*w(106)
      w(109)= - CA*w(149)
      w(106)=w(106) + 1.D0/3.D0*w(109)
      w(106)=w(118)*w(131)*w(106)
      w(97)=w(106) + w(104) + w(97)
      w(97)=w(97)*w(103)
      w(97)=w(101) + w(97)
      w(97)=w(33)*w(97)
      w(101)=w(88)*w(121)*w(131)*w(174)
      w(96)=64.D0*w(101) + 4.D0*w(97) + 8.D0*w(100) + w(96)
      w(96)=as*w(96)
      w(97)=w(323) + w(111)
      w(97)=w(97)*w(160)
      w(100)=w(135) + 3.D0
      w(101)= - z*w(195)
      w(101)=w(101) + w(100)
      w(101)=w(8)*w(101)
      w(97)=w(97) + w(101) + w(120)
      w(97)=w(32)*w(97)
      w(101)=w(102)*w(144)
      w(102)= - 2.D0 - 11.D0/3.D0*z
      w(102)=w(102)*w(110)
      w(102)= - 1.D0 + w(102)
      w(102)=w(8)*w(102)
      w(104)=w(7)*w(124)
      w(106)= - 25.D0 + 218.D0/9.D0*z
      w(106)=z*w(106)
      w(102)=w(104) + w(102) + 2.D0 + w(106)
      w(104)=LL*w(287)
      w(102)=2.D0*w(102) + w(104)
      w(102)=LL*w(102)
      w(104)= - 4.D0 + 17.D0/3.D0*z
      w(104)=w(104)*w(110)
      w(100)=w(255) + w(104) - w(100)
      w(100)=w(77)*w(100)
      w(104)= - w(163) - w(148)
      w(104)=w(104)*w(257)
      w(106)=w(39) + w(83)
      w(109)=w(184)*w(106)
      w(112)= - w(114)*w(178)
      w(113)= - 43.D0 - 400.D0/3.D0*z
      w(113)=z*w(113)
      w(113)= - 14.D0 + w(113)
      w(113)=w(113)*w(206)
      w(114)=w(311)*w(244)
      w(120)= - 10.D0/3.D0 - LL
      w(120)=w(2)*w(120)*w(317)
      w(122)=w(55)*w(146)
      w(124)=w(123)*w(186)
      w(124)=w(124) + w(224)
      w(124)=w(54)*w(124)
      w(123)= - w(123)*w(274)
      w(125)= - 157.D0 + 1588.D0/9.D0*z
      w(125)=z*w(125)
      w(125)= - 112.D0/9.D0*w(1) - 1.D0 + w(125)
      w(97)=w(104) + w(123) + w(100) + w(124) + w(122) + w(120) + w(97)
     &  + w(102) + w(101) + w(114) + w(113) + 1.D0/3.D0*w(125) + w(112)
     &  + w(109)
      w(97)=w(97)*w(119)
      w(100)=w(107) + w(161)
      w(100)=w(8)*w(100)
      w(101)= - w(128) - w(167)
      w(101)=2.D0*w(101) - w(111)
      w(101)=w(101)*w(160)
      w(102)=13.D0 - w(253)
      w(102)=z*w(102)
      w(100)=w(101) + w(102) + w(100)
      w(100)=w(100)*w(192)
      w(101)= - 3.D0 - w(221)
      w(101)=w(8)*w(101)
      w(102)=29.D0 - w(105)
      w(102)=z*w(102)
      w(101)=w(101) - 14.D0 + w(102)
      w(102)= - w(316) + w(117)
      w(102)=LL*w(102)
      w(101)=2.D0*w(101) + w(102)
      w(101)=LL*w(101)
      w(102)=w(203)*w(110)
      w(104)= - w(242)*w(186)
      w(102)=w(104) - 1.D0 + w(102)
      w(102)=w(54)*w(102)
      w(104)= - w(314)*w(272)
      w(105)= - w(330)*w(165)
      w(104)=w(105) + 1.D0 + w(104)
      w(104)=w(104)*w(259)
      w(105)= - w(116)*w(106)
      w(106)=w(207)*w(164)
      w(107)= - w(139)*w(157)
      w(107)= - 8.D0 + w(107)
      w(107)=w(8)*w(107)
      w(109)=w(191)*w(244)
      w(110)=w(242)*w(274)
      w(111)= - 41.D0 + w(189)
      w(111)=z*w(111)
      w(100)=w(110) + w(104) + w(102) + w(100) + w(101) + w(109) + 
     & w(107) + w(106) + 13.D0 + w(111) + w(105)
      w(100)=CF*w(100)
      w(101)=w(21)*w(115)
      w(102)=TF*w(131)*LL**2
      w(97)= - 8.D0/3.D0*w(102) + 8.D0*w(101) + w(100) + w(97)
      w(97)=TF*w(97)
      w(100)=w(131)*w(160)
      w(100)=w(100) - 1.D0
      w(101)=w(249)*w(157)
      w(101)=w(101) + w(100)
      w(101)=w(101)*w(171)
      w(102)=4.D0 - z
      w(102)=z*w(102)
      w(104)=LL*w(272)
      w(102)=w(104) + w(102) - w(147)
      w(102)=CA*w(102)
      w(101)=w(101) + w(102)
      w(101)=w(101)*w(103)
      w(102)= - w(10)*w(234)
      w(99)=w(42)*w(99)
      w(99)=8.D0*w(99) + w(102)
      w(99)=w(121)*w(99)
      w(102)=w(153) - w(132) + w(301)
      w(102)=CA*w(102)
      w(100)= - w(308) - w(100)
      w(100)=w(100)*w(171)
      w(100)=w(100) + w(102)
      w(100)=w(33)*w(100)*w(118)
      w(97)=w(100) + w(97) + w(101) + w(99)
      w(96)=2.D0*w(97) + w(96)
      w(96)=as*w(96)
      w(97)= - TF*w(153)
      w(96)=w(97) + w(96)
      w(96)=as*w(96)
      w(97)=w(108) - w(115)
      w(97)=w(26)*w(97)
      w(99)= - w(173)*w(171)
      w(100)= - CA*w(207)
      w(99)=w(99) + w(100)
      w(99)=w(51)*w(99)
      w(98)=CA*w(98)
      w(98)= - w(108) + w(98)
      w(98)=w(18)*w(98)
      w(97)=384.D0*w(98) + 64.D0*w(99) + 160.D0*w(97)
      w(97)=w(121)*w(97)*as**3
*
      QG = w(96) + w(97)
*
      RETURN
      END
      REAL*8 FUNCTION AGQ(z,nf,as,LL)
*     ------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     complete OME AgQ  
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(89)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
      REAL*8 B4,z4
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z4 = 1.0823232337111381915D0
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      B4=-4.0D0*z2*ln2**2 + 2.0D0/3.0D0*ln2**4 - 13.0D0/2.0D0*z4
     &  + 16.0D0*li4half
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
      w(1)=z**(-1)
      w(2)=Hr1(-1)
      w(3)=Hr1(0)
      w(4)=Hr2(0,-1)
      w(5)=Hr2(0,1)
      w(6)=Hr1(1)
      w(7)=Hr3(0,-1,-1)
      w(8)=Hr3(0,-1,1)
      w(9)=Hr3(0,0,-1)
      w(10)=Hr3(0,0,1)
      w(11)=Hr3(0,1,-1)
      w(12)=Hr3(0,1,1)
      w(13)=Hr2(-1,1)
      w(14)=Hr4(0,0,0,1)
      w(15)=Hr4(0,-1,-1,-1)
      w(16)=Hr4(0,-1,-1,1)
      w(17)=Hr4(0,-1,0,1)
      w(18)=Hr4(0,-1,1,-1)
      w(19)=Hr4(0,-1,1,1)
      w(20)=Hr4(0,0,-1,-1)
      w(21)=Hr4(0,0,-1,1)
      w(22)=Hr4(0,0,0,-1)
      w(23)=Hr5(0,0,0,0,1)
      w(24)=Hr4(0,0,1,-1)
      w(25)=Hr4(0,0,1,1)
      w(26)=Hr4(0,1,-1,-1)
      w(27)=Hr4(0,1,-1,1)
      w(28)=Hr4(0,1,1,-1)
      w(29)=Hr4(0,1,1,1)
      w(30)=z - 2
      w(31)=2*w(1)
      w(32)=w(30) + w(31)
      w(33)=2*w(6)
      w(34)=w(32)*w(33)
      w(35)=20*w(1)
      w(36)=8*z
      w(37)=1063 + w(36)
      w(37)=z*w(37)
      w(37)=w(35) + 644 + w(37)
      w(38)=1.D0/3.D0*w(3)
      w(39)=29*z
      w(40)=14 - w(39)
      w(40)=w(40)*w(38)
      w(37)=w(40) + 1.D0/9.D0*w(37) - w(34)
      w(37)=w(3)*w(37)
      w(40)=w(32)*w(6)
      w(41)=13*z
      w(42)=37 + w(36)
      w(42)=w(42)*w(41)
      w(42)=166 + w(42)
      w(42)= - w(40) + 1.D0/3.D0*w(42) - 278*w(1)
      w(43)=1.D0/3.D0*w(6)
      w(42)=w(42)*w(43)
      w(44)=z + 2
      w(45)=w(44) + w(31)
      w(46)=w(45)*w(3)
      w(47)=2*z
      w(48)=16.D0/3.D0 + z
      w(48)=w(48)*w(47)
      w(48)=193.D0/3.D0*w(1) + 115.D0/3.D0 + w(48)
      w(48)=2*w(48) + w(46)
      w(49)=w(45)*w(2)
      w(48)=1.D0/3.D0*w(48) + 5*w(49)
      w(50)=2*w(2)
      w(48)=w(48)*w(50)
      w(51)=4.D0/3.D0*w(5)
      w(52)=43*z
      w(53)= - 22*w(1) - 62 - w(52)
      w(53)=w(53)*w(51)
      w(54)=82*w(1)
      w(55)=4805.D0/27.D0 + w(36)
      w(55)=z*w(55)
      w(56)=10*w(1)
      w(57)=w(56) - 62 + 17*z
      w(57)=w(4)*w(57)
      w(58)= - 1838 + 521*z
      w(58)=1.D0/3.D0*w(58) + 346*w(1)
      w(58)=z2*w(58)
      w(37)=2.D0/5.D0*w(58) + w(53) + 2.D0/3.D0*w(57) + w(48) + w(37)
     &  + w(42) + w(54) + 8525.D0/27.D0 + w(55)
      w(37)=z2*w(37)
      w(42)=41 + 23.D0/3.D0*z
      w(48)=7*z
      w(42)=w(42)*w(48)
      w(53)=1703.D0/3.D0 + 50*w(1)
      w(53)=w(53)*w(1)
      w(42)=w(53) + w(42) + 479
      w(53)=4.D0/9.D0*z
      w(55)=w(53) + 3
      w(55)=w(55)*z
      w(57)=31.D0/9.D0*w(1)
      w(55)=w(55) - w(57)
      w(58)=w(55)*w(33)
      w(42)=w(58) + 1.D0/27.D0*w(42)
      w(59)=8*w(40)
      w(60)= - 13*w(30) - w(56)
      w(60)=w(3)*w(60)
      w(61)=131.D0/3.D0 + w(36)
      w(61)=z*w(61)
      w(60)=w(60) + w(59) + 101*w(1) + 16 + w(61)
      w(60)=w(3)*w(60)
      w(61)=133*w(1) + 70 + 41*z
      w(62)=1.D0/9.D0*w(61)
      w(63)= - 11.D0/3.D0*w(49) - w(62) + w(46)
      w(63)=w(2)*w(63)
      w(60)=8.D0/3.D0*w(4) + 2.D0/3.D0*w(63) + 1.D0/9.D0*w(60) - w(42)
      w(60)=w(4)*w(60)
      w(63)=16*w(40)
      w(64)=32*z
      w(65)= - 287 - w(64)
      w(65)=z*w(65)
      w(65)= - 331*w(1) + 34 + w(65)
      w(65)=1.D0/3.D0*w(65) - w(63)
      w(66)= - 94 + 25*z
      w(66)=1.D0/9.D0*w(66) + w(31)
      w(67)=2*w(3)
      w(66)=w(66)*w(67)
      w(65)= - 2.D0/3.D0*w(49) + 1.D0/9.D0*w(65) + w(66)
      w(65)=w(9)*w(65)
      w(66)=8*w(3)
      w(68)= - 1 - 5.D0/9.D0*z
      w(68)=w(68)*w(66)
      w(69)=16*z
      w(70)=727.D0/3.D0 + w(69)
      w(70)=z*w(70)
      w(70)= - 115*w(1) + 610.D0/3.D0 + w(70)
      w(68)=w(68) + 1.D0/9.D0*w(70) - w(34)
      w(68)=w(12)*w(68)
      w(60)=w(65) + w(68) + w(60)
      w(65)=64*z
      w(68)=419 + w(65)
      w(68)=z*w(68)
      w(68)= - 202 + w(68)
      w(68)=1.D0/9.D0*w(68) + w(56)
      w(68)=2*w(68) - 29.D0/3.D0*w(40)
      w(68)=w(6)*w(68)
      w(70)= - 565 - w(65)
      w(70)=w(70)*w(47)
      w(70)= - 761*w(1) + 1249 + w(70)
      w(68)=4.D0/9.D0*w(70) + w(68)
      w(68)=w(6)*w(68)
      w(70)=3613.D0/3.D0 + 280*z
      w(70)=z*w(70)
      w(70)= - 43.D0/3.D0*w(1) - 1889.D0/3.D0 + w(70)
      w(68)=8.D0/9.D0*w(70) + w(68)
      w(68)=w(6)*w(68)
      w(70)=w(26)*w(45)
      w(71)= - 37861.D0/3.D0 + 6136*z
      w(71)=z*w(71)
      w(71)= - 167689.D0/9.D0*w(1) + 205262.D0/9.D0 + w(71)
      w(72)=14*w(1)
      w(73)=w(72) - 46 - z
      w(73)=w(29)*w(73)
      w(74)= - 26*w(1) + 202 - 37*z
      w(74)=w(22)*w(74)
      w(68)=w(68) + 8*w(74) - 32*w(70) + 1.D0/3.D0*w(71) + 16*w(73)
      w(70)=209 + w(69)
      w(70)=z*w(70)
      w(70)= - 142 + w(70)
      w(70)=1.D0/3.D0*w(70) - 19*w(1)
      w(70)=w(70)*w(43)
      w(71)=74*z
      w(73)=731.D0/3.D0 - w(71)
      w(73)=z*w(73)
      w(73)= - 1546.D0/3.D0 + w(73)
      w(74)=4*w(1)
      w(75)=361.D0/27.D0 + w(74)
      w(75)=w(1)*w(75)
      w(70)=w(70) + 1.D0/9.D0*w(73) + w(75)
      w(70)=w(6)*w(70)
      w(73)=1787.D0/3.D0 - 3904*z
      w(73)=z*w(73)
      w(73)= - 4648*w(1) - 29839.D0/3.D0 + w(73)
      w(70)=1.D0/27.D0*w(73) + w(70)
      w(73)= - 71 - w(36)
      w(73)=z*w(73)
      w(73)=70*w(1) + 34 + w(73)
      w(73)=1.D0/3.D0*w(73) - w(59)
      w(73)=w(73)*w(33)
      w(75)= - 11 - 8.D0/3.D0*z
      w(75)=z*w(75)
      w(75)=w(34) - 20.D0/3.D0 + w(75)
      w(76)=10 + z
      w(76)=w(3)*w(76)
      w(75)=2*w(75) + w(76)
      w(75)=w(75)*w(38)
      w(76)=1129 + 544*z
      w(76)=z*w(76)
      w(76)=9148 + w(76)
      w(73)=w(75) + 1.D0/9.D0*w(76) + w(73)
      w(73)=w(73)*w(38)
      w(70)=2*w(70) + w(73)
      w(70)=w(3)*w(70)
      w(68)=1.D0/3.D0*w(68) + w(70)
      w(70)=4*z
      w(73)=w(70) - 133.D0/3.D0
      w(73)=w(73)*z
      w(75)= - 80.D0/3.D0 + 13.D0/3.D0*w(1)
      w(73)=w(73) + w(75)
      w(76)= - 2*w(49) - 4*w(46) + w(73)
      w(76)=w(2)*w(76)
      w(77)=4*w(40)
      w(78)= - 101 - 32.D0/3.D0*z
      w(78)=z*w(78)
      w(78)=w(77) + 221.D0/3.D0*w(1) + 38 + w(78)
      w(78)=w(6)*w(78)
      w(71)= - 5083.D0/3.D0 + w(71)
      w(71)=z*w(71)
      w(71)= - 4486.D0/3.D0 + w(71)
      w(79)= - 323.D0/9.D0 - w(74)
      w(79)=w(1)*w(79)
      w(71)=2.D0/3.D0*w(78) + 1.D0/9.D0*w(71) + w(79)
      w(78)=11.D0/3.D0 + w(47)
      w(78)=w(78)*w(70)
      w(78)=w(77) + w(1) + 149.D0/3.D0 + w(78)
      w(79)=2.D0/3.D0*w(1)
      w(80)=w(79) + 2 + 7.D0/9.D0*z
      w(80)=w(3)*w(80)
      w(78)=2.D0/9.D0*w(78) + w(80)
      w(78)=w(78)*w(67)
      w(80)=4.D0/3.D0 + z
      w(80)=7*w(80) + w(79)
      w(51)=w(80)*w(51)
      w(80)=w(56) + 14 + 11*z
      w(80)=w(4)*w(80)
      w(51)=w(51) + 8.D0/9.D0*w(80) + 4.D0/9.D0*w(76) + 1.D0/3.D0*w(71)
     &  + w(78)
      w(51)=w(5)*w(51)
      w(71)= - 389 - 104*z
      w(71)=z*w(71)
      w(54)= - w(54) - 982 + w(71)
      w(71)=23*z
      w(72)= - w(72) - 74 - w(71)
      w(72)=w(72)*w(67)
      w(54)=20*w(49) + w(72) + 1.D0/3.D0*w(54) + w(77)
      w(54)=w(10)*w(54)
      w(72)=4.D0/9.D0*w(49)
      w(76)=w(79) + z
      w(78)= - 10.D0/9.D0 - w(76)
      w(78)=w(78)*w(67)
      w(73)=w(72) - 1.D0/9.D0*w(73) + w(78)
      w(73)=w(11)*w(73)
      w(78)=22.D0/9.D0*w(1) + 178.D0/9.D0 + w(48)
      w(78)=w(14)*w(78)
      w(73)=w(73) + w(78)
      w(46)=44.D0/9.D0*w(49) + 2.D0/9.D0*w(61) - w(46)
      w(46)=w(2)*w(46)
      w(61)=w(45)*w(67)
      w(78)=25 - w(69)
      w(78)=z*w(78)
      w(61)=w(61) - 275*w(1) - 130 + w(78)
      w(61)=w(3)*w(61)
      w(42)=1.D0/3.D0*w(46) + 2*w(42) + 1.D0/27.D0*w(61)
      w(42)=w(3)*w(42)
      w(46)=32.D0/9.D0*w(45)
      w(61)=w(12)*w(46)
      w(42)=w(61) + w(42)
      w(42)=w(42)*w(50)
      w(50)=133.D0/3.D0 + w(70)
      w(50)=z*w(50)
      w(50)=w(50) - w(75)
      w(61)=26.D0/9.D0 - w(76)
      w(61)=w(61)*w(67)
      w(50)=w(72) + 1.D0/9.D0*w(50) + w(61)
      w(50)=w(8)*w(50)
      w(61)= - 176.D0/9.D0*w(15) - 80.D0/9.D0*w(24)
      w(61)=w(45)*w(61)
      w(72)=64.D0/9.D0*w(45)
      w(75)= - w(19) - w(28)
      w(75)=w(72)*w(75)
      w(72)= - w(27)*w(72)
      w(78)=5*z
      w(80)=w(78) + w(56)
      w(81)= - 82 - w(80)
      w(81)=w(21)*w(81)
      w(82)= - w(16)*w(46)
      w(46)= - w(18)*w(46)
      w(83)=128.D0/3.D0*w(25)
      w(84)= - w(76)*w(83)
      w(85)=w(31) + z
      w(86)=w(85) - 5
      w(87)=w(86)*li4half
      w(88)=w(13)*z**2
      w(89)=w(3)*w(88)
      w(37)=w(84) + w(46) + w(75) + 8*w(73) + w(82) + 4.D0/9.D0*w(54)
     &  - 256.D0/3.D0*w(87) + 16.D0/9.D0*w(81) + w(72) + 2.D0/3.D0*
     & w(37) + 2*w(51) + 8*w(50) + w(42) - 64.D0/9.D0*w(89) + 1.D0/3.D0
     & *w(68) + w(61) + 4*w(60)
      w(42)=CA*TF
      w(37)=w(42)*w(37)
      w(46)=5*w(1)
      w(50)=w(70) + w(46) - 5
      w(51)= - w(40) + 2*w(50)
      w(51)=w(51)*w(43)
      w(54)= - w(74) + 4 - 7.D0/3.D0*z
      w(54)=10*w(54) + w(51)
      w(54)=w(6)*w(54)
      w(60)=868*w(1) - 868 + 569*z
      w(54)=4.D0/27.D0*w(60) + w(54)
      w(50)= - w(40) + 2.D0/3.D0*w(50)
      w(60)=z2*w(50)
      w(54)=1.D0/3.D0*w(54) + w(60)
      w(54)=NF*w(54)
      w(41)=w(51) - 58.D0/3.D0*w(1) + 58.D0/3.D0 - w(41)
      w(41)=w(41)*w(33)
      w(51)=1214*w(1) - 1214 + 1129*z
      w(41)=1.D0/27.D0*w(51) + w(41)
      w(51)=2*z2
      w(50)=w(50)*w(51)
      w(41)=w(54) + 1.D0/3.D0*w(41) + w(50)
      w(41)=8.D0/3.D0*w(41)
      w(50)=TF**2
      w(41)=w(50)*w(41)
      w(37)=w(41) + w(37)
      w(41)= - 203.D0/3.D0 - w(65)
      w(41)=z*w(41)
      w(54)=130 - 173*z
      w(54)=w(3)*w(54)
      w(41)=w(54) + w(59) - 560.D0/3.D0*w(1) + 472.D0/3.D0 + w(41)
      w(41)=w(3)*w(41)
      w(54)= - 3817.D0/3.D0 + 224*z
      w(54)=w(54)*w(70)
      w(54)= - 8584.D0/3.D0*w(1) + 18545.D0/3.D0 + w(54)
      w(59)= - 55.D0/3.D0 - w(69)
      w(59)=z*w(59)
      w(59)= - 17*w(40) + 766.D0/3.D0*w(1) - 370.D0/3.D0 + w(59)
      w(59)=w(59)*w(33)
      w(41)=w(41) + 1.D0/3.D0*w(54) + w(59)
      w(53)=w(53) - 3
      w(53)=w(53)*z
      w(53)=w(53) + w(57)
      w(54)=w(53)*w(2)
      w(57)=w(79) + w(30)
      w(59)=16*w(4)
      w(60)= - w(57)*w(59)
      w(61)=5 + w(47)
      w(46)=4*w(61) + w(46)
      w(46)=w(5)*w(46)
      w(61)=w(30)*w(3)
      w(68)=457 - 143*z
      w(68)=1.D0/15.D0*w(68) - 13*w(1)
      w(68)=1.D0/3.D0*w(68) + 12.D0/5.D0*w(61)
      w(68)=z2*w(68)
      w(41)=16*w(68) + 16.D0/9.D0*w(46) + w(60) + 1.D0/9.D0*w(41) - 8*
     & w(54)
      w(41)=w(41)*w(51)
      w(46)=w(57)*w(67)
      w(46)=w(46) - w(55)
      w(46)=w(8)*w(46)
      w(44)=w(79) + w(44)
      w(55)=w(44)*w(67)
      w(53)=w(55) + w(53)
      w(53)=w(11)*w(53)
      w(55)=82 - w(71)
      w(35)=1.D0/3.D0*w(55) - w(35)
      w(35)=1.D0/3.D0*w(35) + 4*w(61)
      w(35)=w(14)*w(35)
      w(35)=w(35) + w(46) + w(53)
      w(46)= - 25 + 16.D0/3.D0*z
      w(46)=w(46)*w(47)
      w(46)=w(34) - 176.D0/3.D0*w(1) + 107 + w(46)
      w(53)=4*w(6)
      w(46)=w(46)*w(53)
      w(53)= - 9229 + 896*z
      w(47)=w(53)*w(47)
      w(47)=1339 + w(47)
      w(46)=1.D0/9.D0*w(47) + w(46)
      w(47)=1051 - w(65)
      w(47)=z*w(47)
      w(47)= - 368 + w(47)
      w(47)=1.D0/9.D0*w(47) + w(63)
      w(52)=14 - w(52)
      w(52)=w(52)*w(38)
      w(47)=2*w(47) + w(52)
      w(47)=w(3)*w(47)
      w(46)=2*w(46) + w(47)
      w(46)=w(46)*w(38)
      w(47)=w(1) + 37
      w(47)=w(47)*w(31)
      w(52)=1.D0/3.D0*z
      w(53)=203*z
      w(55)=1303 - w(53)
      w(55)=w(55)*w(52)
      w(55)= - 494 + w(55)
      w(55)=1.D0/3.D0*w(55) + w(47)
      w(60)= - 146*w(1) + 122 - 97*z
      w(60)=1.D0/9.D0*w(60) + w(34)
      w(60)=w(6)*w(60)
      w(55)=2*w(55) + w(60)
      w(33)=w(55)*w(33)
      w(55)=26941.D0/81.D0 - 88*z
      w(55)=w(55)*w(70)
      w(33)=w(33) + 4016.D0/81.D0*w(1) + 72235.D0/81.D0 + w(55)
      w(33)=4*w(33) + w(46)
      w(33)=w(33)*w(38)
      w(46)=w(1) + 166.D0/9.D0
      w(46)=w(46)*w(79)
      w(55)=8 + 161.D0/27.D0*z
      w(55)=w(55)*z
      w(46)=w(55) + w(46) + 41.D0/3.D0
      w(46)=w(58) + 1.D0/3.D0*w(46)
      w(55)=2.D0/3.D0*z
      w(58)= - 1 - w(55)
      w(55)=w(58)*w(55)
      w(57)=w(3)*w(57)
      w(55)=2.D0/3.D0*w(57) - 76.D0/27.D0*w(1) - 1 + w(55)
      w(55)=w(55)*w(67)
      w(55)=w(55) + w(46)
      w(55)=w(55)*w(59)
      w(57)=58*w(1) - 58 + 35*z
      w(57)=1.D0/3.D0*w(57) - w(40)
      w(43)=w(57)*w(43)
      w(53)= - 3868.D0/3.D0 + w(53)
      w(53)=z*w(53)
      w(53)=3965.D0/3.D0 + w(53)
      w(43)=w(43) + 1.D0/9.D0*w(53) - w(47)
      w(47)=89 - w(70)
      w(47)=z*w(47)
      w(47)=130*w(1) - 430 + w(47)
      w(47)=1.D0/3.D0*w(47) + w(34)
      w(53)= - 56*w(1) - 4 + w(78)
      w(53)=w(3)*w(53)
      w(47)=4*w(47) + w(53)
      w(38)=w(47)*w(38)
      w(38)=2*w(43) + w(38)
      w(43)=w(4)*w(44)
      w(38)= - 8*w(43) + 1.D0/3.D0*w(38) - 4*w(54)
      w(38)=w(5)*w(38)
      w(43)= - 1 + 4.D0/27.D0*z
      w(43)=z*w(43)
      w(43)=w(43) + 31.D0/27.D0*w(1)
      w(43)=w(43)*w(67)
      w(43)=w(43) - w(46)
      w(43)=w(2)*w(43)
      w(43)=16*w(43) + 256.D0/9.D0*w(88)
      w(43)=w(3)*w(43)
      w(44)= - 36170 + 23231.D0/3.D0*z
      w(36)=w(44)*w(36)
      w(36)= - 63052.D0/3.D0*w(1) + 231205 + w(36)
      w(31)= - 5.D0/9.D0*w(30) - w(31)
      w(31)=w(29)*w(31)
      w(44)= - 4 + w(76)
      w(44)=w(22)*w(44)
      w(31)=32*w(44) + 1.D0/243.D0*w(36) + 8*w(31)
      w(36)= - 250*w(1) + 154 - 179*z
      w(36)=2.D0/3.D0*w(36) + 29*w(40)
      w(36)=w(6)*w(36)
      w(44)=661*w(1) - 439 + 236*z
      w(36)=4.D0/3.D0*w(44) + w(36)
      w(36)=w(6)*w(36)
      w(44)=2276*w(1) - 4784 + 2713*z
      w(36)=4.D0/9.D0*w(44) + w(36)
      w(36)=w(6)*w(36)
      w(31)=2*w(31) + 1.D0/27.D0*w(36)
      w(36)=50 + w(39)
      w(34)=w(34) + 1.D0/3.D0*w(36) + w(56)
      w(36)=w(1) + 1.D0/3.D0*w(30)
      w(39)= - w(36)*w(66)
      w(34)=1.D0/3.D0*w(34) + w(39)
      w(34)=w(12)*w(34)
      w(39)= - 205 - w(69)
      w(39)=z*w(39)
      w(39)= - 344*w(1) + 1427 + w(39)
      w(39)=1.D0/3.D0*w(39) - w(77)
      w(44)=88*w(1) - 44 + w(48)
      w(44)=1.D0/9.D0*w(44) - w(61)
      w(44)=w(3)*w(44)
      w(39)=1.D0/9.D0*w(39) + w(44)
      w(39)=w(10)*w(39)
      w(44)= - w(23) + z5
      w(46)=320*w(30)
      w(44)=w(46)*w(44)
      w(46)=5 + 8.D0/9.D0*z
      w(46)=w(46)*w(52)
      w(47)= - 2.D0/9.D0*w(1) + 1 - w(52)
      w(47)=w(3)*w(47)
      w(46)=4*w(47) + 5.D0/3.D0*w(1) + 1 + w(46)
      w(46)=w(9)*w(46)
      w(36)=w(36)*w(83)
      w(31)=w(36) + 16*w(39) + 1024.D0/3.D0*w(87) + 512*w(21) + 64*
     & w(46) + w(41) + 8*w(38) + w(55) + 16.D0/3.D0*w(34) + 2*w(31) + 
     & w(33) + w(43) + w(44) + 32*w(35)
      w(33)=CF*TF
      w(31)=w(31)*w(33)
      w(34)= - 1381 + w(64)
      w(34)=z*w(34)
      w(34)=640*w(1) + 1768 + w(34)
      w(34)=1.D0/3.D0*w(34) - 28*w(40)
      w(35)=17.D0/9.D0*w(1) - 11.D0/9.D0 + w(78)
      w(35)=w(35)*w(66)
      w(34)= - 28.D0/3.D0*w(49) + 1.D0/9.D0*w(34) + w(35)
      w(34)=w(34)*w(42)
      w(35)= - 56.D0/9.D0*NF + 128.D0/9.D0
      w(32)=w(50)*w(32)*w(35)
      w(35)= - 1489 + 400*z
      w(35)=z*w(35)
      w(35)=1076*w(1) - 4168 + w(35)
      w(35)=1.D0/3.D0*w(35) + 44*w(40)
      w(30)=w(30)*w(67)
      w(30)=w(30) - 80.D0/3.D0*w(1) - 98 + 67*z
      w(30)=w(30)*w(67)
      w(30)=1.D0/3.D0*w(35) + w(30)
      w(30)=w(30)*w(33)
      w(30)=1.D0/3.D0*w(30) + w(34) + w(32)
      w(30)=z3*w(30)
      w(32)=w(86)*w(42)
      w(34)=z2*w(32)
      w(35)=w(33)*w(86)
      w(36)= - w(51)*w(35)
      w(34)=w(34) + w(36)
      w(32)= - w(32) + 2*w(35)
      w(35)=ln2**2
      w(32)=w(32)*w(35)
      w(32)=2*w(34) + 1.D0/3.D0*w(32)
      w(32)=w(35)*w(32)
      w(34)= - 10 - w(85)
      w(34)=w(3)*w(34)
      w(34)=22.D0/3.D0*w(49) + w(62) + w(34)
      w(34)=w(7)*w(34)
      w(35)=w(20)*w(45)
      w(34)=w(34) + w(35)
      w(34)=16.D0/3.D0*w(34)
      w(34)=w(42)*w(34)
      w(35)= - 38 - w(80)
      w(35)=w(35)*w(42)
      w(33)=1.D0/9.D0*w(35) + 8*w(33)
      w(33)=w(17)*w(33)
      w(30)=64.D0/3.D0*w(32) + 32*w(33) + 4*w(30) + w(34) + 2*w(37) + 
     & w(31)
*
      AGQ = CF*w(30)*as**3
*
      RETURN
      END      
      REAL*8 FUNCTION GQ(z,nf,as,LL)
*     ------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     OME AgQ  without a_gqQ^(3) 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(65)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL,ln2,z5,li4half
* 
      ln2= 0.69314718055994530942D0
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      z5 = 1.0369277551433699263D0
      li4half = 0.51747906167389938633D0
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 5
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
      w(1)=z**(-1)
      w(2)=Hr1(-1)
      w(3)=Hr1(0)
      w(4)=Hr2(-1,-1)
      w(5)=Hr3(-1,0,-1)
      w(6)=Hr3(-1,0,1)
      w(7)=Hr1(1)
      w(8)=Hr3(0,0,1)
      w(9)=Hr2(0,1)
      w(10)=Hr2(0,-1)
      w(11)=Hr4(0,0,0,1)
      w(12)=Hr3(0,-1,-1)
      w(13)=Hr3(0,0,-1)
      w(14)=Hr5(0,0,0,0,1)
      w(15)=Hr3(0,1,1)
      w(16)=Hr4(0,1,1,1)
      w(17)=8*z
      w(18)=457.D0/3.D0 - 20*z
      w(18)=w(18)*w(17)
      w(18)=176*w(1) + 3407.D0/3.D0 + w(18)
      w(19)=2*w(1)
      w(20)=z - 2
      w(21)=w(19) + w(20)
      w(22)=w(21)*w(7)
      w(23)=35 - 8.D0/3.D0*z
      w(23)=z*w(23)
      w(23)=104.D0/3.D0*w(1) - 55 + w(23)
      w(23)=1.D0/3.D0*w(23) - w(22)
      w(24)=4*w(7)
      w(23)=w(23)*w(24)
      w(18)=1.D0/9.D0*w(18) + w(23)
      w(23)=8*w(22)
      w(25)=7*z
      w(26)=w(20)*w(3)
      w(27)=w(26) - 2 + w(25)
      w(28)=2*w(3)
      w(27)=w(27)*w(28)
      w(29)= - 640 - 973*z
      w(27)=w(27) + 1.D0/3.D0*w(29) - w(23)
      w(29)=1.D0/3.D0*w(3)
      w(27)=w(27)*w(29)
      w(18)=2*w(18) + w(27)
      w(18)=w(3)*w(18)
      w(27)=w(21)*w(24)
      w(30)=56*w(1)
      w(31)=191 - 32*z
      w(31)=z*w(31)
      w(31)= - w(30) + 134 + w(31)
      w(32)=2.D0/3.D0*w(3)
      w(33)=13*z
      w(34)=8 - w(33)
      w(34)=w(34)*w(32)
      w(31)=w(34) + 1.D0/9.D0*w(31) + w(27)
      w(31)=w(31)*w(28)
      w(34)=17*z
      w(35)=16*w(1)
      w(36)= - w(35) + 2 + w(34)
      w(37)= - w(20)*w(28)
      w(36)=1.D0/3.D0*w(36) + w(37)
      w(36)=w(36)*w(28)
      w(37)= - 109 + 32.D0/3.D0*z
      w(37)=z*w(37)
      w(37)= - w(27) - 248.D0/3.D0*w(1) + 184 + w(37)
      w(36)=1.D0/3.D0*w(37) + w(36)
      w(37)=2.D0/3.D0*LL
      w(36)=w(36)*w(37)
      w(38)=17*w(1) - 23 + w(33)
      w(38)=2.D0/3.D0*w(38) + w(22)
      w(38)=w(38)*w(24)
      w(39)= - 1315 + 896.D0/3.D0*z
      w(39)=z*w(39)
      w(39)= - 680.D0/3.D0*w(1) + 1156 + w(39)
      w(38)=1.D0/3.D0*w(39) + w(38)
      w(31)=w(36) + 1.D0/3.D0*w(38) + w(31)
      w(31)=LL*w(31)
      w(36)=2*w(7)
      w(38)=w(21)*w(36)
      w(39)= - w(38) - 86*w(1) + 110 - 67*z
      w(39)=w(7)*w(39)
      w(40)=115*z
      w(41)= - 49*w(1) + 16 + w(40)
      w(39)=4.D0/3.D0*w(41) + w(39)
      w(39)=w(39)*w(36)
      w(41)= - 15955 + 3616*z
      w(41)=z*w(41)
      w(41)=4676*w(1) + 7324 + w(41)
      w(39)=1.D0/3.D0*w(41) + w(39)
      w(18)=w(31) + 1.D0/9.D0*w(39) + w(18)
      w(31)=2*LL
      w(18)=w(18)*w(31)
      w(39)=62*w(1)
      w(41)=4 + w(33)
      w(41)=5*w(41) - w(39)
      w(42)=5*z
      w(43)=4*w(1)
      w(44)= - w(43) - 2 - w(42)
      w(44)=w(44)*w(28)
      w(45)=w(20)*LL
      w(41)= - w(45) + w(44) + 1.D0/3.D0*w(41) + w(27)
      w(41)=LL*w(41)
      w(44)=4*z
      w(46)=224*z
      w(47)=1139 - w(46)
      w(47)=w(47)*w(44)
      w(47)=1976*w(1) - 5107 + w(47)
      w(48)= - 94*w(1) + 118 - 71*z
      w(48)=1.D0/3.D0*w(48) - w(22)
      w(48)=w(48)*w(36)
      w(47)=1.D0/9.D0*w(47) + w(48)
      w(48)=55*z
      w(49)= - 38 + w(48)
      w(49)=w(49)*w(29)
      w(50)= - 319 + 64*z
      w(50)=z*w(50)
      w(50)=112*w(1) - 232 + w(50)
      w(23)=w(49) + 1.D0/9.D0*w(50) - w(23)
      w(23)=w(3)*w(23)
      w(33)=58.D0/5.D0 - w(33)
      w(33)=12.D0/5.D0*w(45) + 1.D0/3.D0*w(33) - 16.D0/5.D0*w(1)
      w(33)=z2*w(33)
      w(23)=8*w(33) + 8.D0/3.D0*w(41) + 1.D0/3.D0*w(47) + w(23)
      w(33)=2*z2
      w(23)=w(23)*w(33)
      w(41)=w(44) - 5
      w(47)=5*w(1)
      w(49)=w(41) + w(47)
      w(50)= - w(22) + 2*w(49)
      w(51)=w(50)*w(7)
      w(52)=2*z
      w(53)= - 73 + 82.D0/3.D0*z
      w(53)=w(53)*w(52)
      w(53)=w(51) - 575.D0/3.D0*w(1) + 245 + w(53)
      w(53)=w(7)*w(53)
      w(54)= - 19801.D0/3.D0 + 5024*z
      w(54)=z*w(54)
      w(54)=3088.D0/3.D0*w(1) - 1312.D0/3.D0 + w(54)
      w(53)=1.D0/3.D0*w(54) + 8*w(53)
      w(54)=z + 1
      w(55)=16*z
      w(56)=w(54)*w(55)
      w(56)= - 41 + w(56)
      w(56)=1.D0/3.D0*w(56) - w(38)
      w(57)= - 22 + 29*z
      w(57)=w(3)*w(57)
      w(56)=4*w(56) + w(57)
      w(56)=w(3)*w(56)
      w(57)=19 - w(52)
      w(57)=w(57)*w(44)
      w(57)=52*w(1) - 125 + w(57)
      w(57)=w(57)*w(24)
      w(58)=2779 - 656*z
      w(58)=z*w(58)
      w(58)= - 2888 + w(58)
      w(57)=1.D0/3.D0*w(58) + w(57)
      w(56)=2*w(57) + w(56)
      w(56)=w(3)*w(56)
      w(53)=2*w(53) + w(56)
      w(53)=w(3)*w(53)
      w(56)=w(2)*w(3)
      w(57)=4.D0/9.D0*z
      w(58)=w(57) - 3
      w(58)=w(58)*z
      w(58)=w(58) + 31.D0/9.D0*w(1)
      w(59)=w(58)*w(56)
      w(60)=w(20) + 2.D0/3.D0*w(1)
      w(61)= - w(60)*w(28)
      w(58)=w(61) + w(58)
      w(58)=w(10)*w(58)
      w(58)=w(59) - w(58)
      w(59)=w(13)*w(60)
      w(58)=128*w(59) - 32*w(58)
      w(58)=LL*w(58)
      w(35)=w(35) + 142 - 269*z
      w(35)=1.D0/3.D0*w(35) + 14*w(26)
      w(35)=w(35)*w(28)
      w(59)=233 - 416.D0/3.D0*z
      w(59)=z*w(59)
      w(35)=w(35) + 28.D0/3.D0*w(22) - 712.D0/9.D0*w(1) - 320 + 1.D0/3.D
     & 0*w(59)
      w(59)=w(26) - 5.D0/3.D0*w(1) - w(41)
      w(59)=LL*w(59)
      w(35)=1.D0/3.D0*w(35) + 16*w(59)
      w(59)=4*z3
      w(35)=w(35)*w(59)
      w(60)=w(14) - z5
      w(60)=64*w(60)
      w(60)=w(20)*w(60)
      w(61)=4.D0/9.D0*w(22)
      w(62)= - w(61) + 40.D0/9.D0*w(1) - 4 + 31.D0/9.D0*z
      w(62)=w(7)*w(62)
      w(40)= - 116*w(1) + 77 - w(40)
      w(40)=2.D0/27.D0*w(40) + w(62)
      w(40)=w(7)*w(40)
      w(62)=440*w(1) - 14 + 193*z
      w(40)=2.D0/81.D0*w(62) + w(40)
      w(24)=w(40)*w(24)
      w(40)=20*w(1)
      w(62)= - w(40) + 16 - 11*z
      w(63)=LL*w(43)
      w(62)=1.D0/3.D0*w(62) + w(63)
      w(62)=w(15)*w(62)
      w(63)=w(43) + w(20)
      w(64)=16.D0/3.D0*w(16)
      w(63)=w(63)*w(64)
      w(65)=322861.D0/2.D0 - 153536.D0/3.D0*z
      w(65)=z*w(65)
      w(65)=163994.D0/3.D0*w(1) - 165035 + w(65)
      w(18)=w(23) + w(35) + w(63) + 16.D0/3.D0*w(62) + w(18) + 1.D0/9.D0
     & *w(53) + 1.D0/81.D0*w(65) + w(24) + w(58) + w(60)
      w(23)=as**3
      w(24)=w(23)*CF
      w(18)=w(18)*w(24)
      w(35)=1 - w(52)
      w(35)=z*w(35)
      w(35)= - 8 + w(35)
      w(53)=z + 4
      w(58)=w(3)*w(53)
      w(35)=w(58) + w(22) + 4.D0/3.D0*w(35) - w(1)
      w(58)=4*w(3)
      w(35)=w(35)*w(58)
      w(60)=2.D0/3.D0*w(22)
      w(62)= - 1 - 1.D0/9.D0*z
      w(62)=w(62)*w(52)
      w(62)=7.D0/3.D0 + w(62)
      w(54)=w(54) + w(1)
      w(63)=w(3)*w(54)
      w(62)=4.D0/3.D0*w(63) + w(60) + 2*w(62) - 35.D0/9.D0*w(1)
      w(62)=w(62)*w(31)
      w(63)=5 + w(44)
      w(63)=z*w(63)
      w(63)=14 + w(63)
      w(63)=1.D0/3.D0*w(63) - 7*w(1)
      w(63)=2*w(63) - w(22)
      w(63)=w(63)*w(36)
      w(34)= - 26 + w(34)
      w(34)=z*w(34)
      w(34)=74 + w(34)
      w(34)=4*w(34) - 191*w(1)
      w(34)=w(62) + w(35) + 1.D0/3.D0*w(34) + w(63)
      w(34)=LL*w(34)
      w(35)= - 4.D0/3.D0 - z
      w(35)=w(35)*w(44)
      w(35)= - 149.D0/3.D0*w(1) - 1.D0/3.D0 + w(35)
      w(62)=19 - w(55)
      w(62)=z*w(62)
      w(62)= - 122 + w(62)
      w(60)=w(60) + 1.D0/3.D0*w(62) + 50*w(1)
      w(60)=w(7)*w(60)
      w(35)=4.D0/3.D0*w(35) + w(60)
      w(35)=w(7)*w(35)
      w(34)=w(35) + w(34)
      w(35)=19*z
      w(60)=w(40) - 20 + w(35)
      w(60)=2.D0/9.D0*w(60) - w(22)
      w(60)=w(7)*w(60)
      w(62)= - 17 - 53*z
      w(62)=w(62)*w(17)
      w(62)= - 356*w(1) - 1387 + w(62)
      w(60)=1.D0/27.D0*w(62) + w(60)
      w(62)=z + 2
      w(32)= - w(62)*w(32)
      w(63)=151 + w(17)
      w(63)=z*w(63)
      w(63)=292 + w(63)
      w(27)=w(32) + 1.D0/3.D0*w(63) + w(27)
      w(27)=w(27)*w(29)
      w(27)=2*w(60) + w(27)
      w(27)=w(3)*w(27)
      w(32)=56*z
      w(60)= - 1811.D0/27.D0 + w(32)
      w(60)=z*w(60)
      w(27)=w(27) - 4043.D0/27.D0*w(1) + 1648.D0/9.D0 + w(60) + 1.D0/3.D
     & 0*w(34)
      w(27)=w(27)*w(31)
      w(34)=1 - w(17)
      w(34)=z*w(34)
      w(34)=w(39) + 76 + w(34)
      w(39)=w(19) + 4 + w(42)
      w(39)=w(39)*w(28)
      w(42)=w(19) + 5 + w(52)
      w(42)=w(42)*w(31)
      w(34)=w(42) + w(39) + 1.D0/3.D0*w(34) - w(22)
      w(34)=LL*w(34)
      w(39)= - 61 + w(17)
      w(39)=z*w(39)
      w(39)=20 + w(39)
      w(42)= - 10 - z
      w(42)=w(3)*w(42)
      w(39)=w(42) + w(38) + 1.D0/3.D0*w(39) + w(43)
      w(39)=w(3)*w(39)
      w(42)= - 7 - w(55)
      w(42)=z*w(42)
      w(42)=98*w(1) - 70 + w(42)
      w(42)=1.D0/3.D0*w(42) + w(22)
      w(42)=w(7)*w(42)
      w(46)= - 505 - w(46)
      w(46)=z*w(46)
      w(46)=224*w(1) - 1687 + w(46)
      w(34)=4*w(34) + w(39) + 1.D0/9.D0*w(46) + w(42)
      w(39)=w(62) + w(19)
      w(42)=w(39)*w(3)
      w(46)=w(39)*LL
      w(55)= - w(42) + 16.D0/3.D0*w(46)
      w(55)=w(2)*w(55)
      w(60)=w(10)*w(39)
      w(55)=w(55) + w(60)
      w(60)=w(1) + z
      w(60)= - 43 - 28*w(60)
      w(60)=z2*w(60)
      w(34)=4.D0/15.D0*w(60) + 1.D0/3.D0*w(34) + 2*w(55)
      w(34)=w(34)*w(33)
      w(55)=w(57) - 1
      w(55)=w(55)*z
      w(55)=w(46) + w(55) + 10.D0/3.D0 + 82.D0/9.D0*w(1)
      w(57)= - 2.D0/3.D0*w(42) + w(55)
      w(56)=w(57)*w(56)
      w(57)=w(19) - 2 + 5.D0/3.D0*z
      w(28)=w(57)*w(28)
      w(28)=w(28) - w(55)
      w(28)=w(10)*w(28)
      w(28)=w(56) + w(28)
      w(39)=64.D0/3.D0*w(39)
      w(55)= - w(12) - w(5)
      w(39)=w(39)*w(55)
      w(42)=w(4)*w(42)
      w(55)= - w(1) - w(20)
      w(55)=w(13)*w(55)
      w(28)=128.D0/3.D0*w(55) + 64.D0/3.D0*w(42) + w(39) + 8*w(28)
      w(28)=LL*w(28)
      w(39)=13 + 4.D0/3.D0*z
      w(25)=w(39)*w(25)
      w(39)= - w(54)*w(58)
      w(25)=w(39) - 14*w(22) + 155.D0/3.D0*w(1) + 74 + w(25)
      w(39)=26.D0/3.D0*w(1) - 20 + 19.D0/3.D0*z
      w(39)=w(39)*w(31)
      w(25)=1.D0/9.D0*w(25) + w(39)
      w(25)=w(25)*w(59)
      w(39)= - 266 - w(48)
      w(39)=w(39)*w(44)
      w(39)=793 + w(39)
      w(39)=2.D0/3.D0*w(39) - 367*w(1)
      w(42)= - 5 - 8.D0/27.D0*z
      w(42)=z*w(42)
      w(42)=w(61) - 124.D0/27.D0*w(1) + 46.D0/9.D0 + w(42)
      w(42)=w(7)*w(42)
      w(48)=313 + 38*z
      w(48)=z*w(48)
      w(48)= - 317 + w(48)
      w(48)=1.D0/27.D0*w(48) + 10*w(1)
      w(42)=2*w(48) + w(42)
      w(42)=w(7)*w(42)
      w(39)=2.D0/27.D0*w(39) + w(42)
      w(39)=w(39)*w(36)
      w(36)= - w(50)*w(36)
      w(42)=242*w(1) - 242 + 181*z
      w(36)=1.D0/3.D0*w(42) + w(36)
      w(36)=w(7)*w(36)
      w(42)= - 2953 - 880*z
      w(42)=z*w(42)
      w(42)= - 3160 + w(42)
      w(36)=2.D0/3.D0*w(36) + 1.D0/27.D0*w(42) + w(40)
      w(40)=8.D0/3.D0*w(41) + w(26)
      w(40)=w(3)*w(40)
      w(41)= - 157 + 131*z
      w(40)=2.D0/3.D0*w(41) + w(40)
      w(40)=w(40)*w(29)
      w(36)=2*w(36) + w(40)
      w(29)=w(36)*w(29)
      w(36)= - 47 - w(17)
      w(36)=z*w(36)
      w(36)= - 128 + w(36)
      w(40)=LL*w(53)
      w(36)=1.D0/3.D0*w(36) + 8*w(40)
      w(36)=w(15)*w(36)
      w(40)=w(53)*w(64)
      w(32)=14165.D0/3.D0 + w(32)
      w(32)=z*w(32)
      w(32)=13937.D0/3.D0*w(1) - 5960 + w(32)
      w(25)=w(34) + w(25) + w(40) + 4.D0/3.D0*w(36) + w(27) + w(29) + 1.
     & D0/27.D0*w(32) + w(39) + w(28)
      w(27)=w(23)*CA
      w(25)=w(25)*w(27)
      w(28)= - w(22) + 4.D0/3.D0*w(49)
      w(28)=w(28)*w(7)
      w(22)= - w(22) + 2.D0/3.D0*w(49)
      w(29)=LL*w(21)
      w(29)=2*w(22) + w(29)
      w(29)=w(29)*w(31)
      w(30)=43*z + w(30) - 56
      w(29)=w(29) + 2.D0/9.D0*w(30) - w(28)
      w(29)=w(29)*as**2
      w(25)=2.D0/3.D0*w(29) + w(25)
      w(29)= - w(21)*w(58)
      w(32)= - w(53)*w(31)
      w(34)= - 29 + w(17)
      w(34)=z*w(34)
      w(34)= - 56 + w(34)
      w(29)=w(32) + 1.D0/3.D0*w(34) + w(29)
      w(29)=w(29)*w(37)
      w(32)=283 + 76*z
      w(32)=z*w(32)
      w(32)=670 + w(32)
      w(34)=z2*w(53)
      w(29)=4.D0/3.D0*w(34) + w(29) + 1.D0/27.D0*w(32) - w(19)
      w(29)=w(29)*w(27)
      w(32)= - 29 + w(52)
      w(32)=w(32)*w(44)
      w(32)= - 104*w(1) + 145 + w(32)
      w(32)=1.D0/3.D0*w(32) + w(38)
      w(34)=3*z
      w(36)= - w(34) + w(43)
      w(36)=w(3)*w(36)
      w(20)=w(1) + 1.D0/3.D0*w(20)
      w(37)= - w(20)*w(31)
      w(32)=w(37) + 1.D0/3.D0*w(32) + w(36)
      w(32)=LL*w(32)
      w(36)=w(47) - 1 - w(44)
      w(36)=w(3)*w(36)
      w(37)= - 29 - w(44)
      w(37)=z*w(37)
      w(37)= - 46*w(1) + 13 + w(37)
      w(36)=2.D0/3.D0*w(37) + w(36)
      w(36)=w(3)*w(36)
      w(37)=419 - 164*z
      w(37)=z*w(37)
      w(37)=575*w(1) - 709 + w(37)
      w(36)=1.D0/9.D0*w(37) + w(36)
      w(20)=w(20)*w(33)
      w(20)=w(20) + 1.D0/3.D0*w(36) + w(32)
      w(20)=w(20)*w(24)
      w(20)=w(29) + 2*w(20)
      w(20)=w(9)*w(20)
      w(18)=8*w(20) + 2*w(25) + w(18)
      w(18)=CF*w(18)
      w(20)= - w(51) + 2.D0/3.D0*w(30)
      w(20)=w(20)*w(7)
      w(25)=41*w(1) - 41 + 31*z
      w(20)= - w(20) + 16.D0/9.D0*w(25)
      w(20)=1.D0/3.D0*w(20)
      w(25)=2.D0/3.D0*w(21)
      w(29)=LL**2*w(25)
      w(30)=68*w(1) - 68 + w(35)
      w(28)=w(29) + 2.D0/9.D0*w(30) + w(28)
      w(28)=LL*w(28)
      w(28)= - w(20) + w(28)
      w(28)=nf*w(28)
      w(29)= - z2*w(22)
      w(30)= - z3*w(21)
      w(29)=w(29) + 2.D0/3.D0*w(30)
      w(30)=nf + 2
      w(29)=w(30)*w(29)
      w(25)=LL*w(25)
      w(22)=w(25) + w(22)
      w(22)=LL*w(22)
      w(21)=31.D0/9.D0*w(21) + w(22)
      w(21)=w(21)*w(31)
      w(20)=w(28) - w(20) + w(21) + w(29)
      w(20)=TF*w(20)
      w(21)=w(6)*CA*w(46)
      w(20)= - 64.D0/3.D0*w(21) + 16.D0/3.D0*w(20)
      w(20)=w(24)*w(20)
      w(21)=5.D0/3.D0 + z
      w(17)=w(21)*w(17)
      w(17)=44*w(1) + 73.D0/3.D0 + w(17)
      w(21)=7 - w(52)
      w(21)=1.D0/3.D0*w(21) - 3*w(1)
      w(21)=2*w(21) + w(26)
      w(21)=w(3)*w(21)
      w(22)=w(26) + w(34) - 10.D0/3.D0*w(1)
      w(22)=w(22)*w(31)
      w(17)=w(22) + 1.D0/3.D0*w(17) + w(21)
      w(17)=w(17)*w(24)
      w(21)= - 10 + z
      w(19)=1.D0/3.D0*w(21) + w(19)
      w(19)=w(19)*w(31)*w(27)
      w(17)=w(19) + w(17)
      w(17)=w(8)*CF*w(17)
      w(19)=26*w(1) - 38 + 37*z
      w(19)= - 6*w(45) + 1.D0/3.D0*w(19) - 4*w(26)
      w(19)=w(11)*w(19)*w(23)*CF**2
      w(17)=w(19) + w(17)
      w(17)=w(18) + 16*w(17) + w(20)
*
      GQ = TF*w(17)
*
      RETURN
      END
      REAL*8 FUNCTION QGL(z,nf,as,LL)
*     -------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     complete OME AqgQ  
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(61)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL
* 
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 4
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
      w(1)=z**(-1)
      w(2)=z**(-1)
      w(3)=Hr1(-1)
      w(4)=Hr1(0)
      w(5)=Hr1(1)
      w(6)=Hr1(0)
      w(7)=Hr1(1)
      w(8)=Hr2(0,-1)
      w(9)=Hr3(0,0,-1)
      w(10)=Hr3(0,0,-1)
      w(11)=Hr4(0,0,0,-1)
      w(12)=Hr4(0,0,0,-1)
      w(13)=Hr4(0,0,0,1)
      w(14)=Hr4(0,0,0,1)
      w(15)=Hr3(0,0,1)
      w(16)=Hr4(0,0,1,1)
      w(17)=Hr2(0,1)
      w(18)=Hr4(0,1,0,1)
      w(19)=Hr3(0,1,1)
      w(20)=Hr4(0,1,1,1)
      w(21)=z - 1.D0
      w(22)=2.D0*z
      w(23)=w(21)*w(22)
      w(23)=w(23) + 1.D0
      w(24)=w(23)*w(17)
      w(25)=2.D0/3.D0*w(24)
      w(26)=4.D0/3.D0*z
      w(27)=LL*z
      w(28)= - w(26) - w(27)
      w(29)=4.D0*LL
      w(28)=w(28)*w(29)
      w(30)=4.D0*z
      w(31)= - 43.D0/9.D0 - w(30)
      w(31)=z*w(31)
      w(31)= - 2.D0 + w(31)
      w(32)=z2*z
      w(28)= - 56.D0/15.D0*w(32) - w(25) + 1.D0/3.D0*w(31) + w(28)
      w(28)=z2*w(28)
      w(31)=w(23)*w(6)
      w(32)=w(23)*LL
      w(33)=14.D0*z
      w(34)=17.D0 - w(33)
      w(34)=z*w(34)
      w(34)= - 10.D0 + w(34)
      w(34)=1.D0/9.D0*w(34) - w(32)
      w(34)=2.D0*w(34) - w(31)
      w(34)=w(6)*w(17)*w(34)
      w(28)=w(28) + w(34)
      w(34)=z + 1.D0
      w(35)=w(34)*w(22)
      w(35)=w(35) + 1.D0
      w(36)=w(35)*LL
      w(37)=2.D0/3.D0*w(6)
      w(38)= - w(35)*w(37)
      w(39)=7.D0*z
      w(40)= - 4.D0 - w(39)
      w(40)=z*w(40)
      w(38)=w(38) + 2.D0/9.D0*w(40) - w(36)
      w(38)=w(10)*w(38)
      w(40)=w(23)*w(20)
      w(41)=2.D0*LL
      w(42)= - 8.D0 + 31.D0/3.D0*z
      w(42)=z*w(42)
      w(42)= - 1.D0 + w(42)
      w(42)=w(42)*w(41)
      w(43)= - 56.D0 + 205.D0/3.D0*z
      w(43)=z*w(43)
      w(42)=w(42) + 14.D0 + w(43)
      w(42)=LL*w(42)
      w(43)= - 1106.D0 + 3791.D0/3.D0*z
      w(43)=z*w(43)
      w(43)=155.D0 + w(43)
      w(42)=1.D0/3.D0*w(43) + w(42)
      w(42)=w(42)*w(41)
      w(43)= - 44012.D0 + 147643.D0/3.D0*z
      w(43)=z*w(43)
      w(43)= - 17368.D0/3.D0*w(1) + 7217.D0 + w(43)
      w(42)=1.D0/27.D0*w(43) + w(42)
      w(43)= - 67.D0/3.D0 + w(30)
      w(43)=w(43)*w(22)
      w(43)=1.D0 + w(43)
      w(44)= - LL*w(30)
      w(43)=1.D0/3.D0*w(43) + w(44)
      w(44)=4.D0*w(17)
      w(43)=w(43)*w(44)
      w(45)=16.D0*z
      w(46)=w(45)*w(21)
      w(47)= - w(14)*w(46)
      w(48)=w(12)*z*w(34)
      w(28)=8.D0/3.D0*w(40) + 16.D0*w(38) - 8.D0*w(13) + w(47) + 64.D0/
     & 3.D0*w(48) + 32.D0/3.D0*w(11) - 160.D0/9.D0*w(9) + 1.D0/3.D0*
     & w(42) + w(43) + 4.D0*w(28)
      w(28)=CA*w(28)
      w(34)=w(34)*w(39)
      w(34)=w(34) + 5.D0
      w(38)=w(36) + 2.D0/3.D0*w(34)
      w(38)=w(38)*LL
      w(42)=191.D0 + 200.D0*z
      w(42)=w(42)*z
      w(42)=w(42) + 112.D0
      w(38)=w(38) + 1.D0/27.D0*w(42)
      w(42)=w(3)*w(38)
      w(43)= - 11831.D0/3.D0 + 962.D0*z
      w(43)=w(43)*w(22)
      w(43)= - 451.D0/3.D0 + w(43)
      w(47)= - 20.D0/9.D0 + z
      w(47)=z*w(47)
      w(48)=w(30) + 1.D0
      w(49)=LL*w(48)
      w(47)= - 1.D0/9.D0*w(49) - 5.D0/9.D0 + w(47)
      w(47)=w(47)*w(41)
      w(49)= - 19.D0 + w(39)
      w(49)=z*w(49)
      w(49)=5.D0 + 34.D0*w(49)
      w(47)=1.D0/27.D0*w(49) + w(47)
      w(47)=w(47)*w(41)
      w(49)=z2*w(26)
      w(42)=w(49) - 4.D0/3.D0*w(42) + 1.D0/81.D0*w(43) + w(47)
      w(42)=CA*w(42)
      w(34)= - 2.D0/9.D0*w(34) - w(36)
      w(34)=w(3)*w(34)
      w(43)=w(3)*w(35)
      w(47)= - 185.D0 - w(45)
      w(47)=z*w(47)
      w(47)=56.D0 + w(47)
      w(49)=1.D0 - 10.D0/3.D0*z
      w(49)=LL*w(49)
      w(43)=1.D0/3.D0*w(5) - 4.D0/9.D0*w(43) + 1.D0/27.D0*w(47) + w(49)
      w(47)=1.D0/3.D0*w(4)
      w(49)=5.D0/3.D0 - w(30)
      w(49)=w(49)*w(47)
      w(43)=2.D0*w(43) + w(49)
      w(43)=w(4)*w(43)
      w(49)= - 8.D0 - z
      w(49)=z*w(49)
      w(49)= - 8.D0*w(27) + 1.D0 + 16.D0/3.D0*w(49)
      w(49)=LL*w(49)
      w(50)= - 125.D0 - 202.D0/9.D0*z
      w(50)=z*w(50)
      w(50)=205.D0/9.D0 + w(50)
      w(34)=w(43) + 4.D0*w(34) + 2.D0/3.D0*w(50) + w(49)
      w(34)=w(47)*CA*w(34)
      w(34)=w(42) + w(34)
      w(34)=w(4)*w(34)
      w(28)=1.D0/3.D0*w(28) + 2.D0*w(34)
      w(34)=w(39)*w(21)
      w(34)=w(34) + 5.D0
      w(42)=w(32) + 2.D0/9.D0*w(34)
      w(43)=w(42)*CF
      w(49)= - 1.D0 + 14.D0/9.D0*z
      w(49)=z*w(49)
      w(49)= - 1.D0/9.D0*w(31) + w(32) + 10.D0/9.D0 + w(49)
      w(49)=CA*w(49)
      w(49)=w(49) - 1.D0/3.D0*w(43)
      w(49)=w(19)*w(49)
      w(50)=5.D0 - 124.D0/27.D0*z
      w(50)=z*w(50)
      w(51)=LL*w(22)
      w(48)=w(4)*w(48)
      w(48)=8.D0/9.D0*w(48) + w(51) + 4.D0/9.D0 + w(50)
      w(48)=CA*w(48)
      w(50)=w(22) - 1.D0/3.D0
      w(50)=w(50)*z
      w(50)=w(50) - 4.D0/3.D0
      w(51)=w(22) - 1.D0
      w(52)=w(51)*w(4)
      w(53)=2.D0*w(50) + w(52)
      w(54)=4.D0*w(4)
      w(53)=w(53)*w(54)
      w(55)=251.D0/9.D0 - 13.D0*z
      w(55)=z*w(55)
      w(55)= - 154.D0/9.D0 + w(55)
      w(53)=w(53) + 2.D0*w(55) + w(32)
      w(55)=1.D0/3.D0*CF
      w(53)=w(53)*w(55)
      w(55)=2.D0*CF
      w(56)=CA + w(55)
      w(56)=w(2)*w(56)
      w(48)=16.D0/27.D0*w(56) + w(53) + w(48)
      w(48)=z3*w(48)
      w(48)=w(49) + w(48)
      w(49)=w(30)*w(21)
      w(49)=w(49) + 5.D0
      w(53)=2.D0/3.D0*w(49) + w(32)
      w(53)=LL*w(53)
      w(53)=28.D0/27.D0*w(34) + w(53)
      w(53)=w(17)*w(53)
      w(34)=w(32) + 2.D0/3.D0*w(34)
      w(56)=w(34)*LL
      w(57)=w(23)*z2
      w(58)=w(21)*z
      w(59)= - 8.D0 - 13.D0*w(58)
      w(59)= - 2.D0/15.D0*w(57) + 14.D0/27.D0*w(59) - w(56)
      w(59)=z2*w(59)
      w(46)=w(46) + 5.D0
      w(46)=w(32) + 1.D0/9.D0*w(46)
      w(60)= - w(17)*w(46)
      w(61)=1.D0/3.D0*w(6)
      w(24)= - w(61)*w(24)
      w(24)=w(60) + w(24)
      w(24)=w(6)*w(24)
      w(24)=w(24) + w(53) + w(59)
      w(53)= - 113.D0 - 250.D0/3.D0*z
      w(53)=w(53)*w(30)
      w(59)= - 43.D0 + 62.D0/3.D0*z
      w(59)=w(59)*w(22)
      w(59)=53.D0 + w(59)
      w(59)=w(59)*w(41)
      w(53)=w(59) + 749.D0 + w(53)
      w(53)=LL*w(53)
      w(59)= - 5557.D0 - 1778.D0/3.D0*z
      w(59)=w(59)*w(22)
      w(59)=12725.D0 + w(59)
      w(53)=1.D0/3.D0*w(59) + w(53)
      w(53)=w(53)*w(29)
      w(59)= - 422543.D0 - 89762.D0/3.D0*z
      w(59)=w(59)*w(22)
      w(59)=904621.D0 + w(59)
      w(53)=1.D0/27.D0*w(59) + w(53)
      w(59)=w(14)*w(23)
      w(24)=32.D0/3.D0*w(40) - 64.D0/3.D0*w(59) + 1.D0/3.D0*w(53) + 32.D
     & 0*w(24)
      w(40)= - 29.D0 - w(39)
      w(40)=w(40)*w(30)
      w(40)=67.D0 + w(40)
      w(53)=w(51)*LL
      w(40)= - 4.D0/5.D0*w(52) + 1.D0/9.D0*w(40) - 5.D0*w(53)
      w(40)=w(4)*w(40)
      w(45)= - 1291.D0/9.D0 + w(45)
      w(45)=w(45)*w(22)
      w(45)=2341.D0/9.D0 + w(45)
      w(51)= - w(51)*w(41)
      w(52)= - 37.D0 - w(33)
      w(52)=z*w(52)
      w(52)=26.D0 + w(52)
      w(51)=1.D0/3.D0*w(52) + w(51)
      w(51)=w(51)*w(29)
      w(40)=w(40) + 1.D0/3.D0*w(45) + w(51)
      w(40)=w(40)*w(47)
      w(26)= - 1.D0 - w(26)
      w(26)=w(26)*w(22)
      w(26)= - 1.D0/3.D0*w(53) + 3.D0 + w(26)
      w(26)=w(26)*w(29)
      w(45)= - 29.D0 + 68.D0/3.D0*z
      w(45)=z*w(45)
      w(26)=w(26) + 101.D0 + 2.D0/3.D0*w(45)
      w(26)=LL*w(26)
      w(45)= - 5141.D0 - 2152.D0*z
      w(45)=z*w(45)
      w(45)=13876.D0 + w(45)
      w(26)=w(40) + 1.D0/81.D0*w(45) + w(26)
      w(26)=w(4)*w(26)
      w(40)=73.D0 + 80.D0*z
      w(40)=w(40)*w(22)
      w(40)=305.D0 + w(40)
      w(29)= - w(50)*w(29)
      w(29)=1.D0/3.D0*w(40) + w(29)
      w(29)=LL*w(29)
      w(40)=32.D0*z
      w(45)=2393.D0 - w(40)
      w(45)=z*w(45)
      w(45)=4145.D0 + w(45)
      w(29)=1.D0/9.D0*w(45) + w(29)
      w(29)=w(29)*w(41)
      w(45)=2296.D0/3.D0 + 155.D0*z
      w(40)=w(45)*w(40)
      w(40)=154967.D0/3.D0 + w(40)
      w(29)=1.D0/27.D0*w(40) + w(29)
      w(26)=1.D0/3.D0*w(29) + w(26)
      w(26)=w(26)*w(54)
      w(24)=1.D0/3.D0*w(24) + w(26)
      w(24)=CF*w(24)
      w(26)=w(35)*w(61)
      w(29)=11.D0 + w(33)
      w(29)=z*w(29)
      w(29)=10.D0 + w(29)
      w(26)=w(26) + 1.D0/9.D0*w(29) + w(36)
      w(26)=w(6)*w(26)
      w(26)=2.D0*w(26) + w(38)
      w(26)=w(8)*CA*w(26)
      w(29)=w(23)*CA
      w(33)=160.D0/9.D0*w(16) + 64.D0/9.D0*w(18)
      w(33)=w(29)*w(33)
      w(35)=13.D0 - w(41)
      w(35)=LL*w(35)
      w(35)= - 50.D0/3.D0 + w(35)
      w(35)=LL*w(35)
      w(35)=667.D0/27.D0 + w(35)
      w(35)=w(35)*w(55)
      w(36)= - 13.D0 - w(41)
      w(36)=LL*w(36)
      w(36)= - 86.D0/3.D0 + w(36)
      w(36)=CA*LL*w(36)
      w(35)=w(36) + w(35)
      w(35)=w(2)*w(35)
      w(24)=32.D0/27.D0*w(35) + 32.D0/3.D0*w(26) + 4.D0*w(28) + w(24)
     &  + 32.D0*w(48) + w(33)
      w(26)=nf*as**3
      w(24)=w(24)*w(26)
      w(28)=2.D0*z2 - w(44)
      w(28)=w(42)*w(28)
      w(33)=w(32) + w(49)
      w(33)=w(33)*w(41)
      w(35)= - 149.D0 + 158.D0*z
      w(35)=z*w(35)
      w(35)=85.D0 + w(35)
      w(35)=1.D0/3.D0*w(35) + w(33)
      w(35)=LL*w(35)
      w(36)= - 1067.D0 + 1184.D0*z
      w(36)=z*w(36)
      w(36)=487.D0 + w(36)
      w(35)=2.D0/27.D0*w(36) + w(35)
      w(36)=w(58)*w(37)
      w(36)=w(36) + w(42)
      w(36)=w(6)*w(36)
      w(30)=7.D0 - w(30)
      w(30)=z*w(30)
      w(30)=2.D0*w(57) + 1.D0 + w(30)
      w(30)=1.D0/3.D0*w(30) + w(36)
      w(30)=w(6)*w(30)
      w(36)=w(19)*w(23)
      w(28)=4.D0/9.D0*w(2) - 4.D0/3.D0*w(36) + w(30) + 1.D0/3.D0*w(35)
     &  + w(28)
      w(28)=CA*w(28)
      w(30)=149.D0 - 140.D0*z
      w(30)=z*w(30)
      w(30)= - 112.D0 + w(30)
      w(30)=1.D0/3.D0*w(30) - w(33)
      w(30)=LL*w(30)
      w(33)=1139.D0 - 1094.D0*z
      w(33)=w(33)*w(22)
      w(33)= - 1613.D0 + w(33)
      w(30)=1.D0/27.D0*w(33) + w(30)
      w(33)=2.D0/9.D0*w(31) + w(46)
      w(33)=w(6)*w(33)
      w(21)=w(21)*w(27)
      w(27)= - 2.D0 - w(58)
      w(21)=7.D0/27.D0*w(27) + w(21)
      w(21)=4.D0*w(21) + w(33)
      w(21)=w(6)*w(21)
      w(21)=1.D0/3.D0*w(30) + w(21)
      w(21)=CF*w(21)
      w(27)=w(23)*CF
      w(30)= - w(29) + 4.D0/3.D0*w(27)
      w(30)=z3*w(30)
      w(21)=4.D0*w(30) + w(21) + w(28)
      w(28)=4.D0*w(26)
      w(21)=w(21)*w(28)
      w(30)=100.D0 - 109.D0*z
      w(30)=w(30)*w(22)
      w(30)= - 85.D0 + w(30)
      w(25)= - 1.D0/3.D0*w(57) + w(25) + 1.D0/27.D0*w(30) - w(56)
      w(30)=w(23)*w(61)
      w(33)=2.D0*w(42) + w(30)
      w(33)=w(6)*w(33)
      w(25)=2.D0*w(25) + w(33)
      w(25)=CA*w(25)
      w(33)=w(34)*w(41)
      w(34)= - 373.D0 + 364.D0*z
      w(34)=z*w(34)
      w(34)=224.D0 + w(34)
      w(33)=1.D0/27.D0*w(34) + w(33)
      w(33)=CF*w(33)
      w(25)=w(25) + w(33)
      w(25)=w(25)*w(26)
      w(30)= - w(30) + w(42)
      w(30)=CA*w(30)
      w(30)=w(30) - w(43)
      w(28)=w(30)*w(28)
      w(27)= - w(29) + w(27)
      w(30)=w(26)*w(7)
      w(27)=w(27)*w(30)
      w(27)=w(28) + 1.D0/3.D0*w(27)
      w(27)=w(7)*w(27)
      w(25)=2.D0*w(25) + 1.D0/3.D0*w(27)
      w(25)=w(7)*w(25)
      w(21)=w(21) + w(25)
      w(21)=w(7)*w(21)
      w(22)= - 5.D0 + w(22)
      w(22)=w(22)*w(39)
      w(22)=10.D0 + w(22)
      w(22)=w(31) + 1.D0/9.D0*w(22) + w(32)
      w(22)=CA*w(22)
      w(23)=w(23)*w(37)
      w(23)=w(23) + w(46)
      w(23)=CF*w(23)
      w(22)=w(22) + w(23)
      w(22)=w(22)*w(26)
      w(23)=w(29)*w(30)
      w(22)=w(22) - 2.D0/3.D0*w(23)
      w(22)=w(15)*w(22)
      w(21)=32.D0/3.D0*w(22) + w(24) + 4.D0/3.D0*w(21)
*
      QGL = TF**2*w(21)
*
      RETURN
      END
      REAL*8 FUNCTION PSL(z,nf,as,LL)
*     -------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     OME AqqQPS  without aPS3 
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 w(30)
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5
      real*8 Hr1,Hr2,Hr3,Hr4,Hr5
      real*8 Hi1,Hi2,Hi3,Hi4,Hi5
      dimension Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),
     $          Hc4(-1:1,-1:1,-1:1,-1:1),
     $          Hc5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),
     $          Hr4(-1:1,-1:1,-1:1,-1:1),
     $          Hr5(-1:1,-1:1,-1:1,-1:1,-1:1)
      dimension Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),
     $          Hi4(-1:1,-1:1,-1:1,-1:1),
     $          Hi5(-1:1,-1:1,-1:1,-1:1,-1:1)
      INTEGER nw
      REAL*8 z,CF,CA,TF,nf,as,z2,z3,LL
* 
      z2 = 1.6449340668482264365D0
      z3 = 1.2020569031595942854D0 
      CF=4.0D0/3.0D0
      TF=1.0D0/2.0D0
      CA=3.0D0
*
      nw = 4
      call hplog5(z,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,-1,1)
*
      w(1)=z**(-1)
      w(2)=z**(-1)
      w(3)=Hr1(0)
      w(4)=Hr1(0)
      w(5)=Hr4(0,0,0,1)
      w(6)=Hr3(0,0,1)
      w(7)=Hr4(0,0,1,1)
      w(8)=Hr4(0,0,1,1)
      w(9)=Hr2(0,1)
      w(10)=Hr3(0,1,1)
      w(11)=Hr4(0,1,1,1)
      w(12)=Hr1(1)
      w(13)= - 7.D0 + 4.D0*z
      w(13)=w(13)*z
      w(13)=w(13) - 13.D0
      w(14)=z + 1.D0
      w(15)=w(14)*LL
      w(15)= - w(15) + 1.D0/3.D0*w(13)
      w(15)=w(15)*LL
      w(16)=19.D0*z
      w(17)=w(16) - 131.D0/3.D0
      w(17)=w(17)*z
      w(17)=w(17) - 203.D0/3.D0
      w(15)=w(15) + 1.D0/9.D0*w(17)
      w(17)=w(9)*w(15)
      w(18)=2.D0*LL
      w(19)=w(18)*w(14)
      w(20)=2.D0*z
      w(21)= - 25.D0/3.D0 + w(20)
      w(21)=z*w(21)
      w(21)= - 34.D0/3.D0 + w(21)
      w(21)=1.D0/3.D0*w(21) - w(19)
      w(21)=w(10)*w(21)
      w(17)=w(17) - w(21)
      w(21)=40.D0/3.D0*w(11) + 16.D0/3.D0*w(5)
      w(21)=w(14)*w(21)
      w(22)=w(20) - 7.D0
      w(22)=w(22)*z
      w(22)=w(22) + 5.D0
      w(23)=1.D0 + 4.D0/3.D0*z
      w(23)=w(23)*z
      w(23)=w(23) - 1.D0 - 4.D0/3.D0*w(2)
      w(24)=w(23)*LL
      w(25)=1.D0/3.D0*w(24) - w(22)
      w(25)=LL*w(25)
      w(26)=197.D0 - 80.D0/3.D0*z
      w(26)=z*w(26)
      w(26)= - 136.D0/3.D0*w(2) - 125.D0 + w(26)
      w(25)=1.D0/9.D0*w(26) + w(25)
      w(25)=LL*w(25)
      w(26)=1280.D0 - 139.D0/3.D0*z
      w(26)=z*w(26)
      w(26)= - 983.D0 + w(26)
      w(26)=2.D0*w(26) - 2351.D0*w(3)
      w(25)= - 8.D0/3.D0*w(7) - 1504.D0/243.D0*w(1) + 1.D0/81.D0*w(26)
     &  + w(25)
      w(26)=w(19) + 29.D0/9.D0*w(14)
      w(27)=w(6)*w(26)
      w(28)=w(8)*z
      w(17)= - 32.D0/3.D0*w(28) + 8.D0*w(27) + 4.D0*w(25) + w(21) - 8.D0
     & *w(17)
      w(17)=nf*w(17)
      w(21)= - 29.D0/3.D0*w(14) - w(19)
      w(21)=LL*w(21)
      w(21)= - 332.D0/27.D0*w(14) + w(21)
      w(25)=2.D0*nf
      w(21)=w(21)*w(25)
      w(27)= - w(26)*w(25)
      w(28)=w(14)*w(4)*nf
      w(29)=1.D0/3.D0*w(28)
      w(27)=w(27) - w(29)
      w(30)=1.D0/3.D0*w(4)
      w(27)=w(27)*w(30)
      w(21)=w(21) + w(27)
      w(21)=w(4)*w(21)
      w(13)= - w(19) + w(13)
      w(13)=LL*w(13)
      w(16)= - 16.D0 + w(16)
      w(16)=z*w(16)
      w(16)= - 40.D0 + w(16)
      w(13)=4.D0/3.D0*w(16) + w(13)
      w(13)=LL*w(13)
      w(16)= - 1577.D0/3.D0 + 220.D0*z
      w(16)=z*w(16)
      w(13)=1.D0/9.D0*w(16) + w(13)
      w(13)=nf*w(13)
      w(13)=4.D0/3.D0*w(13) + w(21)
      w(13)=w(4)*w(13)
      w(16)=2.D0*w(22) - w(24)
      w(16)=LL*w(16)
      w(19)= - 203.D0 - 10.D0*z
      w(19)=z*w(19)
      w(19)=64.D0*w(2) + 149.D0 + w(19)
      w(16)=2.D0/27.D0*w(19) + w(16)
      w(16)=w(16)*w(25)
      w(18)= - 5.D0/9.D0*w(12) + w(18)
      w(18)=w(23)*w(18)
      w(19)=29.D0 - 34.D0/9.D0*z
      w(19)=z*w(19)
      w(19)= - 20.D0/9.D0*w(2) - 23.D0 + w(19)
      w(18)=1.D0/3.D0*w(19) + w(18)
      w(18)=w(12)*nf*w(18)
      w(16)=w(16) + w(18)
      w(16)=w(12)*w(16)
      w(14)=z2*w(14)
      w(14)= - 4.D0/5.D0*w(14) + w(15)
      w(14)=nf*w(14)
      w(15)= - nf*w(26)
      w(15)=w(15) - w(29)
      w(15)=w(4)*w(15)
      w(14)=w(15) + w(14)
      w(14)=z2*w(14)
      w(15)= - 8.D0 - 11.D0*z
      w(15)=w(15)*w(20)
      w(15)=16.D0*w(2) + 17.D0 + w(15)
      w(18)=1.D0/3.D0*nf
      w(15)=w(15)*w(18)
      w(15)=w(15) + 10.D0*w(28)
      w(15)=z3*w(15)
      w(13)=8.D0*w(14) + 2.D0*w(16) + 8.D0/3.D0*w(15) + w(13) + w(17)
*
      PSL = 8.D0/3.D0*TF**2*CF*as**3*w(13)
*
      RETURN
      END
      SUBROUTINE PARAM(INLO,IQ2)
*     --------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Effective parameterization of PDFs; parameter lists 
*     Only the respective choice is transferred by COMMON/PARAME/.   
*     ------------------------------------------------------------------------
*
*     EFFECTIVE MELLIN-SPACE PARAMETERS OF THE PARTON DENSITIES      
*
*******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMMON/PARAME/ XNNS,XANS,XBNS,XCNS,XDNS,
     &               XNSI,XASI,XBSI,XCSI,XDSI,
     &               XNGL,XAGL,XBGL,XCGL,XDGL
*     
      IF(INLO.EQ.1.AND.IQ2.EQ.1) GOTO 1
      IF(INLO.EQ.1.AND.IQ2.EQ.2) GOTO 2
      IF(INLO.EQ.1.AND.IQ2.EQ.3) GOTO 3
      IF(INLO.EQ.2.AND.IQ2.EQ.1) GOTO 4
      IF(INLO.EQ.2.AND.IQ2.EQ.2) GOTO 5
      IF(INLO.EQ.2.AND.IQ2.EQ.3) GOTO 6
      GOTO 101
*
*******************************************************************
1     CONTINUE
*
*     NLO: Q^2 = 30 GeV^2
*
*------------------------------------------------------------------
*
      XNNS =   0.10852D+00    
      XANS =   0.32896D+00    
      XBNS =   0.30695D+01    
      XCNS =   0.87302D+01    
      XDNS =   0.21965D+00    
*------------------------------------------------------------------
      XNSI =   0.60000D+00    
      XASI =  -0.30000D+00    
      XBSI =   0.35000D+01    
      XCSI =   0.50000D+01    
      XDSI =   0.80000D+00    
*------------------------------------------------------------------
      XNGL =   0.11518D+01    
      XAGL =  -0.32230D+00
      XBGL =   0.60445D+01
      XCGL =   0.96178D-01
      XDGL =   0.42125D-03
*------------------------------------------------------------------
      RETURN
2     CONTINUE
*
*     NLO: Q^2 = 100 GeV^2
*
*------------------------------------------------------------------
      XNNS =   0.91250D-01    
      XANS =   0.34190D+00  
      XBNS =   0.32150D+01
      XCNS =   0.97489D+01
      XDNS =   0.16524D+00
*------------------------------------------------------------------
      XNSI =   0.66225D+00  
      XASI =  -0.32101D+00     
      XBSI =   0.38873D+01    
      XCSI =   0.52523D+01      
      XDSI =   0.94456D+00
*------------------------------------------------------------------
      XNGL =   0.23605D+01
      XAGL =  -0.30918D+00
      XBGL =   0.29238D+01 
      XCGL =  -0.11007D+01      
      XDGL =   0.24354D+00
*------------------------------------------------------------------
      RETURN
3     CONTINUE
*
*     NLO: Q^2 = 10000 GeV^2
*
*------------------------------------------------------------------
      XNNS =   0.60626D-01      
      XANS =   0.45319D+00
      XBNS =   0.37446D+01   
      XCNS =   0.13034D+02
      XDNS =  -0.25986D-01   
*------------------------------------------------------------------
      XNSI =   0.69551D+00 
      XASI =  -0.39319D+00
      XBSI =   0.62707D+01
      XCSI =   0.20543D+02
      XDSI =   0.15716D+01
*------------------------------------------------------------------
      XNGL =   0.23605D+01
      XAGL =  -0.30918D+00
      XBGL =   0.29238D+01  
      XCGL =  -0.11007D+01
      XDGL =   0.24354D+00
*------------------------------------------------------------------
      RETURN
*******************************************************************
4     CONTINUE
*
*     NNLO: Q^2 = 30 GeV^2
*
*------------------------------------------------------------------
*
      XNNS =   0.10852D+00    
      XANS =   0.32896D+00    
      XBNS =   0.30695D+01    
      XCNS =   0.87302D+01    
      XDNS =   0.21965D+00    
*------------------------------------------------------------------
      XNSI =   0.60000D+00    
      XASI =  -0.30000D+00    
      XBSI =   0.35000D+01    
      XCSI =   0.50000D+01    
      XDSI =   0.80000D+00    
*------------------------------------------------------------------
      XNGL =   0.11518D+01    
      XAGL =  -0.32230D+00
      XBGL =   0.60445D+01
      XCGL =   0.96178D-01
      XDGL =   0.42125D-03
*------------------------------------------------------------------
      RETURN
5     CONTINUE
*
*     NNLO: Q^2 = 100 GeV^2
*
*------------------------------------------------------------------
      XNNS =   0.90012D-01   
      XANS =   0.37334D+00
      XBNS =   0.32021D+01
      XCNS =   0.97870D+01
      XDNS =   0.12447D+00
*------------------------------------------------------------------
      XNSI =   0.65216D+00
      XASI =  -0.32286D+00  
      XBSI =   0.38537D+01
      XCSI =   0.51846D+01
      XDSI =   0.93115D+00
*------------------------------------------------------------------
      XNGL =   0.24791D+01    
      XAGL =  -0.30370D+00
      XBGL =   0.28748D+01   
      XCGL =  -0.10991D+01
      XDGL =   0.23315D+00   
*------------------------------------------------------------------
      RETURN
6     CONTINUE
*
*     NNLO: Q^2 = 10000 GeV^2
*
*------------------------------------------------------------------
      XNNS =   0.63573D-01 
      XANS =   0.26936D+00
      XBNS =   0.37529D+01 
      XCNS =   0.12672D+02   
      XDNS =   0.18461D+00   
*------------------------------------------------------------------
      XNSI =   0.69074D+00
      XASI =  -0.39337D+00
      XBSI =   0.62209D+01
      XCSI =   0.20164D+02
      XDSI =   0.15624D+01
*------------------------------------------------------------------
      XNGL =   0.19093D+02  
      XAGL =  -0.31917D+00 
      XBGL =   0.35054D+01
      XCGL =  -0.10375D+01
      XDGL =   0.31523D-01   
*------------------------------------------------------------------
      RETURN
*
101   WRITE(6,*) '**** STOP: INLO,IQ2=', INLO,IQ2
      STOP
*
      END
      REAL*8 FUNCTION NSDIS(z)
*     ------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     NS parametreization: z NS(z)
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 z
      REAL*8 XNNS,XANS,XBNS,XCNS,XDNS
      REAL*8 XNSI,XASI,XBSI,XCSI,XDSI
      REAL*8 XNGL,XAGL,XBGL,XCGL,XDGL 
      COMMON/PARAME/ XNNS,XANS,XBNS,XCNS,XDNS,
     &               XNSI,XASI,XBSI,XCSI,XDSI,
     &               XNGL,XAGL,XBGL,XCGL,XDGL
*     
      NSDIS = XNNS*z**XANS*(1.0D0-z)**XBNS*(1.0D0+XCNS*z**XDNS)
*
      RETURN
      END
      REAL*8 FUNCTION SIGDIS(z)
*     ------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Singlet parameterization: z Sigma(z)
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 z,CDIS
      REAL*8 XNNS,XANS,XBNS,XCNS,XDNS
      REAL*8 XNSI,XASI,XBSI,XCSI,XDSI
      REAL*8 XNGL,XAGL,XBGL,XCGL,XDGL 
      INTEGER INLO,IQ2
      COMMON/PARAME/ XNNS,XANS,XBNS,XCNS,XDNS,
     &               XNSI,XASI,XBSI,XCSI,XDSI,
     &               XNGL,XAGL,XBGL,XCGL,XDGL
      COMMON/PILOT/INLO,IQ2
      IF(IQ2.EQ.1.AND.INLO.EQ.2)  CDIS=1.0085561697363055D0
      IF(IQ2.EQ.2.AND.INLO.EQ.2)  CDIS=1.0003257386718929D0
      IF(IQ2.EQ.3.AND.INLO.EQ.2)  CDIS=1.0524360309677239D0 
*     
      SIGDIS = XNSI*z**XASI*(1.0D0-z)**XBSI*(1.0D0+XCSI*z**XDSI)
     &         *CDIS
*
      RETURN
      END
      REAL*8 FUNCTION GLUDIS(z)
*     -------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Gluon parameterization: z G(z)
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      INTEGER INLO,IQ2
      REAL*8 z,CDIS 
      REAL*8 XNNS,XANS,XBNS,XCNS,XDNS
      REAL*8 XNSI,XASI,XBSI,XCSI,XDSI
      REAL*8 XNGL,XAGL,XBGL,XCGL,XDGL 
      COMMON/PARAME/ XNNS,XANS,XBNS,XCNS,XDNS,
     &               XNSI,XASI,XBSI,XCSI,XDSI,
     &               XNGL,XAGL,XBGL,XCGL,XDGL
      COMMON/PILOT/INLO,IQ2
*     
      IF(IQ2.EQ.1.AND.INLO.EQ.2)  CDIS=1.0085561697363055D0
      IF(IQ2.EQ.2.AND.INLO.EQ.2)  CDIS=1.0003257386718929D0
      IF(IQ2.EQ.3.AND.INLO.EQ.2)  CDIS=1.0524360309677239D0 
      GLUDIS = XNGL*z**XAGL*(1.0D0-z)**XBGL*(1.0D0+XCGL*z**XDGL)
     &       *CDIS
*
      RETURN
      END
      REAL*8 FUNCTION MOMNS(N,x,nf,as,LL)
*     -----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 16, 2024
*     Mellin moments for NS
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND
      REAL*8 MOMNSDEL,MOMNSPLU,MOMNSREG
!      REAL*8 NSDEL,ANSDEL,aMOMNSPLU,aMOMNSREG
!      REAL*8 aMOMNSPLU,aMOMNSREG
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B,DEL,A1,B1
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL DAIND
!      EXTERNAL aMOMNSPLU,aMOMNSREG
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0
      DEL=1.0D-14
      A1=A+DEL
      B1=1.0D0-DEL
      KEY=1
      MAX=10000
      EPS=1.0D-10 ! GPS mod for use with dgauss, previously was 1.0D-18
*
      MOMNSDEL=NSDEL(x,nf,as,LL)+ANSDEL(x,nf,as,LL)
*      WRITE(6,*) 'N,x,nf,as,LL=',N,x,nf,as,LL 
*
      MOMNSPLU=DAIND(A,B,aMOMNSPLU,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
      MOMNSREG=DAIND(A,B,aMOMNSREG,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      MOMNS=MOMNSDEL+MOMNSPLU+MOMNSREG 
*
      RETURN
      END
      REAL*8 FUNCTION aMOMNSPLU(y)
*     ----------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMPSL
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z,y
!      REAL*8 NSPLU,ANSPLU1,ANSPLU2
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL NSPLU,ANSPLU
*
      aMOMNSPLU = (y**(NN-1.0D0)-1.0D0)*(NSPLU(y,nf,as,LL)
     &       +   ANSPLU1(y,nf,as,LL) + ANSPLU2(y,nf,as,LL))
* 
      RETURN
      END
      REAL*8 FUNCTION aMOMNSREG(y)
*     ----------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMPSL
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z,y
!      REAL*8 NSREG,ANSREG
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL NSREG,ANSREG
*
      aMOMNSREG = y**(NN-1.0D0)*(NSREG(y,nf,as,LL)
     &            +ANSREG(y,nf,as,LL))
* 
      RETURN
      END
      REAL*8 FUNCTION MOMGG(N,x,nf,as,LL)
*     -----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Mellin moments for AGG
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND
      REAL*8 MOMGGDEL,MOMGGPLU,MOMGGREG
!      REAL*8 GGDEL,AGGDEL,aMOMGGPLU,aMOMGGREG
!      REAL*8 aMOMGGPLU,aMOMGGREG
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL DAIND
!      EXTERNAL aMOMGGPLU,aMOMGGREG
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0
      KEY=2
      MAX=10000
      EPS=1.0D-10 ! GPS mod for use with dgauss, previously was 1.0D-18
*
      MOMGGDEL=GGDEL(x,nf,as,LL)+AGGDEL(x,nf,as,LL)
      MOMGGPLU=DAIND(A,B,aMOMGGPLU,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
      ! GPS redefine EPS to be even less accurate for the following term,
      ! which is especially difficult to integrate when LL is 10, as
      ! in the tests of the VFNS3sub routine
      EPS=1.0D-6
      MOMGGREG=DAIND(A,B,aMOMGGREG,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      MOMGG=MOMGGDEL+MOMGGPLU+MOMGGREG
*
      RETURN
      END
      REAL*8 FUNCTION aMOMGGPLU(y)
*     ----------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMGG
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z,y
!      REAL*8 GGPLU,AGGPLU
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL GGPLU,AGGPLU
*
      aMOMGGPLU = (y**(NN-1.0D0)-1.0D0)*(
     &          GGPLU(y,nf,as,LL)+AGGPLU(y,nf,as,LL))
* 
      RETURN
      END
      REAL*8 FUNCTION aMOMGGREG(y)
*     ----------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMGG
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z,y
      REAL*8 xx,nnf,aas,LLL
!      REAL*8 GGREG,AGGREG0,AGGREG1,AGG
      REAL*8 AGG
      INTEGER NN
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
!      EXTERNAL GGREG,AGGREG0,AGGREG1
      nf=nnf
      LL=LLL
      as=aas
*
      IF(y.LT.0.5D0) AGG=AGGREG0(y,nf,as,LL)
      IF(y.GE.0.5D0) AGG=AGGREG1(y,nf,as,LL)
*
*      LL=10.0D0 !!
*      WRITE(6,*) '>>> LL=',LL
      aMOMGGREG = y**(NN-1.0D0)*(GGREG(y,nf,as,LL)+AGG)
*
      RETURN
      END
      REAL*8 FUNCTION MOMQG(N,x,nf,as,LL)
*     -----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 15, 2024
*     Mellin moments for QG
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND,MQG
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL MQG,DAIND
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0
      KEY=1
      MAX=2000
      EPS=1.0D-8 ! GPS mod for use with dgauss, previously was 1.0D-18
*
      MOMQG=DAIND(A,B,MQG,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      RETURN
      END
      REAL*8 FUNCTION MQG(z)
*     ----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMQG
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
!      REAL*8 x,nf,as,LL,z,QG,aQg3,AQG3NF,RED0,RED1,RED
      REAL*8 x,nf,as,LL,z,RED
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL QG,aQg3,AQG3NF,RED0,RED1
*
*      MQG = z**(NN-1.0D0)*(QG(z,nf,as,LL)+aQg3(z)*as**3)
      IF(z.LT.0.5D0) RED=RED0(z)
      IF(z.GE.0.5D0) RED=RED1(z)
      MQG = z**(NN-1.0D0)*(QG(z,nf,as,LL)
     &                   +(aQg3(z)/2+RED+nf*AQG3NF(z))*as**3)
* 
      RETURN
      END
      REAL*8 FUNCTION MOMPSL(N,x,nf,as,LL)
*     ------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Mellin moments for PSL
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND,MPSL
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL MPSL,DAIND
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0
      KEY=1
      MAX=10000
      EPS=1.0D-10 ! GPS mod for use with dgauss, previously was 1.0D-18
*      WRITE(6,*) 'EPS=',EPS
*
      MOMPSL=DAIND(A,B,MPSL,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      RETURN
      END
      REAL*8 FUNCTION MPSL(z)
*     -----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMPSL
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z!,PSL
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL PSL
*
      MPSL = z**(NN-1.0D0)*PSL(z,nf,as,LL)
* 
      RETURN
      END
      REAL*8 FUNCTION MPS(z)
*     ----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMPS
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z!,PS
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL PS
*
      MPS = z**(NN-1.0D0)*PS(z,nf,as,LL)
* 
      RETURN
      END
      REAL*8 FUNCTION MOMPS(N,x,nf,as,LL)
*     ------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Mellin moments for PS
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND!,MPS
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL MPS,DAIND
!      EXTERNAL DAIND
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0
      KEY=1
      MAX=10000
      EPS=1.0D-10 ! GPS mod for use with dgauss, previously was 1.0D-18
*      WRITE(6,*) 'EPS=',EPS
*
      MOMPS=DAIND(A,B,MPS,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      RETURN
      END
      REAL*8 FUNCTION MOMQGL(N,x,nf,as,LL)
*     ------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Mellin moments for QGL
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND,MQGL
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL MQGL,DAIND
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0
      KEY=1
      MAX=10000
      EPS=1.0D-10 ! GPS mod for use with dgauss, previously was 1.0D-18
*
      MOMQGL=DAIND(A,B,MQGL,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      RETURN
      END
      REAL*8 FUNCTION MQGL(z)
*     -----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMQGL
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z!,QGL
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL QGL
*
      MQGL = z**(NN-1.0D0)*QGL(z,nf,as,LL)
* 
      RETURN
      END
      REAL*8 FUNCTION MOMGQ(N,x,nf,as,LL)
*     -----------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Mellin moments for GQ
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND,MGQ
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL MGQ,DAIND
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0-1.0D-14
      KEY=1
      MAX=10000
      EPS=1.0D-10 ! GPS mod for use with dgauss, previously was 1.0D0-18
*
      MOMGQ=DAIND(A,B,MGQ,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      RETURN
      END
      REAL*8 FUNCTION MGQ(z)
*     ----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMGQ
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z!,GQ
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL GQ
*
      MGQ = z**(NN-1.0D0)*GQ(z,nf,as,LL)
* 
      RETURN
      END
      REAL*8 FUNCTION MOMAPS(N,x,nf,as,LL)
*     ------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Mellin moments for APS
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND,AMPS
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL AMPS,DAIND
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0
      KEY=1
      MAX=10000
      EPS=1.0D-10 ! GPS mod for use with dgauss, previously was 1.0D0-18
*
      MOMAPS=DAIND(A,B,AMPS,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      RETURN
      END
      REAL*8 FUNCTION AMPS(z)
*     -----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     integrand for MOMAPS
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z!,APS1,APS2
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL APS1,APS2
*
      AMPS = z**(NN-1.0D0)*(APS1(z,nf,as,LL)+APS2(z,nf,as,LL))
* 
      RETURN
      END
      REAL*8 FUNCTION MOMAGQ(N,x,nf,as,LL)
*     ------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     Mellin moments for AGQ
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,TA!,DAIND,AMGQ
      REAL*8 xx,nnf,aas,LLL
      REAL*8 EPS,EST,A,B
      INTEGER KEY,MAX,KOUNT,N,NN
!      EXTERNAL AMGQ,DAIND
      COMMON/VAR/ xx,nnf,aas,LLL
      COMMON/MO/ NN
*
      NN=N
      xx=x
      nnf=nf
      aas=as
      LLL=LL
      A=0.0D0
      B=1.0D0-1.0D-14
      KEY=1
      MAX=10000
      EPS=1.0D-10 ! GPS mod for use with dgauss, previously was 1.0D0-18
*
      MOMAGQ=DAIND(A,B,AMGQ,EPS,KEY,MAX,KOUNT,EST)
      WRITE(6,*) 'KOUNT,EST=',KOUNT,EST
* 
      RETURN
      END
      REAL*8 FUNCTION AMGQ(z)
*     -----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     interand for MOMAGQ
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 x,nf,as,LL,z!,AGQ
      INTEGER NN
      COMMON/VAR/ x,nf,as,LL
      COMMON/MO/ NN
!      EXTERNAL AGQ
*
      AMGQ = z**(NN-1.0D0)*AGQ(z,nf,as,LL)
* 
      RETURN
      END
      REAL*8 FUNCTION MOMGGF(A1,A2,A3,A4,A5)
*     --------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 19, 2024
*     Fixed Mellin moments for AGG (2,4,6,8,10): test case
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 Pi,ln2,z3,li4half
      INTEGER A1,A2,A3,A4,A5
      REAL*8 w(9)
*
      Pi =  3.1415926535897932385D0
      ln2 = 0.69314718055994530942D0
      z3 = 1.2020569031595942854D0 
      li4half = 0.51747906167389938633D0
*
      w(1)=Pi**2
      w(2)= - 4 - 476117.D0/72000.D0*w(1)
      w(2)=w(2)*w(1)
      w(3)=ln2**2
      w(3)=w(3) - w(1)
      w(3)=w(3)*ln2
      w(4)=8*w(1)
      w(5)= - w(4) + 4709.D0/200.D0*w(3)
      w(6)=1.D0/3.D0*ln2
      w(5)=w(5)*w(6)
      w(2)=4709.D0/25.D0*li4half - 112996108487.D0/41472000.D0*z3 + 
     & w(5) + 1836968368007939.D0/69984000000.D0 + w(2)
      w(2)=a2*w(2)
      w(5)= - 4 - 6929.D0/3240.D0*w(1)
      w(5)=w(5)*w(1)
      w(7)= - w(4) + 83.D0/9.D0*w(3)
      w(7)=w(7)*w(6)
      w(5)=664.D0/9.D0*li4half - 757289.D0/648.D0*z3 + w(7) + 85678351.D
     & 0/23328.D0 + w(5)
      w(5)=a1*w(5)
      w(7)= - 4 - 1623421.D0/198450.D0*w(1)
      w(7)=w(7)*w(1)
      w(8)= - 4*w(1) + 31739.D0/2205.D0*w(3)
      w(8)=ln2*w(8)
      w(7)=507824.D0/2205.D0*li4half - 138256128195967.D0/45519667200.D0
     & *z3 + 2.D0/3.D0*w(8) + 286361615685032825647.D0/8782450790400000.
     & D0 + w(7)
      w(7)=a3*w(7)
      w(8)= - 4 - 300520733.D0/32659200.D0*w(1)
      w(8)=w(8)*w(1)
      w(9)= - w(4) + 2927861.D0/90720.D0*w(3)
      w(9)=w(9)*w(6)
      w(8)=2927861.D0/11340.D0*li4half - 527047428947803.D0/
     & 175575859200.D0*z3 + w(9) + 2058576456832261300537.D0/
     & 57164344876800000.D0 + w(8)
      w(8)=a4*w(8)
      w(9)= - 4 - 456140737.D0/45738000.D0*w(1)
      w(1)=w(9)*w(1)
      w(3)= - w(4) + 4437799.D0/127050.D0*w(3)
      w(3)=w(3)*w(6)
      w(1)=17751196.D0/63525.D0*li4half - 69798388147227963749.D0/
     & 24927089983488000.D0*z3 + w(3) + 19803154527719672837313969830143
     & .D0/518503034748358066176000000.D0 + w(1)
      w(1)=a5*w(1)
      w(1)=w(1) + w(8) + w(7) + w(2) + w(5)
*
      MOMGGF = 1.0D0 + 1.D0/1125.D0*w(1)
*
      RETURN
      END
      REAL*8 FUNCTION MOMGQF(A1,A2,A3,A4,A5)
*     --------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 19, 2024
*     Fixed Mellin moments for AGQ (2,4,6,8,10): test case
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 Pi,ln2,z3,li4half
      INTEGER A1,A2,A3,A4,A5
      REAL*8 w(5)
*
      Pi =  3.1415926535897932385D0
      ln2 = 0.69314718055994530942D0
      z3 = 1.2020569031595942854D0 
      li4half = 0.51747906167389938633D0
*
      w(1)=Pi**2
      w(2)=ln2**2
      w(1)=w(1) - w(2)
      w(1)=w(1)*w(2)
      w(2)=Pi**4
      w(3)=125309332.D0/2205.D0*z3 + 352.D0/3.D0*w(1) + 4532.D0/45.D0*
     & w(2) - 25344699947995223.D0/163364040000.D0 - 2816*li4half
      w(3)=a3*w(3)
      w(4)=440934308.D0/3025.D0*z3 + 896.D0/3.D0*w(1) + 11536.D0/45.D0*
     & w(2) - 40511005719157764003079.D0/173967779795040000.D0 - 7168*
     & li4half
      w(4)=a5*w(4)
      w(5)=6353843.D0/450.D0*z3 + 88.D0/3.D0*w(1) + 1133.D0/45.D0*w(2)
     &  - 913891259981.D0/16200000.D0 - 704*li4half
      w(5)=a2*w(5)
      w(3)=w(5) + 1.D0/7.D0*w(3) + 1.D0/33.D0*w(4)
      w(4)=83462.D0/9.D0*z3 + 64.D0/3.D0*w(1) + 824.D0/45.D0*w(2) - 
     & 32835733.D0/486.D0 - 512*li4half
      w(4)=a1*w(4)
      w(1)=136071571.D0/5670.D0*z3 + 148.D0/3.D0*w(1) + 3811.D0/90.D0*
     & w(2) - 84740936431197853.D0/1728324864000.D0 - 1184*li4half
      w(1)=a4*w(1)
      w(1)=1.D0/21.D0*w(1) + 1.D0/5.D0*w(3) + w(4)
*
      MOMGQF = 1.D0/10125.D0*w(1)
*
      RETURN
      END
      REAL*8 FUNCTION MOMPSLF(a1,a2,a3,a4,a5)
*     ---------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 19, 2024
*     Fixed Mellin moments for APSL (2,4,6,8,10): test case
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 Pi,ln2,z3,li4half,as,nf,LL
      INTEGER a1,a2,a3,a4,a5
      REAL*8 w(3)
      COMMON/VAR1/ nf,as,LL
*
      Pi =  3.1415926535897932385D0
      ln2 = 0.69314718055994530942D0
      z3 = 1.2020569031595942854D0 
      li4half = 0.51747906167389938633D0
*
      w(1)= - 15857574557.D0/200000.D0*a2 - 7777476991657459.D0/
     & 799014273750.D0*a5 - 531041*a1 - 147811731914839.D0/9219840000.D0
     & *a4
      w(2)=484.D0/25.D0*a2 + 25088.D0/9075.D0*a5 + 128*a1 + 1369.D0/315.
     & D0*a4
      w(2)=z3*w(2)
      w(3)= - 344048264191.D0/1234800.D0 + 1936*z3
      w(3)=a3*w(3)
      w(1)=1.D0/245.D0*w(3) + 1.D0/27.D0*w(1) + w(2)
*
      MOMPSLF =  1.D0/10125.D0*w(1)
*
      RETURN
      END
      REAL*8 FUNCTION MOMQGLF(a1,a2,a3,a4,a5)
*     ---------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 19, 2024
*     Fixed Mellin moments for AQGL (2,4,6,8,10): test case
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 Pi,ln2,z3,li4half,as
      INTEGER a1,a2,a3,a4,a5
      REAL*8 w(2)
      COMMON/as/as
*
      Pi =  3.1415926535897932385D0
      ln2 = 0.69314718055994530942D0
      z3 = 1.2020569031595942854D0 
      li4half = 0.51747906167389938633D0
*
      w(1)=416589745654819474587919.D0/5315682160404000000.D0*a5 +
     & 9262656068352038659.D0/100818950400000.D0*a4 + 7399018176313.D0/
     & 54000000.D0*a2 + 901919.D0/9.D0*a1 + 2085272800153061.D0/
     & 18823840000.D0*a3
      w(2)= - 17723218.D0/831875.D0*a5 - 9099743.D0/378000.D0*a4 -
     & 39809.D0/1250.D0*a2 - 14*a1 - 949003.D0/34300.D0*a3
      w(2)=z3*w(2)
      w(1)=1.D0/20.D0*w(1) + w(2)
*
      MOMQGLF =  1.D0/2025.D0*w(1)
*
      RETURN
      END
      REAL*8 FUNCTION MOMPSF(a1,a2,a3,a4,a5)
*     --------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 19, 2024
*     Fixed Mellin moments for APS (2,4,6,8,10): test case
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 Pi,ln2,z3,li4half,as
      INTEGER a1,a2,a3,a4,a5
      REAL*8 w(6)
      COMMON/as/as
*
      Pi =  3.1415926535897932385D0
      ln2 = 0.69314718055994530942D0
      z3 = 1.2020569031595942854D0 
      li4half = 0.51747906167389938633D0
*
      w(1)=Pi**2
      w(2)=ln2**2
      w(1)=w(1) - w(2)
      w(1)=w(1)*w(2)
      w(2)=Pi**4
      w(3)= - 44288.D0/9.D0*z3 - 32.D0/3.D0*w(1) - 412.D0/45.D0*w(2) +
     & 10582393.D0/162.D0 + 256*li4half
      w(3)=a1*w(3)
      w(4)= - 11.D0/3.D0*w(1) - 1133.D0/360.D0*w(2) + 354488673629.D0/
     & 9720000.D0 + 88*li4half
      w(4)=11*w(4) - 1678937.D0/144.D0*z3
      w(4)=a2*w(4)
      w(5)=18155539.D0/2016.D0*z3 - 484.D0/3.D0*w(1) - 12463.D0/90.D0*
     & w(2) + 33664444401825143.D0/18670176000.D0 + 3872*li4half
      w(5)=a3*w(5)
      w(6)=59857684399.D0/967680.D0*z3 - 1369.D0/12.D0*w(1) - 141007.D0/
     & 1440.D0*w(2) + 216560453347331774509.D0/161310320640000.D0 +
     & 2738*li4half
      w(6)=a4*w(6)
      w(1)=463874074507.D0/190080.D0*z3 - 6272.D0/3.D0*w(1) - 80752.D0/
     & 45.D0*w(2) + 78954810295596095321843.D0/3163050541728000.D0 +
     & 50176*li4half
      w(1)=a5*w(1)
      w(1)=1.D0/9075.D0*w(1) + 1.D0/315.D0*w(6) + 1.D0/245.D0*w(5) +
     & w(3) + 1.D0/25.D0*w(4)

      MOMPSF = 1.D0/10125.D0*w(1)
*
      RETURN
      END
      REAL*8 FUNCTION MOMQGF(a1,a2,a3,a4,a5)
*     --------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 19, 2024
*     Fixed Mellin moments for AQG (2,4,6,8,10): test case
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 Pi,ln2,z3,li4half,as
      INTEGER a1,a2,a3,a4,a5
      REAL*8 w(7)
      COMMON/as/as
*
      Pi =  3.1415926535897932385D0
      ln2 = 0.69314718055994530942D0
      z3 = 1.2020569031595942854D0 
      li4half = 0.51747906167389938633D0
*
      w(1)=Pi**2
      w(2)=4 + 6929.D0/3240.D0*w(1)
      w(2)=w(2)*w(1)
      w(3)=ln2**2
      w(3)=w(3) - w(1)
      w(3)=w(3)*ln2
      w(4)=8*w(1) - 83.D0/9.D0*w(3)
      w(5)=1.D0/3.D0*ln2
      w(4)=w(4)*w(5)
      w(2)= - 664.D0/9.D0*li4half + w(4) + 762329.D0/648.D0*z3 - 
     & 150616519.D0/23328.D0 + w(2)
      w(2)=a1*w(2)
      w(4)=11*w(1)
      w(6)=1 + 636331.D0/453600.D0*w(1)
      w(6)=w(6)*w(4)
      w(7)=22*w(1) - 80147.D0/1260.D0*w(3)
      w(7)=w(7)*w(5)
      w(6)= - 160294.D0/315.D0*li4half + w(7) + 30965132157571.D0/
     & 3251404800.D0*z3 - 38691286402008212737.D0/209105971200000.D0 + 
     & w(6)
      w(6)=a3*w(6)
      w(2)=w(2) + 1.D0/7.D0*w(6)
      w(6)=1 + 2279.D0/1920.D0*w(1)
      w(4)=w(6)*w(4)
      w(6)=22.D0/3.D0*w(1) - 289.D0/16.D0*w(3)
      w(6)=ln2*w(6)
      w(4)=w(6) + 128152427233.D0/16588800.D0*z3 - 3405951881231689.D0/
     & 27993600000.D0 + w(4)
      w(4)=1.D0/3.D0*w(4) - 289.D0/2.D0*li4half
      w(4)=a2*w(4)
      w(6)=37 + 619499579.D0/10886400.D0*w(1)
      w(6)=w(6)*w(1)
      w(6)=2118626955610810711.D0/59929893273600.D0*z3 - 
     & 246504354200960177529497533.D0/312193542153830400000.D0 + w(6)
      w(7)=37*w(1) - 7045343.D0/60480.D0*w(3)
      w(7)=w(7)*w(5)
      w(6)= - 7045343.D0/7560.D0*li4half + 1.D0/2.D0*w(6) + w(7)
      w(6)=a4*w(6)
      w(7)=56 + 10897393.D0/118800.D0*w(1)
      w(7)=w(7)*w(1)
      w(1)=112*w(1) - 11201.D0/30.D0*w(3)
      w(1)=w(1)*w(5)
      w(1)= - 44804.D0/15.D0*li4half + w(1) + 53186528046305136613.D0/
     & 1054766121615360.D0*z3 - 241791674486509332299960601778223.D0/
     & 172385424539713850572800000.D0 + w(7)
      w(1)=a5*w(1)
      w(1)=1.D0/165.D0*w(1) + 1.D0/45.D0*w(6) + 1.D0/3.D0*w(2) + 1.D0/5.
     & D0*w(4)
*
      MOMQGF = 1.D0/375.D0*w(1)
*
      RETURN
      END
      REAL*8 FUNCTION MOMQQNS(a1,a2,a3,a4,a5)
*     ---------------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, March 19, 2024
*     Fixed Mellin moments for AqqQNS (2,4,6,8,10): test case
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL*8 Pi,ln2,z3,li4half,as
      INTEGER a1,a2,a3,a4,a5
      REAL*8 w(6)
      COMMON/as/as
*
      Pi =  3.1415926535897932385D0
      ln2 = 0.69314718055994530942D0
      z3 = 1.2020569031595942854D0 
      li4half = 0.51747906167389938633D0
*
      w(1)=Pi**2
      w(2)=ln2**2
      w(1)=w(1) - w(2)
      w(1)=w(1)*w(2)
      w(2)=Pi**4
      w(3)= - 73027.D0/90.D0*w(2) + 1944399704556653.D0/740880000.D0 + 
     & 22688*li4half
      w(3)= - 2836.D0/9.D0*w(1) + 1.D0/3.D0*w(3) - 17446769.D0/140.D0*
     & z3
      w(3)=a3*w(3)
      w(4)= - 16.D0/3.D0*w(1) - 6721.D0/3.D0*z3 - 206.D0/45.D0*w(2) + 
     & 2661823.D0/243.D0 + 128*li4half
      w(4)=a1*w(4)
      w(5)= - 314.D0/3.D0*w(1) - 5133173.D0/120.D0*z3 - 16171.D0/180.D0
     & *w(2) + 5292070523963.D0/19440000.D0 + 2512*li4half
      w(5)=a2*w(5)
      w(6)= - 9883.D0/3.D0*w(1) - 6342799481.D0/5040.D0*z3 - 1017949.D0/
     & 360.D0*w(2) + 162741536890063389677.D0/17283248640000.D0 + 79064
     & *li4half
      w(6)=a4*w(6)
      w(3)=1.D0/315.D0*w(6) + 1.D0/15.D0*w(5) + 1.D0/35.D0*w(3) + 2.D0/
     & 3.D0*w(4)
      w(1)= - 4822.D0/3.D0*w(1) - 410587711097.D0/693000.D0*z3 - 248333.
     & D0/180.D0*w(2) + 131748303499583029785257.D0/28241522694000000.D0
     &  + 38576*li4half
      w(1)=a5*w(1)
      w(1)=1.D0/5.D0*w(3) + 1.D0/693.D0*w(1)
*
      MOMQQNS = 1.0D0 + 1.D0/675.D0*w(1)
*
      RETURN
      END
      REAL*8 FUNCTION AQG(z,nf,as,LL)
*     -------------------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of aQg3: not yet available. to be supplemented
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 z,nf,as,LL,aQg3
*
*      AQG = aQg3(z)*as**3
       AQG = as**3*0.0D0
*
      RETURN
      END      
      REAL*8 FUNCTION aQg3old(z)
*     -----------------------
*     ------------------------------------------------------------------------
*     Code: J. Bluemlein, February 20, 2024
*     regular part of aQg3: not yet available. to be supplemented
*     without as-factor  
*     ------------------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 z,nf,as,LL
*
      aqg3old=0.0D0
*
      RETURN
      END      
*     
      !include 'hplog5.f'
!include 'daind.f'

      
      
      end module
