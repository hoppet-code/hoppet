**********************************************************************
*
*     File: evolfactn.f
*
*     Returns the evolution factors in N space, according to eq. (20)
*     of the "Notes on Perturbative Evolution".
*
*     Comments: as = alphas/(4*pi)
*
**********************************************************************
*
*



      SUBROUTINE ZFUNC(ZN,X,Q20,Q2,ZFUNCNS,ZFUNCNS24,ZFUNCSG)
*
*     Returns the evolution factor Gamma(N), from Q20 to Q2 
* 
*      - Fixed Flavour Number Scheme (IVFN=0)
*      - Zero Mass-Variable Flavour Number Scheme (IVFN=1)
*
      IMPLICIT none
*
      include "../../commons/alphas.h"
      include "../../commons/evflags.h"
      include "../../commons/expptkin.h"
*
      INTEGER I,J
      INTEGER NFI,NFF,STEP
*
      DOUBLE PRECISION Q20,Q2
      DOUBLE PRECISION X
*
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX EFNNS(4),EFNSG(2,2)
      DOUBLE COMPLEX EFNNS0(4),EFNSG0(2,2),EFNSG1(2,2)
      DOUBLE COMPLEX ZFUNCNS(4),ZFUNCNS24(2),ZFUNCSG(2,2)
*
      DOUBLE COMPLEX EFNNSBQ(4),EFNNSCB(4),EFNNSBT(4),EFNNSTQ(4)
      DOUBLE COMPLEX EFNSGBQ(2,2),EFNSGCB(2,2),EFNSGBT(2,2),EFNSGTQ(2,2)
      DOUBLE COMPLEX CNS(4),CSG(2)
*
*     Allow backward evolution      
*
      IF(Q20.LE.Q2) STEP =  1
      IF(Q20.GT.Q2) STEP = -1
*
*     Determine # of active flavours at the initial (NFI) and 
*     final (NFF) scale.
*
      IF(Q20.GE.Q2TH(6)) THEN
         NFI=6
      ELSEIF(Q20.GE.Q2TH(5)) THEN
         NFI=5
      ELSEIF(Q20.GE.Q2TH(4)) THEN
         NFI=4
      ELSE
         NFI=3
      ENDIF
*
      IF(Q2.GT.Q2TH(6)) THEN
         NFF=6
      ELSEIF(Q2.GT.Q2TH(5)) THEN
         NFF=5
      ELSEIF(Q2.GT.Q2TH(4)) THEN
         NFF=4
      ELSE
         NFF=3
      ENDIF

*
*     Coefficient Functions 
*
      IF(IVFN.EQ.0) NFF=4
      CALL WCOEFF(ZN,NFF,X,Q2,Y,CNS,CSG)
*     
*     Evolution Kernel
*      
      DO I=1,4
         EFNNS0(I) = (0d0,0d0)
      ENDDO
      DO I=1,2
         DO J=1,2
            EFNSG0(I,J) = (0.D0,0.D0)
            EFNSG1(I,J) = (0.D0,0.D0)
            EFNSG(I,J) = (0.D0,0.D0)
         END DO
      END DO
*
*     Fixed Flavour Number Scheme (NF=4)
*
      IF (IVFN.EQ.0) THEN
         NFI = 4
         CALL EVOLFACTN(ZN,Q20,Q2,NFI,EFNNS,EFNSG1)
*     
*     (Zero Mass-) Variable Flavour Number Scheme
*     
      ELSEIF(IVFN.EQ.1) THEN
         IF(NFF.EQ.NFI) THEN
            CALL EVOLFACTN(ZN,Q20,Q2,NFI,EFNNS,EFNSG1)
         ELSE
            CALL EVOLFACTN(ZN,Q20,Q2TH(NFI+STEP),NFI,EFNNS0,EFNSG0)
            DO I=1,4
               EFNNS(I) = EFNNS0(I)
            ENDDO
            CALL MEQUAL_C(EFNSG,EFNSG0,2,2)
 10         NFI = NFI+STEP
            IF(NFI.NE.NFF) THEN 
               CALL EVOLFACTN(ZN,Q2TH(NFI),Q2TH(NFI+STEP),NFI,
     #              EFNNS0,EFNSG0)
            ELSE
               CALL EVOLFACTN(ZN,Q2TH(NFI),Q2,NFI,EFNNS0,EFNSG0)
            ENDIF
            DO I=1,4
               EFNNS(I) = EFNNS(I)*EFNNS0(I)
            ENDDO
            CALL MMULT_C(EFNSG0,2,2,EFNSG,2,2,EFNSG1)
            
            IF(NFI.NE.NFF) THEN
               CALL MEQUAL_C(EFNSG,EFNSG1,2,2)
               GOTO 10         
            ENDIF
         ENDIF
*
*     Evolution of T_24     
*     
         IF(Q2.LE.Q2TH(5))THEN
            ZFUNCNS24(1) = EFNSG1(1,1)
            ZFUNCNS24(2) = EFNSG1(1,2)
         ELSEIF(Q2.LE.Q2TH(6)) THEN
            CALL EVOLFACTN(ZN,Q20,Q2TH(5),4,EFNNSCB,EFNSGCB)
            CALL EVOLFACTN(ZN,Q2TH(5),Q2 ,5,EFNNSBQ,EFNSGBQ)
            ZFUNCNS24(1) = EFNNSBQ(1)*EFNSGCB(1,1)
            ZFUNCNS24(2) = EFNNSBQ(1)*EFNSGCB(1,2)         
         ELSE
            CALL EVOLFACTN(ZN,Q20,Q2TH(5),4,EFNNSCB,EFNSGCB)
            CALL EVOLFACTN(ZN,Q2TH(5),Q2TH(6),5,EFNNSBT,EFNSGBT)
            CALL EVOLFACTN(ZN,Q2TH(6),Q2,6,EFNNSTQ,EFNSGTQ)
            ZFUNCNS24(1) = EFNNSTQ(1)*EFNNSBT(1)*EFNSGCB(1,1)
            ZFUNCNS24(2) = EFNNSTQ(1)*EFNNSBT(1)*EFNSGCB(1,2)         
         ENDIF
      ENDIF
*      
*     Output
*
      DO I=1,4
         ZFUNCNS(I) = CNS(I) * EFNNS(I)
      ENDDO

      DO I=1,2
c         ZFUNCNS24(I)= CNS(1) * ZFUNCNS24(I) 
         ZFUNCNS24(I)= CSG(1) * ZFUNCNS24(I)
      ENDDO
      
      DO I=1,2
         DO J=1,2
            ZFUNCSG(I,J)= CSG(I) * EFNSG1(I,J)
         ENDDO
      ENDDO

      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE EVOLFACTN(ZN,Q2I,Q2F,NF,EFNNS,EFNSG)
*
*     Returns the evolution factors, in N space, from Q20 to Q2 with NF 
*     active flavours, for Non Singlet and Singlet-Gluon.   
*     The PT order is determined through the EVFLAGS common block.
*
      IMPLICIT none
*
      include "../../commons/evscale.h"
      include "../../commons/evflags.h"
      include "../../commons/beta.h"
      include "../../commons/consts.h"
      include "../../commons/alphas.h"
*
      INTEGER NF
      INTEGER I,J,K
      DOUBLE PRECISION Q2I,Q2F,ASI,ASF,T
      DOUBLE PRECISION TMP
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX EFNNS(4),EFNSG(2,2),EFNSGTMP(2,2)
      DOUBLE COMPLEX EXPNS,EXPM,EXPP
      DOUBLE COMPLEX UNS0,UNS1(3),LP,LM
      DOUBLE COMPLEX U(20,2,2)
      double complex EM(2,2),EP(2,2)
!
      DOUBLE COMPLEX U1(2,2)
      DOUBLE COMPLEX L(2,2),LU1(2,2),U1L(2,2)
      DOUBLE COMPLEX USUM(2,2),USUMTMP(2,2),UINV(2,2),USUML(2,2)
      DOUBLE COMPLEX DETINV
*
*     Compute alphas 
*
      IF(Q2I.EQ.Q20) THEN
         ASI = AS0
      ELSEIF(Q2I.EQ.Q2TH(4)) THEN
         ASI = ASC
      ELSEIF(Q2I.EQ.Q2TH(5)) THEN
         ASI = ASB
      ELSEIF(Q2I.EQ.Q2TH(6)) THEN
         ASI = AST
      ENDIF
      
      IF(Q2F.EQ.Q20) THEN
         ASF = AS0
      ELSEIF(Q2F.EQ.Q2TH(4)) THEN
         ASF = ASC
      ELSEIF(Q2F.EQ.Q2TH(5)) THEN
         ASF = ASB
      ELSEIF(Q2F.EQ.Q2TH(6)) THEN
         ASF = AST
      ELSE
         ASF = ASQ
      ENDIF
*     
*     Evaluation of evolution factors
*
*
      T = DLOG(ASF/ASI)

      DO I=1,4
         EFNNS(I) = (0d0,0d0)
      ENDDO
      
      DETINV=0.D0
      DO I=1,2
         DO J=1,2
            UINV(I,J)=(0d0,0d0)
            EFNSG(I,J) = (0d0,0d0)
            EFNSGTMP(I,J) = (0d0,0d0)
         ENDDO
      ENDDO
*
*     U matrices
*
      CALL UMATRIX(ZN,NF,LP,LM,EP,EM,U,UNS0,UNS1)
*
*     LO evolution factor 
* 
*     Non singlet
*
      EXPNS  = EXP ( - UNS0 * T )
*
*     Singlet
*
      EXPM   = EXP ( - LM * T ) 
      EXPP   = EXP ( - LP * T )

      DO I=1,2 
        DO J=1,2 
           L(I,J) = EXPM * EM(I,J) + EXPP * EP(I,J) ! eq.(22) 
        ENDDO 
      ENDDO 

      IF(IPT.EQ.0)THEN
*
*     LO solution
*
         DO I= 1,4
            EFNNS(I) = EXPNS  
         ENDDO
         DO I=1,2 
            DO J=1,2 
               EFNSG(I,J) = L(I,J)                       
            ENDDO 
         ENDDO 

      ELSEIF(IPT.EQ.1)THEN
*     
*     NLO solution
*     
         IF(IMODEV.EQ.0)THEN
*
*     Truncated solution IMODEV=0
*     
            DO I=1,3            
               EFNNS(I)= EXPNS * ((1D0,0D0) + UNS1(I) * (ASF-ASI))
            ENDDO
            EFNNS(4)=EFNNS(3)

            DO I=1,2
               DO J=1,2
                  U1(I,J)=U(1,I,J)
               ENDDO
            ENDDO

            CALL MMULT_C(U1,2,2,L,2,2,U1L)
            CALL MMULT_C(L,2,2,U1,2,2,LU1)
            DO I=1,2 
               DO J=1,2 
                  EFNSG(I,J) = L(I,J) + ASF * U1L(I,J) - ASI * LU1(I,J)
               ENDDO 
            ENDDO 

         ELSEIF(IMODEV.EQ.1)THEN
*     
*     Iterated solution IMODEV=1
*     
            TMP = 0D0
            TMP = DLOG(( 1D0 + B1(NF) * ASF )/(1D0 + B1(NF)*ASI))/B1(NF)
            DO I=1,3
               EFNNS(I) = EXPNS * EXP( TMP * UNS1(I))        
            ENDDO
            EFNNS(4)=EFNNS(3)
*     
            USUM(1,1) = (1d0,0d0)
            USUM(1,2) = (0d0,0d0)
            USUM(2,1) = (0d0,0d0)
            USUM(2,2) = (1d0,0d0)
            USUMTMP(1,1) = (1d0,0d0)
            USUMTMP(1,2) = (0d0,0d0)
            USUMTMP(2,1) = (0d0,0d0)
            USUMTMP(2,2) = (1d0,0d0)
*
            DO K=1,20
               DO I=1,2
                  DO J=1,2
                     USUM(I,J) = USUM(I,J)+ 
     #                    (ASF**(DBLE(K)) * U(K,I,J))
                     USUMTMP(I,J) = USUMTMP(I,J) + 
     #                    (ASI**(DBLE(K)) * U(K,I,J))
                   ENDDO
               ENDDO
            ENDDO
*


            DETINV = 1.D0 / ( USUMTMP(1,1) * USUMTMP(2,2)-
     #           USUMTMP(1,2) * USUMTMP(2,1) )
            UINV(1,1) = DETINV * USUMTMP(2,2)
            UINV(1,2) = - DETINV * USUMTMP(1,2)
            UINV(2,1) = - DETINV * USUMTMP(2,1)
            UINV(2,2) = DETINV * USUMTMP(1,1)
           
            CALL MMULT_C(USUM,2,2,L,2,2,USUML)
            CALL MMULT_C(USUML,2,2,UINV,2,2,EFNSGTMP)   
            
            DO I=1,2
               DO J=1,2
                  EFNSG(I,J)= EFNSGTMP(I,J)
               ENDDO
            ENDDO
         ENDIF
         
      ENDIF

      RETURN
      END
      
*========================================================================

