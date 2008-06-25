*************************************************************
*
*     File: umatrix.f
*     
*     Computes the matrices used for the singlet evolution. 
*     Equation numbers refer to the "Notes on PT evolution"
*
*     Input variables: IPT : perturbative order
*                      ZN   : Mellin variable
*                      NF  : number of flavours
*
*     Functions:
*
*     Output: LP,LM     : LO matrix eigenvalues
*             EP,EM     : LO matrix eigenvectors
*             U(1..20)  : NLO evolution matrices appearing into the 
*                         iterative solution
*             UNS0      : LO non singlet evolution coefficient
*             UNS1      : NLO non singlet evolution coefficient
*
*     MU: 10.03.07
*
*************************************************************
*
      SUBROUTINE UMATRIX(ZN,NF,LP,LM,EP,EM,U,UNS0,UNS1)
*
      IMPLICIT none
*
      include "../../commons/beta.h"
      include "../../commons/evflags.h"
*
      INTEGER NF
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX UNS0,UNS1(3)
      DOUBLE COMPLEX U(20,2,2)
      DOUBLE COMPLEX LP,LM
      DOUBLE COMPLEX EP(2,2),EM(2,2)
*
      INTEGER I,J,K,JJ,KK,II
      DOUBLE COMPLEX QQ0,QG0,GQ0,GG0
      DOUBLE COMPLEX SQ,LDIFF
      DOUBLE COMPLEX DP(2,2),DM(2,2)
      DOUBLE COMPLEX G0(2,2),G1(2,2)
      DOUBLE COMPLEX G0NS(3),G1NS(3)
      DOUBLE COMPLEX R(20,2,2),RT(20,2,2)
      DOUBLE COMPLEX R0(2,2),R1(2,2)
      DOUBLE COMPLEX R1PP(2,2),R1PM(2,2),R1MP(2,2),R1MM(2,2)
      DOUBLE COMPLEX U1PP(2,2),U1PM(2,2),U1MP(2,2),U1MM(2,2)
*
      CALL ANDIM_LO(ZN,NF,G0NS,G0)
*
      UNS0 = G0NS(1) / BETA0(NF)
*
      DO K=1,20
         DO I=1,2
            DO J=1,2
               R(K,I,J)=(0.D0,0.D0)
               RT(K,I,J)=(0.D0,0.D0)
               U(K,I,J)=(0.D0,0.D0)
            ENDDO
         ENDDO
      ENDDO

*
      DO I=1,2
         DO J=1,2
            R0(I,J) = G0(I,J) / BETA0(NF) ! eq. (18)
         END DO
      END DO
*
      QQ0 = R0(1,1)
      QG0 = R0(1,2)
      GQ0 = R0(2,1)
      GG0 = R0(2,2)
*
      SQ = SQRT((QQ0-GG0)**2d0 + 4.d0*QG0*GQ0)   ! eq.(20)
      LP = .5d0*(QQ0 + GG0 + SQ)
      LM = .5d0*(QQ0 + GG0 - SQ)
*
      EM(1,1) = (QQ0 - LP) / (LM - LP)         ! eq.(21)
      EM(1,2)= QG0 / (LM - LP)
      EM(2,1) = GQ0 / (LM - LP)
      EM(2,2) = (GG0 - LP) / (LM - LP)
*
      EP(1,1) = (QQ0 - LM) / (LP - LM)
      EP(1,2) = QG0 / (LP - LM)
      EP(2,1) = GQ0 / (LP - LM)
      EP(2,2) = (GG0 - LM) / (LP - LM)
*
*     NLO coefficients
*

      IF (IPT.GE.1) THEN        
         CALL ANDIM_NLO(ZN,NF,G1NS,G1)
*     
*     NON SINGLET
*
         DO I=1,3        
            UNS1(I)= -(G1NS(I)/BETA0(NF)) + B1(NF)*UNS0
         ENDDO
*     
*     SINGLET
* 
         DO I=1,2
            DO J=1,2
               R(1,I,J)  = G1(I,J)/BETA0(NF) - B1(NF)*R0(I,J) ! eq. (18)
               RT(1,I,J)  = G1(I,J)/BETA0(NF) - B1(NF)*R0(I,J) ! R1_T=R1
               R1(I,J) = R(1,I,J)
           ENDDO
         ENDDO
*     
*     Computation of U1 according to eq.(25)
*     
         CALL MMULT_C(EP,2,2,R1,2,2,DP)
         CALL MMULT_C(EM,2,2,R1,2,2,DM)
*     
         CALL MMULT_C(DP,2,2,EP,2,2,R1PP)
         CALL MMULT_C(DP,2,2,EM,2,2,R1PM)
         CALL MMULT_C(DM,2,2,EP,2,2,R1MP)
         CALL MMULT_C(DM,2,2,EM,2,2,R1MM)
*
         DO I = 1,2
            DO J = 1,2
               U1PP(I,J) = R1PP(I,J)
               U1PM(I,J) = R1PM(I,J)/(LM - LP - 1d0)
               U1MP(I,J) = R1MP(I,J)/(LP - LM - 1d0)
               U1MM(I,J) = R1MM(I,J)
               U(1,I,J) = -(U1MM(I,J)+U1PP(I,J)) + U1PM(I,J) + U1MP(I,J)
            ENDDO
         ENDDO

         IF(IMODEV.EQ.0)THEN
            DO K=2,20
               DO I=1,2
                  DO J=1,2
                     U(K,I,J)=0.D0
                  ENDDO
               ENDDO
            ENDDO

         ELSEIF(IMODEV.EQ.1)THEN
*
*     Computation of U_k (k>=2) at NLO according to eqs (23,26,27)
*
            LDIFF = LM - LP
            
            DO K=2,20
               DO I=1,2
                  DO J=1,2
                     R(K,I,J) = - B1(NF) * R(K-1,I,J) ! eq (27)
                  ENDDO
               ENDDO
            ENDDO

            DO K=2,20
               DO I=1,2
                  DO J=1,2
                     RT(K,I,J) = R(K,I,J)
                     DO JJ=1,K-1
                        DO KK=1,2                      
                           RT(K,I,J) = RT(K,I,J) +
     #                          ( R(JJ,I,KK) * U(K-JJ,KK,J) ) ! eq (26)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
               
               DO I=1,2
                  DO J=1,2
                     DO II=1,2
                        DO JJ=1,2                       
                           U(K,I,J) = U(K,I,J) ! eq(23)
     #                          -(EM(I,II)*RT(K,II,JJ)*EM(JJ,J)/DBLE(K))
     #                          -(EP(I,II)*RT(K,II,JJ)*EP(JJ,J)/DBLE(K))
     #                          -(EM(I,II)*RT(K,II,JJ)* EP(JJ,J)/
     #                          ( DBLE(K) + LDIFF ))
     #                          -(EP(I,II)*RT(K,II,JJ) 
     #                          * EM(JJ,J) / ( DBLE(K) - LDIFF ))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END
