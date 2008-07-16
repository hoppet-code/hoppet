**********************************************************************
*
*     File: wcoeff.f
*
*     Returns the Coefficient Functions in Mellin space for the
*     observable under consideration.
*     
*     The perturbative order the observable
*
*     Comments: 
*         - as = alphas/(4*pi)
*         - Correspondence table OBSERV --> observable:
*                0: PDF
*                1: F2P
*                2: F2D
*                3: F2P/F2D
*                4: RED SIGMA^NC/X e^+p
*                5: RED SIGMA^CC/X e^+p
*                6: RED SIGMA^NC/X e^-p
*                7: RED SIGMA^CC/X e^-p 
*                8: FL
*                9: SIGMA neutrino
*               10: SIGMA antineutrino
*
**********************************************************************
*
      SUBROUTINE WCOEFF(ZN,NF,X,Q2,Y,CNS,CSG)
*
      IMPLICIT none
*
      include "../../commons/alphas.h"
      include "../../commons/evflags.h"
      include "../../commons/observ.h"
*
      INTEGER NF,J
      DOUBLE PRECISION X, Y
*
      DOUBLE COMPLEX ZN
      DOUBLE COMPLEX C2QTMC(0:2),C2GTMC(0:2)
      DOUBLE COMPLEX C3QTMC(0:2),C3GTMC(0:2)
      DOUBLE COMPLEX CLQTMC(0:2),CLGTMC(0:2)
*
      DOUBLE COMPLEX CNS(4),CSG(2)
      DOUBLE COMPLEX C2Q(0:2),C2G(0:2),C3Q(0:2)
      DOUBLE COMPLEX C3G(0:2),CLQ(0:2),CLG(0:2)
      DOUBLE COMPLEX TMC(3)
      DOUBLE PRECISION YM,YP,MN
      DOUBLE PRECISION Q2
*
      MN = .9389d0 
      YM= 1D0 - (1D0 - Y)**2D0
      YP= 1D0 + (1D0 - Y)**2D0 

      CALL C2N(ZN,NF,C2Q,C2G)
      CALL C3N(ZN,NF,C3Q,C3G)
      CALL CLN(ZN,NF,CLQ,CLG)
*     
      IF(OBSERV.EQ.0) THEN
         CNS(1) = (1.d0,0.d0)
         CNS(2) = (1.d0,0.d0)
         CNS(3) = (1.d0,0.d0)
         CNS(4) = (1.d0,0.d0)
         CSG(1) = (1.d0,0.d0)
         CSG(2) = (1.d0,0.d0)
*
      ELSEIF(OBSERV.EQ.1) THEN
         CNS(1) = C2Q(0)
         CNS(2) = C2Q(0)
         CNS(3) = C2Q(0)
         CNS(4) = C2Q(0)
         CSG(1) = C2Q(0)
         CSG(2) = C2G(0)
         IF(IPT.GE.1) THEN
            CNS(1) = CNS(1) + ASQ * C2Q(1)
            CNS(2) = CNS(2) + ASQ * C2Q(1)
            CNS(3) = CNS(3) + ASQ * C2Q(1)
            CNS(4) = CNS(4) + ASQ * C2Q(1)
            CSG(1) = CSG(1) + ASQ * C2Q(1)
            CSG(2) = CSG(2) + ASQ * C2G(1)
         ENDIF   
         IF(ITMC.EQ.1) THEN
            CALL COMPUTE_TMC(ZN,X,Q2,TMC)
            DO J=1,4
               CNS(J) = TMC(1) * CNS(J)
            ENDDO
            DO J=1,2
               CSG(J) = TMC(1) * CSG(J)
            ENDDO
         ENDIF
*     
      ELSEIF(OBSERV.EQ.2) THEN
         CNS(1) = C2Q(0)
         CNS(2) = C2Q(0)
         CNS(3) = C2Q(0)
         CNS(4) = C2Q(0)
         CSG(1) = C2Q(0)
         CSG(2) = C2G(0)
         IF(IPT.GE.1) THEN

            CNS(1) = CNS(1) + ASQ * C2Q(1)
            CNS(2) = CNS(2) + ASQ * C2Q(1)
            CNS(3) = CNS(3) + ASQ * C2Q(1)
            CNS(4) = CNS(4) + ASQ * C2Q(1)
            CSG(1) = CSG(1) + ASQ * C2Q(1)
            CSG(2) = CSG(2) + ASQ * C2G(1)
         ENDIF   
         IF(ITMC.EQ.1) THEN
            CALL COMPUTE_TMC(ZN,X,Q2,TMC)
            DO J=1,4
               CNS(J) = TMC(1) * CNS(J)
            ENDDO
            DO J=1,2
               CSG(J) = TMC(1) * CSG(J)
            ENDDO
         ENDIF
*     
      ELSEIF(OBSERV.EQ.3) THEN

         CNS(1) = C2Q(0)
         CNS(2) = C2Q(0)
         CNS(3) = C2Q(0)
         CNS(4) = C2Q(0)
         CSG(1) = C2Q(0)
         CSG(2) = C2G(0)
         IF(IPT.GE.1) THEN

            CNS(1) = CNS(1) + ASQ * C2Q(1)
            CNS(2) = CNS(2) + ASQ * C2Q(1)
            CNS(3) = CNS(3) + ASQ * C2Q(1)
            CNS(4) = CNS(4) + ASQ * C2Q(1)
            CSG(1) = CSG(1) + ASQ * C2Q(1)
            CSG(2) = CSG(2) + ASQ * C2G(1)
         ENDIF   
         IF(ITMC.EQ.1) THEN
            CALL COMPUTE_TMC(ZN,X,Q2,TMC)
            DO J=1,4
               CNS(J) = TMC(1) * CNS(J)
            ENDDO
            DO J=1,2
               CSG(J) = TMC(1) * CSG(J)
            ENDDO
         ENDIF
         
*     Same coefficient functions for e^+p and e^-p reduced xsec (NC)
      ELSEIF(OBSERV.EQ.4.or.OBSERV.EQ.6) THEN     
         IF(ITMC.EQ.1) THEN
            CALL COMPUTE_TMC(ZN,X,Q2,TMC)

            DO J=0,1
               C2QTMC(J) = C2Q(J) * (TMC(1)-(Y**2D0/YP)*TMC(2))
               C2GTMC(J) = C2G(J) * (TMC(1)-(Y**2D0/YP)*TMC(2))
               C3QTMC(J) = C3Q(J) * TMC(3)
               C3GTMC(J) = C3G(J) * TMC(3)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
            ENDDO
         ELSEIF (ITMC.EQ.0) THEN
            DO J=0,1
               C2QTMC(J) = C2Q(J)
               C2GTMC(J) = C2G(J)
               C3QTMC(J) = C3Q(J)
               C3GTMC(J) = C3G(J)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
            ENDDO
         ENDIF

         CNS(1) = C2QTMC(0) - (Y**2D0/YP) * CLQTMC(0)
         CNS(2) = C2QTMC(0) - (Y**2D0/YP) * CLQTMC(0)
         CNS(3) = (YM/YP) * C3QTMC(0)
         CNS(4) = (YM/YP) * C3QTMC(0)
         CSG(1) = C2QTMC(0) - (Y**2D0/YP) * CLQTMC(0)
         CSG(2) = C2GTMC(0) - (Y**2D0/YP) * CLGTMC(0)
         IF(IPT.GE.1) THEN

            CNS(1) = CNS(1) + ASQ * (C2QTMC(1)-(Y**2D0/YP)*CLQTMC(1))
            CNS(2) = CNS(2) + ASQ * (C2QTMC(1)-(Y**2D0/YP)*CLQTMC(1))
            CNS(3) = CNS(3) + ASQ * (YM/YP)*C3QTMC(1)
            CNS(4) = CNS(4) + ASQ * (YM/YP)*C3QTMC(1)
            CSG(1) = CSG(1) + ASQ * (C2QTMC(1)-(Y**2D0/YP)*CLQTMC(1))
            CSG(2) = CSG(2) + ASQ * (C2GTMC(1)-(Y**2D0/YP)*CLGTMC(1))
         ENDIF   
*
      ELSEIF(OBSERV.EQ.5) THEN

         IF(ITMC.EQ.1) THEN
            CALL COMPUTE_TMC(ZN,X,Q2,TMC)
            DO J=0,1
               C2QTMC(J) = C2Q(J) * (TMC(1)-(Y**2D0/YP)*TMC(2))
               C2GTMC(J) = C2G(J) * (TMC(1)-(Y**2D0/YP)*TMC(2))
               C3QTMC(J) = C3Q(J) * TMC(3)
               C3GTMC(J) = C3G(J) * TMC(3)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
            ENDDO
         ELSEIF (ITMC.EQ.0) THEN
            DO J=0,1
               C2QTMC(J) = C2Q(J)
               C2GTMC(J) = C2G(J)
               C3QTMC(J) = C3Q(J)
               C3GTMC(J) = C3G(J)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
            ENDDO
         ENDIF
         CNS(1) = YM * C3QTMC(0)
         CNS(2) = YM * C3QTMC(0)
         CNS(3) = - YP * C2QTMC(0) + Y**2D0 * CLQTMC(0) 
     #        - YM * C3QTMC(0)
         CNS(4) =  YP * C2QTMC(0) - Y**2D0 * CLQTMC(0) 
     #        - YM * C3QTMC(0)
         CSG(1) = YP * C2QTMC(0) - Y**2D0 *  CLQTMC(0)
         CSG(2) = YP * C2GTMC(0) - Y**2D0 *  CLGTMC(0)
         IF(IPT.GE.1) THEN
            CNS(1) = CNS(1) + ASQ * (YM * C3QTMC(1))
            CNS(2) = CNS(2) + ASQ * (YM * C3QTMC(1))
            CNS(3) = CNS(3) + ASQ * (-YP * C2QTMC(1) +
     #           Y**2D0*CLQTMC(1) - YM * C3QTMC(1))
            CNS(4) = CNS(4) + ASQ * (YP * C2QTMC(1) -Y**2D0*CLQTMC(1) 
     #           - YM * C3QTMC(1))
            CSG(1) = CSG(1) + ASQ *(YP * C2QTMC(1) -Y**2D0*CLQTMC(1))
            CSG(2) = CSG(2) + ASQ *(YP * C2GTMC(1) -Y**2D0*CLGTMC(1)) 
         ENDIF
*
      ELSEIF(OBSERV.EQ.7) THEN

         IF(ITMC.EQ.1) THEN
            CALL COMPUTE_TMC(ZN,X,Q2,TMC)
            DO J=0,1
               C2QTMC(J) = C2Q(J) *  (TMC(1)-(Y**2D0/YP)*TMC(2))
               C2GTMC(J) = C2G(J) *  (TMC(1)-(Y**2D0/YP)*TMC(2))
               C3QTMC(J) = C3Q(J) * TMC(3)
               C3GTMC(J) = C3G(J) * TMC(3)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
            ENDDO
         ELSEIF (ITMC.EQ.0) THEN
            DO J=0,1
               C2QTMC(J) = C2Q(J)
               C2GTMC(J) = C2G(J)
               C3QTMC(J) = C3Q(J)
               C3GTMC(J) = C3G(J)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
            ENDDO
         ENDIF
         CNS(1) = YM * C3QTMC(0)
         CNS(2) = YM * C3QTMC(0)
c         CNS(3) = YP * C2QTMC(0) - Y**2D0 * CLQTMC(0) - YM * C3QTMC(0)
c         CNS(4) = - YP * C2QTMC(0) + Y**2D0 * CLQTMC(0) - YM * C3QTMC(0)
         CNS(3) = YP * C2QTMC(0) - Y**2D0 * CLQTMC(0) + YM * C3QTMC(0)
         CNS(4) = - YP * C2QTMC(0) + Y**2D0 * CLQTMC(0) + YM * C3QTMC(0)
         CSG(1) = YP * C2QTMC(0) - Y**2D0 *  CLQTMC(0)
         CSG(2) = YP * C2GTMC(0) - Y**2D0 *  CLGTMC(0)  
         IF(IPT.GE.1) THEN
            CNS(1) = CNS(1) + ASQ * (YM * C3QTMC(1))
            CNS(2) = CNS(2) + ASQ * (YM * C3QTMC(1))
            CNS(3) = CNS(3) + ASQ * (YP * C2QTMC(1) -Y**2D0*CLQTMC(1) 
     #           - YM * C3QTMC(1))
            CNS(4) = CNS(4) + ASQ *(-YP * C2QTMC(1) +Y**2D0*CLQTMC(1) 
     #           - YM*C3QTMC(1))
            CSG(1) = CSG(1) + ASQ * (YP * C2QTMC(1)-Y**2D0*CLQTMC(1))
            CSG(2) = CSG(2) + ASQ * (YP * C2GTMC(1)-Y**2D0*CLGTMC(1))
         ENDIF

*     Longitudinal structure function
      ELSEIF(OBSERV.EQ.8) THEN     
        
*     Target mass corrections not available for FL
         do J=0,1
            CLQTMC(J) = CLQ(J) 
            CLGTMC(J) = CLG(J) 
         enddo

         CNS(1) = CLQTMC(0)
         CNS(2) = CLQTMC(0)
         CNS(3) = CLQTMC(0)
         CNS(4) = CLQTMC(0)
         CSG(1) = CLQTMC(0)
         CSG(2) = CLGTMC(0)

         IF(IPT.GE.1) THEN
            
            CNS(1) = CNS(1) + ASQ * CLQTMC(1)
            CNS(2) = CNS(2) + ASQ * CLQTMC(1)
            CNS(3) = CNS(3) + ASQ * CLQTMC(1)
            CNS(4) = CNS(4) + ASQ * CLQTMC(1)
            CSG(1) = CSG(1) + ASQ * CLQTMC(1)
            CSG(2) = CSG(2) + ASQ * CLGTMC(1)
            
         ENDIF 
  
**     sigma neutrino 

      ELSEIF(OBSERV.EQ.9) THEN   
         IF(ITMC.EQ.1) THEN
            CALL COMPUTE_TMC(ZN,X,Q2,TMC)
            DO J=0,1
               C2QTMC(J) = C2Q(J)* (TMC(1) 
     #              -(Y**2.D0/(YP -(2.d0*(MN*X*Y)**2D0)/Q2))*TMC(2))
               C2GTMC(J) = C2G(J)* (TMC(1) 
     #              -(Y**2.D0/(YP -(2.d0*(MN*X*Y)**2D0)/Q2))*TMC(2))
               C3QTMC(J) = C3Q(J) * TMC(3)
               C3GTMC(J) = C3G(J) * TMC(3)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
             ENDDO
         ELSEIF (ITMC.EQ.0) THEN
            DO J=0,1
               C2QTMC(J) = C2Q(J)
               C2GTMC(J) = C2G(J)
               C3QTMC(J) = C3Q(J)
               C3GTMC(J) = C3G(J)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
            ENDDO
         ENDIF
         CNS(1) = YM * C3QTMC(0)
         CNS(2) = YM * C3QTMC(0)
         CNS(3) = - (YP -(2.d0*(MN*X*Y)**2D0)/Q2 ) * C2QTMC(0)
     #            + (Y**2D0) * CLQTMC(0) 
     #            + YM * C3QTMC(0)
         CNS(4) = + (YP -(2.D0*(MN*X*Y)**2D0)/Q2 ) * C2QTMC(0)
     #            - (Y**2D0) * CLQTMC(0) 
     #            + YM * C3QTMC(0)
         CSG(1) = (YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2QTMC(0) 
     #        - (Y**2D0) * CLQTMC(0)
         CSG(2) = (YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2GTMC(0) 
     #        - (Y**2D0) * CLGTMC(0)
         IF(IPT.GE.1) THEN
            CNS(1) = CNS(1) + ASQ * (YM * C3QTMC(1))
            CNS(2) = CNS(2) + ASQ * (YM * C3QTMC(1))
            CNS(3) = CNS(3) + ASQ * 
     #           (- (YP -(2.d0*(MN*X*Y)**2D0)/Q2 ) * C2QTMC(1)
     #           + (Y**2D0) * CLQTMC(1) 
     #           + YM * C3QTMC(1))
            CNS(4) = CNS(4) + ASQ *
     #          (+ (YP -(2.D0*(MN*X*Y)**2D0)/Q2 ) * C2QTMC(1)
     #           - (Y**2D0) * CLQTMC(1) 
     #           + YM * C3QTMC(1))
            CSG(1) = CSG(1) + ASQ * 
     #          ( (YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2QTMC(1) 
     #           -(Y**2D0) * CLQTMC(1))
            CSG(2) = CSG(2) + ASQ * 
     #           ((YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2GTMC(1) 
     #           -(Y**2D0) * CLGTMC(1))
         ENDIF
 
**     sigma antineutrino 

      ELSEIF(OBSERV.EQ.10) THEN          
         IF(ITMC.EQ.1) THEN
            CALL COMPUTE_TMC(ZN,X,Q2,TMC)
            DO J=0,1
               C2QTMC(J) = C2Q(J)* (TMC(1) 
     #              -(Y**2.D0/(YP -(2.d0*(MN*X*Y)**2D0)/Q2))*TMC(2))
               C2GTMC(J) = C2G(J)* (TMC(1) 
     #              -(Y**2.D0/(YP -(2.d0*(MN*X*Y)**2D0)/Q2))*TMC(2))
               C3QTMC(J) = C3Q(J) * TMC(3)
               C3GTMC(J) = C3G(J) * TMC(3)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
             ENDDO    
         ELSEIF (ITMC.EQ.0) THEN
            DO J=0,1
               C2QTMC(J) = C2Q(J)
               C2GTMC(J) = C2G(J)
               C3QTMC(J) = C3Q(J)
               C3GTMC(J) = C3G(J)
               CLQTMC(J) = CLQ(J)
               CLGTMC(J) = CLG(J)
            ENDDO
         ENDIF
         CNS(1) = YM * C3QTMC(0)
         CNS(2) = YM * C3QTMC(0)
         CNS(3) = + (YP -(2.d0*(MN*X*Y)**2D0)/Q2 ) * C2QTMC(0)
     #        - (Y**2D0) * CLQTMC(0) 
     #        - YM * C3QTMC(0)
         CNS(4) = - (YP -(2.D0*(MN*X*Y)**2D0)/Q2 ) * C2QTMC(0)
     #        + (Y**2D0) * CLQTMC(0) 
     #        - YM * C3QTMC(0)
         CSG(1) = (YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2QTMC(0) 
     #        - (Y**2D0) * CLQTMC(0)
         CSG(2) = (YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2GTMC(0) 
     #        - (Y**2D0) * CLGTMC(0)
         IF(IPT.GE.1) THEN
            CNS(1) = CNS(1) + ASQ * (YM * C3QTMC(1))
            CNS(2) = CNS(2) + ASQ * (YM * C3QTMC(1))
            CNS(3) = CNS(3) + ASQ * 
     #          (+ (YP -(2.d0*(MN*X*Y)**2D0)/Q2) * C2QTMC(1)
     #           - (Y**2D0) * CLQTMC(1) 
     #           - YM * C3QTMC(1))
            CNS(4) = CNS(4) + ASQ *
     #         ( - (YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2QTMC(1)
     #           + (Y**2D0) * CLQTMC(1) 
     #           - YM * C3QTMC(1))
            CSG(1) = CSG(1) + ASQ * 
     #           ( (YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2QTMC(1) 
     #           - (Y**2D0) * CLQTMC(1))
            CSG(2) = CSG(2) + ASQ * 
     #           ((YP -(2.D0*(MN*X*Y)**2D0)/Q2) * C2GTMC(1) 
     #           - (Y**2D0) * CLGTMC(1))
         ENDIF
*     
      ENDIF
*     
      RETURN
      END
*
