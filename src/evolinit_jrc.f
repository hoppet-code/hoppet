************************************************************* 
* 
*     File: evolinit.f 
*      
*     Initializes the evolution code by filling the COMMON 
*     BLOCKS which contain the constants
* 
*     AG: 17.01.07 
* 
************************************************************* 

      SUBROUTINE EVOLINIT 
*     -------------------
*    
      IMPLICIT none 
*
      include "../../commons/evflags.h"
      include "../../commons/evscale.h"
      include "../../commons/colfact.h"
      include "../../commons/consts.h"
      include "../../commons/beta.h"
      include "../../commons/alphas.h"
*
      INTEGER I, NFF
      DOUBLE PRECISION alphas,AlphaExactBeta,AlphaMZ
*     
      EXTERNAL AlphaExactBeta,AlphaMZ
*
*     Color Factors
*
      CA= 3D0
      CF= 4D0/3D0
      TR= 0.5D0
*
*     Constants
*
      EMC   = 0.5772156649D0
      ZETA2 = 1.644934067D0  
      ZETA3 = 1.2020569031D0
      ZETA4 = 1.0823232337D0
* 
*     Beta function coefficients 
* 
      DO I=3,6         
         BETA0(I) = (33d0-2d0*I)/3d0  
         BETA1(I) = 102d0-38d0/3d0*I 
         BETA2(I) = 2857d0/2d0-5033d0/18d0*I+325d0/54d0*I**2d0 
         B1(I) = BETA1(I)/BETA0(I)
         B2(I) = BETA2(I)/BETA0(I)
      ENDDO
*
      IF(IVFN.EQ.0) THEN
         NFF=4
         CALL FFN(q20,alphas,ipt,AlphaExactBeta,NFF)
c         CALL FFN(q20,alphas,ipt,AlphaMZ,NFF)
         AS0 = ALPHAS/4d0/PI
         CALL FFN(q2TH(4),alphas,ipt,AlphaExactBeta,NFF)
c         CALL FFN(q2TH(4),alphas,ipt,AlphaMZ,NFF)
         ASC = ALPHAS/4d0/PI
         CALL FFN(q2TH(5),alphas,ipt,AlphaExactBeta,NFF)
c         CALL FFN(q2TH(5),alphas,ipt,AlphaMZ,NFF)
         ASB = ALPHAS/4d0/PI
         CALL FFN(q2TH(6),alphas,ipt,AlphaExactBeta,NFF)
c         CALL FFN(q2TH(6),alphas,ipt,AlphaMZ,NFF)
         AST = ALPHAS/4d0/PI
      ELSE
         CALL VFN(q20,alphas,ipt,AlphaExactBeta)
c         CALL VFN(q20,alphas,ipt,AlphaMZ)
         AS0 = ALPHAS/4d0/PI
         CALL VFN(q2TH(4),alphas,ipt,AlphaExactBeta)
c         CALL VFN(q2TH(4),alphas,ipt,AlphaMZ)
         ASC = ALPHAS/4d0/PI
         CALL VFN(q2TH(5),alphas,ipt,AlphaExactBeta)
c         CALL VFN(q2TH(5),alphas,ipt,AlphaMZ)
         ASB = ALPHAS/4d0/PI
         CALL VFN(q2TH(6),alphas,ipt,AlphaExactBeta)
c         CALL VFN(q2TH(6),alphas,ipt,AlphaMZ)
         AST = ALPHAS/4d0/PI
      ENDIF
*
      ASQ = 0.d0 
*            
      RETURN
      END
*
*--------------------------------------------------------------------
