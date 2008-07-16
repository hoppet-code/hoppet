**********************************************************************
*
*     File: wcDIS.f
*
*     Returns the C2, C3, CL Coefficient Functions coefficients 
*     in Mellin space.
*     
*     MU: 18.05.2007
*
**********************************************************************
*
      SUBROUTINE C2N(N,NF,C2Q,C2G)
*
      IMPLICIT none
*
      include "../../commons/colfact.h"
      include "../../commons/consts.h"
*
      INTEGER NF
      DOUBLE COMPLEX N
      DOUBLE COMPLEX C2Q(0:2),C2G(0:2)
      DOUBLE COMPLEX S1,S2,N1,N2,NS,NM
      DOUBLE COMPLEX PSI,DPSI
*
      S1 = EMC + PSI(N+1d0)
      S2 = ZETA2 - DPSI(N+1d0,1)
*     
      NS = N * N
      NM = N - 1d0
      N1 = N + 1d0
      N2 = N + 2d0
*
      C2Q(0)= 1.D0
      C2G(0)= 0.D0
*
      C2Q(1) =  CF * ( 2D0*S1**2D0 - 2D0*S2 + 3D0*S1 -2*S1/(N*N1)
     #         + 3D0/N + 4D0/N1 + 2D0/NS - 9D0 )
      C2G(1)=   DBLE(NF) * (4D0 * TR * ( 4D0/N1 - 4D0/N2 
     #         - (1D0+S1) * (NS+N+2D0)/(N*N1*N2) + 1D0/NS ))

      RETURN
      END

      
      SUBROUTINE C3N(N,NF,C3Q,C3G)
*
      IMPLICIT none
*
      include "../../commons/colfact.h"
      include "../../commons/consts.h"
*
      INTEGER NF
      DOUBLE COMPLEX C3Q(0:2),C3G(0:2)
      DOUBLE COMPLEX N
      DOUBLE COMPLEX S1,S2,N1,N2,NS,NM
      DOUBLE COMPLEX PSI,DPSI
*
      S1 = EMC + PSI(N+1d0)
      S2 = ZETA2 - DPSI(N+1d0,1)
*     
      NS = N * N
      NM = N - 1d0
      N1 = N + 1d0
      N2 = N + 2d0
*
      C3Q(0) = 1D0
      C3G(0) = 0D0
*
      C3Q(1) =  CF * ( 4D0*S1**2D0 - 4D0*S2 + 3D0*S1 -2*S1/(N*N1)
     #        + 3D0/N + 4D0/N1 + 2D0/NS - 9D0 - (4D0*N+2D0)/(N*N1) )
      C3G(1) = 0.D0

      RETURN
      END


      SUBROUTINE CLN(N,NF,CLQ,CLG)
*
      IMPLICIT none
*
      include "../../commons/colfact.h"
      include "../../commons/consts.h"
*
      INTEGER NF
      DOUBLE COMPLEX CLQ(0:2),CLG(0:2)
      DOUBLE COMPLEX N
      DOUBLE COMPLEX S1,S2,N1,N2,NS,NM
      DOUBLE COMPLEX PSI,DPSI

      S1 = EMC + PSI(N+1d0)
      S2 = ZETA2 - DPSI(N+1d0,1)
*     
      NS = N * N
      NM = N - 1d0
      N1 = N + 1d0
      N2 = N + 2d0

      CLQ(0) = 0D0
      CLG(0) = 0D0

      CLQ(1) =  CF * 4D0/N1
      CLG(1) =  DBLE(NF) * 16D0 * TR /(N1*N2)

      RETURN
      END
