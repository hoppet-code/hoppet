*************************************************************
*
*     File: andim_nlo.f
*     
*     Returns the NLO anomalous dimensions in the 
*     a_s=alpha_s/4pi notation
*
*     Input variables: N  : Mellin variable
*                      NF : Number of flavours
*
*     Functions: PSI from data/theo/harmsums/psi.f
*
*     Output: P1NS : Complex valued Non-singlet LO anomalous
*                    dimension in N space
*
*             P1SG(2,2) : matrix of the complex valued LO
*                         singlet anomalous dimensions in 
*                         N space
*
*************************************************************

      SUBROUTINE ANDIM_NLO(N,NF,P1NS,P1SG) 
      IMPLICIT none
*
      include "../../commons/colfact.h"
      include "../../commons/consts.h"
*
      DOUBLE COMPLEX DPSI,PSI, S1, S2
*
* ---------------------------------------------------------------------
*
*     Internal variables
*
      INTEGER I,J
      DOUBLE COMPLEX NS,NT,NFO,NFI,NSI,NSE,NE,NN
      DOUBLE COMPLEX N1,N2,NM,NMS,N1S,N1T,N2S,N2T
      DOUBLE COMPLEX N3,N4,N5,N6
      DOUBLE COMPLEX S11,S12,S13,S14,S15,S16
      DOUBLE COMPLEX SPMOM,SLC,SLV,SSCHLM,SSTR2M,SSTR3M,SSCHLP
      DOUBLE COMPLEX SSTR2P,SSTR3P
      DOUBLE COMPLEX PPSA,PQGA,PGQA,PGGA,PQGB,PGQB,PGGB,PGQC,PGGC
      DOUBLE COMPLEX PNPA,PNMA,PNSB,PNSC
*
* ---------------------------------------------------------------------
*
*     Input variables
*
       DOUBLE COMPLEX N
       INTEGER NF       
*
*     Output variables  
*
       DOUBLE COMPLEX P1NS(3)
       DOUBLE COMPLEX P1SG(2,2)
*
* ---------------------------------------------------------------------
*     

      S1 = EMC + PSI(N+1d0)
      S2 = ZETA2 - DPSI(N+1d0,1)
*     
      NS = N * N
      NT = NS * N
      NFO = NT * N
      NFI = NFO * N
      NSI = NFI * N
      NSE = NSI * N
      NE = NSE * N
      NN = NE * N
*     
      NM = N - 1d0
      N1 = N + 1d0
      N2 = N + 2d0
      NMS = NM * NM
      N1S = N1 * N1
      N1T = N1S * N1
      N2S = N2 * N2
      N2T = N2S * N2

*
* ---------------------------------------------------------------------
*
* ..Analytic continuations of the occuring sums as given in GRV (1990) 
*   (with an improved parametrization of the moments of  Sp(x)/(1+x).)
*

*
       N3 = N + 3D0
       N4 = N + 4D0
       N5 = N + 5D0
       N6 = N + 6D0
       S11 = S1  + 1D0/N1
       S12 = S11 + 1D0/N2
       S13 = S12 + 1D0/N3
       S14 = S13 + 1D0/N4
       S15 = S14 + 1D0/N5
       S16 = S15 + 1D0/N6
       SPMOM = 1.0000D0 * (ZETA2 - S1 / N ) / N  -
     1         0.9992D0 * (ZETA2 - S11/ N1) / N1 +
     2         0.9851D0 * (ZETA2 - S12/ N2) / N2 -
     3         0.9005D0 * (ZETA2 - S13/ N3) / N3 +
     4         0.6621D0 * (ZETA2 - S14/ N4) / N4 -
     5         0.3174D0 * (ZETA2 - S15/ N5) / N5 +
     6         0.0699D0 * (ZETA2 - S16/ N6) / N6  
*
       SLC = - 5D0/8D0 * ZETA3
       SLV = - ZETA2/2D0* (PSI(N1/2D0) - PSI(N/2D0)) 
     1 + S1/NS + SPMOM
       SSCHLM = SLC - SLV
       SSTR2M = ZETA2 - DPSI (N1/2D0,1)
       SSTR3M = 0.5D0 * DPSI (N1/2D0,2) + ZETA3

       SSCHLP = SLC + SLV
       SSTR2P = ZETA2 - DPSI(N2/2D0,1)
       SSTR3P = 0.5D0 * DPSI(N2/2D0,2) + ZETA3
       
*
* ---------------------------------------------------------------------
*
* ..The contributions to P1NS as given in Gonzalez-Arroyo et al. (1979) 
*   (Note that the anomalous dimensions in the literature often differ 
*    from these moments of the splitting functions by factors -1 or -2,
*    in addition to possible different normalizations of the coupling)
*

       PNMA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +
     1      16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2M ) +
     2      64D0* SSCHLM + 24D0* S2 - 3D0 - 8D0* SSTR3M -
     3      8D0* (3D0* NT + NS -1D0) / (NT * N1T) +
     4      16D0* (2D0* NS + 2D0* N +1D0) / (NT * N1T) ) * (-0.5D0)
       PNPA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +
     1      16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2P ) +
     2      64D0* SSCHLP + 24D0* S2 - 3D0 - 8D0* SSTR3P -
     3      8D0* (3D0* NT + NS -1D0) / (NT * N1T) -
     4      16D0* (2D0* NS + 2D0* N +1D0)/(NT * N1T) ) * (-0.5D0)
       
       PNSB = ( S1 * (536D0/9D0 + 8D0* (2D0* N + 1D0) / (NS * N1S)) -
     1      (16D0* S1 + 52D0/3D0- 8D0/(N * N1)) * S2 - 43D0/6D0 -
     2      (151D0* NFO + 263D0* NT + 97D0* NS + 3D0* N + 9D0) *
     3      4D0/ (9D0* NT * N1T) ) * (-0.5D0)
       PNSC = ( -160D0/9D0* S1 + 32D0/3.* S2 + 4D0/3D0 +
     1      16D0*(11D0*NS+5D0*N-3D0)/(9D0* NS * N1S))*(-0.5D0)
       
*     
*     ---------------------------------------------------------------------
*
*     ..The contributions to P1SG as given in Floratos et al. (1981) 
*     ..Pure singlet (PS) and QG
*     
       PPSA = (5d0* NFI + 32d0* NFO + 49d0* NT+38d0* NS + 28d0* N + 8d0) 
     #      / (NM * NT * N1T * N2S) * 2d0
*     
       PQGA = (-2d0* S1 * S1 + 2d0* S2 - 2d0* SSTR2P) 
     #      * (NS + N + 2d0) / (N * N1 * N2) 
     #      + (8d0* S1 * (2d0* N + 3d0)) / (N1S * N2S)
     #      + 2d0* (NN + 6d0* NE + 15d0* NSE + 25d0* NSI + 36d0* NFI
     #      + 85d0* NFO + 128d0* NT + 104d0* NS + 64d0* N + 16d0)
     #      / (NM * NT * N1T * N2T)
       PQGB = (2d0* S1 * S1 - 2d0* S2 + 5d0) * (NS + N + 2d0)
     #      / (N * N1 * N2) - 4d0* S1 / NS
     #      + (11d0* NFO + 26d0* NT + 15d0* NS + 8d0* N + 4d0)
     #      / (NT * N1T * N2) 

*     
*     ---------------------------------------------------------------------
*     
*     ..GQ and GG
*     
       PGQA = (- S1 * S1 + 5d0* S1 - S2) * (NS + N + 2d0) 
     #      / (NM * N * N1)  -  2d0* S1 / N1S
     #      - (12d0* NSI + 30d0* NFI + 43d0* NFO + 28d0* NT - NS
     #      - 12d0* N - 4d0) / (2d0* NM * NT * N1T) 
       PGQB = (S1*S1 + S2 - SSTR2P) * (NS + N + 2d0) / (NM * N * N1)
     #      - S1 * (17d0* NFO + 41d0* NS - 22d0* N - 12d0) 
     #      / (3d0* NMS * NS * N1)
     #      + (109d0* NN + 621d0* NE + 1400d0* NSE + 1678d0* NSI
     #      + 695d0* NFI - 1031d0* NFO - 1304d0* NT - 152d0* NS
     #      + 432d0* N + 144d0) / (9d0* NMS * NT * N1T * N2S)
       PGQC = (S1 - 8d0/3d0) * (NS + N + 2d0) / (NM * N * N1) + 1d0/ N1S
       PGQC = 4d0/3d0* PGQC
*     
       PGGA = - (2d0* NFI + 5d0* NFO + 8d0* NT + 7d0* NS- 2d0* N - 2d0)
     #      * 8d0* S1 / (NMS * NS * N1S * N2S) -  67d0/9d0* S1 + 8d0/3d0
     #      - 4d0* SSTR2P * (NS + N + 1d0) / (NM * N * N1 * N2)
     #      + 2d0* S1 * SSTR2P - 4d0* SSCHLP + 0.5d0 * SSTR3P
     #      + (457d0* NN + 2742d0* NE + 6040d0* NSE + 6098d0* NSI
     #      + 1567d0* NFI - 2344d0* NFO - 1632d0* NT + 560d0* NS
     #      + 1488d0* N + 576d0) / (18d0* NMS * NT * N1T * N2T)
       PGGB = (38d0* NFO + 76d0* NT + 94d0* NS + 56d0* N + 12d0) *(-2d0)
     #      / (9d0* NM * NS * N1S * N2)  +  20d0/9d0* S1  -  4d0/3d0
       PGGC = (2d0* NSI + 4d0* NFI + NFO - 10d0* NT - 5d0* NS - 4d0* N
     #      - 4d0) * (-2d0) / (NM * NT * N1T * N2)  -  1d0

*
*     Output to the array
*

       DO I=1,2
          DO J=1,2
             P1SG(I,J)=0d0
          ENDDO
       ENDDO
       
       do I=1,3
          P1NS(I)=0d0
       enddo
       
      
*     
*     NON SINGLET
*     

*     Plus 
      
       P1NS(1) = CF *((CF-CA/2d0)* PNPA + CA* PNSB + TR*dble(NF)* PNSC)
      
*     Minus=Valence
      
       P1NS(2) = CF *((CF-CA/2d0)* PNMA + CA* PNSB + TR*dble(NF)* PNSC)
       P1NS(3) = P1NS(2)
       
*     
*     SINGLET
*     
      
       P1SG(1,1) = P1NS(1) + TR*dble(NF)*CF*PPSA*4d0
       P1SG(1,2) = TR*dble(NF) * (CA * PQGA + CF * PQGB)*4d0
       P1SG(2,1) = (CF*CF*PGQA + CF*CA*PGQB+TR*dble(NF)*CF*PGQC)*4d0
       P1SG(2,2) = (CA*CA*PGGA + TR*dble(NF)*(CA*PGGB+CF*PGGC))*4d0

*
* ---------------------------------------------------------------------
*
*
       RETURN
       END
*
* =================================================================av==
