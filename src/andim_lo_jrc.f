*************************************************************
*
*     File: andim_lo.f
*     
*     Returns the LO anomalous dimensions in the 
*     a_s=alpha_s/4pi notation
*
*     Input variables: N  : Mellin variable
*                      NF : Number of flavours
*
*     Functions: PSI from data/theo/harmsums/psi.f
*
*     Output: P0NS : Complex valued Non-singlet LO anomalous
*                    dimension in N space
*
*             P0SG(2,2) : matrix of the complex valued LO
*                         singlet anomalous dimensions in 
*                         N space
*
*     AG: 17.01.07
*
*************************************************************

       SUBROUTINE ANDIM_LO(N,NF,P0NS,P0SG)
*
       IMPLICIT none
*
       include "../../commons/colfact.h"
       include "../../commons/consts.h"
*
       DOUBLE COMPLEX PSI, S1
*
* ---------------------------------------------------------------------
*
*     Internal variables
*
       DOUBLE COMPLEX NS,N1,N2,NM
       DOUBLE COMPLEX PQQA,PQGA,PGQA,PGGA,PGGB
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
       DOUBLE COMPLEX P0NS 
       DOUBLE COMPLEX P0SG (2,2)
*
* ---------------------------------------------------------------------
*
       NS = N * N
       N1 = N + cmplx(1d0,0d0)
       N2 = N + (2d0,0d0)
       NM = N - (1d0,0d0)

       !write(6,*) NS,N1,CF,CA,NF
*       
       S1 = cmplx(EMC,0d0) + PSI(N1)

*     
       PQQA = (3d0,0d0) - 4d0* S1 + 2d0/(N * N1)
       PQGA = 4d0* (NS + N + 2d0) / (N * N1 * N2)
       PGQA = 2d0 * (NS + N + 2d0) / (N * N1 * NM)
       PGGA = 11d0/3D0 - 4d0* S1 + 4d0/(N * NM) + 4d0/(N1 * N2) 
       PGGB = - 4d0/3D0
*     
*     Output to the array
*     

       P0NS      = CF * PQQA
*
       P0SG(1,1) = CF * PQQA
       P0SG(1,2) = TR * dble(NF) * PQGA
       P0SG(2,1) = CF * PGQA
       P0SG(2,2) = CA * PGGA + TR *dble(NF) * PGGB
*
       RETURN
       END
*
*----------------------------------------------------------------------

