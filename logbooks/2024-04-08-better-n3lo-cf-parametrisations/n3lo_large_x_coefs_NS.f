! This file contains the large x expansions of the non-singlet piece of
! the C2, C3, and CL coefficient functions. We retain all terms
! log[1-x]^i where i>0. Note that this is the regular piece, so not plus
! distribution and we also only include the charged current contribution
! as the piece multiplying fl011 is not divergent. Equations are taken
! from 0812.4168 and hep-ph/0504242 with some corrections. Please cite
! arXiv:2404.XXXXX if you use this file.
!
! Expressions are exact for QCD colour factors but everything is
! truncated at double precision.
! Arguments are x, DL = log(x), NF = number of active flavours
! DL1 is computed and is log(1-x)    
      FUNCTION C2NP3A_large_x (X, DL, NF)
      IMPLICIT REAL*8 (A - Z)
      INTEGER NF
      DOUBLE PRECISION X, DL, DL1
      
      DL1 = DL1VAL(X, DL)
      
      C2NP3A_large_x = 5894.634952596363d0*DL1 + 2319.655820717043d0
     $     *DL1**2 -1787.0418273217867d0*DL1**3 + 348.44444444444446d0
     $     *DL1**4 -18.962962962962962d0*DL1**5 +(1199.690656381538d0
     $     *DL1 -787.5420539087113d0*DL1**2 +146.69958847736626d0*DL1
     $     **3 -7.901234567901234d0*DL1**4) *NF +(-65.15652656374832d0
     $     *DL1 +14.617283950617283d0*DL1 **2 -0.7901234567901234d0*DL1
     $     **3)*NF**2
      
      RETURN
      END FUNCTION
      
      FUNCTION C3NM3A_large_x (X, DL, NF)
      IMPLICIT REAL*8 (A - Z)
      INTEGER NF
      
      DL1 = DL1VAL(X, DL)
      
      C3NM3A_large_x = 6889.89378009207d0*DL1 + 1581.208087862925d0*DL1
     $     **2 -1609.6420585153533d0*DL1**3 + 329.48148148148147d0*DL1
     $     **4 -18.962962962962962d0*DL1**5 +(859.380948706157d0*DL1 -
     $     675.1899801124716d0*DL1**2 +134.0576131687243d0*DL1**3 -
     $     7.901234567901234d0*DL1**4)*NF +(-50.144180884735974d0*DL1 +
     $     12.246913580246913d0*DL1**2 -0.7901234567901234d0*DL1**3)*NF
     $     **2
      
      RETURN
      END FUNCTION
      
      FUNCTION CLNP3A_large_x (X, DL, NF)
      IMPLICIT REAL*8 (A - Z)
      INTEGER NF
      
      DL1 = DL1VAL(X, DL)
      
      CLNP3A_large_x = -995.2588274957099d0*DL1 + 738.4477328541179d0
     $     *DL1**2 -177.39976880643366d0*DL1**3 + 18.962962962962962d0
     $     *DL1**4+(340.3097076753812d0*DL1 - 112.35207379623978d0*DL1
     $     **2+12.641975308641975d0*DL1**3)*NF +(-15.012345679012345d0
     $     *DL1 +2.3703703703703702d0*DL1**2)*NF**2
      
      RETURN
      END FUNCTION
      

* ..For X values close to 1, use the series expansion close 
*   to 1-X=0 instead of full value, for numerical convergence
*     
      FUNCTION DL1VAL (X, DL)
      IMPLICIT REAL*8 (A - Z)
      DOUBLE PRECISION D2, D24, D2880
      PARAMETER (D2=0.5D0, D24=1.0D0/24.0D0, D2880=1.0D0/2880.0D0)
      
      IF (ABS(DL).LT.1D-4) THEN
         DL1VAL = LOG(-DL) +  DL*D2 + DL**2*D24 - DL**4*D2880
      ELSE
         DL1VAL = LOG(1.0D0 - X)
      ENDIF
      
      RETURN
      END FUNCTION
      
