!  $Id: xpij2n.f90,v 1.3 2004/06/01 09:30:21 salam Exp $ 
!  Automatically generated from f77 file, with addition of "d0"
!  and the placement inside a module.
module xpij2n
character(len=*), parameter :: name_xpij2 = "xpij2n"
contains
!                                                                       
! ..File: xpij2n.f   UPDATE 7/2000                                      
!                                                                       
!                                                                       
! ..Parametrisation of the 3-loop MS(bar) splitting functions P_ij^(2)  
!    for the evolution of unpolarized singlet partons at mu_r = mu_f.   
!    The expansion parameter is alpha_s/(4 pi).                         
!                                                                       
! ..The two sets providing the error estimate are called via  IMOD = 1  
!    and  IMOD = 2.0d0  Any other value of IMOD invokes their average.     
!                                                                       
! ..The distributions (in the mathematical sense) are given as in eq.   
!   (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.0d0   
!   The name-endings A, B, and C of the functions below correspond to   
!   the kernel superscripts [2], [3], and [1] in that equation.         
!                                                                       
! ..The results are based on the lowest six even-integer moments        
!    calculated by Larin et al. and  Retey and Vermaseren, supplemented 
!    by the leading small-x terms obtained by Catani and Hautmann,      
!    Fadin and Lipatov, and Camici and Ciafaloni.                       
!                                                                       
! ..Reference: W.L. van Neerven and A. Vogt,                            
!              hep-ph/9907472  (non-singlet, needed for P_qq),          
!              hep-ph/0006154, and  hep-ph/0007362                      
!    It is appropriate to cite also the sources of information, i.e,    
!    refs. [2-9] of hep-ph/0007362.0d0                                     
!                                                                       
!                                                                       
! ===================================================================== 
!                                                                       
!                                                                       
! ..This is the (regular) pure-singlet splitting functions P_ps^(2).    
!    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+. 
!                                                                       
       FUNCTION P2PSA (Y, NF, IMOD) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       DIMENSION CL3(2) 
       INTEGER IMOD, NF 
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2PSA1 = (229.497d0 * DL1 + 722.99d0 * Y**2 - 2678.77d0 + 560.20d0/Y)    &
     &          * (1.0d0-Y) - 2008.61d0 * DL - 998.15d0 * DL**2                
       G2PSA2 = (-73.845d0 * DL1**2 - 305.988d0 * DL1 - 2063.19d0 * Y         &
     &          + 387.95d0/Y ) * (1.0d0-Y) - 1999.35d0 * Y*DL + 732.68d0 * DL    
!                                                                       
       IF (IMOD .EQ. 1) THEN 
         G2PSA = G2PSA1 
       ELSE IF (IMOD .EQ. 2) THEN 
         G2PSA = G2PSA2 
       ELSE 
         G2PSA = 0.5d0 * (G2PSA1 + G2PSA2) 
       END IF 
!                                                                       
       G2PSA = G2PSA + 3584.0d0/(27.0d0* Y) * DL                              &
     &         + NF * ( (7.282d0 * DL1 + 38.779d0 * Y**2 - 32.022d0 * Y       &
     &                  + 6.252d0 - 1.767d0 / Y ) * (1.0d0-Y) -7.453d0 * DL**2 ) 
       P2PSA = - NF * G2PSA 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the gluon->quark splitting functions P_qg^(2).              
!                                                                       
       FUNCTION P2QGA (Y, NF, IMOD) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       DIMENSION CL3(2) 
       INTEGER IMOD, NF 
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2QGA1 =   31.830d0 * DL1**3 - 1252.267d0 * DL1 - 1722.47d0      &
     &          - 1999.89d0 * Y - 1223.43d0 * DL**2 + 1334.61d0 / Y     
       G2QGA2 = - 19.428d0 * DL1**4 - 159.833d0*DL1**3-309.384d0*DL1**2 &
     &          - 2631.00d0 * (1.0d0-Y) + 67.25d0 * DL**2 + 776.793d0 / Y
!                                                                       
       IF (IMOD .EQ. 1) THEN 
         G2QGA = G2QGA1 
       ELSE IF (IMOD .EQ. 2) THEN 
         G2QGA = G2QGA2 
       ELSE 
         G2QGA = 0.5d0 * (G2QGA1 + G2QGA2) 
       END IF 
!                                                                       
       G2QGA = G2QGA + 896.0d0/ (3.0d0*Y) * DL                            &
     &         + NF * ( 0.9085d0 * DL1**2 + 35.803d0 * DL1 + 128.023d0  &
     &                  - 200.929d0*(1.0d0-Y)-40.542d0*DL-3.284d0 / Y)   
       P2QGA = - NF * G2QGA 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the quark->gluon splitting functions P_gq^(2).              
!                                                                       
       FUNCTION P2GQA (Y, NF, IMOD) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMOD, NF 
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2GQA1 = -13.1212d0 * DL1**4 - 126.665d0 * DL1**3 - 308.536d0 * DL1**2 &
     &          - 361.21d0 + 2113.45d0 * DL + 17.965d0 * DL/Y                 &
     &   + NF* (- 2.4427d0 * DL1**4 - 27.763d0 * DL1**3 - 80.548d0 * DL1**2   &
     &          +227.135d0 + 151.04d0 * DL**2 - 65.91d0  * DL/Y )             
       G2GQA2 = 4.5107d0 * DL1**4 + 66.618d0 * DL1**3 + 231.535d0 * DL1**2    &
     &          + 1224.22d0* (1.0d0-Y) - 240.08d0 * DL**2 - 379.60d0 * (DL+4.0d0)/Y &
     &   + NF* (1.4028d0 * DL1**4 + 11.638d0 * DL1**3 - 164.963d0 * DL1       &
     &          + 1066.78d0* (1.0d0-Y) + 182.08d0 * DL**2 - 138.54d0 * (DL+2.0d0)/Y)
!                                                                       
       IF (IMOD .EQ. 1) THEN 
         G2GQA = G2GQA1 
       ELSE IF (IMOD .EQ. 2) THEN 
         G2GQA = G2GQA2 
       ELSE 
         G2GQA = 0.5d0 * (G2GQA1 + G2GQA2) 
       END IF 
!                                                                       
       G2GQA = G2GQA                                                    &
     &         + NF**2* (- 1.9361d0 * DL1**2 - 11.178d0 * DL1 - 11.632d0      &
     &                   + 15.145d0 * (1.0d0-Y) - 3.354d0 * DL + 2.133d0 / Y )   
       P2GQA = - G2GQA 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the regular piece of the gg splitting function P_gg^(2).    
!                                                                       
       FUNCTION P2GGA (Y, NF, IMOD) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMOD, NF 
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2GGA1 = 732.715d0 * DL1**2 + 20640.069d0 * Y + 15428.58d0 * (1.0d0-Y**2) &
     &          + 15213.6d0 * DL**2 - 16700.88d0 * 1.0d0/Y                     &
     &   + NF* (425.708d0 * DL1 - 914.548d0  * Y**2 + 1122.86d0               &
     &          + 444.21d0 * DL**2 - 376.98d0 * 1.0d0/Y )                      
       G2GGA2 = - 3748.934d0 * DL1 + 35974.45d0 * (1.0d0+Y**2) - 60879.62d0 * Y  &
     &          - 2002.96d0 * DL**2 - 9762.09d0 * 1.0d0/Y                      &
     &   + NF* (- 62.630d0 * DL1**2 - 801.90d0 - 1891.40d0 * DL               &
     &          - 813.78d0 * DL**2 - 1.360d0 * 1.0d0/Y )                       
!                                                                       
       IF (IMOD .EQ. 1) THEN 
         G2GGA = G2GGA1 
       ELSE IF (IMOD .EQ. 2) THEN 
         G2GGA = G2GGA2 
       ELSE 
         G2GGA = 0.5d0 * (G2GGA1 + G2GGA2) 
       END IF 
!                                                                       
       G2GGA = G2GGA - (2675.85d0 + NF * 157.18d0) * DL / Y                 &
     &         + NF**2 * (- 37.6417d0 * Y**2 + 72.926d0 * Y - 32.349d0        &
     &                    + 0.991d0 * DL**2 - 2.818d0 / Y )                 
       P2GGA = - G2GGA 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the singular piece of the gg splitting function P_gg^(2).   
!                                                                       
       FUNCTION P2GGB (Y, NF, IMOD) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMOD, NF 
       D1 = 1.0d0/(1.0d0-Y) 
!                                                                       
       G2GGB1 = (-2626.38d0 + NF * 415.71d0) * D1 
       G2GGB2 = (-2678.22d0 + NF * 412.00d0) * D1 
!                                                                       
       IF (IMOD .EQ. 1) THEN 
         G2GGB = G2GGB1 
       ELSE IF (IMOD .EQ. 2) THEN 
         G2GGB = G2GGB2 
       ELSE 
         G2GGB = 0.5d0 * (G2GGB1 + G2GGB2) 
       END IF 
!                                                                       
       G2GGB =  G2GGB + NF**2 * 16.0d0/9.0d0 * D1 
       P2GGB = -G2GGB 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the 'local' piece of the gg splitting function P_gg^(2).    
!                                                                       
       FUNCTION P2GGC (Y, NF, IMOD) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMOD, NF 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2GGC1 = (-2626.38d0 + NF * 415.71d0)* DL1 - 4424.168d0 + NF * 548.569d0 
       G2GGC2 = (-2678.22d0 + NF * 412.00d0)* DL1 - 4590.570d0 + NF * 534.951d0 
!                                                                       
       IF (IMOD .EQ. 1) THEN 
         G2GGC = G2GGC1 
       ELSE IF (IMOD .EQ. 2) THEN 
         G2GGC = G2GGC2 
       ELSE 
         G2GGC = 0.5d0 * (G2GGC1 + G2GGC2) 
       END IF 
!                                                                       
       G2GGC =  G2GGC + NF**2 * (16.0d0/9.0d0 * DL1 - 6.4882d0) 
       P2GGC = -G2GGC 
!                                                                       
       RETURN 
      END FUNCTION
end module xpij2n
