!  $Id: xpns2n.f90,v 1.3 2004/06/01 09:30:28 salam Exp $ 
!  Automatically generated from f77 file, with addition of "d0"
!  and the placement inside a module.
module xpns2n
character(len=*), parameter :: name_xpns2 = "xpns2n"
contains
!                                                                       
! ..File: xpns2n.f   UPDATE 7/2000                                      
!                                                                       
!                                                                       
! ..Parametrization of the 3-loop MS(bar) splitting functions P_NS^(2)  
!    for the evolution of unpolarized non-singlet partons, mu_r = mu_f. 
!    The expansion parameter is alpha_s/(4 pi).                         
!                                                                       
! ..The two sets spanning the error estimate are called via  IMODN = 1  
!    and  IMODN = 2.0d0  Any other value of IMODN invokes their average.   
!                                                                       
! ..The distributions (in the mathematical sense) are given as in eq.   
!    (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.0d0  
!    The name-endings A, B, and C of the functions below correspond to  
!    the kernel superscripts [2], [3], and [1] in that equation.        
!                                                                       
! ..The P^+ results are based on the lowest six even-integer moments    
!    calculated by Larin et al. and  Retey and Vermaseren. The P^- and  
!    P^S = P^V - P^-) approximations use the seven lowest odd-integer   
!    moments computed also by Retey and Vermaseren. Also used for P^+   
!    and P^- are the leading small-x terms (ln^4 x) as obtained by      
!    Blumlein and Vogt. The exact N_f^2 term has been derived by Gracey.
!                                                                       
! ..Reference: W.L. van Neerven and A. Vogt,                            
!              hep-ph/9907472 and hep-ph/0007362                        
!    It is appropriate to cite also the above-mentioned sources of      
!    information, i.e, ref. [2-4] and [6,7] of hep-ph/0007362.0d0          
!                                                                       
!                                                                       
! ===================================================================== 
!                                                                       
!                                                                       
! ..This is the regular piece of P^(2)+.                                
!                                                                       
       FUNCTION P2NSPA (Y, NF, IMODN) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMODN, NF 
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2NSPA1 = - 1047.590d0 * DL1 + 843.884d0 * Y**2  + 98.65d0 * (1.0d0-Y)    &
     &           + 33.71d0 * DL**2  - 1.580d0 * DL**3 * (DL + 4.0d0)           &
     &    + NF* (- 9.649d0 * DL1**2 - 406.171d0 * Y**2  - 32.218d0 * (1.0d0-Y)   &
     &           - 5.976d0 * DL**2  - 1.60d0 * DL**3 )                      
       G2NSPA2 = 147.692d0 * DL1**2 + 2602.738d0 * Y**2 + 170.11d0            &
     &           - 148.47d0 * DL    - 1.580d0 * DL**3 * (DL - 4.0d0)           &
     &     + NF* ( 89.941d0 * DL1   - 218.482d0 * Y**2  - 9.623d0             &
     &           - 0.910d0* DL**2   + 1.60d0 * DL**3 )                      
!                                                                       
       IF (IMODN .EQ. 1) THEN 
         G2NSPA = G2NSPA1 
       ELSE IF (IMODN .EQ. 2) THEN 
         G2NSPA = G2NSPA2 
       ELSE 
         G2NSPA = 0.5d0 * (G2NSPA1 + G2NSPA2) 
       END IF 
!                                                                       
       G2NSPA = G2NSPA                                                  &
     &     + NF**2* (- 32.0d0* Y*DL/(1.0d0-Y) * (3.0d0* DL + 10.0d0) - 64.0d0          &
     &               - (48.0d0* DL**2 + 352.0d0* DL + 384.0d0) * (1.0d0-Y) )/81.0d0  
       P2NSPA = -G2NSPA 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
! ..This is the regular piece of P^(2)-.                                
!                                                                       
       FUNCTION P2NSMA (Y, NF, IMODN) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMODN, NF 
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2NSMA1 =  157.387d0 * DL1**2 + 2741.42d0 * Y**2 + 490.43d0 * (1.0d0-Y)   &
     &           - 67.0d0   * DL**2  - 10.005d0 * DL**3 - 1.432d0 * DL**4     &
     &    + NF *(- 17.989d0 * DL1**2 - 355.636d0 * Y**2 + 73.407d0 *(1.0d0-Y)*DL1&
     &           - 11.491d0 * DL**2  - 1.928d0 * DL**3 )                    
       G2NSMA2 = -115.099d0 * DL1**2 - 1581.05d0 * DL1  - 267.33d0  * (1.0d0-Y)  &
     &           + 127.65d0 * DL**2  + 25.22d0 * DL**3  - 1.432d0 * DL**4     &
     &    + NF *(- 11.999d0 * DL1**2 - 397.546d0 * Y**2 - 41.949d0 * (1.0d0-Y)   &
     &           + 1.477d0  * DL**2  + 0.538d0 * DL**3 )                    
!                                                                       
       IF (IMODN .EQ. 1) THEN 
         G2NSMA = G2NSMA1 
       ELSE IF (IMODN .EQ. 2) THEN 
         G2NSMA = G2NSMA2 
       ELSE 
         G2NSMA = 0.5d0 * (G2NSMA1 + G2NSMA2) 
       END IF 
!                                                                       
       G2NSMA =  G2NSMA                                                 &
     &     + NF**2* (- 32.0d0* Y*DL/(1.0d0-Y) * (3.0d0* DL + 10.0d0) - 64.0d0          &
     &               - (48.0d0* DL**2 + 352.0d0* DL + 384.0d0) * (1.0d0-Y) )/81.0d0  
       P2NSMA = -G2NSMA 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
! ..This is the singular NS piece for the `+' case.                     
!                                                                       
       FUNCTION P2NSPB (Y, NF, IMODN) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMODN, NF 
       D1 = 1.0d0/(1.0d0-Y) 
!                                                                       
       G2NSB1 = (-1183.762d0 + NF * 183.148d0) * D1 
       G2NSB2 = (-1182.774d0 + NF * 183.931d0) * D1 
!                                                                       
       IF (IMODN .EQ. 1) THEN 
         G2NSB = G2NSB1 
       ELSE IF (IMODN .EQ. 2) THEN 
         G2NSB = G2NSB2 
       ELSE 
         G2NSB = 0.5d0 * (G2NSB1 + G2NSB2) 
       END IF 
!                                                                       
       G2NSB  =  G2NSB + NF**2 * 64.0d0/81.0d0 * D1 
       P2NSPB = -G2NSB 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
! ..This is the singular NS piece for the `-' case.                     
!                                                                       
       FUNCTION P2NSMB (Y, NF, IMODN) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       DIMENSION CDM0(2), CDM1(2) 
       INTEGER IMODN, NF 
       DATA CDM0 / -1185.229d0, -1174.348d0 / 
       DATA CDM1 /  184.765d0,  183.718d0 / 
       D1 = 1.0d0/(1.0d0-Y) 
!                                                                       
       G2NSB1 = (-1185.229d0 + NF * 184.765d0) * D1 
       G2NSB2 = (-1174.348d0 + NF * 183.718d0) * D1 
!                                                                       
       IF (IMODN .EQ. 1) THEN 
         G2NSB = G2NSB1 
       ELSE IF (IMODN .EQ. 2) THEN 
         G2NSB = G2NSB2 
       ELSE 
         G2NSB = 0.5d0 * (G2NSB1 + G2NSB2) 
       END IF 
!                                                                       
       G2NSB  =  G2NSB + NF**2 * 64.0d0/81.0d0 * D1 
       P2NSMB = -G2NSB 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
! ..This is the 'local' NS+ piece piece for the `+' case.               
!                                                                       
       FUNCTION P2NSPC (Y, NF, IMODN) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMODN, NF 
       DATA Z2, Z3 / 1.644934067d0, 1.202056903d0 / 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2NSC1 = (-1183.762d0 + NF* 183.148d0)* DL1 - 1347.032d0 + NF* 174.402d0 
       G2NSC2 = (-1182.774d0 + NF* 183.931d0)* DL1 - 1351.088d0 + NF* 178.208d0 
!                                                                       
       IF (IMODN .EQ. 1) THEN 
         G2NSC = G2NSC1 
       ELSE IF (IMODN .EQ. 2) THEN 
         G2NSC = G2NSC2 
       ELSE 
         G2NSC = 0.5d0 * (G2NSC1 + G2NSC2) 
       END IF 
!                                                                       
       G2NSC  =  G2NSC                                                  &
     &    + NF**2 * (16.0d0* DL1 + (51.0d0 + 48.0d0* Z3 - 80.0d0* Z2)) * 4.0d0/81.0d0   
       P2NSPC = -G2NSC 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
! ..This is the 'local' NS+ piece piece for the `-' case.               
!                                                                       
       FUNCTION P2NSMC (Y, NF, IMODN) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMODN, NF 
       DATA Z2, Z3 / 1.644934067d0, 1.202056903d0 / 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       G2NSC1 = (-1185.229d0 + NF* 184.765d0)* DL1 - 1365.458d0 + NF* 184.289d0 
       G2NSC2 = (-1174.348d0 + NF* 183.718d0)* DL1 - 1286.799d0 + NF* 177.762d0 
!                                                                       
       IF (IMODN .EQ. 1) THEN 
         G2NSC = G2NSC1 
       ELSE IF (IMODN .EQ. 2) THEN 
         G2NSC = G2NSC2 
       ELSE 
         G2NSC = 0.5d0 * (G2NSC1 + G2NSC2) 
       END IF 
!                                                                       
       G2NSC  =  G2NSC                                                  &
     &    + NF**2 * (16.0d0* DL1 + (51.0d0 + 48.0d0* Z3 - 80.0d0* Z2)) * 4.0d0/81.0d0   
       P2NSMC = -G2NSC 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
! ..This is P^(2)S, the difference of P^(2)V and P^(2)-.                
!                                                                       
       FUNCTION P2NSSA (Y, NF, IMODN) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER IMODN, NF 
       DL  = LOG (Y) 
!                                                                       
       G2NSSA1 = (  1441.57d0 * Y**2 - 12603.59d0 * Y + 15450.01d0) * (1.0d0-Y)  &
     &            - 7876.93d0 * Y*DL**2 + 4260.29d0 * DL + 229.27d0 * DL**2   &
     &            - 4.4075d0 * DL**3                                      
       G2NSSA2 = (  704.67d0 * Y**3 - 3310.32d0  * Y**2 - 2144.81d0  * Y      &
     &            + 244.68d0) * (1.0d0-Y) - 4490.81d0 * Y**2 * DL              &
     &            - 42.875d0 * DL + 11.0165d0 * DL**3                       
!                                                                       
       IF (IMODN .EQ. 1) THEN 
         G2NSSA = G2NSSA1 
       ELSE IF (IMODN .EQ. 2) THEN 
         G2NSSA = G2NSSA2 
       ELSE 
         G2NSSA = 0.5d0 * (G2NSSA1 + G2NSSA2) 
       END IF 
!                                                                       
       P2NSSA  = - NF * G2NSSA 
!                                                                       
       RETURN 
      END FUNCTION
end module xpns2n
