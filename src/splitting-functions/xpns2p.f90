!  $Id: xpns2p.f90,v 1.1 2004/06/01 09:30:31 salam Exp $ 
!  Automatically generated from f77 file, with addition of "d0"
!  and the placement inside a module.
module xpns2p
character(len=*), parameter :: name_xpns2 = "xpns2p"
contains
!                                                                       
! ..File: xpns2p.f                                                      
!                                                                       
!                           __                                          
! ..The parametrized 3-loop MS non-singlet splitting functions P^(2)    
!    for the evolution of unpolarized partons densities, mu_r = mu_f.   
!    The expansion parameter is alpha_s/(4 pi).                         
!                                                                       
! ..The distributions (in the mathematical sense) are given as in eq.   
!    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.0d0  
!    The name-endings A, B, and C of the functions below correspond to  
!    the kernel superscripts [2], [3], and [1] in that equation.        
!                                                                       
!  ..The relative accuracy of these parametrizations, as well as of     
!    the convolution results, is better than one part in thousand.      
!                                                                       
! ..References: S. Moch, J. Vermaseren and A. Vogt,                     
!               hep-ph/0209100 = Nucl. Phys. B646 (2002) 181,           
!               hep-ph/0403192 (submitted to Nucl. Phys. B)             
!                                                                       
! ===================================================================== 
!                                                                       
!                                                                       
! ..This is the regular piece of P2_NS+.  The rational coefficients are 
!    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.0d0 
!    The N_f^2 part is exact and was first determined in N-space by     
!    J.A. Gracey in Phys. Lett. B322 (1994) 141.0d0                        
!                                                                       
       FUNCTION P2NSPA (Y, NF) 
       IMPLICIT REAL*8 (A - Z) 
       INTEGER NF 
!                                                                       
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
       D81 = 1.0d0/81.0d0 
!                                                                       
       P2NSPA =   1641.1d0 - 3135.0d0* Y + 243.6d0 * Y**2 - 522.1d0 * Y**3       &
     &            + 128.0d0*D81 * DL**4 + 2400.0d0*D81 * DL**3                &
     &            + 294.9d0 * DL**2 + 1258.0d0* DL                           &
     &            + 714.1d0 * DL1 + DL*DL1 * (563.9d0 + 256.8d0 * DL)         &
     &        + NF * ( -197.0d0 + 381.1d0 * Y + 72.94d0 * Y**2 + 44.79d0 * Y**3 &
     &            - 192.0d0*D81 * DL**3  - 2608.0d0*D81 * DL**2 - 152.6d0 * DL  &
     &            - 5120.0d0*D81 * DL1 - 56.66d0 * DL*DL1 - 1.497d0 * Y*DL**3 )&
     &        + NF**2 * ( 32.0d0* Y*DL/(1.0d0-Y) * (3.0d0* DL + 10.0d0) + 64.0d0       &
     &            + (48.0d0* DL**2 + 352.0d0* DL + 384.0d0) * (1.0d0-Y) ) * D81     
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the regular piece of P2_NS-.  The rational coefficients are 
!    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.0d0 
!    The N_f^2 part is exact (and identical to that of P2_NS+).         
!                                                                       
       FUNCTION P2NSMA (Y, NF) 
       IMPLICIT REAL*8 (A - Z) 
       INTEGER NF 
!                                                                       
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
       D81 = 1.0d0/81.0d0 
!                                                                       
       P2NSMA =   1860.2d0 - 3505.0d0* Y + 297.0d0 * Y**2 - 433.2d0 * Y**3       &
     &            + 116.0d0*D81 * DL**4 + 2880.0d0*D81 * DL**3                &
     &            + 399.2d0 * DL**2 + 1465.2d0 * DL                         &
     &            + 714.1d0 * DL1 + DL*DL1 * (684.0d0 + 251.2d0 * DL)         &
     &        + NF * ( -216.62d0 + 406.5d0 * Y + 77.89d0 * Y**2 + 34.76d0 * Y**3&
     &            - 256.0d0*D81 * DL**3  - 3216.0d0*D81 * DL**2 - 172.69d0 * DL &
     &            - 5120.0d0*D81 * DL1 - 65.43d0 * DL*DL1 - 1.136d0 * Y*DL**3 )&
     &        + NF**2 * ( 32.0d0* Y*DL/(1.0d0-Y) * (3.0d0* DL + 10.0d0) + 64.0d0       &
     &            + (48.0d0* DL**2 + 352.0d0* DL + 384.0d0) * (1.0d0-Y) ) * D81     
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the singular piece of both P2_NS+ and P2_NS-. It is exact   
!    up to the truncation of the irrational coefficients.               
!                                                                       
       FUNCTION P2NSB (Y, NF) 
       IMPLICIT REAL*8 (A-Z) 
       INTEGER NF 
!                                                                       
       P2NSB = ( 1174.898d0 - NF * 183.187d0 - NF**2 * 64.0d0/81.0d0 ) / (1.0d0-Y) 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the 'local' piece of P2_NS+. The coefficients of delta(1-x) 
!    have been partly shifted relative to the exact (truncated) values. 
!                                                                       
       FUNCTION P2NSPC (Y, NF) 
       IMPLICIT REAL*8 (A - Z) 
       INTEGER NF 
!                                                                       
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       P2NSPC =       1174.898d0 * DL1 + 1295.624d0 - 0.24d0                  &
     &        - NF * ( 183.187d0 * DL1 + 173.938d0 - 0.011d0 )                &
     &        + NF**2 * ( - 64.0d0/81.0d0 * DL1 + 1.13067d0 )                 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the 'local' piece of P2_NS-. The coefficients of delta(1-x) 
!    have been partly shifted relative to the exact (truncated) values. 
!                                                                       
       FUNCTION P2NSMC (Y, NF) 
       IMPLICIT REAL*8 (A - Z) 
       INTEGER NF 
!                                                                       
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       P2NSMC =       1174.898d0 * DL1 + 1295.624d0 - 0.154d0                 &
     &        - NF * ( 183.187d0 * DL1 + 173.938d0  - 0.005d0 )               &
     &        + NF**2 * ( - 64.0d0/81.0d0 * DL1 + 1.13067d0 )                 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is P2_NSS, the difference of P2_NSV and P2_NS-.                
!                                                                       
       FUNCTION P2NSSA (Y, NF) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER NF 
!                                                                       
       D27 = 1.0d0/27.0d0 
       DL  = LOG (Y) 
       Y1  = 1.0d0- Y 
       DL1 = LOG (Y1) 
!                                                                       
       P2NSSA = Y1* ( 151.49d0 + 44.51d0 * Y - 43.12d0 * Y**2 + 4.820d0 * Y**3 )&
     &          + 40.0d0*D27 * DL**4 - 80.0d0*D27 * DL**3 + 6.892d0 * DL**2     &
     &          + 178.04d0 * DL + DL*DL1 * ( - 173.1d0 + 46.18d0 * DL )       &
     &          + Y1*DL1 * ( - 163.9d0 / Y - 7.208d0 * Y )                  
!                                                                       
       P2NSSA  = NF * P2NSSA 
!                                                                       
       RETURN 
      END FUNCTION
end module xpns2p
