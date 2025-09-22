!  $Id: xpij2p.f90,v 1.1 2004/06/01 09:30:23 salam Exp $ 
!  Automatically generated from f77 file, with addition of "d0"
!  and the placement inside a module.
module xpij2p
character(len=*), parameter :: name_xpij2 = "xpij2p"
contains
!                                                                       
! ..File: xpij2p.f                                                      
!                                                                       
!                           __                                          
! ..The parametrized 3-loop MS singlet splitting functions P^(2) for    
!    the evolution of unpol. singlet parton densities at mu_r = mu_f.   
!    The expansion parameter is alpha_s/(4 pi).                         
!                                                                       
! ..The distributions (in the mathematical sense) are given as in eq.   
!   (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.0d0   
!   The name-endings A, B, and C of the functions below correspond to   
!   the kernel superscripts [2], [3], and [1] in that equation.         
!                                                                       
! ..The relative accuracy of these parametrisations, as well as of      
!    the convolution results, is better than one part in thousand.      
                                                                        
! ..The coefficients of 1/(1-x)_+, (ln x)/x and 1/x are exact (up       
!    to a truncation of irrational coefficients).  Furthermore all      
!    coefficients written as fractions (e.g., 160.0d0/27.0d0) are exact.    
!    The other terms at x < 1 have fitted to the exact results for x    
!    between 10^-6 and 1 - 10^-6.0d0  The coefficient of delta(1-x) of     
!    P_gg^(2) have been slightly adjusted using the second moments.     
!                                                                       
! ..References: S. Moch, J. Vermaseren and A. Vogt,                     
!               hep-ph/0403192 (to appear in Nucl. Phys. B)             
!               A. Vogt, S. Moch and J. Vermaseren,                     
!               hep-ph/0404111 (submitted to Nucl. Phys. B)             
!                                                                       
! ===================================================================== 
!                                                                       
!                                                                       
! ..The (regular) pure-singlet splitting functions P_ps^(2).            
!    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+. 
!    A parametrization of the latter is provided in the file  xpns2p.f. 
                                                                        
       FUNCTION P2PSA (Y, NF) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER NF 
!                                                                       
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       P2PS1 = - 3584.0d0/(27.0d0*Y) * DL - 506.0d0/ Y + 160.0d0/27.0d0 * DL**4   &
     &         - 400.0d0/9.0d0 * DL**3 + 131.4d0 * DL**2 - 661.6d0 * DL         &
     &         - 5.926d0  * DL1**3 - 9.751d0 * DL1**2 - 72.11d0 * DL1         &
     &         + 177.4d0 + 392.9d0 * Y - 101.4d0 * Y**2 - 57.04d0 * DL*DL1      
       P2PS2 =   256.0d0/(81.0d0*Y) + 32.0d0/27.0d0 * DL**3 + 17.89d0 * DL**2       &
     &         + 61.75d0 * DL + 1.778d0 * DL1**2 + 5.944d0 * DL1 + 100.1d0      &
     &         - 125.2d0 * Y + 49.26d0 * Y**2 - 12.59d0 * Y**3                &
     &         - 1.889d0 * DL*DL1                                         
!                                                                       
       P2PSA = (1.0d0-Y) * NF * ( P2PS1 + NF * P2PS2 ) 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..The gluon->quark splitting functions P_qg^(2).                      
!                                                                       
       FUNCTION P2QGA (Y, NF) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER  NF 
!                                                                       
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       P2QG1 = - 896.0d0/(3.0d0*Y) * DL - 1268.3d0 / Y + 536.0d0/27.0d0 * DL**4   &
     &         - 44.0d0/3.0d0 * DL**3 + 881.5d0 * DL**2 + 424.9d0 * DL          &
     &         + 100.0d0/27.0d0 * DL1**4 - 70.0d0/9.0d0 * DL1**3                &
     &         - 120.5d0 * DL1**2 + 104.42d0 * DL1                          &
     &         + 2522.0d0 - 3316.0d0* Y + 2126.0d0* Y**2                         &
     &         + DL*DL1 * (1823.0d0 - 25.22d0 * DL) - 252.5d0 * Y*DL**3        
       P2QG2 =   1112.0d0/(243.0d0*Y) - 16.0d0/9.0d0 * DL**4                    &
     &         - 376.0d0/27.0d0 * DL**3 - 90.8d0 * DL**2 - 254.0d0 * DL         &
     &         + 20.0d0/27.0d0 * DL1**3 + 200.0d0/27.0d0 * DL1**2 - 5.496d0 * DL1 &
     &         - 252.0d0  + 158.0d0 * Y + 145.4d0 * Y**2 - 139.28d0 * Y**3      &
     &         - DL*DL1 * ( 53.09d0  + 80.616d0 * DL) - 98.07d0 * Y*DL**2     &
     &         + 11.70d0 * Y*DL**3                                        
!                                                                       
       P2QGA = NF * ( P2QG1 + NF * P2QG2 ) 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..The quark->gluon splitting functions P_gq^(2).  P2GQ2 is exact.     
!                                                                       
       FUNCTION P2GQA (Y, NF) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER NF 
!                                                                       
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       P2GQ0 =   1189.3d0 * DL/Y + 6163.1d0 / Y - 4288.0d0/81.0d0 * DL**4       &
     &         + 1568.0d0/9.0d0 * DL**3 - 1794.0d0 * DL**2 + 4033.0d0 * DL        &
     &         + 400.0d0/81.0d0 * DL1**4 + 2200.0d0/27.0d0 * DL1**3             &
     &         + 606.3d0 * DL1**2 + 2193.0d0* DL1                            &
     &         - 4307.0d0 + 489.3d0 * Y + 1452.0d0* Y**2 + 146.0d0 * Y**3         &
     &         - 447.3d0 * DL**2*DL1 - 972.9d0 * Y*DL**2                    
       P2GQ1 =   71.082d0 * DL/Y  - 46.41d0 / Y + 128.0d0/27.0d0 * DL**4        &
     &         + 704/81.0d0 * DL**3 + 20.39d0  * DL**2 + 174.8d0 * DL        &
     &         - 400.0d0/81.0d0 * DL1**3 - 68.069d0 * DL1**2 - 296.7d0 * DL1    &
     &         - 183.8d0 + 33.35d0 * Y - 277.9d0 * Y**2 + 108.6d0 * Y*DL**2     &
     &         - 49.68d0 * DL*DL1                                         
       P2GQ2 = ( 64.0d0 * ( - 1.0d0/Y + 1.0d0 + 2.0d0* Y)                           &
     &         + 320.0d0* DL1 * ( 1.0d0/Y - 1.0d0 + 0.8d0 * Y)                     &
     &         + 96.0d0* DL1**2 * ( 1.0d0/Y - 1.0d0 + 0.5d0 * Y) ) / 27.0d0         
!                                                                       
       P2GQA = ( P2GQ0 + NF * (P2GQ1 + NF * P2GQ2) ) 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..The regular piece of the gluon-gluon splitting function P_gg^(2).   
!                                                                       
       FUNCTION P2GGA (Y, NF) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER NF 
!                                                                       
       DL  = LOG (Y) 
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       P2GGA0 = 2675.8d0 * DL/Y + 14214.0d0/ Y - 144.0d0 * DL**4 + 72.0d0 * DL**3  &
     &          - 7471.0d0 * DL**2 + 274.4d0 * DL + 3589.0d0 * DL1 - 20852.0d0     &
     &          + 3968.0d0* Y - 3363.0d0 * Y**2 + 4848.0d0 * Y**3                &
     &          + DL*DL1 * ( 7305.0d0 + 8757.0d0 * DL )                       
       P2GGA1 = 157.27d0 * DL/Y + 182.96d0 / Y + 512.0d0/27.0d0 * DL**4         &
     &          + 832.0d0/9.0d0 * DL**3 + 491.3d0 * DL**2 + 1541.0d0 * DL        &
     &          - 320.0d0 * DL1 - 350.2d0 + 755.7d0 * Y - 713.8d0 * Y**2        &
     &          + 559.3d0 * Y**3 + DL*DL1 * ( 26.15d0 - 808.7d0 * DL )        
       P2GGA2 = - 680.0d0/(243.0d0 * Y) - 32.0d0/27.0d0 * DL**3 + 9.680d0 * DL**2 &
     &          - 3.422d0 * DL - 13.878d0 + 153.4d0 * Y - 187.7d0 * Y**2        &
     &          + 52.75d0 * Y**3 - DL*DL1 * (115.6d0 - 85.25d0* Y + 63.23d0* DL)
!                                                                       
       P2GGA = P2GGA0 + NF * ( P2GGA1 + NF * P2GGA2 ) 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..The singular piece of the gluon-gluon splitting function P_gg^(2).  
!                                                                       
       FUNCTION P2GGB (Y, NF) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER NF 
!                                                                       
       P2GGB = ( 2643.521d0 - NF * 412.172d0 - NF**2 * 16.0d0/9.0d0 ) / ( 1.0d0-Y) 
!                                                                       
       RETURN 
      END FUNCTION
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..The 'local' piece of the gluon-gluon splitting function P_gg^(2).   
!                                                                       
       FUNCTION P2GGC (Y, NF) 
!                                                                       
       IMPLICIT REAL*8 (A-Z) 
       INTEGER NF 
!                                                                       
       DL1 = LOG (1.0d0-Y) 
!                                                                       
       P2GGC =       2643.521d0 * DL1 + 4425.448d0 + 0.446d0                  &
     &       - NF * ( 412.172d0 * DL1 +  528.720d0 + 0.003d0 )                &
     &       + NF**2 * ( - 16.0d0/9.0d0 * DL1 + 6.4630d0)                     
!                                                                       
       RETURN 
      END FUNCTION
end module xpij2p
