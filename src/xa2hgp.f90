!  $Id: xa2hgp.f90,v 1.1 2002/03/06 17:12:12 gsalam Exp $ 
!  Automatically generated from f77 file, with addition of "d0"
!  From avogt@nikhef.nl Thu Feb  7 15:34:35 2002                        
!  Date: Wed, 6 Feb 2002 18:23:05 +0100 (MET)                           
!  From: Andreas Vogt <avogt@nikhef.nl>                                 
!  To: G. Salam <Gavin.Salam@cern.ch>                                   
!  Cc: A. Vogt <avogt@nikhef.nl>                                        
!  Subject: the subroutine                                              
!                                                                       
!                                                                       
! ..File: xa2hgp.f                                                      
!                                                                       
!                                                                       
! ..Calculation of the alpha_s^2 heavy-quark singlet operator matrix    
!    element (OME) a^2_Hg(x) in the MS(bar) scheme for mu_f = m_H via   
!    a compact parametrization involving only logarithms.               
!    The coupling constant is normalized as  a_s = alpha_s/(4*pi).      
!                                                                       
! ..This quantity, presented in Appendix B of M. Buza, Y. Matiounine,   
!    J. Smith and W.L. van Neerven, Eur. Phys. J. C1 (1998) 301 (BSMN), 
!    is required for the N_f matching of the NNLO parton densities.     
!                                                                       
!  ..The distributions (in the mathematical sense) are given as in eq.  
!    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.  
!    The name-endings A, B, and C of the functions below correspond to  
!    the kernel superscripts [2], [3], and [1] in that equation.        
!                                                                       
!  ..The relative accuracy of the OME, as well as of its convolutions   
!    with gluon distributions, amounts to a few thousandth.             
!                                                                       
!  ..Reference: not yet known                                           
!               hep-ph/yymmnnn                                          
!  ..The user should also cite the original calculation by BSMN.        
!                                                                       
!                                                                       
! ===================================================================== 
!                                                                       
!                                                                       
! ..This is the regular piece.                                          
!                                                                       
       FUNCTION A2HGA (Y) 
       IMPLICIT REAL*8 (A-Z) 
!                                                                       
       DL  = LOG (Y) 
       DL1 = LOG (1.-Y) 
!                                                                       
       A2HGA = - 24.89d0 / Y - 187.8d0 + 249.6d0 * Y - 146.8d0 * DL**2 * DL1    &
     &         - 1.556d0 * DL**3  - 3.292d0 * DL**2  - 93.68d0 * DL           &
     &         - 1.111d0 * DL1**3 - 0.400d0 * DL1**2 - 2.770d0 * DL1          
!                                                                       
       RETURN 
      END                                           
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
!                                                                       
! ..This is the 'local' piece, which has no counterpart in W. van       
!    Neerven's program, as it does not exist in the exact expressions.  
!                                                                       
       FUNCTION A2HGC (Y) 
       IMPLICIT REAL*8 (A-Z) 
!                                                                       
       A2HGC = - 0.006d0 
!                                                                       
       RETURN 
      END                                           
