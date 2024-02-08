      MODULE XCLNS2P
      CONTAINS
*
* ..File: xclns2p.f    FL_NS
*
*
* ..Calculation of the 2-loop x-space MS(bar) coefficient functions 
*    for FL via compact parametrizations involving only logarithms.
*    Non-singlet, mu_r = mu_f = Q. Expansion parameter: alpha_s/(4 pi).
*
*  ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to 
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of the coefficient functions, as well as of 
*    the convolution results, amounts to a few thousandth.
*
*  ..Reference: W.L. van Neerven and A. Vogt, 
*               hep-ph/9907472 = Nucl. Phys. B568 (2000) 263
*  ..The user should also cite the original calculation,
*      J. Sanchez Guillen et al, Nucl. Phys. B353 (1991) 337.
*    
* 
* =====================================================================
*
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to CLNSP+C2NSM in W. van Neerven's program. The 
*    8 numerical coefficients are fitted to his results, using x values 
*    between 10^-6 and 1-10^-6. 
*
       FUNCTION CLNN2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       DOUBLE PRECISION D16O27
       PARAMETER (D16O27 = 16.0D0/27.0D0)
       INTEGER NF
*
       DL  = LOG (Y)
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
*
       CLNN2A = - 40.41D0 + 97.48D0 * Y + (26.56D0 * Y - 0.031D0) * DL
     $      **2 - 14.85D0 * DL + 13.62D0 * DL1**2 - 55.79D0 * DL1 -
     $      150.5D0 * DL * DL1 + NF * D16O27 * ( 6.D0* Y*DL1 - 12.D0* Y
     $      *DL - 25.D0* Y + 6.D0)
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular non-singlet piece for the odd-moment (CC) FL, 
*    corresponding to CLNSP-CLNSN in WvN's program. The 8 numerical 
*    coefficients the NF^0 piece are fitted to his results. The NF part
*    is the same as in CLNN2A.
*
       FUNCTION CLNC2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       DOUBLE PRECISION D16O27
       PARAMETER (D16O27 = 16.0D0/27.0D0)
       INTEGER NF
*
       DL  = LOG (Y)
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
*
       CLNC2A = - 52.27D0 + 100.8D0 * Y + (23.29D0 * Y - 0.043D0) * DL
     $      **2 - 22.21D0 * DL + 13.30D0 * DL1**2 - 59.12D0 * DL1 -
     $      141.7D0 * DL * DL1 + NF * D16O27 * ( 6.D0* Y*DL1 -
     $      12.D0* Y*DL - 25.D0* Y + 6.D0)
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. FL, with no counterpart 
*    in WvN's program, as it does not exist in the exact expressions.
*    The value is fixed from the lowest integer moments.
*
       FUNCTION CLNN2C (Y)
       IMPLICIT REAL*8 (A-Z)
*
       CLNN2C = -0.164D0
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the CC FL, with no counterpart 
*    in WvN's program, as it does not exist in the exact expressions.  
*    The value is fixed from the lowest integer moments.
*
       FUNCTION CLNC2C (Y)
       IMPLICIT REAL*8 (A-Z)
*
       CLNC2C = -0.150D0
*
       RETURN
       END FUNCTION
*
* =================================================================av==

      END MODULE XCLNS2P
