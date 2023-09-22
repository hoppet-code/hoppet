      module xc3ns2p
      contains
*     
* ..File: xc3ns2p.f    F3_NS
*
*
* ..Calculation of the 2-loop x-space MS(bar) coefficient functions 
*    for xF3 via compact parametrizations involving only logarithms.
*    mu_r = mu_f = Q. Expansion parameter: alpha_s/(4 pi).
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
*     E.B. Zijlstra and W.L. van Neerven, Phys. Lett. B297 (1992) 377.
*
* 
* =====================================================================
*
*                                                               __
* ..This is the regular non-singlet piece for the sum F3(nu)+F3(nu), 
*    corresponding to C3NSP-C3NSN in W. van Neerven's program. The 
*    (9+8) numerical coefficients are fitted to his results, using x 
*    values between 10^-6 and 1-10^-6. 
*     
      FUNCTION C3NM2A (Y, NF)
      IMPLICIT REAL*8 (A-Z)
      INTEGER NF
*     
      DL  = LOG (Y)
      OMY = 1.0D0 - Y
      DL1 = LOG (OMY)
*     
      C3NM2A = - 206.1D0 - 576.8D0 * Y - 3.922D0 * DL**3 - 33.31D0 * DL
     $     **2 - 67.60D0 * DL - 15.20D0 * DL1**3 + 94.61D0 * DL1**2 -
     $     409.6D0 * DL1 - 147.9D0 * DL * DL1**2 + NF * ( - 6.337D0 -
     $     14.97D0 * Y + 2.207D0 * DL**2 + 8.683D0 * DL + 0.042D0 * DL1
     $     **3 - 0.808D0 * DL1**2 + 25.00D0 * DL1 + 9.684D0 * DL * DL1 )
     $     
*     
      RETURN
      END function 
*
* ---------------------------------------------------------------------
*
*                                                                 __
* ..This is the regular non-singlet piece for the diff. F3(nu)-F3(nu), 
*    corresponding to C3NSP+C3NSN in WvN's program. For the NF^0 piece
*    7 numerical coefficients are fitted to his results, the ones of
*    ln^3(1-y) and ln^2(1-y) are taken over from C3NM2A. The NF piece
*    is also the same as in C3NM2A.
*
       FUNCTION C3NP2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
*
       C3NP2A = - 242.9D0 - 467.2D0 * Y - 3.049D0 * DL**3 - 30.14D0 * DL
     $      **2 - 79.14D0 * DL - 15.20D0 * DL1**3 + 94.61D0 * DL1**2 -
     $      396.1D0 * DL1 - 92.43D0 * DL * DL1**2 + NF * ( - 6.337D0 -
     $      14.97D0 * Y + 2.207D0 * DL**2 + 8.683D0 * DL + 0.042D0 * DL1
     $      **3 - 0.808D0 * DL1**2  + 25.00D0 * DL1 + 9.684D0 * DL * DL1
     $      )     
*
       RETURN
       END function 
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular NS piece, denoted by SOFT2 in WvN's program. 
*    It is the same for all F2 and F3 cases. The numerical coefficients 
*    are exact, but truncated.
*
       FUNCTION C3NS2B (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
       DM  = 1.D0/OMY
*
       C3NS2B = + 14.2222D0 * DL1**3 - 61.3333D0 * DL1**2 - 31.105D0 *
     $      DL1 + 188.64D0 + NF * ( 1.77778D0 * DL1**2 - 8.5926D0 * DL1
     $      + 6.3489D0 ) 
       C3NS2B = DM * C3NS2B
*
       RETURN
       END function 
*
* ---------------------------------------------------------------------
*
*                                        __
* ..This is the 'local' NS piece for the nu+nu F3, denoted by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (- 0.104 + 0.013 NF) using the lowest moments.
*
       FUNCTION C3NM2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
*
       C3NM2C = + 3.55555D0 * DL1**4 - 20.4444D0 * DL1**3 - 15.5525D0 *
     $      DL1**2 + 188.64D0 * DL1 - 338.531D0 - 0.104D0 + NF *
     $      (0.592593D0 * DL1**3 - 4.2963D0 * DL1**2 + 6.3489D0 * DL1 +
     $      46.844D0 + 0.013D0)
*
       RETURN
       END function 
*
* ---------------------------------------------------------------------
*
*                                        __
* ..This is the 'local' NS piece for the nu-nu F3, also given by COR2 in
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (- 0.152 + 0.013 NF) using the lowest moments.
*
       FUNCTION C3NP2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
*
       C3NP2C = + 3.55555D0 * DL1**4 - 20.4444D0 * DL1**3 - 15.5525D0 *
     $      DL1**2 + 188.64D0 * DL1 - 338.531D0 - 0.152D0 + NF *
     $      (0.592593D0 * DL1**3 - 4.2963D0 * DL1**2 + 6.3489D0 * DL1 +
     $      46.844D0 + 0.013D0)
*
       RETURN
       END function 
*
* =================================================================av==
      end module xc3ns2p
