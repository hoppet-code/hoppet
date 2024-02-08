      MODULE XC2NS2P
      CONTAINS
*     
* ..File: xc2ns2p.f    F2_NS
*
*
* ..Calculation of the 2-loop x-space MS(bar) coefficient functions 
*    for F2 via compact parametrizations involving only logarithms.
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
*     E.B. Zijlstra and W.L. van Neerven, Phys. Lett. B272 (1991) 127.
*
* 
* =====================================================================
*
*
* ..This is the regular non-singlet piece for the electromagnetic F2, 
*    corresponding to C2NSP+C2NSN in W. van Neerven's program. The 
*    (10+8) numerical coefficients are fitted to his results, using x 
*    values between 10^-6 and 1-10^-6. 
*
       FUNCTION C2NN2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
*
       C2NN2A = 
     1          - 69.59D0 - 1008.D0* Y
     2          - 2.835D0 * DL**3 - 17.08D0 * DL**2 + 5.986D0 * DL 
     3          - 17.19D0 * DL1**3 + 71.08D0 * DL1**2 - 660.7D0 * DL1
     4          - 174.8D0 * DL * DL1**2 + 95.09D0 * DL**2 * DL1
     5        + NF * ( - 5.691D0 - 37.91D0 * Y 
     6          + 2.244D0 * DL**2 + 5.770D0 * DL 
     7          - 1.707D0 * DL1**2  + 22.95D0 * DL1
     8          + 3.036D0 * DL**2 * DL1 + 17.97D0 * DL * DL1 )     
*
       RETURN
       END FUNCTION C2NN2A
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular non-singlet piece for the odd-moment (CC) F2, 
*    corresponding to C2NSP-C2NSN in WvN's program. For the NF^0 piece
*    8 numerical coefficients are fitted to his results, the ones of
*    ln^3(1-y) and ln^2(1-y) are taken over from C3NN2A. The NF piece
*    is also the same as in C3NN2A.
*
       FUNCTION C2NC2A (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
*
       C2NC2A = 
     1          - 84.18D0 - 1010.D0* Y
     2          - 3.748D0 * DL**3 - 19.56D0 * DL**2 - 1.235D0 * DL 
     3          - 17.19D0 * DL1**3 + 71.08D0 * DL1**2 - 663.0D0 * DL1
     4          - 192.4D0 * DL * DL1**2 + 80.41D0 * DL**2 * DL1
     5        + NF * ( - 5.691D0 - 37.91D0 * Y 
     6          + 2.244D0 * DL**2 + 5.770D0 * DL 
     7          - 1.707D0 * DL1**2  + 22.95D0 * DL1
     8          + 3.036D0 * DL**2 * DL1 + 17.97D0 * DL * DL1 )     
*
       RETURN
       END FUNCTION C2NC2A
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular NS piece, denoted by SOFT2 in WvN's program. 
*    It is the same for all F2 and F3 cases. The numerical coefficients 
*    are exact, but truncated.
*
       FUNCTION C2NS2B (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
       DM  = 1.D0/(OMY)
*
       C2NS2B = + 14.2222D0 * DL1**3 - 61.3333D0 * DL1**2 - 31.105D0 *
     $      DL1 + 188.64D0 + NF * ( 1.77778D0 * DL1**2 - 8.5926D0 * DL1
     $      + 6.3489D0 ) 
       C2NS2B = DM * C2NS2B
*
       RETURN
       END FUNCTION C2NS2B
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the e.m. F2, denoted by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one (from the delta-function) is 
*    slightly adjusted (+ 0.485 - 0.0035 NF) using the lowest moments.
*
       FUNCTION C2NN2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       OMY = 1.0D0 - Y
       DL1 = LOG (OMY)
*
       C2NN2C = + 3.55555D0 * DL1**4 - 20.4444D0 * DL1**3 - 15.5525D0 *
     $      DL1**2 + 188.64D0 * DL1 - 338.531D0 + 0.485D0 + NF *
     $      (0.592593D0 * DL1**3 - 4.2963D0 * DL1**2 + 6.3489D0 * DL1 +
     $      46.844D0 - 0.0035D0)
*
       RETURN
       END FUNCTION C2NN2C
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' NS piece for the CC F2, also given by COR2 in 
*    WvN's program. The numerical coefficients of the logs are exact,
*    but truncated, the constant one is adjusted (- 0.2652 - 0.0035 NF) 
*    using the lowest moments.
*
       FUNCTION C2NC2C (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       OMY = 1.0D0 - Y 
       DL1 = LOG (OMY)
*
       C2NC2C = + 3.55555D0 * DL1**4 - 20.4444D0 * DL1**3 - 15.5525D0 *
     $      DL1**2 + 188.64D0 * DL1 - 338.531D0 + 0.537D0 + NF *
     $      (0.592593D0 * DL1**3 - 4.2963D0 * DL1**2 + 6.3489D0 * DL1 +
     $      46.844D0 - 0.0035D0)
*
       RETURN
       END FUNCTION C2NC2C

      end module xc2ns2p

*     
* =================================================================av==
