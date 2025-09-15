      MODULE XCDIFF3P
      USE XC2NS3P
      CONTAINS
*     
* ..File: xcdiff3p.f    Fi_NS,  i = 2,3,L  (even-N - odd-N)
*
*
* ..Approximations for the differences between the even-N and odd-N 
*    based three-loop coefficient functions for the structure functions
*    F_L, F_2 and F_3 in charged-current deep-inelastic scattering.
*    The expansion parameter is normalized as  a_s = alpha_s/(4 pi).
*
* ..The two sets providing the error estimate are called via  IMOD = 1
*    and  IMOD = 2.  Any other value of IMOD invokes their average.
*    The results are based on the first five even/odd-integer moments
*    as derived in S. Moch and M. Rogal, arXiv:0704:1740 [hep-ph] and
*    complete results mentined below.
*
* ..Reference: M. Rogal, S. Moch and A. Vogt, arXiv:0708.3731 [hep-ph]
*
* ..The functions given here should be used together with the complete 
*    results for the even-N based F_L, F_2 and the odd-N based F_3 
*    respectively presented by S. Moch, J. Vermaseren and A. Vogt in 
*    hep-ph/0411112, hep-ph/0504242 and hep-ph/0608307.
*
* ===================================================================== 
*
* ..F2
*
       FUNCTION C2Q3DF (Y, DL, NF, IMOD)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, NF
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
*
       C2Q30A = (54.478 * DL1**2 + 304.60 * DL1 + 691.68 * Y) * Y1 
     1          + 179.14 * DL - 0.1826 * DL**3
       C2Q30B = - (13.378 * DL1**2 + 97.60 * DL1 + 118.12 * Y) * Y1 
     1          - 91.196 * DL**2 - 0.4644 * DL**5
*
       C2Q31A = (20.822 * Y**2 - 282.10 * (1.+Y/2.)) * Y1
     1          - 285.58 * Y*DL - 112.30 * DL + 3.587 * DL**3
       C2Q31B = (4.522 * DL1 + 447.88 * (1.+Y/2.)) * Y1
     1          + 514.02 * Y*DL + 147.05 * DL + 7.386 * DL**2
*
       IF (IMOD .EQ. 1) THEN
         C2Q3DF = C2Q30A + NF * C2Q31A
       ELSE IF (IMOD .EQ. 2) THEN
         C2Q3DF = C2Q30B + NF * C2Q31B
       ELSE
         C2Q3DF = 0.5 * (C2Q30A + C2Q30B + NF * (C2Q31A + C2Q31B))
       END IF
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
* ..F3
*
       FUNCTION C3Q3DF (Y, DL, NF, IMOD)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, NF
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
*
       C3Q30A = - (46.72 * DL1**2 + 267.26 * DL1 + 719.49 * Y) * Y1
     1          - 171.98 * DL + 9.470 * DL**3
       C3Q30B = (3.216 * DL1**2 + 44.50 * DL1 - 34.588) * Y1
     1          + 98.719 * DL**2 + 2.6208 * DL**5
*
       C3Q31A = (0.8489 * DL1 + 67.928 * (1.+0.5*Y)) * Y1
     1          + 97.922 * Y*DL - 17.070 * DL**2 - 3.132 * DL**3
       C3Q31B = - (0.186 * DL1 + 61.102 * (1.+Y)) * Y1
     1          - 122.51 * Y*DL + 10.914 * DL**2 + 2.748 * DL**3
*
       IF (IMOD .EQ. 1) THEN
         C3Q3DF = C3Q30A + NF * C3Q31A
       ELSE IF (IMOD .EQ. 2) THEN
         C3Q3DF = C3Q30B + NF * C3Q31B
       ELSE
         C3Q3DF = 0.5 * (C3Q30A + C3Q30B + NF * (C3Q31A + C3Q31B))
       END IF
*
       RETURN
       END FUNCTION
*
* ---------------------------------------------------------------------
*
* ..FL
*
       FUNCTION CLQ3DF (Y, DL, NF, IMOD)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, NF
*
       Y1  = Y1VAL(Y, DL)
       DL1 = DL1VAL(Y, DL)
*
       CLQ30A = - (495.49 * Y**2 + 906.86) * Y1**2
     1          - 983.23 * Y*Y1*DL + 53.706 * DL**2 + 5.3059 * DL**3
       CLQ30B = (78.306 * DL1 + 6.3838 * Y) *Y1**2
     1          + 20.809 * Y*Y1*DL - 114.47 * DL**2 - 22.222 * DL**3
*
       CLQ31A = (29.95 * Y**3 - 59.087 * Y**2 + 379.91) * Y1**2
     1          - 273.042 * Y*DL**2 + 71.482 * Y1*DL
       CLQ31B = (12.532 * DL1 + 141.99 * Y**2 - 250.62 * Y)* Y1**2
     1          - (153.586 * Y - 0.6569) * Y1*DL
*
       IF (IMOD .EQ. 1) THEN
         CLQ3DF = CLQ30A + NF * CLQ31A
       ELSE IF (IMOD .EQ. 2) THEN
         CLQ3DF = CLQ30B + NF * CLQ31B
       ELSE
         CLQ3DF = 0.5 * (CLQ30A + CLQ30B + NF * (CLQ31A + CLQ31B))
       END IF
*
       RETURN
       END FUNCTION
*
* =================================================================av==
      END MODULE XCDIFF3P
