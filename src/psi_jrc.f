* =====================================================================
*
*     Taken from Vogt Pegasus code (hep-ph/0408244)
*
* ..The complex psi function,  PZI(Z),  and its m'th derivatives, 
*    DPZI(Z,M),  calculated from the asymtotic expansions. The 
*    functional equations are used for  |Im(Z)| < 10  to shift the 
*    argument to  Re(Z) >= 10  before applying the expansions.
*
* =====================================================================
*
*
       FUNCTION PSI (Z)
*
       IMPLICIT DOUBLE COMPLEX (A-Z)
       SUB = 0D0
       ZZ = Z
*
* ..Shift of the argument using the functional equation
*
       IF (DABS(DIMAG(ZZ)) .LT. 10.D0) THEN
*
  1      CONTINUE
         IF (DBLE(ZZ) .LT. 10.D0) THEN
           SUB = SUB - 1d0/ ZZ
           ZZ = ZZ + 1d0
           GOTO 1
         END IF
*
       END IF
*
* ..Use of the asymtotic expansion (at the shifted argument)
*     Abramowitz-Stegun eq. 6.3.18 plus one term
*
       RZ = 1d0/ ZZ
       DZ = RZ * RZ
       PSI = SUB + LOG(ZZ) - 0.5d0 * RZ - DZ/5040D0 * ( 420d0+ DZ *
     1       ( - 42d0 + DZ * (20d0 - 21d0 * DZ) ) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
       FUNCTION DPSI (Z,M)
*
       IMPLICIT DOUBLE COMPLEX (A-Z)
       INTEGER M, K1, K2
       SUB = 0D0
       ZZ = Z
*
* ..Shift of the argument using the functional equations
*
       IF (DABS(DIMAG(ZZ)) .LT. 10.D0) THEN
*
  1      CONTINUE
         SUBM = -1d0/ZZ
         DO 10 K1 = 1, M
           SUBM = - SUBM * K1 / ZZ
 10      CONTINUE
*
         IF (DBLE(ZZ) .LT. 10.D0) THEN
           SUB = SUB + SUBM
           ZZ = ZZ + 1d0
           GOTO 1
         END IF
*
       END IF
*
* ..Expansion (Bernoulli) coefficients for the first derivative
*
       A1 =  1D0
       A2 =  1d0/2D0
       A3 =  1d0/6D0
       A4 = -1d0/30D0
       A5 =  1d0/42D0
       A6 = -1d0/30D0
       A7 =  5d0/66D0
*
* ..Expansion coefficients for the higher derivatives
*
       IF (M .EQ. 1) GO TO 2
       DO 11 K2 = 2, M
         A1 = A1 * (dble(K2)-1d0)
         A2 = A2 *  dble(K2)
         A3 = A3 * (dble(K2)+1d0)
         A4 = A4 * (dble(K2)+3d0)
         A5 = A5 * (dble(K2)+5d0)
         A6 = A6 * (dble(K2)+7d0)
         A7 = A7 * (dble(K2)+9d0)
  11   CONTINUE
  2    CONTINUE 
*
* ..Use of the asymtotic expansion (at the shifted argument)
*
       RZ = 1d0/ ZZ
       DZ = RZ * RZ
       DPSI = SUB + (-1)**(M+1d0) * RZ**M * ( A1 + RZ * (A2 + RZ * 
     1        (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) )
*
       RETURN
       END
*
* =================================================================av==
