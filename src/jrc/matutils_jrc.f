************************************************************************
*                                                        	       *
*      Program: matutils.f                                             *
*							  	       *
*      A collection of Fortran utilities to deal with matrices:        *
*        - multiplication (mmult_c, mmult_r). 			       *
*        - equating (mequal_c,mequal_r)                                *
*                                                                      *
*      The ending _# refers to the type of matrix entries:             *
*         C=complex, R=real                                            *     
*		 	       					       *
*      Last Modified: 18.01.07, AG			 	       *
*	        						       *
************************************************************************
*
      SUBROUTINE MMULT_C(A,ROWSA,COLSA,B,ROWSB,COLSB,C)
*     ----------------------------------------------------------
*
      IMPLICIT none
*
      INTEGER ROWSA,COLSA,ROWSB,COLSB
*
      COMPLEX*16 A(ROWSA,COLSA), B(ROWSB,COLSB), C(ROWSA,COLSB)
*
      INTEGER I,J,K
*
*     Check that the matrices we would like to multiply have the 
*     correct dimensions
*
      IF (COLSA .NE. ROWSB) THEN
          WRITE(*,*) 'Coglione... dimensioni delle matrici sbagliate!'
          RETURN
      ELSE
*
*     Initialize the output matrix to zero
*
        DO I = 1,ROWSA
          DO J = 1,COLSB
	    C(I,J) = (0.0d0,0.0d0)
          ENDDO
        ENDDO

*
*     Perform the multiplication according to: 
*                          C[i,j] = Sum_k(A[i,k]*B[k,i])
*
        DO I = 1,ROWSA
          DO J = 1,COLSB
	    DO K = 1,COLSA
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)
	    ENDDO
	  ENDDO
        ENDDO
*
      ENDIF
      RETURN
      END
*-----------------------------------------------------------------------
*
      SUBROUTINE MMULT_R(A,ROWSA,COLSA,B,ROWSB,COLSB,C)
*     ----------------------------------------------------------
*
      IMPLICIT none
*
      INTEGER ROWSA,COLSA,ROWSB,COLSB
*
      REAL*8 A(ROWSA,COLSA), B(ROWSB,COLSB), C(ROWSA,COLSB)
*
      INTEGER I,J,K
*
*     Check that the matrices we would like to multiply have the 
*     correct dimensions
*
      IF (COLSA .NE. ROWSB) THEN
          WRITE(*,*) 'Coglione... dimensioni delle matrici sbagliate!'
          RETURN
      ELSE
*
*     Initialize the output matrix to zero
*
        DO I = 1,ROWSA
          DO J = 1,COLSB
	    C(I,J) = 0.0d0
          ENDDO
        ENDDO

*
*     Perform the multiplication according to: 
*                          C[i,j] = Sum_k(A[i,k]*B[k,i])
*
        DO I = 1,ROWSA
          DO J = 1,COLSB
	    DO K = 1,COLSA
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)
	    ENDDO
	  ENDDO
        ENDDO
*
      ENDIF
      RETURN
      END
*----------------------------------------------------------------------
*
      SUBROUTINE MEQUAL_C(A,B,I,J)
*
      IMPLICIT none
*
      INTEGER I,J
      INTEGER K,L
      DOUBLE COMPLEX A(I,J), B(I,J) 
*
      DO K=1,I
         DO L=1,J
            A(K,L) = B(K,L)
         ENDDO
      ENDDO
      RETURN
      END
*
*----------------------------------------------------------------------
*
      SUBROUTINE MEQUAL_R(A,B,I,J)
*
      IMPLICIT none
*
      INTEGER I,J
      INTEGER K,L
      DOUBLE PRECISION A(I,J), B(I,J) 
*
      DO K=1,I
         DO L=1,J
            A(K,L) = B(K,L)
         ENDDO
      ENDDO
      RETURN
      END
*-----------------------------------------------------------------------
