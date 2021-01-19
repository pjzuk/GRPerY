C Procedures here are most computationally costly.
C For better performance it is advised to replace them with parallel
C implementations.


C  lapack based matrix inversion

      SUBROUTINE MATREV(A,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 A(NN,NN)
      INTEGER INFO,I,J
      CHARACTER*1 UPLO
      PARAMETER (UPLO = 'U')

      CALL DPOTRF(UPLO,NN,A,NN,INFO)
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
      IF(INFO.GT.0) THEN
      WRITE(*,*) 'k=',INFO
      WRITE(*,*) 'The leading minor of order k is not'
      WRITE(*,*) 'positive definite, and the factorization could not be'
      WRITE(*,*) 'completed.'
      ENDIF
      CALL DPOTRI(UPLO,NN,A,NN,INFO)

      DO I=2,NN
       DO J=1,I-1
        A(I,J)=A(J,I)
       ENDDO
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

C lapack based Cholesky decomposition

      SUBROUTINE CHOLESKY(B,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 B(NN,NN)
      CHARACTER*1 UPLO
      PARAMETER (UPLO = 'L')
      INTEGER I,J,INFO

      CALL DPOTRF(UPLO,NN,B,NN,INFO)
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
      IF(INFO.GT.0) THEN
      WRITE(*,*) 'k=',INFO
      WRITE(*,*) 'The leading minor of order k is not'
      WRITE(*,*) 'positive definite, and the factorization could not be'
      WRITE(*,*) 'completed.'
      STOP
      ENDIF

C     RETURN     !!!  The strictly upper triangular part of A is not referenced.
      DO I=1,NN-1
      DO J=I+1,NN
       B(I,J)=0.D0
      ENDDO
      ENDDO

      RETURN
      END

 
