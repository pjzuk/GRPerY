      SUBROUTINE FORCE_PER(FEXT,CONFL)
      USE SIZE        ! NN
      USE CONFIG      ! CONFL
      USE FORCE_PAR   ! k_PAR, L0_PAR, A_PAR
      IMPLICIT NONE
      REAL*8 FEXT(3*NN),CONFL(3,NN)
      REAL*8 FS(3*NN),FB(3*NN),FLJ(3*NN)    ! STRECHING and BENDING FORCES
      REAL*8 L(1:NN+1),LINV(1:NN+1),T(3,NN+2)
      REAL*8 R,L0,T0(3),FIJ(6),FIJK(9)
      INTEGER I,J,K,M
      REAL*8 B1(3),B2(3)
      REAL*8 CS,S,TH,DF
      REAL*8 a11,a12,a22
      REAL*8 ELONGMAX,LJ_SIGMAI

      T=0.D0
      L   =0.D0
      LINV=0.D0
      DO I=1,NN-1
       T(:,I)=CONFL(:,I)-CONFL(:,I+1)
       CALL PER_SKEW_CORRECTION(T(:,I),L(I))      
        IF (L(I).GT.0.D0) THEN
         T(:,I)=T(:,I)/L(I)
        ENDIF
         LINV(I)=1.D0/L(I)
      ENDDO

      FLJ=0.D0
 
      DO I=1,NN-1
       DO J=I+1,NN
        T0=CONFL(:,I)-CONFL(:,J)
        CALL PER_SKEW_CORRECTION(T0,L0)
        CALL FORCE_LJ(T0,FIJ)
        FLJ(3*(I-1)+1:3*I) = FLJ(3*(I-1)+1:3*I) + FIJ(1:3)
        FLJ(3*(J-1)+1:3*J) = FLJ(3*(J-1)+1:3*J) + FIJ(4:6)
       ENDDO
      ENDDO

CC TOTAL FORCE -------------------------------------------------

      FEXT=FLJ

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE FORCE_LJ(R,F)
      USE FORCE_PAR   ! LJ_SIGMA, LJ_EPS 
      IMPLICIT NONE
      REAL*8 R(3),T0(3),L0,F(6)        
      INTEGER I,J

      T0=R
      L0=SQRT(SUM(T0**2))
      IF(L0.LE.LJ_CUT) THEN
       T0=T0/L0
       F(1:3) =
     *    +4.D0*LJ_EPS/LJ_SIGMA*(
     *      12.D0*(LJ_SIGMA/L0)**13
     *     -6.D0*(LJ_SIGMA/L0)**7
     *     )*T0
       F(4:6)=
     *    -4.D0*LJ_EPS/LJ_SIGMA*(
     *      12.D0*(LJ_SIGMA/L0)**13
     *     -6.D0*(LJ_SIGMA/L0)**7
     *     )*T0
      ELSE
       F=0.D0
      ENDIF

      RETURN
      END



