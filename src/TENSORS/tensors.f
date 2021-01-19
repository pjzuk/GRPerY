C     4.12.2014 procedury do tensorow


      SUBROUTINE INIT_U()
      USE TENSORS
      IMPLICIT NONE

      U=0.D0

      U(1,1)=1.D0
      U(2,2)=1.D0
      U(3,3)=1.D0

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE INIT_EPS()
      USE TENSORS
      IMPLICIT NONE

      EPS=0.D0

      EPS(1,2,3) = 1.D0
      EPS(1,3,2) = -1.D0
      EPS(2,1,3) = -1.D0
      EPS(2,3,1) = 1.D0
      EPS(3,1,2) = 1.D0
      EPS(3,2,1) = -1.D0

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE INIT_Y2()
      USE TENSORS
      IMPLICIT NONE

      Y2=0.D0

      Y2(1,1,2)=1.D0/SQRT(2.D0)
      Y2(1,2,1)=1.D0/SQRT(2.D0)

      Y2(2,2,3)=-1.D0/SQRT(2.D0)
      Y2(2,3,2)=-1.D0/SQRT(2.D0)

      Y2(3,1,1)=-1.D0/SQRT(6.D0)
      Y2(3,2,2)=-1.D0/SQRT(6.D0)
      Y2(3,3,3)=SQRT(2.D0/3.D0)

      Y2(4,1,3)=-1.D0/SQRT(2.D0)
      Y2(4,3,1)=-1.D0/SQRT(2.D0)

      Y2(5,1,1)=1.D0/SQRT(2.D0)
      Y2(5,2,2)=-1.D0/SQRT(2.D0)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE INIT_II()
      USE TENSORS
      IMPLICIT NONE

      II=0.D0

      II(1,1,1,1)=0.6666666666666666
      II(1,1,2,2)=-0.3333333333333333
      II(1,1,3,3)=-0.3333333333333333
      II(1,2,1,2)=0.5
      II(1,2,2,1)=0.5
      II(1,3,1,3)=0.5
      II(1,3,3,1)=0.5
      II(2,1,1,2)=0.5
      II(2,1,2,1)=0.5
      II(2,2,1,1)=-0.3333333333333333
      II(2,2,2,2)=0.6666666666666666
      II(2,2,3,3)=-0.3333333333333333
      II(2,3,2,3)=0.5
      II(2,3,3,2)=0.5
      II(3,1,1,3)=0.5
      II(3,1,3,1)=0.5
      II(3,2,2,3)=0.5
      II(3,2,3,2)=0.5
      II(3,3,1,1)=-0.3333333333333333
      II(3,3,2,2)=-0.3333333333333333
      II(3,3,3,3)=0.6666666666666666

      RETURN
      END

     
C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE CALC_RR(RR,RW)
      IMPLICIT NONE
      REAL*8 RR(3,3),RW(3)
      INTEGER I,J

      DO I=1,3
       DO J=1,3
        RR(I,J)=RW(I)*RW(J)
       ENDDO
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CALC_EPSR(EPSR,RW)
      IMPLICIT NONE
      REAL*8 EPSR(3,3),RW(3)
      INTEGER I,J

      EPSR=0.D0

      EPSR(1,2)= RW(3)
      EPSR(2,3)= RW(1)
      EPSR(3,1)= RW(2)
      
      EPSR(2,1)=-RW(3)
      EPSR(3,2)=-RW(1)
      EPSR(1,3)=-RW(2)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CALC_UR(t3UR,RW)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 t3UR(3,3,3),RW(3)
      INTEGER I

      t3UR=0.D0

      DO I=1,3
       t3UR(I,I,1:3)= RW
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CALC_RRR(t3RRR,RW)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 t3RRR(3,3,3),RW(3)
      INTEGER I,J,K

      t3RRR=0.D0

      DO I=1,3
       DO J=1,3
        DO K=1,3
         t3RRR(I,J,K)= RW(I)*RW(J)*RW(K)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CALC_EPSRR(t3EPSRR,RW)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 t3EPSRR(3,3,3),EPSR(3,3),RW(3)
      INTEGER I,J,K

      t3EPSRR=0.D0

      CALL CALC_EPSR(EPSR,RW)

      DO I=1,3
       DO J=1,3
        DO K=1,3
         t3EPSRR(I,J,K) = EPSR(I,J)*RW(K)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CALC_t4D0(t4D0,RR)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 t4D0(3,3,3,3),RR(3,3),TMP(3,3)
      INTEGER I,J

      TMP = RR-U/3.D0

      CALL t2at2b_t4(t4D0,TMP,TMP)
      t4D0 = 3.D0/2.D0*t4D0

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CALC_t4D1(t4D1,RR)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 t4D1(3,3,3,3),RR(3,3),TMP(3,3,3,3)
      INTEGER I,J

      t4D1=0.D0

C     DELTAjlRiRk
      CALL t2at2b_t4(TMP,U,RR)
      CALL t4trans(TMP,1,4)
      t4D1=t4D1+TMP
C     DELTAikRjRl
      CALL t2at2b_t4(TMP,U,RR)
      CALL t4trans(TMP,2,3)
      t4D1=t4D1+TMP
C     DELTAjkRiRl
      CALL t2at2b_t4(TMP,U,RR)
      CALL t4trans(TMP,1,3)
      t4D1=t4D1+TMP
C     DELTAilRjRk
      CALL t2at2b_t4(TMP,U,RR)
      CALL t4trans(TMP,2,4)
      t4D1=t4D1+TMP
C     RiRjRkRl
      CALL t2at2b_t4(TMP,RR,RR)
      t4D1=t4D1 -4.D0*TMP

      t4D1=0.5D0*t4D1

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CALC_t4D2(t4D2,RR)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 t4D2(3,3,3,3),RR(3,3),TMP(3,3,3,3)
      INTEGER I,J

      t4D2=0.D0

C     DELTAikDELTAjl
      CALL t2at2b_t4(TMP,U,U)
      CALL t4trans(TMP,2,3)
      t4D2=t4D2+TMP
C     DELTAjkDELTAil
      CALL t2at2b_t4(TMP,U,U)
      CALL t4trans(TMP,1,3)
      t4D2=t4D2+TMP
C     DELTAijDELTAkl
      CALL t2at2b_t4(TMP,U,U)
      t4D2=t4D2-TMP
C     DELTAklRiRj
      CALL t2at2b_t4(TMP,RR,U)
      t4D2=t4D2+TMP
C     DELTAijRkRl
      CALL t2at2b_t4(TMP,U,RR)
      t4D2=t4D2+TMP
C     RiRjRkRl
      CALL t2at2b_t4(TMP,RR,RR)
      t4D2=t4D2-3.D0*TMP

      t4D2 = 0.5D0*t4D2

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE t2at2b_t4(t4,t2a,t2b)
      IMPLICIT NONE
      REAL*8 t4(3,3,3,3),t2a(3,3),t2b(3,3)
      INTEGER I,J,K,L

      t4=0.D0

      DO I=1,3
       DO J=1,3
        DO K=1,3
         DO L=1,3
          t4(I,J,K,L)=t2a(I,J)*t2b(K,L)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

C     transposes two indices of 3x3x3x3 matrix i1,i2, rule i1<i2

      SUBROUTINE t4trans(t4,i1,i2)
      IMPLICIT NONE
      REAL*8 t4(3,3,3,3),t4old(3,3,3,3)
      INTEGER I,J,K,L,i1,i2

      t4old=t4 
     
      IF ((i1.EQ.1).AND.(i2.EQ.2)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(J,I,K,L)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.1).AND.(i2.EQ.3)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(K,J,I,L)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.1).AND.(i2.EQ.4)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(L,J,K,I)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.2).AND.(i2.EQ.3)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(I,K,J,L)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.2).AND.(i2.EQ.4)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(I,L,K,J)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ELSE IF ((i1.EQ.3).AND.(i2.EQ.4)) THEN
       DO I=1,3
        DO J=1,3
         DO K=1,3
          DO L=1,3
           t4(I,J,K,L)=t4old(I,J,L,K)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF
 
      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE mulT2aT4T2b(t2a,t4,t2b,a)
      IMPLICIT NONE
      REAL*8 t4(3,3,3,3),t2a(3,3),t2b(3,3),a
      INTEGER I,J,K,L

      a=0.D0

      DO I=1,3
       DO J=1,3
        DO K=1,3
         DO L=1,3
          a = a + t4(I,J,K,L)*t2a(I,J)*t2b(K,L)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE mulT3T2(t3,t2,v)
      IMPLICIT NONE
      REAL*8 t3(3,3,3),t2(3,3),v(3)
      INTEGER I,J

      v=0.D0

      DO I=1,3
       DO J=1,3
          v = v + t3(1:3,I,J)*t2(I,J)
       ENDDO
      ENDDO

      RETURN
      END

     
