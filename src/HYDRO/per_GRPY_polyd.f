
C EVALUATES THE 11NNx11NN ROTNE-PRAGER MOBILITY MATRIX.
C DISPLAYS BLOKS:
C PP  6NNx6NN
C PQ  6NNx5NN
C QQ  5NNx5NN


      SUBROUTINE GRPERY_MOB(APP,APQ,AQQ,CONF,RADII,NN,LATTICE,EWS_ERR)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN),ARR(6*NN,6*NN)
      REAL*8 CONF(3,NN),RADII(NN),LATTICE(3,3),EWS_ERR

      CALL GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN,LATTICE,EWS_ERR)

      ARR = APP
      CALL INVFRI_TO_FRI(ARR,APQ,AQQ,NN)
      CALL FRI_TO_MOB_RED(APP,APQ,AQQ,NN)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE GRPERY_INV_FRI(APP,APQ,AQQ,CONF,RADII,NN,
     *               LATTICE,EWS_ERR)
      USE TENSORS
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN),CONF(3,NN)
      REAL*8 RADII(NN),aaI,aaJ,LATTICE(3,3)
      REAL*8 A1PP(6,6),A1PQ(6,5),A1QQ(5,5)
      REAL*8 C1PP(3,3),C1PQ(3,5),C1QQ(5,5)
      REAL*8 R(3),DMIN,EWS_ERR
      INTEGER I,J,M

      CALL CALC_LATTICE_INVERSE(LATTICE,EWS_ERR)

      APP=0.D0
      APQ=0.D0
      AQQ=0.D0

      DO I=1,NN

       CALL PER_ROTNE_PRAGER_TT_SELF(C1PP,RADII(I))
       C1PP = C1PP + 1.D0/(6.D0*RADII(I))*U 
       APP(6*(I-1)+1:6*(I-1)+3,6*(I-1)+1:6*(I-1)+3) = C1PP

       CALL PER_ROTNE_PRAGER_RR_SELF(C1PP,RADII(I))
       C1PP = C1PP + 1.D0/(8.D0*RADII(I)**3)*U 
       APP(6*(I-1)+4:6*(I-1)+6,6*(I-1)+4:6*(I-1)+6) = C1PP

       CALL PER_ROTNE_PRAGER_DD_SELF_Y2(C1QQ,RADII(I))
       DO J=1,5
        C1QQ(J,J) = C1QQ(J,J) + 3.D0/(20.D0*RADII(I)**3) 
       ENDDO            
       AQQ(5*(I-1)+1:5*(I-1)+5,5*(I-1)+1:5*(I-1)+5) = C1QQ

      ENDDO

      DO I=1,NN-1           
       DO J=I+1,NN

        aaI=RADII(I)
        aaJ=RADII(J)
        R=CONF(:,I)-CONF(:,J)
        CALL PER_SKEW_CORRECTION(R,DMIN)

        A1PP = 0.D0
        CALL PER_ROTNE_PRAGER_TT_IJ(C1PP,R,aaI,aaJ)
        A1PP(1:3,1:3) = C1PP
        CALL PER_ROTNE_PRAGER_RT_IJ(C1PP,R,aaI,aaJ)
        A1PP(4:6,1:3) = C1PP
        CALL PER_ROTNE_PRAGER_RT_IJ(C1PP,-R,aaJ,aaI)
        A1PP(1:3,4:6) = -C1PP
        CALL PER_ROTNE_PRAGER_RR_IJ(C1PP,R,aaI,aaJ)
        A1PP(4:6,4:6) = C1PP
        APP( 6*(I-1)+1:6*I , 6*(J-1)+1:6*J )=A1PP
        APP( 6*(J-1)+1:6*J , 6*(I-1)+1:6*I )=TRANSPOSE(A1PP)

        CALL PER_ROTNE_PRAGER_TD_IJ_Y2(C1PQ,R,aaI,aaJ)
        A1PQ(1:3,1:5) = C1PQ
        CALL PER_ROTNE_PRAGER_RD_IJ_Y2(C1PQ,R,aaI,aaJ)
        A1PQ(4:6,1:5) = C1PQ
        APQ(6*(I-1)+1:6*(I-1)+6,5*(J-1)+1:5*(J-1)+5) = A1PQ
        CALL PER_ROTNE_PRAGER_TD_IJ_Y2(C1PQ,-R,aaI,aaJ)
        A1PQ(1:3,1:5) = C1PQ
        CALL PER_ROTNE_PRAGER_RD_IJ_Y2(C1PQ,-R,aaI,aaJ)
        A1PQ(4:6,1:5) = C1PQ
        APQ(6*(J-1)+1:6*(J-1)+6,5*(I-1)+1:5*(I-1)+5) = A1PQ

        CALL PER_ROTNE_PRAGER_DD_IJ_Y2(C1QQ,R,aaI,aaJ)              
        AQQ(5*(I-1)+1:5*(I-1)+5,5*(J-1)+1:5*(J-1)+5) = C1QQ
        AQQ(5*(J-1)+1:5*(J-1)+5,5*(I-1)+1:5*(I-1)+5) = C1QQ 
       
       ENDDO
      ENDDO

C  IF PARTICLES OVERLAP WE SUBTRACT OSEEN AND ADD YAMAKAWA

      DO I=1,NN-1           
       DO J=I+1,NN

        aaI=RADII(I)
        aaJ=RADII(J)
        R=CONF(:,I)-CONF(:,J)
        CALL PER_SKEW_CORRECTION(R,DMIN)

        IF(DMIN<=(aaI+aaJ)) THEN
 
         CALL CORRECTION_TT_IJ(C1PP,R,aaI,aaJ)              
         APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+1:6*(J-1)+3) =
     * APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+1:6*(J-1)+3) + C1PP
         APP(6*(J-1)+1:6*(J-1)+3,6*(I-1)+1:6*(I-1)+3) =
     * APP(6*(J-1)+1:6*(J-1)+3,6*(I-1)+1:6*(I-1)+3) +  C1PP
         
         CALL CORRECTION_RR_IJ(C1PP,R,aaI,aaJ)              
         APP(6*(I-1)+4:6*(I-1)+6,6*(J-1)+4:6*(J-1)+6) =
     * APP(6*(I-1)+4:6*(I-1)+6,6*(J-1)+4:6*(J-1)+6) + C1PP
         APP(6*(J-1)+4:6*(J-1)+6,6*(I-1)+4:6*(I-1)+6) =
     * APP(6*(J-1)+4:6*(J-1)+6,6*(I-1)+4:6*(I-1)+6) + C1PP

         CALL CORRECTION_RT_IJ(C1PP,R,aaI,aaJ)              
         APP(6*(I-1)+4:6*(I-1)+6,6*(J-1)+1:6*(J-1)+3) =        ! RT IJ
     * APP(6*(I-1)+4:6*(I-1)+6,6*(J-1)+1:6*(J-1)+3) + C1PP     ! RT IJ
         APP(6*(J-1)+1:6*(J-1)+3,6*(I-1)+4:6*(I-1)+6) =        ! TR JI
     * APP(6*(J-1)+1:6*(J-1)+3,6*(I-1)+4:6*(I-1)+6) - C1PP     ! TR JI
   
         CALL CORRECTION_RT_IJ(C1PP,-R,aaJ,aaI)              
         APP(6*(J-1)+4:6*(J-1)+6,6*(I-1)+1:6*(I-1)+3) =        ! RT JI
     * APP(6*(J-1)+4:6*(J-1)+6,6*(I-1)+1:6*(I-1)+3) + C1PP     ! RT JI
         APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+4:6*(J-1)+6) =        ! TR IJ
     * APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+4:6*(J-1)+6) - C1PP     ! TR IJ
  
         CALL CORRECTION_DD_IJ_Y2(C1QQ,R,aaI,aaJ)              
         AQQ(5*(I-1)+1:5*(I-1)+5,5*(J-1)+1:5*(J-1)+5) =        ! DD IJ
     * AQQ(5*(I-1)+1:5*(I-1)+5,5*(J-1)+1:5*(J-1)+5) + C1QQ     ! DD IJ
         AQQ(5*(J-1)+1:5*(J-1)+5,5*(I-1)+1:5*(I-1)+5) =        ! DD JI
     * AQQ(5*(J-1)+1:5*(J-1)+5,5*(I-1)+1:5*(I-1)+5) + C1QQ     ! DD JI
   
         CALL CORRECTION_TD_IJ_Y2(C1PQ,R,aaI,aaJ)
         APQ(6*(I-1)+1:6*(I-1)+3,5*(J-1)+1:5*(J-1)+5) =        ! TD IJ
     * APQ(6*(I-1)+1:6*(I-1)+3,5*(J-1)+1:5*(J-1)+5) - C1PQ     ! TD IJ

         CALL CORRECTION_TD_IJ_Y2(C1PQ,-R,aaJ,aaI)
         APQ(6*(J-1)+1:6*(J-1)+3,5*(I-1)+1:5*(I-1)+5) =        ! TD JI
     * APQ(6*(J-1)+1:6*(J-1)+3,5*(I-1)+1:5*(I-1)+5) - C1PQ     ! TD JI
   
         CALL CORRECTION_RD_IJ_Y2(C1PQ,R,aaI,aaJ)
         APQ(6*(I-1)+4:6*(I-1)+6,5*(J-1)+1:5*(J-1)+5) =        ! RD IJ
     * APQ(6*(I-1)+4:6*(I-1)+6,5*(J-1)+1:5*(J-1)+5) - C1PQ     ! RD IJ
 
         CALL CORRECTION_RD_IJ_Y2(C1PQ,-R,aaJ,aaI)
         APQ(6*(J-1)+4:6*(J-1)+6,5*(I-1)+1:5*(I-1)+5) =        ! RD JI
     * APQ(6*(J-1)+4:6*(J-1)+6,5*(I-1)+1:5*(I-1)+5) - C1PQ     ! RD JI
 
        ENDIF

       ENDDO
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CALC_LATTICE_INVERSE(LATTICE,EWS_ERR)
      USE LATTICE_BASE   ! LB(3)
      USE LATTICE_SKEW
      USE FORCE_PAR      ! GAMMA
      IMPLICIT NONE
      REAL*8 MM(3,3),LE(3),LATTICE(3,3)  ! MM = METRIC MATRIX, LE = EIGENVALUS
      REAL*8 RNR,RNI,EWS_ERR
      REAL*8 LMIN
      REAL*8 TSHEAR,PI
      SAVE TSHEAR,PI
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN      
       INI=.FALSE.
       PI=ACOS(-1.D0)
      ENDIF

      LR = LATTICE

      CALL MAT_INV3(LI,LR,LV)
      LV=ABS(LV)
      SIG=(LV)**(1.D0/3.D0)/SQRT(2.D0*PI)
      LOGERR=-LOG(EWS_ERR)
      
      MM=MATMUL(TRANSPOSE(LR),LR)      
      CALL EIGEN(MM,LE,3)
      LMIN=MINVAL(LE)
      RNR=2.D0*LOGERR*SIG**2/LMIN
      NR=SQRT(RNR)+1

      LMIN=1.D0/MAXVAL(LE)
      RNI=LOGERR/(2.D0*PI**2*SIG**2*LMIN)
      NI=SQRT(RNI)+1

      RETURN
      END



C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE EIGEN(A,D,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 A(NN,NN)
      REAL*8 D(NN),E(NN-1),TAU(NN-1),WORK(2*NN)
      INTEGER INFO

      CALL DSYTRD('U', NN, A, NN, D, E, TAU, WORK, 2*NN, INFO )
      IF(INFO.NE.0) WRITE(*,*) 'EIGEN ERROR'

      CALL DORGTR('U', NN, A, NN, TAU, WORK, 2*NN, INFO )
      IF(INFO.NE.0) WRITE(*,*) 'EIGEN ERROR'

      CALL DSTEQR('V', NN, D, E, A, NN, WORK, INFO )
      IF(INFO.NE.0) WRITE(*,*) 'EIGEN ERROR'

      RETURN
      END
      
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE MAT_INV3(LI,LR,DET_LR)
      IMPLICIT NONE
      REAL*8 LI(3,3),LR(3,3),DET_LR
      
      det_LR=LR(1,1)*LR(2,2)*LR(3,3)+LR(1,2)*LR(2,3)*LR(3,1)+
     *        LR(1,3)*LR(2,1)*LR(3,2)-LR(1,3)*LR(2,2)*LR(3,1)-      
     *        LR(2,3)*LR(3,2)*LR(1,1)-LR(3,3)*LR(1,2)*LR(2,1)

      LI(1,1)=(1.D0/det_LR)*(LR(2,2)*LR(3,3)-LR(2,3)*LR(3,2))
      LI(1,2)=(1.D0/det_LR)*(LR(1,3)*LR(3,2)-LR(1,2)*LR(3,3))
      LI(1,3)=(1.D0/det_LR)*(LR(1,2)*LR(2,3)-LR(1,3)*LR(2,2))

      LI(2,1)=(1.D0/det_LR)*(LR(2,3)*LR(3,1)-LR(2,1)*LR(3,3))
      LI(2,2)=(1.D0/det_LR)*(LR(1,1)*LR(3,3)-LR(1,3)*LR(3,1))
      LI(2,3)=(1.D0/det_LR)*(LR(1,3)*LR(2,1)-LR(1,1)*LR(2,3))

      LI(3,1)=(1.D0/det_LR)*(LR(2,1)*LR(3,2)-LR(2,2)*LR(3,1))
      LI(3,2)=(1.D0/det_LR)*(LR(1,2)*LR(3,1)-LR(1,1)*LR(3,2))
      LI(3,3)=(1.D0/det_LR)*(LR(1,1)*LR(2,2)-LR(1,2)*LR(2,1))
      
      RETURN
      END
                  
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE DISTANCE_PER_3D(CONF)
      USE SIZE       ! NN,LLG
      USE LATTICE_SKEW 
      IMPLICIT NONE
      REAL*8 CONF(3,NN)
      REAL*8 R(3),DIST
      INTEGER I,J

      DO I=1,NN-1
      DO J=I+1,NN
      R=CONF(:,I)-CONF(:,J)      
      CALL PER_SKEW_CORRECTION(R,DIST)      
      IF( DIST.LE.1.D0 ) THEN
       WRITE(*,*) 'WARNING! Periodic distance between particles'
       WRITE(*,*) I,' and ',J,' is less or equal 1.0'
      ENDIF
      ENDDO
      ENDDO
      
      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_SKEW_CORRECTION(R,DMIN)
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 R(3)
      REAL*8 RC(3),RL(3),RMIN(3),X(3),M(3),DMIN,D
      INTEGER I,M1,M2,M3
      
      X=MATMUL(LI,R)

      X=MOD(X,1.D0)

      DO I=1,3
      IF(X(I).LT.0.D0) X(I)=X(I)+1.D0
      ENDDO

      RC=MATMUL(LR,X)
      DMIN=SUM(RC**2)
      RMIN=RC

      DO M1=0,1
      DO M2=0,1
      DO M3=0,1
      IF(M1==0.AND.M2==0.AND.M3==0) CYCLE
       M(1)=M1
       M(2)=M2
       M(3)=M3
       RL=MATMUL(LR,M)
       D=SUM( (RC-RL)**2 )
       IF(D.LT.DMIN) THEN
        DMIN=D
        RMIN=RC-RL
       ENDIF
      ENDDO
      ENDDO
      ENDDO

      R=RMIN
      DMIN=SQRT(DMIN)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CORRECTION_TT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3),ARP(3,3),AY(3,3)
      REAL*8 aaI,aaJ           ! radius

      A1=0.D0
      ARP=0.D0
      AY=0.D0

      CALL ROTNE_PRAGER_TT_IJ(ARP,R,aaI,aaJ)
      CALL YAMAKAWA_TT_IJ(AY,R,aaI,aaJ)

      A1=AY-ARP

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CORRECTION_RR_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3),ARP(3,3),AY(3,3)
      REAL*8 aaI,aaJ           ! radius

      A1=0.D0
      ARP=0.D0
      AY=0.D0

      CALL ROTNE_PRAGER_RR_IJ(ARP,R)
      CALL YAMAKAWA_RR_IJ(AY,R,aaI,aaJ)

      A1=AY-ARP

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CORRECTION_RT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3),ARP(3,3),AY(3,3)
      REAL*8 aaI,aaJ           ! radius

      A1=0.D0
      ARP=0.D0
      AY=0.D0

      CALL ROTNE_PRAGER_RT_IJ(ARP,R)
      CALL YAMAKAWA_RT_IJ(AY,R,aaI,aaJ)

      A1=AY-ARP

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CORRECTION_DD_IJ_Y2(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1(5,5),R(3),BRP(5,5),BY(5,5)
      REAL*8 aaI,aaJ           ! radius

      B1=0.D0
      BRP=0.D0
      BY=0.D0

      CALL ROTNE_PRAGER_DD_IJ_Y2(BRP,R,aaI,aaJ)
      CALL YAMAKAWA_DD_IJ_Y2(BY,R,aaI,aaJ)

      B1=BY-BRP

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CORRECTION_TD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,5),R(3),CRP(3,5),CY(3,5)
      REAL*8 aaI,aaJ           ! radius

      C1=0.D0
      CRP=0.D0
      CY=0.D0

      CALL ROTNE_PRAGER_TD_IJ_Y2(CRP,R,aaI,aaJ)
      CALL YAMAKAWA_TD_IJ_Y2(CY,R,aaI,aaJ)

      C1=CY-CRP

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE CORRECTION_RD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,5),R(3),CRP(3,5),CY(3,5)
      REAL*8 aaI,aaJ           ! radius

      C1=0.D0
      CRP=0.D0
      CY=0.D0

      CALL ROTNE_PRAGER_RD_IJ_Y2(CRP,R,aaI,aaJ)
      CALL YAMAKAWA_RD_IJ_Y2(CY,R,aaI,aaJ)

      C1=CY-CRP

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_TT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST3,RW(3),RR(3,3)
      REAL*8 aaI,aaJ           ! radius
      REAL*8 P18,P23a2
      PARAMETER(P18=1.D0/8.D0)
      P23a2=(aaI**2 + aaJ**2)/3.D0

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST3=DIST**3

      CALL CALC_RR(RR,RW)

      A1=P18*( (1.D0/DIST+P23a2/DIST3)*U + (1/DIST-3.D0*P23a2/DIST3)*RR)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_TT_SELF(A1,aaI)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 A1(3,3),A2A(3,3),A2B(3,3)
      REAL*8 DIST,DIST3,RW(3),RR(3,3)
      REAL*8 aaI           ! radius
      REAL*8 PI,PRR,PRE
      INTEGER NX,NY,NZ
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      A1=0.D0

C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_TT_IJ(A2A,RN,aaI,aaI)
           A2A = A2A*ERFC(DIST/SIG/SQRT(2.D0))

           RW=RN/DIST
           CALL CALC_RR(RR,RW)
           PI =(aaI**2 + aaI**2)*
     *         (1.D0/(6.D0*SIG**2) + 1.D0/(3.D0*DIST**2))
           PRR=((aaI**2 + aaI**2)*
     *        (DIST**2/(6.D0*SIG**4) - 1.D0/(3.D0*SIG**2)
     *          - 1.D0/DIST**2  )
     *        + 1.D0)
           PRE=EXP(-DIST**2/(2.D0*SIG**2))/
     *         (4.D0*SQRT(2*P)*SIG)
           A2B=(PI*U + PRR*RR)*PRE
           A1 = A1 + A2A + A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      A1 = A1 - (SIG**2)*P/(2.D0*LV)*U
     * + (aaI**2/(9.D0*SIG**2)-1)/(4.D0*SQRT(2*P)*SIG)*U

C INVERSE LATTICE
      A2A = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST
           CALL CALC_RR(RR,RW)
           PI = 1.D0
           PRR = (1.D0 + (SIG**2)*(DIST**2)/2.D0)
           PRE = (P/LV)*
     *          (1.D0/(DIST**2) - (aaI**2+aaI**2)/6.D0)
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
           A2B = PRE*(PI*U - PRR*RR)
           A2A = A2A + A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      A1 = A1 + A2A

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_TT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 A1(3,3),A2A(3,3),A2B(3,3)
      REAL*8 DIST,DIST3,RW(3),RR(3,3),R(3)
      REAL*8 aaI,aaJ           ! radius
      REAL*8 PI,PRR,PRE
      INTEGER NX,NY,NZ
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      A1=0.D0

C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=R+MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_TT_IJ(A2A,RN,aaI,aaJ)
           A2A = A2A*ERFC(DIST/SIG/SQRT(2.D0))
           RW=RN/DIST
           CALL CALC_RR(RR,RW)
           PI =(aaI**2 + aaJ**2)*
     *         (1.D0/(6.D0*SIG**2) + 1.D0/(3.D0*DIST**2))
           PRR=((aaI**2 + aaJ**2)*
     *        (DIST**2/(6.D0*SIG**4) - 1.D0/(3.D0*SIG**2)
     *          - 1.D0/DIST**2  )
     *        + 1.D0)
           PRE=EXP(-DIST**2/(2.D0*SIG**2))/
     *         (4.D0*SQRT(2*P)*SIG)
           A2B=(PI*U + PRR*RR)*PRE
           A1 = A1 + A2A + A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      A1 = A1 - (SIG**2)*P/(2.D0*LV)*U

C INVERSE LATTICE
      A2A = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST
           CALL CALC_RR(RR,RW)
           PI = 1.D0
           PRR = (1.D0 + (SIG**2)*(DIST**2)/2.D0)
           PRE = (P/LV)*
     *          (1.D0/(DIST**2) - (aaI**2+aaJ**2)/6.D0)
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
     *           *COS(SUM(R*KN))
           A2B = PRE*(PI*U - PRR*RR)
           A2A = A2A + A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      A1 = A1 + A2A

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_TT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST2,DIST3,RW(3),RR(3,3)
      REAL*8 aaI,aaJ            ! radius
      REAL*8 M0TT          ! tt single sphere mobility 
      PARAMETER(M0TT=1.D0/6.D0)
      REAL*8 PRE,PU,PRR,aaIaaJ

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_RR(RR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=M0TT/aaJ*U
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        PRE=1.D0/(6*32*aaI*aaJ)
        PU=16*(aaI+aaJ)-(aaIaaJ**2+3*DIST2)**2/DIST3
        PRR=3*(aaIaaJ**2 - DIST2)**2/DIST3
        A1=PRE*(PU*U+PRR*RR)
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=M0TT/aaI*U
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        PRE=1.D0/(6*32*aaI*aaJ)
        PU=16*(aaI+aaJ)-(aaIaaJ**2+3*DIST2)**2/DIST3
        PRR=3*(aaIaaJ**2 - DIST2)**2/DIST3
        A1=PRE*(PU*U+PRR*RR)
       ENDIF       
      ENDIF

      RETURN
      END
      
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_RR_IJ(A1,R)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST3,RW(3),RR(3,3)
      REAL*8 P16
      PARAMETER(P16=1.D0/16.D0)

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST3=DIST**3

      CALL CALC_RR(RR,RW)

      A1=-P16*(U/DIST3-3.D0*RR/DIST3)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_RR_SELF(A1,aaI)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 A1(3,3),A2A(3,3),A2B(3,3)
      REAL*8 DIST,DIST3,RW(3),RR(3,3)
      REAL*8 aaI           ! radius
      REAL*8 PI,PRR,PRE
      INTEGER NX,NY,NZ
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      A1=0.D0

C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_RR_IJ(A2A,RN)
           A2A = A2A*ERFC(DIST/SIG/SQRT(2.D0))

           RW=RN/DIST
           CALL CALC_RR(RR,RW)
           PI =(1.D0/(4.D0*SIG**2) + 1.D0/(2.D0*DIST**2))
           PRR= (
     *          DIST**2/(4.D0*SIG**4) - 1.D0/(2.D0*SIG**2)
     *          - 3.D0/(2.D0*DIST**2))  
           PRE=EXP(-DIST**2/(2.D0*SIG**2))/
     *         (4.D0*SQRT(2*P)*SIG)
           A2B=(PI*U + PRR*RR)*PRE
           A1 = A1 + A2A - A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      A1 = A1 - 1.D0/(48.D0*SQRT(2.D0*P)*SIG**3)*U

C INVERSE LATTICE
      A2A = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST
           CALL CALC_RR(RR,RW)
           PI = 1.D0
           PRR = (1.D0 + (SIG**2)*(DIST**2)/2.D0)
           PRE = (P/(4.D0*LV))
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
           A2B = PRE*(PI*U - PRR*RR)
           A2A = A2A + A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      A1 = A1 + A2A

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_RR_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 A1(3,3),A2A(3,3),A2B(3,3)
      REAL*8 DIST,DIST3,RW(3),RR(3,3),R(3)
      REAL*8 aaI,aaJ           ! radius
      REAL*8 PI,PRR,PRE
      INTEGER NX,NY,NZ
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      A1=0.D0

C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=R+MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_RR_IJ(A2A,RN)
           A2A = A2A*ERFC(DIST/SIG/SQRT(2.D0))

           RW=RN/DIST
           CALL CALC_RR(RR,RW)
           PI =(1.D0/(4.D0*SIG**2) + 1.D0/(2.D0*DIST**2))
           PRR= (
     *          DIST**2/(4.D0*SIG**4) - 1.D0/(2.D0*SIG**2)
     *          - 3.D0/(2.D0*DIST**2))  
           PRE=EXP(-DIST**2/(2.D0*SIG**2))/
     *         (4.D0*SQRT(2*P)*SIG)
           A2B=(PI*U + PRR*RR)*PRE
           A1 = A1 + A2A - A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO

C INVERSE LATTICE
      A2A = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST
           CALL CALC_RR(RR,RW)
           PI = 1.D0
           PRR = (1.D0 + (SIG**2)*(DIST**2)/2.D0)
           PRE = (P/(4.D0*LV))
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
     *           *COS(SUM(R*KN))
           A2B = PRE*(PI*U - PRR*RR)
           A2A = A2A + A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      A1 = A1 + A2A

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_RR_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST2,DIST3,RW(3),RR(3,3)
      REAL*8 aaI,aaJ            ! radius
      REAL*8 M0RR          ! rr single sphere mobility 
      PARAMETER(M0RR=1.D0/8.D0)
      REAL*8 PRE,PU,PRR,aaIaaJ,aaI2,aaJ2,aaI3,aaJ3

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_RR(RR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=M0RR/aaJ**3*U
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        aaI2=aaI**2
        aaJ2=aaJ**2
        aaI3=aaI**3
        aaJ3=aaJ**3
        PRE=1.D0/(8*64*aaI3*aaJ3)

        PU=  5*DIST3 - 27*DIST*(aaI2+aaJ2) 
     *     + 32*(aaI3+aaJ3)
     *     - 9*(aaI2-aaJ2)**2/DIST 
     *     - aaIaaJ**4*(aaI2 + 4*aaI*aaJ + aaJ2)/DIST3

        PRR=3*(aaIaaJ**2 - DIST2)**2
     *    *((aaI2 + 4*aaI*aaJ + aaJ2)-DIST2)/DIST3

        A1=PRE*(PU*U+PRR*RR)
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=M0RR/aaI**3*U
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        aaI2=aaI**2
        aaJ2=aaJ**2
*-----------------------------------------------------------
        aaI3=aaI**3
        aaJ3=aaJ**3
*-----------------------------------------------------------
        PRE=1.D0/(8*64*(aaI**3)*(aaJ**3))

        PU=  5*DIST3 - 27*DIST*(aaI2+aaJ2) 
     *     + 32*(aaI3+aaJ3)
     *     - 9*(aaI2-aaJ2)**2/DIST 
     *     - aaIaaJ**4*(aaI2 + 4*aaI*aaJ + aaJ2)/DIST3

        PRR=3*(aaIaaJ**2 - DIST2)**2
     *    *((aaI2 + 4*aaI*aaJ + aaJ2)-DIST2)/DIST3

        A1=PRE*(PU*U+PRR*RR)
       ENDIF       
      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_RT_IJ(A1,R)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST2,RW(3),EPSR(3,3)
      REAL*8 P18
      PARAMETER(P18=1.D0/8.D0)
      
      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST2=DIST**2
      
      CALL CALC_EPSR(EPSR,RW)

      A1=P18*EPSR/DIST2

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_RT_SELF(A1,aaI)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 A1(3,3),A2A(3,3),A2B(3,3)
      REAL*8 DIST,DIST3,RW(3),EPSR(3,3)
      REAL*8 aaI           ! radius
      REAL*8 PI,PRR,PRE
      INTEGER NX,NY,NZ
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      A1=0.D0

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_RT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 A1(3,3),A2A(3,3),A2B(3,3)
      REAL*8 DIST,DIST3,RW(3),EPSR(3,3),R(3)
      REAL*8 aaI,aaJ           ! radius
      REAL*8 PEPSR,PRE
      INTEGER NX,NY,NZ
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      A1=0.D0

C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=R+MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_RT_IJ(A2A,RN)
           A2A = A2A*ERFC(DIST/SIG/SQRT(2.D0))

           RW=RN/DIST
           CALL CALC_EPSR(EPSR,RW)
           PEPSR = 1.D0/DIST
           PRE=EXP(-DIST**2/(2.D0*SIG**2))/
     *         (4.D0*SQRT(2*P)*SIG)
           A2B= (PEPSR*EPSR)*PRE
           A1 = A1 + A2A + A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO

C INVERSE LATTICE
      A2A = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST
           CALL CALC_EPSR(EPSR,RW)
           PEPSR =1.D0/DIST
           PRE = (P/(2.D0*LV))
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
     *           *SIN(SUM(R*KN))
           A2B= (PEPSR*EPSR)*PRE
           A2A = A2A + A2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      A1 = A1 + A2A

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_RT_IJ(A1,R,aaI,aaJ)
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,RW(3),EPSR(3,3)
      REAL*8 aaI,aaJ            ! radius
      REAL*8 PRE,PEPSR

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_EPSR(EPSR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=0.D0
       ELSE
        PRE=1.D0/(16*8*aaI**3*aaJ)
        PEPSR=((aaI-aaJ)+DIST)**2
     *   *(aaJ**2+2*aaJ*(aaI+DIST)-3*(aaI-DIST)**2)
     *   /DIST**2
        A1=PRE*PEPSR*EPSR
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=DIST/(8*(aaI**3))*EPSR
       ELSE
        PRE=1.D0/(16*8*aaI**3*aaJ)
        PEPSR=((aaI-aaJ)+DIST)**2
     *   *(aaJ**2+2*aaJ*(aaI+DIST)-3*(aaI-DIST)**2)
     *   /DIST**2
        A1=PRE*PEPSR*EPSR
       ENDIF       
      ENDIF
      
      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_DD_IJ_Y2(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1_CART(3,3,3,3),B1(5,5),R(3),a
      REAL*8 aaI,aaJ           ! radius
      INTEGER I,J

      CALL ROTNE_PRAGER_DD_IJ(B1_CART,R,aaI,aaJ)

      DO I=1,5
       DO J=1,5
        CALL mulT2aT4T2b(Y2(I,1:3,1:3),B1_CART,Y2(J,1:3,1:3),a) 
        B1(I,J)=a
       ENDDO
      ENDDO

      RETURN
      END
     
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_DD_IJ(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1(3,3,3,3),t4D1(3,3,3,3),t4D2(3,3,3,3),t4D0(3,3,3,3)
      REAL*8 DIST,DIST2,DIST5,RW(3),RR(3,3),R(3)
      REAL*8 aaI,aaJ,aaI2,aaJ2          ! radius

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST5=DIST**5
      DIST2=DIST**2
      aaI2=aaI**2
      aaJ2=aaJ**2

      CALL CALC_RR(RR,RW)

      CALL CALC_t4D0(t4D0,RR)
      CALL CALC_t4D1(t4D1,RR)
      CALL CALC_t4D2(t4D2,RR)
      t4D2 = t4D2 - t4D1

      B1=0.D0

      B1=B1 +3.D0
     *      *(6.D0*(aaI2+aaJ2)-5.D0*DIST2)
     *      /(20.D0*DIST5)
     *      *t4D0

      B1=B1 -3.D0
     *      *(8.D0*(aaI2+aaJ2)-5.D0*DIST2)
     *      /(40.D0*DIST5)
     *      *t4D1

      B1=B1 +3.D0
     *      *(aaI2+aaJ2)
     *      /(20.D0*DIST5)
     *      *t4D2

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_DD_SELF_Y2(B1,aaI)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1_CART(3,3,3,3),B1(5,5),R(3),a
      REAL*8 aaI           ! radius
      INTEGER I,J

      CALL PER_ROTNE_PRAGER_DD_SELF(B1_CART,aaI)
      
      DO I=1,5
       DO J=1,5
        CALL mulT2aT4T2b(Y2(I,1:3,1:3),B1_CART,Y2(J,1:3,1:3),a) 
        B1(I,J)=a
       ENDDO
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_DD_SELF(B1,aaI)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 B1(3,3,3,3),t4D1(3,3,3,3),t4D2(3,3,3,3),t4D0(3,3,3,3)
      REAL*8 B2A(3,3,3,3),B2B(3,3,3,3),B2C(3,3,3,3),B1Y(5,5)
      REAL*8 DIST,DIST2,DIST5,RW(3),RR(3,3),R(3)
      REAL*8 aaI          ! radius
      REAL*8 PD0,PD1,PD2,PRE
      INTEGER NX,NY,NZ,I,J
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      B1=0.D0
      B2A=0.D0
      B2B=0.D0
C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_DD_IJ(B2C,RN,aaI,aaI)
           B2A = B2A + B2C*ERFC(DIST/SIG/SQRT(2.D0))

           RW=RN/DIST
           CALL CALC_RR(RR,RW)
           CALL CALC_t4D0(t4D0,RR)
           CALL CALC_t4D1(t4D1,RR)
           CALL CALC_t4D2(t4D2,RR)
           t4D2 = t4D2 - t4D1

           PD0 =(aaI**2 + aaI**2)*
     *  (
     *   - DIST**4/(15.D0*SIG**8) + 4.D0*DIST**2/(15.D0*SIG**6)
     *   + 2.D0/(5.D0*SIG**4) + 12.D0/(5.D0*(SIG**2)*DIST**2)
     *   + 36.D0/(5.D0*DIST**4)
     *  )
     *   - (2.D0*DIST**2)/(3.D0*SIG**4) - 2.D0/SIG**2 - 6.D0/DIST**2   

           PD1=(aaI**2 + aaI**2)/10.D0*
     *  (     
     *     DIST**2/SIG**6 - 4.D0/SIG**4
     *   - 16.D0/((SIG**2)*DIST**2) - 48.D0/DIST**4
     *  )
     *   + 3.D0/DIST**2 + 1.D0/SIG**2

           PD2=2.D0/5.D0*(aaI**2 + aaI**2)*
     *  (
     *   1.D0/((SIG**2)*DIST**2) + 3.D0/DIST**4
     *  )

           PRE=EXP(-DIST**2/(2.D0*SIG**2))/
     *         (4.D0*SQRT(2*P)*SIG)

           B2B = B2B + (PD0*t4D0 + PD1*t4D1 + PD2*t4D2)*PRE

         ENDIF
        ENDDO
       ENDDO
      ENDDO
      B1 = B1 + B2A + B2B

      B1 = B1 - aaI**2/(25.D0*SQRT(2.D0*P)*SIG**5)*II

C INVERSE LATTICE
      B2C = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST

           CALL CALC_RR(RR,RW)
           CALL CALC_t4D0(t4D0,RR)
           CALL CALC_t4D1(t4D1,RR)

           PD0 = -2.D0/3.D0*((DIST**2)*SIG**2)
           PD1 = 1.D0

           PRE = (P/(20.D0*LV))
     *      *(10.D0 - (aaI**2 + aaI**2)*DIST**2)
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
           B2B = PRE*(PD0*t4D0 + PD1*t4D1)
           B2C = B2C + B2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      B1 = B1 + B2C

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_DD_IJ_Y2(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1_CART(3,3,3,3),B1(5,5),R(3),a
      REAL*8 aaI,aaJ           ! radius
      INTEGER I,J

      CALL PER_ROTNE_PRAGER_DD_IJ(B1_CART,R,aaI,aaJ)

      DO I=1,5
       DO J=1,5
        CALL mulT2aT4T2b(Y2(I,1:3,1:3),B1_CART,Y2(J,1:3,1:3),a) 
        B1(I,J)=a
       ENDDO
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_DD_IJ(B1,R,aaI,aaJ)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 B1(3,3,3,3),t4D1(3,3,3,3),t4D2(3,3,3,3),t4D0(3,3,3,3)
      REAL*8 B2A(3,3,3,3),B2B(3,3,3,3),B2C(3,3,3,3)
      REAL*8 DIST,DIST2,DIST5,RW(3),RR(3,3),R(3)
      REAL*8 aaI,aaJ          ! radius
      REAL*8 PD0,PD1,PD2,PRE
      INTEGER NX,NY,NZ,I,J
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      B1=0.D0
      B2A=0.D0
      B2B=0.D0
C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=R+MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_DD_IJ(B2C,RN,aaI,aaJ)
           B2A = B2A + B2C*ERFC(DIST/SIG/SQRT(2.D0))

           RW=RN/DIST
           CALL CALC_RR(RR,RW)
           CALL CALC_t4D0(t4D0,RR)
           CALL CALC_t4D1(t4D1,RR)
           CALL CALC_t4D2(t4D2,RR)
           t4D2 = t4D2 - t4D1

           PD0 =(aaI**2 + aaJ**2)*
     *  (
     *   - DIST**4/(15.D0*SIG**8) + 4.D0*DIST**2/(15.D0*SIG**6)
     *   + 2.D0/(5.D0*SIG**4) + 12.D0/(5.D0*(SIG**2)*DIST**2)
     *   + 36.D0/(5.D0*DIST**4)
     *  )
     *   - (2.D0*DIST**2)/(3.D0*SIG**4) - 2.D0/SIG**2 - 6.D0/DIST**2   

           PD1=(aaI**2 + aaJ**2)/10.D0*
     *  (     
     *     DIST**2/SIG**6 - 4.D0/SIG**4
     *   - 16.D0/((SIG**2)*DIST**2) - 48.D0/DIST**4
     *  )
     *   + 3.D0/DIST**2 + 1.D0/SIG**2

           PD2=2.D0/5.D0*(aaI**2 + aaJ**2)*
     *  (
     *   1.D0/((SIG**2)*DIST**2) + 3.D0/DIST**4
     *  )

           PRE=EXP(-DIST**2/(2.D0*SIG**2))/
     *         (4.D0*SQRT(2*P)*SIG)

           B2B = B2B + (PD0*t4D0 + PD1*t4D1 + PD2*t4D2)*PRE
          ENDIF
        ENDDO
       ENDDO
      ENDDO

      B1 = B1 + B2A + B2B

C INVERSE LATTICE
      B2C = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST

           CALL CALC_RR(RR,RW)
           CALL CALC_t4D0(t4D0,RR)
           CALL CALC_t4D1(t4D1,RR)

           PD0 = -2.D0/3.D0*((DIST**2)*SIG**2)
           PD1 = 1.D0

           PRE = (P/(20.D0*LV))
     *      *(10.D0 - (aaI**2 + aaJ**2)*DIST**2)      
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
     *           *COS(SUM(R*KN))
           B2B = PRE*(PD0*t4D0 + PD1*t4D1)
           B2C = B2C + B2B
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      B1 = B1 + B2C

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_DD_IJ_Y2(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1_CART(3,3,3,3),B1(5,5),R(3),a
      REAL*8 aaI,aaJ           ! radius
      INTEGER I,J

      CALL YAMAKAWA_DD_IJ(B1_CART,R,aaI,aaJ)

      DO I=1,5
       DO J=1,5
        CALL mulT2aT4T2b(Y2(I,1:3,1:3),B1_CART,Y2(J,1:3,1:3),a) 
        B1(I,J)=a
       ENDDO
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_DD_IJ(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1(3,3,3,3),t4D1(3,3,3,3),t4D2(3,3,3,3),t4D0(3,3,3,3)
      REAL*8 DIST,DIST3,DIST5,RW(3),RR(3,3),R(3)
      REAL*8 aaI,aaJ,aaI2,aaJ2,aaI3,aaJ3

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST5=DIST**5
      DIST3=DIST**3
      aaI2=aaI**2
      aaJ2=aaJ**2
      aaI3=aaI**3
      aaJ3=aaJ**3

      CALL CALC_RR(RR,RW)

      CALL CALC_t4D0(t4D0,RR)
      CALL CALC_t4D1(t4D1,RR)
      CALL CALC_t4D2(t4D2,RR)
      t4D2 = t4D2 - t4D1

      B1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        B1= 3.D0/(aaJ**3*20.D0)*(t4D0+t4D1+t4D2)
       ELSE

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3)
     *     *(
     *     +3.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     -10.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3
     *     +32.D0*(aaI3+aaJ3)
     *     -30.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *      )*t4D0

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3)
     *     *(
     *     -2.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     +5.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3
     *     -15.D0*(aaI2-aaJ2)**2/DIST
     *     +32.D0*(aaI3+aaJ3)
     *     -25.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *     )*t4D1

      B1=B1+3.D0/(2560.D0*aaI3*aaJ3)
     *     *(
     *     +(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     -30.D0*(aaI2-aaJ2)**2/DIST
     *     +64.D0*(aaI3+aaJ3)
     *     -40.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *     )*t4D2

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        B1= 3.D0/(aaI**3*20.D0)*(t4D0+t4D1+t4D2)
       ELSE

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3)
     *     *(
     *     +3.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     -10.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3
     *     +32.D0*(aaI3+aaJ3)
     *     -30.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *      )*t4D0

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3)
     *     *(
     *     -2.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     +5.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3
     *     -15.D0*(aaI2-aaJ2)**2/DIST
     *     +32.D0*(aaI3+aaJ3)
     *     -25.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *     )*t4D1

      B1=B1+3.D0/(2560.D0*aaI3*aaJ3)
     *     *(
     *     +(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     -30.D0*(aaI2-aaJ2)**2/DIST
     *     +64.D0*(aaI3+aaJ3)
     *     -40.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *     )*t4D2

       ENDIF       
      ENDIF

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_TD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3UR(3,3,3),t3RRR(3,3,3)
      REAL*8 DIST,DIST4,RW(3),R(3)
      REAL*8 aaI,aaJ,aaI2,aaJ2          ! radius

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST4=DIST**4
      aaI2=aaI**2
      aaJ2=aaJ**2

      CALL CALC_UR(t3UR,RW)
      CALL CALC_RRR(t3RRR,RW)

      C1=0.D0

      C1= ( - 2.D0*(5.D0*aaI2+3.D0*aaJ2)/(5.D0*DIST4)*t3UR
     *      + (5.D0*aaI2+3.D0*aaJ2-3.D0*DIST**2)/DIST4*t3RRR
     *    )/(8.D0)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_TD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ           ! radius
      INTEGER I

      CALL ROTNE_PRAGER_TD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_TD_SELF_Y2(C1,aaI)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ           ! radius
      INTEGER I

      C1=0.D0

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_TD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3UR(3,3,3),t3RRR(3,3,3)
      REAL*8 C2A(3,3,3),C2B(3,3,3),C2C(3,3,3)
      REAL*8 DIST,DIST4,RW(3),R(3)
      REAL*8 aaI,aaJ           ! radius
      REAL*8 PP0,PP1,PRE
      INTEGER NX,NY,NZ
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      C1=0.D0
      C2A=0.D0
      C2B=0.D0
C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=R+MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_TD_IJ(C2C,RN,aaI,aaJ)
           C2A = C2A - C2C*ERFC(DIST/SIG/SQRT(2.D0))

           RW=RN/DIST
           CALL CALC_UR(t3UR,RW)
           CALL CALC_RRR(t3RRR,RW)
           PP0 = 2.D0/15.D0
     *  * ( 5.D0*aaI**2 + 3.D0*aaJ**2 )
     *  * (1.D0/((SIG**2)*DIST) + 3.D0/(DIST**3))  
           PP1 = ( 5.D0*aaI**2 + 3.D0*aaJ**2 )
     *  * ( DIST**3/(30.D0*SIG**6) - DIST/(15.D0*SIG**4)
     *   -1.D0/(3.D0*(SIG**2)*DIST) - 1.D0/(DIST**3) )
     *  + DIST/SIG**2 + 3.D0/DIST
           PRE =EXP(-DIST**2/(2.D0*SIG**2))/
     *         (4.D0*SQRT(2*P)*SIG)
           C2B = C2B + (PP0*t3UR + PP1*t3RRR)*PRE
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      C1 = C1 + C2A + C2B

C INVERSE LATTICE
      C2C = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST
           CALL CALC_UR(t3UR,RW)
           CALL CALC_RRR(t3RRR,RW)
           PP0 = 2.D0  
           PP1 = - ( 2.D0 + (SIG**2)*DIST**2)
           PRE = (P/(2.D0*LV))/DIST
     *      *(1.D0 - (5.D0*aaI**2 + 3.D0*aaJ**2)/30.D0*DIST**2)      
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
     *           *SIN(SUM(R*KN))
           C2C = C2C + (PP0*t3UR + PP1*t3RRR)*PRE
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      C1 = C1 + C2C

      RETURN
      END

***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_TD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ           ! radius
      INTEGER I

      CALL PER_ROTNE_PRAGER_TD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE YAMAKAWA_TD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3UR(3,3,3),t3RRR(3,3,3)
      REAL*8 DIST,DIST2,DIST4,RW(3),R(3),PRE,PUR,PRRR
      REAL*8 aaI,aaJ         ! radius

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST2=DIST**2
      DIST4=DIST**4

      CALL CALC_UR(t3UR,RW)
      CALL CALC_RRR(t3RRR,RW)

      C1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        C1=-3.D0*DIST/(20.D0*aaJ**3)*t3UR
       ELSE

        PRE = 1.D0/(aaI*aaJ**3)

        PUR = (10.D0*DIST2 - 24.D0*aaI*DIST
     *       + 15.D0*(aaI**2-aaJ**2) 
     *       - (aaI-aaJ)**5*(aaI+5.D0*aaJ)/DIST4 )/320.D0

        PRRR = ((aaI-aaJ)**2-DIST2)**2
     *        *((aaI-aaJ)*(aaI+5*aaJ)-DIST2)
     *        /(128.D0*DIST4)

        C1 = PRE*(PUR*t3UR+PRRR*t3RRR)

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        C1=0.D0
       ELSE

        PRE = 1.D0/(aaI*aaJ**3)

        PUR = (10.D0*DIST2 - 24.D0*aaI*DIST
     *       + 15.D0*(aaI**2-aaJ**2) 
     *       - (aaI-aaJ)**5*(aaI+5.D0*aaJ)/DIST4 )/320.D0

        PRRR = ((aaI-aaJ)**2-DIST2)**2
     *        *((aaI-aaJ)*(aaI+5*aaJ)-DIST2)
     *        /(128.D0*DIST4)

        C1 = PRE*(PUR*t3UR+PRRR*t3RRR)

       ENDIF       
      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_TD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ           ! radius
      INTEGER I

      CALL YAMAKAWA_TD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_RD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3EPSRR(3,3,3)
      REAL*8 DIST,RW(3),R(3)
      REAL*8 aaI,aaJ          ! radius

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_EPSRR(t3EPSRR,RW)

      C1=0.D0

      C1= - 3.D0/(8.D0*DIST**3)*t3EPSRR

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_RD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ           ! radius
      INTEGER I

      CALL ROTNE_PRAGER_RD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_RD_SELF_Y2(C1,aaI)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ           ! radius
      INTEGER I

      C1=0.D0

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_RD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      USE LATTICE_SKEW
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3EPSRR(3,3,3)
      REAL*8 C2A(3,3,3),C2B(3,3,3),C2C(3,3,3)
      REAL*8 DIST,DIST4,RW(3),R(3)
      REAL*8 aaI,aaJ           ! radius
      REAL*8 PEPSRR,PRE
      INTEGER NX,NY,NZ
      REAL*8 RN(3),KN(3),NL(3)
      REAL*8 P,YY
      SAVE   P
      LOGICAL :: INI=.TRUE.
      SAVE INI

      IF(INI) THEN
       INI=.FALSE.
       P=ACOS(-1.D0)
      ENDIF

      C1=0.D0
      C2A=0.D0
      C2B=0.D0
C REAL LATTICE
      DO NX=-NR,NR
       DO NY=-NR,NR
        DO NZ=-NR,NR
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          RN=R+MATMUL(LR,NL)
          YY=DOT_PRODUCT(RN,RN)/SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(RN**2) )
           CALL ROTNE_PRAGER_RD_IJ(C2C,RN,aaI,aaJ)
           C2A = C2A - C2C*ERFC(DIST/SIG/SQRT(2.D0))

           RW=RN/DIST
           CALL CALC_EPSRR(t3EPSRR,RW)
           PEPSRR = 2.D0
     *  * (1.D0/(8.D0*SIG**2) + 3.D0/(8.D0*DIST**2))  
           PRE =EXP(-DIST**2/(2.D0*SIG**2))/
     *         (SQRT(2*P)*SIG)
           C2B = C2B + PEPSRR*t3EPSRR*PRE
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      C1 = C1 + C2A + C2B 

C INVERSE LATTICE
      C2C = 0.D0
      DO NX=-NI,NI
       DO NY=-NI,NI
        DO NZ=-NI,NI
         IF(NX.EQ.0.AND.NY.EQ.0.AND.NZ.EQ.0) CYCLE
          NL(1)=NX
          NL(2)=NY
          NL(3)=NZ
          KN=2*P*MATMUL(NL,LI)
          YY=DOT_PRODUCT(KN,KN)*SIG**2/2.D0
          IF(YY.LE.LOGERR) THEN
           DIST=SQRT( SUM(KN**2) )
           RW=KN/DIST
           CALL CALC_EPSRR(t3EPSRR,RW)
           PEPSRR = -2.D0
           PRE = (P/(4.D0*LV))
     *           *EXP(-(DIST**2)*(SIG**2)/2.D0)
     *           *COS(SUM(R*KN))
           C2C = C2C + PEPSRR*t3EPSRR*PRE
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      C1 = C1 + C2C

      RETURN
      END

***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PER_ROTNE_PRAGER_RD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ           ! radius
      INTEGER I

      CALL PER_ROTNE_PRAGER_RD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END



C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_RD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3EPSRR(3,3,3)
      REAL*8 DIST,DIST2,RW(3),R(3)
      REAL*8 aaI,aaJ         ! radius

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST2=DIST**2

      CALL CALC_EPSRR(t3EPSRR,RW)

      C1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        C1=0.D0
       ELSE

        C1 = -3.D0/(256.D0*(aaI**3)*(aaJ**3)*DIST**3)
     *       *((aaI-aaJ)**2 - DIST2)**2
     *       *((aaI**2 + 4*aaI*aaJ + aaJ**2)-DIST2)
     *       *t3EPSRR

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        C1=0.D0
       ELSE

        C1 = -3.D0/(256.D0*(aaI**3)*(aaJ**3)*DIST**3)
     *       *((aaI-aaJ)**2 - DIST2)**2
     *       *((aaI**2 + 4*aaI*aaJ + aaJ**2)-DIST2)
     *       *t3EPSRR

       ENDIF       
      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_RD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ           ! radius
      INTEGER I

      CALL YAMAKAWA_RD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END

C**********************************************************
C**********************************************************
C**********************************************************


      SUBROUTINE MOB_TO_FRI(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN)

      CALL MATREV(APP,6*NN)

      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ + MATMUL(AQP,APQ)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE FRI_TO_MOB(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN)

      CALL MATREV(APP,6*NN)
      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ - MATMUL(AQP,APQ)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE FRI_TO_MOB_RED(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN)

      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ - MATMUL(AQP,APQ)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE INVFRI_TO_FRI(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN),AP(11*NN,11*NN)

      AP(11*NN,11*NN) = 0.D0

      AP(1:6*NN,1:6*NN) = APP
      AP(1:6*NN,6*NN+1:11*NN) = APQ 
      AP(6*NN+1:11*NN,1:6*NN) = TRANSPOSE(APQ)
      AP(6*NN+1:11*NN,6*NN+1:11*NN) = AQQ 

      CALL MATREV(APP,6*NN)
      CALL MATREV(AP,11*NN)

      APP = AP(1:6*NN,1:6*NN)
      APQ = AP(1:6*NN,6*NN+1:11*NN)
      AQQ = AP(6*NN+1:11*NN,6*NN+1:11*NN)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE APQtoAQP(TPQ,TQP,NN)
      IMPLICIT NONE
      REAL*8 TPQ(6*NN,5*NN),TQP(5*NN,6*NN)
      INTEGER I,J,NN

      TQP = 0.D0

      DO I=1,6*NN
       DO J=1,5*NN
        TQP(J,I) = -TPQ(I,J)
       ENDDO
      ENDDO

      RETURN
      END


 
