      SUBROUTINE STEPPER_NOBROWN(CONF,RADII,DT,T)    ! RETURNS NEW T=T+DT
      USE SIZE           ! NN
      USE LATTICE_SKEW   ! ERR
      USE FORCE_PAR      ! GAMMA
      IMPLICIT NONE
      REAL*8 CONF(3*NN),RADII(NN),DT,T
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 A(3*NN,3*NN),C(3,NN),LATTICE(3,3),EWS_ERR
      REAL*8 CONF_1(3*NN),DCONF_1(3*NN),DCONF_2(3*NN)
      REAL*8 F(3*NN),VM(3*NN),VS(3*NN),VC(3*NN),VB(3*NN),G(3*NN)
      INTEGER I,J
                  
C EVALUATING DCONF_1 and CONF_1 starting in CONF ----
      
      CALL SET_LATTICE_SKEW(T,LATTICE)
      EWS_ERR = LINERR

      CALL GRPERY_MOB(APP,APQ,AQQ,CONF,RADII,NN,LATTICE,EWS_ERR)

      A=0.D0
      DO I=1,NN
       DO J=1,NN
        A(3*(I-1)+1:3*I,3*(J-1)+1:3*J) = 
     * APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+1:6*(J-1)+3)
       ENDDO
      ENDDO

      C=0.D0
      DO I=1,NN
       DO J=1,NN
        C(1:3,I) = C(1:3,I) -
     * APQ(6*(I-1)+1:6*(I-1)+3,5*(J-1)+4)/SQRT(2.D0)
       ENDDO
      ENDDO

      CALL FORCE_PER(F,CONF)

      VM=MATMUL(A,F)    ! MOBILITY*FORCE

      VS=0.D0
      VC=0.D0
      DO I=1,NN
       VS(3*(I-1)+1)=CONF(3*(I-1)+3)        ! CONVECTION FAXEN
       VC(3*(I-1)+1:3*(I-1)+3) = C(1:3,I)   ! CONVECTION PARTICLE PAIRWISE INTERACTIONS
      ENDDO
      VS=GAMMA*VS
      VC=GAMMA*VC

      DCONF_1=DT*(VM+VS+VC) 

      CONF_1=CONF+DCONF_1

C EVALUATING DCONF_2 starting in CONF_1 ----

      T=T+DT
      CALL SET_LATTICE_SKEW(T,LATTICE)

      CALL GRPERY_MOB(APP,APQ,AQQ,CONF_1,RADII,NN,LATTICE,EWS_ERR)

      A=0.D0
      DO I=1,NN
       DO J=1,NN
        A(3*(I-1)+1:3*I,3*(J-1)+1:3*J) = 
     * APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+1:6*(J-1)+3)
       ENDDO
      ENDDO

      C=0.D0
      DO I=1,NN
       DO J=1,NN
        C(1:3,I) = C(1:3,I) -
     * APQ(6*(I-1)+1:6*(I-1)+3,5*(J-1)+4)/SQRT(2.D0)
       ENDDO
      ENDDO

      CALL FORCE_PER(F,CONF_1)

      VM=MATMUL(A,F)    ! MOBILITY*FORCE

      VS=0.D0
      VC=0.D0
      DO I=1,NN
       VS(3*(I-1)+1)=CONF_1(3*(I-1)+3)        ! CONVECTION FAXEN
       VC(3*(I-1)+1:3*(I-1)+3) = C(1:3,I)     ! CONVECTION PARTICLE PAIRWISE INTERACTIONS
      ENDDO
      VS=GAMMA*VS
      VC=GAMMA*VC

      DCONF_2=DT*(VM+VS+VC)

C EVALUATING THE FINAL CONF -------------

      CONF=CONF + 0.5D0*(DCONF_1+DCONF_2)

      RETURN
      END

