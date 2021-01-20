      SUBROUTINE INITIAL_CONF(T)
      USE SIZE      ! NN
      USE LATTICE_BASE ! LB
      USE CONFIG    ! CONF,POLY_LEN
      IMPLICIT NONE
      INTEGER I
      REAL*8 T
      LOGICAL THERE

      INQUIRE( FILE='initial_config.dat', EXIST=THERE ) 
       IF ( THERE ) THEN
       OPEN(11,FILE='initial_config.dat')
       READ(11,*) NN

       ALLOCATE( CONF(3,NN) )
       ALLOCATE( RADII(NN) )

       READ(11,*) LB
       READ(11,*) T
       DO I=1,NN
        READ(11,*) CONF(1,I),CONF(2,I),CONF(3,I),RADII(I)
       ENDDO
       CLOSE(11)
      ELSE
       WRITE(0,*) "error: no initial_config.dat file "
     *   // "present in this folder"
       WRITE(0,*) "        Please provide one"
       CALL EXIT()
      ENDIF

      RETURN
      END     
 
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE PERIOD_CORRECTION(R,LB,GAMMA,T)
      IMPLICIT NONE
      REAL*8 R(3),LB(3),GAMMA,T
      INTEGER I

C   CORRECTIONS FOR PERIOD

      IF(R(3).GT.LB(3)) THEN
       R(3) = R(3) - LB(3)
       R(1) = MOD(R(1) - LB(3)*GAMMA*T,LB(1))
      ENDIF
      IF(R(3).LT.0) THEN
       R(3) = R(3) + LB(3)
       R(1) = MOD(R(1) + LB(3)*GAMMA*T,LB(1))
      ENDIF

      IF(R(1).GT.LB(1)) THEN
       R(1)=R(1)-LB(1)
      ELSE IF(R(1).LT.0) THEN
       R(1)=R(1)+LB(1)
      ENDIF

      IF(R(2).GT.LB(2)) THEN
       R(2)=R(2)-LB(2)
      ELSE IF(R(2).LT.0) THEN
       R(2)=R(2)+LB(2)
      ENDIF

      RETURN
      END



