C WRITES THE CURRENT CONFIGURATION ON FILE

      PROGRAM MAIN
      USE SIZE           ! NN 
      USE CONFIG         ! CONF,RADII
      USE LATTICE_BASE   ! LB
      USE LATTICE_SKEW   ! LB
      USE FORCE_PAR      ! GAMMA,SIGMA
      INTEGER FROM_FILE, J, SEEDSIZE 
      INTEGER,ALLOCATABLE :: SEED(:) 
      REAL*8 T, DT, DR(3), END_TIME, BASE_DT, INP
      LOGICAL BROWN, THERE

      CALL INIT_U
      CALL INIT_EPS
      CALL INIT_Y2
      CALL INIT_II

      INQUIRE( FILE='control_file.dat', EXIST=THERE ) 
       IF ( THERE ) THEN
       OPEN(10,FILE='control_file.dat')
       READ(10,*) GAMMA
       READ(10,*) LINERR 
       READ(10,*) BROWN
       IF ( BROWN ) THEN
        CALL RANDOM_SEED(SIZE=SEEDSIZE)
        ALLOCATE(SEED(SEEDSIZE))
        READ(10,*) SEED
        CALL RANDOM_SEED(PUT=SEED)
       ENDIF
       READ(10,*) END_TIME
       READ(10,*) BASE_DT
       READ(10,*) K_WRITE
       READ(10,*) LJ_EPS
       READ(10,*) LJ_SIGMA
       READ(10,*) LJ_CUT
       CLOSE(10)
      ELSE
       WRITE(0,*) "error: no control_file.dat file "
     *   // "present in this folder"
       WRITE(0,*) "        Please provide one"
       CALL EXIT()
      ENDIF
      
      DT = BASE_DT
      K=0

      CALL INITIAL_CONF(T)

      WRITE(*,*) 'N: ',NN,' dt: ',BASE_DT,' shear: ',GAMMA,
     + ' X: ',LB(1),' Y: ',LB(2),
     + ' Z: ',LB(3)

      CALL INIT_VMD_WRITE(RADII,DT)

      DO
       DO I = 1, NN
        CALL PERIOD_CORRECTION(CONF(:,I),LB,GAMMA,T)
       ENDDO

       IF (mod(K,K_WRITE) .eq. 0) then
        CALL VMD_WRITE(CONF,T,K)
       ENDIF

       IF ( BROWN ) THEN
        CALL STEPPER_BROWN(CONF,RADII,DT,T)
       ELSE
        CALL STEPPER_NOBROWN(CONF,RADII,DT,T)
       ENDIF

       IF (T.GT.END_TIME) THEN
        EXIT
       ENDIF
          
       K=K+1

      ENDDO

      CALL END_VMD_WRITE()

      CLOSE(31)

      STOP
      END


C***********************************************************
C***********************************************************
C***********************************************************

C VTF file for vmd reader header

      SUBROUTINE INIT_VMD_WRITE(RADII,DT)
      USE SIZE
      USE LATTICE_BASE
      USE FORCE_PAR
      IMPLICIT NONE
      REAL*8  RADII(NN),DT
      INTEGER I

      OPEN(31,FILE='output.VTF') ! xyz file for vmd

      WRITE(31,*) '# ',NN*1,' ',L0_PAR,' ',K_PAR,' ',LJ_SIGMA,
     *           ' ',LJ_EPS
      WRITE(31,*) '# ',LB(1),' ',LB(2),' ',LB(3),' ',DT,' ',GAMMA
      DO I=1,NN
       WRITE(31,*) 'atom ',I-1,'radius',RADII(I),'name H'
      ENDDO

      WRITE(31,*) ''
      WRITE(31,*) 'pbc',LB(1),LB(2),LB(3)
      WRITE(31,*)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

C VTF file for vmd reader timestep write

      SUBROUTINE VMD_WRITE(CONF,T,K)
      USE SIZE
      IMPLICIT NONE
      REAL(8) CONF(3,NN),T
      INTEGER I,K

      WRITE(31,*) ""
      WRITE(31,*) "timestep"
      WRITE(31,*) "#",T,K
      DO I = 1, NN
       WRITE(31,*) CONF(:,I)
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

C close VTF file for vmd reader

      SUBROUTINE END_VMD_WRITE()
      IMPLICIT NONE

      CLOSE(31) ! xyz file for vmd

      RETURN
      END

      
