C $Header: /u/gcmpack/MITgcm/pkg/iceberg/iceberg_init_fixed.F,v 1.5 2011/12/28 20:45:05 jmc Exp $
C $Name: checkpoint65m $

#include "ICEBERG_OPTIONS.h"

      SUBROUTINE ICEBERG_INIT_FIXED( myThid )
C     *============================================================*
C     | SUBROUTINE ICEBERG_INIT_FIXED
C     | o Routine to initialize ICEBERG parameters and variables.
C     *============================================================*
C     | Initialize ICEBERG parameters and variables.
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEBERG.h"

C     === Routine arguments ===
C     myThid -  Number of this instance of ICEBERG_INIT_FIXED
      INTEGER myThid

#ifdef ALLOW_ICEBERG
C     === Local variables ===
C     I,J,K,bi,bj - Loop counters
      INTEGER I, J, K, bi, bj

C      Read iceberg mask
       IF (ICEBERGmaskFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGmaskFile,icebergMask,1,0,myThid)
       ENDIF

C      Read iceberg numbers file
       IF (ICEBERGmaskNumsFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGmaskNumsFile,icebergMaskNums,1,0,myThid)
       ENDIF
	   
C      Read iceberg num per cell file
       IF (ICEBERGnumPerCellFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGnumPerCellFile,icebergNumBergs,1,0,myThid)
       ENDIF
	   
C      Read drift mask
       IF (ICEBERGdriftFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGdriftFile,driftMask,1,0,myThid)
       ENDIF
	   
C      Read barrier mask
       IF (ICEBERGbarrierFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL
     &     (ICEBERGbarrierFile,barrierMask,1,0,myThid)
       ENDIF
	   
C      Read iceberg area
       IF (ICEBERGareaFile .NE. ' ' ) THEN
        CALL READ_REC_XYZ_RL
     &     (ICEBERGareaFile,icebergArea3D,1,0,myThid)
       ENDIF
	   
C      Read open fraction file
       IF ( ICEBERGopenFracFile .NE. ' ' ) THEN
        CALL READ_REC_XYZ_RL
     &     (ICEBERGopenFracFile,openFraction,1,0,myThid)  
	 
C          Partially block cells filled by icebergs
           DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
             DO J = 1, sNy
              DO I = 1, sNx
			  
C             If we want to make icebergs physical barriers
               IF ( barrierMask(I,J,bi,bj) .GT. 0. _d 0 ) THEN
C              Loop through water column
                DO K = 1, Nr
                  hFacC(I,J,K,bi,bj) = openFraction(I,J,K,bi,bj)
                ENDDO
               ENDIF
			   
              ENDDO
             ENDDO
            ENDDO
           ENDDO
       ENDIF
	   


#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       CALL ICEBERG_DIAGNOSTICS_INIT(myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_ICEBERG */

      RETURN
      END
