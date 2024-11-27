#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
CBOP
C !ROUTINE: ICEBERG_DIAGNOSTICS_DRAG

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_DIAGNOSTICS_DRAG(
     I        uFld, vFld,
     I        bi, bj, myIter, myThid )

C !DESCRIPTION:
C Finish to calculate iceberg stress and fill-up corresponding diagnostics
C CURRENTLY NOT USED
C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"

C !INPUT PARAMETERS: ===================================================
C  uFld           :: velocity, zonal component
C  vFld           :: velocity, meridional component
C  bi,bj          :: tile indices
C  myIter         :: current iteration number
C  myThid         :: my Thread Id number
      _RL uFld   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL vFld   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER bi, bj
      INTEGER myIter, myThid

#ifdef ALLOW_ICEBERG
#ifdef ALLOW_DIAGNOSTICS
C !LOCAL VARIABLES : ====================================================
C  i,j            :: loop indices
C  k              :: vertical level
      INTEGER i, j, k
CEOP

      IF ( selectImplicitDrag.EQ.2 ) THEN
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx+1,sNx+OLx
         DO K = 1,Nr
           icebergDragU(i,j,k,bi,bj) = -icebergDragU(i,j,k,bi,bj)
     &                              * uFld(i,j,k)
         ENDDO
        ENDDO
       ENDDO
       DO j=1-OLy+1,sNy+OLy
        DO i=1-OLx,sNx+OLx
         DO K = 1,Nr
           icebergDragV(i,j,k,bi,bj) = -icebergDragV(i,j,k,bi,bj)
     &                              * vFld(i,j,k)
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      CALL DIAGNOSTICS_FILL_RS( icebergDragU, 'BRG_TauX',
     &                          0, Nr, 1, bi, bj, myThid )
      CALL DIAGNOSTICS_FILL_RS( icebergDragV, 'BRG_TauY',
     &                          0, Nr, 1, bi, bj, myThid )

#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_ICEBERG */

      RETURN
      END
