#include "ICEBERG_OPTIONS.h"

CBOP
C !ROUTINE: ICEBERG_V_DRAG_COEFF

C !INTERFACE: ==========================================================
      SUBROUTINE ICEBERG_V_DRAG_COEFF(
     I        bi, bj, k, inp_KE,
     I        uFld, vFld, kappaRV,
     U        KE,
     O        cDrag,
     I        myIter, myThid )

C !DESCRIPTION:
C Calculates the drag coefficient due to friction and the no-slip condition
C at the bottom of the iceberg (in analogy to bottom drag)
C such as the iceberg stress: tauy_{ice} = -Cd * V_{top} * rUnit2mass ;

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "ICEBERG_SIZE.h"
#include "ICEBERG.h"

C !INPUT PARAMETERS: ===================================================
C  bi,bj          :: tile indices
C  k              :: vertical level to process
C  inp_KE         :: =T : KE is provided as input ; =F : to compute here
C  uFld           :: velocity, zonal component
C  vFld           :: velocity, meridional component
C  kappaRV        :: vertical viscosity
C  KE             :: Kinetic energy (input when inp_KE = T)
C  myIter         :: current iteration number
C  myThid         :: my Thread Id number
      INTEGER bi, bj, k
      LOGICAL inp_KE
      _RL uFld   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL vFld   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL KE     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL kappaRV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr+1)
      INTEGER myIter, myThid

C !OUTPUT PARAMETERS: ==================================================
C  KE             :: Kinetic energy (output when inp_KE = F)
C  cDrag          :: bottom drag coefficient
      _RL cDrag  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)

#ifdef ALLOW_ICEBERG
C !LOCAL VARIABLES : ====================================================
C  i,j                  :: loop indices
C  Kp1                  :: =k+1 for k<Nr, =Nr for k>=Nr
      INTEGER i, j
      INTEGER kUpC, kTop
      _RL viscFac, dragFac, vSq
      _RL rDrCkp1
CEOP

C-  No-slip BCs impose a drag at top
      IF ( usingZCoords ) THEN
       kTop    = 1
       kUpC    = k
c      dragFac = mass2rUnit*rhoConst
       dragFac = 1. _d 0
      ELSE
       kTop    = Nr
       kUpC    = k+1
       dragFac = mass2rUnit*rhoConst
      ENDIF
      rDrCkp1 = recip_drC(kUpC)
CML      IF (k.EQ.kTop) rDrCkp1=recip_drF(k)
      viscFac = 0. _d 0
      IF (brg_NoSlip) viscFac = 2. _d 0

C--   Initialise drag-coeff
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
         cDrag(i,j) = 0. _d 0
       ENDDO
      ENDDO

C--   Friction at the bottom of iceberg (no-slip BC)
      IF ( brg_NoSlip ) THEN
C-    friction accounts for true distance (including hFac) to the surface
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx,sNx+OLx-1
         IF ( k.LE.
     &       MAX( kBergBotC(i,j-1,bi,bj),kBergBotC(i,j,bi,bj) ) ) THEN
           cDrag(i,j) = cDrag(i,j)
     &                + kappaRV(i,j,kUpC)*rDrCkp1*viscFac
     &                        * _recip_hFacS(i,j,k,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

C--   Add Linear drag contribution:
      IF ( brg_DragLinear.NE.zeroRL ) THEN
       DO j=1-OLy+1,sNy+OLy
        DO i=1-OLx,sNx+OLx
         IF ( k.LE.
     &       MAX( kBergBotC(i,j-1,bi,bj),kBergBotC(i,j,bi,bj) ) ) THEN
          cDrag(i,j) = cDrag(i,j) + brg_DragLinear*dragFac
         ENDIF
        ENDDO
       ENDDO
      ENDIF

C--   Add quadratic drag
      IF ( brg_SelectDragQuad.EQ.0 ) THEN
        IF ( .NOT.inp_KE ) THEN
          DO j=1-OLy,sNy+OLy-1
           DO i=1-OLx,sNx+OLx-1
            KE(i,j) = 0.25*(
     &          ( uFld( i , j )*uFld( i , j )*_hFacW(i,j,k,bi,bj)
     &           +uFld(i+1, j )*uFld(i+1, j )*_hFacW(i+1,j,k,bi,bj) )
     &        + ( vFld( i , j )*vFld( i , j )*_hFacS(i,j,k,bi,bj)
     &           +vFld( i ,j+1)*vFld( i ,j+1)*_hFacS(i,j+1,k,bi,bj) )
     &                     )*_recip_hFacC(i,j,k,bi,bj)
           ENDDO
          ENDDO
        ENDIF
C-    average grid-cell-center KE to get velocity norm @ U.pt
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx,sNx+OLx-1
          vSq = 0. _d 0
          IF ( k.LE.
     &       MAX( kBergBotC(i,j-1,bi,bj),kBergBotC(i,j,bi,bj) ) ) THEN
           vSq = KE(i,j)+KE(i,j-1)
          ENDIF
          IF ( vSq.GT.zeroRL ) THEN
           cDrag(i,j) = cDrag(i,j)  !only bottom drag
     &                + brgDragQuadFld(i,j,bi,bj)*SQRT(vSq)*dragFac
     &            + brg_DragForm * dragFac * SQRT(vSq) * 
     &                   recip_dyF(i,j,bi,bj) * drF(k) *
     &                   (1 - _hFacC(i,j,k,bi,bj**2)
          ENDIF
         ENDDO
        ENDDO
      ELSEIF ( brg_SelectDragQuad.EQ.1 ) THEN
C-    calculate locally velocity norm @ U.pt (local U & 4 V averaged)
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx,sNx+OLx-1
          vSq = 0. _d 0
          IF ( k.LE.
     &       MAX( kBergBotC(i,j-1,bi,bj),kBergBotC(i,j,bi,bj) ) ) THEN
           vSq = vFld(i,j)*vFld(i,j)
     &       + ( (uFld( i ,j-1)*uFld( i ,j-1)*hFacW( i ,j-1,k,bi,bj)
     &           +uFld( i , j )*uFld( i , j )*hFacW( i , j ,k,bi,bj))
     &         + (uFld(i+1,j-1)*uFld(i+1,j-1)*hFacW(i+1,j-1,k,bi,bj)
     &           +uFld(i+1, j )*uFld(i+1, j )*hFacW(i+1, j ,k,bi,bj))
     &         )*recip_hFacS(i,j,k,bi,bj)*0.25 _d 0
          ENDIF
          IF ( vSq.GT.zeroRL ) THEN
           cDrag(i,j) = cDrag(i,j)
     &                + brgDragQuadFld(i,j,bi,bj)*SQRT(vSq)*dragFac
     &            + brg_DragForm * dragFac * SQRT(vSq) * 
     &                   recip_dyF(i,j,bi,bj) * drF(k) *
     &                   (1 - _hFacC(i,j,k,bi,bj**2)
          ENDIF
        ENDDO
       ENDDO
      ELSEIF ( brg_SelectDragQuad.EQ.2 ) THEN
C-    same as above but using wet-point method to average 4 V
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx,sNx+OLx-1
          vSq = 0. _d 0
          IF ( k.LE.
     &       MAX( kBergBotC(i,j-1,bi,bj),kBergBotC(i,j,bi,bj) ) ) THEN
           vSq = ( hFacW( i ,j-1,k,bi,bj) + hFacW( i , j ,k,bi,bj) )
     &         + ( hFacW(i+1,j-1,k,bi,bj) + hFacW(i+1, j ,k,bi,bj) )
           IF ( vSq.GT.zeroRL ) THEN
            vSq = vFld(i,j)*vFld(i,j)
     &        +( (uFld( i ,j-1)*uFld( i ,j-1)*hFacW( i ,j-1,k,bi,bj)
     &           +uFld( i , j )*uFld( i , j )*hFacW( i , j ,k,bi,bj))
     &         + (uFld(i+1,j-1)*uFld(i+1,j-1)*hFacW(i+1,j-1,k,bi,bj)
     &           +uFld(i+1, j )*uFld(i+1, j )*hFacW(i+1, j ,k,bi,bj))
     &         )/vSq
           ELSE
            vSq = vFld(i,j)*vFld(i,j)
           ENDIF
          ENDIF
          IF ( vSq.GT.zeroRL ) THEN
           cDrag(i,j) = cDrag(i,j)
     &                + brgDragQuadFld(i,j,bi,bj)*SQRT(vSq)*dragFac
     &            + brg_DragForm * dragFac * SQRT(vSq) * 
     &                   recip_dyF(i,j,bi,bj) * drF(k) *
     &                   (1 - _hFacC(i,j,k,bi,bj**2)
          ENDIF
        ENDDO
       ENDDO
      ENDIF

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       IF ( selectImplicitDrag.EQ.0 ) THEN !this is default
C-     Explicit case: diagnose directly the Iceberg stress, only implemented case
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           icebergDragV(i,j,k,bi,bj) = icebergDragV(i,j,k,bi,bj)
     &                              - cDrag(i,j)*vFld(i,j)*rUnit2mass
          ENDDO
         ENDDO
       ENDIF
       CALL DIAGNOSTICS_FILL_RS( icebergDragV, 'BRG_TauY',
     &                          0, Nr, 1, bi, bj, myThid )
      ENDIF !useDiagnostics+ 
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_ICEBERG */

      RETURN
      END
