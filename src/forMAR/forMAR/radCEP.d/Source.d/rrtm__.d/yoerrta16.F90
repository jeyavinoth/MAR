MODULE YOERRTA16


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA16* - RRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NG16 = 2

REAL_B :: FRACREFA(NG16,9)

REAL_B :: KA(9,5,13,NG16) ,ABSA(585,NG16)
REAL_B :: SELFREF(10,NG16)
REAL_B :: STRRAT

!     -----------------------------------------------------------------
! EQUIVALENCE Instruction is suppressed and
! EQUIVALENCE             is FORCED in RRTM_CMBGB16   (HG, 13-DEC-2003)
! EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))                           ! (OLD)
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL     
! SELFREF : REAL     
! STRRAT  : REAL
!     -----------------------------------------------------------------
END MODULE YOERRTA16
