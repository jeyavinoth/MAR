MODULE YOERRTO4


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTO4* - RRTM ORIGINAL COEFFICIENTS FOR INTERVAL 4
!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     -----------------------------------------------------------------

INTEGER_M, PARAMETER :: NO4  = 16

REAL_B :: FRACREFAO(NO4,9)  ,FRACREFBO(NO4,6)
REAL_B :: KAO(9,5,13,NO4)
REAL_B :: KBO(6,5,13:59,NO4)
REAL_B :: SELFREFO(10,NO4)

!     ------------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! FRACREFB: REAL
! KA      : REAL     
! KB      : REAL     
! SELFREF : REAL     
!     -----------------------------------------------------------------
END MODULE YOERRTO4
