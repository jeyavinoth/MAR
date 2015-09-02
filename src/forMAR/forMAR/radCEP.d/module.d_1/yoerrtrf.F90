MODULE YOERRTRF


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTRF* - RRTM REFERENCE ATMOSPHERE
!     -----------------------------------------------------------------

REAL_B , DIMENSION(59) :: PREF
REAL_B , DIMENSION(59) :: PREFLOG
REAL_B , DIMENSION(59) :: TREF

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! PREF   :  REAL    
! PREFLOG: REAL
! TREF   : REAL
!     -----------------------------------------------------------------
END MODULE YOERRTRF
