MODULE YOERIP


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Real time related variables for ECMWF radiation: updated UPDTIER

!     RIP0M  : I0 WEIGHTED BY THE DISTANCE EARTH-SUN

!     RCODECM: COSINE OF THE DECLINATION
!     RSIDECM:   SINE OF THE DECLINATION

!     RCOVSRM: COSINE OF TRUE SOLAR TIME
!     RSIVSRM:   SINE OF TRUE SOLAR TIME

REAL_B :: RIP0M
REAL_B :: RCODECM
REAL_B :: RSIDECM
REAL_B :: RCOVSRM
REAL_B :: RSIVSRM
!     ------------------------------------------------------------------
END MODULE YOERIP
