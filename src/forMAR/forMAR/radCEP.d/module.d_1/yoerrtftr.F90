MODULE YOERRTFTR

#include "tsmbkind.h"

USE PARRRTM


IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

INTEGER_M :: NGC(JPBAND)
INTEGER_M :: NGS(JPBAND)
INTEGER_M :: NGN(JPGPT)
INTEGER_M :: NGB(JPGPT)

INTEGER_M :: NGM(JPG*JPBAND)
REAL_B ::    WT(JPG)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
!  NGC   : INTEGER :
!  NGS   : INTEGER :
!  NGN   : INTEGER :
!  NGB   : INTEGER :
!  NGM   : INTEGER :
!  WT    : REAL    :
!    -------------------------------------------------------------------
END MODULE YOERRTFTR


