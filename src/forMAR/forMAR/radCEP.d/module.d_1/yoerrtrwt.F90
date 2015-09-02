MODULE YOERRTRWT

#include "tsmbkind.h"

USE PARRRTM


IMPLICIT NONE

SAVE

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------

REAL_B :: FREFA  (JPGPT,13)
REAL_B :: FREFB  (JPGPT,6)
REAL_B :: FREFADF(JPGPT,13)
REAL_B :: FREFBDF(JPGPT,6)
REAL_B :: RWGT   (JPG*JPBAND)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! FREFA  :  REAL   :
! FREFB  :  REAL   :
! FREFADF:  REAL   :
! FREFBDF:  REAL   :
! RWT    :  REAL   : 
!    -------------------------------------------------------------------
END MODULE YOERRTRWT
