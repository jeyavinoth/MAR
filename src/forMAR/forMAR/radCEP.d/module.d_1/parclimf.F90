MODULE PARCLIMF

#include "tsmbkind.h"

USE PARCLI


IMPLICIT NONE

SAVE

!------------------------------------------------------------------------------
!     JPNCHS : Number of *N108* fields (depending on JPNMOI).
!     JPNCHF : Number of final  fields (depending on JPNMOI).

INTEGER_M, PARAMETER :: JPNCHS=JPNMOS*JPNMOI+JPNFIS
INTEGER_M, PARAMETER :: JPNCHF=JPNMOF*JPNMOI+JPNFIF

END MODULE PARCLIMF
