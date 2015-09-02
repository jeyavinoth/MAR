SUBROUTINE RADOZC ( KIDIA , KFDIA , KLON , KTDIA , KLEV &
                  &, KRINT , KDLON , KSHIFT &
                  &, PAPRS , PGEMU &
                  &, POZON                 )

!***********************************************************************
! CAUTION: THIS ROUTINE WORKS ONLY ON A NON-ROTATED, UNSTRETCHED GRID
!***********************************************************************

!**** *RADOZC* - COMPUTES DISTRIBUTION OF OZONE FROM CLIMATOLOGY

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *RADOZC* FROM *RADINT*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     J.-J. MORCRETTE  E.C.M.W.F.    95/01/25

!     MODIFICATIONS.
!     --------------
!          D.SALMOND   ECMWF    99-06-14   Optimisation

!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOMCST   , ONLY : R
USE YOEOZOC  , ONLY : RSINC    ,ROZT     ,RPROC


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KDLON
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KRINT
INTEGER_M :: KSHIFT
INTEGER_M :: KTDIA




!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

REAL_B :: PAPRS(KLON,KLEV+1),PGEMU(KLON)

REAL_B :: POZON(KDLON,KLEV)

!     ----------------------------------------------------------------- 

!*       0.2   LOCAL ARRAYS.
!              -------------

REAL_B :: ZOZLT(0:35) , ZOZON(KDLON,KLEV+1)
REAL_B :: ZRRR(0:34)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IL, INLA, JC, JK, JL

!     LOCAL REAL SCALARS
REAL_B :: ZPMR, ZSILAT, ZSIN

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

!*         1.     LATITUDE INDEX WITHIN OZONE CLIMATOLOGY
!                 ---------------------------------------



ZSIN=PGEMU(KIDIA)
INLA=0
ZSILAT=-9999._JPRB
DO JL=18,1,-1
  IF (ZSIN <= RSINC(JL+1).AND.ZSIN >= RSINC(JL)) THEN
    INLA=JL
  ENDIF
ENDDO
IF (INLA == 0) THEN
!  CALL ABOR1(' Problem with lat. interpolation in radozc!')
  PRINT *,' Problem with lat. interpolation in radozc!'
ENDIF
ZSILAT=(ZSIN-RSINC(INLA))/(RSINC(INLA+1)-RSINC(INLA))

!     ------------------------------------------------------------------


!*         2.     LATITUDE INTERPOLATED FIELD
!                 ---------------------------


IF(INLA == 18.OR.INLA == 1) THEN
  DO JC=0,35
    ZOZLT(JC)=ROZT(INLA,JC)
  ENDDO
ELSE
  DO JC=0,35
    ZOZLT(JC)=ROZT(INLA,JC)+ZSILAT*(ROZT(INLA+1,JC)-ROZT(INLA,JC))
  ENDDO
ENDIF


!     ------------------------------------------------------------------

!*         3.     VERTICAL INTERPOLATION 
!                 ----------------------

DO JC=0,34
  ZRRR(JC)=(_ONE_/(RPROC(JC)-RPROC(JC+1)))*(ZOZLT(JC)-ZOZLT(JC+1))
ENDDO


IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  DO JC=0,34
    DO JK=1,KLEV+1
    ZPMR=PAPRS(JL,JK)
      IF(ZPMR >= RPROC(JC).AND.ZPMR < RPROC(JC+1)) &
&        ZOZON(IL,JK)=ZOZLT(JC+1)+(ZPMR-RPROC(JC+1))*ZRRR(JC)
    ENDDO
  ENDDO
ENDDO

IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  DO JK=1,KLEV+1
    ZPMR=PAPRS(JL,JK)
    ZPMR=PAPRS(JL,JK)
    IF(ZPMR >= RPROC(35)) ZOZON(IL,JK)=ZOZLT(35)
  ENDDO
ENDDO

! INTEGRATION IN THE VERTICAL:
IL=KSHIFT
DO JL=KIDIA,KFDIA,KRINT
  IL=IL+1
  DO JK=1,KLEV
    POZON(IL,JK)=(PAPRS(JL,JK+1)-PAPRS(JL,JK))&
     &*(ZOZON(IL,JK)+ZOZON(IL,JK+1))*_HALF_
  ENDDO
ENDDO

!     -----------------------------------------------------------

RETURN
END SUBROUTINE RADOZC
