!***************************************************************************
SUBROUTINE RRTM_CMBGB10
!***************************************************************************

!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE PARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT 

USE YOERRTO10, ONLY : KAO     ,KBO      ,FRACREFAO   ,FRACREFBO
USE YOERRTA10, ONLY : KA      ,KB       ,FRACREFA    ,FRACREFB             &
             &      , ABSA    ,ABSB     ,NG10
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JP, JT
INTEGER_M :: MEQ, NEQ                ! To force equivalence, HG, 13-DEC-2003

!     LOCAL REAL SCALARS
REAL_B :: SUMF1, SUMF2, SUMK


DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(10)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(9)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+144)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(10)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(9)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+144)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(10)
  SUMF1= _ZERO_
  SUMF2= _ZERO_
  DO IPR = 1, NGN(NGS(9)+IGC)
    IPRSM = IPRSM + 1


    SUMF1= SUMF1+ FRACREFAO(IPRSM)
    SUMF2= SUMF2+ FRACREFBO(IPRSM)
  ENDDO


  FRACREFA(IGC) = SUMF1
  FRACREFB(IGC) = SUMF2
ENDDO

DO IGC = 1,NGC(10)


  FREFA(NGS(9)+IGC,1) = FRACREFA(IGC)
  FREFB(NGS(9)+IGC,1) = FRACREFB(IGC)
ENDDO


! +--Force the equivalence: BEGIN (HG, 13-DEC-2003)
! +  ============================

! +--ABSA
! +  ^^^^
         JT  = 0
         JP  = 1
         IGC = 1
      DO NEQ=1,NG10
      DO MEQ=1,65
             JT =  JT  + 1
       IF  ( JT == 5   + 1 )                                        THEN
             JT =  1
             JP =  JP  + 1
        IF ( JP == 13  + 1 )                                        THEN
             JP =  1
             IGC=  IGC + 1
        END IF
       END IF
             ABSA(MEQ,NEQ) = KA(JT,JP,IGC)
      ENDDO
      ENDDO

! +--ABSB
! +  ^^^^
         JT  = 0
         JP  = 13
         IGC = 1
      DO NEQ=1,NG10
      DO MEQ=1,235
             JT =  JT  + 1
       IF  ( JT == 5   + 1 )                                        THEN
             JT =  1
             JP =  JP  + 1
        IF ( JP == 59  + 1 )                                        THEN
             JP =  13
             IGC=  IGC + 1
        END IF
       END IF
             ABSB(MEQ,NEQ) = KB(JT,JP,IGC)
      ENDDO
      ENDDO

! +--Force the equivalence: END   (HG, 13-DEC-2003)
! +  ==========================


RETURN
END SUBROUTINE RRTM_CMBGB10
