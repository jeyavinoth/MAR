!***************************************************************************
SUBROUTINE RRTM_CMBGB12
!***************************************************************************

!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE PARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT 

USE YOERRTO12, ONLY : KAO     ,SELFREFO   ,FRACREFAO
USE YOERRTA12, ONLY : KA      ,SELFREF    ,FRACREFA                        &
             &      , ABSA    ,NG12
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JN, JP, JT
INTEGER_M :: MEQ, NEQ                ! To force equivalence, HG, 13-DEC-2003

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK


DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(12)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(11)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+176)
        ENDDO

        KA(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(12)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(11)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+176)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(12)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(11)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

DO JP = 1,9
  DO IGC = 1,NGC(12)

    FREFA(NGS(11)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,8
  DO IGC = 1,NGC(12)


    FREFADF(NGS(11)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO


! +--Force the equivalence: BEGIN (HG, 13-DEC-2003)
! +  ============================

! +--ABSA
! +  ^^^^
         JN  = 0
         JT  = 1
         JP  = 1
         IGC = 1
      DO NEQ=1,NG12
      DO MEQ=1,585
             JN =  JN  + 1
      IF   ( JN == 9   + 1)                                         THEN
             JN =  1
             JT =  JT  + 1
       IF  ( JT == 5   + 1 )                                        THEN
             JT =  1
             JP =  JP  + 1
        IF ( JP == 13  + 1 )                                        THEN
             JP =  1
             IGC=  IGC + 1
        END IF
       END IF
      END IF
             ABSA(MEQ,NEQ) = KA(JN,JT,JP,IGC)
      ENDDO
      ENDDO

! +--Force the equivalence: END   (HG, 13-DEC-2003)
! +  ==========================


RETURN
END SUBROUTINE RRTM_CMBGB12
