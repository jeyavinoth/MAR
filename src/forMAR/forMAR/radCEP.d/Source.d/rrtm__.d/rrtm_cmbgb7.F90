!***************************************************************************
SUBROUTINE RRTM_CMBGB7
!***************************************************************************

!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE PARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT 

USE YOERRTO7 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO, ABSCO2O
USE YOERRTA7 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
             &        FRACREFB,ABSCO2                          ,& 
             &        ABSA    ,ABSB    ,NG7
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JN, JP, JT
INTEGER_M :: MEQ, NEQ                    ! To force equivalence, HG, 13-DEC-2003

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK


DO JN = 1,9
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(7)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(6)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+96)
        ENDDO

        KA(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO
DO JT = 1,5
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(7)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(6)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+96)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(7)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(6)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+96)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,9
  IPRSM = 0
  DO IGC = 1,NGC(7)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(6)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(7)
  SUMF = _ZERO_
  SUMK = _ZERO_
  DO IPR = 1, NGN(NGS(6)+IGC)
    IPRSM = IPRSM + 1


    SUMF = SUMF + FRACREFBO(IPRSM)
    SUMK = SUMK + ABSCO2O(IPRSM)*RWGT(IPRSM+96)
  ENDDO


  FRACREFB(IGC) = SUMF
  ABSCO2(IGC) = SUMK
ENDDO

DO JP = 1,9
  DO IGC = 1,NGC(7)

    FREFA(NGS(6)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,8
  DO IGC = 1,NGC(7)


    FREFADF(NGS(6)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO IGC = 1,NGC(7)

  FREFB(NGS(6)+IGC,1) = FRACREFB(IGC)
ENDDO

! +--Force the equivalence: BEGIN (HG, 13-DEC-2003)
! +  ============================

! +--ABSA
! +  ^^^^
         JN  = 0
         JT  = 1
         JP  = 1
         IGC = 1
      DO NEQ=1,NG7
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

! +--ABSB
! +  ^^^^
         JN  = 0
         JP  = 13
         IGC = 1
      DO NEQ=1,NG7
      DO MEQ=1,235
             JN =  JN  + 1
      IF   ( JN == 5   + 1)                                         THEN
             JN =  1
             JP =  JP  + 1
        IF ( JP == 59  + 1 )                                        THEN
             JP =  13
             IGC=  IGC + 1
        END IF
      END IF
             ABSB(MEQ,NEQ) = KB(JN,JP,IGC)
      ENDDO
      ENDDO

! +--Force the equivalence: END   (HG, 13-DEC-2003)
! +  ==========================


RETURN
END SUBROUTINE RRTM_CMBGB7
