!***************************************************************************
SUBROUTINE RRTM_CMBGB2
!***************************************************************************

!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE PARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT 

USE YOERRTO2 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO  ,FORREFO
USE YOERRTA2 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
             &        FRACREFB,FORREF  ,REFPARAM               ,&
             &        ABSA    ,ABSB    ,NG2
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JP, JT
INTEGER_M :: MEQ, NEQ                    ! To force equivalence, HG, 13-DEC-2003

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK


DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(2)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(1)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+16)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(2)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(1)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM+16)
      ENDDO
!               KBC(JT,JP,IGC) = SUMK
      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(2)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(1)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+16)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

DO JP = 1,13
  IPRSM = 0
  DO IGC = 1,NGC(2)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(1)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFAO(IPRSM,JP)
    ENDDO

    FRACREFA(IGC,JP) = SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(2)
  SUMK = _ZERO_
  SUMF = _ZERO_
  DO IPR = 1, NGN(NGS(1)+IGC)
    IPRSM = IPRSM + 1


    SUMK = SUMK + FORREFO(IPRSM)*RWGT(IPRSM+16)
    SUMF = SUMF + FRACREFBO(IPRSM)
  ENDDO


  FORREF(IGC) = SUMK
  FRACREFB(IGC) = SUMF
ENDDO

DO JP = 1,13
  DO IGC = 1,NGC(2)

    FREFA(NGS(1)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 2,13
  DO IGC = 1,NGC(2)


    FREFADF(NGS(1)+IGC,JP) = FRACREFA(IGC,JP-1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO IGC = 1,NGC(2)

  FREFB(NGS(1)+IGC,1) = FRACREFB(IGC)
ENDDO


! +--Force the equivalence: BEGIN (HG, 13-DEC-2003)
! +  ============================

! +--ABSA
! +  ^^^^
         JT  = 0
         JP  = 1
         IGC = 1
      DO NEQ=1,NG2
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
      DO NEQ=1,NG2
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
END SUBROUTINE RRTM_CMBGB2
