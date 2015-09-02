!***************************************************************************
SUBROUTINE RRTM_CMBGB3
!***************************************************************************

!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE PARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT 

USE YOERRTO3 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO  ,FORREFO    ,ABSN2OAO   ,ABSN2OBO
USE YOERRTA3 , ONLY : KA      ,KB      ,SELFREF    ,FRACREFA   ,&
             &        FRACREFB,FORREF  ,ABSN2OA    ,ABSN2OB    ,&
             &        ABSA    ,ABSB    ,NG3
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JN, JP, JT
INTEGER_M :: MEQ, NEQ                    ! To force equivalence, HG, 13-DEC-2003

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK, SUMK1, SUMK2, SUMK3


DO JN = 1,10
  DO JT = 1,5
    DO JP = 1,13
      IPRSM = 0
      DO IGC = 1,NGC(3)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(2)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KAO(JN,JT,JP,IPRSM)*RWGT(IPRSM+32)
        ENDDO

        KA(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO
DO JN = 1,5
  DO JT = 1,5
    DO JP = 13,59
      IPRSM = 0
      DO IGC = 1,NGC(3)
        SUMK = _ZERO_
        DO IPR = 1, NGN(NGS(2)+IGC)
          IPRSM = IPRSM + 1

          SUMK = SUMK + KBO(JN,JT,JP,IPRSM)*RWGT(IPRSM+32)
        ENDDO

        KB(JN,JT,JP,IGC) = SUMK
      ENDDO
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(3)
    SUMK = _ZERO_
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(2)+IGC)
      IPRSM = IPRSM + 1


      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+32)
      SUMF = SUMF + FRACREFAO(IPRSM,JT)
    ENDDO


    SELFREF(JT,IGC) = SUMK
    FRACREFA(IGC,JT) = SUMF
  ENDDO
ENDDO

DO JP = 1,5
  IPRSM = 0
  DO IGC = 1,NGC(3)
    SUMF = _ZERO_
    DO IPR = 1, NGN(NGS(2)+IGC)
      IPRSM = IPRSM + 1

      SUMF = SUMF + FRACREFBO(IPRSM,JP)
    ENDDO

    FRACREFB(IGC,JP) = SUMF
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(3)
  SUMK1= _ZERO_
  SUMK2= _ZERO_
  SUMK3= _ZERO_
  DO IPR = 1, NGN(NGS(2)+IGC)
    IPRSM = IPRSM + 1



    SUMK1= SUMK1+ FORREFO(IPRSM)*RWGT(IPRSM+32)
    SUMK2= SUMK2+ ABSN2OAO(IPRSM)*RWGT(IPRSM+32)
    SUMK3= SUMK3+ ABSN2OBO(IPRSM)*RWGT(IPRSM+32)
  ENDDO



  FORREF(IGC) = SUMK1
  ABSN2OA(IGC) = SUMK2
  ABSN2OB(IGC) = SUMK3
ENDDO

DO JP = 1,10
  DO IGC = 1,NGC(3)

    FREFA(NGS(2)+IGC,JP) = FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,9
  DO IGC = 1,NGC(3)


    FREFADF(NGS(2)+IGC,JP) = FRACREFA(IGC,JP+1) -FRACREFA(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,5
  DO IGC = 1,NGC(3)

    FREFB(NGS(2)+IGC,JP) = FRACREFB(IGC,JP)
  ENDDO
ENDDO
DO JP = 1,4
  DO IGC = 1,NGC(3)


    FREFBDF(NGS(2)+IGC,JP) = FRACREFB(IGC,JP+1) -FRACREFB(IGC,JP)
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
      DO NEQ=1,NG3
      DO MEQ=1,650
             JN =  JN  + 1
      IF   ( JN == 10  + 1)                                         THEN
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
         JT  = 1
         JP  = 13
         IGC = 1
      DO NEQ=1,NG3
      DO MEQ=1,1175
             JN =  JN  + 1
      IF   ( JN == 5   + 1)                                         THEN
             JN =  1
             JT =  JT  + 1
       IF  ( JT == 5   + 1 )                                        THEN
             JT =  1
             JP =  JP  + 1
        IF ( JP == 59  + 1 )                                        THEN
             JP =  13
             IGC=  IGC + 1
        END IF
       END IF
      END IF
             ABSB(MEQ,NEQ) = KB(JN,JT,JP,IGC)
      ENDDO
      ENDDO

! +--Force the equivalence: END   (HG, 13-DEC-2003)
! +  ==========================


RETURN
END SUBROUTINE RRTM_CMBGB3
