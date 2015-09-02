!***************************************************************************
SUBROUTINE RRTM_CMBGB6
!***************************************************************************

!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE PARRRTM  , ONLY : JPBAND  ,JPG        ,JPXSEC     ,JPGPT 

USE YOERRTO6 , ONLY : KAO     ,SELFREFO   ,FRACREFAO  ,&
           &ABSCO2O ,CFC11ADJO,CFC12O
USE YOERRTA6 , ONLY : KA      ,SELFREF    ,FRACREFA   ,&
             &        ABSCO2  ,CFC11ADJ   ,CFC12      ,&
             &        ABSA    ,NG6
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JN, JP, JT
INTEGER_M :: MEQ, NEQ                    ! To force equivalence, HG, 13-DEC-2003

!     LOCAL REAL SCALARS
REAL_B :: SUMF, SUMK, SUMK1, SUMK2, SUMK3


DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(6)
      SUMK = _ZERO_
      DO IPR = 1, NGN(NGS(5)+IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM+80)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(6)
    SUMK = _ZERO_
    DO IPR = 1, NGN(NGS(5)+IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM+80)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(6)
  SUMF = _ZERO_
  SUMK1= _ZERO_
  SUMK2= _ZERO_
  SUMK3= _ZERO_
  DO IPR = 1, NGN(NGS(5)+IGC)
    IPRSM = IPRSM + 1




    SUMF = SUMF + FRACREFAO(IPRSM)
    SUMK1= SUMK1+ ABSCO2O(IPRSM)*RWGT(IPRSM+80)
    SUMK2= SUMK2+ CFC11ADJO(IPRSM)*RWGT(IPRSM+80)
    SUMK3= SUMK3+ CFC12O(IPRSM)*RWGT(IPRSM+80)
  ENDDO




  FRACREFA(IGC) = SUMF
  ABSCO2(IGC) = SUMK1
  CFC11ADJ(IGC) = SUMK2
  CFC12(IGC) = SUMK3
ENDDO

DO IGC = 1,NGC(6)


  FREFA(NGS(5)+IGC,1) = FRACREFA(IGC)
  FREFB(NGS(5)+IGC,1) = FRACREFA(IGC)
ENDDO


! +--Force the equivalence: BEGIN (HG, 13-DEC-2003)
! +  ============================

! +--ABSA
! +  ^^^^
         JN  = 0
         JT  = 1
         JP  = 1
      DO NEQ=1,NG6
      DO MEQ=1,65
             JN =  JN  + 1
      IF   ( JN == 5   + 1)                                         THEN
             JN =  1
             JT =  JT  + 1
       IF  ( JT == 13  + 1 )                                        THEN
             JT =  1
             JP =  JP  + 1
       END IF
      END IF
             ABSA(MEQ,NEQ) = KA(JN,JT,JP)
      ENDDO
      ENDDO

! +--Force the equivalence: END   (HG, 13-DEC-2003)
! +  ==========================


RETURN
END SUBROUTINE RRTM_CMBGB6
