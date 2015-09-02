!***************************************************************************
SUBROUTINE RRTM_CMBGB1
!***************************************************************************

!  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 16 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
!  in arrays FRACREFA and FRACREFB are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTM.

!  BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!***************************************************************************
!  INSTRUCTION EQUIVALENCE SUPPRESSED (H. Gallée, LGGE, 15 décembre 2003)
!***************************************************************************

! Parameters
#include "tsmbkind.h"

USE PARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT 

USE YOERRTO1 , ONLY : KAO , KBO  , SELFREFO, FORREFO, FRACREFAO,FRACREFBO
USE YOERRTA1 , ONLY : KA  , KB   , SELFREF , FORREF , FRACREFA ,FRACREFB      &
             &      , ABSA, ABSB , NG1
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT

IMPLICIT NONE

!     LOCAL INTEGER SCALARS
INTEGER_M :: IGC, IPR, IPRSM, JP, JT
INTEGER_M :: MEQ, NEQ                    ! To force equivalence, HG, 13-DEC-2003

!     LOCAL REAL SCALARS
REAL_B :: SUMF1, SUMF2, SUMK


DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(1)
      SUMK = _ZERO_
      DO IPR = 1, NGN(IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KAO(JT,JP,IPRSM)*RWGT(IPRSM)
      ENDDO

      KA(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(1)
      SUMK = _ZERO_
      DO IPR = 1, NGN(IGC)
        IPRSM = IPRSM + 1

        SUMK = SUMK + KBO(JT,JP,IPRSM)*RWGT(IPRSM)
      ENDDO

      KB(JT,JP,IGC) = SUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(1)
    SUMK = _ZERO_
    DO IPR = 1, NGN(IGC)
      IPRSM = IPRSM + 1

      SUMK = SUMK + SELFREFO(JT,IPRSM)*RWGT(IPRSM)
    ENDDO

    SELFREF(JT,IGC) = SUMK
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(1)
  SUMK = _ZERO_
  SUMF1 = _ZERO_
  SUMF2 = _ZERO_
  DO IPR = 1, NGN(IGC)
    IPRSM = IPRSM + 1



    SUMK = SUMK + FORREFO(IPRSM)*RWGT(IPRSM)
    SUMF1= SUMF1+ FRACREFAO(IPRSM)
    SUMF2= SUMF2+ FRACREFBO(IPRSM)
  ENDDO



  FORREF(IGC) = SUMK
  FRACREFA(IGC) = SUMF1
  FRACREFB(IGC) = SUMF2
ENDDO

DO IGC = 1,NGC(1)


  FREFA(IGC,1) = FRACREFA(IGC)
  FREFB(IGC,1) = FRACREFB(IGC)
ENDDO


! +--Force the equivalence: BEGIN (HG, 13-DEC-2003)
! +  ============================

! +--ABSA
! +  ^^^^
         JT  = 0
         JP  = 1
         IGC = 1
      DO NEQ=1,NG1
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
      DO NEQ=1,NG1
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
END SUBROUTINE RRTM_CMBGB1
