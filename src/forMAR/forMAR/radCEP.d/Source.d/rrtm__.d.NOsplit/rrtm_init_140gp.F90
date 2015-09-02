!***************************************************************************
SUBROUTINE RRTM_INIT_140GP
!***************************************************************************
!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

! Parameters
#include "tsmbkind.h"

USE PARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT
USE YOERRTWN , ONLY : NG       ,NSPA     ,NSPB
USE YOERRTFTR, ONLY : NGC      ,NGS      ,NGN      ,NGB       ,NGM     , WT
! Output
USE YOERRTBG2, ONLY : CORR1    ,CORR2
USE YOERRTRWT, ONLY : FREFA    ,FREFB    ,FREFADF  ,FREFBDF   ,RWGT

! Local

IMPLICIT NONE
REAL_B :: WTSM(JPG)

!     LOCAL INTEGER SCALARS
INTEGER_M :: I, IBND, IG, IGC, IGCSM, IND, IPR, IPRSM, IPT

!     LOCAL REAL SCALARS
REAL_B :: FP, RTFP, WTSUM


!  Calculate lookup tables for functions needed in routine TAUMOL (TAUGB2)
CORR1(0) = _ONE_
CORR1(200) = _ONE_
CORR2(0) = _ONE_
CORR2(200) = _ONE_
DO I = 1,199
  FP = 0.005_JPRB*REAL(I)
  RTFP = SQRT(FP)
  CORR1(I) = RTFP/FP
  CORR2(I) = (_ONE_-RTFP)/(_ONE_-FP)
ENDDO

!  Perform g-point reduction from 16 per band (256 total points) to
!  a band dependant number (140 total points) for all absorption
!  coefficient input data and Planck fraction input data.
!  Compute relative weighting for new g-point combinations.

IGCSM = 0
DO IBND = 1,JPBAND
  IPRSM = 0
  IF (NGC(IBND) < 16) THEN
    DO IGC = 1,NGC(IBND)
      IGCSM = IGCSM + 1
      WTSUM = _ZERO_
      DO IPR = 1, NGN(IGCSM)
        IPRSM = IPRSM + 1
        WTSUM = WTSUM + WT(IPRSM)
      ENDDO
      WTSM(IGC) = WTSUM
    ENDDO
    DO IG = 1,NG(IBND)
      IND = (IBND-1)*16 + IG
      RWGT(IND) = WT(IG)/WTSM(NGM(IND))
    ENDDO
  ELSE
    DO IG = 1,NG(IBND)
      IGCSM = IGCSM + 1
      IND = (IBND-1)*16 + IG
      RWGT(IND) = _ONE_
    ENDDO
  ENDIF
ENDDO

!  Initialize arrays for combined Planck fraction data.

DO IPT = 1,13
  DO IPR = 1, JPGPT
    FREFA(IPR,IPT) = _ZERO_
    FREFADF(IPR,IPT) = _ZERO_
  ENDDO
ENDDO
DO IPT = 1,6
  DO IPR = 1, JPGPT
    FREFB(IPR,IPT) = _ZERO_
    FREFBDF(IPR,IPT) = _ZERO_
  ENDDO
ENDDO

!  Reduce g-points for relevant data in each LW spectral band.

CALL RRTM_CMBGB1
CALL RRTM_CMBGB2
CALL RRTM_CMBGB3
CALL RRTM_CMBGB4
CALL RRTM_CMBGB5
CALL RRTM_CMBGB6
CALL RRTM_CMBGB7
CALL RRTM_CMBGB8
CALL RRTM_CMBGB9
CALL RRTM_CMBGB10
CALL RRTM_CMBGB11
CALL RRTM_CMBGB12
CALL RRTM_CMBGB13
CALL RRTM_CMBGB14
CALL RRTM_CMBGB15
CALL RRTM_CMBGB16

RETURN
END SUBROUTINE RRTM_INIT_140GP
