SUBROUTINE SUCST(KULOUT,KDAT,KSSS,KPRINTLEV)

!**** *SUCST * - Routine to initialize the constants of the model.

!     Purpose.
!     --------
!           Initialize and print the common YOMCST + initialize
!         date and time of YOMRIP.

!**   Interface.
!     ----------
!        *CALL* *SUCST (..)

!        Explicit arguments :
!        --------------------

!        KULOUT  - logical unit for the output
!        KDAT    - date in the form AAAAMMDD
!        KSSS    - number of seconds in the day
!        KPRINTLEV - printing level

!        Implicit arguments :
!        --------------------
!        COMMON YOMCST
!        COMMON YOMRIP

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        Additions : 90-07-30 (J.-F. Geleyn)
!                    91-11-15 (M. Deque)
!                    96-08-12 M.Hamrud - Reduce printing
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE YOMCST   , ONLY : RPI      ,RCLUM    ,RHPLA    ,RKBOL    ,&
            &RNAVO    ,RDAY     ,REA      ,REPSM    ,RSIYEA   ,&
            &RSIDAY   ,ROMEGA   ,RA       ,RG       ,R1SA     ,&
            &RSIGMA   ,RI0      ,R        ,RMD      ,RMV      ,&
            &RMO3     ,RD       ,RV       ,RCPD     ,RCPV     ,&
            &RCVD     ,RCVV     ,RKAPPA   ,RETV     ,RCW      ,&
            &RCS      ,RLVTT    ,RLSTT    ,RLVZER   ,RLSZER   ,&
            &RLMLT    ,RTT      ,RATM     ,RDT      ,RESTT    ,&
            &RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
            &RGAMS    ,RALPD    ,RBETD    ,RGAMD
USE YOMRIP   , ONLY : RTIMST   ,RTIMTR

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KDAT
INTEGER_M :: KPRINTLEV
INTEGER_M :: KSSS
INTEGER_M :: KULOUT


!     LOCAL INTEGER SCALARS
INTEGER_M :: IA, ID, IDAT, IM, ISSS, J

!     LOCAL REAL SCALARS
REAL_B :: ZDE, ZET, ZJU, ZRS, ZRSREL, ZTETA, ZTI


#include "fctast.h"
#include "fcttrm.h"
#include "fcttim.h"
!      -----------------------------------------------------------------

!*       1.    DEFINE FUNDAMENTAL CONSTANTS.
!              -----------------------------


RPI=_TWO_*ASIN(_ONE_)
RCLUM=299792458._JPRB
RHPLA=6.6260755E-34_JPRB
RKBOL=1.380658E-23_JPRB
RNAVO=6.0221367E+23_JPRB

!     ------------------------------------------------------------------

!*       2.    DEFINE ASTRONOMICAL CONSTANTS.
!              ------------------------------

RDAY=86400._JPRB
REA=149597870000._JPRB
REPSM=0.409093_JPRB

RSIYEA=365.25_JPRB*RDAY*_TWO_*RPI/6.283076_JPRB
RSIDAY=RDAY/(_ONE_+RDAY/RSIYEA)
ROMEGA=_TWO_*RPI/RSIDAY

IDAT=KDAT
ISSS=KSSS
ID=NDD(IDAT)
IM=NMM(IDAT)
IA=NCCAA(IDAT)
ZJU=RJUDAT(IA,IM,ID)
ZTI=RTIME(IA,IM,ID,ISSS)
RTIMST=ZTI
RTIMTR=ZTI
ZTETA=RTETA(ZTI)
ZRS=RRS(ZTETA)
ZDE=RDS(ZTETA)
ZET=RET(ZTETA)
ZRSREL=ZRS/REA

!     ------------------------------------------------------------------

!*       3.    DEFINE GEOIDE.
!              --------------

RG=9.80665_JPRB
RA=6371229._JPRB
R1SA=REAL(_ONE_/REAL(RA,KIND(_ONE_)),KIND(R1SA))

!     ------------------------------------------------------------------

!*       4.    DEFINE RADIATION CONSTANTS.
!              ---------------------------

RSIGMA=_TWO_ * RPI**5 * RKBOL**4 /(15._JPRB* RCLUM**2 * RHPLA**3)
RSIGMA= 5.670373E-8_JPRB ! W m-2 K-4

RI0=1360.8_JPRB ! XF XF XF

!#EU RI0=RI0*0.85

!     ------------------------------------------------------------------

!*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
!              ------------------------------------------

R=RNAVO*RKBOL
RMD=28.9644_JPRB
RMV=18.0153_JPRB
RMO3=47.9942_JPRB
RD=1000._JPRB*R/RMD
RV=1000._JPRB*R/RMV
RCPD=3.5_JPRB*RD
RCVD=RCPD-RD
RCPV=4._JPRB *RV
RCVV=RCPV-RV
RKAPPA=RD/RCPD
RETV=RV/RD-_ONE_

!     ------------------------------------------------------------------

!*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
!              ---------------------------------------------

RCW=4218._JPRB

!     ------------------------------------------------------------------

!*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
!              --------------------------------------------

RCS=2106._JPRB

!     ------------------------------------------------------------------

!*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
!              ----------------------------------------------------

RTT=273.16_JPRB
RDT=11.82_JPRB
RLVTT=2.5008E+6_JPRB
RLSTT=2.8345E+6_JPRB
RLVZER=RLVTT+RTT*(RCW-RCPV)
RLSZER=RLSTT+RTT*(RCS-RCPV)
RLMLT=RLSTT-RLVTT
RATM=100000._JPRB

!     ------------------------------------------------------------------

!*       9.    SATURATED VAPOUR PRESSURE.
!              --------------------------

RESTT=611.14_JPRB
RGAMW=(RCW-RCPV)/RV
RBETW=RLVTT/RV+RGAMW*RTT
RALPW=LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
RGAMS=(RCS-RCPV)/RV
RBETS=RLSTT/RV+RGAMS*RTT
RALPS=LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
RGAMD=RGAMS-RGAMW
RBETD=RBETS-RBETW
RALPD=RALPS-RALPW

!     ------------------------------------------------------------------

!*      10.    PRINTS

IF (KPRINTLEV >= 1) THEN
  WRITE(KULOUT,'(''0*** Constants of the ICM   ***'')')
  WRITE(KULOUT,'('' *** Fundamental constants ***'')')
  WRITE(KULOUT,'(''           PI = '',E13.7,'' -'')')RPI
  WRITE(KULOUT,'(''            c = '',E13.7,''m s-1'')')RCLUM
  WRITE(KULOUT,'(''            h = '',E13.7,''J s'')')RHPLA
  WRITE(KULOUT,'(''            K = '',E13.7,''J K-1'')')RKBOL
  WRITE(KULOUT,'(''            N = '',E13.7,''mol-1'')')RNAVO
  WRITE(KULOUT,'('' *** Astronomical constants ***'')')
  WRITE(KULOUT,'(''          day = '',E13.7,'' s'')')RDAY
  WRITE(KULOUT,'('' half g. axis = '',E13.7,'' m'')')REA
  WRITE(KULOUT,'('' mean anomaly = '',E13.7,'' -'')')REPSM
  WRITE(KULOUT,'('' sideral year = '',E13.7,'' s'')')RSIYEA
  WRITE(KULOUT,'(''  sideral day = '',E13.7,'' s'')')RSIDAY
  WRITE(KULOUT,'(''        omega = '',E13.7,'' s-1'')')ROMEGA

  WRITE(KULOUT,'('' The initial date of the run is :'')')
  WRITE(KULOUT,'(1X,I8,1X,I5,5X,I4,1X,I2,1X,I2)')IDAT,ISSS,IA,IM,ID
  WRITE(KULOUT,'('' The Julian date is : '',F11.2)') ZJU
  WRITE(KULOUT,'('' Time of the model  : '',F15.2,'' s'')')ZTI
  WRITE(KULOUT,'('' Distance Earth-Sun : '',E13.7,'' m'')')ZRS
  WRITE(KULOUT,'('' Relative Dist. E-S : '',E13.7,'' m'')')ZRSREL
  WRITE(KULOUT,'('' Declination        : '',F12.5)') ZDE
  WRITE(KULOUT,'('' Eq. of time        : '',F12.5,'' s'')')ZET
  WRITE(KULOUT,'('' ***         Geoide         ***'')')
  WRITE(KULOUT,'(''      Gravity = '',E13.7,'' m s-2'')')RG
  WRITE(KULOUT,'('' Earth radius = '',E13.7,'' m'')')RA
  WRITE(KULOUT,'('' Inverse E.R. = '',E13.7,'' m'')')R1SA
  WRITE(KULOUT,'('' ***        Radiation       ***'')')
  WRITE(KULOUT,'('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4'')')  RSIGMA
  WRITE(KULOUT,'('' Solar const. = '',E13.7,'' W m-2'')')RI0
  WRITE(KULOUT,'('' *** Thermodynamic, gas     ***'')')
  WRITE(KULOUT,'('' Perfect gas  = '',e13.7)') R
  WRITE(KULOUT,'('' Dry air mass = '',e13.7)') RMD
  WRITE(KULOUT,'('' Vapour  mass = '',e13.7)') RMV
  WRITE(KULOUT,'('' Ozone   mass = '',e13.7)') RMO3
  WRITE(KULOUT,'('' Dry air cst. = '',e13.7)') RD
  WRITE(KULOUT,'('' Vapour  cst. = '',e13.7)') RV
  WRITE(KULOUT,'(''         Cpd  = '',e13.7)') RCPD
  WRITE(KULOUT,'(''         Cvd  = '',e13.7)') RCVD
  WRITE(KULOUT,'(''         Cpv  = '',e13.7)') RCPV
  WRITE(KULOUT,'(''         Cvv  = '',e13.7)') RCVV
  WRITE(KULOUT,'(''      Rd/Cpd  = '',e13.7)') RKAPPA
  WRITE(KULOUT,'(''     Rv/Rd-1  = '',e13.7)') RETV
  WRITE(KULOUT,'('' *** Thermodynamic, liquid  ***'')')
  WRITE(KULOUT,'(''         Cw   = '',E13.7)') RCW
  WRITE(KULOUT,'('' *** thermodynamic, solid   ***'')')
  WRITE(KULOUT,'(''         Cs   = '',E13.7)') RCS
  WRITE(KULOUT,'('' *** Thermodynamic, trans.  ***'')')
  WRITE(KULOUT,'('' Fusion point  = '',E13.7)') RTT
  WRITE(KULOUT,'('' RTT-Tx(ew-ei) = '',E13.7)') RDT
  WRITE(KULOUT,'(''        RLvTt  = '',E13.7)') RLVTT
  WRITE(KULOUT,'(''        RLsTt  = '',E13.7)') RLSTT
  WRITE(KULOUT,'(''        RLv0   = '',E13.7)') RLVZER
  WRITE(KULOUT,'(''        RLs0   = '',E13.7)') RLSZER
  WRITE(KULOUT,'(''        RLMlt  = '',E13.7)') RLMLT
  WRITE(KULOUT,'('' Normal press. = '',E13.7)') RATM
  WRITE(KULOUT,'('' Latent heat :  '')')
  WRITE(KULOUT,'(10(1X,E10.4))') (10._JPRB*J,J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (RLV(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (RLS(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'('' *** Thermodynamic, satur.  ***'')')
  WRITE(KULOUT,'('' Fusion point = '',E13.7)') RTT
  WRITE(KULOUT,'(''      es(Tt)  = '',e13.7)') RESTT
  WRITE(KULOUT,'('' es(T) :  '')')
  WRITE(KULOUT,'(10(1X,E10.4))') (10._JPRB*J,J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (ESW(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (ESS(RTT+10._JPRB*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E10.4))') (ES (RTT+10._JPRB*J),J=-4,4)
ENDIF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SUCST






