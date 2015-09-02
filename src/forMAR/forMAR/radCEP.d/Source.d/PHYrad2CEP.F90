       subroutine PHYrad2CEP(rklonr,rklevr,rn_Aer,yrmmdd,irhhss   &
      &                     ,rDiST ,rAlbe ,rp_CEP,rphCEP,rcdCEP   &
      &                     ,emsCEP,lsmCEP,cszCEP,rLaCEP,rLoCEP   &
      &                     ,rAeCEP,rO3CEP,rqvCEP,rqiCEP,rqwCEP   &
      &                     ,rswCEP,rqrCEP,rtaCEP,rtsCEP          &
      &                     ,rFIRnc,rFIRnt,rFSOnc,rFSOnt,rFSOst   &
      &                     ,rCD_OD,rCDtOD,rAe_OD,rAetOD,YYYY     &
! #DB &                                                 ,k2ii,k2jj&
      &)

! +------------------------------------------------------------------------+
! | MAR PHYSICS                                            17-11-2004  MAR |
! |                                                                        |
! |                                                                        |
! |   SubRoutine PHYrad2CEP     interfaces MAR        with the    new      |
! |              ECMWF Solar/Infrared      Radiative  Transfer Scheme      |
! |                                                                        |
! |                                                                        |
! |   Content:   CALL of - ECMWF Code initializing the Radiation Transfert |
! |                      - ECMWF                       Radiation Transfert |
! |                                                                        |
! |                                                                        |
! |   ECMWF Code Source:  J.-J. Morcrette, 28 nov 2002                     |
! |                                                                        |
! |                                                                        |
! +------------------------------------------------------------------------+

#include "tsmbkind.h"

! +--Global Variables (ECMWF)
! +  ========================

USE PARRTM1D , ONLY : JP_LON   ,JP_IDIA  ,JP_FDIA  ,JP_TDIA  ,&
 &            JP_LEV ,JP_LW    ,JP_SW    ,JP_NUA   ,JP_MODE  ,&
 &            JP_AER ,JP_LEVP1
USE YOMCST   , ONLY : RD       ,RG       ,RTT      ,RSIGMA   ,&
 &            RCPD   ,RPI      ,RDAY     ,REA      ,RI0      ,&
 &            REPSM  ,RMD      ,RKBOL    ,RNAVO    ,R        ,&
 &            RLVTT  ,RLSTT
USE YOERAD   , ONLY : NSW      ,NTSW     ,NRADFR             ,&
 &            LRRTM  ,LINHOM   ,LOIFUEC  ,LTEMPDS  ,LOWASYF  ,&
 &            LOWHSSS,LONEWSW  ,LNEWAER  ,LHVOLCA            ,&
 &            NRADIP ,NRADLP   ,NOZOCL                       ,&
 &            NICEOPT,NLIQOPT  ,NOVLP    ,NHOWINH  ,RMINICE
USE YOEAERD  , ONLY : CVDAES   ,CVDAEL   ,CVDAEU   ,CVDAED            ,&
 &            RCTRBGA,RCVOBGA  ,RCSTBGA  ,RCAEOPS  ,RCAEOPL  ,RCAEOPU ,&
 &            RCAEOPD,RCTRPT   ,RCAEADK  ,RCAEADM  ,RCAEROS
USE YOERDI   , ONLY : RCH4     ,RN2O     ,RO3      ,RCFC11   ,&
 &                                                  RCFC12
USE YOERDU   , ONLY : RCDAY    ,R10E     ,DIFF     ,REPLOG   ,&
 &            REPSC  ,REPSCO   ,REPSCQ   ,REPSCT   ,REPSCW   ,&
 &            NTRAER
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &            R4IES  ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &            RALVDCP,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU



      IMPLICIT NONE



! +--INPUT (from MAR/ECMWF Interface)
! +  -----

      integer rklonr                       ! nb de pts de grilles
      integer rklevr,kj                    ! nb de niveaux
      integer rn_Aer,nAe                   ! nb d'a?rosols
! #DB integer k2ii(rklonr),kio
! #DB integer k2jj(rklonr),kjo
! #DB integer jkjllw      ,lijio

      integer yrmmdd,YYYY                  ! Date   in the form  yyyyMMdd
      integer irhhss                       ! Number of seconds in the day

      real    rDiST                        ! Distance Soleil-Terre     [UA]
      real    rLaCEP(rklonr)               ! Latitude              [radian]
      real    rLoCEP(rklonr)               ! Longitude             [radian]
      real    rAlbe (rklonr)               ! Albedo     de la Surface
      real    emsCEP(rklonr)               ! Emissivit? de la Surface
      real    lsmCEP(rklonr)               ! Land/Sea Mask    (1.=land)  
      real    cszCEP(rklonr)               ! cos(Solar zenithal Dist.)
      real    rp_CEP(rklonr,rklevr)        ! Pressure (Layer)          [Pa]
      real    rphCEP(rklonr,rklevr+1)      ! Pressure (Layer Interface)[Pa]
      real    rcdCEP(rklonr,rklevr)        ! Cloud Fraction (dropplets)
      real    rAeCEP(rklonr,rn_Aer,rklevr) ! Aerosol Concentr.
      real    rO3CEP(rklonr,rklevr)        ! O3      Concentr.
      real    rqvCEP(rklonr,rklevr)        ! Vapor   Concentr.      [kg/kg]
      real    rqiCEP(rklonr,rklevr)        ! Cryst.  Concentr.      [kg/kg]
      real    rqwCEP(rklonr,rklevr)        ! Droppl. Concentr.      [kg/kg]
      real    rswCEP(rklonr,rklevr)        ! Saturat. % Water       [kg/kg]
      real    rqrCEP(rklonr,rklevr)        ! Drops   Concentr.      [kg/kg]
      real    rtaCEP(rklonr,rklevr)        ! Temperat.(Layer)           [K]
      real    rtsCEP(rklonr)               ! Temperat.(Surface+Air)     [K]


! +--OUTPUT
! +  ------

      real    rFIRnc(rklonr,rklevr+1)      ! CLEAR-SKY LW NET FLUXES [W/m2]
      real    rFIRnt(rklonr,rklevr+1)      ! TOTAL-SKY LW NET FLUXES [W/m2]
      real    rFSOnc(rklonr,rklevr+1)      ! CLEAR-SKY SW NET FLUXES [W/m2]
      real    rFSOnt(rklonr,rklevr+1)      ! TOTAL-SKY SW NET FLUXES [W/m2]
      real    rFSOst(rklonr)               ! TOTAL-SKY SW SRF DOWN F [W/m2]
      real    rCD_OD(rklonr,rklevr)        ! Cloud Optical Depth, 1st Intrv.
      real    rCDtOD(rklonr)               ! Cloud Optical Depth, 1st Intrv.
                                           !  (Vertical Integral)
      real    rAe_OD(rklonr,rklevr)        ! All Aerosol Optical Depth
      real    rAetOD(rklonr)               ! All Aerosol Optical Depth
                                           !  (Vertical Integral)


! +--INTERNAL VARIABLES
! +  ------------------

! +--For Use in radlsw
! +  ^^^^^^^^^^^^^^^^^
      REAL_B :: PGELAM5(JP_LON)
      REAL_B ::  PGEMU5(JP_LON)
      REAL_B ::  PSLON5(JP_LON)
      REAL_B ::  PCLON5(JP_LON)
      REAL_B ::  ZOZON5(JP_LON,JP_LEV)
      REAL_B ::   ZAER5(JP_LON,JP_AER,JP_LEV)
!
!
      INTEGER_M :: KIDIA ,KFDIA ,KTDIA ,KLON  ,KLEV  
      INTEGER_M :: KMODE ,KAER  ,KSW

      INTEGER_M :: KBOX  ,NBOX
      INTEGER_M :: NDUMP ,ILWRAD
!
      REAL_B ::  PRII05

      REAL_B ::   PAER5(JP_LON,JP_AER,JP_LEV) ! Aerosol Optical Depth
      REAL_B ::  PALBD5(JP_LON,JP_SW)  
      REAL_B ::  PALBP5(JP_LON,JP_SW)
!
      REAL_B ::   PAPH5(JP_LON,JP_LEVP1) 
      REAL_B ::    PAP5(JP_LON,JP_LEV)
!
      REAL_B ::  PCCO25
      REAL_B ::  PCLFR5(JP_LON,JP_LEV) 
      REAL_B ::    PDP5(JP_LON,JP_LEV)
      REAL_B ::  PEMIS5(JP_LON)
      REAL_B ::  PEMIW5(JP_LON)
      REAL_B ::   PLSM5(JP_LON)
      REAL_B ::   PMU05(JP_LON)
      REAL_B ::  POZON5(JP_LON,JP_LEV)
      REAL_B ::     PQ5(JP_LON,JP_LEV)
!
      REAL_B ::  PQIWP5(JP_LON,JP_LEV) 
      REAL_B ::  PQLWP5(JP_LON,JP_LEV)        ! Dropplets    Concentration
      REAL_B ::  PSQIW5(JP_LON,JP_LEV)        ! Ice Crystals Concentration
      REAL_B ::  PSQLW5(JP_LON,JP_LEV)        ! 
      REAL_B ::    PQS5(JP_LON,JP_LEV)
      REAL_B :: PQRAIN5(JP_LON,JP_LEV)
      REAL_B :: PRAINT5(JP_LON,JP_LEV)
      REAL_B :: PRLVRI5(JP_LON,JP_LEV)
      REAL_B :: PRLVRL5(JP_LON,JP_LEV)
      REAL_B ::    PTH5(JP_LON,JP_LEVP1) 
      REAL_B ::     PT5(JP_LON,JP_LEV) 
      REAL_B ::    PTS5(JP_LON)
      REAL_B ::  PNBAS5(JP_LON)        
      REAL_B ::  PNTOP5(JP_LON)
!
      REAL_B ::  PEMIT5(JP_LON)
      REAL_B ::   PFCT5(JP_LON,JP_LEVP1)
      REAL_B ::   PFLT5(JP_LON,JP_LEVP1)
      REAL_B ::   PFCS5(JP_LON,JP_LEVP1)
      REAL_B ::   PFLS5(JP_LON,JP_LEVP1)
      REAL_B :: PFRSOD5(JP_LON)
      REAL_B ::  PSUDU5(JP_LON)
      REAL_B ::  PUVDF5(JP_LON)
      REAL_B ::  PPARF5(JP_LON)
      REAL_B ::  PFDCT5(JP_LON,JP_LEVP1)
      REAL_B ::  PFUCT5(JP_LON,JP_LEVP1)
      REAL_B ::  PFDLT5(JP_LON,JP_LEVP1)
      REAL_B ::  PFULT5(JP_LON,JP_LEVP1)
      REAL_B ::  PFDCS5(JP_LON,JP_LEVP1)
      REAL_B ::  PFUCS5(JP_LON,JP_LEVP1)
      REAL_B ::  PFDLS5(JP_LON,JP_LEVP1)
      REAL_B ::  PFULS5(JP_LON,JP_LEVP1)
!
      REAL_B ::  ZTAU5 (JP_LON,JP_SW,JP_LEV)  ! Cloud Optical Depth
      REAL_B ::  ZTAUI5(JP_LON)               ! Cloud Optical Depth (vert.int.)
!
      REAL_B ::  ASWBOX(JP_LON,100), OLRBOX(JP_LON,100)
      REAL_B ::  SLWBOX(JP_LON,100), SSWBOX(JP_LON,100)
      REAL_B ::  TAUBOX(JP_LON,100), CLDBOX(JP_LON,100,JP_LEV)

! +--For Use in SUCLD
! +  ^^^^^^^^^^^^^^^^
      REAL_B ::  ZETA(JP_LEV)      ,  ZETAH(JP_LEVP1)


! +--For Use in SUOVLP
! +  ^^^^^^^^^^^^^^^^^
      REAL_B :: ZTVIR              ,   ZFACT
      REAL_B ::   ZAZ(JP_LEV)      ,    ZAZH(JP_LEVP1)



! +--Local  Variables
! +  ================

      LOGICAL  RADini
      LOGICAL  RADin2
      common/c_RADini/RADini  , RADin2

      REAL_B ::       ZEPAER                               ! Generic [O3]
      REAL_B ::       RGAMMAS                              ! CLD Fract. Param.
      common/rbRADini/ZEPAER,RGAMMAS

      REAL_B ::       RTIMTR  , ZTHETOZ , ZANGOZC

      INTEGER_M ::    JL      , JK      , JAER    ,JNU
      INTEGER_M ::    KULOUT  , NINDAT  , NSSSSS  ,KPRTLEV
      INTEGER_M ::    IYR     , MONTH   , IDAY    ,IMINUT
      INTEGER_M ::    KPRINT  , SKSHIFT

      CHARACTER*5 :: RCP_CMIP5
      common/c_RCP_CMIP5/RCP_CMIP5


! +--Load External Functions
! +  =======================
#include "fctast.h"
#include "fcttim.h"
#include "fcttre.h"


      KPRTLEV= 1
      KPRINT = 1
      SKSHIFT= 0

! #DB kio = 58
! #DB kjo = 50


! +--Time Base
! +  =========

      NINDAT = yrmmdd
      NSSSSS = irhhss
      IYR    =   NAA(NINDAT)
      MONTH  =   NMM(NINDAT)
      IDAY   =   NDD(NINDAT)
      RTIMTR = RTIME(IYR,MONTH,IDAY,NSSSSS)
      IMINUT =   INT(FLOAT(NSSSSS)/60.)



! +--Basic Initialization
! +  ====================

! +--Dimensions (auxiliary variables)
! +  --------------------------------

        KIDIA  = 1                       ! DO'NT CHANGE
        KFDIA  = JP_LON                  ! Nb Columns
        KTDIA  = 1                       !
        KLON   = rklonr                  ! Nb Columns
        KLEV   = rklevr                  ! Nb Levels
        KMODE  = JP_MODE                 ! Used in Planck Fcts Specification
        KAER   = JP_AER                  !

! +--Nb of Solar Spectral Intervals
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        KSW    =  JP_SW       ! SW Nb of Spectral Intervals (max is JP_SW=6)
        NSW    =  JP_SW       ! SW Nb of Spectral Intervals (max is JP_SW=6)
        NTSW   =  JP_SW       ! SW Nb of Spectral Intervals (max is JP_SW=6)

        KBOX   = 0                       !                                      \VER
        NBOX   = 1                       !                                      \VER

        ILWRAD = 1       ! 0: Morcrette,     1991 operational before 20000627
                         ! 1: Mlawer et al., 1997 now ECMWF-operational
                         ! 2: Morcrette,     1991 original as in ERA-15'

        NDUMP  = 3       ! No Print
!       NDUMP  = 2       ! 1D Results
!       NDUMP  = 1       ! Debug
!       NDUMP  = 0       ! ALL
        KULOUT = 6       ! Output Device for SUCST


! +--Verification of the Dimensions
! +  ------------------------------

      IF (.NOT.RADini)                                              THEN

        write(6,*) 'INITIALISATION OF ECMWF RADIATIVE TRANSFERT: BEGIN'


        IF    (rklonr/=JP_LON)                                      THEN
          write(6,6001) rklonr,JP_LON
 6001     format(' @!#& BAD SET-UP of dimensions (rklonr,JP_LON);  = ('&
     &                                            ,i3,',',i3,')')
          STOP
        END IF
        IF    (rklevr/=JP_LEV)                                      THEN
          write(6,6002) rklevr,JP_LEV
 6002     format(' @!#& BAD SET-UP of dimensions (rklevr,JP_LEV);  = ('&
     &                                            ,i3,',',i3,')')
          STOP
        END IF
        IF    (rn_Aer/=JP_AER)                                      THEN
          write(6,6003) rn_Aer,JP_AER
 6003     format(' @!#& BAD SET-UP of dimensions (rn_Aer,JP_AER);  = ('&
     &                                            ,i3,',',i3,')')
          STOP
        END IF


! +--Mathematical Constants
! +  ----------------------

        REPLOG = 1.E-12                  ! Minimum Logarithm Argument
        REPSC  = 1.E-12
        REPSCO = 1.E-12
        REPSCQ = 1.E-12
        REPSCT = 1.E-12
        REPSCW = 1.E-12


! +--Switches (general)
! +  ------------------

        LONEWSW= .TRUE.  ! .TRUE. SWSN radiative routine is     used
       IF (ILWRAD.EQ.1)     THEN
        LRRTM  = .TRUE.  ! .TRUE. RRTM radiative routine is     used
       ELSE
        LRRTM  = .FALSE. ! .FALSE.RRTM radiative routine is NOT used
       END IF
        LTEMPDS= .FALSE. ! .TRUE. ALLOWS FOR SURF. T DISCONTIN. IN RAD.COMPUT.

        NTRAER = 19      ! NUMBER OF TRANSMISSION FUNCTIONS  W OR W/O AEROSOLS


! +--Switches (Clouds Optical Properties)
! +  ------------------------------------

        LINHOM = .FALSE. ! Tiedke (1995) correct. factor (0.7) of tau not used
        NHOWINH=  2      ! Tau correction factor:           exp(-(sig/tau)^2)
                         !  (used if LINHOM = .TRUE.) 
        LOIFUEC= .FALSE. ! .FALSE. IF ICE   CLOUDS AS EBERT-CURRY (LW & SW)
        LOWASYF= .FALSE. ! .FALSE. IF WATER CLOUDS AS FOUQUART         (SW)
        LOWHSSS= .FALSE. ! .FALSE. IF WATER CLOUDS AS SMITH-SHI   (LW)
        NRADIP =  3      ! Ice    effective Radius:
                         !  0   fixed        40 microns
                         !  1   f(T)   40 - 130 microns
                         !  2   f(T)   30 -  60 microns Jakob-Klein
                         !  3   f(T,IWC)                Sun-Rikus,     1999
        RMINICE=  15     ! Minimum Diameter for Ice Particles (micronm)
                         ! Needed only if   NRADIP = 3
                         ! (see von Walden et al., 2003 (Oct) Tab. 2 p.1393)
        NRADLP =  0      ! Liquid effective Radius: f(Pressure)              !+CA+! (0>2)
                         !  0 effective radius - liquid as f(Pressure)
                         !  1 fixed 10 microns over land, 13 over ocean
                         !  2 computed from LWC          Martin et al, 1994
        NLIQOPT=  1      ! Cloud Optical Properties (Water): 1=ECMWF Operat.
                         !  0  LW: Smith-Shi,   1992; SW: Fouquart,    1987
                         !  1  LW: Savijarvi,   1997; SW: Slingo  ,    1989
                         !  2  LW: Lindner,Li,  2000; SW: Slingo  ,    1989
!XF
! WARNING: NLIQOPT= 2 increases SWD and LWD but MAR is more unstable.

        NICEOPT=  2      ! Cloud Optical Properties (Ice)  : 1=ECMWF Operat. !+CA+! (2>0)

                         !  0  LW: Smith,Shi  , 1992; SW: Ebert-Curry, 1992
                         !  1  LW: Ebert,Curry, 1992; SW: Ebert-Curry, 1992
                         !  2  LW: Fu,Liou    , 1993; SW: Fu,Liou    , 1993
                         !  3  LW: Fu et al.  , 1998; SW: Fu         , 1996
        NOVLP  =  2      ! CLOUD OVERLAP CONFIGURATION:
                         !  1=MRN, 2=MAX, 3=RAN, 4=Hogan


! +--Switches (Aerosols/O3)
! +  ----------------------
        LNEWAER= .TRUE.  ! Climatology of Aerosols: TEGEN ET AL. 1997 / GADS
        LHVOLCA= .TRUE.  ! .TRUE. IF GISS HISTORY OF VOLCANIC AEROSOLS IS ON
!       NOZOCL = -1      ! TESTS the vertical quadrature       (NO absorber)
!       NOZOCL =  0      ! whatever is read for O3 as input profile
!       NOZOCL =  1      ! OLD         ECMWF O3 climatology and    aerosols
        NOZOCL =  2      ! Fortuin-Langematz O3 climatology and    aerosols
!       NOZOCL =  3      ! OLD         ECMWF O3 climatology and NO aerosols
!       NOZOCL =  4      ! Fortuin-Langematz O3 climatology and NO aerosols
 

! +--BASIC CONSTANTS
! +  ---------------

!            *****
        CALL SUCST (KULOUT, NINDAT, NSSSSS, KPRTLEV) ! Initialize common YOMCST
!            *****                                   !  (Basic Constants)
                                                     ! Initialize common YOMRIP
                                                     !  (only date and time)

! +--YOENCST - THERMODYNAMIC TRANSITION OF PHASE
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        RTT=273.16                                   !

! +--YOETHF  - DERIVED CONSTANTS SPECIFIC TO ECMWF THERMODYNAMICS
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        RTWAT=RTT                                    ! 
        RTICE=RTT-23.                                !

        RGAMMAS=0.02

! +--YOERDU  - CONTROL, PARAMETERS AND SECURITY IN RADIATION
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        RCDAY  = RDAY * RG / RCPD                    !
        R10E   = 0.4342945                           ! DECIMAL /  NATURAL
                                                     ! LOG.        FACTOR
        DIFF   = 1.66                                ! DIFFUSIVITY FACTOR

        ZEPAER = 1.E-12                              ! Generic [O3]


! +--Space/Time Independant Coefficients
! +  -----------------------------------

        CALL SURDI(YYYY,RCP_CMIP5)                    ! ECMWF Surface Albedo, Emissivity
        CALL SULWN            ! Initialize common YOELW (new LW Coeff.)
        CALL SUOLW            ! Initialize common YOELW (old LW Coeff.)

! +--Initialization routine for RRTM
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        CALL SURRTAB          ! AER'S RRTM LW RADIATION
        CALL SURRTPK          ! Initialize common YOERRTWN 
                              !  (k-coefficients in spectral intervals)
        CALL SURRTRF          ! Initialize common YOERRTRF
                              !  (RRTM Reference Atmosphere)
        CALL SURRTFTR         ! Initialize common YOERRTRF

!            RRTM routine     ! BAND   [cm-1] ! low           ! high
!       ----------------------+---------------+---------------+---------
        CALL RRTM_KGB1        !  1:   10- 250 ! H2O           ! H2O
        CALL RRTM_KGB2        !  2:  250- 500 ! H2O           ! H2O
        CALL RRTM_KGB3        !  3:  500- 630 ! H2O,CO2       ! H2O,CO2
        CALL RRTM_KGB4        !  4:  630- 700 ! H2O,CO2       ! O3,CO2
        CALL RRTM_KGB5        !  5:  700- 820 ! H2O,CO2       ! O3,CO2
        CALL RRTM_KGB6        !  6:  820- 980 ! H2O           ! nothing
        CALL RRTM_KGB7        !  7:  980-1080 ! H2O,O3        ! O3
        CALL RRTM_KGB8        !  8: 1080-1180 ! (i.e.>~300mb) ! O3
                              !               ! H2O           !
        CALL RRTM_KGB9        !  9: 1180-1390 ! H2O,CH4       ! CH4
        CALL RRTM_KGB10       ! 10: 1390-1480 ! H2O           ! H2O
        CALL RRTM_KGB11       ! 11: 1480-1800 ! H2O           ! H2O
        CALL RRTM_KGB12       ! 12: 1800-2080 ! H2O,CO2       ! nothing
        CALL RRTM_KGB13       ! 13: 2080-2250 ! H2O,N2O       ! nothing
        CALL RRTM_KGB14       ! 14: 2250-2380 ! CO2           ! CO2
        CALL RRTM_KGB15       ! 15: 2380-2600 ! N2O,CO2       ! nothing
        CALL RRTM_KGB16       ! 16: 2600-3000 ! H2O,CH4       ! nothing

! +--Reduce absorption coefficient data from 256 to 140 g-points
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        CALL RRTM_INIT_140GP

! +--Initialization routine for SW (6 spectral interval resolution)
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        CALL SUSWN ( NTSW  , KSW  )      ! Initialize common YOESW


        radINI=.TRUE.

        write(6,*) 'INITIALISATION OF ECMWF RADIATIVE TRANSFERT: END  '


      END IF                             ! End of Basic Initialization


! +--Radiation: Global (Time dependant) Parameters
! +  =============================================

        PRII05 = RI0/(rDiST*rDiST)       ! INSOLATION
        PCCO25 = 360.E-06*44./29.        ! CONCENTRATION IN CO2 (PA/PA)



! +--Surface Properties
! +  ==================

! +--Surface Albedo
! +  --------------

  DO jnu=1,KSW
      DO jl=1,KLON
        PALBD5(JL,JNU)    = rAlbe(jl)
        PALBP5(JL,JNU)    = rAlbe(jl)
      ENDDO
  ENDDO


! +--Surface Emissivity
! +  ------------------

      DO jl=1,KLON
        PEMIS5(JL)   =emsCEP(jl)
        PEMIW5(JL)   =emsCEP(jl)


! +--Land/sea Mask
! +  -------------

        PLSM5 (JL)   =lsmCEP(jl)


! +--Cosine (Solar zenithal Distance)
! +  --------------------------------

        PMU05 (JL)   =cszCEP(jl)
      ENDDO



! +--Atmospheric Thermodynamics (Time and Space dependant)
! +  =====================================================

! +--Pressure
! +  --------

     JK=1+KLEV
      DO jl=1,KLON
        PAPH5 (JL,JK)     = rphCEP(jl,jk)
      ENDDO
  DO JK=1,KLEV
      DO JL=1,KLON
        PAPH5 (JL,JK)     = rp_CEP(jl,jk)
        PAP5  (JL,JK)     = rphCEP(jl,jk)
        PDP5  (JL,JK)     = rphCEP(jl,jk+1)-rphCEP(jl,jk)


! +--Water Species      Distributions
! +  --------------------------------

        PQ5    (JL,JK)    = rqvCEP(jl,jk)
        PQIWP5 (JL,JK)    = max(0.,rqiCEP(jl,jk)-1.E-9)
        PQLWP5 (JL,JK)    = max(0.,rqwCEP(jl,jk)-1.E-9)
        PQS5   (JL,JK)    = rswCEP(jl,jk)
        PQRAIN5(JL,JK)    = rqrCEP(jl,jk)
        PRAINT5(JL,JK)    = 0.                   !                              \VER
        PRLVRI5(JL,JK)    = 0.                   ! e-mail J.-J.M. 20031203
        PSQIW5 (JL,JK)    = 1.                   !    exp(-PRLVRI5(JL,JK))
        PRLVRL5(JL,JK)    = 0.                   ! e-mail J.-J.M. 20031203
        PSQLW5 (JL,JK)    = 1.                   !    exp(-PRLVRL5(JL,JK))


! +--Cloud Fraction
! +  --------------

!XF
        PCLFR5(JL,JK)     = rcdCEP(jl,jk)                    ! dropplets
        PCLFR5(JL,JK)     = (PQIWP5(JL,JK)+PQLWP5(JL,JK)) &  ! ECMWF Paramet.
     &                     /(RGAMMAS      *  PQS5(JL,JK))    !  (VERY Crude)
        PCLFR5(JL,JK)     = min( _ONE_    ,PCLFR5(JL,JK))    !
        PCLFR5(JL,JK)     = max(1.0E-3    ,PCLFR5(JL,JK)) &  ! no small values
     &      *max(0.,sign(1.0,rqiCEP(JL,JK)+rqwCEP(JL,JK)-2.E-9))
        rcdCEP(jl,jk)     = PCLFR5(JL,JK)                    !
!       write(4,4) JL,JK,PQLWP5(JL,JK),PQIWP5(JL,JK)      &
!    &                  ,  PQS5(JL,JK)                    &
!    &                  ,rcdCEP(jl,jk),PCLFR5(JL,JK)
!4      format(2i6,' Cloud Liq.W.= ', f9.6,5x             &
!    &            ,' Cloud Sol.W.= ', f9.6,5x             &
!    &            ,' Satur.W.Vap.= ', f9.6,5x             &
!    &            ,' Cloud Fract.= ',2f9.6)


! +--Temperature        Distribution
! +  -------------------------------

        PT5    (JL,JK)    = rtaCEP(jl,jk)
        PTH5   (JL,JK)    =(rtaCEP(jl,jk)+rtaCEP(jl,max(1,jk-1)))*0.5
      END DO
  END DO
      DO JL=1,KLON
        PTH5   (JL,KLEV+1)= rtsCEP(jl)
        PTS5   (JL)       = rtsCEP(jl)
      END DO


! +--Convective Layer
! +  ----------------

      DO JL=1,KLON
        PNBAS5 (JL)       = 1.                !                                 \VER
        PNTOP5 (JL)       = 1.                !                                 \VER
      END DO



! +--Initialization (Climatologies, Time independant)
! +  ================================================

      IF (.NOT.RADin2)                                              THEN

! +--Aerosols Radiative Characteristics (YOEAER)
! +  ----------------------------------

!            *******
        CALL SUAERL                  ! Aerosols LW Radiative Charact.
        CALL SUAERSN ( NTSW , KSW )  ! Aerosols SW Radiative Charact.
!            *******


! +--Aerosols Optical Thickness Horizontal Distribution (model grid
! +  --------------------------------------------------  independant)

!            ********
        CALL SUAERH                  ! 
        CALL SUECAEBC                ! BLACK CARBON (URBAN/FOREST FIRE ORIGIN)
        CALL SUECAEOR                ! ORGANIC-TYPE
        CALL SUECAESD                ! SOIL-DUST                       ORIGIN
        CALL SUECAESS                ! SEA -SALT                       ORIGIN
        CALL SUECAESU                ! SULFATE-TYPE
!            ********


! +--Clouds (YOECLD)
! +  ---------------

        DO jk=1,klev
          ZETA(jk) =PAP5(1,jk) /PAPH5(1,klev+1)
        ENDDO
        DO jk=1,klev+1
          ZETAH(jk)=PAPH5(1,jk)/PAPH5(1,klev+1)
        ENDDO

!            *****
        CALL SUCLD  ( KLEV  , ZETA ) 
!            *****


! +--Cloud Optical Parameters SW/LW (all parameterizations)
! +  ------------------------------------------------------

!            *******
        CALL SUCLOPN ( NTSW , KSW , KLEV )  ! Initialize YOECLOP
!            *******


! +--Radar Reflectivity
! +  ------------------

         ZAZH(KLEV+1)=     0.
          ZAZ(1)     =100000.
           JL=KIDIA
        DO jk=KLEV,2,-1
          ZTVIR      =       PT5(JL,jk)   /(1.-0.608*PQ5(jl,jk))
          ZFACT      = LOG(PAPH5(jl,jk+1))-    LOG(PAPH5(jl,jk))
          ZAZH(jk)   =      ZAZH(jk+1)    + R *    ZTVIR/(RMD*RG)*ZFACT
           ZAZ(jk)   = 0.5*(ZAZH(jk+1)+ZAZH(jk))*1000.
        END DO

!            ******
        CALL SUOVLP  ( KLEV , ZAZ )         ! Initialize ALPHA1 (%radar refl.)
!            ******                         !  (Hogan & Illingsworth, 1999)


! +--NO Absorber
! +  -----------

        IF (NOZOCL.EQ.-1) THEN
                RCH4             = 1.E-18
                RN2O             = 1.E-18
                RO3              = 1.E-18
                RCFC11           = 1.E-18
                RCFC12           = 1.E-18
                PCCO25           = 1.E-18
          DO jk=1,klev
              DO jl=KIDIA,KFDIA
                POZON5(JL,JK)    = 0.
              ENDDO
          ENDDO
          DO JK=1,KLEV
            DO JAER=1,KAER
              DO JL=KIDIA,KFDIA
                PAER5(JL,JAER,JK)= ZEPAER
              ENDDO
            ENDDO
          ENDDO
        END IF

      END IF                                            ! End Initialization



! +--Initialization (Climatologies, Time   dependant)
! +  ================================================

! +--Aerosols Optical Thickness Horizontal Distribution (model grid
! +  --------------------------------------------------  independant)

!            *******
        CALL SUECAEC ( NINDAT, IMINUT )     ! TEGEN ET AL. (1997, JGR 102, 
!            *******                        !               pp23895-23915)
 
! +--Aerosols Optical Thickness Vertical   Distribution
! +  --------------------------------------------------
!            *******
        CALL SUAERV&
                 & ( KLEV   ,ZETAH                                   &
                 & , CVDAES ,CVDAEL ,CVDAEU ,CVDAED                  &
                 & , RCTRBGA,RCVOBGA,RCSTBGA,RCAEOPS,RCAEOPL,RCAEOPU &
                 & , RCAEOPD,RCTRPT ,RCAEADK,RCAEADM,RCAEROS         &
                 & )
!            *******

! +--O3
! +  --   

! +--Fortuin-Langematz O3 climatology
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!            *******
        CALL SUECOZC ( NINDAT , IMINUT )
!            *******

! +--ECMWF   Geleyn    O3 Climatology
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ZTHETOZ=RTETA( RTIMTR  )
        ZANGOZC=  REL( ZTHETOZ ) - 1.7535

!            *******
        CALL SUECOZO ( ZANGOZC )
!            *******



! +--Interpolation on the MAR Grid
! +  =============================

        DO JL=1,KLON
          PGELAM5(JL)=rLoCEP(JL)
           PGEMU5(JL)=SIN(rLaCEP(JL))
           PCLON5(JL)=COS(rLoCEP(JL))
           PSLON5(JL)=SIN(rLoCEP(JL))
        END DO

!              ******
         CALL  RADACA                                       &
          &( KIDIA , KLON   , KLON  , KTDIA , KLEV          &
          &, PAPH5 , PGELAM5, PGEMU5, PCLON5, PSLON5, PTH5  &
          &, ZAER5 , ZOZON5                                 &
          &  )
!              ******


! +--OLD         ECMWF O3 Climatology
! +  --------------------------------

        IF (NOZOCL.EQ.1 .OR. NOZOCL.EQ.3)                           THEN
          DO JK=1,KLEV
            DO JL=KIDIA,KFDIA
              POZON5(JL,JK)      = ZOZON5(JL,JK)
            ENDDO
          ENDDO
        END IF


! +--FORTUIN LANGEMATZ O3 Climatology
! +  --------------------------------

        IF (NOZOCL.EQ.2 .OR. NOZOCL.EQ.4)                           THEN

!              ******
          CALL RADOZC ( KIDIA ,KLON  , KLON   , KTDIA, KLEV            &
                     &, KPRINT,KLON  , SKSHIFT, PAPH5, PGEMU5, ZOZON5)
!              ******

          DO JK=1,KLEV
            DO JL=KIDIA,KFDIA
              POZON5(JL,JK)      = ZOZON5(JL,JK)
            ENDDO
          ENDDO
 
        END IF

! +--AEROSOLS
! +  --------

        IF (NOZOCL.EQ.1 .OR. NOZOCL.EQ.2)                           THEN
          DO jk=1,klev
            DO jaer=1,KAER
              DO jl=KIDIA,KFDIA
                PAER5(JL,JAER,JK)= ZAER5(JL,JAER,JK)
              ENDDO
            ENDDO
          ENDDO
        END IF


! +--NO AEROSOLS
! +  -----------

        IF (NOZOCL.GT.2)                                            THEN
          DO jk=1,klev
            DO jaer=1,KAER
              DO jl=KIDIA,KFDIA
                PAER5(JL,JAER,JK)= ZEPAER
              ENDDO
            ENDDO
          ENDDO
        END IF


! +--SECURITY CHECK ON AEROSOL AMOUNTS
! +  ---------------------------------

        DO JK=1,KLEV
          DO JAER=1,KAER
            DO JL=KIDIA,KFDIA
              PAER5(JL,JAER,JK)=MAX(ZEPAER,PAER5(JL,JAER,JK))
            ENDDO
          ENDDO
        ENDDO



! +--Transmission to MAR Variables
! +  =============================

          DO JK=1,KLEV
            DO JL=1,KLON
              rO3CEP(jl,jk)      = ZOZON5(JL,JK)
            END DO
          END DO

        DO JK=1,KLEV
          DO JAER=1,KAER
            DO JL=1,KLON
              rAeCEP(jl,jaer,jk) =  PAER5(JL,JAER,JK)
            END DO
          END DO
        END DO



! +--Solar and IR Transfer through the Atmosphere
! +  ============================================

! +         ***********
            CALL RADLSW                                                 &
     &    ( KIDIA , KFDIA , KLON  , KTDIA , KLEV   , KMODE , KAER,      &
     &      KBOX  , NBOX                                                &
     &    , NDUMP , ILWRAD                                              &
     &    , PRII05                                                      &
     &    , PAER5 , PALBD5, PALBP5, PAPH5 , PAP5                        &
     &    , PCCO25, PCLFR5, PDP5  , PEMIS5, PEMIW5 , PLSM5 ,  PMU05,    &
     &      POZON5                                                      &
     &    , PQ5   , PQIWP5, PQLWP5, PSQIW5, PSQLW5 , PQS5  ,  PQRAIN5,  &
     &      PRAINT5                                                     &
     &    , PRLVRI5,PRLVRL5,PTH5  , PT5   , PTS5   , PNBAS5,  PNTOP5    &
     &    , PEMIT5, PFCT5 , PFLT5 , PFCS5 , PFLS5  , PFRSOD5, PSUDU5,   &
     &      PUVDF5, PPARF5                                              &
     &    , PFDCT5, PFUCT5, PFDLT5, PFULT5, PFDCS5 , PFUCS5,  PFDLS5,   &
     &      PFULS5                                                      &
     &    , ZTAU5 , ZTAUI5                                              &
     &    , ASWBOX, OLRBOX, SLWBOX, SSWBOX, TAUBOX , CLDBOX             &
! #DB&                                                     ,  k2ii,k2jj &
     &    )
! +         ***********

! +--Radiative Fluxes   Distributions
! +  ================================

  DO JK=1,KLEV+1
      DO JL=1,KLON
        rFIRnc (jl,jk)    = PFCT5 (JL,JK)  ! CLEAR-SKY LW NET FLUXES
        rFIRnt (jl,jk)    = PFLT5 (JL,JK)  ! TOTAL-SKY LW NET FLUXES
        rFSOnc (jl,jk)    = PFCS5 (JL,JK)  ! TOTAL-SKY SW NET FLUXES
        rFSOnt (jl,jk)    = PFLS5 (JL,JK)  ! TOTAL-SKY SW NET FLUXES
      END DO
  END DO


! +--Cloud Optical Depth
! +  ===================

  DO JK=1,KLEV
        kj              =KLEV + 1    -JK   ! 
      DO JL=1,KLON
        rCD_OD (jl,kj)  =ZTAU5(JL,  1,JK)  ! Cloud Optical Depth
       DO nAe=1,rn_Aer
        rAe_OD (jl,jk)  =rAe_OD (jl  ,jk) &! Aeros.Optical Depth
     &                  +PAER5(JL,nAe,JK)  !
       END DO
      END DO
  END DO

      DO JL=1,KLON
        rAetOD (jl)     = 0.
      END DO
      DO JL=1,KLON
        rCDtOD (jl)     = ZTAUI5(JL)       ! Cloud Optical Depth (vert.integr., 
                                           !                      1st interval)
        DO JK=1,KLEV
        rAetOD (jl)     = rAetOD (jl)   &  ! Aeros.Optical Depth
     &                  + rAe_OD (jl,jk)   !
        END DO

! +--SURFACE RADIATIVE CHARACTERISTICS (SW)
! +  ======================================

        rFSOst (jl)     = PFRSOD5(JL)      ! TOTAL-SKY SRF SW DOWNWARD FLUX
!       ?      (jl)     =  PSUDU5(JL)      ! SOLAR RADIANCE IN SUN'S DIRECT.
!       ?      (jl)     =  PUVDF5(JL)      ! SURFAC.DOWNWARD U.V. RADIATION \VER
!       ?      (jl)     =  PPARF5(JL)      ! PHOTOSYNTHET. ACTIVE RADIATION \VER


! +--SURFACE RADIATIVE CHARACTERISTICS (LW)
! +  ======================================
!       emsCEP (jl)     =  PEMIT5(JL)      ! TOTAL         LW EMISSIVITY    \VER

      END DO


! +--OUTPUT
! +  ======

! #DB     jkjllw=0
! #DB DO JL=1,KLON
! #DB     lijio =0
! #DB DO JK=1,KLEV
!!      write(6,*) k2ii(JL),k2jj(JL),jl,jk,' rFIRnt: ', PFLT5(jl,jk)
!!      write(6,*) k2ii(JL),k2jj(JL),jl,jk,' rFSOnt: ', PFLS5(jl,jk)
! #DB   IF ( PFLT5(jl,jk).GT. 500..OR. PFLS5(jl,jk).GT. 500. .OR.                                                    &
! #DB &      PFLT5(jl,jk).LT.-500..OR. PFLS5(jl,jk).LT.-500. .OR.                                                    &
! #DB &      (k2ii(jl).EQ.kio.AND.k2jj(jl).EQ.kjo)) lijio=1
! #DB END DO
! #DB IF   (lijio.EQ.1)                                                                                        THEN
! #DB   DO JK=1,KLEV
! #DB     IF (mod(jkjllw,20).EQ.0) write(6,600)
! #DB                                      600 format('IN   PHYrad2CEP: Radiative Fluxes ',/                         &
! #DB &                                               ,'    i    j   JL   JK',9x,'Ta',9x,'Qv',9x,'Qi',9x,'Qw'        &
! #DB &                                               ,9x,'O3',8x,'CLD',8x,'COD',8x,'AOD',8x,'SOn',8x,'IRn')
! #DB     jkjllw=jkjllw+1
! #DB     write(6,601) k2ii(JL),k2jj(JL),JL,JK,rtaCEP(jl,jk),rqvCEP(jl,jk),rqiCEP(jl,jk),rqwCEP(jl,jk),rO3CEP(jl,jk) &
! #DB &                                       ,rcdCEP(jl,jk),rCD_OD(jl,kj),rAe_OD(jl,jk),rFSOnt(jl,jk),rFIRnt(jl,jk)
! #DB             601  format(4i5,10e11.3)
! #DB   END DO
! #DB ENDIF
!
! #DB END DO

      return
      end
