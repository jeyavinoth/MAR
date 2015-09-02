C +************************************************************************+
C |                                                                        |
C |                                                                        |
C |               ************************************                     |
C |               *                                  *                     |
C |               *    N    E    S    T    O    R    *                     |
C |               *                                  *                     |
C |               ************************************                     |
C |                                                                        |
C |                                                                        |
C |     NESTing Organization for the preparation of                        |
C |             meteorological and surface fields in Regional models       |
C |                                                                        |
C |                 \__ _                    ____ /                        |
C |               \_/     ****              /    \                         |
C |               / \    ******            /    / \                        |
C |              /   \  ******            |    /   |                       |
C |                    ....                \  /   / UCL-IAG                |
C |                   .....                 \/___/ LGGE                    |
C |                  .....                  / LTHE                         |
C |                                                                        |
C |                                                                        |
C |        Institut d'Astronomie et de Geophysique Georges Lemaitre        |
C |                                                                        |
C |                  Universite catholique de Louvain                      |
C |                       Chemin du Cyclotron, 2                           |
C |                   1348 Louvain-la-Neuve - BELGIUM                      |
C |                                                                        |
C |                  - - - - - - - - - - - - - - - - -                     |
C |                                                                        |
C |                  L. G. G. E.   -    G R E N O B L E                    |
C |                                                                        |
C |    Laboratoire de Glaciologie et de Géophysique de l'Environnement     |
C |                       Rue Molière, 54 - BP 96                          |
C |                    38402 St-Martin d'Hères CEDEX                       |
C |                                                                        |
C |                                                                        |
C |                  - - - - - - - - - - - - - - - - -                     |
C |                                                                        |
C |                  L. T. H. E.   -    G R E N O B L E                    |
C |                                                                        |
C |   Laboratoire d'Etude des Transferts en Hydrologie et Environnement    |
C |                   Domaine Universitaire - BP 53                        |
C |                    Rue de la Piscine 1023-1025                         |
C |                  38041 Grenoble Cedex 9 - FRANCE                       |
C |                                                                        |
C +************************************************************************+
C |                                                                        |
C |  NESTOR 4.1.5                              Date : 20 June  2004        |
C |  ------------                              ------                      |
C |                                                                        |
C |  Development :                                                         |
C |                                                                        |
C |      Olivier Brasseur  (brasseur@oma.be):                              |
C |         General structure, development (several components).           |
C |                                                                        |
C |      Hubert Gallée     (gallee@lgge.obs.ujf-grenoble.fr):              |
C |         General structure, MAR team manager.                           |
C |                                                                        |
C |      Philippe Marbaix  (marbaix@astr.ucl.ac.be):                       |
C |         Development (1st version, atmospheric data interpolation, LMDz)|
C |                                                                        |
C |      Xavier Fettweis   (fettweis@astr.ucl.ac.be)                       |
C |         Maintainer: version NESTOR > 4.0.0                             |
C |                                                                        |
C +************************************************************************+


      PROGRAM NESTOR


C +************************************************************************+


      IMPLICIT NONE


C +---LSC and NST domain dimensions
C +   -----------------------------

      INCLUDE 'NSTdim.inc'


C +---LSC,INT,NST,SND variables
C +   -------------------------

      INCLUDE 'CTRvar.inc'
      INCLUDE 'LSCvar.inc'
      INCLUDE 'INTvar.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'LSCmod.inc'
      INCLUDE 'SNDvar.inc'
      INCLUDE 'WGEvar.inc'
      INCLUDE 'DESvar.inc'
      INCLUDE 'NESTOR.inc'      
            
C +---Soil and surface data files locations
C +   -------------------------------------
      INCLUDE 'LOCfil.inc'

C +---Local variables
C +   ---------------

      INTEGER  VARSIZE
      EXTERNAL VARSIZE

      LOGICAL  Vtrue

C +---Data
C +   ----

      DATA Vtrue   /  .true.  /

C +---Soil and surface data files locations
C +   -------------------------------------
C +   (code to actually set the paths)
      INCLUDE 'LOCset.inc'

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---START MESSAGE
C +   =============
      

      write(6,*)
      write(6,*)
      write(6,*) '    ***********************************************'
      write(6,*) '    *                                             *'
      write(6,*) '    *         N    E    S    T    O    R          *'
      write(6,*) '    *                                             *'
      write(6,*) '    *   NESTing Organization for the preparation  *'
      write(6,*) '    *   of meteorological and surface fields in   *'
      write(6,*) '    *             3-D Regional models.            *'
      write(6,*) '    *                                             *'
      write(6,*) '    *          Rain disagregation models          *'
      write(6,*) '    *                                             *'
      write(6,*) '    *         Wind gust estimate methods          *'
      write(6,*) '    *                                             *'
      write(6,*) '    ***********************************************'
      write(6,*)
      write(6,*) '              ---   Version  4.1.8   ---           '
      write(6,*) '              ---     20/06/2004     ---           '
      write(6,*)
      write(6,*)


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +   Reading of control file
C +   =======================

      OPEN (unit=20,status='old',file='NSTing.ctr')

      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(i4)')    SELECT
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(A3)')    LABLio
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(A60)')   NSTdir
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(A3)')    LSCmod
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(A3)')    NSTmod
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(A3)')    region
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(i4,3(x,i2))') RUNiyr,RUNmma,RUNjda,RUNjhu
      read (20,'(i10,x,i2)')                 DURjda,DURjhu
      read (20,'(i13)')                             FORjhu
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    LoutDA
      read (20,'(l4)')    ASCfor 
      read (20,'(l4)')    LoutLS 
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    SPHgrd
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*)         HORint
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*)         VERint
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    TOPetopo
      read (20,'(a4)')    TOP30 
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    TOPcst
      read (20,'(l4)')    TOPcstLSC
      read (20,'(l4)')    TOPdomLSC
      read (20,'(l4)')    TOPcst0
      read (20,'(l4)')    TOPfilt
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    CORzz6
      read (20,'(l4)')    CORsurf
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    RUGdat
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    VEGdat
      read (20,'(l4)')    VEGcor
      read (20,'(l4)')    VEGbel
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    NDV1km
      read (20,'(l4)')    NDV8km
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    SVTmod
      read (20,*)         SVTwet
      read (20,'(l4)')    SVTlsc
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    SSTrey
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    SNDing
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(A57)')   SNDfil
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(l4)')    CLDcor
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(i4)')    DESGpr
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,'(i4)')    WGEmet
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -
      read (20,*) !- - - - - - - - - - - - - - - - - -

      CLOSE(unit=20)

                             GTOPO30=.false.
      IF (TOP30 .eq. 'T   ') GTOPO30=.true.
      
C +        ******
      CALL WARNms
C +        ******


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      WRITE(6,*) "NESTOR characteristics"
      WRITE(6,*) "~~~~~~~~~~~~~~~~~~~~~~"

C +---MODELS
C +   ======

      WRITE(6,'(A,$)') ' Forcing fields  (LSC) : '

      IF (LSCmod.eq.'ECM') WRITE(6,*) 'ECMWF - Hybrid levels'
      IF (LSCmod.eq.'E15') WRITE(6,*) 'ERA-15 - Hybrid levels'
      IF (LSCmod.eq.'E40') WRITE(6,*) 'ERA-40 - Hybrid levels'
      IF (LSCmod.eq.'E20') WRITE(6,*) 'ERA-20C - Hybrid levels'
      IF (LSCmod.eq.'EIN') WRITE(6,*) 'ERA-Interim - Hybrid levels'
      IF (LSCmod.eq.'ECP') WRITE(6,*) 'ECMWF - Pressure levels'
      IF (LSCmod.eq.'MAR') WRITE(6,*) 'MAR'
      IF (LSCmod.eq.'LMD') WRITE(6,*) 'LMD - phys. out. version (old)'
      IF (LSCmod.eq.'LMz') WRITE(6,*) 'LMDz - standard LSC input file'
      IF (LSCmod.eq.'NCP') WRITE(6,*) 'NCEP'
      IF (LSCmod.eq.'ALA') WRITE(6,*) 'ALADIN'
      IF (LSCmod.eq.'CAN') WRITE(6,*) 'CanESM2 (CMIP5)'
      IF (LSCmod.eq.'CM3') WRITE(6,*) 'HadCM3 (ICE2SEA)'
      IF (LSCmod.eq.'EM5') WRITE(6,*) 'ECHAM5 (ICE2SEA)'
      IF (LSCmod.eq.'NOR') WRITE(6,*) 'NorESM1 (CMIP5)'
      IF (LSCmod.eq.'CSI') WRITE(6,*) 'CSIRO-Mk3 (CMIP5)'
      IF (LSCmod.eq.'BCC') WRITE(6,*) 'BCC-CSM1-1 (CMIP5)'
      IF (LSCmod.eq.'20C') WRITE(6,*) 
     .    '20th Century Reanalysis V2 (NOAA-CIRES)'
      IF (LSCmod.eq.'MIR') WRITE(6,*) 'MIROC5 (CMIP5)'
      IF (LSCmod.eq.'CM5') WRITE(6,*) 'CNRM-CM5 (CMIP5)'
      IF (LSCmod.eq.'AC3') WRITE(6,*) 'ACCESS1-3 (CMIP5)'
      IF (LSCmod.eq.'NC1') WRITE(6,*) 'NCEP-NCAPv1'
      IF (LSCmod.eq.'NC2') WRITE(6,*) 'NCEP-NCAPv2'


      WRITE(6,'(A,$)') ' Nested model    (NST) : '

      IF (NSTmod.eq.'MAR') WRITE(6,*) 'MAR'
      IF (NSTmod.eq.'M2D') WRITE(6,*) 'MAR - 2D version'
      IF (NSTmod.eq.'GRA') WRITE(6,*) 'GRADS'
      IF (NSTmod.eq.'CPL') WRITE(6,*) 'SVAT - Coupling'
      
      WRITE(6,*)

      IF (LSCmod.eq.'MAR'.or.LSCmod.eq.'LMD') THEN
       REGgrd=.false.      ! Non-regular input grid (lat/long)
C      Note: might not be necessary for LMD ! 
C      The current LMz input = ok for zoom & staggered grid
      ELSE
       REGgrd=.true.       ! Regular input grid (lat/long)
      ENDIF

      IF (LSCmod.EQ.'LMz'.OR.LSCmod.EQ.'LMD') THEN 
         M30d=.true.
         WRITE(6,*) 'LMD standard model: months have 30 days'
      ELSE
         M30d=.false.
      ENDIF

                           f28d=.false. 
      IF (LSCmod.EQ.'CM3') M30d=.true.
      IF (LSCmod.EQ.'CAN') f28d=.true.
      IF (LSCmod.EQ.'NOR') f28d=.true.
      IF (LSCmod.EQ.'MIR') f28d=.true.
      IF (LSCmod.EQ.'CSI') f28d=.true.
      IF (LSCmod.EQ.'BCC') f28d=.true.

      Sinclair=.false.
      Alpert  =.false.
      IF (DESGpr.eq.1) Sinclair=.true.
      IF (DESGpr.eq.2) Alpert  =.true.

C +---REGION
C +   ======
   
      WRITE(6,'(A,$)') ' Region                : '

      IF (region.eq.'GRD') then
          WRITE(6,*) 'GRD (Greenland)'
          USRreg='GRD'
      ENDIF
      IF (region.eq.'ANT') then
          WRITE(6,*) 'ANT (Antarctic)'
          region='GRD'
          USRreg='ANT'
      ENDIF
      IF (region.eq.'EUR') WRITE(6,*) 'EUR (Europe)'
      IF (region.eq.'AFW')                                          THEN
                           WRITE(6,*) 'AFW (West Africa)'
          IF (mw.NE.3)     STOP '#@!&  mw badly specified / must be = 3'
      END IF
      IF (region.eq.' NO') WRITE(6,*) 'No specified'

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CONTROL OF DATES
C +   ================   


C +---Initial, incremental, and final dates of run
C +   --------------------------------------------

C +        ******
      CALL DATcnv (RUNiyr,RUNmma,RUNjda,RUNjhu,DATini,Vtrue)
C +        ******
C +...Initial date of run

      DAT_dt=FORjhu
C +...Time interval

      DATfin=DATini+24*DURjda+DURjhu
C +...End date of run


C +---Number of steps
C +   ---------------

      DATstp=(DATfin-DATini)/DAT_dt + 1
      IF (SNDing) DATstp=1


C +---Initialisation
C +   --------------

      DATtim=DATini


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---HORIZONTAL GRID IN NESTED MODEL
C +   ===============================


C +---Modele Atmospherique Regional (MAR)
C +   -----------------------------------

      IF ((NSTmod.eq.'MAR'.or.NSTmod.eq.'M2D').and.(SELECT.ne.2)) THEN

C +         ******
       CALL MARhgd
C +         ******

      ENDIF


C +---GRADS graphic output
C +   --------------------

      IF ((NSTmod.eq.'GRA').and.(SELECT.ne.2)) THEN

C +         ******
       CALL GRAhgd
C +         ******

      ENDIF


C +---Hydrology - meteo coupling
C +   --------------------------

      IF ((NSTmod.eq.'CPL').and.(SELECT.ne.2)) THEN

C +         ******
       CALL CPLhgd
C +         ******

      ENDIF


C +---Rain disagregation model (RDM)
C +   ------------------------------

      IF (SELECT.eq.2) THEN

C +         ******
       CALL DEShgd 
C +         ******

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---TOPOGRAPHY SOURCES
C +   ==================
C +
C +---Bamber data set (1 km)/Racmo2 topo
C +   ----------------------------------
      IF (USRreg.eq."ANT") THEN
C +            ******
          CALL USRant('Bamber') !+ 'Bamber'/'Racmo2' +!
C +            ******
C +
C +---ETOPO data set (1/5 minutes)
C +   ----------------------------
      ELSE IF (TOPetopo) THEN
C +         ******

c       CALL ETOPOg ! ETOPO 5min
       CALL ETOPO1 ! ETOPO 1min
       CALL ICEmsk 

C +         ******
C +
C +---GTOPO data set (30 secondes)
C +   ----------------------------
      ELSE IF (GTOPO30)  THEN
C +         ******
       CALL GTOP30
C +         ******
      ENDIF 


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CORRECTION OF PRESCRIBED TOPOGRAPHY
C +   ===================================


      IF (SELECT.eq.1) THEN


C +---Border of constant NST topography at boundaries
C +   -----------------------------------------------

       IF (TOPcst) THEN

        TOPopt=1
C +          ******
        CALL TOPcor (TOPopt)
C +          ******

       ENDIF


C +---Zero topography in the relaxation zone
C +   --------------------------------------

       IF (TOPcst0) THEN

        TOPopt=4
C +          ******
        CALL TOPcor (TOPopt)
C +          ******

       ENDIF


C +---Topography filtering (2D and 3D)
C +   --------------------------------

       IF (TOPfilt) THEN

        TOPopt=5
C +          ******
        CALL TOPcor (TOPopt)
C +          ******

       ENDIF


C +---Other options related to LSC topography
C +   ---------------------------------------

C +...Options TOPcstLSC (TOPopt=2) and TOPdomLSC (TOPopt=3)
C +...will be treated in NSTint subroutine since they require
C +...the specification of the LSC topography.
C +...If used, they will be followed by a filtering (if required)

      ENDIF

  
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---VEGETATION CHACTERISTICS
C +   ========================


C +---Global land cover - IGBP (only Africa and Europe)
C +   -------------------------------------------------

      IF (VEGdat) THEN

C +         ******
       CALL GLOveg
C +         ******

       ENDIF


C +---Vegetation cover in Europe (Corine)
C +   -----------------------------------

      IF (VEGcor) THEN

C +         ******
       CALL CORveg
C +         ******

      ENDIF


C +---Vegetation cover in Belgium
C +   ---------------------------

      IF (VEGbel) THEN

C +         ******
       CALL BELveg
C +         ******

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---SOIL CHACTERISTICS
C +   ==================


C +---FAO soil types classification
C +   -----------------------------

C +        ******
c     CALL FAOsol
C +        ******

C +---GSWP soil types classification
C +   ------------------------------

C +        ******
      CALL GSWPsl
C +        ******

C +---Surface characteristics
C +   -----------------------

C +        ******
      CALL SOLdom
C +        ******

  
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---CORRECTED VEGETATION FRACTION WITH MAX. NDVI INDEX
C +   ==================================================


C +---NDVI index (1 or 8 km resolution)
C +   ---------------------------------

      IF ((VEGdat.or.VEGcor).and.(NDV1km.or.NDV8km)) THEN

C +         ******
       CALL GLOfrc
C +         ******

       ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO iter=1,DATstp    ! TEMPORAL ITERATION

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---SELECT A LSC DATA FILE
C +   ======================


      IF (.not.SNDing.or.SVTmod) THEN

C +         ******
       CALL LSCinp
C +         ******

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---MAIN TREATMENT
C +   ==============


C +---1. Interpolation forcing fields for NST model
C +   ---------------------------------------------

      IF (SELECT.eq.1.and.(.not.SNDing.or.SVTmod)) THEN

C +         ******
       CALL NSTint
C +         ******

       IF (SSTrey) THEN

C +          ******
        CALL SSTint
C +          ******

       ENDIF

      ENDIF  ! {SELECT.eq.1.and.(.not.SNDing)}


C +---2. Meteorological initialisation based on sounding
C +   --------------------------------------------------

      IF (SELECT.eq.1.and.SNDing) THEN

C +... Old soundings from IRM
C +         ******
c      CALL SOUNDg
C +         ******

C +... New soundings availables from WEB sites
C +         ******
       CALL SNDweb
C +         ******

      ENDIF  ! {SELECT.eq.1.and.SNDing}


C +---3. Rain disagregation models
C +   ----------------------------

      IF (SELECT.eq.2) THEN

C +         ******
       CALL PRCdes
C +         ******

      ENDIF


C +---4. Wind Gust Estimate Methods
C +   -----------------------------

      IF (SELECT.eq.3) THEN

C +         ******
       CALL NSTint
C +         ******
c      CALL WGustE (WGEmet,LSCmod,NST__u,NST__v,NST__w,NST__t,
c    .              NST_st,NST_qv,NST_zz,NST_sh,NST_z0,NST__p,
c    .              NSTtke,NSTuts,NST_qt,WGEest,WGElwb,WGEupb,
c    .              WGEsta,WNDsta,WGE_zi)
C +         ******

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CORRECTION OF SURFACE CHARACTERISTICS
C +   =====================================


      IF (iter.eq.1.and.SELECT.eq.1) THEN

C +         ******
       CALL SOLdom
C +         ******

      ENDIF

C +...Note : this call is useful only if NSTint subroutine has
C +   modified the surface types (specified in NSTsol variable)

  
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CORRECTION OF FIELDS IN THE SURFACE LAYER
C +   =========================================


      IF (CORsurf.and.(.not.SNDing).and.SELECT.eq.1) THEN

C +         ******
       CALL SL_cor
C +         ******

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---PROGNOSTIC AND ADDITIONAL VARIABLES FOR SVAT MODEL
C +   ==================================================


      IF (SVTmod.and.SELECT.eq.1) THEN

C +         ******
       CALL SVTpar
C +         ******

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CORRECTED GREEN LEAF FRACTION WITH NDVI INDEX
C +   =============================================


C +---NDVI index (1 or 8-km resolution)
C +   ---------------------------------

      IF (VEGdat.and.(NDV1km.or.NDV8km).and.SELECT.eq.1) THEN

C +         ******
       CALL GLOglf
C +         ******

       ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---OUTPUT FILES
C +   ============


C +---Modele Atmospherique Regional (MAR)
C +   -----------------------------------

      IF ((NSTmod.eq.'MAR'.or.NSTmod.eq.'M2D'.or.NSTmod.eq.'CPL').and.
     .    SELECT.le.2) THEN

C +         ******
       CALL MARout
C +         ******

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---GRAPHIC (NetCDF) FILES
C +   ======================


C +---Standard fields
C +   ---------------

      IF (LoutLS.and.SELECT.eq.1) THEN

C +         ******
       CALL NSTout
C +         ******

      ENDIF


C +---Precipitation fields (only for disagregation)
C +   ---------------------------------------------

      IF (SELECT.eq.2) THEN

C +         ******
       CALL PRCout
C +         ******

      ENDIF


C +---Wind gust estimates and wind distribution
C +   -----------------------------------------

      IF (SELECT.eq.3) THEN

C +         ******
       CALL WGEout
C +         ******

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DATtim=DATtim+DAT_dt

      ENDDO               ! TEMPORAL ITERATION

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +        ******
      CALL DATcnv (RUNiyr,RUNmma,RUNjda,RUNjhu,DATtim,Vtrue)
C +        ******      

      open (unit=1,status='replace',file='NST.OK')
      write (1,111) RUNjda,RUNmma,RUNiyr,RUNjhu  
  111 format("NESTOR execution stopped normaly the",i2,"/",i2,"/",i4,
     .       " at ",i2,"h.")
      close(1)  


      write(6,*)
      write(6,*)
      write(6,*) '    ***********************************************'
      write(6,*)
      write(6,*) '            END    OF    N  E  S  T  O  R          '
      write(6,*)
      write(6,*)


      END
