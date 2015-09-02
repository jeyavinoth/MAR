      subroutine SBCnew

C +------------------------------------------------------------------------+
C | ____________________________                                           |
C | MAR OUTPUT   Generic Routine                      Mon  6-Feb-2012  MAR |
C | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                           |
C |                                                                        |
C |                                                    vvvvvvvvvvv         |
C | CAUTION:    (THIS    ROUTINE MAY BE ADAPTED / see "DEVELOPMENT" below) |
C |                                                    ^^^^^^^^^^^         |
C |                                                                        |
C | CAUTION:    (dTiANI, dTiAWS MUST BE ADAPTED / see  data         below) |
C |                                                                        |
C |                                                                        |
C | DOMAIN    OUTPUT                           (in  Subroutine ANI_nc     )|
C |              Dynamics      Statistics                                  |
C |              Precipitation Statistics                                  |
C |              Surface       Statistics                                  |
C |              Wind Gust     Estimates                                   |
C |           !! PREDEFINES TIME INTERVAL      (see       data dTiANI     )|
C |                                                                        |
C | DOMAIN    OUTPUT                           (in  Subroutine MIP_nc     )|
C |              Zonal         Averages        (AMMA-MIP)                  |
C |                                                                        |
C | AWS       OUTPUT                           (in  Subroutine AWS_nc     )|
C |           !! PREDEFINES TIME INTERVAL      (see       data dTiAWS     )|
C |           !! PREDEFINES AWS LOCATIONS      (see block data AWS_nc_DATA)|
C |                                                                        |
C | AWS CHARACTERISTICS ARE OUTPUT ON ASCII        FILE SBCnew.AWS         |
C |                                   ferret shell FILE SBCnew.AWS.JNL     |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C | DEVELOPMENT: New AWS's characteristics may be implemented:             |
C | ^^^^^^^^^^^  (See last part of the block data AWS_nc_DATA)             |
C |              (containing the DATA  of AWS characteristics)             |
C |              (                     organised in 7 columns)             |
C |              ( - Short Name:             1st      column )             |
C |              ( - Full  Name in 2 arrays: 2th, 3th columns)             |
C |              ( - Geographic Coordinates: 4th--6th columns)             |
C |              ( - OUTPUT Status:          7th      column )             |
C |              (                  0: NO     OUTPUT         )             |
C |              (                  1: netcdf OUTPUT possible)             |
C |              (                  2: ascii  OUTPUT possible)             |
C |                                                                        |
C |              The CRITERION for determining the Grid Point Coordinates  |
C |                            corresponding to the AWS Location           |
C |                            MAY BE MODIFIED (search string "CRITERION") |
C |                                                                        |
C |              BE SURE TO CONSEQUENTLY MODIFY THE parameter n_AWS        |
C |                                             IN ALL SUBROUTINES         |
C |                                                                        |
C |              PRECISION: Conversion from real   to real   under vi      |
C |                         :.,$g/real  /s//real*8/g                       |
C |                         Conversion from real   to real   under vi      |
C |                         :.,$g/real[*]8/s//real  /g                     |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARCTR.inc'
      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
      include 'MAR_LB.inc'

      include 'MAR_DY.inc'
      include 'MAR_HY.inc'
      include 'MAR_CA.inc'
      include 'MAR_TE.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MAR_IB.inc'

      logical           WRIlog,wriAWS,wriMIP
      common/SBCnew_log/WRIlog,wriAWS,wriMIP


C +--OUTPUT Time Interval
C +  --------------------

      REAL              dTiANI,dTiMIP,dTiAWS
      common/SBCnew_rea/dTiANI,dTiMIP,dTiAWS


C +--For OUTPUT
C +  ----------

      integer           iprANI,nprANI,nboANI
      common/ANI_nc_arg/iprANI,nprANI,nboANI

      integer           iprMIP,nprMIP,nboMIP
      common/MIP_nc_arg/iprMIP,nprMIP,nboMIP

      integer           iprAWS,nprAWS,nboAWS
      common/AWS_nc_arg/iprAWS,nprAWS,nboAWS

      real          Tair_1(mx,my,2),TNois1(mx,my)
      real          Tair_2(mx,my,2),TNois2(mx,my)
      common/TNoise/Tair_1         ,TNois1
     .             ,Tair_2         ,TNois2


C +--Local   Variables
C +  =================

      integer     isn   ,isl   ,misl_2,nisl_2,isnh
      integer     n     ,l


C +--Snow Pack
C +  ---------

      real        deglon,deglat
      real        dzWE
      real        dz_dSN(-nsno:1)
      REAL        SnowZZ( nsno+1)
      REAL        rosNEW
      REAL        Ti_AAv,dTi   ,arg_zz,t_year,arg_tt,TiClim


C +--For Snow Pack Initialization
C +  ----------------------------

      integer     iteSVA
      integer     iyrSVA,mmaSVA,jdaSVA,jhuSVA
      integer     IOiSVA(iptx) ,IOjSVA(iptx) 
      integer     ifrSVA(imx,jmx,nvx)
      real        CaSSVA(imx,jmx,nvx)
      real        CaWSVA(imx,jmx,nvx)
      real        psvSVA(imx,jmx,nvx)
      real        psgSVA(imx,jmx,nvx)
      real        TgrSVA(imx,jmx,nvx)
      real        evpSVA(imx,jmx)
      real        draSVA(imx,jmx)
      real        rofSVA(imx,jmx)

      logical     MARsno
      DATA        MARsno/.false./

      logical     UPDful
      DATA        UPDful/.true./


C +--Snow Model Initialisation
C +  =========================

      IF (itexpe.eq.       0)                                       THEN

       include 'SBCnew.SNO'


C +--Read DATA from a previous simulation
C +  ------------------------------------

       IF      (MARsno)                                             THEN
          OPEN (unit=33,status='old'    ,form='unformatted',
     .                                   file='MARsvt.INI')

            READ(33) iteSVA
            READ(33) iyrSVA,mmaSVA,jdaSVA,jhuSVA
              write(6,600)iteSVA,iyrSVA,mmaSVA,jdaSVA,jhuSVA
 600          format(
     .        /,' vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',
     .           'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',
     .           'vvvvvvvvvvvvvvvvvvvvvvvvvvv'
     .        /,' ++WARNING++ MARsvt from former simulation (after',i9,
     .           ' iterations, on',i3,'-',i2,'-',i4,' at',i3,'h00m00s',
     .           ') is used to initialize Snow Pack'
     .        /,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',
     .           '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',
     .           '^^^^^^^^^^^^^^^^^^^^^^^^^^^')

            READ(33) IOiSVA  ! IO i Index
            READ(33) IOjSVA  ! IO j Index
            READ(33) isolTV  ! Soil       Type Index
            READ(33) iwafTV  ! =0 ==> no Water Flux
                             ! =1 ==> free Drainage
            READ(33) AlbSTV  ! Dry Soil       Albedo
            READ(33) ivegTV  ! Vegetation Type Index
            READ(33) ifrSVA  ! Vegetation Class Coverage
            READ(33) alaiTV  !       Leaf Area Index                [-]
            READ(33) glf_TV  ! Green Leaf Fraction                  [-]
            READ(33) TsolTV  ! Soil Temperature                     [K]
            READ(33) eta_TV  ! Soil Moisture                    [m3/m3]


            READ(33) CaSSVA  ! Canopy Intercepted Snow  Content[m w.e.]
            READ(33) CaWSVA  ! Canopy Intercepted Water Content [kg/m2]
            READ(33) psvSVA  ! Vegetation   Hydraulic Potential     [m]
            READ(33) psgSVA  ! Ground       Hydraulic Potential     [m]
            READ(33) TvegTV  ! Skin      Vegetation Temperature     [K]
            READ(33) TgrSVA  ! Skin            Soil Temperature     [K]

            READ(33) evpSVA  ! Total Evapotranspiration [mm w.e.]
            READ(33) draSVA  ! Drainage            Flow    [mm/s]
            READ(33) rofSVA  ! Integrated Drainage Flow      [mm]

C +--Snow Pack Characteristics
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~
            READ(33) mskSNo  ! Snow/Ice Type Index            [-]
            READ(33) nssSNo  ! Nb Snow and     Ice Layers
            READ(33) issSNo  ! Nb Superimposed Ice Layers
            READ(33) nisSNo  ! Nb              Ice Layers
            READ(33) nhsSNo  ! Snow History                   [-]
            READ(33) dzsSNo  ! Snow Layers  Thickness         [m]
            READ(33) rosSNo  ! Snow Volumic Mass          [kg/m3]
            READ(33) wasSNo  ! Snow Water   Content       [kg/kg]
            READ(33) tisSNo  ! Snow Temperature               [K]
            READ(33) g1sSNo  ! Snow Dendricity / Sphericity   [-]
            READ(33) g2sSNo  ! Snow Sphericity / Size         [-] [0.0001 m]
            READ(33) agsSNo  ! Snow Age                     [day]
            READ(33) snohSN  ! Snow Buffer    Layer [kg/m2], [mm]
            READ(33) SWaSNo  ! Snow Surficial Water           [m]
            READ(33) zWE0SN  ! Initial      Snow Thickn.[mm w.e.]
            READ(33) zWE_SN  ! Current      Snow Thickn.[mm w.e.]
            READ(33) zWEcSN  ! Non-erodible Snow Thickn.[mm w.e.]
            READ(33) SaltSN  ! u*_th                        [m/s]
            READ(33) SLussl  ! u*_s*                  [kg/kg.m/s]
            READ(33) blowSN  ! NEW  MAX Erosion           [kg/m2]
            READ(33) WEq_SN  ! Added Snow Amount         [m w.e.]
c #ZM       READ(33) SLn_z0  ! Z0_momentum    (instantaneous) [m]
c #ZM       READ(33) SLn_r0  ! Z0_scalar      (instantaneous) [m]
c #BS       READ(33) SLn_b0  ! Z0_erosion     (instantaneous) [m]
c #SZ       READ(33) Z0SaBS  ! Z0      (Sastrugi Height)      [m]

          CLOSE(unit=33)

          DO j=1,my
          DO i=1,mx
          DO n=1,nsx
                TsrfSL(i,j,n)      =  TvegTV(i,j,n)        ! Surface Temperat.
          END DO
          END DO
          END DO


C +--Prescribe Snow Pack Characteristics
C +  -----------------------------------

       ELSE

          WRIlog = .false.

C +--Surface Air Temperature
C +  ~~~~~~~~~~~~~~~~~~~~~~~
          DO j=1,my
          DO i=1,mx
                tairSL(i,j)         = pktaDY(i,j,mz)       !            [K]
     .                  *(pstDYn(i,j)+ptopDY)**cap         !
          END DO
          END DO

C +--Vertical Discretization
C +  ~~~~~~~~~~~~~~~~~~~~~~~
          DO isl= 0,-nsno+1,-1
            misl_2 =     -mod(isl,2)
            nisl_2 =         -isl/2
            dz_dSN(isl)=(((1-misl_2) * 0.001
     .                    +  misl_2  * 0.003) * 10**(nisl_2)) * 2.
            dz_dSN(isl)=max(dz_dSN(isl),0.01)
            dz_dSN(isl)=min(dz_dSN(isl),0.20)
          END DO

          DO j=1,my
          DO i=1,mx
          DO n=1,nsx

C +--Over the Ocean
C +  ~~~~~~~~~~~~~~
            IF (isolSL(i,j).le.2)                     THEN ! NO Snow on Sea-Ice
                nssSNo(i,j,n)      =   0                   ! Nb Snow/Ice Layers
              DO isn = 1,nsno  
                nhsSNo(i,j,n,isn)  =   0                   !            [-]
                dzsSNo(i,j,n,isn)  =   0.                  !            [m]
                rosSNo(i,j,n,isn)  =   0.                  !        [kg/m3]
                wasSNo(i,j,n,isn)  =   0.                  !        [kg/kg]
                tisSNo(i,j,n,isn)  =   0.                  !            [K]
                g1sSNo(i,j,n,isn)  =   0.                  ! [-]        [-]
                g2sSNo(i,j,n,isn)  =   0.                  ! [-] [0.0001 m]
                agsSNo(i,j,n,isn)  =  20.                  !          [day]
              END DO
              DO isn = 1,llx
                TsolTV(i,j,n,isn)  =  SST_LB(i,j)          !            [K]
                eta_TV(i,j,n,isn)  =   1.
              END DO

C +--Over the the Ice Sheet
C +  ~~~~~~~~~~~~~~~~~~~~~~
            ELSE 
     .      IF (isolSL(i,j).EQ.3)                     THEN ! Dry Snow   assumed
                nssSNo(i,j,n)      =   nsno                ! Nb Snow/Ice Layers
                SnowZZ(nsno+1)     =   0.
              DO isn = nsno,1,-1
                dzsSNo(i,j,n,isn)  =  dz_dSN(isn-nsno)     !            [m]
                rosSNo(i,j,n,isn)  =  SnowZZ(      isn+1)  !            [m]
     .                             *  10.                  !(+10kg/m3 each m)  
c #DC.                             +  350.00               ! Dome C Measurem.
     .                             +  250.00               ! Dome F Measurem.
                rosSNo(i,j,n,isn)=    rosSNo(i,j,n,isn)    !
     .                        + 200.*(4000.-sh(i,j))/4000. ! To get 450. at D1

                deglat           =    GElatr(i,j)   /degrad
                deglon           =    GElonh(i,j)   *15.0

!  NANSEN    BLUE ICE AREA
!  ~~~~~~~~~~~~~~~~~~~~~~~
               IF  (deglat .LT.-74.50-0.3*max(0.,deglon-162.) .AND.
     .              deglat .GT.-75.00-0.3*max(0.,deglon-162.) .AND.
     .              deglat .GT.-75.20                         .AND.
     .              deglon .GT.160.00                         .AND.
     .              deglon .LT.164.00                         .AND.
     .              sh(i,j).LT.1000.0                              )
     .          rosSNo(i,j,n,isn)=900.

!  PRIESTLEY BLUE ICE AREA
!  ~~~~~~~~~~~~~~~~~~~~~~~
               IF  (deglat .LT.-73.25-0.3*max(0.,deglon-160.) .AND.
     .              deglat .GT.-73.55-0.3*max(0.,deglon-160.) .AND.
!    .              deglat .GT.-75.20                         .AND.
     .              deglon .GT.161.00                         .AND.
     .              deglon .LT.164.00                         .AND.
     .              sh(i,j).GT. 500.0                         .AND.
     .              sh(i,j).LT.2100.0                              )
     .          rosSNo(i,j,n,isn)=900.


                rosSNo(i,j,n,isn)=min(rosSNo(i,j,n,isn)    !
     .                               ,900.00           )   !
                wasSNo(i,j,n,isn)  =   0.                  !        [kg/kg]
                Ti_AAv             =  250.-40.*sh(i,j)/3200.
                dTi                =       40.
                arg_zz             = (SnowZZ(      isn+1)
     .                           +0.5*dzsSNo(i,j,n,isn)  )
     .                               /2.2
                t_year             =  0.0
                arg_tt             =  2.0 *pi *t_year / 365.
                TiClim             =  Ti_AAv 
     .                         + 0.5*dTi*exp(-arg_zz)
     .                            *cos(arg_tt-arg_zz)
                tisSNo(i,j,n,isn)= min(tairSL(i,j),TiClim) !            [K]
                g1sSNo(i,j,n,isn)  =  99.                  ! [-]        [-]
                g2sSNo(i,j,n,isn)  =   3.                  ! [-] [0.0001 m]
                agsSNo(i,j,n,isn)  =  20.                  !          [day]
                SnowZZ(      isn)  =  SnowZZ(      isn+1)  !
     .                             +  dzsSNo(i,j,n,isn)    !
                IF                   (SnowZZ(      isn).LE.0.05)    THEN
                nhsSNo(i,j,n,isn)  =   0                   ! (MAY BE 2) [-]
                ELSE                                       !
                nhsSNo(i,j,n,isn)  =   2                   !  MAY BE 2! [-]
            !   nhsSNo(i,j,n,isn)  =   2   =>              !
            !   Glazed Snow // Snow previously in Contact with Water  (crude)
            !  (Glazed Snow is generated by Dissipation of Kinetic Energy 
            !                                           of Blown Snow Grains)
                ENDIF
              END DO
              IF (.NOT.WRIlog)                                      THEN
              WRIlog = .true.
              write(6,6001) i,j,n,(             isn ,isn=1,20)
     .                           ,(SnowZZ(      isn),isn=1,20)
     .                           ,(dzsSNo(i,j,n,isn),isn=1,20)
     .                           ,(rosSNo(i,j,n,isn),isn=1,20)
     .                           ,(tisSNo(i,j,n,isn),isn=1,20)
     .                           ,(G1sSNo(i,j,n,isn),isn=1,20)
     .                           ,(G2sSNo(i,j,n,isn),isn=1,20)
     .                           ,(nhsSNo(i,j,n,isn),isn=1,20)
 6001         format(/,'SNOW MODEL INITIALIZATION:'
     .               /,'=========================='
     .               /,3i4,20i6,/,'  Z ',8x,20f6.3,/,'  dz',8x,20f6.3
     .                         ,/,'  ro',8x,20f6.0,/,'  Ts',8x,20f6.1
     .                         ,/,'  G1',8x,20f6.1,/,'  G2',8x,20f6.1
     .                         ,/,'  nh',8x,20i6  ,/)
              END IF
              DO isn = 1,llx
                TsolTV(i,j,n,isn)  =  tairSL(i,j)          !            [K]
                eta_TV(i,j,n,isn)  =   0.
              END DO
                isolTV(i,j)        =  12
                iwafTV(i,j)        =   0
                AlbSTV(i,j)        =   0.55
                ivegTV(i,j,n)      =   0
                alaiTV(i,j,n)      =   0.
                glf_TV(i,j,n)      =   0.
            END IF

C +--Upper Limit of the non erodible Snow (nhsSNo .GT. 1)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               isnh =   0
            DO isn=  nsno,1,-1
               isnh =isnh + isn*min(nhsSNo(i,j,n,isn)-1,1)*max(0,1-isnh)
            ENDDO
               zWE0SN(i,j,n) =    0. 
               zWEcSN(i,j,n) =    0. 
            DO isn=1,nsno
               dzWE          =      dzsSNo(i,j,n,isn) *rosSNo(i,j,n,isn)
               zWE0SN(i,j,n) =      zWE0SN(i,j,n)     +  dzWE
               zWEcSN(i,j,n) =      zWEcSN(i,j,n)     +  dzWE
     .                         *max(0,min(1,isnh+1-isn))
            END DO

C +--General
C +  ~~~~~~~
                TsrfSL(i,j,n)      =  tairSL(i,j)          ! Surface Temperat.
                TvegTV(i,j,n)      =  tairSL(i,j)          ! Vegetat.Temperat.

                issSNo(i,j,n)      =   0                   ! Nb Supr.Ice Layer
                nisSNo(i,j,n)      =   0                   ! Nb      Ice Layer
                snohSN(i,j,n)      =   0.                  ! Snow Buffer Layer
                SWaSNo(i,j,n)      =   0.                  ! Snow Surfic.Water
          END DO
          END DO
          END DO

        END IF

      END IF


C +--Update the SNOW Pack over the Ice Sheet
C +  =======================================

      IF (jhurGE.eq.0.AND.minuGE.eq.0.AND.jsecGE.eq.0)              THEN

       IF(UPDful)                                                   THEN

C +          *******
        CALL UPDsnow
C +          *******

       ELSE
        DO j=1,my
        DO i=1,mx
          IF   (isolSL(i,j).EQ.3)                                   THEN
            DO n=1,nsx
                SnowZZ(nssSNo(i,j,n)+1) = 0.
            DO l=      nssSNo(i,j,n),1,-1
                SnowZZ(l) =       SnowZZ(l+1) 
     .        + dzsSNo(i,j,n, l) *rosSNo(i,j,n, l)
            END DO
                SnowZZ(1) =       SnowZZ(1) / ro_Wat
            IF (SnowZZ(1).LT.0.4)                                   THEN
                rosNEW          = rosSNo(i,j,n, 1) +   5.0 !  + 1.0m Snow
                                                           !=>+ 0.5m *10kg/m3/m
                WEq_SN(i,j,n)   = WEq_SN(i,j,n)    +rosNEW !    1.0m *rosNEW
                zWE0SN(i,j,n)   = zWE0SN(i,j,n)    +rosNEW ! New Initial SWE
                SnowZZ(1)       = dzsSNo(i,j,n, 1) *rosSNo(i,j,n,1)
     .                                             +rosNEW !    1.0m *rosNEW
                rosSNo(i,j,n,1) = rosSNo(i,j,n, 1) *dzsSNo(i,j,n,1)
     .                          + rosNEW                   !    1.0m *rosNEW
                dzsSNo(i,j,n,1) = dzsSNo(i,j,n, 1) +1.0000
                rosSNo(i,j,n,1) = rosSNo(i,j,n, 1) /dzsSNo(i,j,n,1)
                g1sSNo(i,j,n,1) =     99.
                g2sSNo(i,j,n,1) =      3.
                nhsSNo(i,j,n,1) =      2                   !            [-]

                wei0IB(i,j,n)   = wei0IB(i,j,n)    +rosNEW *1.e3 /ro_Wat
            END IF
            END DO
          END IF
        END DO
        END DO

       END IF
      END IF


C +--Initialization of the Temperature Numerical Noise Diagnostics
C +  =============================================================

      IF (iterun.EQ.0)                                              THEN

        DO n=1,2
        DO j=1,my
        DO i=1,mx
          TAir_1(i,j,n) = TairDY(i,j,mz)
          TAir_2(i,j,n) = TairDY(i,j, 2)
        END DO
        END DO
        END DO


C +--OPEN      NetCDF FILE
C +  =====================


C +--ANImation NetCDF File
C +  ---------------------



          wriMIP = .true.
          wriAWS = .true.
          dTiANI = 1800.0     ! Time Interval between IO on ANI.nc
          dTiMIP =86400.0     !                             MIP.nc
          dTiAWS =  360.1     !                             AWS.nc
                              ! (the minimum sampling time step
                              !  when t  [hr since 15 jan 1901]
                              !  has only one decimal digit   )

          include 'SBCnew.dt'

          nboANI = dTiANI          /dt
C +...    nboANI : Number of Time steps between IO
c ###   IF    (mod(dTiANI,dt).NE.0)                                 THEN
c ###              dTiANI = nboANI *dt
c ###   END IF

          iprANI =                   1
          nprANI = nterun / nboANI + 1
C +...    nprANI : Number of                    IO


C +--AMMA-MIP  NetCDF File
C +  ---------------------

              nboMIP = dTiMIP / dt     ! Number of Time steps between IO
          IF (nboMIP.EQ.0)                                          THEN
              nboMIP = 1
              dTiMIP = dt
          END IF

          iprMIP =                   1
          nprMIP = nterun / nboMIP + 1
C +...    nprMIP : Number of                    IO


C +--AWS       NetCDF File
C +  ---------------------

              nboAWS = dTiAWS / dt     ! Number of Time steps between IO
          IF (nboAWS.EQ.0)                                          THEN
              nboAWS = 1
              dTiAWS = dt
          END IF

          iprAWS =                   1
          nprAWS = nterun / nboAWS + 1
C +...    nprAWS : Number of                    IO

C +                        ******
                      call ANI_nc
          IF (wriMIP) call MIP_nc
                      call AWSloc
          IF (wriAWS) call AWS_nc
C +                        ******

      ELSE


C +--Temperature Numerical Noise Diagnostics
C +  =======================================

        DO j=1,my
        DO i=1,mx
          TNois1(i,j)   = abs(0.5* ( TAir_1(i,j,2)+TairDY(i,j,mz) )
     .                       -       TAir_1(i,j,1)                  )
          TAir_1(i,j,2) = TAir_1(i,j, 1)
          TAir_1(i,j,1) = TairDY(i,j,mz)
          TNois2(i,j)   = abs(0.5* ( TAir_2(i,j,2)+TairDY(i,j, 2) )
     .                       -       TAir_2(i,j,1)                  )
          TAir_2(i,j,2) = TAir_2(i,j, 1)
          TAir_2(i,j,1) = TairDY(i,j, 2)
        END DO
        END DO


C +--OUTPUT on NetCDF FILE
C +  =====================

C +--ANImation NetCDF File
C +  ---------------------

        IF       (mod(iterun,nboANI).eq.0)                        THEN
                      iprANI=iprANI  +  1

C +            ******
          call ANI_nc
C +            ******

        END IF


C +--MIP       NetCDF File
C +  ---------------------


        IF                                     (wriMIP)           THEN
C +            ******
          call MIPadd
C +            ******

          IF     (mod(iterun,nboMIP).eq.0)                        THEN
                      iprMIP=iprMIP  +  1

C +            ******
          call MIP_nc
C +            ******

          END IF
        END IF


C +--AWS       NetCDF File
C +  ---------------------

        IF       (mod(iterun,nboAWS).eq.0 .AND. wriAWS)           THEN
                      iprAWS=iprAWS  +  1

C +            ******
          call AWS_nc
C +            ******

        END IF
      END IF


C +--Numerical Validation
C +  ====================

      include "SBCnew.VER"

      return
      end 



C +------------------------------------------------------------------------+
      subroutine UPDsnow
C +------------------------------------------------------------------------+
C | MAR SURFACE  XF                                                    MAR |
C +------------------------------------------------------------------------+
C +
      implicit none
C +
C +--General Variables
C +  =================

      include 'MARphy.inc'
      include 'MARCTR.inc'
      include 'MAR_SV.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc' 
      include 'MAR_GE.inc'
      include 'MAR_DY.inc'
      include 'MAR_LB.inc'
      include 'MAR_SL.inc'
      include 'MAR_SN.inc'
      include 'MAR_BS.inc'
      include 'MAR_IO.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MAR_IB.inc'

C +--Local   Variables
C +  =================

      integer                 n
 
      real                    dz_old

      real      ,parameter :: maxlimit = 15.0 ! Maximum Snow Height
      real      ,parameter :: minlimit =  8.0 ! Minimum Snow Height

      write(*,*)  ' '
      write(*,9) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE
    9 format(' UPDsnow call:',i5,4i3)
      write(*,*)  ' '

C +   2. checking of zn1IB
C +   --------------------
             
      do i=1,mx ; do j=1,my ; if (isolSL(i,j) == 3) then 
      do n=1,max(1,nsx)

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       if (SIh_IB(i,j,n) > maxlimit ) then
       
        dz_old          =                    dzsSNo(i,j,n,1)
        dzsSNo(i,j,n,1) = max(0.10,0.5    *  dzsSNo(i,j,n,1))
         
        wei0IB(i,j,n)   = wei0IB(i,j,n)   +  rosSNo(i,j,n,1)
     .                  *(dzsSNo(i,j,n,1) -  dz_old)  
     .                  * 1.e3 / ro_Wat 

c #WR   write(*,*)  ' '
c #WR   write(*,10) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE,i,j
   10   format(' Update (minus) of snow pack at',
     .           i5,4i3,' for (',i3,','i3,')')
c #WR   write(*,*)  ' '

       end if

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       if (SIh_IB(i,j,n) < minlimit ) then
 
c #WR   write(*,*)  ' '
c #WR   write(*,11) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE,i,j
   11   format(' Update (add) of snow pack at',
     .           i5,4i3,' for (',i3,','i3,')')
c #WR   write(*,*)  ' '

         dz_old          =                    dzsSNo(i,j,n,1)
         dzsSNo(i,j,n,1) = min(4.0,2.0     *  dzsSNo(i,j,n,1))

         wei0IB(i,j,n)   = wei0IB(i,j,n)   +  rosSNo(i,j,n,1)
     .                   *(dzsSNo(i,j,n,1) -  dz_old)
     .                   * 1.e3 / ro_Wat
                  
       end if
      
      enddo; endif; enddo ; enddo
      end subroutine


      subroutine WGustE

C +-------------------------------------------------------------------+
C |                                                                   |
C | MAR GUSTS                                         22-07-2004  MAR |
C |   SubRoutine WGustE computes diagnostic Wind Gust Estimates       |
C |                                                                   |
C +-------------------------------------------------------------------+
C |                                                                   |
C | This routine aims at estimating wind gusts. It also includes the  |
C | computation of a bounding interval around the estimate which aims |
C | to contain with a high probability observed gusts.                |
C |                                                                   |
C | Ref.  : Brasseur O., MWR, 129, 5-25.                              |
C | ^^^^^^^                                                           |
C |                                                                   |
C | Input : - uairDY : U-wind                                         |
C | ^^^^^^^ - vairDY : V-wind                                         |
C |         - wairDY : W-wind                                         |
C |         - tairDY : REAL temperature                               |
C |         - tairSL : surface temperature                            |
C |         - qvDY   : specific humidity                              |
C |         - qwHY   : cloud dropplets                                |
C |         - qiHY   : ice crystals                                   |
C |         - qrHY   : rain                                           |
C |         - qsHY   : snow                                           |
C |         - SLuts  : surface heat flux                              |
C |         - ect_TE : turbulent kinetic energy                       |
C |         - zzDY   : level heights                                  |
C |         - sh     : surface elevation                              |
C |         - sigma  : sigma levels                                   |
C |         - pstDY  : pressure depth                                 |
C |         - ptopDY : pressure at the top of the model               |
C |                                                                   |
C | Output: - SL_wge : gust estimate                          (m/s)   |
C | ^^^^^^^ - SLlwge : lower bound of the bounding interval   (m/s)   |
C |         - SLuwge : upper bound of the bounding interval   (m/s)   |
C |                                                                   |
C +-------------------------------------------------------------------+


       IMPLICIT NONE


C +---General Variables
C +   -----------------

       INCLUDE "MARdim.inc"
       INCLUDE "MARgrd.inc"
       INCLUDE "MAR_DY.inc"
       INCLUDE "MAR_HY.inc"
       INCLUDE "MAR_TE.inc"
       INCLUDE "MAR_SL.inc"
       INCLUDE "MAR_WK.inc"


C +---Local variables
C +   ---------------

       INTEGER l,lev_up(mx,my),lev_dw(mx,my),Top_BL(mx,my),kzi

       REAL    dtkemin   , ra      ,    cp,gravit,   cap
       REAL    coeffmin  , srf_TKE , rzero


C +---Data
C +   ----


       DATA dtkemin  /    1.e-5    /
       DATA coeffmin /    0.01     /
       DATA cap      /    0.28586  /
       DATA ra       /  287.       /
       DATA cp       / 1004.       /
       DATA gravit   /    9.81     /
       DATA rzero    /    0.0      /


C +   ====================================================


C +   ====================================================



C +---Compute Virtual and Equivalent Potential Temperature
C +   ----------------------------------------------------

      DO   k=1,mz
        DO j=1,my
        DO i=1,mx
          WKxyz1(i,j,k)=tairDY(i,j,k)
     .           *(100. / (pstDY(i,j)*sigma(k)+ptopDY))**cap
     .           *(1.+0.608*qvDY(i,j,k)-qwHY(i,j,k)-qiHY(i,j,k)
     .                                 -qrHY(i,j,k)-qsHY(i,j,k))
C +       ^^^ Virtual Potential Temperature


C +---Compute wind norm
C +   -----------------

          WKxyz2(i,j,k)=SQRT(uairDY(i,j,k)*uairDY(i,j,k)
     .                      +vairDY(i,j,k)*vairDY(i,j,k)
     .                      +wairDY(i,j,k)*wairDY(i,j,k)*1.e-4)

        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}
      ENDDO                                                       ! {Loop on k}


C +---Compute Integrated Buoyancy
C +   ---------------------------

           k=  mz
        DO j=1,my
        DO i=1,mx
          WKxyz8(i,j,mzz) =
     .       MAX(0.0,
     .          (0.5*WKxyz1(i,j,k  )+0.5*WKxyz1(i,j,k -1)
     .          -0.5*WKxyz1(i,j,mz )-0.5*WKxyz1(i,j,mz-1))*gravit
     .         /(0.5*WKxyz1(i,j,mz )+0.5*WKxyz1(i,j,mz-1))*0.5
     .         *(    gplvDY(i,j,k-1)-    gplvDY(i,j,k   ))/gravit)! Int. Buoy.
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}

           k=  mz - 1
        DO j=1,my
        DO i=1,mx
          WKxyz8(i,j,mzz) =     
     .       MAX(0.0,WKxyz8(i,j,mzz)
     .         -(0.5*WKxyz1(i,j,k  )+0.5*WKxyz1(i,j,k -1)
     .          -0.5*WKxyz1(i,j,mz )-0.5*WKxyz1(i,j,mz-1))*gravit
     .         /(0.5*WKxyz1(i,j,mz )+0.5*WKxyz1(i,j,mz-1))*0.5
     .         *(    gplvDY(i,j,k-1)-    gplvDY(i,j,k   ))/gravit)! Int. Buoy.
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}

      DO   k=  mz,2,-1
        DO j=1,my
        DO i=1,mx
          WKxyz8(i,j,k) =
     .               WKxyz8(i,j,k+1)
     .         -(0.5*WKxyz1(i,j,k  )+0.5*WKxyz1(i,j,k -1)
     .          -0.5*WKxyz1(i,j,mz )-0.5*WKxyz1(i,j,mz-1))*gravit
     .         /(0.5*WKxyz1(i,j,mz )+0.5*WKxyz1(i,j,mz-1))*0.5
     .         *(    gplvDY(i,j,k-1)-    gplvDY(i,j,k   ))/gravit ! Int. Buoy.
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}
      ENDDO                                                       ! {Loop on k}


C +---Filtering of turbulent kinetic energy
C +   -------------------------------------

        DO j=1,my
        DO i=1,mx
          WKxyz4(i,j, 1) = ect_TE(i,j, 1)                         ! tmpECT TBC
          WKxyz4(i,j,mz) = ect_TE(i,j,mz)                         ! tmpECT SBC
          WKxyz5(i,j, 1) = ect_TE(i,j, 1)                         ! filECT TBC
          WKxyz5(i,j,mz) = ect_TE(i,j,mz)                         ! filECT SBC
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}

      DO   k=   2,mz-1
        DO j=jp11,my1
        DO i=ip11,mx1
          WKxyz4(i,j,k)  =
     .               (4. * ect_TE(i  ,j  ,k  )
     .               +2. * ect_TE(i-1,j  ,k  )+2.*ect_TE(i+1,j  ,k  )
     .               +2. * ect_TE(i  ,j-1,k  )+2.*ect_TE(i  ,j+1,k  )
     .               +1. * ect_TE(i-1,j-1,k  )+1.*ect_TE(i-1,j+1,k  )
     .               +1. * ect_TE(i+1,j-1,k  )+1.*ect_TE(i+1,j+1,k  ))
     .              /16.
           WKxy1(i,j)    =(gplvDY(i  ,j  ,k-1)-   gplvDY(i  ,j  ,k  )) 
     .                    /gravit
           WKxy2(i,j)    =(gplvDY(i  ,j  ,k  )-   gplvDY(i  ,j  ,k+1)) 
     .                    /gravit
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}

        DO j=jp11,my1
        DO i=ip11,mx1
          WKxyz5(i,j,k) = 
     .    0.25*(2.                                     *WKxyz4(i,j,k  )
     .         +2.*(WKxy2(i,j)/(WKxy1(i,j)+WKxy2(i,j)))*WKxyz4(i,j,k-1)
     .         +2.*(WKxy1(i,j)/(WKxy1(i,j)+WKxy2(i,j)))*WKxyz4(i,j,k+1))
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}
      ENDDO                                                       ! {Loop on k}


C +---Determination of mean tke below level k
C +   - - - - - - - - - - - - - - - - - - - -

        DO j=jp11,my1
        DO i=ip11,mx1
          WKxy1(i,j)=0.
          WKxy2(i,j)=0.
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}

      DO k=mz-1,2,-1
        DO j=jp11,my1
        DO i=ip11,mx1
         WKxy1(i,j)=WKxy1(i,j)+(gplvDY(i,j,k)+gplvDY(i,j,k+1))*0.5
     .                        *(gplvDY(i,j,k)-gplvDY(i,j,k+1))
     .                        * WKxyz5(i,j,k)
         WKxy2(i,j)=WKxy2(i,j)+(gplvDY(i,j,k)+gplvDY(i,j,k+1))*0.5
     .                        *(gplvDY(i,j,k)-gplvDY(i,j,k+1))
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}

        DO j=jp11,my1
        DO i=ip11,mx1
          WKxyz6(i,j,k) = WKxy1(i,j) / WKxy2(i,j)                 ! mean_tke
        ENDDO                                                     ! {Loop on i}
        ENDDO                                                     ! {Loop on j}
      ENDDO                                                       ! {Loop on k}


C +---Compute wstar
C +   -------------

c #WW   DO j=jp11,my1
c #WW   DO i=ip11,mx1
c #WW     WKxy3(i,j)=MAX(zero,(-gravit*zi__TE(i,j)/290.
c #WW.                                * SLuts(i,j))**(1./3.))
c #WW   ENDDO
c #WW   ENDDO


C +---Evaluation of Gust Wind Speed
C +   -----------------------------

C +---Initial value
C +   - - - - - - -

        DO j=jp11,my1
        DO i=ip11,mx1
          SL_wge(i,j)=MAX(SL_wge(i,j),WKxyz2(i,j,mz))
          SLuwge(i,j)=MAX(SLuwge(i,j),WKxyz2(i,j,mz))
          SLlwge(i,j)=MAX(SLlwge(i,j),WKxyz2(i,j,mz))

          lev_dw(i,j)=1
          lev_up(i,j)=mz
          Top_BL(i,j)=mz
           WKxy5(i,j)=   (WKxyz5(i,j,mz)+WKxyz5(i,j,mz-1))*0.5
        ENDDO                                                    ! {Loop on i}
        ENDDO                                                    ! {Loop on j}


      DO   k=  mz-1,2,-1

C +---Determination of the Top of Boundary Layer
C +   - - - - - - - - - - - - - - - - - - - - -

        DO j=jp11,my1
        DO i=ip11,mx1
          WKxy4(i,j) =                                           ! TKE_min
     .                (coeffmin                                  ! coeff
     .          + 0.1*((gplvDY(i,j,k )-gplvDY(i,j,mz+1))         !
     .                 /gravit        -2000.            )/1000.) !
     .               *   WKxy5(i,j)
          WKxy6(i,j) = (WKxyz5(i,j,k) +WKxyz5(i,j,k-1 )) *0.5    ! Local TKE
          WKxy7(i,j) =                                           ! ENERGY_low
     .                  2.5/11.       *WKxyz5(i,j,k)             ! Source
     .                                +WKxyz8(i,j,k)             ! Sink
          WKxy8(i,j) =                                           ! Estimate
     .              MAX(WKxyz6(i,j,k) ,WKxyz5(i,j,k)*2.5/11.)    ! Source
     .               +  WKxyz8(i,j,k)                            ! Sink

        ENDDO                                                    ! {Loop on i}
        ENDDO                                                    ! {Loop on j}

        DO j=jp11,my1
        DO i=ip11,mx1
          IF (Top_BL(i,j)    .eq.mz         .and.
     .        WKxyz5(i,j,k+1).gt.WKxy4(i,j) .and.
     .        WKxyz5(i,j,k  ).le.WKxy4(i,j)     ) Top_BL(i,j)=k
        ENDDO                                                    ! {Loop on i}
        ENDDO                                                    ! {Loop on j}


C +---Upper bound on Gust Wind Speed
C +   - - - - - - - - - - - - - - - -

        DO j=jp11,my1
        DO i=ip11,mx1
          IF  ( WKxy6(i,j).gt.WKxy4(i,j)   .and.
     .         Top_BL(i,j).eq.mz           .and.
     .         SLuwge(i,j).lt.WKxyz2(i,j,k)     )                   THEN

               SLuwge(i,j) =  MAX(SLuwge(i,j),WKxyz2(i,j,k))
C +            ^^^            Max Wind Speed  (m/s)
               lev_up(i,j) =  k
C +            ^^^ Level  of  Max Wind Speed  (m/s)

          ENDIF
        ENDDO                                                    ! {Loop on i}
        ENDDO                                                    ! {Loop on j}


C +---Lower bound on Wind Gust
C +   - - - - - - - - - - - -

        DO j=jp11,my1
        DO i=ip11,mx1
          IF (WKxy7(i,j).ge.0.)                          THEN ! ENERGY_low > 0.
             SLlwge(i,j) = MAX(SLlwge(i,j),
     .                         WKxyz2(i,j,k))
C +          ^^^ Min Wind Speed  (m/s)
             lev_dw(i,j) = k
          ENDIF
        ENDDO                                                 ! {Loop on i}
        ENDDO                                                 ! {Loop on j}


C +---Estimate of Wind Gust
C +   - - - - - - - - - - -

        DO j=jp11,my1
        DO i=ip11,mx1
          IF (WKxy8(i,j).ge.0.)                          THEN ! ENERGY_low > 0.
             SL_wge(i,j) = MAX(SL_wge(i,j),
     .                        (WKxyz2(i,j,k)+WKxyz2(i,j,k-1))*0.5)
          ENDIF
        ENDDO                                                 ! {Loop on i}
        ENDDO                                                 ! {Loop on j}
      ENDDO                                                   ! {Loop on k}

        DO j=1,my
        DO i=1,mx
           WKxy1(i,j)     = 0.
           WKxy2(i,j)     = 0.
           WKxy3(i,j)     = 0.
           WKxy4(i,j)     = 0.
           WKxy5(i,j)     = 0.
           WKxy6(i,j)     = 0.
           WKxy7(i,j)     = 0.
           WKxy8(i,j)     = 0.
          WKxyz8(i,j,mzz) = 0.
        ENDDO                                                 ! {Loop on i}
        ENDDO                                                 ! {Loop on j}

      DO   k=1,mz
        DO j=1,my
        DO i=1,mx
          WKxyz1(i,j,k)   = 0.
          WKxyz2(i,j,k)   = 0.
          WKxyz4(i,j,k)   = 0.
          WKxyz5(i,j,k)   = 0.
          WKxyz6(i,j,k)   = 0.
          WKxyz8(i,j,k)   = 0.
        ENDDO                                                 ! {Loop on i}
        ENDDO                                                 ! {Loop on j}
      ENDDO                                                   ! {Loop on k}

C +   ====================================================


C +   ====================================================

      return
      end


      subroutine ANI_nc

C +------------------------------------------------------------------------+
C |                                                                        |
C | MAR OUTPUT                                             05-07-2007  MAR |
C |   SubRoutine ANI_nc is used to write the main Model Variables          |
C |                                      on  a NetCDF file                 |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT: iprANI: Current time step    number                           |
C |   ^^^^^^         (starting from iprANI=1, which => new file creation)  |
C |          nprANI: Total  'time slices' number (max value of iprANI)     |
C |                                                                        |
C |   OUTPUT: NetCDF File adapted to IDL Graphic Software                  |
C |   ^^^^^^                                                               |
C |                                                                        |
C |   CAUTION: 1) This Routine requires the usual NetCDF library,          |
C |   ^^^^^^^^    and the complementary access library  'libUN.a'          |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_DY.inc'
      include 'MAR_TE.inc'
      include 'MAR_TU.inc'

      include 'MAR_HY.inc'
      include 'MAR_RA.inc'
c #TC include 'MAR_TC.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'

      include 'MAR_WK.inc'

      include 'MAR_IO.inc'

      integer           iprANI,nprANI,nboANI
      common/ANI_nc_arg/iprANI,nprANI,nboANI


C +--Local   Variables
C +  =================

      integer    Lfnam,     Ltit,     Luni,     Lnam,     Llnam
      PARAMETER (Lfnam= 40, Ltit= 90, Luni= 31, Lnam= 13, Llnam=50)
C +...Length of char strings 

      CHARACTER*(Lfnam)  fnamNC
      common/ANI_nc_loc/ fnamNC
C +...                   fnamNC: To retain file name.

      real               snow0(mx,my),rain0(mx,my)
      common/Out2rr_loc/ snow0       ,rain0
C +...                   snow0 : Integrated Snow over Previous Time Interval
C +...                   rain0 : Integrated Rain over Previous Time Interval

      integer    NdimNC
      PARAMETER (NdimNC = 5)
C +...Number of defined spatial dimensions (exact)

      integer    MXdim
      PARAMETER (MXdim = 4500)
C +...Maximum Number of all dims: recorded Time Steps
C +   and also maximum of spatial grid points for each direction. 

      integer    MX_var
      PARAMETER (MX_var = 80)
C +...Maximum Number of Variables 

      integer    NattNC
      PARAMETER (NattNC = 2)
C +...Number of REAL attributes given to all variables

      INTEGER           RCODE

      integer           jourNC(MXdim)
      integer           moisNC(MXdim)
      real              yearNC(MXdim)
      real              dateNC(MXdim)
      real              timeNC(MXdim)
      common/OUT2nc_r/  yearNC,dateNC
      real              VALdim(MXdim,0:NdimNC)
      integer           nDFdim(      0:NdimNC)
      common/OUT2nc_d/  nDFdim
      integer           NvatNC(NattNC)
      CHARACTER*(Lnam)  NAMdim(      0:NdimNC)
      CHARACTER*(Luni)  UNIdim(      0:NdimNC)
      CHARACTER*(Lnam)  SdimNC(4,MX_var)       
      CHARACTER*(Luni)  unitNC(MX_var)
      CHARACTER*(Lnam)  nameNC(MX_var)
      CHARACTER*(Llnam) lnamNC(MX_var)
      CHARACTER*(Ltit ) tit_NC
      CHARACTER*(Lnam)  NAMrat(NattNC)
c #TC CHARACTER*9   labelc
      CHARACTER*120 tmpINP

      integer   mz_SBL
      parameter(mz_SBL=1)
      real      ua_SBL(mx,my,mz_SBL),va_SBL(mx,my,mz_SBL)
      real      Ta_SBL(mx,my,mz_SBL),zz_SBL(mx,my,mz_SBL)
      real      qa_SBL(mx,my,mz_SBL)

      real          TNois1(mx,my),TNois2(mx,my)
      real          Tair_1(mx,my,2)
      real          Tair_2(mx,my,2)

      common/TNoise/TNois1       ,TNois2
     .             ,Tair_1       ,Tair_2

      real          thicks0      ,thicks1

      integer   n1000 ,n100a ,n100  ,n10_a ,n10   ,n1    
      integer   m10   ,       jd10  ,jd1
      integer   MMXstp,it    ,mois  ,mill  ,iu
      integer   itotNC,NtotNC,ID__nc
      real      starta(1)
      REAL      starti,DayLen,optwa ,optia ,rhodzk


C +--NetCDF File Initialization
C +  ==========================

      IF (iprANI.eq.1) THEN

          n1000 = 1 +     iyrrGE/1000
          n100a =     mod(iyrrGE,1000)
          n100  = 1 +     n100a /100
          n10_a =     mod(n100a ,100)
          n10   = 1 +     n10_a /10
          n1    = 1 + mod(n10_a ,10)
          m10   = 1 +     mmarGE/10
          m1    = 1 + mod(mmarGE,10)
          jd10  = 1 +     jdarGE/10
          jd1   = 1 + mod(jdarGE,10)


C +--Output File Label
C +  -----------------

        fnamNC = 'ANI.'
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)
     .         // '.' // explIO
     .         // '.nc    '


C +--Output Title
C +  ------------

        tit_NC = 'MAR'
     .         // ' - Exp: ' // explIO
     .         // ' - '
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)


C +--Create File / Write Constants
C +  -----------------------------
        MMXstp = MXdim
C +...  To check array bounds... silently

C +--Time Variable (hour)
C +  ~~~~~~~~~~~~~~~~~~~~

C +...  To define a NetCDF dimension (size, name, unit):
c _UL   nDFdim(0)= nprANI
        nDFdim(0)= 0
        NAMdim(0)= 'time'
        UNIdim(0)= 'HOURS since 1901-01-15 00:00:00'

C +...  Check temporary arrays: large enough ?
        IF (nprANI.gt.MMXstp)
     &  STOP '*** ANI_nc - ERROR : MXdim to low ***'

            starti = jhurGE + minuGE/60.d0 + jsecGE/3600.d0
C +...      starti : Starting Time (= current time in the day)

            starta(1) = (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                      +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                      + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                      + njybGE(mmarGE)            ! (including Leap Day)
     .                  *max(0,1-mod(iyrrGE,4))         !
     .                      + jdarGE     -1      )*  24 !
     .                  +jhurGE                         !
     .                + (minuGE *60 +jsecGE      )/3600.!

        DO it = 1,nprANI
              timeNC(it)   = starti    + (it-1) * nboANI  *dt / 3600.
C +...                                         nboANI: #iter between output

              VALdim(it,0) = starta(1) + (it-1) * nboANI  *dt / 3600.
C +...        VALdim(  ,0) : values of the dimension # 0 (time) 

C +--Time Variable (date)
C +  ~~~~~~~~~~~~~~~~~~~~
              dateNC(it) =          timeNC(it)
              jourNC(it) = jdarGE + timeNC(it) / 24.d0
        END DO
                  mois       =  mmarGE
                  mill       =  iyrrGE
        DO it = 1,nprANI
          IF     (jourNC(it).gt.njmoGE(mois))                     THEN
            DO iu=it,nprANI
                  jourNC(iu) =  jourNC(iu) - njmoGE(mois)
            END DO
                  mois       =  mois + 1
              IF (mois.gt.12)                                     THEN
                  mois       =         1
                  mill       =  mill + 1
              END IF
          END IF
                  moisNC(it) =  mois
                  yearNC(it) =  mill

          IF     (dateNC(it).gt.24.d0-epsi)                       THEN
                  DayLen     =  24.d0
            DO iu=it,nprANI
                  dateNC(iu) = mod(dateNC(iu),DayLen)
            END DO
          END IF
        END DO

        DO it = 1,nprANI
              dateNC(it) =  dateNC(it)
     .             + 1.d+2 *jourNC(it)
     .             + 1.d+4 *moisNC(it)
        END DO


C +--Define horizontal spatial dimensions :    
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +...  Check temporary arrays: large enough ?
        IF (    mx .gt.MMXstp.or.my.gt.MMXstp
     &      .or.mzz.gt.MMXstp.or.mw.gt.MMXstp)
     &    STOP '*** ANI_nc - ERROR : MXdim to low ***'

C +...To define NetCDF dimensions (size, name, unit):

        DO i = 1, mx
          VALdim(i,1) = xxkm(i)
        END DO
          nDFdim(1)   = mx
          NAMdim(1)   = 'x'
          UNIdim(1)   = 'km'

        DO j = 1, my
          VALdim(j,2) = yykm(j)
        END DO
          nDFdim(2)   = my
          NAMdim(2)   = 'y'
          UNIdim(2)   = 'km'

        do k = 1, mz_SBL
          VALdim(k,3) = sigma(mzz-k)
        enddo
          nDFdim(3)   =  mz_SBL
          NAMdim(3)   = 'level'
          UNIdim(3)   = '[sigma]'
C +...    For levels k

        do k = 1, mz_SBL
          VALdim(k,4) = sigmid(mzz-k+1)
        enddo
          nDFdim(4)   =  mz_SBL
          NAMdim(4)   = 'level2'
          UNIdim(4)   = '[sigma]'
C +...    For levels k+1/2

        do k = 1, mw
          VALdim(k,5) = k 
        enddo
          nDFdim(5)   = mw
          NAMdim(5)   = 'sector'
          UNIdim(5)   = '[index]'
C +...    For Surface Sectors 

C +--Variable's Choice (Table MARgou.dat)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        OPEN(unit=10,status='unknown',file='MARgou.dat')

        itotNC = 0
 980    CONTINUE
          READ (10,'(A120)',end=990) tmpINP
          IF (tmpINP(1:4).eq.'    ')                                THEN
            itotNC = itotNC + 1
            READ (tmpINP,'(4x,5A9,A12,A50)')
     &          nameNC(itotNC)  ,SdimNC(1,itotNC),SdimNC(2,itotNC),
     &          SdimNC(3,itotNC),SdimNC(4,itotNC),
     &          unitNC(itotNC)  ,lnamNC(itotNC)
C +...          nameNC: Name
C +             SdimNC: Names of Selected Dimensions (max.4/variable) 
C +             unitNC: Units
C +             lnamNC: Long_name, a description of the variable

c #TC       IF (nameNC(itotNC).eq.'qxTC     '.and.nkWri.ge.1)       THEN
c #TC           nameNC(itotNC) =  namTC(1)
c #TC        IF                                  (nkWri.gt.1)       THEN
c #TC           itot = itotNC
c #TC         DO  n=2,nkWri
c #TC           itot = itot    +  1
c #TC           nameNC(itot)   =  namTC(n)
c #TC          DO m=1,4
c #TC           SdimNC(m,itot) = SdimNC(m,itotNC)
c #TC          END DO
c #TC           unitNC(itot)   = unitNC(itotNC)
c #TC           lnamNC(itot)   = lnamNC(itotNC)
c #TC         END DO
c #TC           itotNC         = itot
c #TC        ENDIF
c #TC       ENDIF

          ENDIF
        GOTO 980
 990    CONTINUE

        CLOSE(unit=10)

        NtotNC = itotNC 
C +...  NtotNC : Total number of variables writen in NetCDF file.

C +--List of NetCDF attributes given to all variables:
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +...  The "actual_range" is the (min,max)
C +     of all data for each variable:
        NAMrat(1) = 'actual_range'
        NvatNC(1) = 2

C +...  The "[var]_range" is NOT of attribute type,
C +     it is a true variable containing the (min,max) for
C +     each level, for 4D (space+time) variables only
C +     (automatic handling by UN library;
C +      must be the LAST attribute)
        NAMrat(NattNC) = '[var]_range'
        NvatNC(NattNC) = 2

C +--Automatic Generation of the NetCDF File Structure
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +     **************
        CALL UNscreate (fnamNC, tit_NC,
     &                  NdimNC, nDFdim, MXdim , NAMdim, UNIdim, VALdim,
     &                  MX_var, NtotNC, nameNC, SdimNC, unitNC, lnamNC,
     &                  NattNC, NAMrat, NvatNC,
     &                  ID__nc) 
C +     **************


C +--Write Time - Constants
C +  ~~~~~~~~~~~~~~~~~~~~~~
        DO j=1,my
        DO i=1,mx
          Wkxy1(i,j) =  GElonh(i,j) * 15.d0
C +...    Conversion: Hour->degrees

          WKxy2(i,j) =  GElatr(i,j) / degrad
C +...    Conversion: rad ->degrees

          WKxy3(i,j) =  isolSL(i,j)
C +...    Conversion to REAL type (integer not allowed)


C +--Crude re-initialization of previous precipitation
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          rain0(i,j) =  rainHY(i,j)
          snow0(i,j) =  snowHY(i,j)
        END DO
        END DO


C +       ************
          CALL UNwrite (ID__nc, 'lon   ', 1  , mx    , my, 1 , Wkxy1)
          CALL UNwrite (ID__nc, 'lat   ', 1  , mx    , my, 1 , Wkxy2)
          CALL UNwrite (ID__nc, 'sh    ', 1  , mx    , my, 1 , sh)
          CALL UNwrite (ID__nc, 'isol  ', 1  , mx    , my, 1 , Wkxy3)
C +       ************


C +--Re-Open file if already created.
C +  ================================


      ELSE

C +     ************
        CALL UNwopen (fnamNC,ID__nc)
C +     ************

      END IF


C +--Write Time-dependent variables:
C +  ===============================


C +--UNLIMITED Time Dimension
C +  ------------------------

      IF (nDFdim(0).eq.0)                         THEN !
          starta(1) =  (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
     .               + (minuGE *60 +jsecGE      )/3600.!

C +     ************
        CALL UNwrite (ID__nc, 'time   ', iprANI,  1,  1,  1, starta(1))
C +     ************

      END IF


C +     ************
        CALL UNwrite (ID__nc, 'date   ',iprANI, 1, 1, 1, dateNC(iprANI))
        CALL UNwrite (ID__nc, 'year   ',iprANI, 1, 1, 1, yearNC(iprANI))
C +     ************


C +--Wind Gust Estimate
C +  ------------------

      DO j=1,my
      DO i=1,mx
         SL_wge(i,j)=0.
         SLlwge(i,j)=0.
         SLuwge(i,j)=0.
      END DO
      END DO

C +           ******
         call WGustE
C +           ******


C +--Cloud Optical Depth
C +  -------------------

      do j=1,my
      do i=1,mx
           optwa  = zero
C +...     optwa  : liquid water path (kg/m2) (droplets)

           optia  = zero
c +...     optia  : liquid water path (kg/m2) (crystals)

         DO k = mzabso+1,mz
           rhodzk = (pstDY(i,j)*sigma(k)+ptopDY)
     .            / (ra*tairDY(i,j,k)*(1.+.608*qvDY(i,j,k)))
     .            * (gpmiDY(i,j,k)-gpmiDY(i,j,k+1))
C +...     rhodzk : (rho / 1000) * (dz * gravit)

           optwa  = optwa + rhodzk * qwHY(i,j,k)
           optia  = optia + rhodzk * qiHY(i,j,k)

         END DO

           WKxy2(i,j)  = 1.5 * ( optwa / 20.d-6
     .                         + optia / 40.d-6 ) *grvinv

      enddo
      enddo


C +--Dynamics, Precipitation
C +  -----------------------

      DO j=1,my
      DO i=1,mx
       DO k=1,mz_SBL
        ua_SBL(i,j,k) = uairDY(i,j,mzz-k)
        va_SBL(i,j,k) = vairDY(i,j,mzz-k)
        Ta_SBL(i,j,k) = tairDY(i,j,mzz-k)
        qa_SBL(i,j,k) =   qvDY(i,j,mzz-k)
        zz_SBL(i,j,k) = gplvDY(i,j,mzz-k)*grvinv
       END DO
         WKxy3(i,j)   = rainHY(i,j) - rain0(i,j)
         rain0(i,j)   = rainHY(i,j)
         WKxy4(i,j)   = snowHY(i,j) - snow0(i,j)
         snow0(i,j)   = snowHY(i,j)
      END DO
      END DO


C +--Snow Pack Vertically Integrated Water
C +  -------------------------------------

      do j=1,my
      do i=1,mx

                WKxy1(i,j)  = 0.
                thicks0     = 0.
                thicks1     = 0.

        k = nsno
 1001   CONTINUE
        k = k -1
        IF (k.EQ.0)                                         GO TO   1000

        IF     (rosSNo(i,j,1,k).GT.0. .AND. dzsSNo(i,j,1,k).GT.0.)  THEN
                thicks1     = thicks1     + dzsSNo(i,j,1,k)
          IF   (thicks1.LT.1.0)                                     THEN
                WKxy1(i,j)  = WKxy1(i,j)  + rosSNo(i,j,1,k)
     .                                     *wasSNo(i,j,1,k)
     .                                     *dzsSNo(i,j,1,k)
          ELSE
            IF (thicks0.LT.1.0)                                     THEN
                WKxy1(i,j)  = WKxy1(i,j)  + rosSNo(i,j,1,k)
     .                                     *wasSNo(i,j,1,k)
     .                                *(1.0-thicks0)     
            END IF
          END IF

                thicks0     = thicks1

        END IF

        GO TO 1001
 1000   CONTINUE

      enddo
      enddo

C +   ************
      CALL UNwrite (ID__nc, 'pstar  ', iprANI, mx, my, 1      , pstDY )
      CALL UNwrite (ID__nc, 'zz_SBL ', iprANI, mx, my, mz_SBL , zz_SBL)
      CALL UNwrite (ID__nc, 'Ta_SBL ', iprANI, mx, my, mz_SBL , Ta_SBL)
      CALL UNwrite (ID__nc, 'qa_SBL ', iprANI, mx, my, mz_SBL , qa_SBL)
      CALL UNwrite (ID__nc, 'ua_SBL ', iprANI, mx, my, mz_SBL , ua_SBL)
      CALL UNwrite (ID__nc, 'va_SBL ', iprANI, mx, my, mz_SBL , va_SBL)
      CALL UNwrite (ID__nc, 'Ta_SRF ', iprANI, mx, my, 1      , TairSL)
      CALL UNwrite (ID__nc, 'H__Sen ', iprANI, mx, my, 1      , hsenSL)
      CALL UNwrite (ID__nc, 'H__Lat ', iprANI, mx, my, 1      , hlatSL)
      CALL UNwrite (ID__nc, 'SL_wge ', iprANI, mx, my, 1      , SL_wge)
      CALL UNwrite (ID__nc, 'SLlwge ', iprANI, mx, my, 1      , SLlwge)
      CALL UNwrite (ID__nc, 'SLuwge ', iprANI, mx, my, 1      , SLuwge)
      CALL UNwrite (ID__nc, 'RadOLR ', iprANI, mx, my, 1      , RAdOLR)
      CALL UNwrite (ID__nc, 'Clouds ', iprANI, mx, my, 1      , cld_SL)
c #EE CALL UNwrite (ID__nc, 'OptDep ', iprANI, mx, my, 1      , RAcdtO)
      CALL UNwrite (ID__nc, 'OptDep ', iprANI, mx, my, 1      , WKxy2 )
      CALL UNwrite (ID__nc, 'RadSol ', iprANI, mx, my, 1      , RAdsol)
      CALL UNwrite (ID__nc, 'Rad_IR ', iprANI, mx, my, 1      , RAD_ir)
      CALL UNwrite (ID__nc, 'Albedo ', iprANI, mx, my, 1      , albeSL)
      CALL UNwrite (ID__nc, 'Rain   ', iprANI, mx, my, 1      , WKxy3 )
      CALL UNwrite (ID__nc, 'Snow   ', iprANI, mx, my, 1      , WKxy4 )
      CALL UNwrite (ID__nc, 'mosaic ', iprANI, mx, my, mw     , SLsrfl)
      CALL UNwrite (ID__nc, 'T_Surf ', iprANI, mx, my, mw     , TsrfSL)
      CALL UNwrite (ID__nc, 'uustar ', iprANI, mx, my, mw     , SLuusl)
      CALL UNwrite (ID__nc, 'utstar ', iprANI, mx, my, mw     , SLutsl)
      CALL UNwrite (ID__nc, 'uqstar ', iprANI, mx, my, mw     , SLuqsl)
      CALL UNwrite (ID__nc, 'usstar ', iprANI, mx, my, mw     , SLussl)
      CALL UNwrite (ID__nc, 'IntROF ', iprANI, mx, my, 1      , runoTV)
      CALL UNwrite (ID__nc, 'IntEVA ', iprANI, mx, my, 1      , evapTV)
      CALL UNwrite (ID__nc, 'VIntWA ', iprANI, mx, my, 1      , WKxy1 )
      CALL UNwrite (ID__nc, 'TNoisS ', iprANI, mx, my, 1      , TNois1)
      CALL UNwrite (ID__nc, 'TNoisT ', iprANI, mx, my, 1      , TNois2)
C +   ************


C +--That 's all, folks: NetCDF File Closure
C +  =======================================

C +   ************
      CALL UNclose (ID__nc)
C +   ************


C +--Work Arrays Reset
C +  =================

      do j=1,my
      do i=1,mx
        WKxy1(i,j)   =0.0
        WKxy2(i,j)   =0.0
        WKxy3(i,j)   =0.0
        WKxy4(i,j)   =0.0
      enddo
      enddo

      return
      end


      subroutine MIPadd

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                             12-03-2007  MAR |
C |   SubRoutine MIPadd is used to update selected Model Variables         |
C |                            for write  on  a NetCDF file                |
C +------------------------------------------------------------------------+


C +--General Variables
C +  =================

      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_DY.inc'
      include 'MAR_TE.inc'
      include 'MAR_TU.inc'

      include 'MAR_HY.inc'
      include 'MAR_RA.inc'
c #TC include 'MAR_TC.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'

      include 'MAR_WK.inc'

      include 'MAR_IO.inc'

      integer           iprMIP,nprMIP,nboMIP
      common/MIP_nc_arg/iprMIP,nprMIP,nboMIP


C +--Variables for daily & zonal averages
C +  ------------------------------------

      integer            mxz
      parameter         (mxz=1)
      real               Ta            ,deglon
      real               pslDAY(mxz,my),pr_DAY(mxz,my),za_DAY(mxz,my,mz)
      real               TmnDAY(mxz,my),TmxDAY(mxz,my),tasDAY(mxz,my)
      real               uasDAY(mxz,my),vasDAY(mxz,my)
      real               husDAY(mxz,my),HL_DAY(mxz,my),HS_DAY(mxz,my)
      real               lwdDAY(mxz,my),lwuDAY(mxz,my),olrDAY(mxz,my)
      real               swdDAY(mxz,my),swuDAY(mxz,my)
      real               ua_DAY(mxz,my,mz)            ,va_DAY(mxz,my,mz)
      real               Ta_DAY(mxz,my,mz)        
      real               qv_DAY(mxz,my,mz)            ,cl_DAY(mxz,my,mz)
      real               nb_itd        ,nb_lon(    my),zonmsk(mx ,my)
      common/MIP_DAY_av/ pslDAY        ,pr_DAY        ,za_DAY
     .                  ,TmnDAY        ,TmxDAY        ,tasDAY    
     .                  ,uasDAY        ,vasDAY    
     .                  ,husDAY        ,HL_DAY        ,HS_DAY    
     .                  ,lwdDAY        ,lwuDAY        ,olrDAY    
     .                  ,swdDAY        ,swuDAY    
     .                  ,ua_DAY                       ,va_DAY       
     .                  ,Ta_DAY       
     .                  ,qv_DAY                       ,cl_DAY       
     .                  ,nb_itd        ,nb_lon        ,zonmsk



C +--Incrementation
C +  --------------

          nb_itd        = nb_itd      +                 1.
      DO j=1,my
          Ta            =                               0.
      DO i=1,mx
          pslDAY(1,j)   = pslDAY(1,j)   +  pstDY(i,j)       *zonmsk(i,j)
          Ta            =     Ta        + tairDY(i,j,mz)    *zonmsk(i,j)
          tasDAY(1,j)   = tasDAY(1,j)   + tairDY(i,j,mz)    *zonmsk(i,j)
          uasDAY(1,j)   = uasDAY(1,j)   + uairDY(i,j,mz)    *zonmsk(i,j)
          vasDAY(1,j)   = vasDAY(1,j)   + vairDY(i,j,mz)    *zonmsk(i,j)
          husDAY(1,j)   = husDAY(1,j)   +   qvDY(i,j,mz)    *zonmsk(i,j)
          HL_DAY(1,j)   = HL_DAY(1,j)   + hlatSL(i,j)       *zonmsk(i,j)
          HS_DAY(1,j)   = HS_DAY(1,j)   + hsenSL(i,j)       *zonmsk(i,j)
          lwdDAY(1,j)   = lwdDAY(1,j)   + RAd_ir(i,j)       *zonmsk(i,j)
          lwuDAY(1,j)   = lwuDAY(1,j)   + firmSL(i,j)       *zonmsk(i,j)
          olrDAY(1,j)   = olrDAY(1,j)   + RAdOLR(i,j)       *zonmsk(i,j)
          swdDAY(1,j)   = swdDAY(1,j)   + RAdsol(i,j)       *zonmsk(i,j)
          swuDAY(1,j)   = swuDAY(1,j)   + RAdsol(i,j)       *zonmsk(i,j)
     .                                  * albeSL(i,j)
        DO k=1,mz
          ua_DAY(1,j,k) = ua_DAY(1,j,k) + uairDY(i,j,mzz-k) *zonmsk(i,j)
          va_DAY(1,j,k) = va_DAY(1,j,k) + vairDY(i,j,mzz-k) *zonmsk(i,j)
          Ta_DAY(1,j,k) = Ta_DAY(1,j,k) + tairDY(i,j,mzz-k) *zonmsk(i,j)
          za_DAY(1,j,k) = za_DAY(1,j,k) + gplvDY(i,j,mzz-k) *zonmsk(i,j)
     .                                  * grvinv
          qv_DAY(1,j,k) = qv_DAY(1,j,k) +   qvDY(i,j,mzz-k) *zonmsk(i,j)
          cl_DAY(1,j,k) = cl_DAY(1,j,k) + cfraHY(i,j,mzz-k) *zonmsk(i,j)
        END DO
      END DO
          Ta          =                   Ta         /       nb_lon(  j)
          TmnDAY(1,j)   = min(TmnDAY(1,j),Ta)
          TmxDAY(1,j)   = max(TmxDAY(1,j),Ta)
      END DO

      return
      end


      subroutine MIP_nc

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                             15-03-2007  MAR |
C |   SubRoutine MIP_nc is used to write selected Model Variables          |
C |                                      on  a NetCDF file                 |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT: iprMIP: Current time step    number                           |
C |   ^^^^^^         (starting from iprANI=1, which => new file creation)  |
C |          nprMIP: Total  'time slices' number (max value of iprMIP)     |
C |                                                                        |
C |   OUTPUT: NetCDF File adapted to IDL Graphic Software                  |
C |   ^^^^^^                                                               |
C |                                                                        |
C |   CAUTION: 1) This Routine requires the usual NetCDF library,          |
C |   ^^^^^^^^    and the complementary access library  'libUN.a'          |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_DY.inc'
      include 'MAR_TE.inc'
      include 'MAR_TU.inc'

      include 'MAR_HY.inc'
      include 'MAR_RA.inc'
c #TC include 'MAR_TC.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'

      include 'MAR_WK.inc'

      include 'MAR_IO.inc'

      integer           iprMI1
      integer           iprMIP,nprMIP,nboMIP
      common/MIP_nc_arg/iprMIP,nprMIP,nboMIP


C +--Variables for daily & zonal averages
C +  ------------------------------------

      integer            mxz
      parameter         (mxz=1)
      real               Ta            ,deglon        ,sh_ZON(mxz,my)
      real                                             solZON(mxz,my)
      real               pslDAY(mxz,my),pr_DAY(mxz,my),za_DAY(mxz,my,mz)
      real               TmnDAY(mxz,my),TmxDAY(mxz,my),tasDAY(mxz,my)
      real               uasDAY(mxz,my),vasDAY(mxz,my)
      real               husDAY(mxz,my),HL_DAY(mxz,my),HS_DAY(mxz,my)
      real               lwdDAY(mxz,my),lwuDAY(mxz,my),olrDAY(mxz,my)
      real               swdDAY(mxz,my),swuDAY(mxz,my)
      real               ua_DAY(mxz,my,mz)            ,va_DAY(mxz,my,mz)
      real               Ta_DAY(mxz,my,mz)        
      real               qv_DAY(mxz,my,mz)            ,cl_DAY(mxz,my,mz)
      real               nb_itd        ,nb_lon(    my),zonmsk(mx ,my)
      common/MIP_DAY_av/ pslDAY        ,pr_DAY        ,za_DAY
     .                  ,TmnDAY        ,TmxDAY        ,tasDAY    
     .                  ,uasDAY        ,vasDAY    
     .                  ,husDAY        ,HL_DAY        ,HS_DAY    
     .                  ,lwdDAY        ,lwuDAY        ,olrDAY    
     .                  ,swdDAY        ,swuDAY    
     .                  ,ua_DAY                       ,va_DAY       
     .                  ,Ta_DAY       
     .                  ,qv_DAY                       ,cl_DAY       
     .                  ,nb_itd        ,nb_lon        ,zonmsk


C +--Local   Variables
C +  =================

      integer    Lfnam,     Ltit,     Luni,     Lnam,     Llnam
      PARAMETER (Lfnam= 40, Ltit= 90, Luni= 31, Lnam= 13, Llnam=50)
C +...Length of char strings 

      CHARACTER*(Lfnam)  fnamNC
      common/MIP_nc_loc/ fnamNC
C +...                   fnamNC: To retain file name.

      real               snow0(mx,my),rain0(mx,my)
      common/MIPorr_loc/ snow0       ,rain0
C +...                   snow0 : Integrated Snow over Previous Time Interval
C +...                   rain0 : Integrated Rain over Previous Time Interval

      integer    NdimNC
      PARAMETER (NdimNC = 5)
C +...Number of defined spatial dimensions (exact)

      integer    MXdim
      PARAMETER (MXdim = 4500)
C +...Maximum Number of all dims: recorded Time Steps
C +   and also maximum of spatial grid points for each direction. 

      integer    MX_var
      PARAMETER (MX_var = 80)
C +...Maximum Number of Variables 

      integer    NattNC
      PARAMETER (NattNC = 2)
C +...Number of REAL attributes given to all variables

      INTEGER           RCODE

      integer           jourNC(MXdim)
      integer           moisNC(MXdim)
      real              yearNC(MXdim)
      real              dateNC(MXdim)
      real              timeNC(MXdim)
      common/OUT3nc_r/  yearNC,dateNC
      real              VALdim(MXdim,0:NdimNC)
      integer           nDFdim(      0:NdimNC)
      common/OUT3nc_d/  nDFdim
      integer           NvatNC(NattNC)
      CHARACTER*(Lnam)  NAMdim(      0:NdimNC)
      CHARACTER*(Luni)  UNIdim(      0:NdimNC)
      CHARACTER*(Lnam)  SdimNC(4,MX_var)       
      CHARACTER*(Luni)  unitNC(MX_var)
      CHARACTER*(Lnam)  nameNC(MX_var)
      CHARACTER*(Llnam) lnamNC(MX_var)
      CHARACTER*(Ltit ) tit_NC
      CHARACTER*(Lnam)  NAMrat(NattNC)
c #TC CHARACTER*9   labelc
      CHARACTER*120 tmpINP

      integer   mz_MIP
      parameter(mz_MIP=mz)
      real      ua_MIP(mx,my,mz_MIP),va_MIP(mx,my,mz_MIP)
      real      Ta_MIP(mx,my,mz_MIP),zz_MIP(mx,my,mz_MIP)

      integer   n1000 ,n100a ,n100  ,n10_a ,n10   ,n1    
      integer   m10   ,       jd10  ,jd1
      integer   MMXstp,it    ,mois  ,mill  ,iu
      integer   itotNC,NtotNC,ID__nc
      real      starta(1)
      REAL      starti,DayLen,optwa ,optia ,rhodzk


C +--NetCDF File Initialization
C +  ==========================

      IF (iprMIP.eq.1) THEN

          n1000 = 1 +     iyrrGE/1000
          n100a =     mod(iyrrGE,1000)
          n100  = 1 +     n100a /100
          n10_a =     mod(n100a ,100)
          n10   = 1 +     n10_a /10
          n1    = 1 + mod(n10_a ,10)
          m10   = 1 +     mmarGE/10
          m1    = 1 + mod(mmarGE,10)
          jd10  = 1 +     jdarGE/10
          jd1   = 1 + mod(jdarGE,10)


C +--Output File Label
C +  -----------------

        fnamNC = 'MIP.'
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)
     .         // '.' // explIO
     .         // '.nc    '


C +--Output Title
C +  ------------

        tit_NC = 'MAR'
     .         // ' - Exp: ' // explIO
     .         // ' - '
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)


C +--Create File / Write Constants
C +  -----------------------------
        MMXstp = MXdim
C +...  To check array bounds... silently

C +--Time Variable (hour)
C +  ~~~~~~~~~~~~~~~~~~~~

C +...  To define a NetCDF dimension (size, name, unit):
c _UL   nDFdim(0)= nprMIP
        nDFdim(0)= 0
        NAMdim(0)= 'time'
c #HH   UNIdim(0)= 'HOURS since 1901-01-15 00:00:00'
        UNIdim(0)=  'DAYS since 2000-01-01 00:00:00'

C +...  Check temporary arrays: large enough ?
        IF (nprMIP.gt.MMXstp)
     &  STOP '*** ANI_nc - ERROR : MXdim to low ***'

            starti = jhurGE + minuGE/60.d0 + jsecGE/3600.d0
C +...      starti : Starting Time (= current time in the day)

c #HH       starta(1) = (351+(iyrrGE  -1902) *365        ! Nb Days before iyrrGE
c #HH.                      +(iyrrGE  -1901) /  4        ! Nb Leap Years
c #HH.                      + njyrGE(mmarGE)             ! Nb Days before mmarGE
c #HH.                      + njybGE(mmarGE)             ! (including Leap Day)
c #HH.                  *max(0,1-mod(iyrrGE,4))          !
c #HH.                      + jdarGE     -1      )*  24  !
c #HH.                  +jhurGE                          !
c #HH.                + (minuGE *60 +jsecGE      )/3600. !

            starta(1) =((351+(iyrrGE  -2001) *365 +15-1  ! Nb Days before iyrrGE
     .                      +(iyrrGE  -2000) /  4        ! Nb Leap Years
     .                      + njyrGE(mmarGE)             ! Nb Days before mmarGE
     .                      + njybGE(mmarGE)             ! (including Leap Day)
     .                  *max(0,1-mod(iyrrGE,4))          !
     .                      + jdarGE     -1      )*24    !
     .                  +jhurGE                  )/24.   !
     .                + (minuGE *60 +jsecGE      )/86.4e3!
     .                - 0.5                              !

        DO it = 1,nprMIP
              timeNC(it)   = starti    + (it-1) * nboMIP  *dt / 3600.
C +...                                         nboMIP: #iter between output

              VALdim(it,0) = starta(1) + (it-1) * nboMIP  *dt / 3600.
C +...        VALdim(  ,0) : values of the dimension # 0 (time) 

C +--Time Variable (date)
C +  ~~~~~~~~~~~~~~~~~~~~
              dateNC(it) =          timeNC(it)
              jourNC(it) = jdarGE + timeNC(it) / 24.d0
        END DO
                  mois       =  mmarGE
                  mill       =  iyrrGE
        DO it = 1,nprMIP
          IF     (jourNC(it).gt.njmoGE(mois))                     THEN
            DO iu=it,nprMIP
                  jourNC(iu) =  jourNC(iu) - njmoGE(mois)
            END DO
                  mois       =  mois + 1
              IF (mois.gt.12)                                     THEN
                  mois       =         1
                  mill       =  mill + 1
              END IF
          END IF
                  moisNC(it) =  mois
                  yearNC(it) =  mill

          IF     (dateNC(it).gt.24.d0-epsi)                       THEN
                  DayLen     =  24.d0
            DO iu=it,nprMIP
                  dateNC(iu) = mod(dateNC(iu),DayLen)
            END DO
          END IF
        END DO

        DO it = 1,nprMIP
              dateNC(it) =  dateNC(it)
     .             + 1.d+2 *jourNC(it)
     .             + 1.d+4 *moisNC(it)
        END DO


C +--Define horizontal spatial dimensions :    
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +...  Check temporary arrays: large enough ?
        IF (    mx .gt.MMXstp.or.my.gt.MMXstp
     &      .or.mzz.gt.MMXstp.or.mw.gt.MMXstp)
     &    STOP '*** ANI_nc - ERROR : MXdim to low ***'

C +...To define NetCDF dimensions (size, name, unit):

        DO i = 1, mxz
          VALdim(i,1) = 0.
        END DO
          nDFdim(1)   = mxz
          NAMdim(1)   = 'x'
c #geo    NAMdim(1)   = 'lon'
          UNIdim(1)   = 'degrees_east'
c #deg    UNIdim(1)   = 'deg'

        DO j = 1, my
c #km     VALdim(j,2) = yykm(j)
          VALdim(j,2) = GElatr(imez,j) / degrad
        END DO
          nDFdim(2)   = my
          NAMdim(2)   = 'y'
c #geo    NAMdim(2)   = 'lat'
c #km     UNIdim(2)   = 'km'
          UNIdim(2)   = 'degrees_north'
c #deg    UNIdim(2)   = 'deg'

        do k = 1, mz_MIP
          VALdim(k,3) = sigma(mzz-k)
        enddo
          nDFdim(3)   =  mz_MIP
          NAMdim(3)   = 'level'
          UNIdim(3)   = '[sigma]'
C +...    For levels k

        do k = 1, mz_MIP
          VALdim(k,4) = sigmid(mzz-k+1)
        enddo
          nDFdim(4)   =  mz_MIP
          NAMdim(4)   = 'level2'
          UNIdim(4)   = '[sigma]'
C +...    For levels k+1/2

        do k = 1, mw
          VALdim(k,5) = k 
        enddo
          nDFdim(5)   = mw
          NAMdim(5)   = 'sector'
          UNIdim(5)   = '[index]'
C +...    For Surface Sectors 

C +--Variable's Choice (Table MIPvou.dat)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        OPEN(unit=10,status='unknown',file='MIPvou.dat')

        itotNC = 0
 980    CONTINUE
          READ (10,'(A120)',end=990) tmpINP
          IF (tmpINP(1:4).eq.'    ')                                THEN
            itotNC = itotNC + 1
            READ (tmpINP,'(4x,5A9,A12,A50)')
     &          nameNC(itotNC)  ,SdimNC(1,itotNC),SdimNC(2,itotNC),
     &          SdimNC(3,itotNC),SdimNC(4,itotNC),
     &          unitNC(itotNC)  ,lnamNC(itotNC)
C +...          nameNC: Name
C +             SdimNC: Names of Selected Dimensions (max.4/variable) 
C +             unitNC: Units
C +             lnamNC: Long_name, a description of the variable

c #TC       IF (nameNC(itotNC).eq.'qxTC     '.and.nkWri.ge.1)       THEN
c #TC           nameNC(itotNC) =  namTC(1)
c #TC        IF                                  (nkWri.gt.1)       THEN
c #TC           itot = itotNC
c #TC         DO  n=2,nkWri
c #TC           itot = itot    +  1
c #TC           nameNC(itot)   =  namTC(n)
c #TC          DO m=1,4
c #TC           SdimNC(m,itot) = SdimNC(m,itotNC)
c #TC          END DO
c #TC           unitNC(itot)   = unitNC(itotNC)
c #TC           lnamNC(itot)   = lnamNC(itotNC)
c #TC         END DO
c #TC           itotNC         = itot
c #TC        ENDIF
c #TC       ENDIF

          ENDIF
        GOTO 980
 990    CONTINUE

        CLOSE(unit=10)

        NtotNC = itotNC 
C +...  NtotNC : Total number of variables writen in NetCDF file.

C +--List of NetCDF attributes given to all variables:
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +...  The "actual_range" is the (min,max)
C +     of all data for each variable:
        NAMrat(1) = 'actual_range'
        NvatNC(1) = 2

C +...  The "[var]_range" is NOT of attribute type,
C +     it is a true variable containing the (min,max) for
C +     each level, for 4D (space+time) variables only
C +     (automatic handling by UN library;
C +      must be the LAST attribute)
        NAMrat(NattNC) = '[var]_range'
        NvatNC(NattNC) = 2

C +--Automatic Generation of the NetCDF File Structure
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +     **************
        CALL UNscreate (fnamNC, tit_NC,
     &                  NdimNC, nDFdim, MXdim , NAMdim, UNIdim, VALdim,
     &                  MX_var, NtotNC, nameNC, SdimNC, unitNC, lnamNC,
     &                  NattNC, NAMrat, NvatNC,
     &                  ID__nc) 
C +     **************


C +--Write Time - Constants
C +  ~~~~~~~~~~~~~~~~~~~~~~
        DO j=1,my
        DO i=1,mx
          Wkxy1(i,j) =  GElonh(i,j) * 15.d0
C +...    Conversion: Hour->degrees

          WKxy2(i,j) =  GElatr(i,j) / degrad
C +...    Conversion: rad ->degrees

          WKxy3(i,j) =  isolSL(i,j)
C +...    Conversion to REAL type (integer not allowed)


C +--Crude re-initialization of previous precipitation
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          rain0(i,j) =  rainHY(i,j)
          snow0(i,j) =  snowHY(i,j)
        END DO
        END DO


C +       ************
c #mxy    CALL UNwrite (ID__nc, 'lon   ', 1  , mx    , my, 1 , Wkxy1)
c #mxy    CALL UNwrite (ID__nc, 'lat   ', 1  , mx    , my, 1 , Wkxy2)
c #mxy    CALL UNwrite (ID__nc, 'sh    ', 1  , mx    , my, 1 , sh)
c #mxy    CALL UNwrite (ID__nc, 'isol  ', 1  , mx    , my, 1 , Wkxy3)
C +       ************


C +--Initialization / Zonal Mask
C +  ---------------------------

              nb_itd      =                             0.
        DO j=1,my
              nb_lon(  j) =                             0.
        DO i=1,mx
              deglon      =        GElonh(i,j)       * 15.
          IF (deglon.GT.180.)      deglon  =  deglon -360.
          IF (deglon.GE.-10. .AND. deglon.LE. 10.)                  THEN
              zonmsk(i,j) =                             1.
              nb_lon(  j) =        nb_lon(  j)       +  1.
          ELSE
              zonmsk(i,j) =                             0
          END IF
        END DO
              nb_lon(  j) =    max(nb_lon(  j),epsi)

              sh_ZON(1,j) =                             0.
              solZON(1,j) =                        isolSL(imez,j)
        DO i=1,mx
              sh_ZON(1,j) = sh_ZON(1,j) + zonmsk(i,j) *sh(i   ,j)
        END DO
              sh_ZON(1,j) = sh_ZON(1,j) / nb_lon(  j)
        END DO

C +     ************
        CALL UNwrite (ID__nc, 'sh     ', 1 , mxz, my, 1 , sh_ZON)  ! Surface Altitude [m]
        CALL UNwrite (ID__nc, 'isol   ', 1 , mxz, my, 1 , solZON)  ! Soil    Type     [-]
C +     ************


C +--Set / Reset
C +  -----------

            nb_itd        =      0.
        DO j=1,my
            pslDAY(1,j)   =      0.
          DO i=1,mx
             rain0(i,j)   = rainHY(i,j)
             snow0(i,j)   = snowHY(i,j)
          END DO
            TmnDAY(1,j)   =   6000.
            TmxDAY(1,j)   =   -300.
            tasDAY(1,j)   =      0.
            uasDAY(1,j)   =      0.
            vasDAY(1,j)   =      0.
            husDAY(1,j)   =      0.
            HL_DAY(1,j)   =      0.
            HS_DAY(1,j)   =      0.
            lwdDAY(1,j)   =      0.
            lwuDAY(1,j)   =      0.
            olrDAY(1,j)   =      0.
            swdDAY(1,j)   =      0.
            swuDAY(1,j)   =      0.
          DO k=1,mz
            ua_DAY(1,j,k) =      0.
            va_DAY(1,j,k) =      0.
            Ta_DAY(1,j,k) =      0.
            za_DAY(1,j,k) =      0.
            qv_DAY(1,j,k) =      0.
            cl_DAY(1,j,k) =      0.
          END DO
        END DO


C +--Re-Open file if already created.
C +  ================================


      ELSE

C +     ************
        CALL UNwopen (fnamNC,ID__nc)
C +     ************


C +--Write Time-dependent variables:
C +  ===============================


C +--UNLIMITED Time Dimension
C +  ------------------------

          iprMI1 = iprMIP - 1

        IF (nDFdim(0).eq.0)                         THEN !
c #HH       starta(1) = (351+(iyrrGE  -1902) *365        ! Nb Days before iyrrGE
c #HH.                      +(iyrrGE  -1901) /  4        ! Nb Leap Years
c #HH.                      + njyrGE(mmarGE)             ! Nb Days before mmarGE
c #HH.                      + njybGE(mmarGE)             ! (including Leap Day)
c #HH.                  *max(0,1-mod(iyrrGE,4))          !
c #HH.                      + jdarGE     -1      )*24    !
c #HH.                  +jhurGE                          !
c #HH.                + (minuGE *60 +jsecGE      )/3600. !

            starta(1) =((351+(iyrrGE  -2001) *365 +15-1  ! Nb Days before iyrrGE
     .                      +(iyrrGE  -2000) /  4        ! Nb Leap Years
     .                      + njyrGE(mmarGE)             ! Nb Days before mmarGE
     .                      + njybGE(mmarGE)             ! (including Leap Day)
     .                  *max(0,1-mod(iyrrGE,4))          !
     .                      + jdarGE     -1      )*24.   !
     .                  +jhurGE                  )/24.   !
     .                + (minuGE *60 +jsecGE      )/86.4e3!
     .                - 0.5                              !

C +       ************
          CALL UNwrite (ID__nc, 'time   ',iprMI1, 1,1,1, starta(1))
C +       ************

        END IF

C +       ************
          CALL UNwrite (ID__nc, 'date   ',iprMI1, 1,1,1, dateNC(iprMIP))
          CALL UNwrite (ID__nc, 'year   ',iprMI1, 1,1,1, yearNC(iprMIP))
C +       ************


C +--Daily Average 
C +  -------------

        DO j=1,my
            pslDAY(1,j)   = pslDAY(1,j) * 1000. /(nb_itd *nb_lon(  j))
            pr_DAY(1,j)   =                  0.
        DO i=1,mx
            pr_DAY(1,j)   = pr_DAY(1,j)
     .    +(rainHY(i,j)   -  rain0(i,j)
     .     +snowHY(i,j)   -  snow0(i,j))* zonmsk(i,j) / 86.4
        END DO
            pr_DAY(1,j)   = pr_DAY(1,j)         /         nb_lon(  j)
            tasDAY(1,j)   = tasDAY(1,j)         /(nb_itd *nb_lon(  j))
            uasDAY(1,j)   = uasDAY(1,j)         /(nb_itd *nb_lon(  j))
            vasDAY(1,j)   = vasDAY(1,j)         /(nb_itd *nb_lon(  j))
            husDAY(1,j)   = husDAY(1,j)         /(nb_itd *nb_lon(  j))
            HL_DAY(1,j)   = HL_DAY(1,j)         /(nb_itd *nb_lon(  j))
            HS_DAY(1,j)   = HS_DAY(1,j)         /(nb_itd *nb_lon(  j))
            lwdDAY(1,j)   = lwdDAY(1,j)         /(nb_itd *nb_lon(  j))
            lwuDAY(1,j)   = lwuDAY(1,j)         /(nb_itd *nb_lon(  j))
            olrDAY(1,j)   = olrDAY(1,j)         /(nb_itd *nb_lon(  j))
            swdDAY(1,j)   = swdDAY(1,j)         /(nb_itd *nb_lon(  j))
            swuDAY(1,j)   = swuDAY(1,j)         /(nb_itd *nb_lon(  j))
          DO k=1,mz
            za_DAY(1,j,k) = za_DAY(1,j,k)       /(nb_itd *nb_lon(  j))
            Ta_DAY(1,j,k) = Ta_DAY(1,j,k)       /(nb_itd *nb_lon(  j))
            ua_DAY(1,j,k) = ua_DAY(1,j,k)       /(nb_itd *nb_lon(  j))
            va_DAY(1,j,k) = va_DAY(1,j,k)       /(nb_itd *nb_lon(  j))
            qv_DAY(1,j,k) = qv_DAY(1,j,k)       /(nb_itd *nb_lon(  j))
            cl_DAY(1,j,k) = cl_DAY(1,j,k)       /(nb_itd *nb_lon(  j))
          END DO
        END DO


C +     ************
        CALL UNwrite (ID__nc, 'psl    ', iprMI1, mxz, my, 1 , pslDAY)  ! Surface Pressure                     [Pa]
        CALL UNwrite (ID__nc, 'pr     ', iprMI1, mxz, my, 1 , pr_DAY)  ! Precipitation                   [kg/m2/s]
        CALL UNwrite (ID__nc, 'tasmin ', iprMI1, mxz, my, 1 , TmnDAY)  ! Minimum daily surface air temperature [K]
        CALL UNwrite (ID__nc, 'tasmax ', iprMI1, mxz, my, 1 , TmxDAY)  ! Maximum daily surface air temperature [K]
        CALL UNwrite (ID__nc, 'tas    ', iprMI1, mxz, my, 1 , TasDAY)  ! Surface Air Temperature               [K]
        CALL UNwrite (ID__nc, 'huss   ', iprMI1, mxz, my, 1 , husDAY)  ! Surface air specific humidity     [kg/kg]
        CALL UNwrite (ID__nc, 'hfls   ', iprMI1, mxz, my, 1 , HL_DAY)  ! Surface Latent   Heat Flux         [W/m2]
        CALL UNwrite (ID__nc, 'hfss   ', iprMI1, mxz, my, 1 , HS_DAY)  ! Surface Sensible Heat Flux         [W/m2]
        CALL UNwrite (ID__nc, 'rlds   ', iprMI1, mxz, my, 1 , lwdDAY)  ! Surface Downwelling LW Radiation   [W/m2]
        CALL UNwrite (ID__nc, 'rlus   ', iprMI1, mxz, my, 1 , lwuDAY)  ! Surface   Upwelling LW Radiation   [W/m2]
        CALL UNwrite (ID__nc, 'rsds   ', iprMI1, mxz, my, 1 , swdDAY)  ! Surface Downwelling SW Radiation   [W/m2]
        CALL UNwrite (ID__nc, 'rsus   ', iprMI1, mxz, my, 1 , swuDAY)  ! Surface   Upwelling SW Radiation   [W/m2]
        CALL UNwrite (ID__nc, 'uas    ', iprMI1, mxz, my, 1 , uasDAY)  ! Zonal      Surface Wind Speed       [m/s]
        CALL UNwrite (ID__nc, 'vas    ', iprMI1, mxz, my, 1 , vasDAY)  ! Meridional Surface Wind Speed       [m/s]
        CALL UNwrite (ID__nc, 'rlut   ', iprMI1, mxz, my, 1 , olrDAY)  ! Outgoing      Longwave Radiation   [W/m2]
        CALL UNwrite (ID__nc, 'za     ', iprMI1, mxz, my, mz, za_DAY)  ! Altitude                              [m]
        CALL UNwrite (ID__nc, 'ta     ', iprMI1, mxz, my, mz, Ta_DAY)  ! Temperature                           [K]
        CALL UNwrite (ID__nc, 'ua     ', iprMI1, mxz, my, mz, ua_DAY)  ! Zonal      Wind Speed               [m/s]
        CALL UNwrite (ID__nc, 'va     ', iprMI1, mxz, my, mz, va_DAY)  ! Meridional Wind Speed               [m/s]
        CALL UNwrite (ID__nc, 'hus    ', iprMI1, mxz, my, mz, qv_DAY)  ! Specific   Humidity               [kg/kg]
        CALL UNwrite (ID__nc, 'cl     ', iprMI1, mxz, my, mz, cl_DAY)  ! Cloud Fraction                        [%]
C +     ************


C +--Set / Reset
C +  -----------

            nb_itd        =      0.
        DO j=1,my
            pslDAY(1,j)   =      0.
          DO i=1,mx
             rain0(i,j)   = rainHY(i,j)
             snow0(i,j)   = snowHY(i,j)
          END DO
            TmnDAY(1,j)   =   6000.
            TmxDAY(1,j)   =   -300.
            tasDAY(1,j)   =      0.
            uasDAY(1,j)   =      0.
            vasDAY(1,j)   =      0.
            husDAY(1,j)   =      0.
            HL_DAY(1,j)   =      0.
            HS_DAY(1,j)   =      0.
            lwdDAY(1,j)   =      0.
            lwuDAY(1,j)   =      0.
            olrDAY(1,j)   =      0.
            swdDAY(1,j)   =      0.
            swuDAY(1,j)   =      0.
          DO k=1,mz
            ua_DAY(1,j,k) =      0.
            va_DAY(1,j,k) =      0.
            Ta_DAY(1,j,k) =      0.
            za_DAY(1,j,k) =      0.
            qv_DAY(1,j,k) =      0.
            cl_DAY(1,j,k) =      0.
          END DO
        END DO


      END IF


C +--That 's all, folks: NetCDF File Closure
C +  =======================================

C +   ************
      CALL UNclose (ID__nc)
C +   ************


C +--Work Arrays Reset
C +  =================

      do j=1,my
      do i=1,mx
        WKxy1(i,j)   =0.0
        WKxy2(i,j)   =0.0
        WKxy3(i,j)   =0.0
        WKxy4(i,j)   =0.0
      enddo
      enddo

      return
      end


      subroutine AWSloc

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                             21-09-2007  MAR |
C |   SubRoutine AWSloc searches AWS and Manned Stations on the Model Grid |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MAR_TV.inc'

      include 'MAR_IO.inc'


C +--AWS
C +  ---

      integer                  n_AWS ,      n     ,      nn
      parameter               (n_AWS=333)
      integer            AWSio(n_AWS),AWS_i(n_AWS),AWS_j(n_AWS)
      integer            nnAWS
      REAL               AWSla(n_AWS),AWSlo(n_AWS),AWS_z(n_AWS)
      REAL               AWS_x(n_AWS),AWS_y(n_AWS),AWSdd
      character*6        AWS_0(n_AWS)
      character*8                     AWS_1(n_AWS),AWS_2(n_AWS)

      common /AWS_nc_INT/AWSio       ,AWS_i       ,AWS_j
     .                  ,nnAWS
      common /AWS_nc_REA/AWSla       ,AWSlo       ,AWS_z
      common /AWS_nc_CH6/AWS_0
      common /AWS_nc_CH8/             AWS_1       ,AWS_2


      integer            AWSij(iptx)  ,ij

      REAL               distmn,Pt_la ,Pt_lo ,dd_la ,dd_lo ,Dista 
      REAL               dshmin,dsh
      integer            i__min,j__min,in_min,jn_min,l     ,kkk
      real               dd_min(4)
      integer            ii_min(4)
      integer            jj_min(4)
      logical            AWSver
      logical            SamAlt
      logical            RddAWS
      logical            MAR_nc_GEO


C +--AWS (from AMRC DATA treatment, stored in MAR_nc.GEO)
C +  ----------------------------------------------------

      character*16 stat
      real         xsta     ,ysta     ,zsta     
      real         xMAR     ,yMAR     
      integer      ista     ,jsta
      real         a_Ny                              ! angle of N direction
                                                     ! from     y direction


C +--DATA
C +  ====

      data      AWSver/.false./    ! OUTPUT of the 4 AWS closest grid pnts
      data      SamAlt/.false./    ! AWS Model Grid Point Refinement SWITCH
      data      RddAWS/.false./
      data  MAR_nc_GEO/.true./



C +--Find the closest Grid Point  in the MAR Domain 
C +  ----------------------------------------------

      IF                                    (MAR_nc_GEO)            THEN
        open( unit=30,status='unknown',file='MAR_nc.GEO')
        rewind     30
        open( unit=31,status='unknown',file='MAR_nc.GEOn')
        rewind     31
        write(    6,205        )
 205    format(/,'MAR-grid AWS-coordinates / MAR_nc.GEO') 
 200    CONTINUE
        read (   30,201,END=202) stat    ,xsta    ,ysta    
     .                                   ,zsta    ,a_Ny    
        write(   31,201        ) stat    ,xsta    ,ysta    
     .                                   ,zsta    ,a_Ny    
 201    format(a16,2f8.2,2f8.1,$)
              xMAR     = 9999.
              yMAR     = 9999.
              distmn   =         mx*mx +my*my
              distmn   = dx*sqrt(distmn)
              xsta     = xsta        * degrad 
          IF (ysta     .LT. 0.)
     .        ysta     = ysta        +    360.
              ysta     = ysta        * degrad 
              i__min   = 0
              j__min   = 0
          DO j=1,my
          DO i=1,mx
              Pt_la    = GElatr(i,j) 
              Pt_lo    = GElonh(i,j)
          IF (Pt_lo    .LT. 0.)
     .        Pt_lo    = Pt_lo       +     24.
              Pt_lo    = Pt_lo       * hourad
              dd_la    = earthr      *              (xsta    -Pt_la)
              dd_lo    = earthr      * cos(Pt_la) * (ysta    -Pt_lo)
              Dista    = sqrt(dd_la*dd_la+dd_lo*dd_lo)
          IF (Dista.lt.distmn)                                      THEN
              distmn   = Dista
              i__min   = i
              j__min   = j
          END IF
          END DO
          END DO


C +--A Grid Point has been found in the MAR Domain
C +  ---------------------------------------------

        write(6,204) stat,i__min,j__min,1.e-3*distmn
 204    format(a16,2i6,f8.1)
        IF (distmn.LT.2.*dx .AND. i__min.GE. 1 .AND. i__min.LE.mx      
     .                      .AND. j__min.GE. 1 .AND. j__min.LE.my)  THEN

C +--(x,y) Coordinates of the AWS (kriging)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 xMAR     =  0.
                 yMAR     =  0.
                 AWSdd    =  0.
            DO j=max(1,j__min-1),min(j__min+1,my)
            DO i=max(1,i__min-1),min(i__min+1,mx)
                 Pt_la    = GElatr(i,j) 
                 Pt_lo    = GElonh(i,j)
             IF (Pt_lo    .LT. 0.)
     .           Pt_lo    = Pt_lo       +     24.
                 Pt_lo    = Pt_lo       * hourad
                 dd_la    = earthr      *              (xsta    -Pt_la)
                 dd_lo    = earthr      * cos(Pt_la) * (ysta    -Pt_lo)
                 Dista    =      max(epsi,sqrt(dd_la*dd_la+dd_lo*dd_lo))
                 xMAR     =  xMAR       + dx*(i-imez) / (Dista * Dista)
                 yMAR     =  yMAR       + dx*(j-jmez) / (Dista * Dista)
                 AWSdd    =  AWSdd      + 1.          / (Dista * Dista)
            ENDDO
            ENDDO
                 xMAR     =  xMAR       * 1.0e-3      /  AWSdd
                 yMAR     =  yMAR       * 1.0e-3      /  AWSdd
        END IF

        write(   31,203        )          xMAR    ,yMAR    
     .                                   ,i__min  ,j__min
     .                                   ,imez    ,jmez  
 203    format(2f8.1,4i4)
        GO TO 200
 202    CONTINUE
        write(    6,206        )
 206    format(1x)

        close(unit=30)
        close(unit=31)

      END IF


C +--Geographic and Grid Point Coordinates for IO
C +  ============================================

      open( unit=30,status='unknown',file='SBCnew.AWS')
      rewind     30
      open( unit=31,status='unknown',file='SBCnew.AWS.JNL')
      rewind     31
      IF (AWSver)                                                   THEN
      open( unit=32,status='unknown',file='SBCnew.AWS.VER')
      rewind     32
      END IF
           write( 6,6000) 
           write(30,6000) 
 6000      format('AWS    Station         | Latit. | Longit.|'
     .     ,                     ' x [km] | y [km] | Altit.||'
     .     ,           ' Grid pt.| x [km] | y [km] |'
     .     ,                     ' Latit. | Longit.| Altit.||'
     .     ,             ' D(AWS)|')
           write( 6,6002) 
           write(30,6002) 
 6002      format('-----------------------+--------+--------+'
     .     ,                     '--------+--------+-------++'
     .     ,           '---------+--------+--------+'
     .     ,                     '--------+--------+-------++'
     .     ,             '-------+')


C +--Find  the 4 closest Grid Points  in MAR Domain 
C +  (AWS having a 0 PRiority are *NOT* considered)
C +  ----------------------------------------------

      DO n=1,n_AWS
        DO k=1,4
              dd_min(k)=         mx*mx +my*my
              dd_min(k)= dx*sqrt(dd_min(k))
              ii_min(k)= 0
              jj_min(k)= 0
        ENDDO
        IF   (AWSio(n) .GT. 0 )                                     THEN
              AWSla(n) = AWSla(n)    * degrad 
          IF (AWSlo(n) .LT. 0.)
     .        AWSlo(n) = AWSlo(n)    +    360.
              AWSlo(n) = AWSlo(n)    * degrad 
          DO j=1,my
          DO i=1,mx
              Pt_la    = GElatr(i,j) 
              Pt_lo    = GElonh(i,j)
          IF (Pt_lo    .LT. 0.)
     .        Pt_lo    = Pt_lo       +     24.
              Pt_lo    = Pt_lo       * hourad
              dd_la    = earthr      *              (AWSla(n)-Pt_la)
              dd_lo    =                             AWSlo(n)-Pt_lo
          IF (dd_lo    .GT. 2.0*pi)
     .        dd_lo    = dd_lo       - 2.0*pi
              dd_lo    = earthr      * cos(AWSla(n))*dd_lo
              Dista    = sqrt(dd_la*dd_la+dd_lo*dd_lo)
          IF (Dista .lt. dd_min(1))                                 THEN
              dd_min(1)= Dista
              ii_min(1)= i
              jj_min(1)= j
          ELSE
     .    IF (Dista .lt. dd_min(2))                                 THEN
              dd_min(2)= Dista
              ii_min(2)= i
              jj_min(2)= j
          ELSE
     .    IF (Dista .lt. dd_min(3))                                 THEN
              dd_min(3)= Dista
              ii_min(3)= i
              jj_min(3)= j
          ELSE
     .    IF (Dista .lt. dd_min(4))                                 THEN
              dd_min(4)= Dista
              ii_min(4)= i
              jj_min(4)= j
          END IF

          END DO
          END DO

              i__min   = ii_min(1)
              j__min   = jj_min(1)
         IF  (i__min.GT.0.AND.j__min.GT.0)                          THEN
             kkk=0 
          DO k=1,4
             kkk=(ii_min(k)-ii_min(1))*(ii_min(k)-ii_min(1))
     .          +(jj_min(k)-jj_min(1))*(jj_min(k)-jj_min(1))
          ENDDO
              distmn   = dd_min(1)
          IF (AWSver)
     .        write(32,6660)  
     .                AWS_0(n),kkk
     .              ,       (ii_min(k),k=1,4)
     .              ,       (jj_min(k),k=1,4)
     .              , AWSla(n)                          /degrad 
     .              ,(GElatr(ii_min(k),jj_min(k))       /degrad,k=1,4)
     .              , AWSlo(n)                          /degrad 
     .              ,(GElonh(ii_min(k),jj_min(k))*hourad/degrad,k=1,4)
     .              ,(1.e-3 *dd_min(k),                         k=1,4)
 6660         format(/,a6,i12,4i12,/,18x,4i12,2(/,13x,5(f8.2,4x))
     .                                          /,25x,5(f8.2,4x))
         END IF
        ELSE
              i__min   = 0
              j__min   = 0
        END IF


C +--A Grid Point has been found in the MAR Domain
C +  ---------------------------------------------

        IF   (i__min.GT. 1 .AND. i__min.LT.mx      
     .   .AND.j__min.GT. 1 .AND. j__min.LT.my)                      THEN

C +--MAR Coordinates of the AWS
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~
              AWS_x(n) = 0.
              AWS_y(n) = 0.
              AWSdd    = 0.
          DO k=1,4
              AWS_x(n) = AWS_x(n) + dx*         (ii_min(k)-imez) 
     .                                /(max(epsi,dd_min(k)*dd_min(k)))
              AWS_y(n) = AWS_y(n) + dx*         (jj_min(k)-jmez)
     .                                /(max(epsi,dd_min(k)*dd_min(k)))
              AWSdd    = AWSdd    + 1./(max(epsi,dd_min(k)*dd_min(k)))
          ENDDO 
              AWS_x(n) = AWS_x(n) * 1.0e-3      /  AWSdd
              AWS_y(n) = AWS_y(n) * 1.0e-3      /  AWSdd

C +--CRITERION - CRITERION - CRITERION - CRITERION - CRITERION - CRITERION 
C +  Determine  a MAR Grid Point with the closest Altitude in the Vicinity
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (SamAlt.AND.abs(sh(i__min,j__min)-AWS_z(n)).GT.100.)   THEN
                 dshmin   = 1.e6
            DO k=1,4
                 i        = ii_min(k)
                 j        = jj_min(k)
                 Pt_la    = GElatr(i,j) 
                 Pt_lo    = GElonh(i,j) * hourad
                 dd_la    = earthr      *              (AWSla(n)-Pt_la)
                 dd_lo    = earthr      * cos(Pt_la) * (AWSlo(n)-Pt_lo)
                 Dista    = sqrt(dd_la*dd_la+dd_lo*dd_lo)
                 dsh      =  abs(sh(i,j)    -AWS_z(n)   )
              IF(isolSL(i,j).GE.3                     .AND.
     .           dsh        .LT.dshmin                .AND.
     .           Dista      .LT.dx*1.5                     )        THEN
                 in_min   = i
                 jn_min   = j
                 dshmin   = dsh
                 distmn   = Dista
              END IF
            ENDDO
                 i__min   = in_min
                 j__min   = jn_min
          ENDIF

              AWS_i(n) = i__min
              AWS_j(n) = j__min
                  i    = i__min
                  j    = j__min

C +--Summary Table
C +  ~~~~~~~~~~~~~
           write( 6,6001) AWS_0(n)
     .                   ,AWS_1(n)       ,AWS_2(n)
     .                   ,AWSla(n)/degrad,AWSlo(n)/degrad
     .                   ,AWS_x(n)       ,AWS_y(n)
     .                   ,AWS_z(n)
     .                   ,       i__min,j__min
     .                   ,   dx*(i__min-imez)     * 1.e-3
     .                   ,   dx*(j__min-jmez)     * 1.e-3
     .                   ,GElatr(i__min,j__min)   /degrad
     .                   ,GElonh(i__min,j__min)   *hourad/degrad
     .                   ,    sh(i__min,j__min),    1.e-3*distmn
           write(30,6001) AWS_0(n)
     .                   ,AWS_1(n)       ,AWS_2(n)
     .                   ,AWSla(n)/degrad,AWSlo(n)/degrad
     .                   ,AWS_x(n)       ,AWS_y(n)
     .                   ,AWS_z(n)
     .                   ,       i__min,j__min
     .                   ,   dx*(i__min-imez)     * 1.e-3
     .                   ,   dx*(j__min-jmez)     * 1.e-3
     .                   ,GElatr(i__min,j__min)   /degrad
     .                   ,GElonh(i__min,j__min)   *hourad/degrad
     .                   ,    sh(i__min,j__min),    1.e-3*distmn
 6001      format(a6,1x,2a8, '|',2(f7.2,' |'),2(f7.1,' |'),f6.0,' ||',
     .                  2i4,' |',2(f7.1,' |'),
     .                           2(f7.2,' |'),  f6.0,' ||',f6.1,' |')
           write(31,3100)  AWS_x(n)
     .                   , AWS_y(n) + 1.e-3*dx*real(4*mx)/real(160)
     .                   , AWS_0(n)
     .                   , AWS_x(n)
     .                   , AWS_y(n)
     .                   ,   dx*(i     -imez)     * 1.e-3
     .                   ,   dx*(j     -jmez)     * 1.e-3
c #OU.                   , AWS_1(n)
c #OU.                   ,   dx*(i__min-imez)     * 1.e-3
c #OU.                   ,   dx*(j__min-jmez)     * 1.e-3
 3100      format('LABEL ',2(f8.2,','),'0,0,.10 ',a6
     .         ,/,'LABEL ',2(f8.2,','),'0,0,.12 @P7x'
     .         ,/,'LABEL ',2(f8.2,','),'0,0,.15 @P5.'
c #OU.         ,/,'sp echo " "           >>dsSWEq.LIST'
c #OU.         ,/,'sp echo "',a8,  'SMB" >>dsSWEq.LIST'
c #OU.         ,/,'set region/x=',f8.2,'/y=',f8.2
c #OU.         ,/,'list/append/nohead/file=dsSWEq.LIST '
c #OU.           ,                        'd_BSWE,dnBSWE'
     .                                                   )


C +--Stations which are not in the Model Domain are discarded
C +  --------------------------------------------------------
        ELSE
               AWSio(n)  = 0
               AWS_i(n)  = 0
               AWS_j(n)  = 0
        END IF


C +--Redondant Grid Points
C +  ---------------------
        
        IF (RddAWS)                                             THEN
          DO nn=1,n-1
            IF ((AWS_i(n).EQ.AWS_i(nn)).AND.
     .          (AWS_j(n).EQ.AWS_j(nn))     )                   THEN
                 AWSio(n)  = 0
                 AWS_i(n)  = 0
                 AWS_j(n)  = 0
            END IF
          ENDDO
        END IF

      END DO


             nnAWS       = 0
      DO n=1,n_AWS
        IF  (AWSio(n).GT.0)                                     THEN
             nnAWS       = nnAWS + 1
             AWSio(nnAWS)= AWSio(n)
             AWS_i(nnAWS)= AWS_i(n)
             AWS_j(nnAWS)= AWS_j(n)
             AWS_x(nnAWS)= AWS_x(n)
             AWS_y(nnAWS)= AWS_y(n)
             AWSla(nnAWS)= AWSla(n)
             AWSlo(nnAWS)= AWSlo(n)
             AWS_z(nnAWS)= AWS_z(n)
             AWS_0(nnAWS)= AWS_0(n)
             AWS_1(nnAWS)= AWS_1(n)
             AWS_2(nnAWS)= AWS_2(n)
        ENDIF
      ENDDO


C +--Modification of the Indices for SISVAT ASCII OUTPUT
C +  ---------------------------------------------------

                ij =       0
      DO n=1,nnAWS
          IF (AWSio (n).EQ.2)                                   THEN
                ij = ij  + 1
            IF (ij.LE.iptx)                                     THEN
              AWSij (ij) =       n
              IOi_TV(ij) = AWS_i(n)
              IOj_TV(ij) = AWS_j(n)
            END IF
          END IF
      ENDDO
      IF (ij.LT.iptx)                                           THEN
        DO n=ij+1,iptx
              AWSij (n)  = AWSij(ij)
              IOi_TV(n)  = AWS_i(ij)
              IOj_TV(n)  = AWS_j(ij)
        ENDDO
      ENDIF


           write( 6,6002) 
           write(30,6002) 

           write( 6,*)  ' '
           write(30,*)  ' '

        DO ij=1,iptx
           write( 6,6003) AWS_1(AWSij(ij)), AWS_2(AWSij(ij))
     .                  ,IOi_TV(      ij ),IOj_TV(      ij )
 6003      format('SISVAT OUTPUT for AWS ',2a8,
     .            ' at Pt (',2i4,')')
        ENDDO
                  
           write( 6,*)  ' '
           write(30,*)  ' '

      close(unit=30)
      close(unit=31)
      IF (AWSver)                                                   THEN
      close(unit=32)
      END IF

      return
      end


      subroutine AWS_nc

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                         Wed 12-03-2008  MAR |
C |   SubRoutine AWS_nc is used to write an AWS Model  Variables           |
C |                                      on a   NetCDF File                |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT: iprAWS: Current time step    number                           |
C |   ^^^^^^         (starting from iprAWS=1, which => new file creation)  |
C |          nprAWS: Total  'time slices' number (max value of iprAWS)     |
C |                                                                        |
C |   OUTPUT: NetCDF File adapted to IDL Graphic Software                  |
C |   ^^^^^^                                                               |
C |                                                                        |
C |   CAUTION: 1) This Routine requires the usual NetCDF library,          |
C |   ^^^^^^^^    and the complementary access library  'libUN.a'          |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_DY.inc'
      include 'MAR_TE.inc'
      include 'MAR_TU.inc'

      include 'MAR_HY.inc'
      include 'MAR_RA.inc'
c #TC include 'MAR_TC.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MARdSV.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MAR_BS.inc'

      include 'MAR_WK.inc'

      include 'MAR_IO.inc'

      integer           iprAWS,nprAWS,nboAWS
      common/AWS_nc_arg/iprAWS,nprAWS,nboAWS


C +--AWS
C +  ---

      integer                  n_AWS,  n
      parameter               (n_AWS=333)
      integer            AWSio(n_AWS),AWS_i(n_AWS),AWS_j(n_AWS)
      integer            nnAWS
      REAL               AWSla(n_AWS),AWSlo(n_AWS),AWS_z(n_AWS)
      character*6        AWS_0(n_AWS)
      character*8                     AWS_1(n_AWS),AWS_2(n_AWS)

      common /AWS_nc_INT/AWSio       ,AWS_i       ,AWS_j
     .                  ,nnAWS
      common /AWS_nc_REA/AWSla       ,AWSlo       ,AWS_z
      common /AWS_nc_CH6/AWS_0
      common /AWS_nc_CH8/             AWS_1       ,AWS_2


C +--Local   Variables
C +  =================

      integer            AWSij(iptx) ,ij

      REAL               distmn,Pt_la ,Pt_lo ,dd_la ,dd_lo ,Dista 
      REAL               dshmin,dsh
      integer            i__min,j__min,in_min,jn_min,l
      logical            SamAlt,INTERP,RESET

      integer            ii           ,jj           ,kzi   ,kmx
      REAL               uu           ,vv           ,deglon,deglat
      REAL               FF(1)        ,DD(1)
      integer            nsnsno, ns
      parameter         (nsnsno= nsno+nsol+1)

      REAL                zzDY(mz)
      REAL               WK__1(1)     ,WK__2(1)     ,WK__3(1)
      REAL               WK__4(1)     ,WK__5(1)     ,WK__6(1)
      REAL               WK_z0(1,1,mz),WK_z1(1,1,mz),WK_z2(1,1,mz)  
      REAL               WK_z3(1,1,mz),WK_z4(1,1,mz),WK_z5(1,1,mz)
      REAL               WK_z6(1,1,mz),WK_z7(1,1,mz),WK_z8(1,1,mz)
      REAL               WK_z9(1,1,mz)
      REAL               WK_n1(1,1,nsnsno),WK_n2(1,1,nsnsno)
      REAL               WK_n3(1,1,nsnsno),WK_n4(1,1,nsnsno)
      REAL               WK_n5(1,1,nsnsno),WK_n6(1,1,nsnsno)
      REAL               WK_n7(1,1,nsnsno),WK_n8(1,1,nsnsno)
      REAL               WK_n9(1,1,nsnsno)

      REAL               VV_SBL,zvvSBL,TT_SBL,zttSBL,Kz_dz
      REAL               LMO   ,zetv_i,zetT_i,Psi__i,Psih_i
      REAL               xx2Psi,xx1Psi,yy2Psi,yy1Psi
      REAL               FF_AWS(1)
      REAL               TT_AWS(1)
      REAL               zvvAWS,zttAWS
      data               zvvAWS,zttAWS/3.0,3.0/
      REAL               ectmin,ectkTE,epskTE
      data               ectmin/1.e-6/


C +--NetCDF  Control Variables
C +  -------------------------

      integer    Lfnam,     Ltit,     Luni,     Lnam,     Llnam
      PARAMETER (Lfnam= 40, Ltit= 90, Luni= 31, Lnam= 13, Llnam=50)
C +...Length of char strings 

      CHARACTER*(Lfnam)  fnamNC
      common/AWS_ncl_ii/ fnamNC
C +...                   fnamNC: To retain file name.

      real               snowb(mx,my),rainb(mx,my)
      common/AWS_ncl_rr/ snowb       ,rainb
C +...                   snowb : Integrated Snow over Previous Time Interval
C +...                   rainb : Integrated Rain over Previous Time Interval

      integer    NdimNC
      PARAMETER (NdimNC = 5)
C +...Number of defined spatial dimensions (exact)

      integer    MXdim
      PARAMETER (MXdim = 15000)
C +...Maximum Number of all dims: recorded Time Steps
C +   and also maximum of spatial grid points for each direction. 

      integer    MX_var
      PARAMETER (MX_var = 80)
C +...Maximum Number of Variables 

      integer    NattNC
      PARAMETER (NattNC = 2)
C +...Number of REAL attributes given to all variables

      INTEGER           RCODE

      integer           jourNC(MXdim)
      integer           moisNC(MXdim)
      real              yearNC(MXdim)
      real              dateNC(MXdim)
      real              timeNC(MXdim)
      common/AWS_nc_r/  yearNC,dateNC
      real              VALdim(MXdim,0:NdimNC)
      integer           nDFdim(      0:NdimNC)
      common/AWS_nc_d/  nDFdim
      integer           NvatNC(NattNC)
      CHARACTER*(Lnam)  NAMdim(      0:NdimNC)
      CHARACTER*(Luni)  UNIdim(      0:NdimNC)
      CHARACTER*(Lnam)  SdimNC(4,MX_var)       
      CHARACTER*(Luni)  unitNC(MX_var)
      CHARACTER*(Lnam)  nameNC(MX_var)
      CHARACTER*(Llnam) lnamNC(MX_var)
      CHARACTER*(Ltit ) tit_NC
      CHARACTER*(Lnam)  NAMrat(NattNC)
c #TC CHARACTER*9   labelc
      CHARACTER*120 tmpINP

      integer   n1000 ,n100a ,n100  ,n10_a ,n10   ,n1    
      integer   m10   ,       jd10  ,jd1
      integer   MMXstp,it    ,mois  ,mill  ,iu
      integer   itotNC,NtotNC,ID__nc
      real      starta(1)
      REAL      starti,DayLen,optwa ,optia ,rhodzk


C +--NetCDF File Initialization
C +  ==========================

C +--Common      Initialization
C +  --------------------------

      IF (iprAWS.eq.1)                                              THEN

C +--Output File Label
C +  ~~~~~~~~~~~~~~~~~
        n1000  = 1 +     iyrrGE/1000
        n100a  =     mod(iyrrGE,1000)
        n100   = 1 +     n100a /100
        n10_a  =     mod(n100a ,100)
        n10    = 1 +     n10_a /10
        n1     = 1 + mod(n10_a ,10)
        m10    = 1 +     mmarGE/10
        m1     = 1 + mod(mmarGE,10)
        jd10   = 1 +     jdarGE/10
        jd1    = 1 + mod(jdarGE,10)

        fnamNC = 'AWS___.'
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)
     .         // '.' // explIO
     .         // '.nc    '

C +--Variable's Choice (Table AWSvar.dat)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OPEN(unit=10,status='unknown',file='AWSvou.dat')

        itotNC = 0
 980    CONTINUE
          READ (10,'(A120)',end=990) tmpINP
          IF (tmpINP(1:4).eq.'    ')                                THEN
            itotNC = itotNC + 1
            READ (tmpINP,'(4x,5A9,A12,A50)')
     &          nameNC(itotNC)  ,SdimNC(1,itotNC),SdimNC(2,itotNC),
     &          SdimNC(3,itotNC),SdimNC(4,itotNC),
     &          unitNC(itotNC)  ,lnamNC(itotNC)
C +...          nameNC: Name
C +             SdimNC: Names of Selected Dimensions (max.4/variable) 
C +             unitNC: Units
C +             lnamNC: Long_name, a description of the variable

          ENDIF
        GOTO 980
 990    CONTINUE

        CLOSE(unit=10)

        NtotNC = itotNC 
C +...  NtotNC : Total number of variables writen in NetCDF file.

C +--List of NetCDF Attributes given to all Variables
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +...  The "actual_range" is the (min,max) of all data for each variable:
        NAMrat(1) = 'actual_range'
        NvatNC(1) = 2

C +...  The "[var]_range"  is NOT of attribute type,
C +     it is a true variable containing the (min,max)
C +     for each level, for 4D (space+time) variables only
C +     (automatic handling by UN library; must be the LAST attribute)
        NAMrat(NattNC) = '[var]_range'
        NvatNC(NattNC) = 2

C +--Array Bounds Check Variable
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MMXstp = MXdim

C +--Time Variable (Hour)
C +  ~~~~~~~~~~~~~~~~~~~~
C +...Temporary Arrays Check: large enough?
        IF (nprAWS.gt.MMXstp)
     &  STOP '*** AWS_nc - ERROR : MXdim to low ***'

C +...NetCDF Dimensions (Size, Name, Unit):
c _UL   nDFdim(0)= nprAWS
        nDFdim(0)= 0
        NAMdim(0)= 'time'
        UNIdim(0)= 'HOURS since 1901-01-15 00:00:00'

           starti = jhurGE + minuGE/60.d0 + jsecGE/3600.d0
C +...     starti : Starting Time (= current time in the day)

           starta(1) = (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
     .               + (minuGE *60 +jsecGE      )/3600.!

        DO it = 1,nprAWS
              timeNC(it)   = starti    + (it-1) * nboAWS  *dt / 3600.
C +...                                         nboAWS: #iter between output

              VALdim(it,0) = starta(1) + (it-1) * nboAWS  *dt / 3600.
C +...        VALdim(  ,0) : values of the dimension # 0 (time) 

C +--Time Variable (Date)
C +  ~~~~~~~~~~~~~~~~~~~~
              dateNC(it) =          timeNC(it)
              jourNC(it) = jdarGE + timeNC(it) / 24.d0
        END DO
                  mois       =  mmarGE
                  mill       =  iyrrGE
        DO it = 1,nprAWS
          IF     (jourNC(it).gt.njmoGE(mois))                     THEN
            DO iu=it,nprAWS
                  jourNC(iu) =  jourNC(iu) - njmoGE(mois)
            END DO
                  mois       =  mois + 1
              IF (mois.gt.12)                                     THEN
                  mois       =         1
                  mill       =  mill + 1
              END IF
          END IF
                  moisNC(it) =  mois
                  yearNC(it) =  mill

          IF     (dateNC(it).gt.24.d0-epsi)                       THEN
                  DayLen     =  24.d0
            DO iu=it,nprAWS
                  dateNC(iu) = mod(dateNC(iu),DayLen)
            END DO
          END IF
        END DO

        DO it = 1,nprAWS
              dateNC(it) =  dateNC(it)
     .             + 1.d+2 *jourNC(it)
     .             + 1.d+4 *moisNC(it)
        END DO


C +--Vertical Variables
C +  ~~~~~~~~~~~~~~~~~~
C +...Temporary Arrays Check: large enough?
        IF (nsnsno.gt.MMXstp.OR.mz.gt.MMXstp
     .                      .OR.mw.gt.MMXstp)   
     .  STOP '*** AWS_nc - ERROR : MXdim to low ***'

C +...NetCDF Dimensions (Size, Name, Unit):
          nDFdim(1)   =  1
          NAMdim(1)   = 'x'
          UNIdim(1)   = 'deg'

          nDFdim(2)   =  1   
          NAMdim(2)   = 'y'
          UNIdim(2)   = 'deg'

        do k = 1, mz
          VALdim(k,3) =  sigma(k)
        enddo
          nDFdim(3)   =  mz
          NAMdim(3)   = 'level'
          UNIdim(3)   = '[index]'
C +...    For atmospheric layers k

        do k = 1, nsnsno
          VALdim(k,4) =  k
        enddo
          nDFdim(4)   =  nsnsno
          NAMdim(4)   = 'level2'
          UNIdim(4)   = '[index]'
C +...    For soil        layers k

        do k = 1, mw
          VALdim(k,5) = k 
        enddo
          nDFdim(5)   =  mw
          NAMdim(5)   = 'sector'
          UNIdim(5)   = '[index]'
C +...    For Surface Sectors 

      END IF


C +--Specific    Initialization
C +  --------------------------

      DO n = 1,nnAWS

          ii          = AWS_i(n)
          jj          = AWS_j(n)

          fnamNC(1:6) = AWS_0(n)

        IF (iprAWS.eq.1)                                            THEN

          tit_NC      = '(A)WS ' // AWS_1(n) // AWS_2(n) // ' / '
     .               // ' OUTPUT of MAR (Modele Atmospherique Regional)'
     .               // ' / Experiment ' // explIO //  '  '
C +...                   1234567890123456789012345678901234567890123456
c #VER    write(6,600) tit_NC
  600     format(/,a90,/,9('1234567890'))

          k=1
          l=  Ltit
 601      CONTINUE
              IF      (tit_NC(k  :k+1).EQ.'  ')                     THEN
                       tit_NC(k  :l-1) =  tit_NC(k+1:l)
                       tit_NC(l  :l  ) =  ' '
                              k  =k-1
                              l  =l-1
              END IF
                              k  =k+1
          IF (k.LT.l)                                          GO TO 601

c #VER    write(6,600) tit_NC

          VALdim(1,1) = AWSlo(n)/degrad
          VALdim(1,2) = AWSla(n)/degrad

C +--Automatic Generation of the NetCDF File Structure
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +       **************
          CALL UNscreate(fnamNC, tit_NC,
     &                   NdimNC, nDFdim, MXdim , NAMdim, UNIdim, VALdim,
     &                   MX_var, NtotNC, nameNC, SdimNC, unitNC, lnamNC,
     &                   NattNC, NAMrat, NvatNC,
     &                   ID__nc) 
C +       **************


C +--Write Time and Geographic Coordinates
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          WK__1 =  GElonh(ii,jj) * 15.000
C +...    Conversion: Hour->degrees

          WK__2 =  GElatr(ii,jj) / degrad
C +...    Conversion: rad ->degrees

          WK__3 =      sh(ii,jj)

          WK__4 =  isolSL(ii,jj)
C +...    Conversion to REAL type (integer not allowed)

C +       ************
          CALL UNwrite (ID__nc, 'sh_AWS', 1  , 1  , 1 , 1 , AWS_z(n))
          CALL UNwrite (ID__nc, 'lonMAR', 1  , 1  , 1 , 1 , WK__1)
          CALL UNwrite (ID__nc, 'latMAR', 1  , 1  , 1 , 1 , WK__2)
          CALL UNwrite (ID__nc, 'sh_MAR', 1  , 1  , 1 , 1 , WK__3)
          CALL UNwrite (ID__nc, 'solTyp', 1  , 1  , 1 , 1 , WK__4)
C +       ************


C +--Re-Open file if already created.
C +  ================================


        ELSE

C +       ************
          CALL UNwopen (fnamNC,ID__nc)
C +       ************

        END IF


C +--Write Time-dependent variables:
C +  ===============================

C +--UNLIMITED Time Dimension
C +  ------------------------

        IF (nDFdim(0).eq.0)                         THEN !
            starta(1) =  (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                       +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                       + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                       + njybGE(mmarGE)            ! (including Leap Day)
     .                   *max(0,1-mod(iyrrGE,4))         !
     .                       + jdarGE     -1      )*  24 !
     .                   +jhurGE                         !
     .                 + (minuGE *60 +jsecGE      )/3600.!

C +       ************
          CALL UNwrite (ID__nc, 'time   ',iprAWS, 1, 1, 1, starta(1))
C +       ************

        END IF


C +     ************
        CALL UNwrite (ID__nc, 'date   ',iprAWS, 1, 1, 1, dateNC(iprAWS))
        CALL UNwrite (ID__nc, 'year   ',iprAWS, 1, 1, 1, yearNC(iprAWS))
C +     ************


C +--Radiative Transfert
C +  -------------------

          WK__1    = RAdOLR(ii,jj)
          WK__2    = RAd_ir(ii,jj)
          WK__3    = RAdsol(ii,jj)
          WK__4    = RAcdtO(ii,jj)
          WK__5    = cld_SL(ii,jj)
          WK__6    = RAertO(ii,jj)

C +     ************
        CALL UNwrite (ID__nc, 'OLR    ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'DownLW ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'DownSW ', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'CloudOD', iprAWS, 1, 1, 1      , WK__4)
        CALL UNwrite (ID__nc, 'CloudFR', iprAWS, 1, 1, 1      , WK__5)
        CALL UNwrite (ID__nc, 'AerosOD', iprAWS, 1, 1, 1      , WK__6)
C +     ************


C +--Dynamics, Microphysics, Turbulence
C +  ----------------------------------

        DO k=1,mz
          WK_z0(1,1,k) = gplvDY(ii,jj,k)*grvinv
          WK_z1(1,1,k) = uairDY(ii,jj,k)
          WK_z2(1,1,k) = vairDY(ii,jj,k)
          WK_z3(1,1,k) = tairDY(ii,jj,k)
          WK_z4(1,1,k) =   qvDY(ii,jj,k)         *1000.
          WK_z5(1,1,k) =   qiHY(ii,jj,k)         *1000.
          WK_z6(1,1,k) =   qwHY(ii,jj,k)         *1000.
          WK_z7(1,1,k) =   qsHY(ii,jj,k)         *1000.
          WK_z8(1,1,k) =   qrHY(ii,jj,k)         *1000.
          WK_z9(1,1,k) = ccniHY(ii,jj,k)
        END DO
          uu       =  WK_z1(1,1,mz)
          vv       =  WK_z2(1,1,mz)


C +     ************
c #VER         write (6,*)    'ZZ'
        CALL UNwrite (ID__nc, 'ZZ     ', iprAWS,  1, 1,mz     , WK_z0)
c #VER         write (6,*)    'UU'
        CALL UNwrite (ID__nc, 'UU     ', iprAWS,  1, 1,mz     , WK_z1)
c #VER         write (6,*)    'VV'
        CALL UNwrite (ID__nc, 'VV     ', iprAWS,  1, 1,mz     , WK_z2)
        CALL UNwrite (ID__nc, 'TT     ', iprAWS,  1, 1,mz     , WK_z3)
        CALL UNwrite (ID__nc, 'QQ     ', iprAWS,  1, 1,mz     , WK_z4)
        CALL UNwrite (ID__nc, 'QI     ', iprAWS,  1, 1,mz     , WK_z5)
        CALL UNwrite (ID__nc, 'QW     ', iprAWS,  1, 1,mz     , WK_z6)
        CALL UNwrite (ID__nc, 'QS     ', iprAWS,  1, 1,mz     , WK_z7)
        CALL UNwrite (ID__nc, 'QR     ', iprAWS,  1, 1,mz     , WK_z8)
        CALL UNwrite (ID__nc, 'ICN    ', iprAWS,  1, 1,mz     , WK_z9)
C +     ************


C +--Turbulence
C +  ----------

        DO k=1,mz
          WK_z4(1,1,k) = ect_TE(ii,jj,k)
          WK_z5(1,1,k) = eps_TE(ii,jj,k)
          WK_z0(1,1,k) = WK_z0(1,1,k) - sh(ii,jj)
        END DO

C +     ************
        CALL UNwrite (ID__nc, 'TKE    ', iprAWS,  1, 1,mz     , WK_z4)
        CALL UNwrite (ID__nc, 'epsTKE ', iprAWS,  1, 1,mz     , WK_z5)
C +     ************


C +--Seeing
C +  ------

C +Height of the PSL and SSL (Primary and Secondary Seeing Layer)
C +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +     ***********
        CALL PBLtop (WK_z4,WK_z0,WK__3,WK__4)
C +     ***********

C +     ************
        CALL UNwrite (ID__nc, 'h__PSL ', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'h__SSL ', iprAWS, 1, 1, 1      , WK__4)
C +     ************

C +CN2 (Refractive Index Structure Coefficient)
C +~~~ (Masciadri & al. 1999, Astron.Astrophys.Suppl.Ser. 137, 185-202)

        DO k=1,mz
          WKxyz1(ii,jj,k)   = pktaDY(ii,jj,k)*3.7222
          WKxyz6(ii,jj,k)   =(uairDY(ii,jj,k)*uairDY(ii,jj,k)            !  V^(5/3)
     .                      + vairDY(ii,jj,k)*vairDY(ii,jj,k))**(5./6.)  !
            zzDY(      k)   = gplvDY(ii,jj,k)*grvinv
          WKxyz8(ii,jj,k)   =   zzDY(      k)
        END DO
          WKxyz6(ii,jj,mzz) =      0.
          WKxyz8(ii,jj,mzz) =     sh(ii,jj)

        DO k=1,mz1
          ectkTE = max(ectmin,ect_TE(ii,jj,k))
          epskTE = max(1.e-08,eps_TE(ii,jj,k))
          WKxyz5(ii,jj,k) =   ect_TE(ii,jj,k)
          WKxyz2(ii,jj,k) =(((WKxyz1(ii,jj,k)-WKxyz1(ii,jj,k+1))
     .                      /(WKxyz8(ii,jj,k)-WKxyz8(ii,jj,k+1)))**2)
     .              *(0.066 * ectkTE         *ectkTE         
     .                      / epskTE                         )
     .              * 1.6   /(epskTE          ** 0.33        )
     .              *(8.e-4 *  pstDY(ii,jj)  *(sigma(k)+sigma(k+1))*0.5
     .               *2.0   /(tairDY(ii,jj,k)+tairDY(ii,jj,k+1))**2)**2
          WKxyz3(ii,jj,k) =   WKxyz2(ii,jj,k)
     .                      *(WKxyz8(ii,jj,k)-WKxyz8(ii,jj,k+1))
        END DO
            
          WKxyz7(ii,jj,1)=(3.5 *zzDY(      1)-  zzDY(      2) *0.5 )*0.5
        DO k=2,mz
          WKxyz7(ii,jj,k) =   ( zzDY(      k)+  zzDY(      k-1)    )*0.5
        END DO
          WKxyz7(ii,jj,mzz) = ( zzDY(     mz)+    sh(ii,jj)        )*0.5

           k=  mz
          WKxyz2(ii,jj,k) =   WKxyz2(ii,jj,k-1)
          WKxyz3(ii,jj,k) =   WKxyz2(ii,jj,k)
     .                      *(WKxyz8(ii,jj,k)-WKxyz8(ii,jj,k+1))

C +CN2 vertical integral
C +~~~~~~~~~~~~~~~~~~~~~

          WKxy1(ii,jj) = 0.
          WKxy2(ii,jj) = 0.
          WKxy3(ii,jj) = 0.
          WKxy4(ii,jj) = 0.
        DO k=mzabso,mz
          WKxyz4(ii,jj,k)=0.
        ENDDO

        DO k=mzabso ,mz
          WKxyz4(ii,jj,k) = WKxyz4(ii,jj,k-1)+ WKxyz3(ii,jj,k)
        ENDDO

C +Seeing at 20 m above the surface
C +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         kmx=mzabso
        DO k=mzabso+1,mz
          IF (zzDY(      k) -   sh(ii,jj) .GT. 20.)   kmx = k
        ENDDO
        DO k=mzabso,kmx,-1
            WKxy1(ii,jj) = WKxy1(ii,jj) + WKxyz3(ii,jj,k)
        ENDDO
            k=      kmx
            WKxy1(ii,jj) = WKxy1(ii,jj) + WKxyz2(ii,jj,k)
     .                                  *(WKxyz8(ii,jj,k)-sh(ii,jj)-20.)

C +Seeing at 30 m above the surface
C +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         kmx=mzabso
        DO k=mzabso+1,mz
          IF (zzDY(      k) -   sh(ii,jj) .GT. 30.)   kmx = k
        ENDDO
        DO k=mzabso,kmx,-1
            WKxy2(ii,jj) = WKxy2(ii,jj) + WKxyz3(ii,jj,k)
        ENDDO
            k=      kmx
            WKxy2(ii,jj) = WKxy2(ii,jj) + WKxyz2(ii,jj,k)
     .                                  *(WKxyz8(ii,jj,k)-sh(ii,jj)-30.)

C +CN2*V^(5/3) vertical integral
C +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             WKxy6(ii,jj) = 0.
        DO k=     1,mz
             WKxy6(ii,jj) = WKxy6(ii,jj) 
     .                   + WKxyz3(ii,jj,k)   * WKxyz6(ii,jj,k)
        ENDDO

C +Fried Parameter/Seeing
C +~~~~~~~~~~~~~~~~~~~~~~
           IF (WKxy1(ii,jj).GT.0.)                                  THEN
            WKxy1(ii,jj) = (1.e-12/(0.423*157.904356*WKxy1(ii,jj)))**0.6
            WKxy1(ii,jj) = (180.*3600./3.1416)*0.49e-6/WKxy1(ii,jj)
           ENDIF

           IF (WKxy2(ii,jj).GT.0.)                                  THEN
            WKxy2(ii,jj) = (1.e-12/(0.423*157.904356*WKxy2(ii,jj)))**0.6
            WKxy2(ii,jj) = (180.*3600./3.1416)*0.49e-6/WKxy2(ii,jj)
           ENDIF

           IF (WKxy4(ii,jj).GT.0.)                                  THEN
            WKxy4(ii,jj) = (1.e-12/(0.423*157.904356*WKxy4(ii,jj)))**0.6
            WKxy4(ii,jj) = (180.*3600./3.1416)*0.49e-6/WKxy4(ii,jj)
           ENDIF

           DO k=mzabso,mz
           IF (WKxyz4(ii,jj,k).GT.0.)                               THEN
           WKxyz4(ii,jj,k)=(1.e-12/(0.423*157.904356*WKxyz4(ii,jj,k)))
     .                                                             **0.6
           WKxyz4(ii,jj,k)=(180.*3600./3.1416)*0.49e-6/WKxyz4(ii,jj,k)
           ENDIF
           ENDDO


        DO k=1,mz
          WK_z4(1,1,k) = WKxyz4(ii,jj,k)
        END DO

C +     ************
        CALL UNwrite (ID__nc, 'Seeing ', iprAWS,  1, 1,mz     , WK_z4)
C +     ************

          WK__1    =  WKxy1(ii,jj)
          WK__2    =  WKxy2(ii,jj)
          WK__4    =  57.*((0.5e-6)**(6./5.))/(WKxy6(ii,jj)**0.6)

C +     ************
        CALL UNwrite (ID__nc, 'S__20m ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'S__30m ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'CoherT ', iprAWS, 1, 1, 1      , WK__4)
C +     ************


C +--Precipitation
C +  -------------

          WK__1    =(snowHY(ii,jj) - snowb (ii,jj)) *1000.
          WK__2    =(rainHY(ii,jj) - rainb (ii,jj)) *1000.

C +     ************
        CALL UNwrite (ID__nc, 'Snow   ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'Rain   ', iprAWS, 1, 1, 1      , WK__2)
C +     ************


C +--SBL
C +  ---

C +--Momentum Turbulent Fluxes
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~
           k=  mz
          WK_z4(1,1,k) =  TUkvm(ii,jj,k)
          Kz_dz        = -TUkvm(ii,jj,k) /(WK_z0(1,1,mz)- sh(ii,jj))
          WK_z5(1,1,k) =  Kz_dz          * WK_z1(1,1,mz)
          WK_z6(1,1,k) =  Kz_dz          * WK_z2(1,1,mz)
        DO k=1,mz1
          WK_z4(1,1,k) =  TUkvm(ii,jj,k)
          Kz_dz        = -TUkvm(ii,jj,k) /(WK_z0(1,1, k)-WK_z0(1,1,k+1))
          WK_z5(1,1,k) =  Kz_dz          *(WK_z1(1,1, k)-WK_z1(1,1,k+1))
          WK_z6(1,1,k) =  Kz_dz          *(WK_z2(1,1, k)-WK_z2(1,1,k+1))
        ENDDO

C +     ************
        CALL UNwrite (ID__nc, 'Kz     ', iprAWS,  1, 1,mz     , WK_z4)
        CALL UNwrite (ID__nc, 'UpWp   ', iprAWS,  1, 1,mz     , WK_z5)
        CALL UNwrite (ID__nc, 'VpWp   ', iprAWS,  1, 1,mz     , WK_z6)
C +     ************

C +--Wind Gusts
C +  ~~~~~~~~~~
          WK__1    = SL_wge(ii,jj) 
          WK__2    = SLlwge(ii,jj) 
          WK__3    = SLuwge(ii,jj)

C +     ************
        CALL UNwrite (ID__nc, 'WGE    ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'WGE_LB ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'WGE_UB ', iprAWS, 1, 1, 1      , WK__3)
C +     ************

C +--3m Meteorological Variables
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VV_SBL=              ssvSL(ii,jj,mz)                     ! Refer. Level Wind
          zvvSBL=    grvinv * gplvDY(ii,jj,mz)-    sh(ii,jj)       ! Refer. Level Height
          TT_SBL=             tairDY(ii,jj,mz)                     ! Refer. Level Temper.
          zttSBL=    grvinv * gplvDY(ii,jj,mz)-    sh(ii,jj)       ! Refer. Level Height
          LMO   =             SLlmol(ii,jj,1)                      ! Monin-Obukhov Length
                                                                   !
        IF (LMO.GT.0)                                         THEN ! STABLE   Situation
          zetv_i=                         zvvAWS/max(epsi,LMO)     ! Stability Limit
          zetv_i=               zetv_i-   zvvSBL/max(epsi,LMO)     !
          zetT_i=                         zttAWS/max(epsi,LMO)     ! zeta MAX = 4.28
          zetT_i=               zetT_i-   zttSBL/max(epsi,LMO)     !
          Psi__i=  -6.0*sign(1.,zetv_i)*min(4.28,abs(zetv_i))      !
          Psih_i=  -6.0*sign(1.,zetT_i)*min(4.28,abs(zetT_i))      !
        ELSE                                                       ! UNSTABLE Situation
          xx2Psi=sqrt(sqrt(1.-15.*zvvAWS       /min(-epsi,LMO)))   !
          xx1Psi=sqrt(sqrt(1.-15.*      zvvSBL /min(-epsi,LMO)))   !
          Psi__i=2. * log(0.5*(1.+xx2Psi))                         ! 
     .          +     log(0.5*(1.+xx2Psi*xx2Psi))                  ! 
     .          -2. *atan(        xx2Psi        )                  ! 
     .          -2. * log(0.5*(1.+xx1Psi))                         ! 
     .          -     log(0.5*(1.+xx1Psi*xx1Psi))                  ! 
     .          +2. *atan(        xx1Psi        )                  ! 
          yy2Psi=     sqrt(1.- 9.* zttAWS      /min(-epsi,LMO))    !
          yy1Psi=     sqrt(1.- 9.*      zttSBL /min(-epsi,LMO))    !
          Psih_i=     log(0.5*(1.+yy2Psi))                         ! 
     .          -     log(0.5*(1.+yy1Psi))                         !
        END IF                                                     !
          FF_AWS     = VV_SBL                                      ! Businger, 1973, 
     .               + SLuusl(ii,jj,1)*(log(zvvAWS/zvvSBL)-Psi__i) ! WkShop on mi-mto,
     .                                / 0.4                        !       (3.7) p.77
          FF_AWS     =              max(0.0,FF_AWS)                !
          TT_AWS     = TT_SBL                                      ! Businger, 1973, 
     .          + 0.74*SLutsl(ii,jj,1)*(log(zttAWS/zttSBL)-Psih_i) ! WkShop on mi-mto,
     .       /max(epsi,SLuusl(ii,jj,1))                            !       (3.8) p.77

C +     ************
        CALL UNwrite (ID__nc, 'FF_AWS ', iprAWS, 1, 1, 1      ,FF_AWS)
        CALL UNwrite (ID__nc, 'TT_AWS ', iprAWS, 1, 1, 1      ,TT_AWS)
C +     ************

C +   Wind Direction and Wind Speed 
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        deglon = GElonh(ii,jj) * 15.
        deglat = GElatr(ii,jj) / degrad

C +          ******
        call DDD_FF(uu    ,vv    ,deglon,deglat,xxkm(ii),yykm(jj)
     .             ,GEddxx,GEtrue,GElat0,GElon0,FF(1)   ,DD(1)   )
C +          ******

C +--Surface Thermodynamics
C +  ~~~~~~~~~~~~~~~~~~~~~~
          WK__3    =(pstDy (ii,jj)     + ptopDY)  * 10.
          WK__4    = TairSL(ii,jj)
          WK__5    = TsrfSL(ii,jj,1)

C +     ************
        CALL UNwrite (ID__nc, 'FF     ', iprAWS, 1, 1, 1      , FF   )
        CALL UNwrite (ID__nc, 'DD     ', iprAWS, 1, 1, 1      , DD   )
        CALL UNwrite (ID__nc, 'ppSurf ', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'TaSurf ', iprAWS, 1, 1, 1      , WK__4)
        CALL UNwrite (ID__nc, 'TTSurf ', iprAWS, 1, 1, 1      , WK__5)
C +     ************

C +--SBL Turbulence
C +  ~~~~~~~~~~~~~~
          WK__1    = SLuuSL(ii,jj,1)
          WK__2    = SLutSL(ii,jj,1)
          WK__3    = SLuqSL(ii,jj,1)
          WK__4    = SLusSL(ii,jj,1)
          WK__5    = SaltSN(ii,jj,1)

C +     ************
        CALL UNwrite (ID__nc, 'u_star ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'uT_star', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'uq_star', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'us_star', iprAWS, 1, 1, 1      , WK__4)
        CALL UNwrite (ID__nc, 'us_BSth', iprAWS, 1, 1, 1      , WK__5)
C +     ************

          WK__1    = SL_z0( ii,jj,1) * 1.e3
          WK__2    = SL_r0( ii,jj,1) * 1.e3
          WK__3    = Z0emBS(ii,jj,1) * 1.e3

C +     ************
        CALL UNwrite (ID__nc, 'z0_m   ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'z0_h   ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'z0_eff ', iprAWS, 1, 1, 1      , WK__3)
C +     ************


C +--Soil, Ice and Snow
C +  ------------------

          WK__1    = albeSL(ii,jj)
          WK__2    = nssSNo(ii,jj,1)
          WK__3    = nisSNo(ii,jj,1)
          WK__4    = issSNo(ii,jj,1)

C +     ************
        CALL UNwrite (ID__nc, 'Albedo ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'n_SNOW ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'n__ICE ', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'nSUPER ', iprAWS, 1, 1, 1      , WK__4)
C +     ************

          DO k=1,nsnsno
          WK_n1(1,1,k) = 0.
          WK_n2(1,1,k) = 0.
          WK_n3(1,1,k) = 0.
          WK_n4(1,1,k) = 0.
          WK_n5(1,1,k) = 0.
          WK_n6(1,1,k) = 0.
          WK_n7(1,1,k) = 0.
          WK_n8(1,1,k) = 0.
          END DO

          WK__4    = 0.
        IF          (nssSNo(ii,jj,1).GT.0)              THEN
          ns       = nssSNo(ii,jj,1)
          DO k=1,ns
          WK_n1(1,1,k) = tisSNo(ii,jj,1,    ns+1-k)
          WK_n2(1,1,k) = dzsSNo(ii,jj,1,    ns+1-k)
          WK_n3(1,1,k) = rosSNo(ii,jj,1,    ns+1-k)
          WK_n4(1,1,k) = wasSNo(ii,jj,1,    ns+1-k)
          WK_n5(1,1,k) = g1sSNo(ii,jj,1,    ns+1-k)
          WK_n6(1,1,k) = g2sSNo(ii,jj,1,    ns+1-k)
          WK_n7(1,1,k) = agsSNo(ii,jj,1,    ns+1-k)
          WK_n8(1,1,k) = nhsSNo(ii,jj,1,    ns+1-k)
          WK__4        = WK__4 
     .                 + rosSNo(ii,jj,1,         k)
     .                  *wasSNo(ii,jj,1,         k)
     .                 + dzsSNo(ii,jj,1,         k)
          END DO
        ELSE
          ns = 0
        END IF

          DO k=ns+1,ns+llx
          WK_n1(1,1,k) = TsolTV(ii,jj,1,llx+ns+1-k)
          WK_n2(1,1,k) = dz_dSV(            ns+1-k)
          WK_n3(1,1,k) = ro_Ice
          WK_n4(1,1,k) = eta_TV(ii,jj,1,llx+ns+1-k)
          END DO


C +     ************
        CALL UNwrite (ID__nc, 'T__SIS ', iprAWS, 1, 1, nsnsno , WK_n1)
        CALL UNwrite (ID__nc, 'dz_SIS ', iprAWS, 1, 1, nsnsno , WK_n2)
        CALL UNwrite (ID__nc, 'ro_SIS ', iprAWS, 1, 1, nsnsno , WK_n3)
        CALL UNwrite (ID__nc, 'wa_SIS ', iprAWS, 1, 1, nsnsno , WK_n4)
        CALL UNwrite (ID__nc, 'G1SNOW ', iprAWS, 1, 1, nsnsno , WK_n5)
        CALL UNwrite (ID__nc, 'G2SNOW ', iprAWS, 1, 1, nsnsno , WK_n6)
        CALL UNwrite (ID__nc, 'AgSNOW ', iprAWS, 1, 1, nsnsno , WK_n7)
        CALL UNwrite (ID__nc, 'HiSNOW ', iprAWS, 1, 1, nsnsno , WK_n8)
C +     ************


C +--Water Budget
C +  ------------

          WK__1    = runoTV(ii,jj)
          WK__2    = evapTV(ii,jj)
          WK__3    = SWaSNo(ii,jj,1)

C +     ************
        CALL UNwrite (ID__nc, 'IntROFF', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'IntEVAP', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'SurfWAT', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'VIntWAT', iprAWS, 1, 1, 1      , WK__4)
C +     ************


C +--That 's all, folks: NetCDF File Closure
C +  =======================================

C +     ************
        CALL UNclose (ID__nc)
C +     ************

      ENDDO


C +--Work Arrays Reset
C +  =================

      do j=1,my
      do i=1,mx
        rainb(i,j) = rainHY(i,j)
        snowb(i,j) = snowHY(i,j)
      enddo
      enddo


      return
      end


      subroutine PBLtop(TKE_1D,HHH_1D,h__PSL,h__SSL)

C +---------------------------------------------------------------------------+
C |                                                               07-APR-2005 |
C |   subroutine PBLtop computes the height of the Primary   Seeing Layer PSL |
C |                                                Secondary Seeing Layer SSL |
C |                                                                           |
C |   INPUT:   TKE_1D: Turbulent Kinetic Energy                       [m2/s2] |
C |            HHH_1D: Height above the Surface                           [m] |
C |                                                                           |
C |   OUTPUT:  h__PSL: Height of the Primary   Seeing Layer               [m] |
C |            h__SSL: Height of the Secondary Seeing Layer               [m] |
C |                                                                           |
C +---------------------------------------------------------------------------+


      IMPLICIT NONE

      include "MARdim.inc"

      real     TKE_1D(mz)
      real     HHH_1D(mz)
      real     h__PSL
      real     h__SSL

      real     TKEmin
      real     TKEtop
      
      integer  k     ,kmx   ,kzi

      logical  RESET
      logical  INTERP


      DATA     TKEmin/1.e-6/


C +--Height of the Primary   Seeing Layer (PSL)
C +  ==========================================

C +--Search the lowest TKE maximum
C +  -----------------------------

            k           =              mz
            TKEtop      =  0.01*TKE_1D(k)
 1001 CONTINUE
            k           =              k-1
      IF   (k          .LE.     mzabso     )                  GO TO 1000
      IF   (TKE_1D(k)  .LT.     TKE_1D(k+1).AND.
     .      TKE_1D(k+1).GT.     TKEmin*3.00)                  GO TO 1000
C +                                    3.00     = 1/2 order of magnitude
C +        (in order to only detect a significant maximum)
                                                              GO TO 1001
 1000 CONTINUE
            kmx         =              k+1
            TKEtop      =  0.01*TKE_1D(kmx)
            TKEtop      =   max(TKEmin*1.50,TKEtop)
C +                                    1.50     = 1/4 order of magnitude


C +--Search (from above) the lowest TKE minimum above the lowest TKE maximum
C +  ------ (This mimimum may be                                ) ----------
C +         (either a  TRUE         minimum  => INTERP = .FALSE.)
C +         (    or an arbitrary small value => INTERP = .TRUE. )
C +         -----------------------------------------------------

C +--Index  of the layer containing the minimum
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            kzi         =       mzabso
      DO k= mzabso,kmx
        IF (TKE_1D(k) .LT.      TKEtop          .OR.
     .      TKE_1D(k) .LT.      TKE_1D(k-1)*0.3     )               THEN
            kzi         =              k
         IF(TKE_1D(k) .LT.      TKEtop              )               THEN
            INTERP      =      .TRUE.
         ELSE
            INTERP      =      .FALSE.
         END IF
        END IF
      ENDDO

C +--Height of the                      minimum
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            k           =       kzi     
      IF   (kzi       .LE.      mzabso+1)                           THEN
            h__PSL      =       HHH_1D(mz)
      ELSE
        IF (INTERP)                                                 THEN
            h__PSL      =       HHH_1D(k+1)
     .                        +(HHH_1D(k)  -HHH_1D(k+1))
     .                        *(TKEtop     -TKE_1D(k+1))
     .                        /(TKE_1D(k)  -TKE_1D(k+1))
        ELSE
            h__PSL      =       HHH_1D(k)
        END IF
      END IF

            h__PSL      =   min(h__PSL     ,HHH_1D(1)  )
            h__PSL      =   max(HHH_1D(mz) ,h__PSL     )


C +--Height of the Secondary Seeing Layer (SSL)
C +  ==========================================

            RESET       =      .TRUE.


C +--Search the        TKE minimum above the Primary Seeing Layer (PSL)
C +  (necessary if the TKE has decreased below the minimum value)
C +  ------------------------------------------------------------------

      IF   (INTERP)                                                 THEN
            k           =       kzi + 1 
 1011   CONTINUE
            k           =       k-1
        IF (k         .LE.      mzabso     )                  GO TO 1010
        IF (TKE_1D(k) .LT.      TKE_1D(k+1))                  GO TO 1011
 1010   CONTINUE
      ELSE
            k           =       kzi
      END IF


C +--Search the first  TKE maximum above the Primary Seeing Layer (PSL)
C +  ------------------------------------------------------------------

            kmx         =       kzi
            k           =       k+1
 1021   CONTINUE
            k           =       k-1
      IF   (k            .LE.   mzabso)                       GO TO 1020
        IF (TKE_1D(k)    .GT.   TKE_1D(k-1)     .AND.
     .      TKE_1D(k)    .GT.   TKE_1D(k+1)     .AND.
     .      TKE_1D(k)    .GT.   TKEmin*3.0           )              THEN
C +                                    3.0      = 1/2 order of magnitude
C +        (in order to only detect a significant maximum)


C +--Define the TKE at the SSL top from the largest maximum in the SSL
C +  (thus examine the remaining upper part of the atmospheric column)
C +  -----------------------------------------------------------------

         IF(RESET)                                                  THEN
            RESET       =      .FALSE. ! indicates TKEtop is initialized
            TKEtop      =       0.00   !
         END IF
         IF(TKEtop       .LT.   TKE_1D(k)  *0.01)                   THEN
            TKEtop      =       TKE_1D(k)  *0.01
            kmx         =              k
         END IF
        END IF
                                                              GO TO 1021
 1020 CONTINUE
            TKEtop      =   max(TKEmin*3.0 ,TKEtop)
C +                                    3.0      = 1/2 order of magnitude


C +--Search (from above) the SSL top            above the SSL    TKE maximum
C +  ------ (This         may be                                ) ----------
C +         (either a  TRUE         minimum  => INTERP = .FALSE.)
C +         (    or an arbitrary small value => INTERP = .TRUE. )
C +         -----------------------------------------------------

            kzi         =       mzabso
      DO k= kmx,mzabso,-1     
        IF (TKE_1D(k) .GT.      TKEtop)
     .      kzi         =       k
      ENDDO

            k           =       kzi   -1
      IF   (kzi       .LE.      mzabso+1)                           THEN
            h__SSL      =       HHH_1D(mz)
      ELSE
        IF (INTERP)                                                 THEN
            h__SSL      =       HHH_1D(k+1)
     .                        +(HHH_1D(k)  -HHH_1D(k+1))
     .                        *(TKEtop     -TKE_1D(k+1))
     .                        /(TKE_1D(k)  -TKE_1D(k+1))
        ELSE
            h__SSL      =       HHH_1D(k)
        END IF
      END IF

            h__SSL      =   min(h__SSL     ,HHH_1D(1)  )
            h__SSL      =   max(h__PSL     ,h__SSL     )

      RETURN
      END


      subroutine DDD_FF(u_wind,v_wind,lon   ,lat   ,x  ,y
     .                 ,Dir__X,trulat,lat__0,lon__0,FFF,DDD)

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                             26-10-2004  MAR |
C |   subroutine DDD_FF computes the wind speed and direction              |
C |                              from stereographic MAR wind components    |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


      REAL     u_wind,v_wind
      REAL     lon   ,lat
      REAL     x     ,y   
      REAL     Dir__X,trulat
      REAL     lat__0,lon__0
      REAL     cosbet,sinbet
      REAL     up    ,vp
      REAL     sinla0,cosla0
      REAL     coslon,sinlon
      REAL     coslat,sinlat
      REAL     denomi
      REAL     cosalp,sinalp
      REAL     uu    ,vv    ,magnif
      REAL     FFF   ,DDD


C +--Restauration Coordonn?es ds un plan o? (x,y)=(E,N) au Centre de Projection

      cosbet=cos((Dir__X-90.)*3.1415/180.)
      sinbet=sin((Dir__X-90.)*3.1415/180.)
      up    =     u_wind *cosbet + v_wind *sinbet
      vp    =    -u_wind *sinbet + v_wind *cosbet

C +--Coordonn?es G?ographiques de chaque Point

      sinla0=sin( lat__0     *3.1415/180.)
      cosla0=cos( lat__0     *3.1415/180.)

      coslon=cos((lon-lon__0)*3.1415/180.) ! lat,lon: coord.   g?ograph. d'un pt
      sinlon=sin((lon-lon__0)*3.1415/180.) !   x,  y: coord.st?r?ograph. d'un pt

      sinlat=(sinla0*(4*6371*6371-x*x-y*y)+4*6371*y*cosla0)
     .      /        (4*6371*6371+x*x+y*y)
      coslat=sqrt(1-sinlat*sinlat)

C +--Rotation locale par rapport au syst?me d'axes (E,N)

      denomi= 1.+sinla0*sinlat+cosla0*coslat *coslon
      cosalp=(cosla0*coslat+(1+sinla0*sinlat)*coslon)/denomi
      sinalp=(0            -(  sinlat+sinla0)*sinlon)/denomi

      uu=up*cosalp-vp*sinalp
      vv=up*sinalp+vp*cosalp

C +--Direction du vecteur Vent

          DDD=270.-(180./3.1415)*atan2(vv,uu)
      IF (DDD     .LT.        0.)                                   THEN
          DDD    = DDD    + 360.
      ENDIF
      IF (DDD     .GT.      360.)                                   THEN
          DDD    = DDD    - 360.
      ENDIF

C +--Vitesse   du Vecteur Vent

      magnif=(1+cos(trulat*3.14159/180.))/(1-sin(lat *3.1415/180.))
      FFF   =  sqrt(u_wind*u_wind+v_wind*v_wind)

      RETURN
      END
      subroutine uv2fd(ruu,rvv,rFF,rDD,iu,ju)

C +------------------------------------------------------------------------+
C |                                                                        |
C | MAR OUTPUT                                             26-10-2004  MAR |
C |   SubRoutine uv2fd  transforms Wind(u,v) into Wind(ff,dd)              |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'


C +--Local   Variables
C +  =================

      integer            iu   ,ju
      REAL               uugeo,ruu,rFF
      REAL               vvgeo,rvv,rDD

      REAL               conv
      parameter         (conv  = 15.0*3.141592/180.0)          ! Conversion
                                                               ! Hour -> Radian

C +--From (u,v) to (ff,dd)
C +  =====================

        IF   (ruu   .NE.0.0 .OR.  rvv   .NE.0.0)                    THEN

              uugeo    = (GElonh(iu+1,ju) - GElonh(iu,ju))*conv/dx
     .                   *earthr       *cos(GElatr(iu,ju))*    ruu
     .                 + (GElonh(iu,ju+1) - GElonh(iu,ju))*conv/dx
     .                   *earthr       *cos(GElatr(iu,ju))*    rvv
     
              vvgeo    = (GElatr(iu+1,ju) - GElatr(iu,ju))     /dx
     .                   *earthr                          *    ruu
     .                 + (GElatr(iu,ju+1) - GElatr(iu,ju))     /dx
     .                   *earthr                          *    rvv

              rDD       =  0.0

          IF (uugeo    .GT. 0.0 .and. vvgeo    .GE.0.0)  
     .        rDD       =    3.0*pi/2.0 - atan(vvgeo    /uugeo    )
          IF (uugeo    .LT. 0.0 .and. vvgeo    .GE.0.0)  
     .        rDD       =    5.0*Pi/2.0  -atan(vvgeo    /uugeo    )
          IF (uugeo    .LT. 0.0 .and. vvgeo    .LE.0.0)  
     .        rDD       =    5.0*Pi/2.0  -atan(vvgeo    /uugeo    )
          IF (uugeo    .GT. 0.0 .and. vvgeo    .LE.0.0)  
     .        rDD       =    3.0*Pi/2.0  -atan(vvgeo    /uugeo    )
          IF (uugeo    .EQ. 0.0 .and. vvgeo    .GE.0.0) 
     .        rDD       =        Pi
          IF (uugeo    .EQ. 0.0 .and. vvgeo    .LE.0.0) 
     .        rDD       =        0.0
          IF (rDD       .GT. 2.0*Pi)
     .        rDD       = rDD-2.0*Pi

              rDD       = rDD        / degrad

          if( rDD       .LT. 0.0 )       
     .        rDD       = rDD        + 360.0
     
              rDD       = max(0.0,min(360.0,rDD))
     
              rFF      = sqrt(ruu*ruu + rvv*rvv)
       
        ELSE
              rFF     = 0.0 
              rDD     = 0.0
        END IF 

      return
      end


      block data AWS_nc_DATA

C +----------------------------------------------------------------------------+
C |                                                                            |
C | MAR OUTPUT   Generic Routine                          Mon  6-Feb-2012  MAR |
C |   Manned and Automatic Weather Stations (AWS) Geographic Coordinates       |
C |                                                                            |
C +----------------------------------------------------------------------------+


C +--General Variables
C +  =================

      integer                  n_AWS,  n
      parameter               (n_AWS=333)
      integer            AWSio(n_AWS),AWS_i(n_AWS),AWS_j(n_AWS)
      integer            nnAWS
      REAL               AWSla(n_AWS),AWSlo(n_AWS),AWS_z(n_AWS)
      character*6        AWS_0(n_AWS)
      character*8                     AWS_1(n_AWS),AWS_2(n_AWS)

      common /AWS_nc_INT/AWSio       ,AWS_i       ,AWS_j
     .                  ,nnAWS
      common /AWS_nc_REA/AWSla       ,AWSlo       ,AWS_z
      common /AWS_nc_CH6/AWS_0
      common /AWS_nc_CH8/             AWS_1       ,AWS_2


C +--DATA
C +  ====

C +--ANT
C +  ---

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=001,094)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'CapSpe' ,'Cape Spe','ncer    ', -77.97, 167.55,   30., 1,  !   1
     .     'HerbA1' ,'Herbie A','lley    ', -78.10, 166.67,   30., 1,  !   2
     .     'Harry_' ,'Harry   ','        ', -83.00,-121.39,  945., 1,  !   3
     .     'CapBir' ,'Cape Bir','d       ', -77.22, 166.44,  100., 1,  !   4
     .     'Butler' ,'Butler I','sland   ', -72.21, 299.84,   91., 1,  !   5
     .     'Byrd__' ,'Byrd    ','        ', -80.01,-119.40, 1530., 1,  !   6
     .     'Dome-F' ,'Dome F  ','        ', -77.31,  39.70, 3810., 1,  !   7
     .     'Manuel' ,'Manuela ','        ', -74.95, 163.69,   78., 1,  !   8
     .     'Marble' ,'Marble P','oint    ', -77.44, 163.75,  120., 1,  !   9
     .     'Whitlo' ,'Whitlock','        ', -76.14, 168.39,  275., 1,  !  10
     .     'Lettau' ,'Lettau  ','        ', -82.52,-174.45,   30., 1,  !  11
     .     'PortMa' ,'Port Mar','tin     ', -66.82, 141.40,   39., 1,  !  12
     .     'PengPt' ,'Penguin ','Point   ', -67.62, 146.18,   30., 1,  !  13
     .     'Gill_1' ,'Gill    ','        ', -79.99,-178.61,   25., 1,  !  14
     .     'Schwer' ,'Schwerdt','feger   ', -79.90, 169.97,   60., 1,  !  15
     .     'D10___' ,'D-10    ','        ', -66.71, 139.83,  243., 1,  !  16
     .     'Elaine' ,'Elaine  ','        ', -83.13, 174.17,   60., 1,  !  17
     .     'Ski_Hi' ,'Ski Hi  ','        ', -74.79, 289.51, 1395., 1,  !  18
     .     'Relay_' ,'Relay St','        ', -74.02,  43.06, 3353., 1,  !  19
     .     'Linda_' ,'Linda   ','        ', -78.46, 168.38,   50., 1,  !  20
     .     'Uranus' ,'Uranus G','lacier  ', -71.43, 291.07,  780., 1,  !  21
     .     'MADISO' ,'MADISON ','        ',  43.08, -89.38,    0., 1,  !  22
     .     'Doug__' ,'Doug    ','        ', -82.32,-113.24, 1430., 1,  !  23
     .     'BonaPt' ,'Bonapart','e Point ', -64.78, 295.93,    8., 1,  !  24
     .     'Nico__' ,'Nico    ','        ', -89.00,  89.67, 2935., 1,  !  25
     .     'Limbrt' ,'Limbert ','        ', -75.42, 300.15,   40., 1,  !  26
     .     'Larsen' ,'Larsen I','ce      ', -66.95, 299.10,   17., 1,  !  27
     .     'Wndlss' ,'Wndlss B','t       ', -77.73, 167.70,   60., 1,  !  28
     .     'Ferrel' ,'Ferrell ','        ', -77.91, 170.82,   45., 1,  !  29
     .     'Kirkwd' ,'Kirkwood','        ', -68.34, 290.99,   30., 1,  !  30
     .     'Dismal' ,'Dismal I','s       ', -68.09, 291.18,   10., 1,  !  31
     .     'Marily' ,'Marilyn ','        ', -79.95, 165.13,   75., 1,  !  32
     .     'MinnaB' ,'Minna Bl','uff     ', -78.55, 166.66,  920., 1,  !  33
     .     'PegasS' ,'Pegasus ','South   ', -77.99, 166.58,   10., 1,  !  34
     .     'SipleD' ,'Siple Do','me      ', -81.66,-148.77,  620., 1,  !  35
     .     'Sutton' ,'Sutton  ','        ', -67.10, 141.40,  871., 1,  !  36
     .     'RacerR' ,'Racer Ro','ck      ', -64.07, 298.39,   17., 1,  !  37
     .     'YoungI' ,'Young Is','        ', -66.20, 162.30,   30., 1,  !  38
     .     'MtSipl' ,'Mount Si','ple     ', -73.20,-127.05,  230., 1,  !  39
     .     'PossIs' ,'Poss Is ','        ', -71.89, 171.21,   30., 1,  !  40
     .     'Henry_' ,'Henry   ','        ', -89.01,  -1.02, 2755., 1,  !  41
     .     'D47___' ,'D-47    ','        ', -67.40, 138.73, 1560., 2,  !  42
     .     'D57___' ,'D-57    ','        ', -68.30, 137.87, 2105., 2,  !  43
     .     'CapDen' ,'Cape Den','ison    ', -67.01, 142.66,   31., 1,  !  44
     .     'DomeC2' ,'Dome C2 ','        ', -75.12, 123.37, 3250., 1,  !  45
     .     'SwiBnk' ,'Swithinb','ank     ', -81.20,-126.17,  945., 1,  !  46
     .     'PegasN' ,'Pegasus ','North   ', -77.95, 166.50,    8., 1,  !  47
     .     'Theres' ,'Theresa ','        ', -84.60,-115.81, 1463., 1,  !  48
     .     'Mizuho' ,'Mizuho  ','        ', -70.70,  44.29, 2260., 1,  !  49
     .     'LaurII' ,'Laurie I','I       ', -77.55, 170.81,   30., 1,  !  50
     .     'Elizab' ,'Elizabet','h       ', -82.61,-137.08,  519., 1,  !  51
     .     'Briana' ,'Brianna ','        ', -83.89,-134.15,  526., 1,  !  52
     .     'Erin__' ,'Erin    ','        ', -84.90, 231.17,  990., 1,  !  53
     .     'WillF1' ,'Willie F','ield    ', -77.87, 167.02,   20., 1,  !  54
     .     'Young2' ,'Young Is','land    ', -62.23, 162.28,   30., 1,  !  55
     .     'CleanA' ,'Clean Ai','r       ', -90.00,   0.00, 2835., 1,  !  56
     .     'OdellG' ,'Odell Gl','acier   ', -76.63, 160.05, 1637., 1,  !  57
     .     'HerbA2' ,'Herbie A','lley    ', -77.97, 167.54,   24., 1,  !  58
     .     'SkyBlu' ,'Sky Blu ','        ', -74.79, 288.51, 1395., 1,  !  59
     .     'A028-A' ,'A028-A  ','        ', -67.59, 112.22, 1622., 1,  !  60
     .     'A028__' ,'A028    ','        ', -67.59, 112.22, 1622., 1,  !  61
     .     'GC41__' ,'GC41    ','        ', -70.40, 111.26, 2791., 1,  !  62
     .     'GC46__' ,'GC46    ','        ', -70.86, 109.84, 3096., 1,  !  63
     .     'GF08__' ,'GF08    ','        ', -67.51, 102.18, 2123., 1,  !  64
     .     'LawDom' ,'Law Dome','        ', -65.27, 112.74, 1376., 1,  !  65
     .     'DDU___' ,'DDU     ','        ', -66.67, 140.02,   42., 2,  !  66
     .     'Mawson' ,'Mawson  ','        ', -67.60,  62.87,   10., 2,  !  67
     .     'Casey_' ,'Casey   ','        ', -66.28, 110.52,   40., 2,  !  68
     .     'McMurd' ,'McMurdo ','(Fogle) ', -77.82, 166.75,  202., 2,  !  69
     .     'Mirny_' ,'Mirny   ','        ', -66.33,  93.01,   30., 2,  !  70
     .     'Vostok' ,'Vostok  ','        ', -78.45, 106.84, 3471., 2,  !  71
     .     'Alison' ,'Allison ','        ', -89.88, 300.00, 2835., 0,  !  72
     .     'Bowers' ,'Bowers  ','        ', -85.20, 163.40, 2090., 0,  !  73
     .     'D80___' ,'D80     ','        ', -70.02, 134.72, 2500., 0,  !  74
     .     'Dollem' ,'Dolleman',' Island ', -70.58, 299.08,  396., 0,  !  75
c #C1.     'DomeC1' ,'Dome  C ','        ', -74.65, 124.40, 3232., 0,  !  76
     .     'Dome-A' ,'Dome  A ','China   ', -81.00,  77.00, 4100., 0,  !  76
     .     'DomeCA' ,'Dome  C ','AMRC    ', -74.50, 123.00, 3280., 0,  !  77
     .     'DomeCE' ,'Dome  C ','EPICA   ', -75.11, 123.32, 3232., 2,  !  78
     .     'Eneid_' ,'Eneid   ','(TNB)   ', -74.41, 164.06,   88., 0,  !  79
     .     'Gill_2' ,'Gill    ','        ', -80.00, 181.00,   55., 0,  !  80
     .     'Maning' ,'Manning ','        ', -78.80, 166.80,   65., 0,  !  81
     .     'Martha' ,'Martha  ','        ', -78.31,-172.50,   42., 0,  !  82
     .     'Patrik' ,'Patrick ','        ', -89.88,  45.00, 2835., 0,  !  83
     .     'RidgeB' ,'Ridge B ','        ', -77.08,  94.92, 3400., 1,  !  84
     .     'Tiffan' ,'Tiffany ','        ', -77.95, 168.17,   25., 0,  !  85
     .     'Whitlo' ,'Whitlok ','        ', -76.24, 168.66,  274., 0,  !  86
     .     'WindlB' ,'Windless',' Bight  ', -77.70, 167.70,   40., 0,  !  87
     .     'GPS2__' ,'GPS2    ','        ', -74.61, 157.38, 1804., 0,  !  88
     .     '31Dpt_' ,'31Dpt   ','        ', -74.06, 155.93, 2041., 0,  !  89
     .     'M2____' ,'M2      ','        ', -74.80, 151.10, 2272., 0,  !  90
     .     'MidPt2' ,'MdPt2   ','        ', -75.53, 145.92, 2454., 0,  !  91
     .     'D2____' ,'D2      ','        ', -75.65, 140.48, 2482., 0,  !  92
     .     'D4____' ,'D4      ','        ', -75.60, 135.83, 2793., 0,  !  93
     .     'D6____' ,'D6      ','        ', -75.44, 129.63, 3038., 0/  !  94


      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n= 95,100)
     .  /  'Rother' ,'Rothera ','        ', -67.50, 291.90,   16., 1,  !  95
     .     'Primav' ,'Primaver','a       ', -64.20, 259.00,   50., 1,  !  96
     .     'O_Higg' ,'O Higgin','s       ', -63.30, 302.10,   10., 1,  !  97
     .     'Bellin' ,'Bellings','hausen  ', -62.20, 301.10,   16., 1,  !  98
     .     'Adelai' ,'Adelaide','        ', -67.80, 302.10,   26., 1,  !  99
     .     'SanMar' ,'San Mart','in      ', -68.10, 292.90,    4., 1/  ! 100


C +--AFW
C +  ---

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=101,175)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'Bamako' ,' Bamako ','    Mali',  12.32,  -7.57,  381., 2,  ! 101
     .     'Tombou' ,' Tombouc','tou Mali',  16.43,  -3.00,  264., 2,  ! 102
     .     'NiameA' ,' Niamey-','Aero  NI',  13.48,   2.16,  223., 2,  ! 103
     .     'Abidjn' ,' Abidjan',' Cote Iv',   5.15,  -3.56,    8., 2,  ! 104
     .     'Tamanr' ,' Tamanra','sset  AL',  22.78,   5.51, 1377., 2,  ! 105
     .     'Dakar_' ,' Dakar  ',' Senegal',  14.44, -17.30,   27., 2,  ! 106
     .     'Adiake' ,' Adiake ',' Cote Iv',    5.3,   -3.3,    7., 0,  ! 107
     .     'Bondok' ,' Bondouk',' Cote Iv',  08.03,  -2.47,  370., 0,  ! 108
     .     'Bouake' ,' Bouake ',' Cote Iv',  07.44,  -5.04,  376., 0,  ! 109
     .     'MAN___' ,' Man    ',' Cote Iv',  07.23,  -7.31,  340., 0,  ! 110
     .     'Korhog' ,' Korhogo',' Cote Iv',  09.25,  -5.37,  381., 0,  ! 111
     .     'Sassan' ,' Sassand',' Cote Iv',  04.57,  -6.05,   66., 0,  ! 112
     .     'Tabou_' ,' Tabou  ',' Cote Iv',  04.25,  -7.22,   21., 0,  ! 113
     .     'Dimbok' ,' Dimbokr',' Cote Iv',  06.39,  -4.42,   92., 0,  ! 114
     .     'Odienn' ,' Odienn?',' Cote Iv',  09.30,  -7.34,  421., 0,  ! 115
     .     'Bobodi' ,' Bobodio',' Burkfas',  11.10,  -4.19,  460., 0,  ! 116
     .     'Ouaga_' ,' Ouaga  ',' Burkfas',  12.21,  -1.31,  306., 0,  ! 117
     .     'Ouahig' ,' Ouahigo',' Burkfas',  13.34,  -2.25,  336., 0,  ! 118
     .     'Fadago' ,' Fadagou',' Burkfas',  12.02,   0.22,  309., 0,  ! 119
     .     'Conakr' ,' Conakry',' Guinee ',  09.34, -13.37,   26., 2,  ! 120
     .     'Labe__' ,' Labe   ',' Guinee ',  11.19, -12.18, 1026., 0,  ! 121
     .     'Nzere_' ,' Nzere  ',' Guinee ',   7.44,  -8.50,  470., 0,  ! 122
     .     'Siguir' ,' Siguiri',' Guinee ',  11.26,  -9.10,  366., 0,  ! 123
     .     'Bissau' ,' Bissau ',' Gbissau',  11.53, -15.39,   26., 2,  ! 124
     .     'StLoui' ,' Stlouis',' Senegal',  16.03, -16.26,   04., 0,  ! 125
     .     'Matam_' ,' Matam  ',' Senegal',  15.39, -13.15,   17., 0,  ! 126
     .     'Tambac' ,' Tambaco',' Senegal',  13.46, -13.41,   50., 0,  ! 127
     .     'Kolda_' ,' Kolda  ',' Senegal',  12.53, -14.58,   10., 0,  ! 128
     .     'Ziguin' ,' Ziguinc',' Senegal',  12.33, -16.16,   23., 0,  ! 129
     .     'Diourb' ,' Diourbe',' Senegal',  14.39, -16.14,    9., 0,  ! 130
     .     'Kedouk' ,' Kedouko',' Senegal',  12.34, -12.13,  167., 0,  ! 131
     .     'Lungi_' ,' Lungi  ',' Sierral',   8.37, -13.12,   27., 0,  ! 132
     .     'Robert' ,' Robertf',' Liberia',   6.15, -10.21,   18., 0,  ! 133
     .     'Nouakc' ,' Nouakch',' Maurita',  18.06, -15.57,    3., 0,  ! 134
     .     'Nouadi' ,' Nouadib',' Maurita',  20.56, -17.02,    3., 0,  ! 135
     .     'Zouera' ,' Zouerat',' Maurita',  22.45, -12.29,  343., 0,  ! 136
     .     'Nema__' ,' Nema   ',' Maurita',  16.36,  -7.16,  269., 0,  ! 137
     .     'Atar__' ,' Atar   ',' Maurita',  20.31, -13.04,  224., 0,  ! 138
     .     'Kayes_' ,' Kayes  ',' Mali   ',  14.26, -11.26,   47., 0,  ! 139
     .     'Kidal_' ,' Kidal  ',' Mali   ',  18.26,   1.21,  459., 0,  ! 140
     .     'Mopti_' ,' Mopti  ',' Mali   ',  14.31,  -4.06,  272., 0,  ! 141
     .     'Sikaso' ,' Sikasso',' Mali   ',  11.21,  -5.41,  375., 0,  ! 142
     .     'Accra_' ,' Accra  ',' Ghana  ',   5.36,  -0.10,   69., 0,  ! 143
     .     'Kumasi' ,' Kumassi',' Ghana  ',   6.43,  -1.36,  293., 2,  ! 144
     .     'Tamal_' ,' Tamal  ',' Ghana  ',   9.30,  -0.51,  173., 0,  ! 145
     .     'WAG___' ,' WA     ',' Ghana  ',  10.03,  -0.30,  323., 0,  ! 146
     .     'Wenchi' ,' Wenchi ',' Ghana  ',   7.43,  -2.06,  340., 0,  ! 147
     .     'Ada___' ,' Ada    ',' Ghana  ',   5.47,  -0.38,    7., 0,  ! 148
     .     'Abuja_' ,' Abuja  ',' Nigeria',   9.15,   7.00,  344., 0,  ! 149
     .     'Lagos_' ,' Lagos  ',' Nigeria',   6.35,   3.20,   38., 2,  ! 150
     .     'Maidug' ,' Maidugu',' Nigeria',  11.51,  13.05,  354., 0,  ! 151
     .     'Cotonu' ,' Cotonou',' Benin  ',   6.36,   2.44,    9., 2,  ! 152
     .     'Paraku' ,' Parakou',' Benin  ',   9.21,   2.37,  393., 2,  ! 153
     .     'Kandi_' ,' Kandi  ',' Benin  ',  11.08,   2.56,  292., 0,  ! 154
     .     'Natiti' ,' Natitin',' Benin  ',  10.19,   1.23,  461., 0,  ! 155
     .     'Lome__' ,' Lome   ','  Togo  ',   6.10,   1.15,   25., 0,  ! 156
     .     'Atakpa' ,' Atakpam','  Togo  ',   7.35,   1.07,  402., 0,  ! 157
     .     'Sokode' ,' Sokode ','  Togo  ',   8.59,   1.09,  387., 0,  ! 158
     .     'Dapaon' ,' Dapaon ','  Togo  ',  10.52,   0.15,  330., 0,  ! 159
     .     'Faya__' ,' Faya   ','  Tchad ',  18.00,  19.10,  234., 0,  ! 160
     .     'Moundu' ,' Moundou','  Tchad ',   8.37,  16.10,  422., 0,  ! 161
     .     'Ndjame' ,' Ndjamen','a Tchad ',  12.08,  15.02,  295., 2,  ! 162
     .     'Sarh__' ,' Sarh   ','  Tchad ',   9.09,  18.23,  365., 0,  ! 163
     .     'Douala' ,' Douala ',' Camerou',   4.00,   9.44,    9., 2,  ! 164
     .     'Garoua' ,' Garoua ',' Camerou',   9.20,  13.23,  244., 0,  ! 165
     .     'Malabo' ,' Malabo ',' Camerou',   3.49,   8.46,   56., 0,  ! 166
     .     'Tindou' ,' Tindouf',' Algerie',  27.40,  -8.08,  431., 0,  ! 167
     .     'Agadir' ,' Agadir ',' Maroc  ',  30.20,  -9.24,   74., 0,  ! 168
     .     'Villac' ,'Villacis',' Sah Occ',  23.42, -15.52,   10., 0,  ! 169
     .     'Daloa_' ,'Daloa   ',' Cote Iv',   6.52,  -6.28,  276., 2,  ! 170 
     .     'Gagnoa' ,'Gagnoa  ',' Cote Iv',   6.08,  -5.57,  205., 2,  ! 171
     .     'SanPed' ,'San Pedr',' Cote Iv',   4.45,  -6.39,   31., 2,  ! 172
     .     'Tabou_' ,'Tabou   ',' Cote Iv',   4.25,  -7.22,   20., 2,  ! 173
     .     'Yamous' ,'Yamousso',' Cote Iv',   6.54,  -5.21,  196., 2,  ! 174
     .     'Sal___' ,'Cap Vert',' Cote Iv',  16.44, -22.57,   54., 2/  ! 175


C +--ANT (again)
C +  ---

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=176,188)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'AGO-A8' ,'AGO-A8  ','        ', -84.36, -23.86, 2103., 1,  ! 176
     .     'Davis_' ,'Davis   ','        ', -68.35,  77.59,   18., 2,  ! 177
     .     'Kenton' ,'Kenton  ','        ', -72.28, -38.82, 3185., 1,  ! 178
     .     'SantCI' ,'Santa Cl','aus Isla', -64.96, -65.67,   25., 1,  ! 179
     .     'ScottI' ,'Scott Is','land    ', -67.37,-179.97,   30., 1,  ! 180
     .     'ScottB' ,'Scott Ba','se      ', -77.51, 166.45,   94., 2,  ! 181
     .     'Ski-Hi' ,'Ski-Hi  ','        ', -74.98, -70.77, 1395., 1,  ! 182
     .     'TNB___' ,'Terra No','va Bay  ', -74.42, 164.06,   80., 1,  ! 183
     .     'Theres' ,'Theresa ','        ', -84.60,-115.81, 1463., 1,  ! 184
     .     'WhiteO' ,'White Ou','t       ', -77.87, 168.16,   30., 1,  ! 185
     .     'WhiteI' ,'White Is','land    ', -78.09, 168.01,   30., 1,  ! 186
     .     'Kohnen' ,'Kohnen  ','        ', -75.00,   0.07, 2892., 1,  ! 187
     .     'ElizBe' ,'Princess',' Elizab ', -71.90,  23.33, 1390., 2/  ! 188


C +--ANT (Balises Glacioclim Samba, Antarctique)
C +  ---

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=189,279)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     . /'B____1' ,'Balise S','AMBA   1',-66.69592,139.8985,    0., 1,  ! 189
     .  'B____2' ,'Balise S','AMBA   2',-66.69861,139.8915,    0., 1,  ! 190
     .  'B____3' ,'Balise S','AMBA   3',-66.69849,139.8801,    0., 1,  ! 191
     .  'B____4' ,'Balise S','AMBA   4',-66.69848,139.8690,    0., 1,  ! 192
     .  'B____5' ,'Balise S','AMBA   5',-66.69933,139.8545,    0., 1,  ! 193
     .  'B____6' ,'Balise S','AMBA   6',-66.70067,139.8398,    0., 1,  ! 194
     .  'B____7' ,'Balise S','AMBA   7',-66.70199,139.8236,    0., 1,  ! 195
     .  'B____8' ,'Balise S','AMBA   8',-66.70349,139.8088,    0., 1,  ! 196
     .  'B____9' ,'Balise S','AMBA   9',-66.70551,139.8021,    0., 1,  ! 197
     .  'B___10' ,'Balise S','AMBA  10',-66.70689,139.7838,    0., 1,  ! 198
     .  'B___11' ,'Balise S','AMBA  11',-66.70894,139.7746,    0., 1,  ! 199
     .  'B___12' ,'Balise S','AMBA  12',-66.71090,139.7628,    0., 1,  ! 200
     .  'B___13' ,'Balise S','AMBA  13',-66.71597,139.7400,    0., 1,  ! 201
     .  'B___14' ,'Balise S','AMBA  14',-66.72215,139.7237,    0., 1,  ! 202
     .  'B___15' ,'Balise S','AMBA  15',-66.72684,139.6999,    0., 1,  ! 203
     .  'B___16' ,'Balise S','AMBA  16',-66.73164,139.6770,    0., 1,  ! 204
     .  'B___17' ,'Balise S','AMBA  17',-66.73676,139.6557,    0., 1,  ! 205
     .  'B___18' ,'Balise S','AMBA  18',-66.74189,139.6336,    0., 1,  ! 206
     .  'B___19' ,'Balise S','AMBA  19',-66.74656,139.6087,    0., 1,  ! 207
     .  'B___20' ,'Balise S','AMBA  20',-66.75196,139.5914,    0., 1,  ! 208
     .  'B___21' ,'Balise S','AMBA  21',-66.76015,139.5783,    0., 1,  ! 209
     .  'B___22' ,'Balise S','AMBA  22',-66.76746,139.5677,    0., 1,  ! 210
     .  'B___23' ,'Balise S','AMBA  23',-66.77496,139.5568,    0., 1,  ! 211
     .  'B___24' ,'Balise S','AMBA  24',-66.78358,139.5459,    0., 1,  ! 212
     .  'B___25' ,'Balise S','AMBA  25',-66.79026,139.5321,    0., 1,  ! 213
     .  'B___26' ,'Balise S','AMBA  26',-66.79795,139.5204,    0., 1,  ! 214
     .  'B___27' ,'Balise S','AMBA  27',-66.80594,139.5093,    0., 1,  ! 215
     .  'B___28' ,'Balise S','AMBA  28',-66.81308,139.4968,    0., 1,  ! 216
     .  'B___29' ,'Balise S','AMBA  29',-66.82068,139.4857,    0., 1,  ! 217
     .  'B___30' ,'Balise S','AMBA  30',-66.82845,139.4750,    0., 1,  ! 218
     .  'B___31' ,'Balise S','AMBA  31',-66.83596,139.4623,    0., 1,  ! 219
     .  'B___32' ,'Balise S','AMBA  32',-66.84335,139.4495,    0., 1,  ! 220
     .  'B___33' ,'Balise S','AMBA  33',-66.85136,139.4383,    0., 1,  ! 221
     .  'B___34' ,'Balise S','AMBA  34',-66.85864,139.4251,    0., 1,  ! 222
     .  'B___35' ,'Balise S','AMBA  35',-66.86616,139.4135,    0., 1,  ! 223
     .  'B___36' ,'Balise S','AMBA  36',-66.87332,139.4004,    0., 1,  ! 224
     .  'B___37' ,'Balise S','AMBA  37',-66.88036,139.3867,    0., 1,  ! 225
     .  'B___38' ,'Balise S','AMBA  38',-66.88902,139.3759,    0., 1,  ! 226
     .  'B___39' ,'Balise S','AMBA  39',-66.89655,139.3638,    0., 1,  ! 227
     .  'B___40' ,'Balise S','AMBA  40',-66.90511,139.3541,    0., 1,  ! 228
     .  'B___41' ,'Balise S','AMBA  41',-66.91162,139.3414,    0., 1,  ! 229
     .  'B___42' ,'Balise S','AMBA  42',-66.92023,139.3318,    0., 1,  ! 230
     .  'B___43' ,'Balise S','AMBA  43',-66.92771,139.3174,    0., 1,  ! 231
     .  'B___44' ,'Balise S','AMBA  44',-66.93586,139.3041,    0., 1,  ! 232
     .  'B___45' ,'Balise S','AMBA  45',-66.94995,139.2799,    0., 1,  ! 233
     .  'B___46' ,'Balise S','AMBA  46',-66.96581,139.2567,    0., 1,  ! 234
     .  'B___47' ,'Balise S','AMBA  47',-66.98064,139.2319,    0., 1,  ! 235
     .  'B___48' ,'Balise S','AMBA  48',-66.99638,139.2094,    0., 1,  ! 236
     .  'B___49' ,'Balise S','AMBA  49',-67.01013,139.1807,    0., 1,  ! 237
     .  'B___50' ,'Balise S','AMBA  50',-67.02596,139.1579,    0., 1,  ! 238
     .  'B___51' ,'Balise S','AMBA  51',-67.04180,139.1366,    0., 1,  ! 239
     .  'B___52' ,'Balise S','AMBA  52',-67.05727,139.1129,    0., 1,  ! 240
     .  'B___53' ,'Balise S','AMBA  53',-67.07246,139.0880,    0., 1,  ! 241
     .  'B___54' ,'Balise S','AMBA  54',-67.08774,139.0636,    0., 1,  ! 242
     .  'B___55' ,'Balise S','AMBA  55',-67.10892,139.0396,    0., 1,  ! 243
     .  'B___56' ,'Balise S','AMBA  56',-67.13006,139.0164,    0., 1,  ! 244
     .  'B___57' ,'Balise S','AMBA  57',-67.15158,138.9922,    0., 1,  ! 245
     .  'B___58' ,'Balise S','AMBA  58',-67.17299,138.9688,    0., 1,  ! 246
     .  'B___59' ,'Balise S','AMBA  59',-67.19441,138.9449,    0., 1,  ! 247
     .  'B___60' ,'Balise S','AMBA  60',-67.21529,138.9217,    0., 1,  ! 248
     .  'B___61' ,'Balise S','AMBA  61',-67.23701,138.8970,    0., 1,  ! 249
     .  'B___62' ,'Balise S','AMBA  62',-67.25777,138.8741,    0., 1,  ! 250
     .  'B___63' ,'Balise S','AMBA  63',-67.27934,138.8516,    0., 1,  ! 251
     .  'B___64' ,'Balise S','AMBA  64',-67.29985,138.8272,    0., 1,  ! 252
     .  'B___65' ,'Balise S','AMBA  65',-67.32085,138.8031,    0., 1,  ! 253
     .  'B___66' ,'Balise S','AMBA  66',-67.34199,138.7810,    0., 1,  ! 254
     .  'B___67' ,'Balise S','AMBA  67',-67.36324,138.7575,    0., 1,  ! 255
     .  'B___68' ,'Balise S','AMBA  68',-67.38330,138.7339,    0., 1,  ! 256
     .  'B___69' ,'Balise S','AMBA  69',-67.39895,138.6845,    0., 1,  ! 257
     .  'B___70' ,'Balise S','AMBA  70',-67.41148,138.6308,    0., 1,  ! 258
     .  'B___71' ,'Balise S','AMBA  71',-67.42337,138.5791,    0., 1,  ! 259
     .  'B___72' ,'Balise S','AMBA  72',-67.43371,138.5271,    0., 1,  ! 260
     .  'B___73' ,'Balise S','AMBA  73',-67.44520,138.4770,    0., 1,  ! 261
     .  'B___74' ,'Balise S','AMBA  74',-67.45652,138.4266,    0., 1,  ! 262
     .  'B___75' ,'Balise S','AMBA  75',-67.46742,138.3745,    0., 1,  ! 263
     .  'B___76' ,'Balise S','AMBA  76',-67.48021,138.3261,    0., 1,  ! 264
     .  'B___77' ,'Balise S','AMBA  77',-67.49033,138.2753,    0., 1,  ! 265
     .  'B___78' ,'Balise S','AMBA  78',-67.50198,138.2254,    0., 1,  ! 266
     .  'B___79' ,'Balise S','AMBA  79',-67.51268,138.1755,    0., 1,  ! 267
     .  'B___80' ,'Balise S','AMBA  80',-67.52428,138.1233,    0., 1,  ! 268
     .  'B___81' ,'Balise S','AMBA  81',-67.53552,138.0752,    0., 1,  ! 269
     .  'B___82' ,'Balise S','AMBA  82',-67.54633,138.0226,    0., 1,  ! 270
     .  'B___83' ,'Balise S','AMBA  83',-67.55824,137.9738,    0., 1,  ! 271
     .  'B___84' ,'Balise S','AMBA  84',-67.56938,137.9211,    0., 1,  ! 272
     .  'B___85' ,'Balise S','AMBA  85',-67.57919,137.8671,    0., 1,  ! 273
     .  'B___86' ,'Balise S','AMBA  86',-67.59188,137.8209,    0., 1,  ! 274
     .  'B___87' ,'Balise S','AMBA  87',-67.60346,137.7698,    0., 1,  ! 275
     .  'B___88' ,'Balise S','AMBA  88',-67.61475,137.7186,    0., 1,  ! 276
     .  'B___89' ,'Balise S','AMBA  89',-67.62682,137.6693,    0., 1,  ! 277
     .  'B___90' ,'Balise S','AMBA  90',-67.64194,137.5968,    0., 1,  ! 278
     .  'B___91' ,'Balise S','AMBA  91',-67.66007,137.5244,    0., 1/  ! 279


C +--DDU-Dome-C (source table 5 these Philippe Ciais,
C +  ----------    ALT. table III Pettre al,1986,J.Glac: D03 J18 L24 K33 D59)

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=280,295)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'D13___' ,'D13     ','        ',-66.71, 139.78,  135., 1,  ! 280
     .     'D18___' ,'D18     ','        ',-66.72, 139.77,  421., 1,  ! 281
     .     'D24___' ,'D24     ','        ',-66.75, 139.59,  565., 1,  ! 282
     .     'K33___' ,'K33     ','        ',-66.79, 139.44,  690., 1,  ! 283
     .     'D59___' ,'D59     ','        ',-68.34, 137.32, 2228., 1,  ! 284
     .     'R80___' ,'R80     ','        ',-70.03, 134.84, 2487., 1,  ! 285
     .     'D19___' ,'D19     ','        ',-66.73, 139.69,  440., 1,  ! 286
     .     'D20___' ,'D20     ','        ',-66.74, 139.67,  470., 1,  ! 287
     .     'D21___' ,'D21     ','        ',-66.74, 139.65,  500., 1,  ! 288
     .     'D22___' ,'D22     ','        ',-66.74, 139.63,  523., 1,  ! 289
     .     'D23___' ,'D23     ','        ',-66.75, 139.61,  556., 1,  ! 290
     .     'D32___' ,'D32     ','        ',-66.79, 139.44,  662., 1,  ! 291
     .     'D34___' ,'D34     ','        ',-66.80, 139.41,  709., 1,  ! 292
     .     'D36___' ,'D36     ','        ',-66.81, 139.38,  722., 1,  ! 293
     .     'E39___' ,'E39     ','        ',-66.82, 139.32,  759., 1,  ! 294
     .     'D40___' ,'D40     ','        ',-66.83, 139.31,  783., 1/  ! 295


C +--TNB
C +  ---

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=296,310)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'Alfa__' ,'TNB R_W ','       0',-74.686,164.097,    0., 1,  ! 296
     .     'Bravo_' ,'TNB R_W ','       0',-74.686,164.097,    0., 1,  ! 297
     .     'Eneide' ,'TNB NFH ','WMO89662',-74.696,164.092,   92., 1,  ! 298
     .     'Giulia' ,'Plateau ','WMO89648',-75.536,145.859, 2509., 1,  ! 299
     .     'Irene_' ,'Plateau ','WMO89646',-71.652,148.656, 2094., 1,  ! 300
     .     'Itase_' ,'Plateau ','       0',-72.812,159.298, 2230., 1,  ! 301
     .     'Lola__' ,'TNB N_H ','       0',-74.135,163.431, 1621., 1,  ! 302
     .     'Lucia_' ,'Nansen  ','       0',-74.950,161.771,10000., 1,  ! 303
     .     'Maria_' ,'TNB B_P ','       0',-74.626,164.011,  355., 1,  ! 304
     .     'Minni_' ,'TNB B_P ','       0',-74.633,163.850,  166., 1,  ! 305
     .     'Modest' ,'TNB N_H ','WMO89659',-73.639,160.646, 1924., 1,  ! 306
     .     'Paola_' ,'Reeves  ','       0',-72.767,159.039, 2385., 1,  ! 307
     .     'Pengui' ,'TNB NE  ','       0',-74.333,165.133,   30., 1,  ! 308
     .     'SofiaB' ,'Reeves  ','       0',-75.612,158.590, 1552., 1,  ! 309
     .     'Larsen' ,'Flowcp  ','       0',-74.950,161.767, 1350., 1/  ! 310


C +--EUR (first occurrence, source: http://weather.uwyo.edu/upperair/sounding.html)
C +  ---

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=311,322)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'Uccle_' ,'Uccle   ','        ',  50.80,   4.35,  100., 1,  ! 311
     .     'Trappe' ,'Trappes ','        ',  48.76,   2.00,  168., 1,  ! 312
     .     'Brest_' ,'Brest   ','        ',  48.45,  -4.41,   95., 1,  ! 313
     .     'Bordea' ,'Bordeaux','Merignac',  44.83,  -0.68,   61., 1,  ! 314
     .     'Nimes_' ,'Nimes Co','urbessac',  43.86,   4.40,   62., 1,  ! 315
     .     'Cuneo_' ,'Cuneo Le','valdigi ',  44.53,   7.61, 0000., 2,  ! 316
     .     'Milano' ,'Milano  ','        ',  45.83,   9.28,  103., 1,  ! 317
     .     'Payern' ,'Payerne ','        ',  46.81,   6.95,  491., 1,  ! 318
     .     'Idar_O' ,'Idar Obe','rstein  ',  49.70,   7.33,  377., 1,  ! 319
     .     'Essen_' ,'Essen   ','        ',  51.40,   6.96,  153., 1,  ! 320
     .     'DeBilt' ,'De Bilt ','        ',  52.10,   5.18,    4., 1,  ! 321
     .     'Herstm' ,'Herstmon','ceux    ',  59.90,   0.32,   52., 2/  ! 322
C +         |          |         |             |       |        |  |
C +         |          |         |             |       |        |  v_
C +         |          |         |             |       |        |  OUTPUT = 0: All    OUTPUT are    prohibited 
C +         |          |         |             |       |        |  OUTPUT = 1: netcdf OUTPUT decided in AWSloc
C +         |          |         |             |       |        |  OUTPUT = 2: netcdf OUTPUT decided in AWSloc
C +         |          |         |             |       |        |              ASCII  OUTPUT decided in AWSloc
C +         |          |         |             |       |        v                           (see #WV in SISVAT)
C +         |          |         |             |       v        ALTITUDE of the Station
C +         |          |         |             v       LONGITUDE         of the Station
C +         |          v         v             LATITUDE                  of the Station
C +         v          ATTRIBUTE of the Station, will be written in a title of the corresponding netcdf file
C +         LABEL of the Station,    will be used as the first 3 characters of the corresponding netcdf file


C +--ANT (Ice Cores)
C +  ---

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=323,327)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'TalosD' ,'TALDICE ','        ', -72.78, 159.07, 2315., 2,  ! 323
     .     'Taylor' ,'Taylor D','ome     ', -77.78, 158.72, 2365., 2,  ! 324
     .   'D1____' ,'D1  AWS ','        ',-66.6904219,139.905025,10.,2, ! 325
     .   'D3____' ,'D3  AWS ','        ',-66.6934221,139.895031,10.,2, ! 326
     .   'D17___' ,'D17 AWS ','        ',-66.7212943,139.716433,30.,2/ ! 327


C +--PAPRIKA 
C +  -------

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=328,n_AWS)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'PYR___' ,'Pyramide','        ',  27.96,  86.81, 5035., 2,  ! 328
     .     'PHE___' ,'Pheriche','        ',  27.9 ,  86.82, 4260., 2,  ! 329
     .     'NAM___' ,'Namche  ','        ',  27.8 ,  86.71, 3570., 2,  ! 330
     .     'LUK___' ,'Agromt L','ukla    ',  27.7 ,  86.72, 2660., 2,  ! 331
     .     'Id1006' ,'Gumthang','        ',  27.87,  85.87, 2000., 2,  ! 332
     .     'Id1009' ,'Chautara','        ',  27.78,  85.72, 1660., 2/  ! 333





C +  *******
C +--CAUTION: DO'NT FORGET TO MODIFY THE parameter n_AWS in AWSloc       IF YOU ADD NEW STATIONS!
C +  *******                                                AWS_nc
C +                                                         AWS_nc_DATA 
 
      end
