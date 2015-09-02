      subroutine SBCnew

C +------------------------------------------------------------------------+
C | MAR INPUT      SURFACE BOUNDARY LAYER                  08-03-2005  MAR |
C |                                                                        |
C | Purpose: Water Mass Balance                                            |
C |          Free  Response of the Temperature to the Forcing              |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--Global Variables
C +  ================

      include 'MARCTR.inc'
      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_RA.inc'
      include 'MAR_DY.inc'
      include 'MAR_HY.inc'
      include 'MAR_CA.inc'
      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MARdSV.inc'
      include 'MAR0SV.inc'
      include 'MAR_TV.inc'

      include 'MAR_LB.inc'

      integer           ipr2nc,npr2nc,n_iter
      common/Out2nc_arg/ipr2nc,npr2nc,n_iter


C +--Saved  Variables
C +  ================

      INTEGER ipt         ,i_1D        ,j_1D
      INTEGER IOi_1D(iptx),IOj_1D(iptx)
C +...        IO Grid                     Indices

      INTEGER isol1D(  1,  1)
C +...        isol1D:     Soil       Type Index

      INTEGER iWaF1D(  1,  1)
C +...       (iWaF1D=0 ==> no Water Flux;
C +           iWaF1D=1 ==> free drainage)

      REAL    AlbS1D(  1,  1)
C +...        AlbS1D: Dry Soil       Albedo

      INTEGER iveg1D(  1,  1,nvx)
C +...        iveg1D:     Vegetation Type Index

      INTEGER ifra1D(  1,  1,nvx)
C +...        ifra1D:     Vegetation Class Coverage
C +                               (3 Class, Last One is Open Water)

      REAL    alai1D(  1,  1,nvx)
C +...        alai1D:           Leaf Area Index

      REAL    glf_1D(  1,  1,nvx)
C +...        glf_1D:     Green Leaf Fraction

      REAL    Tsol1D(  1,  1,nvx,llx)
C +...        Tsol1D: Layer Soil       Temperature

      REAL    eta_1D(  1,  1,nvx,llx)
C +           eta_1D: Soil                Moisture Content

      REAL    CaSn1D(  1,  1,nvx)
C +...        CaSn1D: Canopy     Intercepted Snow  Content

      REAL    CaWa1D(  1,  1,nvx)
C +...        CaWa1D: Canopy     Intercepted Water Content

      REAL    psiv1D(  1,  1,nvx)
C +           psiv1D: Vegetation         Hydraulic Potential

      REAL    psig1D(  1,  1,nvx)
C +           psig1D: Soil               Hydraulic Potential

      REAL    Tveg1D(  1,  1,nvx)
C +...        Tveg1D: Skin  Vegetation Temperature

      REAL    Tgrd1D(  1,  1,nvx)
C +...        Tgrd1D: Skin  Soil       Temperature

      REAL    evap1D(  1,  1)
C +           evap1D: Time Integrated Evapotranspiration

      REAL    drai1D(  1,  1)
C +           drai1D: Time Integrated Drainage     Flow

      REAL    runo1D(  1,  1)
C +           runo1D: Time Integrated (Sub)surface Flow


C +--Local  Variables
C +  ================

      real              SrfOLD,Water0,dT_UpW,AfwLat,z__msl(mz)
      common/SBCnew_rea/SrfOLD,Water0,dT_UpW,AfwLat

      integer           nperpe,ihs   ,ihn
      common/SBCnew_int/nperpe,ihs   ,ihn

      integer           No_Day,n     ,l     ,lm1   ,lp1   ,isl
      real              SSTlat,radMed,arglat
      real              SrfINP,SrfNEW,SrfOUT,BLeft ,BRight
      real              WaterF,Water1,Waterd,rain_C,rain_S
      real              waterC,depthS

      logical           obsSST,fixSST,zzzzav,UpWell,fix_LB,freeLB
      real              Tsolav,Height,Raylei,rfrict
      real              TimeNC


C +--DATA
C +  ====

      data              obsSST/.true./
      data              fixSST/.false./
      data              zzzzav/.true./
      data              UpWell/.true./
      data              fix_LB/.false./
      data              freeLB/.true./
      data              Raylei/0.100e-6/
      data              TimeNC/10800.e0/
C +...                  TimeNC:Time Interval between IO on NetCDF File

      include          'SBCnew_1D.ctr'



C +--Prescribed SST
C +  ==============

      include          'SBCnewSST.ctr'

C +                     ******
      if (obsSST)  call SST_in
C +                     ******


C +--SST "Upwelling" Anomaly Extent
C +  ==============================

      AfwLat = 6. * degrad

      IF (iterun.eq.0)                                            THEN ! CTR


C +--SISVAT Forcing   File
C +  =====================

          open (unit=33,status='unknown',file='SISVAT_IN.DAT',
     .          form='unformatted')
          rewind     33


C +--Animation NetCDF File
C +  =====================

          n_iter = TimeNC / dt
C +...    n_iter : Number of Time steps between IO

          ipr2nc =                   1
          npr2nc = nterun / n_iter + 1
C +...    npr2nc : Number of                    IO

          dt_Loc = dt

C +       ***********
          call Out2nc
C +       ***********

      ELSE                                                             ! CTR
        IF       (mod(iterun,n_iter).eq.0)                        THEN ! CTR
                      ipr2nc=ipr2nc  +  1

          dt_Loc = dt

C +       ***********
          call Out2nc
C +       ***********

        END IF                                                         ! CTR
      END IF                                                           ! CTR


C +--SISVAT 1D Init.  File
C +  =====================

      IF (iterun.eq.1)                                            THEN ! CTR

          write(33) maskSL(i_1D,j_1D)
          write(33)     sh(i_1D,j_1D)   

          open (unit=11,status='unknown',file='SISVAT_1D.DAT',
     .          form='unformatted')
          rewind     11

          WRITE(11) itexpe
          WRITE(11) iyrrGE,mmarGE,jdarGE,jhurGE

          DO ipt=1 ,iptx
             IOi_1D(ipt)     = i_1D
             IOj_1D(ipt)     = j_1D
          END DO
             isol1D(1,1)     = isolTV(i_1D,j_1D)
             iWaF1D(1,1)     = iWaFTV(i_1D,j_1D)
             AlbS1D(1,1)     = AlbSTV(i_1D,j_1D)
          DO n  =1 ,nvx
             iveg1D(1,1,n)   = ivegTV(i_1D,j_1D,n)
             ifra1D(1,1,n)   = ifraTV(i_1D,j_1D,n)
             alai1D(1,1,n)   = alaiTV(i_1D,j_1D,n)
             glf_1D(1,1,n)   = glf_TV(i_1D,j_1D,n)
          DO l  =1 ,llx
             Tsol1D(1,1,n,l) = TsolTV(i_1D,j_1D,n,l)
             eta_1D(1,1,n,l) = eta_TV(i_1D,j_1D,n,l)
          END DO
             CaSn1D(1,1,n)   = CaSnTV(i_1D,j_1D,n)
             CaWa1D(1,1,n)   = CaWaTV(i_1D,j_1D,n)
             psiv1D(1,1,n)   = psivTV(i_1D,j_1D,n)
             psig1D(1,1,n)   = psigTV(i_1D,j_1D,n)
             Tveg1D(1,1,n)   = TvegTV(i_1D,j_1D,n)
             Tgrd1D(1,1,n)   = TgrdTV(i_1D,j_1D,n)
          END DO
             evap1D(1,1)     = evapTV(i_1D,j_1D)
             drai1D(1,1)     = draiTV(i_1D,j_1D)
             runo1D(1,1)     = runoTV(i_1D,j_1D)

          WRITE(11) IOi_1D  ! IO i Index
          WRITE(11) IOj_1D  ! IO j Index
          WRITE(11) isol1D  ! Soil       Type Index
          WRITE(11) iWaF1D  ! =0 ==> no Water Flux
C +                         ! =1 ==> free Drainage
          WRITE(11) AlbS1D  ! Dry Soil       Albedo
          WRITE(11) iveg1D  ! Vegetation Type Index
          WRITE(11) ifra1D  ! Vegetation Class Coverage
          WRITE(11) alai1D  !       Leaf Area Index                [-]
          WRITE(11) glf_1D  ! Green Leaf Fraction                  [-]
          WRITE(11) Tsol1D  ! Soil Temperature                     [K]
          WRITE(11) eta_1D  ! Soil Moisture                    [kg/kg]
C +                         !
          WRITE(11) CaSn1D  ! Canopy Intercepted Snow  Content   [mWE]
          WRITE(11) CaWa1D  ! Canopy Intercepted Water Content [kg/m2]
          WRITE(11) psiv1D  ! Vegetation   Hydraulic Potential     [m]
          WRITE(11) psig1D  ! Ground       Hydraulic Potential     [m]
          WRITE(11) Tveg1D  ! Skin      Vegetation Temperature     [K]
          WRITE(11) Tgrd1D  ! Skin            Soil Temperature     [K]
C +                         !
          WRITE(11) evap1D  ! Total Evapotranspiration    [mmWE]
          WRITE(11) drai1D  ! Drainage            Flow    [mm/s]
          WRITE(11) runo1D  ! Integrated Drainage Flow      [mm]

          close(unit=11)

             write(6,600) 
 600         format(/,'    z     |   dz     |   eta       |',
     .              /,'   [m]    |   [m]    |   [-]       |',
     .              /,'----------+----------+-------------+')
             waterC = 0.
             depthS = 0.
          DO l  =1 ,llx
             waterC = waterC + eta_1D(1,1,1,l)*deptTV(l)
             depthS = depthS +                 deptTV(l)
             write(6,601) depthS,deptTV(l),eta_1D(1,1,1,l)
 601         format(2(f9.3,' |'),f12.6,' |')
          END DO
             waterC = waterC * 1.e3
             write(6,602) i_1D,isolSL(i_1D,j_1D),waterC
 602         format('----------+----------+-------------+',
     .           //,'Grid Point',i6,' (',i1,')',
     .                '  Water Content:',f9.3,' mmWE',/)
      END IF                                                           ! CTR
C +
C +
C +--Water Conservation
C +  ==================
C +
       IF (iterun.eq.     0)
     .   open( unit=31,file='Sbcnew.HYD',status='unknown')
C +
C +
C +--Contribution from Surface Fluxes
C +  --------------------------------
C +
          SrfINP = 0.0d+0
          SrfNEW = 0.0d+0
          j=jmez
       DO i=mx/8+1,7*mx/8
          SrfINP = SrfINP + hlatSL(i,j) * dx * dt     / Lv_H2O
          SrfNEW = SrfNEW + rainHY(i,j) * dx * 1.0d+3
       END DO
          SrfOUT = SrfNEW - SrfOLD
          SrfOLD = SrfNEW
C +
C +
C +--Contribution from Lateral Fluxes
C +  --------------------------------
C +
          BLeft  = 0.0d+0
          BRight = 0.0d+0
          j=jmez
          i=mx/8+1
       DO k=1,mz
          BLeft  = BLeft  
     .      + dt * 1.0d+3 * pstDYn(i,j)    *dsigm1(k) *uairDY(i,j,k)
     .           *           (qvDY(i,j,k)  +qwHY(i,j,k) +qrHY(i,j,k)
     .                                     +qiHY(i,j,k) +qsHY(i,j,k))
       END DO
          i=7*mx/8
       DO k=1,mz
          BRight = BRight 
     .      + dt * 1.0d+3 * pstDYn(i,j)    *dsigm1(k) *uairDY(i,j,k)
     .           *           (qvDY(i,j,k)  +qwHY(i,j,k) +qrHY(i,j,k)
     .                                     +qiHY(i,j,k) +qsHY(i,j,k))
       END DO
C +
C +
C +--Total Contribution
C +  ------------------
C +
          WaterF = BLeft  - BRight + SrfINP - SrfOUT
C +
C +
C +--Variation
C +  ---------
C +
          Water1 = 0.0d+0
          j=jmez
       DO i=mx/8+1,7*mx/8
       DO k=1,mz 
          Water1 = Water1 
     .      +     1.0d+3 * pstDYn(i,j)    *dsigm1(k) *dx
     .           *           (qvDY(i,j,k)  +qwHY(i,j,k) +qrHY(i,j,k)
     .                                     +qiHY(i,j,k) +qsHY(i,j,k))
       END DO
       END DO
       IF(itexpe.eq.    0)  Water0 = Water1
          Waterd = Water1 - Water0
          Water0 = Water1
C +
C +
C +--OUTPUT
C +  ------
C +
       IF (jsecGE.eq.0)                                            THEN
         IF               (minuGE.eq.0)                            THEN
             write(31,310)
 310         format(
     .       /, '-------------------+',
     .       1('-------+-------+-------+',2('-------+')),
     .       9('-------+'),
     .       /, 'JJMMMYYYY HH.mm.ss |',
     .       1(' Tair  | HSENS | HLATE | RR CA | RR HY |'),
     .         ' Av CA | Av HY | Water | Wat_D | Wat_F |',
     .         ' Wat_L | Wat_R | WatSI | WatSO |',
     .       /, '                   |',
     .       1(' [K]   | [W/m2]| [W/m2]|',2(' [mm]  |')),
     .       9(' [mm]  |'),
     .       /, ' ------------------+',
     .       1('-------+-------+-------+',2('-------+')),
     .       9('-------+'))
         END IF
C +
              rain_C = 0.0d0
              rain_S = 0.0d0
            j=jmez
         DO i=mx/8+1,7*mx/8
              rain_C = rain_C+            rainCA(i,j)
              rain_S = rain_S+rainHY(i,j)-rainCA(i,j)
         END DO
              rain_C = rain_C*4/(3*mx)
              rain_S = rain_S*4/(3*mx)
C +
             write(31,311)jdarGE,labmGE(mmarGE),iyrrGE,jhlrGE(imez,jmez)
     .                   ,minuGE,jsecGE,    
     .       (tairDY(i,j,mz)  , hsenSL(i,j),hlatSL(i,j),
     .        rainCA(i,j)*1.d3,(rainHY(i,j)-rainCA(i,j))*1.d3,
     .               i=4*mx/8,4*mx/8),
     .        rain_C     *1.d3, rain_S                  *1.d3,
     .        Water1/(dx*mx*3/4),Waterd/(dx*3*mx/4),WaterF/(dx*3*mx/4),
     .        BLeft /(dx*mx*3/4),BRight/(dx*3*mx/4),
     .        SrfINP/(dx*mx*3/4),SrfOUT/(dx*3*mx/4)
 311         format(i2,a3,i4,i3,'.',i2,'.',i2,' |',
     .       1(3(f6.1,' |'),2(f6.2,' |')),2(f6.2,' |'),
     .           f6.1,' |' ,7(f6.3,' |'))
       END IF
C +
C +
C +--1D SISVAT Forcing
C +  -----------------
C +
         write(33) uairDY(i_1D,j_1D,mz),vairDY(i_1D,j_1D,mz)
     .            ,tairDY(i_1D,j_1D,mz),  qvDY(i_1D,j_1D,mz)
     .            ,gplvDY(i_1D,j_1D,mz),pstDYn(i_1D,j_1D)
     .            ,rainHY(i_1D,j_1D)   ,czenGE(i_1D,j_1D)
     .            ,RAdsol(i_1D,j_1D)   ,RAd_ir(i_1D,j_1D)
     .            ,cld_SL(i_1D,j_1D)   ,RAd_ir(i_1D,j_1D)
C +
       IF (iterun.eq.nterun)                                        THEN
         close(unit=31)
         close(unit=33)
       END IF
C +
C +
C +--Perpetual Mode
C +  ==============
C +
       IF  (iterun.eq.0)                                            THEN
         open (unit=30,status='old',file='SBCnew.ctr')
         rewind     30
         read (     30,300) nperpe
 300     format(i6)
         read (     30,301) dT_UpW
 301     format(f6.3)
         read (     30,302) freeLB
         read (     30,302) fix_LB
 302     format(l6)
         close(unit=30)
       END IF
       IF  (jhurGE   .eq.    0.and.minuGE.eq.0.and.jsecGE.eq.0)     THEN
         IF(itexpe*dt.ge.86400.and.nperpe.gt.0                )     THEN
            jda0GE = jda0GE-1
            nperpe = nperpe-1
         END IF
       END IF
C +
C +
C +--SST
C +  ===
C +
      IF (fixSST)                                                   THEN
            No_Day          = itexpe  *dt /86400.
            radMed          = 36.
C +...      radMed          : Mediterranean Sea Southern Border (36 deg N)
C +
        DO j=1,my
        DO i=1,mx
          IF (No_Day.le.2 .AND.       isolSL(i,j) .le. 2)           THEN
              arglat        =(min(abs(GElatr(i,j) /degrad),radMed))
     .                                            *0.5*pi /radMed
              SSTlat        =         sst_SL      -9.0*sin(arglat)
C +...                        9 deg C less over the Mediterranean Sea
C + 
            DO n=1,nvx
            DO l=1,llx
              TsolTV(i,j,n,l) = TsolTV(i,j,n,l)
     .                + (SSTlat-TsolTV(i,j,n,l))*dt/86400.
C +
            END DO
            END DO
          END IF
        END DO
        END DO
      END IF
C +
C +--Height Average (over the ocean only)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (zzzzav.AND.iterun.gt.1)                                   THEN
        DO j=1,my
        DO i=1,mx
          IF    (isolSL(i,j).le.2)                                  THEN
            DO n=1,nvx
                 Tsolav = 0.
                 Height = 0.
                 isl= -nsol
                 l  = 1-isl
                 lp1=  -isl
                 Tsolav = Tsolav
     .                  +(TsolTV(i,j,n,l  ) *dz78SV(isl)
     .                   +TsolTV(i,j,n,lp1) *dz_8SV(isl))
                 Height = Height
     .                  +(                   dz78SV(isl)
     .                   +                   dz_8SV(isl))
C +
              DO isl= -nsol+1,-1
                 l  = 1-isl
                 lp1=  -isl
                 lm1= 2-isl
                 Tsolav = Tsolav
     .                  +(TsolTV(i,j,n,l  ) *dz34SV(isl)
     .                  +(TsolTV(i,j,n,lm1)
     .                   +TsolTV(i,j,n,lp1))*dz_8SV(isl))
                 Height = Height
     .                  +(                   dz34SV(isl)
     .                  + 2.                *dz_8SV(isl))
              END DO
C +
                 isl=     0
                 l  = 1-isl
                 lp1=  -isl
                 lm1= 2-isl
                 Tsolav = Tsolav
     .                  +(TsolTV(i,j,n,l  ) *dz78SV(isl)
     .                  + TsolTV(i,j,n,lm1) *dz_8SV(isl))
                 Height = Height
     .                  +(                   dz78SV(isl)
     .                  +                    dz_8SV(isl))
C +
                 Tsolav = Tsolav            /Height
C +
              DO l=1,llx
                 TsolTV(i,j,n,l) = Tsolav
              END DO
            END DO
          END IF
        END DO
        END DO
      END IF
C +
C +--Upwelling
C +  ~~~~~~~~~
      IF (UpWell)                                                   THEN
        DO j=1,my
        DO i=1,mx
          IF (isolSL(i,j).le.2)                                     THEN
            DO n=1,nvx
            DO l=1,llx
              TsolTV(i,j,n,l) = TsolTV(i,j,n,l) 
     .      - exp(-max(zero,abs(GElatr(i,j))-AfwLat))*dT_UpW*dt/86400.
            END DO
            END DO
          END IF
        END DO
        END DO
      END IF
C +
C +
C +--LB Air Temperature, Pressure Depth Prescription
C +  ===============================================
C +
      IF (fix_LB.AND.iterun.eq.0)                                   THEN
        i   = 0
 110    CONTINUE
        i   = i + 1
        IF (sst_LB(i,1).lt.epsi)                               GO TO 110
        ihs = i
C +
        i   = mx + 1
 120    CONTINUE
        i   = i - 1
        IF (sst_LB(i,1).lt.epsi)                               GO TO 120
        ihn = i
        write(6,160) ihs,sst_LB(ihs,1),ihn,sst_LB(ihn,1)
 160    format(/,' SBCnew: LB temperatures are prescribed from SST',
     .         /,'    SST(',i3,') =',f9.3,
     .         /,'    SST(',i3,') =',f9.3)
        DO k=1,mz
          z__msl(k) = gplvDY(1,1,k)*grvinv
        ENDDO
      END IF
C +
      IF (fix_LB.AND.iterun.gt.1)                                   THEN
C + 
         DO k=1,mz
         DO j=1,my
         DO i=1,n7mxLB
           SSTlat          =(sst_LB(ihs,j)-0.0065*z__msl(k))/pkDY(i,j,k)
           vaxgLB(i,j,k,4) = vaxgLB(i,j,k,4)
     .                    + (SSTlat        -vaxgLB(i,j,k,4))*dt/0.864e6
C +...                                                         10 days
C +
         END DO
         END DO
         END DO
C +
         DO k=1,mz
         DO j=1,my
         DO i=mx-n6mxLB,mx
           SSTlat          =(sst_LB(ihn,j)-0.0065*z__msl(k))/pkDY(i,j,k)
           vaxdLB(i,j,k,4) = vaxdLB(i,j,k,4)
     .                    + (SSTlat        -vaxdLB(i,j,k,4))*dt/0.864e6
         END DO
         END DO
         END DO
C +
      END IF
C +
C +
C +--LB Air Temperature, Pressure Depth Adjustment
C +  =============================================
C +
       IF (freeLB)                                                  THEN
         DO k=1,mz
         DO j=1,my
         DO i=1,n7mxLB
           vaxgLB(i,j,k,4) = pktaDY(  2,j,k)
     .             -3.e-5 *dt /pkDY(  i,j,k)      ! Extratropical Eddies
CC#        vaxgLB(i,j,k,5) = pstdYn(  2,j)        ! Mass blows up
         END DO
         END DO
         END DO
C +
         DO k=1,mz
         DO j=1,my
         DO i=mx-n6mxLB,mx
           vaxdLB(i,j,k,4) = pktaDY(mx1,j,k)
     .             -3.e-5 *dt /pkDY(  i,j,k)      ! Extratropical Eddies
CC#        vaxdLB(i,j,k,5) = pstdYn(mx1,j)        ! Mass blows up
         END DO
         END DO
         END DO
       END IF
C +
C +
C +--Rayleigh Friction
C +  =================
C +
      DO k=1,mzabso
        rfrict          = Raylei * dt * (sigma(mzabso)-sigma(k))
     .                                / (sigma(mzabso)-sigma(1))
      DO j=1,my
      DO i=1,mx
        uairDY(i,j,k)   = uairDY(i,j,k) -rfrict*uairDY(i,j,k)*dt
        vairDY(i,j,k)   = vairDY(i,j,k) -rfrict*vairDY(i,j,k)*dt
      END DO
      END DO
      END DO
C +
C +
C +--OUTPUT
C +  ======
C +
C +--SISVAT Forcing
C +  --------------
C +

C +
C +
C +--SISVAT Soil/Vegetation
C +  ----------------------
C +
      IF (iterun.EQ.1)                                              THEN
        DO i=1,mx
        write(6,6000) 
     .     i,xxkm(i)  ,     sh(i,1)  ,
     .     isolSL(i,1),
     .     isolTV(i,1),(ifraTV(i,1,n),n=1,mw)
 6000   format((i4,2f6.0,2i3,3i3,15('|',3i2)))
        write(6,6001)  (ivegTV(i,1,n),n=1,mw)
 6001   format((4x,12x  ,6x ,3i3,15('|',3i2)))
        END DO
      END IF
C +
      return
      end 


      subroutine Out2nc

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                             16-02-2005  MAR |
C |   SubRoutine Out2nc is used to write the main Model Variables          |
C |                                      on  a NetCDF file                 |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT: ipr2nc: Current time step    number                           |
C |   ^^^^^^         (starting from ipr2nc=1, which => new file creation)  |
C |          npr2nc: Total  'time slices' number (max value of ipr2nc)     |
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
c #NH include 'MAR_NH.inc'
c #OL include 'MAR_OL.inc'

      include 'MAR_TE.inc'
      include 'MAR_TU.inc'

      include 'MAR_HY.inc'
      include 'MAR_CA.inc'
      include 'MAR_RA.inc'

      include "SBCnew.hTC"

      include 'MAR_SL.inc'
c #PO include 'MAR_PO.inc'
c #IB include 'MAR_IB.inc'
      include 'MAR_SV.inc'
      include 'MAR_TV.inc'

      include 'MAR_WK.inc'

      include 'MAR_IO.inc'

      integer           ipr2nc,npr2nc,n_iter
      common/Out2nc_arg/ipr2nc,npr2nc,n_iter


C +--Local   Variables
C +  =================

      integer    Lfnam,     Ltit,     Luni,     Lnam,     Llnam
      PARAMETER (Lfnam= 40, Ltit= 90, Luni= 31, Lnam= 13, Llnam=50)
C +...Length of char strings 

      CHARACTER*(Lfnam)  fnamNC
      common/Out2nc_loc/ fnamNC
C +...                   fnamNC: To retain file name.

      real               rrCV0(mx,my),rrtt0(mx,my)
      common/Out2rr_loc/ rrCV0       ,rrtt0
C +...                   rrCV0 : Convective Rain over Previous Time Interval
C +                      rrtt0 : Total      Rain over Previous Time Interval

      real               T_Soil(mx,my,mw)
      real               q_Soil(mx,my,mw)

      integer    NdimNC
      PARAMETER (NdimNC = 5)
C +...Number of defined spatial dimensions (exact)

      integer    MXdim
      PARAMETER (MXdim =500000)
C +...Maximum Number of all dims: recorded Time Steps
C +   and also maximum of spatial grid points for each direction. 

      integer    MX_var
      PARAMETER (MX_var = 80)
C +...Maximum Number of Variables 
C +
      integer    NattNC
      PARAMETER (NattNC = 2)
C +...Number of real attributes given to all variables
C +
      INTEGER           RCODE
C +
      integer           jourNC(MXdim)
      integer           moisNC(MXdim)
      real              yearNC(MXdim)
      real              dateNC(MXdim)
      real              timeNC(MXdim)
      real              VALdim(MXdim,0:NdimNC)
      integer           nDFdim(      0:NdimNC)
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
C +
      REAL     starti,DayLen,optwa ,optia ,rhodzk
      real     starta
      integer  n1000 ,n100a ,n100  ,n10_a ,n10   ,n1    ,m10
      integer                              jd10  ,jd1
      integer  MMXstp,it    ,mois  ,mill  ,iu    ,n
      integer  itotNC,NtotNC,ID__nc
c #TC integer  itot  
C +
C +
C +--NetCDF File Initialization
C +  ==========================
C +
      IF (ipr2nc.eq.1) THEN
C +
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
C +
C +
C +--Output File Label
C +  -----------------
C +
        fnamNC = 'ANI.'
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)
     .         // '.' // explIO
     .         // '.nc    '
C +
C +
C +--Output Title
C +  ------------
C +
        tit_NC = 'MAR'
     .         // ' - Exp: ' // explIO
     .         // ' - '
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)
C +
C +
C +--Create File / Write Constants
C +  -----------------------------
        MMXstp = MXdim
C +...  To check array bounds... silently
C +
C +--Time Variable (hour)
C +  ~~~~~~~~~~~~~~~~~~~~
C +
C +...  To define a NetCDF dimension (size, name, unit):
c _UL   nDFdim(0)= npr2nc
        nDFdim(0)= 0
        NAMdim(0)= 'time'
        UNIdim(0)= 'HOURS since 1901-01-15 00:00:00'
C +
C +...  Check temporary arrays: large enough ?
        IF (npr2nc.gt.MMXstp)
     &  STOP '*** Out2nc - ERROR : MXdim to low ***'
C +
              starti = jhurGE + minuGE/60.d0 + jsecGE/3600.d0
C +...        starti : Starting Time (= current time in the day)
C +
              starta = (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
C +
        DO it = 1,npr2nc
              timeNC(it)   = starti + (it-1) * n_iter  *dt / 3600.
C +...                                         n_iter: #iter between output
C +
              VALdim(it,0) = starta + (it-1) * n_iter  *dt / 3600.
C +...        VALdim(  ,0) : values of the dimension # 0 (time) 

C +--Time Variable (date)
C +  ~~~~~~~~~~~~~~~~~~~~
              dateNC(it) =          timeNC(it)
              jourNC(it) = jdarGE + timeNC(it) / 24.d0
        END DO
                  mois       =  mmarGE
                  mill       =  iyrrGE
        DO it = 1,npr2nc
          IF     (jourNC(it).gt.njmoGE(mois))                     THEN ! CTR
            DO iu=it,npr2nc
                  jourNC(iu) =  jourNC(iu) - njmoGE(mois)
            END DO
                  mois       =  mois + 1
              IF (mois.gt.12)                                     THEN ! CTR
                  mois       =         1
                  mill       =  mill + 1
              END IF                                                   ! CTR
          END IF                                                       ! CTR
                  moisNC(it) =  mois
                  yearNC(it) =  mill
C +
          IF     (dateNC(it).gt.24.d0-epsi)                       THEN ! CTR
                  DayLen     =  24.d0
            DO iu=it,npr2nc
                  dateNC(iu) = mod(dateNC(iu),DayLen)
            END DO
          END IF                                                       ! CTR
        END DO
C +
        DO it = 1,npr2nc
              dateNC(it) =  dateNC(it)
     .             + 1.d+2 *jourNC(it)
     .             + 1.d+4 *moisNC(it)
        END DO
C +

C     Define horizontal spatial dimensions :    
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +
C +...  Check temporary arrays: large enough ?
        IF (    mx .gt.MMXstp.or.my.gt.MMXstp
     &      .or.mzz.gt.MMXstp.or.mw.gt.MMXstp)
     &    STOP '*** Out2nc - ERROR : MXdim to low ***'
C +
C +...To define NetCDF dimensions (size, name, unit):
C +
        DO i = 1, mx
          VALdim(i,1) = xxkm(i)
        END DO
          nDFdim(1)   = mx
          NAMdim(1)   = 'x'
          UNIdim(1)   = 'km'
C +
        DO j = 1, my
          VALdim(j,2) = yykm(j)
        END DO
          nDFdim(2)   = my
          NAMdim(2)   = 'y'
          UNIdim(2)   = 'km'
C +
        do k = 1, mz
          VALdim(k,3) = sigma(k)
        enddo
          nDFdim(3)   = mz
          NAMdim(3)   = 'level'
          UNIdim(3)   = '[sigma]'
C +...    For levels k
C +
        do k = 1, mz
          VALdim(k,4) = sigmid(k+1)
        enddo
          nDFdim(4)   = mz
          NAMdim(4)   = 'level2'
          UNIdim(4)   = '[sigma]'
C +...    For levels k+1/2
C +
        do k = 1, mw
          VALdim(k,5) = k 
        enddo
          nDFdim(5)   = mw
          NAMdim(5)   = 'sector'
          UNIdim(5)   = '[index]'
C +...    For Surface Sectors 
C +
C +--Variable's Choice (Table MARgou.dat)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +
        OPEN(unit=10,status='unknown',file='MARgou.dat')
C +
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
C +
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
C +
          ENDIF
        GOTO 980
 990    CONTINUE
C +
        CLOSE(unit=10)
C +
        NtotNC = itotNC 
C +...  NtotNC : Total number of variables writen in NetCDF file.
C +
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
C +
C +--Automatic Generation of the NetCDF File Structure
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +
C +     **************
        CALL UNscreate (fnamNC,tit_NC,
     &                  NdimNC, nDFdim, MXdim , NAMdim, UNIdim, VALdim,
     &                  MX_var, NtotNC, nameNC, SdimNC, unitNC, lnamNC,
     &                  NattNC, NAMrat, NvatNC,
     &                  ID__nc) 
C +     **************
C +
C +
C +--Write Time - Constants
C +  ~~~~~~~~~~~~~~~~~~~~~~
        DO j=1,my
        DO i=1,mx
          Wkxy1(i,j) =  GElonh(i,j) * 15.d0
C +...    Conversion: Hour->degrees
C +
          WKxy2(i,j) =  GElatr(i,j) / degrad
C +...    Conversion: rad ->degrees
C +
          WKxy3(i,j) =  isolSL(i,j)
C +...    Conversion to REAL type (integer not allowed)
C +
        END DO
        END DO
C +
        IF (nDFdim(0).ne.0)                                         THEN
C +
C +       ************
          CALL UNwrite (ID__nc, 'year  ', 1  , npr2nc, 1 , 1 , yearNC)
          CALL UNwrite (ID__nc, 'date  ', 1  , npr2nc, 1 , 1 , dateNC)
C +       ************
C +
        END IF
C +
C +       ************
          CALL UNwrite (ID__nc, 'lon   ', 1  , mx    , my, 1 , Wkxy1)
          CALL UNwrite (ID__nc, 'lat   ', 1  , mx    , my, 1 , Wkxy2)
          CALL UNwrite (ID__nc, 'sh    ', 1  , mx    , my, 1 , sh)
          CALL UNwrite (ID__nc, 'isol  ', 1  , mx    , my, 1 , Wkxy3)
          CALL UNwrite (ID__nc, 'S_Frac', 1  , mx    , my, mw, SLsrfl)
C +       ************
C +
C +--Re-Open file if already created.
C +  ================================
C +
C +
      ELSE
C +
C +     ************
        CALL UNwopen (fnamNC,ID__nc)
C +     ************
C +
      END IF
C +
C +
C +--Write Time-dependent variables:
C +  ===============================
C +
C +
C +--UNLIMITED Time Dimension
C +  ------------------------
C +
      IF (nDFdim(0).eq.0)                         THEN !
              starta = (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
C +
C +     ************
        CALL UNwrite (ID__nc, 'time   ', ipr2nc,  1,  1,  1, starta)
C +     ************
C +
      END IF
C +
C +
C +--Dynamics
C +  --------
C +
      do j=1,my
      do i=1,mx
        WKxy1(i,j) = uairDY(i,j,mz)
        WKxy2(i,j) = vairDY(i,j,mz)
        WKxy3(i,j) = TairDY(i,j,mz)
        WKxy4(i,j) =   qvDY(i,j,mz)
        WKxy6(i,j) = rainCA(i,j) - rrCV0(i,j)
        rrCV0(i,j) = rainCA(i,j)
        WKxy7(i,j) = rainHY(i,j) - rrtt0(i,j)
        rrtt0(i,j) = rainHY(i,j)
        DO n=1,mw
          T_Soil(i,j,n) = TsolTV(i,j,n,1)
          q_Soil(i,j,n) = eta_TV(i,j,n,1)
        ENDDO
      enddo
      enddo
C +
C +   ************
      CALL UNwrite (ID__nc, 'ua_SBL ', ipr2nc, mx, my, 1 , WKxy1 )
      CALL UNwrite (ID__nc, 'va_SBL ', ipr2nc, mx, my, 1 , WKxy2 )
      CALL UNwrite (ID__nc, 'Ta_SBL ', ipr2nc, mx, my, 1 , WKxy3 )
      CALL UNwrite (ID__nc, 'qa_SBL ', ipr2nc, mx, my, 1 , WKxy4 )
      CALL UNwrite (ID__nc, 'Ts_SRF ', ipr2nc, mx, my, 1 , TairSL)
      CALL UNwrite (ID__nc, 'RadSol ', ipr2nc, mx, my, 1 , RAdsol)
      CALL UNwrite (ID__nc, 'Sol_SL ', ipr2nc, mx, my, 1 , Sol_SL)
      CALL UNwrite (ID__nc, 'alb0SL ', ipr2nc, mx, my, 1 , alb0SL)
      CALL UNwrite (ID__nc, 'albeSL ', ipr2nc, mx, my, 1 , albeSL)
      CALL UNwrite (ID__nc, 'albsTV ', ipr2nc, mx, my, 1 , albsTV)
      CALL UNwrite (ID__nc, 'Rad_IR ', ipr2nc, mx, my, 1 , RAD_ir)
      CALL UNwrite (ID__nc, 'hsenSL ', ipr2nc, mx, my, 1 , hsenSL)
      CALL UNwrite (ID__nc, 'hlatSL ', ipr2nc, mx, my, 1 , hlatSL)
      CALL UNwrite (ID__nc, 'uqstar ', ipr2nc, mx, my, mw, SLuqsl)
      CALL UNwrite (ID__nc, 'RainCV ', ipr2nc, mx, my, 1 , WKxy6 )
      CALL UNwrite (ID__nc, 'Raintt ', ipr2nc, mx, my, 1 , WKxy7 )
      CALL UNwrite (ID__nc, 'DDraft ', ipr2nc, mx, my, 1 , WKxy8 )
      CALL UNwrite (ID__nc, 'T_Soil ', ipr2nc, mx, my, mw, T_Soil)
      CALL UNwrite (ID__nc, 'q_Soil ', ipr2nc, mx, my, mw, q_Soil)
C +   ************
C +
C +
C +--Work Arrays Reset
C +  -----------------
C +
      do j=1,my
      do i=1,mx
        WKxy1(i,j)   =0.0
        WKxy2(i,j)   =0.0
        WKxy3(i,j)   =0.0
        WKxy4(i,j)   =0.0
        WKxy5(i,j)   =0.0
        WKxy6(i,j)   =0.0
        WKxy7(i,j)   =0.0
        WKxy8(i,j)   =0.0
      enddo
      enddo
C +
C +
C +--Insolation
C +  ==========
C +
      do j=1,my
      do i=1,mx
           WKxy5(i,j)  = czenGE(i,j) * rsunGE
      enddo
      enddo
C +
C +
C +--Cloud Optical Depth
C +  ===================
C +
      do j=1,my
      do i=1,mx
           optwa  = zero
C +...     optwa  : liquid water path (kg/m2) (droplets)
C +
           optia  = zero
c +...     optia  : liquid water path (kg/m2) (crystals)
C +
         DO k = mzabso+1,mz
           rhodzk = (pstDY(i,j)*sigma(k)+ptopDY)
     .            / (ra*tairDY(i,j,k)*(1.+.608*qvDY(i,j,k)))
     .            * (gpmiDY(i,j,k)-gpmiDY(i,j,k+1))
C +...     rhodzk : (rho / 1000) * (dz * gravit)
C +
           optwa  = optwa + rhodzk * qwHY(i,j,k) 
           optia  = optia + rhodzk * qiHY(i,j,k)
c +
         END DO
C +
           WKxy6(i,j)  = 1.5 * ( optwa / 20.d-6 
     .                         + optia / 40.d-6 ) *grvinv
C +
      enddo
      enddo
C +
C +
C +--Vertically Integrated Water Vapor
C +  =================================
C +
      do j=1,my
      do i=1,mx
           WKxy7(i,j)  = 0.
         DO k = mzabso+1,mz
           WKxy7(i,j)  = WKxy7(i,j) + qvDY(i,j,k) *dsigm1(k)
         END DO
           WKxy7(i,j)  = WKxy7(i,j) * pstDYn(i,j) /gravit
      enddo
      enddo
C +
C +
C +   ************
      CALL UNwrite (ID__nc, 'RadOLR ', ipr2nc, mx, my, 1 , RAdOLR)
      CALL UNwrite (ID__nc, 'RSoTop ', ipr2nc, mx, my, 1 , WKxy5 )
      CALL UNwrite (ID__nc, 'OptDep ', ipr2nc, mx, my, 1 , WKxy6 )
      CALL UNwrite (ID__nc, 'WaterV ', ipr2nc, mx, my, 1 , WKxy7 )
C +   ************
C +
C +
C +--Work Arrays Reset
C +  -----------------
C +
      do j=1,my
      do i=1,mx
        WKxy1(i,j)   =0.0
        WKxy2(i,j)   =0.0
        WKxy3(i,j)   =0.0
        WKxy4(i,j)   =0.0
        WKxy5(i,j)   =0.0
        WKxy6(i,j)   =0.0
        WKxy7(i,j)   =0.0
        WKxy8(i,j)   =0.0
      enddo
      enddo
C +
      include "SBCnew.fTC"
C +
C +
C +--That 's all, folks: NetCDF File Closure
C +  =======================================
C +
C +   ************
      CALL UNclose (ID__nc)
C +   ************
C +
      return
      end
      subroutine SST_in
C +
C +------------------------------------------------------------------------+
C | MAR INPUT    SST                                        4-06-2002  MAR |
C |   SubRoutine SST_in initializes MAR         SSTs                       |
C |                     verifies    MARsst.CLI  EOF                        |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   OUTPUT: newlbc    : (0,1) ==> (NO new LBC ,new LBC)                  |
C |   ^^^^^^^                                                              |
C |                                                                        |
C |   OUTPUT (via common block)                                            |
C |   ^^^^^^  sst_LB       : Current                    SSTs               |
C |           sst1LB       : Previous Nesting Time Step SSTs               |
C |           sst2LB       : Next     Nesting Time Step SSTs               |
C |           tim1LB,tim2LB: LBC Nesting      Times  n, n+1                |
C |                                                                        |
C |   CAUTION: It is assumed that tim1LB and tim2LB do not change when the |
C |   ^^^^^^^^ Variables are reassigned after the dynamical Initialization |
C |            (Reassignation => itexpe:= nham => itimar:= itimar-nham*dt) |
C |                                                                        |
C |                                                                        |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
      IMPLICIT NONE
C +
C +
C +--Global Variables
C +  ================
C +
      include 'MARphy.inc'
C +
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
C +
      include 'MAR_LB.inc'
C +
      integer  newlbc
C +
      integer           isolCL(mx)
      common  /SST_in_i/isolCL
C +
      logical           SSTlin
      common  /SST_in_l/SSTlin
C +
      real              xxkmCL(mx),sh__CL(mx)
      common  /SST_in_r/xxkmCL,sh__CL
C +
C +
C +--Local  Variables
C +  ================
C +
      integer  itimar,iv_ilb
      real     rate
C +
C +
C +--Current Time
C +  ============
C +
          itimar=(((iyrrGE-1950)*365        + (iyrrGE-1949)/4
C +...    % Before Present(1950)               Nb of Leap Years % 1950
     .             +njyrGE(mmarGE)
     .             +njybGE(mmarGE)*max(0,1-mod(iyrrGE,4))
     .             +jdarGE-1                             )*24
     .            + jhurGE                                   )*3600
     .           +  minuGE                                         *  60
     .           +  jsecGE
C +
C +
C +--First Opening
C +  =============
C +
      IF (.NOT.SSTlin)                                            THEN ! CTR
CC#            SSTlin = .true.
       open (unit=32,status='old',form='unformatted',file='MARsst.CLI')! CTR
       rewind     32
         read    (32)         isolCL
         read    (32)         xxkmCL
         read    (32)         sh__CL
 101     CONTINUE
         read    (32,END=320) iyr_LB,mma_LB,jda_LB,jhu_LB,sst2LB
C +
             tim2LB=(((iyr_LB-1950)*365      +   (iyr_LB-1949)/4
     .                +njyrGE(mma_LB)
     .                +njybGE(mma_LB)*max(0,1-mod(iyr_LB,4))
     .                +jda_LB-1                             )*24
     .               + jhu_LB                                   )*3600
C +
         IF (tim2LB.gt.itimar)                               GO TO 100 ! CTR
             tim1LB      = tim2LB
         DO j=1,my
         DO i=1,mx
             sst1LB(i,j) = sst2LB(i,j)
         END DO
         END DO
             jdh_LB      =      1
         GO TO 101
 320     CONTINUE
             jdh_LB      =      0
 100     CONTINUE
C +
          write(6,6002)jda_LB,labmGE(mma_LB),iyr_LB,
     .                 jhu_LB,jdh_LB,               tim2LB
C +
       close(unit=32)
      END IF
C +
C +
C +--New LBC
C +  =======
C +
      IF (itimar.gt.   tim2LB)                                    THEN ! CTR
C +
          tim1LB =     tim2LB
C +
          write(6,6001)jda_LB,labmGE(mma_LB),iyr_LB,
     .                 jhu_LB,                      tim1LB,
     .                 jdarGE,labmGE(mmarGE),iyrrGE,
     .                 jhurGE,minuGE,        jsecGE,itimar
 6001     format(/, '  1st SST /',i3,'-',a3,'-',i4,i3,' ',2x,'/',2x,
     .              '   t =',i12,'s A.P.',
     .           /, '  Current /',i3,'-',a3,'-',i4,i3,':',i2,':',i2,
     .              '   t =',i12)
C +
       IF (jdh_LB.eq.0)jdh_LB = -1
       open (unit=32,status='old',form='unformatted',file='MARsst.CLI')! CTR
       rewind     32
         read    (32)         isolCL
         read    (32)         xxkmCL
         read    (32)         sh__CL
 11    CONTINUE
       IF (jdh_LB.le.0)                                       GO TO 10 ! CTR
C +
C +
C +--LBC at nesting time step n
C +  --------------------------
C +
         DO j=1,my
         DO i=1,mx
           sst1LB(i,j)      = sst2LB(i,j)     
         END DO
         END DO
C +
C +
C +--LBC at nesting time step n+1
C +  ----------------------------
C +
         read    (32,END=321) iyr_LB,mma_LB,jda_LB,jhu_LB,sst2LB
C +
          tim2LB=(((iyr_LB-1950)*365      +   (iyr_LB-1949)/4
     .             +njyrGE(mma_LB)
     .             +njybGE(mma_LB)*max(0,1-mod(iyr_LB,4))
     .             +jda_LB-1                             )*24
     .            + jhu_LB                                   )*3600
C +
       IF(itimar.gt.tim2LB)                                   GO TO 11 ! CTR
C +
          write(6,6002)jda_LB,labmGE(mma_LB),iyr_LB,
     .                 jhu_LB,jdh_LB,               tim2LB
 6002     format(   '  2nd LBC /',i3,'-',a3,'-',i4,i3,' ',2x,'/(',i1,
     .              ')  t =',i12)
C +
 10    CONTINUE
       GO TO 12
 321   CONTINUE
          jdh_LB=0
 12    CONTINUE
       close(unit=32)
C +
      ELSE                                                             ! CTR
c #WR     write(6,6003)jdarGE,labmGE(mmarGE),iyrrGE,
c #WR.                 jhurGE,minuGE,        jsecGE,itimar
 6003     format(   '  Current /',i3,'-',a3,'-',i4,i3,':',i2,':',i2,
     .              '   t =',i12,'s A.P.')
      END IF                                                           ! CTR
C +
C +
C +--Time Interpolation
C +  ==================
C +
      IF            (itimar.le.tim2LB  .and.   tim1LB.lt.tim2LB)  THEN ! CTR
C +
        rate = float(itimar  - tim1LB) / float(tim2LB  - tim1LB)
C +
          DO j=1,my
          DO i=1,mx
            sst_LB(i,j)     =sst1LB(i,j)      +
     .     (sst2LB(i,j)     -sst1LB(i,j)     )*rate
          END DO
          END DO
C +
        newlbc = 1
C +
      ELSE                                                             ! CTR
        newlbc = 0
      END IF                                                           ! CTR
C +
C +
C +--OUTPUT (1st Opening)
C +  ====================
C +
      IF (.NOT.SSTlin)                                            THEN ! CTR
               SSTlin = .true.
          write(6,6003)jdarGE,labmGE(mmarGE),iyrrGE,
     .                 jhurGE,minuGE,        jsecGE,itimar
          write(6,*) ' '
          write(6,*) '1st End:  SST_in'
      END IF
C +
      return
      end
