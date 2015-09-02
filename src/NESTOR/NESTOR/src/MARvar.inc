

C +---SPECIFIC VARIABLES FOR MAR
C +---(Modele Atmospherique Regional)
C +---===============================

C +---Parameters
C +   ----------

      INTEGER mzz
      PARAMETER (mzz=mz+1)


C +---Variables of MARgrd.inc
C +   -----------------------

      REAL    dx,dy,dt,dtquil,tequil
C +...        dx    : horizontal grid size (x direction)
C +...        dy    : horizontal grid size (y direction)
C +...        dt    : time step
C +...        dtquil: Time Step for 1-D Initialisation (s)
C +...        tequil: Time Span for 1-D Initialisation (s)

      INTEGER imez,jmez
C +...        imez  : reference grid point (x direction)
C +...        jmez  : reference grid point (y direction)

      REAL    sigma(mz),sh(mx,my)
C +...        sigma : independant variable sigma on sigma levels (k)
C +...        sh    : surface height (m)

c #CO COMMON / MARgrd / dx,dy,dt,dtquil,tequil,imez,jmez,sigma,sh


C +---Variables of MAR_GE.inc
C +   -----------------------

      REAL    GElat0,GElon0,GEddxx
C +...        GElat0: Latitude  (Degrees)
C +...        GElon0: Longitude (Degrees)
C +...        GEddxx: Direction x-axis

c #CO COMMON / MAR_GE1 / GElat0,GElon0
      COMMON / MAR_GE2 / GEddxx

C +---Variables of MAR_DY.inc
C +   -----------------------

      INTEGER iyrDYN,mmaDYN,jdaDYN,jhuDYN
C +...        iyrDYN: Year
C +...        mmaDYN: Month
C +...        jdaDYN: Day
C +...        jhuDYN: Hour  (UT)

      REAL    pstDY (mx,my)    ,ptopDY          ,
     .        fcorDY(mx,my)
c ###.       ,uairDY(mx,my,mz) ,vairDY(mx,my,mz),
c ###.       ,pktaDY(mx,my,mzz),tairDY(mx,my,mz),
c ###.       ,qvDY  (mx,my,mz)
C +...        uairDY: x-wind speed component                       (m/s)
C +...        vairDY: y-wind speed component                       (m/s)
C +...        pktaDY: potential temperature divided by 100.[kPa]**(R/Cp)
C +...        tairDY: real      temperature                          (K)
C +...        qvDY  :            Specific Humidity               (kg/kg)
C +...        qsatDY: Saturation Specific Humidity               (kg/kg)
C +...        pstDY1: Model Pressure Depth at INITIAL  Time Step   (kPa)
C +...        pstDY : Model Pressure Depth at current  Time Step   (kPa)
C +...        ptopDY:       Pressure       at Model Top            (kPa)
C +...        fcorDY: Coriolis Parameter                           (s-1)

c #CO COMMON / MAR_DY / iyrDYN,mmaDYN,jdaDYN,jhuDYN,uairDY,vairDY,
c #CO.                  pktaDY,tairDY,pkDY  ,qvDY  ,pstDY ,ptopDY,
c #CO.                  fcorDY


C +---Variables of MAR_SL.inc
C +   -----------------------

      INTEGER isolSL(mx,my),maskSL(mx,my),nSLsrf(mx,my)
C +...        isolSL : Surface Type : 1 -> open ocean
C +...                                2 -> glacier + ice sheet + snow
C +...                                3 -> sea ice            (+ snow)
C +...                                4 -> soil               (+ snow)
C +...                                5 -> soil + vegetation
C +...        maskSL : Land--Sea Mask 0 -> Continent
C +...                                1 -> Ocean
C +...        nSLsrf : Number of  Sectors in a Grid Box

      REAL    zs_SL,zn_SL,zl_SL,cs2SL,sst_SL,dtagSL,
     .        wk0SL,wx0SL,w20SL,wg0SL
C +...        zs_SL : Typical Sea  Roughness Length (m)
C +...        zn_SL : Typical Snow Roughness Length (m)
C +...        zl_SL : Typical Land Roughness Length (m)
C +...        cs2SL : Soil Temperature Variation Time Scale (s)
C +...                 (usually 86400 s, i.e. diurnal cycle)
C +...        sst_SL: Sea Surface Temperature
C +...        dtagSL: Air-Surface Temperature Difference    (K)
C +...        w**SL : Initial soil humidity variables

      REAL    tairSL(mx,my)   ,
     .        tsrfSL(mx,my,mw),qvapSL(mx,my)   ,
     .        alb0SL(mx,my)   ,albeSL(mx,my)   ,
     .        albsSL(mx,my)   ,eps0SL(mx,my)   ,
     .        SLsrfl(mx,my,mw),SL_z0(mx,my,mw) ,
     .        SL_r0 (mx,my,mw),SLlmo(mx,my)    ,
     .        SLlmol(mx,my,mw),SLuus(mx,my)    ,
     .        SLuusl(mx,my,mw),SaltSL(mx,my)   ,
     .        virSL (mx,my)   ,fracSL          ,
     .        SLuts (mx,my)   ,SLutsl(mx,my,mw),
     .        SLuqs (mx,my)   ,SLuqsl(mx,my,mw),
     .        ch0SL (mx,my)   ,roseSL(mx,my)   ,
     .        raerSL(mx,my)   ,rsurSL(mx,my)   ,
     .        hmelSL(mx,my)   ,ro_SL0(mx,my)   ,
     .        ro_SL (mx,my)   ,d1_SL (mx,my)   ,
     .        t2_SL (mx,my)   ,w2_SL (mx,my)   ,
     .        wg_SL (mx,my)   ,hsnoSL(mx,my)
C +...        tairSL : Extrapolation of the sounding tempature to the surface.(K)
C +...        tsrfSL : Surface Temperature                                    (K)
C +...        qvapSL : specific humidity close to the surface             (kg/kg)
C +...        alb0SL : Background Surface Albedo
C +...        albeSL :            Surface Albedo
C +...        albsSL : Underlaying   Soil Albedo
C +...        eps0SL : Surface IR Emissivity
C +...        SLsrfl : Normalized Sector  Area                         (m/s)
C +...        SL_z0  : Roughness     Length for Momentum               (m/s)
C +...        SL_r0  : Roughness     Length for Heat                   (m/s)
C +...        SLlmo  : Monin-Obukhov Length                (Average)   (m/s)
C +...        SLlmol : Monin-Obukhov Length                            (m/s)
C +...        SLuus  : Friction Velocity                   (Average)   (m/s)
C +...        SLuusl : Friction Velocity                               (m/s)
C +...        SaltSL : Friction Velocity Saltation   Threshold         (m/s)
C +...        virSL  : Air Loading for SBL Parameterization                           
C +...        fracSL : Fractional Time used in Blowing Snow Surface Flux Computation
C +...        SLuts  : Surface Potential Temperature Turbulent Flux  (Average) (K.m/s)
C +...        SLutsl : Surface Potential Temperature Turbulent Flux            (K.m/s)
C +...        SLuqs  : Surface Specific Humidity     Turbulent Flux  (Av.) (kg/kg.m/s)
C +...        SLuqsl : Surface Specific Humidity     Turbulent Flux        (kg/kg.m/s)
C +...        ch0SL  : Bulk  Aerodynamic Coefficient Air/Surface Humidity Flux
C +...        roseSL : Depletion of water vapor in the surface layer           (kg/kg)
C +...                 due to deposition of dew or rime on ground
C +...        raerSL :`Bulk' Stomatal Resistance (Thom & Oliver, 1977, p. 347)
C +...        rsurSL : Aerodynamic Resistance
C +...        hmelSL : cumulative snowmelt height (m water equivalent)
C +...        hsnoSL : cumulative snow accumulation height (m water equivalent)
C +...        ro_SL  : rhos               (Surface Density)                    (kg/m3)
C +...        ro_SL0 : rhos       (Initial Surface Density)                    (kg/m3)
C +...        d1_SL  : rhos * cs *(Depth diurnal Wave)                        (J/m2/K)
C +...        t2_SL  : Soil Deep Layers Temperature                                (K)
C +...        dtgSL  : Soil Temperature Variation during time interval dt          (K)
C +...                 Is renewed every 6 minutes
C +...        wg_SL and w2_SL Adimensional Numbers measuring Soil Water Content 
C +...        wg_SL  : ... near the     Surface
C +...        w2_SL  : ... over a large Soil Thickness 

c #CO COMMON / MAR_SL / isolSL,maskSL,nSLsrf,zs_SL,zn_SL,zl_SL,
c #CO.                  cs2SL,sst_SL,dtagSL,wk0SL,wx0SL,w20SL,
c #CO.                  wg0SL,tairSL,tsrfSL,qvapSL,alb0SL,albeSL,
c #CO.                  albsSL,eps0SL,SLsrfl,SL_z0,SL_r0,SLlmo,
c #CO.                  SLlmol,SLuus,SLuusl,SaltSL,virSL,fracSL,
c #CO.                  SLuts,SLutsl,SLuqs,SLuqsl,ch0SL,roseSL,
c #CO.                  raerSL,rsurSL,hmelSL,ro_SL0,ro_SL,d1_SL,
c #CO.                  t2_SL,w2_SL,wg_SL


C +---Variables of MAR_TV.inc
C +   -----------------------

      INTEGER isolTV(mx,my)
C +...        isolTV : Texture Type :  2 -> loamy sand
C +...                                 5 -> loam
C +...                                11 -> clay


C +---Variables of MAR_TU.inc
C +   -----------------------

      REAL    TUkhff,TUkhmx
C +...        TUkhff: Horiz.vKar**2 (horizontal diffusion)
C +...        TUkhmx: Upper sponge

c #CO COMMON / MAR_TU / TUkhff,TUkhmx


C +---Variables of MAR_LB.inc
C +   -----------------------

      INTEGER iyr_LB,mma_LB,jda_LB,jhu_LB,jdh_LB
C +...        iyr_LB: Year
C +...        mma_LB: Month
C +...        jda_LB: Day
C +...        jhu_LB: Hour  (UT)
C +...        jdh_LB: Time Interval before next GCM/NWP LBC (hour)
C +...        jdh_LB=0 ==> NO further GCM/NWP LBC available

      REAL    vaxgLB(1:n7,my,mz,5),vaxdLB(mx-n6:mx ,my,mz,5),
     .        vayiLB(mx,1:n7,mz,5),vaysLB(mx,my-n6:my ,mz,5)
C +...        vaXX : large scale values of relevant dependant variables 
C +...          ^X=(x->x axis border, y->y axis border)
C +...           ^X=(g->x small, d->x large, b->y small, h->y large)

c #CO COMMON / MAR_LB / iyr_LB,mma_LB,jda_LB,jhu_LB,jdh_LB,
c #CO.                  vaxgLB,vaxdLB,vayiLB,vaysLB


C +---Variables of MAR_PO.inc
C +   -----------------------

      REAL    uocnPO(mx,my),vocnPO(mx,my),aPOlyn(mx,my)
C +...        uocnPO: Oceanic Current (prescribed, x-direction)
C +...        vocnPO: Oceanic Current (prescribed, y-direction)
C +...        aPOlyn: Initial (observed) lead fraction

c #CO COMMON / MAR_PO / uocnPO,vocnPO,aPOlyn


C +---Variables of MAR_IO.inc
C +   -----------------------

      INTEGER igrdIO(5),jgrdIO(5)
C +...        igrdIO: i (x-direc.) Index Ref. Grid Point (for detailed Output)
C +...        jgrdIO: j (y-direc.) Index Ref. Grid Point (for detailed Output)

      CHARACTER*3 explIO

c #CO COMMON / MAR_IO / igrdIO,jgrdIO,explIO


C +---Variables of MAR_FI.inc
C +   -----------------------

      REAL    FIslot,FIslou,FIslop,FIkhmn
C +...        FIslot: Implicit Filter Parameter (Temperature)
C +...        FIslou: Implicit Filter Parameter (Wind Speed)
C +...        FIslop: Implicit Filter Parameter (Pressure)
C +...        FIkhmn: Horizontal Diffusion Coefficient

c #CO COMMON / MAR_FI / FIslot,FIslou,FIslop,FIkhmn


C +---Other variables
C +   ---------------

      INTEGER   maptyp
      INTEGER*8 itexpe
C +...          maptyp: Projection Type
C +...          itexpe: Iteration number

      REAL    zmin,aavu,bbvu,ccvu
C +...        zmin  : height of the first model level from the surface (m)
C +...        aavu  : grid parameter (geometric progression)
C +...        bbvu  : grid parameter (geometric progression)
C +...        ccvu  : grid parameter (geometric progression)

      LOGICAL vertic

c #CO COMMON / MARoth / itexpe,zmin,aavu,bbvu,ccvu,vertic
      COMMON / MARmap / maptyp

