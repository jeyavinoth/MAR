
C +  ==========
C +--MAR_SL.inc
C +  ==========

      real           Ta2mSL(mx,my,mw),TminSL(mx,my,mw),TmaxSL(mx,my,mw)
      real           Ta3mSL(mx,my,mw),V03mSL(mx,my,mw),V10mSL(mx,my,mw)
      common /sheltr/Ta2mSL          ,TminSL          ,TmaxSL
     .              ,Ta3mSL          ,V03mSL          ,V10mSL
C +... Ta2mSL: interpolated 2 meter air temperature
C +    TminSL: minimum      2 meter air temperature 
C +    TmaxSL: maximum      2 meter air temperature

      real           tairSL(mx,my),dtagSL,tsrfSL(mx,my,mw)
      real           qvapSL(mx,my)
      common /srfdyn/tairSL       ,dtagSL,tsrfSL
     .              ,qvapSL
C +... tairSL: extrapolation of the sounding tempature to the surface.    (K)
C +    dtagSL: Air-Surface Temperature Difference                         (K)
C +    tsrfSL:     Surface Temperature                                    (K)
C +    qvapSL: specific humidity close to the surface                 (kg/kg)

      integer          isolSL(mx,my),maskSL(mx,my)
      common /srftyp/  isolSL       ,maskSL
C +...isolSL : Surface Type : 1 -> open ocean
C +                           2 -> glacier + ice sheet + snow
C +                           3 -> sea ice            (+ snow)
C +                           4 -> soil               (+ snow)
C +                           5 -> soil + vegetation
C +   maskSL : Land--Sea Mask 0 -> Continent
C +                           1 -> Ocean

      real             sol_SL(mx,my),firdSL(mx,my),firmSL(mx,my)
      real             hsenSL(mx,my),hlatSL(mx,my)
      real             fmelSL(mx,my),freeSL(mx,my),hbalSL(mx,my)
      common /srfbal/  sol_SL       ,firdSL       ,firmSL
     .              ,  hsenSL       ,hlatSL
     .              ,  fmelSL       ,freeSL       ,hbalSL
C +...                 sol_SL : Absorbed Solar    Heat Flux   (W/m2)
C +                    firdSL : Downward IR       Heat Flux   (W/m2)
C +                    firmSL : Upward   IR       Heat Flux   (W/m2)
C +                    hsenSL : Downward Sensible Heat Flux   (W/m2)
C +                    hlatSL : Evapora. Latent   Heat Flux   (W/m2)
C +                    fmelSL : Melting  Latent   Heat Flux   (W/m2)
C +                    freeSL : Freezing Latent   Heat Flux   (W/m2)
C +                    hbalSL : Total Surface     Heat Loss   (W/m2)

      real             cld_SL(mx,my),clduSL(mx,my),
     .                 cldmSL(mx,my),clddSL(mx,my)
      common /cldalb/  cld_SL,clduSL,cldmSL,clddSL
C +...                 cld_SL: Total Cloudiness above the Surface

      real             alb0SL  (mx,my),albeSL(mx,my),albsSL(mx,my)
     .                ,albxSL  (mx,my,mw)
     .                ,alb0SL_1(mx,my)
     .                ,alb0SL_2(mx,my)
     .                ,alb0SL_3(mx,my)  
      common /srfalb/  alb0SL         ,albeSL       ,albsSL
     .                ,albxSL
     .                ,alb0SL_1       ,alb0SL_2     ,alb0SL_3  
c +...alb0SL  : Background Surface Albedo
C +   albeSL  : Grid Mesh  Surface Albedo
C +   albxSL  : Mosaic     Surface Albedo
C +   albsSL  : Underlaying   Soil Albedo
C +   alb0SL_1: Background Surface Albedo for first  spectral band [0.3-0.8um]
C +   alb0SL_2: Background Surface Albedo for second spectral band [0.8-1.5um]
C +   alb0SL_3: Background Surface Albedo for third  spectral band [1.5-2.8um]          

      real             tgirSL(mx,my),eps0SL(mx,my)
      common /srf_IR/  tgirSL       ,eps0SL
C +...tgirSL : Effective Surface Radiative Temperature
C +   eps0SL : Surface IR Emissivity

      integer         nSLsrf(mx,my)
      common /srfrac/ nSLsrf
C +...nSLsrf : Number of  Sectors in a Grid Box                          (m/s)

      integer         ntavSL
      parameter      (ntavSL=20)
      real            SLsrfl(mx,my,mw)
      real            SLzoro(mx,my,mw)
      real             SL_z0(mx,my,mw),SLn_z0(mx,my,mw,ntavSL)
      real             SL_r0(mx,my,mw),SLn_b0(mx,my,mw,ntavSL)
      real                             SLn_r0(mx,my,mw,ntavSL)
      real             SLlmo(mx,my),SLlmol(mx,my,mw)
      real             SLuus(mx,my),SLuusl(mx,my,mw),duusSL(mx,my)
      real             SLuts(mx,my),SLutsl(mx,my,mw),SLdSdT(mx,my,mw)
      real            dutsSL(mx,my),cutsSL(mx,my,mw)
      real             SLuqs(mx,my),SLuqsl(mx,my,mw),SLdLdT(mx,my,mw)
      real            SaltSL(mx,my),SLussl(mx,my,mw),SLubsl(mx,my,mw)
      real             virSL(mx,my),fracSL
      common /srftur/ SLsrfl, SL_z0,SL_r0           ,SLzoro
     .              , SLn_z0       ,SLn_r0          ,SLn_b0
     .              ,  SLlmo       ,SLlmol
     .              ,  SLuus       ,SLuusl          ,duusSL
     .              ,  SLuts       ,SLutsl          ,SLdSdT
     .              , dutsSL       ,cutsSL
     .              ,  SLuqs       ,SLuqsl          ,SLdLdT
     .              , SaltSL       ,SLussl          ,SLubsl
     .              ,  virSL       ,fracSL
C +...SLsrfl : Normalized Sector  Area                                   (m/s)
C +   SLzoro : Roughness     Length for Momentum             (Orography) (m/s)
C +    SL_z0 : Roughness     Length for Momentum                         (m/s)
C +    SL_r0 : Roughness     Length for Heat                             (m/s)
C +    SLlmo : Monin-Obukhov Length                          (Average)   (m/s)
C +   SLlmol : Monin-Obukhov Length                                      (m/s)
C +    SLuus : Friction Velocity                             (Average)   (m/s)
C +   SLuusl : Friction Velocity                                         (m/s)
C +   duusSL : Correct. Friction Velocity                    (Average) (m2/s2)
C +    SLuts : Surface  Potential Temperature Turbulent Flux (Average) (K.m/s)
C +   SLutsl : Surface  Potential Temperature Turbulent Flux           (K.m/s)
C +   SLdSdT : Sensible Turbulent Heat Flux   T-Derivative            (W/m2/K)
C +   dutsSL : Correct. Potential Temperature Turbulent Flux (Average) (K.m/s)
C +   cutsSL : Correct. Potential Temperature Turbulent Flux           (K.m/s)
C +    SLuqs : Surface  Specific  Humidity    Turbulent Flux (Av.) (kg/kg.m/s)
C +   SLuqsl : Surface  Specific  Humidity    Turbulent Flux       (kg/kg.m/s)
C +   SLdLdT : Latent   Turbulent Heat Flux   T-Derivative            (W/m2/K)
C +   SaltSL : Friction Velocity  Saltation   Threshold                  (m/s)
C +   SLussl : Surface  Blowing   Snow        Turbulent Flux       (kg/kg.m/s)
C +   SLubsl : Surface  Blowing   Dust        Turbulent Flux       (kg/kg.m/s)
C +    virSL : Air Loading for SBL Parameterization                           
C +   fracSL : Fractional Time used in Blowing Snow Surface Flux Computation

C +--u*, u*T*, u*s* Time Moving Averages
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer                       ntaver
      parameter                    (ntaver=04) !09 in the old MAR version
      real          u_0aSL(mx,my,mw,ntaver)
      real          uT0aSL(mx,my,mw,ntaver)
      real          us0aSL(mx,my,mw,ntaver)
      common/srfmem/u_0aSL
     .             ,uT0aSL
     .             ,us0aSL

C +--V,  dT(a-s)    Time Moving Averages
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real          V_0aSL(mx,my,   ntaver)
      real          dT0aSL(mx,my,mw,ntaver)
      common/srfmem/V_0aSL
     .             ,dT0aSL

      real            cdnSL(mx,my,mw),cdmSL(mx,my,mw),cdhSL(mx,my,mw)
      real            zs_SL          ,zn_SL          ,zl_SL
      common /srf_z0/ cdnSL          ,cdmSL          ,cdhSL 
     .              , zs_SL          ,zn_SL          ,zl_SL
C +...cdnSL  : Square  Root of Momentum Neutral Drag Coefficient           (-)
C +   cdmSL  : Contribution of Momentum      to Drag Coefficient           (-)
C +   cdhSL  : Contribution of Heat          to Drag Coefficient           (-)
C +   zs_SL  : Typical Sea  Roughness Length                               (m)
C +   zn_SL  : Typical Snow Roughness Length                               (m)
C +   zl_SL  : Typical Land Roughness Length                               (m)

C +--u*, u*T*, u*s* Time Moving Averages
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real          aeCdSL(mx,my)          ! Aerodynamic Conductance
      real          Kv__SL(mx,my)          ! Equivalent  Turbulent Diffusivity
      real          pktaSL(mx,my)          ! SBC         Turbulence Pot. Temp.
      common/srfIMP/aeCdSL
     .             ,Kv__SL
     .             ,pktaSL

C +--Wind Speed
C +  ~~~~~~~~~~
      real            ssvSL(mx,my,mz)
      common /windno/ ssvSL
C +...ssvSL  : Wind Speed Norm                                           (m/s)

      real             ch0SL(mx,my),roseSL(mx,my)
      real            raerSL(mx,my),rsurSL(mx,my)
      common /srfhum/  ch0SL       ,roseSL,
     .                raerSL       ,rsurSL
C +... ch0SL : Bulk  Aerodynamic Coefficient Air/Surface Humidity Flux
C +   roseSL : Depletion of water vapor in the surface layer           (kg/kg)
C +            due to deposition of dew or rime on ground
C +   raerSL :`Bulk' Stomatal Resistance (Thom & Oliver, 1977, p. 347)
C +   rsurSL : Aerodynamic Resistance

      real            precSL(mx,my)
      common /precip/ precSL
C +...precSL : Rain 

      integer                      nLimit      !
      parameter                   (nLimit=10)  !
      real            WV__SL(mx,my,nLimit)     ! Water Vapor Flux Limitor
      common /WV_Lim/ WV__SL                   !

      real            bsnoSL(mx,my),snobSL(mx,my)
      real            hsnoSL(mx,my),snohSL(mx,my)
      common /snoacc/ bsnoSL       ,snobSL
     .              , hsnoSL       ,snohSL
C +...bsnoSL : cumulative snow       erosion height                        (m)
C +   snobSL :            snow       erosion height                        (m)
C +   hsnoSL : cumulative snow  accumulation height                        (m)
C +   snohSL :            snow precipitation height                        (m)

      real            hmelSL(mx,my)
      common /snomel/ hmelSL
C +...hmelSL : cumulative snowmelt           height       (m water equivalent)

      real            ro_SL0(mx,my),ro_SL(mx,my),d1_SL(mx,my),cs2SL
      real            t2SLtn,t2SLtd,t2_SL(mx,my),dtgSL(mx,my,mw)
      common /srfFrT/ ro_SL0       ,ro_SL       ,d1_SL       ,cs2SL
     .              , t2SLtn,t2SLtd,t2_SL       ,dtgSL
C +...ro_SL  : rhos               (Surface Density)                    (kg/m3)
C +   ro_SL0 : rhos       (Initial Surface Density)                    (kg/m3)
C +   d1_SL  : rhos * cs *(Depth diurnal Wave)                        (J/m2/K)
C +   cs2SL  : Soil Temperature Variation Time Scale                       (s)
C +           (usually 86400 s, i.e. diurnal cycle)
C +   t2SLtn,t2SLtd: time constants for t2_SL time integration
C +   t2_SL  : Soil Deep Layers Temperature                                (K)
C +   dtgSL  : Soil Temperature Variation during time interval dt          (K)
C +            Is renewed every 6 minutes

      logical         qsolSL
      common /srf_lo/ qsolSL

      real            hwatSL(mx,my)
      real             w20SL,w2_SL(mx,my),wg0SL,wg_SL(mx,my)
      real             wk0SL,wk_SL(mx,my),wx0SL,wx_SL(mx,my)
      common /srfFrH/ hwatSL
     .              ,  w20SL,w2_SL       ,wg0SL,wg_SL
     .              ,  wk0SL,wk_SL       ,wx0SL,wx_SL
C +...hwatSL ........ Rain Accumulation                                     (m)
C +   wg_SL and w2_SL Adimensional Numbers measuring Soil Water Content 
C +   wg_SL ......... near the     Surface
C +   w2_SL ......... over a large Soil Thickness 

      real            sst_SL
      common /srfOcn/ sst_SL
C +...sst_SL ........ Sea Surface Temperature

      real            SL_wge(mx,my),SLlwge(mx,my),SLuwge(mx,my)
      common /srfwge/ SL_wge       ,SLlwge       ,SLuwge
C +...SL_wge : wind gust estimate
C +   SLlwge : lower bound of the bounding interval around the estimate
C +   SLuwge : upper bound of the bounding interval around the estimate
