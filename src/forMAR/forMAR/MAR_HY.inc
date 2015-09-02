
C +  ==========
C +--MAR_HY.inc
C +  ==========

      logical           turnHY
      common /hydlog/   turnHY
C +...    turnHY: set = T when rhcrHY is reached => microphysics turned on

      real              cminHY
      common /hy_dat/   cminHY
C +...cminHY : cloud fraction (minimum)

      real              rhcrHY,tim_HY
      common /hydctr/   rhcrHY,tim_HY          
C +       rhcrHY: critical relative humidity value for condensation 
C +       tim_HY: time at which microphysics are taken into account

      real           cfraHY(mx,my,mz)
      common /hydfra/cfraHY
C +...cfraHY : cloud fraction

      real           ccnwHY(mx,my,mz),ccniHY(mx,my,mz)
      common /hydccn/ccnwHY          ,ccniHY
C +...ccnwHY : cloud dropplets    number concentration   (Nb/m3)
C +   ccniHY : cloud ice crystals number concentration   (Nb/m3)

      real             qwHY(mx,my,mz),  qiHY(mx,my,mz )   
     .                ,qrHY(mx,my,mz),  qsHY(mx,my,mzz)
C +HG.                ,qgHY(mx,my,mz)
      common /hydmet/  qwHY          ,  qiHY
     .                ,qrHY          ,  qsHY
C +HG.                ,qgHY
C +...qwHY   : cloud dropplets    concentration          (kg/kg)
C +   qiHY   : cloud ice crystals concentration          (kg/kg)
C +   qrHY   : rain               concentration          (kg/kg)
C +   qsHY   : snow flakes        concentration          (kg/kg)
C +   qgHY   : graupels           concentration          (kg/kg)
C +        (translucent snow pellets encased in a thin layer of ice)

      real           rainHY(mx,my),snowHY(mx,my),crysHY(mx,my)
      real           rai0HY(mx,my),sno0HY(mx,my),sfa0HY(mx,my)
      real           qsrfHY(mx,my),uss_HY(mx,my)
      common /hydsrf/rainHY       ,snowHY       ,crysHY       ,
     .               rai0HY       ,sno0HY       ,sfa0HY       ,
     .               qsrfHY       ,uss_HY
C +...rainHY : integrated precipited        rain 
C +   rai0HY : integrated precipited        rain  (previous time step)
C +   snowHY : integrated precipited/eroded snow 
C +   sno0HY : integrated precipited/eroded snow  (previous time step)
C +   sfa0HY : integrated precipited        snow  (previous time step)
C +   qsrfHY : Blowing Snow Concentration  (0.325 m above the surface)
C +   uss_HY : Blowing Snow Turbulent Flux (0.325 m above the surface)

      real           hlatHY(mx,my,mz), dqiHY(mx,my,mz),dqwHY(mx,my,mz)
      real           Hcd_HY(mx,my)   ,Tcd_HY(mx,my)  ,zcd_HY(mx,my)
      real           Hsb_HY(mx,my)   ,Tsb_HY(mx,my)  ,zsb_HY(mx,my)
      common /hyddyn/hlatHY          , dqiHY          ,dqwHY
     .              ,Hcd_HY          ,Tcd_HY         ,zcd_HY
     .              ,Hsb_HY          ,Tsb_HY         ,zsb_HY
C +...hlatHY : latent heat release [K/s] &  BLOWN * SUBLIMATION [m w.e.]
C +                                         (for k=1)
C +    dqiHY :  Sublimation(-) /   Deposition(+)                [m w.e.]
C +    dqwHY : Vaporisation(-) / Condensation(+)                [m w.e.]
C +   Hcd_HY : latent heat release                             [mm w.e.]
C +   Tcd_HY : latent heat release weighted Air Temperature          [K]
C +   Hsb_HY : latent heat absorb.                             [mm w.e.]
C +   Tsb_HY : latent heat absorb. weighted Air Temperature          [K]

      integer        icntHY
      common /hydint/icntHY
