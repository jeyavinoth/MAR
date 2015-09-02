
C +  ==========
C +--MAR_RA.inc
C +  ==========

      real           pktRAd(mx,my,mz),RAdsol(mx,my),
     .               RAd_ir(mx,my)   ,RAdOLR(mx,my), RadOSR(mx,my),
     .                tviRA(mx,my)   ,IRsoil(mx,my,mw)
      real           RAcd_O(mx,my,mz),RAcdtO(mx,my) 
      real           RAer_O(mx,my,mz),RAertO(mx,my)

      common/radmar/ pktRAd          ,RAdsol 
     .             , RAd_ir          ,RAdOLR, RadOSR 
     .             ,  tviRA          ,IRsoil          
     .             , RAcd_O          ,RAcdtO          
     .             , RAer_O          ,RAertO          
C +...               pktRAd: Radiation Divergence ==> Diabatic Heating
C +                  RAdsol: Surface Downward  Solar Heat Flux
C +                  RAd_ir: Surface Downward  IR    Heat Flux
C +                  RAdOLR: Outgoing  Longwave Radiation
C +                  RAdOSR: Outgoing Shortwave Radiation
C +                   tviRA: Radiating Surface Temperature    
C +                  IRsoil: Soil     upward   IR    
C +                  RAcd_O, RAcdtO: Cloud Optical Depth [-]
C +                  RAer_O, RAertO: Aeros.Optical Depth [-]

      real          pmb_RA(mzir),tairRA(mzir),qvapRA(mzir),
     .              cld_RA(mzir),qlwpRA(mzir),toptRA(mzir),
     .              wo3_RA(mzir),dwo3RA(mzir),htngRA(mzir)
      common/radiat/pmb_RA      ,tairRA      ,qvapRA      ,
     .              cld_RA      ,qlwpRA      ,toptRA      ,
     .              wo3_RA      ,dwo3RA      ,htngRA      
C +...              pmb_RA: Pressure                                   (hPa)
C +                 tairRA: Temperature                                  (K)
C +                 qvapRA: Specific Humidity                        (kg/kg)
C +                 cld_RA: Cloudiness              at level k
C +                 qlwpRA: Cloud Liquid Water Path   (level k)            
C +                 toptRA: Cloud Optical Depth       (level k)            
C +                 wo3_RA: Ozone Content                           (cm STP)
C +                 dwo3RA: Ozone Content          (Divergence)                
C +                 htngRA: Radiative Heating                          (K/s)
