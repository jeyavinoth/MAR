C +
C +  ==========
C +--MAR_TV.inc
C +  ==========
C +
      INTEGER    nvx   ,llx       ,iptx                                               
      PARAMETER (nvx=mw,llx=nsol+1,iptx=5)                                            
C +                                                                             
      integer    imx   ,jmx                                                     
      PARAMETER (imx=mx,jmx=my)                                                 
C +                                                                             
      real      deptTV(0:llx)                                                   
C +...          deptTV:     Soil Level Depth                                    
C +                                                                             
      real      dep2TV(0:llx)                                                   
C +...          dep2TV:     Soil Layer Depth                                    
C +                                                                             
      real      slopTV(imx,jmx)
C +...          slopTV:     Surface    Slope
C +
      real      AlbSTV(imx,jmx)                                                 
C +...          AlbSTV: Dry Soil       Albedo                                   
C +                                                                             
      real      alaiTV(imx,jmx,nvx)                                             
C +...          alaiTV:           Leaf Area Index                               
C +                                                                             
      real      glf_TV(imx,jmx,nvx)                                             
C +...          glf_TV:     Green Leaf Fraction                                 
C +                                                                             
      real      CaWaTV(imx,jmx,nvx)                                             
C +...          CaWaTV: Canopy     Intercepted Water Content                    
C +                                                                             
      real      CaSnTV(imx,jmx,nvx)                                             
C +...          CaSnTV: Canopy     Intercepted Snow  Content                    
C +                                                                             
      real      TvegTV(imx,jmx,nvx)                                             
C +...          TvegTV: Skin  Vegetation Temperature                            
C +                                                                             
      real      TgrdTV(imx,jmx,nvx)                                             
C +...          TgrdTV: Skin  Soil       Temperature                            
C +                                                                             
      real      TsolTV(imx,jmx,nvx,llx)                                         
C +...          TsolTV: Layer Soil       Temperature                            
C +                                                                             
      real      eta_TV(imx,jmx,nvx,llx)                                         
C +             eta_TV: Soil                Moisture Content                    
C +                                                                             
      real      psigTV(imx,jmx,nvx)                                             
C +             psigTV: Soil               Hydraulic Potential                  
C +                                                                             
      real      psivTV(imx,jmx,nvx)                                             
C +             psivTV: Vegetation         Hydraulic Potential                  
C +                                                                             
      real      evapTV(imx,jmx)                                                 
C +             evapTV: Time Integrated Evapotranspiration                      
C +                                                                             
      real      runoTV(imx,jmx)                                                 
C +             runoTV: Time Integrated (Sub)surface Flow                       
C +                                                                             
      real      draiTV(imx,jmx)                                                 
C +             draiTV: Time Integrated Drainage     Flow                       
C +                                                                             
      INTEGER iWaFTV(imx,jmx)                                                   
C +...       (iWaFTV=0 ==> no Water Flux;                                       
C +           iWaFTV=1 ==> free drainage)                                       
C +                                                                             
      INTEGER ivegTV(imx,jmx,nvx)                                               
C +...        ivegTV:     Vegetation Type Index                                 
C +                                                                             
      INTEGER isolTV(imx,jmx)                                                   
C +...        isolTV:     Soil       Type Index                                 
C +                                                                             
      real    ifraTV(imx,jmx,nvx)                                               
C +...        ifraTV:     Vegetation Class Coverage                             
C +                               (3 Class, Last One is Open Water)             
C +                                                                             
      INTEGER IOi_TV(iptx),IOj_TV(iptx)                                         
C +...        IO Grid                     Indices                               
C +                                                                             
      integer       itx   ,ivg                                                  
C +                                                                             
      common/rsvaTV/AlbSTV,alaiTV,glf_TV,CaWaTV,CaSnTV,runoTV,draiTV,
     .              evapTV,TvegTV,TgrdTV,TsolTV,eta_TV,psigTV,psivTV,                  
     .              slopTV,deptTV,dep2TV,ifraTV                                               
      common/isvaTV/iWaFTV,ivegTV,isolTV,IOi_TV,IOj_TV,
     .              itx   ,ivg                                                  
