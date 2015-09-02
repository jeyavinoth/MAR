      real           dtPO,  silfPO,sw00PO,si00PO,aPOlyn
      common/polyda/ dtPO,  silfPO,
     .               sw00PO,si00PO,
     .               aPOlyn(mx,my)
C +...                 dtPO        : Polynya Model Time Step
C +                  silfPO        : Sea Ice Solidification: Latent Heat
C +                  sw00PO, si00PO: water, sea-ice salinities (TN Bay)
C +                           (Bromwich and Kurtz,   1984, JGR, p.3568;
C +                            Cavalieri and Martin, 1985,      p. 248)
C +...               aPOlyn        : Initial (observed) Lead Fraction
C +
      integer        iPO1,iPO2,jPO1,jPO2,
     .               iPO3,iPO4,jPO3,jPO4
      common/polyni/ iPO1,iPO2,jPO1,jPO2,
     .               iPO3,iPO4,jPO3,jPO4
C +...               iPO1,iPO2,iPO3,iPO4 define output area (x-axis)
C +                  jPO1,jPO2,jPO3,jPO4 define output area (y-axis)
C +
      real           uocnPO       ,vocnPO       ,focnPO,
     .               albwPO       ,hatmPO       ,
     .               swsaPO       ,tsswPO       ,cwswPO,
     .               awmxPO
      common/polocn/ uocnPO(mx,my),vocnPO(mx,my),focnPO(mx,my),
     .               albwPO(mx,my),hatmPO(mx,my),
     .               swsaPO(mx,my),tsswPO(mx,my),cwswPO(mx,my),
     .               awmxPO(mx,my)
C +...               uocnPO: Oceanic Current   (prescribed,    x-direction)
C +                  vocnPO: Oceanic Current   (prescribed,    y-direction)
C +                  focnPO: Oceanic Heat Flux (in the oceanic mixed layer)
C +                  albwPO: Open Water Surface Albedo
C +                  hatmPO: Heat Losses from Open Water Area
C +                  swsaPO: Oceanic Salinity  (in the oceanic mixed layer)
C +                  tsswPO: Sea Water Freezing Temperature 
C +                  cwswPO: Sea Water Heat Capacity     
C +                  awmxPO: Max Fraction of Open Water
C +
      real           uicePO       ,vicePO,
     .               ufraPO       ,vfraPO       ,hfraPO,
     .               hicePO       ,aicePO,
     .               hiavPO       ,vgriPO
      common/police/ uicePO(mx,my),vicePO(mx,my),
     .               ufraPO(mx,my),vfraPO(mx,my),hfraPO(mx,my),
     .               hicePO(mx,my),aicePO(mx,my),
     .               hiavPO(mx,my),vgriPO(mx,my)
C +                  uicePO, vicePO: Sea-Ice Velocity
C +                  ufraPO, vfraPO: Fraizil Velocity
C +                  hfraPO        : Fraizil Thickness
C +                  hicePO        : Sea-Ice Thickness (Old   Ice)
C +                  aicePO        : Sea-Ice Fraction  (Old   Ice)
C +                  hiavPO        : Sea-Ice Thickness (Average)   
C +                  vgriPO        : Grease Ice Volume (New   Ice)
C +
