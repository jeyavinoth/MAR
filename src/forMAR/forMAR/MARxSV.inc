C +...BEGIN    "MARxSV.inc"


C +--SISVAT INPUT        Variables
C +  -----------------------------

      integer   LSmask(klonv)                 ! Land-Sea   Mask
      integer   isotSV(klonv)                 ! Soil       Type
      integer   iWaFSV(klonv)                 ! Soil       Drainage:(1,0)=(y,n)
      integer   ivgtSV(klonv)                 ! Vegetation Type

      real      coszSV(klonv)                 ! Cosine of Sun zenithal Angle
      real      sol_SV(klonv)                 ! Downward  Solar    Radiation
      real      IRd_SV(klonv)                 ! Downward  Longwave Radiation

      real      drr_SV(klonv)                 ! Rain  Intensity       [kg/m2/s]
      real      dsn_SV(klonv)                 ! Snow  Intensity       [kg/m2/s]
      real      dsnbSV(klonv)                 ! Idem, fraction, from Drift  [-]
      real      esnbSV(klonv)                 ! Idem, fraction, from Drift  [-]
      real      dbs_SV(klonv)                 ! Drift Amount            [kg/m2]
      real      BrosSV(klonv)                 ! Buffer Snow Layer Density
      real      BG1sSV(klonv)                 ! Buffer Snow Layer Dendr/Sphe[-]
      real      BG2sSV(klonv)                 ! Buffer Snow Layer Spher/Size[-] [0.0001 m]
      real      dz0_SV(klonv)                 ! dz0(Sastrugi dh)            [m]

      real      cld_SV(klonv)                 ! Cloudiness (seen from SBL)
      real      za__SV(klonv)                 ! SBL Height
      real      VV__SV(klonv)                 !(SBL Top)  Wind Velocity
      real      VV10SV(klonv)                 ! 10-m      Wind Velocity
      real      VVs_SV(klonv)                 !(Sastr,V)  Relevance
      real      RRs_SV(klonv)                 !(Sastr,V)  Counter
      real      DDs_SV(klonv)                 !(Sastr,V)  Angle
      real      TaT_SV(klonv)                 ! SBL Top   Temperature
      real      ExnrSV(klonv)                 ! Exner     Potential
      real      dSdTSV(klonv)                 ! Sensible Heat Flux T Derivat.
      real      dLdTSV(klonv)                 ! Latent   Heat Flux T Derivat.
      real      rhT_SV(klonv)                 ! SBL Top   Air  Density
      real      QaT_SV(klonv)                 ! SBL Top   Specific Humidity
      real      dQa_SV(klonv)                 ! SBL Flux  Limitation of Qa
      real      qsnoSV(klonv)                 ! SBL Mean  Snow       Content

      real      LAI0SV(klonv)                 ! Nominal Leaf Area Index
      real      glf0SV(klonv)                 ! Green   Leaf Fraction

      real      alb0SV(klonv)                 ! Soil    Albedo
      real      slopSV(klonv)                 ! Snow/Ice/Soil-Water Surf. Slope

      real      zSBLSV                        ! SBL Height (Initial Value)
      real      dt__SV                        ! Time Step 
      character*18     daHost                 ! Date Host Model


C +--SISVAT INPUT/OUTPUT Variables
C +  -----------------------------

      integer   isnoSV(klonv)                 ! Nb of Ice/Snow Layers
      integer   ispiSV(klonv)                 ! Uppermost superimposed ice
      integer   iiceSV(klonv)                 ! Nb of Ice      Layers
      integer   istoSV(klonv,    0:nsno)      ! Snow Layer     History

      real      alb_SV(klonv)                 ! Surface-Canopy Albedo
      real      alb1SV(klonv)                 ! Snow Albedo
      real      alb2SV(klonv)                 ! Snow Albedo
      real      alb3SV(klonv)                 ! Snow Albedo
      real      emi_SV(klonv)                 ! Surface-Canopy Emissivity
      real      IRs_SV(klonv)                 ! Soil           IR Flux
      real      LMO_SV(klonv)                 ! Monin-Obukhov  Scale
      real      us__SV(klonv)                 ! Friction       Velocity
      real      uts_SV(klonv)                 ! Temperature Turbulent Scale
      real      cutsSV(klonv)                 ! Temperature Turbulent Scale C.
      real      uqs_SV(klonv)                 ! Spec.Humid. Turbulent Scale
      real      uss_SV(klonv)                 ! Blow.Snow   Turbulent Scale
      real      usthSV(klonv)                 ! Blowing Snow Erosion Thresh.
      real      rCDmSV(klonv)                 ! Square  Root Contribut. Drag_m
      real      rCDhSV(klonv)                 ! Square  Root Contribut. Drag_h
      real      Z0m_SV(klonv)                 ! Momentum     Roughness Length
      real      Z0mmSV(klonv)                 !  z0(Momentum,    Time Mean) [m]
      real      Z0mnSV(klonv)                 !  z0(Momentum,    instanta.) [m]
      real      Z0roSV(klonv)                 ! Subgrid Topo Roughness Length
      real      Z0SaSV(klonv)                 !  z0(Sastrugi  h)            [m]
      real      Z0e_SV(klonv)                 !  z0(Snow eroded)            [m]
      real      Z0emSV(klonv)                 !  z0(Snow eroded, Time Mean) [m]
      real      Z0enSV(klonv)                 !  z0(Snow eroded, instanta.) [m]
      real      Z0h_SV(klonv)                 ! Heat         Roughness Length
      real      Z0hmSV(klonv)                 !  z0(Heat,        Time Mean) [m]
      real      Z0hnSV(klonv)                 !  z0(Heat,        instanta.) [m]

      real      snCaSV(klonv)                 ! Canopy  Snow   Thickness
      real      rrCaSV(klonv)                 ! Canopy  Water  Content
      real      psivSV(klonv)                 ! Leaf    Water  Potential 
      real      TvegSV(klonv)                 ! Vegetation     Temperature

      real      TsisSV(klonv,-nsol:nsno)      ! Snow/Ice/Soil-Water Temperature
      real      ro__SV(klonv,-nsol:nsno)      ! Snow/Ice/Soil-Water VolumicMass
      real      eta_SV(klonv,-nsol:nsno)      ! Snow/Ice/Soil     Water Content
      real      G1snSV(klonv,    0:nsno)      ! Snow Dendricity/Sphericity
      real      G2snSV(klonv,    0:nsno)      ! Snow Sphericity/Size
      real      dzsnSV(klonv,    0:nsno)      ! Snow Layer  Thickness
      real      agsnSV(klonv,    0:nsno)      ! Snow Age
      real      BufsSV(klonv)                 ! Snow Buffer Layer
      real      rusnSV(klonv)                 ! Surficial   Water
      real      SWf_SV(klonv)                 ! Normalized  Decay
      real      SWS_SV(klonv)                 ! Surficial Water Status
      real      HFraSV(klonv)                 ! Frazil      Thickness

      real      zWE_SV(klonv)                 ! Current   Snow Thickness [mmWE]
      real      zWEcSV(klonv)                 ! Compacted Snow Thickness [mmWE]
      real      wem_SV(klonv)                 ! Only Melting             [mmWE]
      real      wer_SV(klonv)                 ! Refreezing               [mmWE]
      real      wes_SV(klonv)                 ! Sublimation              [mmWE]


C +--SISVAT OUTPUT       Variables
C +  -----------------------------

      integer   no__SV(nb_wri)                ! OUTPUT file Unit Number
      integer   i___SV(nb_wri)                ! OUTPUT point   i Coordinate
      integer   j___SV(nb_wri)                ! OUTPUT point   j Coordinate
      integer   n___SV(nb_wri)                ! OUTPUT point   n Coordinate
      integer   lwriSV(klonv)                 ! OUTPUT point vec Index

      integer   ii__SV(klonv),iwr_SV          ! WORK   point   i Coordinate
      integer   jj__SV(klonv),jwr_SV          ! WORK   point   j Coordinate
      integer   nn__SV(klonv),nwr_SV          ! WORK   point   n Coordinate

      real      IRu_SV(klonv)                 ! UPward    IR Flux (effective)
      real      hSalSV(klonv)                 ! Saltating Layer Height
      real      qSalSV(klonv)                 ! Saltating Snow  Concentration
      real      RnofSV(klonv)                 ! RunOFF    Intensity

      common/xSISVAT_C/daHost
      common/xSISVAT_I/LSmask,ivgtSV,isotSV,iWaFSV
     .                ,isnoSV,ispiSV,iiceSV,istoSV
     .                ,no__SV,i___SV,j___SV,n___SV,lwriSV
     .                ,ii__SV,jj__SV,nn__SV
     .                ,iwr_SV,jwr_SV,nwr_SV
      common/xSISVAT_R/coszSV,sol_SV,IRd_SV
     .                ,drr_SV,dsn_SV,dsnbSV,esnbSV,dbs_SV
     .                ,BrosSV,BG1sSV,BG2sSV,dz0_SV
     .                ,cld_SV,za__SV,VV__SV,VV10SV,TaT_SV,ExnrSV
     .                ,VVs_SV,RRs_SV,DDs_SV
     .                ,dSdTSV,dLdTSV,rhT_SV,QaT_SV,qsnoSV
     .                ,LAI0SV,glf0SV,alb0SV,slopSV,zSBLSV,dt__SV
     .                ,alb_SV,emi_SV,IRs_SV
     .                ,LMO_SV,us__SV,uts_SV,cutsSV,uqs_SV,uss_SV,usthSV
     .                ,rCDmSV,rCDhSV
     .                ,Z0m_SV,Z0mnSV,Z0mmSV,Z0roSV,Z0SaSV
     .                ,Z0e_SV,Z0enSV,Z0emSV
     .                ,Z0h_SV,Z0hnSV,Z0hmSV,snCaSV,rrCaSV
     .                ,psivSV,TvegSV,TsisSV,ro__SV,eta_SV
     .                ,G1snSV,G2snSV,dzsnSV,agsnSV,BufsSV,rusnSV,SWf_SV
     .                ,SWS_SV,HFraSV,IRu_SV,hSalSV,qSalSV,RnofSV
     .                ,zWE_SV,zWEcSV,wem_SV,wer_SV,wes_SV
     .                ,alb1SV,alb2SV,alb3SV

C +...END      "MARxSV.inc"
