      integer   NLaysv(klonv)                 ! New   Snow     Layer   Switch
      integer   i_thin(klonv)                 ! Index of the thinest Layer
      integer   LIndsv(klonv)                 ! Contiguous Layer relative Index
      common  /ySISVAT_I/
     .          NLaysv,i_thin,LIndsv
C +
      real      albisv(klonv)                 ! Integrated Surface Albedo
      real      albssv(klonv)                 ! Soil               Albedo [-]
      real      SoCasv(klonv)                 ! Canopy  Absorbed Solar Radiat.
      real      SoSosv(klonv)                 ! Surface Absorbed Solar Radiat.
      real      IRv_sv(klonv)                 ! Vegetation IR Flux  [W/m2]
      real      Evg_sv(klonv)                 ! Emissivity of Vegetation+Snow
      real      Eso_sv(klonv)                 ! Emissivity of       Soil+Snow
      real      tau_sv(klonv)                 ! Transmited Radiation Fraction
      real      rrMxsv(klonv)                 ! Canopy Maximum Intercepted Rain
      real      LAIesv(klonv)                 ! effective LAI for transpirati.
      real      LAI_sv(klonv)                 ! corrected LAI in case of snow
      real      glf_sv(klonv)                 ! Green  Leaf Fraction
      real      Sigmsv(klonv)                 ! Canopy Ventilation  Factor
      real      HSv_sv(klonv)                 ! Sensible Heat Flux  [W/m2]
      real      HLv_sv(klonv)                 ! Latent   Heat Flux  [W/m2]
      real      HSs_sv(klonv)                 ! Sensible Heat Flux (t)
      real      HLs_sv(klonv)                 ! Latent   Heat Flux (t)
      real      sqrCm0(klonv)                 ! in Neutral Drag Coef.Moment.
      real      sqrCh0(klonv)                 ! in Neutral Drag Coef.Heat
      real      Lx_H2O(klonv)                 ! Latent Heat of Vaporiz./Sublim.
      real      ram_sv(klonv)                 ! Aerodyn.Resistance (Moment.)
      real      rah_sv(klonv)                 ! Aerodyn.Resistance (Heat)
      real      Fh__sv(klonv)                 ! Stability Function
      real      dFh_sv(klonv)                 ! Stability Function (Deriv.)
      real      Evp_sv(klonv)                 ! Evaporation        [kg/m2]
      real      EvT_sv(klonv)                 ! Evapotranspiration [kg/m2]
      real      LSdzsv(klonv)                 ! Land/Sea Vert. Discretiz. Fact.
      real      Tsrfsv(klonv)                 ! Surface    Temperature
      real      sEX_sv(klonv,-nsol:nsno+1)    ! Verticaly Integr.Extinct.Coef.
      real      zzsnsv(klonv,    0:nsno)      ! Snow  Pack Thickness      [m]
      real      psi_sv(klonv,-nsol:0   )      ! Soil   Water        Potential
      real      Khydsv(klonv,-nsol:0   )      ! Soil   Hydraulic    Conductiv.
      real      Rootsv(klonv,-nsol:0)         ! Root Water Pump      [kg/m2/s]
      real      EExcsv(klonv)                 ! Energy in Excess, current
      common  /ySISVAT_R/
     .          albisv,albssv,SoCasv,SoSosv,IRv_sv,Evg_sv,Eso_sv,tau_sv
     .         ,rrMxsv,LAIesv,LAI_sv,glf_sv
     .         ,Sigmsv,HSv_sv,HLv_sv,HSs_sv,HLs_sv
     .         ,sqrCm0,sqrCh0,Lx_H2O,ram_sv,rah_sv,Evp_sv,EvT_sv
     .         ,LSdzsv,Tsrfsv,sEX_sv,zzsnsv
     .         ,psi_sv,Khydsv,Rootsv,EExcsv
