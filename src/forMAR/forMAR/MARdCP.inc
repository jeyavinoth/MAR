C +..BEGIN    "MARdCP.inc"
C +--Col de Porte specific Constants
C +  ===============================
C +
      logical         ColPrt                ! Col de Porte Switch
      common/ColPrt_l/ColPrt
C +
C +
C +--Fractions of total solar irradiances in 3 spectral intervals
C +  ------------------------------------------------------------
C +
      real            Dr_1SN,Dr_2SN,Dr_3SN  ! Direct  Radiation
      real            Df_1SN,Df_2SN,Df_3SN  ! Diffuse Radiation, Clear  Sky
      real            Dfc1SN,Dfc2SN,Dfc3SN  ! Diffuse Radiation, Cloudy Sky
      common/ColPrt_r/Dr_1SN,Dr_2SN,Dr_3SN
     .               ,Df_1SN,Df_2SN,Df_3SN
     .               ,Dfc1SN,Dfc2SN,Dfc3SN
C +
      real            DirSol,DifSol,TotSol,Clouds
      common/ColPrt_v/DirSol,DifSol,TotSol,Clouds
C +
C +..END      "MARdCP.inc"
