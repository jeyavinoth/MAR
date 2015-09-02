C +...KLON:  defined in MARdim.inc
C +...KLEV:  defined in MARdim.inc
C +
      integer Kstep1(KLON)             ! convective counter
      integer K_CbT1(KLON)             ! cloud top  level
      integer K_CbB1(KLON)             ! cloud base level
      common /CVAgen_MNH_int/Kstep1,K_CbT1,K_CbB1
C +
      integer                  KTCCH0
      parameter(               KTCCH0=1)
      REAL    P_CH_0(KLON,KLEV,KTCCH0)
      REAL    PdCH_1(KLON,KLEV,KTCCH0)
      common /CVAgen_MNH_chi/P_CH_0,PdCH_1
C +
      REAL    PdTa_1(KLON,KLEV)
      REAL    PdQa_1(KLON,KLEV)
      REAL    PdQw_1(KLON,KLEV)
      REAL    PdQi_1(KLON,KLEV)
      REAL    Pdrr_1(KLON)
      REAL    Pdss_1(KLON)
      REAL    PuMF_1(KLON,KLEV)        !   Upward        Mass Flux
      REAL    PdMF_1(KLON,KLEV)        ! Downward        Mass Flux
      REAL    Pfrr_1(KLON,KLEV)        ! Liquid Precipitation Flux
      REAL    Pfss_1(KLON,KLEV)        !  Solid Precipitation Flux
      REAL    Pcape1(KLON)             !   CAPE [J/kg]
      common /CVAgen_MNH_rea/PdTa_1,PdQa_1,PdQw_1,PdQi_1,Pdrr_1,Pdss_1,
     .                       PuMF_1,PdMF_1,Pfrr_1,Pfss_1,Pcape1
