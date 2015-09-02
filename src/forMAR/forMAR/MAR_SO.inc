      logical       sodata
      common/solar0/sodata
C +
      real          ocldSO,
     .                FSOm      ,  FSOd      ,  FSOn      , dFSOn,
     .              htngSO
      common/solar1/ocldSO(mzir),
     .                FSOm(mzir),  FSOd(mzir),  FSOn(mzir), dFSOn(mzir),
     .              htngSO(mzir)
C +...                ocldSO: Simple Scattering Cloud Albedo           
C +                   FSOm  : Upward   SOLAR Flux                     (W/m2)
C +                   FSOd  : Downward SOLAR Flux                     (W/m2)
C +                   FSOn  : Net      SOLAR Flux                     (W/m2)
C +                   htngSO:          SOLAR Heating                   (K/s)
C +
