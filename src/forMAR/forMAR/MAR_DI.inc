      real           pkttDI       ,qvtoDI
     .              ,uairDI       ,vairDI
     .              ,ect_DI       ,eps_DI
      common/dirtbc/ pkttDI(mx,my),qvtoDI(mx,my)
     .              ,uairDI(mx,my),vairDI(mx,my)
     .              ,ect_DI(mx,my),eps_DI(mx,my)
C +... pkttDI: potential temperature at the model top divided by p00**(R/Cp)
C +    qvtoDI: specific humidity     at the model top                 (kg/kg)
C +
