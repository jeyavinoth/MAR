C +
C +  ==========
C +--MAR_NH.inc
C +  ==========
C +
      real           ua0_NH(mx,my,mz),va0_NH(mx,my,mz),wa0_NH(mx,my,mz)
      real           wairNH(mx,my,mz),pairNH(mx,my,mz)
      real           dt_ONH,dtxONH,dtyONH
     .              ,gamaNH,  c2NH,  taNH
      common/nonhyd/ ua0_NH          ,va0_NH          ,wa0_NH
     .              ,wairNH          ,pairNH
     .              ,dt_ONH,dtxONH,dtyONH
     .              ,gamaNH,  c2NH,  taNH
C +... ua0_NH : Horizontal (x-Direction) Wind Speed  (m/s) / Previous Time Step
C +    va0_NH : Horizontal (y-Direction) Wind Speed  (m/s) / Previous Time Step
C +    wa0_NH : Vertical   (z-Direction) Wind Speed  (m/s) / Previous Time Step
C +    wairNH : Vertical Non-Hydrostatic Wind Speed  (m/s)
C +    pairNH :          Non-Hydrostatic Pressure    (kPa)
C +    dt_ONH : Non-Hydrostatic    Time  Step          (s)
C +    dtxONH = dt_ONH / dx                          (s/m)
C +    dtyONH = dt_ONH / dy                          (s/m)
C +    gamaNH : Cp / Cv
C +      c2NH : Prescribed Squared Sound Speed     (m2/s2)
C +      taNH : Reference  Surface Air   Temperature   (K)
C +
