      logical       IRdata
      common/infra0/IRdata
C +
      real          ecldIR      ,epsIR       ,  
     .                FIRm      ,  FIRd      ,  FIRn,
     .              htngIR
      common/infra1/ecldIR(mzir),epsIR (mzir),  
     .                FIRm(mzir),  FIRd(mzir),  FIRn(mzir),
     .              htngIR(mzir)
C +...                ecldIR:        Cloud IR Emissivity        
C +                   epsIR : Local  Cloud IR Emissivity               
C +                   FIRm  : Upward   IR    Flux                     (W/m2)
C +                   FIRd  : Downward IR    Flux                     (W/m2)
C +                   FIRn  : Net      IR    Flux                     (W/m2)
C +                   htngIR:          IR    Heating                   (K/s)
