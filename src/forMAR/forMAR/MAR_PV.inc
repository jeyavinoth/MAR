      real            adugPV,adubPV       ,advgPV       ,advbPV       ,
     .                                     afgzPV    ,   afpsPV    ,
     .                ulscPV    ,          vlscPV    
      common /PV_aux/ adugPV,adubPV(mx,mz),advgPV(mx,mz),advbPV(mx,mz),
     .                                     afgzPV(mx),   afpsPV(mx),
     .                ulscPV(mz),          vlscPV(mz)
C +...................adugPV : factor of u(sounding) ==> ug 
C +                   adubPV : factor of ug          ==> u(large scale) 
C +                   advgPV : factor of ug**2       ==> vg 
C +                   advbPV : Initial Value of vL
C +                   adiaPV : Topographic Contribution  to PGF
C +                   ulscPV : Large Scale Value of u at Sounding Point 
C +                                                  and Previous Time Step
C +                   vlscPV : Large Scale Value of v at Sounding Point 
C +                                                  and Previous Time Step
