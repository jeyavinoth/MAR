C +...BEGIN    "MAR0SV.inc"
      integer         islpSV(-nsol:0)
      integer         isnpSV(      nsno)
      integer         islmSV(-nsol:0)
      common  /SoI0SV/islpSV,isnpSV,islmSV
      integer         nkhy
      parameter      (nkhy=50)
      real            Implic,Explic
      real            dzmiSV(-nsol:0)          ! dz_(i-1/2)
      real            dzi_SV(-nsol:0)          ! dz_(i-1)/(dz_(i)+dz_(i-1))
      real            dziiSV(-nsol:0)          ! dz_(i)  /(dz_(i)+dz_(i-1))
      real            dtz_SV(-nsol:0)          ! dt / dz
      real            dz78SV(-nsol:0)          ! 7/8 (dz)
      real            dz34SV(-nsol:0)          ! 3/4 (dz)
      real            dz_8SV(-nsol:0)          ! 1/8 (dz)
      real            dzAvSV(-nsol:0)          ! 1/8dz_(-1)+3/4dz+1/8dz_(+1)
      real            OcndSV                   ! Swab Ocean / Soil Ratio 
      real            RF__SV( 0:nvgt,-nsol:0)  ! Root Fraction
      real            rocsSV( 0:nsot)          ! Soil Contribution to (ro c)_s
      real            etamSV( 0:nsot)          ! Soil Minimum Humidity
      real            s1__SV( 0:nsot)          ! ... X eta**( b+2), DR97(3.36)
      real            s2__SV( 0:nsot)          ! ... X eta**(2b+3), DR97(3.35)
      real            aKdtSV( 0:nsot, 0:nkhy)  ! Khyd=a*eta+b: a * dt
      real            bKdtSV( 0:nsot, 0:nkhy)  ! Khyd=a*eta+b: b * dt
      common  /SoR0SV/Implic,Explic,dzmiSV,dzi_SV,dziiSV,
     .                dtz_SV,dz78SV,dz34SV,dz_8SV,dzAvSV,OcndSV,
     .                RF__SV,rocsSV,etamSV,s1__SV,s2__SV,aKdtSV,bKdtSV
C +...END      "MAR0SV.inc"
