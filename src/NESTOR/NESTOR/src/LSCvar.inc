
C +---Input model
C +   -----------

      CHARACTER*3 LSCmod
      LOGICAL     SPHgrd,REGgrd
    
C +...LSCmod : acronym of the large-scale model
C +...SPHgrd : stereographic grid
C +...REGgrd : regular grid (in latitude and longitude) 


C +---Horizontal and vertical grid
C +   ----------------------------

      REAL LSC1Dx(ni   ),LSC1Dy(   nj),LSC__z(nk),
     .                   LSC1Vy(  njv),
     .     LSC__x(ni,nj),LSC__y(ni,nj),LSCgdz(nk)

C +...LSC__x : X-coordinates
C +...LSC1Dx : X-coordinates (for regular input grid)
C +...LSC__y : Y-coordinates
C +...LSC1Dy : Y-coordinates (for regular input grid)
C +...LSC1Vy : Y-coordinates (for regular grid, wind)
C +...LSC__z : Z-coordinates (hybrid)
C +...LSCgdz : Vertical grid of the LSC model


C +---2-D surface variables
C +   ---------------------

      REAL LSC_st(ni,nj),LSCdst(ni,nj),LSC_sw(ni,nj),
     .     LSCdsw(ni,nj),LSC_sp(ni,nj),LSC_dt(ni,nj),
     .     LSC_sh(ni,nj),LSC_pr(ni,nj),LSCppr(ni,nj),
     .     LSC_sn(ni,nj),LSCtcc(ni,nj),LSCuts(ni,nj),
     .     LSCsic(ni,nj),LSCsst(ni,nj),LSClsm(ni,nj)

C +...LSC_st : soil or sea surface temperature
C +...LSCdst : deep soil temperature
C +...LSC_sw : soil wetness
C +...LSCdsw : deep soil wetness
C +...LSC_sp : surface pressure
C +...LSC_dt : temperature diff. between 1st lev. and surf.
C +...LSC_sh : surface elevation
C +...LSC_pr : rain precipitation at the current  time step
C +...LSCppr : rain precipitation at the previous time step
C +...LSC_sn : snow precipitation
C +...LSCtcc : total cloud cover
C +...LSCuts : surface heat flux
C +...LSCsic : Sea Ice Fraction
C +...LSCsst : Sea Surface Temperature
C +...LSClsm : Land Sea Mask


C +---3-D atmospheric variables (storred on 1 level = 2D)
C +   -------------------------

      REAL LSC__p(ni,nj),
     .     LSC__u(ni,nj),LSC__v(ni,njv),LSC__w(ni,nj),
     .     LSC_pt(ni,nj),LSC__t(ni,nj) ,LSC_qv(ni,nj),
     .     LSCtke(ni,nj),LSC_qt(ni,nj) ,LSCtmp(ni,nj),
     .     LSCtm2(ni,nj),LSCgv1(ni,njv),LSCgv2(ni,njv),
     .     LSC_rh(ni,nj)

C +...LSC__p : pressure
C +...LSC__u : U-wind
C +...LSC__v : V-wind
C +...LSC__w : W-wind
C +...LSC_pt : potential temperature
C +...LSC__t : real temperature
C +...LSC_qv : specific humidity
C +...LSCtke : turbulent kinetic energy
C +...LSC_qt : total cloud water
C +...LSCtmp : temporary array
C +...LSCtm2 : temporary array


      COMMON/LSCvar_c/LSCmod
      
      COMMON/LSCvar_l/SPHgrd,REGgrd
      
      COMMON/LSCvar_r/LSC1Dx,LSC1Dy,LSC__z,LSC1Vy,LSC__x,LSC__y,
     .                LSCgdz,LSC_st,LSCdst,LSC_sw,LSCdsw,LSC_sp,
     .                LSC_dt,LSC_sh,LSC_pr,LSCppr,LSC_sn,LSCtcc,
     .                LSCuts,LSC__p,LSC__u,LSC__v,LSC__w,LSC_pt,
     .                LSC__t,LSC_qv,LSCtke,LSC_qt,LSCtmp,LSCtm2,
     .                LSCgv1,LSCgv2,LSCsic,LSCsst,LSC_rh,LSClsm
   
