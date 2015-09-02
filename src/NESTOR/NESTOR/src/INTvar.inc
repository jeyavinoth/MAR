
C +---Horizontal and vertical grid
C +   ----------------------------

      REAL INT__x(mx,my),INT__y(mx,my),INT__z(mx,my,nk)

C +...INT__x : X-coordinates
C +...INT__y : Y-coordinates
C +...INT__z : Z-coordinates (hybrid)


C +---2-D surface variables
C +   ---------------------

      REAL INT_st(mx,my),INTdst(mx,my),INT_sw(mx,my),
     .     INTdsw(mx,my),INT_sp(mx,my),INT_dt(mx,my),
     .     INT_sh(mx,my),INT_pr(mx,my),INT_sn(mx,my),
     .     INTsol(mx,my),INTtcc(mx,my),INTuts(mx,my),
     .     INTsic(mx,my),INTsst(mx,my)

C +...INT_st : soil or sea surface temperature
C +...INTdst : deep soil temperature
C +...INT_sw : soil wetness
C +...INTdsw : deep soil wetness
C +...INT_sp : surface pressure
C +...INT_dt : temperature diff. between 1st lev and surf.
C +...INT_sh : surface elevation
C +...INT_pr : rain precipitation
C +...INT_sn : snow precipitation
C +...INTsol : soil types (water,ice,snow,land,...)
C +...INTtcc : total cloud cover
C +...INTuts : surface heat flux
C +...INTsic : Sea-Ice Fraction
C +...INTsic : Sea Surface Temperature

C +---3-D atmospheric variables
C +   -------------------------

      REAL INT__u(mx,my,nk+1),INT_pt(mx,my,nk+1),INT__p(mx,my,nk),
     .     INT__v(mx,my,nk+1),INT_qv(mx,my,nk+1),INTtke(mx,my,nk+1),
     .     INT_qt(mx,my,nk+1),INT_rh(mx,my,nk+1)

C +...INT__u : U-wind
C +...INT__v : V-wind
C +...INT_pt : potential temperature
C +...INT_qv : specific humidity
C +...INT_rh : relative humidity
C +...INT__p : pressure
C +...INTtke : turbulent kinetic energy
C +...INT_qt : total cloud water content


C +---Temporary arrays
C +   ----------------

      REAL INT1Dz(nk+1),INT1Dp(nk+1),INT1Du(nk+1),INT1Dv(nk+1),
     .     INT1Dt(nk+1),INT1Dq(nk+1),INT1De(nk+1),INT1Dh(nk+1)


      common/INTvar_r/INT__x,INT__y,INT__z,INT_st,INTdst,INT_sw,
     .                INTdsw,INT_sp,INT_dt,INT_sh,INT_pr,INT_sn,
     .                INTsol,INTtcc,INTuts,INT__u,INT_pt,INT__p,
     .                INT__v,INT_qv,INTtke,INT_qt,INTsic,INT_rh


