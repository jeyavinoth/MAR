
C +---Nested model
C +   ------------

      CHARACTER*3 NSTmod

C +...NSTmod : acronym of the nested model


C +---Horizontal and vertical grid
C +   ----------------------------

      REAL NST__x(mx,my),NST__y(mx,my),
     .     NSTgdx(mx)   ,NSTgdy(my)   ,NSTgdz(mz)      ,
     .     NST_dx       

      LOGICAL NSTinc(mx,my,nbdom)

C +...NST__x : X-coordinates (longitude)
C +...NST__y : Y-coordinates (latitude)
C +...NSTgdx : simple X-grid (specific to each model)
C +...NSTgdy : simple Y-grid (specific to each model)
C +...NSTgdz : simple Z-grid (specific to each model)
C +...NST_dx : horizontal resolution for regular grid
C +...NSTinc : Check the location (continent) of a grid point


C +---2-D surface variables
C +   ---------------------

      REAL NST_st(mx,my),NSTdst(mx,my),NST_sp(mx,my),
     .     NST_sh(mx,my),NST_pr(mx,my),NST_sn(mx,my),
     .     NST_z0(mx,my),NST_r0(mx,my),NST_d1(mx,my),
     .     NSTalb(mx,my),NSTeps(mx,my),NSTres(mx,my),
     .     NSTch0(mx,my),NSTIpr(mx,my),NSTewc(mx,my),
     .     NSTsst(mx,my),NSTuts(mx,my),NSTpr1(mx,my),
     .     NSTpr2(mx,my),NSTpr3(mx,my),NSTsic(mx,my)

      INTEGER NSTveg(mx,my,nvx),NSTsvt(mx,my,nvx),
     .        NSTiwf(mx,my    ),NSTsol(mx,my    ),
     .        NSTtex(mx,my    )

      REAL    NSTlai(mx,my,nvx),NSTglf(mx,my,nvx),
     .        NSTvfr(mx,my,nvx),NSTsfr(mx,my,nvx),
     .        NSTdsa(mx,my    ),NSTlmx(mx,my,nvx),
     .        NSTfrc(mx,my    ),NSTdv1(mx,my    ),
     .        NSTdv2(mx,my    ),NSTndv(mx,my    ),
     .        NSTice(mx,my    )

C +...NST_st : soil or sea surface temperature
C +...NSTsst : sea surface temperature (Reynolds)
C +...NSTdst : deep soil temperature
C +...NST_sp : surface pressure
C +...NST_sh : surface elevation
C +...NSTpr1 : desaggregated precipitation (without conservation)
C +...NSTpr2 : desaggregated precipitation (with global conservation)
C +...NSTpr3 : desaggregated precipitation (with local and global conservation)
C +...NSTsic : Sea-Ice   Fraction
C +...NSTIpr : rain precipitation (non desaggregated)
C +...NSTewc : equivalent water content (water vapor)
C +...NST_sn : snow precipitation
C +...NSTsol : soil types (water,ice,snow,land,...)
C +...NST_z0 : roughness length for momentum
C +...NST_r0 : roughness length for heat
C +...NSTtex : soil texture
C +...NST_d1 : surface heat capacity (Deardorff, 1978)
C +...NSTalb : surface albedo
C +...NSTeps : surface IR emissivity
C +...NSTres : aerodynamic resistance
C +...NSTch0 : bulk aerodynamic coefficient air/surface
C +...         humidity flux
C +...NSTveg : vegetation type (IGBP classification)
C +...NSTvfr : fraction of vegetation grid cells (IGBP)
C +...NSTsvt : vegetation type (SVAT classification)
C +...NSTsfr : fraction of vegetation grid cells (SVAT)
C +...NSTfrc : fraction of vegetation cover (from NDVI)
C +...NSTndv : NDVI index
C +...NSTdv1 : minimum NDVI index
C +...NSTdv2 : maximum NDVI index
C +...NSTlai : leaf area index
C +...NSTglf : green leaf fraction
C +...NSTdsa : dry soil albedo
C +...NSTiwf : 0=no water flux, 1=free drainage
C +...NSTuts : surface heat flux


C +---2.5-D surface variables
C +   -----------------------

      REAL NST_ts(mx,my,nvx,nsl),NST_sw(mx,my,nvx,nsl)

C +...NST_ts : soil temperature
C +...NST_sw : soil moisture content


C +---3-D atmospheric variables
C +   -------------------------
      REAL NST_hp(mx,my,mz)
C +...NST_hp : Local vertic coordinate for interpolation

      REAL NST__u(mx,my,mz),NST__v(mx,my,mz),NST__w(mx,my,mz),
     .     NST_pt(mx,my,mz),NST__t(mx,my,mz),NST_qv(mx,my,mz),
     .     NST_zz(mx,my,mz),NST__p(mx,my,mz),NSTtke(mx,my,mz),
     .     NST_qt(mx,my,mz),NSTtmp(mx,my,mz+1),NST_rh(mx,my,mz)

C +...NST__u : U-wind
C +...NST__v : V-wind
C +...NST__w : W-wind
C +...NST_pt : potential temperature
C +...NST__t : real temperature
C +...NST_qv : specific humidity
C +...NST_rh : relative humidity
C +...NST_zz : geopotential height
C +...NST__p : pressure at each level
C +...NSTtke : turbulent kinetic energy
C +...NST_qt : total cloud water content
C +...NSTtmp : temporary array


      COMMON/NSTvar_c/NSTmod
      
      COMMON/NSTvar_i/NSTveg,NSTsvt,NSTiwf,NSTsol,
     .                NSTtex
      
      COMMON/NSTvar_l/NSTinc
      
      COMMON/NSTvar_r/NST__x,NST__y,NSTgdx,NSTgdy,NSTgdz,NST_dx,
     .                NST_st,NSTdst,NST_sp,NST_sh,NST_pr,NST_sn,
     .                NST_z0,NST_r0,NST_d1,NSTalb,NSTeps,NSTres,
     .                NSTch0,NSTIpr,NSTewc,NSTsst,NSTuts,NSTpr1,
     .                NSTpr2,NSTpr3,NSTlai,NSTglf,NSTdsa,NSTlmx,
     .                NSTfrc,NSTdv1,NSTdv2,NSTndv,NST_ts,NST_sw,
     .                NST_hp,NST__u,NST__v,NST__w,NST_pt,NST__t,
     .                NST_qv,NST_zz,NST__p,NSTtke,NST_qt,NSTtmp,
     .                NSTsic,NST_rh,NSTvfr,NSTsfr,NSTice
 
