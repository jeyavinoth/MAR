      block data PHYdat
C +
C +------------------------------------------------------------------------+
C | MAR PHYSIC                                              8-07-1996  MAR |
C |   Block Data PHYdat is used to define physical constants used in   MAR |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
C +--General Variables
C +  =================
C +
c #DP include 'MAR_r8.inc'
C +
      include 'MARphy.inc'
C +
      include 'MARdim.inc'
      include 'MAR_GE.inc'
      include 'MARSND.inc'
C +
C +
C +--Numeric   Constants
C +  ===================
C +
      data  izr /0     /
      data  iun /1     /
      data zero /0.0d+0/
      data demi /0.5d+0/
      data unun /1.0d+0/
      data epsi /1.0d-6/
      data thous/1.0d+3/
C +
C +
C +--Earth     Constants
C +  ===================
C +
      data earthr/6371.229d+3/
C +...     earthr: Earth Radius
C +
      data earthv/ 729.217d-7/
C +...     earthv: Earth Angular Velocity
C +
      data gravit/   9.81d0/
C +...     gravit: Gravity Acceleration
C +
C +
C +--Dynamical Constants
C +  ===================
C +
      data akmol /   1.35d-5/
C +...     akmol : Air Viscosity
C +
      data vonkar/   0.40d0/
C +...     vonkar: von Karman constant
C +
C +
C +--Thermodynamical Constants (Atmosphere)
C +  ======================================
C +
      data ra    / 287.d0/
C +...     ra    : Perfect Gas Law  Constant (J/kg/K)
C +
      data cp    /1004.d0/
C +...     cp    : Air Specific Heat         (J/kg/K)
C +
      data cap   /   0.285856574d0/
C +...     cap   = R / Cp
C +
      data pcap  /   3.730037070d0/
C +...     pcap  = 100 ** (R / Cp)
C +
      data epsq  /   3.00d-6/
C +...     epsq  : Minimum Water Vapor Content
C +
      data h2olv /   2.5000d+6/
C +...     h2olv : Latent Heat of Vaporisation for Water (J/kg)
C +
      data h2ols /   2.8336d+6/
C +...     h2ols : Latent Heat of Sublimation  for Ice   (J/kg)
C +
      data stefan/   5.67d-8/
C +...     stefan: Stefan-Bolstzman Constant (W/m2/K4)
C +
C +
C +--Thermodynamical Constants (Sea)
C +  ===============================
C +
      data tfrwat/ 271.20d0/
C +...     tfrwat: Sea-Water freezing Temperature  (K)
C +                (definition: cfr. data in .main.)
C +
      data siclf / 302.00d6/
C +...     siclf : Sea-Ice Heat of Fusion       (J/m3)
C +
      data cdsice/   2.04d0/
C +...     cdsice: Conductivity of  Sea Ice    (W/m/K)
C +
      data hic0  /   1.00d0/
C +...     hic0  : Sea-Ice Initial Thickness       (m)
C +                (Ross Sea Typical Value)
C +
      data rowat /1000.00d0/
C +...     rowat : Density of Water            (kg/m3)
C +
      data cwwat /4186.00d0/
C +...     cwwat : Heat Capacity of Water     (J/kg/K)
C +
      data fracoh/   0.75d0/
C +...     fracoh=   1.00 - 0.25
C +***Hibler (1984):   25% of cooling => Oceanic Heat Flux (ANTARCTIC Ocean)
C +
C +
C +--Thermodynamical Constants (Ice)
C +  ===============================
C +
      data roice  /   920.d0/
C +...     roice:  Density of Pure Ice          (kg/m3)
C +
      data cdice  /   2.51d0/
C +...     cdice:  Conductivity of pure Ice     (W/m/K)
C +
C +
C +--Thermodynamical Constants (Snow)
C +  ================================
C +
      data  csnow/2105.00d+0/
C +...      csnow:  Heat Capacity of Snow                        [J/kg/K]
C +***             (see Loth et al. 1993, JGR 98 D6, 2.2.2 2e para p.10453)
C +
      data tfsnow/ 273.16d+0/
C +...     tfsnow: Snow      melting  Temperature  (K)
C +
      data r0sno /   5.00d+1/
C +...     r0sno :  Fresh  Snow Density        (kg/m3)
C +
      data blsno /   3.30d+2/
C +...     blsno :  Blowed Snow Density        (kg/m3)
C +
      data snolf /   3.34d+5/
C +...     snolf : Latent  Heat of Fusion/ Snow (J/kg)
C +
C +
C +--Standard Atmosphere (Winter / Middle Latitudes)
C +  ===============================================
C +
      data iSND / 1 /
      data jSND / 1 /
C +
      data zSND / -143.,    0., 1001., 1993., 2992., 3993., 4994.,
     .            5983., 6978., 7988., 8984., 9970.,10968.,11975.,
     .           12966.,13949.,14945.,15932.,16950.,17900.,18914.,
     .           19884.,20894.,21933.,22985.,23799.,24990.,29928.,
     .           35068.,38589.,46673.,49408.,49583.,49761.,49948.,
     .           50132.,50324.,50521.,50723.,50930.,51143.       ,
     .            -143.,    0., 1001., 1993., 2992., 3993., 4994.,
     .            5983., 6978., 7988., 8984., 9970.,10968.,11975.,
     .           12966.,13949.,14945.,15932.,16950.,17900.,18914.,
     .           19884.,20894.,21933.,22985.,23799.,24990.,29928.,
     .           35068.,38589.,46673.,49408.,49583.,49761.,49948.,
     .           50132.,50324.,50521.,50723.,50930.,51143.       /
C +
      data tSND / 277.0, 272.2, 268.7, 265.2, 261.7, 255.7, 249.7,
     .            243.7, 237.7, 231.7, 225.7, 219.7, 219.2, 218.7,
     .            218.2, 217.7, 217.2, 216.7, 216.2, 215.7, 215.2,
     .            215.2, 215.2, 215.2, 215.2, 215.2, 215.2, 217.4,
     .            227.8, 243.2, 258.5, 265.7, 265.7, 265.7, 265.7,
     .            265.7, 265.7, 265.7, 265.7, 265.7, 265.7       ,
     .            277.0, 272.2, 268.7, 265.2, 261.7, 255.7, 249.7,
     .            243.7, 237.7, 231.7, 225.7, 219.7, 219.2, 218.7,
     .            218.2, 217.7, 217.2, 216.7, 216.2, 215.7, 215.2,
     .            215.2, 215.2, 215.2, 215.2, 215.2, 215.2, 217.4,
     .            227.8, 243.2, 258.5, 265.7, 265.7, 265.7, 265.7,
     .            265.7, 265.7, 265.7, 265.7, 265.7, 265.7       /
C +
      data qSND /.27e-2,.27e-2,.21e-2,.17e-2,.13e-2,.80e-3,.51e-3,
     .           .32e-3,.14e-3,.67e-4,.35e-4,.18e-4,.20e-4,.20e-4,
     .           .70e-5,.45e-5,.40e-5,.39e-5,.40e-5,.42e-5,.48e-5,
     .           .51e-5,.68e-5,.81e-5,.10e-4,.13e-4,.17e-4,.20e-4,
     .           .14e-4,.10e-4,.14e-4,.69e-5,.70e-5,.72e-5,.74e-5,
     .           .75e-5,.77e-5,.79e-5,.81e-5,.83e-5,.86e-5       ,
     .           .27e-2,.27e-2,.21e-2,.17e-2,.13e-2,.80e-3,.51e-3,
     .           .32e-3,.14e-3,.67e-4,.35e-4,.18e-4,.20e-4,.20e-4,
     .           .70e-5,.45e-5,.40e-5,.39e-5,.40e-5,.42e-5,.48e-5,
     .           .51e-5,.68e-5,.81e-5,.10e-4,.13e-4,.17e-4,.20e-4,
     .           .14e-4,.10e-4,.14e-4,.69e-5,.70e-5,.72e-5,.74e-5,
     .           .75e-5,.77e-5,.79e-5,.81e-5,.83e-5,.86e-5       /
C +
      data pSND / 1036., 1018.,  897.,  790.,  694.,  608.,  531.,
     .             463.,  402.,  347.,  299.,  257.,  220.,  188.,
     .             161.,  138.,  118.,  101.,   86.,   74.,   63.,
     .              54.,   46.,   39.,   33.,   29.,   24.,   11.,
     .               5.,    3.,    1.,   0.7, 0.685, 0.669, 0.654,
     .            0.637, 0.622, 0.606, 0.591, 0.576, 0.560       ,
     .            1036., 1018.,  897.,  790.,  694.,  608.,  531.,
     .             463.,  402.,  347.,  299.,  257.,  220.,  188.,
     .             161.,  138.,  118.,  101.,   86.,   74.,   63.,
     .              54.,   46.,   39.,   33.,   29.,   24.,   11.,
     .               5.,    3.,    1.,   0.7, 0.685, 0.669, 0.654,
     .            0.637, 0.622, 0.606, 0.591, 0.576, 0.560       /
C +
C +
C +--Time Constants
C +  ==============
C +
      data (njmoGE(n),n=0,12)
     .     /0,31,28,31,30, 31, 30, 31, 31, 30, 31, 30, 31/
      data (njyrGE(n),n=0,12)
     .     /0, 0,31,59,90,120,151,181,212,243,273,304,334/
C +...      njmoGE: Nb of Days in each Month of the Year
C +         njyrGE: Nb of Days since   Begin of the Year,
C +                            before  Current Month
C +
      data (labmGE(n),n=0,12)
     .     /'---','Jan','Feb','Mar','Apr','May','Jun',
     .            'Jul','Aug','Sep','Oct','Nov','Dec'/
C +
C +
C +--Characters
C +  ==========
C +
      data labnum/'0','1','2','3','4','5','6','7','8','9'/
C +
      end
C +
      subroutine GRDgeo(maptyp)
C +
C +------------------------------------------------------------------------+
C | MAR GRID                                               17-07-1996  MAR |
C |   SubRoutine GRDgeo computes the Latitudes, Longitudes and             |
C |                              the Time Zone of each Grid Point          |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |    INPUT: imez,jmez    : Indices of the MAR Domain Center              |
C |    ^^^^^^ GEddxx       : (2-D): x-Axis      Direction                  |
C |                          (3-D): South-North Direction along            |
C |                                 90E, 180E, 270E or 360E Meridians      |
C |           GElat0       : Latitude  of (0,0) in  MAR              (deg) |
C |           GElon0       : Longitude of (0,0) in  MAR              (deg) |
C |                                                                        |
C |   HYPOT.: maptyp = 0   : Polar   Stereogr. Project. (SOUTH HEMISPHERE) |
C |   ^^^^^^^ maptyp = 1   : Oblique Stereogr. Project. (ALL    LATITUDES) |
C |           maptyp = 2   : Lambert Comformal Project. (ALL LAT, 3D only) |
C |                                                                        |
C |   OUTPUT: GElatr(mx,my): Latitude  of the (x,y) MAR coordinate   (rad) |
C |   ^^^^^^^ GElonh(mx,my): Longitude of the (x,y) MAR coordinate     (h) |
C |           itizGE(mx,my): Time Zone                                     |
C |           fcorDY(mx,my): Coriolis Parameter (Variable/only 3-D Domain) |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
C +--General Variables
C +  =================
C +
c #DP include 'MAR_r8.inc'
C +
      include 'MARphy.inc'
C +
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
C +
      include 'MAR_DY.inc'
C +
C +
C +--GEOGRAPHIC Coordinates
C +  ======================
C +
C +
C +--1-D and 2-D Cases
C +  -----------------
C +
      if (mmy.eq.1) then                                               ! CTR
C +
         argrot = (GEddxx-90.d0)*degrad
         cosrot =  cos(argrot)
         sinrot =  sin(argrot)
C +
        do 21 j=1,my
        do 21 i=1,mx
         xxmar = cosrot*(i-imez)*dx + sinrot*(j-jmez)*dx
         yymar = cosrot*(j-jmez)*dx - sinrot*(i-imez)*dx
C +
C +      ***********
         call GRDstr(xxmar,yymar,GElon0,GElat0,GElon,GElat)
C +      ***********
C +
         GElatr(i,j) =  GElat
         GElonh(i,j) =  GElon
C +
 21     continue
C +
C +
C +--3-D         Cases
C +  -----------------
C +
      else                                                             ! CTR
C +
C +- ANTARCTICA (Polar   Stereographic Projection is assumed)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if (maptyp.eq.0) then                                           ! CTR
        ddista = earthr * 2.d0 * tan((45.d0+GElat0*0.5d0)*degrad)
        xdista = ddista        * cos((90.d0-GElon0)      *degrad)
        ydista = ddista        * sin((90.d0-GElon0)      *degrad)
        do 31 j=1,my
        do 31 i=1,mx
        if (abs(GEddxx- 90.d0).lt.epsi) then
         xxmar = (i-imez)*dx
         yymar = (j-jmez)*dy
        end if
        if (abs(GEddxx       ).lt.epsi) then
         xxmar = (j-jmez)*dy
         yymar =-(i-imez)*dx
        end if
        if (abs(GEddxx-270.d0).lt.epsi) then
         xxmar =-(i-imez)*dx
         yymar =-(j-jmez)*dy
        end if
        if (abs(GEddxx-180.d0).lt.epsi) then
         xxmar =-(j-jmez)*dy
         yymar = (i-imez)*dx
        end if
C +
         xxmar = xxmar + xdista
         yymar = yymar + ydista
C +
            ddista      = sqrt(xxmar*xxmar+yymar*yymar)
            GElatr(i,j) =-0.5d0*pi +2.d0*atan(ddista*0.5d0/earthr)
        if(abs(xxmar).gt.zero) then
            GElonh(i,j) = atan(yymar/xxmar)
         if   (xxmar.lt.zero)
     .      GElonh(i,j) = GElonh(i,j) + pi
C +
            GElonh(i,j) =   0.5d0 * pi - GElonh(i,j)
         if(GElonh(i,j).gt.         pi)
     .      GElonh(i,j) =  -2.0d0 * pi + GElonh(i,j)
         if(GElonh(i,j).lt.        -pi)
     .      GElonh(i,j) =   2.0d0 * pi + GElonh(i,j)
C +
        else
         if   (yymar.gt.zero) then
            GElonh(i,j) =   0.0d0
         else
            GElonh(i,j) =           pi
         end if
        end if
C +...  transformation stereographic coordinates (center = South Pole)
C +                 -> spherical     coordinates
C +
            GElonh(i,j) =                GElonh(i,j)   / hourad
C +...                     Conversion:   radian       -> Hour
C +
 31     continue
C +
       end if                                                          ! CTR
C +
C +- Oblique Stereographic Projection 
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if (maptyp.eq.1) then                                           ! CTR
C +
        do 32 j=1,my
        do 32 i=1,mx
         xxmar = (i-imez)*dx
         yymar = (j-jmez)*dy
C +
C +      ***********
         call GRDstr(xxmar,yymar,GElon0,GElat0,GElon,GElat)
C +      ***********
C +
         GElatr(i,j) =  GElat
         GElonh(i,j) =  GElon
C +
 32     continue
C +
C +- Lambert Comformal Projection (2 std parallels)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       else if (maptyp.eq.2) then                                      ! CTR
C +
C +     ***********
        call GRDlam
C +     ***********
C +
       end if                                                          ! CTR

C +
C +
C +--Coriolis Parameter
C +  ==================
C +
      do 33 j=1,my
      do 33 i=1,mx
      fcorDY(i,j) = 2.d0*earthv*sin(GElatr(i,j))
C +...fcorDY      : Coriolis Parameter
C +
 33   continue
C +
      end if                                                           ! CTR
C +
C +
C +--Time Zone
C +  =========
C +
      do 300 j=1,my
      do 300 i=1,mx
          itizGE(i,j)    =    GElonh(i,j)
      if (itizGE(i,j).gt. 12) itizGE(i,j) = itizGE(i,j)-24
      if (itizGE(i,j).lt.-12) itizGE(i,j) = itizGE(i,j)+24
 300  continue
C +
C +
C +--OUTPUT
C +  ======
C +
          i1   = imez  - 50
          i2   = imez  + 50
          j1   = jmez  - 50
          j2   = jmez  + 50
          i1   = max(i1, 1)
          i2   = min(i2,mx)
          j1   = max(j1, 1)
          j2   = min(j2,my)
          id10 = 1 + min(mx-1,10)
          jd10 = 1 + min(my-1,10)
C +
c         write(4,994)(i,i=i1,i2,id10)
c994      format(/,' LATITUDES / LONGITUDES / TOPOGRAPHY:  x ->  y ^ ',
c    .           /,' ===================================',/,9x,13i9)
c         do 999 j=j2,j1,-jd10
c         write(4,995)j,(GElatr(i,j)/degrad,i=i1,i2,id10)
c995      format(  i9,11f9.3)
c         write(4,996)  (GElonh(i,j)*1.5d+1,i=i1,i2,id10)
c996      format(  9x,11f9.3)
c         write(4,997)  (sh    (i,j)*1.0d-3,i=i1,i2,id10)
c997      format(  9x,11f9.3)
c         write(4,998)  (itizGE(i,j),       i=i1,i2,id10)
c998      format(  9x,11i9  )
c999      continue
C +
      return
      end
C +
      subroutine GRDstr(xxmar,yymar,GElon0,GElat0,GElon,GElat)
C +
C +------------------------------------------------------------------------+
C | MAR GRID                                                8-03-1996  MAR |
C |   SubRoutine GRDstr computes the Latitudes, Longitudes                 |
C |                              of a MAR Domain Grid Point                |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   METHOD: Inverse Stereographic Oblique Projection                     |
C |   ^^^^^^^                                                              |
C |                                                                        |
C |   REFERENCE: F. Pearson, Map projection methods, CRC Press, 1990.      |
C |   ^^^^^^^^^^                                                           |
C |                                                                        |
C |   INPUT:  xxmar,yymar  : MAR        Coordinates                        |
C |   ^^^^^^  GElon0,GElat0: Geographic Coordinates of MAR Domain Center   |
C |                          (3-D): South-North Direction along            |
C |                                 90E, 180E, 270E or 360E Meridians      |
C |                                                                        |
C |   OUTPUT: GElat        : Latitude  of the MAR grid point      (radian) |
C |   ^^^^^^^ GElon        : Longitude of the MAR grid point        (hour) |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
C +--General Variables
C +  =================
C +
C #DP include 'MAR_r8.inc'
C +
      include 'MARphy.inc'
C +
      include 'MARdim.inc'
      include 'MARgrd.inc'
C +
C +--local  Parameters
C +  =================
C +
      pidemi= pi / 2.0d0
C +
      CphiP = cos(degrad*GElat0)
      SphiP = sin(degrad*GElat0)
C +
C +
C +--Coordinates relative to a Pole set to the Domain Center
C +  =======================================================
C +
C +
C +--Relative Longitude -OBLlon (0 <= OBLlon < 2pi)
C +  ----------------------------------------------
C +
      if      (xxmar.gt.0.)                      then
         OBLlon = pidemi - atan(yymar/xxmar)
      else if (xxmar.eq.0. .and. yymar.lt.0.)    then
         OBLlon = pi
      else if (xxmar.lt.0.)                      then
         OBLlon = 3.0d0*pidemi - atan(yymar/xxmar)
      else if (xxmar.eq.0. .and. yymar.ge.0.)    then
         OBLlon = 0.d0
      end if
C +
C +
C +--Relative  Latitude  OBLlat
C +  --------------------------
C +
      ddista = sqrt ( xxmar*xxmar + yymar*yymar )
      OBLlat = 0.5d0*pi - 2.d0*atan(ddista*0.5d0/earthr)
C +
C +
C +--Coordinates Change (OBLlon,OBLlat) -> (GElon,GElat)
C +                   / (rotation, Pearson p.57)
C +  ===================================================
C +
C +
C +--Latitude (radians)
C +  ------------------
C +
      Sphi = SphiP * sin(OBLlat) + CphiP * cos(OBLlat) * cos(OBLlon)
      GElat= asin(Sphi)
C +
C +
C +--Longitude  (hours)
C +  ------------------
C +
C +--dGElon = GElon - GElon0  (-pi < dGElon <= pi)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      denomi =  CphiP * tan (OBLlat) - SphiP * cos(OBLlon)
C +
      if (OBLlon.gt.epsi          .and. OBLlon.lt.(pi-epsi))     then
C +
C +--1) OBLlon in trigonometric quadrant 1 or 4 ("right"):
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        dGElon = atan(sin(OBLlon)/denomi)
        if (dGElon.lt.0.d0) then
            dGElon = dGElon + pi
C +...      Go to Quadrant 1 by adding       180 degrees
        end if
C +
C +--2) OBLlon is in trigonometric quadrant 2or3 ("left "):
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      else if (OBLlon.gt.(pi+epsi).and. OBLlon.lt.(2.0*pi-epsi)) then
C +
        dGElon = atan(sin(OBLlon)/denomi)
        if (dGElon.gt.0.d0) then
            dGElon = dGElon - pi
C +...      Go to Quadrant 2 by substracting 180 degrees
        end if
C +
      else if (OBLlon.le.epsi .or. OBLlon.ge.(2.0*pi-epsi))      then
C +
C +--3) OBLlon = 0 -> dGElon = 0 or pi :
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ((pidemi-OBLlat) .gt. (pidemi-degrad*GElat0) ) then
C +...    North pole crossed ==> add 180 degrees to Longitude
          dGElon = pi
        else
          dGElon = 0.d0
        end if
C +
      else if (OBLlon.ge.(pi-epsi) .and. OBLlon.le.(pi+epsi))    then
C +
C +--4) OBLlon = pi -> dGElon = 0 or pi :
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ((pidemi-OBLlat) .gt. (pidemi+degrad*GElat0) ) then
C +...    South pole crossed ==> add 180 degrees to Longitude
          dGElon = pi
        else
          dGElon = 0.d0
        end if
      end if
C +
C +--Longitude (hours)
C +  ~~~~~~~~~
      GElon= (dGElon + GElon0 * degrad) / hourad
C +
      return
      end
C +
C +
      subroutine GRDlam
C +
C +------------------------------------------------------------------------+
C | MAR GRID                                               20-10-1997  MAR |
C |   SubRoutine GRDlam computes the Latitudes, Longitudes                 |
C |                              of a MAR Domain Grid Point                |
C |                     using Lambert projection                           |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   METHOD: Inverse Lambert Comformal Projection                         |
C |   ^^^^^^^       (conical, two standard parrallels)                     |
C |                                                                        |
C |   REFERENCE: F. Pearson, Map projection methods, CRC Press, 1990.      |
C |   ^^^^^^^^^^                                                           |
C |                                                                        |
C |   INPUT:                                                               |
C |   ^^^^^^  GElon0,GElat0: Geographic Coordinates of MAR Domain Center   |
C |                          (both in degree!)                             |
C |                                                                        |
C |   OUTPUT: GElatr(mx,my): Latitude  of MAR grid points         (radian) |
C |   ^^^^^^^ GElonh(mx,my): Longitude of MAR grid points           (hour) |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
C +
C +--General Variables
C +  =================
C +
c #DP include 'MAR_r8.inc'
C +
      include 'MARphy.inc'
C +
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
C +
      include 'MAR_DY.inc'

C +
C +--Local constants:
C +  ================
      rayter=6371229.0
      pi = acos (-1.)
      toDEG = 180./pi
C +
C +--Domain center (radiant) and size:
C +  ---------------------------------
      Rlat0 = degrad*GElat0
 
      RlatSz = float(my) * dy / earthr
 
C +--"True latitudes" phi1 and phi2 :
C +  --------------------------------
      phi1 = Rlat0 - RlatSz / 3.
      phi2 = Rlat0 + RlatSz / 3.
 
C +--Constants rK (sin phi0) and psi
C +  -------------------------------
      xx = cos (phi1) / cos (phi2)
      yy = tan (pi /4. - phi1 /2.) / tan (pi /4. - phi2 /2.)
      rK= log (xx) / log (yy)
      psi = rayter*cos(phi1) / (rK*(tan(pi/4.-phi1 /2.))**rK)

C +--y distance from center to pole
C +  ------------------------------
      delty = psi * (tan(pi/4. - Rlat0/2.))**rK
 
C +--Main loop over grid points.
C +  ===========================
      do j = 1, my     
      do i = 1, mx
  
C +--  Mar coordinate
C +    --------------
        xxmar = (i-imez)*dx
        yymar = (j-jmez)*dy

C +--  Transformation to pole-centered xP,yP
C +    -------------------------------------
        xxP = delty - yymar
        yyP = xxmar

C +--  Coordinate change : to polar.
C +    -----------------------------
        pol_r = SQRT (xxP**2. + yyP**2.)
        theta = ATAN (yyP/xxP)
 
C +--  Compute longitude (hour)
C +    ------------------------
        GElonh(i,j) = (GElon0 + theta/rK / degrad) / 15.

C +--  Compute latitude (radian)
C +    -------------------------
        GElatr(i,j) = (pi/2.)-2.*ATAN((pol_r/psi)**(1./rK))

      end do
      end do

      return
      end
      subroutine GRDsig(zmin,aavu,bbvu,ccvu,vertic)
C +
C +------------------------------------------------------------------------+
C | MAR GRID                                               12-06-1996  MAR |
C |   SubRoutine GRDsig is used to initialize the vertical grid            |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   ASSUMPTION: Sigma is calculated from initial level height amsl       |
C |   ^^^^^^^^^^^                     assumig that T(msl) = SST            |
C |                                                dT/dz  = -0.0065 K/m    |
C |                                                p_s    = 100     hPa    |
C |                                                                        |
C |   INPUT  : zmin           : Height above Surface / 1st Sigma Level (m) |
C |   ^^^^^^^^ aavu,bbvu,ccvu : Vertical Discretization Parameters         |
C |            vertic         : Logical Variable caracteris.vertic.discris.|
C |            sigpar(10)     : Parish Model Vertical Discretisation       |
C |                                                                        |
C |   OUTPUT : Variable  which is  initialized is:                         |
C |   ^^^^^^^^  sigma(mz): Independant Variable (Normalized Pressure)      |
C |                                                                        |
C |   OPTIONS: #PA  Parish Model Vertical Discretisation                   |
C |   ^^^^^^^^ #SA  Regular      Vertical Discretisation                   |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
C +--General Variables
C +  =================
C +
c #DP include 'MAR_r8.inc'
C +
      include 'MARphy.inc'
C +
      include 'MARdim.inc'
      include 'MARgrd.inc'
C +
      include 'MARSND.inc'
C +
      include 'MAR_SL.inc'
      include 'MAR_TU.inc'
C +
      include 'MAR_IO.inc'
C +
      logical vertic
C +
C +
C +--Local   Variables
C +  =================
C +
      dimension   zpbl(mz)
C +
c #HE dimension sighei(mz)
c #PA dimension sigpar(mz)
C +
C +
C +--DATA
C +  ====
C +
      data      ps_sig/101.3d0/
c #HE data      sighei/
c #HE.       /0.10015,0.19077,0.27276,0.34695,0.41409,0.47483,0.52979,
c #HE.0.57952,0.62452,0.66524,0.70208,0.73542,0.76558,0.79288,0.81757,
c #HE.0.83992,0.86014,0.87844,0.89499,0.90997,0.92352,0.93579,0.94688,
c #HE.0.95692,0.96601,0.97423,0.98167,0.98840,0.99111/
C +***          sighei: DNMI model Vertical Discretisat.  (Heinemann 1996)
C +
c #PA data      sigpar/0.100,0.350,0.600,0.800,0.900,
c #PA.                 0.930,0.950,0.970,0.985,0.996/
C +***          sigpar: Vertical Discretisation of Parish Model
C +                     (Bromwich, Du and Parish 1994 MWR 122 No 7 p.1418)
      data lsf / 1/
C +
      data ga0 / 0.0065d0/
C +...     ga0 : Standard Atmospheric Lapse Rate
C +
      lvg = 0
C +...lvg : set to 1 if |Vg(sounding)| .ne. 0 anywhere
C +
C +
C +--Entry Checking Point
C +  ====================
C +
      if (IO_loc.ge.2)
     . open(unit=21,status='unknown',file='Out/GRID.info')

      if (IO_loc.ge.2) write(21,999)
 999  format(//,'   --- Initialisation / GRDsig ---')
C +
C +
C +--Temperature Vertical Profile
C +  ============================
C +
      ga  = ga0
C +
      if (IO_loc.ge.2) write(21,1)ga,sst_SL,ps_sig,gravit,ra
 1    format(/,'  dT/dz  =',f8.5,' K/m',
     .       /,'  SST    =',f8.2,' K',
     .       /,'  ps_sig =',f8.2,' kPa',
     .       /,'  gravit =',f8.2,' m/s2',
     .       /,'  ra     =',f8.2,' J/kg/K')
C +
C +
C +--Sigma Levels
C +  ============
C +
C +- 1) Coarse Resolution of the Surface Layer
C +  -----------------------------------------
C +
      if (.not.vertic) then
C +
C +    aa = 0.5
C +    bb = 1.5
C +    cc =-1.0
C +... Reference : E. Richard, these, 1991, p.29
C +
       vu =       0.0d0
       do 11 k=1,mz
       vu = vu  + 1.0d0/dble(mzz)
       sigma(k) = aavu*vu + bbvu*vu*vu + ccvu*vu*vu*vu
C +
c #HE  sigma(k) = sighei(k)
c #PA  sigma(k) = sigpar(k)
C +... Vertical Discretisation of Parish Model is used
C +
       if (abs(ga).gt.1.d-5) then
        zpbl(k) =-(   sst_SL  /ga) *   ((1.d0+(sigma(k)-1.d0)
     .                                 *(1.d2/ps_sig))
     .                                 **(ra*ga/gravit)-1.d0)
       else
        if (IO_loc.ge.2.and.k.eq.1) write(21,110)
 110    format(/,'  t(z)   = CONSTANT')
        zpbl(k) =-(ra*sst_SL  /gravit ) *log((unun+(sigma(k)-unun)
     .                                 *(1.d2/ps_sig)))
       end if
 11    continue
C +
C +
C +- 2) Fine   Resolution of the Surface Layer
C +  -----------------------------------------
C +
      else
C +
       ga     =max(ga,epsi)
C +
       zpbl(1)=     zmin
       zpbl(2)=2.d0*zmin
C +... Last option (#SA) is used when testing radiative routines.
C +
       dzz    =0.d0
       do 120 k=3,mz
       rz     =zmin*aavu **(k-1)
       rzb    =ccvu*bbvu **(k-1)
       if (TUkhmx.gt.0.d0) then
       zpbl(k)=rzb   *rz /(rz + rzb  )
       else
       zpbl(k)=       rz
       end if
C +
       zzo    = zpbl(k)
       zpbl(k)= max(zpbl(k),zpbl(k-1)+zpbl(2))
       dzz    = max(zpbl(k)-zzo,      zero   ) + dzz
 120   continue
       do 12 k=1,mz
       kk=mz+1-k
       sigma(kk)=(ps_sig/100.d0)
     .         *((1.0d0-ga*zpbl(k)/sst_SL)**(gravit/(ga*ra))-1.d0)+1.d0
C +... sigma(kk): the fine resolution of the surface layer is computed
C +           using a geometric progression
C +
 12    continue
      end if
C +
C +
C +--Output
C +  ======
C +
      if (IO_loc.ge.2) then
       write(21,130)(sigma(k)   ,k=1,mz   )
 130   format(/,'  Sigma    Levels :',/,(1x,15f8.4))
       write(21,131)(zpbl(k),k=mz,1,-1)
 131   format(/,'  Altitude Levels :',/,(1x,15f8.1))
      end if
C +
      if (IO_loc.ge.2)
     . close(unit=21)

      return
      end


      subroutine DYNgpo
C +
C +------------------------------------------------------------------------+
C | MAR DYNAMICS   FAST                                    11-03-1996  MAR |
C |   SubRoutine DYNgpo contains the Integration of Hydrostatic Relation   |
C |                          and the Computation of Real Temperature t [K] |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |    INPUT: pktaDY(mx,my,mzz): Reduced Potential Temperature             |
C |    ^^^^^^ gplvDY(mx,my,mzz): Surface Geopotential (i.e. for k=mzz)     |
C |            virDY(mx,my,mz) : Air Loading by water vapor & hydrometeors |
C |                                                                        |
C |   OUTPUT:   pkDY(mx,my,mz) : Exner Potential                           |
C |   ^^^^^^^ tairDY(mx,my,mz) : Temperature                           [K] |
C |           gplvDY(mx,my,mzz):    Geopotential                   [m2/s2] |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
C +--General Variables
C +  =================
C +
c #DP include 'MAR_r8.inc'
C +
      include 'MARphy.inc'
C +
      include 'MARdim.inc'
      include 'MARgrd.inc'
C +
      include 'MAR_DY.inc'
C +
      include 'LSM_WK.inc'
C +
C +
C +--EXNER Potential and Temperature
C +  ===============================
C +
      sigmid(1)   = 0.d0
      sigmid(mzz) = 1.d0
      do k=2,mz
       sigmid(k)  =( sigma(k)   +  sigma(k-1)) /2.d0
      enddo
C +
      ab=0.5d0*(1.d0-sigmid(mz))/(1.d0-sigmid(mz1))
C +
      do 1 k=1,mz
      do 1 j=1,my
      do 1 i=1,mx
         pkDY(i,j,k) = exp(cap *log(pstDYn(i,j)*sigma(k)+ptopDY))
       tairDY(i,j,k) = pktaDY(i,j,k) *pkDY(i,j,k)
 1    continue
C +
      do 11 k=1,mz
      do 11 j=1,my
      do 11 i=1,mx
       WKxyz1(i,j,k) = cp * pkDY(i,j,k)
 11   continue
C +
      do 12 k=1,mz1
      do 12 j=1,my
      do 12 i=1,mx
       WKxyz2(i,j,k) = cp * pkDY(i,j,k+1)
 12   continue
C +
             k=  mz
      do 120 j=1,my
      do 120 i=1,mx
       WKxyz2(i,j,k) = cp *
     .                 exp(cap *log(pstDYn(i,j)         +ptopDY))
 120  continue
C +
C +
C +--Integration of the Hydrostatic Equation
C +  =======================================
C +
      do 20 j=1,my
      do 20 i=1,mx
       gplvDY(i,j,mz)=gplvDY(i,j,mzz)+(WKxyz2(i,j,mz)-WKxyz1(i,j,mz))
     .    *((1.d0+ab)*pktaDY(i,j,mz) *(1.d0+virDY(i,j,mz ))
     .           -ab *pktaDY(i,j,mz1)*(1.d0+virDY(i,j,mz1)))
 20   continue
C +
      do 2  k=mz1,1,-1
      do 2  j=1,my
      do 2  i=1,mx
       gplvDY(i,j,k) =gplvDY(i,j,k+1)+(WKxyz2(i,j,k )-WKxyz1(i,j,k ))
     .              *(pktaDY(i,j,k)  *(1.d0+virDY(i,j,k))
     .               +pktaDY(i,j,k+1)*(1.d0+virDY(i,j,k+1)))*0.5d0
C +
 2    continue
C +
C +
C +--Work Arrays Reset
C +  =================
C +
      do 30 k=1,mz
      do 30 j=1,my
      do 30 i=1,mx
      WKxyz1(i,j,k) = 0.d0
      WKxyz2(i,j,k) = 0.d0
 30   continue
C +
      return
      end
