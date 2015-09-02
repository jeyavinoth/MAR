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
C +
      include 'LSMphy.inc'
      include 'NSTdim.inc'
C +
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
