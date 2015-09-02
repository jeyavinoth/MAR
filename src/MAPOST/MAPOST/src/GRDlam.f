      subroutine GRDlam(xxmar,yymar,dy,GElon0,GElat0,GElon,GElat)
C +
C +------------------------------------------------------------------------+
C | MAR GRID                                               20-10-1997  MAR |
C |   SubRoutine GRDlam computes the Latitudes, Longitudes                 |
C |                              of a MAR Domain Grid Point                |
C |                     using Lambert projection                           |
C |                                                                        |
C |   WARNING: modified version for MAPOST (only 1 grid pt = more simple)  |
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
C |   OUTPUT: GElat(mx,my): Latitude  of MAR grid points          (radian) |
C |   ^^^^^^^ GElon(mx,my): Longitude of MAR grid points            (hour) |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
C +
C +--General Variables
C +  =================
C +
C +...* MAR dimensions :
      include 'NSTdim.inc'

      include 'LSMphy.inc'
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
      psi = earthr*cos(phi1) / (rK*(tan(pi/4.-phi1 /2.))**rK)

C +--y distance from center to pole
C +  ------------------------------
      delty = psi * (tan(pi/4. - Rlat0/2.))**rK
 
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
        GElon = (GElon0 + theta/rK / degrad) / 15.

C +--  Compute latitude (radian)
C +    -------------------------
        GElat = (pi/2.)-2.*ATAN((pol_r/psi)**(1./rK))

      return
      end
