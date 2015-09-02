C   +-------------------------------------------------------------------+
C   +  Subroutine MARhgd                            19-06-2004  NESTING +
C   +-------------------------------------------------------------------+
C   +                                                                   +
C   + Input : Parameters from MARgrd.ctr                                +
C   + ^^^^^^^                                                           +
C   +                                                                   +
C   + Output: Creation of the horizontal grid of MAR                    +
C   + ^^^^^^^ Variables : NST__x(mx,my) and NST__y(mx,my)  (long./lat.) +
C   +                     NSTgdx(mx)    and NSTgdy(my)     (Lambert)    +
C   +                     NST_dx (horizontal resolution)                +
C   +                                                                   +
C   +-------------------------------------------------------------------+


      SUBROUTINE MARhgd


      IMPLICIT NONE


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'MARvar.inc'

C +---Local variables
C +   ---------------

      INTEGER i,j

      REAL degrad,MinLon,MaxLon,MinLat,MaxLat,DEGresol,argrot
                         
C +---Constants
C +   ---------

      DATA degrad /  1.745329252d-2/
C +...     degrad : Conversion Factor: Radian --> Degrees


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---READING OF GRID PARAMETERS IN MARgrd.ctr
C +   ========================================

      OPEN (unit=51,status='old',file='MARgrd.ctr')

       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) maptyp
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) GElon0
       read (51,*) imez  
       read (51,*) GElat0
       read (51,*) jmez  
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) dx
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) GEddxx
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) ptopDY
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) zmin
       read (51,*) aavu
       read (51,*) bbvu
       read (51,*) ccvu
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,'(l4)') vertic
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) sst_SL
       read (51,*) !- - - - - - - - - - - - - - - - - -

      CLOSE(unit=51)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---HORIZONTAL RESOLUTION
C +   =====================

      dx = dx * 1000.
      dy = dx

      NST_dx=dx


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CREATION OF HORIZONTAL 3-D MAR GRID
C +   ===================================


      IF (NSTmod.eq.'MAR') THEN


C +---Choice of map projection
C +   ------------------------

       write(6,'(A,$)') ' Map projection        :'

       IF (maptyp.eq.0)
     .  write(6,*) ' Polar Stereographic Projection (Antarctica)'
       IF (maptyp.eq.1)
     .  write(6,*) ' Stereographic Oblique Projection'
       IF (maptyp.eq.2)
     .  write(6,*) ' Lambert Conformal, 2 Std. Par. Projection'
       write(6,*)


C +---Domain reference grid point
C +   ---------------------------

       IF (imez.le.0.or.imez.gt.mx) imez = mx/2
       IF (jmez.le.0.or.jmez.gt.my) jmez = my/2


C +---Simple grid (Lambert coordinates)
C +   ---------------------------------

       DO i=1,mx
        NSTgdx(i)=(i-imez)*dx/1000.
       ENDDO
 
       DO j=1,my
        NSTgdy(j)=(j-jmez)*dy/1000.
       ENDDO


C +---Compute map projection
C +   ----------------------

C +         ******
       CALL GRDgeo (maptyp,imez,jmez,dx,dy,GElon0,GElat0,
     .              GEddxx,NST__x,NST__y)
C +         ******


C +---Convertion to degree units
C +   --------------------------

       DO j=1,my
       DO i=1,mx
        NST__x(i,j) =  NST__x(i,j) * 15.d0
C +...  Conversion: Hour->degrees
        NST__y(i,j) =  NST__y(i,j) / degrad
C +...  Conversion: rad ->degrees
       ENDDO
       ENDDO


      ENDIF  !  {NSTmod.eq.'MAR'}


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CREATION OF HORIZONTAL 2-D MAR GRID
C +   ===================================


      IF (NSTmod.eq.'M2D') THEN


C +---Simple grid (Lambert coordinates)
C +   ---------------------------------

       DO i=1,mx
        NSTgdx(i)=(i-imez)*dx/1000.
       ENDDO
 
       DO j=1,my
        NSTgdy(j)=(j-jmez)*dy/1000.
       ENDDO


C +---Compute map projection
C +   ----------------------

       DEGresol = 111.111111 * ABS(COS(GElat0*degrad))

       argrot = (90.d0-GEddxx)*degrad

       DO j=1,my
       DO i=1,mx
        NST__x(i,j) =  GElon0 
     .              + (NSTgdx(i) / DEGresol * COS(argrot))
     .              + (NSTgdy(j) / DEGresol * SIN(argrot))
        NST__y(i,j) =  GElat0
     .              + (NSTgdx(i) / DEGresol * SIN(argrot))
     .              + (NSTgdy(j) / DEGresol * COS(argrot))
C +...  Conversion: km -> degrees
       ENDDO
       ENDDO


      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Compute horizontal extent of the horizontal domain
C +   ==================================================

      MinLon = NST__x(1,1)
      MaxLon = NST__x(1,1)
      MinLat = NST__y(1,1)
      MaxLat = NST__y(1,1)
      DO j=1,my
      DO i=1,mx
       MinLon = MIN(NST__x(i,j),MinLon)
       MaxLon = MAX(NST__x(i,j),MaxLon)
       MinLat = MIN(NST__y(i,j),MinLat)
       MaxLat = MAX(NST__y(i,j),MaxLat)
      ENDDO
      ENDDO


C +---Print the characteristics of the grid
C +   =====================================

      write(6,200) mx,my,dx/1000.,GEddxx,MinLon,MaxLon,
     .             MinLat,MaxLat
200   format(' Grid points           : ',i4,' * ',i4,/,
     .       ' Horizontal resolution : ',f7.0,' km.',/,
     .       ' Domain orientation    : ',f7.0,' deg.',/,
     .       ' MAR longitude between : ',f7.2,' and ',f7.2,/,
     .       ' MAR latitude  between : ',f7.2,' and ',f7.2,/)

       write(6,300) mz,ptopDY
300    format(' Number of grid points : ',i4,/,
     .        ' Pressure at the top   : ',f9.4,' kPa.')
       write(6,310) zmin, aavu, bbvu, ccvu
310    format(' First level height    : ', f6.1,/,
     .        ' aavu, bbvu, ccvu      : ',(f6.1,', ',f6.1,', ',f6.1),/)


C +---nvx = mw ?
C +   ==========

      If(nvx.ne.mw)then
       write(6,201) nvx,mw
201    format(' WARNING -- nvx(',i1,') ne mw(',i1,') -- WARNING',/) 
      endif

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      subroutine GRDgeo(maptyp,imez,jmez,dx,dy,GElon0,GElat0,
     .                  GEddxx,GElonh,GElatr)
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
      implicit none
C +
C +
C +--General Variables
C +  =================
C +
      include 'NSTdim.inc'
C +
      integer i,j,mmx,mmy,maptyp,imez,jmez
C +
      real degrad,pi,argrot,cosrot,sinrot,xxmar,yymar,epsi,
     .     ddista,xdista,ydista,zero,dx,dy,earthr,hourad

      real GEddxx,GElonh(mx,my),GElatr(mx,my),GElon0,GElat0,
     .     GElon,GElat
C +
C +
C +--Some initialisations
C +  --------------------
C +
      mmx=mx
      mmy=my
C +
      pi    = 3.141592653589793238462643d0
      degrad= pi / 180.d0
      hourad= pi / 12.d0
      epsi  = 1.0d-6
      zero  = 0.d0
      earthr= 6371.229d+3
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
C +
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
C +- EPSG Polar Stereographic transformation Variant B 
C +-    (http://www.epsg.org/guides/docs/G7-2.pdf)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +
       if (maptyp.eq.0) then                                           ! CTR
         do j=1,my
         do i=1,mx
C +
           xxmar  = (i-imez)*dx/1000.
           yymar  = (j-jmez)*dy/1000.
C +
C +        ***********
           call StereoSouth(xxmar,yymar,GEddxx,GElon,GElat)
C +        ***********
C +
           GElonh(i,j) =  GElon / 15.
C +...     Conversion: degrees->hour
           GElatr(i,j) =  GElat * degrad
C +...     Conversion: rad ->degrees
C +
         enddo
         enddo
C +
       end if
C +
C C +
C +- Oblique Stereographic Projection 
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if (maptyp.eq.1) then                                           ! CTR
C +
        do 32 j=1,my
        do 32 i=1,mx
C +
         argrot = (GEddxx-90.d0)*degrad
         cosrot = cos(argrot)
         sinrot = sin(argrot)
         xxmar  = cosrot*(i-imez)*dx+sinrot*(j-jmez)*dy
         yymar  = cosrot*(j-jmez)*dy-sinrot*(i-imez)*dx
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
        call GRDlam (dx,dy,imez,jmez,GElon0,GElat0,GElatr,GElonh)
C +     ***********
C +
       end if                                                          ! CTR
C +
      end if                                                           ! CTR
C +
      return
      end

! +----------------------------------------------------------------------+
      Subroutine StereoSouth (E,N,GEddxx,lon,lat)
! |  Compute the lon, lat from Oblique Stereographic Projection          |
! |  Written by Cecile Agosta                                 17-05-10   |
! |  EPSG Polar Stereographic transformation Variant B                   |
! |  (http://www.epsg.org/guides/docs/G7-2.pdf)                          |
! |  Equivalent to EPSG 3031 (WGS-84 ellipsoid)                          |
! +----------------------------------------------------------------------+
! |                                                                      |
! | INPUT :  E      : Stereo coordinate on the East  (X, km)             |
! | ^^^^^^^  N      : Stereo coordinate on the North (Y, km)             |
! |          GEddxx : Longitude of X axis (=GEddxx, 90 = East, clockwise)|
! |         [lat true = 71 S]                                            |
! |                                                                      |
! | OUTPUT : lon    : longitude (deg)                                    |
! | ^^^^^^^  lat    : latitude  (deg)                                    |
! |                                                                      |
! +----------------------------------------------------------------------+
      IMPLICIT NONE

      INCLUDE 'NSTdim.inc'

! +-- General Variables
! +   -----------------
      Real,INTENT(in ) :: E,N,GEddxx
      Real,INTENT(out) :: lon,lat

! +-- Local Variables
! +   ---------------
      Real ddista

! +-- Constants
! +   ---------
      Real  aa,ex,pi,degrad,latF,FE,FN,tF,mF,k0,t,rho,khi,lon0
      Real  trulat

      aa     = 6378.1370         ! aa (km) = demi grand axe
      ex     = 0.081819190842621 ! excentricity WGS-84 : 0.081 819 190 842 622 0.081 819 190 842 621
      trulat = -71.              ! Latitude of standard parallel, 71 S for ESPG 3031
      pi     = 4. * atan(1.)
      degrad = pi / 180.

      latF = trulat*degrad
      lon0 = (GEddxx-90.)*degrad
      
      FE = 0. !False Easting
      FN = 0. !False Northing

! +-  Polar Stereographic Projection
! +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! +----------------------------------------------------------------------+
      tF = tan (pi/4 + latF/2) /
     .    ( (1 + ex*sin(latF)) / (1 - ex*sin(latF)) )**(ex/2)
      mF = cos(latF) / (1 - ex**2 * sin(latF)**2)**0.5
      k0 = mF   *( (1+ex)**(1+ex) * (1-ex)**(1-ex) )**0.5 / (2*tF)

      rho = ( (E-FE)**2 + (N-FN)**2 )**0.5
      t   = rho * ( (1+ex)**(1+ex) * (1-ex)**(1-ex) )**0.5 / (2*aa*k0)
      khi = 2*atan(t) - pi/2

      lat = khi 
     .      + (  ex**2/2   +  5*ex**4/24  +     ex**6/12 + 13*ex**8/360)
     .        *sin(2*khi)
     .      + (7*ex**4/48  + 29*ex**6/240 + 811*ex**8/11520            )
     .        *sin(4*khi)  
     .      + (7*ex**6/120 + 81*ex**8/1120                             )
     .        *sin(6*khi)  
     .      + (            4279*ex**8/161280                           )
     .        *sin(8*khi)

      if      (E-FE .eq. 0. .and. N-FN .ge. 0) then
        lon = lon0
      else if (E-FE .eq. 0. .and. N-FN .le. 0) then
        lon = lon0 + pi
      else
        lon = lon0 + atan2(E-FE,N-FN)
      endif

      lat = lat / degrad
      lon = lon / degrad

      Return
      End Subroutine StereoSouth
C +------------------------------------------------------------------------+

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
      implicit none
C +
C +
C +--General Variables
C +  =================
C +
      include 'NSTdim.inc'
C +
      real pidemi,pi,CphiP,SphiP,degrad,OBLlon,xxmar,yymar,Sphi,
     .     denomi,dGElon,epsi,OBLlat,ddista,earthr,hourad
C +
      real GElon0,GElat0,GElon,GElat
C +
C +
C +--local  Parameters
C +  =================
C +
      pi    = 3.141592653589793238462643d0
      degrad= pi / 180.d0
      hourad= pi / 12.d0
      epsi  = 1.0d-6
      earthr= 6371.229d+3
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


      subroutine GRDlam (dx,dy,imez,jmez,GElon0,GElat0,GElatr,GElonh)
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
      implicit none
C +
C +
C +--General Variables
C +  =================
C +
      include 'NSTdim.inc'
C +
      integer i,j,imez,jmez
C +
      real rayter,pi,Rlat0,RlatSz,phi1,phi2,xx,yy,rK,psi,
     .     delty,xxmar,yymar,xxP,yyP,pol_r,theta,GElon0,GElat0,
     .     GElonh(mx,my),GElatr(mx,my),dx,dy,earthr,degrad
C +
C +
C +--Local constants:
C +  ----------------
      rayter=6371229.0
      pi    = 3.141592653589793238462643d0
      degrad= pi / 180.d0
      earthr= 6371.229d+3
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
