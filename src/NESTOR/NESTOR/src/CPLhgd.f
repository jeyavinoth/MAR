C   +-------------------------------------------------------------------+
C   +  Subroutine CPLhgd                          January 2002  NESTING +
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


      SUBROUTINE CPLhgd


      IMPLICIT NONE


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'

C +---Local variables
C +   ---------------

      INTEGER i,j,i1,i2,i3,ilat_min,ilat_max,ilon_min,
     .        ilon_max,mmx,mmy,imez,jmez

      REAL degrad,MinLon,MaxLon,MinLat,MaxLat,long1,lati1,
     .     long2,lati2,long3,lati3,long4,lati4

      REAL long_tmp(4),lati_tmp(4),lati_min,lati_max,
     .     long_min,long_max,auxlong,auxlati,vlon(4),vlat(4),
     .     dx,dy,delta_lon,delta_lat,center_lon,center_lat


C +---Constants
C +   ---------

      DATA degrad /  1.745329252d-2/
C +...     degrad : Conversion Factor: Radian --> Degrees


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Check validity of array dimensions
C +   ==================================

      mmx = mx
      mmy = my

      IF (mmx.ne.3.or.mmy.ne.3) THEN
       WRITE(6,*) ' '
       WRITE(6,*) 'CPL (SVAT coupling) is valid only when mx = 3 and'
       WRITE(6,*) 'my = 3 in NSTdim.inc.'
       WRITE(6,*) 'Please modify these parameters and rerun NESTOR.'
       WRITE(6,*) ' '
       WRITE(6,*) 'STOP in CPLhgd.f'
       WRITE(6,*) ' '
       STOP
      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---READING OF GRID PARAMETERS IN CPLgrd.ctr
C +   ========================================

      OPEN (unit=51,status='old',file='CPLgrd.ctr')

       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) long1
       read (51,*) lati1
       read (51,*) long2
       read (51,*) lati2
       read (51,*) long3
       read (51,*) lati3
       read (51,*) long4
       read (51,*) lati4
ccccc  read (51,*) !- - - - - - - - - - - - - - - - - -
ccccc  read (51,*) ptopDY
ccccc  read (51,*) !- - - - - - - - - - - - - - - - - -
ccccc  read (51,*) zmin
ccccc  read (51,*) aavu
ccccc  read (51,*) bbvu
ccccc  read (51,*) ccvu
ccccc  read (51,*) !- - - - - - - - - - - - - - - - - -
ccccc  read (51,'(l4)') vertic
ccccc  read (51,*) !- - - - - - - - - - - - - - - - - -
ccccc  read (51,*) sst_SL
ccccc  read (51,*) !- - - - - - - - - - - - - - - - - -

      CLOSE(unit=51)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Ordering input grid points
C +   --------------------------

C +...Temporary arrays

      long_tmp(1) = long1
      long_tmp(2) = long2
      long_tmp(3) = long3
      long_tmp(4) = long4
      lati_tmp(1) = lati1
      lati_tmp(2) = lati2
      lati_tmp(3) = lati3
      lati_tmp(4) = lati4

C +...Search for minimum and maximum latitudes

      lati_min = lati_tmp(1)
      ilat_min = 1
      lati_max = lati_tmp(1)
      ilat_max = 1
      DO i=2,4
       IF (lati_tmp(i).lt.lati_min) THEN
        lati_min = lati_tmp(i)
        ilat_min = i
       ENDIF
       IF (lati_tmp(i).gt.lati_max) THEN
        lati_max = lati_tmp(i)
        ilat_max = i
       ENDIF
      ENDDO

C +...Grid point of minimum latitude becomes the first grid point

      auxlong = long_tmp(1)
      auxlati = lati_tmp(1)
      long_tmp(1) = long_tmp(ilat_min)
      lati_tmp(1) = lati_tmp(ilat_min)
      long_tmp(ilat_min) = auxlong
      lati_tmp(ilat_min) = auxlati 

C +...Grid point of maximum latitude becomes the fourth grid point

      auxlong = long_tmp(4)
      auxlati = lati_tmp(4)
      long_tmp(4) = long_tmp(ilat_max)
      lati_tmp(4) = lati_tmp(ilat_max)
      long_tmp(ilat_max) = auxlong
      lati_tmp(ilat_max) = auxlati

C +...Search for the next minimum latitude

      lati_min = lati_tmp(2)
      ilat_min = 2
      DO i=3,4
       IF (lati_tmp(i).lt.lati_min) THEN
        lati_min = lati_tmp(i)
        ilat_min = i
       ENDIF
      ENDDO

C +...Grid point of next minimum latitude becomes grid point 2

      auxlong = long_tmp(2)
      auxlati = lati_tmp(2)
      long_tmp(2) = long_tmp(ilat_min)
      lati_tmp(2) = lati_tmp(ilat_min)
      long_tmp(ilat_min) = auxlong
      lati_tmp(ilat_min) = auxlati 

C +...Exchange grid points 1 and 2 if necessary (depending on their long.)

      IF (long_tmp(2).lt.long_tmp(1)) THEN
       auxlong     = long_tmp(2)
       auxlati     = lati_tmp(2)
       long_tmp(2) = long_tmp(1)
       lati_tmp(2) = lati_tmp(1)
       long_tmp(1) = auxlong
       lati_tmp(1) = auxlati 
      ENDIF

C +...Search for the next maximum latitude

      lati_max = lati_tmp(1)
      ilat_max = 1
      DO i=2,3
       IF (lati_tmp(i).gt.lati_max) THEN
        lati_max = lati_tmp(i)
        ilat_max = i
       ENDIF
      ENDDO

C +...Grid point of next maximum latitude becomes grid point 3

      auxlong = long_tmp(3)
      auxlati = lati_tmp(3)
      long_tmp(3) = long_tmp(ilat_max)
      lati_tmp(3) = lati_tmp(ilat_max)
      long_tmp(ilat_max) = auxlong
      lati_tmp(ilat_max) = auxlati 

C +...Exchange grid points 3 and 4 if necessary (depending on their long.)

      IF (long_tmp(3).lt.long_tmp(4)) THEN
       auxlong     = long_tmp(4)
       auxlati     = lati_tmp(4)
       long_tmp(4) = long_tmp(3)
       lati_tmp(4) = lati_tmp(3)
       long_tmp(3) = auxlong
       lati_tmp(3) = auxlati 
      ENDIF

C +...Store ordered lat/long

      DO i=1,4
       vlon(i) = long_tmp(i)
       vlat(i) = lati_tmp(i)
      ENDDO

C +...Center of the cell

      center_lon = 0.25*(vlon(1)+vlon(2)+vlon(3)+vlon(4))
      center_lat = 0.25*(vlat(1)+vlat(2)+vlat(3)+vlat(4))


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---HORIZONTAL RESOLUTION
C +   =====================

      delta_lon = ABS(vlon(2)-vlon(1)+vlon(3)-vlon(4))
      delta_lat = ABS(vlat(4)-vlat(1)+vlat(3)-vlat(2))

      dx = 0.5*delta_lon*111111.0
      dy = 0.5*delta_lat*111111.0 * ABS(COS(vlat(1)/degrad))

      NST_dx=MAX(dx,dy)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CREATION OF HORIZONTAL MAR GRID
C +   ===============================


C +---Domain reference grid point
C +   ---------------------------

      IF (imez.lt.0.or.imez.gt.mx) imez = 2
      IF (jmez.lt.0.or.jmez.gt.my) jmez = 2


C +---Simple grid (Lambert coordinates)
C +   ---------------------------------

      DO i=1,mx
       NSTgdx(i)=(i-imez)*dx/1000.
      ENDDO

      DO j=1,my
       NSTgdy(j)=(j-jmez)*dy/1000.
      ENDDO


C +---Create NST grid
C +   ---------------

      i1 = 1
      i2 = 2
      i3 = 3

      NST__x(i2,i2) = center_lon
      NST__y(i2,i2) = center_lat

      NST__x(i1,i1) = NST__x(i2,i2) - 2.0*ABS(center_lon-vlon(1))
      NST__y(i1,i1) = NST__y(i2,i2) - 2.0*ABS(center_lat-vlat(1))
      NST__x(i3,i1) = NST__x(i2,i2) + 2.0*ABS(center_lon-vlon(2))
      NST__y(i3,i1) = NST__y(i2,i2) - 2.0*ABS(center_lat-vlat(2))
      NST__x(i1,i3) = NST__x(i2,i2) - 2.0*ABS(center_lon-vlon(4))
      NST__y(i1,i3) = NST__y(i2,i2) + 2.0*ABS(center_lat-vlat(4))
      NST__x(i3,i3) = NST__x(i2,i2) + 2.0*ABS(center_lon-vlon(3))
      NST__y(i3,i3) = NST__y(i2,i2) + 2.0*ABS(center_lat-vlat(3))
      NST__x(i2,i1) = 0.5*(vlon(1)+vlon(2))
      NST__y(i2,i1) = NST__y(i2,i2) 
     .              - 1.0*ABS(2.0*center_lat-vlat(1)-vlat(2))
      NST__x(i2,i3) = 0.5*(vlon(3)+vlon(4))
      NST__y(i2,i3) = NST__y(i2,i2) 
     .              + 1.0*ABS(2.0*center_lat-vlat(3)-vlat(4))
      NST__x(i1,i2) = NST__x(i2,i2) 
     .              - 1.0*ABS(2.0*center_lon-vlon(1)-vlon(4))
      NST__y(i1,i2) = 0.5*(vlat(1)+vlat(4))
      NST__x(i3,i2) = NST__x(i2,i2) 
     .              + 1.0*ABS(2.0*center_lon-vlon(2)-vlon(3))
      NST__y(i3,i2) = 0.5*(vlat(2)+vlat(3))


C +---Compute horizontal extent of the horizontal domain
C +   --------------------------------------------------

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


C +---Print the characteristics of the horizontal grid
C +   ------------------------------------------------

      write(6,*) 'Horizontal CPL grid created'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,200) mx,my,dx/1000.,dy/1000.,
     .             MinLon,MaxLon,MinLat,MaxLat
200   format(' Grid points           : ',i4,' * ',i4,/,
     .       ' Horizontal resol. (X) : ',f7.1,' km.',/,
     .       ' Horizontal resol. (Y) : ',f7.1,' km.',/,
     .       ' MAR longitude between : ',f7.2,' and ',f7.2,/,
     .       ' MAR latitude  between : ',f7.2,' and ',f7.2,/)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END

