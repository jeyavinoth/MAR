C   +-------------------------------------------------------------------+
C   |  Subroutine ETOPOg                               June 03  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : Grid : NST__x and NST__y (longitude and latitude, degrees)|
C   | ^^^^^^^ ETOPO data set, resolution: 5 minutes                     |
C   |                                                                   |
C   | Output: NST_sh: surface elevation                                 |
C   | ^^^^^^^ NSTsol: land (4) / sea (1) mask                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE ETOPOg


      IMPLICIT NONE


C +---Netcdf specifications
C +   ---------------------

      INCLUDE 'NetCDF.inc'


C +---NST variables
C +   -------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'LOCfil.inc'


C +---Local variables
C +   ---------------

      INTEGER      i,j,mlon,mlat,elat,midlon,mlon1,mlat1,
     .             nbchar,TOPOmx,TOPOmy

      PARAMETER   (mlon  =  4320)
      PARAMETER   (mlat  =  2160)
      PARAMETER   (mlon1 =mlon+1)
      PARAMETER   (mlat1 =mlat+1)
C +...Size of full ETOPO file

      PARAMETER   (elat  =  1300)
      PARAMETER   (TOPOmx=  mlon)
      PARAMETER   (TOPOmy=  elat)
C +...Size of sub-domain (latitude only)

      INTEGER     DIMS(1),TOPO_ID,LAT_ID,LON_ID,sol,start(3),
     .            count(3),i1lon,i2lon,i1lat,i2lat,imlon,imlat,
     .            irien,ncid,Rcode,TOPOid

      INTEGER*2   etopo(mlon,elat)

      REAL*8      etopo_lon(mlon), etopo_lat(mlat)

      REAL        topo_lon(mlon),topo_lat(mlat),size_lon,
     .            TOPlon(TOPOmx),TOPlat(TOPOmy),size_lat,
     .            TOPsh(TOPOmx,TOPOmy),tmpTOP(TOPOmx,TOPOmy),
     .            tmp_in(0:TOPOmx+1,0:TOPOmy+1),MINlon,MINlat,
     .            MAXlon,MAXlat,center_lon,center_lat,AUXlon,AUXlat

      LOGICAL     Vtrue

      CHARACTER*80 ETOPOfile


C +---Data
C +   ----

      DATA start  / 1,1,1/
      DATA Vtrue  /.true./


C +---Opening and reading of ETOPO data file
C +   ======================================

      write(6,*) 'Topography : ETOPO data set (5 minutes)'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*)


C +---Open NetCDF file
C +   ----------------

      nbchar=1

      DO i=1,60
       IF (ETOPO_dir(i:i).ne.' ') nbchar=i
      ENDDO

      ETOPOfile = ETOPO_dir(1:nbchar) // 'etopo.nc'
      ncid = NCOPN(ETOPOfile,NCNOWRIT,Rcode)


C +---Find out the id of the variables
C +   --------------------------------

      LON_ID =NCVID(ncid,'lon' ,Rcode)
      LAT_ID =NCVID(ncid,'lat' ,Rcode)
      TOPO_ID=NCVID(ncid,'topo',Rcode)


C +---Read latitudes and longitudes of ETOPO
C +   --------------------------------------

C +...! etopo_lon and _lat are real*8 !

      start(1)=1
      count(1)=mlon
C +        *****
      CALL NCVGT(ncid,LON_ID,start,count,etopo_lon,Rcode)
C +        *****
      midlon = mlon/2
      DO i=1,midlon
        topo_lon(i+midlon) = etopo_lon(i)
      END DO
      DO i=midlon+1,mlon
        topo_lon(i-midlon) = etopo_lon(i) - 360.
      END DO
C +...etopo_lon : from    0 to 360 deg.
C +...topo_lon  : from -180 to 180 deg.

      start(1)=1
      count(1)=mlat
C +        *****
      CALL NCVGT(ncid,LAT_ID,start,count,etopo_lat,Rcode)
C +        *****
      DO j=1,mlat
        topo_lat(j) = etopo_lat(mlat1-j)
      END DO
C +...etopo_lat : from  90 to -90 deg.
C +...topo_lat  : from -90 to  90 deg.


C +---Compute the extension of the sub-domain to be read
C     --------------------------------------------------

      AUXlon = NST__x(1,1)
      AUXlat = NST__y(1,1)
C +        ******
      CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +        ******
      MINlon = AUXlon
      MAXlon = AUXlon
      MINlat = AUXlat
      MAXlat = AUXlat
      DO j=1,my
      DO i=1,mx
       AUXlon = NST__x(i,j)
       AUXlat = NST__y(i,j)
C +         ****** 
       CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +         ******
       MINlon = min(AUXlon,MINlon)
       MAXlon = max(AUXlon,MAXlon)
       MINlat = min(AUXlat,MINlat)
       MAXlat = max(AUXlat,MAXlat)
      ENDDO
      ENDDO


C +---Define extraction zone
C +   ----------------------

C +        ******
      CALL SEARCH (topo_lon,mlon,MINlon,i1lon,irien)
      CALL SEARCH (topo_lon,mlon,MAXlon,irien,i2lon)
C +        ******
      imlon = i2lon - i1lon + 1
C +        ******
      CALL SEARCH (topo_lat,mlat,MINlat,i1lat,irien)
      CALL SEARCH (topo_lat,mlat,MAXlat,irien,i2lat)
C +        ******
      imlat = i2lat - i1lat + 1

      IF (imlat.ge.elat) THEN
       write(*,*) 'Extent of the simulation domain in latitude'
       write(*,*) 'is too large. Please choose a larger value '
       write(*,*) 'for the elat parameter.             - STOP '
       STOP
      ENDIF

      i1lat = i1lat + (i2lat-i1lat)/2 - elat/2
      i1lat = MAX(1,i1lat)
      i2lat = i1lat + elat - 1
      IF (i2lat.gt.mlat) THEN
       i2lat= mlat
       i1lat= i2lat - elat + 1
      ENDIF


C +---Read values of the variables for the sub-domain
C +   -----------------------------------------------

      start(1)=1
      start(2)=mlat-i2lat+1
      count(1)=mlon
      count(2)=elat

C +        *****
      CALL NCVGT(ncid,TOPO_ID,start,count,etopo,Rcode)
C +        *****

      DO i=1,midlon
        DO j=1,elat
        TOPsh(i+midlon,j) = etopo(i,elat+1-j)
        END DO
      END DO
      DO i=midlon+1,mlon
        DO j=1,elat
        TOPsh(i-midlon,j) = etopo(i,elat+1-j)
        END DO
      END DO

C +        ******
      CALL NCCLOS (ncid,Rcode)
C +        ******

      DO i=1,TOPOmx
       TOPlon(i)=topo_lon(i)
      ENDDO

      DO j=1,TOPOmy
       TOPlat(j)=topo_lat(i1lat-1+j)
      ENDDO


C +---Interpolation of topography to the NST grid
C +   -------------------------------------------

C +        ******
      CALL INTbil (TOPOmx,TOPOmy,TOPlon,TOPlat,TOPsh ,Vtrue ,
     .             mx    ,my    ,NST__x,NST__y,NST_sh,tmp_in)
C +        ******


C +---Distinction between land and sea (further refined)
C +   --------------------------------

      DO j=1,my
      DO i=1,mx

       IF (NST_sh(i,j).lt.0.01) THEN
        NSTsol(i,j)=1
       ELSE
        NSTsol(i,j)=4
       ENDIF

      ENDDO
      ENDDO

C +---Special topo for Greenland Simulation
C +   -------------------------------------

      IF (region.eq."GRD") call USRgrd ('ETOPOg')

C +---No atmosphere below sea level...
C +   --------------------------------

      DO j=1,my
      DO i=1,mx

       IF (NST_sh(i,j).lt.0.0) THEN
        NST_sh(i,j)= 0.0
       ENDIF

      ENDDO
      ENDDO


      RETURN
      END


C +   ********************************************************************


      subroutine SEARCH (xx,ni,xreq,KLO,KHI)

      REAL xx(ni)

      KLO=1
      KHI=ni
 1    IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(xx(K).GT.xreq)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

      RETURN
      END


