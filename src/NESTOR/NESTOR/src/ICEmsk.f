C   +-------------------------------------------------------------------+
C   |  Subroutine ICEmsk                                Dec 12  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : Grid : NST__x and NST__y (longitude and latitude, degrees)|
C   | ^^^^^^^ ETOPO data set, resolution: 1 minutes                     |
C   |                                                                   |
C   | Output: NST_sh: surface elevation                                 |
C   | ^^^^^^^ NSTsol: land (4) / sea (1) mask                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE ICEmsk


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
      INCLUDE 'MARvar.inc'
 

C +---Local variables
C +   ---------------

      INTEGER      i,j,mlon,mlat,elat,midlon,mlon1,mlat1,
     .             nbchar,TOPOmx,TOPOmy

      PARAMETER   (mlon  =  1440)
      PARAMETER   (mlat  =  720)
      PARAMETER   (mlon1 =mlon+1)
      PARAMETER   (mlat1 =mlat+1)
C +...Size of full ETOPO file

      PARAMETER   (elat  =  360)
      PARAMETER   (TOPOmx=  mlon)
      PARAMETER   (TOPOmy=  elat)
C +...Size of sub-domain (latitude only)

      INTEGER     DIMS(1),TOPO_ID,LAT_ID,LON_ID,sol,start(3),
     .            count(3),i1lon,i2lon,i1lat,i2lat,imlon,imlat,
     .            irien,ncid,Rcode,TOPOid

      REAL*4      etopo(mlon,elat)

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


C +---Opening and reading of ICEmask data file
C +   ======================================

      write(6,*) 'Topography : ICEmask data set (4 minutes)'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*)


C +---Open NetCDF file
C +   ----------------

      nbchar=1

      DO i=1,60
       IF (ETOPO_dir(i:i).ne.' ') nbchar=i
      ENDDO

      ETOPOfile = 'input/ICEmask/ICEmask.nc'
      ncid = NCOPN(ETOPOfile,NCNOWRIT,Rcode)


C +---Find out the id of the variables
C +   --------------------------------

      LON_ID =NCVID(ncid,'LON' ,Rcode)
      LAT_ID =NCVID(ncid,'LAT' ,Rcode)
      TOPO_ID=NCVID(ncid,'ICE',Rcode)


C +---Read latitudes and longitudes of ETOPO
C +   --------------------------------------

C +...! etopo_lon and _lat are real*8 !

      start(1)=1
      count(1)=mlon
C +        *****
      CALL NCVGT(ncid,LON_ID,start,count,etopo_lon,Rcode)
C +        *****
      DO i=1,mlon
        topo_lon(i) = etopo_lon(i)-360.
      END DO
C +...topo_lon  : from -180 to 180 deg.

      start(1)=1
      count(1)=mlat
C +        *****
      CALL NCVGT(ncid,LAT_ID,start,count,etopo_lat,Rcode)
C +        *****
      DO j=1,mlat
        topo_lat(j) = etopo_lat(j)
      END DO
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
      start(2)=max(1,i1lat-1)
      count(1)=mlon
      count(2)=elat

C +        *****
      CALL NCVGT(ncid,TOPO_ID,start,count,etopo,Rcode)
C +        *****

      DO i=1,mlon
        DO j=1,elat
        TOPsh(i,j) = etopo(i,j)
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
      CALL INTsimple (TOPOmx,TOPOmy,TOPlon,TOPlat,TOPsh ,Vtrue ,
     .                mx    ,my    ,NST__x,NST__y,NSTice,tmp_in)
C +        ******


C +---Distinction between land and sea (further refined)
C +   --------------------------------

      DO j=1,my
      DO i=1,mx

       IF (NSTsol(i,j)<=2) THEN
        NSTice(i,j)=0
       ENDIF

      ENDDO
      ENDDO

C +---Special topo for Greenland Simulation
C +   -------------------------------------

      IF (region.eq."GRD") call USRgrd ('ETOPOg')


      RETURN
      END

!---------------------------------------------------------------------




