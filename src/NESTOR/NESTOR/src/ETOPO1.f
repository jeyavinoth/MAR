C   +-------------------------------------------------------------------+
C   |  Subroutine ETOPO1                                Dec 12  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : Grid : NST__x and NST__y (longitude and latitude, degrees)|
C   | ^^^^^^^ ETOPO data set, resolution: 1 minutes                     |
C   |                                                                   |
C   | Output: NST_sh: surface elevation                                 |
C   | ^^^^^^^ NSTsol: land (4) / sea (1) mask                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE ETOPO1


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

      PARAMETER   (mlon  =  21601)
      PARAMETER   (mlat  =  10801)
      PARAMETER   (mlon1 =mlon+1)
      PARAMETER   (mlat1 =mlat+1)
C +...Size of full ETOPO file

      PARAMETER   (elat  =  4000)
      PARAMETER   (TOPOmx=  mlon)
      PARAMETER   (TOPOmy=  elat)
C +...Size of sub-domain (latitude only)

      INTEGER     DIMS(1),TOPO_ID,LAT_ID,LON_ID,sol,start(3),
     .            count(3),i1lon,i2lon,i1lat,i2lat,imlon,imlat,
     .            irien,ncid,Rcode,TOPOid

      INTEGER*4   etopo(mlon,elat)

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

      write(6,*) 'Topography : ETOPO1 data set (1 minutes)'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*)


C +---Open NetCDF file
C +   ----------------

      nbchar=1

      DO i=1,60
       IF (ETOPO_dir(i:i).ne.' ') nbchar=i
      ENDDO

      ETOPOfile = ETOPO_dir(1:nbchar-1) // '1/etopo1.nc'
      ncid = NCOPN(ETOPOfile,NCNOWRIT,Rcode)


C +---Find out the id of the variables
C +   --------------------------------

      LON_ID =NCVID(ncid,'x' ,Rcode)
      LAT_ID =NCVID(ncid,'y' ,Rcode)
      TOPO_ID=NCVID(ncid,'z',Rcode)


C +---Read latitudes and longitudes of ETOPO
C +   --------------------------------------

C +...! etopo_lon and _lat are real*8 !

      start(1)=1
      count(1)=mlon
C +        *****
      CALL NCVGT(ncid,LON_ID,start,count,etopo_lon,Rcode)
C +        *****
      DO i=1,mlon
        topo_lon(i) = etopo_lon(i)
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
     .                mx    ,my    ,NST__x,NST__y,NST_sh,tmp_in)
C +        ******


C +---Distinction between land and sea (further refined)
C +   --------------------------------

      DO j=1,my
      DO i=1,mx

       IF (NST_sh(i,j).lt.0.01) THEN
        NSTsol(i,j)=1
        NSTice(i,j)=0
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

!---------------------------------------------------------------------


      SUBROUTINE INTsimple (ni,nj,grd_Ix,grd_Iy,var_I,SPHgrd,
     .                      mx,my,grd_Ox,grd_Oy,var_O,tmp_in)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INTEGER    ns,i,j,ii,jj,p,q,is,ind0,ind1,nsamp,LocDim,
     .           ni,nj,mx,my,mmx,mmy,icent1,jcent1,icent2,jcent2,
     .           i1(4),j1(4),i2

      PARAMETER (ns = 5) ! Number of sampling points
      PARAMETER (LocDim=21601) ! Dim. of local 1D arrays

      REAL x,y,tmp,tmp2,x0,x1,y0,y1,epsi,AUXlon,MINlon,MAXlon,
     .     AUXlat,MINlat,MAXlat,dist_O,dist_I,AUXlo1,AUXlo2,
     .     AUXla1,AUXla2,dx,dy,degrad,ns2,tmp3

      REAL grd_Ix(ni),grd_Iy(nj),grd_Ox(mx,my),grd_Oy(mx,my),
     .     tmp_in(0:ni+1,0:nj+1),tmp_Ix(0:LocDim+1),samOx(ns),
     .     samOy(ns),tmp_Iy(0:LocDim+1),var_I(ni,nj),var_O(mx,my)

      LOGICAL SPHgrd,cyclic,npole,spole


C +---Data
C +   ----

      DATA epsi   / 1.d-4          /
      DATA degrad / 1.745329252d-2 /


C +---Check dimensions of temporary arrays
C +   ====================================

      IF (ni.gt.LocDim .or. nj.gt.LocDim) THEN
        WRITE(6,*) 'INTsimple - fatal error: dimension',LocDim
        WRITE(6,*) 'Please change LocDim   -   STOP'
        STOP
      ENDIF


C +---Check if the sampling technique is required
C +   ===========================================

      mmx = mx
      mmy = my
      
      dx    =(grd_Ix(ni/2)-grd_Ix(ni/2-1))*111111.
     .                   *COS(grd_Iy(nj/2)*degrad)
      dy    =(grd_Iy(nj/2)-grd_Iy(nj/2-1))*111111.
      dist_I=max(dx,dy)

      icent1=MAX(1,mx/2)
      icent2=MAX(1,mx/2-1)
      jcent1=MAX(1,my/2)
      jcent2=MAX(1,my/2-1)
      IF (mmx.eq.2) icent1=2
      IF (mmy.eq.2) jcent1=2
      
      AUXlo1=grd_Ox(icent1,jcent1)
CWARNINGXla1=grd_Oy(icent1,icent1)
      AUXla1=grd_Oy(icent1,jcent1)
      AUXlo2=grd_Ox(icent2,jcent2)
      AUXla2=grd_Oy(icent2,jcent2)

C +        ******
      CALL SPHERC (SPHgrd,AUXlo1,AUXla1)
      CALL SPHERC (SPHgrd,AUXlo2,AUXla2)
C +        ******

      dx    =(AUXlo1-AUXlo2)*111111.*COS(AUXla1*degrad)
      IF (mmx.le.1) dx = 1000.
      dy    =(AUXla1-AUXla2)*111111.
      IF (mmy.le.1) dy = 1000.
      dist_O=max(dx,dy)

      nsamp=1
      ns2  =  max(2.,(dist_O/dist_I))

      if(ns2==1) then
       print *,"WARNING: in INTsimple dist_O< dist_I!!"
      endif
       
C +---Coordinates indexes inversion (if necessary)
C +   ============================================


C +---Storage in temporary arrays
C +   ---------------------------

      DO jj=1,nj
      DO ii=1,ni 
       tmp_in(ii,jj)=var_I(ii,jj)
      ENDDO
      ENDDO

      DO ii=1,ni
       tmp_Ix(ii)=grd_Ix(ii)
      ENDDO

      DO jj=1,nj
       tmp_Iy(jj)=grd_Iy(jj)
      ENDDO


C +---Revert grd_Ix (1) <--> grd_Ix (n), ... ?
C +   ----------------------------------------
      
      IF (grd_Ix(ni).lt.grd_Ix(1)) THEN     
       DO ii=1,ni   
        DO jj=1,nj                       
         tmp_in(ii,jj)=var_I(ni-ii+1, jj)
        ENDDO
        tmp_Ix(ii)=grd_Ix(ni-ii+1) 
       ENDDO
      ENDIF


C +---Revert grd_Iy (1) <--> grd_Iy (n), ... ?
C +   ----------------------------------------
      
      IF (grd_Iy(nj).lt.grd_Iy(1)) THEN     
       DO jj=1,nj   
        DO ii=1,ni                       
         tmp_in(ii,jj)=var_I(ii,nj-jj+1)
        ENDDO
        tmp_Iy(jj)=grd_Iy(nj-jj+1)
       ENDDO
      ENDIF


C +---Extended coordinates in longitude and latitude
C +   ==============================================


C +---Check validity of longitude
C +   ---------------------------

      IF (SPHgrd) THEN
       IF ((tmp_Ix(1).lt.(-180.)).or.(tmp_Ix(ni).gt.180.)) THEN
        WRITE(6,*) 'Longitudes of data fields are not between'
        WRITE(6,*) '-180 and +180 deg. (as required by INTsimple)'
        WRITE(6,*) 'but rather between : '
        WRITE(6,*) tmp_Ix(1),tmp_Ix(ni)
        WRITE(6,*) '--- STOP in INTsimple ---'
        STOP
       ENDIF
      ENDIF


C +---Extended left/right boundaries (longitude if SPHgrd)
C +   ----------------------------------------------------

      tmp_Ix(0)   =2.*tmp_Ix( 1)-tmp_Ix(2)
      tmp_Ix(ni+1)=2.*tmp_Ix(ni)-tmp_Ix(ni-1) 


C +---Extended bottom/top boundaries (latitude if SPHgrd)
C +   ---------------------------------------------------

      tmp_Iy(0)   =2.*tmp_Iy( 1)-tmp_Iy(2)
      tmp_Iy(nj+1)=2.*tmp_Iy(nj)-tmp_Iy(nj-1)


C +---Define the cyclic field in longitude
C +   ------------------------------------

      IF (SPHgrd) THEN     ! Stereographic coordinates

       ind0=-1
       ind1=-1
 
       AUXlon=tmp_Ix(0)+360.
       DO i=1,ni
        IF (ABS(AUXlon-tmp_Ix(i)).lt.epsi) ind0=i
       ENDDO
 
       AUXlon=tmp_Ix(ni+1)-360.
       DO i=1,ni
        IF (ABS(AUXlon-tmp_Ix(i)).lt.epsi) ind1=i
       ENDDO
 
       IF (ind0.gt.(-1).and.ind1.gt.(-1)) THEN
        cyclic=.true.
       ELSE
        cyclic=.false.
        ind0=ni
        ind1= 1
       ENDIF
 
       IF (ABS(tmp_Ix(ni+1)-180.).lt.epsi) tmp_Ix(ni+1)=180.+epsi

      ELSE                 ! Non spherical coordinates

       ind0=ni
       ind1= 1

      ENDIF

      DO j=1,nj
       tmp_in(   0,j)=tmp_in(ind0,j)
       tmp_in(ni+1,j)=tmp_in(ind1,j)
      ENDDO
      
C +---Define extra lower and upper boundaries (latitude)
C +   --------------------------------------------------

      IF (SPHgrd) THEN     ! Stereographic coordinates

       IF (tmp_Iy(0).lt.(-90.))
     .  tmp_Iy(0)=MIN(-90.,tmp_Iy(1)-epsi)

       IF (tmp_Iy(nj+1).gt.90.)      
     .  tmp_Iy(nj+1)=MAX(90.,tmp_Iy(nj)+epsi)

       spole=.false.
       npole=.false.

       IF (tmp_Iy(0).le.(-90.)) spole=.true.
       IF (tmp_Iy(nj+1).ge.90.) npole=.true.

      ENDIF

      DO i=0,ni+1
       tmp_in(i,   0)=tmp_in(i, 1)
       tmp_in(i,nj+1)=tmp_in(i,nj)
      ENDDO
      

C +---Check the extension of the sub-domain to be read
C     ================================================

      AUXlon = grd_Ox(1,1)
      AUXlat = grd_Oy(1,1)
C +        ******
      CALL SPHERC (SPHgrd,AUXlon,AUXlat)
C +        ******
      MINlon = AUXlon
      MAXlon = AUXlon
      MINlat = AUXlat
      MAXlat = AUXlat

      DO j=1,my
      DO i=1,mx
       AUXlon = grd_Ox(i,j)
       AUXlat = grd_Oy(i,j)
C +         ******
       CALL SPHERC (SPHgrd,AUXlon,AUXlat)
C +         ******
   
       MINlon = min(AUXlon,MINlon)
       MAXlon = max(AUXlon,MAXlon)
       MINlat = min(AUXlat,MINlat)
       MAXlat = max(AUXlat,MAXlat)
      ENDDO
      ENDDO

      IF ((tmp_Ix(   0).gt.MINlon) .or.
     .    (tmp_Ix(ni+1).le.MAXlon) .or.
     .    (tmp_Iy(   0).gt.MINlat) .or.
     .    (tmp_Iy(nj+1).le.MAXlat)) THEN
       WRITE(6,*) 'Output domain is not (fully) included in'
       WRITE(6,*) 'the input domain.'
       WRITE(6,*) 'Input  domain :'
       WRITE(6,*) tmp_Ix(0),tmp_Ix(ni+1),tmp_Iy(0),tmp_Iy(nj+1)
       WRITE(6,*) 'Output domain :'
       WRITE(6,*) MINlon,MAXlon,MINlat,MAXlat
       WRITE(6,*) '--- STOP in INTsimple ---'
      ENDIF


C +---Bi-linear interpolation
C +   =======================


C +---Some initialisations
C +   --------------------

      p=0
      q=0
      
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO i=1,mx   ! LOOP on output grid-points : BEGIN
      DO j=1,my  
          
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Define sampling point positions
C +   -------------------------------


        DO is=1,nsamp   ! Boundaries : no sampling
         samOx(is)=grd_Ox(i,j)
         samOy(is)=grd_Oy(i,j)
        ENDDO


       tmp2=0.0       ! Initialisation of sum of sampled values

       DO is=1,nsamp  ! Loop on the sampling points: BEGIN

        x=samOx(is)
        y=samOy(is)


C +---Check the range of latitude and longitude
C +   -----------------------------------------

C +          ******
        CALL SPHERC (SPHgrd,x,y)
C +          ******


C +---Search for the bottom-left corner of the surrounding mesh
C +   ---------------------------------------------------------

C +...This simple method accounts for the fact that two successive
C +...requests usually refer to neighbour points


C +---Search for dimension 1 value
C +   ----------------------------

        IF (tmp_Ix(p).le.x) THEN   ! Search upwards
         DO WHILE (tmp_Ix(p+1).lt.x)
          p=p+1
         ENDDO
        ELSE                       ! Search downwards
         DO WHILE (tmp_Ix(p).gt.x)
          p=p-1
         ENDDO
        ENDIF


C +---Search for dimension 2 value
C +   ----------------------------

        IF (tmp_Iy(q).le.y) THEN  ! Search upwards
         DO WHILE (tmp_Iy(q+1).lt.y)
          q=q+1
         ENDDO
        ELSE                      ! Search downwards
         DO WHILE (tmp_Iy(q).gt.y)
          q=q-1
         ENDDO
        ENDIF


C +---Check the validity of p/q indexes
C +   ---------------------------------

        IF ((p.lt.     0).or.(q.lt.     0).or.
     .      (p.gt.(ni+1)).or.(q.gt.(nj+1))) THEN
         WRITE (6,*) 'Inconsistency between input and output'
         WRITE (6,*) 'domains.'
         WRITE (6,*) 'p and q = ',p,q
         WRITE (6,*) '--- STOP in INTsimple ---'
         STOP
        ENDIF


C +---Linear interpolation
C +   --------------------

        tmp2=0 ; tmp3=0

        do ii=nint(-1*ns2/2.),nint(ns2/2.) 
        do jj=nint(-1*ns2/2.),nint(ns2/2.)  

         x0=min(ni,max(1,p+ii))  
         y0=min(nj,max(1,q+jj))  

                                           tmp = 1.
        !if(max(abs(ii),abs(jj))>= ns2/2.) tmp = 2/3.

         tmp2 = tmp2 + tmp_in(x0,y0) * tmp
         tmp3 = tmp3 + tmp  
        enddo
        enddo

       ENDDO          ! LOOP on the sampling points: END


C +---Output value given by the average of the samplings
C +   --------------------------------------------------

       var_O(i,j)=tmp2/tmp3


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDDO
      ENDDO           ! Loop on output grid-points : END

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c #MS IF (cyclic) WRITE(6,*) 'INTbil-info: cyclic boundaries'
c #MS IF (npole ) WRITE(6,*) 'INTbil-info: North pole included'
c #MS IF (spole ) WRITE(6,*) 'INTbil-info: South pole included'


      RETURN
      END


