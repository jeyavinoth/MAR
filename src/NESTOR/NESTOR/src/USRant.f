C   +-------------------------------------------------------------------+
C   |  Subroutine USRant          June 2003 (HG) Oct 2014 (CA)  NESTING |
C   +-------------------------------------------------------------------+
C   | USRant adapt NESTOR to Antarctic region                           |
C   | See bellow a new version (AntTOPO) with Bamber (2009)             |
C   |                                                                   |
C   | Input : - subnam : Name of the subroutine                         |
C   | ^^^^^^^            where USRant is called                         |
C   |                                                                   |
C   | Maintainer : Hubert Gallée, Cécile Agosta                         |
C   | ^^^^^^^^^^^^                                                      |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE USRant (subnam)


      IMPLICIT NONE

C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'LSCvar.inc'

C +---local variables
C +   ---------------

      CHARACTER*6 subnam
  
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Topography for Antarctic    
C +   ========================

      IF (subnam.eq.'NESTOR') call ANTOPO  

      IF (subnam.eq.'Bamber') call AntTOPO  !+ CA +!

      If (subnam.eq.'Racmo2') call AntRACMO !+ CA +!

      END SUBROUTINE USRant

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
C   +-------------------------------------------------------------------+
      Subroutine AntRACMO
C   +-------------------------------------------------------------------+

C +---Netcdf specifications
C +   ---------------------
      INCLUDE 'NetCDF.inc'
      
C +---NST variables
C +   -------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'LOCfil.inc'
      
      Real NST_dh
      
C +...Reading dem file        
      Integer fID
      Character*80 dem_file
      Integer dem_mx,dem_my
      Parameter(dem_mx=110,dem_my=91)
      Integer dem_xmin,dem_ymin,dem_xmax,dem_ymax
      Parameter(dem_xmin=-2600,dem_ymin=-2200)
      Parameter(dem_xmax= 2850,dem_ymax= 2300)
      Character*80   var_units
      Real           dem_sh(dem_mx,dem_mx),dem_msk(dem_mx,dem_mx)
      
      Integer i,j,ii,jj

      NST_dh  = NST_dx/1000.
      
      write(6,*) 'Topography : RACMO-27 interpolated on 50km MAR grid'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    
      dem_file = trim(PFXdir)//'AntTOPO/RACMO-ANT27onMAR-50km.sh.nc'
      write(6,*) "> Read ", dem_file,"  ..."
      write(6,*)
      
      CALL UNropen (dem_file,fID,title)
      CALL UNsread(fID,'sh',1,1,1,1,
     .                    dem_mx,dem_my,1,
     .                    var_units,dem_sh)
     
      CALL UNsread(fID,'msk',1,1,1,1,
     .                    dem_mx,dem_my,1,
     .                    var_units,dem_msk)

      Do j=1,my
      Do i=1,mx
        If ( NSTgdx(i).ge.dem_xmin .and. NSTgdy(j).ge.dem_ymin
     .  .and.NSTgdx(i).le.dem_xmax .and. NSTgdy(j).le.dem_ymax) then
          ii = int((NSTgdx(i)-dem_xmin)/NST_dh) +1
          jj = int((NSTgdy(j)-dem_ymin)/NST_dh) +1
C           print*,i,j,ii,jj
          NST_sh(i,j)=dem_sh (ii,jj)
          NSTice(i,j)=dem_msk(ii,jj)*100.
        else
          NST_sh(i,j)=0
          NSTice(i,j)=0.
        EndIf
        ! +---No atmosphere below sea level
        ! +   -----------------------------
        NST_sh(i,j) = max(NST_sh(i,j),0.)
        If (NSTice(i,j).gt.50.) then
          NSTsol(i,j) = 3
        else
          NSTsol(i,j) = 1
        EndIf
      EndDo
      EndDo
      
      End Subroutine AntRACMO
C   +-------------------------------------------------------------------+


C   +-------------------------------------------------------------------+
C   |  Subroutine ANTOPO                          October 2002  NESTING |
C   |             ANTARTIC TOPOGRAPHY specific Assimilation             |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : NST__x : longitude (degree) of the NST grid               |
C   | ^^^^^^^ NST__y : latitude  (degree) of the NST grid               |
C   |                                                                   |
C   | Output: NST_sh: surface elevation                                 |
C   | ^^^^^^^ NSTsol: land (4) / sea (1) mask                           |
C   |                                                                   |
C   | Method: Divide each nested Grid Cell in elementary Meshes         |
C   | ^^^^^^^ much small than either the nested or the DTM Mesh         |
C   |         Compute geographic Coordinates of each elementary Mesh    |
C   |         Compute Distance of this Mesh to 4 closest DTM Meshes     |
C   |         Compute Height of this Mesh by kriging                    |
C   |         Nested Grid Cell Height is the Average                    |
C   |                                 of the elementary Meshes Heights  |
C   |                                                                   |
C   | DATA Source: Radarsat Antarctic Mapping Project Digital           |
C   | ^^^^^^^^^^^^ Digital  Elevation Model   Version 2                 |
C   |              ftp site : sidads.colorado.edu                       |
C   |              directory: /pub/DATASETS/RAMP/DEM_V2/1KM/ASCII       |
C   |              file     : ramp1kmdem_wgsosu_v2.txt.gz               |
C   | Reference:   Liu, H., K. Jezek, B. Li, and Z. Zhao. 2001.         |
C   | ^^^^^^^^^^   Radarsat Antarctic Mapping Project                   |
C   |              Digital Elevation Model Version 2.                   |
C   |              Boulder, CO: National Snow and Ice Data Center.      |
C   |              Digital media.                                       |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE ANTOPO


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'MARvar.inc'
      INCLUDE 'LOCfil.inc'
      INCLUDE 'NetCDF.inc'

      INTEGER   nxdata,nydata ,nwdata,nndata
      PARAMETER(nxdata=13670204) 
      INTEGER   ii,i,j

      REAL      ANT_sh(nxdata)
      REAL      ANTlat(nxdata),ANTlon(nxdata)

      REAL      LOC_sh(100000),ddxx          ,ddll       ,dd_min
      REAL      LOC__x(100000),LOC__y(100000)

      REAL      pi           ,degrad       ,t_grad       ,t__rad
      REAL      phi          ,angl
      REAL      x_st_i       ,y_st_j       ,x_st_n       ,y_st_n
      REAL                                  x_st_d       ,y_st_d
      REAL      di_lon       ,di_lat       ,dj_lon       ,dj_lat
      REAL      inilon       ,inilat
      REAL      curlon       ,curlat       ,earthr
      REAL      minLON       ,minLAT       ,difLON       ,difLAT
      INTEGER   inx          ,jny
      INTEGER   in           ,jn

      LOGICAL   dd00
      INTEGER   nn           ,iinn(4)      ,jjnn(4)      ,n0
      REAL      x0           ,y0           ,ddxxi        ,ddxxj
      REAL      xx           ,yy           ,xxii         ,yyii
      REAL      dd           ,ddnn(4)                    ,hh


      DATA      earthr/6371.e3/           ! Earth Radius
      DATA      t_grad/ 360.  /           ! 


      pi     =  acos( -1.)
      degrad =  pi / 180.
      t__rad =  pi *   2.
      nndata =         0

 
C +---INPUT
C +   -----

      IF (TOP30.eq.'a   ')                                        THEN
      write(6,*)  'ANTOPO: INPUT File rampxkmdem_wgsosu_v2.dat OPENING'
      open (unit=1,status='old',file='rampxkmdem_wgsosu_v2.dat')
      rewind     1
        DO  j=1,my
        DO  i=1,mx
        read(1,100) NST_sh(i,j),NSTsol(i,j)
 100    format(e15.6,i15)
        END DO
        END DO
      write(6,*)  'ANTOPO: INPUT File rampxkmdem_wgsosu_v2.dat IN'
      close(unit=1)
      ELSE
      write(6,*)  'ANTOPO: INPUT File ramp1kmdem_wgsosu_v2.bin OPENING'
      open (unit=1,status='old',file='ramp1kmdem_wgsosu_v2.bin',
     .      form='unformatted')
      rewind     1
           ii= 0
 1001   CONTINUE
           ii=1+ii
c #DO   DO ii=1,nxdata
        read(1,end=1000) ANTlat(ii),ANTlon(ii),ANT_sh(ii)
        GO TO 1001
 1000   CONTINUE
c #DO   END DO
                nwdata = ii
      close(unit=1)
      write(6,*)  'ANTOPO: INPUT File ramp1kmdem_wgsosu_v2.bin IN'
      write(6,*)  '        Nb DATA: ',nwdata
      write(6,*)  '                 '


C +---Finest Resolution
C +   -----------------

            ddll= abs(NST__x(2,2)-NST__x(2,1))
        IF (ddll .GT. t_grad) 
     .      ddll=ddll-t_grad
            ddxx = 
     .         ((ddll*cos(NST__y(2,2)*degrad     )*degrad)**2
     .         +(     abs(NST__y(2,2)-NST__y(2,1))*degrad)**2)
            ddxx =   sqrt(ddxx)      *earthr
            inx  =        ddxx       *1.e-3
            jny  =        ddxx       *1.e-3
            write(6,600)  ddxx       *1.e-3,inx
 600        format(8x,'dx =',f9.3,'km   inx =',i6)


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO j=1,my                               ! Loop on NST grid points
      DO i=1,mx

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Interpolation/Average
C +   ---------------------
          ddll   = earthr * cos(NST__y(i,j)*degrad)
          angl   = 0.5*pi -     NST__x(i,j)*degrad
          x_st_i = ddll   * cos(angl)
          y_st_j = ddll   * sin(angl)


C +---Interpolation/Average
C +   ---------------------

C +---Relevant Points
C +   ~~~~~~~~~~~~~~~
          nndata = nndata + 1
          nydata = 0
          dd_min = 1.e15 
        DO ii=1,nwdata
          phi    =              ANTlat(ii) *degrad
          ddll   = earthr * cos(phi)
          angl   = 0.5*pi -     ANTlon(ii) *degrad
          x_st_n = ddll   * cos(angl)
          y_st_n = ddll   * sin(angl)
          x_st_d =(x_st_n - x_st_i)
          y_st_d =(y_st_n - y_st_j)
          ddll   = x_st_d * x_st_d 
     .           + y_st_d * y_st_d
          ddll   = sqrt(ddll)
          dd_min =  min(ddll,dd_min)
          IF (ddll .LT. max(0.71*ddxx,4.e3))                      THEN
              nydata         = nydata + 1
              LOC__x(nydata) = x_st_n
              LOC__y(nydata) = y_st_n
              LOC_sh(nydata) = ANT_sh(ii)
          END IF
        ENDDO
        IF (mod(nndata,20).eq.1)
     .    write(6,601) NST__x(i,j),x_st_i*1.e-3,
     .                 NST__y(i,j),y_st_j*1.e-3,
     .                 dd_min*1.e-3,nydata
 601      format(9x,'LON,LAT =',2(f12.3,'  (',f9.3,')'),
     .              '  (',f9.3,')',i8,' Points')


C +---Elementary Meshes Characteristics
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        NST_sh(i,j) = 0.

          x0 = x_st_i - 0.5 *ddxx
          y0 = y_st_j - 0.5 *ddxx
          ddxxi = ddxx/inx
          ddxxj = ddxx/jny 
        DO jn=1,jny
        DO in=1,inx
          xx    = x0 + in*ddxxi
          yy    = y0 + jn*ddxxj

C +---Distances to the current elementary Mesh
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            DO nn=1,4
                  ddnn(nn)=nn*earthr     ! arbitrary large distance
            END DO
          DO ii=1,nydata
            xxii = xx - LOC__x(ii)
            yyii = yy - LOC__y(ii)
            dd   = xxii*xxii + yyii*yyii
            dd   = sqrt(dd)

C +---Closest DTM points
C +   ~~~~~~~~~~~~~~~~~~
                  dd00     =.true.
            DO nn=1,4
              IF (dd.lt.ddnn(nn).AND.dd00)                        THEN
                IF (nn.lt.4)                                      THEN
                  DO n0=4,nn+1,-1
                  ddnn(n0) = ddnn(n0-1)
                  iinn(n0) = iinn(n0-1)
                  END DO
                END IF
                  ddnn(nn) = dd
                  iinn(nn) = ii
                  dd00     =.false.
              END IF
            END DO
          ENDDO


C +---Kriging
C +   ~~~~~~~
                  dd = 0.
                  hh = 0.
            DO nn=1,4
            IF   (ddnn(nn).LT.1500.)                           THEN
                  hh = hh 
     .               + LOC_sh(iinn(nn))/ddnn(nn)
                  dd = dd    +      1. /ddnn(nn)
            END IF
            END DO
            IF   (dd      .GT.   0.)
     .            hh = hh                       /dd
                  NST_sh(i,j) = NST_sh(i,j)     +hh
        ENDDO
        ENDDO

C +---Nested Grid Cell Average
C +   ~~~~~~~~~~~~~~~~~~~~~~~~
                  NST_sh(i,j) = NST_sh(i,j) / (inx*jny)
        IF (mod(nndata,20).eq.1)
     .    write(6,602)   i,j,   NST_sh(i,j)
 602      format(9x,i3,i4,f14.3)


C +---Distinction between land and sea (further refined)
C +   --------------------------------

       IF (NST_sh(i,j).lt.0.01) THEN
        NSTsol(i,j)=1
       ELSE
        NSTsol(i,j)=3
       ENDIF


C +---No atmosphere below sea level...
C +   --------------------------------

       IF (NST_sh(i,j).lt.0.0) THEN
        NST_sh(i,j)= 0.0
       ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDDO                                   ! Loop on NST grid points
      ENDDO

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---OUTPUT
C +   ------

      write(6,*) 'ANTOPO: OUTPUT File rampxkmdem_wgsosu_v2.dat OPENING'
      open (unit=2,status='new',file='rampxkmdem_wgsosu_v2.dat')
      rewind     2
        DO  j=1,my
        DO  i=1,mx
        write(2,100) NST_sh(i,j),NSTsol(i,j)
        END DO
        END DO
      close(unit=1)
      write(6,*) 'ANTOPO: OUTPUT File rampxkmdem_wgsosu_v2.dat OUT'
      write(6,*) ' '
      write(6,*) ' '


      END IF


      RETURN
      END SUBROUTINE ANTOPO
      
C    +-------------------------------------------------------------------+
      SUBROUTINE AntTOPO
C    | ANTARTIC TOPOGRAPHY specific Assimilation - Bamber 2009           |
C    +-------------------------------------------------------------------+
C    | Cecile Agosta, October 2012                                       |
C    | Input : NST__x : longitude (degree) of the NST grid               |
C    | ^^^^^^^ NST__y : latitude  (degree) of the NST grid               |
C    |                                                                   |
C    | Output: NST_sh: surface elevation                                 |
C    | ^^^^^^^ NSTsol: land (4) / ice(3) / sea (1) mask                  |
C    |         NSTice: percentage of grid cell covered by ice [0-100]    |
C    |                                                                   |
C    | Method:                                                           |
C    | ^^^^^^^                                                           |
C    |                                                                   |
C    | Data  :                                                           |
C    | ^^^^^^^                                                           |
C    | File name        : krigged_dem_nsidc.bin
C    | Dataset Title    : Antarctic 1 km Digital Elevation Model (DEM)   |
C    |                    from Combined ERS-1 Radar                      |
C    |                    and ICESat Laser Satellite Altimetry           |
C    | Dataset Creator  : Bamber, Jonathan L., Jose Luis Gomez-Dans,     |
C    |                    and Jennifer A. Griggs                         |
C    | Dataset Release Place: Boulder, Colorado USA                      |
C    | Dataset Publisher: NSIDC > National Snow and Ice Data Center      |
C    | Online Resource  : http://nsidc.org/data/nsidc-0422.html          |
C    |                                                                   |
C    | Reference:  Bamber, J.L., J.L. Gomez-Dans, and J.A. Griggs. 2009. |
C    | ^^^^^^^^^^  A New 1 km Digital Elevation Model of the Antarctic   |
C    |             Derived from Combined Satellite Radar and Laser Data  |
C    |             Ð Part 1: Data and Methods.The Cryosphere,3,101-111.  |
C    |                                                                   |
C    |             Griggs, J. A. and J. L. Bamber. 2009.                 |
C    |             A New 1 km Digital Elevation Model of Antarctica      | 
C    |             Derived from Combined Radar and Laser Data Ð Part 2:  |
C    |             Validation and Error Estimates.The Cryosphere,3,113Ð123.
C    |                                                                   |
C    | Metadata :  krigged_dem_nsidc.bin.hdr                             |
C    | ^^^^^^^^^^                                                        |
C    |    description = { File Imported into ENVI.}                      |
C    |    samples = 5601                                                 |
C    |    lines   = 5601                                                 |
C    |    bands   = 1                                                    |
C    |    header offset = 0                                              |
C    |    file type = ENVI Standard                                      |
C    |    data type = 4                                                  |
C    |    interleave = bsq                                               |
C    |    sensor type = Unknown                                          |
C    |    byte order = 0                                                 |
C    |    map info = {Polar Stereographic South, 1, 1,                   |
C    |           -2800500., 2800500., 1000., 1000., WGS-84, units=Meters}|
C    |    projection info = {31, 6378137.0, 6356752.3, -71., 0.,         |
C    |         0.0, 0.0, WGS-84, Polar Stereographic South, units=Meters}|
C    |                                                                   |
C    +-------------------------------------------------------------------+

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
      Integer      i,j,in,jn,ii,jj,maptyp
      Integer      dem_ii,dem_jj,di0,dj0,ndh
      Real         dem_max,dem_min,dem_dh
      Real         dem_x0,dem_y0
      Real         lon,lat,GEddxx,xx,yy,sh,ice
      Real         nul
C +...Size of full dem file
      Integer      dem_mx,dem_my
      Parameter   (dem_mx=  5601)
      Parameter   (dem_my=  5601)
C +...dem variables      
      Real*4       dem_sh_1km(dem_mx,dem_my)
C +...Reading dem file
      Character*80 dem_file
      Integer      recl      ! Record length
C +...Local NST variables
      Real         NST_dh
      Real         NST_xx(mx),NST_yy(my)

C +---Opening and reading Bamber 2009 data file
C +   =========================================

      write(6,*) 'Topography : Bamber (2009) Antarctic 1km DEM'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

C +---Open BSQ binary file
C +   --------------------
      dem_file = trim(PFXdir)//'AntTOPO/krigged_dem_nsidc.bin'
      write(6,*) "> Read ", dem_file,"  ..."
      write(6,*)
      Inquire(iolength=recl) dem_sh_1km
      Open (unit=10,status='old',file=dem_file,form='unformatted',
     .      access='direct', recl=recl)
      Read (10,rec=1) dem_sh_1km
      Close(10)
C +---Characteristics of the DEM grid
C +   -------------------------------
      dem_x0  = -2800.5
      dem_y0  =  2800.5
      dem_dh  =  1.
      
C +---Verify than the map projection used for MAR 
C +   is the same than in the DEM
C +   ========================================

      OPEN (unit=51,status='old',file='MARgrd.ctr')

       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) maptyp
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) nul
       read (51,*) nul  
       read (51,*) nul
       read (51,*) nul 
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) nul
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) GEddxx
       read (51,*) !- - - - - - - - - - - - - - - - - -

      CLOSE(unit=51)
      
      If (maptyp.ne.0) then
        Print*, "++++++++++++++++++++++++++++++++++"
        Print*, "Routine USRant.f ----> Warning !!!"
        Print*, "For Antarctica,                   "
        Print*, "choose map type = 0 in MARgrd.ctr "
        Print*, "(Standard stereo south projection)"
        Print*, "NESTOR stopped NOW !!!            "
        Print*, "++++++++++++++++++++++++++++++++++"
        STOP
      EndIf
      
      If (GEddxx.eq.90) then
        Do j=1,my
          NST_yy(j) = NSTgdy(j)
        EndDo
        Do i=1,mx
          NST_xx(i) = NSTgdx(i)
        EndDo
      else
        Do j=1,my
        Do i=1,mx
          lon = NST__x(i,j)
          lat = NST__y(i,j)
          call StereoSouth_inverse(lon,lat,90.,xx,yy)
          NST_xx(i) = xx
          NST_yy(j) = yy
        EndDo
        EndDo
      EndIf
      
C +---Average the DEM on the MAR grid
C +   ===============================
      NST_dh  = NST_dx/1000.
      ndh     = nint(NST_dh/dem_dh)
      di0     = int((ndh+1)/2.)+1
      dj0     = int((ndh+1)/2.)+1
      dem_min = -2800.5 + NST_dh/2.
      dem_max =  2800.5 - NST_dh/2.
      Do j=1,my
      Do i=1,mx
        sh  = 0.
        ice = 0.
        If (     NST_xx(i).lt.dem_min 
     .      .or. NST_xx(i).gt.dem_max
     .      .or. NST_yy(j).lt.dem_min 
     .      .or. NST_yy(j).gt.dem_max) then
          NST_sh(i,j) = 0.
          NSTice(i,j) = 0.
          NSTsol(i,j) = 1
        else
          dem_ii = nint(NST_xx(i)  - dem_x0 + 1)
          dem_jj = nint(dem_y0 - NST_yy(j)  + 1)
          Do jn=1,ndh
          Do in=1,ndh
            ii = dem_ii - di0 + in
            jj = dem_jj - dj0 + jn
            If (dem_sh_1km(ii,jj).ge.-990.) then
              ice = ice + 1
              sh  = sh + dem_sh_1km(ii,jj)
            EndIf
          EndDo
          EndDo
          NST_sh(i,j) = sh  / (ndh*ndh)
          ! +---No atmosphere below sea level
          ! +   -----------------------------
          NST_sh(i,j) = max(NST_sh(i,j),0.)
          NSTice(i,j) = ice / (ndh*ndh) *100.
          If (NSTice(i,j).gt.50.) then
            NSTsol(i,j) = 3
          else
            NSTsol(i,j) = 1
          EndIf
        EndIf
      EndDo
      EndDo

      Return
      END SUBROUTINE AntTOPO

! +----------------------------------------------------------------------+
      Subroutine StereoSouth_inverse (lon,lat,lonE,E,N)
! |  Compute Oblique Stereographic Projection from lon,lat coordinates   |
! |  Written by Cecile Agosta                                 17-05-10   |
! |  EPSG Polar Stereographic transformation Variant B                   |
! |  (http://www.epsg.org/guides/docs/G7-2.pdf)                          |
! |  Equivalent to EPSG 3031 (WGS-84 ellipsoid)                          |
! +----------------------------------------------------------------------+
! |                                                                      |
! | INPUT : lon     : Longitude (deg)                                    |
! | ^^^^^^^ lat     : Latitude  (deg)                                    |
! |         lon0    : Longitude of X axis (90 = East, clockwise)         |
! |        [lat true = 71 S]                                             |
! |                                                                      |
! | OUTPUT : E      : Stereo coordinate on the East  (X, km)             |
! | ^^^^^^^^ N      : Stereo coordinate on the North (Y, km)             |
! |                                                                      |
! +----------------------------------------------------------------------+
      IMPLICIT NONE

      INCLUDE 'NSTdim.inc'

! +-- General Variables
! +   -----------------
      Real,INTENT(in ) :: lon,lat,lonE
      Real,INTENT(out) :: E,N

! +-- Local Variables
! +   ---------------
      Real  costru,ddista

! +-- Constants
! +   ---------
      Real  aa,ex,pi,degrad,latF,FE,FN,tF,mF,k0,t,rho,lonrad,latrad
      Real  lon0,trulat
      
      aa     = 6378.1370         ! aa (km) = demi grand axe
      ex     = 0.081819190842621 ! excentricity WGS-84 : 0.081 819 190 842 622 0.081 819 190 842 621
      trulat = -71.              ! Latitude of standard parallel, 71 S for ESPG 3031
      pi     = 4. * atan(1.)
      degrad = pi / 180.

      latF   = trulat*degrad ! Latitude of standard parallel, 71 for ESPG 3031  
      lon0   = (lonE-90)*degrad
      lonrad = lon *degrad
      latrad = lat *degrad      
      
      FE = 0. !False Easting
      FN = 0. !False Northing

! +
! +- Polar Stereographic Projection
! +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! +
      tF  = tan (pi/4 + latF  /2) / 
     .      ( (1 + ex*sin(latF  )) / (1 - ex*sin(latF  )) )**(ex/2)
      mF  = cos(latF) / (1 - ex**2 * sin(latF)**2)**0.5
      k0  = mF*( (1+ex)**(1+ex) * (1-ex)**(1-ex) )**0.5 / (2*tF)
      
      t   = tan (pi/4 + latrad/2) / 
     .      ( (1 + ex*sin(latrad)) / (1 - ex*sin(latrad)) )**(ex/2)
      rho = 2*aa*k0*t / ( (1+ex)**(1+ex) * (1-ex)**(1-ex) )**0.5
      
      E = FE + rho*sin (lonrad - lon0)
      N = FN + rho*cos (lonrad - lon0)

      RETURN
      End Subroutine StereoSouth_inverse
C +------------------------------------------------------------------------+
