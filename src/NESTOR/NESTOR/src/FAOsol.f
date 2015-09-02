C   +-------------------------------------------------------------------+
C   |  Subroutine FAOsol                         February 2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : - NST__x : NST grid, longitude (degrees)                  |
C   | ^^^^^^^ - NST__y : NST grid, latitude  (degrees)                  |
C   |                                                                   |
C   | Output: - NSTdsa : soil albedo                                    |
C   | ^^^^^^^ - NSTtex : soil texture (fine, medium, rough)             |
C   |                                                                   |
C   | This routine reads the files :                                    |
C   | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                      |
C   | - TEXUNIT.ASC : soil type characteristics given in the FAO        |
C   |                 classification ;                                  |
C   | - FAO_SOIL.nc : soil types over the globe                         |
C   |                 --> domain : -180 -> 180 (lon.), -90  -> 90 (lat) |
C   |                 --> resolution : 5 minutes (8 to 10 km)           |
C   |                                                                   |
C   | FAO soil data file contains 9 fields :                            |
C   | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                              |
C   | - Soil Texture Map :                                              |
C   |   ----------------                                                |
C   |   . 1 : % of cell containing rough  texture ;                     |
C   |   . 2 : % of cell containing medium texture ;                     |
C   |   . 3 : % of cell containing fine   texture ;                     |
C   |   . 4 : % of cell containing dominant texture ;                   |
C   |                                                                   |
C   | - Soil Moisture Storage Capacity Map :                            |
C   |   ----------------------------------                              |
C   |   . 5 : phase ;                                                   |
C   |   . 6 : % phase ;                                                 |
C   |                                                                   |
C   | - Soil Map :                                                      |
C   |   --------                                                        |
C   |   . 7 : slope < 8 % ;                                             |
C   |   . 8 : 8 % < slope < 30 % ;                                      |
C   |   . 9 : slope > 30 % .                                            |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE FAOsol


      IMPLICIT NONE


C +---Netcdf specifications
C +   ---------------------

      INCLUDE 'NetCDF.inc'


C +---General and local variables
C +   ---------------------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'LOCfil.inc'

      INTEGER i,j,num,FAO_nx,FAO_ny,nline,ntype,ix,iy,i2,j2,compt,
     .        stype,ncl,FAO_ID,i1lat,i2lat,i1lon,i2lon,imlat,imlon,
     .        start(3),count(3),nbi,nbj,ncid,Rcode,FAOmaxy,nbchar

      REAL coarfr,amedfr,finefr,phase,plat,prphase,ondule,pente,dx,
     .     dy,FAO_px,FAO_py,FAOres,pi,degrad,AUXlon,AUXlat,MINlat,
     .     MAXlat,FAO_dx,FAO_dy,dx1,dx2,dy1,dy2
    
      CHARACTER*20 name

      CHARACTER*80 FAO_text,FAO_file


C +---General parameters of FAO file
C +   ------------------------------

      PARAMETER(FAO_nx=4320,FAO_ny=1300,FAOmaxy=2160)
C +...          ^ size of the large-scale arrays  

      PARAMETER(ntype=7000,nline=4931)
C +...          ^ number of soil types in FAO classification
C +...                     ^ number of lines to be read in EUROSOIL

      PARAMETER(FAOres=1./12.)
C +...          ^ horizontal resolution of the data

      PARAMETER(FAO_px=-180.,FAO_py=-90.)
C +...          ^ initial grid point in longitude and latitude

      INTEGER*2 FAOsoil(FAO_nx,FAO_ny)

      INTEGER FAOtext(0:3),texture(ntype)

      REAL    FAOlon(FAO_nx),FAOlat(FAO_ny)

      LOGICAL Vtrue


C +---Data
C +   ----

      DATA start  / 1,1,1/
      DATA count  / 0,0,0/
      DATA Vtrue  /.true./

      pi    = 3.141592653589793238462643d0
      degrad= pi / 180.d0


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Read the FAO CLASSIFICATION OF SOIL TYPES
C +   =========================================

      nbchar=1

      DO i=1,60
       IF (FAO_dir(i:i).ne.' ') nbchar=i
      ENDDO

      FAO_text = FAO_dir(1:nbchar) // 'TEXUNIT.ASC'
      OPEN(unit=10,status='old',file=FAO_text)

      DO i=1,ntype
       texture(i)=0.
      ENDDO

      DO i=1,nline

       READ(10,100) num,name,coarfr,amedfr,finefr,phase,prphase,
     .              plat,ondule,pente

100    FORMAT(i4,1x,a20,3f6.0,f4.0,f9.3,3f8.0)

C +    COARSE TEXTURE --> SAND --> num.= 1
       IF (coarfr.ge.amedfr.and.coarfr.ge.finefr) texture(num)=1
C +    MEDIUM TEXTURE --> LOAM --> num.= 2
       IF (amedfr.ge.coarfr.and.amedfr.ge.finefr) texture(num)=2
C +    FINE   TEXTURE --> CLAY --> num.= 3
       IF (finefr.ge.coarfr.and.finefr.ge.amedfr) texture(num)=3

      ENDDO

      CLOSE(unit=10)


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Read FAO SOIL TYPES
C +   ===================


      write(6,*) 'Soil types : FAO data set (5 minutes)'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*)


C +---Open NetCDF file
C +   ----------------

      FAO_file = FAO_dir(1:nbchar) // 'FAO_SOIL.nc'

      ncid   = NCOPN(FAO_file,NCNOWRIT,Rcode)
      FAO_ID = NCVID(ncid,'FAO_SOIL'  ,Rcode)


C +---Compute the extension of the sub-domain to be read
C     --------------------------------------------------

      AUXlon = NST__x(1,1)
      AUXlat = NST__y(1,1)
C +        ******
      CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +        ******
      MINlat = AUXlat
      MAXlat = AUXlat
      DO j=1,my
      DO i=1,mx
       AUXlon = NST__x(i,j)
       AUXlat = NST__y(i,j)
C +         ****** 
       CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +         ******
       MINlat = min(AUXlat,MINlat)
       MAXlat = max(AUXlat,MAXlat)
      ENDDO
      ENDDO


C +---Define extraction zone
C +   ----------------------

      i1lat = AINT ( (MINlat-(-90.))/FAOres )
      i2lat = AINT ( (MAXlat-(-90.))/FAOres )
      imlat = i2lat - i1lat + 1

      i1lon = 1
      i2lon = FAO_nx
      imlon = i2lon - i1lon + 1

      IF (imlat.ge.FAO_ny) THEN
       write(*,*) 'Extent of the simulation domain in latitude'
       write(*,*) 'is too large. Please choose a larger value '
       write(*,*) 'for the FAO_ny parameter.           - STOP '
       STOP
      ENDIF

      i1lat = i1lat + (i2lat-i1lat)/2 - FAO_ny/2
      i2lat = i1lat + FAO_ny - 1

      IF (i1lat.lt.1) THEN
       i1lat=1
       i2lat=FAO_ny
      ENDIF

      IF (i2lat.gt.FAOmaxy) THEN
       i1lat=FAOmaxy-FAO_ny+1
       i2lat=FAOmaxy
      ENDIF


C +---Define the FAO grid
C +   -------------------

      DO i=i1lon,i2lon
       FAOlon(i)=FAO_px+real(i-1)*FAOres+FAOres/2.
      ENDDO

      DO j=i1lat,i2lat
       FAOlat(j-i1lat+1)=FAO_py+real(j-1)*FAOres+FAOres/2.
      ENDDO


C +---Read values of the variables for the sub-domain
C +   -----------------------------------------------

      start(1)=i1lon
      start(2)=i1lat
      count(1)=FAO_nx
      count(2)=FAO_ny
      CALL NCVGT(ncid,FAO_ID,start,count,FAOsoil,Rcode)


C +---Close Netcdf data file
C +   ----------------------

      CALL NCCLOS(ncid,Rcode)


C +   ###################################################################

      DO j=1,my    ! Loop for each NST grid point
      DO i=1,mx    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C +   ###################################################################


C +---FAO grid  --->  NST grid
C +   ========================


C +---Nearest FAO grid point
C +   ----------------------

       AUXlon = NST__x(i,j)
       AUXlat = NST__y(i,j)
C +         ****** 
       CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +         ******
       ix=NINT((AUXlon-FAOlon(1))/FAOres)+1
       iy=NINT((AUXlat-FAOlat(1))/FAOres)+1


C +---Correction to avoid limits in the representation of numbers
C +   -----------------------------------------------------------

       IF (ix.eq.0) THEN
        IF (((AUXlon-FAOlon(1))/FAOres).gt.(-0.51)) ix=1
       ENDIF

       IF (iy.eq.0) THEN
        IF (((AUXlat-FAOlat(1))/FAOres).gt.(-0.51)) iy=1
       ENDIF

       IF (ix.eq.(FAO_nx+1)) THEN
        IF (((AUXlon-FAOlon(1))/FAOres).lt.(real(FAO_nx-1)+0.51))
     .   ix=FAO_nx
       ENDIF

       IF (iy.eq.(FAO_ny+1)) THEN
        IF (((AUXlat-FAOlat(1))/FAOres).lt.(real(FAO_ny-1)+0.51))
     .   iy=FAO_ny
       ENDIF


C +---Check if NST latitudes/longitudes are valid
C +   -------------------------------------------

       IF ( (ix.lt.1).or.(ix.gt.FAO_nx).or.
     .      (iy.lt.1).or.(iy.gt.FAO_ny) ) THEN
        WRITE (6,*) 'NST grid point (i,j)= ',i,j,' out of global'
        WRITE (6,*) 'FAO domain.'
        WRITE (6,*) 'STOP in FAOsol.'
        STOP
       ENDIF


C +---Compute the resolution of the considered NST cell
C +   -------------------------------------------------

       IF (NST_dx.gt.0.01) THEN

        dx =NST_dx
        dy =NST_dx

       ELSE

        dx1=(NST__x(i,j)-NST__x(i-1,j))*111111.
     .                                 *COS(NST__y(i,j)*degrad)
        dx2=(NST__y(i,j)-NST__y(i-1,j))*111111.
        dx =SQRT(dx1*dx1+dx2*dx2)

        dy1=(NST__x(i,j)-NST__x(i,j-1))*111111.
     .                                 *COS(NST__y(i,j)*degrad)
        dy2=(NST__y(i,j)-NST__y(i,j-1))*111111.
        dy =SQRT(dy1*dy1+dy2*dy2)

       ENDIF


C +---Define the data points to be read around (i_cent,j_cent)
C +   --------------------------------------------------------

       FAO_dx = FAOres*111111.*COS(NST__y(i,j)*degrad)
       FAO_dy = FAOres*111111.

       nbi=NINT(dx/FAO_dx/2.)-1
       nbj=NINT(dy/FAO_dy/2.)-1

       nbi=MAX(nbi,0)
       nbj=MAX(nbj,0)


C +---Examine all FAO grid points contained in a given NST cell
C +   ---------------------------------------------------------

       FAOtext(0)=0
       FAOtext(1)=0
       FAOtext(2)=0
       FAOtext(3)=0

       DO j2=MAX(1,iy-nbj),MIN(FAO_ny,iy+nbj)
       DO i2=MAX(1,ix-nbi),MIN(FAO_nx,ix+nbi)

        IF (FAOsoil(i2,j2).ne.0) THEN
         stype=texture(FAOsoil(i2,j2))
        ELSE
         stype=0
        ENDIF

C +...  stype = 0 : water
C +...          1 : coarse texture
C +...          2 : medium texture
C +...          3 : fine   texture

        IF (stype.ge.0.and.stype.le.3)
     .   FAOtext(stype)=FAOtext(stype)+1

       ENDDO
       ENDDO


C +---Determination of the dominant soil type
C +   ---------------------------------------

       compt=0
       stype=2

       DO ncl=0,3
        IF (FAOtext(ncl).gt.compt) THEN
         compt=FAOtext(ncl)
         stype=ncl
        ENDIF
       ENDDO


C +---Output : NSTsol and NSTtex
C +   --------------------------

C +---Sea and Sea Ice
 
       IF (NSTsol(i,j).le.2 ) then
        NSTtex(i,j) = 0
        NSTdsa(i,j) = 0.15 
        IF (region.eq."GRD".or.region.eq."ANT") NSTdsa(i,j) = 0.20
       ENDIF

C +---Snow - Ice
 
       IF (NSTsol(i,j).eq.3 ) then
        NSTtex(i,j) = 3
        NSTdsa(i,j) = 0.85 
       ENDIF

C +---Soil - Tundra

       IF (NSTsol(i,j).ge.4) THEN

        IF (stype.ne.0) THEN
         NSTtex(i,j)=stype
        ELSE
         NSTtex(i,j)=2
        ENDIF

        IF (NSTtex(i,j).eq.1) THEN
         NSTdsa(i,j) = 0.40
C +...   Dry Quartz Sand (Deardorff 1978 JGR p.1891)
        ELSE IF (NSTtex(i,j).eq.3) THEN
         NSTdsa(i,j) = 0.15
C +...   Clay Pasture    (Deardorff 1978 JGR p.1891)
        ELSE
         NSTdsa(i,j) = 0.25
C +...   O'Neill average (Deardorff 1978 JGR p.1891)
        ENDIF

       ENDIF

C +---Special Texture for Polar Simulation
C +   ------------------------------------

       IF (region.eq."GRD".or.region.eq."ANT") THEN
       
        if (i.eq.2.and.j.eq.2) then
         write(6,*) "Special Texture for Polar region" 
         write(6,*)
        endif
 
        IF (NSTsol(i,j).le.2 ) NSTtex(i,j) = 0 ! sea
        IF (NSTsol(i,j).eq.3 ) NSTtex(i,j) = 3 ! ice
        IF (NSTsol(i,j).ge.4 ) NSTtex(i,j) = 2 ! tundra              
       
       ENDIF
    
C +   ###################################################################
    
      ENDDO        !  Loop for i (NST grid)
      ENDDO        !  Loop for j (NST grid)
    
C +   ###################################################################


      RETURN
      END
