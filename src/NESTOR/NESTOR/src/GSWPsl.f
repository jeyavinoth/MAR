C   +-------------------------------------------------------------------+
C   |  Subroutine GSWPsl                            April 2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : - NST__x : NST grid, longitude (degrees)                  |
C   | ^^^^^^^ - NST__y : NST grid, latitude  (degrees)                  |
C   |                                                                   |
C   | Output: - NSTdsa : soil albedo                                    |
C   | ^^^^^^^ - NSTtex : soil texture (fine, medium, rough)             |
C   |                                                                   |
C   |           from GSWP data set (http://grads.iges.org/gswp/)        |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE GSWPsl


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
      
      INTEGER error,itmp1,itmp2,Nlon,Nlat,ii,jj,i,j,inc,Nlon1,Nlat1
      INTEGER iimin,jjmin 

      PARAMETER(Nlon=360,Nlat=150,Nlon1=320,Nlat1=320)

      REAL    rtmp1,rtmp2,distmin,
     .        GSWPalb(NLON,NLAT),GSWPsol(NLON,NLAT),GSWPmsk(NLON,NLAT),
     .        GSWPlon(NLON),GSWPlat(NLAT),
     .        NOAAalb(NLON1,NLAT1),NOAAlon(NLON1),NOAAlat(NLAT1)
      
C +---1. Read cdf file for soil parameters : texture and bare soil albedo
C +   ===================================================================

C +---1.1 GSWP (GSWP-SOIL.nc)
C     -----------------------

      error = nf_open('input/SOIL/GSWP-SOIL.nc',nf_nowrite,inc)

      IF (error.ne.nf_noerr) THEN
       write(6,*) '+++++++++++++++++++++++++++++++++'
       write(6,*) 'Routine GSWPsl.f -----> Warning !!!'
       write(6,*) 'File GSWP-SOIL.nc not provided'
       write(6,*) 'Check the directory input/SOIL/'
       write(6,*) 'NESTOR stopped NOW !!!'
       write(6,*) '+++++++++++++++++++++++++++++++++'
       stop
      ENDIF
      
      error = nf_inq_varid(inc,'LON'        ,itmp1)
      error = nf_get_var_real(inc           ,itmp1,GSWPlon)
      error = nf_inq_varid(inc,'LAT'        ,itmp1)
      error = nf_get_var_real(inc           ,itmp1,GSWPlat)
      error = nf_inq_varid(inc,'ALBEDO_SOIL',itmp1)
      error = nf_get_var_real(inc           ,itmp1,GSWPalb)
      error = nf_inq_varid(inc,'SOILCLASS'  ,itmp1)
      error = nf_get_var_real(inc           ,itmp1,GSWPsol)
      error = nf_inq_varid(inc,'LANDMASK'   ,itmp1)
      error = nf_get_var_real(inc           ,itmp1,GSWPmsk)
      error = nf_close(inc)

C +---1.2 NOAA (AFRmax-alb.nc)
C     ------------------------

      error = nf_open('input/SOIL/AFRmax-alb.nc',nf_nowrite,inc)
      
      IF (region.eq.'AFW'.and.error.ne.nf_noerr) THEN
       write(6,*) '+++++++++++++++++++++++++++++++++'
       write(6,*) 'Routine GSWPsl.f -----> Warning !!!'
       write(6,*) 'File AFRmax-alb.nc not provided'
       write(6,*) 'Check the directory input/SOIL/'
       write(6,*) 'NESTOR stopped NOW !!!'
       write(6,*) '+++++++++++++++++++++++++++++++++'
       stop
      ENDIF
      
      error = nf_inq_varid(inc,'lon',itmp1)
      error = nf_get_var_real(inc,   itmp1,NOAAlon)
      error = nf_inq_varid(inc,'lat',itmp1)
      error = nf_get_var_real(inc,   itmp1,NOAAlat)
      error = nf_inq_varid(inc,'alb',itmp1)
      error = nf_get_var_real(inc,   itmp1,NOAAalb)
      error = nf_close(inc)      

C +---2. GSWP grid  --->  NST grid
C +   ============================

      DO j=1,my    ! Loop for each NST grid point
      DO i=1,mx    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C +---2.1 Sea and Sea Ice
C     -------------------
 
       IF (NSTsol(i,j).le.2) then
        NSTtex(i,j) = 0
        NSTdsa(i,j) = 0.15 
        IF (region.eq."GRD".or.region.eq."ANT") NSTdsa(i,j) = 0.20
       ENDIF

C +---2.2 Snow - Ice
C     --------------   
 
       IF (NSTsol(i,j).eq.3) then
        NSTtex(i,j) = 3
        NSTdsa(i,j) = 0.85 
       ENDIF

C +---2.3 Soil - Tundra
C     -----------------
 
       IF (NSTsol(i,j).ge.4) then

        itmp1=0
        rtmp1=0
        rtmp2=0
        
        distmin=10000 
   
        DO ii = 1,NLON
         DO jj = 1,NLAT

          if(abs(NST__x(i,j)-GSWPlon(ii))+
     .       abs(NST__y(i,j)-GSWPlat(jj))<distmin.and.
     .       GSWPmsk(ii,jj)              .ne.0.0) then

            distmin=
     .       abs(NST__x(i,j)-GSWPlon(ii))+
     .       abs(NST__y(i,j)-GSWPlat(jj))
            iimin  = ii
            jjmin  = jj
            itmp1  = 1             + itmp1 
           endif

c          if(abs(NST__x(i,j)-GSWPlon(ii)).le.0.6 .and.
c     .       abs(NST__y(i,j)-GSWPlat(jj)).le.0.6 .and.
c     .       GSWPmsk(ii,jj)              .ne.0.0) then
c                              
c          ! GSWP resolution = 1°
c
c           rtmp1 = GSWPalb(ii,jj) + rtmp1   
c           rtmp2 = GSWPsol(ii,jj) + rtmp2  
c           itmp1 = 1              + itmp1
c
c          ENDIF
 
         ENDDO
        ENDDO
 
        IF (itmp1.gt.0) THEN
         NSTdsa(i,j) = REAL(rtmp1/itmp1)
         NSTtex(i,j) = INT (rtmp2/itmp1)
         NSTdsa(i,j) = GSWPalb(iimin,jjmin)
         NSTtex(i,j) = GSWPsol(iimin,jjmin)
        ELSE
         NSTtex(i,j) = 5
         NSTdsa(i,j) = 0.25
        ENDIF

       ENDIF

C +---2.3.1 Special Albedo for AFW simulation

C XF Jan 2014: albedo too low in the "congo" basin

c       IF (region.eq.'AFW'.and.NSTsol(i,j).ge.4) THEN
c
c        itmp1=0
c        rtmp1=0
c
c        DO ii = 1,NLON1
c         DO jj = 1,NLAT1
c 
c          if(abs(NST__x(i,j)-NOAAlon(ii))  .le.  0.30.and.
c     .       abs(NST__y(i,j)-NOAAlat(jj))  .le.  0.30.and.
c     .                       NOAAalb(ii,jj).ne.-99.0 ) then
c 
c          ! NOAA resolution = 0.25 deg
c 
c           rtmp1 = NOAAalb(ii,jj) + rtmp1
c           itmp1 = 1              + itmp1
c 
c          ENDIF
c 
c         ENDDO
c        ENDDO
c 
c        IF (itmp1.gt.0) THEN
c         NSTdsa(i,j) = REAL(rtmp1/itmp1)/100.
c        ELSE
c         NSTdsa(i,j) = 0.25
c        ENDIF
c
c       ENDIF

        IF (region.eq.'AFW'.and.NSTsol(i,j).ge.4) THEN
        NSTdsa(i,j) = max(0.1,min(0.45,NSTdsa(i,j))) 
        ENDIF

C +---2.3.2 Special Texture/Albedo for GRD Simulation

       IF (region.eq."GRD".and.NSTsol(i,j).ge.4) THEN
        NSTdsa(i,j) = 0.25 
        NSTtex(i,j) = 2
       ENDIF 

C +---2.3.3 Max/Min of Texture/Albedo 

       IF (NSTsol(i,j).ge.4) THEN
        NSTdsa(i,j) = max(0.05,min(0.5,NSTdsa(i,j))) 
        NSTalb(i,j) =                  NSTdsa(i,j)
        NSTtex(i,j) = max(1   ,min(12 ,NSTtex(i,j)))
       ENDIF
    
      ENDDO        !  Loop for i (NST grid)
      ENDDO        !  Loop for j (NST grid)

      RETURN
      END
