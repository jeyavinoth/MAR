C   +-------------------------------------------------------------------+
C   |  Subroutine BETOPO                        September 2002  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : NST__x : longitude (degree) of the NST grid               |
C   | ^^^^^^^ NST__y : latitude  (degree) of the NST grid               |
C   |                                                                   |
C   | Output: NST_sh: surface elevation                                 |
C   | ^^^^^^^ NSTsol: land (4) / sea (1) mask                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE BETOPO (NST__x,NST__y,NST_sh,NSTsol)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'LOCfil.inc'
      INCLUDE 'NetCDF.inc'

      INTEGER i,j,n,nbchar,i_cent,j_cent,i1,
     .        i2,j1,j2,TOP_id,Rcode,start(3),count(3),nsamp,
     .        k,l,G_nx,G_ny,ii,jj,fID,N_lon,N_lat

      INTEGER*2 TOPhgt(1,1),zero

      REAL    G_reso,degrad,dx,dy,dx1,dx2,dy1,dy2,altsum,
     .        G_dx,G_dy,AUXlon,AUXlat,AUXlo1,AUXlo2,AUXla1,
     .        AUXla2,G_lon1,G_lon2,G_lat1,G_lat2

      REAL    NST__x(mx,my),NST__y(mx,my),NST_sh(mx,my)

      INTEGER NSTsol(mx,my)

      LOGICAL Vtrue


C +---Data
C +   ----

      DATA degrad  / 1.745329252d-2     /
      DATA start   / 1,1,1              /
      DATA count   / 0,0,0              /
      DATA zero    / 0                  /
      DATA Vtrue   / .true.             /


C +---Info about data file
C +   --------------------

      G_lon1 =  2.50     !  South-West corner (longitude)
      G_lat1 = 49.25     !  South-West corner (latitude )

      G_lon2 =  6.50     !  North-East corner (longitude)
      G_lat2 = 51.50     !  North-East corner (latitude )

      G_reso = 1./3600.  !  Resolution : about 30 m

      N_lon  = 14401     !  Number of grid points in longitude
      N_lat  =  8101     !  Number of grid points in latitude


C +---Screen message
C +   --------------

      write(6,*) 'Topography of Belgium (30 m resolution)'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'


C +---Input directory
C +   ---------------

      nbchar = 1

      DO i=1,60
       IF (BTOPO_dir(i:i).ne.' ') nbchar=i
      ENDDO


C +---Open the Netcdf data file
C +   -------------------------

      fID=NCOPN(BTOPO_dir(1:nbchar)//'B_TOPO.nc',NCNOWRIT,Rcode)


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO j=1,my     ! Loop on horizontal NST grid points
      DO i=1,mx

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Choice of the appropriate data file
C +   -----------------------------------

       AUXlon = NST__x(i,j)
       AUXlat = NST__y(i,j)
C +         ******
       CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +         ******


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Check if the grid point is located in the data domain
C +   -----------------------------------------------------

       IF ((AUXlon.ge.G_lon1).and.(AUXlon.le.G_lon2).and.
     .     (AUXlat.ge.G_lat1).and.(AUXlat.le.G_lat2)) THEN

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Search for the closest point in data file
C +   -----------------------------------------

        i_cent=NINT((AUXlon-REAL(G_lon1))/G_reso)+1
        j_cent=NINT((AUXlat-REAL(G_lat1))/G_reso)+1


C +---Compute the resolution of the considered NST cell
C +   -------------------------------------------------

        ii = MAX(i,2)
        jj = MAX(j,2)

        AUXlo1 = NST__x(ii  ,jj  )
        AUXla1 = NST__y(ii  ,jj  )
        AUXlo2 = NST__x(ii-1,jj-1)
        AUXla2 = NST__y(ii-1,jj-1)
C +          ******
        CALL SPHERC (Vtrue,AUXlo1,AUXla1)
        CALL SPHERC (Vtrue,AUXlo2,AUXla2)
C +          ******

        dx=ABS(AUXlo1-AUXlo2)*111111.*COS(AUXla1*degrad)
        dy=ABS(AUXla1-AUXla2)*111111.


C +---Define the data points to be read around (i_cent,j_cent)
C +   --------------------------------------------------------

        G_dx = G_reso*111111.*COS(AUXla1*degrad)
        G_dy = G_reso*111111.

        G_nx=NINT(dx/G_dx/2.)-1
        G_ny=NINT(dy/G_dy/2.)-1

        G_nx=MAX(G_nx,0)
        G_ny=MAX(G_ny,0)

        i1=i_cent-G_nx
        i2=i_cent+G_nx
        j1=j_cent-G_ny
        j2=j_cent+G_ny

        i1=MAX(i1,1)
        i2=MIN(i2,N_lon)
        j1=MAX(j1,1)
        j2=MIN(j2,N_lat)


C +---Read subset of data
C +   -------------------

        TOP_id=NCVID(fID,'topo',Rcode)

        nsamp =0
        altsum=0.

        DO l=j1,j2   ! Loop on data grid points
        DO k=i1,i2   ! contained in the (i,j) NST cell

         start(1)=k
         start(2)=l
         count(1)=1
         count(2)=1

         CALL NCVGT(fID,TOP_id,start,count,TOPhgt,Rcode)

         IF (TOPhgt(1,1).ge.0.0.and.TOPhgt(1,1).le.2000.0) THEN
          altsum=altsum+MAX(zero,TOPhgt(1,1))
          nsamp =nsamp+1
         ENDIF

        ENDDO
        ENDDO


C +---Final computation of the topography at (i,j) location
C +   -----------------------------------------------------

        IF (nsamp.ne.0) THEN
         NST_sh(i,j)=altsum/real(nsamp)
        ENDIF


C +---Distinction between land and sea (further refined)
C +   --------------------------------

        IF (NST_sh(i,j).lt.0.01) THEN
         NSTsol(i,j)=1
        ELSE
         NSTsol(i,j)=4
        ENDIF


C +---No atmosphere below sea level...
C +   --------------------------------

        IF (NST_sh(i,j).lt.0.0) THEN
         NST_sh(i,j)= 0.0
        ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ENDIF   ! Grid point (i,j) is in Belgium

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDDO    ! Loop on NST grid points
      ENDDO

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
