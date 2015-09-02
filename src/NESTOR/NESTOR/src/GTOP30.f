C   +-------------------------------------------------------------------+
C   |  Subroutine GTOP30                             June 2003  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : NST__x : longitude (degree) of the NST grid               |
C   | ^^^^^^^ NST__y : latitude  (degree) of the NST grid               |
C   |                                                                   |
C   | Output: NST_sh: surface elevation                                 |
C   | ^^^^^^^ NSTsol: land (4) / sea (1) mask                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE GTOP30


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'LOCfil.inc'
      INCLUDE 'NetCDF.inc'

      INTEGER NBfile,i,j,n,ns,nbchar,nfile,i_cent,j_cent,i1,
     .        i2,j1,j2,TOP_id,Rcode,start(3),count(3),nsamp,
     .        k,l,G_nx,G_ny,num,ii,jj

      PARAMETER (NBfile= 33)

      INTEGER G_lon1(NBfile),G_lon2(NBfile),G_lat1(NBfile),
     .        G_lat2(NBfile),Gridx (NBfile),Gridy (NBfile),
     .        fID   (NBfile),GTOPOnum(NBfile)

      INTEGER*2 TOPhgt(1,1),zero

      REAL MINlon,MAXlon,MINlat,MAXlat,G_reso,degrad,dx,dy,
     .     dx1,dx2,dy1,dy2,altsum,G_dx,G_dy,AUXlon,AUXlat,
     .     AUXlo1,AUXlo2,AUXla1,AUXla2

      CHARACTER*7  NAMfil   (NBfile)
      CHARACTER*60 GTOPOdir
      CHARACTER*80 GTOPOfile(NBfile)

      LOGICAL Vtrue


C +---Data
C +   ----

      DATA G_reso  / 0.00833333333333d0 /
      DATA degrad  / 1.745329252d-2     /
      DATA start   / 1,1,1              /
      DATA count   / 0,0,0              /
      DATA zero    / 0                  /
      DATA Vtrue   / .true.             /


C +---Screen message
C +   --------------

      write(6,*) 'Topography : GTOPO data set (30 secondes)'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'


C +---Input directory
C +   ---------------

      nbchar  =1
      GTOPOdir=GTOPO_dir

      DO i=1,60
       IF (GTOPOdir(i:i).ne.' ') nbchar=i
      ENDDO


C +---Data : informations about input data files
C +   ------------------------------------------

      OPEN (unit=25,status='old',
     .      file=GTOPOdir(1:nbchar)//'File_List')

      DO n=1,6
       READ (25,*)
      ENDDO

      DO n=1,NBfile
       READ (25,250) NAMfil(n),G_lon1(n),G_lon2(n),G_lat1(n),
     .               G_lat2(n),Gridy (n),Gridx (n)
250    FORMAT (A7,5X,I5,5X,I5,6X,I4,6X,I4,6X,I5,5X,I5)
      ENDDO

      CLOSE (unit=25)


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


C +---Selection of the data files including the NST domain
C +   ----------------------------------------------------

      write(6,*) 'Files to be read : '
      ns=0

      DO n=1,NBfile

       GTOPOfile(n)= GTOPOdir(1:nbchar)//NAMfil(n)//'.nc'
       
C +... Check if point (1,1) is included in simulation domain
       IF ( ((REAL(G_lon1(n)).ge.MINlon) .and.
     .       (REAL(G_lon1(n)).le.MAXlon) .and.
     .       (REAL(G_lat1(n)).ge.MINlat) .and.
     .       (REAL(G_lat1(n)).le.MAXlat))     .or.
C +... Check if point (1,2) is included in simulation domain
     .      ((REAL(G_lon1(n)).ge.MINlon) .and.
     .       (REAL(G_lon1(n)).le.MAXlon) .and.
     .       (REAL(G_lat2(n)).ge.MINlat) .and.
     .       (REAL(G_lat2(n)).le.MAXlat))     .or.
C +... Check if point (2,1) is included in simulation domain
     .      ((REAL(G_lon2(n)).ge.MINlon) .and.
     .       (REAL(G_lon2(n)).le.MAXlon) .and.
     .       (REAL(G_lat1(n)).ge.MINlat) .and.
     .       (REAL(G_lat1(n)).le.MAXlat))     .or.
C +... Check if point (2,2) is included in simulation domain
     .      ((REAL(G_lon2(n)).ge.MINlon) .and.
     .       (REAL(G_lon2(n)).le.MAXlon) .and.
     .       (REAL(G_lat2(n)).ge.MINlat) .and.
     .       (REAL(G_lat2(n)).le.MAXlat))     .or.
C +... Check if point (1,1) is included in GTOPO domain
     .      ((MINlon.ge.REAL(G_lon1(n))) .and.
     .       (MINlon.le.REAL(G_lon2(n))) .and.
     .       (MINlat.ge.REAL(G_lat1(n))) .and.
     .       (MINlat.le.REAL(G_lat2(n))))    .or.
C +... Check if point (1,2) is included in GTOPO domain
     .      ((MINlon.ge.REAL(G_lon1(n))) .and.
     .       (MINlon.le.REAL(G_lon2(n))) .and.
     .       (MAXlat.ge.REAL(G_lat1(n))) .and.
     .       (MAXlat.le.REAL(G_lat2(n))))    .or.
C +... Check if point (2,1) is included in GTOPO domain
     .      ((MAXlon.ge.REAL(G_lon1(n))) .and.
     .       (MAXlon.le.REAL(G_lon2(n))) .and.
     .       (MINlat.ge.REAL(G_lat1(n))) .and.
     .       (MINlat.le.REAL(G_lat2(n))))    .or.
C +... Check if point (2,2) is included in GTOPO domain
     .      ((MAXlon.ge.REAL(G_lon1(n))) .and.
     .       (MAXlon.le.REAL(G_lon2(n))) .and.
     .       (MAXlat.ge.REAL(G_lat1(n))) .and.
     .       (MAXlat.le.REAL(G_lat2(n))))    .or.
C +... Check if the simulation domain is fully included in the file
     .      ((REAL(G_lon1(n)).ge.MINlon) .and.      
     .       (REAL(G_lon2(n)).le.MAXlon) .and.      
     .       (REAL(G_lat1(n)).le.MINlat) .and.      
     .       (REAL(G_lat2(n)).ge.MAXlat))    .or.
     .      ((REAL(G_lon1(n)).le.MINlon) .and.
     .       (REAL(G_lon2(n)).ge.MAXlon) .and.
     .       (REAL(G_lat1(n)).le.MINlat) .and.
     .       (REAL(G_lat2(n)).ge.MAXlat)) )   THEN
        ns           = ns+1
        GTOPOnum (ns)= n
        write(6,*)     GTOPOfile(n)
       ENDIF

      ENDDO

      write(6,*)


C +---Open required Netcdf data files
C +   -------------------------------

      DO n=1,ns
       num=GTOPOnum(n)
       fID(num)=NCOPN(GTOPOfile(num),NCNOWRIT,Rcode)
      ENDDO


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO j=1,my     ! Loop on horizontal NST grid points
      DO i=1,mx

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Choice of the appropriate data file
C +   -----------------------------------

       nfile=0

       AUXlon = NST__x(i,j)
       AUXlat = NST__y(i,j)
C +         ******
       CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +         ******

       DO n=1,ns
        num = GTOPOnum(n)
        IF ( (AUXlon.ge.REAL(G_lon1(num))) .and.
     .       (AUXlon.lt.REAL(G_lon2(num))) .and.
     .       (AUXlat.ge.REAL(G_lat1(num))) .and.
     .       (AUXlat.lt.REAL(G_lat2(num))) ) nfile=num
       ENDDO

       IF (nfile.eq.0) THEN
        WRITE (6,*) 'No file found for the NST grid point ',i,j
        WRITE (6,*) 'i.e. lat= ',NST__y(i,j),', long= ',NST__x(i,j)
        WRITE (6,*) 'STOP in GTOP30.'
        STOP
       ENDIF


C +---Search for the closest point in data file
C +   -----------------------------------------

       i_cent=NINT((AUXlon-REAL(G_lon1(nfile)))/G_reso)+1
       j_cent=NINT((AUXlat-REAL(G_lat1(nfile)))/G_reso)+1


C +---Compute the resolution of the considered NST cell
C +   -------------------------------------------------

       ii=MAX(i,2)
       jj=MAX(j,2)

       AUXlo1 = NST__x(ii  ,jj  )
       AUXla1 = NST__y(ii  ,jj  )
       AUXlo2 = NST__x(ii-1,jj-1)
       AUXla2 = NST__y(ii-1,jj-1)
C +         ******
       CALL SPHERC (Vtrue,AUXlo1,AUXla1)
       CALL SPHERC (Vtrue,AUXlo2,AUXla2)
C +         ******

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
       i2=MIN(i2,Gridx(nfile))
       j1=MAX(j1,1)
       j2=MIN(j2,Gridy(nfile))


C +---Read subset of data
C +   -------------------

       TOP_id=NCVID(fID(nfile),'topo',Rcode)

       nsamp =0
       altsum=0.

       DO l=j1,j2   ! Loop on data grid points
       DO k=i1,i2   ! contained in the (i,j) NST cell

        start(1)=k
        start(2)=l
        count(1)=1
        count(2)=1

        CALL NCVGT(fID(nfile),TOP_id,start,count,TOPhgt,Rcode)

        IF (TOPhgt(1,1).ne.9999) THEN
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

      ENDDO    ! Loop on NST grid points
      ENDDO

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
      RETURN
      END
