C   +-------------------------------------------------------------------+
C   |  Subroutine GLOfrc                      18 March    2009  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | GLOfrc read NDVI index (-> max. fract. of vegetation) over Africa |
C   | and Europe.                                                       |
C   |                                                                   |
C   | Input : - NST__x, NST__y : NST grid coordinates (lat./long.)      |
C   | ^^^^^^^ - NSTsol : soil type                                      |
C   |         - NSTveg : vegetation type (IGBP classification)          |
C   |         - NSTvfr : fraction of vegetation in the grid cell (IGBP) |
C   |         - NSTsvt : vegetation type (SVAT classification)          |
C   |         - NSTsfr : fraction of vegetation in the grid cell (SVAT) |
C   |                                                                   |
C   | Output: - NSTveg : vegetation type (IGBP classification)          |
C   | ^^^^^^^ - NSTvfr : fraction of vegetation in the grid cell (IGBP) |
C   |         - NSTsvt : vegetation type (SVATclassification)           |
C   |         - NSTsfr : fraction of vegetation in the grid cell (SVAT) |
C   |         - NSTfrc : fraction of vegetation cover (from NDVI index) |
C   |         - NSTdv1 : minimum NDVI index over a period of one year   |
C   |         - NSTdv2 : maximum NDVI index over a period of one year   |
C   |                                                                   |
C   | Remark: Note that NSTveg = -1 (IGBP) or NSTsvt = 0 (SVAT) corres- |
C   | ^^^^^^^ pond to bare soil (no vegetation).                        |
C   |         NSTvfr and NSTsfr give vegetation fraction in % (integer) |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE GLOfrc

      IMPLICIT none


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'LOCfil.inc'
      INCLUDE 'NetCDF.inc'
      INCLUDE 'NESTOR.inc'

C +---Local variables
C +   ---------------

      INTEGER nbchar,
     .        i,j,k,l,ii,jj,size_X(nbdom),size_Y(nbdom),
     .        i1,i2,j1,j2,i_cent,j_cent,G_nx,G_ny,totvfr,
     .        ncid(nbdom),start(3),count(3),lmin,frac_ini,
     .        vegtmp,frctmp,ndv1ID(nbdom),ndv2ID(nbdom),
     .        EURcid,EUR1ID,EUR2ID,AFR_size_X,AFR_size_Y,
     .        AFRcid,EUR_size_X,EUR_size_Y,Rcode,frac_max,vauxID,
     .        AFR1ID,AFR2ID,lmax,idom,EUidom,AFidom,error,iauxID,
     .        NAidom,SAidom,NAMcid,NAM1ID,NAM2ID,NAM_size_X,
     .        NAM_size_Y,SAMcid,SAM1ID,SAM2ID,SAM_size_X,
     .        SAM_size_Y,int_1,int_2,ii1,ii2,jj1,jj2,mmx,mmy

      INTEGER int_3,first

      INTEGER*2 val1,val2

      REAL    AUXlo1,AUXla1,AUXlo2,AUXla2,G_reso(nbdom),
     .        AUXlon,AUXlat,G_lon1(nbdom),G_lat1(nbdom),
     .        aux1,aux2,aux3,cmpt1,cmpt2,iAVndv1,iAVndv2,
     .        AFR_G_lon1,AFR_G_lat1,degrad,
     .        AFR_G_reso,AFR_G_lon2,AFR_G_lat2,EUR_G_lat1,
     .        EUR_G_lon1,EUR_G_reso,EUR_G_lon2,EUR_G_lat2,
     .        NAM_G_lon1,NAM_G_lat1,NAM_G_reso,NAM_G_lon2,
     .        NAM_G_lat2,SAM_G_lon1,SAM_G_lat1,SAM_G_reso,
     .        SAM_G_lon2,SAM_G_lat2,Rval1,Rval2,
     .        dx,dy,G_dx,G_dy,VEGfrc,AVndv1,AVndv2,VEGaux

      LOGICAL Vtrue,Vfalse,AFRdom,EURdom,NAMdom,SAMdom,NDVclim

      CHARACTER*2  nustri(0:99)
      CHARACTER*60 EURndvdir,AFRndvdir,NAMndvdir,SAMndvdir
      CHARACTER*80 EURndv_file,AFRndv_file,NAMndv_file,
     .             SAMndv_file


C +---Data
C +   ----

      DATA degrad  / 1.745329252d-2     /
      DATA Vtrue   / .true.             /
      DATA Vfalse  / .false.            /

      DATA (nustri(i),i=0,99)
     .     /'00','01','02','03','04','05','06','07','08','09',
     .      '10','11','12','13','14','15','16','17','18','19',
     .      '20','21','22','23','24','25','26','27','28','29',
     .      '30','31','32','33','34','35','36','37','38','39',
     .      '40','41','42','43','44','45','46','47','48','49',
     .      '50','51','52','53','54','55','56','57','58','59',
     .      '60','61','62','63','64','65','66','67','68','69',
     .      '70','71','72','73','74','75','76','77','78','79',
     .      '80','81','82','83','84','85','86','87','88','89',
     .      '90','91','92','93','94','95','96','97','98','99'/


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Grid parameters
C +   ===============


       EUidom = 1  ! Europe
       AFidom = 2  ! Africa
c #AM  NAidom = 3  ! North America
c #AM  SAidom = 4  ! South America
 

C +---1-KM RESOLUTION DATA
C +   = = = = = = = = = = =


      IF (NDV1km) THEN


C +---Europe
C +   ------

       EUR_G_lat1= 0.1200000E+02
       EUR_G_lon1=-0.1100000E+02
       EUR_G_reso= 0.0100000E+00
       EUR_size_X= 8500
       EUR_size_Y= 6000
       EUR_G_lon2= EUR_G_lon1+REAL(EUR_size_X)*EUR_G_reso
       EUR_G_lat2= EUR_G_lat1+REAL(EUR_size_Y)*EUR_G_reso


C +---Africa
C +   ------

       AFR_G_lat1=-0.3500000E+02
       AFR_G_lon1=-0.2000000E+02
       AFR_G_reso= 0.0100000E+00
       AFR_size_X= 8000
       AFR_size_Y= 7500
       AFR_G_lon2= AFR_G_lon1+REAL(AFR_size_X)*AFR_G_reso
       AFR_G_lat2= AFR_G_lat1+REAL(AFR_size_Y)*AFR_G_reso


C +---North America
C +   -------------

c #AM  NAM_G_lat1= 0.1200000E+02
c #AM  NAM_G_lon1=-0.1100000E+02
c #AM  NAM_G_reso= 0.0100000E+00
c #AM  NAM_size_X= 8500
c #AM  NAM_size_Y= 6000
c #AM  NAM_G_lon2= NAM_G_lon1+REAL(NAM_size_X)*NAM_G_reso
c #AM  NAM_G_lat2= NAM_G_lat1+REAL(NAM_size_Y)*NAM_G_reso


C +---South America
C +   -------------

c #AM  SAM_G_lat1= 0.1200000E+02
c #AM  SAM_G_lon1=-0.1100000E+02
c #AM  SAM_G_reso= 0.0100000E+00
c #AM  SAM_size_X= 8500
c #AM  SAM_size_Y= 6000
c #AM  SAM_G_lon2= SAM_G_lon1+REAL(SAM_size_X)*SAM_G_reso
c #AM  SAM_G_lat2= SAM_G_lat1+REAL(SAM_size_Y)*SAM_G_reso


      ENDIF  ! (NDV1km)


C +---8-KM RESOLUTION DATA
C +   = = = = = = = = = = =


      IF (NDV8km) THEN


C +---Africa
C +   ------

       AFR_G_lat1=-0.380000000000E+02
       AFR_G_lon1=-0.200000000000E+02
       AFR_G_reso= 0.083333333333E+00
       AFR_size_X= 984
       AFR_size_Y= 924
       AFR_G_lon2= AFR_G_lon1+REAL(AFR_size_X)*AFR_G_reso
       AFR_G_lat2= AFR_G_lat1+REAL(AFR_size_Y)*AFR_G_reso


      ENDIF  ! (NDV8km)


C +---Select grid parameters
C +   ----------------------

      AFRdom=.false.
      EURdom=.false.
      NAMdom=.false.
      SAMdom=.false.

      DO j=1,my
      DO i=1,mx

       IF (NSTinc(i,j,AFidom)) AFRdom=.true.
       IF (NSTinc(i,j,EUidom)) EURdom=.true.
c #AM  IF (NSTinc(i,j,NAidom)) NAMdom=.true.
c #AM  IF (NSTinc(i,j,SAidom)) SAMdom=.true.

      ENDDO
      ENDDO


C +---Screen message
C +   ==============

      IF (AFRdom) THEN
       write(6,*) 'NDVI Min/Max index over Africa'
      ENDIF

      IF (EURdom) THEN
       write(6,*) 'NDVI Min/Max index over Europe'
      ENDIF

      IF (NAMdom) THEN
       write(6,*) 'NDVI Min/Max index over N. America'
      ENDIF

      IF (SAMdom) THEN
       write(6,*) 'NDVI Min/Max index over S. America'
      ENDIF

      IF (AFRdom.or.EURdom.or.NAMdom.or.SAMdom) THEN
       write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       write(6,*)
      ENDIF

      IF ((.not.EURdom).and.(.not.AFRdom).and.
     .    (.not.NAMdom).and.(.not.SAMdom)) THEN
       write(6,*) '***************'
       write(6,*) '*** CAUTION ***'
       write(6,*) '***************'
       write(6,*)
       write(6,*) 'No NDVI index available for this domain !!!'
       write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       write(6,*)
       GOTO 990
      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Select MIN/MAX NDVI file
C +   ========================


C +---Year and month ?
C +   ----------------

C +        ******
      CALL DATcnv (RUNiyr,RUNmma,RUNjda,RUNjhu,DATtim,Vfalse)
C +        ******


      IF (NDV1km) THEN


C +....Africa

       nbchar=1
       AFRndvdir=AFRndv_dir

       DO i=1,60
        IF (AFRndvdir(i:i).ne.' ') nbchar=i
       ENDDO

       AFRndv_file=AFRndvdir(1:nbchar)//'AFRndv.nc'

C +... Europe

       nbchar=1
       EURndvdir=EURndv_dir

       DO i=1,60
        IF (EURndvdir(i:i).ne.' ') nbchar=i
       ENDDO

       EURndv_file=EURndvdir(1:nbchar)//'EURndv.nc'


C +... North America

c #AM  nbchar=1
c #AM  NAMndvdir=NAMndv_dir

c #AM  DO i=1,60
c #AM   IF (NAMndvdir(i:i).ne.' ') nbchar=i
c #AM  ENDDO

c #AM  NAMndv_file=NAMndvdir(1:nbchar)//'NAMndv.nc'


C +... South America

c #AM  nbchar=1
c #AM  SAMndvdir=SAMndv_dir

c #AM  DO i=1,60
c #AM   IF (SAMndvdir(i:i).ne.' ') nbchar=i
c #AM  ENDDO

c #AM  SAMndv_file=SAMndvdir(1:nbchar)//'SAMndv.nc'


      ENDIF  ! (NDV1km)


      IF (NDV8km) THEN


C +....Africa

       nbchar=1
       AFRndvdir=AFRndv8dir

       DO i=1,60
        IF (AFRndvdir(i:i).ne.' ') nbchar=i
       ENDDO

       int_1 = RUNiyr/100
       int_2 = RUNiyr - (int_1*100)

       AFRndv_file=AFRndvdir(1:nbchar)//'AFRndv.'
     .                                //nustri(int_1 )
     .                                //nustri(int_2 )
     .                                //'.nc'

      ENDIF  ! (NDV8km)


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Open Netcdf data file - NDVI index
C +   ==================================

      IF (AFRdom) THEN
C +             *****
       AFRcid = NCOPN(AFRndv_file,NCNOWRIT,Rcode)
       AFR1ID = NCVID(AFRcid,'NDVImin',Rcode)
       AFR2ID = NCVID(AFRcid,'NDVImax',Rcode)
C +             *****
      ENDIF

      IF (EURdom.and..not.NDV8km) THEN
C +             *****
       EURcid = NCOPN(EURndv_file,NCNOWRIT,Rcode)
       EUR1ID = NCVID(EURcid,'NDVImin',Rcode)
       EUR2ID = NCVID(EURcid,'NDVImax',Rcode)
C +             *****
      ENDIF

c #AM IF (NAMdom.and..not.NDV8km) THEN
C +             *****
c #AM  NAMcid = NCOPN(NAMndv_file,NCNOWRIT,Rcode)
c #AM  NAM1ID = NCVID(NAMcid,'NDVImin',Rcode)
c #AM  NAM2ID = NCVID(NAMcid,'NDVImax',Rcode)
C +             *****
c #AM ENDIF

c #AM IF (SAMdom.and..not.NDV8km) THEN
C +             *****
c #AM  SAMcid = NCOPN(SAMndv_file,NCNOWRIT,Rcode)
c #AM  SAM1ID = NCVID(SAMcid,'NDVImin',Rcode)
c #AM  SAM2ID = NCVID(SAMcid,'NDVImax',Rcode)
C +             *****
c #AM ENDIF


C +---Initialisation of fraction of vegetation cover
C +   ==============================================

      DO j=1,my
      DO i=1,mx
       NSTfrc(i,j)=0.
       NSTdv1(i,j)=0.
       NSTdv2(i,j)=0.
      ENDDO
      ENDDO


C +---Select domains (Africa and/or Europe)
C +   =====================================

C +---idom = 1 : Europe
C +   -----------------

      G_lat1(EUidom)=EUR_G_lat1
      G_lon1(EUidom)=EUR_G_lon1
      G_reso(EUidom)=EUR_G_reso
      size_X(EUidom)=EUR_size_X
      size_Y(EUidom)=EUR_size_Y
      ndv1ID(EUidom)=EUR1ID
      ndv2ID(EUidom)=EUR2ID
      ncid  (EUidom)=EURcid


C +---idom = 2 : Africa
C +   -----------------

      G_lat1(AFidom)=AFR_G_lat1
      G_lon1(AFidom)=AFR_G_lon1
      G_reso(AFidom)=AFR_G_reso
      size_X(AFidom)=AFR_size_X
      size_Y(AFidom)=AFR_size_Y
      ndv1ID(AFidom)=AFR1ID
      ndv2ID(AFidom)=AFR2ID
      ncid  (AFidom)=AFRcid


C +---idom = 3 : North America
C +   ------------------------

c #AM G_lat1(NAidom)=NAM_G_lat1
c #AM G_lon1(NAidom)=NAM_G_lon1
c #AM G_reso(NAidom)=NAM_G_reso
c #AM size_X(NAidom)=NAM_size_X
c #AM size_Y(NAidom)=NAM_size_Y
c #AM ndv1ID(NAidom)=NAM1ID
c #AM ndv2ID(NAidom)=NAM2ID
c #AM ncid  (NAidom)=NAMcid


C +---idom = 4 : South America
C +   ------------------------

c #AM G_lat1(SAidom)=SAM_G_lat1
c #AM G_lon1(SAidom)=SAM_G_lon1
c #AM G_reso(SAidom)=SAM_G_reso
c #AM size_X(SAidom)=SAM_size_X
c #AM size_Y(SAidom)=SAM_size_Y
c #AM ndv1ID(SAidom)=SAM1ID
c #AM ndv2ID(SAidom)=SAM2ID
c #AM ncid  (SAidom)=SAMcid


C +---Search for MIN/MAX values of NDVI
C +   =================================

      NDVmin = 1000.
      NDVmax = 0.
                             NDVclim=.false.
      IF (Region .eq. "AFW") NDVclim=.true.

      IF (NDVclim) THEN

       int_3  = NCOPN('./input/NDVI08/maxNDVI83-92.nc',NCNOWRIT,Rcode)
       iauxID = NCVID(int_3,'NDVImax',Rcode)
c #BUG error  = nf_open('./input/NDVI08/maxNDVI83-92.nc',iauxID)
c #BUG error  = nf_inq_varid(iauxID,'NDVImax',vauxID)
       NDVmin(idom) = 135.
       NDVmax(idom) = 240.

      ELSE

      DO idom=1,nbdom

c       IF (idom.eq.2.and.AFRdom) THEN

c        DO l=1,size_Y(idom)
c        DO k=1,size_X(idom)
 
c         start(1)=k
c         start(2)=l
c         count(1)=1
c         count(2)=1

C +           *****
c         CALL NCVGT (ncid(idom),ndv1ID(idom),start,count,val1,Rcode)
c         CALL NCVGT (ncid(idom),ndv2ID(idom),start,count,val2,Rcode)
C +           *****

c         Rval1  = REAL(val1)
c         IF (Rval1.gt.2.) NDVmin(idom) = MIN(NDVmin(idom),Rval1)
C +...            ^^^ value (in input file) corresponding to bare soil
C +...                e.g. : desert (Sahara) or city
 
c         Rval2  = REAL(val2)
c         NDVmax(idom) = MAX(NDVmax(idom),Rval2)
C +...            ^^^ value (in input file) corresponding to dense vegetation
C +...                e.g. : tropical forest

c        ENDDO   ! k=1,size_X(idom)
c        ENDDO   ! l=1,size_Y(idom)

c       ENDIF   ! Condition on idom, EURdom and AFRdom

       IF ((idom.eq.2.and.AFRdom)) THEN
        NDVmin(idom) = 49.
        NDVmax(idom) = 185.
       ENDIF

       IF ((idom.eq.1.and.EURdom)) THEN
        NDVmin(idom) = 98.
        NDVmax(idom) = 185.
       ENDIF

      ENDDO    ! idom=1,nbdom
      ENDIF

      NDVmin = NDVmin * 1.10
C +...Corrected NDVmin


C +---Treatment of each NST grid point (except boundaries)
C +   ====================================================

      mmx = mx
      mmy = my
      ii1 = MIN(2,mmx)
      ii2 = MAX(1,mmx-1)
      IF (mmx.eq.1) THEN
       ii1 = 1
       ii2 = 1
      ENDIF
      jj1 = MIN(2,mmy)
      jj2 = MAX(1,mmy-1)
      IF (mmy.eq.1) THEN
       jj1 = 1
       jj2 = 1
      ENDIF


      DO j=jj1,jj2  ! Loop on NST grid points
      DO i=ii1,ii2  ! -----------------------


C +---Initialisation of NDVI variables
C +   ================================

       AVndv1=0.
       AVndv2=0.


C +---Location of NST grid cell in the input data grid
C +   ================================================

       DO idom=1,nbdom

       IF (NSTinc(i,j,idom).and.(G_reso(idom).gt.0.)) THEN


C +---Search for the closest point in data file
C +   -----------------------------------------

        AUXlon = NST__x(i,j)
        AUXlat = NST__y(i,j)
C +          ******
        CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +          ******

        i_cent=NINT((AUXlon-G_lon1(idom))/G_reso(idom))+1
        j_cent=NINT((AUXlat-G_lat1(idom))/G_reso(idom))+1


C +---Check if (i,j) belong to the data domain
C +   ----------------------------------------

        IF (i_cent.lt.3 .or. i_cent.gt.(size_X(idom)-2) .or.
     .      j_cent.lt.3 .or. j_cent.gt.(size_Y(idom)-2)) GOTO 900


C +---Compute the resolution of the considered NST cell
C +   -------------------------------------------------

        dx = 0.0
        dy = 0.0

        IF (mmx.ne.1.and.mmy.ne.1) THEN
         ii=MAX(i,2)
         jj=MAX(j,2)
         AUXlo1 = NST__x(ii  ,jj  )
         AUXlo2 = NST__x(ii-1,jj-1)
         AUXla1 = NST__y(ii  ,jj  )
         AUXla2 = NST__y(ii-1,jj-1)
        ELSE
         IF (mmx.ne.1) THEN
          ii=MAX(i,2)
          jj=1
          AUXlo1 = NST__x(ii  ,jj)
          AUXlo2 = NST__x(ii-1,jj)
          AUXla1 = NST__y(ii  ,jj)
          AUXla2 = NST__y(ii-1,jj)
         ELSE
          ii=1
          jj=1
          AUXlo1 = NST__x(ii,jj)
          AUXlo2 = NST__x(ii,jj)
          AUXla1 = NST__y(ii,jj)
          AUXla2 = NST__y(ii,jj)
         ENDIF
        
        ENDIF
        
C +          ******
        CALL SPHERC (Vtrue,AUXlo1,AUXla1)
        CALL SPHERC (Vtrue,AUXlo2,AUXla2)
C +          ******
        dx=ABS(AUXlo1-AUXlo2)*111111.*COS(AUXla1*degrad)
        dy=ABS(AUXla1-AUXla2)*111111.


C +---Define the data points to be read around (i_cent,j_cent)
C +   --------------------------------------------------------

        G_dx = G_reso(idom)*111111.*COS(AUXla1*degrad)
        G_dy = G_reso(idom)*111111.

        G_nx=NINT(dx/G_dx/2.)-1
        G_ny=NINT(dy/G_dy/2.)-1

        G_nx=MAX(G_nx,0)
        G_ny=MAX(G_ny,0)

        first=0.
1000    continue        

        IF (mmx.eq.1) G_nx=0
        IF (mmy.eq.1) G_ny=0

        i1=i_cent-G_nx
        i2=i_cent+G_nx
        j1=j_cent-G_ny
        j2=j_cent+G_ny

        i1=MAX(i1,1)
        i2=MIN(i2,size_X(idom))
        j1=MAX(j1,1)
        j2=MIN(j2,size_Y(idom))

C +---Initialisation of temporary NDVI variables
C +   ==========================================

        iAVndv1=0.
        iAVndv2=0.
        cmpt1  =0.
        cmpt2  =0.


C +---Reading of input data
C +   =====================

        IF (NDVclim) THEN

         DO l=j1,j2
         DO k=i1,i2

          start(1)=k
          start(2)=l
          count(1)=1
          count(2)=1

C +           *****
          CALL NCVGT (int_3,iauxID,start,count,val2,Rcode)
C +           *****
          IF (NDV8km) THEN
           IF (val2.gt.100.) THEN
            iAVndv2=iAVndv2+REAL(val2)
            cmpt2  =cmpt2  +1.
           ENDIF
          ENDIF
c +
         ENDDO
         ENDDO
c +
         IF (cmpt2.gt.0.) THEN
          iAVndv2 = iAVndv2 / cmpt2
         ELSE
          iAVndv2 = 0.
         ENDIF

        ELSE
c +
         DO l=j1,j2
         DO k=i1,i2

          start(1)=k
          start(2)=l
          count(1)=1
          count(2)=1

C +           *****
          CALL NCVGT (ncid(idom),ndv1ID(idom),start,count,val1,Rcode)
          CALL NCVGT (ncid(idom),ndv2ID(idom),start,count,val2,Rcode)
C +           *****

          IF (NDV1km) THEN
           IF (val1.gt.100) THEN
            iAVndv1=iAVndv1+REAL(val1)
            cmpt1  =cmpt1  +1.
           ENDIF
           IF (val2.gt.100) THEN
            iAVndv2=iAVndv2+REAL(val2)
            cmpt2  =cmpt2  +1.
           ENDIF
          ENDIF

          IF (NDV8km) THEN
           IF (val1.gt.1) THEN
            iAVndv1=iAVndv1+REAL(val1)
            cmpt1  =cmpt1  +1.
           ENDIF
           IF (val2.gt.1) THEN
            iAVndv2=iAVndv2+REAL(val2)
            cmpt2  =cmpt2  +1.
           ENDIF
          ENDIF

         ENDDO
         ENDDO

         IF (cmpt1.gt.0.) THEN
          iAVndv1 = iAVndv1 / cmpt1
         ELSE
          iAVndv1 = 0.
         ENDIF

         IF (cmpt2.gt.0.) THEN
          iAVndv2 = iAVndv2 / cmpt2
         ELSE
          iAVndv2 = 0.
         ENDIF
C +
        ENDIF ! Condition on NDVclim


        IF (iAVndv1.ge.AVndv1 .and. iAVndv2.ge.AVndv2) THEN
         AVndv1=iAVndv1
         AVndv2=iAVndv2
        ENDIF

        IF (NSTsol(i,j).eq.4.and.AVndv1.eq.0.and.
     .                           AVndv2.eq.0.and.first.eq.0) THEN
         G_nx=G_nx*2.
         G_ny=G_ny*2.
         print *,"WARNING: no NDVImin data for (DOM,i,j):"
     .          ,idom,i,j
         first=1
         goto 1000
        ENDIF

       ENDIF  ! Condition on NSTinc

       ENDDO  ! Loop on idom

C +---Compute normalized NDVI index
C +   =============================

       DO idom=1,nbdom
        IF (NSTinc(i,j,idom)) THEN

         NSTdv1(i,j) = AVndv1
         NSTdv1(i,j) = MIN(NSTdv1(i,j),NDVmax(idom))
         NSTdv1(i,j) = MAX(NSTdv1(i,j),NDVmin(idom))

         NSTdv2(i,j) = AVndv2
         NSTdv2(i,j) = MIN(NSTdv2(i,j),NDVmax(idom))
         NSTdv2(i,j) = MAX(NSTdv2(i,j),NDVmin(idom))
      
C +---Estimate of vegetation cover
C +   ============================

         VEGaux      = 
     .   (NSTdv2(i,j) - NDVmin(idom)) / (NDVmax(idom) - NDVmin(idom))
        ENDIF

C +---Exclusion grid cell dominated by water, ice or snow
C +   ===================================================

        IF(NSTsol(i,j).le.3) THEN
         NSTdv1(i,j) = 0
         NSTdv2(i,j) = 0
        ENDIF

       ENDDO

       VEGfrc = VEGaux*(2.-VEGaux-EXP(-2.5*VEGaux))
       VEGfrc = MAX(VEGfrc,0.01)

       IF (NSTsol(i,j).le.3) THEN
        NSTfrc(i,j) = 0.
       ELSE
        NSTfrc(i,j) = VEGfrc
       ENDIF


C +---Modification of fractions of IGBP vegetation cover
C +   ==================================================

       IF (NSTsol(i,j).ge.4) THEN

C +...  Search for the less dominant class
        frac_ini=100
        DO l=1,nvx
         IF (NSTvfr(i,j,l).lt.frac_ini) THEN
          lmin    =l
          frac_ini=NSTvfr(i,j,l)
         ENDIF
        ENDDO

C +...  Attribution of bare soil type (convention = -1)
        NSTveg(i,j,lmin)=-1
        NSTvfr(i,j,lmin)=NINT((1.-VEGfrc)*100.)

C +...  Normalization of NSTvfr
        totvfr=0
        DO l=1,nvx
         totvfr=totvfr+NSTvfr(i,j,l)
        ENDDO
        totvfr=totvfr-NSTvfr(i,j,lmin)
        IF (totvfr.ne.0) THEN
         DO l=1,nvx
          IF (l.ne.lmin) THEN
           aux1         =REAL(NSTvfr(i,j,l))
           aux2         =REAL(totvfr)
           aux3         =REAL(NSTvfr(i,j,lmin))
           NSTvfr(i,j,l)=NINT(aux1/aux2*(100.-aux3))
          ENDIF
         ENDDO
        ELSE
         DO l=1,nvx
          NSTvfr(i,j,l) =0
         ENDDO
         NSTvfr(i,j,lmin)=100
        ENDIF

C +...  Reordering of vegetation types
        lmax    =1
        frac_max=0
        DO l=1,nvx
         IF (NSTvfr(i,j,l).gt.frac_max) THEN
          lmax    =l
          frac_max=NSTvfr(i,j,l)
         ENDIF
        ENDDO
        vegtmp=NSTveg(i,j,lmax)
        frctmp=NSTvfr(i,j,lmax)
        NSTveg(i,j,lmax)=NSTveg(i,j,1)
        NSTvfr(i,j,lmax)=NSTvfr(i,j,1)
        NSTveg(i,j,1)   =vegtmp
        NSTvfr(i,j,1)   =frctmp

       ENDIF


C +---Modification of fractions of SVAT vegetation cover
C +   ==================================================

       IF (NSTsol(i,j).ge.4) THEN

C +...  Attribution of bare soil type
        NSTsvt(i,j,nvx)=0
        NSTsfr(i,j,nvx)=NINT((1.-VEGfrc)*100.)

C +...  Normalization of NSTsfr
        totvfr=0
        DO l=1,nvx
         totvfr=totvfr+NSTsfr(i,j,l)
        ENDDO
        totvfr=totvfr-NSTsfr(i,j,nvx)
        IF (totvfr.ne.0) THEN
         DO l=1,nvx-1
          aux1         =REAL(NSTsfr(i,j,l))
          aux2         =REAL(totvfr)
          aux3         =REAL(NSTsfr(i,j,nvx))
          NSTsfr(i,j,l)=NINT(aux1/aux2*(100.-aux3))
         ENDDO
        ELSE
         DO l=1,nvx
          NSTsfr(i,j,l) =0
         ENDDO
         NSTsfr(i,j,nvx)=100
        ENDIF

C +...  Reordering of vegetation types
ccccc   lmax    =1
ccccc   frac_max=0
ccccc   DO l=1,nvx-1
ccccc    IF (NSTsfr(i,j,l).gt.frac_max) THEN
ccccc     lmax    =l
ccccc     frac_max=NSTsfr(i,j,l)
ccccc    ENDIF
ccccc   ENDDO
ccccc   vegtmp=NSTsvt(i,j,lmax)
ccccc   frctmp=NSTsfr(i,j,lmax)
ccccc   NSTsvt(i,j,lmax)=NSTsvt(i,j,1)
ccccc   NSTsfr(i,j,lmax)=NSTsfr(i,j,1)
ccccc   NSTsvt(i,j,1)   =vegtmp
ccccc   NSTsfr(i,j,1)   =frctmp

       ENDIF


900    CONTINUE


      ENDDO   ! Loop on NST grid points
      ENDDO   ! -----------------------

C +---Close Netcdf data file
C +   ======================

      IF (AFRdom) THEN
C +         ******
       CALL NCCLOS(AFRcid,Rcode)
C +         ******
      ENDIF

      IF (EURdom.and..not.NDV8km) THEN
C +         ******
       CALL NCCLOS(EURcid,Rcode)
C +         ******
      ENDIF

c #AM IF (NAMdom.and..not.NDV8km) THEN
C +         ******
c #AM  CALL NCCLOS(NAMcid,Rcode)
C +         ******
c #AM ENDIF

c #AM IF (SAMdom.and..not.NDV8km) THEN
C +         ******
c #AM  CALL NCCLOS(SAMcid,Rcode)
C +         ******
c #AM ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


990   CONTINUE

      RETURN
      END
