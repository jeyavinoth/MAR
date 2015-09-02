C   +-------------------------------------------------------------------+
C   |  Subroutine GLOglf                      18 March    2009  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | GLOglf read NDVI index over Africa and Europe to determine the    |
C   | green leaf fraction.                                              |
C   |                                                                   |
C   | Input : - NST__x, NST__y : NST grid coordinates (lat./long.)      |
C   | ^^^^^^^ - NSTsol : soil type                                      |
C   |         - NSTdv1 : minimum NDVI index                             |
C   |         - NSTdv2 : maximum NDVI index                             |
C   |                                                                   |
C   | Output: - NSTglf : green leaf fraction (from NDVI index)          |
C   | ^^^^^^^ - NSTlai : leaf area index (from NDVI index)              |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE GLOglf

      IMPLICIT none


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'LOCfil.inc'
      INCLUDE 'NetCDF.inc'

      LOGICAL            INIwri
      common/GLOglf_log/ INIwri


C +---Local variables
C +   ---------------

      INTEGER nbchar,i,j,k,l,ii,jj,size_X(nbdom),size_Y(nbdom),
     .        j1,j2,i_cent,j_cent,G_nx,G_ny,ncid(nbdom),
     .        start(3),count(3),EUR_ID,AFR_ID,AFR_size_X,
     .        AFR_size_Y,EUR_size_X,EUR_size_Y,AFRcid,
     .        NAM_size_X,NAM_size_Y,SAM_size_X,SAM_size_Y,
     .        i1,i2,EURcid,Rcode,ndviID(nbdom),idom,idat,
     .        EUidom,AFidom,NAidom,SAidom,NAM_ID,SAM_ID,
     .        SAMcid,int_1,int_2,RUNdec,ndat,GLOiyr,GLOmma,
     .        GLOjda,GLOjhu,GLOdec,ii1,ii2,jj1,jj2,mmx,mmy

      INTEGER*2 Ival_NDVI,first

      INTEGER*4 GLOtim

      REAL    AUXlo1,AUXla1,AUXlo2,AUXla2,dx,dy,degrad,
     .        G_dx,G_dy,G_lon1(nbdom),G_lat1(nbdom),G_lon2,
     .        G_lat2,aux1,aux2,aux3,cmpt,iAVndvi,zero,unun,
     .        NDVrea,AFR_G_lon1,AFR_G_lat1,
     .        AFR_G_reso,AFR_G_lon2,AFR_G_lat2,EUR_G_lat1,
     .        EUR_G_lon1,EUR_G_reso,EUR_G_lon2,EUR_G_lat2,
     .        NAM_G_lat1,NAM_G_lon1,NAM_G_reso,NAM_G_lon2,
     .        NAM_G_lat2,SAM_G_lat1,SAM_G_lon1,SAM_G_reso,
     .        SAM_G_lon2,SAM_G_lat2,LAImax,GLFmax,kval,
     .        AUXlon,AUXlat,G_reso(nbdom),AVndvi,aux,alpha,
     .        weight(5),shiftDAY,tmp_z0(mx,my),val_NDVI,
     .        laiaux,ndvaux,alpha1,beta,alpha2,raux

      LOGICAL Vtrue,Vfalse,AFRdom,EURdom,NAMdom,SAMdom

      CHARACTER*2  nustri(0:99)
      CHARACTER*60 AFRndvdir,EURndvdir,NAMndvdir,SAMndvdir
      CHARACTER*80 AFRndv_file,EURndv_file,NAMndv_file,
     .             SAMndv_file


C +---Data
C +   ----

      DATA degrad  / 1.745329252d-2 /
      DATA zero    / 0.             /
      DATA unun    / 1.             /
      DATA Vtrue   / .true.         /
      DATA Vfalse  / .false.        /

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

       AFR_G_lat1=-0.088894023001E+02
       AFR_G_lon1=-0.249166007340E+02
       AFR_G_reso= 0.083333333333E+00
       AFR_size_X= 612
       AFR_size_Y= 587
       AFR_G_lon2= AFR_G_lon1+REAL(AFR_size_X)*AFR_G_reso
       AFR_G_lat2= AFR_G_lat1+REAL(AFR_size_Y)*AFR_G_reso

C +---Europe
C +   ------

       EUR_G_lat1=-0.380000000000E+02
       EUR_G_lon1=-0.200000000000E+02
       EUR_G_reso= 0.083333333333E+00
       EUR_size_X= 984
       EUR_size_Y= 924
       EUR_G_lon2= EUR_G_lon1+REAL(EUR_size_X)*EUR_G_reso
       EUR_G_lat2= EUR_G_lat1+REAL(EUR_size_Y)*EUR_G_reso


      ENDIF  ! (NDV8km)


C +---Select grid parameters
C +   ----------------------

      AFRdom=.false.
      EURdom=.false.
      NAMdom=.false.
      SAMdom=.false.

C     Redefine the domain in order to get:
C     - African  continent: use of the new set of NDVI (only Africa)
C     - European continent: use of the old set of NDVI

      DO j=1,my
      DO i=1,mx
       AUXlon=NST__x(i,j)
       AUXlat=NST__y(i,j)
       IF (AUXlat.le.36.0) THEN
         NSTinc(i,j,EUidom)=.false. !but still missing Crete
         NSTinc(i,j,AFidom)=.true.  !African continent
                                    !but still missing northern Algeria
       ELSE
         IF (AUXlon.le.-15.0) THEN
           NSTinc(i,j,EUidom)=.false. !still in Africa as Europe dataset
           NSTinc(i,j,AFidom)=.true.  !doesn't extend further than 20W
         ELSE
           NSTinc(i,j,EUidom)=.true.  !    in Europe
           NSTinc(i,j,AFidom)=.false. !not in Africa
         ENDIF
       ENDIF
       IF (AUXlat.le.36.0 .and.
     .     AUXlon.ge. 0.0 .and. AUXlon.le.12.0) THEN  !Northern Algeria
         NSTinc(i,j,EUidom)=.false. !not in Europe
         NSTinc(i,j,AFidom)=.true.  !    in Africa
       ENDIF
       IF (AUXlat.ge.34.0 .and.
     .     AUXlon.ge.22.0 .and. AUXlon.le.27.0) THEN  !Crete
         NSTinc(i,j,EUidom)=.true.  !    in Europe
         NSTinc(i,j,AFidom)=.false. !not in Africa
       ENDIF
      ENDDO
      ENDDO


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
       write(6,*) 'Green leaf fraction from NDVI index over Africa'
      ENDIF

      IF (EURdom) THEN
       write(6,*) 'Green leaf fraction from NDVI index over Europe'
      ENDIF

      IF (NAMdom) THEN
       write(6,*) 'Green leaf fraction from NDVI index over N. America'
      ENDIF

      IF (SAMdom) THEN
       write(6,*) 'Green leaf fraction from NDVI index over S. America'
      ENDIF

      IF (AFRdom.or.EURdom.or.NAMdom.or.SAMdom) THEN
       write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
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

C +---Initialisation
C +   ==============

      DO j=1,my
      DO i=1,mx
       NSTndv(i,j) = 0.
      ENDDO
      ENDDO

      DO k=1,nvx
      DO j=1,my
      DO i=1,mx
       NSTglf(i,j,k) = 0.
      ENDDO
      ENDDO
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      IF (NDV8km) THEN
       ndat = 5
      ELSE
       ndat = 1
      ENDIF


      DO idat = 1,ndat
       

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Select a monthly file
C +   =====================


C +---Year, month and decade ?
C +   ------------------------

C +         ******
       CALL DATcnv (RUNiyr,RUNmma,RUNjda,RUNjhu,DATtim,Vfalse)
C +         ******

       IF (RUNjda.lt.11) THEN
        RUNdec = 1
       ELSE
        IF (RUNjda.lt.21) THEN
         RUNdec = 11
        ELSE
         RUNdec = 21
        ENDIF
       ENDIF


       IF (NDV1km) THEN


C +---File name
C +   ---------

C +.... Africa

        nbchar=1
        AFRndvdir=AFRndv_dir

        DO i=1,60
         IF (AFRndvdir(i:i).ne.' ') nbchar=i
        ENDDO

        AFRndv_file=AFRndvdir(1:nbchar)//'AFRn'
     .                                 //nustri(RUNmma)
     .                                 //'.nc'

C +.... Europe

        nbchar=1
        EURndvdir=EURndv_dir

        DO i=1,60
         IF (EURndvdir(i:i).ne.' ') nbchar=i
        ENDDO

        EURndv_file=EURndvdir(1:nbchar)//'EURn'
     .                                 //nustri(RUNmma)
     .                                 //'.nc'

C +.... North America

c #AM   nbchar=1
c #AM   NAMndvdir=NAMndv_dir

c #AM   DO i=1,60
c #AM    IF (NAMndvdir(i:i).ne.' ') nbchar=i
c #AM   ENDDO

c #AM   NAMndv_file=NAMndvdir(1:nbchar)//'NAMn'
c #AM.                                 //nustri(RUNmma)
c #AM.                                 //'.nc'

C +.... South America

c #AM   nbchar=1
c #AM   SAMndvdir=SAMndv_dir

c #AM   DO i=1,60
c #AM    IF (SAMndvdir(i:i).ne.' ') nbchar=i
c #AM   ENDDO

c #AM   SAMndv_file=SAMndvdir(1:nbchar)//'SAMn'
c #AM.                                 //nustri(RUNmma)
c #AM.                                 //'.nc'


       ENDIF  ! (NDV1km)


       IF (NDV8km) THEN

        shiftDAY = (REAL(RUNjda) - REAL(RUNdec) - 5.0)   ! day
     .           +  REAL(RUNjhu)/24.                     ! hour

        IF (RUNiyr.ge.1983) THEN

         IF (idat.le.2) THEN
          IF (RUNmma.ge.2) THEN
           GLOtim = DATtim + NINT(24.*(-shiftDAY - 10.0*REAL(3-idat)))
          ELSE
           GLOtim = DATtim
          ENDIF
         ENDIF
         IF (idat.ge.3) THEN
          GLOtim = DATtim + NINT(24.*(-shiftDAY + 10.0*REAL(idat-3)))
         ENDIF

        ELSE

         write(6,*) 'No NDVI files available before 1983.'
         write(6,*) 'Please select NDVI database at 1-km resolution'
         write(6,*) 'in NSTing.ctr input file.'
         write(6,*)
         write(6,*) 'STOP in GLOglf.f'
         write(6,*)
         STOP

        ENDIF

C +          ******
        CALL DATcnv (GLOiyr,GLOmma,GLOjda,GLOjhu,GLOtim,Vfalse)
C +          ******

        IF (GLOjda.lt.11) THEN
         GLOdec = 1
        ELSE
         IF (GLOjda.lt.21) THEN
          GLOdec = 11
         ELSE
         GLOdec = 21
         ENDIF
        ENDIF


C +---Compute weights for time interpolation of green leaf fraction
C +   -------------------------------------------------------------

        IF (idat.eq.1) THEN
         weight(3) = 10. / 30.
         IF (shiftDAY.le.0.) THEN
          weight(5) = 0.0
          weight(2) = 10. / 30.
          weight(1) = ABS(shiftDAY) / 30.
          weight(4) = 10. / 30. - weight(1)
         ELSE
          weight(1) = 0.0
          weight(4) = 10. / 30.
          weight(5) = ABS(shiftDAY) / 30.
          weight(2) = 10. / 30. - weight(5)
         ENDIF
        ENDIF

C +      ----X---------X---------X---------X---------X---------X
C +          1         2         3         4         5
C +                                  ^^ current date
C +
C +     - X represents the beginning of each decade
C +     - The current day is assumed to be in decade 3 (between 3 and 4)
C +     - shiftDAY is the difference between the current date and the center
C +       of the decade


C +---File name
C +   ---------

C +.... Africa

        nbchar=1
        AFRndvdir=AFRndv8dir

        DO i=1,60
         IF (AFRndvdir(i:i).ne.' ') nbchar=i
        ENDDO

        int_1 = GLOiyr/100
        int_2 = GLOiyr - (int_1*100)

C +     Monthly data files
ccccc   AFRndv_file=AFRndvdir(1:nbchar)//'AFRn.'
ccccc.                                 //nustri(int_1 )
C +     10-days data files
        AFRndv_file=AFRndvdir(1:nbchar)//'ndvi_mar'
     .                                 //nustri(int_2 )
     .                                 //nustri(GLOmma)
     .                                 //nustri(GLOdec)
     .                                 //'.nc'

C +.... Europe

        nbchar=1
        EURndvdir=EURndv8dir

        DO i=1,60
         IF (EURndvdir(i:i).ne.' ') nbchar=i
        ENDDO

        EURndv_file=EURndvdir(1:nbchar)//'avhrrpf.ndvi.1ntfgl.'
     .                                 //nustri(int_2)
     .                                 //nustri(GLOmma)
     .                                 //nustri(GLOdec)
     .                                 //'.nc'

       ENDIF  ! (NDV8km)


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Open Netcdf data file - NDVI index
C +   ==================================

       IF (AFRdom) THEN
        WRITE(6,*) 'Open file : ',AFRndv_file
C +              *****
        AFRcid = NCOPN(AFRndv_file,NCNOWRIT,Rcode)
        AFR_ID = NCVID(AFRcid,'NDVI',Rcode)
C +              *****
       ENDIF

c       IF (EURdom.and..not.NDV8km) THEN
       IF (EURdom) THEN
        WRITE(6,*) 'Open file : ',EURndv_file
C +              *****
        EURcid = NCOPN(EURndv_file,NCNOWRIT,Rcode)
        EUR_ID = NCVID(EURcid,'NDVI',Rcode)
C +              *****
       ENDIF

c #AM  IF (NAMdom.and..not.NDV8km) THEN
c #AM   WRITE(6,*) 'Open file : ',NAMndv_file
C +              *****
c #AM   NAMcid = NCOPN(NAMndv_file,NCNOWRIT,Rcode)
c #AM   NAM_ID = NCVID(NAMcid,'NDVI',Rcode)
C +              *****
c #AM  ENDIF

c #AM  IF (SAMdom.and..not.NDV8km) THEN
c #AM   WRITE(6,*) 'Open file : ',SAMndv_file
C +              *****
c #AM   SAMcid = NCOPN(SAMndv_file,NCNOWRIT,Rcode)
c #AM   SAM_ID = NCVID(SAMcid,'NDVI',Rcode)
C +              *****
c #AM  ENDIF


C +---Select domains (Africa and/or Europe)
C +   =====================================

C +---idom = 1 : Europe
C +   -----------------

       G_lat1(EUidom)=EUR_G_lat1
       G_lon1(EUidom)=EUR_G_lon1
       G_reso(EUidom)=EUR_G_reso
       size_X(EUidom)=EUR_size_X
       size_Y(EUidom)=EUR_size_Y
       ndviID(EUidom)=EUR_ID
       ncid  (EUidom)=EURcid


C +---idom = 2 : Africa
C +   -----------------

       G_lat1(AFidom)=AFR_G_lat1
       G_lon1(AFidom)=AFR_G_lon1
       G_reso(AFidom)=AFR_G_reso
       size_X(AFidom)=AFR_size_X
       size_Y(AFidom)=AFR_size_Y
       ndviID(AFidom)=AFR_ID
       ncid  (AFidom)=AFRcid


C +---idom = 3 : North America
C +   ------------------------

c #AM  G_lat1(NAidom)=NAM_G_lat1
c #AM  G_lon1(NAidom)=NAM_G_lon1
c #AM  G_reso(NAidom)=NAM_G_reso
c #AM  size_X(NAidom)=NAM_size_X
c #AM  size_Y(NAidom)=NAM_size_Y
c #AM  ndviID(NAidom)=NAM_ID
c #AM  ncid  (NAidom)=NAMcid


C +---idom = 4 : South America
C +   ------------------------

c #AM  G_lat1(SAidom)=SAM_G_lat1
c #AM  G_lon1(SAidom)=SAM_G_lon1
c #AM  G_reso(SAidom)=SAM_G_reso
c #AM  size_X(SAidom)=SAM_size_X
c #AM  size_Y(SAidom)=SAM_size_Y
c #AM  ndviID(SAidom)=SAM_ID
c #AM  ncid  (SAidom)=SAMcid


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


C +---Initialisation of temporary NDVI variables
C +   ==========================================

        AVndvi=0.


C +---Location of NST grid cell in the input data grid
C +   ================================================

        DO idom=1,nbdom

        IF (NSTinc(i,j,idom).and.(G_reso(idom).gt.0.)) THEN


C +---Search for the closest point in data file
C +   -----------------------------------------

         AUXlon = NST__x(i,j)
         AUXlat = NST__y(i,j)
C +           ******
         CALL SPHERC (Vtrue,AUXlon,AUXlat)
C +           ******
 
         i_cent=NINT((AUXlon-G_lon1(idom))/G_reso(idom))+1
         j_cent=NINT((AUXlat-G_lat1(idom))/G_reso(idom))+1


C +---Check if (i,j) belong to the data domain
C +   ----------------------------------------

         IF (.NOT.INIwri)                                           THEN
         IF (i_cent.lt.3 .or. i_cent.gt.(size_X(idom)-2) .or.
     .       j_cent.lt.3 .or. j_cent.gt.(size_Y(idom)-2)) THEN
c           write(6,899) i,j,idom
c 899       format('GLOglf.f: (',i3,',',i3,') does not belong to the'
c     .            ' data domain (idom=',i1,')')
           GOTO 900
         ENDIF
         ENDIF


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
        
C +           ******
         CALL SPHERC (Vtrue,AUXlo1,AUXla1)
         CALL SPHERC (Vtrue,AUXlo2,AUXla2)
C +           ******
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

         first=0
1000     continue

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

         iAVndvi=0.
         cmpt   =0.


C +---Reading of input data
C +   =====================

         DO l=j1,j2
         DO k=i1,i2

          start(1)=k
          start(2)=l
          start(3)=1
          count(1)=1
          count(2)=1
          count(3)=1

C +            *****
c        IF (idom.eq.2 .and. NSTinc(i,j,idom)) THEN !read as REAL (NDVI anom.)
c        CALL NCVGT(ncid(idom),ndviID(idom),start,count,val_NDVI,Rcode)
c        ENDIF ! BUG BUG !!!!
c        IF (idom.eq.1 .and. NSTinc(i,j,idom)) THEN !read as INTEGER*2
         CALL NCVGT(ncid(idom),ndviID(idom),start,count,Ival_NDVI,Rcode)
         val_NDVI = REAL(Ival_NDVI)
c        ENDIF
C +            *****

          IF (NDV1km) THEN
           IF (val_NDVI.gt.100.) THEN
            iAVndvi = iAVndvi + val_NDVI
            cmpt    = cmpt + 1.
           ENDIF
          ENDIF
 
          IF (NDV8km) THEN
           IF (val_NDVI.gt.20.) THEN
            iAVndvi = iAVndvi + val_NDVI
            cmpt    = cmpt + 1.
           ENDIF
          ENDIF
 
         ENDDO
         ENDDO

         IF (cmpt.gt.0.) THEN
          iAVndvi = iAVndvi / cmpt
         ELSE
          iAVndvi = 0.
         ENDIF

         IF (iAVndvi.gt.AVndvi) AVndvi=iAVndvi

         IF (NSTsol(i,j).eq.4.and.AVndvi.eq.0.and.first.eq.0) THEN
          G_nx=G_nx*2.
          G_ny=G_ny*2.
          first=1 
          goto 1000
         ENDIF

        ENDIF  ! Condition on NSTinc

        ENDDO  ! Loop on idom


C +---Compute normalized NDVI index
C +   =============================

       DO idom=1,nbdom
        IF (NSTinc(i,j,idom)) THEN

         NDVrea = AVndvi
         NDVrea = MIN(NDVrea,NDVmax(idom))
         NDVrea = MAX(NDVrea,NDVmin(idom))
      
        ENDIF
       ENDDO

C +---Exclusion grid cell dominated by water, ice or snow
C +   ===================================================

        IF (NSTsol(i,j).le.3) THEN
         NDVrea = 0.
        ENDIF 

C +---Time interpolation of NDVI
C +   ==========================

        IF (ndat.gt.1) THEN
         NSTndv(i,j) = NSTndv(i,j) + weight(idat)*NDVrea
        ELSE
         NSTndv(i,j) = NDVrea
        ENDIF


900     CONTINUE


       ENDDO   ! Loop on NST grid points
       ENDDO   ! -----------------------


C +---Close Netcdf data file
C +   ======================

       IF (AFRdom) THEN
C +          ******
        CALL NCCLOS(AFRcid,Rcode)
C +          ******
       ENDIF
 
c       IF (EURdom.and..not.NDV8km) THEN
       IF (EURdom) THEN
C +          ******
        CALL NCCLOS(EURcid,Rcode)
C +          ******
       ENDIF

c #AM  IF (NAMdom.and..not.NDV8km) THEN
C +          ******
c #AM   CALL NCCLOS(NAMcid,Rcode)
C +          ******
c #AM  ENDIF

c #AM  IF (SAMdom.and..not.NDV8km) THEN
C +          ******
c #AM   CALL NCCLOS(SAMcid,Rcode)
C +          ******
c #AM  ENDIF


      ENDDO    ! {idat=1,ndat}
       


      mmx = mx
      mmy = my
      ii1 = MIN(2,mmx)
      ii2 = MAX(1,mmx-1)
      jj1 = MIN(2,mmy)
      jj2 = MAX(1,mmy-1)

 
      DO j=jj1,jj2  ! Loop on NST grid points
      DO i=ii1,ii2  ! -----------------------


C +---Estimate of green leaf fraction
C +   ===============================

       IF (NSTdv1(i,j).ge.0.   .and. 
     .     NSTdv1(i,j).le.256. .and.
     .     NSTdv2(i,j).ge.0.   .and.    
     .     NSTdv2(i,j).le.256. .and.
     .     NSTdv1(i,j).ne.NSTdv2(i,j)) THEN
 
        DO k=1,nvx
         DO idom=1,nbdom
          IF (NSTinc(i,j,idom)) THEN
          NSTglf(i,j,k)= (NSTndv(i,j)-NDVmin(idom))
     .                 / (NSTdv2(i,j)-NDVmin(idom))
          ENDIF
         ENDDO
         NSTglf(i,j,k) = MIN(unun,NSTglf(i,j,k))
         NSTglf(i,j,k) = MAX(zero,NSTglf(i,j,k))

        ENDDO

        IF (region.eq."AFW" .OR. region.eq."WAF") THEN
         DO k=1,nvx
          NSTglf(i,j,k) = 1.0
         ENDDO
        ENDIF

       ELSE

        DO k=1,nvx
         NSTglf(i,j,k) = 1.0
        ENDDO
  
       ENDIF

       IF (NSTsol(i,j).le.3 ) THEN
       
        DO k=1,nvx
         NSTglf(i,j,k) = 0.0
         NSTlai(i,j,k) = 0.0
        ENDDO

       ENDIF
 
C +---Compute leaf area index
C +   =======================

       IF (NSTsol(i,j).ge.4) THEN
       
c       DO l=1,nvx
c        LAImax = NSTlmx(i,j,l)
c        GLFmax = 1.0
c        kval   = 0.5
c        alpha  = (1.0-EXP(-kval*LAImax)) / GLFmax
c        aux    = MIN(0.999,alpha*NSTglf(i,j,l))
c        NSTlai(i,j,l) = -LOG(1.0-aux) / kval
c        NSTlai(i,j,l) = MIN(NSTlai(i,j,l),LAImax)
c       ENDDO
        DO idom=1,nbdom
         IF (NSTinc(i,j,idom)) THEN
         ndvaux = (NSTndv(i,j) -NDVmin(idom))
     .          / (NDVmax(idom)-NDVmin(idom))
         ENDIF
        ENDDO
        ndvaux = MAX(ndvaux,0.)
        ndvaux = MIN(ndvaux,0.99)
        kval = 0.15
        laiaux = -LOG(1.0-ndvaux) / kval
        IF  (NSTsfr(i,j,2).eq.0 ) THEN
             NSTlai(i,j,2) =  0.
         IF (NSTsfr(i,j,1).eq.0)  THEN
             NSTlai(i,j,1) =  0.
         ELSE
             alpha1        = REAL(NSTsfr(i,j,1))/100.
             NSTlai(i,j,1) =      laiaux/alpha1
         ENDIF
        ELSE
         IF (NSTsfr(i,j,1).eq.0)  THEN
             NSTlai(i,j,1) =  0.
             alpha2        = REAL(NSTsfr(i,j,2))/100.
             NSTlai(i,j,2) =      laiaux/alpha2
         ELSE

C +vv Hubert Gall
c #@X        beta          =      NSTlmx(i,j,1) /NSTlmx(i,j,2)  ! BOUM !
c #@X        alpha1        = REAL(NSTsfr(i,j,1))/100.
c #@X        alpha2        = REAL(NSTsfr(i,j,2))/100.
c #@X        NSTlai(i,j,1) = beta*laiaux/(alpha2+beta*alpha1)
c #@X        NSTlai(i,j,2) =      laiaux/(alpha2+beta*alpha1)
c #@X        write(6,6000) i,j   ,NSTlmx(i,j,1),NSTlmx(i,j,2),beta
 6000        format(2i3,3e15.3)
             alpha1        = NSTlmx(i,j,1)*REAL(NSTsfr(i,j,1))/100.
             alpha2        = NSTlmx(i,j,2)*REAL(NSTsfr(i,j,2))/100.
             beta          = laiaux       /    (alpha1 + alpha2)
             NSTlai(i,j,1) = NSTlmx(i,j,1)*beta
             NSTlai(i,j,2) = NSTlmx(i,j,2)*beta
C +^^ Hubert Gall

         ENDIF
        ENDIF
             raux     = 0.01*REAL(NSTsfr(i,j,1) +NSTsfr(i,j,2))
     .       / NSTfrc(i,j)
        DO l = 1,nvx-1
             NSTlai(i,j,l) =      NSTlai(i,j,l) *raux 
        ENDDO
             NSTlai(i,j,nvx) =    0.0   ! Emilie Vanvyve, e-mail 22-8-2006
        
       ENDIF

      ENDDO   ! Loop on NST grid points
      ENDDO   ! -----------------------




C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


990   CONTINUE

      INIwri=.TRUE.

      RETURN
      END
