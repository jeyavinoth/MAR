C   +-------------------------------------------------------------------+
C   |  Subroutine GLOveg                            March 2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | GLOveg read vegetation classification (IGBP) only for Africa and  |
C   | Europe.                                                           |
C   |                                                                   |
C   | Input : - NST__x, NST__y : NST grid coordinates (lat./long.)      |
C   | ^^^^^^^ - NST_sh : surface elevation                              |
C   |         - NSTsol : soil type                                      |
C   |                                                                   |
C   | Output: - NSTveg : vegetation type (IGBP classification)          |
C   | ^^^^^^^ - NSTvfr : fraction of vegetation in the grid cell (IGBP) |
C   |         - NSTsvt : vegetation type (SVATclassification)           |
C   |         - NSTsfr : fraction of vegetation in the grid cell (SVAT) |
C   |         - NSTlai : leaf area index                                |
C   |                                                                   |
C   | Remark: Note that NSTveg = -1 (IGBP) or NSTsvt = 0 (SVAT) corres- |
C   | ^^^^^^^ pond to bare soil (no vegetation).                        |
C   |         NSTvfr and NSTsfr give vegetation fraction in % (integer) |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE GLOveg

      IMPLICIT none


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'LOCfil.inc'
      INCLUDE 'NetCDF.inc'
      INCLUDE 'NESTOR.inc'

C +---Vegetation classes
C +   ------------------

      INTEGER nsvat,nigbp
      PARAMETER (nsvat=12)
      PARAMETER (nigbp=17)


C +---Local variables
C +   ---------------

      INTEGER*2 val_IGBP

      INTEGER i,j,k,l,ii,jj,size_X(nbdom),size_Y(nbdom),nbchar,
     .        i1,i2,j1,j2,i_cent,j_cent,G_nx,G_ny,totvfr,Rcode,
     .        AFR_ID,ncid(nbdom),start(3),count(3),AFR_size_X,
     .        VEG_ID(nbdom),AFR_size_Y,EUR_size_X,EUR_size_Y,
     .        EUR_ID,EURcid,AFRcid,idom,frac_itot,EUidom,AFidom,
     .        ii1,ii2,jj1,jj2,mmx,mmy,ido,idomi,s0,s1,
     .        US_size_X,US_size_Y,US_ID,UScid,USidom

      REAL    AUXlo1,AUXla1,AUXlo2,AUXla2,dx,dy,degrad,frac_tot,
     .        AUXlon,AUXlat,G_dx,G_dy,G_lon1(nbdom),G_lat1(nbdom),
     .        aux,aux1,aux2,aux3,AFR_G_lon1,AFR_G_lat1,AFR_G_reso,
     .        cmpt,AFR_G_lon2,AFR_G_lat2,EUR_G_lat1,EUR_G_lon1,
     .        icmpt,EUR_G_reso,EUR_G_lon2,EUR_G_lat2,G_reso(nbdom),
     .        NSTmsk(mx,my),tmp_z0(mx,my),
     .        US_G_reso,US_G_lon1,US_G_lon2,US_G_lat1,US_G_lat2

      INTEGER svat_class(nvx),WK_tmp(nigbp)

      REAL    SVAT(0:nsvat),IGBP(nigbp),convert(nigbp,0:nsvat),
     .        svat_frac (nvx),iIGBP(nigbp),igbp_z0(nigbp)

      LOGICAL Vtrue,AFRdom,EURdom,USdom

      CHARACTER*80 AFRveg_file,EURveg_file,USveg_file


C +---Data
C +   ----

      DATA degrad  / 1.745329252d-2     /
      DATA Vtrue   / .true.             /


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Grid parameters
C +   ===============


C +---Europe
C +   ------

      EUidom = 1

      EUR_G_lat1= 0.1200000E+02
      EUR_G_lon1=-0.1100000E+02
      EUR_G_reso= 0.0100000E+00
      EUR_size_X= 8500         
      EUR_size_Y= 6000         
      EUR_G_lon2= EUR_G_lon1+REAL(EUR_size_X)*EUR_G_reso
      EUR_G_lat2= EUR_G_lat1+REAL(EUR_size_Y)*EUR_G_reso


C +---Africa
C +   ------

      AFidom = 2

      AFR_G_lat1=-0.3500000E+02
      AFR_G_lon1=-0.2000000E+02
      AFR_G_reso= 0.0100000E+00
      AFR_size_X= 8000         
      AFR_size_Y= 7500
      AFR_G_lon2= AFR_G_lon1+REAL(AFR_size_X)*AFR_G_reso
c     AFR_G_lat2= AFR_G_lat1+REAL(AFR_size_Y)*AFR_G_reso
      AFR_G_lat2= 0.3800000E+02

C +-- JJ's adding VEGE data for US
 
      USidom = 3

      US_G_lat1=-0.9999999E+01
      US_G_lon1=-0.1311400E+03
      US_G_reso= 0.0100000E+00
      US_size_X= 5001         
      US_size_Y= 6000         
      US_G_lon2= US_G_lon1+REAL(US_size_X)*US_G_reso
      US_G_lat2= US_G_lat1+REAL(US_size_Y)*US_G_reso
c      US_G_lat2= 0.500000E+02
     
      write(*,*) "US Bounding Vegetation Box" 
      write(*,*) US_G_lat1, US_G_lat2, US_G_lon1, US_G_lon2
      write(*,*) 


C +---Select grid parameters
C +   ----------------------

      AFRdom=.false.
      EURdom=.false.
      USdom =.false.

      DO j=1,my
      DO i=1,mx

       AUXlon=NST__x(i,j)
       AUXlat=NST__y(i,j)

       IF (AUXlon.gt.AFR_G_lon1.and.AUXlon.lt.AFR_G_lon2.and.
     .     AUXlat.gt.AFR_G_lat1.and.AUXlat.lt.AFR_G_lat2) THEN
        AFRdom            =.true.
        NSTinc(i,j,AFidom)=.true.
       ELSE
        NSTinc(i,j,AFidom)=.false.
       ENDIF

       IF (AUXlon.gt.EUR_G_lon1.and.AUXlon.lt.EUR_G_lon2.and.
     .     AUXlat.gt.EUR_G_lat1.and.AUXlat.lt.EUR_G_lat2) THEN
        EURdom            =.true.
        NSTinc(i,j,EUidom)=.true.
       ELSE
        NSTinc(i,j,EUidom)=.false.
       ENDIF

c     added US check for domain criteria - JJ

       IF (AUXlon.gt.US_G_lon1.and.AUXlon.lt.US_G_lon2.and.
     .     AUXlat.gt.US_G_lat1.and.AUXlat.lt.US_G_lat2) THEN
        USdom            =.true.
        NSTinc(i,j,USidom)=.true.
       ELSE
        NSTinc(i,j,USidom)=.false.
       ENDIF
       
      ENDDO
      ENDDO


C +---Screen message
C +   ==============

c      IF (Region.eq."GRD") GOTO 990

      IF (AFRdom) THEN
       write(6,*) 'Global land cover (IGBP) over Africa'
      ENDIF

      IF (EURdom) THEN
       write(6,*) 'Global land cover (IGBP) over Europe'
      ENDIF

      IF (USdom) THEN
       write(6,*) 'Global land cover (IGBP) over United States'
      ENDIF


      IF (AFRdom.or.EURdom.or.USdom) THEN
       write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      ENDIF

      IF ((.not.EURdom).and.(.not.AFRdom).and.(.not.USdom)) THEN
       write(6,*) '***************'
       write(6,*) '*** CAUTION ***'
       write(6,*) '***************'
       write(6,*)
       write(6,*) 'No Global land cover available for this domain !!!'
       write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       write(6,*)
       GOTO 990
      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Convertion table : IGBP  -> SVAT classification
C +   ===============================================


C +---Initialisation
C +   --------------

      DO k=1,nigbp
      DO l=0,nsvat
       convert(k,l)=0.
      ENDDO
      ENDDO


C +---Convertion table
C +   ----------------

C +...1. Evergreen Needleleaf Forest
C +   ------------------------------
      convert( 1,12)=70.   ! needleleaf high
      convert( 1,11)=30.   ! needleleaf medium
      igbp_z0( 1   )=0.94

C +...2. Evergreen Broadleaf Forest
C +   -----------------------------
      convert( 2, 9)=70.   ! broadleaf high
      convert( 2, 8)=30.   ! broadleaf medium
      igbp_z0( 2   )=0.94

C +...3. Deciduous Needleleaf Forest
C +   ------------------------------
      convert( 3,11)=70.   ! needleleaf medium
      convert( 3,12)=30.   ! needleleaf high
      igbp_z0( 3   )=0.86

C +...4. Deciduous Broadleaf Forest
C +   -----------------------------
      convert( 4, 8)=70.   ! broadleaf medium
      convert( 4, 9)=30.   ! broadleaf high
      igbp_z0( 4   )=0.86

C +...5. Mixed Forest
C +   ---------------
      convert( 5, 7)=10.   ! broadleaf  low
      convert( 5, 8)=20.   ! broadleaf  medium
      convert( 5, 9)=20.   ! broadleaf  high
      convert( 5,10)=10.   ! needleleaf low
      convert( 5,11)=20.   ! needleleaf medium
      convert( 5,12)=20.   ! needleleaf high
      igbp_z0( 5   )=0.76

C +...6. Closed Shrublands
C +   --------------------
      convert( 6, 7)=60.   ! broadleaf low
      convert( 6, 8)=40.   ! broadleaf medium
      igbp_z0( 6   )=0.44

C +...7. Open Shrublands
C +   ------------------
      convert( 7, 5)=30.   ! grass medium
      convert( 7, 7)=40.   ! broadleaf low
      convert( 7, 8)=30.   ! broadleaf medium
      igbp_z0( 7   )=0.33

C +...8. Woody Savannas
C +   -----------------
      convert( 8, 5)=30.   ! grass medium
      convert( 8, 8)=35.   ! broadleaf medium
      convert( 8, 9)=35.   ! broadleaf high
      igbp_z0( 8   )=0.64

C +...9. Savannas
C +   -----------
      convert( 9, 6)=60.   ! grass high
      convert( 9, 8)=40.   ! broadleaf medium
      igbp_z0( 9   )=0.38

C +...10. Grasslands
C +   --------------
      convert(10, 4)=70.   ! grass low
      convert(10, 5)=30.   ! grass medium
      igbp_z0(10   )=0.016

C +...11. Permanent Wetlands
C +   ----------------------
      convert(11, 4)=20.   ! grass low
      convert(11, 5)=50.   ! grass medium
      convert(11, 6)=30.   ! grass high
      igbp_z0(11   )=0.047

C +...12. Croplands
C +   -------------
      convert(12, 1)=20.   ! crops low
      convert(12, 2)=30.   ! crops medium
      convert(12, 3)=20.   ! crops high
      convert(12, 0)=30.   ! barren soil
      igbp_z0(12   )=0.047

C +...13. Urban and Built-Up
C +   ----------------------
      convert(13, 9)=100.  ! broadleaf high
      igbp_z0(13   )=1.0

C +...14. Cropland/Natural Vegetation Mosaic
C +   --------------------------------------
      convert(14, 1)=20.   ! crops low
      convert(14, 2)=20.   ! crops medium
      convert(14, 3)=20.   ! crops high
      convert(14, 5)=20.   ! grass medium
      convert(14, 7)=20.   ! broadleaf low
      igbp_z0(14   )=0.074

C +...15. Snow and Ice
C +   ----------------
C     If dominant, NSTsol is set to 2 (ice) or 3 (snow)
C     depending on the height.
      igbp_z0(15   )=0.001

C +...16. Barren or Sparsely Vegetated
C +   --------------------------------
      convert(16, 4)=20.   ! grass low
      convert(16, 7)=5.    ! broadleaf low
      convert(16, 0)=75.   ! barren soil
      igbp_z0(16   )=0.022

C +...17. Water Bodies
C +   ----------------
C     If dominant, NSTsol is set to 1
      igbp_z0(17   )=0.001


C +---Correction
C +   ----------

      if (nvx.le.4) then

       DO k=1,nigbp
  
        DO s0=1,10,3
  
         ! S0=1 [crop]     , S0=2 [grass]
         ! S0=3 [broadleaf], S0=4 [needleleaf]

                                                s1=s0+1

         if(convert(k,s0  ).gt.convert(k,s0+1).and.
     .      convert(k,s0  ).gt.convert(k,s0+2)) s1=s0 

         if(convert(k,s0+2).gt.convert(k,s0).and.
     .      convert(k,s0+2).gt.convert(k,s0+1)) s1=s0+2 

         convert(k,s1)=convert(k,s0)+convert(k,s0+1)+convert(k,s0+2)
 
        ENDDO
  
       ENDDO

      endif

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Open Netcdf data file(s) - IGBP Classification
C +   ==============================================

      IF (AFRdom) THEN
       nbchar=1
       DO i=1,60
        IF (AFRveg_dir(i:i).ne.' ') nbchar=i
       ENDDO
       AFRveg_file = AFRveg_dir(1:nbchar) // 'AFRveg_IGBP.nc'
C +             *****
       AFRcid = NCOPN(AFRveg_file,NCNOWRIT,Rcode)
       AFR_ID = NCVID(AFRcid,'IGBP',Rcode)
C +             *****
      ENDIF

      IF (EURdom) THEN
       nbchar=1
       DO i=1,60
        IF (EURveg_dir(i:i).ne.' ') nbchar=i
       ENDDO
       EURveg_file = EURveg_dir(1:nbchar) // 'EURveg_IGBP.nc'
C +             *****
       EURcid = NCOPN(EURveg_file,NCNOWRIT,Rcode)
       EUR_ID = NCVID(EURcid,'IGBP',Rcode)
C +             *****
      ENDIF

      IF (USdom) THEN
       nbchar=1
       DO i=1,60
        IF (USveg_dir(i:i).ne.' ') nbchar=i
       ENDDO
       USveg_file = USveg_dir(1:nbchar) // 'USveg_IGBP.nc'
C +             *****
       UScid = NCOPN(USveg_file,NCNOWRIT,Rcode)
       US_ID = NCVID(UScid,'IGBP',Rcode)
C +             *****
      ENDIF


C +---Select domains (Africa and/or Europe)
C +   =====================================

C +---idom = 1 : Europe
C +   -----------------

      G_lat1(EUidom)=EUR_G_lat1
      G_lon1(EUidom)=EUR_G_lon1
      G_reso(EUidom)=EUR_G_reso
      size_X(EUidom)=EUR_size_X
      size_Y(EUidom)=EUR_size_Y
      VEG_ID(EUidom)=EUR_ID
      ncid  (EUidom)=EURcid


C +---idom = 2 : Africa
C +   -----------------

      G_lat1(AFidom)=AFR_G_lat1
      G_lon1(AFidom)=AFR_G_lon1
      G_reso(AFidom)=AFR_G_reso
      size_X(AFidom)=AFR_size_X
      size_Y(AFidom)=AFR_size_Y
      VEG_ID(AFidom)=AFR_ID
      ncid  (AFidom)=AFRcid

C +---idom = 3 : United States - JJ 
C +   -----------------

      G_lat1(USidom)=US_G_lat1
      G_lon1(USidom)=US_G_lon1
      G_reso(USidom)=US_G_reso
      size_X(USidom)=US_size_X
      size_Y(USidom)=US_size_Y
      VEG_ID(USidom)=US_ID
      ncid  (USidom)=UScid


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

C +---Initialisation of IGBP variables
C +   ================================

       cmpt=0.
       DO k=1,nigbp
        IGBP(k)=0.
       ENDDO


C +---Location of NST grid cell in the input data grid
C +   ================================================

       DO idom=1,nbdom

       IF (NSTinc(i,j,idom)) THEN


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

        i1=i_cent-G_nx
        i2=i_cent+G_nx
        j1=j_cent-G_ny
        j2=j_cent+G_ny

        i1=MAX(i1,1)
        i2=MIN(i2,size_X(idom))
        j1=MAX(j1,1)
        j2=MIN(j2,size_Y(idom))


C +---Initialisation of temporary IGBP variables
C +   ==========================================

        icmpt=0.
        DO k=1,nigbp
         iIGBP(k)=0.
        ENDDO

        tmp_z0(i,j) = 0.


C +---Reading of input data
C +   =====================

        DO l=j1,j2
        DO k=i1,i2

         start(1)=k
         start(2)=l
         count(1)=1
         count(2)=1
 
C +           *****
         CALL NCVGT (ncid(idom),VEG_ID(idom),start,count,val_IGBP,Rcode)
C +           *****

         IF ((val_IGBP.gt.0).and.(val_IGBP.le.nigbp)) THEN
          icmpt          =icmpt          +1.
          iIGBP(val_IGBP)=iIGBP(val_IGBP)+1.
          tmp_z0(i,j)    =tmp_z0(i,j)    +igbp_z0(val_IGBP)
         ENDIF

        ENDDO
        ENDDO

        aux = 0.
        DO k=1,nigbp-1
         aux = aux + REAL(iIGBP(k))
        ENDDO

        IF (aux.gt.0.) THEN
         cmpt       =icmpt
         NSTmsk(i,j)=idom
         DO k=1,nigbp
          IGBP(k)=iIGBP(k)
         ENDDO
         tmp_z0(i,j) = tmp_z0(i,j)/aux
        ENDIF

       ENDIF  ! NSTinc

       ENDDO  ! Loop on idom


C +---Particular case 1 : water area
C +   ==============================

       IF (IGBP(17).gt.(cmpt/2.).and.NST_sh(i,j).le.300) THEN
        cmpt       =0.
        NSTsol(i,j)=1
       ENDIF

C +---Particular case 2 : dominant ice/snow
C +   =====================================

       IF (IGBP(15).gt.(cmpt/2.)) THEN 
        cmpt       =0. 
        NSTsol(i,j)=3
        write(6,*)
        write(6,*)
     .  'WARNING (GLOveg.f): snow/ice imposed for grid point ',i,j
        write(6,*)
     .  '                    You must initialise SISVAT snow model !!'
       ENDIF

C +---Particular case 3 : dominant land
C +   =================================

       IF (cmpt.GT.0.1E-10.and.NSTsol(i,j).le.2) NSTsol(i,j)=4

C +    **************************
       IF (NSTsol(i,j).ge.4) THEN  ! Continental areas
C +    **************************


C +---Initialisation of surface variables
C +   ===================================

        NSTsvt(i,j,1)  =  6
        NSTsfr(i,j,1)  =100
        NSTveg(i,j,1)  =  9
        NSTvfr(i,j,1)  =100
        DO k=1,nvx
         NSTlai(i,j,k) =  2.0
         NSTglf(i,j,k) =  1.0
        ENDDO

C +---Dominant IGBP classes
C +   =====================

C +...  Initialization
        DO k=1,nigbp
         WK_tmp(k)=IGBP(k)
        ENDDO

C +...  Search for dominant classes
        DO l=1,nvx
         DO k=1,nigbp
          IF (WK_tmp(k).gt.NSTvfr(i,j,l)) THEN
           NSTvfr(i,j,l)=WK_tmp(k)
           NSTveg(i,j,l)=k
           WK_tmp(k)    =0
          ENDIF
         ENDDO
        ENDDO

C +...  Normalization of NSTvfr
        totvfr=0
        DO l=1,nvx
         totvfr=totvfr+NSTvfr(i,j,l)
        ENDDO
        IF (totvfr.ne.0) THEN
         DO l=1,nvx
          aux1         =REAL(NSTvfr(i,j,l))
          aux2         =REAL(totvfr)
          NSTvfr(i,j,l)=NINT(aux1/aux2*100.)
         ENDDO
        ENDIF


C +---Convertion of IGBP to SVAT classification
C +   =========================================

        IF (cmpt.gt.0.) THEN

C +...   initialisation
C +      ~~~~~~~~~~~~~~
         DO l=0,nsvat
          SVAT(l)=0.
         ENDDO

C +...   convertion to SVAT
C +      ~~~~~~~~~~~~~~~~~~
         DO k=1,nigbp
          DO l=0,nsvat
           SVAT(l)=SVAT(l)+convert(k,l)*IGBP(k)/cmpt
          ENDDO
         ENDDO

C        SVAT(l) is the fraction covered by class l
C        for each SVAT class

C +...   retain the (nvx-1) dominant classes
C +      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         DO ido=1,nvx
           idomi=0
           DO l=0,nsvat
             IF (SVAT(l).GT.SVAT(idomi)) THEN
               idomi=l
             ENDIF
           ENDDO
         svat_class(ido)=idomi
         svat_frac (ido)=SVAT(idomi)
         SVAT(idomi)=0.0
         ENDDO

C +...   class (nvx) was reserved for barren soil
C +      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c        svat_class(nvx) = 0
c        svat_frac (nvx) = SVAT(0)

C +...   normalizing the three dominant fractions
C +      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         frac_tot=0.
         DO l=1,nvx
          frac_tot=frac_tot+svat_frac(l)
         ENDDO
         IF (frac_tot.ne.0.) THEN
          DO l=1,nvx
           svat_frac(l)=svat_frac(l)/frac_tot
          ENDDO
         ENDIF

C +...   attribute classes and fractions to NST variables
C +      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         DO k=1,nvx
          NSTsvt(i,j,k)=     svat_class(k)
          NSTsfr(i,j,k)=NINT(svat_frac (k)*100.)
         ENDDO

        ENDIF


C +---Final check of soil fractions
C +   =============================

        frac_itot=0
        DO l=1,nvx
         frac_itot=frac_itot+NSTsfr(i,j,l)
        ENDDO

        IF (frac_itot.le.0) THEN   ! Imposed bare soil
         NSTsvt(i,j,nvx)=  0
         NSTsfr(i,j,nvx)=100
         NSTveg(i,j,nvx)= -1
         NSTvfr(i,j,nvx)=100
         DO k=1,nvx-1
          NSTsvt(i,j,k)=0
          NSTsfr(i,j,k)=0
          NSTveg(i,j,k)=0
          NSTvfr(i,j,k)=0
         ENDDO
         write(6,*) 'Warning : bare soil imposed for grid point ',i,j
     .              ,frac_itot
        ENDIF


C +---Roughness length
C +   ================

        IF (RUGdat) THEN
         NST_z0(i,j) = tmp_z0(i,j)
C +      NST_r0(i,j) = 0.1*NST_z0(i,j)
        ENDIF


C +---Define max leaf area index
C +   ==========================

        DO l=1,nvx

         IF (NSTsvt(i,j,l).eq. 0) NSTlmx(i,j,l) = 0.0
         IF (NSTsvt(i,j,l).eq. 1) NSTlmx(i,j,l) = 0.6
         IF (NSTsvt(i,j,l).eq. 2) NSTlmx(i,j,l) = 0.9
         IF (NSTsvt(i,j,l).eq. 3) NSTlmx(i,j,l) = 1.2
         IF (NSTsvt(i,j,l).eq. 4) NSTlmx(i,j,l) = 0.7
         IF (NSTsvt(i,j,l).eq. 5) NSTlmx(i,j,l) = 1.4
         IF (NSTsvt(i,j,l).eq. 6) NSTlmx(i,j,l) = 2.0
         IF (NSTsvt(i,j,l).eq. 7.or.NSTsvt(i,j,l).eq.10)
     .    NSTlmx(i,j,l) = 3.0
         IF (NSTsvt(i,j,l).eq. 8.or.NSTsvt(i,j,l).eq.11)
     .    NSTlmx(i,j,l) = 4.5
         IF (NSTsvt(i,j,l).eq. 9.or.NSTsvt(i,j,l).eq.12)
     .    NSTlmx(i,j,l) = 6.0

         NSTlai(i,j,l) = NSTlmx(i,j,l)
         NSTglf(i,j,l) = 1.0

        ENDDO


C +    ****
       ELSE   ! Ocean, ice, snow
C +    ****


        NSTsvt(i,j,nvx)=  0
        NSTsfr(i,j,nvx)=100
        NSTveg(i,j,nvx)= -1
        NSTvfr(i,j,nvx)=100
        DO l=1,nvx
         NSTlai(i,j,l) = 0.0
         NSTglf(i,j,l) = 0.0
        ENDDO


C +    *****
       ENDIF  ! Continental areas
C +    *****

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

      IF (EURdom) THEN
C +         ******
       CALL NCCLOS(EURcid,Rcode)
C +         ******
      ENDIF

      IF (USdom) THEN
C +         ******
       CALL NCCLOS(UScid,Rcode)
C +         ******
      ENDIF



C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


990   CONTINUE

      IF (region.eq."GRD") call USRgrd('GLOveg') ! Greenland
      IF (region.eq."EUR") call USReur('GLOveg') ! Iceland

      write(6,*)


      RETURN
      END
