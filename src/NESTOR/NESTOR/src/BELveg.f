C   +-------------------------------------------------------------------+
C   |  Subroutine BELveg                             August 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : - NST__x, NST__y : horizontal grid of NST model           |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: - NSTveg : vegetation type (IGBP classification)          |
C   | ^^^^^^^ - NSTvfr : fraction of vegetation in the grid cell (IGBP) |
C   |                                                                   |
C   | Source : L. VAN DER AUWERA (IRM - Pollution Dep.)                 |
C   | ^^^^^^                                                            |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Explanation given with the data (L. Van der Auwera) :             |
C   | ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~               |
C   | VOLGENDE FILE BEVAT 'LAND USE' GEGEVENS VAN BELGIE                |  
C   | DE X,Y COORDINATEN ZIJN UITGEDRUKT IN DEZE LAMBERT-COORDINATEN    |
C   | DE INDELING IS VOOR ELKE LIJN:                                    |
C   | VOLGNR,MATRIXPLAATS,STAFKAARTNR,ONDERVERDELINGSNR V. D. STAFKAART,|
C   |   X-LAMBERT COORD,Y-LAMBERT COORD,RUWHEIDSCODE,LAMBDA,PHI,        |
C   |   CODE(HOOGTE INFORMATIEBRON),COMMENTAAR                          |
C   | MATRIXPLAATS=(Y-1)*300+X  ,PUNT (1,1) IS LINKSONDER,RIJ PER RIJ   |
C   | LAMBDA,PHI IN GRADEN EN DECIMALE GRAADINDELING                    |
C   | CODE = 0 HOOGTE EN RUWHEIDSCODE GENOMEN VANAF BELGISCHE STAFKAART |
C   |          1/25000, ALLE PUNTEN BINNEN BELGIE GELEGEN               |
C   |        1 GEINTERPOLEERDE WAARDEN UIT NEDERLANDSE GEGEVENS         |
C   |                                                      (J.WIERINGA) |
C   |        2 GEINTERPOLEERDE WAARDEN UIT ALPEX GEGEVENS               |
C   |        3 GEINTERPOLEERDE WAARDEN UIT DUITSE GEGEVENS (H.SCHMIDT   |
C   |                                                          D.WETT.) |
C   |        4 GEINTERPOLEERDE WAARDEN UIT DUITSE + ALPEX COMBINATIES   |
C   |        5 WAARDEN OP DE NOORDZEE: hoogte = 0 en ruwheidscode = Z   |
C   |  VERKLARING AANGENOMEN RUWHEIDSCODE = (CF. J.WIERINGA)            |
C   |  CODE - OMSCHRIJVING               (MOGELIJKE)RUWHEIDSLENGTE(CM)  |
C   |     Z - ZEE                                     0.03              |
C   |     M - MEER,WATER                              0.6               |
C   |     R - RIET,MOERAS,DRAS                        1.5               |
C   |     P - POLDER,ZAND,WEIDE,GRAS                  7.0               |
C   |     D - DUIN,HEIDE,LAAG KREUPELHOUT            10.0               |
C   |     A - AKKERLAND,BEBOUWD LAND                 17.0               |
C   |     W - WEGEN,KANALEN,BOMEN RIJEN              24.0               |
C   |     G - BOOMGAARDEN,BOSJES                     35.0               |
C   |     B - BOS,LOOFWOUD                           75.0               |
C   |     N - NAALDWOUD                             100.0               |
C   |     H - HUIZEN,DORPEN                         112.0               |
C   |     S - STEDEN,HOOGBOUW                       160.0               |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE BELveg


      IMPLICIT NONE


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'LOCfil.inc'

C +---Local variables
C +   ---------------

      INTEGER i,j,k,l,dim_x,dim_y,px,py,nbx,nby,count,nigbp,
     .        aux1,aux2,aux3,aux4,aux5,TMPuse

      REAL    tmp_lon,tmp_lat,dxi,dyi,dist,dist_max,dist_best,
     .        distx,disty,degrad,VEGtot,TMPfrc

      CHARACTER*1 soil_id

      PARAMETER (dim_x=310,dim_y=260)
      PARAMETER (nigbp=17)

      INTEGER IBGPcl(dim_x,dim_y),VEGcls(mx,my,nvx)
      REAL    VEGlon(dim_x,dim_y),VEGlat(dim_x,dim_y),
     .        VEGfrc(mx,my,3),IGBPfr(nigbp)

      CHARACTER*80 BELveg_file


C +---Data
C +   ----

      DATA degrad  / 1.745329252d-2     /
 

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Initialization
C +   ==============

      DO j=1,my
      DO i=1,mx

C +    Standard vegetation : 60 % natural prairies
C +    ^^^^^^^^^^^^^^^^^^^   20 % agricultural crop
C +                          20 % evergreen forest

       DO k=1,nvx
        VEGcls(i,j,k)=0
        VEGfrc(i,j,k)=0.
       ENDDO

      ENDDO
      ENDDO

     
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     
C +---Read the data file : LANDUSE_Belg
C +   =================================


C +---Initialisation
C +   --------------

      DO j=1,dim_y
      DO i=1,dim_x
       IBGPcl(i,j)=10
      ENDDO
      ENDDO


C +---Data File
C +   ---------

      BELveg_file = BELveg_dir // 'BELveg_IRM.asc'
      OPEN (unit=10,status='old',file=BELveg_file)
      REWIND 10

      DO k=1,37836

       READ (10,100) aux1,aux2,aux3,aux4,i,j,soil_id,
     .               tmp_lon,tmp_lat,aux5

100    FORMAT (2i6,2i3,2i4,5x,a1,f8.5,f9.5,i2)


C +---Longitude and latitude
C +   ----------------------

       VEGlon(i,j) = tmp_lon
       VEGlat(i,j) = tmp_lat


C +---Convertion of land use to IGBP classification
C +   ---------------------------------------------

       IF (soil_id .eq. 'Z')  IBGPcl(i,j) =  17  ! Water
       IF (soil_id .eq. 'M')  IBGPcl(i,j) =  17  ! Water
       IF (soil_id .eq. 'R')  IBGPcl(i,j) =  11  ! Permanent wetland
       IF (soil_id .eq. 'P')  IBGPcl(i,j) =  10  ! Grassland
       IF (soil_id .eq. 'D')  IBGPcl(i,j) =  16  ! Barren or sparsely veg.
       IF (soil_id .eq. 'A')  IBGPcl(i,j) =  12  ! Croplands
       IF (soil_id .eq. 'W')  IBGPcl(i,j) =  7   ! Open shrublands
       IF (soil_id .eq. 'G')  IBGPcl(i,j) =  6   ! Closed shrublands
       IF (soil_id .eq. 'B')  IBGPcl(i,j) =  4   ! Deciduous broadleaf
       IF (soil_id .eq. 'N')  IBGPcl(i,j) =  1   ! Evergreen needleleaf
       IF (soil_id .eq. 'H')  IBGPcl(i,j) =  13  ! Urban and built-up
       IF (soil_id .eq. 'S')  IBGPcl(i,j) =  13  ! Urban and built-up

      ENDDO

      CLOSE(unit=10)


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Average of vegetation types for each grid cell
C +   ==============================================


      DO j=1,my
      DO i=1,mx

       dxi=ABS(NST__x(i,j)-NST__x(i-1,j-1))*111111.
     .                     *COS(degrad*NST__y(i,j))     
       dyi=ABS(NST__y(i,j)-NST__y(i-1,j-1))*111111.

       dist_max =SQRT(dxi*dxi+dyi*dyi)


C +---Seraching for the data grid point closest to the NST grid point
C +   ---------------------------------------------------------------

       dist_best=dist_max
       px       =0
       py       =0

       DO l=1,dim_y
       DO k=1,dim_x

        distx=(NST__x(i,j)-VEGlon(k,l))*111111.
     .                 *COS(degrad*NST__y(i,j))
        disty=(NST__y(i,j)-VEGlat(k,l))*111111.
        dist =SQRT(distx*distx+disty*disty)

        IF (dist.lt.dist_best) THEN
         dist_best=dist
         px       =k
         py       =l
        ENDIF 

       ENDDO
       ENDDO

       IF (dist_best.gt.dist_max) THEN
         px=0
         py=0
       ENDIF


C +---Compute an average for each vegetation types in the grid cell
C +   -------------------------------------------------------------

       IF (px.ne.0.and.py.ne.0.and.NSTsol(i,j).ge.4) THEN

C +---Initialization
C +   --------------

        nbx  =NINT(dxi/1000.)/2
        nby  =NINT(dyi/1000.)/2
        count=0

        DO k=1,nigbp
         IGBPfr(k)=0.
        ENDDO


C +---Vegetation types in the grid cell
C +   ---------------------------------

        DO l=MAX(py-nby,1),MIN(py+nby,dim_y)
        DO k=MAX(px-nbx,1),MIN(px+nbx,dim_x)
         count=count+1
         IGBPfr(IBGPcl(k,l))=IGBPfr(IBGPcl(k,l))+1.
        ENDDO
        ENDDO


C +---Percentages of each vegetation types
C +   ------------------------------------

        DO k=1,nigbp
         IF (count.gt.0) THEN
          IGBPfr(k)=IGBPfr(k)/MAX(1.,REAL(count))
         ELSE
          IGBPfr(k)=0.
         ENDIF
        ENDDO


C +---Retain three more important vegetation types
C +   --------------------------------------------

        DO k=1,nvx
         VEGcls(i,j,k)=0
        ENDDO
 
        DO k=1,nvx

         TMPuse=0
         TMPfrc=0.

         DO l=1,nigbp
          IF (l.ne.13 .and. l.ne.17 .and. IGBPfr(l).gt.TMPfrc) THEN
           TMPfrc=IGBPfr(l)
           TMPuse=l
          ENDIF
         ENDDO

         VEGcls(i,j,k)=TMPuse
         VEGfrc(i,j,k)=TMPfrc

        ENDDO


C +---Fraction of the three retained vegetation types
C +   -----------------------------------------------

        VEGtot=0
        DO k=1,nvx
         VEGtot=VEGtot+VEGcls(i,j,k)
        ENDDO

        DO k=1,nvx
         IF (VEGtot.gt.0.2) THEN
          VEGfrc(i,j,k)=REAL(VEGcls(i,j,k))/REAL(VEGtot)*100.
         ELSE
          VEGfrc(i,j,k)=0.
         ENDIF
        ENDDO


C +---Use of SVAT model ?
C +   -------------------

        IF (VEGtot.eq.0) THEN 
         NSTsol(i,j)=4   ! No
        ELSE
         NSTsol(i,j)=5   ! Yes
        ENDIF


       ENDIF

      ENDDO
      ENDDO


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Attribution of vegetation types to some NST grid points
C +   =======================================================

      DO j=1,my
      DO i=1,mx

       IF (NSTsol(i,j).eq.5) THEN
        DO k=1,nvx
         NSTveg(i,j,k)=VEGcls (i,j,k)
         NSTvfr(i,j,k)=NINT(VEGfrc(i,j,k))
        ENDDO
       ELSE
        DO k=1,nvx
         NSTveg(i,j,k)=0
         NSTvfr(i,j,k)=0.
        ENDDO
       ENDIF
 
      ENDDO
      ENDDO


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
