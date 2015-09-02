C   +-------------------------------------------------------------------+
C   |  File contents:                                                   |
C   |     INThor                                                        |
C   |     INThor1                                                       |
C   |     INTbic                                                        |
C   |     INTbil                                                        |
C   |     INTlin                                                        |
C   |     SPLINE, SPLINT, SPLIE2, SPLIN2 (from Numerical Recipies)      |
C   +-------------------------------------------------------------------+


C   +-------------------------------------------------------------------+
C   |  Subroutine INThor                               Dec. 95  NESTING |
C   |                                                  (Rev 2002 may)   |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Horizontal interpolation from LSC grid to NST grid distribute     |
C   | tasks to bicubic, linear... routines according to the "intype"    |
C   | variable (1= bilin, 3= bicub).                                    |
C   | Note that this routine uses the dimensions specified in NSTdim.inc|
C   | The bilinear interpolation is able to treat cyclic domains, or    |
C   | domains including the South/North pole.                           |
C   |                                                                   |
C   | Input : intype         : requested interpolation type             |
C   | ^^^^^^^ grd_Ix (ni,nj) : input grid points position x(i,j)        |
C   |         grd_Iy (ni,nj) : input grid points position y(i,j)        |
C   |         var_I  (ni,nj) : input field values                       |
C   |         grd_Ox (mx,my) : output grid positions x(i,j)             |
C   |         grd_Oy (mx,my) : output grid positions y(i,j)             |
C   |         SPHgrd (T/F)   : if true, spherical coordinates  for      |
C   |                          input fields                             |
C   |                                                                   |
C   | Output: var_O  (mx,my) : output field values                      |
C   | ^^^^^^^ pos_Ox (mx,my) : retained posit.for non-regular grid(long)|
C   |         pos_Oy (mx,my) : retained posit.for non-regular grid (lat)|
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE INThor (intype,grd_Ix,grd_Iy,var_I,
     .                   SPHgrd,grd_Ox,grd_Oy,var_O,
     .                   REGgrd,pos_Ox,pos_Oy)
 


C +---LSC and NST domain dimensions
C +   -----------------------------

      include 'NSTdim.inc'


C +---Local variables
C +   ---------------

      INTEGER intype,i,j

      INTEGER pos_Ox(mx,my),pos_Oy(mx,my)

      REAL grd_Ix(ni,nj),grd_Iy(ni,nj),var_I(ni,nj),
     .     grd_Ox(mx,my),grd_Oy(mx,my),var_O(mx,my)

      LOGICAL SPHgrd,REGgrd


C +---Temporary arrays
C +   ----------------

      REAL tmp_I2a(ni,nj),tmp1in(ni,nj),tmp2in(0:ni+1,0:nj+1),
     .     grd1Ix(ni),grd1Iy(nj)

C +---Interpolation
C +   -------------

      IF (REGgrd) THEN            ! Regular input grid

       DO i=1,ni
        grd1Ix(i)=grd_Ix(i,1)
       ENDDO

       DO j=1,nj
        grd1Iy(j)=grd_Iy(1,j)
       ENDDO

       IF (intype.EQ.1) THEN      ! Bilinear interpolation

C +           ******
         CALL INTbil (ni,nj,grd1Ix,grd1Iy,var_I,SPHgrd,
     .                mx,my,grd_Ox,grd_Oy,var_O,tmp2in)
C +           ******

       ELSE IF (intype.EQ.3) THEN ! Bicubic interpolation

C +           ******
         CALL INTbic (tmp_I2a,tmp1in,
     .                ni,nj,grd1Ix,grd1Iy,var_I,
     .                mx,my,grd_Ox,grd_Oy,var_O)
C +           ******

       ENDIF

      ELSE                        ! Non-regular input grid

C +         ******
       CALL INTnrg2 (ni,nj,grd_Ix,grd_Iy,var_I,
     .              mx,my,grd_Ox,grd_Oy,var_O,
     .              pos_Ox,pos_Oy)
C +         ******

      ENDIF


      RETURN
      END


C   +-------------------------------------------------------------------+
C   |  Subroutine INThorV                                               |
C   |                                                  (Rev 2002 june)  |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | A variant of INThor adapted to staggered grids                    |
C   |   which have a different lat size, namely njv instead of nj       |
C   |   (this should not exist if INTERp is rewriten in f90 !)          |
C   |                                                                   |
C   | Input : intype         : requested interpolation type             |
C   | ^^^^^^^ grd_Ix (ni,njv) : input grid points position x(i,j)       |
C   |         grd_Iy (ni,njv) : input grid points position y(i,j)       |
C   |         var_I  (ni,njv) : input field values                      |
C   |         grd_Ox (mx,my) : output grid positions x(i,j)             |
C   |         grd_Oy (mx,my) : output grid positions y(i,j)             |
C   |         SPHgrd (T/F)   : if true, spherical coordinates for       |
C   |                          input fields                             |
C   |                                                                   |
C   | Output: var_O  (mx,my) : output field values                      |
C   | ^^^^^^^ pos_Ox (mx,my) : retained posit.for non-regular grid(long)|
C   |         pos_Oy (mx,my) : retained posit.for non-regular grid (lat)|
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE INThorV(intype,grd_Ix,grd_Iy,var_I,
     .                   SPHgrd,grd_Ox,grd_Oy,var_O,
     .                   REGgrd,pos_Ox,pos_Oy)
 


C +---LSC and NST domain dimensions
C +   -----------------------------

      include 'NSTdim.inc'


C +---Local variables
C +   ---------------

      INTEGER intype,i,j

      INTEGER pos_Ox(mx,my),pos_Oy(mx,my)

      REAL grd_Ix(ni,njv),grd_Iy(ni,njv),var_I(ni,njv),
     .     grd_Ox(mx,my) ,grd_Oy(mx,my) ,var_O(mx,my)

      LOGICAL SPHgrd,REGgrd


C +---Temporary arrays
C +   ----------------

      REAL tmp_I2a(ni,njv),tmp1in(ni,njv),tmp2in(0:ni+1,0:njv+1),
     .     grd1Ix(ni),grd1Iy(njv)

C +---Interpolation
C +   -------------

      IF (REGgrd) THEN            ! Regular input grid

       DO i=1,ni
        grd1Ix(i)=grd_Ix(i,1)
       ENDDO

       DO j=1,njv
        grd1Iy(j)=grd_Iy(1,j)
       ENDDO

       IF (intype.EQ.1) THEN      ! Bilinear interpolation

C +           ******
         CALL INTbil (ni,njv,grd1Ix,grd1Iy,var_I,SPHgrd,
     .                mx,my,grd_Ox,grd_Oy,var_O,tmp2in)
C +           ******

       ELSE IF (intype.EQ.3) THEN ! Bicubic interpolation

C +           ******
         CALL INTbic (tmp_I2a,tmp1in,
     .                ni,njv,grd1Ix,grd1Iy,var_I,
     .                mx,my,grd_Ox,grd_Oy,var_O)
C +           ******

       ENDIF

      ELSE                        ! Non-regular input grid

         write(*,*) 'Non-regular grid  '
         write(*,*) '  + staggered :   '
         write(*,*) 'This is not possible with current NESTOR'
         STOP

      ENDIF


      RETURN
      END

C   +-------------------------------------------------------------------+
C   |  Subroutine INThor1                         February 2000 NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Horizontal interpolation, ...                       distribute    |
C   | tasks to bicubic, linear... routines according to the "intype"    |
C   | variable (1= bilin, 3= bicub).                                    |
C   | Note that this routine uses the dimensions specified in NSTdim.inc|
C   | The bilinear interpolation is able to treat cyclic domains, or    |
C   | domains including the South/North pole.                           |
C   |                                                                   |
C   | Input : intype         : requested interpolation type             |
C   | ^^^^^^^ grd_Ix (mx,my) : input grid points position x(i,j)        |
C   |         grd_Iy (mx,my) : input grid points position y(i,j)        |
C   |         var_I  (mx,my) : input field values                       |
C   |         grd_Ox ( 1, 1) : output grid positions x(i,j)             |
C   |         grd_Oy ( 1, 1) : output grid positions y(i,j)             |
C   |         SPHgrd (T/F)   : if true, spherical coordinates for       |
C   |                          input fields                             |
C   |                                                                   |
C   | Output: var_O  ( 1, 1) : output field values                      |
C   | ^^^^^^^ pos_Ox ( 1, 1) : retained posit.for non-regular grid(long)|
C   |         pos_Oy ( 1, 1) : retained posit.for non-regular grid (lat)|
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE INThor1 (intype,grd_Ix,grd_Iy,var_I,
     .                    SPHgrd,grdsOx,grdsOy,varsO,
     .                    REGgrd,possOx,possOy)
 


C +---LSC and NST domain dimensions
C +   -----------------------------

      include 'NSTdim.inc'


C +---Local variables
C +   ---------------

      INTEGER intype,nx,ny,i,j

      PARAMETER (nx=1)
      PARAMETER (ny=1)

      INTEGER pos_Ox(nx,ny),pos_Oy(nx,ny),possOx,possOy

      REAL grd_Ix(mx,my),grd_Iy(mx,my),var_I(mx,my),
     .     grd_Ox(nx,ny),grd_Oy(nx,ny),var_O(nx,ny),
     .     grdsOx       ,grdsOy       ,varsO

      LOGICAL SPHgrd,REGgrd


C +---Temporary arrays
C +   ----------------

      REAL tmp_I2a(mx,my),tmp1in(mx,my),tmp2in(0:mx+1,0:my+1),
     .     grd1Ix(mx),grd1Iy(my)


C +---Interpolation
C +   -------------

      grd_Ox(1,1) = grdsOx
      grd_Oy(1,1) = grdsOy
      pos_Ox(1,1) = possOx
      pos_Oy(1,1) = possOy

      IF (REGgrd) THEN            ! Regular input grid

       DO i=1,mx
        grd1Ix(i)=grd_Ix(i,1)
       ENDDO

       DO j=1,my
        grd1Iy(j)=grd_Iy(1,j)
       ENDDO

       IF (intype.EQ.1) THEN      ! Bilinear interpolation

C +           ******
         CALL INTbil (mx,my,grd1Ix,grd1Iy,var_I,SPHgrd,
     .                nx,ny,grd_Ox,grd_Oy,var_O,tmp2in)
C +           ******

       ELSE IF (intype.EQ.3) THEN ! Bicubic interpolation

C +           ******
         CALL INTbic (tmp_I2a,tmp1in,
     .                mx,my,grd1Ix,grd1Iy,var_I,
     .                nx,ny,grd_Ox,grd_Oy,var_O)
C +           ******

       ENDIF

      ELSE                        ! Non-regular input grid

C +         ******
       CALL INTnrg2 (mx,my,grd_Ix,grd_Iy,var_I,
     .              nx,ny,grd_Ox,grd_Oy,var_O,
     .              pos_Ox,pos_Oy)
C +         ******

      ENDIF

      varsO  = var_O (1,1)
      possOx = pos_Ox(1,1)
      possOy = pos_Oy(1,1)


      RETURN
      END


C   +-------------------------------------------------------------------+
C   |  Subroutine INTbic                               Dec. 96  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | This routine computes a natural bicubic spline interpolation of   |
C   | 2D scalar field, from a RECTANGULAR grid to an IRREGULAR grid,    |
C   | i.e.    var_I(i,j) (coarse) ===>  var_O(ii,jj) (fine mesh)        |
C   |                                                                   |
C   | Note 1: The independant variable should be adequatelly chosen     |
C   |         to have the "regularity" in the INPUT grid                |
C   | Note 2: This routine is simple but not efficient; a faster        |
C   |         routine is available for 2 rectangular grids, and         |
C   |         modifications to INTbic possible if necessary.            |
C   |                                                                   |
C   | References : W.H.Pres et al., Numerical recipes,Cambridge UP.1992 |
C   | ^^^^^^^^^^^^ Lancaster P., Curve and Surface fitting, AP. 1986    |
C   |                                                                   |
C   | INPUT : grd_Ix (dim_Ix)         : Input grid points position x(i) |
C   | ^^^^^^^ grd_Iy (dim_Iy)         :   "     "    "       "     y(j) |
C   |         var_I  (dim_Ix, dim_Iy) : Input field values              |
C   |         grd_Ox (dim_Ox, dim_Oy) : Output grid positions x(i,j)    |
C   |         grd_Oy (dim_Ox, dim_Oy) : Output grid positions y(i,j)    |
C   |                                                                   |
C   |         dim_Ix, dim_Iy : (parameter type) array dimensions        |
C   |                          ! = size of the interpolated data field  |
C   |         dim_Ox, dim_Oy : (parameter type) array dimensions        |
C   |                          = size of the used data field            |
C   |                                                                   |
C   | OUTPUT: var_O  (dim_Ox, dim_Oy) : Output field values             |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | TEMPORARY arrays : tmp_I2a (dim_Ix, dim_Iy)                       |
C   | ^^^^^^^^^^^^^^^^^^ tmp_in  (dim_Ix, dim_Iy)                       |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE INTbic(tmp_I2a,tmp_in,
     .                  dim_Ix,dim_Iy,grd_Ix,grd_Iy,var_I,
     .                  dim_Ox,dim_Oy,grd_Ox,grd_Oy,var_O)
 

      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INTEGER dim_Ix,dim_Iy,dim_Ox,dim_Oy,i,j,ii,jj

      REAL grd_Ix(dim_Ix),grd_Iy(dim_Iy),var_I(dim_Ix,dim_Iy),
     .     grd_Ox(dim_Ox,dim_Oy),grd_Oy(dim_Ox,dim_Oy),
     .     var_O (dim_Ox,dim_Oy),tmp_in(dim_Ix,dim_Iy),
     .     tmp_I2a(dim_Ix,dim_Iy),tmp_Ix(500),tmp_Iy(500)
 

 
C +-- Coordinates indexes inversion (if necessary)
C +   ============================================
 
C +...Note : unmodified routines from Num. Recipes
C +...       implies that coords are such as
C +...       x(1) < x(2) <...  ===> revert if necessary


C +---Storage in temporary arrays
C +   ---------------------------
 
      DO jj=1,dim_Iy
      DO ii=1,dim_Ix
       tmp_in(ii,jj)=var_I(ii, jj)
      ENDDO
      ENDDO

      DO ii=1,dim_Ix
       tmp_Ix(ii)=grd_Ix(ii)
      ENDDO

      DO jj=1,dim_Iy
       tmp_Iy(jj)=grd_Iy(jj)
      ENDDO
 

C +---Revert grd_Ix (1) <--> grd_Ix (n), ... ?
C +   ----------------------------------------

      IF (grd_Ix(dim_Ix).lt.grd_Ix(1)) THEN
       DO ii=1,dim_Ix
        DO jj=1,dim_Iy
         tmp_in(ii,jj)=var_I(dim_Ix-ii+1,jj)
        ENDDO
        tmp_Ix(ii)=grd_Ix(dim_Ix-ii+1)
       ENDDO
      ENDIF

 
C +---Revert grd_Iy (1) <--> grd_Iy (n), ... ?
C +   ----------------------------------------
 
      IF (grd_Iy(dim_Iy).lt.grd_Iy(1)) THEN
       DO jj=1,dim_Iy
        DO ii=1,dim_Ix
         tmp_in(ii,jj)=var_I(ii,dim_Iy-jj+1)
        ENDDO
        tmp_Iy(jj)=grd_Iy(dim_Iy-jj+1)
       ENDDO
      ENDIF

 
C +---Construction of 1D splines of rows
C +   ==================================
 
C +        ****** 
      CALL SPLIE2(tmp_Ix,tmp_Iy,tmp_in,dim_Ix,dim_Iy,tmp_I2a)
C +        ****** 

 
C +---Interpolation for each output grid point
C +   ========================================
 
      DO j=1,dim_Oy
      DO i=1,dim_Ox

C +         ****** 
       CALL SPLIN2(tmp_Ix,tmp_Iy,tmp_in,tmp_I2a,dim_Ix,dim_Iy,
     .             grd_Ox(i,j),grd_Oy(i,j),var_O(i,j))
C +         ****** 
 
      ENDDO
      ENDDO
 
 
      RETURN
      END


C   +-------------------------------------------------------------------+
C   |  Subroutine INTbil                            01-07-2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | This routine is a bilinear interpolation of a 2D scalar fields.   |
C   | If the output resolution is lower than input, an average of 5     |
C   | bilinear interpolations is performed, considering 5 sampling      |
C   | points located around the selected point in the output mesh.      |
C   | Note that a specific treatment of latitudes/longitudes is         |
C   | included for input grids using spherical coordinates.             |
C   |                                                                   |
C   | Input : grd_Ix (ni)     : Input grid points position x(i)         |
C   | ^^^^^^^ grd_Iy (nj)     :   "     "    "       "     y(j)         |
C   |         var_I  (ni, nj) : Input field values                      |
C   |         grd_Ox (mx, my) : Output grid positions x(i,j)            |
C   |         grd_Oy (mx, my) : Output grid positions y(i,j)            |
C   |         SPHgrd (T/F)    : If true, spherical coordinates for      |
C   |                           input fields                            |
C   |                                                                   |
C   | Output: var_O  (mx, my) : Output field values                     |
C   | ^^^^^^^                                                           |
C   +-------------------------------------------------------------------+


      SUBROUTINE INTbil (ni,nj,grd_Ix,grd_Iy,var_I,SPHgrd,
     .                   mx,my,grd_Ox,grd_Oy,var_O,tmp_in)


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
     .     AUXla1,AUXla2,dx,dy,degrad

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
        WRITE(6,*) 'INTbil - fatal error: dimension',LocDim
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
      dist_I=SQRT(dx*dx+dy*dy)

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
      dist_O=SQRT(dx*dx+dy*dy)

      IF (dist_I.lt.dist_O) THEN
       nsamp=ns
      ELSE
       nsamp=1
      ENDIF


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
        WRITE(6,*) '-180 and +180 deg. (as required by INTbil)'
        WRITE(6,*) 'but rather between : '
        WRITE(6,*) tmp_Ix(1),tmp_Ix(ni)
        WRITE(6,*) '--- STOP in INTbil ---'
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
       WRITE(6,*) '--- STOP in INTbil ---'
      ENDIF


C +---Bi-linear interpolation
C +   =======================


C +---Some initialisations
C +   --------------------

      p=0
      q=0

      i1(1)=-1  ; j1(1)=+1
      i1(2)=+1  ; j1(2)=-1
      i1(3)=-1  ; j1(3)=-1
      i1(4)=+1  ; j1(4)=+1

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO i=1,mx   ! LOOP on output grid-points : BEGIN
      DO j=1,my  
          
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Define sampling point positions
C +   -------------------------------

       IF (i.ne.1.and.i.ne.mx.and.j.ne.1.and.j.ne.my) THEN

        samOx(1)= grd_Ox(i,j)
        samOy(1)= grd_Oy(i,j)

        do i2=1,4
         samOx(i2+1)=0.6*grd_Ox(i,j)+0.4*grd_Ox(i+i1(i2),j+j1(i2))
         samOy(i2+1)=0.6*grd_Oy(i,j)+0.4*grd_Oy(i+i1(i2),j+j1(i2))

         if(sign(1.,grd_Ox(i,j))              .ne.
     .      sign(1.,grd_Ox(i+i1(i2),j+j1(i2))).and.
     .       abs(   grd_Ox(i,j))              .ge. 170.0 ) then
          samOx(i2+1)=grd_Ox(i,j)
         endif

        enddo

       ELSE
        DO is=1,nsamp   ! Boundaries : no sampling
         samOx(is)=grd_Ox(i,j)
         samOy(is)=grd_Oy(i,j)
        ENDDO
       ENDIF


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
         WRITE (6,*) '--- STOP in INTbil ---'
         STOP
        ENDIF


C +---Linear interpolation
C +   --------------------

        x0 = tmp_Ix(p)
        x1 = tmp_Ix(p+1)
        y0 = tmp_Iy(q)
        y1 = tmp_Iy(q+1)
        tmp=(x-x0)*((y-y0)*tmp_in(p+1,q+1)
     .             +(y1-y)*tmp_in(p+1,q  ))
     .     +(x1-x)*((y-y0)*tmp_in(p  ,q+1)
     .             +(y1-y)*tmp_in(p  ,q  ))
        tmp2=tmp2+tmp/((x1-x0)*(y1-y0))


       ENDDO          ! LOOP on the sampling points: END


C +---Output value given by the average of the samplings
C +   --------------------------------------------------

       var_O(i,j)=tmp2/REAL(nsamp)


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDDO
      ENDDO           ! Loop on output grid-points : END

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c #MS IF (cyclic) WRITE(6,*) 'INTbil-info: cyclic boundaries'
c #MS IF (npole ) WRITE(6,*) 'INTbil-info: North pole included'
c #MS IF (spole ) WRITE(6,*) 'INTbil-info: South pole included'


      RETURN
      END


C   +-------------------------------------------------------------------+
C   |  Subroutine INTnrg                             June 2000  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | This routine is a linear interpolation of a 2D scalar fields from |
C   | a non-regular grid to a non-regular grid.                         |
C   |                                                                   |
C   | Input : grd_Ix (ni, nj) : Input grid points position x(i,j)       |
C   | ^^^^^^^ grd_Iy (ni, nj) :   "     "    "       "     y(i,j)       |
C   |         var_I  (ni, nj) : Input field values                      |
C   |         grd_Ox (mx, my) : Output grid positions x(i,j)            |
C   |         grd_Oy (mx, my) : Output grid positions y(i,j)            |
C   |                                                                   |
C   | Output: var_O  (mx, my) : Output field values                     |
C   | ^^^^^^^                                                           |
C   +-------------------------------------------------------------------+


      SUBROUTINE INTnrg (ni,nj,grd_Ix,grd_Iy,var_I,
     .                   mx,my,grd_Ox,grd_Oy,var_O,
     .                   pos_Ox,pos_Oy)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INTEGER  i,j,k,l,ni,nj,mx,my,k1,k2,l1,l2

      INTEGER pos_Ox(mx,my),pos_Oy(mx,my)

      REAL clon,clat,lat1,lat2,lat3,lat4,lon1,lon2,lon3,lon4,val1,
     .     val2,val3,val4,tmp,dlon12,dlon14,dlat12,dlat14,ilat,ilon,
     .     x1,x2,x3,y1,y2,y3,z1,z2,z3,cfa,cfb,cfc,val123,denom2,
     .     val214,val314,val423,d1,d2,d3,d4,dtot,epsi,eps6,denom1

      REAL grd_Ix(ni,nj),grd_Iy(ni,nj),grd_Ox(mx,my),grd_Oy(mx,my),
     .     var_I (ni,nj),var_O (mx,my)

      LOGICAL includ,found,krig,inc123,inc124

      DATA epsi / 1.d-5 /
      DATA eps6 / 1.d-6 /


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO i=1,mx   ! LOOP on output grid-points : BEGIN
      DO j=1,my  
          
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       var_O(i,j)=0.

       clon=grd_Ox(i,j)
       clat=grd_Oy(i,j)

       IF (pos_Ox(i,j).gt.0.and.pos_Ox(i,j).lt.ni .and.
     .     pos_Oy(i,j).gt.0.and.pos_Oy(i,j).lt.nj) THEN
        k1=pos_Ox(i,j)
        k2=pos_Ox(i,j)
        l1=pos_Oy(i,j)
        l2=pos_Oy(i,j)
        found=.true.
        includ=.true.
       ELSE
        k1=1
        k2=ni-1
        l1=1
        l2=nj-1
        found=.false.
        includ=.false.
       ENDIF

       DO l=l1,l2    ! LOOP on input grid-points : BEGIN
       DO k=k1,k2

        lon1=grd_Ix(k  ,l  )
        lon2=grd_Ix(k+1,l  )
        lon3=grd_Ix(k+1,l+1)
        lon4=grd_Ix(k  ,l+1)

        lat1=grd_Iy(k  ,l  )
        lat2=grd_Iy(k+1,l  )
        lat3=grd_Iy(k+1,l+1)
        lat4=grd_Iy(k  ,l+1)

        val1=var_I (k  ,l  )
        val2=var_I (k+1,l  )
        val3=var_I (k+1,l+1)
        val4=var_I (k  ,l+1)


C +---Rotation of input grid cell ?
C +   =============================

        dlat12=ABS(lat2-lat1)
        dlat14=ABS(lat4-lat1)
        dlon12=ABS(lon2-lon1)
        dlon14=ABS(lon4-lon1)

        IF ((dlat14.lt.dlat12).and.(dlon14.gt.dlon12)) THEN
         tmp =lat1   ! Latitude
         lat1=lat2
         lat2=lat3
         lat3=lat4
         lat4=tmp
         tmp =lon1   ! Longitude
         lon1=lon2
         lon2=lon3
         lon3=lon4
         lon4=tmp
         tmp =val1   ! Values
         val1=val2
         val2=val3
         val3=val4
         val4=tmp
        ENDIF


C +---Invert latitude ?
C +   =================

        IF (lat4.lt.lat1) THEN

         IF (lat3.lt.lat2) THEN
          tmp =lat2   ! Latitude
          lat2=lat3
          lat3=tmp
          tmp =lat1
          lat1=lat4
          lat4=tmp
          tmp =lon2   ! Longitude
          lon2=lon3
          lon3=tmp
          tmp =lon1
          lon1=lon4
          lon4=tmp
          tmp =val2   ! Values
          val2=val3
          val3=tmp
          tmp =val1
          val1=val4
          val4=tmp
         ELSE
          WRITE(6,*) 'Inconsistance in latitude. Cannot be resolved'
          WRITE(6,*) 'by INTnrg subroutine.                --- STOP'
          WRITE(6,*)
          WRITE(6,*) 'Info : ',lat1,lat2,lat3,lat4
          STOP
         ENDIF

        ELSE

         IF (lat3.lt.lat2) THEN

          WRITE(6,*) 'Inconsistance in latitude. Cannot be resolved'
          WRITE(6,*) 'by INTnrg subroutine.                --- STOP'
          WRITE(6,*)
          WRITE(6,*) 'Info : ',lat1,lat2,lat3,lat4
          STOP

         ENDIF

        ENDIF


C +---Invert longitude ?
C +   ==================

        IF (lon2.lt.lon1) THEN

         IF (lon3.lt.lon4) THEN
          tmp =lat3   ! Latitude
          lat3=lat4
          lat4=tmp
          tmp =lat1
          lat1=lat2
          lat2=tmp
          tmp =lon3   ! Longitude
          lon3=lon4
          lon4=tmp
          tmp =lon1
          lon1=lon2
          lon2=tmp
          tmp =val3   ! Values
          val3=val4
          val4=tmp
          tmp =val1
          val1=val2
          val2=tmp
         ELSE
          WRITE(6,*) 'Inconsistance in longitude. Cannot be resolved'
          WRITE(6,*) 'INTnrg subroutine.                    --- STOP'
          WRITE(6,*)
          WRITE(6,*) 'Info : ',lon1,lon2,lon3,lon4
          STOP
         ENDIF

        ELSE

         IF (lon3.lt.lon4) THEN

          WRITE(6,*) 'Inconsistance in longitude. Cannot be resolved'
          WRITE(6,*) 'INTnrg subroutine.                    --- STOP'
          WRITE(6,*)
          WRITE(6,*) 'Info : ',lon1,lon2,lon3,lon4
          STOP

         ENDIF

        ENDIF


C +   At this stage, it is assumed that the input grid cell is
C +   such that :
C +
C +              4----------3
C +              |          |
C +              |          |
C +              |          |
C +              1----------2
C +
C +   with lon1 < lon2, lon4 < lon3
C +        lat1 < lat4, lat2 < lat3


C +---Check if input cell includes output grid point
C +   ==============================================


        IF (.not.found) THEN

         includ=.true.

C +---Segment 1-2
C +   -----------

         ilat=lat1+(clon-lon1)/(lon2-lon1)*(lat2-lat1)
         IF (ilat.gt.clat) includ=.false.

C +---Segment 4-3
C +   -----------

         ilat=lat4+(clon-lon4)/(lon3-lon4)*(lat3-lat4)
         IF (ilat.lt.clat) includ=.false.

C +---Segment 2-3
C +   -----------

         ilon=lon2+(clat-lat2)/(lat3-lat2)*(lon3-lon2)
         IF (ilon.lt.clon) includ=.false.

C +---Segment 1-4
C +   -----------

         ilon=lon1+(clat-lat1)/(lat4-lat1)*(lon4-lon1)
         IF (ilon.gt.clon) includ=.false.

C +---Have found an input cell ?
C +   --------------------------

         found=includ

        ENDIF


C +---Interpolation
C +   =============

        IF (includ.and.found) THEN

C +---Initialization
C +   --------------

         krig  =.false.
         val123= 0.
         val214= 0.
         val314= 0.
         val423= 0.


C +---Retain position on input grid
C +   -----------------------------

         pos_Ox(i,j)=k
         pos_Oy(i,j)=l
 

C +---Check the inclusion in triangle 1-3-4 or 1-2-3
C +   ----------------------------------------------

         ilat=lat1+(clon-lon1)/(lon3-lon1)*(lat3-lat1)
         IF (ilat.gt.clat) THEN
          inc123=.true.      ! in triangle 1-2-3
         ELSE
          inc123=.false.     ! in triangle 1-3-4
         ENDIF


C +---Check the inclusion in triangle 1-2-4 or 2-3-4
C +   ----------------------------------------------

         ilat=lat4+(clon-lon4)/(lon2-lon4)*(lat2-lat4)
         IF (ilat.gt.clat) THEN
          inc124=.true.      ! in triangle 1-2-4
         ELSE
          inc124=.false.     ! in triangle 2-3-4
         ENDIF


C +---Solution : plane 1, 2, 3
C +   ------------------------

         IF (inc123) THEN

          x1 = lon1
          x2 = lon2
          x3 = lon3
          y1 = lat1
          y2 = lat2
          y3 = lat3
          z1 = val1
          z2 = val2
          z3 = val3

          denom1 = (y3-y2)
          denom2 = (x1-x3)+(x3-x2)/(y3-y2)*(y1-y3)

          IF (ABS(denom1).lt.eps6.or.ABS(denom2).lt.eps6) krig=.true.

          cfa = ( (z1-z3) + (z3-z2)/(y3-y2)*(y3-y1) ) / denom2
          cfb = ( cfa*(x2-x3) + (z3-z2) )             / denom1
          cfc = z3 - cfa*x3 - cfb*y3

          val123 = cfa*clon + cfb*clat + cfc

         ENDIF


C +---Solution : plane 2, 1, 4
C +   ------------------------

         IF (inc124) THEN

          x1 = lon2
          x2 = lon1
          x3 = lon4
          y1 = lat2
          y2 = lat1
          y3 = lat4
          z1 = val2
          z2 = val1
          z3 = val4
 
          denom1 = (y3-y2)
          denom2 = (x1-x3)+(x3-x2)/(y3-y2)*(y1-y3)
 
          IF (ABS(denom1).lt.eps6.or.ABS(denom2).lt.eps6) krig=.true.

          cfa = ( (z1-z3) + (z3-z2)/(y3-y2)*(y3-y1) ) / denom2
          cfb = ( cfa*(x2-x3) + (z3-z2) )             / denom1
          cfc = z3 - cfa*x3 - cfb*y3

          val214 = cfa*clon + cfb*clat + cfc

         ENDIF


C +---Solution : plane 3, 1, 4
C +   ------------------------

         IF (.not.inc123) THEN

          x1 = lon3
          x2 = lon1
          x3 = lon4
          y1 = lat3
          y2 = lat1
          y3 = lat4
          z1 = val3
          z2 = val1
          z3 = val4
 
          denom1 = (y3-y2)
          denom2 = (x1-x3)+(x3-x2)/(y3-y2)*(y1-y3)

          IF (ABS(denom1).lt.eps6.or.ABS(denom2).lt.eps6) krig=.true.

          cfa = ( (z1-z3) + (z3-z2)/(y3-y2)*(y3-y1) ) / denom2
          cfb = ( cfa*(x2-x3) + (z3-z2) )             / denom1
          cfc = z3 - cfa*x3 - cfb*y3

          val314 = cfa*clon + cfb*clat + cfc

         ENDIF


C +---Solution : plane 4, 2, 3
C +   ------------------------

         IF (.not.inc124) THEN

          x1 = lon4
          x2 = lon2
          x3 = lon3
          y1 = lat4
          y2 = lat2
          y3 = lat3
          z1 = val4
          z2 = val2
          z3 = val3
 
          denom1 = (y3-y2)
          denom2 = (x1-x3)+(x3-x2)/(y3-y2)*(y1-y3)
 
          IF (ABS(denom1).lt.eps6.or.ABS(denom2).lt.eps6) krig=.true.

          cfa = ( (z1-z3) + (z3-z2)/(y3-y2)*(y3-y1) ) / denom2
          cfb = ( cfa*(x2-x3) + (z3-z2) )             / denom1
          cfc = z3 - cfa*x3 - cfb*y3

          val423 = cfa*clon + cfb*clat + cfc

         ENDIF


C +---Global solution
C +   ---------------

         IF (krig) THEN

          d1=MAX(epsi,(lon1-clon)**2.+(lat1-clat)**2.)
          d2=MAX(epsi,(lon2-clon)**2.+(lat2-clat)**2.)
          d3=MAX(epsi,(lon3-clon)**2.+(lat3-clat)**2.)
          d4=MAX(epsi,(lon4-clon)**2.+(lat4-clat)**2.)
          dtot=1./d1+1./d2+1./d3+1./d4
          var_O(i,j)=(val1/d1+val2/d2+val3/d3+val4/d4)/dtot

         ELSE

          var_O(i,j)=0.5*(val123+val214+val314+val423)

         ENDIF

         includ=.false.

        ENDIF         ! { IF (includ.and.found) }


       ENDDO          ! Loop on input grid-points : END
       ENDDO


       IF (.not.found) THEN
        WRITE(6,*) 'No cell of input grid includes an output grid'
        WRITE(6,*) 'point.                     --- STOP in INTnrg.'
        WRITE(6,*)
        WRITE(6,*) 'Info : ',i,j,clon,clat
        STOP
       ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDDO
      ENDDO           ! Loop on output grid-points : END

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END



C   +-------------------------------------------------------------------+
C   |  Subroutine SPHERC                               July 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | This routine sets longitude between -180 and +180, and latitude   |
C   | between -90 and +90, as required by some interpolation sub-       |
C   | routines.                                                         |
C   |                                                                   |
C   | Input : LONval : longitude   (degree)                             |
C   | ^^^^^^^ LATval : latitude    (degree)                             |
C   |         SPHgrd : If true, LONval and LATval really are spherical  |
C   |                  coordinates                                      |
C   |                                                                   |
C   | Output: LONval : longitude   (degree)                             |
C   | ^^^^^^^ LATval : latitude    (degree)                             |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE SPHERC (SPHgrd,LONval,LATval)


      IMPLICIT NONE

      REAL    LONval,LATval
      LOGICAL SPHgrd


      IF (SPHgrd) THEN

C +---Longitude defined between -180 and +180 degree
C +   ----------------------------------------------

       IF (LONval.ge.  180. ) LONval=LONval-360.
       IF (LONval.lt.(-180.)) LONval=LONval+360.

C +---Latitude defined between -90 and +90 degree
C +   -------------------------------------------

       IF (LATval.ge.   90. ) LATval=LATval-180.
       IF (LATval.lt. (-90.)) LATval=LATval+180.

      ENDIF


      RETURN
      END


C     +--------------------------------------------------------------+
      subroutine INTlin (xx,vv,ni,xreq,outvar)   ! Last modif. : 04/99
C     +--------------------------------------------------------------+

      REAL xx(ni), vv(ni)
      REAL xreq, outvar, fdis
      INTEGER ind, KLO, KHI,ni,k

      KLO=1
      KHI=ni
 1    IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(xx(K).GT.xreq)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      ind=KLO

      fdis = xx(ind+1)-xx(ind)
      outvar= vv(ind)*((xx(ind+1)-xreq)/fdis)                 
     .    + vv(ind+1)*((xreq-xx(ind  ))/fdis)                 

      IF (xreq.LT.xx(ind  )) outvar=vv(ind  )
      IF (xreq.GT.xx(ind+1)) outvar=vv(ind+1)


      RETURN
      END

C     +--------------------------------------------------------------+
      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
C     |     * From numerical recipes (H. Press et al., 1992)         |
C     +--------------------------------------------------------------+
 
      PARAMETER (NN=500)
C +...NN = max value allowed for N (for 1D arrays only -> overdim.)
 
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N)
      DIMENSION YTMP(NN),Y2TMP(NN)
 
      DO 13 J=1,M
 
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
 
C +          ******
        CALL SPLINE(X2A,YTMP,N,1.E30,1.E30,Y2TMP)
C +          ******
C +...  NB : 1.E30 = switch value to select "natural" bicub spline
 
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
 
13    CONTINUE
 
      RETURN
      END
 
C     +--------------------------------------------------------------+
      subroutine SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
C     |     * From numerical recipes (H. Press et al., 1992)         |
C     |       Version grossierement optimisee                        |
C     +--------------------------------------------------------------+
 
      PARAMETER (NN=500)
C +   NN = max value allowed for N (for 1D arrays only -> overdim.)
 
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)
      DIMENSION YYTMP(NN)

      KLO=1
      KHI=N
31    IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(X2A(K).GT.X2) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
          GOTO 31
      ENDIF
      H=X2A(KHI)-X2A(KLO)
      IF (H.EQ.0.) THEN
        write(*,*) 'INTERp/SPLIN2: error (Bad XA input)'
        STOP
      ENDIF

      A=(X2A(KHI)-X2)/H
      B=(X2-X2A(KLO))/H

      DO J=1,M
 
        YYTMP(J)=A*YA(J,KLO)+B*YA(J,KHI)+
     .      ((A**3-A)*Y2A(J,KLO)+(B**3-B)*Y2A(J,KHI))*(H**2)/6.

      ENDDO
 
C +        ******
      CALL SPLINE(X1A,YYTMP,M,1.E30,1.E30,Y2TMP)
C +        ******
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)
C +        ******
 
      RETURN
      END
 
C     +--------------------------------------------------------------+
      subroutine SPLINE(X,Y,N,YP1,YPN,Y2)
C     |     * From numerical recipes (H. Press et al., 1992)         |
C     +--------------------------------------------------------------+
 
      PARAMETER (NMAX=500)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
 
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
 
C     +--------------------------------------------------------------+
      subroutine SPLINT(XA,YA,Y2A,N,X,Y)
C     |     * From numerical recipes (H. Press et al., 1992)         |
C     +--------------------------------------------------------------+
 
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
 
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
 
      RETURN
      END

