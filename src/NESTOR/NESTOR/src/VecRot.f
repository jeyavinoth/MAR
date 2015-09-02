C   +-------------------------------------------------------------------+
C   |  Subroutine VecRot                               Sept 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Rotates vector (wind) following the rotation of the reference     |
C   | frame that occurs on a map projection.                            |
C   | A special treatment is done for grid points close to N/S pole.    |
C   |                                                                   |
C   | References : Map projections and equations of motion suitable     |
C   | ^^^^^^^^^^^^ for mesoscale alpha simulations by the MAR model.    |
C   |                                                                   |
C   | Input : - grd_lon (mx, my) : grid positions lon(i,j)              |
C   | ^^^^^^^ - grd_lat (mx, my) : grid positions lat(i,j)              |
C   |         - dx : mesh size   : the grid is assumed to be square     |
C   |                              and fully regullar                   |
C   |         - var_1   (mx, my) : "x" component of the vector          |
C   |         - var_2   (mx, my) : "y" component of the vector          |
C   |                              (local cartesian on the sphere)      |
C   |                                                                   |
C   | Output: - var_1   (mx, my) : x component of the vector            |
C   | ^^^^^^^ - var_2   (mx, my) : y component of the vector            |
C   |                              (cartesian frame on the map)         |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE VecRot (grd_lon,grd_lat,dx,var_1,var_2)


      IMPLICIT NONE


C +---Dimensions
C +   ----------
 
      INCLUDE 'NSTdim.inc'

 
C +---Input
C +   -----
      REAL grd_lon(mx,my),grd_lat(mx,my)
      REAL dx

 
C +---Input and output
C +   ----------------

      REAL var_1(mx,my),var_2(mx,my)

 
C +---Local variables
C +   ---------------

      INTEGER i,j,ii,jj,iii,jjj,cmpt
      REAL m11,m12,m21,m22,vx,vy,rayter,DtR,dist_sp,dist_np,
     .     dist_dx,d_sp,d_np,aux_v1,aux_v2,dlamx,dphix,dlamy,
     .     dphiy,auxi1_lon,auxi2_lon,auxj1_lon,auxj2_lon,
     .     auxi1_lat,auxi2_lat,auxj1_lat,auxj2_lat


C +---Data
C +   ----

      DATA rayter  / 6371229.0 /
      DtR = ACOS(-1.)/180.


C +---Correction on wind direction (stereog. grid only)
C +   =================================================

      DO j = 1,my
       jj = j
       IF (j.eq.1 ) jj=2
       IF (j.eq.my) jj=my-1

       DO i = 1,mx
	ii = i
	IF (i.eq.1) ii = 2
	IF (i.eq.mx)ii = mx-1

        auxi1_lon = grd_lon(ii+1,jj  ) 
        auxi2_lon = grd_lon(ii-1,jj  ) 
        auxj1_lon = grd_lon(ii  ,jj+1)
        auxj2_lon = grd_lon(ii  ,jj-1)
        auxi1_lat = grd_lat(ii+1,jj  ) 
        auxi2_lat = grd_lat(ii-1,jj  ) 
        auxj1_lat = grd_lat(ii  ,jj+1)
        auxj2_lat = grd_lat(ii  ,jj-1)
        IF ((auxi1_lon-auxi2_lon).gt.180.) auxi2_lon=auxi2_lon+360.
        IF ((auxi2_lon-auxi1_lon).gt.180.) auxi1_lon=auxi1_lon+360.
        IF ((auxj1_lon-auxj2_lon).gt.180.) auxj2_lon=auxj2_lon+360.
        IF ((auxj2_lon-auxj1_lon).gt.180.) auxj1_lon=auxj1_lon+360.
        IF ((auxi1_lat-auxi2_lat).gt. 90.) auxi2_lat=auxi2_lat+180.
        IF ((auxi2_lat-auxi1_lat).gt. 90.) auxi1_lat=auxi1_lat+180.
        IF ((auxj1_lat-auxj2_lat).gt. 90.) auxj2_lat=auxj2_lat+180.
        IF ((auxj2_lat-auxj1_lat).gt. 90.) auxj1_lat=auxj1_lat+180.

C +---Correction for latitude and longitude
C +   -------------------------------------

        dlamx = (auxi1_lon-auxi2_lon)/(2*dx)
        dphix = (auxi1_lat-auxi2_lat)/(2*dx)
        dlamy = (auxj1_lon-auxj2_lon)/(2*dx)
        dphiy = (auxj1_lat-auxj2_lat)/(2*dx)

        m11 =  dphiy * DtR * rayter
        m12 = -dlamy * DtR * rayter * cos(grd_lat(ii,jj) * DtR)
        m21 = -dphix * DtR * rayter
        m22 =  dlamx * DtR * rayter * cos(grd_lat(ii,jj) * DtR)

C +...equivalent to:
C +...  m11 =  dlamx * DtR * rayter * cos(grd_lat(ii,jj) * DtR)
C +...  m12 =  dphix * DtR * rayter
C +...  m21 =  dlamy * DtR * rayter * cos(grd_lat(ii,jj) * DtR)
C +...  m22 =  dphiy * DtR * rayter
C +...
C +...or (simplier):
C +...  m11 =  dphiy * DtR * rayter
C +...  m12 =  dphix * DtR * rayter
C +...  m21 = -dphix * DtR * rayter
C +...  m22 =  dphiy * DtR * rayter

C +---Corrected wind direction
C +   ------------------------

        vx = m11 * var_1(i,j) + m12 * var_2(i,j)
        vy = m21 * var_1(i,j) + m22 * var_2(i,j)
        var_1(i,j) = vx
        var_2(i,j) = vy

       ENDDO
      ENDDO


C +---Special treatment close to North/South pole
C +   ===========================================
 
      DO j = 1,my
       jj = j
       IF (j.eq.1 ) jj=2
       IF (j.eq.my) jj=my-1

       DO i = 1,mx
	ii = i
	IF (i.eq.1) ii = 2
	IF (i.eq.mx)ii = mx-1

        dist_sp=ABS(grd_lat(i,j)+90.) * 111111.1
        dist_np=ABS(grd_lat(i,j)-90.) * 111111.1
        dist_dx=dx

        IF (dist_sp.lt.dist_dx .or. dist_np.lt.dist_dx) THEN

C +...For grid points whose distance with the pole is less than
C +...Dx, no satisfying wind direction can be computed using the
C +...above-mentioned formula. Therefore a "mean wind" is computed
C +...considering the 8 closest grid points around (i,j).

         cmpt  =0
         aux_v1=0.
         aux_v2=0.
         DO jjj=jj-1,jj+1
         DO iii=ii-1,ii+1
          d_sp=ABS(grd_lat(iii,jjj)+90.) * 111111.1
          d_np=ABS(grd_lat(iii,jjj)-90.) * 111111.1
          IF (d_sp.gt.dist_dx .and. d_np.gt.dist_dx) THEN
           aux_v1=aux_v1+var_1(iii,jjj)
           aux_v2=aux_v2+var_2(iii,jjj)
           cmpt  =cmpt+1
          ENDIF
         ENDDO
         ENDDO
         IF (cmpt.gt.0) THEN
          var_1(i,j)=aux_v1/REAL(cmpt)
          var_2(i,j)=aux_v2/REAL(cmpt)
         ENDIF
        ENDIF
 
       ENDDO
      ENDDO


      RETURN
      END


C   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

C   +---------------------------------------------------------------------+
C   | Subroutine VecRot_StereoSouth                     May 2009  C.Agosta |
C   +---------------------------------------------------------------------+
C   | Projection of the wind vectors in the polar south stereo. grid      |
C   | lon0 : lon-Direction (2D runs only ; 90 = East, clockwise)          |
C   | MAR : lon0 = GEddxx = 75                                            |
C   +---------------------------------------------------------------------+

      Subroutine VecRot_StereoSouth(lon0,NSTlon,INT_uu,INT_vv)

      Implicit None

      Include 'NSTdim.inc'
      Real, Intent(in)    :: lon0
      Real, Intent(in)    :: NSTlon(mx,my)
      Real, Intent(inout) :: INT_uu(mx,my)
      Real, Intent(inout) :: INT_vv(mx,my)

      Integer i,j
      Real pi, dr, phi, cphi, sphi, deltaphi
      Real uu, vv

      pi = 4.e0*atan(1.)
      dr = pi/180.

      deltaphi = 90. - lon0

      Do j=1,my
      Do i=1,mx
        phi = (-1.) * (NSTlon(i,j)+deltaphi) * dr
        cphi = cos( phi )
        sphi = sin( phi )
        uu = cphi*INT_uu(i,j) - sphi*INT_vv(i,j)
        vv = sphi*INT_uu(i,j) + cphi*INT_vv(i,j)
        INT_uu(i,j) = uu
        INT_vv(i,j) = vv
      EndDo
      EndDo

      Return
      End Subroutine VecRot_StereoSouth
      
C   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
