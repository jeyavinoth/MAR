C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine VecRotP     +                           01/96    MAR  +
C**  +-------------------------+                                         +
C**  +                                                                   +
C**  +  * Rotates vector (wind) following the rotation of the            +
C**  +    reference frame that occurs on a map projection.               +
C**  +                                                                   +
C**  +  Author     : Marbaix Ph.                                         +
C**  +  Revisions  :                                                     +
C**  +                                                                   +
C**  +  References : Map projections and equations of motion suitable    +
C**  +               for mesoscale alpha simulations by the MAR model.   +
C**  +               (Progress report, ASTR 1996)                        +
C**  +                                                                   +
C**  +  INPUT :                                                          +
C**  +          grd_lon (mx, my) : grid positions lon(i,j)               +
C**  +          grd_lat (mx, my) : grid positions lat(i,j)               +
C**  +          dx : mesh size; the grid is assumed to be square         +
C**  +                    and fully regullar                             +
C**  +          var_1  (mx, my): "x" component of the vector,            +
C**  +          var_2  (mx, my): "y" component of the vector,            +
C**  +                                  >> local cartesian on the sphere +
C**  +                                                                   +
C**  +  OUTPUT: var_1  (mx, my): x component of the vector               +
C**  +          var_2  (mx, my): y component of the vector,              +
C**  +                                 >> cartesian frame on the map.    +
C**  +-------------------------------------------------------------------+
C
        subroutine VecRotP (
     &     grd_lon, grd_lat, dx,
     &     var_1, var_2 )
 
C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'      
 
      INTEGER icheck
 
C     ** input
      REAL grd_lon (mx, my), grd_lat (mx, my)
      REAL dx
 
C     ** input AND output (!)
      REAL var_1  (mx, my), var_2  (mx, my)
 
C     ** local
      REAL dlamx, dphix, dlamy, dphiy
      REAL m11, m12, m21, m22
      REAL vx, vy
      REAL rayter
      REAL DtR
      INTEGER i,j,ii,jj
      PARAMETER (rayter = 6371229.0)
      DtR = ACOS(-1.)/180.

      icheck = 0           !Debugging output level
 
 
      IF (icheck.ge.1) WRITE(*,*) 'VecRotP  : Begin'
 
C     **
      DO j = 1,my
       jj = j
       IF (j.eq.1 ) jj=2
       IF (j.eq.my) jj=my-1
       DO i = 1,mx
	ii = i
	IF (i.eq.1) ii = 2
	IF (i.eq.mx)ii = mx-1
 
C       **
        dlamx = (grd_lon(ii+1,jj  )-grd_lon(ii-1,jj  ))/(2*dx)
        dphix = (grd_lat(ii+1,jj  )-grd_lat(ii-1,jj  ))/(2*dx)
        dlamy = (grd_lon(ii  ,jj+1)-grd_lon(ii  ,jj-1))/(2*dx)
        dphiy = (grd_lat(ii  ,jj+1)-grd_lat(ii  ,jj-1))/(2*dx)
C       **
        m11 =  dphiy*DtR * rayter
        m12 = -dlamy*DtR * rayter * cos(grd_lat(ii,jj) * DtR)
        m21 = -dphix*DtR * rayter
        m22 =  dlamx*DtR * rayter * cos(grd_lat(ii,jj) * DtR)
        IF (icheck.ge.2) WRITE(*,*)
     &    'Determinant: ', m11*m22-m12*m21
C       **
        vx = m11 * var_1(i,j) + m12 * var_2(i,j)
        vy = m21 * var_1(i,j) + m22 * var_2(i,j)
        var_1(i,j) = vx
        var_2(i,j) = vy
 
       ENDDO
      ENDDO
 
      IF (icheck.ge.1) WRITE(*,*) 'VecRotP  : End  '
      RETURN
      END
