C +  +-------------------------------------------------------------------+
C +    File contents:           
C +       INThor
C +       INTbic
C +       INTmean
C +       INTlin
C +       SPLINE, SPLINT, SPLIE2, SPLIN2 (from Numerical Recipies)
C +  +-------------------------------------------------------------------+

C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine INThor      +                           12/95    MAR  +
C**  +-------------------------+                                         +
C**  +  * Horizontal interpolation from LS grid to MAR grid              +
C**  +    distribute tasks to bicubic, linear... routines                +
C**  +    according to the "intype" variable (1= bilin, 3= bicub)        +
C**  +    (! this routine uses the LS and MAR grid dimensions)           +
C**  +-------------------------------------------------------------------+
C**  +  INPUT : intype                : requested interpolation type     +
C**  +          grd_Ix (LSni)         : Input grid points position x(i)  +
C**  +          grd_Iy (LSnj)         :   "     "    "       "     y(j)  +
C**  +          var_I  (LSni, LSnj) : Input field values                 +
C**  +          grd_Ox (mx, my) : Output grid positions x(i,j)           +
C**  +          grd_Oy (mx, my) : Output grid positions y(i,j)           +
C**  +                                                                   +
C**  +                                                                   +
C**  +  OUTPUT: var_O  (mx, my) : Output field values                    +
C**  +                                                                   +
C**  +-------------------------------------------------------------------+
C
      subroutine INThor (intype, 
     &                   grd_Ix, grd_Iy, var_I,
     &                   grd_Ox, grd_Oy, var_O)
 

C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

      REAL wk1_LS (LSni, LSnj), wk2_LS(LSni, LSnj)

C +  Input
C +  ~~~~~
      INTEGER intype
      REAL grd_Ix (LSni),   grd_Iy (LSnj)
      REAL grd_Ox (mx, my), grd_Oy (mx, my)
      REAL var_I  (LSni, LSnj)

C +  Output
C +  ~~~~~~
      REAL var_O  (mx, my)

C +  Create temporary arrays:
C +  ~~~~~~~~~~~~~~~~~~~~~~~~
      REAL tmp_I2a (LSni, LSnj)
      REAL tmp_in  (LSni, LSnj)
      REAL samOx  (mx, my, 5)
      REAL samOy  (mx, my, 5)

C +  Interpolation:
C +  ==============

      IF      (intype.EQ.-1) THEN
C +   * Bilinear interpolation.

        CALL INTblin (tmp_in,
     &       LSni, LSnj, grd_Ix, grd_Iy, var_I,
     &       mx, my, grd_Ox, grd_Oy, var_O)

      ELSE IF (intype.EQ.1) THEN
C +   * The usual "mean-on-the-mesh" bilinear interpolation.

        CALL INTmean (tmp_in, samOx, samOy,
     &       LSni, LSnj, grd_Ix, grd_Iy, var_I,
     &       mx, my, grd_Ox, grd_Oy, var_O)

      ELSE IF (intype.EQ.3) THEN 
C +   * Bicubic interpolation.

        CALL INTbic (tmp_I2a, tmp_in,
     &       LSni, LSnj, grd_Ix, grd_Iy, var_I,
     &       mx, my, grd_Ox, grd_Oy, var_O)

      ENDIF


      RETURN
      END


C**  +-------------------------+-----------------------------------------+
C**  +  Subroutine INTbic      +                           12/95    MAR  +
C**  +-------------------------+                                         +
C**  +  * Natural bicubic spline interpolation of a 2D scalar field,     +
C**  +    from a  RECTANGULAR grid   to an  IRREGULAR grid.              +
C**  +    i.e.    var_I(i,j) (coarse) ===>  var_O(ii,jj) (fine mesh)     +
C**  +                                                                   +
C**  +    (Note 1: The independant variable should be adequatelly        +
C**  +        choosed to have the "regularity" in the INPUT grid)        +
C**  +     Note 2: This routine is simple but not efficient;             +
C**  +        a faster routine is available for 2 rectangular grids,     +
C**  +        and modifications to INTbic possible if necessary.         +
C**  +                                                                   +
C**  +  Authors    : Marbaix Ph. (this version), Wang Y (general study)  +
C**  +  Revisions  :                                                     +
C**  +                                                                   +
C**  +  References : W.H.Pres et al., Numerical recipes,Cambridge UP.1992+
C**  +               Lancaster P., Curve and Surface fitting, AP. 1986   +
C**  +                                                                   +
C**  +  INPUT : grd_Ix (dim_Ix)         : Input grid points position x(i)+
C**  +          grd_Iy (dim_Iy)         :   "     "    "       "     y(j)+
C**  +          var_I  (dim_Ix, dim_Iy) : Input field values             +
C**  +          grd_Ox (dim_Ox, dim_Oy) : Output grid positions x(i,j)   +
C**  +          grd_Oy (dim_Ox, dim_Oy) : Output grid positions y(i,j)   +
C**  +                                                                   +
C**  +          dim_Ix, dim_Iy : (parameter type) array dimensions       +
C**  +                           ! = size of the interpolated data field +
C**  +          dim_Ox, dim_Oy : (parameter type) array dimensions       +
C**  +                           = size of the used data field           +
C**  +                                                                   +
C**  +  OUTPUT: var_O  (dim_Ox, dim_Oy) : Output field values            +
C**  +                                                                   +
C**  +  TEMPORARY arrays : tmp_I2a (dim_Ix, dim_Iy)                      +
C**  +                 and tmp_in  (dim_Ix, dim_Iy) must be provided     +
C**  +-------------------------------------------------------------------+
C
      subroutine INTbic (tmp_I2a, tmp_in,
     &       dim_Ix, dim_Iy, grd_Ix, grd_Iy, var_I,
     &       dim_Ox, dim_Oy, grd_Ox, grd_Oy, var_O)
 
C     ** input
      INTEGER dim_Ix, dim_Iy
      INTEGER dim_Ox, dim_Oy
      REAL grd_Ix (dim_Ix),         grd_Iy (dim_Iy)
      REAL grd_Ox (dim_Ox, dim_Oy), grd_Oy (dim_Ox, dim_Oy)
      REAL var_I  (dim_Ix, dim_Iy)
 
C     ** output
      REAL var_O  (dim_Ox, dim_Oy)
 
C     ** temporary arrays:
      REAL tmp_I2a (dim_Ix, dim_Iy)
      REAL tmp_in  (dim_Ix, dim_Iy)
C     ** local:
      INTEGER i,j,ii,jj
      REAL tmp_Ix(500)
      REAL tmp_Iy(500)
 
      icheck = 0  ! Debuging output level

      IF (icheck.ge.1) WRITE(*,*) 'INTbic: Begin'
 
C +-- 1. Coordinates indexes inversion (if necessary).
C     ------------------------------------------------
 
C     ** Note : unmodified routines from Num. Recipes
C     **        implies that coords are such as
C     **        x(1) < x(2) <...  ===> revert if necessary
C     ** tmp is used to keep input data unmodified.
 
      DO jj = 1,dim_Iy
        DO ii = 1,dim_Ix  !(vectorizable loop)
          tmp_in (ii,jj) = var_I(ii, jj)
        ENDDO
      ENDDO
      DO ii = 1,dim_Ix
        tmp_Ix (ii) = grd_Ix (ii)
      ENDDO
      DO jj = 1,dim_Iy
        tmp_Iy (jj) = grd_Iy (jj)
      ENDDO
 
C     ** Revert grd_Ix (1) <--> grd_Ix (n), ... ?
 
      IF (grd_Ix(dim_Ix).lt.grd_Ix(1)) THEN
        DO ii = 1,dim_Ix
          DO jj = 1,dim_Iy
           tmp_in (ii,jj) = var_I(dim_Ix-ii+1, jj)
          ENDDO
          tmp_Ix (ii) = grd_Ix (dim_Ix-ii+1)
        ENDDO
      IF(icheck.ge.2)WRITE(*,*)'Lon. coord. indexes reverted'
      ENDIF
 
C     ** Revert grd_Iy (1) <--> grd_Iy (n), ... ?
 
      IF (grd_Iy(dim_Iy).lt.grd_Iy(1)) THEN
        DO jj = 1,dim_Iy
          DO ii = 1,dim_Ix
           tmp_in (ii,jj) = var_I(ii,dim_Iy-jj+1)
          ENDDO
          tmp_Iy (jj) = grd_Iy (dim_Iy-jj+1)
        ENDDO
      IF(icheck.ge.2)WRITE(*,*)'Lat. coord. indexes reverted'
      ENDIF
 
C +-- 2. Construction of 1D splines of rows
C     -------------------------------------
 
      CALL SPLIE2(tmp_Ix, tmp_Iy, tmp_in, dim_Ix, dim_Iy, tmp_I2a)
 
C +-- 3. Interpolation for each output grid point.
C     --------------------------------------------
 
      DO j = 1,dim_Oy
       DO i = 1,dim_Ox
 
        CALL SPLIN2(tmp_Ix, tmp_Iy, tmp_in, tmp_I2a, dim_Ix, dim_Iy,
     &              grd_Ox (i,j), grd_Oy (i,j), var_O(i,j) )
 
       ENDDO
      ENDDO
 
 
      IF (icheck.ge.2) WRITE(*,*) 'INTbic: End  '
 
      RETURN
      END

C   +-------------------------------------------------------------------+
C   | Subroutine INTmean                               1996.09  NESTOR |
C   +-------------------------------------------------------------------+
C   |  * modified BILINEAR interpolation of a 2D scalar field:          |
C   |    AVERAGE of 5 BILINEAR interpolations at 5 sampling points      |
C   |    located around the selected point in the output mesh.          |
C   |    (usefull when output resolution is LOWER then input)           |
C   |                                                                   |
C   |  Warning: this version was modified for CYCLIC fields!            |
C   |  ^^^^^^^^ = x bondary is assume cyclic when crossed,              |
C   |           but works only for x = longitude in degrees !           |
C   |                              ^^^^^^^^^^^^^^^^^^^^^^^^^^           |
C   |                                                                   |
C   |  ---Ph. Marbaix                                                   |
C   |                                                                   |
C   |  INPUT : grd_Ix (dim_Ix)         : Input grid points position x(i)|
C   |  ^^^^^^^ grd_Iy (dim_Iy)         :   "     "    "       "     y(j)|
C   |          var_I  (dim_Ix, dim_Iy) : Input field values             |
C   |          grd_Ox (dim_Ox, dim_Oy) : Output grid positions x(i,j)   |
C   |          grd_Oy (dim_Ox, dim_Oy) : Output grid positions y(i,j)   |
C   |                                                                   |
C   |          dim_Ix, dim_Iy : (parameter type) array dimensions       |
C   |                           ! = size of the interpolated data field |
C   |          dim_Ox, dim_Oy : (parameter type) array dimensions       |
C   |                           = size of the used data field           |
C   |                                                                   |
C   |  OUTPUT: var_O  (dim_Ox, dim_Oy) : Output field values            |
C   |  ^^^^^^^                                                          |
C   |  TEMPORARY arrays : tmp_in (dim_Ix, dim_Iy)                       |
C   |  ^^^^^^^^^^^^^^^^^^ samOx  (dim_Ox, dim_Oy, 5)                    |
C   |                 and samOy  (dim_Ox, dim_Oy, 5) must be provided   |
C   +-------------------------------------------------------------------+
      SUBROUTINE INTmean   (tmp_in, samOx, samOy,       
     &       dim_Ix, dim_Iy, grd_Ix, grd_Iy, var_I, 
     &       dim_Ox, dim_Oy, grd_Ox, grd_Oy, var_O)

      IMPLICIT NONE

C     ** number of "sampling points" :
      INTEGER ns
      PARAMETER (ns = 5)

C     ** input
      INTEGER dim_Ix, dim_Iy
      INTEGER dim_Ox, dim_Oy
      REAL grd_Ix (dim_Ix),         grd_Iy (dim_Iy)
      REAL grd_Ox (dim_Ox, dim_Oy), grd_Oy (dim_Ox, dim_Oy)
      REAL var_I  (dim_Ix, dim_Iy)

C     ** output
      REAL var_O  (dim_Ox, dim_Oy)

C     ** temporary arrays:
      REAL tmp_in  (dim_Ix, dim_Iy)
      REAL samOx  (dim_Ox, dim_Oy, 5)
      REAL samOy  (dim_Ox, dim_Oy, 5)
C     ** local:
      INTEGER   LocDim
      PARAMETER(LocDim=1000) ! dim of Local 1D arrays
      INTEGER i,j,ii,jj,p,q,is
      REAL x,y,tmp,tmp2,x0,x1,y0,y1
      REAL tmp_Ix(0:LocDim+1)
      REAL tmp_Iy(0:LocDim+1)
      INTEGER icc(0:LocDim+1), jcc(-1:LocDim+1)
      LOGICAL cyclic
      
      INTEGER icheck
      icheck = 0     ! Debbuging output level
       
      cyclic = .FALSE. 
C +.. This routine sets cyclic to T automatically
C +.. If data boundary is crossed. 

      IF (icheck.ge.1) WRITE(*,*) 'INTmean : Begin'
      IF (dim_Ix.gt.LocDim .or. dim_Iy.gt.LocDim) THEN
        WRITE(*,*) 'INTmean - fatal error: dimension',LocDim
        WRITE(*,*) ' (change LocDim or correct CALL error)'
        STOP
      ENDIF


C +---Coordinates indexes inversion (if necessary).
C +   ---------------------------------------------
C     (tmp is used to keep input data unmodified.)

      DO jj = 1,dim_Iy
        DO ii = 1,dim_Ix 
          tmp_in (ii,jj) = var_I(ii, jj)
        ENDDO
      ENDDO
      DO ii = 1,dim_Ix
        tmp_Ix (ii) = grd_Ix (ii)
      ENDDO
      DO jj = 1,dim_Iy
        tmp_Iy (jj) = grd_Iy (jj)
      ENDDO

C     ** Revert grd_Ix (1) <--> grd_Ix (n), ... ?
      
      IF (grd_Ix(dim_Ix).lt.grd_Ix(1)) THEN     
        DO ii = 1,dim_Ix   
          DO jj = 1,dim_Iy                       
           tmp_in (ii,jj) = var_I(dim_Ix-ii+1, jj)
          ENDDO
          tmp_Ix (ii) = grd_Ix (dim_Ix-ii+1) 
        ENDDO
      IF(icheck.ge.2)WRITE(*,*)'Lon. coord. indexes reverted'
      ENDIF

C     ** Revert grd_Iy (1) <--> grd_Iy (n), ... ?
      
      IF (grd_Iy(dim_Iy).lt.grd_Iy(1)) THEN     
        DO jj = 1,dim_Iy   
          DO ii = 1,dim_Ix                       
           tmp_in (ii,jj) = var_I(ii,dim_Iy-jj+1)
          ENDDO
          tmp_Iy (jj) = grd_Iy (dim_Iy-jj+1)
        ENDDO
      IF(icheck.ge.2)WRITE(*,*)'Lat. coord. indexes reverted'
      ENDIF

C +---Define index conversion if cyclic field:
C     ----------------------------------------
      DO ii = 1,dim_Ix
        icc(ii)= ii
      ENDDO
      DO jj = 1,dim_Iy
        jcc(jj)= jj
      ENDDO
      icc(0)        = dim_Ix
      jcc(0)        = 1
      icc(dim_Ix+1) = 1
      jcc(dim_Iy+1) = dim_Iy
C +..Longitude is cyclic but latitude is "limited"

      tmp_Ix(0)        = 2.*tmp_Ix(1)-tmp_Ix(2)
      tmp_Ix(dim_Ix+1) = 2.*tmp_Ix(dim_Ix)-tmp_Ix(dim_Ix-1) 
C +..Define the "cross boundary" longitude.

      tmp_Iy(0)        = 2.*tmp_Iy(1)-tmp_Iy(2)
      tmp_Iy(dim_Iy+1) = 2.*tmp_Iy(dim_Iy)-tmp_Iy(dim_Iy-1)
C +..Define the "cross boundary" latitude. 

C +---Define the sampling points positions:
C     -------------------------------------

      DO j=2,dim_Oy-1
       DO i=2,dim_Ox-1
         samOx(i,j,1) = grd_Ox(i,j)
         samOy(i,j,1) = grd_Oy(i,j)
         samOx(i,j,2) = 0.6 * grd_Ox(i,j) + 0.4 * grd_Ox(i-1,j+1)
         samOy(i,j,2) = 0.6 * grd_Oy(i,j) + 0.4 * grd_Oy(i-1,j+1)
         samOx(i,j,3) = 0.6 * grd_Ox(i,j) + 0.4 * grd_Ox(i+1,j-1)
         samOy(i,j,3) = 0.6 * grd_Oy(i,j) + 0.4 * grd_Oy(i+1,j-1)
         samOx(i,j,4) = 0.6 * grd_Ox(i,j) + 0.4 * grd_Ox(i-1,j-1)
         samOy(i,j,4) = 0.6 * grd_Oy(i,j) + 0.4 * grd_Oy(i-1,j-1)
         samOx(i,j,5) = 0.6 * grd_Ox(i,j) + 0.4 * grd_Ox(i+1,j+1)
         samOy(i,j,5) = 0.6 * grd_Oy(i,j) + 0.4 * grd_Oy(i+1,j+1)
       END DO
      END DO
C     ** Exterior (edge row)  grid points : no sampling
C     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      DO j=1,dim_Oy
       DO is = 1,ns
         samOx(1,j,is) = grd_Ox(1,j)
         samOy(1,j,is) = grd_Oy(1,j)
         samOx(dim_Ox,j,is) = grd_Ox(dim_Ox,j)
         samOy(dim_Ox,j,is) = grd_Oy(dim_Ox,j)
       END DO
      END DO
      DO i=1,dim_Ox
       DO is = 1,ns
         samOx(i,1,is) = grd_Ox(i,1)
         samOy(i,1,is) = grd_Oy(i,1)
         samOx(i,dim_Oy,is) = grd_Ox(i,dim_Oy)
         samOy(i,dim_Oy,is) = grd_Oy(i,dim_Oy)
       END DO
      END DO


C +---Bi-linear interpolation.
C     ------------------------
C
C     ** Initial values for searching input mesh
C     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      p = 1
      q = 1
      
      DO i=1,dim_Ox    ! LOOP on output grid-points : BEGIN
       DO j=1,dim_Oy  
          
C      ** initialise tmp2 (to sum sampled values):
       tmp2 = 0.0

       DO is = 1, ns  ! Loop on the sampling points: BEGIN
        x = samOx(i,j,is)
        y = samOy(i,j,is)

C       ** search for the bottom-left corner of the surrounding mesh :
C       ** (This simple method accounts for the fact that two
C       **  successive requests usually refer to neighbour points)
C       --------------------------------------------------------------

C +-----search for dimension 1 value
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF (x.ge.180) x = x - 360.
C +     Large Scale data is supposed to use -180 to +180 longitudes...

        IF (tmp_Ix(p).le.x) THEN

C         ** Search upwards
          DO WHILE (tmp_Ix(p+1).lt.x)
            p = p + 1
          END DO
        ELSE
C         ** Search downwards
          DO WHILE (tmp_Ix(p).gt.x)
            p = p - 1
          END DO
        END IF

C +-----search for dimension 2 value
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        IF (tmp_Iy(q).le.y ) THEN
C         ** Search upwards
          DO WHILE (tmp_Iy(q+1).lt.y)
            q = q + 1
          END DO
        ELSE
C         ** Search downwards
          DO WHILE (tmp_Iy(q).gt.y)
            q = q - 1
          END DO
        END IF
        IF (icheck.ge.3) WRITE(*,*) tmp_Ix(p),tmp_Iy(q),x,y

C +-----interpolation
C +     ~~~~~~~~~~~~~

        x0 = tmp_Ix(icc(p))
        x1 = tmp_Ix(icc(p+1))
        y0 = tmp_Iy(jcc(q))
        y1 = tmp_Iy(jcc(q+1))
        tmp=(x-x0)*((y-y0)*tmp_in(icc(p+1),jcc(q+1))
     &             +(y1-y)*tmp_in(icc(p+1),jcc(q  )))
     &     +(x1-x)*((y-y0)*tmp_in(icc(p  ),jcc(q+1))
     &             +(y1-y)*tmp_in(icc(p  ),jcc(q  )))
        tmp2 = tmp2 + tmp / ( (x1-x0)*(y1-y0) )
        IF (icheck.ge.3) WRITE(*,*) tmp , x, y, p, q

        IF (p.lt.1 .OR. q.lt.1 
     &    .OR. p.gt.dim_Ix .OR. q.gt.dim_Iy ) THEN 
          cyclic = .TRUE.
        ENDIF

       END DO         ! LOOP on the sampling points: END

C +----output value = average of the samplings:
C +    ----------------------------------------
       var_O(i,j) = tmp2 / REAL(ns)

      END DO
      END DO          ! Loop on output grid-points : END

      IF (cyclic)WRITE(*,*) 'INTmean-info: cyclic boundaries'

      IF (icheck.ge.2) WRITE(*,*) 'INTmean : End'
      END
C
C
C
C
C     +--------------------------------------------------------------+
      subroutine INTlin (xx,vv,ni,xreq,outvar)
C     +--------------------------------------------------------------+

      REAL xx(ni), vv(ni)
      REAL xreq, outvar, fdis
      INTEGER ind, KLO, KHI

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


      RETURN
      END

C     +--------------------------------------------------------------+
      subroutine SPLIE2(X1A,X2A,YA,M,N,Y2A)
C     +     * From numerical recipes (H. Press et al., 1992)         +
C     +--------------------------------------------------------------+
 
      PARAMETER (NN=500)
C     ** NN = max value allowed for N (for 1D arrays only -> overdim.)
 
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N)
      DIMENSION YTMP(NN),Y2TMP(NN)
 
      DO 13 J=1,M
 
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
 
        CALL SPLINE(X2A,YTMP,N,1.E30,1.E30,Y2TMP)
C       ** NB : 1.E30 = switch value to select "natural" bicub spline
 
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
 
13    CONTINUE
 
      RETURN
      END
 
C     +--------------------------------------------------------------+
      subroutine SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
C     +     * From numerical recipes (H. Press et al., 1992)         +
C     +       Version grossierement optimisee                        +
C     +--------------------------------------------------------------+

      PARAMETER (NN=500)
C +   NN = max value allowed for N (for 1D arrays only -> overdim.)

      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),Y2TMP(NN)
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
C     +     * From numerical recipes (H. Press et al., 1992)         +
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
C     +     * From numerical recipes (H. Press et al., 1992)         +
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

C   +-------------------------------------------------------------------+
C   | Subroutine INTblin                               1999.03  LSM*    |
C   +-------------------------------------------------------------------+
C   |  * BILINEAR interpolation of a 2D scalar field:                   |
C   |    Indentical to INTmean but the sampling points are removed.     |
C   |                                                                   |
C   |  Warning: this version was modified for CYCLIC fields!            |
C   |  ^^^^^^^^ = x bondary is assume cyclic when crossed,              |
C   |           but works only for x = longitude in degrees !           |
C   |                              ^^^^^^^^^^^^^^^^^^^^^^^^^^           |
C   |                                                                   |
C   |  INPUT : grd_Ix (dim_Ix)         : Input grid points position x(i)|
C   |  ^^^^^^^ grd_Iy (dim_Iy)         :   "     "    "       "     y(j)|
C   |          var_I  (dim_Ix, dim_Iy) : Input field values             |
C   |          grd_Ox (dim_Ox, dim_Oy) : Output grid positions x(i,j)   |
C   |          grd_Oy (dim_Ox, dim_Oy) : Output grid positions y(i,j)   |
C   |                                                                   |
C   |          dim_Ix, dim_Iy : (parameter type) array dimensions       |
C   |                           ! = size of the interpolated data field |
C   |          dim_Ox, dim_Oy : (parameter type) array dimensions       |
C   |                           = size of the used data field           |
C   |                                                                   |
C   |  OUTPUT: var_O  (dim_Ox, dim_Oy) : Output field values            |
C   |  ^^^^^^^                                                          |
C   |  TEMPORARY arrays : tmp_in (dim_Ix, dim_Iy)                       |
C   +-------------------------------------------------------------------+
      SUBROUTINE INTblin   (tmp_in, 
     &       dim_Ix, dim_Iy, grd_Ix, grd_Iy, var_I, 
     &       dim_Ox, dim_Oy, grd_Ox, grd_Oy, var_O)

      IMPLICIT NONE

C     ** input
      INTEGER dim_Ix, dim_Iy
      INTEGER dim_Ox, dim_Oy
      REAL grd_Ix (dim_Ix),         grd_Iy (dim_Iy)
      REAL grd_Ox (dim_Ox, dim_Oy), grd_Oy (dim_Ox, dim_Oy)
      REAL var_I  (dim_Ix, dim_Iy)

C     ** output
      REAL var_O  (dim_Ox, dim_Oy)

C     ** temporary arrays:
      REAL tmp_in  (dim_Ix, dim_Iy)
C     ** local:
      INTEGER   LocDim
      PARAMETER(LocDim=4500) ! dim of Local 1D arrays
      INTEGER i,j,ii,jj,p,q
      REAL x,y,tmp,x0,x1,y0,y1
      REAL tmp_Ix(0:LocDim+1)
      REAL tmp_Iy(0:LocDim+1)
      INTEGER icc(0:LocDim+1), jcc(-1:LocDim+1)
      LOGICAL cyclic, loverflow,laterr
      REAL dmax, dmin  

      INTEGER icheck
      icheck = 0     ! Debbuging output level
       
      laterr = .FALSE.
      cyclic = .FALSE. 
C +.. This routine sets cyclic to T automatically
C +.. If data boundary is crossed. 

      IF (icheck.ge.1) WRITE(*,*) 'INTblin : Begin'
      IF (dim_Ix.gt.LocDim .or. dim_Iy.gt.LocDim) THEN
        WRITE(*,*) 'INTblin - fatal error: dimension',LocDim
        WRITE(*,*) ' (change LocDim or correct CALL error)'
        STOP
      ENDIF

C +---Coordinates indexes inversion (if necessary).
C +   ---------------------------------------------
C     (tmp is used to keep input data unmodified.)

      loverflow = .FALSE.
      DO jj = 1,dim_Iy
        DO ii = 1,dim_Ix 
          IF (var_I(ii,jj).LE.1E30) THEN 
            tmp_in (ii,jj) = var_I(ii, jj)
          ELSE
            tmp_in (ii,jj) = 1.E20
            loverflow = .TRUE.
          ENDIF
        ENDDO
      ENDDO
      IF (loverflow) write(*,*) 'INTblin: input overflow'
      DO ii = 1,dim_Ix
        tmp_Ix (ii) = grd_Ix (ii)
      ENDDO
      DO jj = 1,dim_Iy
        tmp_Iy (jj) = grd_Iy (jj)
      ENDDO

C     ** Revert grd_Ix (1) <--> grd_Ix (n), ... ?
      
      IF (grd_Ix(dim_Ix).lt.grd_Ix(1)) THEN     
        DO ii = 1,dim_Ix   
          DO jj = 1,dim_Iy                       
           tmp_in (ii,jj) = var_I(dim_Ix-ii+1, jj)
          ENDDO
          tmp_Ix (ii) = grd_Ix (dim_Ix-ii+1) 
        ENDDO
      IF(icheck.ge.2)WRITE(*,*)'Lon. coord. indexes reverted'
      ENDIF

C     ** Revert grd_Iy (1) <--> grd_Iy (n), ... ?
      
      IF (grd_Iy(dim_Iy).lt.grd_Iy(1)) THEN     
        DO jj = 1,dim_Iy   
          DO ii = 1,dim_Ix                       
           tmp_in (ii,jj) = var_I(ii,dim_Iy-jj+1)
          ENDDO
          tmp_Iy (jj) = grd_Iy (dim_Iy-jj+1)
        ENDDO
      IF(icheck.ge.2)WRITE(*,*)'Lat. coord. indexes reverted'
      ENDIF

C +---Define index conversion if cyclic field:
C     ----------------------------------------
      DO ii = 1,dim_Ix
        icc(ii)= ii
      ENDDO
      DO jj = 1,dim_Iy
        jcc(jj)= jj
      ENDDO
      icc(0)        = dim_Ix
      jcc(0)        = 1
      icc(dim_Ix+1) = 1
      jcc(dim_Iy+1) = dim_Iy
C +..Longitude is cyclic but latitude is "limited"

      tmp_Ix(0)        = 2.*tmp_Ix(1)-tmp_Ix(2)
      tmp_Ix(dim_Ix+1) = 2.*tmp_Ix(dim_Ix)-tmp_Ix(dim_Ix-1) 
C +..Define the "cross boundary" longitude.

      tmp_Iy(0)        = 2.*tmp_Iy(1)-tmp_Iy(2)
      tmp_Iy(dim_Iy+1) = 2.*tmp_Iy(dim_Iy)-tmp_Iy(dim_Iy-1)
C +..Define the "cross boundary" latitude. 

C +---Bi-linear interpolation.
C     ------------------------
C
C     ** Initial values for searching input mesh
C     ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      p = 1
      q = 1
      
      DO i=1,dim_Ox    ! LOOP on output grid-points : BEGIN
       DO j=1,dim_Oy  
          
        x = grd_Ox(i,j)
        y = grd_Oy(i,j)

C       ** search for the bottom-left corner of the surrounding mesh :
C       ** (This simple method accounts for the fact that two
C       **  successive requests usually refer to neighbour points)
C       --------------------------------------------------------------

C +-----search for dimension 1 value
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF (x.ge.180) x = x - 360.
C +     Large Scale data is supposed to use -180 to +180 longitudes...

        IF (tmp_Ix(p).le.x) THEN

C         ** Search upwards
          DO WHILE (tmp_Ix(p+1).lt.x)
            p = p + 1
          END DO
        ELSE
C         ** Search downwards
          DO WHILE (tmp_Ix(p).gt.x)
            p = p - 1
          END DO
        END IF

C +-----search for dimension 2 value
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        IF (tmp_Iy(q).le.y ) THEN
C         ** Search upwards
          DO WHILE (tmp_Iy(q+1).lt.y)
            q = q + 1
          END DO
        ELSE
C         ** Search downwards
          DO WHILE (tmp_Iy(q).gt.y)
            q = q - 1
          END DO
        END IF
        IF (icheck.ge.3) WRITE(*,*) tmp_Ix(p),tmp_Iy(q),x,y

C +-----interpolation
C +     ~~~~~~~~~~~~~

        x0 = tmp_Ix(icc(p))
        x1 = tmp_Ix(icc(p+1))
        y0 = tmp_Iy(jcc(q))
        y1 = tmp_Iy(jcc(q+1))
        tmp=(x-x0)*((y-y0)*tmp_in(icc(p+1),jcc(q+1))
     &             +(y1-y)*tmp_in(icc(p+1),jcc(q  )))
     &     +(x1-x)*((y-y0)*tmp_in(icc(p  ),jcc(q+1))
     &             +(y1-y)*tmp_in(icc(p  ),jcc(q  )))

C +----output value 
C +     - - - - - -
        var_O(i,j) = tmp / ( (x1-x0)*(y1-y0) )

        IF (icheck.ge.3) WRITE(*,*) tmp , x, y, p, q

        IF (p.lt.1 .OR. p.gt.dim_Ix ) THEN 
          cyclic = .TRUE.
        ENDIF
        IF (q.lt.1 .OR. q.gt.dim_Ix ) THEN                      
          laterr = .TRUE.
        ENDIF


      END DO
      END DO          ! Loop on output grid-points : END


      IF (cyclic .OR. laterr) THEN

      IF (cyclic) WRITE(*,*) 'INTblin-info: cyclic boundaries'
      IF (laterr) WRITE(*,*) 'Latitude - error'

      WRITE(*,*) 'Input grid (LS)'

      dmin= grd_Iy(1)
      dmax=dmin
      DO jj = 1,dim_Iy
         dmin = min(grd_Iy(jj),dmin)
         dmax = max(grd_Iy(jj),dmax)
      ENDDO
      WRITE(*,*) 'latitude: ',dmin,dmax

      dmin= grd_Ix(1)
      dmax=dmin
      DO ii = 1,dim_Ix
         dmin = min(grd_Ix(ii),dmin)
         dmax = max(grd_Ix(ii),dmax)
      ENDDO
      WRITE(*,*) 'Longitude: ',dmin,dmax

      WRITE(*,*) 'Output grid (NST,MAR)'

      dmin= grd_Oy(1,1)
      dmax=dmin
      DO jj = 1,dim_Oy
      DO ii = 1,dim_Ox
         dmin = min(grd_Oy(ii,jj),dmin)
         dmax = max(grd_Oy(ii,jj),dmax)
      ENDDO
      ENDDO
      WRITE(*,*) 'Latitude: ',dmin,dmax

      dmin= grd_Ox(1,1)
      dmax=dmin
      DO jj = 1,dim_Oy
      DO ii = 1,dim_Ox
         dmin = min(grd_Ox(ii,jj),dmin)
         dmax = max(grd_Ox(ii,jj),dmax)
      ENDDO
      ENDDO
      WRITE(*,*) 'Longitude: ',dmin,dmax

      ENDIF


      IF (icheck.ge.2) WRITE(*,*) 'INTblin : End'
      END
