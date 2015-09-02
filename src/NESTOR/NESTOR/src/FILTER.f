      subroutine DYNfil_2H (f,eps)
C +
C +------------------------------------------------------------------------+
C | MAR DYNAMICS FILTER (2-D)                              30-12-2000  MAR |
C |   SubRoutine DYNfil_2H  is used to Filter Horizontal Fields in 3D Code |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT:   f(i,j) : variable to be filtered for a particular Level k   |
C |   ^^^^^    eps    : value of the selectivity parameter                 |
C |                                                                        |
C |   OUTPUT:  f(i,j)                                                      |
C |   ^^^^^^                                                               |
C |                                                                        |
C |   LATERAL BOUNDARIES:                                                  |
C |   ^^^^^^^^^^^^^^^^^^^                                                  |
C |   1. The value    of the variable is fixed at the Boundary             |
C |                                                                        |
C |   REFER.:  Raymond and Garder, MWR 116, Jan 1988, p209                 |
C |   ^^^^^^                                                               |
C +------------------------------------------------------------------------+
C +
C +
      IMPLICIT NONE
C +
C +
C +--Global Variables
C +  ================
C +
      include 'NSTdim.inc'
C +
      real    f(mx,my),eps
C +
C +
C +--Local  Variables
C +  ================
C +
      integer i,j,ip11,jp11,mx1,my1,m,m1,m2,mn1,mn2,mx2,my2,mn
      real    aa,bb
C +
      real    x(mx,my)
      real    a1(1:mx) ,b1(1:mx) ,d1(1:mx) ,p1(mx) ,q1(mx) ,wk1(1:mx) 
      real    a2(1:my) ,b2(1:my) ,d2(1:my) ,p2(my) ,q2(my) ,wk2(1:my) 
C +
C +
C +--Initialisation
C +  ==============
C +
      ip11 = 2
      jp11 = 2
      mx1  = mx-1
      my1  = my-1
      mx2  = mx-2
      my2  = my-2
C +
      m  = mx
      m1 = mx1
      m2 = mx2
      mn = my
      mn1= my1
      mn2= my2
C +
C +
C +--1st Matrix Initialisation
C +  -------------------------
C +
        a1( 1) = 0.d0
        b1( 1) = 1.d0
        a1(mx) = 0.d0
        b1(mx) = 1.d0
C +
        aa     =       1.d0-eps
        bb     = 2.d0*(1.d0+eps)
C +
      DO i=ip11,mx1
        a1(i)   = aa
        b1(i)   = bb
      END DO
C +
C +
C +--2th Matrix Initialisation
C +  -------------------------
C +
        a2( 1) = 0.d0
        b2( 1) = 1.d0
        a2(my) = 0.d0
        b2(my) = 1.d0
C +
      DO j=jp11,my1
        a2(j)   = aa
        b2(j)   = bb
      END DO
C +
C +
C +--1st Equations System
C +  ====================
C +
      DO j=jp11,my1
C +
        d1( 1) =f( 1,j-1)    +2.d0*f( 1,j)    +     f( 1,j+1)
        d1(mx) =f(mx,j-1)    +2.d0*f(mx,j)    +     f(mx,j+1)
C +
        DO  i=ip11,mx1
        d1(i  )=f(i-1,j+1)+2.d0*f(i,j+1)+     f(i+1,j+1)+
     &     2.d0*f(i-1,j)  +4.d0*f(i,j)  +2.d0*f(i+1,j)  +
     &          f(i-1,j-1)+2.d0*f(i,j-1)+     f(i+1,j-1)
        END DO
C +
C +     *********
        call tlat(a1,b1,a1,d1,p1,q1,m  ,m  ,wk1)
C +     *********
C +
        DO i=ip11,mx1
          x(i,j) = wk1(i)
        END DO
C +
      END DO
C +
C +
C +--2th Equations System
C +  ====================
C +
      DO i=ip11,mx1
C +
          d2( 1) = f(i, 1)
          d2(my) = f(i,my)
C +
        DO j=jp11,my1
          d2( j) = x(i,j)
        END DO
C +
C +     *********
        call tlat(a2,b2,a2,d2,p2,q2,mn ,mn ,wk2)
C +     *********
C +
        DO j=jp11,my1
          f(i,j) = wk2(j)
        END DO
C +
      END DO
C +
      return
      end
C +
C +
      subroutine tlat(tlat_a,tlat_b,tlat_c,tlat_d,tlat_p,tlat_q,nx,n 
     .               ,tlat_x)
C +
C +------------------------------------------------------------------------+
C | MAR DYNAMICS FILTER                                    20-09-2001  MAR |
C |   SubRoutine tlat  uses the Gaussian Elimination Algorithm             | 
C |    (e.g. Pielke (1984), pp.302--303)                                   |
C |    (needed to solve the implicit scheme developped for filtering)      |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT:   tlat_a,tlat_b,tlat_c: tri-diagional matrix coefficients     |
C |   ^^^^^    tlat_d              : tri-diagional matrix independent term |
C |            tlat_p,tlat_q       : working          variables            |
C |            n                   : dimension of the variables            |
C |            tlat_x              : variable to solve                     |
C |                                                                        |
C |   OUTPUT:  tlat_x                                                      |
C |   ^^^^^^                                                               |
C +------------------------------------------------------------------------+
C +
      IMPLICIT NONE
C +
      integer  nx,n
      real     tlat_a(nx),tlat_b(nx),tlat_c(nx),tlat_d(nx)
      real     tlat_x(nx),tlat_p(nx),tlat_q(nx)
C + 
      integer  k ,l
C +
C +
C +--Forward  Sweep
C +  ==============
C +
          tlat_p(1)= tlat_b(1)
          tlat_q(1)=-tlat_c(1)/tlat_p(1)
        DO k=2,n
          tlat_p(k)= tlat_a(k)*tlat_q(k-1)+tlat_b(k)
          tlat_q(k)=-tlat_c(k)/tlat_p(k)
        END DO
C +
          tlat_x(1)= tlat_d(1)/tlat_p(1)
        DO k=2,n
          tlat_x(k)=(tlat_d(k)-tlat_a(k)  *tlat_x(k-1))/tlat_p(k)
        END DO
C +
C +
C +--Backward Sweep
C +  ==============
C +
        DO l=2,n
          k=n-l+1
          tlat_x(k)=tlat_q(k)*tlat_x(k+1)+tlat_x(k)
        END DO
C +
      return
      end
