C   +-------------------------------------------------------------------+
C   |  Subroutine SL_cor                               Sept 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | This routine is designed to correct the vertical profile of wind, |
C   | temperature and moisture in the surface layer.                    |
C   |                                                                   |
C   | Input : NST_pt : potential temperature     (K)                    |
C   | ^^^^^^^ NST_st : surface temperature       (K)                    |
C   |         NST_qv : specific humidity         (kg/kg)                |
C   |         NST_z0 : roughness length          (m)                    |
C   |         NST__u : U-wind                    (m/s)                  |
C   |         NST__v : V-wind                    (m/s)                  |
C   |         NST_sp : surface pressure          (kPa)                  |
C   |         NST__p : pressure at each level    (kPa)                  |
C   |                                                                   |
C   | Output: Corrected NST__u, NST__v, and NST_pt in the surface layer.|
C   | ^^^^^^^ Optional (#QV) for NST_qv.                                |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE SL_cor


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'

      INTEGER i,j,k,level

      REAL tair(mz),diff,karman,ustar,vstar,SLheight,deltaz,
     .     plevel,ra,gravit,cap


C +---Data
C +   ---- 

      DATA karman    /   0.4 /
      DATA ra        /  287. /
      DATA gravit    /  9.81 /
      DATA SLheight  /   40. /
      DATA cap       /  0.285856574 /

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO j=1,my    ! horizontal loop on grid points
      DO i=1,mx    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Compute real temperature
C +   ========================

      DO k=1,mz
       tair(k)=NST_pt(i,j,k)/(100./NST__p(i,j,k))**cap
      ENDDO


C +---Search for the top level in the surface layer
C +   =============================================

       plevel=NST_sp(i,j)*exp(-gravit/ra/tair(mz)*SLheight)

       diff=100.

       DO k=1,mz
        IF (abs(NST__p(i,j,k)-plevel).lt.diff) THEN
         level=k
         diff =abs(NST__p(i,j,k)-plevel)
        ENDIF
       ENDDO

       level=min(level,mz-1)
       level=max(level,1)


C +---Compute friction velocity in neutral atmosphere
C +   ===============================================

       ustar=karman*NST__u(i,j,level)/log(SLheight/NST_z0(i,j))
       vstar=karman*NST__v(i,j,level)/log(SLheight/NST_z0(i,j))


C +---CORRECTION OF WIND, TEMPERATURE AND SPECIFIC HUMIDITY IN SL
C +   ===========================================================


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       DO k=level+1,mz   ! loop on levels in the surface layer

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        deltaz=-ra*tair(k)/gravit*log(NST__p(i,j,k)/NST_sp(i,j))
C +...  Height above the surface for the level k

        IF (deltaz.lt.SLheight) THEN


C +---LOG vertical profile for velocities in the surface layer
C +   --------------------------------------------------------

         NST__u(i,j,k)=ustar/karman*log(deltaz/NST_z0(i,j))
         NST__v(i,j,k)=vstar/karman*log(deltaz/NST_z0(i,j))


C +---Vertical temperature profile in the surface layer
C +   -------------------------------------------------

         IF (NST_pt(i,j,level).lt.NST_st(i,j)) THEN

C +       ... Constant potential temperature if unstable layer

          NST_pt(i,j,k)=NST_pt(i,j,level)

         ENDIF


C +---Vertical profile of specific humidity (constant Qv)
C +   ---------------------------------------------------

c #QV    NST_qv(i,j,k)=NST_qv(i,j,level)

        ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ENDDO     ! loop on levels in the surface layer

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDDO      ! horizontal loop on grid points
      ENDDO      ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END

