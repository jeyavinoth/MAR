      subroutine SBCnew

C +---------------------------------------------------------------------------+
C | MAR Surface Boundary Conditions                            16-03-2004 MAR |
C |                                                                           |
C |   subroutine SBCnew solves analytically Non Linear Mountain Waves         |
C |                                                                           |
C +---------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARCTR.inc'
      include 'MARphy.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
      include 'MAR_DY.inc'
      include 'MAR_NH.inc'

      REAL              sh0(mx),sh_max           ! Topography                [m]
      REAL              sua2ub                   ! Adjustment Time (MAX)     [s]
      REAL              tua2ub                   ! Adjustment Time (CURRENT) [s]
      REAL              uubuum ,uufact           ! Adjustment Factor
      common/SBCnew_rea/sh0    ,sh_max           !
     .                 ,sua2ub ,tua2ub           !
     .                 ,uubuum ,uufact           !


C +--Local   Variables
C +  =================

      real          hhmaam                       ! hhm X aam
      real          N_BruV                       ! Brunt-Vaisala Frequency [1/s]
      real          Scorer                       ! Scorer        Parameter
      real          L2                           ! 
      REAL          xx(mx)                       ! Distance                  [m]


C +--NH    Variables
C +  ---------------

      integer            mf                      !
      parameter         (mf=mx)                  !
      real          dk                           ! Wave Number Increment   [1/m]
      real          kk(0:mf)                     ! Wave Number             [1/m]
      real          k2L2  ,kkll  ,arg__z         !
      real          zz_m  ,dz_m  ,dz_min         !
      parameter                  (dz_min=20.)    ! Min. Vert.  Grid Spacing
      integer       imi   ,ipi   ,n              !
      real          ddelta(mx,mz,-1:1)           ! StreamLine  Displacement
      real          zz_m_n(mx,mz,-1:1)           ! Altitude of Computations
      real          uAnaNH(mx,mz)                ! Analyt.Hor. Perturb.Veloc.
      real          wAnaNH(mx,mz)                ! Analyt.Vert.Perturb.Veloc.
      common/sAnaNH/uAnaNH,wAnaNH                !

      real          paNH  ,paNHmz,pa_Hmz         ! Statistics


C +--Non-Linear Mountain Wave DATA Parameters
C +  ----------------------------------------

      logical       Ana2MAR                      ! Solution Switch
      real          hhm                          ! Mountain Height           [m]
      real          aam                          ! Mountain Half-Width       [m]
      REAL          uum                          ! Initial Base Wind Speed [m/s]
      real          uub                          ! MAXIMUM Base Wind Speed [m/s]
      real          BVF                          ! Brunt-Vaisala Frequency [1/s]


C +--DATA
C +  ====

      include "SBCnew.dat"                       ! Mountain Wave DATA Parameters


C +   ================                                              ====
      IF (itexpe.GE.0)                                              THEN
C +   ================                                              ====


        sua2ub = 5.0 * aam / uub                 ! Adjustment Time           [s]
        tua2ub = 0.0                             ! (ARPS.4.0 Users Guide, p.312)


C +---Altitude
C +   ========

        sh_max =             0.0
      DO i=1,mx
        sh0(i) =     sh(i,1)
        sh_max = max(sh(i,1),sh_max)
      ENDDO


C +---Parameters of the Non-Hydrostatic Non-Linear Analytical Solution
C +   ================================================================

      N_BruV = BVF                               ! Brunt-Vaisala Frequency [1/s]
      Scorer = N_BruV / ugeoDY(1,1,1)            ! Scorer Parameter
      L2     = Scorer * Scorer                      

      hhmaam = hhm    * aam                      !

      DO i=1,mx
        xx(i)   = 1.0e3* xxkm(i) *sign(1.00,uum)
        kk(i)   = 2.0  * pi /(dx *    (mx +1 -i))
      ENDDO
        kk(0)   = 0.0


C +--Non-Hydrostatic Solution (solution to the   Long's 1953 equation)
C +  ======================== (see ARPS 4.0 User's Guide 13.8.2 p.317)

C +--Streamline Displacement
C +  -----------------------

         j=1
      DO k=1,mz
      DO i=1,mx
          DO n=-1,1
                ddelta(i,k,n)  = 0.
          ENDDO
        DO m= 1,mf
c #MF           kk(m)          = m    *50./mf    ! Regular Distribution
                dk             = kk(m)-kk(m-1)

                k2L2           = kk(m)*kk(m)                  -L2
                zz_m           =           gplvDY(i,j,k)  /gravit
                dz_m           =  min(zz_m-gplvDY(i,j,k+1),dz_min)
                zz_m_n(i,k,-1) =      zz_m-dz_m
                zz_m_n(i,k, 0) =      zz_m
                zz_m_n(i,k, 1) =      zz_m+dz_m

          DO n=-1,1
                arg__z        = sqrt(max(0.,k2L2)) * zz_m_n(i,k,n)  
     .                                              +kk(m)*aam
            IF (kk(m).LE.Scorer)                                    THEN
                kkll          = sqrt(      -k2L2)  * zz_m_n(i,k,n)
                ddelta(i,k,n) = ddelta(i,k,n) + exp(-kk(m)*aam         )
     .                                         *cos( kk(m)*xx(i) + kkll)
     .                                         *     dk
            ELSE
                kkll          = sqrt(       k2L2)  * zz_m_n(i,k,n)  
     .                                             + kk(m)*aam
                ddelta(i,k,n) = ddelta(i,k,n) + exp(-kkll              )
     .                                         *cos( kk(m)*xx(i)       )
     .                                         *     dk
            END IF
          ENDDO

        ENDDO
          DO n=-1,1
                ddelta(i,k,n) = ddelta(i,k,n)  *     hhmaam
          ENDDO
      ENDDO 
      ENDDO 


C +--Streamline Displacement
C +  -----------------------

         j=1
      DO k=1,mz
      DO i=1,mx
        ipi           =           min(i+1,mx)
        imi           =           max(i-1, 1)
        uAnaNH(i,k)   = - uum*(ddelta(i  ,k, 1)-ddelta(i  ,k,-1))
     .                       /(zz_m_n(i  ,k, 1)-zz_m_n(i  ,k,-1))
        wAnaNH(i,k)   = - uum*(ddelta(ipi,k, 0)-ddelta(imi,k, 0))/dx
      ENDDO 
      ENDDO 

C +   ======
      END IF
C +   ======


C +--Solution Forcing
C +  ================

      IF      (uub.LT.uum.AND.tua2ub.LT.sua2ub)                     THEN

          uubuum        =                 (uub/ uum)*(tua2ub/sua2ub)
          uufact        =                 (uub/ uum)*(dt    /sua2ub)
           j=1
        DO k=1,mz
        DO i=1,mx
          tua2ub        = tua2ub         + dt
          uairDY(i,j,k) = uairDY(i,j,k)  * uufact
          ugeoDY(i,j,k) = ugeoDY(i,j,k)  * uufact
        ENDDO 
        ENDDO 

      ELSE IF (Ana2MAR)                                             THEN

           j=1
        DO k=1,mz
        DO i=1,mx
          uairDY(i,j,k) = uairDY(i,j,k) 
     .     + (uAnaNH(i,k)-uairDY(i,j,k)) * 0.1
        ENDDO 
        ENDDO 

      END IF


C +--Statistics
C +  ==========

        paNH   = 0.
        paNHmz = 0.
        pa_Hmz = 0.
      DO j=1,my
      DO i=1,mx
      DO k=1,mz
        paNH   = paNH   + pairNH(i,j,k)  *pstDYn(i,j) *dsigm1(k)
      ENDDO 
        paNHmz = paNHmz + pairNH(i,j,mz) *pstDYn(i,j)
        pa_Hmz = pa_Hmz +                 pstDYn(i,j)
      ENDDO 
      ENDDO 
      IF (mod(itexpe,20).EQ.0)                                      THEN
        write(6,6000)
 6000   format(' hh:mm:ss  hPa Domain',9x,'SBL',' Hydrostatic'
     .                                          '  Wind Speed')
      ENDIF
        paNH   = 10.*paNH   /(mx*my)
        paNHmz = 10.*paNHmz /(mx*my)
        pa_Hmz = 10.*pa_Hmz /(mx*my)
        write(6,6001) jhurGE,minuGE,jsecGE,paNH  ,paNHmz,pa_Hmz
     .               ,uairDY(imez,jmez,1) ,uubuum
 6001   format(i3,2(':',i2),5f12.3)

      return
      end
