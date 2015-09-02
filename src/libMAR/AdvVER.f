      subroutine AdvVER(ff,dtProc,NoPass,iniBAL,makBAL,iterat,i_step
     .                 ,lb)
      
C +----------------------------------------------------------------------------+
C |                                                                            |
C |                                                                            |
C |   subroutine AdvVER verifies the Leap-Frog Advection Scheme                |
C |                                                                            |
C |                                                                            |
C +----------------------------------------------------------------------------+


      IMPLICIT NONE


      include "MARphy.inc"
      include "MARdim.inc"
      include "MARgrd.inc"
      include "MAR_DY.inc"

      real     ff(mx,my,mzz)
      real     dtProc,NoPass
      logical  iniBAL,makBAL
      integer  iterat,i_step
      character*6     lb


C +--Local   Variables
C +  =================

      integer   LD
      parameter(LD  =12)
      integer   ngx ,ndx ,ngy ,ndy ,npx ,npy ,nxy
      parameter(ngx = min(LD+1,mx     ))
      parameter(ndx = max(   1,mx-LD-2))
      parameter(ngy = min(LD+1,my     ))
      parameter(ndy = max(   1,my-LD-2))
      parameter(npx = max(1,mx-LD-LD-2))
      parameter(npy = max(1,my-LD-LD-2))
      parameter(nxy = npx    *npy      )
      integer   ngx1 ,ndx1 ,ngy1 ,ndy1
      integer   ngx2 ,ndx2 ,ngy2 ,ndy2
      parameter(ngx1= max(1,ngx-1   ))
      parameter(ngx2= max(1,ngx-2   ))
      parameter(ndx1= min(  ndx+1,mx))
      parameter(ndx2= min(  ndx+2,mx))
      parameter(ngy1= max(1,ngy-1   ))
      parameter(ngy2= max(1,ngy-2   ))
      parameter(ndy1= min(  ndy+1,my))
      parameter(ndy2= min(  ndy+2,my))
      integer   ngxx ,ndxx ,ngyy ,ndyy
      parameter(ngxx= min(  ngx+1,mx))
      parameter(ndxx= max(1,ndx-1   ))
      parameter(ngyy= min(  ngy+1,my))
      parameter(ndyy= max(1,ndy-1   ))

      integer         npt

      real*8          ff_int(mx,my)
      real*8          ff_BAL
      real*8          ff_LBC
      common/AdvVER_r/ff_BAL,ff_LBC

      real*8          p_star(mx,my)    ! Local       Mass
      real*8          qt(    -2: 2)    !
      real*8          qFag  ,qFad      !
      real*8          qFagi ,qFadi     !
      real*8          qFdg  ,qFdd      !
      real*8          f2_3

      logical         ps_var,iniWRI,zeroLB

      data            ps_var/.TRUE./
      data            iniWRI/.TRUE./
      data            zeroLB/.TRUE./

        IF (iniWRI)                                                 THEN
            iniWRI = .FALSE.
            write(6,6001) ngx2,ngx1,ngx,ngxx,ndxx,ndx,ndx1,ndx2
 6001       format(/,' xBC:',4i6,6x,4i6)
            write(6,6002) ngy2,ngy1,ngy,ngyy,ndyy,ndy,ndy1,ndy2
 6002       format(  ' yBC:',4i6,6x,4i6,/,)
        END IF


C +--Water Mass Increment (vertical Integral)                      [mWE]
C +  ========================================                      =====

            npt = (npx+npx)*min(1,max(0,my-1))
     .          + (npy+npy)*min(1,max(0,mx-1))
            f2_3   = 2./3.
        IF (NoPass.LT.0.)                                           THEN
            ff_BAL = 0.
            ff_LBC = 0.
         IF(zeroLB)                                                 THEN
          DO k=  1,mz
          DO j=  1,my
          DO i=ngx-3,ngx+3
            ff(i,j,k) = ff(ngx,j  ,k)
          ENDDO
          DO i=ndx-3,ndx+3
            ff(i,j,k) = ff(ngx,j  ,k)
          ENDDO
          ENDDO
          DO i=  1,mx
          DO j=ngy-3,ngy+3
            ff(i,j,k) = ff(i  ,ngy,k)
          ENDDO
          DO j=ndy-3,ndy+3
            ff(i,j,k) = ff(i  ,ndy,k)
          ENDDO
          ENDDO
          DO i=ngx-3,ngx+3
          DO j=ngy-3,ngy+3
            ff(i,j,k) = ff(ngx,ngy,k)
          ENDDO
          ENDDO
          DO i=ngx-3,ngx+3
          DO j=ndy-3,ndy+3
            ff(i,j,k) = ff(ngx,ndy,k)
          ENDDO
          ENDDO
          DO i=ndx-3,ndx+3
          DO j=ngy-3,ngy+3
            ff(i,j,k) = ff(ndx,ngy,k)
          ENDDO
          ENDDO
          DO i=ndx-3,ndx+3
          DO j=ndy-3,ndy+3
            ff(i,j,k) = ff(ndx,ndy,k)
          ENDDO
          ENDDO
          ENDDO
         END IF
        END IF

        IF (ps_var)                                                 THEN
          DO j=  1,my 
          DO i=  1,mx 
            p_star(i,j) = pstDYn(i,j)
          ENDDO
          ENDDO
        ELSE
          DO j=  1,my 
          DO i=  1,mx 
            p_star(i,j) = 100.00
          ENDDO
          ENDDO
        END IF


        DO j=ngy,ndy
        DO i=ngx,ndx

            ff_int(i,j)   = 0.
          DO k=1,mz
            ff_int(i,j)   = ff_int(i,j) + ff(i,j,k)*dsigm1(k)
          ENDDO
            ff_int(i,j)   = ff_int(i,j)  *p_star(i,j)*grvinv
            ff_BAL        = ff_BAL       +ff_int(i,j)*NoPass

        ENDDO
        ENDDO


C +--LBC
C +  ===

C +--Advection                                               [mWE m2/s3]
C +  ~~~~~~~~~                                               ~~~~~~~~~~~
        IF (NoPass.LT.0.)                                           THEN
          IF (mx .GT. 1)                                            THEN
            DO j=ngy,ndy
               qFag      = 0.
               qFad      = 0.
               qFagi     = 0.
               qFadi     = 0.
             DO k=1,mz

C +--Water Concentration, x< Boundary
C +  
               qt(-2)    =   ff(ngx2,j,k) * p_star(ngx2,j)
               qt(-1)    =   ff(ngx1,j,k) * p_star(ngx1,j)
               qt( 0)    =   ff(ngx ,j,k) * p_star(ngx ,j)

C +--Water Fluxes,        x< Boundary
C +  
               qFag      =   qFag      +uairDY(ngx ,j,k) *dsigm1(k)
     .                   *f2_3*((qt(-1)-qt( 0)) +0.125*(qt( 0)-qt(-2)))
               qFagi     =   qFagi     +uairDY(ngxx,j,k) *dsigm1(k)
     .                   *f2_3*(                 0.125*(qt( 0)-qt(-1)))

C +--Water Concentration, x> Boundary
C +  
               qt( 0)    =   ff(ndx ,j,k) * p_star(ndx ,j)
               qt( 1)    =   ff(ndx1,j,k) * p_star(ndx1,j)
               qt( 2)    =   ff(ndx2,j,k) * p_star(ndx2,j)

C +--Water Fluxes,        x> Boundary
C +  
               qFad      =   qFad      -uairDY(ndx ,j,k) *dsigm1(k)
     .                   *f2_3*((qt( 1)-qt( 0)) +0.125*(qt( 0)-qt( 2)))
               qFadi     =   qFadi     -uairDY(ndxx,j,k) *dsigm1(k)
     .                   *f2_3*(                 0.125*(qt( 0)-qt( 1)))

             ENDDO

               ff_LBC    = ff_LBC + qFag + qFad + qFagi + qFadi !  [mWE m2/s3]
            ENDDO


C +--Water Concentration, y< Boundary
C +  
           IF (my .GT. 1)                                           THEN
            DO i=ngx,ndx
               qFag      = 0.
               qFad      = 0.
               qFagi     = 0.
               qFadi     = 0.
             DO k=1,mz
               qt(-2)    =   ff(i,ngy2,k) * p_star(i,ngy2)
               qt(-1)    =   ff(i,ngy1,k) * p_star(i,ngy1)
               qt( 0)    =   ff(i,ngy ,k) * p_star(i,ngy )

C +--Water Fluxes,        y< Boundary
C +  
               qFag      =   qFag      +vairDY(i,ngy ,k) *dsigm1(k)
     .                   *f2_3*((qt(-1)-qt( 0)) +0.125*(qt( 0)-qt(-2)))
               qFagi     =   qFagi     +vairDY(i,ngyy,k) *dsigm1(k)
     .                   *f2_3*(                 0.125*(qt( 0)-qt(-1)))

C +--Water Concentration, y> Boundary
C +  
               qt( 0)    =   ff(i,ndy ,k) * p_star(i,ndy )
               qt( 1)    =   ff(i,ndy1,k) * p_star(i,ndy1)
               qt( 2)    =   ff(i,ndy2,k) * p_star(i,ndy2)

C +--Water Fluxes,        y> Boundary
C +  
               qFad      =   qFad      -vairDY(i,ndy ,k) *dsigm1(k)
     .                   *f2_3*((qt( 1)-qt( 0)) +0.125*(qt( 0)-qt( 2)))
               qFadi     =   qFadi     -vairDY(i,ndyy,k) *dsigm1(k)
     .                   *f2_3*(                 0.125*(qt( 0)-qt( 1)))

             ENDDO
               ff_LBC    = ff_LBC + qFag + qFad + qFagi + qFadi !  [mWE m2/s3]
            ENDDO
           END IF
          END IF
        END IF


C +--OUTPUT
C +  ======

        IF (NoPass.GT.0.)                                           THEN
          ff_BAL = ff_BAL                   /    nxy
          ff_LBC = ff_LBC * dtProc * grvinv /(dx*nxy)
          write(6,6000) iterat,i_step,lb,ff_BAL,ff_LBC,ff_BAL-ff_LBC
 6000     format(2i6,3x,a6,': HYDver: ',f21.12,' - ',f21.12,
     .                                         ' = ',f21.12)
        END IF

      return
      end
