      program PHYradDRIV

! +------------------------------------------------------------------------+
! |   program PHYradDRIV drives MAR radiative routines         15-Dec-2003 |
! |                                                                        |
! |                                                                        |
! +------------------------------------------------------------------------+

      IMPLICIT none

      include 'MARphy.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
      include 'MAR_DY.inc'
      include 'MAR_HY.inc'
      include 'MAR_SL.inc'
      include 'MAR_IO.inc'

      integer  n      ,io


! +--SET-UP
! +  ======

! +--"MAR" Constants
!    ---------------

      zero              =    0.0
      unun              =    1.0
      stefan            =    5.67e-8       ! Stefan-Boltzman
      cp                = 1004.            !
      ra                =  287.            !
      cap               =   ra/cp          !
      gravit            =    9.81          ! Gravity Acceleration

      jtRadi            =    1


! +--"MAR" Vertical Discretization
!    -----------------------------

        DO k=1,mz
           sigma(k)     =        float(k)             / float(mzz)
          sigmid(k)     = 0.5 * (float(k)+float(k-1)) / float(mzz)
        ENDDO
          sigmid(  1)   = 0.
          sigmid(mzz)   = 1.
        DO k=1,mz
          dsigm1(  k)   = sigmid(k+1) - sigmid(k)
        ENDDO


! +--"MAR" Variables
!    ---------------

          iyrrGE        = 1993
          mmarGE        =   12
          jdarGE        =   12
          jhurGE        =   12
          jmmMAR        =    0
          jssMAR        =    0

      DO j=1,my
      DO i=1,mx
          GElatr(i,j)   =   20.*3.1415/180.
          GElonh(i,j)   =    0.
      ENDDO
      ENDDO

          IO_loc        =    3
      DO  io = 1,5
          igrdIO(io)    =    1
          jgrdIO(io)    =    1
      ENDDO

          ptopDY        =    0.001         ! MAR     Top Pressure
      DO j=1,my
      DO i=1,mx
        DO k=1,mz
          tairDY(i,j,k) =  288.-mz+k       !         Air Temperature
        ENDDO
          tairSL(i,j)   =  288.            ! Surface Air Temperature
           pstDY(i,j)   =  100.            ! Surface     Pressure
          eps0SL(i,j)   =    1.            ! Surface     Emissivity
          albeSL(i,j)   =    0.8           ! Surface     Albedo (Snow)
        DO n=1,mw
          tsrfSL(i,j,n) =  288.            ! Surface     Temperature
          SLsrfl(i,j,n) =    1./float(mw)  !
        ENDDO
          maskSL(i,j)   =    0             ! Land(0)/Sea(1) MASK
        DO k=mz-4,mz
            qwHY(i,j,k) =    0.00015
        ENDDO
      ENDDO
      ENDDO

      DO j=1,my
      DO i=1,mx
        DO k=1,mz
            pkDY(i,j,k) = exp(cap *log(pstDYn(i,j)*sigma(k)+ptopDY))
        ENDDO
      ENDDO
      ENDDO

!          ******
      call qsat3D
!          ******


      DO j=1,my
      DO i=1,mx
        DO k=1,mz
            qvDY(i,j,k) =    0.8 * qvswDY(i,j,k)  !      Water Vapor
        ENDDO
      ENDDO
      ENDDO


! +--Radiative Transfert
! +  ===================

      open (unit=4,status='unknown',file='PHYradDRIV.out')
      rewind     4

!          **********
      call PHYrad_CEP
!          **********

      close(unit=4)

      stop
      end

      subroutine PHYrad_top(DistST)

! +------------------------------------------------------------------------+
! |   subroutine PHYrad_top generates the necessary insolation parameters  |
! |                                   for Radiative Transfert routines     |
! |                                                                        |
! +------------------------------------------------------------------------+

      IMPLICIT none

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      real DistST

        DistST      = 1.       ! Earth-Sun Distance    [UA]

      DO j=1,my
      DO i=1,mx
        czenGE(i,j) = 0.100    ! Cos Zenithal Distance  [-]
      ENDDO
      ENDDO

      return
      END

      subroutine TIMcor

! +------------------------------------------------------------------------+
! |   subroutine TIMcor                                                    |
! |                                                                        |
! +------------------------------------------------------------------------+

      include 'MARdim.inc'
      include 'MAR_GE.inc'

      return
      end
      subroutine qsat3D
C +
C +------------------------------------------------------------------------+
C | MAR PHYSICS                                            14-11-2002  MAR |
C |   SubRoutine qsat3D computes the Saturation Specific Humidity  (kg/kg) |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT :   TairSL: Surface Air Temperature                        (K) |
C |   ^^^^^^^   TairDY:         Air Temperature                        (K) |
C |              pstDY: Model Pressure Thickness                     (kPa) |
C |                                                                        |
C |   OUTPUT :  qvswDY: Saturation Specific Humidity  over Water   (kg/kg) |
C |   ^^^^^^^   qvsiDY: Saturation Specific Humidity  over Ice     (kg/kg) |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
      IMPLICIT NONE
C +
C +
C +--Global Variables
C +  ================
C +
      include 'MARphy.inc'
C +
      include 'MARdim.inc'
      include 'MARgrd.inc'
C +
      include 'MAR_DY.inc'
      include 'MAR_SL.inc'
      include 'MAR_WK.inc'
C +
C +
C +--Local  Variables
C +  ================
C +
      real     WatIce,ExpWat,ExpWa2,ExpIce
C +
C +
C +--DATA
C +  ====
C +
      data     WatIce/273.16e0/
      data     ExpWat/5.138e0/
      data     ExpWa2/6827.e0/
      data     ExpIce/6150.e0/
C +
C +
C +--Temperature (K) and Pressure (hPa)
C +  ==================================
C +
      DO   k=1,mz
        DO j=1,my
        DO i=1,mx
          WKxyz5(i,j,k)   = tairDY(i,j,k)
          WKxyz6(i,j,k)   = (pstDY(i,j)  * sigma(k) + ptopDY) * 10.0d0
        END DO
        END DO
      END DO
C +
        DO j=1,my
        DO i=1,mx
          WKxyz5(i,j,mzz) = TairSL(i,j)
          WKxyz6(i,j,mzz) = (pstDY(i,j)             + ptopDY) * 10.0d0
        END DO
        END DO
C +
C +
C +--Saturation Vapor Pressure over Ice
C +  ==================================
C +
      DO   k=1,mzz
        DO j=1,my
        DO i=1,mx
          WKxyz7(i,j,k) = 
     .    6.1070d0 * exp (ExpIce*(unun/WatIce-unun/WKxyz5(i,j,k)))
C +...    Dudhia (1989) MWR, (B1) and (B2) p.3103 
C +
          WKxyz8(i,j,k) = .622d0*WKxyz7(i,j,k)
     .  /(WKxyz6(i,j,k) - .378d0*WKxyz7(i,j,k))
C +
C +
C +--Saturation Vapor Pressure over Water
C +  ====================================
C +
          WKxyz7(i,j,k) = 
     .    6.1078d0 * exp (ExpWat*  log(WatIce     /WKxyz5(i,j,k)))
     .             * exp (ExpWa2*(unun/WatIce-unun/WKxyz5(i,j,k)))
C +...    Dudhia (1989) MWR, (B1) and (B2) p.3103 
C +       See also Pielke (1984), p.234 and Stull (1988), p.276 
C +
          qvswDY(i,j,k) = .622d0*WKxyz7(i,j,k)
     .  /(WKxyz6(i,j,k) - .378d0*WKxyz7(i,j,k))
C +...    Saturation Vapor Specific Concentration over Water
C +       (even for temperatures less than freezing point)
C +
C +
C +
C +--Water Phase Discriminator
C +  =========================
C +
          WKxyz7(i,j,k) = max(zero,sign(unun,WKxyz5(i,j,k)-WatIce))
C +...    WKxyz7(i,j,k) =     1    if        Tair     >    273.16
C +                           0    if        Tair     <    273.16
C +
C +
C +--Saturation Vapor Specific Concentration over Ice
C +  ================================================
C +
          qvsiDY(i,j,k) = qvswDY(i,j,k) *      WKxyz7(i,j,k) 
     .                  + WKxyz8(i,j,k) *(unun-WKxyz7(i,j,k))
C +
C +
C +--Work Area Reset
C +  ===============
C +
          WKxyz5(i,j,k) = 0.0d+0
          WKxyz6(i,j,k) = 0.0d+0
          WKxyz7(i,j,k) = 0.0d+0
          WKxyz8(i,j,k) = 0.0d+0
        END DO
        END DO
      END DO
C +
      return
      end      
