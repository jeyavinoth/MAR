C   +-------------------------------------------------------------------+
C   |  Subroutine USRgrd                             July 2012  NESTING |
C   +-------------------------------------------------------------------+
C   | USRgrd adapt NESTOR to Greenland region                           |
C   |                                                                   |
C   | Input : - subnam : Name of the subroutine                         |
C   | ^^^^^^^            where USRgrd is called                         |
C   |                                                                   |
C   | Maintainer : Xavier Fettweis                                      |
C   | ^^^^^^^^^^^^                                                      |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE USRgrd (subnam)


      IMPLICIT NONE

C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'LSCvar.inc'

C +---local variables
C +   ---------------

      CHARACTER*6 subnam

      INTEGER     nsvat,nigbp
      PARAMETER  (nsvat=12)
      PARAMETER  (nigbp=17)

      INTEGER     i,j,k,l,var2(mx,my)
  
      REAL        SVAT(0:nsvat),IGBP(nigbp),convert(nigbp,0:nsvat),
     .            svat_frac (3),iIGBP(nigbp),igbp_z0(nigbp),
     .            tmp1,tmp2,ELA,var1(mx,my), svat_class(3),frac_tot

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Topography for ETOPOg    
C +   =====================

      IF (subnam.eq.'ETOPOg') THEN

      write (*,*) 'Special topo for Greenland Simulation'
      write (*,*)

      call CF_READ2D("input/TOPO/GRD-25km-80x135.cdf"
     .               ,'MSK',1,mx,my, 1,var1)

      do i=1,mx ; do j=1,my
       if(var1(i+0,j+0)==1) then
        NSTsol(i,j)=4
       else
        NSTsol(i,j)=1
       endif
      enddo ; enddo

      call CF_READ2D("input/TOPO/GRD-25km-80x135.cdf "
     .               ,'ICE',1,mx,my, 1,var1)
      
      do i=1,mx ; do j=1,my
       NSTice(i,j)=var1(i+0,j+0)*100.
      enddo     ; enddo

      call CF_READ2D("input/TOPO/GRD-25km-80x135.cdf "
     .               ,'SRF',1,mx,my, 1,var1)

      do i=1,mx ; do j=1,my
       NST_sh(i,j)=var1(i+0,j+0)

      enddo     ; enddo


      ENDIF        

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Soil Type for GLOveg    
C +   ====================~

C +---GRD 1: Initialisation of surface variables
C +   ..........................................

      IF (subnam.eq.'GLOveg') THEN

      write(6,*) 'Global land cover (IGBP) over Greenland Region'      
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

      DO j=2,my-1
      DO i=2,mx-1

       NSTfrc(i,j)     = 0.
       
       DO k=1,nvx
        NSTsvt(i,j,k)  = 0 
        NSTsfr(i,j,k)  = 0
        NSTveg(i,j,k)  = 0
        NSTvfr(i,j,k)  = 0
       END DO
       
       IF (NSTsol(i,j).ge.4) THEN
        
        NSTsvt(i,j,1)  = 4
        NSTsfr(i,j,1)  = 100
        NSTveg(i,j,1)  = 9
        NSTvfr(i,j,1)  = 100
        
        DO k=1,nvx
         NSTlai(i,j,k) = 2.0
         NSTglf(i,j,k) = 1.0
        END DO
      
       ELSE
      
        NSTsvt(i,j,1)  = 0
        NSTsfr(i,j,1)  = 100
        NSTveg(i,j,1)  = -1
        NSTvfr(i,j,1)  = 100
      
        DO k=1,nvx
         NSTlai(i,j,k) = 0.0
         NSTglf(i,j,k) = 0.0
        END DO
      
        NST_z0(i,j)    = 0.0013
      
       ENDIF

      END DO
      END DO

C +---IGBP Surface variables
C +   ======================

      DO j=2,my-1
      DO i=2,mx-1 
      IF (NSTsol(i,j) .ge. 4 ) THEN
      
       DO k=1,nigbp
       DO l=0,nsvat
        SVAT(l)       = 0.
        IGBP(k)       = 0.
        convert(k,l)  = 0.
       ENDDO
       ENDDO
                 
       IF(NST__x(i,j).gt.-43) then ! Equilibrium line (m)
        ELA           = -32759.680d0 + 1001.782d0  * NST__y(i,j)
     .                  - 7.331d0    * NST__y(i,j) * NST__y(i,j)
       ELSE
        ELA           = -23201.445d0 + 746.249d0   * NST__y(i,j)
     .                  - 5.640d0    * NST__y(i,j) * NST__y(i,j)
       END IF

       IF (nvx .eq. 2) THEN
              
       convert(16, 4) = 100.   ! grass low
       convert(16, 0) = 0.     ! barren soil
       igbp_z0(16   ) = 0.022

       convert(15, 4) = 0.     ! grass low
       convert(15, 0) = 100.   ! barren soil
       igbp_z0(15   ) = 0.001       

       if(NSTice(i,j).ge.0) then
        NSTveg(i,j,1) = 15 
        NSTveg(i,j,2) = 16        
        NSTvfr(i,j,1) = NSTice(i,j)
        NSTvfr(i,j,2) = 100.0 - NSTvfr(i,j,1)
       else
        NSTveg(i,j,1) = 15       
        NSTvfr(i,j,1) = 100.0
        NSTvfr(i,j,2) = 0.0      
       end if

        IGBP(15)      = NSTvfr(i,j,1) / 100.0       
        IGBP(16)      = NSTvfr(i,j,2) / 100.0   
                  
       END IF
          
C +...  convertion to SVAT
C +     ~~~~~~~~~~~~~~~~~~
        DO k=1,nigbp
         DO l=0,nsvat
          SVAT(l)=SVAT(l)+convert(k,l)*IGBP(k)
         ENDDO
        ENDDO

C +...  retain the (nvx-1) dominant classes
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO k=2,nvx
         svat_class(k)=1
         svat_frac (k)=SVAT(1)
         DO l=1,nsvat
          IF (svat_frac(k).lt.SVAT(l)) THEN
           svat_class(k)=l
           svat_frac (k)=SVAT(l)
          ENDIF
          SVAT(svat_class(k))=0.
         ENDDO
        ENDDO

C +...  class (nvx) is reserved for barren soil
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        svat_class(1) = 0
        svat_frac (1) = SVAT(0)

C +...  normalizing the three dominant fractions
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        frac_tot=0.
        DO l=1,nvx
         frac_tot=frac_tot+svat_frac(l)
        ENDDO
        IF (frac_tot.ne.0.) THEN
         DO l=1,nvx
          svat_frac(l)=svat_frac(l)/frac_tot
         ENDDO
        ENDIF
C +...  attribute classes and fractions to NST variables
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO k=1,nvx
         NSTsvt(i,j,k)=     svat_class(k)
         NSTsfr(i,j,k)=min(100.,svat_frac (k)*100.)
        ENDDO

        DO l=1,nvx

         IF (NSTsvt(i,j,l).eq. 0) NSTlmx(i,j,l) = 0.0
         IF (NSTsvt(i,j,l).eq. 1) NSTlmx(i,j,l) = 0.6
         IF (NSTsvt(i,j,l).eq. 2) NSTlmx(i,j,l) = 0.9
         IF (NSTsvt(i,j,l).eq. 3) NSTlmx(i,j,l) = 1.2
         IF (NSTsvt(i,j,l).eq. 4) NSTlmx(i,j,l) = 0.7
         IF (NSTsvt(i,j,l).eq. 5) NSTlmx(i,j,l) = 1.4
         IF (NSTsvt(i,j,l).eq. 6) NSTlmx(i,j,l) = 2.0
         IF (NSTsvt(i,j,l).eq. 7.or.NSTsvt(i,j,l).eq.10)
     .    NSTlmx(i,j,l) = 3.0
         IF (NSTsvt(i,j,l).eq. 8.or.NSTsvt(i,j,l).eq.11)
     .    NSTlmx(i,j,l) = 4.5
         IF (NSTsvt(i,j,l).eq. 9.or.NSTsvt(i,j,l).eq.12)
     .    NSTlmx(i,j,l) = 6.0

         NSTlai(i,j,l) = NSTlmx(i,j,l)
         NSTglf(i,j,l) = 1.0

        ENDDO

          
      END IF 
      END DO 
      END DO

      ENDIF

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      END SUBROUTINE
      
