C   +-------------------------------------------------------------------+
C   |  Subroutine USRgrd                           February 04  NESTING |
C   +-------------------------------------------------------------------+
C   | USRgrd adapt NESTOR to Greenland region                           |
C   |                                                                   |
C   | Input : - subnam : Name of the subroutine                         |
C   | ^^^^^^^            where USRgrd is called                         |
C   |                                                                   |
C   | Maintainer : Emilie Vanvyve                                       |
C   | ^^^^^^^^^^^^                                                      |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE USReur (subnam)


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

      INTEGER     nsvat,nigbp,frac_tot
      PARAMETER  (nsvat=12)
      PARAMETER  (nigbp=17)

      INTEGER     i,j,k,l,svat_class(3),ii,jj
  
      REAL        SVAT(0:nsvat),IGBP(nigbp),convert(nigbp,0:nsvat),
     .            svat_frac (3),iIGBP(nigbp),igbp_z0(nigbp),
     .            tmp1,tmp2,ELA,ww  

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF (subnam.eq.'ETOPOg') THEN

      write (*,*) 'Special topo for Mont Rigi'
      write (*,*)

      go to 100      

      ii=46 ; jj=24      

      do i=-1,1 ; do j=-1,1

                          ww=0.05
       if (i==0.or. j==0) ww=0.10
       if (i==0.and.j==0) ww=0.20

       NST_sh(i+ii,j+jj)=NST_sh(i+ii,j+jj)*(1+ww)

      enddo ; enddo

100   continue

      endif


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Soil Type for GLOveg    
C +   ====================~

C +---GRD 1: Initialisation of surface variables
C +   ..........................................

      IF (subnam.eq.'GLOveg') THEN

      go to 200

C +---IGBP Surface variables
C +   ======================

      DO j=2,my-1
      DO i=2,mx-1 

      IF(NSTsol(i,j)==3) NSTsol(i,j)=4

      IF (NSTsol(i,j) .ge. 4    .and.
     .    NST__x(i,j) .ge. -30. .and. 
     .    NST__x(i,j) .le. -10. .and. 
     .    NST__y(i,j) .ge.  60. .and. 
     .    NST__y(i,j) .le.  70.) THEN
      
       DO k=1,nigbp
       DO l=0,nsvat
        SVAT(l)       = 0.
        IGBP(k)       = 0.
        convert(k,l)  = 0.
       ENDDO
       ENDDO
          
       ELA=5000  
 
       If (NST_sh(i,j).ge.ELA) NSTsol(i,j)=3
          
       IF (nvx .eq. 3) THEN
       
       convert( 7, 5) = 30.   ! grass medium
       convert( 7, 7) = 40.   ! broadleaf low    TUNDRA
       convert( 7, 8) = 30.   ! broadleaf medium
       igbp_z0( 7   ) = 0.33
       
       convert(16, 4) = 20.   ! grass low
       convert(16, 7) = 5.    ! broadleaf low  MOUNTAIN
       convert(16, 0) = 75.   ! barren soil
       igbp_z0(16   ) = 0.022
       NSTveg(i,j,nvx)= -1 
       NSTvfr(i,j,nvx)= 0 
       
       if (NST_sh(i,j) .le. ELA) then
        NSTveg(i,j,1) = 7 
        NSTveg(i,j,2) = 16        
        NSTvfr(i,j,1) = 100.0 * (1- NST_sh(i,j)/ ELA)
        NSTvfr(i,j,2) = 100.0 - NSTvfr(i,j,1)
       else
        NSTveg(i,j,1) = 7  
        NSTveg(i,j,2) = 16        
        NSTvfr(i,j,1) = 0.0
        NSTvfr(i,j,2) = 100.0      
       end if
                     
       if (NSTvfr(i,j,2) .gt. NSTvfr(i,j,1)) then
        NSTveg(i,j,1) = 16 
        NSTveg(i,j,2) = 7
        tmp1          = NSTvfr(i,j,1)
        NSTvfr(i,j,1) = NSTvfr(i,j,2)
        NSTvfr(i,j,2) = tmp1
        IGBP(7)       = NSTvfr(i,j,2) / 100.0       
        IGBP(16)      = NSTvfr(i,j,1) / 100.0 
       else
        IGBP(7)       = NSTvfr(i,j,1) / 100.0       
        IGBP(16)      = NSTvfr(i,j,2) / 100.0           
       end if
 
       END IF

       IF (nvx .eq. 2) THEN
              
       convert(16, 4) = 50.   ! grass low
       convert(16, 0) = 50.   ! barren soil
       igbp_z0(16   ) = 0.022

       convert(15, 4) = 10.    ! grass low
       convert(15, 0) = 90.    ! barren soil
       igbp_z0(15   ) = 0.001       

       if(NST_sh(i,j).le.ELA) then
        NSTveg(i,j,1) = 15 
        NSTveg(i,j,2) = 16        
        NSTvfr(i,j,2) = 100.0 * (1- NST_sh(i,j)/ ELA)
        NSTvfr(i,j,1) = 100.0 -     NSTvfr(i,j,2)
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
         NSTsfr(i,j,k)=NINT(svat_frac (k)*100.)
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

      DO j=2,my-1
      DO i=2,mx-1 
      IF (NSTsol(i,j).eq.3) THEN

       DO k=1,nvx
        NSTsvt(i,j,k)  =  0
        NSTsfr(i,j,k)  =  0
        NSTveg(i,j,k)  =  0
        NSTvfr(i,j,k)  =  0
        NSTlai(i,j,k)  =  0
        NSTglf(i,j,k)  =  0
       ENDDO

        NSTsvt(i,j,nvx)=  0
        NSTsfr(i,j,nvx)=100
        NSTveg(i,j,nvx)= -1
        NSTvfr(i,j,nvx)=100

      END IF 
      END DO 
      END DO

200   continue
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      END SUBROUTINE
