C   +-------------------------------------------------------------------+
C   |  Subroutine SVTpar                         February 2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Interpolation of large-scale wetness fields to nested grid.       |
C   | SVTpar completes the current data sets (surface characteristics   |
C   | and prognostic variables) for the SVAT model.                     |
C   |                                                                   |
C   | INPUT  : - I_time: time for which the data is requested           |
C   | ^^^^^^^^ - HORint: horizontal interp. type  (1= bilin, 3= bicub)  |
C   |          - LSCfil: input LSC data file (path+name)                |
C   |          - SPHgrd: true if spherical coordinates for LSC model    |
C   |          - NST__x, NST__y : NST grid coordinates (lat./long.)     |
C   |          - NST_sh: topography in nested model                     |
C   |          - NSTsol: soil type                                      |
C   |          - NSTtex: soil texture over land                         |
C   |          - NST_st: soil or sea surface temperature                |
C   |          - NSTdst: deep soil temperature                          |
C   |          - NSTsvt: vegetation type (SVATclassification)           |
C   |          - NSTsfr: fraction of vegetation in the grid cell (SVAT) |
C   |          - NST__t: real temperature                               |
C   |          - SVTwet: imposed soil moisture in all layers (%)        |
C   |          - SVTlsc: soil wetness computed from ECMWF fields        |
C   |                                                                   |
C   | OUTPUT : - NST_ts: soil temperature      (  K  )                  |
C   | ^^^^^^^^ - NST_sw: soil water content    ( m/s )                  |
C   |          - NSTglf: green leaf fraction                            |
C   |          - NSTiwf: 0=no water flux, 1=free drainage               |
C   |                                                                   | 
C   +-------------------------------------------------------------------+

      SUBROUTINE SVTpar

 
      IMPLICIT NONE


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'NSTvar.inc'      
      INCLUDE 'SNDvar.inc'
      INCLUDE 'LSCvar.inc'

C +---Local variables
C +   ---------------

      INTEGER i,j,k,l,it,fID,ierror,SOLtex,tmpvfr

      INTEGER pos_Ox(mx,my),pos_Oy(mx,my)

      REAL    empty1(1),SW_dry(0:12),SW_wet(0:12),relSW1,relSW2,relSW3,
     .        SW_max,zero,unun,aux1,aux2,aux3,totvfr

      REAL    INT_sw(mx,my),INTdsw(mx,my)

      REAL    LSC1Dxs(ni),LSC1Dys(nj),LSC_sws(ni,nj),LSCdsws(ni,nj),
     .        LSC__xs(ni,nj),LSC__ys(ni,nj),LSC_shs(ni,nj)


      CHARACTER*10  var_units
      CHARACTER*100 LSCtit

C +---Data
C +   ----

      DATA zero    /  0.    /
      DATA unun    /  1.    /
 

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Screen message
C +   --------------

      write(6,*) 'Initialisation of soil prognostic variables'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Dry and nearly saturated water contents
C     =======================================


      SW_dry( 0)=1.04e-5 ! 
      SW_dry( 1)=1.38e-2 ! 
      SW_dry( 2)=1.72e-2 ! values in agreement with those of SVAT
      SW_dry( 3)=3.07e-2 ! 
      SW_dry( 4)=5.32e-2 !
      SW_dry( 5)=4.68e-2 ! 
      SW_dry( 6)=7.08e-2 ! 
      SW_dry( 7)=9.50e-2 ! 
      SW_dry( 8)=0.1173  ! 
      SW_dry( 9)=0.1180  ! 
      SW_dry(10)=0.1526  ! 
      SW_dry(11)=0.1628  ! 
      SW_dry(12)=0.0     ! 

      SW_wet( 0)=1.000   ! 
      SW_wet( 1)=0.395   ! 
      SW_wet( 2)=0.410   !  values in agreement with those of SVAT
      SW_wet( 3)=0.435   ! 
      SW_wet( 4)=0.485   ! 
      SW_wet( 5)=0.451   ! 
      SW_wet( 6)=0.420   ! 
      SW_wet( 7)=0.477   !
      SW_wet( 8)=0.476   ! 
      SW_wet( 9)=0.426   ! 
      SW_wet(10)=0.492   !  
      SW_wet(11)=0.482   !   
      SW_wet(12)=0.001   !  


      SW_max   =0.032    ! Soil water content corresponding 
                         ! to saturation in ERA-15

      IF(LSCmod.eq.'E40'.or.LSCmod.eq.'ECM')      
     .SW_max   =0.47     ! Soil water content corresponding
                         ! to saturation in ERA-40

C http://www.ecmwf.int/products/data/technical/soil/discret_soil_lay.html 

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Completion of surface characteristics data set
C +   ==============================================


      DO j=1,my
      DO i=1,mx

       IF (NSTsol(i,j).ge.4) THEN

        NSTiwf(i,j)=1

        DO k=1,nvx
         NSTglf(i,j,k)=1.
        ENDDO
	
C + ... Check fractions of vegetation

        totvfr=0
        DO l=1,nvx
         totvfr=totvfr+NSTsfr(i,j,l)
        ENDDO
        IF (totvfr.ne.100) THEN
         totvfr=totvfr-NSTsfr(i,j,nvx)
         IF (totvfr.ne.0) THEN
          DO l=2,nvx-1
           aux1         =REAL(NSTsfr(i,j,l))
           aux2         =REAL(totvfr)
           aux3         =REAL(NSTsfr(i,j,nvx))
           NSTsfr(i,j,l)=aux1/aux2*(100.-aux3)
          ENDDO
          tmpvfr=0
          DO l=2,nvx
           tmpvfr=tmpvfr+NSTsfr(i,j,l)
          ENDDO
          NSTsfr(i,j,1) = 100 - tmpvfr
         ELSE
          DO l=1,nvx
           NSTsfr(i,j,l) =0
          ENDDO
          NSTsfr(i,j,nvx)=100
         ENDIF
        ENDIF

       ELSE

        NSTiwf(i,j)=0

        DO k=1,nvx
         NSTglf(i,j,k)=0.
        ENDDO

        DO k=1,nvx
         NSTsfr(i,j,k)=0
         NSTsvt(i,j,k)=0
        ENDDO
        NSTsfr(i,j,1) = 100

       ENDIF

      ENDDO
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      IF (SNDing) THEN

       DO l=1,nsl
       DO k=1,nvx
       DO j=1,my
       DO i=1,mx
        NST_ts(i,j,k,l)=NST__t(1,1,mz)-0.006*NST_sh(i,j)
       ENDDO
       ENDDO
       ENDDO
       ENDDO

      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      IF (.not.SNDing.and.SVTlsc) THEN


C +---Open NetCDF file containing LSC data
C +   ====================================

 
       write(6,*) 'Open file  : ',LSCfil
       write(6,*) 'Time step  : ',I_time

C +         *******
       CALL UNropen (LSCfil,fID,LSCtit)
C +         *******


C +---Time for data extraction
C +   ------------------------

       it = I_time


C +---Horizontal coordinates
C +   ----------------------

       DO j=1,my
       DO i=1,mx
        pos_Ox(i,j)=0
        pos_Oy(i,j)=0
       ENDDO
       ENDDO

       IF (REGgrd) THEN

C +          ******
        CALL UNread (fID,'SH'  ,it,1,bi,bj,ni,nj,1,
     .               LSC1Dxs,LSC1Dys,empty1,var_units,LSC_shs)
C +          ****** 

        DO j=1,nj
        DO i=1,ni
         LSC__xs(i,j)=LSC1Dxs(i)
         LSC__ys(i,j)=LSC1Dys(j)
        ENDDO
        ENDDO

       ELSE

C +          ******
        CALL UNread (fID,'lon' ,it,1,bi,bj,ni,nj,1,
     .               LSC1Dxs,LSC1Dys,empty1,var_units,LSC__xs)
C +          ******
        CALL UNread (fID,'lat' ,it,1,bi,bj,ni,nj,1,
     .              LSC1Dxs,LSC1Dys,empty1,var_units,LSC__ys)
C +          ****** 

       ENDIF


C +---Soil wetness
C +   ------------

       IF (LSCmod.ne.'E40'.and.LSCmod.ne.'ECM') THEN

        write(6,'(A,$)') ' 2-D fields : SWL1'

C +         ******
        CALL UNread (fID,'SWL1' ,it,1,bi,bj,ni,nj,1,
     .               LSC1Dxs,LSC1Dys,empty1,var_units,LSC_sws)
C +         ****** 
       ELSE

        write(6,'(A,$)') ' 2-D fields : SWVL1'

C +         ******
        CALL UNread (fID,'SWVL1' ,it,1,bi,bj,ni,nj,1,
     .               LSC1Dxs,LSC1Dys,empty1,var_units,LSC_sws)


       ENDIF

C +         ******
       CALL INThor (HORint,LSC__xs,LSC__ys,LSC_sws,
     .              SPHgrd,NST__x,NST__y,INT_sw,
     .              REGgrd,pos_Ox,pos_Oy)
C +         ******


C +---Deep soil wetness
C +   -----------------

       IF (LSCmod.ne.'E40'.and.LSCmod.ne.'ECM') THEN

        write(6,'(A,$)') ' - SWL2'
        write(6,*)

C +          ******
        CALL UNread (fID,'SWL2' ,it,1,bi,bj,ni,nj,1,
     .               LSC1Dxs,LSC1Dys,empty1,var_units,LSCdsws)
C +          ****** 
       ELSE

        write(6,'(A,$)') ' - SWVL2'
        write(6,*)

C +          ******
        CALL UNread (fID,'SWVL2' ,it,1,bi,bj,ni,nj,1,
     .               LSC1Dxs,LSC1Dys,empty1,var_units,LSCdsws)
C +          ******

       ENDIF

C +         ******
       CALL INThor (HORint,LSC__xs,LSC__ys,LSCdsws,
     .              SPHgrd,NST__x,NST__y,INTdsw,
     .              REGgrd,pos_Ox,pos_Oy)
C +         ******


C +---Close the NetCDF file
C +   =====================

C +         ******
       CALL NCCLOS (fID,ierror)
C +         ******


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ELSE

       write(6,*) 'Imposed soil wetness'

      ENDIF  ! (.not.SNDing.and.SVTlsc)

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Soil temperature
C +   ================

      IF (.not.SNDing) THEN

       DO l=1,nvx
       DO j=1,my
       DO i=1,mx
        IF (NSTsol(i,j).ge.4) THEN
         DO k=1,nsl
          NST_ts(i,j,l,k)= NSTdst(i,j)
         ENDDO
         NST_ts(i,j,l,1)= NST_st(i,j)
         NST_ts(i,j,l,2)=(NST_st(i,j)+NSTdst(i,j))*0.5
        ELSE
         DO k=1,nsl
          NST_ts(i,j,l,k)= NST_st(i,j)
         ENDDO
        ENDIF
       ENDDO
       ENDDO
       ENDDO

      ENDIF

 
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Soil water content
C +   ==================

      DO l=1,nvx
      DO j=1,my
      DO i=1,mx

       IF (.not.SNDing.and.SVTlsc) THEN
        relSW1=INT_sw(i,j)/SW_max
        relSW3=INTdsw(i,j)/SW_max
        relSW2=0.5*(relSW1+relSW3)
       ELSE
        relSW1=SVTwet/100.
        relSW2=SVTwet/100.
        relSW3=SVTwet/100.
       ENDIF

       relSW1=MAX(zero,relSW1)
       relSW2=MAX(zero,relSW2)
       relSW3=MAX(zero,relSW3)
       relSW1=MIN(unun,relSW1)
       relSW2=MIN(unun,relSW2)
       relSW3=MIN(unun,relSW3)

       IF (NSTsol(i,j).ge.4) THEN

        SOLtex=NSTtex(i,j)

        IF (SOLtex.eq.0) SOLtex=2

        DO k=1,nsl-3
        NST_sw(i,j,l,    k)=SW_dry(SOLtex)
     .                     +relSW1*(SW_wet(SOLtex)-SW_dry(SOLtex))
        ENDDO
        NST_sw(i,j,l,nsl-2)=SW_dry(SOLtex)
     .                     +relSW2*(SW_wet(SOLtex)-SW_dry(SOLtex))
        DO k=nsl-1,nsl
        NST_sw(i,j,l,    k)=SW_dry(SOLtex)
     .                     +relSW3*(SW_wet(SOLtex)-SW_dry(SOLtex))
        ENDDO
	
       ELSE

        DO k=1,nsl
         NST_sw(i,j,l,k)=1.
        ENDDO

       ENDIF

      ENDDO
      ENDDO
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      write(6,*) ' '


      RETURN
      END
