C   +-------------------------------------------------------------------+
C   |  Subroutine NSTint                            20/09/2012  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Interpolation of large-scale data to nested grid.                 |
C   |                                                                   |
C   | INPUT  : I_time: time for which the data is requested             |
C   | ^^^^^^^^ HORint: horizontal interp. type  (1= bilin, 3= bicub)    |
C   |          VERint: vertical   interp. type  (1= lin. , 3= cubic)    |
C   |          NST_sh: topography in nested model                       |
C   |          NST_dx: horizontal resolution                            |
C   |          LSCfil: input LSC data file (path+name)                  |
C   |          NSTmod: nested model (used for the vertical grid)        |
C   |          TOPcstLSC,TOPdomLSC,TOPfilt : LSC corrected topography   |
C   |          SPHgrd: true if spherical coordinates for LSC model      |
C   |          CLDcor: true if parameterized cloud water at boundaries  |
C   |                                                                   |
C   | INPUT FILE: a LSRD - Large Scale Raw Data file (NetCDF)           |
C   | ^^^^^^^^^^^                                                       |
C   |   *file name = {LSCmodel}.YY.MM.DD.nC                             |
C   |                 where YYMMDD = Year, Month, and Day of file begin |
C   |   *'time' variable = Universal time from 0 hour the YYMMDD day    |
C   |                      (unit = DAYS)                                |
C   |   *file contents:                                                 |
C   |                      - - - - - - - + - - + - - - + - - - - -      |
C   |                                      variable    |Unit            |
C   |                                   in atm.| 10m   |                | 
C   |                      - - - - - - - + - - + - - - + - - - - -      |
C   |                      Wind          |U    |U10    |m/s             |
C   |                        "           |V    |V10    |m/s             |
C   |                      Specif. humid.|Q    |Q10    |Kg/Kg           |
C   |                      Temperature   |T    |T10    |K               |
C   |                                    |     |       |                |
C   |                      Pressure      |     |SP     |hPa             |
C   |                      Surf. height  |-    |SH     |m               |
C   |                      - - - - - - - + - - + - - - + - - - - -      |
C   |                                                                   | 
C   | OUTPUT : NST__u: U-wind                       ( m/s )             |
C   | ^^^^^^^^ NST__v: V-wind                       ( m/s )             |
C   |          NST__t: real      temperature        (  K  )             |
C   |          NST_pt: potential temperature        (  K  )             |
C   |          NST_qv: specific humidity            (kg/kg)             |
C   |          NST_sp: surface pressure             ( kPa )             |
C   |          NSTsic: sea-ice fraction             (  -  )             |
C   |          NSTsst: sea surface temperature      (  K  )             |
C   |          NST_st: surface temperature          (  K  )             |
C   |          NSTdst: soil    temperature          (  K  )             |
C   |          NST__p: pressure at each lev.        ( kPa )             |
C   |          NST_zz: levels height                (  m  )             |
C   |          NSTgdz: sigma coordinate                                 |
C   |          NSTsol: soil types (ice taken into account)              |
C   |          NSTtke: turbulent kinetic energy     (m2/s2)             |
C   |          NSTuts: surface heat flux            (K.m/s)             |
C   |          NST_qt: total cloud water            (kg/kg)             |
C   |                                                                   | 
C   +-------------------------------------------------------------------+

      SUBROUTINE NSTint

 
      IMPLICIT NONE


C +---Include files
C +   -------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'CTRvar.inc'
      INCLUDE 'LSCvar.inc'
      INCLUDE 'INTvar.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'MARvar.inc'


C +---Local variables
C +   ---------------

      INTEGER i,j,k,it,fID,ierror,zero,kl,mmx,mmy

      INTEGER pos_Ox(mx,my),pos_Oy(mx,my)

      REAL    exxpo,gamTz,empty1(1),cap,ra,grav,humrel,erzmax,
     .        qsat,cp,h2olv,thr_rh,qwater,qsuppl,qclmax,nul,
     .        auxz,auxu,auxv,auxt,auxq,auxe, getpkt, fcort,
     .        lwblon,upblon,lwblat,upblat,lwb_SH,upb_SH,
     .        lwb_SP,upb_SP,lwbST1,upbST1,lwbST2,upbST2,
     .        lwbSW1,upbSW1,lwb10U,upb10U,lwb10V,upb10V,
     .        lwbTCC,upbTCC,lwb__U,upb__U,lwb__V,upb__V,
     .        lwb__T,upb__T,lwb__Q,upb__Q,lwbTKE,upbTKE,
     .        lwbUTS,upbUTS,lwbSIC,upbSIC,lwbSST,upbSST

      REAL    INtmp1(mx,my),INtmp2(mx,my),INtmp3(mx,my),
     .        NSTpk6(mx,my),NSTpx1(mx,my),NSTlp1(mx,my),
     .        LSC_z6(ni,nj),INT_z6(mx,my),NST_z6(mx,my),
     .        LSCpk1(ni,nj),LSCpx1(ni,nj),LSClp1(ni,nj),
     .        LSC1sp,LSC1_t(nk),LSC1_p(nk+1),LSC1hp(nk+1),LSC1sh,
     .        NST1sp,NST1_t(mz),NST1_p(mz+1),NST1hp(mz+1),NST1sh,
     .        INT1sp,INT1sh,
     .        qv_sat(mz),rhoair(mz),deltaz(mz),EQtemp(mz),
     .        WK1Dq(nk+1),WK1Du(nk+1),WK1Dv(nk+1),WK1Dt(nk+1),
     .        WK1_1D(nk),WK2_1D(nk),WK3_1D(nk),qv_max(mz),
     .        WK1m1D(mz),WK2m1D(mz),WK3m1D(mz),qcloud(mz),
     .        WK1De(nk+1),WK1Dh(nk+1),correction
     
      LOGICAL CORsat,LSCiZp(mx,my), NSTiZp(mx,my), iZterm
      LOGICAL decreaseSIC

      CHARACTER*3   emptyC
      CHARACTER*7   nam_SH,nam_SP,namST1,namST2,namSW1,nam10U,nam10V,
     .              nam__U,nam__V,nam__T,nam__Q,namTCC,namlon,namlat,
     .              namTKE,namUTS,nam_QW,nam_QR,nam_QI,nam_QS,namSIC,
     .              namSST
      CHARACTER*10  var_units
      CHARACTER*100 LSCtit
      INTEGER icheck, ipchk, jpchk
      INTEGER im1,ip1,jm1,jp1
      
 
C +---Physical constants
C +   ------------------

      DATA ra    /  287.     d0/
      DATA cp    / 1004.     d0/
      DATA h2olv /    2.5000d+6/
      DATA cap   /    0.28586d0/
      DATA grav  /    9.81   d0/
      DATA emptyC/   '   '    /
      DATA zero  /    0        /
      DATA nul   /    0.       /

      getpkt= exp(-cap*log(100.))
C +...     getpkt: 1. / (100. (kPa) ** cap)


C +---Debug verbose level (0=silent - 3=flood):
C +   -------------------
      icheck=2
C +---Horizontal point for extended check:
      ipchk =39
      jpchk =22


C +---Initialisation
C +   --------------

      lfirst_LSC = .true.
      lfirst_NST = .true.

      IF (NSTmod.eq.'GRA') TOPdomLSC = .true.

      mmx = mx
      mmy = my


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Open NetCDF file containing LSC data
C +   ====================================

C +        *******
      CALL UNropen (LSCfil,fID,LSCtit)
C +        *******


C +---Time for data extraction
C +   ------------------------

      it = I_time


C +---Screen message
C +   --------------

      write(6,*)  'Horizontal and vertical interpolations'
      write(6,*)  '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF (NSTmod.eq.'GRA') THEN
       write(6,*) 'Output for GRADS : imposed LSC topography'
      ENDIF
      write(6,*)  'Open file  : ',LSCfil
      write(6,*)  'Time step  : ',I_time


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Initialisation of temporary variables
C +   =====================================

      DO k=1,nk
       WK1_1D(k) = 0.0
       WK2_1D(k) = 0.0
       WK3_1D(k) = 0.0
      ENDDO

      DO k=1,mz
       WK1m1D(k) = 0.0
       WK2m1D(k) = 0.0
       WK3m1D(k) = 0.0
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Specification of valid interval for LSC values
C +   ==============================================

      lwblon =   -400.0
      upblon =    400.0
      lwblat =   -100.0
      upblat =    100.0
      lwb_SH =  -1000.0
      upb_SH =  10000.0
      lwb_SP =  10000.0
      upb_SP = 130000.0
      lwbSIC =      0.0
      upbSIC =     10.e20
      lwbSST =    100.0
      upbSST =     10.e20
      lwbST1 =    100.0
      upbST1 =    370.0
      lwbST2 =    100.0
      upbST2 =    370.0
      lwbSW1 =     -0.1
      upbSW1 =      1.0
      lwb10U =   -100.0
      upb10U =    100.0
      lwb10V =   -100.0
      upb10V =    100.0
      lwbTCC =     -0.1
      upbTCC =      1.1
      lwb__U =   -300.0
      upb__U =    300.0
      lwb__V =   -300.0
      upb__V =    300.0
      lwb__T =    100.0
      upb__T =    370.0
      lwb__Q =    -0.01
      upb__Q =      1.0
      lwbTKE =     -0.1
      upbTKE =  10000.0
      lwbUTS =   -100.0
      upbUTS =    100.0

      IF (LSCmod.eq.'MAR') THEN
       lwb_SP =  10.0
       upb_SP = 130.0
      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Define variable names
C +   =====================

      IF (LSCmod.eq.'MAR') THEN
       namlon='lon'
       namlat='lat'
       nam_SH='sh'
       nam_SP='pstar'
       namST1='tairSL'
       namST2='tairSL'
       namSW1='-'   
       nam10U='-'  
       nam10V='-'  
       namTCC='-'  
       nam__U='uairDY'
       nam__V='vairDY'
       nam__T='tairDY'
       nam__Q='qvDY'
       namTKE='ect_TE'
       namUTS='SLutsl'
       nam_QW='qwHY'
       nam_QR='qrHY'
       nam_QS='qsHY'
       nam_QI='qiHY'
      ENDIF
      IF (LSCmod.eq.'LMD') THEN
       namlon='nav_lon'
       namlat='nav_lat'
       nam_SH='phis'
       nam_SP='psol'
       namST1='tsol'
       namST2='tsol'
       namSW1='-'
       nam10U='-'
       nam10V='-'
       namTCC='-'  
       nam__U='vitu'
       nam__V='vitv'
       nam__T='temp'
       nam__Q='ovap'
       namTKE='-'
       namUTS='-'
       nam_QW='-'
       nam_QR='-'
       nam_QS='-'
       nam_QI='-'
      ENDIF 
      IF (LSCmod.eq.'LMz') THEN
       namlon='-'
       namlat='-'
       nam_SH='SH'
       nam_SP='SP'
       namST1='SST'  ! Temporary solution !
       namST2='STL1' ! Temporary solution !
       namSW1='-'    ! changed
       nam10U='-'    ! changed + NB: srfT may be useful ?
       nam10V='-'
       namTCC='-' 
       nam__U='U'
       nam__V='V'
       nam__T='T'
       nam__Q='Q'
       namTKE='-'
       namUTS='-'
       nam_QW='-'    ! not in use yet
       nam_QR='-'
       nam_QS='-'
       nam_QI='-'
      ENDIF 
      IF (LSCmod.eq.'ECM'.or.LSCmod.eq.'E15'.or.LSCmod.eq.'E20'
     ..or.LSCmod.eq.'E40'.or.LSCmod.eq.'EIN') THEN
       LSCmod="E40"
       namlon='lon'
       namlat='lat'
       nam_SH='SH'
       nam_SP='SP'
       namSST='SSTK'
       namSIC='CI'
       namST1='STL1'
       namST2='STL2'
       namSW1='SWVL1'
       if (LSCmod.eq.'E15') namSW1='SWL1'
       nam10U='-'
       nam10V='-'
       namTCC='-'  
       nam__U='U'
       nam__V='V'
       nam__T='T'
       nam__Q='Q'
       namTKE='-'
       namUTS='-'
       nam_QW='CLWC'
       nam_QR='-'
       nam_QS='-'
       nam_QI='CIWC'
      ENDIF
      IF (LSCmod.eq.'CM3'.or.LSCmod.eq.'EM5'.or.LSCmod.eq.'CAN'.or.
     .    LSCmod.eq.'NOR'.or.LSCmod.eq.'CSI'.or.LSCmod.eq.'BCC'.or.
     .    LSCmod.eq.'MIR'.or.LSCmod.eq.'CM5'.or.LSCmod.eq.'AC3') THEN
       namlon='lon'
       namlat='lat'
       nam_SH='SH'
       nam_SP='SP'
       namSST='SST2'
       namSIC='CI'
       namST1='SST1'
       namST2='SST1'
       namSW1='-'
       nam10U='-'
       nam10V='-'
       namTCC='-'  
       nam__U='U'
       nam__V='V'
       nam__T='T'
       nam__Q='Q'
       namTKE='-'
       namUTS='-'
       nam_QW='-'
       nam_QR='-'
       nam_QS='-'
       nam_QI='-'
       LSCmod='GCM'
      ENDIF
      IF (LSCmod.eq.'20C'.or.LSCmod.eq.'NCP'.or.
     .    LSCmod.eq.'NC1'.or.LSCmod.eq.'NC2') THEN
       namlon='lon'
       namlat='lat'
       nam_SH='SH'
       nam_SP='SP'
       namSST='SST2'
       namSIC='CI'
       namST1='SST1'
       namST2='SST1'
       namSW1='-'
       nam10U='-'
       nam10V='-'
       namTCC='-'  
       nam__U='U'
       nam__V='V'
       nam__T='T'
       nam__Q='Q'
       namTKE='-'
       namUTS='-'
       nam_QW='-'
       nam_QR='-'
       nam_QS='-'
       nam_QI='-'
       LSCmod='NCP'
      ENDIF
C http://www.ecmwf.int/products/data/technical/GRIB_tables/table_128.html

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Horizontal coordinates
C +   ======================


      DO j=1,my
      DO i=1,mx
       pos_Ox(i,j)=0
       pos_Oy(i,j)=0
      ENDDO
      ENDDO

      IF (REGgrd) THEN

C +         ******
       CALL UNread (fID,nam_SH,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_sh)
C +         ****** 
       CALL VALchk (nam_SH,ni,nj,LSC_sh,lwb_SH,upb_SH)
C +         ******

       DO j=1,nj
       DO i=1,ni
        LSC__x(i,j)=LSC1Dx(i)
        LSC__y(i,j)=LSC1Dy(j)
       ENDDO
       ENDDO

       IF (LSCmod.eq.'LMD') THEN
        DO j=1,nj
        DO i=1,ni
         LSC_sh(i,j)=LSC_sh(i,j)/9.81
        ENDDO
        ENDDO
       ENDIF
C +    This is strange since 'LMD' sets REGgrd = false ?
C +    (remember this is _not_ the Jan2002 LMD-Z, which is named LMz)

      ELSE

C +         ******
       CALL UNread (fID,namlon,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__x)
C +         ******
       CALL VALchk (namlon,ni,nj,LSC__x,lwblon,upblon)
C +         ******
       CALL UNread (fID,namlat,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__y)
C +         ******
       CALL VALchk (namlat,ni,nj,LSC__y,lwblat,upblat)
C +         ******

      ENDIF
      
      IF (icheck.GE.3) THEN
         write(*,*) 'NSTint: input coordinates:'
         write(*,*) (LSC__x(i,1),i=1,ni)
         write(*,*) (LSC__y(1,j),j=1,nj)
         IF (REGgrd) write(*,*) 'Grid is assumed rectangular'
      ENDIF

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Initialize vertical grid
C +   ========================

      DO k=1,nk
       LSCgdz(k)=0.
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Topography
C +   ==========


C +        ******
      CALL UNread (fID,nam_SH,it,1,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSC_sh)
C +        ****** 
      CALL VALchk (nam_SH,ni,nj,LSC_sh,lwb_SH,upb_SH)
C +        ******

C +        ******
      CALL INThor (HORint,LSC__x,LSC__y,LSC_sh,
     .             SPHgrd,NST__x,NST__y,INT_sh,
     .             REGgrd,pos_Ox,pos_Oy)
C +        ******


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CORRECTION OF PRESCRIBED TOPOGRAPHY WITH LSC TOPOGRAPHY
C +   =======================================================


C +---Imposed LSC topography in the relaxation zone
C +   ---------------------------------------------

      IF (TOPcstLSC) THEN

       TOPopt=2
C +         ******
       CALL TOPcor (TOPopt,NST__x,NST__y,NST_sh,INT_sh)
C +         ******

      ENDIF


C +---Imposed LSC topography in the whole domain
C +   ------------------------------------------

      IF (TOPdomLSC) THEN

       TOPopt=3
C +         ******
       CALL TOPcor (TOPopt,NST__x,NST__y,NST_sh,INT_sh)
C +         ******

      ENDIF


C +---Topography filtering (2D and 3D)
C +   --------------------------------

      IF (TOPdomLSC.and.TOPfilt) THEN

       TOPopt=5
C +         ******
       CALL TOPcor (TOPopt,NST__x,NST__y,NST_sh,INT_sh)
C +         ******

      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Horizontal interpolation of surface fields
C +   ==========================================


C +---Surface Pressure
C +   ----------------

      WRITE(6,'(A,$)') ' 2-D fields : '//nam_SH//'- '//nam_SP

C +        ******
      CALL UNread (fID,nam_SP,it,1,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSC_sp)
C +        ****** 
      CALL VALchk (nam_SP,ni,nj,LSC_sp,lwb_SP,upb_SP)
C +        ******

      IF (LSCmod.ne.'MAR') THEN
C +         ******
       CALL LSuCHG (LSC_sp,1.E-3) !(Change units: Pa-->kPa)
C +         ******
      ENDIF

C +        ******
      CALL INThor (HORint,LSC__x,LSC__y,LSC_sp,
     .             SPHgrd,NST__x,NST__y,INT_sp,
     .             REGgrd,pos_Ox,pos_Oy)
C +        ******


C +---Sea-Ice Fraction
C +   ----------------

      IF (LSCmod.eq.'E40'.or.LSCmod.eq.'ECM'.or.LSCmod.eq.'GCM'
     .                                      .or.LSCmod.eq.'NCP') THEN

        WRITE(6,'(A,$)') '- '//namSIC

C +          ******
        CALL UNread (fID,namSIC,it,1,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSCsic)
C +          ****** 
        CALL VALchk (namSIC,ni,nj,LSCsic,lwbSIC,upbSIC)
C +          ******

C +          ******
        CALL INTmsk (LSCsic)
C +          ******

C +          ******
        CALL INThor (HORint,LSC__x,LSC__y,LSCsic,
     .               SPHgrd,NST__x,NST__y,INTsic,
     .               REGgrd,pos_Ox,pos_Oy)
C +          ******

        DO i=1,mx
        DO j=1,my
          IF (NSTsol(i,j).GE.3)    INTsic(i,j) = 0.
          INTsic(i,j) = max(0.,min(INTsic(i,j),1.0))
        ENDDO
        ENDDO

      END IF


C +---Soil or Sea surface temperature
C +   -------------------------------

      WRITE(6,'(A,$)') '- '//namST1

C +        ******
      CALL UNread (fID,namST1,it,1,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSC_st)
C +        ****** 
      CALL VALchk (namST1,ni,nj,LSC_st,lwbST1,upbST1)
C +        ******

C +        ******
      CALL INThor (HORint,LSC__x,LSC__y,LSC_st,
     .             SPHgrd,NST__x,NST__y,INT_st,
     .             REGgrd,pos_Ox,pos_Oy)
C +        ******

C +---Sea surface temperature
C +   -----------------------


      IF (LSCmod.eq.'E40'.or.LSCmod.eq.'ECM'.or.LSCmod.eq.'GCM'
     .                                      .or.LSCmod.eq.'NCP') THEN

        WRITE(6,'(A,$)') '- '//namSST

C +          ******
        CALL UNread (fID,namSST,it,1,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSCsst)
C +          ****** 
        CALL VALchk (namSST,ni,nj,LSCsst,lwbSST,upbSST)
C +          ******

        IF(LSCmod.eq.'E40'.or.LSCmod.eq.'ECM') THEN
         DO i = 1,ni
         DO j = 1,nj
          if(LSCsst(i,j)>=250..and.LSCsst(i,j)<=350.) then
           !+ ---------------------------------------------------------- +!
           !+ http://www.ecmwf.int/research/ifsdocs/CY28r1/Assimilation/
           !+ Assimilation-14-4.html
           !+ For grid boxes characterized by sea-ice concentrations
           !+ exceeding 20% the SST is set to -1.7 degC.
           !+ ---------------------------------------------------------- +!
           LSCsst(i,j)=LSC_st(i,j)
          endif
         ENDDO
         ENDDO
        ENDIF

C +          ******
        CALL INTmsk (LSCsst)
C +          ******

C +          ******
        CALL INThor (HORint,LSC__x,LSC__y,LSCsst,
     .               SPHgrd,NST__x,NST__y,INTsst,
     .               REGgrd,pos_Ox,pos_Oy)
C +          ******

        DO j = 1,my
        DO i = 1,mx
 
         IF (NSTsol(i,j).le.2) THEN
 
          INT_st(i,j)=INTsst(i,j)  
 
         ENDIF
 
        ENDDO
        ENDDO

      END IF

C +---Soil or Sea temperature
C +   -----------------------

      WRITE(6,'(A,$)') '- '//namST2

C +        ******
      CALL UNread (fID,namST2,it,1,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSCdst)
C +        ****** 
      CALL VALchk (namST2,ni,nj,LSCdst,lwbST2,upbST2)
C +        ******

C +        ******
      CALL INThor (HORint,LSC__x,LSC__y,LSCdst,
     .             SPHgrd,NST__x,NST__y,INTdst,
     .             REGgrd,pos_Ox,pos_Oy)
C +        ******


C +---Total Cloud Cover
C +   -----------------

      IF (namTCC.ne.'-') THEN

       WRITE(6,'(A,$)') '- '//namTCC

C +         ******
       CALL UNread (fID,namTCC,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSCtcc)
C +	    ****** 
       CALL VALchk (namTCC,ni,nj,LSCtcc,lwbTCC,upbTCC)
C +	    ******

C +	    ******
       CALL INThor (HORint,LSC__x,LSC__y,LSCtcc,
     .              SPHgrd,NST__x,NST__y,INTtcc,
     .              REGgrd,pos_Ox,pos_Oy)
C +	    ******

      ENDIF

C +---Temperature difference between 1st atm. level and soil/sea
C +   ----------------------------------------------------------

      WRITE(6,'(A,$)') '- '//nam__T

C +        ******
      CALL UNread (fID,nam__T,it,nk,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSC__t)
C +        ****** 
      CALL VALchk (nam__T,ni,nj,LSC__t,lwb__T,upb__T)
C +        ******

      DO j=1,nj
      DO i=1,ni
       LSC_dt(i,j)=LSC__t(i,j)-LSC_st(i,j)
      ENDDO
      ENDDO

C +        ******
      CALL INThor (HORint,LSC__x,LSC__y,LSC_dt,
     .             SPHgrd,NST__x,NST__y,INT_dt,
     .             REGgrd,pos_Ox,pos_Oy)
C +        ******

      DO j=1,nj
      DO i=1,ni
       LSC_pt(i,j)=LSC__t(i,j)*(100./LSC_sp(i,j))**cap 
      ENDDO
      ENDDO

C +        ******
      CALL INThor (HORint,LSC__x,LSC__y,LSC_pt,
     .             SPHgrd,NST__x,NST__y,INtmp1,
     .             REGgrd,pos_Ox,pos_Oy)
C +        ******

C +        ******
      CALL PUT2D3 (INtmp1,nk,INT_pt)
C +        ******


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Correction of surface pressure according to topography change
C +   =============================================================
 
C +...Computation of a surface pressure adapted to MAR topography,
C +...using 2 simple assumptions : - constant T gradient = gamTz
C +...                             - basic T = 1st level near surface

C +---Constants
C +   ---------

      gamTz = - 6.5E-3
      exxpo = - grav / (gamTz * ra)


C +---Compute surface pressure according to topography changes
C +   --------------------------------------------------------

      DO j = 1,my
      DO i = 1,mx
        NST_sp(i,j)= INT_sp(i,j) 
     .  * (1.+gamTz*(NST_sh(i,j)-INT_sh(i,j))/INT_st(i,j))**exxpo
      END DO
      END DO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



C +---Correction of Sea-Ice Fraction
C +   ==============================

      IF (LSCmod.eq.'E40'.or.LSCmod.eq.'ECM'.or.LSCmod.eq.'GCM'
     .                                      .or.LSCmod.eq.'NCP')THEN
        DO j = 1,my
        DO i = 1,mx
          NSTsic(i,j)= INTsic(i,j)
          NSTsst(i,j)= INTsst(i,j)
        END DO
        END DO
      ELSE  
        DO j = 1,my
        DO i = 1,mx
          NSTsic(i,j)= -1.           ! Feeds MAR with nitroglycerine
          NSTsst(i,j)= -99.9
        END DO
        END DO
      END IF



C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



C +---Correction of soil temperature according to topography change
C +   =============================================================

      DO j = 1,my
      DO i = 1,mx
      
       IF (NST_sh(i,j).gt.1..and.NSTsol(i,j).ge.3) THEN
       
        NST_st(i,j)=INT_pt(i,j,nk)/(100./NST_sp(i,j))**cap
     .             -INT_dt(i,j)
     
C +...  Temperature diff. between 1st level and surface is conserved

       ELSE
       
        NST_st(i,j)=INT_st(i,j)

C +...  No correction for the sea surface temperature
C       Possible correction for SST-Reynolds in SSTint.f

       ENDIF

C +...  No correction for the sea surface temperature

       fcort = gamTz * (NST_sh(i,j)-INT_sh(i,j))
       NSTdst(i,j)= INTdst(i,j) + fcort
       
      ENDDO
      ENDDO
 
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Atmospheric variables at surface:
C +   (bottom boundary for vertic. interpolation)
C +   ===========================================       


C +---10-m U-wind
C +   -----------

      IF (nam10U.ne.'-') THEN
C +   (if 10m wind not available, 0 will be used for interpolation)

       WRITE(6,'(A,$)') '- '//nam10U
C +         ******
       CALL UNread (fID,nam10U,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__u)
C +         ****** 
       CALL VALchk (nam10U,ni,nj,LSC__u,lwb__U,upb__U)
C +         ******

C +         ******
       CALL INThor (HORint,LSC__x,LSC__y,LSC__u,
     .             SPHgrd,NST__x,NST__y,INtmp2,
     .             REGgrd,pos_Ox,pos_Oy)
C +         ******

      ENDIF


C +---10-m V-wind
C +   -----------

      IF (nam10V.NE.'-') THEN
C +   (if 10m wind not available, 0 will be used for interpolation)

       WRITE(6,'(A,$)') '- '//nam10V
C +         ******
       CALL UNread (fID,nam10V,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__v)
C +         ****** 
       CALL VALchk (nam10V,ni,nj,LSC__v,lwb__V,upb__V)
C +         ******

C +         ******
       CALL INThor (HORint,LSC__x,LSC__y,LSC__v,
     .             SPHgrd,NST__x,NST__y,INtmp3,
     .             REGgrd,pos_Ox,pos_Oy)
C +         ******
      ENDIF


C +---Wind vector rotation (according to projection)
C +   --------------------

      IF (nam10U.NE.'-'.AND.nam10V.NE.'-') THEN
        IF (NST_dx.gt.0.01.and.NSTmod.ne.'GRA') then
           if (maptyp.ge.1) then
C +                 ******
               CALL VecRot (NST__x,NST__y,NST_dx,INtmp2,INtmp3)
C +                 ******
           else
C              ->Polar Stereographic Projection (Antarctica)
C +                 *****************
               CALL VecRot_StereoSouth (GEddxx,NST__x,INtmp2,INtmp3)
C +                 *****************
           endif
       ENDIF

C +          ******
        CALL PUT2D3 (INtmp2,nk+1,INT__u)
C +          ******

C +          ******
        CALL PUT2D3 (INtmp3,nk+1,INT__v)
C +          ******
      ENDIF

C +---Potential temperature
C +   ---------------------

      IF (namST1.eq.'-') THEN

       WRITE(6,'(A,$)') '- '//nam__T
C +         ******
       CALL UNread (fID,nam__T,it,nk,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__t)
C +         ****** 
       CALL VALchk (nam__T,ni,nj,LSC__t,lwb__T,upb__T)
C +         ******

      ELSE

       WRITE(6,'(A,$)') '- '//namST1
C +         ******
       CALL UNread (fID,namST1,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__t)
C +         ****** 
       CALL VALchk (namST1,ni,nj,LSC__t,lwbST1,upbST1)
C +         ******

      ENDIF

      DO j=1,nj
      DO i=1,ni
       LSC_pt(i,j)=LSC__t(i,j)*(100./LSC_sp(i,j))**cap 
      ENDDO
      ENDDO

C +        ******
      CALL INThor (HORint,LSC__x,LSC__y,LSC_pt,
     .             SPHgrd,NST__x,NST__y,INtmp1,
     .             REGgrd,pos_Ox,pos_Oy)
C +        ******

C +        ******
      CALL PUT2D3 (INtmp1,nk+1,INT_pt)
C +        ******


C +---Water vapour
C +   ------------

      IF (namSW1.eq.'-') THEN

       WRITE(6,*)
       WRITE(6,'(A,$)') '              '//nam__Q
C +         ******
       CALL UNread (fID,nam__Q,it,nk,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_qv)
C +         ****** 
       CALL VALchk (nam__Q,ni,nj,LSC_qv,lwb__Q,upb__Q)
C +         ******

      ELSE

       WRITE(6,*)
       WRITE(6,'(A,$)') '              '//namSW1
C +         ******
       CALL UNread (fID,namSW1,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_qv)
C +         ****** 

       IF (LSCmod.eq.'E40'.or.LSCmod.eq.'ECM') THEN
        do i=1,ni
        do j=1,nj
c        LSC_qv(i,j) = max(0.,LSC_qv(i,j))/0.47
c    .               * qsat(LSC__t(i,j),LSC_sp(i,j))
         LSC_qv(i,j) = LSC_qv(i,j) * 0.07  
        enddo      
        enddo
       ENDIF

c http://www.ecmwf.int/products/data/technical/soil/discret_soil_lay.html

C +         ****** 
       CALL VALchk (namSW1,ni,nj,LSC_qv,lwbSW1,upbSW1)
C +         ******

      ENDIF

C +        ******
      CALL INThor (HORint,LSC__x,LSC__y,LSC_qv,
     .             SPHgrd,NST__x,NST__y,INtmp1,
     .             REGgrd,pos_Ox,pos_Oy)
C +        ******

C +        ******
      CALL PUT2D3 (INtmp1,nk+1,INT_qv)
C +        ******


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Reading and horizontal interpolation (+rotation) 
C +   (for each atm. prognostic variable and each level)
C +   ==================================================

      IF (CLDcor .AND. nam_QW.NE.'-') THEN
         WRITE(6,*)
         WRITE(6,'(A)') ' LSC Cloud water will be added to Qv' 
      ENDIF

      WRITE(6,*)
      WRITE(6,'(A,$)') ' 3-D fields :' 

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO k = nk,1,-1   !*BEGIN LOOP on vertical levels
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       WRITE(6,'(I3,$)') k
       IF (MOD(nk-k+1,20).eq.0) THEN
        WRITE(6,*)
        WRITE(6,'(A,$)') '             '
       ENDIF


C +----U-Wind
C +    ------

C +         ******
       CALL UNread (fID,nam__U,it,k,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__u)
C +         ****** 
       CALL VALchk (nam__U,ni,nj,LSC__u,lwb__U,upb__U)
C +         ******

       IF (REGgrd) THEN
C +      Rectangular grids can be staggered, irregular can not:
C +           ******
         CALL STGgrd (LSC1Dx,LSC1Dy,LSCtmp,LSCtm2,ni,nj)
         CALL INThor (HORint,LSCtmp,LSCtm2,LSC__u,
     .                SPHgrd,NST__x,NST__y,INtmp2,
     .                REGgrd,pos_Ox,pos_Oy)
C +           ******
       ELSE      
C +           ******
         CALL INThor (HORint,LSC__x,LSC__y,LSC__u,
     .                SPHgrd,NST__x,NST__y,INtmp2,
     .                REGgrd,pos_Ox,pos_Oy)
C +           ******
       ENDIF
       

C +----V-Wind
C +    ------
 
C +         ******
       CALL UNread (fID,nam__V,it,k,bi,bj,ni,njv,1,
     .              LSC1Dx,LSC1Vy,empty1,var_units,LSC__v)
C +         ****** 
       CALL VALchk (nam__V,ni,njv,LSC__v,lwb__V,upb__V)
C +         ******

       IF (REGgrd) THEN
C +      Rectangular grids can be staggered, irregular can not 
C +      + staggered grids may have a different size => INThorV
C +           ******
       CALL   STGgrd (LSC1Dx,LSC1Vy,LSCgv1,LSCgv2,ni,njv)
       CALL   INThorV(HORint,LSCgv1,LSCgv2,LSC__v,
     .                SPHgrd,NST__x,NST__y,INtmp3,
     .                  REGgrd,pos_Ox,pos_Oy)
C +           ******
       ELSE
C +           ******
       CALL   INThor (HORint,LSC__x,LSC__y,LSC__v,
     .                SPHgrd,NST__x,NST__y,INtmp3,
     .                REGgrd,pos_Ox,pos_Oy)
C +           ******
       ENDIF

C +----Wind vector rotation (according to projection)
C +    --------------------

       IF (NST_dx.gt.0.01.and.NSTmod.ne.'GRA'
     .    .and.mmx.ne.1.and.mmy.ne.1) then
           if (maptyp.ge.1) then
C +                 ******
               CALL VecRot (NST__x,NST__y,NST_dx,INtmp2,INtmp3)
C +                 ******
           else
C              ->Polar Stereographic Projection (Antarctica)
C +                 *****************
               CALL VecRot_StereoSouth (GEddxx,NST__x,INtmp2,INtmp3)
C +                 *****************
           endif
       ENDIF

C +         ******
       CALL PUT2D3 (INtmp2,k,INT__u)
C +         ******

C +         ******
       CALL PUT2D3 (INtmp3,k,INT__v)
C +         ******

 
C +----Water vapour I : read
C +    ----------------------

C +         ******
       CALL UNread (fID,nam__Q,it,k,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_qv)
C +         ****** 
       CALL VALchk (nam__Q,ni,nj,LSC_qv,lwb__Q,upb__Q)
C +         ******

C +----Add cloud water vapour to Qv -> clouds at boundaries
C +    ----------------------------------------------------
C +    Only if cloud water is available in LSC fields:
C +    Note : Qv is added in the LSC variables because this is
C +    somewhat more consistent for the 600 hPa correction, which
C +    compares LSC and interpolated output fields

       IF (CLDcor) THEN
       
       IF (nam_QW.NE.'-') THEN

C +          ******
        CALL UNread (fID,nam_QW,it,k,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSCtmp)
        CALL VALchk (nam_QW,ni,nj,LSCtmp,lwb__Q,upb__Q)
C +          ******

        DO j=1,nj
        DO i=1,ni
         LSC_qv(i,j) = LSC_qv(i,j) + LSCtmp(i,j)
        ENDDO
        ENDDO
       ENDIF

       IF (nam_QI.NE.'-') THEN

C +          ******
        CALL UNread (fID,nam_QI,it,k,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSCtmp)
        CALL VALchk (nam_QI,ni,nj,LSCtmp,lwb__Q,upb__Q)
C +          ******

        DO j=1,nj
        DO i=1,ni
         LSC_qv(i,j) = LSC_qv(i,j) + LSCtmp(i,j)
        ENDDO
        ENDDO
       ENDIF

       ENDIF
       
C!+CA  WARNING, decreases the specific humidity
C!+CA  LSC_qv = LSC_qv *0.93 !CA WARNING WARNING > NOR

c      LSC_qv = LSC_qv *1.02 ! WARNING NCEPv1 over GRD
c      LSC_qv = LSC_qv *1.02 ! WARNING NCEPv2 over GRD

c      correction=1.08*      (1950.-RUNiyr)**1.5/(1950.-1871.+1.)**1.5
c    .           +0.95 * (1.-(1950.-RUNiyr)**1.5/(1950.-1871.+1.)**1.5)
c      if( RUNiyr >= 1950)          correction=0.95
c      correction=min(1.05,max(0.95,correction))
c      LSC_qv = LSC_qv *correction ! WARNING 20CR over GRD

C +----Water vapour II : interpolate / store
C +    -------------------------------------

C +         ******
       CALL INThor (HORint,LSC__x,LSC__y,LSC_qv,
     .              SPHgrd,NST__x,NST__y,INtmp1,
     .              REGgrd,pos_Ox,pos_Oy)
C +         ******

C +         ******
       CALL PUT2D3 (INtmp1,k,INT_qv)
C +         ******


C +--- Pressure (required here for theta and geop. correction)
C +    --------
       
       IF (LSCmod.NE.'LMz') THEN

         DO j=1,nj
         DO i=1,ni

          LSC1sp=LSC_sp(i,j)
          LSC1sh=LSC_sh(i,j)

C +            ******
          CALL VERgrd (LSCmod,emptyC,fID,k,nk,LSC1sp,
     .                 LSC1sh,LSC1_t,LSC1_p,LSC1hp,
     .                 LSCgdz,WK1_1D,WK2_1D,WK3_1D)
C +            ******
C         LSC1_t ne semble ni affecte ni utilise ! Supprimer ?
C +
C +       NOTE: the code is restructured for models which have
C +       arbitrary-pressure levels (not hybrid, or not exacly: LMDz)
C +       => LSC1hp is not usefull here (maintained for compatibility)
          
           LSC__p(i,j)=LSC1_p(k)

         ENDDO
         ENDDO
 
         
         DO j=1,my
         DO i=1,mx

          INT1sp=INT_sp(i,j)
          INT1sh=INT_sh(i,j)

C +            ******
          CALL VERgrd (LSCmod,emptyC,fID,k,nk,INT1sp,
     .                INT1sh,LSC1_t,INT1Dp,INT1Dz,
     .                LSCgdz,WK1_1D,WK2_1D,WK3_1D)
C +            ******

          INT__p(i,j,k)=INT1Dp(k)

         ENDDO
         ENDDO
                 
       ELSE

C +    ** This is necessary for models which have
C +       arbitrary-pressure levels (not hybrid)
C +       Used for LMDZ, for which mid-levels are not exactly hybrid.
C +       It must be done separately because VERgrd does not
C +       know about horizontal postion==> VERgrd might be rewriten.

C +           ******
         CALL UNread (fID,'P',it,k,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSC__p)
         CALL VALchk ('P',ni,nj,LSC__p,0.0,1.2E5)
C +           ******

         IF (var_units(1:2).EQ.'Pa') THEN
C +             ******
           CALL LSuCHG (LSC__p,1.E-3) !(Change units: Pa-->kPa)
C +             ******
         ELSE
           WRITE(*,*) 'NSTint: unknown p units =',var_units
           STOP
         ENDIF

C +           ******
         CALL INThor (1     ,LSC__x,LSC__y,LSC__p,
     .                SPHgrd,NST__x,NST__y,INtmp1,
     .                REGgrd,pos_Ox,pos_Oy)
C +           ******

C +           ******
         CALL PUT2D3 (INtmp1,k,INT__p)
C +           ******

       ENDIF

C +--- Potential temperature
C +    ---------------------

C +         ******
       CALL UNread (fID,nam__T,it,k,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__t)
C +         ****** 
       CALL VALchk (nam__T,ni,nj,LSC__t,lwb__T,upb__T)
C +         ******

       DO j=1,nj
       DO i=1,ni

        LSC_pt(i,j)=LSC__t(i,j)*(100./LSC__p(i,j))**cap 

       ENDDO
       ENDDO

C +         ******
       CALL INThor (HORint,LSC__x,LSC__y,LSC_pt,
     .              SPHgrd,NST__x,NST__y,INtmp1,
     .              REGgrd,pos_Ox,pos_Oy)
C +         ******

C +         ******
       CALL PUT2D3 (INtmp1,k,INT_pt)
C +         ******

C +--- Geopotential of 600 hPa level (for later use)
C +    ---------------------------------------------
C +    (must be done for each k starting from nk->1)

      IF (CORzz6) THEN
C +         ******
       CALL NSTzz6(LSC_pt, LSC_qv, LSC_sh, LSC_sp, LSC__p, k,ni,nj,nk,
     .             LSCpk1, LSCpx1, LSClp1, LSCiZp, iZterm, LSC_z6)
C +         ******
      ENDIF


C +--- Relative Humidity
C +    -----------------

       DO j=1,nj
       DO i=1,ni

        LSC_rh(i,j)=LSC_qv(i,j)/qsat(LSC__t(i,j),LSC__p(i,j))

       ENDDO
       ENDDO

C +         ******
       CALL INThor (HORint,LSC__x,LSC__y,LSC_rh,
     .              SPHgrd,NST__x,NST__y,INtmp1,
     .              REGgrd,pos_Ox,pos_Oy)
C +         ******

C +         ******
       CALL PUT2D3 (INtmp1,k,INT_rh)
C +         ******


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ENDDO       ! END LOOP ON VERTICAL LEVELS
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(6,*)
      WRITE(6,*)

C +---Interpolate Z600 to MAR grid (for later use)
C +   --------------------------------------------
C +   (must be done after calls to NSTzz6 for all k)

      IF (CORzz6) THEN
C +         ******
       CALL INThor (HORint,LSC__x,LSC__y,LSC_z6,
     .             SPHgrd,NST__x,NST__y,INT_z6,
     .             REGgrd,pos_Ox,pos_Oy)
C +         ******
      ENDIF


C +---Vertical grid in the NST model (depend on SP)
C +   ==============================

      DO j=1,my
      DO i=1,mx

       NST1sp=NST_sp(i,j)
       NST1sh=NST_sh(i,j)

C +         ******
       CALL VERgrd (emptyC,NSTmod,fID,zero,mz,NST1sp,
     .              NST1sh,NST1_t,NST1_p,NST1hp,
     .              NSTgdz,WK1m1D,WK2m1D,WK3m1D)
C +         ******

       DO k=1,mz
        NST_hp(i,j,k)=NST1hp(k)
        NST__p(i,j,k)=NST1_p(k)
       ENDDO

      ENDDO
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Vertical interpolation
C +   ======================

C +---Prepare LSC vertical grid and data
C +   ----------------------------------

      IF (VERint.ne.1.and.VERint.ne.3) THEN
       write(6,*) 'CAUTION :'
       write(6,*) 'Horizontal interpolation order incorrectly'
       write(6,*) 'specified. Default is set to linear'
       VERint=1
      ENDIF

      DO j=1,my
      DO i=1,mx

        DO k=1,nk
         INT1Du(k)=INT__u(i,j,k)
         INT1Dv(k)=INT__v(i,j,k)
         INT1Dt(k)=INT_pt(i,j,k)
         INT1Dq(k)=INT_qv(i,j,k)
         INT1Dp(k)=INT__p(i,j,k)
         INT1Dh(k)=INT_rh(i,j,k)
        ENDDO

      INT1sp=INT_sp(i,j)
      INT1Du(nk+1)=0.
      INT1Dv(nk+1)=0.
      INT1Dt(nk+1)=INT_pt(i,j,nk)
      INT1Dq(nk+1)=INT_qv(i,j,nk)
      INT1Dh(nk+1)=INT_rh(i,j,nk)
C +   WARNING: are you erasing the nk+1 values set above ?

      INT1Dp(nk+1)=INT_sp(i,j)
       
C +        ******
      CALL VERhyb (nk, INT1sp, INT1Dp, INT1Dz)
C +        ******


C +---Linear interpolation (default)
C +   ------------------------------

       IF (VERint.eq.1) THEN

        DO k=1,mz
         auxz = NST_hp(i,j,k)
         CALL INTlin(INT1Dz,INT1Du,nk+1,auxz,auxu)
         CALL INTlin(INT1Dz,INT1Dv,nk+1,auxz,auxv)
         CALL INTlin(INT1Dz,INT1Dt,nk+1,auxz,auxt)
         CALL INTlin(INT1Dz,INT1Dq,nk+1,auxz,auxq)
         NST__u(i,j,k) = auxu
         NST__v(i,j,k) = auxv
         NST_pt(i,j,k) = auxt
         NST_qv(i,j,k) = auxq
         CALL INTlin(INT1Dz,INT1Dh,nk+1,auxz,auxq)
         NST_rh(i,j,k) = auxq
        ENDDO

       ENDIF

 
C +---Natural cubic spline (optional)
C +   -------------------------------

       IF (VERint.eq.3.and.NSTmod.ne.'GRA') THEN

        CALL SPLINE(INT1Dz,INT1Du,nk+1,1.E30,1.E30,WK1Du)
        CALL SPLINE(INT1Dz,INT1Dv,nk+1,1.E30,1.E30,WK1Dv)
        CALL SPLINE(INT1Dz,INT1Dt,nk+1,1.E30,1.E30,WK1Dt)
        CALL SPLINE(INT1Dz,INT1Dq,nk+1,1.E30,1.E30,WK1Dq)
        CALL SPLINE(INT1Dz,INT1Dh,nk+1,1.E30,1.E30,WK1Dh)

        DO k =1,mz
         auxz = NST_hp(i,j,k)
         CALL SPLINT(INT1Dz,INT1Du,WK1Du,nk+1,auxz,auxu)
         CALL SPLINT(INT1Dz,INT1Dv,WK1Dv,nk+1,auxz,auxv)
         CALL SPLINT(INT1Dz,INT1Dt,WK1Dt,nk+1,auxz,auxt)
         CALL SPLINT(INT1Dz,INT1Dq,WK1Dq,nk+1,auxz,auxq)
         NST__u(i,j,k) = auxu
         NST__v(i,j,k) = auxv
         NST_pt(i,j,k) = auxt
         NST_qv(i,j,k) = auxq
         CALL SPLINT(INT1Dz,INT1Dh,WK1Dh,nk+1,auxz,auxq)
         NST_rh(i,j,k) = auxq
        ENDDO

       ENDIF


C +---Linear interpolation (pressure coord.) for GRADS output
C +   -------------------------------------------------------

       IF (NSTmod.eq.'GRA') THEN

        IF (i.eq.1.and.j.eq.1) THEN
         DO k=1,mz
          NSTgdz(k) = NST__p(i,j,k)
         ENDDO
        ENDIF

        DO k=1,mz
         auxz = NST__p(i,j,k)
         CALL INTlin(INT1Dp,INT1Du,nk+1,auxz,auxu)
         CALL INTlin(INT1Dp,INT1Dv,nk+1,auxz,auxv)
         CALL INTlin(INT1Dp,INT1Dt,nk+1,auxz,auxt)
         CALL INTlin(INT1Dp,INT1Dq,nk+1,auxz,auxq)
         IF (auxz.le.INT1Dp(nk+1)) THEN
          NST__u(i,j,k) = auxu
          NST__v(i,j,k) = auxv
          NST_pt(i,j,k) = auxt
          NST_qv(i,j,k) = auxq
          NST__t(i,j,k) = NST_pt(i,j,k)*(NST__p(i,j,k)/100.)**cap
         ELSE
          NST__u(i,j,k) = 999.999 ! Missing value
          NST__v(i,j,k) = 999.999 ! Avoid extrapolation
          NST_pt(i,j,k) = 999.999
          NST_qv(i,j,k) = 999.999
          NST__t(i,j,k) = 999.999
         ENDIF
        ENDDO

       ENDIF

 
      ENDDO
      ENDDO
      

C +---Impose stability of layers for the equiv. potential temp.
C +   ---------------------------------------------------------


c     DO j=1,my
c     DO i=1,mx
c      DO k=1,mz
c       NST__t(i,j,k)=NST_pt(i,j,k)/(100./NST__p(i,j,k))**cap
c       EQtemp(k)    =NST_pt(i,j,k)
c    .               *EXP(h2olv*NST_qv(i,j,k)/cp/NST__t(i,j,k))
c      ENDDO
c
c      DO k=mz-1,1,-1
c       IF (EQtemp(k).lt.EQtemp(k+1)) THEN
c        EQtemp(k)    =EQtemp(k+1)
c        NST_pt(i,j,k)=MAX(NST_pt(i,j,k),EQtemp(k))
c        NST__t(i,j,k)=NST_pt(i,j,k)/(100./NST__p(i,j,k))**cap
c        NST_qv(i,j,k)=cp*NST__t(i,j,k)/h2olv
c    .                *LOG(EQtemp(k)/NST_pt(i,j,k))
c       ENDIF
c      ENDDO
c     ENDDO
c     ENDDO
      

C +---Compute real temperature
C +   ------------------------
      DO j=1,my
      DO i=1,mx

       DO k=1,mz
        NST__t(i,j,k)=NST_pt(i,j,k)*exp(cap*log(NST__p(i,j,k)/100.))
       ENDDO

      ENDDO
      ENDDO

C +---Filtering of the surface temperature above sea ice and land
C +   -----------------------------------------------------------
       DO j=1,my
       DO i=1,mx

        ! Filtering of STL1 from the ECMWF reanalysis 

        if(LSCmod.eq.'E40'.or.LSCmod.eq.'ECM'.or.LSCmod.eq.'GCM'
     .                                       .or.LSCmod.eq.'NCP') then

         if (NSTsol(i,j).ge.3) then
          NST_st(i,j) = max(NST_st(i,j),NST__t(i,j,mz)-10.)
          NST_st(i,j) = min(NST_st(i,j),NST__t(i,j,mz)+10.)
         endif

         if (NSTsol(i,j).le.2) then
          NST_st(i,j) = max(NST_st(i,j),NST__t(i,j,mz)-15.)
          NST_st(i,j) = min(NST_st(i,j),NST__t(i,j,mz)+15.)
          NST_st(i,j) =     NSTsic(i,j)  * min(273.15,NST_st(i,j))
     .                + (1.-NSTsic(i,j)) * max(270.15,NST_st(i,j))
         endif
        endif

        if (NSTsol(i,j).ge.3) NSTsst(i,j) = NST_st(i,j)

      ENDDO
      ENDDO
      
C!+CA NOR / decreaseSIC : decrease sea-ice extent > WARNING 
      decreaseSIC=.false.
      if (decreaseSIC) then
       Do k=1,10
       
         Do j=1,my
         Do i=1,mx
           INtmp1(i,j)=0
           im1=max(i-1,1)
           ip1=min(i+1,mx)
           jm1=max(j-1,1)
           jp1=min(j+1,my)
           if (NSTsol(i,j).le.2) then ! sea
             INtmp1(i,j)=NSTsic(i,j)
             if (NSTsol(im1,j).le.2) then
               INtmp1(i,j)=min(INtmp1(i,j),NSTsic(im1,j))
             endif
             if (NSTsol(ip1,j).le.2) then
               INtmp1(i,j)=min(INtmp1(i,j),NSTsic(ip1,j))
             endif
             if (NSTsol(i,jm1).le.2) then
               INtmp1(i,j)=min(INtmp1(i,j),NSTsic(i,jm1))
             endif
             if (NSTsol(i,jp1).le.2) then
               INtmp1(i,j)=min(INtmp1(i,j),NSTsic(i,jp1))
             endif
           endif
         EndDo
         EndDo
         
         Do j=1,my
         Do i=1,mx
           NSTsic(i,j) = INtmp1(i,j)
         EndDo
         EndDo
         
       EndDo
      endif


C +---Correct surface pressure <==> Z600 NST = Z600 LSC    
C +   =================================================
      IF (CORzz6.and.NSTmod.ne.'GRA') THEN
 
C +---Geopotential of 600 hPa level in the NST data
C +   ---------------------------------------------
 
       DO k = mz,1,-1 !(begin at surface)
         DO j = 1,my
         DO i = 1,mx
            INtmp1(i,j)=NST_pt(i,j,k)
            INtmp2(i,j)=NST_qv(i,j,k)
            INtmp3(i,j)=NST__p(i,j,k)
         END DO
         END DO
 
         CALL NSTzz6(INtmp1,INtmp2,NST_sh,NST_sp,INtmp3,k,mx,my,mz,
     .               NSTpk6, NSTpx1, NSTlp1, NSTiZp, iZterm, NST_z6)
 
       ENDDO
 
       IF(icheck.ge.1) THEN
         WRITE(*,*) 'NST surf press at chk pt', NST_sp(ipchk,jpchk)
         WRITE(*,*) 'INT Z600       at chk pt', INT_z6(ipchk,jpchk)
         WRITE(*,*) 'NST Z600       at chk pt', NST_z6(ipchk,jpchk)
         WRITE(*,*)
       ENDIF
 
C +---Correct surface pressure
C +   ------------------------
        DO j = 1,my
        DO i = 1,mx
          NST_sp (i,j)= NST_sp(i,j) * (1.0 +    
     .                 (INT_z6(i,j)-NST_z6(i,j)) * grav
     .               / (ra*NSTpk6(i,j)*exp(cap*log(60.)) ))
 
C +..  From Marbaix(2000), Thesis, chapter 3,
C      but with conserverd real temperature as 
C      explained in footnote 4 (not as proposed in
C      the thesis, this one is better !)
 
        END DO
        END DO
       
C +---Update p in atm (3D)
C +   --------------------

      DO j=1,my
      DO i=1,mx

       NST1sp=NST_sp(i,j)
       NST1sh=NST_sh(i,j)

C +         ******
       CALL VERgrd (emptyC,NSTmod,fID,zero,mz,NST1sp,
     .              NST1sh,NST1_t,NST1_p,NST1hp,
     .              NSTgdz,WK1m1D,WK2m1D,WK3m1D)
C +         ******

       DO k=1,mz
        NST__p(i,j,k)=NST1_p(k)
       ENDDO

      ENDDO
      ENDDO

C +---Update potential temperature
C +   ----------------------------
      DO j=1,my
      DO i=1,mx

       DO k=1,mz
         NST_pt(i,j,k)=NST__t(i,j,k)*exp(cap*log(100./NST__p(i,j,k)))
       ENDDO  
       
      ENDDO
      ENDDO

 
C +---Option: Check correction 
C +   ------------------------
 
       IF(icheck.ge.1) THEN
       DO k = mz,1,-1 !(begin at surface)
         DO j = 1,my
         DO i = 1,mx
            INtmp1(i,j)=NST_pt(i,j,k)
            INtmp2(i,j)=NST_qv(i,j,k)
            INtmp3(i,j)=NST__p(i,j,k)
         END DO
         END DO
 
         CALL NSTzz6(INtmp1,INtmp2,NST_sh,NST_sp,INtmp3,k,mx,my,mz,
     .               NSTpk6, NSTpx1, NSTlp1, NSTiZp, iZterm, NST_z6)
 
         ENDDO
         erzmax=0.0
         DO j = 1,my
         DO i = 1,mx
          erzmax = max(erzmax, abs(INT_z6(i,j)-NST_z6(i,j)) )
         END DO
         END DO
         IF (erzmax.GE.1.0) THEN
          write(*,*) 'WARNING (NSTint): '              
          write(*,*) 'Z600 error remains after correction: ',erzmax
          write(*,*) ' (this should not occur)'
         ENDIF
       ENDIF 
       
       IF(icheck.ge.2) THEN
         WRITE(*,*) 'new surf press at chk pt',NST_sp(ipchk,jpchk)
         WRITE(*,*) 'new Z600       at chk pt',NST_z6(ipchk,jpchk)
         WRITE(*,*) 'Z600 error after correction (control): ',erzmax
         WRITE(*,*)
       ENDIF
       
      ENDIF ! END CORzz6 section


      IF (NSTmod.ne.'GRA') THEN

C +---Remove all sursaturations
C +   -------------------------

      CORsat = .false.

      DO j=1,my
      DO i=1,mx
        
        IF (CORsat) THEN
         DO k=1,mz
          qv_max(    k)=0.999*qsat(NST__t(i,j,k),NST__p(i,j,k))
          NST_qv(i,j,k)=MIN(NST_qv(i,j,k),qv_max(k))
         ENDDO
        ENDIF
         
       ENDDO
       ENDDO
       
       IF (CORsat) THEN
        write(*,*) 'WARNING (NSTint): Sursaturation corr.'
        write(*,*)
       ENDIF

C +---Compute levels height
C +   ---------------------

C +         ******
       CALL VERhyd(NST_pt, NST_qv, NST_sh, NST_sp, NST__p,
     .             getpkt, mx, my, mz, NST_zz)
C +         ******


C +---Compute layer depths
C +   --------------------

       DO j=1,my
       DO i=1,mx

        DO k=1,mz
         qv_sat(k)  = qsat(NST__t(i,j,k),NST__p(i,j,k))
         rhoair(k)  = NST__p(i,j,k)*1000./287./NST__t(i,j,k)
        ENDDO
 
        deltaz(mz)  = 0.5*(NST_zz(i,j,k-1)-NST_sh(i,j))
        DO k=2,mz-1
         deltaz(k)  = 0.5*(NST_zz(i,j,k-1)-NST_zz(i,j,k+1))
        ENDDO
        deltaz(1)   = deltaz(2)

      ENDDO
      ENDDO

C +---Increase of specific humidity to take into account cloud cover
C +   --------------------------------------------------------------

      IF (CLDcor .AND. nam_QW.EQ.'-') THEN
       WRITE(6,'(A)') ' Adding parameterized cloud water to Qv'

       DO j=1,my
       DO i=1,mx

         thr_rh = 0.8
         qclmax = 0.0005
         DO k=mz,1,-1
          qv_sat(k) = qsat(NST__t(i,j,k),NST__p(i,j,k)) 
          qv_max(k) = 0.999*qv_sat(k)
          humrel    = NST_qv(i,j,k) / qv_sat(k)
          qcloud(k) = qclmax * exp(-(1.-humrel)/(1.-thr_rh)*3.)
          qwater    = qcloud(k)*rhoair(k)*deltaz(k)
          DO kl=k,mz
           qsuppl   = qwater/rhoair(kl)/deltaz(kl)
           IF ((NST_qv(i,j,kl)+qsuppl).gt.qv_max(kl)) THEN
            qwater  = qwater - (qv_max(kl)-NST_qv(i,j,kl))
     .                         *rhoair(kl)*deltaz(kl)
            qwater  = MAX(qwater,nul)
            NST_qv(i,j,kl) = max(NST_qv(i,j,kl),qv_max(kl))
           ENDIF
          ENDDO
         ENDDO

        ENDDO
        ENDDO

       ENDIF


C +---Increase of specific humidity to take into the interpolated relative humidity
C +   -----------------------------------------------------------------------------

      DO j=1,my
      DO i=1,mx
        
        DO k=1,mz
         qv_max(    k)=qsat(NST__t(i,j,k),NST__p(i,j,k))
         NST_qv(i,j,k)=max(NST_qv(i,j,k),NST_rh(i,j,k)*qv_max(k))
        ENDDO
        
       ENDDO
       ENDDO

C +---Compute equivalent water content
C +   --------------------------------
       DO j=1,my
       DO i=1,mx

        NSTewc(i,j) = 0.
        DO k=1,mz
         NSTewc(i,j)= NSTewc(i,j) + NST_qv(i,j,k)*rhoair(k)*deltaz(k)
        ENDDO

       ENDDO
       ENDDO


      ENDIF  ! {NSTmod.ne.'GRA'}


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Reset variables
C +   ===============

      DO k=1,nk
       WK1_1D(k)=0.
       WK2_1D(k)=0.
       WK3_1D(k)=0.
      ENDDO

      DO k=1,mz
       WK1m1D(k)=0.
       WK2m1D(k)=0.
       WK3m1D(k)=0.
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +   ++++++++++ USEFUL FOR WIND GUST ESTIMATE METHODS +++++++++++++++
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Additional variables in LSC fields provided by MAR simulations
C +   **************************************************************


      IF (SELECT.eq.3) THEN


C +---Horizontal interpolation of surface fields
C +   ==========================================


C +---Surface heat flux
C +   -----------------

      IF (namUTS.ne.'-') THEN

       WRITE(6,'(A,$)') ' 2-D fields : '//namUTS

C +         ******
       CALL UNread (fID,namUTS,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSCuts)
C +         ****** 
       CALL VALchk (namUTS,ni,nj,LSCuts,lwbUTS,upbUTS)
C +         ******
       CALL INThor (HORint,LSC__x,LSC__y,LSCuts,
     .              SPHgrd,NST__x,NST__y,INTuts,
     .              REGgrd,pos_Ox,pos_Oy)
C +         ******

      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Reading and horizontal interpolation
C +   (for each atm. prognostic variable and each level)
C +   ==================================================

      IF (nam_QW.ne.'-'.or.nam_QR.ne.'-'.or.
     .    nam_QI.ne.'-'.or.nam_QS.ne.'-'.or.
     .    namTKE.ne.'-') THEN

       WRITE(6,*)
       WRITE(6,'(A,$)') ' 3-D fields :' 

       DO k = nk,1,-1   !*BEGIN LOOP on vertical levels

        WRITE(6,'(I3,$)') k
        IF (MOD(nk-k+1,20).eq.0) THEN
         WRITE(6,*)
         WRITE(6,'(A,$)') '             '
        ENDIF


C +----Turbulent kinetic energy
C +    ------------------------

C +          ******
        CALL UNread (fID,namTKE,it,k,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSCtke)
C +          ****** 
        CALL VALchk (namTKE,ni,nj,LSCtke,lwbTKE,upbTKE)
C +          ******
        CALL INThor (HORint,LSC__x,LSC__y,LSCtke,
     .               SPHgrd,NST__x,NST__y,INtmp1,
     .               REGgrd,pos_Ox,pos_Oy)
C +          ******
        CALL PUT2D3 (INtmp1,k,INTtke)
C +          ******


C +----Cloud water content
C +    -------------------

C +          ******
        CALL UNread (fID,nam_QW,it,k,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSC_qt)
C +          ****** 
        CALL VALchk (nam_QW,ni,nj,LSC_qt,lwb__Q,upb__Q)
C +          ******

C +          ******
        CALL UNread (fID,nam_QR,it,k,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSCtmp)
C +          ****** 
        CALL VALchk (nam_QR,ni,nj,LSCtmp,lwb__Q,upb__Q)
C +          ******

        DO j=1,nj
        DO i=1,ni
         LSC_qt(i,j)=LSC_qt(i,j)+LSCtmp(i,j)
        ENDDO
        ENDDO

C +          ******
        CALL UNread (fID,nam_QI,it,k,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSCtmp)
C +          ****** 
        CALL VALchk (nam_QI,ni,nj,LSCtmp,lwb__Q,upb__Q)
C +          ******

        DO j=1,nj
        DO i=1,ni
         LSC_qt(i,j)=LSC_qt(i,j)+LSCtmp(i,j)
        ENDDO
        ENDDO

C +          ******
        CALL UNread (fID,nam_QS,it,k,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSCtmp)
C +          ****** 
        CALL VALchk (nam_QS,ni,nj,LSCtmp,lwb__Q,upb__Q)
C +          ******

        DO j=1,nj
        DO i=1,ni
         LSC_qt(i,j)=LSC_qt(i,j)+LSCtmp(i,j)
        ENDDO
        ENDDO

C +          ******
        CALL INThor (HORint,LSC__x,LSC__y,LSC_qt,
     .               SPHgrd,NST__x,NST__y,INtmp1,
     .               REGgrd,pos_Ox,pos_Oy)
C +          ******
        CALL PUT2D3 (INtmp1,k,INT_qt)
C +          ******


       ENDDO
       WRITE(6,*)
       WRITE(6,*)
      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Vertical grid in the NST model (depend on SP)
C +   ==============================

      DO j=1,my
      DO i=1,mx

       NST1sp=NST_sp(i,j)
       NST1sh=NST_sh(i,j)

C +         ******
       CALL VERgrd (emptyC,NSTmod,fID,zero,mz,NST1sp,
     .              NST1sh,NST1_t,NST1_p,NST1hp,
     .              NSTgdz,WK1m1D,WK2m1D,WK3m1D)
C +         ******

       DO k=1,mz
        NST_hp(i,j,k)=NST1hp(k)
        NST__p(i,j,k)=NST1_p(k)
       ENDDO
C +    Note: less memory may be used if _hp is not put in 3D
C            but you must change the position  of this routine

      ENDDO
      ENDDO

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Vertical interpolation
C +   ======================


      DO j=1,my
      DO i=1,mx

       LSC1sp=INT_sp(i,j)
       LSC1sh=INT_sh(i,j)

C +         ******
       CALL VERgrd (LSCmod,emptyC,fID,zero,nk,LSC1sp,
     .              LSC1sh,LSC1_t,LSC1_p,LSC1hp,
     .              LSCgdz,WK1_1D,WK2_1D,WK3_1D)
C +         ******

       DO k=1,nk

        INT1Dz(k)=LSC1hp(    k)

        INT1De(k)=INTtke(i,j,k)
        INT1Dq(k)=INT_qt(i,j,k)

       ENDDO

       INT1Dz(nk+1)=LSC1hp(    nk+1)

       INT1De(nk+1)=INTtke(i,j,nk)
       INT1Dq(nk+1)=INT_qt(i,j,nk)


C +---Linear interpolation (default)
C +   ------------------------------

       IF (VERint.eq.1) THEN

        DO k=1,mz
         auxz = NST_hp(i,j,k)
         CALL INTlin(INT1Dz,INT1De,nk+1,auxz,auxe)
         CALL INTlin(INT1Dz,INT1Dq,nk+1,auxz,auxq)
         NSTtke(i,j,k) = auxe
         NST_qt(i,j,k) = auxq
        ENDDO

       ENDIF

 
C +---Natural cubic spline (optional)
C +   -------------------------------

       IF (VERint.eq.3) THEN

        CALL SPLINE(INT1Dz,INT1De,nk+1,1.E30,1.E30,WK1De)
        CALL SPLINE(INT1Dz,INT1Dq,nk+1,1.E30,1.E30,WK1Dq)

        DO k =1,mz
         auxz = NST_hp(i,j,k)
         CALL SPLINT(INT1Dz,INT1De,WK1Dt,nk+1,auxz,auxe)
         CALL SPLINT(INT1Dz,INT1Dq,WK1Dq,nk+1,auxz,auxq)
         NSTtke(i,j,k) = auxe
         NST_qt(i,j,k) = auxq
        ENDDO

       ENDIF


      ENDDO
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Reset variables
C +   ===============

      DO k=1,nk
       WK1_1D(k)=0.
       WK2_1D(k)=0.
       WK3_1D(k)=0.
      ENDDO

      DO k=1,mz
       WK1m1D(k)=0.
       WK2m1D(k)=0.
       WK3m1D(k)=0.
      ENDDO


      ENDIF  !  {SELECT.eq.3}


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +   +++++++++ END OF READING ADDITIONAL VARIABLES FOR WGE ++++++++++
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Close the NetCDF file
C +   =====================

C +        *******
      CALL UNclose (fID)
C +        *******


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END


C     +--------------------------------------------------------------+
      SUBROUTINE put2D3 (var2D,lev,var3D)
C     +--------------------------------------------------------------+
 
C +   ** variables dimensions
      include 'NSTdim.inc'
 
      REAL var2D (mx,my)
      REAL var3D (mx,my,nk+1)
 
      DO j=1,my
      DO i=1,mx
       IF (ABS(var2D(i,j)).gt.1.E+29) THEN
C ...   Look for strange values : wrong input file ?
C ...   (side effect of this routine)
        WRITE(*,*) 'LSget - put2D3 control :'
        WRITE(*,*) 'Strange value at i,j,lev ='
        WRITE(*,*) ' ',i,j,lev
        WRITE(*,*) 'Value = ',var2D (i,j)
        STOP
       ENDIF
       var3D (i,j,lev) = var2D (i,j)
      ENDDO
      ENDDO
 
      RETURN
      END

C     +--------------------------------------------------------------+
      SUBROUTINE LSuCHG (var2D,unitfact)
C     +--------------------------------------------------------------+
 
      include 'NSTdim.inc'
 
      REAL var2D (ni,nj)
      REAL unitfact
 
      DO j=1,nj
      DO i=1,ni
       var2D (i,j) = var2D (i,j) * unitfact
      ENDDO
      ENDDO
 
      RETURN
      END


C     +--------------------------------------------------------------+
      FUNCTION qsat(tt,pr)
C     +--------------------------------------------------------------+

C     Function qsat computes the Saturation Specific Humidity (kg/kg)

      DATA r273p1/273.16/

      IF (tt.ge.273.16) THEN

       esat =  6.1078 * exp (5.138*log(r273p1/tt))
     .                * exp (6827.*(1.0/r273p1-1.0/tt))

C +... esat : saturated vapor pressure with respect to water
C +... Dudhia (1989) MWR, (B1) and (B2) p.3103
C +... See also Pielke (1984), p.234 and Stull (1988), p.276

      ELSE

       esat =  6.107  * exp (6150.*(1.0/r273p1-1.0/tt))

C +... esat : saturated vapor pressure with respect to ice
C +... Dudhia (1989) MWR, 1989, (B1) and (B2) p.3103

      ENDIF 

      qsat = 0.622*esat/(10.*pr-0.378*esat)
C +...pr : pressure (kPa) multiplied by 10. -> hPa

      RETURN
      END



C     +--------------------------------------------------------------+
      FUNCTION esat(tt,pr)   ! hPa
C     +--------------------------------------------------------------+

C     Function qsat computes the Saturation Specific Humidity (kg/kg)

      DATA r273p1/273.16/

      IF (tt.ge.273.16) THEN

       esat =  6.1078 * exp (5.138*log(r273p1/tt))
     .                * exp (6827.*(1.0/r273p1-1.0/tt))

C +... esat : saturated vapor pressure with respect to water
C +... Dudhia (1989) MWR, (B1) and (B2) p.3103
C +... See also Pielke (1984), p.234 and Stull (1988), p.276

      ELSE

       esat =  6.107  * exp (6150.*(1.0/r273p1-1.0/tt))

C +... esat : saturated vapor pressure with respect to ice
C +... Dudhia (1989) MWR, 1989, (B1) and (B2) p.3103

      ENDIF 

      RETURN
      END


C     +--------------------------------------------------------------+
      SUBROUTINE VALchk (varname,nx,ny,var,lwb,upb)
C     +--------------------------------------------------------------+

      IMPLICIT NONE
 
      INTEGER     nx,ny,i,j,ierror,ipe,jpe
      REAL        var (nx,ny)
      REAL        lwb,upb
      CHARACTER*7 varname
 
      ierror = 0

      DO j=1,ny
      DO i=1,nx
       IF (var(i,j).lt.lwb .or.
     .     var(i,j).gt.upb) THEN
         ierror = ierror+1 
         ipe    = i
         jpe    = j
       ENDIF
      ENDDO
      ENDDO

      IF (ierror.ge.1) THEN
       write(6,*)
       write(6,*)
       write(6,*) 'The range of values for the variable'
       write(6,*) '  ',varname,' in large-scale fields'
       write(6,*) 'is probably incorrect (out of specified'
       write(6,*) 'bounds). Please check it before running'
       write(6,*) 'Error is found on',ierror,'points'
       write(6,*) 'such as (i,j)=', ipe,jpe
       write(6,*) 'value        =', var(ipe,jpe)
       write(6,*)
       write(6,*) 'NESTOR.                        --- STOP'
       write(6,*)
       STOP
      ENDIF

      RETURN
      END
      
C     +--------------------------------------------------------------+
      SUBROUTINE STGgrd (LSC1Dx,LSC1Dy,LSCtmp,LSCtm2,ni,nj)
C     +
C     + This is merily a patch to allow staggered horizontal grids
C     +  when the grid is "regular", i.e. rectangular in lat/lon
C     +  but possibly "zoomed" the LMD-way.
C     +  Coordinates are then taken from the NetCDF "dimensions",
C     +  so that it will not work for irregular grids.
C     +--------------------------------------------------------------+

      REAL LSC1Dx(ni)   ,LSC1Dy(nj)
      REAL LSCtmp(ni,nj),LSCtm2(ni,nj)

       DO i=1,ni
        LSCtmp(i,1)=LSC1Dx(i)
       ENDDO

       DO j=1,nj
        LSCtm2(1,j)=LSC1Dy(j)
       ENDDO

      RETURN
      END

