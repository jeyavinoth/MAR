C   +-------------------------------------------------------------------+
C   |  Subroutine PRCdes                            21/09/2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Disaggregation of precipitation fields.                           |
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
C   | OUTPUT : Precipitation from desagregation of LSC fields:          |
C   | ^^^^^^^^ - NSTpr1 : no water conservation                         |
C   |          - NSTpr2 : water conservation over global domain         | 
C   |          - NSTpr3 : water conservation at local scale + global    | 
C   |                                                                   | 
C   +-------------------------------------------------------------------+

      SUBROUTINE PRCdes

 
      IMPLICIT NONE


C +---Include files
C +   -------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'CTRvar.inc'
      INCLUDE 'LSCvar.inc'
      INCLUDE 'INTvar.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'DESvar.inc'

C +---Local variables
C +   ---------------

      INTEGER i,j,k,klev,it,fID,ierror,nstep,jstag,
     .        l,izero,kINI,kLCL,MXlist,ind,itloc,
     .        KLO,KHI,KPL,ntSlow,i1,i2,j1,j2,k1,ii,jj,kk,ngrid,
     .        iloc,jloc,nt,ii1,ii2,jj1,jj2,kk1,kk2,ll1,ll2,idiv,
     .        jdiv,nx,TkLCL(mx,my),ic1,ic2,jc1,jc2,mmx,mmy,istag

      PARAMETER (MXlist=10000)
      PARAMETER (nstep =   10)

      INTEGER pos_Ox(mx,my),pos_Oy(mx,my)

      REAL    exxpo,gamTz,empty1(1),cap,ra,rv,grav,qv_max,
     .        qsat,cp,h2olv,zLCL,zETL,zLFC,aBuoyE,pBuoyE,hPBL,
     .        rzero,lambds,valhum,rhmin,NSTdhx,NSTdhy,dz,aux,
     .        dztot,Ubl,Vbl,Wbl,VgradS,omegaS,omega1,height,
     .        TKEsrf,gamma,dst_dx,dst_dy,pr_top,pr_tot,Uint,
     .        dt_max,alpha,factpr,fact_w,rhcrHY,pfreeze,Vint,
     .        prItot,pr_avg,prIavg,ratio,ezs,tzs,xilon,yilat,
     .        valint,possOx,possOy,lifetime,x0lon,y0lat,dist,
     .        sumtmp,tot_wf,Wtot,sigma,degrad,esat,efficiency,
     .        auxz,auxu,auxv,auxw,auxt,auxq,auxe,nbr_pr,RHbl,
     .        lon1,lon2,dlon,lat1,lat2,dlat,fdis,aux1,aux2,
     .        coslat,LSCdhx,LSCdhy,VgradL,omegaL,dt,max_V,eps,
     .        tfreeze,vrain,vsnow,uloc,vloc,xloc,yloc,area,
     .        dxic,dyic,tfall,hgtmax,hgt1,hgt2,dtSlow,deltap,
     .        tloc,val_pr,gravit,xcent,ycent,wmLCL,wLFC,BEneg2,
     .        stbLCL,wlow,w_conv,third,wthres,wLCL,timecv,dur_dt,
     .        deltaz,sflux,wwprim,omegaT,switch_srf,int_buo,w__CNV,
     .        int_rib,int_u,eps3,zbuo,ribc,zLCLkf,auxEff,Effic2,
     .        TIMcor_orog,avwind,compt,buotot,avgU,avgV,TIMcor_deep,
     .        Qbl,rhoCNV,qvsCNV,diffqv,rate_CNV,wturb_max,LSC_dx,LSC_dy

      REAL    INT1Dk(nk+1),INT__w(mx,my,nk+1),INT1Dw(nk+1),
     .        LSCdpr(ni,nj),NSTlpr(mx,my),ratepr(mx,my),
     .        Us_avg(mx,my),Vs_avg(mx,my),Ws_avg(mx,my),
     .        wfct  (nstep),precip(nstep),TMP_pr(mx,my,2),
     .        rateDC(mz),TT_Envr(mx,my,mz),dstLB (mx,my,2),
     .        TaBuoyE(mx,my),TpBuoyE(mx,my),TzLCL(mx,my),
     .        TzLFC(mx,my),TzETL(mx,my),TAdiabT(mx,my,mz),
     .        TMP2Du(mx,my),TMP2Dv(mx,my),tmp_2D(mx,my),
     .        LSC1sp,LSC1_t(nk),LSC1_p(nk+1),LSC1hp(nk+1),LSC1sh,
     .        NST1sp,NST1_t(mz),NST1_p(mz+1),NST1hp(mz+1),NST1sh,
     .        WK1_1D(nk),WK2_1D(nk),WK3_1D(nk),relhum(mz),
     .        WK1m1D(mz),WK2m1D(mz),WK3m1D(mz),vfall(mz),
     .        T_Envr(mz),QvEnvr(mz),H_Envr(mz),AdiabT(mz),
     .        P_Envr(mz),lambda(mz),omegaP(mz),qvsat(mz),
     .        omegaQ(mz),Ffact(mz),timef(mz),rateqv(mz),
     .        dhx1(mx,my),dhx2(mx,my),dhy1(mx,my),dhy2(mx,my),
     .        w_Updr(mz),wturb(mz),zinv(mx,my),tetav(mx,my,mz),
     .        buoy(mx,my,mz),rib(mx,my,mz),zrib(mx,my)

      REAL    LSC_sf(ni,nj),r(mz),rst(mz),rstfreeze,rstsrf,Nfre
      INTEGER kfreeze
   
      LOGICAL compPR,REGnst,STEnst,BOLcel(mx,my),lev_OK(mz),
     .        global_conserv,input_filt,output_filt,conv_adjust,
     .        deep_conv,surf_flux,vdelb_model,local_conserv

      CHARACTER*3   empty
      CHARACTER*6   nam_SH,nam_SP,namST1,namST2,namSW1,nam10U,nam10V,
     .              nam__U,nam__V,nam__T,nam__Q,namTKE,nam__W,nam_PR,
     .              namUTS,nam_SF
      CHARACTER*10  var_units
      CHARACTER*100 LSCtit
 

C +---Physical constants
C +   ------------------

      DATA ra       /  287.     e0   /
      DATA rv       /  461.     e0   /
      DATA cp       / 1004.     e0   /
      DATA gravit   /    9.81   e0   /
      DATA h2olv    /    2.5000e+6   /
      DATA cap      /    0.28586e0   /
      DATA grav     /    9.81   e0   /
      DATA third    /    0.333333333 /
      DATA wthres   /    0.0         /
      DATA empty    /   '   '        /
      DATA izero    /    0           /
      DATA rzero    /    0.          /
      DATA rhcrHY   /    1.00        /
      DATA pr_top   /20000.0         /
      DATA gamma    /    0.7         /         
      DATA degrad   / 1.745329252e-2 /
      DATA lifetime / 2700.          /           
      DATA tfreeze  /  273.15        /
      DATA vrain    /    5.0         /          
      DATA vsnow    /    1.0         /            
      DATA eps3     /    0.001       /
      DATA ribc     /    0.6         /


C +---Initialisation
C +   --------------

      lfirst_LSC = .true.
      lfirst_NST = .true.

      local_conserv  = .true. 
      global_conserv = .true. 
      deep_conv      = .false.
      surf_flux      = .false.
      input_filt     = .false.
      output_filt    = .false.
      vdelb_model    = .true.

      mmx = mx
      mmy = my
      ic1 = MIN(2,mmx)
      ic2 = MAX(1,mmx-1)
      jc1 = MIN(2,mmy)
      jc2 = MAX(1,mmy-1)


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

      write(6,*) 'Horizontal and vertical interpolations'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*) 'Open file  : ',LSCfil
      write(6,*) 'Time step  : ',I_time


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Define variable names
C +   =====================

      IF (LSCmod.eq.'MAR') THEN
       nam_SH='sh'
       nam_SP='pstar'
       namST1='tairSL'
       namST2='tairSL'
       namSW1='-'   
       nam10U='-'  
       nam10V='-'  
       nam__U='uairDY'
       nam__V='vairDY'
       nam__W='wairDY'
       nam__T='tairDY'
       nam__Q='qvDY'
       namTKE='-'
       namUTS='hsenSL'
c      namUTS='hlatSL'
       nam_PR='rainHY'
       nam_SF='snowHY'
       fact_w=0.01  ! factor for NST__w : cm/s -> m/s
      ELSE
       nam_SH='SH'
       nam_SP='SP'
       namST1='STL1'
       namST2='STL2'
       namSW1='SWL1'
       nam10U='10U'
       nam10V='10V'
       nam__U='U'
       nam__V='V'
       nam__W='-'
       nam__T='T'
       nam__Q='Q'
       namTKE='-'
       namUTS='-'
       nam_PR='-'
       nam_SF='-'
       fact_w=1.
      ENDIF
      factpr=1000. ! factor for NST__p and NST_sp : kPa -> Pa


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

       DO j=1,nj
       DO i=1,ni
        LSC__x(i,j)=LSC1Dx(i)
        LSC__y(i,j)=LSC1Dy(j)
       ENDDO
       ENDDO

      ELSE

C +         ******
       CALL UNread (fID,'lon' ,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__x)
C +         ******
       CALL UNread (fID,'lat' ,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__y)
C +         ****** 

      ENDIF


      IF (LSCmod.eq.'MAR') THEN
       LSC_dx = (LSC1Dx(2)-LSC1Dx(1))*1000.
       LSC_dy = (LSC1Dy(2)-LSC1Dy(1))*1000.
      ELSE
       LSC_dx = ABS(LSC__x(ni/2+1,nj/2)-LSC__x(ni/2,nj/2))
     .        * 111111.11111*ABS(COS(LSC__y(ni/2,nj/2)/degrad))
       LSC_dy = ABS(LSC__y(ni/2,nj/2+1)-LSC__y(ni/2,nj/2))
     .        * 111111.11111
      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Initialize vertical grid
C +   ========================

      DO k=1,nk
       LSCgdz(k)=0.
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Initialization of precipitation variable
C +   ========================================


      compPR=.true.

      IF (iter.eq.1) THEN

       IF ((it-1).gt.0.and.nam_PR.ne.'-') THEN

C +          ******
        CALL UNread (fID,nam_PR,it-1,1,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSC_pr)
C +          ****** 

        compPR=.true.

       ELSE

        DO j=1,nj
        DO i=1,ni
         LSC_pr(i,j)=0.
        ENDDO
        ENDDO
        compPR=.false.

       ENDIF

      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Topography
C +   ==========


C +        ******
      CALL UNread (fID,nam_SH,it,1,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSC_sh)
C +        ****** 


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Horizontal interpolation of surface fields
C +   ==========================================


C +---Surface Pressure
C +   ----------------

      WRITE(6,'(A,$)') ' 2-D fields : '//nam_SH//' - '//nam_SP

C +        ******
      CALL UNread (fID,nam_SP,it,1,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSC_sp)
C +        ****** 

      IF (LSCmod.ne.'MAR') THEN
C +         ******
       CALL LSuCHG (LSC_sp,1.E-3) !(Change units: Pa-->KPa)
C +         ******
      ENDIF


C +---Soil or Sea surface temperature
C +   -------------------------------

      WRITE(6,'(A,$)') ' - '//namST1

C +        ******
      CALL UNread (fID,namST1,it,1,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSC_st)
C +        ****** 


C +---Soil or Sea temperature
C +   -----------------------

      WRITE(6,'(A,$)') ' - '//namST2

C +        ******
      CALL UNread (fID,namST2,it,1,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSCdst)
C +        ****** 


C +---Surface heat flux
C +   -----------------

      IF (namUTS.ne.'-') THEN

       WRITE(6,'(A,$)') ' - '//namUTS

C +         ******
       CALL UNread (fID,namUTS,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSCuts)
C +         ****** 

      ELSE

       DO j=1,nj
       DO i=1,ni
        LSCuts(i,j) = 0.
       ENDDO
       ENDDO

      ENDIF


C +---Include ice cover over ocean
C +   ----------------------------

c     DO j=1,my
c     DO i=1,mx
c      IF (NSTsol(i,j).le.2.and.INT_st(i,j).le.271.2)
c    .  NSTsol(i,j)=2
c     ENDDO
c     ENDDO


C +---Include snow cover over land
C +   ----------------------------

c     DO j=1,my
c     DO i=1,mx
c      IF (NSTsol(i,j).ge.3.and.INT_st(i,j).lt.273.2)
c    .  NSTsol(i,j)=3
c     ENDDO
c     ENDDO


C +---Temperature difference between 1st atm. level and soil/sea
C +   ----------------------------------------------------------

      WRITE(6,'(A,$)') ' - '//nam__T

C +        ******
      CALL UNread (fID,nam__T,it,nk,bi,bj,ni,nj,1,
     .             LSC1Dx,LSC1Dy,empty1,var_units,LSC__t)
C +        ****** 

      DO j=1,nj
      DO i=1,ni
       LSC_dt(i,j)=LSC__t(i,j)-LSC_st(i,j)
      ENDDO
      ENDDO

      DO j=1,nj
      DO i=1,ni
       LSC_pt(i,j)=LSC__t(i,j)*(100./LSC_sp(i,j))**cap 
      ENDDO
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Horizontal interpolation of surface variables
C +   =============================================


      DO j=1,my
      DO i=1,mx

       INT_sh(i,j) = auxL2N(i,j)       
     .             * (   auyL2N(i,j) *LSC_sh(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_sh(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j) *LSC_sh(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_sh(iiL2N(i,j)  ,jjL2N(i,j)  ))

       INT_sp(i,j) = auxL2N(i,j)       
     .             * (   auyL2N(i,j) *LSC_sp(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_sp(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j) *LSC_sp(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_sp(iiL2N(i,j)  ,jjL2N(i,j)  ))

       INT_st(i,j) = auxL2N(i,j)       
     .             * (   auyL2N(i,j) *LSC_st(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_st(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j) *LSC_st(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_st(iiL2N(i,j)  ,jjL2N(i,j)  ))

       INTuts(i,j) = auxL2N(i,j)       
     .             * (   auyL2N(i,j) *LSCuts(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSCuts(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j) *LSCuts(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSCuts(iiL2N(i,j)  ,jjL2N(i,j)  ))

       INTdst(i,j) = auxL2N(i,j)       
     .             * (   auyL2N(i,j) *LSCdst(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSCdst(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j) *LSCdst(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSCdst(iiL2N(i,j)  ,jjL2N(i,j)  ))

       INT_dt(i,j) = auxL2N(i,j)       
     .             * (   auyL2N(i,j) *LSC_dt(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_dt(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j) *LSC_dt(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_dt(iiL2N(i,j)  ,jjL2N(i,j)  ))

       INT_pt(i,j,nk) = auxL2N(i,j)       
     .             * (   auyL2N(i,j) *LSC_pt(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_pt(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j) *LSC_pt(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSC_pt(iiL2N(i,j)  ,jjL2N(i,j)  ))

      ENDDO
      ENDDO


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


C +---Correction of soil temperature according to topography change
C +   =============================================================

      DO j = 1,my
      DO i = 1,mx
       IF (NST_sh(i,j).gt.1.) THEN
        NST_st(i,j)=INT_pt(i,j,nk)/(100./INT_sp(i,j))**cap
     .             -INT_dt(i,j)
C +...  Temperature diff. between 1st level and surface is conserved
       ELSE
        NST_st(i,j)=INT_st(i,j)
C +...  No correction for the sea surface temperature
       ENDIF
        NSTdst(i,j)=INTdst(i,j)
        NSTuts(i,j)=INTuts(i,j)
      ENDDO
      ENDDO

 
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Atmospheric variables at surface:
C +   (bottom boundary for vertic. interpolation)
C +   ===========================================       


C +---10-m U-wind
C +   -----------

      IF (nam10U.eq.'-') THEN

       WRITE(6,'(A,$)') ' - '//nam__U
C +         ******
       CALL UNread (fID,nam__U,it,nk,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__u)
C +         ****** 

      ELSE

       WRITE(6,'(A,$)') ' - '//nam10U
C +         ******
       CALL UNread (fID,nam10U,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__u)
C +         ****** 

      ENDIF

      DO j=1,my
      DO i=1,mx
       TMP2Du(i,j)= auxL2N(i,j)       
     .            * (   auyL2N(i,j) *LSC__u(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__u(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .            + (1.0-auxL2N(i,j))
     .            * (   auyL2N(i,j) *LSC__u(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__u(iiL2N(i,j)  ,jjL2N(i,j)  ))
      ENDDO
      ENDDO


C +---10-m V-wind
C +   -----------

      IF (nam10V.eq.'-') THEN

       WRITE(6,'(A,$)') ' - '//nam__V
C +         ******
       CALL UNread (fID,nam__V,it,nk,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__v)
C +         ****** 

      ELSE

       WRITE(6,'(A,$)') ' - '//nam10V
C +         ******
       CALL UNread (fID,nam10V,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__v)
C +         ****** 

      ENDIF

      DO j=1,my
      DO i=1,mx
       TMP2Dv(i,j)= auxL2N(i,j)       
     .            * (   auyL2N(i,j) *LSC__v(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__v(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .            + (1.0-auxL2N(i,j))
     .            * (   auyL2N(i,j) *LSC__v(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__v(iiL2N(i,j)  ,jjL2N(i,j)  ))
      ENDDO
      ENDDO


C +---Wind vector rotation (according to projection)
C +   --------------------

C   IF (NST_dx.gt.0.01)
C +        ******
c   .CALL VecRot (NST__x,NST__y,NST_dx,TMP2Du,TMP2Dv)
C +        ******

C +        ******
      CALL PUT2D3 (TMP2Du,nk+1,INT__u)
C +        ******

C +        ******
      CALL PUT2D3 (TMP2Dv,nk+1,INT__v)
C +        ******


C +---Water vapour
C +   ------------

      IF (namSW1.eq.'-') THEN

       WRITE(6,'(A,$)') ' - '//nam__Q
C +         ******
       CALL UNread (fID,nam__Q,it,nk,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_qv)
C +         ****** 

      ELSE

       WRITE(6,'(A,$)') ' - '//namSW1
C +         ******
       CALL UNread (fID,namSW1,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_qv)
C +         ****** 

      ENDIF

C +        ******
      CALL LSuCHG (LSC_qv,1.)
C +        ******

      DO j=1,my
      DO i=1,mx
       INT_qv(i,j,nk+1)= auxL2N(i,j)       
     .           * (   auyL2N(i,j) *LSC_qv(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .             +(1-auyL2N(i,j))*LSC_qv(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .           + (1.0-auxL2N(i,j))
     .           * (   auyL2N(i,j) *LSC_qv(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .             +(1-auyL2N(i,j))*LSC_qv(iiL2N(i,j)  ,jjL2N(i,j)  ))
      ENDDO
      ENDDO


C +---Potential temperature
C +   ---------------------

      IF (namST1.eq.'-') THEN

       WRITE(6,'(A,$)') ' - '//nam__T
C +         ******
       CALL UNread (fID,nam__T,it,nk,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__t)
C +         ****** 

      ELSE

       WRITE(6,'(A,$)') ' - '//namST1
C +         ******
       CALL UNread (fID,namST1,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__t)
C +         ****** 

      ENDIF

      DO j=1,nj
      DO i=1,ni
       LSC_pt(i,j)=LSC__t(i,j)*(100./LSC_sp(i,j))**cap 
      ENDDO
      ENDDO

      DO j=1,my
      DO i=1,mx
       INT_pt(i,j,nk+1)= auxL2N(i,j)       
     .           * (   auyL2N(i,j) *LSC_pt(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .             +(1-auyL2N(i,j))*LSC_pt(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .           + (1.0-auxL2N(i,j))
     .           * (   auyL2N(i,j) *LSC_pt(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .             +(1-auyL2N(i,j))*LSC_pt(iiL2N(i,j)  ,jjL2N(i,j)  ))

      ENDDO
      ENDDO

C +---Precipitation
C +   -------------

      IF (nam_PR.ne.'-') THEN

       DO j=1,nj
       DO i=1,ni
        LSCppr(i,j)=LSC_pr(i,j)
       ENDDO
       ENDDO

       WRITE(6,'(A,$)') ' - '//nam_PR
C +         ******
       CALL UNread (fID,nam_PR,it,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_pr)
C +         ****** 

       IF (nam_SF.ne.'-') THEN
        
        WRITE(6,'(A,$)') ' - '//nam_SF
C +          ******
        CALL UNread (fID,nam_SF,it,1,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSC_sf)
C +          ******       
        DO j=1,nj
        DO i=1,ni
         LSC_pr(i,j)=LSC_pr(i,j)+LSC_sf(i,j)
        ENDDO
        ENDDO        
       END IF

       IF (compPR) THEN

        DO j=1,nj
        DO i=1,ni
         LSCdpr(i,j)=LSC_pr(i,j)-LSCppr(i,j)
        ENDDO
        ENDDO

       ELSE

        DO j=1,nj
        DO i=1,ni
         LSCdpr(i,j)=0.
        ENDDO
        ENDDO

       ENDIF

       DO j=1,my
       DO i=1,mx
        NSTIpr(i,j)=auxL2N(i,j)       
     .            *(   auyL2N(i,j) *LSCdpr(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .             +(1-auyL2N(i,j))*LSCdpr(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .            +(1.0-auxL2N(i,j))
     .            *(   auyL2N(i,j) *LSCdpr(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .             +(1-auyL2N(i,j))*LSCdpr(iiL2N(i,j)  ,jjL2N(i,j)  ))
       ENDDO
       ENDDO

      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Reading and horizontal interpolation (+rotation) 
C +   (for each atm. prognostic variable and each level)
C +   ==================================================

      WRITE(6,*)
      WRITE(6,'(A,$)') ' 3-D fields :' 

      DO k = nk,1,-1   !*BEGIN LOOP on vertical levels

       WRITE(6,'(I3,$)') k
       IF (k.eq.(nk/2)) THEN
        WRITE(6,*)
        WRITE(6,'(A,$)') '             '
       ENDIF


C +----U-Wind
C +    ------

C +         ******
       CALL UNread (fID,nam__U,it,k,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__u)
C +         ****** 

       DO j=1,my
       DO i=1,mx
        TMP2Du(i,j)=auxL2N(i,j)       
     .             *(   auyL2N(i,j) *LSC__u(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__u(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             +(1.0-auxL2N(i,j))
     .             *(   auyL2N(i,j) *LSC__u(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__u(iiL2N(i,j)  ,jjL2N(i,j)  ))
       ENDDO
       ENDDO


C +----V-Wind
C +    ------
 
C +         ******
       CALL UNread (fID,nam__V,it,k,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__v)
C +         ****** 

       DO j=1,my
       DO i=1,mx
        TMP2Dv(i,j)=auxL2N(i,j)       
     .             *(   auyL2N(i,j) *LSC__v(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__v(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             +(1.0-auxL2N(i,j))
     .             *(   auyL2N(i,j) *LSC__v(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__v(iiL2N(i,j)  ,jjL2N(i,j)  ))
       ENDDO
       ENDDO


C +----Wind vector rotation (according to projection)
C +    --------------------

c      IF (NST_dx.gt.0.01) 
C +         ******
c    . CALL VecRot (NST__x,NST__y,NST_dx,TMP2Du,TMP2Dv)
C +         ******

C +         ******
       CALL PUT2D3 (TMP2Du,k,INT__u)
C +         ******

C +         ******
       CALL PUT2D3 (TMP2Dv,k,INT__v)
C +         ******

 
C +----W-Wind
C +    ------

       IF (nam__W.ne.'-') THEN

C +          ******
        CALL UNread (fID,nam__W,it,k,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSC__w)
C +          ****** 

        DO j=1,my
        DO i=1,mx
         INT__w(i,j,k)=auxL2N(i,j)       
     .             *(   auyL2N(i,j) *LSC__w(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__w(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             +(1.0-auxL2N(i,j))
     .             *(   auyL2N(i,j) *LSC__w(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__w(iiL2N(i,j)  ,jjL2N(i,j)  ))
        ENDDO
        ENDDO

       ELSE

        DO j=1,my
        DO i=1,mx
         INT__w(i,j,k)=0.0
        ENDDO
        ENDDO

       ENDIF


C +----Water vapour
C +    ------------

C +         ******
       CALL UNread (fID,nam__Q,it,k,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_qv)
C +         ****** 

C +         ******
       CALL LSuCHG (LSC_qv,1.)
C +         ******

       DO j=1,my
       DO i=1,mx
        INT_qv(i,j,k)=auxL2N(i,j)       
     .             *(   auyL2N(i,j) *LSC_qv(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC_qv(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             +(1.0-auxL2N(i,j))
     .             *(   auyL2N(i,j) *LSC_qv(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC_qv(iiL2N(i,j)  ,jjL2N(i,j)  ))
       ENDDO
       ENDDO


C +----Turbulent kinetic energy
C +    ------------------------

       IF (namTKE.ne.'-') THEN

C +          ******
        CALL UNread (fID,namTKE,it,k,bi,bj,ni,nj,1,
     .               LSC1Dx,LSC1Dy,empty1,var_units,LSCtke)
C +          ****** 

        DO j=1,my
        DO i=1,mx
         INTtke(i,j,k)=auxL2N(i,j)       
     .              *(   auyL2N(i,j) *LSCtke(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSCtke(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .              +(1.0-auxL2N(i,j))
     .              *(   auyL2N(i,j) *LSCtke(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .               +(1-auyL2N(i,j))*LSCtke(iiL2N(i,j)  ,jjL2N(i,j)  ))
        ENDDO
        ENDDO

       ELSE

        DO j=1,my
        DO i=1,mx
         INTtke(i,j,k)=0.0
        ENDDO
        ENDDO

       ENDIF


C +--- Potential temperature
C +    ---------------------

C +         ******
       CALL UNread (fID,nam__T,it,k,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__t)
C +         ****** 

       DO j=1,nj
       DO i=1,ni

        LSC1sp=LSC_sp(i,j)
        LSC1sh=LSC_sh(i,j)

C +          ******
        CALL VERgrd (LSCmod,empty ,fID,k,nk,LSC1sp,
     .               LSC1sh,LSC1_t,LSC1_p,LSC1hp,
     .               LSCgdz,WK1_1D,WK2_1D,WK3_1D)
C +          ******

        LSC_pt(i,j)=LSC__t(i,j)*(100./LSC1_p(k))**cap 

       ENDDO
       ENDDO

       DO j=1,my
       DO i=1,mx
        INT_pt(i,j,k)=auxL2N(i,j)       
     .             *(   auyL2N(i,j) *LSC_pt(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC_pt(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             +(1.0-auxL2N(i,j))
     .             *(   auyL2N(i,j) *LSC_pt(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC_pt(iiL2N(i,j)  ,jjL2N(i,j)  ))
       ENDDO
       ENDDO


      ENDDO
      WRITE(6,*)
      WRITE(6,*)


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Vertical grid in the NST model (depend on SP)
C +   ==============================

      DO j=1,my
      DO i=1,mx

       NST1sp=NST_sp(i,j)
       NST1sh=NST_sh(i,j)

C +         ******
       CALL VERgrd (empty ,NSTmod,fID,izero,mz,NST1sp,
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

      IF (VERint.ne.1.and.VERint.ne.3) THEN
       write(6,*) 'CAUTION :'
       write(6,*) 'Horizontal interpolation order incorrectly'
       write(6,*) 'specified. Default is set to linear'
       VERint=1
      ENDIF

      DO j=1,my
      DO i=1,mx

       LSC1sp=INT_sp(i,j)
       LSC1sh=INT_sh(i,j)

C +         ******
       CALL VERgrd (LSCmod,empty,fID,izero,nk,LSC1sp,
     .              LSC1sh,LSC1_t,LSC1_p,LSC1hp,
     .              LSCgdz,WK1_1D,WK2_1D,WK3_1D)
C +         ******

       DO k=1,nk

        INT1Dz(k)=LSC1hp(    k)

        INT1Du(k)=INT__u(i,j,k)
        INT1Dv(k)=INT__v(i,j,k)
        INT1Dw(k)=INT__w(i,j,k)
        INT1Dt(k)=INT_pt(i,j,k)
        INT1Dq(k)=INT_qv(i,j,k)
        INT1Dk(k)=INTtke(i,j,k)

       ENDDO

       INT1Dz(nk+1)=LSC1hp(    nk+1)

       INT1Du(nk+1)=0.
       INT1Dv(nk+1)=0.
       INT1Dw(nk+1)=0.
       INT1Dt(nk+1)=INT_pt(i,j,nk)
       INT1Dq(nk+1)=INT_qv(i,j,nk)
       INT1Dk(nk+1)=INTtke(i,j,nk)


C +---Linear interpolation
C +   --------------------

       IF (VERint.eq.1) THEN

        DO k=1,mz

         KLO=1
         KHI=nk+1
 1       IF (KHI-KLO.GT.1) THEN
           KPL=(KHI+KLO)/2
           IF(INT1Dz(KPL).GT.NST_hp(i,j,k))THEN
             KHI=KPL
           ELSE
             KLO=KPL
           ENDIF
         GOTO 1
         ENDIF
         ind=KLO

         fdis = INT1Dz(ind+1)-INT1Dz(ind)

         aux1 = ((INT1Dz(ind+1)-NST_hp(i,j,k))/fdis)
         aux2 = ((NST_hp(i,j,k)-INT1Dz(ind  ))/fdis)

         NST__u(i,j,k) = aux1*INT1Du(ind) + aux2*INT1Du(ind+1)
         NST__v(i,j,k) = aux1*INT1Dv(ind) + aux2*INT1Dv(ind+1)
         NST__w(i,j,k) = aux1*INT1Dw(ind) + aux2*INT1Dw(ind+1)
         NST_pt(i,j,k) = aux1*INT1Dt(ind) + aux2*INT1Dt(ind+1)
         NST_qv(i,j,k) = aux1*INT1Dq(ind) + aux2*INT1Dq(ind+1)
         NSTtke(i,j,k) = aux1*INT1Dk(ind) + aux2*INT1Dk(ind+1)

         IF (NST_hp(i,j,k).LT.INT1Dz(ind  )) THEN
          NST__u(i,j,k) =INT1Du(ind  )
          NST__v(i,j,k) =INT1Dv(ind  )
          NST__w(i,j,k) =INT1Dw(ind  )
          NST_pt(i,j,k) =INT1Dt(ind  )
          NST_qv(i,j,k) =INT1Dq(ind  )
          NSTtke(i,j,k) =INT1Dk(ind  )
         ENDIF
         IF (NST_hp(i,j,k).GT.INT1Dz(ind+1)) THEN
          NST__u(i,j,k) =INT1Du(ind+1)
          NST__v(i,j,k) =INT1Dv(ind+1)
          NST__w(i,j,k) =INT1Dw(ind+1)
          NST_pt(i,j,k) =INT1Dt(ind+1)
          NST_qv(i,j,k) =INT1Dq(ind+1)
          NSTtke(i,j,k) =INT1Dk(ind+1)
         ENDIF

        ENDDO

       ENDIF

 
C +---Compute real temperature and remove all sursaturations
C +   ------------------------------------------------------

       DO k=1,mz
        NST__t(i,j,k)=NST_pt(i,j,k)/(100./NST__p(i,j,k))**cap
        qv_max       =qsat(NST__t(i,j,k),NST__p(i,j,k))
        NST_qv(i,j,k)=MIN(NST_qv(i,j,k),qv_max)
       ENDDO

      ENDDO
      ENDDO


C +---Computation of levels height
C +   ----------------------------

C +        ******
      CALL GEOpot (NST__t,NST_sh,NST_sp,NST__p,NST_zz)
C +        ******


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


C +---Close the NetCDF file
C +   =====================

C +        ******
      CALL NCCLOS (fID,ierror)
C +        ******


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Filtering of thermodynamical fields
C +   ===================================

      IF (input_filt) THEN
      
       eps = 0.10  ! Selectivity parameter

       DO k=1,mz

C +...  1. Temperature
C +     - - - - - - - -

        DO j=1,my
        DO i=1,mx
         tmp_2D(i,j) = NST__t(i,j,k)
        ENDDO
        ENDDO

C +          *********
        CALL DYNfil_2H (tmp_2D,eps)
C +          *********

        DO j=1,my
        DO i=1,mx
         NST__t(i,j,k) = tmp_2D(i,j)
        ENDDO
        ENDDO

C +...  2. Specific humidity
C +     - - - - - - - - - - -

        DO j=1,my
        DO i=1,mx
         tmp_2D(i,j) = NST_qv(i,j,k)
        ENDDO
        ENDDO

C +          *********
        CALL DYNfil_2H (tmp_2D,eps)
C +          *********

        DO j=1,my
        DO i=1,mx
         NST_qv(i,j,k) = tmp_2D(i,j)
        ENDDO
        ENDDO

       ENDDO

      ENDIF
      

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Disaggregation model of precipitation (Sinclair, 1994)
C +   ======================================================


C +---Initialisation
C +   --------------

      DO j=1,my
      DO i=1,mx
       Us_avg(i,j) = 0.d0
       Vs_avg(i,j) = 0.d0
       Ws_avg(i,j) = 0.d0
      ENDDO
      ENDDO

      DO i=1,mx
       NSTgdx(i) = (REAL(i)-0.5)*NST_dx
      ENDDO

      DO j=1,my
       NSTgdy(j) = (REAL(j)-0.5)*NST_dx
      ENDDO


C +   ******************
      IF (Sinclair) THEN
C +   ******************


C +---Screen message
C +   --------------

      WRITE(6,*)  'Rain disaggregation: Sinclair (1994) model'
      WRITE(6,*)  '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF (input_filt)
     . WRITE(6,*) 'Filtering of LSC temp. and moisture fields'
      IF (output_filt)
     . WRITE(6,*) 'Filtering of disaggregated precipitation field'
      IF (local_conserv)
     . WRITE(6,*) 'Local conservation of precipitated water'
      IF (global_conserv)
     . WRITE(6,*) 'Global conservation of precipitated water'
      IF (deep_conv)
     . WRITE(6,*) 'Rough representation of deep convection'
      IF (surf_flux)
     . WRITE(6,*) 'Interactions with surface fluxes'
      WRITE(6,*)


C +---Initialisation
C +   --------------

      DO j=1,my
      DO i=1,mx
       NSTlpr(i,j) = 0.0
      ENDDO
      ENDDO

      IF (surf_flux) THEN
       switch_srf = 1.0
      ELSE
       switch_srf = 0.0
      ENDIF


C +---Compute inversion height
C +   ------------------------

      IF (namTKE.ne.'-') THEN

       DO j=1,my
       DO i=1,mx

        zinv(i,j)=0.
        TKEsrf=0.5*(NSTtke(i,j,mz)+NSTtke(i,j,mz-1))

        DO k=mz-1,2,-1
         IF (NSTtke(i,j,k).gt.TKEsrf.and.NSTtke(i,j,k-1).le.TKEsrf.and.
     .       zinv(i,j).lt.1.d-5) THEN
          IF ((NSTtke(i,j,k)-NSTtke(i,j,k-1)).gt.1.d-5) THEN
           zinv(i,j)=( (TKEsrf-NSTtke(i,j,k-1))*NST_zz(i,j,k  )
     .               + (NSTtke(i,j,k) - TKEsrf)*NST_zz(i,j,k-1) )
     .               / (NSTtke(i,j,k)-NSTtke(i,j,k-1))
     .              -   NST_sh(i,j)
          ELSE
           zinv(i,j)=0.5*(NST_zz(i,j,k)+NST_zz(i,j,k-1)) - NST_sh(i,j)
          ENDIF
         ENDIF
        ENDDO

        IF (zinv(i,j).lt.1.d-5) THEN
         zinv(i,j)=NST_zz(i,j,mz)-NST_sh(i,j)
        ENDIF
 
       ENDDO
       ENDDO


      ELSE      ! {IF (namTKE.ne.'-')}


C +---Compute Virtual and Equivalent Potential Temperature --> tetav
C +   ----------------------------------------------------     =====

       DO k=1,mz
       DO j=1,my
       DO i=1,mx
        tetav(i,j,k)=NST__t(i,j,k) * (100./NST__p(i,j,k))**cap
     .              *(1.+0.608*NST_qv(i,j,k)-NST_qt(i,j,k))
       ENDDO
       ENDDO
       ENDDO


C +---Compute Integrated Buoyancy --> buo
C +   ---------------------------     ===

       DO j=1,my
       DO i=1,mx

        int_buo =0.

        DO k=mz,2,-1
         int_buo =int_buo-(0.5*tetav(i,j,k -1)+0.5*tetav(i,j,k)
     .                    -0.5*tetav(i,j,mz-1)-0.5*tetav(i,j,mz-2))*9.81
     .                   /(tetav(i,j,k-1)+tetav(i,j,k))/2.
     .                   *(NST_zz(i,j,k-1)-NST_zz(i,j,k))
         IF (k.ge.(mz-1)) int_buo=MAX(0.,int_buo)
         buoy(i,j,k)=int_buo
        ENDDO

       ENDDO
       ENDDO


C +---Bulk Richardson Number --> rib
C +   ----------------------     ===

       DO j=1,my
       DO i=1,mx

        int_rib=0.
        int_u=0.

        DO k=mz,2,-1

         int_u=NST__u(i,j,k)*NST__u(i,j,k)
     .        +NST__v(i,j,k)*NST__v(i,j,k)
         int_u=MAX(int_u,eps3)

         int_rib=(NST_zz(i,j,k)-NST_sh(i,j))
     .          *(tetav(i,j,k)-tetav(i,j,mz))
     .          /tetav(i,j,k)*9.81
     .          /int_u

         rib(i,j,k)=int_rib

        ENDDO

       ENDDO
       ENDDO


C +---Boundary layer depth with Bulk Richardson method --> zrib
C +   ------------------------------------------------     ====

       DO j=1,my
       DO i=1,mx

        zrib(i,j)=0.0
        zbuo     =0.0

        DO k=mz-1,2,-1

         IF (buoy(i,j,k).gt.0.1      .and.
     .       buoy(i,j,k).ne.buoy(i,j,k-1)) THEN
          aux1=(buoy(i,j,k)-0.001)/(buoy(i,j,k)-buoy(i,j,k-1))
          IF (aux1.lt.0.) aux1=0.0
          IF (aux1.gt.1.) aux1=1.0
          aux2=1.0-aux1
          zbuo=aux1*NST_zz(i,j,k-1)+aux2*NST_zz(i,j,k)
     .        -NST_sh(i,j)
         ENDIF

         IF (zrib(i,j).eq.0.      .and.
     .      rib(i,j,k+1).lt.ribc .and.
     .      rib(i,j,k).gt.ribc) THEN
          aux1=(rib(i,j,k)-ribc)/(rib(i,j,k)-rib(i,j,k+1))
          IF (aux1.lt.0.) aux1=0.0
          IF (aux1.gt.1.) aux1=1.0
          aux2=1.0-aux1
          zrib(i,j)=aux1*NST_zz(i,j,k+1)+aux2*NST_zz(i,j,k)
     .             -NST_sh(i,j)
         ENDIF

        ENDDO

        IF (zbuo.gt.zrib(i,j)) zrib(i,j)=zbuo

ccccc   zinv(i,j) = MAX( 100.,zrib(i,j))
ccccc   zinv(i,j) = MIN(3000.,zrib(i,j))

        IF (zinv(i,j).lt.1.d-5) THEN
         zinv(i,j)=NST_zz(i,j,mz)-NST_sh(i,j)
        ENDIF
 
       ENDDO
       ENDDO


      ENDIF     ! {IF (namTKE.ne.'-')}


C +---Compute available and potential convective energy
C +   -------------------------------------------------

      DO j=1,my
      DO i=1,mx

       DO k=1,mz
        T_Envr(k)=NST__t(i,j,k)
        QvEnvr(k)=NST_qv(i,j,k)
        H_Envr(k)=NST_zz(i,j,k)-NST_sh(i,j)
        P_Envr(k)=NST__p(i,j,k)
       ENDDO

       kINI=mz
       hPBL=200.  !  Depth of the layer considered for parcel lifting

C +         ******
       CALL CVAdia (kINI,kLCL,hPBL,zLCL,zLFC,zETL,aBuoyE,pBuoyE,
     .              AdiabT,T_Envr,QvEnvr,H_Envr,P_Envr)
C +         ******

       TaBuoyE(i,j) = aBuoyE
       TpBuoyE(i,j) = pBuoyE
       TkLCL  (i,j) = kLCL
       TzLCL  (i,j) = zLCL
       TzLFC  (i,j) = zLFC
       TzETL  (i,j) = zETL
       DO k=1,mz
        TAdiabT(i,j,k) = AdiabT(k)
        TT_Envr(i,j,k) = T_Envr(k)
       ENDDO

      ENDDO
      ENDDO


C +---Mean loop on horizontal grid
C +   ----------------------------


      DO j=jc1,jc2  ! Loop on NST grid points
      DO i=ic1,ic2  ! -----------------------


C +---Initialisation
C +   --------------

       conv_adjust = .false.

       DO k=1,mz
        lev_OK(k)  = .false.
       ENDDO

       kLCL   = MIN(TkLCL(i,j),mz-1)
       aBuoyE = TaBuoyE(i,j)
       pBuoyE = TpBuoyE(i,j)
       zLCL   = TzLCL  (i,j)
       zLFC   = TzLFC  (i,j)
       zETL   = TzETL  (i,j)


C +---Compute relative humidity
C +   -------------------------

       DO k=1,mz
        qvsat (k) = qsat(NST__t(i,j,k),NST__p(i,j,k))
        r(k)      = NST_qv(i,j,k) / max(0.1e-7,1.-NST_qv(i,j,k))
        rst(k)    = qvsat (k)     / max(0.1e-7,1.-qvsat (k))
        relhum(k) = (r(k)/(0.622+r(k))) / max(0.1e-7,rst(k)
     .            / (0.622+rst(k))) 
        relhum(k) = min(1.,relhum(k))
       ENDDO


C +---Compute vertical component of velocity due to thermals
C +   ------------------------------------------------------

      wturb_max = 0.
      DO k=1,mz
       aux1   = NST_zz(i,j,k)-NST_sh(i,j)
       aux2   = aux1/zinv(i,j)
       aux2   = MIN(aux2,1.25)
       sflux  = MAX(rzero,NSTuts(i,j))/cp
       wwprim = 1.8 * (aux2)**(2./3.) * (1.0-0.8*aux2)**2.0
     .        * (gravit/NST__t(i,j,mz)*sflux*zinv(i,j))**(2./3.)
       wwprim = MAX(rzero,wwprim)
       wturb(k) = SQRT(wwprim)
       wturb_max= MAX(wturb_max,wturb(k))
      ENDDO


C +---Compute updraft velocity
C +   ------------------------

       w_conv = 2.0*MAX(rzero,pBuoyE)
       w_conv = SQRT(w_conv)

       DO k=mz,kLCL,-1
        w_Updr(k)=fact_w*NST__w(i,j,k)
       ENDDO
       DO k=kLCL+1,1,-1
        deltaz   =NST_zz(i,j,k)-NST_zz(i,j,k+1)
        w_Updr(k)=w_Updr(k+1)*w_Updr(k+1)
     .           +2.d0*deltaz*gravit*(TAdiabT(i,j,k)-TT_Envr(i,j,k))
     .                                      /(1.5d+0*TT_Envr(i,j,k))
        w_Updr(k)=SQRT(MAX(w_Updr(k),rzero))
       ENDDO


C +---Averaged variables in the boundary layer
C +   ----------------------------------------


C +... Depth
C +    - - -
       hgtmax = 1500. 
       hgtmax = MIN(hgtmax,NST_zz(i,j,kLCL)-NST_sh(i,j))
       hgtmax = MAX(hgtmax,NST_zz(i,j,mz-1)-NST_sh(i,j))

C +... Initialisation
C +    - - - - - - - -

       dztot = 0.
       Ubl   = 0.
       Vbl   = 0.
       Qbl   = 0.
       RHbl  = 0.


C +... Averaged horizontal wind and relative humidity
C +    - - - - - - - - - - - - - - - - - - - - - - - -

       DO k=2,mz 

        hgt1   = MIN(hgtmax,NST_zz(i,j,k-1)-NST_sh(i,j))
        hgt2   = MIN(hgtmax,NST_zz(i,j,k  )-NST_sh(i,j))
        dz     = hgt1-hgt2
        dztot  = dztot + dz

        Ubl    = Ubl  + dz*NST__u(i,j,k)
        Vbl    = Vbl  + dz*NST__v(i,j,k)

        Qbl    = Qbl  + dz*NST_qv(i,j,k)
        RHbl   = RHbl + dz*relhum(k)
       
       ENDDO

       Ubl  = Ubl  / dztot
       Vbl  = Vbl  / dztot
       Qbl  = Qbl  / dztot
       RHbl = RHbl / dztot

       Us_avg(i,j) = Ubl
       Vs_avg(i,j) = Vbl


C +---Slopes
C +   ------

C +...Distances on the sphere
       dst_dx = ACOS( sin(NST__y(i  ,j)*3.1415926/180.)*
     .                sin(NST__y(i-1,j)*3.1415926/180.)
     .              + cos(NST__y(i  ,j)*3.1415926/180.)*
     .                cos(NST__y(i-1,j)*3.1415926/180.)*
     .                cos((NST__x(i,j)-NST__x(i-1,j))*3.1415926/180.) )*
     .          6371229.
       dst_dy = ACOS( sin(NST__y(i,j  )*3.1415926/180.)*
     .                sin(NST__y(i,j-1)*3.1415926/180.)
     .              + cos(NST__y(i,j  )*3.1415926/180.)*
     .                cos(NST__y(i,j-1)*3.1415926/180.)*
     .                cos((NST__x(i,j)-NST__x(i,j-1))*3.1415926/180.) )*
     .          6371229.

       IF (dst_dx.lt.0.) write(90,*) i,j,NST__x(i-2,j),NST__x(i-1,j),
     .                                   NST__x(i  ,j),NST__x(i+1,j)
       IF (dst_dx.lt.0.) write(90,*) auxL2N(i-2,j),auxL2N(i-1,j),
     .                               auxL2N(i  ,j),auxL2N(i+1,j)
       dhx1(i,j) = (NST_sh(i,j)-NST_sh(i-1,j))/dst_dx  ! if U > 0
       dhx2(i,j) = (NST_sh(i+1,j)-NST_sh(i,j))/dst_dx  ! if U < 0
       dhy1(i,j) = (NST_sh(i,j)-NST_sh(i,j-1))/dst_dy  ! if V > 0
       dhy2(i,j) = (NST_sh(i,j+1)-NST_sh(i,j))/dst_dy  ! if V < 0


C +---Vertical velocity (pressure coordinates)
C +   ----------------------------------------

       Nfre =0
       do k =1,mz-1
        dz   = NST_zz(i,j,k) - NST_zz(i,j,k+1)
        dztot= dztot + dz
        Nfre = Nfre  + dz *
     .         sqrt(max(0.0,gravit*
     .        (log(NST_pt(i,j,k))-log(NST_pt(i,j,k+1)))/dz))
       end do
       Nfre=Nfre/dztot
      
       alpha  = ATAN (Vs_avg(i,j)/Us_avg(i,j))
       IF (Us_avg(i,j).lt.0.) alpha=alpha+3.1415927
       
       VgradS = ( MAX(0.,Us_avg(i,j)) * dhx1(i,j) 
     .          + MIN(0.,Us_avg(i,j)) * dhx2(i,j) ) * (COS(alpha))**2.
     .        + ( MAX(0.,Vs_avg(i,j)) * dhy1(i,j) 
     .          + MIN(0.,Vs_avg(i,j)) * dhy2(i,j) ) * (SIN(alpha))**2.

       omegaS  =-gravit*(factpr*NST__p(i,j,mz)/ra/NST__t(i,j,mz))
     .         * VgradS

       DO k=1,mz
        omega1 =-gravit*(factpr*NST__p(i,j,mz)/ra/NST__t(i,j,mz))
     .         * fact_w*NST__w(i,j,k)
        omegaT =-gravit*(factpr*NST__p(i,j,mz)/ra/NST__t(i,j,mz))
     .         * wturb(k)
        deltap = MAX(rzero,factpr*NST__p(i,j,k)-pr_top)

        if (vdelb_model) then 

        omegaP(k)=omega1+omegaS*sin(Nfre*(NST_sp(i,j)-NST__p(i,j,k))
     .           /sqrt(NST__u(i,j,mz)**2+NST__v(i,j,mz)**2))

c +...  VdelB model: http://www.geog.ucsb.edu/~chris/orthorain.pdf
c +...             : Funk et al. (2003), Int. Journal of Clim.,
c                    23:47-66.  
        
	else

        omegaP(k)=omega1+(omegaS+switch_srf*omegaT)
     .                   *(deltap/(factpr*NST_sp(i,j)-pr_top))
     .                    **TAN(gamma*3.1415926/4.)

C +...  Local component of omegaP(k) depends on gamma and pr_top

        endif

       ENDDO


C +---Check triggering of deep convection
C +   -----------------------------------

       conv_adjust = .false.

       IF (deep_conv) THEN

C +...Mean vertical velocity at the LCL = W
        wmLCL = fact_w*NST__w(i,j,kLCL) 
     .        + MAX(rzero,switch_srf*0.2*wturb_max)

C +...Stability at the LCL when a parcel is perturbated
        wLCL   = MAX(rzero,wmLCL)
c #ZF   stbLCL = TAdiabT(i,j,kLCL)-TT_Envr(i,j,kLCL) + wLCL**third

C +...aBuoyE = Available Buoyant Energy (between LCL and ETL)
C +...pBuoyE = Potential Buoyant Energy (between LFC and ETL)
C +...Negative Buoyancy between LCL and LFC (X 2)
        BEneg2 = MAX(pBuoyE-aBuoyE,rzero) * 2.0d+0

C +...Negative Velocity due to Negative Buoyancy between LCL and LFC
        wlow = SQRT(BEneg2)

C +...Estimated vertical velocity at the LFC
        wLFC = omegaP(kLCL)
     .        /(-gravit*(factpr*NST__p(i,j,mz)/ra/NST__t(i,j,mz))) 
     .       - wlow
     .       + switch_srf*0.2*wturb_max
ccccc   wLFC = (wmLCL+VgradS)-wlow

        IF ( wLFC.gt.0.            .and.
     .       wmLCL.gt.wthres       .and.
C +......... ^ Ability of a parcel to reach the LFC

c #ZF.       stbLCL.gt.0.          .and.
C +......... ^ Unstability at the LCL if perturbation

     .       zLCL.lt.(zETL-2000.)  .and.
     .       zLFC.lt.(zETL- 500.)  .and.
     .       zLFC.lt.      5000.   .and.
C +......... ^ Realizibility criterion based on LCL,LFC,ETL heights

     .       zETL.gt.      3000.   .and.
C +......... ^ Eliminate shallow convection --> only deep convection

     .       aBuoyE.gt.     300.       ) THEN
C +......... ^ Positive available buoyant energy (minimum threshold)

         conv_adjust = .true.

        ENDIF

       ENDIF


C +---Rough representation of deep convection
C +   ---------------------------------------

       IF (conv_adjust) THEN

C +...  Assess time of efficient convection

        aux1   = 0.0
        avgU   = 0.0
        avgV   = 0.0
        avwind = 0.0
        DO k=2,mz-1
         aux = NST_zz(i,j,k)-NST_sh(i,j)
         IF (aux.lt.zETL) THEN
          avgU = avgU + 0.5*(NST_zz(i,j,k-1)-NST_zz(i,j,k+1))
     .                     * NST__u(i,j,k)
          avgV = avgV + 0.5*(NST_zz(i,j,k-1)-NST_zz(i,j,k+1))
     .                     * NST__v(i,j,k)
          aux1 = aux1 + 0.5*(NST_zz(i,j,k-1)-NST_zz(i,j,k+1))
         ENDIF
        ENDDO
        IF (aux1.gt.0.1) THEN
         avgU   = avgU / aux1
         avgV   = avgV / aux1
        ELSE
         avgU   = NST__u(i,j,mz/2)
         avgV   = NST__v(i,j,mz/2)
        ENDIF
        avwind = SQRT(avgU*avgU+avgV*avgV)

ccccc   nx = NINT(0.75*avwind * (REAL(DAT_dt)*3600.) / NST_dx)
        nx = NINT(0.5*avwind * (REAL(DAT_dt)*3600.) / NST_dx)

        i1 = MAX( 1,i-nx/2)
        i2 = MIN(mx,i+nx/2)
        j1 = MAX( 1,j-nx/2)
        j2 = MIN(my,j+nx/2)
        IF (avgU.ge.0.0) THEN
         i1=i
        ELSE
         i2=i
        ENDIF
        IF (avgV.ge.0.0) THEN
         j1=j
        ELSE
         j2=j
        ENDIF

        compt  = 0.0
        buotot = 0.0

        DO jj=j1,j2
        DO ii=i1,i2
         IF (TaBuoyE(ii,jj).gt.300.) THEN
          compt  = compt  + 1.0
          buotot = buotot + TaBuoyE(ii,jj)
         ENDIF
        ENDDO
        ENDDO

        IF (compt.gt.0.5) THEN
         TIMcor_deep = compt / REAL((i2-i1+1)*(j2-j1+1))
     .               * (buotot/compt) / 3000.
        ELSE
         TIMcor_deep = 0.05
        ENDIF
        TIMcor_deep = MIN(1.00,TIMcor_deep)
        TIMcor_deep = MAX(0.05,TIMcor_deep)

C +...  Modified vertical velocities

        w__CNV = w_conv
        qvsCNV = NST_qv(i,j,kLCL)
        rhoCNV = factpr*NST__p(i,j,kLCL)/Ra/NST__t(i,j,kLCL)

        DO k=mz,1,-1
         aux = NST_zz(i,j,k)-NST_sh(i,j)
         IF (aux.le.(zLCL+1500.)) THEN
          w__CNV = w_Updr(k)
          qvsCNV = qvsat (k)
          rhoCNV = factpr*NST__p(i,j,k)/Ra/NST__t(i,j,k)
         ENDIF
        ENDDO

C +...  Characteristic time scale for deep convection

ccccc   IF (w_conv.gt.0.0001) THEN
        IF (w__CNV.gt.0.0001) THEN
ccccc    timecv = (zETL-zLCL) / w_conv
         timecv = MIN(zETL,5500.) / w__CNV
C +...   Starting point for downdraft
         timecv = MAX(rzero,timecv)
        ELSE
         timecv = 0.0
        ENDIF
        timecv = MIN(timecv,REAL(DAT_dt)*3600.)
        timecv = MIN(timecv,1800.)               ! Max conv. adj. time
        dur_dt = timecv /  (REAL(DAT_dt)*3600.) * TIMcor_deep

C +...  Compute condensation rate 150 mb above LCL (Fritsch and Chappell, 1980)

        diffqv   = MAX(rzero,Qbl-qvsCNV)
        rate_CNV = rhoCNV * diffqv * w__CNV
     .           * dur_dt * 0.2

       ENDIF


C +---Condensation rates (Haltiner and Williams, 1980)
C +   ------------------------------------------------

       DO k=1,mz


        rhmin    =0.8 ! 0.6  

        valhum   =MAX(rhmin,relhum(k))
        lambda(k)=sqrt((valhum-rhmin)/(rhcrHY-rhmin))
        Ffact(k) = rst(k) *NST__t(i,j,k)
     .           / (factpr*NST__p(i,j,k))
     .           * (h2olv*ra-cp*rv*NST__t(i,j,k))
     .           / (cp*rv*NST__t(i,j,k)**2.+rst(k)*h2olv*h2olv) 

C +.... Condensation rate = lambda(k) * Ffact(k) * omegaP(k)
C +.... if omegaP(k) < 0

       ENDDO


C +---Efficiency factor of orographic precipitation 
C +   (function of low level RH) ------------------
C +   --------------------------

       rhmin =0.8
       valhum=MAX(rhmin,RHbl)
       lambds=((valhum-rhmin)/(rhcrHY-rhmin))**(0.25)

C +---Precipitation efficiency for deep convection
C +   depending on Cloud Bottom ------------------
C +   -------------------------

       IF (conv_adjust) THEN

        zLCLkf=zLCL/(1000.*0.3048)
C +...  ^^^^ LCL height in units of thousands of feet

        auxEff=0.967-0.700*zLCLkf+0.162*zLCLkf*zLCLkf
     .                  -0.01257*zLCLkf*zLCLkf*zLCLkf

        Effic2=1./(1.+auxEff)
C +...  ^^^^ Precipitation efficiency based on LCL zLCLkf
C +          (Zhang and Fritsch, 1986)

        Effic2=max(0.0,Effic2)
        Effic2=min(1.0,Effic2)

       ENDIF


C +---Precipitation rate
C +   ------------------

       IF (conv_adjust) THEN

        DO k=1,mz
         rateqv(k) = 0.0     ! no orographic precipitation
         aux = NST_zz(i,j,k)-NST_sh(i,j)
         IF (aux.ge.zLCL.and.aux.le.5500.) THEN
          rateDC(k) = lambds * Effic2 * rate_CNV
         ELSE
          rateDC(k) = 0.0
         ENDIF
        ENDDO

       ELSE

        DO k=1,mz

         rateqv(k) = lambds*(-omegaP(k))*Ffact(k)   *lambda(k)
c        rateqv(k) = lambds*(-omegaP(k))*Ffact(kLCL)*lambda(k)
 
C +...   Rem.: Ffact(kLCL) = this means that air ascent properties
C +...         corresponds to the characteristics of LCL
         rateDC(k) = 0.0     ! no deep convection
        ENDDO

       ENDIF


C +---Falling speed of precipitating hydrometeors
C +   -------------------------------------------

       pfreeze    = 85.0
       kfreeze    = 1
       IF (NST__t(i,j,mz).lt.tfreeze) pfreeze = NST__p(i,j,mz)

       DO k=1,mz-1
        IF (NST__t(i,j,k).le.tfreeze) THEN
         vfall(k) = vsnow
         kfreeze  = max(kfreeze,k)
         IF (NST__t(i,j,k+1).gt.tfreeze) pfreeze = NST__p(i,j,k+1)
        ELSE
         vfall(k) = vrain 
        ENDIF
       ENDDO
       vfall(mz)  = vrain
   
       DO k=1,mz-1
        r(k)      = NST_qv(i,j,k) / max(0.1e-7,1.-NST_qv(i,j,k))
        rstsrf    = qsat(NST__t(i,j,mz),NST__p(i,j,mz))     
     .            / max(0.1e-7,1.-qsat(NST__t(i,j,mz),NST__p(i,j,mz)))
        rstfreeze = qsat(NST__t(i,j,kfreeze),NST__p(i,j,kfreeze))     
     .            / max(0.1e-7,
     .              1.-qsat(NST__t(i,j,kfreeze),NST__p(i,j,kfreeze)))
        if (k.gt.kfreeze) then
         vfall(k) = 6 * r(k) / rstsrf
        else
         vfall(k) = 1 * r(k) / rstfreeze
        end if
       END DO     
   
       IF (conv_adjust) THEN
       DO k=1,mz
        vfall(k) = 3.0 * vfall(k)
       ENDDO
       ENDIF


C +---Formation time of precipitation (Sinclair, 1994)
C +   ------------------------------------------------

       DO k=1,mz
        timef(k) = 1000. * (0.5 + 1.0/3.1415927
     .                     *ATAN((pfreeze-NST__p(i,j,k))/5.0))
       ENDDO


C +---Horizontal Lagragian advection
C +   ------------------------------

C +... Time step definition
C +    - - - - - - - - - - -

       dtSlow= 300.0
       max_V = 0.00001d0

       DO k=1,mz
        max_V = MAX(max_V,abs(NST__u(i,j,k)))
        max_V = MAX(max_V,abs(NST__v(i,j,k)))
       ENDDO
       
       dt_max = INT(0.5 * MIN(dst_dx,dst_dy) / ABS(max_V))
       
       IF (dt_max.eq.0) THEN
        dt_max = 1
        write(6,*) 'dt_max 0 vers 1'
       ENDIF
       
       dtSlow = MIN(dtSlow,dt_max)


C +... Loop on vertical levels above the LCL 
C +... and for levels where condensed occured
C +    - - - - - - - - - - - - - - - - - - - -

       DO k=2,kLCL
       IF (rateqv(k).gt.0..or.rateDC(k).gt.0.) THEN


C +...  Initial location of precipitation
C +     - - - - - - - - - - - - - - - - -
        
        xloc = NSTgdx(i)    ! On staggered grid
        yloc = NSTgdy(j)
        iloc = i
        jloc = j


C +...  Loop until hydrometeors reach the surface
C +     - - - - - - - - - - - - - - - - - - - - -

        itloc = 0

        DO WHILE (.not.lev_OK(k))

         itloc = itloc + 1
         tloc  = REAL(itloc)*dtSlow


C +...   Horizontal linear interpolation of wind field
C +      - - - - - - - - - - - - - - - - - - - - - - -

         aux1 = (xloc - (REAL(iloc)-1.0)*NST_dx)  ! Non-staggered
     .        / NST_dx
         aux1 = MAX(rzero,aux1)
         aux1 = MIN(1.0,aux1)
         aux2 = (yloc - (REAL(jloc)-1.0)*NST_dx)
     .        / NST_dx
         aux2 = MAX(rzero,aux2)
         aux2 = MIN(1.0,aux2)

         uloc = (1.0-aux2) * ((1.0-aux1)*NST__u(i  ,j  ,k)
     .                       + aux1     *NST__u(i+1,j  ,k))
     .        + aux2       * ((1.0-aux1)*NST__u(i  ,j+1,k)
     .                       + aux1     *NST__u(i+1,j+1,k))

         vloc = (1.0-aux2) * ((1.0-aux1)*NST__v(i  ,j  ,k)
     .                       + aux1     *NST__v(i+1,j  ,k))
     .        + aux2       * ((1.0-aux1)*NST__v(i  ,j+1,k)
     .                       + aux1     *NST__v(i+1,j+1,k))

        
C +...   Advection of the location of precipitation center
C +      - - - - - - - - - - - - - - - - - - - - - - - - -

         xloc = xloc + uloc*dtSlow
         yloc = yloc + vloc*dtSlow


C +...   New location of (xloc,yloc) in the NST grid
C +      - - - - - - - - - - - - - - - - - - - - - -

         iloc = INT((xloc+0.5*NST_dx)/NST_dx)
         jloc = INT((yloc+0.5*NST_dx)/NST_dx)

        
C +...  Fall time of hydrometeors to reach the surface
C +     - - - - - - - - - - - - - - - - - - - - - - - -

         IF (tloc.gt.timef(k)) THEN

          tfall = 0.0
          k1    = MIN(mz-1,k)

          DO kk=k1,mz-1
           tfall = tfall 
     .           + (NST_zz(i,j,kk)-NST_zz(i,j,kk+1))
     .            /(0.5*(vfall(kk)+vfall(kk+1)))
          ENDDO


C +...    IF the hydrometeors reach the surface ...
C +       - - - - - - - - - - - - - - - - - - - - -

          IF (tloc.gt.(tfall+timef(k))) THEN


C +...     Distribution of precipitation over adjacent grid cells
C +        - - - - - - - - - - - - - - - - - - - - - - - - - - - -

           i1 = iloc - 1
           i1 = MAX( 1,i1)
           i2 = iloc + 1
           i2 = MIN(mx,i2)
           j1 = jloc - 1
           j1 = MAX( 1,j1)
           j2 = jloc + 1
           j2 = MIN(my,j2)

           DO jj=j1,j2
           DO ii=i1,i2
            dxic = ABS(xloc-NSTgdx(ii))
            dyic = ABS(yloc-NSTgdy(jj))
            dxic = MIN(dxic,NST_dx)
            dyic = MIN(dyic,NST_dx)
            area = NST_dx*NST_dx
     .           - (dyic*NST_dx+dxic*NST_dx-dxic*dyic)
            tmp_2D(ii,jj) = MAX(rzero,area) / (NST_dx*NST_dx)
           ENDDO
           ENDDO


C +...     Precipitation rates due to hydrometeor fall at level k
C +        - - - - - - - - - - - - - - - - - - - - - - - - - - - -

           deltap = factpr * 0.5 * (NST__p(i,j,k+1)-NST__p(i,j,k-1))
           val_pr = ( rateqv(k) * deltap / gravit 
     .              + rateDC(k) ) 
     .            * (REAL(DAT_dt)*3600.)
C +....... Rem.: -deltap * gravit = deltaz * rhoair

           DO jj=j1,j2
           DO ii=i1,i2
            NSTlpr(ii,jj) = NSTlpr(ii,jj) + tmp_2D(ii,jj)*0.001*val_pr
C +........ Rem.: 0.001 -> kg(water)/m2 converted in m of water (height)
           ENDDO
           ENDDO
           
           lev_OK(k) = .true.


          ELSE  ! tloc.gt.(tfall+timef(k))


C +...     Exit condition if hydrometeors takes too much time to reach
C +...     the surface - - - - - - - - - - - - - - - - - - - - - - - -
C +        - - - - - -

           IF (tfall.gt.6000.)  lev_OK(k) = .true.


          ENDIF ! tloc.gt.(tfall+timef(k))


         ENDIF  ! tloc.gt.timef(k)
       

        ENDDO   ! WHILE loop: until (.not.lev_OK(k))


       ENDIF    ! rateqv(k).gt.0.L
       ENDDO    ! Loop on k=1,kLCL


      ENDDO     ! Loop on i=2,mx-1
      ENDDO     ! Loop on j=2,my-1


C +   *****
      ENDIF ! Sinclair
C +   *****


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Disaggregation model of precipitation (Alpert, 1989)
C +   ====================================================


C +   ****************
      IF (Alpert) THEN
C +   ****************

      
C +---Screen message
C +   --------------

      WRITE(6,*) 'Rain disaggregation: Alpert (1989) model'
      WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF (input_filt)
     . WRITE(6,*) 'Filtering of LSC temp. and moisture fields'
      IF (output_filt)
     . WRITE(6,*) 'Filtering of disaggregated precipitation field'
      IF (local_conserv)
     . WRITE(6,*) 'Local conservation of precipitated water'
      IF (global_conserv)
     . WRITE(6,*) 'Global conservation of precipitated water'
      WRITE(6,*)


      DO j=jc1,jc2
      DO i=ic1,ic2


C +---Mean wind in the boundary layer
C +   -------------------------------

       dztot = 0.
       Ubl   = 0.
       Vbl   = 0.
       Wbl   = 0.
       DO k=1,mz-1
        height = NST_zz(i,j,k)-NST_sh(i,j)
        IF (height.lt.1500) THEN
         dz    = NST_zz(i,j,k)-NST_zz(i,j,k+1)
         dztot = dztot + dz
         Ubl   = Ubl + NST__u(i,j,k)*dz
         Vbl   = Vbl + NST__v(i,j,k)*dz
         Wbl   = Wbl + NST__w(i,j,k)*dz*fact_w
        ENDIF
       ENDDO
       Ubl = Ubl / dztot
       Vbl = Vbl / dztot
       Wbl = Wbl / dztot
       Us_avg(i,j) = Ubl
       Vs_avg(i,j) = Vbl
       Ws_avg(i,j) = Wbl


C +---Slopes
C +   ------

       dst_dx=ABS(NST__x(i,j)-NST__x(i-1,j))
     .       *40000000./360.*ABS(COS(3.1415926*NST__y(i,j)/180.))
       dst_dy=ABS(NST__y(i,j)-NST__y(i,j-1))
     .       *40000000./360.

       IF (Us_avg(i,j).ge.0.) THEN
        NSTdhx = (NST_sh(i,j)-NST_sh(i-1,j))/dst_dx
        LSCdhx = (INT_sh(i,j)-INT_sh(i-1,j))/dst_dx
       ELSE
        NSTdhx = (NST_sh(i+1,j)-NST_sh(i,j))/dst_dx
        LSCdhx = (INT_sh(i+1,j)-INT_sh(i,j))/dst_dx
       ENDIF

       IF (Vs_avg(i,j).ge.0.) THEN
        NSTdhy = (NST_sh(i,j)-NST_sh(i,j-1))/dst_dy
        LSCdhy = (INT_sh(i,j)-INT_sh(i,j-1))/dst_dy
       ELSE
        NSTdhy = (NST_sh(i,j+1)-NST_sh(i,j))/dst_dy
        LSCdhy = (INT_sh(i,j+1)-INT_sh(i,j))/dst_dy
       ENDIF

C +---Local vertical velocity
C +   -----------------------

       alpha  = ATAN (Vs_avg(i,j)/Us_avg(i,j))
       IF (Us_avg(i,j).lt.0.) alpha=alpha+3.1415927
       VgradS=Us_avg(i,j)*NSTdhx*(COS(alpha))**2.
     .       +Vs_avg(i,j)*NSTdhy*(SIN(alpha))**2.


C +---Total vertical velocity
C +   -----------------------

       Wtot = Vgrads + Wbl


C +---Precipitation rate
C +   ------------------

       k           = mz-1
       qvsat(k)    = qsat(NST__t(i,j,k),NST__p(i,j,k))
       relhum(k)   = NST_qv(i,j,k) / qvsat(k)
       ezs         = esat(NST__t(i,j,k),NST__p(i,j,k))
       tzs         = NST__t(i,j,k)
       ratepr(i,j) = 0.622*relhum(k)*ezs*Wtot/ra/tzs

      ENDDO
      ENDDO


C +---Weight function
C +   ---------------

      REGnst = .false.
      STEnst = .true.

      DO j=jc1,jc2
      DO i=ic1,ic2

       x0lon = NST__x(i,j)
       y0lat = NST__y(i,j)

       xilon = x0lon
       yilat = y0lat

       DO l=1,nstep

C +          ******
        CALL INThor1 (HORint,NST__x,NST__y,Us_avg,
     .                STEnst,xilon ,yilat ,Uint  ,
     .                REGnst,possOx,possOy)
        CALL INThor1 (HORint,NST__x,NST__y,Vs_avg,
     .                STEnst,xilon ,yilat ,Vint  ,
     .                REGnst,possOx,possOy)
C +          ******

        IF (l.eq.1) sigma = sqrt(Us_avg(i,j)*Uint
     .                          +Vs_avg(i,j)*Vint) * lifetime

        xilon  = xilon - Uint * lifetime / REAL(nstep)
     .                   /(111111.*COS(yilat*degrad))
        yilat  = yilat - Vint * lifetime / REAL(nstep)
     .                   / 111111.

C +          *******
        CALL INThor1 (HORint,NST__x,NST__y,ratepr,
     .                STEnst,xilon ,yilat ,valint,
     .                REGnst,possOx,possOy)
C +          *******

        dist = SQRT ( ((xilon-x0lon)*111111.*COS(y0lat*degrad))**2.0
     .              + ((yilat-y0lat)*111111.                  )**2.0 )

        wfct  (l) = EXP(-(dist*dist)/(2.0*sigma*sigma))
        precip(l) = valint
      
       ENDDO

       tot_wf = 0.
       sumtmp = 0.
       DO l=1,nstep
        sumtmp = sumtmp + wfct(l)*precip(l)
        tot_wf = tot_wf + wfct(l)
       ENDDO

       efficiency  = 0.5

       NSTlpr(i,j) = efficiency*sumtmp / tot_wf * REAL(DAT_dt) * 3600.
C +                  ^^^^^^^^^^^^^^^^^^^^^^^^^^ precipitation rate
       

      ENDDO
      ENDDO


C +   *****
      ENDIF ! Alpert
C +   *****


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Horizontal filtering of noise in precipitation field
C +   ====================================================

      IF (output_filt) THEN
      
       eps = 0.20    ! Filter Selectivity Parameter

C +         *********
       CALL DYNfil_2H (NSTlpr,eps)
C +         *********

      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Horizontal advection (in order to take into account the time
C +                         for precipitation processes)
C +   ============================================================

ccccc dt    = 1500.
ccccc max_V = 0.d0

ccccc DO j=1,my
ccccc DO i=1,mx
ccccc  max_V = MAX(max_V,Us_avg(i,j))
ccccc  max_V = MAX(max_V,Vs_avg(i,j))
ccccc ENDDO
ccccc ENDDO

ccccc ntSlow = INT(0.5 * MIN(dst_dx,dst_dy) / max_V) + 1
ccccc dtSlow = dt / REAL(ntSlow)

C +        ******
ccccc CALL ADV_2D (ntSlow,dtSlow,dt,dst_dx,dst_dy,Us_avg,Vs_avg,NSTlpr)
C +        ******


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Store non-conserved water fields
C +   ================================

      DO j=1,my
      DO i=1,mx
       NSTpr1(i,j)=NSTlpr(i,j)
      ENDDO
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Water conservation over 4*4 LSC cells
C +   =====================================


      IF (local_conserv) THEN

       idiv = NINT(4.0*LSC_dx/NST_dx)
       jdiv = NINT(4.0*LSC_dy/NST_dx)


C +    **************
       DO ngrid = 1,2
C +    **************


C +---1. Initialisation
C +   -----------------

        DO j=1,my
        DO i=1,mx
         dstLB (i,j,ngrid) = 0.0
         TMP_pr(i,j,ngrid) = 0.0
        ENDDO
        ENDDO


C +---2. Non-staggered (ngrid=1) and staggered (ngrid=2) grids
C +   --------------------------------------------------------

        istag = (ngrid-1)*idiv/2
        jstag = (ngrid-1)*jdiv/2

        DO l=1,(my/jdiv)+ngrid
        DO k=1,(mx/idiv)+ngrid

         nbr_pr = 0.0
         pr_tot = 0.0
         prItot = 0.0

         ii1 = 1+(k-1)*idiv - istag
         ii2 = k*idiv       - istag
         jj1 = 1+(l-1)*jdiv - jstag
         jj2 = l*jdiv       - jstag
         ii1 = MIN(mx,ii1)
         ii1 = MAX( 1,ii1)
         ii2 = MIN(mx,ii2)
         ii2 = MAX( 1,ii2)
         jj1 = MAX( 1,jj1)
         jj1 = MIN(my,jj1)
         jj2 = MAX( 1,jj2)
         jj2 = MIN(my,jj2)

         xcent = 0.25 * (NST__x(ii1,jj1)+NST__x(ii1,jj2)
     .                  +NST__x(ii2,jj1)+NST__x(ii2,jj2))
         ycent = 0.25 * (NST__y(ii1,jj1)+NST__y(ii1,jj2)
     .                  +NST__y(ii2,jj1)+NST__y(ii2,jj2))

         IF (ii1.ne.ii2.and.jj1.ne.jj2) THEN
          DO j=jj1,jj2
          DO i=ii1,ii2
           pr_tot = pr_tot + NSTlpr(i,j)
           prItot = prItot + NSTIpr(i,j)
           nbr_pr = nbr_pr + 1.0
           coslat = ABS(COS(3.1415926*NST__y(i,j)/180.))
           dstLB(i,j,ngrid) = 
     .               SQRT( ((NST__x(i,j)-xcent)*111.111*coslat)**2.0
     .                   + ((NST__y(i,j)-ycent)*111.111       )**2.0)
          ENDDO
          ENDDO
         ENDIF

C +...   Disaggregated precipitation = large-scale precipitation
          IF (pr_tot.le.0.00001) THEN

          DO j=jj1,jj2
          DO i=ii1,ii2
           TMP_pr(i,j,ngrid) = NSTIpr(i,j)
          ENDDO
          ENDDO

         ELSE

          ratio = prItot/pr_tot
 
C +...   Reduce desagregated precipitation
          IF (ratio.le.1.0) THEN
           DO j=jj1,jj2
           DO i=ii1,ii2
            TMP_pr(i,j,ngrid) = ratio * NSTlpr(i,j)
           ENDDO
           ENDDO
C +...   Increase desagregated precipitation with background
C +...   large-scale precipitation
          ELSE
           IF (prItot.gt.0.00001) THEN
            DO j=1,my
            DO i=1,mx
             TMP_pr(i,j,ngrid) = NSTlpr(i,j) + (prItot-pr_tot)/prItot
     .                                                   *NSTIpr(i,j)
            ENDDO
            ENDDO
           ELSE
            DO j=1,my
            DO i=1,mx
             TMP_pr(i,j,ngrid) = NSTlpr(i,j) + (prItot-pr_tot)/nbr_pr
            ENDDO
            ENDDO
           ENDIF
          ENDIF
 
         ENDIF


        ENDDO
        ENDDO


C +    ************************
       ENDDO   !  {ngrid = 1,2}
C +    ************************


C +---3. Combine the fields of staggered and non-staggered grids
C +   ----------------------------------------------------------

       DO j=1,my
       DO i=1,mx
        aux = dstLB(i,j,1)**2.0 + dstLB(i,j,2)**2.0
        IF (aux.gt.0.0001) THEN
         NSTlpr(i,j) = (dstLB(i,j,2)**2.0/aux) * TMP_pr(i,j,1)
     .               + (dstLB(i,j,1)**2.0/aux) * TMP_pr(i,j,2)
        ELSE
         NSTlpr(i,j) = 0.5 * (TMP_pr(i,j,1)+TMP_pr(i,j,2))
        ENDIF
       ENDDO
       ENDDO


      ENDIF  !  {local_conserv}


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---New distribution of precipitation
C +   =================================

C +---1. Fields without local conservation
C +   ------------------------------------

      prItot=0.
      pr_tot=0.

      DO j=1,my
      DO i=1,mx
       prItot=prItot+NSTIpr(i,j)
       pr_tot=pr_tot+NSTpr1(i,j)
      ENDDO
      ENDDO

      prIavg=prItot/REAL(mx)/REAL(my)
      pr_avg=pr_tot/REAL(mx)/REAL(my)

      IF (pr_avg.gt.prIavg) THEN
       ratio=prIavg/pr_avg
      ELSE
       ratio=1.
      ENDIF

      pr_tot=0.

      write(6,*) 'Correction for global conservation : ',ratio
      write(6,*)
      DO j=1,my
      DO i=1,mx
       NST_pr(i,j)=NSTIpr(i,j)+ratio*(NSTpr1(i,j)-pr_avg)
       NST_pr(i,j)=MAX(rzero,NST_pr(i,j))
       pr_tot     =pr_tot+NST_pr(i,j)
      ENDDO
      ENDDO

      IF (pr_tot.le.0.00001) THEN
       DO j=1,my
       DO i=1,mx
        NST_pr(i,j)=NSTIpr(i,j)
        NSTpr2(i,j)=NST_pr(i,j)
       ENDDO
       ENDDO
      ELSE
       ratio=prItot/pr_tot
       DO j=1,my
       DO i=1,mx
        NST_pr(i,j)=NST_pr(i,j)*ratio
        NSTpr2(i,j)=NST_pr(i,j)
       ENDDO
       ENDDO
      ENDIF


C +---2. Fields with local conservation
C +   ---------------------------------

      prItot=0.
      pr_tot=0.

      DO j=1,my
      DO i=1,mx
       prItot=prItot+NSTIpr(i,j)
       pr_tot=pr_tot+NSTlpr(i,j)
      ENDDO
      ENDDO

      prIavg=prItot/REAL(mx)/REAL(my)
      pr_avg=pr_tot/REAL(mx)/REAL(my)

      IF (pr_avg.gt.prIavg) THEN
       ratio=prIavg/pr_avg
      ELSE
       ratio=1.
      ENDIF

      pr_tot=0.

      DO j=1,my
      DO i=1,mx
       NST_pr(i,j)=NSTIpr(i,j)+ratio*(NSTlpr(i,j)-pr_avg)
       NST_pr(i,j)=MAX(rzero,NST_pr(i,j))
       pr_tot     =pr_tot+NST_pr(i,j)
      ENDDO
      ENDDO

      IF (pr_tot.le.0.00001) THEN
       DO j=1,my
       DO i=1,mx
        NST_pr(i,j)=NSTIpr(i,j)
        NSTpr3(i,j)=NST_pr(i,j)
       ENDDO
       ENDDO
      ELSE
       ratio=prItot/pr_tot
       DO j=1,my
       DO i=1,mx
        NST_pr(i,j)=NST_pr(i,j)*ratio
        NSTpr3(i,j)=NST_pr(i,j)
       ENDDO
       ENDDO
      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END


C   +-------------------------------------------------------------------+
C   |  Subroutine CVAdia                             June 2000  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | INPUT  : T_Envr: vertical sounding of (real) Temperature      (K) |
C   | ^^^^^^^^ QvEnvr: vertical sounding of specific humidity   (kg/kg) |
C   |          H_Envr: height of model levels for this sounding     (m) |
C   |                                                                   |
C   | OUTPUT : kLCL  : lifting condensation level                       |
C   | ^^^^^^^^ zLCL  : height of the lifting condensation level     (m) |
C   |          zLFC  : height of the level of free convection       (m) |
C   |          zETL  : height of the equilibrium Temperature level  (m) |
C   |          aBuoyE: available buoyant energy                 (m2/s2) |
C   |          pBuoyE: potential buoyant energy                 (m2/s2) |
C   |          AdiabT: adiabatic Temperature                            |
C   |                                                                   |
C   | REFER. : K. Emanuel, 'Atmospheric Convection', 1994.              |
C   | ^^^^^^^^                                                          |
C   |                                                                   |
C   +-------------------------------------------------------------------+
C +
C +
      subroutine CVAdia(kINI,kLCL,hPBL,zLCL,zLFC,zETL,aBuoyE,pBuoyE,
     .                  AdiabT,T_Envr,QvEnvr,H_Envr,P_Envr)
C +
C +
      IMPLICIT NONE
C +
C +
C +---General Variables
C +   =================
C +
      include 'NSTdim.inc'
C +
C +
C +---Local   Variables
C +   =================
C +
      integer k,kINI,kLCL,int_3
C +
      real T_Envr(mz),QvEnvr(mz),H_Envr(mz),AdiabT(mz),AdiabQ(mz),
     .     TvEnvr(mz),AdiavT(mz),P_Envr(mz)
C +
      real aux,tLCL,zLCL,zETL,zLFC,rw_Sat,rLiqid,T__Sat,qsat,Qv_Sat,
     .     deltaz1,deltaz2,dEnCin,aBuoyE,pBuoyE,T__PBL,Qv_PBL,H__PBL,
     .     Rd,Rv,Cpd,Cpv,Cl,r__PBL,dT__dz,Height,delz,diffqv,hPBL,
     .     gravit,h2olv,zero
C +
C +
C +---Constants
C +   =========
C +
      data Rd      / 287.04    /
      data Rv      / 461.5     /
      data Cpd     / 1005.7    /
      data Cpv     / 1870.     /
      data Cl      / 4190.     /
      data zero    /  0.       /
      data int_3   /  3        /
      data gravit  / 9.81      /
      data h2olv   / 2500000.  /
C +
C +
C +---Computation of Adiabatic Temperature Profile
C +   ============================================
C +
C +
C +---Mean Temperature and moisture of the air in the lower levels
C +   ------------------------------------------------------------
C +
      T__PBL=0.
      Qv_PBL=0.
      H__PBL=0.
      Height=0.
C +
      DO k=2,kINI
C +
        IF (H_Envr(k) .lt.  hPBL) THEN
C +...  ^^^^ Lower levels = 600 m from the surface
          delz  =H_Envr(k-1)-H_Envr(k)
          Height=Height+          delz
          T__PBL=T__PBL+T_Envr(k)*delz
          Qv_PBL=Qv_PBL+QvEnvr(k)*delz
          H__PBL=H__PBL+H_Envr(k)*delz
        ENDIF
C +
      ENDDO
C +
      IF (Height.gt.0.d0) THEN
        T__PBL=T__PBL/Height
        H__PBL=H__PBL/Height
        Qv_PBL=Qv_PBL/Height
C +...  ^^^^ Specific humidity
      ELSE
        T__PBL=T_Envr(kINI)
        Qv_PBL=QvEnvr(kINI)
        H__PBL=H_Envr(kINI)
      ENDIF
      r__PBL=Qv_PBL/(1.-Qv_PBL)
C +...^^^^ Mixing ratio
C +
C +---Lift a parcel from bottom to top
C +   --------------------------------
C +
      kLCL=0
C +
C +---Vertical Gradient of Temperature : Unsaturated Case
C +   ---------------------------------------------------
C +
C +...Dry adiabatic Temperature gradient (K. Emanuel, 1994)
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dT__dz=-(gravit/Cpd)*(1.+r__PBL)/(1.+r__PBL*Cpv/Cpd)
      DO k=mz,kINI    ,-1
        AdiabT(k)=T__PBL+ dT__dz*(H_Envr(k)-H__PBL-H_Envr(mz))
      ENDDO
C +
      DO k=kINI-1,1,-1
C +
        Qv_Sat    =  qsat(AdiabT(k+1),P_Envr(k+1))
        AdiabQ(k+1) = min(AdiabQ(k+1),Qv_Sat)
C +
        IF (kLCL.eq.0.and.k.lt.(mz-1).and.Qv_PBL.ge.Qv_Sat) THEN
          kLCL  =        k+1           !  LCL detected
          tLCL  = AdiabT(k+1)          !  ------------
          zLCL  = H_Envr(k+1)
        ENDIF
C +
        IF (kLCL.ne.0) THEN            ! Above LCL
C +                                    ! ---------
C +.... Moist adiabatic Temperature gradient (K. Emanuel, 1994)
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          T__Sat = AdiabT(k+1)
          rw_Sat = Qv_Sat/(1.-Qv_Sat)
          diffqv = r__PBL-rw_Sat
          rLiqid = max(zero,diffqv)
          dT__dz =-(gravit/Cpd)
     .           *(1.+            r__PBL)/ (1.0+rw_Sat*Cpv/Cpd)
     .           *(1.+h2olv/Rd   *rw_Sat / T__Sat)
     .           /(1.+            rLiqid * Cl/(Cpd+rw_Sat*Cpv)
     .               +h2olv*h2olv*rw_Sat /(Rv  *T__Sat*T__Sat)
     .           *(1.+            rw_Sat * Rv/Rd)/(Cpd+rw_Sat*Cpv))
C +
        ENDIF
C +
        AdiabT(k) = AdiabT(k+1) + dT__dz*(H_Envr(k)-H_Envr(k+1))
C +
      ENDDO
C +
      DO k=1,mz
        TvEnvr(k) = T_Envr(k) 
c #vT.            * (1.+0.61*QvEnvr(k))
        AdiavT(k) = AdiabT(k) 
c #vT.            * (1.+0.61*AdiabQ(k))
      ENDDO
C +
      kLCL=MAX0(kLCL,int_3)
      kLCL=MIN0(kLCL, mz-1)
C +
C +   ----------------------------------------------------------------------
C +
C +---Equilibrium Temperature and Free Convection Levels
C +   ==================================================
C +
C +---Initialisation of search
C +   ------------------------
C +
      zETL=zLCL
      zLFC=zLCL
C +
      DO k=kLCL,1,-1
C +
C +
C +----Equilibrium Temperature level
C +    -----------------------------
C +
        IF (TvEnvr(k+1).le.AdiavT(k+1) .and.
     .      TvEnvr(k)  .ge.AdiavT(k)  ) THEN
C +
          aux=(TvEnvr(k+1)-AdiavT(k+1)+AdiavT(k)-TvEnvr(k))
C +
          IF (abs(aux).gt.1.e-3) THEN      ! ETL detected
            zETL=(H_Envr(k  )*(TvEnvr(k+1)-AdiavT(k+1))
     .           +H_Envr(k+1)*(AdiavT(k  )-TvEnvr(k  )) ) / aux
          ELSE
            zETL=(H_Envr(k)+H_Envr(k+1))*0.5d+0
          ENDIF
C +
        ENDIF
C +
C +
C +----Level of free convection
C +    ------------------------
C +
        IF (TvEnvr(k+1).ge.AdiavT(k+1) .and.
     .      TvEnvr(k)  .le.AdiavT(k)  ) THEN
C +
          aux=(TvEnvr(k+1)-AdiavT(k+1)+AdiavT(k)-TvEnvr(k))
C +
          IF (abs(aux).gt.1.e-3) THEN     ! LFC detected
            zLFC=(H_Envr(k  )*(TvEnvr(k+1)-AdiavT(k+1))
     .           +H_Envr(k+1)*(AdiavT(k  )-TvEnvr(k  )) ) / aux
          ELSE
            zLFC=(H_Envr(k)+H_Envr(k+1))*0.5d+0
          ENDIF
C +
        ENDIF
C +
C +
      ENDDO    !  { k=1,mz-1 }
C +
C +   ----------------------------------------------------------------------
C +
C +---Available / Potential Buoyant Energy
C +   ====================================
C +
      aBuoyE = 0.
      pBuoyE = 0.
C +
      DO k=kLCL+1,2,-1
C +
        aux    =min(zETL,H_Envr(k-1))-max(zLCL,H_Envr(k))
        deltaz1=max(zero,aux)
        deltaz2=0.
        IF (H_Envr(k-1).gt.zLFC) THEN
          aux    =min(zETL,H_Envr(k-1))-max(zLFC,H_Envr(k))
          deltaz2=max(zero,aux)
        ENDIF
C +
C +
C +----Kinetic energy variation (k+1/2)
C +    --------------------------------
C +
        dEnCin=9.81*((AdiavT(k  )-TvEnvr(k  ))/TvEnvr(k  )
     .              +(AdiavT(k-1)-TvEnvr(k-1))/TvEnvr(k-1) ) *0.5d0
        IF (H_Envr(k-1).gt.zLFC)
     .    dEnCin=max(zero,dEnCin)
C +
C +
C +----Available buoyant energy
C +    ------------------------
C +
        aBuoyE=aBuoyE+dEnCin*deltaz1
C +
C +
C +----Potential buoyant energy
C +    ------------------------
C +
        pBuoyE=pBuoyE+dEnCin*deltaz2
C +
      ENDDO
C +
      return
      end
