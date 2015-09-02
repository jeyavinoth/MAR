      block data PHYrad_CEP_DATA
      LOGICAL   RADini
      LOGICAL   RADin2
      common /c_RADini/RADini,RADin2
      data RADini/.FALSE./
      data RADin2/.FALSE./
      end

      subroutine PHYrad_CEP(dST_UA)

! +------------------------------------------------------------------------+
! | MAR PHYSICS                                         Md 11-09-2007  MAR |
! |                                                                        |
! |   SubRoutine PHYrad_CEP  interfaces MAR        with the    new         |
! |              ECMWF Solar/Infrared   Radiative  Transfer Scheme         |
! |                                                                        |
! |   f77 / f90  MAR /ECMWF  Interface                                     |
! |                                                                        |
! |   ECMWF Code Source:  J.-J. Morcrette, 28 nov 2002                     |
! |                                                                        |
! +------------------------------------------------------------------------+


      IMPLICIT NONE


! +--Global Variables (MAR)
! +  ======================

      include 'MARCTR.inc'
      include 'MARphy.inc'

      include 'MARdim.inc'   
      include 'MARgrd.inc'
      include 'MAR_GE.inc'   

      include 'MAR_DY.inc'  

      include 'MAR_HY.inc'  
c #AR include 'MAR_TC.inc'  
      include 'MAR_RA.inc'  

      include 'MAR_SL.inc'    
c #sn include 'MAR_SN.inc'    

      include 'MAR_WK.inc'  
      include 'MAR_IO.inc'  

      REAL     dST_UA                  ! Distance Soleil-Terre [UA]
      REAL     CldFRA(mx,my,mz)        ! Cloud Fraction         [-]
      REAL     RAcldE(mx,my,mz)        ! Cloud Emissivity       [-]
      REAL     RAfnIR(mx,my,mzz)       ! IR  Net Heat Fluxes [W/m2]
      REAL     htngIR(mx,my,mz)        ! IR      Heating      [K/s]
      REAL     htngSO(mx,my,mz)        ! Solar   Heating      [K/s]


! +--Interface  Variables
! +  ====================

      include 'radCEP.inc'


! +--INPUT
! +  -----

      LOGICAL   RADini
      LOGICAL   RADin2
      common /c_RADini/RADini,RADin2

      integer   yymmdd                 ! Date   in the form  yyyyMMdd
      integer   i_hhss                 ! Number of seconds in the day

      real      AlbCEP(klonr)          ! Surface Albedo
      real      pa_CEP(klonr,klevr)    ! Air     Pressure (layer)
      real      pahCEP(klonr,klevr+1)  ! Air     Pressure (layer interface)
      real      fcdCEP(klonr,klevr)    ! Cloud   Fraction (dropplets)
      real      emsCEP(klonr)          ! Surface IR Emissivity
      real      lsmCEP(klonr)          ! Land/Sea Mask: (1.=land  0.=ocean)
      real      cszCEP(klonr)          ! cosine   (solar zenithal Distance)
      real      czeMIN                 ! Minimum accepted for cszCEP
      real      larCEP(klonr)          ! Latitude                  [radian]
      real      lorCEP(klonr)          ! Longitude                 [radian]
      real      ta_CEP(klonr,klevr)    ! Air     Temperature
      real      tasCEP(klonr)          ! Surface Temperature

      real      AerCEP(klonr,nn_aer,klevr)     ! Aerosol Concentration
      real      Ae_MAR(mx,my,nn_aer,mz)        !
      real      O3rCEP(klonr,klevr)            ! O3 Concentration
      real      O3_MAR(mx,my,mz)               !
      common   /MAR_vs_CEP_FIX/ 
     .          Ae_MAR,O3_MAR                  !

      real      cldMAX(mx,my)                  ! Cloud Max  Fraction    [-]
      real      CD_OD1(klonr,klevr)            ! Cloud Optical Depth    [-]
      real      CDtOD1(klonr)                  ! Cloud Optical Depth    [-]
                                               ! (vertically integrated)
      real      Ae_ODa(klonr,klevr)            ! Aeros.Optical Depth    [-]
      real      AetODa(klonr)                  ! Aeros.Optical Depth    [-]
                                               ! (vertically integrated)
      real      qv_CEP(klonr,klevr)            ! Vapor   Concentr.  [kg/kg]
      real      qi_CEP(klonr,klevr)            ! Cryst.  Concentr.  [kg/kg]
      real      qw_CEP(klonr,klevr)            ! Droppl. Concentr.  [kg/kg]
      real      sw_CEP(klonr,klevr)            ! Saturation % water [kg/kg]
      real      qr_CEP(klonr,klevr)            ! Drops   Concentr.  [kg/kg]
      integer   n                              ! 


! +--OUTPUT
! +  ------

      real      FIRn_c(klonr,klevr+1)  ! CLEAR-SKY         LW NET      FLUXES
      real      FIRn_t(klonr,klevr+1)  ! TOTAL             LW NET      FLUXES
      real      FSOn_c(klonr,klevr+1)  ! CLEAR-SKY         SW NET      FLUXES
      real      FSOn_t(klonr,klevr+1)  ! TOTAL             SW NET      FLUXES
      real      FSOs_t(klonr)          ! TOTAL-SKY SURFACE SW DOWNWARD FLUX


      integer   ij0MAX , ij_MAX
      parameter(ij0MAX = mx2    *my2  )
      integer   nbvMAX , nb_MAX
      parameter(nbvMAX = ij0max /klonr)
      integer   klonrb
      parameter(klonrb = ij0max -klonr*nbvMAX)

      integer                  ikl   ,lkl   ,nkl
      integer                  ij    ,nnn   ,nvc
      integer                  k2i(klonr)         ! i index corresp. to kl
      integer                  k2j(klonr)         ! j index corresp. to kl
c #VR integer               ij0ver(mx,my)         ! For Verif. of Vectoriz.
c #VR integer               ij_ver(mx,my)         ! For Verif. of Vectoriz.
c #VR integer           ij2,ijdver(mx,my)         ! For Verif. of Vectoriz.


! +--Surface Albedo
! +  --------------

      real        bsegal   ,albmax   ,albx     ,dalb     ,albu  
      real        czeMAX   ,czrx
      real        siceOK   ,ciceOK   ,zangOK   ,sign_T   ,ColdOK
      real        sign_S   ,snowOK   


! +--OUTPUT
! +  ------

      integer  io
      real     zlevel,pr_atm,qcloud,fcloud
      REAL     heatng                                 ! Total Heating [K/s]


! +--DATA
! +  ====

      data bsegal/2.00e0/      !
      data albmax/0.99e0/      !
      data czeMAX/0.173648178/ ! czeMAX: 80.deg (Segal et al., 1991 JAS)
      data czeMIN/5.00e-6/     !    MIN (Solar Zenithal Distance)


! +--INITIALIZATION
! +  ==============

! #DB open(unit=30,status='unknown',file='PHYrad_CEP.txt')
! #DB rewind    30

      IF (iterun.EQ.0)                                              THEN
            ij_MAX = mx2* my2
        IF      (mod(ij_MAX ,klonr).eq.0)                           THEN
            nb_MAX = ij_MAX /klonr
        ELSE
            nb_MAX = ij_MAX /klonr + 1
        END IF
      END IF

            qcloud = 0.
            fcloud = 0.


! +--Time & Insolation (top of the atmosphere)
! +  =========================================

      yymmdd=min(iyrrGE,2004)*10000+mmarGE*100+jdarGE
      i_hhss=    jhurGE      * 3600+minuGE* 60+jsecGE

! +--Zenith Angle Correction of Snow Albedo
! +  ======================================

      IF  ( mod(itexpe,jtRadi).eq.0)                              THEN ! CTR
        IF(.NOT.VSISVAT)                                          THEN ! CTR

          DO j=jp11,my1
          DO i=ip11,mx1

            siceOK = 1 - min(iabs(isolSL(i,j)-2),iun)
            ciceOK = 1 - min(iabs(isolSL(i,j)-3),iun)

            zangOK =     max(siceOK,ciceOK)

            sign_T =    sign(unun  ,TfSnow-TairSL(i,j))
            ColdOK =     max(zero  ,sign_T)
            zangOK =     max(zangOK,ColdOK)

            sign_S =         zero
c #sn       sign_S =    sign(unun  ,dzSNow(i,j,1)-eps9)
            snowOK =     max(zero  ,sign_S)
            zangOK =     max(zangOK,snowOK)

c #CP       zangOK =         0.0e+0


! +--Snow and/or ice covered surface 
! +  -------------------------------

            albx        =            alb0SL(i,j)
            czrx        = max(czeMAX,czenGE(i,j))
            dalb        = 0.32e0*((bsegal+unun)/(unun+2.e0*bsegal*czrx)
     .                           - unun ) / bsegal
            dalb        = max(dalb,zero)
            albx        =     dalb+alb0SL(i,j)
            albx        = min(albx,albmax)
! +***      Influence of Sun Zenith Angle 
! +         (Segal et al., 1991 JAS 48, p.1025)


! +--Underlying Surface Albedo
! +  -------------------------

            albu        = alb0SL(i,j)


! +--Actual Albedo
! +  -------------

            albeSL(i,j) = zangOK *albx +(1-zangOK) *albu

          END DO
          END DO

        END IF                                                         ! CTR


! +--Effective Radiating Surface Temperature
! +  =======================================

        DO j=jp11,my1
        DO i=ip11,mx1
               WKxy1(i,j) = 0.0
        END DO
        END DO


        DO n=1   ,mw
        DO j=jp11,my1
        DO i=ip11,mx1
               WKxy1(i,j) =  WKxy1(i,j) 
     .      + eps0SL(i,j) * SLsrfl(i,j,n) *tsrfSL(i,j,n)*tsrfSL(i,j,n)
     .                                    *tsrfSL(i,j,n)*tsrfSL(i,j,n)
        END DO
        END DO
        END DO

        DO j=jp11,my1
        DO i=ip11,mx1
               tviRA(i,j) = sqrt( sqrt( WKxy1(i,j) 
     .              +(1.-eps0SL(i,j)) *RAd_ir(i,j)/stefan) )
              cld_SL(i,j) = 0.
              cldMAX(i,j) = 0.
c #VR         ij0ver(i,j) = 0
c #VR         ij_ver(i,j) = 0
c #VR         ijdver(i,j) = 0
        END DO
        END DO


! +--Solar and IR Transfer through the Atmosphere
! +  ============================================

              ij              = 0
c #VR         ij2             = 0


! +--Grid  Point   Dependant Variables --> PHYrad_CEP "Vector"Variables
! +  ------------------------------------------------------------------

        DO    nvc             = 1,nbvMAX

          DO  ikl             = 1,klonr
              ij              =           ij       + 1
               j              = my1 - mod(ij ,my2)
              i               = ip11+    (ij -  1)/my2
              i               =       min(mx1,  i)
              k2i(ikl)        = i
              k2j(ikl)        = j
c #VR         ij0ver(i,j)     = ij0ver(i,j) + 1
c #VR         ijdver(i,j)     = ijdver(i,j) + ij

! +--Geographic Coordinates
! +  ^^^^^^^^^^^^^^^^^^^^^^
              larCEP(ikl)     = GElatr(i,j)
              lorCEP(ikl)     = GElonh(i,j)*hourad
              lorCEP(ikl)     = lorCEP(ikl)
     .       -pi*2.*min(sign(1.,lorCEP(ikl)),0.)

! +--Albedo
! +  ^^^^^^
              AlbCEP(ikl)     = albeSL(i,j)

! +--Surface 
! +  ^^^^^^^
              pahCEP(ikl,mzz) = (pstDY(i,j)*sigmid(mzz)+ptopDY)*1.e3

              emsCEP(ikl)     = eps0SL(i,j)         ! Emissivity
              lsmCEP(ikl) = 1 - maskSL(i,j)         ! Land/sea Mask
              cszCEP(ikl) = max(czenGE(i,j),czeMIN) ! cos(zenith.Dist.)

              tasCEP(ikl)     = tairSL(i,j)

          ENDDO

          DO  lkl=1,mz
          DO  ikl=1,klonr
              i               = k2i(ikl)
              j               = k2j(ikl)

! +--Pressure
! +  ^^^^^^^^
              pahCEP(ikl,lkl) = (pstDY(i,j)*sigmid(lkl)+ptopDY)*1.e3
              pa_CEP(ikl,lkl) = (pstDY(i,j)* sigma(lkl)+ptopDY)*1.e3

! +--Temperature
! +  ^^^^^^^^^^^
              ta_CEP(ikl,lkl) = tairDY(i,j        ,lkl)

! +--Water   Concentration       (qsHY: Von Walden et al. 2003, JAM 42, p.1400)
! +  ^^^^^^^^^^^^^^^^^^^^^       (      Fall  snow 24 mim                     )
!                                (      Blown snow 11 mim, only over ice sheet)
              qv_CEP(ikl,lkl) =   qvDY(i,j        ,lkl)      ! Vapor
              qi_CEP(ikl,lkl) = 0.                           ! Crystals
              qi_CEP(ikl,lkl) =   qiHY(i,j        ,lkl)      !
cXF
     .        + (1.-min(1.,exp((tairDY(i,j,lkl)-273.15)*0.1))) ! XF 258.15
     .        * (                 qsHY(i,j        ,lkl)*0.6  ! Snow(fall)
c #BS.            * max(0,min(3-isolSL(i,j)  ,1))            !
c #BS.                        +   qsHY(i,j        ,lkl)*0.50 ! Snow(fall)
c #BS.            * max(0,min(  isolSL(i,j)-3,1))            !
c #BS.                        +   qsHY(i,j        ,lkl)*0.67 ! Snow(fall+blown)
c #BS.            * max(0,min(  isolSL(i,j)-2,1))            !
c #BS.            * max(0,min(4-isolSL(i,j)  ,1))            !
     .          )
              qw_CEP(ikl,lkl) = 0.                           ! Dropplets
              qw_CEP(ikl,lkl) =       qwHY(i,j,lkl)          !
              sw_CEP(ikl,lkl) = min(qvswDY(i,j,lkl),0.03)    ! Saturation % W
              qr_CEP(ikl,lkl) = 0.                           ! Rain Drops
              qr_CEP(ikl,lkl) =       qrHY(i,j,lkl)          !

! +--Cloud Fraction (liquid water)
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              fcdCEP(ikl,lkl) =     cfraHY(i,j,lkl)

! +--O3      Concentration
! +  ^^^^^^^^^^^^^^^^^^^^^
              O3rCEP(ikl,lkl) =     O3_MAR(i,j,lkl)

          ENDDO
          ENDDO

! +--Aerosol Concentration
! +  ^^^^^^^^^^^^^^^^^^^^^
          DO  nae=1,nn_aer
          DO  lkl=1,mz
          DO  ikl=1,klonr
              AerCEP(ikl,nae,lkl) 
     .                        = Ae_MAR(k2i(ikl),k2j(ikl),nae,lkl)
          END DO
          END DO
          END DO


! +--Radiative Transfert Computation
! +  -------------------------------

! +            **********
          CALL PHYrad2CEP(klonr ,klevr ,nn_aer,yymmdd,i_hhss
     .                   ,dST_UA,AlbCEP,pa_CEP,pahCEP,fcdCEP
     .                   ,emsCEP,lsmCEP,cszCEP,larCEP,lorCEP
     .                   ,AerCEP,O3rCEP,qv_CEP,qi_CEP,qw_CEP
     .                   ,sw_CEP,qr_CEP,ta_CEP,tasCEP
     .                   ,FIRn_c,FIRn_t,FSOn_c,FSOn_t,FSOs_t
     .                   ,CD_OD1,CDtOD1,Ae_ODa,AetODa,iyrrGE
! #DB.                   ,k2i,k2j
     .                   )
! +            **********


! +--Grid  Point   Dependant Variables <-- PHYrad_CEP "Vector"Variables
! +  ------------------------------------------------------------------

          DO  ikl=1,klonr
              i               = k2i(ikl)
              j               = k2j(ikl)
c #VR         ij2             = ij2         + 1
c #VR         ijdver(i,j)     = ijdver(i,j) - ij2
c #VR         ij_ver(i,j)     = ij_ver(i,j) + 1

C +--Surface Cloud/Aerosol Optical Depth
C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAcdtO(i,j)     =    CDtOD1(ikl)
              RAertO(i,j)     =    AetODa(ikl)

C +--Surface Downward Radiative Fluxes
C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAdOLR(i,j)     =-FIRn_t(ikl,1)
              RAdOSR(i,j)     =-FSOn_t(ikl,1)
              RAdsol(i,j)     = FSOs_t(ikl)
              sol_SL(i,j)     = FSOs_t(ikl)*(1.-albeSL(i,j))
              RAd_ir(i,j)     = FIRn_t(ikl,1+klevr)
     .     +  eps0SL(i,j)      *TairSL(i,j)*    TairSL(i,j)
     .                         *TairSL(i,j)*    TairSL(i,j)
     .                         *stefan

C +--Surface IR Net   Radiative Fluxes
C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAfnIR(i,j,mzz) = FIRn_t(ikl,mzz)
          END DO

          DO lkl=1,mz
          DO ikl=1,klonr
              i               = k2i(ikl)
              j               = k2j(ikl)

! +--Atmosph.IR Net   Radiative Fluxes
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAfnIR(i,j,lkl) = FIRn_t(ikl,lkl)

! +--Cloud   Fraction
! +  ^^^^^^^^^^^^^^^^
              cldMAX(i,j) = max(fcdCEP(ikl,lkl),cldMAX(i,j))
              CldFRA(i,j,lkl) = fcdCEP(ikl,lkl)

C +--Atmosph.Cloud/Aerosol Optical Depth
C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAcd_O(i,j,lkl) = CD_OD1(ikl,lkl)
              RAer_O(i,j,lkl) = Ae_ODa(ikl,lkl)

C +--Radiative Heating
C +  ^^^^^^^^^^^^^^^^^
              WKxyz1(i,j,lkl) = -( FIRn_t(ikl,lkl+1)-  FIRn_t(ikl,lkl))
     .                  *  gravit / (cp *1.e3 *pstDY(i,j) *dsigm1(lkl))
              WKxyz2(i,j,lkl) = -( FSOn_t(ikl,lkl+1)-  FSOn_t(ikl,lkl))
     .                  *  gravit / (cp *1.e3 *pstDY(i,j) *dsigm1(lkl))

! +--O3      Concentration
! +  ^^^^^^^^^^^^^^^^^^^^^
              O3_MAR(i,j,lkl) = O3rCEP(ikl,lkl)

          END DO
          END DO

! +--Cloud   Fraction
! +  ^^^^^^^^^^^^^^^^
          DO ikl=1,klonr
              i               = k2i(ikl)
              j               = k2j(ikl)
              cld_SL(i,j)     = cldMAX(i,j)
              clduSL(i,j)     = 0.
              cldmSL(i,j)     = 0.
              clddSL(i,j)     = 0.
          DO lkl=1,mz
           if(pahCEP(ikl,lkl)<44000) 
     .     clduSL(i,j)=max(clduSL(i,j),CldFRA(i,j,lkl)) 
           if(pahCEP(ikl,lkl)>=44000.and.pahCEP(ikl,lkl)<=68000) 
     .     cldmSL(i,j)=max(cldmSL(i,j),CldFRA(i,j,lkl)) 
           if(pahCEP(ikl,lkl)>68000) 
     .     clddSL(i,j)=max(clddSL(i,j),CldFRA(i,j,lkl)) 
          ENDDO
          ENDDO

C +--Radiative Heating
C +  ^^^^^^^^^^^^^^^^^
          DO lkl=1,mz
          DO ikl=1,klonr
              i               = k2i(ikl)
              j               = k2j(ikl)
              pktRAd(i,j,lkl) =  ( WKxyz1(i,j,lkl)  +  WKxyz2(i,j,lkl))
     .                  *   dt  /    pkDY(i,j,lkl)
              htngIR(i,j,lkl) =    WKxyz1(i,j,lkl)  *  86400.
              htngSO(i,j,lkl) =    WKxyz2(i,j,lkl)  *  86400.
          END DO
          END DO

! +--Aerosol Concentration
! +  ^^^^^^^^^^^^^^^^^^^^^
          DO nae=1,nn_aer
          DO lkl=1,mz
          DO ikl=1,klonr
              Ae_MAR(k2i(ikl),k2j(ikl),nae,lkl)
     .                        = AerCEP(ikl,nae,lkl)
          END DO
          END DO
          END DO

        END DO


! +--Grid  Point   Dependant Variables --> PHYrad_CEP "Vector"Variables
! +  ------------------------------------------------------------------

        IF (klonrb.GT.0)                                            THEN

          DO  ikl             = 1,klonr
              ij              =           ij       + 1
               j              = my1 - mod(ij ,my2)
              i               = ip11+    (ij -  1)/my2
              i               =       min(mx1,  i)
              k2i(ikl)        = i
              k2j(ikl)        = j
c #VR         ij0ver(i,j)     = ij0ver(i,j) + 1
c #VR         ijdver(i,j)     = ijdver(i,j) + ij

! +--Geographic Coordinates
! +  ^^^^^^^^^^^^^^^^^^^^^^
              larCEP(ikl)     = GElatr(i,j)
              lorCEP(ikl)     = GElonh(i,j)*hourad
              lorCEP(ikl)     = lorCEP(ikl)
     .       -pi*2.*min(sign(1.,lorCEP(ikl)),0.)

! +--Albedo
! +  ^^^^^^
              AlbCEP(ikl)     = albeSL(i,j)

! +--Surface 
! +  ^^^^^^^
              pahCEP(ikl,mzz) = (pstDY(i,j)*sigmid(mzz)+ptopDY)*1.e3

              emsCEP(ikl)     = eps0SL(i,j)         ! Emissivity
              lsmCEP(ikl) = 1 - maskSL(i,j)         ! Land/sea Mask
              cszCEP(ikl) = max(czenGE(i,j),czeMIN) ! cos(zenith.Dist.)

              tasCEP(ikl)     = tairSL(i,j)

          ENDDO

          DO  lkl=1,mz
          DO  ikl=1,klonr
              i               = k2i(ikl)
              j               = k2j(ikl)

! +--Pressure
! +  ^^^^^^^^
              pahCEP(ikl,lkl) = (pstDY(i,j)*sigmid(lkl)+ptopDY)*1.e3
              pa_CEP(ikl,lkl) = (pstDY(i,j)* sigma(lkl)+ptopDY)*1.e3

! +--Temperature
! +  ^^^^^^^^^^^
              ta_CEP(ikl,lkl) = tairDY(i,j        ,lkl)

! +--Water   Concentration       (qsHY: Von Walden et al. 2003, JAM 42, p.1400)
! +  ^^^^^^^^^^^^^^^^^^^^^       (      Fall  snow 24 mim                     )
!                                (      Blown snow 11 mim, only over ice sheet)
              qv_CEP(ikl,lkl) =   qvDY(i,j        ,lkl)      ! Vapor
              qi_CEP(ikl,lkl) = 0.                           ! Crystals
              qi_CEP(ikl,lkl) =   qiHY(i,j        ,lkl)      !
     .        + (1.-min(1.,exp((tairDY(i,j,lkl)-273.15)*0.1)))
     .        * (                 qsHY(i,j        ,lkl)*0.6  ! Snow(fall)
c #BS.            * max(0,min(3-isolSL(i,j)  ,1))            !
c #BS.                        +   qsHY(i,j        ,lkl)*0.50 ! Snow(fall)
c #BS.            * max(0,min(  isolSL(i,j)-3,1))            !
c #BS.                        +   qsHY(i,j        ,lkl)*0.67 ! Snow(fall+blown)
c #BS.            * max(0,min(  isolSL(i,j)-2,1))            !
c #BS.            * max(0,min(4-isolSL(i,j)  ,1))            !
     .          )
              qw_CEP(ikl,lkl) = 0.                           ! Dropplets
              qw_CEP(ikl,lkl) =   qwHY(i,j        ,lkl)      !
              sw_CEP(ikl,lkl) = qvswDY(i,j        ,lkl)      ! Saturation % W
              qr_CEP(ikl,lkl) = 0.                           ! Rain Drops
              qr_CEP(ikl,lkl) =   qrHY(i,j        ,lkl)      !

! +--Cloud Fraction (liquid water)
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              fcdCEP(ikl,lkl) = cfraHY(i,j        ,lkl)

! +--O3      Concentration
! +  ^^^^^^^^^^^^^^^^^^^^^
              O3rCEP(ikl,lkl) = O3_MAR(i,j        ,lkl)

          ENDDO
          ENDDO

! +--Aerosol Concentration
! +  ^^^^^^^^^^^^^^^^^^^^^
          DO  nae=1,nn_aer
          DO  lkl=1,mz
          DO  ikl=1,klonr
              AerCEP(ikl,nae,lkl) 
     .                        = Ae_MAR(k2i(ikl),k2j(ikl),nae,lkl)
          END DO
          END DO
          END DO


! +--Radiative Transfert Computation
! +  -------------------------------

! +            **********
          CALL PHYrad2CEP(klonr ,klevr ,nn_aer,yymmdd,i_hhss
     .                   ,dST_UA,AlbCEP,pa_CEP,pahCEP,fcdCEP
     .                   ,emsCEP,lsmCEP,cszCEP,larCEP,lorCEP
     .                   ,AerCEP,O3rCEP,qv_CEP,qi_CEP,qw_CEP
     .                   ,sw_CEP,qr_CEP,ta_CEP,tasCEP
     .                   ,FIRn_c,FIRn_t,FSOn_c,FSOn_t,FSOs_t
     .                   ,CD_OD1,CDtOD1,Ae_ODa,AetODa)
! +            **********


! +--Grid  Point   Dependant Variables <-- PHYrad_CEP "Vector"Variables
! +  ------------------------------------------------------------------

          DO  ikl=1,klonrb
              i               = k2i(ikl)
              j               = k2j(ikl)
c #VR         ij2             = ij2         + 1
c #VR         ijdver(i,j)     = ijdver(i,j) - ij2
c #VR         ij_ver(i,j)     = ij_ver(i,j) + 1

C +--Surface Cloud/Aerosol Optical Depth
C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAcdtO(i,j)     =    CDtOD1(ikl)
              RAertO(i,j)     =    AetODa(ikl)

C +--Surface Downward Radiative Fluxes
C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAdOLR(i,j)     =-FIRn_t(ikl,1)
              RAdsol(i,j)     = FSOs_t(ikl)
              sol_SL(i,j)     = FSOs_t(ikl)*(1.-albeSL(i,j))
              RAd_ir(i,j)     = FIRn_t(ikl,1+klevr)
     .     +  eps0SL(i,j)      *TairSL(i,j)*    TairSL(i,j)
     .                         *TairSL(i,j)*    TairSL(i,j)
     .                         *stefan

C +--Surface IR Net   Radiative Fluxes
C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAfnIR(i,j,mzz) = FIRn_t(ikl,mzz)
          END DO

          DO lkl=1,mz
          DO ikl=1,klonrb
              i               = k2i(ikl)
              j               = k2j(ikl)

! +--Atmosph.IR Net   Radiative Fluxes
! +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAfnIR(i,j,lkl) = FIRn_t(ikl,lkl)

! +--Cloud   Fraction
! +  ^^^^^^^^^^^^^^^^
              cldMAX(i,j) = max(fcdCEP(ikl,lkl),cldMAX(i,j))
              CldFRA(i,j,lkl) = fcdCEP(ikl,lkl)

C +--Atmosph.Cloud/Aerosol Optical Depth
C +  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              RAcd_O(i,j,lkl) = CD_OD1(ikl,lkl)
              RAer_O(i,j,lkl) = Ae_ODa(ikl,lkl)

C +--Radiative Heating
C +  ^^^^^^^^^^^^^^^^^
              WKxyz1(i,j,lkl) = -( FIRn_t(ikl,lkl+1)-  FIRn_t(ikl,lkl))
     .                  *  gravit / (cp *1.e3 *pstDY(i,j) *dsigm1(lkl))
              WKxyz2(i,j,lkl) = -( FSOn_t(ikl,lkl+1)-  FSOn_t(ikl,lkl))
     .                  *  gravit / (cp *1.e3 *pstDY(i,j) *dsigm1(lkl))

! +--O3      Concentration
! +  ^^^^^^^^^^^^^^^^^^^^^
              O3_MAR(i,j,lkl) = O3rCEP(ikl,lkl)

          END DO
          END DO

! +--Cloud   Fraction
! +  ^^^^^^^^^^^^^^^^

          DO  ikl=1,klonrb
              i               = k2i(ikl)
              j               = k2j(ikl)
              cld_SL(i,j)     = cldMAX(i,j)
              clduSL(i,j)     = 0.
              cldmSL(i,j)     = 0.
              clddSL(i,j)     = 0.
          DO lkl=1,mz
           if(pahCEP(ikl,lkl)<44000) 
     .     clduSL(i,j)=max(clduSL(i,j),CldFRA(i,j,lkl)) 
           if(pahCEP(ikl,lkl)>=44000.and.pahCEP(ikl,lkl)<=68000) 
     .     cldmSL(i,j)=max(cldmSL(i,j),CldFRA(i,j,lkl)) 
           if(pahCEP(ikl,lkl)>68000) 
     .     clddSL(i,j)=max(clddSL(i,j),CldFRA(i,j,lkl)) 
          ENDDO
          ENDDO

C +--Radiative Heating
C +  ^^^^^^^^^^^^^^^^^
          DO lkl=1,mz
          DO ikl=1,klonrb
              i               = k2i(ikl)
              j               = k2j(ikl)
              pktRAd(i,j,lkl) =  ( WKxyz1(i,j,lkl)  +  WKxyz2(i,j,lkl))
     .                  *   dt  /    pkDY(i,j,lkl)
              htngIR(i,j,lkl) =    WKxyz1(i,j,lkl)  *  86400.
              htngSO(i,j,lkl) =    WKxyz2(i,j,lkl)  *  86400.
          END DO
          END DO

! +--Aerosol Concentration
! +  ^^^^^^^^^^^^^^^^^^^^^
          DO nae=1,nn_aer
          DO lkl=1,mz
          DO ikl=1,klonrb
              Ae_MAR(k2i(ikl),k2j(ikl),nae,lkl)
     .                        = AerCEP(ikl,nae,lkl)
          END DO
          END DO
          END DO

        END IF


! +--Lateral Boundary Conditions for Radiative Variables
! +  ===================================================

        DO k=1,mz
          DO j=1,my
            pktRAd( 1,j,k) = pktRAd(ip11,j,k)
            pktRAd(mx,j,k) = pktRAd( mx1,j,k)
          END DO
          DO i=1,mx
            pktRAd(i, 1,k) = pktRAd(i,jp11,k)
            pktRAd(i,my,k) = pktRAd(i, my1,k)
          END DO
        END DO

        DO k=1,mz
          DO j=1,my
          DO i=1,mx
            WKxyz1(i,j,k) = 0.
            WKxyz2(i,j,k) = 0.
          END DO
          END DO
        END DO


! +--OUTPUT, write
! +  =============

        IF (  jmmMAR.eq.0.and.    jssMAR   .eq.0       .and.
     .      ((IO_loc.ge.2.and.mod(jhurGE,3).eq.0).or.
     .       (IO_loc.ge.3                            ))     )       THEN

          DO  io = 1,min(5,mx)
              i  = igrdIO(io)
              j  = jgrdIO(io)

C +         ***********
            call TIMcor
C +         ***********


            write(4,401)jdplus,mmplus,jhlrGE(i ,j ),minuGE,
     .                                itizGE(i ,j ), i ,j
  401       format(//,' +-- Radiative  Heat Fluxes --+',
     .       i4,'/',i2,i4,':',i2,' h.LT (',i3,')',
     .               '   (i,j) = (',i3,',',i3,')',
     .        //,' | Altitud | Pressur.| Temper.|  Ozone  | W.Vapor|',
     .      ' Clouds | Clouds |  Opt.D | Aer.OD. | So Warm.|',
     .       ' Emiss.| IR Cool.| IR NetF.|',
     .         /,' |   [km]  |   [hPa] |   [K]  | [cmSTP] | [g/kg] |',
     .      ' [g/kg] |   [%]  |   [-]  |   [-]   | [K/day] |',
     .       '  [-]  | [K/day] |  [W/m2] |',
     .         /,' +---------+---------+--------+---------+--------+',
     .              '--------+--------+--------+---------+---------+',
     .       '-------+---------+---------+')

            write(4,404)                           ptopDY*1.e+1 
     .           ,                                 RAdOLR(i,j)
  404       format(
     .                     ' |  SOMMET |',f8.2,' |', 7x ,' |', 8x ,' |'
     .        , 3( 7x ,' |'),    7x ,' |', 8x ,' |', 8x ,' |'
     .        ,                            6x ,' |', 8x ,' |',f8.1,' |'
     .        ,/,' +---------+---------+--------+---------+--------+',
     .              '--------+--------+--------+---------+---------+',
     .       '-------+---------+---------+')
            DO k =1,mz
               zlevel = 1.e-3* grvinv    *gplvDY(i ,j ,k)
               pr_atm = 1.e+1*(pstDY(i,j)* sigma(k)+ptopDY)
               qcloud = 1.e+3*( qwHY(i,j,k)+qiHY(i ,j ,k))
               fcloud = 1.e+2            *CldFRA(i ,j ,k)
               write(4,402)          zlevel       ,pr_atm       
     .           ,tairDY(i,j,k)     ,O3_MAR(i,j,k)
     .           ,  qvDY(i,j,k)*1.e3,qcloud       ,fcloud       
     .           ,RAcd_O(i,j,k)     ,RAer_O(i,j,k),htngSO(i,j,k)
     .           ,RAcldE(i,j,k)                   ,htngIR(i,j,k)     
     .                                            ,RAfnIR(i,j,k)
  402          format( ' | ',  f7.3,' |' ,f8.2,' |',f7.2,' |',e8.2,' |',
     .          2(f7.3,' |'),2(f7.2,' |'),e8.2,' |',f8.4,' |'
     .          ,                         f6.3,' |',f8.4,' |',f8.1,' |')
            END DO

               pr_atm = 1.e+1*(pstDY(i,j)         +ptopDY)

            write(4,403)                           pr_atm       
     .           ,TairSL(i,j)       
     .           ,RAcdtO(i,j)       ,RAertO(i,j)  ,RAdsol(i,j)
     .                                            ,RAfnIR(i,j,mzz)
     .                              ,albeSL(i,j)  ,sol_SL(i,j)
     .                              ,eps0SL(i,j)  ,RAd_ir(i,j)
     .                              ,stefan       * WKxy1(i,j)
  403       format(
     .           ' +---------+---------+--------+---------+--------+',
     .              '--------+--------+--------+---------+---------+',
     .       '-------+---------+---------+',
     .                   /,' | AIR-SOL |',f8.2,' |',f7.2,' |', 8x ,' |'
     .        , 3( 7x ,' |'),   f7.2,' |',e8.2,' |',f8.1,' |'
     .        ,                            6x ,' |', 8x ,' |',f8.1,' |'
     .        ,          /,' |     SOL |', 8x ,' |', 7x ,' |', 8x ,' |'
     .        , 3( 7x ,' |'),    7x ,' |',f8.2,' |',f8.1,' |'
     .        ,                           f6.3,' |',f8.1,' |',f8.1,' |')
          END DO
        END IF


! C +--OUTPUT for Verification
! C +  -----------------------

c #VR   write(6,6000)
 6000   format(/,'Verification of Vectorization: Before CALL')
c #VR   DO j=my,1,-1
c #VR   write(6,6010) (ij0ver(i,j),i=1,mx)
 6010   format(132i1)
c #VR   ENDDO

c #VR   write(6,6001)
 6001   format(/,'Verification of Vectorization: After  CALL')
c #VR   DO j=my,1,-1
c #VR   write(6,6010) (ij_ver(i,j),i=1,mx)
c #VR   ENDDO

c #VR   DO j=1,my
c #VR   DO i=1,mx
c #VR     IF (ijdver(i,j).NE.0 .AND. ij_ver(i,j).NE.1) 
c #VR.      write(6,6002) i,j,ijdver(i,j)
 6002       format(' Vectorization ERROR on',2i4,'   (',i6,')')
c #VR   ENDDO
c #VR   ENDDO


      END IF 


! +--TURN OFF "PHYrad_CEP" Initialization Procedure
! +  ==============================================

      RADin2 = .TRUE.

! #DB close (unit=30)

      return
      end
