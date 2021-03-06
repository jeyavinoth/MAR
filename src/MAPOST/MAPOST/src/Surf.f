C +----------------------------------------------------------------------+
C | MAR post-processing  XF                                 02/02/2002   | 
C |                                                              v.3.0   |
C | Surf                                                                 |
C |                                                                      |
C +----------------------------------------------------------------------+
C
      SUBROUTINE Surf(idRLS, idMAR, itM, itR, istat,
     $                InReg, MARlon, MARlat, jhuCUR, mmaCUR,
     $                MAPtyp, iOUTnc,idt1)

      IMPLICIT NONE
 
C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
      
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'LSMphy.inc'

C +...* Internal values of Surf (retain for next call):
      include 'Surf.inc'
C     test (for new developments):
      REAL tmp3D(mx,my,mz)
      REAL tT2MAR(mx,my), tT2RLS(mx,my)
      SAVE tT2MAR       , tT2RLS

C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM, idRLS, itR, idOBS
      INTEGER InReg(mx,my,nreg)
      REAL MARlon(mx,my), MARlat(mx,my)
      INTEGER istat,jhuCUR,mmaCUR,MAPtyp,iOUTnc,idt1

C +---OUTPUT 
C +   ~~~~~~
C     (none)

C +---LOCAL
C +   ~~~~~
      CHARACTER *10 vUnits, tmp_units
      CHARACTER *100 nfOBS, OBStit
C +-  -Gridded observations size
      INTEGER OBni, OBnj
      PARAMETER(OBni=720,OBnj=360)
      LOGICAL readgOBS
C +-  -Coordonnees:
      REAL LSlon(LSni), LSlat(LSnj), LSlev(LSnk)
      REAL MARx(mx), MARy(my), empty1(1), sigma(mz)
      REAL ptop, CSTpMA(mzz), SIGpMA(mzz)
      REAL CSTp (LSnk1), SIGp (LSnk1), wkzRLS(LSnk)
      REAL OBlon(OBni),OBlat(OBnj)
C +-  -Valeurs lues
      REAL shMAR(mx,my), shRLS(LSni,LSnj), shIRLS(mx,my)
      REAL PlspRLS(LSni,LSnj), Pcp_RLS(LSni,LSnj)
      REAL Psf_RLS(LSni,LSnj), LSM_LS (LSni,LSnj)
      REAL T2_MAR(mx,my), T2_RLS (LSni,LSnj)
      REAL ALB_MAR(mx,my), ALB_RLS (LSni,LSnj) 
      REAL tmpOBS(OBni,OBnj)
C +-  -Intermediaires de calcul
      REAL Rrr_MAR(nreg), Rcp_MAR(nreg),Rsf_MAR(nreg),Rtot_MAR(nreg)
      REAL RrrIRLS(nreg), RcpIRLS(nreg),RsfIRLS(nreg),RtotIRLS(nreg)
      REAL rrIRLS(mx,my), cpIRLS(mx,my),sfIRLS(mx,my),totIRLS(mx,my)
      REAL rrIOBS(mx,my) 
      REAL T2IRLS(mx,my)
      REAL ALBIRLS(mx,my) 
      REAL rT2IRLS(nreg), rT2_MAR(nreg), rTmiMAR(nreg)
      REAL rALBIRLS(nreg), rALB_MAR(nreg) 
      REAL rTmaMAR(nreg)
      REAL TOTmaxMAR
      REAL tmp_in(OBni,OBnj), samOx(mx,my,5), samOy(mx,my,5)

C +-  -Indices, ...
      INTEGER kl, ii, jj, nelem, ireg
      INTEGER ildt, ilstat


C +---Read general MAR data     
C +   ---------------------

      CALL UNread
     &   (idMAR, 'sh', 0,0, 1,1,mx,my,1,
     &    MARx, MARy, empty1, vUnits, shMAR)

C +---Read general RLS data       
C +   ---------------------
      CALL UNread
     &   (idRLS, 'SH', 0,0, 1,1,LSni,LSnj,1,
     &    LSlon, LSlat, empty1, vUnits, shRLS)

      CALL INThor (-1, LSlon , LSlat , shRLS,
     &                 MARlon, MARlat, shIRLS)


C +===Precipitation.
C +   ==============

C +   IF variable is not available, warn and set=0.0:
C +   (NOVAR_WARN level 1 = replace and warn, don't stop)
      CALL UNparam('NOVAR_REPLACE',0.0)
      CALL UNparam('NOVAR_WARNING',1.0)

C +- -RLS/ERA
C +   - - - - - - - - - -
C +   ERA (ECMWF) => from 6 hour to current (??) LSP- 
C         PxxxRLS variables = Period (usually 6 hour) data
C         (! usually more MAR data then RLS=> ok. 
C          more RLS would not be allowed: missing period)

cXF  
      CALL UNsread
     &   (idRLS,'LSP'   ,itR,0, 
     &    1,1,LSni,LSnj,1, vUnits,PlspRLS) !Large Scale

c     DO jj=1,LSnj
c     DO ii=1,LSni
c       PlspRLS(ii,jj) =  0.0  
c     ENDDO
c     ENDDO

      CALL UNsread
     &   (idRLS,'CP'     ,itR,0,
     &    1,1,LSni,LSnj,1, vUnits,Pcp_RLS) !Convective

      CALL UNsread
     &   (idRLS,'SF'     ,itR,0,
     &    1,1,LSni,LSnj,1, vUnits,Psf_RLS) !Snow Fall

C +   *Compute cumulated precipitation:

      CALL Cumul (PlspRLS, LSni, LSnj, istat, lspRLS)
      CALL Cumul (Pcp_RLS, LSni, LSnj, istat, cp_RLS)
      CALL Cumul (Psf_RLS, LSni, LSnj, istat, sf_RLS)

      DO jj=1,LSnj
      DO ii=1,LSni
         rr_RLS(ii,jj) =  lspRLS(ii,jj) + cp_RLS(ii,jj)
     .                 -  sf_RLS(ii,jj)
        tot_RLS(ii,jj) =  lspRLS(ii,jj) + cp_RLS(ii,jj)      
      ENDDO
      ENDDO

C +- -MAR
C +   - - - - - - - - - -
C +   MAR=> Cumulated -> use end of month

      CALL UNsread
     &   (idMAR,'rainHY' ,itM,0,
     &    1,1,mx  ,my  ,1, vUnits,rr_MAR) !Total Rain

      CALL UNsread
     &   (idMAR,'rainCA' ,itM,0,
     &    1,1,mx  ,my  ,1, vUnits,cp_MAR) !Convective
   
      CALL UNsread
     &   (idMAR,'snowHY' ,itM,0,
     &    1,1,mx  ,my  ,1, vUnits,sf_MAR) !Snow Fall
    
       DO jj=1,my
       DO ii=1,mx
         tot_MAR(ii,jj) =  rr_MAR(ii,jj) + sf_MAR(ii,jj)
       ENDDO
       ENDDO

C +. -Memorize initial value for "cumulated rain" in MAR
C     - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (istat.EQ.0) THEN
       DO jj=1,my
       DO ii=1,mx
         rr_iMAR(ii,jj) =  rr_MAR(ii,jj)
       ENDDO
       ENDDO
      ENDIF

      IF (istat.EQ.0) THEN
       DO jj=1,my
       DO ii=1,mx
         cp_iMAR(ii,jj) =  cp_MAR(ii,jj)
       ENDDO
       ENDDO
      ENDIF

      IF (istat.EQ.0) THEN
       DO jj=1,my
       DO ii=1,mx
         sf_iMAR(ii,jj) =  sf_MAR(ii,jj)
       ENDDO
       ENDDO
      ENDIF

      IF (istat.EQ.0) THEN
       TOTmaxMAR=0.0
       DO jj=1,my
       DO ii=1,mx
         tot_iMAR(ii,jj) =  tot_MAR(ii,jj)
         TOTmaxMAR = max(TOTmaxMAR, tot_MAR(ii,jj))
       ENDDO
       ENDDO
       IF (TOTmaxMAR.GE.0.001) THEN
       write(*,*) 'WARNING (Surf.f):'
       write(*,*) '  Total MAR precip is not 0 at run start'
       write(*,*) '  (this may be ok). Max. init. is', TOTmaxMAR
       ENDIF
      ENDIF

C +. -Write instant. values of RegM[RR] (time)
C     - - - - - - - - - - - - - - - - - - - - -
      CALL HORmean (InReg, rr_MAR, Rrr_MAR)
      CALL UNwrite (iOUTnc, 'r_RR_M', idt1,nreg,1,1, Rrr_MAR)

      CALL HORmean (InReg, sf_MAR, Rsf_MAR)
      CALL UNwrite (iOUTnc, 'r_SF_M', idt1,nreg,1,1, Rsf_MAR)

      CALL HORmean (InReg, tot_MAR, Rtot_MAR)
      CALL UNwrite (iOUTnc, 'r_TOT_M', idt1,nreg,1,1, Rtot_MAR)

      CALL HORmean (InReg, cp_MAR, Rcp_MAR)
      CALL UNwrite (iOUTnc, 'r_CP_M', idt1,nreg,1,1, Rcp_MAR)

      CALL INThor (-1, LSlon , LSlat , rr_RLS,
     &                 MARlon, MARlat, rrIRLS)
      CALL HORmean (InReg, rrIRLS, RrrIRLS)
      CALL UNwrite (iOUTnc, 'r_RR_R', idt1,nreg,1,1, RrrIRLS)

      CALL INThor (-1, LSlon , LSlat , sf_RLS,
     &                 MARlon, MARlat, sfIRLS)
      CALL HORmean (InReg, sfIRLS, RsfIRLS)
      CALL UNwrite (iOUTnc, 'r_SF_R', idt1,nreg,1,1, RsfIRLS)

      CALL INThor (-1, LSlon , LSlat , tot_RLS,
     &                 MARlon, MARlat, totIRLS)
      CALL HORmean (InReg, totIRLS, RtotIRLS)
      CALL UNwrite (iOUTnc, 'r_TOT_R', idt1,nreg,1,1, RtotIRLS)

      CALL INThor (-1, LSlon , LSlat , cp_RLS,
     &                 MARlon, MARlat, cpIRLS)
      CALL HORmean (InReg, cpIRLS, RcpIRLS)
      CALL UNwrite (iOUTnc, 'r_CP_R', idt1,nreg,1,1, RcpIRLS)


C +---Last time-step only: get cumulated value + write.
C +   -------------------------------------------------
      IF (istat.EQ.2) THEN

        DO jj=1,my
        DO ii=1,mx
          rr_MAR(ii,jj) =  rr_MAR(ii,jj) - rr_iMAR(ii,jj)
        ENDDO
        ENDDO

        DO jj=1,my
        DO ii=1,mx
          sf_MAR(ii,jj) =  sf_MAR(ii,jj) - sf_iMAR(ii,jj)
        ENDDO
        ENDDO
        
        DO jj=1,my
        DO ii=1,mx       
          cp_MAR(ii,jj) =  cp_MAR(ii,jj) - cp_iMAR(ii,jj)       
        ENDDO
        ENDDO

        DO jj=1,my
        DO ii=1,mx
          tot_MAR(ii,jj) =  tot_MAR(ii,jj) - tot_iMAR(ii,jj)
        ENDDO
        ENDDO

        CALL UNwrite (iOUTnc, 'C_RR_M', 1, mx, my, 1, rr_MAR)
        CALL UNwrite (iOUTnc, 'C_RR_R', 1, mx, my, 1, rrIRLS)                
        CALL UNwrite (iOUTnc, 'C_CP_M', 1, mx, my, 1, CP_MAR)        
        CALL UNwrite (iOUTnc, 'C_CP_R', 1, mx, my, 1, cpIRLS)     
        CALL UNwrite (iOUTnc, 'C_SF_M', 1, mx, my, 1, sf_MAR)        
        CALL UNwrite (iOUTnc, 'C_SF_R', 1, mx, my, 1, sfIRLS)             
        CALL UNwrite (iOUTnc, 'C_TOT_M', 1, mx, my, 1, tot_MAR)        
        CALL UNwrite (iOUTnc, 'C_TOT_R', 1, mx, my, 1, totIRLS)   
        

        CALL TransGrid  (MAPtyp,
     &      MARlon, MARlat, MARx, MARy, rr_MAR,
     &      LSlon, LSlat, rr_TMAR)

        CALL UNwrite (iOUTnc, 'C_RR_LR', 0,LSni,LSnj,1,rr_RLS)
        CALL UNwrite (iOUTnc, 'C_RR_LM', 0,LSni,LSnj,1,rr_TMAR)    

        CALL TransGrid  (MAPtyp,
     &      MARlon, MARlat, MARx, MARy, cp_MAR,
     &      LSlon, LSlat, cp_TMAR)

        CALL UNwrite (iOUTnc, 'C_CP_LR', 0,LSni,LSnj,1,cp_RLS)
        CALL UNwrite (iOUTnc, 'C_CP_LM', 0,LSni,LSnj,1,cp_TMAR)  

        CALL TransGrid  (MAPtyp,
     &      MARlon, MARlat, MARx, MARy, sf_MAR,
     &      LSlon, LSlat, sf_TMAR)

        CALL UNwrite (iOUTnc, 'C_SF_LR', 0,LSni,LSnj,1,sf_RLS)
        CALL UNwrite (iOUTnc, 'C_SF_LM', 0,LSni,LSnj,1,sf_TMAR)    

        CALL TransGrid  (MAPtyp,
     &      MARlon, MARlat, MARx, MARy, tot_MAR,
     &      LSlon, LSlat, tot_TMAR)

        CALL UNwrite (iOUTnc, 'C_TOT_LR', 0,LSni,LSnj,1,tot_RLS)
        CALL UNwrite (iOUTnc, 'C_TOT_LM', 0,LSni,LSnj,1,tot_TMAR) 

      ENDIF

C +===Surface Albedo
C +   ==============

cXF
      CALL UNsread (idRLS,'FAL',itR,0,
     &    1,1,LSni,LSnj,1, vUnits,ALB_RLS) 

      CALL UNsread (idMAR,'albeSL',itM,0, 
     &    1,1,mx  ,my  ,1, vUnits,ALB_MAR) 

C +. -Interpolate RLS values to MAR grid.
C     - - - - - - - - - - - - - - - - - - - - -
      CALL INThor (-1, LSlon , LSlat , ALB_RLS,
     &                 MARlon, MARlat, ALBIRLS)


C +- -Compute & write tM[ALB] & tSD[ALB]
C     - - - - - - - - - - - - - - - - - - -

      CALL HDynSTA2D(iOUTnc,istat,idt1,ALB_MAR,
     &               'tM_al_M','tSD_al_M',tM_al_M,tSD_al_M)
      CALL HDynSTA2D(iOUTnc,istat,idt1,ALBIRLS,
     &               'tM_al_R','tSD_al_R',tM_al_R,tSD_al_R)


C +. -Write instant. values of RegM[alb] (time)
C     - - - - - - - - - - - - - - - - - - - - -
      CALL HORmean (InReg, ALBIRLS, rALBIRLS)
      CALL UNwrite (iOUTnc, 'rM_al_R'  , idt1,nreg,1,1, rALBIRLS)


      CALL HORmean (InReg, ALB_MAR, rALB_MAR)
      CALL UNwrite (iOUTnc, 'rM_al_M'  , idt1,nreg,1,1, rALB_MAR)

C +===2m air temperature.
C +   ===================

C +- -RLS/ERA
C +   - - - - 
      CALL UNsread (idRLS,'2T'   ,itR,0,
     &    1,1,LSni,LSnj,1, vUnits,T2_RLS) !2m Temp       


C +- -MAR
C +   - - 
c     CALL UNsread (idMAR,'Ta2mSL',itM,0,
c    &    1,1,mx  ,my  ,1, vUnits,T2_MAR) !2m Temp
cXF
c     write(*,*) 'WARNING = temporary T2'
c     CALL UNsread (idMAR,'tsrfSL',itM,1,
c    &    1,1,mx  ,my ,3, vUnits,T2_MAR) !2m Temp

      CALL UNsread (idMAR,'tairDY',itM,1,
     &    1,1,mx  ,my , mz, vUnits,tmp3D) !2m Temp

      do ii=1,mx ; do jj=1,my
       T2_MAR(ii,jj) = tmp3D(ii,jj,mz)
      end do ; end do 
      ! On GRD, 1st level = 2m
    

C +   -* T_min and _max are valid only at 6h in current MAR
C +      (may change with versions: complicated algorithm.)
C 
C       De toute facon, Tmin/max ici = pas au point !
C       (notamment : ecrit pour journee suivante ds fichier...)
 
      IF (jhuCUR.EQ.6) THEN

        CALL UNsread (idMAR,'TminSL',itM,0,
     &    1,1,mx  ,my  ,1, vUnits,TmiMAR) !2m min T

        CALL UNsread (idMAR,'TmaxSL',itM,0,
     &    1,1,mx  ,my  ,1, vUnits,TmaMAR) !2m max T
   
      ENDIF


C +. -Interpolate RLS values to MAR grid.      
C     - - - - - - - - - - - - - - - - - - - - -
      CALL INThor (-1, LSlon , LSlat , T2_RLS,
     &                 MARlon, MARlat, T2IRLS)
C +   Optional: attempt to correct the LS value for the
C +   topo difference => go to MAR sh using -6.5 K/km grad:
      DO jj=1,my
      DO ii=1,mx
         T2IRLS(ii,jj) = T2IRLS(ii,jj)
     &                 + gamTz*(shMAR(ii,jj)-shIRLS(ii,jj)) 
      ENDDO
      ENDDO
      IF (istat.EQ.2) THEN
        write(*,*) 'INFO: T air 2m: rough cor. for topo dif. ON'
      ENDIF


C +. -Write instant. values of RegM[T2m] (time)
C     - - - - - - - - - - - - - - - - - - - - -
      CALL HORmean (InReg, T2IRLS, rT2IRLS)
      CALL UNwrite (iOUTnc, 'rT2_R'  , idt1,nreg,1,1, rT2IRLS)


      CALL HORmean (InReg, T2_MAR, rT2_MAR)
      CALL UNwrite (iOUTnc, 'rT2_M'  , idt1,nreg,1,1, rT2_MAR)

      CALL HORmean (InReg, TmiMAR, rTmiMAR)
      CALL UNwrite (iOUTnc, 'rT2minM', idt1,nreg,1,1, rTmiMAR)

      CALL HORmean (InReg, TmaMAR, rTmaMAR)
      CALL UNwrite (iOUTnc, 'rT2maxM', idt1,nreg,1,1, rTmaMAR)
 

C +. -Compute & write mean temperature.            
C     - - - - - - - - - - - - - - - - -
C     (This routine was created for HDyn, but is perfect 
C      also here...)
C
      ildt  = idt1-1
      ilstat= istat
      IF (idt1.LE.2) ilstat = 0
C     ^The first data does not exist, so reset stats to 0
C      after ther first step (cause: MAR calc of T2)

      CALL HDynSTVAL (iOUTnc,ilstat,ildt,rT2IRLS,
     &               'tMrT2_R',tMT2RLS)

      CALL HDynSTVAL (iOUTnc,ilstat,ildt,rT2_MAR,
     &               'tMrT2_M',tMT2MAR)

      CALL HDynSTVAL (iOUTnc,ilstat,ildt,rTmiMAR,
     &               'tMrTmiM',tMTmiMAR)

      CALL HDynSTVAL (iOUTnc,ilstat,ildt,rTmaMAR,
     &               'tMrTmaM',tMTmaMAR)

C     ? not ready ?

      CALL Mean2D  (iOUTnc,istat,idt1,T2_MAR,
     &                 'tM_T2_M',tT2MAR)

      CALL Mean2D  (iOUTnc,istat,idt1,T2_RLS,
     &                 'tM_T2_R',tT2RLS)




C +===Copy LS Land Sea Mask for graphics
C +   ==================================
      IF (istat.EQ.2) THEN

         CALL UNsread (idRLS,'LSM',0,0,
     &       1,1,LSni,LSnj,1, vUnits,LSM_LS)

         CALL UNwrite (iOUTnc,'LSM_L', 0,LSni,LSnj,1,LSM_LS)

      ENDIF
      
C +   Go back to standard treatment of missing variables:
C +   (warn level 2 = standard = stop all)
      CALL UNparam('NOVAR_WARNING',2.0)


C +- -Observations (CRU) 
C +   - - - - - - - - - -
C     (Cumulated / one month)
      readgOBS=.false. 
C      WARNING: under development 
      IF (readgOBS.AND.(istat.EQ.0)) THEN

        nfOBS = 'CRU05_1986.nc'
        WRITE(*,*) 'Writing CRU surface clim for month:', mmaCUR
        CALL UNropen (nfOBS, idOBS, OBStit)
        CALL UNread
     &   (idOBS, 'pre', mmaCUR,0, 1,1,OBni,OBnj,1,
     &    OBlon, OBlat, empty1, vUnits, tmpOBS)

        CALL UNclose(idOBS)

        CALL INTmean (tmp_in, samOx, samOy,
     &       OBni, OBnj, OBlon, OBlat, tmpOBS,
     &       mx, my, MARlon, MARlat, rrIOBS)

        CALL UNwrite (iOUTnc, 'C_TOT_O', 1, mx, my, 1, rrIOBS)

       ENDIF

      RETURN
      END
