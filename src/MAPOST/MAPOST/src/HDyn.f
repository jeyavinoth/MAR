C +----------------------------------------------------------------------+
C | MAR post-processing                                     02/02/2002   | 
C |                                                              3.0.0   |
C | HDyn                                                                 |
C |  Horizontal, mostly "dynamics oriented" maps & stat. values.         |
C +----------------------------------------------------------------------+
C
      SUBROUTINE HDyn(idRLS, idMAR, itM, itR, istat,
     $                InReg, MARlon, MARlat, jhuCUR,
     $                ioutNC,idt1,slpMAR)

      IMPLICIT NONE
 
C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
            
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'globals.inc'
      include 'LSMphy.inc'

C +...* Z500 filter size:
      include 'HDynSDBP.inc'
      
C +...* Internal values of HDyn (retain for next call):
      include 'HDyn.inc'

C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM, idRLS, itR
      INTEGER InReg(mx,my,nreg)
      REAL MARlon(mx,my), MARlat(mx,my)
      INTEGER istat,jhuCUR,ioutNC,idt1

C +---OUTPUT 
C +   ~~~~~~
C     (none)

C +---LOCAL
C +   ~~~~~
      CHARACTER *10 var_units, tmp_units
C +-  -Coordonnees:
      REAL LSlon(LSni), LSlat(LSnj), LSlev(LSnk)
      REAL LSlonU(LSni), LSlatU(LSnj)
      REAL LSlonV(LSni), LSlatV(LSnj)

      REAL MARx(mx), MARy(my), empty1(1), sigma(mz)
      REAL ptop, CSTpMA(mzz), SIGpMA(mzz)
      REAL CSTp (LSnk1), SIGp (LSnk1), wkzRLS(LSnk)
C +-  -Valeurs lues
      REAL pstMAR(mx,my), shMAR(mx,my)
      REAL tairDY(mx,my,mz), qvDY (mx,my,mz), zzDY(mx,my,mz)
      REAL uairDY(mx,my,mz), vairDY(mx,my,mz)
      REAL spRLS(LSni,LSnj), shRLS(LSni,LSnj)
      REAL T_RLS(LSni,LSnj,LSnk), Q_RLS(LSni,LSnj,LSnk)
      REAL U_RLS(LSni,LSnj,LSnk), V_RLS(LSni,LSnj,LSnk)
C +-  -Valeurs Calculees
      REAL plvMAR(mx,my,mz)      , zz_MAR(mx,my,mz) 
      REAL lplMAR(mx,my,mz)
      REAL plvRLS(LSni,LSnj,LSnk), zz_RLS(LSni,LSnj,LSnk)
      REAL lplRLS(LSni,LSnj,LSnk)
      REAL zplMAR(mx,my), zplRLS(LSni,LSnj), zplIRLS(mx,my)
      REAL zplBIA(mx,my), RzplBIA(nreg)
      REAL zplRME(mx,my), RzplRME(nreg)
      REAL slpMAR(mx,my), slpRLS(LSni,LSnj), slpIRLS(mx,my)
      REAL slpBIA(mx,my), slpRME(mx,my)
      REAL RslpBIA(nreg), RslpRME(nreg)
      REAL RslpMAR(nreg), RslpRLS(nreg)
      REAL TplMAR(mx,my), TplRLS(LSni,LSnj), TplIRLS(mx,my)
      REAL RTplMAR(nreg), RTplRLS(nreg)
      REAL QplMAR(mx,my), QplRLS(LSni,LSnj), QplIRLS(mx,my)
      REAL RQplMAR(nreg), RQplRLS(nreg)
      REAL UV_MAR(mx,my,mz), UV_RLS(LSni,LSnj,LSnk)
      REAL UVpMAR(mx,my), UVpRLS(LSni,LSnj), UVpIRLS(mx,my)
      REAL UplMAR (mx,my), UplRLS (LSni,LSnj), UplIRLS (mx,my)
      REAL RUplMAR(nreg), RUplRLS(nreg)      
      REAL VplMAR (mx,my), VplRLS (LSni,LSnj), VplIRLS (mx,my)
      REAL RVplMAR(nreg), RVplRLS(nreg)   
C +-  -Intermediaires de calcul:
      REAL wk3MAR(mx,my,mz)
      REAL wk3MIX(mx,my,LSnk)
      REAL wk3RLS(LSni,LSnj,LSnk)
      REAL exxpo, dx
C +-  -Miscellaneous
      REAL plint, lplint
C +-  -Indices, ...
      INTEGER kl, ii, jj, nelem, ireg, icheck, izslp
      INTEGER ipl

C +---Ouptuts de controle:
      icheck = 0

C +---Read the MAR coordinate and data
C +   --------------------------------

      CALL UNread
     &   (idMAR, 'pstar', itM, 0, 1,1,mx,my,1,
     &    MARx, MARy, empty1, var_units, pstMAR)
     
      ptop   = 0.01
C +...^MAR top pressure (en kPa= MAR)

      CALL UNread
     &   (idMAR, 'sh', 0, 0, 1,1,mx,my,1,
     &    MARx, MARy, empty1, var_units, shMAR)

      CALL UNread
     &   (idMAR, 'tairDY', itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, tairDY)
C +...  NB: FileID VarName time 0=3D subregion #levels

      CALL UNread
     &   (idMAR, 'qvDY', itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, qvDY)
     
C     CALL UNread
C    &   (idMAR, 'zzDY', itM, 0, 1,1,mx,my,mz,
C    &    MARx, MARy, sigma, var_units, zzDY)
C     **Alternative : test d'apres calc dans MAR

      CALL UNread
     &   (idMAR, 'uairDY', itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, uairDY)
      CALL UNread
     &   (idMAR, 'vairDY', itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, vairDY)

C +...Create CSTp/SIGp hybrid coord equivalent to sigma:
      DO kl=1,mz
        CSTpMA(kl) = ptop
        SIGpMA(kl) = sigma(kl) 
      ENDDO
C +...And define surface:
      CSTpMa(mzz) = ptop
      SIGpMa(mzz) = 1.0

C +---Read the RLS coordinate and data
C +   --------------------------------

      CALL UNread
     &      (idRLS,'CSTp', 0, 0, 0, 0,
     &      LSnk    ,1      , 1     ,
     &      LSlev   ,empty1 , empty1,
     &      tmp_units, wkzRLS)
C +...Change units from Pa-->kPa (MAR units)
      DO kl=1,LSnk
         CSTp(kl) = wkzRLS(kl)* 1.E-3
      ENDDO
      CSTp(LSnk1) = 0.0
         
      CALL UNread
     &      (idRLS,'SIGp', 0, 0, 0, 0,
     &      LSnk     ,1      , 1     ,
     &      LSlev    ,empty1 , empty1,
     &      tmp_units, wkzRLS)
      DO kl=1,LSnk
         SIGp(kl) = wkzRLS(kl)
      ENDDO
      SIGp(LSnk1) = 1.0

      CALL UNread
     &         (idRLS, 'SP', itR, 0, 1,1,LSni,LSnj,1,
     &         LSlon, LSlat, empty1, tmp_units, spRLS)
C +...Change units from Pa-->kPa (MAR units)
      DO ii=1,LSni
      DO jj=1,LSnj
        spRLS(ii,jj)= spRLS(ii,jj)* 1.E-3
      ENDDO
      ENDDO
     
      CALL UNread
     &   (idRLS, 'SH', 0, 0, 1,1,LSni,LSnj,1,
     &    LSlon, LSlat, empty1, tmp_units, shRLS)

      CALL UNread
     &   (idRLS, 'T', itR,  0, 1,1,LSni,LSnj,LSnk,
     &    LSlon, LSlat   , LSlev, var_units, T_RLS)
C +...  NB: FileID VarName time level subregion #levels
C +...  (data, still on the RLS grid)

      CALL UNread
     &   (idRLS, 'Q', itR,  0, 1,1,LSni,LSnj,LSnk,
     &    LSlon, LSlat   , LSlev, var_units, Q_RLS)

      CALL UNread
     &   (idRLS, 'U', itR,  0, 1,1,LSni,LSnj,LSnk,
     &    LSlonU, LSlatU , LSlev, var_units, U_RLS)

      CALL UNread
     &   (idRLS, 'V', itR,  0, 1,1,LSni,LSnj,LSnk,
     &    LSlonV, LSlatV , LSlev, var_units, V_RLS)

      IF (icheck.GE.2) WRITE(*,*) 'Reading ok'
      
C +---Compute geopotentials and pressure of levels.
C +   ---------------------------------------------

C +.. MAR
C +   - -

C +...*Compute pressure and'pkt':
C +... (potential temp divided by 100.[kPa]**(R/Cp))
      DO kl = 1, mz
       DO jj  = 1, my
       DO ii  = 1, mx
         plvMAR(ii,jj,kl)= ptop+pstMAR(ii,jj)*sigma(kl)
         wk3MAR(ii,jj,kl)= tairDY(ii,jj,kl)
     &                 /(plvMAR(ii,jj,kl)**cap)
         lplMAR(ii,jj,kl)=log(plvMAR(ii,jj,kl))
       END DO
       END DO
      END DO  
      
      CALL ZZlevs(wk3MAR, qvDY, shMAR, pstMAR, CSTpMA,SIGpMA,
     &            mx,my,mz, zz_MAR)

C +.. RLS
C +   - -

C +...*Compute pressure and 'pkt':
C +... (potential temp divided by 100.[kPa]**(R/Cp))
      DO kl = 1, LSnk
       DO jj  = 1, LSnj
       DO ii  = 1, LSni
         plvRLS(ii,jj,kl)= CSTp(kl)+spRLS(ii,jj)*SIGp(kl)
         wk3RLS(ii,jj,kl)= T_RLS(ii,jj,kl)
     &                 /(plvRLS(ii,jj,kl)**cap)
         lplRLS(ii,jj,kl)=log(plvRLS(ii,jj,kl))
       END DO
       END DO
      END DO      

      CALL ZZlevs(wk3RLS, Q_RLS, shRLS, spRLS, CSTp, SIGp,
     &            LSni,LSnj,LSnk, zz_RLS)

      IF (icheck.GE.2) WRITE(*,*) 'geop and plevs ok'

C +---Compute and write variables.
C +   ============================

C +---MSLP
C +   ----
      exxpo = - grav / (gamTz * ra)

C +- -MAR
C +   - -
C +...*Reference level: (mz-4 in MAR +/- = LSnk-2)
C +...  mz - 6 should be used if lowest MAR-level is at 2m
      
      IF (REGstudy.EQ.'GRN') THEN
        izslp =  mz - 6 ! Greenland
      ELSEIF (REGstudy.EQ.'ANT') THEN
        izslp =  mz - 4 ! Antar
      ELSE
        izslp =  mz - 4
      ENDIF

      
      DO jj=1,my
      DO ii=1,mx
        slpMAR(ii,jj)= 10.*plvMAR(ii,jj,izslp)
     .    *(1.-gamTz*zz_MAR(ii,jj,izslp)/tairDY(ii,jj,izslp))**exxpo
C +...  (hPa)
      ENDDO
      ENDDO

C +- -RLS
C +   - -
C +...*Reference level:
      izslp =  LSnk - 2

      DO jj=1,LSnj
      DO ii=1,LSni
        slpRLS(ii,jj)= 10.*plvRLS(ii,jj,izslp)
     .    *(1.-gamTz*zz_RLS(ii,jj,izslp)/T_RLS(ii,jj,izslp))**exxpo
C +...  (hPa)
      ENDDO
      ENDDO

C +...*Horizontal interp. of RLS to MAR:
      CALL INThor (-1, LSlon , LSlat , slpRLS,
     &                 MARlon, MARlat, slpIRLS)

      IF (icheck.GE.2) WRITE(*,*) 'slp 1stp ok      '

C +- -Compute & write tM[MSLP] & tSD[MSLP]
C     - - - - - - - - - - - - - - - - - - -

      CALL HDynSTA2D(iOUTnc,istat,idt1,slpMAR,
     &               'tM_slpM','tSD_slpM',tMslpM,tSDslpM)
      CALL HDynSTA2D(iOUTnc,istat,idt1,slpIRLS,
     &               'tM_slpR','tSD_slpR',tMslpR,tSDslpR)

C +- -Compute & write tM[RegM[MSLP]] & tM[RMSE[MSLP]]
C     - - - - - - - - - - - - - - - - - - - - - - - -
      DO jj=1,my
      DO ii=1,mx
        slpBIA(ii,jj)= slpMAR(ii,jj) - slpIRLS(ii,jj)
        slpRME(ii,jj)= slpBIA(ii,jj) * slpBIA(ii,jj) 
      ENDDO
      ENDDO

      CALL HORmean (InReg, slpBIA, RslpBIA)
      CALL HORmean (InReg, slpRME, RslpRME)

      DO ireg = 1,nreg
         RslpRME(ireg)= SQRT(RslpRME(ireg))
      END DO

      CALL HDynSTVAL (iOUTnc,istat,idt1,RslpBIA,
     &               'tMrMslp',tMslpBIA)

      CALL HDynSTVAL (iOUTnc,istat,idt1,RslpRME,
     &               'tMsEslp',tMslpRME)

C +. -Write instant. values of [RegM[MSLP]] & [RMSE[MSLP]]
C     - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL HORmean (InReg, slpMAR,  RslpMAR)
      CALL HORmean (InReg, slpIRLS, RslpRLS)

      CALL UNwrite (iOUTnc, 'rM_slpM', idt1, nreg, 1,1, RslpMAR)
      CALL UNwrite (iOUTnc, 'rM_slpR', idt1, nreg, 1,1, RslpRLS)
      CALL UNwrite (iOUTnc, 'RMSEslp', idt1, nreg, 1,1, RslpRME)

      IF (icheck.GE.2) WRITE(*,*) 'slp finished     '

C +== Begin Interpolation on the surface (2m)
C +   =======================================
C +
C +   En construction !!!! (Xavier Fettweis ?)
C +
C +== END Interpolation on the surface (2m)
C +   =======================================
 
C +== Begin Interpolation on pressure levels.
C +   =======================================
      DO ipl = 1, npl

C +..  Linear interpolation requested pressure levels:

       plint  = plevel(ipl)       !(kPa)
       lplint = log(plint) 
       
C +---Geopotential on pressure levels
C +    -------------------------------
C 
C +...Vertical interp:
      CALL INTsver(lplMAR,zz_MAR,mx, my, mz,    lplint,zplMAR)
      CALL INTsver(lplRLS,zz_RLS,LSni,LSnj,LSnk,lplint,zplRLS)

C +...Horizontal interp. of RLS to MAR:
      CALL INThor (-1, LSlon , LSlat , zplRLS,
     &                 MARlon, MARlat, zplIRLS)

C +...Compute & write tM and tSD:
      CALL HDynSTApl(iOUTnc,istat,idt1,ipl,zplMAR,
     &               'tM_Z_M','tSD_Z_M',tMz_M,tSDz_M)
      CALL HDynSTApl(iOUTnc,istat,idt1,ipl,zplIRLS,     
     &               'tM_Z_R','tSD_Z_R',tMz_R,tSDz_R)
     
C +...Compute & write tSDBP[z500]:
C     - - - - - - - - - - - - - -
      IF (plint.EQ.50) THEN
         CALL HDynSDBP(iOUTnc,istat,idt1,jhuCUR,zplMAR,
     &              'SBPZ5_M',ipfMAR,FilMAR,MBPz_M,SDBPz_M)

         CALL HDynSDBP(iOUTnc,istat,idt1,jhuCUR,zplIRLS,
     &              'SBPZ5_R',ipfRLS,FilRLS,MBPz_R,SDBPz_R)

C +- -Compute & write RMSE[Z500]
C     - - - - - - - - - - - - - -
         DO jj=1,my
         DO ii=1,mx
           zplBIA(ii,jj)= zplMAR(ii,jj) - zplIRLS(ii,jj)
           zplRME(ii,jj)= zplBIA(ii,jj) * zplBIA(ii,jj)
         ENDDO
         ENDDO

C #      CALL HORmean (InReg, zplBIA, RzplBIA)
         CALL HORmean (InReg, zplRME, RzplRME)

         DO ireg = 1,nreg
            RzplRME(ireg)= SQRT(RzplRME(ireg))
         END DO

C        CALL HDynSTVAL (iOUTnc,istat,idt1,RzplBIA,
C    &                  'tMrMzpl',tMzplBIA)
     
C        CALL HDynSTVAL (iOUTnc,istat,idt1,RzplRME,
C    &                  'tMsEzpl',tMzplRME)
C #   
C #      ^^Inutile de s'encombrer. A activer si on veut...

C +. -Write instant. values of [RMSE[Z500]]
C     - - - - - - - - - - - - - - - - - - - -
         CALL UNwrite (iOUTnc, 'RMSE_Z5', idt1, nreg, 1,1, RzplRME)

      ENDIF ! (End 500 hPa -only section)

      IF (icheck.GE.2) WRITE(*,*) 'geop 500 ok      '

C +--  Temperature on pressure levels.
C +    -------------------------------
C +.. .Vertical interp:
       CALL INTsver(lplMAR,tairDY,mx, my, mz,    lplint,TplMAR)
       CALL INTsver(lplRLS,T_RLS ,LSni,LSnj,LSnk,lplint,TplRLS)

C +.. .Horizontal interp. of RLS to MAR:
       CALL INThor(-1, LSlon , LSlat , TplRLS,
     &                 MARlon, MARlat, TplIRLS)

C +.. .Compute & write tM and tSD:
       CALL HDynSTApl(iOUTnc,istat,idt1,ipl,TplMAR,
     &               'tM_T_M','tSD_T_M',tMT_M,tSDT_M)
       CALL HDynSTApl(iOUTnc,istat,idt1,ipl,TplIRLS,
     &               'tM_T_R','tSD_T_R',tMT_R,tSDT_R)

C +. -Write instant. values of RegM[T] (time,plevel)
C     - - - - - - - - - - - - - - - - - - - - - - - -
       CALL HORmean (InReg, TplMAR,  RTplMAR)
       CALL HORmean (InReg, TplIRLS, RTplRLS)
       CALL UNlwrite(iOUTnc, 'rM_T_M', idt1, ipl , nreg,1 , RTplMAR)
       CALL UNlwrite(iOUTnc, 'rM_T_R', idt1, ipl , nreg,1 , RTplRLS)
C +..       ^^^^^^^^ FILEid, VARname ,itime, ilev, Ni  ,Nj, var)
C +..      Writes 1 Level in (time,level,reg) 2D+time variable:
C +..      this requires UNlib version > september 2000. 
C

C +---Specific humidity (Qv) on pressure levels.
C +   ------------------------------------------
C +.. .Vertical interp:
       CALL INTsver(lplMAR,qvDY  ,mx, my, mz,    lplint,QplMAR)
       CALL INTsver(lplRLS,Q_RLS ,LSni,LSnj,LSnk,lplint,QplRLS)

C +.. .Horizontal interp. of RLS to MAR:
       CALL INThor(-1, LSlon , LSlat , QplRLS,
     &                MARlon, MARlat, QplIRLS)

C +.. .Compute & write tM and tSD:
       CALL HDynSTApl(iOUTnc,istat,idt1,ipl,QplMAR,
     &               'tM_Q_M','tSD_Q_M',tMQ_M,tSDQ_M)
       CALL HDynSTApl(iOUTnc,istat,idt1,ipl,QplIRLS,
     &               'tM_Q_R','tSD_Q_R',tMQ_R,tSDQ_R)

C +. -Write instant. values of RegM[Q] (time,plevel)
C     - - - - - - - - - - - - - - - - - - - - - - - - 
       CALL HORmean (InReg, QplMAR,  RQplMAR)
       CALL HORmean (InReg, QplIRLS, RQplRLS)
       CALL UNlwrite(iOUTnc, 'rM_Q_M', idt1, ipl , nreg,1 , RQplMAR)
       CALL UNlwrite(iOUTnc, 'rM_Q_R', idt1, ipl , nreg,1 , RQplRLS)
C +..  .See temperature section above for tech. explanation.    

      IF (icheck.GE.2) WRITE(*,*) 'Qv ok            '


C +---Mean wind vector on pressure levels.
C +   ------------------------------------

      CALL INTsver(lplMAR,uairDY,mx, my, mz,    lplint,UplMAR)
      CALL INTsver(lplMAR,vairDY,mx, my, mz,    lplint,VplMAR)

      CALL INTsver(lplRLS,U_RLS ,LSni,LSnj,LSnk,lplint,UplRLS)
      CALL INTsver(lplRLS,V_RLS ,LSni,LSnj,LSnk,lplint,VplRLS)

      IF (icheck.GE.2) WRITE(*,*) 'Wind 2 models read ok'

C +...Horizontal interp. of RLS to MAR:
      CALL INThor(-1, LSlonU, LSlatU, UplRLS,
     &                MARlon, MARlat, UplIRLS)
      CALL INThor(-1, LSlonV, LSlatV, VplRLS,
     &                MARlon, MARlat, VplIRLS)

      IF (icheck.GE.2) WRITE(*,*) 'Wind 2 models Hor int ok'

C +...Rotates wind vector according to projection >MAR grid:
      dx = (MARx(2)-MARx(1)) * 1000.
      CALL VecRotP (
     &     MARlon, MARlat, dx,
     &     UplIRLS, VplIRLS )

      IF (icheck.GE.2) WRITE(*,*) 'Wind rotate        ok'


C +...Compute & write tM and tSD:
C 
      CALL HDynSTApl(iOUTnc,istat,idt1,ipl,UplMAR,
     &               'tM_U_M','tSD_U_M',tMu_M,tSDu_M)
      CALL HDynSTApl(iOUTnc,istat,idt1,ipl,VplMAR,
     &               'tM_V_M','tSD_V_M',tMv_M,tSDv_M)

      IF (icheck.GE.2) WRITE(*,*) 'Wind stats MAR ok, now RLS'

      CALL HDynSTApl(iOUTnc,istat,idt1,ipl,UplIRLS,
     &               'tM_U_R','tSD_U_R',tMu_R,tSDu_R)
      CALL HDynSTApl(iOUTnc,istat,idt1,ipl,VplIRLS,
     &               'tM_V_R','tSD_V_R',tMv_R,tSDv_R)

       CALL HORmean (InReg, UplMAR,  RUplMAR)
       CALL HORmean (InReg, UplIRLS, RUplRLS)
       CALL HORmean (InReg, VplMAR,  RVplMAR)
       CALL HORmean (InReg, VplIRLS, RVplRLS)
       CALL UNlwrite(iOUTnc, 'rM_U_M', idt1, ipl , nreg,1 , RUplMAR)
       CALL UNlwrite(iOUTnc, 'rM_U_R', idt1, ipl , nreg,1 , RUplRLS)
       CALL UNlwrite(iOUTnc, 'rM_V_M', idt1, ipl , nreg,1 , RVplMAR)
       CALL UNlwrite(iOUTnc, 'rM_V_R', idt1, ipl , nreg,1 , RVplRLS)


C +== END Interpolation on pressure levels.
C +   =======================================
      ENDDO

      IF (icheck.GE.2) WRITE(*,*) 'int on p finished'

C +---Wind norm 500
C +   -------------

C +...Wind Norm:
      DO kl=1,mz
      DO jj=1,my
      DO ii=1,mx
        UV_MAR(ii,jj,kl)=SQRT(uairDY(ii,jj,kl)*uairDY(ii,jj,kl)
     &                      +vairDY(ii,jj,kl)*vairDY(ii,jj,kl))
      ENDDO
      ENDDO
      ENDDO

      DO kl=1,LSnk
      DO jj=1,LSnj
      DO ii=1,LSni
        UV_RLS(ii,jj,kl)=SQRT(U_RLS(ii,jj,kl)*U_RLS(ii,jj,kl)
     &                      +V_RLS(ii,jj,kl)*V_RLS(ii,jj,kl))
      ENDDO
      ENDDO
      ENDDO

C +...Pressure interpolation level:    
      lplint = log(50.) !(kPa)

      CALL INTsver(lplMAR,UV_MAR,mx, my, mz,    lplint,UVpMAR)
      CALL INTsver(lplRLS,UV_RLS,LSni,LSnj,LSnk,lplint,UVpRLS)

C +...Horizontal interp. of RLS to MAR:
      CALL INThor(-1, LSlon , LSlat , UVpRLS,
     &                MARlon, MARlat, UVpIRLS)

C +...Compute & write tM and tSD:
      CALL HDynSTA2D(iOUTnc,istat,idt1,UVpMAR,
     &               'tM_UV5M','tSD_UV5M',tMuvM,tSDuvM)
      CALL HDynSTA2D(iOUTnc,istat,idt1,UVpIRLS,
     &               'tM_UV5R','tSD_UV5R',tMuvR,tSDuvR)

     
      RETURN
      END
