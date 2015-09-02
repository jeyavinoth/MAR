C +----------------------------------------------------------------------+
C | MAR post-processing                                          2001    |
C |                                                              v.2.2   |
C | HDyn                                                                 |
C |  Horizontal, mostly "dynamics oriented" maps & stat. values.         |
c
c
c
c   ATTENTION: version provisoire, largment a revoir
c
c
C +----------------------------------------------------------------------+
C
      SUBROUTINE SurSta(idRLS, idMAR, itM, itR, istat,
     $                wmSTA, MARlon, MARlat, jhuCUR,
     $                ioutNC,idt1)

      IMPLICIT NONE
C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
            
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'globals.inc'
      include 'LSMphy.inc'

C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM, idRLS, itR
      REAL    MARlon(mx,my), MARlat(mx,my)
      INTEGER istat,jhuCUR,ioutNC,idt1
      REAL wmSTA(mx,my,nreg)

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
C +-  -For the stations values :
      INTEGER ii_ml, jj_ml, stop_ml, stop_ml_old
      Logical onetime
      REAL Ti3DRLS(mx,my,mz), Ui3DRLS(mx,my,mz)
      REAL Vi3DRLS(mx,my,mz), Qi3DRLS(mx,my,mz) 
      REAL Ui2DRLS(mx,my)   , Vi2DRLS(mx,my)      
      REAL TmlMAR(nreg,nml) , TmlRLS(nreg,nml) 
      REAL QmlMAR(nreg,nml) , QmlRLS(nreg,nml) 
      REAL UmlMAR(nreg,nml) , UmlRLS(nreg,nml) 
      REAL VmlMAR(nreg,nml) , VmlRLS(nreg,nml) 
      REAL SmlMAR(nreg,nml) , SmlRLS(nreg,nml) 
      REAL DmlMAR(nreg,nml) , DmlRLS(nreg,nml) 
      REAL TmlMAR_tmp(nreg) , TmlRLS_tmp(nreg)
      REAL QmlMAR_tmp(nreg) , QmlRLS_tmp(nreg)
      REAL UmlMAR_tmp(nreg) , UmlRLS_tmp(nreg)
      REAL VmlMAR_tmp(nreg) , VmlRLS_tmp(nreg)
      REAL SmlMAR_tmp(nreg) , SmlRLS_tmp(nreg) 
      REAL DmlMAR_tmp(nreg) , DmlRLS_tmp(nreg) 
C +-  -Intermediaires de calcul:
      REAL wk3MAR(mx,my,mz)
      REAL wk3MIX(mx,my,LSnk)
      REAL wk3RLS(LSni,LSnj,LSnk)
      REAL exxpo, dx
C +-  -Miscellaneous
      REAL plint, lplint
C +-  -Indices, ...
      INTEGER kl, kk, ii, jj, nelem, ireg, icheck, izslp
      INTEGER ipl, iml

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


C +== Begin Interpolation on the surface (for the meteo stations).
C +   ============================================================
C +
      onetime = .true.

      DO ireg = 1, nreg          
       
       stop_ml = 0
       ii_ml   = 10
       jj_ml   = 10
              
       DO ii = 1,mx
       DO jj = 1,my
         stop_ml_old = stop_ml
         stop_ml     = stop_ml + wmSTA(ii,jj,ireg)
        if (stop_ml.eq.1.and.stop_ml_old.eq.0) then
         ii_ml = ii
         jj_ml = jj	
        end if	
        if (stop_ml.gt.2) then
         go to 101	
        end if	                 
       END DO
       END DO

       if (stop_ml.eq.0) then
        write (*,*) 'No definition for the region/station n°:', ireg
        go to 101	
       end if	
       
       if (onetime) then              

C +--  RLS var to MAR grid
C +    -------------------

       call Int3D (idRLS, idMAR, itM, itR,
     .        MARlon, MARlat, 'T',Ti3DRLS)
       call Int3D (idRLS, idMAR, itM, itR,
     .        MARlon, MARlat, 'Q',Qi3DRLS)
       call Int3D (idRLS, idMAR, itM, itR,
     .        MARlon, MARlat, 'U',Ui3DRLS)
       call Int3D (idRLS, idMAR, itM, itR,
     .        MARlon, MARlat, 'V',Vi3DRLS)  

       DO kk =1, mz
C +...Rotates wind vector according to projection >MAR grid:
              
       DO ii = 1,mx
       DO jj = 1,my
        Ui2DRLS(ii,jj) = Ui3DRLS(ii,jj,kk)
        Vi2DRLS(ii,jj) = Vi3DRLS(ii,jj,kk)
       END DO
       END DO
       	
       dx = (MARx(2)-MARx(1)) * 1000.
       CALL VecRotP (
     &     MARlon, MARlat, dx,
     &     Ui2DRLS, Vi2DRLS)
       
       DO ii = 1,mx
       DO jj = 1,my
        Ui3DRLS(ii,jj,kk) = Ui2DRLS(ii,jj)
        Vi3DRLS(ii,jj,kk) = Vi2DRLS(ii,jj)
       END DO
       END DO
       
       END DO          
       	
       onetime = .false.
       
       end if

C +--  Metre levels
C +    ------------
                 
       DO iml = 1, nml
              
       plint  =   plvMAR(ii_ml,jj_ml,mz)*10.
     .        *((tairDY(ii_ml,jj_ml,mz)
     .        -   gamTz *(zz_MAR(ii_ml,jj_ml,mz)-shMAR(ii_ml,jj_ml))
     .        +   gamTz * mlevel(iml)) / tairDY(ii_ml,jj_ml,mz))**exxpo
     
       IF(zz_MAR(ii_ml,jj_ml,mz)-shMAR(ii_ml,jj_ml).gt.mlevel(iml))then	
       plint  =   plvMAR(ii_ml,jj_ml,mz)*10.
       write (*,113) zz_MAR(ii_ml,jj_ml,mz)-shMAR(ii_ml,jj_ml),
     .               mlevel(iml),ireg
 113   format('-----------------------------------------------',/,
     .        'First metre level must be GE first level of MAR',/,
     .        'Then, first metre level =', f7.3,' m in place of',/,
     .        f5.2,' m for the station ',i2,/, 
     .        '-----------------------------------------------')
       END IF
       
       lplint = log(plint)       

C +--  Temperature on metre levels.
C +    ----------------------------

       call station(idRLS, idMAR, itM, itR,
     $               MARlon, MARlat, 'tairDY',Ti3DRLS,
     $               ii_ml,jj_ml,lplint,
     $               TmlMAR(ireg,iml), TmlRLS(ireg,iml))  
       
C +---Specific humidity (Qv) on  metre levels.
C +   ------------------------------------------

       call station(idRLS, idMAR, itM, itR,
     $               MARlon, MARlat, 'qvDY',Qi3DRLS,
     $               ii_ml,jj_ml,lplint,
     $               QmlMAR(ireg,iml), QmlRLS(ireg,iml))     

C +---Mean wind vector on metre levels.
C +   ----------------------------------

       call station(idRLS, idMAR, itM, itR,
     $               MARlon, MARlat, "uairDY",Ui3DRLS,
     $               ii_ml,jj_ml,lplint,
     $               UmlMAR(ireg,iml), UmlRLS(ireg,iml))    

       call station(idRLS, idMAR, itM, itR,
     $               MARlon, MARlat, "vairDY",Vi3DRLS,
     $               ii_ml,jj_ml,lplint,
     $               VmlMAR(ireg,iml), VmlRLS(ireg,iml))    

C +---Wind direction and speed according to MAR lat & lon:
C +   ----------------------------------------------------  
          
       call wind(UmlMAR(ireg,iml),VmlMAR(ireg,iml),
     .           ii_ml, jj_ml, dx, MARlat, MARlon,
     .           DmlMAR(ireg,iml),SmlMAR(ireg,iml))      

       call wind(UmlRLS(ireg,iml),VmlRLS(ireg,iml),
     .           ii_ml, jj_ml, dx, MARlat, MARlon,
     .           DmlRLS(ireg,iml),SmlRLS(ireg,iml))      
             
       END DO                  
         
  101 Continue                         
      END DO
      	
       DO iml = 1, nml
       
        DO ireg = 1, nreg
        TmlMAR_tmp(ireg) = TmlMAR(ireg,iml)
        TmlRLS_tmp(ireg) = TmlRLS(ireg,iml)       
        QmlMAR_tmp(ireg) = QmlMAR(ireg,iml)
        QmlRLS_tmp(ireg) = QmlRLS(ireg,iml) 
        UmlMAR_tmp(ireg) = UmlMAR(ireg,iml)
        UmlRLS_tmp(ireg) = UmlRLS(ireg,iml)                
        VmlMAR_tmp(ireg) = VmlMAR(ireg,iml)
        VmlRLS_tmp(ireg) = VmlRLS(ireg,iml)
        SmlMAR_tmp(ireg) = SmlMAR(ireg,iml)
        SmlRLS_tmp(ireg) = SmlRLS(ireg,iml)                 
        DmlMAR_tmp(ireg) = DmlMAR(ireg,iml)
        DmlRLS_tmp(ireg) = DmlRLS(ireg,iml)  
        END DO

        CALL UNlwrite(iOUTnc,'sM_T_M',idt1,iml,nreg,1,TmlMAR_tmp)
        CALL UNlwrite(iOUTnc,'sM_T_R',idt1,iml,nreg,1,TmlRLS_tmp)
        CALL UNlwrite(iOUTnc,'sM_Q_M',idt1,iml,nreg,1,QmlMAR_tmp)
        CALL UNlwrite(iOUTnc,'sM_Q_R',idt1,iml,nreg,1,QmlRLS_tmp)            
        CALL UNlwrite(iOUTnc,'sM_U_M',idt1,iml,nreg,1,UmlMAR_tmp)
        CALL UNlwrite(iOUTnc,'sM_U_R',idt1,iml,nreg,1,UmlRLS_tmp)
        CALL UNlwrite(iOUTnc,'sM_V_M',idt1,iml,nreg,1,VmlMAR_tmp)
        CALL UNlwrite(iOUTnc,'sM_V_R',idt1,iml,nreg,1,VmlRLS_tmp) 
        CALL UNlwrite(iOUTnc,'sM_S_M',idt1,iml,nreg,1,SmlMAR_tmp)
        CALL UNlwrite(iOUTnc,'sM_S_R',idt1,iml,nreg,1,SmlRLS_tmp)
        CALL UNlwrite(iOUTnc,'sM_D_M',idt1,iml,nreg,1,DmlMAR_tmp)
        CALL UNlwrite(iOUTnc,'sM_D_R',idt1,iml,nreg,1,DmlRLS_tmp)   
       END DO

C +
C +== END Interpolation on the surface (2m)
C +   =======================================
 


     
      RETURN
      END
