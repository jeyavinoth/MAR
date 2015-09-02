C +----------------------------------------------------------------------+
C | MAR post-processing                                          2002.02 |
C | Station                                                              |
C | Computes vertical profiles of a RLS variables                        |
c  attention, ce code est incorrect pour donnees LSC (RLS)
c  et se limite a une interp vertic pour MAR.
c
C +----------------------------------------------------------------------+
      SUBROUTINE Station (idRLS, idMAR, itM, itR,
     $           MARlon, MARlat,
     $           VarN_MAR, IN_RLS, ii_in, jj_in,level, 
     $           OUT_MAR, OUT_RLS)

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
      INTEGER idMAR, itM, idRLS, itR, ii_in, jj_in
      CHARACTER*(*) VarN_MAR
      REAL level
      REAL IN_RLS(mx,my,mz)
      REAL MARlon(mx,my), MARlat(mx,my)
      
C +---OUTPUT 
C +   ~~~~~~
      REAL OUT_MAR, OUT_RLS
      
C +---LOCAL
C +   ~~~~~
      CHARACTER *10 var_units, tmp_units
C +-  -Coordonnees:
      REAL LSlon(LSni), LSlat(LSnj), LSlev(LSnk)
      REAL MARx(mx), MARy(my), empty1(1), sigma(mz)
      REAL ptop, ptopPa
      REAL CSTp (LSnk), SIGp (LSnk)
      INTEGER ilv
C +-  -Valeurs lues et intermÈdiaires de calcul:
      REAL wk3MAR(mx,my,mz), wkMAR(mx,my)
      REAL wk3MIX(mx,my,LSnk), wkMIX(mx,my)
      REAL wk3RLS(LSni,LSnj,LSnk), wkRLS(LSni,LSnj)
      REAL pstMAR(mx,my), spIRLS(mx,my), splocRLS
      REAL valINT, pINT,VAR_RLS(mz),VAR_RLS_corr(mz)
      REAL valRLS(LSnk), pRLS(LSnk)
      REAL va3INT(mx,my,mz)
      REAL VlplMAR(mz),VlplRLS(mz)
C +-  -Derniere etape, profils vertic. moyens:
      REAL var_MAR(mz),VP_IRLS(mz),SVP_RLS(mz), VP_FRLS
      REAL VPpMAR(mz),VPpRLS(mz),SVPpRLS(mz)
      REAL SVPpTMP, VP_TMP
      
C +---Read the MAR coordinate and data
C +   --------------------------------

      CALL UNread
     &   (idMAR, 'pstar', itM, 0, 1,1,mx,my,1,
     &    MARx, MARy, empty1, var_units, pstMAR)
     
      ptop   = 0.01
      ptopPa = ptop * 1000.
C +   ...MAR top pressure (en kPa = MAR, et en Pa)

      CALL UNread
     &   (idMAR, VarN_MAR, itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, wk3MAR)
C +...  NB: FileID VarName time 0=3D subregion #levels

C +---Read the RLS coordinate and data
C +   --------------------------------

      CALL UNread
     &   (idRLS,'CSTp', 0, 0, 0, 0,
     &    LSnk    ,1      , 1     ,
     &    LSlev   ,empty1 , empty1,
     &    tmp_units, CSTp)
C +...constant pressure contrib. to coordinate, in Pa

      CALL UNread
     &   (idRLS,'SIGp', 0, 0, 0, 0,
     &    LSnk     ,1      , 1     ,
     &    LSlev    ,empty1 , empty1,
     &    tmp_units, SIGp)

      CALL UNread
     &   (idRLS, 'SP', itR, 0, 1,1,LSni,LSnj,1,
     &    LSlon, LSlat, empty1, tmp_units, wkRLS)
     
      CALL INThor (1,
     &    LSlon , LSlat , wkRLS,
     &    MARlon, MARlat, spIRLS)
     
C +---Mean Pressure Coordinates
C +   -------------------------

c      ce code et ce qui suit vien de VerP. Ca
c       n'est pas approprié ici
c      (il n'y a d'ailleurs pas de pression moyenne !)
c
C +-  -For MAR:
C     - - - - -
       DO ilv = 1, mz
        VPpMAR(ilv)      = 10.*(pstMAR(ii_in,jj_in)*sigma(ilv)+ptop)
        VlplMAR(ilv)= log(VPpMAR(ilv))
C +     ...En hPa.
       ENDDO
 
C +-  -For RLS interpolated to MAR levels:
C     - - - - - - - - - - - - - - - - - - -
C +    ...(RLS est en Pa)
       DO ilv = 1, mz
         VPpRLS(ilv)= 0.01
     &             *((spIRLS(ii_in,jj_in)-ptopPa)*sigma(ilv)+ptopPa)
        VlplRLS(ilv)= log(VPpRLS(ilv))
C +       ...Pres. surf. RLS et niv. MAR, en hPa.
      ENDDO
          
          
C +---Station value
C     -------------
 
      DO ilv = 1, mz
       var_MAR(ilv)= wk3MAR(ii_in,jj_in,ilv)
       var_RLS(ilv)= IN_RLS(ii_in,jj_in,ilv)
      END DO
      	
      DO ilv = 1, mz      	   
       CALL INTlin(VPpRLS, var_RLS, mz,
     &               VPpMAR(ilv), var_RLS_corr(ilv))           
      ENDDO



C +---Vertical int.
C     -------------
        IF (level.lt.VlplMAR(mz)) then
         CALL INTsver(VlplMAR,var_MAR,1, 1, mz,level,Out_MAR)
         CALL INTsver(VlplMAR,var_RLS_corr,1, 1, mz,level,Out_RLS)
c PhM je ne vois pas pourquoi utiliser INTsver
c     INTsver est plus ou moins la m chose que INTlin
c     mais pour une grille horiz complète.
        ELSE
         Out_MAR = var_MAR(mz)
         Out_RLS = var_RLS(mz)
        END if
      END
                         