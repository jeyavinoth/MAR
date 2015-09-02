C +----------------------------------------------------------------------+
C | MAR post-processing                                          v3.0.0  |
C |                                                                      |
C | VerP_STAbias (called by VerP)                                        |
C |                                                                      |
C | Computes MAR-LSC biases over subregions (time&space average)         |
C | Principle: horizontal averaging is done on MAR sigma levels for      |
C |            each dataset.                                             |
C |                                                                      |
C | ATTENTION: vous pourriez ne pas être d'accord avec les details de    |
C |            cette méthode ! C'est un "truc" qu'on peut juger plus ou  |
C |            moins astucieux, mais qui a été testé.                    |
C |            Il y a dans l'ordre:                                      |
C |            .interpolation horizontale LSC->MAR (evident)             |    
C |            .interpolation depuis les niv. verticaux LSC              | 
C |                vers des niv. sigma definis comme dans MAR            |
C |                mais sur base de la p au sol LSC (c'est ici le truc:  |
C |                on espère que c'est mieux que passer de suite à       |
C |                la grille MAR, car la moyenne sur sous-région vient   |
C |                encore après et a tendance à diminuer les écarts de   |
C |                relief MAR-LSC)                                       |
C |            .moyenne sur les sous-regions                             |
C |            .interpolation verticale du resultat "LSC" vers les       | 
C |             pressions MAR: le but est d'obtenir des donnees sur les  | 
C |             mêmes niveaux de pression malgré la petite différence    |  
C |             de pression en surface qui subsiste en moyenne           |
C |             entre LSC et MAR.                                        |
C |                                                                      |
C |  --> Le but de cette routine n'est certainement pas de faire         |
C |  une interpolation verticale "intelligente" ni une désaggrégation    |
C |  des données LSC à l'échelle MAR: c'est plutot une aggrégation des   |
C |  deux jeux de données sur les sous-régions. On ne peut donc pas      |
C |  se baser sur ce code pour faire tout autre chose... sorry !         |
C +----------------------------------------------------------------------+
C
      SUBROUTINE VerP_STAbias (idRLS, idMAR, itM, itR, istat, idt1,
     $           InReg, MARlon, MARlat,
     $           spIRLS, CSTp,SIGp, ptop, VPpRLS, VPpMAR,
     $           VarN_MAR, VarN_RLS, VarN_BI, VP_BIAS, iOUTnc)

      IMPLICIT NONE

C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
      
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'

C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM, idRLS, itR, istat, idt1
      CHARACTER*(*) VarN_RLS, VarN_MAR, varN_BI
      INTEGER InReg(mx,my,nreg), iOUTnc
      REAL MARlon(mx,my), MARlat(mx,my)
      REAL spIRLS(mx,my), CSTp (LSnk), SIGp (LSnk), VPpRLS(mz,nreg)
      REAL ptop, VPpMAR(mz,nreg)

C +---INPUT/OUTPUT 
C +   ~~~~~~~~~~~~
      REAL VP_BIAS(mz,nreg)

C +---LOCAL
C +   ~~~~~
      CHARACTER*10 var_units
C +-  -Coordonnees:
      REAL LSlon(LSni), LSlat(LSnj), LSlev(LSnk)
      REAL MARx(mx), MARy(my), empty1(1), sigma(mz)
      REAL ptopPa
C +-  -Valeurs lues et intermédiaires de calcul:
      REAL wk3MAR(mx,my,mz), wkMAR(mx,my)
      REAL wk3MIX(mx,my,LSnk), wkMIX(mx,my)
      REAL wk3RLS(LSni,LSnj,LSnk), wkRLS(LSni,LSnj)
      REAL splocRLS
      REAL valINT      , pINT
      REAL valRLS(LSnk), pRLS(LSnk)
      REAL va3INT(mx,my,mz)
C +-  -Derniere etape, profils vertic. moyens:
      REAL VP_MAR(mz,nreg),VP_IRLS(mz,nreg),SVP_RLS(mz), VP_FRLS
      REAL SVPpRLS(mz), SVPpTMP, VP_TMP(nreg)
C +-  -Indices, ...
      INTEGER ilv, ii, jj, nelem, ireg, icheck

C +---Ouptuts de controle:
      icheck = 1

      
C +-- Initialise computation of mean values
C +   =====================================
      IF (istat.EQ.0) THEN
        DO ireg= 1,nreg 
        DO ilv = 1, mz
          VP_BIAS(ilv,ireg) = 0.0
        ENDDO
        ENDDO
      ENDIF 

      
C +   Read and compute instantaneous bias
C +   ===================================

C +---Read the MAR and RLS data
C +   -------------------------

      CALL UNread
     &   (idMAR, VarN_MAR, itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, wk3MAR)
C +...  NB: FileID,VarName,time,lev,subreg,#levels

      CALL UNread
     &   (idRLS, VarN_RLS, itR,  0, 1,1,LSni,LSnj,LSnk,
     &    LSlon, LSlat, LSlev, var_units, wk3RLS)
C +...  (data, still on the RLS grid)

      IF (icheck.GE.2) WRITE(*,*) 'Reading ok'


C +---Interpolate the RLS data to MAR grid.
C +   -------------------------------------
      
      
C +   Horizontal:

      DO ilv = 1, LSnk
      DO jj  = 1, LSnj
      DO ii  = 1, LSni
        wkRLS(ii,jj)=wk3RLS(ii,jj,ilv)
      END DO
      END DO
        Call INThor (1,
     &    LSlon , LSlat , wkRLS,
     &    MARlon, MARlat, wkMIX)     
C +...  Horizontal interp. to MAR grid 
C +...  (RLS vertical grid)
      DO jj  = 1, my
      DO ii  = 1, mx
        wk3MIX(ii,jj,ilv)=wkMIX(ii,jj)
      END DO
      END DO
      END DO
 
C +  Vertical:
      PtopPa = Ptop*1000.
 
      DO ii= 1,mx
      DO jj= 1,my

         splocRLS = spIRLS(ii,jj)
C +      ...RLS pressure thickness (MAR hor. grid), Pa

         DO ilv = 1, LSnk  
           pRLS  (ilv)= CSTp(ilv)+splocRLS*SIGp(ilv)
           pRLS  (ilv)= log(pRLS(ilv))
C +        ...(log of) pressure on the RLS levels, Pa
           valRLS(ilv)= wk3MIX(ii,jj,ilv)
C +        ...value of "VarN_RLS" on the RLS levels
         END DO

         DO ilv = 1, mz
           pINT= (splocRLS-ptopPa)*sigma(ilv)+ptopPa
           pINT= log(pINT)
C +        ...RLS pres at MAR-type sigma coord., Pa

           CALL INTlin(pRLS, valRLS,LSnk,
     .                 pINT, valINT)
           va3INT(ii,jj,ilv)= valINT
C +        ...Interpolation to MAR vert. grid.

         END DO

      ENDDO
      ENDDO


C +---Sub-region horizontal averaging.
C     --------------------------------
      DO ilv = 1, mz
         DO jj  = 1, my
         DO ii  = 1, mx
            wkMAR(ii,jj)=wk3MAR(ii,jj,ilv)
         END DO
         END DO
         CALL HORmean (InReg, wkMAR, VP_TMP)
         DO ireg = 1,nreg
           VP_MAR(ilv,ireg)= VP_TMP(ireg)
         ENDDO

         DO jj  = 1, my
         DO ii  = 1, mx
            wkMAR(ii,jj)=va3INT(ii,jj,ilv)
         END DO
         END DO
         CALL HORmean (InReg, wkMAR, VP_TMP)
         DO ireg = 1,nreg
            VP_IRLS(ilv,ireg)= VP_TMP(ireg)
         END DO
      ENDDO

C +---Bias computation and outputs
C     ----------------------------

      DO ireg=1,nreg
      
        DO ilv=1,mz    
          SVP_RLS(ilv)=VP_IRLS(ilv,ireg)
          SVPpRLS(ilv)=VPpRLS(ilv,ireg)
          SVPpRLS(ilv)=log(SVPpRLS(ilv))
C +       ...Get the column for current (ireg) region.
C            (using _log_ p is an option for vert interpolation)
        ENDDO

        DO ilv = 1, mz

         SVPpTMP = VPpMAR(ilv,ireg)
         SVPpTMP = log(SVPpTMP)
C +      ...Get the (log of) MAR levels pressure for current region

          CALL INTlin(SVPpRLS, SVP_RLS, mz,
     &                SVPpTMP, VP_FRLS)
C +     ...Legere interp. lineaire en coord press. = ramener au
C +     ...niveaux de p du MAR du a petite diff p surf. moyenne
C
C       Pour test sans cette interpolation (pas recommande):
C       VP_FRLS = VP_IRLS(ilv)

         VP_BIAS(ilv,ireg)=VP_BIAS(ilv,ireg)+VP_MAR(ilv,ireg)-VP_FRLS

         IF (SVPpTMP.GT.SVPpRLS(mz)
     &   .OR.SVPpTMP.LT.SVPpRLS(1 ) ) THEN
             VP_BIAS(ilv,ireg)= 1.E28
C +         ...do not allow extrapolation (set as missing)
         ENDIF
  
       ENDDO

      ENDDO
      


C +--Final computation & write result (when istat = 2)
C +  =================================================
      IF (istat.EQ.2) THEN

       DO ireg= 1,nreg
       DO ilv = 1, mz
        VP_BIAS(ilv,ireg) = VP_BIAS(ilv,ireg)/idt1
       ENDDO
       ENDDO

       CALL UNwrite (iOUTnc, VarN_BI, 1, mz, nreg, 1, VP_BIAS )

      ENDIF

      RETURN
      END
