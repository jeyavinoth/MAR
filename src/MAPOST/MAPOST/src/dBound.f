C +----------------------------------------------------------------------+
C | MAR post-processing                                           02.2001|
C |                                                              v.1.0   |
C | dBound                                                               |
C |  Computes statistics using the distance to boundary as indep. var.   |
C |  Ex: MKE(dbound) => diagnostic for LBC treatment                     |
C +----------------------------------------------------------------------+
C
      SUBROUTINE dBound(idRLS, idMAR, itR, itM, istat, ioutNC, 
     &                  slpMAR,idt1, iptMKE)

      IMPLICIT NONE
 
C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
      
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'LSMphy.inc'

C +...* Internal values of dBound (retain for next call):
      include 'dBound.inc'

C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM, idRLS, itR
      INTEGER istat,ioutNC,idt1, iptMKE
      REAL slpMAR(mx,my)

C +---OUTPUT 
C +   ~~~~~~
C     (none)

C +---LOCAL
C +   ~~~~~
      CHARACTER *10 var_units, tmp_units
C +-  -Coordonnees:
      REAL LSlon(LSni), LSlat(LSnj), LSlev(LSnk)
      REAL CSTp (LSnk1), SIGp (LSnk1), wkzRLS(LSnk)
      REAL MARx(mx), MARy(my), empty1(1), sigma(mz), sigmid(mz)
      REAL ptop, CSTpMA(mzz), SIGpMA(mzz)
C +-  -Valeurs lues
      REAL uairDY(mx,my,mz), vairDY(mx,my,mz)
      REAL pstar(mx,my)
      REAL U_RLS(LSni,LSnj,LSnk)
      REAL qwHY(mx,my,mz), qiHY(mx,my,mz)
C +-  -Valeurs Calculees
      REAL uLarge (mx,my), vLarge(mx,my)
      REAL uMeso (mx,my)   , vMeso (mx,my)
      REAL grMKE (mx,my), slpMeso(mx,my)
      REAL grCLOU(mx,my), dsig(mz)

C +-  -Intermediaires de calcul:
      REAL tnorm
C +-  -Miscellaneous
      LOGICAL taksqr
C +-  -Indices, ...
      INTEGER kl, ii, jj, icheck

C +---Ouptuts de controle:
      icheck = 0

      IF(icheck.GE.2)WRITE(*,*) 'Begin dBound'

C +---Read the MAR coordinate and data
C +   --------------------------------

      CALL UNread
     &   (idMAR, 'uairDY', itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, uairDY)
      CALL UNread
     &   (idMAR, 'vairDY', itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, vairDY)

C +---Read the LS coordinate and data
C +   --------------------------------
C     CALL UNread
C    &   (idRLS, 'U', itR, 0, 1,1,LSni,LSnj,LSnk,
C    &    LSlon, LSlat, LSlev, var_units, U_RLS)

      
C +---Compute and write variables.
C +   ============================

C +---MKE. (Giorgi 1993, MWR121)
C +   ----
      iptMKE = 2    

      DO kl = 1,mz

C +     Spatial running average => "mesoscale" comps. of u,v
C +
        CALL  MesoFilt (uairDY, uMeso, kl,iptMKE)
        CALL  MesoFilt (vairDY, vMeso, kl,iptMKE)
      
C +     MKE (instantaneous) + recomp u,v filtered => in/outflow 

        DO jj=1,my
        DO ii=1,mx
          grMKE(ii,jj)= 0.5*( (uMeso(ii,jj)*uMeso(ii,jj)) 
     &                       +(vMeso(ii,jj)*vMeso(ii,jj)) )
          uLarge(ii,jj) = uairDY(ii,jj,kl) - uMeso(ii,jj)
          vLarge(ii,jj) = vairDY(ii,jj,kl) - vMeso(ii,jj)
        ENDDO      
        ENDDO

C +     Compute MKE (dbound) + write

        CALL dBndOfSTA(ioutNC,istat, grMKE, uLarge, vLarge,
     &        kl, idt1, 'MKE2',tM_MKE,tMiMKE,tMoMKE,NPiMKE,NPoMKE)

      ENDDO
     
C +---Temporary addition for LBC paper: 2nd MKE output = 3*2+1 pt avera. 

      iptMKE = 3

      DO kl = 1,mz

C +     Spatial running average => "mesoscale" comps. of u,v
C +
        CALL  MesoFilt (uairDY, uMeso, kl,iptMKE)
        CALL  MesoFilt (vairDY, vMeso, kl,iptMKE)

C +     MKE (instantaneous)

        DO jj=1,my
        DO ii=1,mx
          grMKE(ii,jj)= 0.5*( (uMeso(ii,jj)*uMeso(ii,jj))
     &                       +(vMeso(ii,jj)*vMeso(ii,jj)) )
        ENDDO     
        ENDDO

C +     Compute MKE (dbound) + write

        CALL dBoundSTA(ioutNC,istat, grMKE, kl, idt1,
     &                'MKE3_db',tM_MKE3)

      ENDDO

C +---Temporary addition for LBC paper: 2-3 dx noise

      iptMKE = 1

      DO kl = 1,mz

C +     Spatial running average => "mesoscale" comps. of u,v
C +
        CALL  MesoFilt (uairDY, uMeso, kl,iptMKE)
        CALL  MesoFilt (vairDY, vMeso, kl,iptMKE)

C +     MKE (instantaneous)

        DO jj=1,my
        DO ii=1,mx
          grMKE(ii,jj)= 0.5*( (uMeso(ii,jj)*uMeso(ii,jj))
     &                       +(vMeso(ii,jj)*vMeso(ii,jj)) )
        ENDDO
        ENDDO

C +     Compute MKE (dbound) + write

        CALL dBoundSTA(ioutNC,istat, grMKE, kl, idt1,
     &                'MKE1_db',tM_MKE1)
      ENDDO

C +---MSLP
C +   ----
        iptMKE=2  
        CALL MesoF2D (slpMAR, slpMeso, iptMKE)

        DO jj=1,my
        DO ii=1,mx
          slpMeso(ii,jj)= slpMeso(ii,jj)*slpMeso(ii,jj)
        ENDDO
        ENDDO
        
        taksqr = .TRUE.
        CALL dBnd2DSTA(ioutNC,istat, slpMeso, idt1,
     &                'MSLP_db',tM_SPdb,taksqr)

C +---Clouds - integrated column Qw + Qi
C +   -----------------------------------

C +     IF variable is not available, warn and set=0.0:
C +     (NOVAR_WARN level 1 = replace and warn, don't stop)
        CALL UNparam('NOVAR_REPLACE',0.0)
        CALL UNparam('NOVAR_WARNING',1.0)

        CALL UNsread
     &   (idMAR,'pstar' ,itM,0,
     &    1,1,mx,my,1,var_units,pstar)

        CALL UNsread
     &   (idMAR, 'qwHY',itM,0,
     &    1,1,mx,my,mz,var_units,qwHY) !Cloud water

        CALL UNsread
     &   (idMAR,'qiHY' ,itM,0,
     &    1,1,mx,my,mz,var_units,qiHY) !Cloud ice  

        CALL UNsread
     &   (idMAR,'level2',1,1,
     &    1,1,mz,1,1,var_units,sigmid)
C           == sigma at levels k+1/2 (1/2 closer to surface)

C +     Go back to standard treatment of missing variables:
C +     (warn level 2 = standard = stop all)
        CALL UNparam('NOVAR_WARNING',2.0)


C       fmult = 1000./gravit => Kg/m2 en eau 1 Kg = 1mm*1m2 
C       => pour m(water):   1/g           

        DO kl=1,mz-1
         dsig(kl)=sigmid(kl+1)-sigmid(kl)
         DO jj=1,my
         DO ii=1,mx
c         grCLOU(ii,jj)= 
c    &    (qwHY(ii,jj,kl)+qiHY(ii,jj,kl))*pstar(ii,jj)*dsig(kl)/9.81
         ENDDO
         ENDDO
        ENDDO

        taksqr=.FALSE.
c       CALL dBnd2DSTA(ioutNC,istat, grCLOU, idt1,
c    &                'CLOU_db',tM_CLdb,taksqr)
C       write(*,*) 'Can t output CLOU_db for BIOCLIM now...'

      RETURN
      END

C +-------------------------------------------------------------
      SUBROUTINE MesoFilt (var, varMeso, kl, iptMKE)

      IMPLICIT NONE

C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'

      INTEGER ii, jj, ia, ja, kl, iptMKE
      REAL var(mx,my,mz)
      REAL varMeso (mx,my)
      REAL tnorm, tmpFilt

C +...*Note: iptMKE is the /2 size of the running mean area
      tnorm= (iptMKE*2)+1
      tnorm= tnorm*tnorm

      DO jj=1,my
      DO ii=1,mx
        varMeso(ii,jj)=0.0
      ENDDO
      ENDDO     

      DO jj=1+iptMKE,my-iptMKE
      DO ii=1+iptMKE,mx-iptMKE
          tmpFilt = 0.0
          DO ja= jj-iptMKE, jj+iptMKE
          DO ia= ii-iptMKE, ii+iptMKE
           tmpFilt = tmpFilt + var(ia,ja,kl)
          ENDDO
          ENDDO  
          varMeso(ii,jj)= var(ii,jj,kl)-(tmpFilt/tnorm)
      ENDDO
      ENDDO  

      END 

      SUBROUTINE MesoF2D (var, varMeso, iptMKE)

      IMPLICIT NONE

C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'

      INTEGER ii, jj, ia, ja, iptMKE
      REAL var(mx,my)
      REAL varMeso (mx,my)
      REAL tnorm, tmpFilt

C +...*Note: iptMKE is the /2 size of the running mean area
      tnorm= (iptMKE*2)+1
      tnorm= tnorm*tnorm

      DO jj=1,my
      DO ii=1,mx
        varMeso(ii,jj)=0.0
      ENDDO
      ENDDO

      DO jj=1+iptMKE,my-iptMKE
      DO ii=1+iptMKE,mx-iptMKE
          tmpFilt = 0.0
          DO ja= jj-iptMKE, jj+iptMKE
          DO ia= ii-iptMKE, ii+iptMKE
           tmpFilt = tmpFilt + var(ia,ja)
          ENDDO
          ENDDO
          varMeso(ii,jj)= var(ii,jj)-(tmpFilt/tnorm)
      ENDDO
      ENDDO

      END

