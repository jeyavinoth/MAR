C +----------------------------------------------------------------------+
C | MAR post-processing                                                  |
C |                                                              v 3.0.0 |
C |                                                                      |
C |  Computes vertical profiles, by horizontal averaging                 |
C |  on the MAR sigma levels.                                            |
C +----------------------------------------------------------------------+
      SUBROUTINE VerP (idRLS, idMAR, itM, itR, istat,
     $                 InReg, MARlon, MARlat,
     $                 iOUTnc,idt1)

      IMPLICIT NONE

C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
      
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'LSMphy.inc'
      
C +...* Usefull global things
      include 'globals.inc'

C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM, idRLS, itR
      INTEGER InReg(mx,my,nreg)
      REAL MARlon(mx,my), MARlat(mx,my)
      INTEGER istat,iOUTnc,idt1

C +---OUTPUT 
C +   ~~~~~~
C     (none)

C +---LOCAL
C +   ~~~~~
      CHARACTER *10  var_units, tmp_units
C +-  -Coordonnees:
      REAL LSlon(LSni), LSlat(LSnj)
      REAL MARx(mx), MARy(my), empty1(1), sigma(mz)
      REAL ptop, ptopPa
      REAL CSTp (LSnk1), SIGp (LSnk1)
      
C +-  -Valeurs lues et intermédiaires de calcul:
      REAL pstMAR(mx,my), spIRLS(mx,my), splocRLS, wkRLS(LSni,LSnj)
      REAL valINT      , pINT

C +-  -Valeurs calculées
      REAL pM_VP(mz,nreg),VPpOUT(mz,nreg)
      REAL t__VP(mz,nreg), qv_VP(mz,nreg),u__VP(mz,nreg), v__VP(mz,nreg)
      SAVE pM_VP, t__VP, qv_VP, u__VP, v__VP, VPpOUT

C +-  -Derniere etape, profils vertic. moyens:
      REAL VPpMAR(mz,nreg),VPpRLS(mz,nreg)
      REAL pm_reg(nreg)
      
C +-  -Indices, ...
      INTEGER ilv, ii, jj, nelem, ireg, icheck

C +---Ouptuts de controle:
      icheck = 1


C +---Read and compute pressure coordinates for all variables
C     =======================================================
C     (prepare computation of MAR-LSC biases)

C +---Read the MAR coordinate
C +   -----------------------
      CALL UNsread (idMAR,'level',0,0,1,1,
     &              mz,1,1, tmp_units, sigma)

      CALL UNread
     &   (idMAR, 'pstar', itM, 0, 1,1,mx,my,1,
     &    MARx, MARy, empty1, var_units, pstMAR)
     
      ptop   = 0.01
      ptopPa = ptop * 1000.
C +   ...MAR top pressure (en kPa = MAR, et en Pa)
C +   This should be read in files, but as we don't have
C     it until now...


C +---Read the RLS coordinate
C +   -----------------------

      CALL UNsread (idRLS,'CSTp',0,0,1,1,
     &              LSnk,1,1, tmp_units, CSTp)
C +...constant pressure contrib. to coordinate, in Pa

      CALL UNsread (idRLS,'SIGp',0,0,1,1,
     &              LSnk,1,1, tmp_units, SIGp)

      CALL UNread
     &   (idRLS, 'SP', itR, 0, 1,1,LSni,LSnj,1,
     &    LSlon, LSlat, empty1, tmp_units, wkRLS)
     
      CALL INThor (1,
     &    LSlon , LSlat , wkRLS,
     &    MARlon, MARlat, spIRLS)

      IF (icheck.GE.2) WRITE(*,*) 'Pressure reading ok'

C +---Mean Pressure Coordinates
C +   -------------------------

C +-  -For MAR:
C     - - - - -
       CALL HORmean (InReg, pstMAR, pm_reg)
       DO ireg= 1,nreg
        DO ilv = 1, mz
          VPpMAR(ilv,ireg)=10.*(pm_reg(ireg)*sigma(ilv)+ptop)
C +       ...En hPa.
        ENDDO
       ENDDO
 
C +-  -For RLS interpolated to MAR levels:
C     - - - - - - - - - - - - - - - - - - -
       CALL HORmean (InReg, spIRLS, pm_reg)
C +    ...(RLS est en Pa)
       DO ireg= 1,nreg
        DO ilv = 1, mz
          VPpRLS(ilv,ireg)= 0.01
     &          *((pm_reg(ireg)-ptopPa)*sigma(ilv)+ptopPa)
C +       ...Pres. surf. RLS et niv. MAR, en hPa.
        ENDDO
       ENDDO

      IF (icheck.GE.2) WRITE(*,*) 'Mean press coord ok'
 
C +-  Initialise the computation of mean pressure
C +   ===========================================
      IF (istat.EQ.0) THEN
        DO ireg= 1,nreg 
        DO ilv = 1, mz
          VPpOUT (ilv,ireg) = 0.0
        ENDDO
        ENDDO
      ENDIF 


C +-  Compute & write profiles for MAR-RLS biases of the
C      main variables...
C     ==================================================
      
      CALL VerP_STAbias(idRLS, idMAR, itM, itR, istat, idt1,
     $             InReg, MARlon, MARlat,
     $             spIRLS, CSTp,SIGp, ptop, VPpRLS, VPpMAR,
     $             'tairDY', 'T', 't__VP ', t__VP, iOUTnc)

      CALL VerP_STAbias(idRLS, idMAR, itM, itR, istat, idt1,
     $             InReg, MARlon, MARlat,
     $             spIRLS, CSTp,SIGp, ptop, VPpRLS, VPpMAR,
     $             'qvDY'  , 'Q', 'qv_VP ', qv_VP, iOUTnc)

      CALL VerP_STAbias(idRLS, idMAR, itM, itR, istat, idt1,
     $             InReg, MARlon, MARlat,
     $             spIRLS, CSTp,SIGp, ptop, VPpRLS, VPpMAR,
     $             'uairDY', 'U', 'u__VP ', u__VP, iOUTnc)

      CALL VerP_STAbias(idRLS, idMAR, itM, itR, istat, idt1,
     $             InReg, MARlon, MARlat,
     $             spIRLS, CSTp,SIGp, ptop, VPpRLS, VPpMAR,
     $             'vairDY', 'V', 'v__VP ', v__VP, iOUTnc)


C +--Compute and write mean pressure at the levels
C +  ==============================================

      DO ireg= 1,nreg
       DO ilv = 1, mz
         VPpOUT (ilv,ireg)=VPpOUT(ilv,ireg) +VPpMAR(ilv,ireg) 
       ENDDO
      ENDDO

      IF (istat.EQ.2) THEN

        DO ireg= 1,nreg
        DO ilv = 1, mz
         VPpOUT (ilv,ireg) = VPpOUT (ilv,ireg)/idt1
        ENDDO
        ENDDO

        CALL UNwrite (iOUTnc, 'pM_VP', 1, mz, nreg, 1, VPpOUT)

      ENDIF

      END
