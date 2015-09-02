C +----------------------------------------------------------------------+
C | MAR post-processing                                          02/1999 |
C |                                                                      |
C | Vprofi_HYB                                                           |
C |  Computes vertical profiles, by horizontal averaging                 |
C |  over sig/p hybrid coordinate from RLS file                          |
C +----------------------------------------------------------------------+
C
      SUBROUTINE Vprofi_HYB (idRLS,
     $             idMAR, itM, MARisol, VarName, pVP, varVP)

C ds  IMPLICIT NONE
 
C +...* Large Scale input data dimension.
      include 'LSMARIN.inc'
C +...* MAR dimensions :
      include 'MARdim.inc'

C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM, idRLS
      CHARACTER*(*) VarName
      REAL MARisol(mx,my)

C +---OUTPUT
C +   ~~~~~~
      REAL pVP(LSnk), varVP(LSnk)

C +---LOCAL
C +   ~~~~~
      CHARACTER *10 var_units, tmp_units
      REAL CSTp (LSnk), SIGp (LSnk), LSlev(LSnk)
      REAL wk3MAR(mx,my,mz), p2dMAR(mx,my), spMAR
      REAL valMAR(mz), pMAR(mz)
      REAL tvHYB(LSnk), tmpVP
      REAL MARx(mx), MARy(my), empty1(1), sigma(mz)
      REAL pHYB, vintHYB, ptop
      INTEGER ilv, ii, jj, nelem

C +---Read the RLS coordinate
C +   -----------------------

      CALL UNread
     &   (idRLS,'CSTp', 0, 0, 0, 0,
     &    LSnk    ,1      , 1     ,
     &    LSlev   ,empty1 , empty1,
     &    tmp_units, CSTp)
C +...constant pressure contrib. to coordinate, in Pa

      CALL UNread
     &   (idRLS,'SIGp', 0, 0, 0, 0,
     &    LSnk    ,1      , 1     ,
     &    LSlev   ,empty1 , empty1,
     &    tmp_units, SIGp)

C +---Read/average the MAR coordinate
C     Not really the MAR coord, but MAR SLP
C     and the sigma/pres hyb coordinate of RLS
C +   -----------------------------------------

      CALL UNread
     &   (idMAR, 'pstar', itM, 0, 1,1,mx,my,1,
     &    MARx, MARy, empty1, var_units, p2dMAR)

      CALL VPmean (MARisol, p2dMAR, tmpVP)
C +...surface pressure in kPa

      tmpVP = tmpVP * 1000. ! From kPa to Pa.
      DO ilv = 1, LSnk
        pVP(ilv) = pVP(ilv)
     &           + 0.01 * (CSTp(ilv)+tmpVP*SIGp(ilv))
C       ...MAR pressure on the RLS-type levels, hPa
      ENDDO


C +---Read/average the data
C +   ---------------------
        CALL UNread
     &     (idMAR, VarName, itM, 0,  1,1,mx,my ,mz,
     &      MARx, MARy, sigma, var_units, wk3MAR)
C +..   NB: FileID VarName time 0=3D subregion #levels

      DO ilv = 1, LSnk
        tvHYB(ilv) = 0
      END DO
        
      nelem= 0
      DO ii= 10,mx-10
      DO jj= 10,my-10
C      IF (MARisol(ii,jj).GE.3.5) THEN

         ptop = 10.
C +      ...MAR top pressure, Pa

         spMAR = p2DMAR(ii,jj) * 1000.
C +      ...tot. pressure thickness in MAR, Pa

         DO ilv = 1, mz  
           pMAR(ilv)  = (spMAR*sigma(ilv)+ptop)
C +        ...?log MAR pressure on the MAR levels, Pa
           valMAR(ilv)= wk3MAR(ii,jj,ilv)
C +        ...value of "VarName" on the MAR levels

c          valMAR(ilv)= valMAR(ilv)
c    .                * (100000./pMAR(ilv))**(.28586)
c      TEMPORARY: TEMP to POT TEMP

         END DO

         DO ilv = 1, LSnk
           pHYB= (CSTp(ilv)+spMAR*SIGp(ilv))
C +        ...?log MAR pres at RLS-type hybrid crd., Pa

           CALL INTlin(pMAR, valMAR, mz,
     .                 pHYB, vintHYB)

c          vintHYB = vintHYB 
c    .                * (pHYB/100000.)**(.28586)
c      TEMPORARY: POT TEMP to TEMP

           tvHYB(ilv) = tvHYB(ilv) + vintHYB

         END DO
         nelem= nelem + 1

C      ENDIF
      ENDDO
      ENDDO

      DO ilv = 1, LSnk
        varVP(ilv) = varVP(ilv) + tvHYB(ilv) / nelem
      ENDDO

      END
