C +----------------------------------------------------------------------+
C | MAR post-processing                                          01/1999 |
C |                                                                      |
C | Vprofi_MAR                                                           |
C |  Computes vertical profiles, by horizontal averaging                 |
C |  over model (sigma...) levels C                                      |
C +----------------------------------------------------------------------+

      SUBROUTINE Vprofi_MAR (idMAR, itM, MARisol, VarName, pVP, varVP)


C +...* Large Scale input data dimension.
      include 'LSMARIN.inc'
C +...* MAR dimensions :
      include 'MARdim.inc'

C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM
      CHARACTER*(*) VarName

C +---OUTPUT
C +   ~~~~~~
      REAL pVP(mz), varVP(mz)

C +---LOCAL
C +   ~~~~~
      CHARACTER *10 var_units, tmp_units
      REAL sigma(mz)
      REAL wkMAR(mx,my), tmpVP
      REAL MARx(mx), MARy(my), MARz(mz), empty1(1)
      INTEGER ilv

C +---Read/average the data
C +   ---------------------
      DO ilv = 1, mz

        CALL UNread
     &     (idMAR, VarName, itM, ilv, 1,1,mx,my,1,
     &      MARx, MARy, empty1, var_units, wkMAR)
C +..   NB: FileID VarName time level subregion #levels

        CALL VPmean (MARisol, wkMAR, tmpVP)
       
        varVP(ilv) = varVP(ilv) + tmpVP

      ENDDO

C +---Read/average the coordinate
C +   ---------------------------

      CALL UNread
     &   (idMAR, 'pstar', itM, 0, 1,1,mx,my,1,
     &    MARx, MARy, empty1, var_units, wkMAR)

      CALL VPmean (MARisol, wkMAR, tmpVP)
C +...surface pressure in kPa

      CALL UNread
     &   (idMAR,'level', 0, 0, 0, 0,
     &    mz     ,1      , 1     ,
     &    MARz   ,empty1 , empty1,
     &    tmp_units, sigma)

      DO ilv = 1, mz
        pVP(ilv) = pVP(ilv) + 10. * tmpVP * sigma(ilv)
      ENDDO

      END

C +----------------------------------------------------------------------+
C +
      SUBROUTINE VPmean (MARisol, wkMAR, tmpVP)
C +...* MAR dimensions :
      include 'MARdim.inc'
      REAL MARisol(mx,my)     
      REAL wkMAR(mx,my)
      REAL tmpVP
      INTEGER ii, jj, nelem
      
      tmpVP= 0
      nelem= 0
      DO ii= 10,mx-10
      DO jj= 10,my-10
C         IF (MARisol(ii,jj).GE.3.5) THEN
            tmpVP = tmpVP + wkMAR(ii,jj)     
            nelem = nelem + 1
C         ENDIF              
      ENDDO
      ENDDO
      tmpVP= tmpVP / nelem

      END
