C +----------------------------------------------------------------------+
C | MAR post-processing                                          01/1999 |
C |                                                                      |
C | VProfi_RLS                                                           |
C |  Computes vertical profiles, by horizontal averaging                 |
C |  over model (sigma...) levels C                                      |
C +----------------------------------------------------------------------+

      SUBROUTINE VProfi_RLS (MARlon, MARlat, MARisol,
     &                       idRLS, itR, VarName, pVP, varVP)


C +...* Large Scale input data dimension.
      include 'LSMARIN.inc'
C +...* MAR dimensions :
      include 'MARdim.inc'

C +---INPUT
C +   ~~~~~
      INTEGER idRLS, itR
      CHARACTER*(*) VarName

C +---OUTPUT
C +   ~~~~~~
      REAL pVP(LSnk), varVP(LSnk)

C +---LOCAL
C +   ~~~~~
      CHARACTER *10 var_units, tmp_units
      REAL CSTp(LSnk), SIGp(LSnk), LSlev(LSnk)
      REAL wkRLS(LSni,LSnk), wkMAR(mx,my), tmpVP
      REAL LSlon(LSni), LSlat(LSnj), empty1(1)
      INTEGER ilv

C +---Read/average the data
C +   ---------------------
      DO ilv = 1, LSnk

        CALL UNread
     &     (idRLS, VarName, itR, ilv, 1,1,LSni,LSnj,1,
     &      LSlon, LSlat  , empty1, var_units, wkRLS)
C +..   NB: FileID VarName time level subregion #levels

        CALL INThor (1,
     &      LSlon  , LSlat , wkRLS,
     &      MARlon, MARlat , wkMAR)

        CALL VPmean (MARisol, wkMAR, tmpVP)

        varVP(ilv) = varVP(ilv) + tmpVP

      ENDDO

C +---Read/average the coordinate
C +   ---------------------------

      CALL UNread
     &   (idRLS, 'SP', itR, 0, 1,1,LSni,LSnj,1,
     &    LSlon, LSlat, empty1, tmp_units, wkRLS)

      CALL INThor (1,
     &    LSlon  , LSlat , wkRLS,
     &    MARlon, MARlat , wkMAR)

      CALL VPmean (MARisol, wkMAR, tmpVP)
      
C +...surface pressure in Pa 

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

      DO ilv = 1, LSnk
        pVP(ilv) = pVP(ilv)
     &           + 0.01 * (CSTp(ilv)+tmpVP*SIGp(ilv))
C +...pressure in hPa

      ENDDO

      END
