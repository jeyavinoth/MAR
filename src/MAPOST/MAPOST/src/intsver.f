C +---------------------------------------------------------------------+
C | MAR (output)                                           03-1999  MAR |
C |   SubRoutine INTsver handles linear interpol. along the vertical    |
C |                  (e.g. to pressure levels)                          |
C |           It is a "simplified" version of INTver: only 1 out. level |
C +---------------------------------------------------------------------+
C |                                                                     |
C |       xin(mx,my,mz)    = input vertical coordinate (ascending order)|
C |       vin(mx,my,mz)    = input value                                |
C |       alvout           = output coordinate (e.g. pressure level)    |
C |       valout(mx,my)    = output values                              |
C |                                                                     |
C |   NOTES:                                                            |
C |      -The 'MISSING VALUE' = 1.0E30 is returned when the requested   |
C |                             alvout is out of xin range              |
C |      -The routine is best for sequential processing, not vectorial. |
C +---------------------------------------------------------------------+
      SUBROUTINE INTsver (xin,vin,ni,nj,nk,alvout,valout)

      IMPLICIT NONE

      INTEGER ni,nj,nk
      REAL xin(ni,nj,nk),vin(ni,nj,nk)
      REAL alvout,valout(ni,nj)
      INTEGER ind, KLO, KHI
      REAL fdis
      INTEGER ii,jj,km

       DO jj= 1, nj
       DO ii= 1, ni

C +---  Search for the appropriate level in the input values:
        KLO=1
        KHI=nk
 1      IF (KHI-KLO.GT.1) THEN
          km=(KHI+KLO)/2
          IF(xin(ii,jj,km).GT.alvout)THEN
            KHI=km
          ELSE
            KLO=km
          ENDIF
        GOTO 1
        ENDIF
        ind=KLO

        IF (alvout.LE.xin(ii,jj,nk) 
     .    .AND. alvout.GE.xin(ii,jj,1)) THEN

C +---    Linearly interpolate:
          fdis  = xin(ii,jj,ind+1)-xin(ii,jj,ind)
          valout(ii,jj)
     .        = vin(ii,jj,ind)* ((xin(ii,jj,ind+1)-alvout) /fdis)
     .      + vin(ii,jj,ind+1)* ((alvout-xin(ii,jj,ind  )) /fdis)

        ELSE
C +---    Set a missing value:
          valout(ii,jj)= 1.0E30

        ENDIF

       ENDDO
       ENDDO

      RETURN
      END
