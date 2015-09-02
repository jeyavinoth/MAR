      SUBROUTINE DInfo (ilv, TitStr) 
      COMMON /CoInfo/ icheck  
      INTEGER icheck, ilv
      CHARACTER*(*) TitStr
          IF (icheck.GE.ilv) THEN 
             WRITE(6,*) TitStr
          ENDIF
      END

      SUBROUTINE VInfo (ilv, TitStr,Val) 
      COMMON /CoInfo/ icheck 
      INTEGER icheck, ilv
      CHARACTER*(*) TitStr
      REAL Val
          IF (icheck.GE.ilv) THEN
             WRITE(6,*) TitStr, Val
          ENDIF
      END

C     ****************************************************
C     Search:
C     Returns the KLO such that the "vreq" values is in
C     interval [xin(KLO),xin(KLO+1)] 
C           or [xin(KLO+1),xin(KLO)] :
C     in the later case, it means xin(i+1) < xin(i)

      SUBROUTINE SEARCH(xin, nxin, vreq, KLO)
 
      IMPLICIT NONE
      INTEGER nxin
      REAL xin(nxin), vreq
      INTEGER KLO, KHI, km

      IF (xin(nxin).GT.xin(1)) THEN

        KLO=1
        KHI=nxin
 1      IF (KHI-KLO.GT.1) THEN
          km=(KHI+KLO)/2
          IF(xin(km).GT.vreq)THEN
            KHI=km
          ELSE
            KLO=km
          ENDIF
        GOTO 1
        ENDIF

      ELSE

        KLO=1
        KHI=nxin
 2      IF (KHI-KLO.GT.1) THEN
          km=(KHI+KLO)/2
          IF(xin(km).LT.vreq)THEN
            KHI=km
          ELSE
            KLO=km
          ENDIF
        GOTO 2
        ENDIF

      ENDIF

      RETURN
      END
