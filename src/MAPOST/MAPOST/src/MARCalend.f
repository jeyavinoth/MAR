C +************************************************************************+
C | MARCalend XF                                                           |
C |                                                                        |
C |                                                                        |
C +************************************************************************+
C |                                                                        |
C | Version    January 1999                                                |
C | ^^^^^^^^^^^^^^^^^^^^^^^                                                |
C |                                                                        |
C | ---Ph. Marbaix                                                         |
C +************************************************************************+
      SUBROUTINE MARCalend (iyrBEG,mmaBEG,jdaBEG,jhuBEG, iadtime,
     .                      iyrCUR,mmaCUR,jdaCUR,jhuCUR)


      INCLUDE 'globals.inc'

C +---Arguments
C +   ~~~~~~~~~
      INTEGER iyrBEG,mmaBEG,jdaBEG,jhuBEG, iadtime
     .        jhSTP, iyrCUR,mmaCUR,jdaCUR,jhuCUR,ii

C +---Local
C +   ~~~~~
      INTEGER mmaTMP,jdaTMP,jhuTMP

C +---Time Constants
C +   ~~~~~~~~~~~~~~
      INTEGER njmoGE(0:12), njyrGE(0:13)
      data (njmoGE(n),n=0,12)
     .     /0,31,28,31,30, 31, 30, 31, 31, 30, 31, 30, 31/
      data (njyrGE(n),n=0,13)
     .     /0, 0,31,59,90,120,151,181,212,243,273,304,334, 365/
C +...     njmoGE: Nb of Days in each Month of the Year
C +---     njyrGE: Nb of Days since   Begin of the Year,
C +---                        before  Current Month

      IF (LSCmod.EQ.'LMz') THEN
        DO ii=1,12
         njmoGE(ii)=30
        ENDDO
        DO ii=2,13
         njyrGE(ii)=30*(ii-1)
        ENDDO
      ENDIF 

C +---

      IF (mod(iyrBEG,4).EQ.0) THEN
       njmoGE(2)=29
       do ii=3,13
       njyrGE(ii) = njyrGE(ii) + 1
       enddo
      ENDIF

C +---Hour:
      jhuTMP = jhuBEG + iadtime
      jhuCUR = mod (jhuTMP, 24)  

C +---Days from begin of year:
      jdaTMP = njyrGE(mmaBEG) + jdaBEG + (jhuTMP-jhuCUR)/24

C +---Month:
      mmaCUR = mmaBEG
  10  CONTINUE
      IF (njyrGE(mmaCUR+1).LT.jdaTMP) THEN
         mmaCUR = mmaCUR + 1
         GOTO 10
      ENDIF

C +---Day:
      jdaCUR = jdaTMP - njyrGE(mmaCUR)

C +---Year:
      iyrCUR = iyrBEG
C     ...Year change not handled for the moment.
      
      END
