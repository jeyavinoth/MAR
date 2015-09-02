C +-------------------------------------------------------------------+
C | TIME CONVERSION in HOURS                     E. Vanvyve - 2003.04 |
C |                                                                   |
C |   Subroutine CtimeDD computes the nbr of hours since 15.01.1901   |
C |                      0:00 (cfr MAR output) to the current date    |
C +-------------------------------------------------------------------+
C |                                                                   |
C |  INPUT : Cdate  : current date  (format : mmddhh)                 |
C |          iyrCUR : current year  (format : aaaa  )                 |
C |          mmaCUR : current month (format : mm    )                 |
C |          jdaCUR : current day   (format : dd    )                 |
C |          jhuCUR : current hour  (format : hh    )                 |
C |                                                                   |
C |  OUTPUT: Ctime  : current time  (format : nbr of hours since      |
C |                                  15.01.1901 0:00 (cfr MAR output))|
C |                                                                   |
C |  NOTES: changes in MAPOST.F (see end of this file)                |
C |                                                                   |
C +-------------------------------------------------------------------+

      SUBROUTINE CtimeDD (Cdate, iyrCUR, mmaCUR, jdaCUR, jhuCUR, Ctime)
      
      IMPLICIT NONE

C +---Arguments
C +   ---------
      INTEGER, INTENT(IN)  :: iyrCUR, mmaCUR, jdaCUR, jhuCUR
      REAL, INTENT(IN)  :: Cdate
      REAL, INTENT(OUT) :: Ctime

C +---Local variables
C +   ---------------
      INTEGER :: y, n, L004, L100, L400
      INTEGER :: njyrNY(0:13), njyrLY(0:13)
        !njyrN/LY: nb of days since begin of the (Normal/Leap) year,
        !          before current month
      data (njyrLY(n),n=0,13) /0, 0,31,60,91,121,152,182,213,244,274,
     .                            305,335,366/            !leap   year
      data (njyrNY(n),n=0,13) /0, 0,31,59,90,120,151,181,212,243,273,
     .                            304,334,365/            !normal year
     

C +---Existence conditions
C +   --------------------

      n = (mmaCUR*100 + jdaCUR) *100 + jhuCUR
      IF (n .ne. Cdate) THEN
        write(*,*) 'Error in Cdate (n/Cdate) ', n, Cdate
        STOP
      ENDIF

      IF (iyrCUR .le. 1901) THEN
        write(*,*) 'TIMyear.f error : iyrCUR ', iyrCUR,' <= 1901'
        STOP
      ENDIF

C +---Number of days upto iyrCUR-1
C +   ----------------------------

      !15.01.1901 > 31.12.1901
      Ctime=njyrNY(13)-15
      !01.01.1902 > 31.12.iyrCUR-1
      Ctime=Ctime+(iyrCUR-1902)*365  !(-1-1901=-1902)
      !01.01.1902 > 31.12.iyrCUR-1 (number of 29.02)
      n=0
      DO y=1902,iyrCUR-1
        L004 = mod (y, 4  )
        L100 = mod (y, 100)
        L400 = mod (y, 400)
        IF (L004==0 .and. .not.(L100==0 .and. L400/=0)) THEN
          n=n+1
        ENDIF
      ENDDO
      Ctime=Ctime+n
      !01.01.iyrCUR to jdaCUR.mmaCUR.iyrCUR 0:00
      L004 = mod (iyrCUR, 4  )
      L100 = mod (iyrCUR, 100)
      L400 = mod (iyrCUR, 400)
      IF (L004==0 .and. .not.(L100==0 .and. L400/=0)) THEN
        Ctime=Ctime+njyrLY(mmaCUR)+jdaCUR
      ELSE
        Ctime=Ctime+njyrNY(mmaCUR)+jdaCUR
      ENDIF
      !in hours + jhuCUR
      Ctime=Ctime*24+jhuCUR
      
      END

C +-------------------------------------------------------------------+
C | CHANGES TO DO in MAPOST.F :                                       |
C |                                                                   |
C | in  Compute the current date                                      |
C |     ========================                                      |
C |      CALL MARCalend (iyrBEG,mmaBEG,jdaBEG,jhuBEG, iadtime,        |
C |     .                iyrCUR,mmaCUR,jdaCUR,jhuCUR)                 |
C |      Cdate = (mmaCUR*100 + jdaCUR) *100 + jhuCUR                  |
C |   >  CALL CtimeDD (Cdate, iyrCUR, mmaCUR, jdaCUR, jhuCUR, Ctime)  |
C |                                                                   |
C | in  Open Files / MAR files                                        |
C |     ==========                                                    |
C |   <  CALL UNgindx (idMAR, "date", Cdate, Rdate, Fdate, itM)       |
C |   >  CALL UNgindx (idMAR, "time", Ctime, Rtime, Ftime, itM)       |
C |                                                                   |
C |   <  IF (ABS(Rdate-Cdate).GT.0.5) THEN                            |
C |   >  IF (ABS(Rtime-Ctime).GT.0.5) THEN                            |
C |        nMAR = nMAR + 1                                            |
C |        IF (nMAR .GT. ntMAR) THEN                                  |
C |    <      write(*,*) 'MAPOST - error: Requested date is not'      |
C |    <      write(*,*) '  available in MAR input, date=', Cdate     |
C |    >      write(*,*) 'MAPOST - error: Requested time is not'      |
C |    >      write(*,*) '  available in MAR input, time=', Ctime     |
C |                                                                   |
C +-------------------------------------------------------------------+
