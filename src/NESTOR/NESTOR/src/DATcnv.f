C   +-------------------------------------------------------------------+
C   |  Subroutine DATcnv                              May 2001  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : Date : year, month, day, hour                             |
C   | ^^^^^^^ or I_date (integer)                                       |
C   |         datTOint : if true, transformation from date to I_date    |
C   |                    if false,transformation from I_date to date    |
C   |                                                                   |
C   | Output: I_date (integer)                                          |
C   | ^^^^^^^ or Date (year, month, day, hour)                          |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE DATcnv (year,month,day,hour,I_date,datTOint)


      IMPLICIT NONE

      INCLUDE 'LSCmod.inc'

      INTEGER imonth,year,month,day,hour,i,Nmonth(12),Nhour(0:12),
     .        nday,nextnh,iyear

      INTEGER*4 I_date,N_date,nyhour
      INTEGER   y365, y366

      LOGICAL datTOint


C +---Number of days in each month
C +   ----------------------------

      DATA (Nmonth(i),i=1,12)
     .     /31,   ! January
     .      28,   ! February
     .      31,   ! March
     .      30,   ! April
     .      31,   ! May
     .      30,   ! June
     .      31,   ! July
     .      31,   ! August
     .      30,   ! September
     .      31,   ! October
     .      30,   ! November
     .      31/   ! December

      IF (M30d) THEN ! For LMD standard models: 30days/month
         DO i=1,12
            Nmonth(i)=30
         ENDDO
                  y365 = 360
                  y366 = 360
      ELSE
                  y365 = 365
                  y366 = 366
         if(f28d) y366 = 365
      ENDIF


C +   ******************
      IF (datTOint) THEN
C +   ******************


C +---Case of bisextil year
C +   ---------------------
      IF (M30d.eqv..false..and.f28d.eqv..false.) THEN
        IF (mod(year,4).eq.0.and.(mod(year,100).ne.0
     .                        .or.mod(year,400).eq.0)) THEN
         Nmonth(2)=29
        ELSE
         Nmonth(2)=28
        ENDIF
      ENDIF

C +---Convertion in hours
C +   -------------------

      nday    =0
      Nhour(0)=0

      DO i=1,12
       nday    =nday+Nmonth(i)
       Nhour(i)=nday*24
      ENDDO


C +---Hours from year 0 to the considered year
C +   ----------------------------------------

      nyhour = 0

      DO iyear=0,year-1
        IF (mod(iyear,4).eq.0.and.(mod(iyear,100).ne.0
     .                         .or.mod(iyear,400).eq.0)) THEN
        nyhour = nyhour + y366*24
       ELSE
        nyhour = nyhour + y365*24
       ENDIF
      ENDDO

C +---Convert from DATE to I_date
C +   ---------------------------

       I_date=nyhour+Nhour(month-1)+(day-1)*24+hour


C +   ****
      ELSE
C +   ****


C +---Search for year
C +   ---------------

       nyhour = I_date

       iyear  = 0
       nextnh = y366*24

       DO WHILE (nyhour.ge.nextnh)
        nyhour = nyhour - nextnh
        iyear  = iyear  + 1
        IF (mod(iyear,4).eq.0.and.(mod(iyear,100).ne.0
     .                         .or.mod(iyear,400).eq.0)) THEN
         nextnh = y366*24
        ELSE
         nextnh = y365*24
        ENDIF
       ENDDO

       year   = iyear
       N_date = nyhour


C +---Case of bisextil year
C +   ---------------------
      IF (M30d.eqv..false..and.f28d.eqv..false.) THEN
        IF (mod(year,4).eq.0.and.(mod(year,100).ne.0
     .                        .or.mod(year,400).eq.0)) THEN

         Nmonth(2)=29
        ELSE
         Nmonth(2)=28
        ENDIF
      ENDIF



C +---Convertion in hours
C +   -------------------

      nday    =0
      Nhour(0)=0

      DO i=1,12
       nday    =nday+Nmonth(i)
       Nhour(i)=nday*24
      ENDDO


C +---Convert from I_date to DATE
C +   ---------------------------

       imonth=0

       DO i=1,12
        IF ((N_date.ge.Nhour(i-1)).and.(N_date.lt.Nhour(i))) imonth=i
       ENDDO

       print *,N_date,imonth

       IF (imonth.eq.0) THEN
        write(6,*) 'I_date =',N_date,' cannot be converted.'
        write(6,*) 'STOP.'
        STOP
       ENDIF
 
       month=imonth
       day  =AINT(REAL(N_date-Nhour(month-1))/24.) + 1
       hour =N_date-Nhour(month-1)-24*(day-1)


C +   *****
      ENDIF
C +   *****


      RETURN
      END
