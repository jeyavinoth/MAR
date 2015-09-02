C   +-------------------------------------------------------------------+
C   |  Subroutine LSCinp                            April 2001  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : - DATtim : date given in hours from beginning of the year |
C   | ^^^^^^^ - LSCmod : LSC model used for init. and forcing fields    |
C   |                                                                   |
C   | Output: - LSCfil : file to be read for the fields at DATtim       |
C   | ^^^^^^^ - I_time : time corresponding to DATtim in LSCfil         |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE LSCinp


      IMPLICIT NONE


C +---General variables
C +   -----------------

      INCLUDE 'NetCDF.inc'
      INCLUDE 'NSTdim.inc'
      INCLUDE 'LSCvar.inc'
      INCLUDE 'NESTOR.inc'      

C +---Local variables
C +   ---------------

      INTEGER nfile,maxfile
      PARAMETER (maxfile=200)

      INTEGER i,itR,Ierro,FILEid,
     .        DATiyr,DATmma,DATjda,DATjhu

      REAL    Cdate,Rdate,Fdate,yyyy,mm,dd,hh

      CHARACTER*15  namTIM
      CHARACTER*200 LSCfln(maxfile),LSCtitle,xline

      CHARACTER*10 var_units

      LOGICAL Vtrue,Vfalse


C +---Local data
C +   ----------

      DATA Vtrue  /  .true. /
      DATA Vfalse / .false. /


C +---Variable name for date
C +   ----------------------

      IF (LSCmod.eq.'LMD') THEN
       namTIM='time_counter'
      ELSE
       namTIM='date'
      ENDIF


C +---Requested date
C +   --------------

C +        ******
      CALL DATcnv (DATiyr,DATmma,DATjda,DATjhu,DATtim,Vfalse)
C +        ******

      Cdate = (DATmma*100 + DATjda) *100 + DATjhu

      WRITE(6,*) 'Date'
      WRITE(6,*) '~~~~'
      WRITE(6,1000) DATiyr,DATmma,DATjda,DATjhu
1000  FORMAT(' Year,month,day,hour : ',i5,',',i3,',',i3,',',i3)
      WRITE(6,*)


C +---Input file names
C +   ----------------

      OPEN (unit=52,status='old',file='LSCfil.dat')
      i = 0

210   CONTINUE
      READ (52,'(a100)',END=230) xline
      IF (xline.eq.'') GOTO 210
      IF (xline(1:1).eq.' ') THEN
       write(6,*) 'Blank characters in LSCfil.dat. Please remove'
       write(6,*) 'them and restart NESTOR.'
       write(6,*) 'STOP in LSCinp.f'
       STOP
      ENDIF
      i = i + 1
      LSCfln(i) = xline
      GOTO 210

230   CONTINUE
      CLOSE(unit=52)

      nfile = i

      IF (nfile.gt.maxfile) THEN
       write(6,*) 'Increase maxfile in LSCinp.f '
       write(6,*) 'Error   -   STOP '
       STOP
      ENDIF


C +---Search the LSC file for the requested date
C +   ------------------------------------------

      I_time = -1

      DO i=1,nfile

C +         *******
       CALL UNropen (LSCfln(i),FILEid,LSCtitle)
C +         *******
       CALL UNgindx (FILEid,namTIM,Cdate,Rdate,Fdate,itR)
C +         *******

       IF (ABS(Rdate-Cdate).LE.0.5) THEN
        LSCfil = LSCfln(i)
        I_time = itR
       ENDIF

C +         ******
       CALL NCCLOS (FILEid,Ierro)
C +         ******

      ENDDO


C +---Case of no data file found
C +   --------------------------

      IF (I_time.eq.(-1)) THEN

       write(6,*) 'No LSC data file found for the following date :'
       write(6,*) DATiyr,DATmma,DATjda,DATjhu
       write(6,*)
       write(6,*) '              --- STOP in LSCinp ---           '
       write(6,*)
 
       STOP

      ELSE

       CALL UNropen (LSCfil,FILEid,LSCtitle)

       CALL UNsread (FILEid,'YEAR',I_time,1,
     &               I_time,I_time,1,1,1,var_units,yyyy)

       CALL UNsread (FILEid,'MONTH',I_time,1,
     &               I_time,I_time,1,1,1,var_units,MM)

       CALL UNsread (FILEid,'DAY',I_time,1,
     &               I_time,I_time,1,1,1,var_units,DD)      

       CALL UNsread (FILEid,'HOUR',I_time,1,
     &               I_time,I_time,1,1,1,var_units,HH)

       CALL NCCLOS (FILEid,Ierro)

       write(6,'(a14,i5,3i3)') " LSCfile date:",
     &                         int(yyyy),int(MM),int(DD),int(HH)
 

      ENDIF


      RETURN
      END
