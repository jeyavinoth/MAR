C +------------------------------------------------------------------------+
C | MAPOST XF                                                              |
C |   SubRoutine CreateOutNC                                               |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT:                                                               |
C |   ^^^^^^                                                               |
C |                                                                        |
C |   OUTPUT: NetCDF File                                                  |
C |   ^^^^^^                                                               |
C |                                                                        |
C |   CAUTION: 1) This Routine requires the usual NetCDF library,          |
C |   ^^^^^^^^    and a complementary access library called 'libUN.a'      |
C |                                                                        |
C +------------------------------------------------------------------------+
      subroutine CreateOutNC (LoutHOR,LoutBND,LoutVER,LoutMIS,LoutSRF,
     &                        OUTfile, iyrBEG,mmaBEG,jdaBEG,jhuBEG, 
     &                        MARlon, MARlat, MARsig,
     &                        MARx, MARy, RLSlon, RLSlat,
     &                        jhSTP, itlast, ID__nc)

C +
C +
C +--General Variables
C +  =================
C +
C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

C +---MAPOST parameters & variables
C +   -----------------------------
      include 'MAPOST.inc'
      include 'globals.inc'
C +
C +
C +--Routine arguments.
C +  ==================
C +
      LOGICAL LoutHOR,LoutBND,LoutVER,LoutMIS,LoutSRF
      CHARACTER*60 OUTfile
      INTEGER iyrBEG,mmaBEG
      REAL MARlon(mx,my), MARlat(mx,my), MARx(mx), MARy(my)
      REAL MARsig(mz)
      REAL RLSlon(LSni) , RLSlat(LSnj)
      INTEGER ID__nc, jhSTP, itlast

C +--Local   Variables
C +  =================
      common/IDLloc/ fnamNC
C +...               fnamNC: To retain file name.
C +
      PARAMETER (Lfnam= 80, Ltit= 90, Luni= 16+15, Lnam= 13, Llnam=50)
C +...Length of char strings 
C +
      PARAMETER (NdimNC = 9)
C +...Number of defined spatial dimensions (exact)
C +
      PARAMETER (MXdim = 800)
C +...Maximum Number of all dims: recorded Time Steps
C +   and also maximum of spatial grid points for each direction. 
C +
      PARAMETER (MX_var = 120) 
C +...Maximum Number of Variables 
C +
      PARAMETER (NattNC = 2)
C +...Number of real attributes given to all variables
C +
      INTEGER iyr1, iyr2, iyr3, iyr4, mma1, mma2, numtmp
      CHARACTER*1       numstr(0:9)
      REAL              timeNC(MXdim)
      INTEGER           nDFdim(      0:NdimNC)
      REAL              VALdim(MXdim,0:NdimNC)
      INTEGER           NvatNC(NattNC)
      CHARACTER*(Lnam)  NAMdim(      0:NdimNC)
      CHARACTER*(Luni)  UNIdim(      0:NdimNC)
      CHARACTER*(Lnam)  SdimNC(4,MX_var)       
      CHARACTER*(Luni)  unitNC(MX_var)
      CHARACTER*(Lnam)  nameNC(MX_var)
      CHARACTER*(Llnam) lnamNC(MX_var)
      CHARACTER*(Lfnam) fnamNC
      CHARACTER*(Ltit ) tit_NC
      CHARACTER*(Lnam)  NAMrat(NattNC)
      CHARACTER*6   CHRdate
      CHARACTER*120 tmpINP
      INTEGER VARSIZE, ichrsz
      EXTERNAL VARSIZE
cXF
      REAL njmoGE(12),njmbGE(12),njyrGE(12),njybGE(12)
      data  njmoGE
     .     /31,28,31,30, 31, 30, 31, 31, 30, 31, 30, 31/
      data  njmbGE
     .     /0, 1, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0/
      data  njyrGE
     .     /0,31,59,90,120,151,181,212,243,273,304,334/
      data  njybGE
     .     /0, 0, 1, 1,  1,  1,  1,  1,  1,  1,  1,  1/


      DATA numstr/'0','1','2','3','4','5','6','7','8','9'/
C +
      CALL DInfo(2,'OUTmpp: Begin')

C +
C +--NetCDF File Initialization
C +  ============================
C +
C +--Output file
C +  -----------
C +
      iyr1  =  iyrBEG/1000
      numtmp=  iyrBEG-iyr1*1000
      iyr2  =  numtmp/100             
      numtmp=  numtmp-iyr2*100
      iyr3  =  numtmp/10           
      iyr4  =  numtmp-iyr3*10
      mma1  =  mmaBEG/10
      mma2  =  mmaBEG-mma1*10

      CHRdate= numstr(iyr1) // numstr(iyr2)
     .      // numstr(iyr3) // numstr(iyr4)
     .      // numstr(mma1) // numstr(mma2)

      ichrsz = VARSIZE (OUTfile)
C +.. (from libUN)
c     fnamNC = OUTfile(1:ichrsz)//'_'//CHRdate//'.nc'
cXF
      fnamNC = 'output/MPP.'
     .         //CHRdate//'.'
     .         //OUTfile(1:ichrsz)//'.nc'
C +
C +
C +--Output Title
C +  ------------
C +
c      tit_NC = 'MAR post processor output'
c    .        // ' - Exp: ' // CHRdate
cXF
       tit_NC = 'MAPOST'
     .        // ' - Exp: ' // OUTfile(1:ichrsz)
     .        // ' - '      // CHRdate
C +
C +
C +--Create File / Write Constants
C +  -----------------------------
       MMXstp = MXdim
C +...To check array bounds... silently
C +
C +--Time Variable (hour)
C +  ~~~~~~~~~~~~~~~~~~~~
C +...  Define a NetCDF dimension (size, name, unit):
        nDFdim(0)= itlast+1
        NAMdim(0)= 'time'
cXF     UNIdim(0)= 'days'
        UNIdim(0)= 'HOURS since 1901-01-15 00:00:00'
C +
C +...  Check temporary arrays: large enough ?
        IF ((itlast+1).gt.MMXstp)
     &  STOP '** OUTmpp - ERROR : MXdim too low (1)**'
C +
C +
        DO it = 1,itlast+1
cXF        timeNC(it)  = FLOAT((it-1)*jhSTP) 
cXF        VALdim(it,0)= timeNC(it) 
           VALdim(it,0)= real
     .                  ((351  +(iyrBEG  -1902) *365   ! Nb Days before iyrrGE
     .                  +(iyrBEG  -1901) /  4          ! Nb Leap Years
     .                  + njyrGE(mmaBEG)               ! Nb Days before mmarGE
     .                  + njybGE(mmaBEG)               ! (including Leap Day)
     .                  * max(0,1-mod(iyrBEG,4))       !
     .                  + jdaBEG     -1      )*  24    !
     .                  + jhuBEG )                     !
     .                  + (it-1)*jhSTP                 ! hours
           timeNC(it)  = VALdim(it,0)
        ENDDO
C +
C     Define horizontal spatial dimensions :    
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +...Check temporary arrays: large enough ?
      IF (    mx  .gt.MMXstp.or.my  .gt.MMXstp
     &    .or.LSni.gt.MMXstp.or.LSnj.gt.MMXstp
     &    .or.mzz .gt.MMXstp.or.LSnk.gt.MMXstp)
     &  STOP '** OUTgra - ERROR : MXdim too low (2)**'
C +
C +...Define NetCDF dimensions (size, name, unit):
      CALL DInfo(2,'OUTmpp: define spatial dimensions...')
C +
C +- -Dimensions of MAR:
C +   - - - - - - - - - -
      DO i = 1, mx
        VALdim(i,1) = MARx(i)
      END DO
      nDFdim(1)= mx
      NAMdim(1)= 'x'
      UNIdim(1)= 'km'
C +
      DO j = 1, my
        VALdim(j,2) = MARy(j)
      END DO
      nDFdim(2)= my
      NAMdim(2)= 'y'
      UNIdim(2)= 'km'
C +
      do k = 1, mz
        VALdim(k,3) = MARsig(k)
      enddo
      nDFdim(3)= mz
      NAMdim(3)= 'MARlev'
      UNIdim(3)= '[sigma]'
C +
C +- -Dimensions of Large Scale data (RLS)
C +   - - - - - - - - - - - - - - - - - - -
      CALL DInfo(2,'OUTmpp: RLS dimensions...')

      DO i = 1, LSni
        VALdim(i,4) = RLSlon(i)
      END DO
      nDFdim(4)= LSni
      NAMdim(4)= 'LS_lon'
      UNIdim(4)= 'degrees'
C +
      DO j = 1, LSnj
        VALdim(j,5) = RLSlat(j)
      END DO
      nDFdim(5)= LSnj
      NAMdim(5)= 'LS_lat'
      UNIdim(5)= 'degrees'

      do k = 1, LSnk
        VALdim(k,6) = k
      enddo
      nDFdim(6)= LSnk
      NAMdim(6)= 'RLSlev'
      UNIdim(6)= '-'

C +- -Other dimensions:
C +   - - - - - - - - -
      do k = 1, nreg
        VALdim(k,7) = k
      enddo
      nDFdim(7)= nreg
      NAMdim(7)= 'SubReg'
      UNIdim(7)= '-'
C +... For Averaging Sub-regions

      do k = 1, ndb 
        VALdim(k,8) = k
      enddo
      nDFdim(8)= ndb
      NAMdim(8)= 'dbound'
      UNIdim(8)= '-'
C +... For fn of "dbound", dist. to boundary.

      do k = 1, npl
        VALdim(k,9) = plevel(k) * 10. !kPa -> hPa
      enddo
      nDFdim(9)= npl
      NAMdim(9)= 'plevel'
      UNIdim(9)= 'hPa'
C +... For MAPOST pressure levels.

C +
C +--Variable's Choice (Table MPPvou.dat)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL DInfo(2,'OUTmpp: Variables table...')
C +
        OPEN(unit=15,status='unknown',file='MPPvou.dat')
C +
        itotNC = 0
 980    CONTINUE
          READ (15,'(A120)',end=990) tmpINP
          IF ( tmpINP(1:4).eq.'    '  
     &    .OR.(tmpINP(1:4).eq.'H   '.AND.LoutHOR)
     &    .OR.(tmpINP(1:4).eq.'V   '.AND.LoutVER)
     &    .OR.(tmpINP(1:4).eq.'B   '.AND.LoutBND)
     &    .OR.(tmpINP(1:4).eq.'M   '.AND.LoutMIS)
     &    .OR.(tmpINP(1:4).eq.'S   '.AND.LoutSRF) ) THEN

            itotNC = itotNC + 1
            IF (itotNC.GE.MX_var) THEN
              write(*,*) 'OUTmpp: sorry, MX_var too small'
              STOP
            ENDIF
            READ (tmpINP,'(4x,5A9,A12,A50)')
     &          nameNC(itotNC)  ,SdimNC(1,itotNC),SdimNC(2,itotNC),
     &          SdimNC(3,itotNC),SdimNC(4,itotNC),
     &          unitNC(itotNC)  ,lnamNC(itotNC)
C +...          nameNC: Name
C +             SdimNC: Names of Selected Dimensions (max.4/variable) 
C +             unitNC: Units
C +             lnamNC: Long_name, a description of the variable
C +
          ENDIF
        GOTO 980
 990    CONTINUE
C +
        NtotNC = itotNC 
C +...  NtotNC : Total number of variables writen in NetCDF file.
C +
C +--List of NetCDF attributes given to all variables:
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +... The "actual_range" is the (min,max)
C +    of all data for each variable:
      NAMrat(1) = 'actual_range'
      NvatNC(1) = 2

C +... The "[var]_range" is NOT of attribute type,
C +    it is a true variable containing the (min,max) for
C +    each level, for 4D (space+time) variables only
C +    (automatically handled by UN library;
C +     must be the LAST attribute)
      NAMrat(NattNC) = '[var]_range'
      NvatNC(NattNC) = 2
C +
C +--Automatic Generation of the NetCDF File Structure
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL DInfo(2,'Calling UNscreate')
C +
C +     **************
        CALL UNscreate (fnamNC,tit_NC,
     &                  NdimNC, nDFdim, MXdim , NAMdim, UNIdim, VALdim,
     &                  MX_var, NtotNC, nameNC, SdimNC, unitNC, lnamNC,
     &                  NattNC, NAMrat, NvatNC,
     &                  ID__nc) 
C +     **************
C +
C +
C +--Write Time - Constants
C +  ~~~~~~~~~~~~~~~~~~~~~~
C +     ************
        CALL UNwrite (ID__nc, 'lat   ', 1  , mx    , my, 1 , MARlat)
        CALL UNwrite (ID__nc, 'lon   ', 1  , mx    , my, 1 , MARlon)
C +     ************
C +
C +
      return
      end
