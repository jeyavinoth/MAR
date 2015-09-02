C +************************************************************************+
C | MAPOST XF                                                              |
C |                                                                        |
C | Post-processor for climate runs with the Modele Atmospherique Regional |
C +************************************************************************+
C |                                                                        |
C | Version & contributors: see below or run this program                  |
C |                                                                        |
C +************************************************************************+

      PROGRAM MAPOST 

C +************************************************************************+

      IMPLICIT NONE

C +---External function (from libUN.a, finds true size of strings)
C +   -----------------
      INTEGER  VARSIZE
      EXTERNAL VARSIZE

C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

C +---Physics (smaller package with defined variables.)
C +   -------------------------------------------------
      include 'LSMphy.inc'

C +---MAPOST parameters & variables
C +   -----------------------------
      include 'MAPOST.inc'
      include 'globals.inc'
      
C +...*For numerical filter in HDynSDBP:
C +   ----------------------------------
      INCLUDE 'HDynSDBP.inc'

C +---Local variables
C +   ---------------

C +---General:
C +   ~~~~~~~~
      LOGICAL LoutHOR, LoutBND, LoutVER ,LoutMIS, LoutSRF, LoutSta
      INTEGER itexpe,iyrBEG,mmaBEG,jdaBEG,jhuBEG, iadtime,
     .        jdDLN,jhDLN, jhSTP, iyrCUR,mmaCUR,jdaCUR,jhuCUR,
     .        intype, istat, listat, idt, ntt, idt1, itlast,
     .        iError, ireg, ilv, idf, MAPtyp
      REAL    Cdate, Fdate, Rdate, empty1(1)
      CHARACTER*100 xline, nfile
      CHARACTER*80 UNver,NCDFver
      CHARACTER*13 tmp_units
      
C +---Sub-region and station definitions:
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER InReg(mx,my,nreg)
      REAL    RegMsk(mx,my)
      REAL    wmSta (mx,my,nsta)  
      REAL    StaMsk(mx,my)

C +---RLS files
C +   ~~~~~~~~~
      CHARACTER*100 RLSfile(100)
      CHARACTER*80  RLStit
      CHARACTER*3  prefix
      INTEGER ntRLS, nRLS, idRLS, itR
      REAL RLSlon(LSni), RLSlat(LSnj), shRLS(LSni,LSnj)

C +---MAR files
C +   ~~~~~~~~~
      CHARACTER*100 MARfile(10)
      CHARACTER*80  MARtit
      INTEGER ntMAR, nMAR, idMAR, itM
      REAL MARlon(mx,my), MARlat(mx,my), MARisol(mx,my)
      REAL shMAR(mx,my), MARsig(mz), tmp_mz(mz)
      REAL MARx(mx), MARy(my), shIRLS(mx,my)

C +---MAPOST output
C +   ~~~~~~~~~~~~~~~
      CHARACTER*60 OUTfile
      INTEGER iOUTnc, iptMKE

C +---Data exchange between MAPOST components
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     This is a temporary solution, i.e.
C     it is desirable to find better.
      REAL slpMAR(mx,my)


C +---Program testing, "control-ouput" amount:
      INTEGER icheck
      COMMON /CoInfo/ icheck
      icheck = 2

      Write(*,*)('              ---------------------             ')
      Write(*,*)('                   MAPOST 3.1                   ')
      Write(*,*)('              Creation:       Jan99             ')
      Write(*,*)('              Last Major Rev: Feb02             ')
      Write(*,*)('              Last Minor Rev: >may              ')
      Write(*,*)('              ---------------------             ')
      Write(*,*)('                                                ')
      Write(*,*)('  mailto:philippe.marbaix@lmd.jussieu.fr        ')
      Write(*,*)('  (this also include contributions from others) ')
      Write(*,*)('                                                ')


C +---Initialisation
C +   ===============
      CALL DInfo(1,'Initialisation...')

C +   Reading control files
C +   --------------------

      OPEN (unit=51,status='old',file='MAPOST.ctr')
      read (51,*) !- - - - - - - - - - - - - - - - - -
      read (51,*) !- - - - - - - - - - - - - - - - - -
      read (51,'(A3)')   LSCmod
      read (51,'(A3)')   REGmod
      read (51,'(A3)')   REGstudy
      read (51,*) !- - - - - - - - - - - - - - - - - -
      read (51,'(A60)')  OUTfile
      read (51,*) !- - - - - - - - - - - - - - - - - -
      read (51,'(i4,3(1x,i2))')iyrBEG,mmaBEG,jdaBEG,jhuBEG         
      read (51,'(i3,  1x,i2 )')              jdDLN ,jhDLN             
      read (51,'(1i3)')                             jhSTP
      read (51,*) !- - - - - - - - - - - - - - - - - -
      read (51,'(1i3)')                       MAPtyp 
      read (51,*) !- - - - - - - - - - - - - - - - - -
      read (51,'(l4)')    LoutHOR
      read (51,'(l4)')    LoutSta
      read (51,'(l4)')    LoutBND
      read (51,'(l4)')    LoutVER
      read (51,'(l4)')    LoutMIS
      read (51,'(l4)')    LoutSRF
      read (51,*) !- - - - - - - - - - - - - - - - - -
      CLOSE(unit=51)

      IF (LoutBND.AND.(.NOT.LoutHOR)) THEN
         write(*,*) 'WARNING: LBC section requires MSLP'
         write(*,*) 'computed in Horiz/Dyn    '
      ENDIF

C +   Read input files names
C +   ----------------------
 
C +---RLS files names
C     ~ ~ ~ ~ ~ ~ ~ ~
      OPEN (unit=52,status='old',file='LSCfil.dat')
      idf = 0
 210  continue
        read (52,'(a100)',end=230) xline
        idf = idf + 1
        RLSfile(idf) = xline
        goto 210

 230  continue
      close(unit=52)
      ntRLS = idf      
      nRLS  = 1
      idRLS = -1

C +---MAR files names
C     ~ ~ ~ ~ ~ ~ ~ ~
      OPEN (unit=52,status='old',file='MARfil.dat')
      idf = 0
 310  continue
        read (52,'(a100)',end=330) xline
        idf = idf + 1
        MARfile(idf) = xline
        goto 310

 330  continue
      close(unit=52)
      ntMAR = idf
      nMAR  = 1
      idMAR = -1

      write(*,*) 'Number of LS  files (LSCfil.dat):', ntRLS
      write(*,*) 'Number of MAR files (MARfil.dat):', ntMAR


      istat= 0
C +   ^Tell stat. routines to intialise variables.

C +    ---------------------------------------------
C +    NOTE: istat tells what to do:
C +    (0 -> init, 
C +     1 -> continue, 
C +     2 -> terminate work)
C +
C +    this is a method to make the main part
C +    of the code more simple: stat. routines must
C +    intialise the statistical variables internaly, 
C +    and divide by the number of elements or 
C +    do extra computations at the last time-step.

C +   Initalise some global libUN parameters:
C +  ----------------------------------------

C +   Check all variables for anomalous values at reading
C +   (we must first reset params, otherwise it will be
C +    done at first libUN call)
      CALL UNparam('RESET_PARAMS_',0.0   )
      CALL UNparam('READOVER_WARN',1.0E10)

      CALL UNversion(UNver,NCDFver)
      write(*,*) 'NetCDF library version: ', UNver
      write(*,*) 'NetCDF access (libUN) : ', NCDFver

C +--Read time-constant for the models...
C +  ------------------------------------
      CALL DInfo(1,'Read Time-Constants...')


C +-  1) MAR
C +   - - - -
C     Open a file (the first, because nMAR init=1):
      nfile = MARfile(nMAR)
      WRITE(*,*) 'Opening ',nfile
      CALL UNropen (nfile, idMAR, MARtit)

C +  Read Lat and Lon of MAR points
      CALL UNread
     &  (idMAR, 'lon', 0, 0, 1,1, mx,my,1,
     &    MARx, MARy, empty1, tmp_units, MARlon)
      CALL UNread
     &  (idMAR, 'lat', 0, 0, 1,1, mx,my,1,
     &    MARx, MARy, empty1, tmp_units, MARlat)
      CALL UNread
     &  (idMAR,'level',0, 0, 1,1, mz,1,1,
     &   tmp_mz,empty1,empty1,tmp_units, MARsig)

C +  MAR soil types
      CALL UNread
     &  (idMAR, 'isol', 0,0, 1,1, mx,my,1,
     &    MARx, MARy, empty1, tmp_units, MARisol)

      CALL UNread
     &  (idMAR, 'sh'  , 0,0, 1,1, mx,my,1,
     &    MARx, MARy, empty1, tmp_units, shMAR  )

C +   *NB: the file is left open;                  
C     (ok because = the first. 
C      If closed, must set idMAR=-1 again)

C +-  2) Large-scale (RLS)
C +   - - - - - - - - - - -

C +   Open a file (the first):
      nfile = RLSfile(nRLS)
      WRITE(*,*) 'Opening ',nfile
      CALL UNropen (nfile, idRLS, RLStit)

C +   RLS surface height & lon,lat coordinate:
      CALL UNread
     &  (idRLS, 'SH', 0,0, 1,1, LSni,LSnj,1,
     &    RLSlon,RLSlat,empty1, tmp_units, shRLS)

C +   *NB: the file is left open;
C     (ok because = the first.
C      If closed, must set idRLS=-1 again)


C +---Get the definitions of averaging sub-regions.
C +   ---------------------------------------------
      CALL SetRegs (MARisol,shMAR,MARlon,MARlat,InReg,RegMsk,
     .              wmSta,StaMsk)

C +---Perform specific initialisations: 
C +   ---------------------------------
C +.. Z500 filter:
      CALL InitZFilt

C +---Initialise time-loop.
C +   ---------------------
      idt   = 0
C +   ^The main counter for used time-steps. Begins at 0.

      itlast= (jdDLN*24+jhDLN) / jhSTP
C +   ^ the number of the last time-step (total= itlast+1) 

C +---Create the output file.
C +   =======================
C +   Note: currently left open until the end.
C +   However, it might be closed and open with UN.
      CALL DInfo(1,'Create nc output...')

      CALL CreateOutNC (LoutHOR,LoutBND,LoutVER,LoutMIS,LoutSRF,
     &                  OUTfile, iyrBEG,mmaBEG,jdaBEG,jhuBEG, 
     &                  MARlon, MARlat, MARsig,
     &                  MARx, MARy, RLSlon, RLSlat, 
     &                  jhSTP, itlast, iOUTnc)

      CALL DInfo(1,'Begin time loop...')
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +---Begin of time loop.
  10  CONTINUE 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C +---Compute the current date
C +   ========================
      iadtime= idt*jhSTP
C +   ^Time elapsed from initial, in hours.
      idt1   = idt+1
C +   ^The number of the current time step, starting at 1.
      
      CALL MARCalend (iyrBEG,mmaBEG,jdaBEG,jhuBEG, iadtime,
     .                iyrCUR,mmaCUR,jdaCUR,jhuCUR)

      Cdate= (mmaCUR*100 + jdaCUR) *100 + jhuCUR

      WRITE(*,110) iyrCUR,mmaCUR,jdaCUR,jhuCUR
110   FORMAT(' >Reading data for Year: ',
     .       I4,', Month: ',I2,', Day: ',I2,', Hour: ',I2)

C +  -Write the current date in the NetCDF file:
      CALL UNwrite (iOUTnc, 'date', idt1, 1,1,1, Cdate)

C +---Open Files
C +   ==========

C +---RLS files
C +...Open NetCDF Large Scale Raw Data file
120   IF (idRLS.EQ. -1) THEN 
        nfile = RLSfile(nRLS)
         WRITE(*,*) 'Opening ',nfile
        CALL UNropen (nfile, idRLS, RLStit)
      ENDIF
   
C +...inquire about the position of the requested date in file:
      CALL UNgindx (idRLS, "date", Cdate, Rdate, Fdate, itR)
C +...If the requested date was not found, close file and
C     select next file:
      IF (ABS(Rdate-Cdate).GT.0.5) THEN
        nRLS = nRLS + 1
        IF (nRLS .GT. ntRLS) THEN
           write(*,*) 'MAPOST - error: Requested date is not'
           write(*,*) '  available in RLS input, date=', Cdate
           STOP
        ENDIF
        CALL NCCLOS(idRLS, iError)
        idRLS = -1
        GOTO 120
      ENDIF

      IF(icheck.ge.3) THEN
         WRITE(*,*)'RLS Req. / Read date:', Cdate, Rdate
      END IF

C +---MAR files
C +...Open NetCDF MAR file
130   IF (idMAR.EQ. -1) THEN 
        nfile = MARfile(nMAR)
        WRITE(*,*) 'Opening ',nfile
        CALL UNropen (nfile, idMAR, MARtit)
      ENDIF
   
C +...inquire about the position of the requested date in file:
      CALL UNgindx (idMAR, "date", Cdate, Rdate, Fdate, itM)
C +...If the requested date was not found, close file and
C     select next file:
cXF   
      if ( Rdate == 10100.00 .and. Cdate == 123118.0) THEN
           Rdate = Cdate
      end if ! 31 december problem
cXF  
      IF (ABS(Rdate-Cdate).GT.0.5) THEN
        nMAR = nMAR + 1
        IF (nMAR .GT. ntMAR) THEN
           write(*,*) 'MAPOST - error: Requested date is not'
           write(*,*) '  available in MAR input, date=', Cdate
           STOP
        ENDIF
        CALL NCCLOS(idMAR, iError)
        idMAR = -1
        GOTO 130
      ENDIF

      IF(icheck.ge.3) THEN
         WRITE(*,*)'MAR Req. / Read date:', Cdate, Rdate
      END IF


C +---"Horizontal fields / dynamics"
C +   ==============================
      IF (LoutHOR) THEN

        CALL HDyn(idRLS, idMAR, itM, itR, istat,
     $            InReg, MARlon, MARlat, jhuCUR,
     $            iOUTnc,idt1,slpMAR)

      ENDIF
      
C +---Surface Stations observations
C +   ==============================

      IF (LoutSta) THEN
cXF
c       CALL SurSta(idRLS, idMAR, itM, itR, istat,
c    $              wmSta, MARlon, MARlat, jhuCUR,
c    $              iOUTnc,idt1)

      ENDIF

C +---Vertical Profiles
C +   =================
      IF (LoutVER) THEN

        CALL VerP(idRLS, idMAR, itM, itR, istat,
     $             InReg, MARlon, MARlat,
     $             iOUTnc,idt1)

      ENDIF

C +---local impact of boundaries (fn dbound)
C +   ======================================
      IF (LoutBND) THEN
        iptMKE=2 
        CALL dBound(idRLS, idMAR, itR, itM, istat, ioutNC, 
     $              slpMAR, idt1, iptMKE)

      ENDIF

C +---Surface (Precip, T2m)
C +   =====================
      IF (LoutSRF) THEN

        CALL Surf(idRLS, idMAR, itM, itR, istat,
     $            InReg, MARlon, MARlat, jhuCUR, mmaCUR,
     $            MAPtyp, iOUTnc,idt1)

      ENDIF


C +---Get ready for the next time-step...
C +   ===================================
C +   ** Increment time step:
      idt= idt + 1

C +   ** Remember if it was the last time step (see below)
      listat= istat

C +   ** Set "istat", i.e.
C +      tell routines that the time-serie goes on (=1)
C +      or ends (=2) at the NEXT time step
C +      -> routines must compute final values :
C +
      IF (idt.LT.itlast) THEN
        istat = 1
      ELSE
        istat = 2 
      ENDIF

C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +---End of time loop.
      IF (listat .NE. 2) GOTO 10
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      

C +--Write global outputs
C +  ====================
C +..(statistics must be writen inside subroutines, not here).

      CALL UNwrite (iOUTnc, 'isol   ', 1, mx, my, 1, MARisol)
      CALL UNwrite (iOUTnc, 'RegMsk ', 1, mx, my, 1, RegMsk)
      CALL UNwrite (iOUTnc, 'MAR_sh ', 1, mx, my, 1, shMAR)

      CALL INThor  (-1, RLSlon,RLSlat , shRLS,
     &                 MARlon, MARlat, shIRLS)      
      CALL UNwrite (iOUTnc, 'RLS_sh ', 1, mx, my, 1, shIRLS)


C +--NetCDF File Closure
C +  ===================
C +
      CALL UNclose(iOUTnc, iError)


      Write(*,*)('                                            ')
      Write(*,*)(' Execution terminated succesfully (I guess).')
      Write(*,*)(' Ite, missa est !                           ')
      Write(*,*)('                                            ')

      open (unit=1,status='replace',file='MPP.OK')
      write (1,*) "MAPOST OK"
      close(1)

      END
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
