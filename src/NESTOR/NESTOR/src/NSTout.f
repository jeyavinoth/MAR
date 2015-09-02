C   +-------------------------------------------------------------------+
C   |  Subroutine NSTout                          January 2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : Interpolated LSC (large-scale) fields                     |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: Creation of Netcdf file containing most of interpolated   |
C   | ^^^^^^^ meteorological fields (wind, temperature, humidity, ...)  |
C   |         This NetCDF File is adapted to IDL Graphic Software.      |
C   |                                                                   |
C   | Note that this routine requires the usual NetCDF library, and a   |
C   | complementary access to the library named 'libUN.a'.              |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE NSTout


      IMPLICIT NONE


C +---General variables
C +   =================

      include 'NSTdim.inc'
      include 'NSTvar.inc'
      include 'NESTOR.inc'

C +---Parameters
C +   ==========

      INTEGER i,j,k,it,NdimNC,MXdim,MX_var,NattNC,nnsl

      PARAMETER (NdimNC =   5)
C +...Number of defined spatial dimensions (exact)

      PARAMETER (MXdim  = 900)
C +...Maximum Number of all dims: recorded Time Steps
C +   and also maximum of spatial grid points for each direction. 

      PARAMETER (MX_var =  80)
C +...Maximum Number of Variables 

      PARAMETER (NattNC =   2)
C +...Number of real attributes given to all variables


C +---Local variables
C +   ===============

      INTEGER nbchar,ID__nc,Rcode,ipr_nc,NtotNC,npr_nc,
     .        INImma,INIjda,INIjhu,MMXstp,itotNC,year_1,year_2

      INTEGER*4 TMPdat

      INTEGER nDFdim(0:NdimNC),NvatNC(NattNC)

      REAL dateNC(MXdim),VALdim(MXdim,0:NdimNC),
     .     WK2D_1(mx,my),WK2D_2(mx,my),WK2D_3(mx,my),
     .     WK2D_4(mx,my),WK2D_5(mx,my),WK2D_6(mx,my),
     .     WK2D_7(mx,my),WK2D_8(mx,my),timeNC(MXdim),
     .     WK2D_9(mx,my,nvx),WK2D10(mx,my,nvx)

      CHARACTER*2   nustri(0:99)
      CHARACTER*13  NAMdim(0:NdimNC),nameNC(MX_var),
     .              SdimNC(4,MX_var),NAMrat(NattNC)
      CHARACTER*31  UNIdim(0:NdimNC),unitNC(MX_var)
      CHARACTER*17  suffix
      CHARACTER*50  lnamNC(MX_var)
      CHARACTER*90  fnamNC,fnam_U,fnam_V,fnam_T,fnam_Q,
     .              fnamSP,fnamSH
      CHARACTER*100 tit_NC
      CHARACTER*120 tmpINP

      LOGICAL Vfalse,Tferret


C +---Data
C +   ====

      DATA Tferret / .true. /  ! Time base for FERRET graphic tools

      DATA Vfalse  / .false. /

      DATA (nustri(i),i=0,99)
     .     /'00','01','02','03','04','05','06','07','08','09',
     .      '10','11','12','13','14','15','16','17','18','19',
     .      '20','21','22','23','24','25','26','27','28','29',
     .      '30','31','32','33','34','35','36','37','38','39',
     .      '40','41','42','43','44','45','46','47','48','49',
     .      '50','51','52','53','54','55','56','57','58','59',
     .      '60','61','62','63','64','65','66','67','68','69',
     .      '70','71','72','73','74','75','76','77','78','79',
     .      '80','81','82','83','84','85','86','87','88','89',
     .      '90','91','92','93','94','95','96','97','98','99'/

      nnsl = nsl


C +---Initial date of run
C +   ===================

C +         ******
       CALL DATcnv (RUNiyr,INImma,INIjda,INIjhu,DATini,Vfalse)
C +         ******


C +---Output directory
C +   ================

       nbchar=1

       DO i=1,60
        IF (NSTdir(i:i).ne.' ') nbchar=i
       ENDDO


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF (NSTmod.ne.'GRA') THEN

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output file name
C +   ================

      year_1 = INT(AINT(REAL(RUNiyr)/100.))
      year_2 = RUNiyr-(year_1*100)

      fnamNC = NSTdir(1:nbchar)// 'NST' // '.'
     .       // nustri(year_1) // nustri(year_2) // '.' 
     .       // nustri(INImma) // '.'
     .       // nustri(INIjda) // '.' 
     .       // nustri(INIjhu) // '.'
     .       // LABlio         // '.nc'

      write(6,*) 'Graphic output'
      write(6,*) '~~~~~~~~~~~~~~'
      write(6,*) 'File : ',fnamNC
      write(6,*)


C +---NetCDF File Initialization
C +   ==========================

      IF (DATini.eq.DATtim) THEN


C +---Output Title
C +   ------------

       tit_NC = 'NESTOR output'
     .        // ' - Mod: ' // NSTmod
     .        // ' - Exp: ' // LABlio
     .        // ' - '
     .        // nustri(year_1) // nustri(year_2) // '/' 
     .        // nustri(INImma) // '/'
     .        // nustri(INIjda) // '/' 
     .        // nustri(INIjhu)


C +---Number of time steps
C +   --------------------

       npr_nc = DATstp


C +---Check of size of temporary arrays
C +   ---------------------------------

       MMXstp = MXdim

       IF (npr_nc.gt.MMXstp)
     .  STOP '*** NSTout - ERROR : MXdim to low (temporally)***'

       IF (mx.gt.MMXstp.or.my.gt.MMXstp.or.(mz+1).gt.MMXstp)
     .  STOP '*** NSTout - ERROR : MXdim to low (spatially) ***'


C +---Time Variable (date)
C +   --------------------

       DO it=1,npr_nc

        TMPdat = DATini + (it-1)*DAT_dt

        timeNC(it)=REAL(TMPdat - 16664520+15*24)
C +...  HOURS since 1901-01-15 00:00:00

C +          ******
        CALL DATcnv (RUNiyr,RUNmma,RUNjda,RUNjhu,TMPdat,Vfalse)
C +          ******
        dateNC(it)=REAL(RUNjhu)+1.d2*REAL(RUNjda)+1.d4*REAL(RUNmma)

       ENDDO


C +---Define temporal and spatial dimensions    
C +   --------------------------------------

       DO it = 1,npr_nc
        VALdim(it,0) = REAL(INIjhu + (it-1)*DAT_dt)  ! hours
       ENDDO
       IF (Tferret) THEN
        DO it = 1,npr_nc
         VALdim(it,0) = timeNC(it)  ! Hours since 1901-01-15
        ENDDO
       ENDIF
       nDFdim(0)= npr_nc
c #UL  nDFdim(0)= 0
       NAMdim(0)= 'time'
       UNIdim(0)= 'HOURS since 1901-01-15 00:00:00'

       DO i=1,mx
        VALdim(i,1) = NSTgdx(i)
       ENDDO
       nDFdim(1)= mx
       NAMdim(1)= 'x'
       UNIdim(1)= 'km'
 
       DO j=1,my
        VALdim(j,2) = NSTgdy(j)
       ENDDO
       nDFdim(2)= my
       NAMdim(2)= 'y'
       UNIdim(2)= 'km'
 
       DO k=1,mz
        VALdim(k,3) = NSTgdz(k)
       ENDDO
       nDFdim(3)= mz
       NAMdim(3)= 'level'
       UNIdim(3)= '[sigma]'

       DO k = 1,nvx
        VALdim(k,4) = k
       ENDDO
       nDFdim(4)= nvx
       NAMdim(4)= 'sector'
       UNIdim(4)= '[index]'

       DO k = 1,nsl
        VALdim(k,5) = k
       ENDDO
       nDFdim(5)= nsl
       NAMdim(5)= 'soil'
       UNIdim(5)= '[index]'


C +---Variable's Choice (Table LSMvou.dat)
C +   ------------------------------------

       OPEN(unit=15,status='unknown',file='NSTvou.dat')

       itotNC = 0

980    CONTINUE
        READ (15,'(A120)',end=990) tmpINP

        IF (tmpINP(1:4).eq.'    ') THEN 
         itotNC = itotNC + 1
         READ (tmpINP,'(4x,5A9,A12,A25)')
     .        nameNC(itotNC)  ,SdimNC(1,itotNC),SdimNC(2,itotNC),
     .        SdimNC(3,itotNC),SdimNC(4,itotNC),unitNC(itotNC)  ,
     .        lnamNC(itotNC)
C +...        - nameNC: Name
C +...        - SdimNC: Names of Selected Dimensions (max.4/variable) 
C +...        - unitNC: Units
C +...        - lnamNC: Long_name, a description of the variable
        ENDIF

       GOTO 980
990    CONTINUE

       NtotNC = itotNC 
C +... NtotNC : Total number of variables writen in NetCDF file


C +---List of NetCDF attributes given to all variables
C +   ------------------------------------------------

       NAMrat(1) = 'actual_range'
       NvatNC(1) = 2

       NAMrat(2) = 'valid_range'
       NvatNC(2) = 2

C      NAMrat(NattNC) = '[var]_range'
C      NvatNC(NattNC) = 2
C      Used by IDL/INA but probably not by ferret... 
C      (purpose was to create animations)


C +---Automatic Generation of the NetCDF File Structure
C +   -------------------------------------------------

C +          *********
       CALL UNscreate (fnamNC,tit_NC,NdimNC,nDFdim,MXdim ,NAMdim,
     .                 UNIdim,VALdim,MX_var,NtotNC,nameNC,SdimNC,
     .                 unitNC,lnamNC,NattNC,NAMrat,NvatNC,ID__nc)
C +          *********


C +---Write Time - Constants fields
C +   -----------------------------

       DO j=1,my
       DO i=1,mx
        WK2D_3(i,j)=NST_ts(i,j,1,1)
        WK2D_4(i,j)=NST_ts(i,j,1,2)
        WK2D_5(i,j)=NST_sw(i,j,1,1)
        WK2D_6(i,j)=NST_sw(i,j,1,2)
        WK2D_7(i,j)=REAL(NSTsol(i,j))
        WK2D_8(i,j)=REAL(NSTtex(i,j))
        DO k=1,nvx
         WK2D_9(i,j,k)=REAL(NSTveg(i,j,k))
         WK2D10(i,j,k)=REAL(NSTsvt(i,j,k))
        ENDDO
       ENDDO
       ENDDO

C +         *******
       CALL UNwrite (ID__nc, 'DATE   ', 1  , npr_nc, 1 , 1  , dateNC)
       CALL UNwrite (ID__nc, 'LON    ', 1  , mx    , my, 1  , NST__x)
       CALL UNwrite (ID__nc, 'LAT    ', 1  , mx    , my, 1  , NST__y)
       CALL UNwrite (ID__nc, 'SH     ', 1  , mx    , my, 1  , NST_sh)
       CALL UNwrite (ID__nc, 'SOL    ', 1  , mx    , my, 1  , WK2D_7)
       CALL UNwrite (ID__nc, 'TEX    ', 1  , mx    , my, 1  , WK2D_8)
       CALL UNwrite (ID__nc, 'Z0     ', 1  , mx    , my, 1  , NST_z0)
       CALL UNwrite (ID__nc, 'R0     ', 1  , mx    , my, 1  , NST_r0)
       
       CALL UNwrite (ID__nc, 'VEG    ', 1  , mx    , my, nvx, WK2D_9)
       CALL UNwrite (ID__nc, 'SVT    ', 1  , mx    , my, nvx, WK2D10)

       CALL UNwrite (ID__nc, 'ALB    ', 1  , mx    , my, 1  , NSTalb)
       CALL UNwrite (ID__nc, 'DSA    ', 1  , mx    , my, 1  , NSTdsa)
       
       CALL UNwrite (ID__nc, 'DV1    ', 1  , mx    , my, 1  , NSTdv1)
       CALL UNwrite (ID__nc, 'DV2    ', 1  , mx    , my, 1  , NSTdv2)

       CALL UNwrite (ID__nc, 'TS     ', 1  , mx    , my, nsl, NST_ts)
       CALL UNwrite (ID__nc, 'SW     ', 1  , mx    , my, nsl, NST_sw)

       CALL UNwrite (ID__nc, 'RES    ', 1  , mx    , my, 1  , NSTres)
C +         *******

       DO j=1,my
       DO i=1,mx
       DO k=1,nvx
        WK2D_9(i,j,k)=REAL(NSTvfr(i,j,k))
        WK2D10(i,j,k)=REAL(NSTsfr(i,j,k))
       ENDDO
       ENDDO
       ENDDO

C +         *******
       CALL UNwrite (ID__nc, 'FRC    ', 1  , mx    , my, 1  , NSTfrc)
       CALL UNwrite (ID__nc, 'VFR    ', 1  , mx    , my, nvx, WK2D_9)
       CALL UNwrite (ID__nc, 'SFR    ', 1  , mx    , my, nvx, WK2D10)
C +         *******

      ELSE

C +---Re-Open file if already created.
C +   ================================


C +         *******
       CALL UNwopen (fnamNC,ID__nc)
C +         *******

      END IF


C +---Write time-dependent variables
C +   ==============================

      ipr_nc=(DATtim-DATini)/DAT_dt + 1

      DO j=1,my
      DO i=1,mx
       WK2D_1(i,j) = 0.0
       WK2D_2(i,j) = 0.0
       DO k=1,nvx-1
        WK2D_1(i,j) = WK2D_1(i,j) 
     .              + NSTglf(i,j,k)*REAL(NSTsfr(i,j,k))/100.
        WK2D_2(i,j) = WK2D_2(i,j) 
     .              + NSTlai(i,j,k)
     .               *NSTglf(i,j,k)*REAL(NSTsfr(i,j,k))/100.
       ENDDO
      ENDDO
      ENDDO

C +        *******
      CALL UNwrite (ID__nc, 'UU     ', ipr_nc, mx, my, mz , NST__u)
      CALL UNwrite (ID__nc, 'VV     ', ipr_nc, mx, my, mz , NST__v)
      CALL UNwrite (ID__nc, 'TT     ', ipr_nc, mx, my, mz , NST__t)
      CALL UNwrite (ID__nc, 'PT     ', ipr_nc, mx, my, mz , NST_pt)
      CALL UNwrite (ID__nc, 'RH     ', ipr_nc, mx, my, mz , NST_rh)
      CALL UNwrite (ID__nc, 'QQ     ', ipr_nc, mx, my, mz , NST_qv)
      CALL UNwrite (ID__nc, 'ZZ     ', ipr_nc, mx, my, mz , NST_zz)
      CALL UNwrite (ID__nc, 'SP     ', ipr_nc, mx, my, 1  , NST_sp)
      CALL UNwrite (ID__nc, 'ST     ', ipr_nc, mx, my, 1  , NST_st)
      CALL UNwrite (ID__nc, 'SST    ', ipr_nc, mx, my, 1  , NSTsst)
      CALL UNwrite (ID__nc, 'SIC    ', ipr_nc, mx, my, 1  , NSTsic)
      CALL UNwrite (ID__nc, 'EWC    ', ipr_nc, mx, my, 1  , NSTewc)
      CALL UNwrite (ID__nc, 'NDV    ', ipr_nc, mx, my, 1  , NSTndv)
      CALL UNwrite (ID__nc, 'LAI    ', ipr_nc, mx, my, nvx, NSTlai)
      CALL UNwrite (ID__nc, 'GLF    ', ipr_nc, mx, my, nvx, NSTglf)
      CALL UNwrite (ID__nc, 'EFRV   ', ipr_nc, mx, my, 1  , WK2D_1)
      CALL UNwrite (ID__nc, 'ELAI   ', ipr_nc, mx, my, 1  , WK2D_2)
      CALL UNwrite (ID__nc, 'ICE    ', ipr_nc, mx, my, 1  , NSTice)
C +        *******


C +---NetCDF File Closure
C +   ===================

C +        ******
      CALL NCCLOS (ID__nc, Rcode)
C +        ******


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDIF

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF (NSTmod.eq.'GRA') THEN

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Open unformatted files
C +   ======================

      IF (DATini.eq.DATtim) THEN

       year_1 = INT(AINT(REAL(RUNiyr)/100.))
       year_2 = RUNiyr-(year_1*100)

       suffix = nustri(year_1) // nustri(year_2) // '.'
     .       // nustri(INImma) // '.'
     .       // nustri(INIjda) // '.'
     .       // nustri(INIjhu) // '.bin'

       fnam_U = NSTdir(1:nbchar)// 'U__' // suffix
       fnam_V = NSTdir(1:nbchar)// 'V__' // suffix
       fnam_T = NSTdir(1:nbchar)// 'T__' // suffix
       fnam_Q = NSTdir(1:nbchar)// 'Q__' // suffix
       fnamSP = NSTdir(1:nbchar)// 'SP_' // suffix
       fnamSH = NSTdir(1:nbchar)// 'SH_' // suffix

       write(6,*) 'Graphic output'
       write(6,*) '~~~~~~~~~~~~~~'
       write(6,*) 'Files : ',fnam_U
       write(6,*) '        ',fnam_V
       write(6,*) '        ',fnam_T
       write(6,*) '        ',fnam_Q
       write(6,*) '        ',fnamSP
       write(6,*) '        ',fnamSH
       write(6,*)

       OPEN (unit=61,status='unknown',form='unformatted',file=fnam_U)
       OPEN (unit=62,status='unknown',form='unformatted',file=fnam_V)
       OPEN (unit=63,status='unknown',form='unformatted',file=fnam_T)
       OPEN (unit=64,status='unknown',form='unformatted',file=fnam_Q)
       OPEN (unit=65,status='unknown',form='unformatted',file=fnamSP)
       OPEN (unit=66,status='unknown',form='unformatted',file=fnamSH)

      ENDIF


C +---Write variables in output file
C +   ==============================

      write(6,*) 'Graphic output'
      write(6,*) '~~~~~~~~~~~~~~'
      write(6,*) 'Append unformatted files'
      write(6,*)

      WRITE(61) NST__u
      WRITE(62) NST__v
      WRITE(63) NST__t
      WRITE(64) NST_qv
      WRITE(65) NST_sp
      WRITE(66) NST_sh


C +---Close unformatted files
C +   =======================

      IF (DATtim.eq.DATfin) THEN
       CLOSE (unit=61)
       CLOSE (unit=62)
       CLOSE (unit=63)
       CLOSE (unit=64)
       CLOSE (unit=65)
       CLOSE (unit=66)
      ENDIF
       

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDIF

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
