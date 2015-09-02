C   +-------------------------------------------------------------------+
C   |  Subroutine WGEout                        September 2001  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : Wind gust estimates + wind distribution                   |
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


      SUBROUTINE WGEout


      IMPLICIT NONE


C +---General variables
C +   =================

      include 'NSTdim.inc'
      include 'NESTOR.inc'
      include 'NSTvar.inc'
      include 'WGEvar.inc'

C +---Parameters
C +   ==========

      INTEGER i,j,k,l,it,NdimNC,MXdim,MX_var,NattNC

      PARAMETER (NdimNC =   4)
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

      INTEGER nbchar,ID__nc,Rcode,ipr_nc,NtotNC,
     .        INIiyr,INImma,INIjda,INIjhu,MMXstp,npr_nc,
     .        itotNC,year_1,year_2

      INTEGER*4 TMPdat

      INTEGER nDFdim(0:NdimNC),NvatNC(NattNC)

      REAL dateNC(MXdim),VALdim(MXdim,0:NdimNC),
     .     WK2D_1(mx,my)

      CHARACTER*2   nustri(0:99)
      CHARACTER*13  NAMdim(0:NdimNC),nameNC(MX_var),
     .              SdimNC(4,MX_var),NAMrat(NattNC)
      CHARACTER*16  UNIdim(0:NdimNC),unitNC(MX_var)
      CHARACTER*50  lnamNC(MX_var)
      CHARACTER*80  fnamNC
      CHARACTER*90  tit_NC

      LOGICAL Vfalse


C +---Data
C +   ====

      DATA Vfalse / .false. /

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


C +---Output file name
C +   ================

       year_1 = INT(AINT(REAL(RUNiyr)/100.))
       year_2 = RUNiyr-(year_1*100)

       fnamNC = NSTdir(1:nbchar)// NSTmod // '_WGE_'
     .        // nustri(year_1) // nustri(year_2) // '.' 
     .        // nustri(INImma) // '.'
     .        // nustri(INIjda) // '.' 
     .        // nustri(INIjhu) // '.nc'

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
     .        // ' - Exp: ' // NSTmod
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
     .  STOP '*** PRCout - ERROR : MXdim to low (temporally)***'

       IF (mx.gt.MMXstp.or.my.gt.MMXstp.or.(mz+1).gt.MMXstp)
     .  STOP '*** PRCout - ERROR : MXdim to low (spatially) ***'


C +---Time Variable (date)
C +   --------------------

       DO it=1,npr_nc

        TMPdat = DATini + (it-1)*DAT_dt
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
       nDFdim(0)= npr_nc
       NAMdim(0)= 'time'
       UNIdim(0)= 'days'

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
        VALdim(k,3) = 10./3.6*REAL(k)
       ENDDO
       nDFdim(3)= gc
       NAMdim(3)= 'classes'
       UNIdim(3)= '[sigma]'

       DO k = 1,mw
        VALdim(k,4) = REAL(k)
       ENDDO
       nDFdim(4)= mw
       NAMdim(4)= 'sector'
       UNIdim(4)= '[index]'



C +---Variable's Choice
C +   -----------------

       itotNC = 1
       nameNC(  itotNC)='date'
       SdimNC(1,itotNC)='-'
       SdimNC(2,itotNC)='-'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='MMDDHH'
       lnamNC(  itotNC)='Date (MM/DD/HH)'

       itotNC = 2
       nameNC(  itotNC)='lon'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='-'
       unitNC(  itotNC)='degrees'
       lnamNC(  itotNC)='Longitude'

       itotNC = 3
       nameNC(  itotNC)='lat'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='-'
       unitNC(  itotNC)='degrees'
       lnamNC(  itotNC)='Latitude'

       itotNC = 4
       nameNC(  itotNC)='sh'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='-'
       unitNC(  itotNC)='m'
       lnamNC(  itotNC)='Surface height'

       itotNC = 5
       nameNC(  itotNC)='isol'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='-'
       unitNC(  itotNC)='-'
       lnamNC(  itotNC)='Surface type'

       itotNC = 6
       nameNC(  itotNC)='WGEest'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='m/s'
       lnamNC(  itotNC)='Estimated wind gusts'

       itotNC = 7
       nameNC(  itotNC)='WGElwb'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='m/s'
       lnamNC(  itotNC)='Lower bound on wind gusts'

       itotNC = 8
       nameNC(  itotNC)='WGEupb'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='m/s'
       lnamNC(  itotNC)='Upper bound on wind gusts'

       itotNC = 9
       nameNC(  itotNC)='WGE_zi'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='m'
       lnamNC(  itotNC)='Inversion height'

       itotNC = 10
       nameNC(  itotNC)='WNDsta'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='classes'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='-'
       lnamNC(  itotNC)='Wind distribution'

       itotNC = 11
       nameNC(  itotNC)='WGEsta'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='classes'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='-'
       lnamNC(  itotNC)='Gust distribution'

C +...        - nameNC: Name
C +...        - SdimNC: Names of Selected Dimensions (max.4/variable) 
C +...        - unitNC: Units
C +...        - lnamNC: Long_name, a description of the variable

       NtotNC = itotNC 
C +... NtotNC : Total number of variables writen in NetCDF file


C +---List of NetCDF attributes given to all variables
C +   ------------------------------------------------

       NAMrat(1) = 'actual_range'
       NvatNC(1) = 2

       NAMrat(NattNC) = '[var]_range'
       NvatNC(NattNC) = 2


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
        WK2D_1(i,j)=REAL(NSTsol(i,j))
       ENDDO
       ENDDO

C +         *******
       CALL UNwrite (ID__nc, 'date  ', 1  , npr_nc, 1 , 1 , dateNC)
       CALL UNwrite (ID__nc, 'lon   ', 1  , mx    , my, 1 , NST__x)
       CALL UNwrite (ID__nc, 'lat   ', 1  , mx    , my, 1 , NST__y)
       CALL UNwrite (ID__nc, 'sh    ', 1  , mx    , my, 1 , NST_sh)
       CALL UNwrite (ID__nc, 'isol  ', 1  , mx    , my, 1 , WK2D_1)
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

C +        *******
      CALL UNwrite (ID__nc, 'WGEest ', ipr_nc, mx, my, 1 , WGEest)
      CALL UNwrite (ID__nc, 'WGElwb ', ipr_nc, mx, my, 1 , WGElwb)
      CALL UNwrite (ID__nc, 'WGEupb ', ipr_nc, mx, my, 1 , WGEupb)
      CALL UNwrite (ID__nc, 'WGE_zi ', ipr_nc, mx, my, 1 , WGE_zi)
      CALL UNwrite (ID__nc, 'WNDsta ', ipr_nc, mx, my, gc, WNDsta)
      CALL UNwrite (ID__nc, 'WGEsta ', ipr_nc, mx, my, gc, WGEsta)
C +        *******


C +---NetCDF File Closure
C +   ===================

C +        ******
      CALL NCCLOS (ID__nc, Rcode)
C +        ******


      RETURN
      END
