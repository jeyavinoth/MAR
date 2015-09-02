C   +-------------------------------------------------------------------+
C   |  Subroutine PRCout                            April 2003  NESTING |
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


      SUBROUTINE PRCout


      IMPLICIT NONE


C +---General variables
C +   =================

      include 'NSTdim.inc'
      include 'NSTvar.inc'
      include 'NESTOR.inc'

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

      REAL VALdim(MXdim,0:NdimNC),WK2D_1(mx,my),dateNC

      CHARACTER*2   nustri(0:99)
      CHARACTER*13  NAMdim(0:NdimNC),nameNC(MX_var),
     .              SdimNC(4,MX_var),NAMrat(NattNC)
      CHARACTER*31  UNIdim(0:NdimNC),unitNC(MX_var)
      CHARACTER*50  lnamNC(MX_var)
      CHARACTER*80  fnamNC
      CHARACTER*90  tit_NC

      LOGICAL Vfalse

      REAL njmoGE(12),njmbGE(12),njyrGE(12),njybGE(12)
      REAL PR1T(mx,my),PR2T(mx,my),PR3T(mx,my),PR0T(mx,my)

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
     
      data  njmoGE
     .     /31,28,31,30, 31, 30, 31, 31, 30, 31, 30, 31/
      data  njmbGE
     .     /0, 1, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0/
      data  njyrGE
     .     /0,31,59,90,120,151,181,212,243,273,304,334/
      data  njybGE
     .     /0, 0, 1, 1,  1,  1,  1,  1,  1,  1,  1,  1/


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

       fnamNC = NSTdir(1:nbchar)// NSTmod // '_DSG_'
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

       DO j=1,my
       DO i=1,mx
        PR1T(i,j)=0.
        PR2T(i,j)=0.
        PR3T(i,j)=0.
        PR0T(i,j)=0.
       ENDDO
       ENDDO


C +---Output Title
C +   ------------

       tit_NC = 'NESTOR-DSG output'
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


C +---Define temporal and spatial dimensions    
C +   --------------------------------------

       DO it = 1,npr_nc
        VALdim(it,0)    = real
     .                  ((351  +(RUNiyr  -1902) *365   ! Nb Days before iyrrGE
     .                  +(RUNiyr  -1901) /  4          ! Nb Leap Years
     .                  + njyrGE(INImma)               ! Nb Days before mmarGE
     .                  + njybGE(INImma)               ! (including Leap Day)
     .                  * max(0,1-mod(RUNiyr,4))       !
     .                  + INIjda     -1      )*  24    !
     .                  + INIjhu )                     !
     .                  + (it-1)*DAT_dt                ! hours
       ENDDO
       nDFdim(0)= 0 
       NAMdim(0)= 'time'
       UNIdim(0)= 'HOURS since 1901-01-15 00:00:00'

       DO i=1,mx
        VALdim(i,1) = NSTgdx(i)/1000
       ENDDO
       nDFdim(1)= mx
       NAMdim(1)= 'x'
       UNIdim(1)= 'km'
 
       DO j=1,my
        VALdim(j,2) = NSTgdy(j)/1000
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

       DO k = 1,mw
        VALdim(k,4) = k
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
       nameNC(  itotNC)='PR1'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='mm'
       lnamNC(  itotNC)='Dsg. Precip. (no cons.)'

       itotNC = 7
       nameNC(  itotNC)='PR2'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='mm'
       lnamNC(  itotNC)='Dsg. Precip. (gbl. cons.)'

       itotNC = 8
       nameNC(  itotNC)='PR3'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='mm'
       lnamNC(  itotNC)='Dsg. Precip. (lcl. & glb. cons.)'

       itotNC = 9
       nameNC(  itotNC)='PR0'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='-'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='mm'
       lnamNC(  itotNC)='MAR  Precip.'
 
       itotNC = 10
       nameNC(  itotNC)='TT'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='level'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='deg K'
       lnamNC(  itotNC)='MAR tairDY'
       
       itotNC = 11
       nameNC(  itotNC)='QQ'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='level'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='kg/kg'
       lnamNC(  itotNC)='MAR qvDY'
       
       itotNC = 12
       nameNC(  itotNC)='UU'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='level'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='m/s'
       lnamNC(  itotNC)='MAR uairDY'
       
       itotNC = 13
       nameNC(  itotNC)='VV'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='level'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='m/s'
       lnamNC(  itotNC)='MAR vairDY'
       
       itotNC = 14
       nameNC(  itotNC)='WW'
       SdimNC(1,itotNC)='x'
       SdimNC(2,itotNC)='y'
       SdimNC(3,itotNC)='level'
       SdimNC(4,itotNC)='time'
       unitNC(  itotNC)='m/s'
       lnamNC(  itotNC)='MAR wairDY'

C +...        - nameNC: Name
C +...        - SdimNC: Names of Selected Dimensions (max.4/variable) 
C +...        - unitNC: Units
C +...        - lnamNC: Long_name, a description of the variable

       NtotNC = itotNC 
       NtotNC = 9
C +... NtotNC : Total number of variables writen in NetCDF file


C +---List of NetCDF attributes given to all variables
C +   ------------------------------------------------

       NAMrat(1) = 'actual_range'
       NvatNC(1) = 2

       NAMrat(NattNC) = '[var]_range'
       NvatNC(NattNC) = 2


C +---Automatic Generation of the NetCDF File Structure
C +   -------------------------------------------------

C +         *********
       CALL UNscreate (fnamNC,tit_NC,NdimNC,nDFdim,MXdim ,NAMdim,
     .                 UNIdim,VALdim,MX_var,NtotNC,nameNC,SdimNC,
     .                 unitNC,lnamNC,NattNC,NAMrat,NvatNC,ID__nc)
C +         *********


C +---Write Time - Constants fields
C +   -----------------------------

       DO j=1,my
       DO i=1,mx
        WK2D_1(i,j)=REAL(NSTsol(i,j))
       ENDDO
       ENDDO

C +         *******
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


C +---Time Variable (date)
C +   ====================

      ipr_nc=(DATtim-DATini)/DAT_dt + 1
      TMPdat = DATini + (ipr_nc-1)*DAT_dt
C +          ******
      CALL DATcnv (RUNiyr,RUNmma,RUNjda,RUNjhu,TMPdat,Vfalse)
C +          ******
      dateNC=REAL(RUNjhu)+1.d2*REAL(RUNjda)+1.d4*REAL(RUNmma)


C +---Write time-dependent variables
C +   ==============================

      DO j=1,my
      DO i=1,mx
       PR1T(i,j)=NSTpr1(i,j)*1000.+PR1T(i,j)
       PR2T(i,j)=NSTpr2(i,j)*1000.+PR2T(i,j)
       PR3T(i,j)=NSTpr3(i,j)*1000.+PR3T(i,j)
       PR0T(i,j)=NSTIpr(i,j)*1000.+PR0T(i,j)      
      ENDDO
      ENDDO

C +        *******

      CALL UNwrite (ID__nc, 'time', ipr_nc,  1,  1,  1, VALdim(ipr_nc,0)
     .) 
      CALL UNwrite (ID__nc, 'date', ipr_nc,  1,  1,  1, dateNC         )
      CALL UNwrite (ID__nc, 'PR0 ', ipr_nc, mx, my,  1, PR0T           )
      CALL UNwrite (ID__nc, 'PR1 ', ipr_nc, mx, my,  1, PR1T           )
      CALL UNwrite (ID__nc, 'PR2 ', ipr_nc, mx, my,  1, PR2T           )
      CALL UNwrite (ID__nc, 'PR3 ', ipr_nc, mx, my,  1, PR3T           )
      CALL UNwrite (ID__nc, 'TT  ', ipr_nc, mx, my, mz, NST__t         )
      CALL UNwrite (ID__nc, 'QQ  ', ipr_nc, mx, my, mz, NST_qv         )
      CALL UNwrite (ID__nc, 'UU  ', ipr_nc, mx, my, mz, NST__u         )
      CALL UNwrite (ID__nc, 'VV  ', ipr_nc, mx, my, mz, NST__v         )
      CALL UNwrite (ID__nc, 'WW  ', ipr_nc, mx, my, mz, NST__w         ) 
C +        *******

      
C +---NetCDF File Closure
C +   ===================

C +        ******
      CALL NCCLOS (ID__nc, Rcode)
C +        ******


      RETURN
      END
