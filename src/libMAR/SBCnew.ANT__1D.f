      subroutine SBCnew

C +------------------------------------------------------------------------+
C | MAR  INPUT               Snow Characteristics          31-01-2006  MAR |
C |     OUTPUT     Blowing   Snow Characteristics                          |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C | PRECISION: Conversion from real   to real   under vi                   |
C | ^^^^^^^^^  :.,$g/real  /s//real*8/g                                    |
C |            Conversion from real   to real   under vi                   |
C |            :.,$g/real[*]8/s//real  /g                                  |
C |                                                                        |
C |                                                                        |
C |                                                                        |
C |                                                                        |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARCTR.inc'
      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
      include 'MAR_LB.inc'

      include 'MAR_DY.inc'
      include 'MAR_HY.inc'
      include 'MAR_CA.inc'
      include 'MAR_TE.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'

      real              uairLS(mz),vairLS(mz),pktaLS(mz),qvLS(mz)
      real              uair_0(mz),vair_0(mz),pkta_0(mz),qv_0(mz)
      real              dVair                ,dTair
      common/SBCnew_loc/uairLS    ,vairLS    ,pktaLS    ,qvLS
     .                 ,uair_0    ,vair_0    ,pkta_0    ,qv_0

      integer           ipr2nc,npr2nc,n_iter
      common/Out2nc_arg/ipr2nc,npr2nc,n_iter


C +--Local   Variables
C +  =================

      logical     relxLS
      integer     isn   ,isl   ,misl_2,nisl_2,n
      real        dz_dSN(-nsno:1)
      real        Dendri,Spheri,GrainS,Densit,SnowEq,Age
      REAL        TimeNC
      real        argsin,sinarg


C +--DATA
C +  ====

      data        relxLS/ .true./
      data        TimeNC/  60.00/
C +...            TimeNC:Time Interval between IO on NetCDF File


C +--Snow Model Initialisation
C +  =========================

      IF (itexpe.eq.       0)                                     THEN ! CTR

          open (unit=30,status='unknown',file='SBCnew.data')
          rewind     30
                read(30,'(f3.1)') dVair
                read(30,'(f3.1)') dTair
          close(unit=30)


C +--Surface Air Temperature
C +  -----------------------

          INCLUDE 'SBCnew.dat'

          DO j=1,my
          DO i=1,mx
                tairSL(i,j)         = pktaDY(i,j,mz)       !            [K]
     .                  *(pstDYn(i,j)+ptopDY)**cap         !
          END DO
          END DO


C +--Vertical Discretization
C +  -----------------------

          DO isl= 0,-nsno+1,-1
            misl_2 =     -mod(isl,2)
            nisl_2 =         -isl/2
            dz_dSN(isl)=(((1-misl_2) * 0.001
     .                    +  misl_2  * 0.003) * 10**(nisl_2)) * 10.
            dz_dSN(isl)=min(dz_dSN(isl),0.20)
          END DO


C +--Large Scale Forcing
C +  -------------------

        IF (relxLS)                                               THEN ! CTR
          DO k=1,mz
            uairLS(k) = uairDY(imez,jmez,k)
            vairLS(k) = vairDY(imez,jmez,k)
            pktaLS(k) = pktaDY(imez,jmez,k)
              qvLS(k) =   qvDY(imez,jmez,k)
            uair_0(k) = uairDY(imez,jmez,k)
            vair_0(k) = vairDY(imez,jmez,k)
            pkta_0(k) = pktaDY(imez,jmez,k)
              qv_0(k) =   qvDY(imez,jmez,k)
          END DO
          OPEN (unit=30,status='unknown',form='unformatted',
     .          file='SBCnew.DAT')
          REWIND     30
          WRITE     (30) uairLS,vairLS,pktaLS,qvLS
          CLOSE(UNIT=30)
        END IF                                                         ! CTR

          DO j=1,my
          DO i=1,mx
          DO n=1,nsx


C +--Over the Ocean
C +  --------------

            IF (isolSL(i,j).le.2)                                   THEN
                nssSNo(i,j,n)      =   0                   ! Nb Snow/Ice Lay.
              DO isn = 1,nsno  
                nhsSNo(i,j,n,isn)  =   0                   !            [-]
                dzsSNo(i,j,n,isn)  =   0.                  !            [m]
                rosSNo(i,j,n,isn)  =   0.                  !        [kg/m3]
                wasSNo(i,j,n,isn)  =   0.                  !        [kg/kg]
                tisSNo(i,j,n,isn)  =   0.                  !            [K]
                g1sSNo(i,j,n,isn)  =   0.                  ! [-]        [-]
                g2sSNo(i,j,n,isn)  =   0.                  ! [-] [0.0001 m]
                agsSNo(i,j,n,isn)  =   0.                  !          [day]
              END DO
              DO isn = 1,llx
                TsolTV(i,j,n,isn)  =  SST_LB(i,j)          !            [K]
                eta_TV(i,j,n,isn)  =   1.
              END DO


C +--Over the Continent
C +  ------------------

            ELSE
                  nssSNo(i,j,n)      =   nsno                ! Nb Snow/Ice Lay.
              DO isn = 1,nsno                                !
                  nhsSNo(i,j,n,isn)  =   0                   !            [-]
                  dzsSNo(i,j,n,isn)  =  dz_dSN(1-isn)        !            [m]
                  rosSNo(i,j,n,isn)  =  Densit               !        [kg/m3]
                  wasSNo(i,j,n,isn)  =   0.                  !        [kg/kg]
                  tisSNo(i,j,n,isn)  =  tairSL(i,j)          !            [K]
                IF (Dendri.lt.0.)     THEN                   !
                  g1sSNo(i,j,n,isn)  =  Dendri               ! [-]        [-]
                  g2sSNo(i,j,n,isn)  =  Spheri               ! [-]        [-]
                ELSE                                         !
                  g1sSNo(i,j,n,isn)  =  Spheri               ! [-]        [-]
                  g2sSNo(i,j,n,isn)  =  GrainS               ! [-] [0.0001 m]
                END IF                                       !
                  agsSNo(i,j,n,isn)  =  Age                  !          [day]
              END DO
              DO isn = 1,llx
                  TsolTV(i,j,n,isn)  =  tairSL(i,j)          !            [K]
                  eta_TV(i,j,n,isn)  =   0.
              END DO
                  isolTV(i,j)        =  12
                  iwafTV(i,j)        =   0
                  AlbSTV(i,j)        =   0.55
                  ivegTV(i,j,n)      =   0
                  alaiTV(i,j,n)      =   0.
                  glf_TV(i,j,n)      =   0.
            END IF


C +--General
C +  -------

                TsrfSL(i,j,n)      =  tairSL(i,j)          ! Surface Temperat.
                TvegTV(i,j,n)      =  tairSL(i,j)          ! Vegetat.Temperat.

                issSNo(i,j,n)      =   0                   ! Nb Supr.Ice Layer
                nisSNo(i,j,n)      =   0                   ! Nb      Ice Layer
                snohSN(i,j,n)      =   0.                  ! Snow Buffer Layer
                SWaSNo(i,j,n)      =   0.                  ! Snow Surfic.Water
          END DO
                SLsrfl(i,j,1)      =   1.                  ! Surface Fraction
                ifraTV(i,j,1)      = 100                   ! Surface Fraction
          END DO
          END DO
      END IF                                                           ! CTR


C +--IO FILES
C +  ========

      IF (iterun.eq.       0)                                     THEN ! CTR


C +--Large Scale "Climatology"
C +  -------------------------

        IF (relxLS)                                               THEN ! CTR
          OPEN (unit=31,status='OLD',form='unformatted',
     .          file='SBCnew.DAT')
          REWIND     31
          READ      (31) uairLS,vairLS,pktaLS,qvLS
          CLOSE(UNIT=31)
        END IF                                                         ! CTR


C +--Animation NetCDF File
C +  ---------------------

c #VR     open(unit=33,status='new',file='_DUMP_.MAR')
c #VR     rewind    33
          n_iter = TimeNC / dt
C +...    n_iter : Number of Time steps between IO

          ipr2nc =                   1
          npr2nc = nterun / n_iter + 1
C +...    npr2nc : Number of                    IO

C +       ***********
          call Out2nc
C +       ***********

      ELSE                                                             ! CTR
        IF (relxLS)                                               THEN ! CTR
          argsin        = 2.*pi*iterun*dt/(3.*86400.)
          sinarg        =   sin(argsin)
        DO k=1,mz
          uairLS(    k) = uair_0(k) + dVair*sinarg
          ugeoDY(1,1,k) = uairLS(k)
          pktaLS(    k) = pkta_0(k) + dTair*sinarg*3.72
        ENDDO
        END IF                                                         ! CTR
        IF       (mod(iterun,n_iter).eq.0)                        THEN ! CTR
                      ipr2nc=ipr2nc  +  1

C +       ***********
          call Out2nc
C +       ***********

        END IF                                                         ! CTR
      END IF                                                           ! CTR


C +--Relax to  Large Scale Forcing
C +  -----------------------------

      IF (relxLS)                                                 THEN ! CTR
        DO k=1,mz
        DO j=1,my
        DO i=1,mx
          uairDY(i,j,k)=uairDY(i,j,k)+dt*(uairLS(k)-uairDY(i,j,k))/3.6e3
          vairDY(i,j,k)=vairDY(i,j,k)+dt*(vairLS(k)-vairDY(i,j,k))/3.6e3
          pktaDY(i,j,k)=pktaDY(i,j,k)+dt*(pktaLS(k)-pktaDY(i,j,k))/3.6e3
            qvDY(i,j,k)=  qvDY(i,j,k)+dt*(  qvLS(k)-  qvDY(i,j,k))/3.6e3
        END DO
        END DO
        END DO
      END IF                                                           ! CTR


C +--Blowing Snow is advected away
C +  -----------------------------

        DO j=1,my
        DO i=1,mx
        DO k=1,mz
            qsHY(i,j,k)=  0.00
        END DO
          SnowEq       =  0.00
        DO k=1,nssSNo(i,j,1)
          SnowEq       =  SnowEq + rosSNo(i,j,1,k) *dzsSNo(i,j,1,k)
        END DO
          SnowEq       =  SnowEq / ro_Wat
        END DO
        END DO

      write(6,6000) itexpe,jdarGE,labmGE(mmarGE),jhurGE,minuGE,jsecGE,
     .              snowHY(1,1)  ,SnowEq        ,WEq_SN(1,1,1),
     .              rosSNo(1,1,1 ,nssSNo(1,1,1)),SaltSN(1,1,1),
     .              albeSL(1,1)
 6000 format(i8,i6,'-',a3,i6,'h',i2,':',i2,6f12.6)

c #VR IF (itexpe.GE.30.AND.itexpe.LE.30) call DUMP

      return
      end 


      subroutine Out2nc

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                              3-10-2002  MAR |
C |   SubRoutine Out2nc is used to write the main Model Variables          |
C |                                      on  a NetCDF file                 |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT: ipr2nc: Current time step    number                           |
C |   ^^^^^^         (starting from ipr2nc=1, which => new file creation)  |
C |          npr2nc: Total  'time slices' number (max value of ipr2nc)     |
C |                                                                        |
C |   OUTPUT: NetCDF File adapted to IDL Graphic Software                  |
C |   ^^^^^^                                                               |
C |                                                                        |
C |   CAUTION: 1) This Routine requires the usual NetCDF library,          |
C |   ^^^^^^^^    and the complementary access library  'libUN.a'          |
C |                                                                        |
C +------------------------------------------------------------------------+


      implicit none


C +--General Variables
C +  =================

      include 'MARCTR.inc'
      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_DY.inc'
      include 'MAR_TE.inc'
      include 'MAR_TU.inc'

      include 'MAR_HY.inc'
      include 'MAR_RA.inc'
c #TC include 'MAR_TC.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MARsSN.inc'

      include 'MAR_WK.inc'

      include 'MAR_IO.inc'

      integer           ipr2nc,npr2nc,n_iter
      common/Out2nc_arg/ipr2nc,npr2nc,n_iter


C +--Local   Variables
C +  =================

      integer    Lfnam,     Ltit,     Luni,     Lnam,     Llnam
      PARAMETER (Lfnam= 40, Ltit= 90, Luni= 31, Lnam= 13, Llnam=50)
C +...Length of char strings 

      CHARACTER*(Lfnam)  fnamNC
      common/Out2nc_loc/ fnamNC
C +...                   fnamNC: To retain file name.

      real               snow0(mx,my),bsnow0(mx,my),bsnow1(mx,my)
      common/Out2rr_loc/ snow0       ,bsnow0       ,bsnow1
C +...                   snow0 : Snow Precipitated over Previous Time Interval
C +                     bsnow0 : Snow Budget       at   Previous Time Interval
C +                     bsnow1 : Snow Budget       at   Current  Time Interval

      integer    NdimNC
      PARAMETER (NdimNC = 5)
C +...Number of defined spatial dimensions (exact)

      integer    MXdim
      PARAMETER (MXdim  = 600000)
C +...Maximum Number of all dims: recorded Time Steps
C +   and also maximum of spatial grid points for each direction. 

      integer    MX_var
      PARAMETER (MX_var = 80)
C +...Maximum Number of Variables 

      integer    NattNC
      PARAMETER (NattNC = 2)
C +...Number of REAL attributes given to all variables

      INTEGER           RCODE

      integer           jourNC(MXdim)
      integer           moisNC(MXdim)
      real              yearNC(MXdim)
      real              dateNC(MXdim)
      real              timeNC(MXdim)
      real              VALdim(MXdim,0:NdimNC)
      integer           nDFdim(      0:NdimNC)
      common/c_nDFdim/  nDFdim
      integer           NvatNC(NattNC)
      CHARACTER*(Lnam)  NAMdim(      0:NdimNC)
      CHARACTER*(Luni)  UNIdim(      0:NdimNC)
      CHARACTER*(Lnam)  SdimNC(4,MX_var)       
      CHARACTER*(Luni)  unitNC(MX_var)
      CHARACTER*(Lnam)  nameNC(MX_var)
      CHARACTER*(Llnam) lnamNC(MX_var)
      CHARACTER*(Ltit ) tit_NC
      CHARACTER*(Lnam)  NAMrat(NattNC)
c #TC CHARACTER*9   labelc
      CHARACTER*120 tmpINP

      integer   mz_SBL
      parameter(mz_SBL=nsno)
      real      ua_SBL(mx,my,mz_SBL),va_SBL(mx,my,mz_SBL)
      real      Ta_SBL(mx,my,mz_SBL),zz_SBL(mx,my,mz_SBL)
      real      HS_SBL(mx,my,mz_SBL),HL_SBL(mx,my,mz_SBL)
      real      HS_10m(mx,my)       ,HL_10m(mx,my)       ,fac10m

      integer   isn   ,n
      real      T_SNOW(mx,my,nsno),roSNOW(mx,my,nsno)
      real      Z_SNOW(mx,my,nsno),AGSNOW(mx,my,nsno)
      real      G1SNOW(mx,my,nsno),G2SNOW(mx,my,nsno),GGSNOW(mx,my,nsno)

      integer   n1000 ,n100a ,n100  ,n10_a ,n10   ,n1    
      integer   m10   ,       jd10  ,jd1
      integer   MMXstp,it    ,mois  ,mill  ,iu
      integer   itotNC,NtotNC,ID__nc
      real      starta
      REAL      starti,DayLen,optwa ,optia ,rhodzk

      REAL      gsMAXg


C +--DATA
C +  ====

      data      gsMAXg/15./


C +--NetCDF File Initialization
C +  ==========================

      IF (ipr2nc.eq.1) THEN

          n1000 = 1 +     iyrrGE/1000
          n100a =     mod(iyrrGE,1000)
          n100  = 1 +     n100a /100
          n10_a =     mod(n100a ,100)
          n10   = 1 +     n10_a /10
          n1    = 1 + mod(n10_a ,10)
          m10   = 1 +     mmarGE/10
          m1    = 1 + mod(mmarGE,10)
          jd10  = 1 +     jdarGE/10
          jd1   = 1 + mod(jdarGE,10)


C +--Output File Label
C +  -----------------

        fnamNC = 'ANI.'
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)
     .         // '.' // explIO
     .         // '.nc    '


C +--Output Title
C +  ------------

        tit_NC = 'MAR'
     .         // ' - Exp: ' // explIO
     .         // ' - '
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)


C +--Create File / Write Constants
C +  -----------------------------
        MMXstp = MXdim
C +...  To check array bounds... silently

C +--Time Variable (hour)
C +  ~~~~~~~~~~~~~~~~~~~~

C +...  To define a NetCDF dimension (size, name, unit):
c _UL   nDFdim(0)= npr2nc
        nDFdim(0)= 0
        NAMdim(0)= 'time'
        UNIdim(0)= 'HOURS since 1901-01-15 00:00:00'

C +...  Check temporary arrays: large enough ?
        IF (npr2nc.gt.MMXstp)
     &  STOP '*** Out2nc - ERROR : MXdim to low ***'

              starti = jhurGE + minuGE/60.d0 + jsecGE/3600.d0
C +...        starti : Starting Time (= current time in the day)

              starta = (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
     .               + (minuGE *60 +jsecGE      )/3600.!

        DO it = 1,npr2nc
              timeNC(it)   = starti + (it-1) * n_iter  *dt / 3600.
C +...                                         n_iter: #iter between output

              VALdim(it,0) = starta + (it-1) * n_iter  *dt / 3600.
C +...        VALdim(  ,0) : values of the dimension # 0 (time) 

C +--Time Variable (date)
C +  ~~~~~~~~~~~~~~~~~~~~
              dateNC(it) =          timeNC(it)
              jourNC(it) = jdarGE + timeNC(it) / 24.d0
        END DO
                  mois       =  mmarGE
                  mill       =  iyrrGE
        DO it = 1,npr2nc
          IF     (jourNC(it).gt.njmoGE(mois))                     THEN ! CTR
            DO iu=it,npr2nc
                  jourNC(iu) =  jourNC(iu) - njmoGE(mois)
            END DO
                  mois       =  mois + 1
              IF (mois.gt.12)                                     THEN ! CTR
                  mois       =         1
                  mill       =  mill + 1
              END IF                                                   ! CTR
          END IF                                                       ! CTR
                  moisNC(it) =  mois
                  yearNC(it) =  mill

          IF     (dateNC(it).gt.24.d0-epsi)                       THEN ! CTR
                  DayLen     =  24.d0
            DO iu=it,npr2nc
                  dateNC(iu) = mod(dateNC(iu),DayLen)
            END DO
          END IF                                                       ! CTR
        END DO

        DO it = 1,npr2nc
              dateNC(it) =  dateNC(it)
     .             + 1.d+2 *jourNC(it)
     .             + 1.d+4 *moisNC(it)
        END DO


C +   Define horizontal spatial dimensions :    
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +...  Check temporary arrays: large enough ?
        IF (    mx .gt.MMXstp.or.my.gt.MMXstp
     &      .or.mzz.gt.MMXstp.or.mw.gt.MMXstp)
     &    STOP '*** Out2nc - ERROR : MXdim to low ***'

C +...To define NetCDF dimensions (size, name, unit):

        DO i = 1, mx
          VALdim(i,1) = xxkm(i)
        END DO
          nDFdim(1)   = mx
          NAMdim(1)   = 'x'
          UNIdim(1)   = 'km'

        DO j = 1, my
          VALdim(j,2) = yykm(j)
        END DO
          nDFdim(2)   = my
          NAMdim(2)   = 'y'
          UNIdim(2)   = 'km'

        do k = 1, mz_SBL
          VALdim(k,3) = sigma(mzz-k)
        enddo
          nDFdim(3)   =  mz_SBL
          NAMdim(3)   = 'level'
          UNIdim(3)   = '[sigma]'
C +...    For levels k

        do k = 1, mz_SBL
          VALdim(k,4) = sigmid(mzz-k+1)
        enddo
          nDFdim(4)   =  mz_SBL
          NAMdim(4)   = 'level2'
          UNIdim(4)   = '[sigma]'
C +...    For levels k+1/2

        do k = 1, mw
          VALdim(k,5) = k 
        enddo
          nDFdim(5)   = mw
          NAMdim(5)   = 'sector'
          UNIdim(5)   = '[index]'
C +...    For Surface Sectors 

C +--Variable's Choice (Table MARgou.dat)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        OPEN(unit=10,status='unknown',file='MARgou.dat')

        itotNC = 0
 980    CONTINUE
          READ (10,'(A120)',end=990) tmpINP
          IF (tmpINP(1:4).eq.'    ')                                THEN
            itotNC = itotNC + 1
            READ (tmpINP,'(4x,5A9,A12,A50)')
     &          nameNC(itotNC)  ,SdimNC(1,itotNC),SdimNC(2,itotNC),
     &          SdimNC(3,itotNC),SdimNC(4,itotNC),
     &          unitNC(itotNC)  ,lnamNC(itotNC)
C +...          nameNC: Name
C +             SdimNC: Names of Selected Dimensions (max.4/variable) 
C +             unitNC: Units
C +             lnamNC: Long_name, a description of the variable

c #TC       IF (nameNC(itotNC).eq.'qxTC     '.and.nkWri.ge.1)       THEN
c #TC           nameNC(itotNC) =  namTC(1)
c #TC        IF                                  (nkWri.gt.1)       THEN
c #TC           itot = itotNC
c #TC         DO  n=2,nkWri
c #TC           itot = itot    +  1
c #TC           nameNC(itot)   =  namTC(n)
c #TC          DO m=1,4
c #TC           SdimNC(m,itot) = SdimNC(m,itotNC)
c #TC          END DO
c #TC           unitNC(itot)   = unitNC(itotNC)
c #TC           lnamNC(itot)   = lnamNC(itotNC)
c #TC         END DO
c #TC           itotNC         = itot
c #TC        ENDIF
c #TC       ENDIF

          ENDIF
        GOTO 980
 990    CONTINUE

        CLOSE(unit=10)

        NtotNC = itotNC 
C +...  NtotNC : Total number of variables writen in NetCDF file.

C +--List of NetCDF attributes given to all variables:
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +...  The "actual_range" is the (min,max)
C +     of all data for each variable:
        NAMrat(1) = 'actual_range'
        NvatNC(1) = 2

C +...  The "[var]_range" is NOT of attribute type,
C +     it is a true variable containing the (min,max) for
C +     each level, for 4D (space+time) variables only
C +     (automatic handling by UN library;
C +      must be the LAST attribute)
        NAMrat(NattNC) = '[var]_range'
        NvatNC(NattNC) = 2

C +--Automatic Generation of the NetCDF File Structure
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +     **************
        CALL UNscreate (fnamNC, tit_NC,
     &                  NdimNC, nDFdim, MXdim , NAMdim, UNIdim, VALdim,
     &                  MX_var, NtotNC, nameNC, SdimNC, unitNC, lnamNC,
     &                  NattNC, NAMrat, NvatNC,
     &                  ID__nc) 
C +     **************


C +--Write Time - Constants
C +  ~~~~~~~~~~~~~~~~~~~~~~
        DO j=1,my
        DO i=1,mx
          Wkxy1(i,j) =  GElonh(i,j) * 15.d0
C +...    Conversion: Hour->degrees

          WKxy2(i,j) =  GElatr(i,j) / degrad
C +...    Conversion: rad ->degrees

          WKxy3(i,j) =  isolSL(i,j)
C +...    Conversion to REAL type (integer not allowed)

        END DO
        END DO


C +       ************
          CALL UNwrite (ID__nc, 'lon   ', 1  , mx    , my, 1 , Wkxy1)
          CALL UNwrite (ID__nc, 'lat   ', 1  , mx    , my, 1 , Wkxy2)
          CALL UNwrite (ID__nc, 'sh    ', 1  , mx    , my, 1 , sh)
          CALL UNwrite (ID__nc, 'isol  ', 1  , mx    , my, 1 , Wkxy3)
C +       ************

C +--Re-Open file if already created.
C +  ================================


      ELSE

C +     ************
        CALL UNwopen (fnamNC,ID__nc)
C +     ************

      END IF


C +--Write Time-dependent variables:
C +  ===============================


C +--UNLIMITED Time Dimension
C +  ------------------------

      IF (nDFdim(0).eq.0)                         THEN !
              starta = (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
     .               + (minuGE *60 +jsecGE      )/3600.!

C +     ************
        CALL UNwrite (ID__nc, 'time   ',ipr2nc, 1, 1, 1, starta)
C +     ************

      END IF


C +     ************
        CALL UNwrite (ID__nc, 'date   ',ipr2nc, 1, 1, 1, dateNC(ipr2nc))
        CALL UNwrite (ID__nc, 'year   ',ipr2nc, 1, 1, 1, yearNC(ipr2nc))
C +     ************


C +--Dynamics
C +  --------

      DO j=1,my
      DO i=1,mx
        ua_SBL(i,j,1) = 0.
        va_SBL(i,j,1) = 0.
        Ta_SBL(i,j,1) = tairSL(i,j)
        zz_SBL(i,j,1) = 0.
       DO k=2,mz_SBL
        ua_SBL(i,j,k) = uairDY(i,j,mzz-k+1)
        va_SBL(i,j,k) = vairDY(i,j,mzz-k+1)
        Ta_SBL(i,j,k) = tairDY(i,j,mzz-k+1)
        zz_SBL(i,j,k) = gplvDY(i,j,mzz-k+1)*grvinv-sh(i,j)
       END DO


C +--Turbulent Heat Fluxes
C +  ---------------------

        HS_SBL(i,j,1) = -SLuts(i,j)
        HL_SBL(i,j,1) = -SLuqs(i,j)
       DO k=2,mz_SBL-1
        HS_SBL(i,j,k) = -TUkvh(i,j,mzz-k)
     .                *(pktaDY(i,j,mzz-k)-pktaDY(i,j,mzz-k+1))
     .                /(zz_SBL(i,j,k+1)  -zz_SBL(i,j,k)      )
        HL_SBL(i,j,k) = -TUkvh(i,j,mzz-k)
     .                  *(qvDY(i,j,mzz-k)-qvDY(i,j,mzz-k+1))
     .                /(zz_SBL(i,j,k+1)  -zz_SBL(i,j,k)      )
       END DO
       DO k=1,mz_SBL-1
        HS_SBL(i,j,k) = HS_SBL(i,j,k)*cp*pcap
        HL_SBL(i,j,k) = HL_SBL(i,j,k)*Ls_H2O
       END DO

       k = 1
 1000  CONTINUE
       k = k + 1
       IF (zz_SBL(i,j,k).LT.10.)                              GO TO 1000
       fac10m    = (10.-zz_SBL(i,j,k-1))/(zz_SBL(i,j,k)-zz_SBL(i,j,k-1))
       HS_10m(i,j) =    HS_SBL(i,j,k-1)
     .             +    fac10m          *(HS_SBL(i,j,k)-HS_SBL(i,j,k-1))
       HL_10m(i,j) =    HL_SBL(i,j,k-1)
     .             +    fac10m          *(HL_SBL(i,j,k)-HL_SBL(i,j,k-1))


C +--Precipitation
C +  -------------

         WKxy4(i,j)   = snowHY(i,j) - snow0(i,j)
         snow0(i,j)   = snowHY(i,j)


C +--Total   Snow Budget
C +  -------------------

        bsnow1(i,j)   =              snohSN(i,j,1)
     .                + (snowHY(i,j)-sno0HY(i,j  ))  *1.0e+3
       IF        (nssSNo(i,j,1).gt.   0      )                      THEN
       DO k=max(0,nssSNo(i,j,1)),nssSNo(i,j,1)
        bsnow1(i,j)   =  bsnow1(i,j)+rosSNo(i,j,1,k) *dzsSNo(i,j,1,k)
       END DO
       END IF
       DO k=1,mz
        bsnow1(i,j)   = bsnow1(i,j)+ rolvDY(i,j,  k) *  qsHY(i,j,  k)
     .                 *grvinv     *(gpmiDY(i,j,  k) -gpmiDY(i,j,  k+1))
     .                 *1.0e+3
       END DO
       IF (iterun.gt.0)                                             THEN
         WKxy1(i,j)   =(bsnow1(i,j)- bsnow0(i,j)    )*86400./dt
       ELSE
         WKxy1(i,j)   = 0.
       END IF
        bsnow0(i,j)   = bsnow1(i,j)


C +--Blowing Snow Transport
C +  ----------------------

         WKxy2(i,j)   = 0.

       k=0
 5011  CONTINUE
       k=k+1
       IF (                  k         .gt. mz )      GO TO 5010
       IF (grvinv*gplvDY(i,j,k)-sh(i,j).lt.100.)            THEN
         WKxy2(i,j)    = WKxy2(i,j)+ssvSL(i,j,k)*qsHY(i,j,k)
     .                             *pstDY(i,j)*dsigm1(    k)*1.e3*grvinv
       END IF
       GO TO 5011
 5010  CONTINUE

         WKxy5(i,j)   = SaltSN(i,j,1)
         WKxy6(i,j)   = SLuusl(i,j,1)
         WKxy7(i,j)   = SLutsl(i,j,1)
         WKxy8(i,j)   = SLussl(i,j,1)
         WKxy9(i,j)   = agsSNo(i,j,1,nssSNo(i,j,1))

         DO isn=1,nssSNo(i,j,1)
           Z_SNOW(i,j,isn) = Z_SNOW(i,j,  isn-1) 
     .                     + dzsSNo(i,j,1,isn)
           T_SNOW(i,j,isn) = tisSNo(i,j,1,isn)
           roSNOW(i,j,isn) = rosSNo(i,j,1,isn)
           AGSNOW(i,j,isn) = agsSNo(i,j,1,isn)
           G1SNOW(i,j,isn) = G1sSNo(i,j,1,isn)
           IF               (G1sSNo(i,j,1,isn).gt.0.)               THEN
           G2SNOW(i,j,isn) = G2sSNo(i,j,1,isn)
           ELSE
           G2SNOW(i,j,isn) = 0.
           END IF
           IF (g1sSNo(i,j,n,isn).LT.0.)                             THEN
            IF(g2sSNo(i,j,n,isn).LT.50.)                            THEN
               GGSNOW(i,j,isn) = 0.2-0.002*g1sSNo(i,j,1,isn)
            ELSE
               GGSNOW(i,j,isn) = 0.6+0.002*g1sSNo(i,j,1,isn)
            END IF
           ELSE
            IF(g1sSNo(i,j,n,isn).LT.50.)                            THEN
               GGSNOW(i,j,isn) = 
     .                0.2-0.20*(min(gsMAXg,g2sSNo(i,j,1,isn))-3.0)/12.
            ELSE
               GGSNOW(i,j,isn) = 
     .                0.8+0.20*(min(gsMAXg,g2sSNo(i,j,1,isn))-3.0)/12.
            END IF
           END IF
         END DO
         DO isn=nssSNo(i,j,1),nsno
           Z_SNOW(i,j,isn) = Z_SNOW(i,j,  nssSNo(i,j,1))
           T_SNOW(i,j,isn) = tisSNo(i,j,1,nssSNo(i,j,1))
           roSNOW(i,j,isn) = rosSNo(i,j,1,nssSNo(i,j,1))
           AGSNOW(i,j,isn) = AGSNOW(i,j,  nssSNo(i,j,1))
           G1SNOW(i,j,isn) = G1SNOW(i,j,  nssSNo(i,j,1))
           G2SNOW(i,j,isn) = G2SNOW(i,j,  nssSNo(i,j,1))
           GGSNOW(i,j,isn) = GGSNOW(i,j,  nssSNo(i,j,1))
         END DO
      END DO
      END DO
c #WR write(6,5012) jdarGE,labmGE(mmarGE),iyrrGE,jhurGE,minuGE,jsecGE,
c #WR.              WKxy1(1,1)
 5012 format(i6,'-',a3,'-',i4,i3,'h',i2,'.',i2,f12.3,' mmWE')


C +--Cloud Optical Depth
C +  -------------------

      do j=1,my
      do i=1,mx
           optwa  = zero
C +...     optwa  : liquid water path (kg/m2) (droplets)

           optia  = zero
c +...     optia  : liquid water path (kg/m2) (crystals)

         DO k = mzabso+1,mz
           rhodzk = (pstDY(i,j)*sigma(k)+ptopDY)
     .            / (ra*tairDY(i,j,k)*(1.+.608*qvDY(i,j,k)))
     .            * (gpmiDY(i,j,k)-gpmiDY(i,j,k+1))
C +...     rhodzk : (rho / 1000) * (dz * gravit)

           optwa  = optwa + rhodzk * qwHY(i,j,k) 
           optia  = optia + rhodzk * qiHY(i,j,k)

         END DO

           WKxy3(i,j)  = 1.5 * ( optwa / 20.d-6 
     .                         + optia / 40.d-6 ) *grvinv

      enddo
      enddo

C +   ************
      CALL UNwrite (ID__nc, 'pstar  ', ipr2nc, mx, my, 1      , pstDY )
      CALL UNwrite (ID__nc, 'ua_SBL ', ipr2nc, mx, my, mz_SBL , ua_SBL)
      CALL UNwrite (ID__nc, 'va_SBL ', ipr2nc, mx, my, mz_SBL , va_SBL)
      CALL UNwrite (ID__nc, 'Ta_SBL ', ipr2nc, mx, my, mz_SBL , Ta_SBL)
      CALL UNwrite (ID__nc, 'HS_SBL ', ipr2nc, mx, my, mz_SBL , HS_SBL)
      CALL UNwrite (ID__nc, 'HS_10m ', ipr2nc, mx, my, 1      , HS_10m)
      CALL UNwrite (ID__nc, 'HL_SBL ', ipr2nc, mx, my, mz_SBL , HL_SBL)
      CALL UNwrite (ID__nc, 'HL_10m ', ipr2nc, mx, my, 1      , HL_10m)
      CALL UNwrite (ID__nc, 'zz_SBL ', ipr2nc, mx, my, mz_SBL , zz_SBL)
      CALL UNwrite (ID__nc, 'SnowTB ', ipr2nc, mx, my, 1      , WKxy1 )
      CALL UNwrite (ID__nc, 'SnowTR ', ipr2nc, mx, my, 1      , WKxy2 )
      CALL UNwrite (ID__nc, 'RadOLR ', ipr2nc, mx, my, 1      , RAdOLR)
      CALL UNwrite (ID__nc, 'OptDep ', ipr2nc, mx, my, 1      , WKxy3 )
      CALL UNwrite (ID__nc, 'Snow   ', ipr2nc, mx, my, 1      , WKxy4 )
      CALL UNwrite (ID__nc, 'ue_star', ipr2nc, mx, my, 1      , WKxy5 )
      CALL UNwrite (ID__nc, 'uu_star', ipr2nc, mx, my, 1      , WKxy6 )
      CALL UNwrite (ID__nc, 'ut_star', ipr2nc, mx, my, 1      , WKxy7 )
      CALL UNwrite (ID__nc, 'us_star', ipr2nc, mx, my, 1      , WKxy8 )
      CALL UNwrite (ID__nc, 'AgeSnow', ipr2nc, mx, my, 1      , WKxy9 )
      CALL UNwrite (ID__nc, 'Z_SNOW ', ipr2nc, mx, my, nsno   , Z_SNOW)
      CALL UNwrite (ID__nc, 'T_SNOW ', ipr2nc, mx, my, nsno   , T_SNOW)
      CALL UNwrite (ID__nc, 'roSNOW ', ipr2nc, mx, my, nsno   , roSNOW)
      CALL UNwrite (ID__nc, 'AGSNOW ', ipr2nc, mx, my, nsno   , AGSNOW)
      CALL UNwrite (ID__nc, 'GGSNOW ', ipr2nc, mx, my, nsno   , GGSNOW)
      CALL UNwrite (ID__nc, 'G1SNOW ', ipr2nc, mx, my, nsno   , G1SNOW)
      CALL UNwrite (ID__nc, 'G2SNOW ', ipr2nc, mx, my, nsno   , G2SNOW)
C +   ************


C +--That 's all, folks: NetCDF File Closure
C +  =======================================

C +   ***********
      CALL NCCLOS (ID__nc,RCODE)
C +   ***********


C +--Work Arrays Reset
C +  =================

      do j=1,my
      do i=1,mx
        WKxy1(i,j)   =0.0
        WKxy2(i,j)   =0.0
        WKxy3(i,j)   =0.0
        WKxy4(i,j)   =0.0
        WKxy5(i,j)   =0.0
        WKxy6(i,j)   =0.0
        WKxy7(i,j)   =0.0
        WKxy8(i,j)   =0.0
        WKxy9(i,j)   =0.0
      enddo
      enddo

      return
      end
