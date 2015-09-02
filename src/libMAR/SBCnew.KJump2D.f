      subroutine SBCnew

C +------------------------------------------------------------------------+
C | MAR INPUT      SURFACE BOUNDARY LAYER                  17-06-2004  MAR |
C |   SubRoutine SBCnew includes Snow Profile Output                       |
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

      include 'MAR_DY.inc'
      include 'MAR_LB.inc'

      include 'MAR_SV.inc'
      include 'MAR_TV.inc'
      include 'MAR_SL.inc'
      include 'MARsSN.inc'


      logical           WRIlog
      common/SBCnew_log/WRIlog


C +--Local   Variables
C +  =================

      integer           iprANI,nprANI,nt_ANI
      common/ANI_nc_arg/iprANI,nprANI,nt_ANI

      integer           iprAWS,nprAWS,nt_AWS
      common/AWS_nc_arg/iprAWS,nprAWS,nt_AWS

      integer     isl   ,misl_2,nisl_2,n     ,isn,isnh

      real        dz_dSN(-nsno:1)     ,SnowWE(nsno+1)
      real        Dendri,Spheri,GrainS,Densit,Age

      REAL        ANI_dt
      REAL        AWS_dt


C +--DATA
C +  ====

      data        ANI_dt/ 360.0/  ! Sampling Time Interval on ANI NetCDF File
      data        AWS_dt/ 360.0/  ! Sampling Time Interval on AWS NetCDF File


C +--Snow Model Initialisation
C +  =========================

      IF (iterun.eq.       0)                                       THEN


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
     .                    +  misl_2  * 0.003) * 10**(nisl_2)) * 2.
            dz_dSN(isl)=max(dz_dSN(isl),0.01)
            dz_dSN(isl)=min(dz_dSN(isl),0.20)
          END DO


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
                  SnowWE(nsno+1)     =   0.                  !
              DO isn = nsno,1,-1                             !
                  dzsSNo(i,j,n,isn)  =  dz_dSN(isn-nsno)     !            [m]
                  rosSNo(i,j,n,isn)  =  Densit               !        [kg/m3]
                  wasSNo(i,j,n,isn)  =   0.                  !        [kg/kg]
                  tisSNo(i,j,n,isn)  =  tairSL(i,j)          !            [K]
                IF (Dendri.lt.0.)     THEN                   !
                  g1sSNo(i,j,n,isn)  =  Dendri               ! [-]        [-]
                  g2sSNo(i,j,n,isn)  =  Spheri               ! [-]        [-]
                  nhsSNo(i,j,n,isn)  =   0                   !            [-]
                ELSE                                         !
                  g1sSNo(i,j,n,isn)  =  Spheri               ! [-]        [-]
                  g2sSNo(i,j,n,isn)  =  GrainS               ! [-] [0.0001 m]
                  SnowWE(      isn)  =  SnowWE(      isn+1)  !
     .           +dzsSNo(i,j,n,isn)  *  rosSNo(i,j,n,isn)    !
                  IF                   (SnowWE(      isn) .LE. 05.) THEN
                  nhsSNo(i,j,n,isn)  =   0                   ! (MAY BE 2) [-]
                  ELSE                                       !
                  nhsSNo(i,j,n,isn)  =   2                   !  MAY BE 2! [-]
              !   Glazed Snow // Snow previously in Contact with Water  (crude)
              !  (Glazed Snow is generated by Dissipation of Kinetic Energy
              !                                           of Blown Snow Grains)
                  ENDIF
                END IF                                       !
                  agsSNo(i,j,n,isn)  =  Age                  !          [day]
              END DO
              IF (.NOT.WRIlog)                                      THEN
                WRIlog = .true.
                write(6,6000)
 6000           format(/,'SNOW MODEL INITIALIZATION:'
     .                 /,'==========================')
                write(6,6020) i,j,n,isolSL(i,j)
     .                      ,jhurGE,labmGE(mmarGE),minuGE,jsecGE
 6020           format(/,i3,4i8,a8,2i8,
     .                 /,'----+-------+-------+-------+-------+',
     .                                '-------+-------+-------+',
     .                 /,'isn |     Z |    dz |    ro |    Ts |',
     .                                '    G1 |    G2 |    nh |',
     .                 /,'    | [mmWE]|   [m] |[kg/m3]|   [K] |',
     .                                '   [-] |   [-] |   [-] |',
     .                 /,'----+-------+-------+-------+-------+',
     .                                '-------+-------+-------+')
                DO isn=nssSNo(i,j,n),1,-1
                  write(6,6021) isn,SnowWE(isn)
     .                       ,dzsSNo(i,j,n,isn),rosSNo(i,j,n,isn)
     .                       ,tisSNo(i,j,n,isn)
     .                       ,G1sSNo(i,j,n,isn),G2sSNo(i,j,n,isn)
     .                       ,nhsSNo(i,j,n,isn)
 6021             format(i3,' |',f6.0,' |',f6.3,' |',f6.0,' |'
     .                        ,3(f6.1,' |'), i6,' |')
                END DO
                  write(6,6022)
 6022             format(  '----+-------+-------+-------+-------+',
     .                                  '-------+-------+-------+',/)
              END IF
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


C +--Upper Limit of the non erodible Snow (nhsSNo .GT. 1)
C +  ------------------------------------

           isnh =   0
        DO isn=  nsno,1,-1
           isnh =isnh + isn* min(nhsSNo(i,j,n,isn)-1,1)*max(0,1-isnh)
        ENDDO
           zWEcSN(i,j,n) =     0. 
        DO isn=1,nsno
           zWEcSN(i,j,n) =       zWEcSN(i,j,n) 
     .                   +       dzsSNo(i,j,n,isn) *rosSNo(i,j,n,isn)
     .                     * max(0,min(1,isnh+1-isn))
        END DO


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


C +--OPEN  NetCDF FILE
C +  =================

C +--Animation NetCDF File
C +  ---------------------

              nt_ANI = ANI_dt / dt            ! Number of Time steps between IO

              iprANI =                   1
              nprANI = nterun / nt_ANI + 1    ! Number of                    IO

C +       ***********
          call ANI_nc
C +       ***********


C +--AWS       NetCDF File
C +  ---------------------

              nt_AWS = AWS_dt / dt            ! Number of Time steps between IO
          IF (nt_AWS.EQ.0)                                          THEN
              nt_AWS = 1
              AWS_dt = dt
          END IF

              iprAWS =                   1
              nprAWS = nterun / nt_AWS + 1    ! Number of                    IO

C +       ***********
          call AWSloc
          call AWS_nc
C +       ***********

      ELSE


C +--ANImation NetCDF File
C +  ---------------------

        IF       (mod(iterun,nt_ANI).eq.0)                          THEN
                      iprANI=iprANI  +  1

C +       ***********
          call ANI_nc
C +       ***********

        END IF


C +--AWS       NetCDF File
C +  ---------------------

        IF       (mod(iterun,nt_AWS).eq.0)                        THEN
                      iprAWS=iprAWS  +  1

C +       ***********
          call AWS_nc
C +       ***********

        END IF


      END IF


C +--ASCII OUTPUT
C +  ============

                            i=imez+1
                            j=     1
                            n=     1
                  SnowWE(nsno+1)     =   0.                  !
              DO isn = nsno,1,-1                             !
                  SnowWE(      isn)  =  SnowWE(      isn+1)  !
     .           +dzsSNo(i,j,n,isn)  *  rosSNo(i,j,n,isn)    !
              END DO
                write(6,6020) i,j,n,isolSL(i,j)
     .                      ,jhurGE,labmGE(mmarGE),minuGE,jsecGE
              DO isn=nssSNo(i,j,n),1,-1
                write(6,6021) isn,SnowWE(isn)
     .                     ,dzsSNo(i,j,n,isn),rosSNo(i,j,n,isn)
     .                     ,tisSNo(i,j,n,isn)
     .                     ,G1sSNo(i,j,n,isn),G2sSNo(i,j,n,isn)
     .                     ,nhsSNo(i,j,n,isn)
              END DO
                write(6,6022)
                write(6,6009)       SLuusl(i,j,n),     SLutsl(i,j,n)
     .                             ,SaltSN(i,j,n),1.e3* SL_z0(i,j,n)
     .                                           ,1.e3* SL_r0(i,j,n)
 6009           format('  u*  =',f11.3,' m/s   uT* =',f6.3,' K.m/s'
     .               /,'  u*t =',f11.3,' m/s   z0m =',f6.3
     .                                 ,' mm   z0h =',f6.3,' mm')
                write(6,6010)
 6010           format(1x)


C +--GRAPHICS ASCII FILE
C +  ===================

C +       *************
          call Sastrugi
C +       *************

      return
      end 


      subroutine Sastrugi

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                              3-10-2002  MAR |
C |   SubRoutine Sastrugi is used to write the main Model Variables        |
C |                                                                        |
C |   FORMAT: ASCII                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARCTR.inc'
      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_SV.inc'
      include 'MAR_SN.inc'
      include 'MARsSN.inc'


C +--Local   Variables
C +  =================

      real          distxx(mx),depth0(mx),depths(mx)
      common/SBCout/distxx    ,depth0    ,depths

      REAL        facx  ,facy  ,cste  ,countv,depthv
      integer     isgn  ,mark  ,month
      integer     n     ,jda_10,jda_01,jha_10,jha_01,number,ncurv,itimer

      character * 1 filnex
      character * 8  forme, formi,  line
      character *39 titxxx,tityyy,titbla
      character *13 filn_h


C +--DATA
C +  ====

      data facx/ 1.d00/
      data facy/ 1.d00/
      data cste/ 0.d00/
      data isgn/ 1    /
      data mark/-3    /
      data month/1    /
      data titxxx/'             Distance [km]             $'/
      data tityyy/'         Snow Thickness   [cm]         $'/
      data titbla/'                                       $'/


C +--Assignation
C +  ===========

      if (iterun.eq.0     ) then
        open (unit=31,status='unknown',file='MAR_SN.out')
        rewind     31
      end if

      itimer= iterun*dt
      j     = 1

      if (mod(itimer,3600).eq.0) then
        countv    = 0.0d+0
        depthv    = 0.0d+0
        do i=imez+1,mx1
        depths(i) = 0.0d+0
        IF (VSISVAT)                                              THEN ! CTR
          do n=     1,nsno
          depths(i) = depths(i) + dzsSNo(i,j,1,n)
          end do
        ELSE
          do n=     1,nsno
          depths(i) = depths(i) + dzSNow(i,j,n)
          end do
        END IF
        countv    = 1.0d+0    + countv
        depthv    = depthv    + depths(i)
        end do
        depthv    = depthv    / countv

       if    (itimer      .eq.0) then
        do i=imez,mx
        distxx(i) = (i - imez) *dx *1.d-3
        depth0(i) =  depths(i)
        end do
       end if

        do i=imez,mx
        depths(i) = depths(i) - depthv
        end do

        write(31,311) jdaMAR,jhaMAR,(1.d2*depths(i),i=imez+1,mx1)
 311    format(2i3,20f6.1,/,(6x,20f6.1))


C +--GRAPHIC FILE
C +  ============

       filn_h            ='Sastrugi.outp'
       jda_10            =        jdaMAR/10
       jda_01            =    mod(jdaMAR,10)
       filn_h(10:10)     = labnum(jda_10+1)
       filn_h(11:11)     = labnum(jda_01+1)
       jha_10            =        jhaMAR/10
       jha_01            =    mod(jhaMAR,10)
       filn_h(12:12)     = labnum(jha_10+1)
       filn_h(13:13)     = labnum(jha_01+1)
       open (unit=32,file= filn_h           ,status='unknown')
       rewind     32

       write(32,3201) 
 3201  format(' Transect  Dumont-d`Urville/Dome C       0  ',
C +...         12345678901234567890123456789012345678901234
     .      /,'    1 --- read x|y, y|x, x,y (1,2,3)')
C +... TITRE

       number=  mx - imez - 1
       ncurv =  1
       forme ='(2f15.6)'
       formi ='(2i15  )'
       write(32,3202) number,ncurv,forme,formi
 3202  format(' ' ,i4,' --- number of data',
     .      /,'  ',i3,' --- number of curves',
     .      /,'  ',a8,
     .      /,'  ',a8)

C +--OUTPUT
C +  ~~~~~~
       do i=imez+1,mx1
       write(32,3203) distxx(i),1.d2*depths(i)
 3203  format(2f15.6)
       end do

C +--CLOSE OUTPUT FILE
C +  ~~~~~~~~~~~~~~~~~
       write(32,forme) facx,facy
       write(32,forme) cste,cste
       write(32,formi) isgn,isgn
       write(32,formi) mark,mark
                       month=0
       write(32, 3211) month,titxxx
 3211  format(
     .   '    0 --- axe x logarithmique :        (oui,non) = (1,0)',
     . /,'    0 --- axe y logarithmique :        (oui,non) = (1,0)',
     . /,  i5,' --- label axe inferieur : mois : (oui,non) = (1,0)',
     . /,a40)
       write(32, 3212)       tityyy,titbla,titbla
 3212  format(a40)

       filnex = 'n'
       write(32,3213)filnex,filn_h
 3213  format('    ',a1,'   0 --- Fichier Suivant / Type Graphe',
     .      /,'    ',a12)

       write(32,3214) jdarGE,labmGE(mmarGE),jhurGE,jdaMAR,jhaMAR,
     .           1.d2*depthv
 3214  format('    FILE NAME:           Sastrugi.out   ',
     .      /,'    Snow Erosion Dumont-d`Urville/DomeC ',
C +...             1234567890123456789012345678901234567890
     .      /,'    ',i2,'-',a3,2x,i2,'UT',9x,'After',i2,'d',i3,'h',3x,
     .      /,'    Mean Snow Thickness: ',f12.3,' cm    ')

       write(32,3215)
       write(32,3215)
 3215  format('                                            ',
     .      /,'                                            ')
       write(32,3216)
 3216  format('    GENERATION OUTPUT BY  KjTEST/SBCnew.2KJ')
       write(32,3217)
 3217  format(' -------------------------------------------')

       close(unit=32)

      end if

      if (iterun.eq.nterun-1) then
      close(unit=31)
      end if

      return
      end 


      subroutine ANI_nc

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                              3-10-2002  MAR |
C |   SubRoutine ANI_nc is used to write the main Model Variables          |
C |                                      on  a NetCDF file                 |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT: iprANI: Current time step    number                           |
C |   ^^^^^^         (starting from iprANI=1, which => new file creation)  |
C |          nprANI: Total  'time slices' number (max value of iprANI)     |
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
C +
      include 'MAR_DY.inc'
      include 'MAR_TE.inc'
      include 'MAR_TU.inc'
C +
      include 'MAR_HY.inc'
      include 'MAR_RA.inc'
c #TC include 'MAR_TC.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MARsSN.inc'

      include 'MAR_WK.inc'

      include 'MAR_IO.inc'

      integer           iprANI,nprANI,nt_ANI
      common/ANI_nc_arg/iprANI,nprANI,nt_ANI


C +--Local   Variables
C +  =================

      integer    Lfnam,     Ltit,     Luni,     Lnam,     Llnam
      PARAMETER (Lfnam= 40, Ltit= 90, Luni= 31, Lnam= 13, Llnam=50)
C +...Length of char strings 

      CHARACTER*(Lfnam)  fnamNC
      common/ANI_nc_loc/ fnamNC
C +...                   fnamNC: To retain file name.

      real               snow0(mx,my),bsnow0(mx,my),bsnow1(mx,my)
      common/ANI_rr_loc/ snow0       ,bsnow0       ,bsnow1
C +...                   snow0 : Snow Precipitated over Previous Time Interval
C +                     bsnow0 : Snow Budget       at   Previous Time Interval
C +                     bsnow1 : Snow Budget       at   Current  Time Interval

      integer    NdimNC
      PARAMETER (NdimNC = 5)
C +...Number of defined spatial dimensions (exact)

      integer    MXdim
      PARAMETER (MXdim  = 60000)
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
      parameter(mz_SBL=8)
      real      ua_SBL(mx,my,mz_SBL),va_SBL(mx,my,mz_SBL)
      real      Ta_SBL(mx,my,mz_SBL),zz_SBL(mx,my,mz_SBL)

      integer   n1000 ,n100a ,n100  ,n10_a ,n10   ,n1    
      integer   m10   ,       jd10  ,jd1
      integer   MMXstp,it    ,mois  ,mill  ,iu
      integer   itotNC,NtotNC,ID__nc
      real      starta(1)
      REAL      starti,DayLen,optwa ,optia ,rhodzk


C +--NetCDF File Initialization
C +  ==========================

      IF (iprANI.eq.1) THEN

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
c _UL   nDFdim(0)= nprANI
        nDFdim(0)= 0
        NAMdim(0)= 'time'
        UNIdim(0)= 'HOURS since 1981-01-15 00:00:00'

C +...  Check temporary arrays: large enough ?
        IF (nprANI.gt.MMXstp)
     &  STOP '*** ANI_nc - ERROR : MXdim to low ***'

              starti = jhurGE + minuGE/60.d0 + jsecGE/3600.d0
C +...        starti : Starting Time (= current time in the day)

              starta = (351+(iyrrGE  -1982) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1981) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
     .               + (minuGE *60 +jsecGE      )/3600.!

        DO it = 1,nprANI
              timeNC(it)   = starti    + (it-1) * nt_ANI  *dt / 3600.
C +...                                         nt_ANI: #iter between output
C +
              VALdim(it,0) = starta(1) + (it-1) * nt_ANI  *dt / 3600.
C +...        VALdim(  ,0) : values of the dimension # 0 (time) 

C +--Time Variable (date)
C +  ~~~~~~~~~~~~~~~~~~~~
              dateNC(it) =          timeNC(it)
              jourNC(it) = jdarGE + timeNC(it) / 24.d0
        END DO
                  mois       =  mmarGE
                  mill       =  iyrrGE
        DO it = 1,nprANI
          IF     (jourNC(it).gt.njmoGE(mois))                     THEN ! CTR
            DO iu=it,nprANI
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
            DO iu=it,nprANI
                  dateNC(iu) = mod(dateNC(iu),DayLen)
            END DO
          END IF                                                       ! CTR
        END DO

        DO it = 1,nprANI
              dateNC(it) =  dateNC(it)
     .             + 1.d+2 *jourNC(it)
     .             + 1.d+4 *moisNC(it)
        END DO


C +---Define horizontal spatial dimensions :    
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +...  Check temporary arrays: large enough ?
        IF (    mx .gt.MMXstp.or.my.gt.MMXstp
     &      .or.mzz.gt.MMXstp.or.mw.gt.MMXstp)
     &    STOP '*** ANI_nc - ERROR : MXdim to low ***'

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
C +
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
C +
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
              starta = (351+(iyrrGE  -1982) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1981) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
     .               + (minuGE *60 +jsecGE      )/3600.!

C +     ************
        CALL UNwrite (ID__nc, 'time   ',iprANI, 1, 1, 1, starta)
C +     ************

      END IF


C +     ************
        CALL UNwrite (ID__nc, 'date   ',iprANI, 1, 1, 1, dateNC(iprANI))
        CALL UNwrite (ID__nc, 'year   ',iprANI, 1, 1, 1, yearNC(iprANI))
C +     ************


C +--Dynamics, Precipitation
C +  -----------------------

      DO j=1,my
      DO i=1,mx
       DO k=1,mz_SBL
        ua_SBL(i,j,k) = uairDY(i,j,mzz-k)
        va_SBL(i,j,k) = vairDY(i,j,mzz-k)
        Ta_SBL(i,j,k) = tairDY(i,j,mzz-k)
        zz_SBL(i,j,k) = gplvDY(i,j,mzz-k)*grvinv
       END DO
         WKxy4(i,j)   = snowHY(i,j) - snow0(i,j)
         snow0(i,j)   = snowHY(i,j)


C +--Total   Snow Budget
C +  -------------------

        bsnow1(i,j)   =              snohSN(i,j,1)
     .                + (snowHY(i,j)-sno0HY(i,j  ))  *1.0e+3
       IF        (nssSNo(i,j,1).gt.   0      )                      THEN
       DO k=max(0,nssSNo(i,j,1)),nssSNo(i,j,1)
        bsnow1(i,j)   = bsnow1(i,j)+ rosSNo(i,j,1,k) *dzsSNo(i,j,1,k)
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
      CALL UNwrite (ID__nc, 'pstar  ', iprANI, mx, my, 1      , pstDY )
      CALL UNwrite (ID__nc, 'ua_SBL ', iprANI, mx, my, mz_SBL , ua_SBL)
      CALL UNwrite (ID__nc, 'va_SBL ', iprANI, mx, my, mz_SBL , va_SBL)
      CALL UNwrite (ID__nc, 'Ta_SBL ', iprANI, mx, my, mz_SBL , Ta_SBL)
      CALL UNwrite (ID__nc, 'zz_SBL ', iprANI, mx, my, mz_SBL , zz_SBL)
      CALL UNwrite (ID__nc, 'SnowTB ', iprANI, mx, my, 1      , WKxy1 )
      CALL UNwrite (ID__nc, 'SnowTR ', iprANI, mx, my, 1      , WKxy2 )
      CALL UNwrite (ID__nc, 'RadOLR ', iprANI, mx, my, 1      , RAdOLR)
      CALL UNwrite (ID__nc, 'OptDep ', iprANI, mx, my, 1      , WKxy3 )
      CALL UNwrite (ID__nc, 'Snow   ', iprANI, mx, my, 1      , WKxy4 )
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
      enddo
      enddo

      return
      end


      subroutine AWS_nc

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                             23-05-2004  MAR |
C |   SubRoutine AWS_nc is used to write an AWS Model  Variables           |
C |                                      on a   NetCDF File                |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   INPUT: iprAWS: Current time step    number                           |
C |   ^^^^^^         (starting from iprAWS=1, which => new file creation)  |
C |          nprAWS: Total  'time slices' number (max value of iprAWS)     |
C |                                                                        |
C |   OUTPUT: NetCDF File adapted to IDL Graphic Software                  |
C |   ^^^^^^                                                               |
C |                                                                        |
C |   CAUTION: 1) This Routine requires the usual NetCDF library,          |
C |   ^^^^^^^^    and the complementary access library  'libUN.a'          |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

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
      include 'MARdSV.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MAR_BS.inc'

      include 'MAR_WK.inc'

      include 'MAR_IO.inc'

      integer           iprAWS,nprAWS,nt_AWS
      common/AWS_nc_arg/iprAWS,nprAWS,nt_AWS


C +--AWS
C +  ---

      integer                  n_AWS,  n
      parameter               (n_AWS=  6)
      integer            AWSio(n_AWS),AWS_i(n_AWS),AWS_j(n_AWS)
      integer            nnAWS
      REAL               AWSla(n_AWS),AWSlo(n_AWS),AWS_z(n_AWS)
      character*3        AWS_0(n_AWS)
      character*8                     AWS_1(n_AWS),AWS_2(n_AWS)

      common /AWS_nc_INT/AWSio       ,AWS_i       ,AWS_j
     .                  ,nnAWS
      common /AWS_nc_REA/AWSla       ,AWSlo       ,AWS_z
      common /AWS_nc_CH3/AWS_0
      common /AWS_nc_CH8/             AWS_1       ,AWS_2


C +--Local   Variables
C +  =================

      integer            AWSij(iptx) ,ij

      REAL               distmn,Pt_la ,Pt_lo ,dd_la ,dd_lo ,Dista 
      REAL               dshmin,dsh
      integer            i__min,j__min,in_min,jn_min,l
      logical            SamAlt

      integer            ii           ,jj
      REAL               uu           ,vv
      REAL               FF(1)        ,DD(1)
      integer            nsnsno, ns
      parameter         (nsnsno= nsno+nsol+1)

      REAL               WK__1(1)     ,WK__2(1)     ,WK__3(1)
      REAL               WK__4(1)     ,WK__5(1)     ,WK__6(1)
      REAL               WK_z0(mz)    ,WK_z1(mz)    ,WK_z2(mz)  
      REAL               WK_z3(mz)    ,WK_z4(mz)    ,WK_z5(mz)
      REAL               WK_z6(mz)    ,WK_z7(mz)    ,WK_z8(mz)
      REAL               WK_z9(mz)
      REAL               WK_n1(nsnsno),WK_n2(nsnsno),WK_n3(nsnsno)
      REAL               WK_n4(nsnsno),WK_n5(nsnsno),WK_n6(nsnsno)
      REAL               WK_n7(nsnsno),WK_n8(nsnsno),WK_n9(nsnsno)

      REAL               VV_SBL,zvvSBL,TT_SBL,zttSBL,Kz_dz
      REAL               LMO   ,zetv_i,zetT_i,Psi__i,Psih_i
      REAL               xx2Psi,xx1Psi,yy2Psi,yy1Psi
      REAL               FF_AWS(1)
      REAL               TT_AWS(1)
      REAL               zvvAWS,zttAWS
      data               zvvAWS,zttAWS/3.0,3.0/


C +--NetCDF  Control Variables
C +  -------------------------

      integer    Lfnam,     Ltit,     Luni,     Lnam,     Llnam
      PARAMETER (Lfnam= 40, Ltit= 90, Luni= 31, Lnam= 13, Llnam=50)
C +...Length of char strings 

      CHARACTER*(Lfnam)  fnamNC
      common/AWS_ncl_ii/ fnamNC
C +...                   fnamNC: To retain file name.

      real               snowb(mx,my),rainb(mx,my)
      common/AWS_ncl_rr/ snowb       ,rainb
C +...                   snowb : Integrated Snow over Previous Time Interval
C +...                   rainb : Integrated Rain over Previous Time Interval

      integer    NdimNC
      PARAMETER (NdimNC = 5)
C +...Number of defined spatial dimensions (exact)

      integer    MXdim
      PARAMETER (MXdim = 45000)
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
      common/AWS_nc_r/  yearNC,dateNC
      real              VALdim(MXdim,0:NdimNC)
      integer           nDFdim(      0:NdimNC)
      common/AWS_nc_d/  nDFdim
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

      integer   n1000 ,n100a ,n100  ,n10_a ,n10   ,n1    
      integer   m10   ,       jd10  ,jd1
      integer   MMXstp,it    ,mois  ,mill  ,iu
      integer   itotNC,NtotNC,ID__nc
      real      starta(1)
      REAL      starti,DayLen,optwa ,optia ,rhodzk


C +--NetCDF File Initialization
C +  ==========================

C +--Common      Initialization
C +  --------------------------

      IF (iprAWS.eq.1)                                              THEN

C +--Output File Label
C +  ~~~~~~~~~~~~~~~~~
        n1000  = 1 +     iyrrGE/1000
        n100a  =     mod(iyrrGE,1000)
        n100   = 1 +     n100a /100
        n10_a  =     mod(n100a ,100)
        n10    = 1 +     n10_a /10
        n1     = 1 + mod(n10_a ,10)
        m10    = 1 +     mmarGE/10
        m1     = 1 + mod(mmarGE,10)
        jd10   = 1 +     jdarGE/10
        jd1    = 1 + mod(jdarGE,10)

        fnamNC = 'AWS.'
     .         // labnum(n1000) // labnum(n100)
     .         // labnum(  n10) // labnum(  n1)
     .         // labnum(  m10) // labnum(  m1)
     .         // labnum( jd10) // labnum( jd1)
     .         // '.' // explIO
     .         // '.nc    '

C +--Variable's Choice (Table AWSvar.dat)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OPEN(unit=10,status='unknown',file='AWSvou.dat')

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

          ENDIF
        GOTO 980
 990    CONTINUE

        CLOSE(unit=10)

        NtotNC = itotNC 
C +...  NtotNC : Total number of variables writen in NetCDF file.

C +--List of NetCDF Attributes given to all Variables
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +...  The "actual_range" is the (min,max) of all data for each variable:
        NAMrat(1) = 'actual_range'
        NvatNC(1) = 2

C +...  The "[var]_range"  is NOT of attribute type,
C +     it is a true variable containing the (min,max)
C +     for each level, for 4D (space+time) variables only
C +     (automatic handling by UN library; must be the LAST attribute)
        NAMrat(NattNC) = '[var]_range'
        NvatNC(NattNC) = 2

C +--Array Bounds Check Variable
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MMXstp = MXdim

C +--Time Variable (Hour)
C +  ~~~~~~~~~~~~~~~~~~~~
C +...Temporary Arrays Check: large enough?
        IF (nprAWS.gt.MMXstp)
     &  STOP '*** AWS_nc - ERROR : MXdim to low ***'

C +...NetCDF Dimensions (Size, Name, Unit):
c _UL   nDFdim(0)= nprAWS
        nDFdim(0)= 0
        NAMdim(0)= 'time'
        UNIdim(0)= 'HOURS since 1985-01-15 00:00:00'

              starti = jhurGE + minuGE/60.d0 + jsecGE/3600.d0
C +...        starti : Starting Time (= current time in the day)

              starta = (351+(iyrrGE  -1986) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1985) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
     .               + (minuGE *60 +jsecGE      )/3600.!

        DO it = 1,nprAWS
              timeNC(it)   = starti    + (it-1) * nt_AWS  *dt / 3600.
C +...                                         nt_AWS: #iter between output

              VALdim(it,0) = starta(1) + (it-1) * nt_AWS  *dt / 3600.
C +...        VALdim(  ,0) : values of the dimension # 0 (time) 

C +--Time Variable (Date)
C +  ~~~~~~~~~~~~~~~~~~~~
              dateNC(it) =          timeNC(it)
              jourNC(it) = jdarGE + timeNC(it) / 24.d0
        END DO
                  mois       =  mmarGE
                  mill       =  iyrrGE
        DO it = 1,nprAWS
          IF     (jourNC(it).gt.njmoGE(mois))                     THEN
            DO iu=it,nprAWS
                  jourNC(iu) =  jourNC(iu) - njmoGE(mois)
            END DO
                  mois       =  mois + 1
              IF (mois.gt.12)                                     THEN
                  mois       =         1
                  mill       =  mill + 1
              END IF
          END IF
                  moisNC(it) =  mois
                  yearNC(it) =  mill

          IF     (dateNC(it).gt.24.d0-epsi)                       THEN
                  DayLen     =  24.d0
            DO iu=it,nprAWS
                  dateNC(iu) = mod(dateNC(iu),DayLen)
            END DO
          END IF
        END DO

        DO it = 1,nprAWS
              dateNC(it) =  dateNC(it)
     .             + 1.d+2 *jourNC(it)
     .             + 1.d+4 *moisNC(it)
        END DO


C +--Vertical Variables
C +  ~~~~~~~~~~~~~~~~~~
C +...Temporary Arrays Check: large enough?
        IF (nsnsno.gt.MMXstp.OR.mz.gt.MMXstp
     .                      .OR.mw.gt.MMXstp)   
     .  STOP '*** AWS_nc - ERROR : MXdim to low ***'

C +...NetCDF Dimensions (Size, Name, Unit):
          nDFdim(1)   =  1
          NAMdim(1)   = 'x'
          UNIdim(1)   = 'deg'

          nDFdim(2)   =  1   
          NAMdim(2)   = 'y'
          UNIdim(2)   = 'deg'

        do k = 1, mz
          VALdim(k,3) =  sigma(mzz-k)
        enddo
          nDFdim(3)   =  mz
          NAMdim(3)   = 'level'
          UNIdim(3)   = '[index]'
C +...    For atmospheric layers k

        do k = 1, nsnsno
          VALdim(k,4) =  k
        enddo
          nDFdim(4)   =  nsnsno
          NAMdim(4)   = 'level2'
          UNIdim(4)   = '[index]'
C +...    For soil        layers k

        do k = 1, mw
          VALdim(k,5) = k 
        enddo
          nDFdim(5)   =  mw
          NAMdim(5)   = 'sector'
          UNIdim(5)   = '[index]'
C +...    For Surface Sectors 

      END IF


C +--Specific    Initialization
C +  --------------------------

      DO n = 1,nnAWS

          ii          = AWS_i(n)
          jj          = AWS_j(n)

          fnamNC(1:3) = AWS_0(n)

        IF (iprAWS.eq.1)                                            THEN

          tit_NC      = '(A)WS ' // AWS_1(n) // AWS_2(n) // ' / '
     .               // ' OUTPUT of MAR (Modele Atmospherique Regional)'
     .               // ' / Experiment ' // explIO //  '  '
C +...                   1234567890123456789012345678901234567890123456
c #VER    write(6,600) tit_NC
  600     format(/,a90,/,9('1234567890'))

          k=1
          l=  Ltit
 601      CONTINUE
              IF      (tit_NC(k  :k+1).EQ.'  ')                     THEN
                       tit_NC(k  :l-1) =  tit_NC(k+1:l)
                       tit_NC(l  :l  ) =  ' '
                              k  =k-1
                              l  =l-1
              END IF
                              k  =k+1
          IF (k.LT.l)                                          GO TO 601

c #VER    write(6,600) tit_NC

          VALdim(1,1) = AWSlo(n)/degrad
          VALdim(1,2) = AWSla(n)/degrad

C +--Automatic Generation of the NetCDF File Structure
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +       **************
          CALL UNscreate(fnamNC, tit_NC,
     &                   NdimNC, nDFdim, MXdim , NAMdim, UNIdim, VALdim,
     &                   MX_var, NtotNC, nameNC, SdimNC, unitNC, lnamNC,
     &                   NattNC, NAMrat, NvatNC,
     &                   ID__nc) 
C +       **************


C +--Write Time and Geographic Coordinates
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          WK__1 =  GElonh(ii,jj) * 15.000
C +...    Conversion: Hour->degrees

          WK__2 =  GElatr(ii,jj) / degrad
C +...    Conversion: rad ->degrees

          WK__3 =      sh(ii,jj)

          WK__4 =  isolSL(ii,jj)
C +...    Conversion to REAL type (integer not allowed)

C +       ************
          CALL UNwrite (ID__nc, 'sh_AWS', 1  , 1  , 1 , 1 , AWS_z(n))
          CALL UNwrite (ID__nc, 'lonMAR', 1  , 1  , 1 , 1 , WK__1)
          CALL UNwrite (ID__nc, 'latMAR', 1  , 1  , 1 , 1 , WK__2)
          CALL UNwrite (ID__nc, 'sh_MAR', 1  , 1  , 1 , 1 , WK__3)
          CALL UNwrite (ID__nc, 'solTyp', 1  , 1  , 1 , 1 , WK__4)
C +       ************


C +--Re-Open file if already created.
C +  ================================


        ELSE

C +       ************
          CALL UNwopen (fnamNC,ID__nc)
C +       ************

        END IF


C +--Write Time-dependent variables:
C +  ===============================

C +--UNLIMITED Time Dimension
C +  ------------------------

        IF (nDFdim(0).eq.0)                         THEN !
                starta = (351+(iyrrGE  -1986) *365       ! Nb Days before iyrrGE
     .                       +(iyrrGE  -1985) /  4       ! Nb Leap Years
     .                       + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                       + njybGE(mmarGE)            ! (including Leap Day)
     .                   *max(0,1-mod(iyrrGE,4))         !
     .                       + jdarGE     -1      )*  24 !
     .                   +jhurGE                         !
     .                 + (minuGE *60 +jsecGE      )/3600.!

C +       ************
          CALL UNwrite (ID__nc, 'time   ',iprAWS, 1, 1, 1, starta)
C +       ************

        END IF


C +     ************
        CALL UNwrite (ID__nc, 'date   ',iprAWS, 1, 1, 1, dateNC(iprAWS))
        CALL UNwrite (ID__nc, 'year   ',iprAWS, 1, 1, 1, yearNC(iprAWS))
C +     ************


C +--Radiative Transfert
C +  -------------------

          WK__1    = RAdOLR(ii,jj)
          WK__2    = RAd_ir(ii,jj)
          WK__3    = RAdsol(ii,jj)
          WK__4    = RAcdtO(ii,jj)
          WK__5    = cld_SL(ii,jj)
          WK__6    = RAertO(ii,jj)

C +     ************
        CALL UNwrite (ID__nc, 'OLR    ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'DownLW ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'DownSW ', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'CloudOD', iprAWS, 1, 1, 1      , WK__4)
        CALL UNwrite (ID__nc, 'CloudFR', iprAWS, 1, 1, 1      , WK__5)
        CALL UNwrite (ID__nc, 'AerosOD', iprAWS, 1, 1, 1      , WK__6)
C +     ************


C +--Dynamics, Microphysics, Turbulence
C +  ----------------------------------

        DO k=1,mz
          WK_z0(k) = gplvDY(ii,jj,mzz-k)*grvinv
          WK_z1(k) = uairDY(ii,jj,mzz-k)
          WK_z2(k) = vairDY(ii,jj,mzz-k)
          WK_z9(k) = ect_TE(ii,jj,mzz-k)
          WK_z3(k) = tairDY(ii,jj,mzz-k)
          WK_z4(k) =   qvDY(ii,jj,mzz-k)         *1000.
          WK_z5(k) =   qiHY(ii,jj,mzz-k)         *1000.
          WK_z6(k) =   qwHY(ii,jj,mzz-k)         *1000.
          WK_z7(k) =   qsHY(ii,jj,mzz-k)         *1000.
          WK_z8(k) =   qrHY(ii,jj,mzz-k)         *1000.
        END DO
          uu       =  WK_z1(1)
          vv       =  WK_z2(1)


C +     ************
        CALL UNwrite (ID__nc, 'ZZ     ', iprAWS, 1, 1, mz     , WK_z0)
        CALL UNwrite (ID__nc, 'UU     ', iprAWS, 1, 1, mz     , WK_z1)
        CALL UNwrite (ID__nc, 'VV     ', iprAWS, 1, 1, mz     , WK_z2)
        CALL UNwrite (ID__nc, 'TKE    ', iprAWS, 1, 1, mz     , WK_z9)
        CALL UNwrite (ID__nc, 'TT     ', iprAWS, 1, 1, mz     , WK_z3)
        CALL UNwrite (ID__nc, 'QQ     ', iprAWS, 1, 1, mz     , WK_z4)
        CALL UNwrite (ID__nc, 'QI     ', iprAWS, 1, 1, mz     , WK_z5)
        CALL UNwrite (ID__nc, 'QW     ', iprAWS, 1, 1, mz     , WK_z6)
        CALL UNwrite (ID__nc, 'QS     ', iprAWS, 1, 1, mz     , WK_z7)
        CALL UNwrite (ID__nc, 'QR     ', iprAWS, 1, 1, mz     , WK_z8)
C +     ************


C +--Precipitation
C +  -------------

          WK__1    =(snowHY(ii,jj) - snowb (ii,jj)) *1000.
          WK__2    =(rainHY(ii,jj) - rainb (ii,jj)) *1000.

C +     ************
        CALL UNwrite (ID__nc, 'Snow   ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'Rain   ', iprAWS, 1, 1, 1      , WK__2)
C +     ************


C +--SBL
C +  ---

C +--Momentum Turbulent Fluxes
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~
           k=1
          WK_z4(1) =  TUkvm(ii,jj,mzz-k)
          Kz_dz    = -TUkvm(ii,jj,mzz-k) / (WK_z0(1)-   sh(ii,jj))
          WK_z5(k) =  Kz_dz              *  WK_z1(1)
          WK_z6(k) =  Kz_dz              *  WK_z2(1)
        DO k=2,mz
          WK_z4(k) =  TUkvm(ii,jj,mzz-k)
          Kz_dz    = -TUkvm(ii,jj,mzz-k) / (WK_z0(k)-WK_z0(k-1))
          WK_z5(k) =  Kz_dz              * (WK_z1(k)-WK_z1(k-1))
          WK_z6(k) =  Kz_dz              * (WK_z2(k)-WK_z2(k-1))
        ENDDO

C +     ************
        CALL UNwrite (ID__nc, 'Kz     ', iprAWS, 1, 1, mz     , WK_z4)
        CALL UNwrite (ID__nc, 'UpWp   ', iprAWS, 1, 1, mz     , WK_z5)
        CALL UNwrite (ID__nc, 'VpWp   ', iprAWS, 1, 1, mz     , WK_z6)
C +     ************

C +--Wind Gusts
C +  ~~~~~~~~~~
          WK__1    = SL_wge(ii,jj) 
          WK__2    = SLlwge(ii,jj) 
          WK__3    = SLuwge(ii,jj)

C +     ************
        CALL UNwrite (ID__nc, 'WGE    ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'WGE_LB ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'WGE_UB ', iprAWS, 1, 1, 1      , WK__3)
C +     ************

C +--3m Meteorological Variables
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
          VV_SBL=              ssvSL(ii,jj,mz)                     ! Refer. Level Wind
          zvvSBL=    grvinv * gplvDY(ii,jj,mz)-    sh(ii,jj)       ! Refer. Level Height
          TT_SBL=             tairDY(ii,jj,mz)                     ! Refer. Level Temper.
          zttSBL=             gplvDY(ii,jj,mz)-    sh(ii,jj)       ! Refer. Level Height
          LMO   =    grvinv * SLlmol(ii,jj,1)                      ! Monin-Obukhov Lenght
                                                                   !
        IF (LMO.GT.0)                                         THEN ! STABLE   Situation
          zetv_i=                         zvvAWS/max(epsi,LMO)     ! Stability Limit
          zetv_i=               zetv_i-   zvvSBL/max(epsi,LMO)     !
          zetT_i=                         zttAWS/max(epsi,LMO)     ! zeta MAX = 4.28
          zetT_i=               zetT_i-   zttSBL/max(epsi,LMO)     !
          Psi__i=  -6.0*sign(1.,zetv_i)*min(4.28,abs(zetv_i))      !
          Psih_i=  -6.0*sign(1.,zetT_i)*min(4.28,abs(zetT_i))      !
        ELSE                                                       ! UNSTABLE Situation
          xx2Psi=sqrt(sqrt(1.-15.*zvvAWS       /min(-epsi,LMO)))   !
          xx1Psi=sqrt(sqrt(1.-15.*      zvvSBL /min(-epsi,LMO)))   !
          Psi__i=2. * log(0.5*(1.+xx2Psi))                         ! 
     .          +     log(0.5*(1.+xx2Psi*xx2Psi))                  ! 
     .          -2. *atan(        xx2Psi        )                  ! 
     .          -2. * log(0.5*(1.+xx1Psi))                         ! 
     .          -     log(0.5*(1.+xx1Psi*xx1Psi))                  ! 
     .          +2. *atan(        xx1Psi        )                  ! 
          yy2Psi=     sqrt(1.- 9.* zttAWS      /min(-epsi,LMO))    !
          yy1Psi=     sqrt(1.- 9.*      zttSBL /min(-epsi,LMO))    !
          Psih_i=     log(0.5*(1.+yy2Psi))                         ! 
     .          -     log(0.5*(1.+yy1Psi))                         !
        END IF                                                     !
          FF_AWS     = VV_SBL                                      ! Businger, 1973, 
     .               + SLuusl(ii,jj,1)*(log(zvvAWS/zvvSBL)-Psi__i) ! WkShop on mi-mto,
     .                                / 0.4                        !       (3.7) p.77
          FF_AWS     =              max(0.0,FF_AWS)                !
          TT_AWS     = TT_SBL                                      ! Businger, 1973, 
     .          + 0.74*SLutsl(ii,jj,1)*(log(zttAWS/zttSBL)-Psih_i) ! WkShop on mi-mto,
     .       /max(epsi,SLuusl(ii,jj,1))                            !       (3.8) p.77

C +     ************
        CALL UNwrite (ID__nc, 'FF_AWS ', iprAWS, 1, 1, 1      ,FF_AWS)
        CALL UNwrite (ID__nc, 'TT_AWS ', iprAWS, 1, 1, 1      ,TT_AWS)
C +     ************

C +   Wind Direction and Wind Speed 
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +          *****
        call uv2fd(uu,vv,FF(1),DD(1),ii,jj)
C +          *****

C +--Surface Thermodynamics
C +  ~~~~~~~~~~~~~~~~~~~~~~
          WK__3    =(pstDy (ii,jj)     + ptopDY)  * 10.
          WK__4    = TairSL(ii,jj)
          WK__5    = TsrfSL(ii,jj,1)

C +     ************
        CALL UNwrite (ID__nc, 'FF     ', iprAWS, 1, 1, 1      , FF   )
        CALL UNwrite (ID__nc, 'DD     ', iprAWS, 1, 1, 1      , DD   )
        CALL UNwrite (ID__nc, 'ppSurf ', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'TaSurf ', iprAWS, 1, 1, 1      , WK__4)
        CALL UNwrite (ID__nc, 'TTSurf ', iprAWS, 1, 1, 1      , WK__5)
C +     ************

C +--SBL Turbulence
C +  ~~~~~~~~~~~~~~
          WK__1    = SLuuSL(ii,jj,1)
          WK__2    = SLutSL(ii,jj,1)
          WK__3    = SLuqSL(ii,jj,1)
          WK__4    = SLusSL(ii,jj,1)
          WK__5    = SaltSN(ii,jj,1)

C +     ************
        CALL UNwrite (ID__nc, 'u_star ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'uT_star', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'uq_star', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'us_star', iprAWS, 1, 1, 1      , WK__4)
        CALL UNwrite (ID__nc, 'us_BSth', iprAWS, 1, 1, 1      , WK__5)
C +     ************

          WK__1    = SL_z0( ii,jj,1) * 1.e3
          WK__2    = SL_r0( ii,jj,1) * 1.e3
          WK__3    = Z0emBS(ii,jj,1) * 1.e3

C +     ************
        CALL UNwrite (ID__nc, 'z0_m   ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'z0_h   ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'z0_eff ', iprAWS, 1, 1, 1      , WK__3)
C +     ************


C +--Soil, Ice and Snow
C +  ------------------

          WK__1    = albeSL(ii,jj)
          WK__2    = nssSNo(ii,jj,1)
          WK__3    = nisSNo(ii,jj,1)
          WK__4    = issSNo(ii,jj,1)

C +     ************
        CALL UNwrite (ID__nc, 'Albedo ', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'n_SNOW ', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'n__ICE ', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'nSUPER ', iprAWS, 1, 1, 1      , WK__4)
C +     ************

          DO k=1,nsnsno
          WK_n1(k) = 0.
          WK_n2(k) = 0.
          WK_n3(k) = 0.
          WK_n4(k) = 0.
          WK_n5(k) = 0.
          WK_n6(k) = 0.
          WK_n7(k) = 0.
          WK_n8(k) = 0.
          END DO

          WK__4    = 0.
        IF          (nssSNo(ii,jj,1).GT.0)              THEN
          ns       = nssSNo(ii,jj,1)
          DO k=1,ns
          WK_n1(k) = tisSNo(ii,jj,1,    ns+1-k)
          WK_n2(k) = dzsSNo(ii,jj,1,    ns+1-k)
          WK_n3(k) = rosSNo(ii,jj,1,    ns+1-k)
          WK_n4(k) = wasSNo(ii,jj,1,    ns+1-k)
          WK_n5(k) = g1sSNo(ii,jj,1,    ns+1-k)
          WK_n6(k) = g2sSNo(ii,jj,1,    ns+1-k)
          WK_n7(k) = agsSNo(ii,jj,1,    ns+1-k)
          WK_n8(k) = nhsSNo(ii,jj,1,    ns+1-k)
          WK__4    = WK__4 
     .             + rosSNo(ii,jj,1,         k)
     .              *wasSNo(ii,jj,1,         k)
     .             + dzsSNo(ii,jj,1,         k)
          END DO
        ELSE
          ns = 0
        END IF

          DO k=ns+1,ns+llx
          WK_n1(k) = TsolTV(ii,jj,1,llx+ns+1-k)
          WK_n2(k) = dz_dSV(            ns+1-k)
          WK_n3(k) = ro_Ice
          WK_n4(k) = eta_TV(ii,jj,1,llx+ns+1-k)
          END DO


C +     ************
        CALL UNwrite (ID__nc, 'T__SIS ', iprAWS, 1, 1, nsnsno , WK_n1)
        CALL UNwrite (ID__nc, 'dz_SIS ', iprAWS, 1, 1, nsnsno , WK_n2)
        CALL UNwrite (ID__nc, 'ro_SIS ', iprAWS, 1, 1, nsnsno , WK_n3)
        CALL UNwrite (ID__nc, 'wa_SIS ', iprAWS, 1, 1, nsnsno , WK_n4)
        CALL UNwrite (ID__nc, 'G1SNOW ', iprAWS, 1, 1, nsnsno , WK_n5)
        CALL UNwrite (ID__nc, 'G2SNOW ', iprAWS, 1, 1, nsnsno , WK_n6)
        CALL UNwrite (ID__nc, 'AgSNOW ', iprAWS, 1, 1, nsnsno , WK_n7)
        CALL UNwrite (ID__nc, 'HiSNOW ', iprAWS, 1, 1, nsnsno , WK_n8)
C +     ************


C +--Water Budget
C +  ------------

          WK__1    = runoTV(ii,jj)
          WK__2    = evapTV(ii,jj)
          WK__3    = SWaSNo(ii,jj,1)

C +     ************
        CALL UNwrite (ID__nc, 'IntROFF', iprAWS, 1, 1, 1      , WK__1)
        CALL UNwrite (ID__nc, 'IntEVAP', iprAWS, 1, 1, 1      , WK__2)
        CALL UNwrite (ID__nc, 'SurfWAT', iprAWS, 1, 1, 1      , WK__3)
        CALL UNwrite (ID__nc, 'VIntWAT', iprAWS, 1, 1, 1      , WK__4)
C +     ************


C +--That 's all, folks: NetCDF File Closure
C +  =======================================

C +     ***********
        CALL NCCLOS (ID__nc,RCODE)
C +     ***********

      ENDDO


C +--Work Arrays Reset
C +  =================

      do j=1,my
      do i=1,mx
        rainb(i,j) = rainHY(i,j)
        snowb(i,j) = snowHY(i,j)
      enddo
      enddo


      return
      end


      subroutine AWSloc

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                             20-04-2004  MAR |
C |   SubRoutine AWSloc searches AWS and Manned Stations on the Model Grid |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      include 'MAR_SL.inc'
      include 'MAR_SV.inc'
      include 'MAR_TV.inc'

      include 'MAR_IO.inc'


C +--AWS
C +  ---

      integer                  n_AWS,  n
      parameter               (n_AWS=  6)
      integer            AWSio(n_AWS),AWS_i(n_AWS),AWS_j(n_AWS)
      integer            nnAWS
      REAL               AWSla(n_AWS),AWSlo(n_AWS),AWS_z(n_AWS)
      REAL               AWS_x(n_AWS),AWS_y(n_AWS),AWSdd
      character*3        AWS_0(n_AWS)
      character*8                     AWS_1(n_AWS),AWS_2(n_AWS)

      common /AWS_nc_INT/AWSio       ,AWS_i       ,AWS_j
     .                  ,nnAWS
      common /AWS_nc_REA/AWSla       ,AWSlo       ,AWS_z
      common /AWS_nc_CH3/AWS_0
      common /AWS_nc_CH8/             AWS_1       ,AWS_2


      integer            AWSij(iptx)  ,ij

      REAL               distmn,Pt_la ,Pt_lo ,dd_la ,dd_lo ,Dista 
      REAL               dshmin,dsh
      integer            i__min,j__min,in_min,jn_min,l
      logical            SamAlt


C +--DATA
C +  ====

      data      SamAlt/.false./    ! AWS Model Grid Point Refinement SWITCH


C +--Geographic and Grid Point Coordinates for IO
C +  ============================================

      open( unit=30,status='unknown',file='SBCnew.AWS')
      rewind     30
      open( unit=31,status='unknown',file='SBCnew.AWS.JNL')
      rewind     31
           write( 6,6000) 
           write(30,6000) 
 6000      format('    AWS Station | Latit. | Longit.|'
     .     ,                       ' x [km] | y [km] | Altit. ||'
     .     ,             ' Grid pt.| x [km] | y [km] |'
     .     ,                       ' Latit. | Longit.| Altit. ||'
     .     ,            ' D(AWS,pt)|')
           write( 6,6002) 
           write(30,6002) 
 6002      format('----------------+--------+--------+'
     .     ,                       '--------+--------+--------++'
     .     ,             '---------+--------+--------+'
     .     ,                       '--------+--------+--------++'
     .     ,            '----------+')


C +--Find the closest Grid Point in the MAR Domain
C +  ---------------------------------------------

              dd_lo    = 0.
      DO n=1,n_AWS
              distmn   =         mx*mx +my*my
              distmn   = dx*sqrt(distmn)
              AWSla(n) = AWSla(n)    * degrad 
              AWSlo(n) = AWSlo(n)    * degrad 
        DO j=1,my
        DO i=1,mx
              Pt_la    = GElatr(i,j) 
c #3D         Pt_lo    = GElonh(i,j) * hourad
              dd_la    = earthr      *              (AWSla(n)-Pt_la)
c #3D         dd_lo    = earthr      * cos(Pt_la) * (AWSlo(n)-Pt_lo)
              Dista    = sqrt(dd_la*dd_la+dd_lo*dd_lo)
          IF (Dista.lt.distmn)                                      THEN
              distmn   = Dista
              i__min   = i
              j__min   = j
          END IF
c #WR         write(6,6) n,i,j,GElatr(i,j)/degrad,GElonh(i,j)*hourad
c #WR.                                                       /degrad
c #WR.                        , AWSla(n)  /degrad, AWSlo(n)  /degrad
c #WR.                  ,1.e-3*Dista 
 6            format(3i3,4f9.3,f12.3)
        END DO
        END DO


C +--A Grid Point is found in the MAR Domain
C +  ---------------------------------------

        IF (distmn.LT.2.*dx)                                        THEN

C +--(x,y) Coordinates of the AWS
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 AWS_x(n) =  0.
                 AWS_y(n) =  0.
                 AWSdd    =  0.
            DO j=max(1,j__min-1),min(j__min+1,my)
            DO i=max(1,i__min-1),min(i__min+1,mx)
                 Pt_la    = GElatr(i,j) 
c #3D            Pt_lo    = GElonh(i,j) * hourad
                 dd_la    = earthr      *              (AWSla(n)-Pt_la)
c #3D            dd_lo    = earthr      * cos(Pt_la) * (AWSlo(n)-Pt_lo)
                 Dista    =      max(epsi,sqrt(dd_la*dd_la+dd_lo*dd_lo))
                 AWS_x(n) =  AWS_x(n)   + dx*(i-imez) / (Dista * Dista)
                 AWS_y(n) =  AWS_y(n)   + dx*(j-jmez) / (Dista * Dista)
                 AWSdd    =  AWSdd      + 1.          / (Dista * Dista)
            ENDDO
            ENDDO
                 AWS_x(n) =  AWS_x(n)   * 1.0e-3      /  AWSdd
                 AWS_y(n) =  AWS_y(n)   * 1.0e-3      /  AWSdd

C +--CRITERION - CRITERION - CRITERION - CRITERION - CRITERION - CRITERION 
C +  Determine  a MAR Grid Point with the closest Altitude in the Vicinity
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF (SamAlt.AND.abs(sh(i__min,j__min)-AWS_z(n)).GT.100.)   THEN
                 dshmin   = 1.e6
            DO j=max(1,j__min-1),min(j__min+1,my)
            DO i=max(1,i__min-1),min(i__min+1,mx)
                 Pt_la    = GElatr(i,j) 
c #3D            Pt_lo    = GElonh(i,j) * hourad
                 dd_la    = earthr      *              (AWSla(n)-Pt_la)
c #3D            dd_lo    = earthr      * cos(Pt_la) * (AWSlo(n)-Pt_lo)
                 Dista    = sqrt(dd_la*dd_la+dd_lo*dd_lo)
                 dsh      =  abs(sh(i,j)    -AWS_z(n)   )
              IF(isolSL(i,j).GE.3                     .AND.
     .           dsh        .LT.dshmin                .AND.
     .           Dista      .LT.dx*1.5                     )        THEN
                 in_min   = i
                 jn_min   = j
                 dshmin   = dsh
                 distmn   = Dista
              END IF
            ENDDO
            ENDDO
                 i__min   = in_min
                 j__min   = jn_min
                 i        = in_min
                 j        = jn_min
          ENDIF
               AWS_i(n)   = i__min
               AWS_j(n)   = j__min

C +--Summary Table
C +  ~~~~~~~~~~~~~
           write( 6,6001) AWS_1(n)       ,AWS_2(n)
     .                   ,AWSla(n)/degrad,AWSlo(n)/degrad
     .                   ,AWS_x(n)       ,AWS_y(n)
     .                   ,AWS_z(n)
     .                   ,       i__min,j__min
     .                   ,   dx*(i__min-imez)     * 1.e-3
     .                   ,   dx*(j__min-jmez)     * 1.e-3
     .                   ,GElatr(i__min,j__min)   /degrad
     .                   ,GElonh(i__min,j__min)   *hourad/degrad
     .                   ,    sh(i__min,j__min),    1.e-3*distmn
           write(30,6001) AWS_1(n)       ,AWS_2(n)
     .                   ,AWSla(n)/degrad,AWSlo(n)/degrad
     .                   ,AWS_x(n)       ,AWS_y(n)
     .                   ,AWS_z(n)
     .                   ,       i__min,j__min
     .                   ,   dx*(i__min-imez)     * 1.e-3
     .                   ,   dx*(j__min-jmez)     * 1.e-3
     .                   ,GElatr(i__min,j__min)   /degrad
     .                   ,GElonh(i__min,j__min)   *hourad/degrad
     .                   ,    sh(i__min,j__min),    1.e-3*distmn
 6001      format(2a8, '|',2(f7.2,' |'),2(f7.1,' |'),f7.0,' ||',
     .            2i4,' |',2(f7.1,' |'),
     .                     2(f7.2,' |'),  f7.0,' ||',f9.1,' |')
           write(31,3100)    dx*(i__min-imez)     * 1.e-3
     .                   ,   dx*(j__min-jmez)     * 1.e-3
     .                   ,AWS_0(n)
     .                   ,   dx*(i__min-imez)     * 1.e-3
     .                   ,   dx*(j__min-jmez)     * 1.e-3
     .                   ,   AWS_x(n)
     .                   ,   AWS_y(n)
c #OU.                   ,AWS_1(n)
c #OU.                   ,   dx*(i__min-imez)     * 1.e-3
c #OU.                   ,   dx*(j__min-jmez)     * 1.e-3
 3100      format('LABEL ',2(f8.2,','),'0,0,.10 ',a3
     .         ,/,'LABEL ',2(f8.2,','),'0,0,.10 o'
     .         ,/,'LABEL ',2(f8.2,','),'0,0,.10 +'
c #OU.         ,/,'sp echo " "           >>dsSWEq.LIST'
c #OU.         ,/,'sp echo "',a8,  'SMB" >>dsSWEq.LIST'
c #OU.         ,/,'set region/x=',f8.2,'/y=',f8.2
c #OU.         ,/,'list/append/nohead/file=dsSWEq.LIST '
c #OU.           ,                        'd_BSWE,dnBSWE'
     .                                                   )


C +--Stations which are not in the Model Domain are discarded
C +  --------------------------------------------------------
        ELSE
               AWSio(n)  = 0
               AWS_i(n)  = 0
               AWS_j(n)  = 0
        END IF

      END DO
             nnAWS       = 0
      DO n=1,n_AWS
        IF    (AWSio(n).GT.0)                                   THEN
             nnAWS       = nnAWS + 1
             AWSio(nnAWS)= AWSio(n)
             AWS_i(nnAWS)= AWS_i(n)
             AWS_j(nnAWS)= AWS_j(n)
             AWS_x(nnAWS)= AWS_x(n)
             AWS_y(nnAWS)= AWS_y(n)
             AWSla(nnAWS)= AWSla(n)
             AWSlo(nnAWS)= AWSlo(n)
             AWS_z(nnAWS)= AWS_z(n)
             AWS_0(nnAWS)= AWS_0(n)
             AWS_1(nnAWS)= AWS_1(n)
             AWS_2(nnAWS)= AWS_2(n)
        ENDIF
      ENDDO


C +--Modification of the Indices for SISVAT ASCII OUTPUT
C +  ---------------------------------------------------

                ij =       0
      DO n=1,nnAWS
          IF (AWSio (n).EQ.2)                                   THEN
                ij = ij  + 1
            IF (ij.LE.iptx)                                     THEN
              AWSij (ij) =       n
              IOi_TV(ij) = AWS_i(n)
              IOj_TV(ij) = AWS_j(n)
            END IF
          END IF
      ENDDO
      IF (ij.LT.iptx)                                           THEN
        DO n=ij+1,iptx
              AWSij (n)  = AWSij(ij)
              IOi_TV(n)  = AWS_i(ij)
              IOj_TV(n)  = AWS_j(ij)
        ENDDO
      ENDIF


           write( 6,6002) 
           write(30,6002) 

           write( 6,*)  ' '
           write(30,*)  ' '

        DO ij=1,iptx
           write( 6,6003) AWS_1(AWSij(ij)), AWS_2(AWSij(ij))
     .                  ,IOi_TV(      ij ),IOj_TV(      ij )
 6003      format('SISVAT OUTPUT for AWS ',2a8,
     .            ' at Pt (',2i4,')')
        ENDDO
                  
           write( 6,*)  ' '
           write(30,*)  ' '

      close(unit=30)
      close(unit=31)

      return
      end


      subroutine uv2fd(ruu,rvv,rFF,rDD,iu,ju)

C +------------------------------------------------------------------------+
C | MAR OUTPUT                                             04-06-2004  MAR |
C |   SubRoutine uv2fd  transforms Wind(u,v) into Wind(ff,dd)              |
C |                                                                        |
C |   CAUTION: must be modified to take into account the 2-D Domain        |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARphy.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'


C +--Local   Variables
C +  =================

      integer            iu   ,ju
      REAL               uugeo,ruu,rFF
      REAL               vvgeo,rvv,rDD

      REAL               conv
      parameter         (conv  = 15.0*3.141592/180.0)          ! Conversion
                                                               ! Hour -> Radian

C +--From (u,v) to (ff,dd)
C +  =====================

c #     IF   (ruu   .NE.0.0 .OR.  rvv   .NE.0.0)                    THEN

c #           uugeo    = (GElonh(iu+1,ju) - GElonh(iu,ju))*conv/dx
c #  .                   *earthr       *cos(GElatr(iu,ju))*    ruu
c #  .                 + (GElonh(iu,ju+1) - GElonh(iu,ju))*conv/dx
c #  .                   *earthr       *cos(GElatr(iu,ju))*    rvv
     
c #           vvgeo    = (GElatr(iu+1,ju) - GElatr(iu,ju))     /dx
c #  .                   *earthr                          *    ruu
c #  .                 + (GElatr(iu,ju+1) - GElatr(iu,ju))     /dx
c #  .                   *earthr                          *    rvv

              rDD       =  0.0

c #       IF (uugeo    .GT. 0.0 .and. vvgeo    .GE.0.0)  
c #  .        rDD       =    3.0*pi/2.0 - atan(vvgeo    /uugeo    )
c #       IF (uugeo    .LT. 0.0 .and. vvgeo    .GE.0.0)  
c #  .        rDD       =    5.0*Pi/2.0  -atan(vvgeo    /uugeo    )
c #       IF (uugeo    .LT. 0.0 .and. vvgeo    .LE.0.0)  
c #  .        rDD       =    5.0*Pi/2.0  -atan(vvgeo    /uugeo    )
c #       IF (uugeo    .GT. 0.0 .and. vvgeo    .LE.0.0)  
c #  .        rDD       =    3.0*Pi/2.0  -atan(vvgeo    /uugeo    )
c #       IF (uugeo    .EQ. 0.0 .and. vvgeo    .GE.0.0) 
c #  .        rDD       =        Pi
c #       IF (uugeo    .EQ. 0.0 .and. vvgeo    .LE.0.0) 
c #  .        rDD       =        0.0
c #       IF (rDD       .GT. 2.0*Pi)
c #  .        rDD       = rDD-2.0*Pi

c #           rDD       = rDD        / degrad

c #       if( rDD       .LT. 0.0 )       
c #  .        rDD       = rDD        + 360.0
     
c #           rDD       = max(0.0,min(360.0,rDD))
     
              rFF      = sqrt(ruu*ruu + rvv*rvv)
       
c #     ELSE
c #           rFF     = 0.0 
c #           rDD     = 0.0
c #     END IF 

      return
      end


      block data AWS_nc_DATA

C +----------------------------------------------------------------------------+
C |                                                                            |
C | MAR OUTPUT   Generic Routine                               20-05-2004  MAR |
C |   Manned and Automatic Weather Stations (AWS) Geographic Coordinates       |
C |                                                                            |
C +----------------------------------------------------------------------------+


C +--General Variables
C +  =================

      integer                  n_AWS,  n
      parameter               (n_AWS=  6)
      integer            AWSio(n_AWS),AWS_i(n_AWS),AWS_j(n_AWS)
      integer            nnAWS
      REAL               AWSla(n_AWS),AWSlo(n_AWS),AWS_z(n_AWS)
      character*3        AWS_0(n_AWS)
      character*8                     AWS_1(n_AWS),AWS_2(n_AWS)

      common /AWS_nc_INT/AWSio       ,AWS_i       ,AWS_j
     .                  ,nnAWS
      common /AWS_nc_REA/AWSla       ,AWSlo       ,AWS_z
      common /AWS_nc_CH3/AWS_0
      common /AWS_nc_CH8/             AWS_1       ,AWS_2


C +--DATA
C +  ====

C +--ANT
C +  ---

      data (AWS_0(n),AWS_1(n),AWS_2(n)
     .     ,AWSla(n),AWSlo(n),AWS_z(n),AWSio(n),n=001,n_AWS)
C +...      LABel      AWS LABELS            Latit. Longit. Alti. PR
C +                                                                0 => No IO
C +                                                                1 => OUTone
C +                                                                2 => OUTone
C +                                                                     ASCII
     .  /  'DDU'    ,'DDU     ','        ', -66.67, 140.02,   42., 2,  !  01
     .     'D10'    ,'D-10    ','        ', -66.71, 139.83,  243., 2,  !  02
     .     'D47'    ,'D-47    ','        ', -67.40, 138.73, 1560., 2,  !  03
     .     'D57'    ,'D-57    ','        ', -68.30, 137.87, 2105., 2,  !  04
     .     'D80'    ,'D80     ','        ', -70.02, 134.72, 2500., 2,  !  05
     .     'DCA'    ,'Dome  C ','AMRC    ', -74.50, 123.00, 3280., 2/  !  06
C +         |          |         |             |       |        |  |
C +         |          |         |             |       |        |  v_
C +         |          |         |             |       |        |  OUTPUT = 0: All    OUTPUT are    prohibited 
C +         |          |         |             |       |        |  OUTPUT = 1: netcdf OUTPUT decided in AWSloc
C +         |          |         |             |       |        |  OUTPUT = 2: netcdf OUTPUT decided in AWSloc
C +         |          |         |             |       |        |              ASCII  OUTPUT decided in AWSloc
C +         |          |         |             |       |        v                           (see #WV in SISVAT)
C +         |          |         |             |       v        ALTITUDE of the Station
C +         |          |         |             v       LONGITUDE         of the Station
C +         |          v         v             LATITUDE                  of the Station
C +         v          ATTRIBUTE of the Station, will be written in a title of the corresponding netcdf file
C +         LABEL of the Station,    will be used as the first 3 characters of the corresponding netcdf file


C +  *******
C +--CAUTION: DO'NT FORGET TO MODIFY THE parameter n_AWS in AWSloc       IF YOU ADD NEW STATIONS!
C +  *******                                                AWS_nc
C +                                                         AWS_nc_DATA 
 

      end
