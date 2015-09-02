      subroutine Debugg_MAR(debugm)

C +------------------------------------------------------------------------+
C |                                                                        |
C | MAR          Debugg_MAR.f-Profil                    Vd 05-11-2004  MAR |
C |              Verification of MAR Variables                             |
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

      include 'MAR_RA.inc'

      include 'MAR_LB.inc'
      include 'MAR_DY.inc'
      include 'MAR_HY.inc'
      include 'MAR_CA.inc'
      include 'MAR_TE.inc'
      include 'MAR_TU.inc'
      include 'MAR_SV.inc'
      include 'MAR_SL.inc'
      include 'MAR_TV.inc'


C +--Local   Variables
C +  =================

      logical          DUMPOK
      real             dump_v(mx,my,mz)

      character*10     debugm

      real             qsat0D

      logical          Debugg
      common/DebuggLOG/Debugg

      logical          WrFul1
      logical          WrFul2

      real             WVa
      real             Taq(mz),Ta2(mz),Ta1(mz),TCM(mz),CMx(mz)
      common/Debugg_CM/Taq    ,Ta2    ,Ta1    ,TCM    ,CMx
      real             Qsa(mz)
      integer          i_wri  ,j_wri  ,l      ,n

      real             VV_ijk ,VV_MAX ,TH_MIN ,TH_MAX ,TT_MIN ,TT_MAX
      real             QQ_MIN ,QQ_MAX ,TK_MIN ,TK_MAX ,Te_MIN ,Te_MAX
      real             pp_MIN ,pp_MAX ,zi_MIN ,zi_MAX ,HL_MIN ,HL_MAX
      integer          ip_MAX ,jp_MAX ,ip_MIN ,jp_MIN ,iz_MAX ,jz_MAX 
      integer          iz_MIN ,jz_MIN ,iL_MAX ,jL_MAX ,iL_MIN ,jL_MIN
      integer          iV_MAX ,jV_MAX ,kV_MAX 
      integer          iH_MAX ,jH_MAX ,kH_MAX ,iH_MIN ,jH_MIN ,kH_MIN 
      integer          iT_MAX ,jT_MAX ,kT_MAX ,iT_MIN ,jT_MIN ,kT_MIN 
      integer          iQ_MAX ,jQ_MAX ,kQ_MAX ,iQ_MIN ,jQ_MIN ,kQ_MIN 
      integer          iK_MAX ,jK_MAX ,kK_MAX ,iK_MIN ,jK_MIN ,kK_MIN 
      integer          ie_MAX ,je_MAX ,ke_MAX ,ie_MIN ,je_MIN ,ke_MIN 

      data     WrFul1/.true./  !  Grid Pt (i_wri,j_wri): General OUTPUT Switch
      data     WrFul2/.true./  !  Grid Pt (i_wri,j_wri): Surface OUTPUT Switch
      data     DUMPOK/.false./ !  DUMP                                  Switch

      data     Debugg/.false./ !  Auxiliary Variable ** PLEASE DO'NT MODIFY ** 

      data     i_wri/104/
      data     j_wri/20/


C +--VERIFICATION
C +  ============

       VV_MAX = 0.d0
       TH_MIN = 5.d3
       TH_MAX = 0.d0
       TT_MIN = 5.d3
       TT_MAX = 0.d0
       QQ_MIN = 5.d3
       QQ_MAX = 0.d0
       TK_MIN = 5.d3
       TK_MAX =-1.d0
       Te_MIN = 5.d3
       Te_MAX =-1.d0
       pp_MIN = 5.d3
       pp_MAX = 0.d0
       zi_MIN = 5.d3
       zi_MAX =-1.d0
       HL_MIN = 5.d6
       HL_MAX =-5.d6

       DO j=1,my
       DO i=1,mx
        IF (pstDYn(i,j)  .gt.pp_MAX)                                THEN
            pp_MAX =         pstDYn(i,j)
            ip_MAX = i
            jp_MAX = j
        END IF

        IF ( pstDY(i,j)  .lt.pp_MIN)                                THEN
            pp_MIN =         pstDYn(i,j)
            ip_MIN = i
            jp_MIN = j
        END IF

        IF (zi__TE(i,j)  .gt.zi_MAX)                                THEN
            zi_MAX =         zi__TE(i,j)
            iz_MAX = i
            jz_MAX = j
        END IF

        IF (zi__TE(i,j)  .lt.zi_MIN)                                THEN
            zi_MIN =         zi__TE(i,j)
            iz_MIN = i
            jz_MIN = j
        END IF

        IF (HlatSL(i,j)  .gt.HL_MAX)                                THEN
            HL_MAX =         HlatSL(i,j)
            iL_MAX = i
            jL_MAX = j
        END IF

        IF (HlatSL(i,j)  .lt.HL_MIN)                                THEN
            HL_MIN =         HlatSL(i,j)
            iL_MIN = i
            jL_MIN = j
        END IF

       END DO
       END DO

       DO k=1,mz
       DO j=1,my
       DO i=1,mx
            VV_ijk  =   sqrt(uairDY(i,j,k)*uairDY(i,j,k)
     .                      +vairDY(i,j,k)*vairDY(i,j,k))
        IF (VV_ijk       .gt.VV_MAX)                                THEN
            VV_MAX =         VV_ijk
            iV_MAX = i
            jV_MAX = j
            kV_MAX = k
        END IF

        IF (pktaDY(i,j,k).gt.TH_MAX)                                THEN
            TH_MAX =         pktaDY(i,j,k)
            iH_MAX = i
            jH_MAX = j
            kH_MAX = k
        END IF

        IF (pktaDY(i,j,k).lt.TH_MIN)                                THEN
            TH_MIN =         pktaDY(i,j,k)
            iH_MIN = i
            jH_MIN = j
            kH_MIN = k
        END IF

        IF (TairDY(i,j,k).gt.TT_MAX)                                THEN
            TT_MAX =         TairDY(i,j,k)
            iT_MAX = i
            jT_MAX = j
            kT_MAX = k
        END IF

        IF (TairDY(i,j,k).lt.TT_MIN)                                THEN
            TT_MIN =         TairDY(i,j,k)
            iT_MIN = i
            jT_MIN = j
            kT_MIN = k
        END IF

        IF (  qvDY(i,j,k).gt.QQ_MAX)                                THEN
            QQ_MAX =           qvDY(i,j,k)
            iQ_MAX = i
            jQ_MAX = j
            kQ_MAX = k
        END IF

        IF (  qvDY(i,j,k).lt.QQ_MIN)                                THEN
            QQ_MIN =           qvDY(i,j,k)
            iQ_MIN = i
            jQ_MIN = j
            kQ_MIN = k
        END IF

        IF (ect_TE(i,j,k).gt.TK_MAX)                                THEN
            TK_MAX =         ect_TE(i,j,k)
            iK_MAX = i
            jK_MAX = j
            kK_MAX = k
        END IF

        IF (ect_TE(i,j,k).lt.TK_MIN)                                THEN
            TK_MIN =         ect_TE(i,j,k)
            iK_MIN = i
            jK_MIN = j
            kK_MIN = k
        END IF

        IF (eps_TE(i,j,k).gt.Te_MAX)                                THEN
            Te_MAX =         eps_TE(i,j,k)
            ie_MAX = i
            je_MAX = j
            ke_MAX = k
        END IF

        IF (eps_TE(i,j,k).lt.Te_MIN)                                THEN
            Te_MIN =         eps_TE(i,j,k)
            ie_MIN = i
            je_MIN = j
            ke_MIN = k
        END IF

       END DO
       END DO
       END DO

       write(6,6000) iterun,debugm,jdarGE,jhurGE,minuGE,jsecGE
     .              ,iV_MAX,jV_MAX,kV_MAX,VV_MAX
     .              ,ip_MAX,jp_MAX,       pp_MAX
     .              ,ip_MIN,jp_MIN,       pp_MIN
     .              ,iH_MAX,jH_MAX,kH_MAX,TH_MAX*3.73
     .              ,iH_MIN,jH_MIN,kH_MIN,TH_MIN*3.73
     .              ,iT_MAX,jT_MAX,kT_MAX,TT_MAX
     .              ,iT_MIN,jT_MIN,kT_MIN,TT_MIN
     .              ,iQ_MAX,jQ_MAX,kQ_MAX,QQ_MAX*1.d3
     .              ,iQ_MIN,jQ_MIN,kQ_MIN,QQ_MIN*1.d3
     .              ,iL_MAX,jL_MAX,       HL_MAX
     .              ,iL_MIN,jL_MIN,       HL_MIN
     .              ,iK_MAX,jK_MAX,kK_MAX,TK_MAX
     .              ,iK_MIN,jK_MIN,kK_MIN,TK_MIN
     .              ,ie_MAX,je_MAX,ke_MAX,Te_MAX
     .              ,ie_MIN,je_MIN,ke_MIN,Te_MIN
     .              ,iz_MAX,jz_MAX,       zi_MAX
     .              ,iz_MIN,jz_MIN,       zi_MIN
 6000  format(/,i6,3x,a10,3x,4(i2,'.'),'  Vmax',3i3,f15.2
     .       ,/,9x,10('=')
     .       ,/,34x,      '  pmax',2i3,f18.2,8x,'pmin',2i3,f18.2
     .       ,/,34x,      '  THmx',3i3,f15.2,8x,'THmn',3i3,f15.2
     .       ,/,34x,      '  Tmax',3i3,f15.2,8x,'Tmin',3i3,f15.2
     .       ,/,34x,      '  Qmax',3i3,f15.3,8x,'Qmin',3i3,f15.3
     .       ,/,34x,      '  Lmax',2i3,f18.3,8x,'Lmin',2i3,f18.3
     .       ,/,34x,      '  TKmx',3i3,f15.3,8x,'TKmn',3i3,f15.3
     .       ,/,34x,      '  e_mx',3i3,f15.6,8x,'e_mn',3i3,f15.6
     .       ,/,34x,      '  ZImx',2i3,f18.2,8x,'ZImn',2i3,f18.2
     .       )

       IF (WrFul1)                                                  THEN
         i = i_wri
         j = j_wri

            WVa   =   0.
         DO n=1,nLimit
            WVa   = WVa+WV__SL(i,j,n)
         ENDDO
            WVa   = WVa/nLimit

         DO k=1,mz
            Ta2(k)=        Ta1    (k)
            Ta1(k)=        Taq    (k)
            Taq(k)=     pktaDY(i,j,k)*          pkDY(i,j,k)
            Qsa(k)=     qsat0D(Taq(k),sigma(k),pstDY(i,j),ptopDY,0)
         ENDDO
         IF (.NOT.Debugg)                                           THEN
                  Debugg=.true.
           DO k=1,mz
            Ta1(k)=        Taq    (k)
            Ta2(k)=        Taq    (k)
           ENDDO
         END IF
         DO k=1,mz
            TCM(k)= abs(0.5* ( Ta2(k)+Taq  (k) )-Ta1    (k))
            CMx(k)= max(       TCM(k),CMx  (k) )
         ENDDO

         write(6,6001)
 6001    format(/,8x,' Ta',5x,' Qa',5x,'Qsa',
     .            5x,' Qi',5x,' Qw',5x,' Qs',5x,' Qr',
     .            5x,' Ua',5x,' Va',5x,' Wa',
     .            5x,'TKE',5x,'eps',5x,' Kz',5x,'  z',
     .           ' T Comput.Mode')
         write(6,6002)(k, Taq(    k)     ,  qvDY(i,j,k)*1.e3,Qsa(k)*1.e3
     .                ,  qiHY(i,j,k)*1.e3,  qwHY(i,j,k)*1.e3
     .                ,  qsHY(i,j,k)*1.e3,  qrHY(i,j,k)*1.e3
     .                ,uairDY(i,j,k)     ,vairDY(i,j,k)
     .                ,wairDY(i,j,k)
     .                ,ect_TE(i,j,k),eps_TE(i,j,k),   TUkvm(i,j,k)
     .                ,gplvDY(i,j,k)*grvinv            *1.e-3
     .                ,   TCM(    k),   CMx(    k)
     .                ,           k =    1 ,mz)
 6002    format(i3,14f8.3,2f7.3)
         write(6,6003) TairSL(i,j)  ,qvapSL(i,j)
     .           ,1.e3*snowHY(i,j)  
     .           ,1.e3*rainHY(i,j)
     .           ,      SLuus(i,j)  * SLuus(i,j)
 6003    format(3x,2f8.3,24x,2f8.3,16x,f8.3)

        IF (WrFul2)                                                 THEN
         write(6,6100) i,j,isolSL(i,j), jdarGE,mmarGE,
     .                                  jhurGE,minuGE,jsecGE,
     .                                 (SLsrfl(i,j,n  ),n=1,mw)
         write(6,6101)                 (tsrfSL(i,j,n  ),n=1,mw)
         write(6,6102)                  sst_LB(i,j),
     .                                  sst1LB(i,j),
     .                                  sst2LB(i,j)
         write(6,6104)                  tairDY(i,j,mz)
         write(6,6105)                   ssvSL(i,j,mz)
         write(6,6106)                    qvDY(i,j,mz)
         write(6,6107)                  RAdsol(i,j)
         write(6,6108)                  RAD_ir(i,j)
         write(6,6109)                 (IRsoil(i,j,n  ),n=1,mw)
         write(6,6110)                 (SLlmol(i,j,n  ),n=1,mw)
         write(6,6111)                 (SLuusl(i,j,n  ),n=1,mw)
         write(6,6112)                 (SL_z0( i,j,n  ),n=1,mw)
         write(6,6113)                 (SLutsl(i,j,n  ),n=1,mw)
         write(6,6114)                 (SL_r0( i,j,n  ),n=1,mw)
         write(6,6115)                 (SLuqsl(i,j,n  ),n=1,mw),WVa
         write(6,6116)                  maskSL(i,j)
         write(6,6117)                  isolTV(i,j)
         write(6,6118)                  AlbSTV(i,j)
         write(6,6119)                 (ivegTV(i,j,n  ),n=1,mw)
         write(6,6120)                 (alaiTV(i,j,n  ),n=1,mw)
         write(6,6121)                 (glf_TV(i,j,n  ),n=1,mw)
         write(6,6122)                 (TvegTV(i,j,n  ),n=1,mw)
         write(6,6123)                 (CaSnTV(i,j,n  ),n=1,mw)
         write(6,6124)                 (CaWaTV(i,j,n  ),n=1,mw)
         write(6,6125)                 (psivTV(i,j,n  ),n=1,mw)
         write(6,6126)                ((TsolTV(i,j,n,l),l=1,llx),n=1,mw)
         write(6,6127)                ((eta_TV(i,j,n,l),l=1,llx),n=1,mw)
 6100    format(/,' SL characteristics, Grid Point',i4,i3,
     .          /,' Type   = ', i15  ,'   (Time = ',i2,'-',i2,
     .                                              i3,'h',i2,':',i2,')'
     .         ,/,' SrfSL  = ',9e15.6)
 6101    format(  ' Ts     = ',9e15.6)
 6102    format(  ' SST    = ', e15.6,'   (between',e15.6,
     .                                      '  and',e15.6,')')
 6104    format(  ' Details:',
     .          /,' tairDY = ', e15.6)
 6105    format(  '  ssvSL = ', e15.6)
 6106    format(  '   qvDY = ', e15.6)
 6107    format(  ' RAdsol = ', e15.6)
 6108    format(  ' RAD_ir = ', e15.6)
 6109    format(  ' IRsoil = ',9e15.6)
 6110    format(  ' SLlmol = ',9e15.6)
 6111    format(  ' SLuusl = ',9e15.6)
 6112    format(  ' SL_z0  = ',9e15.6)
 6113    format(  ' SLutsl = ',9e15.6)
 6114    format(  ' SL_r0  = ',9e15.6)
 6115    format(  ' SLuqsl = ',2e15.6,'   WVa = ',e15.6)
 6116    format(  ' maskSL = ', i15)
 6117    format(  ' isolTV = ', i15)
 6118    format(  ' AlbSTV = ', e15.6)
 6119    format(  ' ivegTV = ',2i15)
 6120    format(  ' alaiTV = ',9e15.6)
 6121    format(  ' glf_TV = ',9e15.6)
 6122    format(  ' TvegTV = ',9e15.6)
 6123    format(  ' CaSnTV = ',9e15.6)
 6124    format(  ' CaWaTV = ',9e15.6)
 6125    format(  ' psivTV = ',9e15.6)
 6126    format(  ' TsolTV = ',7e15.6,6(/,'          ',7e15.6))
 6127    format(  ' eta_TV = ',7e15.6,6(/,'          ',7e15.6))

        END IF

       END IF


       IF (DUMPOK)                                                  THEN
          DO k=1,mz
          DO j=1,my
          DO i=1,mx
           dump_v(i,j,k) = pktRAd(i,j,k)
          END DO
          END DO
          END DO

C +            **********
          call dump3D_MAR(dump_v,debugm,'pktRAd    ')
C +            **********

       END IF

      return
      end 
      subroutine dump3D_MAR(dump3D,debugm,dump_m)

C +------------------------------------------------------------------------+
C |                                                                        |
C | MAR          dump3D_MAR                             Mc 07-04-2004  MAR |
C |              DUMP for Verification of MAR Variables                    |
C |                                                                        |
C +------------------------------------------------------------------------+


      IMPLICIT NONE


C +--General Variables
C +  =================

      include 'MARCTR.inc'

      include 'MARdim.inc'
      include 'MARgrd.inc'

      real              dump3D(mx,my,mz)
      character*10      debugm,dump_m

      logical           dumpIN
      common/dump3D_log/dumpIN

      IF (.NOT.dumpIN)                                              THEN
        open(unit=80,status='unknown',name='dump3D_MAR.OUT')
        rewind    80
        dumpIN=.true.
      END IF

      DO k=1,mz
      DO j=1,my
      DO i=1,mx
        write(80,800) itexpe,debugm,dump_m,i,j,k,dump3D(i,j,k)
 800    format(i6,3x,a10,3x,a10,3i6,f15.6)
      ENDDO
      ENDDO
      ENDDO

      return
      end
