C   +-------------------------------------------------------------------+
C   |  Subroutine SOUNDg                             June 2002  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : - SNDfil : sounding file                                  |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: - SND_pt : vertical profile of potential temperature      |
C   | ^^^^^^^ - SND_qv :    "        "    "  specific humidity          |
C   |         - SND_ff :    "        "    "  wind velocity (m/s)        |
C   |         - SND_dd :    "        "    "  wind direction             |
C   |           Note that SND_dd is set for the standard trigonometric  |
C   |           system.                                                 |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE SOUNDg

      IMPLICIT NONE


C +---General variables
C +   =================

      include 'NSTdim.inc'
      include 'NESTOR.inc'
      include 'NSTvar.inc'
      include 'SNDvar.inc'
      include 'CTRvar.inc'


C +---Local variables
C +   ===============

      INTEGER i,k,ifin,zero,fID,fact,nlines,index

      REAL    g,ra,cp,delp,gam,dtint,qsat,degrad,cap,pps,ppm,pp,
     .        pp1,ppf,hh,dpsl,WK4_1D(mz),WK5_1D(mz),WK6_1D(mz),
     .        NST1sp,NST1_t(mz),NST1_p(mz+1),NST1_z(mz+1),NST1sh

      CHARACTER*3  empty


C +---DATA
C +   ====

      DATA g       /    9.81        /
      DATA cp      / 1004.0         /
      DATA ra      /  287.0         /
      DATA degrad  / 1.745329252d-2 /
      DATA cap     / 0.28586d0      /
      DATA zero    /    0           /


C +---Initialisation
C +   ==============

      lfirst_LSC = .true.
      lfirst_NST = .true.


C +---Screen message
C +   ==============

      write(6,*) 'Vertical sounding (used instead 3-D fields)'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*) ' '


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Reading of sounding file
C +   ======================== 

c #PR write(6,*) 'Open sounding file : ',SNDfil
c #PR write(6,*) ' Pressure   Temp.  Rel.H.   Wind  Direct.'
c #PR write(6,*) ' --------   -----  ------   ----  -------'

      OPEN (unit=15,status='old',file=SNDfil)
      REWIND 15

       read  (15,222) SNDlat,SNDlon,SNDjda,SNDmma,SNDjhu
222    format(f5.1,4x,f5.1,11x,3i3)

       read  (15,224) SND_pp(1),SND_tt(1),SND_hr(1),SND_ff(1),SND_dd(1)
c #PR  write ( 6,225) SND_pp(1),SND_tt(1),SND_hr(1),SND_ff(1),SND_dd(1)
224    format(f5.0,4f4.0)
225    format(1x,5f8.0)

       read  (15,226) ifin
226    format(i5)

       ifin=MIN(mzm-4,ifin)
 
       DO i=2,ifin
        read (15,224) SND_pp(i),SND_tt(i),SND_hr(i),SND_ff(i),SND_dd(i)
        SND_hr(i)=MAX(SND_hr(i),1.)
c #PR   write( 6,225) SND_pp(i),SND_tt(i),SND_hr(i),SND_ff(i),SND_dd(i)
       ENDDO

      CLOSE (unit=15)


C +---Completion of the sounding
C +   ==========================

      SND_zz(0)=-500.
      SND_zz(1)=   0.
      dtint    =-(g/cp)*(SND_zz(0)-SND_zz(1))
C +...Dry adiabatic Lapse Rate -(g/cp) is assumed

      SND_tt(0)=SND_tt(1) + dtint
      SND_pp(0)=SND_pp(1) * EXP( -g/ra/(273.15+SND_tt(1))
     .                           *(SND_zz(0)-SND_zz(1)) )
      SND_hr(0)=SND_hr(1)

      IF (ifin.lt.mzm-1) THEN
       delp=MAX(SND_pp(ifin)/5.0,SND_pp(ifin)-1.)
     .                            /real(mzm-ifin)
C +... Top Completed Sounding Pressure is 1/5 of SND_pp(ifin)
       DO i=ifin+1,mzm
        SND_pp(i)=SND_pp(i-1)-delp
        SND_tt(i)=SND_tt(ifin)
        SND_hr(i)=SND_hr(ifin)
        SND_ff(i)=SND_ff(ifin)
        SND_dd(i)=SND_dd(ifin)
       ENDDO
      ENDIF 


C +---Conversion into MKS Units System
C +   ================================

      DO i=0,mzm
       SND_tt(i)=SND_tt(i)+273.15
C +... degC -> K
      ENDDO

      DO i=0,mzm
       SND_ff(i)=SND_ff(i)/1.92308
C +... kn -> m/s         1.92308 kn = 1 m/s
       SND_dd(i)=90.-(SND_dd(i)*10.)
C +... Synoptic (RMI sounding) to trigonometric (MAR) convention
       IF (SND_dd(i).ge.180.) THEN
        SND_dd(i)=SND_dd(i)-180.
       ELSE
        SND_dd(i)=SND_dd(i)+180.
       ENDIF
C +... Direction from ---> towards
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Potential temperature
C +   =====================

      DO i=1,mzm
       SND_pt(i)=SND_tt(i)*(SND_pp(1)/SND_pp(i))**cap
      ENDDO


C +---Wind components
C +   ===============

      DO i=1,mzm
       SND_uu(i)=SND_ff(i)*COS(degrad*SND_dd(i))
       SND_vv(i)=SND_ff(i)*SIN(degrad*SND_dd(i))
      ENDDO


C +---Specific Humidity
C +   =================

      DO i=1,mzm
       SND_qv(i)=SND_hr(i)/100.*qsat(SND_tt(i),SND_pp(i)/10.)
       SND_qv(i)=MAX(0.000003,SND_qv(i))
      ENDDO


C +---Underground fields
C +   ==================

      SND_pt(0)=SND_pt(1)
      SND_qv(0)=SND_qv(1)
      SND_hr(0)=SND_hr(1)
      SND_uu(0)=SND_uu(1)
      SND_vv(0)=SND_vv(1)
      SND_ff(0)=0.
      SND_dd(0)=SND_dd(1)


C +---Integration of the Hydrostatic Relation (-> SND_zz)
C +   ===================================================

      SND_zz(1)=0.0

      gam=(SND_tt(1)-SND_tt(0))/(SND_pp(1)-SND_pp(0))
     .   *(SND_pp(1)+SND_pp(0))/2.d0
     .   * g / (ra*(SND_tt(1)+SND_tt(0))/2.d0)

      IF (gam.ne.0.) THEN
       SND_zz(0)=SND_zz(1)+(SND_tt(0)/gam)
     .          *((SND_pp(1)/SND_pp(0))**(ra*gam/g)-1.)
      ELSE
       SND_zz(0)=SND_zz(1)+(ra*SND_tt(0)/g)
     .          *log(SND_pp(1)/SND_pp(0))
      ENDIF 


      DO i=1,mzm-1
       gam=(SND_tt(i+1)-SND_tt(i))/(SND_pp(i+1)-SND_pp(i))
     .    *(SND_pp(i+1)+SND_pp(i))/2.
     .    * g / (ra*(SND_tt(i+1)+SND_tt(i))/2.)
       IF (gam.ne.0.) THEN
        SND_zz(i+1)=SND_zz(i)-(SND_tt(i)/gam)
     .             *((SND_pp(i+1)/SND_pp(i))**(ra*gam/g)-1.)
       ELSE
        SND_zz(i+1)=SND_zz(i)-(ra*SND_tt(i)/g)
     .             *log(SND_pp(i+1)/SND_pp(i))
       ENDIF 
      ENDDO


C +---Pressure hPa -> kPa
C +   ===================

      DO i=0,mzm
       SND_pp(i)=SND_pp(i)/10.
      ENDDO


C +---Computation of hybrid coordinates
C +   =================================

C +---Reference levels for hybrid coordinates
C +   ---------------------------------------

      pp1  = 105.       ! Reference pressure (KPa)
      dpsl = 20.        ! "> boundary layer" (KPa)


C +---Selection of vertical levels
C +   ----------------------------

      pps = SND_pp(1)
      ppm = pps - dpsl
      DO k = 0,mzm
       pp = SND_pp(k)
       hh = pp/pp1
       IF (pp.gt.ppm) THEN
        ppf= (pp-ppm)/(pps-ppm)
        hh = hh + (pp1-pps)/pp1 * ppf * ppf
       ENDIF
       SNDhyb(k) = LOG(hh)
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Vertical grid in the NST model (depend on SP)
C +   ==============================

      NST1sp=SND_pp(1)
      NST1sh=SND_zz(1)

C +        ******
      CALL VERgrd (empty ,NSTmod,fID,zero,mz,NST1sp,
     .             NST1sh,NST1_t,NST1_p,NST1_z,
     .             NSTgdz,WK4_1D,WK5_1D,WK6_1D)
C +        ******

      DO k=1,mz
       NST_hp(1,1,k)=NST1_z(k)
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Vertical interpolation
C +   ======================

      DO i=0,mzm
       TMPhyb(i)=SNDhyb(mzm-i)
       TMP_pt(i)=SND_pt(mzm-i)
       TMP_qv(i)=SND_qv(mzm-i)
       TMP_uu(i)=SND_uu(mzm-i)
       TMP_vv(i)=SND_vv(mzm-i)
      ENDDO

      DO k=1,mz
       CALL INTlin(TMPhyb,TMP_uu,mzm+1,NST_hp(1,1,k),NST__u(1,1,k))
       CALL INTlin(TMPhyb,TMP_vv,mzm+1,NST_hp(1,1,k),NST__v(1,1,k))
       CALL INTlin(TMPhyb,TMP_pt,mzm+1,NST_hp(1,1,k),NST_pt(1,1,k))
       CALL INTlin(TMPhyb,TMP_qv,mzm+1,NST_hp(1,1,k),NST_qv(1,1,k))
      ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Format output for sounding file ("ms" lines)
C +   ============================================

      IF (ifin.gt.(ms-1)) THEN
      
       fact   = (ifin - (ms-1-10)) / 10 + 1
       nlines = (ifin - (ms-1-10)) / fact
       DO k=ms-1,ms-nlines-1,-1
        index     = MIN(mzm,(ms-1-10) + (k-(ms-nlines))*fact)
        SND_pp(k) = SND_pp(index)
        SND_tt(k) = SND_tt(index)
        SND_hr(k) = SND_hr(index)
        SND_qv(k) = SND_qv(index)
        SND_ff(k) = SND_ff(index)
        SND_dd(k) = SND_dd(index)
       ENDDO
       SND_pp(ms) = SND_pp(mzm)
       SND_tt(ms) = SND_tt(mzm)
       SND_hr(ms) = SND_hr(mzm)
       SND_qv(ms) = SND_qv(mzm)
       SND_ff(ms) = SND_ff(mzm)
       SND_dd(ms) = SND_dd(mzm)

      ELSE

       delp=MAX(SND_pp(ifin)/5.0,SND_pp(ifin)-1.)
     .                             /real(ms-ifin)
C +... Top Completed Sounding Pressure is 1/5 of SND_pp(ifin)
       DO i=ifin+1,ms
        SND_pp(i)=SND_pp(i-1)-delp
        SND_tt(i)=SND_tt(ifin)
        SND_hr(i)=SND_hr(ifin)
        SND_qv(i)=SND_qv(ifin)
        SND_ff(i)=SND_ff(ifin)
        SND_dd(i)=SND_dd(ifin)
       ENDDO

      ENDIF


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
