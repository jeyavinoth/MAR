C   +-------------------------------------------------------------------+
C   |  Subroutine MARout                            April 2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : Interpolated LSC (large-scale) fields                     |
C   | ^^^^^^^              or sounding                                  |
C   |                                                                   |
C   | Output: Creation of MARdyn.DAT (initialization)                   |
C   | ^^^^^^^             MARsol.DAT (       "      )                   |
C   |                     MARlbc.DAT (bound. forcing)                   |
C   |                     MARubc.DAT (bound. forcing)                   |
C   |                     MARsic.DAT (bound. forcing / Sea-Ice)         |
C   |                     MARdom.dat (surf. characteristics)            |
C   |         Note that *.DAT file can be written according to ASCII    |
C   |         or Binary format, depending on ASCfor variable.           |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE MARout


      IMPLICIT NONE


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'MARvar.inc'
      INCLUDE 'SNDvar.inc'
      INCLUDE 'CTRvar.inc'


C +---Local variables
C +   ---------------

      INTEGER i,j,k,l,n,ifh,nbchar,losnd,jjmez,
     .        veg_1D(mx),iwf_1D(mx),svt_1D(mx,nvx),sfr_1D(mx,nvx),
     .        isol1D(mx),ii1,ii2,jj1,jj2,tmpveg,mmx,mmy,m1,m2,
     .        ic1,ic2,jc1,jc2,ip11,jp11,mx1,mx2,my1,my2


      REAL    pcap,WK2_1D(mz),zesnd,TMP_dd(0:mzm),compt,compt1
      
      REAL    sst1D(mx),dsa_1D(mx),lai_1D(mx,nvx),SH_1D(mx),
     .        glf_1D(mx,nvx),d1__1D(mx),ts__1D(mx,nvx,nsl),
     .        sw__1D(mx,nvx,nsl),compt2,z0__1D(mx,mw),
     .        r0__1D(mx,mw),ch0_1D(mx),rsur1D(mx),alb01D(mx),
     .        eps01D(mx)

      REAL    uairUB(mx,my,mzabso),vairUB(mx,my,mzabso)
      REAL    pktaUB(mx,my,mzabso)

      CHARACTER*7  cklon
      CHARACTER*10 NSTinfo

      LOGICAL NSTini,NSTfor,NSTend,Vfalse
      
      REAL           MARsig(mz)
      COMMON/cMARvgd/MARsig  
C            See MARvgd.f

C +---Thermodynamical Constants (Atmosphere)
C +   --------------------------------------

      DATA pcap  /   3.730037070d0/
C +...     pcap  = 100 ** (R / Cp)


C +---Data
C +   ----

      DATA Vfalse  / .false. /


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Dates
C +   -----

      NSTini=.false.
      NSTend=.false.
      NSTfor=.true. 

      IF (DATtim.eq.DATini) NSTini=.true.
      IF (DATtim.eq.DATfin) NSTend=.true.


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---READING OF GRID PARAMETERS IN MARgrd.ctr
C +   ========================================

      OPEN (unit=51,status='old',file='MARgrd.ctr')

       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) maptyp
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) GElon0
       read (51,*) imez
       read (51,*) GElat0
       read (51,*) jmez
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) dx
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) GEddxx
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) ptopDY
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) zmin
       read (51,*) aavu
       read (51,*) bbvu
       read (51,*) ccvu
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,'(l4)') vertic
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) sst_SL
       read (51,*) !- - - - - - - - - - - - - - - - - -

      CLOSE(unit=51)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---DATE
C +   ====


      itexpe = 0

C +        ******
      CALL DATcnv (RUNiyr,mmaDYN,jdaDYN,jhuDYN,DATtim,Vfalse)
C +        ******

      iyrDYN=RUNiyr

      IF (DATtim.eq.DATfin) THEN
       jdh_LB=0
      ELSE
       jdh_LB=DAT_dt
      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---PREPARATION OF VARIABLES TO BE WRITTEN
C +   ======================================


C +---Surface characteristics
C +   -----------------------

      DO j=1,my
      DO i=1,mx
       d1_SL (i,j)  =NST_d1(i,j)
       alb0SL(i,j)  =NSTalb(i,j)
       eps0SL(i,j)  =NSTeps(i,j)
       DO n=1,mw
        SL_z0 (i,j,n)=NST_z0(i,j)
        SL_r0 (i,j,n)=NST_r0(i,j)
       ENDDO
       rsurSL(i,j)  =NSTres(i,j)
       ch0SL (i,j)  =NSTch0(i,j)
       ro_SL (i,j)  =0.0
      ENDDO
      ENDDO


C +---Surface layer variables
C +   -----------------------

      DO j=1,my
      DO i=1,mx
       tairSL(i,j)  =NST_st(i,j)
       t2_SL (i,j)  =NSTdst(i,j)
       DO n=1,mw
        tsrfSL(i,j,n)=NST_st(i,j)   ! Bloc Temporaire, a modif 
        SLsrfl(i,j,n)=0.
        SLuusl(i,j,n)=0.
        SLutsl(i,j,n)=0.
       ENDDO
       nSLsrf(i,j)  =1
       SLsrfl(i,j,1)=1.
C +
       qvapSL(i,j)  =1.e-5
       w2_SL (i,j)  =0.15 
       wg_SL (i,j)  =0.15 
       roseSL(i,j)  =0.
       hmelSL(i,j)  =0.
       hsnoSL(i,j)  =0.
       SaltSL(i,j)  =0.
      ENDDO
      ENDDO


C +---Prognostic variables
C +   --------------------

      DO j=1,my
      DO i=1,mx
       DO k=1,mz
        NST_qv(i,j,k)=MAX(1.e-5,NST_qv(i,j,k))
        NSTtmp(i,j,k)=NST_pt(i,j,k)/pcap
       ENDDO
       NSTtmp(i,j,mz+1)=NST_pt(i,j,mz)/pcap
       pstDY (i,j)     =NST_sp(i,j)-ptopDY
      ENDDO
      ENDDO

C +...uairDY <-- NST__u
C +...vairDY <-- NST__v
C +...qvDY   <-- NST_qv
C +...pktaDY <-- NSTtmp


C +---Boundary variables
C +   ------------------

      IF (NSTmod.ne.'CPL') THEN

       DO k=1,mzabso
        DO j=1,my
         DO i=1,mx
           uairUB(i,j,k) = NST__u (i,j,k)
           vairUB(i,j,k) = NST__v (i,j,k)
           pktaUB(i,j,k) = NSTtmp (i,j,k)
         ENDDO
        ENDDO
       ENDDO
       DO k=1,mz
        DO j=1,my
         DO i=1,n7
          vaxgLB(i,j,k,1) = NST__u (i,j,k)
          vaxgLB(i,j,k,2) = NST__v (i,j,k)
          vaxgLB(i,j,k,3) = NST_qv (i,j,k)
          vaxgLB(i,j,k,4) = NSTtmp (i,j,k)
          vaxgLB(i,j,1,5) = pstDY  (i,j)
          vaxgLB(i,j,mz,5)= tsrfSL (i,j,1)
         ENDDO
         DO i=mx-n6,mx
          vaxdLB(i,j,k,1) = NST__u (i,j,k)
          vaxdLB(i,j,k,2) = NST__v (i,j,k)
          vaxdLB(i,j,k,3) = NST_qv (i,j,k)
          vaxdLB(i,j,k,4) = NSTtmp (i,j,k)
          vaxdLB(i,j,1,5) = pstDY  (i,j)
          vaxdLB(i,j,mz,5)= tsrfSL (i,j,1)
         ENDDO
        ENDDO
        DO i=1,mx
         DO j=1,n7
          vayiLB(i,j,k,1) = NST__u (i,j,k)
          vayiLB(i,j,k,2) = NST__v (i,j,k)
          vayiLB(i,j,k,3) = NST_qv (i,j,k)
          vayiLB(i,j,k,4) = NSTtmp (i,j,k)
          vayiLB(i,j,1,5) = pstDY  (i,j)
          vayiLB(i,j,mz,5)= tsrfSL (i,j,1)
         ENDDO
         DO j=my-n6,my
          vaysLB(i,j,k,1) = NST__u (i,j,k)
          vaysLB(i,j,k,2) = NST__v (i,j,k)
          vaysLB(i,j,k,3) = NST_qv (i,j,k)
          vaysLB(i,j,k,4) = NSTtmp (i,j,k)
          vaysLB(i,j,1,5) = pstDY  (i,j)
          vaysLB(i,j,mz,5)= tsrfSL (i,j,1)
         ENDDO
        ENDDO
       ENDDO

      ENDIF


C +---Sounding variables
C +   ------------------

      DO k=0,mzm
       TMP_dd(k)=90.-SND_dd(k)
       IF (TMP_dd(k).lt.  0.) TMP_dd(k)=TMP_dd(k)+360.
       IF (TMP_dd(k).gt.360.) TMP_dd(k)=TMP_dd(k)-360.
      ENDDO


C +---Soil variables
C +   --------------

      DO j=1,my
      DO i=1,mx
       
                                 isolSL(i,j)=NSTsol(i,j)
c      IF (region.eq."AFW")                                         THEN
                                 isolTV(i,j)=NSTtex(i,j)
c      ELSE
c          IF (NSTtex(i,j).eq.1) isolTV(i,j)=2   ! loamy sand
c          IF (NSTtex(i,j).eq.2) isolTV(i,j)=5   ! sand
c          IF (NSTtex(i,j).eq.3) isolTV(i,j)=11  ! clay
c      ENDIF
        
       IF (region.eq."GRD".or.region.eq."ANT")                      THEN
           IF (NSTsol(i,j).le.2) isolTV(i,j)=0
           IF (NSTsol(i,j).eq.4) isolTV(i,j)=11
           IF (NSTsol(i,j).eq.3) isolTV(i,j)=12
C +...     Transform to SVAT (De Ridder) classification

       ENDIF
       
      ENDDO
      ENDDO


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Some constants specific to MAR
C +   ==============================


C +---Deardorff Soil Model Parameters
C +   -------------------------------

      cs2SL  = 86400.0
      w20SL  = 0.15
      wg0SL  = 0.10
      wk0SL  = 0.15
      wx0SL  = 0.20


C +---Typical Roughness Lengths (m) for land, sea, snow
C +   -------------------------------------------------

      zl_SL  = 1.00e-1
      zs_SL  = 1.00e-3
      zn_SL  = 1.00e-4


C +---Inversion surface temperature
C +   -----------------------------

      dtagSL = 0.


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Parameters of the vertical grid
C +   ===============================


C +        ******
      CALL SETsig (mz,zmin,aavu,bbvu,ccvu,ptopDY)
C +        ******

C +        ******
C      CALL GRDsig(mz,zmin,aavu,bbvu,ccvu,vertic,
C     .               sst_SL,TUkhmx,sigma,WK2_1D)
C +        ******
            
      DO k=1,mz
        sigma(k)=MARsig(k)
      ENDDO


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Specifications for horizontal grid and time step
C +   ================================================

      dx = 1000. *dx
      dy = dx
      dt = 4.e-3 *dx

      IF (NSTmod.eq.'M2D') THEN
       mmx   = mx
       mmy   = 1
       ii1   = 1
       ii2   = mx
       jj1   = jmez
       jj2   = jmez
       jjmez = 1
      ELSE
       IF (NSTmod.eq.'CPL') THEN
        mmx   = 1
        mmy   = 1
        ii1   = 2
        ii2   = 2
        jj1   = 2
        jj2   = 2
        jjmez = 1
       ELSE
        mmx   = mx
        mmy   = my
        ii1   = 1
        ii2   = mx
        jj1   = 1
        jj2   = my
        jjmez = jmez
       ENDIF
      ENDIF

      IF (SNDing) THEN
       imez = 1
       jmez = 1
       iyrDYN = SNDyar
       mmaDYN = SNDmma
       jdaDYN = SNDjda
       jhuDYN = SNDjhu
      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---1-D Topography
C +   ==============

      IF (NSTmod.eq.'M2D'.and.mmx.gt.1) THEN

       ic1 = MIN(2,mx)
       ic2 = MAX(1,mx-1)

       jc1 = MIN(2,my)
       jc2 = MAX(1,my-1)

       DO i=1,mx
 
        compt1    = 0.
        compt2    = 0.
        SH_1D(i)  = 0.
        isol1D(i) = 1
 
        DO j=jc1,jc2
         compt1   = compt1   + 1.
         SH_1D(i) = SH_1D(i) + NST_sh(i,j)
         IF (NSTsol(i,j).ge.3) THEN
          compt2  = compt2   + 1.
         ENDIF
        ENDDO
 
        IF (compt1.ge.1.) THEN
         SH_1D(i) = SH_1D(i) / compt1
        ENDIF

        IF (compt2.ge.(my/2)) isol1D(i) = 4
 
        IF (isol1D(i).le.2) THEN
         SH_1D(i) = 0.
        ENDIF

       ENDDO


C +....Topography filtering
C +    --------------------

       IF (TOPfilt) THEN

C +...  First filtering
        DO i=ic1,ic2
         IF (isol1D(i).ge.3) THEN
          SH_1D(i) = (SH_1D(i-2)+SH_1D(i-1)+2.*SH_1D(i)
     .               +SH_1D(i+1)+SH_1D(i+2)) / 6.0
         ENDIF
        ENDDO

C +...  Second filtering
        DO i=ic2,ic1,-1
         IF (isol1D(i).ge.3) THEN
          SH_1D(i) = (SH_1D(i-2)+SH_1D(i-1)+2.*SH_1D(i)
     .               +SH_1D(i+1)+SH_1D(i+2)) / 6.0
         ENDIF
        ENDDO

C +...  Third filtering
        DO i=ic1,ic2
         IF (isol1D(i).ge.3) THEN
          SH_1D(i) = (SH_1D(i-2)+SH_1D(i-1)+2.*SH_1D(i)
     .               +SH_1D(i+1)+SH_1D(i+2)) / 6.0
         ENDIF
        ENDDO

C +...  Fourth filtering
        DO i=ic2,ic1,-1
         IF (isol1D(i).ge.3) THEN
          SH_1D(i) = (SH_1D(i-2)+SH_1D(i-1)+2.*SH_1D(i)
     .               +SH_1D(i+1)+SH_1D(i+2)) / 6.0
         ENDIF
        ENDDO

C +...  Fifth filtering
        DO i=ic1,ic2
         IF (isol1D(i).ge.3) THEN
          SH_1D(i) = (SH_1D(i-1)+2.*SH_1D(i)+SH_1D(i+1)) / 4.0
         ENDIF
        ENDDO

C +...  Sixth filtering
        DO i=ic2,ic1,-1
         IF (isol1D(i).ge.3) THEN
          SH_1D(i) = (SH_1D(i-1)+2.*SH_1D(i)+SH_1D(i+1)) / 4.0
         ENDIF
        ENDDO

       ENDIF


       m1 = MIN(mx,n10)
       DO i=1,m1-1
        SH_1D(i)=SH_1D(m1)
       ENDDO

       m2 = MAX(1,mx-n10+1)
       DO i=m2+1,mx
        SH_1D(i)=SH_1D(m2)
       ENDDO


      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---1-D SST
C +   =======

      IF (NSTmod.eq.'M2D') THEN

       DO i=1,mx

        compt    = 0.
        sst1D(i) = 0.
 
        DO j=1,my
         IF (NSTsol(i,j).le.2) THEN
          IF (NSTsst(i,j).lt.1.) THEN
           compt    = compt    + 1.
           sst1D(i) = sst1D(i) + NST_st(i,j)
          ELSE
           compt    = compt    + 1.
           sst1D(i) = sst1D(i) + NSTsst(i,j)
          ENDIF
         ENDIF
        ENDDO
 
        IF (compt.ge.1.) THEN
         sst1D(i) = sst1D(i) / compt
        ENDIF
 
        IF (isol1D(i).ge.3) THEN
         sst1D(i) = 0.
        ENDIF
 
       ENDDO
 
      ENDIF
 

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---1-D Surface
C +   ===========

      IF (NSTini.and.NSTmod.eq.'M2D'.and.
     .    LoutDA.and.SELECT.eq.1        ) THEN


       DO i=1,mx

        ch0_1D(i) = 0.
        rsur1D(i) = 0.
        alb01D(i) = 0.
        eps01D(i) = 0.
        d1__1D(i) = 0.
        DO k=1,mw
         z0__1D(i,k) = 0.
         r0__1D(i,k) = 0.
        ENDDO

        DO j=1,my
         ch0_1D(i) = ch0_1D(i) + NSTch0(i,j)
         rsur1D(i) = rsur1D(i) + NSTres(i,j)
         alb01D(i) = alb01D(i) + NSTalb(i,j)
         eps01D(i) = eps01D(i) + NSTeps(i,j)
         d1__1D(i) = d1__1D(i) + NST_d1(i,j)
         DO k=1,mw
          z0__1D(i,k) = z0__1D(i,k) + NST_z0(i,j)
          r0__1D(i,k) = r0__1D(i,k) + NST_r0(i,j)
         ENDDO
        ENDDO

        ch0_1D(i) = ch0_1D(i) / REAL(my)
        rsur1D(i) = rsur1D(i) / REAL(my)
        alb01D(i) = alb01D(i) / REAL(my)
        eps01D(i) = eps01D(i) / REAL(my)
        d1__1D(i) = d1__1D(i) / REAL(my)
        DO k=1,mw
         z0__1D(i,k) = z0__1D(i,k) / REAL(my)
         r0__1D(i,k) = r0__1D(i,k) / REAL(my)
        ENDDO

       ENDDO


      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---1-D SISVAT variables
C +   ====================

      IF (NSTini.and.SVTmod.and.
     .    LoutDA.and.SELECT.eq.1.and.
     .    NSTmod.eq.'M2D'           ) THEN


       DO i=1,mx

        compt1    = 0.

        veg_1D(i) = 0
        tmpveg    = 0
        iwf_1D(i) = 0
        dsa_1D(i) = 0.
        DO k=1,nvx
         svt_1D(i,k) = 0
         sfr_1D(i,k) = 0
         lai_1D(i,k) = 0.
         DO l=1,nsl
          ts__1D(i,k,l) = 0.
          sw__1D(i,k,l) = 0.
         ENDDO
        ENDDO

        DO j=1,my
         IF (NSTsol(i,j).ge.3) THEN
          compt1    = compt1    + 1.
          tmpveg    = tmpveg    + NSTtex(i,j)
          iwf_1D(i) = iwf_1D(i) + NSTiwf(i,j)
          dsa_1D(i) = dsa_1D(i) + NSTdsa(i,j)
          DO k=1,nvx
           svt_1D(i,k) = svt_1D(i,k) + NSTsvt(i,j,k)
           sfr_1D(i,k) = sfr_1D(i,k) + NSTsfr(i,j,k)
           lai_1D(i,k) = lai_1D(i,k) + NSTlai(i,j,k)
           DO l=1,nsl
            ts__1D(i,k,l) = ts__1D(i,k,l) + NST_ts(i,j,k,l)
            sw__1D(i,k,l) = sw__1D(i,k,l) + NST_sw(i,j,k,l)
           ENDDO
          ENDDO
         ENDIF
        ENDDO
 
        IF (compt1.ge.1.) THEN
         tmpveg    = NINT (REAL(tmpveg)    / compt1)
         iwf_1D(i) = NINT (REAL(iwf_1D(i)) / compt1)
         dsa_1D(i) = dsa_1D(i) / compt1
         veg_1D(i) = NINT (REAL(veg_1D(i)) / compt1)
         IF (tmpveg.eq.1) veg_1D(i)=2   ! loamy sand
         IF (tmpveg.eq.2) veg_1D(i)=5   ! sand
         IF (tmpveg.eq.3) veg_1D(i)=11  ! clay
         DO k=1,nvx
          svt_1D(i,k) = NINT (REAL(svt_1D(i,k)) / compt1)
          sfr_1D(i,k) = NINT (REAL(sfr_1D(i,k)) / compt1)
          lai_1D(i,k) = lai_1D(i,k) / compt1
         ENDDO
         DO l=1,nsl
         DO k=1,nvx
          ts__1D(i,k,l) = ts__1D(i,k,l) / compt1
          sw__1D(i,k,l) = sw__1D(i,k,l) / compt1
         ENDDO
         ENDDO
        ENDIF

       ENDDO

      ENDIF


      IF (NSTfor.and.SVTmod.and.
     .    LoutDA.and.SELECT.eq.1.and.
     .    NSTmod.eq.'M2D'           ) THEN

       DO i=1,mx

        compt2    = 0.
        DO k=1,nvx
         glf_1D(i,k) = 0.
        ENDDO

        DO j=1,my
         IF (NSTsol(i,j).ge.4) THEN
          DO k=1,nvx
           compt2      = compt2      + 1.
           glf_1D(i,k) = glf_1D(i,k) + NSTlai(i,j,k)*NSTglf(i,j,k)
          ENDDO
         ENDIF
        ENDDO

        IF (compt2.ge.1.) THEN
         DO k=1,nvx
          IF (lai_1D(i,k).gt.0.) THEN
           glf_1D(i,k) = glf_1D(i,k) / compt2 / lai_1D(i,k)
           glf_1D(i,k) = MIN(1.0,glf_1D(i,k))
          ELSE
           glf_1D(i,k) = 0.0
          ENDIF
         ENDDO
        ENDIF

       ENDDO

      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Filter parameters
C +   -----------------

      CALL MARfil(my,dx,dt,FIslot,FIslou,FIslop,
     .                     FIkhmn,TUkhff,TUkhmx)

C     Note (PhM): we give the opportunity to change FIslo* here
C       because standard value is high in comparison to
C       recomendations in 
C       Raymond and Garder, MWR 116, Jan 1988, p209 
C       (suggests 0.0075, while default is 0.05 in MARfil)
C       Note that we do not change  FIlkhm, which
C       is computed in MARfil and used in MAR:TURhor_dyn.f
C       (i.e.: I don't know the reason for changing it
C       with the filter; of course it also smooth horizontal
C       fields, but may be physically based (?) in contrast to
C       the filter, which should only eliminates 2dx)

      IF (NSTfis.GE.0.0001) THEN
         FIslop= NSTfis
         FIslou= FIslop      
         FIslot= FIslop
       ENDIF

      write(6,*) 'Write  files :'
      write(6,*) '~~~~~~~~~~~~~~'


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output directory
C +   ================

      nbchar=1

      DO i=1,60
       IF (NSTdir(i:i).ne.' ') nbchar=i
      ENDDO


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---MAR include file : MARdim.inc
C +   =============================

      IF (NSTini.and.SELECT.eq.1) THEN

       open  (1,status='unknown',file=NSTdir(1:nbchar)
     .                                //'MARdim.inc')
       open  (2,status='unknown',file=NSTdir(1:nbchar)
     .                                //'MARdim.inc_old')
       rewind 1
       rewind 2

       ip11 = MIN(2,mmx)
       jp11 = MIN(2,mmy)
       mx1  = MAX(1,mmx-1)
       my1  = MAX(1,mmy-1)
       mx2  = MAX(1,mmx-2)
       my2  = MAX(1,mmy-2)

       IF (vector) THEN
        cklon = 'mx2*my2'
       ELSE
        cklon = '      1'
       ENDIF

       write(1,300) mmx,mmy,ip11,jp11,mz,mx1,mx2,my1,my2,mzabso,cklon,
     .              n6,n7,mw
       write(2,400) mmx,mmy,ip11,jp11,mz,mx1,mx2,my1,my2,mzabso,cklon

300    format('      integer   mx    ,my   ,ip11  ,jp11',/,
     .        '      parameter(mx=',i4,',my=',i4,
     .                       ',ip11=',i3,',jp11=',i3,')',/,
     .        '      integer   mz   ,mzir1     ,mzir',/,
     .        '      parameter(mz=',i4,',mzir1=mz+1,mzir=mz+2)',/,
     .        'c ... mzir1 may be chosen much larger than mz, ',/,
     .        'c     if the model vertical domain covers a small',/,
     .        'c     part of the air column',/,
     .        'c     ',/,
     .        '      integer   mx1     ,mx2',/,
     .        '      parameter(mx1=',i4,',mx2=',i4,')',/,
     .        '      integer   my1     ,my2     ,myd2',/,
     .        '      parameter(my1=',i4,',my2=',i4,',myd2=1+my/2)',/,
     .        '      integer   mz1     ,mzz',/,
     .        '      parameter(mz1=mz-1,mzz=mz+1)',/,
     .        '      integer   i_2',/,
     .        '      parameter(i_2=mx-mx1+1) ',/,
     .        '      integer   j_2',/,
     .        '      parameter(j_2=my-my1+1) ',/,
     .        '      integer   mzabso      ,mzhyd',/,
     .        '      parameter(mzabso =  ',i2,',mzhyd=mzabso+1)',/,
     .        'c     ',/,
     .        '      integer   klon,        klev',/,
     .        '      parameter(klon=',a7,',klev=mz)',/,
     .        'C +...if #NV removed (NO vectorization)',/,
     .        'C +   then      klon=      1',/,
     .        'C +   ',/,
     .        '      integer   kdlon,       kflev',/,
     .        '      parameter(kdlon=klon  ,kflev=klev)',/,
     .        'C +   ',/,
     .        '      integer   n6   ,n7',/,
     .        '      parameter(n6=',i2,',n7=',i2,')',/,
     .        'C +.. n6 et n7 determine a relaxation zone',
     .               'towards lateral boundaries',/,
     .        'C +   (large scale values of the variables).',/,
     .        'C +   This zone extends over n6-1 points.',/,
     .        'C +   Davies (1976) propose 5 points ',
     .               '(i.e. n6=6 and n7=7)',/,
     .        'C +   ',/,
     .        '      integer   mw',/,
     .        '      parameter(mw=',i3,')',/,
     .        'C +..           mw is the total number of mosaics',/,
     .        'C +   ')

400    format('      integer   mx    ,my   ,ip11  ,jp11',/,
     .        '      parameter(mx=',i4,',my=',i4,
     .                       ',ip11=',i3,',jp11=',i3,')',/,
     .        '      integer   mz   ,mzir1     ,mzir',/,
     .        '      parameter(mz=',i4,',mzir1=mz+1,mzir=mz+2)',/,
     .        'c ... mzir1 may be chosen much larger than mz, ',/,
     .        'c     if the model vertical domain covers a small',/,
     .        'c     part of the air column',/,
     .        'c     ',/,
     .        '      integer   mx1     ,mx2',/,
     .        '      parameter(mx1=',i4,',mx2=',i4,')',/,
     .        '      integer   my1     ,my2     ,myd2',/,
     .        '      parameter(my1=',i4,',my2=',i4,',myd2=1+my/2)',/,
     .        '      integer   mz1     ,mzz',/,
     .        '      parameter(mz1=mz-1,mzz=mz+1)',/,
     .        '      integer   i_2',/,
     .        '      parameter(i_2=mx-mx1+1) ',/,
     .        '      integer   j_2',/,
     .        '      parameter(j_2=my-my1+1) ',/,
     .        '      integer   mzabso      ,mzhyd',/,
     .        '      parameter(mzabso =  ',i2,',mzhyd=mzabso+1)',/,
     .        'c     ',/,
     .        '      integer   klon,        klev',/,
     .        '      parameter(klon=',a7,',klev=mz)',/,
     .        'C +...if #NV removed (NO vectorization)',/,
     .        'C +   then      klon=      1',/,
     .        'C +   ',/,
     .        '      integer   kdlon,       kflev',/,
     .        '      parameter(kdlon=klon  ,kflev=klev)' )

       close (1)
       close (2)

       write(6,*) 'MAR include file      MARdim.inc created'

      ENDIF 


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---MAR include file : MAR_TV.inc
C +   =============================

      IF (NSTini.and.SELECT.eq.1) THEN

       open  (1,status='unknown',file=NSTdir(1:nbchar)
     .                                //'MAR_TV.inc_old')
       rewind 1

       write (1,411) nvx,nsl
411    format('      integer    nvx  ,llx  ,iptx',/,
     .        '      parameter (nvx=',i3,',llx=',i3,',iptx=5)')

       write (1,412)
412    format(
     .  'C +',/,
     .  6x,'integer    imx   ,jmx',/,
     .  6x,'parameter (imx=mx,jmx=my)',/,
     .  'C +',/,
     .  6x,'real      deptTV(0:llx)',/,
     .  'C +...',10x,'deptTV: Soil Level Depth',/,
     .  'C +',/,
     .  6x,'real      dep2TV(0:llx)',/,
     .  'C +...',10x,'dep2TV: Soil Layer Depth',/,
     .  'C +',/,
     .  6x,'real      slopTV(imx,jmx)',/,
     .  'C +...',10x,'slopTV: Surface Slope',/,
     .  'C +',/,
     .  6x,'real      AlbSTV(imx,jmx)',/,
     .  'C +...',10x,'AlbSTV: Dry Soil Albedo',/,
     .  'C +',/,
     .  6x,'real      alaiTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'alaiTV: Leaf Area Index',/,
     .  'C +',/,
     .  6x,'real      glf_TV(imx,jmx,nvx)',/,
     .  'C +...',10x,'glf_TV: Green Leaf Fraction',/,
     .  'C +',/,
     .  6x,'real      CaWaTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'CaWaTV: Canopy Intercepted Water Content',/,
     .  'C +',/,
     .  6x,'real      CaSnTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'CaSnTV: Canopy Intercepted Snow Content',/,
     .  'C +',/,
     .  6x,'real      TvegTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'TvegTV: Skin Vegetation Temperature',/,
     .  'C +',/,
     .  6x,'real      TgrdTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'TgrdTV: Skin Soil Temperature',/,
     .  'C +',/,
     .  6x,'real      TsolTV(imx,jmx,nvx,llx)',/,
     .  'C +...',10x,'TsolTV: Layer Soil Temperature',/,
     .  'C +',/,
     .  6x,'real      eta_TV(imx,jmx,nvx,llx)',/,
     .  'C +...',10x,'eta_TV: Soil Moisture Content',/,
     .  'C +',/,
     .  6x,'real      psigTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'psigTV: Soil Hydraulic Potential',/,
     .  'C +',/,
     .  6x,'real      psivTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'psivTV: Vegetation Hydraulic Potential',/,
     .  'C +',/,
     .  6x,'real      runoTV(imx,jmx)',/,
     .  'C +...',10x,'runoTV: Time Integrated (Sub)surface Flow',/,
     .  'C +',/,
     .  6x,'real      draiTV(imx,jmx)',/,
     .  'C +...',10x,'draiTV: Time Integrated Drainage Flow',/,
     .  'C +',/,
     .  6x,'integer   iWaFTV(imx,jmx)',/,
     .  'C +...',9x,'(iWaFTV=0 ==> no Water Flux;',/,
     .  'C +   ',10x,'iWaFTV=1 ==> free drainage)',/,
     .  'C +',/,
     .  6x,'integer   ivegTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'ivegTV: Vegetation Type Index',/,
     .  'C +',/,
     .  6x,'integer   isolTV(imx,jmx)',/,
     .  'C +...',10x,'isolTV: Soil Type Index',/,
     .  'C +',/,
     .  6x,'integer   ifraTV(imx,jmx,nvx)',/,
     .  'C +...',10x,'ifraTV: Vegetation Class Coverage',/,
     .  'C +   ',10x,'        (3 Class, Last One is Open Water)',/,
     .  'C +',/,
     .  6x,'integer   IOi_TV(iptx),IOj_TV(iptx)',/,
     .  'C +...',10x,'IO Grid Indices',/,
     .  'C +',/,
     .  6x,'integer   itx,ivg',/,
     .  'C +',/,
     .  5x,' common/rsvaTV/AlbSTV,alaiTV,glf_TV,CaWaTV,CaSnTV,',/,
     .  5x,'.              runoTV,draiTV,TvegTV,TgrdTV,TsolTV,',/,
     .  5x,'.              eta_TV,psigTV,psivTV,deptTV,dep2TV ',/,
     .  5x,' common/isvaTV/iWaFTV,ivegTV,isolTV,ifraTV,IOi_TV,',/,
     .  5x,'.              IOj_TV,itx   ,ivg')

       close (1)

       write(6,*) 'MAR include file      MAR_TV.inc created'

      ENDIF 


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---MAR include file : MAR_SV.inc
C +   =============================

      IF (NSTini.and.SELECT.eq.1) THEN

       open  (1,status='unknown',file=NSTdir(1:nbchar)//'MAR_SV.inc')
       open  (2,status='unknown',file=NSTdir(1:nbchar)//'MAR_SV.inc_nv')
       rewind 1
       rewind 2

       write (1,410) nsl-1,nsno,nvx*5
       write (2,409) nsl-1,nsno,nvx*5
410    format('      integer   klonv    ,nsol  ,nsno',/,
     .        '      parameter(klonv=256,nsol=',i3,',nsno=',i4,')',/
     .        '      integer   nb_wri',/,
     .        '      parameter(nb_wri=',i3,')',/)
409    format('      integer   klonv    ,nsol  ,nsno',/,
     .        '      parameter(klonv=  1,nsol=',i3,',nsno=',i4,')',/
     .        '      integer   nb_wri',/,
     .        '      parameter(nb_wri=',i3,')',/)

       close (1)
       close (2)

       write(6,*) 'MAR include file      MAR_SV.inc created'

      ENDIF 


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---MAR include file : MAR_LB.inc
C +   =============================

      IF (NSTini.and.SELECT.eq.1) THEN

       open  (1,status='unknown',file=NSTdir(1:nbchar)
     .                                //'MAR_LB.inc_old')
       rewind 1

       write (1,420) n6,n7

420    format('      integer   n6  ,n7',/,
     .        '      parameter(n6=',i3,',n7=',i3,')',/,
     .        'C +...n6 et n7 determine a relaxation zone towards ',/,
     .        'C +   lateral boundaries ',/,
     .        'C +   (large scale values of the variables). ',/,
     .        'C +   This zone extends over n6-1 points. ',/,
     .        'C +   Davies (1976) propose 5 points ',/,
     .        'C +   (i.e. n6=6 and n7=7)',/,
     .        'C +   ',/,
     .        '      integer iyr_LB,mma_LB,jda_LB,jhu_LB,jdh_LB',/,
     .        '      common/nudite/iyr_LB,mma_LB,jda_LB,jhu_LB,jdh_LB',/,
     .        'C +...      iyr_LB: Year',/,
     .        'C +         mma_LB: Month',/,
     .        'C +         jda_LB: Day',/,
     .        'C +         jhu_LB: Hour  (UT)',/,
     .        'C +         jdh_LB: Time Interval before next ',/,
     .        'C +                 GCM/NWP LBC (hour)',/,
     .        'C +         jdh_LB=0 ==> NO further GCM/NWP ',/,
     .        'C +                      LBC available',/,
     .        '      integer       tim1LB,tim2LB',/,
     .        '      common/nudtim/tim1LB,tim2LB',/,
     .        'C +...      tim1LB: Time of the previous LBC (second)',/,
     .        'C +         tim2LB: Time of the next     LBC (second)',/,
     .        'C +   ',/,
     .        '      integer n40xLB, n50xLB, n5mxLB, n6mxLB, n7mxLB,',/,
     .        '     .        n40yLB, n50yLB, n5myLB, n6myLB, n7myLB ',/,
     .        '      common/nudind/n40xLB, n50xLB, n5mxLB, n6mxLB,  ',/,
     .        '     .              n7mxLB, n40yLB, n50yLB, n5myLB,  ',/,
     .        '     .              n6myLB, n7myLB ',/,
     .        'C +...           ...n6mxLB, n7mxLB, n6myLB, n7myLB...',/,
     .        'C +                 define the effective length of   ',/,
     .        'C +                 the lateral sponge',/,
     .        'C +   ',/,
     .        '      real*4        vaxgLB              ,vaxdLB',/,
     .        '      real*4        v1xgLB              ,v1xdLB',/,
     .        '     .             ,v2xgLB              ,v2xdLB',/,
     .        '      real          tixgLB              ,tixdLB',/,
     .        '      common/nuddax/vaxgLB(1:n7,my,mz,5),',/,
     .        '     .              vaxdLB(mx-n6:mx ,my,mz,5),',/,
     .        '     .              v1xgLB(1:n7,my,mz,5),',/,
     .        '     .              v1xdLB(mx-n6:mx ,my,mz,5),',/,
     .        '     .              v2xgLB(1:n7,my,mz,5),',/,
     .        '     .              v2xdLB(mx-n6:mx ,my,mz,5),',/,
     .        '     .              tixgLB(2:n7,my,mz  ),',/,
     .        '     .              tixdLB(mx-n6:mx1,my,mz  )',/,
     .        'C +  ',/,
     .        '      real*4        vayiLB              ,vaysLB',/,
     .        '      real*4        v1yiLB              ,v1ysLB',/,
     .        '     .             ,v2yiLB              ,v2ysLB',/,
     .        '      real          tiyiLB              ,tiysLB',/,
     .        '      common/nudday/vayiLB(mx,1:n7,mz,5),',/,
     .        '     .              vaysLB(mx,my-n6:my ,mz,5),',/,
     .        '     .              v1yiLB(mx,1:n7,mz,5),',/,
     .        '     .              v1ysLB(mx,my-n6:my ,mz,5),',/,
     .        '     .              v2yiLB(mx,1:n7,mz,5),',/,
     .        '     .              v2ysLB(mx,my-n6:my ,mz,5),',/,
     .        '     .              tiyiLB(mx,2:n7,mz  ),',/,
     .        '     .              tiysLB(mx,my-n6:my1,mz  )',/,
     .        'C +... vaXX : large scale values of relevant ',/,
     .        'C +           dependant variables ',/,
     .        'C +      ^X=(x->x axis border, y->y axis border)',/,
     .        'C +       ^X=(g->x small, d->x large, ',/,
     .        'C +           b->y small, h->y large) ',/,
     .        'C +    tiXXLB : independant term of semi-implicit',/,
     .        'C +             numerical scheme',/,
     .        'C +    ',/,
     .        '      real          wixgLB',/,
     .        '     .             ,wixdLB',/,
     .        '     .             ,wiyiLB',/,
     .        '     .             ,wiysLB',/,
     .        '      common/nuddaw/wixgLB(   2:  n7,   2:  n7)',/,
     .        '     .             ,wixdLB(mx-n6:mx1,mx-n6:mx1)',/,
     .        '     .             ,wiyiLB(   2:  n7,   2:  n7)',/,
     .        '     .             ,wiysLB(my-n6:my1,my-n6:my1)',/,
     .        'C +...              wiXXLB : coefficient used in',/,
     .        'C +                 semi-implicit numerical scheme',/,
     .        'C +   ',/,
     .        '      real          rxLB    ,ryLB',/,
     .        '      common/nuddtk/rxLB(mx),ryLB(my)',/,
     .        'C +...              rXLB   : nudging coefficients',/,
     .        'C +                          of the relaxation zone',/,
     .        'C +   ',/,
     .        '      real*4        sst_LB',/,
     .        '      real*4        sst1LB,sst2LB',/,
     .        '      common/srfbnd/sst_LB(mx,my),',/,
     .        '     .              sst1LB(mx,my),sst2LB(mx,my)',/,
     .        'C +...              sst_LB : external SST' )

       close (1)
       
       write(6,*) 'MAR include file      MAR_LB.inc_old created'

      ENDIF 


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for surface characteristics : MARdom.dat
C +   ===============================================

C +---Time zone
      ifh = NINT(GElon0/15.)
      if (ifh.gt. 12) ifh=ifh-24
      if (ifh.lt.-12) ifh=ifh+24

C +---i Indices (Surface Output)
      igrdIO(1)=  mx/4
      igrdIO(2)=  mx/4
      igrdIO(3)=  mx/2
      igrdIO(4)=3*mx/4
      igrdIO(5)=3*mx/4
      IF (NSTmod.eq.'M2D') THEN
       igrdIO(1)=1*mx/6
       igrdIO(2)=2*mx/6
       igrdIO(3)=3*mx/6
       igrdIO(4)=4*mx/6
       igrdIO(5)=5*mx/6
      ENDIF

C +---j Indices (Surface Output)
      jgrdIO(1)=  my/4
      jgrdIO(2)=3*my/4
      jgrdIO(3)=  my/2
      jgrdIO(4)=  my/4
      jgrdIO(5)=3*my/4

C +---i/j Indices (Surface Output/Afr West)
      IF (abs(GElat0).lt.15.d0.and.
     .    abs(GElon0).lt.15.d0.and.
     .    NSTmod.ne.'M2D'     .and.
     .    NSTmod.ne.'CPL'          ) THEN
       DO i=1,5
        igrdIO(i) =           4   *mx/5
        jgrdIO(i) =   my/2 + (i-1)*my/10
       ENDDO
      ENDIF

      IF (mmx.eq.1) THEN
       DO i=1,5
        igrdIO(i) = 1
       ENDDO
      ENDIF

      IF (mmy.eq.1) THEN
       DO i=1,5
        jgrdIO(i) = 1
       ENDDO
      ENDIF
      
C +---Vertical adjustment time step
      tequil = 0.
      dtquil = dt

      IF (NSTini.and.SELECT.le.2) THEN

       open  (1,status='unknown',file=NSTdir(1:nbchar)//'MARdom.dat')
       rewind 1

       write (1,141)         LABLio,mmx,mmy,mz
141    format(a3,      15x,
     .          'ON (mx,my,mz) = (',i4,' x',i4,' x',i3,' ) ',
     .                    ' Label + GRID of Simulation')

       write (1,1420)        GElat0,GElon0,GEddxx
1420   format(3d13.6,  14x,' Phi,Lam / x-Axis Direction')
       write (1,1425)        mmaDYN,jdaDYN,jhuDYN,ifh,iyrDYN
1425   format(4i4,i4,  33x,' Month:Day:Hour / Time Zone')

       write (1,1426)        imez,jjmez,maptyp
1426   format(3i4,     41x,' x,y Origin/Projection Type')
       write (1,1427)        igrdIO
1427   format(5i4,     33x,' i Indices (Surface Output)')
       write (1,1428)        jgrdIO
1428   format(5i4,     33x,' j Indices (Surface Output)')
       write (1,1429)        2
1429   format( i4,     49x,' Print Amount Parameter')
       write (1,1421)        1,mmx,1
1421   format(3i4,     41x,' i Output Indices  (MARwri)')
       write (1,1422)        1,mmy,1
1422   format(3i4,     41x,' j Output Indices')
       write (1,1423)        1,mz,1
1423   format(3i4,     41x,' k Output Indices')
       write (1,1424)        mz,min(21,mz)
1424   format(2i4,     45x,' Output Parameters / NetCDF')

       write (1,145)         dx,dy,dt
145    format(3d13.6,  14x,' Hor. Grid Dist./ Time Step')
       write (1,1450)        vertic
1450   format(l3,      50x,' Vertical Grid Type Paramet')
       write (1,1455)        ptopDY
1455   format( d13.6,  40x,' Model   Top Pressure (kPa)')
       write (1,1451)        zmin,aavu,bbvu,ccvu
1451   format(4d13.6,     '  Lowest k + 3Vert.Grid Par.')
       write (1,1452)        FIslot,FIslou,FIslop,FIkhmn
1452   format(4d13.6,     '  Filter Selectivity T, u, p')
       write (1,1453)        TUkhff,TUkhmx
1453   format(2d13.6,  27x,' Horiz.vKar**2  / Up.Sponge')
       write (1,1454)        tequil,dtquil
1454   format(2d13.6,  27x,' 1-D Initialis. Time + Step')
       write (1,1456)        zs_SL,zn_SL,zl_SL,cs2SL
1456   format(4d13.6,     '  z0 Par.Sea/Snow/Land-unused')
       write (1,1457)        sst_SL
1457   format( d13.6,  40x,' SST:(for vert. grid only).')
       write (1,1458)        dtagSL
1458   format( d13.6,  40x,' Initial  T(Air)-T(Surface)')
       write (1,1459)        wk0SL,wx0SL,w20SL,wg0SL
1459   format(4d13.6,     '  Initial Soil Humid.Variab.')

       write (1,1430)
1430   format(' SOIL TYPES')
       IF (NSTmod.eq.'M2D') THEN
        write (1,143)  isol1D
       ELSE
        write (1,143)  ((isolSL(i,j),i=ii1,ii2),j=jj1,jj2)
       ENDIF
143    format((10i13))

       write (1,1431)
1431   format(' TOPOGRAPHY')
       IF (NSTmod.eq.'M2D') THEN
        write (1,1432) SH_1D
       ELSE
        write (1,1432) ((NST_sh(i,j),i=ii1,ii2),j=jj1,jj2)
       ENDIF
1432   format((10d13.6))

       write (1,1433)
1433   format(' ROUGHNESS LENGTH (MOMENTUM)')
       IF (NSTmod.eq.'M2D') THEN
        write (1,1432) z0__1D
       ELSE
        write (1,1432) (((SL_z0(i,j,k),i=ii1,ii2),j=jj1,jj2),k=1,mw)
       ENDIF

       write (1,1434)
1434   format(' ROUGHNESS LENGTH (HEAT,HUMIDITY)')
       IF (NSTmod.eq.'M2D') THEN
        write (1,1432) r0__1D
       ELSE
        write (1,1432) (((SL_r0(i,j,k),i=ii1,ii2),j=jj1,jj2),k=1,mw)
       ENDIF

       write (1,1435)
1435   format(' BULK COEFFICIENT      (HUMIDITY)')
       IF (NSTmod.eq.'M2D') THEN
        write (1,1432) ch0_1D
       ELSE
        write (1,1432) ((ch0SL(i,j),i=ii1,ii2),j=jj1,jj2)
       ENDIF

       write (1,1436)
1436   format(' LEAF SURFACE RESISTANCE')
       IF (NSTmod.eq.'M2D') THEN
        write (1,1432) rsur1D
       ELSE
        write (1,1432) ((rsurSL(i,j),i=ii1,ii2),j=jj1,jj2)
       ENDIF

       write (1,1437)
1437   format(' SURFACE ALBEDO')
       IF (NSTmod.eq.'M2D') THEN
        write (1,1432) alb01D
       ELSE
        write (1,1432) ((alb0SL(i,j),i=ii1,ii2),j=jj1,jj2)
       ENDIF

       write (1,1438)
1438   format(' SURFACE EMISSIVITY (IR)')
       IF (NSTmod.eq.'M2D') THEN
        write (1,1432) eps01D
       ELSE
        write (1,1432) ((eps0SL(i,j),i=ii1,ii2),j=jj1,jj2)
       ENDIF

       write (1,1440)
1440   format(' Rhos Cs sqrt(kappas Tau1) (GROUND)')
       IF (NSTmod.eq.'M2D') THEN
        write (1,1432) d1__1D
       ELSE
        write (1,1432) ((d1_SL(i,j),i=ii1,ii2),j=jj1,jj2)
       ENDIF

c #IT  write (1,1443)
1443   format(' INITIAL GROUND TEMPERATURE   ')
c #IT  write (1,1432) (((tsrfSL(i,j,k),i=ii1,ii2),j=jj1,jj2),k=1,mw)

c #IT  write (1,1444)
1444   format(' INITIAL DEEP   TEMPERATURE   ')
c #IT  write (1,1432) ((t2_SL(i,j),i=ii1,ii2),j=jj1,jj2)

c #po  write (1,1447)
1447   format(' OCEANIC CURRENT (x-Direction)')
c #po  write (1,1432) ((uocnPO(i,j),i=ii1,ii2),j=jj1,jj2)

c #po  write (1,1448)
1448   format(' OCEANIC CURRENT (y-Direction)')
c #po  write (1,1432) ((vocnPO(i,j),i=ii1,ii2),j=jj1,jj2)

c #po  write (1,1449)
1449   format(' LEAD CONCENTRATION')
c #po  write (1,1432) ((aPOlyn(i,j),i=ii1,ii2),j=jj1,jj2)

       write (1,1650)
1650   format(' LONGITUDE')
       write (1,1432) ((NST__x(i,j),i=ii1,ii2),j=jj1,jj2)

       write (1,1651)
1651   format(' LATITUDE')
       write (1,1432) ((NST__y(i,j),i=ii1,ii2),j=jj1,jj2)

       write (1,1652)
1652   format(' SIGMA')
       write (1,1432) sigma

       close (1)

       write(6,*) 'Surface charact. file MARdom.dat created'

      ENDIF 


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for sounding profiles : MARsnd.dat
C +   =========================================

      IF (NSTini.and.SNDing.and.SELECT.eq.1) THEN

       OPEN (unit=2,status='unknown',
     .       file=NSTdir(1:nbchar)//'MARsnd.dat')
       REWIND 2

       zesnd = 0.
       losnd = 0

       WRITE (2,146)         SNDyar,SNDmma,SNDjda,SNDjhu
146    FORMAT(4i4,     37x,' Month:Day:Hour of Sounding')

       WRITE (2,1415)        GElat0,GElon0,GEddxx
1415   FORMAT(3d13.6,  14x,' Phi,Lam / x-Axis Direction')
       WRITE (2,1460)        imez,jjmez
1460   FORMAT(2i4,     45x,' Sounding Coordinate in MAR')

       WRITE (2,1461)
1461   FORMAT(' SOUNDING',
     .      /,'  Temperature Spec. Humid.     Altitude     Pressure',
     .        '  Large Scale Wind (ff,dd)')

       WRITE (2,1462)(SND_tt(i),SND_qv(i),SND_zz(i),SND_pp(i)*10.,
     .                SND_ff(i),TMP_dd(i), i=ms,0,-1)
1462   FORMAT((6d13.6))
 
       WRITE (2,1463) zesnd
1463   FORMAT(d13.6,   40x,' Large Scale Vorticity (2D)')

       WRITE (2,1464) losnd
1464   FORMAT( i4,     49x,' 1  =>  Next Sounding Exist')

       CLOSE (unit=2)

       write(6,*) 'Sounding file         MARsnd.dat created'

      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output file for dynamics : MARdyn.DAT
C     =====================================

      IF (NSTini.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'             ) THEN

       IF (ASCfor) THEN

        open (unit=11,status='unknown',
     .                file=NSTdir(1:nbchar)//'MARdyn.DAT')
        rewind     11
        write     (11,*) itexpe,jdh_LB
        write     (11,*) iyrDYN,mmaDYN,jdaDYN,jhuDYN
        write     (11,*) imez,jjmez
        write     (11,*) GElat0,GElon0
        write     (11,*) sigma,ptopDY,dx,dy
        write     (11,*) NST__u
        write     (11,*) NST__v
        write     (11,*) NSTtmp
        write     (11,*) pstDY
        write     (11,*) NST_qv
        write     (11,*) NST_sh
        write     (11,*) pstDY
        write     (11,*) iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
        write     (11,*) vaxgLB,vaxdLB,vayiLB,vaysLB
        write     (11,*) NST_st
        write     (11,*) uairUB,vairUB,pktaUB ! version MAR > 20/02/04
        close(unit=11)

       ELSE

        open (unit=11,status='unknown',form='unformatted',
     .                file=NSTdir(1:nbchar)//'MARdyn.DAT')
        rewind     11
        write     (11) itexpe,jdh_LB
        write     (11) iyrDYN,mmaDYN,jdaDYN,jhuDYN
        write     (11) imez,jjmez
        write     (11) GElat0,GElon0
        write     (11) sigma ,ptopDY,dx,dy
        write     (11) NST__u
        write     (11) NST__v
        write     (11) NSTtmp
        write     (11) pstDY 
        write     (11) NST_qv
        write     (11) NST_sh
        write     (11) pstDY 
        write     (11) iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
        write     (11) vaxgLB,vaxdLB,vayiLB,vaysLB
        write     (11) NST_st
        write     (11) uairUB,vairUB,pktaUB ! version MAR > 20/02/04
        close(unit=11)

       ENDIF

       write(6,*) 'Initialization   file MARdyn.DAT created'

      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for soil and surface layer : MARsol.DAT
C +   ==============================================

      IF (NSTini.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN

       IF (ASCfor) THEN

        open (unit=11,status='unknown',
     .                file=NSTdir(1:nbchar)//'MARsol.DAT')
        rewind     11
        write     (11,*) itexpe
        write     (11,*) iyrDYN,mmaDYN,jdaDYN,jhuDYN
        write     (11,*) nSLsrf
        write     (11,*) SLsrfl
        write     (11,*) tairSL
        write     (11,*) tsrfSL
        write     (11,*) alb0SL,eps0SL
        write     (11,*) SaltSL
        write     (11,*) ro_SL 
        write     (11,*) ro_SL
        write     (11,*) d1_SL
        write     (11,*) t2_SL
        write     (11,*) w2_SL,wg_SL
        write     (11,*) roseSL
        write     (11,*) qvapSL
        write     (11,*) hsnoSL
        write     (11,*) hmelSL
        write     (11,*) SLuusl,SL_z0
        write     (11,*) SLutsl,SL_r0
        close(unit=11)

       ELSE

        open (unit=11,status='unknown',form='unformatted',
     .                file=NSTdir(1:nbchar)//'MARsol.DAT')
        rewind     11
        write     (11) itexpe
        write     (11) iyrDYN,mmaDYN,jdaDYN,jhuDYN
        write     (11) nSLsrf
        write     (11) SLsrfl
        write     (11) tairSL
        write     (11) tsrfSL
        write     (11) alb0SL,eps0SL
        write     (11) SaltSL
        write     (11) ro_SL  
        write     (11) ro_SL 
        write     (11) d1_SL 
        write     (11) t2_SL 
        write     (11) w2_SL ,wg_SL
        write     (11) roseSL
        write     (11) qvapSL
        write     (11) hsnoSL
        write     (11) hmelSL
        write     (11) SLuusl,SL_z0
        write     (11) SLutsl,SL_r0
        close(unit=11)

       ENDIF

       write(6,*) 'Initialization   file MARsol.DAT created'

      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for SVAT model : MARsvt.DAT
C +   ==================================

      IF (NSTini.and.SVTmod.and.
     .    LoutDA.and.SELECT.eq.1) THEN

       IF (ASCfor) THEN

        open (unit=11,status='unknown',
     .        file=NSTdir(1:nbchar)//'MARsvt.DAT')
        rewind     11
        write     (11,*) itexpe
        write     (11,*) iyrDYN,mmaDYN,jdaDYN,jhuDYN
        write     (11,*) igrdIO
        write     (11,*) jgrdIO
        IF (NSTmod.eq.'M2D') THEN
         write    (11,*) veg_1D
         write    (11,*) iwf_1D
         write    (11,*) dsa_1D
         write    (11,*) svt_1D
         write    (11,*) sfr_1D
         write    (11,*) lai_1D
         write    (11,*) glf_1D
         write    (11,*) ts__1D
         write    (11,*) sw__1D
        ELSE
         write    (11,*) ((isolTV(i,j)   ,i=ii1,ii2),j=jj1,jj2)
         write    (11,*) ((NSTiwf(i,j)   ,i=ii1,ii2),j=jj1,jj2)
         write    (11,*) ((NSTdsa(i,j)   ,i=ii1,ii2),j=jj1,jj2)
         write    (11,*) (((NSTsvt(i,j,k),i=ii1,ii2),j=jj1,jj2),
     .                                    k=1,nvx)
         write    (11,*) (((NSTsfr(i,j,k),i=ii1,ii2),j=jj1,jj2),
     .                                    k=1,nvx)
         write    (11,*) (((NSTlai(i,j,k),i=ii1,ii2),j=jj1,jj2),
     .                                    k=1,nvx)
         write    (11,*) (((NSTglf(i,j,k),i=ii1,ii2),j=jj1,jj2),
     .                                    k=1,nvx)
         write    (11,*) ((((NST_ts(i,j,k,l),i=ii1,ii2),j=jj1,jj2),
     .                                       k=1,nvx),l=1,nsl)
         write    (11,*) ((((NST_sw(i,j,k,l),i=ii1,ii2),j=jj1,jj2),
     .                                       k=1,nvx),l=1,nsl)
        ENDIF
        close(unit=11)

       ELSE

        open (unit=11,status='unknown',form='unformatted',
     .        file=NSTdir(1:nbchar)//'MARsvt.DAT')
        rewind     11
        write     (11) itexpe
        write     (11) iyrDYN,mmaDYN,jdaDYN,jhuDYN
        write     (11) igrdIO
        write     (11) jgrdIO
        IF (NSTmod.eq.'M2D') THEN
         write    (11) veg_1D
         write    (11) iwf_1D
         write    (11) dsa_1D
         write    (11) svt_1D
         write    (11) sfr_1D
         write    (11) lai_1D
         write    (11) glf_1D
         write    (11) ts__1D
         write    (11) sw__1D
        ELSE
         write    (11) ((isolTV(i,j)      ,i=ii1,ii2),j=jj1,jj2)
         write    (11) ((NSTiwf(i,j)      ,i=ii1,ii2),j=jj1,jj2)
         write    (11) ((NSTdsa(i,j)      ,i=ii1,ii2),j=jj1,jj2)
         write    (11) (((NSTsvt(i,j,k)   ,i=ii1,ii2),j=jj1,jj2),
     .                                     k=1,nvx)
         write    (11) (((NSTsfr(i,j,k)   ,i=ii1,ii2),j=jj1,jj2),
     .                                     k=1,nvx)
         write    (11) (((NSTlai(i,j,k)   ,i=ii1,ii2),j=jj1,jj2),
     .                                     k=1,nvx)
         write    (11) (((NSTglf(i,j,k)   ,i=ii1,ii2),j=jj1,jj2),
     .                                     k=1,nvx)
         write    (11) ((((NST_ts(i,j,k,l),i=ii1,ii2),j=jj1,jj2),
     .                                     k=1,nvx)  ,l=1,nsl)
         write    (11) ((((NST_sw(i,j,k,l),i=ii1,ii2),j=jj1,jj2),
     .                                     k=1,nvx)  ,l=1,nsl)
        ENDIF
        close(unit=11)

       ENDIF

       write(6,*) 'Initialization   file MARsvt.DAT created'

      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for boundary forcing : MARglf.DAT
C +   ========================================

      IF (NSTini.and.SVTmod.and.
     .    LoutDA.and.SELECT.eq.1) THEN

       NSTinfo = 'NESTOR_3.3'

       IF (ASCfor) THEN
        open (unit=13,status='unknown',
     .                file=NSTdir(1:nbchar)//'MARglf.DAT')
        rewind     13
       ELSE
        open (unit=13,status='unknown',form='unformatted',
     .                file=NSTdir(1:nbchar)//'MARglf.DAT')
        rewind     13
       ENDIF

       write(6,*) 'SVAT evolutive   file MARglf.DAT created'

      ENDIF

      IF (NSTfor.and.SVTmod.and.
     .    LoutDA.and.SELECT.eq.1) THEN
       IF (ASCfor) THEN
        write     (13,*) iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
c #NI.                  ,NSTinfo
        IF (NSTmod.eq.'M2D') THEN
         write    (13,*) glf_1D 
         write    (13,*) lai_1D 
        ELSE
         write    (13,*) (((NSTglf(i,j,k),i=ii1,ii2),j=jj1,jj2),
     .                                               k=1,nvx)
         write    (13,*) (((NSTlai(i,j,k),i=ii1,ii2),j=jj1,jj2),
     .                                               k=1,nvx) 
        ENDIF
       ELSE
        write     (13) iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
c #NI.                ,NSTinfo
        IF (NSTmod.eq.'M2D') THEN
         write    (13)   glf_1D           
         write    (13)   lai_1D           
        ELSE
         write    (13)   (((NSTglf(i,j,k),i=ii1,ii2),j=jj1,jj2),
     .                                               k=1,nvx)
         write    (13)   (((NSTlai(i,j,k),i=ii1,ii2),j=jj1,jj2),
     .                                               k=1,nvx) 

        ENDIF
       ENDIF
       write(6,*) 'SVAT evolutive   file MARglf.DAT appended'
      ENDIF

      IF (NSTend.and.SVTmod.and.
     .    LoutDA.and.SELECT.eq.1) THEN
       CLOSE(unit=13)      
      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for boundary forcing : MARlbc.DAT
C +   ========================================

      IF (NSTini.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN

       IF (ASCfor) THEN
        open (unit=12,status='unknown',
     .                file=NSTdir(1:nbchar)//'MARlbc.DAT')
        rewind     12
       ELSE
        open (unit=12,status='unknown',form='unformatted',
     .                file=NSTdir(1:nbchar)//'MARlbc.DAT')
        rewind     12
       ENDIF

       write(6,*) 'Boundary forcing file MARlbc.DAT created'

      ENDIF

      IF (NSTfor.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN
       IF (ASCfor) THEN
        write     (12,*) iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
        write     (12,*) vaxgLB,vaxdLB,vayiLB,vaysLB
        write     (12,*) NST_st 
       ELSE
        write     (12)   iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
        write     (12)   vaxgLB,vaxdLB,vayiLB,vaysLB
        write     (12)   NST_st 
       ENDIF
       write(6,*) 'Boundary forcing file MARlbc.DAT appended'
      ENDIF

      IF (NSTend.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN
       CLOSE(unit=12)      
      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for boundary forcing : MARubc.DAT (version MAR > 20/02/04)
C +   ========================================

      IF (NSTini.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN

       IF (ASCfor) THEN
        open (unit=17,status='unknown',
     .                file=NSTdir(1:nbchar)//'MARubc.DAT')
        rewind     17
       ELSE
        open (unit=17,status='unknown',form='unformatted',
     .                file=NSTdir(1:nbchar)//'MARubc.DAT')
        rewind     17
       ENDIF

       write(6,*) 'Boundary forcing file MARubc.DAT created'

      ENDIF

      IF (NSTfor.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN
       IF (ASCfor) THEN
        write     (17,*) iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
        write     (17,*) uairUB,vairUB,pktaUB
       ELSE
        write     (17) iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
        write     (17) uairUB,vairUB,pktaUB
       ENDIF
       write(6,*) 'Boundary forcing file MARubc.DAT appended'
      ENDIF

      IF (NSTend.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN
       CLOSE(unit=17)      
      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for boundary forcing : MARsic.DAT (version MAR > 20/02/04)
C +   ========================================

      IF (NSTini                  .and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN

       IF (ASCfor) THEN
        open (unit=16,status='unknown',
     .                file=NSTdir(1:nbchar)//'MARsic.DAT')
        rewind     16
       ELSE
        open (unit=16,status='unknown',form='unformatted',
     .                file=NSTdir(1:nbchar)//'MARsic.DAT')
        rewind     16
       ENDIF

       write(6,*) 'Boundary forcing file MARsic.DAT created'

      ENDIF

      IF (NSTfor.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN
       IF (ASCfor) THEN
        write     (16,*) iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
        write     (16,*) NSTsic 
       ELSE
        write     (16)   iyrDYN,mmaDYN,jdaDYN,jhuDYN,jdh_LB
        write     (16)   NSTsic 
       ENDIF
       write(6,*) 'Boundary forcing file MARsic.DAT appended'
      ENDIF

      IF (NSTend.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.ne.'M2D'         .and.
     .    NSTmod.ne.'CPL'              ) THEN
       CLOSE(unit=16)      
      ENDIF


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Output for sea surface temperature : MARsst.dat
C +   ===============================================

      IF (NSTini.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.eq.'M2D'              ) THEN

       open (unit=14,status='unknown',
     .       file=NSTdir(1:nbchar)//'MARsst.dat')
       rewind 14
       write(6,*) 'Sea surface temp file MARsst.dat created'

      ENDIF

      IF (NSTfor.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.eq.'M2D'              ) THEN
       DO i=1,mx
        write(14,*) sst1D(i)
       ENDDO
       write(6,*) 'Sea surface temp file MARsst.dat appended'
      ENDIF

      IF (NSTend.and.(.not.SNDing).and.
     .    LoutDA.and.SELECT.eq.1  .and.
     .    NSTmod.eq.'M2D'              ) THEN
       CLOSE(unit=14)      
      ENDIF

      write(6,*)


C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
