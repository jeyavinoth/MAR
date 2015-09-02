C   +-------------------------------------------------------------------+
C   |  Subroutine MARvgd                                     01-09-2004 |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Creation of the vertical grid of the MAR model.                   |
C   |                                                                   |
C   | Input : - VGD_sp : surface pressure (kPa)                         |
C   | ^^^^^^^ - fID    : identificator of the Netcdf data file          |
C   |         - nz     : number of vertical levels                      |
C   |         - klev   : if specified, the level at which pressure and  |
C   |                    hybrid coordinate has to be computed           |
C   |         - parameters from MARgrd.ctr                              |
C   |                                                                   |
C   | Output: Creation of the vertical MAR grid given in hybrid         |
C   | ^^^^^^^ coordinates :                                             |
C   |         - VGD_hp(nz+1) : local hybrid coord. for vertic. interp.  |
C   |         - VGD__p(nz+1) : pressure coordinates (kPa)               |
C   |         - VGDgdz(nz  ) : model coordinates (sigma)                |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE MARvgd (LSCmar,klev,nz,fileID,VGD_sp,VGD_hp,
     .                   VGDgdz,VGD__p ,sigma ,WK2_1D,WK3_1D)


      IMPLICIT NONE

      INCLUDE 'CTRvar.inc'
      INCLUDE 'NSTdim.inc'

C +---Local variables
C +   ---------------

      INTEGER k,klev,nz,k1,k2,fileID,maptyp,imez,jmez

      REAL pp1,pps,ppm,dpsl,pp,hh,ppf,GElat0,GElon0,dx,GEddxx,
     .     ptopDY,zmin,aavu,bbvu,ccvu,sst_SL,TUkhmx,empty1(1)
      SAVE ptopDY
     
      REAL VGD_sp,VGD_hp(nz+1),VGDgdz(nz),VGD__p(nz+1),sigma(nz),
     .     WK2_1D(nz),WK3_1D(nz)
     

      LOGICAL LSCmar,vertic

      CHARACTER*10 var_units

      REAL           MARsig(mz)
      COMMON/cMARvgd/MARsig  
C           See MARout.f

      LOGICAL lfirst
      SAVE    lfirst
      DATA lfirst/.true./
      
      IF (mz.NE.nz) THEN 
         write(*,*) 'Wrong #levels in MARvgd ?'
         STOP
      ENDIF
      
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +                          Begin (First call)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CREATION OF SIGMA MAR GRID USING LSC NetCDF FILE
C +   ================================================

      IF (LSCmar.and.lfirst_LSC) THEN


C +---Read SIGMA in NetCDF file
C +   -------------------------

C +         ******
       CALL UNsread (fileID,'level',0,0,
     .               0,0,nz,1,1,var_units,sigma)
C +         ******


C +---Sigma coordinates
C +   -----------------

       DO k=1,nz
        VGDgdz(k)=sigma(k)
        MARsig(k)=sigma(k)
       ENDDO

C +    write(6,*) 'Sigma grid : '
C +    write(6,*) sigma
C +    write(6,*) 'STOP for verification'
C +    STOP
 
 
C +---End first call
C +   --------------

       lfirst_LSC = .false.

      ENDIF


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CREATION OF SIGMA MAR GRID USING PARAMETERS IN MARgrd.ctr
C +   =========================================================

      IF ((.not.LSCmar).and.lfirst) THEN


C +---Read grid parameters in MARgrd.ctr
C +   ----------------------------------

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
        read (51,*) NSTfis
        read (51,*) !- - - - - - - - - - - - - - - - - -


       CLOSE(unit=51)

C +---Sets the standard values of vertical grid parameters
C +   ----------------------------------------------------

C +         ******
       CALL SETsig (nz,zmin,aavu,bbvu,ccvu,ptopDY)
C +         ******


C +---Computation of vertical grid
C +   ----------------------------

C +         ******
       CALL GRDsig(nz,zmin,aavu,bbvu,ccvu,vertic,
     .                sst_SL,TUkhmx,sigma,WK2_1D)
C +         ******

C +---Sigma coordinates
C +   -----------------

       WK3_1D(1)=ptopDY   ! Store ptopDY

       DO k=1,nz
        VGDgdz(k)=sigma(k)
        MARsig(k)=sigma(k)
       ENDDO

C +---End first call
C +   --------------

       lfirst     = .false.

      ELSE
      
        DO k=1,nz
          VGDgdz(k)=MARsig(k)
        ENDDO

      ENDIF
      lfirst_NST = .false.


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +   End first call section
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---HYBRID AND PRESSURE COORDINATES (required by the nesting code)
C +   ===============================


C +---Reference levels for hybrid coordinates
C +   ---------------------------------------

      pp1  = 105.       ! Reference pressure (KPa)
      dpsl = 20.        ! "> boundary layer" (KPa)


C +---Selection of vertical levels
C +   ----------------------------

      IF ((klev.le.0).or.(klev.gt.nz)) THEN
       k1=1
       k2=nz
      ELSE
       k1=1
       k2=klev
      ENDIF


C +---Computation of hybrid coordinates used in vertic. interp.
C +   ---------------------------------------------------------
      
      IF (ptopDY.LT.1.0E-10) THEN
         write(*,*) 'Something is probably going wrong in MARvgd'
         write(*,*) 'ptopDY= ',ptopDY
         write(*,*)
      ENDIF
      
      pps = VGD_sp
      ppm = pps - dpsl
      DO k = k1,k2+1
       IF (k.eq.(nz+1)) THEN
        pp = VGD_sp
       ELSE
        pp = VGDgdz(k)*(VGD_sp-ptopDY) + ptopDY
       ENDIF
       hh = pp/pp1
       IF (pp.gt.ppm) THEN
        ppf= (pp-ppm)/(pps-ppm)
        hh = hh + (pp1-pps)/pp1 * ppf * ppf
       END IF
       VGD_hp(k) = LOG(hh)
       VGD__p(k) = pp
      ENDDO

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      subroutine SETsig (mz,zmin,aavu,bbvu,ccvu,pt)

C +------------------------------------------------------------------------+
C |   SubRoutine SETsig sets the standard values of vert grid parameters   |
C +------------------------------------------------------------------------+
C |   INPUT  : mz, + all other arguments if not equal to '0'               |
C |   ^^^^^^^^                                                             |
C |                                                                        |
C |                                                                        |
C |   OUTPUT : zmin           : Height above Surface / 1st sigma level (m) |
C |   ^^^^^^^^ aavu,bbvu,ccvu : Vertical Discretization Parameters         |
C |            pt             : Pressure at the top of the model (kPa)     |
C |                                                                        |
C +------------------------------------------------------------------------+

      IMPLICIT NONE
     
      INTEGER mz,mmz

      REAL    aavu,bbvu,ccvu,zmin,pt,INzmin,INaavu,INbbvu,INccvu

      mmz    = mz
      INzmin = zmin
      INaavu = aavu
      INbbvu = bbvu
      INccvu = ccvu

C     Prescribed values for vertical grid parameters
C     ----------------------------------------------

      IF (mmz.eq.10) THEN
       zmin = 10.
       aavu = 2.0
       bbvu = 1.25
       ccvu = 2000.
      ENDIF

      IF (mmz.eq.18.or.mmz.eq.19) THEN
       zmin = 10.
       aavu = 2.0
       bbvu = 1.17
       ccvu = 1200.
      ENDIF

      IF (mmz.eq.20.or.mmz.eq.29) THEN
       zmin = 3.
       aavu = 1.8
       bbvu = 1.13
       ccvu = 1000.
      ENDIF

      IF (mmz.eq.30.or.mmz.eq.31) THEN
       zmin = 5.
       aavu = 2.0
       bbvu = 1.11
       ccvu = 900.
      ENDIF

      IF (mmz.eq.40) THEN
       zmin = 5.
       aavu = 2.0
       bbvu = 1.06
       ccvu = 2500.
      ENDIF

      IF (mmz.eq.60) THEN
       zmin = 2.
       aavu = 2.0
       bbvu = 1.10
       ccvu = 70.
      ENDIF

C     Forcing with values given in MARgrd.ctr
C     ---------------------------------------

      IF (INzmin.ne.0.) zmin=INzmin
      IF (INaavu.ne.0.) aavu=INaavu
      IF (INbbvu.ne.0.) bbvu=INbbvu
      IF (INccvu.ne.0.) ccvu=INccvu

C     If insufficient informations ...
C     --------------------------------

      IF (zmin.eq.0..or.aavu.eq.0..or.bbvu.eq.0..or.ccvu.eq.0.) THEN
       WRITE(6,*) 'Chooses other parameters for z-grid Set-Up!'
       WRITE(6,*) 'Program is stopped in SETsig!'
       STOP
      ENDIF
       
      RETURN
      END


      subroutine GRDsig(mz,zmin,aavu,bbvu,ccvu,vertic,
     .                       sst_SL,TUkhmx,sigma,zpbl)
C +
C +------------------------------------------------------------------------+
C | MAR GRID                                               19-06-2004  MAR |
C |   SubRoutine GRDsig is used to initialize the vertical grid            |
C |                                                                        |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   ASSUMPTION: Sigma is calculated from initial level height amsl       |
C |   ^^^^^^^^^^^                     assumig that T(msl) = SST            |
C |                                                dT/dz  = -0.0065 K/m    |
C |                                                p_s    = 100     hPa    |
C |                                                                        |
C |   INPUT  : zmin           : Height above Surface / 1st Sigma Level (m) |
C |   ^^^^^^^^ aavu,bbvu,ccvu : Vertical Discretization Parameters         |
C |            vertic         : Logical Variable caracteris.vertic.discris.|
C |                                                                        |
C |   OUTPUT : Variable  which is  initialized is:                         |
C |   ^^^^^^^^  sigma(mz): Independant Variable (Normalized Pressure)      |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
      implicit none
C +
C +
C +--General Variables
C +  =================
C +
      integer k,kk,mzz,mz
C +
      real zmin,aavu,bbvu,ccvu,zpbl(mz),sigma(mz),ps_sig,ga,
     .     ga0,aa,bb,cc,vu,ra,gravit,unun,sst_SL,dzz,rz,rzb,
     .     TUkhmx,zzo,zero,epsi
C +
      logical vertic

C +
C +--DATA
C +  ====
C +
      data ps_sig / 101.3d0 /
C +
      data ga0    / 0.0065d0/
C +...     ga0 : Standard Atmospheric Lapse Rate
      data ra     / 287.d0  /
      data gravit / 9.81d0  /
      data unun   / 1.d0    /
      data zero   / 0.d0    /
      data epsi   / 1.0d-6  /
    
C +
C +--Initialization
C +  ==============
C +
      mzz=mz+1 
      DO k=1,mz
       zpbl(k)=0.
      ENDDO
C +
C +--Temperature Vertical Profile
C +  ============================
C +
      ga  = ga0
C +
C +--Sigma Levels
C +  ============
C +
C +- 1) Coarse Resolution of the Surface Layer
C +  -----------------------------------------
C +
      if (.not.vertic) then
C +
C +    aa = 0.5
C +    bb = 1.5
C +    cc =-1.0
C +... Reference : E. Richard, these, 1991, p.29
C +
       vu =       0.0d0
       do k=1,mz
       vu = vu  + 1.0d0/dble(mzz)
       sigma(k) = aavu*vu + bbvu*vu*vu + ccvu*vu*vu*vu
C +
       if (abs(ga).gt.1.d-5) then
        zpbl(k) =-(   sst_SL  /ga) *   ((1.d0+(sigma(k)-1.d0)
     .                                 *(1.d2/ps_sig))
     .                                 **(ra*ga/gravit)-1.d0)
       else
        zpbl(k) =-(ra*sst_SL  /gravit ) *log((unun+(sigma(k)-unun)
     .                                 *(1.d2/ps_sig)))
       end if
       enddo
C +
C +
C +- 2) Fine   Resolution of the Surface Layer
C +  -----------------------------------------
C +
      else
C +
       ga     =max(ga,epsi)
C +
       zpbl(1)=     zmin
       zpbl(2)=2.d0*zmin
C +
       dzz    =0.d0
       
       do k=3,mz
        rz     =zmin*aavu **FLOAT(k-1)
        rzb    =ccvu*bbvu **FLOAT(k-1)
        zpbl(k)=rzb   *rz /(rz + rzb  )
C +
        zzo    = zpbl(k)
        zpbl(k)= max(zpbl(k),zpbl(k-1)+zpbl(2))
        dzz    = max(zpbl(k)-zzo,      zero   ) + dzz
       enddo
       
       do k=1,mz
        kk=mz+1-k
       
C +     sigma(kk)=(ps_sig/100.d0) 
C +               Arbitraire et pas utile ˆ mon avis (PhM)
C +  
        sigma(kk)= 1.0d0
     .         *((1.0d0-ga*zpbl(k)/sst_SL)**(gravit/(ga*ra))-1.d0)+1.d0

C +...  sigma(kk): the fine resolution of the surface layer is computed
C +            using a geometric progression
C +
       enddo
      end if
C +
      do k=1,mz
       if (sigma(k).le.0.0)then

        print *, "ERROR in MARvgd.f: sigma < 0."

        do kk=1,mz
         print *,kk,sigma(mz+1-kk),zpbl(kk)
        enddo

        print *, "Change aavu,bbvu,ccvu in MARgrd.ctr or mz"
        print *, "For example: try to decrease ccvu in MARgrd.ctr"
        stop

       endif
      enddo
C +
      write(*,*) 'MAR vertical grid created'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~ '
      write(*,*) ' '

      return
      end
