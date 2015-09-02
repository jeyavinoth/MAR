C   +-------------------------------------------------------------------+
C   |  Subroutine CORveg                             08/2004    NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Derives MAR vegetation and surface type from CORINE               |
C   |     from CORINE land use database  (revised version 2004)         |
C   |                                                                   |
C   | Input : - NST__x, NST__y : horizontal grid of NST model           |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: - NSTsvt : vegetation type (SVAT classification)          |
C   | ^^^^^^^ - NSTsfr : fraction of vegetation in the grid cell (SVAT) |
C   |         - NSTlai : leaf area index                                |
C   |                                                                   |
C   |                                                                   |
C   | Data source :  CORINE land cover (Corine detail level 2 or 3)     |
C   | ^^^^^^^^^^^                                                       |
C   |                                                                   |
C   |                                                                   |
C   |                                                                   |
C   +-------------------------------------------------------------------+
      SUBROUTINE CORveg


      IMPLICIT NONE

C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'LOCfil.inc'
      INCLUDE 'NetCDF.inc'
      INCLUDE 'NESTOR.inc'

C +---Local variables
C +   ---------------

      INTEGER imx,imy,nclass,nsvat
 
C     CORINE grid size
      PARAMETER (imx=18294,imy=18514)
      PARAMETER (nclass=50)
      PARAMETER (nsvat =12)
      
      INTEGER  VARSIZE
      EXTERNAL VARSIZE

      INTEGER ii,jj,kk,ll,Ierror,lomi,
     .        nbchar,iclass,CORcid,nFrNul,
     .        AXX_ID,AXY_ID,frac_itot,AgFlag
      LOGICAL AgNeed

      REAL    FrMwat,FrMfor,FrMnul,
     .        frac_tot,degrad,reqlon,
     .        reqlat,radius,phi0,lam0,C_reso,xlaea,ylaea
     
      REAL    out_X(0:mx,0:my), out_Y(0:mx,0:my),
     .        ainX(imx),ainY(imy)


      INTEGER SVT_class(nvx)
      
      REAL    CORfrc(mx,my,nclass),VEGlon,
     .        convert(nclass,0:nsvat),Curban(nclass),
     .        CORwat(nclass),Cnoveg(nclass),CORtmp(nclass),
     .        COR_z0(nclass),FrAsSVT(0:nsvat),SVT_frac(nvx)
     
      CHARACTER*200 CHnul 

C +---Data
C +   ----

C     CORINE Map projection information
      DATA degrad  /  1.745329252d-2     /
      DATA lam0    /  9.0                /
      DATA phi0    / 48.0                /
      DATA radius  / 6378388.0           /

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Screen message
C +   ==============

      write(6,*) 'CORINE land cover over Europe'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*)


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



C +---Convertion table : CORINE  -> SVAT classification
C +   =================================================


C +---Initialisation
C +   --------------

      DO kk=1,nclass
      DO ll=0,nsvat
       convert(kk,ll)=0.
      ENDDO
      ENDDO


C +---Convertion table
C +   ----------------
      nbchar=VARSIZE(CORveg_dir)
      OPEN (unit=40,status='old',
     .      file=CORveg_dir(1:nbchar)//'CORINEtab.txt')

      READ(40,'(1x)') 
      READ(40,'(1x)') 
      READ(40,'(1x)') 
      READ(40,'(1x)') 
      DO kk=1,nclass
       READ(40,*) iclass,(convert(kk,ll),ll=0,nsvat),Cnoveg(kk),
     .            Curban(kk),CORwat(kk),COR_z0(kk)
       IF(kk.ne.iclass)THEN 
         write(*,*) 'CORveg: error reading CORINEtab.txt '
         STOP
       ENDIF
      ENDDO
      CLOSE (unit=40)

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Initialization
C +   ==============

C      Currently disabled because CORINE is intended to be used
C      "where data is available" only :
C      The initalisation should be performed in GLOveg (IGBP database),
C      and after this CORveg will change the values where
C      there are CORINE data available.

       IF (.NOT.VEGdat) THEN
          write(*,*) 'CORveg: '
          write(*,*) '======= '
          write(*,*) 'This routine should be used only after GLOveg'
          write(*,*) 'Please activate GLOveg (set VEGdat to true)'
          write(*,*) 'or see comments in CORveg ?'
          STOP
        ENDIF

C --    Please delete this commented code in 2006 if still unused...
C --      DO jj=1,my
C --      DO ii=1,mx
C --
C --       DO kk=1,nvx-1
C --        NSTsvt(ii,jj,kk)=0
C --        NSTsfr(ii,jj,kk)=0
C --       ENDDO
C --       
C --       IF (NSTsol(ii,jj).ge.4) THEN ! Continental areas
C --        NSTsvt(ii,jj,nvx)=  6
C --        NSTsfr(ii,jj,nvx)=100
C --        DO kk=1,nvx
C --         NSTlai(ii,jj,kk) = 2.0
C --         NSTglf(ii,jj,kk) = 1.0
C --        ENDDO
C --        
C --       ELSE
C --        NSTsvt(ii,jj,nvx)=  0
C --        NSTsfr(ii,jj,nvx)=100
C --        DO kk=1,nvx
C --         NSTlai(ii,jj,kk) = 0.0
C --         NSTglf(ii,jj,kk) = 0.0
C --        ENDDO
C --       ENDIF
C --
C --      ENDDO
C --      ENDDO

          nFrNul=0

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Find coordinates of the NST grid meshes boundaries on the 
C +   CORINE MAP (Lambert EA projection)
C +   ==========================================================

      DO jj=1,my-1
      DO ii=1,mx-1

       reqlon = (NST__x(ii,jj) + NST__x(ii+1,jj+1)) / 2.
       reqlat = (NST__y(ii,jj) + NST__y(ii+1,jj+1)) / 2.
C +    Approximation valid for "small" meshes of any orientation

C +         ***********
       CALL lamphi2laea (xlaea,ylaea,reqlon*degrad,reqlat*degrad,
     .                   lam0*degrad,phi0*degrad,radius)
C +         ***********

       out_X (ii,jj) = xlaea
       out_Y (ii,jj) = ylaea

      ENDDO
      ENDDO
C
C +   Domain boundaries
      DO ii=1,mx-1
       out_X (ii,0) = 2. * out_X (ii,1)    - out_X (ii,2)
       out_Y (ii,0) = 2. * out_Y (ii,1)    - out_Y (ii,2)
       out_X (ii,my)= 2. * out_X (ii,my-1) - out_X (ii,my-2)
       out_Y (ii,my)= 2. * out_Y (ii,my-1) - out_Y (ii,my-2)
      ENDDO 
      DO jj=0,my
       out_X (0,jj) = 2. * out_X (1,jj)    - out_X (2,jj)
       out_Y (0,jj) = 2. * out_Y (1,jj)    - out_Y (2,jj)
       out_X (mx,jj)= 2. * out_X (mx-1,jj) - out_X (mx-2,jj)
       out_Y (mx,jj)= 2. * out_Y (mx-1,jj) - out_Y (mx-2,jj)
      ENDDO 
      
C +---Read and "upscale" CORINE data
C +   ==============================

     
C +   Open the CORINE data file 
C +   -------------------------

      Ierror=NF_OPEN(CORveg_dir(1:nbchar)//'CORINE250.nc',
     .              NF_NOWRITE,CORcid)
      Ierror=NF_INQ_VARID(CORcid,'x'   ,AXX_ID)
      Ierror=NF_INQ_VARID(CORcid,'y'   ,AXY_ID)
      Ierror=NF_GET_VAR_REAL(CORcid, AXX_ID, ainX)
      Ierror=NF_GET_VAR_REAL(CORcid, AXY_ID, ainY)
      
C +   Read & upscale
C +   --------------

C +        ******
      CALL UPScor (ainX, ainY, imx, imy, nclass,
     .                 CORcid, 'luse', out_X, out_Y, corFRC) 
C +        ******


C +---Close Netcdf data file
C +   ----------------------
      CALL UNCLOSE(CORcid)



C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO jj=1,my
      DO ii=1,mx
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Specific areas
C     ==============

       FrMwat = 0.

C +---Water areas
C +   -----------
        DO kk=38,nclass
         FrMwat = FrMwat + CORfrc(ii,jj,kk)
        ENDDO

C +---Forest areas
C +   ------------
        FrMfor = (CORfrc(ii,jj,23)+CORfrc(ii,jj,24)+CORfrc(ii,jj,25))

        IF (FrMfor.gt.0.8.and.RUGdat) THEN
         FrMfor      = MIN(1.0,FrMfor)
C +  ?   NST_z0(ii,jj) = NST_z0(ii,jj) * (1.0-0.3*(FrMfor-0.8)/0.2)
C +  ?   NST_r0(ii,jj) = 0.1*NST_z0(ii,jj)
        ENDIF

C +---No-data area in CORINE
C     ----------------------
        FrMnul = CORfrc(ii,jj,49)
        IF (FrMnul .GT. 0.5) THEN
C          WRITE (*,*) 'No-data  area in CORINE (i,j,%) :'
C          WRITE (*,*) ii,jj,CORfrc(ii,jj,49)*100.
           nFrNul=nFrNul+1
        ENDIF


C +    **************************************************************
       IF (FrMnul.lt.0.5) THEN  ! Not undefined in CORINE
       IF (FrMwat.lt.0.5) THEN  ! continent
C +    **************************************************************

       NSTsol(ii,jj)=4 ! unclear  ==============================

C +---Convertion of CORINE to SVAT classification
C +   ===========================================


C +...  initialisation
C +     ~~~~~~~~~~~~~~
        DO ll=0,nsvat
         FrAsSVT(ll)=0.
        ENDDO

C +...  convertion to SVAT classes
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO kk=1,nclass
         DO ll=1,nsvat
          FrAsSVT(ll)=FrAsSVT(ll)+convert(kk,ll)
     .                            /100.*CORfrc(ii,jj,kk)
         ENDDO
        ENDDO
        
C +...  If too much classes / nvx -> aggregate some
        AgNeed=.TRUE. ! First call must test if aggr. needed
        AgFlag= 0     ! No aggregation performed here
        
        CALL AgClasses(FrAsSVT,nsvat, 5, 6, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat, 5, 4, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat, 2, 1, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat, 2, 3, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat,11,10, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat,11,12, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat, 8, 7, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat, 8, 9, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat, 5, 2, AgNeed,AgFlag)
        CALL AgClasses(FrAsSVT,nsvat, 8,11, AgNeed,AgFlag)

C +...  retain the nvx dominant classes
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         DO kk=1,nvx
           lomi=0
           DO ll=0,nsvat
             IF (FrAsSVT(ll).GT.FrAsSVT(lomi)) THEN
               lomi=ll
             ENDIF
           ENDDO
         SVT_class(kk)=lomi
         SVT_frac (kk)=FrAsSVT(lomi)
         FrAsSVT(lomi)=0.0
         ENDDO

C +...  normalizing the three dominant fractions
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        frac_tot=0.
        DO ll=1,nvx
         frac_tot=frac_tot+SVT_frac(ll)
        ENDDO
        IF (frac_tot.ne.0.) THEN
         DO ll=1,nvx
          SVT_frac(ll)=SVT_frac(ll)/frac_tot
         ENDDO
        ENDIF

C +...  attribute classes and fractions to NST variables
C +     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO kk=1,nvx
         NSTsvt(ii,jj,kk)=     SVT_class(kk)
         NSTsfr(ii,jj,kk)=NINT(SVT_frac (kk)*100.)
        ENDDO

C +---Final check of soil fractions
C +   =============================

        frac_itot=0
        DO ll=1,nvx
         frac_itot=frac_itot+NSTsfr(ii,jj,ll)
        ENDDO

        IF (frac_itot.le.0) THEN   ! Imposed bare soil
         NSTsvt(ii,jj,nvx)=  0
         NSTsfr(ii,jj,nvx)=100
         DO kk=1,nvx-1
          NSTsvt(ii,jj,kk)=0
          NSTsfr(ii,jj,kk)=0
         ENDDO
         write(6,*) 'Warning : bare soil imposed for grid point ',ii,jj
     .              ,frac_itot
        ENDIF


C +---Define max leaf area index
C +   ==========================

        DO ll=1,nvx

         IF (NSTsvt(ii,jj,ll).eq. 0) NSTlmx(ii,jj,ll) = 0.0
         IF (NSTsvt(ii,jj,ll).eq. 1) NSTlmx(ii,jj,ll) = 0.6
         IF (NSTsvt(ii,jj,ll).eq. 2) NSTlmx(ii,jj,ll) = 0.9
         IF (NSTsvt(ii,jj,ll).eq. 3) NSTlmx(ii,jj,ll) = 1.2
         IF (NSTsvt(ii,jj,ll).eq. 4) NSTlmx(ii,jj,ll) = 0.7
         IF (NSTsvt(ii,jj,ll).eq. 5) NSTlmx(ii,jj,ll) = 1.4
         IF (NSTsvt(ii,jj,ll).eq. 6) NSTlmx(ii,jj,ll) = 2.0
         IF (NSTsvt(ii,jj,ll).eq. 7.or.NSTsvt(ii,jj,ll).eq.10)
     .        NSTlmx(ii,jj,ll) = 3.0
         IF (NSTsvt(ii,jj,ll).eq. 8.or.NSTsvt(ii,jj,ll).eq.11)
     .        NSTlmx(ii,jj,ll) = 4.5
         IF (NSTsvt(ii,jj,ll).eq. 9.or.NSTsvt(ii,jj,ll).eq.12)
     .        NSTlmx(ii,jj,ll) = 6.0

         NSTlai(ii,jj,ll) = NSTlmx(ii,jj,ll)
         NSTglf(ii,jj,ll) = 1.0

        ENDDO


C +    **************************************************************
       ELSE   ! Ocean / lake in CORINE
C +    **************************************************************

        NSTsol(ii,jj) = 1 ! Water

        NSTsvt(ii,jj,nvx)=  0
        NSTsfr(ii,jj,nvx)=100
        DO ll=1,nvx
         NSTlai(ii,jj,ll) = 0.0
         NSTglf(ii,jj,ll) = 0.0
        ENDDO


C +    **************************************************************
       ENDIF ! End of continent / water selection
       ENDIF ! End of "not undefined in CORINE" section
C +    **************************************************************



C +---Roughness 
C +   =========
C +   (computed with all CORINE data, regardless of SISVAT subgrid nvx)


        IF (FrMnul.lt.0.5.and.RUGdat) THEN  ! CORINE data available
          NST_z0(ii,jj)=0.0
          DO kk=1,nclass
           NST_z0(ii,jj)=NST_z0(ii,jj)+COR_z0(kk)*CORfrc(ii,jj,kk)
     .                                 /(1.0-FrMnul)
C +        NST_r0(ii,jj)=0.1*NST_z0(ii,jj)
          ENDDO
        ENDIF
        
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ENDDO
      ENDDO
      
      IF (nFrNul.GT.0) THEN
         write(*,*) 'CORveg (info) : '
         write(*,*) nFrNul,' points with no CORINE data'
         write(*,*)
      ENDIF
      
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END
      

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +   +++++++++++++++++++++  Minor Subroutines ++++++++++++++++++++++++++
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE AgClasses (FrAsSVT, nsvat, 
     .                      iAgCls, iDetCls, AgNeed, AgFlag)
C     
C     Purpose :
C      AgClasses checks that the number of classes whith a significant
C      fraction is not higher then the number of retained classes in
C      SISVAT (nvx). If so, the 2 provided class numbers are
C      "agregated" = the first gets the total fraction of both.
C      Aim of this = avoiding to drop an entire category of vegetation
C      simply because it is divided into "subclasses" which thus have
C      smaller fractions.
      
      INCLUDE 'NSTdim.inc'
      LOGICAL AgNeed
      INTEGER AgFlag, iAgCls, iDetCls, kk, ll, nuclas
      REAL FrAsSVT(0:nsvat), FrMin
      
      PARAMETER (FrMin=0.1)
            
      IF (AgNeed) THEN
      
C +     Count #classes abovre thereshold fraction
        nuclas= 0
        DO ll=0,nsvat
          IF (FrAsSVT(ll).GT.FrMin) THEN
             nuclas= nuclas + 1
          ENDIF
        ENDDO
      
C +     If too much used classes, aggregate
        IF (nuclas.GT.nvx) THEN
           FrAsSVT(iAgCls)= FrAsSVT(iAgCls) + FrAsSVT(iDetCls)
           FrAsSVT(iDetCls)= 0.0
           AgFlag=1
        ELSE
           AgNeed= .FALSE.
        ENDIF


      ENDIF
      END
