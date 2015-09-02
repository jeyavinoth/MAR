C   +-------------------------------------------------------------------+
C   |  Subroutine CPLvgd                           January 2002 Nesting |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Creation of the vertical grid of the MAR model.                   |
C   |                                                                   |
C   | Input : - VGD_sp : surface pressure (kPa)                         |
C   | ^^^^^^^ - fID    : identificator of the Netcdf data file          |
C   |         - nz     : number of vertical levels                      |
C   |         - klev   : if specified, the level at which pressure and  |
C   |                    hybrid coordinate has to be computed           |
C   |         - parameters from CPLgrd.ctr                              |
C   |                                                                   |
C   | Output: Creation of the vertical MAR grid given in hybrid         |
C   | ^^^^^^^ coordinates :                                             |
C   |         - VGD_hp(nz+1) : local hybrid coord. for vertic. interp.  |
C   |         - VGD__p(nz+1) : pressure coordinates (kPa)               |
C   |         - VGDgdz(nz  ) : model coordinates (sigma)                |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE CPLvgd (LSCmar,klev,nz,fileID,VGD_sp,VGD_hp,
     .                   VGDgdz,VGD__p ,sigma ,WK2_1D,WK3_1D)


      IMPLICIT NONE

      INCLUDE 'CTRvar.inc'


C +---Local variables
C +   ---------------

      INTEGER k,klev,nz,k1,k2,fileID,maptyp,imez,jmez

      REAL pp1,pps,ppm,dpsl,pp,hh,ppf,GElat0,GElon0,dx,GEddxx,
     .     ptopDY,zmin,aavu,bbvu,ccvu,sst_SL,TUkhmx,empty1(1),
     .     long1,lati1,long2,lati2,long3,lati3,long4,lati4

      REAL VGD_sp,VGD_hp(nz+1),VGDgdz(nz),VGD__p(nz+1),sigma(nz),
     .     WK2_1D(nz),WK3_1D(nz)

      LOGICAL LSCmar,vertic

      CHARACTER*10 var_units


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +                          Begin (First call)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CREATION OF SIGMA MAR GRID USING LSC NetCDF FILE
C +   ================================================

      IF (LSCmar.and.lfirst_LSC) THEN


C +---Read SIGMA in NetCDF file
C +   -------------------------

C +    write(6,*) 'Caution : please verify if the UNread call'
C +    write(6,*) 'allows to extract the sigma grid.'
 
C +         ******
C +    CALL UNread (fileID,'tairDY',1,0,1,1,1,1,nz,
C +  .              empty1,empty1,sigma,var_units,WK2_1D)
C +         ******

C +    PhM -> pour Olivier:
C +     Si sigma est toujours appel2 "level" dans le fichier,
C +     la solution la plus simple est:
C +
C +         ******
       CALL UNsread (fileID,'level',0,0,
     .               0,0,nz,1,1,var_units,sigma)
C +         ******
C +


C +---Sigma coordinates
C +   -----------------

       DO k=1,nz
        VGDgdz(k)=sigma(k)
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


C +---CREATION OF SIGMA MAR GRID USING PARAMETERS IN CPLgrd.ctr
C +   =========================================================
      
      IF ((.not.LSCmar).and.lfirst_NST) THEN


C +---Read grid parameters in CPLgrd.ctr
C +   ----------------------------------

       OPEN (unit=51,status='old',file='CPLgrd.ctr')

        read (51,*) !- - - - - - - - - - - - - - - - - -
        read (51,*) !- - - - - - - - - - - - - - - - - -
        read (51,*) !- - - - - - - - - - - - - - - - - -
        read (51,*) !- - - - - - - - - - - - - - - - - -
        read (51,*) long1
        read (51,*) lati1
        read (51,*) long2
        read (51,*) lati2
        read (51,*) long3
        read (51,*) lati3
        read (51,*) long4
        read (51,*) lati4
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


C +---Print the characteristics of the vertical grid
C +   ----------------------------------------------

       write(6,*) 'Vertical CPL grid parameters'
       write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       write(6,300) nz,ptopDY
300    format(' Number of grid points : ',i4,/,
     .        ' Pressure at the top   : ',f9.4,' kPa.')
       write(6,310) zmin, aavu, bbvu, ccvu
310    format(' First level height    : ', f6.1,/,
     .        ' aavu, bbvu, ccvu      : ',(f6.1,', ',f6.1,', ',f6.1),/)


C +---Sigma coordinates
C +   -----------------

       WK3_1D(1)=ptopDY   ! Store ptopDY

       DO k=1,nz
        VGDgdz(k)=sigma(k)
       ENDDO


C +---End first call
C +   --------------

       lfirst_NST = .false.


      ENDIF


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

      ptopDY = WK3_1D(1)

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

