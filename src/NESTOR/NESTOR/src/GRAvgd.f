C   +-------------------------------------------------------------------+
C   |  Subroutine GRAvgd                         February 2002  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Vertical grid for GRADS output analysis.                          |
C   |                                                                   |
C   | Input : - GRA_sp : surface pressure                               |
C   | ^^^^^^^ - fID    : identificator of the Netcdf data file          |
C   |         - nk     : number of vertical levels                      |
C   |         - k      : if specified, the level at which pressure and  |
C   |                    hybrid coordinate has to be computed           |
C   |                                                                   |
C   | Output: Vertical grid of the ECMWF model :                        |
C   | ^^^^^^^ - GRA__p : pressure at each level  [kPa]                  |
C   |         - GRA_hp : hybrid coordinates                             |
C   |                                                                   |
C   | Note that vertical coordinates are computed only at given (i,j)   |
C   | horizontal grid point in order to limit memory requirements.      |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE GRAvgd (fID,klev,nk,GRA_sp,GRA__p,GRA_hp,plevel)

 
      IMPLICIT NONE


C +---Local variables
C +   ---------------

      INTEGER fID,k,nk,klev,k1,k2,nkk
 
      REAL pp,ppm,pps,ppf,pp1,dpsl,hh,empty1(1),plevel(nk),
     .     GRA_sp,GRA__p(nk+1),GRA_hp(nk+1)

      CHARACTER*10 var_units
 

C +---Atmospheric levels: pressure levels
C +   -----------------------------------

      nkk = nk
      IF (nkk.ne.12) THEN
       write(6,*)
       write(6,*) 'GRADS output grid is valid only with 12 vertical'
       write(6,*) 'levels. Please set mz=12 in NSTdim.inc'
       write(6,*)
       write(6,*) '--> STOP in GRAvgd.f'
       write(6,*)
       STOP
      ENDIF

      plevel( 1) = 100.
      plevel( 2) = 150.
      plevel( 3) = 200.
      plevel( 4) = 250.
      plevel( 5) = 300.
      plevel( 6) = 400.
      plevel( 7) = 500.
      plevel( 8) = 600.
      plevel( 9) = 700.
      plevel(10) = 850.
      plevel(11) = 925.
      plevel(12) = 1000.


C +---Computation for a given level or all levels ?
C +   ---------------------------------------------

      IF ((klev.le.0).or.(klev.gt.nk)) THEN
       k1=1
       k2=nk
      ELSE
       k1=1
       k2=klev
      ENDIF


C +---Compute pressure at each levels
C +   -------------------------------

      DO k=k1,k2
       GRA__p(k)=plevel(k)/10.  ! (kPa)
ccccc  IF (GRA__p(k).gt.GRA_sp)
ccccc.  GRA__p(k)=GRA_sp-REAL(k)*0.1
      ENDDO

      GRA__p(nk+1)=GRA_sp


C +---Compute hybrid coordinates (required by nesting procedure)
C +   --------------------------

      pp1  = 105.       ! Reference pressure (KPa)
      dpsl = 20.        ! "> boundary layer" (KPa)
C +...Local hybrid coordinate: set parameters:

      pps = GRA_sp
      ppm = pps - dpsl
      DO k = k1,k2+1
       pp = GRA__p(k)
       hh = pp/pp1
       IF (pp.gt.ppm) THEN
        ppf= (pp-ppm)/(pps-ppm)
        hh = hh + (pp1-pps)/pp1 * ppf * ppf
       END IF
       GRA_hp(k) = LOG(hh)
      ENDDO


      RETURN
      END
