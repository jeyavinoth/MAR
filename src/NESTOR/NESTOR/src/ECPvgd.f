C   +-------------------------------------------------------------------+
C   |  Subroutine ECPvgd                               June 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Vertical grid of the ECMWF model.                                 |
C   |                                                                   |
C   | Input : - ECP_sp : surface pressure                               |
C   | ^^^^^^^ - fID    : identificator of the Netcdf data file          |
C   |         - nk     : number of vertical levels                      |
C   |         - k      : if specified, the level at which pressure and  |
C   |                    hybrid coordinate has to be computed           |
C   |                                                                   |
C   | Output: Vertical grid of the ECMWF model :                        |
C   | ^^^^^^^ - ECP__p : pressure at each level  [kPa]                  |
C   |         - ECP__z : hybrid coordinates                             |
C   |                                                                   |
C   | Note that vertical coordinates are computed only at given (i,j)   |
C   | horizontal grid point in order to limit memory requirements.      |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE ECPvgd (fID,klev,nk,ECP_sp,ECP__p,ECP__z,plevel)

 
      IMPLICIT NONE


C +---Local variables
C +   ---------------

      INTEGER fID,k,nk,klev,k1,k2
 
      REAL pp,ppm,pps,ppf,pp1,dpsl,hh,empty1(1),plevel(nk),
     .     ECP_sp,ECP__p(nk+1),ECP__z(nk+1)

      CHARACTER*10 var_units
 

C +---Atmospheric levels: pressure levels
C +   -----------------------------------

C +        ******
      CALL UNread (fID,'p_level',0,0,0,0,nk,1,1,ECP__z,
     .             empty1,empty1,var_units,plevel)
C +        ******


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
       ECP__p(k)=plevel(k)/10.  ! (kPa)
       IF (ECP__p(k).gt.ECP_sp)
     .  ECP__p(k)=ECP_sp-REAL(k)*0.1
      ENDDO

      ECP__p(nk+1)=ECP_sp


C +---Compute hybrid coordinates (required by nesting procedure)
C +   --------------------------

      pp1  = 105.       ! Reference pressure (KPa)
      dpsl = 20.        ! "> boundary layer" (KPa)
C +...Local hybrid coordinate: set parameters:

      pps = ECP_sp
      ppm = pps - dpsl
      DO k = k1,k2+1
       pp = ECP__p(k)
       hh = pp/pp1
       IF (pp.gt.ppm) THEN
        ppf= (pp-ppm)/(pps-ppm)
        hh = hh + (pp1-pps)/pp1 * ppf * ppf
       END IF
       ECP__z(k) = LOG(hh)
      ENDDO


      RETURN
      END
