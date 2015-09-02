C   +-------------------------------------------------------------------+
C   |  Subroutine NCPvgd                         December 2001  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Vertical grid of the NCEP model.                                  |
C   |                                                                   |
C   | Input : - NCP_sp : surface pressure                               |
C   | ^^^^^^^ - fID    : identificator of the Netcdf data file          |
C   |         - nk     : number of vertical levels                      |
C   |         - k      : if specified, the level at which pressure and  |
C   |                    hybrid coordinate has to be computed           |
C   |                                                                   |
C   | Output: Vertical grid of the ECMWF model :                        |
C   | ^^^^^^^ - NCP__p : pressure at each level  [kPa]                  |
C   |         - NCP__z : hybrid coordinates                             |
C   |                                                                   |
C   | Note that vertical coordinates are computed only at given (i,j)   |
C   | horizontal grid point in order to limit memory requirements.      |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE NCPvgd (fID,klev,nk,NCP_sp,NCP__p,NCP__z,plevel)

 
      IMPLICIT NONE


C +---Local variables
C +   ---------------

      INTEGER fID,k,nk,klev,k1,k2
 
      REAL pp,ppm,pps,ppf,pp1,dpsl,hh,empty1(1),plevel(nk),
     .     NCP_sp,NCP__p(nk+1),NCP__z(nk+1)

      CHARACTER*10 var_units
 

C +---Atmospheric levels: pressure levels
C +   -----------------------------------

C +        ******
      CALL UNread (fID,'level',0,0,0,0,nk,1,1,NCP__z,
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
       NCP__p(k)=plevel(k)/10.  ! (kPa)
c       IF (NCP__p(k).gt.NCP_sp)
c     .  NCP__p(k)=NCP_sp-REAL(k)*0.1
      ENDDO

      NCP__p(nk+1)=NCP_sp


C +---Compute hybrid coordinates (required by nesting procedure)
C +   --------------------------

      pp1  = 105.       ! Reference pressure (KPa)
      dpsl = 20.        ! "> boundary layer" (KPa)
C +...Local hybrid coordinate: set parameters:

      pps = NCP_sp
      ppm = pps - dpsl
      DO k = k1,k2+1
       pp = NCP__p(k)
       hh = pp/pp1
       IF (pp.gt.ppm) THEN
        ppf= (pp-ppm)/(pps-ppm)
        hh = hh + (pp1-pps)/pp1 * ppf * ppf
       END IF
       NCP__z(k) = LOG(hh)
      ENDDO

      RETURN
      END
