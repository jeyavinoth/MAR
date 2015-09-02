C   +-------------------------------------------------------------------+
C   |  Subroutine VERhyb                               May 2002 NESTING |
C   +-------------------------------------------------------------------+
C   | Computes local hybrid coordinate used for vertical interpolation. |
C   |                                                                   |
C   | Input :                                                           |
C   | ^^^^^^^ - nk           : number of vertical levels                |
C   |         - LSC_sp       : surface pressure                         |
C   |         - LSC__p(nk+1) : pressure on the levels                   |
C   |                                                                   |
C   | Output:                                                           |
C   | ^^^^^^^                                                           |
C   |         - LSC_hp(nk+1) : local hybrid coord. for vertic. interp.  |
C   |                                                                   |
C   | Note that vertical coordinates are computed only at given (i,j)   |
C   | horizontal grid point in order to limit memory requirements.      |
C   |                                                                   |
C   +-------------------------------------------------------------------+
      SUBROUTINE VERhyb (nk,LSC_sp,LSC__p,LSC_hp)

 
      IMPLICIT NONE

      INTEGER k,nk

      REAL pp,ppm,pps,ppf,pp1,dpsl,hh,
     .     LSC_sp,LSC__p(nk+1),LSC_hp(nk+1) 
     

      pp1  = 105.       ! Reference pressure (KPa)
      dpsl = 20.        ! "> boundary layer" (KPa)
C +...Local hybrid coordinate: set parameters

      pps = LSC_sp
      ppm = pps - dpsl
      DO k = 1,nk+1
       pp = LSC__p(k)
       hh = pp/pp1
       IF (pp.gt.ppm) THEN
        ppf= (pp-ppm)/(pps-ppm)
        hh = hh + (pp1-pps)/pp1 * ppf * ppf
       END IF
       LSC_hp(k) = LOG(hh)
      ENDDO

      RETURN
      END
