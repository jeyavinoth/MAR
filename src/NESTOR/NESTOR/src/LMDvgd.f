C   +-------------------------------------------------------------------+
C   |  Subroutine LMDvgd                          October 2000  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Vertical grid of the LMDZ model.                                  |
C   |                                                                   |
C   | Input : - LMD_sp : surface pressure                               |
C   | ^^^^^^^ - fID    : identificator of the Netcdf data file          |
C   |         - nk     : number of vertical levels                      |
C   |         - k      : if specified, the level at which pressure and  |
C   |                    hybrid coordinate has to be computed           |
C   |                                                                   |
C   | Output: Vertical grid of the LMDZ model :                         |
C   | ^^^^^^^ - LMD__p : pressure at each level  [kPa]                  |
C   |         - LMD__z : hybrid coordinates                             |
C   |                                                                   |
C   | Note that vertical coordinates are computed only at given (i,j)   |
C   | horizontal grid point in order to limit memory requirements.      |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE LMDvgd (fID,klev,nk,LMD_sp,LMD__p,LMD__z,LSC__z,sigma)

 
      IMPLICIT NONE


C +---Local variables
C +   ---------------

      INTEGER fID,k,nk,klev,k1,k2

      REAL pp,ppm,pps,ppf,pp1,dpsl,hh,empty1(1),pres_nivs(nk),
     .     LMD_sp,LMD__p(nk+1),LMD__z(nk+1),LSC__z(nk),sigma(nk)

      CHARACTER*10 var_units
 

C +---Atmospheric levels: pressure levels
C +   -----------------------------------

C +        ******
      CALL UNread (fID,'pres_nivs',0,0,0,0,nk,1,1,LMD__z,
     .             empty1,empty1,var_units,pres_nivs)
      CALL UNread (fID,'pres',1,k,1,1,1,1,nk,
     .             empty1,empty1,LSC__z,var_units,pres_nivs)

C +        ******


C +---Compute sigma levels
C +   --------------------

      DO k=1,nk
       sigma(k)=pres_nivs(k)/LMD_sp
      ENDDO


C +---Computation for a given level or all levels ?
C +   ---------------------------------------------

      IF ((klev.le.0).or.(klev.gt.nk)) THEN
       k1=1
       k2=nk
      ELSE
       k1=1
       k2=klev
      ENDIF


C +---Get pressure at each levels and adapt units
C +   -------------------------------------------

      DO k=k1,k2
       LMD__p(k)=sigma(k)*LMD_sp * 1.E-3  ! (kPa)
      ENDDO

      LMD__p(nk+1)=LMD_sp


C +---Compute hybrid coordinates (required by nesting procedure)
C +   --------------------------

      pp1  = 105.       ! Reference pressure (KPa)
      dpsl = 20.        ! "> boundary layer" (KPa)
C +...Local hybrid coordinate: set parameters:

      pps = LMD_sp
      ppm = pps - dpsl
      DO k = k1,k2+1
       pp = LMD__p(k)
       hh = pp/pp1
       IF (pp.gt.ppm) THEN
        ppf= (pp-ppm)/(pps-ppm)
        hh = hh + (pp1-pps)/pp1 * ppf * ppf
       END IF
       LMD__z(k) = LOG(hh)
      ENDDO


      RETURN
      END
