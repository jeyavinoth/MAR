C   +-------------------------------------------------------------------+
C   |  Subroutine CM3vgd                           January 2002 NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Vertical grid of the hadCM3/ECHAM5/CanESM2/NorESM1 model.         |
C   |                                                                   |
C   | Input : - ECM_sp : surface pressure                               |
C   | ^^^^^^^ - fID    : identificator of the Netcdf data file          |
C   |         - nk     : number of vertical levels                      |
C   |         - klev   : if specified, the level at which pressure and  |
C   |                    hybrid coordinate has to be computed           |
C   |                                                                   |
C   | Output: Vertical grid of the ECMWF model :                        |
C   | ^^^^^^^ - ECM__p : pressure at each level  [kPa]                  |
C   |         - ECM_hp : local hybrid coord. for vertic. interpolation. |
C   |                                                                   |
C   | Note that vertical coordinates are computed only at given (i,j)   |
C   | horizontal grid point in order to limit memory requirements.      |
C   |                                                                   |
C   +-------------------------------------------------------------------+
      SUBROUTINE CM3vgd (fID,klev,ECM_sp,ECM__p,ECM_hp)

 
      IMPLICIT NONE


C +---Local variables
C +   ---------------
      INCLUDE 'NSTdim.inc'

      INTEGER fID,k,klev,k1,k2,k21

      REAL pp,ppm,pps,ppf,pp1,dpsl,hh,empty1(1),CSTp(nk),SIGp(nk),
     .     ECM_sp,ECM__p(nk+1),ECM_hp(nk+1) 
     
      SAVE CSTp, SIGp

      CHARACTER*10 var_units
      
      LOGICAL lfirst
      SAVE    lfirst
      DATA    lfirst/.true./

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +                          Begin (First call)
      IF (lfirst) THEN      
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Atmospheric levels: pressure levels
C +   -----------------------------------

C +          ******
        CALL UNsread (fID,'CSTp',1,1,
     &                1,1,nk,1,1,var_units,CSTp)
C +          ******

C +          ******
        CALL UNsread (fID,'SIGp',1,1,
     &                1,1,nk,1,1,var_units,SIGp)
C +          ******

C +---Adapt units
C +   -----------

        DO k=1,nk
          CSTp(k) = CSTp(k) * 1.E-3  !(Pa-->KPa)
        ENDDO
        
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      lfirst = .false.
      ENDIF                  ! End (First call)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 

C +---Computation for a given level or all levels ?
C +   ---------------------------------------------

      IF ((klev.le.0).or.(klev.gt.(nk+1))) THEN
       k1 =1
       k2 =nk
       k21=nk+1
      ELSE
       k1 =klev
       k2 =klev
       k21=klev
      ENDIF


C +---Compute pressure at each levels
C +   -------------------------------

C +...Pressure in LS atmosphere is such that :
C +...p(level) = CSTp(level) + SIGp(level) * Surf_pressure

      IF (klev.ne.(nk+1)) THEN
        DO k=k1,k2
          ECM__p(k)=CSTp(k)+SIGp(k)*ECM_sp  ! (kPa)
        ENDDO
      ENDIF

      ECM__p(nk+1)=ECM_sp


C +---Compute hybrid coordinates (required by nesting procedure)
C +   --------------------------

      pp1  = 105.       ! Reference pressure (KPa)
      dpsl = 20.        ! "> boundary layer" (KPa)
C +...Local hybrid coordinate: set parameters

      pps = ECM_sp
      ppm = pps - dpsl
      DO k = k1,k21
       pp = ECM__p(k)
       hh = pp/pp1
       IF (pp.gt.ppm) THEN
        ppf= (pp-ppm)/(pps-ppm)
        hh = hh + (pp1-pps)/pp1 * ppf * ppf
       END IF
       ECM_hp(k) = LOG(hh)
      ENDDO

      RETURN
      END
