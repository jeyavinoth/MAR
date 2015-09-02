C   +-------------------------------------------------------------------+
C   |  Subroutine VERgrd                               June 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Creation of the vertical grid of a given model.                   |
C   |                                                                   |
C   | Input : - VGD_sp : surface pressure                               |
C   | ^^^^^^^ - VGD__t : real temperature                               |
C   |         - fID    : identificator of the Netcdf data file          |
C   |         - nz     : number of vertical levels                      |
C   |         - k      : if specified, the level at which pressure and  |
C   |                    hybrid coordinate has to be computed           |
C   |                                                                   |
C   | Output: Vertical grid of the LSC model :                          |
C   | ^^^^^^^ - VGD__p : pressure at each level  [kPa]                  |
C   |         - VGD_hp : local hybrid coord. for vertic. interpolation. |
C   |         - WK2_1D : level heights (if computed)                    |
C   |                                                                   |
C   | Note that vertical coordinates are computed only at given (i,j)   |
C   | horizontal grid point in order to limit memory requirements.      |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE VERgrd (LSCmod,NSTmod,fID,k,nz,VGD_sp,
     .                   VGD_sh,VGD__t,VGD__p,VGD_hp,
     .                   VGDgdz,WK1_1D,WK2_1D,WK3_1D)


 
      IMPLICIT NONE

      INCLUDE 'CTRvar.inc'


C +---Local variables
C +   ===============

      INTEGER fID,k,nz
 
      REAL VGD_sp,VGD_sh,VGD__t(nz),VGD__p(nz+1),VGD_hp(nz+1),
     .     WK1_1D(nz),WK2_1D(nz),WK3_1D(nz),VGDgdz(nz)

      CHARACTER*3 LSCmod,NSTmod

      LOGICAL Vtrue,Vfalse


C +---Data
C +   ----

      DATA Vtrue  /  .true. /
      DATA Vfalse / .false. /
 

C +---Creation of the vertical grid depending on the specified model
C +   ==============================================================


C +---European Center of Medium-Range Forecast (ECMWF)
C +   ------------------------------------------------

C +---Hybrid levels
C +   - - - - - - -

      IF (LSCmod.eq.'ECM' .OR. LSCmod.eq.'LMz'.OR.
     .    LSCmod.eq.'E15' .OR. LSCmod.eq.'E40')

C +         ******
     . CALL ECMvgd (fID,k,VGD_sp,VGD__p,VGD_hp)
C +         ******

      IF (LSCmod.eq.'GCM')

C +         ******
     . CALL CM3vgd (fID,k,VGD_sp,VGD__p,VGD_hp)
C +         ******

C +---Pressure levels
C +   - - - - - - - -

      IF (LSCmod.eq.'ECP')

C +         ******
     . CALL ECPvgd (fID,k,nz,VGD_sp,VGD__p,VGD_hp,WK1_1D)
C +         ******


C +---NCEP analysis
C +   -------------

      IF (LSCmod.eq.'NCP') 

C +         ******
     . CALL NCPvgd (fID,k,nz,VGD_sp,VGD__p,VGD_hp,WK1_1D)
C +         ******


C +---Laboratoire de Meteorologie Dynamique (LMD)
C +   -------------------------------------------

      IF (LSCmod.eq.'LMD') 

C +         ******
     . CALL LMDvgd (fID,k,nz,VGD_sp,VGD__p,VGD_hp,WK1_1D,WK2_1D)
C +         ******


C +---Modele Atmospherique Regional (MAR)
C +   -----------------------------------

      IF (LSCmod.eq.'MAR') 

C +         ******
     . CALL MARvgd (Vtrue ,k,nz,fID,VGD_sp,VGD_hp,
     .              VGDgdz,VGD__p,WK1_1D,WK2_1D,WK3_1D)
C +         ******

      IF (NSTmod.eq.'MAR'.or.NSTmod.eq.'M2D')

C +         ******
     . CALL MARvgd (Vfalse,k,nz,fID,VGD_sp,VGD_hp,
     .              VGDgdz,VGD__p,WK1_1D,WK2_1D,WK3_1D)
C +         ******


C +---GRADS output analysis
C +   ---------------------

      IF (NSTmod.eq.'GRA') 

C +         ******
     . CALL GRAvgd (fID,k,nz,VGD_sp,VGD__p,VGD_hp,WK1_1D)
C +         ******


C +---Output for SVAT coupling
C +   ------------------------

      IF (NSTmod.eq.'CPL') 

C +         ******
     . CALL CPLvgd (Vfalse,k,nz,fID,VGD_sp,VGD_hp,
     .              VGDgdz,VGD__p,WK1_1D,WK2_1D,WK3_1D)
C +         ******


      RETURN
      END
