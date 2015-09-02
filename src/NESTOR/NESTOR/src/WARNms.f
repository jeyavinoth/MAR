C   +-------------------------------------------------------------------+
C   |  Subroutine WARNms                               July 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : Some variables which have to be carrefully specified.     |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: Warnings to prevent inadapted use of the code.            |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE WARNms


      IMPLICIT NONE


C +---General variables
C +   -----------------

      include 'NSTdim.inc'
      include 'LSCvar.inc'
      include 'SNDvar.inc'
      include 'NESTOR.inc'

C +---Horizontal interpolation
C +   ------------------------

      IF (SPHgrd.and.(HORint.gt.1)) THEN
       WRITE(6,*) 'Warning : first order for horizontal interpolation'
       WRITE(6,*) '~~~~~~~   is recommended if the simulation domain '
       WRITE(6,*) '          is cyclic or includes North/South pole. '
       WRITE(6,*)
      ENDIF

C +---Horizontal LSC grid    
C +   -------------------

      IF (LSCmod.NE.'LMz'.AND.nj.NE.njv) THEN
       WRITE(6,*) 'Error   : Grid size (NSTdim.inc) seems wrong.     '
       WRITE(6,*) '~~~~~     For all LSC input except LMDz (LMz),    '
       WRITE(6,*) '          njv must be equal to nj                 '
       WRITE(6,*)
       STOP
      ENDIF


C +---Vertical interpolation
C +   ----------------------

      IF (VERint.eq.3) THEN
       WRITE(6,*) 'Warning : third order for vertical interpolation  '
       WRITE(6,*) '~~~~~~~   is not recommended since it could induce'
       WRITE(6,*) '          strange variations in vertical profiles.'
       WRITE(6,*)
      ENDIF


C +---Correction of 600-hPa geopotential
C +   ----------------------------------

      IF (CORzz6) THEN
       WRITE(6,*) 'Note : 600hPa-based correction activated'
       WRITE(6,*) '~~~~   (NST height = LSC height at 600 hPa)'
       WRITE(6,*)
      ELSE
       WRITE(6,*) 'Note : 600hPa-based correction NOT activated'
       WRITE(6,*) '~~~~   (bad idea, at least if you have'
       WRITE(6,*) '        mountains near the boundaries)'
       WRITE(6,*)
      ENDIF


C +---NDVI Databases
C +   --------------

      IF (NDV1km.and.NDV8km) THEN
       WRITE(6,*) 'NDVI databases : select either 1-km resolution or'
       WRITE(6,*) '~~~~~~~~~~~~~~   8-km resolution in NSTing.ctr file'
       WRITE(6,*)
       WRITE(6,*) 'STOP.'
       WRITE(6,*)
      ENDIF


C +---Prognostic variables of SVAT
C +   ----------------------------

      IF (SNDing.and.SVTlsc) THEN
       WRITE(6,*) 'Warning : sounding and soil wetness estimated '
       WRITE(6,*) '~~~~~~~   from ECMWF fields are not compatible. '
       WRITE(6,*) '          Imposed relative wetness in all layers'
       WRITE(6,*) '          is then considered. '
       WRITE(6,*)
      ENDIF


      RETURN
      END

