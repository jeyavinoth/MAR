C     The purpose of this file is to give MAPOST "modules"
C     access to MAPOST global data whithout needing an
C     explicit argument.
C
C     If you do like to share everything using arguments,
C     you can still do it !
C     These variables should not be modified inside the modules,
C     and they should not be given other names
C     (these things do not make sense for non-MAR code,
C       you probably have no reason to rename)


      CHARACTER*3      LSCmod, REGmod, REGstudy
      common/MODnames/ LSCmod, REGmod, REGstudy
      
      