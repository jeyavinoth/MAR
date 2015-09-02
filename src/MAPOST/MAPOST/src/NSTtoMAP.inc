C     This is to use NESTOR domain sizes in MAPOST
C     --------------------------------------------
C     ni, nj, and nk comes from NSTdim.inc

      INTEGER LSni, LSnj, LSnk
      
      PARAMETER (LSni = ni+bi-1)
      PARAMETER (LSnj = nj+bj-1)
      PARAMETER (LSnk = nk)
      
C     This is because MAPOST can not use the begin subregion
C     indexes bi and bj, at least in the current version.
C
      INTEGER LSnk1
      PARAMETER (LSnk1 = LSnk+1)
C +   .....Note: LSnk1 = nk vertical levels + 10m OR surface values

      INTEGER mzz
      PARAMETER (mzz=mz+1)
