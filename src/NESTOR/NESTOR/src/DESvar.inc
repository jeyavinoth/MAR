
C +---Relation between desagregated grid and input grid
C +   =================================================

      INTEGER iiL2N (mx,my),jjL2N (mx,my)

      REAL    auxL2N(mx,my),auyL2N(mx,my)

      INTEGER fdiv


C +...fdiv   : factor between input and output resolutions
C +...iiL2N  : location of (i,j) NST grid point in LSC grid (X-coord)
C +...jjL2N  : location of (i,j) NST grid point in LSC grid (Y-coord)
C +...auxL2N : interpolation factors between NST and LSC grid (X)
C +...auyL2N : interpolation factors between NST and LSC grid (Y)


      COMMON/DESvar_i/iiL2N , jjL2N,fdiv
      
      COMMON/DESvar_r/auxL2N,auyL2N
