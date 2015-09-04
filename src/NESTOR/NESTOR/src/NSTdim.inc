
C +---NST domain dimensions
C +   ---------------------

      INTEGER mx,my,mz,mzabso,mw,nvx,nsl,nbdom,nsno

      PARAMETER (mx    = 160)   ! X-grid
      PARAMETER (my    = 160)   ! Y-grid
      PARAMETER (mz    =  23)   ! Z-grid 
      PARAMETER (mzabso=   4)   ! Z-grid 

      PARAMETER (nvx   =   2)   ! Sub-division of grid cell (SVAT)
      PARAMETER (mw    =   2)   ! Sub-division of grid cell (Deardorff)

      PARAMETER (nsl   =   7)   ! Soil layers               (SVAT)
      PARAMETER (nsno  =  25)   ! Snow layers               (Snow model)
      PARAMETER (nbdom =   3)   ! Number of continents      ("GLOveg.f")


C +---LSC domain dimensions
C +   ---------------------

      INTEGER ni,nj,njv,nk,bi,bj,isLMz

C     A sub-region of the external large-scale domain is defined in order to 
C     reduced the CPU cost and the memory requirement for the interpolation.

C +---1. SIZE of the SUB-REGION (in grid points)

      PARAMETER (isLMz = 0)

      PARAMETER (ni = 144)
      PARAMETER (nj =  73)
      PARAMETER (njv=  nj-isLMz)
      PARAMETER (nk =  18)

C     Warning:
C     For LMDz, you may use the scalar grid size, nj= size(lat_s) 
C     but in that case, you must set   isLMz=1 (njv = nj-1)                     
C     For all other models, please set isLMz=0 (njv = nj  ) 

C +---2. BEGINNING INDEX of the SUB-REGION

      PARAMETER (bi = 1)
      PARAMETER (bj = 1)

C +---Selector for vectorization of the MAR code
C +   ------------------------------------------

      LOGICAL vector

      PARAMETER (vector = .false.)

C     "vector" is true only if the MAR code is run on vectorial computer



C +---Dimensions of the RELAXATION ZONE towards LATERAL BOUNDARIES
C     ------------------------------------------------------------

      INTEGER n6,n7,n8,n9,n10 
      PARAMETER(n6 = 6)
C     ......... ^ number of grid points of the relaxation zone

      PARAMETER(n10=5) 
C     ......... ^ number of grid points of constant topo. zone

      PARAMETER(n8 = 3)
C     ......... ^ number of grid points of the topography
C                 transition zone (valid if using LS constant
C                 topography at boundaries).

      PARAMETER(n7 = n6+1)
      PARAMETER(n9 =n10-1)


C     Explanation of boundary structure :
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     1. TOPOGRAPHY
C     -------------
C
C        |   Constant  | Transition |  Computation  | Transition |  Constant  |
C        |  topography |    zone    |    domain     |    zone    | topography |
C        |     zone    | (LS -> MAR)|               | (LS -> MAR)|    zone    |
C        ^             ^            ^               ^            ^            ^
C        1    ...     n10  ...  n10+n8+1  ...  mx-n9-n8-1 ...  mx-n9   ...   mx
C
C     2. RELAXATION LSC --> NST
C     -------------------------
C
C        |      Relaxation     |      Computation      |      Relaxation      |
C        |         zone        |        domain         |         zone         |
C        ^                     |                       |                      ^
C        1         ...    n6-1 |         ...           |mx-n6+2   ...        mx

