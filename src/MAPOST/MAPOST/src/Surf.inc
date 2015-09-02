C +-- Internal values of Surf subroutines (retain->next call):
C +   -----------------------------------
C +   * General info about var. names:
C +   tM : time-mean, tSD time Standard Deviation,
C +   _M : MAR, _R : RLS (i.e. ECMWF or other LS, int. to MAR)

C +   Precipitation 
C     - - - - - - -
      REAL lspRLS(LSni,LSnj), rr_RLS(LSni,LSnj)
      REAL cp_RLS(LSni,LSnj), sf_RLS(LSni,LSnj)
      REAL tot_RLS(LSni,LSnj)

      REAL rr_MAR(mx,my), rr_iMAR(mx,my)
      REAL tot_MAR(mx,my), tot_iMAR(mx,my)
      REAL sf_MAR(mx,my), sf_iMAR(mx,my)
      REAL cp_MAR(mx,my), cp_iMAR(mx,my)    

      common /Precip/ lspRLS,rr_RLS ,cp_RLS,sf_RLS,tot_RLS,
     .                rr_MAR,rr_iMAR,cp_MAR,cp_iMAR,
     .                sf_MAR,sf_iMAR,tot_MAR, tot_iMAR 

      REAL rr_TMAR(LSni,LSnj),sf_TMAR(LSni,LSnj),tot_TMAR(LSni,LSnj)
      REAL cp_TMAR(LSni,LSnj)

C +   Surface (2m) Temperature and albedo
C     - - - - - - - - - - - - - - - - - -
      REAL tMT2RLS(nreg),tMT2MAR(nreg)
      REAL tMTmiMAR(nreg),tMTmaMAR(nreg)

      REAL TmiMAR(mx,my), TmaMAR(mx,my)

      REAL tSD_al_M(mx,my),tM_al_M(mx,my)
      REAL tSD_al_R(mx,my),tM_al_R(mx,my) 

      common /TSurf/tMT2RLS,tMT2MAR,tMTmiMAR,tMTmaMAR,
     .              TmiMAR,TmaMAR,tSD_al_M,tM_al_M,
     .              tSD_al_R,tM_al_R
