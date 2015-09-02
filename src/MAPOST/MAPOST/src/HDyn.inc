C +-- Internal values of HDyn subroutines (retain->next call):
C +   -----------------------------------
C +   * General info about var. names:
C +   tM : time-mean, tSD time Standard Deviation,
C +   _M : MAR, _R : RLS (i.e. ECMWF or other LS, int. to MAR)
C +
C +   Version: 3.0 (2002)

C +   z:  geopotential on p-levels:
C     - - - - - - - - - - - - - - - 
      REAL tMz_M(mx,my,npl),tSDz_M(mx,my,npl)
      REAL tMz_R(mx,my,npl),tSDz_R(mx,my,npl)

C +.. On 500 hPa level only:
      REAL SDBPz_R(mx,my), SDBPz_M(mx,my)
      REAL MBPz_R(mx,my), MBPz_M(mx,my)

C +.. Temporary arrays for the time-filter:
      REAL    FilMAR(mx,my,nfl1), FilRLS(mx,my,nfl1)
      INTEGER ipfMAR, ipfRLS

      common /z_plev/ tMz_M, tSDz_M, tMz_R, tSDz_R,
     .                MBPz_R, MBPz_M, SDBPz_R, SDBPz_M,
     .                FilMAR, FilRLS, ipfMAR, ipfRLS

C +   slp: Mean Sea Level Pressure:
C     - - - - - - - - - - - - - - -
      REAL tMslpM(mx,my),tSDslpM(mx,my)
      REAL tMslpR(mx,my),tSDslpR(mx,my)
      REAL tMslpBIA(nreg), tMslpRME(nreg)
 
      common /mslp/ tMslpM,tSDslpM,tMslpR,tSDslpR,tMslpBIA,tMslpRME

C +   Temperature and water vapor:
C     - - - - - - - - - - - - - - -
      REAL tMT_M(mx,my,npl),tSDT_M(mx,my,npl)
      REAL tMT_R(mx,my,npl),tSDT_R(mx,my,npl)

      REAL tMQ_M(mx,my,npl),tSDQ_M(mx,my,npl)
      REAL tMQ_R(mx,my,npl),tSDQ_R(mx,my,npl)

      common /TQplev/ tMT_M, tSDT_M, tMT_R, tSDT_R,
     .                tMQ_M, tSDQ_M, tMQ_R, tSDQ_R

C +   Wind Vector:
C     - - - - - - 
      REAL tMu_M(mx,my,npl),tSDu_M(mx,my,npl)
      REAL tMu_R(mx,my,npl),tSDu_R(mx,my,npl)
      REAL tMv_M(mx,my,npl),tSDv_M(mx,my,npl)
      REAL tMv_R(mx,my,npl),tSDv_R(mx,my,npl)

      common /WindVec/ tMu_M,tSDu_M,tMu_R,tSDu_R,
     .                 tMv_M,tSDv_M,tMv_R,tSDv_R

C +   Wind Norm 500:
C     - - - - - - - -
      REAL tMuvM(mx,my),tSDuvM(mx,my)
      REAL tMuvR(mx,my),tSDuvR(mx,my)

      common /Wind500/ tMuvM, tSDuvM, tMuvR, tSDuvR


