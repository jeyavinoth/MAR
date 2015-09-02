C +------------------------------------------------------------------------+
C |   VERhyd                                      NESTOR/MAPOST January 02 |
C |   Computation of geopotential height by HYDrostatic relation           |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   METHOD:                                                              |
C |   ^^^^^^^                                                              |
C |     Hydrostatic relation is integrated as in MAR (theta * dExner)      |
C |     (small simplification: theta assumed constant between surf-lev1    |
C |      unless a T and q is given for the surface)                        |
C |                                                                        |
C |     This is a "simplified" version of the Zplev routine from LSMARIN   |
C |     (i.e. there is no interpolation to pressure levels, but the method |
C |     is exactly the same).                                              |
C |                                                                        |
C |   INPUT:  ni,nj,nk   : Grid size                                       |
C |   ^^^^^^               Note: nk is 1st level above surface             |
C |                              OR the surface itself                     |
C |                                                                        |
C |                                                                        |
C |     MOD__p(ni,nj,nk) : pressure                                        |
C |     MOD_sp(ni,nj)    : surface pressure; MOD_sp units: kPa.            |
C |                                                                        |
C |     MOD_qv(ni,nj,nk) : Specific Humidity               (kg/kg)         |
C |     MOD_pt(ni,nj,nk) : potential temp.                                 |
C |     getpkt           : If MOD_pt = true potential temp,                |
C |                        set getpkt to 100**(-cap). For mar, set getpkt=1|
C |     MOD_sh(ni,nj)    : surface height (m)                              |
C |                                                                        |
C |   OUTPUT:                                                              |
C |   ^^^^^^^                                                              |
C |     MOD_zz(ni,nj,nk): Computed geopotential                            |
C |                       NOTE: levels with p < sp                         |
C |                             are set to MOD_zz=MOD_sh                   |
C +------------------------------------------------------------------------+
      SUBROUTINE VERhyd(MOD_pt, MOD_qv, MOD_sh, MOD_sp, MOD__p,
     .                  getpkt, ni, nj, nk, MOD_zz)

      IMPLICIT NONE   
      
      INCLUDE 'NSTphy.inc'

C +.. *Input / Output
      INTEGER ni, nj, nk
      REAL MOD__p  (ni, nj, nk)
      REAL MOD_pt  (ni, nj, nk), MOD_qv (ni, nj, nk)
      REAL pktv1, pex1
      REAL MOD_sh  (ni, nj), MOD_sp (ni, nj)
      REAL MOD_zz  (ni, nj, nk)
      REAL getpkt
   
C +.. *Internal
      INTEGER i,j,k
      REAL pex, pktv, zzcalc
      

      DO j=1,nj
      DO i=1,ni

c +..Initialisation phase : compute functions at surface
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C +..   *Exner potential (Cp*p**cap):
       pex1 = cp*exp(cap*log(MOD_sp(i,j)))
C      Note: use of surface pressure is ok for both MAR and 
C            data defined on hybrid or constant p-levels.

       pktv1= MOD_pt(i,j,nk)*getpkt*(1.d0+MOD_qv(i,j,nk)*0.608d0)
C +..  *Assume constant MOD_pt and MOD_qv between surf. - nearest lev.

       zzcalc= MOD_sh(i,j)
C +..  *Begin Z integration at surface.

C +
C +..Compute geopotential (integrate hydrostat. relation)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       DO k=nk,1,-1
 
C +..   *Exner potential (Cp*p**cap):
        pex = cp *exp(cap *log(MOD__p(i,j,k)))

        pktv=  MOD_pt(i,j,k)*getpkt*(1.d0+MOD_qv(i,j,k)*0.608d0)
C      (0.608 -> 0.85 dans MAR;
C         mailto:philippe.marbaix@advalvas.be for info)

        IF (pex1.GT.pex) THEN
          zzcalc= zzcalc + (pex1-pex)
     .        *(pktv1+pktv)*0.5d0/grav
          pex1  = pex
        ENDIF

        pktv1 = pktv

C +..   *output Z of level:
        MOD_zz(i,j,k) = zzcalc 

       ENDDO
      ENDDO
      ENDDO 

      RETURN
      END
