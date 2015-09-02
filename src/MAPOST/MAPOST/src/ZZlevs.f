C +------------------------------------------------------------------------+
C |   ZZlevs                                                      03/1999  |
C |   Computation of geopotential of all levels.                  v1.1r0   |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   METHOD:                                                              |
C |   ^^^^^^^                                                              |
C |     Hydrostatic relation is integrated as in MAR (theta * dExner)      |
C |     (small simplification: theta assumed constant between surf-lev1)   |
C |                                                                        |
C |     This is a "simplified" version of the Zplev routine from NESTOR.   |
C |     (i.e. there is no interpolation to pressure levels, but the method |
C |     is exactly the same).                                              |
C |                                                                        |
C |   INPUT:  ni,nj,nk : Grid size                                         |
C |   ^^^^^^                                                               |
C |     CSTp(nk+1)     : pressure is given in the general form             |
C |     SIGp(nk+1)       pst(i,j)*SIGp(k)+CSTp(k)                          |
C |     pst(ni,nj)       k= nk+1 at surface; pst units: kPa.               |
C |                                                                        |
C |     qv (ni,nj,nk)  : Specific Humidity               (kg/kg)           |
C |     pkt(ni,nj,nk)  : potential temperature divided by 100.[kPa]**(R/Cp)|
C |     sh (ni,nj)     : surface height (m)                                |
C |                                                                        |
C |   OUTPUT:                                                              |
C |   ^^^^^^^                                                              |
C |     zz_out(ni,nj,nk): Computed geopotential                            |
C +------------------------------------------------------------------------+
      SUBROUTINE ZZlevs(pkt, qv, sh, pst, CSTp, SIGp,
     .                  ni, nj, nk, zz_out)

      IMPLICIT NONE   
      
      INCLUDE 'LSMphy.inc'

C +.. *Input / Output
      INTEGER k, ni, nj, nk
      REAL pkt  (ni, nj, nk), qv  (ni, nj, nk)
      REAL pktv1, pex1
      REAL sh   (ni, nj), pst (ni, nj)
      REAL zz_out  (ni, nj, nk)
      REAL CSTp(nk+1), SIGp(nk+1)
   
C +.. *Internal
      INTEGER i,j
      REAL pres, pres1, pex, pktv, zzcalc
      

      DO j=1,nj
      DO i=1,ni

c +..Initialisation phase : compute functions at surface
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       pex1 = cp *exp(cap *log(pst(i,j)+CSTp(nk+1)))
C +..  *Exner potential (Cp*p**cap)
         
       pktv1= pkt(i,j,nk)*(1.d0+qv(i,j,nk)*0.608d0)
C +..  *Assume constant pkt and qv between surf. - nearest lev.
C      (0.608 -> 0.85 dans MAR)

       zzcalc= sh(i,j)
C +..  *Begin Z integration at surface.

C +
C +..Compute geopotential (integrate hydrostat. relation)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       DO k=nk,1,-1
 
        pres= pst(i,j)*SIGp(k)+CSTp(k)

C +..   *Exner potential (Cp*p**cap):
        pex = cp *exp(cap *log(pres))

        pktv=  pkt(i,j,k)*(1.d0+qv(i,j,k)*0.608d0)
C      (0.608 -> 0.85 dans MAR)

        zzcalc= zzcalc + (pex1-pex)
     .        *(pktv1+pktv)*0.5d0/grav

        pktv1 = pktv
        pex1  = pex

C +..   *output Z of level:
        zz_out(i,j,k) = zzcalc 

       ENDDO
      ENDDO
      ENDDO 

      RETURN
      END
