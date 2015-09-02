C +------------------------------------------------------------------------+
C |   NSTzz6                                          NESTOR - January 02  |
C |                                                   (created: 08/97)     |
C |   Computation of geopotential at a given pressure level.               |
C |   => NSTint will correct the NST field to obtain same Z as in LSC      |
C +------------------------------------------------------------------------+
C |                                                                        |
C |   METHOD:                                                              |
C |   ^^^^^^^                                                              |
C |     Must be called for all levels (decreasing from nk) until           |
C |     the requested pressure is reached for every horizontal position.   |
C |     (objective of this: save mem. because T,Q,p are known only in 2D)  |
C |                                                                        |
C |     Hydrostatic relation is integrated as in MAR (theta * dExner)      |
C |     (small simplification: theta assumed constant between surf-lev1)   |
C |                                                                        |
C |   INPUT:  ni,nj,nk : Grid size                                         |
C |   ^^^^^^                                                               |
C |     sp  (ni,nj)    : surface pressure  (kPa)                           |
C |     sh  (ni,nj)    : surface height    (m)                             |
C | on the k-th level:                                                     |
C |     pp  (ni,nj)    : pressure          (kPa)                           |
C |     qv  (ni,nj)    : Specific Humidity (kg/kg)                         |
C |     pt (ni,nj)     : potential temperature (K)                         |
C |                                                                        |
C |   INPUT / OUTPUT (temporary arrays):                                   |
C |   ^^^^^^^^^^^^^^^                                                      |
C |     pktv1(), pex1(),                                                   |
C |            lpres1(): retains informations for the successive calls     |
C |                      (values below current level)                      |
C |     iZp(mx,my)     : .TRUE. after the completition of all iteration    |
C |                      necessary to find the geopotential Zpl at a given |
C |                      horizontal position.                              |
C |     iZterm         : .TRUE. idem, but for completition of all grid pts.|
C |                                                                        |
C |   OUTPUT:                                                              |
C |   ^^^^^^^                                                              |
C |     Zpl            : Computed geopotential (= at the requested         |
C |                      pressure level after all iterations)              |
C +------------------------------------------------------------------------+
      SUBROUTINE NSTzz6(pt, qv, sh, sp, pp,
     .                  k, ni, nj, nk,
     .                  pktv1, pex1, lpres1, iZp, iZterm, Zpl)

      IMPLICIT NONE   

C +.. *Input and/or Output
      INTEGER k, ni, nj, nk
      REAL pt   (ni, nj), qv  (ni, nj)
      REAL sh   (ni, nj), sp  (ni, nj),  pp (ni, nj)
      REAL pktv1(ni, nj), pex1(ni, nj), lpres1(ni,nj)
      REAL Zpl  (ni, nj)
      LOGICAL iZp (ni,nj)
      LOGICAL iZterm

C +.. *Internal
      INTEGER i,j
      REAL pres,  pex, pktv, RefPL
      REAL lpres, lpres2, cpl1, cpl2

C +...*Physical constants
      REAL cp, cap, ra, grav, getpkt

      data ra    / 287.     d0/
C +...     ra    : Perfect Gas Law  Constant (J/kg/K)
      data cp    /1004.d0/
C +...     cp    : Air Specific Heat         (J/kg/K)
      data cap   /   0.28586d0/
      data grav  /   9.81   d0/

      getpkt= exp(-cap*log(100.))
C +...     getpkt: 1. / (100. (kPa) ** cap)


c +..Initialisation phase : compute functions at surface
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (k.EQ.nk) THEN
C +.. *when 1st level is requested

        DO j=1,nj
        DO i=1,ni
          lpres1(i,j)= log(sp(i,j))       
        
          pex1(i,j)  = cp *exp(cap * lpres1(i,j))
C +..     *Exner potential (Cp*p**cap)
         
          pktv1(i,j) = pt(i,j) * getpkt * (1.d0+qv(i,j)*0.608d0)
C +..     *Assume constant pkt and qv between surf. - nearest lev.
C          Please note that 0.608 is the correct coefficient,
C          as you may find from fundamental textbooks such as
C          "Triplet et Roche" (see def of virtual temperature)

          iZp(i,j)   = .FALSE.

          Zpl(i,j)   = sh(i,j)
C +..     *Begin Z integration at surface.

        ENDDO
        ENDDO

        iZterm    = .FALSE.
c +..   *Requested level is not yet reached (everywhere).

      ENDIF
C +
C +..Compute geopotential increment between levels.
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IF (.NOT.iZterm) THEN ! (work not terminated )
      iZterm= .TRUE.        ! (...but might be now )
      DO j=1,nj
      DO i=1,ni
       IF (.NOT.iZp(i,j)) THEN
 
C +..If current level (k) pressure > RefPL   pressure,
C +..compute Z of RefPL. Otherwise compute Z of level.

        RefPL = 60.           ! kPa

        pres  = pp(i,j)
        lpres2= log(pres)

        pktv=  pt(i,j) * getpkt * (1.d0+qv(i,j)*0.608d0)

        IF (pres.GT.RefPL) THEN
          iZterm   = .FALSE.  ! (some work not terminated)
        ELSE
          iZp(i,j) = .TRUE.   ! (work terminated at (i,j) )
          pres     = RefPL
          lpres    = log(pres)
C +..     *Interpolate pkt to final pres:
          cpl2 = (lpres -lpres1(i,j))/(lpres2-lpres1(i,j))
          cpl1 = (lpres2-lpres )/(lpres2-lpres1(i,j))
          pktv = cpl1*pktv1(i,j)+cpl2*pktv
        ENDIF

        pex = cp *exp(cap *log(pres))
C +..   *Exner potential (Cp*p**cap)
 
        Zpl(i,j) = Zpl(i,j) + (pex1 (i,j)-pex )
     .     *(pktv1(i,j)+pktv)*0.5d0/grav

        pktv1 (i,j)  = pktv
        pex1  (i,j)  = pex
        lpres1(i,j)  = lpres2


       ENDIF
      ENDDO
      ENDDO

      ENDIF 

      RETURN
      END
