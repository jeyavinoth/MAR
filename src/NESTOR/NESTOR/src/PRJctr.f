C   +-------------------------------------------------------------------+
C   |  Subroutine lamphi2laea                       April 2001  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Convertion of lat./long. coordinates into Lambert coordinates     |
C   | => Lambert-Azimuthal Equal-Area Projection on Spheroid            |
C   |                                                                   |
C   | Input : - lam, phi : coordinates in latitude and longitude        |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: - xl, yl : Lambert coordinates                            |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE lamphi2laea (xl,yl,lam,phi,lam0,phi0,r)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      REAL  xl,yl,lam,phi,lam0,phi0,r,kprim


C +---Convert to Lambert coordinates
C +   ------------------------------

      kprim = SQRT(2./(1+SIN(phi0)*SIN(phi)
     .                  +COS(phi0)*COS(phi)*COS(lam-lam0)))

      xl = r * kprim *  COS(phi )*SIN(lam-lam0)
      yl = r * kprim * (COS(phi0)*SIN(phi)
     .                 -SIN(phi0)*COS(phi)*COS(lam-lam0))


      RETURN
      END


C   +-------------------------------------------------------------------+
C   |  Subroutine laea2lamphi                       April 2001  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Convertion of Lambert coordinates into lat./long. coordinates     |
C   | => Lambert-Azimuthal Equal-Area Projection on Spheroid            |
C   |                                                                   |
C   | Input : - xl, yl : Lambert coordinates                            |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: - lam, phi : coordinates in latitude and longitude        |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE laea2lamphi (xl,yl,lam,phi,lam0,phi0,r)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      REAL  xl,yl,lam,phi,lam0,phi0,r,EPS10,rho,c

      DATA  EPS10  / 1.E-10 /


C +---Convert Lambert coordinates
C +   ---------------------------

      rho = SQRT(xl*xl+yl*yl)
      c = 2. * ASIN(rho/(2.*r))

      IF (rho.lt.EPS10) THEN
       phi = 0.
       lam = 0.
      ELSE
       phi = ASIN(COS(c)*SIN(phi0)+(yl*SIN(c)*COS(phi0)/rho))
       lam = lam0
     .     + ATAN(xl*SIN(c)/(rho*COS(phi0)*COS(c)-yl*SIN(phi0)*SIN(c)))
      ENDIF


      RETURN
      END


C   +-------------------------------------------------------------------+
C   |  Subroutine lambel2geo                        April 2001  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Convertion of Lambert coordinates into lat./long. coordinates     |
C   | => Lambert-Azimuthal Equal-Area Projection on Spheroid            |
C   |                                                                   |
C   | Input : - lam, phi : coordinates in latitude and longitude        |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: - xl, yl : Lambert coordinates                            |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE lambel2geo (xl,yl,lam,phi)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INTEGER i

      REAL  xl,yl,lam,phi,a,p,phi1,phi2,pi,lambbr,rad,e,gam1,gam2,
     .      e2,w1,w2,z1,z2,n,n1,k,lamdif,gam,teta,rho,dx,dy,dz,
     .      fact


C +---Initialization
C +   --------------

      a      = 6378388.
      p      = 297.
      phi1   = 2990.
      phi2   = 3070.
      pi     = 3.14159265358979
      lambbr = 4.3680
      rad    = pi / 180.

      p      = 1. / p
      e      = SQRT(2. * p - p * p)
      phi1   = phi1 * rad/60.
      phi2   = phi2 * rad/60.
      gam1   = pi/2. - phi1
      gam2   = pi/2. - phi2
      e2     = e/2.
      w1     = SQRT(1.-(e * COS(gam1))**2.)
      w2     = SQRT(1.-(e * COS(gam2))**2.)
      z1     = TAN(gam1/2.) * (((1.+e*COS(gam1)) 
     .                         /(1.-e*COS(gam1)))**e2)
      z2     = TAN(gam2/2.) * (((1.+e*COS(gam2))
     .                         /(1.-e*COS(gam2)))**e2)
      n      = (ALOG10(w2)-ALOG10(w1)+ALOG10(SIN(gam1)) 
     .       - ALOG10(SIN(gam2))) / (ALOG10(z1)-ALOG10(z2))
      n1     = 1./n

      k      = a*SIN(gam1)/(n*w1*(z1**n))

      dx     = xl - 150000.
      dy     = 5400000. - yl
      rho    = SQRT(dx*dx+dy*dy)
      teta   = ATAN(dx/dy)
      lamdif = teta/(n*rad)
      dz     = (rho/k)**n1

      gam    = dz
      DO i=1,3
       fact  = ((1.+e+(1.-e)*gam*gam) / (1.-e+(1.+e)*gam*gam))**e2
       gam   = dz/fact
      ENDDO
      gam    = 2.*ATAN(gam)

      phi    = (pi/2.-gam) / rad
      lam    = lambbr + lamdif


      RETURN
      END


C   +-------------------------------------------------------------------+
C   |  Subroutine geo2lambel                        April 2001  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Convertion of Lambert coordinates into lat./long. coordinates     |
C   | => Lambert-Azimuthal Equal-Area Projection on Spheroid            |
C   |                                                                   |
C   | Input : - xl, yl : Lambert coordinates                            |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: - lam, phi : coordinates in latitude and longitude        |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE geo2lambel (xl,yl,lam,phi)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      REAL  xl,yl,lam,phi,a,p,phi1,phi2,pi,lambbr,rad,e,gam1,gam2,
     .      e2,w1,w2,z1,z2,n,k,expn,lamdif,gam,teta,rho


C +---Initialization
C +   --------------

      a      = 6378388.
      p      = 297.
      phi1   = 2990.
      phi2   = 3070.
      pi     = 3.14159265358979
      lambbr = 4.3680
      rad    = pi / 180.


C +---Convert Lambert coordinates
C +   ---------------------------

      p    = 1./p
      e    = SQRT(2.*p-p*p)
      phi1 = phi1 * rad / 60.
      phi2 = phi2 * rad / 60.
      gam1 = pi/2. - phi1
      gam2 = pi/2. - phi2
      e2   = e/2.
      w1   = SQRT(1.-(e*COS(gam1))**2.)
      w2   = SQRT(1.-(e*COS(gam2))**2.)
      z1   = TAN(gam1 / 2.) 
     .     * (((1.+e*COS(gam1)) / (1.-e*COS(gam1)))**e2)
      z2   = TAN(gam2 / 2.) 
     .     * (((1.+e*COS(gam2)) / (1.-e*COS(gam2)))**e2)
      n    = (ALOG10(w2) - ALOG10(w1) + ALOG10(SIN(gam1)) 
     .     -  ALOG10(SIN(gam2))) / (ALOG10(z1) - ALOG10(z2))

      k    = a * SIN(gam1) / (n * w1 * (z1**n))
      expn = e2*n

      lamdif = lam - lambbr
      gam    = pi / 2. - phi * rad
      teta   = n * lamdif * rad
      rho    = k * (TAN(gam/2.)**n) 
     .       * (((1.+e*COS(gam)) / (1.-e*COS(gam)))**expn)

      xl =  150000. + rho * SIN(teta)
      yl = 5400000. - rho * COS(teta)


      RETURN
      END

