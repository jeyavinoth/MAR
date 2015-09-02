C +...*Physical constants (LSMARIN/LSMAROUT)
C +...
      REAL cp, cap, ra, grav, gamTz
      REAL pi, epsi, earthr, degrad, hourad

      PARAMETER (pi= 3.14159265358979d0)
C +...    The well known value of acos(-1) ;-). 

      PARAMETER (degrad= pi / 180.d0)
      PARAMETER (hourad= pi /  12.d0)
C +...    From deg to rad, and hour to rad.

      PARAMETER (epsi= 1.0d-6)          
C +...    The definition of "little" by Hub.

      PARAMETER (earthr= 6371.229d+3)
C +...    The radius of the little blue ball.

      PARAMETER (ra= 287.d0)
C +...    Perfect Gas Law  Constant (J/kg/K)

      PARAMETER (cp= 1004.d0)
C +...    Air Specific Heat         (J/kg/K)

      PARAMETER (cap= 0.28586d0)
C +...

      PARAMETER (grav= 9.81   d0)
C +...    The way we are falling    (m/s2)

      PARAMETER (gamTz = - 6.5E-3)
C +...    Mean temperature gradient (K/m)

