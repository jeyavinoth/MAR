C +   * Number of points of the numerical filter for Z500:
      INTEGER n2fl, nfl, nfl1
      PARAMETER (n2fl=15, nfl=(n2fl*2+1), nfl1=(nfl+1))

C +   Filco = filter coefficients

      REAL*8 FilCo(nfl)
      
      common /tSDBP/ FilCo
