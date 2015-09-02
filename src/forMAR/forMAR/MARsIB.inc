C +------------------------------------------------------------------------+
C | MARsIB.inc:  Include    Sea-Ice                        19-02-2004  MAR |
C | ==========                                                             |
C |                                                                        |
C |   Commons    MARsIB are used to follow    MAR Sea-Ice    Fraction      |
C |                                                                        |
C +------------------------------------------------------------------------+

      REAL    sicsIB(mx,my)
C +...        sicsIB:     Sea-Ice Fraction          / Time Interpolated

      REAL    sic1sI(mx,my)
C +...        sic1sI:     Sea-Ice Fraction          / Time Step n

      REAL    sic2sI(mx,my)
C +...        sic2sI:     Sea-Ice Fraction          / Time Step n+1

      common /rsvasI/sic1sI,sic2sI,sicsIB

      INTEGER         iyr_sI,mma_sI,jda_sI,jhu_sI,jdh_sI
      common /isvasI1/iyr_sI,mma_sI,jda_sI,jhu_sI,jdh_sI
      INTEGER*8       tim1sI,tim2sI
      common /isvasI2/tim1sI,tim2sI
