C +------------------------------------------------------------------------+
C | MAR Include    SVAT                                    13-05-2002  MAR |
C |   Commons    MAR_VB are used to follow    MAR Vegetation Fraction      |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
      REAL    glf1VB(mx,my,nvx)
C +...        glf1VB:     Vegetation Class Coverage / Time Step n
C +                               (3 Class, Last One is Open Water)
C +
      REAL    LAI1VB(mx,my,nvx)
C +...        LAI1VB:     Leaf Area Index           / Time Step n
C +
C +
      REAL    glf2VB(mx,my,nvx)
C +...        glf2VB:     Vegetation Class Coverage / Time Step n+1
C +                               (3 Class, Last One is Open Water)
C +
      REAL    LAI2VB(mx,my,nvx)
C +...        LAI2VB:     Leaf Area Index           / Time Step n+1
C +
      common /rsvaVB/glf1VB,glf2VB,LAI1VB,LAI2VB
C +
      INTEGER         iyr_VB,mma_VB,jda_VB,jhu_VB,jdh_VB
      common /isvaVB1/iyr_VB,mma_VB,jda_VB,jhu_VB,jdh_VB
      INTEGER*8       tim1VB,tim2VB
      common /isvaVB2/tim1VB,tim2VB
