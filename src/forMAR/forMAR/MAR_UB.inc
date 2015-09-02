C +----------------------------------------------------------------------------+
C | MAR_UB.inc:  Include    Upper Sponge                        2-03-2004  MAR |
C | ==========                                                                 |
C |                                                                            |
C |   Commons    MAR_UB are used to follow Reference State in MAR upper sponge |
C |   Reference: ARPS 4.0 User's Guide, para 6.4.3 p.152                       |
C |                                                                            |
C +----------------------------------------------------------------------------+


      REAL          Ray_UB(      mzabso)

      common/rCstUB/Ray_UB

      REAL          uairUB(mx,my,mzabso),vairUB(mx,my,mzabso)
      REAL          pktaUB(mx,my,mzabso)

      common/rVarUB/uairUB              ,vairUB
     .             ,pktaUB

      REAL          ua1_UB(mx,my,mzabso),va1_UB(mx,my,mzabso)
      REAL          pkt1UB(mx,my,mzabso)
      REAL          ua2_UB(mx,my,mzabso),va2_UB(mx,my,mzabso)
      REAL          pkt2UB(mx,my,mzabso)

      common/rNstUB/ua1_UB              ,va1_UB
     .             ,pkt1UB
     .             ,ua2_UB              ,va2_UB
     .             ,pkt2UB


      INTEGER        iyr_UB,mma_UB,jda_UB,jhu_UB,jdh_UB
      common/iVarUB1/iyr_UB,mma_UB,jda_UB,jhu_UB,jdh_UB
      INTEGER*8      tim1UB,tim2UB
      common/iVarUB2/tim1UB,tim2UB     

