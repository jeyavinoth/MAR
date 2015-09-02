
C +  ==========
C +--MAR_LB.inc
C +  ==========

C +...n6 et n7 determine a relaxation zone towards lateral boundaries 
C +   (large scale values of the variables). 
C +   This zone extends over n6-1 points. 
C +   Davies (1976) propose 5 points (i.e. n6=6 and n7=7)

      integer       iyr_LB,mma_LB,jda_LB,jhu_LB,jdh_LB
      common/nudite/iyr_LB,mma_LB,jda_LB,jhu_LB,jdh_LB
C +...              iyr_LB: Year
C +                 mma_LB: Month
C +                 jda_LB: Day
C +                 jhu_LB: Hour  (UT)
C +                 jdh_LB: Time Interval before next GCM/NWP LBC (hour)
C +                 jdh_LB=0 ==> NO further GCM/NWP LBC available
      integer*8     tim1LB,tim2LB
      common/nudtim/tim1LB,tim2LB
C +...              tim1LB: Time of the previous    LBC         (second)
C +                 tim2LB: Time of the next        LBC         (second)

      integer       n40xLB, n50xLB, n5mxLB, n6mxLB, n7mxLB, 
     .              n40yLB, n50yLB, n5myLB, n6myLB, n7myLB 
      common/nudind/n40xLB, n50xLB, n5mxLB, n6mxLB, n7mxLB, 
     .              n40yLB, n50yLB, n5myLB, n6myLB, n7myLB 
C +...           ...n6mxLB, n7mxLB, n6myLB, n7myLB... 
C +                           define the effective length of the lateral sponge

      REAL          vaxgLB              ,vaxdLB                   
      REAL          v1xgLB              ,v1xdLB                   
     .             ,v2xgLB              ,v2xdLB                   
      real          tixgLB              ,tixdLB                   
      common/nuddax/vaxgLB(1:n7,my,mz,5),vaxdLB(mx-n6:mx ,my,mz,5)
     .             ,v1xgLB(1:n7,my,mz,5),v1xdLB(mx-n6:mx ,my,mz,5)
     .             ,v2xgLB(1:n7,my,mz,5),v2xdLB(mx-n6:mx ,my,mz,5)
     .             ,tixgLB(2:n7,my,mz  ),tixdLB(mx-n6:mx1,my,mz  )

      REAL          vayiLB              ,vaysLB                   
      REAL          v1yiLB              ,v1ysLB                   
     .             ,v2yiLB              ,v2ysLB                   
      real          tiyiLB              ,tiysLB                   
      common/nudday/vayiLB(mx,1:n7,mz,5),vaysLB(mx,my-n6:my ,mz,5)
     .             ,v1yiLB(mx,1:n7,mz,5),v1ysLB(mx,my-n6:my ,mz,5)
     .             ,v2yiLB(mx,1:n7,mz,5),v2ysLB(mx,my-n6:my ,mz,5)
     .             ,tiyiLB(mx,2:n7,mz  ),tiysLB(mx,my-n6:my1,mz  )
C +...              vaXX : large scale values of relevant dependant variables 
C +                   ^X=(x->x axis border, y->y axis border)
C +                    ^X=(g->x small, d->x large, b->y small, h->y large)
C +                 tiXXLB : independant term of semi-implicit numerical scheme

      real          wixgLB 
     .             ,wixdLB
     .             ,wiyiLB
     .             ,wiysLB
      common/nuddaw/wixgLB(   2:  n7,   2:  n7)
     .             ,wixdLB(mx-n6:mx1,mx-n6:mx1)
     .             ,wiyiLB(   2:  n7,   2:  n7)
     .             ,wiysLB(my-n6:my1,my-n6:my1)
C +...              wiXXLB : coefficient used in semi-implicit numerical scheme

      real          rxLB    ,ryLB    
      common/nuddtk/rxLB(mx),ryLB(my)
C +...              rXLB   : nudging coefficients of the relaxation zone

      REAL          sst_LB       ,Tfr_LB
      REAL                        sst1LB       ,sst2LB       
      common/srfbnd/sst_LB(mx,my),sst1LB(mx,my),sst2LB(mx,my),Tfr_LB
C +...              sst_LB : external SST
C +                 Tfr_LB : Freezing Temperature 
C +                         (corrected when using Reynolds data set)

