C +...BEGIN of MAR_WB.inc
      real            dq__WB(mx,my,0:9),timeWB(0:9)
      real            qLB_WB(      0:9),qSB_WB(0:9)
      real            qt__WB(      0:9)
      real            wq__WB(mx,my)    ,uq__WB(mx,my),us__WB(mx,my)
      real            rr_0WB(mx,my)    ,rrc0WB(mx,my)
      real            ss_0WB(mx,my)    ,ssc0WB(mx,my)
      real            drr_WB(mx,my)    ,drrcWB(mx,my)
      real            dss_WB(mx,my)    ,dsscWB(mx,my)
      real            upq_WB(mx,my)    ,vpq_WB(mx,my),cpq_WB(mx,my)
      real            ups_WB(mx,my)    ,vps_WB(mx,my),cps_WB(mx,my)
      real            qinpWB           ,qoutWB
      real            qOv_WB           ,qLv_WB
      real            qLrsWB           ,qOrsWB
      real            qL_0WB           ,qO_0WB
      real            q0__WB           ,qn__WB
      real            dqt_WB           ,timtWB
      common/WaterB_r/dq__WB           ,timeWB
     .               ,qLB_WB           ,qSB_WB
     .               ,qt__WB
     .               ,wq__WB           ,uq__WB       ,us__WB
     .               ,rr_0WB           ,rrc0WB
     .               ,ss_0WB           ,ssc0WB
     .               ,drr_WB           ,drrcWB
     .               ,dss_WB           ,dsscWB
     .               ,upq_WB           ,vpq_WB       ,cpq_WB
     .               ,ups_WB           ,vps_WB       ,cps_WB
     .               ,qinpWB           ,qoutWB
     .               ,qOv_WB           ,qLv_WB
     .               ,qLrsWB           ,qOrsWB
     .               ,qL_0WB           ,qO_0WB
     .               ,q0__WB           ,qn__WB
     .               ,dqt_WB           ,timtWB
C +...END   of MAR_WB.inc
