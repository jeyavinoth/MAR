
C +  ==========
C +--MAR_WK.inc
C +  ==========

      real          WKxyz1(mx,my,mz ),WKxyz2(mx,my,mz ),
     &              WKxyz3(mx,my,mz ),WKxyz4(mx,my,mz ),
     &              WKxyz5(mx,my,mzz),WKxyz6(mx,my,mzz),
     &              WKxyz7(mx,my,mzz),WKxyz8(mx,my,mzz)

      real          WTxyz1(mx,my,mz ),WTxyz2(mx,my,mz ),
     &              WTxyz3(mx,my,mz ),WTxyz4(mx,my,mz ),
     &              WTxyz5(mx,my,mzz),WTxyz6(mx,my,mzz),
     &              WTxyz7(mx,my,mzz),WTxyz8(mx,my,mzz)

      common/wkxyz/ WKxyz1           ,WKxyz2           ,
     &              WKxyz3           ,WKxyz4           ,
     &              WKxyz5           ,WKxyz6           ,
     &              WKxyz7           ,WKxyz8           
      real          WRKxys(mx,my,mw )
      common/wkxys/ WRKxys
      real          WKxy1 (mx,my)    ,WKxy2 (mx,my)    ,WKxy3 (mx,my)
     &            , WKxy4 (mx,my)    ,WKxy5 (mx,my)    ,WKxy6 (mx,my)
     &            , WKxy7 (mx,my)    ,WKxy8 (mx,my)    ,WKxy9 (mx,my)
     &            , WKxy0 (mx,my)
      real          WTxy1 (mx,my)    ,WTxy2 (mx,my)    ,WTxy3 (mx,my)
     &            , WTxy4 (mx,my)    ,WTxy5 (mx,my)    ,WTxy6 (mx,my)
     &            , WTxy7 (mx,my)    ,WTxy8 (mx,my)    ,WTxy9 (mx,my)
     &            , WTxy0 (mx,my)
      common/wkxy/  WKxy1            ,WKxy2            ,WKxy3
     &            , WKxy4            ,WKxy5            ,WKxy6
     &            , WKxy7            ,WKxy8            ,WKxy9
     &            , WKxy0
      real          WKxza (mx,mz)    ,WKxzb (mx,mz)    ,WKxzc (mx,mz)
     &            , WKxzd (mx,mz)    ,WKxzx (mx,mz)
     &            , WKxzp (mx,mz)    ,WKxzq (mx,mz)
      common/wkxz/  WKxza            ,WKxzb            ,WKxzc
     &            , WKxzd            ,WKxzx
     &            , WKxzp            ,WKxzq
C +...wkXXX variables define a work area 
C +         in order to minimize memory requirements
