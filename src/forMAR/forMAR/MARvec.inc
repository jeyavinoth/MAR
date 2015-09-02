      real            vecxa(mx)
      real            vecx1(mx),vecx2(mx),vecx3(mx)
      real            vecx4(mx),vecx5(mx),vecx6(mx)
      real            delta
      common/vectorX/ vecxa
     .               ,vecx1    ,vecx2    ,vecx3    
     .               ,vecx4    ,vecx5    ,vecx6
     .               ,delta
C +
      real            vecy1(my),vecy2(my),vecy3(my)
      common/vectorY/ 
     .                vecy1    ,vecy2    ,vecy3    
C +
      real            dumeps(      mz)
      real            dumy3D(mx,my,mz)
      common/dummy_R/ dumeps,dumy3D
