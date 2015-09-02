
C +  ==========
C +--MAR_BS.inc
C +  ==========

       integer               nn_pro
       parameter            (nn_pro=30)

       real           BS_pro(nn_pro)
       common/ProfBS/ BS_pro
C +...  BS_pro : Normalized Blowing Snow Profile   (Suspension Layer)

       real           z0SaBS(mx,my,mw)
     .               ,z0emBS(mx,my,mw)
       common/SastBS/ z0SaBS
     .               ,z0emBS
C +...  z0SaBS : Z0 due to Sastrugi Height                                  [m]
C +     z0emBS : Z0 due to Snow     Erosion                                 [m]

       real           hSalBS(mx,my),u_stBS(mx,my),wSalBS(mx,my)
       common/SaltBS/ hSalBS       ,u_stBS       ,wSalBS
C +...  hSalBS : Height  above the Surface  (Top of Saltation  Layer)       [m]
C +     u_stBS : Time Averaged Friction Velocity for 
C +              Blowing Snow Surface Turbulent Flux                      [m/s]
C +     wSalBS : Blowing Snow Sedimentation (Top of Saltation  Layer)     [m/s]

       real           FacRBS,FacSBS,FacTBS,FacUBS
       real           g1__BS,g2__BS,sheaBS
       common/StreBS/ FacRBS,FacSBS,FacTBS,FacUBS,
     .                g1__BS,g2__BS,sheaBS

       real           FracBS,ua_0BS(mx,my),va_0BS(mx,my)
     .                      ,VVs_BS(mx,my),DDs_BS(mx,my),RRs_BS(mx,my)
       common/SastBS/ FracBS,ua_0BS       ,va_0BS       
     .                      ,VVs_BS       ,DDs_BS       ,RRs_BS
