
C +  ==========
C +--MAR_DY.inc
C +  ==========

      logical       openmp
      common/dynlog/openmp      

      integer       iyrDYN,mmaDYN,jdaDYN,jhuDYN
      common/dyntim/iyrDYN,mmaDYN,jdaDYN,jhuDYN
C +...              iyrDYN: Year
C +                 mmaDYN: Month
C +                 jdaDYN: Day
C +                 jhuDYN: Hour  (UT)


      real           uairDY(mx,my,mz), vairDY(mx,my,mz)
     .             , pktaDY(mx,my,mzz),pkt_DY(mx,my,mz),tairDY(mx,my,mz)
     .             , pkDY  (mx,my,mz), qvsiDY(mx,my,mzz)
     .             , qvDY  (mx,my,mz), qvswDY(mx,my,mzz)
     .             , pstDY1(mx,my)    
     .             ,opstDY (mx,my)    
     .             , pstDY (mx,my)    
     .             , pstDYn(mx,my) 
     .             , pstDY2(mx,my) 
     .             , ptopDY       
     .             , fcorDY(mx,my)       
     .             , SFm_DY(mx,my)       
      common/dynatm/ uairDY          , vairDY           
     .             , pktaDY          , pkt_DY          ,tairDY           
     .             , pkDY            , qvsiDY           
     .             , qvDY            , qvswDY           
     .             , pstDY1           
     .             ,opstDY            
     .             , pstDY            
     .             , pstDYn        
     .             , pstDY2        
     .             , ptopDY       
     .             , fcorDY              
     .             , SFm_DY              
C +...       uairDY: x-wind speed component                       (m/s)
C +          vairDY: y-wind speed component                       (m/s)
C +...       pktaDY: potential temperature divided by 100.[kPa]**(R/Cp)
C +          pkt_DY: potential temperature divided by 100.[kPa]**(R/Cp)
C +          tairDY: real      temperature                          (K)
C +...       pkDY  : Exner Potential                   = p[kPa]**(R/Cp)
C +...       qvDY  :            Specific Humidity               (kg/kg)
C +          qvsiDY: Saturation Specific Humidity over Ice      (kg/kg)
C +          qvswDY: Saturation Specific Humidity over Water    (kg/kg)
C +          pstDY1: Model Pressure Depth at INITIAL  Time Step   (kPa)
C +         opstDY : Model Pressure Depth at previous Time Step   (kPa)
C +          pstDY : Model Pressure Depth at current  Time Step   (kPa)
C +          pstDYn: Model Pressure Depth at next     Time Step   (kPa)
C +          pstDY2: Model Pressure Depth   (current) (squared)  (kPa2)
C +          ptopDY:       Pressure       at Model Top            (kPa)
C +          fcorDY: Coriolis Parameter                           (s-1)
C +          SFm_DY: Map Scale Factor                               (-)


      real          rolvDY(mx,my,mz),romiDY(mx,my,mz),virDY(mx,my,mz)
      common/rhodyn/rolvDY          ,romiDY          ,virDY          
C +...rolvDY: air specific mass index (i,j,k) at (i,j,k)       (Ton/m3)
C +   romiDY: air specific mass index (i,j,k) at (i,j,k+1/2)   (Ton/m3)
C +    virDY: air loading by water (vapor, liquid and solid)


      real          ubefDY(mx,my,mz),vbefDY(mx,my,mz),
     .              dgzxDY(mx,my,mz),dg1xDY(mx,my,mz),
     .              dgzyDY(mx,my,mz),dg1yDY(mx,my,mz)
      common/dynnum/ubefDY          ,vbefDY          ,
     .              dgzxDY          ,dg1xDY          ,
     .              dgzyDY          ,dg1yDY          
C +... ubefDY: uairDY at Previous Time Step (in the Leap-Frog)
C +    vbefDY: vairDY at Previous Time Step (in the Leap-Frog)
C +    dgzxDY is the geopotential gradient just before current  time step
C +    dg1xDY is the geopotential gradient just before previous time step
C +       ^X=(x->x direction, y->y direction)


      real          wairDY(mx,my,mz) 
     .             ,psigDY(mx,my,mz)
     .             ,wsigDY(mx,my,mz)
     .             ,CFLzDY
      common/dynver/wairDY
     .             ,psigDY          
     .             ,wsigDY          
     .             ,CFLzDY
C +... wairDY :      Vertical Wind Speed (in z     coordinate)    (cm/s)
C +    psigDY : p* X Vertical Wind Speed (in sigma coordinate)
C +    wsigDY :      Vertical Wind Speed (in sigma coordinate)
C +    CFLzDY :      Vertical CFL        (in sigma coordinate)


      real          gplvDY(mx,my,mzz),
     .              gpmiDY(mx,my,mzz)
      common/dyngeo/gplvDY           ,
     .              gpmiDY           
C +... gplvDY: Geopotential of      sigma Levels (= g z)
C +    gpmiDY: Geopotential between sigma Levels (i,j,k)->layer(i,j,k-1/2)


      real          ugeoDY(mx,my,mz),vgeoDY(mx,my,mz)
      common/dynlsc/ugeoDY          ,vgeoDY          
C +... ugeoDY: x-geostrophic wind speed (assumed external)             (m/s)
C +    vgeoDY: y-geostrophic wind speed (assumed external)             (m/s)
