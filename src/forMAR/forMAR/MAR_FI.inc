      real            FIsloq
     .              , FIslot,FIk_st(mz),FIfstt,FIk_ft(mz)
     .              , FIslou,FIk_su(mz),FIfstu,FIk_fu(mz)
     .              , FIslop,           FIfstp,FIk_fp(mz)
     .              , FIkhmn,FIspon(mzabso)
      common /parfil/ FIsloq
     .              , FIslot,FIk_st    ,FIfstt,FIk_ft
     .              , FIslou,FIk_su    ,FIfstu,FIk_fu
     .              , FIslop,           FIfstp,FIk_fp
     .              , FIkhmn,FIspon
C +...                FIsloq           : Implicit Filter Parameter (Slow mPHYS)
C +...                FIslot,FIk_st(mz): Implicit Filter Parameter (Slow Dyn./Temperature)
C +...                FIslou,FIk_su(mz): Implicit Filter Parameter (Slow Dyn./Wind Speed)
C +                   FIfstu,FIk_fu(mz): Implicit Filter Parameter (Fast Dyn./Wind Speed)
C +                   FIslop:            Implicit Filter Parameter (Slow Dyn./Pressure)  
C +                   FIfstp,FIk_fp(mz): Implicit Filter Parameter (Fast Dyn./Pressure)  
C +...                FIkhmn:            Horizontal Diffusion Coefficient
C +                                      equivalent to the Filter Effect on Long Waves
C +                   FIspon:            Implicit Filter Parameter (Top Absorber)
      logical         FIBord
      parameter      (FIBord=.TRUE.)
      real            FacFIk
      parameter      (FacFIk=10.00)                                       !+CA+! (10.00>20.00)
