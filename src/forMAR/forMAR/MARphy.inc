
C +  ==========
C +--MARphy.inc
C +  ==========

      character*1     labnum
      common /labgen/ labnum(10)
C +...                labnum: Alphanumeric Character

      integer         izr,iun
      common /numint/ izr,iun
C +...                izr = 0     
C +                       iun = 1     
      real                zero,demi,unun,epsi,eps9,eps12,third,thous
      real                argmin,argmax
      common /numrea/     zero,demi,unun,epsi,eps9,eps12,third,thous
C +...                    zero = 0.0e+0
C +                            demi = 0.5e+0
C +                                 unun = 1.0e+0
C +                                      epsi = 1.0e-6
C +                                           eps9 = 1.0e-9
C +                                                eps12= 1.0e-12
C +                                                      third=1/3
C +                                                            thous=1.e3
     .                   ,argmin,argmax
C +                       argmin:        Function exp(x) Minimum Argument
C +                              argmax: Function exp(x) Maximum Argument
C +                                      (Min/Max are Machine Dependant)

      real            pi,degrad,hourad
      common /phyphy/ pi,degrad,hourad
C +...                pi    = 3.141592653589793238462643e0
C +                   degrad= pi / 180.d0
C +                   hourad= pi /  12.d0

      real            earthr,earthv,gravit,gravi2,grvinv
      common /phygeo/ earthr,earthv,gravit,gravi2,grvinv
C +...                earthr: Earth Radius                  =6371.d+3    m     
C +...                earthv: Earth Angular Velocity        = 7.29d-5  s-1
C +                   gravit: Earth Gravity Acceleration    = 9.81    m/s2
C +                   gravi2: idem (squared) 
C +                   grvinv: idem (inverse) 

      real            akmol, vonkar, A_Turb, AhTurb, AsTurb, r_Turb
      common /phytur/ akmol, vonkar, A_Turb, AhTurb, AsTurb, r_Turb
C +...                akmol : Air Viscosity                 = 1.35d-5 m2/s
C +                   vonkar: von Karman constant           = 0.4
C +                   A_Turb: Stability  Coefficient Moment = 5.8
C +                   AhTurb: Stability  Coefficient Heat   = 5.4
C +                   AsTurb: Stability  Coefficient BLOW * = 4.0
C +                   r_Turb: Turbulent Diffusivities Ratio K*/Km

      real            RDryAi,ra,cp,CvDrya,racv,cap,pcap,stefan,qv_MIN
      real            RVapor,epsq,Lv_H2O,Ls_H2O,r_LvCp,r_LcCp,r_LsCp
      common /phyatm/ RDryAi,ra,cp,CvDrya,racv,cap,pcap,stefan,qv_MIN
     .               ,RVapor,epsq,Lv_H2O,Ls_H2O,r_LvCp,r_LcCp,r_LsCp
C +...                RDryAi: perfect gas law constant for dry air (287 J/kg/K)
C +                   ra    : perfect gas law constant for dry air (287 J/kg/K)
C +                   cp    : dry air specific heat at constant p (1004 J/kg/K)
C +                   CvDryA: dry air specific heat at constant V ( 717 J/kg/K)
C +                   racv  = ra / cv      = ra /(cp-ra)
C +                   cap   = ra / cp
C +                   pcap  = 100.[kPa]**cap 
C +...                RVapor: perfect gas law constant Water Vapor (461 J/kg/K)
C +                   epsq  : Minimum Water Vapor Content    =    3.00d-6 kg/kg  
C +                   Lv_H2O: Latent Heat Vaporisation/Water = 2500.00d+3 J/kg
C +                   Ls_H2O: Latent Heat Sublimation /Ice   = 2833.60d+3 J/kg
C +                   r_LvCp [Lv=2500000J/kg]/[Cp=1004J/kg/K]= 2490.04    K/[kg/kg]
C +                   r_LvCp [Lc= 333600J/kg]/[Cp=1004J/kg/K]=  332.27    K/[kg/kg]
C +                   r_LvCp [Ls=2833600J/kg]/[Cp=1004J/kg/K]= 2822.31    K/[kg/kg]
C +                   stefan: Stefan Constant                =    5.67d-8 W/m2/K4
C +                   qv_MIN: Minimum Specific Humidity      =    3.00d-6 kg/kg

      real            tfrwat,siclf,cdsice,hic0,ro_Wat,C__Wat,fracoh
      common /physea/ tfrwat,siclf,cdsice,hic0,ro_Wat,C__Wat,fracoh
C +...                tfrwat: Sea-Water freezing Temperature =  271.20d+0 K
C +                   siclf : Sea-Ice Heat of Fusion         =  302.00d+6 J/m3
C +                   cdsice: Sea-Ice Thermal Conductivity   =    2.04d+0 W/mK
C +                   hic0  : Sea-Ice Initial Thickness      =    2.00d+0 m
C +                   ro_Wat: Density       of Water         = 1000.00d+0 kg/m3
C +                   C__Wat: Heat Capacity of Water         = 4186.00d+0 J/kg/K
C +                   fracoh: 25% of cooling = Oceanic Heat Flux 
C +***                        (Hibler 1984)    (ANTARCTIC Ocean)

      real            ro_Ice,cdice
      common /phyice/ ro_Ice,cdice
C +...                ro_Ice:  Density      of Pure Ice       =  920.00d+0 kg/m3
C +                   cdice :  Conductivity of Pure Ice       =    2.51d+0 W/m/K

      real            TfSnow,csnow,r0sno,blsno,Lf_H2O
      common /physno/ TfSnow,csnow,r0sno,blsno,Lf_H2O
C +...                TfSnow:        Snow melting Temperature=  273.15d+0 K
C +                    csnow:Heat Capacity of Snow             2105      J/kg/K
C +                         (Loth et al. 1993, JGR 98 D6, 2.2.2 2e para p.10453)
C +                   r0sno : Fresh  Snow Density                50.00d+0 kg/m3
C +                   blsno : Blowed Snow Density               330.00d+0 kg/m3
C +                   Lf_H2O: Latent Heat of Fusion of Snow  =    3.34d+5 J/kg
