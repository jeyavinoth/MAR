C +..BEGIN    "MARdSV.inc"
C +--SISVAT Global Variables
C +  =======================

      logical       INI_SV               ! Initialisation Switch
      common/SwchSV/INI_SV               ! 
      real          eps_21               ! Arbitrary  very small value
      common/GendSV/eps_21               !


C +--Snow
C +  ----

      integer       istdSV(5)            ! Snow History
      common/SnidSV/istdSV               !

      real          Cn_dSV               ! Snow Heat Capacity          [J/kg/K]
      real          SMndSV               ! Minimum Thickness of new Layers
      real          G1_dSV               ! Conversion 0/99-->0/1 
      real          DDcdSV,DFcdSV,DScdSV ! Snow Grains Optical Diameter [1e-4m]
      real                        ADSdSV ! Snow Grains Actual  Diameter [1e-4m]
      real          So1dSV,So2dSV,So3dSV ! Total Solar Irradiance Fractions [-]
      real          aI1dSV,aI2dSV,aI3dSV ! Bare Ice Albedo                  [-]
      real          ws0dSV               ! Irreducible Water Saturation in Snow
      real          roCdSV               ! Pore Hole Close OFF Density  [kg/m3]
      real          ru_dSV               ! Surficial Water Scale Factor [kg/m2]

      common/Sn_dSV/Cn_dSV,
     .              SMndSV,
     .              G1_dSV,
     .              DDcdSV,DFcdSV,DScdSV,
     .                            ADSdSV,
     .              So1dSV,So2dSV,So3dSV,
     .              aI1dSV,aI2dSV,aI3dSV,
     .              ws0dSV,
     .              roCdSV,
     .              ru_dSV


C +--Ice
C +  ---

      real          CdidSV
      common/IcedSV/CdidSV


C +--Vegetation
C +  ----------

      integer                nvgt
      parameter             (nvgt=12)
      real          DH_dSV(0:nvgt)       ! Displacement            Height   [m]
      real          Z0mdSV(0:nvgt)       ! Roughness  Length for Momentum   [m]
      real          StodSV(0:nvgt)       ! Minimum    Stomatal Resistance [s/m]
      real          rbtdSV(0:nvgt)       ! Root Fraction Beta Coefficient   [-]
      real          PR_dSV(0:nvgt)       ! Internal Plant      Resistance   [s]
      real          pscdSV               ! Critical Leaf Water Potential    [m]
      real          StxdSV               ! maXimum  Stomatal   Resistance [s/m]
      real          LAIdSV               ! maximum  LAI
      common/Vg_dSV/DH_dSV,Z0mdSV,StodSV,rbtdSV,PR_dSV
     .             ,pscdSV,StxdSV,LAIdSV


C +--Soil
C +  ----

      real          rcwdSV               ! Density * Water Specific Heat
      real          dz_dSV(-nsol:0)      ! Vertical  Discretization
      real          zz_dSV               ! Soil      Thickness
      common/SozdSV/rcwdSV,dz_dSV,zz_dSV !

      integer                nsot
      parameter             (nsot=12)
      real          etadSV(0:nsot)       ! Water Content at Saturation  [kg/kg]
      real          psidSV(0:nsot)       ! Water Succion at Saturation      [m]
      real          Ks_dSV(0:nsot)       ! Hydraulic Conductivity
                                         !               at Saturation    [m/s]
      real          bCHdSV(0:nsot)       ! Clapp-Hornberger Coefficient b   [-]
      common/So_dSV/etadSV,psidSV,Ks_dSV,bCHdSV


C +--Water Bodies
C +  ------------

      real          vK_dSV               ! Diffusivity in Water          [m2/s]
      real          TSIdSV               ! Sea-Ice Fraction: SST Scale      [K]
      common/Wa_dSV/vK_dSV,TSIdSV        !


C +..END      "MARdSV.inc"
