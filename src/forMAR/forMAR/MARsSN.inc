
C +  ==========
C +--MARsSN.inc
C +  ==========

      integer   nsx
      parameter(nsx=mw)


c *CL*
      real           sh_int(mx,my,nsx)
      real           tairDY_int(mx,my,nsx), gradTM(mx,my)
      real           qvDY_int(mx,my,nsx)  , gradQM(mx,my)
      common /interp/sh_int, tairDY_int,gradTM
     .                       , qvDY_int,gradQM

      real           mskSNo(mx,my,nsx) !*CL*
      integer        nssSNo(mx,my,nsx),nisSNo(mx,my,nsx)
      integer        issSNo(mx,my,nsx)
      integer        nhsSNo(mx,my,nsx,nsno)
      common /snotyp/nssSNo           ,nisSNo           
     .              ,issSNo
     .              ,nhsSNo
C +...mskSNo: Snow/Ice Type Index     
C +   nssSNo: Number of              Snow Layers
C +   nisSNo: Number of              Ice  Layers
C +   issSNo: Number of superimposed Ice  Layers
C +   nhsSNo: Layer historics: previous occurrence of water or faceted crystal

      real           tisSNo(mx,my,nsx,nsno),wasSNo(mx,my,nsx,nsno)
      real           rosSNo(mx,my,nsx,nsno),dzsSNo(mx,my,nsx,nsno)
      real           g1sSNo(mx,my,nsx,nsno),g2sSNo(mx,my,nsx,nsno)
      real           agsSNo(mx,my,nsx,nsno),SWaSNo(mx,my,nsx)
      real           SWSSNo(mx,my,nsx)   ,SWfSNo(mx,my)
      common /snopro/tisSNo              ,wasSNo
     .              ,rosSNo              ,dzsSNo
     .              ,g1sSNo              ,g2sSNo
     .              ,agsSNo              ,SWaSNo
     .              ,SWSSNo              ,SWfSNo, mskSNo
C +...tisSNo: Soil temperature      (=> in the snow cover )                [K]
C +   wasSNo: Soil humidity content (=> in the snow cover )            [kg/kg]
C +   rosSNo: Snow density                                             [kg/m3]
C +   dzsSNo: Depth of each layer       of the snow cover                  [m]
C +   g1sSNo: Grain Type: -99< g1 < 0 (Dendricity)      0< g1 <99 [Sphericity]
C +   g2sSNo: Grain Type:   0< g2 <99 (Sphericity)  .3mm < g2           [Size]
C +   agsSNo: Snow Age                                                   [day]
C +   SWaSNo: Surficial Water Mass                             [mmWE or kg/m2]
C +   SWSSNo: Surficial Water Status              (0,1 = no freezing,freezing)
C +   SWfSNo: Surficial Water Hydrological Model                           [-]

      real           snohSN(mx,my,nsx)
      real           BrosSN(mx,my,nsx)
      real           BG1sSN(mx,my,nsx)
      real           BG2sSN(mx,my,nsx)
      real           blowSN(mx,my,nsx)
      real           SaltSN(mx,my,nsx)
      real           hSalSN(mx,my,nsx)
      real           zWE0SN(mx,my,nsx)
      real           zWE_SN(mx,my,nsx)
      real           zWEcSN(mx,my,nsx)
      common /snofrc/snohSN,BrosSN,BG1sSN,BG2sSN,blowSN,SaltSN,hSalSN
     .              ,zWE0SN,zWE_SN,zWEcSN
C +...snohSN:       Snow Buffer Layer Thickness                          [mWE]
C +   BrosSN:       Snow Buffer Layer Density                          [kg/m3]
C +   BG1sSN:       Snow Buffer Layer Dendricity / Sphericity              [-]
C +   BG2sSN:       Snow Buffer Layer Sphericity / Size                    [-] [0.1 mm]
C +   blowSN: Blown snow Buffer Layer Thickness                          [mWE]
C +   SaltSN: Friction Velocity Saltation Threshold                      [m/s]
C +   hSalSN: Saltating    Snow Thickness                                [mWE]
C +   zWE0SN: Initial      Snow Thickness                               [mmWE]
C +   zWE_SN: Current      Snow Thickness                               [mmWE]
C +   zWEcSN: Non-erodible Snow Thickness                               [mmWE]

      real           smbsSN(mx,my,nsx)  ! Surface Mass Balance           [mWE]
      real           sihsSN(mx,my,nsx)  ! Superimposed Ice Height          [m]
      real           sshsSN(mx,my,nsx)  ! Slush            Height          [m]
      real           ravsSN(mx,my,nsx)  ! Averaged Density             [kg/m2]
      real           wavsSN(mx,my,nsx)  ! Averaged Water  Content      [kg/kg]
      real           WEq_SN(mx,my,nsx)  ! Added    Snow   Mass          [mmWE]
      common /snobal/smbsSN,sihsSN,sshsSN,ravsSN,wavsSN,WEq_SN
