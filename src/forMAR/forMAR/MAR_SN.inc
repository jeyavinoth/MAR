C +
C +  ==========
C +--MAR_SN.inc
C +  ==========
C +
      integer       nhSNo1,nhSNo2,nhSNo3,nhSNo4,nhSNo5
      common/snoint/nhSNo1,nhSNo2,nhSNo3,nhSNo4,nhSNo5
C +...              nhSNo1,nhSNo2,nhSNo3,nhSNo4,nhSNo5=1,2,3,4,5: Cas Histor.
C +
      real                 h0SNo, zxSNo, 
     .             gSNmax,vg1SNo, dtSNh, taSNo, tfSNo, teSNo, roSNf,
     .              w0SNo, fcSNo,roiwSN, c1SNo, c2SNo, c3SNo,
     .              waSN ,viwaSN,gethSN, epSN1,
     .             aiceSN,asluSN,DiEtSN,DiGfSN,DiFcSN,
     .             den1SN,den2SN,den3SN,den4SN,
     .             sph1SN,sph2SN,sph3SN,sph4SN,
     .             rho1SN,rho2SN,rho3SN,rh10SN,rh11SN,rh12SN,
     .             SNmdz1,SNmdz2,SNmdz3,SNmdz4,
     .             SNmdz5,SNmdz6,SNmdz7,SNmdz8,
     .             SNtyp1,SNtyp2,SNtyp3,SNtyp4,SNtyp5,SNtyp6,
     .             SNdag1,SNdag2,SNdag3,SNdag4,SNdag5,SNdag6
      common/snodat/       h0SNo, zxSNo, 
     .             gSNmax,vg1SNo, dtSNh, taSNo, tfSNo, teSNo, roSNf,
     .              w0SNo, fcSNo,roiwSN, c1SNo, c2SNo, c3SNo,
     .              waSN ,viwaSN,gethSN, epSN1,
     .             aiceSN,asluSN,DiEtSN,DiGfSN,DiFcSN,
     .             den1SN,den2SN,den3SN,den4SN,
     .             sph1SN,sph2SN,sph3SN,sph4SN,
     .             rho1SN,rho2SN,rho3SN,rh10SN,rh11SN,rh12SN,
     .             SNmdz1,SNmdz2,SNmdz3,SNmdz4,
     .             SNmdz5,SNmdz6,SNmdz7,SNmdz8,
     .             SNtyp1,SNtyp2,SNtyp3,SNtyp4,SNtyp5,SNtyp6,
     .             SNdag1,SNdag2,SNdag3,SNdag4,SNdag5,SNdag6
C +...              h0SNo:  Thickness     of the upper Snow Layer        [3cm]
C +                 zxSNo:  Maximum Depth of the       Snow Model       [100m]
C +                gSNmax:  Maximum Grain Size                           [3mm]
C +                vg1SNo:  Conversion 0/99-->0/1 pour Sphericite et Dendricite
C +                 dtSNh= 3600. sec:   Time Step for Settling and Metamorphism
C +                 taSNo=60000, tfSNo=0.3, teSNo=23: Snow Settling Parameters
C +                        (Brun et al. 1989, JGl 35 121,                p.334)
C +                 roSNf:  Density of Fresh snow                  [300 kg/m3]
C +                        (Greuell and Konzelmann 1994,  Glob. Plan. Change 9)
C +                 w0SNo:  Snow Layer Maximum Water Content      [0.04 kg/kg]
C +                        (Rowe et al. 1995, JGR 100 D8, p. 16267, col. 1)
C +                 fcSNo= -7.8 (Male 1980, (96) p.376)
C +                roiwSN:  Density Ratio (Ice / Water)                [0.920]
C +                 c1SNo, c2SNo, c3SNo: Water Run Off Parameters 
C +                        (Zuo and Oerlemans 1995 J.Glacio., in press)
C +                  waSN:  Melt Water Scaling Factor (Slush Albedo)   [200mm]
C +                        (Zuo and Oerlemans 1995 J.Glacio., in press,p.24)
C +                viwaSN:  Kinematic Water Viscosity at 20.C  [0.000001 m2/s]
C +                gethSN:  Geothermal Heat Flux                        [W/m2]
C +                 epSN1:  0.1
C +                aiceSN:  0.55 (Ice   Albedo)
C +                asluSN:  0.45 (Slush Albedo)
C +                DiEtSN:  1.00 (Diametre opt. Etoile     / + Opt. Diameter)
C +                DiGfSN:  3.00 (Diametre opt. Grain fin  / . Opt. Diameter)
C +                DiFcSN:  4.00 (Diametre opt. Face plane / # Opt. Diameter)
C +                den1SN,den2SN,den3SN,den4SN: Dendr.Paramet. / Wind Effect
C +                sph1SN,sph2SN,sph3SN,sph4SN: Spher.Paramet. / Wind Effect
C +                rho1SN,rho2SN,rho3SN,rh10SN,rh11SN,rh12SN --> Snow Density
C +                nhSNo1,nhSNo2,nhSNo3,nhSNo4,nhSNo5=1,2,3,4,5: Cas  Histor.
C +                SNmdz1,2,3,4,5,6,7,8: for aggregation criterium   (    dz)
C +                SNtyp1,2,3,4,5,6    : for aggregation criterium   (* type)
C +                SNdag1,2,3,4,5,6    : for aggregation criterium   (*  age)
C +
      integer   mg0  ,mg   ,mgg     ,mg1     ,mh
      parameter(mg0=0,mg=20,mgg=mg+1,mg1=mg-1,mh=0)

      integer   mg_SN1     ,mg_SN2      ,mg_SN3
      parameter(mg_SN1=mg/5,mg_SN2=mg/10,mg_SN3=mg/5)
C +...          mgs:    Snow  Layers reference Number
C +             mg0:    Ice   Layers reference Number
C +             mg:     Total Layers reference Number
C +
      integer        maskSN
      common /snomsk/maskSN(mx,my)
C +...maskSN: :=0 (Melting occurred)  ;
c              =1 (Fresh Snow==>Blowing Snow Possible)
C +
      integer        nsiiSN,       nsluSN    
      common /snolay/nsiiSN(mx,my),nsluSN(mx,my)  
      integer         nsSNow          ,niSNow       
      common /snolay2/nsSNow(mx,my)   ,niSNow(mx,my)
C +...nsSNow: Number of Snow Layers
C +   niSNow: Number of Ice  Layers
      real           agSNow(mx,my,mg),waSNru(mx,my)  , slwaSN(mx,my)
     .               ,  slopSN(mx,my),snSNow(mx,my)  ,avwatSN(mx,my)
     .               , avdenSN(mx,my),slushSN(mx,my) ,siiceSN(mx,my)
     .               , smbalSN(mx,my),smbalSN0(mx,my)
     .               , znSNow0(mx,my),znSNow_snow(mx,my)
     .               , znSNow(mx,my)
      common /snopar/agSNow          ,waSNru         , slwaSN
     .              ,slopSN          ,snSNow         ,avwatSN
     .              ,avdenSN         ,slushSN        ,siiceSN
     .              ,smbalSN         ,smbalSN0       
     .              ,znSNow0         ,znSNow_snow    
     .              ,znSNow
C +...agSNow: Snow Age                                   
C +   waSNru: Surficial Water Specific Mass                    (mmWE or kg/m2)
C +   slwaSN: Normalized Surficial Water Variation during dtSNo,   fct(slopSN)
C +---------------------------------------------------------------------------
C +   slopSN: Fall Line Slope
C +   itSNow: Numbre of Snow    Layers Used  (Should always be <=mg)
C +   iiSNow: Numbre of Old Ice Layers Used  (Should always be <=15)
C +   snSNow: Snow Amount above Ice                                     (mmWE)
C +   avwatSN:Average Water Content in snow pack (excluding slush)         (%)
C +   avdenSN:Average Density of snow pack (excluding slush)           (kg/m3)
C +   slushSN:Slush Height of Snow Pack                                    (m)
C +   siiceSN:Superimposed Ice Height                                    (m)
C +   smbalSN: Mass Balance                                             (mmWE)
C +   smbalSN0 : Initial Mass Balance                                   (mmWE)
C +   znSNow_snow : thickness of snow/SI pack above itexpe=0 old ice       (m)
C +   znSNow0  : itexpe=0 old ice's thickness                              (m)
C +---------------------------------------------------------------------------
      integer        nhSNow
      common /snohis/nhSNow(mx,my,mg)
C +...nhSNow: Layer historics: previous occurrence of water or faceted crystal
C +
      real           tiSNow               ,waSNow,
     .               roSNow               ,dzSNow,
     .               g1SNow               ,g2SNow
      common /snowzz/tiSNow(mx,my,-mh:mgg),waSNow(mx,my,mg),
     .               roSNow(mx,my,    mg) ,dzSNow(mx,my,mg),
     .               g1SNow(mx,my,    mg) ,g2SNow(mx,my,mg)

c +...tiSNow: Soil temperature      (=> in the snow cover )                (K)
c +   waSNow: Soil humidity content (=> in the snow cover )            (kg/kg)
c +   roSNow: Snow density                                             (kg/m3)
c +   dzSNow: Depth of each layer       of the snow cover                  (m)
c +   g1SNow: Grain Type: -99< g1 < 0 (Dendricity)      0< g1 <99 (Sphericity)
c +   g2SNow: Grain Type:   0< g2 <99 (Sphericity)  .3mm < g2           (Size)
c +   nhSNow: Layer historics: previous occurrence of water or faceted crystal 
C +
      integer       ntSNo         ,ntwaSN       ,nhSNo
      common/snoint/ntSNo         ,ntwaSN       ,nhSNo(    mg)
C +
      real          dtSNo         ,dtwaSN       ,
     .              znSNo         ,dzSNo        ,dzSNr,
     .              roSNo         ,waSNo        ,tiSNo,
     .              g1SNo         ,g2SNo
      common/snoloc/dtSNo         ,dtwaSN       ,
     .              znSNo(-mh:mgg),dzSNo(-mh:mg),dzSNr(    mg),
     .              roSNo(    mg ),waSNo(    mg),tiSNo(-mh:mg),
     .              g1SNo(    mg ),g2SNo(    mg)
C +...              dtSNo: Snow model Time step
C +                           dtwaSN: Snow Water Percolation Time Step
C +                                     ntSNo: Snow model: Nb local Iter.
C +                                               ntwaSN: Water Percol.:idem
C +...              znSNo,    dzSNo: Current   Discretisation
C +                           dzSNr: Reference Discretisation
C +
