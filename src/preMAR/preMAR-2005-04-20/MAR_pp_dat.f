      block data MAR_pp_dat

C +------------------------------------------------------------------------+
C |                                                                    MAR |
C |                                                                        |
C |                                                                        |
C |   Block Data MAR_pp_dat contains the list of MAR Options               |
C |                                                                        |
C |                                                                        |
C |                                                                        |
C +------------------------------------------------------------------------+

      include 'MAR_pp.inc'


      data (labval(inpspe),comval(inpspe),inpspe= 1,96)/
C +          123456789012345678901234567890123456789012345678901234567890
     .'#DP','DOUBLE PRECISION -- DOUBLE PRECISION -- DOUBLE PRECISION ',
     .'#LA','LAPACK LIBRARY is EXTERNAL // LAPACK LIBRARY is EXTERNAL ',
     .'#//','Parallelisation Set Up is activated   (software MPI used)',
     .'#HF','Initialisation of Huang and Lynch 1993   (HAMMING Filter)',

     .'#NH','DYNAMICS: Non-Hydrost. Code  (adapted from Laprise, 1992)',
     .'#nh','DYNAMICS: Non-Hydrost. Code  (Slope            Contribut)',
     .'#DH','DYNAMICS: Non-Hydrost. Code  (Diabatic Heating Contribut)',
     .'#ON','DYNAMICS: Non-Hydrost. Corr. (Weisman &al.1997 MWR p.541)',
     .'#DD','DYNAMICS: Mass Divergence Damper (Skamarock &Klemp, 1992)',
     .'#VN','DYNAMICS: Variable Number of Leap-Frog Loops (Fast Waves)',
     .'#IL','DYNAMICS: PGF: SBL Fraction with  Air = Surface Temperat.',
     .'#GE','DYNAMICS: Geographic Coordinates may be red in MARdom.dat',
     .'#CC','DYNAMICS: Constant Coriolis Parameter = fcorDY(imez,jmez)',
     .'#2Z','DYNAMICS: Zonally Averaged Model (latitude = x-direction)',
     .'#HE','DYNAMICS: DNMI   Model Vertical Discretisation(29 Levels)',
     .'#lm','DYNAMICS: LMDZ   Model Vertical Discretisation(11 Levels)',
     .'#PA','DYNAMICS: Parish Model Vertical Discretisation(10 Levels)',
     .'#PV','DYNAMICS: Large Scale Flow conserves Pot. Vort. (2D ONLY)',
     .'#pv','DYNAMICS: Large Scale Flow conserves Pot. Temp. (2D ONLY)',
     .'#UW','DYNAMICS: Advect.  3rd Accurate in Space Upstream Scheme ',
     .'#UP','DYNAMICS: Vertical Upstr. Advection Scheme (1st accuracy)',
     .'#ZU','DYNAMICS: Vertical Advection: Cubic Spline (4th accurate)',
     .'#ZO','DYNAMICS: Vertical Advection: Cubic Spline (+Open  SrfBC)',
     .'#UR','DYNAMICS: Vertical Advection/ Upper Radiating Bound.Cond.',
     .'#EP','DYNAMICS: Lateral Sponge included in    Horizontal Filter',
     .'#RB','DYNAMICS: Lateral BC: Carpenter(1982) Sommerfeld Modified',
     .'#DA','DYNAMICS: Lateral BC: Davies (1976) BC on Wind // Lat.B. ',
     .'#da','DYNAMICS: Lateral BC: Davies (1976) BC: K, nu  computed. ',
     .'#FB','DYNAMICS: Lateral BC: Fixed in Horizontal Cubic Spline   ',
     .'#OB','DYNAMICS: Lateral BC: Zero Gradient                      ',
     ,'#OG','DYNAMICS: Lateral BC: (Left) NO Nudging if relaxg=.false.',
     .'#ob','DYNAMICS: Lateral BC: Zero Gradient (Subroutine LBC000)  ',
     .'#RF','DYNAMICS: Top BC: Rayleight Friction in the Top Sponge   ',
     .'#Di','DYNAMICS: Top BC: Dirichlet  (fixed)                     ',
     .'#V+','DYNAMICS: Top BC: Von Neuman (prescrib.non zero-gradient)',
     .'#PS','DYNAMICS: Domain Averaged Pressure Thickness   maintained',
     .'#DY','DYNAMICS: OUTPUT: Components  lowest Level Forces Balance',

     .'_PE','DIFFUSION:(%Grad.)   Slope      USE+ _HH or     (_HH #CR)',
     .'#PE','DIFFUSION:(%Deform.) Slope      USE+ #DF or (#DF #DC #CR)',
     .'_HH','DIFFUSION:(%Grad.)   Vert.Cor.  USE+ _PE                 ',
     .'#DF','DIFFUSION:(%Deform.) Vert.Cor.  USE+ #PE or (#PE #DC #CR)',
     .'#DC','DIFFUSION:(%Deform.)            USE+        (#DF #PE #CR)',
     .'#CR','DIFFUSION: Cross Corr.    USE+ (_PE _HH) or (#DF #PE #DC)',

     .'#FE','FILTERING: Digital Filtering of TKE                      ',
     .'#fe','FILTERING: Digital Filtering of TKE  is   not vectorized ',
     .'#FO','FILTERING: Digital Filtering of TKE (zero gradient at LB)',
     .'#KS','FILTERING: Upper Sponge is solved by horizontal filtering',

     .'#BR','TURBULENCE: 2.5 Level  2nd Order  (Brasseur         1997)',
     .'#CA','CONVECTIVE  Adjustment (general                 Set Up)  ',
     .'#cA','CONVECTIVE  Adjustment (no Double Counting      Set Up)  ',
     .'#ca','CONVECTIVE  Adjustment (no Vector               Set Up)NV',
     .'#FC','CONVECTIVE  Adjustment (Fritsch & Chappell 1980 Set Up)  ',
     .'#fc','CONVECTIVE  Adjustment (Fritsch & Chappell 1980 Set Up)NV',
     .'#kf','CONVECTIVE  Adjustment (Kain    & Fritsch  1990 Improvm.)',
     .'#IT','CONVECTIVE  Adjustment (over 5km Adiabatics Starting Pts)',
     .'#AN','CONVECTIVE  Adjustment (Subgrid Mountain Breeze included)',
     .'#WD','CONVECTIVE  Adjustment (Water Detrainment       included)',
     .'#CG','CONVECTIVE  Adjustment (Cloud Glaciation        included)',
     .'#ND','CONVECTIVE  Adjustment (No Precip if LevFSink<LiftCond.L)',
     .'#vT','CONVECTIVE  Adjustment (Virtual Temperature  is computed)',
     .'#PB','CONVECTIVE  Adjustment (Peter Bechtold     2000 Set Up)  ',
     .'#pb','CONVECTIVE  Adjustment (Peter Bechtold     2000 Set Up)NV',
     .'#KE','CONVECTIVE  Adjustment (Emanuel & Zivkovic 1999 Set Up)  ',
     .'#ke','CONVECTIVE  Adjustment (Emanuel & Zivkovic 1999 Set Up)NV',

     .'#TA','TURBULENCE: K-e: Dissipation + Advect.Horiz.TKE Transport',
     .'#TD','TURBULENCE: K-e: Dissipation + Diffus.Horiz.TKE Transport',
     .'#AV','TURBULENCE: K-e: Buoyancy includes      Aerosol Loading  ',
     .'#HR','TURBULENCE: K-e: Huang & Raman              (1991) BLM 55',
     .'#KI','TURBULENCE: K-e: Kitada                     (1987) BLM 41',
     .'#BH','TURBULENCE: K-e: Kitada (modified)           USE with #KI',
     .'#Kl','TURBULENCE: K-l: Therry & Lacarrere         (1983) BLM 25',
     .'#LE','TURBULENCE: K  : Louis                      (1979) BLM 17',
     .'#KC','TURBULENCE: T.K.E.(mz1) := T.K.E.(mz)                    ',
     .'#KA','TURBULENCE: T.K.E. & e(T.K.E.) Filter along the vertical ',
     .'#AM','TURBULENCE: u*   Time Mean (BOX Moving Average)          ',
     .'#AT','TURBULENCE: u*T* Time Mean (BOX Moving Average)          ',
     .'#AS','TURBULENCE: u*s* Time Mean (BOX Moving Average)          ',
     .'#De','TURBULENCE: Top BC: Dirichlet (fixed) (ect_TE and eps_TE)',
     .'#AE','TURBULENCE: Aerosols Erosion / Turbulent Diffusion Coeff.',
     .'#WE','TURBULENCE: T.K.E. OUTPUT on File MAR.TKE                ',

     .'#BU','SBL: Univ.Funct.:    Businger (1973)  USE with _NO OR #NO',
     .'_NO','SBL: Univ.Funct.: NO Noilhan  (1987)  USE with #BU OR #DR',
     .'#NO','SBL: Univ.Funct.:    Noilhan  (1987)  USE with #BU       ',
     .'#DR','SBL: Univ.Funct.:    Dyer     (1974)  USE with _NO       ',
     .'#LP','SBL: Blowing Snow Fric. Veloc. Thr. (Li and Pomeroy 1997)',
     .'#DS','SBL: Blowing Snow SBL   Flux   (analytical Form of dq/dz)',
     .'#OR','SBL: Orography Roughness included from SL_z0 in MARdom   ',
     .'#ZS','SBL: Mom.: Roughn.Length= F(u*) Chamberlain (1983),  Sea ',
     .'#ZN','SBL: Mom.: Roughn.Length= F(u*) Andreas &al.(2004), Snow ',
     .'#RN','SBL: Heat: Roughn.Length= F(u*,z0)  Andreas (1987), Snow ',
     .'#ZM','SBL: M/H   Roughn.Length: Box Moving Average (in Time)   ',
     .'#SB','Surface Boundary: modified externally (from Campain Data)',
     .'#TI','Turbul. Heat Surface Flux: Implicit numerical Scheme     ',
     .'#QE','Turbul. Vap. Surface Flux: Explicit numerical Scheme     ',
     .'#FI','Turbul. Mom. Surface Flux: Implicit numerical Scheme     ',
     .'#BI','Blowing Snow Surface Flux: Implicit numerical Scheme     '/

      data (labval(inpspe),comval(inpspe),inpspe=97,nnppar)/
     .'#OL','TEST:      Linear Mountain Wave: Specific IO    (2D ONLY)',
     .'#OM','TEST: (Non)Linear Mountain Wave: Specific INPUT (2D ONLY)',
     .'#OS','TEST:      Linear Mountain Wave: Specific IO    (2D ONLY)',
     .'#GR','TEST: LBC: 10C/d  Atmos.Cooling at Model Center (2D ONLY)',
     .'#K1','TEST: LBC: Katab. Atmos.Warming                 (1D ONLY)',
     .'#EK','TEST: EKMAN Spiral: Constant Vertical Turbul. Coefficient',
     .'#CL','TEST: Convective Mixed Layer Test         (HS = 100 W/m2)',
     .'#NL','TEST: Nearly   Neutral Layer Test         (HS =   0 W/m2)',

     .'#TC','TraCer   Advection-Diffusion Equation        is turned ON',
     .'#tc','TraCer   Filtering is  not vectorized                    ',
     .'#TO','TraCer   Open Lateral Boundary Conditions on digit.Filter',
     .'#TS','TraCer   Tracer Deposition diagnostic        is turned ON',
     .'#BD','TraCer   Aeolian Erosion  Submodel           is turned ON',
     .'#DV','TraCer   Aeolian Erosion  Submodel: Air Loading by Dust  ',
     .'#CH','Chemical Atmospheric         Model       may be turned ON',
     .'#MV','TraCer   Total Mass          Verification    is turned ON',

     .'#HY','Explicit Cloud MICROPHYSICS              may be turned ON',
     .'#hy','Explicit Cloud MICROPHYSICS: NO Vectorisation Optmization',
     .'#HM','Explicit Cloud MICROPHYSICS: Hallett-Mossop Ice Multipl. ',
     .'#hm','Explicit Cloud MICROPHYSICS: Hallett-Mossop Ice Mult.  NV',
     .'#LI','Explicit Cloud MICROPHYSICS: Lin et al. (1983) Autoconv. ',
     .'#BS','Explicit Cloud MICROPHYSICS: Blow. *(Snow)         Model ',
     .'#SS','Explicit Cloud MICROPHYSICS: Blow. *(Snow)  Linear Model ',
     .'#S0','Explicit Cloud MICROPHYSICS: Blow. *(Byrd)  Linear Model ',
     .'#HV','Explicit Cloud MICROPHYSICS: Air Loading by Hydrometeors ',
     .'#BV','Explicit Cloud MICROPHYSICS: SBL Loading by all Water Sp.',
     .'#bv','Explicit Cloud MICROPHYSICS: SBL Loading not vectorized  ',
     .'#EM','Explicit Cloud MICROPHYSICS: de Montmollin Parameterizat.',
     .'#BW','Explicit Cloud MICROPHYSICS: Blowing Snow Statistics     ',
     .'#b2','Explicit Cloud MICROPHYSICS: Blowing Snow Statistics (II)',
     .'#EV','Explicit Cloud MICROPHYSICS: Snow Erosion Statistics     ',
     .'#HW','Explicit Cloud MICROPHYSICS: OUTPUT of qr,qs, qw,qi on NC',
     .'#EW','Explicit Cloud MICROPHYSICS: OUTPUT (Ener./Mass) (Unit 6)',
     .'#ew','Explicit Cloud MICROPHYSICS: OUTPUT (Ener./Mass) (Unit 6)',
     .'#WH','Explicit Cloud MICROPHYSICS: OUTPUT              (Unit 6)',
     .'#WQ','Explicit Cloud MICROPHYSICS: OUTPUT (Full Verif) (Unit 6)',
     .'#WB','Explicit Cloud MICROPHYSICS: Water Conservation Controled',
     .'#WW','Explicit Cloud MICROPHYSICS: Water Conservation Summary  ',
     .'#VC','Explicit Cloud MICROPHYSICS: Water Conservation is Forced',
     .'#HO','Explicit Cloud MICROPHYSICS: Zero-Gradient Lat.Bound.Cond',

     .'#MR','PHYSICS: MARrad: Solar/Infrared     (Laurent LI set up)  ',
     .'#AZ','PHYSICS: Solar : Direct Radiation:   Surface Slope Impact',
     .'#MM','PHYSICS: Solar : Direct Radiation:   Mountains Mask    ON',
     .'#TR','PHYSICS: Solarn: Clear Sky, without Underlying Reflection',
     .'#EE','PHYSICS: radCEP: ECMWF   routine    (cfr. JJ Morcrette)  ',
     .'#LL','PHYSICS: radLMD: radlwsw routine    (Laurent LI set up)  ',
     .'#ll','PHYSICS: radLMD: radlwsw routine    (Laurent LI set up)NV',
     .'#AR','PHYSICS: radLMD: radlwsw routine Interactive Terr.Aerosol',
     .'#WL','PHYSICS: radLMD: radlwsw routine IO (Laurent LI set up)  ',

     .'#SA','PHYSICS: MAR Code behaves  as a Stand Alone Surface Model',

     .'#FR','Surface Model: Force Restore (Deardorff) at least   is ON',
     .'#WG','Soil Humidity: Force Restore (Deardorff) may be turned ON',

     .'#TV','Soil /Vegetation Model                   may be turned ON',
     .'#GP','Soil /Vegetation Model: LAI, GLF Variations NOT prescrib.',
     .'#LN','Soil /Vegetation Model: LAI(x,y,t) prescribed(MARglf.DAT)',
     .'#SV','Soil /Vegetation Model  (Koen De Ridder) may be turned ON',
     .'#SH','Soil /Vegetation Model: Hapex-Sahel   Vegetation     DATA',
     .'#V1','Soil /Vegetation Model: (KD) Vegetat. IGBP Classification',
     .'#V2','Soil /Vegetation Model: (KD) Vegetat. MAR  Classification',

     .'#GA','SISVAT: Soil Humidity Geometric Average at Layer Interfac',
     .'#GF','SISVAT: Gravitational Saturation Front          turned ON',
     .'#GH','SISVAT: Gravitational Saturation Front - Horton turned ON',
     .'#OP','SISVAT: Interactive Sea Surface Temperature     turned ON',
     .'#op','SISVAT: SST Nudging -->   prescribed values     turned ON',
     .'#IP','SISVAT: Sea-Ice Fraction prescribed from SMMR and SSM/I  ',
     .'#SI','SISVAT: Sea-Ice Fraction calculated from prescribed SST  ',
     .'#MT','SISVAT: Monin-Obukhov Theory is linearized (Garrat schem)',
     .'#SR','SISVAT: traces & OUTPUT a variable among called routines ',
     .'#WV','SISVAT: performs OUTPUT on an ASCII File (1 file each pt)',
     .'#sa','SISVAT: must be pre-processed, except in stand-alone run ',

     .'#ST','EVOLUTIVE SST (Sea Surface Temperature/Swab Ocean)       ',
     .'#RE','PRESCRIB. SST (Sea Surface Temperature/Reynolds DATA Set)',
     .'#PO','POLYNYA Model                            may be turned ON',
     .'#FD','POLYNYA Model: Sea-Ice Velocity is Free Drift            ',
     .'#HA','POLYNYA Model: POLYNYA Surface Energy Balance:  2000 W/m2',
     .'#HI','POLYNYA Model: Hibler (1979) Parameteriz. of Ice Strength',
     .'#CN','POLYNYA Model: Prescription of a Local Avective Time Step',

     .'#SN','SNOW Model                               may be turned ON',
     .'#AB','SNOW Model: Interactive Albedo f(Grain) (Brun et al.1991)',
     .'#AG','SNOW Model: Snow Aging Col de Porte     (Brun et al.1991)',
     .'#DG','SNOW Model: Snow Settling when Melting | Minimum Density ',
     .'#Se','SNOW Model: Energy Conserv. Verific.: Summary, Output    ',
     .'#SE','SNOW Model: Energy Conserv. Verific.: Summary, Output++++',
     .'#SF','SNOW Model: Energy Conserv. Verific.: Forcing, Conduction',
     .'#SW','SNOW Model: Water  Conserv. Verific.: Melting, Freezing  ',
     .'#HS','SNOW Model: Hardened SNOW Pack Initialization            ',
     .'#RU','SNOW Model: Slush:  Internal Run OFF of Water Excess     ',
     .'#GK','SNOW Model: Interactive Albedo (Greuell &Konzelmann 1994)',
     .'#SL','SNOW Model: Interactive Albedo (Zuo     &Oerlemans  1995)',
     .'#SM','SNOW Model: Melting/Freezing Diagnostics                 ',
     .'#SZ','SNOW Model: Z0 Dependance on varying Sastrugi Height     ',
     .'#CP','SNOW Model: For Validation on Col de Porte Data          ',
     .'#GL','SNOW Model: ETH-Camp & Greenland 3D simulations          ',

     .'#CS',' INPUT: Constant Sounding during 1st Hours     (2-D ONLY)',

     .'#IB','OUTPUT: Ice-Sheet Surface Mass Balance  (on MARphy File )',
     .'#ID','OUTPUT: Main Dependant Variables        (on NetCDF File )',
     .'#UL','OUTPUT: Time Dimension is UNLIMITED     (on NetCDF File )',
     .'#T2','OUTPUT: 2-m  Air Temperature            (on NetCDF File )',
     .'#MA','OUTPUT: MesoAnimation                   (on NetCDF Files)',
     .'#W6','OUTPUT, Additional: Simulation Statistics      on MAR.log',
     .'#w6','OUTPUT, Additional: Simulation Statistics (NH) on MAR.log',
     .'#WA','OUTPUT, Additional:                            DYNadv_ver',
     .'#WR','OUTPUT, Additional: INIsnd, infra, SRFmod_sno, SRFmod_pol',

     .'#vL','PORTABILITY: Vectorization enhanced                      ',
     .'#vN','PORTABILITY: Vectorization enhanced: Leap Frog Counter   ',
     .'#vK','PORTABILITY: Vectorization enhanced: TKE                 ',
     .'#vH','PORTABILITY: Vectorization enhanced: Hydrological Cycle  ',
     .'#vB','PORTABILITY: Vectorization enhanced: Blowing  Snow *     ',
     .'#vD','PORTABILITY: Vectorization enhanced: Blowing  Dust .     ',
     .'#vR','PORTABILITY: Vectorization enhanced: Sastrugi Height     ',
     .'#vS','PORTABILITY: Vectorization enhanced: Snow     Model      ',
     .'#vV','PORTABILITY: Vectorization enhanced: SVAT                ',
     .'#vZ','PORTABILITY: Vectorization enhanced: Av.Roughness Length ',
     .'#HP','PORTABILITY: Enables use of own    library on Linux syst.',
     .'#NV','PORTABILITY: Vectorization  is     turned  OFF           ',
     .'   ','                                                         ',
     .'   ','                                                         ',
     .'   ','                                                         ',
     .'   ','                                                         '/


      end
