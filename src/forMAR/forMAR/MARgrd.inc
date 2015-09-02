
C +  ==========
C +--MARgrd.inc
C +  ==========

      real            dt    ,dt_Loc,dtquil,tequil
      real            dtfast,dtAdvH,dtDiff,dtPhys,dtRadi,dthyd
      common /tgrid/  dt    ,dt_Loc,dtquil,tequil
     .             ,  dtfast,dtAdvH,dtDiff,dtPhys,dtRadi,dthyd 
C +...                dt    : Time Step for Slow Dynamics                  (s)
C +                   dt_Loc: Time Step (dummy)                            (s)
C +                   dtquil: Time Step for 1-D Initialisation             (s)
C +                   tequil: Time Span for 1-D Initialisation             (s)
C +...                dtfast: Time Step for Fast Dynamics                  (s)
C +                   dtAdvH: Time Step for Tracer  Advection              (s)
C +                   dtDiff: Time Step for Subgrid Scale Processes        (s)

      integer idt,jdt,ntFast,nt_Mix,ntAdvH,ntDiff,ntPhys,ntRadi,nthyd
      integer         itFast,it_Mix,jtAdvH,jtDiff,jtPhys,jtRadi,jtRadi2
      integer         kssMAR,jdplus,mmplus
      common /tgridi/ idt   ,jdt
     .              , nt_Mix,ntFast,ntAdvH,ntDiff,ntPhys,ntRadi,nthyd
     .              , it_Mix,itFast,jtAdvH,jtDiff,jtPhys,jtRadi,jtRadi2
     .              , kssMAR,jdplus,mmplus

      integer*8       jdaMAR,jhaMAR,jmmMAR,jssMAR,jhaRUN,dt_new,dt_old
      common /tgridi8/jdaMAR,jhaMAR,jmmMAR,jssMAR,jhaRUN,dt_new,dt_old    

C +...                idt   : Time Step for Slow Dynamics                  (s)
C +                   ntFast: Nb Fast Dyn.Tim.Steps for one Mix  Dyn.Tim.Step
C +                   itFast: No Fast Dyn.Tim.Step
C +                   nt_Mix: Nb Mix  Dyn.Tim.Steps for one Slow Dyn.Tim.Step
C +                   it_Mix: No Mix  Dyn.Tim.Step
C +                   ntAdvH: Nb Trac.Adv.Tim.Steps for one Slow Dyn.Tim.Step
C +                   ntDiff: Nb Turbul.  Tim.Steps for one Slow Dyn.Tim.Step
C +                   ntPhys: Nb Srf.Phys.Tim.Steps for one Slow Dyn.Tim.Step
C +                   ntPadi: Nb Rad.Phys.Tim.Steps for one Slow Dyn.Tim.Step
C +...                jtAdvH: Nb Slow Dyn.Tim.Steps for one Trac.Adv.Tim.Step
C +                   jtDiff: Nb Slow Dyn.Tim.Steps for one Turbul.  Tim.Step
C +                   jtPhys: Nb Slow Dyn.Tim.Steps for one Srf.Phys.Tim.Step
C +                   jtPadi: Nb Slow Dyn.Tim.Steps for one Rad.Phys.Tim.Step
C +...                jdaMAR: Nb of simulated Days
C +                   jhaMAR: Nb of simulated Hours   in Day    jdaMAR+1
C +                   jmmMAR: Nb of simulated Minutes in Hour   jhaMAR+1
C +                   jssMAR: Nb of simulated Seconds in Minute jmmMAR+1
C +                   kssMAR: Nb of simul. %  Second
C +                   jdplus: Day   Offset
C +                   mmplus: Month Offset

      logical         center,tsplit,staggr
      common /lgrid/  center,tsplit,staggr
C +...                center: Pressure       Spatial Scheme  centered or not
C +                   tsplit: Advection         Time Scheme  splitted or not
C +                   staggr:                Vertical  Grid staggered or not

      integer         nordps,norhyd,nortra
      common /ngrid/  nordps,norhyd,nortra
C +...                nordps: Pressure       Spatial Scheme Precision (-1 or 4)
C +                   norhyd: Advection H2O  Spatial Scheme Precision (-1 or 4)
C +                   norhyd: Advection Trac.Spatial Scheme Precision (-1 or 4)
      real            fac43,rxy
      common /rgrid/  fac43,rxy
C +...                fac43 =  4.d0 / 3.d0 
C +                   (used in 4th order horizontal difference scheme)
C +                   rxy   =  1.d-3/(mx*my)
      real            dx,dx2,dtx,dxinv,dxinv2
      common /xgrid/  dx,dx2,dtx,dxinv,dxinv2
C +...                dx    : horizontal grid size     (x direction)
C +                   dx2   =       2 dx
C +                   dtx   = dt /    dx
C +                   dxinv = 1  /    dx
C +                   dxinv2= 1  / (2 dx)

      integer         imez, m, m1, m2, m3, m4,i
      integer         im1    ,ip1    ,im2    ,ip2
      integer         im3    ,ip3    ,im4    ,ip4
      integer         mmx,mmx1,mmx2,mmx3,mmx4,mmx5,mmx6
      integer             m0x1,m0x2,m0x3,m0x4,m0x5,m0x6
      common /igrid/  imez, m, m1, m2, m3, m4,i,
     .                im1(mx),ip1(mx),im2(mx),ip2(mx),
     .                im3(mx),ip3(mx),im4(mx),ip4(mx),
     .                mmx,mmx1,mmx2,mmx3,mmx4,mmx5,mmx6,
     .                    m0x1,m0x2,m0x3,m0x4,m0x5,m0x6
C +...                mez   : origin grid point number (x direction)
C +                   im1,2 : max(i-1, 1), max(i-2, 1), etc...
C +                   ip1,2 : min(i+1,mx), min(i+2,mx), etc...

      real            dy,dy2,dty,dyinv,dyinv2
      common /ygrid/  dy,dy2,dty,dyinv,dyinv2
C +...                dy    : horizontal grid size     (y direction)
C +                   dy2   = 2 dy
C +                   dty   = dt / dy
C +                   dyinv = 1  / dy
C +                   dyinv2= 2  / dy

      integer         jmez,mn,mn1,mn2,mn3,mn4,j
      integer         jm1    ,jp1    ,jm2    ,jp2
      integer         jm3    ,jp3    ,jm4    ,jp4
      integer         mmy,mmy1,mmy2,mmy3,mmy4,mmy5,mmy6
      integer             m0y1,m0y2,m0y3,m0y4,m0y5,m0y6
      common /jgrid/  jmez,mn,mn1,mn2,mn3,mn4,j,
     .                jm1(my),jp1(my),jm2(my),jp2(my),
     .                jm3(my),jp3(my),jm4(my),jp4(my),
     .                mmy,mmy1,mmy2,mmy3,mmy4,mmy5,mmy6,
     .                    m0y1,m0y2,m0y3,m0y4,m0y5,m0y6
C +...                jmez  : origin grid point number (y direction)
C +                   jm1,2 : max(j-1, 1), max(j-2, 1), etc...
C +                   jp1,2 : min(j+1,my), min(j+2,my), etc...

      integer         k,mmz,mmz1,mmz2,km1,kp1,km2    
      common /kgrid/  k,mmz,mmz1,mmz2,
     .                km1(mz),kp1(mz),km2(mz)

      real   sigma      ,dsig_1       ,qsig_1
      real               dsig_2
      real   sigmid     ,dsigm1       ,qsigm1
      real               dsigm2       ,qsigm2,zsigma,z__SBL
      common /sgrid/
     .       sigma( mz) ,dsig_1(0:mzz),qsig_1(0:mzz),
     .                   dsig_2(  mzz),
     .       sigmid(mzz),dsigm1(  mz ),qsigm1(  mz ),
     .                   dsigm2(  mz) ,qsigm2(  mz ),
     .      zsigma( mz) ,z__SBL
C +          sigma : independant variable sigma on sigma levels (k)
C +          dsig_1:      difference ds(k+1/2)
C +          qsig_1: 1 / [difference ds(k+1/2) X 2]
C +          dsig_2:      difference ds(k+1/2) X 2
C +          sigmid: independant variable sigma between  levels (k-1/2)
C +          dsigm1:      difference ds(k)
C +          qsigm1: 1 / [difference ds(k)        ]
C +          dsigm2:      difference ds(k)     X 2
C +          qsigm2: 1 / [difference ds(k)     X 2]
C +         zsigma : height of sigma levels
C +...REM.:  THE VERTICAL DISCRETIZATION IS DEFINED in SubRoutine iniver 
C +                                                           and auxgri
      real           xxkm    ,yykm    ,sh
      common/topogr/ xxkm(mx),yykm(my),sh(mx,my)
C +... xxkm  : distance to the x axis origin                             (km) 
C +            this is generally taken as the coast 
C +            when studying land-sea interactions
C +    yykm  : distance to the y axis origin                             (km) 
C +    sh    : surface height                                             (m)
