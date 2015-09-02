      common/parin/ bclap(11),etasat(11),
     .              psisat(11),aksat(11),depth(4),
     .              alai(0:12),rstmin(0:12),
     .              rhof(0:12,4),tauf(0:12,4),height(0:12),
     .              displ(0:12),rootf(0:12,4),rplant(0:12),
     .              wvmax(0:12),taudif(0:12),dlf(0:12),rfr(0:12),
     .              beer0(0:12),sigt(0:12),z0m(0:12),z0h(0:12),
     .              csat1(11),csat2(11),rhocsl(11),
     .              ivgtyp(mx,my,3),isltyp(mx,my),ifrac(mx,my,3),
     .              mscmax(mx,my),mscw(5),nfilea(5),nfileb(5)
C +
C +   surface parameters and derived quantities
C +
      common/radiaf/rrad(0:12,2),trad(0:12,2),gamma(0:12,2),
     .              rrad0(2),trad0(2),anrad(2),anurad(2),
     .              cst1(2),cst2(2),u0(2),dL(2),
     .              sigrad(2),absv(2),absg(2)
C +
C +   vegetation radiative variables (local): most of them
C +   depend on spectral band (1=PAR,2=NIR)
C +
      common/surftp/tv(mx,my,3),psig(mx,my,3),tgr(mx,my,3),
     .              tvir(mx,my),etag(mx,my,3),psiv(mx,my,3),
     .              hstd(4),fvtr(4),wha(mx,my),walea(mx,my),
     .              wrnet(mx,my),wggr(mx,my),wtau(mx,my)
C +
C +   vegetation and soil surface temperature;
C +   soil surface tensions, leaf tension,halstead 
C +   function,transpiration extraction per soil layer
C +
      common/rain/drip(mx,my,3),rainpr(mx,my)
      logical rainwg,openf,writea(5)
      character*11 fnama
      character*15 fnamb
      character*3 ichar,jchar,mschar
C +
C +   water fluxes
C +
      common/soiltemp/tsoil(mx,my,3,4),conth(4),alpha(4)
C +
C +   soil temperatures (4 layers) and auxiliary variables
C +
      common/soilwat/eta(mx,my,3,4),wv(mx,my,3),omega(4),psi(4)
C +
C +   soil water contents (4 layers) and auxiliary variables
C +
      common/inicontrol/inexpe,initial
C +   
C +   inexpe  = initial time step
C +   initial = initial grid point and time step
