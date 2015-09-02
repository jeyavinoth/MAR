      integer   ntrac,nfxTC,kinTC,nPhot,nkWri
      parameter(ntrac=  28)
      parameter(nfxTC=   1)
      parameter(kinTC=  59)
      parameter(nPhot=   5)
      parameter(nkWri=   4)
C +
      logical       ldepTC,locaTC
      common/TracLO/ldepTC,locaTC
C +
      real*8          chTC,fxTC,qvTC,tairTC,suTC,raTC,cochTC
      common/TracDB/  chTC(mx,      ntrac),
     .                fxTC(         nfxTC),
     .                qvTC(mx            ),
     .              tairTC(mx            ),
     .                suTC(mx,      nPhot),
     .                raTC(mx,      kinTC),
     .              cochTC(mx,ntrac,ntrac)
C +
      real            qxTC,qsTC,uqTC,ch0TC,dchTC,efacTC,vdepTC,hhTC
      real            dt_ODE   ,dt2ODE
      common/TracRE/  qxTC(mx,my,mz,ntrac),
     .                qsTC(mx,my,   ntrac),uqTC(mx,my,ntrac),
     .               ch0TC(         ntrac,4),
     .               dchTC(         ntrac,4),
     .              efacTC(mx,my           ),
     .              vdepTC(         ntrac  ),
     .                hhTC(mx              )
C +
      integer         ikTC,krouTC           ,nt_ODE,jt_ODE
      common/TracIN/  ikTC(         nkWri  ),
     .              krouTC(mx,my           ),nt_ODE,jt_ODE
C +
      character *9   namTC,fixTC
      common/TracCH/ namTC(         ntrac),
     .                     fixTC(   nfxTC) 
C +
      real      Unity
      parameter(Unity  =   0.100e+01)
C +...          Unity:[Molec./cm3]
C +
      real      cminTC 
      parameter(cminTC =   0.100e-07)
C +
