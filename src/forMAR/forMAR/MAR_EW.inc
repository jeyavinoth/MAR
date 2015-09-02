C +...Energy and Water Balance Variables
C +
      real          wat0EW       ,wat1EW       ,wat2EW       ,
     .              watfEW       ,
     .              enr0EW       ,enr1EW       ,enr2EW
      common/hyd_EW/wat0EW(mx,my),wat1EW(mx,my),wat2EW(mx,my),
     .              watfEW(mx,my),
     .              enr0EW(mx,my),enr1EW(mx,my),enr2EW(mx,my)
C +...              wat*EW: Total Precipitable  Water  in the  Air Column [m]
C +                 watfEW: Water Flux (Atm. --> Srf.) during 1 Time Step [m]
C +                 enr*EW: Total Energy (Sens. +Lat.) in the  Air Column [m]
C +
      character*20  mphyEW(mx,my),mauxEW
C +   common/hydcEW/mphyEW(mx,my),mauxEW
C +
