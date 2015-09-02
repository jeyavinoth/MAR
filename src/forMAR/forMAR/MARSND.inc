      real           tiSND1,tiSND2      
      common /timsnd/tiSND1,tiSND2      
C +
      integer        iSND,jSND
      common /origin/iSND,jSND
C +    iSND,jSND are the (x,y) coordinates of the sounding grid point
C +
      real            pstSND
      common /modlim/ pstSND
C +...                pstSND: initial model depth (kPa)
      integer         iyrSND,mmaSND,jdaSND,jhuSND,
     .                loSND,nSND
      common /sndpar/ iyrSND,mmaSND,jdaSND,jhuSND,
     .                loSND,nSND
C +...                iyrSND,mmaSND,jdaSND,jhuSND
C +                   Year   Month  Day    HourUT  of  the  Sounding
C +                   nSND      : Sounding No
C +                   loSND     : 1 -> This is NOT the last Sounding
      real            tpSND(0:40,2), tSND(0:40,2),qSND(0:40,2),
     .                 zSND(0:40,2), pSND(0:40,2),
     .                 fSND(0:40,2), dSND(0:40,2),zeSND    (2),
     .                uuSND(0:40  ),vvSND(0:40  )
      common /snddat/ tpSND        , tSND        ,qSND        ,
     .                 zSND        , pSND        ,
     .                 fSND        , dSND        ,zeSND       ,
     .                uuSND        ,vvSND
C +...xSND : sounding characteristics 
C +   ^x=(tp->potential temperature
C +        t->temperature, q->specific humidity, 
C +        z->Altitude, p->pressure, ff,dd-> large scale wind components)
C +                                     ze->relative vorticity
C +       uu,vv->Large Scale Wind Vector in MAR horiz. coord. (iSND,jSND)
