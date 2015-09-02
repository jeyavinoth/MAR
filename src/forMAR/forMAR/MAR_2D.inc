
C +  ==========
C +--MAR_2D.inc
C +  ==========

      integer         indx(klon),indy(klon)

      common /hy_grd/ indx      ,indy

      real            pkta2D(klon,klev  ),tair2D(klon,klev  ),
     .                rolv2D(klon,klev  ), pst2D(klon       ),
     .                gpmi2D(klon,klev+1),gplv2D(klon,klev+1),
     .                uair2D(klon,klev  ),vair2D(klon,klev  ),
     .                  qv2D(klon,klev  ),  pk2D(klon,klev  ),
     .                pst2Dn(klon       ),wair2D(klon,klev  )

c      common /hy_dyn/ pkta2D             ,tair2D             ,
c     .                rolv2D             , pst2D             ,
c     .                gpmi2D             ,gplv2D             ,
c     .                uair2D             ,vair2D             ,
c     .                  qv2D             ,  pk2D             ,
c     .                pst2Dn             ,wair2D


      real              qg2D(klon,klev  ),  qw2D(klon,klev  )
     .              ,   qr2D(klon,klev  ),  qi2D(klon,klev  )
     .              ,   qs2D(klon,klev  ),cfra2D(klon,klev  )
     .              ,  dqi2D(klon,klev  ), dqw2D(klon,klev  )
     .              , qvsw2D(klon,klev+1),qvsi2D(klon,klev+1)
     .              , ccni2D(klon,klev  ),ccnw2D(klon,klev  )
     .              , rain2D(klon       ),snow2D(klon       )
     .              , crys2D(klon       ),hlat2D(klon,klev  )
c #BS.              ,   vs2D(klon,klev  )

c      common /hy_hyd/   qg2D             ,  qw2D
c     .              ,   qr2D             ,  qi2D
c     .              ,   qs2D             ,cfra2D
c     .              ,  dqi2D             , dqw2D
c     .              , qvsw2D             ,qvsi2D
c     .              , ccni2D             ,ccnw2D
c     .              , rain2D             ,snow2D
c     .              , crys2D             ,hlat2D
c #BS.              ,   vs2D


      real            ect_2D(klon,klev),TUkv2D(klon,klev),
     .                zi__2D(klon)

c      common /hy_ect/ ect_2D           ,TUkv2D           ,
c     .                zi__2D      


      real            prec2D(klon),snoh2D(klon),
     .                tsrf2D(klon)

c      common /hy_sfl/ prec2D      ,snoh2D      ,
c     .                tsrf2D      


      real            W2xyz1(klon,klev  ),W2xyz2(klon,klev  ),
     .                W2xyz3(klon,klev  ),W2xyz4(klon,klev  ),
     .                W2xyz5(klon,klev+1),W2xyz6(klon,klev+1),
     .                W2xyz7(klon,klev+1),W2xyz8(klon,klev+1),
     .                W2xyz9(klon,klev  ),W2xyz0(klon,klev  )

c      common /hy_wrk/ W2xyz1             ,W2xyz2             ,
c     .                W2xyz3             ,W2xyz4             ,
c     .                W2xyz5             ,W2xyz6             ,
c     .                W2xyz7             ,W2xyz8             ,
c     .                W2xyz9             ,W2xyz0             


      real            wat01D(klon),wat11D(klon),
     .                wat21D(klon),watf1D(klon),
     .                enr01D(klon),enr11D(klon),
     .                enr21D(klon)

c      common /hy_ew1/ wat01D      ,wat11D      ,
c     .                wat21D      ,watf1D      ,
c     .                enr01D      ,enr11D      ,
c     .                enr21D


      character*20    mphy2D(klon)

c      common /hy_ew2/ mphy2D


      integer         jhlr2D(klon)

c      common /hy__ge/ jhlr2D


      integer         ioutIO(5)

      common /hy__io/ ioutIO
