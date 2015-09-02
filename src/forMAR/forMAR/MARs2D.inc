C +
      real    snob2D,hwat2D,hsno2D,hsen2D,hbal2D,tsub2D,
     .        t2_2D ,dtg2D ,fird2D,firm2D,eps02D,albe2D,
     .        sol_2D,wg_2D ,wk_2D ,w2_2D ,wx_2D ,ssv2D ,
     .        utsl2D,uusl2D,uqsl2D,cdh2D ,ch02D ,raer2D,
     .        rsur2D,r0_2D ,z0_2D ,d1_2D ,qvap2D,lmol2D,
     .        hlsf2D,rose2D
C +
      common /hy_sfl/ snob2D(klon   ),hwat2D(klon     ),
     .                hsno2D(klon   ),hsen2D(klon     ),
     .                hbal2D(klon   ),tsub2D(klon,mw  ),
     .                t2_2D (klon   ),dtg2D (klon,mw  ),
     .                fird2D(klon   ),firm2D(klon     ),
     .                eps02D(klon   ),albe2D(klon     ),
     .                sol_2D(klon   ),wg_2D (klon     ),
     .                wk_2D (klon   ),w2_2D (klon     ),
     .                wx_2D (klon   ),ssv2D (klon,klev),
     .                utsl2D(klon,mw),uusl2D(klon,mw  ),
     .                uqsl2D(klon,mw),cdh2D (klon,mw  ),
     .                ch02D (klon   ),raer2D(klon     ),
     .                rsur2D(klon   ),r0_2D (klon,mw  ),
     .                z0_2D (klon,mw),d1_2D (klon     ),
     .                qvap2D(klon   ),lmol2D(klon,mw  ),
     .                hlsf2D(klon   ),rose2D(klon     )
C +
C +
      real    RAds2D,RAdi2D
C +
      common /hy_rad/ RAds2D(klon),RAdi2D(klon)
C +
