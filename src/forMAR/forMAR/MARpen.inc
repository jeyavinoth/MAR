      real             akhm,
     .               pente1          ,pente2          ,pente3
      common/MARpen/   akhm(mx,my,mz),
     .               pente1(mx,my,mz),pente2(mx,my,mz),pente3(mx,my,mz)
C +...               pente1, pente2 and pente3 are used in the correction
C +                  towards horizontal diffusion (hdiff and vdiff)
C +
      real           slopex,
     .               slopey
      common/MARslo/ slopex(mx,my,mzz),
     .               slopey(mx,my,mzz)
C +...               slopex and slopey are the slope of the sigma surfaces
C +                  in the x and y direction respectively
C +
