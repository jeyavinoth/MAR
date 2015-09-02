
C +  ==========
C +--MARqqm.inc
C +  ==========

      integer   lb
      integer   lgx ,ldx ,lgy ,ldy
      parameter(lb=1)
      parameter(lgx = min(lb+1,mx     ))
      parameter(ldx = max(   1,mx-lb-2))
      parameter(lgy = min(lb+1,my     ))
      parameter(ldy = max(   1,my-lb-2))

      integer   lgx1 ,ldx1 ,lgy1 ,ldy1
      integer   lgx2 ,ldx2 ,lgy2 ,ldy2
      parameter(lgx1= max(1,lgx-1   ))
      parameter(lgx2= max(1,lgx-2   ))
      parameter(ldx1= min(  ldx+1,mx))
      parameter(ldx2= min(  ldx+2,mx))
      parameter(lgy1= max(1,lgy-1   ))
      parameter(lgy2= max(1,lgy-2   ))
      parameter(ldy1= min(  ldy+1,my))
      parameter(ldy2= min(  ldy+2,my))
      integer   lgxx ,ldxx ,lgyy ,ldyy
      parameter(lgxx= min(  lgx+1,mx))
      parameter(ldxx= max(1,ldx-1   ))
      parameter(lgyy= min(  lgy+1,my))
      parameter(ldyy= max(1,ldy-1   ))
      real      f2_3
      parameter(f2_3= 2./3.)
