      integer   mx    ,my   ,ip11  ,jp11
      parameter(mx=  80,my= 140,ip11=  2,jp11=  2)
      integer   mz   ,mzir1     ,mzir
      parameter(mz=  25,mzir1=mz+1,mzir=mz+2)
c ... mzir1 may be chosen much larger than mz, 
c     if the model vertical domain covers a small
c     part of the air column
c     
      integer   mx1     ,mx2
      parameter(mx1=  79,mx2=  78)
      integer   my1     ,my2     ,myd2
      parameter(my1= 139,my2= 138,myd2=1+my/2)
      integer   mz1     ,mzz
      parameter(mz1=mz-1,mzz=mz+1)
      integer   i_2
      parameter(i_2=mx-mx1+1) 
      integer   j_2
      parameter(j_2=my-my1+1) 
      integer   mzabso      ,mzhyd
      parameter(mzabso =   6,mzhyd=mzabso+1)
c     
      integer   klon,        klev
      parameter(klon=      1,klev=mz)
C +...if #NV removed (NO vectorization)
C +   then      klon=      1
C +   
      integer   kdlon,       kflev
      parameter(kdlon=klon  ,kflev=klev)
C +   
      integer   n6   ,n7
      parameter(n6= 6,n7= 7)
C +.. n6 et n7 determine a relaxation zonetowards lateral boundaries
C +   (large scale values of the variables).
C +   This zone extends over n6-1 points.
C +   Davies (1976) propose 5 points (i.e. n6=6 and n7=7)
C +   
      integer   mw
      parameter(mw=  2)
C +..           mw is the total number of mosaics
C +   
