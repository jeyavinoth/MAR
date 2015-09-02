C +---Number of sub-regions:
C +   ~~~~~~~~~~~~~~~~~~~~~~
      integer nreg
      parameter (nreg=14)
      
C +---Number of observation stations:
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nsta
      parameter (nsta=12)

C +---Max distance to boundary :
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer ndb
      parameter (ndb= 0.5*min(mx,my)-2)

C +---Definition of pressure levels (! kPa) :
C +   (where output is calculated)
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer npl
      parameter (npl= 16 )
      real plevel(npl)
      data plevel /85., 75., 70., 65., 60., 55., 50., 45., 
     .             40., 35., 30., 25., 20., 15., 10.,  5./

C +---Definition of metre levels (! m) :
C +   (where output is calculated)
C +   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nml
      parameter (nml= 4)
      real mlevel(nml)
      data mlevel /2.,4.,8.,10./ 
