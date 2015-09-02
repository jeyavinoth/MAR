C +---------------------------------------------------------------------+
C | MAR (output)                                           03-1999  MAR |
C |   SubRoutine Wind.f computes the wind speed/direction from the      |
C |   the x-wind et y-wind MAR composantes                              |
C +---------------------------------------------------------------------+
C |                                                                     |
C |        dx  = MAR resojjution                                         |
C |        ii  = MAR coord.                                             |
C |        jj  = MAR coord                                              |
C |    MARjjat  = MAR jjat                                                |
C |    MARjjon  = MAR jjon                                                |
C |                                                                     |
C |       uin  = x-wind composante                                      |
C |       vin  = x-wind composante                                      |
C |       dout = wind direction                                         |
C |       sout = wind speed                                             |
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE Wind (uin,vin,ii,jj,dx,MARlat, MARlon,dout,sout)
      
      IMPLICIT NONE
      
C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
            
            
      INTEGER  ii , jj
      REAL     uin, vin, dout, sout
      REAL     uu , vv
      REAL     a  , pi , conv, dx
      REAL     MARlat(mx,my), MARlon(mx,my)
          
      a    = 1000.0 * 6371.229 ! radius of the Earth [m]
      pi   = acos(-1.0)        ! 3.141592654 
      conv = pi/180.0d0        ! conversion
      
      uu    = (MARlon(ii+1,jj) - MARlon(ii,jj))*conv/dx*a*
     .         cos(MARlat(ii,jj)*conv) * uin +
     .        (MARlon(ii,jj+1) - MARlon(ii,jj))*conv/dx*a*
     .         cos(MARlat(ii,jj)*conv) * vin

      vv    = (MARlat(ii+1,jj) - MARlat(ii,jj))*conv/dx*a*
     .         uin +                   
     .        (MARlat(ii,jj+1) - MARlat(ii,jj))*conv/dx*a*
     .         vin       

      dout = 0.0d0

      IF (uu.gt.0.0.and.vv.ge.0.0)
     .          dout = 3.0*pi/2.0-atan(vv/uu)
      IF (uu.lt.0.0.and.vv.ge.0.0)  
     .          dout = 5.0*Pi/2.0-atan(vv/uu)
      IF (uu.lt.0.0.and.vv.le.0.0)  
     .          dout = 5.0*Pi/2.0-atan(vv/uu)
      IF (uu.gt.0.0.and.vv.le.0.0) 
     .          dout = 3.0*Pi/2.0-atan(vv/uu)
      IF (uu.eq.0.0.and.vv.ge.0.0)
     .          dout = Pi
      IF (uu.eq.0.0.and.vv.le.0.0) 
     .          dout = 0.0
      IF (dout.gt.2*Pi)         
     .          dout = dout - 2*Pi

	 dout = dout/ conv
      sout = sqrt(uin*uin+vin*vin)
              
      RETURN
      END 