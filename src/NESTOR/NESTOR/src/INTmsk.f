C   +-------------------------------------------------------------------+
C   |  Subroutine INTmsk                            01-07-2004  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Extrapolation of sea large-scale data to land to reduce the       |
C   | problem of the fjord.
C   |                                                                   |
C   +-------------------------------------------------------------------+

      SUBROUTINE INTmsk(var1)

 
      IMPLICIT NONE

C +---Include files
C +   -------------
      
      include "NSTdim.inc"

C +---Local variables
C +   ---------------      
      
      integer i,j,k,l,tt
      real    var1(ni,nj),var2(ni,nj)
      real    nbr,valmax
      
C +---Variable max value
C +   -------------------
      
      valmax=9.9e20
      
C +---Extrapolation to land
C +   ---------------------

      do i=1,ni
      do j=1,nj
       var2(i,j)=var1(i,j)
      enddo
      enddo
     
      do i=2,ni-1
      do j=2,nj-1      
      
       nbr=0.
      
       if(var1(i,j).ge.valmax) then
             
        do k=-1,1,1
        do l=-1,1,1
         if(var1(i+k,j+l).le.valmax) then       
          if(var2(i,j).ge.valmax) var2(i,j)=0.0
          var2(i,j)=var2(i,j)+var1(i+k,j+l)
          nbr=nbr+1
         endif
        enddo
        enddo

        var2(i,j)=var2(i,j)/max(1.,nbr)
       endif

      enddo
      enddo


      do i=1,ni
      do j=1,nj
       var1(i,j)=var2(i,j)
      enddo
      enddo


      do tt=1,5

       do i=1,ni
       do j=1,nj
        var2(i,j)=var1(i,j)
       enddo
       enddo

       do i=2,ni-1
       do j=2,nj-1     
      
        nbr=0.

        if(var1(i,j).ge.valmax) then
             
         do k=-1,1,1
         do l=-1,1,1
          if(var2(i+k,j+l).le.valmax) then
           if(k.ne.0.or.l.ne.0) then
           if(var1(i,j).ge.valmax) var1(i,j)=0.0
           var1(i,j)=var1(i,j)+var2(i+k,j+l)
           nbr=nbr+1
           endif
          endif
         enddo
         enddo

         var1(i,j)=var1(i,j)/max(1.,nbr)

        endif
            
       enddo
       enddo

      enddo

      do i=1,ni
       if(var1(i,1) .ge.valmax) var1(i,1) =var1(i,2)
       if(var1(i,nj).ge.valmax) var1(i,nj)=var1(i,nj-1)
      enddo

      do j=1,nj
       if(var1(1,j) .ge.valmax) var1(1,j) =var1(2,j)
       if(var1(ni,j).ge.valmax) var1(ni,j)=var1(ni-1,j)
      enddo
      
      END
      
      
