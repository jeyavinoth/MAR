C   +-------------------------------------------------------------------+
C   |  Subroutine INTnrg2                          31/08/2004   NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | This routine is a linear interpolation of a 2D scalar fields from |
C   | a non-regular grid to a non-regular grid.                         |
C   |                                                                   |
C   | Input : grd_Ix (ni, nj) : Input grid points position x(i,j)       |
C   | ^^^^^^^ grd_Iy (ni, nj) :   "     "    "       "     y(i,j)       |
C   |         var_I  (ni, nj) : Input field values                      |
C   |         grd_Ox (mx, my) : Output grid positions x(i,j)            |
C   |         grd_Oy (mx, my) : Output grid positions y(i,j)            |
C   |                                                                   |
C   | Output: var_O  (mx, my) : Output field values                     |
C   | ^^^^^^^                                                           |
C   +-------------------------------------------------------------------+


      SUBROUTINE INTnrg2 (ni2,nj2,grd_Ix,grd_Iy,var_I,
     .                    mx2,my2,grd_Ox,grd_Oy,var_O,
     .                    pos_Ox,pos_Oy)


      IMPLICIT NONE

      include "NSTdim.inc"
      include "NESTOR.inc"

C +---General and local variables
C +   ---------------------------

      INTEGER  i,j,k,l,ni2,nj2,mx2,my2,k1,k2,l1,l2,ii,jj

      INTEGER  pos_Ox(mx,my),pos_Oy(mx,my)

      INTEGER  ii_min(mx,my),ii_max(mx,my) 
      INTEGER  jj_min(mx,my),jj_max(mx,my)  

      REAL grd_Ix(ni,nj),grd_Iy(ni,nj),grd_Ox(mx,my),grd_Oy(mx,my),
     .     var_I(ni,nj) ,var_O (mx,my),int_O (mx,my),nbr_meshes,epsi,
     .     delta_lat(mx,my),delta_lon(mx,my),xx

      DATA epsi / 0.1 /

      common/INTnrg2_I/ii_min,ii_max,jj_min,jj_max
      common/INTnrg2_I/delta_lat,delta_lon

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(I_time.le.1.and.pos_Ox(mx,my).eq.0)then

      ii_min = ni ; ii_max = 1
      jj_min = nj ; jj_max = 1

      DO k=1,mx 
      DO l=1,my 

       k1=max(1,min(mx,k+1))
       k2=max(1,min(mx,k-1))

       l1=max(1,min(my,l+1))
       l2=max(1,min(my,l-1))

       delta_lat(k,l)=max(abs(grd_Iy(k,l1)-grd_Iy(k,l)),
     .                    abs(grd_Iy(k,l2)-grd_Iy(k,l)),
     .                    abs(grd_Iy(k1,l)-grd_Iy(k,l)),
     .                    abs(grd_Iy(k2,l)-grd_Iy(k,l)))

       delta_lon(k,l)=max(abs(grd_Iy(k,l1)-grd_Iy(k,l)),
     .                    abs(grd_Iy(k,l2)-grd_Iy(k,l)),
     .                    abs(grd_Iy(k1,l)-grd_Iy(k,l)),
     .                    abs(grd_Iy(k2,l)-grd_Iy(k,l)))

       int_O(k,l)=0.
       
       xx=0.

       do while (int_O(k,l).le.0.)

        delta_lat(k,l)=delta_lat(k,l)*(1+xx)
        delta_lon(k,l)=delta_lon(k,l)*(1+xx)
        xx=xx+0.05

        do i = 1, ni 
        do j = 1, nj
         IF(abs(grd_Ox(k,l)-grd_Ix(i,j)).le.delta_lon(k,l).and. 
     .      abs(grd_Oy(k,l)-grd_Iy(i,j)).le.delta_lat(k,l))then
           int_O(k,l) = int_O(k,l)+1
          jj_min(k,l) = min(j,jj_min(k,l))
          ii_min(k,l) = min(i,ii_min(k,l))  
          jj_max(k,l) = max(j,jj_max(k,l))
          ii_max(k,l) = max(i,ii_max(k,l)) 
         ENDIF
        end do
        end do

        if(xx.eq.5) then
        WRITE(6,*) 'No cell of input grid includes an output grid'
        WRITE(6,*) 'point.                     --- STOP in INTnrg.'
        stop
        endif

       enddo

       pos_Ox(k,l)=jj_min(k,l)

      ENDDO
      ENDDO

      endif          

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO k=1,mx 
      DO l=1,my 

       int_O(k,l)=0.
       var_O(k,l)=0.

       do i = ii_min(k,l),ii_max(k,l)
       do j = jj_min(k,l),jj_max(k,l) 
    
        IF(abs(grd_Ox(k,l)-grd_Ix(i,j)).le.delta_lon(k,l).and. 
     .      abs(grd_Oy(k,l)-grd_Iy(i,j)).le.delta_lat(k,l))then
           int_O(k,l) = int_O(k,l)+1
           var_O(k,l) = var_O(k,l)+var_I(i,j)
         ENDIF
       
       end do
       end do
  
       var_O(k,l)=var_O(k,l)/int_O(k,l)

      ENDDO
      ENDDO

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      RETURN
      END
