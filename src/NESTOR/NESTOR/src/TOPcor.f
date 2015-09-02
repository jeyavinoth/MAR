C   +-------------------------------------------------------------------+
C   |  Subroutine TOPcor                               19/07/04 NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input  : NST_sh : Topography prescribed from data sets            |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output : NST_sh : Corrected topography                            |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Options: values given to TOPopt :                                 |
C   | ^^^^^^^  1 = Border of constant NST topography at boundaries      |
C   |          2 = Imposed LSC topography in the const. border at bound.|
C   |          3 = Imposed LSC topography in the whole domain           |
C   |          4 = Zero topography in the constant border               |
C   |          5 = Topography filtering (2D and 3D)                     |
C   |          Note that these options can be combined.                 |
C   |                                                                   |
C   |                                                                   |
C   | Explanations of boundary structure : see parameters in NSTdim.inc |
C   | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                               |
C   |                                                                   |
C   | 1. TOPOGRAPHY                                                     |
C   | -------------                                                     |
C   |  Constant  | Transition |  Computation  | Transition |  Constant  |
C   | topography |    zone    |    domain     |    zone    | topography |
C   |    zone    | (LS -> MAR)|               | (LS -> MAR)|    zone    |
C   ^            ^            ^               ^            ^            ^
C   1   ...     n10  ...  n10+n8+1  ...  mx-n9-n8-1 ...  mx-n9   ...   mx
C   |                                                                   |
C   | 2. RELAXATION LSC --> NST                                         |
C   | -------------------------                                         |
C   |     Relaxation     |      Computation      |      Relaxation      |
C   |        zone        |        domain         |         zone         |
C   ^                    ^                       ^                      ^
C   1        ...        n7         ...         mx-n6        ...        mx
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE TOPcor (TOPopt)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'INTvar.inc'

      INTEGER i,j,mmx,mmy,n88,n9x,n9y,n10x,n10y,TOPopt,ind,ii2,jj2

      REAL    TMP_sh (mx,my),aux  


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Define temporary variables related to array sizes
C +   -------------------------------------------------

      n9x  = MIN(n9,mx-1)
      n9y  = MIN(n9,my-1)
      n10x = MIN(n10,mx)
      n10y = MIN(n10,my)

      mmx  = mx
      mmy  = my

      IF (mmx.eq.1) THEN
       n9x  = 0
       n10x = 1
      ENDIF

      IF (mmy.eq.1) THEN
       n9y  = 0
       n10y = 1
      ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---OPTION 1
C +   ********


C +---Border of constant NST topography at boundaries
C     ===============================================

C +...Topography in the relaxation zone of the mesoscale domain 
C +...is given by topography from data sets and has constant 
C +...value in this region.

      IF (TOPopt.eq.1) THEN

       WRITE(6,*) 'Border of constant NST topography'
       WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       WRITE(6,*)

       DO j=1,my
        DO i=1,n9x
         NST_sh(   i  ,j) = NST_sh(   n10x,j)
         NST_sh(mx-i+1,j) = NST_sh(mx- n9x,j)
        END DO
       END DO

       DO j=1,n9y  
        DO i=1,mx
         NST_sh(i,   j  ) = NST_sh(i,  n10y)
         NST_sh(i,my-j+1) = NST_sh(i,my-n9y)
        END DO
       END DO

      ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---OPTION 2
C +   ********


C +---Imposed LSC topography in the relaxation zone at the boundaries
C +   ===============================================================

C +...Topography in the relaxation zone of the mesoscale domain is
C +...given by the LSC topography. A transition zone (size=n8)
C +...avoids strong changes of the surface elevation.

      IF (TOPopt.eq.2) THEN

       WRITE(6,*) 'Imposed LSC topography in the relax. zone'
       WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       WRITE(6,*)


C +---Transition area
C +   ---------------

       n88=MAX(1,n8)+1

       IF (mmx.eq.1.or.mmy.eq.1) n88=0


C +---Bottom and top boundaries (corners included)
C +   --------------------------------------------
 
       DO i=n10x,mx-n9x

        DO j=1,n10y
         NST_sh(i,j) = INT_sh(i,n10y)
        ENDDO

        DO j=my-n9y,my
         NST_sh(i,j) = INT_sh(i,my-n9y)
        ENDDO

       ENDDO


C +---Left and right boundaries (corners not included)
C +   ------------------------------------------------

       DO j=1,my
     
        ind=min(j,my-n9y)
        ind=max(ind,n10y)

        DO i=1,n10x
         NST_sh(i,j) = INT_sh(n10x,ind)
        ENDDO

        DO i=mx-n9x,mx
         NST_sh(i,j) = INT_sh(mx-n9x,ind)
        ENDDO

       ENDDO


C +---Treatment of transition area between large-scale
C +   topography and mesoscale topography ------------
C +   -----------------------------------

       DO i=n10x+1,mx-n9x-1

        ind = min(i,mx-n9x-n88)
        ind = max(ind,n10x+n88 )

        DO j=n10y+1,n10y+n88-1
         aux = real(j-n10y)/real(n88)
         IF (i.ge.j.and.(mx-i).ge.j)
     .    NST_sh(i,j) = (1.-aux)*NST_sh(i,n10y)
     .                     +aux *NST_sh(ind,n10y+n88)
        ENDDO

        DO j=my-n9y-n88+1,my-n9y-1
         aux = real(j-(my-n9y-n88))/real(n88)
         IF (i.le.j.and.(mx-i).le.j)
     .    NST_sh(i,j) = (1.-aux)*NST_sh(i,my-n9y-n88)
     .                      +aux*NST_sh(ind,my-n9y)
        ENDDO

       ENDDO


       DO j=n10y+1,my-n9y-1

        ind = min(j,my-n9y-n88)
        ind = max(ind,n10y+n88)

        DO i=n10x+1,n10x+n88-1
         aux = real(i-n10x)/real(n88)
         IF (i.le.j.and.(mx-i).ge.j)
     .    NST_sh(i,j) = (1.-aux)*NST_sh(n10x,j)
     .                     +aux *NST_sh(n10x+n88,ind)
        ENDDO

        DO i=mx-n9x-n88+1,mx-n9x-1
         aux = real(i-(mx-n9x-n88))/real(n88)
         IF (i.ge.j.and.(mx-i).le.j)
     .    NST_sh(i,j) = (1.-aux)*NST_sh(mx-n9x-n88,ind)
     .                      +aux*NST_sh(mx-n9x,j)
        ENDDO

       ENDDO

      ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---OPTION 3
C +   ********


C +---Imposed LSC topography in the whole domain
C +   ==========================================

C +...LSC topography is imposed for each NST grid point.

      IF (TOPopt.eq.3) THEN

       WRITE(6,*) 'Imposed LSC topography in the whole domain'
       WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       WRITE(6,*)

C +---NST topography = LSC topography
C +   -------------------------------
 
       DO j=1,my
       DO i=1,mx
        NST_sh(i,j) = INT_sh(i,j)
       ENDDO
       ENDDO

      ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---OPTION 4
C +   ********


C +---Zero topography in the relaxation zone
C +   ======================================

C +...Topography in the relaxation zone of the NST domain is set to
C +...zero (mean sea level). A transition zone (size=n8) avoids strong
C +...changes between the relaxation area and the computation domain.

      IF (TOPopt.eq.4) THEN

       WRITE(6,*) 'Zero topography in the relaxation zone'
       WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       WRITE(6,*)

C +---Transition area
C +   ---------------

       n88=MAX(1,n8)+1

C +---Bottom and top boundaries (corners included)
C +   --------------------------------------------
 
       DO i=n10x,mx-n9x

        DO j=1,n10y
         NST_sh(i,j) = 0.
        ENDDO

        DO j=my-n9y,my
         NST_sh(i,j) = 0.
        ENDDO

       ENDDO


C +---Left and right boundaries (corners not included)
C +   ------------------------------------------------

       DO j=1,my
     
        ind=min(j,my-n9y)
        ind=max(ind,n10y)

        DO i=1,n10x
         NST_sh(i,j) = 0.
        ENDDO

        DO i=mx-n9x,mx
         NST_sh(i,j) = 0.
        ENDDO

       ENDDO


C +---Treatment of transition area between large-scale
C +   topography and mesoscale topography ------------
C +   -----------------------------------

       DO i=n10x+1,mx-n9x-1

        ind = min(i,mx-n9x-n88)
        ind = max(ind,n10x+n88 )

        DO j=n10y+1,n10y+n88-1
         aux = real(j-n10y)/real(n88)
         IF (i.ge.j.and.(mx-i).ge.j)
     .    NST_sh(i,j) = (1.-aux)*NST_sh(i,n10y)
     .                     +aux *NST_sh(ind,n10y+n88)
        ENDDO

        DO j=my-n9y-n88+1,my-n9y-1
         aux = real(j-(my-n9y-n88))/real(n88)
         IF (i.le.j.and.(mx-i).le.j)
     .    NST_sh(i,j) = (1.-aux)*NST_sh(i,my-n9y-n88)
     .                      +aux*NST_sh(ind,my-n9y)
        ENDDO

       ENDDO


       DO j=n10y+1,my-n9y-1

        ind = min(j,my-n9y-n88)
        ind = max(ind,n10y+n88)

        DO i=n10x+1,n10x+n88-1
         aux = real(i-n10x)/real(n88)
         IF (i.le.j.and.(mx-i).ge.j)
     .    NST_sh(i,j) = (1.-aux)*NST_sh(n10x,j)
     .                     +aux *NST_sh(n10x+n88,ind)
        ENDDO

        DO i=mx-n9x-n88+1,mx-n9x-1
         aux = real(i-(mx-n9x-n88))/real(n88)
         IF (i.ge.j.and.(mx-i).le.j)
     .    NST_sh(i,j) = (1.-aux)*NST_sh(mx-n9x-n88,ind)
     .                      +aux*NST_sh(mx-n9x,j)
        ENDDO

       ENDDO

      ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---OPTION 5
C +   ********


C +---Topography filtering (2D and 3D)
C +   ================================

      IF (TOPopt.eq.5) THEN

       WRITE(6,*) 'Topography filtering'
       WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~'
       WRITE(6,*)

       IF (mmx.gt.1) THEN

        IF (mmy.gt.1) THEN

         DO j = 2,mmy-1
         DO i = 2,mmx-1
         TMP_sh(i,j)=NST_sh(i-1,j+1)+2.*NST_sh(i,j+1)+   NST_sh(i+1,j+1)
     .           +2.*NST_sh(i-1,j  )+4.*NST_sh(i,j  )+2.*NST_sh(i+1,j  )
     .           +   NST_sh(i-1,j-1)+2.*NST_sh(i,j-1)+   NST_sh(i+1,j-1)
         ENDDO
         ENDDO

         DO j = 2,mmy-1
         DO i = 2,mmx-1
          NST_sh (i,j) = TMP_sh(i,j) / 16.
         ENDDO
         ENDDO

         jj2 = 2
         DO i=1,mx
          NST_sh (i, 1) = NST_sh(i,  jj2)
          NST_sh (i,my) = NST_sh(i,mmy-1)
         ENDDO

         ii2 = 2
         DO j=1,my
          NST_sh ( 1,j) = NST_sh(  ii2,j)
          NST_sh (mx,j) = NST_sh(mmx-1,j)
         ENDDO

        ELSE

         j = 1
 
         DO i=2,mmx-1
          TMP_sh(i,j)=NST_sh(i-1,j)+2.*NST_sh(i,j)+NST_sh(i+1,j)
         ENDDO

         DO i=2,mmx-1
          NST_sh(i,j)=TMP_sh(i,j) / 4.
         ENDDO

        ENDIF

       ENDIF

      ENDIF


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
