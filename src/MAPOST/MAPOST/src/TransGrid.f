C   +-------------------------------------------------------------------+
C   | Subroutine TransGrid                                     1999.03  |
C   |                                                          v1.0a2   |
C   +-------------------------------------------------------------------+
C   |  * "Transfer" data with approx. "flux" conservation from          |
C   |     the MAR grid to the Large Scale (LS) Lat/Lon grid (rectang.). |
C   |     Usage: get precip. for a (usually more coarse) Lat/lon grid   |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   |  INPUT :                                                          |
C   |  ^^^^^^^ MARval  (mx, my) : Input MAR values                      |
C   |          MARlon  (mx, my) : MAR (input) grid positions lon(i,j)   |
C   |          MARlat  (mx, my) : MAR (input) grid positions lat(i,j)   | 
C   |          MARx    (mx)     : MAR (input) grid positions x(i)       |
C   |          MARy    (my)     : MAR (input) grid positions x(i)       |
C   |          LSlon   (LSni)   : Output grid positions lon(i)          |
C   |          LSlat   (LSnj)   :   "     "    "        lat(j)          |
C   |                                                                   |
C   |  OUTPUT: LSval(LSni, LSnj): Output field values                   |
C   |  ^^^^^^^                                                          |
C   +-------------------------------------------------------------------+
      SUBROUTINE TransGrid  (MAPtyp,
     &       MARlon, MARlat, MARx, MARy, MARval,
     &       LSlon, LSlat, LSval)

      IMPLICIT NONE

C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

C +---Physics (smaller package with defined variables.)
C +   -------------------------------------------------
      include 'LSMphy.inc'

C +---MAPOST parameters & variables
C +   -----------------------------
      include 'MAPOST.inc'
      include 'globals.inc'
      

C     ** input
      INTEGER MAPtyp 
      REAL LSlon (LSni), LSlat (LSnj)
      REAL MARlon (mx, my), MARlat (mx, my)
      REAL MARx(mx), MARy(my), MARval (mx, my)

C     ** output
      REAL LSval  (LSni, LSnj)

C     ** local:
      INTEGER BSni, BSnj
      PARAMETER (BSni=LSni+1, BSnj=LSnj+1)
C     ^^Number of meshes boundaries.
      REAL BLSlat (BSnj)
      REAL BLSlon (BSni)
      REAL tweig  (LSni, LSnj)
      REAL GElon0, GElat0, SGlon, SGlat
      REAL MARxm, MARym, dy
      REAL MARxL, MARyL, MARxS, MARyS,SGxx,SGyy

      INTEGER ii,jj,iLS,jLS,ix0,jy0,iSi,jSj,nSn,imexcl
      REAL rSn
      
      CALL DInfo(2,'TransGrid: Begin')
       
C +---Set coordinate for the mesh boundaries.
C +   (this routine is valid for both ascending
C +    and descending indexes / data).
C +   !See info/example below to understand !
C +   -----------------------------------------

      DO ii = 2,LSni
        BLSlon(ii) = (LSlon(ii-1)+LSlon(ii))/2.
      ENDDO
      BLSlon(1)    = BLSlon(2)     *2.-BLSlon(3)
      BLSlon(BSni) = BLSlon(BSni-1)*2.-BLSlon(BSni-2)

      DO jj = 2,LSnj
        BLSlat(jj) = (LSlat(jj-1)+LSlat(jj))/2.
      ENDDO
      BLSlat(1)    = BLSlat(2)     *2.-BLSlat(3)
      BLSlat(BSnj) = BLSlat(BSnj-1)*2.-BLSlat(BSnj-2)

C +---Initialisation
C +   --------------
      CALL DInfo(2,'Initialisation')

      DO jLS = 1,LSnj
      DO iLS = 1,LSni
        LSval (iLS,jLS) = 0 
        tweig (iLS,jLS) = 0
      ENDDO
      ENDDO

C +---Number of "integration" Sub-Grid points (/direction)
C +   - - - - - - - - - - - - - - - - - - - - - - - - - - -
      nSn = 7
      rSn = FLOAT(nSn)

C +---Get data for the geographical projection MAR->sphere
C +   - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL DInfo(2,'Projection') 
C
C +-  Find the center of the MAR grid in lat/lon => project:
C
      CALL SEARCH (MARx, mx, 0.0, ix0)
      CALL SEARCH (MARy, my, 0.0, jy0)
      GElon0= MARlon(ix0, jy0)
      GElat0= MARlat(ix0, jy0)
C +   ^(both in degrees = ok for GRDxxx routine)

C +-  Check that we use the correct projection and parameters:
C     !Note that MARx is in km while GRDxxx use meters !
C 
      MARxm = MARx(mx) * 1.E3
      MARym = MARy(my) * 1.E3
      dy = (MARy(2)-MARy(1)) * 1.E3
      IF (MAPtyp.LE.1) THEN
        CALL GRDstr(MARxm,MARym,GElon0,GElat0,SGlon,SGlat)
C            ^^^^^^
      ELSE
        CALL GRDlam(MARxm,MARym,dy,GElon0,GElat0,SGlon,SGlat)
C            ^^^^^^
      ENDIF
      SGlon= SGlon * 15. !(360/24 = hour to deg)
      SGlat= SGlat / degrad
cXF
      IF (ABS(SGlon-MARlon(mx,my)) .GT. 0.01 .OR.
     .    ABS(SGlat-MARlat(mx,my)) .GT. 0.01)THEN
       !WRITE(*,*) 'TransGrid: Error - invalid Geo Proj.'
       !WRITE(*,*) 'MARlon/lat:',MARlon(mx,my),MARlat(mx,my)
       !WRITE(*,*) 'This Proj :',SGlon        ,SGlat      
       !WRITE(*,*) 'Grid Cent :',ix0,jy0
       !STOP
      ENDIF
 
C +---Add values on the LS grid: Main loop over MAR grid.
C +   ---------------------------------------------------
C     (This routine can not handle the boundary row,
C     and, moreover, near boundary results are not relevant.
C     Therefore, "imexcl" rows are excluded).
C
      CALL DInfo(2,'Main Loop...')

      imexcl= 3
      DO jj = 1+imexcl,my-imexcl 
      DO ii = 1+imexcl,mx-imexcl
 
C +-     -*Size and bottom-left corner of MAR mesh:
         MARxL  = (MARx(ii-1)+MARx(ii)) /2.0
         MARyL  = (MARy(jj-1)+MARy(jj)) /2.0
         MARxS  = (MARx(ii+1)-MARx(ii-1))/ 2.0
         MARyS  = (MARy(jj+1)-MARy(jj-1))/ 2.0

C +-     -*Begin Loop over "integration" sub-grid (SG): 
         DO jSj = 1, nSn
         DO iSi = 1, nSn 

C +-       -*MAR coordinate of the current SG point: 
C           !Note that MARx is in km while GRDxxx use meters !

           SGxx = (MARxL + MARxS * (FLOAT(iSi)-0.5)/rSn )*1.E3
           SGyy = (MARyL + MARyS * (FLOAT(jSj)-0.5)/rSn )*1.E3

C +-       -*Spherical coordinate of this SG element:
           IF (MAPtyp.LE.1) THEN
             CALL GRDstr(SGxx,SGyy,GElon0,GElat0,SGlon,SGlat)
C                 ^^^^^^
           ELSE
             CALL GRDlam(SGxx,SGyy,dy,GElon0,GElat0,SGlon,SGlat)
C                 ^^^^^^
           ENDIF
           SGlon= SGlon * 15. !(360/24 = hour to deg)
           SGlat= SGlat / degrad 

C +-       -*Search for the LS mesh which includes this SG:
           CALL SEARCH(BLSlon,BSni,SGlon,iLS)
           CALL SEARCH(BLSlat,BSnj,SGlat,jLS)
C
C +-       -*Increment the out variable and the weight:
           LSval(iLS,jLS) = LSval(iLS,jLS) + MARval(ii,jj)
           tweig(iLS,jLS) = tweig(iLS,jLS) + 1.0

C +        --Info - example:
C 
C          BLSlat(jj)=    61        60        59     ...
C              jj    =     1         2         3
C             jLS    =     |    1    |    2    |  ... 
C          (eg. for a data at lat=60.3, 
C           SEARCH returns the "LOWEST" index i.e. jj=1; 
C           Therefore, jLS=1 and the data is attributed to
C           LSval(...,1)  )
C    
           
         ENDDO
         ENDDO

      ENDDO
      ENDDO
C 
C +---End computation: divide results / weight table.
C +   -----------------------------------------------
C     NB: LS mesh which does not include at least 1/2 mar
C         mesh with valid results are excluded !!
C
      CALL DInfo(2,'End computation')

      DO jLS = 1,LSnj
      DO iLS = 1,LSni
        IF (tweig(iLS,jLS).GT.(rSn*rSn/2.)) THEN
          LSval(iLS,jLS)= LSval (iLS,jLS) / tweig (iLS,jLS) 
        ELSE 
          LSval(iLS,jLS)= 1.E30
        ENDIF
      ENDDO
      ENDDO

      CALL DInfo(2,'TransGrid: End.')
      
      END
