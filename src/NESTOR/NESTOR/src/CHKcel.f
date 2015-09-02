C   +-------------------------------------------------------------------+
C   |  Subroutine CHKcel                             June 2000  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Search for all NST cells contained in a LSC grid cell.            |
C   |                                                                   |
C   | Input : LSC__x (ni, nj) : Input grid points position x(i,j)       |
C   | ^^^^^^^ LSC__y (ni, nj) :   "     "    "       "     y(i,j)       |
C   |         NST__x (mx, my) : Output grid positions x(i,j)            |
C   |         NST__y (mx, my) : Output grid positions y(i,j)            |
C   |         kk,ll           : selection of the LSC cell               |
C   |                                                                   |
C   | Output: icell  (mx, my) : i-index of NST cells in the LSC cell    |
C   | ^^^^^^^ jcell  (mx, my) : j-index of NST cells in the LSC cell    |
C   |         nlist           : number  of NST cells in the LSC cell    |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE CHKcel (ni,nj,LSC__x,LSC__y,mx,my,NST__x,NST__y,
     .                            kk,ll,MXlist,icell,jcell,nlist)


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INTEGER  i,j,k,l,kk,ll,ni,nj,mx,my,ilist,nlist,MXlist

      INTEGER icell(MXlist),jcell(MXlist)

      REAL clon,clat,lat1,lat2,lat3,lat4,lon1,lon2,lon3,lon4,val1,
     .     val2,val3,val4,tmp,dlon12,dlon14,dlat12,dlat14,ilat,ilon

      REAL LSC__x(ni,nj),LSC__y(ni,nj),NST__x(mx,my),NST__y(mx,my)

      LOGICAL includ


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO l=ll,ll    ! LOOP on LSC grid-points : BEGIN
      DO k=kk,kk   
          
C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       lon1=0.5*(LSC__x(k,l)+LSC__x(k-1,l-1))
       lon2=0.5*(LSC__x(k,l)+LSC__x(k+1,l-1))
       lon3=0.5*(LSC__x(k,l)+LSC__x(k+1,l+1))
       lon4=0.5*(LSC__x(k,l)+LSC__x(k-1,l+1))

       lat1=0.5*(LSC__y(k,l)+LSC__y(k-1,l-1))
       lat2=0.5*(LSC__y(k,l)+LSC__y(k+1,l-1))
       lat3=0.5*(LSC__y(k,l)+LSC__y(k+1,l+1))
       lat4=0.5*(LSC__y(k,l)+LSC__y(k-1,l+1))


C +---Rotation of input grid cell ?
C +   =============================

       dlat12=ABS(lat2-lat1)
       dlat14=ABS(lat4-lat1)
       dlon12=ABS(lon2-lon1)
       dlon14=ABS(lon4-lon1)

       IF ((dlat14.lt.dlat12).and.(dlon14.gt.dlon12)) THEN
        tmp =lat1   ! Latitude
        lat1=lat2
        lat2=lat3
        lat3=lat4
        lat4=tmp
        tmp =lon1   ! Longitude
        lon1=lon2
        lon2=lon3
        lon3=lon4
        lon4=tmp
        tmp =val1   ! Values
        val1=val2
        val2=val3
        val3=val4
        val4=tmp
       ENDIF


C +---Invert latitude ?
C +   =================

       IF (lat4.lt.lat1) THEN

        IF (lat3.lt.lat2) THEN
         tmp =lat2   ! Latitude
         lat2=lat3
         lat3=tmp
         tmp =lat1
         lat1=lat4
         lat4=tmp
         tmp =lon2   ! Longitude
         lon2=lon3
         lon3=tmp
         tmp =lon1
         lon1=lon4
         lon4=tmp
         tmp =val2   ! Values
         val2=val3
         val3=tmp
         tmp =val1
         val1=val4
         val4=tmp
        ELSE
         WRITE(6,*) 'Inconsistance in latitude. Cannot be resolved'
         WRITE(6,*) 'by CHKcel subroutine.                --- STOP'
         WRITE(6,*)
         WRITE(6,*) 'Info : ',lat1,lat2,lat3,lat4
         STOP
        ENDIF

       ELSE

        IF (lat3.lt.lat2) THEN

         WRITE(6,*) 'Inconsistance in latitude. Cannot be resolved'
         WRITE(6,*) 'by CHKcel subroutine.                --- STOP'
         WRITE(6,*)
         WRITE(6,*) 'Info : ',lat1,lat2,lat3,lat4
         STOP

        ENDIF

       ENDIF


C +---Invert longitude ?
C +   ==================

       IF (lon2.lt.lon1) THEN

        IF (lon3.lt.lon4) THEN
         tmp =lat3   ! Latitude
         lat3=lat4
         lat4=tmp
         tmp =lat1
         lat1=lat2
         lat2=tmp
         tmp =lon3   ! Longitude
         lon3=lon4
         lon4=tmp
         tmp =lon1
         lon1=lon2
         lon2=tmp
         tmp =val3   ! Values
         val3=val4
         val4=tmp
         tmp =val1
         val1=val2
         val2=tmp
        ELSE
         WRITE(6,*) 'Inconsistance in longitude. Cannot be resolved'
         WRITE(6,*) 'CHKcel subroutine.                    --- STOP'
         WRITE(6,*)
         WRITE(6,*) 'Info : ',lon1,lon2,lon3,lon4
         STOP
        ENDIF

       ELSE

        IF (lon3.lt.lon4) THEN

         WRITE(6,*) 'Inconsistance in longitude. Cannot be resolved'
         WRITE(6,*) 'CHKcel subroutine.                    --- STOP'
         WRITE(6,*)
         WRITE(6,*) 'Info : ',lon1,lon2,lon3,lon4
         STOP

        ENDIF

       ENDIF


C +   At this stage, it is assumed that the input grid cell is
C +   such that :
C +
C +              4----------3
C +              |          |
C +              |          |
C +              |          |
C +              1----------2
C +
C +   with lon1 < lon2, lon4 < lon3
C +        lat1 < lat4, lat2 < lat3


C +---Initialization of list of cells included in LSC cell
C +   ====================================================

       nlist=0
       DO ilist=1,MXlist
        icell(ilist)=0
        jcell(ilist)=0
       ENDDO


C +---Check if input cell includes output grid point
C +   ==============================================


       DO j=1,my    ! LOOP on NST grid-points : BEGIN
       DO i=1,mx

        clon=NST__x(i,j)
        clat=NST__y(i,j)

        includ=.true.

C +---Segment 1-2
C +   -----------

        ilat=lat1+(clon-lon1)/(lon2-lon1)*(lat2-lat1)
        IF (ilat.gt.clat) includ=.false.

C +---Segment 4-3
C +   -----------

        ilat=lat4+(clon-lon4)/(lon3-lon4)*(lat3-lat4)
        IF (ilat.lt.clat) includ=.false.

C +---Segment 2-3
C +   -----------

        ilon=lon2+(clat-lat2)/(lat3-lat2)*(lon3-lon2)
        IF (ilon.lt.clon) includ=.false.

C +---Segment 1-4
C +   -----------

        ilon=lon1+(clat-lat1)/(lat4-lat1)*(lon4-lon1)
        IF (ilon.gt.clon) includ=.false.


C +---Complete list of cells
C +   ======================

        IF (includ) THEN
         nlist=nlist+1
         IF (nlist.gt.MXlist) THEN
          WRITE(6,*) 'The size of the icell and jcell variables has to'
          WRITE(6,*) 'be increased. Please modify MXlist in PRCdes'
          WRITE(6,*) 'subroutine.                             --- STOP'
          WRITE(6,*)
          STOP
         ENDIF
         icell(nlist)=i
         jcell(nlist)=j
        ENDIF


       ENDDO          ! Loop on NST grid-points : END
       ENDDO


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDDO
      ENDDO           ! Loop on LSC grid-points : END

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END

