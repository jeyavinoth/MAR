C   +-------------------------------------------------------------------+
C   |  Subroutine SSTint                             20/06/2004 NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Interpolation of large-scale data to nested grid.                 |
C   |                                                                   |
C   | INPUT  : DATtim: time for which the data is requested             |
C   | ^^^^^^^^ HORint: horizontal interp. type  (1= bilin, 3= bicub)    |
C   |          NST__x: horizontal NST grid (longitude)                  |
C   |          NST__y: horizontal NST grid (latitude )                  |
C   |          NSTsol: soil types                                       |
C   |                                                                   |
C   | INPUT FILE: Global sea surface temperature (Reynolds)             |
C   | ^^^^^^^^^^^                                                       |
C   |                                                                   | 
C   | OUTPUT : NSTsst: Sea surface temperature       ( K )              |
C   | ^^^^^^^^ NST_st: Corrected surface temperature ( K )              |
C   |                                                                   | 
C   +-------------------------------------------------------------------+

      SUBROUTINE SSTint

 
      IMPLICIT NONE


C +---Include files
C +   -------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'      
      INCLUDE 'LOCfil.inc'
      INCLUDE 'NetCDF.inc'


C +---Local variables
C +   ---------------

      INTEGER i,j,k,it1,ierror,nxs,nys,nbchar,Rcode,
     .        nts,SSTcid,lon_ID,lat_ID,tim_ID,sst_ID,msk_ID,it2

      PARAMETER (nxs=  360)
      PARAMETER (nys=  180)
      PARAMETER (nts= 1010)

      INTEGER start1(1),start2(2),start3(3),
     .        count1(1),count2(2),count3(3)

      INTEGER*2 LSCmsk(nxs,nys),LSCtmp(nxs,nys)

      REAL    LSC1Dx(nxs),LSC1Dy(nys),LSCsst(nxs,nys),temp(nxs,nys),
     .        tmp_I2a(nxs,nys),tmp1in(nxs,nys),tmp2in(0:nxs+1,0:nys+1),
     .        SST__1(mx,my),SST__2(mx,my),aux1,aux2

      REAL*8  LSC1Dt(nts),julDAY

      LOGICAL SPHgrd,Vfalse

      CHARACTER*7   namlon,namlat,namSST,namMSK,namTIM

      DATA Vfalse / .false. /

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Open NetCDF file containing LSC data
C +   ====================================

      nbchar = 1

      DO i=1,60
       IF (SST_dir(i:i).ne.' ') nbchar=i
      ENDDO

C +            *****
      SSTcid = NCOPN(SST_dir(1:nbchar) // 
     .              'sst_1981-2001_hanning.nc',NCNOWRIT,Rcode)
C +            *****

      SPHgrd = .true.


C +---Screen message
C +   --------------

      write(6,*) 'Sea surface temperature of Reynolds'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Define variable names in Netcdf file
C +   ====================================

      namlon='lon'
      namlat='lat'
      namTIM='time'
      namMSK='mask'
      namSST='sst'

      lon_ID = NCVID(SSTcid,namlon,Rcode)
      lat_ID = NCVID(SSTcid,namlat,Rcode)
      tim_ID = NCVID(SSTcid,namTIM,Rcode)
      msk_ID = NCVID(SSTcid,namMSK,Rcode)
      sst_ID = NCVID(SSTcid,namSST,Rcode)

      start1(1)=1
      start2(1)=1
      start2(2)=1
      start3(1)=1
      start3(2)=1
      start3(3)=1
      count1(1)=nxs
      count2(1)=nxs
      count2(2)=nys
      count3(1)=nxs
      count3(2)=nys
      count3(3)=nts


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Horizontal coordinates
C +   ======================


      count1(1)=nxs
C +        *****
      CALL NCVGT (SSTcid,lon_ID,start1,count1,LSC1Dx,Rcode)
C +        *****

      count1(1)=nys
C +        *****
      CALL NCVGT (SSTcid,lat_ID,start1,count1,LSC1Dy,Rcode)
C +        *****

      count1(1)=nts
C +        *****
      CALL NCVGT (SSTcid,tim_ID,start1,count1,LSC1Dt,Rcode)
C +        ***** 

C +        ***** 
      CALL NCVGT (SSTcid,msk_ID,start2,count2,LSCmsk,Rcode)
C +        ***** 


C +---Longitudes between -180. and +180.
C +   ----------------------------------

      DO i=1,nxs
       LSC1Dx(i)=LSC1Dx(i)-180.
      ENDDO

      DO j=1,nys
      DO i=1,nxs
       temp(i,j)=LSCmsk(i,j)
      ENDDO
      ENDDO

      DO j=1,nys
      DO i=1,nxs/2
       LSCmsk(      i,j)=temp(nxs/2+i,j)
       LSCmsk(nxs/2+i,j)=temp(      i,j)
      ENDDO
      ENDDO


C +---Correction of the time base of input data
C +   -----------------------------------------

      DO i=1,nts
       LSC1Dt(i) = LSC1Dt(i) + 3.5
      ENDDO


C +---Time for data extraction
C +   ------------------------

C +        ******
      CALL DATcnv (RUNiyr,RUNmma,RUNjda,RUNjhu,DATtim,Vfalse)
C +        ******

      julDAY = REAL(DATtim)/24.0 - 366.0 - 13.0

      it1 = 0
      it2 = 0

      DO i=2,nts
       IF (LSC1Dt(i-1).le.julDAY .and.
     .     LSC1Dt(i  ).gt.julDAY) THEN
        it1 = i-1
        it2 = i
       ENDIF
      ENDDO

      IF (it1.eq.0.or.it2.eq.0) THEN
       WRITE(6,2400) ' No SST available for the date : ',
     .                RUNiyr,RUNmma,RUNjda,RUNjhu
2400   FORMAT (a33,i6,3i4)
      ENDIF

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Horizontal interpolation of sea surface temperature
C +   ===================================================


      IF (it1.ne.0.and.it2.ne.0) THEN

       write(6,*) 'Time steps : ',it1,it2
       write(6,*) 'Data file  : ',SST_dir(1:nbchar) //
     .            'sst_1981-2001_hanning.nc'

C +... First date (for temporal interpolation)

       start3(3) = it1
       count3(3) = 1

C +         *****
       CALL NCVGT (SSTcid,sst_ID,start3,count3,LSCtmp,Rcode)
C +         ***** 
 
       DO j=1,nys
       DO i=1,nxs
        temp(i,j)=LSCtmp(i,j)
       ENDDO
       ENDDO

       DO j=1,nys
       DO i=1,nxs/2
        LSCtmp(      i,j)=temp(nxs/2+i,j)
        LSCtmp(nxs/2+i,j)=temp(      i,j)
       ENDDO
       ENDDO

       DO j=1,nys
       DO i=1,nxs
        LSCsst(i,j)=0.01*REAL(LSCtmp(i,j))+273.15
       ENDDO
       ENDDO

       IF (HORint.EQ.1) THEN      ! Bilinear interpolation

C +           ******
         CALL INTbil (nxs,nys,LSC1Dx,LSC1Dy,LSCsst,SPHgrd,
     .                mx ,my ,NST__x,NST__y,SST__1,tmp2in)
C +           ******

       ELSE IF (HORint.EQ.3) THEN ! Bicubic interpolation

C +           ******
         CALL INTbic (tmp_I2a,tmp1in,
     .                nxs,nys,LSC1Dx,LSC1Dy,LSCsst,
     .                mx ,my ,NST__x,NST__y,SST__1)
C +           ******

       ENDIF


C +... Second date (for temporal interpolation)

       start3(3) = it2
       count3(3) = 1

C +         *****
       CALL NCVGT (SSTcid,sst_ID,start3,count3,LSCtmp,Rcode)
C +         ***** 
 
       DO j=1,nys
       DO i=1,nxs
        temp(i,j)=LSCtmp(i,j)
       ENDDO
       ENDDO

       DO j=1,nys
       DO i=1,nxs/2
        LSCtmp(      i,j)=temp(nxs/2+i,j)
        LSCtmp(nxs/2+i,j)=temp(      i,j)
       ENDDO
       ENDDO

       DO j=1,nys
       DO i=1,nxs
        LSCsst(i,j)=0.01*REAL(LSCtmp(i,j))+273.15
       ENDDO
       ENDDO

       IF (HORint.EQ.1) THEN      ! Bilinear interpolation

C +           ******
         CALL INTbil (nxs,nys,LSC1Dx,LSC1Dy,LSCsst,SPHgrd,
     .                mx ,my ,NST__x,NST__y,SST__2,tmp2in)
C +           ******

       ELSE IF (HORint.EQ.3) THEN ! Bicubic interpolation

C +           ******
         CALL INTbic (tmp_I2a,tmp1in,
     .                nxs,nys,LSC1Dx,LSC1Dy,LSCsst,
     .                mx ,my ,NST__x,NST__y,SST__2)
C +           ******

       ENDIF


C +... Linear interpolation between SST__1 and SST__2

       aux1 = (LSC1Dt(it2)-julDAY)/(LSC1Dt(it2)-LSC1Dt(it1))
       aux2 = 1.0 - aux1

       DO j=1,my
       DO i=1,mx
        NSTsst(i,j) = aux1*SST__1(i,j) + aux2*SST__2(i,j)
       ENDDO
       ENDDO


C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Correction of soil temperature over ocean
C +   =========================================

       DO j = 1,my
       DO i = 1,mx
        IF (NSTsol(i,j).le.2.and.NSTsst(i,j).ge.273.15) THEN
         NST_st(i,j)=NSTsst(i,j)
        ENDIF
       ENDDO
       ENDDO

      ENDIF  ! (it.gt.0.and.it.le.nts)

 
C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---Close the NetCDF file
C +   =====================

C +        ******
      CALL NCCLOS (SSTcid,ierror)
C +        ******

      write(6,*) 

C +   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
