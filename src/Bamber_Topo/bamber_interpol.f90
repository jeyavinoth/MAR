! ATTENTION: a filtering (24/36 in TOPOcor.f) 
! of the ETOPO based topo is needed before!!! 

program interpol_Bamber_NESTOR

implicit none

include '/usr/local/include/netcdf.inc'
include 'bamber.inc'


integer      :: ncid,NCcode,var_ID,ID__nc,nbr_var
character*90 :: fnamNC,tit_NC
real         :: dimval(500)


integer,parameter :: nx    =  mx*2
integer,parameter :: ny    =  my*2

real   d1,d2,area(mx,my),lon1,lon2,lat1,lat2

!--------------------------------------------------------------------------------

! NESTOR

call CF_READ2D(NSTFILE,'LON',1,mx,my,1,m_lon)
call CF_READ2D(NSTFILE,'LAT',1,mx,my,1,m_lat)
call CF_READ2D(NSTFILE,'SH' ,1,mx,my,1,m_sh)
call CF_READ2D(NSTFILE,'SOL',1,mx,my,1,m_sol)

!--------------------------------------------------------------------------------

! Bamber

fnamNC='Greenland_bedrock_topography_and_geometry_JGriggs.nc'

call CF_READ2D(fnamNC,'LON',1,bx,by,1,b_lon)
call CF_READ2D(fnamNC,'LAT',1,bx,by,1,b_lat)

call CF_READ2D(fnamNC,'BED',1,bx,by,1,b_bed)
call CF_READ2D(fnamNC,'SRF',1,bx,by,1,b_srf)

call CF_READ2D(fnamNC,'ICE',1,bx,by,1,b_ice)
call CF_READ2D(fnamNC,'MSK',1,bx,by,1,b_msk)

!--------------------------------------------------------------------------------

! GISM05

!fnamNC='gism05.nc'
!call CF_READ2D(fnamNC,'ICE_GISM',1,mx,my,1,m_gism)

!--------------------------------------------------------------------------------

	!---------------------------------------
	!creation of the binary mask from NESTOR
        !---------------------------------------

 do i=1,mx ; do j=1,my
  if(m_sol(i,j) == 1.) then
   m_msk(i,j) = 0.
  else
   m_msk(i,j) = 1.
  endif
 enddo ; enddo

!--------------------------------------------------------------------------------

 B25_cpt=0

 !WARNING x-1

 open(unit=210,name="bamber_index_"//XXkm//".dat",status="old",form='unformatted')
 do j=mmy,my-mmy+1 ; do i=mmx,mx-mmx+1
   read(210) x,y,B25_cpt(x,y),(B25_i(x,y,h),h=1,mh), &
	    (B25_j(x,y,h),h=1,mh),(B25_dist(x,y,h),h=1,mh)
   
   if(x /= i) then
    print *,"ERROR 1",x,i
    stop
   endif

   if(y /= j) then
    print *,"ERROR 2",y,j
    stop
   endif

!   if(i==20.and.j==69) then
!    print *,m_lon(i,j),m_lat(i,j)
!    print *,"---------------------------------"
!    do h=1,B25_cpt(x,y)
!     print *,B25_i(x,y,h),B25_j(x,y,h),b_lon(B25_i(x,y,h),B25_j(x,y,h)),b_lat(B25_i(x,y,h),B25_j(x,y,h)),B25_dist(x,y,h)
!    enddo
!    stop
!   endif

  enddo ; enddo
 close(210)


!--------------------------------------------------------------------------------
 
! interpolation

 do k=1,bx ; do l=1,by
  if(b_srf(k,l) < 0.) b_srf(k,l) = 0.
 enddo ; enddo

 call interpolation(m_sh,b_srf,m_srf_Bam,P)

 call interpolation(m_sh,b_bed,m_bed_Bam,P)

!--------------------------------------------------------------------------------


 do k=1,bx ; do l=1,by
                                                         b_tmp(k,l) =  0
  if(b_msk(k,l) ==1.or.b_msk(k,l) ==2.or.b_msk(k,l) ==4) b_tmp(k,l) = -1.
 enddo ; enddo

 call interpolation(m_msk,b_tmp,m_msk_Bam,P)

 do k=1,bx ; do l=1,by
                                          b_tmp(k,l) = 0
  if(b_msk(k,l) == 2..or.b_msk(k,l) == 4) b_tmp(k,l) = 1.
 enddo ; enddo

 call interpolation(0*m_msk,b_tmp,m_ice_Bam,P)

!--------------------------------------------------------------------------------

 do i=1,mx ; do j=1,my
  if(m_msk_Bam(i,j)<-0.01) then
   m_ice_Bam(i,j)=min(0.99999,max(0.00001,m_ice_Bam(i,j)))
  endif
 enddo ; enddo

!--------------------------------------------------------------------------------

 do i=1,mx ; do j=1,my

  !---------
  !Spitzberg
  !---------
  if(m_lon(i,j) > 6.) then
   if(m_msk(i,j) == 1) then
    m_ice_Bam(i,j)=0.
   endif
  endif

  !-------

  !Islande
  !-------
  if(m_lon(i,j) > -27. .and. m_lat(i,j) < 67.) then
   if(m_msk(i,j) == 1) then
    m_ice_Bam(i,j)=0.
   endif
  endif

  !-------------
  !Baffin Island
  !-------------
  if(m_lon(i,j) < -58. .and. m_lat(i,j) < 70.) then
   if(m_msk(i,j) == 1) then
    if (m_srf_Bam(i,j) >=1000 ) then
     m_ice_Bam(i,j)=1.
    else
     m_ice_Bam(i,j)=0.
    endif
   endif
  endif

  !----------------
  !Ellesmere Island
  !----------------
  if((m_lon(i,j) < -73  .and. m_lat(i,j) > 74.) .or.	&
     (m_lon(i,j) < -69. .and. m_lat(i,j) > 79.) .or.	&
     (m_lon(i,j) < -62.5.and. m_lat(i,j) > 81.)) then
   if(m_msk(i,j) == 1) then
    if (m_srf_Bam(i,j) >=1000 ) then
     m_ice_Bam(i,j)=1.
    else
     m_ice_Bam(i,j)=0.
    endif
   endif
  endif

  if(m_msk_Bam(i,j)<=-0.50)                                           m_msk_Bam(i,j)=1
  if(m_msk_Bam(i,j)<0.0.and.m_msk_Bam(i,j)>-0.50)                     m_msk_Bam(i,j)=0

 enddo     ; enddo
 
!--------------------------------------------------------------------------------

! ice sheet margin - 100 km
!---------------------------

! m_ice_Bam = 1 if margin > 100km

tmp = max(3,nint(100./range))

 do i=tmp+1,mx-tmp-1 ; do j=tmp+1,my-tmp-1
  if(m_ice_Bam(i,j)>0.99) then
   cpt=0
   do k=-tmp,tmp  ; do l=-tmp,tmp
     if(m_ice_Bam(i+k,j+l)<0.99) cpt=cpt+1
   enddo      ; enddo
   if (cpt==0) m_ice_Bam(i,j)=1
  endif   
 enddo       ; enddo

!--------------------------------------------------------------------------------

! fitering of the sinks
!----------------------

do i=2,mx-1 ; do j=2,my-1
  TMP_sh(i,j)=   m_srf_Bam(i-1,j+1)+        2.*m_srf_Bam(i,j+1)+   m_srf_Bam(i+1,j+1) &
             +2.*m_srf_Bam(i-1,j  )+filtering *m_srf_Bam(i,j  )+2.*m_srf_Bam(i+1,j  ) &
             +   m_srf_Bam(i-1,j-1)+        2.*m_srf_Bam(i,j-1)+   m_srf_Bam(i+1,j-1)
enddo     ; enddo

!do i=2,mx-1 ; do j=2,my-1
! m_srf_Bam (i,j) = TMP_sh(i,j) / (12.+filtering)
!enddo     ; enddo

!--------------------------------------------------------------------------------

! fitering of the sinks
!----------------------

do m=1,1000000

 cpt=0 ; tmp=0
 do i=2,mx-1 ; do j=2,my-1

  if(m_ice_Bam(i,j)<0.90) then

    cpt=0
    if(m_srf_Bam(i+1,j  )<m_srf_Bam(i,j)) cpt=cpt+1
    if(m_srf_Bam(i-1,j  )<m_srf_Bam(i,j)) cpt=cpt+1  
    if(m_srf_Bam(i  ,j+1)<m_srf_Bam(i,j)) cpt=cpt+1
    if(m_srf_Bam(i  ,j-1)<m_srf_Bam(i,j)) cpt=cpt+1  

    tmp1=min(m_srf_Bam(i+1,j),m_srf_Bam(i-1,j),m_srf_Bam(i,j+1),m_srf_Bam(i,j-1))   

    if(cpt==0.and.tmp1>0) then   
      print *,m,i,j,m_srf_Bam(i,j),m_msk_Bam(i,j)
      m_srf_Bam(i,j)=tmp1+2
      tmp=tmp+1
    endif
  endif
 enddo     ; enddo

 if (tmp==0) goto 1000

enddo

1000 continue

!--------------------------------------------------------------------------------

! fitering of the sea pixels
!---------------------------

 call filt_sea_mask(m_msk_Bam,m_msk_Bam)

 do i=1,mx ; do j=1,my
  if(m_msk_Bam(i,j)==0) m_ice_bam(i,j)=0
 enddo     ; enddo

 m_bed_Bam=min(m_bed_Bam,m_srf_Bam)

!--------------------------------------------------------------------------------

d1=0
d2=0
area=0

do i=2,mx-2 ; do j=2,my-2

  d1=0 ; d2=0

  do x=-1,1,2
   lon1=m_lon(i,  j)* pi/180.
   lon2=m_lon(i+x,j)* pi/180
   lat1=m_lat(i,  j)* pi/180
   lat2=m_lat(i+x,j)* pi/180
   d1=d1+distance(lon2,lat2,lon1,lat1)/2.
  enddo

  do x=-1,1,2
   lon1=m_lon(i,j)* pi/180.
   lon2=m_lon(i,j+x)* pi/180
   lat1=m_lat(i,j)* pi/180
   lat2=m_lat(i,j+x)* pi/180
   d2=d2+distance(lon2,lat2,lon1,lat1)/2.
  enddo

  area(i,j)=d1*d2

enddo       ; enddo

!--------------------------------------------------------------------------------
 
fnamNC=filen

call CF_INI_FILE(fnamNC,XXkm)

call CF_CREATE_DIM("time","-",1,1)

DO i=1,mx
 dimval(i)= (i-imez)*range
ENDDO
call CF_CREATE_DIM("x","km",mx,dimval)


DO j=1,my
 dimval(j)= (j-jmez)*range
ENDDO
call CF_CREATE_DIM("y","km",my,dimval)

call CF_CREATE_VAR("LON","Longitude","degrees","-","x","y","-" )
call CF_CREATE_VAR("LAT","Latitude" ,"degrees","-","x","y","-" )
call CF_CREATE_VAR("MSK","Soil Mask"      ,"-","-","x","y","-" )
call CF_CREATE_VAR("ICE","Bamber (2012) Ice Mask"       ,"-","-","x","y","-" )
call CF_CREATE_VAR("SRF" ,"Surface height (without filtering)" ,"m","-","x","y","-" )
call CF_CREATE_VAR("BED","Bedrock height" ,"m","-","x","y","-" )
call CF_CREATE_VAR("AREA","Pixel area" ,"km^2","-","x","y","-" )
call CF_CREATE_VAR("ICE_GISM","GISM05 Ice Mask"       ,"-","-","x","y","-" )

call CF_CREATE_FILE(fnamNC)

call CF_write (fnamNC, 'LON '   , 1  , mx   , my, 1 , m_lon)
call CF_write (fnamNC, 'LAT '   , 1  , mx   , my, 1 , m_lat)
call CF_write (fnamNC, 'MSK'    , 1  , mx   , my, 1 , m_msk_Bam)
call CF_write (fnamNC, 'ICE'    , 1  , mx   , my, 1 , m_ice_Bam)
call CF_write (fnamNC, 'ICE_GISM'    , 1  , mx   , my, 1 , m_gism)
call CF_write (fnamNC, 'BED'    , 1  , mx   , my, 1 , m_bed_Bam)
call CF_write (fnamNC, 'SRF'    , 1  , mx   , my, 1 , m_srf_Bam)
call CF_write (fnamNC, 'AREA'   , 1  , mx   , my, 1 , area)

!--------------------------------------------------------------------------------

end program


!--------------------------------------------------------------------------------

function distance(lon2,lat2,lon1,lat1)	

implicit none

include 'bamber.inc'

 real :: lon1,lat1
 real :: lon2,lat2
 real :: dlat,dlon,a,c

! distance = acos(sin(lat1) * sin(lat2) + &
!	    cos(lat1) * cos(lat2) * cos(lon2 - lon1)) * R

    dlat = (lat2-lat1)
    dlon = (lon2-lon1)
       a =  sin(dLat/2.) * sin(dLat/2.) + cos(lat1) * cos(lat2) * sin(dLon/2) * sin(dLon/2)
       c =  2.d0 * atan2(sqrt(a), sqrt(1-a))
distance =  R * c

end function


!--------------------------------------------------------------------------------------

subroutine interpolation(m_var,b_var,new_var,exp)

implicit none

include 'bamber.inc'


 real    :: m_var   (mx,my)
 real    :: b_var   (bx,by)
 real    :: new_var (mx,my)

 real	 :: numerat (mx,my)			!numerateur de l interpolation
 real    :: denomin (mx,my)			!denominateur de l interpolation
 real	 :: inv_dist(mx,my)

 real    :: exp					!puissance de la distance


 cpt      = 0.
 dist     = 0.
 inv_dist = 0.
 numerat  = 0.
 denomin  = 0.

 do i=1,mx ; do j=1,my

  cpt = B25_cpt(i,j)

  if(cpt == 0) then

   new_var(i,j) = m_var(i,j)

  else

   do h=1,cpt

    x    = B25_i   (i,j,h)
    y    = B25_j   (i,j,h)
    dist = B25_dist(i,j,h)

    if(dist < 1.) dist = 1.

    inv_dist(i,j) = 1./(dist**exp)
    numerat (i,j) = numerat(i,j) + (inv_dist(i,j) * b_var(x,y))
    denomin (i,j) = denomin(i,j) +  inv_dist(i,j)

   enddo

   new_var(i,j) = numerat(i,j)/denomin(i,j)

  !if(new_var(i,j) < 0.) new_var(i,j) = 0.

  endif

  !---------
  !Spitzberg
  !---------
  if(m_lon(i,j) > 6.) then
   new_var(i,j) = m_var(i,j)
  endif

  !-------
  !Islande
  !-------
  if(m_lon(i,j) > -27. .and. m_lat(i,j) < 67.) then
   new_var(i,j) = m_var(i,j)
  endif

  !-------------
  !Baffin Island
  !-------------
  if(m_lon(i,j) < -58. .and. m_lat(i,j) < 70.) then
   new_var(i,j) = m_var(i,j)
  endif

  !----------------
  !Ellesmere Island
  !----------------
  if((m_lon(i,j) < -74. .and. m_lat(i,j) > 74. ) .or.	&
     (m_lon(i,j) < -68.5.and. m_lat(i,j) > 79.5) .or.	&
     (m_lon(i,j) < -66  .and. m_lat(i,j) > 80.8) .or.	&
     (m_lon(i,j) < -64  .and. m_lat(i,j) > 81.2) .or.	&
     (m_lon(i,j) < -61. .and. m_lat(i,j) > 81.8)) then
   new_var(i,j) = m_var(i,j)
  endif

 enddo ; enddo

end subroutine


!--------------------------------------------------------------------------------------

!filtering binary mask after interpolation

subroutine filt_sea_mask(mask,new_mask)

implicit none

include 'bamber.inc'


 real :: mask    (mx,my)
 real :: new_mask(mx,my)


 new_mask = mask

 do i=1,mx,mx-1 ; do j=2,my-1

  cpt = 0

  do l=j-1,j+1
   cpt = cpt + new_mask(i,l)
  enddo

  if(cpt == 2  .and. new_mask(i,j) == 0) new_mask(i,j) = 1

  if(cpt == 1  .and. new_mask(i,j) == 1) new_mask(i,j) = 0

 enddo ; enddo

! do i=2,mx-1 ; do j=2,my-1

!  cpt = 0

!  do k=i-1,i+1 ; do l=j-1,j+1
!   cpt = cpt + new_mask(k,l)
!  enddo ; enddo

!  if(cpt >= 8  .and. new_mask(i,j) == 0) new_mask(i,j) = 1

 !if(cpt <= 1  .and. new_mask(i,j) == 1) new_mask(i,j) = 0

! enddo ; enddo


end subroutine

!--------------------------------------------------------------------------------------



