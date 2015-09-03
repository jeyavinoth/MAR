program indice_Bamber_NESTOR

implicit none

include '/usr/local/include/netcdf.inc'
include 'bamber.inc'


integer :: ncid,NCcode,var_ID,ID__nc,nbr_var
character*90 :: fnamNC,tit_NC

integer  :: k1,k2,l1,l2,x1,y1
real     :: b_lat_min,b_lat_max,dist2
real     :: b_lon_min,b_lon_max


!--------------------------------------------------------------------------------

! NESTOR

call CF_READ2D(NSTFILE,'LON',1,mx,my,1,m_lon)
call CF_READ2D(NSTFILE,'LAT',1,mx,my,1,m_lat)
call CF_READ2D(NSTFILE,'SH',1,mx,my,1,m_sh)
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

! lon/lat in radian

 m_lon = m_lon * pi/180.
 m_lat = m_lat * pi/180.

 b_lon = b_lon * pi/180.
 b_lat = b_lat * pi/180.

 b_lat_min = 9 ; b_lat_max=-9
 b_lon_min = 9 ; b_lon_max=-9

 do i=1,bx ; do j=1,by
  b_lat_min=min(b_lat_min,b_lat(i,j))
  b_lat_max=max(b_lat_max,b_lat(i,j))
  b_lon_min=min(b_lon_min,b_lon(i,j))
  b_lon_max=max(b_lon_max,b_lon(i,j))
 enddo     ; enddo

!--------------------------------------------------------------------------------

! search of the Bamber's indexes in each MAR pixel

 open(unit=210,file='bamber_index_'//XXkm//'.dat',status='replace',form='unformatted')


 B25_i   =9999
 B25_j   =9999
 B25_dist=9999
 B25_cpt =0

 k1=0

 do j=mmy,my-mmy+1 ; do i=mmx,mx-mmx+1


  cpt  = 0. !initialisation du compteur

  if(m_lat(i,j)<b_lat_min.or.m_lat(i,j)>b_lat_max.or.m_lon(i,j)<b_lon_min.or.m_lon(i,j)>b_lon_max) goto 1000

  if(k1==0.or.i==mmx.or.i==1) then   
   k1=1 ; k2=bx ; l1=1 ; l2=by
  endif

  do k=k1,k2 ; do l=l1,l2

   dist = distance(b_lon(k,l),b_lat(k,l),m_lon(i,j),m_lat(i,j))

   if(dist <= sqrt(2.*(range/2.)**2)) then

    dist2=dist
  
    do x=-1,1 ; do y=-1,1
     x1=min(mx,max(1,i+x))
     y1=min(my,max(1,j+y))
     dist2=min(dist2,distance(b_lon(k,l),b_lat(k,l),m_lon(x1,y1),m_lat(x1,y1)))
    enddo       ; enddo

    if(dist <= dist2*1.01) then
                 cpt  = cpt + 1.
     bm_dist(i,j,cpt) = dist+10e-4*cpt/mh
     bm_i   (i,j,cpt) = k
     bm_j   (i,j,cpt) = l
    endif


   endif

   if (cpt>=mh) goto 1000

  enddo ; enddo

  1000 continue

  bm_cpt(i,j) = cpt


! ascending sort following the distance

  do x=1,min(int(range*range*1.1),bm_cpt(i,j))

   dist_min = 500 

   if(x == 1) then
    indice (x) = 0
    tmp_min(x) = 0
   else
    indice (x) = indice (x-1)
    tmp_min(x) = tmp_min(x-1)
   endif

   do h=1,bm_cpt(i,j)
    if((bm_dist(i,j,h) < dist_min    .and. bm_dist(i,j,h) > tmp_min(x)) .or.	   &
       (bm_dist(i,j,h) == tmp_min(x) .and. h > indice(x))) then
     tmp      = h
     dist_min = bm_dist(i,j,h)
    endif
   enddo

   indice (x) = tmp
   tmp_min(x) = bm_dist(i,j,tmp)

   B25_dist(i,j,x) = bm_dist(i,j,tmp)
   B25_i   (i,j,x) = bm_i   (i,j,tmp)
   B25_j   (i,j,x) = bm_j   (i,j,tmp)
   B25_cpt (i,j)   = min(int(range*range*1.1),bm_cpt(i,j))

  enddo

  if(bm_cpt(i,j)==0) then 
   k1=0
  else
   k1=max(1 ,int(-3*range)+min(bm_i(i,j,1),bm_i(i,j,bm_cpt(i,j))))
   k2=min(bx,int( 3*range)+max(bm_i(i,j,1),bm_i(i,j,bm_cpt(i,j))))
   l1=max(1 ,int(-3*range)+min(bm_j(i,j,1),bm_j(i,j,bm_cpt(i,j))))
   l2=min(by,int( 3*range)+max(bm_j(i,j,1),bm_j(i,j,bm_cpt(i,j))))
  endif

  write(210)       i,j                      , &
  	 	   B25_cpt (i,j)            , &
  		  (B25_i   (i,j,x),x=1,mh), &
  		  (B25_j   (i,j,x),x=1,mh), &
  		  (B25_dist(i,j,x),x=1,mh)

  print *,i,j,bm_cpt(i,j)

 enddo ; enddo


 close(210)


!--------------------------------------------------------------------------------

end program


!--------------------------------------------------------------------------------

! distance

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


!--------------------------------------------------------------------------------


