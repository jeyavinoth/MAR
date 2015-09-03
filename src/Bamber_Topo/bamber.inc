!--------------------------------------------------------
!	Configuration of parameters for interpolation    
!	of the Bamber's topography and the Box's mask	
!	onto the NESTOR grid at different resolutions	
!--------------------------------------------------------
!--------------------------------------------------------
!	Parameters have to be activated according	
!	to the spatial resolution of the NESTOR		
!--------------------------------------------------------

!----------- 
! Run Label
!-----------

character*24	  :: NSTFILE = 'NST.1979.01.01.00.GRj.nc'

!---------------------------------------
!Configuration of the general parameters
!---------------------------------------

 real,parameter    :: pi   = 3.1415
 real,parameter    :: R    = 6371.
 real,parameter    :: epsi = 10.E-6

 real,parameter    :: filtering = 36. 

 integer :: i,j,k,l,x,y,h,m
 integer :: cpt,tmp

 real    :: distance
 real	 :: dist
 real    :: dist_min

 real    :: tmp1,tmp2,tmp3

 real	 :: P = 1./4.			!power of the distance

 common/param_gen/	P	


!-----------------------------------------------------------------------
!Configuration for interpolation onto the NESTOR grid at xxkm resolution
!-----------------------------------------------------------------------

  character*4       :: XXkm  = '25km'
  character*100     :: filen = 'GRD-25km-80x135.cdf'
  integer,parameter :: mx    = 80 
  integer,parameter :: my    = 135
  real,parameter    :: range = 25 ! km
  integer,parameter :: imez  = 40
  integer,parameter :: jmez  = 70

!  character*6       :: XXkm  = '12_5km'
!  character*100     :: filen = 'GRD-12_5km-140x250.cdf'
!  integer,parameter :: mx    = 140
!  integer,parameter :: my    = 250
!  real,parameter    :: range = 12.5 ! km
!  integer,parameter :: imez  = 73
!  integer,parameter :: jmez  = 120


!  character*4       :: XXkm  = '15km'
!  character*100     :: filen = 'GRD-15km-115x210.cdf'
!  integer,parameter :: mx    = 115
!  integer,parameter :: my    = 210
!  real,parameter    :: range = 15 ! km
!  integer,parameter :: imez  = 60
!  integer,parameter :: jmez  = 100

!  character*4       :: XXkm  = '20km'
!  character*100     :: filen = 'GRD-20km-95x165.cdf'
!  integer,parameter :: mx    = 95
!  integer,parameter :: my    = 165
!  real,parameter    :: range = 20 ! km
!  integer,parameter :: imez  = 50
!  integer,parameter :: jmez  = 80

!  character*4       :: XXkm  = '25km'
!  character*100     :: filen = 'GRD-25km-80x135.cdf'  
!  integer,parameter :: mx    =  80
!  integer,parameter :: my    = 135
!  real,parameter    :: range = 25 ! km
!  integer,parameter :: imez  = 41
!  integer,parameter :: jmez  = 66

!  character*4       :: XXkm  = '30km'
!  character*100     :: filen = 'GRD-30km.cdf'
!  integer,parameter :: mx    = 80
!  integer,parameter :: my    = 140
!  real,parameter    :: range = 30 ! km
!  integer,parameter :: imez  = 37
!  integer,parameter :: jmez  = 70

!  character*6       :: XXkm  = '32_5km'
!  character*100     :: filen = 'GRD-32_5km-65x110.cdf'
!  integer,parameter :: mx    = 65
!  integer,parameter :: my    = 110
!  real,parameter    :: range = 32.5 ! km
!  integer,parameter :: imez  = 34
!  integer,parameter :: jmez  = 55

 
!  character*4       :: XXkm  = '35km'
!  character*100     :: filen = 'GRD-35km-60x100.cdf'
!  integer,parameter :: mx    = 60
!  integer,parameter :: my    = 100
!  real,parameter    :: range = 35 ! km
!  integer,parameter :: imez  = 32
!  integer,parameter :: jmez  = 50

!  character*4       :: XXkm  = '35km'
!  character*100     :: filen = 'GRD-35km-80x115.cdf'
!  integer,parameter :: mx    = 80
!  integer,parameter :: my    = 115
!  real,parameter    :: range = 35 ! km
!  integer,parameter :: imez  = 40 
!  integer,parameter :: jmez  = 55

! ERIK NOBLE's Big 36km grid
!  character*4       :: XXkm  = '36km'
!  character*100     :: filen = 'GRD-36km-79x112.cdf'
!  integer,parameter :: mx    = 79 
!  integer,parameter :: my    = 112
!  real,parameter    :: range = 36 ! km
!  integer,parameter :: imez  = 40
!  integer,parameter :: jmez  = 55

!  character*6       :: XXkm  = '37_5km'
!  character*100     :: filen = 'GRD-37_5km-55x95.cdf'
!  integer,parameter :: mx    = 55
!  integer,parameter :: my    = 95
!  real,parameter    :: range = 37.5 ! km
!  integer,parameter :: imez  = 29
!  integer,parameter :: jmez  = 46

!  Erik, Test 50 km domain
!  character*4       :: XXkm  = '50km'
!  character*100     :: filen = 'GRD-50km.cdf'
!  integer,parameter :: mx    = 60
!  integer,parameter :: my    = 80
!  real,parameter    :: range = 50 ! km
!  integer,parameter :: imez  = 35
!  integer,parameter :: jmez  = 40

  integer,parameter :: mmx   = 5 ! [mmx,mx-mmx+1]
  integer,parameter :: mmy   = 5 ! [mmy,my-mmy+1]
  integer,parameter :: mh    = 1.1 * range * range

!---------------------------------------------------
!Configuration of the general parameters from NESTOR
!---------------------------------------------------

 real :: m_lon (mx,my)
 real :: m_lat (mx,my)

 real :: m_dist(mx,my)			!distance

 real :: m_sh  (mx,my)			!topography from NESTOR
 real :: m_sol (mx,my)			!mask from NESTOR
 real :: m_gism (mx,my)			!mask from GISM05

 real :: m_msk (mx,my)			!binary mask from NESTOR

	common/param_NESTOR/	m_lon,m_lat,m_dist		, &
				m_sh,m_sol,m_msk


!---------------------------------------------------------------
!Configuration of the parameters for the interpolated topography
!---------------------------------------------------------------

 integer :: bm_cpt   (mx,my)		!number of points to be interpolated
 integer :: bm_i     (mx,my,mh)
 integer :: bm_j     (mx,my,mh)		!indices i and j
 real    :: bm_dist  (mx,my,mh)		!distance

 integer :: B25_cpt  (mx,my)
 integer :: B25_i    (mx,my,mh)
 integer :: B25_j    (mx,my,mh)
 real    :: B25_dist (mx,my,mh)

 integer :: indice   (mh)			
 real    :: tmp_min  (mh)	

 real    :: m_srf_Bam (mx,my)		!Bamber topography interpolated onto the NESTOR grid
 real	 :: m_bed_Bam(mx,my)		!interpolated bed

 real	 :: m_msk_Bam(mx,my)		!interpolated sea-land mask
 real	 :: m_ice_Bam(mx,my)		!interpolated ice mask

 real    :: m_sol_tmp(mx,my),tmp_sh(mx,my)
 real    :: m_sol_Bam(mx,my)		!interpolated and corrected sea-land mask

	common/param_Bam_NEST/	bm_cpt,bm_i,bm_j,bm_dist	, &
				B25_cpt,B25_i,B25_j,B25_dist	, &
				indice,tmp_min,m_srf_Bam        , &
				m_msk_Bam,m_sol_tmp,m_sol_Bam



!---------------------------------------------------
!Configuration of the Bamber's topography parameters
!---------------------------------------------------

 integer,parameter :: bx = 1501
 integer,parameter :: by = 2801

 real :: b_lon    (bx,by)
 real :: b_lat    (bx,by)

 real :: b_tmp    (bx,by)

 real :: b_dist   (bx,by)		!distance 1km Bamber to NESTOR

 real :: b_srf    (bx,by)		!surface height
 real :: b_bed    (bx,by)		!bedrock height
 real :: b_thick  (bx,by)		!ice sheet thickness

 real :: b_msk    (bx,by)		!sea-land mask
 real :: b_ice    (bx,by)		!ice mask

	common/param_Bamber/	b_lon,b_lat	, &
				b_dist,b_srf,b_bed,b_thick	, &
				b_msk,b_ice






















