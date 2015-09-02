      Subroutine SBCnew      
C +------------------------------------------------------------------------+
C | MAR SURFACE  XF                                        27-05-2014  MAR |
C |   SubRoutine SBCnew for Greenland/Svalbard 3D simulation               |
C |                                                                        |
C |                         Simulation GRD                                 |
C +------------------------------------------------------------------------+      
C +
      Implicit none

C     1) General Variables
C     ====================

      include 'MARCTR.inc'
      include 'MARphy.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'

      character*10 TypeGL

      CHARACTER*5     :: RCP_CMIP5
      common/c_RCP_CMIP5/RCP_CMIP5
    
      ! Name of simulation

                              TypeGL ='ERR'

      if(mx==175.and.my==300) TypeGL ='GR10km'
      if(mx==306.and.my==312) TypeGL ='GR11km'
      if(mx==140.and.my==250) TypeGL ='GR12_5km' ! Greenland
      if(mx==115.and.my==210) TypeGL ='GR15km' 
      if(mx==95 .and.my==165) TypeGL ='GR20km'  
      if(mx==80 .and.my==135) TypeGL ='GR25km'  
      if(mx==68 .and.my==118) TypeGL ='GR30km' 
      if(mx==65 .and.my==110) TypeGL ='GR32_5km' 
      if(mx==60 .and.my==100) TypeGL ='GR35km' 
      if(mx==55 .and.my== 95) TypeGL ='GR37_5km'


      if(mx==141.and.my==121) TypeGL ='AN50km'   ! Antarctica
      if(mx==175.and.my==151) TypeGL ='AN40km'

      if(mx==85 .and.my== 95) TypeGL ='SVa'      ! Svalbard
      if(mx==130.and.my==165) TypeGL ='AC50km'   ! CORDEX 50km

      if(TypeGL=="ERR".and.iterun<=1) then
       print*,"ERROR : you must define a new domain in usr/SBCnew.f"
       print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      endif

C     2) INIsnow 
C     ==========
   
      IF(itexpe ==  0) call INIsnow(TypeGL)  

C     3) UPDsnow 
C     ==========

      IF ( iterun >  1 .and. mod(iterun,12*3600/int(dt))==0) 
     . call UPDsnow

C     4) OUTone 
C     =========

      if(njmbGE(2)==1.and.
     .  (trim(TypeGL)=="SVa".or.trim(TypeGL)=="GR25km")) then

       IF ( (minuGE >  0 .and. minuGE <=  0 + dt/60.) .or.
     .      (minuGE > 30 .and. minuGE <= 30 + dt/60.) .or.
     .       iterun==nterun)           
     .   call OUTone(TypeGL)      

       IF (iterun==1) call OUTone(TypeGL) 

      endif
C     5) OUTsta 
C     =========

c      if(iterun>=1)then
c       call OUTsta('066','058', 1)
c       call OUTsta('063','053', 2)
c      endif

C     6) FILsnow 
C     ==========

      IF ( iterun >  1 .and. mod(iterun,3600/int(dt))==0) 
     . call FILsnow

C     7) ASSsnow 
C     ==========

!     IF (mmarGE>4.and.mmarGE<10) call ASSsnow
      
      return
      end

C +------------------------------------------------------------------------+
      subroutine ice_sheet_model_coupling
C +------------------------------------------------------------------------+
C +
      implicit none
C +
C +--General Variables
C +  =================

      include 'MARphy.inc'
      include 'MARCTR.inc'
      include 'MAR_SV.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
      include 'MAR_DY.inc'
      include 'MAR_LB.inc'
      include 'MAR_SL.inc'
      include 'MAR_SN.inc'
      include 'MAR_BS.inc'
      include 'MAR_IO.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MAR_IB.inc'
      include 'MARdSV.inc'

      integer    ,parameter ::  mxx1 = 60   ! 25km
      integer    ,parameter ::  myy1 = 110  ! 25km
      integer    ,parameter ::  xx1  = 9    ! 25km
      integer    ,parameter ::  yy1  = 17   ! 25km

c     integer    ,parameter ::  mxx1 = 51   ! 30km
c     integer    ,parameter ::  myy1 = 93   ! 30km
c     integer    ,parameter ::  xx1  = 9    ! 30km
c     integer    ,parameter ::  yy1  = 16   ! 30km

c     integer    ,parameter ::  mxx1 = 42   ! 35km
c     integer    ,parameter ::  myy1 = 78   ! 35km
c     integer    ,parameter ::  xx1  = 9    ! 35km
c     integer    ,parameter ::  yy1  = 15   ! 35km

      real       ,parameter :: convrd = 180.0/3.141592   ! rad  => deg
      real       ,parameter :: convhd = 15.0             ! hour => deg

      real tmp1(mxx1,myy1),file_OK
      real lat1(mx,my),lon1(mx,my),bed(mx,my)
      real newmsk(mx,my),newsh(mx,my),oldsh(mx,my),oldmsk(mx,my) 

      character*100 file_name
    
      logical       file_exists

      write(6,12) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE
   12 format('XF WARNING: call of ice_sheet_model_coupling: ',i5,4i3)

      file_name="GRISLI4MAR.cdf"

      open(unit=10,file="MARscenario.ctr",status="old",err=1000)
      rewind     10
      read(10,*,end=1001) 
      read(10,*,end=1001)
      read(10,'(a100)',end=1001) file_name
 1000 continue
 1001 continue
      close(10)

      INQUIRE(FILE=trim(file_name),EXIST=file_exists)

      if(file_exists) then
       
       lon1=0 ; lat1=0 ; bed=0 ; newsh=0 ; newmsk=0

       call CF_READ2D(trim(file_name),'LON' ,1,mxx1,myy1,1,tmp1)
       do i=1,mxx1 ; do j=1,myy1 
        lon1(i+xx1,j+yy1)=tmp1(i,j)
       enddo ; enddo 

       call CF_READ2D(trim(file_name),'LAT' ,1,mxx1,myy1,1,tmp1)
       do i=1,mxx1 ; do j=1,myy1 
        lat1(i+xx1,j+yy1)=tmp1(i,j)
       enddo ; enddo 
  
       call CF_READ2D(trim(file_name),'BED' ,1,mxx1,myy1,1, tmp1)
       do i=1,mxx1 ; do j=1,myy1 
        bed(i+xx1,j+yy1)=tmp1(i,j)
       enddo ; enddo 

       call CF_READ2D(trim(file_name),'newSH' ,1,mxx1,myy1,1,tmp1)
       do i=1,mxx1 ; do j=1,myy1 
        newsh(i+xx1,j+yy1)=tmp1(i,j)
       enddo ; enddo 

       call CF_READ2D(trim(file_name),'newMSK',1,mxx1,myy1,1,tmp1)
       do i=1,mxx1 ; do j=1,myy1 
        newmsk(i+xx1,j+yy1)=tmp1(i,j)
       enddo ; enddo 

       oldmsk = mskSNo(:,:,1)
       oldsh  = sh 

       do i=1,mx ; do j=1,my

         if(newmsk(i,j)>0)then

         !SH      
         if(abs(newsh(i,j)-sh(i,j))>0.1) then

          sh(i,j)    = newsh(i,j)
          sh(i,j)    = max(sh(i,j),  bed(i,j))

          write(6,13) i,j,oldsh(i,j), sh(i,j)
   13     format(2i4," ice sheet SH:",f7.1,'=>',f7.1) 
         endif

         ! MSK
         if(mskSNo(i,j,1)>0.and.mskSNo(i,j,1)<100.and.
     .    abs(lon1(i,j)-GElonh(i,j)*convhd)<0.01 .and.
     .    abs(lat1(i,j)-GElatr(i,j)*convrd)<0.01 .and.
     .    abs(newmsk(i,j)-mskSNo(i,j,1))>0.01)   then
 
          mskSNo(i,j,1) = newmsk(i,j)
          mskSNo(i,j,1) = max( 0.001,mskSNo(i,j,1))
          mskSNo(i,j,1) = min(99.999,mskSNo(i,j,1))
          ifraTV(i,j,1) =            mskSNo(i,j,1)
          ifraTV(i,j,2) = 100. -     mskSNo(i,j,1)

          write(6,14) i,j,oldmsk(i,j), mskSNo(i,j,1)
   14     format(2i4," ice sheet MSK:",f7.2,'=>',f7.2) 
         endif

        endif

       enddo      ; enddo 

      else

       print *,trim(file_name)//" not found"

      endif

      end subroutine ice_sheet_model_coupling


C +------------------------------------------------------------------------+
      subroutine ASSsnow
C +------------------------------------------------------------------------+
C | MAR SURFACE  XF                                                    MAR |
C +------------------------------------------------------------------------+
C +
      implicit none
C +
C +--General Variables
C +  =================

      include 'MARphy.inc'
      include 'MARCTR.inc'
      include 'MAR_SV.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc' 
      include 'MAR_GE.inc'
      include 'MAR_DY.inc'
      include 'MAR_LB.inc'
      include 'MAR_SL.inc'
      include 'MAR_SN.inc'
      include 'MAR_BS.inc'
      include 'MAR_IO.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MAR_IB.inc'

C +--Local   Variables
C +  =================

      real,parameter :: melt_thrsd = 8.5  ! mmWE/day
      real,parameter :: dzsn_thrsd = 0.05 ! m 

      integer        :: day,day_1st_may,day_current,kksn,n
      integer        :: ASS_up,ASS_do

      character*20   :: filename
      character*4    :: YYYYc

      real           :: melt_current,thrsd,dz,dzsn
    
      real           :: tmp1(60,112),tmp2(60,112)
      real           :: msk_sat(mx,my),melt_sat(mx,my)

C +   SMMR/SSMI data set reading

       
      day_1st_may=122-min(1,mod(iyrrGE,4))

      day_current=njyrGE(mmarGE)+
     .            njybGE(mmarGE)*max(0,1-mod(iyrrGE,4))+jdarGE

      day        =day_current-day_1st_may+1  

      write(YYYYc,'(i4)') iyrrGE

      msk_sat=0 ; melt_sat=1

      filename='MELT_'//YYYYc//'.nc'

      call CF_READ2D(filename,'MSK_SAT',1  ,60,112,1,tmp1)
      do i=1,60 ; do j=1,112 
       msk_sat(i+9,j+20)=tmp1(i,j)
      enddo     ; enddo 

      call CF_READ2D(filename,'MELT02' ,        day   ,60,112,1,tmp1)
      call CF_READ2D(filename,'MELT02' ,min(153,day+1),60,112,1,tmp2)
      do i=1,60 ; do j=1,112
       if(tmp1(i,j)==0)                          melt_sat(i+9,j+20)=1
       if(tmp1(i,j)==1)                          melt_sat(i+9,j+20)=2
       if(tmp1(i,j)==tmp2(i,j).and.tmp1(i,j)==1) melt_sat(i+9,j+20)=3
       if(tmp1(i,j)==tmp2(i,j).and.tmp1(i,j)==0) melt_sat(i+9,j+20)=0
      enddo     ; enddo

      do i=1,mx ; do j=1,my ; n=1
       if (mskSNo(i,j,1)>=90.and.msk_sat(i,j) >= 3) then

        melt_current=-1*(wem_IB(i,j,n) - wem0IB(i,j,n))
 
        dzsn=0 ; k = nssSNo(i,j,n)+1
        do while(dzsn<=dzsn_thrsd.or.k>nssSNo(i,j,n)-2)
         k    = k -1  
         dzsn = dzsn+dzsSNo(i,j,n,k)
         kksn = k 
        enddo

        ! 15hTU = midday

        ASS_up=0
        ASS_do=0

        do k=1,3
         if(melt_sat(i,j) <= 1                              .and.
     .      melt_current  >= 1.1*melt_thrsd*(3.+2.*real(k))/10. .and.
     .      jhurGE        <= 15+(k-1)*5  ) then 
                    ASS_up = k
                    thrsd  = melt_thrsd*(3.+2.*real(k))/10.
         endif
        enddo   ! 15h => 5/10 ; 20h => 7/10 ; 25h => 9/10 

        if(melt_sat(i,j)<=1.and.melt_current> melt_thrsd) ASS_up=0

        do k=1,4
         if(melt_sat(i,j) >=2                     .and.
     .      melt_current  <= 0.9*melt_thrsd*real(k)/4..and.
     .      jhurGE        >= 15+k*2) then 
          ASS_do           = k
          thrsd            = melt_thrsd*real(k)/4.
         endif
        enddo ! 17h => 1/4 ; 19h => 2/4 ; 21h => 3/4 ; 23h => 4/4 

        if(ASS_up>=1) then 

         dz=0
         do k=nssSNo(i,j,n),kksn,-1

          tisSNo(i,j,n,k) = min(tisSNo(i,j,n,k),
     .                      273.15-(dzsn-dz)/dzsn)
          dz              = dzsSNo(i,j,n,k) + dz

         enddo

         write(*,*)  ' ' 
         write(*,11) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE,i,j,
     .               melt_current,melt_thrsd,ASS_up
   11    format(' ASSsnow (up) at',i5,4i3,
     .          ' for (',i3,','i3,') : ',f5.2,'>',f6.2,i2)
         write(*,*)  ' '

        endif

        if(ASS_do>=1)then
 
         dz=0
         do k=nssSNo(i,j,n),kksn,-1

          tisSNo(i,j,n,k) = 273.15 + (dzsn-dz)/dzsn
          dz              = dzsSNo(i,j,n,k) + dz

         enddo

         write(*,*)  ' ' 
         write(*,12) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE,i,j,
     .               melt_current,melt_thrsd,ASS_do
   12    format(' ASSsnow (down) at',i5,4i3,
     .          ' for (',i3,','i3,') : ',f5.2,'<',f6.2,i2)
         write(*,*)  ' '

        endif

       endif
      enddo ; enddo 

      end subroutine ASSsnow

C +------------------------------------------------------------------------+
C +------------------------------------------------------------------------+
      subroutine INIsnow(TypeGL)
C +------------------------------------------------------------------------+
C | MAR SURFACE  XF                                                        |
C |   SubRoutine INIsnow initialises the SNOW MODEL                        |
C | Modified to run with mw = 5                Charlotte Lang   13/03/2015 |
C +------------------------------------------------------------------------+
C +
      implicit none
C +
C +--General Variables
C +  =================

      include 'MARphy.inc'
      include 'MARCTR.inc'
      include 'MAR_SV.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_GE.inc'
      include 'MAR_DY.inc'
      include 'MAR_LB.inc'
      include 'MAR_SL.inc'
      include 'MAR_SN.inc'
      include 'MAR_BS.inc'
      include 'MAR_IO.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MAR_IB.inc'
      include 'MARdSV.inc'

C +--Local   Variables
C +  =================

      real       Profil_10(10),depth,ela,Profil_15(15)
      real       snwae(mx,my,nsno),znsn(mx,my,nsno)
      real       ini_snow(mx,my),ann_temp(mx,my)
      real       g2s,denss
      real       ice_depth,distup
      integer    nbr_layer

      integer    ni,nj,nk,n,isn
      integer    i_sea, j_sea,i_tundra, j_tundra
      integer    i_dry, j_dry,i_abla, j_abla
      integer    i_perco, j_perco, INI     

      real      ,parameter :: convhd = 15.0           ! hour => deg
      real      ,parameter :: convrd = 180.0/3.141592 ! rad  => deg
 
      character*10 TypeGL
      character*99 Filename
      character*3  mxc,myc

      real       ro_ini(mx,my,10),ti_ini(mx,my,10)
      real       g1_ini(mx,my,10),g2_ini(mx,my,10),zn_ini(mx,my)
      real       tmp1(mx,my)

      real       ns1(mx,my)     
      real       ni1(mx,my)     
      real       dz1(mx,my,nsno)
      real       nh1(mx,my,nsno)
      real       g11(mx,my,nsno)
      real       g21(mx,my,nsno)
      real       ti1(mx,my,nsno)
      real       ro1(mx,my,nsno)
      real       wa1(mx,my,nsno)
      
      character*4 iyrrGEc
      
      ! Antarctica
      ! ----------
      ! Snow surface density (~first meter of snow) computed
      ! from RACMO2-27km outputs following Ligtenberg 2011 [CA]
      integer mii,mjj
      parameter (mii=262,mjj=240)
      real lon_ini(mii,mjj),lat_ini(mii,mjj)
      real ros_ini(mii,mjj),tas_ini(mii,mjj),smb_ini(mii,mjj)
      real ros_min_ant,ros_max_ant
      real tas_min_ant,tas_max_ant
      real smb_min_ant,smb_max_ant
      real x_ini,y_ini
      integer i_ini,j_ini
      !+  Interpolated surface snow density
      real nb_pts(mx,my)
      real rosSNo_ini(mx,my),tasSNo_ini(mx,my),smbSNo_ini(mx,my)
      !+  local variables
      real x0,y0,dh
      integer ii,jj
      real ro_i,ro_w,R_gc,E_1,E_2,Ao
      real dz,ro_up,tas,fe_ro,deltaPa,drhodz
      logical    verboseAntarctica
      parameter (verboseAntarctica=.false.)
      integer firstPoint

C +--DATA
C +  ====

      data    Profil_10 
     .        /0.005,0.01,0.02,0.03,0.04,0.05,0.075,0.10,0.25,0.42/

      data    Profil_15
     .        /0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,
     .         0.08,0.09,0.1,0.12,0.14,0.17/

C ............x m in 10 layers


      ann_temp = tfsnow + 48.38 - (0.007924 * sh)
     .         - (0.7512 * (GElatr /  degrad))

      ! Mean Climatological Ann. Temperature

C +--TT srf 
C +  ======

      DO j=1,my ; DO i=1,mx       
       tairSL(i,j) = pktaDY(i,j,mz)
     .             *(pstDYn(i,j)+ptopDY)**cap            
       tairSL(i,j) = min(tairDY(i,j,mz),tairSL(i,j))
       tairSL(i,j) = tairDY(i,j,mz) !+CA+!
      END DO ; END DO 
                
C +--Mask initialisation
C +  ===================

      do i=1,mx ; do j=1,my
       do k=1,mw
       if(isolSL(i,j)<=2)              mskSNo(i,j,k) =0     
       if(isolSL(i,j)>=3 .and. k/=nsx) mskSNo(i,j,k) =ifraTV(i,j,k) !ifraTV = SFR = NSTsfr
       if(isolSL(i,j)>=3 .and. k==nsx) mskSNo(i,j,k) =0.            !ifraTV = SFR = NSTsfr
       enddo
      enddo         ; enddo

c     endif

C +--Surface initialisation
C +  ======================

C     1) Ocean and Sea Ice Points
C     ---------------------------

      dzsSNo=0. ; rosSNo=0. ; g1sSNo=0. ; g2sSNo=0.
      nhsSNo=0. ; tisSNo=0. ; wasSNo=0. ; agsSNo=0.
      nssSNo=0  ; nisSNo=0  ; issSNo=0  ; snohSN=0
      SWaSNo=0. 
                 
      do j=1,my ; do i=1,mx
       
       DO n=1,nsx
      
        TsrfSL(i,j,n)     =  tairSL(i,j) 
        TvegTV(i,j,n)     =  tairSL(i,j) 
  
        if((isolSL(i,j).lt.3)) then

         i_sea            =  i  
         j_sea            =  j
         isolTV(i,j)      =  0
         SLsrfl(i,j,1)    =  1.  
         ifraTV(i,j,1)    =  100  
         AlbSTV(i,j)      =  0.15 ! / 2 in SISVAT 
                                  ! to not have problem when the ice sea melt
         do isn = 1,llx
          TsolTV(i,j,n,isn)=  SST_LB(i,j) 
          eta_TV(i,j,n,isn)=  1.
         end do          
       
        end if       
       
       END DO
          
      end do ; end do

C     2) Snow zone
C     ------------

      INI = 1 

                  write(mxc,'(i3)') mx
       if(mx<100) write(mxc,'(i2)') mx
                  write(myc,'(i3)') my
       if(my<100) write(myc,'(i2)') my
    
      if (TypeGL(1:2) == "GR") then
       INI = 2

       Filename='MARini-'//trim(TypeGL)//'-'//
     .           trim(mxc)//'x'//trim(myc)//'.cdf'

       nbr_layer=10
       call CF_READ3D(TRIM(Filename),'RO',1,mx,my,nbr_layer,ro_ini)
       call CF_READ3D(TRIM(Filename),'TI',1,mx,my,nbr_layer,ti_ini)
       call CF_READ3D(TRIM(Filename),'G1',1,mx,my,nbr_layer,g1_ini)
       call CF_READ3D(TRIM(Filename),'G2',1,mx,my,nbr_layer,g2_ini)

       do i=1,mx ; do j=1,my
        if(mskSNo(i,j,1)>0.and.ro_ini(i,j,1)==0) then
         print *,"ERROR INIsnow",i,j
        endif      
       enddo     ; enddo

      endif

      if(TypeGL(1:2) == "SV") INI = 4

      if (TypeGL == "GR15km" .or. TypeGL == "GR12_5km" .or.
     .    TypeGL == "GR20km" .or. TypeGL == "GR11km"   .or.
     .    TypeGL == "AC50km" .or. TypeGL == "GR10km"   .or.
     .    TypeGL == "SV7.5km".or. TypeGL == "SV5km") then

       INI = 3

       write(iyrrGEc,'(i4)') iyrrGE

       Filename='MARini-'//trim(TypeGL)//'-'//
     .           trim(mxc)//'x'//trim(myc)//'-'//
     .           iyrrGEc//'.cdf'

       if (TypeGL == "AC50km")
     . Filename='MARini-AC50km-130x165-all.cdf'

       print *,"Read of "//trim(Filename)

       call CF_READ3D(TRIM(Filename),'DZ1',1,mx,my,nsno,dz1)
       call CF_READ3D(TRIM(Filename),'NH1',1,mx,my,nsno,nh1)
       call CF_READ3D(TRIM(Filename),'G11',1,mx,my,nsno,g11)
       call CF_READ3D(TRIM(Filename),'G21',1,mx,my,nsno,g21)
       call CF_READ3D(TRIM(Filename),'TI1',1,mx,my,nsno,ti1)
       call CF_READ3D(TRIM(Filename),'RO1',1,mx,my,nsno,ro1)
       call CF_READ3D(TRIM(Filename),'WA1',1,mx,my,nsno,wa1)
       call CF_READ2D(TRIM(Filename),'NS1',1,mx,my,1,ns1)
       call CF_READ2D(TRIM(Filename),'NI1',1,mx,my,1,ni1)

      endif

      if(TypeGL(1:2) == "AN") then
        INI = 5
        
        firstPoint = 0.
        
        Filename='ANT_ini_racmo2-27km.cdf'
        call CF_READ3D(TRIM(Filename), 'lon',1,mii,mjj,1,lon_ini)
        call CF_READ3D(TRIM(Filename), 'lat',1,mii,mjj,1,lat_ini)
        call CF_READ3D(TRIM(Filename),'rhos',1,mii,mjj,1,ros_ini)
        call CF_READ3D(TRIM(Filename),  'st',1,mii,mjj,1,tas_ini)
        call CF_READ3D(TRIM(Filename),'accu',1,mii,mjj,1,smb_ini)
        ros_min_ant = 300.
        ros_max_ant = 470.
        tas_min_ant = -70. + 273.15
        tas_max_ant = 0. + 273.15
        smb_min_ant = -50.
        smb_max_ant = 5000.
        
        dh=dx/1000.
        x0=(1-imez)*dh
        y0=(1-jmez)*dh
        
        nb_pts = 0
        rosSNo_ini = 0.
        tasSNo_ini = 0.
        smbSNo_ini = 0.
        do jj=1,mjj
        do ii=1,mii
c         call StereoSouth_inverse(lon_ini(ii,jj),lat_ini(ii,jj),
c    .                             GEddxx,x_ini,y_ini)
          i_ini = nint((x_ini - x0) / dh)
          j_ini = nint((y_ini - y0) / dh)
          nb_pts(i_ini,j_ini) = nb_pts(i_ini,j_ini) + 1
          !+ Interpolation
          rosSNo_ini(i_ini,j_ini) = rosSNo_ini(i_ini,j_ini) +
     .                              ros_ini(ii,jj)     
          tasSNo_ini(i_ini,j_ini) = tasSNo_ini(i_ini,j_ini) +
     .                              tas_ini(ii,jj)
          smbSNo_ini(i_ini,j_ini) = smbSNo_ini(i_ini,j_ini) +
     .                              smb_ini(ii,jj)
        enddo
        enddo
        do j=1,my
        do i=1,mx
          if (isolSL(i,j)>=3) then
            if (nb_pts(i,j).ge.1) then
              rosSNo_ini(i,j) = rosSNo_ini(i,j) / nb_pts(i,j)
              tasSNo_ini(i,j) = tasSNo_ini(i,j) / nb_pts(i,j)
              smbSNo_ini(i,j) = smbSNo_ini(i,j) / nb_pts(i,j)
              !+ Min-Max values
              rosSNo_ini(i,j) = max(min(rosSNo_ini(i,j),ros_max_ant),
     .                                                  ros_min_ant)
              tasSNo_ini(i,j) = max(min(tasSNo_ini(i,j),tas_max_ant),
     .                                                  tas_min_ant)
              smbSNo_ini(i,j) = max(min(smbSNo_ini(i,j),smb_max_ant),
     .                                                  smb_min_ant)
            endif
          endif
        enddo
        enddo
      endif ! (TypeGL(1:2) == "AN")

      write (6,*)
      write (6,150) jdarGE,mmarGE,iyrrGE,trim(TypeGL(1:7)),
     .              mx,my,mz,mw,nsno,INI
 150  format('INIsnow at ',i2,'-',i2,'-',i4,':',
     .       ' type=',a7,' mx=',i4,' my='  ,i4,
     .       ' mz='  ,i3,' mw=',i3,' nsno=',i3,' INI=',i2)
      write (6,*)


      do j=1,my ; do i=1,mx ; if (isolSL(i,j) >= 3) then

       IF(GElonh(i,j)*convhd>-43) then ! Equilibrium line (m)
        ELA                 = -32759.680   + 1001.782  
     .                      *  GElatr(i,j) * convrd
     .                      -  7.331       * GElatr(i,j) * GElatr(i,j)
     .                      *  convrd      * convrd      
       ELSE
        ELA                 = -23201.445   + 746.249   
     .                      *  GElatr(i,j) * convrd
     .                      -  5.640       * GElatr(i,j) * GElatr(i,j)
     .                      *  convrd      * convrd      
       END IF

        ELA                 = max(650.,ELA)


C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       denss      = max(450.,1700.*(1.-(sh(i,j)/(1.5*ELA))**2))
       g2s        = max(  5.,  20.*(1.-(sh(i,j)/(1.5*ELA))**2))
 
       ! Svalbard
       if(INI==4) then 
         g2s           = 3
         ann_temp(i,j) = 273.15 - 10
       endif
       
       ! Antarctica
       if(INI==5) then
         g2s   = 3
         ! mean annual temperature
         ann_temp(i,j) = tasSNo_ini(i,j)
       endif

      ! firn density: exp      interpollated between 920 and denss
      ! firn temp   : linearly interpollated between TairSL and T annual
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
       DO n=1,nsx 

        if(mskSNo(i,j,n)>0.and.n/=nsx) then
         nbr_layer=10
        else
         nbr_layer=0
        endif

        nisSNo(i,j,n)         =  0.
                    ice_depth = 15.
        depth                 =  ice_depth


        if (nbr_layer> 0 ) then
        do k = 1, nbr_layer

         dzsSNo(i,j,n,k)    =  max(0.01,
     .                         ice_depth*Profil_10(nbr_layer-k+1))

         depth              =  depth - dzsSNo(i,j,n,k)/2. !+CA+! 

         nhsSNo(i,j,n,k)    =   0.   
         g1sSNo(i,j,n,k)    =  99.
         g2sSNo(i,j,n,k)    =  99.
         wasSNo(i,j,n,k)    =   0.
         agsSNo(i,j,n,k)    =   0. 

         tairSL(i,j)        =  min(273.,tairSL(i,j))
         
         distup             =  min(1.,max(0.,1.-depth/ice_depth))

         tisSNo(i,j,n,k)    =  ann_temp(i,j) *(1.-distup**2)
     .                      +  tairSL(i,j)   *    distup**2

         rosSNo(i,j,n,k)    =  denss         *(1.-distup**2)
     .                      +  300.          *    distup**2

         ! Svalbard
         if (INI==4)
     .   rosSNo(i,j,n,k)    = min(920.,max(500.,
     .                        500. + (400.-sh(i,j)))) 

         ! All
         if (k>=nbr_layer-5)
     .   rosSNo(i,j,n,k)    =  min(250.+ 50.*(nbr_layer-k),
     .                         rosSNo(i,j,n,k))

         if(rosSNo(i,j,n,k)>=roCdSV+70) rosSNo(i,j,n,k)=920
         if(rosSNo(i,j,n,k) <roCdSV+70) g2sSNo(i,j,n,k)=g2s

         ! Greenland
         if(INI==2) then 
          rosSNo(i,j,n,k)   =  ro_ini(i,j,nbr_layer-k+1)
          tisSNo(i,j,n,k)   =  ti_ini(i,j,nbr_layer-k+1)*(1.-distup**2)
     .                      +  tairSL(i,j)  *                distup**2
          g1sSNo(i,j,n,k)   =  g1_ini(i,j,nbr_layer-k+1)
          g2sSNo(i,j,n,k)   =  g2_ini(i,j,nbr_layer-k+1)
         endif

         if (INI==3)
     .   tisSNo(i,j,n,k)    = ti1(i,j,k)

         ! All
         rosSNo(i,j,n,k)    =  max(200.,rosSNo(i,j,n,k))
         rosSNo(i,j,n,k)    =  min(920.,rosSNo(i,j,n,k))
         tisSNo(i,j,n,k)    =  min(273.,tisSNo(i,j,n,k))

         depth              =  depth - dzsSNo(i,j,n,k)/2.        
   
         if(rosSNo(i,j,n,k)>=roCdSV+70) nisSNo(i,j,n)  = k  
 
        end do  ! k = 1,nbr_layer

        ! Antarctica
        if (INI==5) then
          if (rosSNo_ini(i,j).eq.470) firstPoint = firstPoint+1
          if(firstPoint.eq.1.and.verboseAntarctica) then
            print*, "[Ant] snowpack initialization"
          endif
          ! surface layer
          rosSNo(i,j,n,nbr_layer) = rosSNo_ini(i,j)
          deltaPa =  0.0  + 
     .               rosSNo_ini(i,j) *
     .               dzsSNo(i,j,n,nbr_layer) / 2. *
     .               gravit * 1e-6  ! overburden pressure (MPa)
          ! lower layers
          if(nbr_layer>1) then
            ro_w = 1000            ! water density (kg m-3)
            ro_i = 917             ! ice density   (kg m-3)
            R_gc = 8.3144621       ! gas constant  (J mol-1 K-1)
            E_1  = 10160           ! activation energy   (J mol-1)
            E_2  = 60000           ! activation energy   (J mol-1)
            Ao   = 2.54*10.**4     ! constante           (MPa-3 s-1)
            tas  = tasSNo_ini(i,j) ! mean surface temperature (K)
            do k = nbr_layer-1,1,-1
              ro_up = rosSNo(i,j,n,k+1)
              dz = (dzsSNo(i,j,n,k+1) + dzsSNo(i,j,n,k)) / 2.
              drhodz = -9999
              deltaPa = deltaPa + ro_up * dz * gravit * 10.**-6
              if(firstPoint.eq.1.and.verboseAntarctica) then
                 print*, "[Ant] i,j,n,k: ",i,j,n,k
              endif
              if (ro_up<=550.) then
                ! Ligtenberg (2011) > Barnola (1991)
                drhodz = 11 *                 ! (m-1)
     .                   ro_up / ro_w *       ! (s.u.)
     .                   (ro_i-ro_up) *       ! (rho)
     .                   exp(-E_1/(R_gc*tas)) ! (s.u.)
                if(firstPoint.eq.1.and.verboseAntarctica) then
                  print*,"[Ant] drdz  : ",drhodz
                  print*,"[Ant] drdz 1: ",11*ro_up*(ro_i-ro_up)
                  print*,"[Ant] drdz 2: ",E_1/(R_gc*tas)
                  print*,"[Ant] drdz 3: ",exp(-E_1/(R_gc*tas))
                endif
              else
                ! Ligtenberg (2011) > Barnola (1991)
                fe_ro = 10.**
     .                 (-37.455 * (ro_up / 1000.)**3 
     .                  +99.743 * (ro_up / 1000.)**2 
     .                  -95.027 * (ro_up / 1000.)    
     .                  +30.673                      )
                drhodz =  Ao *                          !  (MPa-3 s-1)
     .              deltaPa**3 *                        ! *(MPa+3)
     .              ro_up /                             ! *(kg m-3)
     .              (smbSNo_ini(i,j) / (365*24*3600)) * ! /(kg m-2 s-1)
     .              ro_up *                             ! (rho)
     .              fe_ro *                             ! (s.u).
     .              exp(-E_2 / (R_gc * tas))            ! (s.u.)
                if(firstPoint.eq.1.and.verboseAntarctica) then
                  print*,"[Ant] Ao*ro_up: ",Ao*ro_up
                  print*,"[Ant] smb: ",smbSNo_ini(i,j)/(365*24*3600)
                  print*,"[Ant] fe: ",fe_ro
                  print*,"[Ant] dPa: ",deltaPa
                  print*,"[Ant] ddPa: ",ro_up*dz*gravit*1e-6
                  print*,"[Ant] exp: ",exp(-E_2/(R_gc*tas))
                endif
              endif  ! ro_up<=550.
              rosSNo(i,j,n,k) = rosSNo(i,j,n,k+1) + drhodz * dz
              if(firstPoint.eq.1.and.verboseAntarctica) then
                 print*,"[Ant] ro_up: ",ro_up
                 print*,"[Ant] dz: ",dz
                 print*,"[Ant] tas: ",tas
                 print*,"[Ant] drhodz: ",drhodz
                 print*,"[Ant] rosSNo: ",rosSNo(i,j,n,k)
              endif
            enddo  !  k = nbr_layer-1,1,-1
            do k = 1,nbr_layer
               g1sSNo(i,j,n,k) = 99.
               g2sSNo(i,j,n,k) = g2s
               rosSNo(i,j,n,k) = max(200.,rosSNo(i,j,n,k))
               rosSNo(i,j,n,k) = min(600.,rosSNo(i,j,n,k))
               if(rosSNo(i,j,n,k)>=roCdSV) nisSNo(i,j,n)  = k
            enddo
          endif  ! nbr_layer>1
        endif  ! INI==5 (Antarctica)

        endif  ! if (nbr_layer> 0 )


         nssSNo(i,j,n)      =  nbr_layer
         issSNo(i,j,n)      =  0
         snohSN(i,j,n)      =  0  
         SWaSNo(i,j,n)      =  0.
         TsrfSL(i,j,n)      =  tairSL(i,j) 
         TvegTV(i,j,n)      =  tairSL(i,j)

         if (isolSL(i,j)    == 3) then 

          isolTV(i,j)       =  12
          iwafTV(i,j)       =  0 
          ivegTV(i,j,n)     =  0  
          alaiTV(i,j,n)     =  0. 
          glf_TV(i,j,n)     =  0. 
c          SLsrfl(i,j,1)     =  1.  
c          ifraTV(i,j,1)     =  100

          AlbSTV(i,j)       =  0.55

          do isn = 1,llx                                              
          TsolTV(i,j,n,isn) =  max(230.,tisSNo(i,j,n,1))
          eta_TV(i,j,n,isn) =  0.                                    
          end do     

         else

          AlbSTV(i,j)       =  0.25 ! albedo / 2 in SISVAT if satured soil
          do isn = 1,llx
           eta_TV(i,j,n,isn) =  max(0.01  ,eta_TV(i,j,n,isn)) 
           if (nbr_layer > 0) then
           TsolTV(i,j,n,isn) =  max(230.,tisSNo(i,j,n,1))
           eta_TV(i,j,n,isn) =  0.
           endif
          enddo

         end if

        END DO

       if (INI==3) then 
        do n=1,nsx-1
        nssSNo(i,j,n)=ns1(i,j)
        nisSNo(i,j,n)=ni1(i,j)
        do k = 1, nsno
         dzsSNo(i,j,n,k) = dz1(i,j,k)
         nhsSNo(i,j,n,k) = nh1(i,j,k)
         g1sSNo(i,j,n,k) = g11(i,j,k)
         g2sSNo(i,j,n,k) = g21(i,j,k)
         tisSNo(i,j,n,k) = ti1(i,j,k)
         rosSNo(i,j,n,k) = ro1(i,j,k)
         wasSNo(i,j,n,k) = wa1(i,j,k)
        enddo
        enddo         
       endif

      endif ; end do ; end do

C +--Initial snow height
C +  ===================

       do j=1,my ; do i=1,mx ;  do k=1,nvx

         maskSN(i,j) = mskSNo(i,j,1)

         zn0IB(i,j,k)    = 0.0
         mb0IB(i,j,k)    = 0.0 
                    
        if (nssSNo(i,j,k) > 0) then 
             
          znsn(i,j,nsno) = dzsSNo(i,j,k,nsno)
         snwae(i,j,nsno) = rosSNo(i,j,k,nsno)   
     .                   * dzsSNo(i,j,k,nsno) *1.d3
     .                    / ro_Wat
         do nk=nsno-1,1,-1
           znsn(i,j,nk)  = dzsSNo(i,j,k,nk) +   znsn(i,j,nk+1)
          snwae(i,j,nk)  = rosSNo(i,j,k,nk) * dzsSNo(i,j,k,nk)  *1.d3
     .                   / ro_Wat           +  snwae(i,j,nk+1)
         end do
         smbalSN0(i,j)   = snwae (i,j,1)    
     .                   - snwae (i,j,nssSNo(i,j,1)+1)
         znSNow0 (i,j)   = znsn(i,j,1)    
     .                   - znsn(i,j,nssSNo(i,j,1)+1)

         zn0IB(i,j,k)    = znSNow0 (i,j)
         mb0IB(i,j,k)    = smbalSN0(i,j)

       end if

       if (nssSNo(i,j,1) >= 10) then 
        ni = i    ; nj = j 
       endif

      end do ; end do ; enddo

C +--Output
C +  ======
       
      open(unit=111, status='replace', file='INIsnow.out')

      do nk = 1, 1
                 
       write(111,*) 'Coord:', ni,nj
       write(111,401) 
      
       write(111,402)(n,znsn(ni,nj,  n)     ,
     .                dzsSNo(ni,nj,1,n)*1000,tisSNo(ni,nj,1,n),
     .                rosSNo(ni,nj,1,n)     ,wasSNo(ni,nj,1,n),
     .                snwae (ni,nj,  n)     ,agsSNo(ni,nj,1,n), 
     .                zero,       zero      ,g1sSNo(ni,nj,1,n),
     .                g2sSNo(ni,nj,1,n)     ,nhsSNo(ni,nj,1,n),
     .                n=nsno,1,-1)                 
       write(111,*)  'mb0IB    :', mb0IB   (ni,nj,1)
       write(111,*)  'zn0IB    :', zn0IB   (ni,nj,1)
       write(111,*)  'nssSNo   :', nssSNo  (ni,nj,1)       
       write(111,*)  'nisSNo   :', nisSNo  (ni,nj,1) 
       write(111,*)  'SH       :', sh      (ni,nj)           
       write(111,*)  'tairSL   :', tairSL  (ni,nj)       
       write(111,*)  'tsrfSL   :', tsrfSL  (ni,nj,1)
       write(111,*)  't2_SL    :', t2_SL   (ni,nj)
       write(111,*)  'd1_SL    :', d1_SL   (ni,nj)
       write(111,*)  'SL_z0    :', SL_z0   (ni,nj,1)
       write(111,*)  'SL_r0    :', SL_r0   (ni,nj,1) 
       write(111,*)  'SLuusl   :', SLuusl  (ni,nj,1)
       write(111,*)  'SLutsl   :', SLutsl  (ni,nj,1)                   
       write(111,*)  'eps0SL   :', eps0SL  (ni,nj) 
      end do 
                        
      close(111)

 401  format(/,' Internal Characteristics',
     .        /,' ========================',
     .        /,'  n |  z    |  dz   |   T    | rho   |  W    |',
     .       ' z(WE) | Age   | Extin |  UW   | Dendr.| Spher.| Hist. |',
     .        /,'    | [m]   | [mm]  |  [K]   | kg/m3 | kg/kg |',
     .       '  [mm] | [d]   |       | mim/s | /Sphe.| /Size |       |',
     .        /,'----+-------+-------+--------+-------+-------+',
     .       '-------+-------+-------+-------+-------+-------+-------+')
 402  format((i3,' |',f6.2,' |',  f6.1,' |', f7.2,' |',  f6.1,' |',
     .      f6.3,' |',f7.0, '|',  f6.1,' |', f6.3,' |',  f6.2,' |',
     .    2(f6.1,' |'),i4,'   |'))

      return
      end
      
C +------------------------------------------------------------------------+
      subroutine UPDsnow
C +------------------------------------------------------------------------+
C | MAR SURFACE  XF                                                    MAR |
C +------------------------------------------------------------------------+
C +
      implicit none
C +
C +--General Variables
C +  =================

      include 'MARphy.inc'
      include 'MARCTR.inc'
      include 'MAR_SV.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc' 
      include 'MAR_GE.inc'
      include 'MAR_DY.inc'
      include 'MAR_LB.inc'
      include 'MAR_SL.inc'
      include 'MAR_SN.inc'
      include 'MAR_BS.inc'
      include 'MAR_IO.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MAR_IB.inc'

C +--Local   Variables
C +  =================

      integer                 n
 
      real                    dz_old,zn_old,mb_old,ro_old,dz_tot

      real      ,parameter :: maxlimit =21.0  ! Maximum Snow Height
      real      ,parameter :: minlimit =19.0  ! Minimum Snow Height      
      real      ,parameter ::   change = 0.01 ! percentage change (0-1)
 
C +   checking of zn1IB
C +   -----------------
             
      do i=1,mx ; do j=1,my ; do n=1,nsx

       dz_tot= 0.

       DO k=nsno,1,-1
        dz_tot = dzsSNo(i,j,n,k)   + dz_tot
       ENDDO

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       if ( dz_tot > maxlimit ) then
       
        k               = 1

        zn_old          = zn0IB (i,j,n)
        dz_old          = dzsSNo(i,j,n,k)
        mb_old          = dzsSNo(i,j,n,k)*rosSNo(i,j,n,k)
        ro_old          = rosSNo(i,j,n,k)

        if(dz_old>1.01) then
         dzsSNo(i,j,n,k)= dzsSNo(i,j,n,k) - 1.
         rosSNo(i,j,n,k)= max(rosSNo(i,j,n,k+1),rosSNo(i,j,n,k)*
     .                    (1.-change))
        else
         dzsSNo(i,j,n,k)= 0.1*dzsSNo(i,j,n,k)
         rosSNo(i,j,n,k)= rosSNo(i,j,n,k+1)
        endif  

        zn0IB (i,j,n)   = zn0IB (i,j,n)   +
     .                   (dzsSNo(i,j,n,k) - dz_old)
        mb0IB (i,j,n)   = mb0IB (i,j,n)   + rosSNo(i,j,n,k)
     .                  * dzsSNo(i,j,n,k) - mb_old
        wet0IB(i,j,n)   = wet0IB(i,j,n)   + rosSNo(i,j,n,k)
     .                  * dzsSNo(i,j,n,k) - mb_old


        write(*,*)  ' ' 
        write(*,10) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE,i,j,
     .              dz_tot,dz_tot-(zn_old-zn0IB(i,j,n)),
     .              ro_old,rosSNo(i,j,n,k)
   10   format(' UPDsnow (min) at',i5,4i3,
     .         ' for (',i3,','i3,') : ZN ',f5.2,'=>',f5.2,
     .                              ' RO ',f8.1,'=>',f8.1)
        write(*,*)  ' '

       end if

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       if (mskSNo(i,j,n)>0.and.n/=nsx .and. dz_tot > minlimit/2. !*CL* 
     .                                .and. dz_tot < minlimit ) then  !*CL* 

        k               = 1

        zn_old          = zn0IB (i,j,n)

        dz_old          = dzsSNo(i,j,n,k)
        mb_old          = dzsSNo(i,j,n,k)*rosSNo(i,j,n,k)
        ro_old          = rosSNo(i,j,n,k)

        dzsSNo(i,j,n,k) = dzsSNo(i,j,n,k) + 1.
        tisSNo(i,j,n,k) = min(tisSNo(i,j,n,k),272.65)

        if(dz_tot > minlimit -1) 
     .  rosSNo(i,j,n,k) = rosSNo(i,j,n,k)*(1.+change)
        if(rosSNo(i,j,n,k)>880) 
     .  rosSNo(i,j,n,k) = ro_ice

        zn0IB (i,j,n)   = zn0IB (i,j,n)   +
     .                   (dzsSNo(i,j,n,k) - dz_old)
        mb0IB (i,j,n)   = mb0IB (i,j,n)   + rosSNo(i,j,n,k)
     .                  * dzsSNo(i,j,n,k) - mb_old
        wet0IB(i,j,n)   = wet0IB(i,j,n)   + rosSNo(i,j,n,k)
     .                  * dzsSNo(i,j,n,k) - mb_old

        write(*,*)  ' ' 
        write(*,11) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE,i,j,
     .              dz_tot,dz_tot+(zn0IB(i,j,n)-zn_old),
     .              ro_old,rosSNo(i,j,n,k)
   11   format(' UPDsnow (add) at',i5,4i3,
     .         ' for (',i3,','i3,') : ZN ',f5.2,'=>',f5.2,
     .                              ' RO ',f8.1,'=>',f8.1)
        write(*,*)  ' '
                 
       end if
      
      enddo; enddo ; enddo
      end subroutine             

C +------------------------------------------------------------------------+
      subroutine FILsnow
C +------------------------------------------------------------------------+
C | MAR SURFACE  XF                                                    MAR |
C +------------------------------------------------------------------------+
C +
      implicit none
C +
C +--General Variables
C +  =================

      include 'MARphy.inc'
      include 'MARCTR.inc'
      include 'MAR_SV.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc' 
      include 'MAR_GE.inc'
      include 'MAR_DY.inc'
      include 'MAR_LB.inc'
      include 'MAR_SL.inc'
      include 'MAR_SN.inc'
      include 'MAR_BS.inc'
      include 'MAR_IO.inc'
      include 'MAR_TV.inc'
      include 'MARsSN.inc'
      include 'MARdSV.inc'
      include 'MAR_IB.inc'

C +--Local   Variables
C +  =================


      integer                 n,l,filtering
 
      real                    ro_new,ww,g1_new,g2_new,ti_new,al_new
      real                    nbr1,nbr2      
               
      do i=2,mx-1 ; do j=2,my-1 ;do n=1,nsx-1   !*CL*
       if (mskSNo(i,j,n) >= 50) then    !*CL*

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
       nbr1=0   ; nbr2=0
       g1_new=0 ; g2_new=0 ; ro_new=0 ; ti_new=0 ; al_new=0

       do k=-1,1 ; do l=-1,1
        if(mskSNo(i+k,j+l,n) >= 50                     .and.
     .     rosSNo(i+k,j+l,n,max(1,nssSNo(i+k,j+l,n)))<roCdSV) then
                            ww=1     
         if (k==0 .or.l==0) ww=2
         if (k==0.and.l==0) ww=4

         ww     = ww*dzsSNo(i+k,j+l,n,nssSNo(i+k,j+l,n))

         if (g1sSNo(i+k,j+l,n,nssSNo(i+k,j+l,n))>=0) then 
         g1_new = g1_new+g1sSNo(i+k,j+l,n,nssSNo(i+k,j+l,n))*ww
         g2_new = g2_new+g2sSNo(i+k,j+l,n,nssSNo(i+k,j+l,n))*ww
         endif

         ro_new = ro_new+rosSNo(i+k,j+l,n,nssSNo(i+k,j+l,n))*ww
         ti_new = ti_new+tisSNo(i+k,j+l,n,nssSNo(i+k,j+l,n))*ww
         al_new = al_new+albxSL(i+k,j+l,n)                  *ww

         nbr1   = nbr1+ww
         nbr2   = nbr2+1

        endif
       enddo  ; enddo

       g1_new = g1_new/nbr1
       g2_new = g2_new/nbr1
       ro_new = ro_new/nbr1
       ti_new = ti_new/nbr1
       al_new = al_new/nbr1

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                                                    filtering=0

       if(nbr2==9                                 .and.
     .    g2sSNo(i,j,n,nssSNo(i,j,n))>1.1*g2_new  .and.
     .    rosSNo(i,j,n,nssSNo(i,j,n))<roCdSV      .and.
     .    albxSL(i,j,n)              <al_new-0.05 .and.
     .    g2_new                     <50.         .and.
     .    albxSL(i,j,n)              <0.70        .and.
     .    albxSL(i,j,n)              >0.50        .and.
     .    g2sSNo(i,j,n,nssSNo(i,j,n))<75.         ) filtering=1

c       if(nbr2==9.and.albxSL(i,j,n)<al_new-0.1    .and.
c     .    g2sSNo(i,j,n,nssSNo(i,j,n))<80          .and.
c     .    g1sSNo(i,j,n,nssSNo(i,j,n))<80          ) filtering=2

       if(filtering>0) then

        write(*,*)  ' '
        write(*,12) iyrrGE,mmarGE,jdarGE,jhurGE,minuGE,i,j,filtering,
     .              rosSNo(i,j,n,nssSNo(i,j,n)),ro_new,
     .              g1sSNo(i,j,n,nssSNo(i,j,n)),g1_new,
     .              g2sSNo(i,j,n,nssSNo(i,j,n)),g2_new,
     .              albxSL(i,j,n)              ,al_new
   12   format('Filtering',
     .           i5,4i3,' for (',i3,','i3,')',i2,f6.1,'=>',f6.1,
     .           f6.1,'=>',f6.1,',',f7.2,'=>',f7.2,f5.2,'=>',f5.2)
        write(*,*)

        g1sSNo(i,j,n,nssSNo(i,j,n)) = min(g1_new,
     .                                g1sSNo(i,j,n,nssSNo(i,j,n)))
        g2sSNo(i,j,n,nssSNo(i,j,n)) = min(g2_new,
     .                                g2sSNo(i,j,n,nssSNo(i,j,n)))
        tisSNo(i,j,n,nssSNo(i,j,n)) = min(ti_new,
     .                                tisSNo(i,j,n,nssSNo(i,j,n)))
        ro_new                      = min(ro_new,
     .                                rosSNo(i,j,n,nssSNo(i,j,n)))
        dzsSNo(i,j,n,nssSNo(i,j,n)) = dzsSNo(i,j,n,nssSNo(i,j,n))
     .                              * rosSNo(i,j,n,nssSNo(i,j,n))
     .                              / ro_new
        rosSNo(i,j,n,nssSNo(i,j,n)) = ro_new
       endif

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       if(mmarGE<=3.or.mmarGE>=11) then
        tisSNo(i,j,n,1)=min(273.14,tisSNo(i,j,n,1))
       endif

      endif; enddo ; enddo ; enddo

      end subroutine 


C +------------------------------------------------------------------------+
      subroutine OUTone(TypeGL) 
C +------------------------------------------------------------------------+
C | MAR OUTPUT  XF                                                     MAR |
C |                                                                        |
C +------------------------------------------------------------------------+

      IMPLICIT NONE
      
C +--Global Variables
C +  ================

      include 'MARphy.inc'
      include 'MARCTR.inc'
      include 'MARdim.inc'
      include 'MARgrd.inc'
      include 'MAR_DY.inc'
      include 'MAR_GE.inc'
      include 'MAR_SL.inc'
      include 'MAR_HY.inc'
      include 'MAR_RA.inc'
      include 'MAR_SN.inc'
      include 'MAR_SV.inc'
      include 'MARsSN.inc'
      include 'MAR_IB.inc'
      include 'MAR_IO.inc'
      include 'MAR_WK.inc'
      include 'MAR_CA.inc'
      include 'MAR_TV.inc'
      include 'MAR_TE.inc'
      include 'MAR_TU.inc'

C +--Local  Variables
C +  ================
C +
      CHARACTER*(40)         fnamNC_one
      common/OUT_nc_one_loc/ fnamNC_one
C +...                       fnamNC_one: To retain file name.

      integer    NdimNC_one
      PARAMETER (NdimNC_one = 8)
C +...Number of defined spatial dimensions (exact)

      integer    MXdim
      PARAMETER (MXdim = 9000)
C +...Maximum Number of all dims: recorded Time Steps
C +   and also maximum of spatial grid points for each direction. 

      integer    MX_var
      PARAMETER (MX_var = 120)
C +...Maximum Number of Variables 

      integer    NattNC_one
      PARAMETER (NattNC_one = 2)
C +...Number of real attributes given to all variables

c     ------------------------------------------------------------
      integer    , parameter :: ONEsta   = 60 ! nbr, station
cSVa  integer    , parameter :: ONEsta   = 19 ! nbr, station
      integer    , parameter :: ONElev   = mz ! nbr levels (<mz)
      integer    , parameter :: ONEint   = 30 ! interval (min)
c     ------------------------------------------------------------

      integer           x_one(ONEsta)      , y_one(ONEsta)     

      real              one1 (ONEsta)       , one2(ONEsta) 
      real              one3 (ONEsta)       , one4(ONEsta) 
      real              one5 (ONEsta)       , one6(ONEsta)       
      real              one7 (ONEsta,ONElev), one8(ONEsta,ONElev)
      real              one9 (ONEsta,ONElev),one10(ONEsta,ONElev)
      real              one11(ONEsta,  nsno),one12(ONEsta,  nsno)   
      real              one14(ONEsta,  llx) ,one15(ONEsta,  llx)
      real              one13(8) 
      real              uu   (ONElev)       ,vv(ONElev)
      real              WS   (ONEsta,ONElev),WD(ONEsta,ONElev)
      real              RH   (ONEsta,ONElev)  
      
      real              starti,starta 
      real              yearNC_one(MXdim)
      real              dateNC_one(MXdim)
      real              timeNC_one(MXdim)
      real              VALdim(MXdim,0:NdimNC_one)

      integer           njmo,ipr_nc_one
      integer           jourNC_one(MXdim)
      integer           moisNC_one(MXdim)
      integer           NvatNC_one(NattNC_one)      
      integer           nDFdim(      0:NdimNC_one)
          
      CHARACTER*(13)    NAMdim(      0:NdimNC_one)
      CHARACTER*(31)    UNIdim(      0:NdimNC_one)
      CHARACTER*(13)    SdimNC_one(4,MX_var)       
      CHARACTER*(31)    unitNC_one(MX_var)
      CHARACTER*(13)    nameNC_one(MX_var)
      CHARACTER*(50)    lnamNC_one(MX_var)
      CHARACTER*(100)   tit_nc_one
      CHARACTER*(13)    NAMrat(NattNC_one)
      CHARACTER*120     tmpINP
      CHARACTER*20      n_one(ONEsta)
      CHARACTER*2       station
      CHARACTER*10      TypeGL      

      integer           n1000 ,n100a ,n100,n10_a,n10,n1,m10,jd10,jd1
      integer           it    ,mois  ,mill  ,itotNC_one
      integer           NtotNC_one,ID__nc_one
      integer           nbr_day, nbr_output,dt2

      integer           ii,jj,kk,ll,s,n,mm,nn,one,imex,jmex,t       
      real              qsat0D,q,qst,r,rst,epsilon

      logical           first

      common/OUTone_l/  first 
      common/OUTone_i/  nDFdim,ipr_nc_one,x_one,y_one
      common/OUTone_r/  yearNC_one,dateNC_one,timeNC_one
      common/OUTone_c/  n_one

      real             ,parameter :: a    = 
     .                   6371.229 * 1000.0     ! radius of the Earth
      real             ,parameter :: conv = 
     .                   15.0*3.141592/180.0   ! Conversion 
                                               ! hour ==> rad   
 
C     1. Station Location Initialization
C     ==================================
                
      IF (.not. first .and. ipr_nc_one == 0 ) THEN
                first =.true.  
  
       print *, ' ' 
       print *, 'Begin of OUTone initialisation for '//trim(TypeGL)
       print *, '======================================'
       print *, ' '  

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
       IF (trim(TypeGL)=='GR25km') THEN
                                                                  s=  1
        x_one(s)=20  ; y_one(s)=63  ; n_one(s)="Aasiat"         ; s=s+1 ! 
        x_one(s)=28  ; y_one(s)=126 ; n_one(s)="Alert"          ; s=s+1 ! 
        x_one(s)=28  ; y_one(s)=55  ; n_one(s)="Aurora"         ; s=s+1 ! 
        x_one(s)=52  ; y_one(s)=57  ; n_one(s)="Aputiek"        ; s=s+1 ! 
        x_one(s)=43  ; y_one(s)=75  ; n_one(s)="Barber"         ; s=s+1 ! 
        x_one(s)=31  ; y_one(s)=67  ; n_one(s)="Crawford 1"     ; s=s+1 ! 
        x_one(s)=31  ; y_one(s)=66  ; n_one(s)="Crawford 2"     ; s=s+1 ! 
        x_one(s)=65  ; y_one(s)=91  ; n_one(s)="Daneborg"       ; s=s+1 ! 
        x_one(s)=62  ; y_one(s)=101 ; n_one(s)="Dove Bugt"      ; s=s+1 ! (63,101)
        x_one(s)=12  ; y_one(s)=103 ; n_one(s)="Dundas"         ; s=s+1 ! 
        x_one(s)=30  ; y_one(s)=52  ; n_one(s)="Dye-2"          ; s=s+1 ! 
        x_one(s)=34  ; y_one(s)=46  ; n_one(s)="Dye-3"          ; s=s+1 ! 
        x_one(s)=27  ; y_one(s)=67  ; n_one(s)="ETH-camp"       ; s=s+1 ! 
        x_one(s)=22  ; y_one(s)=56  ; n_one(s)="GIMEX-Mast 1"   ; s=s+1 !
        x_one(s)=23  ; y_one(s)=56  ; n_one(s)="GIMEX-Mast 2"   ; s=s+1 !
        x_one(s)=24  ; y_one(s)=56  ; n_one(s)="GIMEX-Mast 5"   ; s=s+1 !
        x_one(s)=25  ; y_one(s)=56  ; n_one(s)="GIMEX-Mast 6"   ; s=s+1 !
        x_one(s)=26  ; y_one(s)=56  ; n_one(s)="GIMEX-Mast 6-9" ; s=s+1 !
        x_one(s)=27  ; y_one(s)=56  ; n_one(s)="GIMEX-Mast 9"   ; s=s+1 ! 
        x_one(s)=21  ; y_one(s)=104 ; n_one(s)="Gits"           ; s=s+1 ! (21,104)
        x_one(s)=43  ; y_one(s)=79  ; n_one(s)="Gisp2 (Summit)" ; s=s+1 !
        x_one(s)=44  ; y_one(s)=79  ; n_one(s)="Grip"           ; s=s+1 !
        x_one(s)=28  ; y_one(s)=121 ; n_one(s)="Hall Land"      ; s=s+1 ! 
        x_one(s)=27  ; y_one(s)=108 ; n_one(s)="Humboldt"       ; s=s+1 ! 
        x_one(s)=36  ; y_one(s)=31  ; n_one(s)="Ikermiuarsuk"   ; s=s+1 ! 
        x_one(s)=24  ; y_one(s)=66  ; n_one(s)="Ilulissat"      ; s=s+1 ! 
        x_one(s)=26  ; y_one(s)=67  ; n_one(s)="Jar1"           ; s=s+1 ! 
        x_one(s)=25  ; y_one(s)=67  ; n_one(s)="Jar2"           ; s=s+1 ! 
        x_one(s)=25  ; y_one(s)=66  ; n_one(s)="Jar3"           ; s=s+1 ! 
        x_one(s)=48  ; y_one(s)=79  ; n_one(s)="Julie"          ; s=s+1 ! 
        x_one(s)=23  ; y_one(s)=56  ; n_one(s)="Kangerdlugssuaq"; s=s+1 !
        x_one(s)=22  ; y_one(s)=30  ; n_one(s)="Kangilinnguit"  ; s=s+1 ! 
        x_one(s)=44  ; y_one(s)=128 ; n_one(s)="Kap M. Jesup "  ; s=s+1 ! 
        x_one(s)=51  ; y_one(s)=66  ; n_one(s)="Kar"            ; s=s+1 ! 
        x_one(s)=43  ; y_one(s)=78  ; n_one(s)="Kenton"         ; s=s+1 ! 
        x_one(s)=40  ; y_one(s)=78  ; n_one(s)="Klinck (Mast10)"; s=s+1 !
        x_one(s)=42  ; y_one(s)=49  ; n_one(s)="Kulu"           ; s=s+1 !
        x_one(s)=52  ; y_one(s)=91  ; n_one(s)="Nasa-E"         ; s=s+1 !
        x_one(s)=36  ; y_one(s)=53  ; n_one(s)="Nasa-SE"        ; s=s+1 !
        x_one(s)=29  ; y_one(s)=86  ; n_one(s)="Nasa-U"         ; s=s+1 !
        x_one(s)=44  ; y_one(s)=84  ; n_one(s)="Matt"           ; s=s+1 !
        x_one(s)=29  ; y_one(s)=29  ; n_one(s)="Narssarssuaq"   ; s=s+1 !
        x_one(s)=30  ; y_one(s)=102 ; n_one(s)="NEEM1"          ; s=s+1 !
        x_one(s)=31  ; y_one(s)=102 ; n_one(s)="NEEM2"          ; s=s+1 !
        x_one(s)=38  ; y_one(s)=90  ; n_one(s)="Ngrip"          ; s=s+1 !
        x_one(s)=56  ; y_one(s)=122 ; n_one(s)="Nord"           ; s=s+1 !
        x_one(s)=19  ; y_one(s)=44  ; n_one(s)="Nuuk (Godthab)" ; s=s+1 ! 
        x_one(s)=21  ; y_one(s)=34  ; n_one(s)="Paamiut"        ; s=s+1 !
        x_one(s)=34  ; y_one(s)=24  ; n_one(s)="Prins Ch"       ; s=s+1 !
        x_one(s)=29  ; y_one(s)=27  ; n_one(s)="Qaqortoq"       ; s=s+1 !
        x_one(s)=33  ; y_one(s)=50  ; n_one(s)="Saddle"         ; s=s+1 !
        x_one(s)=18  ; y_one(s)=48  ; n_one(s)="Sioralik"       ; s=s+1 ! Maniitsoq   
        x_one(s)=18  ; y_one(s)=55  ; n_one(s)="Sisimiut"       ; s=s+1 !
        x_one(s)=31  ; y_one(s)=38  ; n_one(s)="South Dome"     ; s=s+1 !
        x_one(s)=44  ; y_one(s)=78  ; n_one(s)="Summit (Cathy)" ; s=s+1 !
        x_one(s)=46  ; y_one(s)=49  ; n_one(s)="Tasiilaq"       ; s=s+1 !
        x_one(s)=13  ; y_one(s)=104 ; n_one(s)="Thules"         ; s=s+1 !
        x_one(s)=47  ; y_one(s)=104 ; n_one(s)="Tunu-N"         ; s=s+1 !
        x_one(s)=21  ; y_one(s)=82  ; n_one(s)="Upernavik"      ; s=s+1 !
        x_one(s)=67  ; y_one(s)=71  ; n_one(s)="Uunartoq"       ; s=s+1 !
          one=s-1      

        DO ll=1,ONEsta
         y_one(ll)=max(1,y_one(ll)-4) ! WARMING (X,Y) of the stations made for 80x140
        ENDDO

       END IF

       IF (trim(TypeGL)=='SVa') THEN
        n_one(1)='Svalbard Lufthavn'
        x_one(1)=36 ; y_one(1)=44

        n_one(2)='Hopen'
        x_one(2)=60 ; y_one(2)=26

        n_one(3)='NY-Alesund'
        x_one(3)=27 ; y_one(3)=52

        n_one(4)='Sveagruva'
        x_one(4)=36 ; y_one(4)=40

        n_one(5)='Barentsburg'
        x_one(5)=32 ; y_one(5)=42

        n_one(6)='Adventalen'
        x_one(6)=36 ; y_one(6)=44

        n_one(7)='Breinosa (Kho Aurora Station)'
        x_one(7)=36 ; y_one(7)=40

        n_one(8)='Janssonhaugen'
        x_one(8)=36 ; y_one(8)=40

        n_one(9)='Gruvefjellet'
        x_one(9)=34 ; y_one(9)=40

        n_one(10)='Kapp Lee'
        x_one(10)=48 ; y_one(10)=42

        n_one(11)='Rijpfjorden'
        x_one(11)=45 ; y_one(11)=39

        n_one(12)='Kvitoya'
        x_one(12)=67 ; y_one(12)=68

        n_one(13)='Kongsoya'
        x_one(13)=65 ; y_one(13)=53

        n_one(14)='Verlegenhuken'
        x_one(14)=38 ; y_one(14)=64

        n_one(15)='Edgeoya-Kapp Heuglin'
        x_one(15)=53 ; y_one(15)=44

        n_one(16)='Svarttangen'
        x_one(16)=49 ; y_one(16)=36

        n_one(17)='Isfjord Radio'
        x_one(17)=32 ; y_one(17)=42

        n_one(18)='Hornsund'
        x_one(18)=35 ; y_one(18)=30

        n_one(19)='Svalbard LH - Plataberget' 
        x_one(19)=36 ; y_one(19)=43
          one=19
       END IF

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

       if (one/=ONEsta) then
        print *, "ERROR: incorrect ONEsta",one ; stop 
       end if
       
       print *, ' '
       print *, 'Weather Stations MAR coordinates'
       print *, '       for Simulation '//trim(TypeGL)  ; print *, ' ' 
       print *, ' N°  ii  jj Name' 
       print *, '------------------------------'

       DO ll=1,ONEsta
        x_one(ll) = max(1, min(x_one(ll),mx))
        y_one(ll) = max(1, min(y_one(ll),my))
        write(*,'(3i4,x,a20)') ll,x_one(ll),y_one(ll),n_one(ll)        
       END DO
       
       print *, '------------------------------'     

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

C     2. NetCDF File Initialization
C     =============================

          n1000 = 1 +     iyrrGE/1000
          n100a =     mod(iyrrGE,1000)
          n100  = 1 +     n100a /100
          n10_a =     mod(n100a ,100)
          n10   = 1 +     n10_a /10
          n1    = 1 + mod(n10_a ,10)
          m10   = 1 +     mmarGE/10
          m1    = 1 + mod(mmarGE,10)
          jd10  = 1 +     jdarGE/10
          jd1   = 1 + mod(jdarGE,10)

C +   2.1 Output File Label
C +   ---------------------

       fnamNC_one = 'ONE.'
     .        // labnum(n1000) // labnum(n100)
     .        // labnum(  n10) // labnum(  n1)
     .        // labnum(  m10) // labnum(  m1)
     .        // labnum( jd10) // labnum( jd1)
     .        // '.' // explIO
     .        // '.nc    '

C +   2.2 Output Title
C +   ----------------

       tit_nc_one = 'ONE'
     .        // ' - Exp: ' // explIO
     .        // ' - '
     .        // labnum(n1000) // labnum(n100)
     .        // labnum(  n10) // labnum(  n1)
     .        // labnum(  m10) // labnum(  m1)
     .        // labnum( jd10) // labnum( jd1)

C +   2.3 Time Variable (hour)
C +   ------------------------

C      ++++++++++++++++++++++++++++
       ipr_nc_one = 0
       nbr_output = int(real(max(1,int(nterun*dt/86400)))
     .            * 24.0*60.0/real(ONEint)) 
       nDFdim(0)  = nbr_output
       nDFdim(0)  = 0             ! Unlimited
C      ++++++++++++++++++++++++++++	

       NAMdim(0)  = 'time'
       UNIdim(0)  = 'HOURS since 1901-01-15 00:00:00'

       IF (nbr_output.gt.MXdim)
     & STOP '*** OUT_nc_one - ERROR : MXdim to low ***'

       starti = jhurGE + minuGE/  60.0               ! Starting Time
     .                 + jsecGE/3600.0 - dt / 3600.0 

       starta = (351   +(iyrrGE  -1902) *365         ! Nb Days before iyrrGE
     .                 +(iyrrGE  -1901) /  4         ! Nb Leap Years
     .                 + njyrGE(mmarGE)              ! Nb Days before mmarGE
     .                 + njybGE(mmarGE)              ! (including Leap Day)
     .                 * max(0,1-mod(iyrrGE,4))      !
     .                 + jdarGE     -1      )*  24   !
     .                 + jhurGE                      !
     .                 + (minuGE *60 + jsecGE - dt)/3600.

       DO it = 1,nbr_output
             timeNC_one(it)   = starti + (it-1) * ONEint / 60.
             VALdim(it,0)     = starta + (it-1) * ONEint / 60.
             dateNC_one(it)   =          timeNC_one(it)
             jourNC_one(it)   = jdarGE + timeNC_one(it) / 24.
       END DO
                 mois       =  mmarGE
                 mill       =  iyrrGE
       DO it = 1,nbr_output
        IF     (mois      .eq.2           .AND.
     .      mod(mill,4)   .eq.0           )                      THEN 
                 njmo       =  njmoGE(mois) + 1
         ELSE                                                        
                 njmo       =  njmoGE(mois)
         END IF
  
         IF     (jourNC_one(it).gt.njmo        )                 THEN 
           DO t=it,nbr_output
                 jourNC_one(t) =  jourNC_one(t) - njmo
           END DO
                 mois       =  mois + 1
             IF (mois.gt.12)                                     THEN 
                 mois       =         1
                 mill       =  mill + 1
             END IF                                                   
         END IF                                                      
                 moisNC_one(it) =  mois
                 yearNC_one(it) =  mill

         IF     (dateNC_one(it).gt.24.-epsi)                     THEN 
           DO t=it,nbr_output
                 dateNC_one(t) = mod(dateNC_one(t),24.)
           END DO
         END IF                                                      
       END DO

       DO it = 1,nbr_output
              dateNC_one(it) =        dateNC_one(it)
     .                       + 1.d+2 *jourNC_one(it) 
     .                       + 1.d+4 *moisNC_one(it) 
     .                       + 1.d+6 *yearNC_one(it) 
       END DO

C +   2.4 Define horizontal spatial dimensions   
C +   ----------------------------------------

       DO i = 1, mx
        VALdim(i,1) = xxkm(i)
       END DO
       nDFdim(1)= mx     ; NAMdim(1)= 'x'       ; UNIdim(1)= 'km'
      
       DO j = 1, my
        VALdim(j,2) = yykm(j)
       END DO
       nDFdim(2)= my     ; NAMdim(2)= 'y'       ; UNIdim(2)= 'km'

       do k = 1, ONElev
        VALdim(k,3) = sigma(mz-k+1)
       enddo
       nDFdim(3)= ONElev ; NAMdim(3)= 'level'   ; UNIdim(3)= '[sigma]'

       do k = 1, ONEsta
        VALdim(k,4) = k
       enddo
       nDFdim(4)= ONEsta ; NAMdim(4)= 'station' ; UNIdim(4)= '[-]'

       DO k = 1, nsx
         VALdim(k,5) = k
       END DO
       nDFdim(5)= nsx    ; NAMdim(5)= 'sector'  ; UNIdim(5)= '[index]'

       DO k = 1, nsno
         VALdim(k,6) = k
       END DO
       nDFdim(6)= nsno   ; NAMdim(6)= 'snolay'  ; UNIdim(6)= '[index]'

       DO k = 1, nsol+1
         VALdim(k,7) = k
       END DO
       nDFdim(7)= nsol+1 ; NAMdim(7)= 'sollay'  ; UNIdim(7)= '[-]'

       DO k = 1, 8
         VALdim(k,8) = k
       END DO
       nDFdim(8)= 8      ; NAMdim(8)= 'info'    ; UNIdim(8)= '[-]'


C +   2.5 Variable's Choice (Table ONEvou.dat)
C +   ----------------------------------------

       DO ll=1,ONEsta   
     
        if (LL<=9)then
         write(station,'(i1)') ll
        else
         write(station,'(i2)') ll
        endif

         nameNC_one  (ll)="S"//station
         SdimNC_one(1,ll)="info"
         SdimNC_one(2,ll)="-"
         SdimNC_one(3,ll)="-"
         SdimNC_one(4,ll)="-"
         unitNC_one  (ll)="-"
         lnamNC_one  (ll)=n_one(ll)

       ENDDO

       OPEN(unit=10,status='old',file='ONEvou.dat')

       itotNC_one = ONEsta
 980   CONTINUE
         READ (10,'(A120)',end=990) tmpINP
         IF (tmpINP(1:4).eq.'    ')                                THEN 
           itotNC_one = itotNC_one + 1
           READ (tmpINP,'(4x,5A9,A12,A50)')
     .          nameNC_one(itotNC_one),  
     .          SdimNC_one(1,itotNC_one),
     .          SdimNC_one(2,itotNC_one),
     .          SdimNC_one(3,itotNC_one),
     .          SdimNC_one(4,itotNC_one),
     .          unitNC_one  (itotNC_one),
     .          lnamNC_one  (itotNC_one)
         ENDIF
       GOTO 980
 990   CONTINUE

       CLOSE(unit=10)

       NtotNC_one = itotNC_one 
C +... NtotNC_one : Total number of variables writen in NetCDF file.

C +   2.6 List of NetCDF attributes given to all variables
C +   ----------------------------------------------------

       NAMrat(1) = 'actual_range'
       NvatNC_one(1) = 2

       NAMrat(NattNC_one) = '[var]_range'
       NvatNC_one(NattNC_one) = 2

C +   2.7 Automatic Generation of the NetCDF File Structure
C +   -----------------------------------------------------

C +    **************
       CALL UNscreate (fnamNC_one,tit_nc_one,
     &                 NdimNC_one, nDFdim, MXdim , NAMdim, UNIdim, 
     &                 VALdim,
     &                 MX_var, NtotNC_one, nameNC_one, SdimNC_one, 
     &                 unitNC_one, 
     &                 lnamNC_one,
     &                 NattNC_one, NAMrat, NvatNC_one,
     &                 ID__nc_one) 
C +    **************

C +   2.8 Write Time - Constants
C +   --------------------------

       imex = int(mx/2)
       jmex = int(my/2)

       Wkxy1  =  GElonh * 15.  ! Hour->degrees
       WKxy2  =  GElatr / degrad ! rad ->degree
       WKxy3  =  real(isolSL)             ! REAL type
       WKxy4(:,:)  =  real(mskSNo(:,:,1)) ! REAL type

C +    ************
       CALL UNwrite (ID__nc_one,'LON', 1,    mx, my, 1, Wkxy1)
       CALL UNwrite (ID__nc_one,'LAT', 1,    mx, my, 1, Wkxy2)
       CALL UNwrite (ID__nc_one,'SH' , 1,    mx, my, 1, sh)
       CALL UNwrite (ID__nc_one,'SRF', 1,    mx, my, 1, Wkxy3)
       CALL UNwrite (ID__nc_one,'SLO', 1,    mx, my, 1, slopTV)
       CALL UNwrite (ID__nc_one,'MSK', 1,    mx, my, 1, WKxy4)
C +    ************

       open(unit=1000,status='replace',file='OUTone.jnl')  

       do ll = 1,ONEsta 

        one13(1) = x_one(ll)
        one13(2) = y_one(ll)
        one13(3) = Wkxy1(x_one(ll),y_one(ll))
        one13(4) = Wkxy2(x_one(ll),y_one(ll))
        one13(5) = sh   (x_one(ll),y_one(ll))
        one13(6) = Wkxy4(x_one(ll),y_one(ll))
        one13(7) = dx*(x_one(ll)-imez)*1.e-3
        one13(8) = dx*(y_one(ll)-jmez)*1.e-3

        if (LL<=9)then
         write(station,'(i1)') ll
        else
         write(station,'(i2)') ll
        endif
         
        CALL UNwrite (ID__nc_one,"S"//station, 1, 8,  1, 1, one13)

        write(1000,1001) one13(7),one13(8),n_one(ll)
 1001   format('LABEL ',2(f8.2,','),'-1,0,.10 @TR+',a20)     

       end do

C +   2.9 Re-Open file if already created.
C +   -----------------------------------

       write (*,*) ' ' 
       write (*,*) 'End of OUTone initialisation'
       write (*,*) '============================'
       write (*,*) ' ' 

       GOTO 1000

      ELSE

C +    ************
       CALL UNwopen (fnamNC_one,ID__nc_one)
C +    ************

      END IF

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

C     3. Write Time-dependent variables
C     =================================

      ipr_nc_one = ipr_nc_one + 1 ! not use with Fortran 77     

c      write(*,*)  ' '
c      write(*,10) ipr_nc_one,iyrrGE,mmarGE,jdarGE,jhurGE,minuGE
c   10 format(' OUTone call (',i4,'):',i5,4i3)
c      write(*,*)  

      IF (nDFdim(0).eq.0)                         THEN !
                                     dt2 = dt
              if(iterun>=nterun -1)  dt2 = 0
              starta = (351+(iyrrGE  -1902) *365       ! Nb Days before iyrrGE
     .                     +(iyrrGE  -1901) /  4       ! Nb Leap Years
     .                     + njyrGE(mmarGE)            ! Nb Days before mmarGE
     .                     + njybGE(mmarGE)            ! (including Leap Day)
     .                 *max(0,1-mod(iyrrGE,4))         !
     .                     + jdarGE     -1      )*  24 !
     .                 +jhurGE                         !
     .               + (minuGE *60 +jsecGE -dt2 )/3600.!
C +     ************
        CALL UNwrite (ID__nc_one, 'time',ipr_nc_one, 1, 1, 1, starta)
C +     ************
      END IF

C +     ************
        CALL UNwrite (ID__nc_one, 'DATE',ipr_nc_one, 1, 1, 1, 
     .                dateNC_one(ipr_nc_one))
        CALL UNwrite (ID__nc_one, 'year',ipr_nc_one, 1, 1, 1, 
     .                yearNC_one(ipr_nc_one))
C +     ************

C +     ************
        starta=iyrrGE
        CALL UNwrite (ID__nc_one, 'YYYY',ipr_nc_one, 1, 1, 1,starta)
        starta=mmarGE
        CALL UNwrite (ID__nc_one, 'MM',ipr_nc_one, 1, 1, 1,starta)
        starta=jdarGE
        CALL UNwrite (ID__nc_one, 'DD',ipr_nc_one, 1, 1, 1,starta)
        starta=jhurGE
        CALL UNwrite (ID__nc_one, 'HH',ipr_nc_one, 1, 1, 1,starta)
        starta=minuGE
        CALL UNwrite (ID__nc_one, 'MIN',ipr_nc_one, 1, 1, 1,starta)
        starta=jsecGE
        CALL UNwrite (ID__nc_one, 'SS',ipr_nc_one, 1, 1, 1,starta)
C +     ************

        if (ipr_nc_one.eq.1)
     .  CALL UNwrite (ID__nc_one,'SLO', 1,    mx, my, 1, slopTV)

      Do ll=1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll)

C +   3.1 Computation of Relative Humidity
C +   ------------------------------------

       DO kk     = 1, ONElev

       epsilon   = 0.622
       q         = qvDY(ii,jj,mz-kk+1)
       qst       = qsat0D(tairDY(ii,jj,mz-kk+1),
     .                  sigma(mz-kk+1),pstDY(ii,jj),ptopDY,1)
          
       r         = q   / max(epsi,1.-q)
       rst       = qst / max(epsi,1.-qst)

       RH(ll,kk) =  (r/(epsilon+r)) 
     .           / max(epsi,(rst/(epsilon+rst))) * 100.

       RH(ll,kk) = max(0.,min(100.,RH(ll,kk)))

C +   3.2 Computation of wind direction and wind speed 
C +   ------------------------------------------------

        IF(uairDY(ii,jj,mz-kk+1)/=0.0  .and.
     .     vairDY(ii,jj,mz-kk+1)/=0.0) THEN

        uu(kk)   = (GElonh(ii+1,jj) - GElonh(ii,jj))*conv/dx
     .              *a*cos(GElatr(ii,jj))*
     .              uairDY(ii,jj,mz-kk+1) +
     .             (GElonh(ii,jj+1) - GElonh(ii,jj))*conv/dx
     .              *a*cos(GElatr(ii,jj))*
     .              vairDY(ii,jj,mz-kk+1)
     
        vv(kk)   = (GElatr(ii+1,jj) - GElatr(ii,jj))/dx*a
     .              *uairDY(ii,jj,mz-kk+1) +                   
     .             (GElatr(ii,jj+1) - GElatr(ii,jj))/dx*a
     .              *vairDY(ii,jj,mz-kk+1)       

        WD(ll,kk) = 0.0

        IF (uu(kk)> 0.0.and.vv(kk)>=0.0)  
     .  WD(ll,kk) = 3.0*pi/2.0-atan(vv(kk)/uu(kk))
        IF (uu(kk)< 0.0.and.vv(kk)>=0.0)  
     .  WD(ll,kk) = 5.0*Pi/2.0-atan(vv(kk)/uu(kk))
        IF (uu(kk)< 0.0.and.vv(kk)<=0.0)  
     .  WD(ll,kk) = 5.0*Pi/2.0-atan(vv(kk)/uu(kk))
        IF (uu(kk)> 0.0.and.vv(kk)<=0.0)  
     .  WD(ll,kk) = 3.0*Pi/2.0-atan(vv(kk)/uu(kk))
        IF (uu(kk)==0.0.and.vv(kk)>=0.0) 
     .  WD(ll,kk) = Pi
        IF (uu(kk)==0.0.and.vv(kk)<=0.0) 
     .  WD(ll,kk) = 0.0
        IF (WD(ll,kk)>2*Pi)
     .  WD(ll,kk) = WD(ll,kk) - 2*Pi

        WD(ll,kk) = WD(ll,kk) * 180.0/Pi

        if( WD(ll,kk) < 0.0 )       
     .  WD(ll,kk) = WD(ll,kk) + 360.0
     
        WD(ll,kk) = max(0.0,min(360.0,WD(ll,kk)))
     
        WS(ll,kk) = sqrt(uairDY(ii,jj,mz-kk+1)*uairDY(ii,jj,mz-kk+1)+
     .                   vairDY(ii,jj,mz-kk+1)*vairDY(ii,jj,mz-kk+1))
       
        ELSE
        WS(ll,kk) = 0.0 ; WD(ll,kk)   = 0.0
        END IF 
       
       END DO

      END DO

C +   3.3 Output 
C +   ----------

      DO ll = 1,ONEsta  ; ii = x_one(ll) ; jj = y_one(ll)       
       do kk=1,ONElev
        one7(ll,kk) = tairDY(ii,jj,mz-kk+1) - 273.15
        one8(ll,kk) = uairDY(ii,jj,mz-kk+1)
        one9(ll,kk) = vairDY(ii,jj,mz-kk+1)
       one10(ll,kk) =   qvDY(ii,jj,mz-kk+1) * 1000.
       end do
      END DO

C +   ************
      CALL UNwrite(ID__nc_one,'TT', ipr_nc_one, ONEsta, ONElev, 1,one7)
      CALL UNwrite(ID__nc_one,'UU', ipr_nc_one, ONEsta, ONElev, 1,one8)
      CALL UNwrite(ID__nc_one,'VV', ipr_nc_one, ONEsta, ONElev, 1,one9)
      CALL UNwrite(ID__nc_one,'WS', ipr_nc_one, ONEsta, ONElev, 1,WS)
      CALL UNwrite(ID__nc_one,'WD', ipr_nc_one, ONEsta, ONElev, 1,WD)
      CALL UNwrite(ID__nc_one,'QQ', ipr_nc_one, ONEsta, ONElev, 1,one10)
      CALL UNwrite(ID__nc_one,'RH', ipr_nc_one, ONEsta, ONElev, 1,RH) 
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta  ; ii = x_one(ll) ; jj = y_one(ll)       
       do kk=1,ONElev
        one7(ll,kk) = TUkvm(ii,jj,mz-kk+1)
        one8(ll,kk) = ect_TE(ii,jj,mz-kk+1)
        one9(ll,kk) = gplvDY(ii,jj,mz-kk+1)*grvinv
       one10(ll,kk) = wairDY(ii,jj,mz-kk+1)
       end do
      END DO

C +   ************
      CALL UNwrite(ID__nc_one,'KZ' ,ipr_nc_one, ONEsta, ONElev, 1, one7)
      CALL UNwrite(ID__nc_one,'TKE',ipr_nc_one, ONEsta, ONElev, 1, one8)
      CALL UNwrite(ID__nc_one,'ZZ' ,ipr_nc_one, ONEsta, ONElev, 1, one9)
      CALL UNwrite(ID__nc_one,'WW' ,ipr_nc_one, ONEsta, ONElev, 1,one10)
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta  ; ii = x_one(ll) ; jj = y_one(ll)       
       do kk=1,ONElev
        one7(ll,kk) = QIHY(ii,jj,mz-kk+1)   * 1000.
        one8(ll,kk) = QWHY(ii,jj,mz-kk+1)   * 1000.
        one9(ll,kk) = QSHY(ii,jj,mz-kk+1)   * 1000.
       one10(ll,kk) = QRHY(ii,jj,mz-kk+1)   * 1000.
       end do
      END DO

C +   ************
      CALL UNwrite(ID__nc_one,'QI', ipr_nc_one, ONEsta, ONElev, 1,one7)
      CALL UNwrite(ID__nc_one,'QW', ipr_nc_one, ONEsta, ONElev, 1,one8)
      CALL UNwrite(ID__nc_one,'QS', ipr_nc_one, ONEsta, ONElev, 1,one9)
      CALL UNwrite(ID__nc_one,'QR', ipr_nc_one, ONEsta, ONElev, 1,one10)
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll) 
       one1(ll)     = tairSL(ii,jj)         - 273.15
       one2(ll)     =  pstDY(ii,jj)         * 10. 
       one3(ll)     = rainHY(ii,jj)         * 1000.
       one4(ll)     = snowHY(ii,jj)         * 1000.
       one5(ll)     = rainCA(ii,jj)         * 1000.
      END DO
      
C +   ************      
      CALL UNwrite(ID__nc_one,'ST', ipr_nc_one, ONEsta, 1, 1, one1)     
      CALL UNwrite(ID__nc_one,'SP', ipr_nc_one, ONEsta, 1, 1, one2) 
      CALL UNwrite(ID__nc_one,'RR', ipr_nc_one, ONEsta, 1, 1, one3)
      CALL UNwrite(ID__nc_one,'SF', ipr_nc_one, ONEsta, 1, 1, one4)
      CALL UNwrite(ID__nc_one,'CP', ipr_nc_one, ONEsta, 1, 1, one5)
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll) 
       one1(ll)     = Radsol(ii,jj)
       one2(ll)     = Rad_IR(ii,jj)
       one3(ll)     = hsenSL(ii,jj)
       one4(ll)     = hlatSL(ii,jj)
       one5(ll)     = firmSL(ii,jj)
      END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'SWD',ipr_nc_one, ONEsta, 1, 1, one1)   
      CALL UNwrite(ID__nc_one,'LWD',ipr_nc_one, ONEsta, 1, 1, one2) 
      CALL UNwrite(ID__nc_one,'LWU',ipr_nc_one, ONEsta, 1, 1, one5)
      CALL UNwrite(ID__nc_one,'SHF',ipr_nc_one, ONEsta, 1, 1, one3)
      CALL UNwrite(ID__nc_one,'LHF',ipr_nc_one, ONEsta, 1, 1, one4)
C +   ************ 

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll) 
       one1(ll)     = albeSL(ii,jj)
       one2(ll)     = cld_SL(ii,jj)
       one3(ll)     = runoTV(ii,jj)
       one4(ll)     = evapTV(ii,jj)
      END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'AL' ,ipr_nc_one, ONEsta, 1, 1, one1)
      CALL UNwrite(ID__nc_one,'CC' ,ipr_nc_one, ONEsta, 1, 1, one2)
      CALL UNwrite(ID__nc_one,'RU' ,ipr_nc_one, ONEsta, 1, 1, one3)
      CALL UNwrite(ID__nc_one,'EV' ,ipr_nc_one, ONEsta, 1, 1, one4)
C +   ************ 

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll)
 
        one1(ll)  = 0.0 ; one2(ll)  = 0.0 ; one3(ll) = 0.0

        DO k = mzabso+1,mz
         one3(ll) = (    pstDY(ii,jj)  * sigma(k)+ptopDY)
     .            / (ra*tairDY(ii,jj,k)*(1.+.608   *qvDY(ii,jj,k)))
     .            * (   gpmiDY(ii,jj,k)-          gpmiDY(ii,jj,k+1)) 
         one1(ll) = one1(ll) + one3(ll) * qwHY(ii,jj,k)
         one2(ll) = one2(ll) + one3(ll) * qiHY(ii,jj,k)
 
        END DO
 
        one4(ll)  = 1.5 * ( one1(ll) / 20.d-6
     .                    + one2(ll) / 40.d-6 ) *grvinv

      END DO

C +   ************
      CALL UNwrite(ID__nc_one,'COD',ipr_nc_one, ONEsta, 1, 1, one4)
C +   ************ 

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll) 
       one1(ll)     = SLlmo(ii,jj)
       one2(ll)     = SLuus(ii,jj)
       one3(ll)     = SLuts(ii,jj)
       one4(ll)     = SLuqs(ii,jj)
      END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'LMO',ipr_nc_one, ONEsta, 1, 1, one1)   
      CALL UNwrite(ID__nc_one,'UUS',ipr_nc_one, ONEsta, 1, 1, one2)  
      CALL UNwrite(ID__nc_one,'UTS',ipr_nc_one, ONEsta, 1, 1, one3)
      CALL UNwrite(ID__nc_one,'UQS',ipr_nc_one, ONEsta, 1, 1, one4)
C +   ************ 

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      SL_wge=0.0 ; SLlwge=0.0 ; SLuwge=0.0

C +   ***********
      call WGustE
C +   *********** 

      DO ll = 1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll) 
       one1(ll)     = SL_wge(ii,jj)
       one2(ll)     = SLlwge(ii,jj)
       one3(ll)     = SLuwge(ii,jj)
      END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'WGE',ipr_nc_one, ONEsta, 1, 1, one1)   
      CALL UNwrite(ID__nc_one,'WGL',ipr_nc_one, ONEsta, 1, 1, one2)  
      CALL UNwrite(ID__nc_one,'WGU',ipr_nc_one, ONEsta, 1, 1, one3)
C +   ************ 

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll) 
       one1(ll)     = SL_z0 (ii,jj,1)
       one2(ll)     = SL_r0 (ii,jj,1)
       one3(ll)     = SWaSNo(ii,jj,1)
      END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'Z0' ,ipr_nc_one, ONEsta, 1, 1, one1) 
      CALL UNwrite(ID__nc_one,'R0' ,ipr_nc_one, ONEsta, 1, 1, one2)
      CALL UNwrite(ID__nc_one,'SWA',ipr_nc_one, ONEsta, 1, 1, one3)
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO k=1,nsx
      DO j=1,my
      DO i=1,mx
          zn1IB (i,j,k) = 1.
          zn2IB (i,j,k) = 0.
          zn3IB (i,j,k) = 0.
          wet_IB(i,j,k) = 0.
      END DO
      END DO
      END DO
 
      DO kk=nsno,1,-1
      DO k=1,nsx
      DO j=1,my
      DO i=1,mx
          zn3IB (i,j,k) = dzsSNo(i,j,k,kk)   + zn3IB(i,j,k)
          wet_IB(i,j,k) = rosSNo(i,j,k,kk)   * dzsSNo(i,j,k,kk)
     .  		* 1.d3/ ro_Wat       *(1.+0.*wasSNo(i,j,k,kk))
     .  		+ wet_IB(i,j,k)
          zn1IB (i,j,k) = zn1IB(i,j,k)
     .  		* max(zero,sign(unun,
     .  		  ro_ice-20.-rosSNo(i,j,k,kk)))
          zn2IB (i,j,k) = dzsSNo(i,j,k,kk)   * zn1IB(i,j,k)
     .  				     + zn2IB(i,j,k)
      END DO
      END DO
      END DO
      END DO
 
      DO k=1,nsx
      DO j=1,my
      DO i=1,mx
          wet_IB(i,j,k)  = wet_IB(i,j,k) +     SWaSNo(i,j,k)
          zn2IB (i,j,k)  = zn2IB (i,j,k) * (1.- zn1IB(i,j,k))
          zn1IB (i,j,k)  = zn3IB (i,j,k)      - zn0IB(i,j,k)
          mbIB  (i,j,k)  = wet_IB(i,j,k)      - mb0IB(i,j,k)
      END DO
      END DO
      END DO

      DO ll = 1,ONEsta ; ii = x_one(ll) ; jj = y_one(ll) 
       one1(ll)     = zn1IB(ii,jj,1)
       one2(ll)     = zn2IB(ii,jj,1)
       one3(ll)     = zn3IB(ii,jj,1)
       one4(ll)     = mbIB (ii,jj,1)
      END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'ZN' ,ipr_nc_one, ONEsta, 1, 1, one1) 
      CALL UNwrite(ID__nc_one,'ZN1',ipr_nc_one, ONEsta, 1, 1, one1) 
      CALL UNwrite(ID__nc_one,'ZN2',ipr_nc_one, ONEsta, 1, 1, one2) 
      CALL UNwrite(ID__nc_one,'ZN3',ipr_nc_one, ONEsta, 1, 1, one3) 
      CALL UNwrite(ID__nc_one,'MB' ,ipr_nc_one, ONEsta, 1, 1, one4)
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta  ; DO nn = 1,nsno
         ii = x_one(ll) ;    jj = y_one(ll)
       one11(ll,nn) = dzsSNo(ii,jj,1,nn) 
       one12(ll,nn) = rosSNo(ii,jj,1,nn)        
      END DO ; END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'DZsn',ipr_nc_one, ONEsta,nsno,1,one11)  
      CALL UNwrite(ID__nc_one,'ROsn',ipr_nc_one, ONEsta,nsno,1,one12)
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta  ; DO nn = 1,nsno
         ii = x_one(ll) ;    jj = y_one(ll)
       one11(ll,nn) = g1sSNo(ii,jj,1,nn) 
       one12(ll,nn) = g2sSNo(ii,jj,1,nn)        
      END DO ; END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'G1sn',ipr_nc_one, ONEsta,nsno,1,one11) 
      CALL UNwrite(ID__nc_one,'G2sn',ipr_nc_one, ONEsta,nsno,1,one12)
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta  ; DO nn = 1,nsno
         ii = x_one(ll) ;    jj = y_one(ll)
       one11(ll,nn) = tisSNo(ii,jj,1,nn) - TfSnow 
       one12(ll,nn) = wasSNo(ii,jj,1,nn)        
      END DO ; END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'TIsn',ipr_nc_one, ONEsta,nsno,1,one11) 
      CALL UNwrite(ID__nc_one,'WAsn',ipr_nc_one, ONEsta,nsno,1,one12)
C +   ************

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO ll = 1,ONEsta  ; DO nn = 1,nsol+1
         ii = x_one(ll) ;    jj = y_one(ll)
       one14(ll,nn) = TsolTV(ii,jj,1,nn) - TfSnow 
       one15(ll,nn) = Eta_TV(ii,jj,1,nn)        
      END DO ; END DO

C +   ************      
      CALL UNwrite(ID__nc_one,'SLT',ipr_nc_one, ONEsta,nsol+1,1,one14) 
      CALL UNwrite(ID__nc_one,'SLQ',ipr_nc_one, ONEsta,nsol+1,1,one15)
C +   ************


 1000 continue

C +   3.4 NetCDF File Closure
C +   -----------------------

      IF              (ID__nc_one.ne.-1)                         THEN

C +       ************
          call UNclose(ID__nc_one)
C +       ************

      END IF

      WKxy1= 0 ; WKxy2= 0 ; WKxy3= 0 ; WKxy4= 0 
       

      end subroutine

C +-------------------------------------------------------------------+
      subroutine WGustE
C +-------------------------------------------------------------------+
C |                                                                   |
C | MAR GUSTS                                         14-09-2001  MAR |
C |   SubRoutine WGustE computes diagnostic Wind Gust Estimates       |
C |                                                                   |
C +-------------------------------------------------------------------+
C |                                                                   |
C | This routine aims at estimating wind gusts. It also includes the  |
C | computation of a bounding interval around the estimate which aims |
C | to contain with a high probability observed gusts.                |
C |                                                                   |
C | Ref.  : Brasseur O., MWR, 129, 5-25.                              |
C | ^^^^^^^                                                           |
C |                                                                   |
C | Input : - uairDY : U-wind                                         |
C | ^^^^^^^ - vairDY : V-wind                                         |
C |         - wairDY : W-wind                                         |
C |         - tairDY : REAL temperature                               |
C |         - tairSL : surface temperature                            |
C |         - qvDY   : specific humidity                              |
C |         - qwHY   : cloud dropplets                                |
C |         - qiHY   : ice crystals                                   |
C |         - qrHY   : rain                                           |
C |         - qsHY   : snow                                           |
C |         - SLuts  : surface heat flux                              |
C |         - ect_TE : turbulent kinetic energy                       |
C |         - zzDY   : level heights                                  |
C |         - sh     : surface elevation                              |
C |         - sigma  : sigma levels                                   |
C |         - pstDY  : pressure depth                                 |
C |         - ptopDY : pressure at the top of the model               |
C |                                                                   |
C | Output: - SL_wge : gust estimate                          (m/s)   |
C | ^^^^^^^ - SLlwge : lower bound of the bounding interval   (m/s)   |
C |         - SLuwge : upper bound of the bounding interval   (m/s)   |
C |                                                                   |
C +-------------------------------------------------------------------+


       IMPLICIT NONE


C +---General Variables
C +   -----------------

       INCLUDE "MARdim.inc"
       INCLUDE "MARgrd.inc"
       INCLUDE "MAR_DY.inc"
       INCLUDE "MAR_HY.inc"
       INCLUDE "MAR_TE.inc"
       INCLUDE "MAR_SL.inc"


C +---Local variables
C +   ---------------

       INTEGER l,lev_up,lev_dw,top_bl,kzi

       REAL    int_buo,   local_tke,tke_min,aux0,aux1,  aux2
       REAL    coeff,     dtkemin,  ra,     cp,  gravit,cap
       REAL    ENERGY_low,ENERGY_est
       REAL    coeffmin,  srf_TKE,  wstar,  rzero

       REAL    tetae(mz),   tetav(mz),buoy(mz),normv(mz),filECT(mz)
       REAL    mean_tke(mz),tmpECT(mz)


C +---Data
C +   ----


       DATA dtkemin  /    1.e-5    /
       DATA coeffmin /    0.01     /
       DATA cap      /    0.28586  /
       DATA ra       /  287.       /
       DATA cp       / 1004.       /
       DATA gravit   /    9.81     /
       DATA rzero    /    0.0      /


C +    ====================================================

       DO j=1,my
       DO i=1,mx

C +    ====================================================

        DO k=1,mz

C +---Compute Virtual and Equivalent Potential Temperature
C +   ----------------------------------------------------

         tetav(k)=tairDY(i,j,k)
     .           *(100./(pstDY(i,j)*sigma(k)+ptopDY))**cap
     .           *(1.+0.608*qvDY(i,j,k)-qwHY(i,j,k)-qiHY(i,j,k)
     .                                 -qrHY(i,j,k)-qsHY(i,j,k))
C +     ^^^ Virtual potential temperature


C +---Compute wind norm
C +   -----------------

         normv(k)=SQRT(uairDY(i,j,k)*uairDY(i,j,k)
     .                +vairDY(i,j,k)*vairDY(i,j,k)
     .                +wairDY(i,j,k)*wairDY(i,j,k)/10000.)

        ENDDO


C +---Compute Integrated Buoyancy
C +   ---------------------------

        k = mz
        int_buo=        (0.5*tetav(k )+0.5*tetav(k -1)
     .                   -0.5*tetav(mz)-0.5*tetav(mz-1))*gravit
     .                  /(0.5*tetav(mz)+0.5*tetav(mz-1))/2.
     .                  *(gplvDY(i,j,k-1)-gplvDY(i,j,k))/gravit
        int_buo=MAX(0.,int_buo)
        buoy(k)=int_buo

        k = mz - 1
        int_buo=int_buo-(0.5*tetav(k )+0.5*tetav(k -1)
     .                   -0.5*tetav(mz)-0.5*tetav(mz-1))*gravit
     .                  /(0.5*tetav(mz)+0.5*tetav(mz-1))/2.
     .                  *(gplvDY(i,j,k-1)-gplvDY(i,j,k))/gravit
        int_buo=MAX(0.,int_buo)
        buoy(k)=int_buo

        DO k=mz,2,-1

         int_buo=int_buo-(0.5*tetav(k )+0.5*tetav(k -1)
     .                   -0.5*tetav(mz)-0.5*tetav(mz-1))*gravit
     .                  /(0.5*tetav(mz)+0.5*tetav(mz-1))/2.
     .                  *(gplvDY(i,j,k-1)-gplvDY(i,j,k))/gravit
         buoy(k)=int_buo

        ENDDO


C +---Filtering of turbulent kinetic energy
C +   -------------------------------------

        tmpECT( 1)=ect_TE(i,j, 1)
        tmpECT(mz)=ect_TE(i,j,mz)
        filECT( 1)=ect_TE(i,j, 1)
        filECT(mz)=ect_TE(i,j,mz)

        DO k=2,mz-1
         IF (i.ne.1.and.i.ne.mx.and.j.ne.1.and.j.ne.my) THEN
          tmpECT(k)=(4.*ect_TE(i,j,k)
     .              +2.*ect_TE(i-1,j  ,k)+2.*ect_TE(i+1,j  ,k)
     .              +2.*ect_TE(i  ,j-1,k)+2.*ect_TE(i  ,j+1,k)
     .              +1.*ect_TE(i-1,j-1,k)+1.*ect_TE(i-1,j+1,k)
     .              +1.*ect_TE(i+1,j-1,k)+1.*ect_TE(i+1,j+1,k))/16.
         ENDIF
        ENDDO

        DO k=2,mz-1
         aux1=(gplvDY(i,j,k-1)-gplvDY(i,j,k)) / gravit
         aux2=(gplvDY(i,j,k)-gplvDY(i,j,k+1)) / gravit
         filECT(k)=0.25*(2.                   *tmpECT(k  )
     .                  +2.*(aux2/(aux1+aux2))*tmpECT(k-1)
     .                  +2.*(aux1/(aux1+aux2))*tmpECT(k+1))
        ENDDO


C +---Determination of mean tke below level k
C +   - - - - - - - - - - - - - - - - - - - -

        aux1=0.
        aux2=0.

        DO k=mz-1,2,-1
         aux1=aux1+(gplvDY(i,j,k)+gplvDY(i,j,k+1))/2.
     .            *filECT(k)
     .            *(gplvDY(i,j,k)-gplvDY(i,j,k+1))
         aux2=aux2+(gplvDY(i,j,k)+gplvDY(i,j,k+1))/2.
     .            *(gplvDY(i,j,k)-gplvDY(i,j,k+1))
        ENDDO

        mean_tke(k)=aux1/aux2


C +---Compute wstar
C +   -------------

c #WW   wstar=MAX(zero,(-gravit*zi__TE(i,j)/290.*SLuts(i,j))**(1./3.))


C +---Evaluation of Gust Wind Speed
C +   -----------------------------

C +---Initial value
C +   - - - - - - -

        SL_wge(i,j)=MAX(SL_wge(i,j),normv(mz))
        SLuwge(i,j)=MAX(SLuwge(i,j),normv(mz))
        SLlwge(i,j)=MAX(SLlwge(i,j),normv(mz))

        lev_dw=1
        lev_up=mz
        top_bl=mz


        DO k=mz-1,2,-1

C +---Determination of the Top of Boundary Layer
C +   - - - - - - - - - - - - - - - - - - - - -

         coeff  = coeffmin 
     .          + 0.1*((gplvDY(i,j,k)-gplvDY(i,j,mz+1))
     .                                   /gravit-2000.)/1000.

         tke_min=coeff*(filECT(mz)+filECT(mz-1))*0.5

         IF (top_bl.eq.mz           .and.
     .       filECT(k+1).gt.tke_min .and.
     .       filECT(k  ).le.tke_min     ) top_bl=k


C +---Upper bound on Gust Wind Speed
C +   - - - - - - - - - - - - - - - -

         local_tke=(filECT(k)+filECT(k-1))*0.5

         IF (local_tke.gt.tke_min.and.top_bl.eq.mz) THEN

          IF (SLuwge(i,j).lt.normv(k)) THEN
           SLuwge(i,j)=MAX(SLuwge(i,j),normv(k))
C +        ^^^ Max Wind Speed  (m/s)
           lev_up  =k
C +        ^^^ Level of Max Wind Speed  (m/s)
          ENDIF

         ENDIF


C +---Lower bound on Wind Gust
C +   - - - - - - - - - - - -

         ENERGY_low = 2.5/11.*filECT(k)    ! Source
     .              + buoy(k)              ! Sink
C +
         IF (ENERGY_low.ge.0.) THEN
          SLlwge(i,j)=MAX(normv(k),SLlwge(i,j))
C +       ^^^ Min Wind Speed  (m/s)
          lev_dw  =k
         ENDIF


C +---Estimate of Wind Gust
C +   - - - - - - - - - - -

         ENERGY_est = MAX(mean_tke(k),2.5/11.*filECT(k))  ! Source
     .              + buoy(k)                             ! Sink

         IF (ENERGY_est.ge.0.) THEN
          SL_wge(i,j)=MAX(0.5*(normv(k)+normv(k-1)),SL_wge(i,j))
         ENDIF


        ENDDO   ! {Loop on k}

C +    ====================================================
C +
       ENDDO    ! {Loop on i}
       ENDDO    ! {Loop on j}
C +
C +    ====================================================

      return
      end
