! 
! 
! 
!
!   PPPP     OOOO     M   M     M   M     EEEEE    W E S T
!  P    P   O    O   M M M M   M M M M   E              
!  PPPPP    O    O   M  M  M   M  M  M   EEEE         A F R I C A N
!  P        O    O   M     M   M     M   E                 
!  P         OOOO    M     M   M     M    EEEEE          M O N S O O N
!
!  -                                -               -        -     -
!  P O S T - P R O C E S S I N G    O F    T H E    M A R    M O D E L
!   -                                -               -        -     -
!
!  Version 1.2b
!  12.09.2005 (01.03.2006)
!
!
!
!
! This program post-processes the output files of the model MAR.
!
! This program concatenates a given number of monthly files (DYN*.nc,
! MAR*.nc, ANI*.nc, and RAD*.nc) produced by MAR into one single file.
! Each variable (see inside the code) is extracted and transformed into
! a daily value.
! The maximum number of each kind of monthly files is 12 (therefore
! maximum one year can be concatened). All files must correspond to the
! same year.
! The output file is named 'WAXXX.YYYY.M1-M2.nc' where:
!  - XXX  is the run name
!  - YYYY is the year to post-processe
!  - M1   is the month at which to start the post-processing
!  - M2   is the month at which to end the post-processing
!                                            (as specified in POMME.ctr)
! The post-processing of XXX starts on the 01/M1/YYYY and ends the last
! day of M2/YYYY.
!


      program POMME


!
! OUTLINE
! _______
!
! 0.   Declaration of variables and initialisation of parameters
! I.   Open/create files
!        I.1.  Open input files
!        I.2.  Create output file
! II.  Read/write dimensions and time-independent variables 
!        II.1. Dimensions ----------------------+
!        II.2. Time-independent variables       | n-month loop to
!                - sigma                        | initialise the output
!                - longitude and latitude       | file
!                - topography and surface kind -+
! III. Read/write variables (daily values)
!        III.1. Temporal dimension -------------+
!        III.2. 6-hourly, 3-D data              |
!                - GreenL                       |
!                - CAPE                         |
!                - pstar                        |
!                - pressure* (4-D)              |
!        III.3. 6-hourly, 4-D data              |
!                - zonal wind                   |
!                - meridional wind              |
!                - geopotential                 | 1-day loop to read/
!                - air temperature              | write the variables
!                - specific humidity            | in daily mean
!                - moist static energy*         |
!        III.4. 30-min, 3-D data                |
!                - albedo                       |
!                - surface temperature          |
!                - latent heat                  |
!                - sensible heat                |
!                - evapotranspiration           |
!                - rainfall                     |
!        III.5. 30-min, 4-D data                |
!                - soil moisture ---------------+
! IV.  Close/save files
!        IV.1.  Close input files
!        IV.2.  Save output file
! V.   Annexes: subroutines and functions
!        - Subroutine handle_leap_years
!        - Subroutine handle_err
!        - Subroutine title
!        - Function num_2_string
!
! * computed variable


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                      +
! 0. DECLARATION OF VARIABLES AND INITIALISATION OF PARAMETERS         +
!                                                                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!use the Fortran 90 netCDF interface
!/usr/local/pub/include/NETCDF.mod copied in the local directory


      use NETCDF


!no implicit variables (all variables must be declared)


      implicit none


!include netCDF declaration (netCDF functions & co)


      include '/usr/local/pub/include/netcdf-3.4/netcdf.inc'
               !f90 -c POMME.f 
               !f90 -o POMME POMME.o -L/usr/local/pub/lib64 -lnetcdf


!include variables specific to this program


      include 'POMME.inc'


!data


      data month /'01','02','03','04','05','06','07','08','09','10',
     &            '11','12'/
      data month_name /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',
     &                 'Sep','Oct','Nov','Dec'/


!title


      t1 = 0
      t2 = 0
      t3 = 0

      write(*,*)
      write(*,2) ' ---------------------------------------------------',
     &           '------------------- '
      write(*,2) '|                                                   ',
     &           '                   |'
      write(*,2) '|  PPPP     OOOO     M   M     M   M     EEEEE    ',
     &           'W E S T              |'
      write(*,2) '| P    P   O    O   M M M M   M M M M   E',
     &           '                              |'
      write(*,2) '| PPPPP    O    O   M  M  M   M  M  M   EEEE      ',
     &           '   A F R I C A N     |'
      write(*,2) '| P        O    O   M     M   M     M   E',
     &           '                              |'
      write(*,2) '| P         OOOO    M     M   M     M    EEEEE    ',
     &           '      M O N S O O N  |'
      write(*,2) '|                                                   ',
     &           '                   |'
      write(*,2) '| -                                -              ',
     &           ' -        -     -    |'
      write(*,2) '| P O S T - P R O C E S S I N G    O F    T H E   ',
     &           ' M A R    M O D E L  |'
      write(*,2) '| -                                -              ',
     &           ' -        -     -    |'
      write(*,2) '|                                                   ',
     &           '                   |'
      write(*,2) '| March 2006 - Version 1.2b',
     &           '                                            |'
      write(*,2) '|                                                   ',
     &           '                   |'
      write(*,2) ' ---------------------------------------------------',
     &           '------------------- '
      write(*,*)

 2    format(a,a)


!user's parameters                                               

      open(7, file='POMME.ctr', status='old', access='sequential',
     &          form='formatted', err=7999, iostat=nerr)

      read(7,*)
      read(7,*)
      read(7,*)
      read(7,'(a3)') run
      read(7,'(a4)') year
      read(7,'(a2)') cmonth_slot(1)
      backspace 7
      read(7,'(i2)') month_slot(1)
      read(7,'(a2)') cmonth_slot(2)
      backspace 7
      read(7,'(i2)') month_slot(2)

      close(7)

!display informations


      write(*,1) ' >> ','POST-PROCESSING OF THE RUN ', run,
     &           ' FROM ', cmonth_slot(1), '/', year,
     &           ' TO ',   cmonth_slot(2), '/', year, ' <<'
 1    format(a,5x,a,a3,a,a2,a1,a4,a,a2,a1,a4,5x,a)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                      +
! I. OPEN/CREATE FILES                                                 +
!                                                                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      call title(1, 'OPEN/CREATE FILES', t1, t2, t3)


!______________________________________________________________________
!                                                                      |
! I.1. OPEN INPUT FILES                                                |
!______________________________________________________________________|


!define time


      call handle_leap_years(year, day(1:12), cday, daycum(1:12))

      day(0)    = 0  !so that the first month is ok when writing the
      daycum(0) = 0  !temporal dimension


!define file kind


      c3tmp(1) = 'MAR'
      c3tmp(2) = 'DYN'
      c3tmp(3) = 'RAD'
      c3tmp(4) = 'ANI'


!define input file names, open files, get their netCDF ID


      write(*,'(/a)') 'open the input files (in input/) : '
      do n = 1, 4  !loop on file kinds
        write(*,'(/$)') 
        do m = month_slot(1), month_slot(2)  !loop on months
          ctmp50 = 'input/'//c3tmp(n)//'.'//run//'.'//year//'.'
     &             //month(m)//'.01-'//cday(m)
     &             //'.nc'
          status = nf90_open(ctmp50, nf90_nowrite, input_id(m,n))
                       if (status .ne. nf_noerr) call handle_err(status)
          write(*,'(1x,a,$)') trim(ctmp50(7:50))
        end do
      end do


!______________________________________________________________________
!                                                                      |
! I.2. CREATE OUTPUT FILE                                              |
!______________________________________________________________________|


!define and create the output file, get its netCDF ID


      write(*,'(//a/)') 'create the output file (in output/):'

      ctmp50 = 'WA'//run//'.'
     &         //year//'.'//cmonth_slot(1)//'-'//cmonth_slot(2)
     &         //'.nc'
      write(*,'(1x,a)') ctmp50  !e.g.: WAF86.1987.01-03.nc  

      !netCDF: define mode

      status = nf90_create(ctmp50, nf_clobber, output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !netCDF: add title

      ctmp50 = 'MAR exp: F86 - '
     &         //cmonth_slot(1)//'-'//cmonth_slot(2)//'/'//year
      status = nf90_put_att(output_id, NF90_GLOBAL,
     &                      'title', trim(ctmp50))
                       if (status .ne. nf_noerr) call handle_err(status)

      !netCDF: quit define mode and enter data mode

      status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                      +
! II. READ/WRITE DIMENSIONS AND TIME-INDEPENDENT VARIABLES             +
!                                                                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      call title(1,'READ/WRITE DIMENSIONS & TIME-INDEPENDENT VARIABLES',
     &           t1, t2, t3)


! This part defines the common setup of the netCF output file.
! The spatial dimensions are defined as well as the temporal dimension.
! Notice that the latter is handled a first time here as the unlimited
! dimension is required for the grid definition of the netCDF file. It
! will be handled again later on for any other time step.
! Then to the output file are added some common variables, time
! independent: sigma, longitude, latitude, topography, surface type.


! Note: temporal dimension in netCDF files
! -----
! There are 2 methods to write the temporal dimension of a netCDF file.
!
! 1) The first one is to write all values of the time axis in the netCDF
!    file in one go. Then each variable is written time step(s) by time
!    step(s). It works as long as the variable position along the time
!    axis is correctly given (with a start(:) array).
!    e.g.: status = nf90_put_var( file_id, time_dimension_var_id,
!                                 time_values(1:time_dimension_length) )
!          start = (/1, 1, time_index, 0 /)          !for a xyt variable
!          status = nf90_put_var( file_id, var_id, values(1:x,1:y),
!                                 start(1:3) )
!
! 2) The second method is to write the time axis step(s) by step(s) in
!    the same time as the variable values are written step(s) by step(s)
!    also. It imports again to provide the nf90_put_var function the 
!    position of these values along the time axis with the help of a
!    start(:) array.
!    e.g.: start = (/ time_index, 0, 0, 0 /)
!          count = (/ 1, 0, 0, 0 /)
!          status = nf90_put_var( file_id, time_dimension_var_id,
!                                 time_value(1), start(1), count(1) )
!          start = (/1, 1, time_index, 0 /)          !for a xyt variable
!          status = nf90_put_var( file_id, var_id, values(1:x,1:y),
!                                 start(1:3) )


      write(*,*)

      m = month_slot(1)

      first_month = .true.

      write(*,'(a3,1x,a4,a1/a)')  month_name(m),year,':',  '---------'


!______________________________________________________________________
!                                                                      |
! II.1. DIMENSIONS                                                     |
!______________________________________________________________________|


! Are extracted here the following dimensions: time, x, y, level (sigma
! levels), and level2 (ground levels).
! The output file gets these dimensions:
! 0 time:  unlimited, daily step (1 unlimited dimension is compulsory)
! 1 x:     1 -> 124 (a 5-grid point buffer zone is removed all around)
! 2 y:     1 -> 124 (a 5-grid point buffer zone is removed all around)
! 3 z_pbl: 1 -> 10 (= sigma levels from 30 to 40)
! 4 z:     1 -> 40 (= all sigma levels)
! 5 depth: 1 -> 7


!     ..................................................................
!     read dimensions


      !get dimension IDs
                                        !1: MAR netCDF file
      status = nf90_inq_dimid(input_id(m,1), 'time', dim_ID(0))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_inq_dimid(input_id(m,1), 'x', dim_ID(1))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_inq_dimid(input_id(m,1), 'y', dim_ID(2))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_inq_dimid(input_id(m,1), 'level', dim_ID(3))
                       if (status .ne. nf_noerr) call handle_err(status)
                                        !4: ANI netCDF file
      status = nf90_inq_dimid(input_id(m,4), 'level2', dim_ID(4))
                       if (status .ne. nf_noerr) call handle_err(status)

      !get dimension names and lengths, get dimension attributes

      att_units = ' '
      name      = ' '
      i = 1
      do n = 0, 4
        if (n .eq. 4) i = 4  !take the right input file for 'level2'
        status = nf90_inquire_dimension(input_id(m,i), dim_ID(n),
     &                                  name(n), dim_len(n))
                       if (status .ne. nf_noerr) call handle_err(status)
        status = nf90_get_att(input_id(m,i), dim_ID(n), 'units',
     &                        att_units(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      end do

      !get dimension values

      values1d = 0.

      !dim_len(0) = 1                !t: only first value required
      dim_len(1) = dim_len(1) - 10  !x (buffer zone: 2 * 5 grid meshes)
      dim_len(2) = dim_len(2) - 10  !y (buffer zone: 2 * 5 grid meshes)

      start = (/ 6, 6, 1, 1 /)                                  !xyzs
      count = (/ dim_len(1),dim_len(2),dim_len(3),dim_len(4) /) 

      !if (first_month) then
      !  n = 0
      !  status = nf90_get_var(input_id(m,1),dim_ID(n), values1d
      !&                      (n,1:dim_len(n)), start(3:3), start(3:3))
      !                      !values1d(n,1),1,1) doesn't work
      !                if (status .ne. nf_noerr) call handle_err(status)
      !end if

      i = 1
      do n = 1, 4
        if (n .eq. 4) i = 4
        status = nf90_get_var(input_id(m,i), dim_ID(n), values1d
     &                        (n,1:count(n)), start(n:n), count(n:n))
                       if (status .ne. nf_noerr) call handle_err(status)
      end do

      !time: compute time dimension with a daily step

      !if (first_month) then
      !  values1d(0,2:) = 0.
      !  do i = 2, daycum(m+2)-daycum(m-1)
      !    values1d(0,i) = values1d(0,i-1) + real(24)
      !  end do
      !end if


!     ..................................................................
!     define dimensions: (0) time, (1) x, (2) y, (3) z_pbl, (4) z,
!                        (5) depth


      !dimension dimensions

      dim_len(5) = dim_len(4)     !depth: 1 -> 7
      dim_len(4) = dim_len(3)     !z:     1 -> 40 (= all sigma levels)
      dim_len(3) = dim_len(3)-30  !z_pbl: 1 -> 10 (= sigma levels 30:40)
      dim_len(2) = dim_len(2)     !y:     1 -> 124
      dim_len(1) = dim_len(1)     !x:     1 -> 124
      !if (first_month) then       !t:     1 -> day per season
      !  dim_len(0) = daycum(m+2)-daycum(m-1)
      !end if

      name(0) = 'time'
      name(1) = 'x'
      name(2) = 'y'
      name(3) = 'z_pbl'
      name(4) = 'z'
      name(5) = 'depth'
      att_units(5) = att_units(4)  !dimensions z and depth shifted (+1)
      att_units(4) = att_units(3)  !so shift units attribute accordingly

      !vertical levels: shift z (3 > 4), depth (4 > 5), create z_pbl (3)

      do k = 1, dim_len(5)
        values1d(5,k) = values1d(4,k)     !depth now is the dimension #5
      end do
      values1d(4,:) = 0.
      do k = 1, dim_len(4)
        values1d(4,k) = values1d(3,k)     !z     now is the dimension #4
      end do
      values1d(3,:) = 0.
      do k = 1, 10
        values1d(3,k) = values1d(4,k+30)  !z_pbl now is the dimension #3
      end do

      !netCDF: enter define mode

      status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define dimensions

      status = nf90_def_dim(output_id, name(0), nf90_unlimited, 
     &                      dim_ID(0))
                       if (status .ne. nf_noerr) call handle_err(status)
      do n = 1, 5
        status = nf90_def_dim(output_id, name(n), dim_len(n),
     &                        dim_ID(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      end do

      !define associated variables + dimension attributes (units, range)

      do n = 0, 5
        status = nf90_def_var(output_id, name(n), NF90_FLOAT, 
     &                        dim_ID(n), i5tmp(n))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, i5tmp(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                              !att_units(n)(1:len(att_units(n)))
                       if (status .ne. nf_noerr) call handle_err(status)
        if (n .ne. 0) then 
          rmin = values1d(n,1)
          rmax = values1d(n,1)
          do i = 1, dim_len(n)
            rmin = min(values1d(n,i),rmin)
            rmax = max(values1d(n,i),rmax)
          end do
          att_range(n,1) = rmin
          att_range(n,2) = rmax
          status = nf90_put_att(output_id, i5tmp(n), 'actual_range', 
     &                          att_range(n,1:2))
                       if (status .ne. nf_noerr) call handle_err(status)
        end if
      end do

      !netCDF: quit define mode and enter data mode

      status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)


!     ..................................................................
!     write dimensions


      !!time
      !
      !if (first_month) then
      !  n = 0
      !  start = (/ 1, 0, 0, 0 /)
      !  count = (/ dim_len(n), 0, 0, 0 /)
      !  status = nf90_put_var(output_id, dim_ID(n),
      !&                        values1d(n,1:dim_len(0)))
      !                if (status .ne. nf_noerr) call handle_err(status)
      !  write(*,9211) name(n), values1d(n,1), values1d(n,dim_len(n)),
      !&                trim(att_units(n)), 1,dim_len(n)
      !end if

      !x, y, z_pbl, z, depth

      do n = 1, 5
        status = nf90_put_var(output_id, dim_ID(n), 
     &                        values1d(n,1:dim_len(n)))
                       if (status .ne. nf_noerr) call handle_err(status)
        if (n .eq. 1) write(*,9211) name(n), values1d(n,1),
     &                values1d(n,dim_len(n)), att_units(n), 1,dim_len(n)
        if (n .eq. 2) write(*,9211) name(n), values1d(n,1),
     &                values1d(n,dim_len(n)), att_units(n), 1,dim_len(n)
        if (n .ge. 3) write(*,9211) name(n), values1d(n,1),
     &                values1d(n,dim_len(n)), att_units(n), 1,dim_len(n)
      end do

      write(*,*)


!     ..................................................................
!     formats


 9211 format(1x,a8,f10.3,' /',f9.3,1x,a33,1x,'(',i1,':',i3,')')


!______________________________________________________________________
!                                                                      |
! II.2. TIME-INDEPENDENT VARIABLES                                     |
!______________________________________________________________________|


! We handle here the following common and time independent variables:
!   0: sigma (z)
!   1: lon   (xy)
!   2: lat   (xy)
!   3: sh    (xy)
!   4: isol  (xy)


      name(0) = 'sigma'
      name(1) = 'lon'
      name(2) = 'lat'
      name(3) = 'sh'
      name(4) = 'isol'

!     ...........
      do n = 0, 4  !loop on variables
!     ...........

      write(*,'(1x,a8,$)') name(n)


!     ..................................................................
!     read variables


      !get variable IDs

      if (n .ne. 0) then
        status = nf90_inq_varid(input_id(m,1), name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      end if

      !get variable attributes (units, long_name, actual_range) + name

      if (n. ne. 0) then
        att_units(n)    = ' '
        att_longname(n) = ' '
        att_range(n,1)  = 0.
        att_range(n,2)  = 0.
        status = nf90_get_att(input_id(m,1), var_id(n), 'units',
     &                        att_units(n))
                       if (status .ne. nf_noerr) call handle_err(status)
        status = nf90_get_att(input_id(m,1), var_id(n), 'long_name',
     &                        att_longname(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      else
        att_longname(0) = 'Normalised pressure'
        att_units(0)    = att_units(4)   !from dimension definition
        att_range(0,1)  = att_range(4,1) !and be cautious as 4 will be
        att_range(0,2)  = att_range(4,2) !re-attributed to 'isol'
      end if

      !get variable values

      if (n .ne. 0) then
        values2d = 0.
        start = (/ 6, 6, 0, 0 /)  !xy
        count = (/ dim_len(1), dim_len(2), 0, 0 /)

        status = nf90_get_var(input_id(m,1), var_id(n),
     &                        values2d(1:count(1),1:count(2)),
     &                        start(1:2), count(1:2))
                       if (status .ne. nf_noerr) call handle_err(status)
        rmin = values2d(1,1)
        rmax = values2d(1,1)
        do j = 1, dim_len(2)
        do i = 1, dim_len(1)
          rmin = min(values2d(i,j),rmin)
          rmax = max(values2d(i,j),rmax)
        end do
        end do
        att_range(n,1) = rmin
        att_range(n,2) = rmax
      else
        do k = 1, dim_len(4)        !sigma created from the
          sigma(k) = values1d(4,k)  !dimension(4) values (i.e. z)
        end do
      end if


!     ..................................................................
!     define variables


      !netCDF: enter define mode

      status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

      count = (/ dim_ID(1), dim_ID(2), dim_ID(4), 0 /)  !xy, z

      if (n .ne. 0) then
        status = nf90_def_var(output_id, name(n), NF90_FLOAT,
     &                        count(1:2), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      else
        status = nf90_def_var(output_id, name(n), NF90_FLOAT,
     &                        count(3), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      end if

      !define variable attributes

      ctmp50 = att_units(n)
      status = nf90_put_att(output_id, var_id(n), 'units', 
     &                      ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
      ctmp50 = att_longname(n)
      status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                      ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_put_att(output_id, var_id(n), 'actual_range', 
     &                      att_range(n,1:2))
                       if (status .ne. nf_noerr) call handle_err(status)
 
      !netCDF: quit define mode and enter data mode

      status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)


!     ..................................................................
!     write variables


      if (n .ne. 0) then
        status = nf90_put_var(output_id, var_id(n), 
     &                        values2d(1:dim_len(1),1:dim_len(2)))
                       if (status .ne. nf_noerr) call handle_err(status)
      else
        status = nf90_put_var(output_id, var_id(n),
     &                        sigma(1:dim_len(4)))
                       if (status .ne. nf_noerr) call handle_err(status)
      end if

      if (n. eq. 0) then
        write(*,9221) att_range(n,1:2), att_units(n), dim_len(4)
      else
        write(*,9222) att_range(n,1:2), att_units(n), dim_len(1),
     &               dim_len(2)
      end if

!     ...........
      end do     !loop on variables
!     ...........

      write(*,*)


!     ..................................................................
!     formats


      !z, z_pbl, depth
 9221 format(f10.3,' /',f10.3,1x,a28,'(   ,   ,',i2,')')
      !xy
 9222 format(f10.3,' /',f10.3,1x,a28,'(',i3,',',i3,',  )')
      !xyz
 9223 format(f10.3,' /',f10.3,1x,a28,'(',i3,',',i3,',',i2,')')


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                      +
! III. READ/WRITE VARIABLES (DAILY VALUES)                             +
!                                                                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      call title(1, 'READ/WRITE VARIABLES (DAILY VALUES)', t1, t2, t3)


! This part handles all variables needed in the output file: dynamics,
! transport, surface variables. A loop on months is performed to add
! each day consecutively.
! For each day the time axis is computed on basis of the first time step
! of the first handled month.
! The variables are read by day (4 or 48 time steps, depending on
! whether the variables are 6-hourly or 30-min.), transformed into daily
! values, and written into the output file.
! The temporal dimension is therefore handled by day.
!
! Important notes:
! - The MAR and DYN netCDF input files contains one excedentary time
!   step (1st day of the next month, 0h). Anyway, we decided not to take
!   the first time step of each month as it represents the 6 past hour
!   (last day of the previous month) .
! - The ANI netCDF intput file contains a bug for the evapotranspiration
!   and rainfall at the 1st time step (1st time step = cumulated value
!   since simulation began no reset to 0 in SBCnew.f of its cumulated
!   value). Therefore this value is systematically reset to 0.


      time_index = 0  !positional index for the writing of time 
                      !dependent values such as time axis and variables

!     -----------------------------------
      DO m = month_slot(1), month_slot(2)  !loop on months
!     -----------------------------------

      first_month = .false.
      if (m.eq.month_slot(1)) first_month = .true.

!     ----------------
      DO d = 1, day(m)  !loop on days
!     ----------------        

      time_index = time_index + 1
      
      if (first_month .and. d.ne.1) first_month = .false.

      write(*,'(/a1,i2,a1,1x,i2,1x,a3,1x,a4,a1,$)') '[',time_index,']',
     &                                     d, month_name(m), year, ':'
      !write(*,'(/i2,1x,a3,1x,a4)') d, month_name(m), year


!______________________________________________________________________
!                                                                      |
! III.1. Temporal dimension                                            |
!______________________________________________________________________|


!     ..................................................................
!     read dimension


      if (first_month) then

        !get dimension value

        timeaxis(1) = 0.

        start = (/ 1, 0, 0, 0 /)
        count = (/ 1, 0, 0, 0 /)
        status = nf90_get_var(input_id(m,1), dim_ID(0), timeaxis(1:1),
     &                        start(1:1), count(1:1))
       !status = nf90_get_var(input_id(m,1),dim_ID(0),values1d(0,1),1,1)
                       if (status .ne. nf_noerr) call handle_err(status)

        !get dimension units

        status = nf90_get_att(input_id(m,1), dim_ID(0), 'units',
     &                        att_units(0))
                       if (status .ne. nf_noerr) call handle_err(status)

      else

        !get dimension value: compute the daily step

        timeaxis(1) = timeaxis(1) + real(24)

        !get dimension units

        status = nf90_get_att(input_id(m,1), dim_ID(0), 'units',
     &                        att_units(0))
                       if (status .ne. nf_noerr) call handle_err(status)
      end if


!     ..................................................................
!     write dimension


      start = (/ time_index, 0, 0, 0 /)
      count = (/ 1, 0, 0, 0 /)
      status = nf90_put_var(output_id, dim_ID(0), timeaxis(1:1),
     &                      start(1:1), count(1:1))
                       if (status .ne. nf_noerr) call handle_err(status)
      write(*,'(f9.1,1x,a/a)') timeaxis(1), trim(att_units(0)),
     &                         '-----------------'


!______________________________________________________________________
!                                                                      |
! III.2. 6-hourly, 3-D data                                            |
!______________________________________________________________________|


! Time steps of retrieval: t
! ------------------------
! Data are retrieved from t to t+3 for each day, the daily mean is
! computed on basis of these 4 values.
! All months except January with start on the 1st, 18h:
!   day (d) 1           2           3           4           5         
!   hour    00 06 12 18 00 06 12 18 00 06 12 18 00 06 12 18 00 06 12 18
!   num      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
!               t-------->  t-------->  t-------->  t-------->  t------
!  => t = d * 4 - 2
! In January, start on the 1st, 18h:
!   day (d) 1           2           3           4           5     
!   hour    00 06 12 18 00 06 12 18 00 06 12 18 00 06 12 18 00 06 12 18
!   num               1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
!                           t-------->  t-------->  t---------> t------
!  => t = d * 4 - 5 = (d * 4 - 2) - 3

      t = d * 4 - 2

! Read/written variables: GreenL, CAPE, pstar (MAR netCDF files)
! -----------------------
! pstar must be the last one to be read so that it can be kept in memory
! and be used with sigma for the computation of a variable 'pressure' as
! explained below.
!
! Note: pstar is 6-hourly in MAR netCDF file; its units will be changed
!       from kPa to hPa just after its retrieval.

! From the normalised pressure (sigma), we compute the pressure levels
! in hPa thanks to:
!
!   sigma = (pressure - pt) / p*   where  p: pressure (kPa)
!                                        pt: pressure at the top of the
!                                            atmosphere (0.01 kPa)
!                                        p*: ps - pt: pressure thickness
!                                        ps: pressure at the surface
!
!   pressure = sigma * p* + pt = sigma * (ps - pt) + pt
!
! In the model MAR, pstar = ps - pt. Therefore, the pressure is:
!
!     pressure(x,y,z,t) = [ sigma(z) * pstar(x,y,t) + 0.01 ] * 10.
!
! (* 10.) to get hPa instead of kPa.


      name            = ''
      att_units       = ''
      att_longname    = ''
      
      name(0)         = 'GreenL'
      name(1)         = 'CAPE'
      name(2)         = 'pstar'
      name(3)         = 'pressure'


!
! 0-2 GreenL, CAPE, pstar
!________________________________________________


!     ...........
      do n = 0, 2  !loop on variables
!     ...........

      write(*,'(1x,a8,$)') name(n)


!     ..................................................................
!     read variables


      !get variable IDs

      status = nf90_inq_varid(input_id(m,1), name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)
                       
      !get variable attributes (units, long_name)
      
      status = nf90_get_att(input_id(m,1), var_id(n), 'units',
     &                      att_units(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_get_att(input_id(m,1), var_id(n), 'long_name',
     &                      att_longname(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      
      att_units(2)    = 'hPa'
      att_longname(2) = 'Pressure thickness'  !orthograph correction
      
      !get variable values

      values4d = 0.
      start = (/ 6, 6, t, 0 /)
      count = (/ dim_len(1), dim_len(2), 4, 0 /)

      status = nf90_get_var(input_id(m,1), var_id(n), values4d
     &                      (1:count(1),1:count(2),1,1:count(3)),
     &                      start(1:3), count(1:3))
                       if (status .ne. nf_noerr) call handle_err(status)

      !compute daily values

      daily4d = 0.
      daily4d(:,:,1) = ( values4d(:,:,1,1)+values4d(:,:,1,2)
     &                  +values4d(:,:,1,3)+values4d(:,:,1,4) ) / 4.

      !compute pstar in hPa

      if (n .eq. 2) daily4d(:,:,1) = daily4d(:,:,1) * 10.

      !compute min/max values

      rmin = daily4d(1,1,1)
      rmax = daily4d(1,1,1)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        rmin = min(daily4d(i,j,1),rmin)
        rmax = max(daily4d(i,j,1),rmax)
      end do
      end do
      att_range(n,1) = rmin
      att_range(n,2) = rmax


!     ..................................................................
!     define variables


      if (first_month) then

      !netCDF: enter define mode

        status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

        count = (/ dim_ID(1), dim_ID(2), dim_ID(0), 0 /)
        status = nf90_def_var(output_id, name(n), NF90_FLOAT,
     &                        count(1:3), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variable attributes (units, long name, range)

        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, var_id(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_longname(n)
        status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)

      !netCDF: quit define mode and enter data mode

        status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      else

      !get variable IDs

        status = nf90_inq_varid(output_id, name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      end if  !.not. first_month


!     ..................................................................
!     write variables


      start = (/ 1, 1, time_index, 0 /)
      status = nf90_put_var(output_id, var_id(n), daily4d(
     &                      1:dim_len(1),1:dim_len(2),1), start(1:3))
                       if (status .ne. nf_noerr) call handle_err(status)

      write(*,9222) att_range(n,1:2),att_units(n), dim_len(1),dim_len(2)


!     ...........
      end do     !loop on variables
!     ...........


!
! 3 pressure
!________________________________________________


      n = 3
      
      write(*,'(1x,a8,$)') name(n)


!     ..................................................................
!     read variables


      !get variable attributes (units, long_name)
      
      att_units(n)    = 'hPa'
      att_longname(n) = 'Pressure'

      !compute pressure
      
      daily4d2 = 0.  !pressure
      do k=1, dim_len(4)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        daily4d2(i,j,k) = sigma(k) * daily4d(i,j,1) + 0.1  !in hPa
      end do
      end do
      end do

      !compute min/max values

      rmin = daily4d2(1,1,1)
      rmax = daily4d2(1,1,1)
      do k=1, dim_len(4)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        rmin = min(daily4d(i,j,k),rmin)
        rmax = max(daily4d(i,j,k),rmax)
      end do
      end do
      end do
      att_range(n,1) = rmin
      att_range(n,2) = rmax


!     ..................................................................
!     define variables


      if (first_month) then

      !netCDF: enter define mode

        status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

        count = (/ dim_ID(1), dim_ID(2), dim_ID(4), dim_ID(0) /)
        status = nf90_def_var(output_id, name(n), NF90_FLOAT, count,
     &                        var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variable attributes (units, long name, range)

        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, var_id(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_longname(n)
        status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)

      !netCDF: quit define mode and enter data mode

        status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      else

      !get variable IDs

        status = nf90_inq_varid(output_id, name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      end if  !.not. first_month


!     ..................................................................
!     write variables


      start = (/ 1, 1, 1, time_index /)
      status = nf90_put_var(output_id, var_id(n), daily4d2(
     &                      1:dim_len(1),1:dim_len(2),1:dim_len(4)),
     &                      start)
                       if (status .ne. nf_noerr) call handle_err(status)
      write(*,9223) att_range(n,1:2), att_units(n),
     &              dim_len(1), dim_len(2), dim_len(4)


!______________________________________________________________________
!                                                                      |
! III.3. 6-hourly, 4-D data                                            |
!______________________________________________________________________|


! Time steps of retrieval: t
! ------------------------
! Data are retrieved from t to t+3 for each day, the daily mean is
! computed on basis of these 4 values.
! All months except January with start on the 1st, 18h:
!   day (d) 1           2           3           4           5         
!   hour    00 06 12 18 00 06 12 18 00 06 12 18 00 06 12 18 00 06 12 18
!   num      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
!               t-------->  t-------->  t-------->  t-------->  t------
!  => t = d * 4 - 2
! In January, start on the 1st, 18h:
!   day (d) 1           2           3           4           5     
!   hour    00 06 12 18 00 06 12 18 00 06 12 18 00 06 12 18 00 06 12 18
!   num               1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
!                           t-------->  t-------->  t---------> t------
!  => t = d * 4 - 5 = (d * 4 - 2) - 3

      t = d * 4 - 2

! Read/written variables: zonal wind (uairDY), meridional wind (vairDY),
! ----------------------- air temperature (tairDY), specific humidity
!                         (qvDY)(DYN netCDF files),
!                         geopotential (zzDY) (MAR netCDF files)
!
! Note: winds are extracted for the whole atmosphere (useful to spot the
!       TEJ and AEJ) instead of only the boundary layer.
!
! The MSE is established from the geopotential, the air temperature, and
! the specific humidity thanks to:
!
!   MSE =  g Z + cp T + Lv q
!
! where g  is the gravitational acceleration                (9.81 m/s^2)
!       Lv is the latent heat of vaporisation         (2.5008 10^6 J/kg)
!       cp is the specific heat at constant pressure for air
!                                                   (1004.708845 J/kg K)
!       Z  is the geopotential      (i.e.   zzDY)
!       T  is the air temperature   (i.e. tairDY)
!       q  is the specific humidity (i.e.   qvDY)


      name            = ''
      att_units       = ''
      att_longname    = ''
      
      name(0)         = 'uairDY'  !all the vertical
      name(1)         = 'vairDY'  !all the vertical
      name(2)         = 'zzDY'    !all the vertical
      name(3)         = 'tairDY'  !all the vertical
      name(4)         = 'qvDY'    !all the vertical
      name(5)         = 'MSE'     !all the vertical


!
! 0-1 zonal and meridional winds
!________________________________________________


!     ...........
      do n = 0, 1  !loop on variables
!     ...........

      write(*,'(1x,a8,$)') name(n)


!     ..................................................................
!     read variables


      !get variable IDs

      status = nf90_inq_varid(input_id(m,2), name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !get variable attributes (units, long_name)

      status = nf90_get_att(input_id(m,2), var_id(n), 'units',
     &                      att_units(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_get_att(input_id(m,2), var_id(n), 'long_name',
     &                      att_longname(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !get variable values

      values4d = 0.
      start = (/ 6, 6, 1, t /)  !xyzt  !start = (/ 6, 6, 30, t /)
      count = (/ dim_len(1), dim_len(2), dim_len(4), 4 /)

      status = nf90_get_var(input_id(m,2), var_id(n), values4d(1:
     &                      count(1),1:count(2),1:count(3),1:count(4)),
     &                      start, count)
                       if (status .ne. nf_noerr) call handle_err(status)

      !compute daily mean

      daily4d = 0.
      daily4d(:,:,:) = ( values4d(:,:,:,1)+values4d(:,:,:,2)
     &                  +values4d(:,:,:,3)+values4d(:,:,:,4) ) / 4.

      !compute min/max values

      rmin = daily4d(1,1,1)
      rmax = daily4d(1,1,1)
      do k=1, dim_len(4)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        rmin = min(daily4d(i,j,k),rmin)
        rmax = max(daily4d(i,j,k),rmax)
      end do
      end do
      end do
      att_range(n,1) = rmin
      att_range(n,2) = rmax


!     ..................................................................
!     define variables


      if (first_month) then

      !netCDF: enter define mode

        status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

        count = (/ dim_ID(1), dim_ID(2), dim_ID(4), dim_ID(0) /)   !xyzt

        status = nf90_def_var(output_id, name(n), NF90_FLOAT, count,
     &                        var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variable attributes

        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, var_id(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_longname(n)
        status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
 
      !netCDF: quit define mode and enter data mode

        status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      else

      !get variable IDs

        status = nf90_inq_varid(output_id, name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      end if  !.not. first_month


!     ..................................................................
!     write variables


      start = (/ 1, 1, 1, time_index /)
      status = nf90_put_var(output_id, var_id(n), daily4d
     &                      (1:dim_len(1),1:dim_len(2),1:dim_len(4)),
     &                      start)
                       if (status .ne. nf_noerr) call handle_err(status)

      write(*,9223) att_range(n,1:2), att_units(n),
     &              dim_len(1), dim_len(2), dim_len(4)

!     ...........
      end do     !loop on variables
!     ...........


!
! 2-4 zzDY, tairDY, qvDY
!________________________________________________


!     ...........
      do n = 2, 4  !loop on variables
!     ...........

      write(*,'(1x,a8,$)') name(n)


!     ..................................................................
!     read variables


      p = 2  !input file type
      if (n .eq. 2) p = 1
      
      !get variable IDs

      status = nf90_inq_varid(input_id(m,p), name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !get variable attributes (units, long_name)

      status = nf90_get_att(input_id(m,p), var_id(n), 'units',
     &                      att_units(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_get_att(input_id(m,p), var_id(n), 'long_name',
     &                      att_longname(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !get variable values

      values4d = 0.
      start = (/ 6, 6, 1, t /)    !xyzt
      count = (/ dim_len(1), dim_len(2), dim_len(4), 4 /)

      status = nf90_get_var(input_id(m,p), var_id(n), values4d(1:
     &                      count(1),1:count(2),1:count(3),1:count(4)),
     &                      start, count)
                       if (status .ne. nf_noerr) call handle_err(status)

      !compute daily values

      daily4d = 0.
      daily4d(:,:,:) = ( values4d(:,:,:,1)+values4d(:,:,:,2)
     &                  +values4d(:,:,:,3)+values4d(:,:,:,4) ) / 4.

      !compute min/max values

      rmin = daily4d(1,1,1)
      rmax = daily4d(1,1,1)
      do k=1, dim_len(4)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        rmin = min(daily4d(i,j,k),rmin)
        rmax = max(daily4d(i,j,k),rmax)
      end do
      end do
      end do
      att_range(n,1) = rmin
      att_range(n,2) = rmax
      
      !compute MSE
      
      if      (n .eq. 2) then  !geopotential
        daily4d2 = 0.
        do k=1, dim_len(4)
        do j=1, dim_len(2)
        do i=1, dim_len(1)
          daily4d2(i,j,k) = 9.81 * daily4d(i,j,k)
        end do
        end do
        end do
      else if (n .eq. 3) then  !temperature
        do k=1, dim_len(4)
        do j=1, dim_len(2)
        do i=1, dim_len(1)
          daily4d2(i,j,k) = daily4d2(i,j,k) + 1004.708845*daily4d(i,j,k)
        end do
        end do
        end do
      
      else if (n .eq. 4) then  !specific humidity
        do k=1, dim_len(4)
        do j=1, dim_len(2)
        do i=1, dim_len(1)
          daily4d2(i,j,k) = daily4d2(i,j,k) + 2500800 * daily4d(i,j,k)
        end do
        end do
        end do
      
      end if


!     ..................................................................
!     define variables


      if (first_month) then

      !netCDF: enter define mode

        status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

        count = (/ dim_ID(1), dim_ID(2), dim_ID(4), dim_ID(0) /)  !xyzt

        status = nf90_def_var(output_id, name(n), NF90_FLOAT, count,
     &                        var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variable attributes

        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, var_id(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_longname(n)
        status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
 
      !netCDF: quit define mode and enter data mode

        status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      else

      !get variable IDs

        status = nf90_inq_varid(output_id, name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      end if  !.not. first_month


!     ..................................................................
!     write variables


      start = (/ 1, 1, 1, time_index /)
      status = nf90_put_var(output_id, var_id(n), daily4d
     &                      (1:dim_len(1),1:dim_len(2),1:dim_len(4)),
     &                      start)
                       if (status .ne. nf_noerr) call handle_err(status)

      write(*,9223) att_range(n,1:2), att_units(n),
     &              dim_len(1), dim_len(2), dim_len(4)

!     ...........
      end do     !loop on variables
!     ...........


!
! 5 moist static energy
!________________________________________________


      n = 5

      write(*,'(1x,a8,$)') name(n)


!     ..................................................................
!     compute variable


      !get variable attributes (units, long_name)

      att_longname(n) = 'Moist Static Energy'
      att_units(n)    = 'J/kg'

      !compute min/max values

      rmin = daily4d2(1,1,1)
      rmax = daily4d2(1,1,1)
      do k=1, dim_len(4)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        rmin = min(daily4d2(i,j,k),rmin)
        rmax = max(daily4d2(i,j,k),rmax)
      end do
      end do
      end do
      att_range(n,1) = rmin
      att_range(n,2) = rmax


!     ..................................................................
!     define variables


      if (first_month) then

      !netCDF: enter define mode

        status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

        count = (/ dim_ID(1), dim_ID(2), dim_ID(4), dim_ID(0) /)  !xyzt

        status = nf90_def_var(output_id, name(n), NF90_FLOAT, count,
     &                        var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variable attributes

        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, var_id(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_longname(n)
        status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
 
      !netCDF: quit define mode and enter data mode

        status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      else

      !get variable IDs

        status = nf90_inq_varid(output_id, name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      end if  !.not. first_month


!     ..................................................................
!     write variables


      start = (/ 1, 1, 1, time_index /)
      status = nf90_put_var(output_id, var_id(n), daily4d2
     &                      (1:dim_len(1),1:dim_len(2),1:dim_len(4)),
     &                      start)
                       if (status .ne. nf_noerr) call handle_err(status)

      write(*,9223) att_range(n,1:2), att_units(n),
     &              dim_len(1), dim_len(2), dim_len(4)


!______________________________________________________________________
!                                                                      |
! III.4. 30-min, 3-D data                                              |
!______________________________________________________________________|


! Time steps of retrieval: t
! ------------------------
! Data are retrieved from t to t+47 for each day, the daily mean is
! computed on basis of these 48 values.
! All months except January with start on the 1st, 18h:
!   day (d) 1                         2                          3
!   hour    00    01        23        00        01     23        00
!   min     00 30 60 90...1380 1410 1440 1470 1500...2820 2850 2880 2910
!   num      1  2  3  4.....47   48   49   50   51.....95   96   97   98
!               t---------------------->   t---------------------->   t-
!  => t = (d - 1) * 48 + 2 = d * 48 - 46
! In January, start on the 1st, 18h:
!   day (d) 1                              2                          3
!   hour    00.....18.....19.....23        00        01     23        00
!   min     00...1080...1140...1380 1410 1440 1470 1500...2820 2850 2880
!   num             1      3     11   12   13   14   15.....59   60   61
!                                               t---------------------->
!  => t = (d - 1) * 48 + 14 = d * 48 - 34 = (d * 48 - 46) + 12
!
! Note:
! -----
! The last time step of the 30-min data files is the last day 18h, which
! means the 48th value for the daily mean is missing. As that one is the
! first value of the 30-min data file for the next month and as that
! value is wrong (no reset in SBCnew.f, see previously), we assume it
! to be 0 in the computation of the daily value. [tt]

      t = d * 48 - 46
      tt = 48
      if (d.eq.day(m)) tt = 47

! Read/written variables: albedo (albeSL), surface temperature (tairSL),
! ----------------------- latent heat (hlatSL), sensible heat (hsenSL)
!                         (RAD netCDF files), rainfall (Raintt),
!                         evapotranspiration (Evapot) (ANI netCDF files)
! Notes:
! - tairSL is renamed tsrfSL (mistake in MAR, this is the surface
!   temperature and not the air temperature at the ground).
! - due to the mistake of the 1st time step of the rainfall and the
!   evapotranspiration in the ANI netCDF files, the 1st time step is
!   systematically reset to 0.


      name            = ''
      att_units       = ''
      att_longname    = ''
      
      name(0)         = 'albeSL'
      name(1)         = 'tairSL'
      name(2)         = 'hlatSL'
      name(3)         = 'hsenSL'
      name(4)         = 'Raintt'
      name(5)         = 'Evapot'


!
! 0-3 albeSL, tairSL, hlatSL, hsenSL
!________________________________________________


!     ...........
      do n = 0, 3  !loop on variables
!     ...........

      name(1) = 'tsrfSL'
      write(*,'(1x,a8,$)') name(n)
      name(1) = 'tairSL'


!     ..................................................................
!     read variables


      !get variable IDs

      status = nf90_inq_varid(input_id(m,3), name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)
                       
      !get variable attributes (units, long_name)

      status = nf90_get_att(input_id(m,3), var_id(n), 'units',
     &                      att_units(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      status = nf90_get_att(input_id(m,3), var_id(n), 'long_name',
     &                      att_longname(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      att_units(0) = '-'

      !get variable values

      values4d = 0.
      start = (/ 6, 6, t, 0 /)
      count = (/ dim_len(1), dim_len(2), tt, 0 /)

      status = nf90_get_var(input_id(m,3), var_id(n), values4d
     &                      (1:count(1),1:count(2),1,1:count(3)),
     &                      start(1:3), count(1:3))
                       if (status .ne. nf_noerr) call handle_err(status)

      !compute daily values

      daily4d = 0.
      do i = 1, tt
        daily4d(:,:,1) = daily4d(:,:,1) + values4d(:,:,1,i)
      end do
      daily4d(:,:,1) = daily4d(:,:,1) / 48.

      !compute min/max values

      rmin = daily4d(1,1,1)
      rmax = daily4d(1,1,1)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        rmin = min(daily4d(i,j,1),rmin)
        rmax = max(daily4d(i,j,1),rmax)
      end do
      end do
      att_range(n,1) = rmin
      att_range(n,2) = rmax


!     ..................................................................
!     define variables

      name(1) = 'tsrfSL'

      if (first_month) then

      !netCDF: enter define mode

        status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

        count = (/ dim_ID(1), dim_ID(2), dim_ID(0), 0 /)
        status = nf90_def_var(output_id, name(n), NF90_FLOAT,
     &                        count(1:3), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variable attributes (units, long name, range)

        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, var_id(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_longname(n)
        status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)

      !netCDF: quit define mode and enter data mode

        status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      else

      !get variable IDs

        status = nf90_inq_varid(output_id, name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      end if  !.not. first_month


!     ..................................................................
!     write variables


      start = (/ 1, 1, time_index, 0 /)
      status = nf90_put_var(output_id, var_id(n), daily4d(
     &                      1:dim_len(1),1:dim_len(2),1), start(1:3))
                       if (status .ne. nf_noerr) call handle_err(status)
      write(*,9222) att_range(n,1:2),att_units(n), dim_len(1),dim_len(2)


!     ...........
      end do     !loop on variables
!     ...........


!
! 4-5 Raintt, Evapot
!________________________________________________


! See explanations of the part III about processing of the 1st time step
! of each month for both variables.


!     ...........
      do n = 4, 5  !loop on variables
!     ...........

      write(*,'(1x,a8,$)') name(n)


!     ..................................................................
!     read variables


      !get variable IDs

      status = nf90_inq_varid(input_id(m,4), name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)
                       
      !get variable attributes (units, long_name)
      
      status = nf90_get_att(input_id(m,4), var_id(n), 'units',
     &                      att_units(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_get_att(input_id(m,4), var_id(n), 'long_name',
     &                      att_longname(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !get variable values

      values4d = 0.
      start = (/ 6, 6, t, 0 /)
      count = (/ dim_len(1), dim_len(2), tt, 0 /)

      status = nf90_get_var(input_id(m,4), var_id(n), values4d
     &                      (1:count(1),1:count(2),1,1:count(3)),
     &                      start(1:3), count(1:3))
                       if (status .ne. nf_noerr) call handle_err(status)

      !compute daily values
      
      daily4d = 0.
      do i = 1, tt  !sum of all values
        daily4d(:,:,1) = daily4d(:,:,1) + values4d(:,:,1,i)
      end do

      !compute min/max values

      rmin = daily4d(1,1,1)
      rmax = daily4d(1,1,1)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        rmin = min(daily4d(i,j,1),rmin)
        rmax = max(daily4d(i,j,1),rmax)
      end do
      end do
      att_range(n,1) = rmin
      att_range(n,2) = rmax


!     ..................................................................
!     define variables


      if (first_month) then

      !netCDF: enter define mode

        status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

        count = (/ dim_ID(1), dim_ID(2), dim_ID(0), 0 /)
        status = nf90_def_var(output_id, name(n), NF90_FLOAT,
     &                        count(1:3), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variable attributes (units, long name, range)

        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, var_id(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_longname(n)
        status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)

      !netCDF: quit define mode and enter data mode

        status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      else

      !get variable IDs

        status = nf90_inq_varid(output_id, name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      end if  !.not. first_month


!     ..................................................................
!     write variables


      start = (/ 1, 1, time_index, 0 /)
      status = nf90_put_var(output_id, var_id(n), daily4d(
     &                      1:dim_len(1),1:dim_len(2),1), start(1:3))
                       if (status .ne. nf_noerr) call handle_err(status)
      write(*,9222) att_range(n,1:2),att_units(n), dim_len(1),dim_len(2)

!     ...........
      end do     !loop on variables
!     ...........


!______________________________________________________________________
!                                                                      |
! III.5. 30-min, 4-D data                                              |
!______________________________________________________________________|


! Time steps of retrieval: t
! ------------------------
! Data are retrieved from t to t+47 for each day, the daily mean is
! computed on basis of these 48 values.
! All months except January with start on the 1st, 18h:
!   day (d) 1                         2                          3
!   hour    00    01        23        00        01     23        00
!   min     00 30 60 90...1380 1410 1440 1470 1500...2820 2850 2880 2910
!   num      1  2  3  4.....47   48   49   50   51.....95   96   97   98
!               t---------------------->   t---------------------->   t-
!  => t = (d - 1) * 48 + 2 = d * 48 - 46
! In January, start on the 1st, 18h:
!   day (d) 1                              2                          3
!   hour    00.....18.....19.....23        00        01     23        00
!   min     00...1080...1140...1380 1410 1440 1470 1500...2820 2850 2880
!   num             1      3     11   12   13   14   15.....59   60   61
!                                               t---------------------->
!  => t = (d - 1) * 48 + 14 = d * 48 - 34 = (d * 48 - 46) + 12
!
! Note:
! -----
! The last time step of the 30-min data files is the last day 18h, which
! means the 48th value for the daily mean is missing. As that one is the
! first value of the 30-min data file for the next month and as that
! value is wrong (no reset in SBCnew.f, see previously), we assume it
! to be 0 in the computation of the daily value. [tt]

      t = d * 48 - 46
      tt = 48
      if (d.eq.day(m)) tt = 47

! Read/written variable: soil moisture (eta_TV) (ANI netCDF files)
! ----------------------


      name            = ''
      att_units       = ''
      att_longname    = ''
      
      name(0)         = 'eta_TV'


!
! 0 soil moisture
!________________________________________________


      n = 0

      write(*,'(1x,a8,$)') name(n)


!     ..................................................................
!     read variables


      !get variable IDs

      status = nf90_inq_varid(input_id(m,4), name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !get variable attributes (units, long_name)

      status = nf90_get_att(input_id(m,4), var_id(n), 'units',
     &                      att_units(n))
                       if (status .ne. nf_noerr) call handle_err(status)
      status = nf90_get_att(input_id(m,4), var_id(n), 'long_name',
     &                      att_longname(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !get variable values

      values4d = 0.
      start = (/ 6, 6, 1, t /)    !xydt
      count = (/ dim_len(1), dim_len(2), dim_len(5), tt /)

      status = nf90_get_var(input_id(m,4), var_id(n), values4d(1:
     &                      count(1),1:count(2),1:count(3),1:count(4)),
     &                      start, count)
                       if (status .ne. nf_noerr) call handle_err(status)

      !compute daily values

      daily4d = 0.
      do i = 1, tt
        daily4d(:,:,:) = daily4d(:,:,:) + values4d(:,:,:,i)
      end do
      daily4d(:,:,:) = daily4d(:,:,:) / 48.
      !call daily_values(values4d2, count, 'ave', 30, daily4d)

      !compute min/max values

      rmin = daily4d(1,1,1)
      rmax = daily4d(1,1,1)
      do k=1, dim_len(5)
      do j=1, dim_len(2)
      do i=1, dim_len(1)
        rmin = min(daily4d(i,j,k),rmin)
        rmax = max(daily4d(i,j,k),rmax)
      end do
      end do
      end do
      att_range(n,1) = rmin
      att_range(n,2) = rmax

     
!     ..................................................................
!     define variables


      if (first_month) then

      !netCDF: enter define mode

        status = nf90_redef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variables

        count = (/ dim_ID(1), dim_ID(2), dim_ID(5), dim_ID(0) /)  !xyzt

        status = nf90_def_var(output_id, name(n), NF90_FLOAT, count,
     &                        var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      !define variable attributes

        ctmp50 = att_units(n)
        status = nf90_put_att(output_id, var_id(n), 'units', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
        ctmp50 = att_longname(n)
        status = nf90_put_att(output_id, var_id(n), 'long_name', 
     &                        ctmp50(1:len(ctmp50)))
                       if (status .ne. nf_noerr) call handle_err(status)
 
      !netCDF: quit define mode and enter data mode

        status = nf90_enddef(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)

      else

      !get variable IDs

        status = nf90_inq_varid(output_id, name(n), var_id(n))
                       if (status .ne. nf_noerr) call handle_err(status)

      end if  !.not. first_month


!     ..................................................................
!     write variables


      start = (/ 1, 1, 1, time_index /)
      status = nf90_put_var(output_id, var_id(n), daily4d
     &                      (1:dim_len(1),1:dim_len(2),1:dim_len(5)),
     &                      start)
                       if (status .ne. nf_noerr) call handle_err(status)

      write(*,9223) att_range(n,1:2), att_units(n),
     &             dim_len(1), dim_len(2), dim_len(5)


!     ------------
 997  END DO      !loop on days
!     ------------

!     -----------------------------------
 998  END DO                             !loop on months
!     -----------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                      +
! IV. CLOSE/SAVE FILES                                                 +
!                                                                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


 301  call title(1, 'CLOSE/SAVE FILES', t1, t2, t3)


!______________________________________________________________________
!                                                                      |
! IV.1. CLOSE INPUT FILES                                              |
!______________________________________________________________________|


      write(*,'(/a)') 'close input files'

      do n = 1, 4
      do m = month_slot(1), month_slot(2)
        status = nf90_close(input_id(m,n))
                       if (status .ne. nf_noerr) call handle_err(status)
      end do
      end do


!______________________________________________________________________
!                                                                      |
! IV.2. SAVE OUTPUT FILES                                              |
!______________________________________________________________________|


      write(*,'(a)') 'save the output file'

      status = nf90_close(output_id)
                       if (status .ne. nf_noerr) call handle_err(status)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                      +
! END                                                                  +
!                                                                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      go to 999


!error messages


 7999 write(*,*) 'POMME.ctr: ', nerr


!end of the programm


 999  stop

      end program POMME




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                      +
! V. ANNEXES: SUBROUTINES AND FUNCTIONS                                +
!                                                                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! This part contains the code of subroutines and functions used in the
! program POMME.
! These are:
!   - subroutine handle_leap_years (intyear, intday, charday, intdaycum)
!   - subroutine handle_err (status)
!   - subroutine title (tit_level, tit_text, tit1_num,tit2_num,tit3_num)
!   - function num_2_string (intnum)


!     __________________________________________________________________

      subroutine handle_leap_years (intyear, intday, charday, intdaycum)

!     This subroutine checks whether the year 'intyear' is a normal or a
!     leap year. It accordingly adapts the number of days per month
!     (integer array + character array) and the cumulated number of
!     days.
!
!     December 2005, E. Vanvyve

      !declaration

      implicit none

      integer(4)                      :: intyear, c
      integer(4),       dimension(12) :: intday, intdaycum
      character(len=2), dimension(12) :: charday

      intday  = (/  31,  28,  31,  30,  31,  30,  31,  31,  30,
     &              31,  30,  31 /)
      charday = (/ '31','28','31','30','31','30','31','31','30',
     &             '31','30','31' /)

      if (mod(intyear,4) .eq. 0) then
        if (mod(intyear,100) .eq. 0) then
          if (mod(intyear,400) .eq. 0) then
            intday(2)  = 29
            charday(2) = '29'
          end if
        end if
      end if

      intdaycum(1) = intday(1)
      do c = 2, 12
        intdaycum(c) = intdaycum(c-1) + intday(c)
      end do
      
      return

      end subroutine handle_leap_years

!     __________________________________________________________________

      subroutine handle_err(status)

!     This subroutine handles netCDF error codes (status) returned as
!     an integer by netCDF functions. It has to be called if nf_noerr is
!     not equal to 0 (the contrary meaning the netCDF function was
!     successful): if (nf_noerr .ne. 0) call handle_err(status)
!
!     This subroutine should be implemented with the netCDF function
!     'nf_strerror(status)', but this latter one refuses to work.
!
!     It uses the function 'num_2_string'.
!
!     September 2005, E. Vanvyve

      !declaration

      implicit none

      integer(4), intent(in) :: status
      character(len=40) :: message
      character(len=40) :: message_err
      character(len=6)  :: num_2_string


      !convert status into a text message in case of error

      select case (status)
        case (-33)   ; message = 'netCDF: badid (-33)'
        case (-35)   ; message = 'netCDF: exist (-35)'
        case (-36)   ; message = 'netCDF: inval (-36)'
        case (-37)   ; message = 'netCDF: perm (-37)'
        case (-38)   ; message = 'netCDF: notindefine (-38)'
        case (-39)   ; message = 'netCDF: indefine (-39)'
        case (-40)   ; message = 'netCDF: invalcoords (-40)'
        case (-41)   ; message = 'netCDF: maxdims (-41)'
        case (-42)   ; message = 'netCDF: nameinuse (-42)'
        case (-43)   ; message = 'netCDF: notatt (-43)'
        case (-44)   ; message = 'netCDF: maxatts (-44)'
        case (-45)   ; message = 'netCDF: badtype (-45)'
        case (-46)   ; message = 'netCDF: baddim (-46)'
        case (-47)   ; message = 'netCDF: unlimpos (-47)'
        case (-48)   ; message = 'netCDF: maxvars (-48)'
        case (-49)   ; message = 'netCDF: notvar (-49)'
        case (-50)   ; message = 'netCDF: global (-50)'
        case (-51)   ; message = 'netCDF: notnc (-51)'
        case (-52)   ; message = 'netCDF: sts (-52)'
        case (-53)   ; message = 'netCDF: maxname (-53)'
        case (-54)   ; message = 'netCDF: unlimit (-54)'
        case (-55)   ; message = 'netCDF: norecvars (-55)'
        case (-56)   ; message = 'netCDF: char (-56)'
        case (-57)   ; message = 'netCDF: edge (-57)'
        case (-58)   ; message = 'netCDF: stride (-58)'
        case (-59)   ; message = 'netCDF: badname (-59)'
        case (-60)   ; message = 'netCDF: range (-60)'
        case default ; message = 'netCDF: unknown error ('//
     &                            trim(num_2_string(status))//')'
      end select

      message_err = '?%@#! '//message(1:len(message))

      write(*,*) message_err

      stop

      return

      end subroutine handle_err

!     __________________________________________________________________

      subroutine title(tit_level, tit_text,
     &                 tit1_num, tit2_num, tit3_num)

!     This subroutine prints the string 'tit_text' as a title of level
!     'tit_level'.
!     Variables tit*_num must be initialised to 0 before the first call
!     to this subroutine. These are intent(inout) variables in order to
!     keep the title numbering.
!
!     October 2005, EV

      !declaration

      implicit none

      character(len=*), intent(in) :: tit_text
      integer(4), intent(in)       :: tit_level
      integer(4), intent(inout)    :: tit1_num, tit2_num, tit3_num
      character(len=1), dimension(10)  :: numstr
      character(len=4), dimension(10)  :: romstr
      character(len=1), dimension(10)  :: letstr
      data numstr /'1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ','10'/
      data romstr /'I   ','II  ','III ','IV  ','V   ','VI  ','VII ',
     &             'VIII','IX  ','X   '/
      data letstr /'a','b','c','d','e','f','g','h','i','j'/

      !title level and numbering
      
      if (tit_level .eq. 1) then
        tit1_num = tit1_num + 1
        tit2_num = 1
        tit3_num = 1
      end if
 
      !display titles

      if      (tit_level .eq. 1) then  !III. TEXTE
          write(*,*)
          write(*,'(2a)') '-------------------------------------------',
     &                    '-----------------------------'
          write(*,'(3a)') trim(romstr(tit1_num)), '. ', tit_text
          write(*,'(2a)') '-------------------------------------------',
     &                    '-----------------------------'

      else if (tit_level .eq. 2) then  !III.1. Texte
          write(*,*)
          write(*,'(5a)') trim(romstr(tit1_num)), '.', numstr(tit2_num),
     &                    '. ', tit_text
          write(*,'(1a)') '----------------------------------------'
          write(*,*)

      else if (tit_level .eq. 3) then  !III.1.a. Texte
          write(*,*)
          write(*,'(3a)') letstr(tit3_num), '. ', tit_text
          write(*,*)
          !write(*,*) repeat('.',len(letstr(tit3_num))+2+len(tit_text))

      end if

      !title level and numbering
      
      if      (tit_level .eq. 2) then
        tit1_num = tit1_num
        tit2_num = tit2_num + 1
        tit3_num = 1
      else if (tit_level .eq. 3) then
        tit3_num = tit3_num + 1
      end if

      return

      end subroutine title

!     __________________________________________________________________

      character(len=6) function num_2_string (intnum)

!     This function returns an integer between -99999 and 99999 as a
!     string of 6-character length (6 because of the minus sign for
!     negative numbers), left-adjusted.
!     E.g.: num_2_string(-79) returns '-79   '.
!
!     The function name must be declared in the main program:
!     character(len=6) :: num_2_string
!
!     September 2005, J. Langridge &  E. Vanvyve

      !declaration

      implicit none

      character(len=1), dimension(0:9) :: numstr
      character(len=1) :: minus
      integer(4) :: intnum
      integer(4) :: tensthous, thous, hunds, tens, sings, tmp
      data numstr/'0','1','2','3','4','5','6','7','8','9'/

      !check 'intnum' relevance

      if (intnum .lt. -99999 .or. intnum .gt. 99999) then
        write(*,*) 'intnum: ', intnum
        stop '?%@#! num_2_string'
      end if

      !handle negative numbers

      if (intnum .lt. 0) then
        minus = '-'
        intnum = intnum * (-1)
      end if

      !get each number separately

      tensthous = intnum / 10000
      tmp       = intnum - tensthous * 10000
      thous     = tmp / 1000
      tmp       = tmp - thous * 1000
      hunds     = tmp / 100
      tmp       = tmp - hunds * 100
      tens      = tmp / 10
      sings     = tmp - tens * 10

      !convert into a string

      if (intnum .ge. 0)     then
        num_2_string = numstr(sings)
      end if
      if (intnum .ge. 10)    then
        num_2_string = numstr(tens)//num_2_string(1:1)
      end if
      if (intnum .ge. 100)   then
        num_2_string = numstr(hunds)//num_2_string(1:2)
      end if
      if (intnum .ge. 1000)  then
        num_2_string = numstr(thous)//num_2_string(1:3)
      end if
      if (intnum .ge. 10000)  then
        num_2_string = numstr(tensthous)//num_2_string(1:4)
      end if

      if (minus .eq. '-') then  !negative numbers
        num_2_string = minus//num_2_string(1:len(num_2_string))
      end if

      return

      end function num_2_string

