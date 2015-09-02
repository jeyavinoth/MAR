! 
! 
! 
!
!     PPPP     OOOO     M   M     M   M     EEEEE    W E S T
!    P    P   O    O   M M M M   M M M M   E              
!    PPPPP    O    O   M  M  M   M  M  M   EEEE         A F R I C A N
!    P        O    O   M     M   M     M   E                 
!    P         OOOO    M     M   M     M    EEEEE          M O N S O O N
! 
!    -                                -               -        -     -
!    P O S T - P R O C E S S I N G    O F    T H E    M A R    M O D E L
!    -                                -               -        -     -
!
!    Version 1.2
!    12.09.2005 (17.01.2006)
!
!
!
!
! This file contains all variable declarations for the program POMME.f.
!


!function num_2_string declaration


      character(len=6) :: num_2_string


!files


      integer(4)                  :: output_id !output_id
      integer(4), dimension(12,4) :: input_id  !input_id(month,file)


!dimensions


      integer(4), dimension(0:5) :: dim_ID
      integer(4), dimension(0:5) :: dim_len


!variables


      integer(4), dimension(0:4)        :: var_id

      real(8), dimension(40)            :: sigma
      real(8), dimension(1)             :: timeaxis   !t

      real(8), dimension(0:5,124)       :: values1d   !x, y, z, d
      real(8), dimension(124,124)       :: values2d   !xy
      real(8), dimension(124,124,40,48) :: values4d   !xyt, xyzt
      real(8), dimension(124,124,40)    :: daily4d    !xy, xyz
      real(8), dimension(124,124,40)    :: daily4d2


!dummy


      integer(4)               :: status, nerr
      integer(4), dimension(4) :: start, count

      character(len=8),  dimension(0:5)     :: name
      character(len=50), dimension(0:5)     :: att_units, att_longname
      real(8),           dimension(0:5,1:2) :: att_range

      character(len=3)                  :: run
      character(len=4)                  :: year
      character(len=2), dimension(2)    :: cmonth_slot
      character(len=2), dimension(12)   :: month, cday
      character(len=3), dimension(12)   :: month_name
      integer(4),       dimension(0:12) :: day, daycum
      integer(4),       dimension(2)    :: month_slot
      logical                           :: first_month

      integer(4)                     :: i, j, k, m, n, o, p, t, tt, d,
     &                                  t1, t2, t3, time_index
      integer(4), dimension(0:5)     :: i5tmp
      real(8)                        :: rmin, rmax
      character(len=3), dimension(4) :: c3tmp
      character(len=50)              :: ctmp50

