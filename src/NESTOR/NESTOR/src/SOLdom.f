C   +-------------------------------------------------------------------+
C   |  Subroutine SOLdom                             June 2003  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : NSTsol : surface type (ocean, ice, snow, land)            |
C   | ^^^^^^^ NSTtex : soil texture (fine, medium, coarse)              |
C   |                                                                   |
C   | Output: NST_d1 : surface heat capacity                            |
C   | ^^^^^^^ NSTalb : surface albedo                                   |
C   |         NSTeps : surface IR emissivity                            |
C   |         NST_z0 : roughness length (momentum)                      |
C   |         NST_r0 : roughness length (heat)                          |
C   |         NSTres : aerodynamic resistance                           |
C   |         NSTch0 : bulk aerodynamic coefficient air/surface         |
C   |                  humidity flux                                    |
C   |                                                                   |
C   +-------------------------------------------------------------------+


      SUBROUTINE SOLdom


      IMPLICIT NONE


C +---General and local variables
C +   ---------------------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'NESTOR.inc'

      INTEGER i,j

      REAL    zl_SL,zs_SL,zn_SL,argLAT,sinLAT,Amn,Rmin,zero,Shelfb,
     .        Rmax,rsurSL,SH1,SH2,minz0

C +---Data
C +   ----

      zero = 0.0
      Amn  = 0.1
      Rmin = 200.
      Rmax = 900.

      SH1  = 1500.  ! Up to this height   : normal vegetation
      SH2  = 2000.  ! Between SH1 and SH2 : decrease of vegetation


C +---Screen message
C +   --------------

      WRITE(6,*) 'Specification of surface characteristics'
      WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      WRITE(6,*)


C +---SURFACE CHARACTERISTICS
C +   =======================


C +---Typical Roughness Lengths (m) for land, sea, snow
C +   -------------------------------------------------

      zl_SL  = 1.00e-1
      zs_SL  = 1.00e-3
      zn_SL  = 1.00e-4
      minz0  = 1.00e-4


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CORRECTION OF PRESCRIBED SURFACE TYPE
C +   =====================================


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO j=1,my
      DO i=1,mx

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       GO TO (100,200,300,400,400) NSTsol(i,j)

C +---1. Ocean
C +   --------    

100    NST_d1(i,j) = 0.
       NSTalb(i,j) = 0.15
       NSTeps(i,j) = 0.97
       IF (NST_z0(i,j).lt.minz0) THEN
        NST_z0(i,j) = zs_SL
       ENDIF
       NST_r0(i,j) = 0.1*zs_SL
       NSTch0(i,j) = 0.00132
       NSTres(i,j) = 0.0
C +
       IF(region.eq."GRD".or.region.eq."ANT") THEN
       NST_d1(i,j) = 2.09d+8
       NSTalb(i,j) = 0.10d0
       ENDIF
C +
       GOTO 500

C +---2. Sea Ice
C +   ----------

200    NST_d1(i,j) = 1.05d+5
       NSTalb(i,j) = 0.85d00
       NSTeps(i,j) = 0.97d00
       IF (NST_z0(i,j).lt.minz0) THEN
        NST_z0(i,j) = zn_SL
       ENDIF
       NST_r0(i,j) = 0.1*zn_SL
       NSTch0(i,j) = 0.0021
C +... (Kondo and Yamazaki, 1990, JAM 29, p.376)
       NSTres(i,j) = 0.0

       IF(region.eq."GRD".or.region.eq."ANT") THEN
c      NSTalb(i,j) = 0.70d0
       NSTalb(i,j) = 0.10d0 ! To no have problem
                            !  when snow ice melt
       ENDIF
C +
       GOTO 500

C +---3. Snow Field
C +   -------------

300    NST_d1(i,j) = 1.05e+5
       NSTalb(i,j) = 0.85           
       NSTeps(i,j) = 0.97
       IF (NST_z0(i,j).lt.minz0) THEN
        NST_z0(i,j) = zn_SL
       ENDIF
       NST_r0(i,j) = 0.1*zn_SL
       NSTch0(i,j) = 0.0021
C +... (Kondo and Yamazaki, 1990, JAM 29, p.376)
       NSTres(i,j) = 0.0
       GOTO 500


C +---4. Continent  
C +   ------------

400    CONTINUE 

       IF (NSTtex(i,j).eq.1) THEN
        NST_d1(i,j) = 1.65e+5
        NSTalb(i,j) = 0.40
C +...  Dry Quartz Sand (Deardorff 1978 JGR p.1891)
       ELSE IF (NSTtex(i,j).eq.3) THEN
        NST_d1(i,j) = 7.55e+5
        NSTalb(i,j) = 0.15
C +...  Clay Pasture    (Deardorff 1978 JGR p.1891)
       ELSE
        NST_d1(i,j) = 2.88e+5
        NSTalb(i,j) = 0.25
C +...  O'Neill average (Deardorff 1978 JGR p.1891)
       ENDIF
C +
       argLAT      = 0.0628 *  NST__y(i,j)
       sinLAT      = max(zero,sin(argLAT))
C +
       NSTeps(i,j) = 0.97
C +
       NST_r0(i,j) = zl_SL
       NST_r0(i,j) =-0.9d-1*sinLAT + zl_SL
       NST_r0(i,j) = min(zl_SL,NST_r0(i,j))
       NST_r0(i,j) = 0.1000 *  NST_r0(i,j)
C +
       IF (NST_z0(i,j).lt.minz0) THEN
        NST_z0(i,j) = 10.0 * NST_r0(i,j)
       ENDIF
C +
       NSTch0(i,j) = 0.0025
C +
       NSTalb(i,j) = 3.0d-1 *  sinLAT       +0.5d-1
       NSTalb(i,j) = max(Amn  ,NSTalb(i,j))
C +
       NSTres(i,j) = 5.0d+3 *  sinLAT       -1.6d+3
       NSTres(i,j) = max(Rmin  ,NSTres(i,j))
       IF (NST_sh(i,j).le.SH1) THEN
        rsurSL  = Rmin
       ELSE
        IF (NST_sh(i,j).gt.SH2) THEN
         rsurSL = Rmax
        ELSE
         rsurSL = Rmin
     .          + (NST_sh(i,j)-SH1)/(SH2-SH1)*(Rmax-Rmin)
        ENDIF
       ENDIF
       NSTres(i,j) = max(rsurSL,NSTres(i,j))
C +
       IF(region.eq."GRD".or.region.eq."ANT") THEN
       NST_d1(i,j) = 1.74d+5
       NSTalb(i,j) = 0.2d0                   
       ENDIF
C +
       GOTO 500

500    CONTINUE


C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENDDO
      ENDDO

C +   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
