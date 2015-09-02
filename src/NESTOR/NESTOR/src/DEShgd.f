C   +-------------------------------------------------------------------+
C   +  Subroutine DEShgd                            21/09/2004  NESTING +
C   +-------------------------------------------------------------------+
C   +                                                                   +
C   + Input : Parameters from MARgrd.ctr                                +
C   + ^^^^^^^                                                           +
C   +                                                                   +
C   + Output: Creation of the horizontal grid of MAR                    +
C   + ^^^^^^^ Variables : NST__x(mx,my) and NST__y(mx,my)  (long./lat.) +
C   +                     NSTgdx(mx)    and NSTgdy(my)     (Lambert)    +
C   +                     NST_dx (horizontal resolution)                +
C   +                                                                   +
C   +-------------------------------------------------------------------+


      SUBROUTINE DEShgd


      IMPLICIT NONE


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'NSTvar.inc'
      INCLUDE 'LSCvar.inc'      
      INCLUDE 'DESvar.inc'
      INCLUDE 'MARvar.inc'

C +---Local variables
C +   ---------------

      INTEGER i,j,fID,iloc,jloc

      REAL degrad,MinLon,MaxLon,MinLat,MaxLat,lwblon,upblon,
     .     lwblat,upblat,empty1(1),dist,distmin,LSC_dx,LSC_dy

      CHARACTER*7   namlon,namlat,nam_SH
      CHARACTER*10  var_units
      CHARACTER*100 LSCtit

                         
C +---Constants
C +   ---------

      DATA degrad /  1.745329252d-2/
C +...     degrad : Conversion Factor: Radian --> Degrees


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---READING OF GRID PARAMETERS IN MARgrd.ctr
C +   ========================================

      OPEN (unit=51,status='old',file='MARgrd.ctr')

       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) maptyp
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) GElon0
       read (51,*) imez  
       read (51,*) GElat0
       read (51,*) jmez  
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) dx
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) GEddxx
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) ptopDY
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) zmin
       read (51,*) aavu
       read (51,*) bbvu
       read (51,*) ccvu
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,'(l4)') vertic
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) sst_SL
       read (51,*) !- - - - - - - - - - - - - - - - - -

      CLOSE(unit=51)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---HORIZONTAL RESOLUTION
C +   =====================

      dx = dx * 1000.
      dy = dx

      NST_dx=dx


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CREATION OF HORIZONTAL MAR GRID
C +   ===============================


C +---Domain reference grid point
C +   ---------------------------

      IF (imez.le.0.or.imez.gt.mx) imez = mx/2
      IF (jmez.le.0.or.jmez.gt.my) jmez = my/2


C +---Compute interpolated horizontal grid
C +   ------------------------------------

      OPEN (unit=52,status='old',file='LSCfil.dat')
      READ (52,'(a100)',END=230) LSCfil
      GOTO 240

230   write(6,*) 'No file found in LSCfil.dat.'
      write(6,*) 'STOP in DEShgd.f'
      STOP

240   CONTINUE
      CLOSE(unit=52)

C +        *******
      CALL UNropen (LSCfil,fID,LSCtit)
C +        *******


C +---Screen message
C +   --------------

      write(6,*) 'Map projection: interpolation from LSC grid'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*) 'Open file  : ',LSCfil


C +---Read LSC (input) grid
C +   ---------------------

      lwblon =   -400.0
      upblon =    400.0
      lwblat =   -100.0
      upblat =    100.0

      IF (LSCmod.eq.'MAR') THEN
       namlon='lon'
       namlat='lat'
       nam_SH='sh'
      ELSE
       IF (LSCmod.eq.'LMD') THEN
        namlon='nav_lon'
        namlat='nav_lat'
        nam_SH='phis'
       ELSE
        IF (LSCmod.eq.'NCP') THEN
         namlon='lon'
         namlat='lat'
         nam_SH='SH'
        ELSE
         namlon='lon'
         namlat='lat'
         nam_SH='SH'
        ENDIF
       ENDIF
      ENDIF

      IF (REGgrd) THEN

C +         ******
       CALL UNread (fID,nam_SH,1,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC_sh)
C +         ******

       DO j=1,nj
       DO i=1,ni
        LSC__x(i,j)=LSC1Dx(i)
        LSC__y(i,j)=LSC1Dy(j)
       ENDDO
       ENDDO

      ELSE

C +         ******
       CALL UNread (fID,namlon,1,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__x)
C +         ******
       CALL VALchk (namlon,ni,nj,LSC__x,lwblon,upblon)
C +         ******
       CALL UNread (fID,namlat,1,1,bi,bj,ni,nj,1,
     .              LSC1Dx,LSC1Dy,empty1,var_units,LSC__y)
C +         ******
       CALL VALchk (namlat,ni,nj,LSC__y,lwblat,upblat)
C +         ******

      ENDIF


      IF (LSCmod.eq.'MAR') THEN
       LSC_dx = (LSC1Dx(2)-LSC1Dx(1))*1000.
       LSC_dy = (LSC1Dy(2)-LSC1Dy(1))*1000.
      ELSE
       LSC_dx = ABS(LSC__x(ni/2+1,nj/2)-LSC__x(ni/2,nj/2))
     .        * 111111.11111*ABS(COS(LSC__y(ni/2,nj/2)/degrad))
       LSC_dy = ABS(LSC__y(ni/2,nj/2+1)-LSC__y(ni/2,nj/2))
     .        * 111111.11111
      ENDIF


C +---Close the NetCDF file
C +   =====================

C +        *******
      CALL UNclose (fID)
C +        *******


C +---Factor between LSC and NST resolutions
C +   --------------------------------------

      fdiv = NINT(MAX(LSC_dx,LSC_dy)/NST_dx)


C +---Compute new horizontal resolution
C +   ---------------------------------

      NST_dx = LSC_dx / REAL(fdiv)


C +---Simple grid (Lambert coordinates)
C +   ---------------------------------

      DO i=1,mx
       NSTgdx(i)=(i-imez)*LSC_dx/REAL(fdiv)/1000.
      ENDDO

      DO j=1,my
       NSTgdy(j)=(j-jmez)*LSC_dy/REAL(fdiv)/1000.
      ENDDO


C +---Compute NST grid
C +   ----------------

      iloc = -1
      jloc = -1
      distmin = 1.d+30

      DO j=1,nj
      DO i=1,ni
       dist = SQRT( (LSC__x(i,j)-GElon0)**2.
     .      +       (LSC__y(i,j)-GElat0)**2. )
       IF (dist.lt.distmin) THEN
        distmin = dist
        iloc = i     ! This point will correspond to imez in NST
        jloc = j     ! This point will correspond to jmez in NST
       ENDIF
      ENDDO
      ENDDO

      IF (iloc.eq.-1.or.jloc.eq.-1) THEN
       write(6,*) 'The center of NST grid is not included in LSC grid.'
       write(6,*) 'STOP in DEShgd.f'
       STOP
      ENDIF


C +---Define transfer of index
C +   ------------------------

      DO j=1,my
      DO i=1,mx
       iiL2N(i,j) = ((i-imez+(imez*fdiv))/fdiv) - imez + iloc
       jjL2N(i,j) = ((j-jmez+(jmez*fdiv))/fdiv) - jmez + jloc
       IF (iiL2N(i,j).lt.1.or.iiL2N(i,j).gt.ni.or.
     .     jjL2N(i,j).lt.1.or.jjL2N(i,j).gt.nj) THEN
        write(6,*) 'The NST grid is outside the LSC grid.'
        write(6,*) 'Please check the definition of the NST domain.'
        write(6,*) i,j,iiL2N(i,j),jjL2N(i,j)
        write(6,*) 'STOP in DEShgd.f'
        STOP
       ENDIF
      ENDDO
      ENDDO


C +---Create NST grid
C +   ---------------

      DO j=1,my
      DO i=1,mx

       auxL2N(i,j) = REAL(MOD(i-imez+(imez*fdiv),fdiv)) / REAL(fdiv)
       auyL2N(i,j) = REAL(MOD(j-jmez+(jmez*fdiv),fdiv)) / REAL(fdiv)

       NST__x(i,j) = auxL2N(i,j)
     .             * (   auyL2N(i,j)*LSC__x(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__x(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j)*LSC__x(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__x(iiL2N(i,j)  ,jjL2N(i,j)  ))
                                       
       NST__y(i,j) = auxL2N(i,j)
     .             * (   auyL2N(i,j)*LSC__y(iiL2N(i,j)+1,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__y(iiL2N(i,j)+1,jjL2N(i,j)  ))
     .             + (1.0-auxL2N(i,j))
     .             * (   auyL2N(i,j)*LSC__y(iiL2N(i,j)  ,jjL2N(i,j)+1)
     .              +(1-auyL2N(i,j))*LSC__y(iiL2N(i,j)  ,jjL2N(i,j)  ))

      ENDDO
      ENDDO


C +---Compute horizontal extent of the horizontal domain
C +   --------------------------------------------------

      MinLon = NST__x(1,1)
      MaxLon = NST__x(1,1)
      MinLat = NST__y(1,1)
      MaxLat = NST__y(1,1)
      DO j=1,my
      DO i=1,mx
       MinLon = MIN(NST__x(i,j),MinLon)
       MaxLon = MAX(NST__x(i,j),MaxLon)
       MinLat = MIN(NST__y(i,j),MinLat)
       MaxLat = MAX(NST__y(i,j),MaxLat)
      ENDDO
      ENDDO


C +---Print the characteristics of the horizontal grid
C +   ------------------------------------------------

      write(6,*) 'Horizontal MAR grid created'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,200) mx,my,LSC_dx/REAL(fdiv)/1000.,
     .             LSC_dy/REAL(fdiv)/1000.,GEddxx,
     .             MinLon,MaxLon,MinLat,MaxLat
200   format(' Grid points           : ',i4,' * ',i4,/,
     .       ' Horizontal resol. (X) : ',f7.1,' km.',/,
     .       ' Horizontal resol. (Y) : ',f7.1,' km.',/,
     .       ' Domain orientation    : ',f7.0,' deg.',/,
     .       ' MAR longitude between : ',f7.2,' and ',f7.2,/,
     .       ' MAR latitude  between : ',f7.2,' and ',f7.2,/)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END

