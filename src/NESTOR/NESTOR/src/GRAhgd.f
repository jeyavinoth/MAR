C   +-------------------------------------------------------------------+
C   +  Subroutine GRAhgd                         February 2002  NESTING +
C   +-------------------------------------------------------------------+
C   +                                                                   +
C   + Input : Parameters from MARgrd.ctr                                +
C   + ^^^^^^^                                                           +
C   +                                                                   +
C   + Output: Creation of the horizontal grid for GRADS                 +
C   + ^^^^^^^ Variables : NST__x(mx,my) and NST__y(mx,my)  (long./lat.) +
C   +                     NSTgdx(mx)    and NSTgdy(my)     (Lambert)    +
C   +                     NST_dx (horizontal resolution)                +
C   +                                                                   +
C   +-------------------------------------------------------------------+


      SUBROUTINE GRAhgd


      IMPLICIT NONE


C +---General variables
C +   ---------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NESTOR.inc'
      INCLUDE 'NSTvar.inc' 
      INCLUDE 'LSCvar.inc'

C +---Local variables
C +   ---------------

      INTEGER i,j,fID,iloc,jloc,vmx1,vmx2,vmy1,vmy2,imez,jmez

      REAL degrad,MinLon,MaxLon,MinLat,MaxLat,lwblon,upblon,
     .     lwblat,upblat,empty1(1),dist,distmin,GElon0,GElat0,
     .     dx,dy,resol

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

      OPEN (unit=51,status='old',file='GRAgrd.ctr')

       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) GElon0
       read (51,*) imez  
       read (51,*) GElat0
       read (51,*) jmez  
       read (51,*) !- - - - - - - - - - - - - - - - - -
       read (51,*) resol
       read (51,*) !- - - - - - - - - - - - - - - - - -

      CLOSE(unit=51)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---HORIZONTAL RESOLUTION (degree)
C +   =====================

      dx = resol*111.111111
      dy = dx

      NST_dx=dx


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C +---CREATION OF HORIZONTAL MAR GRID
C +   ===============================


C +---Domain reference grid point
C +   ---------------------------

      IF (imez.lt.0.or.imez.gt.mx) imez = mx/2
      IF (jmez.lt.0.or.jmez.gt.my) jmez = my/2


C +---Create NST grid
C +   ---------------

      DO i=1,mx
       NSTgdx(i)=GElon0+(i-imez)*resol
      ENDDO

      DO j=1,my
       NSTgdy(j)=GElat0+(j-jmez)*resol
      ENDDO

      DO j=1,my
      DO i=1,mx
       NST__x(i,j) = NSTgdx(i)
       NST__y(i,j) = NSTgdy(j)
      ENDDO
      ENDDO


C +---Open LSC file in order to verify the inclusion in LSC grid
C +   ----------------------------------------------------------

      OPEN (unit=52,status='old',file='LSCfil.dat')
      READ (52,'(a100)',END=230) LSCfil
      GOTO 240

230   write(6,*) 'No file found in LSCfil.dat.'
      write(6,*) 'STOP in DEShgd.f'
      STOP

240   CONTINUE
      CLOSE(unit=52)


C +---Screen message
C +   --------------

      write(6,*) 'Map projection: regular grid included in LSC grid'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,*) 'Open file  : ',LSCfil


C +        *******
      CALL UNropen (LSCfil,fID,LSCtit)
C +        *******


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


C +---Close the NetCDF file
C +   =====================

C +        *******
      CALL UNclose (fID)
C +        *******


C +---Verify the inclusion of NST grid in LSC grid
C +   --------------------------------------------

      MinLon = LSC__x( 1,1)
      MaxLon = LSC__x(ni,1)
      MinLat = LSC__y(1, 1)
      MaxLat = LSC__y(1,nj)
      DO i=1,ni
       MinLat = MAX(LSC__y(i, 1),MinLat)
       MaxLat = MIN(LSC__y(i,nj),MaxLat)
      ENDDO
      DO j=1,nj
       MinLon = MAX(LSC__x( 1,j),MinLon)
       MaxLon = MIN(LSC__x(ni,j),MaxLon)
      ENDDO

      IF (GElon0.lt.MinLon .or.
     .    GElon0.gt.MaxLon .or.
     .    GElat0.lt.MinLat .or.
     .    GElat0.gt.MaxLat) THEN
       write(6,*)
       write(6,*) 'The center of the NST grid is not included'
       write(6,*) 'in the LSC grid. Please check and modify'
       write(6,*) 'the GRAgrd.ctr file.'
       write(6,*)
       write(6,*) '--> STOP in GRAhgd.f'
       write(6,*)
       STOP
      ENDIF

      IF (MinLon.gt.NSTgdx( 1) .or.
     .    MaxLon.lt.NSTgdx(mx) .or.
     .    MinLat.gt.NSTgdy( 1) .or.
     .    MaxLat.lt.NSTgdy(my)) THEN
       vmx1 = INT((GElon0-MinLon)/resol)
       vmx2 = INT((MaxLon-GElon0)/resol)
       vmy1 = INT((GElat0-MinLat)/resol)
       vmy2 = INT((MaxLat-GElat0)/resol)
       write(6,*)
       write(6,*) 'NST grid is not fully included in LSC grid'
       write(6,*)
       write(6,*) 'Characteristics of the LSC grid :'
       write(6,*) '- MinLon = ',MinLon
       write(6,*) '- MaxLon = ',MaxLon
       write(6,*) '- MinLat = ',MinLat
       write(6,*) '- MaxLat = ',MaxLat
       write(6,*)
       write(6,*) 'Characteristics of the NST grid :'
       write(6,*) '- MinLon = ',NSTgdx(1)
       write(6,*) '- MaxLon = ',NSTgdx(mx)
       write(6,*) '- MinLat = ',NSTgdy(1)
       write(6,*) '- MaxLat = ',NSTgdy(my)
       write(6,*)
       write(6,*) 'Please try the following parameters : '
       write(6,*)
       write(6,*) 'mx   (in NSTdim.inc) = ',vmx1+vmx2
       write(6,*) 'my   (in NSTdim.inc) = ',vmy1+vmy2
       write(6,*) 'imez (in GRAgrd.ctr) = ',vmx1
       write(6,*) 'jmez (in GRAgrd.ctr) = ',vmy1
       write(6,*)
       write(6,*) '--> STOP in GRAhgd.f'
       write(6,*)
       STOP
      ENDIF


C +---Print the characteristics of the horizontal grid
C +   ------------------------------------------------

      write(6,*) 'Horizontal regular grid created'
      write(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(6,200) mx,my,resol,
     .             MinLon,MaxLon,MinLat,MaxLat
200   format(' Grid points             : ',i4,' * ',i4,/,
     .       ' Horizontal resolution   : ',f7.2,' degree',/,
     .       ' GRADS longitude between : ',f7.2,' and ',f7.2,/,
     .       ' GRADS latitude  between : ',f7.2,' and ',f7.2,/)


C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      RETURN
      END
