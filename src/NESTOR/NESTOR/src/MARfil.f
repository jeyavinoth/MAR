C   +-------------------------------------------------------------------+
C   |  Subroutine MARfil                               June 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | Input : my,dx,dt : grid size, horizontal resolution and time step |
C   | ^^^^^^^                                                           |
C   |                                                                   |
C   | Output: - deltat : Implicit Filter Parameter (Temperature)        |
C   | ^^^^^^^ - deltau : Implicit Filter Parameter (Wind Speed)         |
C   |         - deltap : Implicit Filter Parameter (Pressure)           |
C   |         - akhdel : Horizontal Diffusion Coefficient               |
C   |         - akhfac : Horiz.vKar**2 (horizontal diffusion)           |
C   |         - akhmax : Upper sponge                                   |
C   |                                                                   |
C   +-------------------------------------------------------------------+

      subroutine MARfil(my,dx,dt,deltat,deltau,deltap,
     .                           akhdel,akhfac,akhmax)

      implicit none

      integer my,mmy

      real dx,dy,dt,deltat,deltau,deltap,akhdel,akhfac,akhmax,akh

      mmy = my
      dy  = dx


C + Horizontal Diffusion Parameters
C + -------------------------------

      if (mmy.eq.1) then
       call rdelta(dx,dt,deltat,akhdel,akh)
      else
       call sdelta(dx,dt,deltat,akhdel,akh)
      end if

      deltau = deltat
      deltap = deltat
      akhfac = 0.16
      akhmax = 0.1*dx*dx/dt
C...  Vertical upper sponge
      
      return
      end


C +------------------------------------------------------------------------+
C |   SubRoutine rdelta is used to define the horizontal 1-D filter in MAR |
C +------------------------------------------------------------------------+

      subroutine rdelta(dx,dt,delta,akhdel,akh)

      implicit none

      character*1 schema

      integer i

      real pi,al,dx,dt,delta,akhdel,akh,ak,ckx,amorf,akhp,clx,
     .     slx,r1,amin,amax,alpha,beta,r,amor,alph,akhmn

      data schema/'E'/

      data pi    / 3.1415926567 /
           delta = 0.05

      do i=3,21,3

       al    = i * dx

       ak    = dx *2. *pi /al
       ckx   = cos(ak)

       amorf = (ckx + 1) / ((1-delta) *ckx + 1 + delta)

       akhp  = -(dx*dx/(dt *ak *ak)) *alog(amorf)

       clx   = (cos(ak/2))**2
       slx   = (sin(ak/2))**2
       r1    =  dx**2/dt
 
       amin  = 1.e30
       amax  = 0.

       if (schema.eq.'I') then
        alpha = 0.25
       else
        alpha = 1.
       end if

       akh   =  r1 *delta /(4. *(clx + delta *alpha *slx))

       beta  = 1. - alpha
       r     = akh  / r1
       amor  = (1 -4. *alpha *r *slx) / (1. +4. *beta * r *slx)
       amin  = amor
       alph  = alpha
       akhmn = akh
       alpha   = alph
       amor  = amin
       akh   = akhmn

      enddo

      akhdel = akh

      return
      end


C +------------------------------------------------------------------------+
C |   SubRoutine sdelta is used to define the horizontal 2-D filter in MAR |
C +------------------------------------------------------------------------+

      subroutine sdelta(dx,dt,delta,akhdel,akh)

      implicit none

      character*1 schema

      integer i

      real pi,al,dx,dt,delta,akhdel,akh,ak,ckx,amorf,akhp,clx,
     .     slx,r1,amin,amax,alpha,beta,r,amor,alph,akhmn,tlx,
     .     dtlx

      data schema/'E'/

      data pi    / 3.1415926567 /
           delta = 0.05

      do i=3,15,3

       al   = i * dx

       ak   = dx *2. *pi /al
       ckx  = cos(ak)

       amorf   = 1 / (1 + (delta*(1-ckx*ckx)
     .             + delta*delta*(ckx-1)*(ckx-1))
     .                  / ( (ckx+1)*(ckx+1) ) )

       clx   = (cos(ak/2))**2
       slx   = (sin(ak/2))**2
       tlx   = slx / clx
       dtlx  = 1 + delta *tlx
       r1    =  dx**2/dt

       amin  = 1.e30
       amax  = 0.
 
       if (schema.eq.'I') then
        alpha = 0.25
       else
        alpha = 1.
       end if

       akh   =  r1 *delta *dtlx /(4. *(clx + alpha *delta *slx *dtlx))

       beta  = 1. - alpha
       r     = akh  / r1
       amor  = (1 -4. *alpha *r *slx) / (1. +4. *beta * r *slx)
       amin  = amor
       alph  = alpha
       akhmn = akh
       alpha = alph
       amor  = amin
       akh   = akhmn

      enddo

      akhdel = akh

      return
      end
