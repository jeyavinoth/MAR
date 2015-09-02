C   +-------------------------------------------------------------------+
C   |  Subroutine GEOpot                               June 99  NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   | This subroutine contains the integration of the hydrostatic       |
C   | relation.                                                         |
C   |                                                                   |
C   | Input : - NST__t : real temperature     (K)                       |
C   | ^^^^^^^ - NST_sh : surface height       (m)                       |
C   |         - NST_sp : surface pressure     (kPa)                     |
C   |         - NST__p : pressure             (kPa)                     |
C   |                                                                   |
C   | Output: - NST_zz : levels height                                  |
C   | ^^^^^^^                                                           |
C   +-------------------------------------------------------------------+

      SUBROUTINE GEOpot (NST__t,NST_sh,NST_sp,NST__p,NST_zz)


      IMPLICIT NONE


C +--General Variables
C +  =================

      INCLUDE 'NSTdim.inc'

      INTEGER i,j,k

      REAL    ab,ra,gravit

      REAL    NST__t(mx,my,mz),NST__p(mx,my,mz),NST_zz(mx,my,mz),
     .        NST_sh(mx,my   ),NST_sp(mx,my   )


C +--Constants
C +  =========

      DATA ra     / 287.d0 /
C +...     ra     : Perfect Gas Law  Constant (J/kg/K)

      DATA gravit / 9.81d0 /
C +...     gravit : Gravity constant

C     'WARNING - GEOpot: you may consider using VERhyd'

C +---Integration of the Hydrostatic Equation
C +   =======================================

      DO j=1,my
      DO i=1,mx

       NST_zz(i,j,mz)=NST_sh(i,j)
     .               +(ra/gravit)*NST__t(i,j,mz)
     .               *LOG(NST_sp(i,j)/NST__p(i,j,mz))

       DO k=mz-1,1,-1

        NST_zz(i,j,k)=NST_zz(i,j,k+1) 
     .               +(ra/gravit)*0.5*(NST__t(i,j,k)+NST__t(i,j,k+1))
     .               *LOG(NST__p(i,j,k+1)/NST__p(i,j,k))

       ENDDO

C +... z1 = z0 - (RT/g) ln(p1/p0)
C +...    = z0 + (RT/g) ln(p0/p1)
      
      ENDDO 
      ENDDO 


      RETURN
      END
