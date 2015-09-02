C +---------------------------------------------------------------------+
C | MAPOST                                              xx-03-1999  MAR |
C |   SubRoutine InitZFilt.f                                            |
C +---------------------------------------------------------------------+
C | Get values for the Z500 filter (see HDynSDBP, and FltCfs.dat)    |
C +---------------------------------------------------------------------+
      SUBROUTINE InitZFilt                      

      IMPLICIT NONE

C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'      
 
C +...* Numerical filter: 
      INCLUDE 'HDynSDBP.inc'

      INTEGER ifl
      REAL*8 C_low, C_bp, C_high
      
      OPEN (unit=52,status='old',file='FltCfs.dat')
      DO ifl = 1, 10 
         READ(52,*)
      ENDDO
 
C +...* Read the filter coefs. from the table:
C +...  (coefs are symetrical, center index = n2fl+1  
C        nfl= 2*n2fl+1; currently, n2fl = 15 => 31 coefs)
C
      DO ifl = 0, n2fl
         READ(52,'(3F14.10)') C_low, C_bp, C_high
         FilCo(n2fl+1-ifl) = C_bp 
         FilCo(n2fl+1+ifl) = C_bp 
      ENDDO

C #   DO ifl = 1, nfl
C #     write(*,*) FilCo(ifl)
C #   ENDDO
C  
C #   write(*,*) ' << WARNING >> ? pass filter active'
         
      CLOSE(unit=52)

      RETURN
      END
