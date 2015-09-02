C +---------------------------------------------------------------------+
C | MAPOST                                                      03-1999 |
C |   SubRoutine Cumul                                                  |
C |   Add "Period" (output time-step) to compute "cumulated" values.    |
C |                                      (e.g. precipitation)           |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE Cumul (Pvar,ni,nj,istat,Cvar)

      IMPLICIT NONE
 
C +...* Input:
      INTEGER istat, ni, nj
      REAL Pvar(ni,nj) 

C +...* Input/Output:
      REAL Cvar(ni,nj)

C +...* Indexes...
      INTEGER ii,jj    



      IF (istat.EQ.0) THEN

C +--Initialise Cumulated value.
C +  ===========================

        DO jj=1,nj
        DO ii=1,ni
           Cvar(ii,jj)= 0.0
        ENDDO
        ENDDO
C       *Note: First data is not added (start at 0)

      ELSE

C +--Update Cumulated value.
C +  =======================

        DO jj=1,nj
        DO ii=1,ni
         Cvar(ii,jj)= Cvar(ii,jj) + Pvar(ii,jj)
        ENDDO
        ENDDO

      ENDIF 

      RETURN
      END
