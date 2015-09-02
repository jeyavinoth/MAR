C +---------------------------------------------------------------------+
C | MAPOST                                                     3.1.0    |
C |   SubRoutine Mean2D                                                 |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE    Mean2D(iOUTnc,istat,idt1,var,
     &                    tM_nam,tMEAN)

      IMPLICIT NONE
 
C +---LS and MAR domain dimensions :
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

C +---MAPOST parameters & variables
      include 'MAPOST.inc'
      include 'globals.inc'
      
C +...* Input:
      INTEGER iOUTnc, istat, idt1
      REAL var(mx,my)
      CHARACTER*(*) tM_nam

C +...* Input/Output:
      REAL tMEAN(mx,my)

C +...* Indexes...
      INTEGER ii,jj,lv

C     write (*,*) 'Mean2D for ',tM_nam, istat
 
C +--Reset statistics
C +  ================

      IF (istat.EQ.0) THEN
         DO jj= 1, my
         DO ii= 1, mx
             tMEAN(ii,jj)= 0.0
         ENDDO
         ENDDO         
      ENDIF 

C +--Update statistics for each valid data
C +  =====================================
      DO jj= 1,my
      DO ii= 1,mx
         IF (var(ii,jj).LT.1.E25.AND.tMEAN(ii,jj).LT.1.E25)THEN
             tMEAN(ii,jj)= tMEAN(ii,jj)+ var(ii,jj)
         ELSE
             tMEAN(ii,jj)= 1.E30
C +          ..Missing Value
         ENDIF

      ENDDO
      ENDDO

C +--Final computation & write result (when istat = 2)
C +  =================================================
      IF (istat.EQ.2) THEN

        DO jj= 1,my
        DO ii= 1,mx
         IF (tMEAN(ii,jj).LT.1.E25)THEN
           tMEAN(ii,jj)= tMEAN(ii,jj)/REAL(idt1)
         ENDIF
        ENDDO
        ENDDO
        CALL UNwrite (iOUTnc, tM_nam,  0, mx, my, 1, tMEAN)
      ENDIF

      RETURN
      END
