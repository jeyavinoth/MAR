C +---------------------------------------------------------------------+
C | MAPOST                                                     3.0.0    |
C |   SubRoutine HdynSTA2D                                              |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE HDynSTA2D(iOUTnc,istat,idt1,var,
     &                    tM_nam,tSD_nam,tMEAN,t__SD)

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
      CHARACTER*(*) tM_nam,tSD_nam

C +...* Input/Output:
      REAL tMEAN(mx,my),t__SD(mx,my)

C +...* Indexes...
      INTEGER ii,jj,lv

C     write (*,*) 'HDynSTA2D for ',tM_nam, istat
 
C +--Reset statistics
C +  ================

      IF (istat.EQ.0) THEN
         DO jj= 1, my
         DO ii= 1, mx
             tMEAN(ii,jj)= 0.0
             t__SD(ii,jj)= 0.0
         ENDDO
         ENDDO         
      ENDIF 

C +--Update statistics for each valid data
C +  =====================================
      DO jj= 1,my
      DO ii= 1,mx
         IF (var(ii,jj).LT.1.E25.AND.tMEAN(ii,jj).LT.1.E25)THEN
             tMEAN(ii,jj)= tMEAN(ii,jj)+ var(ii,jj)
             t__SD(ii,jj)= t__SD(ii,jj)+ var(ii,jj)*var(ii,jj)
         ELSE
             tMEAN(ii,jj)= 1.E30
             t__SD(ii,jj)= 1.E30
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
           t__SD(ii,jj)=      t__SD(ii,jj)/REAL(idt1)
     &                - tMEAN(ii,jj)*tMEAN(ii,jj) 

           IF(t__SD(ii,jj).GE.0.0)THEN             
             t__SD(ii,jj) = SQRT(t__SD(ii,jj))
           ELSE
             write(*,*) 'Anomalie :' 
             write(*,*) t__SD(ii,jj), ii,jj, idt1
             write(*,*) 'Variable en cause: ', tSD_nam
           ENDIF

         ENDIF
        ENDDO
        ENDDO
        CALL UNwrite (iOUTnc, tM_nam,  0, mx, my, 1, tMEAN)
        CALL UNwrite (iOUTnc, tSD_nam, 0, mx, my, 1, t__SD)
      ENDIF

      RETURN
      END
