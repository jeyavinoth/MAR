C +---------------------------------------------------------------------+
C | MAPOST                                                    3.0.0     |
C |   SubRoutine HdynSTApl                                              |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE HDynSTApl(iOUTnc,istat,idt1,ipl,Var,
     &                    tM_nam,tSD_nam,tMEAN,t__SD)

      IMPLICIT NONE
 
C +---LS and MAR domain dimensions :
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

C +---MAPOST parameters & variables
      include 'MAPOST.inc'
      include 'globals.inc'
      

C +...* Input:
      INTEGER iOUTnc, istat, idt1,ipl
      REAL Var(mx,my)
      CHARACTER*(*) tM_nam,tSD_nam

C +...* Input/Output:
      REAL tMEAN(mx,my,npl),t__SD(mx,my,npl)

C +...* Indexes...
      INTEGER ii,jj,lv

C +--Reset statistics
C +  ================

      IF (istat.EQ.0) THEN
         DO jj= 1, my
         DO ii= 1, mx
             tMEAN(ii,jj,ipl)= 0.0
             t__SD(ii,jj,ipl)= 0.0
         ENDDO
         ENDDO         
      ENDIF 

C +--Update statistics for each valid data
C +  =====================================
      DO jj= 1,my
      DO ii= 1,mx
         IF(   (ABS(var(ii,jj)).LT.1.E15)
     .         .AND.(tMEAN(ii,jj,ipl).LT.1.E15)  ) THEN
             tMEAN(ii,jj,ipl)= tMEAN(ii,jj,ipl)+ var(ii,jj)
             t__SD(ii,jj,ipl)= t__SD(ii,jj,ipl)+ var(ii,jj)*var(ii,jj)
         ELSE
             tMEAN(ii,jj,ipl)= 1.E30
             t__SD(ii,jj,ipl)= 1.E30
C +          ..Missing Value
         ENDIF

      ENDDO
      ENDDO

C +--Final computation & write result (when istat = 2)
C +  =================================================
      IF (istat.EQ.2) THEN
        DO jj= 1,my
        DO ii= 1,mx
         IF (tMEAN(ii,jj,ipl).LT.1.E10)THEN

           tMEAN(ii,jj,ipl)= tMEAN(ii,jj,ipl)/REAL(idt1)
           t__SD(ii,jj,ipl)=      t__SD(ii,jj,ipl)/REAL(idt1)
     &                - tMEAN(ii,jj,ipl)*tMEAN(ii,jj,ipl) 

           IF(t__SD(ii,jj,ipl).GT.0.0)THEN             
              t__SD(ii,jj,ipl) = SQRT(t__SD(ii,jj,ipl))
           ELSE
              write(*,*) 'Anomalie :' 
              write(*,*) t__SD(ii,jj,ipl),ii,jj,idt1,ipl
              write(*,*) 'Variable en cause: ', tSD_nam
           ENDIF

         ENDIF
        ENDDO
        ENDDO
        IF (ipl.EQ.npl) THEN
           CALL UNwrite (iOUTnc, tM_nam, 0, mx, my, npl, tMEAN)
           CALL UNwrite (iOUTnc, tSD_nam,0, mx, my, npl, t__SD)
           write(*,*) 'HdynSTApl completed for ',tM_nam
        ENDIF
      ENDIF

      RETURN
      END
