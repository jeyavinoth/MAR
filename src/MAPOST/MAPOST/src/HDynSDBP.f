C +---------------------------------------------------------------------+
C | MAPOST                                                      3.0.0   |
C |   SubRoutine HDynSDBP                                               |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE HDynSDBP(iOUTnc,istat,idt1,jhuCUR,Var,
     &                    SDBP_nam,ipfl,FilBf,MBPv,SDBPv)

      IMPLICIT NONE
 
C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

C +---MAPOST parameters & variables
C +   -----------------------------
      include 'MAPOST.inc'
      include 'globals.inc'
      

C +...* Numerical filter: 
      INCLUDE 'HDynSDBP.inc'

C +...* Input:
      INTEGER iOUTnc, istat, idt1, jhuCUR, ipfl
      REAL Var(mx,my), FilBf(mx,my,nfl1)
      CHARACTER*(*) SDBP_nam

C +...* Input/Output:
      REAL MBPv(mx,my),SDBPv(mx,my)
      
C +...* Local:
      REAL FilVal
      
C +...* Indexes...
      INTEGER ii,jj,ifl, ntfv
      
      
C +--Reset statistics
C +  ================

      IF (istat.EQ.0) THEN
         DO jj= 1, my
         DO ii= 1, mx
           MBPv (ii,jj)= 0.0
           SDBPv(ii,jj)= 0.0
           DO ifl = 1, nfl1
             FilBf(ii,jj,ifl)= 0.0
           ENDDO
C +        *Initialise the filter buffer 
C +         (Only ifl=nfl+1 is important)
         ENDDO
         ENDDO
         ipfl = 0         
      ENDIF 

C +--Time-filter the data, compute SDBP
C +  ==================================
C +   Filter is prepared for 12 hour interval, so  
C +   enter at appropriate time only: 
      IF (jhuCUR.EQ.0 .OR. jhuCUR.EQ.12) THEN

C +..   Count #data loaded in the filter:
        ipfl= ipfl+1

        DO jj= 1,my
        DO ii= 1,mx
          IF (var(ii,jj).LT.1.E20.AND.MBPv(ii,jj).LT.1.E29)THEN
          
C +..      * Load data in the filter-buffer:
           DO ifl = 1, nfl
             FilBf(ii,jj,ifl)= FilBf(ii,jj,ifl+1)
     .         +  FilCo (ifl) * var(ii,jj)           
           ENDDO
           
C +..      * After nfl steps, filtered data is avail.to calc tSD:
           IF (ipfl.GE.nfl) THEN
              FilVal= FilBf(ii,jj,1)
C +..         ^^Filtered value (at current step - nfl/2)       
              MBPv(ii,jj) = MBPv(ii,jj) + FilVal
              SDBPv(ii,jj)= SDBPv(ii,jj)+ FilVal*FilVal           
           ENDIF 
             
          ELSE
           MBPv(ii,jj) = 1.E30
           SDBPv(ii,jj)= 1.E30
C +        ..Missing Value
          ENDIF
        ENDDO
        ENDDO
      
      ENDIF

C +--Final computation & write result (when istat = 2)
C +  =================================================
      IF (istat.EQ.2) THEN
        ntfv = ipfl+1-nfl
C +..   ^^Number of avail. values after filt.: ipfl+1-nfl
        IF (ntfv.LT.1) THEN 
          write(*,*) 'HDynSDBP : not enough input data.'
          GOTO 999
        ENDIF
        DO jj= 1,my
        DO ii= 1,mx
        
         IF (MBPv(ii,jj).LT.1.E29)THEN
           MBPv (ii,jj)= MBPv(ii,jj)/REAL(ntfv)
           SDBPv(ii,jj)= SDBPv(ii,jj)/REAL(ntfv)
     &                  -MBPv(ii,jj)*MBPv(ii,jj) 

           IF(SDBPv(ii,jj).GT.0.0)THEN             
             SDBPv(ii,jj) = SQRT(SDBPv(ii,jj))
           ELSE
             write(*,*) 'Anomalie dans HDynSDBP :' 
             write(*,*) SDBPv(ii,jj), ii,jj, ntfv
           ENDIF
         ENDIF

        ENDDO
        ENDDO
        CALL UNwrite (iOUTnc, SDBP_nam,  0, mx, my, 1, SDBPv)
      ENDIF

 999  RETURN
      END
