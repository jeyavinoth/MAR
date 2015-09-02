C +---------------------------------------------------------------------+
C | MAPOST                                                    05-1999 |
C |   SubRoutine dBnd2DSTA                                              |
C |   Computes the "same dist to boundary" average, + time average,     |
C |   and writes result                                                 |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE dBnd2DSTA (iOUTnc,istat, var, idt1,
     &                      dbNam,tM_var,taksqr)


      IMPLICIT NONE
 
C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
            
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'globals.inc'

C +...* Input:
      INTEGER iOUTnc, istat, idt1
      REAL var(mx,my)   
      CHARACTER*(*) dbNam
      LOGICAL taksqr

C +...* Input:
      REAL tM_var(ndb)

C +...* Indexes...
      INTEGER idb,ipe 


C +--Reset statistics
C +  ================

      IF (istat.EQ.0) THEN
        DO idb = 1,ndb
           tM_var(idb)= 0.0
        END DO
      ENDIF 

C +--Update statistics.
C +  ==================
C +
C +  Perimeter is devided in 2x2 lines of same lenght:
C +      mx-1 and my-1:
C +
C +  (1)<---------->(mx-1)  (1)    top
C +  (2)                     ^
C +   ^                      +
C +   +                      +
C +   +                   (my-1)
C + (my) (2)<------------->(mx)    bottom
C +
      DO idb = 1,ndb
         
C +      * Top and bottom contributions:
         DO ipe = idb, mx-idb
            tM_var(idb)= tM_var(idb)
     &        + var(ipe,idb) + var(ipe+1,my-idb+1)
         ENDDO
C +      * Left and Right contributions:
         DO ipe = idb, my-idb
            tM_var(idb)= tM_var(idb)
     &        + var(idb,ipe+1) + var(mx-idb+1,ipe)
         ENDDO
      
      ENDDO

C +--Final computation & write result (when istat = 2)
C +  =================================================
      IF (istat.EQ.2) THEN

C +   *Devide by the perimeter AND #time steps (idt1):
C
        DO idb = 1,ndb
             tM_var(idb)= tM_var(idb)
     &        / REAL ((mx+my-4*idb+2)*2*idt1)
             IF (taksqr .AND. tM_var(idb).GE.0.) THEN 
                tM_var(idb)= SQRT(tM_var(idb))
             ENDIF
        ENDDO

        CALL UNwrite (iOUTnc, dbNam,  0, ndb, 1, 1, tM_var)
C +     ...all data written at the last call to this routine
      ENDIF

      RETURN
      END
