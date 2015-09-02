C +---------------------------------------------------------------------+
C | MAPOST                                                    3.0.0     |
C |   SubRoutine dBoundSTA                                              |
C |   Computes the "same dist to boundary" average, + time average,     |
C |   and writes result                                                 |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE dBoundSTA (iOUTnc,istat, var, kl, idt1,
     &                      dbNam,tM_var)


      IMPLICIT NONE
      
C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

C +---MAPOST parameters & variables
C +   -----------------------------
      include 'MAPOST.inc'
      include 'globals.inc'
      

C +...* Input:
      INTEGER iOUTnc, istat, idt1, kl
      REAL var(mx,my)   
      CHARACTER*(*) dbNam

C +...* Input:
      REAL tM_var(ndb, mz)

C +...* Indexes...
      INTEGER idb,ipe 


C +--Reset statistics
C +  ================

      IF (istat.EQ.0) THEN
        DO idb = 1,ndb
           tM_var(idb,kl)= 0.0
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
            tM_var(idb,kl)= tM_var(idb,kl)
     &        + var(ipe,idb) + var(ipe+1,my-idb+1)
         ENDDO
C +      * Left and Right contributions:
         DO ipe = idb, my-idb
            tM_var(idb,kl)= tM_var(idb,kl)
     &        + var(idb,ipe+1) + var(mx-idb+1,ipe)
         ENDDO
      
      ENDDO

C +--Final computation & write result (when istat = 2)
C +  =================================================
      IF (istat.EQ.2) THEN

C +   *Devide by the perimeter AND #time steps (idt1):
C
        DO idb = 1,ndb
             tM_var(idb,kl)= tM_var(idb,kl)
     &        / REAL ((mx+my-4*idb+2)*2*idt1)
        ENDDO

        IF (kl.EQ.mz) THEN
          CALL UNwrite (iOUTnc, dbNam,  0, ndb, mz , 1, tM_var)
        ENDIF
C +     ...all data written at the last call to this routine
      ENDIF

      RETURN
      END
