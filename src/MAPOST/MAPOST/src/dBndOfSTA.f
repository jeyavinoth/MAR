C +---------------------------------------------------------------------+
C | MAPOST                                                      02-2001 |
C |   SubRoutine dBndOfSTA                                              |
C |   Computes the "same dist to boundary" average, + time average,     |
C |   and writes result. This version modified => separe OutFlow points |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE dBndOfSTA (iOUTnc,istat, var, uLarge, vLarge,
     &          kl, idt1, dbNam,tM_var,tM_InF,tM_OutF,NP_InF,NP_OutF)


      IMPLICIT NONE
 
      INTEGER  VARSIZE
      EXTERNAL VARSIZE

C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
            
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'globals.inc'

C +...* Input:
      INTEGER iOUTnc, istat, idt1, kl
      REAL var(mx,my), uLarge(mx,my), vLarge(mx,my)
      CHARACTER*(*) dbNam
      CHARACTER*10, NamVar, NamInF, NamOutF

C +...* Input/output:
      REAL tM_var(ndb, mz), tM_InF(ndb, mz), tM_OutF(ndb, mz)
      INTEGER NP_InF(ndb,mz),NP_OutF(ndb, mz)

C +...* Indexes...
      INTEGER idb,ipe 


C +--Reset statistics
C +  ================

      IF (istat.EQ.0) THEN
        DO idb = 1,ndb
           tM_var (idb,kl)= 0.0
           tM_InF(idb,kl)= 0.0
           tM_OutF(idb,kl)= 0.0
           NP_InF (idb,kl)= 0
           NP_OutF(idb,kl)= 0
        END DO
      ENDIF 

C +--Update statistics.
C +  ==================
C +
C +
C +  (my)(2)---------->(mx-1)(my)    top
C +   +                       +
C +   +                       +
C +   +                       +
C +   +                       +
C +  (1) (2)---------->(mx-1)(1 )  bottom
C +
      DO idb = 1,ndb
         
C +      * Bottom and Top contributions:
         DO ipe = idb+1, mx-idb
            tM_var(idb,kl)= tM_var(idb,kl)
     &        + var(ipe,idb) + var(ipe,my-idb+1)

            IF (vLarge(ipe,3).LT.0.0)      THEN
              tM_OutF(idb,kl)=tM_OutF(idb,kl)+var(ipe,idb)
              NP_OutF(idb,kl)=NP_OutF(idb,kl)+ 1
            ELSE
              tM_InF(idb,kl)=tM_InF(idb,kl)+var(ipe,idb)
              NP_InF(idb,kl)=NP_InF(idb,kl)+ 1
            ENDIF

            IF (vLarge(ipe,my-2).GT.0.0) THEN
              tM_OutF(idb,kl)=tM_OutF(idb,kl)+var(ipe,my-idb+1)
              NP_OutF(idb,kl)=NP_OutF(idb,kl)+ 1
            ELSE
              tM_InF(idb,kl)=tM_InF(idb,kl)+var(ipe,my-idb+1)
              NP_InF(idb,kl)=NP_InF(idb,kl)+ 1
            ENDIF

         ENDDO

C +      * Left and Right contributions:
         DO ipe = idb, my-idb+1

            tM_var(idb,kl)= tM_var(idb,kl)
     &        + var(idb,ipe) + var(mx-idb+1,ipe)

            If (uLarge(3,ipe).LT.0.0)      THEN 
              tM_OutF(idb,kl)=tM_OutF(idb,kl)+var(idb,ipe)
              NP_OutF(idb,kl)=NP_OutF(idb,kl)+ 1
            ELSE
              tM_InF(idb,kl)=tM_InF(idb,kl)+var(idb,ipe)
              NP_InF(idb,kl)=NP_InF(idb,kl)+ 1
            ENDIF

            If (uLarge(mx-2,ipe).GT.0.0) THEN
              tM_OutF(idb,kl)=tM_OutF(idb,kl)+var(mx-idb+1,ipe)
              NP_OutF(idb,kl)=NP_OutF(idb,kl)+ 1
            ELSE
              tM_InF(idb,kl)=tM_InF(idb,kl)+var(mx-idb+1,ipe)
              NP_InF(idb,kl)=NP_InF(idb,kl)+ 1
            ENDIF

         ENDDO
      
      ENDDO

C +--Final computation & write result (when istat = 2)
C +  =================================================
      IF (istat.EQ.2) THEN

        DO idb = 1,ndb
          tM_var(idb,kl)=tM_var(idb,kl)/REAL((mx+my-4*idb+2)*2*idt1)
          tM_InF (idb,kl)= tM_InF (idb,kl)/REAL (NP_InF(idb,kl))
          tM_OutF(idb,kl)= tM_OutF(idb,kl)/REAL(NP_OutF(idb,kl))
        ENDDO

        NamVar = dbNam(1:VARSIZE(dbNam))//'_db'
        NamInF = dbNam(1:VARSIZE(dbNam))//'idb'
        NamOutF= dbNam(1:VARSIZE(dbNam))//'odb'
        IF (kl.EQ.mz) THEN
          CALL UNwrite (iOUTnc, NamVar, 0,ndb,mz,1, tM_var)
          CALL UNwrite (iOUTnc, NamInF ,0,ndb,mz,1, tM_InF )
          CALL UNwrite (iOUTnc, NamOutF,0,ndb,mz,1, tM_OutF)
        ENDIF
C +     ...all data written at the last call to this routine
      ENDIF

      RETURN
      END
