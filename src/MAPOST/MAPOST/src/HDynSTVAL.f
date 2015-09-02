C +---------------------------------------------------------------------+
C | MAPOST                                                    03-1999 |
C |   SubRoutine HDynSTVAL                                              |
C |   Computes time average and write outputs for the "stat values"     |
C |            ! This routine is called by HDyn AND Surf !              |
C +---------------------------------------------------------------------+
C |     istat = Updating Level for statistics:                          |
C |                     (0 -> init, 1-> continue, 2->terminate work)    |   
C |                                                                     |
C +---------------------------------------------------------------------+
      SUBROUTINE HDynSTVAL (iOUTnc,istat,idt1,var,
     &               tM_nam,tM_var)


      IMPLICIT NONE
 
C +---LS and MAR domain dimensions :
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

C +---MAPOST parameters & variables
      include 'MAPOST.inc'
      include 'globals.inc'
      


C +...* Input:
      INTEGER iOUTnc, istat, idt1
      REAL var(nreg)
      CHARACTER*(*) tM_nam

C +...* Input/Output:
      REAL tM_var(nreg)

C +...* Indexes...
      INTEGER ireg     


C +--Reset statistics
C +  ================

      IF (istat.EQ.0) THEN

        DO ireg = 1,nreg
           tM_var(ireg)= 0.0
        END DO

      ENDIF 

C +--Update statistics for each valid data
C +  =====================================
      DO ireg = 1,nreg

         IF (var(ireg).LT.1.E20.AND.tM_var(ireg).LT.1.E20)THEN
             tM_var(ireg)= tM_var(ireg)+ var(ireg)
         ELSE
             tM_var(ireg)= 1.E30
C +          ..Missing Value
         ENDIF

      ENDDO

C +--Final computation & write result (when istat = 2)
C +  =================================================
      IF (istat.EQ.2) THEN

        DO ireg = 1,nreg
         IF (tM_var(ireg).LT.1.E25)THEN
           tM_var(ireg)= tM_var(ireg)/REAL(idt1)
         ENDIF
        ENDDO

        CALL UNwrite (iOUTnc, tM_nam,  0, nreg, 1, 1, tM_var)

      ENDIF

      RETURN
      END
