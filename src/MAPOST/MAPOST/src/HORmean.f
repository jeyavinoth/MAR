      SUBROUTINE HORmean (InReg, varMAR, varVP)

C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
            
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'globals.inc'

      INTEGER ii, jj, nelem
      REAL varMAR(mx,my)
      REAL varVP(nreg)
      INTEGER InReg(mx,my,nreg)

      DO ireg=1,nreg
        varVP(ireg)= 0
        nelem      = 0

        DO jj= 1,my
        DO ii= 1,mx
           IF (varMAR(ii,jj) .LT. 1.E25) THEN
             varVP(ireg)=varVP(ireg)+varMAR(ii,jj)*InReg(ii,jj,ireg)
             nelem=nelem+InReg(ii,jj,ireg)
           ENDIF
        ENDDO
        ENDDO
        IF (nelem.GT.1) THEN
           varVP(ireg)= varVP(ireg)/nelem
        ENDIF
      ENDDO

      END
