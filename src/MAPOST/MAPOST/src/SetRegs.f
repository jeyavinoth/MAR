C +----------------------------------------------------------------------+
C | MAR post-processing                                       3.1.0      |
C |                                                                      |
C | SetRegs                                                              |
C |  Defines the averaging sub-regions.                                  |
C |  (InReg(ii,jj,ireg)=1 if point (ii,jj) belongs to reg # ireg)        |
C | WARNING: InReg is an INTEGER                                         |
C |                                                                      |
C |  ireg = 10-12=> internally defined                                   |
C |                 11 : isol = 1,2,3,4                                  |
C |                 12 : isol = 4                                        |
C |                 13 : isol = 1,2                                      |
C |                 14 : isol = 3                                        |
C |                                                                      |
C |  Under development:                                                  |
C |     definition of surface observation stations                       |
C |    wmSta (mx,my,nsta)  ! Weight of MAR value -> this station         |
C | WARNING: wmSTA is a REAL                                             |
C |                                                                      |
C +----------------------------------------------------------------------+
C
C +
      SUBROUTINE SetRegs (MARisol,shMAR,MARlon,MARlat,InReg,RegMsk,
     .                    wmSTA,StaMsk)
      
      IMPLICIT NONE
      
C +---LS and MAR domain dimensions :
C +   -----------------------------
      INCLUDE 'NSTdim.inc'
      INCLUDE 'NSTtoMAP.inc'

      include 'globals.inc'
      include 'MAPOST.inc'
      
      REAL MARisol(mx,my), shMAR(mx,my),MARlon(mx,my), MARlat(mx,my)
      REAL aLat, aLon, bLat, bLon
      INTEGER InReg(mx,my,nreg), isol(mx,my)
      REAL RegMsk(mx,my)
      INTEGER ii, jj, iexcl, iSpecReg, ista, ireg
      CHARACTER*21 INPstr
      CHARACTER*80 RegDFile
      INTEGER isolCur      ! Currently allowed isol
      REAL    shMin, shMax ! Currently allowed min,max sh

C +   For the stations:
      REAL wmSta (mx,my,nsta)  ! Weight of MAR value -> this station
      REAL StaMsk(mx,my)
      
      iexcl   = 7 
C +...Number of excluded near-boundary points.

      isolCur = 0          ! All isol types allowed
      shMin   = -1000.
      shMax   = 10000.
C +...Default characteristics of allowed points

      iSpecReg= 11 - 1
C +...Sets the number of the first automatic. def. region 

      DO ireg=1,nreg
        DO jj= 1,my
        DO ii= 1,mx
          InReg(ii,jj,ireg)= 0
        ENDDO
        ENDDO
      ENDDO

      DO jj= 1,my
      DO ii= 1,mx
         RegMsk(ii,jj)= 0
         isol  (ii,jj)= MARisol(ii,jj)
      ENDDO
      ENDDO

      ireg= iSpecReg + 1
C +...All points except excluded boundary
        DO jj= iexcl,my-iexcl
        DO ii= iexcl,mx-iexcl
          InReg(ii,jj,ireg)= 1
        ENDDO
        ENDDO

      ireg= iSpecReg + 2
C +...All LAND or TUNDRA points except excluded boundary
        DO jj= iexcl,my-iexcl
        DO ii= iexcl,mx-iexcl
          IF (isol(ii,jj).GT.3) THEN
            InReg(ii,jj,ireg)= 1
            write(*,*) 'SetRegs',ii,jj
          ENDIF
        ENDDO
        ENDDO

      ireg= iSpecReg + 3
C +...All OCEAN/SEA-ICE points except excluded boundary
        DO jj= iexcl,my-iexcl
        DO ii= iexcl,mx-iexcl
          IF (isol(ii,jj).LT.3) THEN
            InReg(ii,jj,ireg)= 1
          ENDIF
        ENDDO
        ENDDO

      ireg= iSpecReg + 4
C +...All ICE SHEET points except excluded boundary
        DO jj= iexcl,my-iexcl
        DO ii= iexcl,mx-iexcl
          IF (isol(ii,jj).EQ.3) THEN
            InReg(ii,jj,ireg)= 1
          ENDIF
        ENDDO
        ENDDO

C +--Read sub-region definitions from file.
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C +
        RegDFile = 'RegDefs/'//REGstudy//'.RegD'
        OPEN(unit=51,status='unknown',file=RegDFile)
C +
        ireg= 0    
 980    CONTINUE
          READ (51,'(A21)',end=1010) INPstr
          IF (INPstr(1:1).EQ.'>') THEN
          
            IF (INPstr(2:7).EQ.'REGNUM') THEN
            
              READ (INPstr,'(8x,I4)') ireg
              WRITE(*,*) 'Reading def. for subregion ', ireg
              IF (ireg.GT.nreg) THEN 
                 STOP
              ENDIF
              READ (51,'(A21)') INPstr
              READ (INPstr,*) aLat, aLon
              READ (51,'(A21)') INPstr
              READ (INPstr,*) bLat, bLon
              DO jj= iexcl,my-iexcl
              DO ii= iexcl,mx-iexcl
                IF  (aLat.LE.MARlat(ii,jj)
     &          .AND.bLat.GE.MARlat(ii,jj)
     &          .AND.aLon.LE.MARlon(ii,jj)
     &          .AND.bLon.GE.MARlon(ii,jj) 
     &          .AND.shMIN.LE.shMAR(ii,jj) 
     &          .AND.shMax.GE.shMAR(ii,jj)   
     &          .AND.(isolCUR.EQ.MARisol(ii,jj)
     &                .OR.isolCUR.EQ.0)
     &          ) THEN   
                  InReg(ii,jj,ireg)= 1
                  IF (isolCur.EQ.0) THEN
                    RegMsk(ii,jj)    = ireg
                  ENDIF
                ENDIF
              ENDDO
              ENDDO
              
            ELSE IF (INPstr(2:8).EQ.'STATION') THEN
            
              READ (INPstr,'(9x,I4)') ireg
              WRITE(*,*) 'Reading def. for station ', ireg 
              IF (ireg.GT.nreg) THEN 
                 STOP
              ENDIF
              READ (51,'(A21)') INPstr
              READ (INPstr,*) ii, jj ! Coordinates MAR
              InReg(ii,jj,ireg)= 1
              RegMsk(ii,jj)    = ireg

            ELSE IF (INPstr(2:11).EQ.'WM-STATION') THEN

C             Proposition for observation stations
C              corresponding to a weighted average of MAR points
              READ (INPstr,'(9x,I4)') ista
              WRITE(*,*) 'Reading def. for wm-station ', ista
              IF (ista.GT.nsta) THEN
                 STOP
              ENDIF
              READ (51,'(A21)') INPstr
              READ (INPstr,*) ii, jj ! Coordinates MAR
C                 on pourrait inclure une recherche du point
C                 le plus proche en lat/lon et eventuellement
C                 l'usage de plusieurs pts MAR avec poids/distance
              wmSta (ii,jj,ista)= 1.0
              StaMsk(ii,jj)     = ista

            ELSE IF (INPstr(2:7).EQ.'SOLTYP') THEN
              READ (INPstr,'(8x,I4)') isolCur 

            ELSE IF (INPstr(2:8).EQ.'ALTIMIN') THEN
              READ (INPstr,'(9x,F8.1)') shMin

            ELSE IF (INPstr(2:8).EQ.'ALTIMAX') THEN
              READ (INPstr,'(9x,F8.1)') shMax

            ELSE
              WRITE(*,*) 'WARNING: '
              WRITE(*,*) '  unknown line in your .RegD file:' 
              WRITE(*,*) INPstr
              WRITE(*,*)        
            ENDIF
            
          ENDIF
        GOTO 980
      
 1010   CONTINUE
        CLOSE(unit=51) 
      END
