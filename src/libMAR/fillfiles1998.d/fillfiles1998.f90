PROGRAM origin

  IMPLICIT NONE

  INCLUDE "netcdf.inc"

  CHARACTER(LEN=100)     :: filename
  INTEGER, PARAMETER     :: nt = 100
  INTEGER, DIMENSION(nt) :: vector
  INTEGER                :: fidM, iret, varIDdate, varIDyear, annee, datelen, dimIDdate
  INTEGER                :: i, j, k, firstdate, ladate, nmois, nmoisp
  CHARACTER(2)           :: JJ(1:3)
  CHARACTER(2)           :: MM(1:12)
  DATA JJ/"01","12","23"/
  DATA MM/"01","02","03","04","05","06","07","08","09","10","11","12"/

  annee  = 1998 ! *** A MODIFIER A LA MAIN ***

  ! Boucle sur les fichiers
  DO i=1,12
     DO j=1,3
        write(*,*) " "
        write(*,*) "--- NEXT ---"
        write(*,*) " "
        filename = "/pc179/work500g/gential/pansementMARnc/MARANT.40km._BS/DSG_MAR.1998"//MM(i)//JJ(j)//".ANT.nc"
        write(*,*) filename

        ! The function NF_OPEN opens an existing netCDF dataset for access.
        ! Ouvrir MAR en �criture
        iret = NF_OPEN(TRIM(filename),NF_WRITE,fidM)
        CALL erreur(iret,.TRUE.,"open")

        ! The function NF_INQ_VARID returns the ID of a netCDF variable, given its
        ! name.
        iret = NF_INQ_VARID(fidM,"date",varIDdate)
        CALL erreur(iret,.TRUE.,"inq varid")

        iret = NF_INQ_VARID(fidM,"year",varIDyear)
        CALL erreur(iret,.TRUE.,"inq varid")

        ! Dimension temporelle
        iret = NF_INQ_DIMID(fidM,"time",dimIDdate)
        CALL erreur(iret,.TRUE.,"inq dimid")
        iret = NF_INQ_DIMLEN(fidM,dimIDdate,datelen)
        CALL erreur(iret,.TRUE.,"inquire dimension")

        ! The function NF_REDEF puts an open netCDF dataset into define mode, so
        ! dimensions, variables, and attributes can be added or renamed and attributes
        ! can be deleted.
        iret = NF_REDEF(fidM)
        CALL erreur(iret,.TRUE.,"redef")

        ! The function NF_ENDDEF takes an open netCDF dataset out of define mode. The
        ! changes made to the netCDF dataset while it was in define mode are checked
        ! and committed to disk if no problems occurred. Non-record variables may be
        ! initialized to a "fill value" as well (see NF_SET_FILL). The netCDF dataset
        ! is then placed in data mode, so variable data can be read or written.
        iret = NF_ENDDEF(fidM)
        CALL erreur(iret,.TRUE.,"enddef")

        ! Modification du champ year de MAR
        vector(:) = annee
        iret      = NF_PUT_VAR_INT(fidM,varIDyear,vector)
        CALL erreur(iret,.TRUE.,"put var int")

        ! Modification du champ date de MAR
        ! Rappel sur la date : 71306 signifie le 13 juillet à 6 h 
        iret   = NF_GET_VAR1_INT(fidM,varIDdate,1,firstdate)
        ladate = firstdate
        write(*,*) 'datelen   :', datelen
        write(*,*) 'firstdate :', firstdate
        DO k=2,datelen                          ! On ne modifie pas la 1re valeur
           ladate = ladate + 6                  ! Sorties toutes les 6 h
           IF (mod(ladate,100).EQ.24) THEN      ! On tombe sur 24 h
              ladate = ladate - 24 + 100        ! J à 24 h => J+1 à 0 h
              IF ((k.EQ.datelen).AND.(mod(ladate,10000).GE.2900)) THEN ! Fin de mois
                 nmois  = ladate / 10000        ! Quotient de la division euclidienne : mois qui s'achève
                 nmoisp = nmois + 1             ! Mois M+1
                 ladate = nmoisp * 10000 + 0100 ! Le 1er du mois M+1
                 write(*,*) "mois++"
              ENDIF
           ENDIF
           iret = NF_PUT_VAR1_INT(fidM,varIDdate,k,ladate)
           CALL erreur(iret,.TRUE.,"put var1 int")
        ENDDO

        ! The function NF_CLOSE closes an open netCDF dataset. If the dataset is in
        ! define mode, NF_ENDDEF will be called before closing. (In this case, if
        ! NF_ENDDEF returns an error, NF_ABORT will automatically be called to restore
        ! the dataset to the consistent state before define mode was last entered.)
        ! After an open netCDF dataset is closed, its netCDF ID may be reassigned to
        ! the next netCDF dataset that is opened or created.
        iret = NF_CLOSE(fidM)
        CALL erreur(iret,.TRUE.,"close")

     ENDDO
  ENDDO


END PROGRAM origin

SUBROUTINE erreur(iret, lstop, chaine)
  ! pour les messages d'erreur
  INTEGER, INTENT(in)                     :: iret
  LOGICAL, INTENT(in)                     :: lstop
  CHARACTER(LEN=*), INTENT(in)            :: chaine
  !
  CHARACTER(LEN=256)                      :: message
  !
  INCLUDE "netcdf.inc"
  !
  IF ( iret .NE. 0 ) THEN
    WRITE(*,*) 'ROUTINE: ', TRIM(chaine)
    WRITE(*,*) 'ERREUR: ', iret
    message=NF_STRERROR(iret)
    WRITE(*,*) 'CA VEUT DIRE:',TRIM(message)
    IF ( lstop ) STOP
  ENDIF
  !
END SUBROUTINE erreur
