C     =========
C     libUN.inc
C     =========

      INCLUDE 'NetCDF.inc'
 
      CHARACTER*50   CF_institution 
      PARAMETER     (CF_institution="ULg (Xavier Fettweis)") 

      CHARACTER*10   CF_libUN_version 
      PARAMETER     (CF_libUN_version="2005.04.08") 

      INTEGER        CF_dimmaxlen
      PARAMETER     (CF_dimmaxlen=99999) ! Maximum dim/axes length

      INTEGER        CF_dimmaxnbr
      PARAMETER     (CF_dimmaxnbr=20)    ! Nbr Maximum of dim/axes                  


      INTEGER        CF_varmaxnbr       
      PARAMETER     (CF_varmaxnbr=300)   ! Nbr maximum of variables

      INTEGER        CF_attnbr       
      PARAMETER     (CF_attnbr   =1)     ! nbr of attibutes

      INTEGER        CF_dim(0:CF_dimmaxnbr),CF_attnum(CF_attnbr)
      INTEGER        CF_varnbrtot,CF_dimnbrtot,CF_fileidopened

      REAL           CF_dimval(CF_dimmaxlen,0:CF_dimmaxnbr)

      CHARACTER*13   CF_dimnam(0:CF_dimmaxnbr)
      CHARACTER*13   CF_varnam(CF_varmaxnbr)
      CHARACTER*13   CF_varnamdim(4,CF_varmaxnbr)
      CHARACTER*13   CF_attnam(CF_attnbr)
      CHARACTER*31   CF_dimnamuni(0:CF_dimmaxnbr)
      CHARACTER*31   CF_varnamuni(CF_varmaxnbr)
      CHARACTER*50   CF_vardes(CF_varmaxnbr)
      CHARACTER*200  CF_filenam,CF_filetit,CF_filenamopened

      COMMON/CF_FILE_CREATi/CF_dim,CF_attnum,CF_fileidopened, 
     .                      CF_varnbrtot,CF_dimnbrtot
                     
      COMMON/CF_FILE_CREATr/CF_dimval

      COMMON/CF_FILE_CREATc/CF_dimnam,CF_varnam,CF_varnamdim,
     .                      CF_attnam,CF_dimnamuni,CF_varnamuni,
     .                      CF_vardes,CF_filenam,CF_filetit,
     .                      CF_filenamopened

C     =================================================================

      INTEGER         iVarWarn,ireadchk
      REAL            VarRepl,vReadMin,vReadMax,
     .                vMissVal,ValRange(2)

      COMMON/UNparams/VarRepl,iVarWarn,ireadchk,vReadMin,vReadMax,
     .                vMissVal,ValRange

