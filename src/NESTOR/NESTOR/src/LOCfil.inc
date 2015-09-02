C +---Files locations
C +   ===============

      CHARACTER*60 ETOPO_dir ,GTOPO_dir ,FAO_dir   ,BELveg_dir,
     .             AFRveg_dir,EURveg_dir,NAMveg_dir,SAMveg_dir,
     .             AFRndv_dir,EURndv_dir,NAMndv_dir,SAMndv_dir,
     .             AFRndv8dir,EURndv8dir,NAMndv8dir,SAMndv8dir,
     .             CORveg_dir,SST_dir,BTOPO_dir,PFXdir

      INTEGER PFXsiz

      common /LOCfil/
     .             ETOPO_dir ,GTOPO_dir ,FAO_dir   ,BELveg_dir,
     .             AFRveg_dir,EURveg_dir,NAMveg_dir,SAMveg_dir,
     .             AFRndv_dir,EURndv_dir,NAMndv_dir,SAMndv_dir,
     .             AFRndv8dir,EURndv8dir,NAMndv8dir,SAMndv8dir,
     .             CORveg_dir,SST_dir,BTOPO_dir,PFXdir,PFXsiz
