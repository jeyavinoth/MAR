C +---Prefix for all directories (e.g. your work directory)
C +   -----------------------------------------------------
 
      PFXdir = './input/'

      PFXsiz = VARSIZE(PFXdir)


C +---Sea Surface Temperature (Reynolds - 1 deg)
C +   ------------------------------------------

      SST_dir = PFXdir(1:PFXsiz)//'SST/'


C +---ETOPO topography - 5 minutes (Netcdf file)
C +   ------------------------------------------

      ETOPO_dir = PFXdir(1:PFXsiz)//'ETOPO/'


C +---GTOPO topography - 30 secondes (Netcdf files)
C +   ---------------------------------------------

      GTOPO_dir = PFXdir(1:PFXsiz)//'GTOPO/'


C +---BTOPO topography - 30-m resol. over Belgium
C +   -------------------------------------------

      BTOPO_dir = PFXdir(1:PFXsiz)//'BTOPO/'


C +---FAO SOIL TYPES (Ascii and Netcdf files)
C +   ---------------------------------------

      FAO_dir = PFXdir(1:PFXsiz)//'FAO/'


C +---VEGETATION - IGBP Classification
C +   --------------------------------

C +...AFRICA
      AFRveg_dir = PFXdir(1:PFXsiz)//'VEGE/'
C +...EUROPE
      EURveg_dir = PFXdir(1:PFXsiz)//'VEGE/'
C +...BELGIUM
      BELveg_dir = PFXdir(1:PFXsiz)//'VEGE/'
C +...NORTH AMERICA
      NAMveg_dir = PFXdir(1:PFXsiz)//'VEGE/'
C +...SOUTH AMERICA
      SAMveg_dir = PFXdir(1:PFXsiz)//'VEGE/'
C +...United States
      USveg_dir = PFXdir(1:PFXsiz)//'VEGE/'

C +---VEGETATION - CORINE Classification - EUROPE
C +   -------------------------------------------

      CORveg_dir = PFXdir(1:PFXsiz)//'CORINE/'


C +---VEGETATION - 1-km NDVI index (April 1992 to March 1993)
C +   -------------------------------------------------------

C +...AFRICA
      AFRndv_dir = PFXdir(1:PFXsiz)//'NDVI01/'
C +...EUROPE
      EURndv_dir = PFXdir(1:PFXsiz)//'NDVI01/'
C +...NORTH AMERICA
      NAMndv_dir = PFXdir(1:PFXsiz)//'NDVI01/'
C +...SOUTH AMERICA
      SAMndv_dir = PFXdir(1:PFXsiz)//'NDVI01/'


C +---VEGETATION - 8-km NDVI index (??? to ???)
C +   -----------------------------------------

C +...AFRICA
      AFRndv8dir = PFXdir(1:PFXsiz)//'NDVI08/'
C +...EUROPE
      EURndv8dir = PFXdir(1:PFXsiz)//'NDVI08/'
C +...NORTH AMERICA
      NAMndv8dir = PFXdir(1:PFXsiz)//'NDVI08/'
C +...SOUTH AMERICA
      SAMndv8dir = PFXdir(1:PFXsiz)//'NDVI08/'


