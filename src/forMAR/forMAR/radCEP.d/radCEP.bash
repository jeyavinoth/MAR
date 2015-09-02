#!/bin/bash

# set echo
  rm -f      *~


  echo " "
  echo "USAGE: radCEP.sh ${1} ${2} ${3} ${4}, with: Vector   Length                                               : "${1}
  echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        Nb  of   Levels                                               : "${2}

  echo "PWD: `pwd`"

if expr ${2} \> 600 
then
   echo "NO COMPILATION: Please increase the dimension of CETA    in module_1/yoecld.F90"
   echo "                                                 CVDAES  in module_1/yoeaerd.F90"
   echo "                                                 CVDAEL  in module_1/yoeaerd.F90"
   echo "                                                 CVDAEU  in module_1/yoeaerd.F90"
   echo "                                                 CVDAED  in module_1/yoeaerd.F90"
   echo "                                                 RA1OVLP in module_1/yoeovlp.F90"
else

   CODE=`echo forMAR`
   SPLT=`echo SPLIT`

if [ ${SPLT} = SPLIT ]; then
   echo " "
   echo "** Splited Routines will be compiled **"
   echo " "
else
   cd Source.d
 if [ -d su____.d.NOsplit ];then
      mv su____.d         su____.d.__SPLIT
      mv su____.d.NOsplit su____.d
 fi
 if [ -d rrtm__.d.NOsplit ];then
      mv rrtm__.d         rrtm__.d.__SPLIT
      mv rrtm__.d.NOsplit rrtm__.d
 fi
   cd ../.
fi


                        COMPIFUL="OK"
                         COMPILE="OK"

                         OPTIONS="$opt $openmp" # Optimisation (avec openmp)
                         OPTION2="$opt $openmp" # Optimisation (avec openmp)
    [ $1 -gt 1      ] && OPTIONS="$opt"         # sans openmp

    [ ${#op2} -gt 1 ] && OPTION2="$op2 $openmp" # Optimisation (sans vectorisation)
                        COMPILER=$foc 

#     OPTIONS="-w -zero -vec_report0 -static -mp1 -ip -O3 -xSSE4.1 -traceback"
#     OPTION2="-w -zero -vec_report0 -static -mp1 -ip -O3 -xSSE4.1 -traceback"

 [ ${#COMPILER} -eq 0 ] && echo "ERROR: COMPILER" && exit

    echo " =========================="
    echo " Compiler OPTIONS (general): "${OPTIONS}
    echo " Compiler OPTIONS (RRTM)   : "${OPTION2}
    echo " =========================="
    echo " "


  echo " "
  echo "Dimension: Set-UP"
  echo "^^^^^^^^^^^^^^^^^"
  echo " "

  rm -f                                                                              radCEP.inc
  echo '      integer   klonr       ,klevr   ,nn_aer,nae'                         >> radCEP.inc
  echo '      parameter(klonr='${1}',klevr='${2}',nn_aer=6)'                      >> radCEP.inc

  rm -f                                                                              parrtm1d.F90
  echo 'MODULE PARRTM1D'                                                          >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo '#include "tsmbkind.h"'                                                    >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo 'IMPLICIT NONE'                                                            >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo 'SAVE'                                                                     >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo '!     ------------------------------------------------------------------' >> parrtm1d.F90
  echo '!     Parameters for 1-D radiation only computations from operational'    >> parrtm1d.F90
  echo '!      library routines'                                                  >> parrtm1d.F90
  echo '!     991007  JJMorcrette'                                                >> parrtm1d.F90
  echo '!     ------------------------------------------------------------------' >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_LON  = '${1}                                   >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_IDIA = 1'                                      >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_FDIA = JP_LON-JP_IDIA+1'                       >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_TDIA = 1'                                      >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo '!-- standard tropical INTEGER_M, PARAMETER :: JP_LEV  = 64'               >> parrtm1d.F90
  echo '!-- ATEX              INTEGER_M, PARAMETER :: JP_LEV  = 83'               >> parrtm1d.F90
  echo '!-- BOMEX             INTEGER_M, PARAMETER :: JP_LEV  = 84'               >> parrtm1d.F90
  echo '!-- OPEN_CELLS        INTEGER_M, PARAMETER :: JP_LEV  = 63'               >> parrtm1d.F90
  echo '!-- GATE_A,B,C        INTEGER_M, PARAMETER :: JP_LEV  = 46'               >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_LEV  = '${2}                                   >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_LW   = 6'                                      >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_SW   = 6'                                      >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_NUA  = 24'                                     >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_MODE = 1'                                      >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_AER  = 6'                                      >> parrtm1d.F90
  echo 'INTEGER_M, PARAMETER :: JP_LEVP1= JP_LEV+1'                               >> parrtm1d.F90
  echo ''                                                                         >> parrtm1d.F90
  echo '!     ------------------------------------------------------------------' >> parrtm1d.F90
  echo 'END MODULE PARRTM1D'                                                      >> parrtm1d.F90


  rm -f                                                                              parrrtm.F90
  echo 'MODULE PARRRTM'                                                           >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo '#include "tsmbkind.h"'                                                    >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo 'IMPLICIT NONE'                                                            >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo 'SAVE'                                                                     >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo '!     ------------------------------------------------------------------' >> parrrtm.F90
  echo '!     Parameters relevant to AER s RRTM-LW radiation scheme'              >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo '!     980714  JJMorcrette'                                                >> parrrtm.F90
  echo '!     ------------------------------------------------------------------' >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo '!-- standard tropical INTEGER_M, PARAMETER :: JPLAY  = 64'                >> parrrtm.F90
  echo '!-- ATEX              INTEGER_M, PARAMETER :: JPLAY  = 83'                >> parrrtm.F90
  echo '!-- BOMEX             INTEGER_M, PARAMETER :: JPLAY  = 84'                >> parrrtm.F90
  echo '!-- OPEN_CELLS        INTEGER_M, PARAMETER :: JPLAY  = 63'                >> parrrtm.F90
  echo '!__ GATE_A,B,C        INTEGER_M, PARAMETER :: JPLAY  = 46'                >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: JPLAY  = '${2}                                    >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: JPG    = 16'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: JPBAND = 16'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: JPXSEC = 4'                                       >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: JPINPX = 35'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: JPGPT  = 140'                                     >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
if [ ${CODE} = ORIGINAL ]; then
  echo 'INTEGER_M, PARAMETER :: NG1  = 8'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG2  = 14'                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG3  = 16'                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG4  = 14'                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG5  = 16'                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG6  = 8'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG7  = 12'                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG8  = 8'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG9  = 12'                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG10 = 6'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG11 = 8'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG12 = 8'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG13 = 4'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG14 = 2'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG15 = 2'                                         >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NG16 = 2'                                         >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
fi
  echo 'INTEGER_M, PARAMETER :: NGS1  = 8'                                        >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS2  = 22'                                       >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS3  = 38'                                       >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS4  = 52'                                       >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS5  = 68'                                       >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS6  = 76'                                       >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS7  = 88'                                       >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS8  = 96'                                       >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS9  = 108'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS10 = 114'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS11 = 122'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS12 = 130'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS13 = 134'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS14 = 136'                                      >> parrrtm.F90
  echo 'INTEGER_M, PARAMETER :: NGS15 = 138'                                      >> parrrtm.F90
  echo ' '                                                                        >> parrrtm.F90
  echo '!     ------------------------------------------------------------------' >> parrrtm.F90
  echo 'END MODULE PARRRTM'                                                       >> parrrtm.F90


  cp  -fp                                                                            parrtm1d.F90 module.d_1/.
  cp  -fp                                                                            parrrtm.F90  module.d_1/.

  rm -rf OBJECT.d/ SOURCE.d/ WorkArea/ &>/dev/null

  echo " "
  echo "START radCEP COMPILATION"
  echo "^^^^^^^^^^^^^^^^^^^^^^^^"
  echo " "

if [ ${COMPILE} == "OK" ]; then

  echo "OPTIONS: "${OPTIONS}


 if [ ${COMPIFUL} == "OK" ]; then
  mkdir              SOURCE.d
  mkdir              OBJECT.d

  echo " "
  echo "AUXILIARY MODULES"                                             # L'ordre des compilations est important
  echo "^^^^^^^^^^^^^^^^^"
  echo " "
  mkdir                                WorkArea
  cp                   -p include/*.h  WorkArea/.
  cp                   -p include/*.h  SOURCE.d/.
  cp                   -p module.d_0/* WorkArea/.
  cp                   -p module.d_0/* SOURCE.d/.
  cp                   -p radCEP.inc   SOURCE.d/.
  cd                                   WorkArea
  for file in *.F90 ; do
   echo "${COMPILER} ${OPTIONS} -c $file"
         ${COMPILER} ${OPTIONS} -c $file
  done
  cp                   -p *.o       ../OBJECT.d/.
  cp                   -p *.mod     ../OBJECT.d/.
  cd                                ../.
  rm                  -rf              WorkArea

  

  echo " "
  echo "MAIN MODULES"                                                  # L'ordre des compilations est important
  echo "^^^^^^^^^^^^"
  echo " "
  mkdir                                WorkArea
  cp                   -p include/*.h  WorkArea/.
  cp                   -p OBJECT.d/*   WorkArea/.
  cp                   -p module.d_1/* WorkArea/.
  cp                   -p module.d_1/* SOURCE.d/.
  cd                                   WorkArea
  for file in *.F90 ; do
   echo "${COMPILER} ${OPTIONS} -c $file"
         ${COMPILER} ${OPTIONS} -c $file
  done
  cp                   -p *.o       ../OBJECT.d/.
  cp                   -p *.mod     ../OBJECT.d/.
  cd                                ../.
  rm                  -rf              WorkArea

  

  echo " "
  echo "SOURCES"                                                       # L'ordre des compilations est important
  echo "^^^^^^^"                      
  echo " "
  echo    "S0    Scheme: If needed tape a blank character in order to continue"
  echo    "~~~~~~~~~~~~"
# read  blank
  mkdir                                       WorkArea
  cp                   -p include/*.h         WorkArea/.
  cp                   -p OBJECT.d/*          WorkArea/.
  cp                   -p Source.d/sw____.d/* WorkArea/.
  cp                   -p Source.d/sw____.d/* SOURCE.d/.
  cd                                          WorkArea
  for file in sw*.F90 ; do
   echo "${COMPILER} ${OPTIONS} -c $file"
         ${COMPILER} ${OPTIONS} -c $file
  done
  cp                   -p *.o              ../OBJECT.d/.
  cp                   -p *.mod            ../OBJECT.d/.
  cd                                       ../.
  rm                  -rf                     WorkArea


  echo " "
  echo    "IR    Scheme: If needed tape a blank character in order to continue"
  echo    "~~~~~~~~~~~~"
# read  blank
  mkdir                                       WorkArea
  cp                   -p include/*.h         WorkArea/.
  cp                   -p OBJECT.d/*          WorkArea/.
  cp                   -p Source.d/olw___.d/* WorkArea/.
  cp                   -p Source.d/olw___.d/* SOURCE.d/.
  cd                                          WorkArea
  for file in olw*.F90 ; do
   echo "${COMPILER} ${OPTIONS} -c $file"
         ${COMPILER} ${OPTIONS} -c $file
  done
  cp                   -p *.o              ../OBJECT.d/.
  cp                   -p *.mod            ../OBJECT.d/.
  cd                                       ../.
  rm                  -rf                     WorkArea


  echo " "
  echo    "IR n  Scheme: If needed tape a blank character in order to continue"
  echo    "~~~~~~~~~~~~"
# read  blank
  mkdir                                       WorkArea
  cp                   -p include/*.h         WorkArea/.
  cp                   -p OBJECT.d/*          WorkArea/.
  cp                   -p Source.d/lw____.d/* WorkArea/.
  cp                   -p Source.d/lw____.d/* SOURCE.d/.
  cd                                          WorkArea
  for file in lw*.F90 ; do
   echo "${COMPILER} ${OPTIONS} -c $file"
         ${COMPILER} ${OPTIONS} -c $file
  done
  cp                   -p *.o              ../OBJECT.d/.
  cp                   -p *.mod            ../OBJECT.d/.
  cd                                       ../.
  rm                  -rf                     WorkArea


  echo " "
  echo    "S0/IR DATA:   If needed tape a blank character in order to continue"
  echo    "~~~~~~~~~~"
# read  blank
  mkdir                                       WorkArea
  cp                   -p include/*.h         WorkArea/.
  cp                   -p OBJECT.d/*          WorkArea/.
  cp                   -p Source.d/su____.d/* WorkArea/.
  cp                   -p Source.d/su____.d/* SOURCE.d/.
  cd                                          WorkArea
  for file in su*.F90 ; do
   echo "${COMPILER} ${OPTIONS} -c $file"
         ${COMPILER} ${OPTIONS} -c $file
  done
  cp                   -p *.o              ../OBJECT.d/.
  cp                   -p *.mod            ../OBJECT.d/.
  cd                                       ../.
  rm                  -rf                     WorkArea


  echo " "
  echo    "RRTM  Scheme: If needed tape a blank character in order to continue"
  echo    "~~~~~~~~~~~~"
# read  blank
  mkdir                                       WorkArea
  cp                   -p include/*.h         WorkArea/.
  cp                   -p OBJECT.d/*          WorkArea/.
  cp                   -p Source.d/rrtm__.d/* WorkArea/.
  cp                   -p Source.d/rrtm__.d/* SOURCE.d/.
  cd                                          WorkArea
  for file in rr*.F90 ; do
   echo "${COMPILER} ${OPTION2} -c $file"
         ${COMPILER} ${OPTION2} -c $file
  done
  cp                     -p *.o              ../OBJECT.d/.
  cp                     -p *.mod            ../OBJECT.d/.
  cd                                         ../.
  rm                    -rf                     WorkArea


  echo " "
  echo    "Miscellanea:  If needed tape a blank character in order to continue"
  echo    "~~~~~~~~~~~"
# read  blank

  mkdir                                       WorkArea
  cp                   -p include/*.h         WorkArea/.
  cp                   -p OBJECT.d/*          WorkArea/.
  cp                   -p Source.d/Divers.d/* WorkArea/.
  cp                   -p Source.d/Divers.d/* SOURCE.d/.
  cd                                          WorkArea
  for file in *.F90 ; do
   echo "${COMPILER} ${OPTIONS} -c $file"
         ${COMPILER} ${OPTIONS} -c $file
  done
  cp                   -p *.o              ../OBJECT.d/.
  cp                   -p *.mod            ../OBJECT.d/.
  cd                                       ../.
  rm                  -rf                     WorkArea
 fi

  echo " "
  echo    "PHYrad2CEP:   If needed tape a blank character in order to continue"
  echo    "~~~~~~~~~~"
# read  blank
  mkdir                                           WorkArea
  cp                   -p include/*               WorkArea/.
  cp                   -p radCEP.inc              WorkArea/.
  cp                   -p radCEP.inc              SOURCE.d/.
  cp                   -p OBJECT.d/*              WorkArea/.
  cp                   -p Source.d/PHYrad2CEP.F90 WorkArea/.
  cp                   -p Source.d/PHYrad2CEP.F90 SOURCE.d/.
  cd                                              WorkArea
  echo "${COMPILER} ${OPTIONS} -c          PHYrad2CEP.F90"
        ${COMPILER} ${OPTIONS} -c          PHYrad2CEP.F90                           # radCEP / MAR interface routine
  cp                   -p *.inc                ../OBJECT.d/.
  cp                   -p *.h                  ../OBJECT.d/.
  cp                   -p *.o                  ../OBJECT.d/.
  cp                   -p *.mod                ../OBJECT.d/.
  cd                                           ../.
  rm                  -rf                         WorkArea


  if [ -d ../radCEP_${2}.d${COMPILER} ]; then
   rm -fr ../radCEP_${2}.d${COMPILER}
  fi
  cp -p radCEP.inc    OBJECT.d/.
  cp -p radCEP.inc    SOURCE.d/.
  cp -p parrtm1d.F90  OBJECT.d/.
  cp -p parrtm1d.F90  SOURCE.d/.
  cp -p parrrtm.F90   OBJECT.d/.
  cp -p parrrtm.F90   SOURCE.d/.
 if [ -f radCEP.log ]; then
  cp -p radCEP.log    OBJECT.d/.
 fi
  mv                  OBJECT.d                radCEP_${2}.d${COMPILER}
  echo " "
  echo "END od radCEP COMPILATION"
  echo "^^^^^^^^^^^^^^^^^^^^^^^^^"
  echo " "

else
  echo " "
  mkdir                                           SOURCE.d
  cp                   -p include/*.h             SOURCE.d/.
  cp                   -p module.d_0/*            SOURCE.d/.
  cp                   -p module.d_1/*            SOURCE.d/.
  cp                   -p radCEP.inc              SOURCE.d/.
  cp                   -p Source.d/sw____.d/*     SOURCE.d/.
  cp                   -p Source.d/olw___.d/*     SOURCE.d/.
  cp                   -p Source.d/lw____.d/*     SOURCE.d/.
  cp                   -p Source.d/su____.d/*     SOURCE.d/.
  cp                   -p Source.d/rrtm__.d/*     SOURCE.d/.
  cp                   -p Source.d/Divers.d/*     SOURCE.d/.
  cp                   -p Source.d/PHYrad2CEP.F90 SOURCE.d/.
# echo "Bad Compilation Options ("${4}") ==> the COMPILATION of ithe ECMWF Radiation package is not performed"
fi

fi
