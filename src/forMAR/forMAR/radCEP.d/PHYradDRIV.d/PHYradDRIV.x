#!/bin/sh

        EXE=`echo RRR`
# echo -n "RUN,GRAPH (RRR,GGG): "
# read  EXE

#       fco=`echo f90`
#       fco=`echo pgf90`
        fco=`echo ifc`
 if [ ${fco} = f90 ] || [ ${fco} = pgf90 ]; then
        OPT=`echo -C  -g`
 fi
 if [ ${fco} = ifc ]; then
        OPT=`echo -d2 -g`
 fi



if [ ${EXE} = GGG ];then
  echo " "
  echo "GRAPH"
  echo "====="
  echo " "

  rm   -f   ferret.jnl* *~ 
  ls   -lst      *.JNL
  echo -n "Chooses the graphic: "
  read  GRAPH


  echo "go "${GRAPH}".JNL" >& ${GRAPH}.jnl
  echo "exit"              >> ${GRAPH}.jnl
  ferret  < ${GRAPH}.jnl  
  Fprint -o ${GRAPH}.ps   -p portrait -l cps ${GRAPH}.plt  
  gv        ${GRAPH}.ps  
  mv        ${GRAPH}*.plt                    ${GRAPH}.plt  

  rm     -f ferret.jnl* *~
fi



if [ ${EXE} = RRR ];then
  echo " "
  echo "RUN"
  echo "==="
  echo " "

  rm   -f               *~ 
        PROGR=`echo PHYradDRIV`
# ls   -lst      *.f
# echo -n "Chooses the program: "
# read  PROGR

  cp -p  ../include/*.inc .
  cp -pf    ../MARdim.inc .
 if [ -f PHYrad2CEP.o ]; then
  echo  "PHYrad2CEP.o exists"
 else
  cp -p  ../../radCEP.dDB/* .
 fi

# ls  -1 *.DATA   >& ${PROGR}.LIST
  ${fco} ${OPT}   -o ${PROGR} ${PROGR}.f *.o
 if [ ${fco} = f90 ]; then
  totalview          ${PROGR}
 fi
 if [ ${fco} = pgf90 ] || [ ${fco} = ifc ]; then
                   ./${PROGR}
 fi
  rm  -f             ${PROGR}  work.pc*
fi



  echo " "
  echo "DATE and TIME"
  echo "============="
  echo " "

  datey=`date +%y`
  datem=`date +%m`
  dated=`date +%d`
  dateH=`date +%H`
  dateM=`date +%M`
  min_0=`echo 0`
  min_1=`echo 1`
  min_2=`echo 2`
  min_3=`echo 3`
  min_4=`echo 4`
  min_5=`echo 5`
  min10=`expr $dateM / 10`


# INSERT APROPRIATE TOUCH HERE, BEFORE DEFINING min20


if [ ${EXE} = GGG ]; then
  dateT=`echo $datem$dated$dateH$min10$min_0$datey`
  touch $dateT                            ${GRAPH}.JNL  
  dateT=`echo $datem$dated$dateH$min10$min_0$datey`
  touch $dateT                            ${GRAPH}.jnl  
  dateT=`echo $datem$dated$dateH$min10$min_5$datey`
  touch $dateT                            ${GRAPH}.ps .
fi


if [ ${EXE} = RRR ]; then
 if [ -f ${PROGR}.x ]; then
  dateT=`echo $datem$dated$dateH$min10$min_0$datey`
  touch $dateT                            ${PROGR}.x
 fi
  dateT=`echo $datem$dated$dateH$min10$min_0$datey`
  touch $dateT                            ${PROGR}.f   .
 if [ -f ${PROGR}.dat ]; then
  dateT=`echo $datem$dated$dateH$min10$min_5$datey`
  touch $dateT                            ${PROGR}.dat .
 fi
 if [ -f ${PROGR}.out ]; then
  dateT=`echo $datem$dated$dateH$min10$min_5$datey`
  touch $dateT                            ${PROGR}.out .
 fi
fi


  min20=`expr ${min10} + 1`
      if expr ${min20} = 6 
      then
                min20=0
                dateH=`expr ${dateH} + 1`
         echo ' dateH=     '${dateH}
      if expr ${dateH} \< 10 
      then
                dateH=`echo ${min_0}${dateH}`
      fi
      fi


if [ ${EXE} = RRR ] && [ -f ${PROGR}.LIST ]; then
  dateT=`echo $datem$dated$dateH$min20$min_0$datey`
  touch $dateT                            ${PROGR}.LIST .
fi
