#!/bin/sh
################################################################################
#
  datem=`date +%m`
  dated=`date +%d`
  dateH=`date +%H`
  dateM=`date +%M`
  datey=`date +%y`
  zero0=`echo 0`
  unun1=`echo 1`
  deux2=`echo 2`
  cinq5=`echo 5`
  six_6=`echo 6`
  neuf9=`echo 9`
#
################################################################################
#
f90 -o MAR_pp -O4 MAR_pp.f MAR_pp_dat.f
cp  -p MAR_pp     MAR_pp-Alpha
#
################# MAR_pp-Alpha #################################################
#
                 dateM=`echo $deux2$zero0`
                 dated=`expr $dated - 1`
  if       expr $dated \< 10
  then
                 dated=`echo $zero0$dated`
    if     expr $dated \<  1
    then
                 dated=`echo 30`
                 datem=`expr $datem - 1`
      if   expr $datem \< 10
      then
                 datem=`echo $zero0$datem`
        if expr $datem \<  1
        then
                 datem=`echo 12`
                 datey=`echo $datey - 1`
        fi
      fi
    fi
  fi
#
            dateT=`echo $datem$dated$dateH$dateM$datey`
#
     echo        'MAR_pp-Alpha'
     echo  $dateT
     touch $dateT MAR_pp-Alpha MAR_pp .
# 
################################################################################
#
################# MAR_pp-Alpha.x ###############################################
#
                 dateM=`echo $unun1$cinq5`
                 dated=`expr $dated - 1`
  if       expr $dated \< 10
  then
                 dated=`echo $zero0$dated`
    if     expr $dated \<  1
    then
                 dated=`echo 30`
                 datem=`expr $datem - 1`
      if   expr $datem \< 10
      then
                 datem=`echo $zero0$datem`
        if expr $datem \<  1
        then
                 datem=`echo 12`
                 datey=`echo $datey - 1`
        fi
      fi
    fi
  fi
#
            dateT=`echo $datem$dated$dateH$dateM$datey`
#
     echo        'MAR_pp-Alpha.x'
     echo  $dateT
     touch $dateT MAR_pp-Alpha.x MAR_pp-Alpha MAR_pp .
# 
################################################################################
#
################# MAR_pp.f #####################################################
#
                 dateM=`echo $unun1$zero0`
                 dated=`expr $dated - 1`
  if       expr $dated \< 10
  then
                 dated=`echo $zero0$dated`
    if     expr $dated \<  1
    then
                 dated=`echo 30`
                 datem=`expr $datem - 1`
      if   expr $datem \< 10
      then
                 datem=`echo $zero0$datem`
        if expr $datem \<  1
        then
                 datem=`echo 12`
                 datey=`echo $datey - 1`
        fi
      fi
    fi
  fi
#
            dateT=`echo $datem$dated$dateH$dateM$datey`
#
     echo        'MAR_pp.f'
     echo  $dateT
     touch $dateT MAR_pp.f
# 
################################################################################
#
################# MAR_pp_dat.f #################################################
#
                 dateM=`echo $zero0$cinq5`
                 dated=`expr $dated - 1`
  if       expr $dated \< 10
  then
                 dated=`echo $zero0$dated`
    if     expr $dated \<  1
    then
                 dated=`echo 30`
                 datem=`expr $datem - 1`
      if   expr $datem \< 10
      then
                 datem=`echo $zero0$datem`
        if expr $datem \<  1
        then
                 datem=`echo 12`
                 datey=`echo $datey - 1`
        fi
      fi
    fi
  fi
#
            dateT=`echo $datem$dated$dateH$dateM$datey`
#
     echo        'MAR_pp_dat.f'
     echo  $dateT
     touch $dateT MAR_pp_dat.f
# 
################################################################################
#
################# MAR_pp.inc ###################################################
#
                 dateM=`echo $zero0$cinq5`
                 dated=`expr $dated - 1`
  if       expr $dated \< 10
  then
                 dated=`echo $zero0$dated`
    if     expr $dated \<  1
    then
                 dated=`echo 30`
                 datem=`expr $datem - 1`
      if   expr $datem \< 10
      then
                 datem=`echo $zero0$datem`
        if expr $datem \<  1
        then
                 datem=`echo 12`
                 datey=`echo $datey - 1`
        fi
      fi
    fi
  fi
#
            dateT=`echo $datem$dated$dateH$dateM$datey`
#
     echo        'MAR_pp.inc'
     echo  $dateT
     touch $dateT MAR_pp.inc
# 
################################################################################
#
################# MAR_pp.inp ###################################################
#
                 dateM=`echo $zero0$zero0`
                 dated=`expr $dated - 1`
  if       expr $dated \< 10
  then
                 dated=`echo $zero0$dated`
    if     expr $dated \<  1
    then
                 dated=`echo 30`
                 datem=`expr $datem - 1`
      if   expr $datem \< 10
      then
                 datem=`echo $zero0$datem`
        if expr $datem \<  1
        then
                 datem=`echo 12`
                 datey=`echo $datey - 1`
        fi
      fi
    fi
  fi
#
            dateT=`echo $datem$dated$dateH$dateM$datey`
#
     echo        'MAR_pp.inp'
     echo  $dateT
     touch $dateT MAR_pp.inp.NEW
# 
################################################################################
