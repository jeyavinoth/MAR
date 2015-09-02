#!/bin/ksh
# 
#  v1.2 - 20.01.2006 (16.09.2005)
#
# NAME
#
#   POMME.x - Executes NCO commands to achieve the regional means on POMME
#             output files (in batch on Rhodes)
#                                                                               
# SYNOPSIS                                                                      
#                                                                               
#   <nihil>
#                                                                               
# OPERANDS                                                                      
#                                                                               
#   <nihil>
#                                                                               
# DESCRIPTION                                                                   
#                                                                               
#   <nihil>
#

echo `date`
#----------------------------------------
# NCO: regional mean
#----------------------------------------

for region in "WAM" "WSA" "CSA" "ESA" "CSU" "CGU" "NIA" "OUE" "GOU" ; do

  #define region borders

  case $region in
  (WAM) imin=19 && imax=117 && jmin=35 && jmax=75 && tmp="WAM domain" ;;
  (WSA) imin=19 && imax=40  && jmin=58 && jmax=72 && tmp="Western Sahel" ;;
  (CSA) imin=41 && imax=95  && jmin=58 && jmax=72 && tmp="Central Sahel" ;;
  (ESA) imin=96 && imax=117 && jmin=58 && jmax=72 && tmp="Eastern Sahel" ;;
  (CSU) imin=41 && imax=95  && jmin=44 && jmax=57 && tmp="Central Sudan" ;;
  (CGU) imin=41 && imax=95  && jmin=35 && jmax=43 && tmp="Central Guinea" ;;
  (NIA) imin=72 && imax=76  && jmin=57 && jmax=60 && tmp="Niamey area" ;;
  (OUE) imin=72 && imax=76  && jmin=46 && jmax=49 && tmp="Oueme catchment" ;;
  (GOU) imin=63 && imax=65  && jmin=61 && jmax=70 && tmp="Gourma" ;;
  (*)   DAMNED "unknown region" && exit ;;
  esac

  # +------+------+------+
  # | i1j3 | i2j3 | i3j3 |
  # +------+------+------+
  # | i1j2 | i2j2 | i3j2 |  WAM region: to be cut for reasons of memory
  # +------+------+------+  restriction (the mean can't be computed in one time)
  # | i1j1 | i2j1 | i3j1 | 
  # +------+------+------+
  case $region in (WAM)
    imin11=$imin         && imax11=$(( imin11+(imax-imin)/3 ))
    imin21=$((imax11+1)) && imax21=$(( imin21+(imax-imin)/3 ))
    imin31=$((imax21+1)) && imax31=$imax
    imin12=$imin11       && imax12=$imax11
    imin22=$imin21       && imax22=$imax21
    imin32=$imin31       && imax32=$imax31
    imin13=$imin11       && imax13=$imax11
    imin23=$imin21       && imax23=$imax21
    imin33=$imin31       && imax33=$imax31
    jmin11=$jmin         && jmax11=$(( jmin11+(jmax-jmin)/3 ))
    jmin12=$((jmax11+1)) && jmax12=$(( jmin12+(jmax-jmin)/3 ))
    jmin13=$((jmax12+1)) && jmax13=$jmax
    jmin21=$jmin11       && jmax21=$jmax11
    jmin22=$jmin12       && jmax22=$jmax12
    jmin23=$jmin13       && jmax23=$jmax13
    jmin31=$jmin11       && jmax31=$jmax11
    jmin32=$jmin12       && jmax32=$jmax12
    jmin33=$jmin13       && jmax33=$jmax13
    #echo "i  $imin11:$imax11 $imin21:$imax21 $imin31:$imax31"
    #echo "   $imin12:$imax12 $imin22:$imax22 $imin32:$imax32"
    #echo "   $imin13:$imax13 $imin23:$imax23 $imin33:$imax33"
    #echo "j  $jmax13:$jmin13 $jmax23:$jmin23 $jmax33:$jmin33"
    #echo "   $jmax12:$jmin12 $jmax22:$jmin22 $jmax32:$jmin32"
    #echo "   $jmax11:$jmin11 $jmax21:$jmin21 $jmax31:$jmin31"
    ;;
  esac

  #input/output file names

  input="$WLDreg$runnam.$yrIF.${tmoI}-$tmoF.nc"
  output=${input%$runnam*}$runnam.$region.${input#??$runnam\.}
                                                  #e.g.: WAF86.CGU.1987.10-12.nc

  #title
  title="MAR exp: $runnam - $tmp"

  #NCO execution
 
  echo
  echo "region: $tmp - i=${imin}:$imax j=${jmin}:$jmax [$output]"
  echo
  set -x
  case $region in
  (WAM) #regional mean => nco failure => ncks to extract + ncwa to average then
       #mean 1
  ncwa --fortran --overwrite -a x,y -d x,$((imin11-1)),$((imax11-1)) -d y,$((jmin11-1)),$((jmax11-1)) $input tmp1.nc
  [ $? -eq 0 ] && nci=1
       #mean 2
  ncwa --fortran --overwrite -a x,y -d x,$((imin12-1)),$((imax12-1)) -d y,$((jmin12-1)),$((jmax12-1)) $input tmp2.nc
  [ $? -eq 0 ] && nci=$((nci+1))
       #mean 3
  ncwa --fortran --overwrite -a x,y -d x,$((imin13-1)),$((imax13-1)) -d y,$((jmin13-1)),$((jmax13-1)) $input tmp3.nc
  [ $? -eq 0 ] && nci=$((nci+1))
       #mean 4
  ncwa --fortran --overwrite -a x,y -d x,$((imin21-1)),$((imax21-1)) -d y,$((jmin21-1)),$((jmax21-1)) $input tmp4.nc
  [ $? -eq 0 ] && nci=$((nci+1))
       #mean 5
  ncwa --fortran --overwrite -a x,y -d x,$((imin22-1)),$((imax22-1)) -d y,$((jmin22-1)),$((jmax22-1)) $input tmp5.nc
  [ $? -eq 0 ] && nci=$((nci+1))
       #mean 6
  ncwa --fortran --overwrite -a x,y -d x,$((imin23-1)),$((imax23-1)) -d y,$((jmin23-1)),$((jmax23-1)) $input tmp6.nc
  [ $? -eq 0 ] && nci=$((nci+1))
       #mean 7
  ncwa --fortran --overwrite -a x,y -d x,$((imin31-1)),$((imax31-1)) -d y,$((jmin31-1)),$((jmax31-1)) $input tmp7.nc
  [ $? -eq 0 ] && nci=$((nci+1))
       #mean 8
  ncwa --fortran --overwrite -a x,y -d x,$((imin32-1)),$((imax32-1)) -d y,$((jmin32-1)),$((jmax32-1)) $input tmp8.nc
  [ $? -eq 0 ] && nci=$((nci+1))
       #mean 9
  ncwa --fortran --overwrite -a x,y -d x,$((imin33-1)),$((imax33-1)) -d y,$((jmin33-1)),$((jmax33-1)) $input tmp9.nc
  [ $? -eq 0 ] && nci=$((nci+1))
       #somme(mean i)/9
  if [ $nci -eq 9 ] ; then
    ncea --fortran --overwrite -y avg tmp[123456789].nc $output
    [ $? -eq 0 ] && nci=$((nci+1))
  fi
  [ $nci -eq 10 ] && rm -f tmp[123456789].nc
  [ $nci -eq 10 ] && ncatted -a title,global,m,c,"$title" $output
  set +x
  ;;
  (*)
  ncwa --fortran --overwrite -a x,y \
       -d x,$((imin-1)),$((imax-1)) \
       -d y,$((jmin-1)),$((jmax-1)) \
       $input $output
  [ $? -eq 0 ] && ncatted -a title,global,m,c,"$title" $output
  set +x
  ;;
  esac

done  #loop on regions

echo `date`
echo
