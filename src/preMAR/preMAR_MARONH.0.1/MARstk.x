#!/bin/sh
  datey=`date +%y`
  datem=`date +%m`
  dated=`date +%d`
  dateH=`date +%H`
  min05=`echo 05`
  min10=`echo 10`
  min15=`echo 15`
  min20=`echo 20`
  min25=`echo 25`
  heure=`expr $dateH - 6`
  if  expr $heure \< 10
  then
            zero0=`echo 0`
            heure=`echo $zero0$heure`
  fi
#
f77 -o MARstk -O3 -finit-local-zero MARstk.f
rm  -f                              MARstk.out
       MARstk
#
  dateT=`echo $datem$dated$heure$min05$datey`
         echo $dateT
        touch $dateT                MARstk.f
  dateT=`echo $datem$dated$heure$min10$datey`
         echo $dateT
        touch $dateT                MARstk.x
# dateT=`echo $datem$dated$heure$min15$datey`
#        echo $dateT
#       touch $dateT                MARstk
        rm -f                       MARstk
  dateT=`echo $datem$dated$heure$min20$datey`
         echo $dateT
        touch $dateT                MARstk.out
