!OCL SCALAR
SUBROUTINE RRTM_KGB5_00

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gall�e   , LGGE  (splitting)

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE YOERRTO5 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO, CCL4O
USE YOERRTA5 , ONLY : STRRAT1   ,STRRAT2

!     ------------------------------------------------------------------


IMPLICIT NONE
CCL4O( :) = (/&
     &26.1407_JPRB,  53.9776_JPRB,  63.8085_JPRB,  36.1701_JPRB,&
     &15.4099_JPRB, 10.23116_JPRB,  4.82948_JPRB,  5.03836_JPRB,&
     &1.75558_JPRB,  _ZERO_     ,  _ZERO_     ,  _ZERO_     ,&
     &_ZERO_     ,  _ZERO_     ,  _ZERO_     ,  _ZERO_      /)

STRRAT1 = 90.4894_JPRB
STRRAT2 = 0.900502_JPRB

!     ------------------------------------------------------------------

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

SELFREFO( :, 1) = (/&
&1.27664E-01_JPRB, 1.09296E-01_JPRB, 9.35703E-02_JPRB, 8.01075E-02_JPRB, 6.85818E-02_JPRB,&
&5.87143E-02_JPRB, 5.02666E-02_JPRB, 4.30343E-02_JPRB, 3.68426E-02_JPRB, 3.15417E-02_JPRB/)
SELFREFO( :, 2) = (/&
&1.39620E-01_JPRB, 1.20381E-01_JPRB, 1.03793E-01_JPRB, 8.94908E-02_JPRB, 7.71595E-02_JPRB,&
&6.65273E-02_JPRB, 5.73601E-02_JPRB, 4.94562E-02_JPRB, 4.26414E-02_JPRB, 3.67656E-02_JPRB/)
SELFREFO( :, 3) = (/&
&1.42628E-01_JPRB, 1.23043E-01_JPRB, 1.06148E-01_JPRB, 9.15726E-02_JPRB, 7.89986E-02_JPRB,&
&6.81511E-02_JPRB, 5.87931E-02_JPRB, 5.07201E-02_JPRB, 4.37556E-02_JPRB, 3.77474E-02_JPRB/)
SELFREFO( :, 4) = (/&
&1.53569E-01_JPRB, 1.33143E-01_JPRB, 1.15435E-01_JPRB, 1.00082E-01_JPRB, 8.67706E-02_JPRB,&
&7.52299E-02_JPRB, 6.52241E-02_JPRB, 5.65491E-02_JPRB, 4.90279E-02_JPRB, 4.25070E-02_JPRB/)
SELFREFO( :, 5) = (/&
&1.70491E-01_JPRB, 1.46448E-01_JPRB, 1.25796E-01_JPRB, 1.08056E-01_JPRB, 9.28182E-02_JPRB,&
&7.97290E-02_JPRB, 6.84856E-02_JPRB, 5.88278E-02_JPRB, 5.05319E-02_JPRB, 4.34059E-02_JPRB/)
SELFREFO( :, 6) = (/&
&1.76394E-01_JPRB, 1.51432E-01_JPRB, 1.30003E-01_JPRB, 1.11606E-01_JPRB, 9.58127E-02_JPRB,&
&8.22542E-02_JPRB, 7.06144E-02_JPRB, 6.06217E-02_JPRB, 5.20431E-02_JPRB, 4.46784E-02_JPRB/)
SELFREFO( :, 7) = (/&
&1.85706E-01_JPRB, 1.59172E-01_JPRB, 1.36429E-01_JPRB, 1.16936E-01_JPRB, 1.00228E-01_JPRB,&
&8.59068E-02_JPRB, 7.36322E-02_JPRB, 6.31114E-02_JPRB, 5.40939E-02_JPRB, 4.63648E-02_JPRB/)
SELFREFO( :, 8) = (/&
&1.88647E-01_JPRB, 1.61657E-01_JPRB, 1.38529E-01_JPRB, 1.18710E-01_JPRB, 1.01726E-01_JPRB,&
&8.71722E-02_JPRB, 7.47005E-02_JPRB, 6.40132E-02_JPRB, 5.48549E-02_JPRB, 4.70068E-02_JPRB/)
SELFREFO( :, 9) = (/&
&1.90074E-01_JPRB, 1.62793E-01_JPRB, 1.39427E-01_JPRB, 1.19415E-01_JPRB, 1.02275E-01_JPRB,&
&8.75959E-02_JPRB, 7.50233E-02_JPRB, 6.42552E-02_JPRB, 5.50327E-02_JPRB, 4.71338E-02_JPRB/)
SELFREFO( :,10) = (/&
&1.94769E-01_JPRB, 1.66338E-01_JPRB, 1.42057E-01_JPRB, 1.21320E-01_JPRB, 1.03611E-01_JPRB,&
&8.84863E-02_JPRB, 7.55696E-02_JPRB, 6.45384E-02_JPRB, 5.51175E-02_JPRB, 4.70718E-02_JPRB/)
SELFREFO( :,11) = (/&
&1.90624E-01_JPRB, 1.64229E-01_JPRB, 1.41488E-01_JPRB, 1.21896E-01_JPRB, 1.05017E-01_JPRB,&
&9.04757E-02_JPRB, 7.79475E-02_JPRB, 6.71542E-02_JPRB, 5.78554E-02_JPRB, 4.98442E-02_JPRB/)
SELFREFO( :,12) = (/&
&1.90502E-01_JPRB, 1.64025E-01_JPRB, 1.41228E-01_JPRB, 1.21599E-01_JPRB, 1.04699E-01_JPRB,&
&9.01472E-02_JPRB, 7.76181E-02_JPRB, 6.68303E-02_JPRB, 5.75419E-02_JPRB, 4.95444E-02_JPRB/)
SELFREFO( :,13) = (/&
&1.86786E-01_JPRB, 1.61636E-01_JPRB, 1.39872E-01_JPRB, 1.21039E-01_JPRB, 1.04741E-01_JPRB,&
&9.06380E-02_JPRB, 7.84338E-02_JPRB, 6.78729E-02_JPRB, 5.87340E-02_JPRB, 5.08256E-02_JPRB/)
SELFREFO( :,14) = (/&
&1.99149E-01_JPRB, 1.71475E-01_JPRB, 1.47646E-01_JPRB, 1.27129E-01_JPRB, 1.09462E-01_JPRB,&
&9.42512E-02_JPRB, 8.11538E-02_JPRB, 6.98764E-02_JPRB, 6.01662E-02_JPRB, 5.18053E-02_JPRB/)
SELFREFO( :,15) = (/&
&2.02676E-01_JPRB, 1.73701E-01_JPRB, 1.48869E-01_JPRB, 1.27587E-01_JPRB, 1.09347E-01_JPRB,&
&9.37144E-02_JPRB, 8.03170E-02_JPRB, 6.88348E-02_JPRB, 5.89941E-02_JPRB, 5.05603E-02_JPRB/)
SELFREFO( :,16) = (/&
&1.99865E-01_JPRB, 1.72699E-01_JPRB, 1.49225E-01_JPRB, 1.28942E-01_JPRB, 1.11416E-01_JPRB,&
&9.62721E-02_JPRB, 8.31866E-02_JPRB, 7.18797E-02_JPRB, 6.21097E-02_JPRB, 5.36676E-02_JPRB/)      

FRACREFAO( :, 1) = (/&
!     From P = 387.6 mb.
    &0.13966499_JPRB,0.14138900_JPRB,0.13763399_JPRB,0.13076700_JPRB,&
    &0.12299100_JPRB,0.10747700_JPRB,0.08942000_JPRB,0.06769200_JPRB,&
    &0.04587610_JPRB,0.00501173_JPRB,0.00415809_JPRB,0.00328398_JPRB,&
    &0.00240015_JPRB,0.00156222_JPRB,0.00059104_JPRB,0.00008323_JPRB/)
FRACREFAO( :, 2) = (/&
    &0.13958199_JPRB,0.14332899_JPRB,0.13785399_JPRB,0.13205400_JPRB,&
    &0.12199700_JPRB,0.10679600_JPRB,0.08861080_JPRB,0.06712320_JPRB,&
    &0.04556030_JPRB,0.00500863_JPRB,0.00416315_JPRB,0.00328629_JPRB,&
    &0.00240023_JPRB,0.00156220_JPRB,0.00059104_JPRB,0.00008323_JPRB/)
FRACREFAO( :, 3) = (/&
    &0.13907100_JPRB,0.14250501_JPRB,0.13889600_JPRB,0.13297300_JPRB,&
    &0.12218700_JPRB,0.10683800_JPRB,0.08839260_JPRB,0.06677310_JPRB,&
    &0.04538570_JPRB,0.00495402_JPRB,0.00409863_JPRB,0.00328219_JPRB,&
    &0.00240805_JPRB,0.00156266_JPRB,0.00059104_JPRB,0.00008323_JPRB/)
FRACREFAO( :, 4) = (/&
    &0.13867700_JPRB,0.14190100_JPRB,0.13932300_JPRB,0.13327099_JPRB,&
    &0.12280800_JPRB,0.10692500_JPRB,0.08844510_JPRB,0.06658510_JPRB,&
    &0.04519340_JPRB,0.00492276_JPRB,0.00408832_JPRB,0.00323856_JPRB,&
    &0.00239289_JPRB,0.00155698_JPRB,0.00059104_JPRB,0.00008323_JPRB/)
FRACREFAO( :, 5) = (/&
    &0.13845000_JPRB,0.14158800_JPRB,0.13929300_JPRB,0.13295600_JPRB,&
    &0.12348300_JPRB,0.10736700_JPRB,0.08859480_JPRB,0.06650610_JPRB,&
    &0.04498230_JPRB,0.00491335_JPRB,0.00406968_JPRB,0.00322901_JPRB,&
    &0.00234666_JPRB,0.00155235_JPRB,0.00058813_JPRB,0.00008323_JPRB/)
FRACREFAO( :, 6) = (/&
    &0.13837101_JPRB,0.14113200_JPRB,0.13930500_JPRB,0.13283101_JPRB,&
    &0.12349200_JPRB,0.10796400_JPRB,0.08890490_JPRB,0.06646480_JPRB,&
    &0.04485990_JPRB,0.00489554_JPRB,0.00405264_JPRB,0.00320313_JPRB,&
    &0.00234742_JPRB,0.00151159_JPRB,0.00058438_JPRB,0.00008253_JPRB/)
FRACREFAO( :, 7) = (/&
    &0.13834500_JPRB,0.14093500_JPRB,0.13896500_JPRB,0.13262001_JPRB,&
    &0.12326900_JPRB,0.10828900_JPRB,0.08950050_JPRB,0.06674610_JPRB,&
    &0.04476560_JPRB,0.00489624_JPRB,0.00400962_JPRB,0.00317423_JPRB,&
    &0.00233479_JPRB,0.00148249_JPRB,0.00058590_JPRB,0.00008253_JPRB/)
FRACREFAO( :, 8) = (/&
    &0.13831300_JPRB,0.14069000_JPRB,0.13871400_JPRB,0.13247600_JPRB,&
    &0.12251400_JPRB,0.10831300_JPRB,0.08977090_JPRB,0.06776920_JPRB,&
    &0.04498390_JPRB,0.00484111_JPRB,0.00398948_JPRB,0.00316069_JPRB,&
    &0.00229741_JPRB,0.00150104_JPRB,0.00058608_JPRB,0.00008253_JPRB/)
FRACREFAO( :, 9) = (/&
    &0.14027201_JPRB,0.14420401_JPRB,0.14215700_JPRB,0.13446601_JPRB,&
    &0.12303700_JPRB,0.10596100_JPRB,0.08650370_JPRB,0.06409570_JPRB,&
    &0.04312310_JPRB,0.00471110_JPRB,0.00393954_JPRB,0.00310850_JPRB,&
    &0.00229588_JPRB,0.00146366_JPRB,0.00058194_JPRB,0.00008253_JPRB/)

FRACREFBO( :, 1) = (/&
!     From P = 1.17 mb.
    &0.14339100_JPRB,0.14358699_JPRB,0.13935301_JPRB,0.13306700_JPRB,&
    &0.12135700_JPRB,0.10590600_JPRB,0.08688240_JPRB,0.06553220_JPRB,&
    &0.04446740_JPRB,0.00483580_JPRB,0.00399413_JPRB,0.00316225_JPRB,&
    &0.00233007_JPRB,0.00149135_JPRB,0.00056246_JPRB,0.00008059_JPRB/)
FRACREFBO( :, 2) = (/&
    &0.14330500_JPRB,0.14430299_JPRB,0.14053699_JPRB,0.13355300_JPRB,&
    &0.12151200_JPRB,0.10529100_JPRB,0.08627630_JPRB,0.06505230_JPRB,&
    &0.04385850_JPRB,0.00476555_JPRB,0.00395010_JPRB,0.00313878_JPRB,&
    &0.00232273_JPRB,0.00149354_JPRB,0.00056246_JPRB,0.00008059_JPRB/)
FRACREFBO( :, 3) = (/&
    &0.14328399_JPRB,0.14442700_JPRB,0.14078601_JPRB,0.13390100_JPRB,&
    &0.12132600_JPRB,0.10510600_JPRB,0.08613660_JPRB,0.06494630_JPRB,&
    &0.04381310_JPRB,0.00475378_JPRB,0.00394166_JPRB,0.00313076_JPRB,&
    &0.00231235_JPRB,0.00149159_JPRB,0.00056301_JPRB,0.00008059_JPRB/)
FRACREFBO( :, 4) = (/&
    &0.14326900_JPRB,0.14453100_JPRB,0.14114200_JPRB,0.13397101_JPRB,&
    &0.12127200_JPRB,0.10493400_JPRB,0.08601380_JPRB,0.06483360_JPRB,&
    &0.04378900_JPRB,0.00474655_JPRB,0.00393549_JPRB,0.00312583_JPRB,&
    &0.00230686_JPRB,0.00148433_JPRB,0.00056502_JPRB,0.00008059_JPRB/)
FRACREFBO( :, 5) = (/&
    &0.14328900_JPRB,0.14532700_JPRB,0.14179000_JPRB,0.13384600_JPRB,&
    &0.12093700_JPRB,0.10461500_JPRB,0.08573010_JPRB,0.06461340_JPRB,&
    &0.04366570_JPRB,0.00473087_JPRB,0.00392539_JPRB,0.00311238_JPRB,&
    &0.00229865_JPRB,0.00147572_JPRB,0.00056517_JPRB,0.00007939_JPRB/)

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB5_00
