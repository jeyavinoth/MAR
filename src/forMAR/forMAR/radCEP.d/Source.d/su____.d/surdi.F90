SUBROUTINE SURDI(IYR,RCP)


!**** *SURDI*   - INITIALIZE COMMON YOERDI CONTROLLING RADINT

!     PURPOSE.
!     --------
!           INITIALIZE YOERDI, THE COMMON THAT CONTROLS THE
!           RADIATION INTERFACE

!**   INTERFACE.
!     ----------
!        CALL *SURDI* FROM *SURAD*
!              ------        -----

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOERDI

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS MODEL

!     AUTHOR.
!     -------
!        Original  JEAN-JACQUES MORCRETTE  *ECMWF*
!        Modified   P. Viterbo   99-03-26    Tiling of the land surface

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 88-12-15
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE YOERDI   , ONLY : RRAE     ,RALBSEAD ,&
            &RALBICEVS_AR,RALBICENI_AR,RALBICEVS_AN,RALBICENI_AN,&
            &RALBSFO  ,REMISD   ,REMISL   ,REMISN   ,REMISS   ,&
            &RCARDI   ,RCH4     ,RN2O     ,RO3      ,RCFC11   ,&
            &RCFC12   ,REPALB   ,REPCLC   ,REPH2O   ,RSUNDUR

IMPLICIT NONE


!     LOCAL REAL SCALARS
REAL_B :: ZAIRMWG, ZC11MWG, ZC12MWG, ZCH4MWG, ZCO2MWG, ZN2OMWG, ZO3MWG
INTEGER_M :: JM, IM,IYR
character*5 :: RCP
REAL :: TMP_RCP(35)

!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

RRAE = 0.1277E-02_JPRB

!* Threshold for computing sunshine duration (W/m2)
RSUNDUR=120._JPRB

!*  Ocean surface albedo for diffuse radiation (Taylor et al., 1997)      
RALBSEAD = 0.06_JPRB
!*  For sea ice, monthly values are based on Ebert and Curry, 1993, Table 2.
!   We take dry snow albedo as the representative value for non-summer
!   months, and bare sea-ice as the representative value for summer
!   months. The values for Antarctic are shifted six-months.
!*  Sea ice surf. albedo for visible rad. (snow covered; Ebert and Curry, 1993)
RALBICEVS_AR(1:12) = (/0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB,&
                      &0.975_JPRB,0.876_JPRB,0.778_JPRB,0.778_JPRB,&
                      &0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB/)
!*  Sea ice surf. albedo for near-IR rad. (snow covered; Ebert and Curry, 1993)
RALBICENI_AR(1:12) = (/0.664_JPRB,0.664_JPRB,0.664_JPRB,0.664_JPRB,&
                      &0.664_JPRB,0.476_JPRB,0.288_JPRB,0.288_JPRB,&
                      &0.664_JPRB,0.664_JPRB,0.664_JPRB,0.664_JPRB/)
DO JM=1,12
  IM=MOD(JM+5,12)+1
  RALBICEVS_AN(JM)=RALBICEVS_AR(IM)
  RALBICENI_AN(JM)=RALBICENI_AR(IM)
ENDDO
!*  Snow albedo in the presence of high vegetation
RALBSFO=0.15_JPRB

!- sea surface emissivity and other surfaces outside the window region     
REMISS  = 0.99_JPRB
!- snow window emissivity      
REMISN  = 0.98_JPRB
!- land window emissivity      
REMISL  = 0.96_JPRB
!- desert window emissivity (lower bound when dry)      
REMISD  = 0.93_JPRB

ZAIRMWG = 28.970_JPRB
ZCO2MWG = 44.011_JPRB
ZCH4MWG = 16.043_JPRB
ZN2OMWG = 44.013_JPRB
ZO3MWG  = 47.9982_JPRB
ZC11MWG = 137.3686_JPRB
ZC12MWG = 120.9140_JPRB


!*  Concentration of the various trace gases (IPCC/SACC values for 1990)
!        CO2         CH4        N2O        CFC11       CFC12
!      353ppmv     1.72ppmv   310ppbv     280pptv     484pptv


RCARDI  = 353.E-06_JPRB*ZCO2MWG/ZAIRMWG
RCH4    = 1.72E-06_JPRB*ZCH4MWG/ZAIRMWG
RN2O    = 310.E-09_JPRB*ZN2OMWG/ZAIRMWG
RO3     =   1.E-06_JPRB*ZO3MWG /ZAIRMWG
RCFC11  = 280.E-12_JPRB*ZC11MWG/ZAIRMWG
RCFC12  = 484.E-12_JPRB*ZC12MWG/ZAIRMWG

OPEN(unit=10,status="old",file=TRIM(RCP)//"_CMIP5.dat")
DO JM=1,3
 read(10,*)
ENDDO
DO IM=1,IYR-1765
 read(10,*)
ENDDO
READ(10,*) JM,(tmp_RCP(IM),IM=1,35)
CLOSE(10)

if(JM==IYR) then
 RCARDI  = tmp_RCP(1)      *1.E-06_JPRB*ZCO2MWG/ZAIRMWG
 RCH4    = tmp_RCP(4)/1000.*1.E-06_JPRB*ZCH4MWG/ZAIRMWG
 RN2O    = tmp_RCP(5)      *1.E-09_JPRB*ZN2OMWG/ZAIRMWG
 RCFC11  = tmp_RCP(20)     *1.E-12_JPRB*ZC11MWG/ZAIRMWG
 RCFC12  = tmp_RCP(21)     *1.E-12_JPRB*ZC12MWG/ZAIRMWG
else
 print *,"ERROR in "//TRIM(RCP)//"_CMIP5.dat" 
 stop
endif

1000 continue

print *, "Concentration of the various trace gases for "//TRIM(RCP)
print *, "Year =",IYR
print *, "CO2  =",RCARDI
print *, "CH4  =",RCH4
print *, "N2O  =",RN2O
print *, "CFC11=",RCFC11
print *, "CFC12=",RCFC12

REPCLC=1.E-12_JPRB
REPH2O=1.E-12_JPRB
REPALB=1.E-12_JPRB

!     -----------------------------------------------------------------

RETURN
END SUBROUTINE SURDI
