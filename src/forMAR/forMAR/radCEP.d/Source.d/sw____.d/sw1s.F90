SUBROUTINE SW1S &
 &( KIDIA , KFDIA , KLON , KLEV , KAER , KNU &
 &, PAER  , PALBD , PALBP, PCG  , PCLD , PCLEAR &
 &, PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD  &
 &, PFD   , PFU   , PCD  , PCU  , PSUDU1 &
 &)

!**** *SW1S* - SHORTWAVE RADIATION, FIRST SPECTRAL INTERVAL

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW1S* IS CALLED FROM *SW*.


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES QUANTITIES FOR THE CLEAR-SKY FRACTION OF THE
!     COLUMN
!          2. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
!     CONTINUUM SCATTERING
!          3. MULTIPLY BY OZONE TRANSMISSION FUNCTION

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWTT*, *SWUVO3*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
!        96-01-15   J.-J. MORCRETTE    SW in nsw SPECTRAL INTERVALS 
!        990128     JJMorcrette        sunshine duration
!        99-05-25   JJMorcrette        Revised aerosols
!        00-12-18   JJMorcrette        6 spectral intervals

!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE YOESW    , ONLY : RRAY     ,RSUN
USE YOERAD   , ONLY : NSW


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KAER
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KKIND
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KNU



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PAER(KLON,6,KLEV)&
  &,  PALBD(KLON,NSW)      , PALBP(KLON,NSW)&
  &,  PCG(KLON,NSW,KLEV)   , PCLD(KLON,KLEV) &
  &,  PCLEAR(KLON)&
  &,  PDSIG(KLON,KLEV)&
  &,  POMEGA(KLON,NSW,KLEV), POZ(KLON,KLEV)&
  &,  PRMU(KLON)           , PSEC(KLON)&
  &,  PTAU(KLON,NSW,KLEV)  , PUD(KLON,5,KLEV+1)

REAL_B :: PFD(KLON,KLEV+1) , PFU(KLON,KLEV+1)&
  &,  PCD(KLON,KLEV+1)     , PCU(KLON,KLEV+1)&
  &,  PSUDU1(KLON)

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

INTEGER_M :: IIND6(6), IIND4(4)

REAL_B :: ZCGAZ(KLON,KLEV)&
  &,  ZDIFF(KLON)        , ZDIRF(KLON)        &
  &,  ZDIFT(KLON)        , ZDIRT(KLON)        &
  &,  ZPIZAZ(KLON,KLEV)&
  &,  ZRAYL(KLON), ZRAY1(KLON,KLEV+1), ZRAY2(KLON,KLEV+1)&
  &,  ZREFZ(KLON,2,KLEV+1)&
  &,  ZRJ(KLON,6,KLEV+1), ZRJ0(KLON,6,KLEV+1)&
  &,  ZRK(KLON,6,KLEV+1), ZRK0(KLON,6,KLEV+1)&
  &,  ZRMUE(KLON,KLEV+1), ZRMU0(KLON,KLEV+1)&
  &,  ZR6(KLON,6)       , ZR4(KLON,4)&
  &,  ZTAUAZ(KLON,KLEV)&
  &,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)&
  &,  ZTRCLD(KLON)      , ZTRCLR(KLON)&
  &,  ZW6(KLON,6)       , ZW4(KLON,4), ZO(KLON,2) ,ZT(KLON,2) 

!     LOCAL INTEGER SCALARS
INTEGER_M :: IKL, IKM1, JAJ, JK, JL


!     ------------------------------------------------------------------

!*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
!                 ----------------------- ------------------


!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------


DO JL = KIDIA,KFDIA
  ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)&
   &* (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)&
   &* (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))
ENDDO
!print *,'SW1S After Rayleigh'


!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------


!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------


CALL SWCLR &
  &( KIDIA  , KFDIA , KLON  , KLEV , KAER , KNU &
  &, PAER   , PALBP , PDSIG , ZRAYL, PSEC &
  &, ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
  &, ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
  &)
!print *,'SW1S After SWCLR'


!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------


CALL SWR &
  &( KIDIA ,KFDIA ,KLON  ,KLEV  , KNU &
  &, PALBD ,PCG   ,PCLD  ,POMEGA, PSEC , PTAU &
  &, ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 , ZREFZ, ZRJ  ,ZRK , ZRMUE &
  &, ZTAUAZ,ZTRA1 ,ZTRA2 ,ZTRCLD &
  &)
!print *,'SW1S After SWR'


!     ------------------------------------------------------------------

!*         3.    OZONE ABSORPTION
!                ----------------

IF (NSW <= 4) THEN

!*         3.1   TWO OR FOUR SPECTRAL INTERVALS
!                ------------------------------

  IIND6(1)=1
  IIND6(2)=2
  IIND6(3)=3
  IIND6(4)=1
  IIND6(5)=2
  IIND6(6)=3


!*         3.1.1  DOWNWARD FLUXES
!                 ---------------


  JAJ = 2

  DO JL = KIDIA,KFDIA
    ZW6(JL,1)=_ZERO_
    ZW6(JL,2)=_ZERO_
    ZW6(JL,3)=_ZERO_
    ZW6(JL,4)=_ZERO_
    ZW6(JL,5)=_ZERO_
    ZW6(JL,6)=_ZERO_
    PFD(JL,KLEV+1)=((_ONE_-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
     &+ PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)
    PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
  ENDDO
  DO JK = 1 , KLEV
    IKL = KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZW6(JL,1)=ZW6(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW6(JL,2)=ZW6(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
      ZW6(JL,3)=ZW6(JL,3)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW6(JL,4)=ZW6(JL,4)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW6(JL,5)=ZW6(JL,5)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
      ZW6(JL,6)=ZW6(JL,6)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ENDDO
    
    KKIND=6
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, KKIND &
      &, IIND6 &
      &, ZW6  &
      &, ZR6                          )

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR6(JL,1)*ZR6(JL,2)*ZR6(JL,3)*ZRJ(JL,JAJ,IKL)
      ZDIRF(JL) = ZR6(JL,4)*ZR6(JL,5)*ZR6(JL,6)*ZRJ0(JL,JAJ,IKL)
      PFD(JL,IKL) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&
       &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
      PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    ZDIFT(JL) = ZR6(JL,1)*ZR6(JL,2)*ZR6(JL,3)*ZTRCLD(JL)
    ZDIRT(JL) = ZR6(JL,4)*ZR6(JL,5)*ZR6(JL,6)*ZTRCLR(JL)
    PSUDU1(JL) = ((_ONE_-PCLEAR(JL)) * ZDIFT(JL)&
     &+PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)
  ENDDO


!*         3.1.2  UPWARD FLUXES
!                 -------------


  DO JL = KIDIA,KFDIA
    PFU(JL,1) = ((_ONE_-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)&
     &+ PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))&
     &* RSUN(KNU)
    PCU(JL,1) = ZDIRF(JL) * PALBP(JL,KNU) * RSUN(KNU)
  ENDDO

  DO JK = 2 , KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZW6(JL,1)=ZW6(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW6(JL,2)=ZW6(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW6(JL,3)=ZW6(JL,3)+POZ(JL,  IKM1)*1.66_JPRB
      ZW6(JL,4)=ZW6(JL,4)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW6(JL,5)=ZW6(JL,5)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW6(JL,6)=ZW6(JL,6)+POZ(JL,  IKM1)*1.66_JPRB
    ENDDO
    
    KKIND=6
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, KKIND &
      &, IIND6 &
      &, ZW6  &
      &, ZR6                          )
  
    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR6(JL,1)*ZR6(JL,2)*ZR6(JL,3)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR6(JL,4)*ZR6(JL,5)*ZR6(JL,6)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&
       &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
      PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO




ELSE IF (NSW == 6) THEN
!print *,'SW1S ozone 6SI'

!*         3.2   SIX SPECTRAL INTERVALS
!                ----------------------

  IIND4(1)=1
  IIND4(2)=2
  IIND4(3)=1
  IIND4(4)=2


!*         3.2,1  DOWNWARD FLUXES
!                 ---------------


  JAJ = 2

  DO JL = KIDIA,KFDIA
    ZW4(JL,1)=_ZERO_
    ZW4(JL,2)=_ZERO_
    ZW4(JL,3)=_ZERO_
    ZW4(JL,4)=_ZERO_
  
    ZO(JL,1)=_ZERO_
    ZO(JL,2)=_ZERO_
    PFD(JL,KLEV+1)=((_ONE_-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
      &+ PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)
    PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
  ENDDO
  DO JK = 1 , KLEV
    IKL = KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZW4(JL,1)=ZW4(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW4(JL,2)=ZW4(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
      ZW4(JL,3)=ZW4(JL,3)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW4(JL,4)=ZW4(JL,4)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
    
      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ENDDO
 
    KKIND=4
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, KKIND &
      &, IIND4 &
      &, ZW4  &
      &, ZR4  &
      & )
!    print *,'SW1S after SWTT1 JK=',JK  

    KKIND=2
    CALL SWUVO3 ( KIDIA, KFDIA, KLON, KNU, KKIND &
      &, ZO  &
      &, ZT  &
      & )
!    print *,'SW1S after SWUVO3 JK=',JK

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR4(JL,1)*ZR4(JL,2)*ZT(JL,1)*ZRJ(JL,JAJ,IKL)
      ZDIRF(JL) = ZR4(JL,3)*ZR4(JL,4)*ZT(JL,2)*ZRJ0(JL,JAJ,IKL)
      PFD(JL,IKL) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&
        &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
      PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    ZDIFT(JL) = ZR4(JL,1)*ZR4(JL,2)*ZT(JL,1)*ZTRCLD(JL)
    ZDIRT(JL) = ZR4(JL,3)*ZR4(JL,4)*ZT(JL,2)*ZTRCLR(JL)
    PSUDU1(JL) = ((_ONE_-PCLEAR(JL)) * ZDIFT(JL)&
      &+PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)
  ENDDO


!*         3.2.2  UPWARD FLUXES
!                 -------------


  DO JL = KIDIA,KFDIA
    PFU(JL,1) = ((_ONE_-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)&
      &+ PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))&
      &* RSUN(KNU)
    PCU(JL,1) = ZDIRF(JL) * PALBP(JL,KNU) * RSUN(KNU)
  ENDDO

  DO JK = 2 , KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZW4(JL,1)=ZW4(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW4(JL,2)=ZW4(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW4(JL,3)=ZW4(JL,3)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW4(JL,4)=ZW4(JL,4)+PUD(JL,2,IKM1)*1.66_JPRB
      
      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKM1)*1.66_JPRB
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKM1)*1.66_JPRB
    ENDDO

    KKIND=4
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, KKIND &
      &, IIND4 &
      &, ZW4  &
      &, ZR4  &
      & )

    KKIND=2
    CALL SWUVO3 ( KIDIA, KFDIA, KLON, KNU, KKIND &
      &, ZO  &
      &, ZT  &
      & )

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR4(JL,1)*ZR4(JL,2)*ZT(JL,1)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR4(JL,3)*ZR4(JL,4)*ZT(JL,2)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&
        &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
      PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO
  
END IF  

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SW1S
