SUBROUTINE COL2BOX &
 & (KIDIA, KFDIA, KLON, KLEV, KBOX, KOVLP &
 & , PCLFR, PCLBX &
 & )
!
!* Subdivide a column of cloud parameters in a set of homogeneous boxes
!
!     from C.Jakob and S.A. Klein
!
!-----------------------------------------------------------------------

#include "tsmbkind.h"

INTEGER_M :: KIDIA
INTEGER_M :: KFDIA
INTEGER_M :: KLON
INTEGER_M :: KLEV
INTEGER_M :: KBOX
INTEGER_M :: KOVLP

REAL_B :: PCLFR(KLON,KLEV)
REAL_B :: PCLBX(KLON,100,KLEV)
            
!-- local
      
INTEGER_M :: IABOX(KLON,KBOX), IABOXM1(KLON,KBOX), IABOXINT(KLON,KBOX)
INTEGER_M :: IBOXTYPE1(KLON), IBOXTYPE2(KLON), IBOXTYPE3(KLON)
INTEGER_M :: ISUMBOX(KBOX), ISUMBOXM1(KBOX)

REAL_B :: ZTCC(KLON)      
       
                  
ZBOXWIDTH=1./FLOAT(KBOX)
ZAMIN =1.E-03
ZEPSEC=1.E-06
      
DO JB=1,KBOX
  DO JL=KIDIA,KFDIA
    IABOXINT(JL,JB)=0
    IABOXM1 (JL,JB)=0
    IABOX   (JL,JB)=0
  END DO
END DO 
DO JL=KIDIA,KFDIA
  ZTCC(JL)     =_ZERO_
  ISUMBOX(JL)  =_ZERO_
  ISUMBOXM1(JL)=_ZERO_
END DO                             
                                            
                                                               
DO JK=1,KLEV
!
  IF (JK.GT.1) THEN
    DO JB=1,KBOX
      DO JL=KIDIA,KFDIA
        IABOXM1(JL,JB)=IABOX(JL,JB)
        ISUMBOXM1(JL)=ISUMBOXM1(JL)+IABOX(JL,JB)
        IABOX(JL,JB)=0
      END DO
    END DO
  END IF
      
  DO JL=KIDIA,KFDIA
    ITCCM1=NINT(REAL(KBOX)*ZTCC(JL))           
    IF (ZTCC(JL).GT.ZAMIN .AND. ZTCC(JL).LT._HALF_*ZBOXWIDTH) THEN
      ITCCM1=1
    END IF
!
!-- various cloud overlap assumptions
!
    IF (JK.GT.1) THEN 
!
!-- maximum-random
!       
      IF (KOVLP.EQ.1) THEN
        ZTCC(JL) = _ONE_ - ( (_ONE_-ZTCC(JL)) &
         &      *(_ONE_ -MAX( PCLFR(JL,JK)  , PCLFR(JL,JK-1))) &
         &      /(_ONE_ -MIN( PCLFR(JL,JK-1), _ONE_-ZEPSEC))  )
!
!-- maximum
!       
      ELSE IF (KOVLP.EQ.2) THEN
        ZTCC(JL)=MAX(ZTCC(JL),PCLFR(JL,JK))
!
!-- random
!       
      ELSE IF (KOVLP.EQ.3) THEN 
        ZTCC(JL)=_ONE_-(_ONE_-ZTCC(JL))*(_ONE_-PCLFR(JL,JK))
      END IF
!
    ELSE
      ZTCC(JL)=PCLFR(JL,JK)
    END IF
!
    ITCC=NINT(REAL(KBOX)*ZTCC(JL))
    IF (ZTCC(JL).GT.ZAMIN .AND. ZTCC(JL).LT. _HALF_*ZBOXWIDTH) THEN
      ITCC=1
    END IF
    IAM1=ISUMBOXM1(JL)
    IA=NINT(REAL(KBOX)*PCLFR(JL,JK))
    IF (PCLFR(JL,JK).GT.ZAMIN & 
      &       .AND. PCLFR(JL,JK).LT. _HALF_*ZBOXWIDTH) THEN       
      IA=1
    END IF
!     
    IBOXTYPE1(JL)=ITCC-ITCCM1
!     IF (KOVLP.NE.3) THEN
    IBOXTYPE2(JL)=MIN( IAM1, IA-IBOXTYPE1(JL))
!     ELSE
!       IBOXTYPE2(JL)=NINT( FLOAT(IAM1)*FLOAT(IA-IBOXTYPE1(JL))
!     &                    /MAX(FLOAT(ITCCM1), ZEPSEC) )
!     END IF       
    IBOXTYPE3(JL)=IA - IBOXTYPE1(JL)-IBOXTYPE2(JL)
  END DO
!
  DO JB=1,KBOX
    DO JL=KIDIA,KFDIA
      IF (IABOXINT(JL,JB).EQ.0) THEN
        IF (IBOXTYPE1(JL).GT.0) THEN
          IABOX(JL,JB)=1
          IABOXINT(JL,JB)=1
          IBOXTYPE1(JL)=IBOXTYPE1(JL)-1
        END IF
      ELSE
        IF (IABOXM1(JL,JB).EQ.1) THEN
          IF (IBOXTYPE2(JL).GT.0) THEN
            IABOX(JL,JB)=1
            IBOXTYPE2(JL)=IBOXTYPE2(JL)-1
          END IF
        ELSE
          IF (IBOXTYPE3(JL).GT.0) THEN
            IABOX(JL,JB)=1
            IBOXTYPE3(JL)=IBOXTYPE3(JL)-1
          END IF
        END IF
      END IF
    END DO
  END DO
!     
  DO JB=1,KBOX
    DO JL=KIDIA,KFDIA
      IF (JB.EQ.1) THEN
        IBOXTYPE1(JL)=IBOXTYPE1(JL)+IBOXTYPE2(JL)+IBOXTYPE3(JL)
      END IF
      IF (IABOX(JL,JB).EQ.0 .AND. IBOXTYPE1(JL).GT.0) THEN
        IABOX(JL,JB)=1
        IBOXTYPE1(JL)=IBOXTYPE1(JL)-1
      END IF
      ISUMBOX(JL)=ISUMBOX(JL)+IABOX(JL,JB)
    END DO             
  END DO
  DO JL=KIDIA,KFDIA
    if (JK.GE.21) THEN
      PRINT 9001,(IABOX(JL,JB),JB=1,KBOX)
    end if  
    DO JB=1,KBOX
      PCLBX(JL,JB,JK)=FLOAT(IABOX(JL,JB))
    END DO  
9001 FORMAT(1X,100I1)        
  END DO
      
END DO
      
RETURN
END SUBROUTINE COL2BOX
