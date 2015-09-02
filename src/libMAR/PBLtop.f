      subroutine PBLtop(TKE_1D,HHH_1D,h__PSL,h__SSL)

C +---------------------------------------------------------------------------+
C |                                                               07-APR-2005 |
C |   subroutine PBLtop computes the height of the Primary   Seeing Layer PSL |
C |                                                Secondary Seeing Layer SSL |
C |                                                                           |
C |   INPUT:   TKE_1D: Turbulent Kinetic Energy                       [m2/s2] |
C |            HHH_1D: Height above the Surface                           [m] |
C |                                                                           |
C |   OUTPUT:  h__PSL: Height of the Primary   Seeing Layer               [m] |
C |            h__SSL: Height of the Secondary Seeing Layer               [m] |
C |                                                                           |
C +---------------------------------------------------------------------------+


      IMPLICIT NONE

      include "See_nc.inc"

      real     TKE_1D(mz)
      real     HHH_1D(mz)
      real     h__PSL
      real     h__SSL

      real     TKEmin
      real     TKEtop
      
      integer  k     ,kmx   ,kzi

      logical  RESET
      logical  INTERP


      DATA     TKEmin/1.e-6/


C +--Height of the Primary   Seeing Layer (PSL)
C +  ==========================================

C +--Search the lowest TKE maximum
C +  -----------------------------

            k           =              mz
            TKEtop      =  0.01*TKE_1D(k)
 1001 CONTINUE
            k           =              k-1
      IF   (k          .LE.     mzabso     )                  GO TO 1000
      IF   (TKE_1D(k)  .LT.     TKE_1D(k+1).AND.
     .      TKE_1D(k+1).GT.     TKEmin*3.00)                  GO TO 1000
C +                                    3.00     = 1/2 order of magnitude
C +        (in order to only detect a significant maximum)
                                                              GO TO 1001
 1000 CONTINUE
            kmx         =              k+1
            TKEtop      =  0.01*TKE_1D(kmx)
            TKEtop      =   max(TKEmin*1.50,TKEtop)
C +                                    1.50     = 1/4 order of magnitude


C +--Search (from above) the lowest TKE minimum above the lowest TKE maximum
C +  ------ (This mimimum may be                                ) ----------
C +         (either a  TRUE         minimum  => INTERP = .FALSE.)
C +         (    or an arbitrary small value => INTERP = .TRUE. )
C +         -----------------------------------------------------

C +--Index  of the layer containing the minimum
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            kzi         =       mzabso
      DO k= mzabso,kmx
        IF (TKE_1D(k) .LT.      TKEtop          .OR.
     .      TKE_1D(k) .LT.      TKE_1D(k-1)*0.3     )               THEN
            kzi         =              k
         IF(TKE_1D(k) .LT.      TKEtop              )               THEN
            INTERP      =      .TRUE.
         ELSE
            INTERP      =      .FALSE.
         END IF
        END IF
      ENDDO

C +--Height of the                      minimum
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            k           =       kzi     
      IF   (kzi       .LE.      mzabso+1)                           THEN
            h__PSL      =       HHH_1D(mz)
      ELSE
        IF (INTERP)                                                 THEN
            h__PSL      =       HHH_1D(k+1)
     .                        +(HHH_1D(k)  -HHH_1D(k+1))
     .                        *(TKEtop     -TKE_1D(k+1))
     .                        /(TKE_1D(k)  -TKE_1D(k+1))
        ELSE
            h__PSL      =       HHH_1D(k)
        END IF
      END IF

            h__PSL      =   min(h__PSL     ,HHH_1D(1)  )
            h__PSL      =   max(HHH_1D(mz) ,h__PSL     )


C +--Height of the Secondary Seeing Layer (SSL)
C +  ==========================================

            RESET       =      .TRUE.


C +--Search the        TKE minimum above the Primary Seeing Layer (PSL)
C +  (necessary if the TKE has decreased below the minimum value)
C +  ------------------------------------------------------------------

      IF   (INTERP)                                                 THEN
            k           =       kzi + 1 
 1011   CONTINUE
            k           =       k-1
        IF (k         .LE.      mzabso     )                  GO TO 1010
        IF (TKE_1D(k) .LT.      TKE_1D(k+1))                  GO TO 1011
 1010   CONTINUE
      ELSE
            k           =       kzi
      END IF


C +--Search the first  TKE maximum above the Primary Seeing Layer (PSL)
C +  ------------------------------------------------------------------

            kmx         =       kzi
            k           =       k+1
 1021   CONTINUE
            k           =       k-1
      IF   (k            .LE.   mzabso)                       GO TO 1020
        IF (TKE_1D(k)    .GT.   TKE_1D(k-1)     .AND.
     .      TKE_1D(k)    .GT.   TKE_1D(k+1)     .AND.
     .      TKE_1D(k)    .GT.   TKEmin*3.0           )              THEN
C +                                    3.0      = 1/2 order of magnitude
C +        (in order to only detect a significant maximum)


C +--Define the TKE at the SSL top from the largest maximum in the SSL
C +  (thus examine the remaining upper part of the atmospheric column)
C +  -----------------------------------------------------------------

         IF(RESET)                                                  THEN
            RESET       =      .FALSE. ! indicates TKEtop is initialized
            TKEtop      =       0.00   !
         END IF
         IF(TKEtop       .LT.   TKE_1D(k)  *0.01)                   THEN
            TKEtop      =       TKE_1D(k)  *0.01
            kmx         =              k
         END IF
        END IF
                                                              GO TO 1021
 1020 CONTINUE
            TKEtop      =   max(TKEmin*3.0 ,TKEtop)
C +                                    3.0      = 1/2 order of magnitude


C +--Search (from above) the SSL top            above the SSL    TKE maximum
C +  ------ (This         may be                                ) ----------
C +         (either a  TRUE         minimum  => INTERP = .FALSE.)
C +         (    or an arbitrary small value => INTERP = .TRUE. )
C +         -----------------------------------------------------

            kzi         =       mzabso
      DO k= kmx,mzabso,-1     
        IF (TKE_1D(k) .GT.      TKEtop)
     .      kzi         =       k
      ENDDO

            k           =       kzi   -1
      IF   (kzi       .LE.      mzabso+1)                           THEN
            h__SSL      =       HHH_1D(mz)
      ELSE
        IF (INTERP)                                                 THEN
            h__SSL      =       HHH_1D(k+1)
     .                        +(HHH_1D(k)  -HHH_1D(k+1))
     .                        *(TKEtop     -TKE_1D(k+1))
     .                        /(TKE_1D(k)  -TKE_1D(k+1))
        ELSE
            h__SSL      =       HHH_1D(k)
        END IF
      END IF

            h__SSL      =   min(h__SSL     ,HHH_1D(1)  )
            h__SSL      =   max(h__PSL     ,h__SSL     )

      RETURN
      END
