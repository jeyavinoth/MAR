!OCL SCALAR
SUBROUTINE RRTM_KGB13

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB13_00
      call RRTM_KGB13_A1
      call RRTM_KGB13_A2

!     ------------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB13
