!OCL SCALAR
SUBROUTINE RRTM_KGB15

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB15_00
      call RRTM_KGB15_A1
      call RRTM_KGB15_A2

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB15
