!OCL SCALAR
SUBROUTINE RRTM_KGB3

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB3_00
      call RRTM_KGB3_A1
      call RRTM_KGB3_A2
      call RRTM_KGB3_B1
      call RRTM_KGB3_B2
      call RRTM_KGB3_B3
      call RRTM_KGB3_B4

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB3
