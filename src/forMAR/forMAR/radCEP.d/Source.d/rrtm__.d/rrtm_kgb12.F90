!OCL SCALAR
SUBROUTINE RRTM_KGB12

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB12_00
      call RRTM_KGB12_A1
      call RRTM_KGB12_A2

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB12
