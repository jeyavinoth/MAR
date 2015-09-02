!OCL SCALAR
SUBROUTINE RRTM_KGB5

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB5_00
      call RRTM_KGB5_A1
      call RRTM_KGB5_A2
      call RRTM_KGB5_B1
      call RRTM_KGB5_B2
      call RRTM_KGB5_B3
      call RRTM_KGB5_B4

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB5
