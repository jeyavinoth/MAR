!OCL SCALAR
SUBROUTINE RRTM_KGB4

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB4_00
      call RRTM_KGB4_A1
      call RRTM_KGB4_A2
      call RRTM_KGB4_B1
      call RRTM_KGB4_B2
      call RRTM_KGB4_B3
      call RRTM_KGB4_B4

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB4
