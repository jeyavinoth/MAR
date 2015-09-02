!OCL SCALAR
SUBROUTINE RRTM_KGB7

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB7_00
      call RRTM_KGB7_A1
      call RRTM_KGB7_A2
      call RRTM_KGB7_BB

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB7
