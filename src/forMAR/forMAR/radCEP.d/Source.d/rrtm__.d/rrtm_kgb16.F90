!OCL SCALAR
SUBROUTINE RRTM_KGB16

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB16_00
      call RRTM_KGB16_A1
      call RRTM_KGB16_A2

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB16
