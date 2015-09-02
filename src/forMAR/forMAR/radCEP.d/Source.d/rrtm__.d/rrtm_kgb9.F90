!OCL SCALAR
SUBROUTINE RRTM_KGB9

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE  (splitting)

!     ------------------------------------------------------------------

      call RRTM_KGB9_00
      call RRTM_KGB9_A1
      call RRTM_KGB9_A2
      call RRTM_KGB9_BB

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB9
