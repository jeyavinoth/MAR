!OCL SCALAR
SUBROUTINE RRTM_KGB1

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gallée   , LGGE 

!     ------------------------------------------------------------------

      call RRTM_KGB1_01
      call RRTM_KGB1_02

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB1
