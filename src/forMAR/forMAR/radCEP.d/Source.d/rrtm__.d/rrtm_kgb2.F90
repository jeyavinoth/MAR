!OCL SCALAR
SUBROUTINE RRTM_KGB2

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by  HGallée   , LGGE 

!     ------------------------------------------------------------------

      call RRTM_KGB2_01
      call RRTM_KGB2_02

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB2
