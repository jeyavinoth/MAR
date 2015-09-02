SUBROUTINE SUECAESS

!**   OPTICAL THICKNESS OF AEROSOLS OF SEA-SALT ORIGIN

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SUECAESS* FROM *SUECRAD*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     J.-J. MORCRETTE  E.C.M.W.F.    98/12/21

!     MODIFICATIONS.
!     --------------
!     H. GALLEE        L.G.G.E.      04/01/15:  split for the NEC SX5

!-----------------------------------------------------------------------

      call SUECAESS_01
      call SUECAESS_02
      call SUECAESS_03
      call SUECAESS_04
      call SUECAESS_05
      call SUECAESS_06
      call SUECAESS_07
      call SUECAESS_08
      call SUECAESS_09
      call SUECAESS_10
      call SUECAESS_11
      call SUECAESS_12

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE SUECAESS
