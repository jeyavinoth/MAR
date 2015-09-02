SUBROUTINE SUECAESU

!**   OPTICAL THICKNESS OF SULFATE-TYPE AEROSOLS

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SUECAESU* FROM *SUECRAD*

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

      call SUECAESU_01
      call SUECAESU_02
      call SUECAESU_03
      call SUECAESU_04
      call SUECAESU_05
      call SUECAESU_06
      call SUECAESU_07
      call SUECAESU_08
      call SUECAESU_09
      call SUECAESU_10
      call SUECAESU_11
      call SUECAESU_12

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE SUECAESU
