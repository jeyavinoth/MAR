SUBROUTINE SUECAESD

!**   OPTICAL THICKNESS OF AEROSOLS OF SOIL-DUST ORIGIN

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SUECAESD* FROM *SUECRAD*

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

      call SUECAESD_01
      call SUECAESD_02
      call SUECAESD_03
      call SUECAESD_04
      call SUECAESD_05
      call SUECAESD_06
      call SUECAESD_07
      call SUECAESD_08
      call SUECAESD_09
      call SUECAESD_10
      call SUECAESD_11
      call SUECAESD_12

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE SUECAESD
