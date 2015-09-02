SUBROUTINE SUECAEOR

!**   OPTICAL THICKNESS OF ORGANIC-TYPE AEROSOLS

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SUECAEOR* FROM *SUECRAD*

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

      call SUECAEOR_01
      call SUECAEOR_02
      call SUECAEOR_03
      call SUECAEOR_04
      call SUECAEOR_05
      call SUECAEOR_06
      call SUECAEOR_07
      call SUECAEOR_08
      call SUECAEOR_09
      call SUECAEOR_10
      call SUECAEOR_11
      call SUECAEOR_12

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE SUECAEOR
