SUBROUTINE SUECAEBC

!**   OPTICAL THICKNESS OF BLACK CARBON AEROSOLS (URBAN/FOREST FIRE ORIGIN)

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SUECAEBC* FROM *SUECRAD*

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

      call suecaebc_01
      call suecaebc_02
      call suecaebc_03
      call suecaebc_04
      call suecaebc_05
      call suecaebc_06
      call suecaebc_07
      call suecaebc_08
      call suecaebc_09
      call suecaebc_10
      call suecaebc_11
      call suecaebc_12

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE SUECAEBC
