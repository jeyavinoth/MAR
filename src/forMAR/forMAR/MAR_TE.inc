
C +  ==========
C +--MAR_TE.inc
C +  ==========

      logical       ini_KA_TE
      common/tkelog/ini_KA_TE

      real          ect_TE          ,zi__TE
      common/tketke/ect_TE(mx,my,mz),zi__TE(mx,my)
C +...              ect_TE(i,j,k) : turbulent kinetic energy (i,j,k+1/2) (m2/s2)
C +                 zi__TE(i,j)   : Inversion Height         (i,j)       (m)

      real          eps_TE          ,edt_TE
      common/tkeeps/eps_TE(mx,my,mz),edt_TE(mx,my)
C +...              eps_TE(i,j,k) : dissipation of T.K.E.    (i,j,k+1/2) (m2/s3)
C +                 edt_TE(i,j)   : minimum dissipation time (i,j)       (s)

      real          tranTE
      common/tketra/tranTE(mx,my,mz)
C +...              tranTE(i,j,k) : transport   of T.K.E.    (i,j,k+1/2) (m2/s3)

