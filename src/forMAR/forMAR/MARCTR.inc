
C +  ==========
C +--MARCTR.inc
C +  ==========

      logical           reaVAR    ! Input INI: Prev.Dyn.Simulat. (MAR, GCM)
      logical           reaLBC    ! Input LBC: Prev.Dyn.Simulat. (MAR, GCM)
      logical           safVAR    ! Full Output on Saving Files   MARxxx.DAT
      logical           sALONE    ! 
      logical           geoNST    !
      logical           conmas    ! Mass       Conserv. Constraint Init.Switch
      logical           potvor    ! P.-Vortic. Conserv. Constraint Init.Switch
      logical           hamfil    ! Initial Filtered Fields (Time, Hamming)
      logical           brocam    ! Brown and Campana  Time      Scheme Switch
      logical           LFrBAK    ! Leap-Frog Backward Advection Scheme Switch
      logical           openLB    ! 
      logical           sommlb    ! 
      logical           FirstC    ! 
      logical           qqmass    ! Mass       Conservation             Switch
      logical           turhor    ! Horizontal Diffusion (Smagorinsky)  Switch
      logical           SBLitr    ! 
      logical           tur_25    ! 
      logical           convec    ! Convective  Adjustment              Switch
      logical           MFLX_d    ! Convective  Adjustment (deep)       Switch
      logical           MFLX_s    ! Convective  Adjustment (shallow)    Switch
      logical           micphy    ! Cloud       Microphysics            Switch
      logical           fracld    ! Fractional  Cloudiness              Switch
      logical           chimod    ! Atmospheric Chemical Model          Switch
      logical           physic    ! Atmospheric/Surface Physics         Switch
      logical           polmod    ! Interactive Polynya    is           Switch
      logical           snomod    ! Interactive Snow Model is           Switch
      logical           BloMod    ! 
      logical           vegmod    ! Interactive SVAT                    Switch
      logical           VSISVAT   ! 
      logical           no_vec    ! Scalar     (NO Vectorization)       Switch

      common/CTRL_L_MAR/reaVAR,reaLBC,safVAR,sALONE,geoNST
     .                 ,conmas,potvor
     .                 ,hamfil
     .                 ,brocam,LFrBAK
     .                 ,openLB,sommlb,FirstC,qqmass
     .                 ,turhor,SBLitr,tur_25
     .                 ,convec,MFLX_d,MFLX_s,micphy,fracld
     .                 ,chimod
     .                 ,physic,polmod,snomod,BloMod,vegmod,VSISVAT
     .                 ,no_vec

      integer*4         itexpe
      integer           iterun,nterun,nbhour,itConv
      integer           iboucl,nboucl,nprint,ipr_nc,npr_nc
      integer           maptyp
      integer           log_1D
      common/CTRL_I_MAR/iterun,nterun,nbhour,itConv
     .                 ,iboucl,nboucl,nprint,ipr_nc,npr_nc
     .                 ,maptyp
     .                 ,log_1D
      common/CTRL_I_MAR2/itexpe

      real              Robert
      real              rrmin ,rrmax
      real              rxbase,rxfact
      real              fxlead
      real              tMFLXd        ! d(time) between 2 deep convection CALL
      real              aMFLXd,aMFLXs ! adjustment times (deep, shallow)
      common/CTRL_R_MAR/Robert
     .                 ,rrmin ,rrmax
     .                 ,rxbase,rxfact
     .                 ,fxlead
     .                 ,tMFLXd
     .                 ,aMFLXd,aMFLXs

      character*16      fnam
      common/CTRL_C_MAR/fnam
