
!  ==========
!  MAR_TU.inc
!  ==========

      real             TUkvm(mx,my,mz),TUkvh(mx,my,mz),TUmin(mz)
     .               ,TU_Pra(mx,my,mz)
      common /turver/  TUkvm          ,TUkvh          ,TUmin    ,TU_Pra
!                      TUkvm : Vertical Diffusion Coefficient (Momentum)  [m2/s]
!                      TUkvm : Vertical Diffusion Coefficient (Heat)      [m2/s]
!                      TUmin : Vertical Diffusion Coefficient (Minimum)   [m2/s]
!                     TU_Pra : Prandtl  Number (TUkvm/TUkvh)                 [-]

      real             TUkhx(mx,my,mz),TUkhff,TUkhmx
     .                ,TUkhy(mx,my,mz)
     .                ,TUspon(mzabso)
      common /turlat/  TUkhx          ,TUkhff,TUkhmx
     .                ,TUkhy
     .                ,TUspon
!                      TUkhx : Horizont.Diffusion Coefficient (x-direct.) [m2/s]
!                      TUkhy : Horizont.Diffusion Coefficient (y-direct.) [m2/s]
!                      TUspon: Horizont.Diffusion Coefficient (Top Sponge)[m2/s]
