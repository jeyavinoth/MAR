      real          CUspxh    ,CUspxb
      common /cubx/ CUspxh(mx),CUspxb(mx)
      real          CUspyh    ,CUspyb
      common /cuby/ CUspyh(my),CUspyb(my)
C +...Xh and Xb are cubic spline auxiliary variables
C +   used for HORIZONTAL advection (blank->x,w->y)
C +
      real          CUspzh    ,CUspzb
      common /cubz/ CUspzh(mz),CUspzb(mz)
C +...Xh and Xb are cubic spline auxiliary variables 
C +   used for VERTICAL   advection (routine cubver)
