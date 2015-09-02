C +	VARIABLES USED IN THE SECOND ORDER TURBULENCE CLOSURE (LEVEL 2.5)
C +	*****************************************************************
C +
C +
	real    uu_TE(mx,my,mz),vv_TE(mx,my,mz),ww_TE(mx,my,mz),
C +	........^ uu_TE    : x-velocity variance (k+1/2)
C +	........^ vv_TE    : y-velocity variance (k+1/2)
C +	........^ ww_TE    : z-velocity variance (k+1/2)
C +
     .		ml_TE(mx,my,mz)
C +	........^ ml_TE    : mixing length (k+1/2)
C +
c #WT.         ,mix_TE(mx,my,mz),shearT(mx,my,mz),buoyTE(mx,my,mz)
C +     ........^ output variables
C +
C +
        common  / turbul / uu_TE,vv_TE,ww_TE,ml_TE
c #WT.                    ,mix_TE,shearT,buoyTE
