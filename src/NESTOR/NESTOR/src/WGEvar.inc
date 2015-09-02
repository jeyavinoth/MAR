
C +---Number of intervals for statistics of wind distribution
C +   -------------------------------------------------------

      INTEGER gc

      PARAMETER (gc = 14)


C +---General variables
C +   -----------------

      REAL WGEest(mx,my   ),WGEupb(mx,my   ),WGElwb(mx,my   ),
     .     WGEsta(mx,my,gc),WNDsta(mx,my,gc),WGE_zi(mx,my   )

C +...WGEest : wind gust estimate
C +...WGElwb : lower bound on gust estimates
C +...WGEupb : upper bound on gust estimates
C +...WGEsta : gust distribution
C +...WNDsta : wind distribution
C +...WGE_zi : inversion height


      COMMON/WGEvar_r/WGEest,WGEupb,WGElwb,WGEsta,WNDsta,WGE_zi 
