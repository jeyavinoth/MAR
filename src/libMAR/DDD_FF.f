      subroutine DDD_FF(u_wind,v_wind,lon   ,lat   ,x  ,y
     .                 ,Dir__X,trulat,lat__0,lon__0,FFF,DDD)

C +----------------------------------------------------------------------------+
C |                                                                            |
C |   subroutine DDD_FF computes the wind speed and direction                  |
C |                              from stereographic MAR wind components        |
C |                                                                            |
C +----------------------------------------------------------------------------+


      IMPLICIT NONE


      REAL     u_wind,v_wind
      REAL     lon   ,lat
      REAL     x     ,y   
      REAL     Dir__X,trulat
      REAL     lat__0,lon__0
      REAL     cosbet,sinbet
      REAL     up    ,vp
      REAL     sinla0,cosla0
      REAL     coslon,sinlon
      REAL     coslat,sinlat
      REAL     denomi
      REAL     cosalp,sinalp
      REAL     uu    ,vv    ,magnif
      REAL     FFF   ,DDD


C +--Restauration Coordonnées ds un plan où (x,y)=(E,N) au Centre de Projection

      cosbet=cos((Dir__X-90.)*3.1415/180.)
      sinbet=sin((Dir__X-90.)*3.1415/180.)
      up    =     u_wind *cosbet + v_wind *sinbet
      vp    =    -u_wind *sinbet + v_wind *cosbet

C +--Coordonnées Géographiques de chaque Point

      sinla0=sin( lat__0     *3.1415/180.)
      cosla0=cos( lat__0     *3.1415/180.)

      coslon=cos((lon-lon__0)*3.1415/180.) ! lat,lon: coord.   géograph. d'un pt
      sinlon=sin((lon-lon__0)*3.1415/180.) !   x,  y: coord.stéréograph. d'un pt

      sinlat=(sinla0*(4*6371*6371-x*x-y*y)+4*6371*y*cosla0)
     .      /        (4*6371*6371+x*x+y*y)
      coslat=sqrt(1-sinlat*sinlat)

C +--Rotation locale par rapport au système d'axes (E,N)

      denomi= 1.+sinla0*sinlat+cosla0*coslat *coslon
      cosalp=(cosla0*coslat+(1+sinla0*sinlat)*coslon)/denomi
      sinalp=(0            -(  sinlat+sinla0)*sinlon)/denomi

      uu=up*cosalp-vp*sinalp
      vv=up*sinalp+vp*cosalp

C +--Direction du vecteur Vent

          DDD=270.-(180./3.1415)*atan2(vv,uu)
      IF (DDD     .LT.        0.)                                   THEN
          DDD    = DDD    + 360.
      ENDIF
      IF (DDD     .GT.      360.)                                   THEN
          DDD    = DDD    - 360.
      ENDIF

C +--Vitesse   du Vecteur Vent

      magnif=(1+cos(trulat*3.14159/180.))/(1-sin(lat *3.1415/180.))
      FFF   =  sqrt(u_wind*u_wind+v_wind*v_wind)

      RETURN
      END
