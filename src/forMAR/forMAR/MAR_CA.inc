C +     VARIABLES USED IN THE CONVECTIVE ADJUSTMENT SCHEME
C +     **************************************************
C +
C +
        integer adj_CA(mx,my),
C +     ........^ Variable designed in order to suppress the vertical
C +     ........  advection when the convective adjustement scheme is
C +     ........  active. This variable gives the number of remaining
C +     ........  time steps needed to remove the available buoyant
C +     ........  energy.
C +
     .           int_CA
C +     ........ ^ Time interval (transformed in number of iterations)
C +     ........   between two checks of convective unstability
C +
        common / adjCA / adj_CA,int_CA
C +
C +
        real      capeCA(mx,my),
C +     ......... ^ CAPE (Convective Available Potential Energy [J/kg])
C +
     .            dpktCA(mx,my,mz),
C +     ......... ^ Changes for Potential Temperature Vertical Profile
C +        
     .            dqv_CA(mx,my,mz),
C +     ......... ^ Changes for Water Vapor           Vertical Profile
C +
     .            dqw_CA(mx,my,mz),
C +     ......... ^ Cloud Droplets  Detrainment in the Environment
C +
     .            dqi_CA(mx,my,mz),
C +     ......... ^ Ice   Crystals  Detrainment in the Environment
C +
     .            drr_CA(mx,my),
C +     ......... ^ Changes in      Precipitation  reaching the Ground
C +
     .            dss_CA(mx,my),
C +     ......... ^ Changes in Snow Precipitation  reaching the Ground
C +
     .            dsn_CA(mx,my),
C +     ......... ^ Snow Melting                                  
C +
     .            rainCA(mx,my),
C +     ......... ^ Cumulated Convective       Precipitation
C +
     .            snowCA(mx,my),
C +     ......... ^ Cumulated Convective Solid Precipitation
C +
     .            tau_CA(mx,my)
C +     ......... ^ Characteristic Time Scale for Convection
C +
        common / tenCA / capeCA,dpktCA,dqv_CA,dqw_CA,dqi_CA,
     .                   drr_CA,dss_CA,dsn_CA,rainCA,snowCA,tau_CA
C +
