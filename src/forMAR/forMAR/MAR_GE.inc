
C +  ==========
C +--MAR_GE.inc
C +  ==========

      integer*8        ou2sGE

      character*3      labmGE
      common /MARgea/  labmGE(0:12)
C +...  labmGE ......  Month Label 

      integer          iyr0GE,mma0GE,jda0GE,jhu0GE,jhl0GE,
     .                 iyrrGE,mmarGE,jdarGE,jhurGE,jhlrGE(mx,my),
     .                        minuGE,jsecGE,ksecGE,itizGE(mx,my),
     .                 njyrGE(0:12), njmoGE(0:12), nhyrGE,
     .                 njybGE(0:12), njmbGE(0:12)
      common /MARgei/  iyr0GE,mma0GE,jda0GE,jhu0GE,jhl0GE,
     .                 iyrrGE,mmarGE,jdarGE,jhurGE,jhlrGE       ,
     .                        minuGE,jsecGE,ksecGE,itizGE       ,
     .                 njyrGE      , njmoGE      , nhyrGE,
     .                 njybGE      , njmbGE
C +...  iyr0GE ......  Initial Year               
C +     iyrrGE ......  Current Year               
C +     mma0GE ......  Initial Month              
C +     mmarGE ......  Current Month              
C +     jda0GE ......  Initial Day                 
C +     jdarGE ......  Current Day                 
C +     jhu0GE ......  Initial Hour   (UT)  
C +     jhurGE ......  Current Hour   (UT)  
C +     jhl0GE ......  Initial Hour   (LT)  
C +     jhlrGE ......  Current Hour   (LT)  
C +     minuGE ......          Minute     
C +     jsecGE ......          second     
C +     ksecGE ......  Cent.of second     
C +     itizGE ......  Time Zone (numero du fuseau horaire)

C +...  njyrGE ......  Number of Days since   Begin         of the Year 
C +                                   before  Current Month
C +     njybGE ......  Leap Year Correction  to current Day of the Year
C +     njmoGE ......  Number of Days in      each    Month of the Year
C +     njmbGE ......  Leap Year Correction   to each Month of the Year
C +     nhyrGE ......  Number of Hours                      in one Year 

      real             GElat0,       GElon0,       GEddxx, 
     .                 GElatr       ,GElonh       ,GEtrue,
     .                 slatGE       ,clatGE       
      common /MARge1/  GElat0,       GElon0,       GEddxx, 
     .                 GElatr(mx,my),GElonh(mx,my),GEtrue,
     .                 slatGE(mx,my),clatGE(mx,my)
C +...  GElat0 ......  Latitude                              (Degrees)
C +     GElon0 ......  Longitude                             (Degrees)
C +     GElatr ......  Latitude                              (Radians)
C +     GElonh ......  Longitude                               (Hours)
C +...  GEtrue ......  On Oblique Stereographic Projection: 
C +                    Relative CoLatitude where distances   (Degrees)
C +                    on Proj.plane = distances on sphere
C +     GEddxx ......  Direction x-axis                      (Degrees)
C +     slatGE ......    Sine of Latitude
C +     clatGE ......  Cosine of Latitude

      integer                       n_azim
      parameter                    (n_azim=8)
      real             czenGE       ,               cz0_GE,
     .                 omenGE       ,               slopGE       ,       
     .                 cmntGE              ,daziGE, czmnGE       ,
     .                 cverGE       ,               rsunGE, 
     .                 tlsrGE       ,       tlssGE, tl__GE
      common /MARge2/  czenGE(mx,my),               cz0_GE,
     .                 omenGE(mx,my),               slopGE(mx,my),       
     .                 cmntGE(mx,my,n_azim),daziGE, czmnGE(mx,my),
     .                 cverGE(mx,my),               rsunGE, 
     .                 tlsrGE       ,       tlssGE, tl__GE
C +...  czenGE ......  Cosine  of Solar Zenith Angle
C +     cz0_GE ......  Cosine  of Solar Zenith Angle 
C +                   (Minimum Value for  solari    call)
C +     omenGE ......  Azimuth of Fall  Line                 (Radians)
C +     slopGE ......  Cosine  of Fall  Line   Angle
C +     cmntGE ......  Cosine  of Solar Zenith Angle 
C +                   (Minimum Value from Mountains Mask/
C +                    All Directions                   )
C +     daziGE ......  Unit Angle    for  Mountains Mask     (Radians)
C +     czmnGE ......  Cosine  of Solar Zenith Angle 
C +                   (Minimum Value from Mountains Mask)
C +     cverGE ......  Cosine  of Solar Normal Angle
C +     rsunGE ......  Insolation (Top of the Atmosphere)       (W/m2)
C +     tlsrGE ......  Sun Rise Local Time /  (iSND,jSND)       (hour)
C +     tlssGE ......  Sun Set  Local Time /  (iSND,jSND)       (hour)
C +     tl__GE ......  Current  Local Time /  (iSND,jSND)       (hour)

