      character *3    explIO
      common /outlab/ explIO
C +...                explIO: Experiment Label
C +
      integer         igrdIO,
     .                jgrdIO,
     .                IO_gen,IO_loc,
     .                mxw1IO,mxw2IO,ixw_IO,
     .                myw1IO,myw2IO,iyw_IO,
     .                mzw1IO,mzw2IO,izw_IO,
     .                kkatIO,kmidIO
      common /outlpr/ igrdIO(5),
     .                jgrdIO(5),
     .                IO_gen,IO_loc,
     .                mxw1IO,mxw2IO,ixw_IO,
     .                myw1IO,myw2IO,iyw_IO,
     .                mzw1IO,mzw2IO,izw_IO,
     .                kkatIO,kmidIO
C +...igrdIO: i (x-direc.) Index Ref. Grid Point (for detailed Output)
C +   jgrdIO: j (y-direc.) Index Ref. Grid Point (for detailed Output)
C +
C +...IO_loc = 1 : Initialisation: General: sigma, z_amsl, (ug,vg)
C +            2 : Initialisation: Details: Soundings   (Ref.)
C +                Insolation              (PHYrad_top)       (Synop)
C +                Radiative      Transfer (PHYrad_int) (Ref. ,Synop)
C +                Turbulent      Transfer (TURtke_gen) (Ref.!,Synop)
C +                Polynyas                (SRFmod_pol) (RES. ,Synop)
C +                Snow Pack               (SRFmod_sno) (Ref. ,Synop)
C +                Surface Energy Balance  (SRFfrm_xxx) (Ref. ,Synop)
C +                Surface Mass   Balance  (OUTice)           (Synop)
C +            3 : Blowing Snow            (HYDblo)     (Ref.!,Synop)
C +                Cloud   mPhy            (HYDmic)     (Ref. ,Synop)
C +                Polynyas                (SRFmod_pol)       (Synop)
C +            7 : Surface Turbulent Flux  (TURsbl)
C +
C +   mxw1IO,mxw2IO,ixw_IO: x-domain to be printed (from mxw1 to mxw2 by ixw)
C +   myw1IO,myw2IO,iyw_IO: y-domain to be printed (from myw1 to myw2 by iyw)
C +   mzw1IO,mzw2IO,izw_IO: z-domain to be printed (from mzw1 to mzw2 by izw)
C +   kkatIO              : Low Level scaned in the meso animation 
C +   kmidIO              : Mid Level scaned in the meso animation
