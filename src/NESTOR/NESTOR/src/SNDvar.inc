
C +---Sounding condition
C +   ------------------

      LOGICAL SNDing


C +---Ascii file containing the sounding
C +   ----------------------------------

      CHARACTER*57 SNDfil 


C +---Sounding location and time
C +   --------------------------

      INTEGER SNDyar,SNDmma,SNDjda,SNDjhu

C +...SNDyar : year  of the sounding
C +...SNDmma : month "       "
C +...SNDjda : day   "       "
C +...SNDjhu : hour  "       "

      REAL    SNDlat,SNDlon

C +...SNDlat,SNDlon : latitude and longitude


C +---Meteorological variables
C +   ------------------------

      INTEGER    mzm
      PARAMETER (mzm=100)

C +...Number of vertical levels of the sounding (read)

      INTEGER    ms
      PARAMETER (ms = 40)

C +...Number of vertical levels of the sounding (written)

      REAL SND_pp(0:mzm),SND_zz(0:mzm),
     .     SND_tt(0:mzm),SND_pt(0:mzm),
     .     SND_hr(0:mzm),SND_qv(0:mzm),
     .     SND_ff(0:mzm),SND_dd(0:mzm),
     .     SND_uu(0:mzm),SND_vv(0:mzm),
     .     SNDhyb(0:mzm)

C +...SND_pp : pressure levels      (kPa)
C +...SND_zz : level heights          (m)
C +...SND_tt : real temperature       (K)
C +...SND_pt : potential temperature  (K)
C +...SND_hr : relative humidity      (%)
C +...SND_qv : specific humidity  (kg/kg)
C +...SND_ff : wind velocity        (m/s)
C +...SND_dd : wind direction       (deg)
C +...SND_uu : U-component of wind  (m/s)
C +...SND_vv : V-component of wind  (m/s)
C +...SNDhyb : hybrid coordinates (interpolation)


C +---Temporary arrays
C +   ---------------
      REAL TMP_pt(0:mzm),TMP_qv(0:mzm),
     .     TMP_uu(0:mzm),TMP_vv(0:mzm),
     .     TMPhyb(0:mzm)


      COMMON/SNDvar_c/SNDfil 
      
      COMMON/SNDvar_i/SNDyar,SNDmma,SNDjda,SNDjhu
      
      COMMON/SNDvar_l/SNDing
      
      COMMON/SNDvar_r/SNDlat,SNDlon,SND_pp,SND_zz,SND_tt,SND_pt,
     .                SND_hr,SND_qv,SND_ff,SND_dd,SND_uu,SND_vv,
     .                SNDhyb  
           
