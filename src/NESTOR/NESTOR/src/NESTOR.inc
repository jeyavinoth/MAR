
C +   NESTOR.inc
C +   ----------

      CHARACTER*3      LABLio
	  CHARACTER*3      region !+ GRD for polar config          +!
	  CHARACTER*3      USRreg !+ user region = GRD,ANT,KER,... +!
      CHARACTER*4      TOP30
      CHARACTER*60     NSTdir
      CHARACTER*100    LSCfil
      
      INTEGER          RUNiyr,RUNmma,RUNjda,RUNjhu,DESGpr,I_time,
     .                 DURjda,DURjhu,FORjhu,HORint,VERint,TOPopt,
     .                 DATstp,SELECT,WGEmet,iter
      
      INTEGER*4        DATini,DAT_dt,DATfin,DATtim
 
      REAL             SVTwet,NDVmin(nbdom),NDVmax(nbdom)

      LOGICAL          LoutDA,LoutLS,VEGdat,VEGbel,TOPBE,
     .                 TOPetopo,GTOPO30,TOPcst,TOPcstLSC,TOPdomLSC,
     .                 TOPcst0,TOPfilt,CORzz6,CORsurf,SVTmod,
     .                 ASCfor,CLDcor,Sinclair,Alpert,VEGcor,NDV1km,
     .                 NDV8km,SSTrey,SVTlsc,RUGdat
  
    
      COMMON/NESTOR_c /LABLio,TOP30,NSTdir,LSCfil,region,USRreg   
          
      COMMON/NESTOR_i1/RUNiyr,RUNmma,RUNjda,RUNjhu,DESGpr,I_time,
     .                 DURjda,DURjhu,FORjhu,HORint,VERint,TOPopt,
     .                 DATstp,SELECT,WGEmet,iter
      
      COMMON/NESTOR_i2/DATini,DAT_dt,DATfin,DATtim
      
      COMMON/NESTOR_l /LoutDA,LoutLS,VEGdat,VEGbel,TOPBE,
     .                 TOPetopo,GTOPO30,TOPcst,TOPcstLSC,TOPdomLSC,
     .                 TOPcst0,TOPfilt,CORzz6,CORsurf,SVTmod,
     .                 ASCfor,CLDcor,Sinclair,Alpert,VEGcor,NDV1km,
     .                 NDV8km,SSTrey,SVTlsc,RUGdat
      
      COMMON/NESTOR_r /NDVmin,NDVmax,SVTwet            
