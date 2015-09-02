C +----------------------------------------------------------------------+
C | MAR post-processing                                          2002.02 |
C | Int3D                                                                |
C | RLS 3D var to MAR grid                                               |
c 
c   code incorrect car basŽ sur VerP, qui ne peut s'appliquer ici
C +----------------------------------------------------------------------+

      SUBROUTINE Int3D (idRLS, idMAR, itM, itR,
     $           MARlon, MARlat, VarN_RLS,OUT_RLS)
     
      IMPLICIT NONE
      
C +...* dimensions :
      include 'NSTdim.inc'
      include 'NSTtoMAP.inc'
            
C +...* MAPOST specific (nreg...)
      include 'MAPOST.inc'
      include 'globals.inc'


C +---INPUT
C +   ~~~~~
      INTEGER idMAR, itM, idRLS, itR
      CHARACTER*(1) VarN_RLS
      REAL MARlon(mx,my), MARlat(mx,my)
      
C +---OUTPUT 
C +   ~~~~~~
      REAL OUT_RLS(mx,my,mz)
      
C +---LOCAL
C +   ~~~~~
      CHARACTER *10 var_units, tmp_units
C +-  -Coordonnees:
      REAL LSlon(LSni), LSlat(LSnj), LSlev(LSnk)
      REAL MARx(mx), MARy(my), empty1(1), sigma(mz)
      REAL ptop, ptopPa
      REAL CSTp (LSnk), SIGp (LSnk)
C +-  -Valeurs lues et intermédiaires de calcul:
      REAL wk3MAR(mx,my,mz)
      REAL wk3MIX(mx,my,LSnk), wkMIX(mx,my)
      REAL wk3RLS(LSni,LSnj,LSnk), wkRLS(LSni,LSnj)
      REAL spIRLS(mx,my), splocRLS
      REAL valINT      , pINT
      REAL valRLS(LSnk), pRLS(LSnk)
      REAL va3INT(mx,my,mz)
      INTEGER ilv,ii,jj,kk

C +---Read the MAR coordinate and data
C +   --------------------------------

      CALL UNread
     &   (idMAR, 'tairDY', itM, 0, 1,1,mx,my,mz,
     &    MARx, MARy, sigma, var_units, wk3MAR)  
C +...  NB: FileID VarName time 0=3D subregion #levels
     
      ptop   = 0.01
      ptopPa = ptop * 1000.
C +   ...MAR top pressure (en kPa = MAR, et en Pa)
      
C +---Read the RLS coordinate and data
C +   --------------------------------

      CALL UNread
     &   (idRLS,'CSTp', 0, 0, 0, 0,
     &    LSnk    ,1      , 1     ,
     &    LSlev   ,empty1 , empty1,
     &    tmp_units, CSTp)
C +...constant pressure contrib. to coordinate, in Pa

      CALL UNread
     &   (idRLS,'SIGp', 0, 0, 0, 0,
     &    LSnk     ,1      , 1     ,
     &    LSlev    ,empty1 , empty1,
     &    tmp_units, SIGp)

      CALL UNread
     &   (idRLS, 'SP', itR, 0, 1,1,LSni,LSnj,1,
     &    LSlon, LSlat, empty1, tmp_units, wkRLS)
     
      CALL INThor (1,
     &    LSlon , LSlat , wkRLS,
     &    MARlon, MARlat, spIRLS)

      CALL UNread
     &   (idRLS, VarN_RLS, itR,  0, 1,1,LSni,LSnj,LSnk,
     &    LSlon, LSlat   , LSlev, var_units, wk3RLS)   
     
  
C +---Interpolate the RLS data to MAR grid.
C +   -------------------------------------             
      DO ilv = 1, LSnk
      DO jj  = 1, LSnj
      DO ii  = 1, LSni
        wkRLS(ii,jj)=wk3RLS(ii,jj,ilv)
      END DO
      END DO
        Call INThor (1,
     &    LSlon , LSlat , wkRLS,
     &    MARlon, MARlat, wkMIX)     
C +...  Horizontal interp. to MAR grid 
C +...  (RLS vertical grid)
      DO jj  = 1, my
      DO ii  = 1, mx
        wk3MIX(ii,jj,ilv)=wkMIX(ii,jj)
      END DO
      END DO
      END DO
 
      DO ii= 1,mx
      DO jj= 1,my

         splocRLS = spIRLS(ii,jj)
C +      ...RLS pressure thickness (MAR hor. grid), Pa

         DO ilv = 1, LSnk  
           pRLS  (ilv)= CSTp(ilv)+splocRLS*SIGp(ilv)
C +        ...pressure on the RLS levels, Pa
           valRLS(ilv)= wk3MIX(ii,jj,ilv)
C +        ...value of "VarN_RLS" on the RLS levels
         END DO


         DO ilv = 1, mz
           pINT= (splocRLS-ptopPa)*sigma(ilv)+ptopPa
C +        ...RLS pres at MAR-type sigma coord., Pa

           CALL INTlin(pRLS, valRLS,LSnk,
     .                 pINT, valINT)
           Out_RLS(ii,jj,ilv)= valINT
C +        ...Interpolation to MAR vert. grid.

         END DO
         
      ENDDO
      ENDDO               
      END      