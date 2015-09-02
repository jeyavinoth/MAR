C +-- Internal values of dBound subroutines (retain->next call):
C +   -------------------------------------

      REAL tM_MKE(ndb,mz), tMiMKE(ndb,mz), tMoMKE(ndb,mz)
      INTEGER NPiMKE(ndb,mz),NPoMKE(ndb,mz)
      REAL tM_MKE3(ndb,mz)
      REAL tM_MKE1(ndb,mz)

      REAL tM_SPdb(ndb), tM_CLdb(ndb)


      common /dBNDcom/tM_MKE,tMiMKE,tMoMKE,NPiMKE,NPoMKE,
     &                tM_MKE3,tM_MKE1,tM_SPdb,tM_CLdb
