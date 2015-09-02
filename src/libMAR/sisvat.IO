          IF (.NOT.opnlog)                                          THEN
                   opnlog=.TRUE.
              open(unit=33,status='unknown',file='SISVAT.IO__95')
              open(unit=34,status='unknown',file='SISVAT.IO_105')
              rewind    33
              rewind    34
          ENDIF

          IF (ii__SV(ikl).EQ. 95.AND.
     .        jj__SV(ikl).EQ.  1.AND.
     .        nn__SV(ikl).EQ.  2)                                   THEN
            IF (mod(itexpe,20).EQ.0)                                THEN
              write(33,6100)
 6100         format(' Iter.    ii',7x,'u*   sqrrZ0',5x,'Z0e_',5x,'Z0_L'
     .                                              ,5x,'Z0_*',5x,'Z0s_'
     .                          ,5x,'Z0em',5x,'Z0sm',5x,'Z0m_',5x,'Z0h_'
     .                          ,5x,'Z(*)',4x,'f_eff',6x,'u*t')
            END IF
              write(33,6101)itexpe,ii__SV(ikl),      us__SV(ikl),sqrrZ0
     .      ,1.e3*Z0e_SV(ikl),1.e3*Z0mLnd     , 1.e3*Z0m_Sn
     .      ,1.e3*Z0s_SV(ikl)
     .      ,1.e3*Z0emSV(ikl),1.e3*Z0smSV(ikl), 1.e3*Z0m_SV(ikl)
     .      ,1.e3*Z0h_SV(ikl),     z_snsv(ikl),f_eff,usthSV(ikl)
 6101         format(2i6,13f9.3)
          END IF
 
          IF (ii__SV(ikl).EQ.105.AND.
     .        jj__SV(ikl).EQ.  1.AND.
     .        nn__SV(ikl).EQ.  1)                                   THEN
            IF (mod(itexpe,20).EQ.0)                                THEN
              write(34,6100)
            END IF
              write(34,6101)itexpe,ii__SV(ikl),      us__SV(ikl),sqrrZ0
     .      ,1.e3*Z0e_SV(ikl),1.e3*Z0mLnd     , 1.e3*Z0m_Sn
     .      ,1.e3*Z0s_SV(ikl)
     .      ,1.e3*Z0emSV(ikl),1.e3*Z0smSV(ikl), 1.e3*Z0m_SV(ikl)
     .      ,1.e3*Z0h_SV(ikl),     z_snsv(ikl),f_eff,usthSV(ikl)
          END IF
