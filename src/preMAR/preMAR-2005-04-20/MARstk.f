      program MARstk
C +
C +------------------------------------------------------------------------+
C | MAR Verification                                       25-09-2001  MAR |
C |     STAK Variables Analysis                                            |
C |                                                                        |
C |                                                                        |
C |                                                                        |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
      IMPLICIT NONE
C +
      integer      n1,n2,i1,i2,ic,jc,kc,lc,is,js,ls,lv,nv
      integer      ilins0,ilinsb,no_var
C +
      character*80 MARlin(100000),curlin
      character*10 sublin,sub000 ,subbak
      character* 6 MARvar(1000)  ,curvar
C +
      open( unit=1,status='old',file='MAR___.f')
      rewind     1
      open( unit=2,status='old',file='MARstk.L')
      rewind     2
      open( unit=3,status='new',file='MARstk.out')
      rewind     3
C +
C +
C +--INPUT
C +  =====
C +
      n1 = 0
 100  CONTINUE
      n1 = n1 + 1
      read(1,1010,END=101) MARlin(n1)
 1010 format(a80)
      GO TO  100
 101  CONTINUE
      n1 = n1 - 1
C +
      n2 = 0
 200  CONTINUE
      n2 = n2 + 1
      read(2,2010,END=201) MARvar(n2)
 2010 format(a6)
      GO TO  200
 201  CONTINUE
      n2 = n2 - 1
C +
      write(6,250) n1,n2
 250  format(/,'END of  INPUT:',i6,' MAR lines',i6,' STAK Variables')
C +
              nv         = 0            ! Nb of effective STAK Variables
C +
C +
C +--OUTPUT
C +  ======
C +
      DO      i2=1,n2
          write(3,6000) MARvar(i2)
 6000     format(/,'Variable ',a6)
          curvar =      MARvar(i2)
C +
              jc=0
 210          CONTINUE
              jc=jc+1
              IF (jc           .gt. 6 ) GO TO 211
              IF (curvar(jc:jc).eq.' ') GO TO 211
                                        GO TO 210
 211          CONTINUE
              lc=jc-1
C +
              is         = 0
              ls         = 0            ! ls: switch:=1  when 1st occurence
              lv         = 0            ! lv: switch:=1  when 2nd occurence
              subbak     ='          '
              ilins0     = 0
              sublin     ='          '
              ilinsb     = 0
C +
        DO    i1=1,n1
              curlin     = MARlin(i1)
          DO  ic=7,30
            IF  (curlin(ic:ic+9).eq.'subroutine'.AND.
     .           curlin( 1:   5).eq.'     '          )             THEN 
              is     = is + 1
              subbak = sublin
              sublin = curlin(ic+11:ic+20)
C +
              jc=0
 300          CONTINUE
              jc=jc+1
              IF (jc           .gt.10 ) GO TO 302
              IF (sublin(jc:jc).eq.'(') GO TO 301
                                        GO TO 300
 301          CONTINUE
              IF (jc           .le.10 ) THEN 
                DO kc=jc,10
                  sublin(kc:kc) =  ' ' 
                END DO
              END IF
 302          CONTINUE
C +
              ilinsb = i1
            END IF
          END DO
          DO  ic=7,30
            IF  (curlin(ic:ic+7).eq.'function'  .AND.
     .           curlin( 1:   5).eq.'     '          )             THEN 
              is     = is + 1
              subbak = sublin
              sublin = curlin(ic+ 8:ic+17)
C +
              jc=0
 310          CONTINUE
              jc=jc+1
              IF (jc           .gt.10 ) GO TO 312
              IF (sublin(jc:jc).eq.'(') GO TO 311
                                        GO TO 310
 311          CONTINUE
              IF (jc           .le.10 ) THEN 
                DO kc=jc,10
                  sublin(kc:kc) =  ' ' 
                END DO
              END IF
 312          CONTINUE
C +
              ilinsb = i1
            END IF
          END DO
                                                           no_var = 0
          DO  ic=7,75
            IF ( curlin(ic:ic)      .eq.'!'              ) no_var = 1
            IF ( curlin(ic:ic+lc-1 ).eq.curvar(1:lc).AND.
     .          (curlin(ic- 1:ic- 1).eq.' '.OR.
     .           curlin(ic- 1:ic- 1).eq.','.OR.
     .           curlin(ic- 1:ic- 1).eq.'+'.OR.
     .           curlin(ic- 1:ic- 1).eq.'-'.OR.
     .           curlin(ic- 1:ic- 1).eq.'*'.OR.
     .           curlin(ic- 1:ic- 1).eq.'/'.OR.
     .           curlin(ic- 1:ic- 1).eq.'('.OR.
     .           curlin(ic- 1:ic- 1).eq.')'    )    .AND.
     .          (curlin(ic+lc:ic+lc).eq.' '.OR.
     .           curlin(ic+lc:ic+lc).eq.','.OR.
     .           curlin(ic+lc:ic+lc).eq.'+'.OR.
     .           curlin(ic+lc:ic+lc).eq.'-'.OR.
     .           curlin(ic+lc:ic+lc).eq.'*'.OR.
     .           curlin(ic+lc:ic+lc).eq.'/'.OR.
     .           curlin(ic+lc:ic+lc).eq.'('.OR.
     .           curlin(ic+lc:ic+lc).eq.')'    )    .AND.
     .           curlin( 1:       5).eq.'     '     .AND.
     .           no_var             .eq. 0               )         THEN
              IF (ls.eq. 0)                                        THEN
                  js =  is          ! js: routine No when 1st occurence
                  ls =   1          ! ls: switch:=1  when 1st occurence
                  sub000 = sublin   ! sub000: subrout.lab.1st occurence
                  ilins0 = ilinsb   ! ilins0: Line Number subrou.sub000
               IF(ilins0.eq.0)        sub000='.MAIN.    '
              END IF
              IF (is.ne.js)                                        THEN 
               IF(lv.eq. 0)                                        THEN
                  lv =   1          ! lv: switch:=1  when 2nd occurence
                  nv =   1 + nv     ! nv: Nb of effective STAK Variables
                  write(3,6001) ilins0,sub000
 6001             format('IN (',i6,')',a10)
               END IF
                  write(3,6002) ilinsb,sublin,i1,curlin
 6002             format('IN (',i6,')',a10,
     .                     ' (line',i6,')',a80)
              END IF
            END IF
          END DO
        END DO
      END DO
C +
                  write(3,6003) nv
                  write(6,6003) nv
 6003             format(/,'END of OUTPUT:',i22,
     .                    ' STAK Variables (effective)')
C +
      close(unit=1)
      close(unit=2)
      close(unit=3)
C +
      end
