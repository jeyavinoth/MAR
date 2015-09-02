      program dump3D_MAR

C +----------------------------------------------------------------------------+
C |   program dump3D_MAR output the maximum difference for a MAR variable      |
C |                             step by step                                   |
C |                                                                            |
C +----------------------------------------------------------------------------+


      IMPLICIT     NONE


      integer      i1    ,j1    ,k1
      integer      i2    ,j2    ,k2
      integer      i__max,j__max,k__max
      real         var1  ,var2
      real         v1_max,v2_max,difmax
      character*32 titr_1,titr_2,titbak


C +--INITIALIZATION
C +  ==============

      titbak = '                                '
      v1_max =  0.
      v2_max =  0.
      i__max =  0
      j__max =  0
      k__max =  0


C +--INPUT
C +  =====

      open (unit=1,status='old',file='dump3D_MAR.1')
      rewind     1
      open (unit=2,status='old',file='dump3D_MAR.2')
      rewind     2
      open (unit=3,status='NEW',file='dump3D_MAR.dif')
      rewind     3

 1000 CONTINUE
      read(1,1010,end=1001) titr_1,i1,j1,k1,var1
      read(2,1010,end=1001) titr_2,i2,j2,k2,var2
 1010 format(a32,3i6,f15.6)


C +--Verification
C +  ============

      IF (i1.NE.i2.OR.j1.NE.j2.OR.k1.NE.k2.OR.titr_1.NE.titr_2)     THEN
          write(6,60)i1,j1,k1,titr_1,i2,j2,k2,titr_2
 60       format('The structures of the OUTPUT FILES are different:',
     .            3i6,a32,' .NE. ',3i6,a32)
          stop 
      END IF


C +--OUTPUT
C +  ======

      IF (titr_1 .NE. titbak)                                       THEN
          write(3,61) titbak,i__max,j__max,k__max,v1_max,v2_max,difmax
          write(6,61) titbak,i__max,j__max,k__max,v1_max,v2_max,difmax
 61       format('Difference for ',a32,' at ',3i6,
     .           ' is |',e15.6,'  -  ',e15.6,'|  =  ',e15.6)
          difmax =         0.
          v1_max =         0.
          v2_max =         0.
          i__max =         0
          j__max =         0
          k__max =         0
          titbak =    titr_1
      END IF
      IF (difmax.LT.abs(var1-var2))                                 THEN
          difmax =  abs(var1-var2)
          i__max =        i1
          j__max =        j1
          k__max =        k1
          v1_max =      var1
          v2_max =      var2
      END IF

      GO TO 1000
 1001 CONTINUE


      close(unit=1)
      close(unit=2)
      close(unit=3)


      stop
      end
