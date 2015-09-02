C   +-------------------------------------------------------------------+
C   |  Subroutine UPScor                             08/2004    NESTING |
C   +-------------------------------------------------------------------+
C   |                                                                   |
C   |  "Upscaling" from 250m CORINE land-cover data to MAR              |
C   |                                                                   |
C   |  Output :                                                         |
C   |     CORfrc(i,j,class  Fraction of MAR mesh covered by type "class"|
C   |                                                                   |
C   |  Method :                                                         |
C   |     MAR mesh is approximated by a quadrilateral in the input grid,|
C   |     all points in this quadrilateral are counted in the fraction  |
C   |     Note: this is a fairly general upscaling technique;           |
C   |     some modest changes may provide accurate upscaling for other  |
C   |     fields.                                                       |
C   +-------------------------------------------------------------------+
      SUBROUTINE UPScor (ainX, ainY, imx, imy, nclass,
     .                   in_FID, inVNAM, out_X, out_Y, CORfrc) 

      IMPLICIT NONE


C +---General variables
C +   -----------------

      INCLUDE 'NSTdim.inc'
      INCLUDE 'NetCDF.inc'

      INTEGER  VARSIZE
      EXTERNAL VARSIZE

C +---Arguments
C +   ---------------
      INTEGER imx,imy, nclass, in_FID
      CHARACTER(*) inVNAM
      REAL ainX(imx), ainY(imy)
      REAL out_X(0:mx,0:my), out_Y(0:mx,0:my)
      REAL CORfrc(mx,my,nclass)
      
C +---Local variables
C +   ---------------
      INTEGER ii,jj,ki,kn,kc,Ierror,nsubFR,itmp
      INTEGER iin,jin,iinL,iinH,jinL,jinH,jin1,jin2,idint
      INTEGER idxsPt(2), icoVAL, corVID
      INTEGER iP(4), jP(4)
      REAL   outPX(4), outPY(4)
      REAL   outTOP,outBOT
      REAL   outLFT,outRHT,pti(4)
      REAL*8 RXin, RYin
      REAL   linY
      REAL   AAA,BBB,XXX,XX1,XX2,YY1,YY2,YMIN,YMAX
      INTEGER meshFR(50)
      
      IF(nclass.GT.50) THEN
        STOP
      ENDIF
      
C     Input grid resolution
C     (this simplifies the problem a bit, but you may generalize !)
C +   ---------------------
      RXin = ainX(2)-ainX(1)
      RYin = ainY(2)-ainY(1)
      
C +   Get NetCDF variable ID
      itmp =VARSIZE(inVNAM) 
      Ierror=NF_INQ_VARID(in_FID, inVNAM(1:itmp), corVID)

C +   ===============================================================
C     MAIN LOOP : on output grid points
      DO jj=1,my
      DO ii=1,mx
    
C     Initialisations :
      DO kc=1,nclass
        CORfrc(ii,jj,kc) =0.0
        meshFR(kc)=0
      ENDDO
      nsubfr=0

C       We work in an output mesh defined as:
C            1         2
C  
C            4         3   (quadrilateral)
C
C       The position of the point are assumed to follow simple rules,
C       but a more general case may be implemented here...
C       (should search e.g. "the points which have 
C        at least two under them..." to construct the quadrilateral)
  
        iP(1)=ii-1
        iP(2)=ii
        iP(3)=ii
        iP(4)=ii-1
        
        jP(1)=jj
        jP(2)=jj
        jP(3)=jj-1
        jP(4)=jj-1
        DO ki=1,4
          outPX(ki)=out_X(iP(ki),jP(ki))
          outPY(ki)=out_Y(iP(ki),jP(ki))
        ENDDO
        
C       Find first + end line pos. indexes in input grid (jinL,jinH)
C +     ------------------------------------------------------------

        outTOP = max(outPY(1),outPY(2))
        outBOT = min(outPY(4),outPY(3))

        jin1 = 1 + FLOOR ( (outBOT - ainY(1) + 0.01D01) / RYin )
        jin2 = 1 + FLOOR ( (outTOP - ainY(1) - 0.01D01) / RYin )
        jinL = min(jin1,jin2) + 1 ! inside pts, no double counting
        jinH = max(jin1,jin2) 
C       This is because corine indexes are from N to S,
C       while corine coordinates are (appropriately) from S to N

        IF(jinH.GT.imy.OR.jinL.LT.1)THEN
           write(*,*)'UPScor: ERROR - NST dom out of CORINE grid'
        ENDIF
        
C +     --------------------------=======----------------------------
C       Loop on lines in the input grid
        DO jin= jinL,jinH
          linY = ainY(jin)
        
C +       Search the intersections of output mesh / input lines
C         (4 segments of the output quadrilateral mesh)
C         -----------------------------------------------------
          idint=0   ! index of the intersection found
          DO ki=1,4
            kn=mod(ki,4)+1
            YY1=outPY(ki)
            YY2=outPY(kn)
            XX1=outPX(ki)
            XX2=outPX(kn)
            YMIN=min(YY1,YY2)
            YMAX=max(YY1,YY2)
            IF (YMIN.LE.linY.AND.linY.LE.YMAX.AND.YY2.NE.YY1) THEN
            idint=idint+1
              IF (XX2.EQ.XX1) THEN
                pti(idint)=XX1
                write(*,*) 'EQ pts'
              ELSE 
                AAA=(YY2-YY1)/(XX2-XX1)
                BBB= YY1- AAA*XX1
                XXX=(linY - BBB) / AAA
                pti(idint)=XXX
              ENDIF
            ENDIF
          ENDDO
          outLFT= min(pti(1),pti(2))
          outRHT= max(pti(1),pti(2))
          IF(idint.NE.2)THEN
           IF(idint.LT.2)THEN
             write(*,*) 'UPScor : internal error;'
             write(*,*) '(number of line intersections: ',idint,')'
             write(*,*) (outPX(ki),outPY(ki),ki=1,4)
             write(*,*) linY,jin1,jin2,jin
             write(*,*) ainY(jin-1),ainY(jin),ainY(jin+1)
             STOP
           ELSE
             outLFT= min(outLFT,pti(3))
             outRHT= max(outRHT,pti(3))
             IF(idint.GT.3)THEN
               write(*,*) 'UPScor : WARNING - something strange'
               write(*,*) '(4 intersect. = 2 peaks ?', idint
               write(*,*) XX1,YY1,XX2,YY2
               write(*,*) linY,jin1,jin2,jin
               write(*,*) 'Please check UPScor / CORINE'
             ENDIF
           ENDIF
          ENDIF

C         Find first + end index along line in input grid (iinL,iinH)
C +       -----------------------------------------------------------
          iinL = 2 + FLOOR ( (outLFT - ainX(1)) / RXin )
          iinH = 1 + FLOOR ( (outRHT - ainX(1)) / RXin )

C         Loop on points in the input grid line
C         -------------------------------------
          DO iin= iinL,iinH
        
C         Read data and update class fractions
C         ------------------------------------  
          idxsPt(1)=iin
          idxsPt(2)=jin
          Ierror=NF_GET_VAR1_INT(in_FID, corVID, idxsPt, kc)
          IF (Ierror.NE.0)  THEN
             write(*,*) 'UPScor: CORINE reading error'
             write(*,*) '  Req point was ',iin,jin
          ENDIF
          IF (kc.GT.nclass) THEN 
             write(*,*) 'Error: CORINE nclass / file (UPScor)'
             write(*,*) 'i,j,read val,nclass: ',iin,jin,kc,nclass
             STOP
          ENDIF
          nsubFR     =nsubFR+1
          meshFR(kc) =meshFR(kc)+1
                  
          ENDDO
        ENDDO
C       Loop on lines in the input grid (END)
C +     --------------------------=======----------------------------

        DO kc=1,nclass
          CORfrc(ii,jj,kc) = FLOAT(meshFR(kc)) / FLOAT(nsubFR)
        ENDDO

      
C     MAIN LOOP : on output grid points (END)
      ENDDO
      ENDDO
C +   ===============================================================


      RETURN
      END
