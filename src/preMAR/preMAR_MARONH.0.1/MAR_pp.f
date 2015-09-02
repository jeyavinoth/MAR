C +------------------------------------------------------------------------+
C |                                                                    MAR |
C |   MAR_pp: Preprocessor of                                              |
C |   ^^^^^^^ Modifications to be performed: see c ###                     |
C |                                                                        |
C |   Version 15 August    2004                                            |
C |   ^^^^^^^^^^^^^^^^^^^^^^^^^                                            |
C |                                                                        |
C |   OPTIONS: #IO Additional Input from stdin                             |
C |   ^^^^^^^^ #DL Nb of lines "ntl" is specified externally               |
C |            #DB Intrinsic functions modification                        |
C |                in case of double precision                             |
C |            #CJ "Court-Jus"  pour MARSIENS avertis                      |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
C +
      program MAR_pp
C +
C +
C +--Declaration of the Variables
C +  ============================
C +
      include  'MAR_pp.inc'
C +
      parameter           (nlbmar=300,nnsplu=nnppar+30)
      character*  3 labmar(nlbmar)
      character*  3 speval(nnppar)
      character* 76 spever
      character* 76 speMAR(nnsplu)
      character* 57 commod(nnppar)
C +
      character*  3 underl
      character*  3 labtot(600)
      character*  5 form2
c #DB character*  8 chafun(600)
      character* 11 statusMAR
      character* 12 fout 
      character* 15 snam(600)
      character* 15 fnam(600)
      character* 15 unam(600)
      character* 15 saux ,uaux ,vaux
      character* 80 line ,lin1 ,buffer(400)
C +
      common/dialab/labtot
      common/diaiii/iiitot
C +
      dimension     inam(600)
      dimension     jnam(600)
      dimension     lnam(600)
      dimension     mnam(600)
      dimension     nnam(600)
      dimension     imas(600)
      dimension     jmas(600)
      dimension     mask(100000)
c #DB dimension     linfun(600)
C +
C +
C +--DATA
C +  ====
C +
      logical public,double,lapack
      logical micphy,polmod,snomod,blomod,extfrc
      data    public/.false./
      data           double/.false./
      data                  lapack/.false./
      data    micphy/.false./
      data           polmod/.false./
      data                  snomod/.false./
      data                         blomod/.false./
      data                                extfrc/.false./
C +
      data ntl/80000/
C +..      ntl:nombre total de lignes
C +
C +
C +--INITIAL Values
C +  ==============
C +
      isf = 0
      ift = 0
      isu = 0
      lev = 0
      nfu = 0
C +
      ns  = 0
C +
      ilbmar = 0
      logspe = 0
C +...logspe = 1 or 0: decide if MAR_pp.inp is User Friendly or not
C +
      no_wri = 1
C +
      nnspev = 0
      logMAR = 1
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--IO START 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
C +--Standard Control Input
C +  ======================
C +
      open (unit=3,status='unknown',file='MAR_pp.inp')
      rewind    (3)
C +
                                            nnpspe = 1
 5000 continue
        read    (3,501,end=5001)            spever
 501    format(a76)
C +
C +
C +--MAR_pp.inp is     User Friendly
C +  -------------------------------
C +
        IF   (spever(1:10).eq.'+---------')
     .        logspe = 1
        IF   (spever(1:10).eq.'+-T-------')                         THEN
              logspe = 1
              public =.true.
        END IF
        IF   (logspe      .eq. 1)                                   THEN
          IF (spever(2: 2).eq.'X')                                  THEN
              speval(nnpspe)=               spever( 8:10)
              commod(nnpspe)=               spever(16:72)
C +
C +
C +--MAR_pp.inp is NOT User Friendly
C +  -------------------------------
C +
          END IF
        ELSE
              speval(nnpspe)=               spever( 1: 3)
        END IF
C +
C +
C +--Special Attention to External Forcing
C +  -------------------------------------
C +
        IF   (speval(nnpspe).eq.'#SB')      extfrc = .true. 
C +
C +
C +--Special Attention to Double Precision
C +  -------------------------------------
C +
        IF   (speval(nnpspe).eq.'#DP')      double = .true. 
C +
C +
C +--Special Attention to lapack Library
C +  -----------------------------------
C +
        IF   (speval(nnpspe).eq.'#LA')      lapack = .true. 
C +
                     nnpspe =               nnpspe + 1
C +
C +
C +--Write MAR_pp.inp status
C +  -----------------------
C +
        IF   (nnpspe      .eq. 1 .AND.
     .        no_wri      .eq. 1 )                                  THEN
              no_wri       =   0
          IF                               (public)                 THEN
            statusMAR = 'the PUBLIC '
          ELSE
            statusMAR = ' a RUNABLE '
          END IF
C +
            write   (6,506)                 statusMAR
 506        format(/,' MAR  PREPROCESSING: GENERATION of ',a11,
     .               ' VERSION')
 555        continue
C +
        END IF
C +
C +
C +--Labels Comments for preprocessed MAR code
C +  -----------------------------------------
C +
        IF     (logspe       .eq.          1 .AND.
     .          logMAR       .eq.          1 )                      THEN
          IF   (spever( 1:10).ne.'         1')                      THEN
c #CC       IF (spever( 2: 2).eq.'X'         .OR.  .NOT.public)     THEN
                nnspev        =  nnspev  + 1
                spever( 1: 3) = 'C |'
                spever(76:76) =   '|'
                speMAR(nnspev)=  spever
c #CC       END IF
          ELSE
                logMAR        =  0
          END IF
c #VR     write(6,5995) nnspev,speMAR(nnspev)
 5995     format(i6,3x,a76)
        END IF
C +
C +
C +--END of MAR_pp.inp
C +  -----------------
C +
      goto 5000
 5001 continue
C +
      close(unit=3)
C +
              nnpspe        =               nnpspe - 1
C +
C +
C +--Labels Comments Assignation
C +  ---------------------------
C +
      IF     (logspe      .eq. 0)                                   THEN
              spever(1:3) =          speval(1)
        IF   (spever(1:1).ne.'#'.and.spever(1:1).ne.'_')            THEN
              nnpspe = 0
        ELSE 
          DO  inpspe = 1,nnpspe
              inplab = 1
 5010     continue
          IF (       inplab .gt.nnppar        )               go to 5012
          IF (speval(inpspe).eq.labval(inplab))               go to 5011
              inplab =                 inplab + 1
          go to 5010
 5011     continue
              commod(inpspe)  = comval(inplab)
 5012     continue
          END DO
        END IF
      END IF
C +
      IF (nnpspe.gt.0)                                              THEN
         write(6,5006)(speval(inpspe),commod(inpspe),inpspe = 1,nnpspe)
 5006    format(/,' STANDARD TREATMENT:',/,(1x,a3,': ',a57))
      END IF 
C +
C +
C +--Special  Control Input
C +  ======================
C +
c #IO write(6, 1)
 1    format(/,' MAR (Modele Atmospherique Regional) PreProcessing',
     .       /,' +++++++++++++++++++++++++++++++++++++++++++++++++')
c #IO write(6,11)
 11   format(' Include Cloud Microphysics (F,T): ',$)
c #IO read (5,115) micphy
c #IO write(6,115) micphy
 115  format(l1)
C +
c #IO write(6,12)
 12   format(' Include Polynya      Model (F,T): ',$)
c #IO read (5,115) polmod
c #IO write(6,115) polmod
C +
c #IO write(6,13)
 13   format(' Include Snow         Model (F,T): ',$)
c #IO read (5,115) snomod
c #IO write(6,115) snomod
C +
c #IO write(6,14)
 14   format(' Include Blowing Snow Model (F,T): ',$)
c #IO read (5,115) blomod
c #IO write(6,115) blomod
C +
c #IO write(6,15)
 15   format(' DOUBLE  Precision          (F,T): ',$)
c #IO read (5,115) double
c #IO write(6,115) double
C +
                   low = 0
c #IO write(6,16)
 16   format(' Full Code OUTPUT only      (0,1): ',$)
c #IO read (5,165) low
c #IO write(6,165) low
 165  format(i1)
C +
C +
C +--OUTPUT
C +  ======
C +
      open (unit=1,status='old',    file='MAR___.FOR')
      rewind    (1)
C +
      open (unit=2,status='unknown',file='MAR___.pp1')
      rewind    (2)
C +
      open (unit=8,status='unknown',file='MAR___.out')
      rewind    (8)
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--CODE MODIFICATIONS 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +
C +   * Begin Loop on input lines; il = current line number :
      il = 0
 1000 continue
C +
c #DL do 100 il=1,ntl
C +...Do Loop specified by data
C +
      read(1,101,end=1001)line
      il = il + 1
 101  format(a80)
      IF (lapack.and.line(24:38).eq.'FUNCTION SLAMCH') il = il - 1
      IF (lapack.and.line(24:38).eq.'FUNCTION SLAMCH') go to 1001
C +                                  456789012345678
C +
C +
C +--TOO Long Lines
C +  ==============
C +
                              icm=1
      do 102 ic=1,80
      if (line(ic:ic).ne.' ') icm=ic
 102  continue
      if (icm.eq.80) write(6,103)il,line
 103  format(i6,' | ',a80)
C +
C +
C +--OUTPUT Format
C +  =============
C +
C +   ***********
      call forwri(form2,icm)
C +   ***********
C +
      if (line( 1: 3).eq.'c _'.or.
     .    line( 1: 3).eq.'c #') then
C +
C +
C +--MAR undefined LABELS Registration
C +  =================================
C +
           llbmar=       0
      do   inppar=1,nnppar
       if   (line(3:5).eq.labval(inppar)) llbmar = 1
      end do
C +
      if  (llbmar.eq.0) then
       if (ilbmar.eq.0) then
           llbmar=       1
           ilbmar=ilbmar+1
           labmar(ilbmar) = line(3:5)
       else
        do jlbmar=1,ilbmar
         if (line(3:5).eq.labmar(jlbmar)) llbmar = 1
        end do
       end if
      end if
C +
      if (llbmar.eq.0) then
           ilbmar=ilbmar+1
           labmar(ilbmar) = line(3:5)
      end if
C +
C +
C +--Standard LABELS
C +  ===============
C +
      do 502 inpspe=1,nnpspe
          call spval(line,speval(inpspe))
 502  continue
C +
C +
C +--Special  LABELS
C +  ===============
C +
        if  (micphy) call spval(line,'#HY')
        if  (blomod) call spval(line,'#BS')
        if  (polmod) call spval(line,'#PO')
        if  (snomod) call spval(line,'#SN')
C +
      end if
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--MAR       SUBROUTINES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
        ic =  6
        los=  0
        kos=  0
C +
 20     continue
C +
        if (line(ic:ic+9).eq.'subroutine'.or.los.eq.1) then
         lev = 1
         if(los.eq.0) then 
c #WR     write(6,form2)line(1:icm)
          isf= isf+1
          ic = ic +10
          los= 1
         end if
         if(line(ic:ic  ).ne.' ')then
           if(kos.eq.0) ic1= ic
           kos= 1
 210      continue
           ic = ic+1
          if(line(ic:ic  ).ne.' '.and.
     .       line(ic:ic  ).ne.'('.and.
     .            ic      .ne.icm+1)     go to 210
           ic2= ic-1
           ics= ic-ic1
           saux(1:ics) = line(ic1:ic2)
           snam(isf)   = saux
           inam(isf)   = ics
           lnam(isf)   = il
           imas(isf)   = 0
           jmas(isf)   = 1
c #WR      write(6,212)  il,isf,ics,ic1,ic2,icm
 212       format(i6,' subroutine no',i3,3x,i3,'=',i3,'->',i3,' | =',i3)
c #WR      call forwri(form1,ics)
cc         write(6,*) saux(1:ics)
           if (saux(1:7).eq.'rel_rea') ilf=il
C +...     Debut de la derniere Routine
           los= 0
         end if
        end if
C +
        if (ic.ge.70) go to 21
         ic = ic + 1
        go to 20
 21     continue
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--MAR       FUNCTIONS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
        ic =  6
        lof=  0
        kof=  0
C +
 30     continue
C +
        if ((line(1:3).ne.'C |'.and.line(ic:ic+7).eq.'function').or.
     .       lof.eq.1) then
         lev = 1
         if(lof.eq.0) then 
c #WR      write(6,form2)line(1:icm)
           isf= isf+1
           ift= ift+1
           ic = ic +8
           lof= 1
         end if
         if(line(ic:ic  ).ne.' ')then
           if(kof.eq.0) ic1= ic
           kof= 1
 310      continue
          ic = ic+1
          if(line(ic:ic  ).ne.' '.and.
     .       line(ic:ic  ).ne.'('.and.
     .            ic      .ne.icm+1)     go to 310
           ic2= ic-1
           ics= ic-ic1
           saux(1:ics) = line(ic1:ic2)
           snam(isf)   = saux
           fnam(ift)   = saux
           inam(isf)   = ics
           lnam(isf)   = il
           imas(isf)   = 1
           jmas(isf)   = 1
c #WR      write(6,312)  il,isf,ics,ic1,ic2,icm
 312       format(i6,' function   no',i3,3x,i3,'=',i3,'->',i3,' | =',i3)
c #WR      call forwri(form1,ics)
c #WR      write(6,form1)saux(1:ics)
           lof= 0
         end if
        end if
C +
        if (ic.ge.72) go to 31
         ic = ic + 1
        go to 30
 31     continue
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--ACTIVATED SUBROUTINES
C +--ACTIVATED FUNCTIONS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
       if (lev.eq.0) then 
           lflag=1
       else
           lflag=0
       end if
C +
       if (line( 1: 3).ne.'c _'.and.
     .     line( 1: 3).ne.'c #'.and.
     .     lflag      .eq.   1) then
C +
C +        **********
           call utile(line,unam,form2,il,icm,isu,jnam,mnam)
C +        **********
C +
       ELSE
C +
C +
C +--Generation  of  
C +             the list of non used routines
C +  ----------------------------------------
C +
                                          lc = 0
           DO ic=7,66
             IF (line(ic:ic+3).eq.'call') lc = 1
           END DO
           IF(lc.eq.1)                                              THEN
                      is  =       0
 200         CONTINUE
                      is  =  is + 1
             IF (     is .gt.ns)                               GO TO 202
             IF (nnam(is).gt.il)                               GO TO 201
             GO TO 200
 201         CONTINUE
               DO js = ns,is,-1
                 nnam(js+1) = nnam(js)
               END DO
                 nnam(is)   =      il
                      ns    = ns  + 1
             GO TO 203
 202         CONTINUE
                      ns    = ns  + 1
                 nnam(ns)   = il
C +...           nnam     : call line number of a non activated subroutine
C +
 203         CONTINUE
c #VR        write(6,5998) ns,(nnam(is),is=1,ns)
 5998        format(' A',i6,10i6,/,(8x,10i6))
           END IF
C +
       END If
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--DOUBLE PRECISION : 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
      if (double) then
C +
C +       **********
          call spval(line,'#DP')
C +       **********
C +
      end if
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--DOUBLE PRECISION : INTRINSINC FUNCTIONS MODIFICATION 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +
C +
c #DB if (double) then
c #DB   icha = 6
 1080   continue
c #DB   icha = icha + 1
C +
c #DB   if (icha+5.gt.icm) go to 1085
c #DB   if (line(icha:icha+5).eq.'alog10') then
c #DB    nfu  = nfu  + 1
c #DB    linfun(nfu) = il
c #DB    chafun(nfu) = line(icha:icha+5)
c #DB    line(icha  :icha+5) =   'dlog10'
c #DB   end if
 1085   continue
C +
c #DB   if (icha+4.gt.icm)                go to 1084
c #DB   if (line(icha:icha+4).eq.'amax1') then
c #DB    nfu  = nfu  + 1
c #DB    linfun(nfu) = il
c #DB    chafun(nfu) = line(icha:icha+4)
c #DB    line(icha  :icha+4) =   'dmax1'
c #DB   end if
C +
c #DB   if (line(icha:icha+4).eq.'amin1') then
c #DB    nfu  = nfu  + 1
c #DB    linfun(nfu) = il
c #DB    chafun(nfu) = line(icha:icha+4)
c #DB    line(icha  :icha+4) =   'dmin1'
c #DB   end if
 1084   continue
C +
c #DB   if (icha+3.gt.icm)               go to 1083
c #DB   if (line(icha:icha+3).eq.'ifix') then
c #wr    write(6,1082)icm,line
 1082    format(i4,a80)
c #DB    nfu  = nfu  + 1
c #DB    linfun(nfu) = il
c #DB    chafun(nfu) = line(icha:icha+3)
c #DB    ichamx =    icm +2
c #DB    lind   =   ' '
c #DB    lind(     1:icha-1) =    line(     1:icha-1)
c #DB    lind(icha  :icha+5) =   'jidint'
c #DB    lind(icha+6:ichamx) =    line(icha+4:icm)
c #DB    line                =    lind
c #DB    icm    =    ichamx
c #wr    write(6,1082)icm,line
c #DB    call forwri(form2,icm)
c #DB   end if
C +
c #DB   if (line(icha:icha+3).eq.'alog') then
c #DB    nfu  = nfu  + 1
c #DB    linfun(nfu) = il
c #DB    chafun(nfu) = line(icha:icha+3)
c #DB    line(icha  :icha+3) =   'dlog'
c #DB   end if
 1083   continue
C +
c #DB   if (icha.lt.icm) go to 1080
C +
c #DB end if
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--FIRST OUTPUT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
      write(2,form2)line( 1:icm)
 100  continue
C +...For Do Loop, when ntl specified by data
C +
C +   * End of the loop on input lines.
      goto 1000
 1001 continue
C +   * Total number of lines found in input file = current line No :
      ntl = il
C +
      close(1)
      close(2)
C +
      if (ilbmar.gt.0) then
        write(8,1086)
 1086   format(/,' Labels MAR repris sans commentaires',
     .         /,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',
     .         /,' Numero | Label   ',
     .         /,' -------+---------')
        write(8,1087)(jlbmar,labmar(jlbmar),jlbmar=1,ilbmar)
 1087   format(i7,' | ',a8)
      end if
C +
c #DB if (double) then
c #DB   write(8,1088)
 1088   format(/,' Fonctions Intrinseques Modifiees (Double Precision)',
     .         /,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',
     .         /,' Numero | Fonction',
     .         /,' -------+---------')
c #DB   write(8,1089)(linfun(ifu),chafun(ifu),ifu=1,nfu)
 1089   format(i7,' | ',a8)
c #DB end if
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--SECOND LEVEL ACTIVATED SUBROUTINES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
      write(6,45)
 45   format(/,' PREPROCESSING of the SubRoutines:')
 41   continue
      isud = isu - isuo
      write(6,46)isud
 46   format(16x,'!',i4,' calls')
      isuo = isu
C +
      DO is = 1,isf
      IF (jmas(is).eq.1)                                            THEN
                                          lflag = 0
          saux = snam(is)
C +
        DO iv = 1,isu
          uaux = unam(iv)
          maux = jnam(iv)
          IF (saux(1:maux).eq.uaux(1:maux).and.
     .        inam(is)    .eq.jnam(iv)    ) lflag = 1
        END DO
C +
        IF (lflag      .eq.   1)                                    THEN
          jmas(is)=0
C +...    Au passage suivant, on ne cherchera plus 
C +       une sous-routine appelee dans la routine "is"
C +
          il1=lnam(is)
          il2=lnam(is+1)-1
C +...    1e et derniere ligne d'une routine particuliere is
C +
          open(unit=1,status='old',    file='MAR___.pp1')
          rewind(1)
C +
          DO  il=1,il1-1
            read(1,101)line
          END DO
          DO  il=il1,il2
            read(1,101)line
            IF (line( 1: 3).ne.'c _'.and.
     .          line( 1: 3).ne.'c #')                               THEN
C +
C +           **********
              call utile(line,unam,form2,il,icm,isu,jnam,mnam)
C +           **********
C +
C +--Elimination of effectively used routines 
C +        from the list of non used routines
C +  ----------------------------------------
C +
              js = 0
              ls = 0
 420          CONTINUE
              js = js + 1
                IF (nnam(js) .eq.mnam(isu))                         THEN
                  DO     ks = js+1,ns
                    nnam(ks)=    nnam(ks+1)
                  END DO
                    nnam(ns)=    0
                         ns = ns-1
                         ls =    1
                END IF
              IF (js.ge.ns   .or.ls.eq.1  )                    GO TO 421
              GO TO 420
 421          CONTINUE
            ELSE
C +
C +
C +--Generation  of  
C +             the list of non used routines
C +  ----------------------------------------
C +
                                             lc = 0
              DO ic=7,66
                IF (line(ic:ic+3).eq.'call') lc = 1
              END DO
              IF(lc.eq.1)                                           THEN
                         ms  =       0
 430            CONTINUE
                         ms  =  ms + 1
                IF (     ms .gt.ns)                            GO TO 432
                IF (nnam(ms).gt.il)                            GO TO 431
                GO TO 430
 431            CONTINUE
                  DO js = ns,ms,-1
                    nnam(js+1) = nnam(js)
                  END DO
                    nnam(ms)   =      il
                         ns    = ns  + 1
                GO TO 433
 432            CONTINUE
                         ns    = ns  + 1
                    nnam(ns)   = il
C +...              nnam     : call line number 
C +                            of a non activated subroutine
C +
 433            CONTINUE
c #VR           write(6,5999) ns,(nnam(ks),ks=1,ns)
 5999           format(' B',i6,10i6,/,(8x,10i6))
              END IF
C +
            END IF
          END DO
C +
          close(1)
C +
        END IF
      END IF
C +
      END DO
C +
      if (isu.eq.isuo) go to 40
      go to 41
 40   continue
C +
C +
C +--PreProcessing Statistics
C +  ------------------------
C +
      write(8,1040)
 1040 format(/,' Sous-routines et fonctions presentes dans le code',
     .       /,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',
     .       /,5x,2x,'Routine',7x,'1e ligne')
      do 1041 is=1,isf
      saux = snam(is)
      vaux = saux(1:inam(is))
      write(8,1046)is,imas(is),vaux,lnam(is)
 1046 format(i3,i2,2x,a15,i10)
 1041 continue
C +
      write(8,1042)
 1042 format(/,' Sous-routines et fonctions utilisees dans le code',
     .       /,' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',
     .       /,5x,2x,'Routine',7x,'1e ligne')
      do 1043 is=1,isu
      saux = unam(is)
      vaux = saux(1:jnam(is))
      write(8,1046)is,jmas(is),vaux,mnam(is)
 1043 continue
C +
      if (low.eq.0) then
C +
      write(8,104) ntl
      write(6,104) ntl
 104  format(/,'                FIN LECTURE, OUTPUT complet : ',
     .         'MAR___.pp1',
     .       /,'                Nombre total de lignes      : ',i5,/,1x)
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--PREPROCESSING: SECOND PHASE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
      do 400 il = 1,lnam(1)-1
      mask(il)=1
C +...masque du programme principal
C +
 400  continue
C +
      do 401 iw = 1,isf-1
      lgw=imas(iw) 
C +...Une fonction est reecrite d'office (dans ce cas imas=1)
C +
      il1=lnam(iw)
      il2=lnam(iw+1)-1
C +...1e et derniere ligne d'une routine particuliere iw
C +
      saux = snam(iw)
C +
      do 402 iv = 1,isu
      uaux = unam(iv)
      maux = jnam(iv)
      if (saux(1:maux).eq.uaux(1:maux).and.
     .    inam(iw)    .eq.jnam(iv)    ) lgw = 1
 402  continue
C +
      if (lgw.eq.1) then
       saux= snam(iw)
       vaux= saux(1:inam(iw))
       write(8,4026) iw,lnam(iw)
 4026  format(' Routine no',i3,' Ligne no',i5,' ',$)
       write(8,4027) vaux
 4027  format(a15)
C +...La routine est-elle utilisee ?
C +
      end if
C +
      do 403 il = il1,il2
      mask(il)=lgw
 403  continue
C +
 401  continue
C +
      il1=lnam(isf)
      il2=ntl
      do 404 il=il1,il2
      mask(il)=1
 404  continue 
C +...La derniere routine (en principe rel_rea) est ecrite
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--SECOND OUTPUT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
c #VR write(6,6000) (nnam(is),is=1,ns)
 6000 format(15i6)
C +
      kbuf = 0
C +...kbuf = 0 ===> Line is     entered in the buffer
C +          1              NOT entered in the buffer
C +
      lbuf = 0
C +...lbuf = 0 ===> Buffer      written and cleared
C +          1                  incremented 
C +          2                  cleared
C +
      mbuf = 0
C +...mbuf = 1 ===> Buffer Line indexed
C +
      nbuf = 1
C +...nbuf:         Buffer Content
C +
      ncm  = 0
      lcr  = 0
C +
      open(unit=1,status='old',file='MAR___.pp1')
      rewind(1)
C +
      open(unit=2,status='new',file='MAR___.for')
      rewind(2)
C +
      open(unit=4,status='new',file='MAR___.rmv')
      rewind(4)
C +
C +
C +--Line number for (non) effective routine call
C +  --------------------------------------------
C +
              ilw=0
C +
                        is=1
              ils= mnam(is)
                        ns=1
              nls= nnam(ns)
      DO      il = 1 , ntl
C +
        IF   (il.gt.ils)                                            THEN
              is = 1 +  is
              ils= mnam(is)
c #VR         write(6,*)'ils ',ils,is
          IF (ils.eq.0)  ils = ntl + 1
C +...      active subroutine call line number
C +
        END IF
C +
        IF   (il.gt.nls +2)                                         THEN
C +...                  +2 to memorize nls until 2nd "***" is passed
C +
              ns = 1 +  ns
              nls= nnam(ns)
c #VR         write(6,*)'nls ',nls,ns
          IF (nls.eq.0)  nls = ntl + 10
C +...  non active subroutine call line number
C +
        END IF
C +
C +
C +--INPUT from intermediary file
C +  ----------------------------
C +
        read(1,101)line
C +
C +--Search Ends of Active Branchings
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   logq      =   0
        DO ic=7,60
              IF ( line(1:1).eq.' '           .AND. 
     .            (line(ic:ic+5).eq.'END IF'  .OR.
     .             line(ic:ic+5).eq.'end if'  .OR.
     .             line(ic:ic+7).eq.'CONTINUE'.OR.
     .             line(ic:ic+7).eq.'continue'    ))                THEN
                   logq      =   1
              END IF
        END DO
C +
C +--Buffer activation
C +  ~~~~~~~~~~~~~~~~~
        IF        (mbuf     .eq. 1            )                     THEN
                   underl    =   line(6:8)
C +...  Previous Title Level ('===','---', or '~~~')
C +
c #VR              write(6,6004) il,mbuf-1,line
        END  IF
C +
        IF        (line(1:4).eq.'C +-'   .AND.
     .             line(5:9).ne.'-----'       )                     THEN
                   lbuf      =   1 + lbuf
                   mbuf      =   1
c #VR              write(6,6004) il,mbuf  ,line
 6004              format(2i6,3x,a80)
        ELSE IF   (line(1:3).eq.'###'     .OR.
     .             line(1:3).eq.'C |'     .OR.
     .             line(1:5).eq.'C +++'   .OR.
     .             line(1:6).eq.'      '  .OR.
     .             logq     .eq. 1        .OR.
     .             il       .eq. ils - 1      )                     THEN
C +...             effective routine located on line il.eq.ils
C +
                   kbuf      =   0
                   mbuf      =   0
C +
C +--Decide to clear the buffer
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~
                   lbuf      =   0
          IF      (logq     .eq. 1        .AND.
     .             lbuf     .ge. 1            )
     .             lbuf      =   2
C +
        ELSE
                   mbuf      =   0
        END  IF
C +
C +--Decide to delay the clearing of the buffer
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF        (mbuf     .eq. 0       .AND.
     .             lbuf     .eq. 2            )                     THEN
          IF      (line(6:8).ne. underl       )
     .             lbuf = lbuf - 1
c #VR     write(6,6005)il,lbuf,buffer(2),underl,line
 6005     format(2i6,3x,a80,/,9x,a3,3x,a80)
        END IF
C +
C +--Decide to stop  the loading  of the buffer (non used routines)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF       ( line(1:3).eq.'c #'     .OR.
     .            (il       .ge. nls - 1 .AND.
     .             il       .le. nls + 2      ))                    THEN
                   kbuf      =   1
c #VR              write(6,6002)il,nls,line
 6002              format(2i6,a80)
        ELSE
                   kbuf      =   0
        END IF
C +
C +--Decide to stop  the loading  of the buffer (specific comments)
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF        (line(1:5).eq.'C | #')                            THEN
                   kbuf = 1
c #VR              write(6,6003) il,line
        ELSE IF   (line(1:3).eq.'C |')                              THEN
                                           kcom = 0
          DO ic=5,60
              IF  (line(ic:ic).eq.'#'.OR.
     .             line(ic:ic).eq.'_'    )                          THEN
                DO nlab = 1,nnppar
                  IF (line(ic:ic+2).eq.labval(nlab))
     .                                     kcom = 1
                END DO
              END IF
          END DO
              IF  (kcom         .eq. 1     )                        THEN
c #VR              write(6,6003) il,line
 6003              format(' Comment',i6,' :  ',a80)
                   kbuf = 1
              END IF
        END IF
C +
C +--"Court-Jus"
C +  ~~~~~~~~~~~
        IF        (.NOT.public)                                     THEN
                   kbuf      =   0
                   lbuf      =   0
                   mbuf      =   0
        END IF
C +
C +--Final preprocessing
C +  ~~~~~~~~~~~~~~~~~~~
        IF        (line(1:3).eq.'###')                              THEN
          IF(public)                                                THEN
                   line(1:3) =  'c #'
          ELSE
                   line(1:5) =  '     '
C +...    Here preprocessed Labels in spval are removed
C +
          END IF
        END IF
C +
C +
C +--DOUBLE PRECISION : MODIFICATION DES DECLARATIONS REAL       
C +  -----------------------------------------------------
C +
        if (double) then
C +
C +       **********
          call spdou(line,il)
C +       **********
C +
        end if
C +
C +
C +--OUTPUT
C +  ------
C +
C +--Include Comments on Preprocessing options
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF    (public)                                              THEN
          IF  (line(10:27).eq.'MAIN       OPTIONS')                 THEN
            DO inspev =          1,nnspev-1
               write(2,501) speMAR(inspev)
            END DO
          END IF
        END IF
C +
C +--Buffer Update
C +  ~~~~~~~~~~~~~
        IF  (kbuf.eq.0)                                             THEN
             buffer(1)    = line
C +
C +--Clear the buffer
C +  ~~~~~~~~~~~~~~~~
          IF         (lbuf.eq.2.AND.mbuf.eq.0  )                    THEN
            IF       (line(1:6).eq.'      '.OR.
     .                logq     .eq. 1          )                    THEN
              IF     (nbuf     .ge. 2          )                    THEN
              DO ibuf=nbuf,2,-1
               write(4,480)   buffer(ibuf)
              END DO
              END IF
               lbuf=0
               nbuf=1
            ELSE
              IF     (nbuf     .ge. 3        )                      THEN
              DO ibuf=nbuf,3,-1
               write(4,480)   buffer(ibuf)
              END DO
              END IF
               lbuf=1
               nbuf=2
            END IF
          END IF
C +
C +--Load  the buffer and prepare next line
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF         (lbuf.ge.1)                                    THEN
               nbuf=1+nbuf
            IF(nbuf.gt.400)   write(6,6001)il,line
 6001       format(' WARNING: Buffer is full, Line',i6,3x,a80)
            DO ibuf=  nbuf,2,-1
               buffer(ibuf) = buffer(ibuf-1)
            END DO
C +
C +--Write down and clear the buffer
C +  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ELSE IF    (lbuf.eq.0)                                    THEN
               lin1         = buffer(1)
            DO ibuf=  nbuf,1,-1
               line         = buffer(ibuf)
C +
                                         icm=1
              DO         ic=1,80
                IF (line(ic:ic).ne.' ')  icm=ic
              END DO
              IF (icm.eq.80) write(6,103)il,line
C +...        Too Long Lines
C +
C +             ***********
                call forwri(form2,icm)
C +             ***********
C +
              IF  (icm.eq.3 .and. line(1:icm).eq.'C +')             THEN
                   ncm =  1 + ncm
              ELSE
                   ncm =  0
              END IF
C +
c #CO         IF  (lin1(1:5).eq.'C +++'.AND.line(1:5).ne.'C +++')   THEN
              IF  (lin1(1:5).eq.'C +++'.AND.line(1:5).ne.'C +++'
     .                                 .AND.line(1:1).ne.' ' 
     .                                 .AND.nbuf     .ne. 1     )   THEN
                   lcr =  1
              ELSE
                   lcr =  0
              END IF
C +
              IF  (mask(il).eq.1 .and. ncm.le.2 .and. lcr .eq. 0)   THEN
                   write(2,form2) line(1:icm)
                   ilw =  1 + ilw
              END IF
C +
            END DO
C +
               nbuf=  1
C +...         nbuf=  1: reset the buffer
C +
          ELSE
               write(4,form2) line(1:icm)
          END IF
C +
        ELSE
               write(4,480)   line
 480           format(a80)
        END IF
C +
      END DO
C +
      close(1)
      close(2)
      close(3)
C +
      else
      ilw = ntl
C +
      end if
C +
                  fout     =   'MAR___.for  '
      write(6,500)fout,ilw
      write(8,500)fout,ilw
 500  format(/,'                Fichier OUTPUT ',a12,' : ',
     .                                          i5,' lignes')
C +
      write(8,88)(labtot(iii),iii=1,iiitot)
 88   format(/,' Labels Enleves :',
     .       /,' ================',20(/,10(3x,a3)))
      close(8)
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--MAR_xy.inc Files preprocessing, when double precision is required
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + 
C +
      IF    (double)                                              THEN ! CTR
C +
C +     ***********
        CALL fildou('MARCTR.inc')
        CALL fildou('MARphy.inc')
        CALL fildou('MARgrd.inc')
        CALL fildou('MAR_GE.inc')
        CALL fildou('MARSND.inc')
        CALL fildou('MAR_RA.inc')
        CALL fildou('MAR_DY.inc')
        CALL fildou('MAR_NH.inc')
        CALL fildou('MAR_CU.inc')
        CALL fildou('MAR_LB.inc')
        CALL fildou('MAR_DI.inc')
        CALL fildou('MAR_FI.inc')
        CALL fildou('MAR_HY.inc')
        CALL fildou('MAR_EW.inc')
        CALL fildou('MAR_CA.inc')
        CALL fildou('MAR_BS.inc')
        CALL fildou('MAR_SN.inc')
        CALL fildou('MARdCP.inc')
        CALL fildou('MAR_IB.inc')
        CALL fildou('MAR_TC.inc')
        CALL fildou('MAR_SL.inc')
        CALL fildou('MAR_PO.inc')
        CALL fildou('MAR_TV.inc')
        CALL fildou('MAR_VB.inc')
        CALL fildou('MAR0SV.inc')
        CALL fildou('MARdSV.inc')
        CALL fildou('MARxSV.inc')
        CALL fildou('MAR_WK.inc')
        CALL fildou('MAR_IO.inc')
C +     ***********
C +
C +     ***********
        CALL fildou('MAR_IR.inc')
        CALL fildou('MAR_SO.inc')
        CALL fildou('MARpen.inc')
        CALL fildou('MAR_BR.inc')
        CALL fildou('MAR_TU.inc')
        CALL fildou('MAR_TE.inc')
        CALL fildou('MAR_OL.inc')
        CALL fildou('MAR_PV.inc')
C +     ***********
C +
C +     ***********
        CALL fildou('MARvec.inc')
C +     ***********
C +
C +     ***********
        CALL fildou('MAR_2D.inc')
        CALL fildou('MARc2D.inc')
C +     ***********
C +
        IF (extfrc)                                               THEN ! CTR
          write(6,*) '*SBCnew.for is pre-processed********************'
C +
C +       ***********
          CALL fildou('SBCnew.for')
C +       ***********
C +
        END IF                                                         ! CTR
      END IF                                                           ! CTR
C +
C +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +--NEW MAR_pp.inc MAR_pp_dat.f 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C +
C +
      IF (public)                                                   THEN
          open( unit=11,status='new',file='MAR_pp.PUB.inc')
          rewind     11
            write(11,1101)               nnpspe
 1101       format(6x,'parameter           (nnppar=',i3,')'
     .           /,6x,'character*  3 labval(nnppar)'
     .           /,6x,'character* 57 comval(nnppar)'
     .           /,6x,'common/Option/labval,comval')
          close(unit=11)
C +
          open( unit=12,status='new',file='MAR_pp.PUB.dat')
          rewind     12
            write(12,1200)
 1200       format('      block data MAR_pp_dat'
     .           /,'C +'
     .           /,'C +------------------------------------------',
     .             '------------------------------+'
     .           /,'C |                                          ',
     .             '                          MAR |'
     .           /,'C |   Block Data MAR_pp_dat contains the list',
     .             ' of MAR Options               |'
     .           /,'C |                                          ',
     .             '                              |'
     .           /,'C +------------------------------------------',
     .             '------------------------------+'
     .           /,'C +'
     .           /,'      include "MAR_pp.inc"'
     .           /,'C +'
     .           /,'      data (labval(inpspe),comval(inpspe),',
     .             'inpspe= 1,nnppar)/'
     .           /,'C +          123456789012345678901234567890',
     .                          '123456789012345678901234567890')
            write(12,1201)
     .           (speval(inpspe),commod(inpspe),inpspe = 1,nnpspe-1)
 1201       format(5x,'.','"',a3,'","',a57,'"',',')
            write(12,1202)
     .            speval(nnpspe),commod(nnpspe)
 1202       format(5x,'.','"',a3,'","',a57,'"/',
     .           /,'C +',
     .           /,5x,' end')
C +
          close(unit=12)
      END IF
C +
      stop
      end
      subroutine spval(line,lab)
C +
C +------------------------------------------------------------------------+
C |   subroutine spval activates option lab (i.e., "#XY")                  |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
      character* 3 lab,labtot(600)
      character*80 line
C +
      common/dialab/labtot
      common/diaiii/iiitot
C +
      if (line(3:5).eq.lab) then
          line(1:3)='###'
      end if
C +
      if (iiitot.eq.0) then
        iiitot         = 1
        labtot(iiitot) = lab
      else
                               iiiold = 0
       do 10 iii=1,iiitot
       if (lab.eq.labtot(iii)) iiiold = 1
 10    continue
       if (iiiold.eq.0) then
        iiitot         = iiitot + 1
        labtot(iiitot) = lab
       end if
      end if
C +
      return
      end
      subroutine spdou(line,il)
C +
C +------------------------------------------------------------------------+
C |   subroutine spdou modifies statements  real into real*8               |
C |                                         Real into Real*8               |
C |                                                                        |
C +------------------------------------------------------------------------+
C +
      implicit     none
      character*80 line
      integer      il,ic
C +
      IF     (line(1:6).eq.'      ')                              THEN ! CTR
        ic = 6
 1000   CONTINUE
        ic = ic+1
        IF   (line(ic:ic  ).eq.' '        )                 GO TO 1000
        IF   (line(ic:ic+3).ne.'real'.AND.
     .        line(ic:ic+3).ne.'Real')                      GO TO 1020
c #WR     write(6,6000) line
 6000     format(a80)
          IF (line(ic+4:ic+6).ne.'   ')                           THEN
            write(6,1060) il,line
 1060       format(' line',i6,' is not correctly programmed',
     .                        ' for preprocessing double precision',
     .           /,  12x, a80)
          ELSE
              line(ic+4:ic+5) = '*8'
          END IF
c #WR     write(6,6000) line
 1020   CONTINUE
      END IF                                                           ! CTR
C +
      return
      end
      subroutine forwri(form2,icm)
      character*  1 ch(0:9)
      character*  5 form2
      data ch/'0','1','2','3','4','5','6','7','8','9'/
C +
      form2(1:2)='(a'
      if (icm.lt.10) then
       form2(3:3)=ch(icm)
       form2(4:5)=') '
      else
       ic1 =     icm/10
       ic0 = mod(icm,10)
       form2(3:3)=ch(ic1)
       form2(4:4)=ch(ic0)
       form2(5:5)=')'
      end if
C +
      return
      end
      subroutine utile(line,unam,form2,il,icm,isu,jnam,mnam)
C +
C +------------------------------------------------------------------------+
c |--SOUS-ROUTINES UTILISEES DANS LE CODE                                  |
C +------------------------------------------------------------------------+
C +
      character*  5 form2
      character* 15 unam(600)
      character* 15 saux 
      character* 80 line
C +
      dimension     jnam(600)
      dimension     mnam(600)
C + 
         ic =  6
         loc=  0
         koc=  0
C +
 40     continue
C +
        if (line(ic:ic+3).eq.'call'.or.loc.eq.1) then
         if(loc.eq.0) then 
c #WR      write(6,form2)line(1:icm)
           isu= isu+1
           ic = ic +4
           loc= 1
         end if
         if(line(ic:ic  ).ne.' ')then
          if(koc.eq.0) ic1= ic
           koc= 1
 410      continue
          ic = ic+1
          if(line(ic:ic  ).ne.' '.and.
     .       line(ic:ic  ).ne.'('.and.
     .            ic      .ne.icm+1)     go to 410
           ic2= ic-1
           ics= ic-ic1
           saux(1:ics) = line(ic1:ic2)
           unam(isu)   = saux
           jnam(isu)   = ics
           mnam(isu)   = il
c #WR      write(6,412)  il,isu,ics,ic1,ic2,icm
 412       format(i6,' call       no',i3,3x,i3,'=',i3,'->',i3,' | =',i3)
c #WR      call forwri(form1,ics)
c #WR      write(6,form1)saux(1:ics)
           loc= 0
         end if
        end if
C +
        if (ic.ge.72) go to 41
         ic = ic + 1
        go to 40
 41     continue
C +
       return
       end
      subroutine fildou(incfil)
C +
C +------------------------------------------------------------------------+
C |   subroutine fildou preprocesses double precision in file MAR_xy.inc   |
C |                                                                        |
C |   INPUT: incfil: label of the MAR_xy.inc file                          |
C +------------------------------------------------------------------------+
C +
      implicit     none
C +
      character*10 incfil
      character*80 inclin(2000),auxlin
      integer      il,jl
C +
      open( unit=21,status='old'    ,file=incfil)
      rewind     21
            il =      0
 100        CONTINUE
            il = il + 1
            read(21,101,end=102)          inclin(il)
 101        format(a80)
            GO TO 100
 102        CONTINUE
      close(unit=21)
C +
      open( unit=22,status='unknown',file=incfil)
      rewind     22
            DO   jl = 1,il-1
              auxlin                     =inclin(jl)
C +
C +           **********
              call spdou(auxlin,jl)
C +           **********
C +
              write(22,101)               auxlin
            END DO
      close(unit=22)
C +
      return
      end
