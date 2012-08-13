C 13 Aug 12 - DGF - padd common block QMMM2
C 21 Apr 10 - NA,KK - allow choosing minimizer, RESTRAIN-POSITION changes
C  1 May 03 - CHC - modify TINKIN to concoct a LINKGE default
C  9 MAR 00 - CHC - FIX for parallel run
C 19 OCT 98 - CHC - Modified for parallel run
C  6 MAY 98 - JRS - READS TINKER CARDS FROM THE INPUT FILE
C                   $TINXYZ: Tinker .xyz file contents
C                   $TINKEY: Tinker .key file contents
C
c  ##################################################################
c  #                                                                #
c  #     tinkin: Read $TINXYZ deck containing TINKER .xyz file      #
c  #                                                                #
c  ##################################################################
c
c-mws-      subroutine chkqmm(ir,iw)
c-mws-c
c-mws-c
c-mws-      LOGICAL master,GOPARR,DSKWRK,MASWRK,mmonly,qmmm
c-mws-      INTEGER me,nproc,ibtyp,iptim,IR,IW
c-mws-c
c-mws-      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c-mws-      common /tinopt/ mmonly, qmmm
c-mws-      qmmm=.true.
c-mws-      mmonly=.false.
c-mws-c
c-mws-      CALL SEQREW(IR)
c-mws-      CALL FNDGRP(IR,' $LINK  ',JEOF)
c-mws-      if (jeof.eq.1) then
c-mws-         qmmm=.false.
c-mws-      end if
c-mws-      if (maswrk .and. qmmm) write(IW,*)
c-mws-     *        '---- QMMM procedure is ON ----'
c-mws-c
c-mws-      return
c-mws-      end
c
      subroutine tinkin(ir,iw)
c
c
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c    first, look for the TINKER .xyz file
c
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $TINXYZ',JEOF)
c      if (maswrk) write(IW,*) 'done check for $TINXYZ '
      if (jeof.eq.1) then
         if (maswrk) then
            write(iw,*) '    *** ERROR ***'
            write(iw,*) ' $TINXYZ CARD NOT FOUND '
            write(iw,*) ' CHECK SPACING AND SPELLING  '
          endif
          call abrt
      else
c
c  -- initialize tinker arrays
c
      call initial
c
c  -- read tinker xyz file and echo to the output file
c
      call readxyz (ir,iw)
      if (maswrk) call prtxyz (iw)
c
c  -- look for tinker key stuff
c
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $TINKEY',JEOF)
      if (jeof.eq.1) then
         if (maswrk)   write(iw,*) '    *** ERROR ***'
         if (maswrk)   write(iw,*) ' $TINKEY CARD NOT FOUND '
         if (maswrk)   write(iw,*) ' CHECK SPACING AND SPELLING  '
         call abrt
      end if
c
c  -- read tinker key file and echo to the output file
c  -- call control to set printing levels
c
C  currently, mm part is not parallelized.
        call getkey(ir,iw)
        call control
c
c  -- call mechanic, last bit in Tinker initialization sequence c
      call mechanic
c
      if (maswrk) then
        write(iw,*) '     *** TINKER INTIALIZATION COMPLETE ***'
        write(iw,*)
      endif
c
      end if
c
      return
      end
c
C
C
c  ##################################################################
c  #                                                                #
c  #     toptin: Read T(inker)opt(ions) for tcng optimization, if   #
c  #              specified in an (optional) $TINOPT deck           #
c  #                                                                #
c  ##################################################################
c
      subroutine toptin
c      subroutine toptin(CTMODE,CTMETH,GRDMIN)
c
c     Tinker newton optimization grabs options MODE, METHOD (CHAR*6)
c     and GRDMIN from interative I/O. Need to move this to GAMESS
c     input file.
c
c     - Disconnect: GAMESS NAMEIO (user friendly) reads HOLLERITHS
c     This sub uses NAMEIO to read HOLLERITH variables, and uses these
c     to define the CHARACTER variables for TINKER
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'sizes.i'
c
ckk -minimize-      PARAMETER (NNAM=3)
      PARAMETER (NNAM=8)
      DIMENSION QNAM(NNAM),KQNAM(NNAM)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL MMONLY,QMMM
      CHARACTER*6 CTMODE, CTMETH
      DOUBLE PRECISION GRDMIN
C
c --  COMMONS FROM GAMESS
c
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
cjrs
ckk -minimize- begin
      logical frzmm
      DOUBLE PRECISION NEWTN, MINIMZ, MQNOPT, MINMET
      COMMON /NWTOPT/ OPTPRG,MINMET,GRDMIN,maxhess,frzmm,CTMODE,CTMETH
ckk -minimize- end
      COMMON /TINOPT/ mparti,MMONLY,QMMM
c
      DOUBLE PRECISION ICCG,NEWTON,NONE
C
      DATA TOPTMZ /8HTOPTMZ  /
ckk -fmoimomm-      DATA QNAM /8HTMODE   ,8HTMETH   ,8HGRDMIN  /
ckk -minimize-      DATA KQNAM/5,5,3/
      DATA QNAM /8HOPTPRG  ,8HTMODE   ,8HTMETH   ,8HMINMET  ,
     *           8HGRDMIN  ,8HFRZMM   ,8HMAXHES  ,8HMODPAR  /
ckk -minimize-      DATA KQNAM/5,5,3/
      DATA KQNAM/5,5,5,5,3,0,1,1/
C
ckk -minimize- begin
C     Acceptable entries for OPTPRG
      DATA NEWTN /8HNEWTN   /,MINIMZ /8HMINIMZ  /
ckk -minimize- end
C
C     Acceptable entries for TMODE and TMETH
C
      DATA AUTO  /8HAUTO    /, NEWTON/8HNEWTON  /
      DATA TNCG  /8HTNCG    /, DTNCG /8HDTNCG   /
      DATA NONE  /8HNONE    /, DIAG  /8HDIAG    /
      DATA BLOCK /8HBLOCK   /, SSOR  /8HSSOR    /
      DATA ICCG  /8HICCG    /
ckk -minimize- begin
C     Acceptable entries for MINMET
c     DATA SDOPT /8HSD      /, FROPT /8HFR      /
c     DATA PROPT /8HPR      /, HSOPT /8HHS      /
c     DATA POWOPT/8HPOW     /, MQNOPT/8HMQN     /
      DATA MQNOPT/8HMQN     /
ckk -minimize- end
C
C     Default values
C
      TMODE=DTNCG
      TMETH=ICCG
      GRDMIN=0.0001D0
ckk -fmoimomm- begin
      frzmm=.false.
ckk -fmoimomm- end
ckk -minimize- begin
      OPTPRG=NEWTN
      MINMET=MQNOPT
ckk -minimize- end
      maxhess=0
      mparti=0
c     bit-additive
c     1 do not parallelise
C
C     ----- READ NAMELIST -$TOPTMZ -----
c
      JRET = 0
      CALL NAMEIO(IR,JRET,TOPTMZ,NNAM,QNAM,KQNAM,
ckk -fmoimomm-     *     TMODE,TMETH,GRDMIN,0,
ckk -minimize     *     TMODE,TMETH,GRDMIN,frzmm,
ckk -minimize     *     0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     *     OPTPRG,TEMODE,TMETH,MINMET,GRDMIN,
     *     frzmm,maxhess,mparti,0,  0,0,0,0,0,  0,0,0,0,0,
     *     0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     *     0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     *     0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0)
C
      NERR=0
c
      IF(JRET .EQ. 0) THEN
ckk -minimize begin
c        IF(MASWRK) WRITE(IW,900) TMODE,TMETH,GRDMIN
c      ENDIF
c
c      IF(JRET .EQ. 1) THEN
c        IF(MASWRK) WRITE(IW,901) TMODE,TMETH,GRDMIN
        if(maswrk) then
            if(optprg .eq. newtn) WRITE(IW,900) OPTPRG,TMODE,TMETH,
     *                                         GRDMIN
            if(optprg .eq. minimz) WRITE(IW,903) OPTPRG,MINMET,GRDMIN
            WRITE(IW,910) maxhess,mparti 
        endif
ckk -minimize- end
      ENDIF
C
      IF(JRET .GE. 2) THEN
         IF(MASWRK) WRITE(IW,*) 'TROUBLE INTERPRETING $TOPTMZ '
         NERR=NERR+1
         CALL ABRT
      END IF
ckk- fmoimomm- begin
      if(frzmm) then
         if(maswrk) write(iw,902)
      endif
ckk -fmoimomm- end
      if(maxhess.eq.0) then
         maxhess=10000
      else if(maxhess.eq.-1) then
c        special value to use all of the 2e integral buffer
c        other nagative values mean set maxhess=abs(maxhess) using this buffer. 
c     else
c        Use the dynamic memory in GAMESS 
      endif
c     write(6,*) 'maxhess=',maxhess
c
c
c   Now we assign values for Tinker
c
      IF(TMODE .EQ. AUTO)   CTMODE='auto  '
      IF(TMODE .EQ. NEWTON) CTMODE='newton'
      IF(TMODE .EQ. TNCG)   CTMODE='tncg  '
      IF(TMODE .EQ. DTNCG)  CTMODE='dtncg '
c
      IF(TMETH .EQ. AUTO)   CTMETH='auto  '
      IF(TMETH .EQ. NONE)   CTMETH='none  '
      IF(TMETH .EQ. DIAG)   CTMETH='diag  '
      IF(TMETH .EQ. BLOCK)  CTMETH='block '
      IF(TMETH .EQ. SSOR)   CTMETH='ssor  '
      IF(TMETH .EQ. ICCG)   CTMETH='iccg  '
c
c
ckk -minimize- begin
c  900 FORMAT(15X,'USING VALUES FOR TINKER OPTIMIZATION '/,
c     #      23X,' FROM $TOPTMZ CARD  ',/
c     #      10x,'MODE= ',A8,5X,'METHOD= ',A8,'GRDMIN= ',F10.6)
  900 FORMAT(15X,'USING VALUES FOR TINKER OPTIMIZATION '/,
     #      23X,' FROM $TOPTMZ CARD  ',/
     #      10x,'OPTPRG= ',A8,'   MODE= ',A8,' METHOD= ',
     #      A8,' GRDMIN= ',F10.6)
ckk -minimize- end
  901 FORMAT(/,15X,' *** NO $TOPTMZ CARD FOUND ****'/,
     #      5X,'USING DEFAULT VALUES FOR TINKER OPTIMIZATION ',/
     #      10x,'MODE= ',A8,5X,'METHOD= ',A8,'GRDMIN= ',F10.6)
  902 FORMAT(/,15X,' *** ALL MM ATOMS ARE FREEZED **',/)
ckk -minimize- begin
  903 FORMAT(15X,'USING VALUES FOR TINKER OPTIMIZATION '/,
     #      23X,' FROM $TOPTMZ CARD  ',/
     #      10x,'OPTPRG= ',A8,' MINOPT= ',A8,' GRDMIN= ',F10.6)
ckk -minimize- end
  910 FORMAT(10X,'MAXHES= ',I8,' MODPAR= ',I8)
C
        return
        end
c   This routine changed by Cheol, June 2003.
c   With this modification, user no more needs to input LINKGE
c   definition, which has been quite troublesome.
c
c  ##################################################################
c  #                                                                #
c  #     linkin: Read $LINK card for info on how to link the QM     #
c  #             and MM regions of a QM/MM optimization             #
c  #                                                                #
c  ##################################################################
c
      SUBROUTINE LINKIN(mode)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      include 'sizes.i'
      include 'usage.i'
      include 'atoms.i'
      include 'couple.i'
casa -restrain- begin
      include 'restrn.i'
casa -restrain- end
c
ckk -fmoimomm- begin
ckk     PARAMETER (NNAM=3)
      PARAMETER (NNAM=4)
ckk -fmoimomm- end
      PARAMETER (MAXLNK=100, MAXR1=2000,zero=0.0d+00, one=1.0d+00)
c
      LOGICAL GOPARR,DSKWRK,MASWRK,mmonly,qmmm
      LOGICAL IMOMM,SIMOMM,flag
      INTEGER TAG,CLASS,ATOMIC,VALENCE
      REAL*8 MASS
      CHARACTER*10 NAME, xyzname
      CHARACTER*20 STORY
      double precision link,blqm,blmm

      DIMENSION QNAM(NNAM),KQNAM(NNAM)
      DIMENSION ITMPLNK(3*MAXLNK)
ckk -fmoimomm- begin
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY
c     ibastmp(i), ibastmp(i+1) stores seq. # of atoms and basis set #, 
c     respectively.
      dimension ibastmp(2*MAXR1)
ckk -fmoimomm- end
C
C -- MAXIMUM SIZES and COMMONS FOR QMMM LINKING
c
      COMMON /QMMM1/ IMOMM,SIMOMM,NPAIR,NSEQ
      COMMON /QMMM2/ IQMATM(MAXR1),ibasfmo(MAXR1)
      COMMON /QMMM3/ LINKge(3*MAXLNK),blqm(MAXLNK),blmm(MAXLNK)
      COMMON /ATMTYP/ MASS(MAXATM),TAG(MAXATM),CLASS(MAXATM),
     *                ATOMIC(MAXATM),VALENCE(MAXATM),NAME(MAXATM),
     *                STORY(MAXATM), xyzname(maxatm)
C
c --  COMMONS FROM GAMESS
c
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /tinopt/ mparti,mmonly, qmmm
      COMMON /ZMTALT/ NZMAT2,NZVAR2,NVAR2,NZMTRD,ICOORD
C
      DATA LINK /8HLINK    /
ckk -fmoimomm- begin
ckk      DATA QNAM /8HIMOMM   ,8HSIMOMM  ,8HIQMATM  /
ckk      DATA KQNAM/0,0,1/
      DATA QNAM /8HIMOMM   ,8HSIMOMM  ,8HIQMATM  ,8HIBAS    /
      DATA KQNAM/0,0,1,1/
ckk -fmoimomm- end
c
ckk -fmoimomm- begin
ckk      KQNAM(3) = 10*100+1
      KQNAM(3) = 10*MAXR1+1
      KQNAM(4) = 10*2*MAXR1+1
ckk -fmoimomm- end
c
C    --- INITIALIZE SOME VALUES
c
      imomm = .false.
      simomm = .false.
      LINKGE(1)=0
      do i=1,3*maxlnk
         linkge(i) = 0
      enddo
      do i=1,maxlnk
         blqm(i) = zero
         blmm(i) = zero
      enddo
c
c     user has to specify at least some QM atoms below!
c
      do i=1,maxr1
         iqmatm(i) = 0
      enddo
ckk -fmoimomm- begin
      if(nfg.ne.0) then
         do i=1,2*maxr1
            ibastmp(i)=0
         enddo
      endif
ckk -fmoimomm- end
C
C     ----- READ NAMELIST $LINK -----
c
      JRET = 0
      CALL NAMEIO(IR,JRET,LINK,NNAM,QNAM,KQNAM,
ckk -fmoimomm- begin
ckk     *     IMOMM,SIMOMM,IQMATM,0,
     *     IMOMM,SIMOMM,IQMATM,IBAStmp,
ckk -fmoimomm- end
     *     0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     *     0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     *     0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     *     0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0)
C
      if (jret.eq.0) then
         qmmm  = .true.
         mmonly= .false.
         if (maswrk .and. qmmm) write(IW,*)
     *        '---- QMMM procedure is ON ----'
      end if
ckk -fmoimomm-      IF (JRET.EQ.1) THEN
ckk -fmoimomm-         qmmm  = .false.
ckk -fmoimomm-         return
ckk -fmoimomm-      end if
      IF (JRET.EQ.1) THEN
         qmmm  = .false.
         return
      end if
      IF (JRET.GE.2) THEN
         IF(MASWRK) WRITE(IW,*) 'TROUBLE INTERPRETING $LINK '
         CALL ABORT
      END IF
c
c         the argument mode=0 is just to find out if this group exists,
c         in a later call, it is fully parsed (after -tinkin- called)
c
      if(mode.eq.0) return
c
C     Check IQMATM
c     Note that currently the maximum number of IQMATM is set to 100.
c     If you have more than 100 QM atoms, you must increase MAXR1.
c
ckk -iqmatm-
      ichk=0
      if(iqmatm(1).ne.0) then
         do i=1,MAXR1
            if(iqmatm(i).lt.0) ichk=ichk+1
         enddo
      endif
      if(ichk.gt.0) call iqmgind
ckk -iqmatm-
      nseq=0
      DO I=1,MAXR1
         IF(iqmatm(I) .NE. 0) nseq=nseq+1
      enddo
      IF (Nseq .eq.  0) THEN
        IF (MASWRK) WRITE(IW,910)
        CALL ABORT
      ENDIF
ckk -fmoimomm- begin
      if(nfg.ne.0) then
         ntemp=0
         do i=1,MAXR1
            ibasfmo(i)=1
            if(ibastmp(2*(i-1)+1).ne.0) ntemp=ntemp+1
         enddo
         maxqma=0
         do i=1,nseq
            if(iqmatm(i).gt.maxqma) maxqma=iqmatm(i)
         enddo
         do i=1,2*MAXR1,2
            if(ibastmp(i).ne.0) then
               if(ibastmp(i).gt.maxqma) then
                   write(iw,980) i,ibastmp(i),maxqma
                   call abrt
               endif
               ifound=1
               do j=1,nseq
                  if(ibastmp(i).eq.iqmatm(j)) then
                     ibasfmo(j)=ibastmp(i+1)
                     ifound=0
                  endif
               
               enddo
               if(ifound.eq.1) then
                  write(iw,990) i,ibastmp(i)
                  call abrt
               endif
            endif
         enddo
      endif
ckk -fmoimomm- end
C  Check IMOMM and SIMOMM, only one can be chosen, but one must be.  
c
      If ((simomm.or.imomm)  .and.  .not.(simomm.and.imomm)) then
        continue 
      else
        IF (MASWRK) WRITE(IW,915)
        CALL ABRT
      endif
c
C     Check Linkge
c
      if (MASWRK) write(IW,935)
c
C     The program generates LINKGE definition on the basis of
c     IQMATM
c     First value of a pair is the QM atom sequence of TINKER part.
c     Second value of a pair is the MM atom sequence of TINKER part.
C
      lkgidx=0
      do i=1,nseq
        do j=1,n12(iqmatm(i))
           flag=.false.
           do k=1,nseq
              if (i12(j,iqmatm(i)).eq.iqmatm(k)) then
                 flag=.true.
              endif
           enddo
           if (.not.flag) then
              lkgidx=lkgidx+1
              linkge(lkgidx*2-1)=iqmatm(i)
              linkge(lkgidx*2)=i12(j,iqmatm(i))
           endif
        enddo
      enddo
c
C     write out the linkage pairs
c
      do i=1,lkgidx
         if (MASWRK) write(IW,945) i,linkge(i*2-1),linkge(i*2)
      enddo
c
c     Count the number of pairs
c     Redefine LINKGE values.
c     IF icoord.eq.5 (TINKER coordinate)
c        First value of a pair is the QM atom sequence of TINKER region
c     ELSE IF iccord.eq.-1 (UNIQUE)
c        First value of a pair is the QM atom sequence of QM region
c     ENDIF
C     IF icoord.eq.5 then
c        First value of a pair will be redefined as the QM atom sequence
c        of QM region in TIN2GMS routine
c     ENDIF
c     Second value of a pair is the H atom sequence of QM region
c     Third value of a pair is the MM atom sequence of TINKER region
c
      npair=0
      DO I=1,MAXLNK
         IF(linkge(I*2) .NE. 0) npair=npair+1
      enddo
ckk -fmoimomm- begin
      IF (imomm .and. Npair .eq. 0) then
        IF (MASWRK) WRITE(IW,920)
        CALL ABORT
      else
ckk -fmoimomm- begin
       if(npair.gt.0) then
ckk -fmoimomm- end
         DO I=1,NPAIR
            LR1=LINKGE(2*I-1)
            LR2=LINKGE(2*I  )
            if (icoord.eq.5) then
               ITMPLNK(I*3-2)=LR1
            elseif (icoord.eq.-1) then
               DO J=1,NSEQ
                  IF (LR1.EQ.IQMATM(J)) ITMPLNK(I*3-2)=J
               enddo
            endif
            ITMPLNK(I*3-1)=I+NSEQ
            ITMPLNK(I*3  )=LR2
C
C        Construct BLQM and BLMM
c
            zan1 = atomic(lr1)
            zan2 = atomic(lr2)
            CALL GTDIST(zan1, ONE,RDIST)
            BLQM(I)=RDIST
            call GTDIST(zan1,zan2,RDIST)
            BLMM(I)=RDIST
         enddo
         DO I=1,NPAIR*3
            LINKGE(I)=ITMPLNK(I)
         enddo
ckk -fmoimomm- begin
        endif
ckk -fmoimomm- end
      endif
c
casa  -restrain- begin
      if (nfg .eq. 0) then
C
C     Generate inactive atoms in MM calculations automatically..
C
          if (maswrk) write(iw,970)
          do i=1,nseq
             use(iqmatm(i))=.false.
             nuse=nuse-1
          enddo
          if (maswrk) write(iw,940) (iqmatm(i),i=1,nseq)
      else
          do i=1,nseq
          if (.not.use(iqmatm(i))) then
              use(iqmatm(i))=.true.
              nuse=nuse+1
          end if
          end do
          npfixm = npfix
      end if
casa  -restrain- end
c     IMOMM requires to inactivate some more atoms in MM region
      if (imomm) then
         do i=1,npair
              use(linkge(i*3))=.false.
              nuse=nuse-1
         enddo
      endif
C
c     Write LINK options to output file
c
      if (maswrk) then
         if(simomm) then
            write(iw,930) 'SIMOMM'
         elseif (imomm)  then
            write(iw,930) 'IMOMM'
         endif
         write(iw,936)
         write(iw,940) (iqmatm(i),i=1,nseq)
         if (imomm) then
            write(iw,950)
            do i=1,npair
               write(iw,960)  i,blqm(i),blmm(i)
            enddo
         endif
      end if
C
      Return
c
  910 FORMAT(15X,'*** ERROR IN LINK SPECIFICATION *** ',/
     *   15X,' QMMM REQUIRES at least one QM atom')
  915 FORMAT(15X,'*** ERROR IN LINK SPECIFICATION *** ',/
     *   15X,' QMMM REQUIRES to specify SIMOMM or IMOMM')
  920 FORMAT(15X,'*** ERROR IN LINK SPECIFICATION *** ',/
     *   15X,' QMMM REQUIRES at least one Linkge')
  930 FORMAT(/10x,a6,' LINKING SELECTED ',/)
  935 FORMAT(/10X,'LINKING ATOM PAIRS'/10X,18(1H-)/
     *   1X,10X,'NO',4X,'QM-ATOM',4X,'MM-ATOM')
  936 format (10x,'ALL QM ATOMS',/10x,12(1H-))
  940 FORMAT(6(7x,i5))
  945 FORMAT(5x,i8,3x,i8,3x,i8)
  950 format(/10x,'BRIDGE QM and MM BOND LENGTH FOR IMOMM RUN',/
     *   10x,42(1H-)/,15x,'NO',6x,'QM REGION',6x,'MM REGION')
  960 FORMAT(9x,i8,8x,f6.4,9x,f6.4)
  970 format(/10X,'INACTIVE ATOMS IN MM REGION'/10X,27(1H-))
ckk -fmoimomm- begin
  980 format(/10X,'ibas(i) exceeds nseq. i=',i5,' ibas(i)=',i5,
     *             ' max qm atom number ',i5/)
  990 format(/10X,'ibas(i) atom is not found in QM atoms. i=',i5,
     *            ' ibas(i)=',i5,/)
ckk -fmoimomm- end
      END
ckk -fmoimomm- added a new subroutine for indat data generation fo fmo-imomm
      subroutine fmommind(maxig,indat,indatg)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL IMOMM,SIMOMM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
C -- MAXIMUM SIZES and COMMONS FOR QMMM LINKING
      PARAMETER (MAXATM=12000, MAXLNK=100, MAXR1=2000)
      double precision blqm,blmm
      COMMON /QMMM1/ IMOMM,SIMOMM,NPAIR,NSEQ
      COMMON /QMMM2/ IQMATM(MAXR1),ibasfmo(MAXR1)
      COMMON /QMMM3/ LINKge(3*MAXLNK),blqm(MAXLNK),blmm(MAXLNK)
c -fmoimomm- begin
      dimension indat(*),indatg(*)
c     note: itemp is used as scratch here
      dimension itemp(MAXATM)
c -fmoimomm- end
c
c     process Gaussian-like INDAT, indicated by indat(1)=0 (therefore, skip
c     indat(1))
c
ckkmod         do i=1,MAXATM
ckkmod            itemp(i)=0
ckkmod         enddo
cdeb
         ifg=1
         nifg=0
         natot=0
         i=1
         iqma=0
  100    continue
           i=i+1
           if(i.ge.maxig) goto 200
           now=indat(i)
           if(now.eq.0) then
             if(ifg.eq.nfg) goto 200
             if(nifg.eq.0) then
               if(maswrk) write(iw,9011) ifg
               call abrt
             endif
             ifg=ifg+1
             natot=natot+nifg
             nifg=0
             goto 100
           endif
           if(indat(i+1).lt.0) then
             i=i+1
             next=abs(indat(i))
           else
             next=now
           endif
           do j=now,next
             nifg=nifg+1
c            indatg(j)=ifg
c            check j is qmatom ?
             iqma=0
             do k=1,nseq
ckkmod                if(j.eq.iqmatm(k)) iqmchk=0
                if(j.eq.iqmatm(k)) iqma=k
             enddo
             if(iqma.eq.0) then
                if(maswrk) write(iw,9020) now,next,ifg
                call abrt
             endif
ckkmod             iqma=iqma+1
             indatg(iqma)=ifg
ckkmod             itemp(j)=ifg
           enddo
         goto 100
  200    continue
c     link atoms
         if(npair.gt.0) then
            do i=1,npair
               lr1=linkge(3*i-2)
               ifgt=indatg(lr1)
               if(ifgt.eq.0) then
                  write(iw,9030) lr1
                  call abrt
               else
                  indatg(nseq+i)=ifgt
               endif
            enddo
         endif
         natot=natot+nifg+npair
         if(ifg.ne.nfg.or.natot.ne.natfmo) then
           if(maswrk) write(iw,9000) ifg,nfg,natot,natfmo
           if(maswrk) write(iw,9010) nseq,npair 
           call abrt
         endif
         call icopy(natfmo,indatg,1,indat,1)
c
      return
 9000 format(/1x,'Bad indat: nfg(indat,nfg)=',2I5,
     *           ' natfmo(indat,fmoxyz)=',2I5,/)
 9010 format(/1x,'nseq=',i5,' npair=',i5)
 9011 format(/1x,'No atoms in fragment ifg=',i5)
 9020 format(/1x,'Wrong indat for QM atoms in FMO-IMOMM: now,next,
     *             ifg ',3i5)  
 9030 format(/1x,'No QM atom linked to LR1=',i5)
      end
      subroutine fmommbon(indat,iabdfg,jabdfg)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL IMOMM,SIMOMM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
C -- MAXIMUM SIZES and COMMONS FOR QMMM LINKING
      PARAMETER (MAXR1=2000)
      double precision blqm,blmm
      COMMON /QMMM1/ IMOMM,SIMOMM,NPAIR,NSEQ
      COMMON /QMMM2/ IQMATM(MAXR1),ibasfmo(MAXR1)
c -fmoimomm- begin
      dimension indat(*),iabdfg(*),jabdfg(*)
c -fmoimomm- end
c
c     renumber bda and baa atoms
         if(nbdfg .gt. 0) then
            do i=1,nbdfg
               iatnew=0
               jatnew=0
               do j=1,nseq
                  if(iqmatm(j).eq.abs(iabdfg(i))) then
                     iatnew=j
                  endif
                  if(iqmatm(j).eq.abs(jabdfg(i))) then
                     jatnew=j
                  endif
               enddo
               if(iatnew.ne.0 .and. jatnew.ne.0) then
                  iabdfg(i)=-iatnew
                  jabdfg(i)=jatnew
               else
                  write(iw,9040) iabdfg(i),jabdfg(i)
                  call abrt
               endif
           enddo
c       check data
           do i=1,nbdfg
             if(indat(abs(iabdfg(i))).eq.indat(jabdfg(i))) then
               if(maswrk) write(iw,9050) i,iabdfg(i),
     *                                jabdfg(i),indat(abs(jabdfg(i)))
               call abrt
             endif
           enddo
         endif
c
      return
 9040 format(/1x,'Bad FMOBND data. iabdfg, jabdfg',2i5)
 9050 format(/1x,'Check entree',I5,
     *          ': intrafragment fractioned bonds are not allowed.',
     *       /1x,2I5,' are in fragment',I5,/)
      end
c
      subroutine iqmgind
      PARAMETER (MAXR1=2000)
      COMMON /QMMM2/ IQMATM(MAXR1),ibasfmo(MAXR1)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
c
      ndat=0
      do i=1,MAXR1
         if(iqmatm(i).ne.0) ndat=ndat+1
      enddo
      do i=1,ndat
         iqmatm(MAXR1-ndat+i)=iqmatm(i)
      enddo
c
      i=MAXR1-ndat
      idat=0
   10 i=i+1
      if(iqmatm(i).eq.0) go to 100
      if(iqmatm(i).gt.0) then
         idat=idat+1
         if(idat.gt.MAXR1) then
            write(iw,9000) MAXR1
            call abrt
         endif
         iqmatm(idat)=iqmatm(i)
      else
         if(-iqmatm(i).lt.iqmatm(idat)) then
            write(iw,9100) iqmatm(idat),iqmatm(i)
            call abrt
         endif
         nd=-iqmatm(i)-iqmatm(idat)
         ini=iqmatm(idat)
         do j=1,nd
            idat=idat+1
            if(idat.gt.MAXR1) then
               write(iw,9000) MAXR1
               call abrt
            endif
            iqmatm(idat)=ini+j
         enddo
      endif
      if(i.lt.MAXR1) go to 10
  100 continue
      do i=idat+1,MAXR1
         iqmatm(i)=0
      enddo
c
      return
 9000 format(10x,' error: too many QM atoms. MAXR1=',i8)
 9100 format(10x,' error in iqmatm. to value is smaller than from.
     . iqmatm(from) and iqmatm(to) are ',2i6) 
      end
