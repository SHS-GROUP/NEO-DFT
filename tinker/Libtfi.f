C  9 MAR 00 - CHC - FIX for parallel run
C 14 Oct 98 - CHC - change title -> ttitle
C                   getkey : modified for parallel run
c  7 May 98 - JRS - getprm: Change error handling if no 
c                   parameter file found. Use SEQOPN, SEQREW,
c                   and SEQCLO instead of open, rewind, and close
c  7 May 98 - JRS - getkey: modified to read GAMESS input file 
c  6 May 98 - JRS - initial: call command commented out
c  5 May 98 - JRS - final: renamed to finalt
c
c
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine fatal  --  terminate the program abnormally  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "fatal" terminates execution due to a user request, a severe
c     error or some other nonstandard condition
c
c
      subroutine fatal
      implicit none
      include 'iounit.i'
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     print a final warning message, then quit
c
      if (maswrk) write (iout,10)
   10 format (/,' TINKER is Unable to Continue; Terminating the',
     &           ' Current Calculation',/)
      stop
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine field  --  get the potential energy functions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "field" sets the force field potential energy functions from
c     a parameter file and modifications specified in a keyfile
c
c
      subroutine field
      implicit none
      include 'sizes.i'
      include 'angpot.i'
      include 'bndpot.i'
      include 'chgpot.i'
      include 'keys.i'
      include 'math.i'
      include 'polpot.i'
      include 'potent.i'
      include 'rxnpot.i'
      include 'torpot.i'
      include 'urypot.i'
      include 'vdwpot.i'
      integer i
      character*80 record
c
c
c     set the default values for the active potentials
c
      use_bond = .true.
      use_angle = .true.
      use_strbnd = .true.
      use_urey = .true.
      use_angang = .true.
      use_opbend = .true.
      use_improp = .true.
      use_imptor = .true.
      use_tors = .true.
      use_strtor = .true.
      use_tortor = .false.
      use_vdw = .true.
      use_charge = .true.
      use_chgdpl = .true.
      use_dipole = .true.
      use_mpole = .true.
      use_polar = .true.
      use_rxnfld = .false.
      use_solv = .true.
      use_geom = .true.
      use_extra = .true.
c
c     set default values for force field control parameters
c
      bndtyp = 'TAYLOR'
      bndunit = 1.0d0
      cbnd = 0.0d0
      qbnd = 0.0d0
      angtyp = 'HARMONIC'
      angunit = 1.0d0 / radian**2
      cang = 0.0d0
      qang = 0.0d0
      pang = 0.0d0
      sang = 0.0d0
      aaunit = 1.0d0 / radian**2
      opbunit = 1.0d0 / radian**2
      stbnunit = 1.0d0
      ureyunit = 1.0d0
      torsunit = 1.0d0
      storunit = 1.0d0
      vdwtyp = 'LENNARD-JONES'
      aterm = 0.0d0
      bterm = 0.0d0
      cterm = 0.0d0
      vdw12use = 1
      vdw13use = 1
      vdw14use = -1
      vdwscale = 1.0d0
      radrule = 'ARITHMETIC'
      radtyp = 'R-MIN'
      radsiz = 'RADIUS'
      epsrule = 'GEOMETRIC'
      gausstyp = 'NONE'
      ngauss = 0
      dielec = 1.0d0
      chg12use = 1
      chg13use = 1
      chg14use = -1
      chgscale = 1.0d0
      neutnbr = .false.
      neutcut = .false.
      poltyp = 'MUTUAL'
      poleps = 0.000001d0
      pgamma = 1.0d0
      pradius = 1.662d0
      rfsize = 1000000.0d0
      rfbulkd = 80.0d0
      rfterms = 1
c
c     read the potential energy force field parameter file
c
      call getprm
c
c     check keywords for potential function control parameters
c
      do i = 1, nkey
         record = keyline(i)
         call prmkey (record)
      end do
      return
      end
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine final  --  final actions before program exit  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "final" performs any final program actions, prints a status
c     message, and then pauses if necessary to avoid killing the
c     execution window
c
c
cjrs
      subroutine finalt
cjrs
      implicit none
      include 'inform.i'
      include 'iounit.i'
c
c
c     print a final status message before exiting TINKER
c
      if (debug) then
         write (iout,10)
   10    format (/,' TINKER is Exiting following Normal Termination',
     &              ' of the Program',/)
      end if
c
c     may need a pause to avoid killing the execution window
c
c     pause
      return
      end
c
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function freeunit  --  gets an unopened logical unit  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "freeunit" finds an unopened Fortran I/O unit and returns
c     its numerical value from 1 to 99; the units already assigned
c     to "input" and "iout" (usually 5 and 6) are skipped since
c     they have special meaning as the default I/O units
c
c
      function freeunit ()
      implicit none
      include 'iounit.i'
      integer freeunit
      logical used
c
c
c     try each logical unit until an unopened one is found
c
      freeunit = 0
      used = .true.
      dowhile (used)
         freeunit = freeunit + 1
         if (freeunit.ne.input .and. freeunit.ne.iout) then
            if (freeunit .gt. 99) then
               write (iout,10)
   10          format (' FREEUNIT  --  No Available Fortran I/O Units')
               call fatal
            end if
            inquire (unit=freeunit,opened=used)
         end if
      end do
      return
      end
c
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine getime  --  get elapsed CPU time in seconds  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "getime" gets elapsed CPU time in seconds for an interval
c
c
      subroutine getime (elapsed)
      implicit none
      include 'chrono.i'
      real*8 elapsed
c
c
c     elapsed time for the interval is the current total CPU
c     time minus the total time at the start of the interval
c
      call clock (elapsed)
      elapsed = elapsed - cputim
      return
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getint  --  get internal coordinate structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getint" asks for an internal coordinate file name, then reads
c     the internal coordinates and computes Cartesian coordinates
c
c
      subroutine getint
      implicit none
      include 'iounit.i'
      include 'output.i'
      integer izmt
      integer freeunit
      character*60 intfile
      logical exist
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (intfile,exist)
      if (exist) then
         call basefile (intfile)
         call suffix (intfile,'int')
         call version (intfile,'old')
         inquire (file=intfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      dowhile (.not. exist)
         write (iout,10)
   10    format (/,' Enter Internal Coordinate File Name :  ',$)
         read (input,20)  intfile
   20    format (a60)
         call basefile (intfile)
         call suffix (intfile,'int')
         call version (intfile,'old')
         inquire (file=intfile,exist=exist)
      end do
c
c     first open and then read the internal coordinates file
c
      coordtype = 'internal'
      izmt = freeunit ()
      open (unit=izmt,file=intfile,status='old')
      rewind (unit=izmt)
      call readint (izmt)
      close (unit=izmt)
c
c     convert internal to Cartesian coordinates
c
C-MWS-change spelling below from connect to avoid confusion with
C-MWS-the TCP/IP socket subroutine of the same name.
c
      call connct
      call makexyz
      return
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getkey  --  find and store contents of keyfile ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getkey" finds a valid keyfile and stores its contents as
c     line images for subsequent keyword parameter searching
c
c
      subroutine getkey (ikey,iw)
      implicit double precision (a-h,o-z)
      include 'sizes.i'
      include 'argue.i'
      include 'files.i'
      include 'iounit.i'
      include 'keys.i'
      integer i,ifnd,ieof,ikey,next,length
      integer cmsg(80),freeunit,trimtext
      character*20 keyword
      character*60 keyfile,string
      character*80 temrec,record,comment
      logical flag, exist,header
      INTEGER me,master,nproc,ibtyp,iptim,iflag
      CHARACTER*8 GRPNAM
C
C
      LOGICAL LOG,FOUND,GOPARR,DSKWRK,MASWRK,TDSKWRK
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c
      nkey = 0
      flag=.true.
      dowhile (flag)
      IF (MASWRK) THEN
         read (ikey,20,end=40)  record
         temrec=record
         call upcase(temrec)
         iflag=index(temrec,'END')
         write(6,20) record
   20    format (a80)
C
C     ----- CORRECT INPUT GROUP FOUND ----
C
         IF (GOPARR) THEN
            DO 200 I=1,80
               CMSG(I) = ICHAR(record(I:I))
  200       CONTINUE
         END IF
      else
         if (goparr) then
         endif
C CHC check this out
      ENDIF
c
      if(goparr) call ddi_bcast(250,'I',iflag,1,master)
      if(iflag.ne.0) flag=.false.
      IF (GOPARR) CALL DDI_BCAST(251,'I',CMSG,80,MASTER)
      IF (.NOT.MASWRK) THEN
         DO 220 I=1,80
            record(I:I) = CHAR(CMSG(I))
  220    CONTINUE
      END IF
         length = trimtext (record)
         if (length .ne. 0) then
            nkey = nkey + 1
            keyline(nkey) = record
         end if
         if (nkey .ge. maxkey) then
            if (maswrk) write (iw,30)
   30       format (/,' GETKEY  --  Keyfile Too Large;',
     &                 ' Increase MAXKEY')
            call fatal
         end if
      end do
   40    continue
c
c     check for comment lines to be echoed to the output
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ECHO ') then
            comment = record(next:80)
            length = trimtext (comment)
            if (header) then
               header = .false.
               if (maswrk) write (iw,50)
   50          format ()
            end if
            if (length .eq. 0) then
               if (maswrk) write (iw,60)
   60          format ()
            else
               if (maswrk) write (iw,70)  comment(1:length)
   70          format (a)
            end if
         end if
      end do
      return
      end
c
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine getmol2  --  get a Sybyl MOL2 format file  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "getmol2" asks for a Sybyl MOL2 molecule file name,
c     then reads the coordinates from the file
c
c
      subroutine getmol2
      implicit none
      include 'files.i'
      include 'iounit.i'
      integer isyb
      integer freeunit
      character*60 sybylfile
      logical exist
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (sybylfile,exist)
      if (exist) then
         call basefile (sybylfile)
         call suffix (sybylfile,'mol2')
         call version (sybylfile,'old')
         inquire (file=sybylfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      dowhile (.not. exist)
         write (iout,10)
   10    format (/,' Enter a Sybyl MOL2 File Name :  ',$)
         read (input,20)  sybylfile
   20    format (a60)
         call basefile (sybylfile)
         call suffix (sybylfile,'mol2')
         call version (sybylfile,'old')
         inquire (file=sybylfile,exist=exist)
      end do
c
c     first open and then read the Sybyl MOL2 coordinates file
c
      isyb = freeunit ()
      open (unit=isyb,file=sybylfile,status='old')
      rewind (unit=isyb)
      call readmol2 (isyb)
      close (unit=isyb)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getnumb  --  extract an integer from a string  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getnumb" searchs an input string from left to right for an
c     integer and puts the numeric value in "number"; returns zero
c     with "next" unchanged if no integer value is found
c
c     string    input character string to be searched
c     number    output with the first integer in the string
c     next      input with first position of search string;
c                  output with position following the number
c
c
      subroutine getnumb (string,number,next)
      implicit none
      integer i,j,next,trimtext,length
      integer first,last,initial,final
      integer number,digit,place(10)
      character*1 letter
      character*(*) string
      logical negate,numeral
      data place  / 1, 10, 100, 1000, 10000, 100000, 1000000,
     &              10000000, 100000000, 1000000000 /
c
c
c     initialize number and get the input text string length
c
      number = 0
      negate = .false.
      numeral = .false.
      length = trimtext(string(next:))
c
c     move through the string one character at a time,
c     searching for the first run of numeric characters
c
      first = 1
      last = 0
      initial = next
      final = next + length - 1
      do i = initial, final
         letter = string(i:i)
         if (letter.ge.'0' .and. letter.le.'9') then
            if (.not. numeral) then
               numeral = .true.
               first = i
            end if
            if (i .eq. final) then
               last = final
               next = i + 1
            end if
         else if (letter.eq.'-' .and. .not.negate) then
            negate = .true.
         else if (numeral) then
            if (letter.eq.' ' .or. letter.eq.',' .or. 
     &          letter.eq.';' .or. letter.eq.'_') then
               last = i - 1
               next = i
               goto 10
            end if
         else if (letter.ne.' ' .or. negate) then
            goto 10
         end if
      end do
   10 continue
c
c     trim the actual number if it is too big to return
c
      last = min(last,first+9)
c
c     convert the text numeral into an integer number
c
      j = 0
      do i = last, first, -1
         j = j + 1
         if (string(i:i) .eq. '0') then
            digit = 0
         else if (string(i:i) .eq. '1') then
            digit = 1
         else if (string(i:i) .eq. '2') then
            digit = 2
         else if (string(i:i) .eq. '3') then
            digit = 3
         else if (string(i:i) .eq. '4') then
            digit = 4
         else if (string(i:i) .eq. '5') then
            digit = 5
         else if (string(i:i) .eq. '6') then
            digit = 6
         else if (string(i:i) .eq. '7') then
            digit = 7
         else if (string(i:i) .eq. '8') then
            digit = 8
         else if (string(i:i) .eq. '9') then
            digit = 9
         end if
         number = number + digit * place(j)
      end do
      if (negate)  number = -number
      return
      end
c
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine getpdb  --  get a Protein Data Bank file  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "getpdb" asks for a Protein Data Bank file name,
c     then reads in the coordinates file
c
c
      subroutine getpdb
      implicit none
      include 'iounit.i'
      integer ipdb
      integer freeunit
      character*60 pdbfile
      logical exist
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (pdbfile,exist)
      if (exist) then
         call basefile (pdbfile)
         call suffix (pdbfile,'pdb')
         call version (pdbfile,'old')
         inquire (file=pdbfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      dowhile (.not. exist)
         write (iout,10)
   10    format (/,' Enter Protein Data Bank File Name :  ',$)
         read (input,20)  pdbfile
   20    format (a60)
         call basefile (pdbfile)
         call suffix (pdbfile,'pdb')
         call version (pdbfile,'old')
         inquire (file=pdbfile,exist=exist)
      end do
c
c     first open and then read the PDB coordinates file
c
      ipdb = freeunit ()
      open (unit=ipdb,file=pdbfile,status='old')
      rewind (unit=ipdb)
      call readpdb (ipdb)
      close (unit=ipdb)
      return
      end
c
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine getprm  --  get force field parameter file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "getprm" finds the potential energy parameter file
c     and then opens and reads the parameters
c
c
      subroutine getprm
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'iounit.i'
      include 'keys.i'
      integer i,iexist,iflag,iprm,next,freeunit
      character*4 none
      character*20 keyword
      character*60 prmfile
      character*80 record
      logical exist,useprm
      INTEGER me,master,nproc,ibtyp,iptim
      CHARACTER*8 GRPNAM
C
C
      LOGICAL LOG,FOUND,GOPARR,DSKWRK,MASWRK,TDSKWRK
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     set the default name for the parameter file
c
      useprm = .true.
      iflag=0
      prmfile = filename(1:leng)//'.prm'
c
c     search the keyword list for the parameter filename
c
      IF (maswrk) then
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'PARAMETERS ') then
            call gettext (record,prmfile,next)
         end if
      end do
c
c     check existence of default or specified parameter file
c
      call suffix (prmfile,'prm')
      call version (prmfile,'old')
      inquire (file=prmfile,exist=exist)
c
c     test for user specified absence of a parameter file
c
      if (exist) then
         iexist=0
      else
         iexist=1
         none = prmfile(1:4)
         call upcase (none)
         if (none .eq. 'NONE') then
            exist = .true.
            iexist=0
            useprm = .false.
            iflag=1
         end if
      end if
c
      endif
      IF (GOPARR) CALL DDI_BCAST(7777,'I',iexist,1,MASTER)                       
c
      if (iexist.eq.1) then
       if (maswrk) then
         WRITE(6,*) '  '
         WRITE(6,*) '          *** ERROR IN GETPRM ***  '
         WRITE(6,*) '  SPECIFIED MM PARAMETER FILE DOES NOT EXIST  '
         WRITE(6,*) '        CHECK INPUT FILE AND RUN AGAIN   '
        endif
         CALL FATAL
      ENDIF
c
      IF (GOPARR) CALL DDI_BCAST(7777,'I',iflag,1,MASTER)
C
c cjrs 
      call initprm
      if (iflag.eq.0) then
         iprm = freeunit ()
c
cjrs need to use seqopn, seqrew, and seqclo for parallel compatability 
c for now we use open until I figure out seqopn
c
        if (maswrk) then
         open (unit=iprm,file=prmfile,status='old')
         rewind (unit=iprm)
        endif
c
          call readprm (iprm)
        if (maswrk) close (unit=iprm)
c
      end if
      return
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getref  --  get structure from reference area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getref" copies structure information from the reference area
c     into the standard variables for the current system structure
c
c
      subroutine getref
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'refer.i'
      include 'titles.i'
      integer i,j
c
c
c     retrieve the filename and title line for the structure
c
      filename = reffile
      leng = refleng
      ttitle = reftitle
      ltitle = refltitle
c
c     retrieve the coordinates, type and connectivity of each atom
c
      n = nref
      do i = 1, n
         name(i) = refnam(i)
         x(i) = xref(i)
         y(i) = yref(i)
         z(i) = zref(i)
         type(i) = reftyp(i)
         n12(i) = n12ref(i)
         do j = 1, n12(i)
            i12(j,i) = i12ref(j,i)
         end do
      end do
      return
      end
c
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine getstring  --  extract single quoted string  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "getstring" searchs for a quoted text string within an input
c     character string; the region between the first and second
c     quotes is returned as the "text"; if the actual text is too
c     long, only the first part is returned
c
c     string    input character string to be searched
c     text      the quoted text found in the input string
c     next      input with first position of search string;
c                  output with position following text
c
c
      subroutine getstring (string,text,next)
      implicit none
      integer i,j,len,length,size,next
      integer extent,first,last,initial,final
      character*1 letter
      character*(*) string,text
c
c
c     get the length of input string and output text
c
      length = len(string(next:))
      size = len(text)
c
c     move through the string one character at a time,
c     searching for the quoted text string characters
c
      first = 1
      last = 0
      initial = next
      final = next + length - 1
      do i = initial, final
         letter = string(i:i)
         if (letter .eq. '"') then
            first = i + 1
            do j = first, final
               letter = string(j:j)
               if (letter .eq. '"') then
                  last = j - 1
                  next = j + 1
                  goto 10
               end if
            end do
         end if
      end do
   10 continue
c
c     trim the actual word if it is too long to return
c
      extent = last - first + 1
      final = first + size - 1
      if (extent .gt. size)  last = final
c
c     transfer the text into the return string
c
      j = 0
      do i = first, last
         j = j + 1
         text(j:j) = string(i:i)
      end do
      do i = last+1, final
         j = j + 1
         text(j:j) = ' '
      end do
      return
      end
c
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine gettext  --  extract text from a string  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "gettext" searchs an input string for the first string of
c     non-blank characters; the region from a non-blank character
c     to the first blank space is returned as "text"; if the
c     actual text is too long, only the first part is returned
c
c     string    input character string to be searched
c     text      output with the first text string found
c     next      input with first position of search string;
c                  output with position following text
c
c
      subroutine gettext (string,text,next)
      implicit none
      integer i,j,len,length,size,next
      integer first,last,extent,initial,final
      character*(*) string,text
c
c
c     get the length of input string and output text
c
      length = len(string(next:))
      size = len(text)
c
c     move through the string one character at a time,
c     searching for the first non-blank character
c
      first = 1
      last = 0
      initial = next
      final = next + length - 1
      do i = initial, final
         if (string(i:i) .gt. ' ') then
            first = i
            do j = i+1, final
               if (string(j:j) .le. ' ') then
                  last = j - 1
                  next = j
                  goto 10
               end if
            end do
         end if
      end do
   10 continue
c
c     trim the actual text if it is too long to return
c
      extent = next - first
      final = first + size - 1
      if (extent .gt. size)  last = final
c
c     transfer the text into the return string
c
      j = 0
      do i = first, last
         j = j + 1
         text(j:j) = string(i:i)
      end do
      do i = next, final
         j = j + 1
         text(j:j) = ' '
      end do
      return
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getword  --  extract first word from a string  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getword" searchs an input string for the first alphabetic
c     character (A-Z or a-z); the region from this first character
c     to the first blank space or comma is returned as a "word"; if
c     the actual word is too long, only the first part is returned
c
c     string    input character string to be searched
c     word      output with the first word in the string
c     next      input with first position of search string;
c                  output with position following word
c
c
      subroutine getword (string,word,next)
      implicit none
      integer i,j,len,length,size,next
      integer first,last,extent,initial,final
      character*1 letter
      character*(*) string,word
c
c
c     get the length of input string and output word
c
      length = len(string(next:))
      size = len(word)
c
c     move through the string one character at a time,
c     searching for the first alphabetic character
c
      first = 1
      last = 0
      initial = next
      final = next + length - 1
      do i = initial, final
         letter = string(i:i)
         if ((letter.ge.'A' .and. letter.le.'Z') .or.
     &       (letter.ge.'a' .and. letter.le.'z')) then
            first = i
            do j = i+1, final
               if (string(j:j).le.' ' .or. string(j:j).eq.',') then
                  last = j - 1
                  next = j
                  goto 10
               end if
            end do
         end if
      end do
   10 continue
c
c     trim the actual word if it is too long to return
c
      extent = next - first
      final = first + size - 1
      if (extent .gt. size)  last = final
c
c     transfer the word into the return string
c
      j = 0
      do i = first, last
         j = j + 1
         word(j:j) = string(i:i)
      end do
      do i = next, final
         j = j + 1
         word(j:j) = ' '
      end do
c
c     skip over the next character when it is a comma
c
      if (string(next:next) .eq. ',')  next = next + 1
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine getxyz  --  get Cartesian coordinate structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "getxyz" asks for a Cartesian coordinate file name,
c     then reads in the coordinates file
c
c
      subroutine getxyz
      implicit none
      include 'iounit.i'
      include 'output.i'
      integer ixyz
      integer freeunit
      character*60 xyzfile
      logical exist
      integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (xyzfile,exist)
      if (exist) then
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz')
         call version (xyzfile,'old')
         inquire (file=xyzfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      dowhile (.not. exist)
         write (iout,10)
   10    format (/,' Enter Cartesian Coordinate File Name :  ',$)
         read (input,20)  xyzfile
   20    format (a60)
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz')
         call version (xyzfile,'old')
         inquire (file=xyzfile,exist=exist)
      end do
c
c     first open and then read the Cartesian coordinates file
c
      coordtype = 'cartesian'
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      rewind (unit=ixyz)
      call readxyz (ixyz,iw)
      close (unit=ixyz)
      return
      end
c
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient  --  find energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine gradient (energy,derivs)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inter.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'virial.i'
      include 'warp.i'
      integer i,j
      real*8 energy,derivs(3,maxatm)
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      e14 = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     zero out each of the first derivative components
c
      do i = 1, n
         do j = 1, 3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deaa(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            deid(j,i) = 0.0d0
            deit(j,i) = 0.0d0
            det(j,i) = 0.0d0
            debt(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            de14(j,i) = 0.0d0
            dec(j,i) = 0.0d0
            decd(j,i) = 0.0d0
            ded(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
            der(j,i) = 0.0d0
            des(j,i) = 0.0d0
            deg(j,i) = 0.0d0
            dex(j,i) = 0.0d0
         end do
      end do
c
c     maintain any periodic boundary conditions
c
      if (use_bounds)  call bounds
c
c     zero out the virial and the intermolecular energy
c
      if (isobaric) then
         virx = 0.0d0
         viry = 0.0d0
         virz = 0.0d0
      end if
      einter = 0.0d0
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call piscf
c
c     call the local geometry energy and gradient routines
c
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_strtor)  call estrtor1
c     if (use_tortor)  call etortor1
c
c     call the van der Waals energy and gradient routines
c
      if (use_vdw) then
         if (use_lights) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj5
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck5
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb5
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal5
         else
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
         end if
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
      end if
c
c     call the electrostatic energy and gradient routines
c
      if (use_charge) then
         if (use_deform) then
            call echarge7
         else if (use_lights) then
            call echarge5
         else
            call echarge1
         end if
      end if
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
      if (use_mpole .or. use_polar)  call empole1
      if (use_rxnfld)  call erxnfld1
c
c     call any miscellaneous energy and gradient routines
c
      if (use_solv)  call esolv1
      if (use_geom)  call egeom1
      if (use_extra)  call extra1
c
c     sum up to get the total energy and first derivatives
c
      energy = eb + ea + eba + eub + eaa + eopb + eid
     &            + eit + et + ebt + ett + ev + e14 + ec
     &            + ecd + ed + em + ep + er + es + eg + ex
      do i = 1, n
         do j = 1, 3
            derivs(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
     &                       + deub(j,i) + deaa(j,i) + deopb(j,i)
     &                       + deid(j,i) + deit(j,i) + det(j,i)
     &                       + debt(j,i) + dett(j,i) + dev(j,i)
     &                       + de14(j,i) + dec(j,i) + decd(j,i)
     &                       + ded(j,i) + dem(j,i) + dep(j,i)
     &                       + der(j,i) + des(j,i) + deg(j,i)
     &                       + dex(j,i)
         end do
      end do
      return
      end
c
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine gradrgd  --  energy & gradient of rigid body  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "gradrgd" calls subroutines to calculate the potential energy
c     and first derivatives with respect to rigid body coordinates
c
c
      subroutine gradrgd (energy,derivs)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      include 'rigid.i'
      integer i,j,k
      integer init,stop
      real*8 energy,derivs(6,maxgrp)
      real*8 xcm,ycm,zcm
      real*8 xterm,yterm,zterm
      real*8 phi,cphi,sphi
      real*8 theta,ctheta,stheta
      real*8 ephi(3),etheta(3),epsi(3)
      real*8 g(3,maxatm),tau(3)
c
c
c     zero out the total of rigid body derivative components
c
      do i = 1, maxgrp
         do j = 1, 6
            derivs(j,i) = 0.0d0
         end do
      end do
c
c     calculate the energy and Cartesian first derivatives
c
      call gradient (energy,g)
c
c     compute the rigid body gradient components for each group
c
      do i = 1, ngrp
         init = igrp(1,i)
         stop = igrp(2,i)
         xcm = rbc(1,i)
         ycm = rbc(2,i)
         zcm = rbc(3,i)
         phi = rbc(4,i)
         theta = rbc(5,i)
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
c
c     get unit vectors along the phi, theta and psi rotation axes
c
         ephi(1) = 0.0d0
         ephi(2) = 0.0d0
         ephi(3) = 1.0d0
         etheta(1) = -sphi
         etheta(2) = cphi
         etheta(3) = 0.0d0
         epsi(1) = ctheta * cphi
         epsi(2) = ctheta * sphi
         epsi(3) = -stheta
c
c     first, get the rigid body gradients for translations
c
         do j = init, stop
            k = kgrp(j)
            derivs(1,i) = derivs(1,i) + g(1,k)
            derivs(2,i) = derivs(2,i) + g(2,k)
            derivs(3,i) = derivs(3,i) + g(3,k)
         end do
c
c     accumulate the moment arm along each axis of rotation
c
         do j = 1, 3
            tau(j) = 0.0d0
         end do
         do j = init, stop
            k = kgrp(j)
            xterm = x(k) - xcm
            yterm = y(k) - ycm
            zterm = z(k) - zcm
            tau(1) = tau(1) + yterm*g(3,k) - zterm*g(2,k)
            tau(2) = tau(2) + zterm*g(1,k) - xterm*g(3,k)
            tau(3) = tau(3) + xterm*g(2,k) - yterm*g(1,k)
         end do
c
c     now, set the rigid body gradients for rotations
c
         do j = 1, 3
            derivs(4,i) = derivs(4,i) + tau(j)*ephi(j)
            derivs(5,i) = derivs(5,i) + tau(j)*etheta(j)
            derivs(6,i) = derivs(6,i) + tau(j)*epsi(j)
         end do
      end do
      return
      end
c
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine gradrot  --  energy and torsional derivs  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "gradrot" calls subroutines to calculate the potential
c     energy and its torsional first derivatives
c
c
      subroutine gradrot (energy,derivs)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'domega.i'
      include 'omega.i'
      include 'potent.i'
      include 'rotate.i'
      integer i,j,k,base,partner
      real*8 xatom,yatom,zatom
      real*8 xdist,ydist,zdist
      real*8 xterm,yterm,zterm
      real*8 energy,g(3,maxatm)
      real*8 norm,derivs(maxrot)
c
c
c     zero out individual components of torsional derivatives
c
      do i = 1, nomega
         teb(i) = 0.0d0
         tea(i) = 0.0d0
         teba(i) = 0.0d0
         teub(i) = 0.0d0
         teaa(i) = 0.0d0
         teopb(i) = 0.0d0
         teid(i) = 0.0d0
         teit(i) = 0.0d0
         tet(i) = 0.0d0
         tebt(i) = 0.0d0
         tett(i) = 0.0d0
         tev(i) = 0.0d0
         te14(i) = 0.0d0
         tec(i) = 0.0d0
         tecd(i) = 0.0d0
         ted(i) = 0.0d0
         tem(i) = 0.0d0
         tep(i) = 0.0d0
         ter(i) = 0.0d0
         tes(i) = 0.0d0
         teg(i) = 0.0d0
         tex(i) = 0.0d0
      end do
c
c     calculate the energy and Cartesian first derivatives
c
      call gradient (energy,g)
c
c     transform Cartesian derivatives to torsional space
c
      do i = 1, nomega
         base = iomega(1,i)
         partner = iomega(2,i)
         call rotlist (base,partner)
         xdist = x(base) - x(partner)
         ydist = y(base) - y(partner)
         zdist = z(base) - z(partner)
         norm = sqrt(xdist**2 + ydist**2 + zdist**2)
         xdist = xdist / norm
         ydist = ydist / norm
         zdist = zdist / norm
         do j = 1, nrot
            k = rot(j)
            xatom = x(k) - x(base)
            yatom = y(k) - y(base)
            zatom = z(k) - z(base)
            xterm = ydist*zatom - zdist*yatom
            yterm = zdist*xatom - xdist*zatom
            zterm = xdist*yatom - ydist*xatom
            teb(i) = teb(i) + deb(1,k)*xterm + deb(2,k)*yterm
     &                              + deb(3,k)*zterm
            tea(i) = tea(i) + dea(1,k)*xterm + dea(2,k)*yterm
     &                              + dea(3,k)*zterm
            teba(i) = teba(i) + deba(1,k)*xterm + deba(2,k)*yterm
     &                              + deba(3,k)*zterm
            teub(i) = teub(i) + deub(1,k)*xterm + deub(2,k)*yterm
     &                              + deub(3,k)*zterm
            teaa(i) = teaa(i) + deaa(1,k)*xterm + deaa(2,k)*yterm
     &                              + deaa(3,k)*zterm
            teopb(i) = teopb(i) + deopb(1,k)*xterm + deopb(2,k)*yterm
     &                              + deopb(3,k)*zterm
            teid(i) = teid(i) + deid(1,k)*xterm + deid(2,k)*yterm
     &                              + deid(3,k)*zterm
            teit(i) = teit(i) + deit(1,k)*xterm + deit(2,k)*yterm
     &                              + deit(3,k)*zterm
            tet(i) = tet(i) + det(1,k)*xterm + det(2,k)*yterm
     &                              + det(3,k)*zterm
            tebt(i) = tebt(i) + debt(1,k)*xterm + debt(2,k)*yterm
     &                              + debt(3,k)*zterm
            tett(i) = tett(i) + dett(1,k)*xterm + dett(2,k)*yterm
     &                              + dett(3,k)*zterm
            tev(i) = tev(i) + dev(1,k)*xterm + dev(2,k)*yterm
     &                              + dev(3,k)*zterm
            te14(i) = te14(i) + de14(1,k)*xterm + de14(2,k)*yterm
     &                              + de14(3,k)*zterm
            tec(i) = tec(i) + dec(1,k)*xterm + dec(2,k)*yterm
     &                              + dec(3,k)*zterm
            tecd(i) = tecd(i) + decd(1,k)*xterm + decd(2,k)*yterm
     &                              + decd(3,k)*zterm
            ted(i) = ted(i) + ded(1,k)*xterm + ded(2,k)*yterm
     &                              + ded(3,k)*zterm
            tem(i) = tem(i) + dem(1,k)*xterm + dem(2,k)*yterm
     &                              + dem(3,k)*zterm
            tep(i) = tep(i) + dep(1,k)*xterm + dep(2,k)*yterm
     &                              + dep(3,k)*zterm
            ter(i) = ter(i) + der(1,k)*xterm + der(2,k)*yterm
     &                              + der(3,k)*zterm
            tes(i) = tes(i) + des(1,k)*xterm + des(2,k)*yterm
     &                              + des(3,k)*zterm
            teg(i) = teg(i) + deg(1,k)*xterm + deg(2,k)*yterm
     &                              + deg(3,k)*zterm
            tex(i) = tex(i) + dex(1,k)*xterm + dex(2,k)*yterm
     &                              + dex(3,k)*zterm
         end do
      end do
c
c     sum up to give the total torsional first derivatives
c
      do i = 1, nomega
         derivs(i) = teb(i) + tea(i) + teba(i) + teub(i) + teaa(i)
     &                  + teopb(i) + teid(i) + teit(i) + tet(i)
     &                  + tebt(i) + tett(i) + tev(i) + te14(i)
     &                  + tec(i) + tecd(i) + ted(i) + tem(i)
     &                  + tep(i) + ter(i) + tes(i) + teg(i) + tex(i)
      end do
      return
      end
c
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine groups  --  group membership of atom set  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "groups" tests a set of atoms to see if all are members
c     of a single atom group or a pair of atom groups; if so,
c     then the correct intra- or intergroup weight is assigned
c
c
      subroutine groups (useset,weight,nset,ia,ib,ic,id,ie)
      implicit none
      include 'sizes.i'
      include 'group.i'
      integer ia,ib,ic,id,ie,nset
      integer iga,igb,igc,igd,ige
      integer gmax,gmin
      real*8 weight
      logical useset
c
c
c     check group membership for a set containing one atom
c
      if (nset .eq. 1) then
         iga = grplist(ia)
         weight = wgrp(iga,iga)
c
c     check group membership for a set containing two atoms
c
      else if (nset .eq. 2) then
         iga = grplist(ia)
         igb = grplist(ib)
         weight = wgrp(iga,igb)
c
c     check group membership for a set containing three atoms
c
      else if (nset .eq. 3) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         if (iga.eq.igb .or. igb.eq.igc) then
            weight = wgrp(iga,igc)
         else if (iga .eq. igc) then
            weight = wgrp(iga,igb)
         else
            weight = 0.0d0
         end if
c
c     check group membership for a set containing four atoms
c
      else if (nset .eq. 4) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         igd = grplist(id)
         gmin = min(iga,igb,igc,igd)
         gmax = max(iga,igb,igc,igd)
         if ((iga.ne.gmin .and. iga.ne.gmax) .or.
     &       (igb.ne.gmin .and. igb.ne.gmax) .or.
     &       (igc.ne.gmin .and. igc.ne.gmax) .or.
     &       (igd.ne.gmin .and. igd.ne.gmax)) then
            weight = 0.0d0
         else
            weight = wgrp(gmin,gmax)
         end if
c
c     check group membership for a set containing five atoms
c
      else if (nset .eq. 5) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         igd = grplist(id)
         ige = grplist(ie)
         gmin = min(iga,igb,igc,igd,ige)
         gmax = max(iga,igb,igc,igd,ige)
         if ((iga.ne.gmin .and. iga.ne.gmax) .or.
     &       (igb.ne.gmin .and. igb.ne.gmax) .or.
     &       (igc.ne.gmin .and. igc.ne.gmax) .or.
     &       (igd.ne.gmin .and. igd.ne.gmax) .or.
     &       (ige.ne.gmin .and. ige.ne.gmax)) then
            weight = 0.0d0
         else
            weight = wgrp(gmin,gmax)
         end if
      end if
c
c     interaction will be used if its group has nonzero weight
c
      if (weight .eq. 0.0d0) then
         useset = .false.
      else
         useset = .true.
      end if
      return
      end
c
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine gyrate  --  compute the radius of gyration  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "gyrate" computes the radius of gyration of a molecular
c     system from its atomic coordinates
c
c
      subroutine gyrate (rg)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i
      real*8 rg,xc,yc,zc
c
c
c     find the centroid of the atomic coordinates
c
      xc = 0.0d0
      yc = 0.0d0
      zc = 0.0d0
      do i = 1, n
         xc = xc + x(i)
         yc = yc + y(i)
         zc = zc + z(i)
      end do
      xc = xc / dble(n)
      yc = yc / dble(n)
      zc = zc / dble(n)
c
c     compute and print out the radius of gyration
c
      rg = 0.0d0
      do i = 1, n
         rg = rg + (x(i)-xc)**2 + (y(i)-yc)**2 + (z(i)-zc)**2
      end do
      rg = sqrt(rg/dble(n))
      return
      end
c
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine hessian  --  atom-by-atom Hessian elements  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "hessian" calls subroutines to calculate the Hessian elements
c     for each atom in turn with respect to Cartesian coordinates
c
c
      subroutine hessian (h,hinit,hstop,hindex,hdiag)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'hescut.i'
      include 'hessn.i'
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer i,j,k,ii,nhess
      integer hinit(3,maxatm),hstop(3,maxatm)
c
      LOGICAL GOPARR,DSKWRK,MASWRK,frzmm
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW,maxhess
      CHARACTER*6 CTMODE, CTMETH
      real*8 OPTPRG,MINMET,GRDMIN
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /NWTOPT/ OPTPRG,MINMET,GRDMIN,maxhess,frzmm,CTMODE,CTMETH
c
      integer hindex(maxhess)
      real*8 percent,filled,rdn,hmax
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      real*8 hdiag(3,maxatm),h(maxhess)
      logical keep(maxatm)
c
      logical first
      data first/.true./
      save first
c
c     zero out total number of indexed Hessian elements
c
      nhess = 0
      do i = 1, n
         do j = 1, 3
            hinit(j,i) = 1
            hstop(j,i) = 0
            hdiag(j,i) = 0.0d0
         end do
      end do
c
c     maintain any periodic boundary conditions
c
      if (use_bounds)  call bounds
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call piscf
c
c     set Born radii for use with GB/SA solvation
c
      if (use_gbsa)  call born
c
c     calculate the "reduced" atomic coordinates
c
      if (use_vdw) then
         do i = 1, n
            ii = ired(i)
            rdn = kred(i)
            xred(i) = rdn*(x(i)-x(ii)) + x(ii)
            yred(i) = rdn*(y(i)-y(ii)) + y(ii)
            zred(i) = rdn*(z(i)-z(ii)) + z(ii)
         end do
      end if
c
c     zero out the Hessian elements for the current atom
c
      do i = 1, n
         if (use(i)) then
            do k = 1, n
               do j = 1, 3
                  hessx(j,k) = 0.0d0
                  hessy(j,k) = 0.0d0
                  hessz(j,k) = 0.0d0
               end do
            end do
c
c     call the local geometry Hessian component routines
c
            if (use_bond)  call ebond2 (i)
            if (use_angle)  call eangle2 (i)
            if (use_strbnd)  call estrbnd2 (i)
            if (use_urey)  call eurey2 (i)
            if (use_angang)  call eangang2 (i)
            if (use_opbend)  call eopbend2 (i)
            if (use_improp)  call eimprop2 (i)
            if (use_imptor)  call eimptor2 (i)
            if (use_tors)  call etors2 (i)
            if (use_strtor)  call estrtor2 (i)
c           if (use_tortor)  call etortor2 (i)
c
c     call the van der Waals Hessian component routines
c
            if (use_vdw) then
               if (vdwtyp .eq. 'LENNARD-JONES') then
                  call elj2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'BUCKINGHAM') then
                  call ebuck2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'MM3-HBOND') then
                  call emm3hb2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
                  call ehal2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'GAUSSIAN') then
                  call egauss2 (i,xred,yred,zred)
               end if
            end if
c
c     call the electrostatic Hessian component routines
c
            if (use_charge) then
               if (use_deform) then
                  call echarge8 (i)
               else
                  call echarge2 (i)
               end if
            end if
            if (use_chgdpl)  call echgdpl2 (i)
            if (use_dipole)  call edipole2 (i)
            if (use_mpole .or. use_polar)   call empole2 (i)
            if (use_rxnfld)   call erxnfld2 (i)
c
c     call any miscellaneous Hessian component routines
c
            if (use_solv)  call esolv2 (i)
            if (use_geom)  call egeom2 (i)
            if (use_extra)  call extra2 (i)
c
c     set the diagonal Hessian matrix elements
c
            hdiag(1,i) = hdiag(1,i) + hessx(1,i)
            hdiag(2,i) = hdiag(2,i) + hessy(2,i)
            hdiag(3,i) = hdiag(3,i) + hessz(3,i)
c
c     search each 3x3 block to see which blocks will be kept
c
            do k = i+1, n
               keep(k) = .false.
               if (use(k)) then
                  hmax = max(abs(hessx(1,k)),abs(hessx(2,k)),
     &                       abs(hessx(3,k)),abs(hessy(1,k)),
     &                       abs(hessy(2,k)),abs(hessy(3,k)),
     &                       abs(hessz(1,k)),abs(hessz(2,k)),
     &                       abs(hessz(3,k)))
                  if (hmax .ge. hesscut)  keep(k) = .true.
               end if
            end do
c
c     copy selected off-diagonal Hessian elements for current
c     atom into an indexed master list of Hessian elements;
c     if any elements of 3x3 block are kept, keep them all
c
            hinit(1,i) = nhess + 1
            do j = 2, 3
               nhess = nhess + 1
               hindex(nhess) = 3*i + j - 3
               h(nhess) = hessx(j,i)
            end do
            do k = i+1, n
               if (keep(k)) then
                  do j = 1, 3
                     nhess = nhess + 1
                     hindex(nhess) = 3*k + j - 3
                     h(nhess) = hessx(j,k)
                  end do
               end if
            end do
            hstop(1,i) = nhess
            hinit(2,i) = nhess + 1
            nhess = nhess + 1
            hindex(nhess) = 3*i
            h(nhess) = hessy(3,i)
            do k = i+1, n
               if (keep(k)) then
                  do j = 1, 3
                     nhess = nhess + 1
                     hindex(nhess) = 3*k + j - 3
                     h(nhess) = hessy(j,k)
                  end do
               end if
            end do
            hstop(2,i) = nhess
            hinit(3,i) = nhess + 1
            do k = i+1, n
               if (keep(k)) then
                  do j = 1, 3
                     nhess = nhess + 1
                     hindex(nhess) = 3*k + j - 3
                     h(nhess) = hessz(j,k)
                  end do
               end if
            end do
            hstop(3,i) = nhess
c
c     check for storage of too many Hessian elements
c
            if (nhess .gt. maxhess) then
               if (maswrk) write (iout,10)  i,n,nhess,maxhess,hesscut
   10          format (/,' Current atom :',i7,6x,'Total atoms :',i7,
     &                 /,' Current required Hessian storage : ',i12,
     &                 /,' Maximum allowed Hessian storage :  ',i12,
     &                 /,' Minimum significant Hessian value :',f12.5)
               if (maswrk) write (iout,20)
   20          format (/,' HESSIAN  --  Increase MAXHES in $TOPTMZ,',
     &                   ' and/or HESSCUT.')
               call fatal
            end if
         end if
      end do
c
c     print message telling how much storage was finally used
c
      if (first  .or.  verbose) then
         first=.false.
         percent = 100.0d0 * dble(nhess)/dble(3*n*(3*n-1)/2)
         filled = 100.0d0 * dble(nhess)/dble(maxhess)
         if (maswrk) write (iout,30)  nhess,maxhess,percent,filled
   30    format (1x,'HESSIAN --',i11,' elements used, MAXHES=',i11,
     &              ' in $TOPTMZ input.'/
     &           1x,f8.2,'% Off-diag H',f8.2,'% of available Storage')
      end if
      return
      end
c
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine hessrgd  --  rigid body Hessian elements  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "hessrgd" computes the numerical Hessian elements with
c     respect to rigid body coordinates; either the full matrix
c     or just the diagonal can be calculated; the full matrix
c     needs 6*ngroup+1 gradient evaluations while the diagonal
c     requires just two gradient calls
c
c
      subroutine hessrgd (mode,hrigid)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      include 'rigid.i'
      integer maxrgd
      parameter (maxrgd=6*maxgrp)
      integer i,j,k,m,nvar
      real*8 eps,old(6,maxgrp)
      real*8 e,g(6,maxgrp),g0(6,maxgrp)
      real*8 hrigid(maxrgd,maxrgd)
      character*4 mode
c
c
c     calculate base values for the rigid body gradient
c
      eps = 0.0001d0
      call gradrgd (e,g0)
c     
c     compute one-sided numerical Hessian from gradient values;
c     set off-diagonal elements to the average symmetric value
c
      if (mode .eq. 'full') then
         nvar = 6 * ngrp
         do i = 1, nvar
            j = (i-1)/6 + 1
            k = mod(i-1,6) + 1
            old(k,j) = rbc(k,j)
            rbc(k,j) = rbc(k,j) + eps
            call rigidxyz
            call gradrgd (e,g)
            do m = 1, nvar
               hrigid(m,i) = (g(k,j) - g0(k,j)) / eps
            end do
            do m = 1, i-1
               hrigid(m,i) = 0.5d0 * (hrigid(m,i)+hrigid(i,m))
               hrigid(i,m) = hrigid(m,i)
            end do
            rbc(k,j) = old(k,j)
         end do
c
c     compute numerical Hessian diagonal from gradient values
c
      else if (mode .eq. 'diag') then
         do i = 1, ngrp
            do j = 1, 6
               old(j,i) = rbc(j,i)
               rbc(j,i) = rbc(j,i) + eps
            end do
         end do
         call rigidxyz
         call gradrot (e,g)
         nvar = 0
         do i = 1, ngrp
            do j = 1, 6
               nvar = nvar + 1
               hrigid(nvar,nvar) = (g(j,i)-g0(j,i)) / eps
               rbc(j,i) = old(j,i)
            end do
         end do
      end if
c
c     restore the Cartesian coordinates to original values
c
      call rigidxyz
      return
      end
c
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine hessrot  --  torsional Hessian elements  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "hessrot" computes the numerical Hessian elements with
c     respect to torsional angles; either the full matrix or
c     just the diagonal can be calculated; the full matrix
c     needs nomega+1 gradient evaluations while the diagonal
c     requires just two gradient calls
c
c
      subroutine hessrot (mode,hrot)
      implicit none
      include 'sizes.i'
      include 'omega.i'
      include 'math.i'
      include 'zcoord.i'
      integer i,j,line
      real*8 eps,old(maxrot)
      real*8 e,g(maxrot),g0(maxrot)
      real*8 hrot(maxrot,maxrot)
      character*4 mode
c
c
c     calculate base values for the torsional gradient
c
      eps = 0.0001d0
      call gradrot (e,g0)
c
c     compute one-sided numerical Hessian from gradient values;
c     set off-diagonal elements to the average symmetric value
c
      if (mode .eq. 'full') then
         do i = 1, nomega
            line = zline(i)
            old(i) = ztors(line)
            ztors(line) = ztors(line) + radian*eps
            call makexyz
            call gradrot (e,g)
            ztors(line) = old(i)
            do j = 1, nomega
               hrot(j,i) = (g(j)-g0(j)) / eps
            end do
            do j = 1, i-1
               hrot(j,i) = 0.5d0 * (hrot(j,i)+hrot(i,j))
               hrot(i,j) = hrot(j,i)
            end do
         end do
c
c     compute numerical Hessian diagonal from gradient values
c
      else if (mode .eq. 'diag') then
         do i = 1, nomega
            line = zline(i)
            old(i) = ztors(line)
            ztors(line) = ztors(line) + radian*eps
         end do
         call makexyz
         call gradrot (e,g)
         do i = 1, nomega
            hrot(i,i) = (g(i)-g0(i)) / eps
            line = zline(i)
            ztors(line) = old(i)
         end do
      end if
c
c     restore the Cartesian coordinates to original values
c
      call makexyz
      return
      end
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine hybrid  --  set parameters for hybrid system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "hybrid" constructs the hybrid hamiltonian for a specified
c     initial state, final state and mutation parameter "lambda"
c
c
      subroutine hybrid
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'mutant.i'
c
c
c     set the potential energy parameters for hybrid atoms
c
      if (nhybrid .ne. 0) then
         write (iout,10)  lambda
   10    format (/,' Lambda Coupling Parameter for FEP :',f12.3)
         call hatom
         call hbond
         call hangle
         call hstrbnd
         call himptor
         call htors
         call hstrtor
         call hvdw
         call hcharge
         call hdipole
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine hatom  --  assign hybrid atom parameters  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "hatom" assigns a new atom type to each hybrid site
c
c
      subroutine hatom
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'mutant.i'
      integer i,k,it,it0,it1,ntype
c
c
c     find the total number of atom types currently used;
c     exclude the "HYB" types so that they can be reused
c
      do i = 1, maxtyp
         if (symbol(i).eq.'   ' .or. symbol(i).eq.'HYB') then
            ntype = i - 1
            goto 10
         end if
      end do
   10 continue
c
c     stop if there are too many atom types required
c
      if (maxtyp .lt. ntype+nhybrid) then
         abort = .true.
         write (iout,20)
   20    format (' HATOM  --  Too many Sites to be Altered;',
     &              ' Increase MAXTYP')
      end if
c
c     create a new atom type for each of the hybrid atoms
c
      do i = 1, nhybrid
         k = ihybrid(i)
         it = ntype + i
         it0 = type0(i)
         it1 = type1(i)
         symbol(it) = 'HYB'
         atmnum(it) = 0
         weight(it) = lambda*weight(it1) + (1.0d0-lambda)*weight(it0)
         ligand(it) = 0
         describe(it) = 'Hybrid Atom Type    '
         type(k) = it
         name(k) = symbol(it)
         atomic(k) = atmnum(it)
         mass(k) = weight(it)
         valence(k) = ligand(it)
         story(k) = describe(it)
      end do
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine hbond  --  find hybrid bond parameters  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "hbond" constructs the hybrid bonded interaction
c     parameters given an initial state, final state and
c     mutation parameter "lambda"
c
c
      subroutine hbond
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'iounit.i'
      include 'inform.i'
      include 'kbonds.i'
      include 'mutant.i'
      integer i,j,k,ia,ib
      integer ita,itb,size
      real*8 bk0,bk1,bl0,bl1
      character*3 pa,pb
      character*6 pt
      logical header
c
c
c     assign the hybrid parameters for individual bonds
c
      header = .true.
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         if (alter(ia) .or. alter(ib)) then
            ita = class(ia)
            itb = class(ib)
c
c     find the bond parameters for the initial state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
            end do
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            bk0 = 0.0d0
            bl0 = 0.0d0
            do j = 1, maxnb
               if (kb(j) .eq. pt) then
                  bk0 = fcon(j)
                  bl0 = blen(j)
                  goto 10
               end if
            end do
   10       continue
c
c     find the bond parameters for the final state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = class1(j)
               if (k .eq. ib)  itb = class1(j)
            end do
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            bk1 = 0.0d0
            bl1 = 0.0d0
            do j = 1, maxnb
               if (kb(j) .eq. pt) then
                  bk1 = fcon(j)
                  bl1 = blen(j)
                  goto 20
               end if
            end do
   20       continue
c
c     form the hybrid parameters for the current bond
c
            if (bl0 .eq. 0.0d0)  bl0 = bl1
            if (bl1 .eq. 0.0d0)  bl1 = bl0
            bk(i) = lambda*bk1 + (1.0d0-lambda)*bk0
            bl(i) = lambda*bl1 + (1.0d0-lambda)*bl0
            if (verbose) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Hybrid Bond Stretching Parameters :',
     &                    //,6x,'Atom Numbers',9x,'KS',7x,'Length',/)
               end if
               write (iout,40)  ia,ib,bk(i),bl(i)
   40          format (6x,2i5,f14.3,f12.4)
            end if
         end if
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine hangle  --  find hybrid angle parameters  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "hangle" constructs the hybrid angle bending parameters
c     given an initial state, final state and "lambda" value
c
c
      subroutine hangle
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'inform.i'
      include 'kangs.i'
      include 'mutant.i'
      integer i,j,k,ia,ib,ic
      integer ita,itb,itc,size
      real*8 acon0,acon1,anat0,anat1
      character*3 pa,pb,pc
      character*9 pt
      logical header
c
c
c     assign the hybrid parameters for individual angles
c
      header = .true.
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (alter(ia) .or. alter(ib) .or. alter(ic)) then
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
c
c     find the angle parameters for the initial state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
            end do
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            acon0 = 0.0d0
            anat0 = 0.0d0
            do j = 1, maxna
               if (ka(j) .eq. pt) then
                  acon0 = con(j)
                  anat0 = ang(1,j)
                  goto 10
               end if
            end do
   10       continue
c
c     find the angle parameters for the final state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = class1(j)
               if (k .eq. ib)  itb = class1(j)
               if (k .eq. ic)  itc = class1(j)
            end do
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            acon1 = 0.0d0
            anat1 = 0.0d0
            do j = 1, maxna
               if (ka(j) .eq. pt) then
                  acon1 = con(j)
                  anat1 = ang(1,j)
                  goto 20
               end if
            end do
   20       continue
c
c     form the hybrid parameters for the current angle
c
            if (anat0 .eq. 0.0d0)  anat0 = anat1
            if (anat1 .eq. 0.0d0)  anat1 = anat0
            acon(i) = lambda*acon1 + (1.0d0-lambda)*acon0
            anat(i) = lambda*anat1 + (1.0d0-lambda)*anat0
            if (verbose) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Hybrid Angle Bending Parameters :',
     &                    //,6x,'Atom Numbers',9x,'KB',8x,'Angle',/)
               end if
               write (iout,40)  ia,ib,ic,acon(i),anat(i)
   40          format (3x,3i5,2f12.3)
            end if
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine hstrbnd  --  hybrid stretch-bend parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "hstrbnd" constructs the hybrid stretch-bend parameters
c     given an initial state, final state and "lambda" value
c
c
      subroutine hstrbnd
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'inform.i'
      include 'katoms.i'
      include 'kstbnd.i'
      include 'mutant.i'
      include 'strbnd.i'
      integer i,j,k,ia,ib,ic
      integer ita,itb,itc,nsb
      integer nba,nba0,nba1
      integer nbc,nbc0,nbc1
      real*8 ksb0,ksb1
      logical header
c
c
c     assign hybrid parameters for the stretch-bend sites
c
      header = .true.
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (alter(ia) .or. alter(ib) .or. alter(ic)) then
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            do j = 1, n12(ib)
               if (i12(j,ib) .eq. ia)  nba = bndlist(j,ib)
               if (i12(j,ib) .eq. ic)  nbc = bndlist(j,ib)
            end do
c
c     find the stretch-bend parameters for the initial state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
            end do
            nsb = 1
            nba0 = nba
            nbc0 = nbc
            if (atmnum(ita) .le. 1) then
               nsb = nsb + 1
               if (stbn(3,itb) .eq. 0.0d0)  nba0 = 0
            end if
            if (atmnum(itc) .le. 1) then
               nsb = nsb + 1
               if (stbn(3,itb) .eq. 0.0d0)  nbc0 = 0
            end if
            ksb0 = stbn(nsb,itb)
c
c     find the stretch-bend parameters for the final state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = class1(j)
               if (k .eq. ib)  itb = class1(j)
               if (k .eq. ic)  itc = class1(j)
            end do
            nsb = 1
            nba1 = nba
            nbc1 = nbc
            if (atmnum(ita) .le. 1) then
               nsb = nsb + 1
               if (stbn(3,itb) .eq. 0.0d0)  nba1 = 0
            end if
            if (atmnum(itc) .le. 1) then
               nsb = nsb + 1
               if (stbn(3,itb) .eq. 0.0d0)  nbc1 = 0
            end if
            ksb1 = stbn(nsb,itb)
c
c     form hybrid parameters for the current stretch-bend
c
            ksb(i) = lambda*ksb1 + (1.0d0-lambda)*ksb0
            if (ksb(i) .eq. 0.0d0) then
               if (isb(1,i) .ne. 0) then
                  nstrbnd = nstrbnd - 1
                  isb(1,i) = 0
               end if
            else
               if (isb(1,i) .ne. i) then
                  nstrbnd = nstrbnd + 1
                  isb(1,i) = i
               end if
               isb(2,i) = max(nba0,nba1)
               isb(3,i) = max(nbc0,nbc1)
               if (verbose) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Hybrid Stretch-Bend Parameters :',
     &                       //,6x,'Atom Numbers',9x,'KSB',/)
                  end if
                  write (iout,20)  ia,ib,ic,ksb(i)
   20             format (3x,3i5,f12.3)
               end if
            end if
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine himptor  --  find hybrid improper torsions  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "himptor" constructs the hybrid improper torsional
c     parameters given an initial state, final state and
c     mutation parameter "lambda"
c
c
      subroutine himptor
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'inform.i'
      include 'imptor.i'
      include 'kitors.i'
      include 'mutant.i'
      integer i,j,k,ia,ib,ic,id
      integer ita,itb,itc,itd
      integer it1,it2,it3,size
      real*8 v1_0,v2_0,v3_0,s1_0,s2_0,s3_0
      real*8 v1_1,v2_1,v3_1,s1_1,s2_1,s3_1
      character*3 pb,p1,p2,p3
      character*12 pt
      logical header,used
c
c
c     construct hybrid improper torsion parameters
c
      header = .true.
      do i = 1, n
         if (n12(i) .eq. 3) then
            used = .false.
            ia = i12(1,i)
            ib = i
            ic = i12(2,i)
            id = i12(3,i)
            if (alter(ia) .or. alter(ib) .or.
     &          alter(ic) .or. alter(id)) then
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
c
c     find improper torsion parameters for the initial state
c
               do j = 1, nhybrid
                  k = ihybrid(j)
                  if (k .eq. ia)  ita = class0(j)
                  if (k .eq. ib)  itb = class0(j)
                  if (k .eq. ic)  itc = class0(j)
                  if (k .eq. id)  itd = class0(j)
               end do
               it1 = min(ita,itc,itd)
               it2 = max(min(ita,itc),min(ita,itd),min(itc,itd))
               it3 = max(ita,itc,itd)
               size = 3
               call numeral (itb,pb,size)
               call numeral (it1,p1,size)
               call numeral (it2,p2,size)
               call numeral (it3,p3,size)
               pt = pb//p1//p2//p3
               v1_0 = 0.0d0
               s1_0 = 0.0d0
               v2_0 = 0.0d0
               s2_0 = 0.0d0
               v3_0 = 0.0d0
               s3_0 = 0.0d0
               do j = 1, maxnti
                  if (kti(j) .eq. pt) then
                     used = .true.
                     v1_0 = ti1(1,j)
                     s1_0 = ti1(2,j)
                     v2_0 = ti2(1,j)
                     s2_0 = ti2(2,j)
                     v3_0 = ti3(1,j)
                     s3_0 = ti3(2,j)
                     goto 10
                  end if
               end do
   10          continue
c
c     find improper torsion parameters for the final state
c
               do j = 1, nhybrid
                  k = ihybrid(j)
                  if (k .eq. ia)  ita = class1(j)
                  if (k .eq. ib)  itb = class1(j)
                  if (k .eq. ic)  itc = class1(j)
                  if (k .eq. id)  itd = class1(j)
               end do
               it1 = min(ita,itc,itd)
               it2 = max(min(ita,itc),min(ita,itd),min(itc,itd))
               it3 = max(ita,itc,itd)
               size = 3
               call numeral (itb,pb,size)
               call numeral (it1,p1,size)
               call numeral (it2,p2,size)
               call numeral (it3,p3,size)
               pt = pb//p1//p2//p3
               v1_1 = 0.0d0
               s1_1 = 0.0d0
               v2_1 = 0.0d0
               s2_1 = 0.0d0
               v3_1 = 0.0d0
               s3_1 = 0.0d0
               do j = 1, maxnti
                  if (kti(j) .eq. pt) then
                     used = .true.
                     v1_1 = ti1(1,j)
                     s1_1 = ti1(2,j)
                     v2_1 = ti2(1,j)
                     s2_1 = ti2(2,j)
                     v3_1 = ti3(1,j)
                     s3_1 = ti3(2,j)
                     goto 20
                  end if
               end do
   20          continue
            end if
c
c     form hybrid parameters for the current improper torsion
c
            if (used) then
               do j = 1, nitors
                  if (iitors(2,j) .eq. ib) then
                     k = j
                     goto 30
                  end if
               end do
               nitors = nitors + 1
               k = nitors
   30          continue
               iitors(1,k) = ia
               iitors(2,k) = ib
               iitors(3,k) = ic
               iitors(4,k) = id
               if (s1_0 .eq. 0.0d0)  s1_0 = s1_1
               if (s2_0 .eq. 0.0d0)  s2_0 = s2_1
               if (s3_0 .eq. 0.0d0)  s3_0 = s3_1
               if (s1_1 .eq. 0.0d0)  s1_1 = s1_0
               if (s2_1 .eq. 0.0d0)  s2_1 = s2_0
               if (s3_1 .eq. 0.0d0)  s3_1 = s3_0
               itors1(1,k) = lambda*v1_1 + (1.0d0-lambda)*v1_0
               itors1(2,k) = lambda*s1_1 + (1.0d0-lambda)*s1_0
               itors2(1,k) = lambda*v2_1 + (1.0d0-lambda)*v2_0
               itors2(2,k) = lambda*s2_1 + (1.0d0-lambda)*s2_0
               itors3(1,k) = lambda*v3_1 + (1.0d0-lambda)*v3_0
               itors3(2,k) = lambda*s3_1 + (1.0d0-lambda)*s3_0
            end if
            if (verbose .and. used) then
               if (header) then
                  header = .false.
                  write (iout,40)
   40             format (/,' Hybrid Improper Torsional Parameters :',
     &                    //,6x,'Atom Numbers',16x,'KIT1',13x,'KIT2',
     &                       13x,'KIT3',/)
               end if
               write (iout,50)  ia,ib,ic,id,itors1(1,k),itors1(2,k),
     &                          itors2(1,k),itors2(2,k),
     &                          itors3(1,k),itors3(2,k)
   50          format (1x,4i5,4x,3(f10.4,f7.3))
            end if
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine htors  --  find hybrid torsion parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "htors" constructs the hybrid torsional parameters for
c     a given initial state, final state, and "lambda" value
c
c
      subroutine htors
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'ktorsn.i'
      include 'mutant.i'
      include 'tors.i'
      integer i,j,k,ia,ib,ic,id
      integer ita,itb,itc,itd,size
      real*8  v1_0,s1_0,v2_0,s2_0,v3_0,s3_0
      real*8  v1_1,s1_1,v2_1,s2_1,v3_1,s3_1
      character*3 pa,pb,pc,pd
      character*12 pt
      logical header
c
c
c     construct hybrid torsional parameters
c
      header = .true.
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         if (alter(ia) .or. alter(ib) .or.
     &       alter(ic) .or. alter(id)) then
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
c
c     find the torsion parameters for the initial state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
               if (k .eq. id)  itd = class0(j)
            end do
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            if (itb .lt. itc) then
               pt = pa//pb//pc//pd
            else if (itc .lt. itb) then
               pt = pd//pc//pb//pa
            else if (ita .le. itd) then
               pt = pa//pb//pc//pd
            else if (itd .lt. ita) then
               pt = pd//pc//pb//pa
            end if
            v1_0 = 0.0d0
            s1_0 = 0.0d0
            v2_0 = 0.0d0
            s2_0 = 0.0d0
            v3_0 = 0.0d0
            s3_0 = 0.0d0
            do j = 1, maxnt
               if (kt(j) .eq. pt) then
                  v1_0 = t1(1,j)
                  s1_0 = t1(2,j)
                  v2_0 = t2(1,j)
                  s2_0 = t2(2,j)
                  v3_0 = t3(1,j)
                  s3_0 = t3(2,j)
                  goto 10
               end if
            end do
   10       continue
c
c     find the torsion parameters for the final state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = class1(j)
               if (k .eq. ib)  itb = class1(j)
               if (k .eq. ic)  itc = class1(j)
               if (k .eq. id)  itd = class1(j)
            end do
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            if (itb .lt. itc) then
               pt = pa//pb//pc//pd
            else if (itc .lt. itb) then
               pt = pd//pc//pb//pa
            else if (ita .le. itd) then
               pt = pa//pb//pc//pd
            else if (itd .lt. ita) then
               pt = pd//pc//pb//pa
            end if
            v1_1 = 0.0d0
            s1_1 = 0.0d0
            v2_1 = 0.0d0
            s2_1 = 0.0d0
            v3_1 = 0.0d0
            s3_1 = 0.0d0
            do j = 1, maxnt
               if (kt(j) .eq. pt) then
                  v1_1 = t1(1,j)
                  s1_1 = t1(2,j)
                  v2_1 = t2(1,j)
                  s2_1 = t2(2,j)
                  v3_1 = t3(1,j)
                  s3_1 = t3(2,j)
                  goto 20
               end if
            end do
   20       continue
c
c     form the hybrid parameters for the current torsion
c
            if (s1_0 .eq. 0.0d0)  s1_0 = s1_1
            if (s2_0 .eq. 0.0d0)  s2_0 = s2_1
            if (s3_0 .eq. 0.0d0)  s3_0 = s3_1
            if (s1_1 .eq. 0.0d0)  s1_1 = s1_0
            if (s2_1 .eq. 0.0d0)  s2_1 = s2_0
            if (s3_1 .eq. 0.0d0)  s3_1 = s3_0
            tors1(1,i) = lambda*v1_1 + (1.0d0-lambda)*v1_0
            tors1(2,i) = lambda*s1_1 + (1.0d0-lambda)*s1_0
            tors2(1,i) = lambda*v2_1 + (1.0d0-lambda)*v2_0
            tors2(2,i) = lambda*s2_1 + (1.0d0-lambda)*s2_0
            tors3(1,i) = lambda*v3_1 + (1.0d0-lambda)*v3_0
            tors3(2,i) = lambda*s3_1 + (1.0d0-lambda)*s3_0
            if (verbose) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Hybrid Torsional Parameters :',
     &                    //,6x,'Atom Numbers',17x,'KT1',14x,'KT2',
     &                       14x,'KT3',/)
               end if
               write (iout,40)  ia,ib,ic,id,tors1(1,i),tors1(2,i),
     &                          tors2(1,i),tors2(2,i),
     &                          tors3(1,i),tors3(2,i)
   40          format (1x,4i5,4x,3(f10.4,f7.3))
            end if
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine hstrtor  --  hybrid stretch-torsion terms  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "hstrtor" constructs the hybrid stretch-torsion parameters
c     given an initial state, final state and "lambda" value
c
c
      subroutine hstrtor
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'ksttor.i'
      include 'mutant.i'
      include 'strtor.i'
      include 'tors.i'
      integer i,j,k,itb,itc
      integer ia,ib,ic,id,size
      character*3 pb,pc
      character*6 pt
      real*8 kst0(3),kst1(3)
      logical header
c
c
c     assign hybrid parameters for the stretch-torsion sites
c
      header = .true.
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         if (alter(ib) .or. alter(ic)) then
            itb = class(ib)
            itc = class(ic)
c
c     find the stretch-torsion parameters for the initial state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
            end do
            size = 3
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (itb .le. itc) then
               pt = pb//pc
            else
               pt = pc//pb
            end if
            do k = 1, 3
               kst0(k) = 0.0d0
            end do
            do j = 1, maxnbt
               if (kbt(j) .eq. pt) then
                  do k = 1, 3
                     kst0(k) = btcon(k,j)
                  end do
                  goto 10
               end if
            end do
   10       continue
c
c     find the stretch-torsion parameters for the final state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ib)  itb = class1(j)
               if (k .eq. ic)  itc = class1(j)
            end do
            size = 3
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (itb .le. itc) then
               pt = pb//pc
            else
               pt = pc//pb
            end if
            do k = 1, 3
               kst1(k) = 0.0d0
            end do
            do j = 1, maxnbt
               if (kbt(j) .eq. pt) then
                  do k = 1, 3
                     kst1(k) = btcon(k,j)
                  end do
                  goto 20
               end if
            end do
   20       continue
c
c     form hybrid parameters for the current stretch-torsion
c
            do j = 1, 3
               kst(j,i) = lambda*kst1(j) + (1.0d0-lambda)*kst0(j)
            end do
            if (kst(1,i).eq.0.0d0 .and. kst(2,i).eq.0.0d0
     &                  .and. kst(3,i).eq.0.0d0) then
               if (ist(1,i) .ne. 0) then
                  nstrtor = nstrtor - 1
                  ist(1,i) = 0
               end if
            else
               if (ist(1,i) .ne. i) then
                  nstrtor = nstrtor + 1
                  ist(1,i) = i
                  do j = 1, n12(ib)
                     if (i12(j,ib) .eq. ic) then
                        ist(2,i) = bndlist(j,ib)
                        goto 30
                     end if
                  end do
   30             continue
               end if
               if (verbose) then
                  if (header) then
                     header = .false.
                     write (iout,40)
   40                format (/,' Hybrid Stretch-Torsion Parameters :',
     &                       //,6x,'Atom Numbers',13x,'KST1',8x,'KST2',
     &                          8x,'KST3',/)
                  end if
                  write (iout,50)  ia,ib,ic,id,(kst(j,i),j=1,3)
   50             format (3x,4i5,3f12.3)
               end if
            end if
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine hvdw  --  hybrid van der Waals parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "hvdw" constructs the hybrid van der Waals interaction
c     parameters given an initial state, final state and
c     value of the mutation parameter "lambda"
c
c
      subroutine hvdw
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kvdws.i'
      include 'math.i'
      include 'mutant.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,it,it0,it1
      real*8 radius,rd,ep
      real*8 srad(maxclass)
      real*8 seps(maxclass)
      logical header
c
c
c     assign the hybrid van der Waals parameters
c
      do j = 1, nhybrid
         i = ihybrid(j)
         it = class(i)
         it0 = class0(j)
         it1 = class1(j)
         rad(it) = lambda*rad(it1) + (1.0d0-lambda)*rad(it0)
         eps(it) = lambda*eps(it1) + (1.0d0-lambda)*eps(it0)
      end do
c
c     get the square roots of the vdw radii and well depths
c
      do i = 1, maxclass
         srad(i) = sqrt(rad(i))
         seps(i) = sqrt(eps(i))
      end do
c
c     use combination rules to set pairwise vdw radii sums
c
      do j = 1, nhybrid
         i = ihybrid(j)
         it = class(i)
         do k = 1, maxclass
            if (rad(it).eq.0.0d0 .and. rad(k).eq.0.0d0) then
               rd = 0.0d0
            else if (radrule(1:10) .eq. 'ARITHMETIC') then
               rd = rad(it) + rad(k)
            else if (radrule(1:9) .eq. 'GEOMETRIC') then
               rd = 2.0d0*(srad(it)*srad(k))
            else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
               rd = 2.0d0*(rad(it)**3+rad(k)**3)/(rad(it)**2+rad(k)**2)
            else
               rd = rad(it) + rad(k)
            end if
            radmin(it,k) = rd
            radmin(k,it) = rd
         end do
      end do
c
c     use combination rules to set pairwise well depths
c
      do j = 1, nhybrid
         i = ihybrid(j)
         it = class(i)
         do k = 1, maxclass
            if (eps(it).eq.0.0d0 .and. eps(k).eq.0.0d0) then
               ep = 0.0d0
            else if (epsrule(1:10) .eq. 'ARITHMETIC') then
               ep = 0.5d0 * (eps(it) + eps(k))
            else if (epsrule(1:9) .eq. 'GEOMETRIC') then
               ep = seps(it) * seps(k)
            else if (epsrule(1:8) .eq. 'HARMONIC') then
               ep = 2.0d0 * (eps(it)*eps(k)) / (eps(it)+eps(k))
            else if (epsrule(1:3) .eq. 'HHG') then
               ep = 4.0d0 * (eps(it)*eps(k)) / (seps(it)+seps(k))**2
            else
               ep = seps(it) * seps(k)
            end if
            epsilon(it,k) = ep
            epsilon(k,it) = ep
         end do
      end do
c
c     print the van der Waals parameters for hybrid atoms
c
      header = .true.
      do j = 1, nhybrid
         if (verbose) then
            if (header) then
               header = .false.
               write (iout,10)
   10          format (/,' Hybrid van der Waals Parameters :',
     &                 //,7x,'Atom Number    Radius     Epsilon')
            end if
            radius = rad(it)
            if (radsiz .eq. 'DIAMETER')  radius = 2.0d0 * radius
            if (radtyp .eq. 'SIGMA')  radius = radius / twosix
            write (iout,20)  i,radius,eps(it)
   20       format (6x,i8,f14.4,f12.4)
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine hcharge  --  find hybrid charge parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "hcharge" constructs the hybrid charge interaction energy
c     parameters given an initial state, final state and "lambda"
c
c
      subroutine hcharge
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kchrge.i'
      include 'mutant.i'
      integer i,j,k,it,it0,it1
      real*8 chg0,chg1,hybchg
      logical header,used
c
c
c     assign the hybrid parameters for atomic charges
c
      header = .true.
      do j = 1, nhybrid
         used = .false.
         i = ihybrid(j)
         it = type(i)
         it0 = type0(j)
         it1 = type1(j)
         chg0 = chg(it0)
         chg1 = chg(it1)
         hybchg = lambda*chg1 + (1.0d0-lambda)*chg0
         do k = 1, nion
            if (iion(k) .eq. i) then
               used = .true.
               pchg(k) = hybchg
               goto 10
            end if
         end do
         if (chg0.ne.0.0d0 .or. chg1.ne.0.0d0) then
            used = .true.
            nion = nion + 1
            iion(nion) = i
            kion(nion) = i
            pchg(nion) = hybchg
         end if
   10    continue
         if (verbose .and. used) then
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Hybrid Atomic Partial Charge Parameters :',
     &                 //,7x,'Atom Number',7x,'Charge',/)
            end if
            write (iout,30)  i,hybchg
   30       format (6x,i8,5x,f12.3)
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine hdipole  --  find hybrid dipole parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "hdipole" constructs the hybrid dipole interaction energy
c     parameters given an initial state, final state and "lambda"
c
c
      subroutine hdipole
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bond.i'
      include 'dipole.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kdipol.i'
      include 'mutant.i'
      integer i,j,k,ia,ib
      integer ita,itb,size
      character*3 pa,pb
      character*6 blank,pt
      real*8 dpl0,dpl1,hybdpl
      real*8 pos0,pos1,hybpos
      logical header,used
      data blank  / '      ' /
c
c
c     assign the hybrid parameters for bond dipoles
c
      header = .true.
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         if (alter(ia) .or. alter(ib)) then
            ita = type(ia)
            itb = type(ib)
c
c     find the dipole parameters for the initial state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = type0(j)
               if (k .eq. ib)  itb = type0(j)
            end do
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            dpl0 = 0.0d0
            pos0 = 0.5d0
            do j = 1, maxnd
               if (kd(j) .eq. blank)  goto 10
               if (kd(j) .eq. pt) then
                  if (ita .le. itb) then
                     dpl0 = bdpl(j)
                     pos0 = sdpl(j)
                  else
                     dpl0 = -bdpl(j)
                     pos0 = 1.0d0 - sdpl(j)
                  end if
               end if
            end do
   10       continue
c
c     find the dipole parameters for the final state
c
            do j = 1, nhybrid
               k = ihybrid(j)
               if (k .eq. ia)  ita = type1(j)
               if (k .eq. ib)  itb = type1(j)
            end do
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            dpl1 = 0.0d0
            pos1 = 0.5d0
            do j = 1, maxnd
               if (kd(j) .eq. blank)  goto 20
               if (kd(j) .eq. pt) then
                  if (ita .le. itb) then
                     dpl1 = bdpl(j)
                     pos1 = sdpl(j)
                  else
                     dpl1 = -bdpl(j)
                     pos1 = 1.0d0 - sdpl(j)
                  end if
               end if
            end do
   20       continue
c
c     form the hybrid parameters for the current dipole
c
            hybdpl = lambda*dpl1 + (1.0d0-lambda)*dpl0
            hybpos = lambda*pos1 + (1.0d0-lambda)*pos0
            used = .false.
            do j = 1, ndipole
               if ((idpl(1,j).eq.ia .and. idpl(2,j).eq.ib) .or.
     &             (idpl(1,j).eq.ib .and. idpl(2,j).eq.ia)) then
                  idpl(1,j) = ia
                  idpl(2,j) = ib
                  bdpl(j) = hybdpl
                  sdpl(j) = hybpos
                  used = .true.
                  goto 30
               end if
            end do
            if (hybdpl .ne. 0.0d0) then
               ndipole = ndipole + 1
               idpl(1,ndipole) = ia
               idpl(2,ndipole) = ib
               bdpl(ndipole) = hybdpl
               sdpl(ndipole) = hybpos
               used = .true.
            end if
   30       continue
            if (verbose .and. used) then
               if (header) then
                  header = .false.
                  write (iout,40)
   40             format (/,' Hybrid Bond Dipole Moment Parameters :',
     &                    //,6x,'Atom Numbers',7x,'Moment',
     &                       7x,'Position',/)
               end if
               write (iout,50)  ia,ib,hybdpl,hybpos
   50          format (6x,2i5,2f15.3)
            end if
         end if
      end do
      return
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine image  --  pairwise distance of minimum image  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in the same or neighboring periodic boxes and
c     converts to the components of the minimum image distance
c
c
      subroutine image (xr,yr,zr,i)
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'cell.i'
      integer i
      real*8 xr,yr,zr
      real*8 xmove,ymove,zmove
      real*8 xfrac,yfrac,zfrac
c
c
c     compute the distance to translate along each cell axis
c
      if (i .eq. 0) then
         xmove = 0.0d0
         ymove = 0.0d0
         zmove = 0.0d0
      else
         xmove = icell(1,i) * xbox
         ymove = icell(2,i) * ybox
         zmove = icell(3,i) * zbox
      end if
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = xr + xmove
         dowhile (xr .gt. xcell2)
            xr = xr - xcell
         end do
         dowhile (xr .lt. -xcell2)
            xr = xr + xcell
         end do
         yr = yr + ymove
         dowhile (yr .gt. ycell2)
            yr = yr - ycell
         end do
         dowhile (yr .lt. -ycell2)
            yr = yr + ycell
         end do
         zr = zr + zmove
         dowhile (zr .gt. zcell2)
            zr = zr - zcell
         end do
         dowhile (zr .lt. -zcell2)
            zr = zr + zcell
         end do
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (monoclinic) then
         zfrac = zr / beta_sin
         xfrac = xr - zfrac*beta_cos
         xfrac = xfrac + xmove
         dowhile (xfrac .gt. xcell2)
            xfrac = xfrac - xcell
         end do
         dowhile (xfrac .lt. -xcell2)
            xfrac = xfrac + xcell
         end do
         yr = yr + ymove
         dowhile (yr .gt. ycell2)
            yr = yr - ycell
         end do
         dowhile (yr .lt. -ycell2)
            yr = yr + ycell
         end do
         zfrac = zfrac + zmove
         dowhile (zfrac .gt. zcell2)
            zfrac = zfrac - zcell
         end do
         dowhile (zfrac .lt. -zcell2)
            zfrac = zfrac + zcell
         end do
         xr = xfrac + zfrac*beta_cos
         zr = zfrac * beta_sin
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (triclinic) then
         zfrac = zr / gamma_term
         yfrac = (yr - zfrac*beta_term) / gamma_sin
         xfrac = xr - yfrac*gamma_cos - zfrac*beta_cos
         xfrac = xfrac + xmove
         dowhile (xfrac .gt. xcell2)
            xfrac = xfrac - xcell
         end do
         dowhile (xfrac .lt. -xcell2)
            xfrac = xfrac + xcell
         end do
         yfrac = yfrac + ymove
         dowhile (yfrac .gt. ycell2)
            yfrac = yfrac - ycell
         end do
         dowhile (yfrac .lt. -ycell2)
            yfrac = yfrac + ycell
         end do
         zfrac = zfrac + zmove
         dowhile (zfrac .gt. zcell2)
            zfrac = zfrac - zcell
         end do
         dowhile (zfrac .lt. -zcell2)
            zfrac = zfrac + zcell
         end do
         xr = xfrac + yfrac*gamma_cos + zfrac*beta_cos
         yr = yfrac*gamma_sin + zfrac*beta_term
         zr = zfrac * gamma_term
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
      else if (octahedron) then
         dowhile (xr .gt. xbox2)
            xr = xr - xbox
         end do
         dowhile (xr .lt. -xbox2)
            xr = xr + xbox
         end do
         dowhile (yr .gt. ybox2)
            yr = yr - ybox
         end do
         dowhile (yr .lt. -ybox2)
            yr = yr + ybox
         end do
         dowhile (zr .gt. zbox2)
            zr = zr - zbox
         end do
         dowhile (zr .lt. -zbox2)
            zr = zr + zbox
         end do
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end
c
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine impose  --  superimpose two coordinate sets  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "impose" performs the least squares best superposition
c     of two atomic coordinate sets via a quaternion method;
c     upon return, the first coordinate set is unchanged while
c     the second set is translated and rotated to give best fit;
c     the final root mean square fit is returned in "rmsvalue"
c
c
      subroutine impose (n1,x1,y1,z1,n2,x2,y2,z2,rmsvalue)
      implicit none
      include 'sizes.i'
      include 'align.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,n1,n2
      real*8 x1(maxatm),y1(maxatm),z1(maxatm)
      real*8 x2(maxatm),y2(maxatm),z2(maxatm)
      real*8 xmid,ymid,zmid
      real*8 rmsvalue,rmsfit
c
c
c     superimpose the full structures if not specified
c
      if (nfit .eq. 0) then
         nfit = min(n1,n2)
         do i = 1, nfit
            ifit(1,i) = i
            ifit(2,i) = i
            wfit(i) = 1.0d0
         end do
      end if
c
c     if the weights are all zero, set them to unity
c
      do i = 1, nfit
         if (wfit(i) .ne. 0.0d0)  goto 10
      end do
      do i = 1, nfit
         wfit(i) = 1.0d0
      end do
   10 continue
c
c     find the rms fit of input coordinates
c
      if (verbose) then
         rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2)
         write (iout,20)  rmsvalue
   20    format (/,' IMPOSE  --  Input Coordinates',12x,f12.6)
      end if
c
c     superimpose the centroids of active atom pairs
c
      call center (n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid)
      if (verbose) then
         rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2)
         write (iout,30)  rmsvalue
   30    format (' IMPOSE  --  After Translation',12x,f12.6)
      end if
c
c     use a quaternion method to achieve the superposition
c
      call quatfit (n1,x1,y1,z1,n2,x2,y2,z2)
      rmsvalue = rmsfit (x1,y1,z1,x2,y2,z2)
      if (verbose) then
         write (iout,40)  rmsvalue
   40    format (' IMPOSE  --  After Rotation',15x,f12.6)
      end if
c
c     translate both coordinate sets so as to return
c     the first set to its original position
c
      do i = 1, n1
         x1(i) = x1(i) + xmid
         y1(i) = y1(i) + ymid
         z1(i) = z1(i) + zmid
      end do
      do i = 1, n2
         x2(i) = x2(i) + xmid
         y2(i) = y2(i) + ymid
         z2(i) = z2(i) + zmid
      end do
      return
      end
c
c
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce --computes the induced dipole moment  ##
c     ##                                                          ##
c     ##############################################################
c
c     "induce" computes the induced dipole moment at each
c     polarizable site due to direct or mutual polarization
c
c
      subroutine induce
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,m,skip(maxatm)
      integer ii,iz,ix,kk,kz,kx
      integer iter,maxiter
      real*8 eps,norm
      real*8 xr,yr,zr,taper
      real*8 r,r2,r3,r4,r5
      real*8 rpi(13),rpk(13)
      real*8 udir(3,maxatm)
      real*8 uold(3,maxatm)
      real*8 fieldi(3),fieldk(3)
      logical iuse,kuse
c
c
c     zero out induced dipoles and count the polarizable sites
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the switching function coefficients
c
      call switch ('CHARGE')
c
c     compute the direct induced dipole moment at each atom
c
      do ii = 1, npole-1
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * chg12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * chg13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * chg14use
         end do
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
            if (iuse .or. kuse) then
               if (skip(k) .ne. i) then
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  if (use_image)  call image (xr,yr,zr,0)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     call udirect (ii,kk,xr,yr,zr,r,r2,rpi,
     &                                rpk,fieldi,fieldk)
                     if (skip(k) .eq. -i) then
                        do j = 1, 3
                           fieldi(j) = fieldi(j) / chgscale
                           fieldk(j) = fieldk(j) / chgscale
                        end do
                     end if
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        do j = 1, 3
                           fieldi(j) = fieldi(j) * taper
                           fieldk(j) = fieldk(j) * taper
                        end do
                     end if
                     do j = 1, 3
                        uind(j,ii) = uind(j,ii) + fieldi(j)
                        uind(j,kk) = uind(j,kk) + fieldk(j)
                     end do                 
                  end if
               end if
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do ii = 1, npole
            i = ipole(ii)
            iz = zaxis(ii)
            ix = xaxis(ii)
            iuse = (use(i) .or. use(iz) .or. use(ix))
            do j = 1, maxpole
               rpi(j) = rpole(j,ii)
            end do
            do kk = ii, npole
               k = ipole(kk)
               kz = zaxis(kk)
               kx = xaxis(kk)
               kuse = (use(k) .or. use(kz) .or. use(kx))
               if (iuse .or. kuse) then
                  do m = 1, ncell
                     xr = x(k) - x(i)
                     yr = y(k) - y(i)
                     zr = z(k) - z(i)
                     call image (xr,yr,zr,m)
                     r2 = xr*xr + yr*yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        do j = 1, maxpole
                           rpk(j) = rpole(j,kk)
                        end do
                        call udirect (ii,kk,xr,yr,zr,r,r2,rpi,
     &                                   rpk,fieldi,fieldk)
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r2 * r3
                           taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                           do j = 1, 3
                              fieldi(j) = fieldi(j) * taper
                              fieldk(j) = fieldk(j) * taper
                           end do
                        end if
                        do j = 1, 3
                           uind(j,ii) = uind(j,ii) + fieldi(j)
                           if (ii .ne. kk)
     &                        uind(j,kk) = uind(j,kk) + fieldk(j)
                        end do                 
                     end if
                  end do
               end if
            end do
         end do
      end if
c
c     convert total field at each atom to induced dipole moment
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = polarize(i) * uind(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         maxiter = 100
         iter = 0
         eps = 1.0d0
         do i = 1, npole
            do j = 1, 3
               udir(j,i) = uind(j,i)
            end do
         end do
c
c     compute mutual induced dipole moments by an iterative method
c
         dowhile (iter.lt.maxiter .and. eps.gt.poleps)
            do ii = 1, npole
               do j = 1, 3
                  uold(j,ii) = uind(j,ii)
                  uind(j,ii) = 0.0d0
               end do
            end do
            do ii = 1, npole-1
               i = ipole(ii)
               iz = zaxis(ii)
               ix = xaxis(ii)
               iuse = (use(i) .or. use(iz) .or. use(ix))
               do j = 1, n12(i)
                  skip(i12(j,i)) = i * chg12use
               end do
               do j = 1, n13(i)
                  skip(i13(j,i)) = i * chg13use
               end do
               do j = 1, n14(i)
                  skip(i14(j,i)) = i * chg14use
               end do
               do kk = ii+1, npole
                  k = ipole(kk)
                  kz = zaxis(kk)
                  kx = xaxis(kk)
                  kuse = (use(k) .or. use(kz) .or. use(kx))
                  if (iuse .or. kuse) then
                     if (skip(k) .ne. i) then
                        xr = x(k) - x(i)
                        yr = y(k) - y(i)
                        zr = z(k) - z(i)
                        if (use_image)  call image (xr,yr,zr,0)
                        r2 = xr*xr + yr*yr + zr*zr
                        if (r2 .le. off2) then
                           r = sqrt(r2)
                           call umutual (ii,kk,xr,yr,zr,r,r2,uold(1,ii),
     &                                      uold(1,kk),fieldi,fieldk)
                           if (skip(k) .eq. -i) then
                              do j = 1, 3
                                 fieldi(j) = fieldi(j) / chgscale
                                 fieldk(j) = fieldk(j) / chgscale
                              end do
                           end if
                           if (r2 .gt. cut2) then
                              r3 = r2 * r
                              r4 = r2 * r2
                              r5 = r2 * r3
                              taper = c5*r5 + c4*r4 + c3*r3
     &                                   + c2*r2 + c1*r + c0
                              do j = 1, 3
                                 fieldi(j) = fieldi(j) * taper
                                 fieldk(j) = fieldk(j) * taper
                              end do
                           end if
                           do j = 1, 3
                              uind(j,ii) = uind(j,ii) + fieldi(j)
                              uind(j,kk) = uind(j,kk) + fieldk(j)
                           end do                 
                        end if
                     end if
                  end if
               end do
            end do
c
c     periodic boundary for large cutoffs via replicates method
c
            if (use_replica) then
               do ii = 1, npole
                  i = ipole(ii)
                  iz = zaxis(ii)
                  ix = xaxis(ii)
                  iuse = (use(i) .or. use(iz) .or. use(ix))
                  do kk = ii, npole
                     k = ipole(kk)
                     kz = zaxis(kk)
                     kx = xaxis(kk)
                     kuse = (use(k) .or. use(kz) .or. use(kx))
                     if (iuse .or. kuse) then
                        do m = 1, ncell
                           xr = x(k) - x(i)
                           yr = y(k) - y(i)
                           zr = z(k) - z(i)
                           call image (xr,yr,zr,m)
                           r2 = xr*xr + yr*yr + zr*zr
                           if (r2 .le. off2) then
                              r = sqrt(r2)
                              call umutual (ii,kk,xr,yr,zr,r,r2,
     &                                      uold(1,ii),uold(1,kk),
     &                                         fieldi,fieldk)
                              if (r2 .gt. cut2) then
                                 r3 = r2 * r
                                 r4 = r2 * r2
                                 r5 = r2 * r3
                                 taper = c5*r5 + c4*r4 + c3*r3
     &                                      + c2*r2 + c1*r + c0
                                 do j = 1, 3
                                    fieldi(j) = fieldi(j) * taper
                                    fieldk(j) = fieldk(j) * taper
                                 end do
                              end if
                              do j = 1, 3
                                 uind(j,ii) = uind(j,ii) + fieldi(j)
                                 if (ii .ne. kk)
     &                              uind(j,kk) = uind(j,kk) + fieldk(j)
                              end do                 
                           end if
                        end do
                     end if
                  end do
               end do
            end if
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            eps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = udir(j,i) + polarize(i)*uind(j,i)
                  eps = eps + (uind(j,i)-uold(j,i))**2
               end do
            end do
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
         end do
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     print out a list of the final induced dipole moments
c
      if (debug) then
         write (iout,40)
   40    format (/,' Induced Dipole Moments (Debyes) :',
     &           //,4x,'Atom',14x,'X',11x,'Y',11x,'Z',9x,'Total'/)
         do i = 1, npole
            if (polarize(i) .ne. 0.0d0) then
               k = ipole(i)
               norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
               write (iout,50)  k,(debye*uind(j,i),j=1,3),debye*norm
   50          format (i8,5x,4f12.4)
            end if
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine udirect  --  field for direct induced dipoles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "udirect" evaluates the electric field at a polarizable atom
c     due to permanent atomic multipoles at a second atom, and vice
c     versa, for use in computation of direct induced dipole moments
c
c
      subroutine udirect (ii,kk,xr,yr,zr,r,r2,rpi,rpk,fieldi,fieldk)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      integer i,ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,damp
      real*8 rr3,rr5,rr7
      real*8 xr5,yr5,xyz
      real*8 xr7,yr7,rr53
      real*8 t2,t3,t4,t5,t6,t7,t8,t9,t10,t11
      real*8 t12,t13,t14,t15,t16,t17,t18,t19
      real*8 i0,i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 k0,k1,k2,k3,k4,k5,k6,k7,k8,k9
      real*8 i9i4,i9i7,k4k9,k7k9
      real*8 rpi(13),rpk(13)
      real*8 fieldi(3),fieldk(3)
      logical iquad,kquad
c
c
c     set extent of the multipole expansion at each site
c
      iquad = (polsiz(ii) .ge. 13)
      kquad = (polsiz(kk) .ge. 13)
c
c     compute the first order T2 matrix elements
c
      rr3 = -1.0d0 / (r2*r)
      t2 = xr * rr3
      t3 = yr * rr3
      t4 = zr * rr3
c
c     compute the second order T2 matrix elements
c
      rr5 = -3.0d0 * rr3 / r2
      xr5 = xr * rr5
      yr5 = yr * rr5
      t5 = rr3 + xr5*xr
      t6 = xr5 * yr
      t7 = xr5 * zr
      t8 = rr3 + yr5*yr
      t9 = yr5 * zr
      t10 = -t5 - t8
c
c     compute the third order T2 matrix elements
c
      if (iquad .or. kquad) then
         rr7 = -5.0d0 * rr5 / r2
         xyz = xr * yr * zr
         xr7 = xr * xr * rr7
         yr7 = yr * yr * rr7
         rr53 = 3.0d0 * rr5
         t11 = xr * (xr7+rr53)
         t12 = yr * (xr7+rr5)
         t13 = zr * (xr7+rr5)
         t14 = xr * (yr7+rr5)
         t15 = xyz * rr7
         t16 = -t11 - t14
         t17 = yr * (yr7+rr53)
         t18 = zr * (yr7+rr5)
         t19 = -t12 - t17
      end if
c
c     compute the field at site k due to multipoles at site i
c
      i0 = rpi(1)
      i1 = rpi(2)
      i2 = rpi(3)
      i3 = rpi(4)
      fieldk(1) = -i0*t2 + i1*t5 + i2*t6 + i3*t7
      fieldk(2) = -i0*t3 + i1*t6 + i2*t8 + i3*t9
      fieldk(3) = -i0*t4 + i1*t7 + i2*t9 + i3*t10
      if (iquad) then
         i4 = rpi(5)
         i5 = rpi(6)
         i6 = rpi(7)
         i7 = rpi(9)
         i8 = rpi(10)
         i9 = rpi(13)
         i9i4 = i9 - i4
         i9i7 = i9 - i7
         fieldk(1) = fieldk(1) + i9i4*t11 + i9i7*t14
     &                  - 2.0d0*(i5*t12 + i6*t13 + i8*t15)
         fieldk(2) = fieldk(2) + i9i4*t12 + i9i7*t17
     &                  - 2.0d0*(i5*t14 + i6*t15 + i8*t18)
         fieldk(3) = fieldk(3) + i9i4*t13 + i9i7*t18
     &                  - 2.0d0*(i5*t15 + i6*t16 + i8*t19)
      end if
c
c     compute the field at site i due to multipoles at site k
c
      k0 = rpk(1)
      k1 = rpk(2)
      k2 = rpk(3)
      k3 = rpk(4)
      fieldi(1) = k0*t2 + k1*t5 + k2*t6 + k3*t7
      fieldi(2) = k0*t3 + k1*t6 + k2*t8 + k3*t9
      fieldi(3) = k0*t4 + k1*t7 + k2*t9 + k3*t10
      if (kquad) then
         k4 = rpk(5)
         k5 = rpk(6)
         k6 = rpk(7)
         k7 = rpk(9)
         k8 = rpk(10)
         k9 = rpk(13)
         k4k9 = k4 - k9
         k7k9 = k7 - k9
         fieldi(1) = fieldi(1) + k4k9*t11 + k7k9*t14
     &                  + 2.0d0*(k5*t12 + k6*t13 + k8*t15)
         fieldi(2) = fieldi(2) + k4k9*t12 + k7k9*t17
     &                  + 2.0d0*(k5*t14 + k6*t15 + k8*t18)
         fieldi(3) = fieldi(3) + k4k9*t13 + k7k9*t18
     &                  + 2.0d0*(k5*t15 + k6*t16 + k8*t19)
      end if
c
c     apply a damping factor to reduce the field at short range
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            do i = 1, 3
               fieldi(i) = fieldi(i) * damp
               fieldk(i) = fieldk(i) * damp
            end do
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine umutual  --  field for mutual induced dipoles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "umutual" evaluates the electric field at a polarizable atom
c     due to the induced atomic dipoles at a second atom, and vice
c     versa, for use in computation of mutual induced dipole moments
c
c
      subroutine umutual (ii,kk,xr,yr,zr,r,r2,upi,upk,fieldi,fieldk)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      integer i,ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,damp
      real*8 rr3,rr5,xr5,yr5
      real*8 upi(3),upk(3)
      real*8 fieldi(3),fieldk(3)
      real*8 u1,u2,u3,v1,v2,v3
      real*8 t5,t6,t7,t8,t9,t10
c
c
c     set the current values of the induced dipoles
c
      u1 = upi(1)
      u2 = upi(2)
      u3 = upi(3)
      v1 = upk(1)
      v2 = upk(2)
      v3 = upk(3)
c
c     compute the second order T2 matrix elements
c
      rr3 = -1.0d0 / (r2*r)
      rr5 = -3.0d0 * rr3 / r2
      xr5 = xr * rr5
      yr5 = yr * rr5
      t5 = rr3 + xr5*xr
      t6 = xr5 * yr
      t7 = xr5 * zr
      t8 = rr3 + yr5*yr
      t9 = yr5 * zr
      t10 = -t5 - t8
c
c     compute the field at site i due to site k and vice versa
c
      fieldi(1) = v1*t5 + v2*t6 + v3*t7
      fieldi(2) = v1*t6 + v2*t8 + v3*t9
      fieldi(3) = v1*t7 + v2*t9 + v3*t10
      fieldk(1) = u1*t5 + u2*t6 + u3*t7
      fieldk(2) = u1*t6 + u2*t8 + u3*t9
      fieldk(3) = u1*t7 + u2*t9 + u3*t10
c
c     apply a damping factor to reduce the field at short range
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            do i = 1, 3
               fieldi(i) = fieldi(i) * damp
               fieldk(i) = fieldk(i) * damp
            end do
         end if
      end if
      return
      end
c
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine inertia  --  principal moments of inertia  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "inertia" computes the principal moments of inertia for the
c     system, and optionally translates the center of mass to the
c     origin and rotates the principal axes onto the global axes
c
c        mode = 1     print the moments and principal axes
c        mode = 2     move coordinates to standard orientation
c        mode = 3     perform both of the above operations
c
c     literature reference:
c
c     Herbert Goldstein, "Classical Mechanics, 2nd Edition",
c     Addison-Wesley, Reading, MA, 1980; see the Euler angle
c     xyz convention in Appendix B
c
c
      subroutine inertia (mode)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'iounit.i'
      include 'math.i'
      integer i,j,k,mode
      real*8 weigh,total,dot
      real*8 xcm,ycm,zcm
      real*8 xx,xy,xz,yy,yz,zz
      real*8 xterm,yterm,zterm
      real*8 moment(3),vec(3,3)
      real*8 work1(3),work2(3)
      real*8 tensor(3,3),a(3,3)
      logical print,move
c
c
c     decide upon the type of output desired
c
      print = .false.
      move = .false.
      if (mode.eq.1 .or. mode.eq.3)  print = .true.
      if (mode.eq.2 .or. mode.eq.3)  move = .true.
c
c     compute the position of the center of mass
c
      total = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      do i = 1, n
         weigh = mass(i)
         total = total + weigh
         xcm = xcm + x(i)*weigh
         ycm = ycm + y(i)*weigh
         zcm = zcm + z(i)*weigh
      end do
      xcm = xcm / total
      ycm = ycm / total
      zcm = zcm / total
c
c     compute and then diagonalize the inertial tensor
c
      xx = 0.0d0
      xy = 0.0d0
      xz = 0.0d0
      yy = 0.0d0
      yz = 0.0d0
      zz = 0.0d0
      do i = 1, n
         weigh = mass(i)
         xterm = x(i) - xcm
         yterm = y(i) - ycm
         zterm = z(i) - zcm
         xx = xx + xterm*xterm*weigh
         xy = xy + xterm*yterm*weigh
         xz = xz + xterm*zterm*weigh
         yy = yy + yterm*yterm*weigh
         yz = yz + yterm*zterm*weigh
         zz = zz + zterm*zterm*weigh
      end do
      tensor(1,1) = yy + zz
      tensor(2,1) = -xy
      tensor(3,1) = -xz
      tensor(1,2) = -xy
      tensor(2,2) = xx + zz
      tensor(3,2) = -yz
      tensor(1,3) = -xz
      tensor(2,3) = -yz
      tensor(3,3) = xx + yy
      call tnk_jacobi (3,3,tensor,moment,vec,work1,work2)
c
c     select the direction for each principal moment axis
c
      do i = 1, 2
         do j = 1, n
            xterm = vec(1,i) * (x(j)-xcm)
            yterm = vec(2,i) * (y(j)-ycm)
            zterm = vec(3,i) * (z(j)-zcm)
            dot = xterm + yterm + zterm
            if (dot .lt. 0.0d0) then
               do k = 1, 3
                  vec(k,i) = -vec(k,i)
               end do
            end if
            if (dot .ne. 0.0d0)  goto 10
         end do
   10    continue
      end do
c
c     moment axes must give a right-handed coordinate system
c
      xterm = vec(1,1) * (vec(2,2)*vec(3,3)-vec(2,3)*vec(3,2))
      yterm = vec(2,1) * (vec(1,3)*vec(3,2)-vec(1,2)*vec(3,3))
      zterm = vec(3,1) * (vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2))
      dot = xterm + yterm + zterm
      if (dot .lt. 0.0d0) then
         do j = 1, 3
            vec(j,3) = -vec(j,3)
         end do
      end if
c
c     print the moments of inertia and the principal axes
c
      if (print) then
         write (iout,20)
   20    format (/,' Moments of Inertia and Principal Axes :',
     &           //,9x,'Moment',11x,'X-, Y- and Z-Components of Axis')
         write (iout,30)  (moment(i),vec(1,i),vec(2,i),vec(3,i),i=1,3)
   30    format (3(/,f16.4,5x,3f12.6))
      end if
c
c     principal moment axes form rows of Euler rotation matrix
c
      if (move) then
         do i = 1, 3
            do j = 1, 3
               a(i,j) = vec(j,i)
            end do
         end do
c
c     translate to origin, then apply Euler rotation matrix
c
         do i = 1, n
            xterm = x(i) - xcm
            yterm = y(i) - ycm
            zterm = z(i) - zcm
            x(i) = a(1,1)*xterm + a(1,2)*yterm + a(1,3)*zterm
            y(i) = a(2,1)*xterm + a(2,2)*yterm + a(2,3)*zterm
            z(i) = a(3,1)*xterm + a(3,2)*yterm + a(3,3)*zterm
         end do
      end if
      return
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine initial  --  initial values and program setup  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "initial" sets up original values for some parameters
c     and variables that might not otherwise get initialized
c
c
      subroutine initial
      implicit none
      include 'sizes.i'
      include 'align.i'
      include 'argue.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cell.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'linmin.i'
      include 'minima.i'
      include 'molcul.i'
      include 'mutant.i'
      include 'output.i'
      include 'pdb.i'
      include 'precis.i'
      include 'scales.i'
      include 'sequen.i'
      include 'warp.i'
      include 'zclose.i'
      real*8 precise
c
c
c     number of atoms used in superposition
c
      nfit = 0
c
c     number of command line arguments
c
      narg = 0
c
c     number of atoms in the system
c
      n = 0
c
c     flags for temperature and pressure baths
c
      isothermal = .false.
      isobaric = .false.
c
c     flags for periodic boundaries
c
      use_bounds = .false.
      use_image = .false.
      use_replica = .false.
c
c     number of unit cell replicates
c
      ncell = 0
c
c     highest numbered previous cycle file
c
      nprior = 0
c
c     information levels within the program
c
      verbose = .false.
      debug = .false.
      abort = .false.
c
c     default input/output unit numbers
c
      input = 5
      iout = 6
c
c     number of lines in the keyfile
c
      nkey = 0
c
c     type of coordinates file
c
      coordtype = 'none'
c
c     number of molecules in the system
c
      nmol = 0
c
c     number of hybrid (mutant) atoms in the system
c
      nhybrid = 0
c
c     number of atoms in Protein Data Bank ATOM records
c
      npdb = 0
c
c     generic parameters used by all optimizations
c
      fctmin = 0.0d0
      maxiter = 0
      nextiter = 0
      iprint = -1
      iwrite = -1
c
c     generic parameters used during line search
c
      cappa = 0.0d0
      stpmin = 0.0d0
      stpmax = 0.0d0
      angmax = 0.0d0
      intmax = 0
c
c     flag to show setting of optimization scale factors
c
      set_scale = .false.
c
c     number of residues and chains in biopolymer sequence
c
      nseq = 0
      nchain = 0
c
c     flag for Gaussian density annealing algorithm
c
      use_gda = .false.
c
c     number of bonds added or deleted from Z-matrix
c
      nadd = 0
      ndel = 0
c
c     display program info and copyright notice
c

      call promo

c
c     determine machine tolerence values
c
      tiny = precise (1)
      small = precise (2)
      huge = precise (3)
c
c     get any command line arguments to the program
c
cjrs  no command line arguments used when run under GAMESS
c
c      call command
cjrs
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine initprm  --  initialize force field parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "initprm" completely initializes a force field by setting
c     all values to zero prior to reading in a parameter set
c
c
      subroutine initprm
      implicit none
      include 'sizes.i'
      include 'fields.i'
      include 'kanang.i'
      include 'kangs.i'
      include 'katoms.i'
      include 'kbonds.i'
      include 'kchrge.i'
      include 'kdipol.i'
      include 'khbond.i'
      include 'kiprop.i'
      include 'kitors.i'
      include 'kmulti.i'
      include 'kopbnd.i'
      include 'korbs.i'
      include 'kpolr.i'
      include 'kstbnd.i'
      include 'ksttor.i'
      include 'ktorsn.i'
      include 'kurybr.i'
      include 'kvdws.i'
      include 'kvdwpr.i'
      integer i,j
c
c
c     initialize strings of parameter atom types and classes
c
      do i = 1, maxnvp
         kvpr(i) = '      '
      end do
      do i = 1, maxnhb
         khb(i) = '      '
      end do
      do i = 1, maxnb
         kb(i) = '      '
      end do
      do i = 1, maxnb5
         kb5(i) = '      '
      end do
      do i = 1, maxnb4
         kb4(i) = '      '
      end do
      do i = 1, maxnb3
         kb3(i) = '      '
      end do
      do i = 1, maxna
         ka(i) = '         '
      end do
      do i = 1, maxna5
         ka5(i) = '         '
      end do
      do i = 1, maxna4
         ka4(i) = '         '
      end do
      do i = 1, maxna3
         ka3(i) = '         '
      end do
      do i = 1, maxnopb
         kaopb(i) = '      '
      end do
      do i = 1, maxnu
         ku(i) = '         '
      end do
      do i = 1, maxndi
         kdi(i) = '            '
      end do
      do i = 1, maxnti
         kti(i) = '            '
      end do
      do i = 1, maxnt
         kt(i) = '            '
      end do
      do i = 1, maxnt5
         kt5(i) = '            '
      end do
      do i = 1, maxnt4
         kt4(i) = '            '
      end do
      do i = 1, maxnbt
         kbt(i) = '      '
      end do
      do i = 1, maxnd
         kd(i) = '      '
      end do
      do i = 1, maxnd5
         kd5(i) = '      '
      end do
      do i = 1, maxnd4
         kd4(i) = '      '
      end do
      do i = 1, maxnd3
         kd3(i) = '      '
      end do
      do i = 1, maxnmp
         kmp(i) = '         '
      end do
      do i = 1, maxnpi
         kpi(i) = '      '
      end do
c
c     initialize some of the force field parameters
c
      forcefield = '                    '
      do i = 1, maxtyp
         symbol(i) = '   '
         atmcls(i) = 0
         atmnum(i) = 0
         weight(i) = 0.0d0
         ligand(i) = 0
         describe(i) = '                    '
         chg(i) = 0.0d0
         polr(i) = 0.0d0
      end do
      do i = 1, maxclass
         rad(i) = 0.0d0
         eps(i) = 0.0d0
         reduct(i) = 0.0d0
         do j = 1, 3
            stbn(j,i) = 0.0d0
            anan(j,i) = 0.0d0
         end do
         electron(i) = 0.0d0
         ionize(i) = 0.0d0
         repulse(i) = 0.0d0
      end do
      do i = 1, maxbio
         biotyp(i) = 0
      end do
      return
      end
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine initrot  --  set bonds for dihedral rotation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "initrot" asks for torsional angles which are to be rotated
c     in subsequent computation, it will automatically locate all
c     rotatable single bonds if desired; assumes that an appropriate
c     internal coordinates file has already been read in
c
c
      subroutine initrot
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'math.i'
      include 'omega.i'
      include 'ring.i'
      include 'zcoord.i'
      integer i,j,k,j1,j2,mode
      integer bond1,bond2
      integer attach1,attach2
      integer nfixed,ifixed(2,maxbnd)
      integer list(maxatm)
      character*80 record,string
      logical exist,query
      logical rotate,rotcheck
c
c
c     initialize the number of rotatable torsional angles
c
      nomega = 0
c
c     choose automatic or manual selection of torsional angles
c
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Selection of Torsional Angles for Rotation :',
     &           //,'    0  - Automatic Selection of Torsional Angles',
     &            /,'    1  - Manual Selection of Angles to Rotate',
     &            /,'    2  - Manual Selection of Angles to Freeze',
     &           //,' Enter the Method of Choice [0] :  ',$)
         read (input,30)  mode
   30    format (i10)
      end if
      if (mode.ne.1 .and. mode.ne.2)  mode = 0
c
c     manual selection of the torsional angles to be rotated
c
      if (mode .eq. 1) then
         dowhile (nomega .lt. maxrot)
            nomega = nomega + 1
            j1 = 0
            j2 = 0
            write (iout,40)  nomega
   40       format (/,' Enter Atoms in Rotatable Bond',i5,' :  ',$)
            read (input,50)  record
   50       format (a80)
            read (record,*,err=80,end=80)  j1,j2
            if (j1.eq.0 .and. j2.eq.0)  goto 80
            do i = 4, n
               if (iz(4,i) .eq. 0) then
                  bond1 = iz(1,i)
                  bond2 = iz(2,i)
                  attach1 = n12(bond1)
                  attach2 = n12(bond2)
                  if (attach1.gt.1 .and. attach2.gt.1) then
                     if ((bond1.eq.j1 .and. bond2.eq.j2) .or.
     &                   (bond1.eq.j2 .and. bond2.eq.j1)) then
                        if (rotcheck(bond1,bond2)) then
                           iomega(1,nomega) = bond1
                           iomega(2,nomega) = bond2
                           dihed(nomega) = ztors(i) / radian
                           zline(nomega) = i
                           goto 70
                        end if
                     end if
                  end if
               end if
            end do
            nomega = nomega - 1
            write (iout,60)  j1,j2
   60       format (/,' INITROT  --  Bond between Atoms',2i6,
     &                 ' is not Rotatable')
   70       continue
         end do
   80    continue
         nomega = nomega - 1
      end if
c
c     manual selection of the torsional angles to be frozen
c
      nfixed = 0
      if (mode .eq. 2) then
         do i = 1, maxrot
            ifixed(1,i) = 0
            ifixed(2,i) = 0
            write (iout,90)  i
   90       format (/,' Enter Atoms in Frozen Bond',i5,' :  ',$)
            read (input,100)  record
  100       format (a80)
            read (record,*,err=110,end=110)  ifixed(1,i),ifixed(2,i)
            if (ifixed(1,i).eq.0 .and. ifixed(2,i).eq.0)  goto 110
            nfixed = nfixed + 1
         end do
  110    continue
      end if
c
c     perform the automatic selection of torsional angles to rotate
c
      if (mode.eq.0 .or. mode.eq.2) then
         do i = 4, n
            if (iz(4,i) .eq. 0) then
               rotate = .true.
               bond1 = iz(1,i)
               bond2 = iz(2,i)
c
c     do not rotate a bond if either bonded atom is univalent
c
               attach1 = n12(bond1)
               attach2 = n12(bond2)
               if (attach1.le.1 .or. attach2.le.1)  rotate = .false.
c
c     do not rotate if both atoms of bond are in a small ring
c
               do j = 1, n
                  list(j) = 0
               end do
               do j = 1, nring3
                  do k = 1, 3
                     list(iring3(k,j)) = j
                  end do
                  if (list(bond1).eq.j .and. list(bond2).eq.j) then
                     rotate = .false.
                     goto 120
                  end if
               end do
               do j = 1, n
                  list(j) = 0
               end do
               do j = 1, nring4
                  do k = 1, 4
                     list(iring4(k,j)) = j
                  end do
                  if (list(bond1).eq.j .and. list(bond2).eq.j) then
                     rotate = .false.
                     goto 120
                  end if
               end do
               do j = 1, n
                  list(j) = 0
               end do
               do j = 1, nring5
                  do k = 1, 5
                     list(iring5(k,j)) = j
                  end do
                  if (list(bond1).eq.j .and. list(bond2).eq.j) then
                     rotate = .false.
                     goto 120
                  end if
               end do
               do j = 1, n
                  list(j) = 0
               end do
               do j = 1, nring6
                  do k = 1, 6
                     list(iring6(k,j)) = j
                  end do
                  if (list(bond1).eq.j .and. list(bond2).eq.j) then
                     rotate = .false.
                     goto 120
                  end if
               end do
  120          continue
c
c     do not rotate bonds explicitly frozen by the user
c
               if (mode.eq.2 .and. rotate) then
                  do j = 1, nfixed
                     j1 = ifixed(1,j)
                     j2 = ifixed(2,j)
                     if ((bond1.eq.j1 .and. bond2.eq.j2) .or.
     &                   (bond1.eq.j2 .and. bond2.eq.j1)) then
                        rotate = .false.
                        goto 130
                     end if
                  end do
               end if
  130          continue
c
c     do not rotate bonds with inactive atoms on both sides
c
               if (rotate) then
                  if (.not. rotcheck(bond1,bond2))  rotate = .false.
               end if
c
c     check for possible duplication of rotatable bonds
c
               if (rotate) then
                  do j = 1, nomega
                     j1 = iomega(1,j)
                     j2 = iomega(2,j)
                     if ((bond1.eq.j1 .and. bond2.eq.j2) .or.
     &                   (bond1.eq.j2 .and. bond2.eq.j1)) then
                        write (iout,140)  bond1,bond2
  140                   format (/,' INITROT  --  Rotation about',2i6,
     &                             ' occurs more than once in Z-matrix')
                        call fatal
                     end if
                  end do
                  nomega = nomega + 1
                  iomega(1,nomega) = bond1
                  iomega(2,nomega) = bond2
                  dihed(nomega) = ztors(i) / radian
                  zline(nomega) = i
               end if
            end if
         end do
      end if
c
c     write out the number of rotatable torsions to be used
c
      if (nomega .eq. 0) then
         write (iout,150)
  150    format (/,' INITROT  --  No Torsions for Subsequent',
     &              ' Computation')
         call fatal
      else if (nomega .gt. maxrot) then
         write (iout,160)
  160    format (/,' INITROT  --  Too many Torsions;',
     &              ' Increase MAXROT')
         call fatal
      end if
      write (iout,170)  nomega
  170 format (/,' Number of Torsions Used in Derivative',
     &           ' Computation :',i6)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function rotcheck  --  check for fixed atoms across bond  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotcheck" tests a specified candidate rotatable bond for
c     the disallowed case where inactive atoms are found on both
c     sides of the candidate bond
c
c
      function rotcheck (base,partner)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'rotate.i'
      include 'usage.i'
      integer i,base,partner
      logical rotcheck,value
      logical list(maxatm)
c
c
c     initialize status and find atoms on short side of the bond
c
      value = .true.
      call rotlist (base,partner)
c
c     rotation is allowed if all atoms on one side are active
c
      do i = 1, nrot
         if (.not. use(rot(i))) then
            value = .false.
            goto 10
         end if
      end do
   10 continue
c
c     if short side had inactive atoms, check the other side
c
      if (.not. value) then
         do i = 1, n
            list(i) = .true.
         end do
         do i = 1, nrot
            list(rot(i)) = .false.
         end do
         do i = 1, n
            if (list(i) .and. .not.use(i))  goto 20
         end do
         value = .true.
   20    continue
      end if
c
c     set the final return value of the function
c
      rotcheck = value
      return
      end
c
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine insert  --  insert atom into coordinates list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "insert" adds the specified atom to the Cartesian
c     coordinates list and shifts the remaining atoms
c
c
      subroutine insert (iatom)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,iatom
c
c
c     increase by one the total number of atoms
c
      n = n + 1
c
c     shift the atom coordinates, types and connectivities
c
      do i = n, iatom+1, -1
         name(i) = name(i-1)
         x(i) = x(i-1)
         y(i) = y(i-1)
         z(i) = z(i-1)
         type(i) = type(i-1)
         n12(i) = n12(i-1)
         do j = 1, n12(i)
            i12(j,i) = i12(j,i-1)
         end do
      end do
c
c     put new atom at the origin with a big atom type number
c
      name(iatom) = 'NEW'
      x(iatom) = 0.0d0
      y(iatom) = 0.0d0
      z(iatom) = 0.0d0
      type(iatom) = maxtyp
      n12(iatom) = 0
c
c     shift the connected atom lists to allow the insertion
c
      do i = 1, n
         do j = 1, n12(i)
            if (i12(j,i) .ge. iatom) then
               i12(j,i) = i12(j,i) + 1
            end if
         end do
      end do
c
c     write a message to describe the atom insertion
c
      if (debug) then
         write (iout,10)  iatom
  10     format (' INSERT  --  Inserting Atom Number :',i8)
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine invert  --  gauss-jordan matrix inversion  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "invert" inverts a matrix using the Gauss-Jordan method
c
c       n    logical dimension of the matrix to be inverted
c       np   physical dimension of the matrix storage area
c       a    matrix to invert; contains inverse on exit
c
c
      subroutine invert (n,np,a)
      implicit none
      include 'iounit.i'
      integer maxinv
      parameter (maxinv=100)
      integer i,j,k,n,np,icol,irow
      integer ipivot(maxinv)
      integer indxc(maxinv)
      integer indxr(maxinv)
      real*8 big,temp,pivot
      real*8 a(np,np)
c
c
c     check to see if the matrix is too large to handle
c
      if (n .gt. maxinv) then
         write (iout,10)
   10    format (/,' INVERT  --  Matrix Too Large; Increase MAXINV')
         call fatal
      end if
c
c     perform matrix inversion via the Gauss-Jordan algorithm
c
      do i = 1, n
         ipivot(i) = 0
      end do
      do i = 1, n
         big = 0.0d0
         do j = 1, n
            if (ipivot(j) .ne. 1) then
               do k = 1, n
                  if (ipivot(k) .eq. 0) then
                     if (abs(a(j,k)) .ge. big) then
                        big = abs(a(j,k))
                        irow = j
                        icol = k
                     end if
                  else if (ipivot(k) .gt. 1) then
                     write (iout,20)
   20                format (/,' INVERT  --  Cannot Invert',
     &                          ' a Singular Matrix')
                     call fatal
                  end if
               end do
            end if
         end do
         ipivot(icol) = ipivot(icol) + 1
         if (irow .ne. icol) then
            do j = 1, n
               temp = a(irow,j)
               a(irow,j) = a(icol,j)
               a(icol,j) = temp
            end do
         end if
         indxr(i) = irow
         indxc(i) = icol
         if (a(icol,icol) .eq. 0.0d0) then
            write (iout,30)
   30       format (/,' INVERT  --  Cannot Invert a Singular Matrix')
            call fatal
         end if
         pivot = a(icol,icol)
         a(icol,icol) = 1.0d0
         do j = 1, n
            a(icol,j) = a(icol,j) / pivot
         end do
         do j = 1, n
            if (j .ne. icol) then
               temp = a(j,icol)
               a(j,icol) = 0.0d0
               do k = 1, n
                  a(j,k) = a(j,k) - a(icol,k)*temp
               end do
            end if
         end do
      end do
      do i = n, 1, -1
         if (indxr(i) .ne. indxc(i)) then
            do k = 1, n
               temp = a(k,indxr(i))
               a(k,indxr(i)) = a(k,indxc(i))
               a(k,indxc(i)) = temp
            end do
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine gdastat  --  compute GDA trajectory averages  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine gdastat (nstep,beta,xx,status)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'math.i'
      include 'warp.i'
      integer maxgda
      parameter (maxgda=4*maxatm)
      integer i,nstep,nvar
      real*8 beta,xx(maxgda)
      real*8 e,energy,rg,m2ave
      character*7 status
c
c
c     translate optimization parameters to coordinates and M2's
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
      do i = 1, n
         nvar = nvar + 1
         m2(i) = abs(xx(nvar))
      end do
c
c     get some info about the current integration step
c
      e = energy ()
      call gyrate (rg)
      m2ave = 0.0d0
      do i = 1, n
         m2ave = m2ave + m2(i)
      end do
      m2ave = m2ave / dble(n)
      write (iout,10)  nstep,log(beta)/logten,e,rg,
     &                 log(m2ave)/logten,status
   10 format (i6,2x,4f13.4,6x,a7)
c
c     save the current coordinates to a disk file
c
      call writeout (xx,nstep)
      return
      end

