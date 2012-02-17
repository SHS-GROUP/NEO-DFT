      PROGRAM GCISELCI 
c     --------------------
c
c     PROGRAM GCISELCI.F
c     ******************
c
c     SELECTS GENERAL CI-CONFIGURATIONS containing
c     up to sextuple excitations from the reference (closed-shell)
c     from the data in CI-SDT wavefunction
c
c     RESULT:  LISTOUT FILE listGCI  to be named file.gci
c     GAMESS SYNTAX: gms -file37 /full/name/of/file.gci file
c     -l file.log -q queue 
c
c     REQUIRES:
c     FILE.F37 = GCILIST
c     FILE.F12 = CIVECTR 
c
c-----------------------------------------------------------------
c     ****************** 
c     PROGRAM gciselCI.f
c     ****************** 
c     May 7, 2007 
c     May 25, 2007 
c     January 7, 2010
c     January 15, 2010
c
c     WRITTEN BY LAIMIS BYTAUTAS AND KLAUS RUEDENBERG
c
c     [Procedure is described in references:
c-----------------------------------------------------------------
c     1. L. BYTAUTAS, K. RUEDENBERG, CHEM. PHYS. 356, p.64 (2009).
c     2. L. BYTAUTAS, K. RUEDENBERG, 
c        ACS SYMPOSIUM SER. vol. 958, p.103 (2007). 
c
c     additional codes are from:
c     JOE IVANIC (SEE GCI files in GAMESS) 
c-----------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c     ----------------------------------
      integer posdet

      character*100 listout
      character*100 listoutCI_NNf1
      character*100 listoutCI_HNO
      character*100 listoutCI_NCCN
      character*100 list_NCCN
      character*100 listCI
      character*100 listGCI

      parameter (mxloop=50000)
      parameter (mxorb = 100)
      parameter (mxneb = 1000)
      parameter (mxlg = 50)

      parameter (mxprod = 1000000)
      parameter (mxsprod = 2400000)
      parameter (mxsprod1 = 100000)

      parameter (mxCI = 1000000)
      parameter (mxstring = 2000000)
      parameter (mxna = 200)
      parameter (mxCIsd = 200000)
      parameter (mxstr = 100)

c---------------------------------------------------
      dimension CIsd(mxCIsd)
      dimension ICONAsd(mxCIsd),ICONBsd(mxCIsd)
      dimension IEXCITsd(mxCIsd)

      dimension IACONsd(mxorb),IBCONsd(mxorb)
      dimension ICONsd(mxorb)
c     
      dimension IACON(mxorb),IBCON(mxorb)
      dimension ICON(mxorb)
      dimension ICON2(mxorb)
      dimension ICON3(mxorb)
      dimension ICON4(mxorb)
      dimension ICON5(mxorb)
      dimension ICON6(mxorb)
c     
      DIMENSION ISCF(mxorb),IVIR(mxorb),MATE(mxorb)
      DIMENSION MOSCF(mxorb),MOVIR(mxorb)
      DIMENSION IOCC(mxorb)
      DIMENSION IOCCV(mxorb)
c     
      DIMENSION CISUM3sd(mxsprod1),IDET3sd(mxsprod1)
      DIMENSION CISUMDsd(mxsprod1),IDETDsd(mxsprod1)
      DIMENSION CIDsd(mxsprod1)
      DIMENSION CI3sd(mxsprod1)
c     
      DIMENSION CISUM1sd(mxsprod1),IDET1sd(mxsprod1)
      DIMENSION DETDsd(mxsprod1),NUMPRO1(mxsprod1)
      DIMENSION DET1sd(mxsprod1),NUMPRO2(mxsprod1)
      DIMENSION DET3sd(mxsprod1)
      DIMENSION CISUM3(mxsprod1),IDET3(mxsprod1)

      DIMENSION CIH(mxsprod),CIHsuma(mxsprod)
      DIMENSION LBH(mxsprod),LBHold(mxsprod)

      DIMENSION CI3EST(mxsprod1)
      DIMENSION CI3(mxsprod1),LB33(mxsprod1)
      DIMENSION CI3suma(mxsprod1)
      DIMENSION CITsuma(mxsprod1)

      DIMENSION LB3(mxsprod1)
      DIMENSION LB2(mxsprod1)
      DIMENSION LB22(mxsprod1)

      DIMENSION CIDEST(mxsprod1)
      DIMENSION CI2EST(mxsprod1)
      DIMENSION CI2suma(mxsprod1)
      DIMENSION CIDsuma(mxsprod1)

      DIMENSION NSCF1(mxorb),NVIR2(mxorb) 
      DIMENSION IDOLD1(500),IDOLD2(500),ISCFA(mxorb)
      DIMENSION IVIRA(mxorb)

      DIMENSION ICONOS(mxorb)
c
C-----------------------------------------------------
C     compiling with 64-bit comliler
c     
C     f77 -o file.x -q64 -qintsize=8 file.f
C     file.x < input > output &
C     listout to be placed at /scr/bytautas/file.F37
C
C     Submitting a job on specific mashine
C
C     qsub -q nsfchem -eo OUT script
C-----------------------------------------------------
      dimension ifa(mxorb*mxorb)
      DIMENSION NEL(mxorb),NEL1(mxorb),ID(mxorb)
      DIMENSION IS(mxorb),IALP(mxorb)
      DIMENSION IBET(mxorb)
      DIMENSION IRREP(mxorb)

      DIMENSION KTAB(8)
      dimension iwrk(43)
      DIMENSION ISYMA(mxstr)
      DIMENSION ISYMB(mxstr)
      DIMENSION IPOSA(mxstr)
      DIMENSION IPOSB(mxstr)

c     *******
c     NOTE!!!
c     *******
c     mxstring needs to be less
c     than ICOUNTS(symmetry-reduced)
c    
      DIMENSION ISTRA(mxstring)
      DIMENSION ISTRB(mxstring)

C--------------------------------------------

       ICOUNT = 0
       ICOUNTS = 0

C      -------------------------------
       READ (5,10) NTOT 
 10    FORMAT (I5)
       READ (5,11) NOCC 
 11    FORMAT (I5)
       READ (5,20) ISYMG 
 20    FORMAT (I5)
       READ (5,21) IRRSYM 
 21    FORMAT (I5)

       write (6,*) '-----------------------------------------'
       write (6,*) ' A PRIORI PROCEDURE FOR SELECTING'
       write (6,*) ' THE CONFIGURATIONAL LIVEWOOD '
       write (6,*) ' TO GENERATE A COMPACT WAVE-FUNC.'
       write (6,*) ' FOR GAMESS RUN USING GCI CODE '
       write (6,*) '-----------------------------------------'
       write (6,*) 'PROGRAM gciselCI.f WRITTEN BY'
       write (6,*) 'LAIMIS BYTAUTAS AND KLAUS RUEDENBERG'
       write (6,*) ' '
       write (6,*) 'REFERENCE:'
       write (6,*) '(1) L. BYTAUTAS, K. RUEDENBERG, ' 
       write (6,*) '    CHEM. PHYS. 356, p.64 (2009).'
       write (6,*) '-----------------------------------------'
       write (6,*) 'no active MOs:             NTOT= ',NTOT
       write (6,*) 'no doubly occ. MOs in ref: NOCC= ',NOCC
       write (6,*) 'symmetry group is              = ',ISYMG
       write (6,*) 'irrep of the group is          = ',IRRSYM
       write (6,*) '-----------------------------------------'

  22   AAA=0
       IF (ISYMG.GT.1) GO TO 23 
       NIRR = 2
       GO TO 27
  23   AAA=0 
       IF (ISYMG.GT.2) GO TO 24 
       NIRR = 4
       GO TO 27
  24   AAA=0
       NIRR = 8
  27   AAA=0

       READ (5,*) (IRREP(J), J=1,NTOT) 
       write (6,*) 'irrep symmetries of active MOs'
       WRITE (6,*) (IRREP(J), J=1,NTOT) 

       idsym=ISYMG
       isym1=IRRSYM

       write (6,*) 'Point group symmetry',idsym
       write (6,*) 'total irrep symmetry',isym1

       READ (5,*) NKEEP3 
       READ (5,*) NKEEP4
       READ (5,*) NKEEP5
       READ (5,*) NKEEP6

       NKEEP4new = NKEEP4
       NKEEP5new = NKEEP5
       NKEEP6new = NKEEP6

       write (6,*) 'MAXIMUM # of x=3 SP   NKEEP3=',NKEEP3 
       write (6,*) 'MAXIMUM # of x=4 SP   NKEEP4=',NKEEP4 
       write (6,*) 'MAXIMUM # of x=5 SP   NKEEP5=',NKEEP5 
       write (6,*) 'MAXIMUM # of x=6 SP   NKEEP6=',NKEEP6 

c      Percentage of triple SP by weight
c      "OK ratio 3"

       read (5,*) OKRATIO3 
       SMALLrat3 = 1.0D0 - OKRATIO3
       write (6,*) '-'
       write (6,*) 'Norm3 Fraction x=3 to include =',OKRATIO3 
       write (6,*) 'Norm3 Fraction: omitted x=3   =',SMALLrat3 
       write (6,*) '-'

       read (5,*) IDIV4
       read (5,*) IDIV5
       read (5,*) IDIV6
       write (6,*) 'ISEL2 = TOTAL # of double SPACE-PRODUCTS'
       write (6,*) 'x=4 START CYCLE at (ISEL2*3) / ',IDIV4
       write (6,*) 'x=5 START CYCLE at  (ISEL2)  / ',IDIV5
       write (6,*) 'x=6 START CYCLE at  (ISEL2)  / ',IDIV6

       read (5,*) TNORM4 
       read (5,*) OKNORM4 
       read (5,*) ICYCMAX4 
       read (5,*) MDEL4 

       read (5,*) TNORM5 
       read (5,*) OKNORM5 
       read (5,*) ICYCMAX5 
       read (5,*) MDEL5 

       read (5,*) TNORM6 
       read (5,*) OKNORM6 
       read (5,*) ICYCMAX6 
       read (5,*) MDEL6 

       write (6,*) '-'
       write (6,*) '4-uple fraction selected =',TNORM4 
       write (6,*) '4-uple convergence       =',OKNORM4 
       write (6,*) '4-uple MAXIMUM no cycles =',ICYCMAX4 
       write (6,*) '4-uple CYCLE STEPSIZE    =',MDEL4 
       write (6,*) '-'
       write (6,*) '5-uple fraction selected =',TNORM5 
       write (6,*) '5-uple convergence       =',OKNORM5 
       write (6,*) '5-uple MAXIMUM no cycles =',ICYCMAX5 
       write (6,*) '5-uple CYCLE STEPSIZE    =',MDEL5 
       write (6,*) '-'
       write (6,*) '6-uple fraction selected =',TNORM6 
       write (6,*) '6-uple convergence       =',OKNORM6 
       write (6,*) '6-uple MAXIMUM no cycles =',ICYCMAX6 
       write (6,*) '6-uple CYCLE STEPSIZE    =',MDEL6 
       write (6,*) '-'

       read (5,*) ITRIP0
       read (5,*) IQUAD
       read (5,*) IQUINT
       read (5,*) ISIXT

       write (6,*) '---------------------------------------'
       write (6,*) 'THE CONVENTION:  (yes=1, no=0)'
       write (6,*) 'IS THIS PRELIMINARY CI-SDT RUN?',ITRIP0
       write (6,*) 'Are x=4 excitations included?  ',IQUAD
       write (6,*) 'Are x=5 excitations included?  ',IQUINT
       write (6,*) 'Are x=6 excitations included?  ',ISIXT
       write (6,*) '---------------------------------------'

       TNORMAL = 0.0D0
       ALLNORMAL = 0.0D0

C*    -------------------------------------
      call gtab(idsym,isym1,ktab,
     *   iwrk(1),iwrk(4),iwrk(7),iwrk(10))
c
       write (6,*) 'ktab array'
       WRITE (6,*) (KTAB(J), J=1,NIRR) 
c     -------------------------------------

       NVIRT = NTOT - NOCC 
       MALP = NOCC 
       MBET = NOCC
       NORB = NTOT

C*    ---------------------------------------------
      CALL ibinom(ifa,NTOT,NALP,NBET,NALPS,NBETS)
C*    ---------------------------------------------
c      Total # of S.P. for S, D, T, Q, 5, 6
c      excitations :
c      LIMITAS1
c      LIMITAS2
c      LIMITAS3
c      LIMITAS4
c      LIMITAS5
c      LIMITAS6
c      --------
       DO JAM=1,NTOT
       ICONOS(JAM)=0
       ENDDO
c      *
c      ALL excitations considered together
c      *
       write (6,*) 'Compute total # of SP=LIMITAS '
       LIMITAS = 0

       IEXCIT = 1 
       CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT,ICONOS,LIMITAS1)
       write (6,*) ' single-excitation SP count: ',LIMITAS1

       IEXCIT = 2 
       CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT,ICONOS,LIMITAS2)
       write (6,*) ' double-excitation SP count: ',LIMITAS2

       IEXCIT = 3 
       CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT,ICONOS,LIMITAS3)
       write (6,*) ' triple-excitation SP count: ',LIMITAS3

       IEXCIT = 4 
       CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT,ICONOS,LIMITAS4)
       write (6,*) ' quadruple-excitation SP count: ',LIMITAS4

       IEXCIT = 5 
       CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT,ICONOS,LIMITAS5)
       write (6,*) ' 5-tuple-excitation SP count: ',LIMITAS5

       IEXCIT = 6 
       CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT,ICONOS,LIMITAS6)
       write (6,*) ' 6-tuple-excitation SP count: ',LIMITAS6

       LIMITAS = LIMITAS + LIMITAS1 + LIMITAS2 + LIMITAS3
       LIMITAS = LIMITAS + LIMITAS4 + LIMITAS5 + LIMITAS6
       LIMITAS = LIMITAS + 1 

       IEXCIT1 = 1
       IEXCIT2 = 2 
       IEXCIT3 = 3 
       IEXCIT4 = 4 
       IEXCIT5 = 5 
       IEXCIT6 = 6 

       write (6,*) '-'
       write (6,*) 'TOTAL SP COUNT =',LIMITAS
       write (6,*) 'WARNING: LIMITAS sould be  <  mxstring=',mxstring
       write (6,*) '*'
       write (6,*) '(Q-,5-,6-) mxsprod =',mxsprod
       write (6,*) '(S-,D-,T-) mxsprod1=',mxsprod1
       write (6,*) '-'
       IF (LIMITAS6.LE.mxsprod) GO TO 7349
       write (6,*) 'PROBLEM with 6-tuple SP limit'
       write (6,*) 'TOTAL # of 6-tuple SP is',LIMITAS6
       write (6,*) 'THE PRESENT dimension is',mxsprod
       write (6,*) 'THUS SET LIMITAS6 = mxsprod '
       LIMITAS6 = mxsprod
c***** STOP
 7349  A=0
       IF (LIMITAS.LE.mxstring) GO TO 7350
       write (6,*) 'Problem with SP total dimension'
       write (6,*) 'THE TOTAL    # of SPs is',LIMITAS
       write (6,*) 'THE PRESENT dimension is',mxstring
       write (6,*) 'WILL SET LIMITAS = mxstring '
       LIMITAS = mxstring
c***** STOP
 7350  A=0

c      Use SDT-CI calculation data (*.12 and *.F37) 

       NA = NOCC
       NB = NOCC

       write (6,*) 'NORB  =',NORB
       write (6,*) 'NTOT  =',NTOT
       write (6,*) 'NOCC  =',NOCC
       write (6,*) 'NVIRT =',NVIRT
       write (6,*) '# alfa electrons NA=',NA
       write (6,*) '# beta electrons NB=',NB



       IF (ITRIP0.EQ.1) GO TO 7777



c      NSTATE = number of states
c      NGCI   = number of determinants


c      -------------------------------------------------
       WRITE (6,*) 'SDT-CI-vector follows:' 
C*     READ-IN CI-vector from  FILE.F12(CIvectorfile)
c      -------------------------------------------------
       open (unit=22,file='FILENAME.civec',status='unknown',
     *      access='sequential',form='unformatted')
       READ (22) NSTATEsd,NGCIsd
       READ (22) (CIsd(I), I=1,NGCIsd)
       close (22)
       write (6,*) 'SDT-CI: NSTATEsd, NGCIsd'
       write (6,*) 'No of states, No of determinants'
       write (6,*) NSTATEsd,NGCIsd

C      ---------------------------------------------------
C      READ-IN determinants-alfa/beta strings from  FILE.F37
C      ---------------------------------------------------
c      CI-SDT gcilist
       open (unit=26,file='FILENAME.gci',status='unknown',
     *      access='sequential',form='unformatted')
       READ (26) NGCIsd
       READ (26) (ICONAsd(I), I=1,NGCIsd)
       READ (26) (ICONBsd(I), I=1,NGCIsd)
       close (26)
c      *
c       write (6,*) 'SDT-CI: space-products'
c       write (6,*) 'alfa string,  beta string  positions'
c       DO JIM=1,10
c       write (6,*) JIM,ICONAsd(JIM),ICONBsd(JIM) 
c       ENDDO

c*     use SDT-wavefunction for S, D, and T-dets
c      Loop over SDT-CI determinants
c      *
       CISUMsd = 0.0D0
       NDOUBLE = LIMITAS2  
       NTRIPLE = LIMITAS3 

       DO JJJ=1,NDOUBLE
       CISUM1sd(JJJ)=0.0D0
       CISUMDsd(JJJ)=0.0D0
       IDET1sd(JJJ) = 0
       IDETDsd(JJJ) = 0
       DET1sd(JJJ) = 0.0D0
       DETDsd(JJJ) = 0.0D0
       ENDDO

       DO JJJ=1,NTRIPLE
       CISUM3sd(JJJ)=0.0D0
       IDET3sd(JJJ) = 0
       DET3sd(JJJ) = 0.0D0
       ENDDO

       write (6,*) '# of double-SPs=',NDOUBLE
       write (6,*) '# of triple-SPs=',NTRIPLE

       ISINGLE = 0
       IDOUBLE = 0
       ITRIPLE = 0
       IEXTRA = 0

       IMAXDETsd = NGCIsd 

       write (6,*) 'MAXIMUM # of determinants in SDT-CI:',NGCIsd
c      DETERMINE CONFIGURATIONS FROM STRING POSITIONS
c      Go over SDT-CI determinants:

c      DO loop cycle for translating a- & b- string
c      position in the list into electron occupations
c      (2,1,0)
c      Every determinant (a- & b- string) will
c      be considered in DO loop


       DO 5350 K=1,IMAXDETsd
c      *
c      *
       CIsd(K) = dabs(CIsd(K))
       IPOSAsd = ICONAsd(K)
       IPOSBsd = ICONBsd(K)
c      *
       DO I=1,NA
       IACONsd(I) = I
       ENDDO
c      
       DO I=1,NB
       IBCONsd(I) = I
       ENDDO
c      
       IXAsd = IPOSAsd - 1
       IXBsd = IPOSBsd - 1

       IF (IXAsd.LT.1) GO TO 5402
c      "alfa electrons"
c      -
       DO I=1,IXAsd
       call advanc(IACONsd,NA,NORB)
       ENDDO
c      *
 5402  A=0.0
       IF (IXBsd.LT.1) GO TO 5409
c      "beta" electrons
c      -
       DO I=1,IXBsd
       call advanc(IBCONsd,NB,NORB)
       ENDDO
c      *
 5409  A=0.0

c------------------------------------------------
c      ORBITALS OCCUPIED BY ALFA/BETA ELECTRONS:
c       write (6,*) 'SDT-CI: alfa electron MOs'
c       write (6,*) (IACONsd(i), i=1,NA)
c       write (6,*) 'SDT-CI: beta electron MOs'
c       write (6,*) (IBCONsd(i), i=1,NB)
c------------------------------------------------

       DO I=1,NORB
c      ---------------
c      *

       IADD = 0

c      *
       DO J=1,NA
c      ---------------
       IF (IACONsd(J).EQ.I) GO TO 5425
c      -------------------------------
       ENDDO
       GO TO 5430

 5425  A=0.0
c      -----------------
       IADD = IADD + 1
c      -----------------
c      *
 5430  A=0.0
c      *
       DO JJ=1,NB
c      --------------
       IF (IBCONsd(JJ).EQ.I) GO TO 5440
c      --------------
       ENDDO
       GO TO 5445
c
 5440  A=0.0
c      -----------------
       IADD = IADD + 1
c      -----------------
c      *
 5445  A=0.0
c      *
c      --------------------
       ICONsd(I) = IADD
c      --------------------
c      *
       ENDDO

C--------------------------------------------------
c      CONFIGURATION HAS BEEN IDENTIFIED:
C      ICONsd(i), i=1,NORB
c      ICOUNT=ISCF  - number of SCFMOs = reference
c      IVIR=ICOUNTS - number of VIR MOs
c--------------------------------------------------

c     Assign the excitation type x: 
c     "Space-Products (SP) will be identified
c     and their weights and importances
c     will be deduced from CI-SDT determinants
C     NOTE: usually there are several dets
c     (differing by their spin-couplings) contributing 
c     to a given SP
c     -
c     Assign the excitation type x: 
c     -
      CALL INFO(NTOT,NOCC,NVIRT,ICONsd,
     *NSCF,ISCF,IOCC,NVIR,IVIR,IOCCV,IEXCIT) 
c     -
      IEXCITsd(K) = IEXCIT 
c     -
      IF (IEXCITsd(K).EQ.0) GO TO 5731
      IF (IEXCITsd(K).NE.1) GO TO 5725

c     "IEXCIT=1"
      ISINGLE = ISINGLE + 1

      IEXCsd = IEXCITsd(K) 
c     -
c     here will label the space-product by assigning
c     a number "NSPRO1sd" ("1"-for singles)
c     IDET1sd = # of dets in the given SP
c
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCsd,ICONsd,NSPRO1sd)
c     *
      CISUM1sd(NSPRO1sd) = CISUM1sd(NSPRO1sd) + (CIsd(K)**2)
      IDET1sd(NSPRO1sd) = IDET1sd(NSPRO1sd) + 1 
      DET1sd(NSPRO1sd) = DET1sd(NSPRO1sd) + 1.0D0
c     *
c*    write (6,*) (ICONsd(I), I=1,NORB),'CIsd=',CIsd(K)
c*    write (6,*) 'NSPRO1sd=',NSPRO1sd
c     *
      GO TO 5731


 5725 A=0
      IF (IEXCITsd(K).NE.2) GO TO 5727
c     ---------------------------------
c     IEXCIT=2
      IDOUBLE = IDOUBLE + 1
c
      IEXCsd = IEXCITsd(K) 

c     here will label the space-product by assigning
c     a number "NSPRODsd" ("D"-for doubles)
c     IDETDsd = # of dets in the given SP

      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCsd,ICONsd,NSPRODsd)
c     *
      CISUMDsd(NSPRODsd) = CISUMDsd(NSPRODsd) + (CIsd(K)**2)
      IDETDsd(NSPRODsd) = IDETDsd(NSPRODsd) + 1 
      DETDsd(NSPRODsd) = DETDsd(NSPRODsd) + 1.0D0
c     *
c*    write (6,*) (ICONsd(I), I=1,NORB),'CIsd=',CIsd(K)
c*    write (6,*) 'NSPRODsd=',NSPRODsd
      GO TO 5731


 5727 A=0
      IF (IEXCITsd(K).NE.3) GO TO 5729
c     ---------------------------------
c     IEXCIT=3
      ITRIPLE = ITRIPLE + 1
c
      IEXCsd = IEXCITsd(K) 

c     here will label the space-product by assigning
c     a number "NSPRO3sd" ("3"-for triples)
c     IDET3sd = # of dets in the given SP

      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCsd,ICONsd,NSPRO3sd)
c     *
      CISUM3sd(NSPRO3sd) = CISUM3sd(NSPRO3sd) + (CIsd(K)**2)
      IDET3sd(NSPRO3sd) = IDET3sd(NSPRO3sd) + 1 
      DET3sd(NSPRO3sd) = DET3sd(NSPRO3sd) + 1.0D0
c     *
c     write (6,*) (ICONsd(I), I=1,NORB),'CIsd=',CIsd(K)
c     write (6,*) 'NSPRO3sd=',NSPRO3sd
c     write (6,*) 'IEXCsd=',IEXCsd
      GO TO 5731

 5729 A=0
      IEXTRA = IEXTRA + 1
      write (6,*) 'excitations >3, PROGRAM WILL STOP'
      write (6,*) (ICONsd(I), I=1,NORB),'CIsd=',CIsd(K)
      write (6,*) 'IEXCITsd(K)=',IEXCITsd(K)
      STOP
 5731 A=0
c      --------------------------------
       CISUMsd = CISUMsd + (CIsd(K)**2)
c      --------------------------------
c      *
 5350  CONTINUE

       write (6,*) '------------------------------------------'
       write (6,*) 'SUMMARY OF EXCITATIONS IN TERMS OF DETS:'
       write (6,*) '------------------------------------------'
       write (6,*) 'ISINGLE=',ISINGLE
       write (6,*) 'IDOUBLE=',IDOUBLE
       write (6,*) 'ITRIPLE=',ITRIPLE
       write (6,*) 'IEXTRA=',IEXTRA
       write (6,*) 'SUM OF (c**2) in SDT-CI: ',CISUMsd
       write (6,*) 'END of SDT-CI cycle'
       write (6,*) '------------------------------------------'


       ZEROsmall = 0.00000000000001D0
       ZERObig = 0.00000001D0


c      THIS STEP DETERMINES WEIGHTS & IMPORTANCES
c      for SINGLES, DOUBLES & TRIPLES:

c      *singles*
c      ----------
       DO JJJ=1,NDOUBLE
c***** DO JJJ=1,LIMITAS1
c      -----------------
       IF (DET1sd(JJJ).LT.ZERObig) GO TO 5600
       IF (CISUM1sd(JJJ).LT.ZEROsmall) GO TO 5600
       CISUM1sd(JJJ) = CISUM1sd(JJJ) / DET1sd(JJJ)
 5600  A=0
       ENDDO


       ISEL2 = 0
c      ----------
c      *doubles*
c      DO loop over ALL POSSIBLE DOUBLE SPs
c      ----------
       DO JJJ=1,NDOUBLE
c***** DO JJJ=1,LIMITAS2
c      -----------------
       IF (DETDsd(JJJ).LT.ZERObig) GO TO 5610
       IF (CISUMDsd(JJJ).LT.ZEROsmall) GO TO 5610
c      *
       CIDsd(JJJ) = CISUMDsd(JJJ) 
       CISUMDsd(JJJ) = CISUMDsd(JJJ) / DETDsd(JJJ)
c      ***
       ISEL2 = ISEL2 + 1
       LB22(ISEL2) = JJJ
       CIDsuma(ISEL2) = CIDsd(JJJ)
       CIDEST(ISEL2) = CISUMDsd(JJJ)
 5610  A=0
       ENDDO

c     FINISHED estimating SPs, now SORT IT OUT:
c     -------------------------------------
      call SORT02(ISEL2,CIDEST,CIDsuma,LB22)
c     -------------------------------------

      DO JJJ=1,ISEL2
      J = ISEL2 + 1 - JJJ
      CI2EST(JJJ) = CIDEST(J)
      CI2suma(JJJ) = CIDsuma(J)
      LB2(JJJ) = LB22(J)
      ENDDO


      CI2sumTOT = 0.0D0
      DO J=1,ISEL2
      CI2sumTOT = CI2sumTOT + CI2suma(J)
      ENDDO


      write (6,*) '*'
      write (6,*) 'NUMBER OF SELECTED DOUBLE SPACE-PRODs=',ISEL2
      write (6,*) 'TOTAL SUM of DOUBLE-SP c**2 =',CI2sumTOT
      write (6,*) '*'


       MAX4IM2 = 0
       IF (MAX4IM2.GT.1) GO TO 777
       MAX4IM2 = (ISEL2*3)/IDIV4 
       MAX5IM2 = ISEL2/IDIV5 
       MAX6SP2 = ISEL2/IDIV6 
       write (6,*) 'new starting DOUBLES-values:' 
       write (6,*) 'FOR QUADRUPLES: MAX4IM2 =',MAX4IM2 
       write (6,*) 'FOR QUINTUPLES: MAX5IM2 =',MAX5IM2 
       write (6,*) 'FOR SEXTUPLES:  MAX6SP2 =',MAX6SP2 
 777   A=0

       write (6,*) 'OPEN-ENDED CYCLE for excitations x=4,5,6'
       write (6,*) '4-uple open-ended cycle'
       write (6,*) '             starts at Double-SP =',MAX4IM2 
       write (6,*) '-'
       write (6,*) '5-uple open-ended cycle '
       write (6,*) '             starts at Double-SP =',MAX5IM2 
       write (6,*) '-'
       write (6,*) '6-uple open-ended cycle '
       write (6,*) '             starts at Double-SP =',MAX6SP2 


       ISEL3 = 0
c      ----------
c      *triples*
c      ----------
       DO JJJ=1,NTRIPLE
       IF (DET3sd(JJJ).LT.ZERObig) GO TO 5620
       IF (CISUM3sd(JJJ).LT.ZEROsmall) GO TO 5620
       CI3coef = CISUM3sd(JJJ)
       CI3sd(JJJ) = CISUM3sd(JJJ) / DET3sd(JJJ)
       ISEL3 = ISEL3 + 1
       LB3(ISEL3) = JJJ
       CI3EST(ISEL3) = CI3sd(JJJ)
       CI3suma(ISEL3) = CI3coef
 5620  A=0
       ENDDO





c     ORDER THE TRIPLE SPs by the IMPORTANCE:
c     -----------------------------------------
      call SORT02(ISEL3,CI3EST,CI3suma,LB3)
c     -----------------------------------------

      DO JJJ=1,ISEL3
c     -----------------
      J = ISEL3 + 1 - JJJ
      CI3(JJJ) = CI3EST(J)
      CITsuma(JJJ) = CI3suma(J)
      LB33(JJJ) = LB3(J)
      ENDDO

      CI3sumTOT = 0.0D0

      DO J=1,ISEL3
      CI3EST(J) = CI3(J)
      CI3suma(J) = CITsuma(J)
      LB3(J) = LB33(J)
      CI3sumTOT = CI3sumTOT + CI3suma(J)
      ENDDO


      write (6,*) 'Triple non-zero sp-prod ISEL3=',ISEL3
      write (6,*) 'TOTAL triple c**2 sum =',CI3sumTOT
c     *




      IF (NKEEP3.LE.ISEL3) GO TO 5057
c     -------------------------------
      write (6,*) 'NOTE: NKEEP3=ISEL3'
      NKEEP3 = ISEL3
c     -------------------------------
 5057 A=0 


      CI3def = CI3sumTOT


      write (6,*) 'deficiency estimate for T '
c     ----------------------------------------
      DO J=1,NKEEP3
      CI3def = CI3def - CI3suma(J)
      ENDDO

      RELDEF3 = CI3def / CI3sumTOT
      
      write (6,*) '-----------------------------------------------'
      write (6,*) '                            ISEL3 =',ISEL3
      write (6,*) '                           NKEEP3 =',NKEEP3
      write (6,*) 'TOTAL Triple-deficiency  CI3def   =',CI3def
      write (6,*) 'RELATIVE Triple-deficiency RELDEF3=',RELDEF3
      write (6,*) '-----------------------------------------------'



c           **************************************
c     ***** SELECTED PERCENT-TYPE TRIPLE SELECTION ********
c           **************************************
      CI3def = CI3sumTOT
      CI3norm = 0.0D0
      CI3norm2 = 0.0D0
      ICUT3 = 0

      write (6,*) 'T-selection based on norm3'


      DO J=1,ISEL3
c     -----------------
c     *
      ICUT3 = ICUT3 + 1
      CI3norm2 = CI3norm2 + CI3suma(J)
      RATIO3 = (CI3sumTOT - CI3norm2)
      RATIO3 = RATIO3 / CI3sumTOT 
      IF (RATIO3.LE.SMALLrat3) GO TO 5095
      ENDDO

      write (6,*) 'Will stop. T-ple selection: '
      write (6,*) 'cannot assign the percentage !!! '
      STOP
 5095 A=0
c     ***


      write (6,*) '-----------------------------------------'
      write (6,*) 'NORM3 FRACTION OF EXCLUDED SP3 =',RATIO3
      write (6,*) '                         ICUT3 =',ICUT3
      write (6,*) '-----------------------------------------'



c     --------------
      NKEEP3 = ICUT3
c     --------------

c     ----------------
c     IMPOSE THIS AT THE END:
c     Select all SP3:
c     NKEEP3 = ISEL3
c     --------------

      write (6,*) '3-ple: Finally: out of  ISEL3 =',ISEL3
      write (6,*) '3-ple: keep            NKEEP3 =',NKEEP3
      write (6,*) '---------------------------------------------'

      CI3def = CI3sumTOT
      DO J=1,NKEEP3
      CI3def = CI3def - CI3suma(J)
      CI3norm = CI3norm + CI3suma(J)
      ENDDO
      RELDEF3 = CI3def / CI3sumTOT
      REL3norm = CI3norm / CI3sumTOT

      write (6,*) 'Triples will have         SP3 = ',NKEEP3
      write (6,*) '---------------------------------------------'
      write (6,*) '3-norm                CI3norm =',CI3norm
      write (6,*) 'RELATIVE 3-norm      REL3norm =',REL3norm
      write (6,*) '---------------------------------------------'
      write (6,*) '3-deficiency        CI3def    =',CI3def
      write (6,*) 'RELATIVE 3-deficiency RELDEF3 =',RELDEF3
      write (6,*) '---------------------------------------------'
      write (6,*) 'S,D,T space products done from SDT-CI'
      write (6,*) '---------------------------------------------'



c**** write (6,*) 'Now proceed with the SDTQ56-CI selection scheme'
c      THE FOLLOWING CODE GENERATES SPACE-PRODUCTS
c      FOR HIGHER EXCITATIONS: x=4,5 and 6
c      MAIN DO LOOP:

 7777  A=0



C      one-determinant (SCF) :
C      -----------------------
       DO 30 I=1,NOCC
       NEL(I) = 2
 30    CONTINUE
       DO 32 I=NOCC+1,NTOT
       NEL(I) = 0
 32    CONTINUE
c*     WRITE (6,*) (NEL(I), I=1,NTOT) 
       

       ICOUNT = ICOUNT + 1
       ICNT = 1
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
C*write (6,*) 'positions of alpha/beta det-strings in the list'
C*write (6,*) IPOSA(ICOUNT),IPOSB(ICOUNT) 

       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 

C*write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
C*write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)

       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 33
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
C*     WRITE (6,*) (NEL(KLL), KLL=1,NTOT) 
 33    AAA=0



       write (6,*) 'single excitations start here'
C      single excitations:
C      -------------------
       DO 35 I=1,NOCC
c      --------------
       DO 36 J=1,NOCC
       NEL(J)= 2
  36   CONTINUE
       NSING = NOCC+1-I
       NEL(NSING) = 1

c      virtual
       DO 39 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 39    CONTINUE
       DO 40 MM=NOCC+1,NTOT
       NEL(MM) = 1 

c***** WRITE (6,*) (NEL(L), L=1,NTOT) 

c      --------------------
       ICOUNT = ICOUNT + 1
       ICNT = 1
C*     --------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
C*      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
C*     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
C*     write (6,*) 'positions of alpha/beta det-strings in the list'
C*     write (6,*) IPOSA(ICOUNT),IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 

C* write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
C* write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)

       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 44 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 44    AAA=0
 43    AAA=0
       NEL(MM) = 0
 40    CONTINUE
 37    AAA=0
 35    CONTINUE

       write (6,*) 'SINGLE-EX. SPs symmetry-unconstrained',ICOUNT
       write (6,*) 'SINGLE-EX. SPs symmetry-constrained',ICOUNTS


C      double excitations:
       write (6,*) 'DIAGONAL-double excitations start here'
       DO 45 I=1,NOCC
c      -------------------
c      doubly-occupied-space
       DO 47 J=1,NOCC
       NEL(J)= 2
  47   CONTINUE
       NSING = NOCC+1-I
       NEL(NSING) = 0 

c      virtual-space
       DO 49 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 49    CONTINUE
       DO 50 MM=NOCC+1,NTOT
       NEL(MM) = 2 


       ICOUNT = ICOUNT + 1
       ICNT = 1
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
C*    CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
C*   *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)

C*     write (6,*) 'positions of alpha/beta det-strings in the list'
C*     write (6,*) IPOSA(ICOUNT),IPOSB(ICOUNT) 

       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
C* write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
C* write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 54 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 54    AAA=0
 53    AAA=0
       NEL(MM) = 0
 50    CONTINUE


c      virtual
       DO 55 JJJ=1,NVIRT
c      ------------------
       DO 56 KKK=NOCC+1,NTOT
       NEL(KKK) = 0
  56   CONTINUE
       M = NTOT+1-JJJ
       NEL(M) = 1
       DO 60 II=1,M-NOCC-1
c      -------------------
       MMM = M - II
       NEL(MMM) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
C*CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
C**NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
C*     write (6,*) positions of alpha/beta det-strings in the list'
C*     write (6,*) IPOSA(ICOUNT),IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
C*       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
C*       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 64 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 64    AAA=0
 63    AAA=0
       NEL(M-II) = 0
 60    CONTINUE

 59    AAA=0.0
       NEL(M) = 0
 55    CONTINUE

 48    AAA=0.0
 45    CONTINUE



c      --------------------------
c      S and D diagonal (MOPAIRS)
c      --------------------------
       WRITE (6,*) 'S[1] and D[0] ICOUNT=',ICOUNT
       write (6,*) 'mixed double excitations start here'
C      double excitations-mixed:
C      -------------------


       DO 275 I=1,NOCC
       DO 277 J=1,NOCC
       NEL(J)= 2
 277   CONTINUE
       M1 = NOCC+1-I 
       NEL(M1) = 1 
       DO 278 K1=1,M1-1
       M2 = M1 - K1
       NEL(M2)= 1

c      virtual-space
       DO 280 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 280   CONTINUE
       DO 285 MM=NOCC+1,NTOT
       NEL(MM) = 2 

       ICOUNT = ICOUNT + 1
       ICNT = 1
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
C*  CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
C*  *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
C*       WRITE (6,291) 
C* 291   FORMAT ('positions of alpha/beta det-strings in the list')
C*       WRITE (6,292) IPOSA(ICOUNT),IPOSB(ICOUNT) 
C* 292   FORMAT (I5,2X,I5)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
C*       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
C*       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 295 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 295   AAA=0
 288   AAA=0
       NEL(MM) = 0
 285   CONTINUE



       DO 300 JJJ=1,NVIRT
       DO 301 KKK=NOCC+1,NTOT
       NEL(KKK) = 0
 301   CONTINUE
       M = NTOT+1-JJJ
       NEL(M) = 1
       DO 305 II=1,M-NOCC-1
       MM = M - II
       NEL(MM) = 1

       ICOUNT = ICOUNT + 1
       ICNT = 1
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
C*      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
C*     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
C*      write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
C*      write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 352 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 352   AAA=0
 308   AAA=0
       NEL(MM) = 0
 305   CONTINUE
 304   AAA=0.0
       NEL(M) = 0
 300   CONTINUE
       NEL(M1-K1) = 2
 355   AAA=0.0
 278   CONTINUE
       NEL(M1) = 2
 354   AAA=0.0
 275   CONTINUE
       write (6,*) 'DOUBLE-EXC. SPs symm-unconstrained:',ICOUNT
       write (6,*) 'DOUBLE-EXC. SPs symm-constrained:',ICOUNTS
 359   ALFA=0








       write (6,*) 'triple excitations start here'
C      -------------------
C      TRIPLE EXCITATIONS:
C      -------------------
       ICOUNT3 = 0


c------------------------------------------
c      warning:
c      *******
c      ITRIP0 = 1 ( no T-selection)
c      ITRIP0 = 0 (selection NKEEP3)
c------------------------------------------

       write (6,*) 'T-case: [1,0]'
C      occupied MO space: CASE-A)  [1,0]
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO 363 KKZ=1,NOCC
       NXZ = NOCC + 1 - KKZ
       NEL(NXZ) = 0
       DO 367 KK1Z=1,NOCC-1
       M1Z = NOCC + 1 - KK1Z 
       IF (M1Z.GT.NXZ) GO TO 371 
       IF (M1Z.LE.NXZ) GO TO 372 
 371   AAA=0
       L1Z = M1Z
       GO TO 373 
 372   AAA=0
       L1Z = M1Z - 1
 373   AAA=0
       NEL(L1Z) = 1

c      Virtual MOs
c      CASE-1 [2,1]
       DO 380 K = NOCC+1,NTOT
       NEL(K) = 0
 380   CONTINUE


       DO 385 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2
       DO 389 KK1=1,NVIRT-1
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 390 
       IF (M1.LE.NX) GO TO 391
 390   AAA=0
       L1 = M1
       GO TO 392
 391   AAA=0
       L1 = M1 - 1
 392   AAA=0
       NEL(L1) = 1


       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT3 = ICOUNT3 + 1

       IF (ITRIP0.EQ.1) GO TO 397 
       CALL FILTER3(NKEEP3,LB3,ICOUNT3,IPUT)
       IF (IPUT.EQ.0) GO TO 396
 397   A=0
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
c      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 396 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 396   AAA=0
       NEL(L1) = 0
 389   CONTINUE
       NEL(NX) = 0
 385   CONTINUE


c      virtual MOs
c      CASE-2) [1,1,1]
c
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO 398 KK = 1,NVIRT-2
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 402 KK1 = 1,M1-NOCC-2
       M2 = M1 - KK1
       NEL(M2) = 1

       DO 406 KK2 = 1,M2-NOCC-1
       M3 = M2 - KK2
       NEL(M3) = 1

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT3 = ICOUNT3 + 1
       IF (ITRIP0.EQ.1) GO TO 409 
       CALL FILTER3(NKEEP3,LB3,ICOUNT3,IPUT)
       IF (IPUT.EQ.0) GO TO 410 
  409  A=0
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
c      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 410 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 410   AAA=0
       NEL(M3) = 0
 406   CONTINUE
       NEL(M2) = 0
 402   CONTINUE
       NEL(M1) = 0
 398   CONTINUE
c      ---------------------
c      end of occupied MOs: CASE-A)
c      ---------------------
       NEL(L1Z) = 2
 412   AAA=0.0 
 367   CONTINUE
       NEL(NXZ) = 2
 411   AAA=0.0
 363   CONTINUE




       write (6,*) 'T-case: [1,1,1]'
c      occupied MOs CASE-B) [1,1,1]
c      ----------------------------
       DO 414 I= 1,NOCC
       NEL(I) = 2
 414   CONTINUE
       DO 415 KKW = 1,NOCC-2
       M1W = NOCC+1-KKW
       NEL(M1W) = 1
       DO 418 KK1W = 1,M1W-2
       M2W = M1W - KK1W
       NEL(M2W) = 1
       DO 421 KK2W = 1,M2W-1
       M3W = M2W - KK2W
       NEL(M3W) = 1

c      Virtual MOs
c      CASE-1 [2,1]
c      ---------------------- 
       DO 424 K = NOCC+1,NTOT
       NEL(K) = 0
 424   CONTINUE
       DO 425 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2
       DO 430 KK1=1,NVIRT-1
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 431 
       IF (M1.LE.NX) GO TO 432
 431   AAA=0
       L1 = M1
       GO TO 433
 432   AAA=0
       L1 = M1 - 1
 433   AAA=0
       NEL(L1) = 1

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT3 = ICOUNT3 + 1
       IF (ITRIP0.EQ.1) GO TO 439 
       CALL FILTER3(NKEEP3,LB3,ICOUNT3,IPUT)
       IF (IPUT.EQ.0) GO TO 440 
 439   A=0
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
c      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 440 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 440   AAA=0
       NEL(L1) = 0
 430   CONTINUE
       NEL(NX) = 0
 425   CONTINUE

c      virtual MOs
c      CASE-2) [1,1,1]
c      ---------------------- 
       DO 441 K = NOCC+1,NTOT
       NEL(K) = 0
 441   CONTINUE
       DO 442 KK = 1,NVIRT-2
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 446 KK1 = 1,M1-NOCC-2
       M2 = M1 - KK1
       NEL(M2) = 1
       DO 450 KK2 = 1,M2-NOCC-1
       M3 = M2 - KK2
       NEL(M3) = 1

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT3 = ICOUNT3 + 1
       IF (ITRIP0.EQ.1) GO TO 453 
       CALL FILTER3(NKEEP3,LB3,ICOUNT3,IPUT)
       IF (IPUT.EQ.0) GO TO 454
 453   A=0
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
c      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 454 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 454   AAA=0
       NEL(M3) = 0
 450   CONTINUE
       NEL(M2) = 0
 446   CONTINUE
       NEL(M1) = 0
 442   CONTINUE
c      -------------------------------
c      end of occupied MOs CASE-B) [1,1,1]
c      -------------------------------
       NEL(M3W) = 2 
 457   AAA=0.0
 421   CONTINUE
       NEL(M2W) = 2 
 456   AAA=0.0
 418   CONTINUE
       NEL(M1W) = 2 
 455   AAA=0.0
 415   CONTINUE

       write (6,*) '------------------------------------------------'
       write (6,*) 'T-t:                           NKEEP3 =',NKEEP3
       write (6,*) 'T-t:                           ICOUNT =',ICOUNT
       write (6,*) 'T-t:                          ICOUNT3 =',ICOUNT3
       write (6,*) 'T-t: symmetry-reduced strings ICOUNTS =',ICOUNTS
       write (6,*) 'CI-SDT done'
       write (6,*) '------------------------------------------------'



c      ---------------------------
       IF (ITRIP0.EQ.1) GO TO 1001 
c      ---------------------------
       IF (IQUAD.EQ.0) GO TO 1000 
c      ---------------------------


  490  AAA=0.0


       write (6,*) '*'
       write (6,*) 'quadruple excitations start here'
       write (6,*) '*'

c**** MAXSP4 = 250000 
c---------------------------------------------------------
c              DEFINITIONS:
c
c     MAXSP4 = # of SP4 components
c---------------------------------------------------------


      MAXSP4 = 3 * NKEEP4new 
      SMALL4 = 0.0000000000005D0

      write (6,*) '---'
      write (6,*) 'maximum components SP4 =',MAXSP4
      write (6,*) '        4-tuple SMALL4 =',SMALL4
      write (6,*) '---'


c     ANGLE(2-2): ISTEP3 = 1 
c     **********************
      ISTEP3 = 1 
      STEP3 = 1.0D0 
c     -----------
c     *
      write (6,*) '---'
      write (6,*) '4-tuple selection: lines-angle'
      write (6,*) 'ISTEP3 =',ISTEP3
      write (6,*) ' STEP3 =',STEP3
      write (6,*) '---'
c     ------------------------------------------------------
c      TNORM4  = CHOSEN FRACTION of QUADRUPLE SPs 
c      OKNORM4 = LIVE-WOOD SELECTION CONVERGENCE CRITERION 
c     ------------------------------------------------------


      write (6,*) '-------------------------------------------'
      write (6,*) 'THE QUADRUPLE SP "LIVE-WOOD" SELECTION PROCEDURE:'
      write (6,*) 'THE CONVERGENCE CRITERION IS '
      write (6,*) 'NORM(i)/NORM(i+1) =',OKNORM4
      write (6,*) 'OUT OF THE MOST IMPORTANT QUADRUPLE SP' 
      write (6,*) 'ONLY THIS FRACTION WILL BE INCLUDED =',TNORM4
      write (6,*) '-------------------------------------------'


c------------------------------------------------
c     OPEN-ENDED APPROACH IMPLEMENTED FOR
c     Q-UPLES WITH ENLARGING THE TRIANGLE
c     TO IMPROVE THE SELECTION PROCEDURE FOR SP4
c
c     HERE: MDEL4(INCREMENT STEPS)
c     ICYCMAX4-1 = No of INCREMENT CYCLES
c------------------------------------------------



      CI4oldTOT = 0.0D0 


      IFIX4 = 1

      DO JJJ=1,LIMITAS4
c     -----------------
c**** CI4(JJJ) = 0.0D0
c**** CI4suma(JJJ) = 0.0D0
      CIH(JJJ) = 0.0D0
      CIHsuma(JJJ) = 0.0D0
      LBH(JJJ) = 0
      LBHold(JJJ) = 0
      ENDDO

      DO 491 ICYC4 = 1,ICYCMAX4
c     --------------------------
      write (6,*) '*'
      write (6,*) '*'
      write (6,*) '*'
      write (6,*) '                   ********************************'
      write (6,*) '                           ICYC4 =',ICYC4
      write (6,*) '                   ********************************'

      NKEEP4 = NKEEP4new
      INUMB4 = 0
      ISEL4 = 0
      ITOT4 = 0


      IF (ICYC4.EQ.1) GO TO 492 
      MAX4IM2 = MAX4IM2 + MDEL4
 492  A=0


      DO JJJ=1,LIMITAS4
c     -----------------
      CIH(JJJ) = 0.0D0
      CIHsuma(JJJ) = 0.0D0
      LBH(JJJ) = 0
      ENDDO



      IF (MAX4IM2.LE.ISEL2) GO TO 4057
      write (6,*) '---'
      write (6,*) 'MAX4IM2 =',MAX4IM2 
      write (6,*) 'NOTE: 4-tuple: MAX4IM2 Exceeded ISEL2=',ISEL2
      write (6,*) '---'
C--------------------------------------------------------------
C     THIS MEANS THAT GOING ALONG THE DIAGONAL, ONE
c     WILL ENCOUNTER REGIONS ALONG THE LINE WITH NO
c     D-VALUES EXISTING, HOWEVER EVENTUALLY THEY
c     MUST START REAPPEARING (AT HIGH J values) SINCE
c     D*D non-zero values form a square, not a triangle
c     For very high D-values the square will be enclosed,
c     and no further contributions will occur.
c     LAIMIS, JANUARY 16, 2007
c--------------------------------------------------------------
 4057 A=0

      EMAX4IM2 = MAX4IM2 + 1
      EMAX4IM2 = EMAX4IM2 - (1.0D0/STEP3)
      write (6,*) 'Q-LIVEWOOD EST. begins at Double-SP =',EMAX4IM2



      JD1 = 1 
c     *
c-------------------------------------------------
c     NUMBER OF D-values will go from 1 to MAX4IM2
c-------------------------------------------------

      DO 4013 JD = JD1,MAX4IM2

c*    *** number of lines ***
c     -----------------------
      DO 4014 LIN = 1,ISTEP3
c     -------------------------
c     *
      KDD = JD + 1

c     *** on each line units ***
      DO 4015 JT = 1,JD 
      KTT = ((JT-1)*ISTEP3) + 1
      KTT = KTT + LIN - 1

      IF (KTT.LE.ISEL2) GO TO 4017
      GO TO 4044 
 4017 A=0
      KDD = KDD - 1

      IF (KDD.GT.ISEL2) GO TO 4040 
      IF (KDD.GE.1) GO TO 4018
      write (6,*) 'WARNING: will stop KDD =',KDD 
      STOP
 4018 A=0

      KD1 = KDD 
      KD2 = KTT 

      I2TRIA = KD1 
      I22TRIA = KD2 
      IDTRIANG = JD
      ILIN = LIN
      ITOT4 = ITOT4 + 1 
      LAB2 = LB2(KD1)
      LAB3 = LB2(KD2)

      CALL RLABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,
     *ICON2,LAB2)

      CALL RLABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,
     *ICON3,LAB3)

      CALL NUM555(NTOT,NOCC,NVIRT,ICON2,ICON3,
     *ICON4,IPUT)

c     ---------------------------
      IF (IPUT.EQ.0) GO TO 4040 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT4,ICON4,
     *IDNUM4)
c     ----------------------------------
c      IDNUM4 = Q-string position
c     ----------------------------------
      INUMB4 = INUMB4 + 1
c     --------------------------
      IF (ISEL4.LE.NKEEP4new) GO TO 4052
c     ----------------------------------
      write (6,*) '---'
      write (6,*) 'WARNING (x=4):'
      write (6,*) 'ISEL4 > NKEEP4new'
      write (6,*) '                ISEL4 =',ISEL4
      write (6,*) '       SP4 NKEEP4new  =',NKEEP4new
      write (6,*) 'MUST INCREASE DIM = NKEEP4'
      write (6,*) '---'
      GO TO 4050
 4052 A=0
c     --------------------------------
      IF (INUMB4.LE.MAXSP4) GO TO 4053
c     --------------------------------
      write (6,*) '---'
      write (6,*) 'WARNING (x=4):'
      write (6,*) '4-tuple: components    IMUMB4 =',INUMB4
      write (6,*) 'exceeded maximum comp. MAXSP4 =',MAXSP4
      write (6,*) '       MAXSP4 = 3 * NKEEP4new'
      write (6,*) 'MUST INCREASE NKEEP4 = NKEEP4new'
      write (6,*) 'x=4: cycle ICYC4 now ends'
      write (6,*) '---'
      GO TO 4050
 4053 A=0


      PD4 = (CI2EST(KD1) * CI2EST(KD2))
      PD4suma = (CI2suma(KD1) * CI2suma(KD2))


      IF (INUMB4.EQ.1) GO TO 4035

c     check if IDNUM4 has occured before:
c     -----------------------------------
      DO ITT=1,ISEL4
      IHIT = ITT
      IF (IDNUM4.NE.LBH(ITT)) GO TO 4030 
      CIH(IHIT) = CIH(IHIT) + PD4 
      CIHsuma(IHIT) = CIHsuma(IHIT) + PD4suma 
      GO TO 4040
 4030 A=0
      ENDDO
 4035 A=0

c     (4-tuple is new):
c     -----------------
      ISEL4 = ISEL4 + 1
      LBH(ISEL4) = IDNUM4
      CIH(ISEL4) = CIH(ISEL4) + PD4 
      CIHsuma(ISEL4) = CIHsuma(ISEL4) + PD4suma 
c-------------------------------------------
c*****LAIMIS:ISEL4<mxsprod
c-------------------------------------------
      IF (ISEL4.LT.mxsprod) GO TO 4040
      write (6,*) 'PROBLEM: ISEL4 is ',ISEL4
      write (6,*) '       mxsprod is ',mxsprod
      write (6,*) 'NOW PROGRAM WILL STOP '
      STOP
 4040 A=0

c     did KD2=KTT reach the end at ISEL2?
      IF (KTT.EQ.ISEL2) GO TO 4044
 4015 CONTINUE
c     *** (UNITS)
 4044 A=0
      IF (KDD.LT.1) GO TO 4051
 4014 CONTINUE
c     --------
c     *** (LINES)
 4051 A=0
 4013 CONTINUE
c     --------
c     (do loop starting at new D1 doubles)


 4050 A=0


c     ***************************
c     4-tuple selection completed
c     ***************************
      write (6,*) '-------------------------------------------'
      write (6,*) '                           ICYC4 =',ICYC4 
      write (6,*) '   Total (D*D) products   ITOT4  =',ITOT4
      write (6,*) 'Number of SP4 components  INUMB4 =',INUMB4
      write (6,*) 'Number of distinct SP4     ISEL4 =',ISEL4
      write (6,*) '                          NKEEP4 =',NKEEP4
      write (6,*) '-------------------------------------------'
c**** write (6,*) 'Last D1 :',I2TRIA 
c**** write (6,*) 'Last D2 :',I22TRIA 
c     -------------------------------
      IF (ISEL4.GE.NKEEP4) GO TO 4079
c     -------------------------------
      NKEEP4 = ISEL4

      write (6,*) '                      PUT NKEEP4 = ISEL4'
      write (6,*) '                          NKEEP4 =',NKEEP4
 4079 A=0


c            -------------------
c     ****** ORDER 4-tuples now: *****
c            -------------------
      IF (ISEL4.GE.IFIX4) GO TO 4081 
      write (6,*) 'ISEL4 is too small'
      STOP
 4081 A=0
c     --------------------------
      IF (ICYC4.GT.1) GO TO 4082 
c     --------------------------
      IFIX4 = ISEL4
      SUM4old = 0.0D0
      DO JJA = 1,IFIX4
      LBHold(JJA) = LBH(JJA)
      SUM4old = SUM4old + CIHsuma(JJA)
      ENDDO
 4082 A=0



      NOORDER4 = 0 
c     -----------------------------
      IF (NOORDER4.EQ.1) GO TO 4088
c     -----------------------------
c     4-excitations:
c     -----------------------------------
      call SORT2(ISEL4,CIH,CIHsuma,LBH)
c     -----------------------------------


      DO JJJ=1,ISEL4
c     -----------------
      J = ISEL4 + 1 - JJJ
      IF (JJJ.GE.J) GO TO 5089
      A1 = CIH(JJJ)
      A2 = CIH(J)
      B1 = CIHsuma(JJJ)
      B2 = CIHsuma(J)
      IC1 = LBH(JJJ)
      IC2 = LBH(J) 

      CIH(JJJ) = A2
      CIH(J) = A1
      CIHsuma(JJJ) = B2
      CIHsuma(J) = B1
      LBH(JJJ) = IC2 
      LBH(J) = IC1 
      ENDDO
 5089 A=0

      CI4newTOT = 0.0D0 
      DO J=1,ISEL4
      CI4newTOT = CI4newTOT + CIHsuma(J)
      ENDDO

      write (6,*) 'NEW ordered set for 4-tuples:' 



      CI4def = CI4newTOT
      CI4norm = 0.0D0
      CI4norm2 = 0.0D0
      ICUT4 = 0
      DO J=1,ISEL4
c     -----------------
      ICUT4 = ICUT4 + 1
      CI4norm2 = CI4norm2 + CIHsuma(J)
      RATIO4 = CI4norm2
      RATIO4 = RATIO4 / CI4newTOT 
      QCUT = CIHsuma(J)
      IF (RATIO4.GE.TNORM4) GO TO 5098
      ENDDO
      write (6,*) 'Will stop. SP4 norm4 has problem'
      write (6,*) 'Triangle is too small !!! '
      STOP
 5098 A=0


c     --------------
      NKEEP4 = ICUT4
c     --------------
      write (6,*) '                           ----------'
      write (6,*) '                           CRITERION-1'
      write (6,*) '                           ----------'
      write (6,*) '                           ICYC4 =',ICYC4 
      write (6,*) 'out of selected            ISEL4 =',ISEL4
      write (6,*) 'will be kept only         NKEEP4 =',NKEEP4
      write (6,*) '                           ICUT4 =',ICUT4
      write (6,*) 'Norm-4 ratio   is         RATIO4 =',RATIO4
      write (6,*) 'Norm-4 target  is         TNORM4 =',TNORM4


      CI4def = CI4newTOT

      DO J=1,NKEEP4
c     ----------------------------------
      CI4def = CI4def - CIHsuma(J)
      CI4norm = CI4norm + CIHsuma(J)
      ENDDO
      RELDEF4 = CI4def / CI4newTOT
      REL4norm = CI4norm / CI4newTOT
      RAT4oldnew = CI4oldTOT / CI4newTOT
c     ----------------------------------
      write (6,*) 'Total 4-norm             CI4norm =',CI4norm
      write (6,*) 'Total 4-deficiency      CI4def   =',CI4def
      write (6,*) 'Relative 4-deficency    RELDEF4  =',RELDEF4
      write (6,*) '                           ----------'
      write (6,*) '                           CRITERION-2'
      write (6,*) '                           ----------'
      write (6,*) '                           ICYC4 =',ICYC4 
      write (6,*) 'Total 4-norm (new)     CI4newTOT =',CI4newTOT
      write (6,*) 'Total 4-norm (old)     CI4oldTOT =',CI4oldTOT
      write (6,*) 'The ratio of these    RAT4oldnew =',RAT4oldnew
      write (6,*) 'The target is          OKNORM4   =',OKNORM4



c     -------------------
c     CRITERION-2 is used:
c     -------------------
      IF (RAT4oldnew.GE.OKNORM4) GO TO 4087 
      CI4oldTOT = CI4newTOT

c     -------------------------------
 491  CONTINUE
c     DO LOOP FOR QUADRUPLES CONTINUE
c     -------------------------------

c     Enlarging triangles (D*D) for 4-uple selection
c     -------------
      write (6,*) 'Q-uples: will stop now'
      write (6,*) 'Convergence of normalization4 not achieved'
      write (6,*) 'actual ratio old/new =',RAT4oldnew
      write (6,*) 'target is            =',OKNORM4
      STOP



 4087 A=0
      write (6,*) 'Q-uples: CRITERION2 is satisfied'
      write (6,*) 'Q-uples: (i)/(i+1) ratio is greater than OKNORM4'

      TNORMAL = TNORMAL + CI4norm
      ALLNORMAL = ALLNORMAL + CI4newTOT

      write (6,*) '-------------------------------------'
      write (6,*) 'For SP4 will keep    NKEEP4 =',NKEEP4
      write (6,*) 'out of total          ISEL4 =',ISEL4
      write (6,*) '4-uples: for NKEEP4 TNORMAL =',TNORMAL
      write (6,*) 'at 4-uples: ALLNORMAL(4,5,6)=',ALLNORMAL
      write (6,*) '-------------------------------------'







 4088 A=0
      write (6,*) 'SELECTION OF QUADRUPLE LIVE-WOOD BEGINS'
      write (6,*) 'WILL INCLUDE MOST IMPORTANT SP4=',NKEEP4



      ICYCMIN = 2
      ICYCMAX = 2



      DO 500 ICYC=ICYCMIN,ICYCMAX
c     ----------------------------
c     ***


C     -----------------------
C     QUADRUPLE-EXCITATIONS:
C     -----------------------
       write (6,*) 'Q-selection is step =',ICYC  
       write (6,*) 'At SDT-CI level there are:'  
       write (6,*) 'Symmetry-reduced strings =',ICOUNTS  

       ICOUNT4 = 0

       write (6,*) 'Q-case: OCCU=[0,0]'
C      occupied MO space: CASE-A)  [0,0]
c      *********************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO 505 I=1,NOCC-1
c      ------------------
       M = NOCC+1-I
       NEL(M) = 0
       DO 510 J=1,M-1
c      ------------------
       MM = M - J
       NEL(MM) = 0


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      CASE-1) [2,2]
c      *************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO 525 KK=1,NVIRT-1
       MK = NTOT+1-KK
       NEL(MK) = 2

       DO 530 KKK=1,NVIRT-KK
       MMK = MK - KKK
       NEL(MMK) = 2

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
c      * 
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 535
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 535 
      ICOUNTS = ICOUNTS + 1
      ISTRA(ICOUNTS) = IPOSA(ICNT)
      ISTRB(ICOUNTS) = IPOSB(ICNT)
 535  AAA=0
       NEL(MMK) = 0
 530   CONTINUE
       NEL(MK) = 0
 525   CONTINUE


c      Virtual MOs
c      CASE-2 [2,1,1]
c      ------------------      
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO 545 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2
       DO 550 KK1=1,NVIRT-2
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 551 
       IF (M1.LE.NX) GO TO 552
 551   AAA=0
       L1 = M1
       GO TO 553
 552   AAA=0
       L1 = M1 - 1
 553   AAA=0
       NEL(L1) = 1
       DO 560 KK2=1,M1-NOCC-2    
       M2 = M1 - KK2
       IF (M2.GT.NX) GO TO 561 
       IF (M2.LE.NX) GO TO 562 
 561   AAA=0
       L2 = M2
       GO TO 563
 562   AAA=0
       L2 = M2 - 1
 563   AAA=0
       NEL(L2) = 1

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 568
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
c      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c      write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c      write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 568 
      ICOUNTS = ICOUNTS + 1
      ISTRA(ICOUNTS) = IPOSA(ICNT)
      ISTRB(ICOUNTS) = IPOSB(ICNT)
 568  AAA=0
       NEL(L2) = 0
 560   CONTINUE
       NEL(L1) = 0
 550   CONTINUE
       NEL(NX) = 0
 545   CONTINUE

c      virtual MOs
c      CASE-3) [1,1,1,1]
c      *
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO 570 KK = 1,NVIRT-3
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 575 KK1 = 1,M1-NOCC-3
       M2 = M1 - KK1
       NEL(M2) = 1

       DO 580 KK2 = 1,M2-NOCC-2
       M3 = M2 - KK2
       NEL(M3) = 1

       DO 585 KK3 = 1,M3-NOCC-1
       M4 = M3 - KK3
       NEL(M4) = 1

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 590
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 

       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 590 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 590   AAA=0
       NEL(M4) = 0
 585   CONTINUE
       NEL(M3) = 0
 580   CONTINUE
       NEL(M2) = 0
 575   CONTINUE
       NEL(M1) = 0
 570   CONTINUE
c      ------------------------------------
c      end of occupied MOs: CASE-A) [0,0]
c      ------------------------------------
       NEL(MM) = 2
 510   CONTINUE
       NEL(M) = 2
 505   CONTINUE
 595   AAA=0.0


       write (6,*) 'Q-case: OCCU=[0,1,1]'
c      occupied MOs: CASE-B) [0,1,1]
c      *
       DO 605 I= 1,NOCC
       NEL(I) = 2
 605   CONTINUE
       DO 620 KKZ=1,NOCC
       NXZ = NOCC + 1 - KKZ
       NEL(NXZ) = 0
       DO 627 KK1Z=1,NOCC-2
       M1Z = NOCC + 1 - KK1Z 
       IF (M1Z.GT.NXZ) GO TO 631 
       IF (M1Z.LE.NXZ) GO TO 632 
 631   AAA=0
       L1Z = M1Z
       GO TO 633 
 632   AAA=0
       L1Z = M1Z - 1
 633   AAA=0
       NEL(L1Z) = 1
       DO 640 KK2Z=1,M1Z-2    
       M2Z = M1Z - KK2Z
       IF (M2Z.GT.NXZ) GO TO 641 
       IF (M2Z.LE.NXZ) GO TO 642 
 641   AAA=0
       L2Z = M2Z
       GO TO 643
 642   AAA=0
       L2Z = M2Z - 1
 643   AAA=0
       NEL(L2Z) = 1


c      virtual space
c      CASE-1) [2,2]
c      --------------
       DO 650 K = NOCC+1,NTOT
       NEL(K) = 0
 650   CONTINUE
       DO 655 KK=1,NVIRT-1
       MK = NTOT+1-KK
       NEL(MK) = 2
       DO 660 KKK=1,NVIRT-KK
       MMK = MK - KKK
       NEL(MMK) = 2

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 665 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 665 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 665   AAA=0
       NEL(MMK) = 0
 660   CONTINUE
       NEL(MK) = 0
 655   CONTINUE


c      Virtual MOs
c      CASE-2 [2,1,1]
c      *
       DO 670 K = NOCC+1,NTOT
       NEL(K) = 0
 670   CONTINUE
       DO 675 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2
       DO 680 KK1=1,NVIRT-2
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 681 
       IF (M1.LE.NX) GO TO 682
 681   AAA=0
       L1 = M1
       GO TO 683
 682   AAA=0
       L1 = M1 - 1
 683   AAA=0
       NEL(L1) = 1



       DO 690 KK2=1,M1-NOCC-2    
c      ------------------------
       M2 = M1 - KK2
       IF (M2.GT.NX) GO TO 691 
       IF (M2.LE.NX) GO TO 692 
 691   AAA=0
       L2 = M2
       GO TO 693
 692   AAA=0
       L2 = M2 - 1
 693   AAA=0
       NEL(L2) = 1

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
c      * 
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 698 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
c      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c      write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c      write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 698 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 698   AAA=0
       NEL(L2) = 0
 690   CONTINUE
       NEL(L1) = 0
 680   CONTINUE
       NEL(NX) = 0
 675   CONTINUE

c      virtual MOs
c      CASE-3) [1,1,1,1]
c      ------------------ 
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO 700 KK = 1,NVIRT-3
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 705 KK1 = 1,M1-NOCC-3
       M2 = M1 - KK1
       NEL(M2) = 1
       DO 710 KK2 = 1,M2-NOCC-2
       M3 = M2 - KK2
       NEL(M3) = 1
       DO 715 KK3 = 1,M3-NOCC-1
       M4 = M3 - KK3
       NEL(M4) = 1

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 719 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
c      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 719 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 719   AAA=0
       NEL(M4) = 0
 715   CONTINUE
       NEL(M3) = 0
 710   CONTINUE
       NEL(M2) = 0
 705   CONTINUE
       NEL(M1) = 0
 700   CONTINUE
c      -----------------------------
c      end of occupied MOs: CASE-B)
c      -----------------------------
       NEL(L2Z) = 2
 640   CONTINUE
       NEL(L1Z) = 2
 627   CONTINUE
       NEL(NXZ) = 2
 620   CONTINUE
 750   AAA=0


       write (6,*) 'Q-case: OCCU=[1,1,1,1]'
c      C. occupied MOs CASE-C: [1,1,1,1]
c      ------------------------------------
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO 805 KKW = 1,NOCC-3
       M1W = NOCC+1-KKW
       NEL(M1W) = 1
       DO 810 KK1W = 1,M1W-3
       M2W = M1W - KK1W
       NEL(M2W) = 1
       DO 815 KK2W = 1,M2W-2
       M3W = M2W - KK2W
       NEL(M3W) = 1
       DO 820 KK3W = 1,M3W-1
       M4W = M3W - KK3W
       NEL(M4W) = 1

c      virtual space
c      Q: CASE-1) [2,2]
c      -----------------------
       DO 850 K = NOCC+1,NTOT
       NEL(K) = 0
 850   CONTINUE
       DO 855 KK=1,NVIRT-1
       MK = NTOT+1-KK
       NEL(MK) = 2
       DO 860 KKK=1,NVIRT-KK
       MMK = MK - KKK
       NEL(MMK) = 2
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 865 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 865 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 865   AAA=0
       NEL(MMK) = 0
 860   CONTINUE
       NEL(MK) = 0
 855   CONTINUE

c      Virtual MOs
c      Q: CASE-2 [2,1,1]
c      ----------------------- 
       DO 870 K = NOCC+1,NTOT
       NEL(K) = 0
 870   CONTINUE
       DO 875 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2
       DO 880 KK1=1,NVIRT-2
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 881 
       IF (M1.LE.NX) GO TO 882
 881   AAA=0
       L1 = M1
       GO TO 883
 882   AAA=0
       L1 = M1 - 1
 883   AAA=0
       NEL(L1) = 1
       DO 890 KK2=1,M1-NOCC-2    
       M2 = M1 - KK2
       IF (M2.GT.NX) GO TO 891 
       IF (M2.LE.NX) GO TO 892 
 891   AAA=0
       L2 = M2
       GO TO 893
 892   AAA=0
       L2 = M2 - 1
 893   AAA=0
       NEL(L2) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 898 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c       write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c       write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 898 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 898   AAA=0
       NEL(L2) = 0
 890   CONTINUE
       NEL(L1) = 0
 880   CONTINUE
       NEL(NX) = 0
 875   CONTINUE

c      virtual MOs
c      Q: CASE-3) [1,1,1,1]
c      -------------------- 
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO 900 KK = 1,NVIRT-3
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 905 KK1 = 1,M1-NOCC-3
       M2 = M1 - KK1
       NEL(M2) = 1
       DO 910 KK2 = 1,M2-NOCC-2
       M3 = M2 - KK2
       NEL(M3) = 1
       DO 915 KK3 = 1,M3-NOCC-1
       M4 = M3 - KK3
       NEL(M4) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT4 = ICOUNT4 + 1
       CALL FILTER(NKEEP4,LBH,ICOUNT4,IPUT)
       IF (IPUT.EQ.0) GO TO 919 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
C*    -----------------------------------
c      CALL WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
C*    -----------------------------------
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
c      WRITE (6,*) 'position of alpha string',IPOSA(ICOUNT) 
c      WRITE (6,*) 'position of beta  string',IPOSB(ICOUNT) 
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
c      write (6,*) 'alfa string symmetry(irrep)',ISYMA(ICOUNT)
c      write (6,*) 'beta string symmetry(irrep)',ISYMB(ICOUNT)
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 919 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 919   AAA=0
       NEL(M4) = 0
 915   CONTINUE
       NEL(M3) = 0
 910   CONTINUE
       NEL(M2) = 0
 905   CONTINUE
       NEL(M1) = 0
 900   CONTINUE
c      ---------------------------------------
c      end of occupied MOs CASE-C: Q[1,1,1,1]
c      ---------------------------------------
       NEL(M4W) = 2 
 820   CONTINUE
       NEL(M3W) = 2 
 815   CONTINUE
       NEL(M2W) = 2 
 810   CONTINUE
       NEL(M1W) = 2 
 805   CONTINUE
 930   A=0
c      -----------------
c*     ICYC-loop
c      -----------------
c      Q-excitations
c      -----------------
 500   CONTINUE 
c      ---------------
       NQUADRO = ISEL4
c      ---------------
       write (6,*) '-------------------------------------'
       write (6,*) 'THE QUADRUPLE SP LIVE-WOOD SELECTION COMPLETED'
       write (6,*) 'CI-SDTQ only'
       write (6,*) '# Q-sp NQUADRO=',NQUADRO
       write (6,*) '# Q-sp ISEL4=',ISEL4
       write (6,*) 'Selection is complete for NKEEP4=',NKEEP4
       write (6,*) 'Symmetry reduced strings ICOUNTS =',ICOUNTS
       write (6,*) '-------------------------------------'
       write (6,*) '*'





c      QUINTUPLE (x=5) EXCITATIONS START HERE:
c      --------------------------------------
       IF (IQUINT.EQ.0) GO TO 1000
c      --------------------------------------
       write (6,*) '*'
       write (6,*) '5-tuple: quintuple excitations start here'
       write (6,*) '*'
c      ----------------------------------------------
c      (1) estimate 5-tuple SP up to (ISEL5 = NKEEP5)
c      ----------------------------------------------
c      (2) make the configurational selection;
c      ----------------------------------------------
c      CI5(j)
c      CI5suma(j)
c      LB5(j),         j=1,ISEL5
c      ----------------------------------------------
c     ESTIMATE: CI5sumTOT
c     -------------------
      CI5sumTOT = 0.0D0

c**** CI5sumTOT = (CI2sumTOT * CI3sumTOT)
c**** write (6,*) 'ESTIMATED 5-tuple NORM. DEFICIENCY:'
c**** write (6,*) 'CI5sumTOT (sum2 * sum3) =',CI5sumTOT



      Dsum = 0.0D0
      Tsum = 0.0D0
      ITcount = 0
      JJD = 2 
      PR1 = CI2suma(JJD)*CI3suma(1)
      PR2 = CI3suma(JJD)*CI2suma(1)

c----------LAIMIS-2007 selection of istep3-------
c      write (6,*) '---'
c      write (6,*) 'J =',JJD
c      write (6,*) 'K =',JJD
c      write (6,*) 'PR1(D(J)*T(1)) =',PR1
c      write (6,*) 'PR2(D(1)*T(K)) =',PR2
c      write (6,*) '---'

      DO IT = JJD,ISEL3
c     ----------------
c     ***
      ITcount = ITcount + 1 
      ITfinal = IT

      PR2 = CI3suma(IT)*CI2suma(1)

      IF (PR2.GT.PR1) GO TO 1011
      GO TO 1012

 1011 A=0
      PR2old = PR2
      Kold = IT
c     ---------------
      ENDDO


c----------LAIMIS-2007 selection of istep3-------


 1012 A=0










c     ----------------------- 
      MAXSP5 = 3 * NKEEP5new 
      SMALL5 = 0.0000000000002D0
c     ----------------------- 
      write (6,*) '---'
      write (6,*) 'maximum comp. SP5 MAXSP5 =',MAXSP5
      write (6,*) '5-tuple SMALL5 =',SMALL5
      write (6,*) '-------------------------------------------------'
      write (6,*) 'CONVERGENCE CRITERIUM: internal   TNORM5 =',TNORM5
      write (6,*) 'CONVERGENCE CRITERIUM: external  OKNORM5 =',OKNORM5
      write (6,*) '-------------------------------------------------'
c--------------------------------------------------
c     OPEN-ENDED APPROACH IMPLEMENTED FOR
c     5-TUPLES WITH ENLARGING THE TRIANGLE
c     TO IMPROVE THE SELECTION PROCEDURE FOR SP5
c
c     HERE: MDEL5(INCREMENT STEPS)
c     ICYCMAX5-1 = No of INCREMENT CYCLES
c-----------------------------------------
c      THE CONVENTION FOR ANGLE DEFINITION:
c      -----------------------------------


      ISTEP3 = ITfinal
      ISTEP3 = 6 
      STEP3 = ISTEP3

c     -----------
      write (6,*) '---'
      write (6,*) '5-tuple selection: lines-angle'
      write (6,*) 'THE ANGLE for (D*T) triangle'
      write (6,*) 'SCREENING STARTS FROM: D(2)*T(1)'
      write (6,*) 'TO: D(1)*T(ISTEP3+1)'
      write (6,*) 'NOTE FOR 5-tuples ISTEP3 is fixed'
      write (6,*) 'HERE IS USED: ISTEP3(FIXED) =',ISTEP3
      write (6,*) ' STEP3 =',STEP3
      write (6,*) ' SUGGESTED ANGLE (ITfinal) =',ITfinal
      write (6,*) '---'
c     -----------------

      CI5oldTOT = 0.0D0 

c     ------------
      IFIX5 = 1
c     ------------
      DO JJJ=1,LIMITAS5
c     -----------------
      CIH(JJJ) = 0.0D0
      CIHsuma(JJJ) = 0.0D0
      LBH(JJJ) = 0
      LBHold(JJJ) = 0
      ENDDO
      DO 1055 ICYC5 = 1,ICYCMAX5
c     --------------------------
      NKEEP5 = NKEEP5new
      INUMB5 = 0
      ISEL5 = 0
      ITOT5 = 0
      IF (ICYC5.EQ.1) GO TO 1056
      MAX5IM2 = MAX5IM2 + MDEL5
c     ------------------------- 
 1056 A=0
      DO JJJ=1,LIMITAS5
c     -----------------
      CIH(JJJ) = 0.0D0
      CIHsuma(JJJ) = 0.0D0
      LBH(JJJ) = 0
      ENDDO

      write (6,*) '*'
      write (6,*) '*'
      write (6,*) '*'
      write (6,*) '                   ********************************'
      write (6,*) '                           ICYC5 =',ICYC5
      write (6,*) '                   ********************************'
      write (6,*) '5-tuple estimation starts at doubles(SP) =',MAX5IM2 


      IF (MAX5IM2.LE.ISEL2) GO TO 1057
c     ***
      write (6,*) '---'
      write (6,*) '                     WARNING'
      write (6,*) 'x=5: MAX5IM2 Exceeded doubles ISEL2=',ISEL2
      write (6,*) '---'
c     ***
 1057 A=0
c     ***


      EMAX5IM2 = MAX5IM2 + 1
      EMAX5IM2 = EMAX5IM2 - (1.0D0/STEP3)
c*****write (6,*) '                   Exact D-SP: EMAX5IM2 =',EMAX5IM2


      JD1 = 1 
c     *
c     start with doubles
      DO 1013 JD = JD1,MAX5IM2
c     ------------------------
c*    *** number of lines ***
c     -----------------------

      DO 1014 LIN = 1,ISTEP3
c     -------------------------
c     *
      KDD = JD + 1
c     ------------------------------------------------ 
c                        IMPORTANT !!!
c     For SP5 selection procedure only the
c     all SP3 configurations with ISEL3 will be used
c     Although the NKEEP3 might be less than ISEL3
c     ------------------------------------------------ 
c     *** on each line units ***
c     --------------------------
      DO 1015 JT = 1,JD 
c     -------------------------
c     *
      KTT = ((JT-1)*ISTEP3) + 1
      KTT = KTT + LIN - 1
c     -----------------------------
      IF (KTT.LE.ISEL3) GO TO 1017
c     -----------------------------
      GO TO 1044 
 1017 A=0
      KDD = KDD - 1
      IF (KDD.GT.ISEL2) GO TO 1040 
      IF (KDD.GE.1) GO TO 1018
      write (6,*) 'will stop KDD =',KDD 
      STOP
 1018 A=0
      KD = KDD 
      KT = KTT 
      I2TRIA = KD 
      I3TRIA = KT 
      IDTRIANG = JD
      ILIN = LIN
      ITOT5 = ITOT5 + 1 
      LAB2 = LB2(KD)
      LAB3 = LB3(KT)

      CALL RLABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,
     *ICON2,LAB2)
      CALL RLABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT3,
     *ICON3,LAB3)
      CALL NUM555(NTOT,NOCC,NVIRT,ICON2,ICON3,
     *ICON5,IPUT)

      IF (IPUT.EQ.0) GO TO 1040 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT5,ICON5,
     *IDNUM5)
      INUMB5 = INUMB5 + 1
      IF (ISEL5.LE.NKEEP5new) GO TO 1052
      write (6,*) '---'
      write (6,*) 'WARNING(x=5)'
      write (6,*) 'distinct SP5           ISEL5=',ISEL5
      write (6,*) 'exceeded max SP5  NKEEP5new =',NKEEP5new
      write (6,*) '    MUST INCREASE DIM = NKEEP5'
      write (6,*) 'x=5: cycle ICYC5 now ends'
      write (6,*) '---'
      GO TO 1050
 1052 A=0
      IF (INUMB5.LE.MAXSP5) GO TO 1051
      write (6,*) '---'
      write (6,*) 'WARNING(x=5)'
      write (6,*) '5-t: there are components INUMB5 =',INUMB5
      write (6,*) 'exceeded max allowed      MAXSP5 =',MAXSP5
      write (6,*) '         # components = 3 * NKEEP5'
      write (6,*) '    MUST INCREASE DIM = NKEEP5'
      write (6,*) 'x=5: cycle ICYC5 now ends'
      write (6,*) '---'
      GO TO 1050
 1051 A=0

      PD5 = (CI2EST(KD) * CI3EST(KT))
      PD5suma = (CI2suma(KD) * CI3suma(KT))
      IF (INUMB5.EQ.1) GO TO 1035

c     --------------------------------
c     check if IDNUM5 has occured before:
c     -----------------------------------
      DO ITT=1,ISEL5
      IHIT = ITT
      IF (IDNUM5.NE.LBH(ITT)) GO TO 1030 
      CIH(IHIT) = CIH(IHIT) + PD5 
      CIHsuma(IHIT) = CIHsuma(IHIT) + PD5suma 
      GO TO 1040
 1030 A=0
      ENDDO
 1035 A=0
c     -----------------
c     (5-tuple is new):
c     -----------------
      ISEL5 = ISEL5 + 1
      LBH(ISEL5) = IDNUM5
      CIH(ISEL5) = CIH(ISEL5) + PD5 
      CIHsuma(ISEL5) = CIHsuma(ISEL5) + PD5suma 

c-------------------------------------------
c*****LAIMIS:ISEL5<mxsprod
c-------------------------------------------
      IF (ISEL5.LT.mxsprod) GO TO 1040
      write (6,*) 'PROBLEM: ISEL5 is ',ISEL5
      write (6,*) '       mxsprod is ',mxsprod
      write (6,*) 'NOW PROGRAM WILL STOP '
      write (6,*) 'MUST INCREASE DIM = mxsprod'
      STOP
 1040 A=0
      IF (KTT.GT.ISEL3) GO TO 1044

 1015 CONTINUE
c     *** (UNITS)
 1044 A=0
c     -----------------------------
      IF (KDD.LT.1) GO TO 1050
c     -----------------------------

 1014 CONTINUE
c     --------
c     *** (LINES)

 1013 CONTINUE
c     --------
c     (doubles)

 1050 A=0

c**** write (6,*) '---'
c**** write (6,*) '           ICYC5 =',ICYC5
c**** write (6,*) '5-TUPLE CYCLE INFO:'
c**** write (6,*) 'Final values: KD =',KD
c**** write (6,*) 'Final values: KT =',KT
c**** write (6,*) '---'
c****     5-tuple selection completed

      write (6,*) '--------------------------------------------------'
      write (6,*) '                               ICYC5 =',ICYC5
      write (6,*) '         Total (D*T) products  ITOT5 =',ITOT5
      write (6,*) '               SP5 components INUMB5 =',INUMB5
      write (6,*) '              Distinct SP5 are ISEL5 =',ISEL5
      write (6,*) '                              NKEEP5 =',NKEEP5
      write (6,*) '--------------------------------------------------'

c      write (6,*) 'Started triangle at doubles=',IDTRIANG 
c      write (6,*) 'JA,LB5(JA),CI5(JA),CI5suma(JA)'
c      write (6,*) 'JA,LBH(JA),CIH(JA),CIHsuma(JA)'
c      write (6,*) '---------------------------'
c      DO JA=1,10
c     ------------
c*****write (6,1078) JA,LB5(JA),CI5(JA),CI5suma(JA)
c      write (6,1078) JA,LBH(JA),CIH(JA),CIHsuma(JA)
c 1078 FORMAT (I7,1X,I7,3X,E12.6,2X,E12.6)
c      ENDDO
c-----------------------------------------------------------------


c     *  IF ISEL5 < NKEEP5  *
      IF (ISEL5.GE.NKEEP5) GO TO 1079
c     -------------------------------
      NKEEP5 = ISEL5
      write (6,*) '                        NOTE: ISEL5 < NKEEP5'
      write (6,*) '                             NKEEP5 =',NKEEP5
 1079 A=0


c     ****** ORDER 5-tuples now: *****
c     ---------------------------------
      IF (ISEL5.GE.IFIX5) GO TO 1080
      write (6,*) 'Number of distinct SP5 ISEL5 is too small'
      STOP
 1080 A=0
      IF (ICYC5.GT.1) GO TO 1081
      IFIX5 = ISEL5
      SUM5old = 0.0D0
      DO JJA = 1,IFIX5
      LBHold(JJA) = LBH(JJA)
      SUM5old = SUM5old + CIHsuma(JJA)
      ENDDO
 1081 A=0.0
      NOORDER5 = 0 
      IF (NOORDER5.EQ.1) GO TO 1088
c     -----------------------------------
      call SORT2(ISEL5,CIH,CIHsuma,LBH)
c     -----------------------------------
      DO JJJ=1,ISEL5
      J = ISEL5 + 1 - JJJ
      IF (JJJ.GE.J) GO TO 5077
      A1 = CIH(JJJ)
      A2 = CIH(J)
      B1 = CIHsuma(JJJ)
      B2 = CIHsuma(J)
      IC1 = LBH(JJJ)
      IC2 = LBH(J) 
      CIH(JJJ) = A2
      CIH(J) = A1
      CIHsuma(JJJ) = B2
      CIHsuma(J) = B1
      LBH(JJJ) = IC2 
      LBH(J) = IC1 
      ENDDO
 5077 A=0
      CI5newTOT = 0.0D0 
      DO J=1,ISEL5
      CI5newTOT = CI5newTOT + CIHsuma(J)
      ENDDO




      CI5def = CI5newTOT
      CI5norm = 0.0D0
      CI5norm2 = 0.0D0
      ICUT5 = 0
c     -------------------
c     * CRITERION-1 *
c     -------------------
      ICUT5 = 0
      DO J=1,ISEL5
      ICUT5 = ICUT5 + 1
      CI5norm2 = CI5norm2 + CIHsuma(J)
      RATIO5 = CI5norm2 / CI5newTOT 
      CI5last = CIHsuma(J)
      IF (RATIO5.GE.TNORM5) GO TO 2098
      ENDDO
      write (6,*) 'Will stop. SP5 norm5 has problem'
      write (6,*) 'Triangle is too small !!! '
      STOP
 2098 A=0
      NKEEP5 = ICUT5

      CI5def = CI5newTOT

      DO J=1,NKEEP5
      CI5def = CI5def - CIHsuma(J)
      CI5norm = CI5norm + CIHsuma(J)
      ENDDO

      RELDEF5 = CI5def / CI5newTOT
      REL5norm = CI5norm / CI5newTOT
      write (6,*) '                           ------------'
      write (6,*) '                           CRITERION-1'
      write (6,*) '                           ------------'
      write (6,*) '                           ICYC5 =',ICYC5 
      write (6,*) '                   out of  ISEL5 =',ISEL5
      write (6,*) '                           ICUT5 =',ICUT5
      write (6,*) '             will keep    NKEEP5 =',NKEEP5
      write (6,*) '     5-norm for NKEEP5   CI5norm =',CI5norm
      write (6,*) '       5-relative norm    RATIO5 =',RATIO5
      write (6,*) '       5-relative target  TNORM5 =',TNORM5
      write (6,*) '         5-deficiency     CI5def =',CI5def
      write (6,*) 'RELATIVE 5-deficiency  RELDEF5 =',RELDEF5
c     ----------------------------------------------------------


c----------------------------------------------------------
c     Following option
c     takes all SP5:
c     --------------
c     -------------------
c     * CRITERION-2 *
c     -------------------

      RAT5oldnew = CI5oldTOT / CI5newTOT

      write (6,*) '                           ------------'
      write (6,*) '                           CRITERION-2'
      write (6,*) '                           ------------'
      write (6,*) '                           ICYC5 =',ICYC5 
      write (6,*) 'total 5-norm (new)     CI5newTOT =',CI5newTOT
      write (6,*) 'total 5-norm (old)     CI5oldTOT =',CI5oldTOT
      write (6,*) 'The ratio of old/new             =',RAT5oldnew
      write (6,*) 'The target is                    =',OKNORM5
      IF (RAT5oldnew.GE.OKNORM5) GO TO 1088 
c     -------------------------------------
c     IF THE CONVERGENCE IS SATISFIED GO TO
c     CRITERION-1 TO SELECT THE SPECIFIED 
c     NORMALIZATION PERCENTAGE
c     -------------------------------------
      CI5oldTOT = CI5newTOT

 1055 CONTINUE
c     DO LOOP ENDS HERE
c
c     Enlarged D*T triangle for 5-tuples
      write (6,*) 'Will stop now: 5-tuple selection incomplete'
      write (6,*) 'Neeed more cycles(5) to succeed'
      STOP
 1088  A=0



      TNORMAL = TNORMAL + CI5norm
      ALLNORMAL = ALLNORMAL + CI5newTOT

      write (6,*) 'CRITERION-1: keep only THESE SP5 =',NKEEP5
      write (6,*) 'CRITERION-2: treshold between last triangles OK' 
      write (6,*) '5-: for NKEEP5    TNORMAL(4,5,6) =',TNORMAL
      write (6,*) 'at 5: for ISEL5  ALLNORMAL(4,5,6)=',ALLNORMAL
      write (6,*) '    '
      write (6,*) '---------------------------------------------- '
      write (6,*) 'FINAL STEP for 5-TUPLES: Start the selection '
      write (6,*) '   INCLUDING ONLY THESE SP5 =',NKEEP5




       ICYCMIN = 2
       ICYCMAX = 2



       DO 1099 ICYC=ICYCMIN,ICYCMAX
c      ----------------------------
c      ***



c      ----------------------
c      5-excitations: 
c      ----------------------
       write (6,*) '5-configuration selection is STEP = ',ICYC
       write (6,*) '      At SDTQ-CI level there are:'
       write (6,*) 'Symmetry-reduced strings ICOUNTS =',ICOUNTS

       ICOUNT5 = 0

c      SCF-SPACE:
c      ----------
       write (6,*) '--------------------------'
       write (6,*) 'Quintuple: OCCU=[1,0,0] starts '
C      occupied MO space: CASE  [1,0,0]
c      ---------------------------------------- 
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       MAX1 = NOCC-1 
       DO 1500 I1=1,MAX1
       N1 = NOCC + 1 - I1 
       NEL(N1) = 0
       MAX2 = N1-1
       DO 1505 I2=1,MAX2
       N2 = N1 - I2
       NEL(N2) = 0
       MAX3 = NOCC-2
       DO 1510 I3=1,MAX3
       NN3 = NOCC + 1 - I3
       IF (I3.EQ.1) GO TO 1511
       NN3 = LL3 - 1
 1511  A=0
       call SHIFT2(N1,N2,NN3,LL3)
       IF (LL3.EQ.0) GO TO 1512 
       NEL(LL3) = 1

c      virtual space
c      5: 'VIRT [2,2,1] starts'
c      CASE) VIRT [2,2,1]
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-1
       DO K1=1,LIM1
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 
       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
       NX2 = NX1 - K2
       NEL(NX2) = 2 
       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1535
       M1 = L1 - 1
 1535  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1538 
       NEL(L1) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1540
c      -------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1540 
      ICOUNTS = ICOUNTS + 1
      ISTRA(ICOUNTS) = IPOSA(ICNT)
      ISTRB(ICOUNTS) = IPOSB(ICNT)
 1540 AAA=0
       NEL(L1) = 0
       ENDDO
 1538  A=0
       NEL(NX2) = 0
       ENDDO
       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'


c      virtual space
c      CASE) 5: VIRT [2,1,1,1]
c      ----------------------- 
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO K=1,NVIRT
       NX = NTOT+1-K
       NEL(NX) = 2
       LIM1 = NVIRT-3
       DO KK1=1,LIM1
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1545
       M1 = L1 - 1
 1545  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1
       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1546
       M2 = L2 - 1
 1546  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1
       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1547
       M3 = L3 - 1
 1547  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1550
       NEL(L3) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1551
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1551 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 1551  AAA=0
       NEL(L3)=0
       ENDDO
 1550  A=0
       NEL(L2)=0
       ENDDO
       NEL(L1)=0
       ENDDO
       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'

c      virtual space
c      CASE) VIRT [1,1,1,1,1]
c      ---------------------- 
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-4
       DO KK1=1,LIM1
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1
       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
       M2 = M1 - KK2
       NEL(M2) = 1
       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
       M3 = M2 - KK3
       NEL(M3) = 1
       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
       M4 = M3 - KK4
       NEL(M4) = 1
       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1555
       NEL(M5) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1556
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1556 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 1556  AAA=0
       NEL(M5)=0
       ENDDO
 1555  A=0
       NEL(M4)=0
       ENDDO
       NEL(M3)=0
       ENDDO
       NEL(M2)=0
       ENDDO
       NEL(M1)=0
       ENDDO
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
c      -----------------------------
c      end of occupied MOs: CASE [1,0,0]
c      -----------------------------
       NEL(LL3) = 2 
 1510  CONTINUE
 1512  A=0
       NEL(N2) = 2 
 1505  CONTINUE
       NEL(N1) = 2 
 1500  CONTINUE

       write (6,*) ' Quintuple OCCU=[1,0,0] done'
       write (6,*) '      5-t: total count ICOUNT5 =',ICOUNT5
       write (6,*) '     total string count ICOUNT =',ICOUNT
       write (6,*) ' symmetry reduced string count =',ICOUNTS
       write (6,*) '--------------------------------------------'

 1600  A=0


cLAIMIS
c      SCF-SPACE:
       write (6,*) 'Quintuple OCCU=[1,1,1,0] starts'
C      occupied MO space: CASE  [1,1,1,0]
c      -------------------------------------------- 
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       MAX1 = NOCC
       DO 1605 I1=1,NOCC
       N1 = NOCC + 1 - I1
       NEL(N1) = 0
       MAX2 = NOCC-3
       DO 1610 I2=1,MAX2
       NN2 = NOCC + 1 - I2
       IF (I2.EQ.1) GO TO 1611
       NN2 = LL2 - 1
 1611  A=0
       call SHIFT1(N1,NN2,LL2)
       NEL(LL2) = 1
       MAX3 = NN2 - 2 
       DO 1615 I3=1,MAX3
       NN3 = LL2 - I3
       IF (I3.EQ.1) GO TO 1616
       NN3 = LL3 - 1
 1616  A=0
       call SHIFT1(N1,NN3,LL3)
       NEL(LL3) = 1
       MAX4 = NN3 - 1 
       DO 1620 I4=1,MAX4
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 1621
       NN4 = LL4 - 1
 1621  A=0
       call SHIFT1(N1,NN4,LL4)
       IF (LL4.EQ.0) GO TO 1622 
       NEL(LL4) = 1

c      virtual space
c      CASE 5: VIRT [2,2,1]
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-1
       DO K1=1,LIM1
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 
       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
       NX2 = NX1 - K2
       NEL(NX2) = 2 
       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1635
       M1 = L1 - 1
 1635  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1638 
       NEL(L1) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1640
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1640 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 1640  AAA=0
       NEL(L1) = 0
       ENDDO
 1638  A=0
       NEL(NX2) = 0
       ENDDO
       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'

c      virtual space
c      CASE) 5: VIRT [2,1,1,1]
c      ------------------------ 
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO K=1,NVIRT
       NX = NTOT+1-K
       NEL(NX) = 2
       LIM1 = NVIRT-3
       DO KK1=1,LIM1
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1645
       M1 = L1 - 1
 1645  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1
       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1646
       M2 = L2 - 1
 1646  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1
       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1647
       M3 = L3 - 1
 1647  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1650
       NEL(L3) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1651
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1651 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 1651  AAA=0
       NEL(L3)=0
       ENDDO
 1650  A=0
       NEL(L2)=0
       ENDDO
       NEL(L1)=0
       ENDDO
       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'

c      virtual space
c      CASE) VIRT [1,1,1,1,1]
c      ---------------------- 
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-4
       DO KK1=1,LIM1
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1
       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
       M2 = M1 - KK2
       NEL(M2) = 1
       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
       M3 = M2 - KK3
       NEL(M3) = 1
       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
       M4 = M3 - KK4
       NEL(M4) = 1
       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1655
       NEL(M5) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1656
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1656 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 1656  AAA=0
       NEL(M5)=0
       ENDDO
 1655  A=0
       NEL(M4)=0
       ENDDO
       NEL(M3)=0
       ENDDO
       NEL(M2)=0
       ENDDO
       NEL(M1)=0
       ENDDO
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
       NEL(LL4)=2
 1620  CONTINUE
 1622  A=0
       NEL(LL3)=2
 1615  CONTINUE
       NEL(LL2)=2
 1610  CONTINUE
       NEL(N1)=2
 1605  CONTINUE

       write (6,*) 'Quintuple OCCU=[1,1,1,0] done'
       write (6,*) '        5-t: total count ICOUNT5 =',ICOUNT5
       write (6,*) '              total string count =',ICOUNT
       write (6,*) '   symmetry reduced string count =',ICOUNTS
       write (6,*) '----------------------------------------------'
 1700  A=0


cLAIMIS
       write (6,*) 'Quintuple OCCU=[1,1,1,1,1] starts'
c      --------------------------------------- 
c      SCF-SPACE:
C      occupied MO space: CASE  [1,1,1,1,1]
c      --------------------------------------- 
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       MAX1 = NOCC-4
       DO 1805 I1=1,MAX1
       N1 = NOCC + 1 - I1 
       NEL(N1) = 1
       MAX2 = N1 - 4 
       DO 1810 I2=1,MAX2
       N2 = N1 - I2
       NEL(N2) = 1
       MAX3 = N2 - 3 
       DO 1815 I3=1,MAX3
       N3 = N2 - I3
       NEL(N3) = 1
       MAX4 = N3 - 2 
       DO 1820 I4=1,MAX4
       N4 = N3 - I4
       NEL(N4) = 1
       MAX5 = N4 - 1 
       DO 1825 I5=1,MAX5
       N5 = N4 - I5
       IF (N5.GT.0) GO TO 1831 
       write (6,*) 'WILL STOP because N5=',N5
       STOP
 1831  A=0
       NEL(N5) = 1

c      virtual space
c      CASE) VIRT [2,2,1]
c      ------------------ 
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-1
       DO K1=1,LIM1
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 
       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
       NX2 = NX1 - K2
       NEL(NX2) = 2 
       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1835
       M1 = L1 - 1
 1835  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1838 
       NEL(L1) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1840
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1840 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 1840  AAA=0
       NEL(L1) = 0
       ENDDO
 1838  A=0
       NEL(NX2) = 0
       ENDDO
       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'


c      virtual space
c      CASE) VIRT [2,1,1,1]
c      -------------------- 
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO K=1,NVIRT
       NX = NTOT+1-K
       NEL(NX) = 2
       LIM1 = NVIRT-3
       DO KK1=1,LIM1
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1845
       M1 = L1 - 1
 1845  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1
       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1846
       M2 = L2 - 1
 1846  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1
       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1847
       M3 = L3 - 1
 1847  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1850
       NEL(L3) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1851
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1851 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 1851  AAA=0
       NEL(L3)=0
       ENDDO
 1850  A=0
       NEL(L2)=0
       ENDDO
       NEL(L1)=0
       ENDDO
       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'


c      virtual space
c      CASE) VIRT [1,1,1,1,1]
c      **********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-4
       DO KK1=1,LIM1
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1
       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
       M2 = M1 - KK2
       NEL(M2) = 1
       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
       M3 = M2 - KK3
       NEL(M3) = 1
       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
       M4 = M3 - KK4
       NEL(M4) = 1
       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1855
       NEL(M5) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT5 = ICOUNT5 + 1
      call FILTER(NKEEP5,LBH,ICOUNT5,IPUT)
      IF (IPUT.EQ.0) GO TO 1856
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 1856 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 1856  AAA=0
       NEL(M5)=0
       ENDDO
 1855  A=0
       NEL(M4)=0
       ENDDO
       NEL(M3)=0
       ENDDO
       NEL(M2)=0
       ENDDO
       NEL(M1)=0
       ENDDO
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
       NEL(N5)=2
 1825  CONTINUE
       NEL(N4)=2
 1820  CONTINUE
       NEL(N3)=2
 1815  CONTINUE
       NEL(N2)=2
 1810  CONTINUE
       NEL(N1)=2
 1805  CONTINUE
 1111  A=0

 1099  CONTINUE
c      -------------------------
c      ICYC-loop (5-excitation)
c      -------------------------

       write (6,*) 'Quintuple OCCU=[1,1,1,1,1] done'
       write (6,*) 'Quintuple 5-excitations done'
       write (6,*) '       5-t: total count ICOUNT5 =',ICOUNT5
       write (6,*) '             total string count =',ICOUNT
       write (6,*) '  symmetry reduced string count =',ICOUNTS
c      -------------------------------------------------------------
       write (6,*) '                   SDTQ5-CI done'
       write (6,*) '------------------------------------------------'
c      *

c      --------------------------------------
       IF (ISIXT.EQ.0) GO TO 1000
c      --------------------------------------


 3000  A=0

       write (6,*) '-'
       write (6,*) '6-excitations: '
       write (6,*) '-'


      MAXSP6 = 3 * NKEEP6new 
      SMALL6 = 0.0000000000005D0

      write (6,*) '-'
      write (6,*) '        maximum SP6-components =',MAXSP6
      write (6,*) '                6-tuple SMALL6 =',SMALL6
      write (6,*) '-'
      write (6,*) '----------------------------------------------'
      write (6,*) 'CRITERION-1:        TNORM6 =',TNORM6
      write (6,*) 'CRITERION-2:       OKNORM6 =',OKNORM6
      write (6,*) '----------------------------------------------'

      INUMB6 = 0
      ISEL6 = 0
      ITOT6 = 0



c     ANGLE(2-2): ISTEP3 = 1 
      write (6,*) '              6-tuple selection starts here'

      CI6oldTOT = 0.0D0

c     ------------
      IFIX6 = 1 
c     ------------

      DO JJJ=1,LIMITAS6
      CIH(JJJ) = 0.0D0
      CIHsuma(JJJ) = 0.0D0
      LBH(JJJ) = 0
      LBHold(JJJ) = 0
      ENDDO



      DO 6058 ICYC6 = 1,ICYCMAX6
c     --------------------------
c     *
      write (6,*) '*'
      write (6,*) '*'
      write (6,*) '*'
      write (6,*) '                          ***************' 
      write (6,*) '                          ICYC6 =',ICYC6
      write (6,*) '                          ***************' 

      NKEEP6 = NKEEP6new
      INUMB6 = 0
      ISEL6 = 0
      ITOT6 = 0

      DO JJJ=1,LIMITAS6
      CIH(JJJ) = 0.0D0
      CIHsuma(JJJ) = 0.0D0
      LBH(JJJ) = 0
      ENDDO

      IF (ICYC6.EQ.1) GO TO 5059
c     -------------------------- 
      MAX6SP2 = MAX6SP2 + MDEL6
c     -------------------------- 
 5059 A=0

      write (6,*) '6-t: estimation starts at Doubles =',MAX6SP2 
      IF (MAX6SP2.LE.ISEL2) GO TO 6057
      write (6,*) '                   NOTE THAT MAX6SP2 > ISEL2 '
      write (6,*) '                             MAX6SP2 =',MAX6SP2 
      write (6,*) '                     MAX6SP2 > ISEL2 =',ISEL2
 6057 A=0


c     Total number of pyramidal slices:
c     ---------------------------------
      DO 6013 ICYC1 = 1,MAX6SP2
c     --------------------------
c     *
      N1 = 0
      N2 = 2
      N3 = ICYC1 + 1
      N123 = ICYC1 + 2
c*    *** Third (N3) variable ***
c     -----------------------
      DN123 = N123
      REALN3 = DN123 / 3.0D0 
      IN3 = N123 / 3 
      RN3 = IN3
      IF (RN3.LT.REALN3) GO TO 6066
c     -----------------------------
      LIM3 = IN3
      GO TO 6067
 6066 A=0
      LIM3 = IN3 + 1 
 6067 A=0
c*    *** Third (N3) variable ***
c     -----------------------
      DO 6014 ICYC2 = 1,ICYC1
      N3 = N3 - 1
      N2 = N123 - N3 - 1
      N1 = 1
      IF (N3.LT.LIM3) GO TO 6055 
c     quitting ... new slice 
c     -----------------------------
      IF (N3.GT.ISEL2) GO TO 6044 
c     skipping the N3 value...
c     -----------------------------
      IF (N3.LT.1) GO TO 6055 
c     quitting ... new slice 
c     -----------------------------
      ICYCM3 = N1 + N2
      ICYCM3 = ICYCM3 / 2
      N1 = N1 - 1
      N2 = N2 + 1
c     *** First two (N1 and N2) variables ***
c     -----------------------------------
      DO 6015 ICYC3 = 1,ICYCM3
      N1 = N1 + 1
      N2 = N2 - 1
      IF (N2.GT.ISEL2) GO TO 6040 
c     skipping the N2 value...new N1, N2
c     ----------------------------------
      IF (N2.GT.N3) GO TO 6040 
c     skip the N2 ...new (N1, N2)
c     -----------------------------
      IF (N1.GT.N2) GO TO 6044 
      IF (N2.LT.1) GO TO 6044 
c     ------------------------------
c     quitting the loop...next N3
c     ------------------------------
      KD1 = N1 
      KD2 = N2 
      KD3 = N3
c     Multiplying factor for D*D*D
      CALL MULTIPLY6(N1,N2,N3,FC6)
c     --------------------------------
      I1X = KD1 
      I2Y = KD2 
      I3Z = KD3 
      ITOT6 = ITOT6 + 1 
      LAB2 = LB2(KD1)
      LAB3 = LB2(KD2)
      LAB4 = LB2(KD3)

      CALL RLABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,
     *ICON2,LAB2)
      CALL RLABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,
     *ICON3,LAB3)
      CALL RLABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,
     *ICON4,LAB4)
      CALL NUM666(NTOT,NOCC,NVIRT,ICON2,ICON3,
     *ICON4,ICON6,IPUT)
      IF (IPUT.EQ.0) GO TO 6040 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT6,ICON6,
     *IDNUM6)
c     --------------------------
      INUMB6 = INUMB6 + 1
      IF (ISEL6.LE.NKEEP6new) GO TO 6052
c     --------------------------------

      write (6,*) 'WARNING:  SP6 exceeded limits'
      write (6,*) 'x=6:                          ISEL6 =',ISEL6
      write (6,*) ' maximum SP6              NKEEP6new =',NKEEP6new
      write (6,*) 'MUST INCREASE DIM = NKEEP6'
      write (6,*) 'x=6: cycle ICYC6 now ends'
      write (6,*) '-'
      GO TO 6050
 6052 A=0
c     *

c     --------------------------------
      IF (INUMB6.LE.MAXSP6) GO TO 6051
c     --------------------------------
      write (6,*) '---'
      write (6,*) 'WARNING: SP6 exceeded maximum-components '
      write (6,*) '6-tuple: numb. of components INUMB6 =',INUMB6
      write (6,*) '# maximum components         MAXSP6 =',MAXSP6
      write (6,*) '         # components = 3 * NKEEP6'
      write (6,*) '    MUST INCREASE DIM = NKEEP6'
      write (6,*) 'x=6: cycle ICYC6 now ends'
      write (6,*) '---'
      GO TO 6050
 6051 A=0

      PD6 = CI2EST(KD1)*CI2EST(KD2)*CI2EST(KD3)
      PD6suma = CI2suma(KD1)*CI2suma(KD2)*CI2suma(KD3)
      PD6 = PD6 * FC6 
      PD6suma = PD6suma * FC6 

      IF (INUMB6.EQ.1) GO TO 6035
c     check if IDNUM6 has occured before:
      DO ITT=1,ISEL6
      IHIT = ITT
      IF (IDNUM6.NE.LBH(ITT)) GO TO 6030 
      CIH(IHIT) = CIH(IHIT) + PD6 
      CIHsuma(IHIT) = CIHsuma(IHIT) + PD6suma 
      GO TO 6040
 6030 A=0
      ENDDO
 6035 A=0

c     (6-tuple is new):
      ISEL6 = ISEL6 + 1
      LBH(ISEL6) = IDNUM6
      CIH(ISEL6) = CIH(ISEL6) + PD6 
      CIHsuma(ISEL6) = CIHsuma(ISEL6) + PD6suma 
c-------------------------------------------
c*****LAIMIS:ISEL6<mxsprod
c-------------------------------------------
      IF (ISEL6.LT.mxsprod) GO TO 6040
      write (6,*) 'PROBLEM: ISEL6 is ',ISEL6
      write (6,*) '       mxsprod is ',mxsprod
      write (6,*) 'NOW PROGRAM WILL STOP '
      STOP
 6040 A=0
      IF (N1.GT.N2) GO TO 6044

 6015 CONTINUE
c     --------
c     *** (N1 and N2 variables)

 6044 A=0
 6014 CONTINUE
c     --------
c     *** (N3 variable)
 6055 A=0
 6013 CONTINUE
c     --------
c     (Number of pyramidal slices)


 6050 A=0
c*****write (6,*) '         Final values:KD1,KD2,KD3:',KD1,KD2,KD3 


c     6-tuple selection completed
      write (6,*) '-------------------------------------------------'
      write (6,*) '6-tuple selection completed'
      write (6,*) '         Total (D*D*D) products  ITOT6 =',ITOT6
      write (6,*) ' number of SP6 components       INUMB6 =',INUMB6
      write (6,*) '        Distinct SP6             ISEL6 =',ISEL6
      write (6,*) '                                NKEEP6 =',NKEEP6
      write (6,*) '-------------------------------------------------'
      IF (ISEL6.GE.NKEEP6) GO TO 6079
      NKEEP6 = ISEL6
      write (6,*) '   NOTE ISEL6 < NKEEP6: WILL SET NKEEP6=ISEL6'
      write (6,*) '                                NKEEP6 =',NKEEP6
 6079 A=0


c            -------------------
c     ****** ORDER 6-tuples now: *****
c            -------------------
      IF (ISEL6.GE.IFIX6) GO TO 6080
      write (6,*) 'ISEL6 is too small'
      STOP
 6080 A=0
      IF (ICYC6.GT.1) GO TO 6081
      IFIX6 = ISEL6
      SUM6old = 0.0D0
      DO JJA = 1,IFIX6
      LBHold(JJA) = LBH(JJA)
      SUM6old = SUM6old + CIHsuma(JJA)
      ENDDO
 6081 A=0.0
c     -----------------------------------
      call SORT2(ISEL6,CIH,CIHsuma,LBH)
c     -----------------------------------
      DO JJJ=1,ISEL6
      J = ISEL6 + 1 - JJJ
      IF (JJJ.GE.J) GO TO 6089
      A1 = CIH(JJJ)
      A2 = CIH(J)
      B1 = CIHsuma(JJJ)
      B2 = CIHsuma(J)
      IC1 = LBH(JJJ)
      IC2 = LBH(J) 

      CIH(JJJ) = A2
      CIH(J) = A1
      CIHsuma(JJJ) = B2
      CIHsuma(J) = B1
      LBH(JJJ) = IC2 
      LBH(J) = IC1 
      ENDDO
 6089 A=0
      CI6newTOT = 0.0D0 
      DO J=1,ISEL6
      CI6newTOT = CI6newTOT + CIHsuma(J)
      ENDDO
      write (6,*) 'NEW ordered set for 6-tuples:' 

      CI6def = CI6newTOT
      CI6norm = 0.0D0
      CI6norm2 = 0.0D0
      ICUT6 = 0
      ICUT6 = 0

      DO J=1,ISEL6
      ICUT6 = ICUT6 + 1
      CI6norm2 = CI6norm2 + CIHsuma(J)
      RATIO6 = CI6norm2 / CI6newTOT 
      CI6last = CIHsuma(J)
      IF (RATIO6.GE.TNORM6) GO TO 6098
      ENDDO
      write (6,*) 'Will stop. SP6 norm6 has problem'
      write (6,*) 'Triangle is too small !!! '
      STOP
 6098 A=0

      NKEEP6 = ICUT6

      CI6def = CI6newTOT

      DO J=1,NKEEP6
      CI6def = CI6def - CIHsuma(J)
      CI6norm = CI6norm + CIHsuma(J)
      ENDDO

      RELDEF6 = CI6def / CI6newTOT
      REL6norm = CI6norm / CI6newTOT


      write (6,*) '                       ----------------'
      write (6,*) '                         CRITERION-1'
      write (6,*) '                       ----------------'
      write (6,*) '                         ICYC6 = ',ICYC6
      write (6,*) ' Number of SP6:          ISEL6 =',ISEL6
      write (6,*) '                         ICUT6 =',NKEEP6
      write (6,*) ' will keep              NKEEP6 =',NKEEP6
      write (6,*) 'The ratio (internal)    RATIO6 =',RATIO6
      write (6,*) 'The target is           TNORM6 =',TNORM6
      write (6,*) 'total 6-norm in NKEEP6 CI6norm =',CI6norm
      write (6,*) '6-deficiency            CI6def =',CI6def
      write (6,*) 'RELATIVE 6-deficiency  RELDEF6 =',RELDEF6

      RAT6oldnew = CI6oldTOT / CI6newTOT

      write (6,*) '                       ----------------'
      write (6,*) '                         CRITERION-2'
      write (6,*) '                       ----------------'
      write (6,*) '                         ICYC6 = ',ICYC6
      write (6,*) 'Actual ratio old/new           =',RAT6oldnew
      write (6,*) 'The target  is                 =',OKNORM6
      write (6,*) 'Total Norm-6 (new)   CI6newTOT =',CI6newTOT 
      write (6,*) 'Total Norm-6 (old)   CI6oldTOT =',CI6oldTOT 
      IF (RAT6oldnew.GE.OKNORM6) GO TO 6088 
      CI6oldTOT = CI6newTOT
 6058 CONTINUE
c     Enlarged pyramid
c     ----------------
      write (6,*) 'PROGRAM WILL STOP!'
      write (6,*) 'Desired accuracy not yet achieved...'
      STOP

 6088 A=0
      write (6,*) '6-tuples CRITERION-2 is satisfied'
      TNORMAL = TNORMAL + CI6norm
      ALLNORMAL = ALLNORMAL + CI6newTOT
      DIFFER = ALLNORMAL - TNORMAL
      RELDEFIC = DIFFER / ALLNORMAL

      write (6,*) '--------------------------------------------'
      write (6,*) 'CRITERION-1: keep only SPs NKEEP6 =',NKEEP6
      write (6,*) 'SELECTED NORM-SUM: TNORMAL(4,5,6) =',TNORMAL
      write (6,*) 'TOTAL NORM-SUM    ALLNORMAL(4,5,6)=',ALLNORMAL
      write (6,*) '--------------------------------------------'
      write (6,*) '   RELATIVE-DEFICIENCY(4,5,6)     =',RELDEFIC
      write (6,*) '--------------------------------------------'
      write (6,*) '       At SDTQ5-CI level there are:'
      write (6,*) 'Symmetry-reduced strings ICOUNTS =',ICOUNTS
      write (6,*) '--------------------------------------------'


       ICYCMAX = 1 

c***   DO 2099 ICYC=1,ICYCMAX
c      ----------------------


       write (6,*) '6-tuple SP selection is STEP =',ICYC
       ICOUNT6 = 0

       write (6,*) 'Sixtuple OCCU=[0,0,0] starts'
c     ----------------------
c     SIXTUPLE-EXCITATIONS:
c     ----------------------
C      occupied MO space: CASE  [0,0,0]
c     ---------------------------------- 
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO 3105 I=1,NOCC-2
       M = NOCC+1-I
       NEL(M) = 0
       DO 3110 J=1,M-2
       MM = M - J
       NEL(MM) = 0
       DO 3115 JJ=1,MM-1
       MMM = MM - JJ
       NEL(MMM) = 0

c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      CASE) VIRT[2,2,2]
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
       MK = NTOT+1-KK
       NEL(MK) = 2
       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
       MMK = MK - KKK
       NEL(MMK) = 2
       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
       MMK2 = MMK - KKK2
       NEL(MMK2) = 2
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1
      call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3232 
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3232 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3232  AAA=0
       NEL(MMK2) = 0
       ENDDO
       NEL(MMK) = 0
       ENDDO
       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'

c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      CASE) VIRT [2,2,1,1]
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-1
       DO K1=1,LIM1
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 
       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
       NX2 = NX1 - K2
       NEL(NX2) = 2 
       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3235
       M1 = L1 - 1
 3235  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1
       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3236
       LIM4 = LIM4 - 1
 3236  A=0
       DO KK2=1,LIM4
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3237
       M2 = L2 - 1
 3237  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3238 
       NEL(L2) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3240 
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3240 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3240  AAA=0
       NEL(L2) = 0
       ENDDO
 3238  A=0
       NEL(L1) = 0
       ENDDO
       NEL(NX2) = 0
       ENDDO
       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       DO K=1,NVIRT
       NX = NTOT+1-K
       NEL(NX) = 2
       LIM1 = NVIRT-4
       DO KK1=1,LIM1
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3245
       M1 = L1 - 1
 3245  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1
       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3246
       M2 = L2 - 1
 3246  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1
       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3247
       M3 = L3 - 1
 3247  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1
       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3248
       M4 = L4 - 1
 3248  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3250
       NEL(L4) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3251 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3251 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3251  AAA=0
       NEL(L4)=0
       ENDDO
 3250  A=0
       NEL(L3)=0
       ENDDO
       NEL(L2)=0
       ENDDO
       NEL(L1)=0
       ENDDO
       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'

c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-5
       DO KK1=1,LIM1
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1
       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
       M2 = M1 - KK2
       NEL(M2) = 1
       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
       M3 = M2 - KK3
       NEL(M3) = 1
       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
       M4 = M3 - KK4
       NEL(M4) = 1
       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
       M5 = M4 - KK5
       NEL(M5) = 1
       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3255
       NEL(M6) = 1
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3256 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3256 
       ICOUNTS = ICOUNTS + 1
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3256  AAA=0
       NEL(M6)=0
       ENDDO
 3255  A=0
       NEL(M5)=0
       ENDDO
       NEL(M4)=0
       ENDDO
       NEL(M3)=0
       ENDDO
       NEL(M2)=0
       ENDDO
       NEL(M1)=0
       ENDDO
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
c      -----------------------------
c      end of occupied MOs: CASE [0,0,0]
c      -----------------------------
       NEL(MMM) = 2
 3115  CONTINUE
       NEL(MM) = 2
 3110  CONTINUE
       NEL(M) = 2
 3105  CONTINUE
       write (6,*) 'Sixtuple OCCU=[0,0,0] done'
 3400  A=0


c      SCF-SPACE:
       write (6,*) 'Sixtuple OCCU[1,1,0,0] starts'
C      occupied MO space: CASE  [1,1,0,0]
c      ***********************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC-1 
       DO 3500 I1=1,MAX1
c      ------------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 0



       MAX2 = N1-1
       DO 3505 I2=1,MAX2
c      -------------------
c      *
       N2 = N1 - I2
       NEL(N2) = 0



       MAX3 = NOCC-3
       DO 3510 I3=1,MAX3
c      ------------------
c      *
       NN3 = NOCC + 1 - I3
       IF (I3.EQ.1) GO TO 3511
       NN3 = LL3 - 1
 3511  A=0
       call SHIFT2(N1,N2,NN3,LL3)
       NEL(LL3) = 1




       MAX4 = NN3 - 1  
 3512  A=0
       DO 3515 I4=1,MAX4
c      ------------------
c      *
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 3516
       NN4 = LL4 - 1
 3516  A=0
       call SHIFT2(N1,N2,NN4,LL4)
       IF (LL4.EQ.0) GO TO 3517
       NEL(LL4) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      --------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3532 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3532 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3532  AAA=0
c      ----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3535
       M1 = L1 - 1
 3535  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3536
       LIM4 = LIM4 - 1
 3536  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3537
       M2 = L2 - 1
 3537  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3538 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3540 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3540 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3540  AAA=0
c      ----------------------------------
       NEL(L2) = 0
       ENDDO
 3538  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3545
       M1 = L1 - 1
 3545  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3546
       M2 = L2 - 1
 3546  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3547
       M3 = L3 - 1
 3547  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3548
       M4 = L4 - 1
 3548  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3550
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3551 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3551 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3551  AAA=0
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3550  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3555
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3556 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3556 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3556  AAA=0
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3555  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
C*     -----------------------------------

       NEL(LL4)=2
 3515  CONTINUE
 3517  A=0

 
       NEL(LL3)=2
 3510  CONTINUE


       NEL(N2)=2
 3505  CONTINUE


       NEL(N1)=2
 3500  CONTINUE
c      ***

       write (6,*) 'Sixtuple OCCU[1,1,0,0] done'






c      ********
 3600  A=0
c      ********


c      SCF-SPACE:
       write (6,*) 'Sixtuple OCCU[1,1,1,1,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,1,1,1,0]
c      ***********************************

       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO

       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC
       DO 3605 I1=1,NOCC
c      -----------------
c      *
       N1 = NOCC + 1 - I1
       NEL(N1) = 0


       MAX2 = NOCC-4
       DO 3610 I2=1,MAX2
c      ------------------
c      *
       NN2 = NOCC + 1 - I2
       IF (I2.EQ.1) GO TO 3611
       NN2 = LL2 - 1
 3611  A=0
       call SHIFT1(N1,NN2,LL2)
       NEL(LL2) = 1



       MAX3 = NN2 - 3
       DO 3615 I3=1,MAX3
c      -----------------
c      *
       NN3 = LL2 - I3
       IF (I3.EQ.1) GO TO 3616
       NN3 = LL3 - 1
 3616  A=0
       call SHIFT1(N1,NN3,LL3)
       NEL(LL3) = 1




       MAX4 = NN3 - 2
       DO 3620 I4=1,MAX4
c      -----------------
c      *
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 3621
       NN4 = LL4 - 1
 3621  A=0
       call SHIFT1(N1,NN4,LL4)
       NEL(LL4) = 1


       MAX5 = NN4 - 1
       DO 3625 I5=1,MAX5
c      ------------------
c      *
       NN5 = LL4 - I5
       IF (I5.EQ.1) GO TO 3626
       NN5 = LL5 - 1
 3626  A=0
       call SHIFT1(N1,NN5,LL5)
       IF (LL5.EQ.0) GO TO 3627
       NEL(LL5) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      --------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3632 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3632 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3632  AAA=0
c      ----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3635
       M1 = L1 - 1
 3635  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3636
       LIM4 = LIM4 - 1
 3636  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3637
       M2 = L2 - 1
 3637  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3638 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3640 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3640 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3640  AAA=0
c      ----------------------------------
       NEL(L2) = 0
       ENDDO
 3638  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3645
       M1 = L1 - 1
 3645  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3646
       M2 = L2 - 1
 3646  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3647
       M3 = L3 - 1
 3647  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3648
       M4 = L4 - 1
 3648  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3650
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3651 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3651 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3651  AAA=0
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3650  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3655
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3656 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3656 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3656  AAA=0
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3655  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
C*     -----------------------------------

       NEL(LL5)=2
 3625  CONTINUE
 3627  A=0


       NEL(LL4)=2 
 3620  CONTINUE


       NEL(LL3)=2
 3615  CONTINUE


       NEL(LL2)=2
 3610  CONTINUE


       NEL(N1)=2
 3605  CONTINUE
c      ******
       write (6,*) 'Sixtuple OCCU[1,1,1,1,0] done'






c      ********
 3800  A=0
c      ********




c      SCF-SPACE:
       write (6,*) 'Sixtuple OCCU[1,1,1,1,1,1] starts'
c
c      **************************************
C      occupied MO space: CASE  [1,1,1,1,1,1]
c      **************************************

       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       MAX1 = NOCC-5
       DO 3805 I1=1,MAX1
c      -----------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 1

       MAX2 = N1 - 5 
       DO 3810 I2=1,MAX2
c      -----------------
c      *
       N2 = N1 - I2
       NEL(N2) = 1

       MAX3 = N2 - 4
       DO 3815 I3=1,MAX3
c      ------------------
c      *
       N3 = N2 - I3
       NEL(N3) = 1

       MAX4 = N3 - 3
       DO 3820 I4=1,MAX4
c      -----------------
c      *
       N4 = N3 - I4
       NEL(N4) = 1

       MAX5 = N4 - 2 
       DO 3825 I5=1,MAX5
c      -----------------
c      *
       N5 = N4 - I5
       NEL(N5) = 1

       MAX6 = N5 - 1
       DO 3830 I6=1,MAX6
c      ------------------
c      *
       N6 = N5 - I6
       IF (N6.GT.0) GO TO 3831 
       write (6,*) 'N6=',N6
       STOP
 3831  A=0
       NEL(N6) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      *************************
c      --------------------
       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3832 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3832 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3832  AAA=0
c      ----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      ***********************************


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3835
       M1 = L1 - 1
 3835  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3836
       LIM4 = LIM4 - 1
 3836  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3837
       M2 = L2 - 1
 3837  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3838 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3840 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3840 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3840  AAA=0
c      ----------------------------------
       NEL(L2) = 0
       ENDDO
 3838  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3845
       M1 = L1 - 1
 3845  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3846
       M2 = L2 - 1
 3846  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3847
       M3 = L3 - 1
 3847  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3848
       M4 = L4 - 1
 3848  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3850
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3851 
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3851 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3851  AAA=0
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3850  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3855
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
       ICNT = 1
       ICOUNT6 = ICOUNT6 + 1

c******call FILTER(NKEEP6,LB6,ICOUNT6,IPUT)
       call FILTER(NKEEP6,LBH,ICOUNT6,IPUT)
       IF (IPUT.EQ.0) GO TO 3856
c      -------------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
      CALL IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
      CALL SORT(NALP,IALP)
      CALL SORT(NBET,IBET)
       IPOSA(ICNT) = posdet(NTOT,NALP,IALP,ifa)
       IPOSB(ICNT) = posdet(NTOT,NBET,IBET,ifa)
       call getsym1(IALP,NTOT,NALP,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMA(ICNT) = isym 
       call getsym1(IBET,NTOT,NBET,IRREP,idsym,isym,
     * iwrk(1),iwrk(4),iwrk(7),iwrk(10))
       ISYMB(ICNT) = isym 
       IF (KTAB(ISYMA(ICNT)).NE.ISYMB(ICNT)) GO TO 3856 
c      **************************
c
c      ---------------------
       ICOUNTS = ICOUNTS + 1
c      ---------------------
       ISTRA(ICOUNTS) = IPOSA(ICNT)
       ISTRB(ICOUNTS) = IPOSB(ICNT)
 3856  AAA=0
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3855  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *************************************
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'



       NEL(N6)=2
 3830  CONTINUE


       NEL(N5)=2
 3825  CONTINUE


       NEL(N4)=2
 3820  CONTINUE


       NEL(N3)=2
 3815  CONTINUE


       NEL(N2)=2
 3810  CONTINUE


       NEL(N1)=2
 3805  CONTINUE




       write (6,*) 'Sixtuple OCCU[1,1,1,1,1,1] done'
       write (6,*) '-------------------------'





c     --------------
c     6-excitations:
c     --------------
c      ***
 2111  A=0
c      ***




c      ****************
C***** 2099  CONTINUE
c      ****************
c      -------------------------
c      ICYC-loop (6-excitation)
c      -------------------------
     

       write (6,*) 'Sixtuple-excitations  done'
       write (6,*) 'CI-SDTQ56 GCI selection  done'
       write (6,*) '------------------------------------'
       write (6,*) ' total strings =',ICOUNT
       write (6,*) ' symmetry reduced strings =',ICOUNTS
       write (6,*) '------------------------------------'



c*******************END OF SDTQ56-excitations************


 1000  AAA=0
      DIFFER = ALLNORMAL - TNORMAL
      RELDEFIC = DIFFER / ALLNORMAL


      write (6,*) '                 FINAL RESULT: '
      write (6,*) '--------------------------------------------'
      write (6,*) 'SELECTED NORM-SUM: TNORMAL(4,5,6) =',TNORMAL
      write (6,*) 'TOTAL NORM-SUM:  ALLNORMAL(4,5,6) =',ALLNORMAL
      write (6,*) '--------------------------------------------'
      write (6,*) '    RELATIVE-DEFICIENCY(4,5,6)    =',RELDEFIC
      write (6,*) '--------------------------------------------'
      write (6,*) '*'
      write (6,*) '*'
      write (6,*) ' *****************************************'
      WRITE (6,*) ' TOTAL SUMMARY (COMPACT CI-SDTQ56) FOR GCI'
      write (6,*) ' *****************************************'
      write (6,*) '*'
      write (6,*) '--------------------------------------------'
      write (6,*) 'RELATIVE-DEFICIENCY(4,5,6) =',RELDEFIC
      write (6,*) '--------------------------------------------'
      write (6,*) '*'
      write (6,*) 'Triples SP3 kept    NKEEP3 =',NKEEP3
      write (6,*) 'out of               ISEL3 =',ISEL3
      write (6,*) '*'
      write (6,*) 'Quadruple SP4 kept  NKEEP4 =',NKEEP4
      write (6,*) 'out of               ISEL4 =',ISEL4
      write (6,*) '*'
      write (6,*) 'Pentuple SP5 kept   NKEEP5 =',NKEEP5
      write (6,*) 'out of               ISEL5 =',ISEL5
      write (6,*) '*'
      write (6,*) 'Sextuple SP6 kept   NKEEP6 =',NKEEP6
      write (6,*) 'out of               ISEL6 =',ISEL6
      write (6,*) '--------------------------------------------'

 1001 A=0
      WRITE (6,*) 'all space-products ICOUNT=',ICOUNT
      WRITE (6,*) 'symmetry-modified ICOUNTS=',ICOUNTS
      write (6,*) '--------------------------------------------'

c       WRITE (6,*) (ISTRA(L), L=1,ICOUNTS)
c       WRITE (6,*) (ISTRB(L), L=1,ICOUNTS)


C*     Writing-out the string positions to the FILE
C      ---------------------------------------------------
       open (unit=11,file='listGCI',status='unknown',
     *      access='sequential',form='unformatted')


C*     'alpha-string array ISTRA'
C*      'beta-string array ISTRB '

       WRITE (11) ICOUNTS
       WRITE (11) (ISTRA(L), L=1,ICOUNTS)
       WRITE (11) (ISTRB(L), L=1,ICOUNTS)
c
       close (11)



C----------------------------------------------------------
C
       STOP
       END










c     -----------------------------------------------
      subroutine MATCH6(IFIX6,LB6old,LB6,isame,perc6)
c     -----------------------------------------------
       implicit double precision(a-h,o-z)
c
       parameter (mxsprod=2400000)
       DIMENSION LB6(mxsprod)
       DIMENSION LB6old(mxsprod)
c      ------------------------

       isame = 0


       DO J=1,IFIX6
c      ------------


       DO JJ=1,IFIX6
c      -------------
       IF (LB6(JJ).NE.LB6old(J)) GO TO 50
       isame = isame + 1
       GO TO 100
  50   A=0
c      ***
c      -----
       ENDDO


 100   A=0
c      ***
c      -----
       ENDDO

       same = isame
       fix6 = IFIX6

       perc6 = same / fix6

C      -----------
       RETURN
       END







c      ----------------------------------------------------
c
       SUBROUTINE FILTER(NKEEP4,LB4,ICOUNT4,IPUT)
c
c      ----------------------------------------------------
       implicit double precision(a-h,o-z)
c
       parameter (mxsprod=2400000)
       DIMENSION LB4(mxsprod)
c
c      ------------------------------
c      NKEEP4 = # of space-products kept
c      LB4    = the list of s-p labels that
c              correspond to these s-p
c      ICOUNT4 = label to be checked-out
c      *** 
c      IPUT = 0 (no match)
c      IPUT = 1 (match)
c      --------------------------------


       IPUT = 0


       DO 100 J = 1,NKEEP4
c      --------------------
c      *
       IF (LB4(J).NE.ICOUNT4) GO TO 150
c      --------------------------------
c      match: 
c      --------------
       IPUT = 1
c      --------------
       GO TO 200

 150   A=0
c      ***
 100   CONTINUE


 200   A=0
C      -----------
       RETURN
       END




c      ----------------------------------------------------
c
       SUBROUTINE FILTER3(NKEEP4,LB4,ICOUNT4,IPUT)
c
c      ----------------------------------------------------
       implicit double precision(a-h,o-z)
c
       parameter (mxsprod1=100000)
       DIMENSION LB4(mxsprod1)
c
c      ------------------------------
c      NKEEP4 = # of space-products kept
c      LB4    = the list of s-p labels that
c              correspond to these s-p
c      ICOUNT4 = label to be checked-out
c      *** 
c      IPUT = 0 (no match)
c      IPUT = 1 (match)
c      --------------------------------


       IPUT = 0


       DO 100 J = 1,NKEEP4
c      --------------------
c      *
       IF (LB4(J).NE.ICOUNT4) GO TO 150
c      --------------------------------
c      match: 
c      --------------
       IPUT = 1
c      --------------
       GO TO 200

 150   A=0
c      ***
 100   CONTINUE


 200   A=0
C      -----------
       RETURN
       END




c      ----------------------------------------------------
c
       SUBROUTINE SHIFT1(NX,M1,L1)
c
c      ----------------------------------------------------
       implicit double precision(a-h,o-z)
C
c      ***
       IF (M1.GT.NX) GO TO 551 
       IF (M1.EQ.NX) GO TO 552
       IF (M1.LT.NX) GO TO 551
c      ***
 551   AAA=0
       L1 = M1
       GO TO 553

 552   AAA=0
       L1 = M1 - 1
 553   AAA=0
c      ***
C      -----------
       RETURN
       END







c      ----------------------------------------------------

       SUBROUTINE SHIFT2(NX1,NX2,M1,L1)

c      ----------------------------------------------------
       implicit double precision(a-h,o-z)

c      ***
       IF (M1.GT.NX1) GO TO 551 
       IF (M1.EQ.NX1) GO TO 552
       IF (M1.LT.NX1) GO TO 553
 551   AAA=0
       L1 = M1
       GO TO 560

 553   A=0
       L1 = M1
       IF (M1.NE.NX2) GO TO 560 
       L1 = M1 - 1




 552   AAA=0
       M1 = M1 - 1
       IF (M1.EQ.NX2) GO TO 555
       L1 = M1 
       GO TO 560


 555   AAA=0
       L1 = M1 - 1 


 560   AAA=0
c      ***
C      -----------
       RETURN
       END







C      -------------
C*     SUBROUTINES:
C      -------------


      SUBROUTINE IDENTIFY(ID,IS,NEL,NTOT,ND,NS,IALP,IBET,
     *NALP,NBET,MBET,MALP)
C     ----------------------------------------------
C     This subroutine identifies singly and doubly
C     occupied MOs and stores them as an array 
C     ID(icount,j) and IS(icount,k)
C     -----------------------------------------------------------------
      implicit double precision(a-h,o-z)
      parameter (mxorb = 100)
      DIMENSION NEL(mxorb),ID(mxorb),IS(mxorb),IALP(mxorb),IBET(mxorb)

       NS = 0
       ND = 0
       DO 40 J=1,NTOT
       IF (NEL(J).EQ.0) GO TO 50 
       IF (NEL(J).GT.1) GO TO 45 

C      singly
       NS = NS + 1
       IS(NS) = J
       GO TO 50 

   45  AAA=0.0
C      doubly
       ND = ND + 1
       ID(ND) = J
   50  AAA=0.0
   40  CONTINUE


C      ----------------------------
C*     count ALP and BET electrons
C      ----------------------------
       NALP = 0
       NBET = 0
       DO 60 J=1,ND
       NALP = NALP + 1
       NBET = NBET + 1
       IALP(NALP) = ID(J)
       IBET(NBET) = ID(J)
  60   CONTINUE



       DO 70 J=1,NS
       IF (NALP.EQ.MALP) GO TO 75 
       NALP = NALP + 1
       IALP(NALP) = IS(J)
       GO TO 77 
 75    AAA=0.0
       IF (NBET.EQ.MBET) GO TO 77 
       NBET = NBET + 1
       IBET(NBET) = IS(J)
 77    AAA=0.0
  70   CONTINUE


C      -----------
       RETURN
       END












C*------------------------------------------------------------
c      SUBROUTINE WRITE(ID,IS,NEL,NTOT,ND,NS,ICOUNT,IALP,IBET,
c     *NALP,NBET,MBET,MALP)
c      implicit double precision(a-h,o-z)
c      parameter (mxorb = 100)
c      DIMENSION NEL(mxorb),ID(mxorb),IS(mxorb)
c      DIMENSION IALP(mxorb),IBET(mxorb)
C--------------------------------------------------------------
c       WRITE (6,25) ICOUNT
c  25   FORMAT ('space-product = ',I5)
c       WRITE (6,30) ND
c  30   FORMAT ('doubly occupied MOs',I5)
c       WRITE (6,55) (ID(I), I=1,ND) 
c  55   FORMAT (9I2)
c       WRITE (6,60) NS
c  60   FORMAT ('singly occupied MOs',I5)
c       WRITE (6,96) (IS(I), I=1,NS) 
c  96   FORMAT (9I2)
c       WRITE (6,130) NALP
c 130   FORMAT ('alpha electr',I5)
c       WRITE (6,155) (IALP(I), I=1,NALP) 
c 155   FORMAT (9I2)
c       WRITE (6,160) NBET 
c 160   FORMAT ('beta electrons',I5)
c       WRITE (6,196) (IBET(I), I=1,NBET) 
c 196   FORMAT (9I2)
C       WRITE (6,200)
C 200   FORMAT ('-------------')
C      -----------
c       RETURN
c       END
C*------------------------------------------------------------







C*    -----------------------------------
      SUBROUTINE SORT(N,IARR)
C*    -----------------------------------
      implicit double precision(a-h,o-z)
C
C     sort an array IARR(N) into ascending numerical
C     order
C     ----------------------------------------------
      parameter (mxorb = 100)
      DIMENSION IARR(mxorb)

      DO 12 J=2,N
      IA = IARR(J)
      DO 11 I=J-1,1,-1
      IF (IARR(I).LE.IA) GO TO 10
      IARR(I+1) = IARR(I)
  11  CONTINUE
      I=0
  10  IARR(I+1) = IA
  12  CONTINUE
      RETURN
      END





c
c     -------------------------------------------------------------
      subroutine ibinom(ifa,n,na,nb,nalp,nblp)
c     -------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
      integer ifa(0:n,0:n)
c
c     Returns all binomial numbers (i,j) for i=0,n and j=0,i in fa .
c     The binomial number (i,j) is stored in ifa(i,j)
c

c     ----------------------------------------------------
C     ifa  = binomial number specifying alpha/beta string 
c     n    = number of active MOs
c     na   = number of alpha electrons 
c     nb   = number of beta electrons 
c     nalp = number of alpha strings
c     nblp = number of beta strings
C     ----------------------------------------------------
C
      do 13 ii=0,n
	 ifa(ii,0)  = 1
	 ifa(ii,ii) = 1
   13 continue
c
      do 113 iy = 1, n
	 do 114 ix = 1, (iy-1)
	    ifa(iy,ix) = ifa(iy-1,ix-1) + ifa(iy-1,ix)
  114    continue
  113 continue
c
      nalp = ifa(n,na)
      nblp = ifa(n,nb)
c
      return
      end
c







c     -------------------------------------------
      integer function posdet(nact,noel,con,ifa)
c     -------------------------------------------
      implicit double precision(a-h,o-z)
      integer con(noel),els,pos1,b,a,pos2,k,i,j
      dimension ifa(0:nact,0:nact)
c
c     ----------------------------------------------
C     nact = number of active MOs
c     noel = number of electrons (alpha or beta)
c     con  = (IALP/IBET MOs label) for alpha or beta MO
c     ---------------------------------------------


      pos1 = 0
      posdet = 1
      do 33 i=1,noel
         do 55 j=pos1+1,con(i)-1
            posdet = posdet + ifa(nact-j,noel-i)
   55    continue
         pos1 = con(i)
   33 continue
c
      return
      end
c


c
c     -----------------------------------------------------
      subroutine gtab(idsym,isym1,itab,iele,ista,iscr,icha)
c     -----------------------------------------------------
c
c     Routine to return table such that i x itab(i) = isym1
c     where i is an irreducible representation.
c     idsym  specifies the point group.
c     isym1 is the desired symmetry
c     itab = returned table explained above.
c 
c     iele = scratch integer array of length 3
c     ista = scratch integer array of length 3
c     iscr = scratch integer array of length 3
c     icha = scratch integer array of length 34      
c     
c     
c
c     CONVENTION FOR idsym,isym1, and itab
c
c     Point group  idsym  Irred rep isym1  Sym operations used
c   -----------------------------------------------------------
c        Ci          1       Ag       1        i
c                            Au       2
c    
c        Cs          1       A'       1      (sigma)h
c                            A''      2      
c       
c        C2          1       A        1       C2
c                            B        2
c
c        D2          2       A        1       C2(z)
c                            B1       2       C2(y)
c                            B2       3     
c                            B3       4
c
c        C2v         2       A1       1       C2
c                            A2       2       (sigma)v(xz)
c                            B1       3  
c                            B2       4
c
c        C2h         2       Ag       1       i
c                            Bg       2       (sigma)h
c                            Bu       3
c                            Au       4
c
c        D2h         3       Ag       1       (sigma)(xy)
c                            B1g      2       (sigma)(xz)
c                            B2g      3       (sigma)(yz)
c                            B3g      4 
c                            Au       5
c                            B1u      6
c                            B2u      7
c                            B3u      8
c
c
      implicit double precision(a-h,o-z)
      dimension iele(3),ista(3),icha(34)
      dimension iscr(3)
      dimension itab(*)
      call getdata(iele,ista,icha)
c      data (iele(i),i=1,3) /2,4,8/
c      data (ista(i),i=1,3) /1,3,11/
c      data (icha(i),i=1,34) /1,-1,1,1,1,-1,-1,1,-1,-1,
c     *   1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,
c     *   1,-1,1,1,1,-1/
c
      ist = ista(idsym)
      iel = iele(idsym)
      call gtab1(icha(ist),iel,idsym,isym1,itab,iscr)
      return
      end
c


c     ----------------------------------------------------
      subroutine gtab1(icha,iel,idi,isym1,itab,iscr)
c     ----------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension iscr(3)
      dimension icha(idi,iel)
      dimension itab(iel)
c
      do 34 ii=1,iel
         do 45 jj=1,iel
	    do 77 kk=1,idi
	       iscr(kk) = icha(kk,ii)*icha(kk,jj)
	       if (iscr(kk).ne.icha(kk,isym1)) goto 45
   77       continue
	    itab(ii) = jj
	    goto 34
   45    continue
   34 continue
c
      return
      end








c
c     -------------------------------------------------
      subroutine gmul(idsym,imul,iele,ista,iscr,icha)
c     -------------------------------------------------
c
c     Routine to return multiplication table ixj = imul(i,j)
c     where i,j are irreducible representations.
c     idsym  specifies the point group.
c
c     CONVENTION FOR idsym,i,j and imul
c
c     Point group  idsym  Irred rep  i,j  Sym operation used
c   -----------------------------------------------------------
c        Ci          1       Ag       1        i
c                            Au       2
c    
c        Cs          1       A'       1      (sigma)h
c                            A''      2      
c       
c        C2          1       A        1       C2
c                            B        2
c
c        D2          2       A        1       C2(z)
c                            B1       2       C2(y)
c                            B2       3     
c                            B3       4
c
c        C2v         2       A1       1       C2
c                            A2       2       (sigma)v(xz)
c                            B1       3  
c                            B2       4
c
c        C2h         2       Ag       1       i
c                            Bg       2       (sigma)h
c                            Bu       3
c                            Au       4
c
c        D2h         3       Ag       1       (sigma)(xy)
c                            B1g      2       (sigma)(xz)
c                            B2g      3       (sigma)(yz)
c                            B3g      4 
c                            Au       5
c                            B1u      6
c                            B2u      7
c                            B3u      8
c
c
      implicit double precision(a-h,o-z)
      dimension iele(3),ista(3),icha(34)
      dimension iscr(3)
      dimension imul(*)
      call getdata(iele,ista,icha)
c      data (iele(i),i=1,3) /2,4,8/
c      data (ista(i),i=1,3) /1,3,11/
c      data (icha(i),i=1,34) /1,-1,1,1,1,-1,-1,1,-1,-1,
c     *   1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,
c     *   1,-1,1,1,1,-1/
c
      ist = ista(idsym)
      iel = iele(idsym)
      call gmul1(icha(ist),iel,idsym,imul,iscr)
      return
      end
c
c     ----------------------------------------------------
      subroutine gmul1(icha,iel,idi,imul,iscr)
c     ----------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension iscr(idi)
      dimension icha(idi,iel)
      dimension imul(iel,iel)
c
      do 34 ii=1,iel
         do 45 jj=1,iel
	    do 77 kk=1,idi
	       iscr(kk) = icha(kk,ii)*icha(kk,jj)
   77       continue
	    do 88 kl=1,iel
	       do 99 lk=1,idi
		  if (iscr(lk).ne.icha(lk,kl)) goto 88
   99          continue
	       imul(ii,jj) = kl
	       goto 45
   88       continue 
   45    continue
   34 continue
c
      return
      end
c
c     -------------------------------------------------
      subroutine getsym1(icon,nact,nele,ibo,idsym,isym,
     *    iele,ista,iscr,icha)
c     -------------------------------------------------
c
c     Routine to return symmetry for a single spin space function.
c     icon(i) contains orbital occupied by electron i.
c     nact    No. of orbitals.
c     nele    No. of electrons.
c     ibo(i) contains symmetry of orbital i.
c     idsym  specifies the point group.
c     isym   returns the symmetry(irreducible rep) of the icon.
c
c     CONVENTION FOR IBO, IDSYM AND ISYM.
c
c     Point group  idsym  Irred rep  isym  Sym operation used
c   -----------------------------------------------------------
c        Ci          1       Ag       1        i
c                            Au       2
c    
c        Cs          1       A'       1      (sigma)h
c                            A''      2      
c       
c        C2          1       A        1       C2
c                            B        2
c
c        D2          2       A        1       C2(z)
c                            B1       2       C2(y)
c                            B2       3     
c                            B3       4
c
c        C2v         2       A1       1       C2
c                            A2       2       (sigma)v(xz)
c                            B1       3  
c                            B2       4
c
c        C2h         2       Ag       1       i
c                            Bg       2       (sigma)h
c                            Bu       3
c                            Au       4
c
c        D2h         3       Ag       1       (sigma)(xy)
c                            B1g      2       (sigma)(xz)
c                            B2g      3       (sigma)(yz)
c                            B3g      4 
c                            Au       5
c                            B1u      6
c                            B2u      7
c                            B3u      8
c
c
      implicit double precision(a-h,o-z)
      dimension iele(3),ista(3),icha(34)
      dimension iscr(3)
      dimension icon(nele),ibo(nact)
      call getdata(iele,ista,icha)
c      data (iele(i),i=1,3) /2,4,8/
c      data (ista(i),i=1,3) /1,3,11/
c      data (icha(i),i=1,34) /1,-1,1,1,1,-1,-1,1,-1,-1,
c     *   1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,
c     *   1,-1,1,1,1,-1/
c
      ist = ista(idsym)
      iel = iele(idsym)
      call sym(icon,nact,nele,ibo,icha(ist),idsym,iel,isym)
      return
      end
c
c     ----------------------------------------------------
      subroutine sym(icon,nact,nele,ibo,icha,idi,iel,isym)
c     ----------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension ibo(nact),icon(nele)
      dimension iscr(3)
      dimension icha(idi,iel)
c
      do 7 kk=1,idi
	 iscr(kk) = 1
    7 continue
      do 13 ii=1,nele
	 ia = icon(ii)
	 do 20 jj=1,idi
	    iscr(jj) = iscr(jj)*icha(jj,ibo(ia))
   20    continue
   13 continue
c
      do 56 ii=1,iel
	 do 89 jj=1,idi
	    if (iscr(jj).ne.icha(jj,ii)) goto 56
   89    continue
	 isym = ii
	 return
   56 continue
c
      write(6,*) 'Element not identified'
      return
      end
c     ----------------------------------
      subroutine getdata(iele,ista,icha)
c     ----------------------------------
      implicit double precision(a-h,o-z)
      dimension iele(3),ista(3),icha(34)
      iele(1) = 2
      iele(2) = 4
      iele(3) = 8
c
      ista(1) = 1
      ista(2) = 3
      ista(3) = 11
c
      icha(1) = 1
      icha(2) = -1
      icha(3) = 1
      icha(4) = 1
      icha(5) = 1
      icha(6) = -1
      icha(7) = -1
      icha(8) = 1
      icha(9) = -1
      icha(10) = -1
      icha(11) = 1
      icha(12) = 1
      icha(13) = 1
      icha(14) = 1
      icha(15) = -1
      icha(16) = -1
      icha(17) = -1
      icha(18) = 1
      icha(19) = -1
      icha(20) = -1
      icha(21) = -1
      icha(22) = 1
      icha(23) = -1
      icha(24) = -1
      icha(25) = -1
      icha(26) = -1
      icha(27) = 1
      icha(28) = 1
      icha(29) = 1
      icha(30) = -1
      icha(31) = 1
      icha(32) = 1
      icha(33) = 1
      icha(34) = -1
      return
      end
c
c     ------------------------
      subroutine prisym(idsym)
c     ------------------------
      implicit double precision(a-h,o-z)
      if (idsym.eq.0) then
	 write(6,*) 'Have specified C1 symmetry'
      endif
      if (idsym.eq.1) then
      write(6,*) ' Definitions of irreducible representations: ',
     *'        1=(1),        2=(-1)'
      endif
      if (idsym.eq.2) then
       write(6,*) ' Definitions of irreducible representations: ',
     *'        1=(11),       2=(1-1)'
         write(6,*) '                                            ',
     *'        3=(-11),      4=(-11)'
      endif
      if (idsym.eq.3) then
      write(6,*) ' Definitions of irreducible representations: ',
     *'        1=(111),      2=(1-1-1)' 
	 write(6,*) '                                            ',
     *'        3=(-11-1),    4=(-1-11)'
	 write(6,*) '                                            ',
     *'        5=(-1-1-1),   6=(-111)'
	 write(6,*) '                                            ',
     *'        7=(1-11),     8=(11-1)'
      endif
      return
      end
c
c     ---------------------------------------------
      subroutine advanc(con,nele,norb)
c     ---------------------------------------------
      implicit double precision(a-h,o-z)
      integer con(*)
c
      if (con(nele).eq.norb) then 
	 do 50 i=nele-1,1,-1
	    if (con(i+1)-con(i).gt.1) then
	       con(i) = con(i) + 1
	       do 40 j=i+1,nele
	          con(j) = con(j-1) + 1 
   40          continue
	       return
            endif
   50    continue
      endif
c
	 con(nele) = con(nele)+1
c
      return
      end
c
cJAKAL








c
c     -------------------------
      subroutine sortq(n,arr,ipica,ipicb)
c     -------------------------
      integer n,M,NSTACK
      integer ia,ib,iat,ibt
      double precision arr(n)
      integer ipica(n),ipicb(n)
      parameter (M=7,NSTACK=50)
c
c         Sorts an array arr(1:n) into ascending numberical order
c         using the Quicksort algorithm.  n is input; arr is
c         replaced on output by its sorted rearrangement.
c         Parameters: M is the size of the subarrays sorted by
c         straight insertion and NSTACK is the required
c         auxiliary storage.
c
      integer i,ir,j,jstack,k,ll,istack(NSTACK)
      double precision a,temp
      jstack=0
      ll=1
      ir=n
    1 if (ir-ll.lt.M) then
         do 12 j=ll+1,ir
            a=arr(j)
            ia=ipica(j)
            ib=ipicb(j)
            do 11 i=j-1,ll,-1
               if (arr(i).le.a) goto 2
               arr(i+1) = arr(i)
               ipica(i+1)=ipica(i)
               ipicb(i+1)=ipicb(i)
   11       continue
            i=ll-1
    2       arr(i+1)=a
            ipica(i+1)=ia
            ipicb(i+1)=ib
   12    continue
         if(jstack.eq.0)return
         ir=istack(jstack)
         ll=istack(jstack-1)
         jstack=jstack-2
      else
         k=(ll+ir)/2
         temp=arr(k)
         iat=ipica(k)
         ibt=ipicb(k)
         arr(k)=arr(ll+1)
         ipica(k)=ipica(ll+1)
         ipicb(k)=ipicb(ll+1)

         arr(ll+1)=temp
         ipica(ll+1)=iat
         ipicb(ll+1)=ibt
         if (arr(ll).gt.arr(ir)) then
            temp=arr(ll)
            iat=ipica(ll)
            ibt=ipicb(ll)
            arr(ll)=arr(ir)
            ipica(ll)=ipica(ir)
            ipicb(ll)=ipicb(ir)
            arr(ir)=temp
            ipica(ir)=iat
            ipicb(ir)=ibt
         endif
         if (arr(ll+1).gt.arr(ir)) then
            temp=arr(ll+1)
            iat=ipica(ll+1)
            ibt=ipicb(ll+1)
            arr(ll+1)=arr(ir)
            ipica(ll+1)=ipica(ir)
            ipicb(ll+1)=ipicb(ir)
            arr(ir)=temp
            ipica(ir)=iat
            ipicb(ir)=ibt
         endif
         if (arr(ll).gt.arr(ll+1)) then
            temp=arr(ll)
            iat=ipica(ll)
            ibt=ipicb(ll)
            arr(ll)=arr(ll+1)
            ipica(ll)=ipica(ll+1)
            ipicb(ll)=ipicb(ll+1)
            arr(ll+1)=temp
            ipica(ll+1)=iat
            ipicb(ll+1)=ibt
         endif
         i=ll+1
         j=ir
         a=arr(ll+1)
         ia=ipica(ll+1)
         ib=ipicb(ll+1)

    3    continue
            i=i+1
         if (arr(i).lt.a) goto 3
    4    continue
            j=j-1
         if (arr(j).gt.a) goto 4
         if (j.lt.i) goto 5
         temp=arr(i)
         iat=ipica(i)
         ibt=ipicb(i)
         arr(i) = arr(j)
         ipica(i)=ipica(j)
         ipicb(i)=ipicb(j)
         arr(j) = temp
         ipica(j)=iat
         ipicb(j)=ibt
         goto 3
    5    arr(ll+1)=arr(j)
         ipica(ll+1)=ipica(j)
         ipicb(ll+1)=ipicb(j)
         arr(j) = a
         ipica(j)=ia
         ipicb(j)=ib
         jstack = jstack+2
         if (jstack.gt.NSTACK) pause 'NSTACK too small in sort'
         if (ir-i+1.ge.j-ll) then
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=ll
            ll=i
         endif
      endif

      goto 1
      end






c
C      -------------
C*     SUBROUTINES:
C      -------------



c     ---------------------------------------------------------
c     ***
      SUBROUTINE RLABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT,ICON,NSPRO)
c
c     ---------------------------------------------------------
c     This subroutine
c     ---------------
c     
c     *
c     * returns  ICON(j=1,NTOT)
c     * given NSPRO (=Number of space product )

c    ************************************
c        RLABELSDTQ:
c      ---------------
c      This code labels SDTQ56-space products in a 
c      systematic way
c      --------------
c

      implicit double precision(a-h,o-z)
c     -------------------------------------------

      parameter (mxorb = 100)
c
      parameter (mxsprod = 2400000)
      parameter (mxsprod1 = 100000)
      parameter (mxCI = 1000000)
      parameter (mxstring = 2000000)
      parameter (mxna = 200)



      DIMENSION ICON(mxorb),NEL(mxorb)




C--------------------------------------------
C*     doubly and singly occupied MOs
C      ------------------------------

       ICOUNT = 0






c      ---------------------------
       IF (IEXCIT.EQ.1) GO TO 33
c      ---------------------------
       IF (IEXCIT.EQ.2) GO TO 44 
c      ---------------------------
       IF (IEXCIT.EQ.3) GO TO 350 
c      ----------------------------
       IF (IEXCIT.EQ.4) GO TO 490 
c      ----------------------------
       IF (IEXCIT.EQ.5) GO TO 1100 
c      ----------------------------
       IF (IEXCIT.EQ.6) GO TO 3000 
c      ----------------------------



       write (6,*) 'subr.RLABELSDTQ:  no excitations '
       STOP




C***** code to generate space-products
C      -------------------------------
c



C      one-determinant (SCF) :
C      -----------------------
       DO 30 I=1,NOCC
       NEL(I) = 2
 30    CONTINUE
       DO 32 I=NOCC+1,NTOT
       NEL(I) = 0
 32    CONTINUE
C*
c***** WRITE (6,*) (NEL(I), I=1,NTOT) 
c***** ICOUNT = ICOUNT + 1
C*     -------------------





  33   A=0
c      ***
c***** 'single excitations '
C      -------------------


       DO 35 I=1,NOCC
c      --------------
       DO 36 J=1,NOCC
       NEL(J)= 2
  36   CONTINUE

       NSING = NOCC+1-I
       NEL(NSING) = 1



c      ^^^^^^^^^^^^^^^^^^^^
c      virtual
c      ^^^^^^^^^^^^^^^^^^^^
       DO 39 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 39    CONTINUE


       DO 40 MM=NOCC+1,NTOT
c      --------------------
       NEL(MM) = 1 
c
c***** WRITE (6,*) (NEL(IX), IX=1,NTOT) 
c      --------------------
       ICOUNT = ICOUNT + 1
C*     --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 43
       GO TO 5000
  43   A=0
c      *****************

       NEL(MM) = 0
 40    CONTINUE
c      -------------
 35    CONTINUE


c      ----------------------------------------------
c****  write (6,*) 'subr. LABELSDTQ: singlet got here-all'
       GO TO 5000
c      ----------








  44   A=0
c      ***
c***** 'DIAGONAL-double excitations'

       DO 45 I=1,NOCC
c      -------------------

       DO 47 J=1,NOCC
       NEL(J)= 2
  47   CONTINUE

       NSING = NOCC+1-I
       NEL(NSING) = 0 



c      ^^^^^^^^^^^^^^^^^^^^
c      virtual
c      ^^^^^^^^^^^^^^^^^^^^
       DO 49 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 49    CONTINUE

       DO 50 MM=NOCC+1,NTOT
c      --------------------
       NEL(MM) = 2 
c
c***** WRITE (6,*) (NEL(L), L=1,NTOT) 
c      -------------------
       ICOUNT = ICOUNT + 1
C*     -------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 51 
       GO TO 5000
  51   A=0
c      *****************

       NEL(MM) = 0
 50    CONTINUE




c      ^^^^^^^^^^^^^^^^^
c      virtual
c      ^^^^^^^^^^^^^^^^^
       DO 55 JJJ=1,NVIRT
c      ------------------
       DO 56 KKK=NOCC+1,NTOT
       NEL(KKK) = 0
  56   CONTINUE
       M = NTOT+1-JJJ
       NEL(M) = 1



       DO 60 II=1,M-NOCC-1
c      -------------------
       MMM = M - II
c      ------------
       NEL(MMM) = 1

c
c****  WRITE (6,*) (NEL(L), L=1,NTOT) 
c      -------------------
       ICOUNT = ICOUNT + 1
C*     -------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 61 
       GO TO 5000
  61   A=0
c      *****************

       NEL(M-II) = 0
 60    CONTINUE

       NEL(M) = 0
 55    CONTINUE


c      --------------
 45    CONTINUE



c      -------------------------
c***** 'mixed double excitations '
C      -------------------------

       DO 275 I=1,NOCC
C      ---------------------

       DO 277 J=1,NOCC
       NEL(J)= 2
 277   CONTINUE

       M1 = NOCC+1-I 
c      -------------
       NEL(M1) = 1 


       DO 278 K1=1,M1-1
C      ----------------
       M2 = M1 - K1
c      *
       NEL(M2)= 1





c      ^^^^^^^^^^^^^^^^^^^^^
c      virtual
c      ^^^^^^^^^^^^^^^^^^^^^
       DO 280 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 280   CONTINUE

       DO 285 MM=NOCC+1,NTOT
c      ---------------------
       NEL(MM) = 2 
c
c***** WRITE (6,*) (NEL(L), L=1,NTOT) 
C*     --------------------
       ICOUNT = ICOUNT + 1
C*     --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 281 
       GO TO 5000
 281   A=0
c      *****************


       NEL(MM) = 0
c      --------
 285   CONTINUE





c      ^^^^^^^^^
c      virtual
c      ^^^^^^^^^
       DO 300 JJJ=1,NVIRT
C      ------------------

       DO 301 KKK=NOCC+1,NTOT
       NEL(KKK) = 0
 301   CONTINUE
       M = NTOT+1-JJJ
c      ------------------------
       NEL(M) = 1


       DO 305 II=1,M-NOCC-1
c      --------------------
       MM = M - II
c      -----------
       NEL(MM) = 1

c
c***** WRITE (6,*) (NEL(L), L=1,NTOT) 
c      -------------------
       ICOUNT = ICOUNT + 1
C*     -------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 303
       GO TO 5000
 303   A=0
c      *****************


       NEL(MM) = 0
 305   CONTINUE

       NEL(M) = 0
 300   CONTINUE


c      ***************
c      SCF-MOs
c      ---------------
       NEL(M1-K1) = 2
 278   CONTINUE


       NEL(M1) = 2
 275   CONTINUE

c      ----------------------------------------------
c      ----------------------------------------------
c***** write (6,*) 'subr. LABELSDTQ: doubles got here-all'
       GO TO 5000
c      ----------






 350   A=0
c      ******
c*     'triple excitations start here'
C      -------------------

C****** 'T-case:[1,0] 
c      -----------------
       DO 362 I= 1,NOCC
       NEL(I) = 2
 362   CONTINUE


       DO 363 KKZ=1,NOCC
c      -------------------
       NXZ = NOCC + 1 - KKZ
c      *
       NEL(NXZ) = 0


       DO 367 KK1Z=1,NOCC-1
c      ----------------------
       M1Z = NOCC + 1 - KK1Z 
c      *
       IF (M1Z.GT.NXZ) GO TO 371 
       IF (M1Z.LE.NXZ) GO TO 372 
 371   AAA=0
       L1Z = M1Z
       GO TO 373 
 372   AAA=0
       L1Z = M1Z - 1
 373   AAA=0
c--------------------
       NEL(L1Z) = 1



c      ^^^^^^^^^^^^^^^^
c      Virtual MOs
c      ^^^^^^^^^^^^^^^^
c      ****************
c      CASE-1 [2,1]
c      ****************
       DO 380 K = NOCC+1,NTOT
       NEL(K) = 0
 380   CONTINUE


       DO 385 KK=1,NVIRT
c      -------------------
       NX = NTOT + 1 - KK
c      ------------------------
       NEL(NX) = 2



       DO 389 KK1=1,NVIRT-1
c      ----------------------
       M1 = NTOT + 1 - KK1 
c      *
       IF (M1.GT.NX) GO TO 390 
       IF (M1.LE.NX) GO TO 391
 390   AAA=0
       L1 = M1
       GO TO 392
 391   AAA=0
       L1 = M1 - 1
 392   AAA=0
c      ------------------------
       NEL(L1) = 1

c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c   *  write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 400 
       GO TO 5000
 400   A=0
c      *************


 395   AAA=0.0
       NEL(L1) = 0
c      ---------
 389   CONTINUE


 388   AAA=0.0
       NEL(NX) = 0
c      ---------
 385   CONTINUE





c      ^^^^^^^^^^^^^^^^^
c      virtual MOs
c      ^^^^^^^^^^^^^^^^^
c      *****************
c      CASE-2) [1,1,1]
c      *****************

       DO 397 K = NOCC+1,NTOT
       NEL(K) = 0
 397   CONTINUE


       DO 398 KK = 1,NVIRT-2
c      ---------------------
       M1 = NTOT+1-KK
c      ------------------------
       NEL(M1) = 1


       DO 402 KK1 = 1,M1-NOCC-2
c      -------------------------
       M2 = M1 - KK1
c      ------------------------
       NEL(M2) = 1


       DO 406 KK2 = 1,M2-NOCC-1
c      -------------------------
       M3 = M2 - KK2
c      ------------------------
       NEL(M3) = 1


c      ----------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 410
       GO TO 5000
 410   A=0
c      *************


 409   AAA=0.0
       NEL(M3) = 0
c      ----------
 406   CONTINUE


 405   AAA=0.0
       NEL(M2) = 0
c      ----------
 402   CONTINUE


 401   AAA=0.0
       NEL(M1) = 0
c      ----------
 398   CONTINUE


 381   A=0.0
c      ---------------


c      ---------------------
c      occupied MOs: CASE-A)
c      ---------------------
 412   AAA=0.0 
       NEL(L1Z) = 2
 367   CONTINUE


 411   AAA=0.0
       NEL(NXZ) = 2
 363   CONTINUE
c      -------------
c      *******
 416   AAA=0.0
c      *******





c*     write (6,*) 'T-case: [1,1,1]'
c      ********************************
c      occupied MOs CASE-B) [1,1,1]
c      ********************************
       DO 414 I= 1,NOCC
       NEL(I) = 2
 414   CONTINUE


       DO 415 KKW = 1,NOCC-2
c      ---------------------
       M1W = NOCC+1-KKW
c      --------------------------
       NEL(M1W) = 1



       DO 418 KK1W = 1,M1W-2
c      -------------------------
       M2W = M1W - KK1W
c      -------------
       NEL(M2W) = 1



       DO 421 KK2W = 1,M2W-1
c      -------------------------
       M3W = M2W - KK2W
c      --------------------------
       NEL(M3W) = 1





c      ^^^^^^^^^^^^^^^^
c      Virtual MOs
c      ^^^^^^^^^^^^^^^^
c      ****************
c      CASE-1 [2,1]
c      ****************
       DO 424 K = NOCC+1,NTOT
       NEL(K) = 0
 424   CONTINUE


       DO 425 KK=1,NVIRT
c      -------------------
       NX = NTOT + 1 - KK
c      -----------
       NEL(NX) = 2



       DO 430 KK1=1,NVIRT-1
c      ----------------------
       M1 = NTOT + 1 - KK1 
c
       IF (M1.GT.NX) GO TO 431 
       IF (M1.LE.NX) GO TO 432
 431   AAA=0
       L1 = M1
       GO TO 433
 432   AAA=0
       L1 = M1 - 1
 433   AAA=0
c      -----------
       NEL(L1) = 1


c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 440
       GO TO 5000
 440   A=0
c      *************


 439   AAA=0.0
       NEL(L1) = 0
c      ---------
 430   CONTINUE


 438   AAA=0.0
       NEL(NX) = 0
c      ---------
 425   CONTINUE






c      ^^^^^^^^^^^^^^^^^
c      virtual MOs
c      ^^^^^^^^^^^^^^^^^
c      *****************
c      CASE-2) [1,1,1]
c      *****************
       DO 441 K = NOCC+1,NTOT
       NEL(K) = 0
 441   CONTINUE


       DO 442 KK = 1,NVIRT-2
c      ---------------------
       M1 = NTOT+1-KK
c      -------------
       NEL(M1) = 1


       DO 446 KK1 = 1,M1-NOCC-2
c      -------------------------
       M2 = M1 - KK1
c      ------------
       NEL(M2) = 1



       DO 450 KK2 = 1,M2-NOCC-1
c      -------------------------
       M3 = M2 - KK2
c      -----------
       NEL(M3) = 1



c      ----------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 451 
       GO TO 5000
 451   A=0
c      *************

 453   AAA=0.0
       NEL(M3) = 0
c      ----------
 450   CONTINUE


 449   AAA=0.0
       NEL(M2) = 0
c      ----------
 446   CONTINUE


 445   AAA=0.0
       NEL(M1) = 0
c      ----------
 442   CONTINUE


c      -------------------------------
c      occupied MOs CASE-B) [1,1,1]
c      -------------------------------
 457   AAA=0.0
       NEL(M3W) = 2 
c      ----------
 421   CONTINUE


 456   AAA=0.0
       NEL(M2W) = 2 
c      ----------
 418   CONTINUE


 455   AAA=0.0
       NEL(M1W) = 2 
c      ----------
 415   CONTINUE

c      ----------------------------------------------
c      ----------------------------------------------
c***** write (6,*) 'subr. LABELSDTQ: triples got here-all'
       GO TO 5000
c      ----------









  490  AAA=0.0
c*     'quadruple excitations start here'
c      **********************************


c*     write (6,*) 'Q-case: [0,0]'
c      ---------------------------

       DO 503 I= 1,NOCC
       NEL(I) = 2
 503   CONTINUE


       DO 505 I=1,NOCC-1
c      ------------------
       M = NOCC+1-I
c      ------------
       NEL(M) = 0


       DO 510 J=1,M-1
c      ------------------
       MM = M - J
c      -----------
       NEL(MM) = 0




c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *************
c      CASE-1) [2,2]
c      *************

       DO 520 K = NOCC+1,NTOT
       NEL(K) = 0
 520   CONTINUE


       DO 525 KK=1,NVIRT-1
c      -------------------
       MK = NTOT+1-KK
c      -----------
       NEL(MK) = 2


       DO 530 KKK=1,NVIRT-KK
c      ----------------------
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2


c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 535 
       GO TO 5000
 535   A=0
c      *************

 533   AAA=0.0
       NEL(MMK) = 0
 530   CONTINUE


 528   AAA=0.0
       NEL(MK) = 0
 525   CONTINUE




c      ^^^^^^^^^^^^^^^^
c      Virtual MOs
c      ^^^^^^^^^^^^^^^^
c      ****************
c      CASE-2 [2,1,1]
c      ****************

       DO 540 K = NOCC+1,NTOT
       NEL(K) = 0
 540   CONTINUE


       DO 545 KK=1,NVIRT
c      -------------------
       NX = NTOT + 1 - KK
c      -----------
       NEL(NX) = 2



       DO 550 KK1=1,NVIRT-2
c      ----------------------
       M1 = NTOT + 1 - KK1 
c      ----------------------
       IF (M1.GT.NX) GO TO 551 
       IF (M1.LE.NX) GO TO 552
 551   AAA=0
       L1 = M1
       GO TO 553
 552   AAA=0
       L1 = M1 - 1
 553   AAA=0
c      -----------
       NEL(L1) = 1



       DO 560 KK2=1,M1-NOCC-2    
c      ------------------------
       M2 = M1 - KK2
c      ------------------------
       IF (M2.GT.NX) GO TO 561 
       IF (M2.LE.NX) GO TO 562 
 561   AAA=0
       L2 = M2
       GO TO 563
 562   AAA=0
       L2 = M2 - 1
 563   AAA=0
c      -----------
       NEL(L2) = 1


c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 565
       GO TO 5000
 565   A=0
c      *************

 566   AAA=0.0
       NEL(L2) = 0
c      ---------
 560   CONTINUE



 556   AAA=0.0
       NEL(L1) = 0
c      ---------
 550   CONTINUE



 548   AAA=0.0
       NEL(NX) = 0
c      ---------
 545   CONTINUE





c      ^^^^^^^^^^^^^^^^^
c      virtual MOs
c      ^^^^^^^^^^^^^^^^^
c      *****************
c      CASE-3) [1,1,1,1]
c      *****************
       DO 569 K = NOCC+1,NTOT
       NEL(K) = 0
 569   CONTINUE


       DO 570 KK = 1,NVIRT-3
c      ---------------------
       M1 = NTOT+1-KK
c      -----------
       NEL(M1) = 1



       DO 575 KK1 = 1,M1-NOCC-3
c      -------------------------
       M2 = M1 - KK1
c      -----------
       NEL(M2) = 1



       DO 580 KK2 = 1,M2-NOCC-2
c      -------------------------
       M3 = M2 - KK2
c      -----------
       NEL(M3) = 1



       DO 585 KK3 = 1,M3-NOCC-1
c      -------------------------
       M4 = M3 - KK3
c      -----------
       NEL(M4) = 1


c      ----------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 590
       GO TO 5000
 590   A=0
c      *************

 588   AAA=0.0
       NEL(M4) = 0
c      ----------
 585   CONTINUE


 583   AAA=0.0
       NEL(M3) = 0
c      ----------
 580   CONTINUE


 578   AAA=0.0
       NEL(M2) = 0
c      ----------
 575   CONTINUE


 573   AAA=0.0
       NEL(M1) = 0
c      ----------
 570   CONTINUE


c      -----------------------------
c      occupied MOs: CASE-A) [0,0]
c      -----------------------------
 592   AAA=0.0
       NEL(MM) = 2
 510   CONTINUE

 591   AAA=0.0 
       NEL(M) = 2
 505   CONTINUE







 595   AAA=0.0
C*     write (6,*) 'Q-case: [0,1,1]'
c      *****************************
c      occupied MOs: CASE-B) [0,1,1]
c      *****************************
       DO 605 I= 1,NOCC
       NEL(I) = 2
 605   CONTINUE


       DO 620 KKZ=1,NOCC
c      -------------------
       NXZ = NOCC + 1 - KKZ
c      ------------
       NEL(NXZ) = 0



       DO 627 KK1Z=1,NOCC-2
c      ----------------------
       M1Z = NOCC + 1 - KK1Z 
c      ----------------------
       IF (M1Z.GT.NXZ) GO TO 631 
       IF (M1Z.LE.NXZ) GO TO 632 
 631   AAA=0
       L1Z = M1Z
       GO TO 633 
 632   AAA=0
       L1Z = M1Z - 1
 633   AAA=0
c      ------------
       NEL(L1Z) = 1



       DO 640 KK2Z=1,M1Z-2    
c      ------------------------
       M2Z = M1Z - KK2Z

       IF (M2Z.GT.NXZ) GO TO 641 
       IF (M2Z.LE.NXZ) GO TO 642 
 641   AAA=0
       L2Z = M2Z
       GO TO 643
 642   AAA=0
       L2Z = M2Z - 1
 643   AAA=0
c      --------------
       NEL(L2Z) = 1



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *************
c      CASE-1) [2,2]
c      *************
c      -----------------------
       DO 650 K = NOCC+1,NTOT
       NEL(K) = 0
 650   CONTINUE


       DO 655 KK=1,NVIRT-1
c      -------------------
       MK = NTOT+1-KK
c      -----------
       NEL(MK) = 2


       DO 660 KKK=1,NVIRT-KK
c      ----------------------
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2


c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c*
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 665 
       GO TO 5000
 665   A=0
c      *************
c      
 663   AAA=0.0
       NEL(MMK) = 0
 660   CONTINUE


 658   AAA=0.0
       NEL(MK) = 0
 655   CONTINUE





c      ^^^^^^^^^^^^^^^^
c      Virtual MOs
c      ^^^^^^^^^^^^^^^^
c      ****************
c      CASE-2 [2,1,1]
c      ****************
       DO 670 K = NOCC+1,NTOT
       NEL(K) = 0
 670   CONTINUE


       DO 675 KK=1,NVIRT
c      -------------------
       NX = NTOT + 1 - KK
c      -----------
       NEL(NX) = 2


       DO 680 KK1=1,NVIRT-2
c      ----------------------
       M1 = NTOT + 1 - KK1 
c      ----------------------
       IF (M1.GT.NX) GO TO 681 
       IF (M1.LE.NX) GO TO 682
 681   AAA=0
       L1 = M1
       GO TO 683
 682   AAA=0
       L1 = M1 - 1
 683   AAA=0
c      -----------
       NEL(L1) = 1



       DO 690 KK2=1,M1-NOCC-2    
c      ------------------------
       M2 = M1 - KK2
c      ------------------------
       IF (M2.GT.NX) GO TO 691 
       IF (M2.LE.NX) GO TO 692 
 691   AAA=0
       L2 = M2
       GO TO 693
 692   AAA=0
       L2 = M2 - 1
 693   AAA=0
c      -----------
       NEL(L2) = 1

c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 698
       GO TO 5000
 698   A=0
c      *************

 696   AAA=0.0
       NEL(L2) = 0
c      ---------
 690   CONTINUE


 686   AAA=0.0
       NEL(L1) = 0
c      ---------
 680   CONTINUE


 678   AAA=0.0
       NEL(NX) = 0
c      ---------
 675   CONTINUE





c      ^^^^^^^^^^^^^^^^^
c      virtual MOs
c      ^^^^^^^^^^^^^^^^^
c      *****************
c      CASE-3) [1,1,1,1]
c      *****************

       DO 699 K = NOCC+1,NTOT
       NEL(K) = 0
 699   CONTINUE


       DO 700 KK = 1,NVIRT-3
c      ---------------------
       M1 = NTOT+1-KK
c      -----------
       NEL(M1) = 1


       DO 705 KK1 = 1,M1-NOCC-3
c      -------------------------
       M2 = M1 - KK1
c      -----------
       NEL(M2) = 1


       DO 710 KK2 = 1,M2-NOCC-2
c      -------------------------
       M3 = M2 - KK2
c      -----------
       NEL(M3) = 1


       DO 715 KK3 = 1,M3-NOCC-1
c      -------------------------
       M4 = M3 - KK3
c      -----------
       NEL(M4) = 1

c      ----------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 719 
       GO TO 5000
 719   A=0
c      *************

 718   AAA=0.0
       NEL(M4) = 0
c      ----------
 715   CONTINUE


 713   AAA=0.0
       NEL(M3) = 0
c      ----------
 710   CONTINUE


 708   AAA=0.0
       NEL(M2) = 0
c      ----------
 705   CONTINUE


 703   AAA=0.0
       NEL(M1) = 0
c      ----------
 700   CONTINUE

c      ---------------------
c      occupied MOs: CASE-B)
c      ---------------------
 722   AAA=0.0
       NEL(L2Z) = 2
 640   CONTINUE


 721   AAA=0.0
       NEL(L1Z) = 2
 627   CONTINUE


 720   AAA=0.0
       NEL(NXZ) = 2
 620   CONTINUE






 750   AAA=0
c*     write (6,*) 'Q-case: [1,1,1,1]'
c      ********************************
c      occupied MOs CASE-C) [1,1,1,1]
c      ********************************
       DO 803 I= 1,NOCC
       NEL(I) = 2
 803   CONTINUE


       DO 805 KKW = 1,NOCC-3
c      ---------------------
       M1W = NOCC+1-KKW
c      -------------
       NEL(M1W) = 1


       DO 810 KK1W = 1,M1W-3
c      -------------------------
       M2W = M1W - KK1W
c      ------------
       NEL(M2W) = 1


       DO 815 KK2W = 1,M2W-2
c      -------------------------
       M3W = M2W - KK2W
c      ------------
       NEL(M3W) = 1


       DO 820 KK3W = 1,M3W-1
c      -------------------------
       M4W = M3W - KK3W
c      ------------
       NEL(M4W) = 1



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *************
c      CASE-1) [2,2]
c      *************
       DO 850 K = NOCC+1,NTOT
       NEL(K) = 0
 850   CONTINUE


       DO 855 KK=1,NVIRT-1
c      -------------------
       MK = NTOT+1-KK
c      -----------
       NEL(MK) = 2


       DO 860 KKK=1,NVIRT-KK
c      ----------------------
       MMK = MK - KKK
c      -------------
       NEL(MMK) = 2


c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 865 
       GO TO 5000
 865   A=0
c      *************

 863   AAA=0.0
       NEL(MMK) = 0
 860   CONTINUE


 858   AAA=0.0
       NEL(MK) = 0
 855   CONTINUE



c      ^^^^^^^^^^^^^^^^
c      Virtual MOs
c      ^^^^^^^^^^^^^^^^
c      ****************
c      CASE-2 [2,1,1]
c      ****************
       DO 870 K = NOCC+1,NTOT
       NEL(K) = 0
 870   CONTINUE


       DO 875 KK=1,NVIRT
c      -------------------
       NX = NTOT + 1 - KK
c      -----------
       NEL(NX) = 2


       DO 880 KK1=1,NVIRT-2
c      ----------------------
       M1 = NTOT + 1 - KK1 
c      ----------------------
       IF (M1.GT.NX) GO TO 881 
       IF (M1.LE.NX) GO TO 882
 881   AAA=0
       L1 = M1
       GO TO 883
 882   AAA=0
       L1 = M1 - 1
 883   AAA=0
c      -----------
       NEL(L1) = 1


       DO 890 KK2=1,M1-NOCC-2    
c      ------------------------
       M2 = M1 - KK2
c      ------------------------
       IF (M2.GT.NX) GO TO 891 
       IF (M2.LE.NX) GO TO 892 
 891   AAA=0
       L2 = M2
       GO TO 893
 892   AAA=0
       L2 = M2 - 1
 893   AAA=0
c      -----------
       NEL(L2) = 1

c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 898
       GO TO 5000
 898   A=0
c      *************

 896   AAA=0.0
       NEL(L2) = 0
c      ---------
 890   CONTINUE


 886   AAA=0.0
       NEL(L1) = 0
c      ---------
 880   CONTINUE


 878   AAA=0.0
       NEL(NX) = 0
c      ---------
 875   CONTINUE




c      ^^^^^^^^^^^^^^^^^
c      virtual MOs
c      ^^^^^^^^^^^^^^^^^
c      *****************
c      CASE-3) [1,1,1,1]
c      *****************
       DO 899 K = NOCC+1,NTOT
       NEL(K) = 0
 899   CONTINUE


       DO 900 KK = 1,NVIRT-3
c      ---------------------
       M1 = NTOT+1-KK
c      -----------
       NEL(M1) = 1


       DO 905 KK1 = 1,M1-NOCC-3
c      -------------------------
       M2 = M1 - KK1
c      -----------
       NEL(M2) = 1


       DO 910 KK2 = 1,M2-NOCC-2
c      -------------------------
       M3 = M2 - KK2
c      -----------
       NEL(M3) = 1


       DO 915 KK3 = 1,M3-NOCC-1
c      -------------------------
       M4 = M3 - KK3
c      -----------
       NEL(M4) = 1


c      ----------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       IF (ICOUNT.NE.NSPRO) GO TO 919 
       GO TO 5000
 919   A=0
c      *************

 918   AAA=0.0
       NEL(M4) = 0
c      ----------
 915   CONTINUE


 913   AAA=0.0
       NEL(M3) = 0
c      ----------
 910   CONTINUE


 908   AAA=0.0
       NEL(M2) = 0
c      ----------
 905   CONTINUE


 903   AAA=0.0
       NEL(M1) = 0
c      ----------
 900   CONTINUE


c      -------------------------------
c      occupied MOs CASE-C) [1,1,1,1]
c      -------------------------------
 923   AAA=0.0
       NEL(M4W) = 2 
c      ----------
 820   CONTINUE


 922   AAA=0.0
       NEL(M3W) = 2 
c      ----------
 815   CONTINUE


 921   AAA=0.0
       NEL(M2W) = 2 
c      ----------
 810   CONTINUE


 920   AAA=0.0
       NEL(M1W) = 2 
c      ----------
 805   CONTINUE
c      ----------------------------------------------
c      ----------------------------------------------
c***** write (6,*) 'subr. LABELSDTQ: Q- got here-all'
       GO TO 5000
c      ----------





c*****LABELSDTQ: SDTQ56-excitation-labeling starts ***
c     --------------------------------------------


c      ---------------------------------------
c      subr. LABELSDTQ:  QUINTUPLE EXCITATIONS
c      ---------------------------------------
c***** write (6,*) 'quintuple(5) exc. start here'
c***** write (6,*) '5-case: [1,0,0] '
c

c      ***
 1100  A=0
c      ***


c      ****************************
c      occ. MO space: CASE [1,0,0]
c      ****************************
c*
c      SCF-SPACE:
c***** write (6,*) 'Quintuple OCCU[1,0,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,0,0]
c      ***********************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC-1 
       DO 1500 I1=1,MAX1
c      ------------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 0



       MAX2 = N1-1
       DO 1505 I2=1,MAX2
c      -------------------
c      *
       N2 = N1 - I2
       NEL(N2) = 0



       MAX3 = NOCC-2
       DO 1510 I3=1,MAX3
c      ------------------
c      *
       NN3 = NOCC + 1 - I3
       IF (I3.EQ.1) GO TO 1511
       NN3 = LL3 - 1
 1511  A=0
       call SHIFT2(N1,N2,NN3,LL3)
       IF (LL3.EQ.0) GO TO 1512 
       NEL(LL3) = 1



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c
c      write (6,*) 'VIRT [2,2,1] starts'
c      *********************
c      CASE) VIRT [2,2,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 


       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 


       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1535
       M1 = L1 - 1
 1535  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1538 
       NEL(L1) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1540 
       GO TO 5000
 1540  A=0
c      *****************
c      ----------------------------------
       NEL(L1) = 0
       ENDDO
 1538  A=0
c      ***

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'




c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2


       LIM1 = NVIRT-3
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1545
       M1 = L1 - 1
 1545  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1546
       M2 = L2 - 1
 1546  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1547
       M3 = L3 - 1
 1547  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1550
       NEL(L3) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1549
       GO TO 5000
 1549  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L3)=0
       ENDDO
 1550  A=0
c      ****

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1555
       NEL(M5) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1556
       GO TO 5000
 1556  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M5)=0
       ENDDO
c      *****
 1555  A=0

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
C*     ---------------------------------------------


c      -----------------------------
c      occupied MOs: CASE [1,0,0]
c      -----------------------------
       NEL(LL3) = 2 
 1510  CONTINUE
 1512  A=0
c      ***

       NEL(N2) = 2 
 1505  CONTINUE

       NEL(N1) = 2 
 1500  CONTINUE
c      ********





c***** write (6,*) 'Quintuple OCCU[1,0,0] done'
c***** write (6,*) '---------------------------'




c      *****
 1600  A=0
c      *****



cLAIMIS
c      SCF-SPACE:
c***** write (6,*) 'Quintuple OCCU[1,1,1,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,1,1,0]
c      ***********************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC
       DO 1605 I1=1,NOCC
c      -----------------
c      *
       N1 = NOCC + 1 - I1
       NEL(N1) = 0


       MAX2 = NOCC-3
       DO 1610 I2=1,MAX2
c      ------------------
c      *
       NN2 = NOCC + 1 - I2
       IF (I2.EQ.1) GO TO 1611
       NN2 = LL2 - 1
 1611  A=0
       call SHIFT1(N1,NN2,LL2)
       NEL(LL2) = 1



       MAX3 = NN2 - 2 
       DO 1615 I3=1,MAX3
c      -----------------
c      *
       NN3 = LL2 - I3
       IF (I3.EQ.1) GO TO 1616
       NN3 = LL3 - 1
 1616  A=0
       call SHIFT1(N1,NN3,LL3)
       NEL(LL3) = 1



       MAX4 = NN3 - 1 
       DO 1620 I4=1,MAX4
c      -----------------
c      *
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 1621
       NN4 = LL4 - 1
 1621  A=0
       call SHIFT1(N1,NN4,LL4)
       IF (LL4.EQ.0) GO TO 1622 
       NEL(LL4) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1] starts'
c      *********************
c      CASE) VIRT [2,2,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 


       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 


       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1635
       M1 = L1 - 1
 1635  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1638 
       NEL(L1) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1640 
       GO TO 5000
 1640  A=0
c      *****************
c      ----------------------------------
c      ----------------------------------
       NEL(L1) = 0
       ENDDO
 1638  A=0
c      ***

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'




c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2


       LIM1 = NVIRT-3
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1645
       M1 = L1 - 1
 1645  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1646
       M2 = L2 - 1
 1646  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1647
       M3 = L3 - 1
 1647  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1650
       NEL(L3) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1651
       GO TO 5000
 1651  A=0
c      *****************
c      ----------------------------------
       NEL(L3)=0
       ENDDO
 1650  A=0
c      ****

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1655
       NEL(M5) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1656
       GO TO 5000
 1656  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M5)=0
       ENDDO
c      *****
 1655  A=0

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
C*     ---------------------------------------------

       NEL(LL4)=2
 1620  CONTINUE
 1622  A=0
c      ***

       NEL(LL3)=2
 1615  CONTINUE


       NEL(LL2)=2
 1610  CONTINUE


       NEL(N1)=2
 1605  CONTINUE



c***** write (6,*) 'Quintuple OCCU=[1,1,1,0] done'
c***** write (6,*) '---------------------------'







c      *****
 1700  A=0
c      *****




cLAIMIS
c      Quintuple OCCU=[1,1,1,1,1] starts
c      *********************************
c      SCF-SPACE:
c***** write (6,*) 'Quintuple OCCU=[1,1,1,1,1] starts'
c      **************************************
C      occupied MO space: CASE  [1,1,1,1,1]
c      **************************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       MAX1 = NOCC-4
       DO 1805 I1=1,MAX1
c      -----------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 1



       MAX2 = N1 - 4 
       DO 1810 I2=1,MAX2
c      -----------------
c      *
       N2 = N1 - I2
       NEL(N2) = 1


       MAX3 = N2 - 3 
       DO 1815 I3=1,MAX3
c      ------------------
c      *
       N3 = N2 - I3
       NEL(N3) = 1


       MAX4 = N3 - 2 
       DO 1820 I4=1,MAX4
c      -----------------
c      *
       N4 = N3 - I4
       NEL(N4) = 1


       MAX5 = N4 - 1 
       DO 1825 I5=1,MAX5
c      -----------------
c      *
       N5 = N4 - I5
       IF (N5.GT.0) GO TO 1831 
       write (6,*) 'N5=',N5
       STOP
 1831  A=0
       NEL(N5) = 1



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1] starts'
c      *********************
c      CASE) VIRT [2,2,1]
c      *********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 


       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 


       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1835
       M1 = L1 - 1
 1835  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1838 
       NEL(L1) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1840
       GO TO 5000
 1840  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L1) = 0
       ENDDO
 1838  A=0
c      ***

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'




c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2


       LIM1 = NVIRT-3
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1845
       M1 = L1 - 1
 1845  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1846
       M2 = L2 - 1
 1846  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1847
       M3 = L3 - 1
 1847  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1850
       NEL(L3) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1851
       GO TO 5000
 1851  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L3)=0
       ENDDO
 1850  A=0
c      ****

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1855
       NEL(M5) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 1856
       GO TO 5000
 1856  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M5)=0
       ENDDO
c      *****
 1855  A=0

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
C*     ---------------------------------------------

       NEL(N5)=2
 1825  CONTINUE


       NEL(N4)=2
 1820  CONTINUE


       NEL(N3)=2
 1815  CONTINUE


       NEL(N2)=2
 1810  CONTINUE


       NEL(N1)=2
 1805  CONTINUE


c---------------------------------------------------
c***** write (6,*) 'Quintuple OCCU=[1,1,1,1,1] done'
c                   LABELSDTQ:
c***** write (6,*) '---------------------------'
c      ----------------------------------------------
c***** write (6,*) 'subr. LABELSDTQ: 5- got here-all'
       GO TO 5000
c      ----------












c      **********
 3000  A=0
c      **********


c     ----------------------
c***** write (6,*) 'sixtuple exc'
c***** write (6,*) '6-case: [0,0,0]'
c     ----------------------
c*
c*     call CHECKMATE(NOCC,NSCF,MATE,M,MM,MMM,
c*   * MK,MMK,MMK2,ISAY)
c*     IF (ISAY.EQ.0) GO TO 3200 
c      write (6,*) M,MM,MMM,MK,MMK,MMK2
c      write (6,*) ISAY
c
c
c     ----------------------
c     SIXTUPLE-EXCITATIONS:
c     ----------------------
c


c      *******
 3100  A=0
c      *******


c      *********************************
C      occupied MO space: CASE  [0,0,0]
c      *********************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO

       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       DO 3105 I=1,NOCC-2
c      ------------------
       M = NOCC+1-I
c      ------------
       NEL(M) = 0



       DO 3110 J=1,M-2
c      ------------------
       MM = M - J
c      ------------
       NEL(MM) = 0



       DO 3115 JJ=1,MM-1
c      ------------------
       MMM = MM - JJ
c      -------------
       NEL(MMM) = 0



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      --------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3150 
       GO TO 5000
 3150  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3235
       M1 = L1 - 1
 3235  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3236
       LIM4 = LIM4 - 1
 3236  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3237
       M2 = L2 - 1
 3237  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3238 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3240 
       GO TO 5000
 3240  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L2) = 0
       ENDDO
 3238  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3245
       M1 = L1 - 1
 3245  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3246
       M2 = L2 - 1
 3246  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3247
       M3 = L3 - 1
 3247  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3248
       M4 = L4 - 1
 3248  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3250
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3251 
       GO TO 5000
 3251  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3250  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3255
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3256
       GO TO 5000
 3256  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3255  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
C*     -----------------------------------
c      **********************VIRTUAL

c      -----------------------------
c      occupied MOs: CASE [0,0,0]
c      -----------------------------
       NEL(MMM) = 2
 3115  CONTINUE


       NEL(MM) = 2
 3110  CONTINUE


       NEL(M) = 2
 3105  CONTINUE
c      ********

c---------------------------------------------------
c***** write (6,*) 'Sixtuple OCCU=[0,0,0] done'
c***** write (6,*) '---------------------------'







c      *******
 3400  A=0
c      *******




c      SCF-SPACE:
c***** write (6,*) 'Sixtuple OCCU[1,1,0,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,1,0,0]
c      ***********************************

       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO

       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC-1 
       DO 3500 I1=1,MAX1
c      ------------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 0



       MAX2 = N1-1
       DO 3505 I2=1,MAX2
c      -------------------
c      *
       N2 = N1 - I2
       NEL(N2) = 0



       MAX3 = NOCC-3
       DO 3510 I3=1,MAX3
c      ------------------
c      *
       NN3 = NOCC + 1 - I3
       IF (I3.EQ.1) GO TO 3511
       NN3 = LL3 - 1
 3511  A=0
       call SHIFT2(N1,N2,NN3,LL3)
       NEL(LL3) = 1




       MAX4 = NN3 - 1  
 3512  A=0
       DO 3515 I4=1,MAX4
c      ------------------
c      *
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 3516
       NN4 = LL4 - 1
 3516  A=0
       call SHIFT2(N1,N2,NN4,LL4)
       IF (LL4.EQ.0) GO TO 3517
       NEL(LL4) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      --------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3532
       GO TO 5000
 3532  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3535
       M1 = L1 - 1
 3535  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3536
       LIM4 = LIM4 - 1
 3536  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3537
       M2 = L2 - 1
 3537  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3538 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3540
       GO TO 5000
 3540  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L2) = 0
       ENDDO
 3538  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3545
       M1 = L1 - 1
 3545  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3546
       M2 = L2 - 1
 3546  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3547
       M3 = L3 - 1
 3547  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3548
       M4 = L4 - 1
 3548  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3550
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3551
       GO TO 5000
 3551  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3550  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3555
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3556
       GO TO 5000
 3556  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3555  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
C*     -----------------------------------

       NEL(LL4)=2
 3515  CONTINUE
 3517  A=0

 
       NEL(LL3)=2
 3510  CONTINUE


       NEL(N2)=2
 3505  CONTINUE


       NEL(N1)=2
 3500  CONTINUE
c      ***

c--------------------------------------------------
c***** write (6,*) 'Sixtuple OCCU[1,1,0,0] done'
c***** write (6,*) '---------------------------'





c      ********
 3600  A=0
c      ********


c      SCF-SPACE:
c****  write (6,*) 'Sixtuple OCCU[1,1,1,1,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,1,1,1,0]
c      ***********************************

       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO

       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC
       DO 3605 I1=1,NOCC
c      -----------------
c      *
       N1 = NOCC + 1 - I1
       NEL(N1) = 0


       MAX2 = NOCC-4
       DO 3610 I2=1,MAX2
c      ------------------
c      *
       NN2 = NOCC + 1 - I2
       IF (I2.EQ.1) GO TO 3611
       NN2 = LL2 - 1
 3611  A=0
       call SHIFT1(N1,NN2,LL2)
       NEL(LL2) = 1



       MAX3 = NN2 - 3
       DO 3615 I3=1,MAX3
c      -----------------
c      *
       NN3 = LL2 - I3
       IF (I3.EQ.1) GO TO 3616
       NN3 = LL3 - 1
 3616  A=0
       call SHIFT1(N1,NN3,LL3)
       NEL(LL3) = 1




       MAX4 = NN3 - 2
       DO 3620 I4=1,MAX4
c      -----------------
c      *
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 3621
       NN4 = LL4 - 1
 3621  A=0
       call SHIFT1(N1,NN4,LL4)
       NEL(LL4) = 1


       MAX5 = NN4 - 1
       DO 3625 I5=1,MAX5
c      ------------------
c      *
       NN5 = LL4 - I5
       IF (I5.EQ.1) GO TO 3626
       NN5 = LL5 - 1
 3626  A=0
       call SHIFT1(N1,NN5,LL5)
       IF (LL5.EQ.0) GO TO 3627
       NEL(LL5) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      --------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3632
       GO TO 5000
 3632  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3635
       M1 = L1 - 1
 3635  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3636
       LIM4 = LIM4 - 1
 3636  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3637
       M2 = L2 - 1
 3637  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3638 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3640
       GO TO 5000
 3640  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L2) = 0
       ENDDO
 3638  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3645
       M1 = L1 - 1
 3645  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3646
       M2 = L2 - 1
 3646  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3647
       M3 = L3 - 1
 3647  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3648
       M4 = L4 - 1
 3648  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3650
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3651
       GO TO 5000
 3651  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3650  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3655
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3656
       GO TO 5000
 3656  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3655  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
C*     -----------------------------------

       NEL(LL5)=2
 3625  CONTINUE
 3627  A=0


       NEL(LL4)=2 
 3620  CONTINUE


       NEL(LL3)=2
 3615  CONTINUE


       NEL(LL2)=2
 3610  CONTINUE


       NEL(N1)=2
 3605  CONTINUE
c      ******

c---------------------------------------------------------
c***** write (6,*) 'Sixtuple OCCU[1,1,1,1,0] done'
c***** write (6,*) '---------------------------'
c      ----------------------------------------------






c      ********
 3800  A=0
c      ********




c      SCF-SPACE:
c***** write (6,*) 'Sixtuple OCCU[1,1,1,1,1,1] starts'
c
c      **************************************
C      occupied MO space: CASE  [1,1,1,1,1,1]
c      **************************************

       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       MAX1 = NOCC-5
       DO 3805 I1=1,MAX1
c      -----------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 1

       MAX2 = N1 - 5 
       DO 3810 I2=1,MAX2
c      -----------------
c      *
       N2 = N1 - I2
       NEL(N2) = 1

       MAX3 = N2 - 4
       DO 3815 I3=1,MAX3
c      ------------------
c      *
       N3 = N2 - I3
       NEL(N3) = 1

       MAX4 = N3 - 3
       DO 3820 I4=1,MAX4
c      -----------------
c      *
       N4 = N3 - I4
       NEL(N4) = 1

       MAX5 = N4 - 2 
       DO 3825 I5=1,MAX5
c      -----------------
c      *
       N5 = N4 - I5
       NEL(N5) = 1

       MAX6 = N5 - 1
       DO 3830 I6=1,MAX6
c      ------------------
c      *
       N6 = N5 - I6
       IF (N6.GT.0) GO TO 3831 
       write (6,*) 'N6=',N6
       STOP
 3831  A=0
       NEL(N6) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      *************************
c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3832
       GO TO 5000
 3832  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      ***********************************


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3835
       M1 = L1 - 1
 3835  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3836
       LIM4 = LIM4 - 1
 3836  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3837
       M2 = L2 - 1
 3837  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3838 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3840
       GO TO 5000
 3840  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L2) = 0
       ENDDO
 3838  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3845
       M1 = L1 - 1
 3845  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3846
       M2 = L2 - 1
 3846  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3847
       M3 = L3 - 1
 3847  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3848
       M4 = L4 - 1
 3848  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3850
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      ************************
       IF (ICOUNT.NE.NSPRO) GO TO 3851
       GO TO 5000
 3851  A=0
c      ************************
c      ----------------------------------
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3850  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3855
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       IF (ICOUNT.NE.NSPRO) GO TO 3856
       GO TO 5000
 3856  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3855  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *************************************
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
c      *************************************

       NEL(N6)=2
 3830  CONTINUE


       NEL(N5)=2
 3825  CONTINUE


       NEL(N4)=2
 3820  CONTINUE


       NEL(N3)=2
 3815  CONTINUE


       NEL(N2)=2
 3810  CONTINUE


       NEL(N1)=2
 3805  CONTINUE

     

c***** write (6,*) 'Sixtuple OCCU[1,1,1,1,1,1] done'
c***** write (6,*) '------------------------------------'


c      ----------------------------------------------
c***** write (6,*) 'subr. LABELSDTQ: 6- got here-all'
       GO TO 5000
c      ----------







c*******************END OF SDTQ56-excitations************
c                   -------------------------
c******LABELSDTQ:  SDTQ56-excitation-labeling ends ******




c      ******
 5000  AAA=0
c      ******




c      -----------------
       DO JAM=1,NTOT
c      -----------------
c      *
       ICON(JAM) = NEL(JAM)
c      *
       ENDDO




c      ---------------------------------
       IF (ICOUNT.LE.mxsprod) GO TO 5001
c      ---------------------------------
c***   write (6,*) 'subr. RLABELSDTQ: exceeded limits'
c***   write (6,*) 'IEXCIT=',IEXCIT
c***   write (6,*) 'ICOUNT=',ICOUNT
c****  STOP
c      *
 5001  A=0
c      ***



c***************************************************
c      *
c*     write (6,*) 'IEXCIT=',IEXCIT
c*     write (6,*) 'NSPRO=',NSPRO
c*     write (6,*) (ICON(LIL), LIL=1,NTOT) 
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
c*     write (6,*) 'subr. RLABELSDTQ will stop'
c      *****************
c      *
c      -----------------------
c      subr RLABELSDTQ ends
c      -----------------------
       RETURN
       END








c     ---------------------------------------------------------
c     ***
      SUBROUTINE LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT,ICON,NSPRO)
c
c     ---------------------------------------------------------
c     This subroutine
c     ---------------
c     * given  ICON(j=1,16)
c     * returns NSPRO (=Number of space product SD-only)
c         LABELSDTQ
c      This subroutine labels SD-space products in a 
c      systematic way
c      ---------------------------------------------------

      implicit double precision(a-h,o-z)

      parameter (mxorb = 100)
      parameter (mxsprod = 2400000)
      parameter (mxsprod1 = 100000)
      parameter (mxCI = 1000000)
      parameter (mxstring = 2000000)
      parameter (mxna = 200)
      DIMENSION ICON(mxorb),NEL(mxorb)

       ICOUNT = 0
       IF (IEXCIT.EQ.1) GO TO 33
       IF (IEXCIT.EQ.2) GO TO 44 
       IF (IEXCIT.EQ.3) GO TO 350 
       IF (IEXCIT.EQ.4) GO TO 490 
       IF (IEXCIT.EQ.5) GO TO 1100 
       IF (IEXCIT.EQ.6) GO TO 3000 
       write (6,*) 'subr.LABELSDTQ:  no excitations '
       STOP


C***** code to generate space-products


C      one-determinant (SCF) :
C      -----------------------
       DO 30 I=1,NOCC
       NEL(I) = 2
 30    CONTINUE
       DO 32 I=NOCC+1,NTOT
       NEL(I) = 0
 32    CONTINUE

  33   A=0
c      'single excitations '
C      -------------------
       DO 35 I=1,NOCC
       DO 36 J=1,NOCC
       NEL(J)= 2
  36   CONTINUE

       NSING = NOCC+1-I
       NEL(NSING) = 1

c      "virtual space"
c      
       DO 39 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 39    CONTINUE

       DO 40 MM=NOCC+1,NTOT
c      
       NEL(MM) = 1 
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 43
       ENDDO
       GO TO 5000
  43   A=0
       NEL(MM) = 0
 40    CONTINUE
 35    CONTINUE

c      ----------
       GO TO 5000
c      ----------




  44   A=0
c      'DIAGONAL-double excitations'

       DO 45 I=1,NOCC
c      -------------------
       DO 47 J=1,NOCC
       NEL(J)= 2
  47   CONTINUE

       NSING = NOCC+1-I
       NEL(NSING) = 0 

c      "virtual"
       DO 49 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 49    CONTINUE

       DO 50 MM=NOCC+1,NTOT
       NEL(MM) = 2 
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
c      
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 51 
       ENDDO
       GO TO 5000
  51   A=0
       NEL(MM) = 0
 50    CONTINUE

c      "virtual"
       DO 55 JJJ=1,NVIRT
c      
       DO 56 KKK=NOCC+1,NTOT
       NEL(KKK) = 0
  56   CONTINUE
       M = NTOT+1-JJJ
       NEL(M) = 1

       DO 60 II=1,M-NOCC-1
       MMM = M - II
       NEL(MMM) = 1

       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 61 
       ENDDO
       GO TO 5000
  61   A=0
       NEL(M-II) = 0
 60    CONTINUE
       NEL(M) = 0
 55    CONTINUE
 45    CONTINUE


c      'mixed double excitations '
       DO 275 I=1,NOCC
       DO 277 J=1,NOCC
       NEL(J)= 2
 277   CONTINUE

       M1 = NOCC+1-I 
       NEL(M1) = 1 

       DO 278 K1=1,M1-1
       M2 = M1 - K1
       NEL(M2)= 1

c      "virtual"
       DO 280 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 280   CONTINUE
       DO 285 MM=NOCC+1,NTOT
       NEL(MM) = 2 
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 281 
       ENDDO
       GO TO 5000
 281   A=0
       NEL(MM) = 0
 285   CONTINUE

c      "virtual"
       DO 300 JJJ=1,NVIRT
       DO 301 KKK=NOCC+1,NTOT
       NEL(KKK) = 0
 301   CONTINUE
       M = NTOT+1-JJJ
       NEL(M) = 1
       DO 305 II=1,M-NOCC-1
       MM = M - II
       NEL(MM) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 303 
       ENDDO
       GO TO 5000
 303   A=0
       NEL(MM) = 0
 305   CONTINUE
       NEL(M) = 0
 300   CONTINUE
c      "SCF-MOs"
       NEL(M1-K1) = 2
 278   CONTINUE
       NEL(M1) = 2
 275   CONTINUE
       GO TO 5000



 350   A=0
c      'TRIPLE excitations '
C      'T-case:[1,0] 
       DO 362 I= 1,NOCC
       NEL(I) = 2
 362   CONTINUE
       DO 363 KKZ=1,NOCC
       NXZ = NOCC + 1 - KKZ
       NEL(NXZ) = 0
       DO 367 KK1Z=1,NOCC-1
       M1Z = NOCC + 1 - KK1Z 
       IF (M1Z.GT.NXZ) GO TO 371 
       IF (M1Z.LE.NXZ) GO TO 372 
 371   AAA=0
       L1Z = M1Z
       GO TO 373 
 372   AAA=0
       L1Z = M1Z - 1
 373   AAA=0
       NEL(L1Z) = 1

c      "Virtual MOs"
c      TRIPLES: [2,1]
       DO 380 K = NOCC+1,NTOT
       NEL(K) = 0
 380   CONTINUE
       DO 385 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2
       DO 389 KK1=1,NVIRT-1
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 390 
       IF (M1.LE.NX) GO TO 391
 390   AAA=0
       L1 = M1
       GO TO 392
 391   AAA=0
       L1 = M1 - 1
 392   AAA=0
       NEL(L1) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 400 
       ENDDO
       GO TO 5000
 400   A=0
 395   AAA=0.0
       NEL(L1) = 0
 389   CONTINUE
 388   AAA=0.0
       NEL(NX) = 0
 385   CONTINUE

c      "virtual MOs"
c      TRIPLES: [1,1,1]
       DO 397 K = NOCC+1,NTOT
       NEL(K) = 0
 397   CONTINUE
       DO 398 KK = 1,NVIRT-2
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 402 KK1 = 1,M1-NOCC-2
       M2 = M1 - KK1
       NEL(M2) = 1
       DO 406 KK2 = 1,M2-NOCC-1
       M3 = M2 - KK2
       NEL(M3) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 410 
       ENDDO
       GO TO 5000
 410   A=0
 409   AAA=0.0
       NEL(M3) = 0
 406   CONTINUE
 405   AAA=0.0
       NEL(M2) = 0
 402   CONTINUE
 401   AAA=0.0
       NEL(M1) = 0
 398   CONTINUE
 381   A=0.0
 412   AAA=0.0 
       NEL(L1Z) = 2
 367   CONTINUE
 411   AAA=0.0
       NEL(NXZ) = 2
 363   CONTINUE
 416   AAA=0.0


c      "occupied MOs TRIPLES: [1,1,1]"
       DO 414 I= 1,NOCC
       NEL(I) = 2
 414   CONTINUE
       DO 415 KKW = 1,NOCC-2
       M1W = NOCC+1-KKW
       NEL(M1W) = 1
       DO 418 KK1W = 1,M1W-2
       M2W = M1W - KK1W
       NEL(M2W) = 1

       DO 421 KK2W = 1,M2W-1
       M3W = M2W - KK2W
       NEL(M3W) = 1

c      "Virtual MOs TRIPLES: [2,1]"
c      ---------------------------
       DO 424 K = NOCC+1,NTOT
       NEL(K) = 0
 424   CONTINUE
       DO 425 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2

       DO 430 KK1=1,NVIRT-1
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 431 
       IF (M1.LE.NX) GO TO 432
 431   AAA=0
       L1 = M1
       GO TO 433
 432   AAA=0
       L1 = M1 - 1
 433   AAA=0
       NEL(L1) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 440 
       ENDDO
       GO TO 5000
 440   A=0
 439   AAA=0.0
       NEL(L1) = 0
 430   CONTINUE
 438   AAA=0.0
       NEL(NX) = 0
 425   CONTINUE

c      virtual MOs TRIPLES: [1,1,1]"
c      ------------------------------
       DO 441 K = NOCC+1,NTOT
       NEL(K) = 0
 441   CONTINUE
       DO 442 KK = 1,NVIRT-2
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 446 KK1 = 1,M1-NOCC-2
       M2 = M1 - KK1
       NEL(M2) = 1
       DO 450 KK2 = 1,M2-NOCC-1
       M3 = M2 - KK2
       NEL(M3) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 451 
       ENDDO
       GO TO 5000
 451   A=0
 453   A=0
       NEL(M3) = 0
 450   CONTINUE
 449   AAA=0.0
       NEL(M2) = 0
 446   CONTINUE
 445   AAA=0.0
       NEL(M1) = 0
 442   CONTINUE

c      occupied MOs TRIPLES: [1,1,1]
c      -------------------------------
 457   AAA=0.0
       NEL(M3W) = 2 
 421   CONTINUE
 456   AAA=0.0
       NEL(M2W) = 2 
 418   CONTINUE
 455   AAA=0.0
       NEL(M1W) = 2 
 415   CONTINUE
       GO TO 5000
c      ----------



  490  AAA=0.0
c      -------------------------- 
c      'QUADRUPLE EXCITATIONS"
c      -------------------------- 

c*     "OCCUPIED QUADRUPLES: Q-case: [0,0]'
       DO 503 I= 1,NOCC
       NEL(I) = 2
 503   CONTINUE
       DO 505 I=1,NOCC-1
       M = NOCC+1-I
       NEL(M) = 0
       DO 510 J=1,M-1
       MM = M - J
       NEL(MM) = 0

c      "virtual space: QUADRUPLES: [2,2]
       DO 520 K = NOCC+1,NTOT
       NEL(K) = 0
 520   CONTINUE
       DO 525 KK=1,NVIRT-1
       MK = NTOT+1-KK
       NEL(MK) = 2
       DO 530 KKK=1,NVIRT-KK
       MMK = MK - KKK
       NEL(MMK) = 2
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 535 
       ENDDO
       GO TO 5000
 535   A=0
 533   AAA=0.0
       NEL(MMK) = 0
 530   CONTINUE
 528   AAA=0.0
       NEL(MK) = 0
 525   CONTINUE

c      "Virtual MOs QUADRUPLES: [2,1,1]"
       DO 540 K = NOCC+1,NTOT
       NEL(K) = 0
 540   CONTINUE
       DO 545 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2
       DO 550 KK1=1,NVIRT-2
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 551 
       IF (M1.LE.NX) GO TO 552
 551   AAA=0
       L1 = M1
       GO TO 553
 552   AAA=0
       L1 = M1 - 1
 553   AAA=0
       NEL(L1) = 1
       DO 560 KK2=1,M1-NOCC-2    
       M2 = M1 - KK2
       IF (M2.GT.NX) GO TO 561 
       IF (M2.LE.NX) GO TO 562 
 561   AAA=0
       L2 = M2
       GO TO 563
 562   AAA=0
       L2 = M2 - 1
 563   AAA=0
       NEL(L2) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 565
       ENDDO
       GO TO 5000
 565   A=0
 566   AAA=0.0
       NEL(L2) = 0
 560   CONTINUE
 556   AAA=0.0
       NEL(L1) = 0
 550   CONTINUE
 548   AAA=0.0
       NEL(NX) = 0
 545   CONTINUE

c      "virtual MOs QUADRUPLES: [1,1,1,1]"
c      -----------------------------------
       DO 569 K = NOCC+1,NTOT
       NEL(K) = 0
 569   CONTINUE
       DO 570 KK = 1,NVIRT-3
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 575 KK1 = 1,M1-NOCC-3
       M2 = M1 - KK1
       NEL(M2) = 1
       DO 580 KK2 = 1,M2-NOCC-2
       M3 = M2 - KK2
       NEL(M3) = 1
       DO 585 KK3 = 1,M3-NOCC-1
       M4 = M3 - KK3
       NEL(M4) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 590
       ENDDO
       GO TO 5000
 590   A=0
 588   AAA=0.0
       NEL(M4) = 0
 585   CONTINUE
 583   AAA=0.0
       NEL(M3) = 0
 580   CONTINUE
 578   AAA=0.0
       NEL(M2) = 0
 575   CONTINUE
 573   AAA=0.0
       NEL(M1) = 0
 570   CONTINUE
 592   AAA=0.0
       NEL(MM) = 2
 510   CONTINUE
 591   AAA=0.0 
       NEL(M) = 2
 505   CONTINUE


 595   AAA=0.0
c      "Occupied QUADRUPLES: [0,1,1]"
c      ------------------------------
       DO 605 I= 1,NOCC
       NEL(I) = 2
 605   CONTINUE
       DO 620 KKZ=1,NOCC
       NXZ = NOCC + 1 - KKZ
       NEL(NXZ) = 0
       DO 627 KK1Z=1,NOCC-2
       M1Z = NOCC + 1 - KK1Z 
       IF (M1Z.GT.NXZ) GO TO 631 
       IF (M1Z.LE.NXZ) GO TO 632 
 631   AAA=0
       L1Z = M1Z
       GO TO 633 
 632   AAA=0
       L1Z = M1Z - 1
 633   AAA=0
       NEL(L1Z) = 1

       DO 640 KK2Z=1,M1Z-2    
       M2Z = M1Z - KK2Z
       IF (M2Z.GT.NXZ) GO TO 641 
       IF (M2Z.LE.NXZ) GO TO 642 
 641   AAA=0
       L2Z = M2Z
       GO TO 643
 642   AAA=0
       L2Z = M2Z - 1
 643   AAA=0
       NEL(L2Z) = 1

c      "virtual space QUADRUPLES: [2,2]"
c      -----------------------------------
       DO 650 K = NOCC+1,NTOT
       NEL(K) = 0
 650   CONTINUE
       DO 655 KK=1,NVIRT-1
       MK = NTOT+1-KK
       NEL(MK) = 2
       DO 660 KKK=1,NVIRT-KK
       MMK = MK - KKK
       NEL(MMK) = 2
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 665 
       ENDDO
       GO TO 5000
 665   A=0
 663   AAA=0.0
       NEL(MMK) = 0
 660   CONTINUE
 658   AAA=0.0
       NEL(MK) = 0
 655   CONTINUE


c      "Virtual MOs QUADRUPLES: [2,1,1]"
c      ----------------------------------
       DO 670 K = NOCC+1,NTOT
       NEL(K) = 0
 670   CONTINUE
       DO 675 KK=1,NVIRT
       NX = NTOT + 1 - KK
       NEL(NX) = 2
       DO 680 KK1=1,NVIRT-2
       M1 = NTOT + 1 - KK1 
       IF (M1.GT.NX) GO TO 681 
       IF (M1.LE.NX) GO TO 682
 681   AAA=0
       L1 = M1
       GO TO 683
 682   AAA=0
       L1 = M1 - 1
 683   AAA=0
       NEL(L1) = 1
       DO 690 KK2=1,M1-NOCC-2    
       M2 = M1 - KK2
       IF (M2.GT.NX) GO TO 691 
       IF (M2.LE.NX) GO TO 692 
 691   AAA=0
       L2 = M2
       GO TO 693
 692   AAA=0
       L2 = M2 - 1
 693   AAA=0
       NEL(L2) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 698 
       ENDDO
       GO TO 5000
 698   A=0
 696   AAA=0.0
       NEL(L2) = 0
 690   CONTINUE
 686   AAA=0.0
       NEL(L1) = 0
 680   CONTINUE
 678   AAA=0.0
       NEL(NX) = 0
 675   CONTINUE

c      "virtual MOs QUADRUPLES: [1,1,1,1]"
c      --------------------------------------
       DO 699 K = NOCC+1,NTOT
       NEL(K) = 0
 699   CONTINUE
       DO 700 KK = 1,NVIRT-3
       M1 = NTOT+1-KK
       NEL(M1) = 1
       DO 705 KK1 = 1,M1-NOCC-3
       M2 = M1 - KK1
       NEL(M2) = 1
       DO 710 KK2 = 1,M2-NOCC-2
       M3 = M2 - KK2
       NEL(M3) = 1
       DO 715 KK3 = 1,M3-NOCC-1
       M4 = M3 - KK3
       NEL(M4) = 1
       ICOUNT = ICOUNT + 1
       DO JAM=1,NTOT
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 719 
       ENDDO
       GO TO 5000
 719   A=0
 718   AAA=0.0
       NEL(M4) = 0
 715   CONTINUE
 713   AAA=0.0
       NEL(M3) = 0
 710   CONTINUE
 708   AAA=0.0
       NEL(M2) = 0
 705   CONTINUE
 703   AAA=0.0
       NEL(M1) = 0
 700   CONTINUE
 722   AAA=0.0
       NEL(L2Z) = 2
 640   CONTINUE
 721   AAA=0.0
       NEL(L1Z) = 2
 627   CONTINUE
 720   AAA=0.0
       NEL(NXZ) = 2
 620   CONTINUE



 750   AAA=0
c*     write (6,*) 'Q-case: [1,1,1,1]'
c      ********************************
c      occupied MOs CASE-C) [1,1,1,1]
c      ********************************
       DO 803 I= 1,NOCC
       NEL(I) = 2
 803   CONTINUE


       DO 805 KKW = 1,NOCC-3
c      ---------------------
       M1W = NOCC+1-KKW
c      -------------
       NEL(M1W) = 1


       DO 810 KK1W = 1,M1W-3
c      -------------------------
       M2W = M1W - KK1W
c      ------------
       NEL(M2W) = 1


       DO 815 KK2W = 1,M2W-2
c      -------------------------
       M3W = M2W - KK2W
c      ------------
       NEL(M3W) = 1


       DO 820 KK3W = 1,M3W-1
c      -------------------------
       M4W = M3W - KK3W
c      ------------
       NEL(M4W) = 1



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *************
c      CASE-1) [2,2]
c      *************
       DO 850 K = NOCC+1,NTOT
       NEL(K) = 0
 850   CONTINUE


       DO 855 KK=1,NVIRT-1
c      -------------------
       MK = NTOT+1-KK
c      -----------
       NEL(MK) = 2


       DO 860 KKK=1,NVIRT-KK
c      ----------------------
       MMK = MK - KKK
c      -------------
       NEL(MMK) = 2


c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 865 
       ENDDO
       GO TO 5000
 865   A=0
c      *************

 863   AAA=0.0
       NEL(MMK) = 0
 860   CONTINUE


 858   AAA=0.0
       NEL(MK) = 0
 855   CONTINUE



c      ^^^^^^^^^^^^^^^^
c      Virtual MOs
c      ^^^^^^^^^^^^^^^^
c      ****************
c      CASE-2 [2,1,1]
c      ****************
       DO 870 K = NOCC+1,NTOT
       NEL(K) = 0
 870   CONTINUE


       DO 875 KK=1,NVIRT
c      -------------------
       NX = NTOT + 1 - KK
c      -----------
       NEL(NX) = 2


       DO 880 KK1=1,NVIRT-2
c      ----------------------
       M1 = NTOT + 1 - KK1 
c      ----------------------
       IF (M1.GT.NX) GO TO 881 
       IF (M1.LE.NX) GO TO 882
 881   AAA=0
       L1 = M1
       GO TO 883
 882   AAA=0
       L1 = M1 - 1
 883   AAA=0
c      -----------
       NEL(L1) = 1


       DO 890 KK2=1,M1-NOCC-2    
c      ------------------------
       M2 = M1 - KK2
c      ------------------------
       IF (M2.GT.NX) GO TO 891 
       IF (M2.LE.NX) GO TO 892 
 891   AAA=0
       L2 = M2
       GO TO 893
 892   AAA=0
       L2 = M2 - 1
 893   AAA=0
c      -----------
       NEL(L2) = 1

c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 898 
       ENDDO
       GO TO 5000
 898   A=0
c      *************

 896   AAA=0.0
       NEL(L2) = 0
c      ---------
 890   CONTINUE


 886   AAA=0.0
       NEL(L1) = 0
c      ---------
 880   CONTINUE


 878   AAA=0.0
       NEL(NX) = 0
c      ---------
 875   CONTINUE




c      ^^^^^^^^^^^^^^^^^
c      virtual MOs
c      ^^^^^^^^^^^^^^^^^
c      *****************
c      CASE-3) [1,1,1,1]
c      *****************
       DO 899 K = NOCC+1,NTOT
       NEL(K) = 0
 899   CONTINUE


       DO 900 KK = 1,NVIRT-3
c      ---------------------
       M1 = NTOT+1-KK
c      -----------
       NEL(M1) = 1


       DO 905 KK1 = 1,M1-NOCC-3
c      -------------------------
       M2 = M1 - KK1
c      -----------
       NEL(M2) = 1


       DO 910 KK2 = 1,M2-NOCC-2
c      -------------------------
       M3 = M2 - KK2
c      -----------
       NEL(M3) = 1


       DO 915 KK3 = 1,M3-NOCC-1
c      -------------------------
       M4 = M3 - KK3
c      -----------
       NEL(M4) = 1


c      ----------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c*     write (6,*) (NEL(LIL), LIL=1,NTOT) 
C*     -----------------------------------
c      *************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 919 
       ENDDO
       GO TO 5000
 919   A=0
c      *************

 918   AAA=0.0
       NEL(M4) = 0
c      ----------
 915   CONTINUE


 913   AAA=0.0
       NEL(M3) = 0
c      ----------
 910   CONTINUE


 908   AAA=0.0
       NEL(M2) = 0
c      ----------
 905   CONTINUE


 903   AAA=0.0
       NEL(M1) = 0
c      ----------
 900   CONTINUE


c      -------------------------------
c      occupied MOs CASE-C) [1,1,1,1]
c      -------------------------------
 923   AAA=0.0
       NEL(M4W) = 2 
c      ----------
 820   CONTINUE


 922   AAA=0.0
       NEL(M3W) = 2 
c      ----------
 815   CONTINUE


 921   AAA=0.0
       NEL(M2W) = 2 
c      ----------
 810   CONTINUE


 920   AAA=0.0
       NEL(M1W) = 2 
c      ----------
 805   CONTINUE
c      ----------------------------------------------
c      ----------------------------------------------
c***** write (6,*) 'subr. LABELSDTQ: Q- got here-all'
       GO TO 5000
c      ----------





c*****LABELSDTQ: SDTQ56-excitation-labeling starts ***
c     --------------------------------------------


c      ---------------------------------------
c      subr. LABELSDTQ:  QUINTUPLE EXCITATIONS
c      ---------------------------------------
c***** write (6,*) 'quintuple(5) exc. start here'
c***** write (6,*) '5-case: [1,0,0] '
c

c      ***
 1100  A=0
c      ***


c      ****************************
c      occ. MO space: CASE [1,0,0]
c      ****************************
c*
c      SCF-SPACE:
c***** write (6,*) 'Quintuple OCCU[1,0,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,0,0]
c      ***********************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC-1 
       DO 1500 I1=1,MAX1
c      ------------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 0



       MAX2 = N1-1
       DO 1505 I2=1,MAX2
c      -------------------
c      *
       N2 = N1 - I2
       NEL(N2) = 0



       MAX3 = NOCC-2
       DO 1510 I3=1,MAX3
c      ------------------
c      *
       NN3 = NOCC + 1 - I3
       IF (I3.EQ.1) GO TO 1511
       NN3 = LL3 - 1
 1511  A=0
       call SHIFT2(N1,N2,NN3,LL3)
       IF (LL3.EQ.0) GO TO 1512 
       NEL(LL3) = 1



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c
c      write (6,*) 'VIRT [2,2,1] starts'
c      *********************
c      CASE) VIRT [2,2,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 


       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 


       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1535
       M1 = L1 - 1
 1535  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1538 
       NEL(L1) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1540 
       ENDDO
       GO TO 5000
 1540  A=0
c      *****************
c      ----------------------------------
       NEL(L1) = 0
       ENDDO
 1538  A=0
c      ***

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'




c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2


       LIM1 = NVIRT-3
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1545
       M1 = L1 - 1
 1545  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1546
       M2 = L2 - 1
 1546  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1547
       M3 = L3 - 1
 1547  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1550
       NEL(L3) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1549 
       ENDDO
       GO TO 5000
 1549  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L3)=0
       ENDDO
 1550  A=0
c      ****

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1555
       NEL(M5) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1556 
       ENDDO
       GO TO 5000
 1556  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M5)=0
       ENDDO
c      *****
 1555  A=0

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
C*     ---------------------------------------------


c      -----------------------------
c      occupied MOs: CASE [1,0,0]
c      -----------------------------
       NEL(LL3) = 2 
 1510  CONTINUE
 1512  A=0
c      ***

       NEL(N2) = 2 
 1505  CONTINUE

       NEL(N1) = 2 
 1500  CONTINUE
c      ********





c***** write (6,*) 'Quintuple OCCU[1,0,0] done'
c***** write (6,*) '---------------------------'




c      *****
 1600  A=0
c      *****



cLAIMIS
c      SCF-SPACE:
c***** write (6,*) 'Quintuple OCCU[1,1,1,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,1,1,0]
c      ***********************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC
       DO 1605 I1=1,NOCC
c      -----------------
c      *
       N1 = NOCC + 1 - I1
       NEL(N1) = 0


       MAX2 = NOCC-3
       DO 1610 I2=1,MAX2
c      ------------------
c      *
       NN2 = NOCC + 1 - I2
       IF (I2.EQ.1) GO TO 1611
       NN2 = LL2 - 1
 1611  A=0
       call SHIFT1(N1,NN2,LL2)
       NEL(LL2) = 1



       MAX3 = NN2 - 2 
       DO 1615 I3=1,MAX3
c      -----------------
c      *
       NN3 = LL2 - I3
       IF (I3.EQ.1) GO TO 1616
       NN3 = LL3 - 1
 1616  A=0
       call SHIFT1(N1,NN3,LL3)
       NEL(LL3) = 1



       MAX4 = NN3 - 1 
       DO 1620 I4=1,MAX4
c      -----------------
c      *
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 1621
       NN4 = LL4 - 1
 1621  A=0
       call SHIFT1(N1,NN4,LL4)
       IF (LL4.EQ.0) GO TO 1622 
       NEL(LL4) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1] starts'
c      *********************
c      CASE) VIRT [2,2,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 


       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 


       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1635
       M1 = L1 - 1
 1635  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1638 
       NEL(L1) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1640 
       ENDDO
       GO TO 5000
 1640  A=0
c      *****************
c      ----------------------------------
c      ----------------------------------
       NEL(L1) = 0
       ENDDO
 1638  A=0
c      ***

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'




c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2


       LIM1 = NVIRT-3
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1645
       M1 = L1 - 1
 1645  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1646
       M2 = L2 - 1
 1646  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1647
       M3 = L3 - 1
 1647  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1650
       NEL(L3) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1651 
       ENDDO
       GO TO 5000
 1651  A=0
c      *****************
c      ----------------------------------
       NEL(L3)=0
       ENDDO
 1650  A=0
c      ****

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1655
       NEL(M5) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1656
       ENDDO
       GO TO 5000
 1656  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M5)=0
       ENDDO
c      *****
 1655  A=0

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
C*     ---------------------------------------------

       NEL(LL4)=2
 1620  CONTINUE
 1622  A=0
c      ***

       NEL(LL3)=2
 1615  CONTINUE


       NEL(LL2)=2
 1610  CONTINUE


       NEL(N1)=2
 1605  CONTINUE



c***** write (6,*) 'Quintuple OCCU=[1,1,1,0] done'
c***** write (6,*) '---------------------------'







c      *****
 1700  A=0
c      *****




cLAIMIS
c      Quintuple OCCU=[1,1,1,1,1] starts
c      *********************************
c      SCF-SPACE:
c***** write (6,*) 'Quintuple OCCU=[1,1,1,1,1] starts'
c      **************************************
C      occupied MO space: CASE  [1,1,1,1,1]
c      **************************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       MAX1 = NOCC-4
       DO 1805 I1=1,MAX1
c      -----------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 1



       MAX2 = N1 - 4 
       DO 1810 I2=1,MAX2
c      -----------------
c      *
       N2 = N1 - I2
       NEL(N2) = 1


       MAX3 = N2 - 3 
       DO 1815 I3=1,MAX3
c      ------------------
c      *
       N3 = N2 - I3
       NEL(N3) = 1


       MAX4 = N3 - 2 
       DO 1820 I4=1,MAX4
c      -----------------
c      *
       N4 = N3 - I4
       NEL(N4) = 1


       MAX5 = N4 - 1 
       DO 1825 I5=1,MAX5
c      -----------------
c      *
       N5 = N4 - I5
       IF (N5.GT.0) GO TO 1831 
       write (6,*) 'N5=',N5
       STOP
 1831  A=0
       NEL(N5) = 1



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1] starts'
c      *********************
c      CASE) VIRT [2,2,1]
c      *********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 


       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 


       LIM3 = NVIRT - 2 
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1835
       M1 = L1 - 1
 1835  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       IF (L1.EQ.NOCC) GO TO 1838 
       NEL(L1) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1840
       ENDDO
       GO TO 5000
 1840  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L1) = 0
       ENDDO
 1838  A=0
c      ***

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Quintuple=[2,2,1] virt done'




c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2


       LIM1 = NVIRT-3
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 1845
       M1 = L1 - 1
 1845  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-2
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 1846
       M2 = L2 - 1
 1846  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-1 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 1847
       M3 = L3 - 1
 1847  A=0
       call SHIFT1(NX,M3,L3)
       IF (L3.EQ.NOCC) GO TO 1850
       NEL(L3) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1851
       ENDDO
       GO TO 5000
 1851  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L3)=0
       ENDDO
 1850  A=0
c      ****

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Pentuple=[2,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-4
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-3
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-2
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-1
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       IF (M5.EQ.NOCC) GO TO 1855
       NEL(M5) = 1


C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 1856
       ENDDO
       GO TO 5000
 1856  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M5)=0
       ENDDO
c      *****
 1855  A=0

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Quintuple virt[1,1,1,1,1] done'
C*     ---------------------------------------------

       NEL(N5)=2
 1825  CONTINUE


       NEL(N4)=2
 1820  CONTINUE


       NEL(N3)=2
 1815  CONTINUE


       NEL(N2)=2
 1810  CONTINUE


       NEL(N1)=2
 1805  CONTINUE


c---------------------------------------------------
c***** write (6,*) 'Quintuple OCCU=[1,1,1,1,1] done'
c                   LABELSDTQ:
c***** write (6,*) '---------------------------'
c      ----------------------------------------------
c***** write (6,*) 'subr. LABELSDTQ: 5- got here-all'
       GO TO 5000
c      ----------












c      **********
 3000  A=0
c      **********


c     ----------------------
c***** write (6,*) 'sixtuple exc'
c***** write (6,*) '6-case: [0,0,0]'
c     ----------------------
c*
c*     call CHECKMATE(NOCC,NSCF,MATE,M,MM,MMM,
c*   * MK,MMK,MMK2,ISAY)
c*     IF (ISAY.EQ.0) GO TO 3200 
c      write (6,*) M,MM,MMM,MK,MMK,MMK2
c      write (6,*) ISAY
c
c
c     ----------------------
c     SIXTUPLE-EXCITATIONS:
c     ----------------------
c


c      *******
 3100  A=0
c      *******


c      *********************************
C      occupied MO space: CASE  [0,0,0]
c      *********************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO

       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       DO 3105 I=1,NOCC-2
c      ------------------
       M = NOCC+1-I
c      ------------
       NEL(M) = 0



       DO 3110 J=1,M-2
c      ------------------
       MM = M - J
c      ------------
       NEL(MM) = 0



       DO 3115 JJ=1,MM-1
c      ------------------
       MMM = MM - JJ
c      -------------
       NEL(MMM) = 0



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      --------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3150 
       ENDDO
       GO TO 5000
 3150  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3235
       M1 = L1 - 1
 3235  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3236
       LIM4 = LIM4 - 1
 3236  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3237
       M2 = L2 - 1
 3237  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3238 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3240 
       ENDDO
       GO TO 5000
 3240  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L2) = 0
       ENDDO
 3238  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3245
       M1 = L1 - 1
 3245  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3246
       M2 = L2 - 1
 3246  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3247
       M3 = L3 - 1
 3247  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3248
       M4 = L4 - 1
 3248  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3250
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3251
       ENDDO
       GO TO 5000
 3251  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3250  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3255
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3256
       ENDDO
       GO TO 5000
 3256  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3255  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
C*     -----------------------------------
c      **********************VIRTUAL

c      -----------------------------
c      occupied MOs: CASE [0,0,0]
c      -----------------------------
       NEL(MMM) = 2
 3115  CONTINUE


       NEL(MM) = 2
 3110  CONTINUE


       NEL(M) = 2
 3105  CONTINUE
c      ********

c---------------------------------------------------
c***** write (6,*) 'Sixtuple OCCU=[0,0,0] done'
c***** write (6,*) '---------------------------'







c      *******
 3400  A=0
c      *******




c      SCF-SPACE:
c***** write (6,*) 'Sixtuple OCCU[1,1,0,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,1,0,0]
c      ***********************************

       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO

       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC-1 
       DO 3500 I1=1,MAX1
c      ------------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 0



       MAX2 = N1-1
       DO 3505 I2=1,MAX2
c      -------------------
c      *
       N2 = N1 - I2
       NEL(N2) = 0



       MAX3 = NOCC-3
       DO 3510 I3=1,MAX3
c      ------------------
c      *
       NN3 = NOCC + 1 - I3
       IF (I3.EQ.1) GO TO 3511
       NN3 = LL3 - 1
 3511  A=0
       call SHIFT2(N1,N2,NN3,LL3)
       NEL(LL3) = 1




       MAX4 = NN3 - 1  
 3512  A=0
       DO 3515 I4=1,MAX4
c      ------------------
c      *
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 3516
       NN4 = LL4 - 1
 3516  A=0
       call SHIFT2(N1,N2,NN4,LL4)
       IF (LL4.EQ.0) GO TO 3517
       NEL(LL4) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      --------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3532
       ENDDO
       GO TO 5000
 3532  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3535
       M1 = L1 - 1
 3535  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3536
       LIM4 = LIM4 - 1
 3536  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3537
       M2 = L2 - 1
 3537  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3538 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3540
       ENDDO
       GO TO 5000
 3540  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L2) = 0
       ENDDO
 3538  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3545
       M1 = L1 - 1
 3545  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3546
       M2 = L2 - 1
 3546  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3547
       M3 = L3 - 1
 3547  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3548
       M4 = L4 - 1
 3548  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3550
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3551
       ENDDO
       GO TO 5000
 3551  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3550  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3555
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3556
       ENDDO
       GO TO 5000
 3556  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3555  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
C*     -----------------------------------

       NEL(LL4)=2
 3515  CONTINUE
 3517  A=0

 
       NEL(LL3)=2
 3510  CONTINUE


       NEL(N2)=2
 3505  CONTINUE


       NEL(N1)=2
 3500  CONTINUE
c      ***

c--------------------------------------------------
c***** write (6,*) 'Sixtuple OCCU[1,1,0,0] done'
c***** write (6,*) '---------------------------'





c      ********
 3600  A=0
c      ********


c      SCF-SPACE:
c****  write (6,*) 'Sixtuple OCCU[1,1,1,1,0] starts'
c      ***********************************
C      occupied MO space: CASE  [1,1,1,1,0]
c      ***********************************

       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO

       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO


       MAX1 = NOCC
       DO 3605 I1=1,NOCC
c      -----------------
c      *
       N1 = NOCC + 1 - I1
       NEL(N1) = 0


       MAX2 = NOCC-4
       DO 3610 I2=1,MAX2
c      ------------------
c      *
       NN2 = NOCC + 1 - I2
       IF (I2.EQ.1) GO TO 3611
       NN2 = LL2 - 1
 3611  A=0
       call SHIFT1(N1,NN2,LL2)
       NEL(LL2) = 1



       MAX3 = NN2 - 3
       DO 3615 I3=1,MAX3
c      -----------------
c      *
       NN3 = LL2 - I3
       IF (I3.EQ.1) GO TO 3616
       NN3 = LL3 - 1
 3616  A=0
       call SHIFT1(N1,NN3,LL3)
       NEL(LL3) = 1




       MAX4 = NN3 - 2
       DO 3620 I4=1,MAX4
c      -----------------
c      *
       NN4 = LL3 - I4
       IF (I4.EQ.1) GO TO 3621
       NN4 = LL4 - 1
 3621  A=0
       call SHIFT1(N1,NN4,LL4)
       NEL(LL4) = 1


       MAX5 = NN4 - 1
       DO 3625 I5=1,MAX5
c      ------------------
c      *
       NN5 = LL4 - I5
       IF (I5.EQ.1) GO TO 3626
       NN5 = LL5 - 1
 3626  A=0
       call SHIFT1(N1,NN5,LL5)
       IF (LL5.EQ.0) GO TO 3627
       NEL(LL5) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      --------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3632
       ENDDO
       GO TO 5000
 3632  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      write (6,*) 'Sixtuple VIRT [2,2,2] done'


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3635
       M1 = L1 - 1
 3635  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3636
       LIM4 = LIM4 - 1
 3636  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3637
       M2 = L2 - 1
 3637  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3638 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3640
       ENDDO
       GO TO 5000
 3640  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L2) = 0
       ENDDO
 3638  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3645
       M1 = L1 - 1
 3645  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3646
       M2 = L2 - 1
 3646  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3647
       M3 = L3 - 1
 3647  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3648
       M4 = L4 - 1
 3648  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3650
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3651
       ENDDO
       GO TO 5000
 3651  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3650  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3655
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3656
       ENDDO
       GO TO 5000
 3656  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3655  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *****
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
C*     -----------------------------------

       NEL(LL5)=2
 3625  CONTINUE
 3627  A=0
       NEL(LL4)=2 
 3620  CONTINUE
       NEL(LL3)=2
 3615  CONTINUE
       NEL(LL2)=2
 3610  CONTINUE
       NEL(N1)=2
 3605  CONTINUE
c      ******
c---------------------------------------------------------
c***** write (6,*) 'Sixtuple OCCU[1,1,1,1,0] done'
c***** write (6,*) '---------------------------'

 3800  A=0




c      SCF-SPACE:
c***** write (6,*) 'Sixtuple OCCU[1,1,1,1,1,1] starts'
c
c      **************************************
C      occupied MO space: CASE  [1,1,1,1,1,1]
c      **************************************
       DO I= 1,NOCC
       NEL(I) = 2
       ENDDO
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       MAX1 = NOCC-5
       DO 3805 I1=1,MAX1
c      -----------------
c      *
       N1 = NOCC + 1 - I1 
       NEL(N1) = 1

       MAX2 = N1 - 5 
       DO 3810 I2=1,MAX2
c      -----------------
c      *
       N2 = N1 - I2
       NEL(N2) = 1

       MAX3 = N2 - 4
       DO 3815 I3=1,MAX3
c      ------------------
c      *
       N3 = N2 - I3
       NEL(N3) = 1

       MAX4 = N3 - 3
       DO 3820 I4=1,MAX4
c      -----------------
c      *
       N4 = N3 - I4
       NEL(N4) = 1

       MAX5 = N4 - 2 
       DO 3825 I5=1,MAX5
c      -----------------
c      *
       N5 = N4 - I5
       NEL(N5) = 1

       MAX6 = N5 - 1
       DO 3830 I6=1,MAX6
c      ------------------
c      *
       N6 = N5 - I6
       IF (N6.GT.0) GO TO 3831 
       write (6,*) 'N6=',N6
       STOP
 3831  A=0
       NEL(N6) = 1


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      *******************
c      CASE) VIRT[2,2,2]
c      *******************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO
       LIM1 = NVIRT-2
       DO KK=1,LIM1
c      -------------------
c      *
       MK = NTOT+1-KK
c      --------------
       NEL(MK) = 2

       ING = KK-1 
       LIM2 = NVIRT-KK-1
       DO KKK=1,LIM2
c      ----------------------
c      *
       MMK = MK - KKK
c      ------------
       NEL(MMK) = 2

       LIM3 = NVIRT-KKK-1-ING
       DO KKK2=1,LIM3
c      -------------------------
c      *
       MMK2 = MMK - KKK2
c      -----------------
       NEL(MMK2) = 2
c      *************************
c      --------------------
       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3832
       ENDDO
       GO TO 5000
 3832  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(MMK2) = 0
       ENDDO

       NEL(MMK) = 0
       ENDDO

       NEL(MK) = 0
       ENDDO
c      ***********************************


c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) 'VIRT [2,2,1,1] starts'
c      *********************
c      CASE) VIRT [2,2,1,1]
c      *********************
c      -----------------------
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-1
       DO K1=1,LIM1
c      --------------------
c      *
       NX1 = NTOT + 1 - K1
       NEL(NX1) = 2 

       LIM2 = NX1-NOCC-1
       DO K2=1,LIM2
c      -------------------------
c      *
       NX2 = NX1 - K2
       NEL(NX2) = 2 

       LIM3 = NVIRT - 3
       DO KK1=1,LIM3
c      ----------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3835
       M1 = L1 - 1
 3835  A=0
       call SHIFT2(NX1,NX2,M1,L1)
       NEL(L1) = 1


       LIM4 = M1-NOCC-1
       IF (KK1.NE.1) GO TO 3836
       LIM4 = LIM4 - 1
 3836  A=0
       DO KK2=1,LIM4
c      ----------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3837
       M2 = L2 - 1
 3837  A=0
       call SHIFT2(NX1,NX2,M2,L2)
       IF (L2.EQ.NOCC) GO TO 3838 
       NEL(L2) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3840
       ENDDO
       GO TO 5000
 3840  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L2) = 0
       ENDDO
 3838  A=0
c      ***
       NEL(L1) = 0
       ENDDO

       NEL(NX2) = 0
       ENDDO

       NEL(NX1) = 0
       ENDDO
c      write (6,*) 'case Sixtuple[2,2,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [2,1,1,1,1] virt starts'
c      ***********************
c      CASE) VIRT [2,1,1,1,1]
c      ***********************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       DO K=1,NVIRT
c      -----------------
c      *
       NX = NTOT+1-K
       NEL(NX) = 2

       LIM1 = NVIRT-4
       DO KK1=1,LIM1
c      --------------------
c      *
       M1 = NTOT + 1 - KK1
       IF (KK1.EQ.1) GO TO 3845
       M1 = L1 - 1
 3845  A=0
       call SHIFT1(NX,M1,L1)
       NEL(L1) = 1



       LIM2 = M1-NOCC-3
       DO KK2=1,LIM2
c      ---------------------
c      *
       M2 = L1 - KK2
       IF (KK2.EQ.1) GO TO 3846
       M2 = L2 - 1
 3846  A=0
       call SHIFT1(NX,M2,L2)
       NEL(L2) = 1


       LIM3 = M2-NOCC-2 
       DO KK3=1,LIM3
c      ---------------------
c      *
       M3 = L2 - KK3
       IF (KK3.EQ.1) GO TO 3847
       M3 = L3 - 1
 3847  A=0
       call SHIFT1(NX,M3,L3)
       NEL(L3) = 1


       LIM4 = M3-NOCC-1
       DO KK4=1,LIM4
c      ---------------------
c      *
       M4 = L3 - KK4
       IF (KK4.EQ.1) GO TO 3848
       M4 = L4 - 1
 3848  A=0
       call SHIFT1(NX,M4,L4)
       IF (L4.EQ.NOCC) GO TO 3850
       NEL(L4) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3851
       ENDDO
       GO TO 5000
 3851  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(L4)=0
       ENDDO
 3850  A=0
c      ****

       NEL(L3)=0
       ENDDO

       NEL(L2)=0
       ENDDO

       NEL(L1)=0
       ENDDO

       NEL(NX)=0
       ENDDO
c      write (6,*) 'case Sixtuple[2,1,1,1,1] virt done'



c      ^^^^^^^^^^^^^^^^^^^^^^^
c      virtual space
c      ^^^^^^^^^^^^^^^^^^^^^^^
c      write (6,*) ' CASE) [1,1,1,1,1,1] virt starts'
c      ***************************
c      CASE) VIRT [1,1,1,1,1,1]
c      ***************************
       DO K = NOCC+1,NTOT
       NEL(K) = 0
       ENDDO

       LIM1 = NVIRT-5
       DO KK1=1,LIM1
c      -------------------
c      *
       M1 = NTOT + 1 - KK1
       NEL(M1) = 1


       LIM2 = M1-NOCC-5
       DO KK2=1,LIM2
c      -------------------
c      *
       M2 = M1 - KK2
       NEL(M2) = 1


       LIM3 = M2-NOCC-4
       DO KK3=1,LIM3
c      -------------------
c      *
       M3 = M2 - KK3
       NEL(M3) = 1


       LIM4 = M3-NOCC-3
       DO KK4=1,LIM4
c      ------------------
c      *
       M4 = M3 - KK4
       NEL(M4) = 1


       LIM5 = M4-NOCC-2
       DO KK5=1,LIM5
c      -------------------
c      *
       M5 = M4 - KK5
       NEL(M5) = 1


       LIM6 = M5-NOCC-1
       DO KK6=1,LIM6
c      ------------------
c      *
       M6 = M5 - KK6
       IF (M6.EQ.NOCC) GO TO 3855
       NEL(M6) = 1
C*     -----------------------------------

       ICOUNT = ICOUNT + 1
c      --------------------
c      *****************
       DO JAM=1,NTOT
c      -------------
       IF (ICON(JAM).NE.NEL(JAM)) GO TO 3856
       ENDDO
       GO TO 5000
 3856  A=0
c      *****************
c      ----------------------------------
C*     -----------------------------------
       NEL(M6)=0
       ENDDO
c      *****
 3855  A=0

       NEL(M5)=0
       ENDDO

       NEL(M4)=0
       ENDDO

       NEL(M3)=0
       ENDDO

       NEL(M2)=0
       ENDDO

       NEL(M1)=0
       ENDDO
c      *************************************
c      write (6,*) 'Sixtuple virt[1,1,1,1,1,1] done'
c      *************************************

       NEL(N6)=2
 3830  CONTINUE


       NEL(N5)=2
 3825  CONTINUE


       NEL(N4)=2
 3820  CONTINUE


       NEL(N3)=2
 3815  CONTINUE


       NEL(N2)=2
 3810  CONTINUE


       NEL(N1)=2
 3805  CONTINUE

     

c***** write (6,*) 'Sixtuple OCCU[1,1,1,1,1,1] done'
c***** write (6,*) '------------------------------------'


c      ----------------------------------------------
c***** write (6,*) 'subr. LABELSDTQ: 6- got here-all'
       GO TO 5000
c      ----------







c*******************END OF SDTQ56-excitations************
c                   -------------------------
c******LABELSDTQ:  SDTQ56-excitation-labeling ends ******




c      ******
 5000  AAA=0
c      ******
c      -----------------
       NSPRO = ICOUNT
c      -----------------
       IF (NSPRO.LT.mxsprod) GO TO 5001
c      ------------------------------
c**    write (6,*) 'subr. LABELSDTQ: exceeded limits'
c**    write (6,*) 'IEXCIT =',IEXCIT
c**    write (6,*) 'NSPRO  =',NSPRO
c**    write (6,*) 'mxsprod  =',mxsprod
c***** STOP
 5001  A=0
c      ***



c***************************************************
c       IF (NSPRO.NE.2146) GO TO 5002
c       write (6,*) 'IEXCIT=',IEXCIT
c       write (6,*) 'NSPRO=',NSPRO
c       write (6,*) (NEL(LIL), LIL=1,NTOT) 
c       write (6,*) 'subr. LABELSDTQ will stop'
c       STOP
c 5002  A=0
c       ****************************
c        IF (IEXCIT.NE.2) GO TO 6000
c        IF (NSPRO.NE.324) GO TO 6000
c        write (6,*) 'NSPRO=',NSPRO
c        write (6,*) (NEL(LIL), LIL=1,NTOT) 
c        write (6,*) 'subr. LABELSDTQ will stop'
c        STOP
c 6000   A=0
c      *****************
c
c      -----------------
c      subr LABELSDTQ ends
c      -----------------
       RETURN
       END











c---------------------------------------------------------
c     ***
      SUBROUTINE ESTQUAD(NTOT,NOCC,NVIRT,IEXCIT,
     *NDOUBLE,CISUMD,CID,NSCF,ISCF,IOCC,
     *NVIR,IVIR,IOCCV,CQEST,CQADD)
c     ----------------------------------------------------
      implicit double precision(a-h,o-z)
C
      parameter (mxcomb = 100)
      parameter (mxorb = 100)
      parameter (mxdia = 1000)
      parameter (mxneb = 1000)
      parameter (mxsprod = 2400000)
      parameter (mxsprod1 = 100000)

      DIMENSION ISCF(mxorb),IOCC(mxorb),IVIR(mxorb) 
      DIMENSION IOCCV(mxorb)
      DIMENSION CISUMD(mxsprod1) 
      DIMENSION CID(mxsprod1) 
      DIMENSION LB2(mxsprod1)
      DIMENSION CCPRD(mxcomb) 
      DIMENSION ISCFA(mxorb),IVIRA(mxorb) 
      DIMENSION ICON1(mxorb),ICON2(mxorb)
c
      DIMENSION IDOLD1(500),IDOLD2(500)
c----------------------------------------------------------
c     This code gives an estimate for the
c     Q-quadruple space-product
c     as the sum of D-double sp products
c     that are consistent with this Q-sp.
c
c     Given:
c     ------
c     CISUMD(j), j=1,NDOUBLE numerical values of D-sp
c     coefficients that have been obtained before.
c
c     ISCF, IOCC & IVIR,IOCCV : occupied and virtual
c     orbitals with their occupation numbers.
c
c     On return:
c     ----------
c     CQEST = the numerical estimate for the Q-sp.
c----------------------------------------------------------

      IADD = 0

      DO 25 N=1,NSCF
c     --------------
      IF (IOCC(N).LT.1) GO TO 30 
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
      GO TO 31
c     *
 30   A=0
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
 31   A=0
c     *
 25   CONTINUE


      IF (IADD.EQ.4) GO TO 35
      write (6,*) 'inconsistent IOCC'
      write (6,*) 'ISCFA(JON)=',(ISCFA(JON),JON=1,4)
      write (6,*) '--------------------------------'
      write (6,*) 'ISCF(i)=',(ISCF(JON),JON=1,NSCF)
      write (6,*) 'IOCC(i)=',(IOCC(JON),JON=1,NSCF)
      write (6,*) '--------------------------------'
      write (6,*) 'subr. ESTQUAD'
      STOP
 35   A=0


      IADD = 0
c     *
      DO 50 M=1,NVIR
c     ---------------
      IF (IOCCV(M).GT.1) GO TO 55
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
      GO TO 61
 55   A=0
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
 61   A=0
c     *
 50   CONTINUE


      IF (IADD.EQ.4) GO TO 65
c     ***
      write (6,*) 'inconsistent IOCCV'
      write (6,*) 'IVIRA(JON)=',(IVIRA(JON),JON=1,4)
      write (6,*) '--------------------------------'
      write (6,*) 'IVIR(i)=',(IVIR(JON),JON=1,NVIR)
      write (6,*) 'IOCCV(i)=',(IOCCV(JON),JON=1,NVIR)
      write (6,*) 'subr. ESTQUAD'
      STOP
 65   A=0



c     -----------------
c     ISCFA(j), j=1,4
c     IVIRA(k), k=1,4.
c     -----------------

      IEXCIT2 = 2

      CQEST = 0.0D0
      CQADD = 0.0D0

      ISUM = 0
      ITEM = 0
      TEM = 0.0D0

      DO JJJ=1,18
      CCPRD(JJJ) = 0.0D0
      ENDDO


c     *********
c     occupied:
c     *********
      I1 = ISCFA(1)
      I2 = ISCFA(2)
      I3 = ISCFA(3)
      I4 = ISCFA(4)
c     *******
c     virtual(1):
c     ------------
      K1 = IVIRA(1)
      K2 = IVIRA(2)
      K3 = IVIRA(3)
      K4 = IVIRA(4)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
c*      write (6,*) 'in subr. ESTQUAD, idnum CISUMD(idnum)'
c*      write (6,*) IDNUM1,CISUMD(IDNUM1)
c*      write (6,*) IDNUM2,CISUMD(IDNUM2)
c*      write (6,*) '***'
c     *
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
c     virtual(2):
c     ------------
      K1 = IVIRA(3)
      K2 = IVIRA(4)
      K3 = IVIRA(1)
      K4 = IVIRA(2)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
c     ***
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 100
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 100  A=0
c     --------------------------------
c     *******
c     virtual(3):
c     -----------
      K1 = IVIRA(1)
      K2 = IVIRA(3)
      K3 = IVIRA(2)
      K4 = IVIRA(4)
c     *
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
c     ***
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 101
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 101  A=0
c     --------------------------------
c     virtual(4):
c     -----------
      K1 = IVIRA(2)
      K2 = IVIRA(4)
      K3 = IVIRA(1)
      K4 = IVIRA(3)
c     *
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
c     ***
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 102
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 102  A=0
c     --------------------------------
c     *******
c     virtual(5):
c     ------------
      K1 = IVIRA(1)
      K2 = IVIRA(4)
      K3 = IVIRA(2)
      K4 = IVIRA(3)
c     *
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
c     ***
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 103
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 103  A=0
c     -------
c     virtual(6):
c     ------------
      K1 = IVIRA(2)
      K2 = IVIRA(3)
      K3 = IVIRA(1)
      K4 = IVIRA(4)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
c     ***
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 104
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 104  A=0
c     **********
      


c     *********
c     occupied:
c     *********
      I1 = ISCFA(1)
      I2 = ISCFA(3)
      I3 = ISCFA(2)
      I4 = ISCFA(4)
c     **********
c     virtual(7):
c     -----------
      K1 = IVIRA(1)
      K2 = IVIRA(2)
      K3 = IVIRA(3)
      K4 = IVIRA(4)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 105
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 105  A=0
c     --------------------------------
c     virtual(8):
c     -----------
      K1 = IVIRA(3)
      K2 = IVIRA(4)
      K3 = IVIRA(1)
      K4 = IVIRA(2)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
c     ***
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 106
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 106  A=0
c     --------------------------------
c     *******
c     virtual(9):
c     ------------
      K1 = IVIRA(1)
      K2 = IVIRA(3)
      K3 = IVIRA(2)
      K4 = IVIRA(4)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 107
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 107  A=0
c     --------------------------------
c     virtual(10):
c     ------------
      K1 = IVIRA(2)
      K2 = IVIRA(4)
      K3 = IVIRA(1)
      K4 = IVIRA(3)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 108
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 108  A=0
c     --------------------------------
c     virtual(11):
c     ------------
      K1 = IVIRA(1)
      K2 = IVIRA(4)
      K3 = IVIRA(2)
      K4 = IVIRA(3)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 109 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 109  A=0
c     --------------------------------
c     virtual(12):
c     ------------
      K1 = IVIRA(2)
      K2 = IVIRA(3)
      K3 = IVIRA(1)
      K4 = IVIRA(4)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 110 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 110  A=0
c     --------------------------------
c     **********
      
      
c     *********
c     occupied:
c     -------
      I1 = ISCFA(1)
      I2 = ISCFA(4)
      I3 = ISCFA(2)
      I4 = ISCFA(3)
c     **********
c     virtual(13):
c     -------------
      K1 = IVIRA(1)
      K2 = IVIRA(2)
      K3 = IVIRA(3)
      K4 = IVIRA(4)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 111 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 111  A=0
c     --------------------------------
c     virtual(14):
c     -------------
      K1 = IVIRA(3)
      K2 = IVIRA(4)
      K3 = IVIRA(1)
      K4 = IVIRA(2)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 112 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 112  A=0
c     --------------------------------
c     virtual(15):
c     -------------
      K1 = IVIRA(1)
      K2 = IVIRA(3)
      K3 = IVIRA(2)
      K4 = IVIRA(4)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 113 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 113  A=0
c     --------------------------------
c     virtual(16):
c     -------------
      K1 = IVIRA(2)
      K2 = IVIRA(4)
      K3 = IVIRA(1)
      K4 = IVIRA(3)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 114 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 114  A=0
c     --------------------------------
c     virtual(17):
c     ------------
      K1 = IVIRA(1)
      K2 = IVIRA(4)
      K3 = IVIRA(2)
      K4 = IVIRA(3)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 115 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 115  A=0
c     --------------------------------
c     virtual(18):
c     ------------
      K1 = IVIRA(2)
      K2 = IVIRA(3)
      K3 = IVIRA(1)
      K4 = IVIRA(4)
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I3,I4,K3,K4,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 116 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQEST=CQEST+(SGN*CISUMD(IDNUM1)*CISUMD(IDNUM2))
      CQADD=CQADD+(SGN*CID(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------
 116  A=0
c     --------------------------------
c
c     *
c     MAXIMUM product used:
c     -------------------
c*    NN = ITEM
c*    CALL SORTFAST(NN,CCPRD)
c     -------------------
c*    CQEST = CCPRD(NN)
c     *

      IF (ITEM.GT.0) GO TO 200
      write (6,*) 'ESTQUAD: ITEM=',ITEM
      STOP
 200  A=0


c
c     -----------------
      CQEST = CQEST/TEM
c     -----------------
c     ---------------------
      CQEST = DABS(CQEST)
c     --------------------- 
c     ---------------------
      CQEST = SQRT(CQEST)
c     --------------------- 
c     *
c      write (6,*) 'subr. ESTQUAD'
c      write (6,*) 'I1,I2,I3,I4;K1,K2,K3,K4'
c      write (6,*) I1,I2,I3,I4,K1,K2,K3,K4
c      write (6,*) 'CQEST=',CQEST 
c      write (6,*) 'CQADD=',CQADD 
c      write (6,*) 'total orb NTOT=',NTOT
c      write (6,*) 'occ orb NOCC=',NOCC
c      write (6,*) 'virt orb NVIRT=',NVIRT
c      do J=1,ISUM
c      write (6,*) 'IDOLD1 2=',IDOLD1(J),IDOLD2(J)
c      enddo
c     *      
c     ----------
      RETURN
      END








c---------------------------------------------------------
c     ***
      SUBROUTINE EST555(NTOT,NOCC,NVIRT,IEXCIT,
     *NDOUBLE,CISUMD,CID,NTRIPLE,CI3,CI3sum,
     *NSCF,ISCF,IOCC,NVIR,IVIR,IOCCV,C5EST,C5ADD)
c     ----------------------------------------------------
      implicit double precision(a-h,o-z)
C
      parameter (mxorb = 100)
      parameter (mxdia = 1000)
      parameter (mxneb = 1000)
      parameter (mxsprod = 2400000)
      parameter (mxsprod1 = 100000)

      DIMENSION ISCF(mxorb),IOCC(mxorb),IVIR(mxorb) 
      DIMENSION IOCCV(mxorb)

      DIMENSION CISUMD(mxsprod1) 
      DIMENSION CID(mxsprod1) 
      DIMENSION LB2(mxsprod1)

      DIMENSION CI3(mxsprod1) 
      DIMENSION CI3sum(mxsprod1)

      DIMENSION ISCFA(mxorb),IVIRA(mxorb) 
      DIMENSION ICON1(mxorb),ICON2(mxorb)
c
      DIMENSION IDOLD1(500),IDOLD2(500)
c
c----------------------------------------------------------
c     This code gives an estimate for the
c     5-pentuple space-product
c     as the sum of D-double*T-triple sp products
c     that are consistent with this 5-sp.
c
c     Given:
c     ------
c     CISUMD(j), j=1,NDOUBLE numerical values of D-sp
c     coefficients that have been obtained before.
c
c     ISCF, IOCC & IVIR,IOCCV : occupied and virtual
c     orbitals with their occupation numbers.
c
c     On return:
c     ----------
c     C5EST = the numerical estimate for the 5-sp.
c     C5ADD = the norm(c**2) estimate for the 5-sp.
c----------------------------------------------------------





      ICYC = 0

      IADD = 0



      DO 25 N=1,NSCF
c     --------------
      IF (IOCC(N).LT.1) GO TO 30 
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
      GO TO 31
c     *
 30   A=0
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
 31   A=0
c     *
 25   CONTINUE


      IF (IADD.EQ.5) GO TO 35
c     -----------------------
      write (6,*) 'inconsistent IOCC'
      write (6,*) 'ISCFA(JON)=',(ISCFA(JON),JON=1,5)
      write (6,*) '--------------------------------'
      write (6,*) 'ISCF(i)=',(ISCF(JON),JON=1,NSCF)
      write (6,*) 'IOCC(i)=',(IOCC(JON),JON=1,NSCF)
      write (6,*) '--------------------------------'
      write (6,*) 'subr. EST555'
      STOP
 35   A=0


      IADD = 0
c     *
      DO 50 M=1,NVIR
c     ---------------
      IF (IOCCV(M).GT.1) GO TO 55
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
      GO TO 61
 55   A=0
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
 61   A=0
c     *
 50   CONTINUE


      IF (IADD.EQ.5) GO TO 65
c     -----------------------
c     ***
      write (6,*) 'inconsistent IOCCV'
      write (6,*) 'IVIRA(JON)=',(IVIRA(JON),JON=1,5)
      write (6,*) '--------------------------------'
      write (6,*) 'IVIR(i)=',(IVIR(JON),JON=1,NVIR)
      write (6,*) 'IOCCV(i)=',(IOCCV(JON),JON=1,NVIR)
      write (6,*) 'subr. EST555'
      STOP
 65   A=0




c     -----------------
c     ISCFA(j), j=1,5
c     IVIRA(k), k=1,5.
c     -----------------




c**** ITESTING = 1
c     ************


      IEXCIT2 = 2
      IEXCIT3 = 3 



      C5EST = 0.0D0
      C5ADD = 0.0D0



      ISUM = 0
      ITEM = 0
      TEM = 0.0D0






c     ************
      DO 66 M1=1,3
c     *
c     ------------
      M2L = M1 + 1



      DO 67 M2=M2L,4
c     *
c     --------------
      M3L = M2 + 1



      DO 68  M3=M3L,5
c     *
c     --------------







      If (M1.EQ.1) GO TO 70
      If (M1.EQ.2) GO TO 90
      If (M1.EQ.3) GO TO 110
      write (6,*) 'subr. EST555 stop'
      STOP
c     ---------------------- 


c     ------
c     M1 = 1
c     ------
c     *****
 70   A=0.0
c     *****
      IF (M2.LT.3) GO TO 72
      M4 = 2
c     ------
      IF (M2.LT.4) GO TO 80
      M5 = 3
c     -----------
c     [1,4,5, 2,3]
c     -----------
      GO TO 500

c---------------------------
c     *******(M2=3)
 80   A=0.0
c     *******
      IF (M3.EQ.5) GO TO 84
      M5 = 5
c     -----------
c     [1,3,4, 2,5]
c     -----------
      GO TO 500


 84   A=0.0
c     ********(M3=5)
      M5 = 4
c     -----------
c     [1,3,5, 2,4]
c     -----------
      GO TO 500


c---------------------------
c     ***** (M2=2)
 72   A=0.0
c     *****
      IF (M3.LT.4) GO TO 74
      M4 = 3
c     ------ 
      IF (M3.EQ.5) GO TO 76
      M5 = 5
c     ------
c     [1,2,4, 3,5]
c     -----------
      GO TO 500



 76   A=0.0
c     *******
      M5 = 4
c     ------
c     [1,2,5, 3,4]
c     -----------
      GO TO 500



 74   A=0.0
c     ***** (M3=3)
      M4 = 4
      M5 = 5
c     ------
c     [1,2,3, 4,5]
c     ------------
      GO TO 500


c---------------------------

c     ------
c     M1 = 2 
c     ------
c     *****
 90   A=0.0
c     *****
      M4 = 1
c     ------
      IF (M2.LT.4) GO TO 96
      M5 = 3
c     ------
c     [2,4,5, 1,3]
c     ------------
      GO TO 500



 96   A=0.0
c     *****
      IF (M3.LT.5) GO TO 98
      M5 = 4
c     ------
c     [2,3,5, 1,4]
c     -----------
      GO TO 500



 98   A=0.0
c     ******
      M5 = 5
c     ------
c     [2,3,4, 1,5]
c     ------------
      GO TO 500



c------------------------------


c     ------
c     M1 = 3 
c     ------
c     *****
 110  A=0.0
c     *****
c     --------
      M4 = 1
      M5 = 2
c     --------
c     [3,4,5, 1,2]
c     -----------


c     *******
 500  A=0.0
c     *******


c***** OCCUPIED MO SPACE DONE **********






c     *********
c     occupied:
c     *********
c     ---------------
      I1 = ISCFA(M1)
      I2 = ISCFA(M2)
      I3 = ISCFA(M3)
      I4 = ISCFA(M4)
      I5 = ISCFA(M5)
c     ---------------







c     VIRTUAL MO SPACE:
c     *****************

      DO 550 MM1=1,3
c     *
c     ------------------
      MM2L = MM1 + 1


      DO 555 MM2=MM2L,4
c     *
c     ------------------
      MM3L = MM2 + 1


      DO 560  MM3=MM3L,5
c     *
c     -------------------





      If (MM1.EQ.1) GO TO 570
      If (MM1.EQ.2) GO TO 590
      If (MM1.EQ.3) GO TO 650 
      write (6,*) 'subr. EST555 stop'
      STOP
c     ---------------------- 


c     ------
c     MM1 = 1
c     ------
c     *****
 570  A=0.0
c     *****
      IF (MM2.LT.3) GO TO 572
      MM4 = 2
c     ------
      IF (MM2.LT.4) GO TO 580
      MM5 = 3
c     -----------
c     [1,4,5, 2,3]
c     -----------
      GO TO 700


c---------------------------
c     *******(MM2=3)
 580  A=0.0
c     *******
      IF (MM3.EQ.5) GO TO 584
      MM5 = 5
c     -----------
c     [1,3,4, 2,5]
c     -----------
      GO TO 700



 584  A=0.0
c     ********(MM3=5)
      MM5 = 4
c     -----------
c     [1,3,5, 2,4]
c     -----------
      GO TO 700


c---------------------------
c     ***** (MM2=2)
 572  A=0.0
c     *****
      IF (MM3.LT.4) GO TO 574
      MM4 = 3
c     ------ 
      IF (MM3.EQ.5) GO TO 576
      MM5 = 5
c     ------
c     [1,2,4, 3,5]
c     -----------
      GO TO 700



 576  A=0.0
c     *******
      MM5 = 4
c     ------
c     [1,2,5, 3,4]
c     -----------
      GO TO 700



 574  A=0.0
c     ***** (MM3=3)
      MM4 = 4
      MM5 = 5
c     ------
c     [1,2,3, 4,5]
c     ------------
      GO TO 700


c---------------------------



c     ------
c     MM1 = 2 
c     ------
c     *****
 590  A=0.0
c     *****

      MM4 = 1
c     ------
      IF (MM2.LT.4) GO TO 596
      MM5 = 3
c     ------
c     [2,4,5, 1,3]
c     ------------
      GO TO 700



 596  A=0.0
c     *****
      IF (MM3.LT.5) GO TO 598
      MM5 = 4
c     ------
c     [2,3,5, 1,4]
c     -----------
      GO TO 700



 598  A=0.0
c     ******
      MM5 = 5
c     ------
c     [2,3,4, 1,5]
c     ------------
      GO TO 700


c------------------------------

c     ------
c     MM1 = 3 
c     ------
c     *****
 650  A=0.0
c     *****
c     --------
      MM4 = 1
      MM5 = 2
c     --------
c     [3,4,5, 1,2]
c     -----------



c     *******
 700  A=0.0
c     *******



c***** VIRTUAL MO SPACE DONE **********




c     ---------------
      ICYC = ICYC + 1
c     ---------------



c     *******
c     virtual:
c     --------------
      K1 = IVIRA(MM1)
      K2 = IVIRA(MM2)
      K3 = IVIRA(MM3)
      K4 = IVIRA(MM4)
      K5 = IVIRA(MM5)
c     ---------------




c**  write (6,1112) I1,I2,I3,I4,I5,K1,K2,K3,K4,K5
c** 1112 FORMAT (I5,I5,I5,I5,I5,4X,I5,I5,I5,I5,I5)





c**** IF (ITESTING.EQ.1) GO TO 1000
c     *****************************


c     *
      CALL NUMTRIPLE(NTOT,NOCC,NVIRT,I1,I2,I3,K1,K2,K3,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT3,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I4,I5,K4,K5,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *

c***  write (6,*) 'TEST: subr. EST555, idnum CISUMD(idnum)'
c***  write (6,*) 'No   CI3(#)=',IDNUM1,CI3(IDNUM1)
c***  write (6,*) 'No   CIdouble(#)=',IDNUM2,CISUMD(IDNUM2)
c***  write (6,*) 'will stop ***'
c     *
c     ------------------------
      IF (ICYC.EQ.1) GO TO 999
c     ------------------------
c     *
c     ***
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
c     ---------------
      ISUM = ISUM + 1
c     ---------------
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 1000
c     ***
 999  A=0.0
c     ***
      SGN = 1.0D0
c     -----------
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
c     ***
c     -----------------------------------------------
      C5EST = C5EST + (CI3(IDNUM1)*CISUMD(IDNUM2))

      C5ADD = C5ADD + (CI3sum(IDNUM1)*CID(IDNUM2))
c     -----------------------------------------------


c     ***
 1000 A=0
c     ***





 560  CONTINUE
c     --------
 555  CONTINUE
c     --------
 550  CONTINUE
c     --------




 68   CONTINUE
c     --------
 67   CONTINUE
c     --------
 66   CONTINUE
c     --------








c*************************
c     *
c     MAXIMUM product used:
c     -------------------
c*    NN = ITEM
c*    CALL SORTFAST(NN,CCPRD)
c*    CQEST = CCPRD(NN)
c*************************

      IF (ITEM.GT.0) GO TO 2000
      write (6,*) 'stop EST555: ITEM=',ITEM
      STOP
 2000 A=0


c
c     -----------------
      C5EST = C5EST/TEM
c     -----------------
c     ---------------------
      C5EST = DABS(C5EST)
c     --------------------- 
c     ---------------------
      C5EST = SQRT(C5EST)
c     --------------------- 
c     *
c**** write (6,*) 'subr. EST555  for 5-tuples'
c***  write (6,*) 'I1,I2,I3,I4,I5;K1,K2,K3,K4,K5'
c***  write (6,*) I1,I2,I3,I4,I5,K1,K2,K3,K4,K5
c***  write (6,*) 'C5EST=',C5EST 
c***  write (6,*) 'C5ADD=',C5ADD 
c***  write (6,*) 'total orb NTOT=',NTOT
c***  write (6,*) 'occ orb NOCC=',NOCC
c***  write (6,*) 'virt orb NVIRT=',NVIRT
c***  write (6,*) '# of distinct comb. ITEM=',ITEM
c***  do J=1,ISUM
c***  write (6,*) 'IDOLD1 2=',IDOLD1(J),IDOLD2(J)
c***  enddo
c***  write (6,*) 'subr. EST555: ICYC=',ICYC
c***  write (6,*) 'EST555: will stop'
c     *      
c
c     ----------
      RETURN
      END







c---------------------------------------------------------
c     ***
      SUBROUTINE EST666(NTOT,NOCC,NVIRT,IEXCIT,
     *NDOUBLE,CISUMD,CID,NQUADRO,CI4,CI4sum,LB4,
     *NSCF,ISCF,IOCC,NVIR,IVIR,IOCCV,C6EST,C6ADD)
c     ----------------------------------------------------
      implicit double precision(a-h,o-z)
C
      parameter (mxorb = 100)
      parameter (mxdia = 1000)
      parameter (mxneb = 1000)

      parameter (mxsprod = 2400000)
      parameter (mxsprod1 = 100000)

      DIMENSION ISCF(mxorb),IOCC(mxorb),IVIR(mxorb) 
      DIMENSION IOCCV(mxorb)

      DIMENSION CISUMD(mxsprod1) 
      DIMENSION CID(mxsprod1) 
      DIMENSION LB2(mxsprod1)

      DIMENSION CI4(mxsprod) 
      DIMENSION CI4sum(mxsprod)
      DIMENSION LB4(mxsprod)

      DIMENSION ISCFA(mxorb),IVIRA(mxorb) 
      DIMENSION ICON1(mxorb),ICON2(mxorb)
c
      DIMENSION IDOLD1(500),IDOLD2(500)
c
c----------------------------------------------------------
c     This code gives an estimate for the
c     5-pentuple space-product
c     as the sum of D-double*T-triple sp products
c     that are consistent with this 5-sp.
c
c     Given:
c     ------
c     CISUMD(j), j=1,NDOUBLE numerical values of D-sp
c     coefficients that have been obtained before.
c
c     ISCF, IOCC & IVIR,IOCCV : occupied and virtual
c     orbitals with their occupation numbers.
c
c     On return:
c     ----------
c     C6EST = the numerical estimate for the 5-sp.
c     C6ADD = the norm(c**2) estimate for the 5-sp.
c----------------------------------------------------------





      ICYC = 0

      IADD = 0



      DO 25 N=1,NSCF
c     --------------
      IF (IOCC(N).LT.1) GO TO 30 
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
      GO TO 31
c     *
 30   A=0
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
 31   A=0
c     *
 25   CONTINUE


      IF (IADD.EQ.6) GO TO 35
c     -----------------------
      write (6,*) 'inconsistent IOCC'
      write (6,*) 'ISCFA(JON)=',(ISCFA(JON),JON=1,6)
      write (6,*) '--------------------------------'
      write (6,*) 'ISCF(i)=',(ISCF(JON),JON=1,NSCF)
      write (6,*) 'IOCC(i)=',(IOCC(JON),JON=1,NSCF)
      write (6,*) '--------------------------------'
      write (6,*) 'subr. EST666'
      STOP
 35   A=0


      IADD = 0
c     *
      DO 50 M=1,NVIR
c     ---------------
      IF (IOCCV(M).GT.1) GO TO 55
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
      GO TO 61
 55   A=0
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
 61   A=0
c     *
 50   CONTINUE


      IF (IADD.EQ.6) GO TO 65
c     -----------------------
c     ***
      write (6,*) 'inconsistent IOCCV'
      write (6,*) 'IVIRA(JON)=',(IVIRA(JON),JON=1,6)
      write (6,*) '--------------------------------'
      write (6,*) 'IVIR(i)=',(IVIR(JON),JON=1,NVIR)
      write (6,*) 'IOCCV(i)=',(IOCCV(JON),JON=1,NVIR)
      write (6,*) 'subr. EST666'
      STOP
 65   A=0




c     -----------------
c     ISCFA(j), j=1,6
c     IVIRA(k), k=1,6.
c     -----------------

c*    ITESTING = 1
c     ************


      IEXCIT2 = 2
      IEXCIT4 = 4 



      C6EST = 0.0D0
      C6ADD = 0.0D0



      ISUM = 0
      ITEM = 0
      TEM = 0.0D0






c     ************
      DO 66 M1=1,3
c     *
c     ------------
      M2L = M1 + 1



      DO 67 M2=M2L,4
c     *
c     --------------
      M3L = M2 + 1



      DO 68  M3=M3L,5
c     *
c     --------------
      M4L = M3 + 1




      DO 69  M4=M4L,6
c     *
c     ---------------







      If (M1.EQ.1) GO TO 70
      If (M1.EQ.2) GO TO 90
      If (M1.EQ.3) GO TO 110
      write (6,*) 'subr. EST666 stop'
      STOP
c     ------------------------------- 




c     ------
c     M1 = 1
c     ------
c     *****
 70   A=0.0
c     *****
      IF (M2.LT.3) GO TO 72
      M5 = 2
c     ------
      IF (M2.LT.4) GO TO 80
      M6 = 3
c     -------------
c     [1,4,5,6, 2,3]
c     --------------
      GO TO 500

c---------------------------
c     *******(M2=3)
 80   A=0.0
c     *******
      IF (M3.EQ.5) GO TO 84
c     *******(M3=4)
      IF (M4.LT.6) GO TO 82 
      M6 = 5
c     --------------
c     [1,3,4,6, 2,5]
c*    --------------
      GO TO 500


 82   A=0
c     ***
      M6 = 6
c     --------------
c     [1,3,4,5, 2,6]
c     --------------
      GO TO 500      


 84   A=0.0
c     ********(M3=5)
      M6 = 4
c     --------------
c     [1,3,5,6, 2,4]
c     --------------
      GO TO 500

c---------------------------
c---------------------------

c     ***** (M2=2)
 72   A=0.0
c     *****
      IF (M3.LT.4) GO TO 74
      M5 = 3
c     ------ 
      IF (M3.EQ.5) GO TO 76
c     *****(M3=4)
      IF (M4.EQ.6) GO TO 78
c     *****(M4=5)
      M6 = 6
c     -------------
c     [1,2,4,5, 3,6]
c     -------------
      GO TO 500

 78   A=0
c     ***
      M6 = 5
c     ----------------
c     [1,2,4,6, 3,5]
c     ----------------
      GO TO 500


 76   A=0.0
c     *******
      M6 = 4
c     ----------------
c     [1,2,5,6, 3,4]
c     ----------------
      GO TO 500


 74   A=0.0
c     ***** (M3=3)
      IF (M4.LT.5) GO TO 81
c     
      IF (M4.LT.6) GO TO 83
c     *****(M4=6)
      M5 = 4
      M6 = 5
c     ----------------
c     [1,2,3,6, 4,5]
c     ----------------
      GO TO 500

c
 81   A=0
c     *****(M4=4)
      M5 = 5
      M6 = 6
c     --------------
c     [1,2,3,4, 5,6]
c     -------------- 
      GO TO 500
c

 83   A=0
c     ***(M4=5)
      M5 = 4
      M6 = 6
c     ---------------
c     [1,2,3,5, 4,6]
c     ---------------
      GO TO 500


c---------------------------
c---------------------------
c---------------------------

c     ------
c     M1 = 2 
c     ------
c     *****
 90   A=0.0
c     *****
c     ------
      M5 = 1
c     ------
      IF (M2.LT.4) GO TO 96
c     ******(M2=4)
      M6 = 3
c     ---------------
c     [2,4,5,6, 1,3]
c     ---------------
      GO TO 500



 96   A=0.0
c     *****(M2=3)
      IF (M3.LT.5) GO TO 98
      M6 = 4
c     ---------------
c     [2,3,5,6, 1,4]
c     ---------------
      GO TO 500



 98   A=0.0
c     ******(M3=4)
      IF (M4.LT.6) GO TO 99
c     ****** (M4=6)
      M6 = 5
c     ---------------
c     [2,3,4,6, 1,5]
c     ---------------
      GO TO 500



 99   A=0
c     ******(M4=5)
      M6 = 6
c     ---------------
c     [2,3,4,5, 1,6]
c     ---------------
      GO TO 500


c------------------------------
c------------------------------
c------------------------------

c     ------
c     M1 = 3 
c     ------
c     *****
 110  A=0.0
c     *****
c     --------
      M5 = 1
      M6 = 2
c     --------------
c     [3,4,5,6, 1,2]
c     --------------


c     *******
 500  A=0.0
c     *******


c***** OCCUPIED MO SPACE DONE **********






c     *********
c     occupied:
c     *********
c     ---------------
      I1 = ISCFA(M1)
      I2 = ISCFA(M2)
      I3 = ISCFA(M3)
      I4 = ISCFA(M4)
      I5 = ISCFA(M5)
      I6 = ISCFA(M6)
c     ---------------
c     subr. EST 666:
c     ---------------






c     VIRTUAL MO SPACE:
c     *****************

      DO 550 MM1=1,3
c     *
c     ------------------
      MM2L = MM1 + 1


      DO 555 MM2=MM2L,4
c     *
c     ------------------
      MM3L = MM2 + 1


      DO 560  MM3=MM3L,5
c     *
c     -------------------
      MM4L = MM3 + 1


      DO 565  MM4=MM4L,6
c     *
c     -------------------





      If (MM1.EQ.1) GO TO 570
      If (MM1.EQ.2) GO TO 590
      If (MM1.EQ.3) GO TO 650 
      write (6,*) 'subr. EST666 stop'
      STOP
c     ------------------------------- 


c     ------
c     MM1 = 1
c     ------
c     *****
 570  A=0.0
c     *****
      IF (MM2.LT.3) GO TO 572
      MM5 = 2
c     ------
      IF (MM2.LT.4) GO TO 580
      MM6 = 3
c     -------------
c     [1,4,5,6, 2,3]
c     --------------
      GO TO 700


c---------------------------
c     *******(MM2=3)
 580  A=0.0
c     *******
      IF (MM3.EQ.5) GO TO 584
c     *******(MM3=4)
      IF (MM4.LT.6) GO TO 582 
      MM6 = 5
c     --------------
c     [1,3,4,6, 2,5]
c     --------------
      GO TO 700


 582  A=0
c     ***
      MM6 = 6
c     --------------
c     [1,3,4,5, 2,6]
c     --------------
      GO TO 700      


 584  A=0.0
c     ********(MM3=5)
      MM6 = 4
c     --------------
c     [1,3,5,6, 2,4]
c     --------------
      GO TO 700


c---------------------------
c---------------------------

c     ***** (MM2=2)
 572  A=0.0
c     *****
      IF (MM3.LT.4) GO TO 574
      MM5 = 3
c     ------ 
      IF (MM3.EQ.5) GO TO 576
c     *****(MM3=4)
      IF (MM4.EQ.6) GO TO 578
c     *****(MM4=5)
      MM6 = 6
c     -------------
c     [1,2,4,5, 3,6]
c     -------------
      GO TO 700



 578  A=0
c     ***
      MM6 = 5
c     ----------------
c     [1,2,4,6, 3,5]
c     ----------------
      GO TO 700



 576  A=0.0
c     *******
      MM6 = 4
c     ----------------
c     [1,2,5,6, 3,4]
c     ----------------
      GO TO 700



 574  A=0.0
c     ***** (MM3=3)
      IF (MM4.LT.5) GO TO 581
c     
      IF (MM4.LT.6) GO TO 583
c     *****(MM4=6)
      MM5 = 4
      MM6 = 5
c     ----------------
c     [1,2,3,6, 4,5]
c     ----------------
      GO TO 700


c
 581  A=0
c     *****(MM4=4)
      MM5 = 5
      MM6 = 6
c     --------------
c     [1,2,3,4, 5,6]
c     -------------- 
      GO TO 700



 583  A=0
c     ***(MM4=5)
      MM5 = 4
      MM6 = 6
c     ---------------
c     [1,2,3,5, 4,6]
c     ---------------
      GO TO 700



c---------------------------
c---------------------------
c---------------------------


c     ------
c     MM1 = 2 
c     ------
c     *****
 590  A=0.0
c     *****
c     ------
      MM5 = 1
c     ------
      IF (MM2.LT.4) GO TO 596
c     ******(MM2=4)
      MM6 = 3
c     ---------------
c     [2,4,5,6, 1,3]
c     ---------------
      GO TO 700



 596  A=0.0
c     *****(MM2=3)
      IF (MM3.LT.5) GO TO 598
      MM6 = 4
c     ---------------
c     [2,3,5,6, 1,4]
c     ---------------
      GO TO 700



 598  A=0.0
c     ******(MM3=4)
      IF (MM4.LT.6) GO TO 599
c     ****** (MM4=6)
      MM6 = 5
c     ---------------
c     [2,3,4,6, 1,5]
c     ---------------
      GO TO 700



 599  A=0
c     ******(MM4=5)
      MM6 = 6
c     ---------------
c     [2,3,4,5, 1,6]
c     ---------------
      GO TO 700



c------------------------------
c------------------------------
c------------------------------

c     ------
c     MM1 = 3 
c     ------
c     *****
 650  A=0.0
c     *****
c     --------
      MM5 = 1
      MM6 = 2
c     --------------
c     [3,4,5,6, 1,2]
c     --------------




c     *******
 700  A=0.0
c     *******




c***** VIRTUAL MO SPACE DONE **********




c     ---------------
      ICYC = ICYC + 1
c     ---------------



c     *******
c     virtual:
c     --------------
      K1 = IVIRA(MM1)
      K2 = IVIRA(MM2)
      K3 = IVIRA(MM3)
      K4 = IVIRA(MM4)
      K5 = IVIRA(MM5)
      K6 = IVIRA(MM6)
c     ---------------




C**  write (6,1112) I1,I2,I3,I4,I5,I6,K1,K2,K3,K4,K5,K6
c* 1112 FORMAT (I5,I5,I5,I5,I5,I5,4X,I5,I5,I5,I5,I5,I5)

c***      write (6,1112) M1,M2,M3,M4,M5,M6,MM1,MM2,MM3,
c***     *MM4,MM5,MM6
c*** 1112 FORMAT (I5,I5,I5,I5,I5,I5,4X,I5,I5,I5,I5,I5,I5)




c*****IF (ITESTING.EQ.1) GO TO 1000
c     *****************************


c     *
      CALL NUMQUADRO(NTOT,NOCC,NVIRT,I1,I2,I3,I4,
     *K1,K2,K3,K4,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT4,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I5,I6,K5,K6,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *

c***  write (6,*) 'TEST: subr. EST666, idnum CISUMD(idnum)'
c***  write (6,*) 'No   CI4(#)=',IDNUM1,CI4(IDNUM1)
c***  write (6,*) 'No   CIdouble(#)=',IDNUM2,CISUMD(IDNUM2)
c***  write (6,*) 'will stop ***'
c     *
c     ------------------------
      IF (ICYC.EQ.1) GO TO 999
c     ------------------------
c     *
c     ***
      CALL CHECKSAME(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
c     ---------------
      ISUM = ISUM + 1
c     ---------------
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 1000
c     ***
 999  A=0.0
c     ***
      SGN = 1.0D0
c     -----------
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
c     ***
c     ********************************************
      DO JACS=1,NQUADRO
c     -----------------
      IHIT1 = JACS
      ISETI = LB4(JACS)
      IF (IDNUM1.NE.ISETI) GO TO 1222
      GO TO 1225
 1222 A=0
      ENDDO


c     
c     -----------------------
      CI4(IHIT1) = 0.0D0
      CI4sum(IHIT1) = 0.0D0
c     -----------------------


 1225 A=0
c     *********************************************
c**** C6EST = C6EST + (CI4(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------

      C6EST = C6EST + (CI4sum(IHIT1)*CISUMD(IDNUM2))

      C6ADD = C6ADD + (CI4sum(IHIT1)*CID(IDNUM2))
c     -----------------------------------------------


c     ***
 1000 A=0
c     ***





 565  CONTINUE
c     --------
 560  CONTINUE
c     --------
 555  CONTINUE
c     --------
 550  CONTINUE
c     --------






 69   CONTINUE
c     --------
 68   CONTINUE
c     --------
 67   CONTINUE
c     --------
 66   CONTINUE
c     --------










c*************************
c     *
c     MAXIMUM product used:
c     -------------------
c*    NN = ITEM
c*    CALL SORTFAST(NN,CCPRD)
c*    CQEST = CCPRD(NN)
c*************************

      IF (ITEM.GT.0) GO TO 2000
      write (6,*) 'stop EST666: ITEM=',ITEM
      STOP
 2000 A=0


c
c     -----------------
c**** C6EST = C6EST/TEM
c     -----------------
c     ---------------------
      C6EST = DABS(C6EST)
c     --------------------- 
c     ---------------------
      C6EST = SQRT(C6EST)
c     --------------------- 
c     *
c*    write (6,*) 'subr. EST666  for 6-tuples'
c*    write (6,*) 'I1,I2,I3,I4,I5,I6;K1,K2,K3,K4,K5,K6'
c*    write (6,*) I1,I2,I3,I4,I5,I6,K1,K2,K3,K4,K5,K6
c*    write (6,*) 'C6EST=',C6EST 
c*    write (6,*) 'C6ADD=',C6ADD 
c*    write (6,*) 'total orb NTOT=',NTOT
c*    write (6,*) 'occ orb NOCC=',NOCC
c*    write (6,*) 'virt orb NVIRT=',NVIRT
c*    write (6,*) '# of distinct comb. ITEM=',ITEM
c*    do J=1,ISUM
c*    write (6,*) 'IDOLD1 2=',IDOLD1(J),IDOLD2(J)
c*    enddo
c*    write (6,*) 'subr. EST666: ICYC=',ICYC
c*    write (6,*) 'EST666: will stop'
c     *      
c
c     ----------
      RETURN
      END











c---------------------------------------------------------

      SUBROUTINE ESTTRIPLE(NTOT,NOCC,NVIRT,IEXCIT,NDOUBLE,
     *CISUM1,CISUMD,
     *NSCF,ISCF,IOCC,NVIR,IVIR,IOCCV,CTEST,CTADD)

c     ----------------------------------------------------

      implicit double precision(a-h,o-z)
C
      parameter (mxcomb = 100)
      parameter (mxorb = 100)
      parameter (mxdia = 1000)
      parameter (mxneb = 1000)
      parameter (mxsprod = 2400000)
      parameter (mxsprod1 = 100000)


      DIMENSION ISCF(mxorb),IOCC(mxorb),IVIR(mxorb) 
      DIMENSION IOCCV(mxorb)
c     *
      DIMENSION CISUMD(mxsprod1) 
      DIMENSION CISUM1(mxsprod1) 
c     *
      DIMENSION CCPRD(mxcomb) 

      DIMENSION ISCFA(mxorb),IVIRA(mxorb) 

      DIMENSION ICON1(mxorb),ICON2(mxorb)
      DIMENSION IDOLD1(500),IDOLD2(500)

c----------------------------------------------------------
c     This code gives an estimate for the
c     Q-quadruple space-product
c     as the sum of D-double sp products
c     that are consistent with this Q-sp.
c
c     Given:
c     ------
c     CISUM1(j), j-single excitation sp.
c     CISUMD(j), j=1,NDOUBLE numerical values of D-sp
c     coefficients that have been obtained before.
c
c     ISCF, IOCC & IVIR,IOCCV : occupied and virtual
c     orbitals with their occupation numbers.
c
c     On return:
c     ----------
c     CTEST = the numerical estimate for the Q-sp.
c----------------------------------------------------------


      IF (IEXCIT.EQ.3) GO TO 11
c     -------------------------
      write (6,*) 'ESTTRIPLE: excitations?'
      STOP
  11  A=0



      IADD = 0

      DO 25 N=1,NSCF
c     --------------
      IF (IOCC(N).LT.1) GO TO 30 
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
      GO TO 31
c     *
 30   A=0
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
      IADD = IADD + 1
      ISCFA(IADD) = ISCF(N)
 31   A=0
c     *
 25   CONTINUE


      IF (IADD.EQ.3) GO TO 35
c     -------------------------------
      write (6,*) 'inconsistent IOCC'
      write (6,*) 'ISCFA(JON)=',(ISCFA(JON),JON=1,3)
      write (6,*) '--------------------------------'
      write (6,*) 'ISCF(i)=',(ISCF(JON),JON=1,NSCF)
      write (6,*) 'IOCC(i)=',(IOCC(JON),JON=1,NSCF)
      write (6,*) '--------------------------------'
      write (6,*) 'subr. ESTTRIPLE'
      STOP
 35   A=0


      IADD = 0
c     *
      DO 50 M=1,NVIR
c     ---------------
      IF (IOCCV(M).GT.1) GO TO 55
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
      GO TO 61
 55   A=0
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
      IADD = IADD + 1
      IVIRA(IADD) = IVIR(M)
 61   A=0
c     *
 50   CONTINUE


      IF (IADD.EQ.3) GO TO 65
c     ----------------------------------
      write (6,*) 'inconsistent IOCCV'
      write (6,*) 'IVIRA(JON)=',(IVIRA(JON),JON=1,3)
      write (6,*) '--------------------------------'
      write (6,*) 'IVIR(i)=',(IVIR(JON),JON=1,NVIR)
      write (6,*) 'IOCCV(i)=',(IOCCV(JON),JON=1,NVIR)
      write (6,*) 'subr. ESTTRIPLE'
      STOP
 65   A=0


c     -----------------
c     ISCFA(j), j=1,3
c     IVIRA(k), k=1,3.
c     -----------------

      IEXCIT1 = 1 
      IEXCIT2 = 2

      CTEST = 0.0D0
      CTADD = 0.0D0

      ISUM = 0
      ITEM = 0
      TEM = 0.0D0

      DO JJJ=1,9
      CCPRD(JJJ) = 0.0D0
      ENDDO


c     *********
c     occupied:
c     *********
      I1 = ISCFA(1)
      I2 = ISCFA(2)
      I3 = ISCFA(3)
c     *******
c     virtual(1):
c     ------------
      K1 = IVIRA(1)
      K2 = IVIRA(2)
      K3 = IVIRA(3)
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------

c     --------------------------------
c     *******
c     virtual(2):
c     -----------
      K1 = IVIRA(2)
      K2 = IVIRA(1)
      K3 = IVIRA(3)
c     *
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     *
      CALL CHECKSAME2(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 101
c     *
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------
 101  A=0
c     --------------------------------
c     virtual(3):
c     -----------
      K1 = IVIRA(3)
      K2 = IVIRA(1)
      K3 = IVIRA(2)
c     *
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     ***
      CALL CHECKSAME2(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 102
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------
 102  A=0
c     --------------------------------




c     *********
c     occupied:
c     *********
      I1 = ISCFA(2)
      I2 = ISCFA(1)
      I3 = ISCFA(3)
c     **********
c     virtual(4):
c     -----------
      K1 = IVIRA(1)
      K2 = IVIRA(2)
      K3 = IVIRA(3)
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     ***
      CALL CHECKSAME2(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 105
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------
 105  A=0
c     --------------------------------
c     virtual(5):
c     -----------
      K1 = IVIRA(2)
      K2 = IVIRA(1)
      K3 = IVIRA(3)
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     ***
      CALL CHECKSAME2(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 106
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------
 106  A=0
c     --------------------------------
c     *******
c     virtual(6):
c     ------------
      K1 = IVIRA(3)
      K2 = IVIRA(1)
      K3 = IVIRA(2)
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     ***
      CALL CHECKSAME2(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 107
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------
 107  A=0
c     --------------------------------




c     *********
c     occupied:
c     -------
      I1 = ISCFA(3)
      I2 = ISCFA(1)
      I3 = ISCFA(2)
c     **********
c     virtual(7):
c     -------------
      K1 = IVIRA(1)
      K2 = IVIRA(2)
      K3 = IVIRA(3)
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     ***
      CALL CHECKSAME2(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 111 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------
 111  A=0
c     --------------------------------
c     ***
c     virtual(8):
c     -------------
      K1 = IVIRA(2)
      K2 = IVIRA(1)
      K3 = IVIRA(3)
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     ***
      CALL CHECKSAME2(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 112 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------
 112  A=0
c     --------------------------------
c     *
c     virtual(9):
c     -------------
      K1 = IVIRA(3)
      K2 = IVIRA(1)
      K3 = IVIRA(2)
c     *
      CALL NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON1) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT1,ICON1,IDNUM1) 
c     *
      CALL NUMDOUBLE(NTOT,NOCC,NVIRT,I2,I3,K2,K3,ICON2) 
      CALL LABELSDTQ(NTOT,NOCC,NVIRT,IEXCIT2,ICON2,IDNUM2) 
c     ***
      CALL CHECKSAME2(ISUM,IDOLD1,IDOLD2,IDNUM1,IDNUM2,
     *IPUT)
      ISUM = ISUM + 1
      IDOLD1(ISUM) = IDNUM1
      IDOLD2(ISUM) = IDNUM2
      IF (IPUT.NE.1) GO TO 113 
c     ***
      SGN = 1.0D0
      ITEM = ITEM + 1
      TEM = TEM + 1.0D0
      CCPRD(ITEM) = (CISUM1(IDNUM1)*CISUMD(IDNUM2))
      CTEST=CTEST+(SGN*CISUM1(IDNUM1)*CISUMD(IDNUM2))
c     -----------------------------------------------
 113  A=0




c     *
c     ---------------------
c     maximum-product used:
c     ---------------------
c*    NN = ITEM
c*    CALL SORTFAST(NN,CCPRD)
c     -------------------
c*    CTEST = CCPRD(NN)
c     *
      IF (ITEM.GT.0) GO TO 200
      write (6,*) 'ESTTRIP: TEM=',TEM
 200  A=0
c     ---------------------
      CTADD = CTEST
c     ---------------------
c     *
c     ---------------------
      CTEST = CTEST / TEM
c     ---------------------
      CTEST = DABS(CTEST)
c     ---------------------
      CTEST = SQRT(CTEST)
c     --------------------- 
c     *
c*    write (6,*) 'subr. ESTTRIPLE'
c*    write (6,*) 'I1,I2,I3,;K1,K2,K3'
c*    write (6,*) I1,I2,I3,K1,K2,K3
c*    write (6,*) 'CTEST=',CTEST 
c     ----------
      RETURN
      END










c     -------------------------
      SUBROUTINE SORTFAST(N,RA)
c     -------------------------
c     *
      implicit double precision (a-h,o-z)
c     *
      parameter (mxsprod=2400000)
      parameter (mxcomb=10000)
c     *
      DIMENSION RA(mxcomb)
c     *
c---------------------------
      L = N/2 + 1
c
      IR = N

 10   CONTINUE
      IF (L.GT.1) THEN
      L = L - 1
      RRA = RA(L)
      ELSE
c
      RRA = RA(IR)
      RA(IR) = RA(1)
      IR = IR - 1
      IF (IR.EQ.1) THEN
      RA(1) = RRA
      RETURN
c
      ENDIF
c
      ENDIF
c
      I = L
      J = L + L
 20   IF (J.LE.IR) THEN
c
      IF (J.LT.IR) THEN
c
      IF (RA(J).LT.RA(J+1)) J=J+1
      ENDIF
c
      IF (RRA.LT.RA(J)) THEN
      RA(I) = RA(J)
      I=J
      J=J+J
      ELSE
c
      J=IR+1
      ENDIF
c
      GO TO 20
      ENDIF
c
      RA(I) = RRA
      GO TO 10
c
c     -------------
      RETURN
      END




c     ---------------------------------
      SUBROUTINE SORT2(N,RA,RC,IRD)
c     ---------------------------------
c     *
      implicit double precision (a-h,o-z)
c     *
      parameter (mxsprod=2400000)
c     *
      DIMENSION RA(mxsprod)
      DIMENSION IRD(mxsprod)
c**   DIMENSION IRB(mxsprod)
      DIMENSION RC(mxsprod)
c     *
c---------------------------
      L = N/2 + 1
c
      IR = N

 10   CONTINUE
      IF (L.GT.1) THEN
      L = L - 1
      RRA = RA(L)
      RRC = RC(L)
      IRRD = IRD(L)
      ELSE
c
      RRA = RA(IR)
      RRC = RC(IR)
      IRRD = IRD(IR)
c
      RA(IR) = RA(1)
      RC(IR) = RC(1)
      IRD(IR) = IRD(1)
c
      IR = IR - 1

      IF (IR.EQ.1) THEN
      RA(1) = RRA
      RC(1) = RRC
      IRD(1) = IRRD
      RETURN
c
      ENDIF
c
      ENDIF
c
      I = L
      J = L + L
 20   IF (J.LE.IR) THEN
c
      IF (J.LT.IR) THEN
c
      IF (RA(J).LT.RA(J+1)) J=J+1
      ENDIF
c
      IF (RRA.LT.RA(J)) THEN
      RA(I) = RA(J)
      RC(I) = RC(J)
      IRD(I) = IRD(J)
      I=J
      J=J+J
      ELSE
c
      J=IR+1
      ENDIF
c
      GO TO 20
      ENDIF
c
      RA(I) = RRA
      RC(I) = RRC
      IRD(I) = IRRD
      GO TO 10
c
c     -------------
      RETURN
      END







c     ---------------------------------
      SUBROUTINE SORT02(N,RA,RC,IRD)
c     ---------------------------------
c     *
      implicit double precision (a-h,o-z)
c     *
      parameter (mxsprod=2400000)
      parameter (mxsprod1=100000)
c     *
      DIMENSION RA(mxsprod1)
      DIMENSION IRD(mxsprod1)
c**** DIMENSION IRB(mxsprod1)
      DIMENSION RC(mxsprod1)
c     *
c---------------------------
      L = N/2 + 1
c
      IR = N

 10   CONTINUE
      IF (L.GT.1) THEN
      L = L - 1
      RRA = RA(L)
      RRC = RC(L)
      IRRD = IRD(L)
      ELSE
c
      RRA = RA(IR)
      RRC = RC(IR)
      IRRD = IRD(IR)
c
      RA(IR) = RA(1)
      RC(IR) = RC(1)
      IRD(IR) = IRD(1)
c
      IR = IR - 1

      IF (IR.EQ.1) THEN
      RA(1) = RRA
      RC(1) = RRC
      IRD(1) = IRRD
      RETURN
c
      ENDIF
c
      ENDIF
c
      I = L
      J = L + L
 20   IF (J.LE.IR) THEN
c
      IF (J.LT.IR) THEN
c
      IF (RA(J).LT.RA(J+1)) J=J+1
      ENDIF
c
      IF (RRA.LT.RA(J)) THEN
      RA(I) = RA(J)
      RC(I) = RC(J)
      IRD(I) = IRD(J)
      I=J
      J=J+J
      ELSE
c
      J=IR+1
      ENDIF
c
      GO TO 20
      ENDIF
c
      RA(I) = RRA
      RC(I) = RRC
      IRD(I) = IRRD
      GO TO 10
c
c     -------------
      RETURN
      END







c     -----------------------------------
      SUBROUTINE SORT1(N,RA,IRB)
c     -----------------------------------
c     *
      implicit double precision (a-h,o-z)
c     *
c*    parameter (mxsprod=500000)
      parameter (mxsprod=2400000)
c     *
      DIMENSION RA(mxsprod)
      DIMENSION IRB(mxsprod)
c     *
      write (6,*) 'subr. SORT1 '
      STOP

c---------------------------
      L = N/2 + 1
c
      IR = N

 10   CONTINUE
      IF (L.GT.1) THEN
      L = L - 1
      RRA = RA(L)
      IRRB = IRB(L)
      ELSE
c
      RRA = RA(IR)
      IRRB = IRB(IR)
c
      RA(IR) = RA(1)
      IRB(IR) = IRB(1)
c
      IR = IR - 1

      IF (IR.EQ.1) THEN
      RA(1) = RRA
      IRB(1) = IRRB
      RETURN
c
      ENDIF
c
      ENDIF
c
      I = L
      J = L + L
 20   IF (J.LE.IR) THEN
c
      IF (J.LT.IR) THEN
c
      IF (RA(J).LT.RA(J+1)) J=J+1
      ENDIF
c
      IF (RRA.LT.RA(J)) THEN
      RA(I) = RA(J)
      IRB(I) = IRB(J)
      I=J
      J=J+J
      ELSE
c
      J=IR+1
      ENDIF
c
      GO TO 20
      ENDIF
c
      RA(I) = RRA
      IRB(I) = IRRB
      GO TO 10
c
c     -------------
      RETURN
      END











c     -----------------------------------------
      SUBROUTINE INFO(NTOT,NOCC,NVIRT,ICON,
     *NSCF,ISCF,IOCC,NVIR,IVIR,IOCCV,IEXCIT) 
c     -----------------------------------------
c     *
      implicit double precision (a-h,o-z)
      parameter (mxorb = 100)
c     *
      DIMENSION ICON(mxorb),ISCF(mxorb),IOCC(mxorb)
      DIMENSION IVIR(mxorb),IOCCV(mxorb)

      DIMENSION MOSCF(mxorb),MOVIR(mxorb)
c     *
C----------------------------------------
C      ICON(i), i=1,NORB
c
c      ICOUNT = ISCF  - number of SCF MOs
c      IVIR = ICOUNTS - number of VIR MOs
c----------------------------------------

       NORB = NTOT

c      Initialize:
c      ***********
       DO JA=1,NORB
c      ------------
       MOSCF(JA)=0
       ISCF(JA)=0
       IOCC(JA)=0
       MOVIR(JA)=0
       IVIR(JA)=0
       IOCCV(JA)=0
       ENDDO
c      ***


       ICOUNT = 0
       IEXC1 = 0




       DO II=1,NOCC
c ***  -----------

       NMO = 2 - ICON(II)
c      -----------------
       IF (NMO.EQ.2) go to 530 
       go to 532
  530  A=0
       IEXC1 = IEXC1 + NMO
       ICOUNT = ICOUNT + 1
       MOSCF(ICOUNT) = II
       ISCF(ICOUNT) = II
       IOCC(ICOUNT) = 0 
c      *
       go to 540 
c      -----------------------
  532  A=0
       IF (NMO.EQ.1) go to 534 
       go to 536
  534  A=0
       IEXC1 = IEXC1 + NMO
       ICOUNT = ICOUNT + 1
       MOSCF(ICOUNT) = II
       ISCF(ICOUNT) = II
       IOCC(ICOUNT) = 1
c      *
       go to 540 
c      ----------------------
  536  A=0
       IF (NMO.EQ.0) go to 540 
       write (6,*) 'undetermined occupation, stop'
       STOP

  540  A=0
c ***  ---
       ENDDO
c ***  -----------



       ICOUNTS = 0
       IEXC2 = 0


c      ^^^^^^^^^^^^^
c      VIRTUAL SPACE
c      ^^^^^^^^^^^^^
       DO III=NOCC+1,NORB
c ***  ------------------

       NMO = ICON(III)
c      -----------------
       IF (NMO.EQ.2) go to 550 
       go to 552
  550  A=0
       IEXC2 = IEXC2 + NMO 
       ICOUNTS = ICOUNTS + 1
       MOVIR(ICOUNTS) = III
       IVIR(ICOUNTS) = III
       IOCCV(ICOUNTS) = 2
c      *
       go to 560 
c      -----------------------
  552  A=0
       IF (NMO.EQ.1) go to 554 
       go to 556
  554  A=0
       IEXC2 = IEXC2 + NMO 
       ICOUNTS = ICOUNTS + 1
       MOVIR(ICOUNTS) = III
       IVIR(ICOUNTS) = III
       IOCCV(ICOUNTS) = 1
c      *
       go to 560 
c      ----------------------
  556  A=0
       IF (NMO.EQ.0) go to 560 
       write (6,*) 'undetermined occupation, stop'
       STOP

  560  A=0
c ***  ---
       ENDDO




       NSCF = ICOUNT
       NVIR = ICOUNTS




       IF (IEXC1.EQ.IEXC2) go to 565
       write (6,*) 'problem with excitations, stop'
       STOP

 565   A=0
c      ***
       IEXCIT = IEXC1
c
c
c       write (6,*) 'excitation IEXCIT =',IEXCIT
c       write (6,*) 'tot numb SCF-MOs =',NSCF
c       write (6,*) 'tot numb VIR-MOs =',NVIR
c      -------------
c       DO J=1,NSCF
c       write (6,*) 'ISCF/IOCC =',ISCF(J),IOCC(J)
c       ENDDO
c       write (6,*) '----------------------------'
c      -----------------
c       DO J=1,NVIR
c       write (6,*) 'IVIR/IOCCV =',IVIR(J),IOCCV(J)
c       ENDDO
c      ----------------------------
c     -------------
      RETURN
      END







c     --------------------------------------------------

      SUBROUTINE CHECKSAME(ISUM,IDOLD1,IDOLD2,ID1,ID2,IPUT)

c     --------------------------------------------------
c
c      RETURNS:
c      *****
c      IPUT (0 or 1) 
c      ----------------------------------

      implicit double precision(a-h,o-z)
c     -------------------------------------------
      parameter (mxorb = 100)
      parameter (mxCI = 1000000, mxstring=2000000, mxna=200)
c
      DIMENSION IDOLD1(500),IDOLD2(500)
c
c
      IPUT = 0
c     --------


c     ***
      DO 20 I=1,ISUM
c     --------------

      IF (IDOLD1(I).NE.ID1) GO TO 25 
      IF (IDOLD2(I).NE.ID2) GO TO 25
      GO TO 100 
c
 25   A=0
      IF (IDOLD1(I).NE.ID2) GO TO 30
      IF (IDOLD2(I).NE.ID1) GO TO 30 
      GO TO 100

 30   A=0
c     ***
c     --------
 20   CONTINUE


c     If got here, all must be distinct:
c     ----------------------------------
      IPUT = 1
c     -----------
 100  A=0
c     ------------
      RETURN
      END





c     --------------------------------------------------
c
      SUBROUTINE CHECKSAME2(ISUM,IDOLD1,IDOLD2,ID1,ID2,IPUT)
c
c     --------------------------------------------------
c
c      RETURNS:
c      *****
c      IPUT (0 or 1) 
c
c      IDOLD1(j) - single excitation ID-number
c      IDOLD2(j) - double excitation ID-number
c      ----------------------------------

      implicit double precision(a-h,o-z)
c     -------------------------------------------
      parameter (mxorb = 100)
      parameter (mxCI = 1000000, mxstring=2000000, mxna=200)
c
      DIMENSION IDOLD1(500),IDOLD2(500)
c
c
      IPUT = 0


c     ***
      DO 20 I=1,ISUM
c     --------------

      IF (IDOLD1(I).NE.ID1) GO TO 25 
      IF (IDOLD2(I).NE.ID2) GO TO 25
      GO TO 100 
c
 25   A=0
c     ***
c     --------
 20   CONTINUE


c     If got here, all must be distinct:
c     ----------------------------------
      IPUT = 1
c     -----------
 100  A=0
c     ------------
      RETURN
      END







c     -----------------------------------------------------
      SUBROUTINE NUMDOUBLE(NTOT,NOCC,NVIRT,I1,I2,K1,K2,ICON) 
c     -----------------------------------------------------
c     *
c      NUMDOUBLE 
c      ---------------
c      This code turns I1,I2/K1,K2 into ICON(j)
c      --------------------------------------- 
c      RETURNS:
c      *****
c      ICON(j), j=1,16
c      ----------------------------------
      implicit double precision(a-h,o-z)
c     -------------------------------------------
      parameter (mxorb = 100)
c
      DIMENSION ICON(mxorb)
c     -----------------------

c**** CO2:
c**** NOCC = 8
c**** NTOT = 16


      DO 25 IA=1,NOCC
c     ----------------
      ICON(IA) = 2
      IF (I1.NE.IA) GO TO 30 
      ICON(IA) = ICON(IA) - 1
c     *
 30   A=0
      IF (I2.NE.IA) GO TO 35 
      ICON(IA) = ICON(IA) - 1
c     *
 35   A=0
c     ----------------
 25   CONTINUE


      DO 50 IA=NOCC+1,NTOT
c     ---------------------
      ICON(IA) = 0 
      IF (K1.NE.IA) GO TO 60 
      ICON(IA) = ICON(IA) + 1
c     *
 60   A=0
      IF (K2.NE.IA) GO TO 65 
      ICON(IA) = ICON(IA) + 1
c     *
 65   A=0
c     ----------------
 50   CONTINUE
c
c*
c     write (6,*) 'subr NUMDOUBLE '
c     write (6,*) 'I1,I2,K1,K2',I1,I2,K1,K2
c     write (6,*) (ICON(J), J=1,NTOT)
c     write (6,*) '***'
c*
c     --------
      RETURN
      END





c     -----------------------------------------------------
      SUBROUTINE NUMTRIPLE(NTOT,NOCC,NVIRT,I1,I2,I3,
     *K1,K2,K3,ICON) 
c     -----------------------------------------------------
c     *
c     NUMTRIPLE: 
c     ----------
c     *
c     This code turns I1,I2,I3/K1,K2,K3 into ICON(j)
c     ---------------------------------------------- 
c     RETURNS:
c     *****
c     ICON(j), j=1,16
c     ----------------------------------
      implicit double precision(a-h,o-z)
c     ----------------------------------
      parameter (mxorb = 100)
c
      DIMENSION ICON(mxorb)
c     -----------------------

c**** CO2:
c**** NOCC = 8
c**** NTOT = 16


      DO 25 IA=1,NOCC
c     ----------------
      ICON(IA) = 2
      IF (I1.NE.IA) GO TO 30 
      ICON(IA) = ICON(IA) - 1
c     *
 30   A=0
      IF (I2.NE.IA) GO TO 35 
      ICON(IA) = ICON(IA) - 1
c     *
 35   A=0
c     *
      IF (I3.NE.IA) GO TO 45 
      ICON(IA) = ICON(IA) - 1
c     *
 45   A=0
c     ----------------
      ICONAS = ICON(IA)
c     ----------------
c     *
      IF (ICONAS.GE.0) GO TO 49
      write (6,*) 'subr. NUMTRIPLE: problem'
      write (6,*) ' ICON(ia) =',ICONAS
      STOP
c     *
 49   A=0
c     *
c     ----------------
 25   CONTINUE




      DO 50 IA=NOCC+1,NTOT
c     ---------------------
      ICON(IA) = 0 
      IF (K1.NE.IA) GO TO 60 
      ICON(IA) = ICON(IA) + 1
c     *
 60   A=0
      IF (K2.NE.IA) GO TO 65 
      ICON(IA) = ICON(IA) + 1
c     *
 65   A=0
c     *
      IF (K3.NE.IA) GO TO 75 
      ICON(IA) = ICON(IA) + 1
c     *
 75   A=0
c     -----------------
      ICONAS = ICON(IA)
c     -----------------
c     *
      IF (ICONAS.LT.3) GO TO 79
      write (6,*) 'subr. NUMTRIPLE: problem'
      write (6,*) ' ICON(ia) =',ICONAS
      STOP
c     *
 79   A=0
c     *
c     *
c     ----------------
 50   CONTINUE
c
c     *
c***  write (6,*) 'subr NUMTRIPLE '
c***  write (6,*) 'I1,I2,I3,K1,K2,K3',I1,I2,I3,K1,K2,K3
c***  write (6,*) (ICON(J), J=1,NTOT)
c**** write (6,*) 'NUMTRIPLE subr. stop ***'
c
c*
c     --------
      RETURN
      END









c     -----------------------------------------------------
      SUBROUTINE NUMQUADRO(NTOT,NOCC,NVIRT,I1,I2,I3,I4,
     *K1,K2,K3,K4,ICON) 
c     -----------------------------------------------------
c     *
c     NUMQUADRO: 
c     ----------
c     *
c     This code turns I1,I2,I3,I4/K1,K2,K3,K4 into ICON(j)
c     ---------------------------------------------------- 
c     RETURNS:
c     *****
c     ICON(j), j=1,16
c
c
c     ----------------------------------
      implicit double precision(a-h,o-z)
c     ----------------------------------
      parameter (mxorb = 100)
c
      DIMENSION ICON(mxorb)
c     -----------------------

c**** CO2:
c**** NOCC = 8
c**** NTOT = 16


      DO 25 IA=1,NOCC
c     ----------------
      ICON(IA) = 2
      IF (I1.NE.IA) GO TO 30 
      ICON(IA) = ICON(IA) - 1
c     *
 30   A=0
      IF (I2.NE.IA) GO TO 35 
      ICON(IA) = ICON(IA) - 1
c     *
 35   A=0
c     *
      IF (I3.NE.IA) GO TO 45 
      ICON(IA) = ICON(IA) - 1
c     *
 45   A=0
c     *
      IF (I4.NE.IA) GO TO 55 
      ICON(IA) = ICON(IA) - 1
c     *
 55   A=0
c     *
c     ----------------
      ICONAS = ICON(IA)
c     ----------------
c     *
      IF (ICONAS.GE.0) GO TO 49
      write (6,*) 'subr. NUMQUADRO: problem'
      write (6,*) ' ICON(ia) =',ICONAS
      STOP
c     *
 49   A=0
c     *
c     ----------------
 25   CONTINUE




      DO 50 IA=NOCC+1,NTOT
c     ---------------------
      ICON(IA) = 0 
      IF (K1.NE.IA) GO TO 60 
      ICON(IA) = ICON(IA) + 1
c     *
 60   A=0
      IF (K2.NE.IA) GO TO 65 
      ICON(IA) = ICON(IA) + 1
c     *
 65   A=0
c     *
      IF (K3.NE.IA) GO TO 75 
      ICON(IA) = ICON(IA) + 1
c     *
 75   A=0
c     *
      IF (K4.NE.IA) GO TO 85 
      ICON(IA) = ICON(IA) + 1
c     *
 85   A=0
c     *
c     -----------------
      ICONAS = ICON(IA)
c     -----------------
c     *
      IF (ICONAS.LT.3) GO TO 79
      write (6,*) 'subr. NUMQUADRO: problem'
      write (6,*) ' ICON(ia) =',ICONAS
      STOP
c     *
 79   A=0
c     *
c     *
c     ----------------
 50   CONTINUE
c
c     *
c*    write (6,*) 'subr NUMQUADRO '
c*    write (6,*) 'I1,I2,I3,I4,K1,K2,K3,K4'
c*    write (6,*) I1,I2,I3,I4,K1,K2,K3,K4
c*    write (6,*) (ICON(J), J=1,NTOT)
c*    write (6,*) '-------------------------' 
c*    write (6,*) 'NUMquadro subr. stop ***'
c*
c     --------
      RETURN
      END






c     -----------------------------------------------------
      SUBROUTINE NUM555(NTOT,NOCC,NVIRT,ICON2,ICON3,ICON,
     *IPUT) 
c     -----------------------------------------------------
c     *
c     NUMQUADRO: 
c     ----------
c     *
c     This code generates 5-tuple: ICON(j)= ICON2 + ICON3
c     ---------------------------------------------------- 
c     RETURNS:
c     *****
c     ICON(j), j=1,NTOT
c     -----------------
c**** CO2:
c**** NOCC = 8
c**** NTOT = 16

c
c
c     ----------------------------------
      implicit double precision(a-h,o-z)
c     ----------------------------------
c     *
      parameter (mxorb = 100)
c     *
      DIMENSION ICON2(mxorb)
      DIMENSION ICON3(mxorb)
      DIMENSION ICON(mxorb)
c     -----------------------


      DO IA=1,NTOT
c     ----------------
      ICON(IA) = 0
      ENDDO
c     ----------------



      IPUT = 0



c     ************
c     OCCUPIED MO:
c     ************
      DO IA=1,NOCC
c     ----------------
      ICON(IA) = ICON2(IA) + ICON3(IA)
      ENDDO
c     ----------------


      DO IA=1,NOCC
c     ----------------
      ICON(IA) = 4 - ICON(IA)
      ENDDO
c     ----------------


      DO IA=1,NOCC
c     ----------------
      ICON(IA) = 2 - ICON(IA)

      IF (ICON(IA).LT.0) GO TO 100 
      IF (ICON(IA).GT.2) GO TO 100 

      ENDDO
c     ----------------




c     ************
c     VIRTUAL MO:
c     ************
      DO IA=NOCC+1,NTOT
c     ----------------
      ICON(IA) = ICON2(IA) + ICON3(IA)
      ENDDO
c     ----------------



      DO IA=NOCC+1,NTOT
c     ----------------
      IF (ICON(IA).LT.0) GO TO 100 
      IF (ICON(IA).GT.2) GO TO 100 
      ENDDO
c     ----------------

c     --------
      IPUT = 1
c     --------

c     *
c     ***
 100  A=0
c     ***

c     *
c*    write (6,*) 'subr NUM555: '
c*    write (6,*) '************ '
c*    write (6,*) ' IPUT = ',IPUT
c*    write (6,*) (ICON2(J), J=1,NTOT)
c*    write (6,*) (ICON3(J), J=1,NTOT)
c*    write (6,*) (ICON(J), J=1,NTOT)
c*    write (6,*) '-------------------------' 
c*    write (6,*) 'NUM555 subr. stop ***'
c     *
c     *
c     --------
      RETURN
      END





c     -------------------------------------
      SUBROUTINE MULTIPLY6(N1,N2,N3,FACTOR)
c     -------------------------------------
c     **********
c     Deduce factor (binomial) to multiply possibilities: 
c     ----------
c     6-tuple excitations only: 
c     ---------------------------------------------------- 
c     ----------------------------------
      implicit double precision(a-h,o-z)
c     ----------------------------------
c*    parameter (mxorb = 100)
c     ----------------------------------


      IF (N1.EQ.N2) GO TO 50
c     ----------------------
c     N1.ne.N2
      IF (N1.EQ.N3) GO TO 51
c     ***
c     N1.ne.N2
c     N1.ne.N3
      IF (N2.EQ.N3) GO TO 55
c     N2.ne.N3
c     --------------
      FACTOR = 6.0D0
c     --------------
      GO TO 100


 51   A=0
c     ***
c     N1.ne.N2
c     N1=N3
c     --------------
      FACTOR = 3.0D0
c     --------------
      GO TO 100


 55   A=0
c     ***
c     N1.ne.N2
c     N1.ne.N3
c     N2=N3
c     --------------
      FACTOR = 3.0D0
c     --------------
      GO TO 100


 50   A=0
c     ***
      IF (N2.EQ.N3) GO TO 60
c     N1.eq.N2
c     N2.ne.N3
c     --------------
      FACTOR = 3.0D0
c     --------------
      GO TO 100


 60   A=0
c     ***
c     N1=N2=N3
c     --------------
      FACTOR = 1.0D0
c     --------------

c     ***
 100  A=0
c     ***
c     *
c     --------
      RETURN
      END















c     -------------------------------------------------
      SUBROUTINE NUM666(NTOT,NOCC,NVIRT,ICON2,ICON3,
     *ICON4,ICON,IPUT)
c     -----------------------------------------------------
c     NUM666: 
c     ----------
c     This code generates 6-tuple: 
c     ICON(j)= ICON2 + ICON3 + ICON4
c     6-       2-      2-      2-
c     ---------------------------------------------------- 
c     RETURNS:
c     *****
c     ICON(j), j=1,NTOT
c     -----------------
c**** CO2:
c**** NOCC = 8
c**** NTOT = 16
c     ----------------------------------
      implicit double precision(a-h,o-z)
c     ----------------------------------
c     *
      parameter (mxorb = 100)
c     *
      DIMENSION ICON2(mxorb)
      DIMENSION ICON3(mxorb)
      DIMENSION ICON4(mxorb)
      DIMENSION ICON(mxorb)
c     -----------------------

      DO IA=1,NTOT
c     ----------------
      ICON(IA) = 0
      ENDDO
c     ----------------


      IPUT = 0


c     ************
c     OCCUPIED MO:
c     ************
      DO IA=1,NOCC
c     ----------------
      ICON(IA) = ICON2(IA) + ICON3(IA) + ICON4(IA)
      ENDDO
c     ----------------


      DO IA=1,NOCC
c     ----------------
      ICON(IA) = 6 - ICON(IA)
      ENDDO
c     ----------------

      DO IA=1,NOCC
c     ----------------
      ICON(IA) = 2 - ICON(IA)

      IF (ICON(IA).LT.0) GO TO 100 
      IF (ICON(IA).GT.2) GO TO 100 

      ENDDO
c     ----------------



c     ************
c     VIRTUAL MO:
c     ************
      DO IA=NOCC+1,NTOT
c     ----------------
      ICON(IA) = ICON2(IA) + ICON3(IA) + ICON4(IA)
      ENDDO
c     ----------------


      DO IA=NOCC+1,NTOT
c     ----------------
      IF (ICON(IA).LT.0) GO TO 100 
      IF (ICON(IA).GT.2) GO TO 100 
      ENDDO
c     ----------------

c     --------
      IPUT = 1
c     --------

c     *
c     ***
 100  A=0
c     ***
c     *
c*    write (6,*) 'subr NUM666: '
c*    write (6,*) '************ '
c*    write (6,*) ' IPUT = ',IPUT
c*    write (6,*) (ICON2(J), J=1,NTOT)
c*    write (6,*) (ICON3(J), J=1,NTOT)
c*    write (6,*) (ICON4(J), J=1,NTOT)
c*    write (6,*) (ICON(J), J=1,NTOT)
c*    write (6,*) '-------------------------' 
c*    write (6,*) 'NUM666 subr. stop ***'
c     *
c     *
c     --------
      RETURN
      END









c     --------------------------------------------------
      SUBROUTINE NUMSINGLE(NTOT,NOCC,NVIRT,I1,K1,ICON) 
c     --------------------------------------------------
c     *
c      NUMSINGLE 
c      ---------------
c      This code turns I1/K1 into ICON(j)
c      --------------------------------------- 
c      RETURNS:
c      *****
c      ICON(j), j=1,16
c      ----------------------------------
      implicit double precision(a-h,o-z)
c     -------------------------------------------
      parameter (mxorb = 100)
c
      DIMENSION ICON(mxorb)
c     -----------------------

c**** CO2:
c**** NOCC = 8
c**** NTOT = 16


      DO 25 IA=1,NOCC
c     ----------------
      ICON(IA) = 2
      IF (I1.NE.IA) GO TO 30 
      ICON(IA) = ICON(IA) - 1
c     *
 30   A=0
c     ----------------
 25   CONTINUE


      DO 50 IA=NOCC+1,NTOT
c     ---------------------
      ICON(IA) = 0 
      IF (K1.NE.IA) GO TO 60 
      ICON(IA) = ICON(IA) + 1
c     *
 60   A=0
c     ----------------
 50   CONTINUE
c
c     --------
      RETURN
      END






c     --------------------------------------------------
      SUBROUTINE SYSTEMDOUBLE(I1,I2,K1,K2,IDNUM) 
c     --------------------------------------------------
c        NUMDOUBLE 
c      ---------------
c      This code labels doubly excited -space products in a 
c      systematic way and matches i1,i2/k1,k2
c      to this label. 
c
c      RETURNS:
c      *****
c      IDNUM
c      ----------------------------------

      implicit double precision(a-h,o-z)
c     -------------------------------------------
      parameter (mxorb = 100)
      parameter (mxCI = 1000000, mxstring=1000000, mxna=200)
c
c
c     First, make sure that I1,I2 / K1,K2
c     are ordered in increasing order!
c     ************************************

      IT1 = I1 
      IT2 = I2
      IF (I2.LT.I1) GO TO 15
      GO TO 17
  15  A=0
      IT1 = I2 
      IT2 = I1
  17  A=0
c
c
      KT1 = K1
      KT2 = K2
      IF (K2.LT.K1) GO TO 20 
      GO TO 21
  20  A=0
      KT1 = K2 
      KT2 = K1
  21  A=0
c
c     ***
      I1 = IT1
      I2 = IT2
c
      K1 = KT1
      K2 = KT2
c     ***


      NOCC = 8
      NTOT = 16
      ICOUNT = 0


c     ***
      DO 25 II1=1,NOCC
c     ----------------
      DO 30 II2=II1,NOCC
c     ------------------

      DO 50 KK1=NOCC+1,NTOT
c     ---------------------
      DO 55 KK2=KK1,NTOT
c     ---------------------

      ICOUNT = ICOUNT + 1


      IF (II1.NE.I1) GO TO 60 
      IF (II2.NE.I2) GO TO 60 
c
      IF (KK1.NE.K1) GO TO 60 
      IF (KK2.NE.K2) GO TO 60 
c
c     match has occured, write it out
c     -------------------------------
      IDNUM = ICOUNT
c     --------------
      GO TO 100
c     *
  60  A=0
c     *
  55  CONTINUE
  50  CONTINUE

  30  CONTINUE
  25  CONTINUE

 100  A=0
c     --------
      RETURN
      END










c     --------------------------------------------------
c
      SUBROUTINE WRITESD(IDETSD,CISUMSD)
c
c     --------------------------------------------------
c     This subroutine
c     ---------------

c    ************************************
c         WRITESD
c      ---------------
c      This code labels SD-space products in a 
c      systematic way and writes them out
c      ----------------------------------

      implicit double precision(a-h,o-z)
c     -------------------------------------------

      parameter (mxorb = 100)
      parameter (mxsprod1 = 100000)
      parameter (mxCI = 1000000, mxstring=2000000, mxna=200)
c
      DIMENSION IDETSD(mxsprod1),CISUMSD(mxsprod1) 
      DIMENSION NEL(mxorb) 
      DIMENSION NSCF1(mxorb),NVIR2(mxorb) 




C--------------------------------------------
C*     doubly and singly occupied MOs
C      ------------------------------

       ICOUNT = 0

       NOCC = 8
       NVIRT = 8
       NTOT = 16
C      -------------------------------




C***** code to generate space-products
C      -------------------------------
c



C      one-determinant (SCF) :
C      -----------------------
       DO 30 I=1,NOCC
       NEL(I) = 2
 30    CONTINUE
       DO 32 I=NOCC+1,NTOT
       NEL(I) = 0
 32    CONTINUE
C*
c***** WRITE (6,*) (NEL(I), I=1,NTOT) 
c***** ICOUNT = ICOUNT + 1
C*     -------------------





c      ***
c***** write (6,*) 'single excitations '
C      -------------------

       IADD1 = 0
       IADD2 = 0
       DO JAMA=1,NTOT
       NSCF1(JAMA) = 0
       NVIR2(JAMA) = 0
       ENDDO
c      *** 


       DO 35 I=1,NOCC
c      --------------
       DO 36 J=1,NOCC
       NEL(J)= 2
  36   CONTINUE

       NSING = NOCC+1-I
       NEL(NSING) = 1



c      ^^^^^^^^^^^^^^^^^^^^
c      virtual
c      ^^^^^^^^^^^^^^^^^^^^
       DO 39 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 39    CONTINUE



       DO 40 MM=NOCC+1,NTOT
c      --------------------
       NEL(MM) = 1 
c
c***** WRITE (6,*) (NEL(IX), IX=1,NTOT) 
c      --------------------
       ICOUNT = ICOUNT + 1
C*    --------------------
c      ***
       IADD1 = 0
       IADD2 = 0
       DO JAMA=1,NTOT
       NSCF1(JAMA) = 0
       NVIR2(JAMA) = 0
       ENDDO
c      *** 
      DO JO=1,NOCC
      IDEL = 2 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 42
      IADD1 = IADD1 + 1
      NSCF1(IADD1) = JO
 42   A=0  
      ENDDO
c     ***
      DO JO=NOCC+1,NTOT
      IDEL = 0 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 43
      IADD2 = IADD2 + 1
      NVIR2(IADD2) = JO
 43   A=0  
      ENDDO
c     ***
c     ***
      IF (CISUMSD(ICOUNT).LE.0.0D0) GO TO 41
c
      WRITE (6,*) ICOUNT,IDETSD(ICOUNT),
     *CISUMSD(ICOUNT),(NEL(IX),IX=1,NTOT),
     *NSCF1(1),NVIR2(1) 
c     ********************

 41   A=0

       NEL(MM) = 0
 40    CONTINUE
c      -------------
 35    CONTINUE






c***** write (6,*) 'DIAGONAL-double excitations'
c      ***

       IADD1 = 0
       IADD2 = 0
       DO JAMA=1,NTOT
       NSCF1(JAMA) = 0
       NVIR2(JAMA) = 0
       ENDDO



       DO 45 I=1,NOCC
c      -------------------

       DO 47 J=1,NOCC
       NEL(J)= 2
  47   CONTINUE

       NSING = NOCC+1-I
       NEL(NSING) = 0 



       IADD1 = 0
       IADD2 = 0

c      ^^^^^^^^^^^^^^^^^^^^
c      virtual
c      ^^^^^^^^^^^^^^^^^^^^
       DO 49 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 49    CONTINUE

       DO 50 MM=NOCC+1,NTOT
c      --------------------
       NEL(MM) = 2 
c
c***** WRITE (6,*) (NEL(L), L=1,NTOT) 
c      -------------------
       ICOUNT = ICOUNT + 1
C*    --------------------
c      ***
       IADD1 = 0
       IADD2 = 0
       DO JAMA=1,NTOT
       NSCF1(JAMA) = 0
       NVIR2(JAMA) = 0
       ENDDO
c      *** 
      DO JO=1,NOCC
      IDEL = 2 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 52 
      IADD1 = IADD1 + 1
      NSCF1(IADD1) = JO
 52   A=0  
      ENDDO
c     ***
      DO JO=NOCC+1,NTOT
      IDEL = 0 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 53 
      IADD2 = IADD2 + 1
      NVIR2(IADD2) = JO
 53   A=0  
      ENDDO
c     ***
c     ***
      IF (CISUMSD(ICOUNT).LE.0.0D0) GO TO 51 
c
      WRITE (6,*) ICOUNT,IDETSD(ICOUNT),
     *CISUMSD(ICOUNT),(NEL(IX),IX=1,NTOT),
     *NSCF1(1),NVIR2(1) 
c     *************************
c
 51   A=0

       NEL(MM) = 0
 50    CONTINUE






c      ^^^^^^^^^^^^^^^^^
c      virtual
c      ^^^^^^^^^^^^^^^^^
       DO 55 JJJ=1,NVIRT
c      ------------------
       DO 56 KKK=NOCC+1,NTOT
       NEL(KKK) = 0
  56   CONTINUE
       M = NTOT+1-JJJ
       NEL(M) = 1



       DO 60 II=1,M-NOCC-1
c      -------------------
       MMM = M - II
c      ------------
       NEL(MMM) = 1

c
c****  WRITE (6,*) (NEL(L), L=1,NTOT) 
c      -------------------
       ICOUNT = ICOUNT + 1
C*     -------------------
c      ***
       IADD1 = 0
       IADD2 = 0
       DO JAMA=1,NTOT
       NSCF1(JAMA) = 0
       NVIR2(JAMA) = 0
       ENDDO
c      *** 
      DO JO=1,NOCC
      IDEL = 2 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 57 
      IADD1 = IADD1 + 1
      NSCF1(IADD1) = JO
 57   A=0  
      ENDDO
c     ***
      DO JO=NOCC+1,NTOT
      IDEL = 0 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 58 
      IADD2 = IADD2 + 1
      NVIR2(IADD2) = JO
 58   A=0  
      ENDDO
c     ***
c     ***
      IF (CISUMSD(ICOUNT).LE.0.0D0) GO TO 61 
c
      WRITE (6,*) ICOUNT,IDETSD(ICOUNT),
     *CISUMSD(ICOUNT),(NEL(IX),IX=1,NTOT),
     *NSCF1(1),NVIR2(1),NVIR2(2) 
c     **************************

 61   A=0


       NEL(M-II) = 0
 60    CONTINUE

       NEL(M) = 0
 55    CONTINUE




c      --------------
 45    CONTINUE






c***** WRITE (6,*) 'S[1] and D[0] ICOUNT=',ICOUNT
c      ***
c      -------------------------
c***** write (6,*) 'mixed double excitations '
C      -------------------


       DO 275 I=1,NOCC
C      ---------------------

       DO 277 J=1,NOCC
       NEL(J)= 2
 277   CONTINUE

       M1 = NOCC+1-I 
c      -------------
       NEL(M1) = 1 


       DO 278 K1=1,M1-1
C      ----------------
       M2 = M1 - K1
c      *
       NEL(M2)= 1






c      ^^^^^^^^^^^^^^^^^^^^^
c      virtual
c      ^^^^^^^^^^^^^^^^^^^^^
       DO 280 JJ=NOCC+1,NTOT
       NEL(JJ) = 0
 280   CONTINUE

       DO 285 MM=NOCC+1,NTOT
c      ---------------------
       NEL(MM) = 2 
c
c***** WRITE (6,*) (NEL(L), L=1,NTOT) 
C*     --------------------
       ICOUNT = ICOUNT + 1
C*     --------------------
c      ***
       IADD1 = 0
       IADD2 = 0
       DO JAMA=1,NTOT
       NSCF1(JAMA) = 0
       NVIR2(JAMA) = 0
       ENDDO
c      *** 
      DO JO=1,NOCC
      IDEL = 2 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 282 
      IADD1 = IADD1 + 1
      NSCF1(IADD1) = JO
 282  A=0  
      ENDDO
c     ***
      DO JO=NOCC+1,NTOT
      IDEL = 0 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 283 
      IADD2 = IADD2 + 1
      NVIR2(IADD2) = JO
 283  A=0  
      ENDDO
c     ***
c     ***
c
c***** IMPORTANT!!!**************************
      IF (CISUMSD(ICOUNT).LE.0.0D0) GO TO 281 
c
      WRITE (6,*) ICOUNT,IDETSD(ICOUNT),
     *CISUMSD(ICOUNT),(NEL(IX),IX=1,NTOT),
     *NSCF1(1),NSCF1(2),NVIR2(1) 
c     ***************************************

 281  A=0



       NEL(MM) = 0
c      --------
 285   CONTINUE






c      ^^^^^^^^^
c      virtual
c      ^^^^^^^^^
       DO 300 JJJ=1,NVIRT
C      ------------------

       DO 301 KKK=NOCC+1,NTOT
       NEL(KKK) = 0
 301   CONTINUE
       M = NTOT+1-JJJ
c      ------------------------
       NEL(M) = 1


       DO 305 II=1,M-NOCC-1
c      --------------------
       MM = M - II
c      -----------
       NEL(MM) = 1

c
c***** WRITE (6,*) (NEL(L), L=1,NTOT) 
c      -------------------
       ICOUNT = ICOUNT + 1
C*     -------------------
c      ***
       IADD1 = 0
       IADD2 = 0
       DO JAMA=1,NTOT
       NSCF1(JAMA) = 0
       NVIR2(JAMA) = 0
       ENDDO
c      *** 
      DO JO=1,NOCC
      IDEL = 2 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 306 
      IADD1 = IADD1 + 1
      NSCF1(IADD1) = JO
 306  A=0  
      ENDDO
c     ***
      DO JO=NOCC+1,NTOT
      IDEL = 0 - NEL(JO)
      IF (IDEL.EQ.0) GO TO 307 
      IADD2 = IADD2 + 1
      NVIR2(IADD2) = JO
 307  A=0  
      ENDDO
c     ***
c     ***
      IF (CISUMSD(ICOUNT).LE.0.0D0) GO TO 303 
c
      WRITE (6,*) ICOUNT,IDETSD(ICOUNT),
     *CISUMSD(ICOUNT),(NEL(IX),IX=1,NTOT),
     *NSCF1(1),NSCF1(2),NVIR2(1),NVIR2(2) 
c     ***********************************

 303  A=0



       NEL(MM) = 0
 305   CONTINUE

       NEL(M) = 0
 300   CONTINUE


c      ***************
c      SCF-MOs
c      ---------------
       NEL(M1-K1) = 2
 278   CONTINUE


       NEL(M1) = 2
 275   CONTINUE


c      *****************
c      subr WRITESD ends
c      -----------------
       RETURN
       END





