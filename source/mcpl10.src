C 19 SEP 05 - MK  - BASED ON MODPOT WRITTEN BY DIETMAR KRAUSE,
C                   WITH CONTRIBUTIONS FROM S. A. DECKER
C
C*MODULE MCPLIB  *DECK MCPLIB
C
C======================================================================
C THE MODEL CORE POTENTIALS AND ASSOCIATED BASIS SETS WERE DEVELOPED BY
C Y. SAKAI,  E. MIYOSHI,  M. KLOBUKOWSKI, AND S. HUZINAGA,
C     "MODEL POTENTIALS FOR MOLECULAR CALCULATIONS.
C     I. THE SD-MP SET FOR TRANSITION METAL ATOMS SC THROUGH HG",
C     J. COMPUT. CHEM. 8 (1987) 226-255.
C Y. SAKAI, E.  MIYOSHI, M. KLOBUKOWSKI, AND S. HUZINAGA,
C     "MODEL POTENTIALS FOR MOLECULAR CALCULATIONS.
C     II. THE SPD-MP SET FOR TRANSITION METAL ATOMS SC THROUGH HG",
C      J. COMPUT. CHEM. 8 (1987) 256-264.
C Y. SAKAI, E. MIYOSHI, M, KLOBUKOWSKI, AND S. HUZINAGA,
C     "MODEL POTENTIALS FOR MAIN GROUP ELEMENTS", 
C     J. CHEM. PHYS. 106 (1997) 8084-8092.
C E. MIYOSHI, Y. SAKAI, K. TANAKA, AND M. MASAMURA
C     "RELATIVISTIC DSP-MODEL CORE POTENTIALS FOR MAIN GROUP ELEMENTS 
C     IN THE FOURTH, FIFTH, AND SIXTH-ROW AND THEIR APPLICATIONS"
C     J. MOL. STRUCT. (THEOCHEM), 451 (1998) 73-79.
C Y. SAKAI, E. MIYOSHI, AND H. TATEWAKI
C     "MODEL CORE POTENTIALS FOR THE LANTHANIDES"
C     J. MOL. STRUCT. (THEOCHEM), 451 (1998) 143-150.
C======================================================================
C
C     MCPLIB RETURNS VALENCE BASIS-SET INFORMATION
C     OR MODEL-POTENTIALS WITH CORE BASIS-SET INFORMATION
C
      SUBROUTINE MCPLIB(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
C
C     --- MMPMOD: 0  READ VALENCE BASIS-SET
C     ---         1  READ MODEL-POTENTIAL-PARAMETERS
C     ---            AND CORE-SHELL BASIS-SET
C
C     --- NUCZ:    NUCLEAR CHARGE
C     --- MMPTYP:  1  NON-RELATIVISTIC   MCP
C     ---          2  QUASI-RELATIVISTIC MCP
C     ---          3  READ IN MCP FROM $MCP GROUP
C     ---          4  TYPE NOT SPECIFIED
C     ---    W*10**I  FOR I-TH SHELL (I>=1)
C     ---    W: 1,2,3,4,6  S,P,D,F,L - SHELL
C     --- CDFCTR(I) CODED DEFAULT CONTRACTION OF VALENCE-SHELL I
C     --- NOSH     NO.OF (VALENCE OR CORE) SHELLS
C     --- KTYP(I)  TYPE OF SHELL I
C     --- KNG(I)   NO.OF GTFS IN SHELL I
C     --- EX(J)    EXPONENT OF J-TH GTF
C     --- C(J)     CONTRACTION COEFF.OF J-TH GTF
C     --- ZCORE    CORE-CHARGE
C     --- NOAN0    NO.OF A-TERMS FOR N=0
C     --- NOAN1    NO.OF A-TERMS FOR N=1
C     --- AN0(I)   A-PARAMETERS FOR N=0
C     --- ALPN0(I) ALPHA-EXPONENTS FOR N=0
C     --- AN1(I)   A-PARAMETERS FOR N=1
C     --- ALPN1(I) ALPHA-EXPONENTS FOR N=1
C     --- BPAR(I)  B-PARAMETERS
C     --- IERR(I)  ERROR-STATUS   FOR MMPMOD=0:
C     ---          0  OK; 1  NOT MMP-ATOM; 2 MMPTYP WRONG
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AN0(*),ALPN0(*),AN1(*),ALPN1(*)
      DIMENSION BPAR(*),KNG(*),KTYP(*)
      DIMENSION C(*),EX(*)
      CHARACTER*20 CDFCTR(*),SHINFO
C
      DIMENSION ITYP(10)
C
      IZSAVE=NUCZ
      IERR=0
C
C     WRITE(*,*) 'MMPMOD',MMPMOD
C     WRITE(*,*) 'MMPTYP',MMPTYP
C     WRITE(*,*) 'IZSAVE',IZSAVE
C ... MAIN GROUP 03 <= Z <= 10  LI - NE
      IF(IZSAVE.GE.03 .AND. IZSAVE.LE.10) THEN
        CALL MG0310(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,
     3                  BPAR,ITYP,IERR)
      ENDIF
C ... MAIN GROUP 11 <= Z <= 18  NA - AR
      IF(IZSAVE.GE.11 .AND. IZSAVE.LE.18) THEN
        CALL MG1118(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
      ENDIF
C ... MAIN GROUP 19 <= Z <= 36   K - KR
      IF   ( IZSAVE.EQ.19  .OR. IZSAVE.EQ.20
     + .OR. (IZSAVE.GE.31 .AND. IZSAVE.LE.36 ) ) THEN
        CALL MG1936(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
      ENDIF
C ... MAIN GROUP 37 <= Z <= 54  RB - XE
      IF   ( IZSAVE.EQ.37  .OR. IZSAVE.EQ.38
     + .OR. (IZSAVE.GE.49 .AND. IZSAVE.LE.54 ) ) THEN
        CALL MG3754(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
      ENDIF
C ... MAIN GROUP 55 <= Z <= 86  CS - RN
      IF   ( IZSAVE.EQ.55  .OR. IZSAVE.EQ.56
     + .OR. (IZSAVE.GE.81 .AND. IZSAVE.LE.86 ) ) THEN
        CALL MG5586(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
      ENDIF
C ... 1ST-ROW TRANSITION METALS 21 <= 30  V - ZN
      IF(IZSAVE.GE.21 .AND. IZSAVE.LE.30) THEN
        CALL TM2130(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
      ENDIF
C ... 2ND-ROW TRANSITION METALS 39 <= 48 Y - CD
      IF(IZSAVE.GE.39 .AND. IZSAVE.LE.48) THEN
        CALL TM3948(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
      ENDIF
C ... 3RD ROW TRANSITION METALS 71 <= 80 LU - HG
      IF(IZSAVE.GE.71 .AND. IZSAVE.LE.80) THEN
        CALL TM7180(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
      ENDIF
C ... LANTHANIDES 57 <= Z <=70  LA - YB
      IF(IZSAVE.GE.57 .AND. IZSAVE.LE.70) THEN
        CALL LA5770(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
      ENDIF
C
      RETURN
      END

C*MODULE MCPLIB  *DECK MG0310
      SUBROUTINE MG0310(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,BPAR,ITYP,IERR)
C
C     --- MMPMOD: 0  READ VALENCE BASIS-SET
C     ---         1  READ MODEL-POTENTIAL-PARAMETERS
C     ---            AND CORE-SHELL BASIS-SET
C
C     --- NUCZ:    NUCLEAR CHARGE
C     --- MMPTYP:  1  NON-RELATIVISTIC   MOD.POT.
C     ---          2  QUASI-RELATIVISTIC MOD.POT.
C     ---          3  READ IN MOD.POT.FROM $MMP GROUP
C     ---          4  TYPE NOT SPECIFIED
C     ---    W*10**I  FOR I-TH SHELL (I>=1)
C     ---    W: 1,2,3,4,6  S,P,D,F,L - SHELL
C     --- CDFCTR(I) CODED DEFAULT CONTRACTION OF VALENCE-SHELL I
C     --- NOSH     NO.OF (VALENCE OR CORE) SHELLS
C     --- KTYP(I)  TYPE OF SHELL I
C     --- KNG(I)   NO.OF GTFS IN SHELL I
C     --- EX(J)    EXPONENT OF J-TH GTF
C     --- C(J)     CONTRACTION COEFF.OF J-TH GTF
C     --- ZCORE    CORE-CHARGE
C     --- NOAN0    NO.OF A-TERMS FOR N=0
C     --- NOAN1    NO.OF A-TERMS FOR N=1
C     --- AN0(I)   A-PARAMETERS FOR N=0
C     --- ALPN0(I) ALPHA-EXPONENTS FOR N=0
C     --- AN1(I)   A-PARAMETERS FOR N=1
C     --- ALPN1(I) ALPHA-EXPONENTS FOR N=1
C     --- BPAR(I)  B-PARAMETERS
C     --- IERR(I)  ERROR-STATUS   FOR MMPMOD=0:
C     ---          0  OK; 1  NOT MMP-ATOM; 2 MMPTYP WRONG
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AN0(*),ALPN0(*)
      DIMENSION BPAR(*),KNG(*),KTYP(*)
      DIMENSION C(*),EX(*)
      CHARACTER*20 CDFCTR(*),SHINFO
C
      DIMENSION ITYP(10)
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
C
      IZSAVE=NUCZ
      IERR=0
C
      IZSAVE=NUCZ
      IERR=0
C
C     --- ITYP(I) IS THE TYPE (MMPTYP) OF THE I-TH BUILT-IN
C     --- MODEL-POTENTIAL FOR A GIVEN ATOM (IZSAVE)
C
      DO 40 I=1,10
         ITYP(I)=0
  40  CONTINUE
C ... ../TOOLS_MCP/MG0310:LI-NE.MPPAR
C
C     --- 
C
      IF(IZSAVE.EQ.03) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : LI/NR             CONFIGURATION : S01P00       *
C   VALENCE : (2S)              BASIS SET     : ( 4)         *
C          EFFECTIVE CORE CHARGE :  1.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =    -.196324 A.U.
C
C  ***     ORBITAL ENERGY( 2S ) =    -.1963   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  1.0D+00
          SHINFO='2S'
          NOSH= 1
          CDFCTR(1)='31'
             KNG(1)= 4
            KTYP(1)= 1
          EX( 1)= 0.92682241D+01
           C( 1)=-0.53730657D-01
          EX( 2)= 0.10744566D+01
           C( 2)=-0.21120103D+00
          EX( 3)= 0.58427760D-01
           C( 3)= 0.76560650D+00
          EX( 4)= 0.21767120D-01
           C( 4)= 0.30009762D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  1.0D+00
          NOAN0=  2
          NOAN1=  0
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.14040697D+01
          ALPN0( 1)= 0.27184225D+01
            AN0( 2)= 0.20408434D+00
          ALPN0( 2)= 0.10619056D+01
C         --- A-PARAMETERS FOR N=1
C         --- CORE-SHELL INFORMATION 
          NOSH     =  1
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.49554840D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.71057931D+03
           C(  1) = 0.20095330D-02
          EX(  2) = 0.10343675D+03
           C(  2) = 0.14823293D-01
          EX(  3) = 0.24433003D+02
           C(  3) = 0.66404805D-01
          EX(  4) = 0.73089673D+01
           C(  4) = 0.19908007D+00
          EX(  5) = 0.25068936D+01
           C(  5) = 0.38191082D+00
          EX(  6) = 0.94380031D+00
           C(  6) = 0.38846630D+00
          EX(  7) = 0.37097102D+00
           C(  7) = 0.10999668D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.04) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : BE/NR             CONFIGURATION : S02P00       *
C   VALENCE : (2S)              BASIS SET     : ( 4)         *
C          EFFECTIVE CORE CHARGE :  2.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =    -.961811 A.U.
C
C  ***     ORBITAL ENERGY( 2S ) =    -.3093   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  2.0D+00
          SHINFO='2S'
          NOSH= 1
          CDFCTR(1)='31'
             KNG(1)= 4
            KTYP(1)= 1
          EX( 1)= 0.19757152D+02
           C( 1)=-0.53813650D-01
          EX( 2)= 0.24449094D+01
           C( 2)=-0.23218774D+00
          EX( 3)= 0.13959402D+00
           C( 3)= 0.68820327D+00
          EX( 4)= 0.47059460D-01
           C( 4)= 0.38991068D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  2.0D+00
          NOAN0=  2
          NOAN1=  0
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.69860589D+00
          ALPN0( 1)= 0.51889921D+01
            AN0( 2)= 0.10611016D+00
          ALPN0( 2)= 0.20239798D+01
C         --- A-PARAMETERS FOR N=1
C         --- CORE-SHELL INFORMATION 
          NOSH     =  1
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.94653380D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.13116155D+04
           C(  1) = 0.19698940D-02
          EX(  2) = 0.19138605D+03
           C(  2) = 0.14519706D-01
          EX(  3) = 0.45278833D+02
           C(  3) = 0.65342847D-01
          EX(  4) = 0.13597320D+02
           C(  4) = 0.19826939D+00
          EX(  5) = 0.47197013D+01
           C(  5) = 0.38353068D+00
          EX(  6) = 0.18173441D+01
           C(  6) = 0.38450433D+00
          EX(  7) = 0.73401209D+00
           C(  7) = 0.10967690D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.05) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : B/NR              CONFIGURATION : S02P01       *
C   VALENCE : (2S/2P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  3.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =   -2.545778 A.U.
C
C  ***     ORBITAL ENERGY( 2S ) =    -.4947   A.U.
C  ***     ORBITAL ENERGY( 2P ) =    -.3098   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  3.0D+00
          SHINFO='2S/2P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.33505342D+02
           C( 1)=-0.58733114D-01
          EX( 2)= 0.42685415D+01
           C( 2)=-0.25074265D+00
          EX( 3)= 0.25387318D+00
           C( 3)= 0.69362350D+00
          EX( 4)= 0.82675457D-01
           C( 4)= 0.38848910D+00
          EX( 5)= 0.53938986D+01
           C( 5)= 0.43838498D-01
          EX( 6)= 0.10322568D+01
           C( 6)= 0.24924528D+00
          EX( 7)= 0.27603416D+00
           C( 7)= 0.54142626D+00
          EX( 8)= 0.80518700D-01
           C( 8)= 0.39246535D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  3.0D+00
          NOAN0=  2
          NOAN1=  0
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.50623323D+00
          ALPN0( 1)= 0.16035599D+02
            AN0( 2)= 0.93626185D-01
          ALPN0( 2)= 0.23223287D+01
C         --- A-PARAMETERS FOR N=1
C         --- CORE-SHELL INFORMATION 
          NOSH     =  1
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.15390670D+02
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.27878673D+03
           C(  1) = 0.24770105D-01
          EX(  2) = 0.38263296D+02
           C(  2) = 0.16786333D+00
          EX(  3) = 0.85763749D+01
           C(  3) = 0.49775698D+00
          EX(  4) = 0.23609779D+01
           C(  4) = 0.45427324D+00
          EX(  5) = 0.34261815D+01
           C(  5) =-0.12354426D-02
          EX(  6) = 0.29298474D+00
           C(  6) = 0.85544710D-02
          EX(  7) = 0.93751420D-01
           C(  7) = 0.69860096D-02
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.06) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE :  C/NR             CONFIGURATION : S02P02       *
C   VALENCE : (2S/2P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  4.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =   -5.335565 A.U.
C
C  ***     ORBITAL ENERGY( 2S ) =    -.7057   A.U.
C  ***     ORBITAL ENERGY( 2P ) =    -.4331   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  4.0D+00
          SHINFO='2S/2P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.36199578D+02
           C( 1)=-0.82198416D-01
          EX( 2)= 0.51796860D+01
           C( 2)=-0.25857220D+00
          EX( 3)= 0.43654500D+00
           C( 3)= 0.64703251D+00
          EX( 4)= 0.13816150D+00
           C( 4)= 0.45496290D+00
          EX( 5)= 0.84281311D+01
           C( 5)= 0.49676223D-01
          EX( 6)= 0.16548770D+01
           C( 6)= 0.26389808D+00
          EX( 7)= 0.44308380D+00
           C( 7)= 0.54247980D+00
          EX( 8)= 0.12620780D+00
           C( 8)= 0.37961949D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  4.0D+00
          NOAN0=  2
          NOAN1=  0
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.39558949D+00
          ALPN0( 1)= 0.18646050D+02
            AN0( 2)= 0.56441564D-01
          ALPN0( 2)= 0.35121000D+01
C         --- A-PARAMETERS FOR N=1
C         --- CORE-SHELL INFORMATION 
          NOSH     =  1
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.22651036D+02
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.40906987D+03
           C(  1) = 0.24466201D-01
          EX(  2) = 0.56163428D+02
           C(  2) = 0.16673090D+00
          EX(  3) = 0.12632156D+02
           C(  3) = 0.49850161D+00
          EX(  4) = 0.34953356D+01
           C(  4) = 0.45331671D+00
          EX(  5) = 0.54043974D+01
           C(  5) =-0.12974100D-02
          EX(  6) = 0.46071150D+00
           C(  6) = 0.90698702D-02
          EX(  7) = 0.14364200D+00
           C(  7) = 0.74568262D-02
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.07) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE :  N/NR             CONFIGURATION : S02P03       *
C   VALENCE : (2S/2P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  5.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =   -9.682770 A.U.
C
C  ***     ORBITAL ENERGY( 2S ) =    -.9455   A.U.
C  ***     ORBITAL ENERGY( 2P ) =    -.5673   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  5.0D+00
          SHINFO='2S/2P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.50958241D+02
           C( 1)=-0.83392980D-01
          EX( 2)= 0.73942370D+01
           C( 2)=-0.26448986D+00
          EX( 3)= 0.62445100D+00
           C( 3)= 0.64839874D+00
          EX( 4)= 0.19434270D+00
           C( 4)= 0.45511766D+00
          EX( 5)= 0.12036455D+02
           C( 5)= 0.52186787D-01
          EX( 6)= 0.24004170D+01
           C( 6)= 0.27436132D+00
          EX( 7)= 0.64243000D+00
           C( 7)= 0.54251781D+00
          EX( 8)= 0.17985000D+00
           C( 8)= 0.37163722D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  5.0D+00
          NOAN0=  2
          NOAN1=  0
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.31915672D+00
          ALPN0( 1)= 0.27639300D+02
            AN0( 2)= 0.52484011D-01
          ALPN0( 2)= 0.50221000D+01
C         --- A-PARAMETERS FOR N=1
C         --- CORE-SHELL INFORMATION 
          NOSH     =  1
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.31258000D+02
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.56350938D+03
           C(  1) = 0.24226670D-01
          EX(  2) = 0.77387535D+02
           C(  2) = 0.16571007D+00
          EX(  3) = 0.17448549D+02
           C(  3) = 0.49798873D+00
          EX(  4) = 0.48435084D+01
           C(  4) = 0.45110532D+00
          EX(  5) = 0.78209567D+01
           C(  5) = 0.13346499D-02
          EX(  6) = 0.66330180D+00
           C(  6) = 0.94109376D-02
          EX(  7) = 0.20294585D+00
           C(  7) = 0.77795411D-02
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.08) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE :  O/NR             CONFIGURATION : S02P04       *
C   VALENCE : (2S/2P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  6.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =  -15.713564 A.U.
C
C  ***     ORBITAL ENERGY( 2S ) =   -1.2442   A.U.
C  ***     ORBITAL ENERGY( 2P ) =    -.6311   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  6.0D+00
          SHINFO='2S/2P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.68223000D+02
           C( 1)= 0.83697503D-01
          EX( 2)= 0.99992500D+01
           C( 2)= 0.26941882D+00
          EX( 3)= 0.84936800D+00
           C( 3)=-0.65251375D+00
          EX( 4)= 0.26144700D+00
           C( 4)=-0.45246762D+00
          EX( 5)= 0.15515807D+02
           C( 5)= 0.55623825D-01
          EX( 6)= 0.31259460D+01
           C( 6)= 0.29047128D+00
          EX( 7)= 0.82453700D+00
           C( 7)= 0.53956068D+00
          EX( 8)= 0.22113600D+00
           C( 8)= 0.36693760D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  6.0D+00
          NOAN0=  2
          NOAN1=  0
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.28259844D+00
          ALPN0( 1)= 0.40648640D+02
            AN0( 2)= 0.36391764D-01
          ALPN0( 2)= 0.51752300D+01
C         --- A-PARAMETERS FOR N=1
C         --- CORE-SHELL INFORMATION 
          NOSH     =  1
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.41337314D+02
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.74324114D+03
           C(  1) = 0.24133060D-01
          EX(  2) = 0.10208699D+03
           C(  2) = 0.16555330D+00
          EX(  3) = 0.23058965D+02
           C(  3) = 0.49961640D+00
          EX(  4) = 0.64161860D+01
           C(  4) = 0.45179500D+00
          EX(  5) = 0.10678027D+02
           C(  5) =-0.13766600D-02
          EX(  6) = 0.90653200D+00
           C(  6) = 0.97613499D-02
          EX(  7) = 0.27406100D+00
           C(  7) = 0.79927749D-02
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.09) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE :  F/NR             CONFIGURATION : S02P05       *
C   VALENCE : (2S/2P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  7.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =  -23.975349 A.U.
C
C  ***     ORBITAL ENERGY( 2S ) =   -1.5725   A.U.
C  ***     ORBITAL ENERGY( 2P ) =    -.7299   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  7.0D+00
          SHINFO='2S/2P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.87947647D+02
           C( 1)= 0.84592289D-01
          EX( 2)= 0.12986085D+02
           C( 2)= 0.27375227D+00
          EX( 3)= 0.11057444D+01
           C( 3)=-0.65608505D+00
          EX( 4)= 0.33734000D+00
           C( 4)=-0.44982082D+00
          EX( 5)= 0.19731803D+02
           C( 5)= 0.59012237D-01
          EX( 6)= 0.40040800D+01
           C( 6)= 0.30134499D+00
          EX( 7)= 0.10487002D+01
           C( 7)= 0.53782967D+00
          EX( 8)= 0.27462770D+00
           C( 8)= 0.36167421D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  7.0D+00
          NOAN0=  2
          NOAN1=  0
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.24505135D+00
          ALPN0( 1)= 0.52027812D+02
            AN0( 2)= 0.32854788D-01
          ALPN0( 2)= 0.62063539D+01
C         --- A-PARAMETERS FOR N=1
C         --- CORE-SHELL INFORMATION 
          NOSH     =  1
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.52765500D+02
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.94761572D+03
           C(  1) = 0.24017600D-01
          EX(  2) = 0.13018410D+03
           C(  2) = 0.16511310D+00
          EX(  3) = 0.29448853D+02
           C(  3) = 0.49980400D+00
          EX(  4) = 0.82089415D+01
           C(  4) = 0.45149800D+00
          EX(  5) = 0.13975570D+02
           C(  5) =-0.14112000D-02
          EX(  6) = 0.11843901D+01
           C(  6) = 0.10059900D-01
          EX(  7) = 0.35459670D+00
           C(  7) = 0.82025000D-02
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.10) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : NE/NR             CONFIGURATION : S02P06       *
C   VALENCE : (2S/2P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  8.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =  -34.779888 A.U.
C
C  ***     ORBITAL ENERGY( 2S ) =   -1.9304   A.U.
C  ***     ORBITAL ENERGY( 2P ) =    -.8504   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  8.0D+00
          SHINFO='2S/2P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.15329425D+03
           C( 1)=-0.62040289D-01
          EX( 2)= 0.20444605D+02
           C( 2)=-0.28256022D+00
          EX( 3)= 0.12742196D+01
           C( 3)= 0.69518252D+00
          EX( 4)= 0.38758434D+00
           C( 4)= 0.39378657D+00
          EX( 5)= 0.24619842D+02
           C( 5)= 0.61768366D-01
          EX( 6)= 0.50227709D+01
           C( 6)= 0.30746955D+00
          EX( 7)= 0.13106569D+01
           C( 7)= 0.53682128D+00
          EX( 8)= 0.33828849D+00
           C( 8)= 0.35837653D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  8.0D+00
          NOAN0=  2
          NOAN1=  0
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.20441772D+00
          ALPN0( 1)= 0.54076675D+02
            AN0( 2)= 0.26360134D-01
          ALPN0( 2)= 0.76129727D+01
C         --- A-PARAMETERS FOR N=1
C         --- CORE-SHELL INFORMATION 
          NOSH     =  1
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.65544850D+02
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.11758697D+04
           C(  1) = 0.23947436D-01
          EX(  2) = 0.16156796D+03
           C(  2) = 0.16489091D+00
          EX(  3) = 0.36590690D+02
           C(  3) = 0.50015823D+00
          EX(  4) = 0.10212845D+02
           C(  4) = 0.45096867D+00
          EX(  5) = 0.17711633D+02
           C(  5) =-0.14322814D-02
          EX(  6) = 0.14973299D+01
           C(  6) = 0.10256591D-01
          EX(  7) = 0.44467001D+00
           C(  7) = 0.83527640D-02
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C--   IF (MODE.EQ.0) THEN
C--      IERR=1
C--      RETURN
C--   END IF
C
      WRITE(IW,9999) IZSAVE
      CALL ABRT
      STOP
 9999 FORMAT(/1X,'*** NO MOD-POT AVAILABLE FOR ATOM A=',I4)
      END

C*MODULE MCPLIB  *DECK MG1118
      SUBROUTINE MG1118(MMPMOD,NUCZ,MMPTYP,SHINFO,CDFCTR,NOSH,KTYP,KNG,
     2                  EX,C,ZCORE,NOAN0,NOAN1,AN0,ALPN0,AN1,
     3                  ALPN1,BPAR,ITYP,IERR)
C
C     --- MMPMOD: 0  READ VALENCE BASIS-SET
C     ---         1  READ MODEL-POTENTIAL-PARAMETERS
C     ---            AND CORE-SHELL BASIS-SET
C
C     --- NUCZ:    NUCLEAR CHARGE
C     --- MMPTYP:  1  NON-RELATIVISTIC   MOD.POT.
C     ---          2  QUASI-RELATIVISTIC MOD.POT.
C     ---          3  READ IN MOD.POT.FROM $MMP GROUP
C     ---          4  TYPE NOT SPECIFIED
C     ---    W*10**I  FOR I-TH SHELL (I>=1)
C     ---    W: 1,2,3,4,6  S,P,D,F,L - SHELL
C     --- CDFCTR(I) CODED DEFAULT CONTRACTION OF VALENCE-SHELL I
C     --- NOSH     NO.OF (VALENCE OR CORE) SHELLS
C     --- KTYP(I)  TYPE OF SHELL I
C     --- KNG(I)   NO.OF GTFS IN SHELL I
C     --- EX(J)    EXPONENT OF J-TH GTF
C     --- C(J)     CONTRACTION COEFF.OF J-TH GTF
C     --- ZCORE    CORE-CHARGE
C     --- NOAN0    NO.OF A-TERMS FOR N=0
C     --- NOAN1    NO.OF A-TERMS FOR N=1
C     --- AN0(I)   A-PARAMETERS FOR N=0
C     --- ALPN0(I) ALPHA-EXPONENTS FOR N=0
C     --- AN1(I)   A-PARAMETERS FOR N=1
C     --- ALPN1(I) ALPHA-EXPONENTS FOR N=1
C     --- BPAR(I)  B-PARAMETERS
C     --- IERR(I)  ERROR-STATUS   FOR MMPMOD=0:
C     ---          0  OK; 1  NOT MMP-ATOM; 2 MMPTYP WRONG
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AN0(*),ALPN0(*),AN1(*),ALPN1(*)
      DIMENSION BPAR(*),KNG(*),KTYP(*)
      DIMENSION C(*),EX(*)
      CHARACTER*20 CDFCTR(*),SHINFO
C
      DIMENSION ITYP(10)
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
C
      IZSAVE=NUCZ
      IERR=0
C
      IZSAVE=NUCZ
      IERR=0
C
C     --- ITYP(I) IS THE TYPE (MMPTYP) OF THE I-TH BUILT-IN
C     --- MODEL-POTENTIAL FOR A GIVEN ATOM (IZSAVE)
C
      DO 40 I=1,10
         ITYP(I)=0
  40  CONTINUE
C ... ../TOOLS_MCP/MG1118:NA-AR.MPPAR
C
C     --- 
C
      IF(IZSAVE.EQ.11) THEN
        ITYP(1)=112
        IF (MMPTYP.EQ.112) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : NA/NR             CONFIGURATION : P06S01       *
C   VALENCE : (3S/2P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  7.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =  -27.702102 A.U.
C
C  ***     ORBITAL ENERGY( 3S ) =    -.1821   A.U.
C  ***     ORBITAL ENERGY( 2P ) =   -1.5181   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  7.0D+00
          SHINFO='3S/2P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='4'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.45183832D+02
           C( 1)= 0.35464597D-01
          EX( 2)= 0.91956191D+00
           C( 2)=-0.25355735D+00
          EX( 3)= 0.46387630D-01
           C( 3)= 0.85351074D+00
          EX( 4)= 0.16546950D-01
           C( 4)= 0.21154996D+00
          EX( 5)= 0.34506121D+02
           C( 5)= 0.55903870D-01
          EX( 6)= 0.71463353D+01
           C( 6)= 0.28368743D+00
          EX( 7)= 0.19283157D+01
           C( 7)= 0.54056417D+00
          EX( 8)= 0.53448009D+00
           C( 8)= 0.36399836D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  7.0D+00
          NOAN0=  3
          NOAN1=  3
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.24903555D+00
          ALPN0( 1)= 0.30506560D+03
            AN0( 2)= 0.21616059D+00
          ALPN0( 2)= 0.14646747D+02
            AN0( 3)= 0.11485600D+00
          ALPN0( 3)= 0.21368669D+01
C         --- A-PARAMETERS FOR N=1
            AN1( 1)= 0.22725132D+02
          ALPN1( 1)= 0.47184272D+03
            AN1( 2)= 0.64983925D+00
          ALPN1( 2)= 0.20523771D+03
            AN1( 3)= 0.69301799D-01
          ALPN1( 3)= 0.68946093D+00
C         --- CORE-SHELL INFORMATION 
          NOSH     =  2
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.80957000D+02
          KTYP( 2) =  1
           KNG( 2) =  7
          BPAR( 2) = 0.55940000D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.14373868D+04
           C(  1) = 0.23716362D-01
          EX(  2) = 0.19752425D+03
           C(  2) = 0.16365706D+00
          EX(  3) = 0.44777665D+02
           C(  3) = 0.49871140D+00
          EX(  4) = 0.12518309D+02
           C(  4) = 0.45296990D+00
          EX(  5) = 0.21567240D+02
           C(  5) =-0.15394690D-02
          EX(  6) = 0.19100396D+01
           C(  6) = 0.10624124D-01
          EX(  7) = 0.60717102D+00
           C(  7) = 0.85477456D-02
          EX(  8) = 0.14373868D+04
           C(  8) = 0.62947990D-02
          EX(  9) = 0.19752425D+03
           C(  9) = 0.43437874D-01
          EX( 10) = 0.44777665D+02
           C( 10) = 0.13236802D+00
          EX( 11) = 0.12518309D+02
           C( 11) = 0.12022731D+00
          EX( 12) = 0.21567240D+02
           C( 12) = 0.89529738D-01
          EX( 13) = 0.19100396D+01
           C( 13) =-0.61785918D+00
          EX( 14) = 0.60717102D+00
           C( 14) =-0.49710482D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.12) THEN
        ITYP(1)=112
        IF (MMPTYP.EQ.112) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : MG/NR             CONFIGURATION : P06S02       *
C   VALENCE : (3S/2P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  8.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =  -38.216852 A.U.
C
C  ***     ORBITAL ENERGY( 3S ) =    -.2530   A.U.
C  ***     ORBITAL ENERGY( 2P ) =   -2.2822   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  8.0D+00
          SHINFO='3S/2P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='4'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.52378743D+02
           C( 1)= 0.47038043D-01
          EX( 2)= 0.12877140D+01
           C( 2)=-0.31349585D+00
          EX( 3)= 0.87692270D-01
           C( 3)= 0.76695487D+00
          EX( 4)= 0.31454784D-01
           C( 4)= 0.32201052D+00
          EX( 5)= 0.45327095D+02
           C( 5)= 0.50914639D-01
          EX( 6)= 0.94830887D+01
           C( 6)= 0.27406187D+00
          EX( 7)= 0.26137933D+01
           C( 7)= 0.54400336D+00
          EX( 8)= 0.75515353D+00
           C( 8)= 0.36187855D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  8.0D+00
          NOAN0=  3
          NOAN1=  3
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.21170786D+00
          ALPN0( 1)= 0.33718154D+03
            AN0( 2)= 0.17784947D+00
          ALPN0( 2)= 0.15822919D+02
            AN0( 3)= 0.93981276D-01
          ALPN0( 3)= 0.23093307D+01
C         --- A-PARAMETERS FOR N=1
            AN1( 1)= 0.19770101D+02
          ALPN1( 1)= 0.55524754D+03
            AN1( 2)= 0.61954714D+00
          ALPN1( 2)= 0.27979393D+03
            AN1( 3)= 0.56876624D-01
          ALPN1( 3)= 0.76813234D+00
C         --- CORE-SHELL INFORMATION 
          NOSH     =  2
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.98063470D+02
          KTYP( 2) =  1
           KNG( 2) =  7
          BPAR( 2) = 0.75354400D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.17265765D+04
           C(  1) = 0.23503634D-01
          EX(  2) = 0.23731223D+03
           C(  2) = 0.16247869D+00
          EX(  3) = 0.53846608D+02
           C(  3) = 0.49722124D+00
          EX(  4) = 0.15078842D+02
           C(  4) = 0.45499774D+00
          EX(  5) = 0.25865162D+02
           C(  5) =-0.16385139D-02
          EX(  6) = 0.23752961D+01
           C(  6) = 0.11025669D-01
          EX(  7) = 0.79379808D+00
           C(  7) = 0.86412833D-02
          EX(  8) = 0.17265765D+04
           C(  8) = 0.64661553D-02
          EX(  9) = 0.23731223D+03
           C(  9) = 0.44699999D-01
          EX( 10) = 0.53846608D+02
           C( 10) = 0.13679203D+00
          EX( 11) = 0.15078842D+02
           C( 11) = 0.12517579D+00
          EX( 12) = 0.25865162D+02
           C( 12) = 0.92846212D-01
          EX( 13) = 0.23752961D+01
           C( 13) =-0.62476832D+00
          EX( 14) = 0.79379808D+00
           C( 14) =-0.48965738D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.13) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : AL/NR             CONFIGURATION : S02P01       *
C   VALENCE : (3S/3P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  3.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =   -1.879478 A.U.
C
C  ***     ORBITAL ENERGY( 3S ) =    -.3934   A.U.
C  ***     ORBITAL ENERGY( 3P ) =    -.2099   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  3.0D+00
          SHINFO='3S/3P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.61013425D+02
           C( 1)=-0.71088671D-01
          EX( 2)= 0.16677184D+01
           C( 2)= 0.39365648D+00
          EX( 3)= 0.14178169D+00
           C( 3)=-0.78872575D+00
          EX( 4)= 0.50957100D-01
           C( 4)=-0.31891433D+00
          EX( 5)= 0.16738095D+02
           C( 5)= 0.47854721D-01
          EX( 6)= 0.30480537D+01
           C( 6)= 0.12752243D+00
          EX( 7)= 0.17760526D+00
           C( 7)=-0.55703526D+00
          EX( 8)= 0.51052100D-01
           C( 8)=-0.55062242D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  3.0D+00
          NOAN0=  3
          NOAN1=  3
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.11286709D+01
          ALPN0( 1)= 0.13169528D+03
            AN0( 2)= 0.20155823D+01
          ALPN0( 2)= 0.52173240D+01
            AN0( 3)= 0.12157937D+00
          ALPN0( 3)= 0.23314592D+01
C         --- A-PARAMETERS FOR N=1
            AN1( 1)= 0.38022861D+02
          ALPN1( 1)= 0.91422726D+03
            AN1( 2)= 0.14026833D+02
          ALPN1( 2)= 0.20775601D+03
            AN1( 3)= 0.11166382D+00
          ALPN1( 3)= 0.12510797D+01
C         --- CORE-SHELL INFORMATION 
          NOSH     =  3
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.76051300D+02
          KTYP( 2) =  1
           KNG( 2) =  7
          BPAR( 2) = 0.63839100D+01
          KTYP( 3) =  2
           KNG( 3) =  5
          BPAR( 3) = 0.41837900D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.20438076D+04
           C(  1) = 0.23306138D-01
          EX(  2) = 0.28093674D+03
           C(  2) = 0.16139194D+00
          EX(  3) = 0.63789115D+02
           C(  3) = 0.49584335D+00
          EX(  4) = 0.17888624D+02
           C(  4) = 0.45690622D+00
          EX(  5) = 0.30619272D+02
           C(  5) =-0.17238765D-02
          EX(  6) = 0.28868382D+01
           C(  6) = 0.11427572D-01
          EX(  7) = 0.99935005D+00
           C(  7) = 0.86598621D-02
          EX(  8) = 0.20438076D+04
           C(  8) = 0.66039049D-02
          EX(  9) = 0.28093674D+03
           C(  9) = 0.45731176D-01
          EX( 10) = 0.63789115D+02
           C( 10) = 0.14049957D+00
          EX( 11) = 0.17888624D+02
           C( 11) = 0.12946655D+00
          EX( 12) = 0.30619272D+02
           C( 12) = 0.95662365D-01
          EX( 13) = 0.28868382D+01
           C( 13) =-0.63414553D+00
          EX( 14) = 0.99935005D+00
           C( 14) =-0.48055816D+00
          EX( 15) = 0.10896206D+03
           C( 15) = 0.18243123D-01
          EX( 16) = 0.23630353D+02
           C( 16) = 0.12121048D+00
          EX( 17) = 0.71077322D+01
           C( 17) = 0.35717809D+00
          EX( 18) = 0.23660639D+01
           C( 18) = 0.49939064D+00
          EX( 19) = 0.80250185D+00
           C( 19) = 0.23093033D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.14) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : SI/NR             CONFIGURATION : S02P02       *
C   VALENCE : (3S/3P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  4.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =   -3.717317 A.U.
C
C  ***     ORBITAL ENERGY( 3S ) =    -.5399   A.U.
C  ***     ORBITAL ENERGY( 3P ) =    -.2971   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  4.0D+00
          SHINFO='3S/3P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.70678910D+02
           C( 1)= 0.57145105D-01
          EX( 2)= 0.20763480D+01
           C( 2)=-0.42826185D+00
          EX( 3)= 0.20335900D+00
           C( 3)= 0.45774670D+00
          EX( 4)= 0.13816150D+00
           C( 4)= 0.63104883D+00
          EX( 5)= 0.21190060D+02
           C( 5)=-0.51895309D-01
          EX( 6)= 0.39626990D+01
           C( 6)=-0.16851947D+00
          EX( 7)= 0.26056700D+00
           C( 7)= 0.60674124D+00
          EX( 8)= 0.76803329D-01
           C( 8)= 0.49612106D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  4.0D+00
          NOAN0=  3
          NOAN1=  3
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.83001083D+00
          ALPN0( 1)= 0.13498630D+03
            AN0( 2)= 0.14754701D+01
          ALPN0( 2)= 0.43357137D+01
            AN0( 3)=-0.86341681D-01
          ALPN0( 3)= 0.45479023D+01
C         --- A-PARAMETERS FOR N=1
            AN1( 1)= 0.30456884D+02
          ALPN1( 1)= 0.85406377D+03
            AN1( 2)= 0.95568910D+01
          ALPN1( 2)= 0.21197610D+03
            AN1( 3)= 0.86513019D-01
          ALPN1( 3)= 0.99809167D+00
C         --- CORE-SHELL INFORMATION 
          NOSH     =  3
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.82574880D+02
          KTYP( 2) =  1
           KNG( 2) =  7
          BPAR( 2) = 0.73878360D+01
          KTYP( 3) =  2
           KNG( 3) =  5
          BPAR( 3) = 0.51072648D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.23879574D+04
           C(  1) = 0.23137295D-01
          EX(  2) = 0.32828020D+03
           C(  2) = 0.16045018D+00
          EX(  3) = 0.74586123D+02
           C(  3) = 0.49462158D+00
          EX(  4) = 0.20944123D+02
           C(  4) = 0.45858021D+00
          EX(  5) = 0.35822223D+02
           C(  5) =-0.17964583D-02
          EX(  6) = 0.34446205D+01
           C(  6) = 0.11811096D-01
          EX(  7) = 0.12222835D+01
           C(  7) = 0.86220458D-02
          EX(  8) = 0.23879574D+04
           C(  8) = 0.67191789D-02
          EX(  9) = 0.32828020D+03
           C(  9) = 0.46595484D-01
          EX( 10) = 0.74586123D+02
           C( 10) = 0.14364042D+00
          EX( 11) = 0.20944123D+02
           C( 11) = 0.13317384D+00
          EX( 12) = 0.35822223D+02
           C( 12) = 0.98068957D-01
          EX( 13) = 0.34446205D+01
           C( 13) =-0.64476972D+00
          EX( 14) = 0.12222835D+01
           C( 14) =-0.47067892D+00
          EX( 15) = 0.13482278D+03
           C( 15) = 0.17311800D-01
          EX( 16) = 0.29361481D+02
           C( 16) = 0.11697476D+00
          EX( 17) = 0.88906349D+01
           C( 17) = 0.35318628D+00
          EX( 18) = 0.30003493D+01
           C( 18) = 0.50346337D+00
          EX( 19) = 0.10376018D+01
           C( 19) = 0.22874705D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.15) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE :  P/NR             CONFIGURATION : S02P03       *
C   VALENCE : (3S/3P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  5.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =   -6.335286 A.U.
C
C  ***     ORBITAL ENERGY( 3S ) =    -.6964   A.U.
C  ***     ORBITAL ENERGY( 3P ) =    -.3917   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  5.0D+00
          SHINFO='3S/3P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.81256103D+02
           C( 1)= 0.35183989D-01
          EX( 2)= 0.25178780D+01
           C( 2)=-0.43408372D+00
          EX( 3)= 0.27338900D+00
           C( 3)= 0.76857011D+00
          EX( 4)= 0.10006351D+00
           C( 4)= 0.36680263D+00
          EX( 5)= 0.26078700D+02
           C( 5)=-0.45444009D-01
          EX( 6)= 0.49744954D+01
           C( 6)=-0.15889781D+00
          EX( 7)= 0.35323726D+00
           C( 7)= 0.60566870D+00
          EX( 8)= 0.10560026D+00
           C( 8)= 0.49860766D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  5.0D+00
          NOAN0=  3
          NOAN1=  3
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.54765153D+00
          ALPN0( 1)= 0.12699655D+03
            AN0( 2)= 0.11426965D+01
          ALPN0( 2)= 0.51250394D+01
            AN0( 3)=-0.67031639D-01
          ALPN0( 3)= 0.28490414D+01
C         --- A-PARAMETERS FOR N=1
            AN1( 1)= 0.15647806D+02
          ALPN1( 1)= 0.11691293D+04
            AN1( 2)= 0.51897151D+01
          ALPN1( 2)= 0.18439113D+03
            AN1( 3)= 0.47283639D-01
          ALPN1( 3)= 0.13267763D+01
C         --- CORE-SHELL INFORMATION 
          NOSH     =  3
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.87966687D+02
          KTYP( 2) =  1
           KNG( 2) =  7
          BPAR( 2) = 0.82622100D+01
          KTYP( 3) =  2
           KNG( 3) =  5
          BPAR( 3) = 0.59410450D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.27593369D+04
           C(  1) = 0.22988443D-01
          EX(  2) = 0.37937300D+03
           C(  2) = 0.15961617D+00
          EX(  3) = 0.86241959D+02
           C(  3) = 0.49352883D+00
          EX(  4) = 0.24245952D+02
           C(  4) = 0.46007652D+00
          EX(  5) = 0.41467939D+02
           C(  5) =-0.18587396D-02
          EX(  6) = 0.40493273D+01
           C(  6) = 0.12168495D-01
          EX(  7) = 0.14620909D+01
           C(  7) = 0.85539149D-02
          EX(  8) = 0.27593369D+04
           C(  8) = 0.68154743D-02
          EX(  9) = 0.37937300D+03
           C(  9) = 0.47322036D-01
          EX( 10) = 0.86241959D+02
           C( 10) = 0.14631844D+00
          EX( 11) = 0.24245952D+02
           C( 11) = 0.13640070D+00
          EX( 12) = 0.41467939D+02
           C( 12) = 0.10014081D+00
          EX( 13) = 0.40493273D+01
           C( 13) =-0.65558566D+00
          EX( 14) = 0.14620909D+01
           C( 14) =-0.46084775D+00
          EX( 15) = 0.16357377D+03
           C( 15) = 0.16532835D-01
          EX( 16) = 0.35742036D+02
           C( 16) = 0.11332407D+00
          EX( 17) = 0.10881465D+02
           C( 17) = 0.34933115D+00
          EX( 18) = 0.37132893D+01
           C( 18) = 0.50656661D+00
          EX( 19) = 0.13040087D+01
           C( 19) = 0.22820296D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.16) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE :  S/NR             CONFIGURATION : S02P04       *
C   VALENCE : (3S/3P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  6.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =   -9.947196 A.U.
C
C  ***     ORBITAL ENERGY( 3S ) =    -.8795   A.U.
C  ***     ORBITAL ENERGY( 3P ) =    -.4376   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  6.0D+00
          SHINFO='3S/3P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.92693689D+02
           C( 1)= 0.26349969D-01
          EX( 2)= 0.29953266D+01
           C( 2)=-0.46428207D+00
          EX( 3)= 0.35397237D+00
           C( 3)= 0.76997684D+00
          EX( 4)= 0.12962166D+00
           C( 4)= 0.37724028D+00
          EX( 5)= 0.31777309D+02
           C( 5)=-0.56289957D-01
          EX( 6)= 0.61605477D+01
           C( 6)=-0.18475972D+00
          EX( 7)= 0.45462028D+00
           C( 7)= 0.61404369D+00
          EX( 8)= 0.13024690D+00
           C( 8)= 0.49406578D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  6.0D+00
          NOAN0=  3
          NOAN1=  3
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.60336902D+00
          ALPN0( 1)= 0.11559165D+03
            AN0( 2)= 0.95621820D+00
          ALPN0( 2)= 0.59711412D+01
            AN0( 3)=-0.78411783D-01
          ALPN0( 3)= 0.28590902D+01
C         --- A-PARAMETERS FOR N=1
            AN1( 1)= 0.62934707D+01
          ALPN1( 1)= 0.20626801D+04
            AN1( 2)= 0.57199110D+01
          ALPN1( 2)= 0.14871381D+03
            AN1( 3)= 0.38138852D-01
          ALPN1( 3)= 0.16440104D+01
C         --- CORE-SHELL INFORMATION 
          NOSH     =  3
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.96604672D+02
          KTYP( 2) =  1
           KNG( 2) =  7
          BPAR( 2) = 0.94545045D+01
          KTYP( 3) =  2
           KNG( 3) =  5
          BPAR( 3) = 0.70166355D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.31580329D+04
           C(  1) = 0.22855723D-01
          EX(  2) = 0.43422727D+03
           C(  2) = 0.15886937D+00
          EX(  3) = 0.98759464D+02
           C(  3) = 0.49254160D+00
          EX(  4) = 0.27795217D+02
           C(  4) = 0.46142784D+00
          EX(  5) = 0.47548414D+02
           C(  5) =-0.19129356D-02
          EX(  6) = 0.47034445D+01
           C(  6) = 0.12489837D-01
          EX(  7) = 0.17201378D+01
           C(  7) = 0.84759572D-02
          EX(  8) = 0.31580329D+04
           C(  8) = 0.68973859D-02
          EX(  9) = 0.43422727D+03
           C(  9) = 0.47943503D-01
          EX( 10) = 0.98759464D+02
           C( 10) = 0.14863891D+00
          EX( 11) = 0.27795217D+02
           C( 11) = 0.13924942D+00
          EX( 12) = 0.47548414D+02
           C( 12) = 0.10196352D+00
          EX( 13) = 0.47034445D+01
           C( 13) =-0.66573474D+00
          EX( 14) = 0.17201378D+01
           C( 14) =-0.45178646D+00
          EX( 15) = 0.19310309D+03
           C( 15) = 0.16156360D-01
          EX( 16) = 0.42295572D+02
           C( 16) = 0.11184668D+00
          EX( 17) = 0.12926176D+02
           C( 17) = 0.34919431D+00
          EX( 18) = 0.44463082D+01
           C( 18) = 0.50941853D+00
          EX( 19) = 0.15753638D+01
           C( 19) = 0.22319216D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.17) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : CL/NR             CONFIGURATION : S02P05       *
C   VALENCE : (3S/3P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  7.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =  -14.816523 A.U.
C
C  ***     ORBITAL ENERGY( 3S ) =   -1.0893   A.U.
C  ***     ORBITAL ENERGY( 3P ) =    -.5016   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  7.0D+00
          SHINFO='3S/3P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.10497889D+03
           C( 1)= 0.11234322D-01
          EX( 2)= 0.35078200D+01
           C( 2)=-0.52552059D+00
          EX( 3)= 0.44380400D+00
           C( 3)= 0.79309396D+00
          EX( 4)= 0.16275300D+00
           C( 4)= 0.36478934D+00
          EX( 5)= 0.37786850D+02
           C( 5)=-0.67877953D-01
          EX( 6)= 0.74179980D+01
           C( 6)=-0.22000177D+00
          EX( 7)= 0.56676530D+00
           C( 7)= 0.61435686D+00
          EX( 8)= 0.16001640D+00
           C( 8)= 0.49137659D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  7.0D+00
          NOAN0=  3
          NOAN1=  3
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.37923915D+00
          ALPN0( 1)= 0.11205427D+03
            AN0( 2)= 0.81531801D+00
          ALPN0( 2)= 0.68609777D+01
            AN0( 3)=-0.93121468D-01
          ALPN0( 3)= 0.27894290D+01
C         --- A-PARAMETERS FOR N=1
            AN1( 1)= 0.50228932D+01
          ALPN1( 1)= 0.16675256D+04
            AN1( 2)= 0.46924079D+01
          ALPN1( 2)= 0.10706773D+03
            AN1( 3)= 0.14140620D+00
          ALPN1( 3)= 0.35745528D+01
C         --- CORE-SHELL INFORMATION 
          NOSH     =  3
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.10488440D+03
          KTYP( 2) =  1
           KNG( 2) =  7
          BPAR( 2) = 0.10607500D+02
          KTYP( 3) =  2
           KNG( 3) =  5
          BPAR( 3) = 0.80722000D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.35842922D+04
           C(  1) = 0.22734786D-01
          EX(  2) = 0.49287736D+03
           C(  2) = 0.15818514D+00
          EX(  3) = 0.11214646D+03
           C(  3) = 0.49162402D+00
          EX(  4) = 0.31593867D+02
           C(  4) = 0.46268289D+00
          EX(  5) = 0.54066098D+02
           C(  5) =-0.19616783D-02
          EX(  6) = 0.54056440D+01
           C(  6) = 0.12791593D-01
          EX(  7) = 0.19951056D+01
           C(  7) = 0.83992448D-02
          EX(  8) = 0.35842922D+04
           C(  8) = 0.69661986D-02
          EX(  9) = 0.49287736D+03
           C(  9) = 0.48469738D-01
          EX( 10) = 0.11214646D+03
           C( 10) = 0.15063923D+00
          EX( 11) = 0.31593867D+02
           C( 11) = 0.14177133D+00
          EX( 12) = 0.54066098D+02
           C( 12) = 0.10355442D+00
          EX( 13) = 0.54056440D+01
           C( 13) =-0.67525142D+00
          EX( 14) = 0.19951056D+01
           C( 14) =-0.44338510D+00
          EX( 15) = 0.22390776D+03
           C( 15) = 0.15961660D-01
          EX( 16) = 0.49160411D+02
           C( 16) = 0.11122717D+00
          EX( 17) = 0.15077593D+02
           C( 17) = 0.35014972D+00
          EX( 18) = 0.52222618D+01
           C( 18) = 0.51156468D+00
          EX( 19) = 0.18624703D+01
           C( 19) = 0.21760594D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C
C     --- 
C
      IF(IZSAVE.EQ.18) THEN
        ITYP(1)=111
        IF (MMPTYP.EQ.111) THEN
          IF (MMPMOD.EQ.0) THEN
C
C*************************************************************
C     TITLE : AR/NR             CONFIGURATION : S02P06       *
C   VALENCE : (3S/3P)           BASIS SET     : ( 4/ 4)      *
C          EFFECTIVE CORE CHARGE :  8.000000                 *
C*************************************************************
C
C  ***     TOTAL ENERGY         =  -21.579774 A.U.
C
C  ***     ORBITAL ENERGY( 3S ) =   -1.2778   A.U.
C  ***     ORBITAL ENERGY( 3P ) =    -.5910   A.U.
C
C         --- VALENCE-SHELL BASIS-SET
C
          ZCORE = IZSAVE -  8.0D+00
          SHINFO='3S/3P'
          NOSH= 2
          CDFCTR(1)='31'
          CDFCTR(2)='31'
             KNG(1)= 4
            KTYP(1)= 1
             KNG(2)= 4
            KTYP(2)= 2
          EX( 1)= 0.11811521D+03
           C( 1)=-0.24944500D-03
          EX( 2)= 0.40558648D+01
           C( 2)=-0.57508113D+00
          EX( 3)= 0.54293853D+00
           C( 3)= 0.79943743D+00
          EX( 4)= 0.19926562D+00
           C( 4)= 0.36611905D+00
          EX( 5)= 0.44222080D+02
           C( 5)=-0.97624207D-01
          EX( 6)= 0.87705467D+01
           C( 6)=-0.29878046D+00
          EX( 7)= 0.69053744D+00
           C( 7)= 0.62127700D+00
          EX( 8)= 0.19424861D+00
           C( 8)= 0.46517595D+00
          RETURN
        ELSE
C
C         --- CORE-SHELLS AND MODEL-POTENTIAL
C
          ZCORE = IZSAVE -  8.0D+00
          NOAN0=  3
          NOAN1=  3
C         --- A-PARAMETERS FOR N=0
            AN0( 1)= 0.38157243D+00
          ALPN0( 1)= 0.10229371D+03
            AN0( 2)= 0.66740068D+00
          ALPN0( 2)= 0.72540802D+01
            AN0( 3)=-0.60145792D-01
          ALPN0( 3)= 0.35391335D+01
C         --- A-PARAMETERS FOR N=1
            AN1( 1)= 0.38798637D+01
          ALPN1( 1)= 0.19196248D+04
            AN1( 2)= 0.42267664D+01
          ALPN1( 2)= 0.92372859D+02
            AN1( 3)= 0.12529355D+00
          ALPN1( 3)= 0.33077683D+01
C         --- CORE-SHELL INFORMATION 
          NOSH     =  3
          KTYP( 1) =  1
           KNG( 1) =  7
          BPAR( 1) = 0.11861000D+03
          KTYP( 2) =  1
           KNG( 2) =  7
          BPAR( 2) = 0.12322200D+02
          KTYP( 3) =  2
           KNG( 3) =  5
          BPAR( 3) = 0.95714700D+01
C         --- EXPONENTS AND COEFFICIENTS 
          EX(  1) = 0.40378052D+04
           C(  1) = 0.22626543D-01
          EX(  2) = 0.55528005D+03
           C(  2) = 0.15757120D+00
          EX(  3) = 0.12639300D+03
           C(  3) = 0.49079641D+00
          EX(  4) = 0.35639049D+02
           C(  4) = 0.46381464D+00
          EX(  5) = 0.61018190D+02
           C(  5) =-0.20044308D-02
          EX(  6) = 0.61564043D+01
           C(  6) = 0.13063236D-01
          EX(  7) = 0.22871154D+01
           C(  7) = 0.83215279D-02
          EX(  8) = 0.40378052D+04
           C(  8) = 0.70252355D-02
          EX(  9) = 0.55528005D+03
           C(  9) = 0.48923726D-01
          EX( 10) = 0.12639300D+03
           C( 10) = 0.15238565D+00
          EX( 11) = 0.35639049D+02
           C( 11) = 0.14400817D+00
          EX( 12) = 0.61018190D+02
           C( 12) = 0.10495532D+00
          EX( 13) = 0.61564043D+01
           C( 13) =-0.68401272D+00
          EX( 14) = 0.22871154D+01
           C( 14) =-0.43572903D+00
          EX( 15) = 0.26372414D+03
           C( 15) = 0.15147014D-01
          EX( 16) = 0.58007550D+02
           C( 16) = 0.10692360D+00
          EX( 17) = 0.17853637D+02
           C( 17) = 0.34336629D+00
          EX( 18) = 0.62306300D+01
           C( 18) = 0.51296422D+00
          EX( 19) = 0.22498520D+01
           C( 19) = 0.22438203D+00
          RETURN
          END IF
        END IF
        IERR=2
        RETURN
      END IF
C ...
C
C--   IF (MODE.EQ.0) THEN
C--      IERR=1
C--      RETURN
C--   END IF
C
      WRITE(IW,9999) IZSAVE
      CALL ABRT
      STOP
 9999 FORMAT(/1X,'*** NO MOD-POT AVAILABLE FOR ATOM A=',I4)
      END
