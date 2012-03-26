      PROGRAM GET_STO6G
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                  get_gen_sto6g.f - a VB2000 utility
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     This program outputs a GEN basis set file for general STO-6G basis 
C     sets for use by VB2000. It thus allows the use of different exponents 
C     for s and p functions, best atom exponents or indeed any exponents of 
C     your choice. Input is:-
C
C     N (Atom atomic number) followed optionally by " *" for polarisation 
C     d (or p on H) to be added, and the polarisation exponent if different 
C     from the default described below (in format I2,A2,F10.0)
C     List of required exponents (Free format) -
C        1 for H, He; 3 for Li - Ne; 5 for Na-Ar; 8 for K - Kr.
C        If the 8th exponent for D is negative, the D functions
C        are omitted (usefull for K and Ca).
C
C     The polarisation functions are the ones described in GAMESS(US) as 
C     standard.
C
C     Repeated for each required atom until N = 0.
C
C     To compile the program on linux, type:
C
C     gfortran -o get_gen_sto6g.x get_gen_sto6g.f
C
C     at the prompt and then run it as:
C
C     ./get_gen_sto6g.x < data.inp >& output.bas &
C
C     where data.inp is a data file prepared as above and output.bas
C     is an output file to be used in VB2000 as a general basis file.
C
C     Written by Brian J. Duke, December 2011, to get data that matches
C     the Roso VB program used by Richard Harcourt.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      DOUBLE PRECISION EX(8,6),C(8,6),XP(8),POLF1(36),CC,CCEXP
      INTEGER I,J,IATOM,N,M
      LOGICAL LSTAR,L
      CHARACTER*2 ATMSYM,STAR,ORBSYM
      DIMENSION ATMSYM(36),ORBSYM(8)
      DATA ATMSYM/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     1            'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     2            ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co',
     3            'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'/
      DATA STAR/'  '/
      DATA ORBSYM/' S',' S',' P',' S',' P',' S',' P',' D'/
C
C    ---- "COMMON" POLARIZATION SET FROM GAMESS(UK) ----
C    THESE VALUES ARE INSPIRED BY POPLE BASIS SETS FOR THE LIGHTER
C    ELEMENTS, AND BY THE HUZINAGA "GREEN BOOK" FOR THE HEAVIER ATOMS.
C
      DATA POLF1    /1.100D+00,1.100D+00,
     *   0.200D+00,0.400D+00,
     *   0.600D+00,0.800D+00,0.800D+00,0.800D+00,0.800D+00,0.800D+00,
     *   0.175D+00,0.175D+00,
     *   0.325D+00,0.395D+00,0.550D+00,0.650D+00,0.750D+00,0.850D+00,
     *   0.200D+00,0.200D+00,    10*0.0D+00,
     *   0.207D+00,0.246D+00,0.293D+00,0.338D+00,0.389D+00,0.443D+00/
      CALL GETCEX(EX,C)
 1    READ(5,'(I2,A2,F10.0)') IATOM,STAR,CCEXP
      LSTAR=.FALSE.
      IF(STAR.EQ." *") LSTAR=.TRUE.
      IF (IATOM.GT.20.AND.IATOM.LE.31) LSTAR=.FALSE.
C     No polarisation for transition metals
      IF (IATOM.EQ.0) GOTO 3
      IF (IATOM.GT.36) THEN
        WRITE(6,'(" ATOM NO. TOO HIGH. IATOM=",I2)') IATOM
        GOTO 3
      ENDIF
      IF (IATOM.LE.2) THEN
        READ(5,*) XP(1)
        N=1
      ELSE IF (IATOM.LE.10) THEN
        READ(5,*) XP(1),XP(2),XP(3)
        N=3
      ELSE IF (IATOM.LE.18) THEN
        READ(5,*) XP(1),XP(2),XP(3),XP(4),XP(5)
        N=5
      ELSE IF (IATOM.LE.36) THEN
        READ(5,*) XP(1),XP(2),XP(3),XP(4),XP(5),
     &   XP(6),XP(7),XP(8)
        N=8
      ENDIF
      WRITE(6,'("-",A2,"  0")') ATMSYM(IATOM)
      DO J=1,N
        IF (J.NE.8.OR.XP(8).GT.0.0D-00) THEN
          WRITE(6,'("     1    ",A2,9X,"6")') ORBSYM(J)
          DO I=1,6
           WRITE(6,'(I5,1X,F14.7,F15.10)') I,EX(J,I)*XP(J)*XP(J),C(J,I)
          ENDDO
        ENDIF
      ENDDO
      IF(LSTAR) THEN
        CC=1.0D00
        IF(CCEXP.LT.1.0D-6) CCEXP=POLF1(IATOM)
        I=1
        M=8
        IF(IATOM.LE.2) M=3
        WRITE(6,'("     1    ",A2,9X,"1")') ORBSYM(M)
        WRITE(6,'(I5,1X,F14.7,F15.10)') I,CCEXP,CC
      ENDIF
      WRITE(6,'("****")')
      GOTO 1
 3    CONTINUE
      END
      SUBROUTINE GETCEX(EX,C)
      IMPLICIT NONE
      DOUBLE PRECISION EX(8,6),C(8,6)
C  1S
      EX(1,1)=2.310303149D+01 
      C(1,1)=9.163596280D-03
      EX(1,2)=4.235915534D+00 
      C(1,2)=4.936149294D-02
      EX(1,3)=1.185056519D+00 
      C(1,3)=1.685383049D-01
      EX(1,4)=4.070988982D-01 
      C(1,4)=3.705627997D-01
      EX(1,5)=1.580884151D-01 
      C(1,5)=4.164915298D-01
      EX(1,6)=6.510953954D-02 
      C(1,6)=1.303340841D-01
C  2S
      EX(2,1)=2.768496241D+01 
      C(2,1)=-4.151277819D-03
      EX(2,2)=5.077140627D+00 
      C(2,2)=-2.067024148D-02
      EX(2,3)=1.426786050D+00 
      C(2,3)=-5.150303337D-02
      EX(2,4)=2.040335729D-01 
      C(2,4)=3.346271174D-01
      EX(2,5)=9.260298399D-02 
      C(2,5)=5.621061301D-01
      EX(2,6)=4.416183978D-02 
      C(2,6)=1.712994697D-01
C  2P
      EX(3,1)=5.868285913D+00 
      C(3,1)=7.924233646D-03
      EX(3,2)=1.530329631D+00 
      C(3,2)=5.144104825D-02
      EX(3,3)=5.475665231D-01 
      C(3,3)=1.898400060D-01
      EX(3,4)=2.288932733D-01 
      C(3,4)=4.049863191D-01
      EX(3,5)=1.046655969D-01 
      C(3,5)=4.012362861D-01
      EX(3,6)=4.948220127D-02 
      C(3,6)=1.051855189D-01
C 3S
      EX(4,1) = 3.273031938D+00
      C(4,1) = -6.775596947D-03
      EX(4,2) = 9.200611311D-01
      C(4,2) = -5.639325779D-02
      EX(4,3) = 3.593349765D-01
      C(4,3) = -1.587856086D-01
      EX(4,4) = 8.636686991D-02
      C(4,4) = 5.534527651D-01
      EX(4,5) = 4.797373812D-02
      C(4,5) = 5.015351020D-01
      EX(4,6) = 2.724741144D-02
      C(4,6) = 7.223633674D-02
C 3P
      EX(5,1) = 5.077973607D+00
      C(5,1) = -3.329929840D-03
      EX(5,2) = 1.340786940D+00
      C(5,2) = -1.419488340D-02
      EX(5,3) = 2.248434849D-01
      C(5,3) = 1.639395770D-01
      EX(5,4) = 1.131741848D-01
      C(5,4) = 4.485358256D-01
      EX(5,5) = 6.076408893D-02
      C(5,5) = 3.908813050D-01
      EX(5,6) = 3.315424265D-02
      C(5,6) = 7.411456232D-02
C 4S
      EX(6,1) = 3.232838646D+00
      C(6,1) = 1.374817488D-03
      EX(6,2) = 3.605788802D-01
      C(6,2) = -8.666390043D-02
      EX(6,3) = 1.717905487D-01
      C(6,3) = -3.130627309D-01
      EX(6,4) = 5.277666487D-02
      C(6,4) = 7.812787397D-01
      EX(6,5) = 3.163400284D-02
      C(6,5) = 4.389247988D-01
      EX(6,6) = 1.874093091D-02
      C(6,6) = 2.487178756D-02
C 4P
      EX(7,1) = 2.389722618D+00
      C(7,1) = -1.665913575D-03
      EX(7,2) = 7.960947826D-01
      C(7,2) = -1.657464971D-02
      EX(7,3) = 3.415541380D-01
      C(7,3) = -5.958513378D-02
      EX(7,4) = 8.847434525D-02
      C(7,4) = 4.053115554D-01
      EX(7,5) = 4.958248334D-02
      C(7,5) = 5.433958189D-01
      EX(7,6) = 2.816929784D-02
      C(7,6) = 1.204970491D-01
C 3D
      EX(8,1) = 2.488296923D+00
      C(8,1) = 7.283828112D-03
      EX(8,2) = 7.981487853D-01
      C(8,2) = 5.386799363D-02
      EX(8,3) = 3.311327490D-01
      C(8,3) = 2.072139149D-01
      EX(8,4) = 1.559114463D-01
      C(8,4) = 4.266269092D-01
      EX(8,5) = 7.877734732D-02
      C(8,5) = 3.843100204D-01
      EX(8,6) = 4.058484363D-02
      C(8,6) = 8.902827546D-02
C
      RETURN
      END
