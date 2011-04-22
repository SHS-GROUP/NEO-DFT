C 11 AUG 10 - MK  - NEW CODE FOR RESOLUTION OF THE IDENTITY 3-CENTER INTS
C
C*MODULE RIINT   *DECK FORM2C1OV
      SUBROUTINE FORM2C1OV(S,DI,DJ,SS,LLI,LLJ,LI,LJ,NI,NJ)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION S(3,0:LI,0:LJ),DI(NI),DJ(NJ),SS(*)
      INTEGER LLI(3,NI),LLJ(3,NJ)
C
      IJ = 0
      DO J = 1, NJ
         NJX = LLJ(1,J)
         NJY = LLJ(2,J)
         NJZ = LLJ(3,J)
         DO I = 1, NI
            NIX = LLI(1,I)
            NIY = LLI(2,I)
            NIZ = LLI(3,I)
            VAL = S(1,NIX,NJX)*S(2,NIY,NJY)*S(3,NIZ,NJZ)
C      WRITE(6,*) 'VAL=',VAL,S(1,NIX,NJX),S(2,NIY,NJY),S(3,NIZ,NJZ)
            IJ = IJ+1
            SS(IJ) = SS(IJ)+DI(I)*DJ(J)*VAL
         END DO
      END DO
C
      RETURN
      END
C*MODULE RIINT   *DECK FORM2C2E
      SUBROUTINE FORM2C2E(G,DI,DJ,GG,LLI,LLJ,LI,LJ,NI,NJ,NROOTS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION G(3,NROOTS,0:LI,0:LJ),DI(NI),DJ(NJ),GG(*)
      INTEGER LLI(3,NI),LLJ(3,NJ)
C
      PARAMETER (ZERO=0.0D+00,PI252=34.986836655250D+00)
C
      IJ = 0
      DO J = 1, NJ
         NJX = LLJ(1,J)
         NJY = LLJ(2,J)
         NJZ = LLJ(3,J)
         DO I = 1, NI
            NIX = LLI(1,I)
            NIY = LLI(2,I)
            NIZ = LLI(3,I)
            VAL = ZERO
            DO IR = 1, NROOTS
               VAL = VAL+G(1,IR,NIX,NJX)*G(2,IR,NIY,NJY)*G(3,IR,NIZ,NJZ)
            END DO
            IJ = IJ+1
            GG(IJ) = GG(IJ)+PI252*DI(I)*DJ(J)*VAL
         END DO
      END DO
C
      RETURN
      END
C*MODULE RIINT   *DECK FORM3C2EP
      SUBROUTINE FORM3C2EP(G3,DB,DK,GG,LLI,LLJ,LLK,
     *     LI,LJ,LK,NI,NJ,NK,NROOTS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION G3(3,NROOTS,0:LI,0:LJ,0:LK),DB(NI*NJ),DK(NK),
     *     GG(*)
      INTEGER LLI(3,NI),LLJ(3,NJ),LLK(3,NK)
C
      PARAMETER (ZERO=0.0D+00, PI252=34.986836655250D+00)
C
      IJ = 0
      IJK = 0
      DO I = 1, NI
         NIX = LLI(1,I)
         NIY = LLI(2,I)
         NIZ = LLI(3,I)
         DO J = 1, NJ
            NJX = LLJ(1,J)
            NJY = LLJ(2,J)
            NJZ = LLJ(3,J)
            IJ = IJ+1
            DO K = 1, NK
               NKX = LLK(1,K)
               NKY = LLK(2,K)
               NKZ = LLK(3,K)
               VAL = ZERO
               DO IR = 1, NROOTS
                  VAL = VAL+G3(1,IR,NIX,NJX,NKX)*G3(2,IR,NIY,NJY,NKY)
     *                 *G3(3,IR,NIZ,NJZ,NKZ)
               END DO
               IJK = IJK+1
               GG(IJK) = GG(IJK)+PI252*DB(IJ)*DK(K)*VAL
            END DO
         END DO
      END DO
C
      RETURN
      END
C*MODULE RIINT   *DECK GEN2C1OV
      SUBROUTINE GEN2C1OV(CC,EXI,CIS,CIP,CID,CIF,CIG,CIH,CII,
     *     EXJ,CJS,CJP,CJD,CJF,CJG,CJH,CJI,SS,
     *     LI,LJ,NGI,NGJ,MINI,MAXI,MINJ,MAXJ,NORMI,NORMJ)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL NORMI,NORMJ
C
      DOUBLE PRECISION CC(3,2),SS(*)
      DOUBLE PRECISION EXI(NGI),CIS(NGI),CIP(NGI),CID(NGI),CIF(NGI),
     *     CIG(NGI),CIH(NGI),CII(NGI),EXJ(NGJ),CJS(NGJ),CJP(NGJ),
     *     CJD(NGJ),CJF(NGJ),CJG(NGJ),CJH(NGJ),CJI(NGJ)
      DOUBLE PRECISION DI(28),DJ(28),S1(39),S(147)
      DOUBLE PRECISION A1(3),A2(3)
      INTEGER LLI(3,28),LLJ(3,28),LX(84),LY(84),LZ(84)
C
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00, SQRT9=3.0D+00,
     *           SQRT11=3.3166247903553998D+00,
     *           ZERO=0.0D+00, RLN10=2.30258D+00)
      PARAMETER (ITOLOV=20)
C
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,
     *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,
     *            4,   0,   0,   3,   3,   1,   0,   1,   0,   2,
     *            2,   0,   2,   1,   1,
     *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,
     *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,
     *            1,
     *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,
     *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,
     *            0,   3,   3,   2,   1,   2,   1,   2/
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,
     *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,
     *            0,   4,   0,   1,   0,   3,   3,   0,   1,   2,
     *            0,   2,   1,   2,   1,
     *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,
     *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,
     *            2,
     *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,
     *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,
     *            3,   2,   1,   3,   3,   1,   2,   2/
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,
     *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,
     *            0,   0,   4,   0,   1,   0,   1,   3,   3,   0,
     *            2,   2,   1,   1,   2,
     *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,
     *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,
     *            2,
     *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,
     *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,
     *            3,   1,   2,   1,   2,   3,   3,   2/
C
C     COMPUTES 2C-1E OVERLAP INTEGRALS (IJ)
C
C     EXI, EXJ: EXPONENT OF GAUSSIAN FUNCTIONS
C     CSI, CSJ: COEFFICINT OF GAUSSIAN FUNCTIONS
C     LI, LJ:   ANGLAR MOMENTUM
C     NGI,NGJ:  # OF PRIMITIVE GAUSSIAN FUNCTIONS
C     MINI, MINJ:
C     MAXI, MAXJ:
C     ITOL: THRESHOLD OF PRIMITIVE GAUSSIAN FUNCTIONS
C
C      TOL = RLN10*ITOL
      TOL = RLN10*ITOLOV
C
      LTOT = LI+LJ
      NI = MAXI-MINI+1
      NJ = MAXJ-MINJ+1
      A1(1) = CC(1,1)-CC(1,2)
      A1(2) = CC(2,1)-CC(2,2)
      A1(3) = CC(3,1)-CC(3,2)
      RIJ2 = A1(1)*A1(1)+A1(2)*A1(2)+A1(3)*A1(3)
C
C      WRITE(6,*) 'LI=',LI,'LJ=',LJ,'LTOT=',LTOT
C      WRITE(6,*) 'NGI=',NGI,'NGJ=',NGJ
C      WRITE(6,*) 'NI=',NI,'NJ=',NJ
C
      DO I = 1, NI*NJ
         SS(I) = ZERO
      END DO
C
      ICNT = 0
      DO I = MINI, MAXI
         ICNT = ICNT+1
         LLI(1,ICNT) = LX(I)
         LLI(2,ICNT) = LY(I)
         LLI(3,ICNT) = LZ(I)
C         WRITE(6,*) ICNT,(LLI(IA,ICNT),IA=1,3)
      END DO
      JCNT = 0
      DO J = MINJ, MAXJ
         JCNT = JCNT+1
         LLJ(1,JCNT) = LX(J)
         LLJ(2,JCNT) = LY(J)
         LLJ(3,JCNT) = LZ(J)
C         WRITE(6,*) JCNT,(LLJ(IA,JCNT),IA=1,3)
      END DO
C
C     I PRIMITIVE
C
      DO IG = 1, NGI
         AI  = EXI(IG)
C
C        ----- DENSITY FACTOR OF I PRIMITIVE
C
         ICNT = 0
         DO I = MINI, MAXI
            IF(I .EQ. 1) THEN
               DUM1 = CIS(IG)
            ELSE IF(I .LE.  4) THEN
               DUM1 = CIP(IG)
            ELSE IF(I .LE. 10) THEN
               DUM1 = CID(IG)
               IF(I .GE. 8 .AND. NORMI) DUM1 = DUM1*SQRT3
            ELSE IF(I .LE. 20) THEN
               DUM1 = CIF(IG)
               IF(I .GE. 14 .AND. NORMI) THEN
                  IF(I .LE. 19) THEN
                     DUM1 = DUM1*SQRT5
                  ELSE
                     DUM1 = DUM1*SQRT3
                  END IF
               END IF
            ELSE IF(I .LE. 35) THEN
               DUM1 = CIG(IG)
               IF(I .GE. 24 .AND. NORMI) THEN
                  IF(I .LE. 29) THEN
                     DUM1 = DUM1*SQRT7
                  ELSE IF(I .LE. 32) THEN
                     DUM1 = DUM1*SQRT5/SQRT3
                  ELSE
                     DUM1 = DUM1*SQRT3
                  END IF
               END IF
            ELSE IF(I .LE. 56) THEN
               DUM1 = CIH(IG)
               IF(I .GE. 39 .AND. NORMI) THEN
                  IF(I .LE. 44) THEN
                     DUM1 = DUM1*SQRT9
                  ELSE IF(I .LE. 50) THEN
                     DUM1 = DUM1*SQRT7/SQRT3
                  ELSE IF(I .LE. 53) THEN
                     DUM1 = DUM1*SQRT3
                  ELSE
                     DUM1 = DUM1*SQRT5/SQRT3
                  ENDIF
               ENDIF
            ELSE
               DUM1 = CII(IG)
               IF(I .GE. 60 .AND. NORMI) THEN
                  IF(I .LE. 65) THEN
                     DUM1 = DUM1*SQRT11
                  ELSE IF(I .LE. 71) THEN
                     DUM1 = DUM1*SQRT3
                  ELSE IF(I .LE. 74) THEN
                     DUM1 = DUM1*SQRT3
                  ELSE IF(I .LE. 77) THEN
                     DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)
                  ELSE IF(I .LE. 83) THEN
                     DUM1 = DUM1*SQRT5
                  ELSE
                     DUM1 = DUM1*SQRT5/SQRT3
                  ENDIF
               END IF
            ENDIF
            ICNT = ICNT+1
            DI(ICNT) = DUM1
         END DO
C
C     J PRIMITIVE
C
         DO JG = 1, NGJ
            AJ  = EXJ(JG)
            AMU = AI*AJ/(AI+AJ)
            DUM = AMU*RIJ2
            IF(DUM .GT. TOL) GO TO 500
            EXPFAC = EXP(-DUM)
C            WRITE(6,*) 'EXPFAC=',EXPFAC
C
C           ----- DENSITY FACTOR OF J PRIMITIVE
C
            JCNT = 0
            DO J = MINJ, MAXJ
               IF(J .EQ. 1) THEN
                  DUM2 = CJS(JG)
               ELSE IF(J .LE.  4) THEN
                  DUM2 = CJP(JG)
               ELSE IF(J .LE. 10) THEN
                  DUM2 = CJD(JG)
                  IF(J .GE. 8 .AND. NORMJ) DUM2 = DUM2*SQRT3
               ELSE IF(J .LE. 20) THEN
                  DUM2 = CJF(JG)
                  IF(J .GE. 14 .AND. NORMJ) THEN
                     IF(J .LE. 19) THEN
                        DUM2 = DUM2*SQRT5
                     ELSE
                        DUM2 = DUM2*SQRT3
                     END IF
                  END IF
               ELSE IF(J .LE. 35) THEN
                  DUM2 = CJG(JG)
                  IF(J .GE. 24 .AND. NORMJ) THEN
                     IF(J .LE. 29) THEN
                        DUM2 = DUM2*SQRT7
                     ELSE IF(J .LE. 32) THEN
                        DUM2 = DUM2*SQRT5/SQRT3
                     ELSE
                        DUM2 = DUM2*SQRT3
                     END IF
                  END IF
               ELSE IF(J .LE. 56) THEN
                  DUM2 = CJH(JG)
                  IF(J .GE. 39 .AND. NORMJ) THEN
                     IF(J .LE. 44) THEN
                        DUM2 = DUM2*SQRT9
                     ELSE IF(J .LE. 50) THEN
                        DUM2 = DUM2*SQRT7/SQRT3
                     ELSE IF(J .LE. 53) THEN
                        DUM2 = DUM2*SQRT3
                     ELSE
                        DUM2 = DUM2*SQRT5/SQRT3
                     ENDIF
                  END IF
               ELSE
                  DUM2 = CJI(JG)
                  IF(J .GE. 60 .AND. NORMJ) THEN
                     IF(J .LE. 65) THEN
                        DUM2 = DUM2*SQRT11
                     ELSE IF(J .LE. 71) THEN
                        DUM2 = DUM2*SQRT3
                     ELSE IF(J .LE. 74) THEN
                        DUM2 = DUM2*SQRT3
                     ELSE IF(J .LE. 77) THEN
                        DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)
                     ELSE IF(J .LE. 83) THEN
                        DUM2 = DUM2*SQRT5
                     ELSE
                        DUM2 = DUM2*SQRT5/SQRT3
                     ENDIF
                  END IF
               END IF
               JCNT = JCNT+1
               DJ(JCNT) = DUM2
            END DO
C
C           ----- EVALUATE 1D INTEGRAL
C
            CALL XYZ2C1OV(S1,S,A1,A2,AI,AJ,EXPFAC,LI,LJ,LTOT)
C
C           ----- FORM (IJ) OVERLAP
C
            CALL FORM2C1OV(S,DI,DJ,SS,LLI,LLJ,LI,LJ,NI,NJ)
C
 500        CONTINUE
         END DO
      END DO
C
      RETURN
      END
C*MODULE RIINT   *DECK GEN2C2ERI
      SUBROUTINE GEN2C2ERI(CC,EXI,CIS,CIP,CID,CIF,CIG,CIH,CII,
     *     EXJ,CJS,CJP,CJD,CJF,CJG,CJH,CJI,GG,
     *     LI,LJ,NGI,NGJ,MINI,MAXI,MINJ,MAXJ,NORMI,NORMJ)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL NORMI,NORMJ
C
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
C
      DOUBLE PRECISION CC(3,2),GG(*)
      DOUBLE PRECISION EXI(NGI),CIS(NGI),CIP(NGI),CID(NGI),CIF(NGI),
     *     CIG(NGI),CIH(NGI),CII(NGI),EXJ(NGJ),CJS(NGJ),CJP(NGJ),
     *     CJD(NGJ),CJF(NGJ),CJG(NGJ),CJH(NGJ),CJI(NGJ)
      DOUBLE PRECISION DI(28),DJ(28),G(1911)
      DOUBLE PRECISION A1(3,13),B1(13),C1(13),A2(3,13),B2(13),F00(13)
      INTEGER LLI(3,28),LLJ(3,28),LX(84),LY(84),LZ(84)
C
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00,
     *           SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,
     *           ZERO=0.0D+00, ONE=1.0D+00, THREE=3.0D+00)
C
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,
     *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,
     *            4,   0,   0,   3,   3,   1,   0,   1,   0,   2,
     *            2,   0,   2,   1,   1,
     *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,
     *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,
     *            1,
     *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,
     *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,
     *            0,   3,   3,   2,   1,   2,   1,   2/
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,
     *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,
     *            0,   4,   0,   1,   0,   3,   3,   0,   1,   2,
     *            0,   2,   1,   2,   1,
     *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,
     *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,
     *            2,
     *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,
     *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,
     *            3,   2,   1,   3,   3,   1,   2,   2/
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,
     *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,
     *            0,   0,   4,   0,   1,   0,   1,   3,   3,   0,
     *            2,   2,   1,   1,   2,
     *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,
     *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,
     *            2,
     *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,
     *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,
     *            3,   1,   2,   1,   2,   3,   3,   2/
C
C     COMPUTES 2C-2E ERI (I|J) BY GAUSS-RYS METHOD
C
C     EXI, EXJ: EXPONENT OF GAUSSIAN FUNCTIONS
C     CSI, CSJ: COEFFICINT OF GAUSSIAN FUNCTIONS
C     LI, LJ:   ANGLAR MOMENTUM
C     NGI,NGJ:  # OF PRIMITIVE GAUSSIAN FUNCTIONS
C     MINI, MINJ:
C     MAXI, MAXJ:
C     ITOL: THRESHOLD OF PRIMITIVE GAUSSIAN FUNCTIONS
C
      LTOT = LI+LJ
      NROOTS = LTOT/2+1
      NI = MAXI-MINI+1
      NJ = MAXJ-MINJ+1
      RXIJ = CC(1,1)-CC(1,2)
      RYIJ = CC(2,1)-CC(2,2)
      RZIJ = CC(3,1)-CC(3,2)
      RIJ2 = RXIJ**2+RYIJ**2+RZIJ**2
C
C      WRITE(6,*) 'LI=',LI,'LJ=',LJ
C      WRITE(6,*) 'NGI=',NGI,'NGJ=',NGJ
C      WRITE(6,*) 'NI=',NI,'NJ=',NJ
C      WRITE(6,*) 'LTOT=',LTOT,'NROOTS=',NROOTS
C
      DO I = 1, NI*NJ
         GG(I) = ZERO
      END DO
C
      ICNT = 0
      DO I = MINI, MAXI
         ICNT = ICNT+1
         LLI(1,ICNT) = LX(I)
         LLI(2,ICNT) = LY(I)
         LLI(3,ICNT) = LZ(I)
C         WRITE(6,*) ICNT,(LLI(IA,ICNT),IA=1,3)
      END DO
      JCNT = 0
      DO J = MINJ, MAXJ
         JCNT = JCNT+1
         LLJ(1,JCNT) = LX(J)
         LLJ(2,JCNT) = LY(J)
         LLJ(3,JCNT) = LZ(J)
C         WRITE(6,*) JCNT,(LLJ(IA,JCNT),IA=1,3)
      END DO
C
C     I PRIMITIVE
C
      DO IG = 1, NGI
         AI  = EXI(IG)
         AIINV = ONE/AI
C
C        ----- DENSITY FACTOR OF I PRIMITIVE
C
         ICNT = 0
         DO I = MINI, MAXI
            IF(I .EQ. 1) THEN
               DUM1 = CIS(IG)
            ELSE IF(I .LE.  4) THEN
               DUM1 = CIP(IG)
            ELSE IF(I .LE. 10) THEN
               DUM1 = CID(IG)
               IF(I .GE. 8 .AND. NORMI) DUM1 = DUM1*SQRT3
            ELSE IF(I .LE. 20) THEN
               DUM1 = CIF(IG)
               IF(I .GE. 14 .AND. NORMI) THEN
                  IF(I .LE. 19) THEN
                     DUM1 = DUM1*SQRT5
                  ELSE
C                     DUM1 = DUM1*SQRT3
                     DUM1 = DUM1*SQRT3*SQRT5
                  END IF
               END IF
            ELSE IF(I .LE. 35) THEN
               DUM1 = CIG(IG)
               IF(I .GE. 24 .AND. NORMI) THEN
                  IF(I .LE. 29) THEN
                     DUM1 = DUM1*SQRT7
                  ELSE IF(I .LE. 32) THEN
C                     DUM1 = DUM1*SQRT5/SQRT3
                     DUM1 = DUM1*SQRT7*SQRT5/SQRT3
                  ELSE
C                     DUM1 = DUM1*SQRT3
                     DUM1 = DUM1*SQRT7*SQRT5
                  END IF
               END IF
            ELSE IF(I .LE. 56) THEN
               DUM1 = CIH(IG)
               IF(I .GE. 39 .AND. NORMI) THEN
                  IF(I .LE. 44) THEN
                     DUM1 = DUM1*SQRT9
                  ELSE IF(I .LE. 50) THEN
C                     DUM1 = DUM1*SQRT7/SQRT3
                     DUM1 = DUM1*SQRT7*SQRT3
                  ELSE IF(I .LE. 53) THEN
C                     DUM1 = DUM1*SQRT3
                     DUM1 = DUM1*SQRT9*SQRT7
                  ELSE
C                     DUM1 = DUM1*SQRT5/SQRT3
                     DUM1 = DUM1*SQRT7*SQRT5*SQRT3
                  ENDIF
               ENDIF
            ELSE
               DUM1 = CII(IG)
               IF(I .GE. 60 .AND. NORMI) THEN
                  IF(I .LE. 65) THEN
                     DUM1 = DUM1*SQRT11
                  ELSE IF(I .LE. 71) THEN
C                     DUM1 = DUM1*SQRT3
                     DUM1 = DUM1*SQRT11*SQRT3
                  ELSE IF(I .LE. 74) THEN
C                     DUM1 = DUM1*SQRT3
                     DUM1 = DUM1*SQRT11*THREE
                  ELSE IF(I .LE. 77) THEN
C                     DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT3/SQRT5
                  ELSE IF(I .LE. 83) THEN
C                     DUM1 = DUM1*SQRT5
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT3
                  ELSE
C                     DUM1 = DUM1*SQRT5/SQRT3
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT5
                  ENDIF
               END IF
            ENDIF
            ICNT = ICNT+1
            DI(ICNT) = DUM1*AIINV
         END DO
C
C     J PRIMITIVE
C
         DO JG = 1, NGJ
            AJ  = EXJ(JG)
            AJINV = ONE/AJ
C
C           ----- DENSITY FACTOR OF J PRIMITIVE
C
            JCNT = 0
            DO J = MINJ, MAXJ
               IF(J .EQ. 1) THEN
                  DUM2 = CJS(JG)
               ELSE IF(J .LE.  4) THEN
                  DUM2 = CJP(JG)
               ELSE IF(J .LE. 10) THEN
                  DUM2 = CJD(JG)
                  IF(J .GE. 8 .AND. NORMJ) DUM2 = DUM2*SQRT3
               ELSE IF(J .LE. 20) THEN
                  DUM2 = CJF(JG)
                  IF(J .GE. 14 .AND. NORMJ) THEN
                     IF(J .LE. 19) THEN
                        DUM2 = DUM2*SQRT5
                     ELSE
C                        DUM2 = DUM2*SQRT3
                        DUM2 = DUM2*SQRT3*SQRT5
                     END IF
                  END IF
               ELSE IF(J .LE. 35) THEN
                  DUM2 = CJG(JG)
                  IF(J .GE. 24 .AND. NORMJ) THEN
                     IF(J .LE. 29) THEN
                        DUM2 = DUM2*SQRT7
                     ELSE IF(J .LE. 32) THEN
C                        DUM2 = DUM2*SQRT5/SQRT3
                        DUM2 = DUM2*SQRT7*SQRT5/SQRT3
                     ELSE
C                        DUM2 = DUM2*SQRT3
                        DUM2 = DUM2*SQRT7*SQRT5
                     END IF
                  END IF
               ELSE IF(J .LE. 56) THEN
                  DUM2 = CJH(JG)
                  IF(J .GE. 39 .AND. NORMJ) THEN
                     IF(J .LE. 44) THEN
                        DUM2 = DUM2*SQRT9
                     ELSE IF(J .LE. 50) THEN
C                        DUM2 = DUM2*SQRT7/SQRT3
                        DUM2 = DUM2*SQRT7*SQRT3
                     ELSE IF(J .LE. 53) THEN
C                        DUM2 = DUM2*SQRT3
                        DUM2 = DUM2*SQRT9*SQRT7
                     ELSE
C                        DUM2 = DUM2*SQRT5/SQRT3
                        DUM2 = DUM2*SQRT7*SQRT5*SQRT3
                     ENDIF
                  END IF
               ELSE
                  DUM2 = CJI(JG)
                  IF(J .GE. 60 .AND. NORMJ) THEN
                     IF(J .LE. 65) THEN
                        DUM2 = DUM2*SQRT11
                     ELSE IF(J .LE. 71) THEN
C                        DUM2 = DUM2*SQRT3
                        DUM2 = DUM2*SQRT11*SQRT3
                     ELSE IF(J .LE. 74) THEN
C                        DUM2 = DUM2*SQRT3
                        DUM2 = DUM2*SQRT11*THREE
                     ELSE IF(J .LE. 77) THEN
C                        DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)
                        DUM2 = DUM2*SQRT11*SQRT7*SQRT5/SQRT3
                     ELSE IF(J .LE. 83) THEN
C                        DUM2 = DUM2*SQRT5
                        DUM2 = DUM2*SQRT11*SQRT7*SQRT3
                     ELSE
C                        DUM2 = DUM2*SQRT5/SQRT3
                        DUM2 = DUM2*SQRT11*SQRT7*SQRT5
                     ENDIF
                  END IF
               END IF
               JCNT = JCNT+1
               DJ(JCNT) = DUM2*AJINV
            END DO
C
C           ----- GET ROOTS AND WEIGHTS FOR QUADRATURE
C
            RHO = AI*AJ/(AI+AJ)
            XX = RHO*RIJ2
            IF (NROOTS .LE. 3) CALL RT123
            IF (NROOTS .EQ. 4) CALL ROOT4
            IF (NROOTS .EQ. 5) CALL ROOT5
            IF (NROOTS .GE. 6) CALL ROOT6
C
C           ----- EVALUATE 1D INTEGRAL BY R.R.
C
            CALL XYZ2C2E(CC,U,W,G,A1,B1,C1,A2,B2,F00,AI,AJ,LI,LJ,NROOTS)
C
C           ----- FORM (I/J) ERIS
C
            CALL FORM2C2E(G,DI,DJ,GG,LLI,LLJ,LI,LJ,NI,NJ,NROOTS)
C
         END DO
      END DO
C
      RETURN
      END
C*MODULE RIINT   *DECK GEN3C2ERI
      SUBROUTINE GEN3C2ERI(CC,EXI,CIS,CIP,CID,CIF,CIG,CIH,CII,
     *     EXJ,CJS,CJP,CJD,CJF,CJG,CJH,CJI,EXK,CKS,CKP,CKD,CKF,CKG,CKH,
     *     CKI,GG,
     *     LI,LJ,LK,NGI,NGJ,NGK,MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     *     NORMI,NORMJ,NORMK,AKINVV,DK)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL NORMI,NORMJ,NORMK
C
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
C
      PARAMETER (MXLTOT=6)
      DOUBLE PRECISION CC(3,3),GG(*)
      DOUBLE PRECISION EXI(NGI),CIS(NGI),CIP(NGI),CID(NGI),CIF(NGI),
     *     CIG(NGI),CIH(NGI),CII(NGI),EXJ(NGJ),CJS(NGJ),CJP(NGJ),
     *     CJD(NGJ),CJF(NGJ),CJG(NGJ),CJH(NGJ),CJI(NGJ),
     *     EXK(NGK),CKS(NGK),CKP(NGK),CKD(NGK),CKF(NGK),CKG(NGK),
     *     CKH(NGK),CKI(NGK),AKINVV(NGK),DK(28,NGK)
      DOUBLE PRECISION G2(2730),G3(10290),CAB(3,0:MXLTOT),
     *     CNK(28),DB(784),P(3),PAV(3),CPV(3)
      DOUBLE PRECISION A1(3,13),B1(13),C1(13),
     *                 A2(3,13),B2(13),F00(13)
      INTEGER LLI(3,28),LLJ(3,28),LLK(3,28),LX(84),LY(84),LZ(84)
C
      PARAMETER (SQRT3=1.73205080756888D+00, SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00,
     *           SQRT9=3.0D+00,SQRT11=3.3166247903553998D+00,
     *           ZERO=0.0D+00, ONE=1.0D+00, THREE=3.0D+00)
      PARAMETER (ITOLRI=10, RLN10=2.30258D+00)
C
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,
     *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,
     *            4,   0,   0,   3,   3,   1,   0,   1,   0,   2,
     *            2,   0,   2,   1,   1,
     *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,
     *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,
     *            1,
     *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,
     *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,
     *            0,   3,   3,   2,   1,   2,   1,   2/
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,
     *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,
     *            0,   4,   0,   1,   0,   3,   3,   0,   1,   2,
     *            0,   2,   1,   2,   1,
     *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,
     *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,
     *            2,
     *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,
     *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,
     *            3,   2,   1,   3,   3,   1,   2,   2/
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,
     *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,
     *            0,   0,   4,   0,   1,   0,   1,   3,   3,   0,
     *            2,   2,   1,   1,   2,
     *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,
     *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,
     *            2,
     *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,
     *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,
     *            3,   1,   2,   1,   2,   3,   3,   2/
C
C      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,
C     *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,
C     *            4,   3,   3,   2,   2,   2,   1,   1,   1,   1,
C     *            0,   0,   0,   0,   0,
C     *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,
C     *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,
C     *            1,
C     *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,
C     *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,
C     *            0,   3,   3,   2,   1,   2,   1,   2/
C      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,
C     *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,
C     *            0,   1,   0,   2,   1,   0,   3,   2,   1,   0,
C     *            4,   3,   2,   1,   0,
C     *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,
C     *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,
C     *            2,
C     *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,
C     *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,
C     *            3,   2,   1,   3,   3,   1,   2,   2/
C      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,
C     *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,
C     *            0,   0,   1,   0,   1,   2,   0,   1,   2,   3,
C     *            0,   1,   2,   3,   4,
C     *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,
C     *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,
C     *            2,
C     *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,
C     *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,
C     *            3,   1,   2,   1,   2,   3,   3,   2/
C
      TOL = ITOLRI*RLN10
C
C     COMPUTES 3C-2E ERI (IJ|K) BY GAUSS-RYS METHOD
C
C     EXI, EXJ, EXK: EXPONENT OF GAUSSIAN FUNCTIONS
C     CIS, CJS, CKS: COEFFICINT OF GAUSSIAN FUNCTIONS
C     LI, LJ, LK:    ANGLAR MOMENTUM
C     NGI,NGJ,NGK:   NUMBER OF PRIMITIVE GAUSSIAN FUNCTIONS
C     MINI, MINJ, MINK:
C     MAXI, MAXJ, MAXK:
C     ITOL: THRESHOLD OF PRIMITIVE GAUSSIAN FUNCTIONS
C
      LB = LI+LJ
      LTOT = LB+LK
      NROOTS = LTOT/2+1
      NI = MAXI-MINI+1
      NJ = MAXJ-MINJ+1
      NK = MAXK-MINK+1
C
C      WRITE(6,*) 'LI=',LI,'LJ=',LJ,'LK=',LK
C      WRITE(6,*) 'NGI=',NGI,'NGJ=',NGJ,'NGK=',NGK
C      WRITE(6,*) 'NI=',NI,'NJ=',NJ,'NK=',NK
C
      RAB2 = ZERO
      DO IA = 1, 3
         VAL = CC(IA,1)-CC(IA,2)
         RAB2 = RAB2+VAL*VAL
         DO J = 0, LJ
            CAB(IA,J) = VAL**J
         END DO
      END DO
C      WRITE(6,*) 'CAB'
C      DO J = 0,L
C         WRITE(6,'(3E15.8)') (CAB(IA,J),IA=1,3)
C      END DO
C
      DO J = 0, LJ
         JJ = J*(J+1)/2+1
         CNK(JJ) = ONE
         JJ = (J+1)*(J+2)/2
         CNK(JJ) = ONE
         IF(J .GE. 2) THEN
            DO K = 1, J-1
               JJ = J*(J+1)/2+K+1
               JA = (J-1)*J/2+K
               JB = (J-1)*J/2+K+1
               CNK(JJ) = CNK(JA)+CNK(JB)
            END DO
         END IF
      END DO
C      DO J = 0, LJ
C         DO K = 0, J
C            JJ = J*(J+1)/2+K+1
C            WRITE(6,*) 'CNK=',J,K,JJ,CNK(JJ)
C         END DO
C      END DO
C
      DO IJK = 1, NI*NJ*NK
         GG(IJK) = ZERO
      END DO
C
      ICNT = 0
      DO I = MINI, MAXI
         ICNT = ICNT+1
         LLI(1,ICNT) = LX(I)
         LLI(2,ICNT) = LY(I)
         LLI(3,ICNT) = LZ(I)
      END DO
      JCNT = 0
      DO J = MINJ, MAXJ
         JCNT = JCNT+1
         LLJ(1,JCNT) = LX(J)
         LLJ(2,JCNT) = LY(J)
         LLJ(3,JCNT) = LZ(J)
      END DO
      KCNT = 0
      DO K = MINK, MAXK
         KCNT = KCNT+1
         LLK(1,KCNT) = LX(K)
         LLK(2,KCNT) = LY(K)
         LLK(3,KCNT) = LZ(K)
      END DO
C
C     ----- DENSITY FACTOR OF K PRIMITIVE
C
      DO KG = 1, NGK
         AKINVV(KG) = ONE/EXK(KG)
         KCNT = 0
         DO K = MINK, MAXK
            IF(K .EQ. 1) THEN
               DUM4 = CKS(KG)
            ELSE IF(K .LE.  4) THEN
               DUM4 = CKP(KG)
            ELSE IF(K .LE. 10) THEN
               DUM4 = CKD(KG)
               IF(K .GE. 8 .AND. NORMK) DUM4 = DUM4*SQRT3
            ELSE IF(K .LE. 20) THEN
               DUM4 = CKF(KG)
               IF(K .GE. 14 .AND. NORMK) THEN
                  IF(K .LE. 19) THEN
                     DUM4 = DUM4*SQRT5
                  ELSE
C                     DUM4 = DUM4*SQRT3
                     DUM4 = DUM4*SQRT3*SQRT5
                  END IF
               END IF
            ELSE IF(K .LE. 35) THEN
               DUM4 = CKG(KG)
               IF(K .GE. 24 .AND. NORMK) THEN
                  IF(K .LE. 29) THEN
                     DUM4 = DUM4*SQRT7
                  ELSE IF(K .LE. 32) THEN
C                     DUM4 = DUM4*SQRT5/SQRT3
                     DUM4 = DUM4*SQRT7*SQRT5/SQRT3
                  ELSE
C                     DUM4 = DUM4*SQRT3
                     DUM4 = DUM4*SQRT7*SQRT5
                  END IF
               END IF
            ELSE IF(K .LE. 56) THEN
               DUM4 = CKH(KG)
               IF(K .GE. 39 .AND. NORMK) THEN
                  IF(K .LE. 44) THEN
                     DUM4 = DUM4*SQRT9
                  ELSE IF(K .LE. 50) THEN
C                     DUM4 = DUM4*SQRT7/SQRT3
                     DUM4 = DUM4*SQRT7*SQRT3
                  ELSE IF(K .LE. 53) THEN
C                     DUM4 = DUM4*SQRT3
                     DUM4 = DUM4*SQRT9*SQRT7
                  ELSE
C                     DUM4 = DUM4*SQRT5/SQRT3
                     DUM4 = DUM4*SQRT7*SQRT5*SQRT3
                  ENDIF
               ENDIF
            ELSE
               DUM4 = CKI(KG)
               IF(K .GE. 60 .AND. NORMK) THEN
                  IF(K .LE. 65) THEN
                     DUM4 = DUM4*SQRT11
                  ELSE IF(K .LE. 71) THEN
C                     DUM4 = DUM4*SQRT3
                     DUM4 = DUM4*SQRT11*SQRT3
                  ELSE IF(K .LE. 74) THEN
C                     DUM4 = DUM4*SQRT3
                     DUM4 = DUM4*SQRT11*THREE
                  ELSE IF(K .LE. 77) THEN
C                     DUM4 = DUM4*SQRT7/(SQRT5*SQRT3)
                     DUM4 = DUM4*SQRT11*SQRT7*SQRT3/SQRT5
                  ELSE IF(K .LE. 83) THEN
C                     DUM4 = DUM4*SQRT5
                     DUM4 = DUM4*SQRT11*SQRT7*SQRT3
                  ELSE
C                     DUM4 = DUM4*SQRT5/SQRT3
                     DUM4 = DUM4*SQRT11*SQRT7*SQRT5
                  ENDIF
               END IF
            END IF
            KCNT = KCNT+1
            DK(KCNT,KG) = DUM4*AKINVV(KG)
         END DO
      END DO                    !KG
C
C     I PRIMITIVE
C
      DO IG = 1, NGI
         AI = EXI(IG)
         TMP1 = AI*CC(1,1)
         TMP2 = AI*CC(2,1)
         TMP3 = AI*CC(3,1)
         AIRAB = AI*RAB2
C
C        J PRIMITIVE
C
         DO JG = 1, NGJ
            AJ = EXJ(JG)
            AB = AI+AJ
            ABINV = ONE/AB
            DUM = AJ*ABINV*AIRAB
            IF (DUM .GT. TOL) GO TO 500
            EXPB = EXP(-DUM)
            P(1) = (TMP1+AJ*CC(1,2))*ABINV
            P(2) = (TMP2+AJ*CC(2,2))*ABINV
            P(3) = (TMP3+AJ*CC(3,2))*ABINV
            PAV(1) = P(1)-CC(1,1)
            PAV(2) = P(2)-CC(2,1)
            PAV(3) = P(3)-CC(3,1)
            CPV(1) = CC(1,3)-P(1)
            CPV(2) = CC(2,3)-P(2)
            CPV(3) = CC(3,3)-P(3)
            RR = CPV(1)*CPV(1)+CPV(2)*CPV(2)+CPV(3)*CPV(3)
C
C           ----- DENSITY FACTOR OF IJ PRIMITIVE
C
            DUM0 = ABINV*EXPB
            IJ = 0
            DO I = MINI, MAXI
               IF(I .EQ. 1) THEN
                  DUM1 = CIS(IG)
               ELSE IF(I .LE.  4) THEN
                  DUM1 = CIP(IG)
               ELSE IF(I .LE. 10) THEN
                  DUM1 = CID(IG)
                  IF(I .GE. 8 .AND. NORMI) DUM1 = DUM1*SQRT3
               ELSE IF(I .LE. 20) THEN
                  DUM1 = CIF(IG)
                  IF(I .GE. 14 .AND. NORMI) THEN
                     IF(I .LE. 19) THEN
                        DUM1 = DUM1*SQRT5
                     ELSE
C                        DUM1 = DUM1*SQRT3
                        DUM1 = DUM1*SQRT3*SQRT5
                     END IF
                  END IF
               ELSE IF(I .LE. 35) THEN
                  DUM1 = CIG(IG)
                  IF(I .GE. 24 .AND. NORMI) THEN
                     IF(I .LE. 29) THEN
                        DUM1 = DUM1*SQRT7
                     ELSE IF(I .LE. 32) THEN
C                        DUM1 = DUM1*SQRT5/SQRT3
                        DUM1 = DUM1*SQRT7*SQRT5/SQRT3
                     ELSE
C                        DUM1 = DUM1*SQRT3
                        DUM1 = DUM1*SQRT7*SQRT5
                     END IF
                  END IF
               ELSE IF(I .LE. 56) THEN
                  DUM1 = CIH(IG)
                  IF(I .GE. 39 .AND. NORMI) THEN
                     IF(I .LE. 44) THEN
                        DUM1 = DUM1*SQRT9
                     ELSE IF(I .LE. 50) THEN
C                        DUM1 = DUM1*SQRT7/SQRT3
                        DUM1 = DUM1*SQRT7*SQRT3
                     ELSE IF(I .LE. 53) THEN
C                        DUM1 = DUM1*SQRT3
                        DUM1 = DUM1*SQRT9*SQRT7
                     ELSE
C                        DUM1 = DUM1*SQRT5/SQRT3
                        DUM1 = DUM1*SQRT7*SQRT5*SQRT3
                     ENDIF
                  ENDIF
               ELSE
                  DUM1 = CII(IG)
                  IF(I .GE. 60 .AND. NORMI) THEN
                     IF(I .LE. 65) THEN
                        DUM1 = DUM1*SQRT11
                     ELSE IF(I .LE. 71) THEN
C                        DUM1 = DUM1*SQRT3
                        DUM1 = DUM1*SQRT11*SQRT3
                     ELSE IF(I .LE. 74) THEN
C                        DUM1 = DUM1*SQRT3
                        DUM1 = DUM1*SQRT11*THREE
                     ELSE IF(I .LE. 77) THEN
C                        DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)
                        DUM1 = DUM1*SQRT11*SQRT7*SQRT3/SQRT5
                     ELSE IF(I .LE. 83) THEN
C                        DUM1 = DUM1*SQRT5
                        DUM1 = DUM1*SQRT11*SQRT7*SQRT3
                     ELSE
C                        DUM1 = DUM1*SQRT5/SQRT3
                        DUM1 = DUM1*SQRT11*SQRT7*SQRT5
                     ENDIF
                  END IF
               END IF
               DUM2 = DUM0*DUM1
               DO J = MINJ, MAXJ
                  IF(J .EQ. 1) THEN
                     DUM3 = CJS(JG)
                  ELSE IF(J .LE.  4) THEN
                     DUM3 = CJP(JG)
                  ELSE IF(J .LE. 10) THEN
                     DUM3 = CJD(JG)
                     IF(J .GE. 8 .AND. NORMJ) DUM3 = DUM3*SQRT3
                  ELSE IF(J .LE. 20) THEN
                     DUM3 = CJF(JG)
                     IF(J .GE. 14 .AND. NORMJ) THEN
                        IF(J .LE. 19) THEN
                           DUM3 = DUM3*SQRT5
                        ELSE
C                           DUM3 = DUM3*SQRT3
                           DUM3 = DUM3*SQRT3*SQRT5
                        END IF
                     END IF
                  ELSE IF(J .LE. 35) THEN
                     DUM3 = CJG(JG)
                     IF(J .GE. 24 .AND. NORMJ) THEN
                        IF(J .LE. 29) THEN
                           DUM3 = DUM3*SQRT7
                        ELSE IF(J .LE. 32) THEN
C                           DUM3 = DUM3*SQRT5/SQRT3
                           DUM3 = DUM3*SQRT7*SQRT5/SQRT3
                        ELSE
C                           DUM3 = DUM3*SQRT3
                           DUM3 = DUM3*SQRT5*SQRT3
                        END IF
                     END IF
                  ELSE IF(J .LE. 56) THEN
                     DUM3 = CJH(JG)
                     IF(J .GE. 39 .AND. NORMJ) THEN
                        IF(J .LE. 44) THEN
                           DUM3 = DUM3*SQRT9
                        ELSE IF(J .LE. 50) THEN
C                           DUM3 = DUM3*SQRT7/SQRT3
                           DUM3 = DUM3*SQRT7*SQRT3
                        ELSE IF(J .LE. 53) THEN
C                           DUM3 = DUM3*SQRT3
                           DUM3 = DUM3*SQRT9*SQRT7
                        ELSE
C                           DUM3 = DUM3*SQRT5/SQRT3
                           DUM3 = DUM3*SQRT7*SQRT5*SQRT3
                        ENDIF
                     END IF
                  ELSE
                     DUM3 = CJI(JG)
                     IF(J .GE. 60 .AND. NORMJ) THEN
                        IF(J .LE. 65) THEN
                           DUM3 = DUM3*SQRT11
                        ELSE IF(J .LE. 71) THEN
C                           DUM3 = DUM3*SQRT3
                           DUM3 = DUM3*SQRT11*SQRT3
                        ELSE IF(J .LE. 74) THEN
C                           DUM3 = DUM3*SQRT3
                           DUM3 = DUM3*SQRT11*THREE
                        ELSE IF(J .LE. 77) THEN
C                           DUM3 = DUM3*SQRT7/(SQRT5*SQRT3)
                           DUM3 = DUM3*SQRT11*SQRT7*SQRT5/SQRT3
                        ELSE IF(J .LE. 83) THEN
C                           DUM3 = DUM3*SQRT5
                           DUM3 = DUM3*SQRT11*SQRT7*SQRT3
                        ELSE
C                           DUM3 = DUM3*SQRT5/SQRT3
                           DUM3 = DUM3*SQRT11*SQRT7*SQRT5
                        ENDIF
                     END IF
                  END IF
                  IJ = IJ+1
                  DB(IJ) = DUM2*DUM3
               END DO
            END DO
C
C           K PRIMITIVE
C
            DO KG = 1, NGK
               AK = EXK(KG)
C
C              ----- ROOTS AND WEIGHTS FOR QUADRATURE
C
               RHO = AB*AK/(AB+AK)
               XX = RHO*RR
               IF (NROOTS .LE. 3) CALL RT123
               IF (NROOTS .EQ. 4) CALL ROOT4
               IF (NROOTS .EQ. 5) CALL ROOT5
               IF (NROOTS .GE. 6) CALL ROOT6
C
C              ----- FORM 3D-INTEGRALS
C
               CALL XYZ3C2E(U,W,PAV,CPV,G2,A1,B1,C1,A2,B2,F00,
     *              AB,AK,LB,LK,NROOTS)
               CALL XYZHRR3C2E(G2,CAB,G3,LI,LJ,LK,LB,NROOTS)
C
C              ----- FORM 3C-2E ERIS (IJ|K)
C
               CALL FORM3C2EP(G3,DB,DK(1,KG),GG,LLI,LLJ,LLK,
     *              LI,LJ,LK,NI,NJ,NK,NROOTS)
C
            END DO              !KG
 500        CONTINUE
         END DO                 !JG
      END DO                    !IG
C
      RETURN
      END
C*MODULE RIINT   *DECK XYZ2C1OV
      SUBROUTINE XYZ2C1OV(S1,S2,A1,A2,AI,AJ,EXPFAC,LI,LJ,LTOT)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION S1(3,0:LTOT),S2(3,0:LI,0:LJ),A1(3),A2(3)
C
      PARAMETER (HALF=0.5D+00, ONE=1.0D+00, PI12=1.77245385090552D+00)
C
      AIJINV  = ONE/(AI+AJ)
      AIJINV2 = HALF*AIJINV
      S00 = PI12*SQRT(AIJINV)
      VAL = -AJ*AIJINV
C     A1(1) = CC(1,1)-CC(1,2)
C     A1(2) = CC(2,1)-CC(2,2)
C     A1(3) = CC(3,1)-CC(3,2)
      A2(1) = VAL*A1(1)
      A2(2) = VAL*A1(2)
      A2(3) = VAL*A1(3)
C
C     LJ = 0
C
C     (SS)
      S1(1,0) = S00
      S1(2,0) = S00
      S1(3,0) = S00*EXPFAC
      IF (LTOT .EQ. 0) GO TO 1000
C     (PS)
      S1(1,1) = A2(1)*S1(1,0)
      S1(2,1) = A2(2)*S1(2,0)
      S1(3,1) = A2(3)*S1(3,0)
      IF (LTOT .EQ. 1) GO TO 1000
C     (DS), (FS), (GS)...
      DO I = 2, LTOT
         VAL = (I-1)*AIJINV2
         S1(1,I) = A2(1)*S1(1,I-1)+VAL*S1(1,I-2)
         S1(2,I) = A2(2)*S1(2,I-1)+VAL*S1(2,I-2)
         S1(3,I) = A2(3)*S1(3,I-1)+VAL*S1(3,I-2)
      END DO
 1000 CONTINUE
      DO I = 0, LI
         S2(1,I,0) = S1(1,I)
         S2(2,I,0) = S1(2,I)
         S2(3,I,0) = S1(3,I)
      END DO
      IF(LJ .EQ. 0) RETURN
C
C     LJ >= 1
C
      DO J = 1, LJ
         DO I = 0, LTOT-J
            S1(1,I) = S1(1,I+1)+A1(1)*S1(1,I)
            S1(2,I) = S1(2,I+1)+A1(2)*S1(2,I)
            S1(3,I) = S1(3,I+1)+A1(3)*S1(3,I)
         END DO
         DO I = 0, LI
            S2(1,I,J) = S1(1,I)
            S2(2,I,J) = S1(2,I)
            S2(3,I,J) = S1(3,I)
         END DO
      END DO
C
      RETURN
      END
C*MODULE RIINT   *DECK XYZ2C2E
      SUBROUTINE XYZ2C2E(CC,U,W,G,A1,B1,C1,A2,B2,F00,AI,AJ,LI,LJ,NROOTS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION CC(3,2),U(NROOTS),W(NROOTS)
      DOUBLE PRECISION G(3,NROOTS,0:LI,0:LJ)
      DOUBLE PRECISION A1(3,NROOTS),B1(NROOTS),C1(NROOTS),
     *                 A2(3,NROOTS),B2(NROOTS),F00(NROOTS)
C
      PARAMETER (HALF=0.5D+00, ONE=1.0D+00)
C
      AIJINV = ONE/(AI+AJ)
      RHO = AI*AJ*AIJINV
      AIINV  = ONE/AI
      AJINV  = ONE/AJ
      AIINV2 = HALF*AIINV
      AJINV2 = HALF*AJINV
      AIJINV2 = HALF*AIJINV
      AIJINVSQ = SQRT(AIJINV)
      DO IR = 1, NROOTS
         U2 = U(IR)*RHO
         T2 = U2/(RHO+U2)
         VAL1 = RHO*AIINV*T2
         VAL2 = RHO*AJINV*T2
         VAL3 = AIINV2-AIINV2*VAL1
         VAL4 = AJINV2-AJINV2*VAL2
         A1(1,IR) = VAL1*(CC(1,2)-CC(1,1))
         A1(2,IR) = VAL1*(CC(2,2)-CC(2,1))
         A1(3,IR) = VAL1*(CC(3,2)-CC(3,1))
         A2(1,IR) = VAL2*(CC(1,1)-CC(1,2))
         A2(2,IR) = VAL2*(CC(2,1)-CC(2,2))
         A2(3,IR) = VAL2*(CC(3,1)-CC(3,2))
         B1(IR) = VAL3
         B2(IR) = VAL4
         C1(IR) = AIJINV2*T2
         F00(IR) = AIJINVSQ*W(IR)
      END DO
C
C     LJ = 0
C     (S/S)
      DO IR = 1, NROOTS
         G(1,IR,0,0) = ONE
         G(2,IR,0,0) = ONE
         G(3,IR,0,0) = F00(IR)
C         WRITE(6,'(I1,3D20.10)') IR,(G(IA,IR,0,0),IA=1,3)
      END DO
      IF (LI .EQ. 0) GO TO 1000
C     (P/S)
      DO IR = 1, NROOTS
         G(1,IR,1,0) = A1(1,IR)
         G(2,IR,1,0) = A1(2,IR)
         G(3,IR,1,0) = A1(3,IR)*F00(IR)
C         WRITE(6,'(I1,3D20.10)') IR,(G(IA,IR,1,0),IA=1,3)
      END DO
      IF (LI .EQ. 1) GO TO 1000
C     (D/S), (F/S), (G/S)
      DO I = 2, LI
         DO IR = 1, NROOTS
            DO K = 1, 3
               G(K,IR,I,0) = A1(K,IR)*G(K,IR,I-1,0)
     *                  +(I-1)*B1(IR)*G(K,IR,I-2,0)
            END DO
C            WRITE(6,'(I1,3D20.10)') IR,(G(IA,IR,I,0),IA=1,3)
         END DO
      END DO
 1000 CONTINUE
      IF (LJ .EQ. 0) RETURN
C
C     LJ = 1
C     (S/P)
      DO IR = 1, NROOTS
         G(1,IR,0,1) = A2(1,IR)
         G(2,IR,0,1) = A2(2,IR)
         G(3,IR,0,1) = A2(3,IR)*F00(IR)
C         WRITE(6,'(I1,3D20.10)') IR,(G(IA,IR,0,1),IA=1,3)
      END DO
      IF (LI .EQ. 0) GO TO 2000
C     (P/P)
      DO IR = 1, NROOTS
         DO K = 1, 3
            G(K,IR,1,1) = A1(K,IR)*G(K,IR,0,1)+C1(IR)*G(K,IR,0,0)
         END DO
C         WRITE(6,'(I1,3D20.10)') IR,(G(IA,IR,1,1),IA=1,3)
      END DO
      IF (LI .EQ. 1) GO TO 2000
C     (D/P), (F/P), (G/P)
      DO I = 2, LI
         DO IR = 1, NROOTS
            DO K = 1, 3
               G(K,IR,I,1) = A1(K,IR)*G(K,IR,I-1,1)
     *                  +(I-1)*B1(IR)*G(K,IR,I-2,1)+C1(IR)*G(K,IR,I-1,0)
            END DO
C            WRITE(6,'(I1,3D20.10)') IR,(G(IA,IR,I,1),IA=1,3)
         END DO
      END DO
 2000 CONTINUE
      IF(LJ .EQ. 1) RETURN
C
C     LJ >= 2
      DO J = 2, LJ
C        NAI = 0
         DO IR = 1, NROOTS
            DO K = 1, 3
               G(K,IR,0,J) = A2(K,IR)*G(K,IR,0,J-1)
     *                  +(J-1)*B2(IR)*G(K,IR,0,J-2)
            END DO
         END DO
         IF (LI .EQ. 0) GO TO 3000
C        NAI = 1
         DO IR = 1, NROOTS
            DO K = 1, 3
               G(K,IR,1,J) = A1(K,IR)*G(K,IR,0,J)
     *                    +J*C1(IR)*G(K,IR,0,J-1)
            END DO
         END DO
         IF (LI .EQ. 1) GO TO 3000
C        NAI = 2
         DO I = 2, LI
            DO IR = 1, NROOTS
               DO K = 1, 3
                  G(K,IR,I,J) = A1(K,IR)*G(K,IR,I-1,J)
     *                     +(I-1)*B1(IR)*G(K,IR,I-2,J)
     *                       +J*C1(IR)*G(K,IR,I-1,J-1)
               END DO
            END DO
         END DO
 3000    CONTINUE
      END DO
C
      RETURN
      END
C*MODULE RIINT   *DECK XYZ3C2E
      SUBROUTINE XYZ3C2E(U,W,PA,CP,G,A1,B1,C1,A2,B2,F00,
     *     AI,AJ,LI,LJ,NROOTS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION U(NROOTS),W(NROOTS),PA(3),CP(3)
      DOUBLE PRECISION G(3,NROOTS,0:LI,0:LJ)
      DOUBLE PRECISION A1(3,NROOTS),B1(NROOTS),C1(NROOTS),
     *                 A2(3,NROOTS),B2(NROOTS),F00(NROOTS)
C
      PARAMETER (HALF=0.5D+00, ONE=1.0D+00)
C
      AIJINV = ONE/(AI+AJ)
      RHO = AI*AJ*AIJINV
      AIINV  = ONE/AI
      AJINV  = ONE/AJ
      AIINV2 = HALF*AIINV
      AJINV2 = HALF*AJINV
      AIJINV2 = HALF*AIJINV
      AIJINVSQ = SQRT(AIJINV)
      DO IR = 1, NROOTS
         U2 = U(IR)*RHO
         T2 = U2/(RHO+U2)
         VAL1 = RHO*AIINV*T2
         VAL2 = RHO*AJINV*T2
         VAL3 = AIINV2-AIINV2*VAL1
         VAL4 = AJINV2-AJINV2*VAL2
         A1(1,IR) = PA(1)+VAL1*CP(1)
         A1(2,IR) = PA(2)+VAL1*CP(2)
         A1(3,IR) = PA(3)+VAL1*CP(3)
         A2(1,IR) = -VAL2*CP(1)
         A2(2,IR) = -VAL2*CP(2)
         A2(3,IR) = -VAL2*CP(3)
         B1(IR) = VAL3
         B2(IR) = VAL4
         C1(IR) = AIJINV2*T2
         F00(IR) = AIJINVSQ*W(IR)
      END DO
C
C     LJ = 0
C     (S/S)
      DO IR = 1, NROOTS
         G(1,IR,0,0) = ONE
         G(2,IR,0,0) = ONE
         G(3,IR,0,0) = F00(IR)
      END DO
      IF (LI .EQ. 0) GO TO 1000
C     (P/S)
      DO IR = 1, NROOTS
         G(1,IR,1,0) = A1(1,IR)
         G(2,IR,1,0) = A1(2,IR)
         G(3,IR,1,0) = A1(3,IR)*F00(IR)
      END DO
      IF (LI .EQ. 1) GO TO 1000
C     (D/S), (F/S), (G/S)
      DO I = 2, LI
         DO IR = 1, NROOTS
            DO K = 1, 3
               G(K,IR,I,0) = A1(K,IR)*G(K,IR,I-1,0)
     *                  +(I-1)*B1(IR)*G(K,IR,I-2,0)
            END DO
         END DO
      END DO
 1000 CONTINUE
      IF (LJ .EQ. 0) RETURN
C
C     LJ = 1
C     (S/P)
      DO IR = 1, NROOTS
         G(1,IR,0,1) = A2(1,IR)
         G(2,IR,0,1) = A2(2,IR)
         G(3,IR,0,1) = A2(3,IR)*F00(IR)
      END DO
      IF (LI .EQ. 0) GO TO 2000
C     (P/P)
      DO IR = 1, NROOTS
         DO K = 1, 3
            G(K,IR,1,1) = A1(K,IR)*G(K,IR,0,1)+C1(IR)*G(K,IR,0,0)
         END DO
      END DO
      IF (LI .EQ. 1) GO TO 2000
C     (D/P), (F/P), (G/P)
      DO I = 2, LI
         DO IR = 1, NROOTS
            DO K = 1, 3
               G(K,IR,I,1) = A1(K,IR)*G(K,IR,I-1,1)
     *                  +(I-1)*B1(IR)*G(K,IR,I-2,1)+C1(IR)*G(K,IR,I-1,0)
            END DO
         END DO
      END DO
 2000 CONTINUE
      IF(LJ .EQ. 1) RETURN
C
C     LJ >= 2
C     NAI = 0
C     (S/*)
      DO J = 2, LJ
         DO IR = 1, NROOTS
            DO K = 1, 3
               G(K,IR,0,J) = A2(K,IR)*G(K,IR,0,J-1)
     *                  +(J-1)*B2(IR)*G(K,IR,0,J-2)
            END DO
         END DO
      END DO
      IF (LI .EQ. 0) GO TO 3000
C     NAI = 1
C     (P/*)
      DO J = 2, LJ
         DO IR = 1, NROOTS
            DO K = 1, 3
               G(K,IR,1,J) = A1(K,IR)*G(K,IR,0,J)
     *                    +J*C1(IR)*G(K,IR,0,J-1)
            END DO
         END DO
      END DO
      IF (LI .EQ. 1) GO TO 3000
C     NAI = 2
C     (*/*) D>=*
      DO J = 2, LJ
         DO I = 2, LI
            DO IR = 1, NROOTS
               DO K = 1, 3
                  G(K,IR,I,J) = A1(K,IR)*G(K,IR,I-1,J)
     *                     +(I-1)*B1(IR)*G(K,IR,I-2,J)
     *                       +J*C1(IR)*G(K,IR,I-1,J-1)
               END DO
            END DO
         END DO
      END DO
 3000 CONTINUE
C
      RETURN
      END
C*MODULE RIINT   *DECK XYZHRR3C2E
      SUBROUTINE XYZHRR3C2E(G,CAB,G3,LI,LJ,LK,LB,NROOTS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION G(3,NROOTS,0:LB,0:LK),CAB(3,0:LJ),
     *     G3(3,NROOTS,0:LI,0:LJ,0:LK)
C
      DO KK = 0, LK
         DO JJ = 0, LJ-1
            DO II = 0, LI
               DO IR = 1, NROOTS
                  DO IA = 1, 3
                     G3(IA,IR,II,JJ,KK) = G(IA,IR,II,KK)
                  END DO
               END DO
            END DO
            IF(LJ .GT. 0) THEN
               DO II = 0, LI+LJ-JJ-1
CC                  WRITE(6,*) 'LI=',LI,'LJ=',LJ,'II=',II,'JJ=',JJ
                  DO IR = 1, NROOTS
                     DO IA = 1, 3
                        G(IA,IR,II,KK) = G(IA,IR,II+1,KK)
     *                       +CAB(IA,1)*G(IA,IR,II,KK)
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO
      DO KK = 0, LK
         DO II = 0, LI
            DO IR = 1, NROOTS
               DO IA = 1, 3
                  G3(IA,IR,II,LJ,KK) = G(IA,IR,II,KK)
               END DO
            END DO
         END DO
      END DO
C      DO IR = 1, NROOTS
C         DO KK = 0, LK
C            DO JJ = 0, LJ
C               DO II = 0, LI
C                  WRITE(6,'(4I4,3D15.8)') II,JJ,KK,IR,
C     *                 (G3(IA,IR,II,JJ,KK),IA=1,3)
C               END DO
C            END DO
C         END DO
C      END DO
C
      RETURN
      END
