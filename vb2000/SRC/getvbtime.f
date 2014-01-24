      PROGRAM GETVBTIME
C
C     Program to sum VB2000 times for VB2000/GAMESS test set.
C
C     Compile with:-
C     gfortran -o getvbtime.x getvbtime.f
C
C     Put getvbtime.x on your and then run as:-
C     getvbtime.x
C     or run with full path to getvbtime.x, e.g
C     ./getvbtime.x  or
C     ../getvbtime.x or whatever
C
      IMPLICIT NONE
      CHARACTER*13 FSTR
      CHARACTER*40 RECORD
      INTEGER I,IOUT,IOVBOUT
      DOUBLE PRECISION T, TOT, TIME
C
C     To read and total the "TOTAL CPU TIME" from this data:-
C CPU TIME FOR MACROITERATION       14.330
C TOTAL CPU TIME                    14.900
C
C
      IOUT = 6
      IOVBOUT = 5
      TOT = 0.0D00
      DO I=1,25
        TIME = 0.0D00
        IF(I.EQ.1) GOTO 2
        IF(I.EQ.2) GOTO 2
        IF(I.EQ.9) GOTO 2
        FSTR(1:7)="exam-vb"
        WRITE(FSTR(8:9),'(I2.2)') I
        FSTR(10:13)=".log"
        OPEN(UNIT=IOVBOUT,FILE=FSTR,STATUS='OLD',
     1    ACCESS='SEQUENTIAL',FORM='FORMATTED')
 1      READ(IOVBOUT,11,END=3) RECORD
        IF(RECORD(1:27).NE."CPU TIME FOR MACROITERATION") GOTO 1
        READ(IOVBOUT,'(1X,A40)') RECORD
        READ(RECORD(31:40),'(F10.3)') T
        TIME = T + TIME
        GOTO 1
 3      WRITE(IOUT,'(A13,F10.3," secs")') FSTR,TIME
        TOT = TOT + TIME
 2      CONTINUE
      ENDDO
      WRITE(IOUT,'("TOTAL TIME",F10.3," secs")') TOT
 11   FORMAT(1X,A40)
      END
