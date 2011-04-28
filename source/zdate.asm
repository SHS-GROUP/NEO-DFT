*  8 MAY 1983 - MWS - REASSEMBLED
*       THIS ROUTINE WAS TAKEN FROM THE PROGRAM SYSTEM ALIS
*       IT ALLOWS THE IBM VERSION OF GAMESS TO OBTAIN A
*       CHARACTER STRING PACKED WITH THE WALL CLOCK TIME.
*
*   FORTRAN CALL IS 'CALL ZDATE(FIELD)'
*        FIELD IS A 24 BYTE REGION THAT HOLDS ON RETURN
*             >HH:MM:SS DD MMM YY (JJJ)<
ZDATE    CSECT
         ENTRY DATE,YEAR,MONTH,DAYR,DAYMO
         SPACE 3
         SPACE 3
*     ENTRY DATE(FIELD)
DATE     BC    ALWAYS,12(15)
         DC    X'08'
         DC    CL7'DATE'
         STM   14,12,12(13)  SAVE ROUTINE
         BALR  EBASE,0   ESTABLISH BASE REGISTER
         USING *,EBASE
         LR    WORK,13
         LA    13,SAVEMAC
         ST    13,8(WORK)
         ST    WORK,4(13)
         L     OUTADDR,0(1)  GET ADDR OF OUTPUT ARRAY
         SPACE 3
*   GET TIME & DATE-PUT TIME INTO SOURCE STRING
         TIME  DEC
         SRL   0,8                     ZERO OUT LEFT BYTE
         ST    0,SOURCE                PUT TIME INTO SOURCE BYTES 2-4
*    DECODE DATE AND PUT INTO SOURCE STRING
DCDTE    BAL   INTREG2,LEAPTEST    DECODE DATE, PUT INTO SOURCE STRING
         MVC   PATTERND(28),PATTERN    COPY TEXT PATTERN FOR EDITING
         ED    PATTERND(28),SOURCE   EDIT SOURCE INTO PATTERN
         MVC   PATTERND+16(4),NAMONTH  COPY MONTH NAME INTO PATTERN
         MVI   PATTERND+23,C'('       INSERT '(' INTO PATTTERN
         MVI   PATTERND+27,C')'       INSERT ')' INTO PATTTERN
         MVC   0(24,OUTADDR),PATTERND+4 MOVE TEXT INTO OUTPUT ARRAY
         L     13,SAVEMAC+4
         BC    ALWAYS,RETURN
         DROP  EBASE
         SPACE 3
*   GET THE MONTH OF THE YEAR  (BINARY)
MONTH    BC    ALWAYS,12(15)
         DC    X'08'
         DC    CL7'MONTH'
         STM   14,12,12(13)  SAVE ROUTINE
         BALR  EBASE,0   ESTABLISH BASE REGISTER
         USING *,EBASE
         BAL   INTREG1,GOODSTUF
         ST    WORK,0(OUTADDR)
         L     13,SAVEMAC+4
         BC    ALWAYS,RETURN
         DROP  EBASE
         SPACE 3
*   GET VALUE OF THE YEAR (BINARY)
YEAR     BC    ALWAYS,12(15)
         DC    X'08'
         DC    CL7'YEAR'
         STM   14,12,12(13)  SAVE ROUTINE
         BALR  EBASE,0   ESTABLISH BASE REGISTER
         USING *,EBASE
         BAL   INTREG1,GOODSTUF
         A     XYEAR,NINTEENH
         ST    XYEAR,0(OUTADDR)
         L     13,SAVEMAC+4
         BC    ALWAYS,RETURN
         DROP  EBASE
         SPACE 3
*   GET DAY OF YEAR (BINARY)
DAYR     BC    ALWAYS,12(15)
         DC    X'08'
         DC    CL7'DAYR'
         STM   14,12,12(13)  SAVE ROUTINE
         BALR  EBASE,0   ESTABLISH BASE REGISTER
         USING *,EBASE
         BAL   INTREG1,GOODSTUF
         ST    XDAYR,0(OUTADDR)
         L     13,SAVEMAC+4
         BC    ALWAYS,RETURN
         DROP  EBASE
         SPACE 3
*   GET DAY OF MONTH (BINARY)
DAYMO    BC    ALWAYS,12(15)
         DC    X'08'
         DC    CL7'DAYMO'
         STM   14,12,12(13)  SAVE ROUTINE
         BALR  EBASE,0   ESTABLISH BASE REGISTER
         USING *,EBASE
         BAL   INTREG1,GOODSTUF
         ST    0,0(OUTADDR)
         L     13,SAVEMAC+4
         BC    ALWAYS,RETURN
         DROP  EBASE
         SPACE 3
*   SET UP INTERNAL CALL OF 'LEAPTEST'
GOODSTUF BALR  MBASE,0
         USING *,MBASE
         LR    WORK,13
         LA    13,SAVEMAC
         ST    13,8(WORK)
         ST    WORK,4(13)
         L     OUTADDR,0(1)
         TIME  DEC
         BAL   INTREG2,LEAPTEST   BRANCH TO DECODING ROUTINE
         BCR   ALWAYS,INTREG1     RETURN FROM 'GOODSTUF'
         DROP  MBASE
         SPACE 3
*   ROUTINE TO DECODE DATE (MO,DAY,...)
*     AND INSERTS STANDARD FORMAT DATE INTO THE SOURCE STRING
LEAPTEST BALR  LBASE,0
         USING *,LBASE
         LR    WORK,1                  REG 1 CONTAINS 00 YY DD DF
         N     WORK,YRMASK             ISOLATE YEAR IN BYTE 2 OF WORK
         O     WORK,YRSGN              TURN ON SIGN IN BYTE 3
         SRL   WORK,12                 ALIGN ON RIGHT BOUNDARY
         ST    WORK,TEMP+4             STORE PACKED DECIMAL
         LA    0,0                        VALUE OF YEAR (LAST 2 DIGITS)
         ST    0,TEMP                     IN TEMP (8 BYTES)
         CVB   XYEAR,TEMP              LOAD 2 DIGITS OF YEAR IN BINARY
         LR    0,XYEAR                 COPY XYEAR=REG6 TO REG0
         N     0,LEAPMASK              TEST FOR MOD(4) ON YEAR
         BC    NOTZERO,NOLEAP          IF(MOD(YEAR,4).NE.0) GOTO NOLEAP
         MVI   TABLE+11,29             CHANGE DAYS IN FEB TO 29
         SPACE 3
*
NOLEAP   SRL   WORK,4                  PUT 2 YEAR DIGITS IN BYTE 4
         STC   WORK,SOURCE+5           STORE YEAR IN BYTE 6 OF SOURCE
         SPACE 3
*
         N     1,YRSTRIP               DELETE YEAR DIGITS FROM R1
         ST    1,TEMP+4                STORE PACKED DECIMAL
         LA    WORK,0                     VALUE OF DAYOY (3 DIGITS)
         ST    WORK,TEMP                  IN TEMP (8 BYTES)
         STC   1,SOURCE+7              STORE DF IN BYTE 8 OF SOURCE
         SRL   1,8                     SHIFT OUT DF
         STC   1,SOURCE+6              STORE DD IN BYTE 7 OF SOURCE
         CVB   XDAYR,TEMP              LOAD BINARY VALUE OF DAYOY
         LR    0,XDAYR                 COPY DAYOY FROM R6 TO R0
*                                      WORK=R8 WAS ZEROED OUT ABOVE
AGAIN    LA    WORK,4(WORK)            WORK=WORK+4
         S     0,TABLE+0(WORK)         R0=R0-TABLE(WORK)
         BC    POS,AGAIN               IF(R0.GT.0) GO TO AGAIN
MOFOUND  A     0,TABLE+0(WORK)         R0 = DAY OF MONTH
         CVD   0,TEMP                  STORE DAYOM AS PK DEC IN TEMP
         L     1,TEMP+4                RELOAD DAYOM AS PK DEC IN R1
         SRL   1,4                     SHIFT OUT SIGN
         STC   1,SOURCE+4              STORE DAYOM IN BYTE 5 OF SOURCE
         L     1,NAMONTH+0(WORK)       LOAD NAME OF MONTH
         ST    1,NAMONTH               STORE CURRENT MONTH
         SRL   WORK,2                  BINARY VALUE OF MONTH IS IN WORK
         BCR   ALWAYS,INTREG2
         SPACE 3
RETURN   LM    14,12,12(13)  RESTORE ROUTINE
         MVI   12(13),X'FF'
         BCR   ALWAYS,14
*  CONSTANT DEFINITIONS AND STORAGE AREAS
         DS    0D
TEMP     DS    CL8
SAVEMAC  DS    18F
SOURCE   DS    CL8
PATTERND DS    7F
PATTERN  DC    XL4'40402021'
         DC    XL8'20207A20207A2020'
         DC    XL8'2220202240404040'
         DC    XL8'202040222020205D'
HSTRIP   DC    X'FFFFFFF0'
SGNINSRT DC    X'0000000C'
MMSSTOFF DC    X'FF00000F'
YRMASK   DC    X'00FF0000'
YRSGN    DC    X'0000C000'
LEAPMASK DC    X'00000003'
YRSTRIP  DC    X'0000FFFF'
NINTEENH DC    F'1900'
TABLE    DC    F'0'
         DC    F'31'
         DC    F'28'
         DC    F'31'
         DC    F'30'
         DC    F'31'
         DC    F'30'
         DC    F'31'
         DC    F'31'
         DC    F'30'
         DC    F'31'
         DC    F'30'
         DC    F'31'
NAMONTH  DC    CL4'    '
         DC    CL4'JAN '
         DC    CL4'FEB '
         DC    CL4'MAR '
         DC    CL4'APR '
         DC    CL4'MAY '
         DC    CL4'JUN '
         DC    CL4'JUL '
         DC    CL4'AUG '
         DC    CL4'SEP '
         DC    CL4'OCT '
         DC    CL4'NOV '
         DC    CL4'DEC '
ADDRLEAP DC    A(LEAPTEST)   ADDRESS OF LEAPTEST CODE
*
INTREG1  EQU   2
INTREG2  EQU   3
EBASE    EQU   4
XYEAR    EQU   5
XDAYR    EQU   6
OUTADDR  EQU   7
WORK     EQU   8
LBASE    EQU   9
MBASE    EQU   10
*
POS      EQU   2
NEG      EQU   4
NOTZERO  EQU   7
ALWAYS   EQU   15
         END
