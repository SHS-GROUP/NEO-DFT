      PROGRAM CHGLIB
      REAL*8 X
      Character*40 FSTR,INPFILE
      CHARACTER*80 A
      CHARACTER*2 XKEY(10),SYM1,SYM2,C
      DATA XKEY/' s',' x',' y',' z','xx','yy','zz','xy','xz','yz'/
      call getarg(1,FSTR)
      do i=1,40
      if(FSTR(i:i).eq." ") then
         LENFIL=i-1
         goto 1
      endif
      enddo
 1    Continue
      INPFILE=FSTR(1:LENFIL)
      OPEN(UNIT=5,FILE=INPFILE,STATUS='OLD',
     1  ACCESS='SEQUENTIAL',FORM="FORMATTED")
 2    READ(5,'(A80)',END=3) A
      READ(A(11:12),'(A2)') SYM1
      L=0
      DO IKEY=1,10
      IF(SYM1.eq.XKEY(IKEY)) THEN
        L=1
        READ(A(1:22),'(I4,I4,A2,A2,F10.5)') J,K,C,SYM2,X
        WRITE(6,100) J,K,C,SYM2,X
      ENDIF
      ENDDO
      do i=80,1,-1
      if(A(i:i).ne." ") then
         LENFIL=i
         goto 4
      endif
      enddo
 4    Continue
      IF(L.eq.0) WRITE(6,'(80A1)') (A(i:i),i=1,LENFIL)
      GOTO 2
 3    CONTINUE
 100  FORMAT(I4,I3,1X,A2,4X,A2,F10.5)
      END

