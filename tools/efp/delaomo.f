c  99-08-20  add info to the output data file
c  99-06-28  increase number of AOs to 2000
c  99-04-23  add minor corrections in wording
C  99-02-23  add cntmos subroutine which counts mos and aos
C  99-02-10  read-in of filenames from standard i/o made formatted
C  98-07-02  skips part of the program if none of the vecs are to be removed 
C                  same for aos; fixed vector file format (1P added)
C 
C  Program for deleting a vector and/or ao from the $vec containing file
C   V.Kairys
C
C  Limits: 2000 Basis functions, 500 AOs or MOs can be deleted 
C
      PROGRAM DELAOMO_RANGE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER NEWFIL*255,OLDFIL*255,WORD*8
      DIMENSION NPLACE(500),NIFL(2000),EL(5),ARRAY(2000),ARNEW(2000),
     &NIFA(2000),NPLACA(500)
      WRITE(6,*)'Enter name of a file which contains a $VEC group'
      READ(5,1000)OLDFIL
      OPEN(7,FILE=OLDFIL,STATUS='OLD',FORM='FORMATTED')
      WRITE(6,*)'Enter destination file name'
      READ(5,1000)NEWFIL
      OPEN(8,FILE=NEWFIL,STATUS='NEW',FORM='FORMATTED')
c      WRITE(6,*)'How many vectors do you have in the vector file?'
c      READ(5,*)NVEC
c      WRITE(6,*)'What is the total number of basis functions?'
c      READ(5,*)NTOT
      call cntmos(7,nvec,ntot)
C   Calculates the number of lines for one vector   
C   and also # of elements in the last line
      IF(MOD(NTOT,5).EQ.0)THEN
        NLIN=NTOT/5
        NELEMI=5
      ELSE
        NLIN=INT(NTOT/5)+1
        NELEMI=MOD(NTOT,5)
      END IF
C
      WRITE(6,*)'How many MOs do you want to delete? Enter 0 if none.'
      READ(5,*)NMOS
      WRITE(6,*)'How many AOs do you want to delete? Enter 0 if none.'
      READ(5,*)NAOS
      NLEFT=NTOT-NAOS
C   Calculates the number of lines in the new vector file
C   and also # of elements in the last line
      IF(MOD(NLEFT,5).EQ.0)THEN
        NLINL=NLEFT/5
        NELEM=5
      ELSE
        NLINL=INT(NLEFT/5)+1
        NELEM=MOD(NLEFT,5)
      END IF
C
C Initialize nplace and nifl
      DO 1 I=1,2000
1     NPLACE(I)=0
      DO 2 I=1,2000
2     NIFL(I)=0
      if(nmos.gt.0)then
      WRITE(6,*)'Write numbers of the MO-s which you'
      WRITE(6,*)' want to delete. Write all, separated by commas'
      WRITE(6,*)' For example, 1,5,15,27 <Enter>'
      READ(5,*)(NPLACE(I),I=1,NMOS)
      end if
C    if nifl=1 then the vector is to be removed
      DO 9 I=1,NMOS
9       NIFL(NPLACE(I))=1
C
C Initialize nplaca and nifa
      DO 21 I=1,500
21     NPLACA(I)=0
      DO 22 I=1,500
22     NIFA(I)=0
      if(naos.gt.0)then
      WRITE(6,*)'Enter the range of the AOs which you'
      WRITE(6,*)' want to delete. First enter the first member'
      WRITE(6,*)' and then the last, separated by a comma.'
      WRITE(6,*)' Example: 94,145 <Enter>'
      READ(5,*)NFIRST,NLAST
      naos1=nlast-nfirst+1
      if(naos1.ne.naos)then
        write(6,*)' Error: your range contains ',naos1,' AOs'
        goto 116
      end if
      do i=1,naos
         nplaca(i)=nfirst+i-1
      end do
c      READ(5,*)(NPLACA(I),I=1,NAOS)
      end if
C    if nifa=1 then the ao is to be removed
      DO 29 I=1,NAOS
29       NIFA(NPLACA(I))=1
C
3     FORMAT(I2,I3,1P,5E15.8)
C
C
   20  READ(7,9100,END=115)WORD
 9100  FORMAT(A8)
       IF(WORD.NE.' $VEC')GOTO 20
      WRITE(8,9200)OLDFIL,NAOS,NMOS
9200  FORMAT('THIS FILE IS CREATED FROM ',A15,' BY REMOVING ',I3,
     *       ' AO-S AND ',I3,' MO-S')
      WRITE(8,FMT='(5H $VEC)')
      IVC=0
      DO 8 I=1,NVEC 
        IF(NIFL(I).EQ.1)THEN
         DO 11 J=1,NLIN
C  Read but no write!
11          READ(7,3)IN,IA,EL
        ELSE
         IVC=IVC+1
         DO 50 J=1,NLIN-1
50          READ(7,3)IN,IA,(ARRAY((J-1)*5+K),K=1,5)
c just one last line
         READ(7,3)IN,IA,(ARRAY((NLIN-1)*5+K),K=1,NELEMI)
c
c
c
c    redoing array into arnew
         IAR=0
         DO 51 K=1,NTOT 
         IF(NIFA(K).EQ.1)GOTO 52
         IAR=IAR+1
         ARNEW(IAR)=ARRAY(K)
 52      CONTINUE
 51      CONTINUE
         if(iar.ne.nleft)then
           write(6,*)'ERROR: iar.neq.nleft'
         end if
c  writing arnew into a new vec file
         DO 53 IA=1,NLINL-1
53      WRITE(8,3)MOD(IVC,100),IA,(ARNEW((IA-1)*5+II),II=1,5)
c just one line more
        WRITE(8,3)MOD(IVC,100),NLINL,(ARNEW((NLINL-1)*5+II),II=1,NELEM)
        END IF
8     CONTINUE
      WRITE(8,FMT='(5H $END)')
      WRITE(6,*)'Finished...' 
      GOTO 116
115   WRITE(6,*)' $VEC GROUP NOT FOUND'
      GOTO 116
116   CLOSE(7)
      CLOSE(8)
      STOP
1000  FORMAT(A255)
      END
C
      SUBROUTINE CNTMOS(IUNIT,NMOS,NAOS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER CSTRNG*80,WORD*8
      DIMENSION DMY(5)
c
c
   20 READ(IUNIT,9100,END=115)WORD
      IF(WORD.NE.' $VEC')GOTO 20
c
c  scan records
c
      nmos=0
      numlv=0
      line=0
      iaprev=1
      nldone=0
      nfdone=0
      nleft=5
   30 READ (IUNIT,9068,END = 1400) CSTRNG
      line = line + 1
c  since len_trim isn't a standard fortran77 function, we have to do without it
      do  40 i=80,1,-1
          IF (CSTRNG(i:i) .NE. ' ') go to 50
 40   continue
 50   lentrm=i
c
ccc      write(6,*)'line',line,'lentrm',lentrm
      if(nfdone.eq.0 .and. lentrm.gt.10 .and. lentrm.lt.80) then
          nleft=(lentrm-5)/15
          nfdone=1
      end if
c
      if(lentrm.lt.10)then
c  end of vectors has been reached
        linem1=line-1
        if( mod(linem1,numlv) .ne. 0) then
           write(6,*)' ERROR: Vec File is irregular:' 
           write(6,*)' the total number of lines cannot be divided by ' 
           write(6,*)' the number of lines per vector without remainder' 
           go to 116
        end if
        nmos=linem1/numlv
        goto 114
      end if
c
      if(nldone.gt.0)goto 30
c   check out the first two integers in cstrng
      READ (UNIT=CSTRNG,FMT='(I2,I3,5E15.8)') IA,IB,DMY(1)
ccc     write(6,*)'iaprev',iaprev,'ia',ia
      if(iaprev.ne.ia)then
          numlv=line-1
          nldone=1
          naos=(numlv-1)*5+nleft
      end if
      iaprev=ia
      goto 30
cc
114   write(6,*)'THE PROGRAM FOUND ',NMOS,' MOS, ',naos,' AOS'
      GOTO 116
115   WRITE(6,*)' $VEC GROUP NOT FOUND'
      GOTO 116
1400  WRITE(6,*)' EOF FILE ENCOUNTERED'
 116  REWIND(IUNIT)
      RETURN
9068  FORMAT(A80)
9100  FORMAT(A8)
      END
