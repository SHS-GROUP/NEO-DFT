      PROGRAM ADDZEROS
C 
C  Program for adding a zero row to the MO matrix (expanding the basis set)
C  Tailored from delaomo
C   V.Kairys
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER NEWFIL*15,OLDFIL*15,WORD*8
      PARAMETER (ZERO=0.0D+00,NDIM=1500)
      DIMENSION EL(5),ARRAY(NDIM),ARNEW(NDIM),NIFA(NDIM),NPLACA(NDIM)
      WRITE(6,*)'Input vector source file name'
      READ(5,*)OLDFIL
      OPEN(7,FILE=OLDFIL,STATUS='OLD',FORM='FORMATTED')
      WRITE(6,*)'Input vector destination file name'
      READ(5,*)NEWFIL
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
      WRITE(6,*)'How many zeroes - AOs do you want to add?'
      READ(5,*)NAOS
      NLEFT=NTOT+NAOS
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
C Initialize nplaca and nifa
      DO 21 I=1,NDIM
21     NPLACA(I)=0
      DO 22 I=1,NDIM
22     NIFA(I)=0
      WRITE(6,*)' Enter the number of the AO after which ',
     *           'you want the zeroes to be inserted'
      READ(5,*)NWHERE
      do i=1,naos
         nplaca(i)=nwhere+i
      end do
C    if nifa=1 then the ao is to be the additional row of zeroes
      DO 29 I=1,NLEFT
29       NIFA(NPLACA(I))=1
C
3     FORMAT(I2,I3,1P,5E15.8)
C
C
   20  READ(7,9100,END=115)WORD
 9100  FORMAT(A8)
       IF(WORD.NE.' $VEC')GOTO 20
      WRITE(8,FMT='(5H $VEC)')
      DO 8 I=1,NVEC 
         DO 50 J=1,NLIN-1
50          READ(7,3)IN,IA,(ARRAY((J-1)*5+K),K=1,5)
c just one last line
         READ(7,3)IN,IA,(ARRAY((NLIN-1)*5+K),K=1,NELEMI)
c
c
c
c    redoing array into arnew
         K=1
         DO 51 IAR=1,NLEFT 
         IF(NIFA(IAR).EQ.1)THEN
            ARNEW(IAR)=ZERO
         ELSE
            ARNEW(IAR)=ARRAY(K)
            K=K+1
         END IF
 51      CONTINUE
c  writing arnew into a new vec file
         DO 53 IA=1,NLINL-1
53      WRITE(8,3)MOD(I,100),IA,(ARNEW((IA-1)*5+II),II=1,5)
c just one line more
        WRITE(8,3)MOD(I,100),NLINL,(ARNEW((NLINL-1)*5+II),II=1,NELEM)
8     CONTINUE
      WRITE(8,FMT='(5H $END)')
      WRITE(6,*)'Finished...' 
      GOTO 116
115   WRITE(6,*)' $VEC GROUP NOT FOUND'
      GOTO 116
116   CLOSE(7)
      CLOSE(8)
      STOP
      END
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
