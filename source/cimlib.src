CMWSC-WL- 2010.02.18 Copy from njsubr-lw.src
CMWSC -1- 2005.11.17 ---
CMWS      subroutine wnfo_int(nfo,key,nvalu,yes)
CMWS      implicit none
CMWS      integer nfo,nvalu,i,k,l,yes
CMWS      character kline*40,key*(*)
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=200)
CMWS      enddo
CMWS
CMWS 100  write(nfo) '#'//kline(1:40)
CMWS      write(nfo) 'I '
CMWS      write(nfo) nvalu
CMWS      yes=1
CMWS 200  return
CMWS
CMWS      end

CMWSC -2- 2005.11.17 ---
CMWS      subroutine rnfo_int(nfo,key,nvalu,yes)
CMWS      implicit none
CMWS      integer nfo,nvalu,i,k,l,yes
CMWS      character kline*40,ktype*2,line*41,key*(*)
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      rewind(nfo)
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=120) line
CMWS         if (line(1:1).ne.'#') goto 120
CMWS         if (index(line,kline).eq.0) goto 120
CMWS         read(nfo,end=100,err=100) ktype
CMWS         if (ktype(1:1).ne.'I') return
CMWS         read(nfo,end=100,err=100) nvalu
CMWS         yes=1
CMWS         return
CMWS 120  enddo
CMWS
CMWS 100  return
CMWS
CMWS      end

CMWSC -3- 2005.11.18 ---
CMWS      subroutine wnfo_real(nfo,key,rvalu,yes)
CMWS      implicit none
CMWS      integer nfo,i,k,l,yes
CMWS      character kline*40,key*(*)
CMWS      real*4 rvalu
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=200)
CMWS      enddo
CMWS
CMWS 100  write(nfo) '#'//kline(1:40)
CMWS      write(nfo) 'r '
CMWS      write(nfo) rvalu
CMWS      yes=1
CMWS 200  return
CMWS
CMWS      end

CMWSC -4- 2005.11.17 ---
CMWS      subroutine rnfo_real(nfo,key,rvalu,yes)
CMWS      implicit none
CMWS      integer nfo,i,k,l,yes
CMWS      character kline*40,ktype*2,line*41,key*(*)
CMWS      real*4 rvalu
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      rewind(nfo)
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=120) line
CMWS         if (line(1:1).ne.'#') goto 120
CMWS         if (index(line,kline).eq.0) goto 120
CMWS         read(nfo,end=100,err=100) ktype
CMWS         if (ktype(1:1).ne.'r') return
CMWS         read(nfo,end=100,err=100) rvalu
CMWS         yes=1
CMWS         return
CMWS 120  enddo
CMWS
CMWS 100  return
CMWS
CMWS      end

CMWSC -5- 2005.11.18 ---
CMWS      subroutine wnfo_double(nfo,key,dvalu,yes)
CMWS      implicit none
CMWS      integer nfo,i,k,l,yes
CMWS      character kline*40,key*(*)
CMWS      real*8 dvalu
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=200)
CMWS      enddo
CMWS
CMWS 100  write(nfo) '#'//kline(1:40)
CMWS      write(nfo) 'R '
CMWS      write(nfo) dvalu
CMWS      yes=1
CMWS 200  return
CMWS
CMWS      end
CMWS
CMWSC -6- 2005.11.17 ---
CMWS      subroutine rnfo_double(nfo,key,dvalu,yes)
CMWS      implicit none
CMWS      integer nfo,i,k,l,yes
CMWS      character kline*40,ktype*2,line*41,key*(*)
CMWS      real*8 dvalu
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      rewind(nfo)
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=120) line
CMWS         if (line(1:1).ne.'#') goto 120
CMWS         if (index(line,kline).eq.0) goto 120
CMWS         read(nfo,end=100,err=100) ktype
CMWS         if (ktype(1:1).ne.'R') return
CMWS         read(nfo,end=100,err=100) dvalu
CMWS         yes=1
CMWS         return
CMWS 120  enddo
CMWS
CMWS 100  return
CMWS
CMWS      end

CMWSC -a- 2005.11.17 ---
CMWS      subroutine wnfo_iarray(nfo,key,x,y,z,dim,yes)
CMWS      implicit none
CMWS      integer nfo,x,y,z,i,j,k,l,yes
CMWS      integer dim(x,y,z)
CMWS      character kline*40,key*(*)
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=200)
CMWS      enddo
CMWS
CMWS 100  write(nfo) '#'//kline(1:40)
CMWS      write(nfo) 'II'
CMWS      write(nfo) x*y*z
CMWS      write(nfo) (((dim(i,j,k),i=1,x),j=1,y),k=1,z)
CMWS      yes=1
CMWS 200  return
CMWS
CMWS      end

CMWSC -b- 2005.11.17 ---
CMWS      subroutine rnfo_iarray(nfo,key,x,y,z,dim,yes)
CMWS      implicit none
CMWS      integer nfo,x,y,z,i,j,k,l,num,yes
CMWS      integer dim(x,y,z)
CMWS      character kline*40,ktype*2,line*41,key*(*)
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      rewind(nfo)
CMWS      do L=1,1000000
CMWS         read(nfo,end=100,err=120) line
CMWS         if (line(1:1).ne.'#') goto 120
CMWS         if (index(line,kline).eq.0) goto 120
CMWS         read(nfo,end=100,err=100) ktype
CMWS         if (ktype.ne.'II') return
CMWS         read(nfo,end=100,err=100) num
CMWS         if (num.ne.x*y*z) return
CMWS         read(nfo,end=100,err=100) (((dim(i,j,k),i=1,x),j=1,y),k=1,z)
CMWS         yes=1
CMWS         return
CMWS 120  enddo
CMWS
CMWS 100  return
CMWS
CMWS      end


CMWSC -c- 2005.11.17 ---
CMWS      subroutine wnfo_rarray(nfo,key,x,y,z,dim,yes)
CMWS      implicit none
CMWS      integer nfo,x,y,z,i,j,k,l,yes
CMWS      real*4 dim(x,y,z)
CMWS      character kline*40,key*(*)
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=200)
CMWS      enddo
CMWS
CMWS 100  write(nfo) '#'//kline(1:40)
CMWS      write(nfo) 'rr'
CMWS      write(nfo) x*y*z
CMWS      write(nfo) (((dim(i,j,k),i=1,x),j=1,y),k=1,z)
CMWS      yes=1
CMWS 200  return
CMWS
CMWS      end


CMWSC -d- 2005.11.17 ---
CMWS      subroutine rnfo_rarray(nfo,key,x,y,z,dim,yes)
CMWS      implicit none
CMWS      integer nfo,x,y,z,i,j,k,l,num,yes
CMWS      real*4 dim(x,y,z)
CMWS      character kline*40,ktype*2,line*41,key*(*)
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      rewind(nfo)
CMWS      do L=1,1000000
CMWS         read(nfo,end=100,err=120) line
CMWS         if (line(1:1).ne.'#') goto 120
CMWS         if (index(line,kline).eq.0) goto 120
CMWS         read(nfo,end=100,err=100) ktype
CMWS         if (ktype.ne.'rr') return
CMWS         read(nfo,end=100,err=100) num
CMWS         if (num.ne.x*y*z) return
CMWS         read(nfo,end=100,err=100) (((dim(i,j,k),i=1,x),j=1,y),k=1,z)
CMWS         yes=1
CMWS         return
CMWS 120  enddo
CMWS
CMWS 100  return
CMWS
CMWS      end

CMWSC -e- 2005.11.17 ---
CMWS      subroutine wnfo_darray(nfo,key,x,y,z,dim,yes)
CMWS      implicit none
CMWS      integer nfo,x,y,z,i,j,k,l,yes
CMWS      real*8 dim(x,y,z)
CMWS      character kline*40,key*(*)
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      do i=1,1000000
CMWS         read(nfo,end=100,err=200)
CMWS      enddo
CMWS
CMWS 100  write(nfo) '#'//kline(1:40)
CMWS      write(nfo) 'RR'
CMWS      write(nfo) x*y*z
CMWS      write(nfo) (((dim(i,j,k),i=1,x),j=1,y),k=1,z)
CMWS      yes=1
CMWS 200  return
CMWS
CMWS      end


CMWSC -f- 2005.11.17 ---
CMWS      subroutine rnfo_darray(nfo,key,x,y,z,dim,yes)
CMWS      implicit none
CMWS      integer nfo,x,y,z,i,j,k,l,num,yes
CMWS      real*8 dim(x,y,z)
CMWS      character kline*40,ktype*2,line*41,key*(*)
CMWS
CMWS      l=len(key)
CMWS      if (l.GE.40) then
CMWS         kline=key(1:40)
CMWS      else
CMWS         kline(1:l)=key(1:l)
CMWS         do k=l+1,40
CMWS            kline(k:k)=' '
CMWS         enddo
CMWS      endif
CMWS
CMWS      yes=0
CMWS      rewind(nfo)
CMWS      do L=1,1000000
CMWS         read(nfo,end=100,err=120) line
CMWS         if (line(1:1).ne.'#') goto 120
CMWS         if (index(line,kline).eq.0) goto 120
CMWS         read(nfo,end=100,err=100) ktype
CMWS         if (ktype.ne.'RR') return
CMWS         read(nfo,end=100,err=100) num
CMWS         if (num.ne.x*y*z) return
CMWS         read(nfo,end=100,err=100) (((dim(i,j,k),i=1,x),j=1,y),k=1,z)
CMWS         yes=1
CMWS         return
CMWS 120  enddo
CMWS
CMWS 100  return
CMWS
CMWS      end
C
C
      subroutine iwrit(iunit, key, n, M)
      implicit none
      character key*(*)
      integer iunit,n,L,k
      integer M(n)
C
      L=len(key)
      write(iunit,*) key(1:L), n
      write(iunit,'(10I10)') (M(k),k=1,n)
      write(iunit,*) '$END'
      end
      subroutine iwrit8(iunit, key, n, M)
      implicit none
      character key*(*)
      integer iunit,n,L,k
      integer M(n)
C
      L=len(key)
      write(iunit,*) key(1:L), n
      write(iunit,'(8I10)') (M(k),k=1,n)
      write(iunit,*) '$END'
      end
C
      subroutine rwrit(iunit, key, n, A)
      implicit none
      character key*(*)
      integer iunit,n,L,k
      double precision A(n)
C
      L=len(key)
      write(iunit,*) key(1:L), n
      write(iunit,'(1P,5E20.12)') (A(k),k=1,n)
      write(iunit,*) '$END'
      end
      subroutine rwrit8(iunit, key, n, A)
      implicit none
      character key*(*)
      integer iunit,n,L,k
      double precision A(n)
C
      L=len(key)
      write(iunit,*) key(1:L), n
      write(iunit,'(1P,4E20.12)') (A(k),k=1,n)
      write(iunit,*) '$END'
      end
C
C
      subroutine iread(iunit, key, n, M)
      implicit none
      character key*(*),line*200
      integer iunit,n,L,k,iyes,itry
      integer M(n)
C
      do k=1,n
         M(k)=0
      enddo
C
      iyes=0
      itry=0
      L=len(key)
 50   do while(iyes.eq.0)
         read(iunit,'(a)',end=100,err=100) line
         if (index(line,key(1:L)).ne.0) then
            iyes=1
            goto 100
         endif
      enddo
C
 100  itry=itry+1
      if (iyes.eq.0) then
         if (itry.eq.1) then
            rewind(iunit)
            goto 50
         else
            write(*,*) 'ERROR: can not find '//key(1:L)//' in .txt'
            stop
         endif
      else
         read(iunit,'(10I10)') (M(k),k=1,n)
         read(iunit,*,end=999)
      endif

 999  return

      end
      subroutine iread8(iunit, key, n, M)
      implicit none
      character key*(*),line*200
      integer iunit,n,L,k,iyes,itry
      integer M(n)
C
      do k=1,n
         M(k)=0
      enddo
C
      iyes=0
      itry=0
      L=len(key)
 50   do while(iyes.eq.0)
         read(iunit,'(a)',end=100,err=100) line
         if (index(line,key(1:L)).ne.0) then
            iyes=1
            goto 100
         endif
      enddo
C
 100  itry=itry+1
      if (iyes.eq.0) then
         if (itry.eq.1) then
            rewind(iunit)
            goto 50
         else
            write(*,*) 'ERROR: can not find '//key(1:L)//' in .txt'
            stop
         endif
      else
         read(iunit,'(8I10)') (M(k),k=1,n)
         read(iunit,*,end=999)
      endif

 999  return

      end
C
C
      subroutine rread(iunit, key, n, A)
      implicit none
      character key*(*),line*200
      integer iunit,n,L,k,iyes,itry
      double precision A(n)
C
      do k=1,n
         A(k)=0d0
      enddo
C
      iyes=0
      itry=0
      L=len(key)
 50   do while(iyes.eq.0)
         read(iunit,'(a)',end=100,err=100) line
         if (index(line,key(1:L)).ne.0) then
            iyes=1
            goto 100
         endif
      end do
C
 100  itry=itry+1
      if (iyes.eq.0) then
         if (itry.eq.1) then
            rewind(iunit)
            goto 50
         else
            write(*,*) 'ERROR: can not find '//key(1:L)//' in .txt'
            stop
         endif
      else
         read(iunit,'(1P,5E20.12)') (A(k),k=1,n)
         read(iunit,*,end=999)
      endif

 999  return

      end
      subroutine rread8(iunit, key, n, A)
      implicit none
      character key*(*),line*200
      integer iunit,n,L,k,iyes,itry
      double precision A(n)
C
      do k=1,n
         A(k)=0d0
      enddo
C
      iyes=0
      itry=0
      L=len(key)
 50   do while(iyes.eq.0)
         read(iunit,'(a)',end=100,err=100) line
         if (index(line,key(1:L)).ne.0) then
            iyes=1
            goto 100
         endif
      end do
C
 100  itry=itry+1
      if (iyes.eq.0) then
         if (itry.eq.1) then
            rewind(iunit)
            goto 50
         else
            write(*,*) 'ERROR: can not find '//key(1:L)//' in .txt'
            stop
         endif
      else
         read(iunit,'(1P,4E20.12)') (A(k),k=1,n)
         read(iunit,*,end=999)
      endif

 999  return

      end
      subroutine cwrit(iunit, key, n, A)
      implicit none
      character key*(*)
      integer iunit,n,L,k
      character*8 A(n)

      L=len(key)
      write(iunit,*) key(1:L), n
      write(iunit,'(10(2x,A8))') (A(k),k=1,n)
      write(iunit,*) '$END'
      end



      subroutine cread(iunit, key, n, A)
      implicit none
      character key*(*),line*200
      integer iunit,n,L,k,iyes,itry
      character*8 A(n)

      iyes=0
      itry=0
      L=len(key)
 50   do while(iyes.eq.0)
         read(iunit,'(a)',end=100,err=100) line
         if (index(line,key(1:L)).ne.0) then
            iyes=1
            goto 100
         endif
      end do

 100  itry=itry+1
      if (iyes.eq.0) then
         if (itry.eq.1) then
            rewind(iunit)
            goto 50
         else
            write(*,*) 'ERROR: can not find '//key(1:L)//' in .txt'
            stop
         endif
      else
         read(iunit,'(10(2x,A8))') (A(k),k=1,n)
         read(iunit,*,end=999)
      endif
 999  return

      End
