C --- Update 2007.02.02 14:20 by liwei
C --- Update 2008.03.25 13:41 by liwei
      program OMPJOB
      implicit none
      integer hlp,inp,io,tmp,logf
      parameter (hlp=10,inp=11,io=12,tmp=16,logf=26)
      character*100 inname,outname,tmpname, logname,machname, pfxname
      character*100 snfname
      integer i,j,k,l,m,n,k1,k2,k3,k4,exlin,kk1,kk2,jj,kk,ii
c
      integer(kind=4) Wall0,Wall,TIME
      real*8 Tim0,CPUTim
      integer system
      external CPUTim,system,TIME
c
      character line*100,command*500
c
C +++ For OpenMP +++
      integer maxcpus,myrk
      integer(kind=4) np
      integer omp_get_num_procs,omp_get_thread_num
      external omp_get_num_procs,omp_get_thread_num
C +++ End for OpenMP
C
C
C --- For command lines
      integer marg
      parameter (marg=20)
      integer narg,typarg(0:marg)
      character*100 arg(0:marg)
      character suffix*10
C
C --- For parallel ---
      logical lg_in,lg_out,lg
      integer iyes,njobs,irerun,mjobs,merr,merglog
      character(len=100),allocatable:: jobfile(:),outfile(:)
      character(len=500),allocatable:: jobcmd(:)
      character(len=500),allocatable:: postcmd(:)
      integer,allocatable:: jobdone(:)
C
C --- 01. Initial variables of cpu time and wall time
      Tim0=CPUTim(0)
      Wall0=TIME()
C
C --- 02. Read strings after command and online help
      call NJ_cmdline(narg,arg,typarg)
      call NJ_help(narg,arg,typarg,i)
      if (i.ge.0) call ompjob_title(0)
C
C --- 03. Define used file name based on input file
      call LS_infile(narg,arg,typarg,inname,suffix)
      if (suffix.ne.'job       ') inname=trim(inname)//'.job'
      call NJ_trim(inname,k1,k2); pfxname=inname(k1:k2-4)
      call NJ_outfile(narg,arg,typarg,inname,'o','omp',outname)
      call NJ_outfile(narg,arg,typarg,inname,'l','log',logname)
      call NJ_outfile(narg,arg,typarg,inname,' ','snf',snfname)
C
C --- 04. Open needed file
      inquire(file=inname,exist=lg)
      if (lg) then
         open(inp,file=inname,status='old')
      else
         write(*,*) trim(inname)//' is not existed'; stop
      endif
C
      call NJ_argfind(narg,arg,typarg,'append',iyes)
      if (iyes==1) then
         call LS_openf(-io,outname,'formatted',1)
      else
         call LS_openf(io,outname,'formatted',1)
      endif
C
      call NJ_date(io,'Task begin from:')
      call NJ_sysinfo2(io,snfname)
      write(io,*)
C
      call NJ_argvalu(narg,arg,typarg,'np',np)
      call NJ_argfind(narg,arg,typarg,'merge',merglog)
C
C --- 05. Set the number of CPUs for OpenMP parallel
      maxcpus=omp_get_num_procs()
      write(io,'('' max cpus of machine ='',i4)') maxcpus
      if (np==0) then
         np=1
      elseif (np>maxcpus) then
         write(io,'('' max cpus is less than'',i4)') np
         np=maxcpus
      endif
      call omp_set_num_threads(np)
      write(io,'('' omp_set_num_threads ='',i4)') np
      write(io,*)
C
C --- 06. Read the task command and files
      njobs=0
      read(inp,'(a)') line
      read(line,*,end=300,err=300) njobs,irerun
 300  if (njobs<=0) then
         write(io,*) 'Error: No job found! Please check '//trim(inname)
         call ompjob_stop(io,Tim0,Wall0)
      endif
      write(io,'('' The number of jobs = '',i4)') njobs
      if (merglog.ne.0) then
         write(io,*) 'Merge all main outfile to '//trim(logname)
      endif
!     write(io,'('' Retry dead jobs for'',i2,'' times'')') irerun
C
      allocate(jobcmd(njobs),jobfile(njobs),outfile(njobs))
      allocate(postcmd(njobs),jobdone(njobs))
      jobfile=' '
      outfile=' '
      jobcmd=' '
      postcmd=' '
      jobdone=0
C
      mjobs=0
      merr=0
      do i=1,njobs
         read(inp,'(1x,a)')jobfile(i)
         read(inp,'(1x,a)')jobcmd(i)
         read(inp,'(1x,a)')outfile(i)
         call NJ_trim(outfile(i),k1,k2)
         do j=k1,k2
            if (outfile(i)(j:j)==' ') then
               outfile(i)(j:k2)=' '
               exit
            endif
         enddo
         read(inp,'(a)')postcmd(i)
C
         inquire(file=jobfile(i),exist=lg_in)
         inquire(file=outfile(i),exist=lg_out)
         if (lg_in) then
            if (lg_out) then
               jobdone(i)=1
               write(io,'(i4,2x,a)')i,trim(jobcmd(i))//' ! Done'
            else
               mjobs=mjobs+1
               jobdone(i)=0
               write(io,'(i4,2x,a)')i,trim(jobcmd(i))
            endif
         else
            if (lg_out) then
               jobdone(i)=1
               write(io,'(i4,2x,a)')i,trim(jobcmd(i))//' ! Done'
            else
               merr=merr+1
               jobdone(i)=0
               write(io,'(i4,2x,a)')i,trim(jobcmd(i))//' ! No input'
            endif
         endif
      enddo
      write(io,*)
C
      if (mjobs==0) then
         write(io,*) 'No job is needed to be run! Exiting ...'
         goto 9996
      endif
C
      if (merr.ne.0) then
         write(io,*) 'No input files for some jobs! Stopping ...'
         call ompjob_stop(io,Tim0,Wall0)
      endif
C
C --- 07. Enter the OpenMP parallel
      call LS_omprun(io,np,njobs,jobfile,outfile,jobcmd,postcmd,jobdone)
C
      merr=0
      do i=1,njobs
         if (jobdone(i).ne.1) merr=merr+1
      enddo
      if (merr.ne.0) then
         write(io,'('' The number of dead jobs ='',i5)') merr
         write(io,*)
         call ompjob_stop(io,Tim0,Wall0)
      endif
C
C --- 08. Merge all out file to .log
 9996 if (merglog.ne.0) then
         open(logf,file=logname)
         
         do i=1,njobs
            open(tmp,file=trim(jobfile(i)))
            close(tmp,status='delete')
            open(tmp,file=trim(outfile(i)))
            call NJ_copy(tmp,logf)
            close(tmp,status='delete')
         enddo
         close(logf)

         write(io,*) 'Main outfile have been merged to '//trim(logname)
      endif
      write(io,*)
C
C     ----- End Main Body ----------------------------------------------
9999  deallocate(jobfile,outfile,jobcmd,postcmd)
      call ompjob_end(io,Tim0,Wall0)
C
      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_cmdline  --  get command line arguments   ##
c     ##  2005.10.15 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
C
c     arguments are stored at arg(1:narg)
c     arg(0)= initial command; typarg(0)=0
C     default typarg=1
C     For -x (x=A~Z or a-z) typarg=2
C     For -x (x=0~9) typarg=3
C     For -x (x=-) typarg=4
C     For -x (if x is not in "a-z,A-Z,0-9,-" then typarg=-1
c     
      subroutine NJ_cmdline(narg,arg,typarg)
      implicit none
      integer marg
      parameter (marg=20)
      integer narg,typarg(0:marg)
      character*100 arg(0:marg)
      character suffix*10

      
      integer(kind=4) i,iargc,k
      character ch1,ch2
      
      arg=' '
      narg=iargc()
      if (narg.gt.marg) narg=marg
      do i=0,narg
         call getarg (i,arg(i))
      end do
      
      typarg=1
      typarg(0)=0 

      do i=1,narg
         ch1=arg(i)(1:1)
         ch2=arg(i)(2:2)
         if (ch1.eq.'-') then
            k=ichar(ch2)
            if (k.ge.65.and.k.le.90.or.k.ge.97.and.k.le.122) then
               typarg(i)=2    ! -"a-z,A-Z"
            elseif (k.ge.48.and.k.le.57) then
               typarg(i)=3    ! -"0-9"
            elseif (k.eq.45) then
               typarg(i)=4      ! -"-"
            else
               typarg(i)=-1
            endif
         endif
      enddo

      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_help  --  Print help if '-h' or '--h'     ##
c     ##  Recently update on Oct 15, 2005 by Wei Li               ##
c     ##############################################################
c     
      subroutine NJ_help(narg,arg,typarg,hlp)
      implicit none
      integer marg
      parameter (marg=20)
      integer narg,typarg(0:marg)
      character*100 arg(0:marg)
      character suffix*10

      integer i,j,k,hlp
      
      if (narg.eq.0) then
         write(*,*) 'No any string after command!'
         stop 
      endif
      
      i=typarg(1)
      if (i.eq.2.and.arg(1)(2:2).eq.'h') then
         hlp=0
      elseif (i.eq.4) then
         hlp=1
      else
         hlp=-1
      endif
      
      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_infile  --  get input file name from cmd  ##
c     ##  2005.10.15 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
c
      subroutine LS_infile(narg,arg,typarg,inname,suffix)
      implicit none
      integer marg
      parameter (marg=20)
      integer narg,typarg(0:marg)
      character*100 arg(0:marg)
      character suffix*10

      character inname*(*)
      integer i,j0,j1,j,k,k1,k2

      do i=1,narg
         j0=typarg(i-1)
         j1=typarg(i)
         if (j1.eq.1.and.(j0.eq.0.or.j0.eq.1)) then
            do j=len(arg(i)),1,-1
               if (arg(i)(j:j).ne.' ') exit
            enddo
            k=len(inname)
            j=min(j,k)
            inname=arg(i)(1:j)
C
            suffix='          '
            call NJ_trim(inname,k1,k2)
            do k=k2,max(k1,k2-10),-1
               if (inname(k:k)=='.') exit
            enddo
            suffix=inname(k+1:k2)
C
            return
         endif
      enddo

      write(*,*) 'Wrong command for running input file!'
      stop
      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_outfile -- get output file name from cmd  ##
c     ##  2005.10.15 by Wei Li; Update 2005.10.17 by Wei Li       ##
c     ##############################################################
c
      subroutine NJ_outfile(narg,arg,typarg,inname,ch,ch2,outname)
      implicit none
      integer marg
      parameter (marg=20)
      integer narg,typarg(0:marg)
      character*100 arg(0:marg)
      character suffix*10

      character inname*(*),outname*(*),ch,ch2*(*)
      integer i,j0,j1,j,k,l
      
      if (ch==' ') goto 100
      do i=1,narg
         j0=typarg(i-1)
         j1=typarg(i)
         if (j1==1.and.j0==2.and.arg(i-1)(2:2)==ch
     &                   .and.arg(i-1)(3:3)==' ') then
            do j=len(arg(i)),1,-1
               if (arg(i)(j:j).ne.' ') exit
            enddo
            k=len(outname)
            j=min(j,k)
            outname=arg(i)(1:j)
            return
         endif
      enddo

 100  j=len(inname)
      do i=j,1,-1
         if(inname(i:i).eq.'.') exit
      enddo
      
      k=len(outname)
      l=len(ch2) 
      i=min(i,k-l)
      outname=inname(1:i-1)//'.'//ch2
         
      end
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_help  --  Print help if '-h' or '--h'     ##
c     ##  Recently update on Oct 15, 2005 by Wei Li               ##
c     ##############################################################
c
      subroutine NJ_argvalu(narg,arg,typarg,ch,np)
      implicit none
C     
      integer marg 
      parameter (marg=20)
      integer narg,typarg(0:marg),lth,ith
      character*100 arg(0:marg)
      character ch*(*)
      integer i,j,k,np,j0,j1
C     
      np=0
      lth=len(ch)
C     
      do i=1,narg
         j0=typarg(i-1)
         j1=typarg(i)
         if (j1.eq.1.and.j0.eq.2.and.arg(i-1)(2:lth+1).eq.ch) then
            ith=len(arg(i))
            read(arg(i)(1:ith),*) np
            return
         endif
      enddo
C     
      end
      
c
C --- 2006.12.13 Add
      subroutine NJ_argfind(narg,arg,typarg,ch,yes)
      implicit none
C     
      integer marg 
      parameter (marg=20)
      integer narg,typarg(0:marg),lth,ith
      character*100 arg(0:marg)
      character ch*(*)
      integer i,j,k,np,j0,j1,yes
C
      yes=0
      lth=len(ch)
C
      do i=1,narg
         j0=typarg(i)
         if (j0.eq.2.and.arg(i)(2:lth+1).eq.ch) then
            yes=1
            return
         endif
      enddo
C
      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_openf  -- open an file, and backup old    ##
c     ##  2004.12.22 by Wei Li; Update 2005.10.17 by Wei Li       ##
c     ##############################################################
c              
      subroutine LS_openf(io,fname,forms,nblk)
      implicit none
      integer io,i,k1,k2,k3,system,ioo,ierr,iform,nblk
      character fname*(*),ch,run*256,forms*(*)
      parameter (ch='~')
      logical logi
      external system
C     
      ioo=abs(io)
      ierr=0
C
      if (nblk<0.or.nblk>10) nblk=1
C     
      if (ioo==0.or.ioo>10000) then
         write(*,*) 'Error: file unit==0 or >10000'
         ierr=1
         stop
      endif
C     
      call NJ_trim(fname,k1,k2)
      k3=len(forms)
      iform=0
      if (index(forms(1:k3),'unformat').ne.0) iform=1
C
      inquire(file=fname(k1:k2),exist=logi)
C
cc    run='mv '//fname(k1:k2)//' '//fname(k1:k2)//ch
C
      if (logi) then
         if (io>0) then
            i=system(run)
            open(ioo,file=fname(k1:k2),form=forms(1:k3))
         elseif (io<0) then
            open(ioo,file=fname(k1:k2),form=forms(1:k3))
            call LS_toend(ioo,iform)
            do i=1,nblk
               if (iform==0) then
                  write(ioo,*)
               elseif (iform==1) then
                  write(ioo)
               endif
            enddo
         endif
      else
         open(ioo,file=fname(k1:k2),form=forms(1:k3))
      endif
C
      end
C
C --- goto the end of a file
      subroutine LS_toend(io,iform)
      implicit none
      integer i,j,io,iform

      do
        if (iform==0) then
           read(io,*,end=200,err=200)
        elseif (iform==1) then
           read(io,end=200,err=200)
        endif
      enddo

 200  return
      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_copy  --  copy one file into another      ##
c     ##  2005.10.15 by Wei Li; Update 2005.10.19 by Wei Li       ##
c     ##############################################################
c     
      subroutine NJ_copy(file1a,file2)
      implicit none
      integer i,j,file1a,file1,file2,k,mode
      character line*256,fname*100
C
      mode=file1a
      file1=abs(file1a) 
C        
      inquire(file1,name=fname)
      k=0
      rewind(file1)
      do i=100,1,-1
         if (fname(i:i).ne.' ') exit
      enddo 
C           
      if (mode>0) then
         write(file2,'(a)') '[File] '//fname(1:i)
         write(file2,'(72(''=''))')
      endif
      do
         read(file1,'(a)',err=100,end=100) line
         do i=256,1,-1
            if (line(i:i).ne.' ') exit
         enddo
         write(file2,'(a)') line(1:i)
         k=1
      enddo

 100  if (k==0.and.mode>0) then
         write(file2,*) '<-- Warning: No file could be found! -->'
      endif
      if (mode>0) then
         write(file2,'(72(''=''))')
         write(file2,*)
      endif

      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_tmpunit -- Auto find an unused file unit  ##
c     ##  2005.10.18 by Wei Li; Update 2005.10.18 by Wei Li       ##
c     ##############################################################
c
      subroutine NJ_tmpunit(tmp2)
      implicit none
      integer tmp2
      logical lg

      if (tmp2<=0 .or. tmp2> 10000) tmp2=50
100   inquire(unit=tmp2,opened=lg)
      if (lg) then
         tmp2=tmp2+1
         goto 100
      endif

      end
C
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_trim  --  move blank of two sides         ##
c     ##  2005.01.07 by Wei Li; Update 2005.11.01 by Wei Li       ##
c     ##############################################################
c     
      subroutine NJ_trim(line,k1,k2)
      implicit none
      integer k1,k2,i,j
      character line*(*)

      j=len(line)
      if (j<=0) then
         k1=1
         k2=1
         return
      endif
C
      do i=1,j
         if (line(i:i).ne.' ') then
            k1=i; exit
         endif
         if (i==j) then
            k1=1; k2=1 
            return
         endif
      enddo
C
      do i=j,1,-1
         if (line(i:i).ne.' ') then
            k2=i; exit
         endif
      enddo

      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_date -- Print current time for system     ##
c     ##  2004.12.24 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
      
      subroutine NJ_date(io,note)
      implicit none
      integer io,i
      character datim*26,note*(*)
      
      i=len(note)
      call GDate(datim)
      write (io,*) note(1:i)//' '//datim(1:24)
      end
      
      
*Deck GDate
      Subroutine GDate(Date1)
      Implicit Integer(A-Z)
C     
C     This wrapper routine either calls FDate (on bsd systems) or
C     builds the 24-character date in some other way.
C
      Character*(*) Date1
C
C#ifdef IBM_RS6K
C#define GDATE_DONE
C      Character*26 LDate
C      LDate = ' '
C      Junk = GCTime(LDate)
C      Date1 = LDate
C#endif
C#ifndef GDATE_DONE
      Call FDate(Date1)
C#endif  
      If(Len(Date1).gt.24) Date1(25:) = ' '
      Return
      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_cputim -- Print out total job CPU Time    ##
c     ##  2004.12.24 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
c     in parent program def. "real*8 Tim0,CPUTim"; "Tim0=CPUTim(0)" for initial time

      subroutine NJ_cputim(IOut,RefTim)
      Implicit Real*8(A-H,O-Z)

 1000  Format(' CPU time: ',I3,' days ',I2,' hours ',I2,' minutes ',
     $    F4.1,' seconds.')
      
      Time = CPUTim(0) - RefTim
      NDays = (Time / (3600.0d0*24.0d0))
      Time = Time - (NDays*(3600.0d0*24.0d0))
      NHours = (Time / 3600.0d0)
      Time = Time - (NHours*3600.0d0)
      NMin = (Time / 60.0d0)
      Time = Time - (NMin*60.0d0)
      Write(IOut,1000) NDays, NHours, NMin, Time
      Return
      End
      
      function CPUTim(Junk)
      Implicit Real*8(a-h,o-z)
C     
C#ifdef IBM_RS6K
C#define CPUTIM_DONE
C      Integer IT(4), Times
C      IDum = Times(IT)
C      CPUTim = DFloat(IT(1)+IT(2))/100.0d0
C#endif
C#ifdef IBM_PC
C#define CPUTIM_DONE
C      CPUTim = PCTIME ()
C#endif
C#ifndef CPUTIM_DONE
C#ifdef _SGI64_
C      Real*4 TimArray(2), ETime
C#else
      Real TimArray(2), ETime
C#endif
      CPUTim = ETime(TimArray)
C#endif
      Return
      End
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_walltim -- Print total job wall Time      ##
c     ##  2004.12.24 by Wei Li; Update 2006.03.06 by Wei Li       ##
c     ##############################################################
c     in parent program def. "integer(kind=4) Wall0,TIME"; "Wall0=TIME()" for initial time

      subroutine NJ_walltim(io,Wall0)
      Implicit Real(kind=8) (A-H,O-Z)
      integer(kind=4) TIME,Wall,Wall0  !2006.03.06 Add (kind=4); Some other program with walltime may be modified
      external TIME

 1000  Format(' WALL time:',I3,' days ',I2,' hours ',I2,' minutes ',
     $    I4,' seconds.')

      Wall = TIME()-Wall0
      NDays= Wall/(3600*24)
      Wall = Wall-NDays*(3600*24)
      NHours= Wall/3600
      Wall = Wall-NHours*3600
      NMin = Wall/60
      Wall = Wall-NMin*60
      Write(io,1000) NDays, NHours, NMin, Wall

      Return
      
      End
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_sysinfo  --  show host os & directory     ##
c     ##  2005.12.26 by Wei Li; Update 2005.12.26 by Wei Li       ##
c     ##############################################################
c
      subroutine NJ_sysinfo(io)
      implicit none
      integer io,system,i,j,k,itmp
      character line*100,line2*100,line3*100
      external system
C     
      itmp=55
      call NJ_tmpunit(itmp)
C     
C --- 01. Hostname
      i=system('hostname > systeminfo.tmp')
      open(itmp,file='systeminfo.tmp')
      read(itmp,'(a)') line
      do i=100,1,-1
         if (line(i:i).ne.' ') exit
      enddo
      close(itmp) 
C --- 02. Operation system
      j=system('uname -sp > systeminfo.tmp')
      open(itmp,file='systeminfo.tmp')
      read(itmp,'(a)') line2
      do j=100,1,-1
         if (line2(j:j).ne.' ') exit
      enddo
      close(itmp) 
C --- 03. Directory
      k=system('pwd > systeminfo.tmp')
      open(itmp,file='systeminfo.tmp')
      read(itmp,'(a)') line3
      do k=100,1,-1
         if (line3(k:k).ne.' ') exit
      enddo
      close(itmp,status='delete')
C --- Output
      write(io,*) line(1:i)//':'//line3(1:k)//' ('//line2(1:j)//')'
C
      end
C
c     ##############################################################
c     ##  COPYRIGHT (C) 2006 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_sysinfo  --  show host os & directory     ##
c     ##  2005.12.26 by Wei Li; Update 2008.03.17 by Wei Li       ##
c     ##############################################################
c
      subroutine NJ_sysinfo2(io,snfname)
      implicit none
      integer io,system,i,j,k,itmp
      character line*100,line2*100,line3*100,snfname*100
      external system
C
      itmp=55
      call NJ_tmpunit(itmp)
C
C --- 01. Hostname
      i=system('hostname > '//trim(snfname))
      open(itmp,file=trim(snfname))
      read(itmp,'(a)') line
      do i=100,1,-1
         if (line(i:i).ne.' ') exit
      enddo
      close(itmp)
C --- 02. Operation system
      j=system('uname -sp > '//trim(snfname))
      open(itmp,file=trim(snfname))
      read(itmp,'(a)') line2
      do j=100,1,-1
         if (line2(j:j).ne.' ') exit
      enddo
      close(itmp)
C --- 03. Directory
      k=system('pwd > '//trim(snfname))
      open(itmp,file=trim(snfname))
      read(itmp,'(a)') line3
      do k=100,1,-1
         if (line3(k:k).ne.' ') exit
      enddo
      close(itmp,status='delete')
C --- Output
      write(io,*) line(1:i)//':'//line3(1:k)//' ('//line2(1:j)//')'
C
      end
C
C --- Stop the program 2007.01.29 by Wei Li ---
      subroutine ompjob_stop(io,Tim0,Wall0)
      real*8 Tim0
      integer (kind=4) Wall0
      call NJ_date(io,'Task stop on:')
      call NJ_cputim(io,Tim0)
      call NJ_walltim(io,Wall0)
      write(io,*) 'Error Termination of OMPJOB Program!'
      stop
      end
C
C --- End the program 2007.01.29 by Wei Li ---
      subroutine ompjob_end(io,Tim0,Wall0)
      real*8 Tim0
      integer (kind=4) Wall0
      call NJ_date(io,'Task over on:')
      call NJ_cputim(io,Tim0)
      call NJ_walltim(io,Wall0)
      write(io,*) 'Normal Termination of OMPJOB Program!'
      end
C
C --- OpenMP parallel for jobs 2007.01.30 by Wei Li
      subroutine LS_omprun(io,np,njobs,jobfile,outfile,
     &           jobcmd,postcmd,jobdone)
      implicit none
      integer io,njobs,myrk,i,j,k,k1,k2,L,system,omp_get_thread_num
      integer(kind=4) np
      integer jobdone(njobs),ierr,jerr
      character*100 jobfile(njobs),outfile(njobs)
      character*500 jobcmd(njobs),postcmd(njobs)
      external system,omp_get_thread_num
      real(kind=8) slptim
C
*$OMP parallel default(private) shared(io,np,njobs,jobdone,
*$OMP&         jobfile,outfile,jobcmd,postcmd)
      myrk=0
      myrk=omp_get_thread_num()
c
      do i=1,njobs
         j=mod(myrk,4)
         slptim=j/256d0
CC       call sleep(slptim)  ! This line may cause problem, you can remove it
         if (jobdone(i)==0) then
            jobdone(i)=1
C
            write(io,701) trim(jobcmd(i)),myrk
            call flush(io)
            ierr=system(trim(jobcmd(i)))
!           if (ierr==0) then
               jerr=system(trim(postcmd(i)))
!           else
!              jerr=system('mv '//trim(outfile(i))//' '
!    &              //trim(outfile(i))//'.err')
!              write(io,702) trim(jobcmd(i)),myrk
!              jobdone(i) = -1
!           endif
C
            call flush(io)
         endif
      enddo
*$OMP end parallel
!     write(id,*)
c
 701  format(1x,a,' running in thread',i5)
 702  format(1x,a,'  FAILED in thread',i5,' *')
C
      end
C
C
C --- ompjob title 2007.01.29 by Wei Li ---
      subroutine ompjob_title(io)
C
!     write(io,*) 'H==================================================H'
!     write(io,*) 'H  Lower Scaling Quantum Chemistry (LSQC) Program  H'
!     write(io,*) 'H  Copyright (C) 2006 - 2007 Nanjing Univ., China  H'
!     write(io,*) 'H  itcc.nju.edu.cn/lsqc  Email: shuhua@nju.edu.cn  H'
!     write(io,*) 'H==================================================H'
!     write(io,*)
      write(io,*) '          ---------------------------'
      write(io,*) '          O M P J O B -- 7 March 2007'
      write(io,*) '          ---------------------------'
      write(io,*) 'Submit jobs in .job with an OpenMP parallel'
C
      write(io,*)
      write(io,*) 'Input file: exam.job (defined jobs file)'
      write(io,*) 'Output files: exam.omp (general output file)'
      write(io,*) '              exam.log (merged job outfiles)'
      write(io,*) 'Usage: ompjob exam.job [-item value]'
      write(io,*) '  -append: append output file to an old .omp'
      write(io,*) '  -np k: OpenMP parallel with k CPUs'
      write(io,*) '  -merge: merge job outfiles to .log'
      write(io,*) '  -o filename: define .omp file'
      write(io,*) '  -l filename: define .log file'
      write(io,*) 'Example for exam.job:'
      write(io,*) '  (Comments after ''!'')'
      write(io,*) '-----------------------------------------------'
      write(io,*) '3                     ! the number of jobs'
      write(io,*) 'job1.gjf              ! gaussian input file'
      write(io,*) 'g03 job1.gjf          ! run gaussian03'
      write(io,*) 'job1.log              ! gaussian03 output file'
      write(io,*) '                      ! blank line'
      write(io,*) 'job2.inp              ! gamess input file'
      write(io,*) 'rungms job2.inp>& job2.gms   ! run gamess'
      write(io,*) 'job2.gms              ! gamess output file'
      write(io,*) '                      ! blank line'
      write(io,*) 'job3.in               ! molpro input file'
      write(io,*) 'molpro job3.in        ! run molpro'
      write(io,*) 'job3.out              ! molpro output file'
      write(io,*) '                      ! a blank line is needed'
      write(io,*) '-----------------------------------------------'
      write(io,*)
C
      if (io==0) stop
c
      end
C

