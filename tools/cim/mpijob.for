* --- mpijob.for based on mpijob2.for (13 Aug 2010) ---*
* --- MPI script for submit jobs 2006.12.1 by Wei Li --- *
* --- Last modified on Dec. 13, 2006 --- *
* --- Last modified on Apr. 19, 2008 --- * scp even for NFS systems
      program mpijob2
      implicit none
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE), np, my_rank, ierr, irr, jerr
      integer i,j,k,L,m,n, L1,L2,k1,k2,k3,k4,k5,k6,k7,k8, system, task
      integer inp,io,ilg,tmp,logf, len_inp,len_io,len_ilg,len_tmp
      parameter(inp=11, io=12, ilg=13, logf=14, tmp=55)
      character*100 inname, outname, logname, tmpname, hostname, line
      character*100 pfxname, nfsname, stopname
      character command*500,mypath*500,mnpath*500
      integer len_mypath, len_mnpath, slptime, blankcyc, rerun, mpierr
      integer merglog,iyes,mjobs,merr
      parameter(slptime = 1)
      logical lg,lg_nfs,lg_in,lg_out
C
      integer marg
      parameter (marg=20)
      integer narg, typarg(0:marg)
      character*100 arg(0:marg),suffix*10
C
      real*8 Tim0,CPUTim
      integer Wall0,Wall,TIME
      external CPUTim
C
      integer njobs
      integer,allocatable:: jobdone(:)
      character(len=100),allocatable:: jobfile(:),outfile(:),nodefile(:)
      character(len=100),allocatable:: outfile1(:)
      character(len=500),allocatable:: jobcmd(:),node(:),workdir(:)
      character(len=500),allocatable:: postcmd(:)
C
      real(kind=8) ZERO, ONE, TWO, HALF
      parameter (ZERO=0d0, ONE=1d0, TWO=2d0, HALF=0.5d0)
      integer IZERO, IONE, ITWO
      parameter (IZERO=0, IONE=1, ITWO=2)
C
      integer nfs, DEBUG  ! nfs: Network File System. For nfs server nfs=1; else nfs=0
      parameter (DEBUG=0)
C
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
C
      mypath=' '; mnpath=' '; command=' '
      rerun = 0
      mpierr= 0
      nfs = 0
C
C --- 01. Read the number of jobs and allocate array ---
      if (my_rank == 0) then
         Tim0=CPUTim(0)
         Wall0=TIME()
C
         call NJ_cmdline(narg, arg, typarg)
         call NJ_help(narg,arg,typarg,i)
         if (i.ge.0) then
            call mpijob_title(0)
            goto 9998
         endif

         call LS_infile(narg, arg, typarg, inname, suffix)
         if (suffix.ne.'job       ') inname=trim(inname)//'.job'
C
!        call NJ_argfind(narg, arg, typarg,'nfs',nfs)
         call NJ_argfind(narg, arg, typarg,'merge',merglog)
         call NJ_trim(inname, k1, len_inp)
         outname = inname(1:len_inp-4)//'.mpi'
         logname = inname(1:len_inp-4)//'.log'
         tmpname = inname(1:len_inp-4)//'.tmp'
!        nfsname = inname(1:len_inp-4)//'.nfs'
         len_io  = len_inp
         len_ilg = len_inp
         len_tmp = len_inp

!        inquire(file=nfsname,exist=lg)
!        if (lg) nfs=1
         inquire(file=inname,exist=lg)
         if (lg) then
            open(inp,file=inname,status='old')
         else
            write(*,*) trim(inname)//' is not existed'; goto 9998
         endif
C
         call NJ_argfind(narg,arg,typarg,'append',iyes)
         if (iyes==1) then
            call LS_openf(-io,outname,'formatted',1)
         else
            call LS_openf(io,outname,'formatted',0)
         endif
C
         call mpijob_title(io)
         if (nfs.ne.0) then
            write(io,*) '    Using NFS (Network File system)'
         else
            write(io,*) '    Without NFS (Network File system)'
         endif
         if (merglog.ne.0) then
            write(io,*) 'Merge all main outfile to '//trim(logname)
         endif
C
cc       write(io,'('' MPI control file:'',a)') inname(1:len_inp)
cc       write(io,'('' MPI output file: '',a)') outname(1:len_io)
         read(inp,'(a)') line
         njobs=0; rerun=0
         read(line, *, err= 90, end = 90) njobs, rerun
 90      if (rerun.ne.0) then
            write(io,*) '    Rerun error jobs once if needs'
         endif
         write(io,*)
      endif
      call MPI_barrier(MPI_comm_world,ierr)
C
C --- 02. Send filename and allocate array for jobs infomations ---
      call MPI_bcast(len_inp,1,mpi_integer,0,MPI_comm_world,ierr)
      call MPI_bcast(len_io, 1,mpi_integer,0,MPI_comm_world,ierr)
      call MPI_bcast(len_tmp,1,mpi_integer,0,MPI_comm_world,ierr)
C
      call MPI_bcast(rerun,1,mpi_integer,0,MPI_comm_world,ierr)
      call MPI_bcast(nfs,1,mpi_integer,0,MPI_comm_world,ierr)
C
      call MPI_bcast(inname, 100,mpi_character,0,MPI_comm_world,ierr)
      call MPI_bcast(outname,100,mpi_character,0,MPI_comm_world,ierr)
      call MPI_bcast(tmpname,100,mpi_character,0,MPI_comm_world,ierr)
      call MPI_barrier(MPI_comm_world,ierr)
C
      if (my_rank.ne.0) then
         write(line,*) my_rank
         call NJ_trim(line,k1,k2)
         open(io,file=trim(outname)//'.'//line(k1:k2))
         tmpname=trim(tmpname)//'.'//line(k1:k2)
      endif
      call MPI_barrier(MPI_comm_world,ierr)
C
      call MPI_bcast(len_inp,1,mpi_integer,0,MPI_comm_world,ierr)
      call MPI_bcast(njobs,1,mpi_integer,0,MPI_comm_world,ierr)
      allocate(jobdone(njobs),jobfile(njobs),jobcmd(njobs))
      allocate(nodefile(0:np-1),outfile(njobs),outfile1(njobs))
      allocate(postcmd(njobs))
      jobdone=0
      jobfile=' '
      outfile=' '
      outfile1=' '
      jobcmd=' '
      postcmd=' '
      nodefile=' '
C
      do i=0, np-1
         write(line,*) i
         call NJ_trim(line,k1,k2)
         nodefile(i)=trim(inname)//'.'//line(k1:k2)
      enddo
!     call MPI_barrier(MPI_comm_world,ierr)
C
!     if (my_rank == 0) then
!        do i=0, np-1
!           open(tmp,file=nodefile(i))
!           write(tmp,*) IZERO
!           close(tmp)
!        enddo
!     else
      open(tmp,file = nodefile(my_rank))
      write(tmp,*) IZERO
      close(tmp)
!     endif
      call MPI_barrier(MPI_comm_world,ierr)
C
      mjobs=0
      merr=0
C --- 03. Read the number of file names and commands ---
      if (my_rank == 0) then
         write(io,*) 'Task to be run:'
         do i=1,njobs
            read(inp,'(a)')jobfile(i)
            read(inp,'(1x,a)')jobcmd(i)
            read(inp,'(a)')outfile(i)
C ---       2007.02.12 Added for get the first main job out file ---
            outfile1(i)=outfile(i)
            call NJ_trim(outfile1(i),k1,k2)
            do j=k1,k2
               if (outfile1(i)(j:j)==' ') then
                  outfile1(i)(j:k2)=' '
                  exit
               endif
            enddo
C ---
            read(inp,'(1x,a)')postcmd(i)
            write(io,'(i4,2x,a)')i,trim(jobcmd(i))
C
C --- 2007.02.12 Added for checking if the job has been completed before you run 
            inquire(file=jobfile(i),exist=lg_in)
            inquire(file=outfile1(i),exist=lg_out)
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
C --- 2007.02.12 Added for checking if the job has been completed before you run
         if (mjobs==0) then
            write(io,*) 'No job is needed to be run! Exiting ...'
         endif 
C
         if (merr.ne.0) then
            write(io,*) 'No input files for some jobs! Stopping ...'
         endif

      endif
      call MPI_barrier(MPI_comm_world,ierr)
C
      call MPI_bcast(mjobs,1,mpi_integer,0,MPI_comm_world,ierr)
      call MPI_bcast(merr, 1,mpi_integer,0,MPI_comm_world,ierr)
      if (mjobs==0.or.merr.ne.0) goto 9997
C
      call MPI_bcast(jobfile,njobs*100,mpi_character,0,
     &     MPI_comm_world,ierr)
      call MPI_bcast(outfile,njobs*100,mpi_character,0,
     &     MPI_comm_world,ierr)
      call MPI_bcast(jobcmd,njobs*500,mpi_character,0,
     &     MPI_comm_world,ierr)
      call MPI_bcast(postcmd,njobs*500,mpi_character,0,
     &     MPI_comm_world,ierr)
      call MPI_bcast(jobdone,njobs,mpi_integer,0,
     &     MPI_comm_world,ierr)
      call MPI_barrier(MPI_comm_world,ierr)
C
C --- 04. Read the hostname of nodes and work directory ---
      allocate(workdir(0:np-1))
      workdir = ' '
      call NJ_mypath(mypath,len_mypath,tmpname)
      workdir(my_rank)=mypath
      do i=0,np-1
         call MPI_bcast(workdir(i),100,mpi_character,i,
     &        MPI_comm_world,ierr)
      enddo
      if (my_rank == 0) then
         mnpath = mypath
         len_mnpath = len_mypath
      endif
      call MPI_barrier(MPI_comm_world,ierr)
      call MPI_bcast(len_mnpath,1,mpi_integer,0,MPI_comm_world,ierr)
      call MPI_bcast(mnpath,100,mpi_character,0,MPI_comm_world,ierr)

      call MPI_barrier(MPI_comm_world,ierr)
C
C --  Add this for printing Host and Path info
      if (my_rank == 0) then
         write(io,*) 'Host to be used:'
         write(io,'(''   M'',2x,a)') trim(workdir(0))
         do j = 1, np-1
            write(io,'(i4,2x,a)') j, trim(workdir(j))
         enddo
         write(io,*)
      endif
      call MPI_barrier(MPI_comm_world,ierr)
C     goto 9999
C
C ................................................................................
C --- 05. Submtted jobs to idle cpus ---
      if (my_rank == 0) then
         do i = 1, njobs
            if (jobdone(i). ne. 0) cycle
 190        do j = 1, np-1
cc             if (jobdone(i) == 1) exit
               open(tmp,file=nodefile(j))
               read(tmp,*) task
               close(tmp)
cc             write(*,*) i,j,task
               if (task == 0) then
                  jobdone(i) = 1
                  task = i
!!                if (nfs == 0) then
!!                   command='scp '//trim(jobfile(i))//' '
!!   &                  //trim(workdir(j))//' >& '//tmpname
!!                   ierr=system(command)
!!                endif
!                 write(io,'(a)') trim(command)
                  write(line,*) task
                  call NJ_trim(line,k1,k2)
                  command='echo '//line(k1:k2)//' > '//trim(nodefile(j))
                  ierr=system(command)
****              open(tmp,file=nodefile(j))
****              write(tmp,*) task
****              call flush(tmp)
****              close(tmp)
                  write(io,201) i, j
                  call flush(io)
****              call sleep(slptime)
!!!!              if (nfs == 0) then
                     command='scp '//trim(nodefile(j))//' '
     &                  //trim(workdir(j))//' >& '//tmpname
                     ierr=system(command)
!!!!              endif
!                 write(io,'(a)') trim(command)
**                call MPI_BSEND(task,1,mpi_integer,j,70,
**   &                 MPI_comm_world,ierr)
**                write(io,*) 'task,ierr=',task,ierr
                  exit
               endif
            enddo
C
 200        if (jobdone(i) == 0) then
               call sleep(slptime)
               call sleep(slptime)
               call sleep(slptime)
               goto 190
            endif
****        call sleep(slptime)
         enddo
         call sleep(slptime)
         write(io,*)
         task = -999
         do j=1, np-1
            command='echo -999 > '//trim(nodefile(j))
            ierr=system(command)
****        open(tmp,file=nodefile(j))
****        write(tmp,*) task
****        call flush(tmp)
****        close(tmp)
            if (nfs == 0) then
               command='scp '//trim(nodefile(j))//' '
     &            //trim(workdir(j))//' >& '//tmpname
               ierr=system(command)
            endif
!           write(io,'(a)') trim(command)
         enddo
         goto 9999
      else
**       call MPI_RECV(task,1,mpi_integer,0,70,
**   &        MPI_comm_world,status,ierr)
**       write(io,*) 'task,ierr=',task,ierr
         blankcyc=0
 210     open(tmp,file=nodefile(my_rank))
         read(tmp,*) task
         close(tmp)
         write(io,*) 'task =',task
         call flush(io)
         if (task == 0) then
            blankcyc = blankcyc + 1
            if (blankcyc<=np) then    ! 2008.04.07  60 -> 30 -> np
               call sleep(slptime)
               goto 210
            else
               goto 9999
            endif
         elseif (task>0.and.task<=njobs) then
            blankcyc=0
            ierr=system(jobcmd(task))
            write(io,'(a)') trim(jobcmd(task))
            call flush(io)
            if (ierr.ne.0) then
               write(99,*) 'task,myrank,ierr=',task,my_rank,ierr
            endif
C
            if (ierr .ne. 0. and. rerun .ne. 0) then
               call sleep(5)
!              command='mv '//trim(outfile(task))//' '
!    &                 //trim(outfile(task))//'.err'
!              jerr=system(command)
               jerr=system(jobcmd(task))
               write(io,'(''Rerun: '',a)') trim(jobcmd(task))
               if (jerr .ne. 0) then
                  write(io,*) 'Error for run this job!'
                  mpierr=1
               endif
            endif
C
!!          if (nfs == 0) then
!!             command='scp '//trim(outfile(task))//' '
!!   &            //trim(mnpath)//' >& '//tmpname
!!             ierr=system(command)
!!             write(io,'(a)') trim(command)
!!             if (ierr .ne. 0) then
!!                write(io,*) 'Error for scp ...'
!!                mpierr=1
!!             endif
!!          endif
C
****        ierr=system(postcmd(task)) 
!!          L1=index(mypath,':')
!!          L2=index(mnpath,':')
!!          if ((mypath(1:L1).ne.mnpath(1:L2)) 
!!   &            .and.nfs==0.and.mpierr.eq.0)then
!!             ierr=system('rm -rf '//jobfile(task))
!!             ierr=system('rm -rf '//outfile(task))
!!          endif
C
            open(tmp,file=nodefile(my_rank))
            read(tmp,*) task
            close(tmp)
            write(io,*) 'task =',task
            call flush(io)
            if (task<0) goto 9999
            command='echo 0 > '//trim(nodefile(my_rank))
            ierr=system(command)
****        task=0
****        write(tmp,*) task
****        call flush(tmp)
            if (nfs == 0) then
               command='scp '//trim(nodefile(my_rank))//' '
     &            //trim(mnpath)//' >& '//tmpname
               ierr=system(command)
               write(io,'(a)') trim(command)
               if (ierr .ne. 0) then
                  write(io,*) 'Error for scp ...'
                  mpierr=1
               endif
            endif
            goto 210
         else
            goto 9999
         endif
            
      endif
C
 9999 deallocate(workdir)
 9997 call MPI_barrier(MPI_comm_world,ierr)
C
C --- 08. Merge all out file to .log
 9996 if (my_rank == 0) then
        if (merglog.ne.0) then
          open(logf,file=logname)
         
          do i=1,njobs
             open(tmp,file=trim(jobfile(i)))
             close(tmp,status='delete')
             open(tmp,file=trim(outfile1(i)))
             call NJ_copy(tmp,logf)
             close(tmp,status='delete')
          enddo
          close(logf)
      
          write(io,*) 'Main outfile have been merged to '//trim(logname)
        endif
      endif

      call MPI_barrier(MPI_comm_world,ierr)
C
      if (DEBUG == 0) then
         if (my_rank == 0) then
            do i=0,np-1
               ierr = system('rm -rf '//nodefile(i))
            enddo
            write(io,*) 'All jobs have been completed!'
         else
            if (mpierr==0) then
               ierr = system('rm -rf '//nodefile(my_rank))
               close(io,status='delete')
            else
               close(io)
            endif
         endif
      endif
      call MPI_barrier(MPI_comm_world,ierr)
C
      ierr = system('rm -rf '//trim(tmpname))
      deallocate(jobdone,jobfile,jobcmd,nodefile)
      deallocate(outfile,postcmd,outfile1)
C
      call MPI_barrier(MPI_comm_world,ierr)
      if (my_rank == 0) then
         call mpijob_end(io,Tim0,Wall0)
         close(io)
      endif
C
 9998 call MPI_Finalize(ierr)
C
 201  format(1x,'Job',i4,' is running in thread',i4)
C
      end
C

c     ##############################################################
c     ##  COPYRIGHT (C) 2005 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_cmdline  --  get command line arguments   ##
c     ##  2005.10.15 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
c
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
      integer narg, typarg(0:marg)
      character*100 arg(0:marg)

      integer i,iargc,k
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
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2005 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_infile  --  get input file name from cmd  ##
c     ##  2005.10.15 by Wei Li; Update 2005.10.16 by Wei Li       ##
c     ##############################################################
c     
      subroutine LS_infile(narg,arg,typarg,inname,suffix)
      implicit none
      integer marg,ierr
      parameter (marg=20)
      integer narg, typarg(0:marg)
      character*100 arg(0:marg),suffix*10
      
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



c     ##############################################################
c     ##  COPYRIGHT (C) 2005 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_trim  --  move blank of two sides         ##
c     ##  2005.01.07 by Wei Li; Update 2005.11.01 by Wei Li       ##
c     ##############################################################
c
      subroutine NJ_trim(line,k1,k2)
      implicit none
      integer k1,k2,i,j
      character line*(*)

      j=len(line)
      do i=1,j
         if (line(i:i).ne.' ') then
            k1=i; exit
         endif
         if (i==j) then
            k1=1; k2=1
            return
         endif
      enddo

      do i=j,1,-1
         if (line(i:i).ne.' ') then
            k2=i; exit
         endif
      enddo

      end

c     ##############################################################
c     ##  COPYRIGHT (C) 2005 by ITCC @ NJU, All Rights Reserved   ##
c     ##  subroutine NJ_sysinfo  --  show host os & directory     ##
c     ##  2005.12.26 by Wei Li; Update 2005.12.26 by Wei Li       ##
c     ##############################################################
c
      subroutine NJ_mypath(mypath,len_mypath,tmpname)
      implicit none
      integer io,system,i,j,k,itmp,len_mypath
      logical lg
      character line*100,line2*100,line3*100,mypath*(*),tmpname*100
C
      itmp=56
100   inquire(unit=itmp,opened=lg)
      if (lg) then
         itmp=itmp+1
         goto 100
      endif
C
C --- 01. Hostname
      i=system('hostname > '//tmpname)
      open(itmp,file = tmpname)
      read(itmp,'(a)') line
      do i=100,1,-1
         if (line(i:i).ne.' ') exit
      enddo
      close(itmp)
C --- 03. Directory
      k=system('pwd > '//tmpname)
      open(itmp,file=tmpname)
      read(itmp,'(a)') line3
      do k=100,1,-1
         if (line3(k:k).ne.' ') exit
      enddo
C
      close(itmp,status='delete')
C --- Output
      mypath=line(1:i)//':'//line3(1:k)
      len_mypath=i+k+1
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

c     ##############################################################
c     ##  COPYRIGHT (C) 2005 by ITCC @ NJU, All Rights Reserved   ##
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

c     ##############################################################
c     ##  COPYRIGHT (C) 2005 by ITCC @ NJU, All Rights Reserved   ##
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

c     ##############################################################
c     ##  COPYRIGHT (C) 2005 by ITCC @ NJU, All Rights Reserved   ##
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
C
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

C --- End the program 2007.01.29 by Wei Li ---
      subroutine mpijob_end(io,Tim0,Wall0)
      real*8 Tim0
      integer Wall0
      call NJ_date(io,'Task over on:')
      call NJ_cputim(io,Tim0)
      call NJ_walltim(io,Wall0)
      write(io,*) 'Normal Termination of MPIJOB Program!'
      end

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
      integer TIME,Wall,Wall0  !2006.03.06 Add (kind=4); Some other program with walltime may be modified

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

C --- mpijob title 2007.01.29 by Wei Li ---
      subroutine mpijob_title(io)
C
!     write(io,*) 'H==================================================H'
!     write(io,*) 'H  Lower Scaling Quantum Chemistry (LSQC) Program  H'
!     write(io,*) 'H  Copyright (C) 2006 - 2007 Nanjing Univ., China  H'
!     write(io,*) 'H  itcc.nju.edu.cn/lsqc  Email: shuhua@nju.edu.cn  H'
!     write(io,*) 'H==================================================H'
!
!     write(io,*) 'S. Li, W. Li, T. Fang, J. Ma, and Y. Jiang, '
!    &          //'LSQC Program,'
!     write(io,*) 'Version 1.1, Nanjing University, '
!    &          //'Nanjing, 2006.'
!     write(io,*)
C
      write(io,*) '          -------------------------------'
      write(io,*) '          M P I J O B -- 12 February 2007'
      write(io,*) '          -------------------------------'
      write(io,*) 'Submit jobs in .job with an MPI parallel'
C
      if (io>0) then
!        call NJ_date(io,'Task begin from:')
!        call NJ_sysinfo(io)
         write(io,*)
         return
      endif
C
      write(io,*)
      write(io,*) 'Input file: exam.job (defined jobs file)'
      write(io,*) 'Output files: exam.mpi (general output file)'
      write(io,*) '              exam.log (merged job outfiles)'
      write(io,*) 'Usage: mpirun -np k -machinename machinefile '
     &          //'./mpijob exam.job [-item value]'
      write(io,*) '  -append: append output file to an old .mpi'
      write(io,*) '  -np k: MPI parallel with k-1 CPUs'
      write(io,*) '  machinefile: hosts list file'
      write(io,*) '  -nfs: network file system (NFS) supported'
      write(io,*) 'Example for exam.job:'
      write(io,*) '  (Comments after ''!'')'
      write(io,*) '-----------------------------------------------'
      write(io,*) '3                     ! number of jobs'
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
      write(io,*) 'Example for machinefile:'
      write(io,*) '  (Comments after ''!'')'
      write(io,*) '-----------------------------------------------'
      write(io,*) 'node0                 ! main node'
      write(io,*) 'node0                 ! if main node is also used'
      write(io,*) 'node1'
      write(io,*) 'node2'
      write(io,*) 'node3'
      write(io,*) '-----------------------------------------------'
      write(io,*) 'Notes: (IMPORTANT!)'
      write(io,*) '  You should make the same directory in all nodes'
      write(io,*) 'and copy mpijob to each nodes including main node'
      write(io,*)
C
      if (io==0) return
c
      end

