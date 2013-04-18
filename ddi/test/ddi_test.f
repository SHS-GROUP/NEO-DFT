      program ddi_test
      implicit none
      integer  ddi_np,ddi_me,memrepl,memdist,ddi_root
      integer  ddi_nn,node_rank,node_root,smp_np,smp_me
      integer  handle,nrows,ncols,ilo,ihi,jlo,jhi
      integer  i,icorr,j,dlb_counter
      integer  len,to,from,msgid
      integer  na,ma,arrhand,indmin(2),indmax(2)
      integer  ngroups,comm(0:2),ddi_me_glob,ngroups0,mygroup,
     *         ddi_world,ddi_group,ddi_masters,mannod(0:1),
     *         dlb_counter_glob,master_id(0:1),ibuff(100)
      integer  itemp,lensm,sm_handle,locsm
      character*1 addr(1)
      logical extra
      parameter (na=6,ma=10)
      double precision arrchk(na,ma),vmin,vmax,dots(na),arrvec(na)
      double precision x,xcorr,buff(100),exetyp
      double precision temp,valsum,expected
      parameter(ddi_world=0,ddi_group=1,ddi_masters=2)
      data exetyp/8H        /
c
c     ----------  DDI library test program  ----------
c     the best test run would be two processes in each of two nodes.
c     the matrices here are pretty small, so don't try lots of cores.
c
      call ddi_init()
      call ddi_nproc(ddi_np,ddi_me)
      if (ddi_me.eq.0) write(6,3) ddi_np
c         set -extra- to true to see additional output from
c         other compute processes.  this can help with being
c         sure process counters are set right, and DLB is OK.
c         It does lead to a certain amount of output scrambling...
      extra = .false.
      if(ddi_me.eq.0 .or. extra) write(6,4) ddi_me
c
      call ddi_nnode(ddi_nn,node_rank)
      if (ddi_me.eq.0) write(6,9) ddi_nn
c
      call ddi_smp_nproc(smp_np,smp_me)
      if (ddi_me.eq.0  .or.  (extra .and. smp_me.eq.0)) then
         write(6,10) node_rank
         write(6,11) smp_np
      end if
      if(ddi_me.eq.0  .or.  extra) write(6,12) ddi_me,smp_me
c
c     enter the required memory, in megawords
c
      memrepl = 1
      memdist = 1
      call ddi_memory(memrepl,memdist,exetyp)
      if (ddi_me.eq.0) write(6,5) memdist
c
c     1. test global sum
c
      if(ddi_me.eq.0) write(6,*) ' '
      if(ddi_me.eq.0) write(6,*) 'testing all-to-all functionality...'
      if(ddi_me.eq.0) write(6,*) 'test global sums:'
      xcorr = (ddi_np*ddi_np + ddi_np)/2.0d+00
      x = ddi_me+1
      call ddi_gsumf(1,x,1)
      if(ddi_me.eq.0) write(6,2)
     *        ddi_me,'floating point global sum',x,xcorr
      icorr = (ddi_np*ddi_np + ddi_np)/2
      i = ddi_me+1
      call ddi_gsumi(2,i,1)
      if(ddi_me.eq.0) write(6,1) ddi_me,'integer global sum',i,icorr
c
c     2. test broadcast
c
      call ddi_sync(5)
      if (ddi_me.eq.0) write(6,*) 'test broadcasts:'
c
c     choose broadcast root process here
c
      ddi_root = 0
      if (ddi_me.eq.0) write(6,7) ddi_root
c
      x = 1.0d99
      if (ddi_me.eq.ddi_root) x = 3.141d+00
      call ddi_bcast(3,'f',x,1,ddi_root)
      if(ddi_me.eq.0) write(6,2)
     *       ddi_me,'floating point broadcast',x,3.141d+00
      i = 99
      if (ddi_me.eq.ddi_root) i = 3141
      call ddi_bcast(4,'i',i,1,ddi_root)
      if(ddi_me.eq.0) write(6,1) ddi_me,'integer broadcast',i,3141
c
c     3. test distributed data tools
c
      call ddi_sync(6)
      if(ddi_me.eq.0) write(6,*) ' '
      if(ddi_me.eq.0) write(6,*) 'testing distributed data...'
      nrows = 10
      ncols = 10
      call ddi_create(nrows,ncols,handle)
      call ddi_distrib(handle,ddi_me,ilo,ihi,jlo,jhi)
      call ddi_ndistrib(handle,node_rank,ilo,ihi,jlo,jhi)
      if (ddi_me.eq.0  .or.  (extra .and. smp_me.eq.0))
     *       write(6,13) node_rank,jlo,jhi,handle
c
c     fill new array with numbers
c
      do j = jlo, jhi
        do i = ilo, ihi
          buff(i) = (i+j+ddi_me)*1.77245d+00
        end do
        call ddi_put(handle,ilo,ihi,j,j,buff)
      end do
c
c     4. local get,acc
c
      call ddi_sync(7)
      if(ddi_me.eq.0) write(6,*) 'test local get,acc:'
      j = jlo
      call ddi_get(handle,ilo,ihi,j,j,buff)
      if(ddi_me.eq.0) write(6,6) ddi_me,'local get',(buff(i),i=ilo,ihi)
      call ddi_sync(8)
      do i = ilo, ihi
        buff(i) = 3.1415927d+00
      end do
      call ddi_acc(handle,ilo,ihi,j,j,buff)
      call ddi_sync(9)
      call ddi_get(handle,ilo,ihi,j,j,buff)
      if(ddi_me.eq.0) write(6,6)
     *      ddi_me,'local get after acc',(buff(i),i=ilo,ihi)
c
c     5. non-local get,acc
c
      call ddi_sync(10)
      if (ddi_me.eq.0) write(6,*) 'test remote get,acc:'
      j = jhi+2
      if (j.gt.ncols) j = j - ncols 
      call ddi_get(handle,ilo,ihi,j,j,buff)
      if(ddi_me.eq.0) write(6,6)
     *      ddi_me,'remote get',(buff(i),i=ilo,ihi)
      call ddi_sync(11)
      do i = ilo, ihi
        buff(i) = 9.8696d+00
      end do
      call ddi_acc(handle,ilo,ihi,j,j,buff)
      call ddi_sync(12)
      call ddi_get(handle,ilo,ihi,j,j,buff)
      if(ddi_me.eq.0) write(6,6)
     *      ddi_me,'remote get after acc',(buff(i),i=ilo,ihi)
c
      if(ddi_me.eq.0) write(6,27) handle
      call ddi_destroy(handle)
c
c     6. test DLB, both global and 'proc' versions
c
      if(ddi_me.eq.0) write(6,*) ' '
      if(ddi_me.eq.0) write(6,*) 'testing dynamic load balance...'
      dlb_counter = -1
      call ddi_dlbreset
      call ddi_dlbnext(dlb_counter)
      if(ddi_me.eq.0  .or.  extra)
     *   write(6,8) ddi_me,'got dlb task',dlb_counter
c          sum to see if we handed out counters 0, 1, ..., ddi_np-1
      call ddi_gsumi(9801,dlb_counter,1)
      icorr = ((ddi_np-1)*ddi_np)/2
      if(ddi_me.eq.0  .or.  extra) write(6,23) dlb_counter,icorr
c
c        next routine's function is not as well understood as desired!
c        you will note that its documentation is also fairly lame.
c
      if(ddi_me.eq.0) write(6,*) 'process level DLB test:'
      call ddi_procdlb_create(handle)
      call ddi_procdlb_reset(handle)
c        calling it repeatedly increments the counter up to ddi_me
      do i=0,ddi_me
         call ddi_procdlb_next(handle,ddi_me,dlb_counter)
      enddo
      call ddi_procdlb_destroy(handle) 
c        It seems each element of the distributed counter array is a
c        separate counter, so if incremented only by its matching
c        compute process, it is not really 'dynamic'...an answer
c        can therefore be predicted for the sum part of this test.
      call ddi_gsumi(9802,dlb_counter,1)
      icorr = ((ddi_np-1)*ddi_np)/2
      if(ddi_me.eq.0  .or.  extra) write(6,24) dlb_counter,icorr
c
c     7. test point-to-point message-passing
c
      if(ddi_me.eq.0) write(6,*) ' '
      if(ddi_me.eq.0) write(6,*) 'beginning point to point tests...'
      if (ddi_np.gt.1) then
        call ddi_sync(13)
        if (ddi_me.eq.0) write(6,*) 'test send/recv:'
        len = 1*8
        buff(1) = 0.0d+00
        from = 0
        to =   1
c
c       note even if proc 0,1 reside on the same 
c       transfer will be via the network card 
c
        if (ddi_me.eq.from.or.ddi_me.eq.to) then
c
c         synchronous message transfer
c
          if (ddi_me.eq.from) then 
            buff(1) = 3.14159d+00
            call ddi_send(buff,len,to)
            write(6,14) ddi_me,to
          end if
          if (ddi_me.eq.to) then 
            call ddi_recv(buff,len,from)
            write(6,15) ddi_me,buff(1),from
          end if
c
c         asynchronous message transfer
c
          buff(1) = 0.0d+00
          if (ddi_me.eq.from) then 
            buff(1) = 1.0d+00/3.0d+00
            call ddi_isend(buff,len,to,msgid)
            write(6,16) ddi_me,to
c
c           a real application would do lots of work independent of 
c           the contents of the message buffer (buff) here !
c
            call ddi_wait(msgid)
          end if
          if (ddi_me.eq.to) then 
            call ddi_irecv(buff,len,from,msgid)
            write(6,17) ddi_me,from
c
c           a real application would do lots of work independent of 
c           the contents of the message buffer (buff) here !
c
            call ddi_wait(msgid)
            write(6,18) ddi_me,buff(1)
          end if
        end if
        call ddi_sync(14)
      end if
c
c     8. test array operations
c
c     here we fetch the entire na x ma=6x10 matrix to a local array,
c     in order to print it out, and thus see if each step works OK.
c     ten columns fits nicely into a screen width.
c
      if(ddi_me.eq.0) write(6,*) ' '
      if(ddi_me.eq.0) write(6,*) 'beginning array operation tests...'
      call ddi_create(na,ma,arrhand)
c
      call ddi_arr_fill(arrhand,1,na,1,ma,3.3d+00)
      if(ddi_me.eq.0) then
         write(6,*) 'fill matrix with 3.3:'
         call ddi_get(arrhand,1,na,1,ma,arrchk)
         call print2d(arrchk,na,ma)
      endif
c
      call ddi_arr_scale(arrhand,1,na,1,ma,2.0d+00)
      if(ddi_me.eq.0) then
         write(6,*) 'scale entire matrix by 2.0:'
         call ddi_get(arrhand,1,na,1,ma,arrchk)
         call print2d(arrchk,na,ma)
      endif
c
      call ddi_arr_zero(arrhand,1,2,1,6)
      call ddi_arr_fill(arrhand,3,5,5,10,3.0d+00)
      if(ddi_me.eq.0) then
         write(6,*) 'zero patch= rows 1-2, columns 1-6'
         write(6,*) 'fill patch= rows 3-5, columns 5-10 with 3.0'
         call ddi_get(arrhand,1,na,1,ma,arrchk)
         call print2d(arrchk,na,ma)
      endif
c
      call ddi_arr_scale(arrhand,1,5,3,9,2.0d+00)
      call ddi_arr_fill(arrhand,2,2,10,10,33.0d+00)
      if(ddi_me.eq.0) then
         write(6,*) '   scale patch= rows 1-5, columns 3-9 by 2.0'
         write(6,*) 'then set element row 2,   column 10  to 33.0'
         call ddi_get(arrhand,1,na,1,ma,arrchk)
         call print2d(arrchk,na,ma)
      endif
c
c          this min/max is OK, operating on the entire matrix
      call ddi_arr_min(arrhand,1,na,1,ma,vmin,indmin(1),indmin(2))
      call ddi_arr_max(arrhand,1,na,1,ma,vmax,indmax(1),indmax(2))
      if(ddi_me.eq.0) write(6,19) vmin,indmin,vmax,indmax
c
c          The next part does not function properly, except if p=1,
c          when the subpatch is stored on the local compute process.
c          Sometimes it hangs, sometimes it returns incorrect values:
c          hang for 3x2 (3 nodes, 2 c.p.), no hang if 1x1, 2x2, 3x1.
c          The former is quite serious, so we comment this bit out.
c          Note to the future:  This call fails even if moved
c          in front of the previous min/max of the entire matrix.
c---  if(ddi_me.eq.0) write(6,*)
c--- *      'zooming in on patch= rows 1-3, columns 7-10:'
c---  call ddi_arr_min(arrhand,1,3,7,10,vmin,indmin(1),indmin(2))
c---  call ddi_arr_max(arrhand,1,3,7,10,vmax,indmax(1),indmax(2))
c---  if(ddi_me.eq.0) write(6,19) vmin,indmin,vmax,indmax
c
c          this is more general than DAXPY (Y=a*X+Y) since:
c          these are patches instead of vectors,
c          both patches can have a factor,
c          and the answer can go into a third patch.
c          And, it could be three different distributed matrices.
      call ddi_arr_add(arrhand,1,2,1,10,2.0d+00,
     *                 arrhand,3,4,1,10,1.0d+00,
     *                 arrhand,5,6,1,10)
      if(ddi_me.eq.0) then
         write(6,*) 'generalized DAXPY:'
         write(6,*) '  take 2.0 times patch at rows 1-2, cols 1-10'
         write(6,*) '  plus 1.0 times patch at rows 3-4, cols 1-10'
         write(6,*) '             and store at rows 5-6, cols 1-10.'
         call ddi_get(arrhand,1,na,1,ma,arrchk)
         call print2d(arrchk,na,ma)
      endif
c
c         the dot product operation works on patches in a single row,
c         but can do multiple rows, returning a vector in -dots-.
      call ddi_arr_dot(arrhand,1,6,1,10,
     *                 arrhand,1,6,1,10,dots)
      if(ddi_me.eq.0) then
         do i=1,na
            write(6,20) i,dots(i)
         enddo
         write(6,21)
      end if
c     
c          a generalized accumulate, which has a multiplicative factor,
c          and can operate on multiple rows if asked to.
c          Only one process should do the accumulation, of course
c
      do i=1,na
         arrvec(i) = i + 0.1d+00 * i
      enddo
c
      if(ddi_me.eq.0) call ddi_arr_acc(arrhand,2,2,2,7,2.0d+00,arrvec)
c
      if(ddi_me.eq.0) then
         write(6,22) arrvec
         call ddi_get(arrhand,1,na,1,ma,arrchk)
         call print2d(arrchk,na,ma)
      endif
c
c         this is the end of the array tests
c
      call ddi_sync(8001)
      call ddi_destroy(arrhand)
      if(ddi_me.eq.0) write(6,*) 'done with array operation tests'
c
c     9. test group operations
c        this test accomplishes something even if there is only 1 node,
c        but the test is more satisfactory when more nodes exist.
c
c         G1: test groups of 2 nodes (typical production)
c
      if(ddi_me.eq.0) write(6,*) ' '
      if(ddi_me.eq.0) write(6,*) 'beginning group DDI tests...'
      ngroups=(ddi_nn-1)/2+1
      if(ddi_nn.eq.2) ngroups=2
      ddi_me_glob=ddi_me
      call ddi_group_create(ngroups,comm(0),comm(1),comm(2))
      call ddi_scope(comm(ddi_group))
      CALL ddi_ngroup(ngroups0,mygroup)
      call ddi_nproc(ddi_np,ddi_me)
      if(ddi_me_glob.eq.0) then
         write(6,100) ngroups0,comm(0),comm(1),comm(2) 
         write(6,101) mygroup,ddi_np,ddi_me
      endif
      i=ddi_me+1
      call ddi_gsumi(15,i,1)
      if(ddi_me_glob.eq.0) then
         write(6,102) mygroup,i,(ddi_np*ddi_np+ddi_np)/2
      endif
c
c         G2: test distributed matrices with groups 
c
      nrows = 10
      ncols = 10
      call ddi_create(nrows,ncols,handle)
      call ddi_distrib(handle,ddi_me,ilo,ihi,jlo,jhi)
      call ddi_ndistrib(handle,node_rank,ilo,ihi,jlo,jhi)
      if (ddi_me_glob.eq.0) write(6,13) node_rank,jlo,jhi,handle
      do j = jlo, jhi
        do i = ilo, ihi
          buff(i) = (i+j+ddi_me)*1.77245d+00
        end do
      end do
      call ddi_put(handle,ilo,ihi,jlo,jhi,buff)
      call ddi_acc(handle,ilo,ihi,jlo,jhi,buff)
      call ddi_get(handle,ilo,ihi,jlo,jhi,buff)
      call ddi_destroy(handle)
c
c         G3: test manual group division
c
      call ddi_scope(comm(ddi_world))
      call ddi_gdlbreset()
      ngroups=min(2,ddi_nn)
      mannod(0)=max(ddi_nn-1,1)
      mannod(1)=min(ddi_nn-1,1)
      call ddi_group_create_custom(ngroups,mannod,
     *                             comm(0),comm(1),comm(2))
      call ddi_scope(comm(ddi_group))
      CALL ddi_ngroup(ngroups0,mygroup)
      call ddi_nproc(ddi_np,ddi_me)
      if(ddi_me_glob.eq.0) then
         write(6,103) ngroups0,comm(0),comm(1),comm(2) 
         write(6,101) mygroup,ddi_np,ddi_me
      endif
      x = ddi_me+1
      call ddi_gsumf(16,x,1)
      if(ddi_me_glob.eq.0) then
         write(6,104) mygroup,x,(ddi_np*ddi_np+ddi_np)/2.0d+00
      endif
      x = ddi_np 
      call ddi_bcast(17,'f',x,1,ddi_np-1)
      if(ddi_me.ne.0) x = 0
c
c         G4: test global counters at two levels (inter- and intra-group).
c             test master scope.
c
      dlb_counter_glob = -1
      call ddi_gdlbnext(dlb_counter_glob)
      dlb_counter = -1
      call ddi_dlbreset
      call ddi_dlbnext(dlb_counter)
      if(ddi_me_glob.eq.0) then
         write(6,105) mygroup,dlb_counter_glob,dlb_counter
      endif
      call ddi_scope(comm(ddi_masters))
      if(ngroups0.eq.2) master_id(1-mygroup)=0
      master_id(mygroup)=ddi_me_glob
      call ddi_gsumi(18,master_id,2)
      if(ddi_me_glob.eq.0) then
         if(ngroups0.eq.1) write(6,106) master_id(0)
         if(ngroups0.gt.1) write(6,107) (master_id(i),i=0,ngroups0-1)
      endif
      call ddi_scope(comm(ddi_world))
c
c         G5: sanity check: world has not changed?
c         Sum x across groups (accumulated in G3).
c
      call ddi_nproc(ddi_np,ddi_me)
      if(ddi_me.eq.0) then
         write(6,108) ddi_np,ddi_me
      endif
      call ddi_gsumf(19,x,1)
      if(ddi_me.eq.0) then
         write(6,109) x,ddi_np*1.0d+00
      endif
      if(ddi_me.eq.0) write(6,*) 'done with group DDI tests'
c
c     10. test node-based routines, which come in two basic types.
c     10a. collective routines inside the node, and between node masters
c
      if(ddi_me.eq.0) write(6,*) ' '
      if(ddi_me.eq.0) write(6,*) 'beginning SMP tests...'
      call ddi_smp_nproc(smp_np,smp_me)
      node_root = 0
      if(ddi_me.eq.ddi_root) then
         temp=1.0d+00
         itemp=1
      else
         temp=0.0d+00
         itemp=0
      endif
c         first, send 'one' from global master to node masters,
      if(smp_me.eq.0) call ddi_masters_bcast(1011,'F', temp,1,ddi_root)
      if(smp_me.eq.0) call ddi_masters_bcast(1012,'I',itemp,1,ddi_root)
c         send 'one' from node master to other processes inside node,
c         every process should now be holding a 'one'.
      call ddi_smp_bcast(1013,'F', temp,1,node_root)
      call ddi_smp_bcast(1014,'I',itemp,1,node_root)
c         sum the 'one' inside node, should total this node's "smp_np"
      call ddi_smp_gsumf(1015, temp,1)
      call ddi_smp_gsumi(1016,itemp,1)
c         sum each node's processor count across all nodes,
c         which should be the total number of compute processes.
      if(smp_me.eq.0) call ddi_masters_gsumf(1017, temp,1)
      if(smp_me.eq.0) call ddi_masters_gsumi(1018,itemp,1)
      if(ddi_me.eq.0) write(6,201) temp,itemp,ddi_np
c
c     10b. shared memory routines  (node-replicated memory)
c          create/use/destroy sequence
c
      lensm = 1000
      call ddi_smp_create(lensm,sm_handle)
      call ddi_smp_offset(sm_handle,addr,locsm)
      locsm = locsm + 1
      call ddi_smp_sync()
c
      call smp_test_routine(addr(locsm),lensm,valsum)
c
      call ddi_smp_destroy(sm_handle)
c
      if(smp_me.eq.0) call ddi_masters_gsumf(1019,valsum,1)
      valsum = valsum/ddi_nn
      expected = (lensm*lensm+lensm)/2.0d+00
      if(ddi_me.eq.0) write(6,202) lensm,valsum,expected
      if(ddi_me.eq.0) write(6,*) 'done with SMP tests.'
c
c     11. test scatter_accumulate
c         create/use/destroy sequence
c         (currently limited to row vector)
c
      call ddi_create(1,10,handle)
      j = 0
      do i = 1,10
         if(mod(i,ddi_np).eq.ddi_me)then
         j = j+1
         ibuff(j) = i

         buff(j) = i*1.0d+00
         endif
      enddo
      call ddi_scatter_acc(handle,j,ibuff,buff,1.0d+00)
      call ddi_sync(911)
      if(ddi_me.eq.0)call ddi_get(handle,1,1,1,10,buff)
      if(ddi_me.eq.0)write(6,203)buff(1:10)
      call ddi_destroy(handle)
c
      if(ddi_me.eq.0) write(6,*) ' '
      call ddi_finalize()
      stop
c
    1 format(2x,'cp=',1i5,2x,1a30,2x,1i8,' expected:',1i8)
    2 format(2x,'cp=',1i5,2x,1a30,2x,1f8.2,' expected:',1f8.2)
    3 format(2x,'there are ',1i5,' total compute processes (cp)')
    4 format(2x,'I am compute process',1i5)
    9 format(2x,'there are ',1i5,' nodes ')
   10 format(2x,'I am node',1i5)
   11 format(2x,'there are ',1i5,' local compute processes ')
   12 format(2x,'overall process',1i5,
     *     ' has local node compute process rank',1i5)
    5 format(2x,1i5,' megawords of memory is distributed')
    6 format(2x,'cp=',1i5,2x,1a30,':'/(2x,10f7.2))
   27 format(2x,'destroying distributed matrix [',i1,'].')
    7 format(2x,'broadcast root process is',1i5)
    8 format(2x,'cp=',1i5,2x,1a30,2x,1i8)
   23 format(2x,'sum of  global DLB counters=',i8,
     *          ', test passes if =',i8)
   24 format(2x,'sum of process DLB counters=',i8,
     *          ', test passes if =',i8)
   13 format(2x,'Node ',1i4,' stores columns ',1i5,' to ',1i5,
     *' of distributed matrix ',1i2)
   14 format(2x,'proc ',1i4,' sent message to ',1i4)
   15 format(2x,'proc ',1i4,' received ',1f8.5,' from ',1i4,
     *       ' (should be pi)')
   16 format(2x,'proc ',1i4,' posted isend to ',1i4)
   17 format(2x,'proc ',1i4,' posted irecv from ',1i4)
   18 format(2x,'proc ',1i4,' irecv done, got ',1f8.5,
     *       ' (should be 1/3)')
   19 format(2x,'smallest value=',f7.3,' at element',2i3/
     *       2x,' biggest value=',f7.3,' at element',2i3)
   20 format(2x,'dot product of row',i2,', cols 1-10 with itself=',f9.2)
   21 format(2x,'the expected answers for the dot products are'/
     *       2x,'566.28, 1611.72, 624.60, 624.60, 3919.32, 8417.88')
   22 format(2x,'generalized accumulate:'/
     *       2x,' add 2.0 times vector=',6f6.2/
     *       2x,' to rows 2-2, at cols 2-7.  Row 2, Col 7 should= 26.4')
  100 format(2x,'Created',i8,' groups, communicators are:',3i8)
  101 format(2x,'Group',i8,' begs to report:',i8,' nodes ready,',i8,
     *       ' rank speaking.')  
  102 format(2x,'Group',i8,' local summed to',i12,
     *                    ', correct answer:',i12)
  103 format(2x,'Manually created',i8,' groups, communicators are:',3i8)
  104 format('2x,Group',i8,' local summed to',f14.1,', correct answer:',
     *       f14.1)
  105 format(2x,'Group',i8,' got global counter:',i8,
     *          ' and local counter',i8)
  106 format(2x,'Master global rank is:',i8)
  107 format(2x,'Master global ranks are:',2i8)
  108 format(2x,'World scope has',i8,' nodes and my rank is',i8)  
  109 format(2x,'Summed to ',f10.1,', correct answer:',f10.1)  
  201 format(2x,'SMP bcast/gsum tests: for f.p.=',f10.1,', for int=',i8/
     *       2x,'test passes if both values above match compute',
     *          ' process count=',i8)
  202 format(2x,'SMP shared memory region test:'/
     *       2x,'sum of first',i8,' integers inside each node=',f10.1/
     *       2x,'test passes if the value above is ',f10.1)
  203 format(/'The following tests the scatter accumulate function: '/
     *    2x,'the vector increments by 1.00 from 1.00 to 10.00'/,10f7.2)

      end
c
      subroutine print2d(a,n,m)
      implicit double precision(a-h,o-z)
      dimension a(n,m)
      do i=1,n
         write(6,900) (a(i,j),j=1,m)
      enddo
      return
  900 format(10f7.3)
      end
c
      subroutine smp_test_routine(shared,lensm,valsum)
      implicit double precision(a-h,o-z)
      dimension shared(lensm)
      integer smp_np,smp_me
c
      call ddi_smp_nproc(smp_np,smp_me)
      ipcount = smp_me - 1
c
c        set entire shared memory region: shared(i) = i,
c        with only one core in each node initializing each element.
c
      do i=1,lensm
         ipcount=ipcount+1
         if(mod(ipcount,smp_np).eq.0) shared(i) = i
      enddo
c
      call ddi_smp_sync()
c
c        the node's master returns the sum of the first i integers,
c        the caller does nothing with other process' return values.
c
      if(smp_me.eq.0) then
         valsum = 0.0d+00
         do i=1,lensm
            valsum = valsum + shared(i)
         enddo
      else
         valsum = -999.0d+00
      end if
      return
      end
