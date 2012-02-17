c  5 May 98 - JRS - TINKER library routines that start with
c                   a to d Lib(rary)t(inker)ad 
c
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine active  --  set the list of active atoms  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "active" sets the list of atoms that are used during
c     each potential energy function calculation
c
c
      subroutine active
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'usage.i'
      integer i,j,next,center
      integer nfree,nfixed,nsphere
      integer nlist,list(maxatm)
      integer free(maxatm),fixed(maxatm)
      real*8 xcenter,ycenter,zcenter
      real*8 radius,radius2,dist2
      character*20 keyword
      character*80 record,string
      INTEGER me,master,nproc,ibtyp,iptim,iflag
      CHARACTER*8 GRPNAM
C
C
      LOGICAL LOG,FOUND,GOPARR,DSKWRK,MASWRK,TDSKWRK
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     set defaults for the numbers and lists of active atoms
c
      nuse = n
      do i = 1, n
         use(i) = .true.
      end do
      nfree = 0
      nfixed = 0
      do i = 1, n
         free(i) = 0
         fixed(i) = 0
      end do
      nsphere = 0
c
c     get any keywords containing active atom parameters
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
c
c     get any lists of atoms whose coordinates are active
c
         if (keyword(1:7) .eq. 'ACTIVE ') then
            string = record(next:80)
            read (string,*,err=10,end=10)  (free(i),i=nfree+1,n)
   10       continue
            dowhile (free(nfree+1) .ne. 0)
               nfree = nfree + 1
            end do
c
c     get any lists of atoms whose coordinates are inactive
c
         else if (keyword(1:9) .eq. 'INACTIVE ') then
            string = record(next:80)
            read (string,*,err=20,end=20)  (fixed(i),i=nfixed+1,n)
   20       continue
            dowhile (fixed(nfixed+1) .ne. 0)
               nfixed = nfixed + 1
            end do
c
c     get the center and radius of the sphere of active atoms
c
         else if (keyword(1:7) .eq. 'SPHERE ') then
            center = 0
            xcenter = 0.0d0
            ycenter = 0.0d0
            zcenter = 0.0d0
            radius = 0.0d0
            string = record(next:80)
            read (string,*,err=30,end=30)  xcenter,ycenter,
     &                                     zcenter,radius
   30       continue
            if (radius .eq. 0.0d0) then
               read (string,*,err=60,end=60)  center,radius
               xcenter = x(center)
               ycenter = y(center)
               zcenter = z(center)
            end if
            nsphere = nsphere + 1
            if (nsphere .eq. 1) then
               nuse = 0
               do i = 1, n
                  use(i) = .false.
               end do
               if (verbose) then
                  if (maswrk) write (iout,40)
   40             format (/,' Active Site Spheres used to',
     &                        ' Select Active Atoms :',
     &                     //,3x,'Atom Center',11x,'Coordinates',
     &                        12x,'Radius',6x,'# Active Atoms')
               end if
            end if
            radius2 = radius * radius
            do i = 1, n
               if (.not. use(i)) then
                  dist2 = (x(i)-xcenter)**2 + (y(i)-ycenter)**2
     &                            + (z(i)-zcenter)**2
                  if (dist2 .le. radius2) then
                     nuse = nuse + 1
                     use(i) = .true.
                  end if
               end if
            end do
            if (verbose) then
               if (maswrk) write (iout,50)  center,xcenter,ycenter,
     &                          zcenter,radius,nuse
   50          format (2x,i8,6x,3f9.2,2x,f9.2,7x,i8)
            end if
   60       continue
         end if
      end do
c
c     set active atoms to those not on the inactive atom list
c
      i = 1
      dowhile (fixed(i) .ne. 0)
         if (fixed(i) .gt. 0) then
            use(fixed(i)) = .false.
            nuse = nuse - 1
            i = i + 1
         else
            do j = abs(fixed(i)), abs(fixed(i+1))
               use(j) = .false.
               nuse = nuse - 1
            end do
            i = i + 2
         end if
      end do
c
c     set active atoms to only those on the active atom list
c
      i = 1
      dowhile (free(i) .ne. 0)
         if (i .eq. 1) then
            nuse = 0
            do j = 1, n
               use(j) = .false.
            end do
         end if
         if (free(i) .gt. 0) then
            use(free(i)) = .true.
            nuse = nuse + 1
            i = i + 1
         else
            do j = abs(free(i)), abs(free(i+1))
               use(j) = .true.
               nuse = nuse + 1
            end do
            i = i + 2
         end if
      end do
c
c     output the final list of the active atoms
c
      if (debug .and. nuse.gt.0 .and. nuse.lt.n) then
         nlist = 0
         do i = 1, n
            if (use(i)) then
               nlist = nlist + 1
               list(nlist) = i
            end if
         end do
         if (maswrk) write (iout,70)
   70    format (/,' List of Active Atoms for Energy',
     &              ' Calculations :',/)
         if (maswrk) write (iout,80)  (list(i),i=1,nlist)
   80    format (3x,10i7)
      end if
c
c     output the final list of the inactive atoms
c
      if (maswrk) write(iout,100)
  100 format (/,' List of InActive Atoms for Energy',
     &              ' Calculations :')
      if (maswrk) write (iout,150)  (fixed(i),i=1,nfixed)
      if (maswrk) write (iout,*)
  150 format (3x,10i7)
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  function adjacent  --  atom adjacent to specified atom  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "adjacent" finds an atom connected to atom "i1" other than
c     atom "i2"; if no such atom exists, then the closest atom
c     in space is returned; used only by subroutine "makeint"
c
c     variables and parameters :
c
c     mode   whether "makeint" is in manual mode, automatic, etc.
c     more   returned true if there is more than one previously
c              defined atom other than "i2" which is directly
c              connected (adjacent) to atom "i1"
c     iz0    line number of the Z-matrix on which an atom is
c              defined, 0 if not yet defined
c     iz1    line number of the Z-matrix on which the atom used
c              defining the bond length to a given atom is defined
c
c
      function adjacent (i1,i2,mode,more,iz0,iz1)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'zclose.i'
      integer i,j,k,i1,i2,nc,adjacent,mode
      integer ic(maxval),iz0(0:maxatm),iz1(maxatm)
      real*8 dist,short
      logical more
c
c
c     get a list of eligible atoms bonded to the atom of interest
c
      nc = 0
      more = .false.
      do j = 1, n12(i1)
         i = i12(j,i1)
         if (iz0(i).ne.0 .and. i.ne.i2) then
            if (i2 .eq. 0) then
               nc = nc + 1
               ic(nc) = i
            else
               if (iz1(i).eq.i1 .or. iz1(i1).eq.i) then
                  nc = nc + 1
                  ic(nc) = i
               end if
            end if
         end if
      end do
      if (nc .gt. 1)  more = .true.
c
c     if no bonded atom is eligible, use the nearest neighbor
c
      if (nc .eq. 0) then
         adjacent = 0
         if (mode .eq. 1) then
            write (iout,10)  i1
   10       format (' ADJACENT  --  Atom',i6,' not Attached',
     &                 ' to any Prior Atom')
         else
            short = 1000000.0d0
            do i = 1, n
               if (iz0(i).ne.0 .and. i.ne.i1 .and. i.ne.i2) then
                  dist = sqrt((x(i)-x(i1))**2 + (y(i)-y(i1))**2
     &                              + (z(i)-z(i1))**2)
                  if (dist .lt. short) then
                     short = dist
                     adjacent = i
                  end if
               end if
            end do
            if (i2 .eq. 0) then
               ndel = ndel + 1
               idel(1,ndel) = adjacent
               idel(2,ndel) = i1
               write (iout,10)  i1
            end if
         end if
c
c     for automatic mode, always use the first eligible bonded atom
c
      else if (mode .eq. 0) then
         adjacent = ic(1)
c
c     for torsion mode, use an adjacent atom bonded to undefined atoms
c
      else if (mode .eq. 3) then
         adjacent = ic(1)
         do k = 1, nc
            do j = 1, n12(ic(k))
               i = i12(j,ic(k))
               if (iz0(i).ne.0 .and. i.ne.i1) then
                  adjacent = ic(k)
                  goto 20
               end if
            end do
         end do
   20    continue
c
c     if only one directly bonded atom is eligible, then use it
c
      else if (nc .eq. 1) then
         adjacent = ic(1)
         write (iout,30)  ic(1)
   30    format (' Atom',i6,' is the only Connected Atom')
c
c     ask the user which eligible bonded atom to use as adjacent
c
      else
   40    continue
         if (nc .eq. 2) then
            write (iout,50)  (ic(j),j=1,nc)
   50       format (' Choose a Connected Atom (',2i6,') :  ',$)
         else if (nc .eq. 3) then
            write (iout,60)  (ic(j),j=1,nc)
   60       format (' Choose a Connected Atom (',3i6,') :  ',$)
         else if (nc .eq. 4) then
            write (iout,70)  (ic(j),j=1,nc)
   70       format (' Choose a Connected Atom (',4i6,') :  ',$)
         else if (nc .eq. 5) then
            write (iout,80)  (ic(j),j=1,nc)
   80       format (' Choose a Connected Atom (',5i6,') :  ',$)
         else if (nc .eq. 6) then
            write (iout,90)  (ic(j),j=1,nc)
   90       format (' Choose a Connected Atom (',6i6,') :  ',$)
         else if (nc .eq. 7) then
            write (iout,100)  (ic(j),j=1,nc)
  100       format (' Choose a Connected Atom (',7i6,') :  ',$)
         else if (nc .eq. 8) then
            write (iout,110)  (ic(j),j=1,nc)
  110       format (' Choose a Connected Atom (',8i6,') :  ',$)
         end if
         read (input,120,err=40)  adjacent
  120    format (i10)
         if (adjacent .eq. 0) then
            adjacent = ic(1)
         else
            do j = 1, nc
               if (ic(j) .eq. adjacent)  goto 130
            end do
            goto 40
  130       continue
         end if
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine analysis  --  energy components and analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "analysis" calls the subroutines to calculate the potential
c     energy and perform an energy partitioning analysis in terms
c     of type of interaction and/or atom number
c
c
      subroutine analysis (energy)
      implicit none
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'bound.i'
      include 'energi.i'
      include 'moment.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer i
      real*8 energy
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      e14 = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     zero out energy partitioning components for each atom
c
      do i = 1, n
         aeb(i) = 0.0d0
         aea(i) = 0.0d0
         aeba(i) = 0.0d0
         aeub(i) = 0.0d0
         aeaa(i) = 0.0d0
         aeopb(i) = 0.0d0
         aeid(i) = 0.0d0
         aeit(i) = 0.0d0
         aet(i) = 0.0d0
         aebt(i) = 0.0d0
         aett(i) = 0.0d0
         aev(i) = 0.0d0
         ae14(i) = 0.0d0
         aec(i) = 0.0d0
         aecd(i) = 0.0d0
         aed(i) = 0.0d0
         aem(i) = 0.0d0
         aep(i) = 0.0d0
         aer(i) = 0.0d0
         aes(i) = 0.0d0
         aeg(i) = 0.0d0
         aex(i) = 0.0d0
      end do
c
c     zero out the net charge and dipole moment components
c
      netchg = 0.0d0
      xdipole = 0.0d0
      ydipole = 0.0d0
      zdipole = 0.0d0
c
c     maintain any periodic boundary conditions
c
      if (use_bounds)  call bounds
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call piscf
c
c     call the local geometry energy component routines
c
      if (use_bond)  call ebond3
      if (use_angle)  call eangle3
      if (use_strbnd)  call estrbnd3
      if (use_urey)  call eurey3
      if (use_angang)  call eangang3
      if (use_opbend)  call eopbend3
      if (use_improp)  call eimprop3
      if (use_imptor)  call eimptor3
      if (use_tors)  call etors3
      if (use_strtor)  call estrtor3
c     if (use_tortor)  call etortor3
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj3
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck3
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb3
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal3
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss3
      end if
c
c     call the electrostatic energy component routines
c
      if (use_charge) then
         if (use_deform) then
            call echarge9
         else
            call echarge3
         end if
      end if
      if (use_chgdpl)  call echgdpl3
      if (use_dipole)  call edipole3
      if (use_mpole .or. use_polar)  call empole3
      if (use_rxnfld)  call erxnfld3
c
c     call any miscellaneous energy component routines
c
      if (use_solv)  call esolv3
      if (use_geom)  call egeom3
      if (use_extra)  call extra3
c
c     sum up to give the total potential energy
c
      energy = eb + ea + eba + eub + eaa + eopb + eid
     &            + eit + et + ebt + ett + ev + e14 + ec
     &            + ecd + ed + em + ep + er + es + eg + ex
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine angles  --  locate and store bond angles  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "angles" finds the total number of bond angles and stores
c     the atom numbers of the atoms defining each angle; for
c     each angle to a tricoordinate central atom, the third
c     bonded atom is stored for use in out-of-plane bending
c
c
      subroutine angles
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmlst.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      integer i,j,k,iangle
c
c
c     loop over all atoms, storing the atoms in each bond angle
c
      nangle = 0
      do i = 1, n
         iangle = 0
         do j = 2, n12(i)
            do k = 1, j-1
               nangle = nangle + 1
               if (nangle .gt. maxang) then
                  write (iout,10)
   10             format (/,' ANGLES  --  Too many Bond Angles;',
     &                       ' Increase MAXANG')
                  call fatal
               end if
               iangle = iangle + 1
               anglist(iangle,i) = nangle
               iang(1,nangle) = i12(k,i)
               iang(2,nangle) = i
               iang(3,nangle) = i12(j,i)
               iang(4,nangle) = 0
            end do
         end do
c
c     set the out-of-plane atom for angles at trivalent centers
c
         if (n12(i) .eq. 3) then
            iang(4,nangle) = i12(1,i)
            iang(4,nangle-1) = i12(2,i)
            iang(4,nangle-2) = i12(3,i)
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine attach  --  setup of connectivity arrays  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "attach" sets up the lists of 1-3 and 1-4 connectivities
c     starting from the previously determined list of attached
c     atoms (ie, 1-2 connectivity)
c
c
      subroutine attach
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      integer i,j,k,m,p
      integer jj,kk,mm
c
c
c     loop over all atoms finding all the 1-3 relationships;
c     note "n12" and "i12" have already been setup elsewhere
c
      do i = 1, n
         n13(i) = 0
         do j = 1, n12(i)
            jj = i12(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (kk .eq. i)  goto 10
               do m = 1, n12(i)
                  if (kk .eq. i12(m,i))  goto 10
               end do
               n13(i) = n13(i) + 1
               i13(n13(i),i) = kk
   10          continue
            end do
         end do
         call sort (n13(i),i13(1,i))
      end do
c
c     loop over all atoms finding all the 1-4 relationships
c
      do i = 1, n
         n14(i) = 0
         do j = 1, n12(i)
            jj = i12(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               do m = 1, n12(kk)
                  mm = i12(m,kk)
                  if (mm .eq. i)  goto 20
                  do p = 1, n12(i)
                     if (mm .eq. i12(p,i))  goto 20
                  end do
                  do p = 1, n13(i)
                     if (mm .eq. i13(p,i))  goto 20
                  end do
                  n14(i) = n14(i) + 1
                  i14(n14(i),i) = mm
   20             continue
               end do
            end do
         end do
         call sort (n14(i),i14(1,i))
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine basefile  --  get base prefix from a filename  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "basefile" extracts from an input filename the portion
c     consisting of any directory name and the base filename
c
c
      subroutine basefile (string)
      implicit none
      include 'files.i'
      integer i,k,trimtext
      character*1 letter
      character*60 string
c
      INTEGER me,master,nproc,ibtyp,iptim
      integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
c
c
c     store the input filename and find its full length
c
      filename = string
      leng = trimtext (string)
c
c     count the number of characters prior to any extension
c
      k = leng
      do i = 1, leng
         letter = string(i:i)
         if (letter .eq. '/')  k = leng
c        if (letter .eq. '\')  k = leng
         if (ichar(letter) .eq. 92)  k = leng
         if (letter .eq. ']')  k = leng
         if (letter .eq. ':')  k = leng
         if (letter .eq. '~')  k = leng
         if (letter .eq. '.')  k = i - 1
      end do
      leng = min(leng,k)
c
c     find the length of any directory name prefix
c
      k = 0
      do i = leng, 1, -1
         letter = string(i:i)
         if (letter .eq. '/')  k = i
c        if (letter .eq. '\')  k = i
         if (ichar(letter) .eq. 92)  k = i
         if (letter .eq. ']')  k = i
         if (letter .eq. ':')  k = i
         if (letter .eq. '~')  k = i
         if (letter .eq. '.')  k = i
         if (k .ne. 0)  goto 10
      end do
   10 continue
      ldir = k
c
c     read and store the keywords from the keyfile
c
      call getkey(ir,iw)
c
c     get the information level and output style
c
      call control
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine beeman  --  computes a molecular dynamics step  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "beeman" performs a single molecular dynamics time step
c     by means of a Beeman multistep recursion formula; the
c     actual coefficients in the recursions are the "Better
c     Beeman" values due to Bernie Brooks which equalize the
c     expected errors in the potential and kinetic energies
c
c
      subroutine beeman (istep,dt,dt_8,dt2_8,ndump)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,istep,ndump
      real*8 dt,dt_8,dt2_8,e_tot,e_kin,e_pot
      real*8 temp,pres,vol,derivs(3,maxatm)
      real*8 xterm,yterm,zterm
      real*8 x_old(maxatm),y_old(maxatm),z_old(maxatm)
c
c
c     store the current atom positions, then find new atom
c     positions and half-step velocities via Beeman recursion
c
      do i = 1, n
         if (use(i)) then
            x_old(i) = x(i)
            y_old(i) = y(i)
            z_old(i) = z(i)
            xterm = 5.0d0*a(1,i) - a_old(1,i)
            yterm = 5.0d0*a(2,i) - a_old(2,i)
            zterm = 5.0d0*a(3,i) - a_old(3,i)
            x(i) = x(i) + v(1,i)*dt + xterm*dt2_8
            y(i) = y(i) + v(2,i)*dt + yterm*dt2_8
            z(i) = z(i) + v(3,i)*dt + zterm*dt2_8
            v(1,i) = v(1,i) + xterm*dt_8
            v(2,i) = v(2,i) + yterm*dt_8
            v(3,i) = v(3,i) + zterm*dt_8
         end if
      end do
c
c     apply "rattle" to correct atom positions and velocities
c
      if (use_rattle)  call rattle (dt,x_old,y_old,z_old)
c
c     get the potential energy and atomic forces
c
      call gradient (e_pot,derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Beeman recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a_old(j,i) = a(j,i)
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + (3.0d0*a(j,i)+a_old(j,i))*dt_8
            end do
         end if
      end do
c
c     use "rattle" method to get corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     sum the total kinetic energy over all atoms
c
      e_kin = 0.0d0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               e_kin = e_kin + mass(i)*v(j,i)**2
            end do
         end if
      end do
      e_kin = 0.5d0 * e_kin / convert
c
c     determine system temperature and total energy
c
      temp = 2.0d0 * e_kin / (dble(3*nuse-nrat-6) * gasconst)
      e_tot = e_kin + e_pot
c
c     control temperature and pressure via external bath coupling
c
      if (isothermal)  call temper (dt,temp)
      if (isobaric)  call pressure (dt,pres,vol,e_kin)
c
c     compute any averages or statistics for this step
c
      call mdstat (istep,dt,e_tot,e_pot,e_kin,temp,pres,vol,ndump)
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  function bndangle  --  bond angle between three atoms  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "bndangle" finds the value of the bond angle defined
c     by three input atoms
c
c
      function bndangle (ia,ib,ic)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'math.i'
      integer ia,ib,ic
      real*8 bndangle
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rcb2
      real*8 dot,cosine
c
c
c     set default value in case atoms are degenerate
c
      bndangle = 0.0d0
c
c     compute the value in degrees of the bond angle
c
      xab = x(ia) - x(ib)
      yab = y(ia) - y(ib)
      zab = z(ia) - z(ib)
      rab2 = xab*xab + yab*yab + zab*zab
      xcb = x(ic) - x(ib)
      ycb = y(ic) - y(ib)
      zcb = z(ic) - z(ib)
      rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
      if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
         dot = xab*xcb + yab*ycb + zab*zcb
         cosine = dot / sqrt(rab2*rcb2)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         bndangle = radian * acos(cosine)
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  function bndleng  --  bond length between two atoms  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "bndleng" finds the value of the bond length defined
c     by two input atoms
c
c
      function bndleng (ia,ib)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer ia,ib
      real*8 bndleng
      real*8 xab,yab,zab
      real*8 rab2
c
c
c     set default value in case atoms are degenerate
c
      bndleng = 0.0d0
c
c     compute the value in Angstroms of the bond length
c
      xab = x(ia) - x(ib)
      yab = y(ia) - y(ib)
      zab = z(ia) - z(ib)
      rab2 = xab*xab + yab*yab + zab*zab
      if (rab2 .ne. 0.0d0) then
         bndleng = sqrt(rab2)
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine bonds  --  locate and store covalent bonds  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "bonds" finds the total number of covalent bonds and
c     stores the atom numbers of the atoms defining each bond
c
c
      subroutine bonds
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'iounit.i'
      integer i,j,k,m
c
c
c     loop over all atoms, storing the atoms in each bond
c
      nbond = 0
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            if (i .lt. k) then
               nbond = nbond + 1
               if (nbond .gt. maxbnd) then
                  write (iout,10)
   10             format (/,' BONDS  --  Too many Bonds; Increase',
     &                       ' MAXBND')
                  call fatal
               end if
               ibnd(1,nbond) = i
               ibnd(2,nbond) = k
               bndlist(j,i) = nbond
               do m = 1, n12(k)
                  if (i .eq. i12(m,k)) then
                     bndlist(m,k) = nbond
                     goto 20
                  end if
               end do
   20          continue
            end if
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine born  --  computation of Born radii for GB/SA  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "born" computes the Born radius of each atom for use with
c     the Macromodel GB/SA solvation model
c
c     literature references:
c
c     W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson,
c     "A Semianalytical Treatment of Solvation for Molecular
c     Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129
c     (1990)  (see supplimentary material for details)
c
c     G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Pairwise Solute
c     Descreening of Solute Charges from a Dielectric Medium", Chem.
c     Phys. Lett., 246, 122-129 (1995)  (approximate method)
c
c
      subroutine born
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'solute.i'
      integer i,k,ncalls,atn
      real*8 area,rold,t
      real*8 shell,fraction
      real*8 inner,outer
      real*8 offset,tinit,ratio
      real*8 total,r(maxatm)
      real*8 xi,yi,zi,dik
      real*8 ri,rk,si,si2,sk,sk2
      real*8 lik,lik2,uik,uik2,term
      real*8 sum(maxatm),fs(maxatm)
      logical done
      data ncalls  / 0 /
      save ncalls
c
c
c     increment the number of calls to this routine
c
      ncalls = ncalls + 1
      if (reborn.eq.0 .or. ncalls.lt.reborn)  return
      ncalls = 0
c
c     get the GB/SA dielectric offset atomic radii
c
      offset = -0.09d0
      do i = 1, n
         r(i) = rsolv(i) + offset
      end do
c
c     get the Born radii by Still's original method
c
      if (n .le. bornmax) then
         tinit = 0.1d0
         ratio = 1.5d0
         do i = 1, n
            t = tinit
            rold = r(i)
            total = 0.0d0
            done = .false.
            dowhile (.not. done)
               r(i) = r(i) + 0.5d0*t
               call surfatom (i,area,r)
               fraction = area / (4.0d0*pi*r(i)**2)
               if (fraction .lt. 0.99d0) then
                  inner = r(i) - 0.5d0*t
                  outer = inner + t
                  shell = 1.0d0/inner - 1.0d0/outer
                  total = total + fraction*shell
                  r(i) = r(i) + 0.5d0*t
                  t = ratio * t
               else
                  inner = r(i) - 0.5d0*t
                  total = total + 1.0d0/inner
                  done = .true.
               end if
            end do
            rborn(i) = 1.0d0 / total
            r(i) = rold
         end do
c
c     set the radii factors for the approximate pairwise method
c
      else
         do i = 1, n
            sum(i) = 1.0d0 / r(i)
            fs(i) = 0.80d0
            atn = atomic(i)
            if (atn .eq. 1)  fs(i) = 0.85d0
            if (atn .eq. 6)  fs(i) = 0.72d0
            if (atn .eq. 7)  fs(i) = 0.79d0
            if (atn .eq. 8)  fs(i) = 0.85d0
            if (atn .eq. 9)  fs(i) = 0.88d0
            if (atn .eq. 15)  fs(i) = 0.86d0
            if (atn .eq. 16)  fs(i) = 0.96d0
         end do
c
c     get the Born radii via the approximate pairwise method
c
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = r(i)
            si = ri * fs(i)
            si2 = si * si
            do k = i+1, n
               rk = r(k)
               sk = rk * fs(k)
               sk2 = sk * sk
               dik = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
               if (ri .lt. dik+sk) then
                  lik = 1.0d0 / max(ri,dik-sk)
                  uik = 1.0d0 / (dik+sk)
                  lik2 = lik * lik
                  uik2 = uik * uik
                  term = lik - uik + 0.25d0*dik*(uik2-lik2)
     &                      + (0.5d0/dik)*log(uik/lik)
     &                      + (0.25d0*sk2/dik)*(lik2-uik2)
                  sum(i) = sum(i) - 0.5d0*term
               end if
               if (rk .lt. dik+si) then
                  lik = 1.0d0 / max(rk,dik-si)
                  uik = 1.0d0 / (dik+si)
                  lik2 = lik * lik
                  uik2 = uik * uik
                  term = lik - uik + 0.25d0*dik*(uik2-lik2)
     &                      + (0.5d0/dik)*log(uik/lik)
     &                      + (0.25d0*si2/dik)*(lik2-uik2)
                  sum(k) = sum(k) - 0.5d0*term
               end if
            end do
         end do
         do i = 1, n
            rborn(i) = 1.0d0 / sum(i)
            rborn(i) = max(r(i),rborn(i))
         end do
      end if
c
c     write out the Born radius of each atom
c
      if (debug) then
         write (iout,10)
   10    format (/,' Born Radius for each Atom :',/)
         k = 1
         dowhile (k .le. n)
            write (iout,20)  (i,rborn(i),i=k,min(k+4,n))
   20       format (1x,5(i7,f8.3))
            k = k + 5
         end do
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bounds  --  check periodic boundary conditions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bounds" finds the center of mass of each molecule, translates
c     any stray molecules back into the periodic box, and saves the
c     offset of each atom relative to the molecular center of mass
c
c
      subroutine bounds
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'centre.i'
      include 'molcul.i'
      integer i,j,k,init,stop
      real*8 weight,x_mid,y_mid,z_mid
      real*8 x_frac,y_frac,z_frac
c
c
c     locate the center of mass of each molecule
c
      do i = 1, nmol
         init = imol(1,i)
         stop = imol(2,i)
         x_mid = 0.0d0
         y_mid = 0.0d0
         z_mid = 0.0d0
         do j = init, stop
            k = kmol(j)
            weight = mass(k)
            x_mid = x_mid + x(k)*weight
            y_mid = y_mid + y(k)*weight
            z_mid = z_mid + z(k)*weight
         end do
         weight = molmass(i)
         x_mid = x_mid / weight
         y_mid = y_mid / weight
         z_mid = z_mid / weight
c
c     store coordinates of atoms relative to center of mass
c
         do j = init, stop
            k = kmol(j)
            xcm(k) = x(k) - x_mid
            ycm(k) = y(k) - y_mid
            zcm(k) = z(k) - z_mid
         end do
c
c     get fractional coordinates of center of mass
c
         if (orthogonal .or. octahedron) then
            z_frac = z_mid
            y_frac = y_mid
            x_frac = x_mid
         else if (monoclinic) then
            z_frac = z_mid / beta_sin
            y_frac = y_mid
            x_frac = x_mid - z_frac*beta_cos
         else if (triclinic) then
            z_frac = z_mid / gamma_term
            y_frac = (y_mid - z_frac*beta_term) / gamma_sin
            x_frac = x_mid - y_frac*gamma_cos - z_frac*beta_cos
         end if
c
c     translate center of mass into the periodic box
c
         dowhile (x_frac .gt. xbox2)
            x_frac = x_frac - xbox
         end do
         dowhile (x_frac .lt. -xbox2)
            x_frac = x_frac + xbox
         end do
         dowhile (y_frac .gt. ybox2)
            y_frac = y_frac - ybox
         end do
         dowhile (y_frac .lt. -ybox2)
            y_frac = y_frac + ybox
         end do
         dowhile (z_frac .gt. zbox2)
            z_frac = z_frac - zbox
         end do
         dowhile (z_frac .lt. -zbox2)
            z_frac = z_frac + zbox
         end do
c
c     truncated octahedron needs to have corners removed
c
         if (octahedron) then
            if (abs(x_frac)+abs(y_frac)+abs(z_frac) .gt. box34) then
               x_frac = x_frac - sign(xbox2,x_frac)
               y_frac = y_frac - sign(ybox2,y_frac)
               z_frac = z_frac - sign(zbox2,z_frac)
            end if
         end if
c
c     convert fractional center of mass back to Cartesian
c
         if (orthogonal .or. octahedron) then
            x_mid = x_frac
            y_mid = y_frac
            z_mid = z_frac
         else if (monoclinic) then
            x_mid = x_frac + z_frac*beta_cos
            y_mid = y_frac
            z_mid = z_frac * beta_sin
         else if (triclinic) then
            x_mid = x_frac + y_frac*gamma_cos + z_frac*beta_cos
            y_mid = y_frac*gamma_sin + z_frac*beta_term
            z_mid = z_frac * gamma_term
         end if
c
c     translate atomic coordinates based upon center of mass
c
         do j = init, stop
            k = kmol(j)
            x(k) = xcm(k) + x_mid
            y(k) = ycm(k) + y_mid
            z(k) = zcm(k) + z_mid
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine calendar  --  find the current date and time  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "calendar" returns the current time as a set of integer
c     values representing the year, month, day, hour, minute
c     and second; only one of the various machine implementations
c     included should be activated by removing comment characters
c
c
      subroutine calendar (year,month,day,hour,minute,second)
      implicit none
      integer year,month,day
      integer hour,minute,second
c
c
c     Unix machines use the "itime" and "idate" intrinsics,
c     this code works for Sun, DEC, SGI and probably others;
c     also for Intel PC's under Digital Visual Fortran 5.0
c     and for Macintosh under Absoft Pro Fortran 5.0
c
      integer hms(3)
      call itime (hms)
      hour = hms(1)
      minute = hms(2)
      second = hms(3)
      call idate (month,day,year)
c
c     code for IBM RS/6000 under AIX appends an underscore
c
c     integer hms(3)
c     call itime_(hms)
c     hour = hms(1)
c     minute = hms(2)
c     second = hms(3)
c     call idate_(month,day,year)
c
c     code for Hewlett-Packard HP-UX systems and Digital OpenVMS;
c     for HP machines the +E1 compiler option is also required
c
c     real secnds
c     second = int(secnds(0.0))
c     hour = mod(second,3600)
c     second = second - 3600*hour
c     minute = mod(second,60)
c     second = second - 60*minute
c     call idate (month,day,year)
c
c     Intel-based PC's under the Microsoft and Watcom compilers
c
c     integer*2 hrs,mins,secs,hsecs
c     integer*2 yrs,mons,days
c     call gettim (hrs,mins,secs,hsecs)
c     hour = hrs
c     minute = mins
c     second = secs
c     call getdat (yrs,mons,days)
c     year = yrs
c     month = mons
c     day = days
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine center  --  superimpose structure centroids  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "center" moves the weighted centroid of each coordinate
c     set to the origin during least squares superposition
c
c
      subroutine center (n1,x1,y1,z1,n2,x2,y2,z2,xmid,ymid,zmid)
      implicit none
      include 'sizes.i'
      include 'align.i'
      integer i,k,n1,n2
      real*8 xmid,ymid,zmid,weight,norm
      real*8 x1(maxatm),y1(maxatm),z1(maxatm)
      real*8 x2(maxatm),y2(maxatm),z2(maxatm)
c
c
c     find the weighted centroid of the second
c     structure and translate it to the origin
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         k = ifit(2,i)
         weight = wfit(i)
         xmid = xmid + x2(k)*weight
         ymid = ymid + y2(k)*weight
         zmid = zmid + z2(k)*weight
         norm = norm + weight
      end do
      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n2
         x2(i) = x2(i) - xmid
         y2(i) = y2(i) - ymid
         z2(i) = z2(i) - zmid
      end do
c
c     now repeat for the first structure, note
c     that this centroid position gets returned
c
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         k = ifit(1,i)
         weight = wfit(i)
         xmid = xmid + x1(k)*weight
         ymid = ymid + y1(k)*weight
         zmid = zmid + z1(k)*weight
         norm = norm + weight
      end do
      xmid = xmid / norm
      ymid = ymid / norm
      zmid = zmid / norm
      do i = 1, n1
         x1(i) = x1(i) - xmid
         y1(i) = y1(i) - ymid
         z1(i) = z1(i) - zmid
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine cholesky  --  modified Cholesky linear solve  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "cholesky" uses a modified Cholesky method to solve the linear
c     system Ax = b, returning "x" in "b"; "A" is assumed to be a
c     real symmetric positive definite matrix with diagonal and upper
c     triangle stored by rows in "A"; thus the actual size of the
c     passed portion of "A" is nvar*(nvar+1)/2
c
c     literature reference:
c
c     R. S. Martin, G. Peters and J. H. Wilkinson, Numer. Math.,
c     7, 362-383 (1965)
c
c
      subroutine cholesky (nvar,a,b)
      implicit none
      integer nvar,i,j,k,ii,ij,ik,im,jk,jm,ki,kk
      real*8 a(nvar*(nvar+1)/2),b(nvar),r,s,t
c
c
c     use Cholesky to reduce "A" to (L)(D)(L transpose) = A;
c     "L" has a unit diagonal; store 1.0/D on the diagonal of "A"
c
      ii = 1
      do i = 1, nvar
         im = i - 1
         if (i .ne. 1) then
            ij = i
            do j = 1, im
               r = a(ij)
               if (j .ne. 1) then
                  ik = i
                  jk = j
                  jm = j - 1
                  do k = 1, jm
                     r = r - a(ik)*a(jk)
                     ik = nvar - k + ik
                     jk = nvar - k + jk
                  end do
               end if
               a(ij) = r
               ij = nvar - j + ij
            end do
         end if
         r = a(ii)
         if (i .ne. 1) then
            kk = 1
            ik = i
            do k = 1, im
               s = a(ik)
               t = s * a(kk)
               a(ik) = t
               r = r - s*t
               ik = nvar - k + ik
               kk = nvar - k + 1 + kk
            end do
         end if
         a(ii) = 1.0d0 / r
         ii = nvar - i + 1 + ii
      end do
c
c     solve linear equations; first solve Ly = b for y
c
      do i = 1, nvar
         if (i .ne. 1) then
            ik = i
            im = i - 1
            r = b(i)
            do k = 1, im
               r = r - b(k)*a(ik)
               ik = nvar - k + ik
            end do
            b(i) = r
         end if
      end do
c
c     now, solve (D)(L transpose)(x) = y for x
c
      ii = nvar*(nvar+1)/2
      do j = 1, nvar
         i = nvar + 1 - j
         r = b(i) * a(ii)
         if (j .ne. 1) then
            im = i + 1
            ki = ii + 1
            do k = im, nvar
               r = r - a(ki)*b(k)
               ki = ki + 1
            end do
         end if
         b(i) = r
         ii = ii - j - 1
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine clock  --  find elapsed time for current job  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "clock" determines elapsed CPU time in seconds since the
c     start of the job; only one of the implementations should
c     be activated by removing comment characters from the code
c
c
      subroutine clock (seconds)
      implicit none
      real*8 seconds
c
c
c     Unix machines have access to the "etime" intrinsic,
c     this code works for Sun, DEC, SGI and probably others;
c     also for Macintosh under Absoft Pro Fortran 5.0
c
      real etime,times(2)
      seconds = dble(etime(times))
c
c     IBM RS/6000 under AIX appends an underscore to "etime"
c
c     real etime_,times(2)
c     seconds = dble(etime_(times))
c
c     Hewlett-Packard machines under HP-UX only give wall clock
c     time; HP's also require use of the +E1 compiler option
c
c     real secnds,start
c     logical initial
c     save initial,start
c     data initial  / .true. /
c     if (initial) then
c        initial = .false.
c        start = secnds (0.0)
c     end if
c     seconds = secnds (start)
c
c     Digital OpenVMS machines use a special system call
c
c     include '($jpidef)'
c     structure /itmlst3_1item/
c        structure item
c           integer*2 buffer_length
c           integer*2 code
c           integer*4 buffer_address
c           integer*4 retlen_address
c        end structure
c     end structure
c     record /itmlst3_1item/ jpi_list
c     integer*4 cputime
c     jpi_list.item.code = jpi$_cputim
c     jpi_list.item.buffer_length = 4
c     jpi_list.item.buffer_address = %loc(cputime)
c     jpi_list.item.retlen_address = 0
c     call sys$getjpiw (,,,jpi_list,,,)
c     seconds = 0.01d0 * dble(cputime)
c
c     Intel PC's under the Digital Visual Fortran 5.0 compiler
c
c     real time
c     call cpu_time (time)
c     seconds = dble(time)
c
c     Intel-based PC's under the Microsoft and Watcom compilers;
c     note that this returns the wall clock time only, and will
c     fail for timing intervals spaning midnight
c
c     integer*2 hrs,mins,secs,hsecs
c     real*8 start
c     logical initial
c     save initial,start
c     data initial  / .true. /
c     if (initial) then
c        initial = .false.
c        call gettim (hrs,mins,secs,hsecs)
c        start = 3600.0d0*dble(hrs) + 60.0d0*dble(mins)
c    &              + dble(secs) + 0.01d0*dble(hsecs)
c     end if
c     call gettim (hrs,mins,secs,hsecs)
c     seconds = 3600.0d0*dble(hrs) + 60.0d0*dble(mins)
c    &             + dble(secs) + 0.01d0*dble(hsecs) - start
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cluster  --  set user-defined groups of atoms  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cluster" gets the partitioning of the system into groups
c     and stores a list of the group to which each atom belongs
c
c
      subroutine cluster
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      integer i,j,k
      integer inum,jnum
      integer last,next
      integer gnum,ga,gb
      integer nlist,list(maxatm)
      real*8 wg
      character*20 keyword
      character*80 record,string
      logical header,selected
c
c
c     set defaults for the group atom list and weights
c
      use_group = .false.
      selected = .false.
      do i = 1, n
         grplist(i) = 0
      end do
      do i = 1, maxgrp
         igrp(1,i) = 1
         igrp(2,i) = 0
      end do
      do i = 0, maxgrp
         do j = 0, maxgrp
            wgrp(j,i) = 0.0d0
         end do
      end do
c
c     get any keywords containing atom group definitions
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:6) .eq. 'GROUP ') then
            use_group = .true.
            gnum = 0
            do i = 1, n
               list(i) = 0
            end do
            call getnumb (record,gnum,next)
            if (gnum .gt. maxgrp) then
               write (iout,10)  maxgrp
   10          format (/,' CLUSTER  --  The Maximum Group Number of',
     &                    i6,' has been Exceeded')
               call fatal
            end if
            string = record(next:80)
            read (string,*,err=20,end=20)  (list(i),i=1,n)
   20       continue
            i = 1
            dowhile (list(i) .ne. 0)
               if (list(i) .gt. 0) then
                  grplist(list(i)) = gnum
                  i = i + 1
               else
                  do k = abs(list(i)), abs(list(i+1))
                     grplist(k) = gnum
                  end do
                  i = i + 2
               end if
            end do
c
c     get any keywords with weights for group interactions
c
         else if (keyword(1:13) .eq. 'SELECT-GROUP ') then
            selected = .true.
            ga = 0
            gb = 0
            wg = 0.0d0
            string = record(next:80)
            read (string,*,err=30,end=30)  ga,gb,wg
   30       continue
            if (gb .eq. 0)  gb = ga
            if (wg .eq. 0.0d0)  wg = 1.0d0
            wgrp(ga,gb) = wg
            wgrp(gb,ga) = wg
         end if
      end do
c
c     pack atoms of each group into a contiguous indexed list
c
      do i = 1, n
         list(i) = grplist(i)
      end do
      call sort3 (n,list,kgrp)
c
c     find the first and last atom in each nonempty group
c
      last = -27
      ngrp = 0
      igrp(1,ngrp) = 1
      do i = 1, n
         j = list(i)
         if (j .ne. last) then
            igrp(2,ngrp) = i - 1
            ngrp = ngrp + 1
            igrp(1,ngrp) = i
            grpnum(ngrp) = j
            last = j
         end if
      end do
      igrp(2,ngrp) = n
c
c     sort the list of atoms in each group by atom number
c
      do i = 0, maxgrp
         j = igrp(2,i) - igrp(1,i) + 1
         call sort (j,kgrp(igrp(1,i)))
      end do
c
c     renumber the group list and weights to remove empty groups
c
      do i = 1, ngrp
         do j = igrp(1,i), igrp(2,i)
            grplist(kgrp(j)) = i
         end do
      end do
      do i = 1, ngrp
         inum = grpnum(i)
         wgrp(0,i) = wgrp(0,inum)
         wgrp(i,0) = wgrp(inum,0)
         do j = 1, ngrp
            jnum = grpnum(j)
            wgrp(j,i) = wgrp(jnum,inum)
            if (jnum.ne.j .or. inum.ne.i)  wgrp(jnum,inum) = 0.0d0
         end do
      end do
c
c     use only intergroup if no group interactions were selected
c
      if (.not.selected) then
         do i = 0, ngrp
            do j = 0, ngrp
               wgrp(j,i) = 1.0d0
            end do
            wgrp(i,i) = 0.0d0
         end do
      end if
c
c     output the final list of atoms in each group
c
      if (debug .and. use_group) then
         do i = 1, ngrp
            nlist = 0
            do j = igrp(1,i), igrp(2,i)
               nlist = nlist + 1
               list(nlist) = kgrp(j)
            end do
            if (nlist .ne. 0) then
               write (iout,40)  grpnum(i)
   40          format (/,' List of Atoms Contained in Group',i3,' :',/)
               write (iout,50)  (list(j),j=1,nlist)
   50          format (3x,10i7)
            end if
         end do
      end if
c
c     output the weights for intragroup and intergroup interactions
c
      if (verbose .and. use_group) then
         header = .true.
         do i = 1, ngrp
            do j = i, ngrp
               if (wgrp(j,i) .ne. 0.0d0) then
                  if (header) then
                     header = .false.
                     write (iout,60)
   60                format (/,' Active Sets of Intra- and InterGroup',
     &                          ' Interactions :',
     &                       //,11x,'Groups',15x,'Type',14x,'Weight',/)
                  end if
                  if (i .eq. j) then
                     write (iout,70)  grpnum(i),grpnum(j),wgrp(j,i)
   70                format (5x,2i6,12x,'IntraGroup',5x,f12.4)
                  else
                     write (iout,80)  grpnum(i),grpnum(j),wgrp(j,i)
   80                format (5x,2i6,12x,'InterGroup',5x,f12.4)
                  end if
               end if
            end do
         end do
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine column  --  access Hessian elements by column  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "column" takes the off-diagonal Hessian elements stored
c     as sparse rows and sets up indices to allow column access
c
c
      subroutine column (nvar,h_init,h_stop,h_index,
     &                   c_init,c_stop,c_index,c_value,maxhess)
      implicit none
      include 'sizes.i'
      integer i,j,k,m,nvar,maxhess
      integer h_init(maxvar),h_stop(maxvar)
      integer c_init(maxvar),c_stop(maxvar)
      integer h_index(maxhess),c_index(maxhess)
      integer c_value(maxhess)
c
c
c     zero out beginning and end marker for each column
c
      do i = 1, nvar
         c_init(i) = 0
         c_stop(i) = 0
      end do
c
c     count the number of elements in each column
c
      do i = 1, nvar
         do j = h_init(i), h_stop(i)
            k = h_index(j)
            c_stop(k) = c_stop(k) + 1
         end do
      end do
c
c     set each beginning marker just beyond
c     the last element for its column
c
      c_init(1) = c_stop(1) + 1
      do i = 2, nvar
         c_init(i) = c_init(i-1) + c_stop(i)
      end do
c
c     set column index by scanning rows in reverse order
c
      do i = nvar, 1, -1
         do j = h_init(i), h_stop(i)
            k = h_index(j)
            m = c_init(k) - 1
            c_init(k) = m
            c_index(m) = i
            c_value(m) = j
         end do
      end do
c
c     convert "c_stop" from number of elements
c     in column to end marker for the column
c
      do i = 1, nvar
         c_stop(i) = c_init(i) + c_stop(i) - 1
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine command  --  get any command line arguments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "command" uses the standard Unix-like iargc/getarg routines
c     to get the number and values of arguments specified on the
c     command line at program runtime
c
c
      subroutine command
      implicit none
      include 'argue.i'
      integer i,iargc
      character*1 letter
      character*20 blank
      data blank / '                    ' /
c
c
c     initialize command line arguments as blank strings
c
      do i = 0, maxarg
         arg(i) = blank//blank//blank
      end do
c
c     get the number of arguments and store each in a string
c
      narg = iargc ()
      if (narg .gt. maxarg)  narg = maxarg
      do i = 0, narg
         call getarg (i,arg(i))
      end do
c
c     mark the command line options as unuseable for input
c
      listarg(0) = .false.
      do i = 1, narg
         listarg(i) = .true.
      end do
      do i = 1, narg-1
         letter = arg(i)(1:1)
         if (letter .eq. '-') then
            letter = arg(i)(2:2)
            call upcase (letter)
            if (letter.ge.'A' .and. letter.le.'Z') then
               listarg(i) = .false.
               listarg(i+1) = .false.
            end if
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine connct  --  attached atom list from Z-matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "connect" sets up the attached atom arrays
c     starting from a set of internal coordinates
c
c
      subroutine connct
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'zcoord.i'
      include 'zclose.i'
      integer i,j,k,id1,id2
c
c
c     zero out the number of atoms attached to each atom
c
      do i = 1, n
         n12(i) = 0
      end do
c
c     loop over the bonds in the Z-matrix, adding each bond
c     to the attach atom lists unless it is to be removed
c
      do i = 2, n
         k = iz(1,i)
         do j = 1, ndel
            id1 = idel(1,j)
            id2 = idel(2,j)
            if ((i.eq.id1 .and. k.eq.id2) .or.
     &          (i.eq.id2 .and. k.eq.id1))  goto 10
         end do
         n12(i) = n12(i) + 1
         n12(k) = n12(k) + 1
         i12(n12(i),i) = k
         i12(n12(k),k) = i
   10    continue
      end do
c
c     add any extra bonds used to make ring closures
c
      do i = 1, nadd
         do j = 1, 2
            k = iadd(j,i)
            n12(k) = n12(k) + 1
            i12(n12(k),k) = iadd(3-j,i)
         end do
      end do
c
c     sort the attached atom lists into ascending order
c
      do i = 1, n
         call sort (n12(i),i12(1,i))
      end do
      return
      end
c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 1990 by Jay William Ponder                    ##
c     ##  COPYRIGHT (C) 1985 by Scripps Clinic & Research Foundation  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine connolly  --  analytical surface area & volume  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "connolly" uses the algorithms from the AMS/VAM programs of
c     Michael Connolly to compute the analytical molecular surface
c     area and volume of a collection of spherical atoms; thus
c     it implements Fred Richards' molecular surface definition as
c     a set of analytically defined spherical and toroidal polygons
c
c     literature references:
c
c     M. L. Connolly, "Analytical Molecular Surface Calculation",
c     Journal of Applied Crystallography, 16, 548-558 (1983)
c
c     M. L. Connolly, "Computation of Molecular Volume",
c     J. Amer. Chem. Soc., 107, 1118-1124 (1985)
c
c     variables only in the Connolly routines:
c
c     na      number of atoms
c     ntt     number of temporary tori
c     nt      number of tori
c     np      number of probe positions
c     nv      number of vertices
c     nen     number of concave edges
c     nfn     number of concave faces
c     nc      number of circles
c     nep     number of convex edges
c     nfs     number of saddle faces
c     ncy     number of cycles
c     fpncy   number of cycles bounding convex face
c     nfp     number of convex faces
c     cynep   number of convex edges in cycle
c
c     a       atomic coordinates
c     ar      atomic radii
c     pr      probe radius
c     skip    if true, atom is not used
c     nosurf  if true, atom has no free surface
c     afree   atom free of neighbors
c     abur    atom buried
c
c     anbr    begin and end pointers for atoms neighbors
c     nbr     atom numbers of neighbors
c     nbrt    pointer from neighbor to torus
c
c     tta     torus atom numbers
c     ttfe    first edge of each temporary torus
c     ttle    last edge of each temporary torus
c     enext   pointer to next edge of torus
c     ttbur   torus buried
c     ttfree  torus free
c
c     t       torus center
c     tr      torus radius
c     tax     torus axis
c     ta      torus atom numbers
c     tfe     torus first edge
c     tfree   torus free of neighbors
c
c     p       probe coordinates
c     pa      probe atom numbers
c     v       vertex coordinates
c     va      vertex atom number
c     vp      vertex probe number
c     c       circle center
c     cr      circle radius
c     ca      circle atom number
c     ct      circle torus number
c
c     env     concave edge vertex numbers
c     fnen    concave face concave edge numbers
c     epc     convex edge circle number
c     epv     convex edge vertex numbers
c     afe     first convex edge of each atom
c     ale     last convex edge of each atom
c     epnext  pointer to next convex edge of atom
c     fsen    saddle face concave edge numbers
c     fsep    saddle face convex edge numbers
c     cyep    cycle convex edge numbers
c     fpa     atom number of convex face
c     fpcy    convex face cycle numbers
c
c
      subroutine connolly (volume,area,radius,probe,exclude)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'faces.i'
      include 'math.i'
      integer i
      real*8 volume,area,radius(maxatm),probe,exclude
c
c
c     set the probe radius and the number of atoms
c
      pr = probe
      na = n
c
c     set atom coordinates and radii, the excluded buffer
c     radius ("exclude") is added to atomic radii
c
      do i = 1, na
         a(1,i) = x(i)
         a(2,i) = y(i)
         a(3,i) = z(i)
         ar(i) = radius(i)
         if (ar(i) .eq. 0.0d0) then
            skip(i) = .true.
         else
            ar(i) = ar(i) + exclude
            skip(i) = .false.
         end if
      end do
c
c     find the analytical volume and surface area
c
      call neighbor
      call torus
      call place
      call compress
      call saddles
      call contact
      call vam (volume,area)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine neighbor  --  list of neighboring atom pairs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "neighbor" finds all of the neighbors of each atom
c
c     local variables :
c
c     ico      integer cube coordinates
c     icuptr   pointer to next atom in cube
c     comin    minimum atomic coordinates (cube corner)
c     icube    pointer to first atom in list for cube
c     scube    true if cube contains active atoms
c     sscube   true if cube or adjacent cubes have active atoms
c     itnl     temporary neighbor list, before sorting
c
c
      subroutine neighbor
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer maxnbra,maxcube
      parameter (maxnbra=1000)
      parameter (maxcube=40)
      integer i,j,k,m,iptr,i1,j1,k1,iatom,jatom
      integer jnbr,ici,icj,ick,jci,jcj,jck
      integer iuse,jmin,jminbr,jmold
      integer nnbr,nnbra,nbra(maxnbra),itnl(maxnbra)
      integer ico(3,maxatm),icuptr(maxatm)
      integer icube(maxcube,maxcube,maxcube)
      real*8 radmax,width,sum,sumi,dist2,d2,r2
      real*8 vect1,vect2,vect3,comin(3)
      logical scube(maxcube,maxcube,maxcube)
      logical sscube(maxcube,maxcube,maxcube)
c
c
c     ignore all atoms that are completely inside another atom;
c     may give nonsense results if this step is not taken
c
      do i = 1, na-1
         if (.not. skip(i)) then
            do j = i+1, na
               d2 = dist2(a(1,i),a(1,j))
               r2 = (ar(i) - ar(j))**2
               if (.not.skip(j) .and. d2.lt.r2) then
                  if (ar(i) .lt. ar(j)) then
                     skip(i) = .true.
                  else
                     skip(j) = .true.
                  end if
               end if
            end do
         end if
      end do
c
c     check for new coordinate minima and radii maxima
c
      radmax = 0.0d0
      do k = 1, 3
         comin(k) = a(k,1)
      end do
      do i = 1, na
         do k = 1, 3
            if (a(k,i) .lt. comin(k))  comin(k) = a(k,i)
         end do
         if (ar(i) .gt. radmax)  radmax = ar(i)
      end do
c
c     calculate width of cube from maximum
c     atom radius and probe radius
c
      width = 2.0d0 * (radmax+pr)
c
c     set up cube arrays; first the integer coordinate arrays
c
      do i = 1, na
         do k = 1, 3
            ico(k,i) = (a(k,i)-comin(k))/width + 1
            if (ico(k,i) .lt. 1) then
               call error ('Cube Coordinate Too Small')
            else if (ico(k,i) .gt. maxcube) then
               call error ('Cube Coordinate Too Large')
            end if
         end do
      end do
c
c     initialize head pointer and srn=2 arrays
c
      do i = 1, maxcube
         do j = 1, maxcube
            do k = 1, maxcube
               icube(i,j,k) = 0
               scube(i,j,k) = .false.
               sscube(i,j,k) = .false.
            end do
         end do
      end do
c
c     initialize linked list pointers
c
      do i = 1, na
         icuptr(i) = 0
      end do
c
c     set up head and later pointers for each atom
c
      do iatom = 1, na
c
c     skip atoms with surface request numbers of zero
c
         if (skip(iatom))  goto 30
         i = ico(1,iatom)
         j = ico(2,iatom)
         k = ico(3,iatom)
         if (icube(i,j,k) .le. 0) then
c
c     first atom in this cube
c
            icube(i,j,k) = iatom
         else
c
c     add to end of linked list
c
            iptr = icube(i,j,k)
   10       continue
c
c     check for duplicate atoms, turn off one of them
c
            if (dist2(a(1,iatom),a(1,iptr)) .le. 0.0d0) then
               skip(iatom) = .true.
               goto 30
            end if
c
c     move on down the list
c
            if (icuptr(iptr) .le. 0)  goto 20
            iptr = icuptr(iptr)
            goto 10
   20       continue
c
c     store atom number
c
            icuptr(iptr) = iatom
         end if
c
c     check for surfaced atom
c
         if (.not. skip(iatom))  scube(i,j,k) = .true.
   30    continue
      end do
c
c     check if this cube or any adjacent cube has active atoms
c
      do k = 1, maxcube
         do j = 1, maxcube
            do i = 1, maxcube
               if (icube(i,j,k) .ne. 0) then
                  do k1 = max(k-1,1), min(k+1,maxcube)
                     do j1 = max(j-1,1), min(j+1,maxcube)
                        do i1 = max(i-1,1), min(i+1,maxcube)
                           if (scube(i1,j1,k1)) then
                              sscube(i,j,k) = .true.
                           end if
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
      nnbr = 0
c
c     zero pointers for atom and find its cube
c
      do i = 1, na
         nnbra = 0
         nosurf(i) = skip(i)
         anbr(1,i) = 0
         anbr(2,i) = 0
         if (skip(i))  goto 70
         ici = ico(1,i)
         icj = ico(2,i)
         ick = ico(3,i)
c
c     skip iatom if its cube and adjoining
c     cubes contain only blockers
c
         if (.not. sscube(ici,icj,ick))  goto 70
         sumi = 2.0d0*pr + ar(i)
c
c     check iatom cube and adjacent cubes for neighboring atoms
c
         do jck = max(ick-1,1), min(ick+1,maxcube)
            do jcj = max(icj-1,1), min(icj+1,maxcube)
               do jci = max(ici-1,1), min(ici+1,maxcube)
                  j = icube(jci,jcj,jck)
   40             continue
c
c     check for end of linked list for this cube
c
                  if (j .le. 0)  goto 60
                  if (i .eq. j)  goto 50
                  if (skip(j))  goto 50
c
c     distance check
c
                  sum = sumi + ar(j)
                  vect1 = abs(a(1,j) - a(1,i))
                  if (vect1 .ge. sum)  goto 50
                  vect2 = abs(a(2,j) - a(2,i))
                  if (vect2 .ge. sum)  goto 50
                  vect3 = abs(a(3,j) - a(3,i))
                  if (vect3 .ge. sum)  goto 50
                  d2 = vect1**2 + vect2**2 + vect3**2
                  if (d2 .ge. sum**2)  goto 50
c
c     atoms are neighbors, save atom number in temporary array
c
                  if (.not. skip(j))  nosurf(i) = .false.
                  nnbra = nnbra + 1
                  if (nnbra .gt. maxnbra) then
                     call error ('Too many Neighbors for Atom')
                  end if
                  itnl(nnbra) = j
   50             continue
c
c     get number of next atom in cube
c
                  j = icuptr(j)
                  goto 40
   60             continue
               end do
            end do
         end do
         if (nosurf(i))  goto 70
c
c     set up neighbors arrays with jatom in increasing order
c
         jmold = 0
         do iuse = 1, nnbra
            jmin = na + 1
            do jnbr = 1, nnbra
c
c     don't use ones already sorted
c
               if (itnl(jnbr) .gt. jmold) then
                  if (itnl(jnbr) .lt. jmin) then
                     jmin = itnl(jnbr)
                     jminbr = jnbr
                  end if
               end if
            end do
            jmold = jmin
            jnbr = jminbr
            jatom = itnl(jnbr)
            nbra(iuse) = jatom
         end do
c
c     set up pointers to first and last neighbors of atom
c
         if (nnbra .gt. 0) then
            anbr(1,i) = nnbr + 1
            do m = 1, nnbra
               nnbr = nnbr + 1
               if (nnbr .gt. maxnbr) then
                  call error ('Too many Neighboring Atom Pairs')
               end if
               nbr(nnbr) = nbra(m)
            end do
            anbr(2,i) = nnbr
         end if
   70    continue
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine torus  --  position of each temporary torus  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torus" sets a list of all of the temporary torus positions
c     by testing for a torus between each atom and its neighbors
c
c
      subroutine torus
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer ia,ibeg,iend,ja,jn
      real*8 tt(3),ttr,ttax(3)
      logical ttok
c
c
c     no torus is possible if there is only one atom
c
      ntt = 0
      do ia = 1, na
         afree(ia) = .true.
      end do
      if (na .le. 1)  return
c
c     get begin and end pointers to neighbors of this atom
c
      do ia = 1, na
         if (.not. nosurf(ia)) then
            ibeg = anbr(1,ia)
            iend = anbr(2,ia)
c
c     check for no neighbors
c
            if (ibeg .gt. 0) then
               do jn = ibeg, iend
c
c     clear pointer from neighbor to torus
c
                  nbrt(jn) = 0
c
c     get atom number of neighbor
c
                  ja = nbr(jn)
c
c     don't create torus twice
c
                  if (ja .ge. ia) then
c
c     do some solid geometry
c
                     call gettor (ia,ja,ttok,tt,ttr,ttax)
                     if (ttok) then
c
c     we have a temporary torus, set up variables
c
                        ntt = ntt + 1
                        if (ntt .gt. maxtt) then
                           call error ('Too many Temporary Tori')
                        end if
c
c     mark both atoms not free
c
                        afree(ia) = .false.
                        afree(ja) = .false.
                        tta(1,ntt) = ia
                        tta(2,ntt) = ja
c
c     pointer from neighbor to torus
c
                        nbrt(jn) = ntt
c
c     initialize torus as both free and buried
c
                        ttfree(ntt) = .true.
                        ttbur(ntt) = .true.
c
c     clear pointers from torus to first and last concave edges
c
                        ttfe(ntt) = 0
                        ttle(ntt) = 0
                     end if
                  end if
               end do
            end if
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine place  --  locate positions of the probe sites  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "place" finds the probe sites by putting the probe sphere
c     tangent to each triple of neighboring atoms
c
c
      subroutine place
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer maxmut
      parameter (maxmut=500)
      integer k,ke,kv,l,l1,l2,ia,ja,ka
      integer ik,ip,jk,km,la,lm,lkf,itt
      integer nmut,iend,jend,iptr,jptr
      integer mut(maxmut),ikt(maxmut)
      integer jkt(maxmut),lknbr(maxnbr)
      real*8 dist2,d2,det,hij,hijk
      real*8 disnbr(maxmut),sumnbr(maxmut)
      real*8 uij(3),bij(3),bijk(3),tempv(3)
      real*8 uijk(3),aijk(3),pijk(3)
      logical tb,ttok,prbok
c
c
c     no possible placement if there are no temporary tori
c
      np = 0
      nfn = 0
      nen = 0
      nv = 0
      if (ntt .le. 0)  return
c
c     consider each torus in turn
c
      do itt = 1, ntt
c
c     get atom numbers
c
         ia = tta(1,itt)
         ja = tta(2,itt)
c
c     form mutual neighbor list; clear number
c     of mutual neighbors of atoms ia and ja
c
         nmut = 0
c
c     get begin and end pointers for each atom's neighbor list
c
         iptr = anbr(1,ia)
         jptr = anbr(1,ja)
         if (iptr.le.0 .or. jptr.le.0)  goto 130
         iend = anbr(2,ia)
         jend = anbr(2,ja)
c
c     collect mutual neighbors
c
   10    continue
c
c     check for end of loop
c
         if (iptr .gt. iend)  goto 40
         if (jptr .gt. jend)  goto 40
c
c     go move the lagging pointer
c
         if (nbr(iptr) .lt. nbr(jptr))  goto 20
         if (nbr(jptr) .lt. nbr(iptr))  goto 30
c
c     both point at same neighbor
c     one more mutual neighbor
c     save atom number of mutual neighbor
c
         nmut = nmut + 1
         if (nmut .gt. maxmut) then
            call error ('Too many Mutual Neighbors')
         end if
         mut(nmut) = nbr(iptr)
c
c     save pointers to second and third tori
c
         ikt(nmut) = nbrt(iptr)
         jkt(nmut) = nbrt(jptr)
   20    continue
c
c     increment pointer to ia atom neighbors
c
         iptr = iptr + 1
         goto 10
   30    continue
c
c     increment pointer to ja atom neighbors
c
         jptr = jptr + 1
         goto 10
   40    continue
c
c     we have all the mutual neighbors of ia and ja
c     if no mutual neighbors, skip to end of loop
c
         if (nmut .le. 0) then
            ttbur(itt) = .false.
            goto 130
         end if
         call gettor (ia,ja,ttok,bij,hij,uij)
         do km = 1, nmut
            ka = mut(km)
            disnbr(km) = dist2(bij,a(1,ka))
            sumnbr(km) = (pr+ar(ka))**2
c
c     initialize link to next farthest out neighbor
c
            lknbr(km) = 0
         end do
c
c     set up a linked list of neighbors in order of
c     increasing distance from ia-ja torus center
c
         lkf = 1
         if (nmut .le. 1)  goto 70
c
c     put remaining neighbors in linked list at proper position
c
         do l = 2, nmut
            l1 = 0
            l2 = lkf
   50       continue
            if (disnbr(l) .lt. disnbr(l2))  goto 60
            l1 = l2
            l2 = lknbr(l2)
            if (l2 .ne. 0)  goto 50
   60       continue
c
c     add to list
c
            if (l1 .eq. 0) then
               lkf = l
               lknbr(l) = l2
            else
               lknbr(l1) = l
               lknbr(l) = l2
            end if
         end do
   70    continue
c
c     loop thru mutual neighbors
c
         do km = 1, nmut
c
c     get atom number of neighbors
c
            ka = mut(km)
            if (skip(ia) .and. skip(ja) .and. skip(ka))  goto 120
c
c     get tori numbers for neighbor
c
            ik = ikt(km)
            jk = jkt(km)
c
c     possible new triple, do some geometry to
c     retrieve saddle center, axis and radius
c
            call getprb (ia,ja,ka,prbok,tb,bijk,hijk,uijk)
            if (tb) then
               ttbur(itt) = .true.
               ttfree(itt) = .false.
               goto 120
            end if
c
c     no duplicate triples
c
            if (ka .lt. ja)  goto 120
c
c     check whether any possible probe positions
c
            if (.not. prbok)  goto 120
c
c     altitude vector
c
            do k = 1, 3
               aijk(k) = hijk * uijk(k)
            end do
c
c     we try two probe placements
c
            do ip = 1, 2
               do k = 1, 3
                  if (ip .eq. 1) then
                     pijk(k) = bijk(k) + aijk(k)
                  else
                     pijk(k) = bijk(k) - aijk(k)
                  end if
               end do
c
c     mark three tori not free
c
               ttfree(itt) = .false.
               ttfree(ik) = .false.
               ttfree(jk) = .false.
c
c     check for collisions
c
               lm = lkf
   80          continue
               if (lm .le. 0)  goto 100
c
c     get atom number of mutual neighbor
c
               la = mut(lm)
c
c     must not equal third atom
c
               if (la .eq. ka)  goto 90
c
c     compare distance to sum of radii
c
               d2 = dist2(pijk,a(1,la))
               if (d2 .le. sumnbr(lm))  goto 110
   90          continue
               lm = lknbr(lm)
               goto 80
  100          continue
c
c     we have a new probe position
c
               np = np + 1
               if (np .gt. maxp) then
                  call error ('Too many Probe Positions')
               end if
c
c     mark three tori not buried
c
               ttbur(itt) = .false.
               ttbur(ik) = .false.
               ttbur(jk) = .false.
c
c     store probe center
c
               do k = 1, 3
                  p(k,np) = pijk(k)
               end do
c
c     calculate vectors from probe to atom centers
c
               if (nv+3 .gt. maxv)  call error ('Too many Vertices')
               do k = 1, 3
                  v(k,nv+1) = a(k,ia) - p(k,np)
                  v(k,nv+2) = a(k,ja) - p(k,np)
                  v(k,nv+3) = a(k,ka) - p(k,np)
               end do
c
c     calculate determinant of vectors defining triangle
c
               det = v(1,nv+1)*v(2,nv+2)*v(3,nv+3)
     &                  + v(1,nv+2)*v(2,nv+3)*v(3,nv+1)
     &                  + v(1,nv+3)*v(2,nv+1)*v(3,nv+2)
     &                  - v(1,nv+3)*v(2,nv+2)*v(3,nv+1)
     &                  - v(1,nv+2)*v(2,nv+1)*v(3,nv+3)
     &                  - v(1,nv+1)*v(2,nv+3)*v(3,nv+2)
c
c     now add probe coordinates to vertices
c
               do k = 1, 3
                  v(k,nv+1) = p(k,np) + v(k,nv+1)*pr/(ar(ia)+pr)
                  v(k,nv+2) = p(k,np) + v(k,nv+2)*pr/(ar(ja)+pr)
                  v(k,nv+3) = p(k,np) + v(k,nv+3)*pr/(ar(ka)+pr)
               end do
c
c     want the concave face to have counter-clockwise orientation
c
               if (det .gt. 0.0d0) then
c
c     swap second and third vertices
c
                  do k = 1, 3
                     tempv(k) = v(k,nv+2)
                     v(k,nv+2) = v(k,nv+3)
                     v(k,nv+3) = tempv(k)
                  end do
c
c     set up pointers from probe to atoms
c
                  pa(1,np) = ia
                  pa(2,np) = ka
                  pa(3,np) = ja
c
c     set up pointers from vertices to atoms
c
                  va(nv+1) = ia
                  va(nv+2) = ka
                  va(nv+3) = ja
c
c     insert concave edges into linked lists for appropriate tori
c
                  call inedge (nen+1,ik)
                  call inedge (nen+2,jk)
                  call inedge (nen+3,itt)
               else
c
c     similarly, if face already counter clockwise
c
                  pa(1,np) = ia
                  pa(2,np) = ja
                  pa(3,np) = ka
                  va(nv+1) = ia
                  va(nv+2) = ja
                  va(nv+3) = ka
                  call inedge (nen+1,itt)
                  call inedge (nen+2,jk)
                  call inedge (nen+3,ik)
               end if
c
c     set up pointers from vertices to probe
c
               do kv = 1, 3
                  vp(nv+kv) = np
               end do
c
c     set up concave edges and concave face
c
               if (nen+3 .gt. maxen) then
                  call error ('Too many Concave Edges')
               end if
c
c     edges point to vertices
c
               env(1,nen+1) = nv+1
               env(2,nen+1) = nv+2
               env(1,nen+2) = nv+2
               env(2,nen+2) = nv+3
               env(1,nen+3) = nv+3
               env(2,nen+3) = nv+1
               if (nfn+1 .gt. maxfn) then
                  call error ('Too many Concave Faces')
               end if
c
c     face points to edges
c
               do ke = 1, 3
                  fnen(ke,nfn+1) = nen + ke
               end do
c
c     increment counters for number of faces, edges and vertices
c
               nfn = nfn + 1
               nen = nen + 3
               nv = nv + 3
  110          continue
            end do
  120       continue
         end do
  130    continue
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine inedge  --  manage linked list of torus edges  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "inedge" inserts a concave edge into the
c     linked list for its temporary torus
c
c
      subroutine inedge (ien,itt)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer ien,itt,iepen
c
c
c     check for a serious error in the calling arguments
c
      if (ien .le. 0)  call error ('Bad Edge Number in INEDGE')
      if (itt .le. 0)  call error ('Bad Torus Number in INEDGE')
c
c     set beginning of list or add to end
c
      if (ttfe(itt) .eq. 0) then
         ttfe(itt) = ien
         enext(ien) = 0
         ttle(itt) = ien
      else
         iepen = ttle(itt)
         enext(iepen) = ien
         enext(ien) = 0
         ttle(itt) = ien
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine compress  --  condense temporary to final tori  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "compress" transfers only the non-buried tori from
c     the temporary tori arrays to the final tori arrays
c
c
      subroutine compress
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer itt,ia,ja,iptr,ned
      integer ip1,ip2,iv1,iv2
      logical ttok
c
c
c     initialize the number of nonburied tori
c
      nt = 0
      if (ntt .le. 0)  return
c
c     if torus is free, then it is not buried;
c     skip to end of loop if buried torus
c
      do itt = 1, ntt
         if (ttfree(itt))  ttbur(itt) = .false.
         if (.not. ttbur(itt)) then
c
c     first, transfer information
c
            nt = nt + 1
            if (nt .gt. maxt)  call error ('Too many NonBuried Tori')
            ia = tta(1,itt)
            ja = tta(2,itt)
            call gettor (ia,ja,ttok,t(1,nt),tr(nt),tax(1,nt))
            ta(1,nt) = ia
            ta(2,nt) = ja
            tfree(nt) = ttfree(itt)
            tfe(nt) = ttfe(itt)
c
c     special check for inconsistent probes
c
            iptr = tfe(nt)
            ned = 0
            dowhile (iptr .ne. 0)
               ned = ned + 1
               iptr = enext(iptr)
            end do
            if (mod(ned,2) .ne. 0) then
               iptr = tfe(nt)
               dowhile (iptr .ne. 0)
                  iv1 = env(1,iptr)
                  iv2 = env(2,iptr)
                  ip1 = vp(iv1)
                  ip2 = vp(iv2)
                  call error ('Odd Torus for Probes IP1 and IP2')
                  iptr = enext(iptr)
               end do
            end if
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine saddles  --  builds saddle pieces from tori  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "saddles" constructs circles, convex edges and saddle faces
c
c
      subroutine saddles
      implicit none
      include 'sizes.i'
      include 'faces.i'
      include 'math.i'
      integer maxent
      parameter (maxent=500)
      integer k,ia,in,ip,it,iv,itwo
      integer ien,ient,nent,l1,l2,m1,n1
      integer ten(maxent),nxtang(maxent)
      real*8 triple,factor,dtev,dt,atvect(3)
      real*8 teang(maxent),tev(3,maxent)
      logical sdstrt(maxent)
c
c
c     zero the number of circles, convex edges and saddle faces
c
      nc = 0
      nep = 0
      nfs = 0
      do ia = 1, na
         afe(ia) = 0
         ale(ia) = 0
         abur(ia) = .true.
      end do
c
c     no saddle faces if no tori
c
      if (nt .lt. 1)  return
c
c     cycle through tori
c
      do it = 1, nt
         if (skip(ta(1,it)) .and. skip(ta(2,it)))  goto 80
c
c     set up two circles
c
         do in = 1, 2
            ia = ta(in,it)
c
c     mark atom not buried
c
            abur(ia) = .false.
c
c     vector from atom to torus center
c
            do k = 1, 3
               atvect(k) = t(k,it) - a(k,ia)
            end do
            factor = ar(ia) / (ar(ia)+pr)
c
c     one more circle
c
            nc = nc + 1
            if (nc .gt. maxc)  call error ('Too many Circles')
c
c     circle center
c
            do k = 1, 3
               c(k,nc) = a(k,ia) + factor*atvect(k)
            end do
c
c     pointer from circle to atom
c
            ca(nc) = ia
c
c     pointer from circle to torus
c
            ct(nc) = it
c
c     circle radius
c
            cr(nc) = factor * tr(it)
         end do
c
c     skip to special code if free torus
c
         if (tfree(it))  goto 70
c
c     now we collect all the concave edges for this torus;
c     for each concave edge, calculate vector from torus center
c     thru probe center and the angle relative to first such vector
c
c     clear the number of concave edges for torus
c
         nent = 0
c
c     pointer to start of linked list
c
         ien = tfe(it)
   10    continue
c
c     finished if concave edge pointer is zero
c
         if (ien .le. 0)  goto 20
c
c     one more concave edge
c
         nent = nent + 1
         if (nent .gt. maxent) then
            call error ('Too many Edges for Torus')
         end if
c
c     first vertex of edge
c
         iv = env(1,ien)
c
c     probe number of vertex
c
         ip = vp(iv)
         do k = 1, 3
            tev(k,nent) = p(k,ip) - t(k,it)
         end do
         dtev = 0.0d0
         do k = 1, 3
            dtev = dtev + tev(k,nent)**2
         end do
         if (dtev .le. 0.0d0)  call error ('Probe on Torus Axis')
         dtev = sqrt(dtev)
         do k = 1, 3
            tev(k,nent) = tev(k,nent) / dtev
         end do
c
c     store concave edge number
c
         ten(nent) = ien
         if (nent .gt. 1) then
c
c     calculate angle between this vector and first vector
c
            dt = 0.0d0
            do k = 1, 3
               dt = dt + tev(k,1)*tev(k,nent)
            end do
c
c     be careful
c
            if (dt .gt. 1.0d0)  dt = 1.0d0
            if (dt .lt. -1.0d0)  dt = -1.0d0
c
c     store angle
c
            teang(nent) = acos(dt)
c
c     get the sign right
c
            if (triple(tev(1,1),tev(1,nent),tax(1,it)) .lt. 0.0d0) then
               teang(nent) = 2.0d0*pi - teang(nent)
            end if
         else
            teang(1) = 0.0d0
         end if
c
c     saddle face starts with this edge if it points parallel
c     to torus axis vector (which goes from first to second atom)
c
         sdstrt(nent) = (va(iv) .eq. ta(1,it))
c
c     next edge in list
c
         ien = enext(ien)
         goto 10
   20    continue
         if (nent .le. 0) then
            call error ('No Edges for Non-free Torus')
         end if
         itwo = 2
         if (mod(nent,itwo) .ne. 0) then
            call error ('Odd Number of Edges for Torus')
         end if
c
c     set up linked list of concave edges in order
c     of increasing angle around the torus axis;
c     clear second linked (angle-ordered) list pointers
c
         do ient = 1, nent
            nxtang(ient) = 0
         end do
         do ient = 2, nent
c
c     we have an entry to put into linked list
c     search for place to put it
c
            l1 = 0
            l2 = 1
   30       continue
            if (teang(ient) .lt. teang(l2))  goto 40
c
c     not yet, move along
c
            l1 = l2
            l2 = nxtang(l2)
            if (l2 .ne. 0)  goto 30
   40       continue
c
c     we are at end of linked list or between l1 and l2;
c     insert edge
c
            if (l1 .le. 0)  call error ('Logic Error in SADDLES')
            nxtang(l1) = ient
            nxtang(ient) = l2
         end do
c
c     collect pairs of concave edges into saddles
c     create convex edges while you're at it
c
         l1 = 1
   50    continue
         if (l1 .le. 0)  goto 60
c
c     check for start of saddle
c
         if (sdstrt(l1)) then
c
c     one more saddle face
c
            nfs = nfs + 1
            if (nfs .gt. maxfs)  call error ('Too many Saddle Faces')
c
c     get edge number
c
            ien = ten(l1)
c
c     first concave edge of saddle
c
            fsen(1,nfs) = ien
c
c     one more convex edge
c
            nep = nep + 1
            if (nep .gt. maxep)  call error ('Too many Convex Edges')
c
c     first convex edge points to second circle
c
            epc(nep) = nc
c
c     atom circle lies on
c
            ia = ca(nc)
c
c     insert convex edge into linked list for atom
c
            call ipedge (nep,ia)
c
c     first vertex of convex edge is second vertex of concave edge
c
            epv(1,nep) = env(2,ien)
c
c     first convex edge of saddle
c
            fsep(1,nfs) = nep
c
c     one more convex edge
c
            nep = nep + 1
            if (nep .gt. maxep)  call error ('Too many Convex Edges')
c
c     second convex edge points to first circle
c
            epc(nep) = nc - 1
            ia = ca(nc-1)
c
c     insert convex edge into linked list for atom
c
            call ipedge (nep,ia)
c
c     second vertex of second convex edge
c     is first vertex of first concave edge
c
            epv(2,nep) = env(1,ien)
            l1 = nxtang(l1)
c
c     wrap around
c
            if (l1 .le. 0)  l1 = 1
            if (sdstrt(l1)) then
               m1 = nxtang(l1)
               if (m1 .le. 0)  m1 = 1
               if (sdstrt(m1))  call error ('Three Starts in a Row')
               n1 = nxtang(m1)
c
c     the old switcheroo
c
               nxtang(l1) = n1
               nxtang(m1) = l1
               l1 = m1
            end if
            ien = ten(l1)
c
c     second concave edge for saddle face
c
            fsen(2,nfs) = ien
c
c     second vertex of first convex edge is
c     first vertex of second concave edge
c
            epv(2,nep-1) = env(1,ien)
c
c     first vertex of second convex edge is
c     second vertex of second concave edge
c
            epv(1,nep) = env(2,ien)
            fsep(2,nfs) = nep
c
c     quit if we have wrapped around to first edge
c
            if (l1 .eq. 1)  goto 60
         end if
c
c     next concave edge
c
         l1 = nxtang(l1)
         goto 50
   60    continue
         goto 80
c
c     free torus
c
   70    continue
c
c     set up entire circles as convex edges for new saddle surface;
c     one more saddle face
c
         nfs = nfs + 1
         if (nfs .gt. maxfs)  call error ('Too many Saddle Faces')
c
c     no concave edges for saddle
c
         fsen(1,nfs) = 0
         fsen(2,nfs) = 0
c
c     one more convex edge
c
         nep = nep + 1
         ia = ca(nc)
c
c     insert convex edge into linked list for atom
c
         call ipedge (nep,ia)
c
c     no vertices for convex edge
c
         epv(1,nep) = 0
         epv(2,nep) = 0
c
c     pointer from convex edge to second circle
c
         epc(nep) = nc
c
c     first convex edge for saddle face
c
         fsep(1,nfs) = nep
c
c     one more convex edge
c
         nep = nep + 1
         ia = ca(nc-1)
c
c     insert second convex edge into linked list
c
         call ipedge (nep,ia)
c
c     no vertices for convex edge
c
         epv(1,nep) = 0
         epv(2,nep) = 0
c
c     convex edge points to first circle
c
         epc(nep) = nc - 1
c
c     second convex edge for saddle face
c
         fsep(2,nfs) = nep
c
c     buried torus; do nothing with it
c
   80    continue
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine gettor  --  test torus site between two atoms  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "gettor" tests for a possible torus position at the interface
c     between two atoms, and finds the torus radius, center and axis
c
c
      subroutine gettor (ia,ja,ttok,torcen,torad,torax)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer k,ia,ja
      real*8 torcen(3),torad,torax(3)
      real*8 vij(3),uij(3),bij(3)
      real*8 dist2,dij,temp,temp1,temp2
      logical ttok
c
c
c     get the distance between the two atoms
c
      ttok = .false.
      dij = sqrt(dist2(a(1,ia),a(1,ja)))
c
c     find a unit vector along interatomic (torus) axis
c
      do k = 1, 3
         vij(k) = a(k,ja) - a(k,ia)
         uij(k) = vij(k) / dij
      end do
c
c     find coordinates of the center of the torus
c
      temp = 1.0d0 + ((ar(ia)+pr)**2-(ar(ja)+pr)**2)/dij**2
      do k = 1, 3
         bij(k) = a(k,ia) + 0.5d0*vij(k)*temp
      end do
c
c     skip if atoms too far apart (should not happen)
c
      temp1 = (ar(ia)+ar(ja)+2.0d0*pr)**2 - dij**2
      if (temp1 .ge. 0.0d0) then
c
c     skip if one atom is inside the other
c
         temp2 = dij**2 - (ar(ia)-ar(ja))**2
         if (temp2 .ge. 0.0d0) then
c
c     store the torus radius, center and axis
c
            ttok = .true.
            torad = sqrt(temp1*temp2) / (2.0d0*dij)
            do k = 1, 3
               torcen(k) = bij(k)
               torax(k) = uij(k)
            end do
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine getprb  --  test probe site between three atoms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "getprb" tests for a possible probe position at the interface
c     between three neighboring atoms
c
c
      subroutine getprb (ia,ja,ka,prbok,tb,bijk,hijk,uijk)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer k,ia,ja,ka
      real*8 dot,dotijk,dotut,wijk,swijk,fact
      real*8 dist2,dat2,rad,rad2,dba,rip2,hijk
      real*8 rij,rik,uij(3),uik(3),uijk(3),utb(3)
      real*8 tij(3),tik(3),bijk(3),tijik(3)
      logical prbok,tb,tok
c
c
c     initialize, then check torus over atoms "ia" and "ja"
c
      prbok = .false.
      tb = .false.
      call gettor (ia,ja,tok,tij,rij,uij)
      if (.not. tok)  return
      dat2 = dist2(a(1,ka),tij)
      rad2 = (ar(ka)+pr)**2 - rij**2
c
c     if "ka" less than "ja", then all we care about
c     is whether the torus is buried
c
      if (ka .lt. ja) then
         if (rad2 .le. 0.0d0)  return
         if (dat2 .gt. rad2)  return
      end if
      call gettor (ia,ka,tok,tik,rik,uik)
      if (.not. tok)  return
      dotijk = dot(uij,uik)
      if (dotijk .gt. 1.0d0)  dotijk = 1.0d0
      if (dotijk .lt. -1.0d0)  dotijk = -1.0d0
      wijk = acos(dotijk)
      swijk = sin(wijk)
c
c     if the three atoms are colinear, then there is no
c     probe placement; but we still care whether the torus
c     is buried by atom "k"
c
      if (swijk .eq. 0.0d0) then
         tb = (rad2.gt.0.0d0 .and. dat2.le.rad2)
         return
      end if
      call vcross (uij,uik,uijk)
      do k = 1, 3
         uijk(k) = uijk(k) / swijk
      end do
      call vcross (uijk,uij,utb)
      do k = 1, 3
         tijik(k) = tik(k) - tij(k)
      end do
      dotut = dot(uik,tijik)
      fact = dotut / swijk
      do k = 1, 3
         bijk(k) = tij(k) + utb(k)*fact
      end do
      dba = dist2(a(1,ia),bijk)
      rip2 = (ar(ia) + pr)**2
      rad = rip2 - dba
      if (rad .lt. 0.0d0) then
         tb = (rad2.gt.0.0d0 .and. dat2.le.rad2)
      else
         prbok = .true.
         hijk = sqrt(rad)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ipedge  --  manage linked list of convex edges  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ipedge" inserts convex edge into linked list for atom
c
c
      subroutine ipedge (iep,ia)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer iep,ia,iepen
c
c
c     first, check for an error condition
c
      if (iep .le. 0)  call error ('Bad Edge Number in IPEDGE')
      if (ia .le. 0)  call error ('Bad Atom Number in IPEDGE')
c
c     set beginning of list or add to end
c
      if (afe(ia) .eq. 0) then
         afe(ia) = iep
         epnext(iep) = 0
         ale(ia) = iep
      else
         iepen = ale(ia)
         epnext(iepen) = iep
         epnext(iep) = 0
         ale(ia) = iep
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine contact  --  builds exposed contact surfaces  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "contact" constructs the contact surface, cycles and convex faces
c
c
      subroutine contact
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer maxepa,maxcypa
      parameter (maxepa=300)
      parameter (maxcypa=100)
      integer i,k,ia,ia2,it,iep,ic,jc,jcy
      integer nepa,iepa,jepa,ncypa,icya,jcya,kcya
      integer ncyep,icyep,jcyep,ncyold,nused,lookv
      integer aic(maxepa),aia(maxepa),aep(maxepa),av(2,maxepa)
      integer ncyepa(maxcypa),cyepa(mxcyep,maxcypa)
      real*8 anorm,anaa,factor
      real*8 acvect(3,maxepa),aavect(3,maxepa)
      real*8 pole(3),unvect(3),acr(maxepa)
      logical ptincy,epused(maxepa),cycy(maxcypa,maxcypa)
      logical cyused(maxcypa),samef(maxcypa,maxcypa)
c
c
c     zero out the number of cycles and convex faces
c
      ncy = 0
      nfp = 0
c
c     mark all free atoms not buried
c
      do ia = 1, na
         if (afree(ia))  abur(ia) = .false.
      end do
c
c     go through all atoms
c
      do ia = 1, na
         if (skip(ia))  goto 130
c
c     skip to end of loop if buried atom
c
         if (abur(ia))  goto 130
c
c     special code for completely solvent-accessible atom
c
         if (afree(ia))  goto 120
c
c     gather convex edges for atom
c     clear number of convex edges for atom
c
         nepa = 0
c
c     pointer to first edge
c
         iep = afe(ia)
   10    continue
c
c     check whether finished gathering
c
         if (iep .le. 0)  goto 20
c
c     one more edge
c
         nepa = nepa + 1
         if (nepa .gt. maxepa) then
            call error ('Too many Convex Edges for Atom')
         end if
c
c      store vertices of edge
c
         av(1,nepa) = epv(1,iep)
         av(2,nepa) = epv(2,iep)
c
c     store convex edge number
c
         aep(nepa) = iep
         ic = epc(iep)
c
c     store circle number
c
         aic(nepa) = ic
c
c     get neighboring atom
c
         it = ct(ic)
         if (ta(1,it) .eq. ia) then
            ia2 = ta(2,it)
         else
            ia2 = ta(1,it)
         end if
c
c     store other atom number, we might need it sometime
c
         aia(nepa) = ia2
c
c     vector from atom to circle center; also
c     vector from atom to center of neighboring atom
c     sometimes we use one vector, sometimes the other
c
         do k = 1, 3
            acvect(k,nepa) = c(k,ic) - a(k,ia)
            aavect(k,nepa) = a(k,ia2) - a(k,ia)
         end do
c
c     circle radius
c
         acr(nepa) = cr(ic)
c
c     pointer to next edge
c
         iep = epnext(iep)
         goto 10
   20    continue
         if (nepa .le. 0) then
            call error ('No Edges for Non-buried, Non-free Atom')
         end if
c
c
c     form cycles; initialize all the
c     convex edges as not used in cycle
c
         do iepa = 1, nepa
            epused(iepa) = .false.
         end do
c
c     save old number of cycles
c
         ncyold = ncy
         nused = 0
         ncypa = 0
   30    continue
c
c     look for starting edge
c
         do iepa = 1, nepa
            if (.not. epused(iepa))  goto 40
         end do
c
c     cannot find starting edge, finished
c
         goto 80
   40    continue
c
c     pointer to edge
c
         iep = aep(iepa)
c
c     one edge so far for this cycle
c
         ncyep = 1
c
c     one more cycle for atom
c
         ncypa = ncypa + 1
         if (ncypa .gt. maxcypa) then
            call error ('Too many Cycles per Atom')
         end if
c
c     mark edge used in cycle
c
         epused(iepa) = .true.
         nused = nused + 1
c
c     one more cycle for molecule
c
         ncy = ncy + 1
         if (ncy .gt. maxcy)  call error ('Too many Cycles')
c
c     index of edge in atom cycle array
c
         cyepa(ncyep,ncypa) = iepa
c
c     store in molecule cycle array a pointer to edge
c
         cyep(ncyep,ncy) = iep
c
c     second vertex of this edge is the vertex to look
c     for next as the first vertex of another edge
c
         lookv = av(2,iepa)
c
c     if no vertex, this cycle is finished
c
         if (lookv .le. 0)  goto 70
   50    continue
c
c     look for next connected edge
c
         do jepa = 1, nepa
            if (epused(jepa))  goto 60
c
c     check second vertex of iepa versus first vertex of jepa
c
            if (av(1,jepa) .ne. lookv)  goto 60
c
c     edges are connected
c     pointer to edge
c
            iep = aep(jepa)
c
c     one more edge for this cycle
c
            ncyep = ncyep + 1
            if (ncyep .gt. mxcyep) then
               call error ('Too many Edges per Cycle')
            end if
            epused(jepa) = .true.
            nused = nused + 1
c
c     store index in local edge array
c
            cyepa(ncyep,ncypa) = jepa
c
c     store pointer to edge
c
            cyep(ncyep,ncy) = iep
c
c     new vertex to look for
c
            lookv = av(2,jepa)
c
c     if no vertex, this cycle is in trouble
c
            if (lookv .le. 0) then
               call error ('Pointer Error in Cycle')
            end if
            goto 50
   60       continue
         end do
c
c     it better connect to first edge of cycle
c
         if (lookv .ne. av(1,iepa)) then
            call error ('Cycle does not Close')
         end if
   70    continue
c
c     this cycle is finished
c     store number of edges in cycle
c
         ncyepa(ncypa) = ncyep
         cynep(ncy) = ncyep
         if (nused .ge. nepa)  goto 80
c
c     look for more cycles
c
         goto 30
   80    continue
c
c     compare cycles for inside/outside relation;
c     check to see if cycle i is inside cycle j
c
         do icya = 1, ncypa
            do jcya = 1, ncypa
               jcy = ncyold + jcya
c
c     initialize
c
               cycy(icya,jcya) = .true.
c
c     check for same cycle
c
               if (icya .eq. jcya)  goto 90
c
c     if cycle j has two or fewer edges, nothing can
c     lie in its exterior; i is therefore inside j
c
               if (ncyepa(jcya) .le. 2)  goto 90
c
c     if cycles i and j have a pair of edges belonging
c     to the same circle, then they are outside each other
c
               do icyep = 1, ncyepa(icya)
                  iepa = cyepa(icyep,icya)
                  ic = aic(iepa)
                  do jcyep = 1, ncyepa(jcya)
                     jepa = cyepa(jcyep,jcya)
                     jc = aic(jepa)
                     if (ic .eq. jc) then
                        cycy(icya,jcya) = .false.
                        goto 90
                     end if
                  end do
               end do
               iepa = cyepa(1,icya)
               anaa = anorm(aavect(1,iepa))
               factor = ar(ia) / anaa
c
c     north pole and unit vector pointing south
c
               do k = 1, 3
                  pole(k) = factor*aavect(k,iepa) + a(k,ia)
                  unvect(k) = -aavect(k,iepa) / anaa
               end do
               cycy(icya,jcya) = ptincy(pole,unvect,jcy)
   90          continue
            end do
         end do
c
c     group cycles into faces; direct comparison for i and j
c
         do icya = 1, ncypa
            do jcya = 1, ncypa
c
c     tentatively say that cycles i and j bound
c     the same face if they are inside each other
c
               samef(icya,jcya) = (cycy(icya,jcya) .and.
     &                               cycy(jcya,icya))
            end do
         end do
c
c     if i is in exterior of k, and k is in interior of
c     i and j, then i and j do not bound the same face
c
         do icya = 1, ncypa
            do jcya = 1, ncypa
               if (icya .ne. jcya) then
                  do kcya = 1, ncypa
                     if (kcya.ne.icya .and. kcya.ne.jcya) then
                        if (cycy(kcya,icya) .and. cycy(kcya,jcya)
     &                        .and. .not.cycy(icya,kcya)) then
                           samef(icya,jcya) = .false.
                           samef(jcya,icya) = .false.
                        end if
                     end if
                  end do
               end if
            end do
         end do
c
c     fill gaps so that "samef" falls into complete blocks
c
         do icya = 1, ncypa-2
            do jcya = icya+1, ncypa-1
               if (samef(icya,jcya)) then
                  do kcya = jcya+1, ncypa
                     if (samef(jcya,kcya)) then
                        samef(icya,kcya) = .true.
                        samef(kcya,icya) = .true.
                     end if
                  end do
               end if
            end do
         end do
c
c     group cycles belonging to the same face
c
         do icya = 1, ncypa
            cyused(icya) = .false.
         end do
c
c     clear number of cycles used in bounding faces
c
         nused = 0
         do icya = 1, ncypa
c
c     check for already used
c
            if (cyused(icya))  goto 110
c
c     one more convex face
c
            nfp = nfp + 1
            if (nfp .gt. maxfp) then
               call error ('Too many Convex Faces')
            end if
c
c     clear number of cycles for face
c
            fpncy(nfp) = 0
c
c     pointer from face to atom
c
            fpa(nfp) = ia
c
c     look for all other cycles belonging to same face
c
            do jcya = 1, ncypa
c
c     check for cycle already used in another face
c
               if (cyused(jcya))  goto 100
c
c     cycles i and j belonging to same face
c
               if (.not. samef(icya,jcya))  goto 100
c
c     mark cycle used
c
               cyused(jcya) = .true.
               nused = nused + 1
c
c     one more cycle for face
c
               fpncy(nfp) = fpncy(nfp) + 1
               if (fpncy(nfp) .gt. mxfpcy) then
                  call error ('Too many Cycles bounding Convex Face')
               end if
               i = fpncy(nfp)
c
c     store cycle number
c
               fpcy(i,nfp) = ncyold + jcya
c
c     check for finished
c
               if (nused .ge. ncypa)  goto 130
  100          continue
            end do
  110       continue
         end do
c
c     should not fall through end of do loops
c
         call error ('Not all Cycles grouped into Convex Faces')
  120    continue
c
c     one face for free atom; no cycles
c
         nfp = nfp + 1
         if (nfp .gt. maxfp) then
            call error ('Too many Convex Faces')
         end if
         fpa(nfp) = ia
         fpncy(nfp) = 0
  130    continue
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine vam  --  volumes and areas of molecules  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "vam" takes the analytical molecular surface defined
c     as a collection of spherical and toroidal polygons
c     and uses it to compute the volume and surface area
c
c
      subroutine vam (volume,area)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      integer maxdot,maxop,nscale
      parameter (maxdot=1000)
      parameter (maxop=100)
      parameter (nscale=20)
      integer k,ke,ke2,kv,ia,ic,ip,it,ien,iep
      integer ifn,ifp,ifs,iv,iv1,iv2,isc,jfn
      integer ndots,idot,nop,iop,nate,neat,neatmx
      integer ivs(3),ispind(3),ispnd2(3)
      integer nlap(maxfn),ifnop(maxop),enfs(maxen)
      integer fnt(3,maxfn),nspt(3,maxfn)
      real*8 volume,area
      real*8 atmarea(maxatm),alens,vint,vcone,vpyr,vlens
      real*8 hedron,totap,totvp,totas,totvs,totasp,totvsp
      real*8 totan,totvn,alenst,alensn,vlenst,vlensn,prism
      real*8 areap,volp,areas,vols,areasp,volsp,arean,voln
      real*8 depth,triple,dist2,areado,voldo,dot,dota
      real*8 ds2,dij2,dt,dpp,rm,rat,rsc,rho
      real*8 sumsc,sumsig,sumlam,stq,scinc,coran,corvn
      real*8 depths(maxfn),alts(3,maxfn),fncen(3,maxfn)
      real*8 cora(maxfn),corv(maxfn),cenop(3,maxop)
      real*8 sdot(3),dotv(nscale),fnvect(3,3,maxfn)
      real*8 tau(3),ppm(3),xpnt1(3),xpnt2(3)
      real*8 qij(3),qji(3),vects(3,3)
      real*8 vect1(3),vect2(3),vect3(3),vect4(3)
      real*8 vect5(3),vect6(3),vect7(3),vect8(3)
      real*8 upp(3),thetaq(3),sigmaq(3)
      real*8 umq(3),upq(3),uc(3),uq(3),uij(3)
      real*8 dots(3,maxdot),tdots(3,maxdot)
      logical fcins(3,maxfn),fcint(3,maxfn)
      logical cinsp,cintp,usenum,vip(3),ate(maxop)
      logical spindl,alli,allj,anyi,anyj,case1,case2
      logical fntrev(3,maxfn),badav(maxfn),badt(maxfn)
c
c
c     compute the volume of the interior polyhedron
c
      hedron = 0.0d0
      do ifn = 1, nfn
         call measpm (ifn,prism)
         hedron = hedron + prism
      end do
c
c     compute the area and volume due to convex faces
c     as well as the area partitioned among the atoms
c
      totap = 0.0d0
      totvp = 0.0d0
      do ia = 1, na
         atmarea(ia) = 0.0d0
      end do
      do ifp = 1, nfp
         call measfp (ifp,areap,volp)
         ia = fpa(ifp)
         atmarea(ia) = atmarea(ia) + areap
         totap = totap + areap
         totvp = totvp + volp
      end do
c
c     compute the area and volume due to saddle faces
c     as well as the spindle correction value
c
      totas = 0.0d0
      totvs = 0.0d0
      totasp = 0.0d0
      totvsp = 0.0d0
      do ifs = 1, nfs
         do k = 1, 2
            ien = fsen(k,ifs)
            if (ien .gt. 0)  enfs(ien) = ifs
         end do
         call measfs (ifs,areas,vols,areasp,volsp)
         totas = totas + areas
         totvs = totvs + vols
         totasp = totasp + areasp
         totvsp = totvsp + volsp
         if (areas-areasp .lt. 0.0d0) then
            call error ('Negative Area for Saddle Face')
         end if
      end do
c
c     compute the area and volume due to concave faces
c
      totan = 0.0d0
      totvn = 0.0d0
      do ifn = 1, nfn
         call measfn (ifn,arean,voln)
         totan = totan + arean
         totvn = totvn + voln
      end do
c
c     compute the area and volume lens correction values
c
      alenst = 0.0d0
      alensn = 0.0d0
      vlenst = 0.0d0
      vlensn = 0.0d0
      if (pr .le. 0.0d0)  goto 140
      ndots = maxdot
      call gendot (ndots,dots,pr,0.0d0,0.0d0,0.0d0)
      dota = (4.0d0 * pi * pr**2) / ndots
      do ifn = 1, nfn
         nlap(ifn) = 0
         cora(ifn) = 0.0d0
         corv(ifn) = 0.0d0
         badav(ifn) = .false.
         badt(ifn) = .false.
         do k = 1, 3
            nspt(k,ifn) = 0
         end do
         ien = fnen(1,ifn)
         iv = env(1,ien)
         ip = vp(iv)
         depths(ifn) = depth(ip,alts(1,ifn))
         do k = 1, 3
            fncen(k,ifn) = p(k,ip)
         end do
         ia = va(iv)
c
c     get vertices and vectors
c
         do ke = 1, 3
            ien = fnen(ke,ifn)
            ivs(ke) = env(1,ien)
            ia = va(ivs(ke))
            ifs = enfs(ien)
            iep = fsep(1,ifs)
            ic = epc(iep)
            it = ct(ic)
            fnt(ke,ifn) = it
            fntrev(ke,ifn) = (ta(1,it) .ne. ia)
         end do
         do ke = 1, 3
            do k = 1, 3
               vects(k,ke) = v(k,ivs(ke)) - p(k,ip)
            end do
         end do
c
c     calculate normal vectors for the three planes
c     that cut out the geodesic triangle
c
         call vcross (vects(1,1),vects(1,2),fnvect(1,1,ifn))
         call vnorm (fnvect(1,1,ifn),fnvect(1,1,ifn))
         call vcross (vects(1,2),vects(1,3),fnvect(1,2,ifn))
         call vnorm (fnvect(1,2,ifn),fnvect(1,2,ifn))
         call vcross (vects(1,3),vects(1,1),fnvect(1,3,ifn))
         call vnorm (fnvect(1,3,ifn),fnvect(1,3,ifn))
      end do
      do ifn = 1, nfn-1
         do jfn = ifn+1, nfn
            dij2 = dist2(fncen(1,ifn),fncen(1,jfn))
            if (dij2 .gt. 4.0d0*pr**2)  goto 90
            if (depths(ifn).gt.pr .and. depths(jfn).gt.pr)  goto 90
c
c     these two probes may have intersecting surfaces
c
            dpp = sqrt(dist2(fncen(1,ifn),fncen(1,jfn)))
c
c     compute the midpoint
c
            do k = 1, 3
               ppm(k) = (fncen(k,ifn) + fncen(k,jfn)) / 2.0d0
               upp(k) = (fncen(k,jfn) - fncen(k,ifn)) / dpp
            end do
            rm = pr**2 - (dpp/2.0d0)**2
            if (rm .lt. 0.0d0)  rm = 0.0d0
            rm = sqrt(rm)
            rat = dpp / (2.0d0*pr)
            if (rat .gt. 1.0d0)  rat = 1.0d0
            if (rat .lt. -1.0d0)  rat = -1.0d0
            rho = asin(rat)
c
c     use circle-plane intersection routine
c
            alli = .true.
            anyi = .false.
            spindl = .false.
            do k = 1, 3
               ispind(k) = 0
               ispnd2(k) = 0
            end do
            do ke = 1, 3
               thetaq(ke) = 0.0d0
               sigmaq(ke) = 0.0d0
               tau(ke) = 0.0d0
               call cirpln (ppm,rm,upp,fncen(1,ifn),fnvect(1,ke,ifn),
     &                              cinsp,cintp,xpnt1,xpnt2)
               fcins(ke,ifn) = cinsp
               fcint(ke,ifn) = cintp
               if (.not. cinsp)  alli = .false.
               if (cintp)  anyi = .true.
               if (.not. cintp)  goto 10
               it = fnt(ke,ifn)
               if (tr(it) .gt. pr)  goto 10
               do ke2 = 1, 3
                  if (it .eq. fnt(ke2,jfn)) then
                     ispind(ke) = it
                     nspt(ke,ifn) = nspt(ke,ifn) + 1
                     ispnd2(ke2) = it
                     nspt(ke2,jfn) = nspt(ke2,jfn) + 1
                     spindl = .true.
                  end if
               end do
               if (ispind(ke) .eq. 0)  goto 10
c
c     check that the two ways of calculating
c     intersection points match
c
               rat = tr(it) / pr
               if (rat .gt. 1.0d0)  rat = 1.0d0
               if (rat .lt. -1.0d0)  rat = -1.0d0
               thetaq(ke) = acos(rat)
               stq = sin(thetaq(ke))
               if (fntrev(ke,ifn)) then
                  do k = 1, 3
                     uij(k) = -tax(k,it)
                  end do
               else
                  do k = 1, 3
                     uij(k) = tax(k,it)
                  end do
               end if
               do k = 1, 3
                  qij(k) = t(k,it) - stq * pr * uij(k)
                  qji(k) = t(k,it) + stq * pr * uij(k)
               end do
               do k = 1, 3
                  umq(k) = (qij(k) - ppm(k)) / rm
                  upq(k) = (qij(k) - fncen(k,ifn)) / pr
               end do
               call vcross (uij,upp,vect1)
               dt = dot(umq,vect1)
               if (dt .gt. 1.0d0)  dt = 1.0d0
               if (dt .lt. -1.0d0)  dt = -1.0d0
               sigmaq(ke) = acos(dt)
               call vcross (upq,fnvect(1,ke,ifn),vect1)
               call vnorm (vect1,uc)
               call vcross (upp,upq,vect1)
               call vnorm (vect1,uq)
               dt = dot(uc,uq)
               if (dt .gt. 1.0d0)  dt = 1.0d0
               if (dt .lt. -1.0d0)  dt = -1.0d0
               tau(ke) = pi - acos(dt)
   10          continue
            end do
            allj = .true.
            anyj = .false.
            do ke = 1, 3
               call cirpln (ppm,rm,upp,fncen(1,jfn),fnvect(1,ke,jfn),
     &                              cinsp,cintp,xpnt1,xpnt2)
               fcins(ke,jfn) = cinsp
               fcint(ke,jfn) = cintp
               if (.not. cinsp)  allj = .false.
               if (cintp)  anyj = .true.
            end do
            case1 = (alli .and. allj .and. .not.anyi .and. .not.anyj)
            case2 = (anyi .and. anyj .and. spindl)
            if (.not.case1 .and. .not.case2)  goto 90
c
c     this kind of overlap can be handled
c
            nlap(ifn) = nlap(ifn) + 1
            nlap(jfn) = nlap(jfn) + 1
            do ke = 1, 3
               ien = fnen(ke,ifn)
               iv1 = env(1,ien)
               iv2 = env(2,ien)
               do k = 1, 3
                  vect3(k) = v(k,iv1) - fncen(k,ifn)
                  vect4(k) = v(k,iv2) - fncen(k,ifn)
               end do
               do ke2 = 1, 3
                  if (ispind(ke) .eq. ispnd2(ke2))  goto 40
                  if (ispind(ke) .eq. 0)  goto 40
                  call cirpln (fncen(1,ifn),pr,fnvect(1,ke,ifn),
     &                           fncen(1,jfn),fnvect(1,ke2,jfn),
     &                           cinsp,cintp,xpnt1,xpnt2)
                  if (.not. cintp)  goto 40
                  ien = fnen(ke2,jfn)
                  iv1 = env(1,ien)
                  iv2 = env(2,ien)
                  do k = 1, 3
                     vect7(k) = v(k,iv1) - fncen(k,jfn)
                     vect8(k) = v(k,iv2) - fncen(k,jfn)
                  end do
c
c     check whether point lies on spindle arc
c
                  do k = 1, 3
                     vect1(k) = xpnt1(k) - fncen(k,ifn)
                     vect2(k) = xpnt2(k) - fncen(k,ifn)
                     vect5(k) = xpnt1(k) - fncen(k,jfn)
                     vect6(k) = xpnt2(k) - fncen(k,jfn)
                  end do
                  if (triple(vect3,vect1,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 20
                  if (triple(vect1,vect4,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 20
                  if (triple(vect7,vect5,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 20
                  if (triple(vect5,vect8,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 20
                  goto 30
   20             continue
                  if (triple(vect3,vect2,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 40
                  if (triple(vect2,vect4,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 40
                  if (triple(vect7,vect6,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 40
                  if (triple(vect6,vect8,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 40
   30             continue
                  badav(ifn) = .true.
   40             continue
               end do
            end do
            do ke = 1, 3
               ien = fnen(ke,ifn)
               iv1 = env(1,ien)
               iv2 = env(2,ien)
               do k = 1, 3
                  vect3(k) = v(k,iv1) - fncen(k,ifn)
                  vect4(k) = v(k,iv2) - fncen(k,ifn)
               end do
               do ke2 = 1, 3
                  if (ispind(ke) .eq. ispnd2(ke2))  goto 70
                  if (ispnd2(ke2) .eq. 0)  goto 70
                  call cirpln (fncen(1,jfn),pr,fnvect(1,ke2,jfn),
     &                           fncen(1,ifn),fnvect(1,ke,ifn),
     &                           cinsp,cintp,xpnt1,xpnt2)
                  if (.not. cintp)  goto 70
                  ien = fnen(ke2,jfn)
                  iv1 = env(1,ien)
                  iv2 = env(2,ien)
                  do k = 1, 3
                     vect7(k) = v(k,iv1) - fncen(k,jfn)
                     vect8(k) = v(k,iv2) - fncen(k,jfn)
                  end do
c
c     check whether point lies on spindle arc
c
                  do k = 1, 3
                     vect1(k) = xpnt1(k) - fncen(k,ifn)
                     vect2(k) = xpnt2(k) - fncen(k,ifn)
                     vect5(k) = xpnt1(k) - fncen(k,jfn)
                     vect6(k) = xpnt2(k) - fncen(k,jfn)
                  end do
                  if (triple(vect3,vect1,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 50
                  if (triple(vect1,vect4,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 50
                  if (triple(vect7,vect5,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 50
                  if (triple(vect5,vect8,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 50
                  goto 60
   50             continue
                  if (triple(vect3,vect2,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 70
                  if (triple(vect2,vect4,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 70
                  if (triple(vect7,vect6,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 70
                  if (triple(vect6,vect8,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 70
   60             continue
                  badav(jfn) = .true.
   70             continue
               end do
            end do
            sumlam = 0.0d0
            sumsig = 0.0d0
            sumsc = 0.0d0
            do ke = 1, 3
               if (ispind(ke) .ne. 0) then
                  sumlam = sumlam + pi - tau(ke)
                  sumsig = sumsig + sigmaq(ke) - pi
                  sumsc = sumsc + sin(sigmaq(ke))*cos(sigmaq(ke))
               end if
            end do
            alens = 2.0d0 * pr**2 * (pi - sumlam - sin(rho)*(pi+sumsig))
            vint = alens * pr / 3.0d0
            vcone = pr * rm**2 * sin(rho) * (pi+sumsig) / 3.0d0
            vpyr =  pr * rm**2 * sin(rho) * sumsc / 3.0d0
            vlens = vint - vcone + vpyr
            cora(ifn) = cora(ifn) + alens
            cora(jfn) = cora(jfn) + alens
            corv(ifn) = corv(ifn) + vlens
            corv(jfn) = corv(jfn) + vlens
c
c     check for vertex on opposing probe in face
c
            do kv = 1, 3
               vip(kv) = .false.
               ien = fnen(kv,jfn)
               iv = env(1,ien)
               do k = 1, 3
                  vect1(k) = v(k,iv) - fncen(k,ifn)
               end do
               call vnorm (vect1,vect1)
               do ke = 1, 3
                  dt = dot(fnvect(1,ke,ifn),v(1,iv))
                  if (dt .gt. 0.0d0)  goto 80
               end do
               vip(kv) = .true.
   80          continue
            end do
   90       continue
         end do
      end do
      do ifn = 1, nfn
         do ke = 1, 3
            if (nspt(ke,ifn) .gt. 1)  badt(ifn) = .true.
         end do
      end do
      do ifn = 1, nfn
         if (nlap(ifn) .le. 0)  goto 130
c
c     gather all overlapping probes
c
         nop = 0
         do jfn = 1, nfn
            if (ifn .ne. jfn) then
               dij2 = dist2(fncen(1,ifn),fncen(1,jfn))
               if (dij2 .le. 4.0d0*pr**2) then
                  if (depths(jfn) .le. pr) then
                     nop = nop + 1
                     if (nop .gt. maxop) then
                        call error ('NOP Overflow in VAM')
                     end if
                     ifnop(nop) = jfn
                     do k = 1, 3
                        cenop(k,nop) = fncen(k,jfn)
                     end do
                  end if
               end if
            end if
         end do
c
c     numerical calculation of the correction
c
         areado = 0.0d0
         voldo = 0.0d0
         scinc = 1.0d0 / nscale
         do isc = 1, nscale
            rsc = isc - 0.5d0
            dotv(isc) = pr * dota * rsc**2 * scinc**3
         end do
         do iop = 1, nop
            ate(iop) = .false.
         end do
         neatmx = 0
         do idot = 1, ndots
            do ke = 1, 3
               dt = dot(fnvect(1,ke,ifn),dots(1,idot))
               if (dt .gt. 0.0d0)  goto 120
            end do
            do k = 1, 3
               tdots(k,idot) = fncen(k,ifn) + dots(k,idot)
            end do
            do iop = 1, nop
               jfn = ifnop(iop)
               ds2 = dist2(tdots(1,idot),fncen(1,jfn))
               if (ds2 .lt. pr**2) then
                  areado = areado + dota
                  goto 100
               end if
            end do
  100       continue
            do isc = 1, nscale
               rsc = isc - 0.5d0
               do k = 1, 3
                  sdot(k) = fncen(k,ifn) + rsc*scinc*dots(k,idot)
               end do
               neat = 0
               do iop = 1, nop
                  jfn = ifnop(iop)
                  ds2 = dist2(sdot,fncen(1,jfn))
                  if (ds2 .lt. pr**2) then
                     do k = 1, 3
                        vect1(k) = sdot(k) - fncen(k,jfn)
                     end do
                     do ke = 1, 3
                        dt = dot(fnvect(1,ke,jfn),vect1)
                        if (dt .gt. 0.0d0)  goto 110
                     end do
                     neat = neat + 1
                     ate(iop) = .true.
  110                continue
                  end if
               end do
               if (neat .gt. neatmx)  neatmx = neat
               if (neat .gt. 0) then
                  voldo = voldo + dotv(isc) * (neat/(1.0d0+neat))
               end if
            end do
  120       continue
         end do
         coran = areado
         corvn = voldo
         nate = 0
         do iop = 1, nop
            if (ate(iop))  nate = nate + 1
         end do
c
c     use either the analytical or numerical correction
c
         usenum = (nate.gt.nlap(ifn) .or. neatmx.gt.1 .or. badt(ifn))
         if (usenum) then
            cora(ifn) = coran
            corv(ifn) = corvn
            alensn = alensn + cora(ifn)
            vlensn = vlensn + corv(ifn)
         else if (badav(ifn)) then
            corv(ifn) = corvn
            vlensn = vlensn + corv(ifn)
         end if
         alenst = alenst + cora(ifn)
         vlenst = vlenst + corv(ifn)
  130    continue
      end do
  140 continue
c
c     print out the decomposition of the area and volume
c
      if (debug) then
         write (iout,150)
  150    format (/,' Convex Surface Area for each Atom :',/)
         k = 1
         dowhile (k .le. na)
            write (iout,160)  (ia,atmarea(ia),ia=k,min(k+4,na))
  160       format (1x,5(i7,f8.3))
            k = k + 5
         end do
         write (iout,170)  nfp,totap,totvp
  170    format (/,' Convex Faces :',i12,5x,'Area :',f13.3,
     &              4x,'Volume :',f13.3)
         write (iout,180)  nfs,totas,totvs
  180    format (' Saddle Faces :',i12,5x,'Area :',f13.3,
     &              4x,'Volume :',f13.3)
         write (iout,190)  nfn,totan,totvn
  190    format (' Concave Faces :',i11,5x,'Area :',f13.3,
     &              4x,'Volume :',f13.3)
         write (iout,200)  hedron
  200    format (' Buried Polyhedra :',36x,'Volume :',f13.3)
         if (totasp.ne.0.0d0 .or. totvsp.ne.0.0d0 .or.
     &       alenst.ne.0.0d0 .or. vlenst.ne.0.0d0) then
            write (iout,210)  -totasp,-totvsp
  210       format (/,' Spindle Correction :',11x,'Area :',f13.3,
     &                 4x,'Volume :',f13.3)
            write (iout,220)  -alenst-alensn,vlenst-vlensn
  220       format (' Lens Analytical Correction :',3x,'Area :',f13.3,
     &                 4x,'Volume :',f13.3)
         end if
         if (alensn.ne.0.0d0 .or. vlensn.ne.0.0d0) then
            write (iout,230)  alensn,vlensn
  230       format (' Lens Numerical Correction :',4x,'Area :',f13.3,
     &                 4x,'Volume :',f13.3)
         end if
      end if
c
c     finally, compute the total area and total volume
c
      area = totap + totas + totan - totasp - alenst
      volume = totvp + totvs + totvn + hedron - totvsp + vlenst
      return
      end
c
c
c     ######################
c     ##                  ##
c     ##  function depth  ##
c     ##                  ##
c     ######################
c
c
      function depth (ip,alt)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer k,ip,ia1,ia2,ia3
      real*8 depth,dot,alt(3)
      real*8 vect1(3),vect2(3),vect3(3),vect4(3)
c
c
      ia1 = pa(1,ip)
      ia2 = pa(2,ip)
      ia3 = pa(3,ip)
      do k = 1, 3
         vect1(k) = a(k,ia1) - a(k,ia3)
         vect2(k) = a(k,ia2) - a(k,ia3)
         vect3(k) = p(k,ip) - a(k,ia3)
      end do
      call vcross (vect1,vect2,vect4)
      call vnorm (vect4,vect4)
      depth = dot(vect4,vect3)
      do k = 1, 3
         alt(k) = vect4(k)
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine measpm  --  volume of interior polyhedron  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "measpm" computes the volume of a single prism section of
c     the full interior polyhedron
c
c
      subroutine measpm (ifn,prism)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer k,ke,ien,iv,ia,ip,ifn
      real*8 prism,height,pav(3,3)
      real*8 vect1(3),vect2(3),vect3(3)
c
c
      height = 0.0d0
      do ke = 1, 3
         ien = fnen(ke,ifn)
         iv = env(1,ien)
         ia = va(iv)
         height = height + a(3,ia)
         ip = vp(iv)
         do k = 1, 3
            pav(k,ke) = a(k,ia) - p(k,ip)
         end do
      end do
      height = height / 3.0d0
      do k = 1, 3
         vect1(k) = pav(k,2) - pav(k,1)
         vect2(k) = pav(k,3) - pav(k,1)
      end do
      call vcross (vect1,vect2,vect3)
      prism = height * vect3(3) / 2.0d0
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine measfp  ##
c     ##                     ##
c     #########################
c
c
      subroutine measfp (ifp,areap,volp)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      include 'math.i'
      integer k,ke,ifp,iep,ia,ia2,ic,it,iv1,iv2
      integer ncycle,ieuler,icyptr,icy,nedge
      real*8 areap,volp,dot,dt,gauss
      real*8 vecang,angle,geo,pcurve,gcurve
      real*8 vect1(3),vect2(3),acvect(3),aavect(3)
      real*8 tanv(3,2,mxcyep),radial(3,mxcyep)
c
c
      ia = fpa(ifp)
      pcurve = 0.0d0
      gcurve = 0.0d0
      ncycle = fpncy(ifp)
      if (ncycle .gt. 0) then
         ieuler = 2 - ncycle
      else
         ieuler = 2
      end if
      do icyptr = 1, ncycle
         icy = fpcy(icyptr,ifp)
         nedge = cynep(icy)
         do ke = 1, nedge
            iep = cyep(ke,icy)
            ic = epc(iep)
            it = ct(ic)
            if (ia .eq. ta(1,it)) then
               ia2 = ta(2,it)
            else
               ia2 = ta(1,it)
            end if
            do k = 1, 3
               acvect(k) = c(k,ic) - a(k,ia)
               aavect(k) = a(k,ia2) - a(k,ia)
            end do
            call vnorm (aavect,aavect)
            dt = dot(acvect,aavect)
            geo = -dt / (ar(ia)*cr(ic))
            iv1 = epv(1,iep)
            iv2 = epv(2,iep)
            if (iv1.eq.0 .or. iv2.eq.0) then
               angle = 2.0d0 * pi
            else
               do k = 1, 3
                  vect1(k) = v(k,iv1) - c(k,ic)
                  vect2(k) = v(k,iv2) - c(k,ic)
                  radial(k,ke) = v(k,iv1) - a(k,ia)
               end do
               call vnorm (radial(1,ke),radial(1,ke))
               call vcross (vect1,aavect,tanv(1,1,ke))
               call vnorm (tanv(1,1,ke),tanv(1,1,ke))
               call vcross (vect2,aavect,tanv(1,2,ke))
               call vnorm (tanv(1,2,ke),tanv(1,2,ke))
               angle = vecang(vect1,vect2,aavect,-1.0d0)
            end if
            gcurve = gcurve + cr(ic)*angle*geo
            if (nedge .ne. 1) then
               if (ke .gt. 1) then
                  angle = vecang(tanv(1,2,ke-1),tanv(1,1,ke),
     &                               radial(1,ke),1.0d0)
                  if (angle .lt. 0.0d0) then
                     call error ('Negative Angle in MEASFP')
                  end if
                  pcurve = pcurve + angle
               end if
            end if
         end do
         if (nedge .gt. 1) then
            angle = vecang(tanv(1,2,nedge),tanv(1,1,1),
     &                         radial(1,1),1.0d0)
            if (angle .lt. 0.0d0) then
               call error ('Negative Angle in MEASFP')
            end if
            pcurve = pcurve + angle
         end if
      end do
      gauss = 2.0d0*pi*ieuler - pcurve - gcurve
      areap = gauss * (ar(ia)**2)
      volp = areap * ar(ia) / 3.0d0
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine measfs  ##
c     ##                     ##
c     #########################
c
c
      subroutine measfs (ifs,areas,vols,areasp,volsp)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      include 'math.i'
      integer k,ifs,iep,ic,ic1,ic2,it,ia1,ia2,iv1,iv2
      real*8 areas,vols,areasp,volsp,vecang,phi
      real*8 dot,d1,d2,w1,w2,theta1,theta2,rat,thetaq
      real*8 cone1,cone2,term1,term2,term3,spin,volt
      real*8 vect1(3),vect2(3),aavect(3)
      logical cusp
c
c
      iep = fsep(1,ifs)
      ic = epc(iep)
      it = ct(ic)
      ia1 = ta(1,it)
      ia2 = ta(2,it)
      do k = 1, 3
         aavect(k) = a(k,ia2) - a(k,ia1)
      end do
      call vnorm (aavect,aavect)
      iv1 = epv(1,iep)
      iv2 = epv(2,iep)
      if (iv1.eq.0 .or. iv2.eq.0) then
         phi = 2.0d0 * pi
      else
         do k = 1, 3
            vect1(k) = v(k,iv1) - c(k,ic)
            vect2(k) = v(k,iv2) - c(k,ic)
         end do
         phi = vecang(vect1,vect2,aavect,1.0d0)
      end if
      do k = 1, 3
         vect1(k) = a(k,ia1) - t(k,it)
         vect2(k) = a(k,ia2) - t(k,it)
      end do
      d1 = -dot(vect1,aavect)
      d2 = dot(vect2,aavect)
      theta1 = atan2(d1,tr(it))
      theta2 = atan2(d2,tr(it))
c
c     check for cusps
c
      if (tr(it).lt.pr .and. theta1.gt.0.0d0
     &                 .and. theta2.gt.0.0d0) then
         cusp = .true.
         rat = tr(it) / pr
         if (rat .gt. 1.0d0)  rat = 1.0d0
         if (rat .lt. -1.0d0)  rat = -1.0d0
         thetaq = acos(rat)
      else
         cusp = .false.
         thetaq = 0.0d0
         areasp = 0.0d0
         volsp = 0.0d0
      end if
      term1 = tr(it) * pr * (theta1+theta2)
      term2 = (pr**2) * (sin(theta1) + sin(theta2))
      areas = phi * (term1-term2)
      if (cusp) then
         spin = tr(it)*pr*thetaq - pr**2 * sin(thetaq)
         areasp = 2.0d0 * phi * spin
      end if
c
      iep = fsep(1,ifs)
      ic2 = epc(iep)
      iep = fsep(2,ifs)
      ic1 = epc(iep)
      if (ca(ic1) .ne. ia1) then
         call error ('IA1 Inconsistency in MEASFS')
      end if
      do k = 1, 3
         vect1(k) = c(k,ic1) - a(k,ia1)
         vect2(k) = c(k,ic2) - a(k,ia2)
      end do
      w1 = dot(vect1,aavect)
      w2 = -dot(vect2,aavect)
      cone1 = phi * (w1*cr(ic1)**2)/6.0d0
      cone2 = phi * (w2*cr(ic2)**2)/6.0d0
      term1 = (tr(it)**2) * pr * (sin(theta1)+sin(theta2))
      term2 = sin(theta1)*cos(theta1) + theta1
     &          + sin(theta2)*cos(theta2) + theta2
      term2 = tr(it) * (pr**2) * term2
      term3 = sin(theta1)*cos(theta1)**2 + 2.0d0*sin(theta1)
     &          + sin(theta2)*cos(theta2)**2 + 2.0d0*sin(theta2)
      term3 = (pr**3 / 3.0d0) * term3
      volt = (phi/2.0d0) * (term1-term2+term3)
      vols = volt + cone1 + cone2
      if (cusp) then
         term1 = (tr(it)**2) * pr * sin(thetaq)
         term2 = sin(thetaq)*cos(thetaq) + thetaq
         term2 = tr(it) * (pr**2) * term2
         term3 = sin(thetaq)*cos(thetaq)**2 + 2.0d0*sin(thetaq)
         term3 = (pr**3 / 3.0d0) * term3
         volsp = phi * (term1-term2+term3)
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine measfn  ##
c     ##                     ##
c     #########################
c
c
      subroutine measfn (ifn,arean,voln)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      include 'math.i'
      integer k,ke,je,ifn,ien,iv,ia,ip
      real*8 arean,voln,vecang,triple,defect,simplx
      real*8 pvv(3,3),pav(3,3),planev(3,3),angle(3)
c
c
      do ke = 1, 3
         ien = fnen(ke,ifn)
         iv = env(1,ien)
         ia = va(iv)
         ip = vp(iv)
         do k = 1, 3
            pvv(k,ke) = v(k,iv) - p(k,ip)
            pav(k,ke) = a(k,ia) - p(k,ip)
         end do
         if (pr .gt. 0.0d0)  call vnorm (pvv(1,ke),pvv(1,ke))
      end do
      if (pr .le. 0.0d0) then
         arean = 0.0d0
      else
         do ke = 1, 3
            je = ke + 1
            if (je .gt. 3)  je = 1
            call vcross (pvv(1,ke),pvv(1,je),planev(1,ke))
            call vnorm (planev(1,ke),planev(1,ke))
         end do
         do ke = 1, 3
            je = ke - 1
            if (je .lt. 1)  je = 3
            angle(ke) = vecang(planev(1,je),planev(1,ke),
     &                             pvv(1,ke),-1.0d0)
            if (angle(ke) .lt. 0.0d0) then
               call error ('Negative Angle in MEASFN')
            end if
         end do
         defect = 2.0d0*pi - (angle(1)+angle(2)+angle(3))
         arean = (pr**2) * defect
      end if
      simplx = -triple(pav(1,1),pav(1,2),pav(1,3)) / 6.0d0
      voln = simplx - arean*pr/3.0d0
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine projct  ##
c     ##                     ##
c     #########################
c
c
      subroutine projct (pnt,unvect,icy,ia,spv,nedge,fail)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer k,ke,icy,ia,nedge,iep,iv
      real*8 dot,dt,f,polev(3),pnt(3)
      real*8 unvect(3),spv(3,mxcyep)
      logical fail
c
c
      fail = .false.
      nedge = cynep(icy)
      do ke = 1, cynep(icy)
c
c     vertex number (use first vertex of edge)
c
         iep = cyep(ke,icy)
         iv = epv(1,iep)
         if (iv .ne. 0) then
c
c     vector from north pole to vertex
c
            do k = 1, 3
               polev(k) = v(k,iv) - pnt(k)
            end do
c
c     calculate multiplication factor
c
            dt = dot(polev,unvect)
            if (dt .eq. 0.0d0) then
               fail = .true.
               return
            end if
            f = (ar(ia)*2) / dt
            if (f .lt. 1.0d0) then
               fail = .true.
               return
            end if
c
c     projected vertex for this convex edge
c
            do k = 1, 3
               spv(k,ke) = pnt(k) + f*polev(k)
               continue
            end do
         end if
      end do
      return
      end
c
c
c     #######################
c     ##                   ##
c     ##  function ptincy  ##
c     ##                   ##
c     #######################
c
c
      function ptincy (pnt,unvect,icy)
      implicit none
      include 'sizes.i'
      include 'faces.i'
      integer k,ke,icy,iep,ic,it,iatom,iaoth,nedge
      real*8 unvect(3),dot,rotang,totang
      real*8 spv(3,mxcyep),epu(3,mxcyep)
      real*8 pnt(3),acvect(3),cpvect(3)
      logical ptincy,fail
c
c
c     check for eaten by neighbor
c
      do ke = 1, cynep(icy)
         iep = cyep(ke,icy)
         ic = epc(iep)
         it = ct(ic)
         iatom = ca(ic)
         if (ta(1,it) .eq. iatom) then
            iaoth = ta(2,it)
         else
            iaoth = ta(1,it)
         end if
         do k = 1, 3
            acvect(k) = a(k,iaoth) - a(k,iatom)
            cpvect(k) = pnt(k) - c(k,ic)
         end do
         if (dot(acvect,cpvect) .ge. 0.0d0) then
            ptincy = .false.
            return
         end if
      end do
      if (cynep(icy) .le. 2) then
         ptincy = .true.
         return
      end if
      call projct (pnt,unvect,icy,iatom,spv,nedge,fail)
      if (fail) then
         ptincy = .true.
         return
      end if
      call epuclc (spv,nedge,epu)
      totang = rotang(epu,nedge,unvect)
      ptincy = (totang .gt. 0.0d0)
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine epuclc  ##
c     ##                     ##
c     #########################
c
c
      subroutine epuclc (spv,nedge,epu)
      implicit none
      integer k,ke,ke2,le,nedge
      real*8 anorm,epun,spv(3,nedge),epu(3,nedge)
c
c
c     calculate unit vectors along edges
c
      do ke = 1, nedge
c
c     get index of second edge of corner
c
         if (ke .lt. nedge) then
            ke2 = ke + 1
         else
            ke2 = 1
         end if
c
c     unit vector along edge of cycle
c
         do k = 1, 3
            epu(k,ke) = spv(k,ke2) - spv(k,ke)
         end do
         epun = anorm(epu(1,ke))
c        if (epun .le. 0.0d0)  call error ('Null Edge in Cycle')
c
c     normalize
c
         if (epun .gt. 0.0d0) then
            do k = 1, 3
               epu(k,ke) = epu(k,ke) / epun
            end do
         else
            do k = 1, 3
               epu(k,ke) = 0.0d0
            end do
         end if
      end do
c
c     vectors for null edges come from following or preceding edges
c
      do ke = 1, nedge
         if (anorm(epu(1,ke)) .le. 0.0d0) then
            le = ke - 1
            if (le .le. 0)  le = nedge
            do k = 1, 3
               epu(k,ke) = epu(k,le)
            end do
         end if
      end do
      return
      end
c
c
c     #######################
c     ##                   ##
c     ##  function rotang  ##
c     ##                   ##
c     #######################
c
c
      function rotang (epu,nedge,unvect)
      implicit none
      integer ke,nedge
      real*8 rotang,totang,dot,dt,ang
      real*8 epu(3,nedge),unvect(3),crs(3)
c
c
      totang = 0.0d0
c
c     sum angles at vertices of cycle
c
      do ke = 1, nedge
         if (ke .lt. nedge) then
            dt = dot(epu(1,ke),epu(1,ke+1))
            call vcross (epu(1,ke),epu(1,ke+1),crs)
         else
c
c     closing edge of cycle
c
            dt = dot(epu(1,ke),epu(1,1))
            call vcross (epu(1,ke),epu(1,1),crs)
         end if
         if (dt .lt. -1.0d0)  dt = -1.0d0
         if (dt .gt. 1.0d0)  dt = 1.0d0
         ang = acos(dt)
         if (dot(crs,unvect) .gt. 0.0d0)  ang = -ang
c
c     add to total for cycle
c
         totang = totang + ang
      end do
      rotang = totang
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine vcross  ##
c     ##                     ##
c     #########################
c
c
c     "vcross" finds the cross product of two vectors
c
c
      subroutine vcross (x,y,z)
      implicit none
      real*8 x(3),y(3),z(3)
c
c
      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
      return
      end
c
c
c     ####################
c     ##                ##
c     ##  function dot  ##
c     ##                ##
c     ####################
c
c
c     "dot" finds the dot product of two vectors
c
c
      function dot (x,y)
      implicit none
      real*8 dot,x(3),y(3)
c
c
      dot = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function anorm  --  find the length of a vector  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "anorm" finds the norm (length) of a vector; used as a
c     service routine by the Connolly surface area and volume
c     computation
c
c
      function anorm (x)
      implicit none
      real*8 anorm,x(3)
c
c
      anorm = x(1)**2 + x(2)**2 + x(3)**2
      if (anorm .lt. 0.0d0)  anorm = 0.0d0
      anorm = sqrt(anorm)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine vnorm  --  normalize a vector to unit length  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "vnorm" normalizes a vector to unit length; used as a
c     service routine by the Connolly surface area and volume
c     computation
c
c
      subroutine vnorm (x,xn)
      implicit none
      integer k
      real*8 x(3),xn(3),ax,anorm
c
c
      ax = anorm(x)
      do k = 1, 3
         xn(k) = x(k) / ax
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function dist2  --  distance squared between two points  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dist2" finds the distance squared between two points; used
c     as a service routine by the Connolly surface area and volume
c     computation
c
c
      function dist2 (x,y)
      implicit none
      real*8 dist2,x(3),y(3)
c
c
      dist2 = (x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function triple  --  form triple product of three vectors  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "triple" finds the triple product of three vectors; used as
c     a service routine by the Connolly surface area and volume
c     computation
c
c
      function triple (x,y,z)
      implicit none
      real*8 triple,x(3),y(3),z(3),xy(3),dot
c
c
      call vcross (x,y,xy)
      triple = dot(xy,z)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function vecang  --  finds the angle between two vectors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "vecang" finds the angle between two vectors handed with respect
c     to a coordinate axis; returns an angle in the range [0,2*pi]
c
c
      function vecang (v1,v2,axis,hand)
      implicit none
      include 'math.i'
      real*8 vecang,v1(3),v2(3),axis(3),hand
      real*8 angle,a1,a2,a12,dt
      real*8 anorm,dot,triple
c
c
      a1 = anorm(v1)
      a2 = anorm(v2)
      dt = dot(v1,v2)
      a12 = a1 * a2
      if (abs(a12) .ne. 0.0d0)  dt = dt/a12
      if (dt .lt. -1.0d0)  dt = -1.0d0
      if (dt .gt. 1.0d0)  dt = 1.0d0
      angle = acos(dt)
      if (hand*triple(v1,v2,axis) .lt. 0.0d0) then
         vecang = 2.0d0*pi - angle
      else
         vecang = angle
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine cirpln  ##
c     ##                     ##
c     #########################
c
c
      subroutine cirpln (circen,cirrad,cirvec,plncen,plnvec,
     &                        cinsp,cintp,xpnt1,xpnt2)
      implicit none
      integer k
      real*8 anorm,dot,dcp,dir,ratio,rlen
      real*8 circen(3),cirrad,cirvec(3)
      real*8 plncen(3),plnvec(3)
      real*8 xpnt1(3),xpnt2(3),cpvect(3),pnt1(3)
      real*8 vect1(3),vect2(3),uvect1(3),uvect2(3)
      logical cinsp,cintp
c
c
      do k = 1, 3
         cpvect(k) = plncen(k) - circen(k)
      end do
      dcp = dot(cpvect,plnvec)
      cinsp = (dcp .gt. 0.0d0)
      call vcross (plnvec,cirvec,vect1)
      if (anorm(vect1) .gt. 0.0d0) then
         call vnorm (vect1,uvect1)
         call vcross (cirvec,uvect1,vect2)
         if (anorm(vect2) .gt. 0.0d0) then
            call vnorm (vect2,uvect2)
            dir = dot(uvect2,plnvec)
            if (dir .ne. 0.0d0) then
               ratio = dcp / dir
               if (abs(ratio) .le. cirrad) then
                  do k = 1, 3
                     pnt1(k) = circen(k) + ratio*uvect2(k)
                  end do
                  rlen = cirrad**2 - ratio**2
                  if (rlen .lt. 0.0d0)  rlen = 0.0d0
                  rlen = sqrt(rlen)
                  do k = 1, 3
                     xpnt1(k) = pnt1(k) - rlen*uvect1(k)
                     xpnt2(k) = pnt1(k) + rlen*uvect1(k)
                  end do
                  cintp = .true.
                  return
               end if
            end if
         end if
      end if
      cintp = .false.
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine gendot  --  find surface points on unit sphere  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "gendot" finds the coordinates of a specified number of surface
c     points for a sphere with the input radius and coordinate center
c
c
      subroutine gendot (ndots,dots,radius,xcenter,ycenter,zcenter)
      implicit none
      include 'math.i'
      integer i,j,k,ndots
      integer nequat,nvert,nhoriz
      real*8 fi,fj,x,y,z,xy
      real*8 xcenter,ycenter,zcenter
      real*8 radius,dots(3,ndots)
c
c
      nequat = sqrt(pi*dble(ndots))
      nvert = 0.5d0 * nequat
      if (nvert .lt. 1)  nvert = 1
      k = 0
      do i = 0, nvert
         fi = (pi * dble(i)) / dble(nvert)
         z = cos(fi)
         xy = sin(fi)
         nhoriz = nequat * xy
         if (nhoriz .lt. 1)  nhoriz = 1
         do j = 0, nhoriz-1
            fj = (2.0d0 * pi * dble(j)) / dble(nhoriz)
            x = cos(fj) * xy
            y = sin(fj) * xy
            k = k + 1
            dots(1,k) = x*radius + xcenter
            dots(2,k) = y*radius + ycenter
            dots(3,k) = z*radius + zcenter
            if (k .ge. ndots)  goto 10
         end do
      end do
   10 continue
      ndots = k
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine error  --  surface area-volume error message  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "error" is the error handling routine for the Connolly
c     surface area and volume computation
c
c
      subroutine error (string)
      implicit none
      include 'iounit.i'
      integer i,leng,trimtext
      character*60 string
c
c
c     find the length of the message string
c
      leng = trimtext(string)
c
c     write out the error message and quit
c
      write (iout,10)  (string(i:i),i=1,leng)
   10 format (/,' CONNOLLY  --  ',60a1)
      call fatal
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine control  --  set information and output types  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "control" gets initial values for parameters that determine
c     the output style and information level provided by TINKER
c
c
      subroutine control
      implicit none
      include 'sizes.i'
      include 'inform.i'
      include 'keys.i'
      include 'output.i'
      integer i,next
      character*20 keyword
      character*80 record
c
c
c     set default values for information and output variables
c
      verbose = .false.
      debug = .false.
      abort = .false.
      savecycle = .false.
      use_version = .true.
c
c     search keywords for various control parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'VERBOSE ') then
            verbose = .true.
         else if (keyword(1:6) .eq. 'DEBUG ') then
            debug = .true.
            verbose = .true.
         else if (keyword(1:10) .eq. 'SAVECYCLE ') then
            savecycle = .true.
         else if (keyword(1:10) .eq. 'OVERWRITE ') then
            savecycle = .false.
         else if (keyword(1:10) .eq. 'NOVERSION ') then
            use_version = .false.
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cutoffs  --  set distance and Hessian cutoffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cutoffs" initializes and stores the cutoff distance windows,
c     the Hessian matrix cutoff, and the neighbor generation method
c     for pairwise potential functions
c
c
      subroutine cutoffs
      implicit none
      include 'sizes.i'
      include 'cutoff.i'
      include 'hescut.i'
      include 'keys.i'
      integer i,next
      character*20 keyword
      character*80 record,string
      logical truncate
c
c
c     set default values for cutoffs and the switching function
c
      vdwcut = 100000000.0d0
      chgcut = 100000000.0d0
      dplcut = 100000000.0d0
      vdwtaper = 0.90d0
      chgtaper = 0.65d0
      dpltaper = 0.75d0
      truncate = .false.
      hesscut = 0.0d0
c
c     neighbor default is double loop, not method of lights
c
      use_lights = .false.
c
c     search the keywords for various cutoff parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
c
c     get the cutoff radii for potential energy functions
c
         if (keyword(1:7) .eq. 'CUTOFF ') then
            string = record(next:80)
            read (string,*,err=10)  vdwcut
            chgcut = vdwcut
            dplcut = vdwcut
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            string = record(next:80)
            read (string,*,err=10)  vdwcut
         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
            string = record(next:80)
            read (string,*,err=10)  chgcut
         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
            string= record(next:80)
            read (string,*,err=10)  dplcut
c
c     get distance for initialization of energy switching
c
         else if (keyword(1:6) .eq. 'TAPER ') then
            string = record(next:80)
            read (string,*,err=10)  vdwtaper
            chgtaper = vdwtaper
            dpltaper = vdwtaper
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            string = record(next:80)
            read (string,*,err=10)  vdwtaper
         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
            string = record(next:80)
            read (string,*,err=10)  chgtaper
         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
            string= record(next:80)
            read (string,*,err=10)  dpltaper
c
c     get truncation, Hessian minimum and neighbor method
c
         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
            string = record(next:80)
            read (string,*,err=10)  hesscut
         else if (keyword(1:7) .eq. 'LIGHTS ') then
            use_lights = .true.
         end if
   10    continue
      end do
c
c     convert any distance percentages to absolute distances
c
      if (vdwtaper .lt. 1.0d0) then
         vdwtaper = vdwtaper * vdwcut
      end if
      if (chgtaper .lt. 1.0d0) then
         chgtaper = chgtaper * chgcut
      end if
      if (dpltaper .lt. 1.0d0) then
         dpltaper = dpltaper * dplcut
      end if
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper = vdwcut
         chgtaper = chgcut
         dpltaper = dplcut
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine delete  --  remove atom from coordinates list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "delete" removes a specified atom from the Cartesian
c     coordinates list and shifts the remaining atoms
c
c
      subroutine delete (iatom)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,k,m,iatom
c
c
c     reduce by one the total number of atoms
c
      n = n - 1
c
c     shift the atom coordinates, types and connectivities
c
      do i = iatom, n
         name(i) = name(i+1)
         x(i) = x(i+1)
         y(i) = y(i+1)
         z(i) = z(i+1)
         type(i) = type(i+1)
         n12(i) = n12(i+1)
         do j = 1, n12(i)
            i12(j,i) = i12(j,i+1)
         end do
      end do
c
c     remove connections to deleted atom and shift the lists
c
      do i = 1, n
         m = 0
         do j = 1, n12(i)
            if (i12(j,i) .eq. iatom) then
               m = m + 1
               do k = j, n12(i)-1
                  i12(k,i) = i12(k+1,i)
               end do
            end if
         end do
         n12(i) = n12(i) - m
         do j = 1, n12(i)
            if (i12(j,i) .gt. iatom)  i12(j,i) = i12(j,i) - 1
         end do
      end do
c
c     write a message to describe the atom deletion
c
      if (debug) then
         write (iout,10)  iatom
   10    format (' DELETE  --  Deleting Atom Number :',i8)
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine diagq  --  fast matrix diagonalization routine  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "diagq" is a matrix diagonalization routine which is a
c     conglomeration of the classical given, housec, and eigen
c     algorithms with several modifications to increase the
c     efficiency and accuracy
c
c     n         logical dimension of the matrix to be diagonalized
c     np        physical dimension of the matrix storage area
c     nv        number of eigenvalues and eigenvectors desired
c     dd        upper triangle of the matrix to be diagonalized
c     ev        returned with the eigenvalues in ascending order
c     vec       returned with the eigenvectors of the matrix
c     a,b,p,w,
c     ta,tb,y   temporary work vectors of dimension (np+1)
c
c     this routine is adapted from an original program written by
c     Bernie Brooks, DCRT, National Inst. of Health, Bethesda, MD
c
c
      subroutine diagq (n,np,nv,dd,ev,vec,a,b,p,w,ta,tb,y)
      implicit none
      include 'inform.i'
      include 'iounit.i'
      integer i,j,k,m,n,ia,ii,ji,mi,mj,mk
      integer nn,nm,np,nv,nom,nomtch,ntot
      integer ipt,iter,j1,mi1,mj1,mk1
      real*8 alimit,anorm,aroot,bx,del1,delbig,delta
      real*8 elapsed,elim1,elim2,epr,eta,expr,factor
      real*8 f0,gamma,rand1,rootl,rootx,rpow1,rpower
      real*8 s,sgn,sum1,t,temp,theta,theta1,toler
      real*8 trial,u,xkap,xnorm
      real*8 dd(np*(np+1)/2),vec(np,np),ev(np)
      real*8 a(np+1),b(np+1),p(np+1),w(np+1)
      real*8 ta(np+1),tb(np+1),y(np+1)
      logical done
c
c
c     initialization and setup of the problem
c
      if (verbose)  call setime
      eta = 1.0d-16
      theta = 1.0d37
      del1 = eta / 100.0d0
      delta = eta**2 * 100.0d0
      gamma = eta**2 / 100.0d0
      delbig = theta*delta / 1000.0d0
      theta1 = 1000.0d0 / theta
      toler = 100.0d0 * eta
      rpower = 8388608.0d0
      rpow1 = rpower * 0.5d0
      rand1 = rpower - 3.0d0
      factor = 0.0d0
      ntot = n*(n+1) / 2
      do i = 1, ntot
         factor = max(factor,abs(dd(i)))
      end do
      if (factor .eq. 0.0d0)  return
      k = 0
      anorm = 0.0d0
      do i = 1, n
         do j = i, n
            k = k + 1
            u = (dd(k)/factor)**2
            if (i .eq. j)  u = u * 0.5d0
            anorm = anorm + u
         end do
      end do
      anorm = sqrt(anorm+anorm) * factor
      do i = 1, ntot
         dd(i) = dd(i) / anorm
      end do
      if (verbose) then
         call getime (elapsed)
         write (iout,10)  elapsed
   10    format (' DIAGQ  --  Time for Initial Setup :',f13.2)
         call setime
      end if
c
c     perform the tridiagonalization step
c
      nn = n - 1
      mi = 0
      mi1 = n - 1
      do i = 1, nn
         sum1 = 0.0d0
         b(i) = 0.0d0
         ji = i + 1
         ipt = mi + i
         a(i) = dd(ipt)
         ipt = ipt + 1
         bx = dd(ipt)
         if (ji .eq. n) then
            b(i) = bx
            dd(mi+ji) = 0.0d0
            mi = mi + mi1
            mi1 = mi1 - 1
         else
            do j = ji+1, n
               ipt = ipt + 1
               sum1 = sum1 + dd(ipt)*dd(ipt)
            end do
            if (sum1 .gt. gamma) then
               s = sqrt(sum1+bx**2)
               sgn = 1.0d0
               if (bx .lt. 0.0)  sgn = -1.0d0
               temp = abs(bx)
               w(ji) = sqrt(0.5d0*(1.0d0+(temp/s)))
               ipt = mi + ji
               dd(ipt) = w(ji)
               ii = i + 2
               if (ii .le. n) then
                  temp = sgn/(2.0d0*w(ji)*s)
                  do j = ii, n
                     ipt = ipt + 1
                     w(j) = temp * dd(ipt)
                     dd(ipt) = w(j)
                  end do
               end if
               b(i) = -sgn * s
               do j = ji, n
                  p(j) = 0.0d0
               end do
               mk = mi + mi1
               mk1 = mi1 - 1
               do k = ji, n
                  ipt = mk + k
                  do m = k, n
                     bx = dd(ipt)
                     p(k) = p(k) + bx*w(m)
                     if (k .ne. m)  p(m) = p(m) + bx*w(k)
                     ipt = ipt + 1
                  end do
                  mk = mk + mk1
                  mk1 = mk1 - 1
               end do
               xkap = 0.0d0
               do k = ji, n
                  xkap = xkap + w(k)*p(k)
               end do
               do k = ji, n
                  p(k) = p(k) - xkap*w(k)
               end do
               mj = mi + mi1
               mj1 = mi1 - 1
               do j = ji, n
                  do k = j, n
                     expr = p(j)*w(k) + p(k)*w(j)
                     dd(mj+k) = dd(mj+k) - expr - expr
                  end do
                  mj = mj + mj1
                  mj1 = mj1 - 1
               end do
               mi = mi + mi1
               mi1 = mi1 - 1
            end if
         end if
      end do
      if (verbose) then
         call getime (elapsed)
         write (iout,20)  elapsed
   20    format (' DIAGQ  --  Time to Tridiagonalize :',f13.2)
         call setime
      end if
c
c     find the eigenvalues of the matrix
c
      a(n) = dd(mi+n)
      b(n) = 0.0d0
      alimit = 1.0d0
      do i = 1, n
         w(i) = b(i)
         b(i) = b(i) * b(i)
      end do
      do i = 1, nv
         ev(i) = alimit
      end do
      rootl = -alimit
      do i = 1, nv
         rootx = alimit
         do j = i, nv
            rootx = min(rootx,ev(j))
         end do
         ev(i) = rootx
   30    continue
         trial = (rootl+ev(i)) * 0.5d0
         if (trial.ne.rootl .and. trial.ne.ev(i)) then
            nomtch = n
            j = 1
   40       continue
            f0 = a(j) - trial
   50       continue
            if (abs(f0) .ge. theta1) then
               if (f0 .ge. 0.0d0)  nomtch = nomtch-1
               j = j + 1
               if (j .gt. n)  goto 60
               f0 = a(j) - trial - b(j-1)/f0
               goto 50
            end if
            j = j + 2
            nomtch = nomtch - 1
            if (j .le. n)  goto 40
   60       continue
            if (nomtch .lt. i) then
               rootl = trial
            else
               ev(i) = trial
               nom = min(nv,nomtch)
               ev(nom) = trial
            end if
            goto 30
         end if
      end do
      if (verbose) then
         call getime (elapsed)
         write (iout,70)  elapsed
   70    format (' DIAGQ  --  Time for Eigenvalues :',f15.2)
         call setime
      end if
c
c     find the eigenvectors of the matrix
c
      do i = 1, nv
         aroot = ev(i)
         do j = 1, n
            y(j) = 1.0d0
         end do
         if (i.eq.1 .or. abs(ev(i-1)-aroot).ge.toler) then
            ia = 0
         else
            ia = ia + 1
         end if
         elim1 = a(1) - aroot
         elim2 = w(1)
         do j = 1, nn
            if (abs(elim1) .gt. abs(w(j))) then
               ta(j) = elim1
               tb(j) = elim2
               p(j) = 0.0d0
               temp = w(j) / elim1
               elim1 = a(j+1) - aroot - temp*elim2
               elim2 = w(j+1)
            else
               ta(j) = w(j)
               tb(j) = a(j+1) - aroot
               p(j) = w(j+1)
               temp = 1.0d0
               if (abs(w(j)) .gt. theta1)  temp = elim1 / w(j)
               elim1 = elim2 - temp*tb(j)
               elim2 = -temp * w(j+1)
            end if
            b(j) = temp
         end do
         ta(n) = elim1
         tb(n) = 0.0d0
         p(n) = 0.0d0
         p(nn) = 0.0d0
         iter = 1
         if (ia .ne. 0)  goto 100
   80    continue
         m = n + 1
         do j = 1, n
            m = m - 1
            done = .false.
            dowhile (.not. done)
               done = .true.
               if (n-m-1 .lt. 0) then
                  elim1 = y(m)
               else
                  if (n-m-1 .eq. 0) then
                     elim1 = y(m) - y(m+1)*tb(m)
                  else
                     elim1 = y(m) - y(m+1)*tb(m) - y(m+2)*p(m)
                  end if
               end if
               if (abs(elim1) .le. delbig) then
                  temp = ta(m)
                  if (abs(temp) .lt. delta)  temp = delta
                  y(m) = elim1 / temp
               else
                  do k = 1, n
                     y(k) = y(k) / delbig
                  end do
                  done = .false.
               end if
            end do
         end do
         if (iter .eq. 2)  goto 110
         iter = iter + 1
   90    continue
         elim1 = y(1)
         do j = 1, nn
            if (ta(j) .ne. w(j)) then
               y(j) = elim1
               elim1 = y(j+1) - elim1*b(j)
            else
               y(j) = y(j+1)
               elim1 = elim1 - y(j+1)*b(j)
            end if
         end do
         y(n) = elim1
         goto 80
  100    continue
         do j = 1, n
            rand1 = mod(4099.0d0*rand1,rpower)
            y(j) = rand1/rpow1 - 1.0d0
         end do
         goto 80
  110    continue
         if (ia .ne. 0) then
            do j1 = 1, ia
               k = i - j1
               temp = 0.0d0
               do j = 1, n
                  temp = temp + y(j)*vec(j,k)
               end do
               do j = 1, n
                  y(j) = y(j) - temp*vec(j,k)
               end do
            end do
         end if
         if (iter .eq. 1)  goto 90
         elim1 = 0.0d0
         do j = 1, n
            elim1 = max(elim1,abs(y(j)))
         end do
         temp = 0.0d0
         do j = 1, n
            elim2 = y(j) / elim1
            temp = temp + elim2*elim2
         end do
         temp = 1.0d0 / (sqrt(temp)*elim1)
         do j = 1, n
            y(j) = y(j) * temp
            if (abs(y(j)) .lt. del1)  y(j) = 0.0d0
         end do
         do j = 1, n
            vec(j,i) = y(j)
         end do
      end do
      if (verbose) then
         call getime (elapsed)
         write (iout,120)  elapsed
  120    format (' DIAGQ  --  Time for Eigenvectors :',f14.2)
         call setime
      end if
c
c     perform the transformation step
c
      do i = 1, nv
         do j = 1, n
            y(j) = vec(j,i)
         end do
         mk = (n*(n-1))/2 - 3
         mk1 = 3
         nm = n - 2
         do j = 1, nm
            t = 0.0d0
            m = n - j
            do k = m, n
               t = t + dd(mk+k)*y(k)
            end do
            do k = m, n
               epr = t * dd(mk+k)
               y(k) = y(k) - 2.0d0*epr
            end do
            mk = mk - mk1
            mk1 = mk1 + 1
         end do
         t = 0.0d0
         do j = 1, n
            t = t + y(j)*y(j)
         end do
         xnorm = sqrt(t)
         do j = 1, n
            y(j) = y(j) / xnorm
         end do
         do j = 1, n
            vec(j,i) = y(j)
         end do
      end do
      do i = 1, n
         ev(i) = ev(i) * anorm
      end do
      if (verbose) then
         call getime (elapsed)
         write (iout,130)  elapsed
  130    format (' DIAGQ  --  Time for Transformation :',f12.2)
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine diffeq  --  differential equation integration  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "diffeq" performs the numerical integration of an ordinary
c     differential equation using an adaptive stepsize method to
c     solve the corresponding coupled first-order equations of the
c     general form dyi/dx = f(x,y1,...,yn) for yi = y1,...,yn
c
c     variables and parameters :
c
c     nvar      number of coupled first-order differential equations
c     y         contains the values of the dependent variables
c     x1        value of the beginning integration limit
c     x2        value of the ending integration limit
c     eps       relative accuracy required of the integration steps
c     h1        initial guess for the first integration stepsize
c     hmin      minimum allowed integration stepsize
c     nok       number of initially successful integration steps
c     nbad      number of integration steps that required retry
c
c     required external routines :
c
c     derivs    subroutine to find the right-hand side of the
c                  first-order differential equations
c
c
      subroutine diffeq (nvar,y,x1,x2,eps,h1,hmin,nok,nbad,derivs)
      include 'sizes.i'
      include 'iounit.i'
      integer maxgda,maxstep
      real*8 tiny
      parameter (maxgda=4*maxatm)
      parameter (maxstep=1000)
      parameter (tiny=1.0d-30)
      integer i,nvar,nstep,nok,nbad
      real*8 x1,x2,eps,h1,hmin
      real*8 x,h,hdid,hnext
      real*8 y(maxgda),dydx(maxgda),yscal(maxgda)
      character*7 status
      logical terminate
      external derivs
c
c
c     initialize starting limit, step size and status counters
c
      terminate = .false.
      x = x1
      h = sign(h1,x2-x1)
      nstep = 0
      nok = 0
      nbad = 0
c
c     perform a series of individual integration steps
c
      dowhile (.not. terminate)
         call derivs (x,y,dydx)
         do i = 1, nvar
            yscal(i) = abs(y(i)) + abs(h*dydx(i)) + tiny
         end do
c
c     set the final step to stop at the integration limit
c
         if ((x+h-x2)*(x+h-x1) .gt. 0.0d0)  h = x2 - x
c
c     take a Bulirsch-Stoer integration step
c
         call bsstep (nvar,x,dydx,y,h,eps,yscal,hdid,hnext,derivs)
c
c     mark the current step as either good or bad
c
         if (hdid .eq. h) then
            nok = nok + 1
            status = 'Success'
         else
            nbad = nbad + 1
            status = ' Retry '
         end if
c
c     update stepsize and get information about the current step
c
         h = hnext
         nstep = nstep + 1
         call gdastat (nstep,x,y,status)
c
c     test for convergence to the final integration limit
c
         if ((x-x2)*(x2-x1) .ge. 0.0d0) then
            write (iout,10)
   10       format (/,' DIFFEQ  --  Normal Termination',
     &                 ' at Integration Limit')
            terminate = .true.
         end if
c
c     test for a trial stepsize that is too small
c
         if (abs(hnext) .lt. hmin) then
            write (iout,20)
   20       format (/,' DIFFEQ  --  Incomplete Integration',
     &                 ' due to SmallStep')
            terminate = .true.
         end if
c
c     test for too many total integration steps
c
         if (nstep .ge. maxstep) then
            write (iout,30)
   30       format (/,' DIFFEQ  --  Incomplete Integration',
     &                 ' due to IterLimit')
            terminate = .true.
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine bsstep  --  Bulirsch-Stoer integration step  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine bsstep (nv,x,dydx,y,htry,eps,yscal,hdid,hnext,derivs)
      include 'sizes.i'
      include 'iounit.i'
      integer maxgda,kmaxx,imax
      real*8 safe1,safe2,redmax,redmin,tiny,scalmx
      parameter (maxgda=4*maxatm)
      parameter (kmaxx=8)
      parameter (imax=kmaxx+1)
      parameter (safe1=0.25d0)
      parameter (safe2=0.7d0)
      parameter (redmax=1.0d-5)
      parameter (redmin=0.7d0)
      parameter (tiny=1.0d-30)
      parameter (scalmx=0.1d0)
      integer i,iq,k,kk,km,kmax,kopt
      integer nv,nseq(imax)
      real*8 eps,hdid,hnext,htry,x
      real*8 eps1,epsold,errmax,fact,h,red
      real*8 scale,work,wrkmin,xest,xnew
      real*8 dydx(maxgda),y(maxgda),yscal(maxgda)
      real*8 a(imax),alf(kmaxx,kmaxx),err(kmaxx)
      real*8 yerr(maxgda),ysav(maxgda),yseq(maxgda)
      logical first,reduct
      save a,alf,epsold,first,kmax,kopt,nseq,xnew
      external derivs
      data first  / .true. /
      data epsold / -1.0d0 /
      data nseq   / 2,4,6,8,10,12,14,16,18 /
c
c
      if (eps .ne. epsold) then
         hnext = -1.0d29
         xnew = -1.0d29
         eps1 = safe1 * eps
         a(1) = 1.0d0 + dble(nseq(1))
         do k = 1, kmaxx
            a(k+1) = a(k) + dble(nseq(k+1))
         end do
         do iq = 2, kmaxx
            do k = 1, iq-1
               alf(k,iq) = eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.0d0)
     &                                                 *(2*k+1)))
            end do
         end do
         epsold = eps
         do kopt = 2, kmaxx-1
            if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt))  goto 10
         end do
   10    continue
         kmax = kopt
      end if
      h = htry
      do i = 1, nv
         ysav(i) = y(i)
      end do
      if (h.ne.hnext .or. x.ne.xnew) then
         first = .true.
         kopt = kmax
      end if
      reduct = .false.
   20 continue
      do k = 1, kmax
         xnew = x + h
         if (xnew .eq. x) then
            write (iout,30)
   30       format (' BSSTEP  --  Underflow of Step Size')
            call fatal
         end if
         call mmid (nseq(k),h,nv,x,dydx,ysav,yseq,derivs)
         xest = (h/dble(nseq(k)))**2
         call pzextr (k,nv,xest,yseq,y,yerr)
         if (k .ne. 1) then
            errmax = tiny
            do i = 1, nv
               errmax = max(errmax,abs(yerr(i)/yscal(i)))
            end do
            errmax = errmax / eps
            km = k - 1
            err(km) = (errmax/safe1)**(1.0d0/(2*km+1))
         end if
         if (k.ne.1 .and. (k.ge.kopt-1 .or. first)) then
            if (errmax .lt. 1.0d0)  goto 50
            if (k.eq.kmax .or. k.eq.kopt+1) then
               red = safe2 / err(km)
               goto 40
            else if (k .eq. kopt) then
               if (alf(kopt-1,kopt) .lt. err(km)) then
                  red = 1.0d0 / err(km)
                  goto 40
               end if
            else if (kopt .eq. kmax)then
               if (alf(km,kmax-1) .lt. err(km)) then
                  red = alf(km,kmax-1) * safe2 / err(km)
                  goto 40
               endif
            else if (alf(km,kopt) .lt. err(km)) then
               red = alf(km,kopt-1) / err(km)
               goto 40
            end if
         end if
      end do
   40 continue
      red = min(red,redmin)
      red = max(red,redmax)
      h = h * red
      reduct = .true.
      goto 20
   50 continue
      x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.0d35
      do kk = 1, km
         fact = max(err(kk),scalmx)
         work = fact * a(kk+1)
         if (work .lt. wrkmin) then
            scale = fact
            wrkmin = work
            kopt = kk + 1
         end if
      end do
      hnext = h / scale
      if (kopt.ge.k .and. kopt.ne.kmax .and. .not.reduct) then
         fact = max(scale/alf(kopt-1,kopt),scalmx)
         if (a(kopt+1)*fact .le. wrkmin) then
            hnext = h / fact
            kopt = kopt + 1
         end if
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine mmid  --  takes a modified midpoint step  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "mmid" implements a modified midpoint method to advance the
c     integration of a set of first order differential equations
c
c
      subroutine mmid (nstep,htot,nvar,xs,dydx,y,yout,derivs)
      include 'sizes.i'
      integer maxgda
      parameter (maxgda=4*maxatm)
      integer i,k,nstep,nvar
      real*8 xs,x,htot,h,h2,temp
      real*8 y(maxgda),yout(maxgda),dydx(maxgda)
      real*8 ym(maxgda),yn(maxgda)
      external derivs
c
c
c     set substep size based on number of steps to be taken
c
      h = htot / dble(nstep)
      h2 = 2.0d0 * h
c
c     take the first substep and get values at ends of step
c
      do i = 1, nvar
         ym(i) = y(i)
         yn(i) = y(i) + h*dydx(i)
      end do
      x = xs + h
      call derivs (x,yn,yout)
c
c     take the second and subsequent substeps
c
      do k = 2, nstep
         do i = 1, nvar
            temp = ym(i) + h2*yout(i)
            ym(i) = yn(i)
            yn(i) = temp
         end do
         x = x + h
         call derivs (x,yn,yout)
      end do
c
c     complete the update of values for the last substep
c
      do i = 1, nvar
         yout(i) = 0.5d0 * (ym(i)+yn(i)+h*yout(i))
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine pzextr  --  polynomial extrapolation method  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine pzextr (iest,nvar,xest,yest,yz,dy)
      include 'sizes.i'
      integer maxgda,imax
      parameter (maxgda=4*maxatm)
      parameter (imax=13)
      integer i,j,iest,nvar
      real*8 xest,delta,f1,f2,q
      real*8 yz(maxgda),dy(maxgda)
      real*8 yest(maxgda),d(maxgda)
      real*8 x(imax),qcol(maxgda,imax)
      save x,qcol
c
c
      x(iest) = xest
      do j = 1, nvar
         dy(j) = yest(j)
         yz(j) = yest(j)
      end do
      if (iest .eq. 1) then
         do j = 1, nvar
            qcol(j,1) = yest(j)
         end do
      else
         do j = 1, nvar
            d(j) = yest(j)
         end do
         do i = 1, iest-1
            delta = 1.0d0 / (x(iest-i)-xest)
            f1 = xest * delta
            f2 = x(iest-i) * delta
            do j = 1, nvar
               q = qcol(j,i)
               qcol(j,i) = dy(j)
               delta = d(j) - q
               dy(j) = f1 * delta
               d(j) = f2 * delta
               yz(j) = yz(j) + dy(j)
            end do
         end do
         do j = 1, nvar
            qcol(j,iest) = dy(j)
         end do
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  function dihedral  --  torsion angle between four atoms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dihedral" finds the value of the dihedral angle in the
c     range from -180 to +180 degrees defined by four input atoms
c
c
      function dihedral (ia,ib,ic,id)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'math.i'
      integer ia,ib,ic,id
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 rt2,ru2,rtru
      real*8 dihedral,cosine,sign
c
c
c     set default in case atoms are colinear or degenerate
c
      dihedral = 0.0d0
c
c     compute the value in degrees of the dihedral angle
c
      xba = x(ib) - x(ia)
      yba = y(ib) - y(ia)
      zba = z(ib) - z(ia)
      xcb = x(ic) - x(ib)
      ycb = y(ic) - y(ib)
      zcb = z(ic) - z(ib)
      xdc = x(id) - x(ic)
      ydc = y(id) - y(ic)
      zdc = z(id) - z(ic)
      xt = yba*zcb - ycb*zba
      yt = xcb*zba - xba*zcb
      zt = xba*ycb - xcb*yba
      xu = ycb*zdc - ydc*zcb
      yu = xdc*zcb - xcb*zdc
      zu = xcb*ydc - xdc*ycb
      rt2 = xt*xt + yt*yt + zt*zt
      ru2 = xu*xu + yu*yu + zu*zu
      rtru = sqrt(rt2 * ru2)
      if (rtru .ne. 0.0d0) then
         cosine = (xt*xu + yt*yu + zt*zu) / rtru
         cosine = min(1.0d0,max(-1.0d0,cosine))
         dihedral = radian * acos(cosine)
         sign = xba*xu + yba*yu + zba*zu
         if (sign .lt. 0.0d0)  dihedral = -dihedral
      end if
      return
      end
