c   TINKER library routines that start with em to eo
c                  T(inker)lib(rary)emo.f
c
C 21 Apr 10 - KK  - allow choosing minimizer
c  8 MAY 98 - JRS call rotmat changed to call rotmatt
c                 call extra changed to call extrat
c
c DGF thinks that maxgeo usage is very dangerous because maxgeo can be smaller
c than the total number of atoms and then array overflow will occur. 
c To uncomment, use "%s/^cf/  /g" in vi.
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine embed  --  structures via distance geometry  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "embed" is a distance geometry routine patterned after the
c     ideas of Gordon Crippen, Irwin Kuntz and Tim Havel; it takes
c     as input a set of upper and lower bounds on the interpoint
c     distances, chirality constraints and torsional constraints,
c     and attempts to generate a set of coordinates that satisfy
c     the input bounds and constraints
c
c     literature references:
c
c     G. M. Crippen and T. F. Havel, "Distance Geometry and Molecular
c     Conformation", Research Studies Press, Letchworth U.K., 1988,
c     John Wiley and Sons, U.S. distributor
c
c     T. F. Havel, "An Evaluation of Computational Strategies for
c     Use in the Determination of Protein Structure from Distance
c     Constraints obtained by Nuclear Magnetic Resonance", Progress
c     in Biophysics and Molecular Biology, 56, 43-78 (1991)
c
c
cf    subroutine embed
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'files.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    include 'minima.i'
cf    include 'refer.i'
cf    integer maxeigen
cf    parameter (maxeigen=5)
cf    integer i,j,nstep,lext
cf    integer igeo,freeunit
cf    integer maxneg,nneg
cf    integer maxinner,ninner
cf    integer maxouter,nouter
cf    real*8 fctval,grdmin,secs
cf    real*8 dt,temp_start,temp_stop
cf    real*8 rg,rmsorig,rmsflip,mass
cf    real*8 bounds,contact,local,chiral,torsion
cf    real*8 bnderr,vdwerr,locerr,chirer,torser
cf    real*8 derivs(3,maxgeo),matrix(maxgeo,maxgeo)
cf    real*8 evl(maxeigen),evc(maxgeo,maxeigen)
cf    real*8 v(maxvar),a(maxvar)
cf    character*7 errtyp,ext
cf    character*60 title,geofile
cf    logical done,valid,exist,info
ckk -minimize- begin
cf    character*3 method
cf    method='   '
ckk -minimize- end
c
c
c     initialize any chirality constraints, then smooth the
c     bounds via triangle and inverse triangle inequalities;
c     currently these functions are performed by "distgeom"
c
c     call chirin
c     if (debug .or. (verbose .and. n.le.130)) then
c        title = 'Input Distance Bounds :'
c        call grafic (n,maxgeo,bnd,title)
c     end if
c     call geodesic
c     if (debug .or. (verbose .and. n.le.130)) then
c        title = 'Triangle Smoothed Bounds :'
c        call grafic (n,maxgeo,bnd,title)
c     end if
c
c     generate a distance matrix between the upper and
c     lower bounds, then convert to a metric matrix
c
cf    maxinner = 3
cf    maxouter = 3
cf    maxneg = 2
cf    nouter = 0
cf    valid = .false.
cf    dowhile (.not. valid)
cf       ninner = 0
cf       done = .false.
cf       dowhile (.not. done)
cf          ninner = ninner + 1
cf          call dstmat (matrix)
cf          call metric (matrix,nneg)
cf          if (nneg.le.maxneg .or. ninner.eq.maxinner)  done = .true.
cf          compact = 0.0d0
cf       end do
cf       if (verbose .and. nneg.gt.maxneg) then
cf          write (iout,10)  nneg
cf 10       format (/,' EMBED  --  Warning, Using Metric Matrix',
cf   &                  ' with',i4,' Negative Distances')
cf       end if
c
c     find the principle components of metric matrix, then
c     generate the trial Cartesian coordinates
c
cf       nouter = nouter + 1
cf       call eigen (evl,evc,matrix,valid)
cf       call coords (evl,evc)
cf       if (nouter.eq.maxouter .and. .not.valid) then
cf          valid = .true.
cf          if (verbose) then
cf             write (iout,20)
cf 20          format (/,' EMBED  --  Warning, Using Poor Initial',
cf   &                    ' Coordinates')
cf          end if
cf       end if
cf    end do
c
c     compute an index of compaction from the embedded structure
c
cf    call chksize
c
c     use majorization to improve the approximate coordinates
c
cf    do i = 1, n
cf       matrix(i,i) = 0.0d0
cf       do j = i+1, n
cf          matrix(j,i) = matrix(i,j)
cf       end do
cf    end do
cf    call majorize (matrix)
c
c     superimpose embedded structure or enantiomer on reference
c
cf    info = verbose
cf    verbose = .false.
cf    call impose (nref,xref,yref,zref,n,x,y,z,rmsorig)
cf    if (use_invert) then
cf       do i = 1, n
cf          x(i) = -x(i)
cf       end do
cf       call impose (nref,xref,yref,zref,n,x,y,z,rmsflip)
cf       if (rmsorig .lt. rmsflip) then
cf          do i = 1, n
cf             x(i) = -x(i)
cf          end do
cf          call impose (nref,xref,yref,zref,n,x,y,z,rmsorig)
cf       end if
cf       write (iout,30)  rmsorig,rmsflip
cf 30    format (/,' RMS Superposition for Original and',
cf   &              ' Enantiomer : ',2f12.4)
cf    end if
cf    verbose = info
c
c     write out the embedded unrefined atomic coordinate set
c
cf    if (debug) then
cf       i = 0
cf       exist = .true.
cf       dowhile (exist)
cf          i = i + 1
cf          lext = 3
cf          call numeral (i,ext,lext)
cf          geofile = filename(1:leng)//'-embed'//'.'//ext(1:lext)
cf          inquire (file=geofile,exist=exist)
cf       end do
cf       igeo = freeunit ()
cf       open (unit=igeo,file=geofile,status='new')
cf       call prtxyz (igeo)
cf       close (unit=igeo)
cf    end if
c
c     square the bounds for use during structure refinement
c
cf    do i = 1, n
cf       do j = 1, n
cf          bnd(j,i) = bnd(j,i)**2
cf       end do
cf    end do
c
c     minimize the error function via simulated annealing
c
cf    if (verbose)  call setime
cf    if (use_anneal) then
cf       iprint = 0
cf       if (verbose)  iprint = 10
cf       grdmin = 1.0d0
cf       mass = 10000.0d0
cf       do i = 1, 3*n
cf          v(i) = 0.0d0
cf          a(i) = 0.0d0
cf       end do
cf       errtyp = 'final'
ckk -minimize-         call refine (errtyp,fctval,grdmin)
cf       call refine (method,errtyp,fctval,grdmin)
cf       nstep = 1000
cf       dt = 4.0d-14
cf       temp_start = 200.0d0
cf       temp_stop = 200.0d0
cf       call explore (errtyp,nstep,dt,mass,temp_start,temp_stop,v,a)
cf       nstep = 10000
cf       dt = 2.0d-13
cf       temp_start = 200.0d0
cf       temp_stop = 0.0d0
cf       call explore (errtyp,nstep,dt,mass,temp_start,temp_stop,v,a)
cf       grdmin = 0.01d0
ckk -minimize-         call refine (errtyp,fctval,grdmin)
cf       call refine (method,errtyp,fctval,grdmin)
c
c     minimize the error function via nonlinear optimization
c
cf    else
cf       iprint = 0
cf       if (verbose)  iprint = 10
cf       grdmin = 0.01d0
cf       errtyp = 'initial'
ckk -minimize-         call refine (errtyp,fctval,grdmin)
cf       call refine (method,errtyp,fctval,grdmin)
cf       errtyp = 'middle'
ckk -minimize-         call refine (errtyp,fctval,grdmin)
cf       call refine (method,errtyp,fctval,grdmin)
cf       errtyp = 'final'
ckk -minimize-         call refine (errtyp,fctval,grdmin)
cf       call refine (method,errtyp,fctval,grdmin)
cf    end if
cf    if (verbose) then
cf       call getime (secs)
cf       write (iout,40)  secs
cf 40    format (/,' Time Required for Refinement :',10x,f12.2,
cf   &              ' seconds')
cf    end if
c
c     print the final error function and its components
c
cf    bounds = bnderr (derivs)
cf    contact = vdwerr (derivs)
cf    local = locerr (derivs)
cf    chiral = chirer (derivs)
cf    torsion = torser (derivs)
cf    write (iout,50)  fctval,bounds,contact,local,chiral,torsion
cf 50 format (/,' Results of Distance Geometry Embedding :',
cf   &        //,' Final Error Function Value :',10x,f16.4,
cf   &        //,' Distance Restraint Error :',12x,f16.4,
cf   &        /,' Hard Sphere Contact Error :',11x,f16.4,
cf   &        /,' Local Geometry Error :',16x,f16.4,
cf   &        /,' Chirality-Planarity Error :',11x,f16.4,
cf   &        /,' Torsional Restraint Error :',11x,f16.4)
c
c     take the root of the currently squared distance bounds
c
cf    do i = 1, n
cf       do j = 1, n
cf          bnd(j,i) = sqrt(bnd(j,i))
cf       end do
cf    end do
c
c     print the final rms deviations and radius of gyration
c
cf    title = 'after Refinement Protocol :'
cf    call rmserror (title)
cf    call gyrate (rg)
cf    write (iout,60)  rg
cf 60 format (/,' Radius of Gyration of the System :',10x,f16.4)
cf    if (debug .or. (verbose .and. n.le.130)) then
cf       call dmdump (matrix)
cf    end if
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine chirin  --  initializes chirality constraints  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "chirin" determines the target value for each chirality
c     and planarity constraint as the signed volume of the
c     parallelpiped spanned by vectors from a common atom to
c     each of three other atoms
c
c
cf    subroutine chirin
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    integer i,j,ia,ib,ic,id
cf    real*8 xad,yad,zad
cf    real*8 xbd,ybd,zbd
cf    real*8 xcd,ycd,zcd
cf    real*8 c1,c2,c3
c
c
c     compute the signed volume of each parallelpiped;
c     if the defining atoms almost lie in a plane, then
c     set the signed volume to exactly zero
c
cf    do i = 1, nchir
cf       ia = ichir(1,i)
cf       ib = ichir(2,i)
cf       ic = ichir(3,i)
cf       id = ichir(4,i)
cf       xad = x(ia) - x(id)
cf       yad = y(ia) - y(id)
cf       zad = z(ia) - z(id)
cf       xbd = x(ib) - x(id)
cf       ybd = y(ib) - y(id)
cf       zbd = z(ib) - z(id)
cf       xcd = x(ic) - x(id)
cf       ycd = y(ic) - y(id)
cf       zcd = z(ic) - z(id)
cf       c1 = ybd*zcd - zbd*ycd
cf       c2 = ycd*zad - zcd*yad
cf       c3 = yad*zbd - zad*ybd
cf       vchir(i) = xad*c1 + xbd*c2 + xcd*c3
cf       if (abs(vchir(i)) .lt. 1.0d0)  vchir(i) = 0.0d0
cf    end do
c
c     print out the results for each constraint
c
cf    if (verbose) then
cf       if (nchir .ne. 0) then
cf          write (iout,10)
cf 10       format (/,' Chirality and Planarity Constraints :')
cf          write (iout,20)
cf 20       format (/,18x,'Atom Numbers',12x,'Signed Volume',/)
cf       end if
cf       do i = 1, nchir
cf          write (iout,30)  i,(ichir(j,i),j=1,4),vchir(i)
cf 30       format (i6,5x,4i6,5x,f12.4)
cf       end do
cf    end if
cf    return
cf    end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine triangle  --  triangle inequality smoothing  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "triangle" smooths the upper and lower distance bounds via
c     the triangle inequality using a full-matrix variant of the
c     Floyd-Warshall shortest path algorithm; this routine is
c     usually much slower than the sparse matrix shortest path
c     methods in "geodesic" and "trifix", and should be used only
c     for comparison with answers generated by those routines
c
c     literature reference:
c
c     A. W. M. Dress and T. F. Havel, "Shortest-Path Problems and
c     Molecular Conformation", Discrete Applied Mathematics, 19,
c     129-144 (1988)
c
c
cf    subroutine triangle
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'iounit.i'
cf    integer i,j,k
cf    integer ik1,ik2
cf    integer jk1,jk2
cf    real*8 eps
cf    real*8 lij,lik,ljk
cf    real*8 uij,uik,ujk
c
c
c     use full-matrix algorithm to smooth upper and lower bounds
c
cf    eps = 1.0d-10
cf    do k = 1, n
cf       do i = 1, n-1
cf          ik1 = min(i,k)
cf          ik2 = max(i,k)
cf          lik = bnd(ik2,ik1)
cf          uik = bnd(ik1,ik2)
cf          do j = i+1, n
cf             lij = bnd(j,i)
cf             uij = bnd(i,j)
cf             jk1 = min(j,k)
cf             jk2 = max(j,k)
cf             ljk = bnd(jk2,jk1)
cf             ujk = bnd(jk1,jk2)
cf             lij = max(lij,lik-ujk,ljk-uik)
cf             uij = min(uij,uik+ujk)
cf             if (lij-uij .gt. eps) then
cf                write (iout,10)
cf 10             format (/,' TRIANGLE  --  Inconsistent Bounds;',
cf   &                       ' Geometrically Impossible')
cf                write (iout,20)  i,j,lij,uij
cf 20             format (/,' Error at :',6x,2i6,3x,2f9.4)
cf                write (iout,30)  i,k,lik,uik,j,k,ljk,ujk
cf 30             format (/,' Traced to :',5x,2i6,3x,2f9.4,
cf   &                    /,17x,2i6,3x,2f9.4)
cf                call fatal
cf             end if
cf             if (lij-bnd(j,i) .gt. eps) then
cf                write (iout,40)  i,j,bnd(j,i),lij
cf 40             format (' TRIANGLE  --  Altered Lower Bound at',
cf   &                       2x,2i6,3x,f9.4,' -->',f9.4)
cf             end if
cf             if (bnd(i,j)-uij .gt. eps) then
cf                write (iout,50)  i,j,bnd(i,j),uij
cf 50             format (' TRIANGLE  --  Altered Upper Bound at',
cf   &                       2x,2i6,3x,f9.4,' -->',f9.4)
cf             end if
cf             bnd(j,i) = lij
cf             bnd(i,j) = uij
cf          end do
cf       end do
cf    end do
cf    return
cf    end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine geodesic  --  sparse matrix triangle smoothing  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "geodesic" smooths the upper and lower distance bounds via
c     the triangle inequality using a sparse matrix version of a
c     shortest path algorithm
c
c     literature reference:
c
c     G. M. Crippen and T. F. Havel, "Distance Geometry and Molecular
c     Conformation", Research Studies Press, Letchworth U.K., 1988,
c     John Wiley and Sons, U.S. distributor, see section 6-2
c
c
cf    subroutine geodesic
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'restrn.i'
cf    integer maxlist
cf    parameter (maxlist=2*maxfix)
cf    integer i,j,k,nlist
cf    integer start(maxatm),stop(maxatm)
cf    integer list(maxlist),key(maxlist)
cf    real*8 upper(maxatm),lower(maxatm)
c
c
c     build an indexed list of atoms in distance restraints
c
cf    do i = 1, n
cf       start(i) = 0
cf       stop(i) = -1
cf    end do
cf    nlist = 2 * ndfix
cf    do i = 1, ndfix
cf       list(i) = idfix(1,i)
cf       list(i+ndfix) = idfix(2,i)
cf    end do
cf    call sort3 (nlist,list,key)
cf    j = -1
cf    do i = 1, nlist
cf       k = list(i)
cf       if (k .ne. j) then
cf          start(k) = i
cf          j = k
cf       end if
cf    end do
cf    j = -1
cf    do i = nlist, 1, -1
cf       k = list(i)
cf       if (k .ne. j) then
cf          stop(k) = i
cf          j = k
cf       end if
cf    end do
cf    do i = 1, nlist
cf       k = key(i)
cf       if (k .le. ndfix) then
cf          list(i) = idfix(2,k)
cf       else
cf          list(i) = idfix(1,k-ndfix)
cf       end if
cf    end do
c
c     triangle smooth bounds via sparse shortest path method
c
cf    do i = 1, n
cf       call minpath (i,upper,lower,start,stop,list)
cf       do j = i+1, n
cf          bnd(i,j) = upper(j)
cf          bnd(j,i) = max(lower(j),bnd(j,i))
cf       end do
cf    end do
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine minpath  --  triangle smoothed bounds to atom  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "minpath" is a routine for finding the triangle smoothed upper
c     and lower bounds of each atom to a specified root atom using a
c     sparse variant of the Bellman-Ford shortest path algorithm
c
c     literature reference:
c
c     D. P. Bertsekas, "A Simple and Fast Label Correcting Algorithm
c     for Shortest Paths", Networks, 23, 703-709 (1993)
c
c
cf    subroutine minpath (root,upper,lower,start,stop,list)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'couple.i'
cf    include 'disgeo.i'
cf    integer maxlist
cf    parameter (maxlist=2*maxfix)
cf    integer i,j,k,root
cf    integer narc,iarc(maxatm)
cf    integer head,tail,queue(maxatm)
cf    integer start(maxatm),stop(maxatm)
cf    integer list(maxlist)
cf    real*8 big,upper(maxatm)
cf    real*8 small,lower(maxatm)
cf    logical enter,queued(maxatm)
c
c
c     initialize candidate atom queue and the path lengths
c
cf    do i = 1, n
cf       queued(i) = .false.
cf       upper(i) = 1000000.0d0
cf       lower(i) = 0.0d0
cf    end do
c
c     put the root atom into the queue of candidate atoms
c
cf    head = root
cf    tail = root
cf    queue(root) = 0
cf    queued(root) = .true.
cf    upper(root) = 0.0d0
c
c     get the next candidate atom from head of queue
c
cf    dowhile (head .ne. 0)
cf       j = head
cf       queued(j) = .false.
cf       head = queue(head)
c
c     make a list of arcs to the current candidate atom
c
cf       narc = 0
cf       do i = 1, n12(j)
cf          k = i12(i,j)
cf          if (k .ne. root) then
cf             narc = narc + 1
cf             iarc(narc) = k
cf          end if
cf       end do
cf       do i = 1, n13(j)
cf          k = i13(i,j)
cf          if (k .ne. root) then
cf             narc = narc + 1
cf             iarc(narc) = k
cf          end if
cf       end do
cf       do i = 1, n14(j)
cf          k = i14(i,j)
cf          if (k .ne. root) then
cf             narc = narc + 1
cf             iarc(narc) = k
cf          end if
cf       end do
cf       do i = start(j), stop(j)
cf          k = list(i)
cf          if (k .ne. root) then
cf             narc = narc + 1
cf             iarc(narc) = k
cf          end if
cf       end do
c
c     check each arc for alteration of the path length bounds
c
cf       do i = 1, narc
cf          k = iarc(i)
cf          if (k .lt. j) then
cf             big = upper(j) + bnd(k,j)
cf             small = max(bnd(j,k)-upper(j),lower(j)-bnd(k,j))
cf          else
cf             big = upper(j) + bnd(j,k)
cf             small = max(bnd(k,j)-upper(j),lower(j)-bnd(j,k))
cf          end if
cf          enter = .false.
cf          if (upper(k) .gt. big) then
cf             upper(k) = big
cf             if (.not. queued(k))  enter = .true.
cf          end if
cf          if (lower(k) .lt. small) then
cf             lower(k) = small
cf             if (.not. queued(k))  enter = .true.
cf          end if
c
c     enter a new candidate atom at the tail of the queue
c
cf          if (enter) then
cf             queued(k) = .true.
cf             if (head .eq. 0) then
cf                head = k
cf                tail = k
cf                queue(k) = 0
cf             else
cf                queue(tail) = k
cf                queue(k) = 0
cf                tail = k
cf             end if
cf          end if
cf       end do
cf    end do
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine trifix  --  update triangle inequality bounds  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "trifix" rebuilds both the upper and lower distance bound
c     matrices following tightening of one or both of the bounds
c     between a specified pair of atoms, "p" and "q", using a
c     modification of Murchland's shortest path update algorithm
c
c     literature references:
c
c     P. A. Steenbrink, "Optimization of Transport Networks", John
c     Wiley and Sons, Bristol, 1974; see section 7.7
c
c     R. Dionne, "Etude et Extension d'un Algorithme de Murchland",
c     Infor, 16, 132-146 (1978)
c
c
cf    subroutine trifix (p,q)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    integer i,k,p,q,ip,iq,np,nq
cf    integer pt(maxgeo),qt(maxgeo)
cf    real*8 eps,ipmin,ipmax,iqmin,iqmax
cf    real*8 pmin(maxgeo),pmax(maxgeo)
cf    real*8 qmin(maxgeo),qmax(maxgeo)
cf    logical pun(maxgeo),qun(maxgeo)
c
c
c     initialize the set of nodes that may have changed bounds
c
cf    eps = 1.0d-10
cf    np = 0
cf    nq = 0
cf    do i = 1, n
cf       pun(i) = .true.
cf       qun(i) = .true.
cf    end do
c
c     store the upper and lower bounds to "p" and "q"
c
cf    do i = 1, p
cf       pmin(i) = bnd(p,i)
cf       pmax(i) = bnd(i,p)
cf    end do
cf    do i = p+1, n
cf       pmin(i) = bnd(i,p)
cf       pmax(i) = bnd(p,i)
cf    end do
cf    do i = 1, q
cf       qmin(i) = bnd(q,i)
cf       qmax(i) = bnd(i,q)
cf    end do
cf    do i = q+1, n
cf       qmin(i) = bnd(i,q)
cf       qmax(i) = bnd(q,i)
cf    end do
c
c     check for changes in the upper bounds to "p" and "q"
c
cf    do i = 1, n
cf       ipmax = qmax(p) + qmax(i)
cf       if (pmax(i) .gt. ipmax+eps) then
cf          np = np + 1
cf          pt(np) = i
cf          pmax(i) = ipmax
cf          pun(i) = .false.
cf       end if
cf       iqmax = pmax(q) + pmax(i)
cf       if (qmax(i) .gt. iqmax+eps) then
cf          nq = nq + 1
cf          qt(nq) = i
cf          qmax(i) = iqmax
cf          qun(i) = .false.
cf       end if
cf    end do
c
c     for node pairs whose bounds to "p" and "q" have changed,
c     make any needed changes to upper bound of the pair
c
cf    do ip = 1, np
cf       i = pt(ip)
cf       ipmax = pmax(i)
cf       do iq = 1, nq
cf          k = qt(iq)
cf          if (i .lt. k) then
cf             bnd(i,k) = min(bnd(i,k),ipmax+pmax(k))
cf          else
cf             bnd(k,i) = min(bnd(k,i),ipmax+pmax(k))
cf          end if
cf       end do
cf    end do
c
c     check for changes in the lower bounds to "p" and "q"
c
cf    do i = 1, n
cf       ipmin = max(qmin(p)-qmax(i),qmin(i)-qmax(p))
cf       if (pmin(i) .lt. ipmin-eps) then
cf          if (pun(i)) then
cf             np = np + 1
cf             pt(np) = i
cf          end if
cf          pmin(i) = ipmin
cf       end if
cf       iqmin = max(pmin(q)-pmax(i),pmin(i)-pmax(q))
cf       if (qmin(i) .lt. iqmin-eps) then
cf          if (qun(i)) then
cf             nq = nq + 1
cf             qt(nq) = i
cf          end if
cf          qmin(i) = iqmin
cf       end if
cf    end do
c
c     for node pairs whose bounds to "p" and "q" have changed,
c     make any needed changes to lower bound of the pair
c
cf    do ip = 1, np
cf       i = pt(ip)
cf       ipmin = pmin(i)
cf       ipmax = pmax(i)
cf       do iq = 1, nq
cf          k = qt(iq)
cf          if (i .lt. k) then
cf             bnd(k,i) = max(bnd(k,i),ipmin-pmax(k),pmin(k)-ipmax)
cf          else
cf             bnd(i,k) = max(bnd(i,k),ipmin-pmax(k),pmin(k)-ipmax)
cf          end if
cf       end do
cf    end do
c
c     update the upper and lower bounds to "p" and "q"
c
cf    do i = 1, p
cf       bnd(p,i) = pmin(i)
cf       bnd(i,p) = pmax(i)
cf    end do
cf    do i = p+1, n
cf       bnd(i,p) = pmin(i)
cf       bnd(p,i) = pmax(i)
cf    end do
cf    do i = 1, q
cf       bnd(q,i) = qmin(i)
cf       bnd(i,q) = qmax(i)
cf    end do
cf    do i = q+1, n
cf       bnd(i,q) = qmin(i)
cf       bnd(q,i) = qmax(i)
cf    end do
c
c     output the atoms updated and amount of work required
c
cf    if (debug) then
cf       write (iout,10)  p,q,np*nq
cf 10    format (' TRIFIX  --  Bounds Update for Atoms',2i6,
cf   &              ' with',i8,' Searches')
cf    end if
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine grafic  --  schematic graphical matrix output  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "grafic" outputs the upper & lower triangles and diagonal
c     of a square matrix in a schematic form for visual inspection
c
c
cf    subroutine grafic (n,np,a,title)
cf    implicit none
cf    include 'iounit.i'
cf    integer i,j,k,m,n,np
cf    integer maxj,nrow,ndash
cf    integer minrow,maxrow
cf    real*8 a(np,np),big
cf    real*8 amin,dmin,bmin
cf    real*8 amax,dmax,bmax
cf    real*8 v,rcl,scl,tcl
cf    real*8 ca,cb,cc,cd
cf    real*8 cw,cx,cy,cz
cf    character*1 ta,tb,tc,td,te
cf    character*1 dash,digit(0:9)
cf    character*1 symbol(130)
cf    character*60 title
cf    data ta,tb,tc,td,te  / ' ','.','+','X','#' /
cf    data dash   / '-' /
cf    data digit  / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     set bounds of length of print row and write the header
c
cf    minrow = 54
cf    maxrow = 130
cf    ndash = min(max(n,minrow),maxrow)
cf    write (iout,10)  (dash,i=1,ndash)
cf 10 format (/,1x,130a1)
cf    write (iout,20)  title
cf 20 format (/,1x,a60)
c
c     find the maximum and minimum elements of the matrix
c
cf    big = 1.0d6
cf    dmax = -big
cf    dmin = big
cf    amax = -big
cf    amin = big
cf    bmax = -big
cf    bmin = big
cf    do i = 1, n
cf       if (a(i,i) .gt. dmax)  dmax = a(i,i)
cf       if (a(i,i) .lt. dmin)  dmin = a(i,i)
cf       do j = 1, i-1
cf          if (a(j,i) .gt. amax)  amax = a(j,i)
cf          if (a(j,i) .lt. amin)  amin = a(j,i)
cf          if (a(i,j) .gt. bmax)  bmax = a(i,j)
cf          if (a(i,j) .lt. bmin)  bmin = a(i,j)
cf       end do
cf    end do
cf    write (iout,30)  amin,amax,dmin,dmax,bmin,bmax
cf 30 format (/,' Range of Above Diag Elements : ',f13.4,' to',f13.4,
cf   &        /,' Range of Diagonal Elements :   ',f13.4,' to',f13.4,
cf   &        /,' Range of Below Diag Elements : ',f13.4,' to',f13.4)
c
c     now, print out the graphical representation
c
cf    write (iout,40)
cf 40 format (/,' Symbol Magnitude Ordering :',14x,
cf   &          '# > X > + > . > '' ''',/)
cf    rcl = (bmax-bmin) / 5.0d0
cf    scl = (amax-amin) / 5.0d0
cf    tcl = (dmax-dmin) / 9.0d0
cf    if (rcl .eq. 0.0d0)  rcl = 1.0d0
cf    if (scl .eq. 0.0d0)  scl = 1.0d0
cf    if (tcl .eq. 0.0d0)  tcl = 1.0d0
cf    ca = amin + scl
cf    cb = ca + scl
cf    cc = cb + scl
cf    cd = cc + scl
cf    cw = bmin + rcl
cf    cx = cw + rcl
cf    cy = cx + rcl
cf    cz = cy + rcl
cf    do j = 1, n, maxrow
cf       maxj = j + maxrow - 1
cf       if (maxj .gt. n)  maxj = n
cf       nrow = maxj - j + 1
cf       do i = 1, n
cf          do k = j, maxj
cf             m = k - j + 1
cf             if (k .lt. i) then
cf                v = abs(a(i,k))
cf                if (v .le. cw) then
cf                   symbol(m) = ta
cf                else if (v .le. cx) then
cf                   symbol(m) = tb
cf                else if (v .le. cy) then
cf                   symbol(m) = tc
cf                else if (v .le. cz) then
cf                   symbol(m) = td
cf                else
cf                   symbol(m) = te
cf                end if
cf             else if (k .eq. i) then
cf                symbol(m) = digit(nint((a(i,i)-dmin)/tcl))
cf             else if (k .gt. i) then
cf                v = abs(a(i,k))
cf                if (v .le. ca) then
cf                   symbol(m) = ta
cf                else if (v .le. cb) then
cf                   symbol(m) = tb
cf                else if (v .le. cc) then
cf                   symbol(m) = tc
cf                else if (v .le. cd) then
cf                   symbol(m) = td
cf                else
cf                   symbol(m) = te
cf                end if
cf             end if
cf          end do
cf          write (iout,50)  (symbol(k),k=1,nrow)
cf 50       format (1x,130a1)
cf       end do
cf       write (iout,60)  (dash,i=1,ndash)
cf 60    format (/,1x,130a1)
cf       if (maxj .lt. n) then
cf          write (iout,70)
cf 70       format ()
cf       end if
cf    end do
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dstmat  --  choose values for distance matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dstmat" selects a distance matrix containing values between
c     the previously smoothed upper and lower bounds; the distance
c     values are chosen as uniform random fractions, in a triangle
c     correlated fashion, or via various metrization algorithms
c
c
cf    subroutine dstmat (dmx)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    include 'keys.i'
cf    include 'restrn.i'
cf    integer maxpair
cf    parameter (maxpair=maxgeo*(maxgeo-1)/2)
cf    integer i,j,k,m,index,next
cf    integer npart,nmetrize,npair
cf    integer mik,mjk,nik,njk
cf    integer list(maxpair)
cf    real*8 random,normal
cf    real*8 fraction,fracdist
cf    real*8 corr,dmean,dwidth
cf    real*8 denom,swap,delta
cf    real*8 eps,gap,secs
cf    real*8 value(maxpair)
cf    real*8 dmx(maxgeo,maxgeo)
cf    character*8 method
cf    character*20 keyword
cf    character*80 record,string
cf    logical initial
cf    save initial,method,npart,dmean
cf    data initial  / .true. /
c
c
c     initialize the method for distance element selection
c
cf    if (initial) then
cf       initial = .false.
cf       method = 'pairwise'
cf       npart = 0
cf       dmean = 0.0d0
c
c     search each line of the keyword file for options
c
cf       do i = 1, nkey
cf          next = 1
cf          record = keyline(i)
cf          call gettext (record,keyword,next)
cf          call upcase (keyword)
cf          if (keyword(1:15) .eq. 'TRIAL-DISTANCE ') then
cf             call gettext (record,method,next)
cf             call lowcase (method)
cf             call getnumb (record,npart,next)
cf          else if (keyword(1:19) .eq. 'TRIAL-DISTRIBUTION ') then
cf             string = record(next:80)
cf             read (string,*,err=10,end=10)  dmean
cf 10          continue
cf          end if
cf       end do
c
c     set the number of atoms to use in partial metrization
c
cf       if (method .eq. 'havel') then
cf          if (npart.le.0 .or. npart.ge.n-1)  npart = n
cf       else if (method .eq. 'partial') then
cf          if (npart.le.0 .or. npart.ge.n-1)  npart = 4
cf       else if (method .eq. 'pairwise') then
cf          if (npart.le.0 .or. npart.ge.100)  npart = 100
cf       end if
c
c     set the initial distribution for selection of trial distances
c
cf       if (method.eq.'random' .or. method.eq.'partial'
cf   &             .or. method.eq.'pairwise') then
cf          compact = 0.0d0
c           if (dmean .eq. 0.0d0)  dmean = 2.35d0 / log(pathmax)
cf          if (dmean .eq. 0.0d0)  dmean = 1.65d0 / (pathmax)**0.25d0
c           if (dmean .eq. 0.0d0)  dmean = 1.30d0 / (pathmax)**0.20d0
cf       end if
cf    end if
c
c     write out the final choice for distance matrix generation
c
cf    if (verbose) then
cf       call setime
cf       if (method .eq. 'classic') then
cf          write (iout,20)
cf 20       format (/,' Distance Matrix via Uniform Random',
cf   &                 ' Fractions without Metrization :')
cf       else if (method .eq. 'random') then
cf          write (iout,30)
cf 30       format (/,' Distance Matrix Generated via Normal',
cf   &                 ' Fractions without Metrization :')
cf       else if (method .eq. 'tricor') then
cf          write (iout,40)
cf 40       format (/,' Distance Matrix Generated via Triangle',
cf   &                 ' Correlated Fractions :')
cf       else if (method.eq.'havel' .and. npart.lt.n) then
cf          write (iout,50)  npart
cf 50       format (/,' Distance Matrix Generated via',i4,'-Atom',
cf   &                 ' Partial Metrization :')
cf       else if (method .eq. 'havel') then
cf          write (iout,60)
cf 60       format (/,' Distance Matrix Generated via Randomized',
cf   &                 ' Atom-Based Metrization :')
cf       else if (method .eq. 'partial') then
cf          write (iout,70)  npart
cf 70       format (/,' Distance Matrix Generated via',i4,'-Atom',
cf   &                 ' Partial Metrization :')
cf       else if (method.eq.'pairwise' .and. npart.lt.100) then
cf          write (iout,80)  npart
cf 80       format (/,' Distance Matrix Generated via',i3,'%',
cf   &                 ' Random Pairwise Metrization :')
cf       else
cf          write (iout,90)
cf 90       format (/,' Distance Matrix Generated via Randomized',
cf   &                 ' Pairwise Metrization :')
cf       end if
cf    end if
c
c     adjust the distribution for selection of trial distances
c
cf    if (method.eq.'random' .or. method.eq.'partial'
cf   &          .or. method.eq.'pairwise') then
cf       dmean = dmean - 0.02d0*sign(sqrt(abs(compact)),compact)
cf       dwidth = 0.15d0
cf       if (verbose) then
cf          write (iout,100)  dmean,dwidth
cf100       format (/,' Trial Distance Normal Distribution :',
cf   &                 10x,f5.2,' +/-',f5.2)
cf       end if
cf    else
cf       write (iout,110)
cf110    format (/,' Trial Distances selected at Random from',
cf   &              ' Uniform Distribution')
cf    end if
c
c     uniform or normally distributed distances without metrization
c
cf    if (method.eq.'classic' .or. method.eq.'random') then
cf       do i = 1, n
cf          dmx(i,i) = 0.0d0
cf       end do
cf       do i = 1, n-1
cf          do j = i+1, n
cf             if (method .eq. 'classic') then
cf                fraction = random ()
cf             else
cf                fraction = fracdist (dmean,dwidth)
cf             end if
cf             delta = bnd(i,j) - bnd(j,i)
cf             dmx(j,i) = bnd(j,i) + delta*fraction
cf             dmx(i,j) = dmx(j,i)
cf          end do
cf       end do
c
c     Crippen's triangle correlated distance selection
c
cf    else if (method .eq. 'tricor') then
cf       do i = 1, n
cf          dmx(i,i) = 0.0d0
cf       end do
cf       do i = 1, n-1
cf          do j = i+1, n
cf             dmx(j,i) = random ()
cf             dmx(i,j) = dmx(j,i)
cf          end do
cf       end do
cf       do i = 1, n-1
cf          do j = i+1, n
cf             denom = 0.0d0
cf             dmx(i,j) = 0.0d0
cf             do k = 1, n
cf                if (k .ne. i) then
cf                   mik = max(i,k)
cf                   mjk = max(j,k)
cf                   nik = min(i,k)
cf                   njk = min(j,k)
cf                   if (k .eq. j) then
cf                      dmx(i,j) = dmx(i,j) + dmx(j,i)
cf                      denom = denom + 1.0d0
cf                   else if (bnd(njk,mjk) .le. 0.2d0*bnd(nik,mik)) then
cf                      if (i .gt. k)  corr = 0.9d0 * dmx(i,k)
cf                      if (k .gt. i)  corr = 0.9d0 * dmx(k,i)
cf                      dmx(i,j) = dmx(i,j) + corr
cf                      denom = denom + 0.9d0
cf                   else if (bnd(nik,mik) .le. 0.2d0*bnd(njk,mjk)) then
cf                      if (j .gt. k)  corr = 0.9d0 * dmx(j,k)
cf                      if (k .gt. j)  corr = 0.9d0 * dmx(k,j)
cf                      dmx(i,j) = dmx(i,j) + corr
cf                      denom = denom + 0.9d0
cf                   else if (bnd(mik,nik) .ge. 0.9d0*bnd(njk,mjk)) then
cf                      if (j .gt. k)  corr = 0.5d0 * (1.0d0-dmx(j,k))
cf                      if (k .gt. j)  corr = 0.5d0 * (1.0d0-dmx(k,j))
cf                      dmx(i,j) = dmx(i,j) + corr
cf                      denom = denom + 0.5d0
cf                   else if (bnd(mjk,njk) .ge. 0.9d0*bnd(nik,mik)) then
cf                      if (i .gt. k)  corr = 0.5d0 * (1.0d0-dmx(i,k))
cf                      if (k .gt. i)  corr = 0.5d0 * (1.0d0-dmx(k,i))
cf                      dmx(i,j) = dmx(i,j) + corr
cf                      denom = denom + 0.5d0
cf                   end if
cf                end if
cf             end do
cf             dmx(i,j) = dmx(i,j) / denom
cf          end do
cf       end do
cf       do i = 1, n-1
cf          do j = i+1, n
cf             delta = bnd(i,j) - bnd(j,i)
cf             dmx(i,j) = bnd(j,i) + delta*dmx(i,j)
cf             dmx(j,i) = dmx(i,j)
cf          end do
cf       end do
c
c     Havel/XPLOR atom-based metrization over choice of distribution
c
cf    else if (method.eq.'havel' .or. method.eq.'partial') then
cf       do i = 1, n
cf          do j = 1, n
cf             dmx(j,i) = bnd(j,i)
cf          end do
cf       end do
cf       do i = 1, n
cf          value(i) = random ()
cf       end do
cf       call sort2 (n,value,list)
cf       gap = 0.0d0
cf       do i = 1, n-1
cf          k = list(i)
cf          do j = i+1, n
cf             m = list(j)
cf             if (method .eq. 'havel') then
cf                fraction = random ()
cf             else
cf                fraction = fracdist (dmean,dwidth)
cf             end if
cf             delta = abs(bnd(k,m) - bnd(m,k))
cf             if (k .lt. m) then
cf                bnd(k,m) = bnd(m,k) + delta*fraction
cf                bnd(m,k) = bnd(k,m)
cf             else
cf                bnd(k,m) = bnd(k,m) + delta*fraction
cf                bnd(m,k) = bnd(k,m)
cf             end if
cf             if (i .le. npart)  call trifix (k,m)
cf             if (i .gt. npart)  gap = gap + delta
cf          end do
cf       end do
cf       do i = 1, n
cf          do j = 1, n
cf             swap = dmx(j,i)
cf             dmx(j,i) = bnd(j,i)
cf             bnd(j,i) = swap
cf          end do
cf       end do
cf       if (verbose .and. npart.lt.n-1) then
cf          write (iout,120)  gap/dble((n-npart)*(n-npart-1)/2)
cf120       format (/,' Average Bound Gap after Partial Metrization :',
cf   &                 3x,f12.4)
cf       end if
c
c     use partial randomized pairwise distance-based metrization
c
cf    else if (method.eq.'pairwise' .and. npart.le.10) then
cf       npair = n*(n-1) / 2
cf       nmetrize = nint(0.01d0*dble(npart)*dble(npair))
cf       do i = 1, n
cf          do j = 1, n
cf             dmx(j,i) = bnd(j,i)
cf          end do
cf       end do
cf       do i = 1, nmetrize
cf130       continue
cf          k = int(dble(n)*random()) + 1
cf          m = int(dble(n)*random()) + 1
cf          if (bnd(k,m) .eq. bnd(m,k))  goto 130
cf          if (k .gt. m) then
cf             swap = k
cf             k = m
cf             m = swap
cf          end if
cf          fraction = fracdist (dmean,dwidth)
cf          delta = bnd(k,m) - bnd(m,k)
cf          bnd(k,m) = bnd(m,k) + delta*fraction
cf          bnd(m,k) = bnd(k,m)
cf          call trifix (k,m)
cf       end do
cf       gap = 0.0d0
cf       do i = 1, n-1
cf          do j = i, n
cf             delta = bnd(i,j) - bnd(j,i)
cf             if (delta .ne. 0.0d0) then
cf                gap = gap + delta
cf                fraction = dmean + dwidth*normal()
cf                fraction = max(0.0d0,min(1.0d0,fraction))
cf                bnd(i,j) = bnd(j,i) + delta*fraction
cf                bnd(j,i) = bnd(i,j)
cf             end if
cf          end do
cf       end do
cf       do i = 1, n
cf          do j = 1, n
cf             swap = dmx(j,i)
cf             dmx(j,i) = bnd(j,i)
cf             bnd(j,i) = swap
cf          end do
cf       end do
cf       if (verbose .and. nmetrize.lt.npair) then
cf          write (iout,140)  gap/dble(npair-nmetrize)
cf140       format (/,' Average Bound Gap after Partial Metrization :',
cf   &                 3x,f12.4)
cf       end if
c
c     use randomized pairwise distance-based metrization
c
cf    else if (method .eq. 'pairwise') then
cf       npair = n*(n-1) / 2
cf       nmetrize = nint(0.01d0*dble(npart)*dble(npair))
cf       do i = 1, n
cf          do j = 1, n
cf             dmx(j,i) = bnd(j,i)
cf          end do
cf       end do
cf       do i = 1, npair
cf          value(i) = random ()
cf       end do
cf       call sort2 (npair,value,list)
cf       eps = 1.0d-10
cf       gap = 0.0d0
cf       do i = 1, npair
cf          index = list(i)
cf          k = int(0.5d0 * (dble(2*n+1)
cf   &                 - sqrt(dble(4*n*(n-1)-8*index+9))) + eps)
cf          m = n*(1-k) + k*(k+1)/2 + index
cf          fraction = fracdist (dmean,dwidth)
cf          delta = bnd(k,m) - bnd(m,k)
cf          bnd(k,m) = bnd(m,k) + delta*fraction
cf          bnd(m,k) = bnd(k,m)
cf          if (i .le. nmetrize)  call trifix (k,m)
cf          if (i .gt. nmetrize)  gap = gap + delta
cf       end do
cf       do i = 1, n
cf          do j = 1, n
cf             swap = dmx(j,i)
cf             dmx(j,i) = bnd(j,i)
cf             bnd(j,i) = swap
cf          end do
cf       end do
cf       if (verbose .and. nmetrize.lt.npair) then
cf          write (iout,150)  gap/dble(npair-nmetrize)
cf150       format (/,' Average Bound Gap after Partial Metrization :',
cf   &                 3x,f12.4)
cf       end if
cf    end if
c
c     get the time required for distance matrix generation
c
cf    if (verbose) then
cf       call getime (secs)
cf       write (iout,160)  secs
cf160    format (/,' Time Required for Distance Matrix :',5x,f12.2,
cf   &              ' seconds')
cf    end if
cf    return
cf    end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function fracdist  --  random fraction distance value  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fracdist" computes a random fractional distance in the
c     range [0,1] based on a skewed Gaussian distribution with
c     the input mean and standard deviation
c
c
cf    function fracdist (dmean,dwidth)
cf    real*8 dmean,dwidth
cf    real*8 fracdist,normal
c
c
c     compute a random fractional distance between 0 and 1
c
cf    fracdist = normal ()
cf    if (fracdist .ge. 0.0d0) then
cf       if (1.0d0-dmean .ge. 3.0d0*dwidth) then
cf          fracdist = dmean + fracdist*dwidth
cf       else
cf          fracdist = dmean + fracdist*(1.0d0-dmean)/3.0d0
cf       end if
cf    else
cf       if (dmean .ge. 3.0d0*dwidth) then
cf          fracdist = dmean + fracdist*dwidth
cf       else
cf          fracdist = dmean + fracdist*dmean/3.0d0
cf       end if
cf    end if
cf    fracdist = max(0.0d0,min(1.0d0,fracdist))
cf    return
cf    end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine metric  --  computation of the metric matrix  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "metric" takes as input the trial distance matrix and computes
c     the metric matrix of all possible dot products between the atomic
c     vectors and the center of mass using the law of cosines and the
c     following formula for the distances to the center of mass:
c
c        dcm(i)**2 = (1/n) * sum(j=1,n)(dist(i,j)**2)
c                          - (1/n**2) * sum(j<k)(dist(j,k)**2)
c
c     upon output, the metric matrix is stored in the lower triangle
c     plus diagonal of the input trial distance matrix, the upper
c     triangle of the input matrix is unchanged
c
c     literature reference:
c
c     G. M. Crippen and T. F. Havel, "Stable Calculation of Coordinates
c     from Distance Information", Acta Cryst., A34, 282-284 (1978)
c
c
cf    subroutine metric (gmx,nneg)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    integer i,j,nneg
cf    real*8 total,sum,rg
cf    real*8 dsq(maxgeo),dcm(maxgeo)
cf    real*8 gmx(maxgeo,maxgeo)
c
c
c     square and sum trial distances to get radius of gyration
c
cf    total = 0.0d0
cf    do i = 1, n
cf       do j = i, n
cf          gmx(j,i) = gmx(j,i)**2
cf          total = total + gmx(j,i)
cf       end do
cf    end do
cf    total = total / dble(n**2)
cf    rg = sqrt(total)
c
c     sum squared distances from each atom; the center
c     of mass is derived using the formula shown above
c
cf    nneg = 0
cf    do i = 1, n
cf       sum = 0.0d0
cf       do j = 1, i-1
cf          sum = sum + gmx(i,j)
cf       end do
cf       do j = i, n
cf          sum = sum + gmx(j,i)
cf       end do
cf       dsq(i) = sum/dble(n) - total
cf       dcm(i) = sqrt(abs(dsq(i)))
cf       if (dsq(i) .lt. 0.0d0) then
cf          nneg = nneg + 1
cf          dcm(i) = -dcm(i)
cf       end if
cf    end do
c
c     calculate the metric matrix using the law of cosines, and
c     place into the lower triangle of the input distance matrix
c
cf    do i = 1, n
cf       do j = i, n
cf          gmx(j,i) = 0.5d0 * (dsq(i)+dsq(j)-gmx(j,i))
cf       end do
cf    end do
c
c     print out the results of metric matrix computation
c
cf    if (verbose) then
cf       write (iout,10)  rg
cf 10    format (/,' Radius of Gyration of the System :',10x,f16.4)
cf    end if
cf    if (debug .or. (verbose .and. n.lt.130)) then
cf       write (iout,20)
cf 20    format (/,' Atomic Distances to the Center of Mass :',/)
cf       write (iout,30)  (dcm(i),i=1,n)
cf 30    format (6f13.4)
cf    end if
cf    return
cf    end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eigen  --  largest eigenvalues of metric metrix  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eigen" uses the power method to compute the largest eigenvalues
c     and eigenvectors of the metric matrix, "valid" is set true if the
c     first three eigenvalues are positive
c
c
cf    subroutine eigen (evl,evc,gmx,valid)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    integer maxeigen
cf    parameter (maxeigen=5)
cf    integer i,j,neigen
cf    real*8 evl(maxeigen)
cf    real*8 evc(maxgeo,maxeigen)
cf    real*8 gmx(maxgeo,maxgeo)
cf    real*8 secs,work(maxgeo)
cf    logical valid
c
c
c     initialize number of eigenvalues and convergence criteria
c
cf    if (verbose)  call setime
cf    neigen = 3
c
c     compute the largest few eigenvalues via the power method
c
cf    call power (n,maxgeo,neigen,gmx,evl,evc,work)
c
c     check to see if the first three eigenvalues are positive
c
cf    valid = .true.
cf    do i = 1, 3
cf       if (evl(i) .lt. 0.0d0)  valid = .false.
cf    end do
c
c     print out the eigenvalues and their eigenvectors
c
cf    if (verbose) then
cf       write (iout,10)
cf 10    format (/,' Eigenvalues from Metric Matrix :',/)
cf       write (iout,20)  (evl(i),i=1,neigen)
cf 20    format (5f15.4)
cf    end if
cf    if (debug) then
cf       write (iout,30)
cf 30    format (/,' Eigenvectors from Metric Matrix :',/)
cf       do i = 1, n
cf          write (iout,40)  (evc(i,j),j=1,neigen)
cf 40       format (5f15.4)
cf       end do
cf    end if
c
c     get the time required for partial matrix diagonalization
c
cf    if (verbose) then
cf       call getime (secs)
cf       write (iout,50)  secs
cf 50    format (/,' Time Required for Eigenvalues :',9x,f12.2,
cf   &              ' seconds')
cf    end if
cf    return
cf    end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine coords  --  converts eigenvalues to coordinates  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "coords" converts the three principal eigenvalues/vectors from
c     the metric matrix into atomic coordinates, and calls a routine
c     to compute the rms deviation from the bounds
c
c
cf    subroutine coords (evl,evc)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    integer maxeigen
cf    parameter (maxeigen=5)
cf    integer i,j,neigen
cf    real*8 rg,evl(maxeigen)
cf    real*8 evc(maxgeo,maxeigen)
cf    character*60 title
c
c
c     compute coordinates from the largest eigenvalues and vectors
c
cf    neigen = 3
cf    do j = 1, neigen
cf       evl(j) = sqrt(abs(evl(j)))
cf    end do
cf    do j = 1, neigen
cf       do i = 1, n
cf          evc(i,j) = evl(j) * evc(i,j)
cf       end do
cf    end do
c
c     transfer the final coordinates back to atomic vectors
c
cf    do i = 1, n
cf       x(i) = evc(i,1)
cf       y(i) = evc(i,2)
cf       z(i) = evc(i,3)
cf    end do
c
c     find the rms bounds deviations and radius of gyration
c
cf    if (verbose) then
cf       title = 'after Projection to 3-Dimensions :'
cf       call rmserror (title)
cf       call gyrate (rg)
cf       write (iout,10)  rg
cf 10    format (/,' Radius of Gyration of the System :',10x,f16.4)
cf    end if
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine chksize  --  estimate compaction or expansion  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "chksize" computes a measure of overall global structural
c     expansion or compaction from the number of excess upper
c     or lower bounds matrix violations
c
c
cf    subroutine chksize
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'couple.i'
cf    include 'disgeo.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    integer i,k,npair,nskip
cf    integer nlarge,nsmall
cf    integer skip(maxatm)
cf    real*8 xi,yi,zi,distmin
cf    real*8 dstsq,bupsq,blosq
c
c
c     zero out the list of atoms locally connected to each atom
c
cf    nskip = 0
cf    do i = 1, n
cf       skip(i) = 0
cf    end do
c
c     initialize counters, total pair number, and cutoff distance
c
cf    nlarge = 0
cf    nsmall = 0
cf    npair = n*(n-1) / 2
cf    distmin = 0.1d0 * (pathmax-100.0d0)
c
c     count the number of excess upper or lower bound violations
c
cf    do i = 1, n-1
cf       xi = x(i)
cf       yi = y(i)
cf       zi = z(i)
cf       do k = 1, n12(i)
cf          skip(i12(k,i)) = i
cf       end do
cf       do k = 1, n13(i)
cf          skip(i13(k,i)) = i
cf       end do
cf       do k = 1, n14(i)
cf          skip(i14(k,i)) = i
cf       end do
cf       do k = i+1, n
cf          if (skip(k) .eq. i) then
cf             nskip = nskip + 1
c           else if (bnd(i,k) .lt. distmin) then
c              nskip = nskip + 1
cf          else
cf             dstsq = (x(k)-xi)**2 + (y(k)-yi)**2 + (z(k)-zi)**2
cf             bupsq = bnd(i,k)**2
cf             blosq = bnd(k,i)**2
cf             if (dstsq .gt. bupsq) then
cf                nlarge = nlarge + 1
cf             else if (blosq .gt. dstsq) then
cf                nsmall = nsmall + 1
cf             end if
cf          end if
cf       end do
cf    end do
c
c     set the value for the overall index of compaction
c
cf    compact = 100.0d0 * dble(nlarge-nsmall)/dble(npair-nskip)
cf    if (verbose) then
cf       write (iout,10)  compact
cf 10    format (/,' Index of Structure Expansion/Compaction :',
cf   &              7x,f12.4)
cf    end if
cf    return
cf    end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine majorize  --  Guttman transform majorization  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "majorize" refines the projected coordinates by attempting to
c     minimize the least square residual between the trial distance
c     matrix and the distances computed from the coordinates
c
c     literature reference:
c
c     T. F. Havel, "An Evaluation of Computational Strategies for
c     Use in the Determination of Protein Structure from Distance
c     Constraints obtained by Nuclear Magnetic Resonance", Progress
c     in Biophysics and Molecular Biology, 56, 43-78 (1991)
c
c
cf    subroutine majorize (dmx)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    integer i,k,iter,niter,period
cf    real*8 pairs,dn1,dn2,rg,secs
cf    real*8 target,dist,error
cf    real*8 rmserr,average
cf    real*8 xi,yi,zi,b(maxatm)
cf    real*8 dmx(maxgeo,maxgeo)
cf    real*8 xx(maxatm),yy(maxatm),zz(maxatm)
cf    character*60 title
c
c
c     set number of iterations and some other needed values
c
cf    if (verbose)  call setime
cf    niter = 20
cf    period = 5
cf    pairs = dble(n*(n-1)/2)
cf    dn1 = dble(n-1)
cf    dn2 = dble(n*n)
c
c     find the average and rms error from trial distances
c
cf    iter = 0
cf    rmserr = 0.0d0
cf    average = 0.0d0
cf    do i = 1, n-1
cf       xi = x(i)
cf       yi = y(i)
cf       zi = z(i)
cf       do k = i+1, n
cf          target = dmx(k,i)
cf          dist = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
cf          error = dist - target
cf          rmserr = rmserr + error**2
cf          average = average + error/target
cf       end do
cf    end do
cf    rmserr = sqrt(rmserr/pairs)
cf    average = 100.0d0 * average / pairs
c
c     write a header with the initial error values
c
cf    if (verbose) then
cf       write (iout,10)
cf 10    format (/,' Majorization to Trial Distances using',
cf   &              ' Constant Weights :',
cf   &           //,4x,'Iteration',6x,'RMS Error',5x,'Ave % Error',/)
cf       write (iout,20)  iter,rmserr,average
cf 20    format (5x,i5,2f16.4)
cf    end if
c
c     initialize the transformed coordinates for each atom
c
cf    do iter = 1, niter
cf       do i = 1, n
cf          xi = x(i)
cf          yi = y(i)
cf          zi = z(i)
cf          b(i) = 0.0d0
cf          xx(i) = 0.0d0
cf          yy(i) = 0.0d0
cf          zz(i) = 0.0d0
c
c     form a single row of the B matrix assuming unity weights
c
cf          do k = 1, i-1
cf             dist = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
cf             b(k) = -dmx(k,i) / dist
cf             b(i) = b(i) - b(k)
cf          end do
cf          do k = i+1, n
cf             dist = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
cf             b(k) = -dmx(k,i) / dist
cf             b(i) = b(i) - b(k)
cf          end do
c
c     multiply the row of the B matrix by the atomic coordinates
c
cf          do k = 1, n
cf             xx(i) = xx(i) + b(k)*x(k)
cf             yy(i) = yy(i) + b(k)*y(k)
cf             zz(i) = zz(i) + b(k)*z(k)
cf          end do
cf       end do
c
c     move the intermediate values into the coordinate arrays
c
cf       do i = 1, n
cf          x(i) = xx(i)
cf          y(i) = yy(i)
cf          z(i) = zz(i)
cf       end do
c
c     multiply the inverse weight matrix S+ by the coordinates
c
cf       do i = 1, n
cf          xx(i) = (dn1/dn2) * x(i)
cf          yy(i) = (dn1/dn2) * y(i)
cf          zz(i) = (dn1/dn2) * z(i)
cf          do k = 1, i-1
cf             xx(i) = xx(i) - x(k)/dn2
cf             yy(i) = yy(i) - y(k)/dn2
cf             zz(i) = zz(i) - z(k)/dn2
cf          end do
cf          do k = i+1, n
cf             xx(i) = xx(i) - x(k)/dn2
cf             yy(i) = yy(i) - y(k)/dn2
cf             zz(i) = zz(i) - z(k)/dn2
cf          end do
cf       end do
c
c     copy the new coordinates into their permanent arrays
c
cf       do i = 1, n
cf          x(i) = xx(i)
cf          y(i) = yy(i)
cf          z(i) = zz(i)
cf       end do
c
c     find the average and rms error from trial distances
c
cf       rmserr = 0.0d0
cf       average = 0.0d0
cf       do i = 1, n-1
cf          xi = x(i)
cf          yi = y(i)
cf          zi = z(i)
cf          do k = i+1, n
cf             target = dmx(k,i)
cf             dist = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
cf             error = dist - target
cf             rmserr = rmserr + error**2
cf             average = average + error/target
cf          end do
cf       end do
cf       rmserr = sqrt(rmserr/pairs)
cf       average = 100.0d0 * average / pairs
cf       if (verbose .and. mod(iter,period).eq.0) then
cf          write (iout,30)  iter,rmserr,average
cf 30       format (5x,i5,2f16.4)
cf       end if
cf    end do
c
c     find the rms bounds deviations and radius of gyration
c
cf    if (verbose) then
cf       title = 'after Majorization :'
cf       call rmserror (title)
cf       call gyrate (rg)
cf       write (iout,40)  rg
cf 40    format (/,' Radius of Gyration of the System :',10x,f16.4)
c
c     get the time required for the majorization procedure
c
cf       call getime (secs)
cf       write (iout,50)  secs
cf 50    format (/,' Time Required for Majorization :',8x,f12.2,
cf   &              ' seconds')
cf    end if
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine refine  --  minimization of initial embedding  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "refine" performs minimization of the atomic coordinates
c     of an initial crude embedded distance geometry structure versus
c     the bound, chirality, planarity and torsional error functions
c
c
ckk -minimize-      subroutine refine (mode,fctval,grdmin)
cf    subroutine refine (method,mode,fctval,grdmin)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    include 'minima.i'
cf    include 'output.i'
cf    integer i,nvar
cf    real*8 initer,miderr,toterr
cf    real*8 fctval,grdmin,xx(maxvar)
cf    character*7 mode
cf    external initer,miderr,toterr,writeout
ckk -minimize- begin
cf    character*3 method
ckk -minimize- end
c
c
c     translate the atomic coordinates to optimization variables
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       xx(nvar) = x(i)
cf       nvar = nvar + 1
cf       xx(nvar) = y(i)
cf       nvar = nvar + 1
cf       xx(nvar) = z(i)
cf    end do
c
c     set values of parameters needed for optimization
c
cf    coordtype = 'none'
cf    savecycle = .true.
c     grdmin = 0.01d0
cf    maxiter = 2 * nvar
cf    iwrite = 0
c     iprint = 0
c     if (verbose)  iprint = 10
c
c     minimize initially only on the local geometry and torsions,
c     then on local geometry and chirality, torsions, and finally
c     minimize on all distance bounds, chirality and torsions
c
cf    if (mode .eq. 'initial') then
ckk -minimize-         call lmqn (nvar,xx,fctval,grdmin,initer,writeout)
cf       call lmqn (method,nvar,xx,fctval,grdmin,initer,writeout)
cf    else if (mode .eq. 'middle') then
ckk -minimize-         call lmqn (nvar,xx,fctval,grdmin,miderr,writeout)
cf       call lmqn (method,nvar,xx,fctval,grdmin,miderr,writeout)
cf    else if (mode .eq. 'final') then
ckk -minimize-         call lmqn (nvar,xx,fctval,grdmin,toterr,writeout)
cf       call lmqn (method,nvar,xx,fctval,grdmin,toterr,writeout)
cf    end if
c
c     translate optimization variables back to coordinates
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       x(i) = xx(nvar)
cf       nvar = nvar + 1
cf       y(i) = xx(nvar)
cf       nvar = nvar + 1
cf       z(i) = xx(nvar)
cf    end do
cf    return
cf    end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine explore  --  simulated annealing refinement  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "explore" uses simulated annealing on an initial crude
c     embedded distance geoemtry structure to refine versus the
c     bound, chirality, planarity and torsional error functions
c
c
cf    subroutine explore (mode,nstep,dt,mass,temp_start,temp_stop,v,a)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'inform.i'
cf    include 'iounit.i'
cf    include 'math.i'
cf    include 'units.i'
cf    integer i,istep,nstep,nvar,period
cf    real*8 error,total,prior,change
cf    real*8 dt,dt2,dt_2,dt2_2
cf    real*8 xbig,xrms,mass,kinetic
cf    real*8 temp_start,temp_stop
cf    real*8 tau_start,tau_stop
cf    real*8 ratio,sigmoid,scale
cf    real*8 target,temp,tautemp
cf    real*8 initer,miderr,toterr
cf    real*8 xx(maxvar),xmove(maxvar)
cf    real*8 g(maxvar),v(maxvar),a(maxvar)
cf    character*7 mode
c
c
c     set values of the basic simulated annealing parameters
c
c     nstep = 5000
c     dt = 1.0d-13
c     temp_start = 200.0d0
c     temp_stop = 0.0d0
c     mass = 1000.0d0
c
c     translate the atomic coordinates to annealing variables
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       xx(nvar) = x(i)
cf       nvar = nvar + 1
cf       xx(nvar) = y(i)
cf       nvar = nvar + 1
cf       xx(nvar) = z(i)
cf    end do
c
c     initialize the velocities, accelerations and other parameters
c
cf    dt2 = dt * dt
cf    dt_2 = dt / 2.0d0
cf    dt2_2 = dt2 / 2.0d0
cf    period = 100
cf    tau_start = 100.0d0 * dt
cf    tau_stop = 10.0d0 * dt
cf    tautemp = tau_start
c
c     print a header for the simulated annealing protocol
c
cf    write (iout,10)
cf 10 format (/,' Molecular Dynamics Simulated Annealing Refinement :')
cf    write (iout,20)  nstep,1.0d12*dt,log(mass)/logten,
cf   &                 temp_start,temp_stop
cf 20 format (/,' Steps:',i6,3x,'Time/Step:',f6.3,' ps',3x,
cf   &           'LogMass:',f5.2,3x,'Temp:',f6.1,' to',f6.1)
c
c     get the total error and temperature at start of dynamics
c
cf    if (mode .eq. 'initial') then
cf       error = initer (xx,g)
cf    else if (mode .eq. 'middle') then
cf       error = miderr (xx,g)
cf    else if (mode .eq. 'final') then
cf       error = toterr (xx,g)
cf    end if
cf    kinetic = 0.0d0
cf    do i = 1, nvar
cf       kinetic = kinetic + mass*v(i)**2
cf    end do
cf    kinetic = 0.5d0 * kinetic / convert
cf    temp = 2.0d0 * kinetic / (dble(nvar) * gasconst)
cf    total = error + kinetic
cf    prior = total
cf    if (verbose) then
cf       write (iout,30)
cf 30    format (/,' MD Step    E Total   E Potential   E Kinetic',
cf   &              '     Temp    MaxMove   RMS Move',/)
cf       write (iout,40)  0,total,error,kinetic,temp
cf 40    format (i6,2f13.4,f12.4,f11.2)
cf    end if
c
c     find new positions and half-step velocities via Verlet
c
cf    do istep = 1, nstep
cf       xbig = 0.0d0
cf       xrms = 0.0d0
cf       do i = 1, nvar
cf          xmove(i) = v(i)*dt + a(i)*dt2_2
cf          xx(i) = xx(i) + xmove(i)
cf          v(i) = v(i) + a(i)*dt_2
cf          if (abs(xmove(i)) .gt. xbig)  xbig = abs(xmove(i))
cf          xrms = xrms + xmove(i)**2
cf       end do
cf       xrms = sqrt(xrms/dble(nvar))
c
c     get the error function value and gradient
c
cf       if (mode .eq. 'initial') then
cf          error = initer (xx,g)
cf       else if (mode .eq. 'middle') then
cf          error = miderr (xx,g)
cf       else if (mode .eq. 'final') then
cf          error = toterr (xx,g)
cf       end if
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Verlet recursion
c
cf       do i = 1, nvar
cf          a(i) = -convert * g(i) / mass
cf          v(i) = v(i) + a(i)*dt_2
cf       end do
c
c     find the total kinetic energy and system temperature
c
cf       kinetic = 0.0d0
cf       do i = 1, nvar
cf          kinetic = kinetic + mass*v(i)**2
cf       end do
cf       kinetic = 0.5d0 * kinetic / convert
cf       temp = 2.0d0 * kinetic / (dble(nvar) * gasconst)
cf       if (temp .eq. 0.0d0)  temp = 0.1d0
c
c     set target temperature and coupling via a sigmoidal cooling
c
cf       ratio = dble(istep) / dble(nstep)
cf       ratio = sigmoid (3.5d0,ratio)
cf       target = temp_start*(1.0d0-ratio) + temp_stop*ratio
cf       tautemp = tau_start*(1.0d0-ratio) + tau_stop*ratio
c
c     couple to external temperature bath via velocity scaling
c
cf       scale = sqrt(1.0d0 + (dt/tautemp)*(target/temp-1.0d0))
cf       do i = 1, nvar
cf          v(i) = scale * v(i)
cf       end do
c
c     write results for the current annealing step
c
cf       total = error + kinetic
cf       if (verbose .and. mod(istep,period).eq.0) then
cf          write (iout,50)  istep,total,error,kinetic,temp,xbig,xrms
cf 50       format (i6,2f13.4,f12.4,f11.2,2f10.4)
cf       end if
c
c     check the energy change for instability in the dynamics
c
cf       change = total - prior
cf       if (change .gt. dble(n)) then
cf          do i = 1, nvar
cf             xx(i) = xx(i) - xmove(i)
cf          end do
cf          if (verbose .and. mod(istep,period).ne.0) then
cf             write (iout,60)  istep,total,error,kinetic,temp,xbig,xrms
cf 60          format (i6,2f13.4,f12.4,f11.2,2f10.4)
cf          end if
cf          write (iout,70)
cf 70       format (/,' EXPLORE  --  Simulated Annealing Unstable;',
cf   &                 ' Switching to Minimization')
cf          goto 80
cf       end if
cf    end do
c
c     translate annealing variables back to coordinates
c
cf 80 continue
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       x(i) = xx(nvar)
cf       nvar = nvar + 1
cf       y(i) = xx(nvar)
cf       nvar = nvar + 1
cf       z(i) = xx(nvar)
cf    end do
cf    return
cf    end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine rmserror  --  rms bound and restraint error  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "rmserror" computes the maximum absolute deviation and the
c     rms deviation from the distance bounds, and the number and
c     rms value of the distance restraint violations
c
c
cf    subroutine rmserror (title)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'iounit.i'
cf    include 'restrn.i'
cf    integer i,j,k,npair
cf    integer nhierr,nloerr
cf    integer ihi,jhi,ilo,jlo
cf    integer leng,trimtext
cf    real*8 rms,himax,lomax
cf    real*8 dist,hierr,loerr
cf    character*60 title
c
c
c     search all atom pairs for maximal bounds deviations
c
cf    npair = n*(n-1) / 2
cf    nloerr = 0
cf    nhierr = 0
cf    ilo = 0
cf    jlo = 0
cf    ihi = 0
cf    jhi = 0
cf    rms = 0.0d0
cf    lomax = 0.0d0
cf    himax = 0.0d0
cf    do i = 1, n-1
cf       do j = i+1, n
cf          dist = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
cf          dist = sqrt(dist)
cf          hierr = dist - bnd(i,j)
cf          if (hierr .gt. 0.0d0) then
cf             nhierr = nhierr + 1
cf             rms = rms + hierr**2
cf             if (hierr .gt. himax) then
cf                himax = hierr
cf                ihi = i
cf                jhi = j
cf             end if
cf          end if
cf          loerr = bnd(j,i) - dist
cf          if (loerr .gt. 0.0d0) then
cf             nloerr = nloerr + 1
cf             rms = rms + loerr**2
cf             if (loerr .gt. lomax) then
cf                lomax = loerr
cf                ilo = i
cf                jlo = j
cf             end if
cf          end if
cf       end do
cf    end do
cf    rms = sqrt(rms/dble(n*(n-1)/2))
c
c     print the maximal and rms bound deviations
c
cf    leng = trimtext(title)
cf    write (iout,10)  title(1:leng)
cf 10 format (/,' Fit to Bounds ',a)
cf    write (iout,20)  nhierr,npair,nloerr,npair,himax,
cf   &                 ihi,jhi,lomax,ilo,jlo,rms
cf 20 format (/,' Num Upper Bound Violations :',4x,i11,'  of ',i12,
cf   &        /,' Num Lower Bound Violations :',4x,i11,'  of ',i12,
cf   &        /,' Max Upper Bound Violation :',4x,f12.4,'  at ',2i6,
cf   &        /,' Max Lower Bound Violation :',4x,f12.4,'  at ',2i6,
cf   &        /,' RMS Deviation from Bounds :',4x,f12.4)
c
c     search the list of distance restraints for violations
c
cf    if (ndfix .gt. 0) then
cf       nloerr = 0
cf       nhierr = 0
cf       ilo = 0
cf       jlo = 0
cf       ihi = 0
cf       jhi = 0
cf       rms = 0.0d0
cf       himax = 0.0d0
cf       lomax = 0.0d0
cf       do k = 1, ndfix
cf          i = idfix(1,k)
cf          j = idfix(2,k)
cf          dist = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
cf          dist = sqrt(dist)
cf          if (dist .lt. dfix(1,k)) then
cf             nloerr = nloerr + 1
cf             loerr = dfix(1,k) - dist
cf             rms = rms + (dist-dfix(1,k))**2
cf             if (loerr .gt. lomax) then
cf                lomax = loerr
cf                ilo = i
cf                jlo = j
cf             end if
cf          else if (dist .gt. dfix(2,k)) then
cf             nhierr = nhierr + 1
cf             hierr = dist - dfix(2,k)
cf             rms = rms + hierr**2
cf             if (hierr .gt. himax) then
cf                himax = hierr
cf                ihi = i
cf                jhi = j
cf             end if
cf          end if
cf       end do
cf       rms = sqrt(rms/dble(ndfix))
c
c     print total number and rms value of restraint violations
c
cf       write (iout,30)  nhierr,ndfix,nloerr,ndfix,himax,
cf   &                    ihi,jhi,lomax,ilo,jlo,rms
cf 30    format (/,' Num Upper Restraint Violations :',i11,'  of ',i12,
cf   &           /,' Num Lower Restraint Violations :',i11,'  of ',i12,
cf   &           /,' Max Upper Restraint Violation :',f12.4,'  at ',2i6,
cf   &           /,' Max Lower Restraint Violation :',f12.4,'  at ',2i6,
cf   &           /,' RMS Restraint Dist Violation : ',f12.4)
cf    end if
cf    return
cf    end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine dmdump  --  final distance and error matrix  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "dmdump" puts the distance matrix of the final structure
c     into the upper half of a matrix, the distance of each atom
c     to the centroid on the diagonal, and the individual terms
c     of the bounds errors into the lower half of the matrix
c
c
cf    subroutine dmdump (dmd)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'iounit.i'
cf    integer i,j
cf    real*8 sum,dist,dist2,rgsq
cf    real*8 dmd(maxgeo,maxgeo)
cf    character*60 title
c
c
c     store the final distance matrix and bound violations
c
cf    do i = 1, n
cf       dmd(i,i) = 0.0d0
cf    end do
cf    sum = 0.0d0
cf    do i = 1, n-1
cf       do j = i+1, n
cf          dist2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
cf          sum = sum + dist2
cf          dmd(i,i) = dmd(i,i) + dist2
cf          dmd(j,j) = dmd(j,j) + dist2
cf          dist = sqrt(dist2)
cf          dmd(i,j) = dist
cf          if (dist .gt. bnd(i,j)) then
cf             dmd(j,i) = dist - bnd(i,j)
cf          else if (dist .lt. bnd(j,i)) then
cf             dmd(j,i) = bnd(j,i) - dist
cf          else
cf             dmd(j,i) = 0.0d0
cf          end if
cf       end do
cf    end do
c
c     put the distance to the centroid on the diagonal
c
cf    rgsq = sum / dble(n**2)
cf    do i = 1, n
cf       dmd(i,i) = sqrt(dmd(i,i)/dble(n) - rgsq)
cf    end do
c
c     write out the interatomic distance and error matrices
c
cf    title = 'Final Dist Matrix Above; DCM on Diag; Error Below :'
cf    call grafic (n,maxgeo,dmd,title)
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function initer  --  initial error function and gradient  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "initer" is the initial error function and derivatives for
c     a distance geometry embedding; it includes components from
c     the local geometry and torsional constraint errors
c
c
cf    function initer (xx,g)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    integer i,j,nvar
cf    real*8 initer,local,torsion
cf    real*8 locerr,torser
cf    real*8 xx(maxvar),g(maxvar)
cf    real*8 derivs(3,maxgeo)
cf    external locerr,torser
c
c
c     translate optimization parameters to atomic coordinates
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       x(i) = xx(nvar)
cf       nvar = nvar + 1
cf       y(i) = xx(nvar)
cf       nvar = nvar + 1
cf       z(i) = xx(nvar)
cf    end do
c
c     zero out the values of the atomic gradient components
c
cf    do i = 1, n
cf       do j = 1, 3
cf          derivs(j,i) = 0.0d0
cf       end do
cf    end do
c
c     compute the local goemetry and the torsional
c     components of the error function and its gradient
c
cf    local = locerr (derivs)
cf    torsion = torser (derivs)
cf    initer = local + torsion
c
c     store the atomic gradients as the optimization gradient
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       g(nvar) = derivs(1,i)
cf       nvar = nvar + 1
cf       g(nvar) = derivs(2,i)
cf       nvar = nvar + 1
cf       g(nvar) = derivs(3,i)
cf    end do
cf    return
cf    end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function miderr  --  second error function and gradient  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "miderr" is the secondary error function and derivatives
c     for a distance geometry embedding; it includes components
c     from the distance bounds, local geometry, chirality and
c     torsional constraint errors
c
c
cf    function miderr (xx,g)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    integer i,j,nvar
cf    real*8 miderr,bounds,local
cf    real*8 chiral,torsion
cf    real*8 bnderr,locerr
cf    real*8 chirer,torser
cf    real*8 xx(maxvar),g(maxvar)
cf    real*8 derivs(3,maxgeo)
c
c
c     translate optimization parameters to atomic coordinates
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       x(i) = xx(nvar)
cf       nvar = nvar + 1
cf       y(i) = xx(nvar)
cf       nvar = nvar + 1
cf       z(i) = xx(nvar)
cf    end do
c
c     zero out the values of the atomic gradient components
c
cf    do i = 1, n
cf       do j = 1, 3
cf          derivs(j,i) = 0.0d0
cf       end do
cf    end do
c
c     compute the local geometry, chirality and torsional
c     components of the error function and its gradient
c
cf    bounds = bnderr (derivs)
cf    local = locerr (derivs)
cf    chiral = chirer (derivs)
cf    torsion = torser (derivs)
cf    miderr = bounds + local + chiral + torsion
c
c     store the atomic gradients as the optimization gradient
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       g(nvar) = derivs(1,i)
cf       nvar = nvar + 1
cf       g(nvar) = derivs(2,i)
cf       nvar = nvar + 1
cf       g(nvar) = derivs(3,i)
cf    end do
cf    return
cf    end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function toterr  --  total error function and gradient  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "toterr" is the error function and derivatives for a distance
c     geometry embedding; it includes components from the distance
c     bounds, hard sphere contacts, local geometry, chirality and
c     torsional constraint errors
c
c
cf    function toterr (xx,g)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    integer i,j,nvar
cf    real*8 toterr,bounds,contact
cf    real*8 local,chiral,torsion
cf    real*8 bnderr,vdwerr,locerr
cf    real*8 chirer,torser
cf    real*8 xx(maxvar),g(maxvar)
cf    real*8 derivs(3,maxgeo)
cf    external bnderr,vdwerr,chirer,torser
c
c
c     translate optimization parameters to atomic coordinates
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       x(i) = xx(nvar)
cf       nvar = nvar + 1
cf       y(i) = xx(nvar)
cf       nvar = nvar + 1
cf       z(i) = xx(nvar)
cf    end do
c
c     zero out the values of the atomic gradient components
c
cf    do i = 1, n
cf       do j = 1, 3
cf          derivs(j,i) = 0.0d0
cf       end do
cf    end do
c
c     compute the distance bound, vdw, chirality and torsional
c     components of the error function and its gradient
c
cf    bounds = bnderr (derivs)
cf    contact = vdwerr (derivs)
cf    local = locerr (derivs)
cf    chiral = chirer (derivs)
cf    torsion = torser (derivs)
cf    toterr = bounds + contact + local + chiral + torsion
c
c     store the atomic gradients as the optimization gradient
c
cf    nvar = 0
cf    do i = 1, n
cf       nvar = nvar + 1
cf       g(nvar) = derivs(1,i)
cf       nvar = nvar + 1
cf       g(nvar) = derivs(2,i)
cf       nvar = nvar + 1
cf       g(nvar) = derivs(3,i)
cf    end do
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function bnderr  --  computes total distance bound error  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bnderr" is the distance bound error function and derivatives;
c     this version implements the original and Havel's normalized
c     lower bound penalty, the normalized version is preferred when
c     lower bounds are small (as with NMR NOE constraints), the
c     original penalty is needed if large lower bounds are present
c
c
cf    function bnderr (derivs)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'restrn.i'
cf    integer i,j,k
cf    real*8 bnderr,error,cutsq
cf    real*8 scale,chain,term
cf    real*8 gap,buffer,weight
cf    real*8 dx,dy,dz,gx,gy,gz
cf    real*8 dstsq,bupsq,blosq
cf    real*8 derivs(3,maxgeo)
c
c
c     zero out the distance bounds error function
c
cf    bnderr = 0.0d0
cf    scale = 10.0d0
cf    cutsq = 40.0d0
c
c     calculate the pairwise distances between atoms
c
cf    do k = 1, ndfix
cf       i = idfix(1,k)
cf       j = idfix(2,k)
cf       dx = x(i) - x(j)
cf       dy = y(i) - y(j)
cf       dz = z(i) - z(j)
c
c     calculate squared actual distance and bound distances;
c     use of a small "buffer" cleans up the final error count
c
cf       dstsq = dx*dx + dy*dy + dz*dz
cf       gap = dfix(2,k) - dfix(1,k)
cf       buffer = 0.05d0 * min(1.0d0,gap)
cf       blosq = (dfix(1,k) + buffer)**2
cf       bupsq = (dfix(2,k) - buffer)**2
c
c     error and derivatives for upper bound violation
c
cf       if (dstsq .gt. bupsq) then
cf          weight = scale * dfix(3,k)
cf          term = (dstsq-bupsq) / bupsq
cf          chain = 4.0d0 * weight * term / bupsq
cf          error = weight * term**2
cf          gx = dx * chain
cf          gy = dy * chain
cf          gz = dz * chain
cf          bnderr = bnderr + error
cf          derivs(1,i) = derivs(1,i) + gx
cf          derivs(2,i) = derivs(2,i) + gy
cf          derivs(3,i) = derivs(3,i) + gz
cf          derivs(1,j) = derivs(1,j) - gx
cf          derivs(2,j) = derivs(2,j) - gy
cf          derivs(3,j) = derivs(3,j) - gz
c
c     error and derivatives for lower bound violation
c
cf       else if (dstsq .lt. blosq) then
cf          weight = scale * dfix(3,k)
cf          if (blosq .gt. cutsq) then
cf             term = (blosq-dstsq) / dstsq
cf             chain = -4.0d0 * weight * term * (blosq/dstsq**2)
cf          else
cf             term = (blosq-dstsq) / (blosq+dstsq)
cf             chain = -8.0d0 * weight * term * (blosq/(blosq+dstsq)**2)
cf          end if
cf          error = weight * term**2
cf          gx = dx * chain
cf          gy = dy * chain
cf          gz = dz * chain
cf          bnderr = bnderr + error
cf          derivs(1,i) = derivs(1,i) + gx
cf          derivs(2,i) = derivs(2,i) + gy
cf          derivs(3,i) = derivs(3,i) + gz
cf          derivs(1,j) = derivs(1,j) - gx
cf          derivs(2,j) = derivs(2,j) - gy
cf          derivs(3,j) = derivs(3,j) - gz
cf       end if
cf    end do
cf    return
cf    end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function vdwerr  --  computes van der Waals bound error  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "vdwerr" is the hard sphere van der Waals bound error function
c     and derivatives that penalizes close nonbonded contacts,
c     pairwise neighbors are generated via the method of lights
c
c
cf    function vdwerr (derivs)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'couple.i'
cf    include 'disgeo.i'
cf    include 'light.i'
cf    include 'shunt.i'
cf    integer i,j,k,kgy,kgz
cf    integer skip(maxgeo)
cf    integer map(maxlight)
cf    real*8 vdwerr,error
cf    real*8 scale,chain,term
cf    real*8 xi,yi,zi
cf    real*8 dx,dy,dz,gx,gy,gz
cf    real*8 dstsq,blosq
cf    real*8 radi,radsq
cf    real*8 derivs(3,maxgeo)
cf    real*8 xsort(maxlight)
cf    real*8 ysort(maxlight)
cf    real*8 zsort(maxlight)
c
c
c     zero out the distance van der Waals error function
c
cf    vdwerr = 0.0d0
cf    scale = 1.0d0
c
c     set maximum value of vdw radii sum for a pair of atoms
c
cf    off = 3.6d0
c
c     transfer coordinates and zero out atoms to be skipped
c
cf    do i = 1, n
cf       xsort(i) = x(i)
cf       ysort(i) = y(i)
cf       zsort(i) = z(i)
cf       skip(i) = 0
cf    end do
c
c     use the method of lights to generate neighbors
c
cf    call lights (n,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
cf    do i = 1, n
cf       radi = vdwrad(i)
cf       xi = xsort(rgx(i))
cf       yi = ysort(rgy(i))
cf       zi = zsort(rgz(i))
cf       do j = 1, n12(i)
cf          skip(i12(j,i)) = i
cf       end do
cf       do j = 1, n13(i)
cf          skip(i13(j,i)) = i
cf       end do
cf       do j = 1, n14(i)
cf          skip(i14(j,i)) = i
cf       end do
cf       do j = kbx(i)+1, kex(i)
cf          k = locx(j)
cf          if (skip(k) .eq. i)  goto 10
cf          kgy = rgy(k)
cf          if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 10
cf          kgz = rgz(k)
cf          if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 10
c
c     calculate squared distances and bounds
c
cf          dx = xi - xsort(j)
cf          dy = yi - ysort(kgy)
cf          dz = zi - zsort(kgz)
cf          dstsq = dx*dx + dy*dy + dz*dz
cf          radsq = (radi + vdwrad(k))**2
cf          blosq = min(bnd(k,i),bnd(i,k),radsq)
c
c     error and derivatives for lower bound violation
c
cf          if (dstsq .lt. blosq) then
cf             term = (blosq-dstsq) / (blosq+dstsq)
cf             chain = -8.0d0 * scale * term * (blosq/(blosq+dstsq)**2)
cf             error = scale * term**2
cf             gx = dx * chain
cf             gy = dy * chain
cf             gz = dz * chain
cf             vdwerr = vdwerr + error
cf             derivs(1,i) = derivs(1,i) + gx
cf             derivs(2,i) = derivs(2,i) + gy
cf             derivs(3,i) = derivs(3,i) + gz
cf             derivs(1,k) = derivs(1,k) - gx
cf             derivs(2,k) = derivs(2,k) - gy
cf             derivs(3,k) = derivs(3,k) - gz
cf          end if
cf 10       continue
cf       end do
cf    end do
cf    return
cf    end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function locerr  --  computes local geometry error value  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "locerr" is the local geometry error function and derivatives
c     including the 1-2, 1-3 and 1-4 distance bound constraints
c
c
cf    function locerr (derivs)
cf    implicit none
cf    include 'sizes.i'
cf    include 'angle.i'
cf    include 'atoms.i'
cf    include 'bond.i'
cf    include 'disgeo.i'
cf    include 'tors.i'
cf    integer i,ia,ib,ic,id
cf    real*8 locerr,error
cf    real*8 scale,chain,term
cf    real*8 dx,dy,dz,gx,gy,gz
cf    real*8 dstsq,bupsq,blosq
cf    real*8 derivs(3,maxgeo)
c
c
c     zero out the local geometry error function
c
cf    locerr = 0.0d0
cf    scale = 10.0d0
c
c     calculate the bounds error for bond length distances
c
cf    do i = 1, nbond
cf       ia = min(ibnd(1,i),ibnd(2,i))
cf       ib = max(ibnd(1,i),ibnd(2,i))
cf       dx = x(ia) - x(ib)
cf       dy = y(ia) - y(ib)
cf       dz = z(ia) - z(ib)
cf       dstsq = dx*dx + dy*dy + dz*dz
cf       bupsq = bnd(ia,ib)
cf       blosq = bnd(ib,ia)
cf       if (dstsq .gt. bupsq) then
cf          term = (dstsq-bupsq) / bupsq
cf          chain = 4.0d0 * scale * term / bupsq
cf          error = scale * term**2
cf          gx = dx * chain
cf          gy = dy * chain
cf          gz = dz * chain
cf          locerr = locerr + error
cf          derivs(1,ia) = derivs(1,ia) + gx
cf          derivs(2,ia) = derivs(2,ia) + gy
cf          derivs(3,ia) = derivs(3,ia) + gz
cf          derivs(1,ib) = derivs(1,ib) - gx
cf          derivs(2,ib) = derivs(2,ib) - gy
cf          derivs(3,ib) = derivs(3,ib) - gz
cf       else if (dstsq .lt. blosq) then
cf          term = (blosq-dstsq) / (blosq+dstsq)
cf          chain = -8.0d0 * scale * term * (blosq/(blosq+dstsq)**2)
cf          error = scale * term**2
cf          gx = dx * chain
cf          gy = dy * chain
cf          gz = dz * chain
cf          locerr = locerr + error
cf          derivs(1,ia) = derivs(1,ia) + gx
cf          derivs(2,ia) = derivs(2,ia) + gy
cf          derivs(3,ia) = derivs(3,ia) + gz
cf          derivs(1,ib) = derivs(1,ib) - gx
cf          derivs(2,ib) = derivs(2,ib) - gy
cf          derivs(3,ib) = derivs(3,ib) - gz
cf       end if
cf    end do
c
c     calculate the bounds error for the bond angle distances
c
cf    do i = 1, nangle
cf       ia = min(iang(1,i),iang(3,i))
cf       ic = max(iang(1,i),iang(3,i))
cf       dx = x(ia) - x(ic)
cf       dy = y(ia) - y(ic)
cf       dz = z(ia) - z(ic)
cf       dstsq = dx*dx + dy*dy + dz*dz
cf       bupsq = bnd(ia,ic)
cf       blosq = bnd(ic,ia)
cf       if (dstsq .gt. bupsq) then
cf          term = (dstsq-bupsq) / bupsq
cf          chain = 4.0d0 * scale * term / bupsq
cf          error = scale * term**2
cf          gx = dx * chain
cf          gy = dy * chain
cf          gz = dz * chain
cf          locerr = locerr + error
cf          derivs(1,ia) = derivs(1,ia) + gx
cf          derivs(2,ia) = derivs(2,ia) + gy
cf          derivs(3,ia) = derivs(3,ia) + gz
cf          derivs(1,ic) = derivs(1,ic) - gx
cf          derivs(2,ic) = derivs(2,ic) - gy
cf          derivs(3,ic) = derivs(3,ic) - gz
cf       else if (dstsq .lt. blosq) then
cf          term = (blosq-dstsq) / (blosq+dstsq)
cf          chain = -8.0d0 * scale * term * (blosq/(blosq+dstsq)**2)
cf          error = scale * term**2
cf          gx = dx * chain
cf          gy = dy * chain
cf          gz = dz * chain
cf          locerr = locerr + error
cf          derivs(1,ia) = derivs(1,ia) + gx
cf          derivs(2,ia) = derivs(2,ia) + gy
cf          derivs(3,ia) = derivs(3,ia) + gz
cf          derivs(1,ic) = derivs(1,ic) - gx
cf          derivs(2,ic) = derivs(2,ic) - gy
cf          derivs(3,ic) = derivs(3,ic) - gz
cf       end if
cf    end do
c
c     calculate the bounds error for the torsion angle distances
c
cf    do i = 1, ntors
cf       ia = min(itors(1,i),itors(4,i))
cf       id = max(itors(1,i),itors(4,i))
cf       dx = x(ia) - x(id)
cf       dy = y(ia) - y(id)
cf       dz = z(ia) - z(id)
cf       dstsq = dx*dx + dy*dy + dz*dz
cf       bupsq = bnd(ia,id)
cf       blosq = bnd(id,ia)
cf       if (dstsq .gt. bupsq) then
cf          term = (dstsq-bupsq) / bupsq
cf          chain = 4.0d0 * scale * term / bupsq
cf          error = scale * term**2
cf          gx = dx * chain
cf          gy = dy * chain
cf          gz = dz * chain
cf          locerr = locerr + error
cf          derivs(1,ia) = derivs(1,ia) + gx
cf          derivs(2,ia) = derivs(2,ia) + gy
cf          derivs(3,ia) = derivs(3,ia) + gz
cf          derivs(1,id) = derivs(1,id) - gx
cf          derivs(2,id) = derivs(2,id) - gy
cf          derivs(3,id) = derivs(3,id) - gz
cf       else if (dstsq .lt. blosq) then
cf          term = (blosq-dstsq) / (blosq+dstsq)
cf          chain = -8.0d0 * scale * term * (blosq/(blosq+dstsq)**2)
cf          error = scale * term**2
cf          gx = dx * chain
cf          gy = dy * chain
cf          gz = dz * chain
cf          locerr = locerr + error
cf          derivs(1,ia) = derivs(1,ia) + gx
cf          derivs(2,ia) = derivs(2,ia) + gy
cf          derivs(3,ia) = derivs(3,ia) + gz
cf          derivs(1,id) = derivs(1,id) - gx
cf          derivs(2,id) = derivs(2,id) - gy
cf          derivs(3,id) = derivs(3,id) - gz
cf       end if
cf    end do
cf    return
cf    end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function chirer  --  computes chirality error function  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "chirer" computes the chirality error and its derivatives
c     with respect to atomic Cartesian coordinates as a sum the
c     squares of deviations of chiral volumes from target values
c
c
cf    function chirer (derivs)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    integer i,ia,ib,ic,id
cf    real*8 chirer,error,scale
cf    real*8 vol,dv,dedv
cf    real*8 xad,yad,zad
cf    real*8 xbd,ybd,zbd
cf    real*8 xcd,ycd,zcd
cf    real*8 c1,c2,c3
cf    real*8 dedxia,dedyia,dedzia
cf    real*8 dedxib,dedyib,dedzib
cf    real*8 dedxic,dedyic,dedzic
cf    real*8 dedxid,dedyid,dedzid
cf    real*8 derivs(3,maxgeo)
c
c
c     zero the chirality constraint error function
c
cf    chirer = 0.0d0
cf    scale = 0.1d0
c
c     find signed volume value and compute chirality error
c
cf    do i = 1, nchir
cf       ia = ichir(1,i)
cf       ib = ichir(2,i)
cf       ic = ichir(3,i)
cf       id = ichir(4,i)
cf       xad = x(ia) - x(id)
cf       yad = y(ia) - y(id)
cf       zad = z(ia) - z(id)
cf       xbd = x(ib) - x(id)
cf       ybd = y(ib) - y(id)
cf       zbd = z(ib) - z(id)
cf       xcd = x(ic) - x(id)
cf       ycd = y(ic) - y(id)
cf       zcd = z(ic) - z(id)
cf       c1 = ybd*zcd - zbd*ycd
cf       c2 = ycd*zad - zcd*yad
cf       c3 = yad*zbd - zad*ybd
cf       vol = xad*c1 + xbd*c2 + xcd*c3
cf       dv = vol - vchir(i)
cf       error = scale * dv**2
cf       dedv = 2.0d0 * scale * dv
c
c     chain rule terms for first derivative components
c
cf       dedxia = dedv * (ybd*zcd - zbd*ycd)
cf       dedyia = dedv * (zbd*xcd - xbd*zcd)
cf       dedzia = dedv * (xbd*ycd - ybd*xcd)
cf       dedxib = dedv * (zad*ycd - yad*zcd)
cf       dedyib = dedv * (xad*zcd - zad*xcd)
cf       dedzib = dedv * (yad*xcd - xad*ycd)
cf       dedxic = dedv * (yad*zbd - zad*ybd)
cf       dedyic = dedv * (zad*xbd - xad*zbd)
cf       dedzic = dedv * (xad*ybd - yad*xbd)
cf       dedxid = -dedxia - dedxib - dedxic
cf       dedyid = -dedyia - dedyib - dedyic
cf       dedzid = -dedzia - dedzib - dedzic
c
c     increment the chirality constraint error and derivatives
c
cf       chirer = chirer + error
cf       derivs(1,ia) = derivs(1,ia) + dedxia
cf       derivs(2,ia) = derivs(2,ia) + dedyia
cf       derivs(3,ia) = derivs(3,ia) + dedzia
cf       derivs(1,ib) = derivs(1,ib) + dedxib
cf       derivs(2,ib) = derivs(2,ib) + dedyib
cf       derivs(3,ib) = derivs(3,ib) + dedzib
cf       derivs(1,ic) = derivs(1,ic) + dedxic
cf       derivs(2,ic) = derivs(2,ic) + dedyic
cf       derivs(3,ic) = derivs(3,ic) + dedzic
cf       derivs(1,id) = derivs(1,id) + dedxid
cf       derivs(2,id) = derivs(2,id) + dedyid
cf       derivs(3,id) = derivs(3,id) + dedzid
cf    end do
cf    return
cf    end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function torser  --  computes torsional error function  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torser" computes the torsional error function and its first
c     derivatives with respect to the atomic Cartesian coordinates
c     based on the deviation of specified torsional angles from
c     desired values, the contained bond angles are also constrained
c     to avoid a numerical instability
c
c
cf    function torser (derivs)
cf    implicit none
cf    include 'sizes.i'
cf    include 'atoms.i'
cf    include 'disgeo.i'
cf    include 'math.i'
cf    include 'restrn.i'
cf    integer i,ia,ib,ic,id
cf    real*8 torser,derivs(3,maxgeo)
cf    real*8 error,dt,deddt,dedphi
cf    real*8 angle,target,scale
cf    real*8 xia,yia,zia,xib,yib,zib
cf    real*8 xic,yic,zic,xid,yid,zid
cf    real*8 rba,rcb,rdc,dot,cosine,sine
cf    real*8 xp,yp,zp,rp
cf    real*8 terma,termb,termc,termd
cf    real*8 angmax,angmin,cosmax,cosmin
cf    real*8 bamax,bamin,cbmax,cbmin,dcmax,dcmin
cf    real*8 camax,camin,dbmax,dbmin
cf    real*8 xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
cf    real*8 xca,yca,zca,xdb,ydb,zdb
cf    real*8 xt,yt,zt,rt2,xu,yu,zu,ru2
cf    real*8 xtu,ytu,ztu,rtru
cf    real*8 tf1,tf2,t1,t2
cf    real*8 dphidxt,dphidyt,dphidzt
cf    real*8 dphidxu,dphidyu,dphidzu
cf    real*8 dphidxia,dphidyia,dphidzia
cf    real*8 dphidxib,dphidyib,dphidzib
cf    real*8 dphidxic,dphidyic,dphidzic
cf    real*8 dphidxid,dphidyid,dphidzid
cf    real*8 dedxia,dedyia,dedzia
cf    real*8 dedxib,dedyib,dedzib
cf    real*8 dedxic,dedyic,dedzic
cf    real*8 dedxid,dedyid,dedzid
c
c
c     zero the torsional constraint error function
c
cf    torser = 0.0d0
cf    scale = 0.01d0
c
c     compute error value and derivs for torsional constraints
c
cf    do i = 1, ntfix
cf       ia = itfix(1,i)
cf       ib = itfix(2,i)
cf       ic = itfix(3,i)
cf       id = itfix(4,i)
cf       xia = x(ia)
cf       yia = y(ia)
cf       zia = z(ia)
cf       xib = x(ib)
cf       yib = y(ib)
cf       zib = z(ib)
cf       xic = x(ic)
cf       yic = y(ic)
cf       zic = z(ic)
cf       xid = x(id)
cf       yid = y(id)
cf       zid = z(id)
cf       xba = xib - xia
cf       yba = yib - yia
cf       zba = zib - zia
cf       xcb = xic - xib
cf       ycb = yic - yib
cf       zcb = zic - zib
cf       xdc = xid - xic
cf       ydc = yid - yic
cf       zdc = zid - zic
c
c     find the actual distances between the four atoms
c
cf       rba = sqrt(xba*xba + yba*yba + zba*zba)
cf       rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
cf       rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
c
c     get the minimum and minimum values for distances
c
cf       bamax = sqrt(bnd(min(ib,ia),max(ib,ia)))
cf       bamin = sqrt(bnd(max(ib,ia),min(ib,ia)))
cf       cbmax = sqrt(bnd(min(ic,ib),max(ic,ib)))
cf       cbmin = sqrt(bnd(max(ic,ib),min(ic,ib)))
cf       dcmax = sqrt(bnd(min(id,ic),max(id,ic)))
cf       dcmin = sqrt(bnd(max(id,ic),min(id,ic)))
cf       camax = sqrt(bnd(min(ic,ia),max(ic,ia)))
cf       camin = sqrt(bnd(max(ic,ia),min(ic,ia)))
cf       dbmax = sqrt(bnd(min(id,ib),max(id,ib)))
cf       dbmin = sqrt(bnd(max(id,ib),min(id,ib)))
c
c     compute the ia-ib-ic bond angle and any error
c
cf       dot = xba*xcb + yba*ycb + zba*zcb
cf       cosine = -dot / (rba*rcb)
cf       cosine = min(1.0d0,max(-1.0d0,cosine))
cf       angle = radian * acos(cosine)
cf       cosmax = (bamin**2+cbmin**2-camax**2) / (2.0d0*bamin*cbmin)
cf       cosmax = min(1.0d0,max(-1.0d0,cosmax))
cf       angmax = radian * acos(cosmax)
cf       cosmin = (bamax**2+cbmax**2-camin**2) / (2.0d0*bamax*cbmax)
cf       cosmin = min(1.0d0,max(-1.0d0,cosmin))
cf       angmin = radian * acos(cosmin)
cf       if (angle .gt. angmax) then
cf          dt = angle - angmax
cf       else if (angle .lt. angmin) then
cf          dt = angle - angmin
cf       else
cf          dt = 0.0d0
cf       end if
cf       error = scale * dt**2
cf       deddt = 2.0d0 * radian * scale * dt
c
c     compute derivative components for this interaction
c
cf       xp = zcb*yba - ycb*zba
cf       yp = xcb*zba - zcb*xba
cf       zp = ycb*xba - xcb*yba
cf       rp = sqrt(xp*xp + yp*yp + zp*zp)
cf       if (rp .ne. 0.0d0) then
cf          terma = -deddt / (rba*rba*rp)
cf          termc = deddt / (rcb*rcb*rp)
cf          dedxia = terma * (zba*yp-yba*zp)
cf          dedyia = terma * (xba*zp-zba*xp)
cf          dedzia = terma * (yba*xp-xba*yp)
cf          dedxic = termc * (ycb*zp-zcb*yp)
cf          dedyic = termc * (zcb*xp-xcb*zp)
cf          dedzic = termc * (xcb*yp-ycb*xp)
cf          dedxib = -dedxia - dedxic
cf          dedyib = -dedyia - dedyic
cf          dedzib = -dedzia - dedzic
c
c     increment the bond angle constraint error and derivatives
c
cf          torser = torser + error
cf          derivs(1,ia) = derivs(1,ia) + dedxia
cf          derivs(2,ia) = derivs(2,ia) + dedyia
cf          derivs(3,ia) = derivs(3,ia) + dedzia
cf          derivs(1,ib) = derivs(1,ib) + dedxib
cf          derivs(2,ib) = derivs(2,ib) + dedyib
cf          derivs(3,ib) = derivs(3,ib) + dedzib
cf          derivs(1,ic) = derivs(1,ic) + dedxic
cf          derivs(2,ic) = derivs(2,ic) + dedyic
cf          derivs(3,ic) = derivs(3,ic) + dedzic
cf       end if
c
c     compute the ib-ic-id bond angle and any error
c
cf       dot = xdc*xcb + ydc*ycb + zdc*zcb
cf       cosine = -dot / (rdc*rcb)
cf       cosine = min(1.0d0,max(-1.0d0,cosine))
cf       angle = radian * acos(cosine)
cf       cosmax = (dcmin**2+cbmin**2-dbmax**2) / (2.0d0*dcmin*cbmin)
cf       cosmax = min(1.0d0,max(-1.0d0,cosmax))
cf       angmax = radian * acos(cosmax)
cf       cosmin = (dcmax**2+cbmax**2-dbmin**2) / (2.0d0*dcmax*cbmax)
cf       cosmax = min(1.0d0,max(-1.0d0,cosmin))
cf       angmin = radian * acos(cosmin)
cf       if (angle .gt. angmax) then
cf          dt = angle - angmax
cf       else if (angle .lt. angmin) then
cf          dt = angle - angmin
cf       else
cf          dt = 0.0d0
cf       end if
cf       error = scale * dt**2
cf       deddt = 2.0d0 * radian * scale * dt
c
c     compute derivative components for this interaction
c
cf       xp = zdc*ycb - ydc*zcb
cf       yp = xdc*zcb - zdc*xcb
cf       zp = ydc*xcb - xdc*ycb
cf       rp = sqrt(xp*xp + yp*yp + zp*zp)
cf       if (rp .ne. 0.0d0) then
cf          termb = -deddt / (rcb*rcb*rp)
cf          termd = deddt / (rdc*rdc*rp)
cf          dedxib = termb * (zcb*yp-ycb*zp)
cf          dedyib = termb * (xcb*zp-zcb*xp)
cf          dedzib = termb * (ycb*xp-xcb*yp)
cf          dedxid = termd * (ydc*zp-zdc*yp)
cf          dedyid = termd * (zdc*xp-xdc*zp)
cf          dedzid = termd * (xdc*yp-ydc*xp)
cf          dedxic = -dedxib - dedxid
cf          dedyic = -dedyib - dedyid
cf          dedzic = -dedzib - dedzid
c
c     increment the bond angle constraint error and derivatives
c
cf          torser = torser + error
cf          derivs(1,ib) = derivs(1,ib) + dedxib
cf          derivs(2,ib) = derivs(2,ib) + dedyib
cf          derivs(3,ib) = derivs(3,ib) + dedzib
cf          derivs(1,ic) = derivs(1,ic) + dedxic
cf          derivs(2,ic) = derivs(2,ic) + dedyic
cf          derivs(3,ic) = derivs(3,ic) + dedzic
cf          derivs(1,id) = derivs(1,id) + dedxid
cf          derivs(2,id) = derivs(2,id) + dedyid
cf          derivs(3,id) = derivs(3,id) + dedzid
cf       end if
c
c     compute the value of the ia-ib-ic-id torsional angle
c
cf       xt = yba*zcb - ycb*zba
cf       yt = zba*xcb - zcb*xba
cf       zt = xba*ycb - xcb*yba
cf       xu = ycb*zdc - ydc*zcb
cf       yu = zcb*xdc - zdc*xcb
cf       zu = xcb*ydc - xdc*ycb
cf       xtu = yt*zu - yu*zt
cf       ytu = zt*xu - zu*xt
cf       ztu = xt*yu - xu*yt
cf       rt2 = xt*xt + yt*yt + zt*zt
cf       ru2 = xu*xu + yu*yu + zu*zu
cf       rtru = sqrt(rt2 * ru2)
cf       if (rtru .ne. 0.0d0) then
cf          rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
cf          cosine = (xt*xu + yt*yu + zt*zu) / rtru
cf          sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
cf          cosine = min(1.0d0,max(-1.0d0,cosine))
cf          angle = radian * acos(cosine)
cf          if (sine .lt. 0.0d0)  angle = -angle
c
c     calculate the torsional constraint error for this angle
c
cf          tf1 = tfix(1,i)
cf          tf2 = tfix(2,i)
cf          if (angle.gt.tf1 .and. angle.lt.tf2) then
cf             target = angle
cf          else if (angle.gt.tf1 .and. tf1.gt.tf2) then
cf             target = angle
cf          else if (angle.lt.tf2 .and. tf1.gt.tf2) then
cf             target = angle
cf          else
cf             t1 = angle - tf1
cf             t2 = angle - tf2
cf             if (t1 .gt. 180.0d0) then
cf                t1 = t1 - 360.0d0
cf             else if (t1 .lt. -180.0d0) then
cf                t1 = t1 + 360.0d0
cf             end if
cf             if (t2 .gt. 180.0d0) then
cf                t2 = t2 - 360.0d0
cf             else if (t2 .lt. -180.0d0) then
cf                t2 = t2 + 360.0d0
cf             end if
cf             if (abs(t1) .lt. abs(t2)) then
cf                target = tf1
cf             else
cf                target = tf2
cf             end if
cf          end if
cf          dt = angle - target
cf          if (dt .gt. 180.0d0) then
cf             dt = dt - 360.0d0
cf          else if (dt .lt. -180.0d0) then
cf             dt = dt + 360.0d0
cf          end if
cf          error = scale * dt**2
cf          dedphi = 2.0d0 * radian * scale * dt
c
c     abbreviations for first derivative chain rule terms
c
cf          xca = xic - xia
cf          yca = yic - yia
cf          zca = zic - zia
cf          xdb = xid - xib
cf          ydb = yid - yib
cf          zdb = zid - zib
cf          dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
cf          dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
cf          dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
cf          dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
cf          dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
cf          dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     chain rule terms for first derivative components
c
cf          dphidxia = zcb*dphidyt - ycb*dphidzt
cf          dphidyia = xcb*dphidzt - zcb*dphidxt
cf          dphidzia = ycb*dphidxt - xcb*dphidyt
cf          dphidxib = yca*dphidzt - zca*dphidyt
cf   &                    + zdc*dphidyu - ydc*dphidzu
cf          dphidyib = zca*dphidxt - xca*dphidzt
cf   &                    + xdc*dphidzu - zdc*dphidxu
cf          dphidzib = xca*dphidyt - yca*dphidxt
cf   &                    + ydc*dphidxu - xdc*dphidyu
cf          dphidxic = zba*dphidyt - yba*dphidzt
cf   &                    + ydb*dphidzu - zdb*dphidyu
cf          dphidyic = xba*dphidzt - zba*dphidxt
cf   &                    + zdb*dphidxu - xdb*dphidzu
cf          dphidzic = yba*dphidxt - xba*dphidyt
cf   &                    + xdb*dphidyu - ydb*dphidxu
cf          dphidxid = zcb*dphidyu - ycb*dphidzu
cf          dphidyid = xcb*dphidzu - zcb*dphidxu
cf          dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     compute first derivative components for torsion angle
c
cf          dedxia = dedphi * dphidxia
cf          dedyia = dedphi * dphidyia
cf          dedzia = dedphi * dphidzia
cf          dedxib = dedphi * dphidxib
cf          dedyib = dedphi * dphidyib
cf          dedzib = dedphi * dphidzib
cf          dedxic = dedphi * dphidxic
cf          dedyic = dedphi * dphidyic
cf          dedzic = dedphi * dphidzic
cf          dedxid = dedphi * dphidxid
cf          dedyid = dedphi * dphidyid
cf          dedzid = dedphi * dphidzid
c
c     increment the torsional constraint error and derivatives
c
cf          torser = torser + error
cf          derivs(1,ia) = derivs(1,ia) + dedxia
cf          derivs(2,ia) = derivs(2,ia) + dedyia
cf          derivs(3,ia) = derivs(3,ia) + dedzia
cf          derivs(1,ib) = derivs(1,ib) + dedxib
cf          derivs(2,ib) = derivs(2,ib) + dedyib
cf          derivs(3,ib) = derivs(3,ib) + dedzib
cf          derivs(1,ic) = derivs(1,ic) + dedxic
cf          derivs(2,ic) = derivs(2,ic) + dedyic
cf          derivs(3,ic) = derivs(3,ic) + dedzic
cf          derivs(1,id) = derivs(1,id) + dedxid
cf          derivs(2,id) = derivs(2,id) + dedyid
cf          derivs(3,id) = derivs(3,id) + dedzid
cf       end if
cf    end do
cf    return
cf    end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emm3hb  --  MM3 van der Waals and hbond energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emm3hb" calculates the van der Waals and hydrogen bonding
c     interaction energy using the MM3 exp-6 formula with directional
c     charge transfer hydrogen bonding
c
c     literature reference:
c
c     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in
c     the MM3 Force Field. I", Journal of Physical Organic Chemistry,
c     7, 591-609 (1994)
c
c
      subroutine emm3hb
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,ii,kk
      integer ia,ib,ic,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expcut2
      real*8 expterm,expmerge
      real*8 dot,cosine
      real*8 fterm,ideal
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical proceed,iuse
c
c
c     zero out the van der Waals energy contributions
c
      ev = 0.0d0
      e14 = 0.0d0
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find the van der Waals energy via double loop search
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
               if (rik2 .le. off2) then
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i) then
                     eps = eps / vdwscale
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = k
                     else
                        ia = k
                        ib = i12(1,k)
                        ic = i
                     end if
                     xia = x(ia)
                     yia = y(ia)
                     zia = z(ia)
                     xib = x(ib)
                     yib = y(ib)
                     zib = z(ib)
                     xic = x(ic)
                     yic = y(ic)
                     zic = z(ic)
                     xab = xia - xib
                     yab = yia - yib
                     zab = zia - zib
                     rab2 = xab*xab + yab*yab + zab*zab
                     rab = sqrt(rab2)
                     xcb = xic - xib
                     ycb = yic - yib
                     zcb = zic - zib
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     ideal = bl(bndlist(1,ia))
                     fterm = cosine * (rab/ideal)
                  end if
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - fterm*cterm*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (skip(k) .eq. -i) then
                     e14 = e14 + e
                  else
                     ev = ev + e
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  kt = class(k)
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
                  if (rik2 .le. off2) then
                     fterm = 1.0d0
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (radhbnd(kt,it) .ne. 0.0d0) then
                        rv = radhbnd(kt,it)
                        eps = epshbnd(kt,it) / dielec
                        if (atomic(i) .eq. 1) then
                           ia = i
                           ib = i12(1,i)
                           ic = k
                        else
                           ia = k
                           ib = i12(1,k)
                           ic = i
                        end if
                        xia = x(ia)
                        yia = y(ia)
                        zia = z(ia)
                        xib = x(ib)
                        yib = y(ib)
                        zib = z(ib)
                        xic = x(ic)
                        yic = y(ic)
                        zic = z(ic)
                        xab = xia - xib
                        yab = yia - yib
                        zab = zia - zib
                        rab2 = xab*xab + yab*yab + zab*zab
                        rab = sqrt(rab2)
                        xcb = xic - xib
                        ycb = yic - yib
                        zcb = zic - zib
                        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                        xp = ycb*zab - zcb*yab
                        yp = zcb*xab - xcb*zab
                        zp = xcb*yab - ycb*xab
                        rp = sqrt(xp*xp + yp*yp + zp*zp)
                        dot = xab*xcb + yab*ycb + zab*zcb
                        cosine = dot / sqrt(rab2*rcb2)
                        ideal = bl(bndlist(1,ia))
                        fterm = cosine * (rab/ideal)
                     end if
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        expterm = aterm * exp(-bterm/p)
                        e = eps * (expterm - fterm*cterm*p6)
                     else
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik = sqrt(rik2)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy component;
c     interaction of an atom with its own image counts half
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine emm3hb1  --  MM3 vdw & hbond energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "emm3hb1" calculates the van der Waals and hydrogen bonding
c     energy and its first derivatives with respect to Cartesian
c     coordinates using the MM3 exp-6 formula with directional
c     charge transfer hydrogen bonding
c
c     literature reference:
c
c     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in
c     the MM3 Force Field. I", Journal of Physical Organic Chemistry,
c     7, 591-609 (1994)
c
c
      subroutine emm3hb1
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bond.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,k,ii,kk
      integer ia,ib,ic,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv,redk,redkv,rvterm
      real*8 dedx,dedy,dedz,de
      real*8 rik,rik2,rik3,rik4,rik5
      real*8 taper,dtaper
      real*8 expcut,expcut2
      real*8 expterm,expmerge
      real*8 dot,cosine,sine
      real*8 fterm,fcterm,term
      real*8 deddr,ideal,ratio
      real*8 deddt,terma,termc
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical proceed,iuse,use_hb
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      e14 = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
         de14(1,i) = 0.0d0
         de14(2,i) = 0.0d0
         de14(3,i) = 0.0d0
      end do
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find van der Waals energy and derivatives via double loop
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute energy and derivatives for this interaction
c
               if (rik2 .le. off2) then
                  use_hb = .false.
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i) then
                     eps = eps / vdwscale
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     use_hb = .true.
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = k
                     else
                        ia = k
                        ib = i12(1,k)
                        ic = i
                     end if
                     xia = x(ia)
                     yia = y(ia)
                     zia = z(ia)
                     xib = x(ib)
                     yib = y(ib)
                     zib = z(ib)
                     xic = x(ic)
                     yic = y(ic)
                     zic = z(ic)
                     xab = xia - xib
                     yab = yia - yib
                     zab = zia - zib
                     rab2 = xab*xab + yab*yab + zab*zab
                     rab = sqrt(rab2)
                     xcb = xic - xib
                     ycb = yic - yib
                     zcb = zic - zib
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     sine = sqrt(abs(1.0d0-cosine**2))
                     ideal = bl(bndlist(1,ia))
                     ratio = rab / ideal
                     fterm = cosine * ratio
                     deddt = -sine * ratio
                     deddr = cosine / (rab*ideal)
                  end if
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bterm / rv
                     expterm = aterm * exp(-bterm/p)
                     fcterm = fterm * cterm * p6
                     e = eps * (expterm - fcterm)
                     de = eps * (rvterm*expterm+6.0d0*fcterm/rik)
                  else
                     use_hb = .false.
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     if (use_hb) then
                        deddt = deddt * fgrp
                        deddr = deddr * fgrp
                     end if
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  if (skip(k) .eq. -i) then
                     e14 = e14 + e
                     if (i .eq. iv) then
                        de14(1,i) = de14(1,i) + dedx
                        de14(2,i) = de14(2,i) + dedy
                        de14(3,i) = de14(3,i) + dedz
                     else
                        de14(1,i) = de14(1,i) + dedx*redi
                        de14(2,i) = de14(2,i) + dedy*redi
                        de14(3,i) = de14(3,i) + dedz*redi
                        de14(1,iv) = de14(1,iv) + dedx*rediv
                        de14(2,iv) = de14(2,iv) + dedy*rediv
                        de14(3,iv) = de14(3,iv) + dedz*rediv
                     end if
                     if (k .eq. kv) then
                        de14(1,k) = de14(1,k) - dedx
                        de14(2,k) = de14(2,k) - dedy
                        de14(3,k) = de14(3,k) - dedz
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        de14(1,k) = de14(1,k) - dedx*redk
                        de14(2,k) = de14(2,k) - dedy*redk
                        de14(3,k) = de14(3,k) - dedz*redk
                        de14(1,kv) = de14(1,kv) - dedx*redkv
                        de14(2,kv) = de14(2,kv) - dedy*redkv
                        de14(3,kv) = de14(3,kv) - dedz*redkv
                     end if
                  else
                     ev = ev + e
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                     end if
                     if (k .eq. kv) then
                        dev(1,k) = dev(1,k) - dedx
                        dev(2,k) = dev(2,k) - dedy
                        dev(3,k) = dev(3,k) - dedz
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        dev(1,k) = dev(1,k) - dedx*redk
                        dev(2,k) = dev(2,k) - dedy*redk
                        dev(3,k) = dev(3,k) - dedz*redk
                        dev(1,kv) = dev(1,kv) - dedx*redkv
                        dev(2,kv) = dev(2,kv) - dedy*redkv
                        dev(3,kv) = dev(3,kv) - dedz*redkv
                     end if
                  end if
c
c     find the chain rule terms for hydrogen bonding components
c
                  if (use_hb) then
                     term = eps * cterm * p6
                     deddt = deddt * term
                     deddr = deddr * term
                     if (rik2 .gt. cut2) then
                        deddt = deddt * taper
                        deddr = deddr * taper
                     end if
                     terma = deddt / (rab2*rp)
                     termc = -deddt / (rcb2*rp)
                     dedxia = terma * (yab*zp-zab*yp) - deddr*xab
                     dedyia = terma * (zab*xp-xab*zp) - deddr*yab
                     dedzia = terma * (xab*yp-yab*xp) - deddr*zab
                     dedxic = termc * (ycb*zp-zcb*yp)
                     dedyic = termc * (zcb*xp-xcb*zp)
                     dedzic = termc * (xcb*yp-ycb*xp)
                     dedxib = -dedxia - dedxic
                     dedyib = -dedyia - dedyic
                     dedzib = -dedzia - dedzic
c
c     increment the derivatives for directional hydrogen bonding
c
                     dev(1,ia) = dev(1,ia) + dedxia
                     dev(2,ia) = dev(2,ia) + dedyia
                     dev(3,ia) = dev(3,ia) + dedzia
                     dev(1,ib) = dev(1,ib) + dedxib
                     dev(2,ib) = dev(2,ib) + dedyib
                     dev(3,ib) = dev(3,ib) + dedzib
                     dev(1,ic) = dev(1,ic) + dedxic
                     dev(2,ic) = dev(2,ic) + dedyic
                     dev(3,ic) = dev(3,ic) + dedzic
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx + xr*dedx
                     viry = viry + yr*dedy
                     virz = virz + zr*dedz
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  kt = class(k)
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     compute energy and derivatives for this interaction
c
                  if (rik2 .le. off2) then
                     use_hb = .false.
                     fterm = 1.0d0
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (radhbnd(kt,it) .ne. 0.0d0) then
                        use_hb = .true.
                        rv = radhbnd(kt,it)
                        eps = epshbnd(kt,it) / dielec
                        if (atomic(i) .eq. 1) then
                           ia = i
                           ib = i12(1,i)
                           ic = k
                        else
                           ia = k
                           ib = i12(1,k)
                           ic = i
                        end if
                        xia = x(ia)
                        yia = y(ia)
                        zia = z(ia)
                        xib = x(ib)
                        yib = y(ib)
                        zib = z(ib)
                        xic = x(ic)
                        yic = y(ic)
                        zic = z(ic)
                        xab = xia - xib
                        yab = yia - yib
                        zab = zia - zib
                        rab2 = xab*xab + yab*yab + zab*zab
                        rab = sqrt(rab2)
                        xcb = xic - xib
                        ycb = yic - yib
                        zcb = zic - zib
                        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                        xp = ycb*zab - zcb*yab
                        yp = zcb*xab - xcb*zab
                        zp = xcb*yab - ycb*xab
                        rp = sqrt(xp*xp + yp*yp + zp*zp)
                        dot = xab*xcb + yab*ycb + zab*zcb
                        cosine = dot / sqrt(rab2*rcb2)
                        sine = sqrt(abs(1.0d0-cosine**2))
                        ideal = bl(bndlist(1,ia))
                        ratio = rab / ideal
                        fterm = cosine * ratio
                        deddt = -sine * ratio
                        deddr = cosine / (rab*ideal)
                     end if
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     rik = sqrt(rik2)
                     if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        rvterm = -bterm / rv
                        expterm = aterm * exp(-bterm/p)
                        fcterm = fterm * cterm * p6
                        e = eps * (expterm - fcterm)
                        de = eps * (rvterm*expterm+6.0d0*fcterm/rik)
                     else
                        use_hb = .false.
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                        de = -12.0d0 * e / rik
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                              + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        de = e*dtaper + de*taper
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                        if (use_hb) then
                           deddt = deddt * fgrp
                           deddr = deddr * fgrp
                        end if
                     end if
c
c     find the chain rule terms for derivative components
c
                     de = de / rik
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                     end if
                     if (i .ne. k) then
                        if (k .eq. kv) then
                           dev(1,k) = dev(1,k) - dedx
                           dev(2,k) = dev(2,k) - dedy
                           dev(3,k) = dev(3,k) - dedz
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           dev(1,k) = dev(1,k) - dedx*redk
                           dev(2,k) = dev(2,k) - dedy*redk
                           dev(3,k) = dev(3,k) - dedz*redk
                           dev(1,kv) = dev(1,kv) - dedx*redkv
                           dev(2,kv) = dev(2,kv) - dedy*redkv
                           dev(3,kv) = dev(3,kv) - dedz*redkv
                        end if
                     end if
c
c     find the chain rule terms for hydrogen bonding components
c
                     if (use_hb) then
                        term = eps * cterm * p6
                        deddt = deddt * term
                        deddr = deddr * term
                        if (rik2 .gt. cut2) then
                           deddt = deddt * taper
                           deddr = deddr * taper
                        end if
                        terma = deddt / (rab2*rp)
                        termc = -deddt / (rcb2*rp)
                        dedxia = terma * (yab*zp-zab*yp) - deddr*xab
                        dedyia = terma * (zab*xp-xab*zp) - deddr*yab
                        dedzia = terma * (xab*yp-yab*xp) - deddr*zab
                        dedxic = termc * (ycb*zp-zcb*yp)
                        dedyic = termc * (zcb*xp-xcb*zp)
                        dedzic = termc * (xcb*yp-ycb*xp)
                        dedxib = -dedxia - dedxic
                        dedyib = -dedyia - dedyic
                        dedzib = -dedzia - dedzic
c
c     increment the derivatives for directional hydrogen bonding
c
                        dev(1,ia) = dev(1,ia) + dedxia
                        dev(2,ia) = dev(2,ia) + dedyia
                        dev(3,ia) = dev(3,ia) + dedzia
                        dev(1,ib) = dev(1,ib) + dedxib
                        dev(2,ib) = dev(2,ib) + dedyib
                        dev(3,ib) = dev(3,ib) + dedzib
                        dev(1,ic) = dev(1,ic) + dedxic
                        dev(2,ic) = dev(2,ic) + dedyic
                        dev(3,ic) = dev(3,ic) + dedzic
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*dedx
                        viry = viry + yr*dedy
                        virz = virz + zr*dedz
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emm3hb2  --  atom-wise MM3 vdw & hbond Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emm3hb2" calculates the van der Waals and hydrogen bonding
c     second derivatives for a single atom at a time using the MM3
c     exp-6 formula with directional charge transfer hydrogen bonding
c
c     note this version only partially incorporates the directional
c     hydrogen bonding term into the Hessian calculation
c
c     literature reference:
c
c     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in
c     the MM3 Force Field. I", Journal of Physical Organic Chemistry,
c     7, 591-609 (1994)
c
c
      subroutine emm3hb2 (iatom,xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'shunt.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer iatom,i,j,k,ii,kk
      integer ia,ib,ic,iv,kv
      integer it,kt,skip(maxatm)
      integer nuse,use(5),jcell
      real*8 e,de,d2e,fgrp
      real*8 p,p2,p6,p12,eps,rv
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 redi2,rediv2,rediiv
      real*8 redik,redivk,redikv,redivkv
      real*8 rik,rik2,rik3,rik4,rik5
      real*8 taper,dtaper,d2taper
      real*8 d2edx,d2edy,d2edz,term(3,3)
      real*8 expcut,expcut2
      real*8 expterm,expmerge
      real*8 rvterm,rvterm2
      real*8 dot,cosine,sine
      real*8 fterm,fcterm
      real*8 ideal,ratio
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical proceed,use_hb
c
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     check to see if the atom of interest is a vdw site
c
      nuse = 0
      do k = 1, nvdw
         if (ivdw(k) .eq. iatom) then
            nuse = nuse + 1
            use(nuse) = iatom
            goto 10
         end if
      end do
      return
   10 continue
c
c     determine the atoms involved via reduction factors
c
      nuse = 1
      use(nuse) = iatom
      do k = 1, n12(iatom)
         i = i12(k,iatom)
         if (ired(i) .eq. iatom) then
            nuse = nuse + 1
            use(nuse) = i
         end if
      end do
c
c     find van der Waals Hessian elements for involved atoms
c
      do ii = 1, nuse
         i = use(ii)
         iv = ired(i)
         redi = kred(i)
         if (i .ne. iv) then
            rediv = 1.0d0 - redi
            redi2 = redi * redi
            rediv2 = rediv * rediv
            rediiv = redi * rediv
         end if
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         skip(i) = i
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         do kk = 1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute Hessian elements for this interaction
c
               if (rik2 .le. off2) then
                  use_hb = .false.
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i) then
                     eps = eps / vdwscale
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     use_hb = .true.
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = k
                     else
                        ia = k
                        ib = i12(1,k)
                        ic = i
                     end if
                     xia = x(ia)
                     yia = y(ia)
                     zia = z(ia)
                     xib = x(ib)
                     yib = y(ib)
                     zib = z(ib)
                     xic = x(ic)
                     yic = y(ic)
                     zic = z(ic)
                     xab = xia - xib
                     yab = yia - yib
                     zab = zia - zib
                     rab2 = xab*xab + yab*yab + zab*zab
                     rab = sqrt(rab2)
                     xcb = xic - xib
                     ycb = yic - yib
                     zcb = zic - zib
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     sine = sqrt(abs(1.0d0-cosine**2))
                     ideal = bl(bndlist(1,ia))
                     ratio = rab / ideal
                     fterm = cosine * ratio
                  end if
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bterm / rv
                     rvterm2 = rvterm * rvterm
                     expterm = aterm * exp(-bterm/p)
                     fcterm = fterm * cterm * p6
                     e = eps * (expterm - fcterm)
                     de = eps * (rvterm*expterm+6.0d0*fcterm/rik)
                     d2e = eps * (rvterm2*expterm-42.0d0*fcterm/rik2)
                  else
                     use_hb = .false.
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                     d2e = 156.0d0 * e / rik2
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     e = eps * (p12 - 2.0d0 * p6)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     d2taper = 20.0d0*c5*rik3 + 12.0d0*c4*rik2
     &                            + 6.0d0*c3*rik + 2.0d0*c2
                     d2e = e*d2taper + 2.0d0*de*dtaper + d2e*taper
                     de = e*dtaper + de*taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     de = de * fgrp
                     d2e = d2e * fgrp
                  end if
c
c     get chain rule terms for van der Waals Hessian elements
c
                  de = de / rik
                  d2e = (d2e-de) / rik2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  if (i .eq. iatom) then
                     if (i.eq.iv .and. k.eq.kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)
                           hessy(j,i) = hessy(j,i) + term(2,j)
                           hessz(j,i) = hessz(j,i) + term(3,j)
                           hessx(j,k) = hessx(j,k) - term(1,j)
                           hessy(j,k) = hessy(j,k) - term(2,j)
                           hessz(j,k) = hessz(j,k) - term(3,j)
                        end do
                     else if (k .eq. kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                           hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                           hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                           hessx(j,k) = hessx(j,k) - term(1,j)*redi
                           hessy(j,k) = hessy(j,k) - term(2,j)*redi
                           hessz(j,k) = hessz(j,k) - term(3,j)*redi
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                        end do
                     else if (i .eq. iv) then
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)
                           hessy(j,i) = hessy(j,i) + term(2,j)
                           hessz(j,i) = hessz(j,i) + term(3,j)
                           hessx(j,k) = hessx(j,k) - term(1,j)*redk
                           hessy(j,k) = hessy(j,k) - term(2,j)*redk
                           hessz(j,k) = hessz(j,k) - term(3,j)*redk
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redkv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redkv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redkv
                        end do
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        redik = redi * redk
                        redikv = redi * redkv
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                           hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                           hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                           hessx(j,k) = hessx(j,k) - term(1,j)*redik
                           hessy(j,k) = hessy(j,k) - term(2,j)*redik
                           hessz(j,k) = hessz(j,k) - term(3,j)*redik
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redikv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redikv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redikv
                        end do
                     end if
                  else if (iv .eq. iatom) then
                     if (k .eq. kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                           hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                           hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                           hessx(j,k) = hessx(j,k) - term(1,j)*rediv
                           hessy(j,k) = hessy(j,k) - term(2,j)*rediv
                           hessz(j,k) = hessz(j,k) - term(3,j)*rediv
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                        end do
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        redivk = rediv * redk
                        redivkv = rediv * redkv
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                           hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                           hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                           hessx(j,k) = hessx(j,k) - term(1,j)*redivk
                           hessy(j,k) = hessy(j,k) - term(2,j)*redivk
                           hessz(j,k) = hessz(j,k) - term(3,j)*redivk
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redivkv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redivkv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redivkv
                        end do
                     end if
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nuse
         i = use(ii)
         iv = ired(i)
         redi = kred(i)
         if (i .ne. iv) then
            rediv = 1.0d0 - redi
            redi2 = redi * redi
            rediv2 = rediv * rediv
            rediiv = redi * rediv
         end if
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do kk = 1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do jcell = 1, ncell
                  kt = class(k)
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,jcell)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     compute Hessian elements for this interaction
c
                  if (rik2 .le. off2) then
                     use_hb = .false.
                     fterm = 1.0d0
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (radhbnd(kt,it) .ne. 0.0d0) then
                        use_hb = .true.
                        rv = radhbnd(kt,it)
                        eps = epshbnd(kt,it) / dielec
                        if (atomic(i) .eq. 1) then
                           ia = i
                           ib = i12(1,i)
                           ic = k
                        else
                           ia = k
                           ib = i12(1,k)
                           ic = i
                        end if
                        xia = x(ia)
                        yia = y(ia)
                        zia = z(ia)
                        xib = x(ib)
                        yib = y(ib)
                        zib = z(ib)
                        xic = x(ic)
                        yic = y(ic)
                        zic = z(ic)
                        xab = xia - xib
                        yab = yia - yib
                        zab = zia - zib
                        rab2 = xab*xab + yab*yab + zab*zab
                        rab = sqrt(rab2)
                        xcb = xic - xib
                        ycb = yic - yib
                        zcb = zic - zib
                        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                        xp = ycb*zab - zcb*yab
                        yp = zcb*xab - xcb*zab
                        zp = xcb*yab - ycb*xab
                        rp = sqrt(xp*xp + yp*yp + zp*zp)
                        dot = xab*xcb + yab*ycb + zab*zcb
                        cosine = dot / sqrt(rab2*rcb2)
                        sine = sqrt(abs(1.0d0-cosine**2))
                        ideal = bl(bndlist(1,ia))
                        ratio = rab / ideal
                        fterm = cosine * ratio
                     end if
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     rik = sqrt(rik2)
                     if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        rvterm = -bterm / rv
                        rvterm2 = rvterm * rvterm
                        expterm = aterm * exp(-bterm/p)
                        fcterm = fterm * cterm * p6
                        e = eps * (expterm - fcterm)
                        de = eps * (rvterm*expterm+6.0d0*fcterm/rik)
                        d2e = eps * (rvterm2*expterm-42.0d0*fcterm/rik2)
                     else
                        use_hb = .false.
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                        de = -12.0d0 * e / rik
                        d2e = 156.0d0 * e / rik2
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        e = eps * (p12 - 2.0d0 * p6)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        d2taper = 20.0d0*c5*rik3 + 12.0d0*c4*rik2
     &                             + 6.0d0*c3*rik + 2.0d0*c2
                        d2e = e*d2taper + 2.0d0*de*dtaper + d2e*taper
                        de = e*dtaper + de*taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        de = de * fgrp
                        d2e = d2e * fgrp
                     end if
c
c     get chain rule terms for van der Waals Hessian elements
c
                     de = de / rik
                     d2e = (d2e-de) / rik2
                     d2edx = d2e * xr
                     d2edy = d2e * yr
                     d2edz = d2e * zr
                     term(1,1) = d2edx*xr + de
                     term(1,2) = d2edx*yr
                     term(1,3) = d2edx*zr
                     term(2,1) = term(1,2)
                     term(2,2) = d2edy*yr + de
                     term(2,3) = d2edy*zr
                     term(3,1) = term(1,3)
                     term(3,2) = term(2,3)
                     term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                     if (i .eq. iatom) then
                        if (i.eq.iv .and. k.eq.kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)
                              hessy(j,i) = hessy(j,i) + term(2,j)
                              hessz(j,i) = hessz(j,i) + term(3,j)
                              hessx(j,k) = hessx(j,k) - term(1,j)
                              hessy(j,k) = hessy(j,k) - term(2,j)
                              hessz(j,k) = hessz(j,k) - term(3,j)
                           end do
                        else if (k .eq. kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                              hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                              hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                              hessx(j,k) = hessx(j,k) - term(1,j)*redi
                              hessy(j,k) = hessy(j,k) - term(2,j)*redi
                              hessz(j,k) = hessz(j,k) - term(3,j)*redi
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediiv
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediiv
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediiv
                           end do
                        else if (i .eq. iv) then
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)
                              hessy(j,i) = hessy(j,i) + term(2,j)
                              hessz(j,i) = hessz(j,i) + term(3,j)
                              hessx(j,k) = hessx(j,k) - term(1,j)*redk
                              hessy(j,k) = hessy(j,k) - term(2,j)*redk
                              hessz(j,k) = hessz(j,k) - term(3,j)*redk
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redkv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redkv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redkv
                           end do
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           redik = redi * redk
                           redikv = redi * redkv
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                              hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                              hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                              hessx(j,k) = hessx(j,k) - term(1,j)*redik
                              hessy(j,k) = hessy(j,k) - term(2,j)*redik
                              hessz(j,k) = hessz(j,k) - term(3,j)*redik
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediiv
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediiv
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediiv
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redikv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redikv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redikv
                           end do
                        end if
                     else if (iv .eq. iatom) then
                        if (k .eq. kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                              hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                              hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                              hessx(j,k) = hessx(j,k) - term(1,j)*rediv
                              hessy(j,k) = hessy(j,k) - term(2,j)*rediv
                              hessz(j,k) = hessz(j,k) - term(3,j)*rediv
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediv2
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediv2
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediv2
                           end do
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           redivk = rediv * redk
                           redivkv = rediv * redkv
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                              hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                              hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                              hessx(j,k) = hessx(j,k) - term(1,j)*redivk
                              hessy(j,k) = hessy(j,k) - term(2,j)*redivk
                              hessz(j,k) = hessz(j,k) - term(3,j)*redivk
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediv2
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediv2
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediv2
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redivkv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redivkv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redivkv
                           end do
                        end if
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emm3hb3  --  MM3 vdw & hbond energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emm3hb3" calculates the van der Waals and hydrogen bonding
c     interaction energy using MM3 exp-6 formula with directional
c     charge transfer hydrogen bonding; also partitions the energy
c     among the atoms
c
c     literature reference:
c
c     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in
c     the MM3 Force Field. I", Journal of Physical Organic Chemistry,
c     7, 591-609 (1994)
c
c
      subroutine emm3hb3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,ii,kk
      integer ia,ib,ic,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expcut2
      real*8 expterm,expmerge
      real*8 dot,cosine
      real*8 fterm,ideal
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical header,huge,proceed,iuse
c
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ne14 = 0
      ev = 0.0d0
      e14 = 0.0d0
      do i = 1, n
         aev(i) = 0.0d0
         ae14(i) = 0.0d0
      end do
      header = .true.
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find the van der Waals energy via double loop search
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
               if (rik2 .le. off2) then
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i) then
                     eps = eps / vdwscale
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = k
                     else
                        ia = k
                        ib = i12(1,k)
                        ic = i
                     end if
                     xia = x(ia)
                     yia = y(ia)
                     zia = z(ia)
                     xib = x(ib)
                     yib = y(ib)
                     zib = z(ib)
                     xic = x(ic)
                     yic = y(ic)
                     zic = z(ic)
                     xab = xia - xib
                     yab = yia - yib
                     zab = zia - zib
                     rab2 = xab*xab + yab*yab + zab*zab
                     rab = sqrt(rab2)
                     xcb = xic - xib
                     ycb = yic - yib
                     zcb = zic - zib
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     ideal = bl(bndlist(1,ia))
                     fterm = cosine * (rab/ideal)
                  end if
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - fterm*cterm*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (skip(k) .eq. -i) then
                     ne14 = ne14 + 1
                     e14 = e14 + e
                     ae14(i) = ae14(i) + 0.5d0*e
                     ae14(k) = ae14(k) + 0.5d0*e
                  else
                     nev = nev + 1
                     ev = ev + e
                     aev(i) = aev(i) + 0.5d0*e
                     aev(k) = aev(k) + 0.5d0*e
                  end if
c
c     print a warning if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',11x,'Atom Names',
     &                             18x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),
     &                                rv,sqrt(rik2),e
   20                format (' VDW-MM3  ',i5,'-',a3,1x,i5,'-',a3,
     &                          12x,2f10.4,f12.4)
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  kt = class(k)
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
                  if (rik2 .le. off2) then
                     fterm = 1.0d0
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (radhbnd(kt,it) .ne. 0.0d0) then
                        rv = radhbnd(kt,it)
                        eps = epshbnd(kt,it) / dielec
                        if (atomic(i) .eq. 1) then
                           ia = i
                           ib = i12(1,i)
                           ic = k
                        else
                           ia = k
                           ib = i12(1,k)
                           ic = i
                        end if
                        xia = x(ia)
                        yia = y(ia)
                        zia = z(ia)
                        xib = x(ib)
                        yib = y(ib)
                        zib = z(ib)
                        xic = x(ic)
                        yic = y(ic)
                        zic = z(ic)
                        xab = xia - xib
                        yab = yia - yib
                        zab = zia - zib
                        rab2 = xab*xab + yab*yab + zab*zab
                        rab = sqrt(rab2)
                        xcb = xic - xib
                        ycb = yic - yib
                        zcb = zic - zib
                        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                        xp = ycb*zab - zcb*yab
                        yp = zcb*xab - xcb*zab
                        zp = xcb*yab - ycb*xab
                        rp = sqrt(xp*xp + yp*yp + zp*zp)
                        dot = xab*xcb + yab*ycb + zab*zcb
                        cosine = dot / sqrt(rab2*rcb2)
                        ideal = bl(bndlist(1,ia))
                        fterm = cosine * (rab/ideal)
                     end if
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        expterm = aterm * exp(-bterm/p)
                        e = eps * (expterm - fterm*cterm*p6)
                     else
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik = sqrt(rik2)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                                + c2*rik2 + c1*rik + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy component
c
                     nev = nev + 1
                     if (i .eq. k) then
                        ev = ev + 0.5d0*e
                        aev(i) = aev(i) + 0.5d0*e
                     else
                        ev = ev + e
                        aev(i) = aev(i) + 0.5d0*e
                        aev(k) = aev(k) + 0.5d0*e
                     end if
c
c     print a warning if the energy of this interaction is large
c
                     huge = (e .gt. 10.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual van der Waals',
     &                                ' Interactions :',
     &                             //,' Type',11x,'Atom Names',
     &                                18x,'Minimum',4x,'Actual',
     &                                6x,'Energy',/)
                        end if
                        write (iout,40)  i,name(i),k,name(k),
     &                                   rv,sqrt(rik2),e
   40                   format (' VDW-MM3  ',i5,'-',a3,1x,i5,'-',a3,
     &                             '   (XTAL)   ',2f10.4,f12.4)
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emm3hb4  --  MM3 van der Waals and hbond energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emm3hb4" calculates the MM3 exp-6 van der Waals and hydrogen
c     bonding interaction energy using the method of lights to locate
c     neighboring atoms
c
c     literature reference:
c
c     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in
c     the MM3 Force Field. I", Journal of Physical Organic Chemistry,
c     7, 591-609 (1994)
c
c
      subroutine emm3hb4
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'iounit.i'
      include 'light.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,ii,kk
      integer iv,it,kt,ia,ib,ic
      integer kgy,kgz,start,stop
      integer kskip,skip(maxatm)
      integer kmap,kvmap,map(maxlight)
      real*8 e,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 expcut,expcut2
      real*8 expterm,expmerge
      real*8 dot,cosine
      real*8 fterm,ideal
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      real*8 xsort(maxlight),ysort(maxlight),zsort(maxlight)
      logical proceed,iuse,repeat
c
c
c     zero out the van der Waals energy contributions
c
      ev = 0.0d0
      e14 = 0.0d0
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do j = 1, nvdw
         i = ivdw(j)
         iv = ired(i)
         rdn = kred(i)
         xred(j) = rdn*(x(i)-x(iv)) + x(iv)
         yred(j) = rdn*(y(i)-y(iv)) + y(iv)
         zred(j) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nvdw,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            if (kk .le. nvdw) then
               kskip = skip(ivdw(kk))
               if (kskip .eq. i)  goto 20
            else
               kskip = 0
            end if
            kmap = ivdw(map(kk))
            kvmap = ired(kmap)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,kmap,0,0,0)
            if (proceed)  proceed = (iuse .or. use(kmap)
     &                                .or. use(kvmap))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(kmap)
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_image) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
               if (rik2 .le. off2) then
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (kskip .eq. -i) then
                     eps = eps / vdwscale
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = kmap
                        xab = xred(ii) - x(ib)
                        yab = yred(ii) - y(ib)
                        zab = zred(ii) - z(ib)
                        xcb = xab - xr
                        ycb = yab - yr
                        zcb = zab - zr
                     else
                        ia = kmap
                        ib = i12(1,kmap)
                        ic = i
                        xab = xred(kk) - x(ib)
                        yab = yred(kk) - y(ib)
                        zab = zred(kk) - z(ib)
                        xcb = xab + xr
                        ycb = yab + yr
                        zcb = zab + zr
                     end if
                     if (ired(ia) .ne. ia) then
                        xab = xab / kred(ia)
                        yab = yab / kred(ia)
                        zab = zab / kred(ia)
                     end if
                     rab2 = xab*xab + yab*yab + zab*zab
                     rab = sqrt(rab2)
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     ideal = bl(bndlist(1,ia))
                     fterm = cosine * (rab/ideal)
                  end if
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - fterm*cterm*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (kskip .eq. -i) then
                     e14 = e14 + e
                  else
                     ev = ev + e
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine emm3hb5  --  MM3 vdw & hbond energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "emm3hb5" calculates the MM3 exp-6 van der Waals and hydrogen
c     bonding interaction energy and its first derivatives using the
c     method of lights to locate neighboring atoms
c
c     literature reference:
c
c     J.-H. Lii and N. L. Allinger, "Directional Hydrogen Bonding in
c     the MM3 Force Field. I", Journal of Physical Organic Chemistry,
c     7, 591-609 (1994)
c
c
      subroutine emm3hb5
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bond.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'iounit.i'
      include 'light.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,ii,kk
      integer iv,it,kt
      integer ia,ib,ic
      integer kgy,kgz,start,stop
      integer kskip,skip(maxatm)
      integer kmap,kvmap,map(maxlight)
      real*8 e,rv,eps
      real*8 rdn,fgrp
      real*8 p,p2,p6,p12
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv,redk,redkv,rvterm
      real*8 dedx,dedy,dedz,de
      real*8 rik,rik2,rik3,rik4,rik5
      real*8 taper,dtaper
      real*8 expcut,expcut2
      real*8 expterm,expmerge
      real*8 dot,cosine,sine
      real*8 fterm,fcterm,term
      real*8 deddr,ideal,ratio
      real*8 deddt,terma,termc
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab2,rab,rcb2
      real*8 xp,yp,zp,rp
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      real*8 xsort(maxlight),ysort(maxlight),zsort(maxlight)
      logical proceed,iuse,repeat,use_hb
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      e14 = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
         de14(1,i) = 0.0d0
         de14(2,i) = 0.0d0
         de14(3,i) = 0.0d0
      end do
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do j = 1, nvdw
         i = ivdw(j)
         iv = ired(i)
         rdn = kred(i)
         xred(j) = rdn*(x(i)-x(iv)) + x(iv)
         yred(j) = rdn*(y(i)-y(iv)) + y(iv)
         zred(j) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nvdw,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = class(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            if (kk .le. nvdw) then
               kskip = skip(ivdw(kk))
               if (kskip .eq. i)  goto 20
            else
               kskip = 0
            end if
            kmap = ivdw(map(kk))
            kvmap = ired(kmap)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,kmap,0,0,0)
            if (proceed)  proceed = (iuse .or. use(kmap)
     &                                .or. use(kvmap))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(kmap)
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_image) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute energy and derivatives for this interaction
c
               if (rik2 .le. off2) then
                  use_hb = .false.
                  fterm = 1.0d0
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (kskip .eq. -i) then
                     eps = eps / vdwscale
                  else if (radhbnd(kt,it) .ne. 0.0d0) then
                     use_hb = .true.
                     rv = radhbnd(kt,it)
                     eps = epshbnd(kt,it) / dielec
                     if (atomic(i) .eq. 1) then
                        ia = i
                        ib = i12(1,i)
                        ic = kmap
                        xab = xred(ii) - x(ib)
                        yab = yred(ii) - y(ib)
                        zab = zred(ii) - z(ib)
                        xcb = xab - xr
                        ycb = yab - yr
                        zcb = zab - zr
                     else
                        ia = kmap
                        ib = i12(1,kmap)
                        ic = i
                        xab = xred(kk) - x(ib)
                        yab = yred(kk) - y(ib)
                        zab = zred(kk) - z(ib)
                        xcb = xab + xr
                        ycb = yab + yr
                        zcb = zab + zr
                     end if
                     if (ired(ia) .ne. ia) then
                        xab = xab / kred(ia)
                        yab = yab / kred(ia)
                        zab = zab / kred(ia)
                     end if
                     rab2 = xab*xab + yab*yab + zab*zab
                     rab = sqrt(rab2)
                     rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
                     xp = ycb*zab - zcb*yab
                     yp = zcb*xab - xcb*zab
                     zp = xcb*yab - ycb*xab
                     rp = sqrt(xp*xp + yp*yp + zp*zp)
                     dot = xab*xcb + yab*ycb + zab*zcb
                     cosine = dot / sqrt(rab2*rcb2)
                     sine = sqrt(abs(1.0d0-cosine**2))
                     ideal = bl(bndlist(1,ia))
                     ratio = rab / ideal
                     fterm = cosine * ratio
                     deddt = -sine * ratio
                     deddr = cosine / (rab*ideal)
                  end if
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bterm / rv
                     expterm = aterm * exp(-bterm/p)
                     fcterm = fterm * cterm * p6
                     e = eps * (expterm - fcterm)
                     de = eps * (rvterm*expterm+6.0d0*fcterm/rik)
                  else
                     use_hb = .false.
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     if (use_hb) then
                        deddt = deddt * fgrp
                        deddr = deddr * fgrp
                     end if
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  if (kskip .eq. -i) then
                     e14 = e14 + e
                     if (i .eq. iv) then
                        de14(1,i) = de14(1,i) + dedx
                        de14(2,i) = de14(2,i) + dedy
                        de14(3,i) = de14(3,i) + dedz
                     else
                        de14(1,i) = de14(1,i) + dedx*redi
                        de14(2,i) = de14(2,i) + dedy*redi
                        de14(3,i) = de14(3,i) + dedz*redi
                        de14(1,iv) = de14(1,iv) + dedx*rediv
                        de14(2,iv) = de14(2,iv) + dedy*rediv
                        de14(3,iv) = de14(3,iv) + dedz*rediv
                     end if
                     if (kmap .eq. kvmap) then
                        de14(1,kmap) = de14(1,kmap) - dedx
                        de14(2,kmap) = de14(2,kmap) - dedy
                        de14(3,kmap) = de14(3,kmap) - dedz
                     else
                        redk = kred(kmap)
                        redkv = 1.0d0 - redk
                        de14(1,kmap) = de14(1,kmap) - dedx*redk
                        de14(2,kmap) = de14(2,kmap) - dedy*redk
                        de14(3,kmap) = de14(3,kmap) - dedz*redk
                        de14(1,kvmap) = de14(1,kvmap) - dedx*redkv
                        de14(2,kvmap) = de14(2,kvmap) - dedy*redkv
                        de14(3,kvmap) = de14(3,kvmap) - dedz*redkv
                     end if
                  else
                     ev = ev + e
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                     end if
                     if (kmap .eq. kvmap) then
                        dev(1,kmap) = dev(1,kmap) - dedx
                        dev(2,kmap) = dev(2,kmap) - dedy
                        dev(3,kmap) = dev(3,kmap) - dedz
                     else
                        redk = kred(kmap)
                        redkv = 1.0d0 - redk
                        dev(1,kmap) = dev(1,kmap) - dedx*redk
                        dev(2,kmap) = dev(2,kmap) - dedy*redk
                        dev(3,kmap) = dev(3,kmap) - dedz*redk
                        dev(1,kvmap) = dev(1,kvmap) - dedx*redkv
                        dev(2,kvmap) = dev(2,kvmap) - dedy*redkv
                        dev(3,kvmap) = dev(3,kvmap) - dedz*redkv
                     end if
                  end if
c
c     find the chain rule terms for hydrogen bonding components
c
                  if (use_hb) then
                     term = eps * cterm * p6
                     deddt = deddt * term
                     deddr = deddr * term
                     if (rik2 .gt. cut2) then
                        deddt = deddt * taper
                        deddr = deddr * taper
                     end if
                     terma = deddt / (rab2*rp)
                     termc = -deddt / (rcb2*rp)
                     dedxia = terma * (yab*zp-zab*yp) - deddr*xab
                     dedyia = terma * (zab*xp-xab*zp) - deddr*yab
                     dedzia = terma * (xab*yp-yab*xp) - deddr*zab
                     dedxic = termc * (ycb*zp-zcb*yp)
                     dedyic = termc * (zcb*xp-xcb*zp)
                     dedzic = termc * (xcb*yp-ycb*xp)
                     dedxib = -dedxia - dedxic
                     dedyib = -dedyia - dedyic
                     dedzib = -dedzia - dedzic
c
c     increment the derivatives for directional hydrogen bonding
c
                     dev(1,ia) = dev(1,ia) + dedxia
                     dev(2,ia) = dev(2,ia) + dedyia
                     dev(3,ia) = dev(3,ia) + dedzia
                     dev(1,ib) = dev(1,ib) + dedxib
                     dev(2,ib) = dev(2,ib) + dedyib
                     dev(3,ib) = dev(3,ib) + dedzib
                     dev(1,ic) = dev(1,ic) + dedxic
                     dev(2,ic) = dev(2,ic) + dedyic
                     dev(3,ic) = dev(3,ic) + dedzic
                  end if
c
c     increment the total intermolecular energy
c
                  if (kk .le. nvdw) then
                     if (molcule(i) .ne. molcule(kmap)) then
                        einter = einter + e
                     end if
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx + xr*dedx
                     viry = viry + yr*dedy
                     virz = virz + zr*dedz
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine empole  --  multipole & polarization energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "empole" calculates the electrostatic energy due to
c     atomic multipole interactions and dipole polarizability
c
c
      subroutine empole
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,m
      integer ii,iz,ix
      integer kk,kz,kx
      integer skip(maxatm)
      real*8 eik,ei,ek,a(3,3)
      real*8 shift,taper,trans
      real*8 f,fik,fgrp
      real*8 xr,yr,zr,r
      real*8 r2,r3,r4,r5,r6,r7
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      logical proceed,iuse,kuse
c
c
c     zero out the atomic multipole interaction energy
c
      em = 0.0d0
      ep = 0.0d0
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     rotate the multipole components into the global frame
c
      do i = 1, npole
cjrs
         call rotmatt (i,a)
cjrs
         call rotpole (i,a)
      end do
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     calculate the multipole interaction energy term
c
      do ii = 1, npole-1
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * chg12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * chg13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * chg14use
         end do
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, maxpole
                     rpk(j) = rpole(j,kk)
                  end do
                  do j = 1, 3
                     indk(j) = uind(j,kk)
                  end do
                  r = sqrt(r2)
                  call empik (ii,kk,xr,yr,zr,r,rpi,rpk,
     &                          indi,indk,eik,ei,ek)
                  fik = f
                  if (skip(k) .eq. -i)  fik = fik / chgscale
                  eik = fik * eik
                  ei = fik * ei
                  ek = fik * ek
c
c     use shifted energy switching if near the cutoff distance
c
                  fik = fik * rpi(1) * rpk(1)
                  shift = fik / (0.5d0*(off+cut))
                  eik = eik - shift
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     eik = eik * taper + trans
                     ei = ei * taper
                     ek = ek * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     eik = eik * fgrp
                     ei = ei * fgrp
                     ek = ek * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  em = em + eik
                  ep = ep + ei + ek
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do m = 1, ncell
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     do j = 1, 3
                        indk(j) = uind(j,kk)
                     end do
                     r = sqrt(r2)
                     call empik (ii,kk,xr,yr,zr,r,rpi,rpk,
     &                             indi,indk,eik,ei,ek)
                     fik = f
                     eik = fik * eik
                     ei = fik * ei
                     ek = fik * ek
c
c     use shifted energy switching if near the cutoff distance
c
                     fik = fik * rpi(1) * rpk(1)
                     shift = fik / (0.5d0*(off+cut))
                     eik = eik - shift
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                        eik = eik * taper + trans
                        ei = ei * taper
                        ek = ek * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        eik = eik * fgrp
                        ei = ei * fgrp
                        ek = ek * fgrp
                     end if
c
c     increment the overall multipole and polarization energies
c
                     if (i .eq. k) then 
                        eik = 0.5d0 * eik
                        ei = 0.5d0 * ei
                        ek = 0.5d0 * ek
                     end if
                     em = em + eik
                     ep = ep + ei + ek
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empik  --  multipole & polarization pair energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empik" computes the permanent multipole and induced dipole
c     energies between a specified pair of atomic multipole sites
c
c
      subroutine empik (ii,kk,xr,yr,zr,r,rpi,rpk,indi,indk,eik,ei,ek)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      integer ii,kk
      real*8 eik,ei,ek
      real*8 xr,yr,zr,r,damp
      real*8 rr1,rr2,rr3,rr5,rr7,rr9
      real*8 xx,yy,xyz,xry,yrx,zrx,zry
      real*8 xr5,yr5,xr7,yr7
      real*8 xr9,xyr9,yr9,xr9x7,yr9y7
      real*8 rr53,xr73,yr73
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t11
      real*8 t12,t13,t14,t15,t17,t18,t21,t22
      real*8 t23,t24,t25,t27,t28,t31,t32
      real*8 i0,i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 k0,k1,k2,k3,k4,k5,k6,k7,k8,k9
      real*8 u1,u2,u3,v1,v2,v3
      real*8 i3k3,i4i9,i7i9,k4k9,k7k9,i3k6
      real*8 i3k8,i6k3,i8k3,ik49,ik59,ik69
      real*8 ik79,ik89,ik68,i6k6,i8k8,i9k9
      real*8 i3v3,i6v3,i8v3,k3u3,k6u3,k8u3
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      logical iquad,kquad
c
c
c     zero out multipole and polarization energy components
c
      eik = 0.0d0
      ei = 0.0d0
      ek = 0.0d0
c
c     check for presence of quadrupole components at either site
c
      iquad = (polsiz(ii) .ge. 13)
      kquad = (polsiz(kk) .ge. 13)
c
c     set permanent and induced multipoles for first site
c
      i0 = rpi(1)
      i1 = rpi(2)
      i2 = rpi(3)
      i3 = rpi(4)
      i4 = rpi(5)
      i5 = rpi(6)
      i6 = rpi(7)
      i7 = rpi(9)
      i8 = rpi(10)
      i9 = rpi(13)
      u1 = indi(1)
      u2 = indi(2)
      u3 = indi(3)
c
c     set permanent and induced multipoles for second site
c
      k0 = rpk(1)
      k1 = rpk(2)
      k2 = rpk(3)
      k3 = rpk(4)
      k4 = rpk(5)
      k5 = rpk(6)
      k6 = rpk(7)
      k7 = rpk(9)
      k8 = rpk(10)
      k9 = rpk(13)
      v1 = indk(1)
      v2 = indk(2)
      v3 = indk(3)
c
c     compute the zeroth order T2 matrix element
c
      rr1 = 1.0d0 / r
      t1 = rr1
c
c     compute the first order T2 matrix elements
c
      rr2 = rr1 * rr1
      rr3 = rr1 * rr2
      t2 = -xr * rr3
      t3 = -yr * rr3
      t4 = -zr * rr3
c
c     compute the second order T2 matrix elements
c
      rr5 = 3.0d0 * rr3 * rr2
      xr5 = xr * rr5
      yr5 = yr * rr5
      t5 = -rr3 + xr5*xr
      t6 = xr5 * yr
      t7 = xr5 * zr
      t8 = -rr3 + yr5*yr
      t9 = yr5 * zr
c
c     compute the third order T2 matrix elements
c
      if (iquad .or. kquad) then
         rr7 = 5.0d0 * rr5 * rr2
         xx = xr * xr
         yy = yr * yr
         xyz = xr * yr * zr
         xr7 = xx * rr7
         yr7 = yy * rr7
         rr53 = 3.0d0 * rr5
         t11 = xr * (rr53-xr7)
         t12 = yr * (rr5-xr7)
         t13 = zr * (rr5-xr7)
         t14 = xr * (rr5-yr7)
         t15 = -xyz * rr7
         t17 = yr * (rr53-yr7)
         t18 = zr * (rr5-yr7)
c
c     compute the fourth order T2 matrix elements
c
         rr9 = 7.0d0 * rr7 * rr2
         if (xr .eq. 0.0d0) then
            yrx = 0.0d0
            zrx = 0.0d0
         else
            yrx = yr / xr
            zrx = zr / xr
         end if
         if (yr .eq. 0.0d0) then
            xry = 0.0d0
            zry = 0.0d0
         else
            xry = xr / yr
            zry = zr / yr
         end if
         xr9 = xx * xx * rr9
         xyr9 = xx * yy * rr9
         yr9 = yy * yy * rr9
         xr73 = 3.0d0 * xr7
         yr73 = 3.0d0 * yr7
         xr9x7 = xr9 - xr73
         yr9y7 = yr9 - yr73
         t21 = xr9x7 - xr73 + rr53
         t22 = yrx * xr9x7
         t23 = zrx * xr9x7
         t24 = xyr9 - xr7 - yr7 + rr5
         t25 = zry * (xyr9-yr7)
         t27 = xry * yr9y7
         t28 = zrx * (xyr9-xr7)
         t31 = yr9y7 - yr73 + rr53
         t32 = zry * yr9y7
      end if
c
c     get the M-M, M-D and D-D parts of the multipole energy
c
      i3k3 = i3 * k3
      eik = i0*k0*t1 + (i0*k1-i1*k0)*t2 + (i0*k2-i2*k0)*t3
     &         + (i0*k3-i3*k0)*t4 + (i3k3-i1*k1)*t5 - (i1*k2+i2*k1)*t6
     &         - (i1*k3+i3*k1)*t7 + (i3k3-i2*k2)*t8 - (i2*k3+i3*k2)*t9
c
c     get the M-indD and D-indD parts of the polarization energy
c
      i3v3 = i3 * v3
      k3u3 = k3 * u3
      ei = -k0*(u1*t2+u2*t3+u3*t4) + (k3u3-k1*u1)*t5 - (k1*u2+k2*u1)*t6
     &        - (k3*u1+k1*u3)*t7 + (k3u3-k2*u2)*t8 - (k2*u3+k3*u2)*t9
      ek = i0*(v1*t2+v2*t3+v3*t4) + (i3v3-i1*v1)*t5 - (i1*v2+i2*v1)*t6
     &        - (i3*v1+i1*v3)*t7 + (i3v3-i2*v2)*t8 - (i2*v3+i3*v2)*t9
c
c     get the M-Q, D-Q and Q-Dind energies, if necessary
c
      if (kquad) then
         k4k9 = k4 - k9
         k7k9 = k7 - k9
         i3k6 = 2.0d0 * i3 * k6
         i3k8 = 2.0d0 * i3 * k8
         eik = eik + i0*k4k9*t5 + 2.0d0*i0*k5*t6 + 2.0d0*i0*k6*t7
     &             + i0*k7k9*t8 + 2.0d0*i0*k8*t9 + (i3k6-i1*k4k9)*t11
     &             + (i3k8-i2*k4k9-2.0d0*i1*k5)*t12
     &             - (2.0d0*i1*k6+i3*k4k9)*t13
     &             + (i3k6-i1*k7k9-2.0d0*i2*k5)*t14
     &             - 2.0d0*(i1*k8+i2*k6+i3*k5)*t15
     &             + (i3k8-i2*k7k9)*t17
     &             - (2.0d0*i2*k8+i3*k7k9)*t18
         k6u3 = k6 * u3
         k8u3 = k8 * u3
         ei = ei + (-k4k9*u1+2.0d0*k6u3)*t11
     &           + (-k4k9*u2+2.0d0*(k8u3-k5*u1))*t12
     &           + (-k4k9*u3-2.0d0*k6*u1)*t13
     &           + (-k7k9*u1+2.0d0*(k6u3-k5*u2))*t14
     &           - 2.0d0*(k5*u3+k6*u2+k8*u1)*t15
     &           + (-k7k9*u2+2.0d0*k8u3)*t17
     &           + (-k7k9*u3-2.0d0*k8*u2)*t18
      end if
      if (iquad) then
         i4i9 = i4 - i9
         i7i9 = i7 - i9
         i6k3 = -2.0d0 * i6 * k3
         i8k3 = -2.0d0 * i8 * k3
         eik = eik + i4i9*k0*t5 + 2.0d0*i5*k0*t6 + 2.0d0*i6*k0*t7
     &             + i7i9*k0*t8 + 2.0d0*i8*k0*t9 + (i4i9*k1+i6k3)*t11
     &             + (i8k3+i4i9*k2+2.0d0*i5*k1)*t12
     &             + (2.0d0*i6*k1+i4i9*k3)*t13
     &             + (i6k3+i7i9*k1+2.0d0*i5*k2)*t14
     &             + 2.0d0*(i8*k1+i6*k2+i5*k3)*t15
     &             + (i8k3+i7i9*k2)*t17
     &             + (2.0d0*i8*k2+i7i9*k3)*t18
         i6v3 = i6 * v3
         i8v3 = i8 * v3
         ek = ek + (i4i9*v1-2.0d0*i6v3)*t11
     &           + (i4i9*v2+2.0d0*(i5*v1-i8v3))*t12
     &           + (i4i9*v3+2.0d0*i6*v1)*t13
     &           + (i7i9*v1+2.0d0*(i5*v2-i6v3))*t14
     &           + 2.0d0*(i5*v3+i6*v2+i8*v1)*t15
     &           + (i7i9*v2-2.0d0*i8v3)*t17
     &           + (i7i9*v3+2.0d0*i8*v2)*t18
      end if
c
c     get the Q-Q part of the multipole interaction energy
c
      if (iquad .and. kquad) then
         ik49 = i4*k9 + i9*k4
         ik59 = i5*k9 + i9*k5
         ik69 = i6*k9 + i9*k6
         ik79 = i7*k9 + i9*k7
         ik89 = i8*k9 + i9*k8
         ik68 = 2.0d0 * (i6*k8 + i8*k6)
         i6k6 = i6 * k6
         i8k8 = i8 * k8
         i9k9 = i9 * k9
         eik = eik + (i4*k4-ik49+i9k9-4.0d0*i6k6)*t21
     &             + 2.0d0*(i5*k4+i4*k5-ik68-ik59)*t22
     &             + 2.0d0*(i4*k6+i6*k4-ik69)*t23
     &             + (i4*k7+i7*k4-ik49-ik79+2.0d0*i9k9
     &                  +4.0d0*(i5*k5-i6k6-i8k8))*t24
     &             + 2.0d0*(i4*k8+i8*k4+2.0d0*(i5*k6+i6*k5)-ik89)*t25
     &             + 2.0d0*(i5*k7+i7*k5-ik68-ik59)*t27
     &             + 2.0d0*(i6*k7+i7*k6+2.0d0*(i5*k8+i8*k5)-ik69)*t28
     &             + (i7*k7-ik79+i9k9-4.0d0*i8k8)*t31
     &             + 2.0d0*(i7*k8+i8*k7-ik89)*t32
      end if
c
c     apply exponential damping to the polarization term
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            ei = ei * damp
            ek = ek * damp
         end if
      end if
c
c     final polarization energy is half of value computed above
c
      ei = 0.5d0 * ei
      ek = 0.5d0 * ek
      return
      end
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1  --  mpole/polar energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,m,jj
      integer ii,iz,ix
      integer kk,kz,kx
      integer skip(maxatm)
      real*8 eik,ei,ek,de
      real*8 f,fik,fgrp
      real*8 shift,taper,dtaper
      real*8 trans,dtrans
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xiz,yiz,ziz,xix,yix,zix
      real*8 xkz,ykz,zkz,xkx,ykx,zkx
      real*8 r,r2,r3,r4,r5,r6,r7
      real*8 a(3,3),d(3,3,3,3)
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 dm(3),dp(3),utu
      real*8 dmi(3,3),dmk(3,3)
      real*8 dpi(3,3),dpk(3,3)
      logical proceed,iuse,kuse
c
c
c     zero out multipole and polarization energy and derivatives
c
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
         end do
      end do
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     rotate multipole components and get rotation derivatives
c
      do i = 1, npole
cjrs
         call rotmatt (i,a)
cjrs
         call rotpole (i,a)
         call drotmat (i,d)
         do j = 1, 3
            do k = 1, 3
               call drotpole (i,a,d,j,k)
            end do
         end do
      end do
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     compute the multipole interaction energy and derivatives
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * chg12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * chg13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * chg14use
         end do
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, maxpole
                     rpk(j) = rpole(j,kk)
                  end do
                  do j = 1, 3
                     indk(j) = uind(j,kk)
                  end do
                  r = sqrt(r2)
                  call empik1 (ii,kk,xr,yr,zr,r,r2,rpi,
     &                         rpk,indi,indk,eik,ei,ek,
     &                         dm,dmi,dmk,dp,dpi,dpk,utu)
                  fik = f
                  if (skip(k) .eq. -i)  fik = fik / chgscale
                  eik = fik * eik
                  ei = fik * ei
                  ek = fik * ek
                  utu = fik * utu
                  do j = 1, 3
                     dm(j) = fik * dm(j)
                     dp(j) = fik * dp(j)
                     do jj = 1, 3
                        dmi(jj,j) = fik * dmi(jj,j)
                        dmk(jj,j) = fik * dmk(jj,j)
                        dpi(jj,j) = fik * dpi(jj,j)
                        dpk(jj,j) = fik * dpk(jj,j)
                     end do
                  end do
c
c     use shifted energy switching if near the cutoff distance
c
                  fik = fik * rpi(1) * rpk(1)
                  shift = fik / (0.5d0*(off+cut))
                  eik = eik - shift
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                               + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                             + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                     de = eik * dtaper + dtrans
                     eik = eik * taper + trans
                     dm(1) = dm(1)*taper + de*(xr/r)
                     dm(2) = dm(2)*taper + de*(yr/r)
                     dm(3) = dm(3)*taper + de*(zr/r)
                     de = (2.0d0*(ei+ek)+utu) * dtaper
                     ei = ei * taper
                     ek = ek * taper
                     dp(1) = dp(1)*taper + de*(xr/r)
                     dp(2) = dp(2)*taper + de*(yr/r)
                     dp(3) = dp(3)*taper + de*(zr/r)
                     do j = 1, 3
                        do jj = 1, 3
                           dmi(jj,j) = dmi(jj,j) * taper
                           dmk(jj,j) = dmk(jj,j) * taper
                           dpi(jj,j) = dpi(jj,j) * taper
                           dpk(jj,j) = dpk(jj,j) * taper
                        end do
                     end do
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     eik = eik * fgrp
                     ei = ei * fgrp
                     ek = ek * fgrp
                     do j = 1, 3
                        dm(j) = dm(j) * fgrp
                        dp(j) = dp(j) * fgrp
                        do jj = 1, 3
                           dmi(jj,j) = dmi(jj,j) * fgrp
                           dmk(jj,j) = dmk(jj,j) * fgrp
                           dpi(jj,j) = dpi(jj,j) * fgrp
                           dpk(jj,j) = dpk(jj,j) * fgrp
                        end do
                     end do
                  end if
c
c     increment the multipole energy and derivative expressions
c
                  em = em + eik
                  dem(1,i) = dem(1,i) - dm(1) + dmi(1,1)
                  dem(2,i) = dem(2,i) - dm(2) + dmi(1,2)
                  dem(3,i) = dem(3,i) - dm(3) + dmi(1,3)
                  dem(1,iz) = dem(1,iz) + dmi(2,1)
                  dem(2,iz) = dem(2,iz) + dmi(2,2)
                  dem(3,iz) = dem(3,iz) + dmi(2,3)
                  dem(1,ix) = dem(1,ix) + dmi(3,1)
                  dem(2,ix) = dem(2,ix) + dmi(3,2)
                  dem(3,ix) = dem(3,ix) + dmi(3,3)
                  dem(1,k) = dem(1,k) + dm(1) + dmk(1,1)
                  dem(2,k) = dem(2,k) + dm(2) + dmk(1,2)
                  dem(3,k) = dem(3,k) + dm(3) + dmk(1,3)
                  dem(1,kz) = dem(1,kz) + dmk(2,1)
                  dem(2,kz) = dem(2,kz) + dmk(2,2)
                  dem(3,kz) = dem(3,kz) + dmk(2,3)
                  dem(1,kx) = dem(1,kx) + dmk(3,1)
                  dem(2,kx) = dem(2,kx) + dmk(3,2)
                  dem(3,kx) = dem(3,kx) + dmk(3,3)
c
c     increment the polarization energy and derivative expressions
c
                  ep = ep + ei + ek
                  dep(1,i) = dep(1,i) - dp(1) + dpi(1,1)
                  dep(2,i) = dep(2,i) - dp(2) + dpi(1,2)
                  dep(3,i) = dep(3,i) - dp(3) + dpi(1,3)
                  dep(1,iz) = dep(1,iz) + dpi(2,1)
                  dep(2,iz) = dep(2,iz) + dpi(2,2)
                  dep(3,iz) = dep(3,iz) + dpi(2,3)
                  dep(1,ix) = dep(1,ix) + dpi(3,1)
                  dep(2,ix) = dep(2,ix) + dpi(3,2)
                  dep(3,ix) = dep(3,ix) + dpi(3,3)
                  dep(1,k) = dep(1,k) + dp(1) + dpk(1,1)
                  dep(2,k) = dep(2,k) + dp(2) + dpk(1,2)
                  dep(3,k) = dep(3,k) + dp(3) + dpk(1,3)
                  dep(1,kz) = dep(1,kz) + dpk(2,1)
                  dep(2,kz) = dep(2,kz) + dpk(2,2)
                  dep(3,kz) = dep(3,kz) + dpk(2,3)
                  dep(1,kx) = dep(1,kx) + dpk(3,1)
                  dep(2,kx) = dep(2,kx) + dpk(3,2)
                  dep(3,kx) = dep(3,kx) + dpk(3,3)
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + eik + ei + ek
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     xiz = x(iz) - x(i)
                     yiz = y(iz) - y(i)
                     ziz = z(iz) - z(i)
                     xix = x(ix) - x(i)
                     yix = y(ix) - y(i)
                     zix = z(ix) - z(i)
                     xkz = x(kz) - x(k)
                     ykz = y(kz) - y(k)
                     zkz = z(kz) - z(k)
                     xkx = x(kx) - x(k)
                     ykx = y(kx) - y(k)
                     zkx = z(kx) - z(k)
                     virx = virx + xr*(dm(1)+dp(1))
     &                           + xiz*(dmi(2,1)+dpi(2,1))
     &                           + xix*(dmi(3,1)+dpi(3,1))
     &                           + xkz*(dmk(2,1)+dpk(2,1))
     &                           + xkx*(dmk(3,1)+dpk(3,1))
                     viry = viry + yr*(dm(2)+dp(2))
     &                           + yiz*(dmi(2,2)+dpi(2,2))
     &                           + yix*(dmi(3,2)+dpi(3,2))
     &                           + ykz*(dmk(2,2)+dpk(2,2))
     &                           + ykx*(dmk(3,2)+dpk(3,2))
                     virz = virz + zr*(dm(3)+dp(3))
     &                           + ziz*(dmi(2,3)+dpi(2,3))
     &                           + zix*(dmi(3,3)+dpi(3,3))
     &                           + zkz*(dmk(2,3)+dpk(2,3))
     &                           + zkx*(dmk(3,3)+dpk(3,3))
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do m = 1, ncell
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     do j = 1, 3
                        indk(j) = uind(j,kk)
                     end do
                     r = sqrt(r2)
                     call empik1 (ii,kk,xr,yr,zr,r,r2,rpi,
     &                            rpk,indi,indk,eik,ei,ek,
     &                            dm,dmi,dmk,dp,dpi,dpk,utu)
                     fik = f
                     eik = fik * eik
                     ei = fik * ei
                     ek = fik * ek
                     utu = fik * utu
                     do j = 1, 3
                        dm(j) = fik * dm(j)
                        dp(j) = fik * dp(j)
                        do jj = 1, 3
                           dmi(jj,j) = fik * dmi(jj,j)
                           dmk(jj,j) = fik * dmk(jj,j)
                           dpi(jj,j) = fik * dpi(jj,j)
                           dpk(jj,j) = fik * dpk(jj,j)
                        end do
                     end do
c
c     use shifted energy switching if near the cutoff distance
c
                     fik = fik * rpi(1) * rpk(1)
                     shift = fik / (0.5d0*(off+cut))
                     eik = eik - shift
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                        dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                                  + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                                + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                        de = eik * dtaper + dtrans
                        eik = eik * taper + trans
                        dm(1) = dm(1)*taper + de*(xr/r)
                        dm(2) = dm(2)*taper + de*(yr/r)
                        dm(3) = dm(3)*taper + de*(zr/r)
                        de = (2.0d0*(ei+ek)+utu) * dtaper
                        ei = ei * taper
                        ek = ek * taper
                        dp(1) = dp(1)*taper + de*(xr/r)
                        dp(2) = dp(2)*taper + de*(yr/r)
                        dp(3) = dp(3)*taper + de*(zr/r)
                        do j = 1, 3
                           do jj = 1, 3
                              dmi(jj,j) = dmi(jj,j) * taper
                              dmk(jj,j) = dmk(jj,j) * taper
                              dpi(jj,j) = dpi(jj,j) * taper
                              dpk(jj,j) = dpk(jj,j) * taper
                           end do
                        end do
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        eik = eik * fgrp
                        ei = ei * fgrp
                        ek = ek * fgrp
                        do j = 1, 3
                           dm(j) = dm(j) * fgrp
                           dp(j) = dp(j) * fgrp
                           do jj = 1, 3
                              dmi(jj,j) = dmi(jj,j) * fgrp
                              dmk(jj,j) = dmk(jj,j) * fgrp
                              dpi(jj,j) = dpi(jj,j) * fgrp
                              dpk(jj,j) = dpk(jj,j) * fgrp
                           end do
                        end do
                     end if
c
c     increment the multipole energy and derivative expressions
c
                     if (i .eq. k)  eik = 0.5d0 * eik
                     em = em + eik
                     dem(1,i) = dem(1,i) - dm(1) + dmi(1,1)
                     dem(2,i) = dem(2,i) - dm(2) + dmi(1,2)
                     dem(3,i) = dem(3,i) - dm(3) + dmi(1,3)
                     dem(1,iz) = dem(1,iz) + dmi(2,1)
                     dem(2,iz) = dem(2,iz) + dmi(2,2)
                     dem(3,iz) = dem(3,iz) + dmi(2,3)
                     dem(1,ix) = dem(1,ix) + dmi(3,1)
                     dem(2,ix) = dem(2,ix) + dmi(3,2)
                     dem(3,ix) = dem(3,ix) + dmi(3,3)
                     if (i .ne. k) then
                        dem(1,k) = dem(1,k) + dm(1) + dmk(1,1)
                        dem(2,k) = dem(2,k) + dm(2) + dmk(1,2)
                        dem(3,k) = dem(3,k) + dm(3) + dmk(1,3)
                        dem(1,kz) = dem(1,kz) + dmk(2,1)
                        dem(2,kz) = dem(2,kz) + dmk(2,2)
                        dem(3,kz) = dem(3,kz) + dmk(2,3)
                        dem(1,kx) = dem(1,kx) + dmk(3,1)
                        dem(2,kx) = dem(2,kx) + dmk(3,2)
                        dem(3,kx) = dem(3,kx) + dmk(3,3)
                     end if
c
c     increment the polarization energy and derivative expressions
c
                     if ( i .eq. k)  ek = 0.0d0
                     ep = ep + ei + ek
                     dep(1,i) = dep(1,i) - dp(1) + dpi(1,1)
                     dep(2,i) = dep(2,i) - dp(2) + dpi(1,2)
                     dep(3,i) = dep(3,i) - dp(3) + dpi(1,3)
                     dep(1,iz) = dep(1,iz) + dpi(2,1)
                     dep(2,iz) = dep(2,iz) + dpi(2,2)
                     dep(3,iz) = dep(3,iz) + dpi(2,3)
                     dep(1,ix) = dep(1,ix) + dpi(3,1)
                     dep(2,ix) = dep(2,ix) + dpi(3,2)
                     dep(3,ix) = dep(3,ix) + dpi(3,3)
                     if (i .ne. k) then
                        dep(1,k) = dep(1,k) + dp(1) + dpk(1,1)
                        dep(2,k) = dep(2,k) + dp(2) + dpk(1,2)
                        dep(3,k) = dep(3,k) + dp(3) + dpk(1,3)
                        dep(1,kz) = dep(1,kz) + dpk(2,1)
                        dep(2,kz) = dep(2,kz) + dpk(2,2)
                        dep(3,kz) = dep(3,kz) + dpk(2,3)
                        dep(1,kx) = dep(1,kx) + dpk(3,1)
                        dep(2,kx) = dep(2,kx) + dpk(3,2)
                        dep(3,kx) = dep(3,kx) + dpk(3,3)
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + eik + ei + ek
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        xiz = x(iz) - x(i)
                        yiz = y(iz) - y(i)
                        ziz = z(iz) - z(i)
                        xix = x(ix) - x(i)
                        yix = y(ix) - y(i)
                        zix = z(ix) - z(i)
                        xkz = x(kz) - x(k)
                        ykz = y(kz) - y(k)
                        zkz = z(kz) - z(k)
                        xkx = x(kx) - x(k)
                        ykx = y(kx) - y(k)
                        zkx = z(kx) - z(k)
                        virx = virx + xr*(dm(1)+dp(1))
     &                              + xiz*(dmi(2,1)+dpi(2,1))
     &                              + xix*(dmi(3,1)+dpi(3,1))
     &                              + xkz*(dmk(2,1)+dpk(2,1))
     &                              + xkx*(dmk(3,1)+dpk(3,1))
                        viry = viry + yr*(dm(2)+dp(2))
     &                              + yiz*(dmi(2,2)+dpi(2,2))
     &                              + yix*(dmi(3,2)+dpi(3,2))
     &                              + ykz*(dmk(2,2)+dpk(2,2))
     &                              + ykx*(dmk(3,2)+dpk(3,2))
                        virz = virz + zr*(dm(3)+dp(3))
     &                              + ziz*(dmi(2,3)+dpi(2,3))
     &                              + zix*(dmi(3,3)+dpi(3,3))
     &                              + zkz*(dmk(2,3)+dpk(2,3))
     &                              + zkx*(dmk(3,3)+dpk(3,3))
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empik1  --  mpole & polarization pair gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empik1" computes the permanent multipole and induced dipole
c     energies and derivatives between a pair of multipole sites
c
c
      subroutine empik1 (ii,kk,xr,yr,zr,r,r2,rpi,rpk,indi,indk,
     &                    eik,ei,ek,dm,dmi,dmk,dp,dpi,dpk,utu)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      integer i,j,ii,kk
      real*8 eik,ei,ek
      real*8 damp,ddamp,de,term
      real*8 xr,yr,zr,r,r2,r4
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 dm(3),dp(3),utu
      real*8 dmi(3,3),dmk(3,3)
      real*8 dpi(3,3),dpk(3,3)
      real*8 rr1,rr2,rr3,rr5,rr7,rr9,rr11
      real*8 xx,yy,xyz,xry,yrx,zrx,zry
      real*8 xr5,yr5,xr7,yr7
      real*8 xr9,xyr9,yr9,xr9x7,yr9y7
      real*8 rr53,xr73,yr73
      real*8 x2,y2,z2,x4,y4,x2y2,x2r2,y2r2
      real*8 xrr11,yrr11,zrr11,xyzrr11
      real*8 r445,r4225,nnn1,nnn2
      real*8 n945x4,n945y4,n945x2y2
      real*8 i3k3,i6k6,i8k8,i9k9
      real*8 ik49,ik59,ik69,ik79,ik89,ik68
      real*8 k4k9,k7k9,i4i9,i7i9
      real*8 k3u3,i3k6,i3k8,k6u3,k8u3
      real*8 i3v3,i6k3,i8k3,i6v3,i8v3
      real*8 i0,i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 k0,k1,k2,k3,k4,k5,k6,k7,k8,k9
      real*8 u1,u2,u3,v1,v2,v3
      real*8 di1,di2,di3,di4,di5,di6,di7,di8,di9
      real*8 dk1,dk2,dk3,dk4,dk5,dk6,dk7,dk8,dk9
      real*8 w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13
      real*8 t14,t15,t16,t17,t18,t19,t21,t22,t23,t24,t25
      real*8 t26,t27,t28,t29,t31,t32,t33,t36,t37,t38,t39
      real*8 t40,t41,t42,t43,t44,t46,t47,t48,t51,t52,t53
      logical iquad,kquad
c
c
c     zero out the energy and first derivative components
c
      eik = 0.0d0
      ei = 0.0d0
      ek = 0.0d0
      utu = 0.0d0
      do i = 1, 3
         dm(i) = 0.0d0
         dp(i) = 0.0d0
         do j = 1, 3
            dmi(j,i) = 0.0d0
            dmk(j,i) = 0.0d0
            dpi(j,i) = 0.0d0
            dpk(j,i) = 0.0d0
         end do
      end do
c
c     check for presence of quadrupole components at either site
c
      iquad = (polsiz(ii) .ge. 13)
      kquad = (polsiz(kk) .ge. 13)
c
c     set permanent and induced multipoles for first site
c
      i0 = rpi(1)
      i1 = rpi(2)
      i2 = rpi(3)
      i3 = rpi(4)
      i4 = rpi(5)
      i5 = rpi(6)
      i6 = rpi(7)
      i7 = rpi(9)
      i8 = rpi(10)
      i9 = rpi(13)
      u1 = indi(1)
      u2 = indi(2)
      u3 = indi(3)
c
c     set permanent and induced multipoles for second site
c
      k0 = rpk(1)
      k1 = rpk(2)
      k2 = rpk(3)
      k3 = rpk(4)
      k4 = rpk(5)
      k5 = rpk(6)
      k6 = rpk(7)
      k7 = rpk(9)
      k8 = rpk(10)
      k9 = rpk(13)
      v1 = indk(1)
      v2 = indk(2)
      v3 = indk(3)
c
c     compute the zeroth order T2 matrix element
c
      rr1 = 1.0d0 / r
      t1 = rr1
c
c     compute the first order T2 matrix elements
c
      rr2 = rr1 * rr1
      rr3 = rr1 * rr2
      t2 = -xr * rr3
      t3 = -yr * rr3
      t4 = -zr * rr3
c
c     compute the second order T2 matrix elements
c
      rr5 = 3.0d0 * rr3 * rr2
      xr5 = xr * rr5
      yr5 = yr * rr5
      t5 = -rr3 + xr5*xr
      t6 = xr5 * yr
      t7 = xr5 * zr
      t8 = -rr3 + yr5*yr
      t9 = yr5 * zr
      t10 = -t5 - t8
c
c     compute the third order T2 matrix elements
c
      rr7 = 5.0d0 * rr5 * rr2
      xx = xr * xr
      yy = yr * yr
      xyz = xr * yr * zr
      xr7 = xx * rr7
      yr7 = yy * rr7
      rr53 = 3.0d0 * rr5
      t11 = xr * (rr53-xr7)
      t12 = yr * (rr5-xr7)
      t13 = zr * (rr5-xr7)
      t14 = xr * (rr5-yr7)
      t15 = -xyz * rr7
      t16 = -t11 - t14
      t17 = yr * (rr53-yr7)
      t18 = zr * (rr5-yr7)
      t19 = -t12 - t17
c
c     compute the fourth order T2 matrix elements
c
      if (iquad .or. kquad) then
         rr9 = 7.0d0 * rr7 * rr2
         if (xr .eq. 0.0d0) then
            yrx = 0.0d0
            zrx = 0.0d0
         else
            yrx = yr / xr
            zrx = zr / xr
         end if
         if (yr .eq. 0.0d0) then
            xry = 0.0d0
            zry = 0.0d0
         else
            xry = xr / yr
            zry = zr / yr
         end if
         xr9 = xx * xx * rr9
         xyr9 = xx * yy * rr9
         yr9 = yy * yy * rr9
         xr73 = 3.0d0 * xr7
         yr73 = 3.0d0 * yr7
         xr9x7 = xr9 - xr73
         yr9y7 = yr9 - yr73
         t21 = xr9x7 - xr73 + rr53
         t22 = yrx * xr9x7
         t23 = zrx * xr9x7
         t24 = xyr9 - xr7 - yr7 + rr5
         t25 = zry * (xyr9-yr7)
         t26 = -t21 - t24
         t27 = xry * yr9y7
         t28 = zrx * (xyr9-xr7)
         t29 = -t22 - t27
         t31 = yr9y7 - yr73 + rr53
         t32 = zry * yr9y7
         t33 = -t24 - t31
c
c     compute the fifth order T2 matrix elements
c
         r4 = r2 * r2
         rr5 = rr2 * rr3
         rr7 = rr2 * rr5
         rr9 = rr2 * rr7
         rr11 = rr2 * rr9
         x2 = xr * xr
         y2 = yr * yr
         z2 = zr * zr
         x4 = x2 * x2
         y4 = y2 * y2
         x2r2 = x2 * r2
         y2r2 = y2 * r2
         x2y2 = x2 * y2
         xrr11 = xr * rr11
         yrr11 = yr * rr11
         zrr11 = zr * rr11
         xyzrr11 = 315.0d0 * xyz * rr11
         r445 = 45.0d0 * r4
         n945x4 = -945.0d0 * x4
         n945y4 = -945.0d0 * y4
         n945x2y2 = -945.0d0 * x2y2
         nnn1 = n945x4 + 630.0d0*x2r2 - r445
         nnn2 = n945y4 + 630.0d0*y2r2 - r445
         r4225 = 225.0d0 * r4
         t36 = (n945x4+1050.0d0*x2r2-r4225) * xrr11
         t37 = nnn1 * yrr11
         t38 = nnn1 * zrr11
         t39 = (n945x2y2+105.0d0*x2r2+315.0d0*y2r2-r445) * xrr11
         t40 = (r2-3.0d0*x2) * xyzrr11
         t41 = -t36 - t39
         t42 = (n945x2y2+105.0d0*y2r2+315.0d0*x2r2-r445) * yrr11
         t43 = (n945x2y2+105.0d0*(x2r2+y2r2)-15.0d0*r4) * zrr11
         t44 = -t37 - t42
         t46 = nnn2 * xrr11
         t47 = (r2-3.0d0*y2) * xyzrr11
         t48 = -t39 - t46
         t51 = (n945y4+1050.0d0*y2r2-r4225) * yrr11
         t52 = nnn2 * zrr11
         t53 = -t42 - t51
      end if
c
c     get the M-M, M-D and D-D parts of the multipole energy
c
      i3k3 = i3 * k3
      w1 = i0 * k0
      w2 = i0*k1 - i1*k0
      w3 = i0*k2 - i2*k0
      w4 = i0*k3 - i3*k0
      w5 = i3k3 - i1*k1
      w6 = -i1*k2 - i2*k1
      w7 = -i1*k3 - i3*k1
      w8 = i3k3 - i2*k2
      w9 = -i2*k3 - i3*k2
      eik = w1*t1 + w2*t2 + w3*t3 + w4*t4 + w5*t5
     &         + w6*t6 + w7*t7 + w8*t8 + w9*t9
      dm(1) = w1*t2 + w2*t5 + w3*t6 + w4*t7 + w5*t11
     &           + w6*t12 + w7*t13 + w8*t14 + w9*t15
      dm(2) = w1*t3 + w2*t6 + w3*t8 + w4*t9 + w5*t12
     &           + w6*t14 + w7*t15 + w8*t17 + w9*t18
      dm(3) = w1*t4 + w2*t7 + w3*t9 + w4*t10 + w5*t13
     &           + w6*t15 + w7*t16 + w8*t18 + w9*t19
c
c     get the M-Q and D-Q parts of the multipole energy
c
      if (kquad) then
         k4k9 = k4 - k9
         k7k9 = k7 - k9
         i3k6 = 2.0d0 * i3 * k6
         i3k8 = 2.0d0 * i3 * k8
         w1 = i0 * k4k9
         w2 = 2.0d0 * i0 * k5
         w3 = 2.0d0 * i0 * k6
         w4 = i0 * k7k9
         w5 = 2.0d0 * i0 * k8
         w6 = i3k6 - i1*k4k9
         w7 = i3k8 - i2*k4k9 - 2.0d0*i1*k5
         w8 = -2.0d0*i1*k6 - i3*k4k9
         w9 = i3k6 - i1*k7k9 - 2.0d0*i2*k5
         w10 = -2.0d0 * (i1*k8 + i2*k6 + i3*k5)
         w11 = i3k8 - i2*k7k9
         w12 = -2.0d0*i2*k8 - i3*k7k9
         eik = eik + w1*t5 + w2*t6 + w3*t7 + w4*t8
     &            + w5*t9 + w6*t11 + w7*t12 + w8*t13
     &            + w9*t14 + w10*t15 + w11*t17 + w12*t18
         dm(1) = dm(1) + w1*t11 + w2*t12 + w3*t13 + w4*t14
     &              + w5*t15 + w6*t21 + w7*t22 + w8*t23
     &              + w9*t24 + w10*t25 + w11*t27 + w12*t28
         dm(2) = dm(2) + w1*t12 + w2*t14 + w3*t15 + w4*t17
     &              + w5*t18 + w6*t22 + w7*t24 + w8*t25
     &              + w9*t27 + w10*t28 + w11*t31 + w12*t32
         dm(3) = dm(3) + w1*t13 + w2*t15 + w3*t16 + w4*t18
     &              + w5*t19 + w6*t23 + w7*t25 + w8*t26
     &              + w9*t28 + w10*t29 + w11*t32 + w12*t33
      end if
c
c     get the M-Q and D-Q parts of the multipole energy
c
      if (iquad) then
         i4i9 = i4 - i9
         i7i9 = i7 - i9
         i6k3 = -2.0d0 * i6 * k3
         i8k3 = -2.0d0 * i8 * k3
         w1 = i4i9 * k0
         w2 = 2.0d0 * i5 * k0
         w3 = 2.0d0 * i6 * k0
         w4 = i7i9 * k0
         w5 = 2.0d0 * i8 * k0
         w6 = i4i9*k1 + i6k3
         w7 = i8k3 + i4i9*k2 + 2.0d0*i5*k1
         w8 = 2.0d0*i6*k1 + i4i9*k3
         w9 = i6k3 + i7i9*k1 + 2.0d0*i5*k2
         w10 = 2.0d0 * (i8*k1 + i6*k2 + i5*k3)
         w11 = i8k3 + i7i9*k2
         w12 = 2.0d0*i8*k2 + i7i9*k3
         eik = eik + w1*t5 + w2*t6 + w3*t7 + w4*t8
     &            + w5*t9 + w6*t11 + w7*t12 + w8*t13
     &            + w9*t14 + w10*t15 + w11*t17 + w12*t18
         dm(1) = dm(1) + w1*t11 + w2*t12 + w3*t13 + w4*t14
     &              + w5*t15 + w6*t21 + w7*t22 + w8*t23
     &              + w9*t24 + w10*t25 + w11*t27 + w12*t28
         dm(2) = dm(2) + w1*t12 + w2*t14 + w3*t15 + w4*t17
     &              + w5*t18 + w6*t22 + w7*t24 + w8*t25
     &              + w9*t27 + w10*t28 + w11*t31 + w12*t32
         dm(3) = dm(3) + w1*t13 + w2*t15 + w3*t16 + w4*t18
     &              + w5*t19 + w6*t23 + w7*t25 + w8*t26
     &              + w9*t28 + w10*t29 + w11*t32 + w12*t33
      end if
c
c     get the Q-Q part of the multipole interaction energy
c
      if (iquad .and. kquad) then
         ik49 = i4*k9 + i9*k4
         ik59 = i5*k9 + i9*k5
         ik69 = i6*k9 + i9*k6
         ik79 = i7*k9 + i9*k7
         ik89 = i8*k9 + i9*k8
         ik68 = 2.0d0 * (i6*k8 + i8*k6)
         i6k6 = i6 * k6
         i8k8 = i8 * k8
         i9k9 = i9 * k9
         w1 = i4*k4 - ik49 + i9k9 - 4.0d0*i6k6
         w2 = 2.0d0 * (i5*k4 + i4*k5 - ik68 - ik59)
         w3 = 2.0d0 * (i4*k6 + i6*k4 - ik69)
         w4 = i4*k7 + i7*k4 - ik49 - ik79 + 2.0d0*i9k9
     &           +4.0d0*(i5*k5-i6k6-i8k8)
         w5 = 2.0d0 * (i4*k8 + i8*k4 + 2.0d0*(i5*k6+i6*k5) - ik89)
         w6 = 2.0d0 * (i5*k7 + i7*k5 - ik68 - ik59)
         w7 = 2.0d0 * (i6*k7 + i7*k6 + 2.0d0*(i5*k8+i8*k5) - ik69)
         w8 = i7*k7 - ik79 + i9k9 - 4.0d0*i8k8
         w9 = 2.0d0 * (i7*k8 + i8*k7 - ik89)
         eik = eik + w1*t21 + w2*t22 + w3*t23 + w4*t24
     &            + w5*t25 + w6*t27 + w7*t28 + w8*t31 + w9*t32
         dm(1) = dm(1) + w1*t36 + w2*t37 + w3*t38 + w4*t39
     &              + w5*t40 + w6*t42 + w7*t43 + w8*t46 + w9*t47
         dm(2) = dm(2) + w1*t37 + w2*t39 + w3*t40 + w4*t42
     &              + w5*t43 + w6*t46 + w7*t47 + w8*t51 + w9*t52
         dm(3) = dm(3) + w1*t38 + w2*t40 + w3*t41 + w4*t43
     &              + w5*t44 + w6*t47 + w7*t48 + w8*t52 + w9*t53
      end if
c
c     get the (dM2/dx)*T*M1 terms for dipoles at both sites
c
      do i = 1, 3
         do j = 1, 3
            di1 = dpole(2,j,i,ii)
            di2 = dpole(3,j,i,ii)
            di3 = dpole(4,j,i,ii)
            dk1 = dpole(2,j,i,kk)
            dk2 = dpole(3,j,i,kk)
            dk3 = dpole(4,j,i,kk)
            w1 = k3 * di3
            dmi(j,i) = -k0*di1*t2 - k0*di2*t3 - k0*di3*t4
     &                    + (w1-k1*di1)*t5 - (k2*di1+k1*di2)*t6
     &                    - (k1*di3+k3*di1)*t7 + (w1-k2*di2)*t8
     &                    - (k2*di3+k3*di2)*t9
            w1 = dk3 * i3
            dmk(j,i) = dk1*i0*t2 + dk2*i0*t3 + dk3*i0*t4
     &                    + (w1-dk1*i1)*t5 - (dk2*i1+dk1*i2)*t6
     &                    - (dk1*i3+dk3*i1)*t7 + (w1-dk2*i2)*t8
     &                    - (dk2*i3+dk3*i2)*t9
            w1 = v3 * di3
            dpi(j,i) = (w1-v1*di1)*t5 - (v2*di1+v1*di2)*t6
     &                    - (v1*di3+v3*di1)*t7 + (w1-v2*di2)*t8
     &                    - (v3*di2+v2*di3)*t9
            w1 = u3 * dk3
            dpk(j,i) = (w1-u1*dk1)*t5 - (u2*dk1+u1*dk2)*t6
     &                    - (u1*dk3+u3*dk1)*t7 + (w1-u2*dk2)*t8
     &                    - (u3*dk2+u2*dk3)*t9
c
c     get the (dM2/dx)*T*M1 terms for quadrupole at first site
c
            if (iquad) then
               di4 = dpole(5,j,i,ii)
               di5 = dpole(6,j,i,ii)
               di6 = dpole(7,j,i,ii)
               di7 = dpole(9,j,i,ii)
               di8 = dpole(10,j,i,ii)
               di9 = dpole(13,j,i,ii)
               w1 = i9 - i4
               w2 = i9 - i7
               w3 = i6 * dk3
               w4 = i8 * dk3
               dmk(j,i) = dmk(j,i) - (w1*dk1+2.0d0*w3)*t11
     &                       - (w1*dk2+2.0d0*(w4-dk1*i5))*t12
     &                       - (w1*dk3-2.0d0*dk1*i6)*t13
     &                       - (w2*dk1+2.0d0*(w3-dk2*i5))*t14
     &                       + 2.0d0*(dk3*i5+dk2*i6+dk1*i8)*t15
     &                       - (w2*dk2+2.0d0*w4)*t17
     &                       - (w2*dk3-2.0d0*dk2*i8)*t18
               w1 = di9 - di4
               w2 = di9 - di7
               w3 = di6 * k3
               w4 = di8 * k3
               dmi(j,i) = dmi(j,i) - k0*(w1*t5+w2*t8)
     &                       + 2.0d0*k0*(di5*t6+di6*t7+di8*t9)
     &                       - (w1*k1+2.0d0*w3)*t11
     &                       - (w1*k2+2.0d0*(w4-k1*di5))*t12
     &                       - (w1*k3-2.0d0*k1*di6)*t13
     &                       - (w2*k1+2.0d0*(w3-k2*di5))*t14
     &                       + 2.0d0*(k3*di5+k2*di6+k1*di8)*t15
     &                       - (w2*k2+2.0d0*w4)*t17
     &                       - (w2*k3-2.0d0*k2*di8)*t18
c              w1 = di9 - di4
c              w2 = di9 - di7
               w3 = di6 * v3
               w4 = di8 * v3
               dpi(j,i) = dpi(j,i) - (w1*v1+2.0d0*w3)*t11
     &                       - (w1*v2+2.0d0*(w4-v1*di5))*t12
     &                       - (w1*v3-2.0d0*v1*di6)*t13
     &                       - (w2*v1+2.0d0*(w3-v2*di5))*t14
     &                       + 2.0d0*(v3*di5+v2*di6+v1*di8)*t15
     &                       - (w2*v2+2.0d0*w4)*t17
     &                       - (w2*v3-2.0d0*v2*di8)*t18
            end if
c
c     get the (dM2/dx)*T*M1 terms for quadrupole at second site
c
            if (kquad) then
               dk4 = dpole(5,j,i,kk)
               dk5 = dpole(6,j,i,kk)
               dk6 = dpole(7,j,i,kk)
               dk7 = dpole(9,j,i,kk)
               dk8 = dpole(10,j,i,kk)
               dk9 = dpole(13,j,i,kk)
               w1 = k9 - k4
               w2 = k9 - k7
               w3 = k6 * di3
               w4 = k8 * di3
               dmi(j,i) = dmi(j,i)
     &                       + (w1*di1+2.0d0*w3)*t11
     &                       + (w1*di2+2.0d0*(w4-k5*di1))*t12
     &                       + (w1*di3-2.0d0*k6*di1)*t13
     &                       + (w2*di1+2.0d0*(w3-k5*di2))*t14
     &                       - 2.0d0*(k5*di3+k6*di2+k8*di1)*t15
     &                       + (w2*di2+2.0d0*w4)*t17
     &                       + (w2*di3-2.0d0*k8*di2)*t18
               w1 = dk9 - dk4
               w2 = dk9 - dk7
               w3 = dk6 * i3
               w4 = dk8 * i3
               dmk(j,i) = dmk(j,i) - i0*(w1*t5+w2*t8)
     &                       + 2.0d0*i0*(dk5*t6+dk6*t7+dk8*t9)
     &                       + (w1*i1+2.0d0*w3)*t11
     &                       + (w1*i2+2.0d0*(w4-dk5*i1))*t12
     &                       + (w1*i3-2.0d0*dk6*i1)*t13
     &                       + (w2*i1+2.0d0*(w3-dk5*i2))*t14
     &                       - 2.0d0*(dk8*i1+dk5*i3+dk6*i2)*t15
     &                       + (w2*i2+2.0d0*w4)*t17
     &                       + (w2*i3-2.0d0*dk8*i2)*t18
c              w1 = dk9 - dk4
c              w2 = dk9 - dk7
               w3 = dk6 * u3
               w4 = dk8 * u3
               dpk(j,i) = dpk(j,i) + (w1*u1+2.0d0*w3)*t11
     &                       + (w1*u2+2.0d0*(w4-u1*dk5))*t12
     &                       + (w1*u3-2.0d0*u1*dk6)*t13
     &                       + (w2*u1+2.0d0*(w3-u2*dk5))*t14
     &                       - 2.0d0*(u1*dk8+u2*dk6+u3*dk5)*t15
     &                       + (w2*u2+2.0d0*w4)*t17
     &                       + (w2*u3-2.0d0*u2*dk8)*t18
            end if
c
c     get the (dM2/dx)*T*M1 terms for quadrupoles at both sites
c
            if (iquad .and. kquad) then
               w1 = di9 - di4
               w2 = di9 - di7
               w3 = di5*k9 + 2.0d0*(di6*k8+di8*k6)
               w4 = di6 * k6
               w5 = di6 * k9
               w6 = di8 * k8
               w7 = di8 * k9
               dmi(j,i) = dmi(j,i) + (w1*(k9-k4)-4.0d0*w4)*t21
     &                       + 2.0d0*(di5*k4-w1*k5-w3)*t22
     &                       + 2.0d0*(di6*k4-w1*k6-w5)*t23
     &                       + (4.0d0*(di5*k5-w4-w6)-w1*k7-w2*k4
     &                            +(2.0d0*di9-di4-di7)*k9)*t24
     &                       + 2.0d0*(di8*k4-w1*k8-w7
     &                            +2.0d0*(di6*k5+di5*k6))*t25
     &                       + 2.0d0*(di5*k7-w2*k5-w3)*t27
     &                       + 2.0d0*(di6*k7-w2*k6-w5
     &                            +2.0d0*(di8*k5+di5*k8))*t28
     &                       + (w2*(k9-k7)-4.0d0*w6)*t31
     &                       + 2.0d0*(di8*k7-w2*k8-w7)*t32
               w1 = dk9 - dk4
               w2 = dk9 - dk7
               w3 = dk5*i9 + 2.0d0*(dk6*i8+dk8*i6)
               w4 = dk6 * i6
               w5 = dk6 * i9
               w6 = dk8 * i8
               w7 = dk8 * i9
               dmk(j,i) = dmk(j,i) + (w1*(i9-i4)-4.0d0*w4)*t21
     &                       + 2.0d0*(dk5*i4-w1*i5-w3)*t22
     &                       + 2.0d0*(dk6*i4-w1*i6-w5)*t23
     &                       + (4.0d0*(dk5*i5-w4-w6)-w1*i7-w2*i4
     &                            +(2.0d0*dk9-dk4-dk7)*i9)*t24
     &                       + 2.0d0*(dk8*i4-w1*i8-w7
     &                            +2.0d0*(dk6*i5+dk5*i6))*t25
     &                       + 2.0d0*(dk5*i7-w2*i5-w3)*t27
     &                       + 2.0d0*(dk6*i7-w2*i6-w5
     &                            +2.0d0*(dk8*i5+dk5*i8))*t28
     &                       + (w2*(i9-i7)-4.0d0*w6)*t31
     &                       + 2.0d0*(dk8*i7-w2*i8-w7)*t32
            end if
         end do 
      end do
c
c     get indD-M, indD-D and indD-Q polarization energy and gradients
c
      k3u3 = k3 * u3
      w1 = k0 * u1
      w2 = k0 * u2
      w3 = k0 * u3
      w4 = k3u3 - k1*u1
      w5 = k1*u2 + k2*u1
      w6 = k3*u1 + k1*u3
      w7 = k3u3 - k2*u2
      w8 = k2*u3 + k3*u2
      ei = -w1*t2 - w2*t3 - w3*t4 + w4*t5 - w5*t6
     &        - w6*t7 + w7*t8 - w8*t9
      dp(1) = -w1*t5 - w2*t6 - w3*t7 + w4*t11 - w5*t12
     &           - w6*t13 + w7*t14 - w8*t15
      dp(2) = -w1*t6 - w2*t8 - w3*t9 + w4*t12 - w5*t14
     &           - w6*t15 + w7*t17 - w8*t18
      dp(3) = -w1*t7 - w2*t9 - w3*t10 + w4*t13 - w5*t15
     &           - w6*t16 + w7*t18 - w8*t19
      if (kquad) then
         k6u3 = k6 * u3
         k8u3 = k8 * u3
         w1 = -k4k9*u1 + 2.0d0*k6u3
         w2 = -k4k9*u2 + 2.0d0*(k8u3-k5*u1)
         w3 = -k4k9*u3 - 2.0d0*k6*u1
         w4 = -k7k9*u1 + 2.0d0*(k6u3-k5*u2)
         w5 = 2.0d0 * (k5*u3 + k6*u2 + k8*u1)
         w6 = -k7k9*u2 + 2.0d0*k8u3
         w7 = -k7k9*u3 - 2.0d0*k8*u2
         ei = ei + w1*t11 + w2*t12 + w3*t13 + w4*t14
     &           - w5*t15 + w6*t17 + w7*t18
         dp(1) = dp(1) + w1*t21 + w2*t22 + w3*t23
     &              + w4*t24 - w5*t25 + w6*t27 + w7*t28
         dp(2) = dp(2) + w1*t22 + w2*t24 + w3*t25
     &              + w4*t27 - w5*t28 + w6*t31 + w7*t32
         dp(3) = dp(3) + w1*t23 + w2*t25 + w3*t26
     &              + w4*t28 - w5*t29 + w6*t32 + w7*t33
      end if
c
c     get M-indD, D-indD and Q-indD polarization energy and gradients
c
      i3v3 = i3 * v3
      w1 = i0 * v1
      w2 = i0 * v2
      w3 = i0 * v3
      w4 = i3v3 - i1*v1
      w5 = i1*v2 + i2*v1
      w6 = i3*v1 + i1*v3
      w7 = i3v3 - i2*v2
      w8 = i2*v3 + i3*v2
      ek = w1*t2 + w2*t3 + w3*t4 + w4*t5 - w5*t6
     &        - w6*t7 + w7*t8 - w8*t9
      dp(1) = dp(1) + w1*t5 + w2*t6 + w3*t7 + w4*t11 - w5*t12
     &           - w6*t13 + w7*t14 - w8*t15
      dp(2) = dp(2) + w1*t6 + w2*t8 + w3*t9 + w4*t12 - w5*t14
     &           - w6*t15 + w7*t17 - w8*t18
      dp(3) = dp(3) + w1*t7 + w2*t9 + w3*t10 + w4*t13 - w5*t15
     &           - w6*t16 + w7*t18 - w8*t19
      if (iquad) then
         i6v3 = i6 * v3
         i8v3 = i8 * v3
         w1 = i4i9*v1 - 2.0d0*i6v3
         w2 = i4i9*v2 + 2.0d0*(i5*v1-i8v3)
         w3 = i4i9*v3 + 2.0d0*i6*v1
         w4 = i7i9*v1 + 2.0d0*(i5*v2-i6v3)
         w5 = 2.0d0 * (i5*v3 + i6*v2 + i8*v1)
         w6 = i7i9*v2 - 2.0d0*i8v3
         w7 = i7i9*v3 + 2.0d0*i8*v2
         ek = ek + w1*t11 + w2*t12 + w3*t13 + w4*t14
     &           + w5*t15 + w6*t17 + w7*t18
         dp(1) = dp(1) + w1*t21 + w2*t22 + w3*t23
     &              + w4*t24 + w5*t25 + w6*t27 + w7*t28
         dp(2) = dp(2) + w1*t22 + w2*t24 + w3*t25
     &              + w4*t27 + w5*t28 + w6*t31 + w7*t32
         dp(3) = dp(3) + w1*t23 + w2*t25 + w3*t26
     &              + w4*t28 + w5*t29 + w6*t32 + w7*t33
      end if
c
c     compute the mutual polarization induced dipole gradient terms
c
      if (poltyp .eq. 'MUTUAL') then
         w1 = -u1*v1 + u3*v3
         w2 = -u2*v1 - u1*v2
         w3 = -u3*v1 - u1*v3
         w4 = -u2*v2 + u3*v3
         w5 = -u3*v2 - u2*v3
         utu = w1*t5 + w2*t6 + w3*t7 + w4*t8 + w5*t9
         dp(1) = dp(1) + w1*t11 + w2*t12 + w3*t13 + w4*t14 + w5*t15
         dp(2) = dp(2) + w1*t12 + w2*t14 + w3*t15 + w4*t17 + w5*t18
         dp(3) = dp(3) + w1*t13 + w2*t15 + w3*t16 + w4*t18 + w5*t19
      end if
c
c     apply a damping factor to polarization energy and derivatives
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         term = -pgamma * (r/damp)**3
         if (term .gt. -50.0d0) then
            term = exp(term)
            ddamp = (3.0d0*pgamma*r2/damp**3) * term
            damp = 1.0d0 - term
            de = (ei+ek+utu) * ddamp
            ei = ei * damp
            ek = ek * damp
            dp(1) = dp(1)*damp + de*(xr/r)
            dp(2) = dp(2)*damp + de*(yr/r)
            dp(3) = dp(3)*damp + de*(zr/r)
            do i = 1, 3
               do j = 1, 3
                  dpi(j,i) = dpi(j,i) * damp
                  dpk(j,i) = dpk(j,i) * damp
               end do
            end do
         end if
      end if
c
c     final polarization energy is half of value computed above
c
      ei = 0.5d0 * ei
      ek = 0.5d0 * ek
      return
      end
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole2  --  multipole & polarization Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole2" calculates second derivatives of the multipole
c     and dipole polarization energy for a single atom at a time
c
c
      subroutine empole2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'hessn.i'
      integer i,j,k
      real*8 eps,old
      real*8 d0(3,maxatm)
c
c
c     set the stepsize to be used for numerical derivatives
c
      eps = 1.0d-7
c
c     get multipole first derivatives for the base structure
c
      call empole2b (i)
      do k = 1, n
         do j = 1, 3
            d0(j,k) = dem(j,k) + dep(j,k)
         end do
      end do
c
c     find numerical x-components via perturbed structures
c
      old = x(i)
      x(i) = x(i) + eps
      call empole2b (i)
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical y-components via perturbed structures
c
      old = y(i)
      y(i) = y(i) + eps
      call empole2b (i)
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical z-components via perturbed structures
c
      old = z(i)
      z(i) = z(i) + eps
      call empole2b (i)
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole2b  --  mpole & polar Hessian; numerical  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole2b" computes multipole and dipole polarization first
c     derivatives for a single angle with respect to Cartesian
c     coordinates; used to get finite difference second derivatives
c
c     note that since polarization effects are many body, it is not
c     really correct to neglect interactions where "iatom" is not
c     directly involved as a multipole site or local coordinate axis;
c     however, other sites are neglected in this version via the
c     "ipart" and "kpart" checks to quickly get approximate values
c
c
      subroutine empole2b (iatom)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'group.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,m
      integer jj,iatom
      integer ii,iz,ix
      integer kk,kz,kx
      integer skip(maxatm)
      real*8 eik,ei,ek,de
      real*8 f,fik,fgrp
      real*8 shift,taper,dtaper
      real*8 trans,dtrans
      real*8 xi,yi,zi,xr,yr,zr
      real*8 r,r2,r3,r4,r5,r6,r7
      real*8 a(3,3),d(3,3,3,3)
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 dm(3),dp(3),utu
      real*8 dmi(3,3),dmk(3,3)
      real*8 dpi(3,3),dpk(3,3)
      logical proceed,iuse,kuse,ipart,kpart
c
c
c     zero out the multipole and polarization first derivatives
c
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
         end do
      end do
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     rotate multipole components and get rotation derivatives
c
      do i = 1, npole
cjrs
         call rotmatt (i,a)
cjrs
         call rotpole (i,a)
         call drotmat (i,d)
         do j = 1, 3
            do k = 1, 3
               call drotpole (i,a,d,j,k)
            end do
         end do
      end do
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     compute and partition multipole interaction energy
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         iz = zaxis(ii)
         ix = xaxis(ii)
         ipart = (i.eq.iatom .or. iz.eq.iatom .or. ix.eq.iatom)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * chg12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * chg13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * chg14use
         end do
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kpart = (k.eq.iatom .or. kz.eq.iatom .or. kx.eq.iatom)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = ((ipart.or.kpart) .and.
     &                                  (iuse.or.kuse))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, maxpole
                     rpk(j) = rpole(j,kk)
                  end do
                  do j = 1, 3
                     indk(j) = uind(j,kk)
                  end do
                  r = sqrt(r2)
                  call empik1 (ii,kk,xr,yr,zr,r,r2,rpi,
     &                         rpk,indi,indk,eik,ei,ek,
     &                         dm,dmi,dmk,dp,dpi,dpk,utu)
                  fik = f
                  if (skip(k) .eq. -i)  fik = fik / chgscale
                  eik = fik * eik
                  ei = fik * ei
                  ek = fik * ek
                  utu = fik * utu
                  do j = 1, 3
                     dm(j) = fik * dm(j)
                     dp(j) = fik * dp(j)
                     do jj = 1, 3
                        dmi(jj,j) = fik * dmi(jj,j)
                        dmk(jj,j) = fik * dmk(jj,j)
                        dpi(jj,j) = fik * dpi(jj,j)
                        dpk(jj,j) = fik * dpk(jj,j)
                     end do
                  end do
c
c     use shifted energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     fik = fik * rpi(1) * rpk(1)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                     dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                                  + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                              + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                     shift = fik / (0.5d0*(off+cut))
                     eik = eik - shift
                     de = eik * dtaper + dtrans
                     dm(1) = dm(1)*taper + de*(xr/r)
                     dm(2) = dm(2)*taper + de*(yr/r)
                     dm(3) = dm(3)*taper + de*(zr/r)
                     de = (2.0d0*(ei+ek)+utu) * dtaper
                     dp(1) = dp(1)*taper + de*(xr/r)
                     dp(2) = dp(2)*taper + de*(yr/r)
                     dp(3) = dp(3)*taper + de*(zr/r)
                     do j = 1, 3
                        do jj = 1, 3
                           dmi(jj,j) = dmi(jj,j) * taper
                           dmk(jj,j) = dmk(jj,j) * taper
                           dpi(jj,j) = dpi(jj,j) * taper
                           dpk(jj,j) = dpk(jj,j) * taper
                        end do
                     end do
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     do j = 1, 3
                        dm(j) = dm(j) * fgrp
                        dp(j) = dp(j) * fgrp
                        do jj = 1, 3
                           dmi(jj,j) = dmi(jj,j) * fgrp
                           dmk(jj,j) = dmk(jj,j) * fgrp
                           dpi(jj,j) = dpi(jj,j) * fgrp
                           dpk(jj,j) = dpk(jj,j) * fgrp
                        end do
                     end do
                  end if
c
c     increment the multipole first derivative expressions
c
                  dem(1,i) = dem(1,i) - dm(1) + dmi(1,1)
                  dem(2,i) = dem(2,i) - dm(2) + dmi(1,2)
                  dem(3,i) = dem(3,i) - dm(3) + dmi(1,3)
                  dem(1,iz) = dem(1,iz) + dmi(2,1)
                  dem(2,iz) = dem(2,iz) + dmi(2,2)
                  dem(3,iz) = dem(3,iz) + dmi(2,3)
                  dem(1,ix) = dem(1,ix) + dmi(3,1)
                  dem(2,ix) = dem(2,ix) + dmi(3,2)
                  dem(3,ix) = dem(3,ix) + dmi(3,3)
                  dem(1,k) = dem(1,k) + dm(1) + dmk(1,1)
                  dem(2,k) = dem(2,k) + dm(2) + dmk(1,2)
                  dem(3,k) = dem(3,k) + dm(3) + dmk(1,3)
                  dem(1,kz) = dem(1,kz) + dmk(2,1)
                  dem(2,kz) = dem(2,kz) + dmk(2,2)
                  dem(3,kz) = dem(3,kz) + dmk(2,3)
                  dem(1,kx) = dem(1,kx) + dmk(3,1)
                  dem(2,kx) = dem(2,kx) + dmk(3,2)
                  dem(3,kx) = dem(3,kx) + dmk(3,3)
c
c     increment the polarization first derivative expressions
c
                  dep(1,i) = dep(1,i) - dp(1) + dpi(1,1)
                  dep(2,i) = dep(2,i) - dp(2) + dpi(1,2)
                  dep(3,i) = dep(3,i) - dp(3) + dpi(1,3)
                  dep(1,iz) = dep(1,iz) + dpi(2,1)
                  dep(2,iz) = dep(2,iz) + dpi(2,2)
                  dep(3,iz) = dep(3,iz) + dpi(2,3)
                  dep(1,ix) = dep(1,ix) + dpi(3,1)
                  dep(2,ix) = dep(2,ix) + dpi(3,2)
                  dep(3,ix) = dep(3,ix) + dpi(3,3)
                  dep(1,k) = dep(1,k) + dp(1) + dpk(1,1)
                  dep(2,k) = dep(2,k) + dp(2) + dpk(1,2)
                  dep(3,k) = dep(3,k) + dp(3) + dpk(1,3)
                  dep(1,kz) = dep(1,kz) + dpk(2,1)
                  dep(2,kz) = dep(2,kz) + dpk(2,2)
                  dep(3,kz) = dep(3,kz) + dpk(2,3)
                  dep(1,kx) = dep(1,kx) + dpk(3,1)
                  dep(2,kx) = dep(2,kx) + dpk(3,2)
                  dep(3,kx) = dep(3,kx) + dpk(3,3)
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         iz = zaxis(ii)
         ix = xaxis(ii)
         ipart = (i.eq.iatom .or. iz.eq.iatom .or. ix.eq.iatom)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kpart = (k.eq.iatom .or. kz.eq.iatom .or. kx.eq.iatom)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = ((ipart.or.kpart) .and.
     &                                  (iuse.or.kuse))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do m = 1, ncell
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     do j = 1, 3
                        indk(j) = uind(j,kk)
                     end do
                     r = sqrt(r2)
                     call empik1 (ii,kk,xr,yr,zr,r,r2,rpi,
     &                            rpk,indi,indk,eik,ei,ek,
     &                            dm,dmi,dmk,dp,dpi,dpk,utu)
                     fik = f
                     eik = fik * eik
                     ei = fik * ei
                     ek = fik * ek
                     utu = fik * utu
                     do j = 1, 3
                        dm(j) = fik * dm(j)
                        dp(j) = fik * dp(j)
                        do jj = 1, 3
                           dmi(jj,j) = fik * dmi(jj,j)
                           dmk(jj,j) = fik * dmk(jj,j)
                           dpi(jj,j) = fik * dpi(jj,j)
                           dpk(jj,j) = fik * dpk(jj,j)
                        end do
                     end do
c
c     use shifted energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        fik = fik * rpi(1) * rpk(1)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                                 + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                     + f3*r3 + f2*r2 + f1*r + f0)
                        dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                                     + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                                 + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                        shift = fik / (0.5d0*(off+cut))
                        eik = eik - shift
                        de = eik * dtaper + dtrans
                        dm(1) = dm(1)*taper + de*(xr/r)
                        dm(2) = dm(2)*taper + de*(yr/r)
                        dm(3) = dm(3)*taper + de*(zr/r)
                        de = (2.0d0*(ei+ek)+utu) * dtaper
                        dp(1) = dp(1)*taper + de*(xr/r)
                        dp(2) = dp(2)*taper + de*(yr/r)
                        dp(3) = dp(3)*taper + de*(zr/r)
                        do j = 1, 3
                           do jj = 1, 3
                              dmi(jj,j) = dmi(jj,j) * taper
                              dmk(jj,j) = dmk(jj,j) * taper
                              dpi(jj,j) = dpi(jj,j) * taper
                              dpk(jj,j) = dpk(jj,j) * taper
                           end do
                        end do
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        do j = 1, 3
                           dm(j) = dm(j) * fgrp
                           dp(j) = dp(j) * fgrp
                           do jj = 1, 3
                              dmi(jj,j) = dmi(jj,j) * fgrp
                              dmk(jj,j) = dmk(jj,j) * fgrp
                              dpi(jj,j) = dpi(jj,j) * fgrp
                              dpk(jj,j) = dpk(jj,j) * fgrp
                           end do
                        end do
                     end if
c
c     increment the multipole first derivative expressions
c
                     dem(1,i) = dem(1,i) - dm(1) + dmi(1,1)
                     dem(2,i) = dem(2,i) - dm(2) + dmi(1,2)
                     dem(3,i) = dem(3,i) - dm(3) + dmi(1,3)
                     dem(1,iz) = dem(1,iz) + dmi(2,1)
                     dem(2,iz) = dem(2,iz) + dmi(2,2)
                     dem(3,iz) = dem(3,iz) + dmi(2,3)
                     dem(1,ix) = dem(1,ix) + dmi(3,1)
                     dem(2,ix) = dem(2,ix) + dmi(3,2)
                     dem(3,ix) = dem(3,ix) + dmi(3,3)
                     if (i .ne. k) then
                        dem(1,k) = dem(1,k) + dm(1) + dmk(1,1)
                        dem(2,k) = dem(2,k) + dm(2) + dmk(1,2)
                        dem(3,k) = dem(3,k) + dm(3) + dmk(1,3)
                        dem(1,kz) = dem(1,kz) + dmk(2,1)
                        dem(2,kz) = dem(2,kz) + dmk(2,2)
                        dem(3,kz) = dem(3,kz) + dmk(2,3)
                        dem(1,kx) = dem(1,kx) + dmk(3,1)
                        dem(2,kx) = dem(2,kx) + dmk(3,2)
                        dem(3,kx) = dem(3,kx) + dmk(3,3)
                     end if
c
c     increment the polarization first derivative expressions
c
                     dep(1,i) = dep(1,i) - dp(1) + dpi(1,1)
                     dep(2,i) = dep(2,i) - dp(2) + dpi(1,2)
                     dep(3,i) = dep(3,i) - dp(3) + dpi(1,3)
                     dep(1,iz) = dep(1,iz) + dpi(2,1)
                     dep(2,iz) = dep(2,iz) + dpi(2,2)
                     dep(3,iz) = dep(3,iz) + dpi(2,3)
                     dep(1,ix) = dep(1,ix) + dpi(3,1)
                     dep(2,ix) = dep(2,ix) + dpi(3,2)
                     dep(3,ix) = dep(3,ix) + dpi(3,3)
                     if (i .ne. k) then
                        dep(1,k) = dep(1,k) + dp(1) + dpk(1,1)
                        dep(2,k) = dep(2,k) + dp(2) + dpk(1,2)
                        dep(3,k) = dep(3,k) + dp(3) + dpk(1,3)
                        dep(1,kz) = dep(1,kz) + dpk(2,1)
                        dep(2,kz) = dep(2,kz) + dpk(2,2)
                        dep(3,kz) = dep(3,kz) + dpk(2,3)
                        dep(1,kx) = dep(1,kx) + dpk(3,1)
                        dep(2,kx) = dep(2,kx) + dpk(3,2)
                        dep(3,kx) = dep(3,kx) + dpk(3,3)
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole interactions and dipole polarizability,
c     and also partitions the energy among the atoms
c
c
      subroutine empole3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moment.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,m
      integer ii,iz,ix
      integer kk,kz,kx
      integer skip(maxatm)
      real*8 eik,ei,ek,a(3,3)
      real*8 shift,taper,trans
      real*8 f,fik,fgrp
      real*8 xr,yr,zr,r
      real*8 r2,r3,r4,r5,r6,r7
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 weight,xcenter,ycenter,zcenter
      real*8 totchg,xsum,ysum,zsum
      logical header,huge,proceed,iuse,kuse
c
c
c     zero out multipole and polarization energy and partitioning
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
      header = .true.
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     rotate the multipole components into the global frame
c
      do i = 1, npole
cjrs
         call rotmatt (i,a)
cjrs
         call rotpole (i,a)
      end do
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     find the center of mass of the total system
c
      weight = 0.0d0
      xcenter = 0.0d0
      ycenter = 0.0d0
      zcenter = 0.0d0
      do i = 1, n
         weight = weight + mass(i)
         xcenter = xcenter + x(i)*mass(i)
         ycenter = ycenter + y(i)*mass(i)
         zcenter = zcenter + z(i)*mass(i)
      end do
      xcenter = xcenter / weight
      ycenter = ycenter / weight
      zcenter = zcenter / weight
c
c     get net charge and dipole components relative to center of mass
c
      totchg = 0.0d0
      xsum = 0.0d0
      ysum = 0.0d0
      zsum = 0.0d0
      do i = 1, npole
         k = ipole(i)
         totchg = totchg + rpole(1,i)
         xsum = xsum + (x(k)-xcenter)*rpole(1,i)
         ysum = ysum + (y(k)-ycenter)*rpole(1,i)
         zsum = zsum + (z(k)-zcenter)*rpole(1,i)
         xsum = xsum + rpole(2,i) + uind(1,i)
         ysum = ysum + rpole(3,i) + uind(2,i)
         zsum = zsum + rpole(4,i) + uind(3,i)
      end do
      netchg = netchg + totchg
      xdipole = xdipole + debye*xsum
      ydipole = ydipole + debye*ysum
      zdipole = zdipole + debye*zsum
c
c     calculate the multipole interaction energy term
c
      do ii = 1, npole-1
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * chg12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * chg13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * chg14use
         end do
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, maxpole
                     rpk(j) = rpole(j,kk)
                  end do
                  do j = 1, 3
                     indk(j) = uind(j,kk)
                  end do
                  r = sqrt(r2)
                  call empik (ii,kk,xr,yr,zr,r,rpi,rpk,
     &                          indi,indk,eik,ei,ek)
                  fik = f
                  if (skip(k) .eq. -i)  fik = fik / chgscale
                  eik = fik * eik
                  ei = fik * ei
                  ek = fik * ek
c
c     use shifted energy switching if near the cutoff distance
c
                  fik = fik * rpi(1) * rpk(1)
                  shift = fik / (0.5d0*(off+cut))
                  eik = eik - shift
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     eik = eik * taper + trans
                     ei = ei * taper
                     ek = ek * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     eik = eik * fgrp
                     ei = ei * fgrp
                     ek = ek * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  nem = nem + 1
                  em = em + eik
                  aem(i) = aem(i) + 0.5d0*eik
                  aem(k) = aem(k) + 0.5d0*eik
                  nep = nep + 2
                  ep = ep + ei + ek
                  aep(i) = aep(i) + 0.5d0*(ei+ek)
                  aep(k) = aep(k) + 0.5d0*(ei+ek)
c
c     print a warning if the energy of this interaction is large
c
                  huge = (max(abs(eik),abs(ei),abs(ek)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',11x,'Atom Names',
     &                             9x,'Distance',6x,'Energies',
     &                             ' (MPole, Pol1, Pol2)',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),r,eik,ei,ek
   20                format (' M-Pole   ',i5,'-',a3,1x,i5,'-',a3,
     &                          5x,f8.4,3f12.4)
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do m = 1, ncell
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     do j = 1, 3
                        indk(j) = uind(j,kk)
                     end do
                     r = sqrt(r2)
                     call empik (ii,kk,xr,yr,zr,r,rpi,rpk,
     &                             indi,indk,eik,ei,ek)
                     fik = f
                     eik = fik * eik
                     ei = fik * ei
                     ek = fik * ek
c
c     use shifted energy switching if near the cutoff distance
c
                     fik = fik * rpi(1) * rpk(1)
                     shift = fik / (0.5d0*(off+cut))
                     eik = eik - shift
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                        eik = eik * taper + trans
                        ei = ei * taper
                        ek = ek * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        eik = eik * fgrp
                        ei = ei * fgrp
                        ek = ek * fgrp
                     end if
c
c     increment the overall multipole and polarization energies
c
                     nem = nem + 1
                     nep = nep + 2
                     if (i .eq. k) then
                        em = em + 0.5d0*eik
                        aem(i) = aem(i) + 0.5d0*eik
                        ep = ep + ei
                        aep(i) = aep(i) + ei
                     else
                        em = em + eik
                        aem(i) = aem(i) + 0.5d0*eik
                        aem(k) = aem(k) + 0.5d0*eik
                        ep = ep + ei + ek
                        aep(i) = aep(i) + 0.5d0*(ei+ek)
                        aep(k) = aep(k) + 0.5d0*(ei+ek)
                     end if
c
c     print a warning if the energy of this interaction is large
c
                     huge = (max(abs(eik),abs(ei),abs(ek)) .gt. 100.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',11x,'Atom Names',
     &                                9x,'Distance',6x,'Energies',
     &                                ' (MPole, Pol1, Pol2)',/)
                        end if
                        write (iout,40)  i,name(i),k,name(k),r,eik,ei,ek
   40                   format (' M-Pole   ',i5,'-',a3,1x,i5,'-',a3,
     &                             ' (X) ',f8.4,3f12.4)
                     end if
                  end if
               end do
            end if
         end do
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
c     #############################################################
c     ##                                                         ##
c     ##  function energy  --  evaluates energy terms and total  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "energy" calls the subroutines to calculate the potential
c     energy terms and sums up to form the total energy
c
c
      function energy ()
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'warp.i'
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
      if (use_bond)  call ebond
      if (use_angle)  call eangle
      if (use_strbnd)  call estrbnd
      if (use_urey)  call eurey
      if (use_angang)  call eangang
      if (use_opbend)  call eopbend
      if (use_improp)  call eimprop
      if (use_imptor)  call eimptor
      if (use_tors)  call etors
      if (use_strtor)  call estrtor
c     if (use_tortor)  call etortor
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (use_lights) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj4
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck4
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb4
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal4
         else
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
         end if
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss
      end if
c
c     call the electrostatic energy component routines
c
      if (use_charge) then
         if (use_deform) then
            call echarge6
         else if (use_lights) then
            call echarge4
         else
            call echarge
         end if
      end if
      if (use_chgdpl)  call echgdpl
      if (use_dipole)  call edipole
      if (use_mpole .or. use_polar)  call empole
      if (use_rxnfld)  call erxnfld
c
c     call any miscellaneous energy component routines
c
      if (use_solv)  call esolv
      if (use_geom)  call egeom
cjrs
      if (use_extra)  call extrat
cjrs
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
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine eopbend  --  out-of-plane bending energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "eopbend" computes the out-of-plane bend potential energy at
c     trigonal centers as per Allinger's MM2 and MM3 forcefields
c
c
      subroutine eopbend
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'opbend.i'
      include 'usage.i'
      integer i,iopbend
      integer ia,ib,ic,id
      real*8 e,force,angle
      real*8 dot,cosine,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd,rbd2
      real*8 xcd,ycd,zcd
      real*8 xpd,ypd,zpd,rpd2
      real*8 xt,yt,zt,rt2,delta
      logical proceed
c
c
c     zero out the out-of-plane bending energy component
c
      eopb = 0.0d0
c
c     calculate the out-of-plane bending energy term
c
      do iopbend = 1, nopbend
         i = iopb(iopbend)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         force = kopb(iopbend)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
c
c     compute the out-of-plane bending angle
c
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xbd = xib - xid
            ybd = yib - yid
            zbd = zib - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            xt = yad*zcd - zad*ycd
            yt = zad*xcd - xad*zcd
            zt = xad*ycd - yad*xcd
            rt2 = xt*xt + yt*yt + zt*zt
            delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
            xpd = xbd + xt*delta
            ypd = ybd + yt*delta
            zpd = zbd + zt*delta
            rbd2 = xbd*xbd + ybd*ybd + zbd*zbd
            rpd2 = xpd*xpd + ypd*ypd + zpd*zpd
            if (rbd2.ne.0.0d0 .and. rpd2.ne.0.0d0) then
               dot = xbd*xpd + ybd*ypd + zbd*zpd
               cosine = dot / sqrt(rbd2*rpd2)
               cosine = min (1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     find the out-of-plane angle bending energy
c
               dt = angle
               dt2 = dt * dt
               dt3 = dt2 * dt
               dt4 = dt2 * dt2
               e = opbunit * force * dt2
     &                * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total bond angle bending energy
c
               eopb = eopb + e
            end if
         end if
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eopbend1  --  out-of-plane energy and derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eopbend1" computes the out-of-plane bend potential energy
c     and first derivatives at trigonal centers as per Allinger's
c     MM2 and MM3 forcefields
c
c
      subroutine eopbend1
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'bath.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'opbend.i'
      include 'usage.i'
      include 'virial.i'
      integer i,iopbend
      integer ia,ib,ic,id
      real*8 e,force,angle
      real*8 dot,cosine,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 deddt,term
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd,rbd2
      real*8 xcd,ycd,zcd
      real*8 xpd,ypd,zpd,rpd2
      real*8 xt,yt,zt,rt2,ptrt2
      real*8 xm,ym,zm,rm,delta,delta2
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxip,dedyip,dedzip
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xdb,ydb,zdb
      logical proceed
c
c
c     zero out out-of-plane energy and first derivatives
c
      eopb = 0.0d0
      do i = 1, n
         deopb(1,i) = 0.0d0
         deopb(2,i) = 0.0d0
         deopb(3,i) = 0.0d0
      end do
c
c     calculate the out-of-plane bend energy and derivatives
c
      do iopbend = 1, nopbend
         i = iopb(iopbend)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         force = kopb(iopbend)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
c
c     compute the out-of-plane bending angle
c
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xbd = xib - xid
            ybd = yib - yid
            zbd = zib - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            xt = yad*zcd - zad*ycd
            yt = zad*xcd - xad*zcd
            zt = xad*ycd - yad*xcd
            rt2 = xt*xt + yt*yt + zt*zt
            delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
            xpd = xbd + xt*delta
            ypd = ybd + yt*delta
            zpd = zbd + zt*delta
            rbd2 = xbd*xbd + ybd*ybd + zbd*zbd
            rpd2 = xpd*xpd + ypd*ypd + zpd*zpd
            xm = ypd*zbd - zpd*ybd
            ym = zpd*xbd - xpd*zbd
            zm = xpd*ybd - ypd*xbd
            rm = sqrt(xm*xm + ym*ym + zm*zm)
            if (rm .ne. 0.0d0) then
               dot = xbd*xpd + ybd*ypd + zbd*zpd
               cosine = dot / sqrt(rbd2*rpd2)
               cosine = min (1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     find the out-of-plane energy and master chain rule terms
c
               dt = angle
               dt2 = dt * dt
               dt3 = dt2 * dt
               dt4 = dt2 * dt2
               e = opbunit * force * dt2
     &                * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
               deddt = opbunit * force * dt * radian
     &                    * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2
     &                        + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  deddt = deddt * fgrp
               end if
c
c     chain rule terms for central atom and its projection
c
               term = -deddt / (rbd2*rm)
               dedxib = term * (ybd*zm-zbd*ym)
               dedyib = term * (zbd*xm-xbd*zm)
               dedzib = term * (xbd*ym-ybd*xm)
               term = deddt / (rpd2*rm)
               dedxip = term * (ypd*zm-zpd*ym)
               dedyip = term * (zpd*xm-xpd*zm)
               dedzip = term * (xpd*ym-ypd*xm)
c
c     chain rule terms for peripheral plane-defining atoms
c
               delta2 = 2.0d0 * delta
               ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2
               term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd)
               dedxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2
               term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd)
               dedyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2
               term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd)
               dedzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2
               term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad)
               dedxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2
               term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad)
               dedyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2
               term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad)
               dedzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2
c
c     get out-of-plane atom chain rule terms by difference
c
               dedxid = -dedxia - dedxib - dedxic
               dedyid = -dedyia - dedyib - dedyic
               dedzid = -dedzia - dedzib - dedzic
c
c     increment the out-of-plane bend energy and gradient
c
               eopb = eopb + e
               deopb(1,ia) = deopb(1,ia) + dedxia
               deopb(2,ia) = deopb(2,ia) + dedyia
               deopb(3,ia) = deopb(3,ia) + dedzia
               deopb(1,ib) = deopb(1,ib) + dedxib
               deopb(2,ib) = deopb(2,ib) + dedyib
               deopb(3,ib) = deopb(3,ib) + dedzib
               deopb(1,ic) = deopb(1,ic) + dedxic
               deopb(2,ic) = deopb(2,ic) + dedyic
               deopb(3,ic) = deopb(3,ic) + dedzic
               deopb(1,id) = deopb(1,id) + dedxid
               deopb(2,id) = deopb(2,id) + dedyid
               deopb(3,id) = deopb(3,id) + dedzid
c
c     increment the virial for use in pressure computation
c
               if (isobaric) then
                  xab = xia - xib
                  yab = yia - yib
                  zab = zia - zib
                  xcb = xic - xib
                  ycb = yic - yib
                  zcb = zic - zib
                  xdb = xid - xib
                  ydb = yid - yib
                  zdb = zid - zib
                  virx = virx + xab*dedxia + xcb*dedxic + xdb*dedxid
                  viry = viry + yab*dedyia + ycb*dedyic + ydb*dedyid
                  virz = virz + zab*dedzia + zcb*dedzic + zdb*dedzid
               end if
            end if
         end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eopbend2  --  out-of-plane bend Hessian; numer  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eopbend2" calculates second derivatives of the out-of-plane
c     bend energy as per Allinger's MM2 and MM3 forcefields for a
c     single atom using finite difference methods
c
c
      subroutine eopbend2 (i)
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'group.i'
      include 'hessn.i'
      include 'opbend.i'
      integer i,j,k,iopbend
      integer ia,ib,ic,id
      real*8 eps,fgrp,old,term
      real*8 d0(3,maxatm)
      logical proceed
c
c
c     set stepsize for derivatives and default group weight
c
      eps = 1.0d-7
      fgrp = 1.0d0
c
c     compute numerical out-of-plane Hessian for current atom
c
      do iopbend = 1, nopbend
         k = iopb(iopbend)
         ia = iang(1,k)
         ib = iang(2,k)
         ic = iang(3,k)
         id = iang(4,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (i.eq.ia .or. i.eq.ib .or.
     &                              i.eq.ic .or. i.eq.id)
c
c     find first derivatives for the base structure
c
         if (proceed) then
            term = fgrp / eps
            call eopbend2b (iopbend)
            do j = 1, 3
               d0(j,ia) = deopb(j,ia)
               d0(j,ib) = deopb(j,ib)
               d0(j,ic) = deopb(j,ic)
               d0(j,id) = deopb(j,id)
            end do
c
c     find numerical x-components via perturbed structures
c
            old = x(i)
            x(i) = x(i) + eps
            call eopbend2b (iopbend)
            x(i) = old
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + term*(deopb(j,ia)-d0(j,ia))
               hessx(j,ib) = hessx(j,ib) + term*(deopb(j,ib)-d0(j,ib))
               hessx(j,ic) = hessx(j,ic) + term*(deopb(j,ic)-d0(j,ic))
               hessx(j,id) = hessx(j,id) + term*(deopb(j,id)-d0(j,id))
            end do
c
c     find numerical y-components via perturbed structures
c
            old = y(i)
            y(i) = y(i) + eps
            call eopbend2b (iopbend)
            y(i) = old
            do j = 1, 3
               hessy(j,ia) = hessy(j,ia) + term*(deopb(j,ia)-d0(j,ia))
               hessy(j,ib) = hessy(j,ib) + term*(deopb(j,ib)-d0(j,ib))
               hessy(j,ic) = hessy(j,ic) + term*(deopb(j,ic)-d0(j,ic))
               hessy(j,id) = hessy(j,id) + term*(deopb(j,id)-d0(j,id))
            end do
c
c     find numerical z-components via perturbed structures
c
            old = z(i)
            z(i) = z(i) + eps
            call eopbend2b (iopbend)
            z(i) = old
            do j = 1, 3
               hessz(j,ia) = hessz(j,ia) + term*(deopb(j,ia)-d0(j,ia))
               hessz(j,ib) = hessz(j,ib) + term*(deopb(j,ib)-d0(j,ib))
               hessz(j,ic) = hessz(j,ic) + term*(deopb(j,ic)-d0(j,ic))
               hessz(j,id) = hessz(j,id) + term*(deopb(j,id)-d0(j,id))
            end do
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eopbend2b  --  out-of-plane bend derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eopbend2b" calculates out-of-plane bending first derivatives
c     for a single angle with respect to Cartesian coordinates;
c     used in computation of finite difference second derivatives
c
c
      subroutine eopbend2b (i)
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'math.i'
      include 'opbend.i'
      integer i,k,ia,ib,ic,id
      real*8 angle,force
      real*8 dot,cosine
      real*8 dt,dt2,dt3,dt4
      real*8 deddt,term
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd,rbd2
      real*8 xcd,ycd,zcd
      real*8 xpd,ypd,zpd,rpd2
      real*8 xt,yt,zt,rt2,ptrt2
      real*8 xm,ym,zm,rm,delta,delta2
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxip,dedyip,dedzip
c
c
c     set the atom numbers and parameters for this angle
c
      k = iopb(i)
      ia = iang(1,k)
      ib = iang(2,k)
      ic = iang(3,k)
      id = iang(4,k)
      force = kopb(i)
c
c     get the coordinates of the atoms in the angle
c
      xia = x(ia)
      yia = y(ia)
      zia = z(ia)
      xib = x(ib)
      yib = y(ib)
      zib = z(ib)
      xic = x(ic)
      yic = y(ic)
      zic = z(ic)
      xid = x(id)
      yid = y(id)
      zid = z(id)
c
c     zero out the first derivative components
c
      deopb(1,ia) = 0.0d0
      deopb(2,ia) = 0.0d0
      deopb(3,ia) = 0.0d0
      deopb(1,ib) = 0.0d0
      deopb(2,ib) = 0.0d0
      deopb(3,ib) = 0.0d0
      deopb(1,ic) = 0.0d0
      deopb(2,ic) = 0.0d0
      deopb(3,ic) = 0.0d0
      deopb(1,id) = 0.0d0
      deopb(2,id) = 0.0d0
      deopb(3,id) = 0.0d0
c
c     compute the out-of-plane bending angle
c
      xad = xia - xid
      yad = yia - yid
      zad = zia - zid
      xbd = xib - xid
      ybd = yib - yid
      zbd = zib - zid
      xcd = xic - xid
      ycd = yic - yid
      zcd = zic - zid
      xt = yad*zcd - zad*ycd
      yt = zad*xcd - xad*zcd
      zt = xad*ycd - yad*xcd
      rt2 = xt*xt + yt*yt + zt*zt
      delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
      xpd = xbd + xt*delta
      ypd = ybd + yt*delta
      zpd = zbd + zt*delta
      rbd2 = xbd*xbd + ybd*ybd + zbd*zbd
      rpd2 = xpd*xpd + ypd*ypd + zpd*zpd
      xm = ypd*zbd - zpd*ybd
      ym = zpd*xbd - xpd*zbd
      zm = xpd*ybd - ypd*xbd
      rm = sqrt(xm*xm + ym*ym + zm*zm)
      if (rm .ne. 0.0d0) then
         dot = xbd*xpd + ybd*ypd + zbd*zpd
         cosine = dot / sqrt(rbd2*rpd2)
         cosine = min (1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
c
c     get the out-of-plane bending master chain rule terms
c
         dt = angle
         dt2 = dt * dt
         dt3 = dt2 * dt
         dt4 = dt2 * dt2
         deddt = opbunit * force * dt * radian
     &              * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2
     &                  + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)
c
c     chain rule terms for central atom and its projection
c
         term = -deddt / (rbd2*rm)
         dedxib = term * (ybd*zm-zbd*ym)
         dedyib = term * (zbd*xm-xbd*zm)
         dedzib = term * (xbd*ym-ybd*xm)
         term = deddt / (rpd2*rm)
         dedxip = term * (ypd*zm-zpd*ym)
         dedyip = term * (zpd*xm-xpd*zm)
         dedzip = term * (xpd*ym-ypd*xm)
c
c     chain rule terms for peripheral plane-defining atoms
c
         delta2 = 2.0d0 * delta
         ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2
         term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd)
         dedxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2
         term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd)
         dedyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2
         term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd)
         dedzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2
         term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad)
         dedxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2
         term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad)
         dedyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2
         term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad)
         dedzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2
c
c     get out-of-plane atom chain rule terms by difference
c
         dedxid = -dedxia - dedxib - dedxic
         dedyid = -dedyia - dedyib - dedyic
         dedzid = -dedzia - dedzib - dedzic
c
c     set the out-of-plane bending derivatives
c
         deopb(1,ia) = dedxia
         deopb(2,ia) = dedyia
         deopb(3,ia) = dedzia
         deopb(1,ib) = dedxib
         deopb(2,ib) = dedyib
         deopb(3,ib) = dedzib
         deopb(1,ic) = dedxic
         deopb(2,ic) = dedyic
         deopb(3,ic) = dedzic
         deopb(1,id) = dedxid
         deopb(2,id) = dedyid
         deopb(3,id) = dedzid
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
c     ##  subroutine eopbend3  --  out-of-plane bending & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eopbend3" computes the out-of-plane bend potential energy at
c     trigonal centers as per Allinger's MM2 and MM3 forcefields;
c     also partitions the energy among the atoms
c
c
      subroutine eopbend3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'opbend.i'
      include 'usage.i'
      integer i,iopbend
      integer ia,ib,ic,id
      real*8 e,angle,force
      real*8 dot,cosine,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd,rbd2
      real*8 xcd,ycd,zcd
      real*8 xpd,ypd,zpd,rpd2
      real*8 xt,yt,zt,rt2,delta
      logical header,huge,proceed
c
c
c     zero out the out-of-plane bend energy and partitioning
c
      neopb = 0
      eopb = 0.0d0
      do i = 1, n
         aeopb(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the out-of-plane bending energy term
c
      do iopbend = 1, nopbend
         i = iopb(iopbend)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         force = kopb(iopbend)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
c
c     compute the out-of-plane bending angle
c
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xbd = xib - xid
            ybd = yib - yid
            zbd = zib - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            xt = yad*zcd - zad*ycd
            yt = zad*xcd - xad*zcd
            zt = xad*ycd - yad*xcd
            rt2 = xt*xt + yt*yt + zt*zt
            delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
            xpd = xbd + xt*delta
            ypd = ybd + yt*delta
            zpd = zbd + zt*delta
            rbd2 = xbd*xbd + ybd*ybd + zbd*zbd
            rpd2 = xpd*xpd + ypd*ypd + zpd*zpd
            if (rbd2.ne.0.0d0 .and. rpd2.ne.0.0d0) then
               dot = xbd*xpd + ybd*ypd + zbd*zpd
               cosine = dot / sqrt(rbd2*rpd2)
               cosine = min (1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     find the out-of-plane angle bending energy
c
               dt = angle
               dt2 = dt * dt
               dt3 = dt2 * dt
               dt4 = dt2 * dt2
               e = opbunit * force * dt2
     &                * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total bond angle bending energy
c
               neopb = neopb + 1
               eopb = eopb + e
               aeopb(ib) = aeopb(ib) + e
c
c     print a warning if the energy of this angle is large
c
               huge = (e .gt. 2.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Out-of-Plane Bending',
     &                          ' Interactions :',
     &                       //,' Type',11x,'Atom Names',20x,'Ideal',
     &                          4x,'Actual',6x,'Energy',/)
                  end if
                  write (iout,20)  ib,name(ib),id,name(id),0.0d0,angle,e
   20             format (' O-P-Bend ',i5,'-',a3,1x,i5,'-',a3,
     &                       12x,2f10.4,f12.4)
               end if
            end if
         end if
      end do
      return
      end
