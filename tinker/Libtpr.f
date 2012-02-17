C 21 Apr 10 - NA - format changes in printing 
C 12 Dec 08 - MWS - PRTXYZ: drop headers if writing to trajectory file
C 28 Feb 07 - YBG - PRTXYZ: change format to stay within 80 cols
C 15 Nov 06 - YBG - READPRM: check input vs. dimension limits
C  9 MAR 00 - CHC - FIX for parallel run
C 14 OCT 98 - CHC change title -> ttitle
C                 readxyz : modified for parallel run
c  8 MAY 98 - JRS readxyz: redirect input from GAMESS input
c                 prtxyz:  redirect ouput to GAMESS output
c                 rotmat: renamed romatt
c
c DGF thinks that maxres usage is very dangerous because
c no check is done to ensure there is enough elements allocated. 
c However many subroutines using it are not envoked. 
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine pialter  --  modify constants for pisystem  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "pialter" first modifies bond lengths and force constants
c     according to the standard bond slope parameters and the
c     bond order values stored in "pnpl"; also alters some 2-fold
c     torsional parameters based on the bond-order * beta matrix
c
c
      subroutine pialter
      implicit none
      include 'sizes.i'
      include 'bond.i'
      include 'inform.i'
      include 'iounit.i'
      include 'pibond.i'
      include 'pistuf.i'
      include 'piterm.i'
      include 'tors.i'
      integer i,j,k
c
c
c     modify the stretching constants and natural bond lengths
c
      if (debug) then
         write (iout,10)
   10    format ()
      end if
      do i = 1, npibond
         k = pibond(1,i)
         bk(k) = bkpi(i) - kslope(i) * (1.0d0-pnpl(i))
         bl(k) = blpi(i) + lslope(i) * (1.0d0-pnpl(i))
         if (debug) then
            write (iout,20)  ibnd(1,k),ibnd(2,k),bkpi(i),
     &                       blpi(i),bk(k),bl(k)
   20       format (' PIALTER  --  Bond',3x,2i6,5x,
     &                f7.3,f8.4,' -->',f7.3,f8.4)
         end if
      end do
c
c     modify the 2-fold torsional constants across pibonds
c
      if (debug) then
         write (iout,30)
   30    format ()
      end if
      do i = 1, npitors
         j = pitors(1,i)
         k = pitors(2,i)
         tors2(1,j) = pbpl(k) * torspi(i)
         if (debug) then
            write (iout,40)  (itors(k,j),k=1,4),torspi(i),tors2(1,j)
   40       format (' PIALTER  --  Torsion',4i6,3x,f10.3,' -->',f10.3)
         end if
      end do
      if (debug) then
         write (iout,50)
   50    format ()
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
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine pimove  --  translate/rotate bond vector  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "pimove" rotates the vector between atoms "list(1)" and
c     "list(2)" so that atom 1 is at the origin and atom 2 along
c     the x-axis; the atoms defining the respective planes are
c     also moved and their bond lengths normalized
c
c
      subroutine pimove (list,xr,yr,zr)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,j,list(8)
      real*8 xr(8),yr(8),zr(8),xt,yt,zt
      real*8 denom,sine,cosine,xold
c
c
c     translate "list" atoms to place atom 1 at origin
c
      j = list(1)
      xt = x(j)
      yt = y(j)
      zt = z(j)
      do i = 1, 8
         j = list(i)
         xr(i) = x(j) - xt
         yr(i) = y(j) - yt
         zr(i) = z(j) - zt
      end do
c
c     rotate "list" atoms to place atom 2 on the x-axis
c
      denom = sqrt(xr(2)**2 + yr(2)**2)
      if (denom .ne. 0.0d0) then
         sine = yr(2) / denom
         cosine = xr(2) / denom
         do i = 1, 8
            xold = xr(i)
            xr(i) = xr(i)*cosine + yr(i)*sine
            yr(i) = yr(i)*cosine - xold*sine
         end do
      end if
      denom = sqrt(xr(2)**2 + zr(2)**2)
      if (denom .ne. 0.0d0) then
         sine = zr(2) / denom
         cosine = xr(2) / denom
         do i = 1, 8
            xold = xr(i)
            xr(i) = xr(i)*cosine + zr(i)*sine
            zr(i) = zr(i)*cosine - xold*sine
         end do
      end if
c
c     normalize the coordinates of atoms defining the plane
c     for atom 1 (ie, make all these atoms have unit length to
c     atom 1) so that the orbital makes equal angles with the
c     atoms rather than simply being perpendicular to the common
c     plane of the atoms
c
      do i = 3, 5
         if (list(i) .ne. list(1)) then
            denom = sqrt(xr(i)**2+yr(i)**2+zr(i)**2)
            xr(i) = xr(i) / denom
            yr(i) = yr(i) / denom
            zr(i) = zr(i) / denom
         end if
      end do
c
c     normalization of plane defining atoms for atom 2; for the
c     x-coordinate we translate back to the origin, normalize
c     and then retranslate back along the x-axis
c
      do i = 6, 8
         if (list(i) .ne. list(2)) then
            denom = sqrt((xr(i)-xr(2))**2+yr(i)**2+zr(i)**2)
            xr(i) = (xr(i)-xr(2))/denom + xr(2)
            yr(i) = yr(i) / denom
            zr(i) = zr(i) / denom
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine piplane  --  plane perpendicular to orbital  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "piplane" selects the three atoms which specify the plane
c     perpendicular to each p-orbital; the current version will
c     fail in certain situations, including ketenes, allenes,
c     and isolated or adjacent triple bonds
c
c
      subroutine piplane
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'pistuf.i'
      integer i,j,iorb,atmnum,trial
      integer alpha,beta,gamma,attach
      logical done
c
c
c     for each pisystem atom, find a set of atoms which define
c     the p-orbital's plane based on piatom's atomic number and
c     the number and type of attached atoms
c
      do i = 1, norbit
         iorb = iorbit(i)
         attach = n12(iorb)
         atmnum = atomic(iorb)
         done = .false.
c
c     most common case of an atom bonded to three atoms
c
         if (attach .eq. 3) then
            piperp(1,i) = i12(1,iorb)
            piperp(2,i) = i12(2,iorb)
            piperp(3,i) = i12(3,iorb)
            done = .true.
c
c     any non-alkyne atom bonded to exactly two atoms
c
         else if (attach.eq.2 .and. atmnum.ne.6) then
            piperp(1,i) = iorb
            piperp(2,i) = i12(1,iorb)
            piperp(3,i) = i12(2,iorb)
            done = .true.
c
c     atom bonded to four different atoms (usually two lone
c     pairs and two "real" atoms); use the "real" atoms
c
         else if (attach .eq. 4) then
            piperp(1,i) = iorb
            do j = 1, n12(iorb)
               trial = i12(j,iorb)
               if (atomic(trial) .ne. 0) then
                  if (piperp(2,i) .eq. 0) then
                     piperp(2,i) = trial
                  else
                     piperp(3,i) = trial
                     done = .true.
                  end if
               end if
            end do
c
c     "carbonyl"-type oxygen atom, third atom is any atom
c     attached to the "carbonyl carbon"; fails for ketenes
c
         else if (attach.eq.1 .and. atmnum.eq.8) then
            alpha = i12(1,iorb)
            beta = i12(1,alpha)
            if (beta .eq. iorb)  beta = i12(2,alpha)
            piperp(1,i) = iorb
            piperp(2,i) = alpha
            piperp(3,i) = beta
            done = .true.
c
c     an sp nitrogen atom, third atom must be a gamma atom
c
         else if (attach.eq.1 .and. atmnum.eq.7) then
            alpha = i12(1,iorb)
            do j = 1, n12(alpha)
               trial = i12(j,alpha)
               if (trial.ne.iorb .and. listpi(trial) .and.
     &             n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            gamma = i12(1,beta)
            if (gamma .eq. alpha)  gamma = i12(2,beta)
            piperp(3,i) = iorb
            piperp(3,i) = alpha
            piperp(3,i) = gamma
c
c     an sp carbon atom; third atom must be an atom attached
c     to the non-sp piatom bonded to the original carbon
c
         else if (attach.eq.2 .and. atmnum.eq.6) then
            alpha = i12(1,iorb)
            if ((n12(alpha).eq.2 .and. atomic(alpha).eq.6) .or.
     &          (n12(alpha).eq.1 .and. atomic(alpha).eq.7))
     &         alpha = i12(2,iorb)
            do j = 1, n12(iorb)
               trial = i12(j,iorb)
               if (trial.ne.iorb .and. trial.ne.alpha .and.
     &             listpi(trial) .and. n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            do j = 1, n12(alpha)
               trial = i12(j,alpha)
               if (trial.ne.iorb .and. trial.ne.alpha .and.
     &             listpi(trial) .and. n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            gamma = i12(1,beta)
            if (gamma.eq.iorb .or. gamma.eq.alpha)  gamma = i12(2,beta)
            piperp(1,i) = iorb
            piperp(2,i) = alpha
            piperp(3,i) = gamma
         end if
c
c     quit if the p-orbital plane remains undefined
c
         if (.not. done) then
            write (iout,10)  iorb
   10       format(/,' PIPLANE  --  Failure to Define a',
     &                ' p-Orbital Plane for Atom',i6)
            call fatal
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine piscf  --  scf molecular orbital calculation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "piscf" performs an scf molecular orbital calculation for
c     the pisystem using a modified Pariser-Parr-Pople method
c
c
      subroutine piscf
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'pibond.i'
      include 'picalc.i'
      include 'pistuf.i'
      include 'units.i'
      integer i,j,k,m,iter,maxiter
      integer iatn,jatn,iorb,jorb
      real*8 f(maxpi,maxpi),hc(maxpi,maxpi)
      real*8 v(maxpi,maxpi),en(maxpi),ip(maxpi)
      real*8 gamma(maxpi,maxpi),ed(maxpi,maxpi)
      real*8 work1(maxpi),work2(maxpi)
      real*8 delta,converge,p,xij,yij,zij
      real*8 hcii,gii,gij,g11,g11sq,g12,g14
      real*8 povlap(maxpib),rij,erij,brij
      real*8 ovlap,covlap,ioniz_i,ioniz_j,cionize
      real*8 rijsq,hcij,qi,total,totold,ebeta
      real*8 aeth,abnz,ebe,ebb,ble,blb,eebond,bebond
      real*8 s1,s2,gjk,vij,vik,vmj,vmk,xi,xj,xk,xg
      character*6 mode
      logical initial
      save f,initial
      data initial  / .true. /
c
c
c     only needs to be done if pisystem is present
c
      if (norbit .eq. 0)  return
c
c     initialize some constants and parameters:
c
c     mode      planar or nonplanar pi-calculation
c     maxiter   maximum number of scf iterations
c     converge  criterion for scf convergence
c     ebeta     value of resonance integral for ethylene
c     cionize   ionizaton potential of carbon (hartree)
c
      mode = 'Planar'
      maxiter = 50
      converge = 0.00000001d0
      ebeta = -0.0757d0
      cionize = -11.16d0 / evolt
c
c     set the bond energies, alpha values and ideal bond length
c     parameter for carbon-carbon pibond type parameters:
c
c     ebe = equilibrium bond energy in ethylene
c     ebb = equilibrium bond energy in benzene
c     aeth = the P-P-P constant "a" in ethylene
c     abnz = the P-P-P constant "a" in benzene
c     ble = equilibrium bond length in ethylene
c     blb = equilibrium bond length in benzene
c
      ebe = 129.37d0
      ebb = 117.58d0
      aeth = 2.309d0
      abnz = 2.142d0
      ble = 1.338d0
      blb = 1.397d0
c
c     assign empirical one-center Coulomb integrals, and
c     first or second ionization potential depending on
c     whether the orbital contributes one or two electrons
c
      do i = 1, norbit
         gamma(i,i) = em(i)
         ip(i) = w(i) + (1.0d0-q(i))*em(i)
      end do
c
c     calculate two-center repulsion integrals
c     according to Ohno's semi-empirical formula
c
      do i = 1, norbit-1
         iorb = iorbit(i)
         gii = gamma(i,i)
         do j = i+1, norbit
            jorb = iorbit(j)
            g11 = 0.5d0 * (gii+gamma(j,j))
            g11sq = 1.0d0 / g11**2
            xij = x(iorb) - x(jorb)
            yij = y(iorb) - y(jorb)
            zij = z(iorb) - z(jorb)
            rijsq = (xij**2 + yij**2 + zij**2) / bohr**2
            g12 = 1.0d0 / sqrt(rijsq+g11sq)
            gamma(i,j) = g12
            gamma(j,i) = g12
         end do
      end do
c
c     zero out the resonance integral values
c
      do i = 1, norbit
         do j = 1, norbit
            hc(j,i) = 0.0d0
         end do
      end do
c
c     the first term in the sum to find alpha is the first
c     or second ionization potential, then the two-center
c     repulsion integrals are added
c
      do i = 1, norbit
         hcii = ip(i)
         do j = 1, norbit
            if (i .ne. j) then
               hcii = hcii - q(j)*gamma(i,j)
            end if
         end do
         hc(i,i) = hcii
      end do
c
c     get two-center repulsion integrals via Ohno's formula
c
      do k = 1, npibond
         i = pibond(2,k)
         j = pibond(3,k)
         iorb = iorbit(i)
         jorb = iorbit(j)
         iatn = atomic(iorb)
         jatn = atomic(jorb)
         xij = x(iorb) - x(jorb)
         yij = y(iorb) - y(jorb)
         zij = z(iorb) - z(jorb)
         rij = sqrt(xij**2 + yij**2 + zij**2)
         rijsq = rij**2 / bohr**2
         g11 = 0.5d0 * (gamma(i,i)+gamma(j,j))
         g11sq = 1.0d0 / g11**2
         g12 = gamma(i,j)
c
c     compute the bond energy using a Morse potential
c
         erij = aeth * (ble-rij)
         brij = abnz * (blb-rij)
         eebond = (2.0d0*exp(erij)-exp(2.0d0*erij)) * ebe / hartree
         bebond = (2.0d0*exp(brij)-exp(2.0d0*brij)) * ebb / hartree
c
c     compute the carbon-carbon resonance integral using
c     the Whitehead and Lo formula
c
         g14 = 1.0d0 / sqrt(4.0d0*rijsq+g11sq)
         hcij = 1.5d0*(bebond-eebond) - 0.375d0*g11
     &             + (5.0d0/12.0d0)*g12 - g14/24.0d0
c
c     if either atom is non-carbon, then factor the resonance
c     integral by overlap ratio and ionization potential ratio
c
         if (iatn.ne.6 .or. jatn.ne.6) then
            call overlap (iatn,jatn,rij,ovlap)
            call overlap (6,6,rij,covlap)
            hcij = hcij * (ovlap/covlap)
            ioniz_i = ip(i)
            if (q(i) .ne. 1.0d0) then
               if (iatn .eq. 7)  ioniz_i = 0.595d0 * ioniz_i
               if (iatn .eq. 8)  ioniz_i = 0.525d0 * ioniz_i
               if (iatn .eq. 16)  ioniz_i = 0.89d0 * ioniz_i
            end if
            ioniz_j = ip(j)
            if (q(j) .ne. 1.0d0) then
               if (jatn .eq. 7)  ioniz_j = 0.595d0 * ioniz_j
               if (jatn .eq. 8)  ioniz_j = 0.525d0 * ioniz_j
               if (jatn .eq. 16)  ioniz_j = 0.89d0 * ioniz_j
            end if
            hcij = hcij * (ioniz_i+ioniz_j)/(2.0d0*cionize)
         end if
c
c     set symmetric elements to the same value
c
         hc(i,j) = hcij
         hc(j,i) = hcij
      end do
c
c     make an initial guess at the Fock matrix if needed
c
      if (initial) then
         initial = .false.
         do i = 1, norbit
            do j = 1, norbit
               f(j,i) = hc(j,i)
            end do
         end do
         do i = 1, norbit
            f(i,i) = 0.5d0 * ip(i)
         end do
      end if
c
c     now, do the scf-mo computation; note that it needs to
c     be done twice, initially for the planar analog of the
c     actual system; then for the nonplanar (actual) system
c
      dowhile (mode.eq.'Planar' .or. mode.eq.'Nonpln')
         if (mode .eq. 'Nonpln') then
            call pitilt (povlap)
            do k = 1, npibond
               i = pibond(2,k)
               j = pibond(3,k)
               hc(i,j) = hc(i,j) * povlap(k)
               hc(j,i) = hc(i,j)
            end do
         end if
c
c     perform scf iterations until convergence is reached;
c     diagonalize the Fock matrix "f" to get the mo's,
c     then use mo's to form the next "f" matrix assuming
c     zero differential overlap except for one-center
c     exchange repulsions
c
         iter = 0
         delta = 2.0d0 * converge
         dowhile (delta.gt.converge .and. iter.lt.maxiter)
            iter = iter + 1
            call tnk_jacobi (norbit,maxpi,f,en,v,work1,work2)
            do i = 1, norbit
               do j = i, norbit
                  s1 = 0.0d0
                  s2 = 0.0d0
                  gij = gamma(i,j)
                  do k = 1, nfill
                     s2 = s2 - v(i,k)*v(j,k)*gij
                     if (i .eq. j) then
                        do m = 1, norbit
                           s1 = s1 + 2.0d0*gamma(i,m)*v(m,k)**2
                        end do
                     end if
                  end do
                  f(i,j) =  s1 + s2 + hc(i,j)
                  f(j,i) = f(i,j)
               end do
            end do
c
c     calculate the ground state energy, where "xi" sums the
c     molecular core integrals, "xj" sums the molecular coulomb
c     repulsion integrals, "xk" sums the molecular exchange
c     repulsion integrals, and "xg" is sums the nuclear repulsion
c
            xi = 0.0d0
            xj = 0.0d0
            xk = 0.0d0
            xg = 0.0d0
            do i = 1, nfill
               do j = 1, norbit
                  vij = v(j,i)
                  do k = 1, norbit
                     vik = v(k,i)
                     gjk = gamma(j,k)
                     xi = xi + 2.0d0*vij*vik*hc(j,k)
                     do m = 1, nfill
                        vmj = v(j,m)
                        vmk = v(k,m)
                        xj = xj + 2.0d0*vij*vij*vmk*vmk*gjk
                        xk = xk - vij*vmj*vik*vmk*gjk
                     end do
                  end do
               end do
            end do
            do i = 1, norbit-1
               qi = q(i)
               do j = i+1, norbit
                  xg = xg + qi*q(j)*gamma(i,j)
               end do
            end do
            total = xi + xj + xk + xg
            if (iter .ne. 1)  delta = abs(total-totold)
            totold = total
         end do
c
c     print warning if scf-mo iteration did not converge
c
         if (delta .gt. converge) then
            write (iout,10)
   10       format (' PISCF  --  The SCF-MO Iteration',
     &              ' has not reached Self-Consistency')
         end if
c
c     calculate electron densities from filled mo's
c
         do i = 1, norbit
            do j = 1, norbit
               ed(i,j) = 0.0d0
               do k = 1, nfill
                  ed(i,j) = ed(i,j) + 2.0d0*v(i,k)*v(j,k)
               end do
            end do
         end do
c
c     print out results for the scf computation
c
         if (debug) then
            write (iout,20)  mode,iter,total,delta
   20       format (/,1x,a6,' SCF-MO Iteration',i3,
     &              //,' Energy =',f10.4,3x,'Delta =',d12.4)
            write (iout,30)  xi,xj,xk,xg
   30       format (/,' Core Integrals      =',f11.4,
     &              /,' Coulomb Repulsion   =',f11.4,
     &              /,' Exchange Repulsion  =',f11.4,
     &              /,' Nuclear Repulsion   =',f11.4)
            write (iout,40)  (en(i),i=1,norbit)
   40       format (/,' Orbital Energies',/,(8f9.4))
            write (iout,50)
   50       format (/,' Molecular Orbitals')
            do i = 1, norbit
               write (iout,60)  (v(i,j),j=1,norbit)
   60          format (8f9.4)
            end do
            write (iout,70)
   70       format (/,' Fock Matrix')
            do i = 1, norbit
               write (iout,60)  (f(i,j),j=1,norbit)
            end do
            write (iout,80)  (ed(i,i),i=1,norbit)
   80       format (/,' Electron Densities',/,(8f9.4))
            write (iout,90)
   90       format (/,' Density Matrix')
            do i = 1, norbit
               write (iout,60)  (ed(i,j),j=1,norbit)
            end do
            write (iout,100)
  100       format (/,' H-core Matrix')
            do i = 1, norbit
               write (iout,60)  (hc(i,j),j=1,norbit)
            end do
            write (iout,110)
  110       format (/,' Gamma Matrix')
            do i = 1, norbit
               write (iout,60)  (gamma(i,j),j=1,norbit)
            end do
         end if
c
c     now, get the bond orders (compute p and p*b)
c
         if (debug) then
            write (iout,120)
  120       format (/,' PI Bond Orders')
         end if
         do k = 1, npibond
            i = pibond(2,k)
            j = pibond(3,k)
            p = 0.0d0
            do m = 1, nfill
               p = p + 2.0d0*v(i,m)*v(j,m)
            end do
            if (mode .eq. 'Planar') then
               pbpl(k) = p * hc(i,j)/ebeta
            else if (mode .eq. 'Nonpln') then
               pnpl(k) = p
            end if
            if (debug) then
               i = ibnd(1,pibond(1,k))
               j = ibnd(2,pibond(1,k))
               write (iout,130)  i,j,p
  130          format (4x,2i5,3x,f10.4)
            end if
         end do
c
c     if we have done planar calculation, do the nonplanar;
c     when both are complete, alter the pisystem constants
c
         if (mode .eq. 'Planar') then
            mode = 'Nonpln'
         else if (mode .eq. 'Nonpln') then
            mode = '      '
         end if
      end do
c
c     alter torsional and bond constants for pisystem
c
      call pialter
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
c     ##  subroutine pitilt  --  direction cosines for pisystem  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "pitilt" calculates for each pibond the ratio of the
c     actual p-orbital overlap integral to the ideal overlap
c     if the same orbitals were perfectly parallel
c
c
      subroutine pitilt (povlap)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'pistuf.i'
      integer i,j,k,m,iorb,jorb,list(8)
      real*8 ideal,povlap(maxpib)
      real*8 xr(8),yr(8),zr(8),cosine
      real*8 xij,yij,zij,rij,rnorm
      real*8 a1,b1,c1,a2,b2,c2
      real*8 x2,y2,z2,x3,y3,z3
c
c
c     planes defining each p-orbital are in "piperp"; transform
c     coordinates of "iorb", "jorb" and their associated planes
c     to put "iorb" at origin and "jorb" along the x-axis
c
      do k = 1, npibond
         i = pibond(2,k)
         j = pibond(3,k)
         iorb = iorbit(i)
         jorb = iorbit(j)
         list(1) = iorb
         list(2) = jorb
         do m = 1, 3
            list(m+2) = piperp(m,i)
            list(m+5) = piperp(m,j)
         end do
         call pimove (list,xr,yr,zr)
c
c     check for sp-hybridized carbon in current bond;
c     assume perfect overlap for any such pibond
c
         if ((atomic(iorb).eq.6 .and. n12(iorb).eq.2) .or.
     &       (atomic(jorb).eq.6 .and. n12(jorb).eq.2)) then
            povlap(k) = 1.0d0
c
c     find and normalize a vector parallel to first p-orbital
c
         else
            x2 = xr(4) - xr(3)
            y2 = yr(4) - yr(3)
            z2 = zr(4) - zr(3)
            x3 = xr(5) - xr(3)
            y3 = yr(5) - yr(3)
            z3 = zr(5) - zr(3)
            a1 = y2*z3 - y3*z2
            b1 = x3*z2 - x2*z3
            c1 = x2*y3 - x3*y2
            rnorm = sqrt(a1*a1+b1*b1+c1*c1)
            a1 = a1 / rnorm
            b1 = b1 / rnorm
            c1 = c1 / rnorm
c
c     now find vector parallel to the second p-orbital,
c     "a2" changes sign to correspond to internuclear axis
c
            x2 = xr(7) - xr(6)
            y2 = yr(7) - yr(6)
            z2 = zr(7) - zr(6)
            x3 = xr(8) - xr(6)
            y3 = yr(8) - yr(6)
            z3 = zr(8) - zr(6)
            a2 = y2*z3 - y3*z2
            b2 = x3*z2 - x2*z3
            c2 = x2*y3 - x3*y2
            rnorm = sqrt(a2*a2+b2*b2+c2*c2)
            a2 = -a2 / rnorm
            b2 = b2 / rnorm
            c2 = c2 / rnorm
c
c     compute the cosine of the angle between p-orbitals;
c     if more than 90 degrees, reverse one of the vectors
c
            cosine = a1*a2 + b1*b2 + c1*c2
            if (cosine .lt. 0.0d0) then
               a2 = -a2
               b2 = -b2
               c2 = -c2
            end if
c
c     find overlap if the orbitals were perfectly parallel
c
            xij = x(iorb) - x(jorb)
            yij = y(iorb) - y(jorb)
            zij = z(iorb) - z(jorb)
            rij = sqrt(xij**2 + yij**2 + zij**2)
            call overlap (atomic(iorb),atomic(jorb),rij,ideal)
c
c     set ratio of actual to ideal overlap for current pibond
c
            povlap(k) = ideal*a1*a2 + b1*b2 + c1*c2
         end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine power  --  largest eigenvalues by power method  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "power" uses the power method with deflation to compute the
c     few largest eigenvalues and eigenvectors of a symmetric matrix
c
c     n     logical dimension of the matrix to be diagonalized
c     np    physical dimension of the matrix storage area
c     nv    number of largest eigenvalues to be extracted
c     a     input with the matrix to be diagonalized; only
c              the lower triangle and diagonal are required
c     ev    returned with the eigenvalues in descending order
c     vec   returned with the eigenvectors of the matrix
c     work  temporary work vector
c
c
      subroutine power (n,np,nv,a,ev,vec,work)
      implicit none
      include 'iounit.i'
      integer i,j,k,n,np,nv
      integer iter,maxiter
      real*8 random,eps
      real*8 dot1,dot2,ratio
      real*8 ev(nv),vec(np,nv)
      real*8 a(np,np),work(np)
c
c
c     initialize number of iterations and convergence criteria
c
      maxiter = 500
      eps = 1.0d-6
c
c     use identity vector as initial guess for eigenvectors
c
      do j = 1, nv
         do i = 1, n
            vec(i,j) = 1.0d0
         end do
      end do
c
c     find the few largest eigenvalues and eigenvectors
c
      do k = 1, nv
         ev(k) = 0.0d0
         dot1 = 0.0d0
         do i = 1, n
            work(i) = 0.0d0
            do j = 1, i-1
               work(i) = work(i) + a(i,j)*vec(j,k)
            end do
            do j = i, n
               work(i) = work(i) + a(j,i)*vec(j,k)
            end do
            dot1 = dot1 + work(i)**2
         end do
c
c     if in or near null space, use random guess as eigenvector
c
         if (dot1 .le. 100.0d0*eps*dble(n)) then
            do i = 1, n
               work(i) = random ()
            end do
         end if
c
c     find the current eigenvalue by iterating to convergence;
c     first multiply vector by matrix and compute dot products
c
         do iter = 1, maxiter
            dot1 = 0.0d0
            dot2 = 0.0d0
            do i = 1, n
               vec(i,k) = 0.0d0
               do j = 1, i-1
                  vec(i,k) = vec(i,k) + a(i,j)*work(j)
               end do
               do j = i, n
                  vec(i,k) = vec(i,k) + a(j,i)*work(j)
               end do
               dot1 = dot1 + vec(i,k)**2
               dot2 = dot2 + vec(i,k)*work(i)
            end do
c
c     normalize new eigenvector and substitute for old one
c
            ratio = abs((ev(k)-dot2) / dot2)
            ev(k) = dot2
            dot1 = sqrt(dot1)
            do i = 1, n
               vec(i,k) = vec(i,k) / dot1
               work(i) = vec(i,k)
            end do
            if (ratio .lt. eps)  goto 20
         end do
         write (iout,10)  k
   10    format (/,' POWER  --  Eigenvalue',i3,' not Fully Converged')
c
c     eliminate the current eigenvalue from the matrix
c
   20    continue
         do i = 1, n
            do j = i, n
               a(j,i) = a(j,i) - ev(k)*vec(i,k)*vec(j,k)
            end do
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
c     #########################################################
c     ##                                                     ##
c     ##  function precise  --  determine machine precision  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "precise" finds one of three machine precision values
c
c        (1) the smallest positive floating point value
c        (2) the smallest relative floating point spacing
c        (3) the largest relative floating point spacing
c
c
      function precise (i)
      implicit none
      integer i
      real*8 precise,value
      real*8 zero,one,delta
c
c
c     set zero, one and multiplicative factor
c
      zero = 0.0d0
      one = 1.0d0
      delta = 1.1d0
      precise = one
c
c     find the smallest positive floating point value;
c     minimum of 0.23x10-307 is a patch for SGI's
c
      if (i .eq. 1) then
c        dowhile (precise .ne. zero)
         dowhile (precise .ge. 0.23d-307)
            value = precise
            precise = precise / delta
         end do
         precise = value
c
c     find the smallest relative floating point spacing
c
      else if (i .eq. 2) then
         dowhile (one+precise .ne. one)
            value = precise
            precise = precise / delta
         end do
         precise = value
c
c     find the largest relative floating point spacing
c
      else if (i .eq. 3) then
         dowhile (one+precise .ne. precise)
           value = precise
           precise = precise * delta
         end do
         precise = value
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
c     ##  subroutine precond  --  precondition linear CG method  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "precond" solves a simplified version of the Newton equations
c     Ms = r, and uses the result to precondition linear conjugate
c     gradient iterations on the full Newton equations in "solve"
c
c     reference for incomplete Cholesky factorization :
c
c     T. A. Manteuffel, "An Incomplete Factorization Technique
c     for Positive Definite Linear Systems", Mathematics of
c     Computation, 34, 473-497 (1980); the present method is
c     based upon the SICCG(0) method described in this paper
c
c     types of preconditioning methods :
c
c     none     use no preconditioning at all
c     diag     exact Hessian diagonal preconditioning
c     block    3x3 block diagonal preconditioning
c     ssor     symmetric successive over-relaxation
c     iccg     shifted incomplete Cholesky factorization
c
c
      subroutine precond (method,iter,nvar,s,r,h,h_init,
     &                    h_stop,h_index,h_diag, f,c_index,c_value,
     &                    maxhess)
      implicit none
      include 'sizes.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,k,ii,kk,iii,kkk,iter,maxhess
      integer nvar,nblock,ix,iy,iz,icount
      integer h_init(maxvar),h_stop(maxvar)
      integer c_init(maxvar),c_stop(maxvar)
      integer h_index(maxhess)
      integer c_index(maxhess),c_value(maxhess)
      real*8 f(maxhess),f_diag(maxvar)
      real*8 h(maxhess),h_diag(maxvar)
      real*8 r(maxvar),s(maxvar)
      real*8 a(6),b(3)
      real*8 omega,factor,diag(maxvar)
      real*8 maxalpha,alpha,f_i,f_k
      character*6 method
      logical stable
      save f_diag,stable
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     use no preconditioning, using M = identity matrix
c
      if (method .eq. 'none') then
         do i = 1, nvar
            s(i) = r(i)
         end do
      end if
c
c     diagonal preconditioning, using M = abs(Hessian diagonal)
c
      if (method .eq. 'diag') then
         do i = 1, nvar
            s(i) = r(i) / abs(h_diag(i))
         end do
      end if
c
c     block diagonal preconditioning with exact atom blocks
c     (using M = 3x3 blocks from diagonal of full Hessian)
c
      if (method .eq. 'block') then
         nblock = 3
         do i = 1, nvar/3
            iz = 3 * i
            iy = iz - 1
            ix = iz - 2
            a(1) = h_diag(ix)
            if (h_index(h_init(ix)) .eq. iy) then
               a(2) = h(h_init(ix))
            else
               a(2) = 0.0d0
            end if
            if (h_index(h_init(ix)+1) .eq. iz) then
               a(3) = h(h_init(ix)+1)
            else
               a(3) = 0.0d0
            end if
            a(4) = h_diag(iy)
            if (h_index(h_init(iy)) .eq. iz) then
               a(5) = h(h_init(iy))
            else
               a(5) = 0.0d0
            end if
            a(6) = h_diag(iz)
            b(1) = r(ix)
            b(2) = r(iy)
            b(3) = r(iz)
            call cholesky (nblock,a,b)
            s(ix) = b(1)
            s(iy) = b(2)
            s(iz) = b(3)
         end do
      end if
c
c     symmetric successive over-relaxation (SSOR) preconditioning
c     (using M = (D/w+U)T * (D/w)-1 * (D/w+U) with 0 < w < 2)
c
      if (method .eq. 'ssor') then
         omega = 1.0d0
         factor = 2.0d0 - omega
         do i = 1, nvar
            s(i) = r(i) * factor
            diag(i) = h_diag(i) / omega
         end do
         do i = 1, nvar
            s(i) = s(i) / diag(i)
            do j = h_init(i), h_stop(i)
               k = h_index(j)
               s(k) = s(k) - h(j)*s(i)
            end do
         end do
         do i = nvar, 1, -1
            s(i) = s(i) * diag(i)
            do j = h_init(i), h_stop(i)
               k = h_index(j)
               s(i) = s(i) - h(j)*s(k)
            end do
            s(i) = s(i) / diag(i)
         end do
      end if
c
c     factorization phase of incomplete cholesky preconditioning
c
      if (method.eq.'iccg' .and. iter.eq.0) then
         call column (nvar,h_init,h_stop,h_index,
     &                c_init,c_stop,c_index,c_value,maxhess)
         stable = .true.
         icount = 0
         maxalpha = 2.1d0
         alpha = -0.001d0
   10    continue
         if (alpha .le. 0.0d0) then
            alpha = alpha + 0.001d0
         else
            alpha = 2.0d0 * alpha
         end if
         if (alpha .gt. maxalpha) then
            stable = .false.
            if (verbose) then
               if (maswrk) write (iout,20)
   20          format (' PRECOND  --  Incomplete Cholesky is',
     &                 ' Unstable, using Diagonal Method')
            end if
         else
            factor = 1.0d0 + alpha
            do i = 1, nvar
               f_diag(i) = factor * h_diag(i)
               do j = c_init(i), c_stop(i)
                  k = c_index(j)
                  f_i = f(c_value(j))
                  f_diag(i) = f_diag(i) - f_i*f_i*f_diag(k)
                  icount = icount + 1
               end do
               if (f_diag(i) .le. 0.0d0)  goto 10
               if (f_diag(i) .lt. 1.0d-7)  f_diag(i) = 1.0d-7
               f_diag(i) = 1.0d0 / f_diag(i)
               do j = h_init(i), h_stop(i)
                  k = h_index(j)
                  f(j) = h(j)
                  ii = c_init(i)
                  kk = c_init(k)
                  dowhile (ii.le.c_stop(i) .and. kk.le.c_stop(k))
                     iii = c_index(ii)
                     kkk = c_index(kk)
                     if (iii .lt. kkk) then
                        ii = ii + 1
                     else if (kkk .lt. iii) then
                        kk = kk + 1
                     else
                        f_i = f(c_value(ii))
                        f_k = f(c_value(kk))
                        f(j) = f(j) - f_i*f_k*f_diag(iii)
                        ii = ii + 1
                        kk = kk + 1
                        icount = icount + 1
                     end if
                  end do
               end do
            end do
            if (verbose) then
               if (maswrk) write (iout,30)  icount,alpha
   30          format (' PRECOND  --  Incomplete Cholesky',i12,
     &                 ' Operations',f8.3,' Alpha Value')
            end if
         end if
      end if
c
c     solution phase of incomplete cholesky preconditioning
c
      if (method .eq. 'iccg') then
         if (stable) then
            do i = 1, nvar
               s(i) = r(i)
            end do
            do i = 1, nvar
               s(i) = s(i) * f_diag(i)
               do j = h_init(i), h_stop(i)
                  k = h_index(j)
                  s(k) = s(k) - f(j)*s(i)
               end do
            end do
            do i = nvar, 1, -1
               s(i) = s(i) / f_diag(i)
               do j = h_init(i), h_stop(i)
                  k = h_index(j)
                  s(i) = s(i) - f(j)*s(k)
               end do
               s(i) = s(i) * f_diag(i)
            end do
         else
            do i = 1, nvar
               s(i) = r(i) / abs(h_diag(i))
            end do
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
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine pressure  --  maintain constant pressure  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "pressure" uses the internal virial to find the pressure
c     in a periodic box and maintains a constant desired pressure
c     by scaling the coordinates via coupling to an external
c     constant pressure bath
c
c
      subroutine pressure (dt,pres,vol,e_kin)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'boxes.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i
      real*8 dt,e_kin,pres,vol
      real*8 virial,third,scale
c
c
c     find the volume of the periodic box
c
      if (orthogonal) then
         vol = xbox * ybox * zbox
      else if (monoclinic) then
         vol = beta_sin * xbox * ybox * zbox
      else if (triclinic) then
         vol = (gamma_sin*gamma_term) * xbox * ybox * zbox
      else if (octahedron) then
         vol = 0.5d0 * xbox * ybox * zbox
      end if
c
c     compute the value of the internal virial and pressure
c
      virial = virx + viry + virz
      pres = (2.0d0*e_kin - virial) / (3.0d0*vol)
      pres = pres * prescon
c
c     find the scale factor to maintain constant pressure
c
      third = 1.0d0 / 3.0d0
      scale = (1.0d0 + (dt*compress/taupres)*(pres-atmsph))**third
c
c     alter the size of the periodic box
c
      xbox = scale * xbox
      ybox = scale * ybox
      zbox = scale * zbox
      xbox2 = 0.5d0 * xbox
      ybox2 = 0.5d0 * ybox
      zbox2 = 0.5d0 * zbox
      if (octahedron)  box34 = 0.75d0 * xbox
c
c     couple to external pressure bath via atom scaling
c
      do i = 1, n
         if (use(i)) then
            x(i) = scale * x(i)
            y(i) = scale * y(i)
            z(i) = scale * z(i)
         end if
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prmkey  --  interpret force field keywords  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "field" parses a text string to extract keywords related to
c     force field potential energy functional forms and constants
c
c
      subroutine prmkey (text)
      implicit none
      include 'sizes.i'
      include 'angpot.i'
      include 'bndpot.i'
      include 'chgpot.i'
      include 'fields.i'
      include 'polpot.i'
      include 'potent.i'
      include 'rxnpot.i'
      include 'torpot.i'
      include 'urypot.i'
      include 'vdwpot.i'
      integer next
      character*4 value
      character*20 keyword
      character*80 text,record,string
c
c
c     parse the line to extract any possible keyword
c
      record = text
      next = 1
      call upcase (record)
      call gettext (record,keyword,next)
c
c     select the individual force field potential terms
c
      if (keyword(1:9) .eq. 'BONDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_bond = .true.
         if (value .eq. 'NONE')  use_bond = .false.
      else if (keyword(1:10) .eq. 'ANGLETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_angle = .true.
         if (value .eq. 'NONE')  use_angle = .false.
      else if (keyword(1:11) .eq. 'STRBNDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_strbnd = .true.
         if (value .eq. 'NONE')  use_strbnd = .false.
      else if (keyword(1:9) .eq. 'UREYTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_urey = .true.
         if (value .eq. 'NONE')  use_urey = .false.
      else if (keyword(1:11) .eq. 'ANGANGTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_angang = .true.
         if (value .eq. 'NONE')  use_angang = .false.
      else if (keyword(1:11) .eq. 'OPBENDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_opbend = .true.
         if (value .eq. 'NONE')  use_opbend = .false.
      else if (keyword(1:13) .eq. 'IMPROPERTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_improp = .true.
         if (value .eq. 'NONE')  use_improp = .false.
      else if (keyword(1:12) .eq. 'IMPTORSTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_imptor = .true.
         if (value .eq. 'NONE')  use_imptor = .false.
      else if (keyword(1:12) .eq. 'TORSIONTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_tors = .true.
         if (value .eq. 'NONE')  use_tors = .false.
      else if (keyword(1:11) .eq. 'STRTORTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_strtor = .true.
         if (value .eq. 'NONE')  use_strtor = .false.
      else if (keyword(1:11) .eq. 'TORTORTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_tortor = .true.
         if (value .eq. 'NONE')  use_tortor = .false.
      else if (keyword(1:8) .eq. 'VDWTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_vdw = .true.
         if (value .eq. 'NONE')  use_vdw = .false.
      else if (keyword(1:11) .eq. 'CHARGETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_charge = .true.
         if (value .eq. 'NONE')  use_charge = .false.
      else if (keyword(1:11) .eq. 'CHGDPLTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_chgdpl = .true.
         if (value .eq. 'NONE')  use_chgdpl = .false.
      else if (keyword(1:11) .eq. 'DIPOLETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_dipole = .true.
         if (value .eq. 'NONE')  use_dipole = .false.
      else if (keyword(1:10) .eq. 'MPOLETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_mpole = .true.
         if (value .eq. 'NONE')  use_mpole = .false.
      else if (keyword(1:13) .eq. 'POLARIZETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_polar = .true.
         if (value .eq. 'NONE')  use_polar = .false.
      else if (keyword(1:13) .eq. 'RXNFIELDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_rxnfld = .true.
         if (value .eq. 'NONE')  use_rxnfld = .false.
      else if (keyword(1:12) .eq. 'SOLVATETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_solv = .true.
         if (value .eq. 'NONE')  use_solv = .false.
      else if (keyword(1:13) .eq. 'RESTRAINTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_geom = .true.
         if (value .eq. 'NONE')  use_geom = .false.
      else if (keyword(1:10) .eq. 'EXTRATERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_extra = .true.
         if (value .eq. 'NONE')  use_extra = .false.
      end if
c
c     set the name of the force field parameter set
c
      if (keyword(1:11) .eq. 'FORCEFIELD ') then
         call getword (record,forcefield,next)
c
c     definition of the bond stretching potential term
c
      else if (keyword(1:9) .eq. 'BONDTYPE ') then
         call getword (record,bndtyp,next)
      else if (keyword(1:9) .eq. 'BONDUNIT ') then
         string = record(next:80)
         read (string,*,err=10)  bndunit
      else if (keyword(1:11) .eq. 'BOND-CUBIC ') then
         string = record(next:80)
         read (string,*,err=10)  cbnd
      else if (keyword(1:13) .eq. 'BOND-QUARTIC ') then
         string = record(next:80)
         read (string,*,err=10)  qbnd
c
c     definition of the bond angle bending potential term
c
      else if (keyword(1:10) .eq. 'ANGLETYPE ') then
         call getword (record,angtyp,next)
      else if (keyword(1:10) .eq. 'ANGLEUNIT ') then
         string = record(next:80)
         read (string,*,err=10)  angunit
      else if (keyword(1:12) .eq. 'ANGLE-CUBIC ') then
         string = record(next:80)
         read (string,*,err=10)  cang
      else if (keyword(1:14) .eq. 'ANGLE-QUARTIC ') then
         string = record(next:80)
         read (string,*,err=10)  qang
      else if (keyword(1:13) .eq. 'ANGLE-PENTIC ') then
         string = record(next:80)
         read (string,*,err=10)  pang
      else if (keyword(1:13) .eq. 'ANGLE-SEXTIC ') then
         string = record(next:80)
         read (string,*,err=10)  sang
c
c     definitions for other local geometry potential terms
c
      else if (keyword(1:11) .eq. 'STRBNDUNIT ') then
         string = record(next:80)
         read (string,*,err=10)  stbnunit
      else if (keyword(1:9) .eq. 'UREYUNIT ') then
         string = record(next:80)
         read (string,*,err=10)  ureyunit
      else if (keyword(1:11) .eq. 'ANGANGUNIT ') then
         string = record(next:80)
         read (string,*,err=10)  aaunit
      else if (keyword(1:11) .eq. 'OPBENDUNIT ') then
         string = record(next:80)
         read (string,*,err=10)  opbunit
      else if (keyword(1:12) .eq. 'TORSIONUNIT ') then
         string = record(next:80)
         read (string,*,err=10)  torsunit
      else if (keyword(1:11) .eq. 'STRTORUNIT ') then
         string = record(next:80)
         read (string,*,err=10)  storunit
c
c     definition of the van der Waals potential term
c
      else if (keyword(1:8) .eq. 'VDWTYPE ') then
         call getword (record,vdwtyp,next)
      else if (keyword(1:10) .eq. 'A-EXPTERM ') then
         string = record(next:80)
         read (string,*,err=10)  aterm
      else if (keyword(1:10) .eq. 'B-EXPTERM ') then
         string = record(next:80)
         read (string,*,err=10)  bterm
      else if (keyword(1:10) .eq. 'C-EXPTERM ') then
         string = record(next:80)
         read (string,*,err=10)  cterm
      else if (keyword(1:11) .eq. 'VDW-12-USE ') then
         call getword (record,value,next)
         vdw12use = -1
         if (value .eq. 'NONE')  vdw12use = 1
      else if (keyword(1:11) .eq. 'VDW-13-USE ') then
         call getword (record,value,next)
         vdw13use = -1
         if (value .eq. 'NONE')  vdw13use = 1
      else if (keyword(1:11) .eq. 'VDW-14-USE ') then
         call getword (record,value,next)
         vdw14use = -1
         if (value .eq. 'NONE')  vdw14use = 1
      else if (keyword(1:10) .eq. 'VDW-SCALE ') then
         string = record(next:80)
         read (string,*,err=10)  vdwscale
      else if (keyword(1:11) .eq. 'RADIUSTYPE ') then
         call getword (record,radtyp,next)
      else if (keyword(1:11) .eq. 'RADIUSSIZE ') then
         call getword (record,radsiz,next)
      else if (keyword(1:11) .eq. 'RADIUSRULE ') then
         call getword (record,radrule,next)
      else if (keyword(1:12) .eq. 'EPSILONRULE ') then
         call getword (record,epsrule,next)
      else if (keyword(1:14) .eq. 'GAUSSTYPE ') then
         call getword (record,gausstyp,next)
c
c     definition of the electrostatic potential term
c
      else if (keyword(1:11) .eq. 'DIELECTRIC ') then
         string = record(next:80)
         read (string,*,err=10)  dielec
      else if (keyword(1:11) .eq. 'CHG-12-USE ') then
         call getword (record,value,next)
         chg12use = -1
         if (value .eq. 'NONE')  chg12use = 1
      else if (keyword(1:11) .eq. 'CHG-13-USE ') then
         call getword (record,value,next)
         chg13use = -1
         if (value .eq. 'NONE')  chg13use = 1
      else if (keyword(1:11) .eq. 'CHG-14-USE ') then
         call getword (record,value,next)
         chg14use = -1
         if (value .eq. 'NONE')  chg14use = 1
      else if (keyword(1:10) .eq. 'CHG-SCALE ') then
         string = record(next:80)
         read (string,*,err=10)  chgscale
      else if (keyword(1:15) .eq. 'NEUTRAL-GROUPS ') then
         neutcut = .true.
      else if (keyword(1:16) .eq. 'GROUP-NEIGHBORS ') then
         neutnbr = .true.
c
c     definition of the polarization potential term
c
      else if (keyword(1:13) .eq. 'POLARIZATION ') then
         call getword (record,poltyp,next)
      else if (keyword(1:10) .eq. 'POLAR-EPS ') then
         string = record(next:80)
         read (string,*,err=10)  poleps
      else if (keyword(1:11) .eq. 'POLAR-DAMP ') then
         string = record(next:80)
         read (string,*,err=10)  pradius,pgamma
c
c     definition of the reaction field potential term
c
      else if (keyword(1:14) .eq. 'REACTIONFIELD ') then
         string = record(next:80)
         read (string,*,err=10)  rfsize,rfbulkd,rfterms
      end if
c
c     jump directly to the end if any error was detected
c
   10 continue
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine potoff  --  turn off all potential functions  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "potoff" clears the forcefield definition by turning off
c     the use of each of the potential energy functions
c
c
      subroutine potoff
      implicit none
      include 'potent.i'
c
c
c     turn off the use of each of the potential energy functions
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_improp = .false.
      use_imptor = .false.
      use_tors = .false.
      use_strtor = .false.
      use_vdw = .false.
      use_charge = .false.
      use_chgdpl = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_polar = .false.
      use_rxnfld = .false.
      use_solv = .false.
      use_geom = .false.
      use_extra = .false.
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
c     ##  subroutine promo  --  copywrite notice and program info  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "promo" writes a short message containing info about
c     the TINKER program package and the copyright notice
c
c
      subroutine promo
      implicit none
      include 'iounit.i'
      LOGICAL GOPARR,DSKWRK,MASWRK
      integer me,master,nproc,ibtyp,iptim
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     write the short informational message
c
      if (maswrk) write (iout,10)
   10 format (/,' ',78('#'),
     &        /,' ',78('#'),
     &        /,' ##',74x,'##',
     &        /,' ##',13x,'TINKER  ---  Software Tools for',
     &            ' Molecular Design',13x,'##',
     &        /,' ##',74x,'##',
     &        /,' ##',24x,'Version 3.6  February 1998',24x,'##',
     &        /,' ##',74x,'##',
     &        /,' ##',15x,'Copyright (c)  Jay William Ponder',
     &            '  1990-1998',15x,'##',
     &        /,' ##',28x,'All Rights Reserved',27x,'##',
     &        /,' ##',74x,'##',
     &        /' --- last minor change for GAMESS/SIMOMM=April 2010'
     &        /,' ',78('#'),
     &        /,' ',78('#'),/)
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtdyn  --  output of MD restart information  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtdyn" writes out the information needed to restart a
c     molecular dynamics trajectory to an external disk file
c
c
      subroutine prtdyn
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'files.i'
      include 'moldyn.i'
      include 'titles.i'
      integer i,idyn
      integer freeunit
      character*60 dynfile
      logical exist
c
c
c     overwrite an existing restart file or open a new one
c
      idyn = freeunit ()
      dynfile = filename(1:leng)//'.dyn'
      inquire (file=dynfile,exist=exist)
      if (exist) then
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
      else
         open (unit=idyn,file=dynfile,status='new')
      end if
c
c     save the number of atoms and the title string
c
      write (idyn,10)
   10 format (' Number of Atoms and Title :')
      if (ltitle .eq. 0) then
         write (idyn,20)  n
   20    format (i6)
      else
         write (idyn,30)  n,ttitle(1:ltitle)
   30    format (i6,2x,a)
      end if
c
c     save the periodic box edge lengths and angles
c
      write (idyn,40)
   40 format (' Periodic Box Dimensions :')
      write (idyn,50)  xbox,ybox,zbox
   50 format (3d26.16)
      write (idyn,60)  alpha,beta,gamma
   60 format (3d26.16)
c
c     save the atomic positions, velocities and accelerations
c
      write (idyn,70)
   70 format (' Current Atomic Positions :')
      do i = 1, n
         write (idyn,80)  x(i),y(i),z(i)
   80    format (3d26.16)
      end do
      write (idyn,90)
   90 format (' Current Atomic Velocities :')
      do i = 1, n
         write (idyn,100)  v(1,i),v(2,i),v(3,i)
  100    format (3d26.16)
      end do
      write (idyn,110)
  110 format (' Current Atomic Accelerations :')
      do i = 1, n
         write (idyn,120)  a(1,i),a(2,i),a(3,i)
  120    format (3d26.16)
      end do
      write (idyn,130)
  130 format (' Previous Atomic Accelerations :')
      do i = 1, n
         write (idyn,140)  a_old(1,i),a_old(2,i),a_old(3,i)
  140    format (3d26.16)
      end do
      close (unit=idyn)
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prterr  --  output coordinates upon error  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prterr" writes out a set of coordinates to a disk
c     file prior to aborting on a serious error
c
c
      subroutine prterr
      implicit none
      include 'files.i'
      include 'output.i'
      integer ierr,freeunit
      character*60 errorfile
c
c
c     write the current coordinates to a file after an error
c
      ierr = freeunit ()
      errorfile = filename(1:leng)//'.err'
      call version (errorfile,'new')
      open (unit=ierr,file=errorfile,status='new')
      if (coordtype .eq. 'cartesian') then
         call prtxyz (ierr)
      else if (coordtype .eq. 'internal') then
         call prtint (ierr)
      else if (coordtype .eq. 'rigidbody') then
         call prtxyz (ierr)
      end if
      close (unit=ierr)
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
c     ##  subroutine prtint  --  output of internal coordinates  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prtint" writes out a set of Z-matrix internal
c     coordinates to an external disk file
c
c
      subroutine prtint (izmt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'titles.i'
      include 'zclose.i'
      include 'zcoord.i'
      integer i,k,izmt
      character*60 coordfile
      logical opened
c
c
c     open output unit if not already done
c
      inquire (unit=izmt,opened=opened)
      if (.not. opened) then
         coordfile = filename(1:leng)//'.int'
         call version (coordfile,'new')
         open (unit=izmt,file=coordfile,status='new')
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         write (izmt,10)  n
   10    format (i6)
      else
         write (izmt,20)  n,ttitle(1:ltitle)
   20    format (i6,2x,a)
      end if
c
c     first three atoms are special cases
c
      if (n .ge. 1) then
         write (izmt,30)  1,name(1),type(1)
   30    format (i6,2x,a3,i5)
      end if
      if (n .ge. 2) then
         write (izmt,40)  2,name(2),type(2),iz(1,2),zbond(2)
   40    format (i6,2x,a3,i5,i6,f10.5)
      end if
      if (n .ge. 3) then
         write (izmt,50)  3,name(3),type(3),iz(1,3),zbond(3),
     &                    iz(2,3),zang(3)
   50    format (i6,2x,a3,i5,i6,f10.5,i6,f10.4)
      end if
c
c     now, the fourth through final atoms
c
      do i = 4, n
         if (iz(4,i) .eq. 0) then
            dowhile (ztors(i) .lt. -180.0d0)
               ztors(i) = ztors(i) + 360.0d0
            end do
            dowhile (ztors(i) .gt. 180.0d0)
               ztors(i) = ztors(i) - 360.0d0
            end do
         end if
         write (izmt,60)  i,name(i),type(i),iz(1,i),zbond(i),
     &                    iz(2,i),zang(i),iz(3,i),ztors(i),iz(4,i)
   60    format (i6,2x,a3,i5,i6,f10.5,i6,f10.4,i6,f10.4,i6)
      end do
c
c     finally, add or delete bonds as required
c
      if (nadd.ne.0 .or. ndel.ne.0)  write (izmt,70)
   70 format ()
      do i = 1, nadd
         write (izmt,80)  (iadd(k,i),k=1,2)
   80    format (2i6)
      end do
      if (ndel .ne. 0)  write (izmt,90)
   90 format ()
      do i = 1, ndel
         write (izmt,100)  (idel(k,i),k=1,2)
  100    format (2i6)
      end do
      if (.not. opened)  close (unit=izmt)
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
c     ##  program prtmol2  --  output of Sybyl MOL2 coordinate file  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "prtmol2" writes out a set of coordinates in Sybyl MOL2
c     format to an external disk file
c
c
      subroutine prtmol2 (isyb)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'files.i'
      include 'iounit.i'
      include 'titles.i'
      integer i,j,k,isyb
      integer thousand,hundred,tens,ones
      character*1 digit(0:9)
      character*4 atmtyp,number
      character*7 atmnam
      character*60 sybylfile
      logical opened
      data digit / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     open output unit if not already done
c
      inquire (unit=isyb,opened=opened)
      if (.not. opened) then
         sybylfile = filename(1:leng)//'.mol2'
         call version (sybylfile,'new')
         open (unit=isyb,file=sybylfile,status='new')
      end if
c
c     write the molecule record type indicator
c
      write (isyb,10)
   10 format ('@<TRIPOS>MOLECULE')
      if (ltitle .eq. 0) then
         write (isyb,20)
   20    format ('****')
      else
         write (isyb,30)  ttitle(1:ltitle)
   30    format (a)
      end if
      write (isyb,40)  n,nbond,1
   40 format (3i8)
      write (isyb,50)
   50 format ('SMALL')
      write (isyb,60)
   60 format ('NO_CHARGES')
c
c     write the atom record type indicator
c
      write (isyb,70)
   70 format (/,'@<TRIPOS>ATOM')
      do i = 1, n
c
c     set the Sybyl atom_name for the atom
c
         thousand = i / 1000
         hundred = (i - 1000*thousand) / 100
         tens = (i - 1000*thousand - 100*hundred) / 10
         ones = i - 1000*thousand - 100*hundred - 10*tens
         number(1:1) = digit(thousand)
         number(2:2) = digit(hundred)
         number(3:3) = digit(tens)
         number(4:4) = digit(ones)
         if (number(1:1) .eq. '0')  number(1:1)=' '
         if (number(2:2).eq.'0' .and. number(1:1).eq.' ') then
            number(2:2) = ' '
         end if
         if (number(3:3).eq.'0' .and. number(2:2).eq.' ') then
            number(3:3) = ' '
         end if
         atmnam = name(i)//number
         do j = 1, 6
            dowhile (atmnam(j:j) .eq. ' ')
               do k = j, 6
                  atmnam(k:k) = atmnam(k+1:k+1)
               end do
               atmnam(7:7) = '*'
            end do
         end do
         do j = 1, 7
            if (atmnam(j:j) .eq. '*')  atmnam(j:j) = ' '
         end do
c
c     set the Sybyl atom_type for the atom
c
         atmtyp = name(i)//' '
         if (atmtyp .eq. 'C  ') then
            if (n12(i) .eq. 4)  atmtyp = 'C.3 '
            if (n12(i) .eq. 3)  atmtyp = 'C.2 '
            if (n12(i) .eq. 2)  atmtyp = 'C.1 '
         else if (atmtyp .eq. 'N  ') then
            if (n12(i) .ge. 3)  atmtyp = 'N.3 '
            if (n12(i) .eq. 2)  atmtyp = 'N.2 '
            if (n12(i) .eq. 1)  atmtyp = 'N.1 '
         else if (atmtyp .eq. 'N+ ') then
            atmtyp = 'N.4 '
         else if (atmtyp .eq. 'O  ') then
            if (n12(i) .ge. 2)  atmtyp = 'O.3 '
            if (n12(i) .le. 1)  atmtyp = 'O.2 '
         else if (atmtyp .eq. 'O- ') then
            atmtyp = 'O.2 '
         else if (atmtyp .eq. 'S  ') then
            if (n12(i) .ge. 2)  atmtyp = 'S.3 '
            if (n12(i) .le. 1)  atmtyp = 'S.2 '
         else if (atmtyp .eq. 'P  ') then
            atmtyp = 'P.3 '
         else if (atmtyp .eq. 'Lp ') then
            atmtyp = 'LP  '
         end if
         write (isyb,80)  i,atmnam,x(i),y(i),z(i),atmtyp
   80    format (i8,3x,a7,2x,3f12.6,3x,a4)
      end do
c
c     write the bond record type indicator
c
      write (isyb,90)
   90 format (/,'@<TRIPOS>BOND')
      do i = 1, nbond
         write (isyb,100)  i,(ibnd(j,i),j=1,2),1
  100    format (4i8)
      end do
c
c     write the substructure record type indicator
c
      write (isyb,110)
  110 format (/,'@<TRIPOS>SUBSTRUCTURE')
      write (isyb,120)  1,'****',1
  120 format (i8,12x,a4,i8)
      if (.not. opened)  close (unit=isyb)
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
c     ##  subroutine prtpdb  --  output of Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtpdb" writes out a set of Protein Data Bank
c     coordinates to an external disk file
c
c
      subroutine prtpdb (ipdb)
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'pdb.i'
      include 'sequen.i'
      include 'titles.i'
      integer i,k,ipdb
      integer start,stop,resnumb
      integer resid(maxres)
      character*1 chnname
      character*1 chain(maxres)
      character*60 pdbfile
      logical opened
c
c
c     open output unit if not already done
c
      inquire (unit=ipdb,opened=opened)
      if (.not. opened) then
         pdbfile = filename(1:leng)//'.pdb'
         call version (pdbfile,'new')
         open (unit=ipdb,file=pdbfile,status='new')
      end if
c
c     write out the header lines and the title
c
      if (ltitle .eq. 0) then
         write (ipdb,10)
   10    format ('HEADER',/,'COMPND',/,'SOURCE')
      else
         write (ipdb,20)  ttitle(1:ltitle)
   20    format ('HEADER',4x,a,/,'COMPND',/,'SOURCE')
      end if
c
c     find the chain name and chain position for each residue
c
      do i = 1, nchain
         start = ichain(1,i)
         stop = ichain(2,i)
         do k = start, stop
            resid(k) = k - start + 1
            chain(k) = chnnam(i)
         end do
      end do
c
c     next, write the coordinates for each PDB atom
c
      do i = 1, npdb
         if (pdbtyp(i) .eq. 'ATOM  ') then
            resnumb = resid(resnum(i))
            chnname = chain(resnum(i))
         else
            resnumb = resnum(i) + nseq
            chnname = ' '
         end if
         if (resnam(i) .eq. 'CYX')  resnam(i) = 'CYS'
         if (resnam(i) .eq. 'HIP')  resnam(i) = 'HIS'
         write (ipdb,30)  pdbtyp(i),i,atmnam(i),resnam(i),chnname,
     &                    resnumb,xpdb(i),ypdb(i),zpdb(i)
   30    format (a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
      end do
c
c     finally, write any connectivity records for PDB atoms
c
      do i = 1, npdb
         if (npdb12(i) .ne. 0) then
            write (ipdb,40)  i,(ipdb12(k,i),k=1,npdb12(i))
   40       format ('CONECT',5i5)
         end if
      end do
      write (ipdb,50)
   50 format ('END')
c     close (unit=ipdb)
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
c     ##  subroutine prtprm  --  output of force field parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtprm" writes out a formatted listing of the default
c     set of potential energy parameters for a force field
c
c
      subroutine prtprm (itxt)
      implicit none
      include 'sizes.i'
      include 'angpot.i'
      include 'bndpot.i'
      include 'fields.i'
      include 'kanang.i'
      include 'kangs.i'
      include 'katoms.i'
      include 'kbonds.i'
      include 'kchrge.i'
      include 'kdipol.i'
      include 'khbond.i'
      include 'kiprop.i'
      include 'kitors.i'
      include 'kmulti.i'
      include 'kopbnd.i'
      include 'korbs.i'
      include 'kpolr.i'
      include 'kstbnd.i'
      include 'ksttor.i'
      include 'ktorsn.i'
      include 'kurybr.i'
      include 'kvdws.i'
      include 'kvdwpr.i'
      integer i,j,k,itxt,number
      integer k1,k2,k3,k4
      integer fold(6)
      real*8 ampli(6),phase(6)
      character*1 formfeed
      logical exist
c
c
c     set the string value of the formfeed character (Ctrl-L)
c
      formfeed = char(12)
c
c     force field atom type definitions
c
      exist = .false.
      do i = 1, maxtyp
         if (symbol(i) .ne. '   ')  exist = .true.
      end do
      if (exist) then
         write (itxt,10)  forcefield
   10    format (//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,20)
   20    format (//,17x,'Force Field Atom Definitions',
     &           //,52x,'Atomic',4x,'Atomic',
     &           /,7x,'Type',3x,'Class',3x,'Symbol',3x,'Description',
     &              10x,'Number',4x,'Weight',3x,'Valence',/)
         do i = 1, maxtyp
            if (symbol(i) .ne. '   ') then
               write (itxt,30)  i,atmcls(i),symbol(i),describe(i),
     &                          atmnum(i),weight(i),ligand(i)
   30          format (5x,i5,3x,i5,5x,a3,5x,a20,i5,f12.3,i7)
            end if
         end do
      end if
c
c     van der Waals parameters for atom types
c
      exist = .false.
      do i = 1, maxclass
         if (rad(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,40)  formfeed,forcefield
   40    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,50)
   50    format (//,17x,'Van der Waals Parameters',
     &           ///,22x,'Class',7x,'Radius',6x,'Epsilon',
     &                 4x,'Reduction',/)
         k = 0
         do i = 1, maxclass
            if (rad(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,60)  k,i,rad(i),eps(i),reduct(i)
   60          format (13x,i5,5x,i3,2x,3f12.3)
            end if
         end do
      end if
c
c     van der Waals parameters for specific atom pairs
c
      if (kvpr(1) .ne. '      ') then
         write (itxt,70)
   70    format (//,17x,'Van der Waals Parameters for Atom Pairs',
     &           ///,23x,'Classes',7x,'Radii Sum',4x,'Epsilon',/)
         do i = 1, maxnvp
            if (kvpr(i) .eq. '      ')  goto 90
            k1 = number (kvpr(i)(1:3))
            k2 = number (kvpr(i)(4:6))
            write (itxt,80)  i,k1,k2,radpr(i),epspr(i)
   80       format (13x,i5,5x,i3,'-',i3,2x,2f12.3)
         end do
   90    continue
      end if
c
c     hydrogen bonding parameters for specific atom pairs
c
      if (khb(1) .ne. '      ') then
         write (itxt,100)
  100    format (//,17x,'Hydrogen Bonding Parameters for Atom Pairs',
     &           ///,23x,'Classes',7x,'Radii Sum',4x,'Epsilon',/)
         do i = 1, maxnhb
            if (khb(i) .eq. '      ')  goto 120
            k1 = number (khb(i)(1:3))
            k2 = number (khb(i)(4:6))
            write (itxt,110)  i,k1,k2,radhb(i),epshb(i)
  110       format (13x,i5,5x,i3,'-',i3,2x,2f12.3)
         end do
  120    continue
      end if
c
c     bond stretching parameters
c
      if (kb(1) .ne. '      ') then
         write (itxt,130)  formfeed,forcefield
  130    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,140)
  140    format (//,17x,'Bond Stretching Parameters',
     &           ///,23x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb
            if (kb(i) .eq. '      ')  goto 160
            k1 = number (kb(i)(1:3))
            k2 = number (kb(i)(4:6))
            write (itxt,150)  i,k1,k2,fcon(i),blen(i)
  150       format (13x,i5,5x,i3,'-',i3,6x,f12.3,f12.4)
         end do
  160    continue
      end if
c
c     bond stretching parameters for 5-membered rings
c
      if (kb5(1) .ne. '      ') then
         write (itxt,170)
  170    format (//,17x,'5-Membered Ring Stretch Parameters',
     &           ///,23x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb5
            if (kb5(i) .eq. '      ')  goto 190
            k1 = number (kb5(i)(1:3))
            k2 = number (kb5(i)(4:6))
            write (itxt,180)  i,k1,k2,fcon5(i),blen5(i)
  180       format (13x,i5,5x,i3,'-',i3,6x,f12.3,f12.4)
         end do
  190    continue
      end if
c
c     bond stretching parameters for 4-membered rings
c
      if (kb4(1) .ne. '      ') then
         write (itxt,200)
  200    format (//,17x,'4-Membered Ring Stretch Parameters',
     &           ///,23x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb4
            if (kb4(i) .eq. '      ')  goto 220
            k1 = number (kb4(i)(1:3))
            k2 = number (kb4(i)(4:6))
            write (itxt,240)  i,k1,k2,fcon4(i),blen4(i)
  210       format (13x,i5,5x,i3,'-',i3,6x,f12.3,f12.4)
         end do
  220    continue
      end if
c
c     bond stretching parameters for 3-membered rings
c
      if (kb3(1) .ne. '      ') then
         write (itxt,230)
  230    format (//,17x,'3-Membered Ring Stretch Parameters',
     &           ///,23x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb3
            if (kb3(i) .eq. '      ')  goto 250
            k1 = number (kb3(i)(1:3))
            k2 = number (kb3(i)(4:6))
            write (itxt,240)  i,k1,k2,fcon3(i),blen3(i)
  240       format (13x,i5,5x,i3,'-',i3,6x,f12.3,f12.4)
         end do
  250    continue
      end if
c
c     cubic and quartic bond stretching parameters
c
      if (cbnd.ne.0.0d0 .or. qbnd.ne.0.0d0) then
         write (itxt,260)  cbnd,qbnd
  260    format (//,17x,'Higher Order Stretching Constants',
     &           ///,23x,'Cubic',f17.3,/,23x,'Quartic',f15.3)
      end if
c
c     bond angle bending parameters
c
      if (ka(1) .ne. '         ') then
         write (itxt,270)  formfeed,forcefield
  270    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,280)
  280    format (//,17x,'Angle Bending Parameters',
     &           ///,18x,'Classes',10x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,43x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna
            if (ka(i) .eq. '         ')  goto 310
            k1 = number (ka(i)(1:3))
            k2 = number (ka(i)(4:6))
            k3 = number (ka(i)(7:9))
            if (ang(2,i).eq.0.0d0 .and. ang(3,i).eq.0.0d0) then
               write (itxt,290)  i,k1,k2,k3,con(i),ang(1,i)
  290          format (5x,i5,5x,i3,'-',i3,'-',i3,2f12.3)
            else
               write (itxt,300)  i,k1,k2,k3,con(i),(ang(j,i),j=1,3)
  300          format (5x,i5,5x,i3,'-',i3,'-',i3,4f12.3)
            end if
         end do
  310    continue
      end if
c
c     bond angle bending parameters for 5-membered rings
c
      if (ka5(1) .ne. '         ') then
         write (itxt,320)
  320    format (//,17x,'5-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',10x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,43x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna5
            if (ka5(i) .eq. '         ')  goto 350
            k1 = number (ka5(i)(1:3))
            k2 = number (ka5(i)(4:6))
            k3 = number (ka5(i)(7:9))
            if (ang5(2,i).eq.0.0d0 .and. ang5(3,i).eq.0.0d0) then
               write (itxt,330)  i,k1,k2,k3,con5(i),ang5(1,i)
  330          format (5x,i5,5x,i3,'-',i3,'-',i3,2f12.3)
            else
               write (itxt,340)  i,k1,k2,k3,con5(i),(ang5(j,i),j=1,3)
  340          format (5x,i5,5x,i3,'-',i3,'-',i3,4f12.3)
            end if
         end do
  350    continue
      end if
c
c     bond angle bending parameters for 4-membered rings
c
      if (ka4(1) .ne. '         ') then
         write (itxt,360)
  360    format (//,17x,'4-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',10x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,43x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna4
            if (ka4(i) .eq. '         ')  goto 390
            k1 = number (ka4(i)(1:3))
            k2 = number (ka4(i)(4:6))
            k3 = number (ka4(i)(7:9))
            if (ang4(2,i).eq.0.0d0 .and. ang4(3,i).eq.0.0d0) then
               write (itxt,370)  i,k1,k2,k3,con4(i),ang4(1,i)
  370          format (5x,i5,5x,i3,'-',i3,'-',i3,2f12.3)
            else
               write (itxt,380)  i,k1,k2,k3,con4(i),(ang4(j,i),j=1,3)
  380          format (5x,i5,5x,i3,'-',i3,'-',i3,4f12.3)
            end if
         end do
  390    continue
      end if
c
c     bond angle bending parameters for 3-membered rings
c
      if (ka3(1) .ne. '         ') then
         write (itxt,400)
  400    format (//,17x,'3-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',10x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,43x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do  i = 1, maxna3
            if (ka3(i) .eq. '         ')  goto 430
            k1 = number (ka3(i)(1:3))
            k2 = number (ka3(i)(4:6))
            k3 = number (ka3(i)(7:9))
            if (ang3(3,i).eq.0.0d0 .and. ang3(3,i).eq.0.0d0) then
               write (itxt,410)  i,k1,k2,k3,con3(i),ang3(1,i)
  410          format (5x,i5,5x,i3,'-',i3,'-',i3,2f12.3)
            else
               write (itxt,420)  i,k1,k2,k3,con3(i),(ang3(j,i),j=1,3)
  420          format (5x,i5,5x,i3,'-',i3,'-',i3,4f12.3)
            end if
         end do
  430    continue
      end if
c
c     cubic through sextic bond angle bending parameters
c
      if (cang.ne.0.0d0 .or. qang.ne.0.0d0 .or.
     &    pang.ne.0.0d0 .or. sang.ne.0.0d0) then
         write (itxt,440)  cang,qang,pang,sang
  440    format (//,17x,'Higher Order Bending Constants',
     &           ///,21x,'Cubic',d17.3,/,21x,'Quartic',d15.3,
     &           /,21x,'Pentic',d16.3,/,21x,'Sextic',d16.3)
      end if
c
c     stretch-bend parameters
c
      exist = .false.
      do i = 1, maxclass
         do k = 1, 3
            if (stbn(k,i) .ne. 0.0d0)  exist = .true.
         end do
      end do
      if (exist) then
         write (itxt,450)  formfeed,forcefield
  450    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,460)
  460    format (//,17x,'Stretch-Bend Parameters',
     &           ///,21x,'Class',9x,'KSB 1',7x,'KSB 2',7x,'KSB 3',
     &           /,34x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         k = 0
         do i = 1, maxclass
            if (stbn(1,i).ne.0.0d0 .or. stbn(2,i).ne.0.0d0
     &               .or. stbn(3,i).ne.0.0d0) then
               k = k + 1
               write (itxt,470)  k,i,(stbn(j,i),j=1,3)
  470          format (10x,i5,7x,i3,3x,3f12.3)
            end if
         end do
      end if
c
c     out-of-plane bending parameters
c
      if (kaopb(1) .ne. '      ') then
         write (itxt,480)  formfeed,forcefield
  480    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,490)
  490    format (//,17x,'Out-of-Plane Bend Parameters',
     &           ///,23x,'Classes',8x,'KOPB',/)
         do i = 1, maxnopb
            if (kaopb(i) .eq. '      ')  goto 510
            k1 = number (kaopb(i)(1:3))
            k2 = number (kaopb(i)(4:6))
            write (itxt,500)  i,k1,k2,copb(i)
  500       format (10x,i5,8x,i3,'-',i3,f12.3)
         end do
  510    continue
      end if
c
c     Urey-Bradley parameters
c
      if (ku(1) .ne. '         ') then
         write (itxt,520)  formfeed,forcefield
  520    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,530)
  530    format (//,17x,'Urey-Bradley Parameters',
     &           ///,18x,'Classes',10x,'KB',6x,'Distance',/)
         do i = 1, maxnu
            if (ku(i) .eq. '         ')  goto 550
            k1 = number (ku(i)(1:3))
            k2 = number (ku(i)(4:6))
            k3 = number (ku(i)(7:9))
            write (itxt,540)  i,k1,k2,k3,ucon(i),dst13(i)
  540       format (5x,i5,5x,i3,'-',i3,'-',i3,f12.3,f12.4)
         end do
  550    continue
      end if
c
c     angle-angle parameters
c
      exist = .false.
      do i = 1, maxclass
         do k = 1, 3
            if (anan(k,i) .ne. 0.0d0)  exist = .true.
         end do
      end do
      if (exist) then
         write (itxt,560)  formfeed,forcefield
  560    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,570)
  570    format (//,17x,'Angle-Angle Parameters',
     &           ///,21x,'Class',9x,'KAA 1',7x,'KAA 2',7x,'KAA 3',
     &           /,34x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         k = 0
         do i = 1, maxclass
            if (anan(1,i).ne.0.0d0 .or. anan(2,i).ne.0.0d0
     &               .or. anan(3,i).ne.0.0d0) then
               k = k + 1
               write (itxt,580)  k,i,(anan(j,i),j=1,3)
  580          format (10x,i5,7x,i3,3x,3f12.3)
            end if
         end do
      end if
c
c     improper dihedral parameters
c
      if (kdi(1) .ne. '            ') then
         write (itxt,590)  formfeed,forcefield
  590    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,600)
  600    format (//,17x,'Improper Dihedral Parameters',
     &           ///,20x,'Classes',11x,'KID',7x,'Target',/)
         do i = 1, maxndi
            if (kdi(i) .eq. '            ')  goto 620
            k1 = number (kdi(i)(1:3))
            k2 = number (kdi(i)(4:6))
            k3 = number (kdi(i)(7:9))
            k4 = number (kdi(i)(10:12))
            write (itxt,610)  i,k1,k2,k3,k4,dcon(i),tdi(i)
  610       format (5x,i5,5x,i3,'-',i3,'-',i3,'-',i3,f12.3,f12.4)
         end do
  620    continue
      end if
c
c     improper torsional parameters
c
      if (kti(1) .ne. '            ') then
         write (itxt,630)  formfeed,forcefield
  630    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,640)
  640    format (//,17x,'Improper Torsion Parameters',
     &           ///,18x,'Classes',15x,'KTI Values',/)
         do i = 1, maxnti
            if (kti(i) .eq. '            ')  goto 660
            k1 = number (kti(i)(1:3))
            k2 = number (kti(i)(4:6))
            k3 = number (kti(i)(7:9))
            k4 = number (kti(i)(10:12))
            j = 0
            if (ti1(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = ti1(1,i)
               phase(j) = ti1(2,i)
            end if
            if (ti2(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = ti2(1,i)
               phase(j) = ti2(2,i)
            end if
            if (ti3(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = ti3(1,i)
               phase(j) = ti3(2,i)
            end if
            write (itxt,650)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  650       format (3x,i5,5x,i3,'-',i3,'-',i3,'-',i3,2x,3(f8.3,f6.1,i2))
         end do
  660    continue
      end if
c
c     torsional angle parameters
c
      if (kt(1) .ne. '            ') then
         write (itxt,670)  formfeed,forcefield
  670    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,680)
  680    format (//,17x,'Torsional Parameters',
     &           ///,18x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt
            if (kt(i) .eq. '            ')  goto 700
            k1 = number (kt(i)(1:3))
            k2 = number (kt(i)(4:6))
            k3 = number (kt(i)(7:9))
            k4 = number (kt(i)(10:12))
            j = 0
            if (t1(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t1(1,i)
               phase(j) = t1(2,i)
            end if
            if (t2(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t2(1,i)
               phase(j) = t2(2,i)
            end if
            if (t3(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t3(1,i)
               phase(j) = t3(2,i)
            end if
            if (t4(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t4(1,i)
               phase(j) = t4(2,i)
            end if
            if (t5(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t5(1,i)
               phase(j) = t5(2,i)
            end if
            if (t6(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t6(1,i)
               phase(j) = t6(2,i)
            end if
            write (itxt,690)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  690       format (3x,i5,5x,i3,'-',i3,'-',i3,'-',i3,2x,6(f8.3,f6.1,i2))
         end do
  700    continue
      end if
c
c     torsional angle parameters for 5-membered rings
c
      if (kt5(1) .ne. '            ') then
         write (itxt,710)
  710    format (//,17x,'5-Membered Ring Torsion Parameters',
     &           ///,18x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt5
            if (kt5(i) .eq. '            ')  goto 730
            k1 = number (kt5(i)(1:3))
            k2 = number (kt5(i)(4:6))
            k3 = number (kt5(i)(7:9))
            k4 = number (kt5(i)(10:12))
            j = 0
            if (t15(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t15(1,i)
               phase(j) = t15(2,i)
            end if
            if (t25(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t25(1,i)
               phase(j) = t25(2,i)
            end if
            if (t35(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t35(1,i)
               phase(j) = t35(2,i)
            end if
            if (t45(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t45(1,i)
               phase(j) = t45(2,i)
            end if
            if (t55(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t55(1,i)
               phase(j) = t55(2,i)
            end if
            if (t65(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t65(1,i)
               phase(j) = t65(2,i)
            end if
            write (itxt,720)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  720       format (3x,i5,5x,i3,'-',i3,'-',i3,'-',i3,2x,6(f8.3,f6.1,i2))
         end do
  730    continue
      end if
c
c     torsional angle parameters for 4-membered rings
c
      if (kt4(1) .ne. '            ') then
         write (itxt,740)
  740    format (//,17x,'4-Membered Ring Torsion Parameters',
     &           ///,18x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt4
            if (kt4(i) .eq. '            ')  goto 760
            k1 = number (kt4(i)(1:3))
            k2 = number (kt4(i)(4:6))
            k3 = number (kt4(i)(7:9))
            k4 = number (kt4(i)(10:12))
            j = 0
            if (t14(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t14(1,i)
               phase(j) = t14(2,i)
            end if
            if (t24(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t24(1,i)
               phase(j) = t24(2,i)
            end if
            if (t34(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t34(1,i)
               phase(j) = t34(2,i)
            end if
            if (t44(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t44(1,i)
               phase(j) = t44(2,i)
            end if
            if (t54(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t54(1,i)
               phase(j) = t54(2,i)
            end if
            if (t64(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t64(1,i)
               phase(j) = t64(2,i)
            end if
            write (itxt,750)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  750       format (3x,i5,5x,i3,'-',i3,'-',i3,'-',i3,2x,6(f8.3,f6.1,i2))
         end do
  760    continue
      end if
c
c     stretch-torsion parameters
c
      if (kbt(1) .ne. '      ') then
         write (itxt,770)  formfeed,forcefield
  770    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,780)
  780    format (//,17x,'Stretch-Torsion Parameters',
     &           ///,24x,'Classes',14x,'KST1',8x,'KST2',8x,'KST3',/)
         do i = 1, maxnbt
            if (kbt(i) .eq. '      ')  goto 800
            k1 = number (kbt(i)(1:3))
            k2 = number (kbt(i)(4:6))
            write (itxt,790)  i,k1,k2,(btcon(j,i),j=1,3)
  790       format (14x,i5,5x,i3,'-',i3,6x,3f12.3)
         end do
  800    continue
      end if
c
c     atomic partial charge parameters
c
      exist = .false.
      do i = 1, maxtyp
         if (chg(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,810)  formfeed,forcefield
  810    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,820)
  820    format (//,17x,'Atomic Partial Charge Parameters',
     &           ///,28x,'Type',9x,'Partial Chg',/)
         k = 0
         do i = 1, maxtyp
            if (chg(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,830)  k,i,chg(i)
  830          format (16x,i5,7x,i3,6x,f12.3)
            end if
         end do
      end if
c
c     bond dipole moment parameters
c
      if (kd(1) .ne. '      ') then
         write (itxt,840)  formfeed,forcefield
  840    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,850)
  850    format (//,17x,'Bond Dipole Moment Parameters',
     &           ///,26x,'Types',9x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd
            if (kd(i) .eq. '      ')  goto 870
            k1 = number (kd(i)(1:3))
            k2 = number (kd(i)(4:6))
            write (itxt,860)  i,k1,k2,dpl(i),pos(i)
  860       format (14x,i5,5x,i3,'-',i3,6x,2f12.3)
         end do
  870    continue
      end if
c
c     bond dipole moment parameters for 5-membered rings
c
      if (kd5(1) .ne. '      ') then
         write (itxt,880)
  880    format (//,17x,'5-Membered Ring Bond Dipole Parameters',
     &           ///,26x,'Types',9x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd5
            if (kd5(i) .eq. '      ')  goto 900
            k1 = number (kd5(i)(1:3))
            k2 = number (kd5(i)(4:6))
            write (itxt,890)  i,k1,k2,dpl5(i),pos5(i)
  890       format (14x,i5,5x,i3,'-',i3,6x,2f12.3)
         end do
  900    continue
      end if
c
c     bond dipole moment parameters for 4-membered rings
c
      if (kd4(1) .ne. '      ') then
         write (itxt,910)
  910    format (//,17x,'4-Membered Ring Bond Dipole Parameters',
     &           ///,26x,'Types',9x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd4
            if (kd4(i) .eq. '      ')  goto 930
            k1 = number (kd4(i)(1:3))
            k2 = number (kd4(i)(4:6))
            write (itxt,920)  i,k1,k2,dpl4(i),pos4(i)
  920       format (14x,i5,5x,i3,'-',i3,6x,2f12.3)
         end do
  930    continue
      end if
c
c     bond dipole moment parameters for 3-membered rings
c
      if (kd3(1) .ne. '      ') then
         write (itxt,940)
  940    format (//,17x,'3-Membered Ring Bond Dipole Parameters',
     &           ///,26x,'Types',9x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd3
            if (kd3(i) .eq. '      ')  goto 960
            k1 = number (kd3(i)(1:3))
            k2 = number (kd3(i)(4:6))
            write (itxt,950)  i,k1,k2,dpl3(i),pos3(i)
  950       format (14x,i5,5x,i3,'-',i3,6x,2f12.3)
         end do
  960    continue
      end if
c
c     atomic multipole electrostatic parameters
c
      if (kmp(1) .ne. '         ') then
         write (itxt,970)  formfeed,forcefield
  970    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,980)
  980    format (//,17x,'Atomic Multipole Parameters',
     &           ///,17x,'Type',3x,'Axis Types',6x,'Frame',
     &              8x,'Multipoles (M-D-Q)',/)
         do i = 1, maxnmp
            if (kmp(i) .eq. '         ')  goto 1000
            k1 = number (kmp(i)(1:3))
            k2 = number (kmp(i)(4:6))
            k3 = number (kmp(i)(7:9))
            write (itxt,990)  i,k1,k2,k3,mpaxis(i),multip(1,i),
     &                        multip(2,i),multip(3,i),multip(4,i),
     &                        multip(5,i),multip(8,i),multip(9,i),
     &                        multip(11,i),multip(12,i),multip(13,i)
  990       format (5x,i5,7x,i3,5x,i3,2x,i3,5x,a8,10x,f10.5,
     &                 /,46x,3f10.5,/,46x,f10.5,
     &                 /,46x,2f10.5,/,46x,3f10.5)
         end do
 1000    continue
      end if
c
c     atomic dipole polarizability parameters
c
      exist = .false.
      do i = 1, maxtyp
         if (polr(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1010)  formfeed,forcefield
 1010    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1020)
 1020    format (//,17x,'Dipole Polarizability Parameters',
     &           ///,28x,'Type',14x,'Alpha',/)
         k = 0
         do i = 1, maxtyp
            if (polr(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1030)  k,i,polr(i)
 1030          format (16x,i5,7x,i3,6x,f12.3)
            end if
         end do
      end if
c
c     conjugated pisystem atom parameters
c
      exist = .false.
      do i = 1, maxclass
         if (ionize(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1040)  formfeed,forcefield
 1040    format (a1,//,17x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1050)
 1050    format (//,17x,'Conjugated Pisystem Atom Parameters',
     &           ///,21x,'Class',3x,'Electron',3x,
     &                 'Ionization',3x,'Repulsion',/)
         k = 0
         do i = 1, maxclass
            if (ionize(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1060)  k,i,electron(i),ionize(i),repulse(i)
 1060          format (10x,i5,7x,i3,f10.1,2x,2f12.3)
            end if
         end do
      end if
c
c     conjugated pisystem bond parameters
c
      if (kpi(1) .ne. '      ') then
         write (itxt,1070)
 1070    format (//,17x,'Conjugated Pisystem Bond Parameters',
     &           ///,20x,'Classes',8x,'d Force',4x,'d Length',/)
         do i = 1, maxnpi
            if (kpi(i) .eq. '      ')  goto 1090
            k1 = number (kpi(i)(1:3))
            k2 = number (kpi(i)(4:6))
            write (itxt,1080)  i,k1,k2,sslope(i),tslope(i)
 1080       format (10x,i5,5x,i3,'-',i3,3x,f12.3,f12.3)
         end do
 1090    continue
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
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prtseq  --  output of biopolymer sequence  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prtseq" writes out a biopolymer sequence to an external
c     disk file with 50 residues per line and distinct chains
c     separated by blank lines
c
c
      subroutine prtseq (iseq)
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'sequen.i'
      integer i,k,iseq,smax,smin
      integer size,start,stop
      character*1 letter
      character*60 seqfile
      logical opened
c
c
c     open output unit if not already done
c
      inquire (unit=iseq,opened=opened)
      if (.not. opened) then
         seqfile = filename(1:leng)//'.seq'
         call version (seqfile,'new')
         open (unit=iseq,file=seqfile,status='new')
      end if
c
c     write out a 1-letter code sequence file
c
      do i = 1, nchain
         letter = chnnam(i)
         start = ichain(1,i)
         stop = ichain(2,i)
         size = stop - start + 1
         smax = 0
         dowhile (smax .lt. size)
            smin = smax + 1
            smax = smax + 50
            smax = min(smax,size)
            if (i.ne.1 .and. smin.eq.1) then
               write (iseq,10)
   10          format ()
            end if
            write (iseq,20)  letter,smin,(seq(k+start-1),k=smin,smax)
   20       format (3x,a1,i6,1x,5(1x,10a1))
         end do
      end do
      if (.not. opened)  close (unit=iseq)
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
c     ##  subroutine prtxyz  --  output of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtxyz" writes out a set of Cartesian coordinates
c     to an external disk file, unit -ixyz-
c
c
      subroutine prtxyz (ixyz)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'titles.i'
      integer i,k,ixyz
      character*60 coordfile
      logical opened
c
c     output coordinates of all MM atoms, suitable for $TINXYZ.
c     use different formats to unit 4 (trajectory) or 6 (log file)
c
      if(ixyz.eq.6) write(ixyz,8000) 
      if(ixyz.eq.6) write(ixyz,8010) 
      if(ixyz.eq.6) write(ixyz,8000) 
c
      if (ltitle .eq. 0) then
         write (ixyz,10)  n
   10    format (i6)
      else
         write (ixyz,20)  n,ttitle(1:ltitle)
   20    format (i6,2x,a)
      end if
c
c     write the coordinates for each atom
c
      do i = 1, n
casa  -printxyz- begin
c         write(ixyz,8020) i,name(i),x(i),y(i),z(i),type(i),
c     &                    (i12(k,i),k=1,n12(i))
         write(ixyz,8020) i,xyzname(i),x(i),y(i),z(i),type(i),
     &                    (i12(k,i),k=1,n12(i))
      end do
c
      if(ixyz.eq.6) write(ixyz,8000) 
      return
c
 8000 format(1x,76(1h-))
 8010 format(11x,'Cartesian Coordinates of Atoms in Bulk Model (ANGS)')
 8020 format(i5,2x,a4,3f12.6,9i6)
c 8020 format(i5,2x,a3,3f12.6,9i5)
casa  -printxyz end
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
c     ##  subroutine qrfact  --  QR factorization of a matrix  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "qrfact" performs Householder transformations with column
c     pivoting (optional) to compute a QR factorization of the
c     m by n matrix a; the routine determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r; the Householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c               i - (1/u(k))*u*u(transpose)
c
c     where u has zeros in the first k-1 positions
c
c     arguments and variables :
c
c     m        positive integer input variable set to
c                the number of rows of a
c     n        positive integer input variable set to
c                the number of columns of a
c     a        m by n array; on input a contains the matrix for
c                which the QR factorization is to be computed; on
c                output the strict upper trapezoidal part contains
c                the strict upper trapezoidal part of r, the lower
c                trapezoidal part contains a factored form of q
c                (non-trivial elements of u vectors described above)
c     lda      positive integer input variable not less than m
c                which specifies the leading dimension of array a
c     pivot    logical input variable; if pivot is set true,
c                then column pivoting is enforced; if pivot is
c                set false, then no column pivoting is done
c     ipvt     integer output array which defines the permutation
c                matrix p such that a*p = q*r; column j of p is
c                column ipvt(j) of the identity matrix
c     rdiag    an output array of length n which contains
c                the diagonal elements of r
c
c
      subroutine qrfact (m,n,a,lda,pivot,ipvt,rdiag,work)
      implicit none
      integer m,n,lda,ipvt(n)
      integer i,j,k,jmax,minmn,itemp
      real*8 aknorm,temp,a(lda,n),rdiag(n),work(n)
      logical pivot
c
c
c     find the initial column norms and initialize some arrays
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + a(i,j)**2
         end do
         rdiag(j) = sqrt(temp)
         work(j) = rdiag(j)
         if (pivot)  ipvt(j) = j
      end do
c
c     reduce the matrix with Householder transformations
c
      minmn = min(m,n)
      do k = 1, minmn
c
c     bring the column of largest norm into the pivot position
c
         if (pivot) then
            jmax = k
            do j = k, n
               if (rdiag(j) .gt. rdiag(jmax))  jmax = j
            end do
            if (jmax .ne. k) then
               do i = 1, m
                  temp = a(i,k)
                  a(i,k) = a(i,jmax)
                  a(i,jmax) = temp
               end do
               rdiag(jmax) = rdiag(k)
               work(jmax) = work(k)
               itemp = ipvt(k)
               ipvt(k) = ipvt(jmax)
               ipvt(jmax) = itemp
            end if
         end if
c
c     compute the Householder transformation to reduce the
c     k-th column of a to a multiple of the k-th unit vector
c
         aknorm = 0.0d0
         do i = k, m
            aknorm = aknorm + a(i,k)**2
         end do
         aknorm = sqrt(aknorm)
         if (aknorm .ne. 0.0d0) then
            if (a(k,k) .lt. 0.0d0)  aknorm = -aknorm
            do i = k, m
               a(i,k) = a(i,k) / aknorm
            end do
            a(k,k) = a(k,k) + 1.0d0
c
c     apply the transformation to the remaining columns
c     and update the column norms
c
            if (n .ge. k+1) then
               do j = k+1, n
                  temp = 0.0d0
                  do i = k, m
                     temp = temp + a(i,k)*a(i,j)
                  end do
                  temp = temp / a(k,k)
                  do i = k, m
                     a(i,j) = a(i,j) - temp*a(i,k)
                  end do
                  if (pivot .and. rdiag(j).ne.0.0d0) then
                     temp = a(k,j) / rdiag(j)
                     if (abs(temp) .lt. 1.0d0) then
                        rdiag(j) = rdiag(j) * sqrt(1.0d0-temp**2)
                     else
                        temp = 0.0d0
                        do i = k+1, m
                           temp = temp + a(i,j)**2
                        end do
                        rdiag(j) = sqrt(temp)
                        work(j) = rdiag(j)
                     end if
                  end if
               end do
            end if
         end if
         rdiag(k) = -aknorm
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine qrsolv  --  least squares solution of QR factor  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "qrsolv" solves a*x=b and d*x=0 in the least squares sense;
c     normally used in combination with routine "qrfact" to solve
c     least squares problems
c
c     arguments and variables :
c
c     n        number of rows and columns in the matrix r
c     r        an n by n array containing the upper triangular
c                matrix r; on output the full triangle is unaltered,
c                and the strict lower triangle contains the transpose
c                of the strict upper triangular matrix s
c     ldr      leading dimension of r exactly as specified in
c                the dimension statement of the calling program
c     ipvt     vector of length n which defines the permutation
c                matrix p such that a*p = q*r; column j of p is
c                column ipvt(j) of the identity matrix
c     diag     vector of length n containing the diagonal elements
c                of the matrix d
c     qtb      vector of length n containing the first n elements
c                of the vector q(transpose)*b
c     x        vector of length n containing the least squares
c                solution of the systems a*x = b, d*x = 0
c     sdiag    vector of length n containing the diagonal elements
c                of the upper triangular matrix s
c
c
      subroutine qrsolv (n,r,ldr,ipvt,diag,qtb,x,sdiag,work)
      implicit none
      integer i,j,k,jj,nsing,n,ldr,ipvt(n)
      real*8 r(ldr,n),diag(n),qtb(n)
      real*8 x(n),sdiag(n),work(n)
      real*8 sine,cosine,tangent,cotangent,qtbpj,temp
c
c
c     copy r and (q transpose)*b to preserve input and
c     initialize s; in particular, save the diagonal
c     elements of r in x
c
      do j = 1, n-1
         do k = j+1, n
            r(k,j) = r(j,k)
         end do
      end do
      do j = 1, n
         x(j) = r(j,j)
         work(j) = qtb(j)
      end do
c
c     eliminate the diagonal matrix d using a Givens rotation
c
      do j = 1, n
c
c     prepare the row of d to be eliminated, locating
c     the diagonal element using p from the QR factorization
c
         jj = ipvt(j)
         if (diag(jj) .ne. 0.0d0) then
            do k = j, n
               sdiag(k) = 0.0d0
            end do
            sdiag(j) = diag(jj)
c
c     the transformations to eliminate the row of d modify
c     only a single element of (q transpose)*b beyond the
c     first n, which is initially zero
c
            qtbpj = 0.0d0
            do k = j, n
c
c     determine a Givens rotation which eliminates the
c     appropriate element in the current row of d
c
               if (sdiag(k) .ne. 0.0d0) then
                  if (abs(r(k,k)) .lt. abs(sdiag(k))) then
                     cotangent = r(k,k) / sdiag(k)
                     sine = 0.5d0 / sqrt(0.25d0+0.25d0*cotangent**2)
                     cosine = sine * cotangent
                  else
                     tangent = sdiag(k) / r(k,k)
                     cosine = 0.5d0 / sqrt(0.25d0+0.25d0*tangent**2)
                     sine = cosine * tangent
                  end if
c
c     compute the modified diagonal element of r
c     and the modified element of ((q transpose)*b,0)
c
                  r(k,k) = cosine*r(k,k) + sine*sdiag(k)
                  temp = cosine*work(k) + sine*qtbpj
                  qtbpj = -sine*work(k) + cosine*qtbpj
                  work(k) = temp
c
c     accumulate the tranformation in the row of s
c
                  if (n .ge. k+1) then
                     do i = k+1, n
                        temp = cosine*r(i,k) + sine*sdiag(i)
                        sdiag(i) = -sine*r(i,k) + cosine*sdiag(i)
                        r(i,k) = temp
                     end do
                  end if
               end if
            end do
         end if
c
c     store the diagonal element of s and restore
c     the corresponding diagonal element of r
c
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
      end do
c
c     solve the triangular system for z; if the system
c     is singular, then obtain a least squares solution
c
      nsing = n
      do j = 1, n
         if (sdiag(j).eq.0.0d0 .and. nsing.eq.n)  nsing = j - 1
         if (nsing .lt. n)  work(j) = 0.0d0
      end do
      if (nsing .ge. 1) then
         do k = 1, nsing
            j = nsing - k + 1
            temp = 0.0d0
            if (nsing .ge. j+1) then
               do i = j+1, nsing
                  temp = temp + r(i,j)*work(i)
               end do
            end if
            work(j) = (work(j)-temp) / sdiag(j)
         end do
      end if
c
c     permute the components of z back to components of x
c
      do j = 1, n
         k = ipvt(j)
         x(k) = work(j)
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine quatfit  --  quaternion superposition of coords  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "quatfit" uses a quaternion based method to achieve the best
c     fit superposition of two sets of coordinates
c
c     literature reference:
c
c     S. J. Kearsley, "An Algorithm for the Simultaneous Superposition
c     of a Structural Series", Journal of Computational Chemistry,
c     11, 1187-1192 (1990)
c
c     adapted from an original program written by David J. Heisterberg,
c     Ohio Supercomputer Center, Columbus, OH
c
c
      subroutine quatfit (n1,x1,y1,z1,n2,x2,y2,z2)
      implicit none
      include 'sizes.i'
      include 'align.i'
      integer i,i1,i2,n1,n2
      real*8 weight,xrot,yrot,zrot
      real*8 xxyx,xxyy,xxyz,xyyx,xyyy
      real*8 xyyz,xzyx,xzyy,xzyz
      real*8 rot(3,3),temp1(4),temp2(4)
      real*8 q(4),d(4),c(4,4),v(4,4)
      real*8 x1(maxatm),y1(maxatm),z1(maxatm)
      real*8 x2(maxatm),y2(maxatm),z2(maxatm)
c
c
c     build the upper triangle of the quadratic form matrix
c
      xxyx = 0.0d0
      xxyy = 0.0d0
      xxyz = 0.0d0
      xyyx = 0.0d0
      xyyy = 0.0d0
      xyyz = 0.0d0
      xzyx = 0.0d0
      xzyy = 0.0d0
      xzyz = 0.0d0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         weight = wfit(i)
         xxyx = xxyx + weight*x1(i1)*x2(i2)
         xxyy = xxyy + weight*y1(i1)*x2(i2)
         xxyz = xxyz + weight*z1(i1)*x2(i2)
         xyyx = xyyx + weight*x1(i1)*y2(i2)
         xyyy = xyyy + weight*y1(i1)*y2(i2)
         xyyz = xyyz + weight*z1(i1)*y2(i2)
         xzyx = xzyx + weight*x1(i1)*z2(i2)
         xzyy = xzyy + weight*y1(i1)*z2(i2)
         xzyz = xzyz + weight*z1(i1)*z2(i2)
      end do
      c(1,1) = xxyx + xyyy + xzyz
      c(1,2) = xzyy - xyyz
      c(2,2) = xxyx - xyyy - xzyz
      c(1,3) = xxyz - xzyx
      c(2,3) = xxyy + xyyx
      c(3,3) = xyyy - xzyz - xxyx
      c(1,4) = xyyx - xxyy
      c(2,4) = xzyx + xxyz
      c(3,4) = xyyz + xzyy
      c(4,4) = xzyz - xxyx - xyyy
c
c     diagonalize the quadratic form matrix
c
      call tnk_jacobi (4,4,c,d,v,temp1,temp2)
c
c     extract the desired quaternion
c
      q(1) = v(1,4)
      q(2) = v(2,4)
      q(3) = v(3,4)
      q(4) = v(4,4)
c
c     assemble the rotation matrix that superimposes molecules
c
      rot(1,1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
      rot(2,1) = 2.0d0 * (q(2) * q(3) - q(1) * q(4))
      rot(3,1) = 2.0d0 * (q(2) * q(4) + q(1) * q(3))
      rot(1,2) = 2.0d0 * (q(3) * q(2) + q(1) * q(4))
      rot(2,2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
      rot(3,2) = 2.0d0 * (q(3) * q(4) - q(1) * q(2))
      rot(1,3) = 2.0d0 * (q(4) * q(2) - q(1) * q(3))
      rot(2,3) = 2.0d0 * (q(4) * q(3) + q(1) * q(2))
      rot(3,3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
c
c     rotate second molecule to best fit with first molecule
c
      do i = 1, n2
         xrot = x2(i)*rot(1,1) + y2(i)*rot(1,2) + z2(i)*rot(1,3)
         yrot = x2(i)*rot(2,1) + y2(i)*rot(2,2) + z2(i)*rot(2,3)
         zrot = x2(i)*rot(3,1) + y2(i)*rot(3,2) + z2(i)*rot(3,3)
         x2(i) = xrot
         y2(i) = yrot
         z2(i) = zrot
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
c     ##  function random  --  portable random number generator  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "random" generates a random number on [0,1] via a long
c     period generator due to L'Ecuyer with Bays-Durham shuffle
c
c     literature references:
c
c     P. L'Ecuyer, Communications of the ACM, 31, 742-774 (1988)
c
c     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
c     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
c     University Press, 1992, Section 7-1
c
c
      function random ()
      implicit none
      include 'sizes.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      integer im1,ia1,iq1,ir1
      integer im2,ia2,iq2,ir2
      integer big,ntable
      integer imm1,ndiv
      real*8 factor
      parameter (im1=2147483563)
      parameter (ia1=40014)
      parameter (iq1=53668)
      parameter (ir1=12211)
      parameter (im2=2147483399)
      parameter (ia2=40692)
      parameter (iq2=52774)
      parameter (ir2=3791)
      parameter (big=141803398)
      parameter (ntable=32)
      parameter (imm1=im1-1)
      parameter (ndiv=1+imm1/ntable)
      parameter (factor=1.0d0/im1)
      integer i,k,next,seed,seed2
      integer iy,itable(ntable)
      integer year,month,day
      integer hour,minute,second
      real*8 random
      character*20 keyword
      character*80 record,string
      logical initial
      save initial,seed,seed2,iy,itable
      data initial  / .true. /
c
c
c     random number seed is first set to a big number,
c     then incremented by the seconds elapsed this decade
c
      if (initial) then
         initial = .false.
         seed = big
         call calendar (year,month,day,hour,minute,second)
         year = mod(year,10)
         seed = seed + 32140800*year + 2678400*(month-1)
         seed = seed + 86400*(day-1) + 3600*hour
         seed = seed + 60*minute + second
c
c     search the keywords for a random number seed
c
         do i = 1, nkey
            next = 1
            record = keyline(i)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:11) .eq. 'RANDOMSEED ') then
               string = record(next:80)
               read (string,*,err=10)  seed
               seed = max(1,seed)
            end if
   10       continue
         end do
c
c     print the value used for the random number seed
c
         if (verbose) then
            write (iout,20)  seed
   20       format (/,' RANDOM  --  Initialized with SEED of',i12)
         end if
c
c     warm up and then load the shuffling table
c
         seed2 = seed
         do i = ntable+8, 1, -1
            k = seed / iq1
            seed = ia1 * (seed-k*iq1) - k*ir1
            if (seed .lt. 0)  seed = seed + im1
            if (i .le. ntable)  itable(i) = seed
         end do
         iy = itable(1)
      end if
c
c     get a new random number value each call
c
      k = seed / iq1
      seed = ia1*(seed-k*iq1) - k*ir1
      if (seed .lt. 0)  seed = seed + im1
      k = seed2 / iq2
      seed2 = ia2*(seed2-k*iq2) - k*ir2
      if (seed2 .lt. 0)  seed2 = seed2 + im2
      i = 1 + iy/ndiv
      iy = itable(i) - seed2
      itable(i) = seed
      if (iy .lt. 1)  iy = iy + imm1
      random = factor * iy
c
c     print the value of the current random number
c
c     if (debug) then
c        write (iout,30)  random
c  30    format (' RANDOM  --  The Random Number Value is',f12.8)
c     end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function normal  --  random number from normal curve  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "normal" generates a random number from a normal Gaussian
c     distribution with a mean of zero and a variance of one
c
c
      function normal ()
      implicit none
      include 'inform.i'
      include 'iounit.i'
      real*8 random,v1,v2,rsq
      real*8 factor,store,normal
      logical compute
      save compute,store
      data compute  / .true. /
c
c
c     get a pair of random values from the distribution
c
      if (compute) then
   10    continue
         v1 = 2.0d0 * random () - 1.0d0
         v2 = 2.0d0 * random () - 1.0d0
         rsq = v1**2 + v2**2
         if (rsq .ge. 1.0d0)  goto 10
         factor = sqrt(-2.0d0*log(rsq)/rsq)
         store = v1 * factor
         normal = v2 * factor
         compute = .false.
c
c     use the second random value computed at the last call
c
      else
         normal = store
         compute = .true.
      end if
c
c     print the value of the current random number
c
c     if (debug) then
c        write (iout,20)  normal
c  20    format (' NORMAL  --  The Random Number Value is',f12.8)
c     end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ranvec  --  unit vector in random direction  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ranvec" generates a unit vector in 3-dimensional
c     space with uniformly distributed random orientation
c
c     literature references:
c
c     G. Marsaglia, Ann. Math. Stat., 43, 645 (1972)
c
c     R. C. Rapaport, The Art of Molecular Dynamics Simulation,
c     Cambridge University Press, 1995, Appendix A4
c
c
      subroutine ranvec (vector)
      implicit none
      include 'inform.i'
      include 'iounit.i'
      real*8 random,vector(3)
      real*8 x,y,s
c
c
c     get a pair of appropriate components in the plane
c
      s = 2.0d0
      dowhile (s .ge. 1.0d0)
         x = 2.0d0 * random () - 1.0d0
         y = 2.0d0 * random () - 1.0d0
         s = x**2 + y**2
      end do
c
c     construct the 3-dimensional random unit vector
c
      vector(3) = 1.0d0 - 2.0d0*s
      s = 2.0d0 * sqrt(1.0d0 - s)
      vector(2) = s * y
      vector(1) = s * x
c
c     print the components of the random unit vector
c
c     if (debug) then
c        write (iout,10)  vector(1),vector(2),vector(3)
c  10    format (' RANVEC  --  The Random Vector is',3f10.4)
c     end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rattle  --  apply rattle position constraints  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rattle" implements the first portion of the rattle algorithm
c     by correcting atomic positions and half-step velocities to
c     maintain constrained interatomic distances
c
c     Literature Reference:
c
c     H. C. Andersen, "Rattle: A Velocity Version of the Shake
c     Algorithm for Molecular Dynamics Calculations", Journal of
c     Computational Physics, 52, 24-34 (1983)
c
c
      subroutine rattle (dt,xold,yold,zold)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'usage.i'
      integer i,ia,ib,niter,maxiter
      real*8 dt,xr,yr,zr,xo,yo,zo
      real*8 eps,dist2,delta,rma,rmb
      real*8 dot,term,xterm,yterm,zterm
      real*8 xold(maxatm),yold(maxatm),zold(maxatm)
      logical done,moved(maxatm),update(maxatm)
c
c
c     initialize the lists of atoms previously corrected
c
      do i = 1, n
         if (use(i)) then
            moved(i) = .true.
         else
            moved(i) = .false.
         end if
         update(i) = .false.
      end do
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 100
      niter = 0
      done = .false.
      eps = 0.000001d0
c
c     apply the rattle algorithm to correct the atom
c     positions and the half-step velocity values
c
      dowhile (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         do i = 1, nrat
            ia = irat(1,i)
            ib = irat(2,i)
            if (moved(ia) .or. moved(ib)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               if (ratimage(i))  call image (xr,yr,zr,0)
               dist2 = xr**2 + yr**2 + zr**2
               delta = krat(i)**2 - dist2
               if (abs(delta) .gt. eps) then
                  done = .false.
                  update(ia) = .true.
                  update(ib) = .true.
                  xo = xold(ib) - xold(ia)
                  yo = yold(ib) - yold(ia)
                  zo = zold(ib) - zold(ia)
                  if (ratimage(i))  call image (xo,yo,zo,0)
                  dot = xr*xo + yr*yo + zr*zo
                  rma = 1.0d0 / mass(ia)
                  rmb = 1.0d0 / mass(ib)
                  term = delta / (2.0d0 * (rma+rmb) * dot)
                  xterm = xo * term
                  yterm = yo * term
                  zterm = zo * term
                  x(ia) =  x(ia) - xterm * rma
                  y(ia) =  y(ia) - yterm * rma
                  z(ia) =  z(ia) - zterm * rma
                  x(ib) =  x(ib) + xterm * rmb
                  y(ib) =  y(ib) + yterm * rmb
                  z(ib) =  z(ib) + zterm * rmb
                  rma = rma / dt
                  rmb = rmb / dt
                  v(1,ia) = v(1,ia) - xterm * rma
                  v(2,ia) = v(2,ia) - yterm * rma
                  v(3,ia) = v(3,ia) - zterm * rma
                  v(1,ib) = v(1,ib) + xterm * rmb
                  v(2,ib) = v(2,ib) + yterm * rmb
                  v(3,ib) = v(3,ib) + zterm * rmb
               end if
            end if
         end do
         do i = 1, n
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
c
c     write information on the number of iterations needed
c
      if (niter .eq. maxiter) then
         write (iout,10)
   10    format (/,' RATTLE  --  Warning, Position Constraints',
     &              ' not Satisfied')
         call prterr
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' RATTLE   --  Position Constraints met at',i6,
     &              ' Iterations')
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rattle2  --  apply rattle velocity constraints  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rattle2" implements the second portion of the rattle algorithm
c     by correcting the full-step velocities in order to maintain
c     constrained interatomic distances
c
c
      subroutine rattle2 (dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moldyn.i'
      include 'shake.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ia,ib,niter,maxiter
      real*8 xr,yr,zr,xv,yv,zv
      real*8 dt,eps,rma,rmb,dot,vterm
      real*8 term,xterm,yterm,zterm
      logical done,moved(maxatm),update(maxatm)
c
c
c     initialize the lists of atoms previously corrected
c
      do i = 1, n
         if (use(i)) then
            moved(i) = .true.
         else
            moved(i) = .false.
         end if
         update(i) = .false.
      end do
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 100
      niter = 0
      done = .false.
      eps = 0.000001d0 / dt
      vterm = 2.0d0 / (dt * convert)
c
c     apply the rattle algorithm to correct the velocities
c
      dowhile (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         do i = 1, nrat
            ia = irat(1,i)
            ib = irat(2,i)
            if (moved(ia) .or. moved(ib)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               if (ratimage(i))  call image (xr,yr,zr,0)
               xv = v(1,ib) - v(1,ia)
               yv = v(2,ib) - v(2,ia)
               zv = v(3,ib) - v(3,ia)
               dot = xr*xv + yr*yv + zr*zv
               rma = 1.0d0 / mass(ia)
               rmb = 1.0d0 / mass(ib)
               term = -dot / ((rma+rmb) * krat(i)**2)
               if (abs(term) .gt. eps) then
                  done = .false.
                  update(ia) = .true.
                  update(ib) = .true.
                  xterm = xr * term
                  yterm = yr * term
                  zterm = zr * term
                  v(1,ia) = v(1,ia) - xterm * rma
                  v(2,ia) = v(2,ia) - yterm * rma
                  v(3,ia) = v(3,ia) - zterm * rma
                  v(1,ib) = v(1,ib) + xterm * rmb
                  v(2,ib) = v(2,ib) + yterm * rmb
                  v(3,ib) = v(3,ib) + zterm * rmb
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx - xr * xterm * vterm
                     viry = viry - yr * yterm * vterm
                     virz = virz - zr * zterm * vterm
                  end if
               end if
            end if
         end do
         do i = 1, n
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
c
c     write information on the number of iterations needed
c
      if (niter .eq. maxiter) then
         write (iout,10)
   10    format (/,' RATTLE2  --  Warning, Velocity Constraints',
     &              ' not Satisfied')
         call prterr
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' RATTLE2  --  Velocity Constraints met at',i6,
     &              ' Iterations')
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readdyn  --  input of MD restart information  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readdyn" get the positions, velocities and accelerations
c     for a molecular dynamics restart from an external disk file
c
c
      subroutine readdyn (idyn)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'files.i'
      include 'iounit.i'
      include 'moldyn.i'
      integer i,idyn,ndyn
      character*60 dynfile
      character*80 record
      logical exist,opened,quit
c
c
c     open the input file if it has not already been done
c
      inquire (unit=idyn,opened=opened)
      if (.not. opened) then
         dynfile = filename(1:leng)//'.dyn'
         call version (dynfile,'old')
         inquire (file=dynfile,exist=exist)
         if (exist) then
            open (unit=idyn,file=dynfile,status='old')
            rewind (unit=idyn)
         else
            write (iout,10)
   10       format (/,' READDYN  --  Unable to Open Dynamics',
     &                 ' Restart File')
            call fatal
         end if
      end if
c
c     initialize error handling during reading of the file
c
      i = 0
      quit = .true.
c
c     get the number of atoms and check for consistency
c
      read (idyn,20)
   20 format ()
      read (idyn,30)  record
   30 format (a80)
      read (record,*,err=160,end=160)  ndyn
      if (ndyn .ne. n) then
         write (iout,40)
   40    format (/,' READDYN  --  Restart File has Incorrect',
     &              ' Number at Atoms')
         call fatal
      end if
c
c     get the periodic box edge lengths and angles
c
      read (idyn,50)
   50 format ()
      read (idyn,60)  record
   60 format (a80)
      read (record,*,err=160,end=160)  xbox,ybox,zbox
      read (idyn,70)  record
   70 format (a80)
      read (record,*,err=160,end=160)  alpha,beta,gamma
      read (idyn,80)
   80 format ()
c
c     get the atomic positions, velocities and accelerations
c
      quit = .true.
      do i = 1, n
         read (idyn,90)  record
   90    format (a80)
         read (record,*,err=160,end=160)  x(i),y(i),z(i)
      end do
      read (idyn,100)
  100 format ()
      do i = 1, n
         read (idyn,110)  record
  110    format (a80)
         read (record,*,err=160,end=160)  v(1,i),v(2,i),v(3,i)
      end do
      read (idyn,120)
  120 format ()
      do i = 1, n
         read (idyn,130)  record
  130    format (a80)
         read (record,*,err=160,end=160)  a(1,i),a(2,i),a(3,i)
      end do
      read (idyn,140)
  140 format ()
      do i = 1, n
         read (idyn,150)  record
  150    format (a80)
         read (record,*,err=160,end=160)  a_old(1,i),a_old(2,i),
     &                                    a_old(3,i)
      end do
      quit = .false.
  160 continue
      if (.not. opened)  close (unit=idyn)
c
c     report any error in reading the dynamics restart file
c
      if (quit) then
         write (iout,170)  i
  170    format (/,' READDYN  --  Error in Dynamics Restart',
     &              ' File at Atom',i6)
         call fatal
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
c     ##  subroutine readint  --  input of internal coordinates  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "readint" gets a set of Z-matrix internal
c     coordinates from an external file
c
c
      subroutine readint (izmt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      include 'titles.i'
      include 'zclose.i'
      include 'zcoord.i'
      integer i,j,izmt
      integer next,first,last
      integer nexttext,trimtext
      character*60 intfile
      character*80 record,string
      logical exist,opened,quit
c
c
c     open the input file if it has not already been done
c
      inquire (unit=izmt,opened=opened)
      if (.not. opened) then
         intfile = filename(1:leng)//'.int'
         open (unit=izmt,file=intfile,status='old')
         rewind (unit=izmt)
         call version (intfile,'old')
         inquire (file=intfile,exist=exist)
         if (exist) then
            open (unit=izmt,file=intfile,status='old')
            rewind (unit=izmt)
         else
            write (iout,10)
   10       format (/,' READINT  --  Unable to Open Internal',
     &                 ' Coordinates File')
            call fatal
         end if
      end if
c
c     read the title line and get the number of atoms
c
      quit = .true.
      i = 0
      read (izmt,20)  record
   20 format (a80)
      next = 1
      call gettext (record,string,next)
      read (string,*,err=60,end=60)  n
c
c     extract the title and determine its length
c
      string = record(next:80)
      first = nexttext (string)
      last = trimtext (string)
      if (last .eq. 0) then
         ttitle = ' '
         ltitle = 0
      else
         ttitle = string(first:last)
         ltitle = trimtext (ttitle)
      end if
c
c     check for too many total atoms in the file
c
      if (n .gt. maxatm) then
         write (iout,30)  maxatm
   30    format (' READINT  --  The Maximum of',i6,' Atoms',
     &           ' has been Exceeded')
         call fatal
      end if
c
c     read the coordinates and connectivities for each atom
c
      do i = 1, n
         next = 1
         read (izmt,40,end=60)  record
   40    format (a80)
         read (record,*,err=60,end=60)  tag(i)
         call getword (record,name(i),next)
         string = record(next:80)
         read (string,*,err=60,end=50)  type(i),iz(1,i),zbond(i),
     &                                  iz(2,i),zang(i),iz(3,i),
     &                                  ztors(i),iz(4,i)
   50    continue
      end do
      quit = .false.
   60 continue
      if (.not. opened)  close (unit=izmt)
c
c     an error occurred in reading the Z-matrix coordinates
c
      if (quit) then
         write (iout,70)  i
   70    format (' READZ  --  Error in Z-Matrix File at Atom',i6)
         call fatal
      end if
c
c     read in any additional bonds to be added or deleted
c
      nadd = 0
      ndel = 0
      read (izmt,80,end=120)
   80 format ()
      do i = 1, maxatm
         read (izmt,90,end=120)  record
   90    format (a80)
         if (nexttext(record) .eq. 0)  goto 100
         read (record,*,err=100)  (iadd(j,i),j=1,2)
         nadd = nadd + 1
      end do
  100 continue
      do i = 1, maxatm
         read (izmt,110,end=120)  record
  110    format (a80)
         if (nexttext(record) .eq. 0)  goto 120
         read (record,*,err=120)  (idel(j,i),j=1,2)
         ndel = ndel + 1
      end do
  120 continue
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
c     ##  subroutine readmol2  --  input of a Sybyl MOL2 file  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "readmol2" gets a set of Sybyl MOL2 coordinates
c     from an external disk file
c
c
      subroutine readmol2 (isyb)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'iounit.i'
      include 'titles.i'
      integer i,j,k,m,ia,ib,isyb
      integer nbond,number
      integer next,trimtext
      character*10 atmnam
      character*60 sybylfile
      character*80 record,string
      logical exist,opened
c
c
c     open the input file if it has not already been done
c
      inquire (unit=isyb,opened=opened)
      if (.not. opened) then
         sybylfile = filename(1:leng)//'.mol2'
         call version (sybylfile,'old')
         inquire (file=sybylfile,exist=exist)
         if (exist) then
            open (unit=isyb,file=sybylfile,status='old')
            rewind (unit=isyb)
         else
            write (iout,10)
   10       format (/,' READMOL2  --  Unable to Open TRIPOS',
     &                 ' Sybyl MOL2 File')
            call fatal
         end if
      end if
c
c     get title line and get the number of atoms and bonds
c
      do i = 1, 1000000
         read (isyb,20)  record
   20    format (a80)
         next = 1
         call gettext (record,string,next)
         call upcase (string)
         if (string .eq. '@<TRIPOS>MOLECULE') then
            read (isyb,30)  ttitle
   30       format (a60)
            ltitle = trimtext (ttitle)
            read (isyb,40)  record
   40       format (a80)
            read (record,*)  n,nbond
            goto 50
         end if
      end do
   50 continue
c
c     read the atom names and coordinates
c
      do i = 1, 1000000
         read (isyb,60)  record
   60    format (a80)
         next = 1
         call gettext (record,string,next)
         call upcase (string)
         if (string .eq. '@<TRIPOS>ATOM') then
            do j = 1, n
               read (isyb,70)  record
   70          format (a80)
               read (record,*)  number
               next = 1
               call getword (record,atmnam,next)
               read (record(next:80),*)  x(j),y(j),z(j)
               call getword (record,atmnam,next)
               name(j) = atmnam(1:3)
               do k = 1, 3
                  if (atmnam(k:k) .eq. '.') then
                     do m = k, 3
                        name(j)(m:m) = ' '
                     end do
                  end if
               end do
               type(j) = 0
            end do
            goto 80
         end if
      end do
   80 continue
c
c     read the bond list to get attached atom lists
c
      do i = 1, 1000000
         read (isyb,90)  record
   90    format (a80)
         next = 1
         call gettext (record,string,next)
         call upcase (string)
         if (string .eq. '@<TRIPOS>BOND') then
            do j = 1, nbond
               read (isyb,100)  record
  100          format (a80)
               read (record,*)  number,ia,ib
               n12(ia) = n12(ia) + 1
               i12(n12(ia),ia) = ib
               n12(ib) = n12(ib) + 1
               i12(n12(ib),ib) = ia
            end do
            goto 110
         end if
      end do
  110 continue
c
c     for each atom, sort its list of attached atoms
c
      do i = 1, n
         call sort (n12(i),i12(1,i))
      end do
      if (.not. opened)  close (unit=isyb)
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
c     ##  subroutine readpdb  --  input of Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readpdb" gets a set of Protein Data Bank coordinates
c     from an external disk file
c
c
      subroutine readpdb (ipdb)
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'iounit.i'
      include 'pdb.i'
      include 'sequen.i'
      include 'titles.i'
      integer i,ipdb,next,resnumb
      integer index,trimtext
      integer nalt,nins,length
      integer serial,residue,reslast
      real*8 xx,yy,zz
      character*1 altloc,altsym,altlast
      character*1 chain,chainsym,chnlast
      character*1 insert,inslast
      character*1 letter,chnatm(maxatm)
      character*3 resname,namelast
      character*4 atmname,word
      character*6 remark
      character*20 alttyp,chaintyp
      character*20 instyp,instemp
      character*60 pdbfile
      character*80 record,string
      logical exist,opened,model
      intrinsic index
c
c
c     open the input file if it has not already been done
c
      inquire (unit=ipdb,opened=opened)
      if (.not. opened) then
         pdbfile = filename(1:leng)//'.pdb'
         call version (pdbfile,'old')
         inquire (file=pdbfile,exist=exist)
         if (exist) then
            open (unit=ipdb,file=pdbfile,status='old')
            rewind (unit=ipdb)
         else
            write (iout,10)
   10       format (/,' READPDB  --  Unable to Brookhaven Protein',
     &                 ' Databank File')
            call fatal
         end if
      end if
c
c     initialize the model flag, residue counter and residue name
c
      model = .false.
      resnumb = 0
      reslast = 10000
      namelast = '   '
c
c     initialize lists of alternate sites, chains and insertions
c
      nalt = 0
      nchain = 0
      nins = 0
      altlast = '#'
      chnlast = '#'
      inslast = '#'
      alttyp = '                    '
      chaintyp = '####################'
      instyp = '                    '
c
c     extract header information from the Protein Data Bank file
c
      dowhile (.true.)
         read (ipdb,20,end=30) record
   20    format (a80)
         call upcase (record(1:6))
         remark = record(1:6)
         if (remark .eq. 'HEADER') then
            ttitle = record(11:70)
            ltitle = trimtext (ttitle)
            goto 30
         end if
      end do
   30 continue
c
c     scan for alternate locations, multiple chains and inserts
c
      rewind (unit=ipdb)
      dowhile (.true.)
         read (ipdb,40,end=60)  record
   40    format (a80)
         call upcase (record)
         remark = record(1:6)
         string = record(7:80)
         if (remark.eq.'ATOM  ' .or. remark.eq.'HETATM') then
            read (string,50)  altloc,chain,insert
   50       format (10x,a1,4x,a1,4x,a1)
            if (altloc .ne. altlast) then
               if (index(alttyp,altloc) .eq. 0) then
                  nalt = nalt + 1
                  alttyp(nalt:nalt) = altloc
                  altlast = altloc
               end if
            end if
            if (chain .ne. chnlast) then
               if (index(chaintyp,chain) .eq. 0) then
                  nchain = nchain + 1
                  chaintyp(nchain:nchain) = chain
                  chnlast = chain
               end if
            end if
            if (insert .ne. inslast) then
               if (index(instyp,insert) .eq. 0) then
                  nins = nins + 1
                  instyp(nins:nins) = insert
                  inslast = insert
               end if
            end if
         end if
      end do
   60 continue
c
c     find out which set of alternate locations will be used
c
      altsym = ' '
      if (nalt .gt. 0) then
         string(1:3) = '['//alttyp(1:1)//']'
         length = 3
         do i = 2, nalt
            string = string(1:length)//' '//alttyp(i:i)
            length = length + 2
         end do
         write (iout,70)  string(1:length)
   70    format (/,' Choose a Set of Alternate Atom Locations',
     &              ' from (',a,') :  ',$)
         read (input,80)  record
   80    format (a80)
         next = 1
         call gettext (record,altsym,next)
         if (altsym .eq. ' ')  altsym = alttyp(1:1)
         call upcase (altsym)
      end if
c
c     find out which, if any, insert records will be used
c
      if (nins .gt. 0) then
         instemp = '                    '
         string(1:1) = instyp(1:1)
         length = 1
         do i = 2, nins
            string = string(1:length)//' '//instyp(i:i)
            length = length + 2
         end do
         string = string(1:length)//' [ALL] NONE'
         length = length + 11
         write (iout,90)  string(1:length)
   90    format (/,' Enter the Insert Records to include',
     &              ' (',a,') :  ',$)
         read (input,100)  instemp
  100    format (a20)
         call upcase (instemp)
         next = 1
         call getword (instemp,word,next)
         if (word.eq.'    ' .or. word.eq.'ALL ') then
            continue
         else if (word .eq. 'NONE') then
            instyp = '                    '
         else
            instyp = instemp
         end if
      end if
c
c     process every atom of each chain from the Protein Data Bank file
c
      do i = 1, nchain
         rewind (unit=ipdb)
         chainsym = chaintyp(i:i)
         dowhile (.true.)
            read (ipdb,110,end=170) record
  110       format (a80)
            call upcase (record)
            remark = record(1:6)
            if (remark .eq. 'ATOM  ') then
               string = record(7:80)
               read (string,120)  serial,atmname,altloc,resname,
     &                            chain,residue,insert,xx,yy,zz
  120          format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)
               if (chain .ne. chainsym)  goto 130
               if (altloc.ne.' ' .and. altloc.ne.altsym)  goto 130
               if (insert.ne.' ' .and. index(instyp,insert).eq.0)
     &            goto 130
               if (residue.ne.reslast .or. resname.ne.namelast
     &                     .or. insert.ne.inslast) then
                  resnumb = resnumb + 1
                  reslast = residue
                  namelast = resname
                  inslast = insert
               end if
               call fixpdb (resname,atmname)
               npdb = npdb + 1
               xpdb(npdb) = xx
               ypdb(npdb) = yy
               zpdb(npdb) = zz
               pdbtyp(npdb) = remark
               atmnam(npdb) = atmname
               resnam(npdb) = resname
               resnum(npdb) = resnumb
               chnatm(npdb) = chain
  130          continue
            else if (remark .eq. 'HETATM') then
               string = record(7:80)
               read (string,140)  serial,atmname,altloc,resname,
     &                            chain,residue,insert,xx,yy,zz
  140          format (i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3)
               if (chain .ne. chainsym)  goto 150
               if (altloc.ne.' ' .and. altloc.ne.altsym)  goto 150
               if (insert.ne.' ' .and. index(instyp,insert).eq.0)
     &            goto 150
               call fixpdb (resname,atmname)
               npdb = npdb + 1
               xpdb(npdb) = xx
               ypdb(npdb) = yy
               zpdb(npdb) = zz
               pdbtyp(npdb) = remark
               atmnam(npdb) = atmname
               resnam(npdb) = resname
               resnum(npdb) = 0
               chnatm(npdb) = chain
  150          continue
            else if (remark .eq. 'MODEL ') then
               if (model) then
                  write (iout,160)
  160             format (/,' READPDB  --  File contains Multiple',
     &                       ' Models; First one Used')
                  goto 170
               else
                  model = .true.
               end if
            end if
         end do
  170    continue
      end do
c
c     set the total sequence length and chain termini information
c
      nseq = npdb
      nchain = 0
      chnlast = '#'
      do i = 1, npdb
         if (pdbtyp(i) .eq. 'ATOM  ') then
            letter = chnatm(i)
            if (letter .ne. chnlast) then
               nchain = nchain + 1
               ichain(1,nchain) = resnum(i)
               chnnam(nchain) = letter
               chnlast = letter
            else
               ichain(2,nchain) = resnum(i)
            end if
         end if
      end do
      if (.not. opened)  close (unit=ipdb)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine fixpdb  --  correct nonstandard PDB atom names  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fixpdb" corrects problems with the AMBER and CHARMM/XPLOR
c     PDB files by converting atom names to the PDB standards
c
c
      subroutine fixpdb (resname,atmname)
      implicit none
      character*3 resname
      character*4 atmname
c
c
c     convert any generically used unusual names
c
      if (atmname .eq. ' HN ')  atmname = ' H  '
c
c     convert any unusual names in terminal residues
c
      if (atmname .eq. ' HN1')  atmname = '1H  '
      if (atmname .eq. ' HN2')  atmname = '2H  '
      if (atmname .eq. ' HN3')  atmname = '3H  '
      if (atmname .eq. ' OT1')  atmname = ' O  '
      if (atmname .eq. 'OCT1')  atmname = ' O  '
      if (atmname .eq. ' OT ')  atmname = ' OXT'
      if (atmname .eq. ' OT2')  atmname = ' OXT'
      if (atmname .eq. 'OCT2')  atmname = ' OXT'
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         if (atmname .eq. ' HA1')  atmname = '1HA '
         if (atmname .eq. ' HA2')  atmname = '2HA '
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HB3')  atmname = '3HB '
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         if (atmname .eq. 'HG11')  atmname = '1HG1'
         if (atmname .eq. 'HG12')  atmname = '2HG1'
         if (atmname .eq. 'HG13')  atmname = '3HG1'
         if (atmname .eq. 'HG21')  atmname = '1HG2'
         if (atmname .eq. 'HG22')  atmname = '2HG2'
         if (atmname .eq. 'HG23')  atmname = '3HG2'
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. 'HD11')  atmname = '1HD1'
         if (atmname .eq. 'HD12')  atmname = '2HD1'
         if (atmname .eq. 'HD13')  atmname = '3HD1'
         if (atmname .eq. 'HD21')  atmname = '1HD2'
         if (atmname .eq. 'HD22')  atmname = '2HD2'
         if (atmname .eq. 'HD23')  atmname = '3HD2'
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         if (atmname .eq. ' CD ')  atmname = ' CD1'
         if (atmname .eq. 'HG11')  atmname = '1HG1'
         if (atmname .eq. 'HG12')  atmname = '2HG1'
         if (atmname .eq. 'HG21')  atmname = '1HG2'
         if (atmname .eq. 'HG22')  atmname = '2HG2'
         if (atmname .eq. 'HG23')  atmname = '3HG2'
         if (atmname .eq. 'HD11')  atmname = '1HD1'
         if (atmname .eq. ' HD1')  atmname = '1HD1'
         if (atmname .eq. 'HD12')  atmname = '2HD1'
         if (atmname .eq. ' HD2')  atmname = '2HD1'
         if (atmname .eq. 'HD13')  atmname = '3HD1'
         if (atmname .eq. ' HD3')  atmname = '3HD1'
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         if (atmname .eq. ' OG1')  atmname = ' OG '
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = ' HG '
         if (atmname .eq. ' HOG')  atmname = ' HG '
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         if (atmname .eq. ' OG ')  atmname = ' OG1'
         if (atmname .eq. ' CG ')  atmname = ' CG2'
         if (atmname .eq. ' HOG')  atmname = ' HG1'
         if (atmname .eq. 'HG21')  atmname = '1HG2'
         if (atmname .eq. 'HG22')  atmname = '2HG2'
         if (atmname .eq. 'HG23')  atmname = '3HG2'
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         if (atmname .eq. ' SG1')  atmname = ' SG '
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = ' HG '
         if (atmname .eq. ' HSG')  atmname = ' HG '
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HD1')  atmname = '1HD '
         if (atmname .eq. ' HD2')  atmname = '2HD '
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         if (atmname .eq. ' HOH')  atmname = ' HH '
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HNE')  atmname = ' HE1'
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HND')  atmname = ' HD1'
         if (atmname .eq. ' HNE')  atmname = ' HE2'
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HND')  atmname = ' HD1'
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HD ')  atmname = ' HD2'
         if (atmname .eq. ' HE ')  atmname = ' HE1'
         if (atmname .eq. ' HNE')  atmname = ' HE2'
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         if (atmname .eq. ' OD ')  atmname = ' OD1'
         if (atmname .eq. ' ND ')  atmname = ' ND2'
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. 'HD21')  atmname = '1HD2'
         if (atmname .eq. 'HND1')  atmname = '1HD2'
         if (atmname .eq. 'HD22')  atmname = '2HD2'
         if (atmname .eq. 'HND2')  atmname = '2HD2'
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         if (atmname .eq. ' OE ')  atmname = ' OE1'
         if (atmname .eq. ' NE ')  atmname = ' NE2'
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. 'HE21')  atmname = '1HE2'
         if (atmname .eq. 'HNE1')  atmname = '1HE2'
         if (atmname .eq. 'HE22')  atmname = '2HE2'
         if (atmname .eq. 'HNE2')  atmname = '2HE2'
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HE1')  atmname = '1HE '
         if (atmname .eq. ' HE2')  atmname = '2HE '
         if (atmname .eq. ' HE3')  atmname = '3HE '
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HD1')  atmname = '1HD '
         if (atmname .eq. ' HD2')  atmname = '2HD '
         if (atmname .eq. ' HE1')  atmname = '1HE '
         if (atmname .eq. ' HE2')  atmname = '2HE '
         if (atmname .eq. ' HZ1')  atmname = '1HZ '
         if (atmname .eq. 'HNZ1')  atmname = '1HZ '
         if (atmname .eq. ' HZ2')  atmname = '2HZ '
         if (atmname .eq. 'HNZ2')  atmname = '2HZ '
         if (atmname .eq. ' HZ3')  atmname = '3HZ '
         if (atmname .eq. 'HNZ3')  atmname = '3HZ '
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HD1')  atmname = '1HD '
         if (atmname .eq. ' HD2')  atmname = '2HD '
         if (atmname .eq. 'HH11')  atmname = '1HH1'
         if (atmname .eq. 'HN11')  atmname = '1HH1'
         if (atmname .eq. 'HH12')  atmname = '2HH1'
         if (atmname .eq. 'HN12')  atmname = '2HH1'
         if (atmname .eq. 'HH21')  atmname = '1HH2'
         if (atmname .eq. 'HN21')  atmname = '1HH2'
         if (atmname .eq. 'HH22')  atmname = '2HH2'
         if (atmname .eq. 'HN22')  atmname = '2HH2'
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
         if (atmname .eq. ' HD1')  atmname = '1HD '
         if (atmname .eq. ' HD2')  atmname = '2HD '
         if (atmname .eq. ' HE1')  atmname = '1HE '
         if (atmname .eq. 'HNE1')  atmname = '1HE '
         if (atmname .eq. ' HE2')  atmname = '2HE '
         if (atmname .eq. 'HNE2')  atmname = '2HE '
         if (atmname .eq. ' HE3')  atmname = '3HE '
         if (atmname .eq. 'HNE3')  atmname = '3HE '
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         if (atmname .eq. 'HB11')  atmname = '1HB1'
         if (atmname .eq. 'HB12')  atmname = '2HB1'
         if (atmname .eq. 'HB13')  atmname = '3HB1'
         if (atmname .eq. 'HB21')  atmname = '1HB2'
         if (atmname .eq. 'HB22')  atmname = '2HB2'
         if (atmname .eq. 'HB23')  atmname = '3HB2'
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         if (atmname .eq. ' HB1')  atmname = '1HB '
         if (atmname .eq. ' HB2')  atmname = '2HB '
         if (atmname .eq. ' HG1')  atmname = '1HG '
         if (atmname .eq. ' HG2')  atmname = '2HG '
c
c     N-terminal acetyl residue  (ACE)
c
      else if (resname .eq. 'ACE') then
         if (atmname .eq. ' CY ')  atmname = ' C  '
         if (atmname .eq. ' CAY')  atmname = ' CH3'
         if (atmname .eq. ' CA ')  atmname = ' CH3'
         if (atmname .eq. ' OY ')  atmname = ' O  '
         if (atmname .eq. ' HY1')  atmname = '1H  '
         if (atmname .eq. ' HY2')  atmname = '2H  '
         if (atmname .eq. ' HY3')  atmname = '3H  '
         if (atmname .eq. ' H1 ')  atmname = '1H  '
         if (atmname .eq. ' H2 ')  atmname = '2H  '
         if (atmname .eq. ' H3 ')  atmname = '3H  '
c
c     N-terminal formyl residue  (FOR)
c
      else if (resname .eq. 'FOR') then
         if (atmname .eq. ' CY ')  atmname = ' C  '
         if (atmname .eq. ' OY ')  atmname = ' O  '
         if (atmname .eq. ' HY ')  atmname = ' H  '
c
c     C-terminal N-methylamide residue  (NME)
c
      else if (resname .eq. 'NME') then
         if (atmname .eq. ' NT ')  atmname = ' N  '
         if (atmname .eq. ' CT ')  atmname = ' CH3'
         if (atmname .eq. ' CAT')  atmname = ' CH3'
         if (atmname .eq. ' CA ')  atmname = ' CH3'
         if (atmname .eq. ' HNT')  atmname = ' H  '
         if (atmname .eq. ' HT1')  atmname = '1H  '
         if (atmname .eq. ' HT2')  atmname = '2H  '
         if (atmname .eq. ' HT3')  atmname = '3H  '
         if (atmname .eq. ' H1 ')  atmname = '1H  '
         if (atmname .eq. ' H2 ')  atmname = '2H  '
         if (atmname .eq. ' H3 ')  atmname = '3H  '
c
c     C-terminal amide residue  (NH2)
c
      else if (resname .eq. 'NH2') then
         if (atmname .eq. ' NT ')  atmname = ' N  '
         if (atmname .eq. ' HT1')  atmname = '1H  '
         if (atmname .eq. ' HT2')  atmname = '2H  '
         if (atmname .eq. ' H1 ')  atmname = '1H  '
         if (atmname .eq. ' H2 ')  atmname = '2H  '
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
c     ##  subroutine readprm  --  input of force field parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readprm" processes the potential energy parameter file
c     in order to define the default force field parameters
c
c
      subroutine readprm (iprm)
      implicit none
      include 'sizes.i'
      include 'fields.i'
      include 'iounit.i'
      include 'kanang.i'
      include 'kangs.i'
      include 'katoms.i'
      include 'kbonds.i'
      include 'kchrge.i'
      include 'kdipol.i'
      include 'khbond.i'
      include 'kiprop.i'
      include 'kitors.i'
      include 'kmulti.i'
      include 'kopbnd.i'
      include 'korbs.i'
      include 'kpolr.i'
      include 'kstbnd.i'
      include 'ksttor.i'
      include 'ktorsn.i'
      include 'kurybr.i'
      include 'kvdws.i'
      include 'kvdwpr.i'
      integer i,j,iprm,size,next
      integer length,trimtext
      integer ia,ib,ic,id
      integer cls,atn,lig,ft(6)
      integer nvp,nhb,nmp,npi
      integer nb,nb5,nb4,nb3
      integer na,na5,na4,na3
      integer nopb,nu,ndi,nti
      integer nt,nt5,nt4,nbt
      integer nd,nd5,nd4,nd3
      real*8 wght,rd,ep,rdn
      real*8 an1,an2,an3
      real*8 ba1,ba2,ba3
      real*8 aa1,aa2,aa3
      real*8 bt1,bt2,bt3
      real*8 dk,vd,dst,cg,dp,ps,fc
      real*8 bd,el,pol,iz,rp,ss,ts
      real*8 vt(6),st(6),pl(13)
      character*3 pa,pb,pc,pd
      character*8 axt
      character*20 keyword
      character*80 record,string
      logical header
c
      integer nproc,me,master,ibtyp,iptim
      integer kk,iflag,iend,CMSG(80)
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c     initialize the counters for some parameter types
c
      nvp = 0
      nhb = 0
      nb = 0
      nb5 = 0
      nb4 = 0
      nb3 = 0
      na = 0
      na5 = 0
      na4 = 0
      na3 = 0
      nopb = 0
      nu = 0
      ndi = 0
      nti = 0
      nt = 0
      nt5 = 0
      nt4 = 0
      nbt = 0
      nd = 0
      nd5 = 0
      nd4 = 0
      nd3 = 0
      nmp = 0
      npi = 0
casa  -charge2- 
      n_chg2 = 1
c
c     number of characters in an atom number text string
c
      size = 3
c
c     set blank line header before echoed comment lines
c
      header = .true.
c
c     process each line of the parameter file, first
c     extract the keyword at the start of each line
c
      dowhile (.true.)
      iflag=0
c
      IF (MASWRK) THEN
         read (iprm,10,iostat=kk)  record
   10    format (a80)
         if (kk.lt.0) then
            iflag=1
         endif
         IF (GOPARR) THEN
           DO I=1,80
              CMSG(I) = ICHAR(record(I:I))
           enddo
         END IF
      ENDIF
      IF (GOPARR) CALL DDI_BCAST(7777,'I',iflag,1,MASTER)
      if (iflag.eq.1) goto 430
      IF (GOPARR) CALL DDI_BCAST(7777,'I',CMSG,80,MASTER)
      IF (.NOT.MASWRK) THEN
         DO I=1,80
            record(I:I) = CHAR(CMSG(I))
         enddo
      END IF
c
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
c
c     check for a force field modification keyword
c
         call prmkey (record)
c
c     comment line to be echoed to the output
c
         if (keyword(1:5) .eq. 'ECHO ') then
            string = record(next:80)
            length = trimtext (string)
            if (header) then
               header = .false.
               write (iout,20)
   20          format ()
            end if
            if (length .eq. 0) then
               write (iout,30)
   30          format ()
            else
               write (iout,40)  string(1:length)
   40          format (a)
            end if
c
c     atom type definitions and parameters
c
         else if (keyword(1:5) .eq. 'ATOM ') then
            ia = 0
            cls = 0
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,ia,next)
            call getnumb (record,cls,next)
            if (cls .eq. 0)  cls = ia
            atmcls(ia) = cls
            if (ia .ge. maxtyp) then
               write (iout,50)
   50          format (/,' READPRM  --  Too many Atom Types;',
     &                    ' Increase MAXTYP')
               call fatal
            else if (cls .ge. maxclass) then
               write (iout,60)
   60          format (/,' READPRM  --  Too many Atom Classes;',
     &                    ' Increase MAXCLASS')
               call fatal
            end if
            call gettext (record,symbol(ia),next)
            call getstring (record,describe(ia),next)
            string = record(next:80)
            read (string,*,err=70,end=70)  atn,wght,lig
   70       continue
            atmnum(ia) = atn
            weight(ia) = wght
            ligand(ia) = lig
c
c     van der Waals parameters for individual atom types
c
         else if (keyword(1:4) .eq. 'VDW ') then
            ia = 0
            rd = 0.0d0
            ep = 0.0d0
            rdn = 0.0d0
            string = record(next:80)
            read (string,*,err=80,end=80)  ia,rd,ep,rdn
   80       continue
            rad(ia) = rd
            eps(ia) = ep
            reduct(ia) = rdn
c
c     van der Waals parameters for specific atom pairs
c
         else if (keyword(1:6) .eq. 'VDWPR ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:80)
            read (string,*,err=90,end=90)  ia,ib,rd,ep
   90       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nvp = nvp + 1
c
            if(nvp.gt.maxnvp) then
               write(6,*) 'too many VDW parameters'
               call abrt
            end if
c
            if (ia .le. ib) then
               kvpr(nvp) = pa//pb
            else
               kvpr(nvp) = pb//pa
            end if
            radpr(nvp) = rd
            epspr(nvp) = ep
c
c     van der Waals parameters for hydrogen bonding pairs
c
         else if (keyword(1:6) .eq. 'HBOND ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:80)
            read (string,*,err=95,end=95)  ia,ib,rd,ep
   95       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nhb = nhb + 1
c
            if(nhb.gt.maxnhb) then
               write(6,*) 'too many VDW for H-bond parameters'
               call abrt
            end if
c
            if (ia .le. ib) then
               khb(nhb) = pa//pb
            else
               khb(nhb) = pb//pa
            end if
            radhb(nhb) = rd
            epshb(nhb) = ep
c
c     bond stretching parameters
c
         else if (keyword(1:5) .eq. 'BOND ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:80)
            read (string,*,err=100,end=100)  ia,ib,fc,bd
  100       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb = nb + 1
c
            if(nb.gt.maxnb) then
cbb
               write(*,*) nb, maxnb
cbb
               write(6,*) 'too many stretch force constants' , nb, maxnb
               call abrt
            end if
c
            if (ia .le. ib) then
               kb(nb) = pa//pb
            else
               kb(nb) = pb//pa
            end if
            fcon(nb) = fc
            blen(nb) = bd
c
c     bond stretching parameters for 5-membered rings
c
         else if (keyword(1:6) .eq. 'BOND5 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:80)
            read (string,*,err=110,end=110)  ia,ib,fc,bd
  110       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb5 = nb5 + 1
c
            if(nb5.gt.maxnb5) then
               write(6,*) 'too many 5-mem ring stretch force constants'
               call abrt
            end if
c
            if (ia .le. ib) then
               kb5(nb5) = pa//pb
            else
               kb5(nb5) = pb//pa
            end if
            fcon5(nb5) = fc
            blen5(nb5) = bd
c
c     bond stretching parameters for 4-membered rings
c
         else if (keyword(1:6) .eq. 'BOND4 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:80)
            read (string,*,err=120,end=120)  ia,ib,fc,bd
  120       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb4 = nb4 + 1
c
            if(nb4.gt.maxnb4) then
               write(6,*) 'too many 4-mem ring stretch force constants'
               call abrt
            end if
c
            if (ia .le. ib) then
               kb4(nb4) = pa//pb
            else
               kb4(nb4) = pb//pa
            end if
            fcon4(nb4) = fc
            blen4(nb4) = bd
c
c     bond stretching parameters for 3-membered rings
c
         else if (keyword(1:6) .eq. 'BOND3 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:80)
            read (string,*,err=130,end=130)  ia,ib,fc,bd
  130       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb3 = nb3 + 1
c
            if(nb3.gt.maxnb3) then
               write(6,*) 'too many 3-mem ring stretch force constants'
               call abrt
            end if
c
            if (ia .le. ib) then
               kb3(nb3) = pa//pb
            else
               kb3(nb3) = pb//pa
            end if
            fcon3(nb3) = fc
            blen3(nb3) = bd
c
c     bond angle bending parameters
c
         else if (keyword(1:6) .eq. 'ANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            string = record(next:80)
            read (string,*,err=140,end=140)  ia,ib,ic,fc,an1,an2,an3
  140       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na = na + 1
c
            if(na.gt.maxna) then
               write(6,*) 'too many angle force constants, maxna=',maxna
               call abrt
            end if
c
            if (ia .le. ic) then
               ka(na) = pa//pb//pc
            else
               ka(na) = pc//pb//pa
            end if
            con(na) = fc
            ang(1,na) = an1
            ang(2,na) = an2
            ang(3,na) = an3
c
c     angle bending parameters for 5-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE5 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            string = record(next:80)
            read (string,*,err=150,end=150)  ia,ib,ic,fc,an1,an2,an3
  150       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na5 = na5 + 1
c
            if(na5.gt.maxna5) then
               write(6,*) 'too many 5-mem ring angle force constants'
               call abrt
            end if
c
            if (ia .le. ic) then
               ka5(na5) = pa//pb//pc
            else
               ka5(na5) = pc//pb//pa
            end if
            con5(na5) = fc
            ang5(1,na5) = an1
            ang5(2,na5) = an2
            ang5(3,na5) = an3
c
c     angle bending parameters for 4-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE4 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            string = record(next:80)
            read (string,*,err=160,end=160)  ia,ib,ic,fc,an1,an2,an3
  160       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na4 = na4 + 1
c
            if(na4.gt.maxna4) then
               write(6,*) 'too many 4-mem ring angle force constants'
               call abrt
            end if
c
            if (ia .le. ic) then
               ka4(na4) = pa//pb//pc
            else
               ka4(na4) = pc//pb//pa
            end if
            con4(na4) = fc
            ang4(1,na4) = an1
            ang4(2,na4) = an2
            ang4(3,na4) = an3
c
c     angle bending parameters for 3-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE3 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            string = record(next:80)
            read (string,*,err=170,end=170)  ia,ib,ic,fc,an1,an2,an3
  170       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na3 = na3 + 1
c
            if(na3.gt.maxna3) then
               write(6,*) 'too many 3-mem ring angle force constants'
               call abrt
            end if
c
            if (ia .le. ic) then
               ka3(na3) = pa//pb//pc
            else
               ka3(na3) = pc//pb//pa
            end if
            con3(na3) = fc
            ang3(1,na3) = an1
            ang3(2,na3) = an2
            ang3(3,na3) = an3
c
c     stretch-bend parameters
c
         else if (keyword(1:7) .eq. 'STRBND ') then
            ia = 0
            ba1 = 0.0d0
            ba2 = 0.0d0
            ba3 = 0.0d0
            string = record(next:80)
            read (string,*,err=180,end=180)  ia,ba1,ba2,ba3
  180       continue
            stbn(1,ia) = ba1
            stbn(2,ia) = ba2
            stbn(3,ia) = ba3
c
c     out-of-plane bend parameters
c
         else if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            string = record(next:80)
            read (string,*,err=190,end=190)  ia,ib,fc
  190       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nopb = nopb + 1
c
            if(nopb.gt.maxnopb) then
               write(6,*) 'too many OOP bend force constants'
               call abrt
            end if
c
            kaopb(nopb) = pa//pb
            copb(nopb) = fc
c
c     Urey-Bradley parameters
c
         else if (keyword(1:9) .eq. 'UREYBRAD ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            dst = 0.0d0
            string = record(next:80)
            read (string,*,err=200,end=200)  ia,ib,ic,fc,dst
  200       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            nu = nu + 1
c
            if(nu.gt.maxnu) then
               write(6,*) 'too many UREYBRAD force constants'
               call abrt
            end if
c
            if (ia .le. ic) then
               ku(nu) = pa//pb//pc
            else
               ku(nu) = pc//pb//pa
            end if
            ucon(nu) = fc
            dst13(nu) = dst
c
c     angle-angle parameters
c
         else if (keyword(1:7) .eq. 'ANGANG ') then
            ia = 0
            aa1 = 0.0d0
            aa2 = 0.0d0
            aa3 = 0.0d0
            string = record(next:80)
            read (string,*,err=210,end=210)  ia,aa1,aa2,aa3
  210       continue
            anan(1,ia) = aa1
            anan(2,ia) = aa2
            anan(3,ia) = aa3
c
c     improper dihedral parameters
c
         else if (keyword(1:9) .eq. 'IMPROPER ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            dk = 0.0d0
            vd = 0.0d0
            string = record(next:80)
            read (string,*,err=220,end=220)  ia,ib,ic,id,dk,vd
  220       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            ndi = ndi + 1
c
            if(ndi.gt.maxndi) then
               write(6,*) 'too many improper force constants'
               call abrt
            end if
c
            kdi(ndi) = pa//pb//pc//pd
            dcon(ndi) = dk
            tdi(ndi) = vd
c
c     improper torsional parameters
c
         else if (keyword(1:8) .eq. 'IMPTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            string = record(next:80)
            read (string,*,err=230,end=230)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  230       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nti = nti + 1
c
            if(nti.gt.maxnti) then
               write(6,*) 'too many improper-torsion force constants'
               call abrt
            end if
c
            kti(nti) = pa//pb//pc//pd
            call torphase (ft,vt,st)
            ti1(1,nti) = vt(1)
            ti1(2,nti) = st(1)
            ti2(1,nti) = vt(2)
            ti2(2,nti) = st(2)
            ti3(1,nti) = vt(3)
            ti3(2,nti) = st(3)
c
c     torsional parameters
c
         else if (keyword(1:8) .eq. 'TORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            string = record(next:80)
            read (string,*,err=240,end=240)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  240       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt = nt + 1
c
            if(nt.gt.maxnt) then
               write(6,*) 'too many torsion force constants'
               call abrt
            end if
c
            if (ib .lt. ic) then
               kt(nt) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt(nt) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt(nt) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt(nt) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t1(1,nt) = vt(1)
            t1(2,nt) = st(1)
            t2(1,nt) = vt(2)
            t2(2,nt) = st(2)
            t3(1,nt) = vt(3)
            t3(2,nt) = st(3)
            t4(1,nt) = vt(4)
            t4(2,nt) = st(4)
            t5(1,nt) = vt(5)
            t5(2,nt) = st(5)
            t6(1,nt) = vt(6)
            t6(2,nt) = st(6)
c
c     torsional parameters for 5-membered rings
c
         else if (keyword(1:9) .eq. 'TORSION5 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            string = record(next:80)
            read (string,*,err=250,end=250)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  250       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt5 = nt5 + 1

            if(nt5.gt.maxnt5) then
               write(6,*) 'too many 5-mem ring torsion force constants'
               call abrt
            end if
c
            if (ib .lt. ic) then
               kt5(nt5) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt5(nt5) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt5(nt5) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt5(nt5) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t15(1,nt5) = vt(1)
            t15(2,nt5) = st(1)
            t25(1,nt5) = vt(2)
            t25(2,nt5) = st(2)
            t35(1,nt5) = vt(3)
            t35(2,nt5) = st(3)
            t45(1,nt5) = vt(4)
            t45(2,nt5) = st(4)
            t55(1,nt5) = vt(5)
            t55(2,nt5) = st(5)
            t65(1,nt5) = vt(6)
            t65(2,nt5) = st(6)
c
c     torsional parameters for 4-membered rings
c
         else if (keyword(1:9) .eq. 'TORSION4 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            string = record(next:80)
            read (string,*,err=260,end=260)  ia,ib,ic,id,
     &                                       (vt(i),st(i),ft(i),i=1,6)
  260       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt4 = nt4 + 1
c
            if(nt4.gt.maxnt4) then
               write(6,*) 'too many 4-mem ring torsion force constants'
               call abrt
            end if
c
            if (ib .lt. ic) then
               kt4(nt4) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt4(nt4) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt4(nt4) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt4(nt4) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t14(1,nt4) = vt(1)
            t14(2,nt4) = st(1)
            t24(1,nt4) = vt(2)
            t24(2,nt4) = st(2)
            t34(1,nt4) = vt(3)
            t34(2,nt4) = st(3)
            t44(1,nt4) = vt(4)
            t44(2,nt4) = st(4)
            t54(1,nt4) = vt(5)
            t54(2,nt4) = st(5)
            t64(1,nt4) = vt(6)
            t64(2,nt4) = st(6)
c
c     stretch-torsion parameters
c
         else if (keyword(1:8) .eq. 'STRTORS ') then
            ia = 0
            ib = 0
            bt1 = 0.0d0
            bt2 = 0.0d0
            bt3 = 0.0d0
            string = record(next:80)
            read (string,*,err=270,end=270)  ia,ib,bt1,bt2,bt3
  270       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nbt = nbt + 1
c
            if(nbt.gt.maxnbt) then
               write(6,*) 'too many stretch torsion force constants'
               call abrt
            end if
c
            if (ia .le. ib) then
               kbt(nbt) = pa//pb
            else
               kbt(nbt) = pb//pa
            end if
            btcon(1,nbt) = bt1
            btcon(2,nbt) = bt2
            btcon(3,nbt) = bt3
c
c     atomic partial charge parameters
c
         else if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            string = record(next:80)
            read (string,*,err=280,end=280)  ia,cg
  280       continue
            chg(ia) = cg
casa  -charge2- begin
c     atomic partial charge parameters2 by asada
c
         else if (keyword(1:7) .eq. 'CHARGE2 ') then
            ia = 0
            cg = 0.0d0
            call gettext (record,atname(n_chg2),next)
            call upcase (atname(n_chg2))
            string = record(next:120)
            read (string,*,err=351,end=351)  cg
  351       continue
            if (cg .ne. 0)  chg2(n_chg2) = cg
            n_chg2 = n_chg2 + 1
casa  -charge2- end
c
c     bond dipole moment parameters
c
         else if (keyword(1:7) .eq. 'DIPOLE ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:80)
            read (string,*,err=290,end=290)  ia,ib,dp,ps
  290       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nd = nd + 1
c
            if(nd.gt.maxnd) then
               write(6,*) 'too many dipole parameters'
               call abrt
            end if
c
            if (ia .le. ib) then
               kd(nd) = pa//pb
            else
               kd(nd) = pb//pa
            end if
            dpl(nd) = dp
            pos(nd) = ps
c
c     bond dipole moment parameters for 5-membered rings
c
         else if (keyword(1:8) .eq. 'DIPOLE5 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:80)
            read (string,*,err=300,end=300)  ia,ib,dp,ps
  300       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nd5 = nd5 + 1
c
            if(nd5.gt.maxnd5) then
               write(6,*) 'too many 5-mem ring dipole parameters'
               call abrt
            end if
c
            if (ia .le. ib) then
               kd5(nd5) = pa//pb
            else
               kd5(nd5) = pb//pa
            end if
            dpl5(nd5) = dp
            pos5(nd5) = ps
c
c     bond dipole moment parameters for 4-membered rings
c
         else if (keyword(1:8) .eq. 'DIPOLE4 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:80)
            read (string,*,err=310,end=310)  ia,ib,dp,ps
  310       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nd4 = nd4 + 1
c
            if(nd4.gt.maxnd4) then
               write(6,*) 'too many 4-mem ring dipole parameters'
               call abrt
            end if
c
            if (ia .le. ib) then
               kd4(nd4) = pa//pb
            else
               kd4(nd4) = pb//pa
            end if
            dpl4(nd4) = dp
            pos4(nd4) = ps
c
c     bond dipole moment parameters for 3-membered rings
c
         else if (keyword(1:8) .eq. 'DIPOLE3 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:80)
            read (string,*,err=320,end=320)  ia,ib,dp,ps
  320       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nd3 = nd3 + 1
c
            if(nd3.gt.maxnd3) then
               write(6,*) 'too many 3-mem ring dipole parameters'
               call abrt
            end if
c
            if (ia .le. ib) then
               kd3(nd3) = pa//pb
            else
               kd3(nd3) = pb//pa
            end if
            dpl3(nd3) = dp
            pos3(nd3) = ps
c
c     atomic multipole moment parameters
c
         else if (keyword(1:10) .eq. 'MULTIPOLE ') then
            ia = 0
            ib = 0
            ic = 0
            axt = 'Z-then-X'
            do i = 1, 13
               pl(i) = 0.0d0
            end do
            string = record(next:80)
            read (string,*,err=370,end=370)  ia,ib,ic,pl(1)
c
      IF (MASWRK) THEN
            read (iprm,330,iostat=kk)  record
  330       format (a80)                                                        
         if (kk.lt.0) then
            iflag=1
         endif
         IF (GOPARR) THEN
           DO I=1,80
              CMSG(I) = ICHAR(record(I:I))
           enddo
         END IF
      ENDIF
      IF (GOPARR) CALL DDI_BCAST(7777,'I',iflag,1,MASTER)
      if (iflag.eq.1) goto 430
      IF (GOPARR) CALL DDI_BCAST(7777,'I',CMSG,80,MASTER)
      IF (.NOT.MASWRK) THEN
         DO I=1,80
            record(I:I) = CHAR(CMSG(I))
         enddo
      END IF
            read (record,*,err=370,end=370)  pl(2),pl(3),pl(4)
       write(*,*) 'PL ME',pl(2),pl(3),pl(4),me
c
      IF (MASWRK) THEN
            read (iprm,340,iostat=kk)  record
  340       format (a80)                                                        
         if (kk.lt.0) then
            iflag=1
         endif
         IF (GOPARR) THEN
           DO I=1,80
              CMSG(I) = ICHAR(record(I:I))
           enddo
         END IF
      ENDIF
      IF (GOPARR) CALL DDI_BCAST(7777,'I',iflag,1,MASTER)
      if (iflag.eq.1) goto 430
      IF (GOPARR) CALL DDI_BCAST(7777,'I',CMSG,80,MASTER)
      IF (.NOT.MASWRK) THEN
         DO I=1,80
            record(I:I) = CHAR(CMSG(I))
         enddo
      END IF
            read (record,*,err=370,end=370)  pl(5)
c
      IF (MASWRK) THEN
            read (iprm,350,iostat=kk)  record
  350       format (a80)                                                        
         if (kk.lt.0) then
            iflag=1
         endif
         IF (GOPARR) THEN
           DO I=1,80
              CMSG(I) = ICHAR(record(I:I))
           enddo
         END IF
      ENDIF
      IF (GOPARR) CALL DDI_BCAST(7777,'I',iflag,1,MASTER)
      if (iflag.eq.1) goto 430
      IF (GOPARR) CALL DDI_BCAST(7777,'I',CMSG,80,MASTER)
      IF (.NOT.MASWRK) THEN
         DO I=1,80
            record(I:I) = CHAR(CMSG(I))
         enddo
      END IF
            read (record,*,err=370,end=370)  pl(8),pl(9)
c
      IF (MASWRK) THEN
            read (iprm,360,iostat=kk)  record
  360       format (a80)                                                        
         if (kk.lt.0) then
            iflag=1
         endif
         IF (GOPARR) THEN
           DO I=1,80
              CMSG(I) = ICHAR(record(I:I))
           enddo
         END IF
      ENDIF
      IF (GOPARR) CALL DDI_BCAST(7777,'I',iflag,1,MASTER)
      if (iflag.eq.1) goto 430
      IF (GOPARR) CALL DDI_BCAST(7777,'I',CMSG,80,MASTER)
      IF (.NOT.MASWRK) THEN
         DO I=1,80
            record(I:I) = CHAR(CMSG(I))
         enddo
      END IF
            read (record,*,err=370,end=370)  pl(11),pl(12),pl(13)
  370       continue
            if (ic .lt. 0) then
               ic = -ic
               axt = 'Bisector'
            end if
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            nmp = nmp + 1
c
            if(nmp.gt.maxnmp) then
               write(6,*) 'too many multipole parameters'
               call abrt
            end if
c
            kmp(nmp) = pa//pb//pc
            mpaxis(nmp) = axt
            multip(1,nmp) = pl(1)
            multip(2,nmp) = pl(2)
            multip(3,nmp) = pl(3)
            multip(4,nmp) = pl(4)
            multip(5,nmp) = pl(5)
            multip(6,nmp) = pl(8)
            multip(7,nmp) = pl(11)
            multip(8,nmp) = pl(8)
            multip(9,nmp) = pl(9)
            multip(10,nmp) = pl(12)
            multip(11,nmp) = pl(11)
            multip(12,nmp) = pl(12)
            multip(13,nmp) = pl(13)
c
c     atomic dipole polarizability parameters
c
         else if (keyword(1:9) .eq. 'POLARIZE ') then
            ia = 0
            pol = 0.0d0
            string = record(next:80)
            read (string,*,err=380,end=380)  ia,pol
  380       continue
            polr(ia) = pol
c
c     conjugated pisystem atom parameters
c
         else if (keyword(1:7) .eq. 'PIATOM ') then
            ia = 0
            el = 0.0d0
            iz = 0.0d0
            rp = 0.0d0
            string = record(next:80)
            read (string,*,err=390,end=390)  ia,el,iz,rp
  390       continue
            electron(ia) = el
            ionize(ia) = iz
            repulse(ia) = rp
c
c     conjugated pisystem bond parameters
c
         else if (keyword(1:7) .eq. 'PIBOND ') then
            ia = 0
            ib = 0
            ss = 0.0d0
            ts = 0.0d0
            string = record(next:80)
            read (string,*,err=400,end=400)  ia,ib,ss,ts
  400       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            npi = npi + 1
c
            if(npi.gt.maxnpi) then
               write(6,*) 'too many pibond parameters'
               call abrt
            end if
c
            if (ia .le. ib) then
               kpi(npi) = pa//pb
            else
               kpi(npi) = pb//pa
            end if
            sslope(npi) = ss
            tslope(npi) = ts
c
c     biopolymer atom type conversion definitions
c
         else if (keyword(1:8) .eq. 'BIOTYPE ') then
            string = record(next:80)
            read (string,*,err=410,end=410)  ia
            call getword (record,string,next)
            call getstring (record,string,next)
            string = record(next:80)
            read (string,*,err=410,end=410)  ib
            if (ia .ge. maxbio) then
               write (iout,440)
  440          format (/,' READPRM  --  Too many Atom biotypes;',
     &                    ' Increase MAXBIO')
               call fatal
            end if
            biotyp(ia) = ib
  410       continue
         end if
  420    continue
      end do
  430 continue
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
c     ##  subroutine readseq  --  read biopolymer sequence file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "readseq" gets a biopolymer sequence containing one or more
c     separate chains from an external file; all lines containing
c     sequence must begin with the starting sequence number, the
c     actual sequence is read from subsequent nonblank characters
c
c
      subroutine readseq (iseq)
      implicit none
      include 'sizes.i'
c     include 'aminos.i'
      include 'files.i'
      include 'iounit.i'
      include 'sequen.i'
      integer i,j,iseq,number
      integer next,length,trimtext
      character*1 letter
      character*60 seqfile
      character*80 record
      logical exist,opened
      include 'aminos.i'
c
c
c     open the input file if it has not already been done
c
      inquire (unit=iseq,opened=opened)
      if (.not. opened) then
         seqfile = filename(1:leng)//'.seq'
         call version (seqfile,'old')
         inquire (file=seqfile,exist=exist)
         if (exist) then
            open (unit=iseq,file=seqfile,status='old')
            rewind (unit=iseq)
         else
            write (iout,10)
   10       format (/,' READSEQ  --  Unable to Open Biopolymer',
     &                 ' Sequence File')
            call fatal
         end if
      end if
c
c     zero out the number and type of residues
c
      nseq = 0
      nchain = 0
      do i = 1, maxres
         seq(i) = ' '
      end do
c
c     read in the biopolymer sequence file
c
      dowhile (.true.)
         read (iseq,20,err=30,end=30)  record
   20    format (a80)
         length = trimtext (record)
         next = 1
         call gettext (record,letter,next)
         if (letter.ge.'0' .and. letter.le.'9') then
            next = 1
            letter = ' '
         end if
         call getnumb (record,number,next)
         if (number .eq. 1) then
            nchain = nchain + 1
            ichain(1,nchain) = nseq + 1
            chnnam(nchain) = letter
         end if
         do i = next, length
            letter = record(i:i)
            if (letter .ne. ' ') then
               nseq = nseq + 1
               seq(nseq) = letter
            end if
         end do
      end do
   30 continue
c
c     set the last residue in each sequence chain
c
      do i = 1, nchain-1
         ichain(2,i) = ichain(1,i+1) - 1
      end do
      if (nchain .ne. 0)  ichain(2,nchain) = nseq
c
c     find the residue type for each sequence element
c
      do i = 1, nseq
         seqtyp(i) = 0
         do j = 1, maxamino
            if (seq(i) .eq. amino1(j))  seqtyp(i) = j
         end do
      end do
      if (.not. opened)  close (unit=iseq)
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
c     ##  subroutine readxyz  --  input of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readxyz" gets a set of Cartesian coordinates
c     from an external disk file
c
c
      subroutine readxyz (ixyz,iw)
      implicit double precision (a-h,o-z)
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'iounit.i'
      include 'titles.i'
      integer iw,nproc,me,master,ibtyp,iptim,ieof,ierr,i,j,k,m,ixyz
      integer jeof,next,first,last
      integer nexttext,trimtext,ifind
      integer list(maxatm)
      integer CMSG(120)
      character*60 xyzfile
      character*120 record,string
      logical exist,opened,quit,reorder
      LOGICAL GOPARR,DSKWRK,MASWRK
      double precision rfind
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     read the title line and get the number of atoms
c
      quit = .true.
      i = 0
      IF (MASWRK) THEN
         read (ixyz,20)  record
   20    format (a120)                                                                         
         IF (GOPARR) THEN
           CALL DDI_BCAST(7777,'I',0,1,MASTER)
           DO I=1,120
              CMSG(I) = ICHAR(record(I:I))
           enddo
         END IF
      ELSE
         IF (GOPARR)
     *     CALL DDI_BCAST(7777,'I',IEND,1,MASTER)
      ENDIF
      IF (GOPARR) CALL DDI_BCAST(7777,'I',CMSG,120,MASTER)
      IF (.NOT.MASWRK) THEN
         DO I=1,120
            record(I:I) = CHAR(CMSG(I))
         enddo
      END IF
      next = 1
      call gettext (record,string,next)
      read (string,*,err=60,end=60)  n
c
c     extract the title and determine its length
c
      string = record(next:120)
      first = nexttext (string)
      last = trimtext (string)
      if (last .eq. 0) then
         ttitle = ' '
         ltitle = 0
      else
         ttitle = string(first:last)
         ltitle = trimtext (ttitle)
      end if
      CALL SEQREW(Ixyz)
      CALL FNDGRP(Ixyz,' $TINXYZ',JEOF)
c
c  -- read tinker xyz file and echo to the output file
c
      call opncrd(ixyz,-iw)
      ieof = 0
      call rdcard('$TINXYZ1',ieof)
      if (ieof.eq.1) call abrt
C
      n=ifind('NATOM ',ierr)
      if (ierr.ne.0) call abrt
      call gstrng(ttitle,-70)
C
c     check for too many total atoms in the file
c
      if (n .gt. maxatm) then
         if (maswrk) write (iout,30)  maxatm
   30    format (/,' READXYZ  --  The Maximum of',i6,' Atoms',
     &              ' has been Exceeded')
         call abrt
      end if
c
c     read the coordinates and connectivities for each atom
c
      do i = 1, n
         next = 1
         do j = 1, maxval
            i12(j,i) = 0
      end do
      call rdcard('$TINXYZ2',ieof)
      if (ieof.eq.1) go to 55
      tag(i)=ifind('ANUM  ',ierr)
      if (ierr.ne.0) go to 55
      call gstrng(name(i),-10)
casa  -printxyz-
      xyzname(i) = name(i)
      x(i)=rfind('COORD ',ierr)
      if (ierr.ne.0) go to 55
      y(i)=rfind('COORD ',ierr)
      if (ierr.ne.0) go to 55
      z(i)=rfind('COORD ',ierr)
      if (ierr.ne.0) go to 55
      type(i)=ifind('TYPE  ',ierr)
      if (ierr.ne.0) go to 55
      do 50 j=1,maxval
         i12(j,i)=ifind('CONNET',ierr)
         if (ierr.ne.0) go to 55
 50   continue
c
      end do
 55   quit = .false.
   60 continue
cjrs
c      if (.not. opened)  close (unit=ixyz)
cjrs
c
c     an error occurred in reading the coordinate file
c
      if (quit) then
         if (maswrk) write (iout,70)  i
   70    format (/,' READXYZ  --  Error in Coordinate File at Atom',i6)
         call fatal
      end if
c
c     for each atom, count and sort its attached atoms
c
      do i = 1, n
         n12(i) = 0
         do j = maxval, 1, -1
            if (i12(j,i) .ne. 0) then
               n12(i) = j
               goto 80
            end if
         end do
   80    continue
         call sort (n12(i),i12(1,i))
      end do
c
c     check for scrambled atom order and attempt to renumber
c
      reorder = .false.
      do i = 1, n
         list(tag(i)) = i
         if (tag(i) .ne. i)  reorder = .true.
      end do
      if (reorder) then
         if (maswrk) write (iout,90)
   90    format (/,' READXYZ  --  Atom Labels not Sequential,',
     &              ' Attempting to Renumber')
         do i = 1, n
            tag(i) = i
            do j = 1, n12(i)
               i12(j,i) = list(i12(j,i))
            end do
            call sort (n12(i),i12(1,i))
         end do
      end if
c
c     make sure that all connectivities are bidirectional
c
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            do m = 1, n12(k)
               if (i12(m,k) .eq. i)  goto 110
            end do
            if (maswrk) write (iout,100)  k,i
  100       format (/,' READXYZ  --  Check Connection of Atom',
     &                 i6,' to Atom',i6)
            call fatal
  110       continue
         end do
      end do
C
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
c     ##  subroutine restrain  --  initialize geometric restraints  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "restrain" defines any geometric restraint interaction terms
c     to be included in the potential energy calculation
c
c
      subroutine restrain
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'keys.i'
      include 'potent.i'
      include 'restrn.i'
      integer i,ia,ib,ic,id,ip,next
      real*8 xp,yp,zp,p1
      real*8 d1,d2,d3
      real*8 t1,t2,t3
      real*8 bndleng,dihedral
      character*20 keyword
      character*80 record,string
      logical exist
c
c
c     set the default values for the restraint variables
c
      npfix = 0
      ndfix = 0
      ntfix = 0
      depth = 0.0d0
      width = 0.0d0
      use_basin = .false.
      use_wall = .false.
c
c     search the keywords for restraint parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
c
c     specified atom will be restrained to a specified position
c
         if (keyword(1:18) .eq. 'RESTRAIN-POSITION ') then
            xp = 0.0d0
            yp = 0.0d0
            zp = 0.0d0
            p1 = 0.0d0
            string = record(next:80)
            read (string,*,err=10,end=10)  ip,xp,yp,zp,p1
   10       continue
            if (p1 .eq. 0.0d0)  p1 = 100.0d0
            npfix = npfix + 1
            ipfix(npfix) = ip
            xpfix(npfix) = xp
            ypfix(npfix) = yp
            zpfix(npfix) = zp
            pfix(npfix) = p1
            if (npfix .gt. maxfix) then
               write (iout,20)
   20          format (/,' RESTRAIN  --  Too many Position',
     &                    ' Restraints; Increase MAXFIX')
               call fatal
            end if
c
c     distance between two atoms will be restrained to a range
c
         else if (keyword(1:18) .eq. 'RESTRAIN-DISTANCE ') then
            d1 = 0.0d0
            d2 = 0.0d0
            d3 = 0.0d0
            string = record(next:80)
            exist = .false.
            read (string,*,err=30,end=30)  ia,ib,d1
            exist = .true.
   30       continue
            read (string,*,err=40,end=40)  ia,ib,d1,d2,d3
   40       continue
            if (.not. exist)  d1 = bndleng (ia,ib)
            if (d2 .eq. 0.0d0)  d2 = d1
            if (d3 .eq. 0.0d0)  d3 = 100.0d0
            ndfix = ndfix + 1
            idfix(1,ndfix) = ia
            idfix(2,ndfix) = ib
            dfix(1,ndfix) = d1
            dfix(2,ndfix) = d2
            dfix(3,ndfix) = d3
            if (ndfix .gt. maxfix) then
               write (iout,50)
   50          format (/,' RESTRAIN  --  Too many Distance',
     &                    ' Restraints; Increase MAXFIX')
               call fatal
            end if
c
c     dihedral angle of four atoms will be restrained to a range
c
         else if (keyword(1:18).eq.'RESTRAIN-DIHEDRAL ') then
            t1 = 0.0d0
            t2 = 0.0d0
            t3 = 0.0d0
            string = record(next:80)
            exist = .false.
            read (string,*,err=60,end=60)  ia,ib,ic,id,t1
            exist = .true.
   60       continue
            read (string,*,err=70,end=70)  ia,ib,ic,id,t1,t2,t3
            exist = .true.
   70       continue
            if (.not. exist)  t1 = dihedral (ia,ib,ic,id)
            if (t2 .eq. 0.0d0)  t2 = t1
            if (t3 .eq. 0.0d0)  t3 = 1.0d0
            dowhile (t1 .gt. 180.0d0)
               t1 = t1 - 360.0d0
            end do
            dowhile (t1 .lt. -180.0d0)
               t1 = t1 + 360.0d0
            end do
            dowhile (t2 .gt. 180.0d0)
               t2 = t2 - 360.0d0
            end do
            dowhile (t2 .lt. -180.0d0)
               t2 = t2 + 360.0d0
            end do
            ntfix = ntfix + 1
            itfix(1,ntfix) = ia
            itfix(2,ntfix) = ib
            itfix(3,ntfix) = ic
            itfix(4,ntfix) = id
            tfix(1,ntfix) = t1
            tfix(2,ntfix) = t2
            tfix(3,ntfix) = t3
            if (ntfix .gt. maxfix) then
               write (iout,80)
   80          format (/,' RESTRAIN  --  Too many Torsional',
     &                    ' Restraints; Increase MAXFIX')
               call fatal
            end if
c
c     shallow Gaussian basin restraint between all pairs of atoms
c
         else if (keyword(1:6) .eq. 'BASIN ') then
            string = record(next:80)
            read (string,*,err=90,end=90)  depth,width
   90       continue
            use_basin = .true.
            if (depth .eq. 0.0d0)  use_basin = .false.
            if (width .eq. 0.0d0)  use_basin = .false.
            if (depth .gt. 0.0d0)  depth = -depth
c
c     boundary wall to restrain atoms to a spherical droplet
c
         else if (keyword(1:5) .eq. 'WALL ') then
            string = record(next:80)
            read (string,*,err=100,end=100)  rwall
  100       continue
            if (rwall .gt. 0.0d0)  use_wall = .true.
         end if
      end do
c
c     if geometric restraints are not used, turn off the potential
c
      if (npfix.eq.0 .and. ndfix.eq.0 .and. ntfix.eq.0 .and.
     &         .not.use_basin .and. .not.use_wall) then
         use_geom = .false.
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
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine rings  --  locate and store small rings  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "rings" searches the structure for small rings and stores
c     their component atoms; code to remove the reducible rings
c     consisting of smaller rings is commented in this version
c     since reducible rings are needed for parameter assignment
c
c
      subroutine rings
      implicit none
      include 'sizes.i'
      include 'angle.i'
c     include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'ring.i'
      include 'tors.i'
      integer i,j,k,ii,iii,imax
      integer ia,ib,ic,id,ie,ig
c     integer m,list(maxatm)
c     integer list1,list2,list3,list4
c
c
c     zero out the number of small rings in the structure
c
      nring3 = 0
      nring4 = 0
      nring5 = 0
      nring6 = 0
c
c     parse the structure for bonds, angles and torsions
c
      if (nbond .eq. 0)  call bonds
      if (nangle .eq. 0)  call angles
      if (ntors .eq. 0)  call torsions
c
c     search for and store all of the 3-membered rings
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (ib.lt.ia .and. ib.lt.ic) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ic) then
                  nring3 = nring3 + 1
                  if (nring3 .gt. maxring) then
                     write (iout,10)
   10                format (/,' RINGS  --  Too many 3-Membered Rings')
                     call fatal
                  end if
                  iring3(1,nring3) = ia
                  iring3(2,nring3) = ib
                  iring3(3,nring3) = ic
                  goto 20
               end if
            end do
   20       continue
         end if
      end do
c
c     search for and store all of the 4-membered rings
c
c     do i = 1, n
c        list(i) = 0
c     end do
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         imax = max(ia,ib,ic)
         do j = 1, n12(ia)
            id = i12(j,ia)
            if (id .gt. imax) then
               do k = 1, n12(ic)
                  if (i12(k,ic) .eq. id) then
                     nring4 = nring4 + 1
                     if (nring4 .gt. maxring) then
                        write (iout,30)
   30                   format (/,' RINGS  --  Too many',
     &                             ' 4-Membered Rings')
                        call fatal
                     end if
                     iring4(1,nring4) = ia
                     iring4(2,nring4) = ib
                     iring4(3,nring4) = ic
                     iring4(4,nring4) = id
c
c     remove the ring if it is reducible into smaller rings
c
c                    list(ia) = nring4
c                    list(ib) = nring4
c                    list(ic) = nring4
c                    list(id) = nring4
c                    do m = 1, nring3
c                       list1 = list(iring3(1,m))
c                       list2 = list(iring3(2,m))
c                       list3 = list(iring3(3,m))
c                       if (list1.eq.nring4 .and. list2.eq.nring4
c    &                          .and. list3.eq.nring4) then
c                          nring4 = nring4 - 1
c                          list(ia) = 0
c                          list(ib) = 0
c                          list(ic) = 0
c                          list(id) = 0
c                          goto 40
c                       end if
c                    end do
                     goto 40
                  end if
               end do
   40          continue
            end if
         end do
      end do
c
c     search for and store all of the 5-membered rings
c
c     do i = 1, n
c        list(i) = 0
c     end do
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         imax = max(ia,ib,ic,id)
         do j = 1, n12(ia)
            ie = i12(j,ia)
            if (ie .gt. imax) then
               do k = 1, n12(id)
                  if (i12(k,id) .eq. ie) then
                     nring5 = nring5 + 1
                     if (nring5 .gt. maxring) then
                        write (iout,50)
   50                   format (/,' RINGS  --  Too many',
     &                             ' 5-Membered Rings')
                        call fatal
                     end if
                     iring5(1,nring5) = ia
                     iring5(2,nring5) = ib
                     iring5(3,nring5) = ic
                     iring5(4,nring5) = id
                     iring5(5,nring5) = ie
c
c     remove the ring if it is reducible into smaller rings
c
c                    list(ia) = nring5
c                    list(ib) = nring5
c                    list(ic) = nring5
c                    list(id) = nring5
c                    list(ie) = nring5
c                    do m = 1, nring3
c                       list1 = list(iring3(1,m))
c                       list2 = list(iring3(2,m))
c                       list3 = list(iring3(3,m))
c                       if (list1.eq.nring5 .and. list2.eq.nring5
c    &                          .and. list3.eq.nring5) then
c                          nring5 = nring5 - 1
c                          list(ia) = 0
c                          list(ib) = 0
c                          list(ic) = 0
c                          list(id) = 0
c                          list(ie) = 0
c                          goto 60
c                       end if
c                    end do
                     goto 60
                  end if
               end do
   60          continue
            end if
         end do
      end do
c
c     search for and store all of the 6-membered rings
c
c     do i = 1, n
c        list(i) = 0
c     end do
      do i = 1, nangle
         ib = iang(1,i)
         ic = iang(2,i)
         id = iang(3,i)
         do ii = 1, n12(ib)
            ia = i12(ii,ib)
            if (ia.ne.ic .and. ia.ne.id) then
               do iii = 1, n12(id)
                  ie = i12(iii,id)
                  if (ie.ne.ic .and. ie.ne.ib .and. ie.ne.ia) then
                     imax = max(ia,ib,ic,id,ie)
                     do j = 1, n12(ia)
                        ig = i12(j,ia)
                        if (ig .gt. imax) then
                           do k = 1, n12(ie)
                              if (i12(k,ie) .eq. ig) then
                                 nring6 = nring6 + 1
                                 if (nring6 .gt. maxring) then
                                    write (iout,70)
   70                               format (/,' RINGS  --  Too many',
     &                                         ' 6-Membered Rings')
                                    call fatal
                                 end if
                                 iring6(1,nring6) = ia
                                 iring6(2,nring6) = ib
                                 iring6(3,nring6) = ic
                                 iring6(4,nring6) = id
                                 iring6(5,nring6) = ie
                                 iring6(6,nring6) = ig
c
c     remove the ring if it is reducible into smaller rings
c
c                                list(ia) = nring6
c                                list(ib) = nring6
c                                list(ic) = nring6
c                                list(id) = nring6
c                                list(ie) = nring6
c                                list(ig) = nring6
c                                do m = 1, nring3
c                                   list1 = list(iring3(1,m))
c                                   list2 = list(iring3(2,m))
c                                   list3 = list(iring3(3,m))
c                                   if (list1.eq.nring6 .and.
c    &                                  list2.eq.nring6 .and.
c    &                                  list3.eq.nring6) then
c                                      nring6 = nring6 - 1
c                                      list(ia) = 0
c                                      list(ib) = 0
c                                      list(ic) = 0
c                                      list(id) = 0
c                                      list(ie) = 0
c                                      list(ig) = 0
c                                      goto 80
c                                   end if
c                                end do
c                                do m = 1, nring4
c                                   list1 = list(iring4(1,m))
c                                   list2 = list(iring4(2,m))
c                                   list3 = list(iring4(3,m))
c                                   list4 = list(iring4(4,m))
c                                   if (list1.eq.nring6 .and.
c    &                                  list2.eq.nring6 .and.
c    &                                  list3.eq.nring6 .and.
c    &                                  list4.eq.nring6) then
c                                      nring6 = nring6 - 1
c                                      list(ia) = 0
c                                      list(ib) = 0
c                                      list(ic) = 0
c                                      list(id) = 0
c                                      list(ie) = 0
c                                      list(ig) = 0
c                                      goto 80
c                                   end if
c                                end do
   80                            continue
                              end if
                           end do
                        end if
                     end do
                  end if
               end do
            end if
         end do
      end do
c
c     print out lists of the small rings in the structure
c
      if (debug) then
         if (nring3 .gt. 0) then
            write (iout,90)
   90       format (/,' Three-Membered Rings Contained',
     &                 ' in the Structure :',
     &              //,11x,'Ring',14x,'Atoms in Ring',/)
            do i = 1, nring3
               write (iout,100)  i,(iring3(j,i),j=1,3)
  100          format (9x,i5,10x,3i6)
            end do
         end if
         if (nring4 .gt. 0) then
            write (iout,110)
  110       format (/,' Four-Membered Rings Contained',
     &                 ' in the Structure :',
     &              //,11x,'Ring',17x,'Atoms in Ring',/)
            do i = 1, nring4
               write (iout,120)  i,(iring4(j,i),j=1,4)
  120          format (9x,i5,10x,4i6)
            end do
         end if
         if (nring5 .gt. 0) then
            write (iout,130)
  130       format (/,' Five-Membered Rings Contained',
     &                 ' in the Structure :',
     &              //,11x,'Ring',20x,'Atoms in Ring',/)
            do i = 1, nring5
               write (iout,140)  i,(iring5(j,i),j=1,5)
  140          format (9x,i5,10x,5i6)
            end do
         end if
         if (nring6 .gt. 0) then
            write (iout,150)
  150       format (/,' Six-Membered Rings Contained',
     &                 ' in the Structure :',
     &              //,11x,'Ring',23x,'Atoms in Ring',/)
            do i = 1, nring6
               write (iout,160)  i,(iring6(j,i),j=1,6)
  160          format (9x,i5,10x,6i6)
            end do
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
c     ###########################################################
c     ##                                                       ##
c     ##  function rmsfit  --  rms deviation for paired atoms  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "rmsfit" computes the rms fit of two coordinate sets
c
c
      function rmsfit (x1,y1,z1,x2,y2,z2)
      implicit none
      include 'sizes.i'
      include 'align.i'
      integer i,i1,i2
      real*8 rmsfit,rmsterm
      real*8 xr,yr,zr,dist2
      real*8 weight,norm
      real*8 x1(maxatm),y1(maxatm),z1(maxatm)
      real*8 x2(maxatm),y2(maxatm),z2(maxatm)
c
c
c     compute the rms fit over superimposed atom pairs
c
      rmsfit = 0.0d0
      norm = 0.0d0
      do i = 1, nfit
         i1 = ifit(1,i)
         i2 = ifit(2,i)
         weight = wfit(i)
         xr = x1(i1) - x2(i2)
         yr = y1(i1) - y2(i2)
         zr = z1(i1) - z2(i2)
         dist2 = xr**2 + yr**2 + zr**2
         norm = norm + weight
         rmsterm = dist2 * weight
         rmsfit = rmsfit + rmsterm
      end do
      rmsfit = sqrt(rmsfit/norm)
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
c     ##  subroutine rotlist  --  find atoms on one side of a bond  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotlist" generates the minimum list of all the atoms lying
c     to one side of a pair of directly bonded atoms
c
c
      subroutine rotlist (base,partner)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'rotate.i'
      include 'zclose.i'
      integer i,k,ia,ib
      integer base,partner
      integer mark,test,nattach
      integer list(0:maxatm)
      logical bonded
c
c
c     initialize the number of atoms to one side of the bond
c
      nrot = 0
c
c     remove any bonds needed for intramolecular ring closures
c
      do i = 1, nadd
         ia = iadd(1,i)
         ib = iadd(2,i)
         if (molcule(ia) .eq. molcule(ib)) then
            do k = 1, n12(ia)
               if (i12(k,ia) .eq. ib)  i12(k,ia) = 0
            end do
            do k = 1, n12(ib)
               if (i12(k,ib) .eq. ia)  i12(k,ib) = 0
            end do
         end if
      end do
c
c     add any links needed to make intermolecular connections
c
      do i = 1, ndel
         ia = idel(1,i)
         ib = idel(2,i)
         if (molcule(ia) .ne. molcule(ib)) then
            if (n12(ia).eq.maxval .or. n12(ib).eq.maxval) then
               write (iout,10)
   10          format (/,' ROTLIST  --  Maximum Valence Exceeded;',
     &                    ' Increase MAXVAL')
               call fatal
            end if
            n12(ia) = n12(ia) + 1
            i12(n12(ia),ia) = ib
            n12(ib) = n12(ib) + 1
            i12(n12(ib),ib) = ia
         end if
      end do
c
c     check to see if the two atoms are still directly bonded
c
      bonded = .false.
      do i = 1, n12(base)
         if (i12(i,base) .eq. partner)  bonded = .true.
      end do
c
c     make a list of atoms to one side of this pair of atoms,
c     taking note of any rings in which the atom pair resides
c
      if (bonded) then
         list(0) = 1
         do i = 1, n
            rot(i) = 0
         end do
         nrot = 0
         do i = 1, n
            list(i) = 0
         end do
         list(base) = 1
         list(partner) = 1
         nattach = n12(base)
         do i = 1, nattach
            test = i12(i,base)
            if (list(test) .eq. 0) then
               nrot = nrot + 1
               rot(nrot) = test
               list(test) = 1
            end if
         end do
         do i = 1, n
            mark = rot(i)
            if (mark .eq. 0)  goto 20
            nattach = n12(mark)
            if (nattach .gt. 1) then
               do k = 1, nattach
                  test = i12(k,mark)
                  if (list(test) .eq. 0) then
                     nrot = nrot + 1
                     rot(nrot) = test
                     list(test) = 1
                  end if
               end do
            end if
         end do
      end if
   20 continue
c
c     remove links added to make intermolecular connections
c
      do i = 1, ndel
         ia = idel(1,i)
         ib = idel(2,i)
         if (molcule(ia) .ne. molcule(ib)) then
            n12(ia) = n12(ia) - 1
            n12(ib) = n12(ib) - 1
         end if
      end do
c
c     add any bonds required for intramolecular ring closures
c
      do i = 1, nadd
         ia = iadd(1,i)
         ib = iadd(2,i)
         if (molcule(ia) .eq. molcule(ib)) then
            do k = 1, n12(ia)
               if (i12(k,ia) .eq. 0) then
                  i12(k,ia) = ib
                  goto 30
               end if
            end do
   30       continue
            do k = 1, n12(ib)
               if (i12(k,ib) .eq. 0) then
                  i12(k,ib) = ia
                  goto 40
               end if
            end do
   40       continue
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole  --  rotate multipoles to global frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotpole" computes the atomic multipole values in the global
c     coordinate frame by applying a rotation matrix to a set of
c     locally defined multipoles
c
c
      subroutine rotpole (imdq,a)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k,m,imdq
      real*8 m2(3,3),r2(3,3),a(3,3)
c
c
c     monopoles have the same value in any coordinate frame
c
      rpole(1,imdq) = pole(1,imdq)
c
c     rotate the dipoles to the global coordinate frame
c
      do i = 2, 4
         rpole(i,imdq) = 0.0d0
         do j = 2, 4
            rpole(i,imdq) = rpole(i,imdq) + pole(j,imdq)*a(i-1,j-1)
         end do
      end do
c
c     rotate the quadrupoles to the global coordinate frame
c
      k = 5
      do i = 1, 3
         do j = 1, 3
            m2(i,j) = pole(k,imdq)
            r2(i,j) = 0.0d0
            k = k + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            if (j .lt. i) then
               r2(i,j) = r2(j,i)
            else
               do k = 1, 3
                  do m = 1, 3
                     r2(i,j) = r2(i,j) + a(i,k)*a(j,m)*m2(k,m)
                  end do
               end do
            end if
         end do
      end do
      k = 5
      do i = 1, 3
         do j = 1, 3
            rpole(k,imdq) = r2(i,j)
            k = k + 1
         end do
      end do
      return
      end
c
c
c     ###########################
c     ##                       ##
c     ##  subroutine drotpole  ##
c     ##                       ##
c     ###########################
c
c
c     "drotpole" computes the derivatives of the atomic multipoles
c     in the global coordinate frame with respect to motion of the
c     sites that define the local coordinate frame
c
c
      subroutine drotpole (imdq,a,d,p,q)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k,m,nn
      integer imdq,p,q
      real*8 a(3,3),d(3,3,3,3)
      real*8 m2(3,3),dm2(3,3)
c
c
c     derivative for rotation of the monopole term is zero
c
      dpole(1,p,q,imdq) = 0.0d0
c
c     derivative terms for rotation of the dipole components
c
      do i = 2, 4
         dpole(i,p,q,imdq) = 0.0d0
         do j = 2, 4
            dpole(i,p,q,imdq) = dpole(i,p,q,imdq)
     &                              + pole(j,imdq)*d(i-1,j-1,p,q)
         end do
      end do
c
c     derivative terms for rotation of the quadrupole components
c
      nn = 5
      do i = 1, 3
         do j = 1, 3
            m2(i,j) = pole(nn,imdq)
            dm2(i,j) = 0.0d0
            nn = nn + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               do m = 1, 3
                  dm2(i,j) = dm2(i,j) + m2(k,m) *
     &                         (d(i,k,p,q)*a(j,m)+a(i,k)*d(j,m,p,q))
               end do
            end do
          end do
      end do
      nn = 5
      do i = 1, 3
         do j = 1, 3
            dpole(nn,p,q,imdq) = dm2(i,j)
            nn = nn + 1
         end do
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine rotmat  --  find multipole rotation matrix  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "rotmat" find the rotation matrix that converts from the local
c     coordinate system at each multipole site to the global system
c
c
cjrs
      subroutine rotmatt (i,a)
cjrs
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i
      real*8 r,dotxz
      real*8 dx,dy,dz
      real*8 dx1,dy1,dz1
      real*8 dx2,dy2,dz2
      real*8 a(3,3)
c
c
c     rotation matrix elements for z- and x-axes, first z then x
c
      if (polaxe(i) .eq. 'Z-then-X') then
         dx = x(zaxis(i)) - x(ipole(i))
         dy = y(zaxis(i)) - y(ipole(i))
         dz = z(zaxis(i)) - z(ipole(i))
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dx = x(xaxis(i)) - x(ipole(i))
         dy = y(xaxis(i)) - y(ipole(i))
         dz = z(xaxis(i)) - z(ipole(i))
         dotxz = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
         dx = dx - dotxz*a(1,3)
         dy = dy - dotxz*a(2,3)
         dz = dz - dotxz*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
c
c     rotation matrix elements for z- and x-axes, bisector method
c
      else if (polaxe(i) .eq. 'Bisector') then
         dx = x(zaxis(i)) - x(ipole(i))
         dy = y(zaxis(i)) - y(ipole(i))
         dz = z(zaxis(i)) - z(ipole(i))
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx1 = dx / r
         dy1 = dy / r
         dz1 = dz / r
         dx = x(xaxis(i)) - x(ipole(i))
         dy = y(xaxis(i)) - y(ipole(i))
         dz = z(xaxis(i)) - z(ipole(i))
         r = sqrt(dx*dx + dy*dy + dz*dz)
         dx2 = dx / r
         dy2 = dy / r
         dz2 = dz / r
         dx = dx1 + dx2
         dy = dy1 + dy2
         dz = dz1 + dz2
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,3) = dx / r
         a(2,3) = dy / r
         a(3,3) = dz / r
         dotxz = dx2*a(1,3) + dy2*a(2,3) + dz2*a(3,3)
         dx = dx2 - dotxz*a(1,3)
         dy = dy2 - dotxz*a(2,3)
         dz = dz2 - dotxz*a(3,3)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         a(1,1) = dx / r
         a(2,1) = dy / r
         a(3,1) = dz / r
      end if
c
c     finally, find rotation matrix elements for the y-axis
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine drotmat  --  multipole rotation matrix derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "drotmat" finds the derivative rotation matrices that convert
c     multipoles from the local coordinate system to the global system
c
c
      subroutine drotmat (i,d)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      integer i
      real*8 d(3,3,3,3)
c
c
      if (polaxe(i) .eq. 'Z-then-X') then
         call drotmat1 (i,d)
      else if (polaxe(i) .eq. 'Bisector') then
         call drotmat2 (i,d)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine drotmat1  --  Z-then-X local coordinate derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "drotmat1" finds the multipole rotation matrix derivatives
c     for local coordinates defined via the "Z-then-X" method
c
c
      subroutine drotmat1 (i,d)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k
      real*8 d(3,3,3,3)
      real*8 xi,yi,zi,xz,yz,zz,xx,yx,zx
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9
      real*8 t10,t12,t13,t14,t15,t16,t17,t18,t19
      real*8 t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
      real*8 t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
      real*8 t40,t41,t42,t43,t44,t46,t47,t49
      real*8 t50,t52,t58,t59
      real*8 t60,t62,t64,t65,t66,t68,t69
      real*8 t70,t72,t73,t74,t75,t76,t77,t78
      real*8 t80,t81,t86,t91,t93,t95,t97,t98
      real*8 t101,t105,t106,t111,t116,t118
      real*8 t120,t123,t124,t126,t130,t131,t136
      real*8 t141,t143,t145,t146,t148
      real*8 t150,t152,t153,t155,t157,t159
      real*8 t163,t167,t168,t173,t178,t180,t182,t186
      real*8 t190,t191,t196,t201,t203,t205,t209,t213,t214
      real*8 t219,t224,t226,t228,t230,t232,t234,t236,t238
      real*8 t240,t242,t243,t246,t249,t251,t253,t254,t255
      real*8 t256,t259,t263,t265,t266,t271,t273,t276,t278
      real*8 t282,t285,t287,t290,t292,t295,t297,t300,t304
      real*8 t308,t310,t314,t317,t324,t329,t332,t354,t358
      real*8 t362,t366,t370,t374,t381,t388,t395,t416,t418
      real*8 t422,t425,t428,t430,t432,t436,t438,t441,t446
c
c
c     coordinates of main atom and those defining local axes
c
      xi = x(ipole(i))
      yi = y(ipole(i))
      zi = z(ipole(i))
      xz = x(zaxis(i))
      yz = y(zaxis(i))
      zz = z(zaxis(i))
      xx = x(xaxis(i))
      yx = y(xaxis(i))
      zx = z(xaxis(i))
c
c     temporary variables from Maple symbolic algebra derivation
c
      t1 = xz - xi
      t2 = t1 * t1
      t3 = yz - yi
      t4 = t3 * t3
      t5 = zz - zi
      t6 = t5 * t5
      t7 = t2 + t4 + t6
      t8 = sqrt(t7)
      t9 = 1.0d0 / t8
      t10 = t7 * t7
      t12 = t8 / t10
      t13 = t1 * t12
      t14 = 2.0d0 * (xi-xz)
      t15 = t13 * t14
      t16 = 2.0d0 * (yi-yz)
      t17 = t13 * t16
      t18 = 2.0d0 * (zi-zz)
      t19 = t13 * t18
      t20 = t3 * t12
      t21 = t20 * t14
      t22 = t20 * t16
      t23 = t20 * t18
      t24 = t5 * t12
      t25 = t24 * t14
      t26 = t24 * t16
      t27 = t24 * t18
      t28 = 2.0d0 * (xz-xi)
      t29 = t13 * t28
      t30 = 2.0d0 * (yz-yi)
      t31 = t13 * t30
      t32 = 2.0d0 * (zz-zi)
      t33 = t13 * t32
      t34 = t20 * t28
      t35 = t20 * t30
      t36 = t20 * t32
      t37 = t24 * t28
      t38 = t24 * t30
      t39 = t24 * t32
      t40 = t1 * t9
      t41 = xx - xi
      t42 = t41 * t9
      t43 = t41 * t1
      t44 = t12 * t14
      t46 = yx - yi
      t47 = t46 * t3
      t49 = zx - zi
      t50 = t49 * t5
      t52 = -t40 - t42 - 0.5d0*t44*(t43+t47+t50)
      t58 = t9 * (t43+t47+t50)
      t59 = t58 * t9
      t60 = t58 * t1
      t62 = -1.0d0 - t52*t1*t9 + t59 + 0.5d0*t60*t44
      t64 = xx - xi - t60*t9
      t65 = t64 * t64
      t66 = t58 * t3
      t68 = yx - yi - t66*t9
      t69 = t68 * t68
      t70 = t58 * t5
      t72 = zx - zi - t70*t9
      t73 = t72 * t72
      t74 = t65 + t69 + t73
      t75 = sqrt(t74)
      t76 = 1.0d0 / t75
      t77 = t62 * t76
      t78 = t74 * t74
      t80 = t75 / t78
      t81 = t64 * t80
      t86 = -t52*t3*t9 + 0.5d0*t66*t44
      t91 = -t52*t5*t9 + 0.5d0*t70*t44
      t93 = 2.0d0 * (t64*t62 + t68*t86 + t72*t91)
      t95 = t12 * t16
      t97 = t3 * t9
      t98 = t46 * t9
      t101 =  -t97 - t98 - 0.5d0*t95*(t43+t47+t50)
      t105 = -t101*t1*t9 + 0.5d0*t60*t95
      t106 = t105 * t76
      t111 = -1.0d0 - t101*t3*t9 + t59 + 0.5d0*t66*t95
      t116 = -t101*t5*t9 + 0.5d0*t70*t95
      t118 = 2.0d0 * (t64*t105 + t68*t111 + t72*t116)
      t120 = t12 * t18
      t123 = t5 * t9
      t124 = t49 * t9
      t126 = - t123 - t124 - 0.5d0*t120*(t43+t47+t50)
      t130 = -t126*t1*t9 + 0.5d0*t60*t120
      t131 = t130 * t76
      t136 = -t126*t3*t9 + 0.5d0*t66*t120
      t141 = -1.0d0 - t126*t5*t9 + t59 + 0.5d0*t70*t120
      t143 = 2.0d0 * (t64*t130 + t68*t136 + t72*t141)
      t145 = t86 * t76
      t146 = t68 * t80
      t148 = t111 * t76
      t150 = t136 * t76
      t152 = t91 * t76
      t153 = t72 * t80
      t155 = t116 * t76
      t157 = t141 * t76
      t159 = t12 * t28
      t163 = t42 - 0.5d0*t159*(t43+t47+t50)
      t167 = -t163*t1*t9 - t59 + 0.5d0*t60*t159
      t168 = t167 * t76
      t173 = -t163*t3*t9 + 0.5d0*t66*t159
      t178 = -t163*t5*t9 + 0.5d0*t70*t159
      t180 = 2.0d0 * (t64*t167 + t68*t173 + t72*t178)
      t182 = t12 * t30
      t186 = t98 - 0.5d0*t182*(t43+t47+t50)
      t190 = -t186*t1*t9 + 0.5d0*t60*t182
      t191 = t190 * t76
      t196 = -t186*t3*t9 - t59 + 0.5d0*t66*t182
      t201 = -t186*t5*t9 + 0.5d0*t70*t182
      t203 = 2.0d0 * (t64*t190 + t68*t196 + t72*t201)
      t205 = t12 * t32
      t209 = t124 - 0.5d0*t205*(t43+t47+t50)
      t213 = -t209*t1*t9 + 0.5d0*t60*t205
      t214 = t213 * t76
      t219 = -t209*t3*t9 + 0.5d0*t66*t205
      t224 = -t209*t5*t9 - t59 + 0.5d0*t70*t205
      t226 = 2.0d0 * (t64*t213 + t68*t219 + t72*t224)
      t228 = t173 * t76
      t230 = t196 * t76
      t232 = t219 * t76
      t234 = t178 * t76
      t236 = t201 * t76
      t238 = t224 * t76
      t240 = 1.0d0 / t7
      t242 = 1.0d0 - t2*t240
      t243 = t242 * t76
      t246 = t240 * t3
      t249 = t240 * t5
      t251 = 2.0d0 * (t64*t242 - t68*t1*t246 - t72*t1*t249)
      t253 = t1 * t240
      t254 = t3 * t76
      t255 = t253 * t254
      t256 = t64 * t1
      t259 = 1.0d0 - t4*t240
      t263 =  2.0d0 * (t68*t259 - t256*t246 - t72*t3*t249)
      t265 = t5 * t76
      t266 = t253 * t265
      t271 = 1.0d0 - t6*t240
      t273 =  2.0d0 * (t72*t271 - t256*t249 - t68*t3*t249)
      t276 = t259 * t76
      t278 = t246 * t265
      t282 = t271 * t76
      t285 = t97 * t93
      t287 = t72 * t76
      t290 = t123 * t93
      t292 = t68 * t76
      t295 = t97 * t118
      t297 = t287 * t9
      t300 = t123 * t118
      t304 = t97 * t143
      t308 = t123 * t143
      t310 = t292 * t9
      t314 = t64 * t76
      t317 = t40 * t93
      t324 = t40 * t118
      t329 = t314 * t9
      t332 = t40 * t143
      t354 = t97 * t180
      t358 = t123 * t180
      t362 = t97 * t203
      t366 = t123 * t203
      t370 = t97 * t226
      t374 = t123 * t226
      t381 = t40 * t180
      t388 = t40 * t203
      t395 = t40 * t226
      t416 = t97 * t251
      t418 = t123 * t251
      t422 = t97 * t263
      t425 = t123 * t263
      t428 = t97 * t273
      t430 = t6 * t76
      t432 = t123 * t273
      t436 = t2 * t12
      t438 = t40 * t251
      t441 = t40 * t263
      t446 = t40 * t273
c
c     set values for the derivatives of the rotation matrix
c
      d(1,3,1,1) = -t9 - 0.5d0*t15
      d(1,3,1,2) = -0.5d0 * t17
      d(1,3,1,3) = -0.5d0 * t19
      d(2,3,1,1) = -0.5d0 * t21
      d(2,3,1,2) = -t9 - 0.5d0*t22
      d(2,3,1,3) = -0.5d0 * t23
      d(3,3,1,1) = -0.5d0 * t25
      d(3,3,1,2) = -0.5d0 * t26
      d(3,3,1,3) = -t9 - 0.5d0*t27
      d(1,3,2,1) = t9 - 0.5d0*t29
      d(1,3,2,2) = -0.5d0 * t31
      d(1,3,2,3) = -0.5d0 * t33
      d(2,3,2,1) = -0.5d0 * t34
      d(2,3,2,2) = t9 - 0.5d0*t35
      d(2,3,2,3) = -0.5d0 * t36
      d(3,3,2,1) = -0.5d0 * t37
      d(3,3,2,2) = -0.5d0 * t38
      d(3,3,2,3) = t9 - 0.5d0*t39
      d(1,1,1,1) = t77 - 0.5d0*t81*t93
      d(1,1,1,2) = t106 - 0.5d0*t81*t118
      d(1,1,1,3) = t131 - 0.5d0*t81*t143
      d(2,1,1,1) = t145 - 0.5d0*t146*t93
      d(2,1,1,2) = t148 - 0.5d0*t146*t118
      d(2,1,1,3) = t150 - 0.5d0*t146*t143
      d(3,1,1,1) = t152 - 0.5d0*t153*t93
      d(3,1,1,2) = t155 - 0.5d0*t153*t118
      d(3,1,1,3) = t157 - 0.5d0*t153*t143
      d(1,1,2,1) = t168 - 0.5d0*t81*t180
      d(1,1,2,2) = t191 - 0.5d0*t81*t203
      d(1,1,2,3) = t214 - 0.5d0*t81*t226
      d(2,1,2,1) = t228 - 0.5d0*t146*t180
      d(2,1,2,2) = t230 - 0.5d0*t146*t203
      d(2,1,2,3) = t232 - 0.5d0*t146*t226
      d(3,1,2,1) = t234 - 0.5d0*t153*t180
      d(3,1,2,2) = t236 - 0.5d0*t153*t203
      d(3,1,2,3) = t238 - 0.5d0*t153*t226
      d(1,1,3,1) = t243 - 0.5d0*t81*t251
      d(1,1,3,2) = -t255 - 0.5d0*t81*t263
      d(1,1,3,3) = -t266 - 0.5d0*t81*t273
      d(2,1,3,1) = -t255 - 0.5d0*t146*t251
      d(2,1,3,2) = t276 - 0.5d0*t146*t263
      d(2,1,3,3) = -t278 - 0.5d0*t146*t273
      d(3,1,3,1) = -t266 - 0.5d0*t153*t251
      d(3,1,3,2) = -t278 - 0.5d0*t153*t263
      d(3,1,3,3) = t282 - 0.5d0*t153*t273
      d(1,2,1,1) = t152*t97 - 0.5d0*(t153*t285 + t287*t21)
     &                - t145*t123 + 0.5d0*(t146*t290 + t292*t25)
      d(1,2,1,2) = t155*t97 - t297 - 0.5d0*(t153*t295 + t287*t22)
     &                - t148*t123 + 0.5d0*(t146*t300 + t292*t26)
      d(1,2,1,3) = t157*t97 + t310 - 0.5d0*(t153*t304 + t287*t23)
     &                - t150*t123 + 0.5d0*(t146*t308 + t292*t27)
      d(2,2,1,1) = t77*t123 + t297 - 0.5d0*(t81*t290 + t314*t25)
     &                - t152*t40 + 0.5d0*(t153*t317 + t287*t15)
      d(2,2,1,2) = t106*t123 - 0.5d0*(t81*t300 + t314*t26)
     &                - t155*t40 + 0.5d0*(t153*t324 + t287*t17)
      d(2,2,1,3) = t131*t123 - t329 - 0.5d0*(t81*t308 + t314*t27)
     &                - t157*t40 + 0.5d0*(t153*t332 + t287*t19)
      d(3,2,1,1) = t145*t40 - t310 - 0.5d0*(t146*t317 + t292*t15)
     &                - t77*t97 + 0.5d0*(t81*t285 + t314*t21)
      d(3,2,1,2) = t148*t40 + t329 - 0.5d0*(t146*t324 + t292*t17)
     &                - t106*t97 + 0.5d0*(t81*t295 + t314*t22)
      d(3,2,1,3) = t150*t40 - 0.5d0*(t146*t332 + t292*t19)
     &                - t131*t97 + 0.5d0*(t81*t304 + t314*t23)
      d(2,2,1,3) = t131*t123 - t329 - 0.5d0*(t81*t308 + t314*t27)
     &                - t157*t40 + 0.5d0*(t153*t332 + t287*t19)
      d(3,2,1,1) = t145*t40 - t310 - 0.5d0*(t146*t317 + t292*t15)
     &                - t77*t97 + 0.5d0*(t81*t285 + t314*t21)
      d(3,2,1,2) = t148*t40 + t329 - 0.5d0*(t146*t324 + t292*t17)
     &                - t106*t97 + 0.5d0*(t81*t295 + t314*t22)
      d(3,2,1,3) = t150*t40 - 0.5d0*(t146*t332 + t292*t19)
     &                - t131*t97 + 0.5d0*(t81*t304 + t314*t23)
      d(1,2,2,1) = t234*t97 - 0.5d0*(t153*t354 + t287*t34)
     &                - t228*t123 + 0.5d0*(t146*t358 + t292*t37)
      d(1,2,2,2) = t236*t97 + t297 - 0.5d0*(t153*t362 + t287*t35)
     &                - t230*t123 + 0.5d0*(t146*t366 + t292*t38)
      d(1,2,2,3) = t238*t97 - t310 - 0.5d0*(t153*t370 + t287*t36)
     &                - t232*t123 + 0.5d0*(t146*t374 + t292*t39)
      d(2,2,2,1) = t168*t123 - t297 - 0.5d0*(t81*t358 + t314*t37)
     &                - t234*t40 + 0.5d0*(t153*t381 + t287*t29)
      d(2,2,2,2) = t191*t123 - 0.5d0*(t81*t366 + t314*t38)
     &                - t236*t40 + 0.5d0*(t153*t388 + t287*t31)
      d(2,2,2,3) = t214*t123 + t329 - 0.5d0*(t81*t374 + t314*t39)
     &                - t238*t40 + 0.5d0*(t153*t395 + t287*t33)
      d(3,2,2,1) = t228*t40 + t310 - 0.5d0*(t146*t381 + t292*t29)
     &                - t168*t97 + 0.5d0*(t81*t354 + t314*t34)
      d(3,2,2,2) = t230*t40 - t329 - 0.5d0*(t146*t388 + t292*t31)
     &                - t191*t97 + 0.5d0*(t81*t362 + t314*t35)
      d(3,2,2,3) = t232*t40 - 0.5d0*(t146*t395 + t292*t33)
     &                - t214*t97 + 0.5d0*(t81*t370 + t314*t36)
      d(1,2,3,1) = 0.5d0 * (t146*t418 - t153*t416)
      d(1,2,3,2) = -t4*t12*t265 - t276*t123
     &                + 0.5d0*(t146*t425 - t153*t422)
      d(1,2,3,3) = t282*t97 + t20*t430 + 0.5d0*(t146*t432-t153*t428)
      d(2,2,3,1) = t243*t123 + t436*t265 + 0.5d0*(t153*t438-t81*t418)
      d(2,2,3,2) = 0.5d0 * (t153*t441 - t81*t425)
      d(2,2,3,3) = -t13*t430 - t282*t40 + 0.5d0*(t153*t446-t81*t432)
      d(3,2,3,1) = -t436*t254 - t243*t97 + 0.5d0*(t81*t416-t146*t438)
      d(3,2,3,2) = t276*t40 + t13*t4*t76 + 0.5d0*(t81*t422-t146*t441)
      d(3,2,3,3) = 0.5d0 * (t81*t428 - t146*t446)
c
c     some of the derivative matrix values are always zero
c
      do j = 1, 3
         do k = 1, 3
            d(k,3,3,j) = 0.0d0
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine drotmat2  --  bisector local coordinate derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "drotmat2" finds the multipole rotation matrix derivatives
c     for local coordinates defined via the "Bisector" method
c
c
      subroutine drotmat2 (i,d)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k,m
      real*8 d(3,3,3,3)
      real*8 xi,yi,zi,xz,yz,zz,xx,yx,zx
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      real*8 t12,t13,t14,t16,t18,t19,t20
      real*8 t21,t22,t23,t24,t25,t26,t27,t28,t29
      real*8 t31,t32,t33,t35,t36,t37,t38,t39,t40
      real*8 t41,t42,t44,t45,t47,t48,t50
      real*8 t52,t53,t54,t55,t57,t58,t60
      real*8 t62,t65,t66,t67,t68,t70
      real*8 t73,t75,t77,t78,t79
      real*8 t81,t82,t83,t84,t86,t88,t89
      real*8 t91,t92,t93,t94,t95,t96,t97,t98,t99,t100
      real*8 t101,t103,t104,t106,t108,t109,t110
      real*8 t111,t113,t114,t116,t117,t118
      real*8 t121,t122,t123,t124,t126,t129,t130
      real*8 t131,t133,t134,t135,t137,t138,t139
      real*8 t141,t143,t145,t146,t147,t150
      real*8 t151,t154,t157,t160
      real*8 t162,t164,t166,t169,t170
      real*8 t172,t174,t175,t176,t178,t179,t180
      real*8 t182,t183,t184,t185,t186,t187,t188,t190
      real*8 t191,t195,t198,t202,t205,t207,t209,t210
      real*8 t217,t220,t222,t224,t225,t232,t238,t240
      real*8 t242,t249,t253,t255,t256,t262,t269
      real*8 t271,t273,t274,t276,t278,t280
      real*8 t281,t283,t285,t295,t300,t303,t308,t310
      real*8 t311,t316,t317,t320,t325,t328,t330
      real*8 t335,t348,t351,t352,t354,t356,t357
      real*8 t364,t371,t373,t390,t393,t395,t397,t398
      real*8 t405,t412,t414,t416,t418,t420
      real*8 t422,t424,t426,t429
      real*8 t431,t432,t436,t438,t439
      real*8 t443,t448,t453,t458,t464,t465,t469
      real*8 t478,t480,t488,t498,t503,t519
      real*8 t521,t525,t527,t531,t536,t541,t546
      real*8 t552,t556,t565,t567,t575,t585,t590
c
c
c     coordinates of main atom and those defining local axes
c
      xi = x(ipole(i))
      yi = y(ipole(i))
      zi = z(ipole(i))
      xz = x(zaxis(i))
      yz = y(zaxis(i))
      zz = z(zaxis(i))
      xx = x(xaxis(i))
      yx = y(xaxis(i))
      zx = z(xaxis(i))
c
c     temporary variables from Maple symbolic algebra derivation
c
      t1 = xz - xi
      t2 = t1**2
      t3 = yz - yi
      t4 = t3**2
      t5 = zz - zi
      t6 = t5**2
      t7 = t2 + t4 + t6
      t8 = dsqrt(t7)
      t9 = 1.0d0 / t8
      t10 = t7**2
      t12 = t8 / t10
      t13 = t1 * t12
      t14 = 2.0d0 * (xz-xi)
      t16 = t9 - 0.5d0*t13*t14
      t18 = xx - xi
      t19 = t18**2
      t20 = yx - yi
      t21 = t20**2
      t22 = zx - zi
      t23 = t22**2
      t24 = t19 + t21 + t23
      t25 = dsqrt(t24)
      t26 = 1.0d0 / t25
      t27 = t18 * t26
      t28 = t1*t9 + t27
      t29 = t28**2
      t31 = t20 * t26
      t32 = t3*t9 + t31
      t33 = t32**2
      t35 = t22 * t26
      t36 = t5*t9 + t35
      t37 = t36**2
      t38 = t29 + t33 + t37
      t39 = dsqrt(t38)
      t40 = 1.0d0 / t39
      t41 = t16 * t40
      t42 = t38**2
      t44 = t39 / t42
      t45 = t28 * t44
      t47 = t32 * t3
      t48 = t12 * t14
      t50 = t36 * t5
      t52 = 2.0d0*t28*t16 - t47*t48 - t50*t48
      t53 = t45 * t52
      t54 = 2.0d0 * (yz-yi)
      t55 = t54 * t40
      t57 = t28 * t1
      t58 = t12 * t54
      t60 = t3 * t12
      t62 = t9 - 0.5d0*t60*t54
      t65 = -t57*t58 + 2.0d0*t32*t62 - t50*t58
      t66 = t45 * t65
      t67 = 2.0d0 * (zz-zi)
      t68 = t67 * t40
      t70 = t12 * t67
      t73 = t5 * t12
      t75 = t9 - 0.5d0*t73*t67
      t77 = -t57*t70 - t47*t70 + 2.0d0*t36*t75
      t78 = t45 * t77
      t79 = t14 * t40
      t81 = t32 * t44
      t82 = t81 * t52
      t83 = t62 * t40
      t84 = t81 * t65
      t86 = t81 * t77
      t88 = t36 * t44
      t89 = t88 * t52
      t91 = t88 * t65
      t92 = t75 * t40
      t93 = t88 * t77
      t94 = t24**2
      t95 = 1.0d0 / t94
      t96 = t25 * t95
      t97 = t18 * t96
      t98 = 2.0d0 * (xx-xi)
      t99 = t97 * t98
      t100 = t26 - 0.5d0*t99
      t101 = t100 * t40
      t103 = t32 * t20
      t104 = t96 * t98
      t106 = t36 * t22
      t108 = 2.0d0*t28*t100 - t103*t104 - t106*t104
      t109 = t45 * t108
      t110 = 2.0d0 * (yx-yi)
      t111 = t110 * t40
      t113 = t28 * t18
      t114 = t96 * t110
      t116 = t20 * t96
      t117 = t116 * t110
      t118 = t26 - 0.5d0*t117
      t121 = -t113*t114 + 2.0d0*t32*t118 - t106*t114
      t122 = t45 * t121
      t123 = 2.0d0 * (zx-zi)
      t124 = t123 * t40
      t126 = t96 * t123
      t129 = t22 * t96
      t130 = t129 * t123
      t131 = t26 - 0.5d0*t130
      t133 = -t113*t126 - t103*t126 + 2.0d0*t36*t131
      t134 = t45 * t133
      t135 = t98 * t40
      t137 = t81 * t108
      t138 = t118 * t40
      t139 = t81 * t121
      t141 = t81 * t133
      t143 = t88 * t108
      t145 = t88 * t121
      t146 = t131 * t40
      t147 = t88 * t133
      t150 = t31 * t3
      t151 = t48 * t40
      t154 = t35 * t5
      t157 = t27*t41 - 0.5d0*(t27*t53 + t150*t151 + t31*t82
     &                            + t154*t151 + t35*t89)
      t160 = t28 * t40
      t162 = t32 * t40
      t164 = t36 * t40
      t166 = t27*t160 + t31*t162 + t35*t164
      t169 = t166 * t28
      t170 = t44 * t52
      t172 = -t157*t28*t40 - t166*t16*t40 + 0.5d0*t169*t170
      t174 = t27 - t169*t40
      t175 = t174**2
      t176 = t166 * t32
      t178 = t31 - t176*t40
      t179 = t178**2
      t180 = t166 * t36
      t182 = t35 - t180*t40
      t183 = t182**2
      t184 = t175 + t179 + t183
      t185 = dsqrt(t184)
      t186 = 1.0d0 / t185
      t187 = t172 * t186
      t188 = t184**2
      t190 = t185 / t188
      t191 = t174 * t190
      t195 = t166 * t3
      t198 = -t157*t32*t40 + 0.5d0*(t195*t151 + t176*t170)
      t202 = t166 * t5
      t205 = -t157*t36*t40 + 0.5d0*(t202*t151 + t180*t170)
      t207 = 2.0d0 * (t174*t172 + t178*t198 + t182*t205)
      t209 = t27 * t1
      t210 = t58 * t40
      t217 = t31*t83 - 0.5d0*(t209*t210 + t27*t66 + t31*t84
     &                            + t154*t210 + t35*t91)
      t220 = t166 * t1
      t222 = t44 * t65
      t224 = -t217*t28*t40 + 0.5d0*(t220*t210 + t169*t222)
      t225 = t224 * t186
      t232 = -t217*t32*t40 - t166*t62*t40 + 0.5d0*t176*t222
      t238 = -t217*t36*t40 + 0.5d0*(t202*t210 + t180*t222)
      t240 = 2.0d0 * (t174*t224 + t178*t232 + t182*t238)
      t242 = t70 * t40
      t249 = t35*t92 - 0.5d0*(t209*t242 + t27*t78 + t150*t242
     &                              + t31*t86 + t35*t93)
      t253 = t44 * t77
      t255 = -t249*t28*t40 + 0.5d0*(t220*t242 + t169*t253)
      t256 = t255 * t186
      t262 = -t249*t32*t40 + 0.5d0*(t195*t242 + t176*t253)
      t269 = -t249*t36*t40 - t166*t75*t40 + 0.5d0*t180*t253
      t271 = 2.0d0 * (t174*t255 + t178*t262 + t182*t269)
      t273 = t198 * t186
      t274 = t178 * t190
      t276 = t232 * t186
      t278 = t262 * t186
      t280 = t205 * t186
      t281 = t182 * t190
      t283 = t238 * t186
      t285 = t269 * t186
      t295 = t21 * t95
      t300 = t23 * t95
      t303 = t26*t28*t40 + t27*t101
     &          - 0.5d0 * (t97*t160*t98 + t27*t109 + t116*t162*t98
     &                        + t295*t135 + t31*t137 + t300*t135
     &                        + t35*t143 + t129*t164*t98)
      t308 = t44 * t108
      t310 = t26 - t303*t28*t40 - t166*t100*t40
     &          + 0.5d0*(t169*t308 - t99)
      t311 = t310 * t186
      t316 = t166 * t20
      t317 = t104 * t40
      t320 = -t303*t32*t40 + 0.5d0*(t316*t317 + t176*t308 - t116*t98)
      t325 = t166 * t22
      t328 = -t303*t36*t40 + 0.5d0*(t325*t317 + t180*t308 - t129*t98)
      t330 = 2.0d0 * (t174*t310 + t178*t320 + t182*t328)
      t335 = t19 * t95
      t348 = t26*t32*t40 + t31*t138
     &          - 0.5d0 * (t97*t160*t110 + t335*t111 + t27*t122
     &                        + t116*t162*t110 + t129*t164*t110
     &                        + t300*t111 + t35*t145 + t31*t139)
      t351 = t166 * t18
      t352 = t114 * t40
      t354 = t44 * t121
      t356 = -t348*t28*t40 + 0.5d0*(t351*t352 + t169*t354 - t97*t110)
      t357 = t356 * t186
      t364 = t26 - t348*t32*t40 - t166*t118*t40
     &          + 0.5d0*(t176*t354 - t117)
      t371 = -t348*t36*t40 + 0.5d0*(t325*t352 + t180*t354 - t129*t110)
      t373 = 2.0d0 * (t174*t356 + t178*t364 + t182*t371)
      t390 = t26*t36*t40 + t35*t146
     &          - 0.5d0 * (t97*t160*t123 + t335*t124 + t27*t134
     &                        + t116*t162*t123 + t129*t164*t123
     &                        + t295*t124 + t31*t141 + t35*t147)
      t393 = t126 * t40
      t395 = t44 * t133
      t397 = -t390*t28*t40 + 0.5d0*(t351*t393 + t169*t395 - t97*t123)
      t398 = t397 * t186
      t405 = -t390*t32*t40 + 0.5d0*(t316*t393 + t176*t395 - t116*t123)
      t412 = t26 - t390*t36*t40 - t166*t131*t40
     &          + 0.5d0*(t180*t395 - t130)
      t414 = 2.0d0 * (t174*t397 + t178*t405 + t182*t412)
      t416 = t320 * t186
      t418 = t364 * t186
      t420 = t405 * t186
      t422 = t328 * t186
      t424 = t371 * t186
      t426 = t412 * t186
      t429 = t162 * t207
      t431 = t182 * t186
      t432 = t431 * t3
      t436 = t164 * t207
      t438 = t178 * t186
      t439 = t438 * t5
      t443 = t162 * t240
      t448 = t164 * t240
      t453 = t162 * t271
      t458 = t164 * t271
      t464 = t174 * t186
      t465 = t464 * t5
      t469 = t160 * t207
      t478 = t160 * t240
      t480 = t431 * t1
      t488 = t160 * t271
      t498 = t464 * t3
      t503 = t438 * t1
      t519 = t162 * t330
      t521 = t431 * t20
      t525 = t164 * t330
      t527 = t438 * t22
      t531 = t162 * t373
      t536 = t164 * t373
      t541 = t162 * t414
      t546 = t164 * t414
      t552 = t464 * t22
      t556 = t160 * t330
      t565 = t160 * t373
      t567 = t431 * t18
      t575 = t160 * t414
      t585 = t464 * t20
      t590 = t438 * t18
c
c     set values for the derivatives of the rotation matrix
c
      d(1,3,2,1) = t41 - 0.5d0*t53
      d(1,3,2,2) = -0.5d0 * (t13*t55 + t66)
      d(1,3,2,3) = -0.5d0 * (t13*t68 + t78)
      d(2,3,2,1) = -0.5d0 * (t60*t79 + t82)
      d(2,3,2,2) = t83 - 0.5d0*t84
      d(2,3,2,3) = -0.5d0 * (t60*t68 + t86)
      d(3,3,2,1) = -0.5d0 * (t73*t79 + t89)
      d(3,3,2,2) = -0.5d0 * (t73*t55 + t91)
      d(3,3,2,3) = t92 - 0.5d0*t93
      d(1,3,3,1) = t101 - 0.5d0*t109
      d(1,3,3,2) = -0.5d0 * (t97*t111 + t122)
      d(1,3,3,3) = -0.5d0 * (t97*t124 + t134)
      d(2,3,3,1) = -0.5d0 * (t116*t135 + t137)
      d(2,3,3,2) = t138 - 0.5d0*t139
      d(2,3,3,3) = -0.5d0 * (t116*t124 + t141)
      d(3,3,3,1) = -0.5d0 * (t129*t135 + t143)
      d(3,3,3,2) = -0.5d0 * (t129*t111 + t145)
      d(3,3,3,3) = t146 - 0.5d0*t147
      d(1,1,2,1) = t187 - 0.5d0*t191*t207
      d(1,1,2,2) = t225 - 0.5d0*t191*t240
      d(1,1,2,3) = t256 - 0.5d0*t191*t271
      d(2,1,2,1) = t273 - 0.5d0*t274*t207
      d(2,1,2,2) = t276 - 0.5d0*t274*t240
      d(2,1,2,3) = t278 - 0.5d0*t274*t271
      d(3,1,2,1) = t280 - 0.5d0*t281*t207
      d(3,1,2,2) = t283 - 0.5d0*t281*t240
      d(3,1,2,3) = t285 - 0.5d0*t281*t271
      d(1,1,3,1) = t311 - 0.5d0*t191*t330
      d(1,1,3,2) = t357 - 0.5d0*t191*t373
      d(1,1,3,3) = t398 - 0.5d0*t191*t414
      d(2,1,3,1) = t416 - 0.5d0*t274*t330
      d(2,1,3,2) = t418 - 0.5d0*t274*t373
      d(2,1,3,3) = t420 - 0.5d0*t274*t414
      d(3,1,3,1) = t422 - 0.5d0*t281*t330
      d(3,1,3,2) = t424 - 0.5d0*t281*t373
      d(3,1,3,3) = t426 - 0.5d0*t281*t414
      d(1,2,2,1) = t280*t162 - t273*t164
     &                - 0.5d0*(t281*t429 + t432*t151 + t431*t82)
     &                + 0.5d0*(t274*t436 + t439*t151 + t438*t89)
      d(1,2,2,2) = t283*t162 + t431*t83 - t276*t164
     &                - 0.5d0*(t281*t443 + t431*t84)
     &                + 0.5d0*(t274*t448 + t439*t210 + t438*t91)
      d(1,2,2,3) = t285*t162 - t278*t164 - t438*t92
     &                - 0.5d0*(t281*t453 + t432*t242 + t431*t86)
     &                + 0.5d0*(t274*t458 + t438*t93)
      d(2,2,2,1) = t187*t164 - t280*t160 - t431*t41
     &                - 0.5d0*(t191*t436 + t465*t151 + t464*t89)
     &                + 0.5d0*(t281*t469 + t431*t53)
      d(2,2,2,2) = t225*t164 - t283*t160
     &                - 0.5d0*(t191*t448 + t465*t210 + t464*t91)
     &                + 0.5d0*(t281*t478 + t480*t210 + t431*t66)
      d(2,2,2,3) = t256*t164 + t464*t92 - t285*t160
     &                - 0.5d0*(t191*t458 + t464*t93)
     &                + 0.5d0*(t281*t488 + t480*t242 + t431*t78)
      d(3,2,2,1) = t273*t160 + t438*t41 - t187*t162
     &                - 0.5d0*(t274*t469 + t438*t53)
     &                + 0.5d0*(t191*t429 + t498*t151 + t464*t82)
      d(3,2,2,2) = t276*t160 - t225*t162 - t464*t83
     &                - 0.5d0*(t274*t478 + t503*t210 + t438*t66)
     &                + 0.5d0*(t191*t443 + t464*t84)
      d(3,2,2,3) = t278*t160 - t256*t162
     &                - 0.5d0*(t274*t488 + t503*t242 + t438*t78)
     &                + 0.5d0*(t191*t453 + t498*t242 + t464*t86)
      d(1,2,3,1) = t422*t162 - t416*t164
     &                - 0.5d0*(t281*t519 + t521*t317 + t431*t137)
     &                + 0.5d0*(t274*t525 + t527*t317 + t438*t143)
      d(1,2,3,2) = t424*t162 + t431*t138 - t418*t164
     &                - 0.5d0*(t281*t531 + t431*t139)
     &                + 0.5d0*(t274*t536 + t527*t352 + t438*t145)
      d(1,2,3,3) = t426*t162 - t420*t164 - t438*t146
     &                - 0.5d0*(t281*t541 + t521*t393 + t431*t141)
     &                + 0.5d0*(t274*t546 + t438*t147)
      d(2,2,3,1) = t311*t164 - t422*t160 - t431*t101
     &                - 0.5d0*(t191*t525 + t552*t317 + t464*t143)
     &                + 0.5d0*(t281*t556 + t431*t109)
      d(2,2,3,2) = t357*t164 - t424*t160
     &                - 0.5d0*(t191*t536 + t552*t352 + t464*t145)
     &                + 0.5d0*(t281*t565 + t567*t352 + t431*t122)
      d(2,2,3,3) = t398*t164 + t464*t146 - t426*t160
     &                - 0.5d0*(t191*t546 + t464*t147)
     &                + 0.5d0*(t281*t575 + t567*t393 + t431*t134)
      d(3,2,3,1) = t416*t160 + t438*t101 - t311*t162
     &                - 0.5d0*(t274*t556 + t438*t109)
     &                + 0.5d0*(t191*t519 + t585*t317 + t464*t137)
      d(3,2,3,2) = t418*t160 - t357*t162 - t464*t138
     &                - 0.5d0*(t274*t565 + t590*t352 + t438*t122)
     &                + 0.5d0*(t191*t531 + t464*t139)
      d(3,2,3,3) = t420*t160 - t398*t162
     &                - 0.5d0*(t274*t575 + t590*t393 + t438*t134)
     &                + 0.5d0*(t191*t541 + t585*t393 + t464*t141)
c
c     some derivative matrix values are combinations of others
c
      do j = 1, 3
         do k = 1, 3
            do m = 1, 3
               d(m,k,1,j) = -d(m,k,2,j) - d(m,k,3,j)
            end do
         end do
      end do
      return
      end
