c  6 May 10 - DGF - add locase 
c  5 May 98 - JRS - TINKER library routines that start with
c                   s to x Lib(rary)t(inker)sx     
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine search  --  performs line search for minimum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "search" is a line search minimizer based upon parabolic
c     extrapolation and cubic interpolation using both function
c     and gradient values; if forced to search in an uphill
c     direction, return is after the initial step
c
c     variables used by the routine :
c
c     f       function value at the best line search point
c     x       current values of variables during line search
c     g       gradient at the current point during line search
c     p       initial search vector, unchanged by this routine
c     s       scaled search vector at current line search point
c     angle   angle between search and negative gradient vector
c
c     parameters used by the routine :
c
c     cappa    reduction in projected gradient for termination
c     stpmin   minimum allowed line search step size
c     stpmax   maximum allowed line search step size
c     angmax   maximum angle between search and -grad directions
c     intmax   maximum number of interpolations in line search
c
c     status codes upon return :
c
c     Success     normal termination after satisfying "cappa" test
c     ReSearch    normal termination after a reinterpolation
c     WideAngle   large angle between search and -grad directions
c     BadIntpln   unsatisfied "cappa" test after two searches
c     IntplnErr   function value increase or serious gradient error
c
c
      subroutine search (nvar,f,g,x,p,f_move,angle,ncalls,
     &                          fgvalue,status)
      implicit none
      include 'sizes.i'
      include 'linmin.i'
      include 'math.i'
casa
      include 'restrn.i'
      integer   atm
casa
      integer i,nvar,ncalls,intpln
      real*8 x(maxvar),g(maxvar),p(maxvar),s(maxvar)
      real*8 fgvalue,f,f_move,s_norm,g_norm,cosang,angle
      real*8 step,parab,cube,cubstp,sss,ttt
      real*8 f_0,f_1,f_a,f_b,f_c,sg_0,sg_1,sg_a,sg_b,sg_c
      character*9 status,blank
      external fgvalue
      save blank
      data blank  / '         ' /
c
c
c     use default parameters for the line search if needed
c
      if (cappa .eq. 0.0d0)  cappa = 0.1d0
      if (stpmin .eq. 0.0d0)  stpmin = 1.0d-20
      if (stpmax .eq. 0.0d0)  stpmax = 2.0d0
      if (angmax .eq. 0.0d0)  angmax = 180.0d0
      if (intmax .eq. 0)  intmax = 5
c
c     copy the search direction into a new vector
c
      do i = 1, nvar
         s(i) = p(i)
      end do
c
c     compute the length of gradient and search direction
c
      g_norm = 0.0d0
      s_norm = 0.0d0
      do i = 1, nvar
         g_norm = g_norm + g(i)*g(i)
         s_norm = s_norm + s(i)*s(i)
      end do
      g_norm = sqrt(g_norm)
      s_norm = sqrt(s_norm)
c
c     store initial function, then normalize the
c     search vector and find directional gradient
c
      f_0 = f
      sg_0 = 0.0d0
      do i = 1, nvar
         s(i) = s(i) / s_norm
         sg_0 = sg_0 + s(i)*g(i)
      end do
c
c     check the angle between the search direction
c     and the negative gradient vector
c
      cosang = -sg_0 / g_norm
      cosang = min(1.0d0,max(-1.0d0,cosang))
      angle = radian * acos(cosang)
      if (angle .gt. angmax) then
         status = 'WideAngle'
         return
      end if
c
c     set the initial stepsize to the length of the passed
c     search vector, or based on previous function decrease
c
      step = 2.0d0 * abs(f_move/sg_0)
      step = min(step,s_norm)
      if (step .gt. stpmax)  step = stpmax
      if (step .lt. stpmin)  step = stpmin
c
c     beginning of the parabolic extrapolation procedure
c
   10 continue
      intpln = 0
      f_b = f_0
      sg_b = sg_0
c
c     replace last point by latest and take another step
c
   20 continue
      f_a = f_b
      sg_a = sg_b
      do i = 1, nvar
         x(i) = x(i) + step*s(i)
      end do
c
c     get new function and gradient, then test for termination
c
      ncalls = ncalls + 1
      f_b = fgvalue (x,g)
      sg_b = 0.0d0
casa
      do i = 1, npfix
         if (pfix(i) .ge. RESTRAINV) then
            atm = (ipfix(i) - 1) * 3
c       write(iout,*) 'atm and g is',  atm, g(atm+1), g(atm+2), g(atm+3)
            g(atm+1) = 0.0d0
            g(atm+2) = 0.0d0
            g(atm+3) = 0.0d0
         end if
      end do
casa
      do i = 1, nvar
         sg_b = sg_b + s(i)*g(i)
      end do
      if (abs(sg_b/sg_0).le.cappa .and. f_b.lt.f_a) then
         f = f_b
         if (status .eq. blank)  status = ' Success '
         return
      end if
c
c     if slope changes or function increases, start interpolation;
c     if the finite difference curvature is negative double the step,
c     else if  step < parabolic estimate < 4*step  use this estimate,
c     else truncate to step or 4*step, respectively
c
      if (sg_a*sg_b.lt.0.0d0 .or. f_a.lt.f_b)  goto 30
      step = 2.0d0 * step
      if (sg_b .gt. sg_a) then
         parab = (f_a-f_b) / (sg_b-sg_a)
         if (parab .gt. 2.0d0*step)  parab = 2.0d0 * step
         if (parab .lt. 0.5d0*step)  parab = 0.5d0 * step
         step = parab
      end if
      if (step .gt. stpmax)  step = stpmax
      goto 20
c
c     beginning of the cubic interpolation procedure
c
   30 continue
      intpln = intpln + 1
      sss = 3.0d0*(f_b-f_a)/step - sg_a - sg_b
      ttt = sss*sss - sg_a*sg_b
      if (ttt .lt. 0.0d0) then
         f = f_b
         status = 'IntplnErr'
         return
      end if
      ttt = sqrt(ttt)
      cube = step * (sg_b+ttt+sss)/(sg_b-sg_a+2.0d0*ttt)
      if (cube.lt.0.0d0 .or. cube.gt.step) then
         f = f_b
         status = 'IntplnErr'
         return
      end if
      do i = 1, nvar
         x(i) = x(i) - cube*s(i)
      end do
c
c     get new function and gradient, then test for termination
c
      ncalls = ncalls + 1
      f_c = fgvalue (x,g)
      sg_c = 0.0d0
casa
      do i = 1, npfix
         if (pfix(i) .ge. RESTRAINV) then
            atm = (ipfix(i) - 1) * 3
c       write(iout,*) 'atm and g is',  atm, g(atm+1), g(atm+2), g(atm+3)
            g(atm+1) = 0.0d0
            g(atm+2) = 0.0d0
            g(atm+3) = 0.0d0
         end if
      end do
casa
      do i = 1, nvar
         sg_c = sg_c + s(i)*g(i)
      end do
      if (abs(sg_c/sg_0) .le. cappa) then
         f = f_c
         if (status .eq. blank)  status = ' Success '
         return
      end if
c
c     get the next pair of bracketing points by replacing one
c     of the current brackets with the interpolated point
c
      if (f_c.le.f_a .or. f_c.le.f_b) then
         cubstp = min(abs(cube),abs(step-cube))
         if (cubstp.ge.stpmin .and. intpln.lt.intmax) then
c
c     if the current brackets have slopes of opposite sign,
c     then substitute the interpolated point for the bracket
c     point with slope of same sign as the interpolated point
c
            if (sg_a*sg_b .lt. 0.0d0) then
               if (sg_a*sg_c .lt. 0.0d0) then
                  f_b = f_c
                  sg_b = sg_c
                  step = step - cube
               else
                  f_a = f_c
                  sg_a = sg_c
                  step = cube
                  do i = 1, nvar
                     x(i) = x(i) + cube*s(i)
                  end do
               end if
c
c     if current brackets have slope of same sign, then replace
c     the far bracket if the interpolated point has a slope of
c     the opposite sign or a lower function value than the near
c     bracket, otherwise replace the far bracket point
c
            else
               if (sg_a*sg_c.lt.0.0d0 .or. f_a.le.f_c) then
                  f_b = f_c
                  sg_b = sg_c
                  step = step - cube
               else
                  f_a = f_c
                  sg_a = sg_c
                  step = cube
                  do i = 1, nvar
                     x(i) = x(i) + cube*s(i)
                  end do
               end if
            end if
            goto 30
         end if
      end if
c
c     interpolation has failed, reset to best current point
c
      f_1 = min(f_a,f_b,f_c)
      if (f_1 .eq. f_a) then
         sg_1 = sg_a
         do i = 1, nvar
            x(i) = x(i) + (cube-step)*s(i)
         end do
      else if (f_1 .eq. f_b) then
         sg_1 = sg_b
         do i = 1, nvar
            x(i) = x(i) + cube*s(i)
         end do
      else if (f_1 .eq. f_c) then
         sg_1 = sg_c
      end if
c
c     try to restart from best point with smaller stepsize
c
      if (f_1 .gt. f_0) then
         ncalls = ncalls + 1
         f = fgvalue (x,g)
         status = 'IntplnErr'
         return
      end if
      f_0 = f_1
      sg_0 = sg_1
      if (sg_1 .gt. 0.0d0) then
         do i = 1, nvar
            s(i) = -s(i)
         end do
         sg_0 = -sg_1
      end if
      step = max(cube,step-cube) / 10.0d0
      if (step .lt. stpmin)  step = stpmin
c
c     if already restarted once, return with best point
c
      if (status .eq. ' ReSearch') then
         ncalls = ncalls + 1
         f = fgvalue (x,g)
         status = 'BadIntpln'
         return
      else
         status = ' ReSearch'
         goto 10
      end if
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
c     ##  subroutine setime  --  initialize elapsed CPU time clock  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setime" initializes the elapsed interval CPU timer
c
c
      subroutine setime
      implicit none
      include 'chrono.i'
c
c
c     initialize interval at elapsed CPU time for current job
c
      call clock (cputim)
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine shakeup  --  setup holonomic constraints  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "shakeup" initializes any holonomic constraints for use
c     with the rattle algorithm during molecular dynamics
c
c
      subroutine shakeup
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'couple.i'
      include 'keys.i'
      include 'math.i'
      include 'molcul.i'
      include 'shake.i'
      integer i,j,k,next
      integer ia,ib,ic,ita,itb,itc
      real*8 rab,rbc,rac,cosine
      character*9 rattyp
      character*20 keyword
      character*80 record,string
c
c
c     zero out the number of rattle distance constraints
c
      nrat = 0
      use_rattle = .true.
c
c     search each line of the keyword file for options
c
      do k = 1, nkey
         next = 1
         record = keyline(k)
         call upcase (record)
         call gettext (record,keyword,next)
c
c     get the distance constraint types for the rattle method
c
         if (keyword(1:7) .eq. 'RATTLE ') then
            call getword (record,rattyp,next)
c
c     constrain all bond lengths to their ideal values
c
            if (rattyp(1:5) .eq. 'BONDS') then
               do i = 1, nbond
                  ia = ibnd(1,i)
                  ib = ibnd(2,i)
                  nrat = nrat + 1
                  irat(1,nrat) = ia
                  irat(2,nrat) = ib
                  krat(nrat) = bl(i)
               end do
c
c     constrain 1-3 distance across each angle at ideal value
c
            else if (rattyp(1:6) .eq. 'ANGLES') then
               do i = 1, nangle
                  ia = iang(1,i)
                  ib = iang(2,i)
                  ic = iang(3,i)
                  do j = 1, n12(ib)
                     if (i12(j,ib) .eq. ia) then
                        rab = bl(bndlist(j,ib))
                     else if (i12(j,ib) .eq. ic) then
                        rbc = bl(bndlist(j,ib))
                     end if
                  end do
                  cosine = cos(anat(i)/radian)
                  rac = sqrt(rab**2+rbc**2-2.0d0*rab*rbc*cosine)
                  nrat = nrat + 1
                  irat(1,nrat) = ia
                  irat(2,nrat) = ic
                  krat(nrat) = rac
               end do
c
c     fix bond length in diatomics to give a rigid molecule
c
            else if (rattyp(1:8) .eq. 'DIATOMIC') then
               do i = 1, nbond
                  ia = ibnd(1,i)
                  ib = ibnd(2,i)
                  if (n12(ia).eq.1 .and. n12(ib).eq.1) then
                     nrat = nrat + 1
                     irat(1,nrat) = ia
                     irat(2,nrat) = ib
                     krat(nrat) = bl(i)
                  end if
               end do
c
c     fix bonds and angle in triatomics to give a rigid molecule
c
            else if (rattyp(1:9) .eq. 'TRIATOMIC') then
               do i = 1, nangle
                  ia = iang(1,i)
                  ib = iang(2,i)
                  ic = iang(3,i)
                  if (n12(ia)+n12(ib)+n12(ic) .eq. 4) then
                     rab = bl(bndlist(1,ia))
                     rbc = bl(bndlist(1,ic))
                     cosine = cos(anat(i)/radian)
                     rac = sqrt(rab**2+rbc**2-2.0d0*rab*rbc*cosine)
                     nrat = nrat + 1
                     irat(1,nrat) = ia
                     irat(2,nrat) = ib
                     krat(nrat) = rab
                     nrat = nrat + 1
                     irat(1,nrat) = ib
                     irat(2,nrat) = ic
                     krat(nrat) = rbc
                     nrat = nrat + 1
                     irat(1,nrat) = ia
                     irat(2,nrat) = ic
                     krat(nrat) = rac
                  end if
               end do
c
c     fix bonds and angle of each water to give a rigid molecule
c
            else if (rattyp(1:5) .eq. 'WATER') then
               do i = 1, nangle
                  ia = iang(1,i)
                  ib = iang(2,i)
                  ic = iang(3,i)
                  ita = atomic(ia)
                  itb = atomic(ib)
                  itc = atomic(ic)
                  if (ita.eq.1 .and. itb.eq.8 .and. itc.eq.1) then
                     rab = bl(bndlist(1,ia))
                     rbc = bl(bndlist(1,ic))
                     cosine = cos(anat(i)/radian)
                     rac = sqrt(rab**2+rbc**2-2.0d0*rab*rbc*cosine)
                     nrat = nrat + 1
                     irat(1,nrat) = ia
                     irat(2,nrat) = ib
                     krat(nrat) = rab
                     nrat = nrat + 1
                     irat(1,nrat) = ib
                     irat(2,nrat) = ic
                     krat(nrat) = rbc
                     nrat = nrat + 1
                     irat(1,nrat) = ia
                     irat(2,nrat) = ic
                     krat(nrat) = rac
                  end if
               end do
c
c     fix all bonds to hydrogen atoms at their ideal length
c
            else
               do i = 1, nbond
                  ia = ibnd(1,i)
                  ib = ibnd(2,i)
                  ita = atomic(ia)
                  itb = atomic(ib)
                  if (ita.eq.1 .or. itb.eq.1) then
                     nrat = nrat + 1
                     irat(1,nrat) = ia
                     irat(2,nrat) = ib
                     krat(nrat) = bl(i)
                  end if
               end do
            end if
c
c     fix a single specified atom pair at the desired distance
c
         else if (keyword(1:12) .eq. 'RATTLE-BOND ') then
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            rab = 0.0d0
            string = record(next:80)
            read (string,*,err=10,end=10)  rab
   10       continue
            if (rab .eq. 0.0d0) then
               do i = 1, n12(ia)
                  if (i12(i,ia) .eq. ib) then
                     rab = bl(bndlist(i,ia))
                  end if
               end do
            end if
            if (rab .eq. 0.0d0) then
               rab = sqrt((x(ia)-x(ib))**2 + (y(ia)-y(ib))**2
     &                           + (z(ia)-z(ib))**2)
            end if
            nrat = nrat + 1
            irat(1,nrat) = ia
            irat(2,nrat) = ib
            krat(nrat) = rab
         end if
      end do
c
c     set flag to apply minimum image to intermolecular rattles
c
      do i = 1, nrat
         ia = irat(1,nrat)
         ib = irat(2,nrat)
         if (use_bounds .and. (molcule(ia).ne.molcule(ib))) then
            ratimage(i) = .true.
         else
            ratimage(i) = .false.
         end if
      end do
c
c     if no rattle constraints are used, turn off its use
c
      if (nrat .eq. 0)  use_rattle = .false.
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
c     ##  function sigmoid  --  general sigmoidal functional form  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "sigmoid" implements a normalized sigmoidal function on the
c     interval [0,1]; the curves connect (0,0) to (1,1) and have
c     a cooperativity controlled by beta, they approach a straight
c     line as beta -> 0 and get more nonlinear as beta increases
c
c
      function sigmoid (beta,x)
      implicit none
      real*8 beta,x,sigmoid
      real*8 expmax,expmin,expterm
c
c
c     compute the value of the normalized sigmoidal function
c
      if (beta .eq. 0.0d0) then
         sigmoid = x
      else
         expmax = 1.0d0 / (exp(-beta) + 1.0d0)
         expmin = 1.0d0 / (exp(beta) + 1.0d0)
         expterm = 1.0d0 / (exp(beta*(2.0d0*x-1.0d0)) + 1.0d0)
         sigmoid = (expmax - expterm) / (expmax - expmin)
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine smooth  --  set potential smoothing parameter  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "smooth" sets extent of potential surface deformation for use
c     with potential smoothing plus search, the diffusion equation
c     method or Gaussian density annealing
c
c
      subroutine smooth
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'fields.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'warp.i'
      integer i,next
      character*20 keyword
      character*80 record,string
      logical query,exist
c
c
c     set the default value for the deformation parameter
c
      deform = 0.0d0
c
c     set default values for the diffusion coefficients
c
      diffb = 0.000156d0
      diffa = 0.0014d0
      diffid = 0.0225d0
      difft = 0.0225d0
      diffv = 1.0d0
      diffc = 1.0d0
c
c     set a flag to determine use of smoothed potentials
c
      if (forcefield .eq. 'SMOOTH') then
         use_deform = .true.
         query = .true.
      else
         use_deform = .false.
         query = .false.
      end if
c
c     get keyword containing the deformation parameter value
c
      if (query) then
         do i = 1, nkey
            next = 1
            record = keyline(i)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:7) .eq. 'DEFORM ') then
               query = .false.
               string = record(next:80)
               read (string,*,err=10)  deform
            end if
   10       continue
         end do
      end if
c
c     try to get the deformation value from the command line
c
      if (query) then
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=20,end=20)  deform
            query = .false.
         end if
   20    continue
      end if
c
c     ask for the potential surface deformation to be used
c
      if (query) then
         if (use_gda) then
            deform = 200.0d0
            write (iout,30)
   30       format (/,' Enter the Initial Mean Squared Gaussian',
     &                 ' Width [200.0] :  ',$)
         else
            write (iout,40)
   40       format (/,' Enter the Potential Surface Smoothing',
     &                 ' Parameter [0.0] :  ',$)
         end if
         read (input,50)  record
   50    format (a80)
         read (record,*,err=60)  deform
   60    continue
      end if
c
c     get keywords containing diffusion coefficients
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:13) .eq. 'DIFFUSE-BOND ') then
            string = record(next:80)
            read (string,*,err=70)  diffb
         else if (keyword(1:14) .eq. 'DIFFUSE-ANGLE ') then
            string = record(next:80)
            read (string,*,err=70)  diffa
         else if (keyword(1:17) .eq. 'DIFFUSE-IMPROPER ') then
            string = record(next:80)
            read (string,*,err=70)  diffid
         else if (keyword(1:16) .eq. 'DIFFUSE-TORSION ') then
            string = record(next:80)
            read (string,*,err=70)  difft
         else if (keyword(1:12) .eq. 'DIFFUSE-VDW ') then
            string = record(next:80)
            read (string,*,err=70)  diffv
         else if (keyword(1:15) .eq. 'DIFFUSE-CHARGE ') then
            string = record(next:80)
            read (string,*,err=70)  diffc
         end if
   70    continue
      end do
c
c     set second moment of Gaussian on each atom for GDA methods
c
      if (use_gda) then
         do i = 1, n
            m2(i) = deform
         end do
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine solvate  --  macroscopic solvation parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "solvate" assigns macroscopic solvation energy parameters
c     for the Eisenberg-McLachlan ASP, Ooi-Scheraga SASA or 
c     Macromodel GB/SA solvation models
c
c     literature references:
c
c     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
c     Applied to Molecular Dynamics of Proteins in Solution",
c     Protein Science, 1, 227-235 (1992)  (Eisenberg-McLachlan)
c
c     T. Ooi, M. Oobatake, G. Nemethy and H. A. Scheraga, "Accessible
c     Surface Areas as a Measure of the Thermodynamic Parameters of
c     Hydration of Peptides", PNAS, 84, 3086-3090 (1987)  (SASA)
c
c     W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson,
c     "A Semianalytical Treatment of Solvation for Molecular
c     Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129
c     (1990)  (Macromodel GB/SA)
c
c
      subroutine solvate
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'keys.i'
      include 'potent.i'
      include 'solute.i'
      integer i,j,k,next
      integer atmnum,brn
      character*4 solvtyp
      character*8 value
      character*20 keyword
      character*80 record,string
c
c
c     default is to not use the macroscopic solvation term
c
      use_solv = .false.
      use_gbsa = .false.
      nsolv = 0
c
c     set default values for the GB/SA solvation model
c
      reborn = 1
      bornmax = 300
c
c     search keywords for macroscopic solvation commands
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'SOLVATE ') then
            use_solv = .true.
            call gettext (record,value,next)
            call upcase (value)
            if (value(1:5) .eq. 'GBSA ') then
               use_gbsa = .true.
               solvtyp = 'GBSA'
               string = record(next:80)
               read (string,*,err=10,end=10)  brn
               bornmax = brn
   10          continue
            else if (value(1:5) .eq. 'SASA ') then
               solvtyp = 'SASA'
            else
               solvtyp = 'E-M'
            end if
         end if
      end do
c
c     assign the Eisenberg-McLachlan ASP solvation parameters;
c     parameters only available for protein-peptide groups
c
      if (use_solv .and. solvtyp.eq.'E-M') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 6) then
               nsolv = nsolv + 1
               rsolv(i) = 1.9d0
               vsolv(i) = 0.004d0
            else if (atmnum .eq. 7) then
               nsolv = nsolv + 1
               rsolv(i) = 1.7d0
               vsolv(i) = -0.113d0
               if (n12(i) .eq. 4) then
                  vsolv(i) = -0.169d0
               end if
            else if (atmnum .eq. 8) then
               nsolv = nsolv + 1
               rsolv(i) = 1.4d0
               vsolv(i) = -0.113d0
               if (n12(i).eq.1 .and. atomic(i12(1,i)).eq.6) then
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (n12(k).eq.1 .and. atomic(k).eq.8) then
                        vsolv(i) = -0.166d0
                     end if
                  end do
               end if
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 15)  vsolv(i) = -0.140d0
               end do
            else if (atmnum .eq. 15) then
               nsolv = nsolv + 1
               rsolv(i) = 1.9d0
               vsolv(i) = -0.140d0
            else if (atmnum .eq. 16) then
               nsolv = nsolv + 1
               rsolv(i) = 1.8d0
               vsolv(i) = -0.017d0
            else
               rsolv(i) = 0.0d0
               vsolv(i) = 0.0d0
            end if
         end do
c
c     assign the Ooi-Scheraga SASA solvation parameters;
c     parameters only available for protein-peptide groups
c
      else if (use_solv .and. solvtyp.eq.'SASA') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 6) then
               nsolv = nsolv + 1
               rsolv(i) = 2.0d0
               vsolv(i) = 0.008d0
               if (n12(i) .eq. 3) then
                  rsolv(i) = 1.75d0
                  vsolv(i) = -0.008d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 8) then
                        rsolv(i) = 1.55d0
                        vsolv(i) = 0.427d0
                     end if
                  end do
               end if
            else if (atmnum .eq. 7) then
               nsolv = nsolv + 1
               rsolv(i) = 1.55d0
               vsolv(i) = -0.132d0
               if (n12(i) .eq. 4)  vsolv(i) = -1.212d0
            else if (atmnum .eq. 8) then
               nsolv = nsolv + 1
               rsolv(i) = 1.4d0
               if (n12(i) .eq. 1) then
                  vsolv(i) = -0.038d0
                  if (atomic(i12(1,i)) .eq. 6) then
                     do j = 1, n13(i)
                        k = i13(j,i)
                        if (n12(k).eq.1 .and. atomic(k).eq.8) then
                           vsolv(i) = -0.770d0
                        end if
                     end do
                  end if
               else if (n12(i) .eq. 2) then
                  vsolv(i) = -0.172d0
               end if
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 15)  vsolv(i) = -0.717d0
               end do
            else if (atmnum .eq. 15) then
               nsolv = nsolv + 1
               rsolv(i) = 2.1d0
               vsolv(i) = 0.0d0
            else if (atmnum .eq. 16) then
               nsolv = nsolv + 1
               rsolv(i) = 2.0d0
               vsolv(i) = -0.021d0
            else if (atmnum .eq. 17) then
               nsolv = nsolv + 1
               rsolv(i) = 2.0d0
               vsolv(i) = 0.012d0
            else
               rsolv(i) = 0.0d0
               vsolv(i) = 0.0d0
            end if
         end do
c
c     assign the GB/SA parameters for the surface area term
c
      else if (use_solv .and. solvtyp.eq.'GBSA') then
         do i = 1, n
            atmnum = atomic(i)
            vsolv(i) = 0.0072d0
            if (atmnum .eq. 1) then
               nsolv = nsolv + 1
               rsolv(i) = 1.21d0
               k = i12(1,i)
               if (atomic(k) .eq. 7)  rsolv(i) = 1.15d0
               if (atomic(k) .eq. 8)  rsolv(i) = 1.15d0
            else if (atmnum .eq. 3) then
               nsolv = nsolv + 1
               rsolv(i) = 1.432d0
            else if (atmnum .eq. 6) then
               nsolv = nsolv + 1
               rsolv(i) = 1.90d0
            else if (atmnum .eq. 7) then
               nsolv = nsolv + 1
               rsolv(i) = 1.625d0
            else if (atmnum .eq. 8) then
               nsolv = nsolv + 1
               rsolv(i) = 1.535d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.48d0
            else if (atmnum .eq. 9) then
               nsolv = nsolv + 1
               rsolv(i) = 1.47d0
            else if (atmnum .eq. 11) then
               nsolv = nsolv + 1
               rsolv(i) = 1.992d0
            else if (atmnum .eq. 12) then
               nsolv = nsolv + 1
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 14) then
               nsolv = nsolv + 1
               rsolv(i) = 1.80d0
            else if (atmnum .eq. 15) then
               nsolv = nsolv + 1
               rsolv(i) = 1.87d0
            else if (atmnum .eq. 16) then
               nsolv = nsolv + 1
               rsolv(i) = 1.775d0
            else if (atmnum .eq. 17) then
               nsolv = nsolv + 1
               rsolv(i) = 1.735d0
            else if (atmnum .eq. 19) then
               nsolv = nsolv + 1
               rsolv(i) = 2.123d0
            else if (atmnum .eq. 20) then
               nsolv = nsolv + 1
               rsolv(i) = 1.817d0
            else if (atmnum .eq. 35) then
               nsolv = nsolv + 1
               rsolv(i) = 1.90d0
            else if (atmnum .eq. 37) then
               nsolv = nsolv + 1
               rsolv(i) = 2.26d0
            else if (atmnum .eq. 53) then
               nsolv = nsolv + 1
               rsolv(i) = 2.10d0
            else if (atmnum .eq. 55) then
               nsolv = nsolv + 1
               rsolv(i) = 2.507d0
            else if (atmnum .eq. 56) then
               nsolv = nsolv + 1
               rsolv(i) = 2.188d0
            else
               rsolv(i) = 0.0d0
            end if
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
c     ##  subroutine solve  --  linear equation solve via CG method  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "solve" uses the linear conjugate gradient method to find
c     an (approximate) solution to the set of linear equations
c     represented in matrix form by Hp = -g (Newton's equations)
c
c     status codes upon return:
c
c     TruncNewt    convergence to (truncated) Newton criterion
c     NegCurve     termination upon detecting negative curvature
c     OverLimit    maximum number of CG iterations exceeded
c
c
      subroutine solve (mode,method,negtest,nvar,p,x,g,h,
     &                  h_init,h_stop,h_index,h_diag,cycle,
     &                  iter_cg,fg_call,fgvalue,status,
     &                  f,c_index,c_value,maxhess)
      implicit none
      include 'sizes.i'
      include 'output.i'
casa
      include 'iounit.i'
      integer i,j,k,nvar,cycle,iter,iter_cg,fg_call,maxiter,maxhess
      integer h_init(maxvar),h_stop(maxvar),h_index(maxhess)
      real*8 x(maxvar),g(maxvar),p(maxvar),q(maxvar)
      real*8 m(maxvar),r(maxvar),s(maxvar),d(maxvar)
      real*8 h_diag(maxvar),h(maxhess)
      real*8 alpha,beta,delta,sigma
      real*8 f_sigma,x_sigma(maxvar),g_sigma(maxvar)
      real*8 fgvalue,g_norm,g_rms,eps,converge
      real*8 hj,gg,dq,rr,dd,rs,rs_new,r_norm
      real*8 f(maxhess)
      integer c_index(maxhess),c_value(maxhess)
      character*6 mode,method
      character*9 status
      logical negtest
      external fgvalue
c
c
c     transformation using exact Hessian diagonal
c
      if (mode.ne.'dtncg ' .and. method.ne.'none  ') then
         do i = 1, nvar
            m(i) = 1.0d0 / sqrt(abs(h_diag(i)))
         end do
         do i = 1, nvar
            g(i) = g(i) * m(i)
            h_diag(i) = h_diag(i) * m(i) * m(i)
            do j = h_init(i), h_stop(i)
               k = h_index(j)
               h(j) = h(j) * m(i) * m(k)
            end do
         end do
      end if
c
c     setup prior to linear conjugate gradient iterations
c
      iter = 0
      gg = 0.0d0
      do i = 1, nvar
         p(i) = 0.0d0
         r(i) = -g(i)
         gg = gg + g(i)*g(i)
      end do
      g_norm = sqrt(gg)
      call precond (method,iter,nvar,s,r,h,h_init,
     &                 h_stop,h_index,h_diag, f,c_index,c_value,maxhess)
      rs = 0.0d0
      do i = 1, nvar
         d(i) = s(i)
         rs = rs + r(i)*s(i)
      end do
      if (mode .eq. 'newton') then
         eps = 1.0d-10
         maxiter = nvar
      else if (mode.eq.'tncg  ' .or. mode.eq.'dtncg ') then
         delta = 1.0d0
         eps = delta / dble(cycle)
         g_rms = g_norm / sqrt(dble(nvar))
         eps = min(eps,g_rms)
         converge = 1.0d0
         eps = eps**converge
         maxiter = nint(10.0d0*sqrt(dble(nvar)))
      end if
      iter = 1
c
c     evaluate or estimate the matrix-vector product
c
      dowhile (.true.)
         if (mode.eq.'tncg  ' .or. mode.eq.'newton') then
            do i = 1, nvar
               q(i) = 0.0d0
            end do
            do i = 1, nvar
               q(i) = q(i) + h_diag(i)*d(i)
               do j = h_init(i), h_stop(i)
                  k = h_index(j)
                  hj = h(j)
                  q(i) = q(i) + hj*d(k)
                  q(k) = q(k) + hj*d(i)
               end do
            end do
         else if (mode .eq. 'dtncg ') then
            dd = 0.0d0
            do i = 1, nvar
               dd = dd + d(i)*d(i)
            end do
            sigma = 1.0d-7 / sqrt(dd)
            if (coordtype .eq. 'internal ') then
               sigma = 1.0d-4 / sqrt(dd)
            end if
            do i = 1, nvar
               x_sigma(i) = x(i) + sigma*d(i)
            end do
            fg_call = fg_call + 1
            f_sigma = fgvalue (x_sigma,g_sigma)
            do i = 1, nvar
               q(i) = (g_sigma(i)-g(i)) / sigma
            end do
         end if
c
c     check for a direction of negative curvature
c
         dq = 0.0d0
         do i = 1, nvar
            dq = dq + d(i)*q(i)
         end do
         if (negtest) then
            if (dq .le. 0.0d0) then
               if (iter .eq. 1) then
                  do i = 1, nvar
                     p(i) = d(i)
                  end do
               end if
               status = ' NegCurve'
               goto 10
            end if
         end if
c
c     test the truncated Newton termination criterion
c
         alpha = rs / dq
         rr = 0.0d0
         do i = 1, nvar
            p(i) = p(i) + alpha*d(i)
            r(i) = r(i) - alpha*q(i)
            rr = rr + r(i)*r(i)
         end do
         r_norm = sqrt(rr)
         if (r_norm/g_norm .le. eps) then
            status = 'TruncNewt'
            goto 10
         end if
c
c     solve the preconditioning equations
c
         call precond (method,iter,nvar,s,r,h,h_init,
     &                 h_stop,h_index,h_diag,f,c_index,c_value,maxhess)
c
c     update the truncated Newton direction
c
         rs_new = 0.0d0
         do i = 1, nvar
            rs_new = rs_new + r(i)*s(i)
         end do
         beta = rs_new / rs
         rs = rs_new
         do i = 1, nvar
            d(i) = s(i) + beta*d(i)
         end do
c
c     check for overlimit, then begin next iteration
c
         if (iter .ge. maxiter) then
            status = 'OverLimit'
            goto 10
         end if
         iter = iter + 1
      end do
c
c     retransform and increment total iterations, then terminate
c
   10 continue
      if (mode.ne.'dtncg ' .and. method.ne.'none  ') then
         do i = 1, nvar
            p(i) = p(i) * m(i)
            g(i) = g(i) / m(i)
         end do
      end if
      iter_cg = iter_cg + iter
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
c     ##  subroutine sort  --  heapsort of an integer array  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "sort" takes an input list of integers and sorts it
c     into ascending order using the Heapsort algorithm
c
c
      subroutine sort (n,list)
      implicit none
      integer i,j,k,n,index
      integer list(*),lists
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sort2  --  heapsort of real array with keys  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sort2" takes an input list of reals and sorts it
c     into ascending order using the Heapsort algorithm;
c     it also returns a key into the original ordering
c
c
      subroutine sort2 (n,list,key)
      implicit none
      integer i,j,k,n,index
      integer key(*),keys
      real*8 list(*),lists
c
c
c     initialize index into the original ordering
c
      do i = 1, n
         key(i) = i
      end do
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
            keys = key(k)
         else
            lists = list(index)
            keys = key(index)
            list(index) = list(1)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
         key(i) = keys
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine sort3  --  heapsort of integer array with keys  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "sort3" takes an input list of integers and sorts it
c     into ascending order using the Heapsort algorithm;
c     it also returns a key into the original ordering
c
c
      subroutine sort3 (n,list,key)
      implicit none
      integer i,j,k,n,index
      integer list(*),lists
      integer key(*),keys
c
c
c     initialize index into the original ordering
c
      do i = 1, n
         key(i) = i
      end do
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
            keys = key(k)
         else
            lists = list(index)
            keys = key(index)
            list(index) = list(1)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
         key(i) = keys
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine sort4  --  heapsort of integer absolute values  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "sort4" takes an input list of integers and sorts it into
c     ascending absolute value using the Heapsort algorithm
c
c
      subroutine sort4 (n,list)
      implicit none
      integer i,j,k,n,index
      integer list(*),lists
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. index)
            if (j .lt. index) then
               if (abs(list(j)) .lt. abs(list(j+1)))  j = j + 1
            end if
            if (abs(lists) .lt. abs(list(j))) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine sort5  --  heapsort of integer array modulo m  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "sort5" takes an input list of integers and sorts it
c     into ascending order based on each value modulo "m"
c
c
      subroutine sort5 (n,list,m)
      implicit none
      integer i,j,k,m,n,index
      integer jmod,j1mod,smod
      integer list(*),lists
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. index)
            if (j .lt. index) then
               jmod = mod(list(j),m)
               j1mod = mod(list(j+1),m)
               if (jmod .lt. j1mod) then
                  j = j + 1
               else if (jmod.eq.j1mod .and. list(j).lt.list(j+1)) then
                  j = j + 1
               end if
            end if
            smod = mod(lists,m)
            jmod = mod(list(j),m)
            if (smod .lt. jmod) then
               list(i) = list(j)
               i = j
               j = j + j
            else if (smod.eq.jmod .and. lists.lt.list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine sort6  --  heapsort of a text string array  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sort6" takes an input list of character strings and sorts
c     it into alphabetical order using the Heapsort algorithm
c
c
      subroutine sort6 (n,list)
      implicit none
      integer i,j,k,n,index
      character*256 lists
      character*(*) list(*)
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine sort7  --  heapsort of text strings with keys  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "sort7" takes an input list of character strings and sorts it
c     into alphabetical order using the Heapsort algorithm; it also
c     returns a key into the original ordering
c
c
      subroutine sort7 (n,list,key)
      implicit none
      integer i,j,k,n,index
      integer key(*),keys
      character*256 lists
      character*(*) list(*)
c
c
c     initialize index into the original ordering
c
      do i = 1, n
         key(i) = i
      end do
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
            keys = key(k)
         else
            lists = list(index)
            keys = key(index)
            list(index) = list(1)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
         key(i) = keys
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
c     ##  subroutine square  --  nonlinear least squares with bounds  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "square" is a nonlinear least squares routine derived from
c     the IMSL routine BCLSF and More's Minpack routine LMDER; the
c     Jacobian is estimated by finite differences and bounds can
c     be specified for the variables to be refined
c
c     arguments and variables :
c
c     n        number of variables to optimize
c     m        number of residual functions
c     xlo      vector of length n containing the lower bounds for
c                the variables
c     xhi      vector of length n containing the upper bounds for
c                the variables
c     xscale   vector of length n containing the diagonal scaling
c                matrix for the variables
c     xc       vector of length n containing the variable values
c                at the approximate solution
c     fc       vector of length m containing the residuals at the
c                approximate solution
c     fp       vector of length m containing the updated residual
c     xp       vector of length n containing the updated point
c     sc       vector of length n containing the last step taken
c     gc       vector of length n containing an estimate of
c                the gradient at the approximate solution
c     fjac     real m by n matrix containing an estimate of
c                the Jacobian at the approximate solution
c     mdim     leading dimension of fjac exactly as specified in
c                the dimension statement of the calling program
c     iactive  integer vector of length n indicating if xc(i) had
c                to be moved to an upper or lower bound
c     ipvt     vector of length n containing the permutation matrix
c                used in the QR factorization of the Jacobian at
c                the approximate solution
c     stpmax   real scalar containing the maximum allowed step size
c     delta    real scalar containing the trust region radius
c
c     required external routines :
c
c     rsdvalue   subroutine to evaluate residual function values
c     writeout   subroutine to write out info about current status
c
c
      subroutine square (m,n,xlo,xhi,xc,fc,gc,fjac,mdim,
     &                      grdmin,rsdvalue,writeout)
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'keys.i'
      include 'minima.i'
      integer maxlsq,maxrsd
      parameter (maxlsq=50)
      parameter (maxrsd=100)
      integer i,j,k,icode,next
      integer niter,ncalls,nactive
      integer m,n,mdim,nbigstp,ndigit
      integer iactive(maxlsq),ipvt(maxlsq)
      real*8 xc(maxlsq),xlo(maxlsq),xhi(maxlsq)
      real*8 fc(maxrsd),gc(maxlsq),fjac(mdim,maxlsq)
      real*8 xp(maxlsq),fp(maxrsd),ga(maxlsq),gs(maxlsq)
      real*8 sc(maxlsq),sa(maxlsq),ftemp(maxrsd)
      real*8 xscale(maxlsq),xsa(maxlsq)
      real*8 rdiag(maxlsq),qtf(maxlsq),work(maxlsq)
      real*8 amu,delta,epsfcn,fcnorm,fpnorm,gcnorm,ganorm
      real*8 precise,eps,grdmin,stpnorm,stpmax,stpmin
      real*8 rftol,faketol,xtemp,stepsz,sum,temp
      character*20 keyword
      character*80 record
      logical done,first,gauss,bigstp
      external rsdvalue,writeout,precise
c
c
c     check for too many variables or residuals
c
      if (n .gt. maxlsq) then
         write (iout,10)
   10    format (/,' SQUARE  --  Too many Parameters,',
     &              ' Increase value of MAXLSQ')
         return
      else if (m .gt. maxrsd) then
         write (iout,20)
   20    format (/,' SQUARE  --  Too many Residuals,',
     &              ' Increase value of MAXRSD')
         return
      end if
c
c     initialize various counters and status code
c
      niter = 0
      ncalls = 0
      nbigstp = 0
      done = .false.
c
c     setup the default tolerances and parameter values
c
      eps = precise (2)
      ndigit = 10
      if (eps .lt. 10.0d0**(-ndigit)) then
         eps = 10.0d0**(-ndigit)
      end if
      if (maxiter .eq. 0)  maxiter = 100
      if (iprint .lt. 0)  iprint = 1
      if (iwrite .lt. 0)  iwrite = 1
      if (fctmin .eq. 0.0d0)  fctmin = eps
      if (grdmin .eq. 0.0d0)  grdmin = eps**(1.0d0/3.0d0)
      epsfcn = sqrt(eps)
      delta = 0.0d0
      stpmax = 1000.0d0 * sqrt(dble(n))
      stpmin = eps**(2.0d0/3.0d0)
      rftol = eps**(2.0d0/3.0d0)
      faketol = 100.0d0 * eps
c
c     search each line of the keyword file for options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'FCTMIN ') then
            read (record(next:80),*,err=30)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (record(next:80),*,err=30)  maxiter
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (record(next:80),*,err=30)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (record(next:80),*,err=30)  iwrite
         end if
   30    continue
      end do
c
c     check feasibility of variables and use bounds if needed
c
      nactive = 0
      do j = 1, n
         if (xc(j) .lt. xlo(j)) then
            xc(j) = xlo(j)
            iactive(j) = -1
         else if (xc(j) .gt. xhi(j)) then
            xc(j) = xhi(j)
            iactive(j) = 1
         else
            nactive = nactive + 1
            iactive(j) = 0
         end if
      end do
c
c     evaluate the function at the initial point
c
      ncalls = ncalls + 1
      call rsdvalue (m,n,xc,fc)
      fcnorm = 0.0d0
      do i = 1, m
         fcnorm = fcnorm + fc(i)**2
      end do
      fcnorm = 0.5d0 * fcnorm
c
c     evaluate the Jacobian at the initial point by finite
c     differences; replace loop with user routine if desired
c
      do j = 1, n
         stepsz = epsfcn * abs(xc(j))
         if (stepsz .lt. epsfcn)  stepsz = epsfcn
         if (xc(j) .lt. 0.0d0)  stepsz = -stepsz
         xtemp = xc(j)
         xc(j) = xtemp + stepsz
         ncalls = ncalls + 1
         call rsdvalue (m,n,xc,ftemp)
         xc(j) = xtemp
         do i = 1, m
            fjac(i,j) = (ftemp(i)-fc(i)) / stepsz
         end do
      end do
c
c     compute More's adaptive variable scale factors
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + fjac(i,j)**2
         end do
         xscale(j) = sqrt(temp)
         if (xscale(j) .eq. 0.0d0)  xscale(j) = 1.0d0
      end do
c
c     compute the total gradient vector for all variables
c
      do j = 1, n
         gc(j) = 0.0d0
         do i = 1, m
            gc(j) = gc(j) + fjac(i,j)*fc(i)
         end do
      end do
c
c     compute the norm of the scaled total gradient
c     and the scaled gradient for active variables
c
      gcnorm = 0.0d0
      ganorm = 0.0d0
      do j = 1, n
         gs(j) = gc(j) * max(abs(xc(j)),1.0d0/xscale(j))
         gcnorm = gcnorm + gs(j)**2
         if (iactive(j) .eq. 0) then
            ganorm = ganorm + gs(j)**2
         end if
      end do
      gcnorm = sqrt(gcnorm/n)
      if (nactive .ne. 0)  ganorm = sqrt(ganorm/nactive)
c
c     print out information about initial conditions
c
      if (iprint .gt. 0) then
         write (iout,40)
   40    format (/,' Levenberg-Marquardt Nonlinear Least Squares :')
         write (iout,50)
   50    format (/,' LS Iter    F Value      Total G     Active G',
     &             '    N Active   F Calls',/)
         if (max(fcnorm,gcnorm) .lt. 10000000.0d0) then
            write (iout,60)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
   60       format (i6,3f13.4,2i10)
         else
            write (iout,70)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
   70       format (i6,3d13.4,2i10)
         end if
      end if
c
c     write out the parameters, derivatives and residuals
c
      if (iwrite .eq. 1) then
         if (.not. done)  call writeout (niter,xc,gs,m,fc)
      end if
c
c     check stopping criteria at the initial point; test the
c     absolute function value and gradient norm for termination
c
      if (fcnorm .le. fctmin)  return
      if (ganorm .le. grdmin)  return
c
c     start of the main body of least squares iteration
c
   80 continue
      niter = niter + 1
c
c     repack the Jacobian to include only active variables
c
      if (nactive .ne. n) then
         k = 0
         do j = 1, n
            if (iactive(j) .ne. 0) then
               if (k .eq. 0)  k = j
            else
               if (k .ne. 0) then
                  do i = 1, m
                     fjac(i,k) = fjac(i,j)
                  end do
                  k = k + 1
               end if
            end if
         end do
      end if
c
c     repack scale factors and gradient for active variables
c
      k = 0
      do j = 1, n
         if (iactive(j) .eq. 0) then
            k = k + 1
            xsa(k) = xscale(j)
            ga(k) = gc(j)
         end if
      end do
c
c     compute the QR factorization of the Jacobian
c
      call qrfact (m,nactive,fjac,mdim,.true.,ipvt,rdiag,work)
c
c     compute the vector Q(transpose) * residuals
c
      do i = 1, m
         qtf(i) = fc(i)
      end do
      do j = 1, nactive
         if (fjac(j,j) .ne. 0.0d0) then
            sum = 0.0d0
            do i = j, m
               sum = sum + fjac(i,j)*qtf(i)
            end do
            temp = -sum / fjac(j,j)
            do i = j, m
               qtf(i) = qtf(i) + fjac(i,j)*temp
            end do
         end if
         fjac(j,j) = rdiag(j)
      end do
c
c     compute the Levenberg-Marquardt step
c
      icode = 6
      first = .true.
      dowhile (icode .ge. 4)
         call lmstep (nactive,ga,fjac,mdim,ipvt,xsa,qtf,stpmax,
     &                       delta,amu,first,sa,gauss)
c
c     unpack the step vector to include all variables
c
         k = 0
         do i = 1, n
            if (iactive(i) .ne. 0) then
               sc(i) = 0.0d0
            else
               k = k + 1
               sc(i) = sa(k)
            end if
         end do
c
c     check new point and update the trust region
c
         call trust (rsdvalue,m,n,xc,fcnorm,gc,fjac,mdim,ipvt,sc,sa,
     &               xscale,gauss,stpmax,delta,icode,xp,fc,fp,fpnorm,
     &               bigstp,ncalls,xlo,xhi,nactive,stpmin,rftol,faketol)
      end do
      if (icode .eq. 1)  done = .true.
c
c     update to the new variables and residuals
c
      do j = 1, n
         xc(j) = xp(j)
      end do
      do i = 1, m
         fc(i) = fp(i)
      end do
      fcnorm = fpnorm
c
c     update the active vs inactive status of the variables;
c     in a true active set strategy, at most one constraint
c     is added to the active set per iteration
c
      do j = 1, n
         if (iactive(j) .eq. 0) then
            if (abs(xc(j)-xlo(j)) .le. eps) then
               nactive = nactive - 1
               iactive(j) = -1
c              goto 110
            else if (abs(xc(j)-xhi(j)) .le. eps) then
               nactive = nactive - 1
               iactive(j) = 1
c              goto 110
            end if
         end if
      end do
   90 continue
c
c     evaluate the Jacobian at the new point using finite
c     differences; replace loop with user routine if desired
c
      do j = 1, n
         stepsz = epsfcn * max(abs(xc(j)),1.0d0/xscale(j))
         if (xc(j) .lt. 0.0d0)  stepsz = -stepsz
         xtemp = xc(j)
         xc(j) = xtemp + stepsz
         ncalls = ncalls + 1
         call rsdvalue (m,n,xc,ftemp)
         xc(j) = xtemp
         do i = 1, m
            fjac(i,j) = (ftemp(i)-fc(i)) / stepsz
         end do
      end do
c
c     compute More's adaptive variable scale factors
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + fjac(i,j)**2
         end do
         xscale(j) = max(xscale(j),sqrt(temp))
      end do
c
c     compute the total gradient vector for all variables
c
      do j = 1, n
         gc(j) = 0.0d0
         do i = 1, m
            gc(j) = gc(j) + fjac(i,j)*fc(i)
         end do
      end do
c
c     compute the norm of the scaled total gradient
c     and the scaled gradient for active variables
c
      gcnorm = 0.0d0
      ganorm = 0.0d0
      do j = 1, n
         gs(j) = gc(j) * max(abs(xc(j)),1.0d0/xscale(j))
         gcnorm = gcnorm + gs(j)**2
         if (iactive(j) .eq. 0) then
            ganorm = ganorm + gs(j)**2
         end if
      end do
      gcnorm = sqrt(gcnorm/n)
      if (nactive .ne. 0)  ganorm = sqrt(ganorm/nactive)
c
c     print out information about current iteration
c
      if (iprint.ne.0 .and. mod(niter,iprint).eq.0) then
         if (max(fcnorm,gcnorm) .lt. 10000000.0d0) then
            write (iout,100)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
  100       format (i6,3f13.4,2i10)
         else
            write (iout,110)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
  110       format (i6,3d13.4,2i10)
         end if
      end if
c
c     check stopping criteria at the new point; test the absolute
c     function value, the gradient norm and step for termination
c
      if (fcnorm .le. fctmin)  done = .true.
      if (ganorm .le. grdmin)  done = .true.
      stpnorm = 0.0d0
      do j = 1, n
         temp = max(abs(xc(j)),1.0d0/xscale(j))
         stpnorm = stpnorm + (sc(j)/temp)**2
      end do
      stpnorm = sqrt(stpnorm/n)
      if (stpnorm .le. stpmin)  done = .true.
c
c     check for inactive variables that can be made active;
c     in a true active set strategy, variables are released
c     one at a time at a minimum of the current active set
c
c     if (done) then
      if (nactive .ne. n) then
         do j = 1, n
            if (iactive(j).eq.-1 .and. gc(j).lt.0.0d0) then
               nactive = nactive + 1
               iactive(j) = 0
               done = .false.
c                 goto 140
            else if (iactive(j).eq.1 .and. gc(j).gt.0.0d0) then
               nactive = nactive + 1
               iactive(j) = 0
               done = .false.
c                 goto 140
            end if
         end do
  120    continue
      end if
c     end if
c
c     if still done, then normal termination has been achieved
c
      if (done) then
         write (iout,130)
  130    format (/,' SQUARE  --  Normal Termination of Least Squares')
      end if
c
c     check the limit on the number of iterations
c
      if (niter .ge. maxiter) then
         done = .true.
         write (iout,140)
  140    format (/,' SQUARE  --  Maximum Number of Allowed Iterations')
      end if
c
c     check relative function convergence and false convergence
c
      if (icode .eq. 2) then
         done = .true.
         write (iout,150)
  150    format (/,' SQUARE  --  Relative Function Convergence',
     &           //,' Both the scaled actual and predicted',
     &               ' reductions in the function',
     &           /,' are less than or equal to the relative',
     &               ' convergence tolerance')
      else if (icode .eq. 3) then
         done = .true.
         write (iout,160)
  160    format (/,' SQUARE  --  Possible False Convergence',
     &           //,' The iterates appear to be converging to',
     &               ' a noncritical point due',
     &           /,' to bad gradient information, discontinuous',
     &               ' function, or stopping',
     &           /,' tolerances being too tight')
      end if
c
c     check for several consecutive maximum steps taken
c
      if (bigstp) then
         nbigstp = nbigstp + 1
         if (nbigstp .eq. 5) then
            done = .true.
            write (iout,170)
  170       format (/,' SQUARE  --  Five Consecutive Maximum',
     &                  ' Length Steps',
     &              //,' Either the function is unbounded below,',
     &                  ' or has a finite',
     &              /,' asymptote in some direction, or STEPMAX',
     &                  ' is too small')
         end if
      else
         nbigstp = 0
      end if
c
c     write out the parameters, derivatives and residuals
c
      if (iwrite.ne.0 .and. mod(niter,iwrite).eq.0) then
         if (.not. done)  call writeout (niter,xc,gs,m,fc)
      end if
c
c     continue with the next iteration if not finished
c
      if (.not. done)  goto 80
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
c     ##  subroutine suffix  --  test for default file extension  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "suffix" checks a filename for the presence of an
c     extension, and appends an extension if none is found
c
c
      subroutine suffix (filename,extension)
      implicit none
      integer i,leng,lext
      integer last,trimtext
      character*1 letter
      character*(*) filename,extension
      logical exist
c
c
c     get the length of the current filename
c
      leng = trimtext (filename)
      lext = trimtext (extension)
c
c     check for an extension on the current filename
c
      last = leng
      do i = 1, leng
         letter = filename(i:i)
         if (letter .eq. '/')  last = leng
c        if (letter .eq. '\')  last = leng
         if (ichar(letter) .eq. 92)  last = leng
         if (letter .eq. ']')  last = leng
         if (letter .eq. ':')  last = leng
         if (letter .eq. '~')  last = leng
         if (letter .eq. '.')  last = i - 1
      end do
      if (last .ne. leng)  return
c
c     append extension if current name does not exist
c
      exist = .false.
      if (leng .ne. 0)  inquire (file=filename(1:leng),exist=exist)
      if (.not. exist) then
         filename = filename(1:leng)//'.'//extension(1:lext)
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
c     ##  subroutine surface  --  accessible surface area & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "surface" performs an analytical computation of the weighted
c     solvent accessible surface area of each atom and the first
c     derivatives of the area with respect to Cartesian coordinates
c
c     literature references:
c
c     T. J. Richmond, "Solvent Accessible Surface Area and
c     Excluded Volume in Proteins", Journal of Molecular Biology,
c     178, 63-89 (1984)
c
c     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
c     Applied to Molecular Dynamics of Proteins in Solution",
c     Protein Science, 1, 227-235 (1992)
c
c     variables and parameters:
c
c     total    total surface area of the whole structure
c     area     accessible surface area of each atom
c     darea    x,y,z components of the gradient of the area of
c                the molecule with respect to atomic coordinates
c     radius   radii of the individual atoms
c     weight   weight assigned to each atom's area; if set to
c                1.0, return is actual area in square Angstroms
c     probe    radius of the probe sphere
c     sig      tolerance used in the tests for sphere overlaps
c                and for colinearity; if a connectivity error
c                occurs for an atom, change the coordinates of
c                the offending atom slightly
c
c
      subroutine surface (total,area,darea,radius,weight,probe)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'usage.i'
      integer maxarc
      parameter (maxarc=300)
      integer i,j,k,l,m,ii,ib,jb,in,io,ir
      integer mi,ni,narc,ider(maxarc),sign_yder(maxarc)
      integer key(maxarc),intag(maxarc),intag1(maxarc)
      real*8 total,area(maxatm),darea(3,maxatm)
      real*8 probe,radius(maxatm),r(maxatm)
      real*8 wght,weight(maxatm)
      real*8 sig,sigsq,arcsum,dsql,wxl,wxlsq,p,s,v,rcn,cosine
      real*8 axx,axy,axz,ayx,ayy,azx,azy,azz,uxl,uyl,uzl
      real*8 tx,ty,tz,txb,tyb,t2,td,tr2,tr,txr,tyr,tk1,tk2
      real*8 thec,the,t,tb,txk,tyk,tzk,t1,ti,tf,tt,txl,tyl,tzl
      real*8 arclen,exang,xr,yr,zr,rr,rrx2,rrsq,rplus,rminus
      real*8 ccsq,cc,xysq,bgl,bsqk,bsql,bk,gi,gl,pix2,pix4,pid2
      real*8 dax,day,daz,deal,decl,dtkal,dtkcl,dtlal,dtlcl
      real*8 therk,dk,gk,risqk,rik,risql,faca,facb,facc,gaca,gacb
      real*8 ux(maxarc),uy(maxarc),uz(maxarc)
      real*8 xc(maxarc),yc(maxarc),zc(maxarc)
      real*8 xc1(maxarc),yc1(maxarc),zc1(maxarc)
      real*8 dsq(maxarc),bsq(maxarc),b(maxarc)
      real*8 arci(maxarc),arcf(maxarc),ex(maxarc)
      real*8 kent(maxarc),kout(maxarc),lt(maxarc)
      real*8 bg(maxarc),ther(maxarc),ri(maxarc),risq(maxarc)
      real*8 b1(maxarc),dsq1(maxarc),bsq1(maxarc),gr(maxarc)
      logical ltop,skip(maxatm),omit(maxarc),komit
c
c
c     set multiples of pi and the overlap significance value
c
      pix2 = 2.0d0 * pi
      pix4 = 4.0d0 * pi
      pid2 = pi / 2.0d0
      sig = 0.01d0
      sigsq = sig**2
      total = 0.0d0
      do i = 1, maxarc
         ider(i) = 0
         sign_yder(i) = 0
      end do
c
c     zero the area and derivatives, and set the sphere radii
c
      do i = 1, n
         area(i) = 0.0d0
         darea(1,i) = 0.0d0
         darea(2,i) = 0.0d0
         darea(3,i) = 0.0d0
         r(i) = radius(i)
         if (r(i) .ne. 0.0d0)  r(i) = r(i) + probe
      end do
c
c     set the "skip" array to exclude all inactive atoms
c     that do not overlap any of the current active atoms
c
      do i = 1, n
         skip(i) = .true.
      end do
      do i = 1, n
         if (use(i)) then
            xr = x(i)
            yr = y(i)
            zr = z(i)
            rr = r(i)
            do k = 1, n
               rplus = (rr + r(k))**2
               ccsq = (x(k)-xr)**2 + (y(k)-yr)**2 + (z(k)-zr)**2
               if (ccsq .le. rplus)  skip(k) = .false.
            end do
         end if
      end do
c
c     compute the area and derivatives of current "ir" sphere
c
      do ir = 1, n
         if (skip(ir))  goto 160
         io = 0
         jb = 0
         ib = 0
         arclen = 0.0d0
         exang = 0.0d0
         xr = x(ir)
         yr = y(ir)
         zr = z(ir)
         rr = r(ir)
         rrx2 = 2.0d0 * rr
         rrsq = rr**2
         wght = weight(ir)
c
c     test each sphere to see if it overlaps the "ir" sphere
c
         do i = 1, n
            if (i .eq. ir)  goto 20
            rplus = rr + r(i)
            tx = x(i) - xr
            if (abs(tx) .ge. rplus)  goto 20
            ty = y(i) - yr
            if (abs(ty) .ge. rplus)  goto 20
            tz = z(i) - zr
            if (abs(tz) .ge. rplus)  goto 20
c
c     check for overlap of spheres by testing center to
c     center distance against sum and difference of radii
c
            xysq = tx**2 + ty**2
            if (xysq .lt. sigsq) then
               tx = sig
               ty = 0.0d0
               xysq = sigsq
            end if
            ccsq = xysq + tz**2
            cc = sqrt(ccsq)
            if (rplus-cc .le. sig)  goto 20
            rminus = rr - r(i)
c
c     check for a completely buried "ir" sphere
c
            if (cc-abs(rminus) .le. sig) then
               if (rminus .le. 0.0d0)  goto 160
               goto 20
            end if
c
c     calculate overlap parameters between "i" and "ir" sphere
c
            io = io + 1
            xc1(io) = tx
            yc1(io) = ty
            zc1(io) = tz
            dsq1(io) = xysq
            bsq1(io) = ccsq
            b1(io) = cc
            gr(io) = (ccsq+rplus*rminus) / (rrx2*b1(io))
            intag1(io) = i
            if (io .gt. maxarc) then
               write (iout,10)
   10          format (/,' SURFACE  --  Increase the Value of MAXARC')
               call fatal
            end if
   20       continue
         end do
c
c     case where no other spheres overlap the current sphere
c
         if (io .eq. 0) then
            area(ir) = pix4
            goto 150
         end if
c
c     case where only one sphere overlaps the current sphere
c
         if (io .eq. 1) then
            k = 1
            txk = xc1(1)
            tyk = yc1(1)
            tzk = zc1(1)
            bsqk = bsq1(1)
            bk = b1(1)
            intag(1) = intag1(1)
            arcsum = pix2
            ib = ib + 1
            arclen = arclen + gr(k)*arcsum
            in = intag(k)
            t1 = arcsum*rrsq*(bsqk-rrsq+r(in)**2) / (rrx2*bsqk*bk)
            darea(1,ir) = darea(1,ir) - txk*t1*wght
            darea(2,ir) = darea(2,ir) - tyk*t1*wght
            darea(3,ir) = darea(3,ir) - tzk*t1*wght
            darea(1,in) = darea(1,in) + txk*t1*wght
            darea(2,in) = darea(2,in) + tyk*t1*wght
            darea(3,in) = darea(3,in) + tzk*t1*wght
            goto 140
         end if
c
c     general case where more than one sphere intersects the
c     current sphere; sort intersecting spheres by their degree
c     of overlap with the current main sphere
c
         call sort2 (io,gr,key)
         do i = 1, io
            k = key(i)
            intag(i) = intag1(k)
            xc(i) = xc1(k)
            yc(i) = yc1(k)
            zc(i) = zc1(k)
            dsq(i) = dsq1(k)
            b(i) = b1(k)
            bsq(i) = bsq1(k)
            omit(i) = .false.
         end do
c
c     radius of the each circle on the surface of the "ir" sphere
c
         do i = 1, io
            gi = gr(i) * rr
            bg(i) = b(i) * gi
            risq(i) = rrsq - gi**2
            ri(i) = sqrt(risq(i))
            ther(i) = pid2 - asin(min(1.0d0,max(-1.0d0,gr(i))))
         end do
c
c     find boundary of inaccessible area on "ir" sphere
c
         do k = 1, io-1
            if (.not. omit(k)) then
               txk = xc(k)
               tyk = yc(k)
               tzk = zc(k)
               bk = b(k)
               therk = ther(k)
               do j = k+1, io
                  if (omit(j))  goto 50
c
c     check to see if J circle is intersecting K circle;
c     get distance between circle centers and sum of radii
c
                  cc = (txk*xc(j)+tyk*yc(j)+tzk*zc(j))/(bk*b(j))
                  cc = acos(min(1.0d0,max(-1.0d0,cc)))
                  td = therk + ther(j)
c
c     check to see if circles enclose separate regions
c
                  if (cc .ge. td)  goto 50
c
c     check for circle J completely inside circle K
c
                  if (cc+ther(j) .lt. therk)  goto 30
c
c     check for circles essentially parallel
c
                  if (cc .gt. sig)  goto 40
   30             continue
                  omit(j) = .true.
                  goto 50
c
c     check for "ir" sphere completely buried
c
   40             continue
                  if (pix2-cc .le. td)  goto 160
   50             continue
               end do
            end if
         end do
c
c     find T value of circle intersections
c
         do k = 1, io
            if (omit(k))  goto 100
            komit = omit(k)
            omit(k) = .true.
            narc = 0
            ltop = .false.
            txk = xc(k)
            tyk = yc(k)
            tzk = zc(k)
            dk = sqrt(dsq(k))
            bsqk = bsq(k)
            bk = b(k)
            gk = gr(k) * rr
            risqk = risq(k)
            rik = ri(k)
            therk = ther(k)
c
c     rotation matrix elements
c
            t1 = tzk / (bk*dk)
            axx = txk * t1
            axy = tyk * t1
            axz = dk / bk
            ayx = tyk / dk
            ayy = txk / dk
            azx = txk / bk
            azy = tyk / bk
            azz = tzk / bk
            do l = 1, io
               if (.not. omit(l)) then
                  txl = xc(l)
                  tyl = yc(l)
                  tzl = zc(l)
c
c     rotate spheres so K vector colinear with z-axis
c
                  uxl = txl*axx + tyl*axy - tzl*axz
                  uyl = tyl*ayy - txl*ayx
                  uzl = txl*azx + tyl*azy + tzl*azz
                  cosine = min(1.0d0,max(-1.0d0,uzl/b(l)))
                  if (acos(cosine) .lt. therk+ther(l)) then
                     dsql = uxl**2 + uyl**2
                     tb = uzl*gk - bg(l)
                     txb = uxl * tb
                     tyb = uyl * tb
                     td = rik * dsql
                     tr2 = risqk*dsql - tb**2
                     if (tr2 .le. 0.000001d0)  tr2 = 0.000001d0
                     tr = sqrt(tr2)
                     txr = uxl * tr
                     tyr = uyl * tr
c
c     get T values of intersection for K circle
c
                     tb = (txb+tyr) / td
                     tb = min(1.0d0,max(-1.0d0,tb))
                     tk1 = acos(tb)
                     if (tyb-txr .lt. 0.0d0)  tk1 = pix2 - tk1
                     tb = (txb-tyr) / td
                     tb = min(1.0d0,max(-1.0d0,tb))
                     tk2 = acos(tb)
                     if (tyb+txr .lt. 0.0d0)  tk2 = pix2 - tk2
                     thec = (rrsq*uzl-gk*bg(l)) / (rik*ri(l)*b(l))
                     if (abs(thec) .lt. 1.0d0) then
                        the = -acos(thec)
                     else if (thec .ge. 1.0d0) then
                        the = 0.0d0
                     else if (thec .le. -1.0d0) then
                        the = -pi
                     end if
c
c     see if "tk1" is entry or exit point; check t=0 point;
c     "ti" is exit point, "tf" is entry point
c
                     cosine = min(1.0d0,max(-1.0d0,
     &                               (uzl*gk-uxl*rik)/(b(l)*rr)))
                     if ((acos(cosine)-ther(l))*(tk2-tk1)
     &                          .le. 0.0d0) then
                        ti = tk2
                        tf = tk1
                     else
                        ti = tk2
                        tf = tk1
                     end if
                     narc = narc + 1
                     if (narc .ge. maxarc) then
                        write (iout,60)
   60                   format (/,' SURFACE  --  Increase the Value',
     &                             ' of MAXARC')
                        call fatal
                     end if
                     if (tf .le. ti) then
                        arcf(narc) = tf
                        arci(narc) = 0.0d0
                        tf = pix2
                        lt(narc) = l
                        ex(narc) = the
                        ltop = .true.
                        narc = narc + 1
                     end if
                     arcf(narc) = tf
                     arci(narc) = ti
                     lt(narc) = l
                     ex(narc) = the
                     ux(l) = uxl
                     uy(l) = uyl
                     uz(l) = uzl
                  end if
               end if
            end do
            omit(k) = komit
c
c     special case; K circle without intersections
c
            if (narc .le. 0)  goto 80
c
c     general case; sum up arclength and set connectivity code
c
            call sort2 (narc,arci,key)
            arcsum = arci(1)
            mi = key(1)
            t = arcf(mi)
            ni = mi
            if (narc .gt. 1) then
               do j = 2, narc
                  m = key(j)
                  if (t .lt. arci(j)) then
                     arcsum = arcsum + arci(j) - t
                     exang = exang + ex(ni)
                     jb = jb + 1
                     if (jb .ge. maxarc) then
                        write (iout,70)
   70                   format (/,' SURFACE  --  Increase the Value',
     &                             ' of MAXARC')
                        call fatal
                     end if
                     l = lt(ni)
                     ider(l) = ider(l) + 1
                     sign_yder(l) = sign_yder(l) + 1
                     kent(jb) = maxarc*l + k
                     l = lt(m)
                     ider(l) = ider(l) + 1
                     sign_yder(l) = sign_yder(l) - 1
                     kout(jb) = maxarc*k + l
                  end if
                  tt = arcf(m)
                  if (tt .ge. t) then
                     t = tt
                     ni = m
                  end if
               end do
            end if
            arcsum = arcsum + pix2 - t
            if (.not. ltop) then
               exang = exang + ex(ni)
               jb = jb + 1
               l = lt(ni)
               ider(l) = ider(l) + 1
               sign_yder(l) = sign_yder(l) + 1
               kent(jb) = maxarc*l + k
               l = lt(mi)
               ider(l) = ider(l) + 1
               sign_yder(l) = sign_yder(l) - 1
               kout(jb) = maxarc*k + l
            end if
c
c     calculate the surface area derivatives
c
            do l = 1, io
               if (ider(l) .ne. 0) then
                  rcn = ider(l) * rrsq
                  ider(l) = 0
                  uzl = uz(l)
                  gl = gr(l) * rr
                  bgl = bg(l)
                  bsql = bsq(l)
                  risql = risq(l)
                  wxlsq = bsql - uzl**2
                  wxl = sqrt(wxlsq)
                  p = bgl - gk*uzl
                  v = risqk*wxlsq - p**2
                  if (v .le. 0.000001d0)  v = 0.000001d0
                  v = sqrt(v)
                  t1 = rr * (gk*(bgl-bsql)+uzl*(bgl-rrsq))
     &                             / (v*risql*bsql)
                  deal = -wxl*t1
                  decl = -uzl*t1 - rr/v
                  dtkal = (wxlsq-p) / (wxl*v)
                  dtkcl = (uzl-gk) / v
                  s = gk*b(l) - gl*uzl
                  t1 = 2.0d0*gk - uzl
                  t2 = rrsq - bgl
                  dtlal = -(risql*wxlsq*b(l)*t1-s*(wxlsq*t2+risql*bsql))
     &                              / (risql*wxl*bsql*v)
                  dtlcl = -(risql*b(l)*(uzl*t1-bgl)-uzl*t2*s)
     &                              / (risql*bsql*v)
                  gaca = rcn * (deal-(gk*dtkal-gl*dtlal)/rr) / wxl
                  gacb = (gk - uzl*gl/b(l)) * sign_yder(l) * rr / wxlsq
                  sign_yder(l) = 0
                  faca = ux(l)*gaca - uy(l)*gacb
                  facb = uy(l)*gaca + ux(l)*gacb
                  facc = rcn * (decl-(gk*dtkcl-gl*dtlcl)/rr)
                  dax = axx*faca - ayx*facb + azx*facc
                  day = axy*faca + ayy*facb + azy*facc
                  daz = azz*facc - axz*faca
                  in = intag(l)
                  darea(1,ir) = darea(1,ir) + dax*wght
                  darea(2,ir) = darea(2,ir) + day*wght
                  darea(3,ir) = darea(3,ir) + daz*wght
                  darea(1,in) = darea(1,in) - dax*wght
                  darea(2,in) = darea(2,in) - day*wght
                  darea(3,in) = darea(3,in) - daz*wght
               end if
            end do
            goto 90
   80       continue
            arcsum = pix2
            ib = ib + 1
   90       continue
            arclen = arclen + gr(k)*arcsum
            in = intag(k)
            t1 = arcsum*rrsq*(bsqk-rrsq+r(in)**2) / (rrx2*bsqk*bk)
            darea(1,ir) = darea(1,ir) - txk*t1*wght
            darea(2,ir) = darea(2,ir) - tyk*t1*wght
            darea(3,ir) = darea(3,ir) - tzk*t1*wght
            darea(1,in) = darea(1,in) + txk*t1*wght
            darea(2,in) = darea(2,in) + tyk*t1*wght
            darea(3,in) = darea(3,in) + tzk*t1*wght
  100       continue
         end do
         if (arclen .eq. 0.0d0)  goto 160
         if (jb .eq. 0)  goto 140
c
c     find number of independent boundaries; check
c     for connectivity error for the current atom
c
         j = 0
         do k = 1, jb
            if (kout(k) .ne. 0) then
               i = k
  110          continue
               m = kout(i)
               kout(i) = 0
               j = j + 1
               do ii = 1, jb
                  if (m .eq. kent(ii)) then
                     if (ii .eq. k) then
                        ib = ib + 1
                        if (j .eq. jb)  goto 140
                        goto 120
                     end if
                     i = ii
                     goto 110
                  end if
               end do
  120          continue
            end if
         end do
         ib = ib + 1
         write (iout,130)  ir
  130    format (/,' SURFACE  --  Connectivity Error at Atom',i6)
c
c     form the accessible area for the current atom
c
  140    continue
         area(ir) = ib*pix2 + exang + arclen
         area(ir) = mod(area(ir),pix4)
  150    continue
         area(ir) = area(ir) * rrsq * wght
         total = total + area(ir)
  160    continue
      end do
c
c     zero out the area derivatives for the inactive atoms
c
      do i = 1, n
         if (.not. use(i)) then
            darea(1,i) = 0.0d0
            darea(2,i) = 0.0d0
            darea(3,i) = 0.0d0
         end if
      end do
c
c     print out the surface area and derivatives for each atom
c
      if (debug) then
         write (iout,170)
  170    format (/,' Weighted Atomic Surface Areas and Derivatives :',
     &           //,4x,'Atom',7x,'Area Term',10x,'dA/dx',
     &              7x,'dA/dy',7x,'dA/dz',/)
         do i = 1, n
            if (.not. skip(i)) then
               write (iout,180)  i,area(i),(darea(j,i),j=1,3)
  180          format (i8,4x,f12.4,3x,3f12.4)
            end if
         end do
         write (iout,190)  total
  190    format (/,' Total Weighted Surface Area :',5x,f16.4)
      end if
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
c     ##  subroutine surfatom  --  exposed surface area of an atom  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "surfatom" performs an analytical computation of the surface
c     area of a specified atom; a simplified version of "surface"
c
c     literature references:
c
c     T. J. Richmond, "Solvent Accessible Surface Area and
c     Excluded Volume in Proteins", Journal of Molecular Biology,
c     178, 63-89 (1984)
c
c     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
c     Applied to Molecular Dynamics of Proteins in Solution",
c     Protein Science, 1, 227-235 (1992)
c
c     variables and parameters:
c
c     ir       number of atom for which area is desired
c     area     accessible surface area of the atom
c     radius   radii of each of the individual atoms
c
c
      subroutine surfatom (ir,area,r)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'math.i'
      integer maxarc
      parameter (maxarc=300)
      integer i,j,k,m,ii,ib,jb,io,ir,mi,ni,narc
      integer ider(maxarc),sign_yder(maxarc)
      integer key(maxarc),intag(maxarc),intag1(maxarc)
      real*8 area,delta,delta2,arcsum,arclen,exang
      real*8 xr,yr,zr,rr,rrsq,rplus,rminus
      real*8 axx,axy,axz,ayx,ayy
      real*8 azx,azy,azz,uxj,uyj,uzj
      real*8 tx,ty,tz,txb,tyb,td
      real*8 tr2,tr,txr,tyr,tk1,tk2
      real*8 thec,the,t,tb,txk,tyk,tzk
      real*8 t1,ti,tf,tt,txj,tyj,tzj
      real*8 ccsq,cc,xysq,bsqk,bk,cosine
      real*8 dsqj,gi,pix2,therk,dk,gk,risqk,rik
      real*8 r(maxatm),ri(maxarc),risq(maxarc)
      real*8 ux(maxarc),uy(maxarc),uz(maxarc)
      real*8 xc(maxarc),yc(maxarc),zc(maxarc)
      real*8 xc1(maxarc),yc1(maxarc),zc1(maxarc)
      real*8 dsq(maxarc),bsq(maxarc)
      real*8 dsq1(maxarc),bsq1(maxarc)
      real*8 arci(maxarc),arcf(maxarc)
      real*8 ex(maxarc),lt(maxarc),gr(maxarc)
      real*8 b(maxarc),b1(maxarc),bg(maxarc)
      real*8 kent(maxarc),kout(maxarc),ther(maxarc)
      logical top,omit(maxarc)
c
c
c     zero out the surface area for the sphere of interest
c
      area = 0.0d0
      if (r(ir) .eq. 0.0d0)  return
c
c     set the overlap significance and zero out counters
c
      pix2 = 2.0d0 * pi
      delta = 0.01d0
      delta2 = delta * delta
      io = 0
      jb = 0
      ib = 0
      arclen = 0.0d0
      exang = 0.0d0
      do i = 1, maxarc
         ider(i) = 0
         sign_yder(i) = 0
      end do
c
c     store coordinates and radius of the sphere of interest
c
      xr = x(ir)
      yr = y(ir)
      zr = z(ir)
      rr = r(ir)
      rrsq = rr * rr
c
c     test each sphere to see if it overlaps the sphere of interest
c
      do i = 1, n
         if (i.eq.ir .or. r(i).eq.0.0d0)  goto 20
         rplus = rr + r(i)
         tx = x(i) - xr
         if (abs(tx) .ge. rplus)  goto 20
         ty = y(i) - yr
         if (abs(ty) .ge. rplus)  goto 20
         tz = z(i) - zr
         if (abs(tz) .ge. rplus)  goto 20
c
c     check for sphere overlap by testing distance against radii
c
         xysq = tx*tx + ty*ty
         if (xysq .lt. delta2) then
            tx = delta
            ty = 0.0d0
            xysq = delta2
         end if
         ccsq = xysq + tz*tz
         cc = sqrt(ccsq)
         if (rplus-cc .le. delta)  goto 20
         rminus = rr - r(i)
c
c     check to see if sphere of interest is completely buried
c
         if (cc-abs(rminus) .le. delta) then
            if (rminus .le. 0.0d0)  goto 150
            goto 20
         end if
c
c     check for too many overlaps with sphere of interest
c
         if (io .ge. maxarc) then
            write (iout,10)
   10       format (/,' SURFATOM  --  Increase the Value of MAXARC')
            call fatal
         end if
c
c     get overlap between current sphere and sphere of interest
c
         io = io + 1
         xc1(io) = tx
         yc1(io) = ty
         zc1(io) = tz
         dsq1(io) = xysq
         bsq1(io) = ccsq
         b1(io) = cc
         gr(io) = (ccsq+rplus*rminus) / (2.0d0*rr*b1(io))
         intag1(io) = i
         omit(io) = .false.
   20    continue
      end do
c
c     case where no other spheres overlap the sphere of interest
c
      if (io .eq. 0) then
         area = 4.0d0 * pi * rrsq
         return
      end if
c
c     case where only one sphere overlaps the sphere of interest
c
      if (io .eq. 1) then
         area = pix2 * (1.0d0 + gr(1))
         area = mod(area,4.0d0*pi) * rrsq
         return
      end if
c
c     case where many spheres intersect the sphere of interest;
c     sort the intersecting spheres by their degree of overlap
c
      call sort2 (io,gr,key)
      do i = 1, io
         k = key(i)
         intag(i) = intag1(k)
         xc(i) = xc1(k)
         yc(i) = yc1(k)
         zc(i) = zc1(k)
         dsq(i) = dsq1(k)
         b(i) = b1(k)
         bsq(i) = bsq1(k)
      end do
c
c     get radius of each overlap circle on surface of the sphere
c
      do i = 1, io
         gi = gr(i) * rr
         bg(i) = b(i) * gi
         risq(i) = rrsq - gi*gi
         ri(i) = sqrt(risq(i))
         ther(i) = 0.5d0*pi - asin(min(1.0d0,max(-1.0d0,gr(i))))
      end do
c
c     find boundary of inaccessible area on sphere of interest
c
      do k = 1, io-1
         if (.not. omit(k)) then
            txk = xc(k)
            tyk = yc(k)
            tzk = zc(k)
            bk = b(k)
            therk = ther(k)
c
c     check to see if J circle is intersecting K circle;
c     get distance between circle centers and sum of radii
c
            do j = k+1, io
               if (omit(j))  goto 50
               cc = (txk*xc(j)+tyk*yc(j)+tzk*zc(j))/(bk*b(j))
               cc = acos(min(1.0d0,max(-1.0d0,cc)))
               td = therk + ther(j)
c
c     check to see if circles enclose separate regions
c
               if (cc .ge. td)  goto 50
c
c     check for circle J completely inside circle K
c
               if (cc+ther(j) .lt. therk)  goto 30
c
c     check for circles that are essentially parallel
c
               if (cc .gt. delta)  goto 40
   30          continue
               omit(j) = .true.
               goto 50
c
c     check to see if sphere of interest is completely buried
c
   40          continue
               if (pix2-cc .le. td)  goto 150
   50          continue
            end do
         end if
      end do
c
c     find T value of circle intersections
c
      do k = 1, io
         if (omit(k))  goto 100
         omit(k) = .true.
         narc = 0
         top = .false.
         txk = xc(k)
         tyk = yc(k)
         tzk = zc(k)
         dk = sqrt(dsq(k))
         bsqk = bsq(k)
         bk = b(k)
         gk = gr(k) * rr
         risqk = risq(k)
         rik = ri(k)
         therk = ther(k)
c
c     rotation matrix elements
c
         t1 = tzk / (bk*dk)
         axx = txk * t1
         axy = tyk * t1
         axz = dk / bk
         ayx = tyk / dk
         ayy = txk / dk
         azx = txk / bk
         azy = tyk / bk
         azz = tzk / bk
         do j = 1, io
            if (.not. omit(j)) then
               txj = xc(j)
               tyj = yc(j)
               tzj = zc(j)
c
c     rotate spheres so K vector colinear with z-axis
c
               uxj = txj*axx + tyj*axy - tzj*axz
               uyj = tyj*ayy - txj*ayx
               uzj = txj*azx + tyj*azy + tzj*azz
               cosine = min(1.0d0,max(-1.0d0,uzj/b(j)))
               if (acos(cosine) .lt. therk+ther(j)) then
                  dsqj = uxj*uxj + uyj*uyj
                  tb = uzj*gk - bg(j)
                  txb = uxj * tb
                  tyb = uyj * tb
                  td = rik * dsqj
                  tr2 = risqk*dsqj - tb*tb
                  if (tr2 .le. 0.000001d0)  tr2 = 0.000001d0
                  tr = sqrt(tr2)
                  txr = uxj * tr
                  tyr = uyj * tr
c
c     get T values of intersection for K circle
c
                  tb = (txb+tyr) / td
                  tb = min(1.0d0,max(-1.0d0,tb))
                  tk1 = acos(tb)
                  if (tyb-txr .lt. 0.0d0)  tk1 = pix2 - tk1
                  tb = (txb-tyr) / td
                  tb = min(1.0d0,max(-1.0d0,tb))
                  tk2 = acos(tb)
                  if (tyb+txr .lt. 0.0d0)  tk2 = pix2 - tk2
                  thec = (rrsq*uzj-gk*bg(j)) / (rik*ri(j)*b(j))
                  if (abs(thec) .lt. 1.0d0) then
                     the = -acos(thec)
                  else if (thec .ge. 1.0d0) then
                     the = 0.0d0
                  else if (thec .le. -1.0d0) then
                     the = -pi
                  end if
c
c     see if "tk1" is entry or exit point; check t=0 point;
c     "ti" is exit point, "tf" is entry point
c
                  cosine = min(1.0d0,max(-1.0d0,
     &                            (uzj*gk-uxj*rik)/(b(j)*rr)))
                  if ((acos(cosine)-ther(j))*(tk2-tk1) .le. 0.0d0) then
                     ti = tk2
                     tf = tk1
                  else
                     ti = tk2
                     tf = tk1
                  end if
                  narc = narc + 1
                  if (narc .ge. maxarc) then
                     write (iout,60)
   60                format (/,' SURFATOM  --  Increase the Value',
     &                          ' of MAXARC')
                     call fatal
                  end if
                  if (tf .le. ti) then
                     arcf(narc) = tf
                     arci(narc) = 0.0d0
                     tf = pix2
                     lt(narc) = j
                     ex(narc) = the
                     top = .true.
                     narc = narc + 1
                  end if
                  arcf(narc) = tf
                  arci(narc) = ti
                  lt(narc) = j
                  ex(narc) = the
                  ux(j) = uxj
                  uy(j) = uyj
                  uz(j) = uzj
               end if
            end if
         end do
         omit(k) = .false.
c
c     special case; K circle without intersections
c
         if (narc .le. 0)  goto 80
c
c     general case; sum up arclength and set connectivity code
c
         call sort2 (narc,arci,key)
         arcsum = arci(1)
         mi = key(1)
         t = arcf(mi)
         ni = mi
         if (narc .gt. 1) then
            do j = 2, narc
               m = key(j)
               if (t .lt. arci(j)) then
                  arcsum = arcsum + arci(j) - t
                  exang = exang + ex(ni)
                  jb = jb + 1
                  if (jb .ge. maxarc) then
                     write (iout,70)
   70                format (/,' SURFATOM  --  Increase the Value',
     &                          ' of MAXARC')
                     call fatal
                  end if
                  i = lt(ni)
                  ider(i) = ider(i) + 1
                  sign_yder(i) = sign_yder(i) + 1
                  kent(jb) = maxarc*i + k
                  i = lt(m)
                  ider(i) = ider(i) + 1
                  sign_yder(i) = sign_yder(i) - 1
                  kout(jb) = maxarc*k + i
               end if
               tt = arcf(m)
               if (tt .ge. t) then
                  t = tt
                  ni = m
               end if
            end do
         end if
         arcsum = arcsum + pix2 - t
         if (.not. top) then
            exang = exang + ex(ni)
            jb = jb + 1
            i = lt(ni)
            ider(i) = ider(i) + 1
            sign_yder(i) = sign_yder(i) + 1
            kent(jb) = maxarc*i + k
            i = lt(mi)
            ider(i) = ider(i) + 1
            sign_yder(i) = sign_yder(i) - 1
            kout(jb) = maxarc*k + i
         end if
         goto 90
   80    continue
         arcsum = pix2
         ib = ib + 1
   90    continue
         arclen = arclen + gr(k)*arcsum
  100    continue
      end do
      if (arclen .eq. 0.0d0)  goto 150
      if (jb .eq. 0)  goto 140
c
c     find number of independent boundaries
c     and check for connectivity error
c
      j = 0
      do k = 1, jb
         if (kout(k) .ne. 0) then
            i = k
  110       continue
            m = kout(i)
            kout(i) = 0
            j = j + 1
            do ii = 1, jb
               if (m .eq. kent(ii)) then
                  if (ii .eq. k) then
                     ib = ib + 1
                     if (j .eq. jb)  goto 140
                     goto 120
                  end if
                  i = ii
                  goto 110
               end if
            end do
  120       continue
         end if
      end do
      ib = ib + 1
      write (iout,130)  ir
  130 format (/,' SURFATOM  --  Connectivity Error at Atom',i6)
c
c     compute the exposed surface area for the sphere of interest
c
  140 continue
      area = ib*pix2 + exang + arclen
      area = mod(area,4.0d0*pi) * rrsq
  150 continue
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
c     ##  subroutine switch  --  get switching function coefficients  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "switch" sets the coeffcients used by the fifth and seventh
c     order polynomial switching functions for spherical cutoffs
c
c
      subroutine switch (mode)
      implicit none
      include 'cutoff.i'
      include 'shunt.i'
      real*8 denom,term
      real*8 off3,off4,off5,off6,off7
      real*8 cut3,cut4,cut5,cut6,cut7
      character*6 mode
c
c
c     set switching coefficients to zero for truncation cutoffs
c
      c0 = 0.0d0
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
      c4 = 0.0d0
      c5 = 0.0d0
      f0 = 0.0d0
      f1 = 0.0d0
      f2 = 0.0d0
      f3 = 0.0d0
      f4 = 0.0d0
      f5 = 0.0d0
      f6 = 0.0d0
      f7 = 0.0d0
c
c     get the switching window for the current potential type
c
      if (mode(1:3) .eq. 'VDW') then
         off = vdwcut
         cut = vdwtaper
      else if (mode .eq. 'CHARGE') then
         off = chgcut
         cut = chgtaper
      else if (mode .eq. 'CHGDPL') then
         off = sqrt(chgcut*dplcut)
         cut = sqrt(chgtaper*dpltaper)
      else if (mode .eq. 'DIPOLE') then
         off = dplcut
         cut = dpltaper
      end if
c
c     store the powers of the switching window cutoffs
c
      off2 = off * off
      off3 = off2 * off
      off4 = off2 * off2
      off5 = off2 * off3
      off6 = off3 * off3
      off7 = off3 * off4
      cut2 = cut * cut
      cut3 = cut2 * cut
      cut4 = cut2 * cut2
      cut5 = cut2 * cut3
      cut6 = cut3 * cut3
      cut6 = cut3 * cut3
      cut7 = cut3 * cut4
c
c     get 5th degree multiplicative switching function coefficients
c
      if (cut .lt. off) then
         denom = (off-cut)**5
         c0 = off*off2 * (off2-5.0d0*off*cut+10.0d0*cut2) / denom
         c1 = -30.0d0 * off2*cut2 / denom
         c2 = 30.0d0 * (off2*cut+off*cut2) / denom
         c3 = -10.0d0 * (off2+4.0d0*off*cut+cut2) / denom
         c4 = 15.0d0 * (off+cut) / denom
         c5 = -6.0d0 / denom
      end if
c
c     get 7th degree additive switching function coefficients
c
      if (cut.lt.off .and. mode.eq.'CHARGE') then
         term = 9.3d0 * cut*off / (off-cut)
         denom = cut7 - 7.0d0*cut6*off + 21.0d0*cut5*off2
     &              - 35.0d0*cut4*off3 + 35.0d0*cut3*off4
     &              - 21.0d0*cut2*off5 + 7.0d0*cut*off6 - off7
         denom = term * denom
         f0 = cut3*off3 * (-39.0d0*cut+64.0d0*off) / denom
         f1 = cut2*off2
     &           * (117.0d0*cut2-100.0d0*cut*off-192.0d0*off2) / denom
         f2 = cut*off * (-117.0d0*cut3-84.0d0*cut2*off
     &                   +534.0d0*cut*off2+192.0d0*off3) / denom
         f3 = (39.0d0*cut4+212.0d0*cut3*off-450.0d0*cut2*off2
     &            -612.0d0*cut*off3-64.0d0*off4) / denom
         f4 = (-92.0d0*cut3+66.0d0*cut2*off
     &            +684.0d0*cut*off2+217.0d0*off3) / denom
         f5 = (42.0d0*cut2-300.0d0*cut*off-267.0d0*off2) / denom
         f6 = (36.0d0*cut+139.0d0*off) / denom
         f7 = -25.0d0 / denom
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
c     ##  subroutine temper  --  maintain constant temperature  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "temper" maintains a constant desired temperature by scaling
c     the velocities via coupling to an external temperature bath
c
c
      subroutine temper (dt,temp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'moldyn.i'
      include 'usage.i'
      integer i,j
      real*8  temp,dt,scale
c
c
c     find the scale factor to maintain constant temperature
c
      if (temp .eq. 0.0d0)  temp = 0.1d0
      scale = sqrt(1.0d0 + (dt/tautemp)*(kelvin/temp-1.0d0))
c
c     couple to external temperature bath via velocity scaling
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = scale * v(j,i)
            end do
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine tncg  --  truncated newton optimization method  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "tncg" implements the truncated Newton optimization algorithm;
c     a linear conjugate gradient method is used to achieve the
c     (truncated) solution of the Newton equations; special features
c     include use of sparse Hessian matrix, various types of linear
c     equation preconditioning, use of the explicit Hessian or a
c     finite-difference gradient approximation within PCG iteration,
c     and optional use of the exact Newton search directions;
c     under default operation, the algorithm checks for directions
c     of negative curvature during the PCG iteration in order to
c     prevent convergence to a stationary point having negative
c     eigenvalues; if optimization to a saddle point is desired
c     this check can be removed by disabling "negtest"
c
c     literature references:
c
c     J. W. Ponder and F. M Richards, "An Efficient Newton-like
c     Method for Molecular Mechanics Energy Minimization of
c     Large Molecules", Journal of Computational Chemistry,
c     8, 1016-1024 (1987)
c
c     R. S. Dembo and T. Steihaug, "Truncated-Newton Algorithms
c     for Large-Scale Unconstrained Optimization", Mathematical
c     Programming, 26, 190-212 (1983)
c
c     variables and parameters:
c
c     mode       determines optimization method; choice of
c                  Newton's method, truncated Newton, or
c                  truncated Newton with finite differencing
c     method     determines which type of preconditioning will
c                  be used on the Newton equations; choice
c                  of none, diagonal, 3x3 block diagonal,
c                  SSOR or incomplete Cholesky preconditioning
c     nvar       number of parameters in the objective function
c     minimum    upon return contains the best value of the
c                  function found during the optimization
c     f          contains current best value of function
c     x          contains starting point upon input, upon
c                  return contains the best point found
c     g          contains gradient of current best point
c     h          contains the Hessian matrix values in an
c                  indexed linear array
c     h_mode     controls amount of Hessian matrix computed;
c                  either the full matrix, diagonal or none
c     h_init     points to the first Hessian matrix element
c                  associated with each parameter
c     h_stop     points to the last Hessian matrix element
c                  associated with each parameter
c     h_index    contains second parameter involved in each
c                  element of the Hessian array
c     h_diag     contains diagonal of the Hessian matrix
c     p          search direction resulting from pcg iteration
c     f_move     function decrease over last tncg iteration
c     f_old      function value at end of last iteration
c     x_move     rms movement per atom over last tn iteration
c     x_old      parameters value at end of last tn iteration
c     g_norm     Euclidian norm of the gradient vector
c     g_rms      root mean square gradient value
c     fg_call    cumulative number of function/gradient calls
c     grdmin     termination criterion based on RMS gradient
c     iprint     print iteration results every iprint iterations
c     iwrite     call user-supplied output every iwrite iterations
c     newhess    number of iterations between the computation
c                  of new Hessian matrix values
c     negtest    determines whether test for negative curvature
c                  is performed during the PCG iterations
c     maxiter    maximum number of tncg iterations to attempt
c
c     parameters used in the line search:
c
c     cappa      accuarcy of line search control  (0 < cappa < 1)
c     stpmin     minimum allowed line search step size
c     stpmax     maximum allowed line search step size
c     angmax     maximum angle between search and -grad directions
c     intmax     maximum number of interpolations in line search
c
c     required external routines:
c
c     fgvalue    function to evaluate function and gradient values
c     hmatrix    subroutine which evaluates Hessian diagonal
c                  and large off-diagonal matrix elements
c     writeout   subroutine to write out info about current status
c
c
      subroutine tncg (mode,method,nvar,x,minimum,grdmin,
     &                     fgvalue,hmatrix,writeout,h,h_index,
     &                     ff,c_index,c_value,maxhess)
casa
       implicit none
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'sizes.i'
      include 'hescut.i'
      include 'iounit.i'
      include 'keys.i'
      include 'linmin.i'
      include 'math.i'
      include 'minima.i'
      include 'output.i'
      include 'pistuf.i'
      include 'potent.i'
      include 'solute.i'
casa
      include 'restrn.i'
      integer   atm
casa
      integer i,nvar,fg_call,next,maxhess
      integer iter_tn,iter_cg,newhess
      integer h_init(maxvar),h_stop(maxvar)
      integer h_index(maxhess)
      real*8 x(maxvar),x_old(maxvar)
      real*8 g(maxvar),p(maxvar)
      real*8 h(maxhess),h_diag(maxvar)
      real*8 minimum,f,fgvalue,grdmin,angle,rms
      real*8 x_move,f_move,f_old,g_norm,g_rms
      real*8 ff(maxhess)
      integer c_index(maxhess),c_value(maxhess)
      character*4 h_mode
      character*6 mode,method
      character*9 status,info_solve,info_search
      character*20 keyword
      character*80 record
      logical done,negtest
      logical automode,automatic
      external fgvalue,hmatrix,writeout
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     check number of variables and get type of optimization
c
      if (nvar .gt. maxvar) then
         if (maswrk) write (iout,10)
   10    format (' TNCG  --  Too many Parameters,',
     &           ' Increase value of MAXVAR')
         return
      end if
      rms = sqrt(dble(nvar))
      if (coordtype .eq. 'cartesian') then
         rms = rms / sqrt(3.0d0)
      else if (coordtype .eq. 'rigidbody') then
         rms = rms / sqrt(6.0d0)
      end if
c
c     set default parameters for the optimization
c
      if (fctmin .eq. 0.0d0)  fctmin = -1000000.0d0
      if (iwrite .lt. 0)  iwrite = 1
      if (iprint .lt. 0)  iprint = 1
      if (maxiter .eq. 0)  maxiter = 1000
      if (nextiter .eq. 0)  nextiter = 1
      newhess = 1
      done = .false.
      status = '          '
      negtest = .true.
      automode = .false.
      automatic = .false.
      if (mode .eq. 'auto')  automode = .true.
      if (method .eq. 'auto')  automatic = .true.
c
c     set default parameters for the line search
c
      if (cappa .eq. 0.0d0)  cappa = 0.1d0
      if (stpmin .eq. 0.0d0)  stpmin = 1.0d-20
      if (stpmax .eq. 0.0d0)  stpmax = 5.0d0
      if (angmax .eq. 0.0d0)  angmax = 180.0d0
      if (intmax .eq. 0)  intmax = 8
c
c     search each line of the keyword file for options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'FCTMIN ') then
            read (record(next:80),*,err=20)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (record(next:80),*,err=20)  maxiter
         else if (keyword(1:9) .eq. 'NEXTITER ') then
            read (record(next:80),*,err=20)  nextiter
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (record(next:80),*,err=20)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (record(next:80),*,err=20)  iwrite
         else if (keyword(1:8) .eq. 'NEWHESS ') then
            read (record(next:80),*,err=20)  newhess
         else if (keyword(1:12) .eq. 'SADDLEPOINT ') then
            negtest = .false.
         else if (keyword(1:6) .eq. 'CAPPA ') then
            read (record(next:80),*,err=20)  cappa
         else if (keyword(1:8) .eq. 'STEPMIN ') then
            read (record(next:80),*,err=20)  stpmin
         else if (keyword(1:8) .eq. 'STEPMAX ') then
            read (record(next:80),*,err=20)  stpmax
         else if (keyword(1:7) .eq. 'ANGMAX ') then
            read (record(next:80),*,err=20)  angmax
         else if (keyword(1:7) .eq. 'INTMAX ') then
            read (record(next:80),*,err=20)  intmax
         end if
   20    continue
      end do
c
c     initialize iteration counter and set its maximum value
c
      iter_tn = nextiter - 1
      maxiter = iter_tn + maxiter
c
c     print header information about the method used
c
      if (iprint .gt. 0) then
         if (mode .eq. 'newton') then
            if (maswrk) write (iout,30)
   30       format (/,' Full-Newton Conjugate-Gradient',
     &                ' Optimization :')
         else if (mode .eq. 'tncg') then
            if (maswrk) write (iout,40)
   40       format (/,' Truncated-Newton Conjugate-Gradient',
     &                ' Optimization :')
         else if (mode .eq. 'dtncg') then
            if (maswrk) write (iout,50)
   50       format (/,' Finite-Difference Truncated-Newton',
     &                ' Conjugate-Gradient Optimization :')
         else if (mode .eq. 'auto') then
            if (maswrk) write (iout,60)
   60       format (/,' Variable-Mode Truncated-Newton',
     &                ' Conjugate-Gradient Optimization :')
         end if
         if (maswrk) write (iout,70)  mode,method,grdmin
   70    format (/,' Algorithm : ',a6,5x,'Preconditioning : ',a6,5x,
     &             ' RMS Grad : ',d8.2)
         if (maswrk) write (iout,80)
   80    format (/,' TN Iter   F Value       G RMS     F Move  ',
     &             '  X Move   CG Iter   Solve   FG Call',/)
      end if
c
c     compute and print initial function and gradient values
c
      iter_cg = 0
      fg_call = 1
      f = fgvalue (x,g)
      f_old = f
      g_norm = 0.0d0
casa
c     write(6,*) 'wwwf3',npfix,(pfix(i),i=1,npfix),(ipfix(i),i=1,npfix)
      do i = 1, npfix 
         if (pfix(i) .ge. RESTRAINV) then
            atm = (ipfix(i) - 1) * 3
c       write(iout,*) 'atm and g is',  atm, g(atm+1), g(atm+2), g(atm+3)
            g(atm+1) = 0.0d0
            g(atm+2) = 0.0d0
            g(atm+3) = 0.0d0
         end if
      end do
casa
      do i = 1, nvar
         x_old(i) = x(i)
         g_norm = g_norm + g(i)**2
      end do
      g_norm = sqrt(g_norm)
      f_move = 0.5d0 * stpmax * g_norm
      g_rms = g_norm / rms
      if (iprint .gt. 0) then
         if (f.lt.1.0d6 .and. f.gt.-1.0d5 .and. g_rms.lt.1.0d6) then
            if (maswrk) write (iout,90)  iter_tn,f,g_rms,fg_call
   90       format (i6,f12.4,f12.4,41x,i7)
         else
            if (maswrk) write (iout,100)  iter_tn,f,g_rms,fg_call
  100       format (i6,d12.4,d12.4,41x,i7)
         end if
      end if
c
c     check for termination criteria met by initial point
c
      if (g_rms .le. grdmin) then
         done = .true.
         minimum = f
         if (iprint .gt. 0) then
            if (maswrk) write (iout,110)
  110       format (/,' TNCG  --  Normal Termination due to SmallGrad')
         end if
      else if (f .le. fctmin) then
         done = .true.
         minimum = f
         if (iprint .gt. 0) then
            if (maswrk) write (iout,120)
  120       format (/,' TNCG  --  Normal Termination due to SmallFct')
         end if
      else if (iter_tn .ge. maxiter) then
         done = .true.
         minimum = f
         if (iprint .gt. 0) then
            if (maswrk) write (iout,130)
  130       format (/,' TNCG  --  Incomplete Convergence',
     &                 ' due to IterLimit')
         end if
      end if
c
c     beginning of the truncated Newton cycle;
c     increment the major iteration counter
c
      dowhile (.not. done)
         iter_tn = iter_tn + 1
c
c     if pisystem is present, do the molecular orbital
c     computation, then recompute function and gradient
c
         if (norbit .ne. 0) then
            use_orbit = .false.
            if (iter_tn .ne. nextiter) then
               call piscf
               fg_call = fg_call + 1
               f = fgvalue (x,g)
            end if
         end if
c
c     if using GB/SA solvation, update Born radii occasionally
c
         if (use_gbsa) then
            if (mod(iter_tn,5) .eq. 0) then
               reborn = 1
               call born
            end if
            reborn = 0
         end if
c
c     choose the optimization mode based on the gradient value
c
         if (automode) then
            if (g_rms .ge. 3.0d0) then
               mode = 'tncg'
            else
               mode = 'dtncg'
            end if
         end if
c
c     decide on an optimal preconditioning based on the gradient
c
         if (automatic) then
            if (g_rms .ge. 10.0d0) then
               method = 'diag'
               hesscut = 1.0d0
            else if (nvar .lt. 10) then
               method = 'ssor'
               hesscut = 1.0d0
            else if (g_rms .lt. 1.0d0) then
               method = 'iccg'
               hesscut = 0.001d0 * nvar
               if (hesscut .gt. 0.1d0)  hesscut = 0.1d0
            else
               method = 'iccg'
               hesscut = 0.001d0 * nvar
               if (hesscut .gt. 1.0d0)  hesscut = 1.0d0
            end if
         end if
c
c     compute needed portions of the Hessian matrix
c
         h_mode = 'full'
         if (mod(iter_tn-1,newhess) .ne. 0)  h_mode = 'none'
         if (mode.eq.'dtncg' .and. method.eq.'none')  h_mode = 'none'
         if (mode.eq.'dtncg' .and. method.eq.'diag')  h_mode = 'diag'
         call hmatrix (h_mode,x,h,h_init,h_stop,h_index,h_diag,maxhess)
c
c     find the next approximate Newton search direction
c
         call solve (mode,method,negtest,nvar,p,x,g,h,
     &               h_init,h_stop,h_index,h_diag,iter_tn,
     &               iter_cg,fg_call,fgvalue,info_solve, 
     &               ff,c_index,c_value,maxhess)
c
c     perform a line search in the chosen direction
c
         info_search = '         '
         call search (nvar,f,g,x,p,f_move,angle,fg_call,
     &                fgvalue,info_search)
c
c     update variables to reflect this iteration
c
casa
c     write(6,*) 'wwwf4',npfix,(pfix(i),i=1,npfix),(ipfix(i),i=1,npfix)
         do i = 1, npfix
            if (pfix(i) .ge. RESTRAINV) then
               atm = (ipfix(i) - 1) * 3
c       write(iout,*) 'atm and g is',  atm, g(atm+1), g(atm+2), g(atm+3)
               g(atm+1) = 0.0d0
               g(atm+2) = 0.0d0
               g(atm+3) = 0.0d0
            end if
         end do
casa
         f_move = f_old - f
         f_old = f
         x_move = 0.0d0
         g_norm = 0.0d0
         do i = 1, nvar
            x_move = x_move + (x(i)-x_old(i))**2
            x_old(i) = x(i)
            g_norm = g_norm + g(i)**2
         end do
         x_move = sqrt(x_move)
         x_move = x_move / rms
         if (coordtype .eq. 'internal') then
            x_move = x_move * radian
         end if
         g_norm = sqrt(g_norm)
         g_rms = g_norm / rms
c
c     quit if the maximum number of iterations is exceeded
c
         if (iter_tn .ge. maxiter) then
            done = .true.
            status = 'IterLimit'
         end if
c
c     quit if the function value did not change
c
         if (f_move .eq. 0.0d0) then
            done = .true.
            status = 'NoMotion '
         end if
c
c     quit if either of the normal termination tests are met
c
         if (g_rms .le. grdmin) then
            done = .true.
            status = 'SmallGrad'
         else if (f .le. fctmin) then
            done = .true.
            status = 'SmallFct '
         end if
c
c     set warning if the line search encountered problems
c
         if (info_search .ne. ' Success ') then
            info_solve = info_search
         end if
c
c     print the results for the current iteration
c
         if (iprint .gt. 0) then
            if (done .or. mod(iter_tn,iprint) .eq. 0) then
               if (f.lt.1.0d6 .and. f.gt.-1.0d5 .and.
     &             g_rms.lt.1.0d6 .and. f_move.lt.1.0d5) then
                  if (maswrk) write (iout,140)  iter_tn,f,g_rms,f_move,
     &                        x_move,iter_cg,info_solve,fg_call
  140             format (i6,f12.4,f12.4,f11.4,f10.4,i8,3x,a9,i7)
               else
                  if (maswrk) write (iout,150)  iter_tn,f,g_rms,f_move,
     &                        x_move,iter_cg,info_solve,fg_call
  150             format (i6,d12.4,d12.4,d11.4,f10.4,i8,3x,a9,i7)
               end if
            end if
         end if
         if (iwrite .gt. 0) then
            if (done .or. mod(iter_tn,iwrite) .eq. 0) then
               if (maswrk) call writeout (x,iter_tn)
            end if
         end if
c
c     print the reason for terminating the optimization
c
         if (done) then
            minimum = f
            if (iprint .gt. 0) then
               if (g_rms.le.grdmin .or. f.le.fctmin) then
                  if (maswrk) write (iout,160)  status
  160             format (/,' TNCG  --  Normal Termination due to ',a9)
               else
                  if (maswrk) write (iout,170)  status
  170             format (/,' TNCG  --  Incomplete Convergence',
     &                       ' due to ',a9)
               end if
            end if
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine torphase  --  torsional amplitude and phase  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torphase" sets the n-fold amplitude and phase values
c     for each torsion via sorting of the input parameters
c
c
      subroutine torphase (ft,vt,st)
      implicit none
      integer i,k
      integer ft(6)
      real*8 vt(6),st(6)
      real*8 ampli(6),phase(6)
c
c
c     copy the input fold, amplitude and phase angles
c
      do i = 1, 6
         ampli(i) = vt(i)
         phase(i) = st(i)
         vt(i) = 0.0d0
         st(i) = 0.0d0
      end do
c
c     shift the phase angles into the standard range
c
      do i = 1, 6
         dowhile (phase(i) .lt. -180.0d0)
            phase(i) = phase(i) + 360.0d0
         end do
         dowhile (phase(i) .gt. 180.0d0)
            phase(i) = phase(i) - 360.0d0
         end do
      end do
c
c     convert input torsional parameters to storage format
c
      do i = 1, 6
         k = ft(i)
         if (k .eq. 0) then
            goto 10
         else if (k .le. 6) then
            vt(k) = ampli(i)
            st(k) = phase(i)
         end if
      end do
   10 continue
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
c     ##  subroutine torsions  --  locate and store torsions  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "torsions" finds the total number of dihedral angles and
c     the numbers of the four atoms defining each dihedral angle
c
c
      subroutine torsions
      implicit none
      include 'sizes.i'
      include 'bond.i'
      include 'couple.i'
      include 'iounit.i'
      include 'tors.i'
      integer i,j,k,ia,ib,ic,id
c
c
c     loop over all bonds, storing the atoms in each torsion
c
      ntors = 0
      do i = 1, nbond
         ib = ibnd(1,i)
         ic = ibnd(2,i)
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia .ne. ic) then
               do k = 1, n12(ic)
                  id = i12(k,ic)
                  if (id.ne.ib .and. id.ne.ia) then
                     ntors = ntors + 1
                     if (ntors .gt. maxtors) then
                        write (iout,10)
   10                   format (/,' TORSIONS  --  Too many Torsional',
     &                             ' Angles; Increase MAXTORS')
                        call fatal
                     end if
                     itors(1,ntors) = ia
                     itors(2,ntors) = ib
                     itors(3,ntors) = ic
                     itors(4,ntors) = id
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
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  function trimtext  --  find last non-blank character  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "trimtext" finds and returns the location of the last
c     non-blank character before the first null character in
c     an input text string; the function returns zero if no
c     such character is found
c
c
      function trimtext (string)
      implicit none
      integer i,size,len,last,trimtext
      character*1 null
      character*(*) string
c
c
c     move forward through the string, one character
c     at a time, looking for first null character
c
      trimtext = 0
      null = char(0)
      size = len(string)
      last = size
      do i = 1, size
         if (string(i:i) .eq. null) then
            last = i - 1
            goto 10
         end if
      end do
   10 continue
c
c     move backward through the string, one character
c     at a time, looking for first non-blank character
c
      do i = last, 1, -1
         if (string(i:i) .ne. ' ') then
            trimtext = i
            goto 20
         end if
      end do
   20 continue
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
c     ##  subroutine trust  --  update the model trust region  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "trust" updates the model trust region for a nonlinear
c     least squares calculation; this version is based on the
c     ideas found in NL2SOL and in Dennis and Schnabel's book
c
c     arguments and variables :
c
c     m        number of functions
c     n        number of variables
c     xc       vector of length n containing the current iterate
c     fcnorm   real scalar containing the norm of f(xc)
c     gc       vector of length n containing the gradient at xc
c     a        real m by n matrix containing the upper triangular
c                matrix r from the QR factorization of the current
c                Jacobian in the upper triangle
c     lda      leading dimension of A exactly as specified in
c                the dimension statement of the calling program
c     ipvt     integer vector of length n containing the permutation
c                matrix from QR factorization of the Jacobian
c     sc       vector of length n containing the Newton step
c     sa       vector of length n containing current step
c     xscale   vector of length n containing the diagonal
c                scaling matrix for x
c     gauss    logical variable equal is true when the Gauss-Newton
c                step is taken
c     stpmax   maximum allowable step size
c     delta    trust region radius with value retained between calls
c     icode    return code set upon exit
c                0  means xp accepted as next iterate, delta
c                     is trust region for next iteration
c                1  means the algorithm was unable to find a
c                     satisfactory xp sufficiently distinct from xc
c                2  means both the scaled actual and predicted
c                     function reductions are smaller than rftol
c                3  means that false convergence is detected
c                4  means fpnorm is too large, current iteration is
c                     continued with a new, reduced trust region
c                5  means fpnorm is sufficiently small, but the
c                     chance of taking a longer successful step
c                     seems good that the current iteration is to
c                     be continued with a new, doubled trust region
c     xpprev   vector of length n containing the value of xp
c                at the previous call within this iteration
c     fpprev   vector of length m containing f(xpprev)
c     xp       vector of length n containing the new iterate
c     fp       vector of length m containing the functions at xp
c     fpnorm   scalar containing the norm of f(xp)
c     bigstp   logical variable; true, if maximum step length was taken,
c                false  otherwise
c     ncalls   number of function evaluations used
c     xlo      vector of length n containing the lower bounds
c     xhi      vector of length n containing the upper bounds
c     nactive  number of columns in the active Jacobian
c
c     required external routines :
c
c     rsdvalue   subroutine to evaluate residual function values
c
c
      subroutine trust (rsdvalue,m,n,xc,fcnorm,gc,a,lda,ipvt,sc,sa,
     &                  xscale,gauss,stpmax,delta,icode,xp,fc,fp,
     &                  fpnorm,bigstp,ncalls,xlo,xhi,nactive,stpmin,
     &                  rftol,faketol)
      implicit none
      integer maxlsq,maxrsd
      parameter (maxlsq=50)
      parameter (maxrsd=100)
      integer i,j,k,m,n,lda,icode,ncalls,nactive,ipvt(maxlsq)
      real*8 fcnorm,stpmax,delta,fpnorm
      real*8 xc(maxlsq),gc(maxlsq)
      real*8 a(lda,maxlsq),sc(maxlsq),sa(maxlsq)
      real*8 xp(maxlsq),fc(maxrsd),fp(maxrsd)
      real*8 xlo(maxlsq),xhi(maxlsq),xscale(maxlsq)
      real*8 reduce,predict,rellen,slope,tiny
      real*8 stplen,stpmin,rftol,faketol
      real*8 alpha,temp,precise
      real*8 xpprev(maxlsq),fpprev(maxrsd),fpnrmp
      logical gauss,bigstp,feas,ltemp
      save xpprev,fpprev,fpnrmp
      external rsdvalue,precise
c
c
c     set value of alpha, logicals and step length
c
      alpha = 0.0001d0
      bigstp = .false.
      feas = .true.
      stplen = 0.0d0
      do i = 1, n
         stplen = stplen + (xscale(i)*sc(i))**2
      end do
      stplen = sqrt(stplen)
c
c     compute new trial point and new function values
c
      do i = 1, n
         xp(i) = xc(i) + sc(i)
         if (xp(i) .gt. xhi(i)) then
            sc(i) = xhi(i) - xc(i)
            xp(i) = xhi(i)
            feas = .false.
         else if (xp(i) .lt. xlo(i)) then
            sc(i) = xlo(i) - xc(i)
            xp(i) = xlo(i)
            feas = .false.
         end if
      end do
      ncalls = ncalls + 1
      call rsdvalue (m,n,xp,fp)
      fpnorm = 0.0d0
      do i = 1, m
         fpnorm = fpnorm + fp(i)**2
      end do
      fpnorm = 0.5d0 * fpnorm
      reduce = fpnorm - fcnorm
      slope = 0.0d0
      do i = 1, n
         slope = slope + gc(i)*sc(i)
      end do
      if (icode .ne. 5)  fpnrmp = 0.0d0
c
c     internal doubling no good; reset to previous and quit
c
      if (icode.eq.5 .and.
     &     ((fpnorm.ge.fpnrmp).or.(reduce.gt.alpha*slope))) then
         icode = 0
         do i = 1, n
            xp(i) = xpprev(i)
         end do
         do i = 1, m
            fp(i) = fpprev(i)
         end do
         fpnorm = fpnrmp
         delta = 0.5d0 * delta
c
c     fpnorm is too large; the step is unacceptable
c
      else if (reduce .ge. alpha*slope) then
         rellen = 0.0d0
         do i = 1, n
            temp = abs(sc(i))/max(abs(xp(i)),1.0d0/xscale(i))
            rellen = max(rellen,temp)
         end do
c
c     magnitude of (xp-xc) is too small, end the global step
c
         if (rellen .lt. stpmin) then
            icode = 1
            do i = 1, n
               xp(i) = xc(i)
            end do
            do i = 1, m
               fp(i) = fc(i)
            end do
c
c     quadratic interpolation step; reduce delta and continue
c
         else
            icode = 4
            tiny = precise (1)
            if (abs(reduce-slope) .gt. tiny) then
               temp = -slope*stplen / (2.0d0*(reduce-slope))
            else
               temp = -slope*stplen / 2.0d0
            end if
            if (temp .lt. 0.1d0*delta) then
               delta = 0.1d0 * delta
            else if (temp .gt. 0.5d0*delta) then
               delta = 0.5d0 * delta
            else
               delta = temp
            end if
         end if
c
c     fpnorm is sufficiently small; the step is acceptable compute
c     the predicted reduction as predict = g(T)*s + (1/2)*s(T)*h*s
c     with h = p * r**t * r * p**t
c
      else
         predict = slope
         do i = 1, nactive
            k = ipvt(i)
            temp = 0.0d0
            do j = i, nactive
               temp = temp + sa(k)*a(i,j)
            end do
            predict = predict + 0.5d0*temp**2
         end do
         ltemp = (abs(predict-reduce) .le. 0.1d0*abs(reduce))
c
c     if reduce and predict agree to within relative error of 0.1
c     or if negative curvature is indicated, and a longer step is
c     possible and delta has not been decreased this iteration,
c     then double trust region and continue global step
c
         if (icode.ne.4 .and. (ltemp.or.(reduce.le.slope)) .and. feas
     &        .and. .not.gauss .and. (delta.le.0.99d0*stpmax)) then
            icode = 5
            do i = 1, n
               xpprev(i) = xp(i)
            end do
            do i = 1, m
               fpprev(i) = fp(i)
            end do
            fpnrmp = fpnorm
            delta = min(2.0d0*delta,stpmax)
c
c     accept the point; choose new trust region for next iteration
c
         else
            icode = 0
            if (stplen .gt. 0.99d0*stpmax)  bigstp = .true.
            if (reduce .ge. 0.1d0*predict) then
               delta = 0.5d0 * delta
            else if (reduce .le. 0.75d0*predict) then
               delta = min(2.0d0*delta,stpmax)
            end if
         end if
c
c     check relative function convergence and false convergence
c
         if (reduce .le. 2.0d0*predict) then
            if (abs(reduce).le.rftol*abs(fcnorm) .and.
     &          abs(predict).le.rftol*abs(fcnorm)) then
               icode = 2
            end if
         else
            rellen = 0.0d0
            do i = 1, n
               temp = abs(sc(i))/max(abs(xp(i)),1.0d0/xscale(i))
               rellen = max(rellen,temp)
            end do
            if (rellen .lt. faketol)  icode = 3
         end if
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine unitcell  --  read periodic cell parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "unitcell" gets the periodic unitcell sizes and related
c     values from an external keyword file
c
c
      subroutine unitcell
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'keys.i'
      integer i,next
      character*20 keyword
      character*80 record,string
c
c
c     set the default values for the unitcell variables
c
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
      alpha = 0.0d0
      beta = 0.0d0
      gamma = 0.0d0
      octahedron = .false.
      spacegrp = '          '
c
c     get keywords containing crystal lattice dimensions
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:80)
         if (keyword(1:7) .eq. 'A-AXIS ') then
            read (string,*,err=20)  xbox
         else if (keyword(1:7) .eq. 'B-AXIS ') then
            read (string,*,err=20)  ybox
         else if (keyword(1:7) .eq. 'C-AXIS ') then
            read (string,*,err=20)  zbox
         else if (keyword(1:6) .eq. 'ALPHA ') then
            read (string,*,err=20)  alpha
         else if (keyword(1:5) .eq. 'BETA ') then
            read (string,*,err=20)  beta
         else if (keyword(1:6) .eq. 'GAMMA ') then
            read (string,*,err=20)  gamma
         else if (keyword(1:11) .eq. 'OCTAHEDRON ') then
            octahedron = .true.
         else if (keyword(1:11) .eq. 'SPACEGROUP ') then
            call getword (record,spacegrp,next)
         end if
   20    continue
      end do
c
c     assume standard defaults for unspecified cell parameters
c
      if (xbox .ne. 0.0d0) then
         if (ybox .eq. 0.0d0)  ybox = xbox
         if (zbox .eq. 0.0d0)  zbox = xbox
      end if
      if (max(alpha,beta,gamma) .ne. 0.0d0) then
         if (alpha .eq. 0.0d0)  alpha = 90.0d0
         if (beta .eq. 0.0d0)  beta = 90.0d0
         if (gamma .eq. 0.0d0)  gamma = 90.0d0
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
c     ##  subroutine upcase  --  convert string to all upper case  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "upcase" converts a text string to all upper case letters
c
c
      subroutine upcase (string)
      implicit none
      integer i,length,code,ichar
      character*1 char
      character*(*) string
c
c
c     move through the string one character at a time,
c     converting lower case letters to upper case
c
      length = len(string)
      do i = 1, length
         code = ichar(string(i:i))
         if (code.ge.97 .and. code.le.122)
     &      string(i:i) = char(code-32)
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
c     ##  subroutine verlet  --  Verlet molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "verlet" performs a single molecular dynamics time step
c     by means of the velocity Verlet multistep recursion formula
c
c
      subroutine verlet (istep,dt,dt_2,dt2_2,ndump)
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
      real*8 dt,dt_2,dt2_2,e_tot,e_kin,e_pot
      real*8 temp,pres,vol,derivs(3,maxatm)
      real*8 x_old(maxatm),y_old(maxatm),z_old(maxatm)
c
c
c     store the current atom positions, then find new atom
c     positions and half-step velocities via Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            x_old(i) = x(i)
            y_old(i) = y(i)
            z_old(i) = z(i)
            x(i) = x(i) + v(1,i)*dt + a(1,i)*dt2_2
            y(i) = y(i) + v(2,i)*dt + a(2,i)*dt2_2
            z(i) = z(i) + v(3,i)*dt + a(3,i)*dt2_2
            do j = 1, 3
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
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
c     find the full-step velocities using the Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
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
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine version  --  create version number for file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "version" checks the name of a file about to be opened; if
c     if "old" status is passed, the name of the highest current
c     version is returned; if "new" status is passed the filename
c     of the next available unused version is generated
c
c
      subroutine version (filename,status)
      implicit none
      include 'iounit.i'
      include 'output.i'
      integer i,leng,trimtext
      integer thousand,hundred,tens,ones
      character*1 digit(0:9)
      character*3 status
      character*60 filename,oldfile,newfile
      logical exist
      data digit / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     process the filename and status variables
c
      call lowcase (status)
      leng = trimtext (filename)
c
c     no change is needed if the file doesn't exist
c
      exist = .false.
      if (leng .ne. 0)  inquire (file=filename(1:leng),exist=exist)
      if (.not. exist)  return
c
c     set initial values for the current and next versions
c
      newfile = filename
      oldfile = filename
c
c     append an artificial version number to the filename;
c     currently handles up to 10000 versions of a file
c
      if (use_version) then
         i = 1
         dowhile (exist)
            i = i + 1
            oldfile = newfile
            thousand = i / 1000
            hundred = (i - 1000*thousand) / 100
            tens = (i - 1000*thousand - 100*hundred) / 10
            ones = i - 1000*thousand - 100*hundred - 10*tens
            if (thousand .ne. 0) then
               newfile = filename(1:leng)//'_'//digit(thousand)
     &                      //digit(hundred)//digit(tens)//digit(ones)
            else if (hundred .ne. 0) then
               newfile = filename(1:leng)//'_'//digit(hundred)
     &                      //digit(tens)//digit(ones)
            else if (tens .ne. 0) then
               newfile = filename(1:leng)//'_'//digit(tens)//digit(ones)
            else
               newfile = filename(1:leng)//'_'//digit(ones)
            end if
            inquire (file=newfile,exist=exist)
         end do
      end if
c
c     set the file name based on the requested status
c
      if (status .eq. 'old') then
         filename = oldfile
      else if (status .eq. 'new') then
         filename = newfile
         inquire (file=filename,exist=exist)
         if (exist) then
            call nextarg (filename,exist)
            if (exist) then
               inquire (file=filename,exist=exist)
            else
               exist = .true.
            end if
            dowhile (exist)
               write (iout,20)
   20          format (/,' Enter File Name for Coordinate Output :  ',$)
               read (input,30)  filename
   30          format (a60)
               inquire (file=filename,exist=exist)
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
c     ############################################################
c     ##                                                        ##
c     ##  subroutine vmetric  --  variable metric optimization  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "vmetric" is an optimally conditioned variable metric
c     function minimization routine without line searches
c
c     literature references:
c
c     W. C. Davidon, "Optimally Conditioned Optimization Algorithms
c     Without Line Searches", Mathematical Programming, 9, 1-30 (1975)
c
c     D. F. Shanno and K-H. Phua, "Matrix Conditioning and Nonlinear
c     Optimization", Mathematical Programming, 14, 149-16 (1977)
c
c     D. F. Shanno and K-H. Phua, "Numerical Comparison of Several
c     Variable-Metric Algorithms", Journal of Optimization Theory
c     and Applications, 25, 507-518 (1978)
c
c     variables and parameters:
c
c     nvar      number of parameters in the objective function
c     x0        contains starting point upon input, upon return
c                 contains the best point found
c     f0        during optimization contains best current function
c                 value; returns final best function value
c     hguess    initial value for diagonal elements of the "h" matrix
c     eps       maximum length of vector considered infinitely small
c     grdmin    normal exit if rms gradient gets below this value
c     stpmax    maximum length of DeltaX vector in any one iteration
c     fctmin    lower bound on the desired minimum function value
c     maxiter   maximum number of iterations to be attempted
c     iprint    print iteration results every iprint iterations
c     iwrite    call user-supplied output every iwrite iterations
c     ncalls    total number of function/gradient evaluations
c     gnorm     returns with the final value of the gradient norm
c     scale     factor by which the actual function has been multiplied
c
c     required external routines:
c
c     fgvalue    function to evaluate function and gradient values
c     writeout   subroutine to write out info about current status
c
c
      subroutine vmetric (nvar,x0,f0,grdmin,fgvalue,writeout)
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'keys.i'
      include 'linmin.i'
      include 'math.i'
      include 'minima.i'
      include 'output.i'
      include 'potent.i'
      include 'scales.i'
      include 'solute.i'
      integer i,j,nvar,mvar,ncalls,next
      integer niter,nbigang,maxbigang
      real*8 fgvalue,eps,grdmin,precise
      real*8 f0,f,fprime,f0old,f0prime,srchnorm
      real*8 sgangle,sg,snorm,zeta,cosang
      real*8 fmove,xmove,gnorm,grms,rms
      real*8 m2,n2,u2,v,micron,mw,us,qk0
      real*8 a,b,b0,c,alpha,gamma,delta
      real*8 x(maxopt),x0(maxopt),x0old(maxopt)
      real*8 g(maxopt),search(maxopt)
      real*8 s(maxopt),w(maxopt)
      real*8 k(maxopt),k0(maxopt)
      real*8 m(maxopt),n(maxopt),u(maxopt)
      real*8 p(maxopt),q(maxopt),hq(maxopt)
      real*8 h(maxopt,maxopt)
      character*9 status
      character*20 keyword
      character*80 record
      logical restart,done
      save h
      external fgvalue,writeout
c
c
c     initialization and set-up for the optimization
c
      if (nvar .gt. maxopt) then
         write (iout,10)
   10    format (' VMETRIC  --  Too many Parameters,',
     &           ' Increase value of MAXOPT')
         return
      end if
      mvar = nvar
      rms = sqrt(dble(nvar))
      if (coordtype .eq. 'cartesian') then
         rms = rms / sqrt(3.0d0)
      else if (coordtype .eq. 'rigidbody') then
         rms = rms / sqrt(6.0d0)
      end if
      maxbigang = 2
      eps = precise (2)
      restart = .true.
      done = .false.
c
c     set Born radii update frequency for GB/SA solvation
c
      if (use_gbsa) then
         if (reborn .eq. 1) then
            call born
            reborn = 50
         end if
      end if
c
c     set default values for variable scale factors
c
      if (.not. set_scale) then
         do i = 1, nvar
            if (scale(i) .eq. 0.0d0)  scale(i) = 1.0d0
         end do
      end if
c
c     set default parameters for the optimization
c
      if (fctmin .eq. 0.0d0)  fctmin = -1000000.0d0
      if (maxiter .eq. 0)  maxiter = 1000000
      if (nextiter .eq. 0)  nextiter = 1
      if (iprint .lt. 0)  iprint = 1
      if (iwrite .lt. 0)  iwrite = 1
      if (stpmax .eq. 0.0d0)  stpmax = 5.0d0
      if (angmax .eq. 0.0d0)  angmax = 180.0d0
      if (hguess .eq. 0.0d0)  hguess = 0.4d0
c
c     search the keywords for optimization parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'FCTMIN ') then
            read (record(next:80),*,err=20)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (record(next:80),*,err=20)  maxiter
         else if (keyword(1:9) .eq. 'NEXTITER ') then
            read (record(next:80),*,err=20)  nextiter
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (record(next:80),*,err=20)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (record(next:80),*,err=20)  iwrite
         else if (keyword(1:7) .eq. 'HGUESS ') then
            read (record(next:80),*,err=20)  hguess
         else if (keyword(1:8) .eq. 'STEPMAX ') then
            read (record(next:80),*,err=20)  stpmax
         else if (keyword(1:7) .eq. 'ANGMAX ') then
            read (record(next:80),*,err=20)  angmax
         end if
   20    continue
      end do
c
c     get initial function and gradient values, then print header
c
      niter = nextiter - 1
      maxiter = niter + maxiter
      do i = 1, nvar
         x0old(i) = x0(i)
      end do
      ncalls = 1
      f0 = fgvalue (x0,g)
      f0old = f0
      if (iprint .gt. 0) then
         write (iout,30)
   30    format (/,' Optimally Conditioned Variable Metric',
     &             ' Optimization :')
         write (iout,40)
   40    format (/,' VM Iter    F Value       G RMS     F Move',
     &             '    X Move      Angle   FG Call'/)
      end if
c
c     step 1:  selection of the search direction and step
c
c     set the "h" matrix to a diagonal matrix on first iteration or
c     if the angle between the "s" vector and the negative gradient
c     has been bigger than angmax for maxbigang straight iterations
c
      dowhile (.not. done)
         if (restart) then
            do j = 1, mvar
               do i = 1, nvar
                  h(i,j) = 0.0d0
               end do
            end do
            do j = 1, mvar
               h(j,j) = hguess
            end do
            do j = 1, mvar
               k0(j) = 0.0d0
               do i = 1, nvar
                  k0(j) = k0(j) + h(i,j)*g(i)
               end do
               w(j) = k0(j)
            end do
            restart = .false.
         end if
c
c     start the next iteration using either an updated "h"
c     matrix or the "h" matrix from the previous iteration
c
         gnorm = 0.0d0
         grms = 0.0d0
         do i = 1, nvar
            gnorm = gnorm + g(i)**2
            grms = grms + (g(i)*scale(i))**2
         end do
         gnorm = sqrt(gnorm)
         grms = sqrt(grms) / rms
         xmove = 0.0d0
         if (niter .ne. 0) then
            do i = 1, nvar
               xmove = xmove + ((x0(i)-x0old(i))/scale(i))**2
               x0old(i) = x0(i)
            end do
            xmove = sqrt(xmove) / rms
            if (coordtype .eq. 'internal') then
               xmove = radian * xmove
            end if
            fmove = f0old - f0
            f0old = f0
         end if
c
c     print intermediate results for the current iteration
c
         if (iprint .gt. 0) then
            if (niter .eq. 0) then
               if (f0.lt.1.0d7 .and. f0.gt.-1.0d6 .and.
     &                    grms.lt.1.0d6) then
                  write (iout,50)  niter,f0,grms,ncalls
   50             format (i6,f13.4,f12.4,32x,i9)
               else
                  write (iout,60)  niter,f0,grms,ncalls
   60             format (i6,d13.4,d12.4,32x,i9)
               end if
            else if (mod(niter,iprint) .eq. 0) then
               if (f0.lt.1.0d7 .and. f0.gt.-1.0d6 .and.
     &             grms.lt.1.0d6 .and. fmove.lt.1.0d5) then
                  write (iout,70)  niter,f0,grms,fmove,
     &                             xmove,sgangle,ncalls
   70             format (i6,f13.4,f12.4,f11.4,f10.4,f11.4,i9)
               else
                  write (iout,80)  niter,f0,grms,fmove,
     &                             xmove,sgangle,ncalls
   80             format (i6,d13.4,d12.4,d11.4,f10.4,f11.4,i9)
               end if
            end if
         end if
         if (iwrite .gt. 0) then
            if (niter.ne.0 .and. mod(niter,iwrite).eq.0) then
               call writeout (x,niter)
            end if
         end if
c
c     before starting the next iteration, check to see whether
c     the gradient norm, function decrease or iteration limit
c     termination criteria have been satisfied
c
         if (grms.lt.grdmin .or. f0.lt.fctmin
     &          .or. niter.ge.maxiter) then
            if (iwrite .gt. 0) then
               if (mod(niter,iwrite) .ne. 0) then
                  call writeout (x,niter)
               end if
            end if
            if (iprint .gt. 0) then
               if (niter.ne.0 .and. mod(niter,iprint).ne.0) then
                  if (f0.lt.1.0d7 .and. f0.gt.-1.0d6 .and.
     &                grms.lt.1.0d6 .and. fmove.lt.1.0d5) then
                     write (iout,90)  niter,f0,grms,fmove,
     &                                xmove,sgangle,ncalls
   90                format (i6,f13.4,f12.4,f11.4,f10.4,f11.4,i9)
                  else
                     write (iout,100)  niter,f0,grms,fmove,
     &                                xmove,sgangle,ncalls
  100                format (i6,d13.4,d12.4,d11.4,f10.4,f11.4,i9)
                  end if
               end if
               if (niter .ge. maxiter)  status = 'IterLimit'
               if (f0 .lt. fctmin)  status = 'SmallFct '
               if (grms .lt. grdmin)  status = 'SmallGrad'
               if (status .eq. 'IterLimit') then
                  write (iout,110)  status
  110             format (/,' VMETRIC  --  Incomplete Convergence',
     &                       ' due to ',a9)
               else
                  write (iout,120)  status
  120             format (/,' VMETRIC  --  Normal Termination',
     &                      ' due to ',a9)
               end if
            end if
            done = .true.
            goto 170
         end if
c
c     start of the next iteration
c
         niter = niter + 1
         sg = 0.0d0
         snorm = 0.0d0
         do j = 1, mvar
            s(j) = -k0(j)
            snorm = snorm + s(j)**2
            sg = sg - s(j)*g(j)
         end do
         f0prime = -snorm
         snorm = sqrt(snorm)
         cosang = sg / (snorm*gnorm)
         cosang = min(1.0d0,max(-1.0d0,cosang))
         sgangle = radian * acos(cosang)
         if (sgangle .gt. angmax) then
            nbigang = nbigang + 1
         else
            nbigang = 0
         end if
         zeta = 2.0d0
         if (4.0d0*(f0-fctmin) .lt. -f0prime) then
            do j = 1, mvar
               s(j) = -s(j) * (4.0d0*(f0-fctmin)/f0prime)
            end do
            f0prime = -4.0d0 * (f0-fctmin)
         end if
c
c     step 2:  location of the next starting point
c
  130    continue
         do i = 1, nvar
            search(i) = 0.0d0
         end do
         do j = 1, mvar
            do i = 1, nvar
               search(i) = search(i) + h(i,j)*s(j)
            end do
         end do
         srchnorm = 0.0d0
         do i = 1, nvar
            srchnorm = srchnorm + search(i)**2
         end do
         srchnorm = sqrt(srchnorm)
         if (srchnorm .gt. stpmax) then
            do j = 1, mvar
               s(j) = (stpmax/srchnorm) * s(j)
            end do
            do i = 1, nvar
               search(i) = (stpmax/srchnorm) * search(i)
            end do
            f0prime = (stpmax/srchnorm) * f0prime
            zeta = 0.5d0
         end if
c
c     check here to see if we should done without meeting
c     usual termination criteria; quit if -f0prime is too small
c
         if (-f0prime .lt. eps) then
            if (iwrite .gt. 0) then
               if (mod(niter,iwrite) .ne. 0) then
                  call writeout (x,niter)
               end if
            end if
            if (iprint .gt. 0) then
               if (niter.ne.0 .and. mod(niter,iprint).ne.0) then
                  if (f0.lt.1.0d7 .and. f0.gt.-1.0d6 .and.
     &                       grms.lt.1.0d6) then
                     write (iout,140)  niter,f0,grms,0.0,0.0,
     &                                 sgangle,ncalls
  140                format (i6,f13.4,f12.4,f11.4,f10.4,f11.4,i9)
                  else
                     write (iout,150)  niter,f0,grms,0.0,0.0,
     &                                 sgangle,ncalls
  150                format (i6,d13.4,d12.4,f11.4,f10.4,f11.4,i9)
                  end if
               end if
               status = 'SmallMove'
               write (iout,160)  status
  160          format (/,' VMETRIC  --  Incomplete Convergence',
     &                    ' due to ',a9)
            end if
            done = .true.
            goto 170
         end if
         do i = 1, nvar
            x(i) = x0(i) + search(i)
         end do
         ncalls = ncalls + 1
         f = fgvalue (x,g)
         if (f .ge. f0) then
            do j = 1, mvar
               s(j) = 0.5d0 * s(j)
            end do
            f0prime = 0.5d0 * f0prime
            zeta = 0.5d0
            goto 130
         end if
c
c     step 3:  decide whether to update or take another step
c
         do j = 1, mvar
            k(j) = 0.0d0
            do i = 1, nvar
               k(j) = k(j) + h(i,j)*g(i)
            end do
         end do
         fprime = 0.0d0
         do j = 1, mvar
            fprime = fprime + k(j)*s(j)
         end do
         b0 = fprime - f0prime
         do j = 1, mvar
            m(j) = s(j) + k0(j) - k(j)
            k0(j) = k(j)
         end do
         do i = 1, nvar
            x0(i) = x(i)
         end do
         f0 = f
         f0prime = fprime
         if (b0 .lt. eps) then
            do j = 1, mvar
               s(j) = s(j) * zeta
            end do
            f0prime = f0prime * zeta
            goto 130
         end if
c
c     step 4:  check to see if we need to update
c
         if (nbigang .ge. maxbigang) then
            restart = .true.
            goto 170
         end if
         m2 = 0.0d0
         do j = 1, mvar
            m2 = m2 + m(j)**2
         end do
         if (m2 .lt. eps) then
            goto 170
         end if
         v = 0.0d0
         do j = 1, mvar
            v = v + m(j)*s(j)
         end do
         micron = v - m2
         mw = 0.0d0
         do j = 1, mvar
            mw = mw + m(j)*w(j)
         end do
         do j = 1, mvar
            u(j) = w(j) - m(j)*(mw/m2)
         end do
         u2 = 0.0d0
         do j = 1, mvar
            u2 = u2 + u(j)**2
         end do
         if (m2*u2 .ge. eps) then
            us = 0.0d0
            do j = 1, mvar
               us = us + u(j)*s(j)
            end do
            do j = 1, mvar
               n(j) = u(j)*(us/u2)
            end do
            n2 = us * us/u2
         else
            do j = 1, mvar
               n(j) = 0.0d0
            end do
            n2 = 0.0d0
         end if
c
c     step 5:  test inner product of projected s and del-g
c
         b = n2 + micron * v/m2
         if (b .lt. eps) then
            do j = 1, mvar
               n(j) = s(j) - m(j)*(v/m2)
            end do
            n2 = b0 - micron * v/m2
            b = b0
         end if
c
c     step 6:  set "gamma" and "delta" for the update
c
         if (micron*v .ge. m2*n2) then
            gamma = 0.0d0
            delta = sqrt(v/micron)
         else
            a = b - micron
            c = b + v
            gamma = sqrt((1.0d0-micron*v/(m2*n2))/(a*b))
            delta = sqrt(c/a)
            if (c .lt. a) then
               gamma = -gamma
            end if
         end if
c
c     step 7:  perform the update of the "h" matrix
c
         alpha = v + micron*delta + m2*n2*gamma
         do j = 1, mvar
            p(j) = m(j)*(delta-n2*gamma) + n(j)*(gamma*v)
            q(j) = m(j)*((1.0d0+n2*gamma)/alpha)
     &             - n(j)*(gamma * micron/alpha)
            w(j) = m(j)*(n2*(1.0d0+gamma*micron*v)/alpha)
     &             - n(j)*((1.0d0+delta)*micron*v/alpha)
         end do
         qk0 = 0.0d0
         do j = 1, mvar
            qk0 = qk0 + q(j)*k0(j)
         end do
         do j = 1, mvar
            k0(j) = k0(j) + p(j)*qk0
         end do
         do i = 1, nvar
            hq(i) = 0.0d0
         end do
         do j = 1, mvar
            do i = 1, nvar
               hq(i) = hq(i) + h(i,j)*q(j)
            end do
         end do
         do j = 1, mvar
            do i = 1, nvar
               h(i,j) = h(i,j) + hq(i)*p(j)
            end do
         end do
         if (n2 .le. 0.0d0) then
            do j = 1, mvar
               w(j) = k0(j)
            end do
         end if
  170    continue
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
c     ##  subroutine volume  --  excluded volume term via Connolly  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "volume" calculates the excluded volume via the Connolly
c     analytical volume and surface area algorithm
c
c
      subroutine volume (volume_tot,radius,exclude)
      implicit none
      include 'sizes.i'
      real*8 volume_tot,area_tot,probe
      real*8 radius(maxatm),exclude
c
c
c     make call to the volume and surface area routine
c
      probe = 0.0d0
      call connolly (volume_tot,area_tot,radius,probe,exclude)
      return
      end
c
c
c     ################################################################
c     ##  COPYRIGHT (C) 1990 by Craig Kundrot & Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine volume1  --  excluded volume derivs; Cartesian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "volume1" calculates first derivatives of the total excluded
c     volume with respect to the Cartesian coordinates of each atom
c
c     literature reference:
c
c     C. E. Kundrot, J. W. Ponder and F. M. Richards, "Algorithms for
c     Calculating Excluded Volume and Its Derivatives as a Function
c     of Molecular Conformation and Their Use in Energy Minimization",
c     Journal of Computational Chemistry, 12, 402-409 (1991)
c
c
      subroutine volume1 (dex,radius,probe)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'math.i'
      integer maxcube,maxarc
      parameter (maxcube=15)
      parameter (maxarc=300)
      integer i,j,k,m,narc,nx,ny,nz
      integer istart,istop,jstart,jstop,kstart,kstop,mstart
      integer mstop,inov(maxarc),io,ir,in,isum,icube,itemp
      integer cube(2,maxcube,maxcube,maxcube),itab(maxatm)
      real*8 dex(3,maxatm),probe,zstep,radius(maxatm)
      real*8 vdwrad(maxatm),xmin,ymin,zmin,xmax,ymax,zmax
      real*8 aa,bb,temp,theta1,theta2,dtheta,phi_term
      real*8 seg_dx,seg_dy,seg_dz,pre_dx,pre_dy,pre_dz
      real*8 rinsq,rsecn,rsec2n,rdiff,cosine,alpha,beta,ti,tf
      real*8 ztop,ztopshave,zstart,phi1,phi2,cos_phi1,cos_phi2
      real*8 zgrid,rsec2r,rsecr,pix2,rr,rrx2,rrsq,rmax,edge
      real*8 xr,yr,zr,dist2,vdwsum,arci(maxarc),arcf(maxarc)
      real*8 dx(maxarc),dy(maxarc),dsq(maxarc),d(maxarc)
      logical skip(maxatm)
c
c
c     fix the stepsize in the z-direction; this value sets
c     the accuracy of the numerical derivatives; zstep=0.06
c     is a good balance between compute time and accuracy
c
      zstep = 0.0601d0
c
c     initialize minimum and maximum ranges of atoms
c
      pix2 = 2.0d0 * pi
      rmax = 0.0d0
      xmin = x(1)
      xmax = x(1)
      ymin = y(1)
      ymax = y(1)
      zmin = z(1)
      zmax = z(1)
c
c     assign van der Waals radii to the atoms; note that
c     the radii are incremented by the size of the probe;
c     then get the maximum and minimum ranges of atoms
c
      do i = 1, n
         vdwrad(i) = radius(i)
         if (vdwrad(i) .eq. 0.0d0) then
            skip(i) = .true.
         else
            skip(i) = .false.
            vdwrad(i) = vdwrad(i) + probe
            if (vdwrad(i) .gt. rmax)  rmax = vdwrad(i)
            if (x(i) .lt. xmin)  xmin = x(i)
            if (x(i) .gt. xmax)  xmax = x(i)
            if (y(i) .lt. ymin)  ymin = y(i)
            if (y(i) .gt. ymax)  ymax = y(i)
            if (z(i) .lt. zmin)  zmin = z(i)
            if (z(i) .gt. zmax)  zmax = z(i)
         end if
      end do
c
c     load the cubes based on coarse lattice; first of all
c     set edge length to the maximum diameter of any atom
c
      edge = 2.0d0 * rmax
      nx = int((xmax-xmin)/edge) + 1
      ny = int((ymax-ymin)/edge) + 1
      nz = int((zmax-zmin)/edge) + 1
      if (max(nx,ny,nz) .gt. maxcube) then
         write (iout,10)
   10    format (' VOLUME1  --  Increase the value of MAXCUBE')
         call fatal
      end if
c
c     initialize the coarse lattice of cubes
c
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               cube(1,i,j,k) = 0
               cube(2,i,j,k) = 0
            end do
         end do
      end do
c
c     find the number of atoms in each cube
c
      do m = 1, n
         if (.not. skip(m)) then
            i = int((x(m)-xmin)/edge) + 1
            j = int((y(m)-ymin)/edge) + 1
            k = int((z(m)-zmin)/edge) + 1
            cube(1,i,j,k) = cube(1,i,j,k) + 1
         end if
      end do
c
c     determine the highest index in the array "itab" for the
c     atoms that fall into each cube; the first cube that has
c     atoms defines the first index for "itab"; the final index
c     for the atoms in the present cube is the final index of
c     the last cube plus the number of atoms in the present cube
c
      isum = 0
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               icube = cube(1,i,j,k)
               if (icube .ne. 0) then
                  isum = isum + icube
                  cube(2,i,j,k) = isum
               end if
            end do
         end do
      end do
c
c     "cube(2,,,)" now contains a pointer to the array "itab"
c     giving the position of the last entry for the list of
c     atoms in that cube of total number equal to "cube(1,,,)"
c
      do m = 1, n
         if (.not. skip(m)) then
            i = int((x(m)-xmin)/edge) + 1
            j = int((y(m)-ymin)/edge) + 1
            k = int((z(m)-zmin)/edge) + 1
            icube = cube(2,i,j,k)
            itab(icube) = m
            cube(2,i,j,k) = icube - 1
         end if
      end do
c
c     set "cube(2,,,)" to be the starting index in "itab"
c     for atom list of that cube; and "cube(1,,,)" to be
c     the stop index
c
      isum = 0
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               icube = cube(1,i,j,k)
               if (icube .ne. 0) then
                  isum = isum + icube
                  cube(1,i,j,k) = isum
                  cube(2,i,j,k) = cube(2,i,j,k) + 1
               end if
            end do
         end do
      end do
c
c     process in turn each atom from the coordinate list;
c     first select the potential intersecting atoms
c
      do ir = 1, n
         pre_dx = 0.0d0
         pre_dy = 0.0d0
         pre_dz = 0.0d0
         if (skip(ir))  goto 50
         rr = vdwrad(ir)
         rrx2 = 2.0d0 * rr
         rrsq = rr * rr
         xr = x(ir)
         yr = y(ir)
         zr = z(ir)
c
c     find cubes to search for overlaps of current atom
c
         istart = int((xr-xmin)/edge)
         istop = min(istart+2,nx)
         istart = max(istart,1)
         jstart = int((yr-ymin)/edge)
         jstop = min(jstart+2,ny)
         jstart = max(jstart,1)
         kstart = int((zr-zmin)/edge)
         kstop = min(kstart+2,nz)
         kstart = max(kstart,1)
c
c     load all overlapping atoms into "inov"
c
         io = 0
         do i = istart, istop
            do j = jstart, jstop
               do k = kstart, kstop
                  mstart = cube(2,i,j,k)
                  if (mstart .ne. 0) then
                     mstop = cube(1,i,j,k)
                     do m = mstart, mstop
                        in = itab(m)
                        if (in .ne. ir) then
                           io = io + 1
                           if (io .gt. maxarc) then
                              write (iout,20)
   20                         format (' VOLUME1  --  Increase ',
     &                                ' the value of MAXARC')
                              call fatal
                           end if
                           dx(io) = x(in) - xr
                           dy(io) = y(in) - yr
                           dsq(io) = dx(io)**2 + dy(io)**2
                           dist2 = dsq(io) + (z(in)-zr)**2
                           vdwsum = (rr+vdwrad(in))**2
                           if (dist2.gt.vdwsum .or. dist2.eq.0.0d0) then
                              io = io - 1
                           else
                              d(io) = sqrt(dsq(io))
                              inov(io) = in
                           end if
                        end if
                     end do
                  end if
               end do
            end do
         end do
c
c     determine resolution along the z-axis
c
         if (io .ne. 0) then
            ztop = zr + rr
            ztopshave = ztop - zstep
            zgrid = zr - rr
c
c     half of the part not covered by the planes
c
            zgrid = zgrid + 0.5d0*(rrx2-(int(rrx2/zstep)*zstep))
            zstart = zgrid
c
c     section atom spheres perpendicular to the z axis
c
            dowhile (zgrid .le. ztop)
c
c     "rsecr" is radius of circle of intersection
c     of "ir" sphere on the current sphere
c
               rsec2r = rrsq - (zgrid-zr)**2
               if (rsec2r .lt. 0.0d0)  rsec2r = 0.000001d0
               rsecr = sqrt(rsec2r)
               if (zgrid .ge. ztopshave) then
                  cos_phi1 = 1.0d0
                  phi1 = 0.0d0
               else
                  cos_phi1 = (zgrid + 0.5d0*zstep - zr) / rr
                  phi1 = acos(cos_phi1)
               end if
               if (zgrid .eq. zstart) then
                  cos_phi2 = -1.0d0
                  phi2 = pi
               else
                  cos_phi2 = (zgrid - 0.5d0*zstep - zr) / rr
                  phi2 = acos(cos_phi2)
               end if
c
c     check intersections of neighbor circles
c
               narc = 0
               do k = 1, io
                  in = inov(k)
                  rinsq = vdwrad(in)**2
                  rsec2n = rinsq - (zgrid-z(in))**2
                  if (rsec2n .gt. 0.0d0) then
                     rsecn = sqrt(rsec2n)
                     if (d(k) .lt. rsecr+rsecn) then
                        rdiff = rsecr - rsecn
                        if (d(k) .le. abs(rdiff)) then
                           if (rdiff .lt. 0.0d0) then
                              narc = 1
                              arci(narc) = 0.0d0
                              arcf(narc) = pix2
                           end if
                           goto 40
                        end if
                        narc = narc + 1
                        if (narc .gt. maxarc) then
                           write (iout,30)
   30                      format (' VOLUME1  --  Increase',
     &                             ' the value of MAXARC')
                           call fatal
                        end if
c
c     initial and final arc endpoints are found for intersection
c     of "ir" circle with another circle contained in same plane;
c     the initial endpoint of the enclosed arc is stored in "arci",
c     the final endpoint in "arcf"; get "cosine" via law of cosines
c
                        cosine = (dsq(k)+rsec2r-rsec2n)
     &                                     / (2.0d0*d(k)*rsecr)
                        cosine = min(1.0d0,max(-1.0d0,cosine))
c
c     "alpha" is the angle between a line containing either point
c     of intersection and the reference circle center and the
c     line containing both circle centers; "beta" is the angle
c     between the line containing both circle centers and x-axis
c
                        alpha = acos(cosine)
                        beta = atan2(dy(k),dx(k))
                        if (dy(k) .lt. 0.0d0)  beta = beta + pix2
                        ti = beta - alpha
                        tf = beta + alpha
                        if (ti .lt. 0.0d0)  ti = ti + pix2
                        if (tf .gt. pix2)  tf = tf - pix2
                        arci(narc) = ti
c
c     if the arc crosses zero, then it is broken into two segments;
c     the first ends at two pi and the second begins at zero
c
                        if (tf .lt. ti) then
                           arcf(narc) = pix2
                           narc = narc + 1
                           arci(narc) = 0.0d0
                        end if
                        arcf(narc) = tf
   40                   continue
                     end if
                  end if
               end do
c
c     find the pre-area and pre-forces on this section (band),
c     "pre-" means a multiplicative factor is yet to be applied
c
               if (narc .eq. 0) then
                  seg_dz = pix2 * (cos_phi1**2 - cos_phi2**2)
                  pre_dz = pre_dz + seg_dz
               else
c
c     sort the arc endpoint arrays, each with "narc" entries,
c     in order of increasing values of the arguments in "arci"
c
                  k = 1
                  dowhile (k .lt. narc)
                     aa = arci(k)
                     bb = arcf(k)
                     temp = 1000000.0d0
                     do i = k, narc
                        if (arci(i) .le. temp) then
                           temp = arci(i)
                           itemp = i
                        end if
                     end do
                     arci(k) = arci(itemp)
                     arcf(k) = arcf(itemp)
                     arci(itemp) = aa
                     arcf(itemp) = bb
                     k = k + 1
                  end do
c
c     consolidate arcs by removing overlapping arc endpoints
c
                  temp = arcf(1)
                  j = 1
                  do k = 2, narc
                     if (temp .lt. arci(k)) then
                        arcf(j) = temp
                        j = j + 1
                        arci(j) = arci(k)
                        temp = arcf(k)
                     else if (temp .lt. arcf(k)) then
                        temp = arcf(k)
                     end if
                  end do
                  arcf(j) = temp
                  narc = j
                  if (narc .eq. 1) then
                     narc = 2
                     arcf(2) = pix2
                     arci(2) = arcf(1)
                     arcf(1) = arci(1)
                     arci(1) = 0.0d0
                  else
                     temp = arci(1)
                     do k = 1, narc-1
                        arci(k) = arcf(k)
                        arcf(k) = arci(k+1)
                     end do
                     if (temp.eq.0.0d0 .and. arcf(narc).eq.pix2) then
                        narc = narc - 1
                     else
                        arci(narc) = arcf(narc)
                        arcf(narc) = temp
                     end if
                  end if
c
c     compute the numerical pre-derivative values
c
                  do k = 1, narc
                     theta1 = arci(k)
                     theta2 = arcf(k)
                     if (theta2 .ge. theta1) then
                        dtheta = theta2 - theta1
                     else
                        dtheta = (theta2+pix2) - theta1
                     end if
                     phi_term = phi2 - phi1 -
     &                          0.5d0*(sin(2.0d0*phi2)-sin(2.0d0*phi1))
                     seg_dx = (sin(theta2)-sin(theta1)) * phi_term
                     seg_dy = (cos(theta1)-cos(theta2)) * phi_term
                     seg_dz = dtheta * (cos_phi1**2 - cos_phi2**2)
                     pre_dx = pre_dx + seg_dx
                     pre_dy = pre_dy + seg_dy
                     pre_dz = pre_dz + seg_dz
                  end do
               end if
               zgrid = zgrid + zstep
            end do
         end if
   50    continue
         dex(1,ir) = 0.5d0 * rrsq * pre_dx
         dex(2,ir) = 0.5d0 * rrsq * pre_dy
         dex(3,ir) = 0.5d0 * rrsq * pre_dz
      end do
      return
      end
c
c
c     ################################################################
c     ##  COPYRIGHT (C) 1990 by Craig Kundrot & Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine volume2  --  excluded volume Hessian; Cartesian  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "volume2" calculates second derivatives of the total excluded
c     volume with respect to the Cartesian coordinates of the atoms
c
c     literature reference:
c
c     C. E. Kundrot, J. W. Ponder and F. M. Richards, "Algorithms for
c     Calculating Excluded Volume and Its Derivatives as a Function
c     of Molecular Conformation and Their Use in Energy Minimization",
c     Journal of Computational Chemistry, 12, 402-409 (1991)
c
c
      subroutine volume2 (xhess,yhess,zhess,iatom,radius,probe)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'math.i'
      integer maxarc
      parameter (maxarc=300)
      integer i,j,k,m,iatom,narc
      integer nnear,inear(maxarc),iblock,itemp
      integer in,iaa,ibb,idtemp,idfirst,id(0:2)
      integer arciatom(maxarc),arcfatom(maxarc)
      real*8 radius(maxatm),probe,zstep,xr,yr,zr
      real*8 ztop,ztopshave,zstart,aa,bb,temp
      real*8 phi1,phi2,phiold,theta1,theta2,tempf,firsti
      real*8 zgrid,rsec2r,rsecr,pix2,dist2,rcut2
      real*8 vdwrad(maxatm),arci(maxarc),arcf(maxarc)
      real*8 dx(maxarc),dy(maxarc),dsq(maxarc),d(maxarc)
      real*8 rr,rrx2,rrsq,alpha,beta,gamma,ti,tf,ri,s2
      real*8 rinsq,rsecn,rsec2n,b,cosine,delx(2),dely(2),delz(2)
      real*8 cos1,cos2,sin1,sin2,phi_xy,phi_z
      real*8 r_s(2),r_s2(2),r(0:2),r_r(0:2),u(2)
      real*8 dfdtheta(3,2),dthetadx(2,3,0:2),dbetadx(2,2,0:2)
      real*8 dalphdx(2,3,0:2),duds(2),dudr(2),u_term(2)
      real*8 dudx(2,3,0:2),dsdx(2,2,0:2),drdz(2,0:2)
      real*8 xhess(3,maxatm),yhess(3,maxatm),zhess(3,maxatm)
      logical covered
c
c
c     fix the stepsize in the z-direction; this value sets
c     the accuracy of the numerical derivatives; zstep=0.06
c     is a good balance between compute time and accuracy
c
      zstep = 0.0601d0
c
c     zero out the Hessian elements for current atom
c
      do i = 1, n
         do j = 1, 3
            xhess(j,i) = 0.0d0
            yhess(j,i) = 0.0d0
            zhess(j,i) = 0.0d0
         end do
      end do
      if (radius(iatom) .eq. 0.0d0)  return
      pix2 = 2.0d0 * pi
c
c     assign van der Waals radii to the atoms; note that
c     the radii are incremented by the size of the probe
c
      do i = 1, n
         vdwrad(i) = radius(i)
         if (vdwrad(i) .ne. 0.0d0)  vdwrad(i) = vdwrad(i) + probe
      end do
c
c     set the radius and coordinates for current atom
c
      rr = vdwrad(iatom)
      rrx2 = 2.0d0 * rr
      rrsq = rr**2
      xr = x(iatom)
      yr = y(iatom)
      zr = z(iatom)
c
c     select potential intersecting atoms
c
      nnear = 1
      do j = 1, n
         if (j.ne.iatom .and. vdwrad(j).ne.0.0d0) then
            dx(nnear) = x(j) - xr
            dy(nnear) = y(j) - yr
            dsq(nnear) = dx(nnear)**2 + dy(nnear)**2
            dist2 = dsq(nnear) + (z(j)-zr)**2
            rcut2 = (vdwrad(j) + rr)**2
            if (dist2 .lt. rcut2) then
               d(nnear) = sqrt(dsq(nnear))
               inear(nnear) = j
               nnear = nnear + 1
               if (nnear .gt. maxarc) then
                  write (iout,10)
   10             format (' VOLUME2  --  Increase the value of MAXARC')
                  call fatal
               end if
            end if
         end if
      end do
      nnear = nnear - 1
c
c     determine the z resolution
c
      if (nnear .ne. 0) then
         ztop = zr + rr
         ztopshave = ztop - zstep
         zgrid = zr - rr
c
c     half of the part not covered by the planes
c
         zgrid = zgrid + (0.5d0*(rrx2-(int(rrx2/zstep)*zstep)))
         zstart = zgrid
c
c     section atom spheres perpendicular to the z axis
c
         dowhile (zgrid .le. ztop)
c
c     "rsecr" is radius of current atom sphere on the z-plane
c
            rsec2r = rrsq - (zgrid-zr)**2
            if (rsec2r .lt. 0.0d0) then
               rsec2r = 0.000001d0
            end if
            rsecr = sqrt(rsec2r)
            if (zgrid .ge. ztopshave) then
               phi1 = 0.0d0
            else
               phi1 = acos(((zgrid+0.5d0*zstep)-zr) / rr)
            end if
            if (zgrid .eq. zstart) then
               phi2 = pi
            else
               phi2 = phiold
            end if
c
c     check intersections of neighbor circles
c
            k = 0
            narc = 0
            covered = .false.
            dowhile (.not.covered .and. k.lt.nnear
     &                   .and. narc.lt.maxarc)
               k = k + 1
               in = inear(k)
               rinsq = vdwrad(in)**2
               rsec2n = rinsq - (zgrid-z(in))**2
               if (rsec2n .gt. 0.0d0) then
                  rsecn = sqrt(rsec2n)
                  if (d(k) .lt. rsecr+rsecn) then
                     b = rsecr - rsecn
                     if (d(k) .le. abs(b)) then
                        if (b .lt. 0.0d0) then
                           narc = 1
                           arci(narc) = 0.0d0
                           arcf(narc) = pix2
                           arciatom(narc) = in
                           arcfatom(narc) = in
                           covered = .true.
                        end if
                     else
                        narc = narc + 1
                        if (narc .gt. maxarc) then
                           write (iout,20)
   20                      format (' VOLUME2  -- Increase',
     &                             ' the value of MAXARC')
                           call fatal
                        else
c
c     initial and final arc endpoints are found for intersection
c     of "ir" circle with another circle contained in same plane;
c     the initial endpoint of the enclosed arc is stored in "arci",
c     the final endpoint in "arcf"; get "cosine" via law of cosines
c
                           cosine = (dsq(k)+rsec2r-rsec2n) /
     &                                      (2.0d0*d(k)*rsecr)
                           cosine = min(1.0d0,max(-1.0d0,cosine))
c
c     "alpha" is the angle between a line containing either point
c     of intersection and the reference circle center and the
c     line containing both circle centers; "beta" is the angle
c     between the line containing both circle centers and x-axis
c
                           alpha = acos(cosine)
                           if (dx(k) .eq. 0.0d0) then
                              gamma = 0.5d0 * pi
                           else
                              gamma = atan(abs(dy(k)/dx(k)))
                           end if
                           if (dy(k) .gt. 0.0d0) then
                              if (dx(k) .gt. 0.0d0) then
                                 beta = gamma
                              else
                                 beta = pi - gamma
                              end if
                           else
                              if (dx(k) .gt. 0.0d0) then
                                 beta = pix2 - gamma
                              else
                                 beta = pi + gamma
                              end if
                           end if
c
c     finally, the arc endpoints
c
                           ti = beta - alpha
                           tf = beta + alpha
                           if (ti .lt. 0.0d0)  ti = ti + pix2
                           if (tf .gt. pix2)  tf = tf - pix2
                           arci(narc) = ti
                           arciatom(narc) = in
                           arcfatom(narc) = in
                           if (tf .lt. ti) then
                              arcf(narc) = pix2
                              narc = narc + 1
                              arci(narc) = 0.0d0
                              arciatom(narc) = in
                              arcfatom(narc) = in
                           end if
                           arcf(narc) = tf
                        end if
                     end if
                  end if
               end if
            end do
c
c     find the pre- area and pre- forces on this section (band)
c     through sphere "ir"; the "pre-" means a multiplicative
c     factor is yet to be applied
c
            if (narc .ne. 0) then
c
c     general case; sort arc endpoints
c
               k = 1
               dowhile (k .lt. narc)
                  aa = arci(k)
                  bb = arcf(k)
                  iaa = arciatom(k)
                  ibb = arcfatom(k)
                  temp = 10000000.0d0
                  do i = k, narc
                     if (arci(i) .le. temp) then
                        temp = arci(i)
                        itemp = i
                     end if
                  end do
                  arci(k) = arci(itemp)
                  arcf(k) = arcf(itemp)
                  arciatom(k) = arciatom(itemp)
                  arcfatom(k) = arcfatom(itemp)
                  arci(itemp) = aa
                  arcf(itemp) = bb
                  arciatom(itemp) = iaa
                  arcfatom(itemp) = ibb
                  k = k + 1
               end do
c
c     eliminate overlapping arc endpoints;
c     first, consolidate the occluded arcs
c
               m = 1
               tempf = arcf(1)
               idtemp = arcfatom(1)
               do k = 2, narc
                  if (tempf .lt. arci(k)) then
                     arcf(m) = tempf
                     arcfatom(m) = idtemp
                     m = m + 1
                     arci(m) = arci(k)
                     arciatom(m) = arciatom(k)
                     tempf = arcf(k)
                     idtemp = arcfatom(k)
                  else if (tempf .lt. arcf(k)) then
                     tempf = arcf(k)
                     idtemp = arcfatom(k)
                  end if
               end do
               arcf(m) = tempf
               arcfatom(m) = idtemp
               narc = m
c
c     change occluded arcs to accessible arcs
c
               if (narc .eq. 1) then
                  if (arci(1).eq.0.0d0 .and. arcf(1).eq.pix2) then
                     narc = 0
                  else
                     firsti = arci(1)
                     idfirst = arciatom(1)
                     arci(1) = arcf(1)
                     arciatom(1) = arcfatom(1)
                     arcf(1) = firsti + pix2
                     arcfatom(1) = idfirst
                  end if
               else
                  firsti = arci(1)
                  idfirst = arciatom(1)
                  do k = 1, narc-1
                     arci(k) = arcf(k)
                     arciatom(k) = arcfatom(k)
                     arcf(k) = arci(k+1)
                     arcfatom(k) = arciatom(k+1)
                  end do
c
c     check gap between first and last arcs; if the
c     occluded arc crossed zero, then no accessible arc
c
                  if (firsti.eq.0.0d0 .and. arcf(narc).eq.pix2) then
                     narc = narc - 1
                  else
                     arci(narc) = arcf(narc)
                     arciatom(narc) = arcfatom(narc)
                     arcf(narc) = firsti
                     arcfatom(narc) = idfirst
                  end if
               end if
c
c     setup prior to application of chain rule
c
               do k = 1, narc
                  ri = sqrt(rrsq - (zgrid-zr)**2)
                  do i = 1, 2
                     if (i .eq. 1) then
                        id(1) = arciatom(k)
                     else
                        id(2) = arcfatom(k)
                     end if
                     delx(i) = x(id(i)) - xr
                     dely(i) = y(id(i)) - yr
                     delz(i) = zgrid - z(id(i))
                     s2 = delx(i)**2 + dely(i)**2
                     r_s(i) = 1.0d0 / sqrt(s2)
                     r_s2(i) = r_s(i)**2
                     r(i) = sqrt(vdwrad(id(i))**2 - delz(i)**2)
                     r_r(i) = 1.0d0 / r(i)
                     u(i) = (ri**2+s2-r(i)**2) * (0.5d0*r_s(i) / ri)
                  end do
c
c     apply the chain rule repeatedly
c
                  theta1 = arci(k)
                  theta2 = arcf(k)
                  cos1 = cos(theta1)
                  cos2 = cos(theta2)
                  sin1 = sin(theta1)
                  sin2 = sin(theta2)
                  phi_xy = phi2 - phi1 -
     &                     0.5d0*(sin(2.0d0*phi2)-sin(2.0d0*phi1))
                  phi_z = sin(phi2)**2 - sin(phi1)**2
                  phi_xy = 0.5d0 * rrsq * phi_xy
                  phi_z = 0.5d0 * rrsq * phi_z
                  dfdtheta(1,1) = -cos1 * phi_xy
                  dfdtheta(2,1) = -sin1 * phi_xy
                  dfdtheta(3,1) = -phi_z
                  dfdtheta(1,2) =  cos2 * phi_xy
                  dfdtheta(2,2) =  sin2 * phi_xy
                  dfdtheta(3,2) =  phi_z
                  do i = 1, 2
                     dbetadx(i,1,0) = dely(i) * r_s2(i)
                     dbetadx(i,2,0) = -delx(i) * r_s2(i)
                     dbetadx(i,1,i) = -dbetadx(i,1,0)
                     dbetadx(i,2,i) = -dbetadx(i,2,0)
                  end do
                  do i = 1, 2
                     duds(i) = (1.0d0/ri) - (u(i)*r_s(i))
                     dsdx(i,1,i) = delx(i) * r_s(i)
                     dsdx(i,2,i) = dely(i) * r_s(i)
                     dsdx(i,1,0) = -dsdx(i,1,i)
                     dsdx(i,2,0) = -dsdx(i,2,i)
                     dudr(i) = -r(i) * r_s(i) / ri
                     drdz(i,i) = delz(i) * r_r(i)
                     drdz(i,0) = -drdz(i,i)
                  end do
                  do m = 0, 2
                     do i = 1, 2
                        dudx(i,1,m) = duds(i) * dsdx(i,1,m)
                        dudx(i,2,m) = duds(i) * dsdx(i,2,m)
                        dudx(i,3,m) = dudr(i) * drdz(i,m)
                     end do
                  end do
                  do i = 1, 2
                     u_term(i) = -1.0d0 / sqrt(1.0d0-u(i)**2)
                  end do
                  do j = 1, 3
                     do m = 0, 2
                        do i = 1, 2
                           dalphdx(i,j,m) =
     &                        u_term(i) * dudx(i,j,m)
                        end do
                     end do
                  end do
                  do j = 1, 2
                     do m = 0, 2
                        dthetadx(1,j,m) = dbetadx(1,j,m)
     &                                     + dalphdx(1,j,m)
                        dthetadx(2,j,m) = dbetadx(2,j,m)
     &                                     - dalphdx(2,j,m)
                     end do
                  end do
                  do m = 0, 2
                     dthetadx(1,3,m) = dalphdx(1,3,m)
                     dthetadx(2,3,m) = -dalphdx(2,3,m)
                  end do
c
c     partials with respect to coordinates of serial atom id(m)
c
                  id(0) = iatom
                  do m = 0, 2
                     iblock = id(m)
                     do j = 1, 3
                        xhess(j,iblock) = xhess(j,iblock)
     &                     + dfdtheta(1,1) * dthetadx(1,j,m)
     &                     + dfdtheta(1,2) * dthetadx(2,j,m)
                        yhess(j,iblock) = yhess(j,iblock)
     &                     + dfdtheta(2,1) * dthetadx(1,j,m)
     &                     + dfdtheta(2,2) * dthetadx(2,j,m)
                        zhess(j,iblock) = zhess(j,iblock)
     &                     + dfdtheta(3,1) * dthetadx(1,j,m)
     &                     + dfdtheta(3,2) * dthetadx(2,j,m)
                     end do
                  end do
               end do
            end if
            zgrid = zgrid + zstep
            phiold = phi1
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
c     ##  subroutine writeout  --  update the optimization results  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "writeout" is used by each of the optimization routines
c     to save imtermediate atomic coordinates to a disk file
c
c
      subroutine writeout (xx,ncycle)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'math.i'
      include 'omega.i'
      include 'output.i'
      include 'scales.i'
      include 'usage.i'
      include 'zcoord.i'
      integer i,iwrt,freeunit
      integer nvar,ncycle,lext
      real*8 xx(maxvar)
      character*3 status
      character*7 ext
      character*60 coordfile
c
c
c     nothing to do if coordinate type is undefined
c
      if (coordtype .eq. 'none')  return
c
c     check scaling factors for optimization parameters
c
      if (.not. set_scale) then
         set_scale = .true.
         if (coordtype .eq. 'cartesian') then
            do i = 1, 3*n
               scale(i) = 1.0d0
            end do
         else if (coordtype .eq. 'internal') then
            do i = 1, nomega
               scale(i) = 1.0d0
            end do
         end if
      end if
c
c     transform optimization parameters back to coordinates
c
      if (coordtype .eq. 'cartesian') then
         nvar = 0
         do i = 1, n
            if (use(i)) then
               nvar = nvar + 1
               x(i) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               y(i) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               z(i) = xx(nvar) / scale(nvar)
            end if
         end do
      else if (coordtype .eq. 'internal') then
         do i = 1, nomega
            dihed(i) = xx(i) / scale(i)
            ztors(zline(i)) = dihed(i) * radian
         end do
      end if
c
c     get the name to be used for the coordinates file
c
      if (savecycle) then
         lext = 3
         call numeral (ncycle,ext,lext)
         coordfile = filename(1:leng)//'.'//ext(1:lext)
         status = 'new'
      else
         coordfile = outfile
         status = 'old'
      end if
c
c     open and then update the coordinates file
c
      iwrt = freeunit ()
      call version (coordfile,status)
      open (unit=iwrt,file=coordfile,status=status)
      if (status .eq. 'old')  rewind (unit=iwrt)
      if (coordtype .eq. 'cartesian') then
         call prtxyz (iwrt)
      else if (coordtype .eq. 'internal') then
         call prtint (iwrt)
      else if (coordtype .eq. 'rigidbody') then
         call prtxyz (iwrt)
      end if
      close (unit=iwrt)
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
c     ##  subroutine xyzatm  --  single atom internal to Cartesian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "xyzatm" computes the Cartesian coordinates of a single
c     atom from its defining internal coordinate values
c
c
      subroutine xyzatm (i,ia,bond,ib,angle1,ic,angle2,chiral)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      integer i,ia,ib,ic,chiral
      real*8 bond,angle1,angle2
      real*8 eps,rad1,rad2
      real*8 sin1,cos1,sin2,cos2
      real*8 cosine,sine,sine2
      real*8 xab,yab,zab,rab
      real*8 xba,yba,zba,rba
      real*8 xbc,ybc,zbc,rbc
      real*8 xac,yac,zac,rac
      real*8 xt,yt,zt,xu,yu,zu
      real*8 cosb,sinb,cosg,sing
      real*8 xtmp,ztmp,a,b,c
c
c
c     convert angles to radians, and get their sines and cosines
c
      eps = 0.00000001d0
      rad1 = angle1 / radian
      rad2 = angle2 / radian
      sin1 = sin(rad1)
      cos1 = cos(rad1)
      sin2 = sin(rad2)
      cos2 = cos(rad2)
c
c     if no second site given, place the atom at the origin
c
      if (ia .eq. 0) then
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
c
c     if no third site given, place the atom along the z-axis
c
      else if (ib .eq. 0) then
         x(i) = x(ia)
         y(i) = y(ia)
         z(i) = z(ia) + bond
c
c     if no fourth site given, place the atom in the x,z-plane
c
      else if (ic .eq. 0) then
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         rab = sqrt(xab**2 + yab**2 + zab**2)
         xab = xab / rab
         yab = yab / rab
         zab = zab / rab
         cosb = zab
         sinb = sqrt(xab**2 + yab**2)
         if (sinb .eq. 0.0d0) then
            cosg = 1.0d0
            sing = 0.0d0
         else
            cosg = yab / sinb
            sing = xab / sinb
         end if
         xtmp = bond*sin1
         ztmp = rab - bond*cos1
         x(i) = x(ib) + xtmp*cosg + ztmp*sing*sinb
         y(i) = y(ib) - xtmp*sing + ztmp*cosg*sinb
         z(i) = z(ib) + ztmp*cosb
c
c     general case where the second angle is a dihedral angle
c
      else if (chiral .eq. 0) then
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         rab = sqrt(xab**2 + yab**2 + zab**2)
         xab = xab / rab
         yab = yab / rab
         zab = zab / rab
         xbc = x(ib) - x(ic)
         ybc = y(ib) - y(ic)
         zbc = z(ib) - z(ic)
         rbc = sqrt(xbc**2 + ybc**2 + zbc**2)
         xbc = xbc / rbc
         ybc = ybc / rbc
         zbc = zbc / rbc
         xt = zab*ybc - yab*zbc
         yt = xab*zbc - zab*xbc
         zt = yab*xbc - xab*ybc
         cosine = xab*xbc + yab*ybc + zab*zbc
         sine = sqrt(max(1.0d0-cosine**2,eps))
         if (abs(cosine) .ge. 1.0d0) then
            write (iout,10)  i
   10       format (/,' XYZATM  --  Undefined Dihedral',
     &                 ' Angle at Atom',i6)
         end if
         xt = xt / sine
         yt = yt / sine
         zt = zt / sine
         xu = yt*zab - zt*yab
         yu = zt*xab - xt*zab
         zu = xt*yab - yt*xab
         x(i) = x(ia) + bond * (xu*sin1*cos2 + xt*sin1*sin2 - xab*cos1)
         y(i) = y(ia) + bond * (yu*sin1*cos2 + yt*sin1*sin2 - yab*cos1)
         z(i) = z(ia) + bond * (zu*sin1*cos2 + zt*sin1*sin2 - zab*cos1)
c
c     general case where the second angle is a bond angle
c
      else if (abs(chiral) .eq. 1) then
         xba = x(ib) - x(ia)
         yba = y(ib) - y(ia)
         zba = z(ib) - z(ia)
         rba = sqrt(xba**2 + yba**2 + zba**2)
         xba = xba / rba
         yba = yba / rba
         zba = zba / rba
         xac = x(ia) - x(ic)
         yac = y(ia) - y(ic)
         zac = z(ia) - z(ic)
         rac = sqrt(xac**2 + yac**2 + zac**2)
         xac = xac / rac
         yac = yac / rac
         zac = zac / rac
         xt = zba*yac - yba*zac
         yt = xba*zac - zba*xac
         zt = yba*xac - xba*yac
         cosine = xba*xac + yba*yac + zba*zac
         sine2 = max(1.0d0-cosine**2,eps)
         if (abs(cosine) .ge. 1.0d0) then
            write (iout,20)  i
   20       format (/,' XYZATM  --  Defining Atoms Colinear',
     &                 ' at Atom',i6)
         end if
         a = (-cos2 - cosine*cos1) / sine2
         b = (cos1 + cosine*cos2) / sine2
         c = (1.0d0 + a*cos2 - b*cos1) / sine2
         if (c .gt. eps) then
            c = chiral * sqrt(c)
         else if (c .lt. -eps) then
            c = sqrt((a*xac+b*xba)**2 + (a*yac+b*yba)**2
     &                       + (a*zac+b*zba)**2)
            a = a / c
            b = b / c
            c = 0.0d0
            if (debug) then
               write (iout,30)  ia
   30          format (/,' XYZATM  --  Sum of Bond Angles',
     &                    ' Too Large at Atom',i6)
            end if
         else
            c = 0.0d0
         end if
         x(i) = x(ia) + bond * (a*xac + b*xba + c*xt)
         y(i) = y(ia) + bond * (a*yac + b*yba + c*yt)
         z(i) = z(ia) + bond * (a*zac + b*zba + c*zt)
      end if
      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine locase  --  convert string to all lower case  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "locase" converts a text string to all lower case letters
c     cloned from upcase
c
c
      subroutine locase (string)
      implicit none
      integer i,length,code,ichar
      character*1 char
      character*(*) string
c
c
c     move through the string one character at a time,
c     converting upper case letters to lower case
c
      length = len(string)
      do i = 1, length
         code = ichar(string(i:i))
         if (code.ge.65 .and. code.le.90)
     &      string(i:i) = char(code+32)
      end do
      return
      end
c
