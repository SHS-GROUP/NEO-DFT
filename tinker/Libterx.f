c   TINKER library routines that start with er to ex
c                  T(inker)lib(rary)exo.f
c
c  8 MAY 98 - JRS call rotmat changed to call rotmatt
c                 extra: renamed extrat
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  function erf  --  evaluate the standard error function  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "erf" computes a numerical approximation to the value of
c     the error function via a Chebyshev approximation
c
c
      function erf (x)
      integer mode
      real*8 erf,x,result
c
c
c     compute the error function via Chebyshev fitting
c
      mode = 0
      call erfcore (x,result,mode)
      erf = result
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function erfc  --  evaluate the error function compliment  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erfc" computes a numerical approximation to the value of the
c     error function complement via a Chebyshev approximation
c
c
      function erfc (x)
      integer mode
      real*8 erfc,x,result
c
c
c     compute the error function compliment via Chebyshev fitting
c
      mode = 1
      call erfcore (x,result,mode)
      erfc = result
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erfcore  --  erf and erfc via Chebyshev approx  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erfcore" evaluates erf(x) or erfc(x) for a real argument x;
c     this routine is called as the numerical kernel by two simple
c     functions: "erf" and "erfc"; if "calerf" is called with a
c     mode value of 0 it returns erf, a mode of 1 returns erfc;
c     it uses rational functions that theoretically approximate
c     erf(x) and erfc(x) to at least 18 significant decimal digits
c
c     literature reference:
c
c     W. J. Cody, "Rational Chebyshev Approximations for the Error
c     Function", Mathematics of Computation, 631-638, 1969
c
c     machine-dependent constants:
c
c     xsmall   argument below which erf(x) may be represented by
c              2*x/sqrt(pi) and above which x*x won't underflow;
c              a conservative value is the largest machine number
c              X such that 1.0 + X = 1.0 to machine precision
c
c     xbig     largest argument acceptable for erfc; solution to
c              the equation:  W(x) * (1-0.5/x**2) = XMIN, where
c              W(x) = exp(-x*x)/[x*sqrt(pi)]
c
c     original author:
c
c     W. J. Cody
c     Mathematics and Computer Science Division
c     Argonne National Laboratory
c     Argonne, IL 60439
c
c
      subroutine erfcore (arg,result,mode)
      implicit none
      integer i,mode
      real*8 arg,result
      real*8 x,y,ysq,del
      real*8 zero,one,two,four,half
      real*8 sqrpi,thresh,sixten
      real*8 xnum,xden,xsmall,xbig
      real*8 a(5),b(4),c(9),d(8),p(6),q(5)
c
c     mathematical constants
c
      data zero   / 0.0d0 /
      data one    / 1.0d0 /
      data two    / 2.0d0 /
      data four   / 4.0d0 /
      data half   / 0.5d0 /
      data sqrpi  / 5.6418958354775628695d-1 /
      data thresh / 0.46875d0 /
      data sixten / 16.0d0 /
c
c     machine-dependent constants
c
      data xsmall / 1.11d-16 /
      data xbig   / 26.543d0 /
c
c     coefficients for approximation to erf in first interval
c
      data a / 3.16112374387056560d0,  1.13864154151050156d2,
     &         3.77485237685302021d2,  3.20937758913846947d3,
     &         1.85777706184603153d-1 /
      data b / 2.36012909523441209d1,  2.44024637934444173d2,
     &         1.28261652607737228d3,  2.84423683343917062d3 /
c
c     coefficients for approximation to erfc in second interval
c
      data c / 5.64188496988670089d-1, 8.88314979438837594d0,
     &         6.61191906371416295d1,  2.98635138197400131d2,
     &         8.81952221241769090d2,  1.71204761263407058d3,
     &         2.05107837782607147d3,  1.23033935479799725d3,
     &         2.15311535474403846d-8 /
      data d / 1.57449261107098347d1,  1.17693950891312499d2,
     &         5.37181101862009858d2,  1.62138957456669019d3,
     &         3.29079923573345963d3,  4.36261909014324716d3,
     &         3.43936767414372164d3,  1.23033935480374942d3 /
c
c     coefficients for approximation to erfc in third interval
c
      data p / 3.05326634961232344d-1, 3.60344899949804439d-1,
     &         1.25781726111229246d-1, 1.60837851487422766d-2,
     &         6.58749161529837803d-4, 1.63153871373020978d-2 /
      data q / 2.56852019228982242d0,  1.87295284992346047d0,
     &         5.27905102951428412d-1, 6.05183413124413191d-2,
     &         2.33520497626869185d-3 /
c
c
c     store the argument and its absolute value
c
      x = arg
      y = abs(x)
c
c     evaluate erf for |x| less than 0.46875
c
      if (y .le. thresh) then
         ysq = zero
         if (y .gt. xsmall)  ysq = y * y
         xnum = a(5) * ysq
         xden = ysq
         do i = 1, 3
            xnum = (xnum + a(i)) * ysq
            xden = (xden + b(i)) * ysq
         end do
         result = x * (xnum + a(4)) / (xden + b(4))
         if (mode .ne. 0)  result = one - result
c
c     evaluate erfc for 0.46875 <= |x| <= 4.0
c
      else if (y .le. four) then
         xnum = c(9) * y
         xden = y
         do i = 1, 7
            xnum = (xnum + c(i)) * y
            xden = (xden + d(i)) * y
         end do
         result = (xnum + c(8)) / (xden + d(8))
         ysq = aint(y*sixten) / sixten
         del = (y-ysq) * (y+ysq)
         result = exp(-ysq*ysq) * exp(-del) * result
         if (mode .eq. 0) then
            result = (half - result) + half
            if (x .lt. zero)  result = -result
         else
            if (x .lt. zero)  result = two - result
         end if
c
c     evaluate erfc for |x| greater than 4.0
c
      else
         result = zero
         if (y .lt. xbig) then
            ysq = one / (y * y)
            xnum = p(6) * ysq
            xden = ysq
            do i = 1, 4
               xnum = (xnum + p(i)) * ysq
               xden = (xden + q(i)) * ysq
            end do
            result = ysq * (xnum + p(5)) / (xden + q(5))
            result = (sqrpi -  result) / y
            ysq = aint(y*sixten) / sixten
            del = (y-ysq) * (y+ysq)
            result = exp(-ysq*ysq) * exp(-del) * result
         end if
         if (mode .eq. 0) then
            result = (half - result) + half
            if (x .lt. zero)  result = -result
         else
            if (x .lt. zero)  result = two - result
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function erfinv  --  evaluate the error function inverse  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erfinv" evaluates the inverse of the error function erf for
c     a real argument in the range (-1,1) using a rational function
c     approximation followed by cycles of Newton-Raphson correction
c
c     adapted from the pseudocode for the Matlab function of the
c     same name; Matlab, version 4.2c, March 1995
c
c
      function erfinv (y)
      implicit none
      include 'iounit.i'
      include 'math.i'
      real*8 erfinv,erf,x,y,z
      real*8 a(4),b(4),c(4),d(2)
c
c     coefficients for approximation to erfinv in central range
c
      data a /  0.886226899d0, -1.645349621d0,
     &          0.914624893d0, -0.140543331d0 /
      data b / -2.118377725d0,  1.442710462d0,
     &         -0.329097515d0,  0.012229801d0 /
c
c     coefficients for approximation to erfinv near endpoints
c
      data c / -1.970840454d0, -1.624906493d0,
     &          3.429567803d0,  1.641345311d0 /
      data d /  3.543889200d0,  1.637067800d0 /
c
c
c     get an initial estimate for the inverse error function
c
      if (abs(y) .le. 0.7d0) then
         z = y * y
         x = y * (((a(4)*z+a(3))*z+a(2))*z+a(1))
     &              / ((((b(4)*z+b(3))*z+b(2))*z+b(1))*z+1.0d0)
      else if (y.gt.0.7d0 .and. y.lt.1.0d0) then
         z = sqrt(-log((1.0d0-y)/2.0d0))
         x = (((c(4)*z+c(3))*z+c(2))*z+c(1)) / ((d(2)*z+d(1))*z+1.0d0)
      else if (y.lt.-0.7d0 .and. y.gt.-1.0d0) then
         z = sqrt(-log((1.0d0+y)/2.0d0))
         x = -(((c(4)*z+c(3))*z+c(2))*z+c(1)) / ((d(2)*z+d(1))*z+1.0d0)
      else
         write (iout,10)
   10    format (/,' ERFINV  --  Illegal Argument to Inverse',
     &              ' Error Function')
         call fatal
      end if
c
c     use two steps of Newton-Raphson correction to increase accuracy
c
      x = x - (erf(x) - y) / (2.0d0/sqrtpi * exp(-x**2))
      x = x - (erf(x) - y) / (2.0d0/sqrtpi * exp(-x**2))
      erfinv = x
      return
      end
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1996 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine erxnfld  --  reaction field potential energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "erxnfld" calculates the macroscopic reaction field energy
c
c
      subroutine erxnfld
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'mpole.i'
      include 'shunt.i'
      include 'usage.i'
      integer i,j,k
      integer ii,iz,ix,kk,kz,kx
      real*8 eik,a(3,3)
      real*8 xr,yr,zr,r2
      real*8 rpi(13),rpk(13)
      logical iuse,kuse
c
c
c     zero out the macroscopic reaction field energy
c
      er = 0.0d0
c
c     set the switching function coefficients
c
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
c
      call ijk_pt
c
c     calculate the multipole interaction energy term
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, polsiz(ii)
            rpi(j) = rpole(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
            if (iuse .or. kuse) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, polsiz(kk)
                     rpk(j) = rpole(j,kk)
                  end do
                  call erfik (ii,kk,i,k,rpi,rpk,eik)
                  er = er + eik
               end if
            end if
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erfik  --  reaction field energy of site pair   ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine erfik (ii,kk,i,k,rpi,rpk,eik)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'chgpot.i'
      include 'mpole.i'
      include 'rxnfld.i'
      include 'rxnpot.i'
      include 'units.i'
      integer i,k,ii,kk
      integer isiz,ksiz
      integer m,n1,n2,nn
      integer fii,fi,fj
      integer p_s1,p_s2,p_e1,p_e2
      integer ind1_x(13),ind1_y(13),ind1_z(13)
      integer ind2_x(13),ind2_y(13),ind2_z(13)
      real*8 eik,term,ratio,factor
      real*8 xi,yi,zi,xk,yk,zk
      real*8 size,size2,d,d2,ri2,rk2
      real*8 xi2,yi2,zi2,xk2,yk2,zk2
      real*8 rpi(13),rpk(13)
      real*8 rklij(13,13),d1d2
      real*8 m2t2(13)
c
c
c     integer ind_x,ind_y,ind_z
c     real*8 d1d2_r1,d1d2_r2s,d1d2_r2t,d1d2_r2u,d1d2_r2i
c     real*8 d1d2_r2is,d1d2_r2it,d1d2_r2iu,d1d2_r2j
c     real*8 d1d2_r2js,d1d2_r2jt,d1d2_r2ju,d1d2_r2k,d1d2_r
c
c
c     get numbers of the atoms
c
      isiz = polsiz(ii)
      ksiz = polsiz(kk)
      xi = x(i)
      yi = y(i)
      zi = z(i)
      xk = x(k)
      yk = y(k)
      zk = z(k)
      d = xi*xk + yi*yk + zi*zk
      ri2 = xi*xi + yi*yi + zi*zi
      rk2 = xk*xk + yk*yk + zk*zk
c
c     set highest order of multipoles at each site (M=0, D=1, Q=2)
c
      eik = 0.0d0
      n1 = 2
      n2 = 2
      nn = rfterms
      ratio = rfbulkd / dielec
      factor = electric * (1.0d0-ratio)
      if (i .eq. k)  factor = 0.5d0 * factor
      size = 1.0d0 / rfsize
      size2 = size * size
c
c     get the values of the indices
c
      m = (3**(n1+1)-1)/2
      call rindex (n1,m,ind1_x,ind1_y,ind1_z,p_s1,p_e1)
      m = (3**(n2+1)-1)/2
      call rindex (n2,m,ind2_x,ind2_y,ind2_z,p_s2,p_e2)
c
c     initialize the stored matrix element areas
c
      do fi = 1, p_e1
         do fj = 1, p_e2
            b1(fi,fj) = 0.0d0
            b2(fi,fj) = 0.0d0
         end do
      end do
c
c     explicit formula for the 0th summation term
c
      if (nn .ge. 0) then
         eik = size * rpi(1) * rpk(1) / ratio
         size = size * size2
      end if
c     
c     explicit formula for the 1st summation term
c
      if (nn .ge. 1) then
         b2(1,1) = d
         b2(1,2) = xi
         b2(1,3) = yi
         b2(1,4) = zi
         b2(2,1) = xk
         b2(3,1) = yk
         b2(4,1) = zk
         b2(2,2) = 1.0d0
         b2(3,3) = 1.0d0
         b2(4,4) = 1.0d0
         do fi = 1, 4
            m2t2(fi) = 0.0d0
            do fj = 1, 4
               m2t2(fi) = m2t2(fi) + b2(fi,fj)*rpk(fj)
            end do
         end do
         term = 0.0d0
         do fi = 1, 4
            term = term + rpi(fi)*m2t2(fi)
         end do
         term = 2.0d0 * size * term / (2.0d0*ratio+1.0d0)
         eik = eik + term
         size = size * size2
      end if
c     
c     explicit formula for the 2nd summation term
c
      if (nn .ge. 2) then
         b2(1,1) = (3.0d0*d*d-ri2*rk2) * 0.5d0
         b2(1,2) = 3.0d0*xi*d - xk*ri2
         b2(1,3) = 3.0d0*yi*d - yk*ri2
         b2(1,4) = 3.0d0*zi*d - zk*ri2
         b2(1,5) = 3.0d0*xi*xi - ri2
         b2(1,6) = 3.0d0*xi*yi
         b2(1,7) = 3.0d0*xi*zi
         b2(1,8) = b2(1,6)
         b2(1,9) = 3.0d0*yi*yi - ri2
         b2(1,10) = 3.0d0*yi*zi
         b2(1,11) = b2(1,7)
         b2(1,12) = b2(1,10)
         b2(1,13) = 3.0d0*zi*zi - ri2
         b2(2,1) = 3.0d0*xk*d - xi*rk2
         b2(2,2) = 3.0d0*d + xi*xk
         b2(2,3) = 3.0d0*xk*yi - 2.0d0*xi*yk
         b2(2,4) = 3.0d0*zi*xk - 2.0d0*xi*zk     
         b2(2,5) = 4.0d0*xi      
         b2(2,6) = 3.0d0*yi      
         b2(2,7) = 3.0d0*zi      
         b2(2,8) = b2(2,6)    
         b2(2,9) = -2.0d0*xi
         b2(2,11) = b2(2,7)    
         b2(2,13) = b2(2,9)    
         b2(3,1) = 3.0d0*yk*d - yi*rk2
         b2(3,2) = 3.0d0*yk*xi - 2.0d0*yi*xk
         b2(3,3) = 3.0d0*d + yi*yk
         b2(3,4) = 3.0d0*yk*zi - 2.0d0*yi*zk
         b2(3,5) = -2.0d0*yi    
         b2(3,6) = 3.0d0*xi      
         b2(3,8) = b2(3,6)    
         b2(3,9) = 4.0d0*yi 
         b2(3,10) = 3.0d0*zi  
         b2(3,12) = b2(3,10)    
         b2(3,13) = b2(3,5)   
         b2(4,1) = 3.0d0*zk*d - zi*rk2
         b2(4,2) = 3.0d0*zk*xi - 2.0d0*zi*xk
         b2(4,3) = 3.0d0*zk*yi - 2.0d0*zi*yk
         b2(4,4) = 3.0d0*d + zi*zk
         b2(4,5) = -2.0d0*zi
         b2(4,7) = 3.0d0*xi
         b2(4,9) = b2(4,5)
         b2(4,10) = 3.0d0*yi
         b2(4,11) = b2(4,7)
         b2(4,12) = b2(4,10)
         b2(4,13) = 4.0d0*zi
         b2(5,1) = 3.0d0*xk*xk - rk2
         b2(5,2) = 4.0d0*xk
         b2(5,3) = -2.0d0*yk
         b2(5,4) = -2.0d0*zk
         b2(5,5) = 4.0d0
         b2(5,9) = -2.0d0
         b2(5,13) = -2.0d0
         b2(6,1) = 3.0d0*xk*yk
         b2(6,2) = 3.0d0*yk
         b2(6,3) = 3.0d0*xk
         b2(6,6) = 3.0d0
         b2(6,8) = 3.0d0
         b2(7,1) = 3.0d0*xk*zk
         b2(7,2) = 3.0d0*zk 
         b2(7,4) = 3.0d0*xk 
         b2(7,7) = 3.0d0 
         b2(7,11) = 3.0d0 
         b2(8,1) = b2(6,1)
         b2(8,2) = b2(6,2)
         b2(8,3) = b2(6,3)
         b2(8,6) = 3.0d0
         b2(8,8) = 3.0d0
         b2(9,1) = 3.0d0*yk*yk - rk2
         b2(9,2) = -2.0d0*xk 
         b2(9,3) = 4.0d0*yk 
         b2(9,4) = -2.0d0*zk
         b2(9,5) = -2.0d0 
         b2(9,9) = 4.0d0 
         b2(9,13) = -2.0d0
         b2(10,1) = 3.0d0*yk*zk
         b2(10,3) = 3.0d0*zk
         b2(10,4) = 3.0d0*yk
         b2(10,10) = 3.0d0
         b2(10,12) = 3.0d0
         b2(11,1) = b2(7,1)
         b2(11,2) = b2(7,2)
         b2(11,4) = b2(7,4)
         b2(11,7) = 3.0d0
         b2(11,11) = 3.0d0
         b2(12,1) = b2(10,1)
         b2(12,3) = b2(10,3)
         b2(12,4) = b2(10,4)
         b2(12,10) = 3.0d0
         b2(12,12) = 3.0d0
         b2(13,1) = 3.0d0*zk*zk - rk2
         b2(13,2) = -2.0d0*xk
         b2(13,3) = -2.0d0*yk
         b2(13,4) = 4.0d0*zk
         b2(13,5) = -2.0d0
         b2(13,9) = -2.0d0
         b2(13,13) = 4.0d0
         do fi = 1, isiz
            m2t2(fi) = 0.0d0
            do fj = 1, ksiz
               m2t2(fi) = m2t2(fi) + b2(fi,fj)*rpk(fj)
            end do
         end do
         term = 0.0d0
         do fi = 1, isiz
            term = term + rpi(fi)*m2t2(fi)
         end do
         term = 3.0d0 * size * term / (3.0d0*ratio+2.0d0)
         eik = eik + term
         size = size * size2
      end if
c
c     explicit formula for the 3rd summation term
c
      if (nn .ge. 3) then
         d2 = d*d
         xi2 = xi*xi
         yi2 = yi*yi
         zi2 = zi*zi
         xk2 = xk*xk
         yk2 = yk*yk
         zk2 = zk*zk
         b1(1,1) = d*(2.5d0*d2-1.5d0*ri2*rk2)
         b1(1,2) = 7.5d0*d2*xi-3.0d0*xk*ri2*d-1.5d0*xi*ri2*rk2
         b1(1,3) = 7.5d0*d2*yi-3.0d0*yk*ri2*d-1.5d0*yi*ri2*rk2
         b1(1,4) = 7.5d0*d2*zi-3.0d0*zk*ri2*d-1.5d0*zi*ri2*rk2
         b1(1,5) = 15.0d0*d*xi2-3.0d0*ri2*(d+2.0d0*xi*xk)
         b1(1,6) = 15.d0*xi*yi*d - 3.0d0*ri2*(xi*yk+xk*yi)
         b1(1,7) = 15.d0*xi*zi*d - 3.0d0*ri2*(xi*zk+xk*zi)
         b1(1,8) = b1(1,6)
         b1(1,9) = 15.0d0*d*yi2-3.0d0*ri2*(d+2.0d0*yi*yk)
         b1(1,10) = 15.d0*yi*zi*d - 3.0d0*ri2*(yi*zk+yk*zi)
         b1(1,11) = b1(1,7)
         b1(1,12) = b1(1,10)
         b1(1,13) = 15.0d0*d*zi2-3.0d0*ri2*(d+2.0d0*zi*zk)
         b1(2,1) = 7.5d0*d2*xk-3.0d0*xi*rk2*d-1.5d0*xk*ri2*rk2
         b1(2,2) = 7.5d0*d2+9.0d0*xi*xk*d-3.0d0*xi2*rk2-3.0d0*xk2*ri2
     &                -1.5d0*ri2*rk2
         b1(2,3) = 3.0d0*((5.0d0*xk*yi-2.0d0*xi*yk)*d
     &                -xi*yi*rk2-xk*yk*ri2)
         b1(2,4) = 3.0d0*((5.0d0*xk*zi-2.0d0*xi*zk)*d
     &                -xi*zi*rk2-xk*zk*ri2)
         b1(2,5) = 24.0d0*xi*yi*yk + 24.0d0*xi*zi*zk + 18.0d0*xi2*xk
     &                - 9.0d0*xk*yi2  - 9.0d0*xk*zi2      
         b1(2,6) = (8.0d0*yi*xk*xi - 3.0d0*xi2*yk + 4.0d0*yi2*yk
     &                - yk*zi2  + 5.0d0*yi*zi*zk)*3.0d0      
         b1(2,7) = 15.0d0*zi*yi*yk + 12.0d0*zi2*zk - 9.0d0*xi2*zk
     &                - 3.0d0*zk*yi2  + 24.0d0*zi*xk*xi     
         b1(2,8) = b1(2,6)    
         b1(2,9) = - 9.0d0*xi2*xk + 12.0d0*xk*yi2  - 3.0d0*xk*zi2
     &                - 18.0d0*xi*yi*yk - 6.0d0*xi*zi*zk
         b1(2,10) = 15.0d0*zi*xk*yi - 6.0d0*zi*xi*yk - 6.0d0*yi*xi*zk  
         b1(2,11) = b1(2,7)    
         b1(2,12) = b1(2,10)    
         b1(2,13) = - 6.0d0*xi*yi*yk - 9.0d0*xi2*xk - 3.0d0*xk*yi2
     &                 + 12.0d0*xk*zi2  - 18.0d0*xi*zi*zk    
         b1(3,1) = 7.5d0*d2*yk-3.0d0*yi*rk2*d-1.5d0*yk*ri2*rk2
         b1(3,2) = 3.0d0*((5.0d0*xi*yk-2.0d0*xk*yi)*d
     &                -xi*yi*rk2-xk*yk*ri2)
         b1(3,3) = 7.5d0*d2+9.0d0*yi*yk*d-3.0d0*yi2*rk2-3.0d0*yk2*ri2
     &                -1.5d0*ri2*rk2
         b1(3,4) = 3.0d0*((5.0d0*yk*zi-2.0d0*yi*zk)*d
     &                -yi*zi*rk2-yk*zk*ri2)
         b1(3,5) = - 9.0d0*yi2*yk - 6.0d0*yi*zi*zk - 18.0d0*yi*xk*xi
     &                + 12.0d0*xi2*yk - 3.0d0*yk*zi2
         b1(3,6) = 12.0d0*xi2*xk + 15.0d0*xi*zi*zk - 9.0d0*xk*yi2
     &                - 3.0d0*xk*zi2  + 24.0d0*xi*yi*yk
         b1(3,7) = 15.0d0*zi*xi*yk - 6.0d0*yi*xi*zk - 6.0d0*zi*xk*yi     
         b1(3,8) = b1(3,6)    
         b1(3,9) = - 9.0d0*xi2*yk + 18.0d0*yi2*yk - 9.0d0*yk*zi2
     &                + 24.0d0*yi*xk*xi + 24.0d0*yi*zi*zk
         b1(3,10) = 24.0d0*zi*yi*yk - 3.0d0*xi2*zk - 9.0d0*zk*yi2
     &                + 12.0d0*zi2*zk + 15.0d0*zi*xk*xi
         b1(3,11) = b1(3,7)  
         b1(3,12) = b1(3,10)    
         b1(3,13) = - 3.0d0*xi2*yk - 9.0d0*yi2*yk + 12.0d0*yk*zi2
     &                 - 18.0d0*yi*zi*zk - 6.0d0*yi*xk*xi   
         b1(4,1) = 7.5d0*d2*zk-3.0d0*zi*rk2*d-1.5d0*zk*ri2*rk2
         b1(4,2) = 3.0d0*((5.0d0*xi*zk-2.0d0*xk*zi)*d
     &                -xi*zi*rk2-xk*zk*ri2)
         b1(4,3) = 3.0d0*((5.0d0*yi*zk-2.0d0*yk*zi)*d
     &                -yi*zi*rk2-yk*zk*ri2)
         b1(4,4) = 7.5d0*d2+9.0d0*zi*zk*d-3.0d0*zi2*rk2-3.0d0*zk2*ri2
     &                -1.5d0*ri2*rk2
         b1(4,5) = 12.0d0*xi2*zk - 3.0d0*zk*yi2 - 9.0d0*zi2*zk
     &                - 18.0d0*zi*xk*xi - 6.0d0*zi*yi*yk
         b1(4,6) = 15.0d0*yi*xi*zk - 6.0d0*zi*xi*yk - 6.0d0*zi*xk*yi
         b1(4,7) = 24.0d0*xi*zi*zk + 12.0d0*xi2*xk - 3.0d0*xk*yi2
     &                - 9.0d0*xk*zi2  + 15.0d0*xi*yi*yk
         b1(4,8) = b1(4,6)
         b1(4,9) = - 6.0d0*zi*xk*xi - 9.0d0*zi2*zk - 3.0d0*xi2*zk
     &                + 12.0d0*zk*yi2  - 18.0d0*zi*yi*yk
         b1(4,10) = 15.0d0*yi*xk*xi + 12.0d0*yi2*yk - 9.0d0*yk*zi2
     &                + 24.0d0*yi*zi*zk - 3.0d0*xi2*yk
         b1(4,11) = b1(4,7)
         b1(4,12) = b1(4,10)
         b1(4,13) = 24.0d0*zi*xk*xi + 18.0d0*zi2*zk - 9.0d0*xi2*zk
     &                - 9.0d0*zk*yi2  + 24.0d0*zi*yi*yk
         b1(5,1) = 15.0d0*d*xk2-3.0d0*rk2*(d+2.0d0*xi*xk)
         b1(5,2) = 18.0d0*xi*xk2  + 24.0d0*xk*yi*yk + 24.0d0*xk*zi*zk
     &                - 9.0d0*xi*yk2  - 9.0d0*xi*zk2 
         b1(5,3) = 12.0d0*yi*xk2 - 9.0d0*yk2*yi - 3.0d0*yi*zk2
     &                - 18.0d0*xk*xi*yk - 6.0d0*yk*zi*zk
         b1(5,4) = - 9.0d0*zk2*zi - 6.0d0*zk*yi*yk - 18.0d0*xk*xi*zk
     &                + 12.0d0*zi*xk2  - 3.0d0*zi*yk2
         b1(5,5) = 24.0d0*zi*zk + 24.0d0*yi*yk + 36.0d0*xi*xk 
         b1(5,6) = -18.0d0*xi*yk + 24.0d0*yi*xk
         b1(5,7) = -18.0d0*xi*zk + 24.0d0*zi*xk
         b1(5,8) = b1(5,6)
         b1(5,9) = -6.0d0*zi*zk - 18.0d0*yi*yk - 18.0d0*xi*xk
         b1(5,10) = -6.0d0*(yi*zk + zi*yk)
         b1(5,11) = b1(5,7)
         b1(5,12) = b1(5,10)
         b1(5,13) = -6.0d0*yi*yk - 18.0d0*xi*xk - 18.0d0*zi*zk
         b1(6,1) = 15.d0*xk*yk*d - 3.0d0*rk2*(xi*yk+xk*yi)
         b1(6,2) = -9.0d0*yi*xk2 + 12.0d0*yk2*yi - 3.0d0*yi*zk2
     &                + 24.0d0*xk*xi*yk + 15.0d0*yk*zi*zk
         b1(6,3) = 12.0d0*xi*xk2 + 15.0d0*xk*zi*zk - 9.0d0*xi*yk2
     &                - 3.0d0*xi*zk2  + 24.0d0*xk*yi*yk
         b1(6,4) = -6.0d0*xk*yi*zk - 6.0d0*yk*xi*zk + 15.0d0*zi*xk*yk
         b1(6,5) = -18.0d0*yi*xk + 24.0d0*xi*yk
         b1(6,6) = 24.0d0*yi*yk + 24.0d0*xi*xk + 15.0d0*zi*zk
         b1(6,7) = -6.0d0*yi*zk + 15.0d0*zi*yk
         b1(6,8) = b1(6,6)
         b1(6,9) = -18.0d0*xi*yk + 24.0d0*yi*xk
         b1(6,10) = -6.0d0*xi*zk + 15.0d0*zi*xk
         b1(6,11) = b1(6,7)
         b1(6,12) = b1(6,10)
         b1(6,13) = -6.0d0*yi*xk - 6.0d0*xi*yk
         b1(7,1) = 15.d0*xk*zk*d - 3.0d0*rk2*(xi*zk+xk*zi)
         b1(7,2) = 15.0d0*zk*yi*yk + 12.0d0*zk2*zi - 9.0d0*zi*xk2
     &                - 3.0d0*zi*yk2  + 24.0d0*xk*xi*zk 
         b1(7,3) = - 6.0d0*zi*xk*yk - 6.0d0*yk*xi*zk + 15.0d0*xk*yi*zk
         b1(7,4) = 12.0d0*xi*xk2  - 3.0d0*xi*yk2  - 9.0d0*xi*zk2
     &                + 15.0d0*xk*yi*yk + 24.0d0*xk*zi*zk
         b1(7,5) = -18.0d0*zi*xk + 24.0d0*xi*zk
         b1(7,6) = -6.0d0*zi*yk + 15.0d0*yi*zk
         b1(7,7) = 24.0d0*xi*xk + 24.0d0*zi*zk + 15.0d0*yi*yk 
         b1(7,8) = b1(7,6)
         b1(7,9) = -6.0d0*zi*xk - 6.0d0*xi*zk
         b1(7,10) = -6.0d0*xi*yk + 15.0d0*yi*xk 
         b1(7,11) = b1(7,7)
         b1(7,12) = b1(7,10) 
         b1(7,13) = -18.0d0*xi*zk + 24.0d0*zi*xk
         b1(9,1) = 15.0d0*d*yk2-3.0d0*rk2*(d+2.0d0*yi*yk)
         b1(9,2) = -9.0d0*xi*xk2 + 12.0d0*xi*yk2 - 3.0d0*xi*zk2
     &                - 18.0d0*xk*yi*yk - 6.0d0*xk*zi*zk
         b1(9,3) = -9.0d0*yi*xk2  + 18.0d0*yk2*yi - 9.0d0*yi*zk2
     &                + 24.0d0*yk*zi*zk + 24.0d0*xk*xi*yk
         b1(9,4) = 12.0d0*zi*yk2 - 18.0d0*zk*yi*yk - 3.0d0*zi*xk2
     &                - 9.0d0*zk2*zi - 6.0d0*xk*xi*zk
         b1(9,5) = -18.0d0*xi*xk - 6.0d0*zi*zk - 18.0d0*yi*yk
         b1(9,6) = -18.0d0*yi*xk + 24.0d0*xi*yk
         b1(9,7) = -6.0d0*zi*xk - 6.0d0*xi*zk
         b1(9,8) = b1(9,6)
         b1(9,9) = 24.0d0*xi*xk + 24.0d0*zi*zk + 36.0d0*yi*yk 
         b1(9,10) = -18.0d0*yi*zk + 24.0d0*zi*yk
         b1(9,11) = b1(9,7)
         b1(9,12) = b1(9,10)
         b1(9,13) = -18.0d0*yi*yk - 6.0d0*xi*xk - 18.0d0*zi*zk
         b1(10,1) = 15.d0*yk*zk*d - 3.0d0*rk2*(yi*zk+yk*zi)
         b1(10,2) = -6.0d0*zi*xk*yk -6.0d0*xk*yi*zk + 15.0d0*yk*xi*zk
         b1(10,3) = 12.0d0* zk2*zi + 15.0d0*xk*xi*zk - 3.0d0*zi*xk2
     &                 - 9.0d0*zi*yk2  + 24.0d0*zk*yi*yk
         b1(10,4) = 15.0d0*xk*xi*yk + 12.0d0*yk2*yi - 3.0d0*yi*xk2
     &                 - 9.0d0*yi*zk2  + 24.0d0*yk*zi*zk
         b1(10,5) = -6.0d0*yi*zk - 6.0d0*zi*yk
         b1(10,6) = -6.0d0*zi*xk + 15.0d0*xi*zk
         b1(10,7) = -6.0d0*yi*xk + 15.0d0*xi*yk
         b1(10,8) = b1(10,6)
         b1(10,9) = 24.0d0*yi*zk - 18.0d0*zi*yk
         b1(10,10) = 15.0d0*xi*xk + 24.0d0*zi*zk + 24.0d0*yi*yk
         b1(10,11) = b1(10,7)
         b1(10,12) = b1(10,10)
         b1(10,13) = -18.0d0*yi*zk + 24.0d0*zi*yk
         b1(13,1) = 15.0d0*d*zk2-3.0d0*rk2*(d+2.0d0*zi*zk)
         b1(13,2) = 12.0d0*xi*zk2 - 18.0d0*xk*zi*zk - 9.0d0*xi*xk2
     &                 - 3.0d0*xi*yk2 - 6.0d0*xk*yi*yk
         b1(13,3) = 12.0d0*yi*zk2 - 3.0d0*yi*xk2 - 9.0d0*yk2*yi
     &                 - 18.0d0*yk*zi*zk - 6.0d0*xk*xi*yk
         b1(13,4) = -9.0d0*zi*xk2 - 9.0d0*zi*yk2 + 18.0d0*zk2*zi
     &                 + 24.0d0*xk*xi*zk + 24.0d0*zk*yi*yk
         b1(13,5) = -6.0d0*yi*yk - 18.0d0*zi*zk - 18.0d0*xi*xk
         b1(13,6) = -6.0d0*yi*xk - 6.0d0*xi*yk
         b1(13,7) = 24.0d0*xi*zk - 18.0d0*zi*xk
         b1(13,8) = b1(13,6)
         b1(13,9) = -18.0d0*yi*yk - 6.0d0*xi*xk - 18.0d0*zi*zk
         b1(13,10) = 24.0d0*yi*zk - 18.0d0*zi*yk
         b1(13,11) = b1(13,7)
         b1(13,12) = b1(13,10)
         b1(13,13) = 36.0d0*zi*zk + 24.0d0*xi*xk + 24.0d0*yi*yk
         do fi = 1, isiz
            b1(8,fi) = b1(6,fi)
            b1(11,fi) = b1(7,fi)
            b1(12,fi) = b1(10,fi)
         end do
         do fi = 1, isiz
            m2t2(fi) = 0.0d0
            do fj = 1, ksiz
               m2t2(fi) = m2t2(fi) + b1(fi,fj)*rpk(fj)
            end do
         end do
         term = 0.0d0
         do fi = 1, isiz
            term = term + rpi(fi)*m2t2(fi)
         end do
         term = 4.0d0 * size * term / (4.0d0*ratio+3.0d0)
         eik = eik + term
         size = size * size2
      end if
c
c     recursive formulation of 4th through nth summation terms
c
      do fii = 4, nn
         do fi = 1, p_e1
            if (fi .eq. 8) then
               do fj = 1, p_e2
                  rklij(fi,fj) = rklij(6,fj)
               end do
            else if (fi .eq. 11) then
               do fj = 1, p_e2
                  rklij(fi,fj) = rklij(7,fj)
               end do
            else if (fi .eq. 12) then
               do fj = 1, p_e2
                  rklij(fi,fj) = rklij(10,fj)
               end do
            else
               do fj = 1, p_e2
                  if (fj .eq. 8) then
                     rklij(fi,fj) = rklij(fi,6)
                  else if (fj .eq. 11) then
                     rklij(fi,fj) = rklij(fi,7)
                  else if (fj .eq. 12) then
                     rklij(fi,fj) = rklij(fi,10)
                  else
                     rklij(fi,fj) = d1d2 (fii,xi,yi,zi,xk,yk,zk,
     &                                    d,ri2,rk2,ind1_x(fi),
     &                                    ind1_y(fi),ind1_z(fi),
     &                              ind2_x(fj),ind2_y(fj),ind2_z(fj))
                  end if
               end do
            end if
         end do
c
c     update storage of the last two sets of matrix elements
c
         do fi = 1, p_e1
           do fj = 1, p_e2
              b2(fj,fi) = b1(fj,fi)
              b1(fj,fi) = rklij(fj,fi)
           end do
         end do
c
c     compute interaction energy between the two multipole sites
c
         do fi = 1, isiz
            m2t2(fi) = 0.0d0
            do fj = 1, ksiz
               m2t2(fi) = m2t2(fi) + rklij(fi,fj)*rpk(fj)
            end do
         end do
         term = 0.0d0
         do fi = 1, isiz
            term = term + rpi(fi)*m2t2(fi)
         end do
         term = term * size * dble(fii+1)
     &             / (dble(fii+1)*ratio+dble(fii))
         eik = eik + term
         size = size * size2
      end do
      eik = factor * eik
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine rindex  ##
c     ##                     ##
c     #########################
c
c
      subroutine rindex (n,m,ind_x,ind_y,ind_z,p_s,p_e)
      implicit none
      integer i,j,k,n,m,p_s,p_e
      integer ind_x(m),ind_y(m),ind_z(m)
c
c
      p_s = 1
      p_e = 1
c
c
      do i = 1, m
         ind_x(i) = 0
         ind_y(i) = 0
         ind_z(i) = 0
      end do
c
c
      k = 1
      do i = 1, n
         do j = p_s, p_e
            k = k + 1  
            ind_x(k) = ind_x(j) + 1
            ind_y(k) = ind_y(j) 
            ind_z(k) = ind_z(j)
         end do
         do j = p_s, p_e
            k = k + 1  
            ind_x(k) = ind_x(j) 
            ind_y(k) = ind_y(j) + 1
            ind_z(k) = ind_z(j) 
         end do
         do j = p_s, p_e
            k = k + 1  
            ind_x(k) = ind_x(j) 
            ind_y(k) = ind_y(j) 
            ind_z(k) = ind_z(j) + 1
         end do
         p_s = p_e + 1
         p_e = k
      end do
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine ijk_pt  ##
c     ##                     ##
c     #########################
c
c
      subroutine ijk_pt
      implicit none
      include 'rxnfld.i'
      integer i,j,k
c
c
      do i = 0, 5
         do j = 0, 5
            do k = 0, 5
               ijk(i,j,k) = (3**(i+j+k) + 3**(j+k) + 3**k - 1) / 2
            end do
         end do
      end do
      return
      end
c
c
c     #####################
c     ##                 ##
c     ##  function d1d2  ##
c     ##                 ##
c     #####################
c
c
      function d1d2 (n,x1,y1,z1,x2,y2,z2,d,r1sq,r2sq,i,j,k,s,t,u)
      implicit none
      include 'rxnfld.i'
      integer n,i,j,k,s,t,u
      integer is,it,iu,js,jt,ju,ks,kt,ku
      real*8 x1,y1,z1,x2,y2,z2
      real*8 d1d2,d,r1sq,r2sq,f,g
c
c
      if (n.lt.i+j+k  .or. n.lt.s+t+u) then
         d1d2 = 0.0d0
         return
      end if
c
c
      is = i*s
      it = i*t
      iu = i*u
      js = j*s
      jt = j*t
      ju = j*u
      ks = k*s
      kt = k*t
      ku = k*u
c
c
      f = d*b1(ijk(i,j,k),ijk(s,t,u))
      g = r1sq*r2sq*b2(ijk(i,j,k),ijk(s,t,u))
c
c
      if (i .ne. 0) then
         f = f + i*x2*b1(ijk(i-1,j,k),ijk(s,t,u))
         g = g + 2.0d0*i*x1*r2sq*b2(ijk(i-1,j,k),ijk(s,t,u))
         if (i .ne. 1)  g = g + i*(i-1)*r2sq*b2(ijk(i-2,j,k),ijk(s,t,u))
      end if
      if (j .ne. 0) then
         f = f + j*y2*b1(ijk(i,j-1,k),ijk(s,t,u))
         g = g + 2.0d0*j*y1*r2sq*b2(ijk(i,j-1,k),ijk(s,t,u))
         if (j .ne. 1)  g = g + j*(j-1)*r2sq*b2(ijk(i,j-2,k),ijk(s,t,u))
      end if
      if (k .ne. 0) then
         f = f + k*z2*b1(ijk(i,j,k-1),ijk(s,t,u))
         g = g + 2.0d0*k*z1*r2sq*b2(ijk(i,j,k-1),ijk(s,t,u))
         if (k .ne. 1)  g = g + k*(k-1)*r2sq*b2(ijk(i,j,k-2),ijk(s,t,u))
      end if
      if (s .ne. 0) then
         f = f + s*x1*b1(ijk(i,j,k),ijk(s-1,t,u))
         g = g + 2.0d0*s*x2*r1sq*b2(ijk(i,j,k),ijk(s-1,t,u))
         if (s .ne. 1)  g = g + s*(s-1)*r1sq*b2(ijk(i,j,k),ijk(s-2,t,u))
      end if
      if (t .ne. 0) then
         f = f + t*y1*b1(ijk(i,j,k),ijk(s,t-1,u))
         g = g + 2.0d0*t*y2*r1sq*b2(ijk(i,j,k),ijk(s,t-1,u))
         if (t .ne. 1)  g = g + t*(t-1)*r1sq*b2(ijk(i,j,k),ijk(s,t-2,u))
      end if
      if (u .ne. 0) then
         f = f + u*z1*b1(ijk(i,j,k),ijk(s,t,u-1))
         g = g + 2.0d0*u*z2*r1sq*b2(ijk(i,j,k),ijk(s,t,u-1))
         if (u .ne. 1)  g = g + u*(u-1)*r1sq*b2(ijk(i,j,k),ijk(s,t,u-2))
      end if
c
c
      if (is .ne. 0) then
         f = f + is*b1(ijk(i-1,j,k),ijk(s-1,t,u))
         g = g + 4.0d0*is*x1*x2*b2(ijk(i-1,j,k),ijk(s-1,t,u))
         if (i .ne. 1) then
            g = g + 2.0d0*(i-1)*is*x2*b2(ijk(i-2,j,k),ijk(s-1,t,u))
            if (s .ne. 1)
     &         g = g + (i-1)*(s-1)*is*b2(ijk(i-2,j,k),ijk(s-2,t,u))
         end if
         if (s .ne. 1)
     &      g = g + 2.0d0*(s-1)*is*x1*b2(ijk(i-1,j,k),ijk(s-2,t,u))
      end if
      if (jt .ne. 0) then
         f = f + jt*b1(ijk(i,j-1,k),ijk(s,t-1,u))
         g = g + 4.0d0*jt*y1*y2*b2(ijk(i,j-1,k),ijk(s,t-1,u))
         if (j .ne. 1) then
            g = g + 2.0d0*(j-1)*jt*y2*b2(ijk(i,j-2,k),ijk(s,t-1,u))
            if (t .ne. 1)
     &         g = g + (j-1)*(t-1)*jt*b2(ijk(i,j-2,k),ijk(s,t-2,u))
         end if
         if (t .ne. 1)
     &      g = g + 2.0d0*(t-1)*jt*y1*b2(ijk(i,j-1,k),ijk(s,t-2,u))
      end if
      if (ku .ne. 0) then
         f = f + ku*b1(ijk(i,j,k-1),ijk(s,t,u-1))
         g = g + 4.0d0*ku*z1*z2*b2(ijk(i,j,k-1),ijk(s,t,u-1))
         if (k .ne. 1) then
            g = g + 2.0d0*(k-1)*ku*z2*b2(ijk(i,j,k-2),ijk(s,t,u-1))
            if (u .ne. 1)
     &         g = g + (k-1)*(u-1)*ku*b2(ijk(i,j,k-2),ijk(s,t,u-2))
         end if
         if (u .ne. 1)
     &      g = g + 2.0d0*(u-1)*ku*z1*b2(ijk(i,j,k-1),ijk(s,t,u-2))
      end if
      if (it .ne. 0) then
         g = g + 4.0d0*it*x1*y2*b2(ijk(i-1,j,k),ijk(s,t-1,u))
         if (i .ne. 1) then
            g = g + 2.0d0*(i-1)*it*y2*b2(ijk(i-2,j,k),ijk(s,t-1,u))
            if (t .ne. 1)
     &         g = g + (i-1)*(t-1)*it*b2(ijk(i-2,j,k),ijk(s,t-2,u))
         end if
         if (t .ne. 1)
     &      g = g + 2.0d0*(t-1)*it*x1*b2(ijk(i-1,j,k),ijk(s,t-2,u))
      end if
      if (iu .ne. 0) then
         g = g + 4.0d0*iu*x1*z2*b2(ijk(i-1,j,k),ijk(s,t,u-1))
         if (i .ne. 1) then
            g = g + 2.0d0*(i-1)*iu*z2*b2(ijk(i-2,j,k),ijk(s,t,u-1))
            if (u .ne. 1)
     &         g = g + (i-1)*(u-1)*iu*b2(ijk(i-2,j,k),ijk(s,t,u-2))
         end if
         if (u .ne. 1)
     &      g = g + 2.0d0*(u-1)*iu*x1*b2(ijk(i-1,j,k),ijk(s,t,u-2))
      end if
      if (js .ne. 0) then
         g = g + 4.0d0*js*y1*x2*b2(ijk(i,j-1,k),ijk(s-1,t,u))
         if (j .ne. 1) then
            g = g + 2.0d0*(j-1)*js*x2*b2(ijk(i,j-2,k),ijk(s-1,t,u))
            if (s .ne. 1)
     &         g = g + (j-1)*(s-1)*js*b2(ijk(i,j-2,k),ijk(s-2,t,u))
         end if
         if (s .ne. 1)
     &      g = g + 2.0d0*(s-1)*js*y1*b2(ijk(i,j-1,k),ijk(s-2,t,u))
      end if
      if (ju .ne. 0) then
         g = g + 4.0d0*ju*y1*z2*b2(ijk(i,j-1,k),ijk(s,t,u-1))
         if (j .ne. 1) then
            g = g + 2.0d0*(j-1)*ju*z2*b2(ijk(i,j-2,k),ijk(s,t,u-1))
            if (u .ne. 1)
     &         g = g + (j-1)*(u-1)*ju*b2(ijk(i,j-2,k),ijk(s,t,u-2))
         end if
         if (u .ne. 1)
     &      g = g + 2.0d0*(u-1)*ju*y1*b2(ijk(i,j-1,k),ijk(s,t,u-2))
      end if
      if (ks .ne. 0) then
         g = g + 4.0d0*ks*z1*x2*b2(ijk(i,j,k-1),ijk(s-1,t,u))
         if (k .ne. 1) then
            g = g + 2.0d0*(k-1)*ks*x2*b2(ijk(i,j,k-2),ijk(s-1,t,u))
            if (s .ne. 1)
     &         g = g + (k-1)*(s-1)*ks*b2(ijk(i,j,k-2),ijk(s-2,t,u))
         end if
         if (s .ne. 1)
     &      g = g + 2.0d0*(s-1)*ks*z1*b2(ijk(i,j,k-1),ijk(s-2,t,u))
      end if
      if (kt .ne. 0) then
         g = g + 4.0d0*kt*z1*y2*b2(ijk(i,j,k-1),ijk(s,t-1,u))
         if (k .ne. 1) then
            g = g + 2.0d0*(k-1)*kt*y2*b2(ijk(i,j,k-2),ijk(s,t-1,u))
            if (t .ne. 1)
     &         g = g + (k-1)*(t-1)*kt*b2(ijk(i,j,k-2),ijk(s,t-2,u))
         end if
         if (t .ne. 1)
     &      g = g + 2.0d0*(t-1)*kt*z1*b2(ijk(i,j,k-1),ijk(s,t-2,u))
      end if
c
c
      f = dble(2*n-1) * f
      g = dble(n-1) * g
      d1d2 = (f-g) / dble(n)
      return
      end
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1996 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine erxnfld1  --  reaction field energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "erxnfld1" calculates the macroscopic reaction field energy
c     and derivatives with respect to Cartesian coordinates
c
c
      subroutine erxnfld1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'energi.i'
      integer i,j
c
c
c     zero out macroscopic reaction field energy and derivatives
c
      er = 0.0d0
      do i = 1, n
         do j = 1, 3
            der(j,i) = 0.0d0
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1996 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erxnfld2  --  atom-wise reaction field Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erxnfld2" calculates second derivatives of the macroscopic
c     reaction field energy for a single atom at a time
c
c
      subroutine erxnfld2 (i)
      implicit none
      integer i
c
c
      return
      end
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1996 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erxnfld3  --  reaction field energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erxnfld3" calculates the macroscopic reaction field energy,
c     and also partitions the energy among the atoms
c
c
      subroutine erxnfld3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'shunt.i'
      include 'usage.i'
      integer i,j,k
      integer ii,iz,ix,kk,kz,kx
      real*8 eik,xr,yr,zr,r2,r,di,dk
      real*8 a(3,3),rpi(13),rpk(13)
      logical header,huge,iuse,kuse
c
c
c     zero out the reaction field energy and partitioning
c
      ner = 0
      er = 0.0d0
      do i = 1, n
         aer(i) = 0.0d0
      end do
      header = .true.
c
c     set the switching function coefficients
c
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
c
      call ijk_pt
c
c     calculate the reaction field interaction energy term
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, polsiz(ii)
            rpi(j) = rpole(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
            if (iuse .or. kuse) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, polsiz(kk)
                     rpk(j) = rpole(j,kk)
                  end do
                  call erfik (ii,kk,i,k,rpi,rpk,eik)
                  ner = ner + 1
                  er = er + eik
                  aer(i) = aer(i) + 0.5d0*eik
                  aer(k) = aer(k) + 0.5d0*eik
c
c     print a warning if the energy of this interaction is large
c
                  huge = (eik .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Reaction Field',
     &                             ' Interactions :',
     &                          //,' Type',11x,'Atom Names',
     &                             9x,'Dist from Origin',4x,'R(1-2)',
     &                             6x,'Energy',/)
                     end if
                     r = sqrt(r2)
                     di = sqrt(x(i)**2 + y(i)**2 + z(i)**2)
                     dk = sqrt(x(k)**2 + y(k)**2 + z(k)**2)
                     write (iout,20)  i,name(i),k,name(k),di,dk,r,eik
   20                format (' RxnFld   ',i5,'-',a3,1x,i5,'-',
     &                          a3,2x,3f10.4,f12.4)
                  end if
               end if
            end if
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
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine esolv  --  macroscopic solvation energy  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "esolv" calculates the macroscopic solvation energy via
c     either the Eisenberg-McLachlan ASP, Ooi-Scheraga SASA or
c     Macromodel GB/SA solvation model
c
c
      subroutine esolv
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'potent.i'
      include 'solute.i'
      real*8 earea,epol,probe
      real*8 aarea(maxatm),darea(3,maxatm)
c
c
c     zero out the macroscopic solvation energy
c
      es = 0.0d0
c
c     compute the surface area-based solvation energy term
c
      probe = 1.4d0
      call surface (earea,aarea,darea,rsolv,vsolv,probe)
      es = es + earea
c
c     get the generalized Born term for GB/SA solvation
c
      if (use_gbsa) then
         call born
         call egbsa (epol)
         es = es + epol
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine egbsa  --  generalized Born energy for GB/SA  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "egbsa" calculates the generalized Born energy term
c     for the Macromodel GB/SA solvation model
c
c
      subroutine egbsa (epol)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'shunt.i'
      include 'solute.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,ii,kk
      real*8 epol,e,f,fi,fik,fgb,rb2,rm2
      real*8 xi,yi,zi,xr,yr,zr
      real*8 dwater,shift,taper,trans
      real*8 r,r2,r3,r4,r5,r6,r7
      logical iuse
c
c
c     zero out the GB/SA polarization energy
c
      epol = 0.0d0
      if (nion .eq. 0)  return
c
c     set the dielectric constant and energy conversion factor
c
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set cutoff distances and switching function coefficients
c
      call switch ('CHARGE')
c
c     calculate GB/SA electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            if (iuse .or. use(k)) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  fik = fi * pchg(kk)
                  rb2 = rborn(i) * rborn(k)
                  fgb = sqrt(r2 + rb2*exp(-0.25d0*r2/rb2))
                  e = fik / fgb
c
c     use shifted energy switching if near the cutoff distance
c
                  rm2 = (0.5d0 * (off+cut))**2
                  shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  e = e - shift
                  if (r2 .gt. cut2) then
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     e = e * taper + trans
                  end if
                  if (i .eq. k)  e = 0.5d0 * e
                  epol = epol + e
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
c     calculate GB/SA polarization energy with other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            if (iuse .or. use(k)) then
               do j = 1, ncell
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     fik = fi * pchg(kk)
                     rb2 = rborn(i) * rborn(k)
                     fgb = sqrt(r2 + rb2*exp(-0.25d0*r2/rb2))
                     e = fik / fgb
c
c     use shifted energy switching if near the cutoff distance
c
                     rm2 = (0.5d0 * (off+cut))**2
                     shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                     e = e - shift
                     if (r2 .gt. cut2) then
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                        e = e * taper + trans
                     end if
                     if (i .eq. k)  e = 0.5d0 * e
                     epol = epol + e
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
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine esolv1  --  solvation energy and derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "esolv1" calculates the macroscopic solvation energy and
c     first derivatives with respect to Cartesian coordinates
c     using either the Eisenberg-McLachlan ASP, Ooi-Scheraga SASA
c     or Macromodel GB/SA solvation model
c
c
      subroutine esolv1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'energi.i'
      include 'potent.i'
      include 'solute.i'
      integer i,j
      real*8 earea,epol,probe
      real*8 aarea(maxatm),darea(3,maxatm)
      real*8 depol(3,maxatm)
c
c
c     zero out the macroscopic solvation energy and derivatives
c
      es = 0.0d0
      do i = 1, n
         des(1,i) = 0.0d0
         des(2,i) = 0.0d0
         des(3,i) = 0.0d0
      end do
c
c     compute the surface area-based solvation energy term
c
      probe = 1.4d0
      call surface (earea,aarea,darea,rsolv,vsolv,probe)
      es = es + earea
      do i = 1, n
         do j = 1, 3
            des(j,i) = des(j,i) + darea(j,i)
         end do
      end do
c
c     get the generalized Born term for GB/SA solvation
c
      if (use_gbsa) then
         call born
         call egbsa1 (epol,depol)
         es = es + epol
         do i = 1, n
            do j = 1, 3
               des(j,i) = des(j,i) + depol(j,i)
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egbsa1  --  generalized Born energy for GB/SA  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egbsa1" calculates the generalized Born energy and first
c     derivatives with respect to Cartesian coordinates for the
c     Macromodel GB/SA solvation model
c
c
      subroutine egbsa1 (epol,depol)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'solute.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,ii,kk
      real*8 epol,e,f,fi,fik,fgb,rb2,rm2
      real*8 xi,yi,zi,xr,yr,zr
      real*8 dedx,dedy,dedz,de,dwater
      real*8 shift,taper,dtaper,trans,dtrans
      real*8 r,r2,r3,r4,r5,r6,r7,expterm
      real*8 depol(3,maxatm)
      logical iuse
c
c
c     zero out the GB/SA polarization energy and derivatives
c
      epol = 0.0d0
      do i = 1, n
         depol(1,i) = 0.0d0
         depol(2,i) = 0.0d0
         depol(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     set the dielectric constant and energy conversion factor
c
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set cutoff distances and switching function coefficients
c
      call switch ('CHARGE')
c
c     calculate GB/SA electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            if (iuse .or. use(k)) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
                  rb2 = rborn(i) * rborn(k)
                  expterm = exp(-0.25d0*r2/rb2)
                  fgb = sqrt(r2 + rb2*expterm)
                  e = fik / fgb
                  de = -fik * (r-0.25d0*r*expterm) / fgb**3
c
c     use energy switching if near the cutoff distance
c
                  rm2 = (0.5d0 * (off+cut))**2
                  shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  e = e - shift
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
                     de = e*dtaper + de*taper + dtrans
                     e = e * taper + trans
                  end if
c
c     increment the overall energy and derivative expressions
c
                  if (i .eq. k) then
                     e = 0.5d0 * e
                     epol = epol + e
                  else
                     epol = epol + e
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     depol(1,i) = depol(1,i) + dedx
                     depol(2,i) = depol(2,i) + dedy
                     depol(3,i) = depol(3,i) + dedz
                     depol(1,k) = depol(1,k) - dedx
                     depol(2,k) = depol(2,k) - dedy
                     depol(3,k) = depol(3,k) - dedz
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*dedx
                        viry = viry + yr*dedy
                        virz = virz + zr*dedz
                     end if
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
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
c     calculate GB/SA polarization energy with other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            if (iuse .or. use(k)) then
               do j = 1, ncell
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     fik = fi * pchg(kk)
                     rb2 = rborn(i) * rborn(k)
                     expterm = exp(-0.25d0*r2/rb2)
                     fgb = sqrt(r2 + rb2*expterm)
                     e = fik / fgb
                     de = -fik * (r-0.25d0*r*expterm) / fgb**3
c
c     use shifted energy switching if near the cutoff distance
c
                     rm2 = (0.5d0 * (off+cut))**2
                     shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                     e = e - shift
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
                        de = e*dtaper + de*taper + dtrans
                        e = e * taper + trans
                     end if
c
c     form the chain rule terms for derivative expressions
c
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
                     if (i .eq. k)  e = 0.5d0 * e
                     epol = epol + e
                     depol(1,i) = depol(1,i) + dedx
                     depol(2,i) = depol(2,i) + dedy
                     depol(3,i) = depol(3,i) + dedz
                     if (i .ne. k) then
                        depol(1,k) = depol(1,k) - dedx
                        depol(2,k) = depol(2,k) - dedy
                        depol(3,k) = depol(3,k) - dedz
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
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine esolv2  --  atom-by-atom solvation Hessian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "esolv2" calculates second derivatives of the macroscopic
c     solvation energy using either the Eisenberg-McLachlan ASP,
c     Ooi-Scheraga SASA or Macromodel GB/SA solvation model
c
c
      subroutine esolv2 (i)
      implicit none
      include 'sizes.i'
      include 'potent.i'
      include 'solute.i'
      integer i
c     real*8 probe
c     real*8 earea,aarea(maxatm),darea(3,maxatm)
c
c
c     compute the surface area-based solvation energy term
c
c     probe = 1.4d0
c     call surface (earea,aarea,darea,rsolv,vsolv,probe)
c
c     get the generalized Born term for GB/SA solvation
c
      if (use_gbsa)  call egbsa2 (i)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egbsa2  --  atom-wise generalized Born Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egbsa2" calculates second derivatives of the generalized
c     Born energy term for the Macromodel GB/SA solvation model
c
c
      subroutine egbsa2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'hessn.i'
      include 'shunt.i'
      include 'solute.i'
      include 'units.i'
      integer i,j,k,kk,jcell
      real*8 fi,fik,e,de,d2e,d2edx,d2edy,d2edz
      real*8 xi,yi,zi,xr,yr,zr,term(3,3)
      real*8 shift,taper,dtaper,d2taper
      real*8 trans,dtrans,d2trans
      real*8 dwater,rb2,rm2,fgb,expterm
      real*8 r,r2,r3,r4,r5,r6,r7
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            dwater = 78.3d0
            fi = -electric * (1.0d0 - 1.0d0/dwater) * pchg(k)
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
c
c     set cutoff distances and switching function coefficients
c
      call switch ('CHARGE')
c
c     calculate GB/SA polarization energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         if (i .ne. k) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            if (use_image)  call image (xr,yr,zr,0)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               fik = fi * pchg(kk)
c
c     compute chain rule terms for Hessian matrix elements
c
               rb2 = rborn(i) * rborn(k)
               expterm = exp(-0.25d0*r2/rb2)
               fgb = sqrt(r2 + rb2*expterm)
               de = -fik * (r-0.25d0*r*expterm) / fgb**3
               d2e = 1.5d0 * fik * (r-0.25d0*r*expterm)**2 / fgb**5
     &                  - fik * (1.0d0-0.25d0*expterm
     &                           +0.125d0*r2*expterm/rb2) / fgb**3
c
c     use energy switching if near the cutoff distance
c
               if (r2 .gt. cut2) then
                  e = fik / fgb
                  rm2 = (0.5d0 * (off+cut))**2
                  shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  e = e - shift
                  r3 = r2 * r
                  r4 = r2 * r2
                  r5 = r2 * r3
                  r6 = r3 * r3
                  r7 = r3 * r4
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                        + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                  d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                         + 6.0d0*c3*r + 2.0d0*c2
                  trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                            + f3*r3 + f2*r2 + f1*r + f0)
                  dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                            + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                            + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                  d2trans = fik * (42.0d0*f7*r5 + 30.0d0*f6*r4
     &                             + 20.0d0*f5*r3 + 12.0d0*f4*r2
     &                             + 6.0d0*f3*r + 2.0d0*f2)
                  d2e = e*d2taper + 2.0d0*de*dtaper
     &                     + d2e*taper + d2trans
                  de = e*dtaper + de*taper + dtrans
               end if
c
c     now, form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
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
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do kk = 1, nion
         k = iion(kk)
         do jcell = 1, ncell
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call image (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               fik = fi * pchg(kk)
c
c     compute chain rule terms for Hessian matrix elements
c
               rb2 = rborn(i) * rborn(k)
               expterm = exp(-0.25d0*r2/rb2)
               fgb = sqrt(r2 + rb2*expterm)
               de = -fik * (r-0.25d0*r*expterm) / fgb**3
               d2e = 1.5d0 * fik * (r-0.25d0*r*expterm)**2 / fgb**5
     &                  - fik * (1.0d0-0.25d0*expterm
     &                           +0.125d0*r2*expterm/rb2) / fgb**3
c
c     use energy switching if near the cutoff distance
c
               if (r2 .gt. cut2) then
                  e = fik / fgb
                  rm2 = (0.5d0 * (off+cut))**2
                  shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  e = e - shift
                  r3 = r2 * r
                  r4 = r2 * r2
                  r5 = r2 * r3
                  r6 = r3 * r3
                  r7 = r3 * r4
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                        + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                  d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                         + 6.0d0*c3*r + 2.0d0*c2
                  trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                            + f3*r3 + f2*r2 + f1*r + f0)
                  dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                            + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                            + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                  d2trans = fik * (42.0d0*f7*r5 + 30.0d0*f6*r4
     &                             + 20.0d0*f5*r3 + 12.0d0*f4*r2
     &                             + 6.0d0*f3*r + 2.0d0*f2)
                  d2e = e*d2taper + 2.0d0*de*dtaper
     &                     + d2e*taper + d2trans
                  de = e*dtaper + de*taper + dtrans
               end if
c
c     now, form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
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
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
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
c     ############################################################
c     ##                                                        ##
c     ##  subroutine esolv3  --  solvation energy and analysis  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "esolv3" calculates the macroscopic solvation energy using
c     either the Eisenberg-McLachlan ASP, Ooi-Scheraga SASA or
c     Macromodel GB/SA solvation model; also partitions the
c     energy among the atoms
c
c
      subroutine esolv3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'energi.i'
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      include 'solute.i'
      integer i
      real*8 earea,epol,probe
      real*8 aarea(maxatm),darea(3,maxatm)
      real*8 aepol(maxatm)
      logical header,huge
c
c
c     zero out the macroscopic solvation energy and partitioning
c
      nes = 0
      es = 0.0d0
      do i = 1, n
         aes(i) = 0.0d0
      end do
c
c     compute the surface area-based solvation energy term
c
      header = .true.
      probe = 1.4d0
      call surface (earea,aarea,darea,rsolv,vsolv,probe)
      nes = nes + n
      es = es + earea
      do i = 1, n
         aes(i) = aes(i) + aarea(i)
      end do
c
c     get the generalized Born term for GB/SA solvation
c
      if (use_gbsa) then
         call born
         call egbsa3 (epol,aepol)
         es = es + epol
         do i = 1, n
            aes(i) = aes(i) + aepol(i)
         end do
      end if
c
c     print a warning if the energy of any atom is large
c
      do i = 1, n
         huge = (abs(aes(i)) .gt. 25.0d0)
         if (debug .or. (verbose.and.huge)) then
            if (header) then
               header = .false.
               write (iout,10)
   10          format (/,' Individual Atomic Solvation Energy',
     &                    ' Terms :',
     &                 //,' Type',11x,'Atom Name',42x,'Energy',/)
            end if
            write (iout,20)  i,name(i),aes(i)
   20       format (' Solvate',7x,i5,'-',a3,37x,f12.4)
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egbsa3  --  generalized Born energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egbsa3" calculates the generalized Born energy term for
c     the Macromodel GB/SA solvation model; also partitions the
c     energy among the atoms
c
c
      subroutine egbsa3 (epol,aepol)
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'shunt.i'
      include 'solute.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,ii,kk
      real*8 epol,e,f,fi,fik,fgb,rb2,rm2
      real*8 xi,yi,zi,xr,yr,zr
      real*8 dwater,shift,taper,trans
      real*8 r,r2,r3,r4,r5,r6,r7
      real*8 aepol(maxatm)
      logical iuse
c
c
c     zero out the GB/SA polarization energy and partitioning
c
      epol = 0.0d0
      do i = 1, n
         aepol(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     set the dielectric constant and energy conversion factor
c
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set cutoff distances and switching function coefficients
c
      call switch ('CHARGE')
c
c     calculate GB/SA electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            if (iuse .or. use(k)) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  fik = fi * pchg(kk)
                  rb2 = rborn(i) * rborn(k)
                  fgb = sqrt(r2 + rb2*exp(-0.25d0*r2/rb2))
                  e = fik / fgb
c
c     use shifted energy switching if near the cutoff distance
c
                  rm2 = (0.5d0 * (off+cut))**2
                  shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  e = e - shift
                  if (r2 .gt. cut2) then
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     e = e * taper + trans
                  end if
                  nes = nes + 1
                  if (i .eq. k) then
                     epol = epol + 0.5d0*e
                     aepol(i) = aepol(i) + 0.5d0*e
                  else
                     epol = epol + e
                     aepol(i) = aepol(i) + 0.5d0*e
                     aepol(k) = aepol(k) + 0.5d0*e
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
c     calculate GB/SA polarization energy with other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         iuse = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            if (iuse .or. use(k)) then
               do j = 1, ncell
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     fik = fi * pchg(kk)
                     rb2 = rborn(i) * rborn(k)
                     fgb = sqrt(r2 + rb2*exp(-0.25d0*r2/rb2))
                     e = fik / fgb
c
c     use shifted energy switching if near the cutoff distance
c
                     rm2 = (0.5d0 * (off+cut))**2
                     shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                     e = e - shift
                     if (r2 .gt. cut2) then
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                        e = e * taper + trans
                     end if
                     nes = nes + 1
                     if (i .eq. k) then
                        epol = epol + 0.5d0*e
                        aepol(i) = aepol(i) + 0.5d0*e
                     else
                        epol = epol + e
                        aepol(i) = aepol(i) + 0.5d0*e
                        aepol(k) = aepol(k) + 0.5d0*e
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
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine estrbnd  --  stretch-bend cross term energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "estrbnd" calculates the stretch-bend potential energy
c
c
      subroutine estrbnd
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'bond.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'strbnd.i'
      include 'usage.i'
      integer i,j,k,istrbnd
      integer ia,ib,ic
      real*8 e,dr,dt,fgrp
      real*8 angle,force
      real*8 dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 rab,xab,yab,zab
      real*8 rcb,xcb,ycb,zcb
      logical proceed
c
c
c     zero out the stretch-bend cross term energy
c
      eba = 0.0d0
c
c     calculate the stretch-bend energy term
c
      do istrbnd = 1, nstrbnd
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         force = ksb(istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,3,ia,ib,ic,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
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
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            if (rab*rcb .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt = angle - anat(i)
c
c     get the stretch-bend interaction energy
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               dr = rab - bl(j) + rcb - bl(k)
               e = stbnunit * force * dt * dr
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total stretch-bend energy
c
               eba = eba + e
            end if
         end if
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrbnd1   --  stretch-bend energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrbnd1" calculates the stretch-bend potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine estrbnd1
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bond.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'strbnd.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,istrbnd
      integer ia,ib,ic
      real*8 e,dt,dr,fgrp
      real*8 angle,force
      real*8 dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab,rab
      real*8 xcb,ycb,zcb,rcb
      real*8 xp,yp,zp,rp
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxic,ddtdyic,ddtdzic
      real*8 ddrdxia,ddrdyia,ddrdzia
      real*8 ddrdxic,ddrdyic,ddrdzic
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 term,terma,termc
      logical proceed
c
c
c     zero out the energy and first derivative components
c
      eba = 0.0d0
      do i = 1, n
         deba(1,i) = 0.0d0
         deba(2,i) = 0.0d0
         deba(3,i) = 0.0d0
      end do
c
c     calculate the stretch-bend energy and first derivatives
c
      do istrbnd = 1, nstrbnd
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         force = ksb(istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,3,ia,ib,ic,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
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
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            if (rp .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     find chain rule terms for the bond angle deviation
c
               dt = angle - anat(i)
               terma = -radian / (rab*rab*rp)
               termc = radian / (rcb*rcb*rp)
               ddtdxia = terma * (yab*zp-zab*yp)
               ddtdyia = terma * (zab*xp-xab*zp)
               ddtdzia = terma * (xab*yp-yab*xp)
               ddtdxic = termc * (ycb*zp-zcb*yp)
               ddtdyic = termc * (zcb*xp-xcb*zp)
               ddtdzic = termc * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the bond length deviations
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               term = stbnunit * force
               dr = rab - bl(j) + rcb - bl(k)
               ddrdxia = xab / rab
               ddrdyia = yab / rab
               ddrdzia = zab / rab
               ddrdxic = xcb / rcb
               ddrdyic = ycb / rcb
               ddrdzic = zcb / rcb
c
c     scale the interaction based on its group membership
c
               if (use_group)  term = term * fgrp
c
c     get the energy and master chain rule terms for derivatives
c
               e = term * dt * dr
               dedxia = term * (dt*ddrdxia+ddtdxia*dr)
               dedyia = term * (dt*ddrdyia+ddtdyia*dr)
               dedzia = term * (dt*ddrdzia+ddtdzia*dr)
               dedxic = term * (dt*ddrdxic+ddtdxic*dr)
               dedyic = term * (dt*ddrdyic+ddtdyic*dr)
               dedzic = term * (dt*ddrdzic+ddtdzic*dr)
               dedxib = -dedxia - dedxic
               dedyib = -dedyia - dedyic
               dedzib = -dedzia - dedzic
c
c     increment the total stretch-bend energy and derivatives
c
               eba = eba + e
               deba(1,ia) = deba(1,ia) + dedxia
               deba(2,ia) = deba(2,ia) + dedyia
               deba(3,ia) = deba(3,ia) + dedzia
               deba(1,ib) = deba(1,ib) + dedxib
               deba(2,ib) = deba(2,ib) + dedyib
               deba(3,ib) = deba(3,ib) + dedzib
               deba(1,ic) = deba(1,ic) + dedxic
               deba(2,ic) = deba(2,ic) + dedyic
               deba(3,ic) = deba(3,ic) + dedzic
c
c     increment the virial for use in pressure computation
c
               if (isobaric) then
                  virx = virx + xab*dedxia + xcb*dedxic
                  viry = viry + yab*dedyia + ycb*dedyic
                  virz = virz + zab*dedzia + zcb*dedzic
               end if
            end if
         end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine estrbnd2  --  stretch-bend Hessian; analytical  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "estrbnd2" calculates the stretch-bend potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine estrbnd2 (iatom)
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'bond.i'
      include 'group.i'
      include 'hessn.i'
      include 'math.i'
      include 'strbnd.i'
      integer i,j,k,iatom
      integer ia,ib,ic,istrbnd
      real*8 angle,force
      real*8 dot,cosine
      real*8 dt,dr,fgrp
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab,rab
      real*8 xcb,ycb,zcb,rcb
      real*8 xp,yp,zp,rp,rp2
      real*8 term,terma,termc
      real*8 xrab,yrab,zrab,rab2
      real*8 xrcb,yrcb,zrcb,rcb2
      real*8 xabp,yabp,zabp
      real*8 xcbp,ycbp,zcbp
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxib,ddtdyib,ddtdzib
      real*8 ddtdxic,ddtdyic,ddtdzic
      real*8 ddrdxia,ddrdyia,ddrdzia
      real*8 ddrdxib,ddrdyib,ddrdzib
      real*8 ddrdxic,ddrdyic,ddrdzic
      real*8 dtxiaxia,dtxiayia,dtxiazia
      real*8 dtxibxib,dtxibyib,dtxibzib
      real*8 dtxicxic,dtxicyic,dtxiczic
      real*8 dtyiayia,dtyiazia,dtziazia
      real*8 dtyibyib,dtyibzib,dtzibzib
      real*8 dtyicyic,dtyiczic,dtziczic
      real*8 dtxibxia,dtxibyia,dtxibzia
      real*8 dtyibxia,dtyibyia,dtyibzia
      real*8 dtzibxia,dtzibyia,dtzibzia
      real*8 dtxibxic,dtxibyic,dtxibzic
      real*8 dtyibxic,dtyibyic,dtyibzic
      real*8 dtzibxic,dtzibyic,dtzibzic
      real*8 dtxiaxic,dtxiayic,dtxiazic
      real*8 dtyiaxic,dtyiayic,dtyiazic
      real*8 dtziaxic,dtziayic,dtziazic
      real*8 drxiaxia,drxiayia,drxiazia
      real*8 drxibxib,drxibyib,drxibzib
      real*8 drxicxic,drxicyic,drxiczic
      real*8 dryiayia,dryiazia,drziazia
      real*8 dryibyib,dryibzib,drzibzib
      real*8 dryicyic,dryiczic,drziczic
      logical proceed
c
c
c     calculate the stretch-bend Hessian elements
c
      do istrbnd = 1, nstrbnd
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         force = ksb(istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,3,ia,ib,ic,0,0)
         if (proceed)  proceed = (iatom.eq.ia .or. iatom.eq.ib
     &                                 .or. iatom.eq.ic)
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
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            if (rp .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     first derivatives of angle with respect to coordinates
c
               dt = angle - anat(i)
               terma = -radian / (rab*rab*rp)
               termc = radian / (rcb*rcb*rp)
               ddtdxia = terma * (yab*zp-zab*yp)
               ddtdyia = terma * (zab*xp-xab*zp)
               ddtdzia = terma * (xab*yp-yab*xp)
               ddtdxic = termc * (ycb*zp-zcb*yp)
               ddtdyic = termc * (zcb*xp-xcb*zp)
               ddtdzic = termc * (xcb*yp-ycb*xp)
               ddtdxib = -ddtdxia - ddtdxic
               ddtdyib = -ddtdyia - ddtdyic
               ddtdzib = -ddtdzia - ddtdzic
c
c     abbreviations used in defining chain rule terms
c
               rab2 = 2.0d0 / (rab*rab)
               xrab = xab * rab2
               yrab = yab * rab2
               zrab = zab * rab2
               rcb2 = 2.0d0 / (rcb*rcb)
               xrcb = xcb * rcb2
               yrcb = ycb * rcb2
               zrcb = zcb * rcb2
               rp2 = 1.0d0 / (rp*rp)
               xabp = (yab*zp-zab*yp) * rp2
               yabp = (zab*xp-xab*zp) * rp2
               zabp = (xab*yp-yab*xp) * rp2
               xcbp = (ycb*zp-zcb*yp) * rp2
               ycbp = (zcb*xp-xcb*zp) * rp2
               zcbp = (xcb*yp-ycb*xp) * rp2
c
c     second derivatives of angle with respect to coordinates
c
               dtxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab)
               dtxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab)
               dtxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab)
               dtyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab)
               dtyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab)
               dtziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab)
               dtxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb)
               dtxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb)
               dtxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb)
               dtyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb)
               dtyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb)
               dtziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb)
               dtxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp
               dtxiayic = -terma*xab*yab - ddtdxia*yabp
               dtxiazic = -terma*xab*zab - ddtdxia*zabp
               dtyiaxic = -terma*xab*yab - ddtdyia*xabp
               dtyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp
               dtyiazic = -terma*yab*zab - ddtdyia*zabp
               dtziaxic = -terma*xab*zab - ddtdzia*xabp
               dtziayic = -terma*yab*zab - ddtdzia*yabp
               dtziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp
c
c     more angle deviation derivatives resulting from symmetry
c
               dtxibxia = -dtxiaxia - dtxiaxic
               dtxibyia = -dtxiayia - dtyiaxic
               dtxibzia = -dtxiazia - dtziaxic
               dtyibxia = -dtxiayia - dtxiayic
               dtyibyia = -dtyiayia - dtyiayic
               dtyibzia = -dtyiazia - dtziayic
               dtzibxia = -dtxiazia - dtxiazic
               dtzibyia = -dtyiazia - dtyiazic
               dtzibzia = -dtziazia - dtziazic
               dtxibxic = -dtxicxic - dtxiaxic
               dtxibyic = -dtxicyic - dtxiayic
               dtxibzic = -dtxiczic - dtxiazic
               dtyibxic = -dtxicyic - dtyiaxic
               dtyibyic = -dtyicyic - dtyiayic
               dtyibzic = -dtyiczic - dtyiazic
               dtzibxic = -dtxiczic - dtziaxic
               dtzibyic = -dtyiczic - dtziayic
               dtzibzic = -dtziczic - dtziazic
               dtxibxib = -dtxibxia - dtxibxic
               dtxibyib = -dtxibyia - dtxibyic
               dtxibzib = -dtxibzia - dtxibzic
               dtyibyib = -dtyibyia - dtyibyic
               dtyibzib = -dtyibzia - dtyibzic
               dtzibzib = -dtzibzia - dtzibzic
c
c     compute the values of the bond length deviations
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               term = stbnunit * force
               dr = term * (rab - bl(j) + rcb - bl(k))
               terma = term / rab
               termc = term / rcb
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dr = dr * fgrp
                  terma = terma * fgrp
                  termc = termc * fgrp
               end if
c
c     first derivatives of bond length with respect to coordinates
c
               ddrdxia = terma * xab
               ddrdyia = terma * yab
               ddrdzia = terma * zab
               ddrdxic = termc * xcb
               ddrdyic = termc * ycb
               ddrdzic = termc * zcb
               ddrdxib = -ddrdxia - ddrdxic
               ddrdyib = -ddrdyia - ddrdyic
               ddrdzib = -ddrdzia - ddrdzic
c
c     abbreviations used in defining chain rule terms
c
               xab = xab / rab
               yab = yab / rab
               zab = zab / rab
               xcb = xcb / rcb
               ycb = ycb / rcb
               zcb = zcb / rcb
c
c     second derivatives of bond length with respect to coordinates
c
               drxiaxia = terma * (1.0d0-xab**2)
               drxiayia = -terma * xab*yab
               drxiazia = -terma * xab*zab
               dryiayia = terma * (1.0d0-yab**2)
               dryiazia = -terma * yab*zab
               drziazia = terma * (1.0d0-zab**2)
               drxicxic = termc * (1.0d0-xcb**2)
               drxicyic = -termc * xcb*ycb
               drxiczic = -termc * xcb*zcb
               dryicyic = termc * (1.0d0-ycb**2)
               dryiczic = -termc * ycb*zcb
               drziczic = termc * (1.0d0-zcb**2)
               drxibxib = drxiaxia + drxicxic
               drxibyib = drxiayia + drxicyic
               drxibzib = drxiazia + drxiczic
               dryibyib = dryiayia + dryicyic
               dryibzib = dryiazia + dryiczic
               drzibzib = drziazia + drziczic
c
c     increment diagonal and non-diagonal Hessian elements
c
               if (ia .eq. iatom) then
                  hessx(1,ia) = hessx(1,ia) + dt*drxiaxia + dr*dtxiaxia
     &                             + 2.0d0*ddtdxia*ddrdxia
                  hessx(2,ia) = hessx(2,ia) + dt*drxiayia + dr*dtxiayia
     &                             + ddtdxia*ddrdyia + ddtdyia*ddrdxia
                  hessx(3,ia) = hessx(3,ia) + dt*drxiazia + dr*dtxiazia
     &                             + ddtdxia*ddrdzia + ddtdzia*ddrdxia
                  hessy(1,ia) = hessy(1,ia) + dt*drxiayia + dr*dtxiayia
     &                             + ddtdyia*ddrdxia + ddtdxia*ddrdyia
                  hessy(2,ia) = hessy(2,ia) + dt*dryiayia + dr*dtyiayia
     &                             + 2.0d0*ddtdyia*ddrdyia
                  hessy(3,ia) = hessy(3,ia) + dt*dryiazia + dr*dtyiazia
     &                             + ddtdyia*ddrdzia + ddtdzia*ddrdyia
                  hessz(1,ia) = hessz(1,ia) + dt*drxiazia + dr*dtxiazia
     &                             + ddtdzia*ddrdxia + ddtdxia*ddrdzia
                  hessz(2,ia) = hessz(2,ia) + dt*dryiazia + dr*dtyiazia
     &                             + ddtdzia*ddrdyia + ddtdyia*ddrdzia
                  hessz(3,ia) = hessz(3,ia) + dt*drziazia + dr*dtziazia
     &                             + 2.0d0*ddtdzia*ddrdzia
                  hessx(1,ib) = hessx(1,ib) - dt*drxiaxia + dr*dtxibxia
     &                             + ddtdxia*ddrdxib + ddtdxib*ddrdxia
                  hessx(2,ib) = hessx(2,ib) - dt*drxiayia + dr*dtxibyia
     &                             + ddtdxia*ddrdyib + ddtdyib*ddrdxia
                  hessx(3,ib) = hessx(3,ib) - dt*drxiazia + dr*dtxibzia
     &                             + ddtdxia*ddrdzib + ddtdzib*ddrdxia
                  hessy(1,ib) = hessy(1,ib) - dt*drxiayia + dr*dtyibxia
     &                             + ddtdyia*ddrdxib + ddtdxib*ddrdyia
                  hessy(2,ib) = hessy(2,ib) - dt*dryiayia + dr*dtyibyia
     &                             + ddtdyia*ddrdyib + ddtdyib*ddrdyia
                  hessy(3,ib) = hessy(3,ib) - dt*dryiazia + dr*dtyibzia
     &                             + ddtdyia*ddrdzib + ddtdzib*ddrdyia
                  hessz(1,ib) = hessz(1,ib) - dt*drxiazia + dr*dtzibxia
     &                             + ddtdzia*ddrdxib + ddtdxib*ddrdzia
                  hessz(2,ib) = hessz(2,ib) - dt*dryiazia + dr*dtzibyia
     &                             + ddtdzia*ddrdyib + ddtdyib*ddrdzia
                  hessz(3,ib) = hessz(3,ib) - dt*drziazia + dr*dtzibzia
     &                             + ddtdzia*ddrdzib + ddtdzib*ddrdzia
                  hessx(1,ic) = hessx(1,ic) + dr*dtxiaxic
     &                             + ddtdxia*ddrdxic + ddtdxic*ddrdxia
                  hessx(2,ic) = hessx(2,ic) + dr*dtxiayic
     &                             + ddtdxia*ddrdyic + ddtdyic*ddrdxia
                  hessx(3,ic) = hessx(3,ic) + dr*dtxiazic
     &                             + ddtdxia*ddrdzic + ddtdzic*ddrdxia
                  hessy(1,ic) = hessy(1,ic) + dr*dtyiaxic
     &                             + ddtdyia*ddrdxic + ddtdxic*ddrdyia
                  hessy(2,ic) = hessy(2,ic) + dr*dtyiayic
     &                             + ddtdyia*ddrdyic + ddtdyic*ddrdyia
                  hessy(3,ic) = hessy(3,ic) + dr*dtyiazic
     &                             + ddtdyia*ddrdzic + ddtdzic*ddrdyia
                  hessz(1,ic) = hessz(1,ic) + dr*dtziaxic
     &                             + ddtdzia*ddrdxic + ddtdxic*ddrdzia
                  hessz(2,ic) = hessz(2,ic) + dr*dtziayic
     &                             + ddtdzia*ddrdyic + ddtdyic*ddrdzia
                  hessz(3,ic) = hessz(3,ic) + dr*dtziazic
     &                             + ddtdzia*ddrdzic + ddtdzic*ddrdzia
               else if (ib .eq. iatom) then
                  hessx(1,ib) = hessx(1,ib) + dt*drxibxib + dr*dtxibxib
     &                             + 2.0d0*ddtdxib*ddrdxib
                  hessx(2,ib) = hessx(2,ib) + dt*drxibyib + dr*dtxibyib
     &                             + ddtdxib*ddrdyib + ddtdyib*ddrdxib
                  hessx(3,ib) = hessx(3,ib) + dt*drxibzib + dr*dtxibzib
     &                             + ddtdxib*ddrdzib + ddtdzib*ddrdxib
                  hessy(1,ib) = hessy(1,ib) + dt*drxibyib + dr*dtxibyib
     &                             + ddtdyib*ddrdxib + ddtdxib*ddrdyib
                  hessy(2,ib) = hessy(2,ib) + dt*dryibyib + dr*dtyibyib
     &                             + 2.0d0*ddtdyib*ddrdyib
                  hessy(3,ib) = hessy(3,ib) + dt*dryibzib + dr*dtyibzib
     &                             + ddtdyib*ddrdzib + ddtdzib*ddrdyib
                  hessz(1,ib) = hessz(1,ib) + dt*drxibzib + dr*dtxibzib
     &                             + ddtdzib*ddrdxib + ddtdxib*ddrdzib
                  hessz(2,ib) = hessz(2,ib) + dt*dryibzib + dr*dtyibzib
     &                             + ddtdzib*ddrdyib + ddtdyib*ddrdzib
                  hessz(3,ib) = hessz(3,ib) + dt*drzibzib + dr*dtzibzib
     &                             + 2.0d0*ddtdzib*ddrdzib
                  hessx(1,ia) = hessx(1,ia) - dt*drxiaxia + dr*dtxibxia
     &                             + ddtdxib*ddrdxia + ddtdxia*ddrdxib
                  hessx(2,ia) = hessx(2,ia) - dt*drxiayia + dr*dtxibyia
     &                             + ddtdxib*ddrdyia + ddtdyia*ddrdxib
                  hessx(3,ia) = hessx(3,ia) - dt*drxiazia + dr*dtxibzia
     &                             + ddtdxib*ddrdzia + ddtdzia*ddrdxib
                  hessy(1,ia) = hessy(1,ia) - dt*drxiayia + dr*dtyibxia
     &                             + ddtdyib*ddrdxia + ddtdxia*ddrdyib
                  hessy(2,ia) = hessy(2,ia) - dt*dryiayia + dr*dtyibyia
     &                             + ddtdyib*ddrdyia + ddtdyia*ddrdyib
                  hessy(3,ia) = hessy(3,ia) - dt*dryiazia + dr*dtyibzia
     &                             + ddtdyib*ddrdzia + ddtdzia*ddrdyib
                  hessz(1,ia) = hessz(1,ia) - dt*drxiazia + dr*dtzibxia
     &                             + ddtdzib*ddrdxia + ddtdxia*ddrdzib
                  hessz(2,ia) = hessz(2,ia) - dt*dryiazia + dr*dtzibyia
     &                             + ddtdzib*ddrdyia + ddtdyia*ddrdzib
                  hessz(3,ia) = hessz(3,ia) - dt*drziazia + dr*dtzibzia
     &                             + ddtdzib*ddrdzia + ddtdzia*ddrdzib
                  hessx(1,ic) = hessx(1,ic) - dt*drxicxic + dr*dtxibxic
     &                             + ddtdxib*ddrdxic + ddtdxic*ddrdxib
                  hessx(2,ic) = hessx(2,ic) - dt*drxicyic + dr*dtxibyic
     &                             + ddtdxib*ddrdyic + ddtdyic*ddrdxib
                  hessx(3,ic) = hessx(3,ic) - dt*drxiczic + dr*dtxibzic
     &                             + ddtdxib*ddrdzic + ddtdzic*ddrdxib
                  hessy(1,ic) = hessy(1,ic) - dt*drxicyic + dr*dtyibxic
     &                             + ddtdyib*ddrdxic + ddtdxic*ddrdyib
                  hessy(2,ic) = hessy(2,ic) - dt*dryicyic + dr*dtyibyic
     &                             + ddtdyib*ddrdyic + ddtdyic*ddrdyib
                  hessy(3,ic) = hessy(3,ic) - dt*dryiczic + dr*dtyibzic
     &                             + ddtdyib*ddrdzic + ddtdzic*ddrdyib
                  hessz(1,ic) = hessz(1,ic) - dt*drxiczic + dr*dtzibxic
     &                             + ddtdzib*ddrdxic + ddtdxic*ddrdzib
                  hessz(2,ic) = hessz(2,ic) - dt*dryiczic + dr*dtzibyic
     &                             + ddtdzib*ddrdyic + ddtdyic*ddrdzib
                  hessz(3,ic) = hessz(3,ic) - dt*drziczic + dr*dtzibzic
     &                             + ddtdzib*ddrdzic + ddtdzic*ddrdzib
               else if (ic .eq. iatom) then
                  hessx(1,ic) = hessx(1,ic) + dt*drxicxic + dr*dtxicxic
     &                            + 2.0d0*ddtdxic*ddrdxic
                  hessx(2,ic) = hessx(2,ic) + dt*drxicyic + dr*dtxicyic
     &                            + ddtdxic*ddrdyic + ddtdyic*ddrdxic
                  hessx(3,ic) = hessx(3,ic) + dt*drxiczic + dr*dtxiczic
     &                            + ddtdxic*ddrdzic + ddtdzic*ddrdxic
                  hessy(1,ic) = hessy(1,ic) + dt*drxicyic + dr*dtxicyic
     &                            + ddtdyic*ddrdxic + ddtdxic*ddrdyic
                  hessy(2,ic) = hessy(2,ic) + dt*dryicyic + dr*dtyicyic
     &                            + 2.0d0*ddtdyic*ddrdyic
                  hessy(3,ic) = hessy(3,ic) + dt*dryiczic + dr*dtyiczic
     &                            + ddtdyic*ddrdzic + ddtdzic*ddrdyic
                  hessz(1,ic) = hessz(1,ic) + dt*drxiczic + dr*dtxiczic
     &                            + ddtdzic*ddrdxic + ddtdxic*ddrdzic
                  hessz(2,ic) = hessz(2,ic) + dt*dryiczic + dr*dtyiczic
     &                            + ddtdzic*ddrdyic + ddtdyic*ddrdzic
                  hessz(3,ic) = hessz(3,ic) + dt*drziczic + dr*dtziczic
     &                            + 2.0d0*ddtdzic*ddrdzic
                  hessx(1,ib) = hessx(1,ib) - dt*drxicxic + dr*dtxibxic
     &                            + ddtdxic*ddrdxib + ddtdxib*ddrdxic
                  hessx(2,ib) = hessx(2,ib) - dt*drxicyic + dr*dtxibyic
     &                            + ddtdxic*ddrdyib + ddtdyib*ddrdxic
                  hessx(3,ib) = hessx(3,ib) - dt*drxiczic + dr*dtxibzic
     &                            + ddtdxic*ddrdzib + ddtdzib*ddrdxic
                  hessy(1,ib) = hessy(1,ib) - dt*drxicyic + dr*dtyibxic
     &                            + ddtdyic*ddrdxib + ddtdxib*ddrdyic
                  hessy(2,ib) = hessy(2,ib) - dt*dryicyic + dr*dtyibyic
     &                            + ddtdyic*ddrdyib + ddtdyib*ddrdyic
                  hessy(3,ib) = hessy(3,ib) - dt*dryiczic + dr*dtyibzic
     &                            + ddtdyic*ddrdzib + ddtdzib*ddrdyic
                  hessz(1,ib) = hessz(1,ib) - dt*drxiczic + dr*dtzibxic
     &                            + ddtdzic*ddrdxib + ddtdxib*ddrdzic
                  hessz(2,ib) = hessz(2,ib) - dt*dryiczic + dr*dtzibyic
     &                            + ddtdzic*ddrdyib + ddtdyib*ddrdzic
                  hessz(3,ib) = hessz(3,ib) - dt*drziczic + dr*dtzibzic
     &                            + ddtdzic*ddrdzib + ddtdzib*ddrdzic
                  hessx(1,ia) = hessx(1,ia) + dr*dtxiaxic
     &                            + ddtdxic*ddrdxia + ddtdxia*ddrdxic
                  hessx(2,ia) = hessx(2,ia) + dr*dtyiaxic
     &                            + ddtdxic*ddrdyia + ddtdyia*ddrdxic
                  hessx(3,ia) = hessx(3,ia) + dr*dtziaxic
     &                            + ddtdxic*ddrdzia + ddtdzia*ddrdxic
                  hessy(1,ia) = hessy(1,ia) + dr*dtxiayic
     &                            + ddtdyic*ddrdxia + ddtdxia*ddrdyic
                  hessy(2,ia) = hessy(2,ia) + dr*dtyiayic
     &                            + ddtdyic*ddrdyia + ddtdyia*ddrdyic
                  hessy(3,ia) = hessy(3,ia) + dr*dtziayic
     &                            + ddtdyic*ddrdzia + ddtdzia*ddrdyic
                  hessz(1,ia) = hessz(1,ia) + dr*dtxiazic
     &                            + ddtdzic*ddrdxia + ddtdxia*ddrdzic
                  hessz(2,ia) = hessz(2,ia) + dr*dtyiazic
     &                            + ddtdzic*ddrdyia + ddtdyia*ddrdzic
                  hessz(3,ia) = hessz(3,ia) + dr*dtziazic
     &                            + ddtdzic*ddrdzia + ddtdzia*ddrdzic
               end if
            end if
         end if
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
c     ##  subroutine estrbnd3  --  stretch-bend energy & analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "estrbnd3" calculates the stretch-bend potential energy;
c     also partitions the energy among the atoms
c
c
      subroutine estrbnd3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'strbnd.i'
      include 'usage.i'
      integer i,j,k,istrbnd
      integer ia,ib,ic
      real*8 e,dr,dt,fgrp
      real*8 angle,force
      real*8 dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 rab,xab,yab,zab
      real*8 rcb,xcb,ycb,zcb
      logical header,huge,proceed
c
c
c     zero out the energy component and partitioning terms
c
      neba = 0
      eba = 0.0d0
      do i = 1, n
         aeba(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the stretch-bend energy term
c
      do istrbnd = 1, nstrbnd
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         force = ksb(istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,3,ia,ib,ic,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
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
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            if (rab*rcb .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt = angle - anat(i)
c
c     get the stretch-bend interaction energy
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               dr = rab - bl(j) + rcb - bl(k)
               e = stbnunit * force * dt * dr
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total stretch-bend energy
c
               neba = neba + 1
               eba = eba + e
               aeba(ib) = aeba(ib) + e
c
c     print a warning if the energy of this interaction is large
c
               huge = (abs(e) .gt. 2.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Stretch-Bend Cross Term',
     &                          ' Interactions :',
     &                       //,' Type',15x,'Atom Names',16x,'dBond',
     &                          4x,'dAngle',6x,'Energy',/)
                  end if
                  write (iout,20)  ia,name(ia),ib,name(ib),
     &                             ic,name(ic),dr,dt,e
   20             format (' StrBend  ',i5,'-',a3,1x,i5,'-',a3,
     &                       1x,i5,'-',a3,2x,2f10.4,f12.4)
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
c     ##  subroutine estrtor  --  stretch-torsion cross term energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "estrtor" calculates the stretch-torsion potential energy
c
c
      subroutine estrtor
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bond.i'
      include 'energi.i'
      include 'group.i'
      include 'strtor.i'
      include 'torpot.i'
      include 'tors.i'
      include 'usage.i'
      integer i,k,ia,ib,ic,id,istrtor
      real*8 e,rcb,dr,fgrp
      real*8 rt2,ru2,rtru
      real*8 xt,yt,zt,xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3,s1,s2,s3
      real*8 cosine,cosine2,cosine3
      real*8 sine,sine2,sine3
      real*8 phi1,phi2,phi3
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      logical proceed
c
c
c     zero out the stretch-torsion energy
c
      ebt = 0.0d0
c
c     calculate the stretch-torsion interaction energy term
c
      do istrtor = 1, nstrtor
         i = ist(1,istrtor)
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the stretch-torsional parameters for this angle
c
               v1 = kst(1,istrtor)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = kst(2,istrtor)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = kst(3,istrtor)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
c
c     calculate the bond-stretch for the central bond
c
               k = ist(2,istrtor)
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               dr = rcb - bl(k)
c
c     compute the stretch-torsion energy for this angle
c
               e = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total stretch-torsion energy
c
               ebt = ebt + e
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrtor1  --  stretch-torsion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrtor1" calculates the stretch-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
      subroutine estrtor1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bond.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'strtor.i'
      include 'torpot.i'
      include 'tors.i'
      include 'usage.i'
      include 'virial.i'
      integer i,k,ia,ib,ic,id,istrtor
      real*8 e,dedphi,dr,ddr,fgrp
      real*8 rt2,ru2,rtru,rcb
      real*8 ddrdx,ddrdy,ddrdz
      real*8 v1,v2,v3
      real*8 c1,c2,c3,s1,s2,s3
      real*8 cosine,cosine2,cosine3
      real*8 sine,sine2,sine3
      real*8 phi1,phi2,phi3
      real*8 dphi1,dphi2,dphi3
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 xt,yt,zt,xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      logical proceed
c
c
c     zero out the stretch-torsion energy and first derivatives
c
      ebt = 0.0d0
      do i = 1, n
         debt(1,i) = 0.0d0
         debt(2,i) = 0.0d0
         debt(3,i) = 0.0d0
      end do
c
c     calculate the stretch-torsion interaction energy term
c
      do istrtor = 1, nstrtor
         i = ist(1,istrtor)
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the stretch-torsional parameters for this angle
c
               v1 = kst(1,istrtor)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = kst(2,istrtor)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = kst(3,istrtor)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
c
c     calculate the bond-stretch for the central bond
c
               k = ist(2,istrtor)
               dr = rcb - bl(k)
c
c     calculate stretch-torsion energy and chain rule terms
c
               e = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = storunit * dr * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddr = e / (dr * rcb)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedphi = dedphi * fgrp
                  ddr = ddr * fgrp
               end if
c
c     first direivative components for the bond stretch
c
               ddrdx = xcb * ddr
               ddrdy = ycb * ddr
               ddrdz = zcb * ddr
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = yca*dphidzt - zca*dphidyt
     &                       + zdc*dphidyu - ydc*dphidzu
               dphidyib = zca*dphidxt - xca*dphidzt
     &                       + xdc*dphidzu - zdc*dphidxu
               dphidzib = xca*dphidyt - yca*dphidxt
     &                       + ydc*dphidxu - xdc*dphidyu
               dphidxic = zba*dphidyt - yba*dphidzt
     &                       + ydb*dphidzu - zdb*dphidyu
               dphidyic = xba*dphidzt - zba*dphidxt
     &                       + zdb*dphidxu - xdb*dphidzu
               dphidzic = yba*dphidxt - xba*dphidyt
     &                       + xdb*dphidyu - ydb*dphidxu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     compute derivative components for this interaction
c
               dedxia = dedphi*dphidxia
               dedyia = dedphi*dphidyia
               dedzia = dedphi*dphidzia
               dedxib = dedphi*dphidxib - ddrdx
               dedyib = dedphi*dphidyib - ddrdy
               dedzib = dedphi*dphidzib - ddrdz
               dedxic = dedphi*dphidxic + ddrdx
               dedyic = dedphi*dphidyic + ddrdy
               dedzic = dedphi*dphidzic + ddrdz
               dedxid = dedphi*dphidxid
               dedyid = dedphi*dphidyid
               dedzid = dedphi*dphidzid
c
c     increment the stretch-torsion energy and gradient
c
               ebt = ebt + e
               debt(1,ia) = debt(1,ia) + dedxia
               debt(2,ia) = debt(2,ia) + dedyia
               debt(3,ia) = debt(3,ia) + dedzia
               debt(1,ib) = debt(1,ib) + dedxib
               debt(2,ib) = debt(2,ib) + dedyib
               debt(3,ib) = debt(3,ib) + dedzib
               debt(1,ic) = debt(1,ic) + dedxic
               debt(2,ic) = debt(2,ic) + dedyic
               debt(3,ic) = debt(3,ic) + dedzic
               debt(1,id) = debt(1,id) + dedxid
               debt(2,id) = debt(2,id) + dedyid
               debt(3,id) = debt(3,id) + dedzid
c
c     increment the components of the virial
c
               if (isobaric) then
                  virx = virx - xba*dedxia + xdc*dedxid
     &                      + xcb*(dedxic+dedxid)
                  viry = viry - yba*dedyia + ydc*dedyid
     &                      + ycb*(dedyic+dedyid)
                  virz = virz - zba*dedzia + zdc*dedzid
     &                      + zcb*(dedzic+dedzid)
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrtor2  --  stretch-torsion Hessian; analyt  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrtor2" calculates the stretch-torsion potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine estrtor2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bond.i'
      include 'group.i'
      include 'hessn.i'
      include 'strtor.i'
      include 'torpot.i'
      include 'tors.i'
      integer i,j,k,ia,ib,ic,id,istrtor
      real*8 dedphi,d2edphi2,fgrp
      real*8 rt2,ru2,rtru,rcb
      real*8 v1,v2,v3,c1,c2,c3,s1,s2,s3
      real*8 cosine,cosine2,cosine3
      real*8 sine,sine2,sine3
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 phi1,phi2,phi3
      real*8 dphi1,dphi2,dphi3
      real*8 d2phi1,d2phi2,d2phi3
      real*8 dr,ddr,d2dr
      real*8 ddrdx,ddrdy,ddrdz
      real*8 d2drdxx,d2drdyy,d2drdzz
      real*8 d2drdxy,d2drdxz,d2drdyz
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 dphidxibt,dphidyibt,dphidzibt
      real*8 dphidxibu,dphidyibu,dphidzibu
      real*8 dphidxict,dphidyict,dphidzict
      real*8 dphidxicu,dphidyicu,dphidzicu
      real*8 dxia,dyia,dzia,dxib,dyib,dzib
      real*8 dxic,dyic,dzic,dxid,dyid,dzid
      real*8 dxiaxia,dyiayia,dziazia,dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic,dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia,dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic,dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib,dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib,dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic,dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid,dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid,dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic,dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid,dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid,dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid,dzicxid,dzicyid,dziczid
      logical proceed
c
c
c     for each stretch-torsion interaction we first calculate
c     the cosine of the dihedral between ia-ib-ic-id
c
      do istrtor = 1, nstrtor
         j = ist(1,istrtor)
         ia = itors(1,j)
         ib = itors(2,j)
         ic = itors(3,j)
         id = itors(4,j)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (i.eq.ia .or. i.eq.ib .or.
     &                              i.eq.ic .or. i.eq.id)
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the stretch-torsional parameters for this angle
c
               v1 = kst(1,istrtor)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = kst(2,istrtor)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = kst(3,istrtor)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
               d2phi1 = -(cosine*c1 + sine*s1)
               d2phi2 = -4.0d0 * (cosine2*c2 + sine2*s2)
               d2phi3 = -9.0d0 * (cosine3*c3 + sine3*s3)
c
c     calculate the bond-stretch for the central bond
c
               k = ist(2,istrtor)
               dr = rcb - bl(k)
c
c     calculate the stretch-torsion master chain rule terms
c
               dedphi = storunit * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = storunit * dr
     &                       * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3)
               ddr = 1.0d0 / rcb
               d2dr = -storunit * (v1*phi1 + v2*phi2 + v3*phi3) / rcb**3
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dr = d2dr * fgrp
               end if
c
c     first and second derivative components for the bond stretch
c
               ddrdx = xcb * ddr
               ddrdy = ycb * ddr
               ddrdz = zcb * ddr
               d2drdxx = (xcb**2-rcb**2) * d2dr
               d2drdyy = (ycb**2-rcb**2) * d2dr
               d2drdzz = (zcb**2-rcb**2) * d2dr
               d2drdxy = xcb * ycb * d2dr
               d2drdxz = xcb * zcb * d2dr
               d2drdyz = ycb * zcb * d2dr
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
               xycb2 = xcb*xcb + ycb*ycb
               xzcb2 = xcb*xcb + zcb*zcb
               yzcb2 = ycb*ycb + zcb*zcb
               rcbxt = -2.0d0 * rcb * dphidxt
               rcbyt = -2.0d0 * rcb * dphidyt
               rcbzt = -2.0d0 * rcb * dphidzt
               rcbt2 = rcb * rt2
               rcbxu = 2.0d0 * rcb * dphidxu
               rcbyu = 2.0d0 * rcb * dphidyu
               rcbzu = 2.0d0 * rcb * dphidzu
               rcbu2 = rcb * ru2
               dphidxibt = yca*dphidzt - zca*dphidyt
               dphidxibu = zdc*dphidyu - ydc*dphidzu
               dphidyibt = zca*dphidxt - xca*dphidzt
               dphidyibu = xdc*dphidzu - zdc*dphidxu
               dphidzibt = xca*dphidyt - yca*dphidxt
               dphidzibu = ydc*dphidxu - xdc*dphidyu
               dphidxict = zba*dphidyt - yba*dphidzt
               dphidxicu = ydb*dphidzu - zdb*dphidyu
               dphidyict = xba*dphidzt - zba*dphidxt
               dphidyicu = zdb*dphidxu - xdb*dphidzu
               dphidzict = yba*dphidxt - xba*dphidyt
               dphidzicu = xdb*dphidyu - ydb*dphidxu
c
c     intermediate terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = dphidxibt + dphidxibu
               dphidyib = dphidyibt + dphidyibu
               dphidzib = dphidzibt + dphidzibu
               dphidxic = dphidxict + dphidxicu
               dphidyic = dphidyict + dphidyicu
               dphidzic = dphidzict + dphidzicu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     chain rule terms for first derivative components
c
               dxia = dedphi * dphidxia
               dyia = dedphi * dphidyia
               dzia = dedphi * dphidzia
               dxib = dedphi * dphidxib
               dyib = dedphi * dphidyib
               dzib = dedphi * dphidzib
               dxic = dedphi * dphidxic
               dyic = dedphi * dphidyic
               dzic = dedphi * dphidzic
               dxid = dedphi * dphidxid
               dyid = dedphi * dphidyid
               dzid = dedphi * dphidzid
               dedphi = dedphi * dr
c
c     chain rule terms for second derivative components
c
               dxiaxia = rcbxt*dphidxia
               dxiayia = rcbxt*dphidyia - zcb*rcb/rt2
               dxiazia = rcbxt*dphidzia + ycb*rcb/rt2
               dxiaxib = rcbxt*dphidxibt + xcb*(zca*ycb-yca*zcb)/rcbt2
               dxiayib = rcbxt*dphidyibt + dphidzt
     &                      + (xca*zcb*xcb+zca*yzcb2)/rcbt2
               dxiazib = rcbxt*dphidzibt - dphidyt
     &                      - (xca*ycb*xcb+yca*yzcb2)/rcbt2
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2
               dxiayic = rcbxt*dphidyict - dphidzt
     &                      - (xba*zcb*xcb+zba*yzcb2)/rcbt2
               dxiazic = rcbxt*dphidzict + dphidyt
     &                      + (xba*ycb*xcb+yba*yzcb2)/rcbt2
               dxiaxid = 0.0d0
               dxiayid = 0.0d0
               dxiazid = 0.0d0
               dyiayia = rcbyt*dphidyia
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2
               dyiaxib = rcbyt*dphidxibt - dphidzt
     &                      - (yca*zcb*ycb+zca*xzcb2)/rcbt2
               dyiayib = rcbyt*dphidyibt + ycb*(xca*zcb-zca*xcb)/rcbt2
               dyiazib = rcbyt*dphidzibt + dphidxt
     &                      + (yca*xcb*ycb+xca*xzcb2)/rcbt2
               dyiaxic = rcbyt*dphidxict + dphidzt
     &                      + (yba*zcb*ycb+zba*xzcb2)/rcbt2
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2
               dyiazic = rcbyt*dphidzict - dphidxt
     &                      - (yba*xcb*ycb+xba*xzcb2)/rcbt2
               dyiaxid = 0.0d0
               dyiayid = 0.0d0
               dyiazid = 0.0d0
               dziazia = rcbzt*dphidzia
               dziaxib = rcbzt*dphidxibt + dphidyt
     &                      + (zca*ycb*zcb+yca*xycb2)/rcbt2
               dziayib = rcbzt*dphidyibt - dphidxt
     &                      - (zca*xcb*zcb+xca*xycb2)/rcbt2
               dziazib = rcbzt*dphidzibt + zcb*(yca*xcb-xca*ycb)/rcbt2
               dziaxic = rcbzt*dphidxict - dphidyt
     &                      - (zba*ycb*zcb+yba*xycb2)/rcbt2
               dziayic = rcbzt*dphidyict + dphidxt
     &                      + (zba*xcb*zcb+xba*xycb2)/rcbt2
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2
               dziaxid = 0.0d0
               dziayid = 0.0d0
               dziazid = 0.0d0
               dxibxic = -xcb*dphidxib/(rcb*rcb)
     &             - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidxibt/rt2
     &             - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidxibu/ru2
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
     &             - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidxibt/rt2
     &             + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidxibu/ru2
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2
               dxibyid = rcbyu*dphidxibu - dphidzu
     &                      - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
               dxibzid = rcbzu*dphidxibu + dphidyu
     &                      + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
               dyibzib = ycb*dphidzib/(rcb*rcb)
     &             - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &             - 2.0d0*(xt*zca-xca*zt)*dphidzibt/rt2
     &             + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &             + 2.0d0*(xu*zdc-xdc*zu)*dphidzibu/ru2
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu
     &             + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidyibt/rt2
     &             - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidyibu/ru2
               dyibyic = -ycb*dphidyib/(rcb*rcb)
     &             - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidyibt/rt2
     &             - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidyibu/ru2
               dyibxid = rcbxu*dphidyibu + dphidzu
     &                      + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2
               dyibzid = rcbzu*dphidyibu - dphidxu
     &                      - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu
     &             - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidzibt/rt2
     &             + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidzibu/ru2
               dzibzic = -zcb*dphidzib/(rcb*rcb)
     &             - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &             - 2.0d0*(xt*yba-xba*yt)*dphidzibt/rt2
     &             - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &             + 2.0d0*(xu*ydb-xdb*yu)*dphidzibu/ru2
               dzibxid = rcbxu*dphidzibu - dphidyu
     &                      - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
               dzibyid = rcbyu*dphidzibu + dphidxu
     &                      + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
               dxicyid = rcbyu*dphidxicu + dphidzu
     &                      + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
               dxiczid = rcbzu*dphidxicu - dphidyu
     &                      - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
               dyicxid = rcbxu*dphidyicu - dphidzu
     &                      - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
               dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
               dyiczid = rcbzu*dphidyicu + dphidxu
     &                      + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
               dzicxid = rcbxu*dphidzicu + dphidyu
     &                      + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
               dzicyid = rcbyu*dphidzicu - dphidxu
     &                      - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
               dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
               dxidxid = rcbxu*dphidxid
               dxidyid = rcbxu*dphidyid + zcb*rcb/ru2
               dxidzid = rcbxu*dphidzid - ycb*rcb/ru2
               dyidyid = rcbyu*dphidyid
               dyidzid = rcbyu*dphidzid + xcb*rcb/ru2
               dzidzid = rcbzu*dphidzid
c
c     get some second derivative chain rule terms by difference
c
               dxibxib = -dxiaxib - dxibxic - dxibxid
               dxibyib = -dyiaxib - dxibyic - dxibyid
               dxibzib = -dxiazib - dzibxic - dzibxid
               dxibzic = -dziaxib - dxibzib - dxibzid
               dyibyib = -dyiayib - dyibyic - dyibyid
               dyibzic = -dziayib - dyibzib - dyibzid
               dzibzib = -dziazib - dzibzic - dzibzid
               dzibyic = -dyiazib - dyibzib - dzibyid
               dxicxic = -dxiaxic - dxibxic - dxicxid
               dxicyic = -dyiaxic - dyibxic - dxicyid
               dxiczic = -dziaxic - dzibxic - dxiczid
               dyicyic = -dyiayic - dyibyic - dyicyid
               dyiczic = -dziayic - dzibyic - dyiczid
               dziczic = -dziazic - dzibzic - dziczid
c
c     now, increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                              + d2edphi2*dphidxia*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                              + d2edphi2*dphidyia*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                              + d2edphi2*dphidzia*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxia*dphidxib
     &                              - dxia*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                              + d2edphi2*dphidyia*dphidxib
     &                              - dyia*ddrdx
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                              + d2edphi2*dphidzia*dphidxib
     &                              - dzia*ddrdx
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                              + d2edphi2*dphidxia*dphidyib
     &                              - dxia*ddrdy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                              + d2edphi2*dphidyia*dphidyib
     &                              - dyia*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                              + d2edphi2*dphidzia*dphidyib
     &                              - dzia*ddrdy
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                              + d2edphi2*dphidxia*dphidzib
     &                              - dxia*ddrdz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                              + d2edphi2*dphidyia*dphidzib
     &                              - dyia*ddrdz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                              + d2edphi2*dphidzia*dphidzib
     &                              - dzia*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxia*dphidxic
     &                              + dxia*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                              + d2edphi2*dphidyia*dphidxic
     &                              + dyia*ddrdx
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                              + d2edphi2*dphidzia*dphidxic
     &                              + dzia*ddrdx
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                              + d2edphi2*dphidxia*dphidyic
     &                              + dxia*ddrdy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                              + d2edphi2*dphidyia*dphidyic
     &                              + dyia*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                              + d2edphi2*dphidzia*dphidyic
     &                              + dzia*ddrdy
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                              + d2edphi2*dphidxia*dphidzic
     &                              + dxia*ddrdz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                              + d2edphi2*dphidyia*dphidzic
     &                              + dyia*ddrdz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                              + d2edphi2*dphidzia*dphidzic
     &                              + dzia*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxia*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                              + d2edphi2*dphidyia*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                              + d2edphi2*dphidzia*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                              + d2edphi2*dphidxia*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                              + d2edphi2*dphidyia*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                              + d2edphi2*dphidzia*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                              + d2edphi2*dphidxia*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                              + d2edphi2*dphidyia*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                              + d2edphi2*dphidzia*dphidzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                              + d2edphi2*dphidxib*dphidxib
     &                              - 2.0d0*dxib*ddrdx + d2drdxx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
     &                              - dyib*ddrdx - dxib*ddrdy + d2drdxy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
     &                              - dzib*ddrdx - dxib*ddrdz + d2drdxz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
     &                              - dxib*ddrdy - dyib*ddrdx + d2drdxy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                              + d2edphi2*dphidyib*dphidyib
     &                              - 2.0d0*dyib*ddrdy + d2drdyy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
     &                              - dzib*ddrdy - dyib*ddrdz + d2drdyz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
     &                              - dxib*ddrdz - dzib*ddrdx + d2drdxz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
     &                              - dyib*ddrdz - dzib*ddrdy + d2drdyz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                              + d2edphi2*dphidzib*dphidzib
     &                              - 2.0d0*dzib*ddrdz + d2drdzz
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxib*dphidxia
     &                              - dxia*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                              + d2edphi2*dphidyib*dphidxia
     &                              - dxia*ddrdy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                              + d2edphi2*dphidzib*dphidxia
     &                              - dxia*ddrdz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                              + d2edphi2*dphidxib*dphidyia
     &                              - dyia*ddrdx
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                              + d2edphi2*dphidyib*dphidyia
     &                              - dyia*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                              + d2edphi2*dphidzib*dphidyia
     &                              - dyia*ddrdz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                              + d2edphi2*dphidxib*dphidzia
     &                              - dzia*ddrdx
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                              + d2edphi2*dphidyib*dphidzia
     &                              - dzia*ddrdy
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                              + d2edphi2*dphidzib*dphidzia
     &                              - dzia*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                              + d2edphi2*dphidxib*dphidxic
     &                              + (dxib-dxic)*ddrdx - d2drdxx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                              + d2edphi2*dphidyib*dphidxic
     &                              + dyib*ddrdx - dxic*ddrdy - d2drdxy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                              + d2edphi2*dphidzib*dphidxic
     &                              + dzib*ddrdx - dxic*ddrdz - d2drdxz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                              + d2edphi2*dphidxib*dphidyic
     &                              + dxib*ddrdy - dyic*ddrdx - d2drdxy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                              + d2edphi2*dphidyib*dphidyic
     &                              + (dyib-dyic)*ddrdy - d2drdyy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                              + d2edphi2*dphidzib*dphidyic
     &                              + dzib*ddrdy - dyic*ddrdz - d2drdyz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                              + d2edphi2*dphidxib*dphidzic
     &                              + dxib*ddrdz - dzic*ddrdx - d2drdxz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                              + d2edphi2*dphidyib*dphidzic
     &                              + dyib*ddrdz - dzic*ddrdy - d2drdyz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                              + d2edphi2*dphidzib*dphidzic
     &                              + (dzib-dzic)*ddrdz - d2drdzz
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                              + d2edphi2*dphidxib*dphidxid
     &                              - dxid*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                              + d2edphi2*dphidyib*dphidxid
     &                              - dxid*ddrdy
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                              + d2edphi2*dphidzib*dphidxid
     &                              - dxid*ddrdz
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                              + d2edphi2*dphidxib*dphidyid
     &                              - dyid*ddrdx
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                              + d2edphi2*dphidyib*dphidyid
     &                              - dyid*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                              + d2edphi2*dphidzib*dphidyid
     &                              - dyid*ddrdz
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                              + d2edphi2*dphidxib*dphidzid
     &                              - dzid*ddrdx
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                              + d2edphi2*dphidyib*dphidzid
     &                              - dzid*ddrdy
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                              + d2edphi2*dphidzib*dphidzid
     &                              - dzid*ddrdz
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                              + d2edphi2*dphidxic*dphidxic
     &                              + 2.0d0*dxic*ddrdx + d2drdxx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
     &                              + dyic*ddrdx + dxic*ddrdy + d2drdxy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
     &                              + dzic*ddrdx + dxic*ddrdz + d2drdxz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
     &                              + dxic*ddrdy + dyic*ddrdx + d2drdxy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                              + d2edphi2*dphidyic*dphidyic
     &                              + 2.0d0*dyic*ddrdy + d2drdyy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
     &                              + dzic*ddrdy + dyic*ddrdz + d2drdyz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
     &                              + dxic*ddrdz + dzic*ddrdx + d2drdxz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
     &                              + dyic*ddrdz + dzic*ddrdy + d2drdyz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                              + d2edphi2*dphidzic*dphidzic
     &                              + 2.0d0*dzic*ddrdz + d2drdzz
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxic*dphidxia
     &                              + dxia*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                              + d2edphi2*dphidyic*dphidxia
     &                              + dxia*ddrdy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                              + d2edphi2*dphidzic*dphidxia
     &                              + dxia*ddrdz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                              + d2edphi2*dphidxic*dphidyia
     &                              + dyia*ddrdx
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                              + d2edphi2*dphidyic*dphidyia
     &                              + dyia*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                              + d2edphi2*dphidzic*dphidyia
     &                              + dyia*ddrdz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                              + d2edphi2*dphidxic*dphidzia
     &                              + dzia*ddrdx
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                              + d2edphi2*dphidyic*dphidzia
     &                              + dzia*ddrdy
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                              + d2edphi2*dphidzic*dphidzia
     &                              + dzia*ddrdz
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                              + d2edphi2*dphidxic*dphidxib
     &                              - (dxic-dxib)*ddrdx - d2drdxx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                              + d2edphi2*dphidyic*dphidxib
     &                              - dyic*ddrdx + dxib*ddrdy - d2drdxy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                              + d2edphi2*dphidzic*dphidxib
     &                              - dzic*ddrdx + dxib*ddrdz - d2drdxz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                              + d2edphi2*dphidxic*dphidyib
     &                              - dxic*ddrdy + dyib*ddrdx - d2drdxy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                              + d2edphi2*dphidyic*dphidyib
     &                              - (dyic-dyib)*ddrdy - d2drdyy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                              + d2edphi2*dphidzic*dphidyib
     &                              - dzic*ddrdy + dyib*ddrdz - d2drdyz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                              + d2edphi2*dphidxic*dphidzib
     &                              - dxic*ddrdz + dzib*ddrdx - d2drdxz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                              + d2edphi2*dphidyic*dphidzib
     &                              - dyic*ddrdz + dzib*ddrdy - d2drdyz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                              + d2edphi2*dphidzic*dphidzib
     &                              - (dzic-dzib)*ddrdz - d2drdzz
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                              + d2edphi2*dphidxic*dphidxid
     &                              + dxid*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                              + d2edphi2*dphidyic*dphidxid
     &                              + dxid*ddrdy
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                              + d2edphi2*dphidzic*dphidxid
     &                              + dxid*ddrdz
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                              + d2edphi2*dphidxic*dphidyid
     &                              + dyid*ddrdx
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                              + d2edphi2*dphidyic*dphidyid
     &                              + dyid*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                              + d2edphi2*dphidzic*dphidyid
     &                              + dyid*ddrdz
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                              + d2edphi2*dphidxic*dphidzid
     &                              + dzid*ddrdx
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                              + d2edphi2*dphidyic*dphidzid
     &                              + dzid*ddrdy
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                              + d2edphi2*dphidzic*dphidzid
     &                              + dzid*ddrdz
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                              + d2edphi2*dphidxid*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                              + d2edphi2*dphidyid*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                              + d2edphi2*dphidzid*dphidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxid*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                              + d2edphi2*dphidyid*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                              + d2edphi2*dphidzid*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                              + d2edphi2*dphidxid*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                              + d2edphi2*dphidyid*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                              + d2edphi2*dphidzid*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                              + d2edphi2*dphidxid*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                              + d2edphi2*dphidyid*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                              + d2edphi2*dphidzid*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                              + d2edphi2*dphidxid*dphidxib
     &                              - dxid*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                              + d2edphi2*dphidyid*dphidxib
     &                              - dyid*ddrdx
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                              + d2edphi2*dphidzid*dphidxib
     &                              - dzid*ddrdx
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                              + d2edphi2*dphidxid*dphidyib
     &                              - dxid*ddrdy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                              + d2edphi2*dphidyid*dphidyib
     &                              - dyid*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                              + d2edphi2*dphidzid*dphidyib
     &                              - dzid*ddrdy
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                              + d2edphi2*dphidxid*dphidzib
     &                              - dxid*ddrdz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                              + d2edphi2*dphidyid*dphidzib
     &                              - dyid*ddrdz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                              + d2edphi2*dphidzid*dphidzib
     &                              - dzid*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                              + d2edphi2*dphidxid*dphidxic
     &                              + dxid*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                              + d2edphi2*dphidyid*dphidxic
     &                              + dyid*ddrdx
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                              + d2edphi2*dphidzid*dphidxic
     &                              + dzid*ddrdx
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                              + d2edphi2*dphidxid*dphidyic
     &                              + dxid*ddrdy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                              + d2edphi2*dphidyid*dphidyic
     &                              + dyid*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                              + d2edphi2*dphidzid*dphidyic
     &                              + dzid*ddrdy
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                              + d2edphi2*dphidxid*dphidzic
     &                              + dxid*ddrdz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                              + d2edphi2*dphidyid*dphidzic
     &                              + dyid*ddrdz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                              + d2edphi2*dphidzid*dphidzic
     &                              + dzid*ddrdz
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine estrtor3  --  stretch-torsion energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "estrtor3" calculates the stretch-torsion potential energy;
c     also partitions the energy terms among the atoms
c
c
      subroutine estrtor3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'strtor.i'
      include 'torpot.i'
      include 'tors.i'
      include 'usage.i'
      integer i,k,ia,ib,ic,id,istrtor
      real*8 e,dr,angle,fgrp
      real*8 rcb,rt2,ru2,rtru
      real*8 xt,yt,zt,xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3,s1,s2,s3
      real*8 cosine,cosine2,cosine3
      real*8 sine,sine2,sine3
      real*8 phi1,phi2,phi3
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      logical header,huge,proceed
c
c
c     zero out the stretch-torsion energy and partitioning terms
c
      nebt = 0
      ebt = 0.0d0
      do i = 1, n
         aebt(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the stretch-torsion interaction energy term
c
      do istrtor = 1, nstrtor
         i = ist(1,istrtor)
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     set the stretch-torsional parameters for this angle
c
               v1 = kst(1,istrtor)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = kst(2,istrtor)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = kst(3,istrtor)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
c
c     calculate the bond-stretch for the central bond
c
               k = ist(2,istrtor)
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               dr = rcb - bl(k)
c
c     compute the stretch-torsion energy for this angle
c
               e = storunit * dr * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total stretch-torsion energy
c
               nebt = nebt + 1
               ebt = ebt + e
               aebt(ib) = aebt(ib) + 0.5d0*e
               aebt(ic) = aebt(ic) + 0.5d0*e
c
c     print a warning if the energy of this angle is large
c
               huge = (abs(e) .gt. 1.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Stretch-Torsion Cross',
     &                          ' Term Interactions :',
     &                       //,' Type',21x,'Atom Names',20x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                             name(ic),id,name(id),angle,e
   20             format (' StrTors  ',i5,'-',a3,3(1x,i5,'-',a3),
     &                       2x,f10.4,f12.4)
               end if
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
c     ########################################################
c     ##                                                    ##
c     ##  subroutine etors  --  torsional potential energy  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "etors" calculates the torsional potential energies
c
c
      subroutine etors
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'torpot.i'
      include 'tors.i'
      include 'usage.i'
      include 'warp.i'
      integer i,ia,ib,ic,id
      real*8 e,fgrp,time
      real*8 xt,yt,zt,xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 v1,v2,v3,v4,v5,v6
      real*8 c1,c2,c3,c4,c5,c6
      real*8 s1,s2,s3,s4,s5,s6
      real*8 cosine,cosine2,cosine3
      real*8 cosine4,cosine5,cosine6
      real*8 sine,sine2,sine3
      real*8 sine4,sine5,sine6
      real*8 phi1,phi2,phi3
      real*8 phi4,phi5,phi6
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      logical proceed
c
c
c     zero out the torsional potential energy
c
      et = 0.0d0
c
c     scale deformation time by the diffusion coefficient
c
      if (use_deform)  time = difft * deform
c
c     calculate the torsional angle energy term
c
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the torsional parameters for this angle
c
               v1 = tors1(1,i)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = tors2(1,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = tors3(1,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               v4 = tors4(1,i)
               c4 = tors4(3,i)
               s4 = tors4(4,i)
               v5 = tors5(1,i)
               c5 = tors5(3,i)
               s5 = tors5(4,i)
               v6 = tors6(1,i)
               c6 = tors6(3,i)
               s6 = tors6(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               cosine4 = cosine*cosine3 - sine*sine3
               sine4 = cosine*sine3 + sine*cosine3
               cosine5 = cosine*cosine4 - sine*sine4
               sine5 = cosine*sine4 + sine*cosine4
               cosine6 = cosine*cosine5 - sine*sine5
               sine6 = cosine*sine5 + sine*cosine5
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               phi4 = 1.0d0 + (cosine4*c4 + sine4*s4)
               phi5 = 1.0d0 + (cosine5*c5 + sine5*s5)
               phi6 = 1.0d0 + (cosine6*c6 + sine6*s6)
c
c     transform the potential via diffusional smoothing
c
               if (use_deform) then
                  phi1 = phi1 * exp(-time)
                  phi2 = phi2 * exp(-4.0d0*time)
                  phi3 = phi3 * exp(-9.0d0*time)
                  phi4 = phi4 * exp(-16.0d0*time)
                  phi5 = phi5 * exp(-25.0d0*time)
                  phi6 = phi6 * exp(-36.0d0*time)
               end if
c
c     calculate the torsional energy for this angle
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3
     &                            + v4*phi4 + v5*phi5 + v6*phi6)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total torsional angle energy
c
               et = et + e
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine etors1  --  torsional energy & derivatives  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "etors1" calculates torsional potential energy and first
c     derivatives with respect to Cartesian coordinates
c
c
      subroutine etors1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'torpot.i'
      include 'tors.i'
      include 'usage.i'
      include 'virial.i'
      include 'warp.i'
      integer i,ia,ib,ic,id
      real*8 e,dedphi,fgrp,time
      real*8 v1,v2,v3,v4,v5,v6
      real*8 c1,c2,c3,c4,c5,c6
      real*8 s1,s2,s3,s4,s5,s6
      real*8 cosine,cosine2,cosine3
      real*8 cosine4,cosine5,cosine6
      real*8 sine,sine2,sine3
      real*8 sine4,sine5,sine6
      real*8 phi1,phi2,phi3
      real*8 phi4,phi5,phi6
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 dphi1,dphi2,dphi3
      real*8 dphi4,dphi5,dphi6
      real*8 deform1,deform2,deform3
      real*8 deform4,deform5,deform6
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      logical proceed
c
c
c     zero out energy and first derivative components
c
      et = 0.0d0
      do i = 1, n
         det(1,i) = 0.0d0
         det(2,i) = 0.0d0
         det(3,i) = 0.0d0
      end do
c
c     scale deformation time by the diffusion coefficient
c
      if (use_deform)  time = difft * deform
c
c     calculate the torsional angle energy term
c
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the torsional parameters for this angle
c
               v1 = tors1(1,i)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = tors2(1,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = tors3(1,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               v4 = tors4(1,i)
               c4 = tors4(3,i)
               s4 = tors4(4,i)
               v5 = tors5(1,i)
               c5 = tors5(3,i)
               s5 = tors5(4,i)
               v6 = tors6(1,i)
               c6 = tors6(3,i)
               s6 = tors6(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               cosine4 = cosine*cosine3 - sine*sine3
               sine4 = cosine*sine3 + sine*cosine3
               cosine5 = cosine*cosine4 - sine*sine4
               sine5 = cosine*sine4 + sine*cosine4
               cosine6 = cosine*cosine5 - sine*sine5
               sine6 = cosine*sine5 + sine*cosine5
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               phi4 = 1.0d0 + (cosine4*c4 + sine4*s4)
               phi5 = 1.0d0 + (cosine5*c5 + sine5*s5)
               phi6 = 1.0d0 + (cosine6*c6 + sine6*s6)
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
               dphi4 = 4.0d0 * (cosine4*s4 - sine4*c4)
               dphi5 = 5.0d0 * (cosine5*s5 - sine5*c5)
               dphi6 = 6.0d0 * (cosine6*s6 - sine6*c6)
c
c     transform the potential via diffusional smoothing
c
               if (use_deform) then
                  deform1 = exp(-time)
                  deform2 = exp(-4.0d0*time)
                  deform3 = exp(-9.0d0*time)
                  deform4 = exp(-16.0d0*time)
                  deform5 = exp(-25.0d0*time)
                  deform6 = exp(-36.0d0*time)
                  phi1 = phi1 * deform1
                  phi2 = phi2 * deform2
                  phi3 = phi3 * deform3
                  phi4 = phi4 * deform4
                  phi5 = phi5 * deform5
                  phi6 = phi6 * deform6
                  dphi1 = dphi1 * deform1
                  dphi2 = dphi2 * deform2
                  dphi3 = dphi3 * deform3
                  dphi4 = dphi4 * deform4
                  dphi5 = dphi5 * deform5
                  dphi6 = dphi6 * deform6
               end if
c
c     calculate torsional energy and master chain rule term
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3
     &                            + v4*phi4 + v5*phi5 + v6*phi6)
               dedphi = torsunit * (v1*dphi1 + v2*dphi2 + v3*dphi3
     &                                 + v4*dphi4 + v5*dphi5 + v6*dphi6)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedphi = dedphi * fgrp
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = yca*dphidzt - zca*dphidyt
     &                       + zdc*dphidyu - ydc*dphidzu
               dphidyib = zca*dphidxt - xca*dphidzt
     &                       + xdc*dphidzu - zdc*dphidxu
               dphidzib = xca*dphidyt - yca*dphidxt
     &                       + ydc*dphidxu - xdc*dphidyu
               dphidxic = zba*dphidyt - yba*dphidzt
     &                       + ydb*dphidzu - zdb*dphidyu
               dphidyic = xba*dphidzt - zba*dphidxt
     &                       + zdb*dphidxu - xdb*dphidzu
               dphidzic = yba*dphidxt - xba*dphidyt
     &                       + xdb*dphidyu - ydb*dphidxu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     compute first derivative components for this angle
c
               dedxia = dedphi * dphidxia
               dedyia = dedphi * dphidyia
               dedzia = dedphi * dphidzia
               dedxib = dedphi * dphidxib
               dedyib = dedphi * dphidyib
               dedzib = dedphi * dphidzib
               dedxic = dedphi * dphidxic
               dedyic = dedphi * dphidyic
               dedzic = dedphi * dphidzic
               dedxid = dedphi * dphidxid
               dedyid = dedphi * dphidyid
               dedzid = dedphi * dphidzid
c
c     increment the total torsional angle energy and gradient
c
               et = et + e
               det(1,ia) = det(1,ia) + dedxia
               det(2,ia) = det(2,ia) + dedyia
               det(3,ia) = det(3,ia) + dedzia
               det(1,ib) = det(1,ib) + dedxib
               det(2,ib) = det(2,ib) + dedyib
               det(3,ib) = det(3,ib) + dedzib
               det(1,ic) = det(1,ic) + dedxic
               det(2,ic) = det(2,ic) + dedyic
               det(3,ic) = det(3,ic) + dedzic
               det(1,id) = det(1,id) + dedxid
               det(2,id) = det(2,id) + dedyid
               det(3,id) = det(3,id) + dedzid
c
c     increment the components of the virial
c
               if (isobaric) then
                  virx = virx - xba*dedxia + xdc*dedxid
     &                      + xcb*(dedxic+dedxid)
                  viry = viry - yba*dedyia + ydc*dedyid
     &                      + ycb*(dedyic+dedyid)
                  virz = virz - zba*dedzia + zdc*dedzid
     &                      + zcb*(dedzic+dedzid)
               end if
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine etors2  --  atom-by-atom torsional Hessian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "etors2" calculates second derivatives of the torsional
c     energy for a single atom
c
c
      subroutine etors2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      include 'hessn.i'
      include 'torpot.i'
      include 'tors.i'
      include 'warp.i'
      integer i,ia,ib,ic,id,ktors
      real*8 dedphi,d2edphi2
      real*8 fgrp,time
      real*8 v1,v2,v3,v4,v5,v6
      real*8 c1,c2,c3,c4,c5,c6
      real*8 s1,s2,s3,s4,s5,s6
      real*8 cosine,cosine2,cosine3
      real*8 cosine4,cosine5,cosine6
      real*8 sine,sine2,sine3
      real*8 sine4,sine5,sine6
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 dphi1,dphi2,dphi3
      real*8 dphi4,dphi5,dphi6
      real*8 d2phi1,d2phi2,d2phi3
      real*8 d2phi4,d2phi5,d2phi6
      real*8 deform1,deform2,deform3
      real*8 deform4,deform5,deform6
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 dphidxibt,dphidyibt,dphidzibt
      real*8 dphidxibu,dphidyibu,dphidzibu
      real*8 dphidxict,dphidyict,dphidzict
      real*8 dphidxicu,dphidyicu,dphidzicu
      real*8 dxiaxia,dyiayia,dziazia,dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic,dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia,dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic,dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib,dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib,dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic,dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid,dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid,dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic,dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid,dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid,dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid,dzicxid,dzicyid,dziczid
      logical proceed
c
c
c     scale deformation time by the diffusion coefficient
c
c        next line added by Yingbin in 11/2006 to prevent crash
c        on Sun Opteron/Solaris (conditional execution below?).
c        USE_DEFORM is expected to be .FALSE. always.
      time = 0.01d+00
c
      if (use_deform)  time = difft * deform
c
c     calculate the torsional angle energy term
c
      do ktors = 1, ntors
         ia = itors(1,ktors)
         ib = itors(2,ktors)
         ic = itors(3,ktors)
         id = itors(4,ktors)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (i.eq.ia .or. i.eq.ib .or.
     &                              i.eq.ic .or. i.eq.id)
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the torsional parameters for this angle
c
               v1 = tors1(1,ktors)
               c1 = tors1(3,ktors)
               s1 = tors1(4,ktors)
               v2 = tors2(1,ktors)
               c2 = tors2(3,ktors)
               s2 = tors2(4,ktors)
               v3 = tors3(1,ktors)
               c3 = tors3(3,ktors)
               s3 = tors3(4,ktors)
               v4 = tors4(1,ktors)
               c4 = tors4(3,ktors)
               s4 = tors4(4,ktors)
               v5 = tors5(1,ktors)
               c5 = tors5(3,ktors)
               s5 = tors5(4,ktors)
               v6 = tors6(1,ktors)
               c6 = tors6(3,ktors)
               s6 = tors6(4,ktors)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               cosine4 = cosine*cosine3 - sine*sine3
               sine4 = cosine*sine3 + sine*cosine3
               cosine5 = cosine*cosine4 - sine*sine4
               sine5 = cosine*sine4 + sine*cosine4
               cosine6 = cosine*cosine5 - sine*sine5
               sine6 = cosine*sine5 + sine*cosine5
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
               dphi4 = 4.0d0 * (cosine4*s4 - sine4*c4)
               dphi5 = 5.0d0 * (cosine5*s5 - sine5*c5)
               dphi6 = 6.0d0 * (cosine6*s6 - sine6*c6)
               d2phi1 = -(cosine*c1 + sine*s1)
               d2phi2 = -4.0d0 * (cosine2*c2 + sine2*s2)
               d2phi3 = -9.0d0 * (cosine3*c3 + sine3*s3)
               d2phi4 = -16.0d0 * (cosine4*c4 + sine4*s4)
               d2phi5 = -25.0d0 * (cosine5*c5 + sine5*s5)
               d2phi6 = -36.0d0 * (cosine6*c6 + sine6*s6)
c
c     transform the potential via diffusional smoothing
c
               if (use_deform) then
                  deform1 = exp(-time)
                  deform2 = exp(-4.0d0*time)
                  deform3 = exp(-9.0d0*time)
                  deform4 = exp(-16.0d0*time)
                  deform5 = exp(-25.0d0*time)
                  deform6 = exp(-36.0d0*time)
                  dphi1 = dphi1 * deform1
                  dphi2 = dphi2 * deform2
                  dphi3 = dphi3 * deform3
                  dphi4 = dphi4 * deform4
                  dphi5 = dphi5 * deform5
                  dphi6 = dphi6 * deform6
                  d2phi1 = d2phi1 * deform1
                  d2phi2 = d2phi2 * deform2
                  d2phi3 = d2phi3 * deform3
                  d2phi4 = d2phi4 * deform4
                  d2phi5 = d2phi5 * deform5
                  d2phi6 = d2phi6 * deform6
               end if
c
c     calculate the torsional master chain rule terms
c
               dedphi = torsunit * (v1*dphi1 + v2*dphi2 + v3*dphi3
     &                            + v4*dphi4 + v5*dphi5 + v6*dphi6)
               d2edphi2 = torsunit * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3
     &                              + v4*d2phi4 + v5*d2phi5 + v6*d2phi6)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
               xycb2 = xcb*xcb + ycb*ycb
               xzcb2 = xcb*xcb + zcb*zcb
               yzcb2 = ycb*ycb + zcb*zcb
               rcbxt = -2.0d0 * rcb * dphidxt
               rcbyt = -2.0d0 * rcb * dphidyt
               rcbzt = -2.0d0 * rcb * dphidzt
               rcbt2 = rcb * rt2
               rcbxu = 2.0d0 * rcb * dphidxu
               rcbyu = 2.0d0 * rcb * dphidyu
               rcbzu = 2.0d0 * rcb * dphidzu
               rcbu2 = rcb * ru2
               dphidxibt = yca*dphidzt - zca*dphidyt
               dphidxibu = zdc*dphidyu - ydc*dphidzu
               dphidyibt = zca*dphidxt - xca*dphidzt
               dphidyibu = xdc*dphidzu - zdc*dphidxu
               dphidzibt = xca*dphidyt - yca*dphidxt
               dphidzibu = ydc*dphidxu - xdc*dphidyu
               dphidxict = zba*dphidyt - yba*dphidzt
               dphidxicu = ydb*dphidzu - zdb*dphidyu
               dphidyict = xba*dphidzt - zba*dphidxt
               dphidyicu = zdb*dphidxu - xdb*dphidzu
               dphidzict = yba*dphidxt - xba*dphidyt
               dphidzicu = xdb*dphidyu - ydb*dphidxu
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = dphidxibt + dphidxibu
               dphidyib = dphidyibt + dphidyibu
               dphidzib = dphidzibt + dphidzibu
               dphidxic = dphidxict + dphidxicu
               dphidyic = dphidyict + dphidyicu
               dphidzic = dphidzict + dphidzicu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     chain rule terms for second derivative components
c
               dxiaxia = rcbxt*dphidxia
               dxiayia = rcbxt*dphidyia - zcb*rcb/rt2
               dxiazia = rcbxt*dphidzia + ycb*rcb/rt2
               dxiaxib = rcbxt*dphidxibt + xcb*(zca*ycb-yca*zcb)/rcbt2
               dxiayib = rcbxt*dphidyibt + dphidzt
     &                      + (xca*zcb*xcb+zca*yzcb2)/rcbt2
               dxiazib = rcbxt*dphidzibt - dphidyt
     &                      - (xca*ycb*xcb+yca*yzcb2)/rcbt2
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2
               dxiayic = rcbxt*dphidyict - dphidzt
     &                      - (xba*zcb*xcb+zba*yzcb2)/rcbt2
               dxiazic = rcbxt*dphidzict + dphidyt
     &                      + (xba*ycb*xcb+yba*yzcb2)/rcbt2
               dxiaxid = 0.0d0
               dxiayid = 0.0d0
               dxiazid = 0.0d0
               dyiayia = rcbyt*dphidyia
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2
               dyiaxib = rcbyt*dphidxibt - dphidzt
     &                      - (yca*zcb*ycb+zca*xzcb2)/rcbt2
               dyiayib = rcbyt*dphidyibt + ycb*(xca*zcb-zca*xcb)/rcbt2
               dyiazib = rcbyt*dphidzibt + dphidxt
     &                      + (yca*xcb*ycb+xca*xzcb2)/rcbt2
               dyiaxic = rcbyt*dphidxict + dphidzt
     &                      + (yba*zcb*ycb+zba*xzcb2)/rcbt2
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2
               dyiazic = rcbyt*dphidzict - dphidxt
     &                      - (yba*xcb*ycb+xba*xzcb2)/rcbt2
               dyiaxid = 0.0d0
               dyiayid = 0.0d0
               dyiazid = 0.0d0
               dziazia = rcbzt*dphidzia
               dziaxib = rcbzt*dphidxibt + dphidyt
     &                      + (zca*ycb*zcb+yca*xycb2)/rcbt2
               dziayib = rcbzt*dphidyibt - dphidxt
     &                      - (zca*xcb*zcb+xca*xycb2)/rcbt2
               dziazib = rcbzt*dphidzibt + zcb*(yca*xcb-xca*ycb)/rcbt2
               dziaxic = rcbzt*dphidxict - dphidyt
     &                      - (zba*ycb*zcb+yba*xycb2)/rcbt2
               dziayic = rcbzt*dphidyict + dphidxt
     &                      + (zba*xcb*zcb+xba*xycb2)/rcbt2
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2
               dziaxid = 0.0d0
               dziayid = 0.0d0
               dziazid = 0.0d0
               dxibxic = -xcb*dphidxib/(rcb*rcb)
     &             - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidxibt/rt2
     &             - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidxibu/ru2
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
     &             - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidxibt/rt2
     &             + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidxibu/ru2
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2
               dxibyid = rcbyu*dphidxibu - dphidzu
     &                      - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
               dxibzid = rcbzu*dphidxibu + dphidyu
     &                      + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
               dyibzib = ycb*dphidzib/(rcb*rcb)
     &             - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &             - 2.0d0*(xt*zca-xca*zt)*dphidzibt/rt2
     &             + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &             + 2.0d0*(xu*zdc-xdc*zu)*dphidzibu/ru2
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu
     &             + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidyibt/rt2
     &             - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidyibu/ru2
               dyibyic = -ycb*dphidyib/(rcb*rcb)
     &             - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidyibt/rt2
     &             - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidyibu/ru2
               dyibxid = rcbxu*dphidyibu + dphidzu
     &                      + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2
               dyibzid = rcbzu*dphidyibu - dphidxu
     &                      - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu
     &             - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidzibt/rt2
     &             + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidzibu/ru2
               dzibzic = -zcb*dphidzib/(rcb*rcb)
     &             - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &             - 2.0d0*(xt*yba-xba*yt)*dphidzibt/rt2
     &             - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &             + 2.0d0*(xu*ydb-xdb*yu)*dphidzibu/ru2
               dzibxid = rcbxu*dphidzibu - dphidyu
     &                      - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
               dzibyid = rcbyu*dphidzibu + dphidxu
     &                      + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
               dxicyid = rcbyu*dphidxicu + dphidzu
     &                      + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
               dxiczid = rcbzu*dphidxicu - dphidyu
     &                      - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
               dyicxid = rcbxu*dphidyicu - dphidzu
     &                      - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
               dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
               dyiczid = rcbzu*dphidyicu + dphidxu
     &                      + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
               dzicxid = rcbxu*dphidzicu + dphidyu
     &                      + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
               dzicyid = rcbyu*dphidzicu - dphidxu
     &                      - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
               dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
               dxidxid = rcbxu*dphidxid
               dxidyid = rcbxu*dphidyid + zcb*rcb/ru2
               dxidzid = rcbxu*dphidzid - ycb*rcb/ru2
               dyidyid = rcbyu*dphidyid
               dyidzid = rcbyu*dphidzid + xcb*rcb/ru2
               dzidzid = rcbzu*dphidzid
c
c     get some second derivative chain rule terms by difference
c
               dxibxib = -dxiaxib - dxibxic - dxibxid
               dxibyib = -dyiaxib - dxibyic - dxibyid
               dxibzib = -dxiazib - dzibxic - dzibxid
               dxibzic = -dziaxib - dxibzib - dxibzid
               dyibyib = -dyiayib - dyibyic - dyibyid
               dyibzic = -dziayib - dyibzib - dyibzid
               dzibzib = -dziazib - dzibzic - dzibzid
               dzibyic = -dyiazib - dyibzib - dzibyid
               dxicxic = -dxiaxic - dxibxic - dxicxid
               dxicyic = -dyiaxic - dyibxic - dxicyid
               dxiczic = -dziaxic - dzibxic - dxiczid
               dyicyic = -dyiayic - dyibyic - dyicyid
               dyiczic = -dziayic - dzibyic - dyiczid
               dziczic = -dziazic - dzibzic - dziczid
c
c     now, increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                              + d2edphi2*dphidxia*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                              + d2edphi2*dphidyia*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                              + d2edphi2*dphidzia*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxia*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                              + d2edphi2*dphidyia*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                              + d2edphi2*dphidzia*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                              + d2edphi2*dphidxia*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                              + d2edphi2*dphidyia*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                              + d2edphi2*dphidzia*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                              + d2edphi2*dphidxia*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                              + d2edphi2*dphidyia*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                              + d2edphi2*dphidzia*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxia*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                              + d2edphi2*dphidyia*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                              + d2edphi2*dphidzia*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                              + d2edphi2*dphidxia*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                              + d2edphi2*dphidyia*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                              + d2edphi2*dphidzia*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                              + d2edphi2*dphidxia*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                              + d2edphi2*dphidyia*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                              + d2edphi2*dphidzia*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxia*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                              + d2edphi2*dphidyia*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                              + d2edphi2*dphidzia*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                              + d2edphi2*dphidxia*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                              + d2edphi2*dphidyia*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                              + d2edphi2*dphidzia*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                              + d2edphi2*dphidxia*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                              + d2edphi2*dphidyia*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                              + d2edphi2*dphidzia*dphidzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                              + d2edphi2*dphidxib*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                              + d2edphi2*dphidyib*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                              + d2edphi2*dphidzib*dphidzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxib*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                              + d2edphi2*dphidyib*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                              + d2edphi2*dphidzib*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                              + d2edphi2*dphidxib*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                              + d2edphi2*dphidyib*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                              + d2edphi2*dphidzib*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                              + d2edphi2*dphidxib*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                              + d2edphi2*dphidyib*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                              + d2edphi2*dphidzib*dphidzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                              + d2edphi2*dphidxib*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                              + d2edphi2*dphidyib*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                              + d2edphi2*dphidzib*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                              + d2edphi2*dphidxib*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                              + d2edphi2*dphidyib*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                              + d2edphi2*dphidzib*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                              + d2edphi2*dphidxib*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                              + d2edphi2*dphidyib*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                              + d2edphi2*dphidzib*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                              + d2edphi2*dphidxib*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                              + d2edphi2*dphidyib*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                              + d2edphi2*dphidzib*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                              + d2edphi2*dphidxib*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                              + d2edphi2*dphidyib*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                              + d2edphi2*dphidzib*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                              + d2edphi2*dphidxib*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                              + d2edphi2*dphidyib*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                              + d2edphi2*dphidzib*dphidzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                              + d2edphi2*dphidxic*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                              + d2edphi2*dphidyic*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                              + d2edphi2*dphidzic*dphidzic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxic*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                              + d2edphi2*dphidyic*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                              + d2edphi2*dphidzic*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                              + d2edphi2*dphidxic*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                              + d2edphi2*dphidyic*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                              + d2edphi2*dphidzic*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                              + d2edphi2*dphidxic*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                              + d2edphi2*dphidyic*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                              + d2edphi2*dphidzic*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                              + d2edphi2*dphidxic*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                              + d2edphi2*dphidyic*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                              + d2edphi2*dphidzic*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                              + d2edphi2*dphidxic*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                              + d2edphi2*dphidyic*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                              + d2edphi2*dphidzic*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                              + d2edphi2*dphidxic*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                              + d2edphi2*dphidyic*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                              + d2edphi2*dphidzic*dphidzib
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                              + d2edphi2*dphidxic*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                              + d2edphi2*dphidyic*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                              + d2edphi2*dphidzic*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                              + d2edphi2*dphidxic*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                              + d2edphi2*dphidyic*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                              + d2edphi2*dphidzic*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                              + d2edphi2*dphidxic*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                              + d2edphi2*dphidyic*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                              + d2edphi2*dphidzic*dphidzid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                              + d2edphi2*dphidxid*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                              + d2edphi2*dphidyid*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                              + d2edphi2*dphidzid*dphidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxid*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                              + d2edphi2*dphidyid*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                              + d2edphi2*dphidzid*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                              + d2edphi2*dphidxid*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                              + d2edphi2*dphidyid*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                              + d2edphi2*dphidzid*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                              + d2edphi2*dphidxid*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                              + d2edphi2*dphidyid*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                              + d2edphi2*dphidzid*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                              + d2edphi2*dphidxid*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                              + d2edphi2*dphidyid*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                              + d2edphi2*dphidzid*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                              + d2edphi2*dphidxid*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                              + d2edphi2*dphidyid*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                              + d2edphi2*dphidzid*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                              + d2edphi2*dphidxid*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                              + d2edphi2*dphidyid*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                              + d2edphi2*dphidzid*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                              + d2edphi2*dphidxid*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                              + d2edphi2*dphidyid*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                              + d2edphi2*dphidzid*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                              + d2edphi2*dphidxid*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                              + d2edphi2*dphidyid*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                              + d2edphi2*dphidzid*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                              + d2edphi2*dphidxid*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                              + d2edphi2*dphidyid*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                              + d2edphi2*dphidzid*dphidzic
               end if
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
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine etors3  --  torsional energy & analysis  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "etors3" calculates the torsional and stretch-torsion
c     potential energies; also partitions the energy terms
c     among the atoms
c
c
      subroutine etors3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'torpot.i'
      include 'tors.i'
      include 'usage.i'
      include 'warp.i'
      integer i,ia,ib,ic,id
      real*8 e,angle,fgrp,time
      real*8 xt,yt,zt,xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 v1,v2,v3,v4,v5,v6
      real*8 c1,c2,c3,c4,c5,c6
      real*8 s1,s2,s3,s4,s5,s6
      real*8 cosine,cosine2,cosine3
      real*8 cosine4,cosine5,cosine6
      real*8 sine,sine2,sine3
      real*8 sine4,sine5,sine6
      real*8 phi1,phi2,phi3
      real*8 phi4,phi5,phi6
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      logical header,huge,proceed
c
c
c     zero out the torsional energy and partitioning terms
c
      net = 0
      et = 0.0d0
      do i = 1, n
         aet(i) = 0.0d0
      end do
      header = .true.
c
c     scale deformation time by the diffusion coefficient
c
      if (use_deform)  time = difft * deform
c
c     calculate the torsional angle energy term
c
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     set the torsional parameters for this angle
c
               v1 = tors1(1,i)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = tors2(1,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = tors3(1,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               v4 = tors4(1,i)
               c4 = tors4(3,i)
               s4 = tors4(4,i)
               v5 = tors5(1,i)
               c5 = tors5(3,i)
               s5 = tors5(4,i)
               v6 = tors6(1,i)
               c6 = tors6(3,i)
               s6 = tors6(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               cosine4 = cosine*cosine3 - sine*sine3
               sine4 = cosine*sine3 + sine*cosine3
               cosine5 = cosine*cosine4 - sine*sine4
               sine5 = cosine*sine4 + sine*cosine4
               cosine6 = cosine*cosine5 - sine*sine5
               sine6 = cosine*sine5 + sine*cosine5
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               phi4 = 1.0d0 + (cosine4*c4 + sine4*s4)
               phi5 = 1.0d0 + (cosine5*c5 + sine5*s5)
               phi6 = 1.0d0 + (cosine6*c6 + sine6*s6)
c
c     transform the potential via diffusional smoothing
c
               if (use_deform) then
                  phi1 = phi1 * exp(-time)
                  phi2 = phi2 * exp(-4.0d0*time)
                  phi3 = phi3 * exp(-9.0d0*time)
                  phi4 = phi4 * exp(-16.0d0*time)
                  phi5 = phi5 * exp(-25.0d0*time)
                  phi6 = phi6 * exp(-36.0d0*time)
               end if
c
c     calculate the torsional energy for this angle
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3
     &                            + v4*phi4 + v5*phi5 + v6*phi6)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total torsional angle energy
c
               net = net + 1
               et = et + e
               aet(ib) = aet(ib) + 0.5d0*e
               aet(ic) = aet(ic) + 0.5d0*e
c
c     print a warning if the energy of this angle is large
c
               huge = (e .gt. 3.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Torsional Angle',
     &                          ' Interactions :',
     &                       //,' Type',21x,'Atom Names',20x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                             name(ic),id,name(id),angle,e
   20             format (' Torsion  ',i5,'-',a3,3(1x,i5,'-',a3),
     &                       2x,f10.4,f12.4)
               end if
            end if
         end if
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
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine eurey  --  Urey-Bradley potential energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "eurey" calculates the Urey-Bradley 1-3 interaction energy
c
c
      subroutine eurey
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'urey.i'
      include 'urypot.i'
      include 'usage.i'
      integer i,ia,ib
      real*8 e,ideal,force
      real*8 dt,dt2,fgrp
      real*8 xab,yab,zab,rab
      logical proceed
c
c
c     zero out the Urey-Bradley interaction energy
c
      eub = 0.0d0
c
c     calculate the Urey-Bradley 1-3 energy term
c
      do i = 1, nurey
         ia = iury(1,i)
         ib = iury(2,i)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,ia,ib,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
            dt2 = dt * dt
c
c     calculate the Urey-Bradley energy for this interaction
c
            e = ureyunit * force * dt2
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total Urey-Bradley energy
c
            eub = eub + e
         end if
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey1" calculates the Urey-Bradley interaction energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
      subroutine eurey1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'urey.i'
      include 'urypot.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ia,ib
      real*8 e,ideal,force
      real*8 dt,dt2,deddt,fgrp
      real*8 de,dedx,dedy,dedz
      real*8 xab,yab,zab,rab
      logical proceed
c
c
c     zero out the Urey-Bradley energy and first derivatives
c
      eub = 0.0d0
      do i = 1, n
         deub(1,i) = 0.0d0
         deub(2,i) = 0.0d0
         deub(3,i) = 0.0d0
      end do
c
c     calculate the Urey-Bradley 1-3 energy and first derivatives
c
      do i = 1, nurey
         ia = iury(1,i)
         ib = iury(2,i)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,ia,ib,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
            dt2 = dt * dt
            e = ureyunit * force * dt2
            deddt = 2.0d0 * ureyunit * force * dt
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               deddt = deddt * fgrp
            end if
c
c     compute chain rule terms needed for derivatives
c
            de = deddt / rab
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total Urey-Bradley energy and first derivatives
c
            eub = eub + e
            deub(1,ia) = deub(1,ia) + dedx
            deub(2,ia) = deub(2,ia) + dedy
            deub(3,ia) = deub(3,ia) + dedz
            deub(1,ib) = deub(1,ib) - dedx
            deub(2,ib) = deub(2,ib) - dedy
            deub(3,ib) = deub(3,ib) - dedz
c
c     increment the virial for use in pressure computation
c
            if (isobaric) then
               virx = virx + xab*dedx
               viry = viry + yab*dedy
               virz = virz + zab*dedz
            end if
         end if
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey2  --  atom-by-atom Urey-Bradley Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey2" calculates second derivatives of the Urey-Bradley
c     interaction energy for a single atom at a time
c
c
      subroutine eurey2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'urey.i'
      include 'urypot.i'
      integer i,j,ia,ib,iurey
      real*8 ideal,force,fgrp
      real*8 xab,yab,zab,rab,rab2
      real*8 dt,dt2,deddt,d2eddt2
      real*8 term,de,d2e(3,3)
      real*8 termx,termy,termz
      logical proceed
c
c
c     compute the Hessian elements of the Urey-Bradley energy
c
      do iurey = 1, nurey
         ia = iury(1,iurey)
         ib = iury(2,iurey)
         ideal = ul(iurey)
         force = uk(iurey)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,ia,ib,0,0,0)
         if (proceed)  proceed = (i.eq.ia .or. i.eq.ib)
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            if (i .eq. ib) then
               ib = ia
               ia = i
            end if
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            rab2 = xab*xab + yab*yab + zab*zab
            rab = sqrt(rab2)
            dt = rab - ideal
            dt2 = dt * dt
            deddt = 2.0d0 * ureyunit * force * dt
            d2eddt2 = 2.0d0 * ureyunit * force
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               deddt = deddt * fgrp
               d2eddt2 = d2eddt2 * fgrp
            end if
c
c     set the chain rule terms for the Hessian elements
c
            de = deddt / rab
            term = (d2eddt2-de) / rab2
            termx = term * xab
            termy = term * yab
            termz = term * zab
            d2e(1,1) = termx*xab + de
            d2e(1,2) = termx*yab
            d2e(1,3) = termx*zab
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yab + de
            d2e(2,3) = termy*zab
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + d2e(1,j)
               hessy(j,ia) = hessy(j,ia) + d2e(2,j)
               hessz(j,ia) = hessz(j,ia) + d2e(3,j)
               hessx(j,ib) = hessx(j,ib) - d2e(1,j)
               hessy(j,ib) = hessy(j,ib) - d2e(2,j)
               hessz(j,ib) = hessz(j,ib) - d2e(3,j)
            end do
         end if
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine eurey3  --  Urey-Bradley energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "eurey3" calculates the Urey-Bradley energy; also
c     partitions the energy among the atoms
c
c
      subroutine eurey3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'urey.i'
      include 'urypot.i'
      include 'usage.i'
      integer i,ia,ib
      real*8 e,ideal,force
      real*8 dt,dt2,fgrp
      real*8 xab,yab,zab,rab
      logical header,huge,proceed
c
c
c     zero out the Urey-Bradley energy and partitioning terms
c
      neub = 0
      eub = 0.0d0
      do i = 1, n
         aeub(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the Urey-Bradley 1-3 energy term
c
      do i = 1, nurey
         ia = iury(1,i)
         ib = iury(2,i)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,ia,ib,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
            dt2 = dt * dt
c
c     calculate the Urey-Bradley energy for this interaction
c
            e = ureyunit * force * dt2
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total Urey-Bradley energy
c
            neub = neub + 1
            eub = eub + e
            aeub(ia) = aeub(ia) + 0.5d0*e
            aeub(ib) = aeub(ib) + 0.5d0*e
c
c     print a warning if the energy of this interaction is large
c
            huge = (e .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Urey-Bradley Interactions :',
     &                    //,' Type',11x,'Atom Names',20x,'Ideal',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ideal,rab,e
   20          format (' UreyBrad ',i5,'-',a3,1x,i5,'-',a3,
     &                    12x,2f10.4,f12.4)
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
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine extra  --  user defined extra potentials  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "extra" calculates any additional user defined potential
c     energy contribution
c
c
cjrs
      subroutine extrat
cjrs
      implicit none
      include 'energi.i'
c
c
c     zero out the energy due to extra potential terms
c
      ex = 0.0d0
c
c     add any user-defined extra potentials below here
c
c     e = ......
c     ex = ex + e
c
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
c     ##  subroutine extra1  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra1" calculates any additional user defined potential
c     energy contribution and its first derivatives
c
c
      subroutine extra1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'energi.i'
      integer i
c
c
c     zero out the extra energy term and first derivatives
c
      ex = 0.0d0
      do i = 1, n
         dex(1,i) = 0.0d0
         dex(2,i) = 0.0d0
         dex(3,i) = 0.0d0
      end do
c
c     add any user-defined extra potentials and derivatives;
c     also increment intermolecular energy and virial as needed
c
c     e = ......
c     ex = ex + e
c     do i = 1, n
c        dex(1,i) = ......
c        dex(2,i) = ......
c        dex(3,i) = ......
c     end do
c
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
c     ##  subroutine extra2  --  atom-wise user defined Hessian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "extra2" calculates second derivatives of any additional
c     user defined potential energy contribution for a single
c     atom at a time
c
c
      subroutine extra2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'hessn.i'
      integer i
c
c
c     compute the Hessian elements for extra energy terms
c
c     do j = 1, n
c        hessx(1,j) = hessx(1,j) + ......
c        hessy(2,j) = hessy(2,j) + ......
c        hessz(3,j) = hessz(3,j) + ......
c     end do
c
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
c     ##  subroutine extra3  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "extra3" calculates any additional user defined potential
c     contribution and also partitions the energy among the atoms
c
c
      subroutine extra3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'energi.i'
      integer i
c
c
c     zero out the energy due to extra potential terms
c
      nex = 0
      ex = 0.0d0
      do i = 1, n
         aex(i) = 0.0d0
      end do
c
c     add any user-defined extra potentials and derivatives
c
c     e = ......
c     nex = nex + 1
c     ex = ex + e
c     aex(i) = ......
c
      return
      end
