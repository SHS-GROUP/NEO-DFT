c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  minima.i  --  general parameters for minimizations  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     fctmin    value below which function is deemed optimized
c     hguess    initial value for the H-matrix diagonal elements
c     maxiter   maximum number of iterations during optimization
c     nextiter  iteration number to use for the first iteration
c     iprint    iterations between status printing (0=no printing)
c     iwrite    iterations between coordinate dumps (0=no dumps)
c
c
      integer maxiter,nextiter
      integer iprint,iwrite
      real*8 fctmin,hguess
      common /minima/ fctmin,hguess,maxiter,nextiter,iprint,iwrite
