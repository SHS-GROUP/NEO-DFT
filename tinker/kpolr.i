c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  kpolr.i  --  forcefield parameters for polarizability  ##
c     ##                                                         ##
c     #############################################################
c
c
c     polr    dipole polarizability parameters for each atom type
c
c
      real*8 polr
      common /kpolr/ polr(maxtyp)
