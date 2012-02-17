c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  polpot.i  --  specifics of polarization functional form  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     poleps    induced dipole convergence criterion (rms Debyes/atom)
c     pradius   radius of an idealized atom with unit polarizability
c     pgamma    prefactor in exponential polarization damping term
c     poltyp    type of polarization potential (direct or mutual)
c
c
      real*8 poleps,pradius,pgamma
      character*6 poltyp
      common /polpot/ poleps,pradius,pgamma,poltyp
