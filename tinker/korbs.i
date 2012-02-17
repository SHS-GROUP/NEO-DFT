c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  korbs.i  --  forcefield parameters for pi-system orbitals  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnpi    maximum number of pi-system bond parameters
c
c     electron  number of pi-electrons for each atom class
c     ionize    ionization potential for each atom class
c     repulse   repulsion integral value for each atom class
c     sslope    slope for bond stretch vs. pi-bond order
c     tslope    slope for 2-fold torsion vs. pi-bond order
c     kpi       string of atom classes for pi-system bond parameters
c
c
      integer maxnpi
      parameter (maxnpi=100)
      real*8 electron,ionize,repulse,sslope,tslope
      character*6 kpi
      common /korbs/ electron(maxclass),ionize(maxclass),
     &               repulse(maxclass),sslope(maxnpi),tslope(maxnpi),
     &               kpi(maxnpi)
