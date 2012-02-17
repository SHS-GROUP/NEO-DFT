c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  polar.i  --  polarizabilities and induced dipole moments  ##
c     ##                                                            ##
c     ################################################################
c
c
c     polarize  dipole polarizability for each multipole site (Ang**3)
c     pdamp     value of polarizability damping factor for each site
c     uind      induced dipole components at each multipole site
c     npolar    total number of polarizable sites in the system
c
c
      integer npolar
      real*8 polarize,pdamp,uind
      common /polar/ polarize(maxatm),pdamp(maxatm),uind(3,maxatm),
     &               npolar
