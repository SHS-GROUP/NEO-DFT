c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  kstbnd.i  --  forcefield parameters for stretch-bending  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     stbn     stretch-bending parameters for each atom class
c
c
      real*8 stbn
      common /kstbnd/ stbn(3,maxclass)
