c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  rotate.i  --  molecule partitions for rotation of a bond  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nrot     total number of atoms moving when bond rotates
c     rot      atom numbers of atoms moving when bond rotates
c
c
      integer nrot,rot
      common /rotatet/ nrot,rot(maxatm)
