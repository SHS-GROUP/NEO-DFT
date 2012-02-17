c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  usage.i  --  atoms active during energy computation  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     nuse    number of active atoms used in energy calculation
c     use     true if an atom is active, false if inactive
c
c
      integer nuse
      logical use
      common /usage/ nuse,use(maxatm)
