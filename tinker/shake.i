c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  shake.i  --  definition of Shake/Rattle constraints  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     krat         ideal distance value for rattle constraint
c     nrat         number of rattle constraints to be applied
c     irat         atom numbers of atoms in a rattle constraint
c     ratimage     flag to use minimum image for rattle constraint
c     use_rattle   logical flag to set use of rattle contraints
c
c
      integer nrat,irat
      real*8 krat
      logical ratimage,use_rattle
      common /shake/ krat(maxbnd),nrat,irat(2,maxbnd),ratimage(maxbnd),
     &               use_rattle
