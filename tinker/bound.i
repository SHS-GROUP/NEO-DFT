c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  bound.i  --  control of periodic boundary conditions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     use_bounds   flag to use periodic boundary conditions
c     use_image    flag to use images for periodic system
c     use_replica  flag to use replicates for periodic system
c
c
      logical use_bounds,use_image,use_replica
      common /bound/ use_bounds,use_image,use_replica
