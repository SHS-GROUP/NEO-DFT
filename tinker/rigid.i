c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  rigid.i  --  rigid body coordinates for atom groups  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     xrb   rigid body reference x-coordinate for each atom
c     yrb   rigid body reference y-coordinate for each atom
c     zrb   rigid body reference z-coordinate for each atom
c     rbc   current rigid body coordinates for each atom group
c
c
      real*8 xrb,yrb,zrb,rbc
      common /rigid/ xrb(maxatm),yrb(maxatm),zrb(maxatm),rbc(6,maxgrp)
