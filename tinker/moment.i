c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  moment.i  --  components of the net multipole moments  ##
c     ##                                                         ##
c     #############################################################
c
c
c     netchg   net electric charge on the total system
c     xdipole  dipole moment of the system along the x-axis
c     ydipole  dipole moment of the system along the y-axis
c     zdipole  dipole moment of the system along the z-axis
c
c
      real*8 netchg,xdipole,ydipole,zdipole
      common /moment/ netchg,xdipole,ydipole,zdipole
