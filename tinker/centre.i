c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  centre.i  --  coordinates relative to molecule centroid  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xcm     x-offset of atom from molecular center of mass
c     ycm     y-offset of atom from molecular center of mass
c     zcm     z-offset of atom from molecular center of mass
c
c
      real*8 xcm,ycm,zcm
      common /centre/ xcm(maxatm),ycm(maxatm),zcm(maxatm)
