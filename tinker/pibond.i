c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  pibond.i  --  bond orders for a conjugated pi-system  ##
c     ##                                                        ##
c     ############################################################
c
c
c     pbpl    pi-bond orders for bonds in "planar" pi-system
c     pnpl    pi-bond orders for bonds in "nonplanar" pi-system
c
c
      real*8 pbpl,pnpl
      common /pibond/ pbpl(maxpib),pnpl(maxpib)
