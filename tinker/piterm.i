c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  piterm.i  --  bonds and torsions in the current pi-system  ##
c     ##                                                             ##
c     #################################################################
c
c
c     bkpi    bond stretch force constants for pi-bond order of 1.0
c     blpi    ideal bond length values for a pi-bond order of 1.0
c     kslope  rate of force constant decrease with bond order decrease
c     lslope  rate of bond length increase with a bond order decrease
c     torspi  2-fold torsional energy barrier for pi-bond order of 1.0
c
c
      real*8 bkpi,blpi,kslope,lslope,torspi
      common /piterm/ bkpi(maxpib),blpi(maxpib),kslope(maxpib),
     &                lslope(maxpib),torspi(maxpit)
