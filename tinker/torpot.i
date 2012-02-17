c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  torpot.i  --  specifics of torsional functional forms  ##
c     ##                                                         ##
c     #############################################################
c
c
c     torsunit  scale factor for torsional parameter amplitudes
c     storunit  convert stretch-torsion force to kcal/mole/Ang
c
c
      real*8 torsunit,storunit
      common /torpot/ torsunit,storunit
