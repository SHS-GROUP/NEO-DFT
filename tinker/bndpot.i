c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  bndpot.i  --  specifics of bond stretch functional form  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     bndunit  convert bond force constant to kcal/mole/Ang**2
c     cbnd     cubic term in bond stretch potential
c     qbnd     quartic term in bond stretch potential
c     bndtyp   type of bond stretch potential energy function
c
c
      real*8 bndunit,cbnd,qbnd
      character*8 bndtyp
      common /bndpot/ bndunit,cbnd,qbnd,bndtyp
