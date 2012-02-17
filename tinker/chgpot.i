c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  chgpot.i  --  specifics of electrostatics functional form  ##
c     ##                                                             ##
c     #################################################################
c
c
c     dielec    dielectric constant for electrostatic interactions
c     chgscale  factor by which 1-4 electrostatic terms are scaled
c     chg12use  usage of 1-2 electrostatics (0=use, 1=omit, -1=scale)
c     chg13use  usage of 1-3 electrostatics (0=use, 1=omit, -1=scale)
c     chg14use  usage of 1-4 electrostatics (0=use, 1=omit, -1=scale)
c     neutnbr   logical flag governing use of neutral group neighbors
c     neutcut   logical flag governing use of neutral group cutoffs
c
c
      real*8 dielec,chgscale
      integer chg12use,chg13use,chg14use
      logical neutcut,neutnbr
      common /electr/ dielec,chgscale,chg12use,chg13use,chg14use,
     &                neutnbr,neutcut
