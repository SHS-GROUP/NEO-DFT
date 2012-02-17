c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  angle.i  --  bond angles within the current structure  ##
c     ##                                                         ##
c     #############################################################
c
c
c     anat    ideal bond angle values in degrees
c     acon    bond angle force constants (kcal/mole/rad**2)
c     nangle  total number of bond angles in the system
c     iang    numbers of the atoms in each bond angle
c     angin   logical flag to set use of in-plane projected angle
c
c
      integer nangle,iang
      real*8 anat,acon
      logical angin
      common /angle/ anat(maxang),acon(maxang),nangle,iang(4,maxang),
     &               angin(maxang)
