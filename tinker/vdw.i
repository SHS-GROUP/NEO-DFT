c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  vdw.i  --  van der Waals parameters for current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     radmin   minimum energy distance for each atom class pair
c     epsilon  well depth parameter for each atom class pair
c     radhbnd  minimum energy distance for hydrogen bonding pairs
c     epshbnd  well depth parameter for hydrogen bonding pairs
c     kred     value of reduction factor parameter for each atom
c     ired     attached atom from which reduction factor is applied
c     nvdw     total number van der Waals active sites in the system
c     ivdw     number of the atom for each van der Waals active site
c
c
      integer ired,nvdw,ivdw
      real*8 radmin,epsilon,radhbnd,epshbnd,kred
      common /vdw/ radmin(maxclass,maxclass),
     &             epsilon(maxclass,maxclass),
     &             radhbnd(maxclass,maxclass),
     &             epshbnd(maxclass,maxclass),
     &             kred(maxatm),ired(maxatm),nvdw,ivdw(maxatm)
