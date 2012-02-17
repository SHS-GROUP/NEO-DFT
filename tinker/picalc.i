c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  picalc.i  --  orbital energies for conjugated pi-system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     q       number of pi-electrons contributed by each atom
c     w       ionization potential of each pi-system atom
c     em      repulsion integral for each pi-system atom
c     nfill   number of filled pi-system molecular orbitals
c
c
      integer nfill
      real*8 q,w,em
      common /picalc/ q(maxpi),w(maxpi),em(maxpi),nfill
