c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  solute.i  --  surface area-based macroscopic solvation  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     rsolv    atomic radius of each atom for empirical solvation
c     vsolv    atomic solvation parameters (kcal/mole/Ang**2)
c     rborn    Born radius of each atom for GB/SA solvation
c     nsolv    number of atoms with non-zero solvation parameters
c     reborn   number of GB/SA calculations between Born radii updates
c     bornmax  maximum atoms for original Born radii computation
c
c
      integer nsolv,reborn,bornmax
      real*8 rsolv,vsolv,rborn
      common /solute/ rsolv(maxatm),vsolv(maxatm),rborn(maxatm),
     &                nsolv,reborn,bornmax
