c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  mutant.i  --  hybrid atoms for free energy perturbation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     lambda   weighting of initial state in hybrid Hamiltonian
c     nhybrid  number of atoms mutated from initial to final state
c     ihybrid  atomic sites differing in initial and final state
c     type0    atom type of each atom in the initial state system
c     class0   atom class of each atom in the initial state system
c     type1    atom type of each atom in the final state system
c     class1   atom class of each atom in the final state system
c     alter    true if an atom is to be mutated, false otherwise
c
c
      integer nhybrid,ihybrid(maxatm)
      integer type0(maxatm),class0(maxatm)
      integer type1(maxatm),class1(maxatm)
      real*8 lambda
      logical alter(maxatm)
      common /mutant/ lambda,nhybrid,ihybrid,type0,class0,
     &                type1,class1,alter
