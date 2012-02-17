c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kiprop.i  --  forcefield parameters for improper dihedral  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxndi  maximum number of improper dihedral parameter entries
c
c     dcon    force constant parameters for improper dihedrals
c     tdi     ideal dihedral angle values for improper dihedrals
c     kdi     atom types for improper dihedral angle parameters
c
c
      integer maxndi
      parameter (maxndi=200)
      real*8 dcon,tdi
      character*12 kdi
      common /kiprop/ dcon(maxndi),tdi(maxndi),kdi(maxndi)
