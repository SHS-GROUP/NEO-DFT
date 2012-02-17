c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  ksttor.i  --  forcefield parameters for stretch-torsions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnbt  maximum number of stretch-torsion parameter entries
c
c     btcon   force constant parameters for stretch-torsion
c     kbt     string of atom types for bonds in stretch-torsion
c
c
      integer maxnbt
      parameter (maxnbt=300)
      real*8 btcon
      character*6 kbt
      common /ksttor/ btcon(3,maxnbt),kbt(maxnbt)
