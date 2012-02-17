c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  disgeo.i  --  distance geometry bounds and parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     bnd         distance geometry upper and lower bounds matrix
c     vdwrad      hard sphere radii for distance geometry atoms
c     vchir       signed volume values for chirality constraints
c     compact     index of local distance compaction on embedding
c     pathmax     maximum value of upper bound after smoothing
c     nchir       total number of chirality constraints
c     ichir       numbers of atoms in each chirality constraint
c     use_invert  flag to use enantiomer closest to input structure
c     use_anneal  flag to use simulated annealing refinement
c
c
      integer nchir,ichir
      real*8 bnd,vdwrad,vchir
      real*8 compact,pathmax
      logical use_invert,use_anneal
      common /disgeo/ bnd(maxgeo,maxgeo),vdwrad(maxatm),vchir(maxatm),
     &                compact,pathmax,nchir,ichir(4,maxatm),use_invert,
     &                use_anneal
