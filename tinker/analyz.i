c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  analyz.i  --  energy components partitioned over atoms  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     aeb     bond stretch energy partitioned over atoms
c     aea     angle bend energy partitioned over atoms
c     aeba    stretch-bend energy partitioned over atoms
c     aeub    Urey-Bradley energy partitioned over atoms
c     aeaa    angle-angle energy partitioned over atoms
c     aeopb   out-of-plane bend energy partitioned over atoms
c     aeid    improper dihedral energy partitioned over atoms
c     aeit    improper torsion energy partitioned over atoms
c     aet     torsional energy partitioned over atoms
c     aebt    stretch-torsion energy partitioned over atoms
c     aett    torsion-torsion energy partitioned over atoms
c     aev     van der Waals energy partitioned over atoms
c     ae14    1-4 van der Waals energy partitioned over atoms
c     aec     charge-charge energy partitioned over atoms
c     aecd    charge-dipole energy partitioned over atoms
c     aed     dipole-dipole energy partitioned over atoms
c     aem     multipole energy partitioned over atoms
c     aep     polarization energy partitioned over atoms
c     aer     reaction field energy partitioned over atoms
c     aes     solvation energy partitioned over atoms
c     aeg     geometric restraint energy partitioned over atoms
c     aex     extra energy term partitioned over atoms
c
c
      real*8 aeb,aea,aeba,aeub,aeaa,aeopb,aeid,aeit,aet,aebt
      real*8 aett,aev,ae14,aec,aecd,aed,aem,aep,aer,aes,aeg,aex
      common /analyz/ aeb(maxatm),aea(maxatm),aeba(maxatm),
     &                aeub(maxatm),aeaa(maxatm),aeopb(maxatm),
     &                aeid(maxatm),aeit(maxatm),aet(maxatm),
     &                aebt(maxatm),aett(maxatm),aev(maxatm),
     &                ae14(maxatm),aec(maxatm),aecd(maxatm),
     &                aed(maxatm),aem(maxatm),aep(maxatm),
     &                aer(maxatm),aes(maxatm),aeg(maxatm),
     &                aex(maxatm)
