c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  pistuf.i  --  conjugated system in the current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     norbit   total number of pisystem orbitals in the system
c     iorbit   numbers of the atoms containing pisystem orbitals
c     piperp   atoms defining a normal plane to each orbital
c     npibond  total number of bonds affected by the pisystem
c     pibond   bond and piatom numbers for each pisystem bond
c     npitors  total number of torsions affected by the pisystem
c     pitors   torsion and pibond numbers for each pisystem torsion
c     listpi   atom list indicating whether each atom has an orbital
c
c
      integer norbit,iorbit,piperp
      integer npibond,pibond,npitors,pitors
      logical listpi
      common /pistuf/ norbit,iorbit(maxpi),piperp(3,maxpi),npibond,
     &                pibond(3,maxpib),npitors,pitors(2,maxpit),
     &                listpi(maxatm)
