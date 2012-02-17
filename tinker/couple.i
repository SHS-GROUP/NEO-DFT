c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  couple.i  --  near-neighbor atom connectivity lists   ##
c     ##                                                        ##
c     ############################################################
c
c
c     n12     number of atoms directly bonded to each atom
c     i12     atom numbers of atoms 1-2 connected to each atom
c     n13     number of atoms in a 1-3 relation to each atom
c     i13     atom numbers of atoms 1-3 connected to each atom
c     n14     number of atoms in a 1-4 relation to each atom
c     i14     atom numbers of atoms 1-4 connected to each atom
c
c
      integer n12,i12,n13,i13,n14,i14
      common /couple/ n12(maxatm),i12(maxval,maxatm),n13(maxatm),
     &                i13(maxval*(maxval-1),maxatm),n14(maxatm),
     &                i14(maxval*(maxval-1)**2,maxatm)
