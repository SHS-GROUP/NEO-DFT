c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  kbonds.i  --  forcefield parameters for bond stretching  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxnb   maximum number of bond stretch parameter entries
c     maxnb5  maximum number of 5-membered ring bond stretch entries
c     maxnb4  maximum number of 4-membered ring bond stretch entries
c     maxnb3  maximum number of 3-membered ring bond stretch entries
c
c     fcon    force constant parameters for bond stretch
c     fcon5   force constant parameters for 5-ring bond stretch
c     fcon4   force constant parameters for 4-ring bond stretch
c     fcon3   force constant parameters for 3-ring bond stretch
c     blen    bond length parameters for bond stretch
c     blen5   bond length parameters for 5-ring bond stretch
c     blen4   bond length parameters for 4-ring bond stretch
c     blen3   bond length parameters for 3-ring bond stretch
c     kb      string of atom classes for bond stretch parameters
c     kb5     string of atom classes for 5-ring bond stretch parameters
c     kb4     string of atom classes for 4-ring bond stretch parameters
c     kb3     string of atom classes for 3-ring bond stretch parameters
c
c
      integer maxnb,maxnb5,maxnb4,maxnb3
      parameter (maxnb=1000)
      parameter (maxnb5=300)
      parameter (maxnb4=300)
      parameter (maxnb3=300)
      real*8 fcon,fcon5,fcon4,fcon3
      real*8 blen,blen5,blen4,blen3
      character*6 kb,kb5,kb4,kb3
      common /kbonds/ fcon(maxnb),fcon5(maxnb5),fcon4(maxnb4),
     &                fcon3(maxnb3),blen(maxnb),blen5(maxnb5),
     &                blen4(maxnb4),blen3(maxnb3),kb(maxnb),
     &                kb5(maxnb5),kb4(maxnb4),kb3(maxnb3)
