c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  kdipol.i  --  forcefield parameters for bond dipoles  ##
c     ##                                                        ##
c     ############################################################
c
c
c     maxnd   maximum number of bond dipole parameter entries
c     maxnd5  maximum number of 5-membered ring dipole entries
c     maxnd4  maximum number of 4-membered ring dipole entries
c     maxnd3  maximum number of 3-membered ring dipole entries
c
c     dpl     dipole moment parameters for bond dipoles
c     dpl5    dipole moment parameters for 5-ring dipoles
c     dpl4    dipole moment parameters for 4-ring dipoles
c     dpl3    dipole moment parameters for 3-ring dipoles
c     pos     dipole position parameters for bond dipoles
c     pos5    dipole position parameters for 5-ring dipoles
c     pos4    dipole position parameters for 4-ring dipoles
c     pos3    dipole position parameters for 3-ring dipoles
c     kd      string of atom classes for bond dipole parameters
c     kd5     string of atom classes for 5-ring dipole parameters
c     kd4     string of atom classes for 4-ring dipole parameters
c     kd3     string of atom classes for 3-ring dipole parameters
c
c
      integer maxnd,maxnd5,maxnd4,maxnd3
      parameter (maxnd=300)
      parameter (maxnd5=100)
      parameter (maxnd4=100)
      parameter (maxnd3=100)
      real*8 dpl,dpl5,dpl4,dpl3
      real*8 pos,pos5,pos4,pos3
      character*6 kd,kd5,kd4,kd3
      common /kdipol/ dpl(maxnd),dpl5(maxnd5),dpl4(maxnd4),dpl3(maxnd3),
     &                pos(maxnd),pos5(maxnd5),pos4(maxnd4),pos3(maxnd3),
     &                kd(maxnd),kd5(maxnd5),kd4(maxnd4),kd3(maxnd3)
