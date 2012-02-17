c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kopbnd.i  --  forcefield parameters for out-of-plane bend  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxnopb  maximum number of out-of-plane bending entries
c
c     copb     force constant parameters for out-of-plane bending
c     kaopb    string of atom types for out-of-plane bend parameters
c
c
      integer maxnopb
      parameter (maxnopb=200)
      real*8 copb
      character*6 kaopb
      common /kopbnd/ copb(maxnopb),kaopb(maxnopb)
