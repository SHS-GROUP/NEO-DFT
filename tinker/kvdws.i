c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  kvdws.i  --  forcefield parameters for van der Waals terms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     rad     van der Waals radius parameter for each atom class
c     eps     van der Waals well depth parameter for each atom class
c     reduct  van der Waals reduction factor for each atom class
c
c
      real*8 rad,eps,reduct
      common /kvdws/ rad(maxclass),eps(maxclass),reduct(maxclass)
