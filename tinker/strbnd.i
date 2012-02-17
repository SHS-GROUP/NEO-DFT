c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  strbnd.i  --  stretch-bends in the current structure  ##
c     ##                                                        ##
c     ############################################################
c
c
c     ksb      force constant for stretch-bend terms
c     nstrbnd  total number of stretch-bend interactions
c     isb      angle and bond numbers used in stretch-bend
c
c
      integer nstrbnd,isb
      real*8 ksb
      common /strbnd/ ksb(maxang),nstrbnd,isb(3,maxang)
