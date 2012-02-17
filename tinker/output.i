c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  output.i  --  parameters to control coordinate output  ##
c     ##                                                         ##
c     #############################################################
c
c
c     use_version  logical flag governing use of filename versions
c     savecycle    flag to mark whether individual cycles are saved
c     coordtype    selects Cartesian, internal, rigid body or none
c
c
      character*9 coordtype
      logical use_version,savecycle
      common /outputt/ use_version,savecycle,coordtype
