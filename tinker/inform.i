c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  inform.i  --  control flags for I/O and program flow  ##
c     ##                                                        ##
c     ############################################################
c
c
c     verbose  logical flag to turn on extra information
c     debug    logical flag to turn on full debug printing
c     abort    logical flag to stop execution at next chance
c
c
      logical verbose,debug,abort
      common /inform/ verbose,debug,abort
