c  May 1, 2001
c  following Brian Salter-Duke, LISTARG is moved before ARG to
c  produce a common aligned on 8 byte as well as 4 byte machines.
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  argue.i  --  command line arguments at program startup  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxarg   maximum number of command line arguments
c
c     narg     number of command line arguments to the program
c     arg      strings containing the command line arguments
c     listarg  flag to mark available command line arguments
c
c
      integer maxarg
      parameter (maxarg=20)
      integer narg
      character*60 arg
      logical listarg
      common /argue/ narg,listarg(0:maxarg),arg(0:maxarg)
