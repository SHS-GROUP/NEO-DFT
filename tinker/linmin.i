c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  linmin.i  --  parameters for line search minimization   ##
c     ##                                                          ##
c     ##############################################################
c
c
c     cappa   stringency of line search (0=tight < cappa < 1=loose)
c     stpmin  minimum step length in current line search direction
c     stpmax  maximum step length in current line search direction
c     angmax  maximum angle between search direction and -gradient
c     intmax  maximum number of cubic interpolation attempts
c
c
      integer intmax
      real*8 cappa,stpmin,stpmax,angmax
      common /linmin/ cappa,stpmin,stpmax,angmax,intmax
