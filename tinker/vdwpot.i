c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  vdwpot.i  --  specifics of van der Waals functional form  ##
c     ##                                                            ##
c     ################################################################
c
c
c     aterm     value of the "A" constant in exp-6 vdw potential
c     bterm     value of the "B" constant in exp-6 vdw potential
c     cterm     value of the "C" constant in exp-6 vdw potential
c     vdwscale  factor by which 1-4 vdw interactions are scaled
c     igauss    coefficients of Gaussian fit to vdw potential
c     ngauss    number of Gaussians used in fit to vdw potential
c     vdw12use  usage of 1-2 vdw terms (1=omit, -1=use scaled value)
c     vdw13use  usage of 1-3 vdw terms (1=omit, -1=use scaled value)
c     vdw14use  usage of 1-4 vdw terms (1=omit, -1=use scaled value)
c     vdwtyp    type of van der Waals potential energy function
c     radtyp    type of parameter (sigma or R-min) for atomic size
c     radsiz    atomic size provided as radius or diameter
c     radrule   combining rule for atomic size parameters
c     epsrule   combining rule for vdw well depth parameters
c     gausstyp  type of Gaussian fit to van der Waals potential
c
c
      real*8 aterm,bterm,cterm,vdwscale,igauss
      integer ngauss,vdw12use,vdw13use,vdw14use
      character*5 radtyp
      character*8 radsiz,gausstyp
      character*13 vdwtyp
      character*10 radrule,epsrule
      common /vdwpot/ aterm,bterm,cterm,vdwscale,igauss(2,4),ngauss,
     &                vdw12use,vdw13use,vdw14use,vdwtyp,radtyp,radsiz,
     &                radrule,epsrule,gausstyp
