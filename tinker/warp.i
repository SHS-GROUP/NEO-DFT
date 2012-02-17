c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  warp.i  --  parameters for potential surface smoothing  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     m2          second moment of Gaussian representing each atom
c     deform      value of diffusional smoothing deformation parameter
c     diffb       diffusion coefficient for bond stretch potential
c     diffa       diffusion coefficient for angle bend potential
c     diffid      diffusion coefficient for improper dihedral potential
c     difft       diffusion coefficient for torsional potential
c     diffv       diffusion coefficient for van der Waals potential
c     diffc       diffusion coefficient for charge-charge potential
c     use_deform  flag to use diffusion smoothed potential terms
c     use_gda     flag to use Straub's GDA instead of Scheraga's DEM 
c
c
      real*8 m2,deform
      real*8 diffb,diffa,diffid,difft,diffv,diffc
      logical use_deform,use_gda
      common /warp/ m2(maxatm),deform,diffb,diffa,diffid,difft,diffv,
     &              diffc,use_deform,use_gda
