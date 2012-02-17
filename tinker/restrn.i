C 21 Apr 10 - NA  - improvement for RESTRAIN-POSITION optimization 
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  restrn.i  --  definition of the geometrical restraints  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     depth      depth of shallow Gaussian basin restraint
c     width      exponential width coefficient of Gaussian basin
c     rwall      radius of spherical droplet boundary restraint
c     xpfix      x-coordinate target for each restrained position
c     ypfix      y-coordinate target for each restrained position
c     zpfix      z-coordinate target for each restrained position
c     pfix       force constant for each restrained position
c     dfix       target range and force constant for each distance
c     tfix       target range and force constant for each torsion
c     npfix      number of position restraints to be applied
c     ipfix      atom number involved in each position restraint
c     ndfix      number of distance restraints to be applied
c     idfix      atom numbers defining each distance restraint
c     ntfix      number of torsional restraints to be applied
c     itfix      atom numbers defining each torsional restraint
c     use_basin  logical flag governing use of Gaussian basin
c     use_wall   logical flag governing use of droplet boundary
c
casa  -restrain- begin
c     npfixm     number of position restraints for mm calculation
c     nfix       number of position restraints for fmomm calculation
c
c      integer npfix,ipfix,ndfix
      integer npfix,ipfix,ndfix, npfixm, nfix
      real*8 RESTRAINV
      PARAMETER (RESTRAINV=1.0D+8)
casa  -restrain- end
      integer idfix,ntfix,itfix
      real*8 depth,width,rwall
      real*8 xpfix,ypfix,zpfix
      real*8 pfix,dfix,tfix
      logical use_basin,use_wall
      common /restrn/ depth,width,rwall,xpfix(maxfix),ypfix(maxfix),
     &                zpfix(maxfix),pfix(maxfix),dfix(3,maxfix),
     &                tfix(3,maxfix),npfix,ipfix(maxfix),ndfix,
     &                idfix(2,maxfix),ntfix,itfix(4,maxfix),use_basin,
     &                use_wall, npfixm, nfix
c     &                use_wall
