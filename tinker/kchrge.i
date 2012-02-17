C 21 Apr 10 - NA  - add CHARGE2 support
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  kchrge.i  --  forcefield parameters for partial charges  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     chg       partial charge parameters for each atom type
casa  -charge2- begin
c     chg2      partial charge parameters for ligand
c     n_chg2  number of charge2 parameter 
c     atname    atom name for charge2 parameter
c
c
c      real*8 chg
c      common /kchrge/ chg(maxtyp)
      real*8 chg, chg2
      integer n_chg2
      character*10 atname
      common /kchrge/ chg(maxtyp), chg2(maxatm), atname(maxatm), n_chg2
