c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  cutoff.i  --  cutoff distances for energy interactions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     vdwcut       cutoff distance for van der Waals interactions
c     vdwtaper     distance at which van der Waals switching begins
c     chgcut       cutoff distance for charge-charge interactions
c     chgtaper     distance at which charge-charge switching begins
c     dplcut       cutoff distance for dipole-dipole interactions
c     dpltaper     distance at which dipole-dipole switching begins
c     use_lights   flag to use method of lights neighbor generation
c
c
      real*8 vdwcut,vdwtaper,chgcut,chgtaper,dplcut,dpltaper
      logical use_lights
      common /cutoff/ vdwcut,vdwtaper,chgcut,chgtaper,dplcut,dpltaper,
     &                use_lights
