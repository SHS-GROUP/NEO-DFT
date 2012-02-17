c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  angpot.i  --  specifics of bond angle functional form  ##
c     ##                                                         ##
c     #############################################################
c
c
c     angunit   convert angle force constant to kcal/mole/deg**2
c     cang      cubic term in angle bending potential
c     qang      quartic term in angle bending potential
c     pang      quintic term in angle bending potential
c     sang      sextic term in angle bending potential
c     aaunit    convert angle-angle constant to kcal/mole/deg**2
c     opbunit   convert out-plane bend constant to kcal/mole/deg**2
c     stbnunit  convert str-bnd constant to kcal/mole/deg-Ang**2
c     angtyp    type of angle bending potential energy function
c
c
      real*8 angunit,cang,qang,pang,sang
      real*8 aaunit,opbunit,stbnunit
      character*8 angtyp
      common /angpot/ angunit,cang,qang,pang,sang,aaunit,opbunit,
     &                stbnunit,angtyp
