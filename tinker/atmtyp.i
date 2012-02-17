C 21 Apr 10 - NA  - add CHARGE2 support 
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  atmtyp.i  --  atomic properties for each current atom  ##
c     ##                                                         ##
c     #############################################################
c
c
c     mass     atomic weight for each atom in the system
c     tag      integer atom labels from input coordinates file
c     class    atom class number for each atom in the system
c     atomic   atomic number for each atom in the system
c     valence  valence number for each atom in the system
c     name     atom name for each atom in the system
casa  xyzname   atom name for charge2 assignment
c     story    descriptive type for each atom in system
c
c
      integer tag,class,atomic,valence
      real*8 mass
casa   character*10 name
      character*10 name, xyzname
      character*20 story
      common /atmtyp/ mass(maxatm),tag(maxatm),class(maxatm),
     &                atomic(maxatm),valence(maxatm),name(maxatm),
     &                story(maxatm), xyzname(maxatm)
c     &                story(maxatm)
