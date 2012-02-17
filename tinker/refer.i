c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  refer.i  --  storage of reference atomic coordinate set  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xref       reference x-coordinate for each atom in the system
c     yref       reference y-coordinate for each atom in the system
c     zref       reference z-coordinate for each atom in the system
c     nref       total number of atoms in the reference system
c     reftyp     atom type for each atom in the reference system
c     n12ref     number of atoms bonded to each reference atom
c     i12ref     atom numbers of atoms 1-2 connected to each atom
c     refleng    length in characters of the reference filename
c     refltitle  length in characters of the reference title string
c     refnam     atom name for each atom in the reference system
c     reffile    base filename for the reference structure
c     reftitle   title used to describe the reference structure
c
c
      integer nref,reftyp,n12ref,i12ref,refleng,refltitle
      real*8 xref,yref,zref
      character*3 refnam
      character*60 reffile,reftitle
      common /refert/ 
     #               xref(maxatm),yref(maxatm),zref(maxatm),nref,
     &               reftyp(maxatm),n12ref(maxatm),
     &               i12ref(maxval,maxatm),refleng,refltitle,
     &               refnam(maxatm),reffile,reftitle
