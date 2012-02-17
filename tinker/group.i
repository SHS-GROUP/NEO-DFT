c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  group.i  --  partitioning of system into atom groups  ##
c     ##                                                        ##
c     ############################################################
c
c
c     wgrp       energetic weight of a group-group interaction
c     grpnum     original group number for each nonempty group
c     ngrp       total number of atom groups in the system
c     kgrp       contiguous list of the atoms in each group
c     igrp       first and last atom of each group in the list
c     grplist    number of the group to which each atom belongs
c     use_group  flag to use partitioning of system into groups
c
c
      integer grpnum,ngrp,kgrp,igrp,grplist
      real*8 wgrp
      logical use_group
      common /group/ wgrp(0:maxgrp,0:maxgrp),grpnum(maxgrp),ngrp,
     &               kgrp(maxatm),igrp(2,0:maxgrp),grplist(maxatm),
     &               use_group
