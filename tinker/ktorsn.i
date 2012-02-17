c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  ktorsn.i  --  forcefield parameters for torsional angles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnt   maximum number of torsion parameter entries
c     maxnt5  maximum number of 5-membered ring torsion entries
c     maxnt4  maximum number of 4-membered ring torsion entries
c
c     t1      torsional parameters for standard 1-fold rotation
c     t2      torsional parameters for standard 2-fold rotation
c     t3      torsional parameters for standard 3-fold rotation
c     t4      torsional parameters for standard 4-fold rotation
c     t5      torsional parameters for standard 5-fold rotation
c     t6      torsional parameters for standard 6-fold rotation
c     t15     torsional parameters for 1-fold rotation in 5-ring
c     t25     torsional parameters for 2-fold rotation in 5-ring
c     t35     torsional parameters for 3-fold rotation in 5-ring
c     t45     torsional parameters for 4-fold rotation in 5-ring
c     t55     torsional parameters for 5-fold rotation in 5-ring
c     t65     torsional parameters for 6-fold rotation in 5-ring
c     t14     torsional parameters for 1-fold rotation in 4-ring
c     t24     torsional parameters for 2-fold rotation in 4-ring
c     t34     torsional parameters for 3-fold rotation in 4-ring
c     t44     torsional parameters for 4-fold rotation in 4-ring
c     t54     torsional parameters for 5-fold rotation in 4-ring
c     t64     torsional parameters for 6-fold rotation in 4-ring
c     kt      string of atom types for standard torsional parameters
c     kt5     string of atom types for 5-ring torsional parameters
c     kt4     string of atom types for 4-ring torsional parameters
c
c
      integer maxnt,maxnt5,maxnt4
      parameter (maxnt=1500)
c Yingbin Ge 11/15/2006: MAXNT5 was raised from 300 to 600
      parameter (maxnt5=600)
      parameter (maxnt4=100)
      real*8 t1,t2,t3,t4,t5,t6
      real*8 t15,t25,t35,t45,t55,t65
      real*8 t14,t24,t34,t44,t54,t64
      character*12 kt,kt5,kt4
      common /ktorsn/ t1(2,maxnt),t2(2,maxnt),t3(2,maxnt),
     &                t4(2,maxnt),t5(2,maxnt),t6(2,maxnt),
     &                t15(2,maxnt5),t25(2,maxnt5),t35(2,maxnt5),
     &                t45(2,maxnt5),t55(2,maxnt5),t65(2,maxnt5),
     &                t14(2,maxnt4),t24(2,maxnt4),t34(2,maxnt4),
     &                t44(2,maxnt4),t54(2,maxnt4),t64(2,maxnt4),
     &                kt(maxnt),kt5(maxnt5),kt4(maxnt4)
