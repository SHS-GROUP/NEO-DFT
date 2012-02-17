c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  kangs.i  --  forcefield parameters for bond angle bending  ##
c     ##                                                             ##
c     #################################################################
c
c
c     maxna    maximum number of bond angle bend parameter entries
c     maxna5   maximum number of 5-membered ring angle bend entries
c     maxna4   maximum number of 4-membered ring angle bend entries
c     maxna3   maximum number of 3-membered ring angle bend entries
c
c     con      force constant parameters for angle bending
c     con5     force constant parameters for 5-ring angle bends
c     con4     force constant parameters for 4-ring angle bends
c     con3     force constant parameters for 3-ring angle bends
c     ang      bond angle parameters for angle bending
c     ang5     bond angle parameters for 5-ring angle bends
c     ang4     bond angle parameters for 4-ring angle bends
c     ang3     bond angle parameters for 3-ring angle bends
c     ka       string of atom classes for angle bending parameters
c     ka5      string of atom classes for 5-ring angle bend parameters
c     ka4      string of atom classes for 4-ring angle bend parameters
c     ka3      string of atom classes for 3-ring angle bend parameters
c
c
      integer maxna,maxna5,maxna4,maxna3
      parameter (maxna=10000)
      parameter (maxna5=2500)
      parameter (maxna4=1200)
      parameter (maxna3=1200)
      real*8 con,con5,con4,con3
      real*8 ang,ang5,ang4,ang3
      character*9 ka,ka5,ka4,ka3
      common /kangs/ con(maxna),con5(maxna5),con4(maxna4),con3(maxna3),
     &               ang(3,maxna),ang5(3,maxna5),ang4(3,maxna4),
     &               ang3(3,maxna3),ka(maxna),ka5(maxna5),ka4(maxna4),
     &               ka3(maxna3)
