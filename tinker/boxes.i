c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  boxes.i  --  box sizes for periodic boundary conditions  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xbox        length in Angs of a-axis of periodic box
c     ybox        length in Angs of b-axis of periodic box
c     zbox        length in Angs of c-axis of periodic box
c     alpha       angle in degrees between b- and c-axes of box
c     beta        angle in degrees between a- and c-axes of box
c     gamma       angle in degrees between a- and b-axes of box
c     xbox2       half of the a-axis length of periodic box
c     ybox2       half of the b-axis length of periodic box
c     zbox2       half of the c-axis length of periodic box
c     box34       three-fourths axis length of truncated octahedron
c     beta_sin    sine of the beta periodic box angle
c     beta_cos    cosine of the beta periodic box angle
c     gamma_sin   sine of the gamma periodic box angle
c     gamma_cos   cosine of the gamma periodic box angle
c     beta_term   term used in generating triclinic box
c     gamma_term  term used in generating triclinic box
c     orthogonal  flag to mark periodic box as orthogonal
c     monoclinic  flag to mark periodic box as monoclinic
c     triclinic   flag to mark periodic box as triclinic
c     octahedron  flag to mark box as truncated octahedron
c     spacegrp    space group symbol for the unitcell type
c
c
      real*8 xbox,ybox,zbox,alpha,beta,gamma
      real*8 xbox2,ybox2,zbox2,box34
      real*8 beta_sin,beta_cos,gamma_sin
      real*8 gamma_cos,beta_term,gamma_term
      character*10 spacegrp
      logical orthogonal,monoclinic,triclinic,octahedron
      common /boxes/ xbox,ybox,zbox,alpha,beta,gamma,xbox2,ybox2,
     &               zbox2,box34,beta_sin,beta_cos,gamma_sin,
     &               gamma_cos,beta_term,gamma_term,orthogonal,
     &               monoclinic,triclinic,octahedron,spacegrp
