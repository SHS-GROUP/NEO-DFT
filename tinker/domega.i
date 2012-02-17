c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  domega.i  --  derivative components over dihedrals  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     teb     bond stretch derivatives over torsions
c     tea     angle bend derivatives over torsions
c     teba    stretch-bend derivatives over torsions
c     teub    Urey-Bradley derivatives over torsions
c     teaa    angle-angle derivatives over torsions
c     teopb   out-of-plane bend derivatives over torsions
c     teid    improper dihedral derivatives over torsions
c     teit    improper torsion derivatives over torsions
c     tet     torsional derivatives over torsions
c     tebt    stretch-torsion derivatives over torsions
c     tett    torsion-torsion derivatives over torsions
c     tev     van der Waals derivatives over torsions
c     te14    1-4 van der Waals derivatives over torsions
c     tec     charge-charge derivatives over torsions
c     tecd    charge-dipole derivatives over torsions
c     ted     dipole-dipole derivatives over torsions
c     tem     atomic multipole derivatives over torsions
c     tep     polarization derivatives over torsions
c     ter     reaction field derivatives over torsions
c     tes     solvation derivatives over torsions
c     teg     geometric restraint derivatives over torsions
c     tex     extra energy term derivatives over torsions
c
c
      real*8 teb,tea,teba,teub,teaa,teopb,teid,teit,tet,tebt
      real*8 tett,tev,te14,tec,tecd,ted,tem,tep,ter,tes,teg,tex
      common /domega/ teb(maxrot),tea(maxrot),teba(maxrot),
     &                teub(maxrot),teaa(maxrot),teopb(maxrot),
     &                teid(maxrot),teit(maxrot),tet(maxrot),
     &                tebt(maxrot),tett(maxrot),tev(maxrot),
     &                te14(maxrot),tec(maxrot),tecd(maxrot),
     &                ted(maxrot),tem(maxrot),tep(maxrot),
     &                ter(maxrot),tes(maxrot),teg(maxrot),
     &                tex(maxrot)
