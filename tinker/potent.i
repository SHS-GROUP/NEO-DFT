c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  potent.i  --  usage of each potential energy component  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     use_bond    logical flag governing use of bond stretch potential
c     use_angle   logical flag governing use of angle bend potential
c     use_strbnd  logical flag governing use of stretch-bend potential
c     use_urey    logical flag governing use of Urey-Bradley potential
c     use_angang  logical flag governing use of angle-angle cross term
c     use_opbend  logical flag governing use of out-of-plane bend term
c     use_improp  logical flag governing use of improper dihedral term
c     use_imptor  logical flag governing use of improper torsion term
c     use_tors    logical flag governing use of torsional potential
c     use_strtor  logical flag governing use of stretch-torsion term
c     use_tortor  logical flag governing use of torsion-torsion term
c     use_vdw     logical flag governing use of vdw der Waals potential
c     use_charge  logical flag governing use of charge-charge potential
c     use_chgdpl  logical flag governing use of charge-dipole potential
c     use_dipole  logical flag governing use of dipole-dipole potential
c     use_mpole   logical flag governing use of multipole potential
c     use_polar   logical flag governing use of polarization term
c     use_rxnfld  logical flag governing use of reaction field term
c     use_solv    logical flag governing use of surface area solvation
c     use_gbsa    logical flag governing use of GB/SA solvation term
c     use_geom    logical flag governing use of geometric restraints
c     use_extra   logical flag governing use of extra potential term
c     use_orbit   logical flag governing use of pi-system computation
c
c
      logical use_bond,use_angle,use_strbnd,use_urey,use_angang
      logical use_opbend,use_improp,use_imptor,use_tors,use_strtor
      logical use_tortor,use_vdw,use_charge,use_chgdpl,use_dipole
      logical use_mpole,use_polar,use_rxnfld,use_solv,use_gbsa
      logical use_geom,use_extra,use_orbit
      common /potent/ use_bond,use_angle,use_strbnd,use_urey,use_angang,
     &                use_opbend,use_improp,use_imptor,use_tors,
     &                use_strtor,use_tortor,use_vdw,use_charge,
     &                use_chgdpl,use_dipole,use_mpole,use_polar,
     &                use_rxnfld,use_solv,use_gbsa,use_geom,use_extra,
     &                use_orbit
