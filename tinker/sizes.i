c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  sizes.i  --  parameter values to set array dimensions  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sizes.i" sets values for critical array dimensions used
c     throughout the software; these parameters will fix the size
c     of the largest systems that can be handled; values too large
c     for the computer's memory and/or swap space to accomodate
c     will result in poor performance or outright failure
c
c     parameter:      maximum allowed number of:
c
c     maxatm          atoms in the molecular system
c     maxval          atoms directly bonded to an atom
c     maxgrp          user-defined groups of atoms
c     maxtyp          force field atom type definitions
c     maxclass        force field atom class definitions
c     maxkey          lines in the keyword file
c     maxrot          bonds for torsional rotation
c     maxopt          optimization variables (full-matrix)
c     maxvib          vibrational frequencies
c     maxgeo          distance geometry points (unused now, see Libtemo.f) 
c     maxpi           atoms in conjugated pisystem
c     maxcell         unit cells in replicated crystal
c     maxring         3-, 4-, or 5-membered rings
c     maxfix          geometric restraints
c     maxbio          biopolymer atom definitions
c     maxres          residues in the macromolecule
c     maxamino        amino acid residue types
c     maxvar          optimization variables (linear storage)
c     maxbnd          covalent bonds in molecular system
c     maxang          bond angles in molecular system
c     maxtors         dihedral angles in molecular system
c     maxlight        sites for method of lights neighbors
c     maxpib          covalent bonds in pisystem
c     maxpit          dihedrals involving pisystem
c
c
      integer maxatm,maxval,maxgrp,maxtyp,maxclass,maxkey,maxrot
      integer maxopt,maxvib,maxpi,maxcell,maxring
      integer maxfix,maxbio,maxres,maxamino,maxvar,maxbnd,maxang
      integer maxtors,maxlight,maxpib,maxpit
c
c        Jim Shoemaker chose to use a thousand atoms, back in the day.
c                    Achtung!    Achtung!     Achtung!
c        MAXATM in inputb.src must be changed to match this parameter,
c        in the event you decide to change the value just below.
c
      parameter (maxatm=12000)
c
c        coordination number was changed by Yingbin Ge in March 2007,
c        from 4 to 6.
      parameter (maxval=6)
c
      parameter (maxgrp=100)
      parameter (maxtyp=3000)
      parameter (maxclass=500)
      parameter (maxkey=5000)
      parameter (maxrot=500)
c
      parameter (maxopt=500)
      parameter (maxvib=500)
c     parameter (maxgeo=500)
      parameter (maxpi=50)
      parameter (maxcell=500)
c         changed from 500 to 1200 in August 2007 for Luke
      parameter (maxring=1200)
      parameter (maxfix=12000)
      parameter (maxbio=5000)
      parameter (maxres=5000)
      parameter (maxamino=31)
      parameter (maxvar=3*maxatm)
c
c         Luke Roskop's bond counting, from thinking about the diamond
c         lattice, was introduced November 2006.  We are not sure how
c         they count torsions, so we guessed at that.
      parameter (maxbnd=2*maxatm)
      parameter (maxang=4*maxatm)
      parameter (maxtors=8*maxatm)
c         original settings are probably for polypeptides
c---  parameter (maxbnd=6*maxatm/5)
c---  parameter (maxang=12*maxatm/5)
c---  parameter (maxtors=4*maxatm)
c
      parameter (maxlight=8*maxatm)
      parameter (maxpib=6*maxpi/5)
      parameter (maxpit=4*maxpi)
