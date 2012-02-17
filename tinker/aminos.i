c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  aminos.i  --  standard abbreviations for amino acid types  ##
c     ##                                                             ##
c     #################################################################
c
c
c     amino     three-letter abbreviations for amino acids and caps
c     amino1    one-letter abbreviations for amino acids and caps
c
c
      character*1 amino1(maxamino)
      character*3 amino(maxamino)
      data amino  / 'GLY','ALA','VAL','LEU','ILE','SER','THR','CYS',
     &              'CYX','PRO','PHE','TYR','TRP','HIS','HID','HIE',
     &              'ASP','ASN','GLU','GLN','MET','LYS','ARG','ORN',
     &              'AIB','PCA','UNK','FOR','ACE','NH2','NME' /
      data amino1 / 'G','A','V','L','I','S','T','C','C','P','F','Y',
     &              'W','H','U','Z','D','N','E','Q','M','K','R','O',
     &              'B','J','X','f','a','n','m' /
      save amino,amino1
