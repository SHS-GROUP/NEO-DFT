c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  units.i  --  physical constants and unit conversions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     avogadro   Avogadro's number (N) in particles/mole
c     boltzmann  Boltzmann constant (kB) in g*Ang**2/s**2/K/mole
c     gasconst   ideal gas constant (R) in kcal/mole/K
c     lightspd   speed of light in vacuum (c) in cm/sec
c     bohr       conversion from Bohrs to Angstroms
c     joule      conversion from calories to joules
c     evolt      conversion from Hartree to electron-volts
c     hartree    conversion from Hartree to kcal/mole
c     electric   conversion from electron**2/Ang to kcal/mole
c     debye      conversion from electron-Ang to Debyes
c     prescon    conversion from kcal/mole/Ang**3 to Atm
c     convert    conversion from kcal to g*Ang**2/s**2
c
c
      real*8 avogadro,boltzmann,gasconst,lightspd
      real*8 bohr,joule,evolt,hartree
      real*8 electric,debye,prescon,convert
      parameter (avogadro=6.022045d+23)
      parameter (boltzmann=8.3143435d+23)
      parameter (gasconst=1.98717623d-3)
      parameter (lightspd=2.99792458d+10)
      parameter (bohr=0.529177249d0)
      parameter (joule=4.184d0)
      parameter (evolt=27.2107d0)
      parameter (hartree=627.503d0)
      parameter (electric=332.05382d0)
      parameter (debye=4.8033324d0)
      parameter (prescon=6.85695d+4)
      parameter (convert=4.184d+26)
