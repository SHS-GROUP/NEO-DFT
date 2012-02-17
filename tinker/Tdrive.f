c 08 Jun 10 - DGF - dynamic memory allocation for the Hessian
C 21 Apr 10 - NA,KK - allow choosing minimizer, RESTRAIN-POSITION changes
C 24 Sep 01 - JMR - tnewtx: change format statements
C  9 MAR 00 - CHC - Modified output option
c 20 OCT 98 - CHC - Modified for parallel run
c 12 MAY 98 - JRS - NEWTON:  Converted to subroutine NEWTZ
c                            Initialization now done in TINKIN
c                            Command line entry of options removed
c  7 MAY 98 - JRS - ANALYZE: Converted to subroutine ANLYZX 
c                            Initialization now done in TINKIN
c                            Command line entry of options removed
c         CHANGED INTO SUBROUTINES TO BE CALLED FROM GAMESS
C
C     
c
c
c
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program analyze  --  energy partitioning and analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "analyze" computes and displays the total potential; options
c     are provided to partition the energy by atom or by potential
c     function type; parameters used in computing interactions can
c     also be displayed by atom; output of large energy interactions
c     and the total dipole moment and components are available
c
c
cjrs
c      program analyze
cjrs
      subroutine analyze
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'angang.i'
      include 'angle.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'dipole.i'
      include 'energi.i'
      include 'improp.i'
      include 'imptor.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kvdws.i'
      include 'moment.i'
      include 'mpole.i'
      include 'opbend.i'
      include 'polar.i'
      include 'potent.i'
      include 'restrn.i'
      include 'solute.i'
      include 'strbnd.i'
      include 'strtor.i'
      include 'tors.i'
      include 'units.i'
      include 'urey.i'
      include 'vdw.i'
ckk -fmoimomm- begin
      include 'domega.i'
ckk -fmoimomm- end
cjrs added
      include 'deriv.i'
      include 'usage.i'
      real*8 derivs(3,maxatm)
      real*8 gnorm,grms
      integer nvar
cjrs
      integer i,j,k,trimtext
      integer ia,ib,ic,id
      integer list(20),fold(6)
      real*8 energy,rg,moment
      real*8 ampli(6),phase(6)
      character*1 letter
      character*80 record,string
      logical doenergy,doatom,dolarge
      logical doprops,doelect,doparam
      logical exist,header,active(maxatm)
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
ckk -fmoimomm- begin
      real*8 engmm,egmm,tegtmp
      integer MXATM,nd
      PARAMETER (MXATM=2000)
      COMMON /FUNCT / ENGMM,EGMM(3*MXATM)
      COMMON /TGRAD/  TEGTMP(3*MAXATM)
ckk -fmoimomm- end
c
c     set up the structure and mechanics calculation
c
c     set option control flags based desired analysis types
c
      doenergy = .true.
cjrs
c      if (maswrk) write(*,*) 'n on entry to analyze ',n
c
c
      doatom = .false.
      dolarge = .false.
      doprops = .false.
      doelect = .false.
      doparam = .false.
c      call upcase (string)
c      do i = 1, trimtext(string)
c         letter = string(i:i)
c         if (letter .eq. 'E') then
c            doenergy = .true.
c         else if (letter .eq. 'A') then
c            doatom = .true.
c         else if (letter .eq. 'L') then
c            dolarge = .true.
c         else if (letter .eq. 'I') then
c            doprops = .true.
c         else if (letter .eq. 'M') then
c            doelect = .true.
c         else if (letter .eq. 'P') then
c            doparam = .true.
c         end if
c      end do
c
c     get the list of atoms for which output is desired
c
      if (doatom .or. doparam) then
         do i = 1, 20
            list(i) = 0
         end do
         do i = 1, n
            active(i) = .true.
         end do
         i = 1
         dowhile (list(i) .ne. 0)
            if (i .eq. 1) then
               do j = 1, n
                  active(j) = .false.
               end do
            end if
            if (list(i) .gt. 0) then
               active(list(i)) = .true.
               i = i + 1
            else
               do j = abs(list(i)), abs(list(i+1))
                  active(j) = .true.
               end do
               i = i + 2
            end if
         end do
      end if
c
c     if desired, write out some of the individual energy terms
c
      if (dolarge) then
         verbose = .true.
      else
cjrs 
c        verbose = .false.
cjrs
      end if
c      if (maswrk) write(*,*) ' debug  verbose ',debug, verbose
c
c     now make the call to compute the potential energy
c
      if (doenergy .or. doatom .or. dolarge) then
         call analysis (energy)
         if (abs(energy) .lt. 1.0d10) then
            if (maswrk) write (iout,90)  energy
   90       format (/,' Total Potential Energy :',6x,f16.4,' Kcal/mole')
         else
            if (maswrk) write (iout,100)  energy
  100       format (/,' Total Potential Energy :',6x,d16.4,' Kcal/mole')
         end if
      else if (doprops) then
         call analysis (energy)
      end if
c
c     energy partitioning by types of potential energy
c
      if (doenergy) then
         if (maswrk) write (iout,110)
  110    format (/,' Energy Component Breakdown :',9x,'Kcal/mole',
     &             6x,'Interactions'/)
         if (use_bond .and. neb.ne.0) then
            if (maswrk) write (iout,120)  eb,neb
  120       format (' Bond Stretching',15x,f16.4,i14)
         end if
         if (use_angle .and. nea.ne.0) then
            if (maswrk) write (iout,130)  ea,nea
  130       format (' Angle Bending',17x,f16.4,i14)
         end if
         if (use_strbnd .and. neba.ne.0) then
            if (maswrk) write (iout,140)  eba,neba
  140       format (' Stretch-Bend',18x,f16.4,i14)
         end if
         if (use_urey .and. neub.ne.0) then
            if (maswrk) write (iout,150)  eub,neub
  150       format (' Urey-Bradley',18x,f16.4,i14)
         end if
         if (use_angang .and. neaa.ne.0) then
            if (maswrk) write (iout,160)  eaa,neaa
  160       format (' Angle-Angle',19x,f16.4,i14)
         end if
         if (use_opbend .and. neopb.ne.0) then
            if (maswrk) write (iout,170)  eopb,neopb
  170       format (' Out-of-Plane Bend',13x,f16.4,i14)
         end if
         if (use_improp .and. neid.ne.0) then
            if (maswrk) write (iout,180)  eid,neid
  180       format (' Improper Dihedral',13x,f16.4,i14)
         end if
         if (use_imptor .and. neit.ne.0) then
            if (maswrk) write (iout,190)  eit,neit
  190       format (' Improper Torsion',14x,f16.4,i14)
         end if
         if (use_tors .and. net.ne.0) then
            if (maswrk) write (iout,200)  et,net
  200       format (' Torsional Angle',15x,f16.4,i14)
         end if
         if (use_strtor .and. nebt.ne.0) then
            if (maswrk) write (iout,210)  ebt,nebt
  210       format (' Stretch-Torsion',15x,f16.4,i14)
         end if
         if (use_tortor .and. nett.ne.0) then
            if (maswrk) write (iout,220)  ett
  220       format (' Torsion-Torsion',15x,f16.4,i14)
         end if
         if (use_vdw .and. ne14.ne.0) then
            if (maswrk) write (iout,230)  e14,ne14
  230       format (' 1-4 van der Waals',13x,f16.4,i14)
         end if
         if (use_vdw .and. nev.ne.0) then
            if (abs(ev) .lt. 1.0d10) then
               if (maswrk) write (iout,240)  ev,nev
  240          format (' Other van der Waals',11x,f16.4,i14)
            else
               if (maswrk) write (iout,250)  ev,nev
  250          format (' Other van der Waals',11x,d16.4,i14)
            end if
         end if
         if (use_charge .and. nec.ne.0) then
            if (abs(ec) .lt. 1.0d10) then
               if (maswrk) write (iout,260)  ec,nec
  260          format (' Charge-Charge',17x,f16.4,i14)
            else
               if (maswrk) write (iout,270)  ec,nec
  270          format (' Charge-Charge',17x,d16.4,i14)
            end if
         end if
         if (use_chgdpl .and. necd.ne.0) then
            if (abs(ecd) .lt. 1.0d10) then
               if (maswrk) write (iout,280)  ecd,necd
  280          format (' Charge-Dipole',17x,f16.4,i14)
            else
               if (maswrk) write (iout,290)  ecd,necd
  290          format (' Charge-Dipole',17x,d16.4,i14)
            end if
         end if
         if (use_dipole .and. ned.ne.0) then
            if (abs(ed) .lt. 1.0d10) then
               if (maswrk) write (iout,300)  ed,ned
  300          format (' Dipole-Dipole',17x,f16.4,i14)
            else
               if (maswrk) write (iout,310)  ed,ned
  310          format (' Dipole-Dipole',17x,d16.4,i14)
            end if
         end if
         if (use_mpole .and. nem.ne.0) then
            if (abs(em) .lt. 1.0d10) then
               if (maswrk) write (iout,320)  em,nem
  320          format (' Atomic Multipoles',13x,f16.4,i14)
            else
               if (maswrk) write (iout,330)  em,nem
  330          format (' Atomic Multipoles',13x,d16.4,i14)
            end if
         end if
         if (use_polar .and. nep.ne.0) then
            if (abs(ep) .lt. 1.0d10) then
               if (maswrk) write (iout,340)  ep,nep
  340          format (' Polarization',18x,f16.4,i14)
            else
               if (maswrk) write (iout,350)  ep,nep
  350          format (' Polarization',18x,d16.4,i14)
            end if
         end if
         if (use_rxnfld .and. ner.ne.0) then
            if (maswrk) write (iout,360)  er,ner
  360       format (' Reaction Field',16x,f16.4,i14)
         end if
         if (use_solv .and. nes.ne.0) then
            if (maswrk) write (iout,370)  es,nes
  370       format (' Macroscopic Solvation',9x,f16.4,i14)
         end if
         if (use_geom .and. neg.ne.0) then
            if (maswrk) write (iout,380)  eg,neg
  380       format (' Geometric Restraints',10x,f16.4,i14)
         end if
         if (use_extra .and. nex.ne.0) then
            if (maswrk) write (iout,390)  ex,nex
  390       format (' Extra Energy Terms',12x,f16.4,i14)
         end if
      end if
c
c     radius of gyration and moments of inertia
c
      if (doprops) then
         call gyrate (rg)
         if (maswrk) write (iout,400)  rg
  400    format (/,' Radius of Gyration :',14x,f12.4,' Angstroms')
         call inertia (1)
      end if
c
c     total electrical charge and dipole moment
c
      if (doelect) then
         if (maswrk) write (iout,410)  netchg
  410    format (/,' Total Electric Charge :',11x,f12.5,' Electrons')
         moment = sqrt(xdipole**2 + ydipole**2 + zdipole**2)
         if (maswrk) write (iout,420)  moment,xdipole,ydipole,zdipole
  420    format (/,' Total Dipole Moment :',13x,f12.4,' Debyes',
     &           /,4x,'X-Component',20x,f12.4,
     &           /,4x,'Y-Component',20x,f12.4,
     &           /,4x,'Z-Component',20x,f12.4)
         if (dielec .ne. 1.0d0) then
            if (maswrk) write (iout,430)  dielec
  430       format (/,' Dielectric Constant :',13x,f12.4)
            if (maswrk) write (iout,440)  netchg/sqrt(dielec)
  440       format (' Effective Total Charge :',10x,f12.5,' Electrons')
            if (maswrk) write (iout,450)  moment/sqrt(dielec)
  450       format (' Effective Dipole Moment :',9x,f12.4,' Debyes')
         end if
      end if
c
c     energy partitioning over the individual atoms
c
      if (doatom) then
         if (maswrk) write (iout,460)
  460    format (/,' Potential Energy Breakdown over Atoms :')
         if (maswrk) write (iout,470)
  470    format (/,'  Atom',10x,'EB',10x,'EA',9x,'EBA',9x,'EUB',
     &               9x,'EAA',8x,'EOPB',
     &           /,15x,'EID',9x,'EIT',10x,'ET',9x,'EBT',9x,'ETT',
     &               9x,'E14',
     &           /,16x,'EV',10x,'EC',9x,'ECD',10x,'ED',10x,'EM',
     &               8x,'EP',
     &           /,16x,'ER',10x,'ES',10x,'EG',10x,'EX')
         do i = 1, n
            if (active(i)) then
               if (maswrk) write (iout,480)  i,aeb(i),aea(i),aeba(i),
     &                     aeub(i),aeaa(i),aeopb(i),aeid(i),aeit(i),
     &                     aet(i),aebt(i),aett(i),ae14(i),
     &                     aev(i),aec(i),aecd(i),aed(i),aem(i),
     &                     aep(i),aer(i),aes(i),aeg(i),aex(i)
  480          format (/,i6,6f12.4,/,6x,6f12.4,/,6x,6f12.4,/,6x,4f12.4)
            end if
         end do
      end if
c
c     list parameters used for molecular mechanics atom types
c
      if (doparam) then
         header = .true.
         do i = 1, n
            if (active(i)) then
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,490)
  490             format (/,' Atom Type Definition Parameters :',
     &                    //,3x,'Atom',4x,'Symbol',3x,'Type',
     &                       2x,'Class',2x,'Atomic',4x,'Mass',
     &                       3x,'Valence',2x,'Description',/)
               end if
               if (maswrk) write (iout,500) i,name(i),type(i),class(i),
     &                    atomic(i),mass(i),valence(i),story(i)
  500          format (i6,7x,a3,3i7,f11.3,i6,5x,a20)
            end if
         end do
c
c     list parameters used for van der Waals interactions
c
         if (use_vdw) then
            header = .true.
            do i = 1, n
               if (active(i)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,510)
  510                format (/,' Van der Waals Parameters :',
     &                       //,11x,'Atom Number',13x,'Radius',
     &                          3x,'Epsilon',3x,'Reduction',/)
                  end if
                  j = class(i)
                  if (reduct(j) .eq. 0.0d0) then
                     if (maswrk) write (iout,520)  i,i,rad(j),eps(j)
  520                format (i6,4x,i6,15x,2f10.4)
                  else
                     if (maswrk) write (iout,530)  i,i,rad(j),eps(j),
     &                   reduct(j)
  530                format (i6,4x,i6,15x,3f10.4)
                  end if
               end if
            end do
         end if
c
c     list parameters used for bond stretching interactions
c
         if (use_bond) then
            header = .true.
            do i = 1, nbond
               ia = ibnd(1,i)
               ib = ibnd(2,i)
               if (active(ia) .or. active(ib)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,540)
  540                format (/,' Bond Stretching Parameters :',
     &                       //,11x,'Atom Numbers',25x,'KS',
     &                          5x,'Length',/)
                  end if
                  if (maswrk) write (iout,550)  i,ia,ib,bk(i),bl(i)
  550             format (i6,4x,2i6,13x,f16.4,f10.4)
               end if
            end do
         end if
c
c     list parameters used for angle bending interactions
c
         if (use_angle) then
            header = .true.
            do i = 1, nangle
               ia = iang(1,i)
               ib = iang(2,i)
               ic = iang(3,i)
               if (active(ia) .or. active(ib) .or. active(ic)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,560)
  560                format (/,' Angle Bending Parameters :',
     &                       //,14x,'Atom Numbers',22x,'KB',
     &                          6x,'Angle',/)
                  end if
                  if (angin(i)) then
                     if (maswrk) write (iout,570)  i,ia,ib,ic,acon(i),
     &                     anat(i)
  570                format (i6,4x,3i6,7x,f16.4,f10.4,2x,'In-Plane')
                  else
                     if (maswrk) write (iout,580)  i,ia,ib,ic,acon(i),
     &                       anat(i)
  580                format (i6,4x,3i6,7x,f16.4,f10.4)
                  end if
               end if
            end do
         end if
c
c     list parameters used for stretch-bend interactions
c
         if (use_strbnd) then
            header = .true.
            do i = 1, nstrbnd
               j = isb(1,i)
               ia = iang(1,j)
               ib = iang(2,j)
               ic = iang(3,j)
               if (active(ia) .or. active(ib) .or. active(ic)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,590)
  590                format (/,' Stretch-Bend Parameters :',
     &                       //,14x,'Atom Numbers',11x,'KSB',
     &                          6x,'Angle',3x,'Length1',
     &                          3x,'Length2',/)
                  end if
                  if (maswrk) write (iout,600)  i,ia,ib,ic,ksb(i),
     &                        anat(j),bl(isb(2,i)),bl(isb(3,i))
  600             format (i6,4x,3i6,f13.4,3f10.4)
               end if
            end do
         end if
c
c     list parameters used for Urey-Bradley interactions
c
         if (use_urey) then
            header = .true.
            do i = 1, nurey
               ia = iury(1,i)
               ib = iury(2,i)
               if (active(ia) .or. active(ib)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,610)
  610                format (/,' Urey-Bradley Parameters :',
     &                       //,11x,'Atom Numbers',24x,'KUB',
     &                          4x,'Distance',/)
                  end if
                  if (maswrk) write (iout,620)  i,ia,ib,uk(i),ul(i)
  620             format (i6,4x,2i6,13x,f16.4,f10.4)
               end if
            end do
         end if
c
c     list parameters used for out-of-plane bending interactions
c
         if (use_opbend) then
            header = .true.
            do i = 1, nopbend
               j = iopb(i)
               ia = iang(1,j)
               ib = iang(2,j)
               ic = iang(3,j)
               id = iang(4,j)
               if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,630)
  630                format (/,' Out-of-Plane Bending Parameters :',
     &                       //,17x,'Atom Numbers',19x,'KOPB',/)
                  end if
                  if (maswrk) write (iout,640)  i,id,ib,ia,ic,kopb(i)
  640             format (i6,3x,4i6,9x,f10.4)
               end if
            end do
         end if
c
c     list parameters used for improper dihedral interactions
c
         if (use_improp) then
            header = .true.
            do i = 1, niprop
               ia = iiprop(1,i)
               ib = iiprop(2,i)
               ic = iiprop(3,i)
               id = iiprop(4,i)
               if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,650)
  650                format (/,' Improper Dihedral Parameters :',
     &                       //,17x,'Atom Numbers',19x,'KID',
     &                          4x,'Dihedral',/)
                  end if
                  if (maswrk) write (iout,660)  i,ia,ib,ic,id,kprop(i),
     &                  vprop(i)
  660             format (i6,3x,4i6,9x,2f10.4)
               end if
            end do
         end if
c
c     list parameters used for improper torsion interactions
c
         if (use_imptor) then
            header = .true.
            do i = 1, nitors
               ia = iitors(1,i)
               ib = iitors(2,i)
               ic = iitors(3,i)
               id = iitors(4,i)
               if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,670)
  670                format (/,' Improper Torsion Parameters :',
     &                       //,17x,'Atom Numbers',11x,
     &                          'Amplitude, Phase and Periodicity',/)
                  end if
                  j = 0
                  if (itors1(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 1
                     ampli(j) = itors1(1,i)
                     phase(j) = itors1(2,i)
                  end if
                  if (itors2(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 2
                     ampli(j) = itors2(1,i)
                     phase(j) = itors2(2,i)
                  end if
                  if (itors3(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 3
                     ampli(j) = itors3(1,i)
                     phase(j) = itors3(2,i)
                  end if
                  if (j .eq. 0) then
                     if (maswrk) write (iout,680)  i,ia,ib,ic,id
  680                format (i6,3x,4i6)
                  else if (j .eq. 1) then
                     if (maswrk) write (iout,690)  i,ia,ib,ic,id,
     &                                 ampli(1),phase(1),fold(1)
  690                format (i6,3x,4i6,10x,f10.3,f8.1,i4)
                  else if (j .eq. 2) then
                     if (maswrk) write (iout,700)  i,ia,ib,ic,id,
     &                        (ampli(k),phase(k),fold(k),k=1,j)
  700                format (i6,3x,4i6,2x,2(f10.3,f6.1,i4))
                  else
                     if (maswrk) write (iout,710)  i,ia,ib,ic,id,
     &                         (ampli(k),nint(phase(k)),fold(k),k=1,j)
  710                format (i6,3x,4i6,4x,3(f8.3,i4,'/',i1))
                  end if
               end if
            end do
         end if
c
c     list parameters used for torsional interactions
c
         if (use_tors) then
            header = .true.
            do i = 1, ntors
               ia = itors(1,i)
               ib = itors(2,i)
               ic = itors(3,i)
               id = itors(4,i)
               if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,720)
  720                format (/,' Torsional Angle Parameters :',
     &                       //,17x,'Atom Numbers',11x,
     &                          'Amplitude, Phase and Periodicity',/)
                  end if
                  j = 0
                  if (tors1(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 1
                     ampli(j) = tors1(1,i)
                     phase(j) = tors1(2,i)
                  end if
                  if (tors2(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 2
                     ampli(j) = tors2(1,i)
                     phase(j) = tors2(2,i)
                  end if
                  if (tors3(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 3
                     ampli(j) = tors3(1,i)
                     phase(j) = tors3(2,i)
                  end if
                  if (tors4(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 4
                     ampli(j) = tors4(1,i)
                     phase(j) = tors4(2,i)
                  end if
                  if (tors5(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 5
                     ampli(j) = tors5(1,i)
                     phase(j) = tors5(2,i)
                  end if
                  if (tors6(1,i) .ne. 0.0d0) then
                     j = j + 1
                     fold(j) = 6
                     ampli(j) = tors6(1,i)
                     phase(j) = tors6(2,i)
                  end if
                  if (j .eq. 0) then
                     if (maswrk) write (iout,730)  i,ia,ib,ic,id
  730                format (i6,3x,4i6)
                  else
                     if (maswrk) write (iout,740)  i,ia,ib,ic,id,
     &                          (ampli(k),nint(phase(k)),fold(k),k=1,j)
  740                format (i6,3x,4i6,4x,6(f8.3,i4,'/',i1))
                  end if
               end if
            end do
         end if
c
c     list parameters used for stretch-torsion interactions
c
         if (use_strtor) then
            header = .true.
            do i = 1, nstrtor
               j = ist(1,i)
               ia = itors(1,j)
               ib = itors(2,j)
               ic = itors(3,j)
               id = itors(4,j)
               if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,750)
  750                format (/,' Stretch-Torsion Parameters :',
     &                       //,17x,'Atom Numbers',10x,'Length',
     &                          5x,'Torsion Terms',/)
                  end if
                  k = 0
                  if (kst(1,i) .ne. 0.0d0) then
                     k = k + 1
                     ampli(k) = kst(1,i)
                     phase(k) = nint(tors1(2,j))
                  end if
                  if (kst(2,i) .ne. 0.0d0) then
                     k = k + 1
                     ampli(k) = kst(2,i)
                     phase(k) = 2 * nint(tors2(2,j))
                  end if
                  if (kst(3,i) .ne. 0.0d0) then
                     k = k + 1
                     ampli(k) = kst(3,i)
                     phase(k) = 3 * nint(tors3(2,j))
                  end if
                  if (maswrk) write (iout,760)  i,ia,ib,ic,id,
     &                   bl(ist(2,i)),(ampli(j),phase(j),j=1,k)
  760             format (i6,3x,4i6,2x,f10.4,1x,3(f8.3,i3))
               end if
            end do
         end if
c
c     list parameters used for atomic partial charges
c
         if (use_charge .or. use_chgdpl) then
            header = .true.
            do i = 1, nion
               ia = iion(i)
               ib = jion(i)
               ic = kion(i)
               if (active(ia) .or. active(ic)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,770)
  770                format (/,' Atomic Partial Charge Parameters :',
     &                       /,46x,'Neighbor',3x,'Cutoff',
     &                       /,11x,'Atom Number',13x,'Charge',
     &                          7x,'Site',6x,'Site',/)
                  end if
                  if (ia.eq.ib .and. ia.eq.ic) then
                     if (maswrk) write (iout,780)  i,ia,pchg(i)
  780                format (i6,4x,i6,15x,f10.4)
                  else
                     if (maswrk) write (iout,790)  i,ia,pchg(i),ib,ic
  790                format (i6,4x,i6,15x,f10.4,5x,i6,4x,i6)
                  end if
               end if
            end do
         end if
c
c     list parameters used for bond dipole moments
c
         if (use_dipole .or. use_chgdpl) then
            header = .true.
            do i = 1, ndipole
               ia = idpl(1,i)
               ib = idpl(2,i)
               if (active(ia) .or. active(ib)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,800)
  800                format (/,' Bond Dipole Moment Parameters :',
     &                       //,11x,'Atom Numbers',22x,'Dipole',
     &                          3x,'Position',/)
                  end if
                  if (maswrk) write (iout,810)  i,ia,ib,bdpl(i),sdpl(i)
  810             format (i6,4x,2i6,13x,f16.4,f10.4)
               end if
            end do
         end if
c
c     list parameters used for atomic multipole moments
c
         if (use_mpole) then
            header = .true.
            do i = 1, npole
               ia = ipole(i)
               if (active(ia)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,820)
  820                format (/,' Atomic Multipole Parameters :',
     &                       //,11x,'Atom Number',3x,'Local Axes',
     &                          ' Definition',10x,'Multipole Moments',/)
                  end if
                  do j = 2, 4
                     pole(j,i) = pole(j,i) / bohr
                  end do
                  do j = 5, 13
                     pole(j,i) = 3.0d0 * pole(j,i) / bohr**2
                  end do
        if (maswrk) write (iout,830)  i,ia,zaxis(i),xaxis(i),polaxe(i),
     &                              (pole(j,i),j=1,5),pole(8,i),
     &                              pole(9,i),(pole(j,i),j=11,13)
  830             format (i6,4x,i6,6x,i6,1x,i6,3x,a8,4x,f9.5,
     &                       /,50x,3f9.5,/,50x,f9.5,
     &                       /,50x,2f9.5,/,50x,3f9.5)
               end if
            end do
         end if
c
c     list parameters used for dipole polarizability
c
         if (use_polar) then
            header = .true.
            do i = 1, npole
               ia = ipole(i)
               if (active(ia)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,840)
  840                format (/,' Dipole Polarizability Parameters :',
     &                       //,11x,'Atom Number',14x,'Alpha',/)
                  end if
                  if (maswrk) write (iout,850)  i,ia,polarize(i)
  850             format (i6,4x,i6,15x,f10.4)
               end if
            end do
         end if
c
c     list parameters used for empirical solvation
c
         if (use_solv) then
            header = .true.
            do i = 1, n
               if (active(i)) then
                  if (header) then
                     header = .false.
                     if (maswrk) write (iout,860)
  860                format (/,' Empirical Solvation Parameters :',
     &                       //,11x,'Atom Number',13x,'Radius',
     &                          3x,'ASP Value',/)
                  end if
             if (maswrk) write (iout,870)  i,i,rsolv(i),vsolv(i)
  870             format (i6,4x,i6,15x,2f10.4)
               end if
            end do
         end if
      end if
c
c     perform any final tasks before program exit
c
cjrs we don't want to stop the program here
c      call final
cjrs 
c
c
cjrs added
c     write out final function value and gradient
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
         end if
      end do
      nvar = nvar*3
c
      call gradient (energy,derivs)
      gnorm = 0.0d0
cjrs
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               gnorm = gnorm + derivs(j,i)**2
cjrs
            end do
         end if
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(dble(nvar/3))
c
      if (maswrk) write(iout,900)
      do 45 i=1,n
        if (maswrk) write(iout,901) derivs(1,i),derivs(2,i),derivs(3,i)
 45   continue
       
      if (grms .gt. 0.0001d0) then
         if (maswrk) write (iout,40)  energy,grms,gnorm
   40    format (/,'            MM Energy (kcal/mol) : ',f15.4,
     &           /,'            MM RMS Gradient      : ',f15.4,
     &           /,'            MM Gradient Norm     : ',f15.4,/)
      else
         if (maswrk) write (iout,50)  energy,grms,gnorm
   50    format (/,'            MM Energy (kcal/mol) : ',f15.4,
     &           /,'            MM RMS Gradient      : ',d15.4,
     &           /,'            MM Gradient Norm     : ',d15.4,/)
      end if
      if (maswrk) write(iout,902)
c
ckk -fmoimomm- begin
      ENGMM=energy
      nd=1
      do i=1,n
         do  j=1,3
            tegtmp(nd)=derivs(j,i)
            nd=nd+1
         enddo
      enddo
ckk -fmoimomm- end
c
 900  format(25x,' **** TINKER Gradient (kcal/mol/ang) **** ',/,
     #     15x,'dE/dx ',15x,'dE/dy',15x,'dE/dz')
 901   format(10x,f12.8,10x,f12.8,10x,f12.8)
 902   format(/)
          
cjrs
      return      
      end
c
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  program newton  --  perform TNCG Cartesian optimization  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "newton" performs an energy minimization in Cartesian
c     coordinate space using a truncated Newton method
c
c
cjrs
c      program newton
cjrs
ckk -minimize-      subroutine tnewtx(mode,method,grdmin)
      subroutine tnewtx(optprg,mode,method,mzmeth,grdmin,maxhess)
      implicit double precision(a-h,o-z)
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      include 'usage.i'
casa -restrain- begin
      include 'restrn.i'
      include 'potent.i'
casa -restrain- end
cjrs
c -- New COMMON between TOYS and Tinker only
      real*8 teg
!      integer MxAtm
!      real*8 Energy,Eg(*)
      PARAMETER (MXATM=2000)
casa -restrain- begin
      PARAMETER (MAXR1=2000)
C -- MAXIMUM SIZES and COMMONS FOR QMMM LINKING
c
      logical IMOMM,SIMOMM
      COMMON /QMMM1/ IMOMM,SIMOMM,NPAIR,NSEQ
      COMMON /QMMM2/ IQMATM(MAXR1),ibasfmo(MAXR1)
casa -restrain- end
      COMMON /FUNCT / ENERGY,EG(3*MXATM)
      COMMON /TGRAD/ TEG(3*MAXATM)
      integer nd
cjrs
ckk -fmoimomm- begin
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody
ckk -fmoimomm- end
      integer i,j,imin,nvar
      integer next,freeunit
      real*8 gnorm,grms,grdmin
      real*8 newton1,minimum
      real*8 xx(maxvar),derivs(3,maxatm)
      character*1 answer
      character*6 mode,method
ckk -minimize- begin
      character*3 mzmeth
      double precision newtn,minimz
      data newtn /8HNEWTN   /,minimz/8HMINIMZ  /
ckk -minimize- end
      character*60 minfile
      character*80 record,string
      logical exist
      external newton1,newton2,writeout
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim
      integer IR,IW,IP,IS,IPK,IDAF,NAV,ioda
      integer NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
c     integer NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX 
      real*8 XXX 
casa -restrain- begin
      LOGICAL use_old(MAXATM)
      integer nuse_old
casa -restrain- end
c
      COMMON /FMCOM / XXX(1)
c     COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
c
c     Begin by allocating memory
c
c     maxhessa=maxhess
      if(optprg.eq.newtn) then
         if(maxhess.gt.0) then
            CALL VALFM(LOADFM)
            lh=LOADFM+1
            lhindex=lh+maxhess
            lf=lhindex+(maxhess-1)/nwdvar+1
            lcindex=lf+maxhess
            lcvalue=lcindex+(maxhess-1)/nwdvar+1
            last=lcvalue+(maxhess-1)/nwdvar+1
            NEED = LAST- LOADFM -1
            CALL GETFM(NEED)
            if(maswrk) write(iw,*)'Allocating',NEED,' words for Hessian'
         else
            call abrt
c           not implemented yet
c           maxhess=-maxhess
         endif
      endif
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i)
            nvar = nvar + 1
            xx(nvar) = y(i)
            nvar = nvar + 1
            xx(nvar) = z(i)
         end if
      end do
c
c     make the call to the optimization routine
c
ckk -minimize- begin
      if(optprg.eq.newtn) then
         call tncg (mode,method,nvar,xx,minimum,grdmin,
     &              newton1,newton2,writeout,
     &              xxx(lh),xxx(lhindex),
     &              xxx(lf),xxx(lcindex),xxx(lcvalue),maxhess)
      elseif(optprg.eq.minimz) then
c        lmqn expects mzmeth to be lower case!
         call locase(mzmeth)
         call lmqn (mzmeth,nvar,xx,minimum,grdmin,newton1,writeout)
      else
c        error exit
         if(maswrk) write(*,9000) ' error (tnewtx): wrong optprg:',
     &                            optprg
         call abrt
      endif
ckk -minimize- end
c     maxhess=maxhessa
c     if(maxhess.gt.0) then
      if(maxhess.gt.0.and.optprg.eq.newtn) then
         CALL RETFM(NEED)
      endif
c
c     untranslate the final coordinates for active atoms
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     write out final function value and gradient
c
casa -restrain- begin
      if (nfg .gt. 0) then
          nfix = npfix - npfixm
          do i=1, nfix
             ipfix(npfix) = 0
             xpfix(npfix) = 0
             ypfix(npfix) = 0
             zpfix(npfix) = 0
             pfix(npfix) =  0
             npfix = npfix - 1
          end do
          if( npfixm .eq. 0) then
              use_geom = .false.
          end if
c     store the use information
          nuse_old = nuse
          nuse = 0
          do i=1, n
             use_old(i) = use(i)
             use(i) = .true.
             nuse = nuse+1
          end do
c     turn off use option for QM atoms
          do i=1,nseq
             use(iqmatm(i))=.false.
             nuse=nuse-1
          enddo
          if (maswrk) write(iw,940) (iqmatm(i),i=1,nseq)
      end if
casa -restrain- end
      call gradient (minimum,derivs)
      gnorm = 0.0d0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               gnorm = gnorm + derivs(j,i)**2
            end do
         end if
      end do
casa -restrain- begin
      if (nfg .gt. 0) then
          nuse = nuse_old
          do i=1, n
             use(i)= use_old(i)
          enddo
          if (maswrk) write(iw,940) (iqmatm(i),i=1,nseq)
      end if
casa -restrain- end
C  CHC
C      write(*,*)
c      write(*,*) ' Tinker Gradient Components (kcal/mole/ang) '
c      do 200 i=1,n
c 200        write(*,*) derivs(1,i),derivs(2,i),derivs(3,i)
c      write(*,*)
c
c Put Tinker derivs into TEG array to pass back to TOYS
c
      nd=1
      do 800 i=1,n
         do 810  j=1,3
            teg(nd)=derivs(j,i)
 810        nd=nd+1
 800  continue
c             
cjrs
c
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(dble(nvar/3))
C     Adding the two energy
ckk -fmoimomm-      Energy=Energy+minimum/627.5095d+00
      if(nfg.eq.0) then
         Energy=Energy+minimum/627.5095d+00
      else
         Energy=minimum
      endif
ckk -fmoimomm-  if (grms .gt. 0.0001d0) then
         if (maswrk) write (iout,80)  minimum,grms,gnorm
ckk -fmoimomm- begin
         if (maswrk . and . nfg.eq.0) write(iout,81) Energy
   80    format (/,'           MM Energy (kcal/mol)  : ',f20.10,
     &           /,'           MM RMS Gradient       : ',f20.10,
     &           /,'           MM Gradient Norm      : ',f20.10)
   81    format (  '           QM+MM Energy (Hartree): ',f20.10,/)
ckk   else
ckk      if (maswrk) write (iout,90)  minimum,grms,gnorm,Energy
ckk90    format (/,'           MM Energy (kcal/mol)  : ',f20.10,
ckk  &           /,'           MM RMS Gradient       : ',f20.10,
ckk  &           /,'           MM Gradient Norm      : ',f20.10,
ckk  &           /,'           QM+MM Energy (Hartree): ',f20.10,/)
ckk   end if
ckk -fmoimomm- end
c
c     write the final coordinates into a file
c
cjrs redirect this to the GAMESS output file
c
c      imin = freeunit ()
c      open (unit=imin,file=minfile,status='old')
c      rewind (unit=imin)
c
      if (maswrk) call prtxyz (iout)
c
c      close (unit=imin)
c
c     perform any final tasks before program exit
c
c  -- Tinker final renamed finalt
c     call finalt
cjrs  don't forget, this is a subroutine now
c
      return
ckk -minimize-
 9000 format(1x,a32,5x,a6)
casa -restrain- 
  940 FORMAT(6(7x,i4))
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function newton1  --  energy/gradient values for newton  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "newton1" is a service routine that computes the energy
c     and gradient for truncated Newton optimization in Cartesian
c     coordinate space
c
c
      function newton1 (xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      integer i,nvar
      real*8 newton1,e,derivs(3,maxatm)
      real*8 xx(maxvar),g(maxvar)
c
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     compute and store the energy and gradient
c
      call gradient (e,derivs)
      newton1 = e
c
c     store atom gradients as optimization gradient, also
c     translate the coordinates of each active atom; the
c     latter may be needed when using periodic boundaries
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i)
            g(nvar) = derivs(1,i)
            nvar = nvar + 1
            xx(nvar) = y(i)
            g(nvar) = derivs(2,i)
            nvar = nvar + 1
            xx(nvar) = z(i)
            g(nvar) = derivs(3,i)
         end if
      end do
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine newton2  --  Hessian values for newton  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "newton2" is a service routine that computes the sparse
c     matrix Hessian elements for truncated Newton optimization
c     in Cartesian coordinate space
c
c
      subroutine newton2 (mode,xx,h,hinit,hstop,hindex,hdiag,maxhess)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      integer i,j,k,nvar,maxhess
      integer hinit(maxvar),hstop(maxvar)
      integer hindex(maxhess)
      integer hvar(maxvar),huse(maxvar)
      real*8 xx(maxvar),hdiag(maxvar),h(maxhess)
      character*4 mode
c
c
c     translate optimization parameters to atomic coordinates
c
      if (mode .eq. 'none')  return
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     compute and store the Hessian elements
c
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     transform the sparse Hessian to use only active atoms
c
      nvar = 0
      if (nuse .ne. n) then
         do i = 1, n
            k = 3 * (i-1)
            if (use(i)) then
               do j = 1, 3
                  nvar = nvar + 1
                  hvar(nvar) = j + k
                  huse(j+k) = nvar
               end do
            else
               do j = 1, 3
                  huse(j+k) = 0
               end do
            end if
         end do
         do i = 1, nvar
            k = hvar(i)
            hinit(i) = hinit(k)
            hstop(i) = hstop(k)
            hdiag(i) = hdiag(k)
            do j = hinit(i), hstop(i)
               hindex(j) = huse(hindex(j))
            end do
         end do
      end if
c
c     translate the coordinates of each active atom;
c     this may be needed when using periodic boundaries
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i)
            nvar = nvar + 1
            xx(nvar) = y(i)
            nvar = nvar + 1
            xx(nvar) = z(i)
         end if
      end do
      return
      end
