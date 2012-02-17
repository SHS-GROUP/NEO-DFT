      program dk3
      implicit double precision(a-h,o-z)
      character*1 type(4),typ
      character*2 symbol
      character*6 filenm
      character*50 title
      parameter (mxzeta=40, mxao=300)
      dimension e(mxzeta,4),vec(mxao,mxao),nocc(4),nzeta(4)
      parameter (tol=1.0d-7)
c
c     prepare DK3 basis sets for H-Lr from the optimized exponents
c     available from Professor Hirao's web page,
c        http://www.riken.jp/qcl/publications/dk3bs/periodic_table.html
c                 this web link was updated September 2009.
c     The formal literature citation for these basis sets is
c       T.Tsuchiya, M.Abe, T.Nakajima, K.Hirao
c       J.Chem.Phys. 115, 4463-4472(2001).
c
c     extraction program written by Mike Schmidt in March 2004.
c
c     a) download the web page for your element, in its entirety,
c     keeping all the lines at the top, to a file called xx.web.
c     Use the 1 or 2 digit chemical symbol for xx, e.g. Au.web.
c     Tell your browser to use 'text' mode when saving the web page,
c     and if it interprets this as including non-Unix return/newline
c     combinations (take a peek at the file with vi to see if it
c     contains any ^M characters), if so remove them, perhaps by
c             perl -pi -e 's/\r\n?/\n/g' xx.web
c     b) run this program, in mode 1, to make an uncontracted input
c     c) run GAMESS, after inserting the correct input from REFS.DOC
c     into the template input file to produce a wavefunction with
c     spherically averaged density, preserving degenerate p,d,f AOs.
c     d) place the $VEC from the GAMESS run into a file named xx.vec
c     e) run this program, in mode 2, to prepare a generally
c     contracted input file, containing a minimal basis set (MBS).
c     f) run GAMESS, again after inserting the correct spherical
c     averaging information from REFS.DOC.
c     check: the energy from steps d and f should be the same as the
c     web page at Tokyo Daigaku, and in step f, the final orbitals
c     should be a unit matrix.  Some idea of the agreement is
c        S     -398.589026149    U.Tokyo           scftyp=gvb mult=3
c              -398.5890261401   uncontracted
c              -398.5890261400   MBS contraction
c        Ga   -1942.234245329                      scftyp=gvb mult=2
c             -1942.2342464819
c             -1942.2342464820
c        Hg  -19626.258988046                      scftyp=rhf mult=1
c            -19626.2592464115
c            -19626.2592464133
c     g) edit the input file contracted to a MBS to include additional
c     Gaussians, to give the basis set additional variational freedom
c     in molecular calculations.  This basis seems well suited to the
c     development of a valence DZ quality contraction, by means of
c     adding the most diffuse primitive in each symmetry type as a new
c     shell, e.g. if the smallest s function is zeta=0.030, add
c         s 1 ; 1 0.030 1.0
c     In the case of TM elements, you probably want to float two of
c     the d's, in order to have a TZ quality valence d, a wise thing
c     to do even if you are trying for DZ quality otherwise.  TM
c     elements will benefit from even-tempered extensions of the p
c     exponents to cover the valence p orbital, not occupied in the
c     web page basis sets.  Add polarization functions to taste.
c     h) it's time to do molecules!
c
      in=1
      iout=2
      ir=5
      iw=6
c
      write(iw,900)
  900 format('This program runs in two modes:'/
     *       '   1 - read exponents from web page,',
     *       ' prepare uncontracted $data'/
     *       '   2 - read web page plus the $VEC matrix,',
     *       ' prepare contracted $data'/
     *       'please enter mode: ',$)
      read(ir,*) mode
c
      write(iw,910)
  910 format('enter name of element: ',$)
      read(ir,fmt='(a2)') symbol
      iend=2
      if(symbol(2:2).eq.' ') iend=1
      kend=6
      if(symbol(2:2).eq.' ') kend=5
c
      filenm = symbol(1:iend)//'.web'
      open(unit=in, file=filenm(1:kend), status='old',
     *     form='formatted', access='sequential', err=800)
      filenm = symbol(1:iend)//'.out'
      open(unit=iout, file=filenm(1:kend), status='unknown',
     *     form='formatted', access='sequential', err=800)
c
      read(in,8005) title
      read(in,8000)
      read(in,8000)
      read(in,8000) ns,np,nd,nf
      read(in,8000) ncs,ncp,ncd,ncf
      read(in,8000) nos,nop,nod,nof
      read(in,8000) nes,nep,ned,nef
      read(in,8000)
      read(in,8010) energy
      read(in,8000)
 8000 format(34x,i4,3i6)
 8005 format(2x,a50)
 8010 format(6x,f20.10)
 8020 format(4x,e16.7)
c
      ne = 2*ncs + 6*ncp + 10*ncd + 14*ncf + nes + nep + ned + nef
      znuc = ne
      nbf = ns + 3*np + 6*nd + 10*nf
      if(nbf.gt.mxao) then
         write(iw,*) 'please reset MXAO in program to',nbf
         stop
      end if
      type(1)='s'
      type(2)='p'
      type(3)='d'
      type(4)='f'
      nzeta(1)=ns
      nzeta(2)=np
      nzeta(3)=nd
      nzeta(4)=nf
      nocc(1)=ncs+nos
      nocc(2)=ncp+nop
      nocc(3)=ncd+nod
      nocc(4)=ncf+nof
c
      write(iout,720) energy,title,symbol,znuc
c
                       nirrep=1
      if(nocc(2).gt.0) nirrep=2
      if(nocc(3).gt.0) nirrep=3
      if(nocc(4).gt.0) nirrep=4
c
      do i=1,nirrep
         read(in,8000)
         read(in,8000)
         do j=1,nzeta(i)
            read(in,8020) e(j,i)
         enddo
      enddo
c
c         dump a $data group for an uncontracted calculation
c
      if(mode.eq.1) then
         do i=1,nirrep
            do j=1,nzeta(i)
               write(iout,730) type(i),e(j,i)
            enddo
         enddo
         go to 790
      endif
  920 format('the uncontracted input is now in file ',a6)
c
      if(mode.ne.2) then
         write(iw,*) 'can''t understand mode=',mode
         stop
      end if
c
c         read a $vec group prepared by executing a mode=1 job,
c         and create a general contraction-type $data
c
      close(unit=in, status='keep')
      filenm = symbol(1:iend)//'.vec'
      open(unit=in, file=filenm(1:kend), status='old',
     *     form='formatted', access='sequential', err=800)
      nbf = ns + 3*np + 6*nd + 10*nf
      nao = ns + 3*np + 5*nd +  7*nf
      call rdgms(in,vec,nao,nbf,mxao)
c
      ms =    ns
      mp =  3*np
      md =  6*nd
      mf = 10*nf
      needs = nocc(1)
      needp = nocc(2)
      needd = nocc(3)
      needf = nocc(4)
c
      do 290 iao=1,nao
c
c           assume clean $VEC given, so that if we find a non-zero
c           s,p,d,f, that's exactly what it is.  Furthermore, since
c           the first input should have been run in D2h, we assume
c           that we get clean px,py,pz (etc) subspecies too.  And
c           finally, we assume spherical harmonics are used, so there
c           is no s contamination in d, or p in f orbitals.  So we
c           can just look for the first non-zero coefficient, and
c           let this define the entire AO's symmetry type.
c
         do ibf=1,nbf
            if(abs(vec(ibf,iao)).gt.tol) go to 210
         enddo
         stop 'entire AO is zero?'
c
  210    continue
         if(ibf.le.ms)                                  typ='s'
         if(ibf.gt.ms        .and.  ibf.le.ms+mp)       typ='p'
         if(ibf.gt.ms+mp     .and.  ibf.le.ms+mp+md)    typ='d'
         if(ibf.gt.ms+mp+md  .and.  ibf.le.ms+mp+md+mf) typ='f'
c
         if(typ.eq.'s'  .and.  needs.gt.0) then
            needs=needs-1
            write(iout,750) typ,ns
            do i=1,ns
               write(iout,760) i,e(i,1),vec(i,iao)
            enddo
         end if
c
c            looking specifically for p-x subspecies
c
         if(typ.eq.'p'  .and.  needp.gt.0  .and.  ibf.eq.ms+1) then
            needp=needp-1
            write(iout,750) typ,np
            do i=1,np
               write(iout,760) i,e(i,2),vec(ms+1+3*(i-1),iao)
            enddo
         end if
c
c            d irrep of Kh resolves to 2Ag+B1g+B2g+B3g in D2h
c            looking specifically for d-xy subspecies (B1g symmery)
c
         if(typ.eq.'d'  .and.  needd.gt.0  .and.  ibf.eq.ms+mp+4) then
            needd=needd-1
            write(iout,750) typ,nd
            do i=1,nd
               write(iout,760) i,e(i,3),vec(ms+mp+4+6*(i-1),iao)
            enddo
         end if
c
c            f irrep of Kh resolves to Au+2B1u+2B2u+2B3u in D2h
c            looking specifically for f-xyz subspecies (Au symmetry)
c
         if(typ.eq.'f'  .and.  needf.gt.0 .and. ibf.eq.ms+mp+md+10) then
            needf=needf-1
            write(iout,750) typ,nf
            do i=1,nf
               write(iout,760) i,e(i,4),vec(ms+mp+md+10+10*(i-1),iao)
            enddo
         end if
  290 continue
c
c        all done, end the $data group, and exit.
c
  790 continue
      write(iout,740)
      filenm = symbol(1:iend)//'.out'
      write(iw,920) filenm
      stop
c
  800 continue
      write(6,990) filenm(1:6)
  990 format('error opening file ',a6)
      stop
c
  720 format('!       energy on Todai''s web page=',f20.10/
     *       ' $contrl scftyp=xxx mult=xxx runtyp=energy',
     *       ' relwfn=dk ispher=1 $end'/
     *       ' $system kdiag=3 mwords=5 $end'/
     *       ' $relwfn norder=3 modeqr=9 qmttol=1.0d-09 $end'/
     *       ' $data'/a50/'Dnh 2'/' '/
     *       a2,2x,f4.1,'    0.0 0.0 0.0')
  730 format(3x,a1,' 1 ; 1 ',1p,e13.7,' 1.0')
  740 format(' '/' $end')
  750 format(2x,a,i4)
  760 format(2x,i4,1x,1p,e13.7,1x,e15.8)
      end
c
      SUBROUTINE rdgms(in,V,M,N,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*5 keywd
      DIMENSION V(NDIM,*)
C
C     ----- READ IN $VEC -----
c
c     position to the $VEC group
C
  100 continue
      read(in,9000,end=800) keywd
      if(keywd.eq.' $VEC'  .or.  keywd.eq.' $vec') go to 200
      go to 100
c
  200 continue
      DO 220 J = 1,M
      IC = 0
      MAX = 0
  210 MIN = MAX+1
      MAX = MAX+5
      IC = IC+1
      IF (MAX .GT. N) MAX = N
      MODJ=MOD(J,100)
      READ(in,9008) MODJ,IC,(V(I,J),I = MIN,MAX)
      IF (MAX .LT. N) GO TO 210
  220 CONTINUE
      RETURN
c
  800 continue
      write(6,*) 'No $VEC card was found in your input file.'
      write(6,*) 'please check it.'
      stop
c
 9000 format(A5)
 9008 FORMAT(I2,I3,1P,5E15.8)
      END
