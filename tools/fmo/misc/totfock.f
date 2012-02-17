c----------------------------------------
c      Molecular orbitals from FMO/F 
c         (and FMO/FX, FMO/XF)
c           totfock program 
c       last updated Oct 28 2009  
c   to be used to process output from
c        GAMESS FMO 3.3 and up
c        coded by D.G.Fedorov
c          RICS, AIST, Japan
c----------------------------------------
c
c     Usage: run FMO in GAMESS with $fmoprp mofock=7 $end .
c     If it was a GDDI parallel run, concatenate all punch files
c     (file naming depends upon the details in rungms, you should make sure
c     that all punch files are copied back to your home directory!)
c     cat job.dat job.F07.* > joball.dat
c     If you fail to collect punch files from all groups, nonsense will ensue
c     when running totfock.x . 
c     Do not pipeline to totfock.x : it will work, but VERY slowly (why?!).
c     Compile this program (see comptot). Because of the comment below
c     numerious problems shall arise, no solution ever available! 
c     To run, provide a concatenated punch file:
c     totfock.x <joball.dat >joball.res
c
c     This program is written in a most atrocious mixture of F77 and F90,
c     truly abhorable to the bone to be used in classroom as a stunning example
c     of an explamplary poor programming.
c
      program totfock
      implicit double precision (A-H,O-Z)
      character*17 aread,apatt
      character*1  jobz,range
      logical doe,dumpth
      integer, dimension(:), allocatable :: l,ian,ianb,ic,isuppz,iwork,
     *                                      iworks,lpart,map,ianbp
      double precision, dimension(:), allocatable :: a,t,h,E,V,s,SE,SV,
     *                                            work,works,apart,bflab
c
      apatt= 'TOTAL FOCK MATRIX'
c     ir=5
c     iw=6
c     jobz='N'
      jobz='V'
      range='A'
      vl=0 
      vu=0
      il=0
      iu=0
c     For computing orbitals
      abstol=1.0D-10
      ispher=1
      qmttol=1.0d-06
c     If the total energy is to be computed (requires jobz.eq.'V').
      jl=0
      ju=0
c     For printing orbitals; if 0,0, nothing is printed, if -1,1, all is.
c
      ndimv=1
      ndimv2=1
  100 continue
        READ(*,9000,END=120,ERR=100) aread,natfmo,ichfmo,nefmo,mulfmo,
     *                         l0fmo,l1fmo,nbody,nlayer,nfg,maxl1,modpar
        if(aread.eq.apatt) then
          write(*,9010) natfmo,ichfmo,nefmo,mulfmo,l0fmo,l1fmo,nbody,
     *                  nlayer,nfg,maxl1,modpar
          if(jobz.eq.'V') ndimv=l1fmo
          dumpth=iand(modpar,4).ne.0
          doe=jobz.eq.'V'.and.dumpth
c         modpar:
c         1 set if dumping Fock from FMO enabled.
c         2 if exchange is added to Fock post factum.
c         4 if T+H (kinetic+1e Hamiltonian) are dumped, needed to do FMO/F.
          l2fmo=(l1fmo*l1fmo+l1fmo)/2
          if(doe) ndimv2=l2fmo
          lwork=26*l1fmo
          liwork=10*l1fmo
          lworks=l1fmo*l1fmo+6*l1fmo+1
          liworks=5*l1fmo+3
          maxl1d=maxl1*nbody
          maxl2d=(maxl1d*maxl1d+maxl1d)/2
          nafmo=nefmo/2
c         i.e., assuming singlet RHF
c         open shells can be added using mulfmo/nefmo.
          n=l1fmo
          n2=l2fmo
          allocate (a(n2),t(ndimv2),h(ndimv2),E(n),V(ndimv*ndimv),s(n2),
     *              SE(n),SV(n*n),l(n),ian(natfmo),ianb(n),ic(n),
     *              isuppz(2*n),work(lwork),iwork(liwork),works(lworks),
     *              iworks(liworks),bflab(ndimv),lpart(maxl1d),
     *              map(maxl1d),apart(maxl2d),ianbp(maxl1d))
c
          call diag(l1fmo, l2fmo, jobz, range, vl, vu, il, iu, abstol,
     *              ispher, qmttol, doe, dumpth, ndimv, ndimv2, lwork,
     *              liwork,lworks,liworks,jl,ju,nbody,nfg,maxl1d,
     *              maxl2d,nafmo,natfmo,e1,ef,a,t,h,E,V,s,SE,SV,l,ian,
     *              ianb,ic,isuppz,work,iwork,works,iworks,bflab,lpart,
     *              map,apart,ianbp)
c         deallocate(apart,map,lpart,iworks,works,iwork,work,isuppz,ic,
c    *               l,SV,SE,s,V,E,h,t,a)
          deallocate(a,t,h,E,V,s,SE,SV,l,ian,ianb,ic,isuppz,work,iwork,
     *               works,iworks,bflab,lpart,map,apart,ianbp)
          goto 120
        endif
      goto 100
  120 continue
 9000 format(1x,A17,I8,I4,I9,I3,2I9,I2,I2,I6,I6,I3)
 9010 format(1x,'FOUND FOCK MATRIX',I8,I4,I9,I3,2I9,I2,I2,I6,I6,I3)
      end

c call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z,
c ldz, isuppz, work, lwork, iwork, liwork, info)
c     dspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info)

      subroutine diag(n, n2, jobz, range, vl, vu, il, iu, abstol,ispher,
     *                qmttol, doe, dumpth, ndimv, ndimv2, lwork,liwork,
     *                lworks,liworks,jl,ju,nbody,nfg,maxl1,maxl2,nafmo,
     *                natfmo,e1,ef,a,t,h,E,V,s,SE,SV,l,ian,ianb,ic,
     *                isuppz,work,iwork,works,iworks,bflab,lpart,map,
     *                apart,ianbp)
      implicit double precision (A-H,O-Z)
      parameter (MAXL=6,ZERO=0.0D+00,ONE=1.0D+00,toeV=27.21138386D+00)
      character*120 str 
      character*20 aread,apatt,apatt2
      character*1  jobz,range
      logical doe,dumpth
      dimension a(n2),t(ndimv2),h(ndimv2),E(n),V(ndimv,ndimv),s(n2),
     *          SE(n),SV(n,n),l(n),ian(natfmo),ianb(n),ic(n),isuppz(2*n)
     *         ,work(lwork),iwork(liwork),works(lworks),iworks(liworks),
     *          bflab(ndimv)
      dimension lpart(maxl1),map(maxl1),apart(maxl2),ianbp(maxl1)
      dimension ln(0:MAXL),lm(0:MAXL),ts(1),tp(3,3),td(6,6),tf(10,10),
     *          tg(15,15),th(21,21),ti(28,28),r(28*28,0:MAXL)
c     save a,h,E,V,s,SE,SV,l,ic,isuppz,work,iwork,works,iworks,
c    *     lpart,map,apart
      DATA ln/1,3,6,10,15,21,28/
      DATA lm/1,3,5,7,9,11,13/
      DATA ts/1/
      DATA tp/1,0,0, 0,1,0, 0,0,1/
      DATA td/-0.5D+00,-0.5D+00,1, 0,0,0,
     *         0.866025403784439D+00,-0.866025403784439D+00,0,0,0,0,
     *         0,0,0,1,0,0,
     *         0,0,0,0,1,0,
     *         0,0,0,0,0,1,
     *         0,0,0,0,0,0/
      data tf/0.790569415042095D+00,0,0,0,0,
     *        -1.06066017177982D+00,0,0,0,0,
     *        0,0,0,0,0.866025403784438D+00,
     *        0,-0.866025403784438D+00,0,0,0,
     *        -0.612372435695794D+00,0,0,0,0,
     *        -0.273861278752583D+00,0,1.09544511501033D+00,0,0,
     *        0,0,1,0,-0.670820393249937D+00,
     *        0,-0.670820393249937D+00,0,0,0,
     *        0,-0.612372435695794D+00,0,-0.273861278752583D+00,0,
     *        0,0,0,1.09544511501033D+00,0,
     *        0,0,0,0,0, 0,0,0,0,1,
     *        0,-0.790569415042095D+00,0,1.06066017177982D+00,0,
     *        0,0,0,0,0,
     *        0,0,0,0,0, 0,0,0,0,0,
     *        0,0,0,0,0, 0,0,0,0,0,
     *        0,0,0,0,0, 0,0,0,0,0/
c
c     a(1,1)=0
c     a(2,1)=1
c     a(2,2)=0
c
c     First, assemble matrices a (F), s (S) and h (1e H)
c
      call dcopy(n2,0.0D+00,0,a,1)
      call dcopy(n2,0.0D+00,0,s,1)
      call dcopy(ndimv2,0.0D+00,0,t,1)
      call dcopy(ndimv2,0.0D+00,0,h,1)
      do i=1,n
        l(i)=-1
      enddo
      call dcopy(ln(0)*ln(0),ts,1,r(1,0),1)
      call dcopy(ln(1)*ln(1),tp,1,r(1,1),1)
      call dcopy(ln(2)*ln(2),td,1,r(1,2),1)
      call dcopy(ln(3)*ln(3),tf,1,r(1,3),1)
      call dcopy(ln(4)*ln(4),tg,1,r(1,4),1)
      call dcopy(ln(5)*ln(5),th,1,r(1,5),1)
      call dcopy(ln(6)*ln(6),ti,1,r(1,6),1)
      apatt='FRAGMENT FOCK MATRIX'
      apatt2='TOTAL FOCK NUCLEI='
      nblock=0
  100 continue
        READ(*,9000,END=120,ERR=100) str 
c       write(*,*) 'read',str,'end'
        READ(str,9005,END=200,ERR=200) aread,icurlay,icurfg,jcurfg,
     *                                 kcurfg,l1,ifac
c       if came here, then read was succesful and a new n-mer was found
        goto 210
  200   continue
c       Try looking for the other string 
        READ(str,9300,END=100,ERR=100) aread,enucr0
        if(aread.eq.apatt2) then
          enucr=enucr0
          write(*,9310) enucr
          read(*,9040,END=120,ERR=100) (ian(i),i=1,natfmo)
        endif 
c
        goto 100
  210   continue
        l2=(l1*l1+l1)/2
        if(aread.eq.apatt) then
          read(*,9020,END=120,ERR=100) (map(i),i=1,l1)
          if(jcurfg.eq.0.and.kcurfg.eq.0) then
            read(*,9025,END=120,ERR=100) (lpart(i),i=1,l1)
            read(*,9315,END=120,ERR=100) (ianbp(i),i=1,l1)
          endif
          write(*,9010) icurlay,icurfg,jcurfg,kcurfg,l1,ifac
          if(ifac.ne.0) nblock=nblock+1
          nm=1
c         The order is: S
          if(ifac.ne.0) then
            nm=2
c           The order is: S+F
            if(dumpth) nm=4
c           The order is: S+T+H+F
          endif
c         loop (e.g.) over S (m=1), T (m=2), H (m=3) and F (m=4)
          do m=1,nm
            read(*,9030,END=120,ERR=100) (apart(i),i=1,l2)
            loop=0
            do i=1,l1
              im=map(i)
              if(jcurfg.eq.0.and.kcurfg.eq.0) then
                l(im)=lpart(i)
                ianb(im)=ianbp(i)
              endif
              do j=1,i
                loop=loop+1
                jm=map(j)
                if(im.ge.jm) then
                  ig=im
                  jg=jm
                else
                  ig=jm
                  jg=im
                endif
c               ig=max(im,jm) 
c               jg=min(im,jm) 
c               write(6,*) 'www',loop,ig,jg
                loopg=(ig*ig-ig)/2+jg
                if(m.eq.1) then
                  s(loopg)=apart(loop)
c                 overwrite S (not accumulate!)
                else if(m.eq.nm) then
                  a(loopg)=a(loopg)+ifac*apart(loop)
c                 Accumulate F
c                 a(ig,jg)=a(ig,jg)+ifac*apart(loop)
                else if(m.eq.2) then
                  t(loopg)=apart(loop)
c                 overwrite t (not accumulate!)
                else if(m.eq.3) then
                  h(loopg)=apart(loop)
c                 overwrite h (not accumulate!)
                endif
c               if(ig.ne.jg) a(ig,jg)=0
c               if(ig.gt.n) a(ig,jg)=0
              enddo
            enddo
          enddo
        endif
      goto 100
  120 continue
      nmonomer=nfg
      ndimer=(nfg*nfg-nfg)/2
      ntrimer=(nfg*nfg*nfg-3*nfg*nfg+2*nfg)/6
      ntotal=nmonomer
      if(nbody.gt.1) ntotal=ntotal+ndimer
      if(nbody.gt.2) ntotal=ntotal+ntrimer
      mblock=nblock-nmonomer
      mtotal=ntotal-nmonomer
      if(mtotal.eq.0) then 
        sparse=0
      else
        sparse=1.0D+02-(mblock*1.0D+02)/mtotal
      endif
      write(*,9100) nblock,ntotal,sparse
c     write(*,*) 'L values for the basis set are:' 
c     write(*,9025) l 
c     write(*,*) 'ian' 
c     write(*,9025) ian 
c     write(*,*) 'ianb' 
c     write(*,9025) ianb 
c     write(*,*) 'Fock matrix is'
c     write(*,9030) a
c     write(*,*) 'Overlaps are'
c     write(*,9030) s
      if(doe) then 
c       write(*,*) '1e ints are'
c       write(*,9030) h
      endif
      l1=n
      nsalc=l1
      if(jobz.eq.'V') call setlab(natfmo,l1,ln,l,ian,ianb,bflab)
      if(ispher.eq.1) call sphers(n,n2,l,s,ln,lm,r,nsalc,ic)
      l0=nsalc
c     write(*,9030) (s(i),i=1,(l0*l0+l0)/2)
      call dspevd('V', 'U', l0, s, SE, SV, n,
     *            works, lworks, iworks, liworks, info)
      neig=l0
      if(info.eq.0) then
        write(*,*) 'Diagonalisation of S was successful.'
        write(*,*) neig,' S eigenvalues found, matrix size=',n
        write(*,9200) (SE(i),i=1,neig)
      else
        write(*,*) 'Error occurred (see dspevd), info=',info
      endif
c     write(*,9030) ((SV(i,j),i=1,n),j=1,l0)
c
c     Construct Q=Xs^(-1/2), where SX=Xs
c
      do i=1,l0
        call dscal(n,1/sqrt(se(i)),sv(1,i),1)
      enddo
        write(*,*) 'Done rooting eigenvalues of S.' 
      if(ispher.eq.1) then
        call spherv(l0,l,sv,n,ln,lm,r)
        write(*,*) 'Done transforming back eignevectors of S.' 
      endif
      call TFTRI(s,a,SV,works(n*n+1),l0,l1,n)
      write(*,*) 'Done transforming F to orthogonal basis.'
      if(jl.eq.-1) jl=1 
      if(ju.eq.-1) ju=l0
c
c     Expand triangular Qt*F*Q (stored in S) to square matrix 
c
c     call dcopy(l0*l0,0.0D+00,0,works,1)
c     loop=0
c     do i=1,l0
c       do j=1,i
c         loop=loop+1
c         works(i+(j-1)*l0)=s(loop)
c       enddo
c     enddo
c     write(*,9030) s
c     write(*,9030) (works(i),i=1,n*n) 
c
c     call dsyevr(jobz, range, 'L', l0, works, l0, vl, vu, il, iu,
c    *            abstol, neig, E, V, ndimv, isuppz, work, lwork, iwork,
c    *            liwork, info)
c     dsyevr produces wrong eignevectors for matrix size more than 13!!
      call dspevd(jobz, 'U', l0, s, E, V, ndimv,
     *            works, lworks, iworks, liworks, info)
      neig=l0
      if(info.eq.0) then
        write(*,*) 'Diagonalisation of F was successful.'
        write(*,*) neig,' F eigenvalues found, matrix size=',l0,ndimv
        write(*,9200) (E(i),i=1,neig)
        if(range.eq.'A') WRITE(*,9210) nafmo,nafmo+1,
     *                                 (e(nafmo+1)-e(nafmo))*toeV
c       do j=1,neig
c         if(jobz.eq.'V') write(*,9030) (v(i,j),i=1,l0)
c       enddo
c     Eigenvectors of F are back-transformed by Q.
      if(jobz.eq.'V') then
        write(*,*) 'Transforming back MOs...'
        CALL DGEMM('N','N',l1,l0,l0,ONE,sv,n,v,ndimv,ZERO,works,l1)
c
c       print orbitals.
c
        if(jl.gt.0.and.ju.gt.0)
     *    call PREVS(bflab,works,E,min(ju-jl+1,neig),l1,l1,jl)
        if(doe) then
          write(*,*) 'Done, computing 1e expectation values...',nafmo
          tkin=0
          e1=0
          ef=0
          nocc=2
          do i=1,min(nafmo,neig)
            call TFTRI(e1i,h,works((i-1)*l1+1),works(n*n+1),1,l1,n)
            call TFTRI(t1i,t,works((i-1)*l1+1),works(n*n+1),1,l1,n)
c           2 is the occupation number
            e1=e1+e1i*nocc
            tkin=tkin+t1i*nocc
            ef=ef+e(i)*nocc
          enddo
          e2=(ef-e1)/2
          etot=(e1+ef)/2+enucr
c         write(*,9300) ef,e1,(ef-e1)/2,(e1+ef)/2
c         Strictly speaking, range should be 'A' for the complete energy,
c         else it would be a partial contribution (may be interesting too?)
          E2 = ETOT - E1 - ENUCR
          VNE = E1 - TKIN
          VNN = ENUCR
          VEE = E2
          VTOT = VNE + VNN + VEE
          VIRIAL = -VTOT/TKIN
          WRITE(*,9320) E1,E2,ENUCR,ETOT
          WRITE(*,9330) VEE,VNE,VNN,VTOT,TKIN,VIRIAL
        endif
      endif
      else
        write(*,*) 'Error occurred (see dsyevr), info=',info
      endif
      return
 9000 format(1x,A120)
 9005 format(A20,I4,3I6,2I6)
c9008 FORMAT(I2,I3,1P,5E15.8)
 9010 format(1x,'FRAGMENT FOCK MATRIX',I4,3I6,2I6)
 9020 FORMAT(10I8)
 9025 FORMAT(40I2)
 9030 FORMAT(4E20.13)
 9040 format(1x,26I3)
 9100 FORMAT(1x,'Read',I9,' blocks, max',I9,', sparsity=',F8.2,' %')
 9200 FORMAT(6F13.6)
 9210 format(/1x,'HOMO=',I8,' , LUMO=',I8,' , gap=',F10.3,' eV',/)
c9300 FORMAT(/1x,'Total energies are:',
c    *       /1x,'   Fock matrix=   ',F20.10,
c    *       /1x,'  one electron=   ',F20.10,
c    *       /1x,'  two electron=   ',F20.10,
c    *       /1x,'total electron=   ',F20.10)
 9300 format(A18,F25.13)
 9310 format(1x,'TOTAL FOCK NUCLEI=',F25.13)
 9315 format(10I8)
 9320 FORMAT(/1X,'               ONE ELECTRON ENERGY =',F24.10/
     *        1X,'               TWO ELECTRON ENERGY =',F24.10/
     *        1X,'          NUCLEAR REPULSION ENERGY =',F24.10/
     *       38X,23(1H-)/
     *        1X,'                      TOTAL ENERGY =',F24.10)
 9330 FORMAT(/1X,'ELECTRON-ELECTRON POTENTIAL ENERGY =',F24.10/
     *        1X,' NUCLEUS-ELECTRON POTENTIAL ENERGY =',F24.10/
     *        1X,'  NUCLEUS-NUCLEUS POTENTIAL ENERGY =',F24.10/
     *       38X,23(1H-)/
     *        1X,'            TOTAL POTENTIAL ENERGY =',F24.10/
     *        1X,'              TOTAL KINETIC ENERGY =',F24.10/
     *        1X,'                VIRIAL RATIO (V/T) =',F24.10)
      end
      SUBROUTINE TFTRI(H,F,T,WRK,M,N,LDT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION H(*),F(*),T(LDT,M),WRK(N)
      PARAMETER (MXROWS=5,zero=0.0D+00)
C
C     ----- TRANSFORM THE TRIANGULAR MATRIX F USING VECTORS T -----
C                      H = T-DAGGER * F * T
C     THE ORDER OF THE TRIANGULAR MATRICES H AND F ARE M AND N.
C
      IJ = 0
      DO 310 J = 1,M,MXROWS
         JJMAX = MIN(M,J+MXROWS-1)
C
C             FIRST CALCULATE T-DAGGER TIMES -F-, A ROW AT A TIME
C
         DO 300 JJ=J,JJMAX
            IK = 0
            DO 140 I = 1,N
               IM1 = I-1
               DUM = ZERO
               TDUM = T(I,JJ)
               IF (IM1.GT.0) THEN
                  DO 100 K = 1,IM1
                     IK = IK+1
                     WRK(K) = WRK(K)+F(IK)*TDUM
                     DUM = DUM+F(IK)*T(K,JJ)
  100             CONTINUE
               END IF
               IK = IK+1
               WRK(I) = DUM+F(IK)*TDUM
  140       CONTINUE
C
C             THEN TAKE THAT ROW TIMES EVERY COLUMN IN -T-
C
            DO 200 I = 1,JJ
               IJ = IJ+1
               HIJ = DDOT(N,T(1,I),1,WRK,1)
               H(IJ)=HIJ
  200       CONTINUE
  300    CONTINUE
  310 CONTINUE
C
      RETURN
      END
      SUBROUTINE sphers(n,n2,l,s,ln,lm,t,nsalc,ic)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter (MAXL=6,ZERO=0.0D+00,ONE=1.0D+00)
      dimension l(n),s(n2),ln(0:MAXL),lm(0:MAXL),t(28*28,0:MAXL)
      dimension sb1(28,28),sb2(28,28)
      dimension ic(n)
c
c     Transform S to the spherical basis, 
c     S'=Vt * S * V, where
c     V is the SALC matrix, block diagonal. Each AO forms a block
c     (one block for all components, such as px, py and pz).
c     V is not stored, instead blocks are taken from T indexed by l.
c     S' overwrites S.
c
      nsalc=0
      i=1
  100 continue
        li=l(i)
        lni=ln(li)
        lmi=lm(li)
        j=1
  200   continue
        lj=l(j)
        lnj=ln(lj)
        if(li.gt.1.or.lj.gt.1) then
c       otherwise T is a unit matrix
c         write(*,*) 'trans',i,j,li,lj,lni,lnj
          do im=1,lni
            loopi=i-1+im
            loopj=(loopi*loopi-loopi)/2+j-1
            maxj=lnj
            if(i.eq.j) maxj=im
            do jm=1,maxj
              loopj=loopj+1
              sb1(im,jm)=s(loopj)
c             write(6,*) im,jm,loopj,s(loopj)
            enddo
            if(i.eq.j) then
              loopj=(loopi*loopi-loopi)/2+j-1
              do jm=1,im-1
                loopj=loopj+1
                sb1(jm,im)=s(loopj)
c             write(6,*) jm,im,loopj,s(loopj)
              enddo
            endif
          enddo
c         if(i.eq.j) write(*,*) ((sb1(ii,jj),ii=1,lni),jj=1,lnj)
c         write(*,*) (t(ii,li),ii=1,lni*lni)
c         write(*,*) (t(ii,lj),ii=1,lnj*lnj)
c     DGEMM should be replaced by DCOPY for unit matrices
      CALL DGEMM('T','N',lni,lnj,lni,ONE,t(1,li),lni,sb1,28,ZERO,sb2,28)
c         write(*,*) ((sb2(ii,jj),ii=1,lni),jj=1,lnj)
      CALL DGEMM('N','N',lni,lnj,lnj,ONE,sb2,28,t(1,lj),lnj,ZERO,sb1,28)
c         write(*,*) ((sb1(ii,jj),ii=1,lni),jj=1,lnj)
c         Put the transformed block back into S 
          do im=1,lni
            loopi=i-1+im
            loopj=(loopi*loopi-loopi)/2+j-1
            maxj=lnj
            if(i.eq.j) maxj=im
            do jm=1,maxj
              loopj=loopj+1
              s(loopj)=sb1(im,jm)
            enddo
          enddo
        endif
        j=j+lnj
        if(j.le.i) goto 200
        do im=1,lmi
          ic(i+im-1)=0
          nsalc=nsalc+1
        enddo
        do im=lmi+1,lni
          ic(i+im-1)=1
        enddo
        i=i+lni
      if(i.le.n) goto 100
c     write(*,*) '15th row',(s(105+ii),ii=1,15)
c     write(*,*)'15th col',(s(((14+ii)*(14+ii)+(14+ii))/2+15),ii=1,n-14)
c
c     Now purge zero rows and columns (spherical contaminants)
c
c     write(*,9030) s
      write(*,*) 'Purging',n-nsalc,' contaminants...'
c     write(*,*) (ic(ii),ii=1,n) 
      loop1=0
      loop2=0
      do i=1,n
        do j=1,i
          loop1=loop1+1
          if(ic(i)+ic(j).eq.0) then
            loop2=loop2+1
            s(loop2)=s(loop1)
          endif
        enddo
      enddo
c
      RETURN
c9030 FORMAT(5E16.9)
      END
      SUBROUTINE spherv(l0,l,sv,n,ln,lm,t)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter (MAXL=6,ZERO=0.0D+00,ONE=1.0D+00)
      dimension l(n),sv(n,l0),ln(0:MAXL),lm(0:MAXL)
      dimension t(28*28,0:MAXL),sb(28,28)
c
c     Transform back S eigenvectors (in-place)
c     X = V * X'
c     V is the SALC matrix, block diagonal. Each AO forms a block
c     (one block for all components, such as px, py and pz).
c     V is not stored, instead blocks are taken from T indexed by l.
c
c     Trudge to find the last block 
c     write(*,9030) ((sv(i,j),j=1,l0),i=1,l0)
      i=1
      i0=1
  100 continue
        li=l(i)
        lni=ln(li)
        lmi=lm(li)
        i=i+lni
        i0=i0+lmi
      if(i0.le.l0) goto 100
c
  110 continue
c       This assumes that spherical components are consequent
        lip=l(i-1)
        i=i-ln(lip)
        i0=i0-lm(lip)
c       if(i0.lt.1) goto 300
        li=l(i)
        lni=ln(li)
        lmi=lm(li)
        j=1
  200   continue
          lj=l(j)
          lmj=lm(lj)
          if(li.gt.1) then
c           otherwise T is a unit matrix
c           write(*,*) 'trans',i,j,li,lj,lni,lmj
            CALL DGEMM('N','N',lni,lmj,lmi,ONE,t(1,li),lni,sv(i0,j),n,
     *                 ZERO,sb,28)
            do ii=lni-1,0,-1 
              do jj=j,j+lmj-1
                sv(i+ii,jj)=sb(ii+1,jj-j+1)
              enddo
            enddo
          else
c           in this case lni=lmi and the V block is a unit matrix
            if(i.ne.i0) then
              do ii=lmi-1,0,-1
                do jj=j,j+lmj-1
                  sv(i+ii,jj)=sv(i0+ii,jj)
                enddo
              enddo
            endif
c           write(*,*) 'dcopy',i0,i,j,li,lj,lni,lmi,lmj,n
          endif
          j=j+lmj
        if(j.le.l0) goto 200
c     goto 110
      if(i.gt.1) goto 110
c 300 continue
c     write(*,9030) ((t(i+(j-1)*6,2),j=1,l0),i=1,n)
c     write(*,*) 'Q matrix' 
c     do j=1,l0
c     write(*,9030) (sv(i,j),i=1,n)
c     enddo
c     write(iw,8030) ((q(i,j),j=l0,l0),i=1,l1)
c
      RETURN
c9030 FORMAT(6E13.5)
      END
C*MODULE MTHLIB  *DECK PREVS
      SUBROUTINE PREVS(bflab,V,E,NMO,NAO,LDV,ISTMO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(LDV,NMO),E(NMO),bflab(nao)
C
C     ----- PRINT OUT EIGENDATA, WITH MO SYMMETRY LABELS -----
C     THE ROWS ARE LABELED WITH THE BASIS FUNCTION NAMES.
C
      MAX = 5
      IMAX = ISTMO-1
C 
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. NMO) IMAX = NMO
      WRITE (*,9008) 
      if(NAO.le.99999) then
        WRITE (*,9028) (I,     I=IMIN,IMAX)
        WRITE (*,9068) (E(I),  I=IMIN,IMAX)
        WRITE (*,9078) ('A   ',I=IMIN,IMAX)
        DO J = 1,NAO
          WRITE (*,9048) J,BFLAB(J),(V(J,I),I = IMIN,IMAX)
        enddo 
      else
        WRITE (*,9029) (I,     I=IMIN,IMAX)
        WRITE (*,9069) (E(I),  I=IMIN,IMAX)
        WRITE (*,9079) ('A   ',I=IMIN,IMAX)
        DO J = 1,NAO
          WRITE (*,9049) J,BFLAB(J),(V(J,I),I = IMIN,IMAX)
        enddo 
      endif
      IF (IMAX .LT. NMO) GO TO 100
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(15X,10(4X,I4,3X))
 9048 FORMAT(I5,2X,A8,10F11.6)
 9068 FORMAT(15X,10F11.4)
 9078 FORMAT(16X,10(5X,A4,2X))
 9029 FORMAT(18X,10(4X,I4,3X))
 9049 FORMAT(I8,2X,A8,10F11.6)
 9069 FORMAT(18X,10F11.4)
 9079 FORMAT(19X,10(5X,A4,2X))
      END
      subroutine setlab(natfmo,l1,ln,l,ian,ianb,bflab)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter (MAXL=6)
      CHARACTER*6 BFNAM(35+49)
      CHARACTER*8 BFL  
      dimension l(l1),ian(natfmo),ianb(l1),bflab(l1),ln(0:MAXL),
     *          loff(0:MAXL)
      CHARACTER*2 ATMLAB(104)
      data BFNAM/'  s ','  x ','  y ','  z ',
     *           ' xx ',' yy ',' zz ',' xy ',' xz ',' yz ',
     *           ' xxx',' yyy',' zzz',' xxy',' xxz',
     *           ' yyx',' yyz',' zzx',' zzy',' xyz',
     *           'xxxx','yyyy','zzzz','xxxy','xxxz',
     *           'yyyx','yyyz','zzzx','zzzy','xxyy',
     *           'xxzz','yyzz','xxyz','yyxz','zzxy',
     *           ' xxxxx',' yyyyy',' zzzzz',' xxxxy',' xxxxz',
     *           ' yyyyx',' yyyyz',' zzzzx',' zzzzy',' xxxyy',
     *           ' xxxzz',' yyyxx',' yyyzz',' zzzxx',' zzzyy',
     *           ' xxxyz',' yyyxz',' zzzxy',' xxyyz',' xxzzy',
     *           ' yyzzx',
     *           '    x6','    y6','    z6','   x5y','   x5z',
     *           '   y5x','   y5z','   z5x','   z5y','  x4y2',
     *           '  x4z2','  y4x2','  y4z2','  z4x2','  z4y2',
     *           '  x4yz','  y4xz','  z4xy','  x3y3','  x3z3',
     *           '  y3z3',' x3y2z',' x3z2y',' y3x2z',' y3z2x',
     *           ' z3x2y',' z3y2x','x2y2z2'/
C
C     104 TRUE ELEMENTS
C
      DATA ATMLAB/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     *            'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     *            'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     *            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     *            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     *            'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     *            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     *            'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     *            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     *            'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     *            'Md','No','Lr','Rf'/
      data loff/1,2,5,11,21,36,57/
c
      i=1
  100 continue
        li=l(i)
        lni=ln(li)
        loffi=loff(li)
        do j=1,lni
          IF(li.LE.4) THEN
            iat=ianb(i)
            WRITE(UNIT=BFL,FMT='(A2,I2,A4)') 
     *        atmlab(ian(iat)),MOD(IAT,100),BFNAM(loffi+j-1)
          ELSE
            WRITE(UNIT=BFL,FMT='(A2,A6)') 
     *        atmlab(ian(iat)),BFNAM(loffi+j-1)
          END IF
          READ (UNIT=BFL,FMT='(A8)') BFLAB(i)
          i=i+1
        enddo
      if(i.le.l1) goto 100
      RETURN
      END
