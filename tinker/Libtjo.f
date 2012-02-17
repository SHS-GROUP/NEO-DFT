C 21 Apr 10 - NA,KK - add CHARGE2 support and allow choosing minimizer
C  9 MAR 00 - CHC - FIX for parallel run
C 14 Oct 98 - CHC - change title -> ttitle
c  8 May 98 - JRS - slater: renamed to slatert
c                   Note: slater bundled with overlap in Tinker
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine tnk_jacobi  --  jacobi matrix diagonalization  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "jacobi" performs a matrix diagonalization of a real
c     symmetric matrix by the method of Jacobi rotations
c
c     n    logical dimension of the matrix to be diagonalized
c     np   physical dimension of the matrix storage area
c     a    input with the matrix to be diagonalized; only
c             the upper triangle and diagonal are required
c     d    returned with the eigenvalues in ascending order
c     v    returned with the eigenvectors of the matrix
c     b    temporary work vector
c     z    temporary work vector
c
c
      subroutine tnk_jacobi (n,np,a,d,v,b,z)
      implicit none
      include 'iounit.i'
      integer i,j,k,ip,iq,n,np,nrot,maxrot
      real*8 sm,tresh,s,c,t,theta,tau,h,g,p
      real*8 a(np,np),d(np),v(np,np),b(np),z(np)
c
c
c     setup and initialization
c
      maxrot = 100
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
c
c     perform the jacobi rotations
c
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip))
     &                    .and. abs(d(iq))+g.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c = 1.0d0 / sqrt(1.0d0+t**2)
                  s = t * c
                  tau = s / (1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do
c
c     print warning if not converged
c
   10 continue
      if (nrot .eq. maxrot) then
         write (iout,20)
   20    format (/,' JACOBI  --  Matrix Diagonalization not Converged')
      end if
c
c     sort the eigenvalues and vectors
c
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kangang  --  angle-angle parameter assignment  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kangang" assigns the parameters for angle-angle cross term
c     interactions and processes new or changed parameter values
c
c
      subroutine kangang
      implicit none
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'kanang.i'
      include 'keys.i'
      include 'potent.i'
      integer i,j,k,m,next
      integer it,ia,ic
      integer nang,jang,kang
      real*8 fa,faa,aak(3)
      character*20 keyword
      character*80 record,string
      logical header
c
c
c     process keywords containing angle-angle parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'ANGANG ') then
            it = 0
            do j = 1, 3
               aak(j) = 0.0d0
            end do
            string = record(next:80)
            read (string,*,err=10,end=10)  it,(aak(j),j=1,3)
   10       continue
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Angle-Angle Parameters :',
     &                 //,5x,'Atom Class',8x,'K(AA) 1',4x,'K(AA) 2',
     &                    4x,'K(AA) 3',/)
            end if
            write (iout,30)  it,(aak(j),j=1,3)
   30       format (9x,i3,7x,3f11.3)
            do j = 1, 3
               anan(j,it) = aak(j)
            end do
         end if
      end do
c
c     assign the angle-angle parameters for each angle pair
c
      nangang = 0
      do i = 1, n
         nang = n12(i) * (n12(i)-1) / 2
         it = class(i)
         do j = 1, nang-1
            jang = anglist(j,i)
            ia = iang(1,jang)
            ic = iang(3,jang)
            m = 1
            if (atomic(ia) .le. 1)  m = m + 1
            if (atomic(ic) .le. 1)  m = m + 1
            fa = anan(m,it)
            do k = j+1, nang
               kang = anglist(k,i)
               ia = iang(1,kang)
               ic = iang(3,kang)
               m = 1
               if (atomic(ia) .le. 1)  m = m + 1
               if (atomic(ic) .le. 1)  m = m + 1
               faa = fa * anan(m,it)
               if (faa .ne. 0.0d0) then
                  nangang = nangang + 1
                  iaa(1,nangang) = jang
                  iaa(2,nangang) = kang
                  kaa(nangang) = faa
               end if
            end do
         end do
      end do
c
c     if no angle-angle terms are used, turn off the potential
c
      if (nangang .eq. 0)  use_angang = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kangle  --  angle bend parameter assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kangle" assigns the force constants and ideal angles for
c     the bond angles; also processes new or changed parameters
c
c
      subroutine kangle
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kangs.i'
      include 'keys.i'
      include 'potent.i'
      include 'ring.i'
      integer i,j,k,iring,next,size
      integer ia,ib,ic,id,ita,itb,itc
      integer jen,ih,nh,minatomic
      real*8 fc,an1,an2,an3
      character*3 pa,pb,pc
      character*9 pt
      character*20 keyword
      character*80 record,string
      logical header,use_ring3
      logical use_ring4,use_ring5
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing bond angle bending parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:6) .eq. 'ANGLE ')  iring = 0
         if (keyword(1:7) .eq. 'ANGLE5 ')  iring = 5
         if (keyword(1:7) .eq. 'ANGLE4 ')  iring = 4
         if (keyword(1:7) .eq. 'ANGLE3 ')  iring = 3
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            jen = 0
            string = record(next:80)
            read (string,*,err=10,end=10)  ia,ib,ic,fc,an1,an2,an3
   10       continue
            if (an2.ne.0.0d0 .or. an3.ne.0.0d0)  jen = 1
            if (header) then
               header = .false.
               if (maswrk) write (iout,20)
   20          format (/,' Additional Angle Bending Parameters :',
     &                 //,5x,'Atom Classes',9x,'K(B)',7x,'Angle',/)
            end if
            if (iring .eq. 0) then
               if (jen .eq. 0) then
                  if (maswrk) write (iout,30)  ia,ib,ic,fc,an1
   30             format (4x,3i4,2x,2f12.3)
               else if (an1 .ne. 0.0d0) then
                  if (maswrk) write (iout,40)  ia,ib,ic,fc,an1
   40             format (4x,3i4,2x,2f12.3,3x,'0-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  if (maswrk) write (iout,50)  ia,ib,ic,fc,an2
   50             format (4x,3i4,2x,2f12.3,3x,'1-H''s')
               end if
               if (an3 .ne. 0.0d0) then
                  if (maswrk) write (iout,60)  ia,ib,ic,fc,an3
   60             format (4x,3i4,2x,2f12.3,3x,'2-H''s')
               end if
            else if (iring .eq. 5) then
               if (jen .eq. 0) then
                  if (maswrk) write (iout,70)  ia,ib,ic,fc,an1
   70             format (4x,3i4,2x,2f12.3,3x,'5-Ring')
               else if (an1 .ne. 0.0d0) then
                  if (maswrk) write (iout,80)  ia,ib,ic,fc,an1
   80             format (4x,3i4,2x,2f12.3,3x,'5-Ring',3x,'0-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  if (maswrk) write (iout,90)  ia,ib,ic,fc,an2
   90             format (4x,3i4,2x,2f12.3,3x,'5-Ring',3x,'1-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  if (maswrk) write (iout,100)  ia,ib,ic,fc,an3
  100             format (4x,3i4,2x,2f12.3,3x,'5-Ring',3x,'2-H''s')
               end if
            else if (iring .eq. 4) then
               if (jen .eq. 0) then
                  if (maswrk) write (iout,110)  ia,ib,ic,fc,an1
  110             format (4x,3i4,2x,2f12.3,3x,'4-Ring')
               else if (an1 .ne. 0.0d0) then
                  if (maswrk) write (iout,120)  ia,ib,ic,fc,an1
  120             format (4x,3i4,2x,2f12.3,3x,'4-Ring',3x,'0-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  if (maswrk) write (iout,130)  ia,ib,ic,fc,an2
  130             format (4x,3i4,2x,2f12.3,3x,'4-Ring',3x,'1-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  if (maswrk) write (iout,140)  ia,ib,ic,fc,an3
  140             format (4x,3i4,2x,2f12.3,3x,'4-Ring',3x,'2-H''s')
               end if
            else if (iring .eq. 3) then
               if (jen .eq. 0) then
                  if (maswrk) write (iout,150)  ia,ib,ic,fc,an1
  150             format (4x,3i4,2x,2f12.3,3x,'3-Ring')
               else if (an1 .ne. 0.0d0) then
                  if (maswrk) write (iout,160)  ia,ib,ic,fc,an1
  160             format (4x,3i4,2x,2f12.3,3x,'3-Ring',3x,'0-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  if (maswrk) write (iout,170)  ia,ib,ic,fc,an2
  170             format (4x,3i4,2x,2f12.3,3x,'3-Ring',3x,'1-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  if (maswrk) write (iout,180)  ia,ib,ic,fc,an3
  180             format (4x,3i4,2x,2f12.3,3x,'3-Ring',3x,'2-H''s')
               end if
            end if
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            if (iring .eq. 0) then
               do j = 1, maxna
                  if (ka(j).eq.'         ' .or. ka(j).eq.pt) then
                     ka(j) = pt
                     con(j) = fc
                     ang(1,j) = an1
                     ang(2,j) = an2
                     ang(3,j) = an3
                     goto 230
                  end if
               end do
               if (maswrk) write (iout,190)
  190          format (/,' KANGLE  --  Too many Bond Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 5) then
               do j = 1, maxna5
                  if (ka5(j).eq.'         ' .or. ka5(j).eq.pt) then
                     ka5(j) = pt
                     con5(j) = fc
                     ang5(1,j) = an1
                     ang5(2,j) = an2
                     ang5(3,j) = an3
                     goto 230
                  end if
               end do
               if (maswrk) write (iout,200)
  200          format (/,' KANGLE  --  Too many 5-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 4) then
               do j = 1, maxna4
                  if (ka4(j).eq.'         ' .or. ka4(j).eq.pt) then
                     ka4(j) = pt
                     con4(j) = fc
                     ang4(1,j) = an1
                     ang4(2,j) = an2
                     ang4(3,j) = an3
                     goto 230
                  end if
               end do
               if (maswrk) write (iout,210)
  210          format (/,' KANGLE  --  Too many 4-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 3) then
               do j = 1, maxna3
                  if (ka3(j).eq.'         ' .or. ka3(j).eq.pt) then
                     ka3(j) = pt
                     con3(j) = fc
                     ang3(1,j) = an1
                     ang3(2,j) = an2
                     ang3(3,j) = an3
                     goto 230
                  end if
               end do
               if (maswrk) write (iout,220)
  220          format (/,' KANGLE  --  Too many 3-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            end if
  230       continue
         end if
      end do
c
c     set generic parameters for use with any number of hydrogens
c
      do j = 1, maxna
         if (ka(j) .eq. '         ')  goto 240
         if (ang(2,j).eq.0.0d0 .and. ang(3,j).eq.0.0d0) then
            ang(2,j) = ang(1,j)
            ang(3,j) = ang(1,j)
         end if
      end do
  240 continue
      do j = 1, maxna5
         if (ka5(j) .eq. '         ')  goto 250
         if (ang5(2,j).eq.0.0d0 .and. ang5(3,j).eq.0.0d0) then
            ang5(2,j) = ang5(1,j)
            ang5(3,j) = ang5(1,j)
         end if
      end do
  250 continue
      do j = 1, maxna4
         if (ka4(j) .eq. '         ')  goto 260
         if (ang4(2,j).eq.0.0d0 .and. ang4(3,j).eq.0.0d0) then
            ang4(2,j) = ang4(1,j)
            ang4(3,j) = ang4(1,j)
         end if
      end do
  260 continue
      do j = 1, maxna3
         if (ka3(j) .eq. '         ')  goto 270
         if (ang3(2,j).eq.0.0d0 .and. ang3(3,j).eq.0.0d0) then
            ang3(2,j) = ang3(1,j)
            ang3(3,j) = ang3(1,j)
         end if
      end do
  270 continue
c
c     use small rings if present and parameters are available
c
      if (nring3.eq.0 .or. ka3(1).eq.'         ') then
         use_ring3 = .false.
      else
         use_ring3 = .true.
      end if
      if (nring4.eq.0 .or. ka4(1).eq.'         ') then
         use_ring4 = .false.
      else
         use_ring4 = .true.
      end if
      if (nring5.eq.0 .or. ka5(1).eq.'         ') then
         use_ring5 = .false.
      else
         use_ring5 = .true.
      end if
c
c     loop over all angles assigning the necessary constants
c
      header = .true.
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         minatomic = min(atomic(ia),atomic(ib),atomic(ic))
         size = 3
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pt = pa//pb//pc
         else
            pt = pc//pb//pa
         end if
         acon(i) = 0.0d0
         anat(i) = 0.0d0
c
c     count number of non-angle hydrogens on the central atom
c
         nh = 1
         do j = 1, n12(ib)
            ih = i12(j,ib)
            if (ih.ne.ia .and. ih.ne.ic .and. atomic(ih).eq.1)
     &         nh = nh + 1
         end do
c
c     make a check for bond angles inside small rings
c
         iring = 0
         if (n12(ia).eq.1 .or. n12(ic).eq.1)  goto 280
         if (use_ring3) then
            do j = 1, n12(ia)
               if (ic .eq. i12(j,ia)) then
                  iring = 3
                  goto 280
               end if
            end do
         end if
         if (use_ring4) then
            do j = 1, n12(ia)
               id = i12(j,ia)
               if (ib .ne. id) then
                  do k = 1, n12(ic)
                     if (id .eq. i12(k,ic)) then
                        iring = 4
                        goto 280
                     end if
                  end do
               end if
            end do
         end if
         if (use_ring5) then
            do j = 1, n12(ia)
               id = i12(j,ia)
               if (ib .ne. id) then
                  do k = 1, n13(ic)
                     if (id .eq. i13(k,ic)) then
                        iring = 5
                        goto 280
                     end if
                  end do
               end if
            end do
         end if
  280    continue
c
c     assign angle bending parameters for each bond angle
c
         if (iring .eq. 0) then
            do j = 1, maxna
               if (ka(j).eq.pt .and. ang(nh,j).ne.0.0d0) then
                  acon(i) = con(j)
                  anat(i) = ang(nh,j)
                  goto 310
               end if
            end do
            if (minatomic .ne. 0) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,290)
  290             format (/,' Undefined Angle Bending Parameters :',
     &                    //,' Type',18x,'Atom Names',19x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,300)  ia,name(ia),ib,name(ib),
     &                           ic,name(ic),ita,itb,itc
  300          format (' Angle',6x,3(i6,'-',a3),7x,3i5)
            end if
  310       continue
c
c     assign bending parameters for 5-membered ring angles
c
         else if (iring .eq. 5) then
            do j = 1, maxna5
               if (ka5(j).eq.pt .and. ang5(nh,j).ne.0.0d0) then
                  acon(i) = con5(j)
                  anat(i) = ang5(nh,j)
                  goto 340
               end if
            end do
            if (minatomic .ne. 0) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,320)
  320             format (/,' Undefined Angle Bending Parameters :',
     &                    //,' Type',18x,'Atom Names',19x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,330)  ia,name(ia),ib,name(ib),
     &                           ic,name(ic),ita,itb,itc
  330          format (' 5-Ring',5x,3(i6,'-',a3),7x,3i5)
            end if
  340       continue
c
c     assign bending parameters for 4-membered ring angles
c
         else if (iring .eq. 4) then
            do j = 1, maxna4
               if (ka4(j).eq.pt .and. ang4(nh,j).ne.0.0d0) then
                  acon(i) = con4(j)
                  anat(i) = ang4(nh,j)
                  goto 370
               end if
            end do
            if (minatomic .ne. 0) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,350)
  350             format (/,' Undefined Angle Bending Parameters :',
     &                    //,' Type',18x,'Atom Names',19x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,360)  ia,name(ia),ib,name(ib),
     &                           ic,name(ic),ita,itb,itc
  360          format (' 4-Ring',5x,3(i6,'-',a3),7x,3i5)
            end if
  370       continue
c
c     assign bending parameters for 3-membered ring angles
c
         else if (iring .eq. 3) then
            do j = 1, maxna3
               if (ka3(j).eq.pt .and. ang3(nh,j).ne.0.0d0) then
                  acon(i) = con3(j)
                  anat(i) = ang3(nh,j)
                  goto 400
               end if
            end do
            if (minatomic .ne. 0) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,380)
  380             format (/,' Undefined Angle Bending Parameters :',
     &                    //,' Type',18x,'Atom Names',19x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,390)  ia,name(ia),ib,name(ib),
     &                           ic,name(ic),ita,itb,itc
  390          format (' 3-Ring',5x,3(i6,'-',a3),7x,3i5)
            end if
  400       continue
         end if
      end do
c
c     if no angle bend terms are used, turn off the potential
c
      if (nangle .eq. 0)  use_angle = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine katom  --  atom type parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "katom" assigns an atom type definitions to each atom in
c     the structure and processes any new or changed values
c
c
      subroutine katom
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'keys.i'
      integer i,k,next
      integer cls,atn,lig
      real*8 wght
      character*3 symb
      character*20 keyword,notice
      character*80 record,string
      logical header
      INTEGER me,master,nproc,ibtyp,iptim,iflag
      CHARACTER*8 GRPNAM
C
C
      LOGICAL LOG,FOUND,GOPARR,DSKWRK,MASWRK,TDSKWRK
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing atom type parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            k = 0
            cls = 0
            symb = ' '
            notice = ' '
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,k,next)
            call getnumb (record,cls,next)
            if (cls .eq. 0)  cls = k
            atmcls(k) = cls
            call gettext (record,symb,next)
            call getstring (record,notice,next)
            string = record(next:80)
            read (string,*,err=50,end=50)  atn,wght,lig
            if (k.ge.1 .and. k.le.maxtyp) then
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,20)
   20             format (/,' Additional Atom Type Parameters :',
     &                    //,4x,'Type  Class  Symbol',4x,'Description',
     &                       12x,'Atomic #',4x,'Mass',3x,'Valence',/)
               end if
               symbol(k) = symb
               describe(k) = notice
               atmnum(k) = atn
               weight(k) = wght
               ligand(k) = lig
              if (maswrk) write (iout,30) k,cls,symb,notice,atn,wght,lig
   30          format (2i7,5x,a3,5x,a20,i8,f12.3,i6)
            else if (k .ge. maxtyp) then
               if (maswrk) write (iout,40)
   40          format (/,' KATOM   --  Too many Atom Types;',
     &                    ' Increase MAXTYP')
               abort = .true.
            end if
   50       continue
         end if
      end do
c
c     transfer atom type values to individual atoms
c
      do i = 1, n
         k = type(i)
         if (k .eq. 0) then
            class(i) = 0
            atomic(i) = 0
            mass(i) = 0.0d0
            valence(i) = 0
            story(i) = 'Undefined Atom Type '
         else
casa -charge2-
            xyzname(i) = name(i)
            if (symbol(k) .ne. '   ')  name(i) = symbol(k)
            class(i) = atmcls(k)
            atomic(i) = atmnum(k)
            mass(i) = weight(k)
            valence(i) = ligand(k)
            story(i) = describe(k)
         end if
      end do
c
c     process keywords containing atom types for specific atoms
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            k = 0
            symb = ' '
            notice = ' '
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,k,next)
            call getnumb (record,cls,next)
            call gettext (record,symb,next)
            call getstring (record,notice,next)
            string = record(next:80)
            read (string,*,err=90,end=90)  atn,wght,lig
            if (k.lt.0 .and. k.ge.-n) then
               if (header) then
                  header = .false.
                  if (maswrk)  write (iout,70)
   70             format (/,' Additional Atom Types for',
     &                       ' Specific Atoms :',
     &                    //,4x,'Type  Class  Symbol',4x,'Description',
     &                       12x,'Atomic #',4x,'Mass',3x,'Valence',/)
               end if
               k = -k
               if (cls .eq. 0)  cls = k
               class(k) = cls
               name(k) = symb
               story(k) = notice
               atomic(k) = atn
               mass(k) = wght
               valence(k) = lig
              if (maswrk) write (iout,80) k,cls,symb,notice,atn,wght,lig
   80          format (2i7,5x,a3,5x,a20,i8,f12.3,i6)
            end if
   90       continue
         end if
      end do
c
c     check for undefined atom types in the molecule
c
      header = .true.
      do i = 1, n
         k = type(i)
         if (k.lt.1 .or. k.gt.maxtyp) then
            abort = .true.
            if (header) then
               header = .false.
               if (maswrk)  write (iout,100)
  100          format (/,' Undefined or Illegal Atom Types :',
     &                 //,' Type',10x,'Atom Number',5x,'Atom Type',/)
            end if
            if (maswrk)  write (iout,110)  i,k
  110       format (' Atom',12x,i5,10x,i5)
         end if
      end do
c
c     check the number of atoms attached to each atom
c
      header = .true.
      do i = 1, n
         if (n12(i) .ne. valence(i)) then
            if (header) then
               header = .false.
               if (maswrk)  write (iout,120)
  120          format (/,' Atoms with an Unusual Number of Attached',
     &                    ' Atoms :',
     &                 //,' Type',11x,'Atom Name',6x,'Atom Type',7x,
     &                    'Expected',4x,'Found',/)
            end if
            if (maswrk) write (iout,130) i,name(i),type(i),valence(i),
     *                         n12(i)
  130       format (' Valence',7x,i5,'-',a3,8x,i5,10x,i5,5x,i5)
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine kbond  --  bond parameter assignment  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "kbond" assigns a force constant and ideal bond length
c     to each bond in the structure and processes any new or
c     changed parameter values
c
c
      subroutine kbond
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kbonds.i'
      include 'keys.i'
      include 'potent.i'
      include 'ring.i'
      integer i,j,k,m
      integer ia,ib,ic,id
      integer ita,itb
      integer iring,size,next
      real*8 fc,bd
      character*3 pa,pb
      character*6 pt
      character*20 keyword
      character*80 record,string
      logical header,use_ring3
      logical use_ring4,use_ring5
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing bond stretching parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:5) .eq. 'BOND ')  iring = 0
         if (keyword(1:6) .eq. 'BOND5 ')  iring = 5
         if (keyword(1:6) .eq. 'BOND4 ')  iring = 4
         if (keyword(1:6) .eq. 'BOND3 ')  iring = 3
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:80)
            read (string,*,err=10,end=10)  ia,ib,fc,bd
   10       continue
            if (header) then
               header = .false.
               if (maswrk) write (iout,20)
   20          format (/,' Additional Bond Stretching Parameters :',
     &                 //,5x,'Atom Classes',9x,'K(S)',6x,'Length',/)
            end if
            if (iring .eq. 0) then
               if (maswrk) write (iout,30)  ia,ib,fc,bd
   30          format (6x,2i4,4x,f12.3,f12.4)
            else if (iring .eq. 5) then
               if (maswrk) write (iout,40)  ia,ib,fc,bd
   40          format (6x,2i4,4x,f12.3,f12.4,3x,'5-Ring')
            else if (iring .eq. 4) then
               if (maswrk) write (iout,50)  ia,ib,fc,bd
   50          format (6x,2i4,4x,f12.3,f12.4,3x,'4-Ring')
            else if (iring .eq. 3) then
               if (maswrk) write (iout,60)  ia,ib,fc,bd
   60          format (6x,2i4,4x,f12.3,f12.4,3x,'3-Ring')
            end if
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            if (iring .eq. 0) then
               do j = 1, maxnb
                  if (kb(j).eq.'      ' .or. kb(j).eq.pt) then
                     kb(j) = pt
                     fcon(j) = fc
                     blen(j) = bd
                     goto 80
                  end if
               end do
               if (maswrk) write (iout,70)
   70          format (/,' KBOND   --  Too many Bond Stretching',
     &                       ' Parameters')
               abort = .true.
   80          continue
            else if (iring .eq. 5) then
               do j = 1, maxnb5
                  if (kb5(j).eq.'      ' .or. kb5(j).eq.pt) then
                     kb5(j) = pt
                     fcon5(j) = fc
                     blen5(j) = bd
                     goto 100
                  end if
               end do
               if (maswrk) write (iout,90)
   90          format (/,' KBOND   --  Too many 5-Ring Stretching',
     &                       ' Parameters')
               abort = .true.
  100          continue
            else if (iring .eq. 4) then
               do j = 1, maxnb4
                  if (kb4(j).eq.'      ' .or. kb4(j).eq.pt) then
                     kb4(j) = pt
                     fcon4(j) = fc
                     blen4(j) = bd
                     goto 120
                  end if
               end do
               if (maswrk) write (iout,110)
  110          format (/,' KBOND   --  Too many 4-Ring Stretching',
     &                       ' Parameters')
               abort = .true.
  120          continue
            else if (iring .eq. 3) then
               do j = 1, maxnb3
                  if (kb3(j).eq.'      ' .or. kb3(j).eq.pt) then
                     kb3(j) = pt
                     fcon3(j) = fc
                     blen3(j) = bd
                     goto 140
                  end if
               end do
               if (maswrk) write (iout,130)
  130          format (/,' KBOND   --  Too many 3-Ring Stretching',
     &                       ' Parameters')
               abort = .true.
  140          continue
            end if
         end if
      end do
c
c     use small rings if present and parameters are available
c
      if (nring3.eq.0 .or. kb3(1).eq.'      ') then
         use_ring3 = .false.
      else
         use_ring3 = .true.
      end if
      if (nring4.eq.0 .or. kb4(1).eq.'      ') then
         use_ring4 = .false.
      else
         use_ring4 = .true.
      end if
      if (nring5.eq.0 .or. kb5(1).eq.'      ') then
         use_ring5 = .false.
      else
         use_ring5 = .true.
      end if
c
c     assign bond length and force constant for each bond
c
      header = .true.
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 3
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt = pa//pb
         else
            pt = pb//pa
         end if
         bk(i) = 0.0d0
         bl(i) = 0.0d0
c
c     make a check for bonds inside small rings
c
         iring = 0
         if (n12(ia).eq.1 .or. n12(ib).eq.1)  goto 150
         if (use_ring3) then
            do j = 1, n12(ia)
               ic = i12(j,ia)
               do k = 1, n12(ib)
                  if (ic .eq. i12(k,ib)) then
                     iring = 3
                     goto 150
                  end if
               end do
            end do
         end if
         if (use_ring4) then
            do j = 1, n12(ia)
               ic = i12(j,ia)
               if (ic .ne. ib) then
                  do k = 1, n12(ib)
                     id = i12(k,ib)
                     if (id .ne. ia) then
                        do m = 1, n12(ic)
                           if (id .eq. i12(m,ic)) then
                              iring = 4
                              goto 150
                           end if
                        end do
                     end if
                  end do
               end if
            end do
         end if
         if (use_ring5) then
            do j = 1, n13(ia)
               ic = i13(j,ia)
               do k = 1, n13(ib)
                  if (ic .eq. i13(k,ib)) then
                     iring = 5
                     goto 150
                  end if
               end do
            end do
         end if
  150    continue
c
c     assign bond stretching parameters for each bond
c
         if (iring .eq. 0) then
            do j = 1, maxnb
               if (kb(j) .eq. pt) then
                  bk(i) = fcon(j)
                  bl(i) = blen(j)
                  goto 180
               end if
            end do
            if (atomic(ia).ne.0 .and. atomic(ib).ne.0) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,160)
  160             format (/,' Undefined Bond Stretching Parameters :',
     &                    //,' Type',13x,'Atom Names',11x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,170)  ia,name(ia),ib,name(ib),
     &                ita,itb
  170          format (' Bond',7x,i6,'-',a3,i6,'-',a3,7x,2i5)
            end if
  180       continue
c
c     assign stretching parameters for 5-membered ring bonds
c
         else if (iring .eq. 5) then
            do j = 1, maxnb5
               if (kb5(j) .eq. pt) then
                  bk(i) = fcon5(j)
                  bl(i) = blen5(j)
                  goto 210
               end if
            end do
            if (atomic(ia).ne.0 .and. atomic(ib).ne.0) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,190)
  190             format (/,' Undefined Bond Stretching Parameters :',
     &                    //,' Type',13x,'Atom Names',11x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,200)  ia,name(ia),ib,name(ib),
     &                ita,itb
  200          format (' 5-Ring',5x,i6,'-',a3,i6,'-',a3,7x,2i5)
            end if
  210       continue
c
c     assign stretching parameters for 4-membered ring bonds
c
         else if (iring .eq. 4) then
            do j = 1, maxnb4
               if (kb4(j) .eq. pt) then
                  bk(i) = fcon4(j)
                  bl(i) = blen4(j)
                  goto 240
               end if
            end do
            if (atomic(ia).ne.0 .and. atomic(ib).ne.0) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,220)
  220             format (/,' Undefined Bond Stretching Parameters :',
     &                    //,' Type',13x,'Atom Names',11x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,230)  ia,name(ia),ib,name(ib),
     &               ita,itb
  230          format (' 4-Ring',5x,i6,'-',a3,i6,'-',a3,7x,2i5)
            end if
  240       continue
c
c     assign stretching parameters for 3-membered ring bonds
c
         else if (iring .eq. 3) then
            do j = 1, maxnb3
               if (kb3(j) .eq. pt) then
                  bk(i) = fcon3(j)
                  bl(i) = blen3(j)
                  goto 270
               end if
            end do
            if (atomic(ia).ne.0 .and. atomic(ib).ne.0) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,250)
  250             format (/,' Undefined Bond Stretching Parameters :',
     &                    //,' Type',13x,'Atom Names',11x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,260)  ia,name(ia),ib,name(ib),
     &                ita,itb
  260          format (' 3-Ring',5x,i6,'-',a3,i6,'-',a3,7x,2i5)
            end if
  270       continue
         end if
      end do
c
c     if no bond stretches are used, turn off the potential
c
      if (nbond .eq. 0)  use_bond = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kcharge  --  assign partial charge parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kcharge" assigns partial charges to the atoms within
c     the structure and processes any new or changed values
c
c
      subroutine kcharge
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kchrge.i'
      include 'keys.i'
      include 'potent.i'
      integer i,j,k,m,ia,next
      integer chglist(maxatm)
      integer nc12(maxatm)
      real*8 cg
      character*20 keyword
      character*80 record,string
      logical header
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing partial charge parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            string = record(next:80)
            read (string,*,err=40,end=40)  ia,cg
            if (ia .gt. 0) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Additional Atomic Partial Charge',
     &                       ' Parameters :',
     &                    //,6x,'Atom Type',9x,'Charge',/)
               end if
               if (ia .le. maxtyp) then
                  chg(ia) = cg
                  write (iout,20)  ia,cg
   20             format (6x,i6,6x,f12.3)
               else
                  write (iout,30)
   30             format (/,' KCHARGE  --  Too many Partial Charge',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
   40       continue
         end if
      end do
c
c     find and store all the atomic partial charges
c
casa  -charge2- begin
      do i = 1, n
         if ( type(i) .ge. 2021 ) then
            do j = 1, n_chg2-1
                if ( xyzname(i) .eq. atname(j) ) then
                    pchg(i) = chg2(j)
                end if
            end do
         else
             pchg(i) = chg(type(i))
         end if
         if (pchg(i) .eq. 0 ) then
             if (maswrk) then
                write(iout, 9000) i, xyzname(i)
 9000           format(/,'WARNING  -- ATOM', 1x, i5, '(', a4, ')', 
     &                      1x, 'has no charge!',/)
             end if
         end if
      end do
casa  -charge2- end
c
c     process keywords containing atom specific partial charges
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            string = record(next:80)
            read (string,*,err=70,end=70)  ia,cg
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Partial Charges for',
     &                       ' Specific Atoms :',
     &                    //,8x,'Atom',12x,'Charge',/)
               end if
               write (iout,60)  ia,cg
   60          format (6x,i6,6x,f12.3)
               pchg(ia) = cg
            end if
   70       continue
         end if
      end do
c
c     remove zero partial charges from the list of charges
c
      nion = 0
      do i = 1, n
         chglist(i) = 0
         if (pchg(i) .ne. 0.0d0) then
            nion = nion + 1
            iion(nion) = i
            jion(nion) = i
            kion(nion) = i
            pchg(nion) = pchg(i)
            chglist(i) = nion
         end if
      end do
c
c     optionally use neutral groups for neighbors and cutoffs
c
      if (neutnbr .or. neutcut) then
         do i = 1, nion
            nc12(i) = 0
            k = iion(i)
            do j = 1, n12(k)
               m = chglist(i12(j,k))
               if (m .ne. 0)  nc12(i) = nc12(i) + 1
            end do
         end do
         do i = 1, nion
            if (nc12(i) .eq. 1) then
               k = iion(i)
               do j = 1, n12(k)
                  m = chglist(i12(j,k))
                  if (m .ne. 0) then
                     if (nc12(m) .gt. 1) then
                        if (neutnbr)  jion(i) = iion(m)
                        if (neutcut)  kion(i) = iion(m)
                     end if
                  end if
               end do
            end if
         end do
      end if
c
c     if there are no partial charges, turn off the potential
c
      if (nion .eq. 0) then
         use_charge = .false.
         use_chgdpl = .false.
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kdipole  --  assign bond dipole parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kdipole" assigns bond dipoles to the bonds within
c     the structure and processes any new or changed values
c
c
      subroutine kdipole
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'dipole.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kdipol.i'
      include 'keys.i'
      include 'potent.i'
      include 'ring.i'
      integer i,j,k,m
      integer ia,ib,ic,id
      integer ita,itb
      integer iring,size,next
      real*8 dp,ps
      character*3 pa,pb
      character*6 pt
      character*20 keyword
      character*80 record,string
      logical header,use_ring3
      logical use_ring4,use_ring5
c
c
c     process keywords containing bond dipole parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:7) .eq. 'DIPOLE ')  iring = 0
         if (keyword(1:8) .eq. 'DIPOLE5 ')  iring = 5
         if (keyword(1:8) .eq. 'DIPOLE4 ')  iring = 4
         if (keyword(1:8) .eq. 'DIPOLE3 ')  iring = 3
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:80)
            read (string,*,err=10,end=10)  ia,ib,dp,ps
   10       continue
            if (ia.gt.0 .and. ib.gt.0) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Bond Dipole Moment ',
     &                       'Parameters :',
     &                    //,6x,'Atom Types',8x,'Moment',
     &                       6x,'Position',/)
               end if
               if (iring .eq. 0) then
                  write (iout,30)  ia,ib,dp,ps
   30             format (6x,2i4,4x,2f12.3)
               else if (iring .eq. 5) then
                  write (iout,40)  ia,ib,dp,ps
   40             format (6x,2i4,4x,2f12.3,3x,'5-Ring')
               else if (iring .eq. 4) then
                  write (iout,50)  ia,ib,dp,ps
   50             format (6x,2i4,4x,2f12.3,3x,'4-Ring')
               else if (iring .eq. 3) then
                  write (iout,60)  ia,ib,dp,ps
   60             format (6x,2i4,4x,2f12.3,3x,'3-Ring')
               end if
               size = 3
               call numeral (ia,pa,size)
               call numeral (ib,pb,size)
               if (ia .le. ib) then
                  pt = pa//pb
               else
                  pt = pb//pa
               end if
               if (iring .eq. 0) then
                  do j = 1, maxnd
                     if (kd(j).eq.'      ' .or. kd(j).eq.pt) then
                        kd(j) = pt
                        if (ia .le. ib) then
                           dpl(j) = dp
                           pos(j) = ps
                        else
                           dpl(j) = -dp
                           pos(j) = 1.0d0 - ps
                        end if
                        goto 110
                     end if
                  end do
                  write (iout,70)
   70             format (/,' KDIPOLE  --  Too many Bond Dipole',
     &                       ' Moment Parameters')
                  abort = .true.
               else if (iring .eq. 5) then
                  do j = 1, maxnd5
                     if (kd5(j).eq.'      ' .or. kd5(j).eq.pt) then
                        kd5(j) = pt
                        if (ia .le. ib) then
                           dpl5(j) = dp
                           pos5(j) = ps
                        else
                           dpl5(j) = -dp
                           pos5(j) = 1.0d0 - ps
                        end if
                        goto 110
                     end if
                  end do
                  write (iout,80)
   80             format (/,' KDIPOLE  --  Too many 5-Ring Bond',
     &                       ' Dipole Parameters')
                  abort = .true.
               else if (iring .eq. 4) then
                  do j = 1, maxnd4
                     if (kd4(j).eq.'      ' .or. kd4(j).eq.pt) then
                        kd4(j) = pt
                        if (ia .le. ib) then
                           dpl4(j) = dp
                           pos4(j) = ps
                        else
                           dpl4(j) = -dp
                           pos4(j) = 1.0d0 - ps
                        end if
                        goto 110
                     end if
                  end do
                  write (iout,90)
   90             format (/,' KDIPOLE  --  Too many 4-Ring Bond',
     &                       ' Dipole Parameters')
                  abort = .true.
               else if (iring .eq. 3) then
                  do j = 1, maxnd3
                     if (kd3(j).eq.'      ' .or. kd3(j).eq.pt) then
                        kd3(j) = pt
                        if (ia .le. ib) then
                           dpl3(j) = dp
                           pos3(j) = ps
                        else
                           dpl3(j) = -dp
                           pos3(j) = 1.0d0 - ps
                        end if
                        goto 110
                     end if
                  end do
                  write (iout,100)
  100             format (/,' KDIPOLE  --  Too many 3-Ring Bond',
     &                       ' Dipole Parameters')
                  abort = .true.
               end if
            end if
  110       continue
         end if
      end do
c
c     use small rings if present and parameters are available
c
      if (nring3.eq.0 .or. kd3(1).eq.'      ') then
         use_ring3 = .false.
      else
         use_ring3 = .true.
      end if
      if (nring4.eq.0 .or. kd4(1).eq.'      ') then
         use_ring4 = .false.
      else
         use_ring4 = .true.
      end if
      if (nring5.eq.0 .or. kd5(1).eq.'      ') then
         use_ring5 = .false.
      else
         use_ring5 = .true.
      end if
c
c     find and store all the bond dipole moments
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = type(ia)
         itb = type(ib)
         size = 3
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt = pa//pb
         else
            pt = pb//pa
         end if
         bdpl(i) = 0.0d0
c
c     make a check for bonds inside small rings
c
         iring = 0
         if (n12(ia).eq.1 .or. n12(ib).eq.1)  goto 120
         if (use_ring3) then
            do j = 1, n12(ia)
               ic = i12(j,ia)
               do k = 1, n12(ib)
                  if (ic .eq. i12(k,ib)) then
                     iring = 3
                     goto 120
                  end if
               end do
            end do
         end if
         if (use_ring4) then
            do j = 1, n12(ia)
               ic = i12(j,ia)
               if (ic .ne. ib) then
                  do k = 1, n12(ib)
                     id = i12(k,ib)
                     if (id .ne. ia) then
                        do m = 1, n12(ic)
                           if (id .eq. i12(m,ic)) then
                              iring = 4
                              goto 120
                           end if
                        end do
                     end if
                  end do
               end if
            end do
         end if
         if (use_ring5) then
            do j = 1, n13(ia)
               ic = i13(j,ia)
               do k = 1, n13(ib)
                  if (ic .eq. i13(k,ib)) then
                     iring = 5
                     goto 120
                  end if
               end do
            end do
         end if
  120    continue
c
c     try to assign bond dipole parameters for the bond
c
         if (iring .eq. 0) then
            do j = 1, maxnd
               if (kd(j) .eq. '      ')  goto 130
               if (kd(j) .eq. pt) then
                  if (ita .le. itb) then
                     idpl(1,i) = ia
                     idpl(2,i) = ib
                  else
                     idpl(1,i) = ib
                     idpl(2,i) = ia
                  end if
                  bdpl(i) = dpl(j)
                  sdpl(i) = pos(j)
               end if
            end do
  130       continue
         else if (iring .eq. 5) then
            do j = 1, maxnd5
               if (kd5(j) .eq. '      ')  goto 140
               if (kd5(j) .eq. pt) then
                  if (ita .le. itb) then
                     idpl(1,i) = ia
                     idpl(2,i) = ib
                  else
                     idpl(1,i) = ib
                     idpl(2,i) = ia
                  end if
                  bdpl(i) = dpl5(j)
                  sdpl(i) = pos5(j)
               end if
            end do
  140       continue
         else if (iring .eq. 4) then
            do j = 1, maxnd4
               if (kd4(j) .eq. '      ')  goto 150
               if (kd4(j) .eq. pt) then
                  if (ita .le. itb) then
                     idpl(1,i) = ia
                     idpl(2,i) = ib
                  else
                     idpl(1,i) = ib
                     idpl(2,i) = ia
                  end if
                  bdpl(i) = dpl4(j)
                  sdpl(i) = pos4(j)
               end if
            end do
  150       continue
         else if (iring .eq. 3) then
            do j = 1, maxnd3
               if (kd3(j) .eq. '      ')  goto 160
               if (kd3(j) .eq. pt) then
                  if (ita .le. itb) then
                     idpl(1,i) = ia
                     idpl(2,i) = ib
                  else
                     idpl(1,i) = ib
                     idpl(2,i) = ia
                  end if
                  bdpl(i) = dpl3(j)
                  sdpl(i) = pos3(j)
               end if
            end do
  160       continue
         end if
      end do
c
c     process keywords containing bond specific bond dipoles
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'DIPOLE ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.0d0
            string = record(next:80)
            read (string,*,err=170,end=170)  ia,ib,dp,ps
  170       continue
            if (ia.lt.0 .and. ib.lt.0) then
               ia = -ia
               ib = -ib
               if (header) then
                  header = .false.
                  write (iout,180)
  180             format (/,' Additional Bond Dipoles for',
     &                       ' Specific Bonds :',
     &                    //,5x,'Bonded Atoms',7x,'Moment',
     &                          6x,'Position',/)
               end if
               do j = 1, n12(ia)
                  if (i12(j,ia) .eq. ib) then
                     k = bndlist(j,ia)
                     if (ps .eq. 0.0d0)  ps = 0.5d0
                     if (idpl(1,k) .eq. ib) then
                        bdpl(k) = dp
                        sdpl(k) = ps
                     else
                        bdpl(k) = -dp
                        sdpl(k) = 1.0d0 - ps
                     end if
                     write (iout,190)  ia,ib,dp,ps
  190                format (4x,i5,' -',i5,2x,2f12.3)
                     goto 200
                  end if
               end do
            end if
  200       continue
         end if
      end do
c
c     remove zero bond dipoles from the list of dipoles
c
      ndipole = 0
      do i = 1, nbond
         if (bdpl(i) .ne. 0.0d0) then
            ndipole = ndipole + 1
            idpl(1,ndipole) = idpl(1,i)
            idpl(2,ndipole) = idpl(2,i)
            bdpl(ndipole) = bdpl(i)
            sdpl(ndipole) = sdpl(i)
         end if
      end do
c
c     if there are no dipole moments, turn off the potential
c
      if (ndipole .eq. 0) then
         use_dipole = .false.
         use_chgdpl = .false.
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine kimprop  --  improper dihedral parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "kimprop" assigns potential parameters to each improper
c     dihedral in the structure and processes any changed values
c
c
      subroutine kimprop
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'improp.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kiprop.i'
      include 'potent.i'
      integer i,j,k,size,next
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      real*8 tk,tv
      character*3 pa,pb,pc,pd
      character*12 blank,pti,pt(6)
      character*20 keyword
      character*80 record,string
      logical header
      data blank  / '            ' /
c
c
c     process keywords containing improper dihedral parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'IMPROPER ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            tk = 0.0d0
            tv = 0.0d0
            string = record(next:80)
            read (string,*,err=10,end=10)  ia,ib,ic,id,tk,tv
   10       continue
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            pti = pa//pb//pc//pd
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Improper Dihedral Parameters :',
     &                 //,5x,'Atom Classes',13x,'KID',11x,'Angle',/)
            end if
            write (iout,30)  ia,ib,ic,id,tk,tv
   30       format (1x,4i4,4x,2f14.4)
            do j = 1, maxndi
               if (kdi(j).eq.blank .or. kdi(j).eq.pti) then
                  kdi(j) = pti
                  dcon(j) = tk
                  tdi(j) = tv
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KIMPROP  --  Too many Improper Dihedral',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     assign improper dihedral parameters for each improper
c     angle by putting parameter values into "iprop" arrays
c
      niprop = 0
      do i = 1, n
         if (n12(i) .eq. 3) then
            ia = i
            ib = i12(1,i)
            ic = i12(2,i)
            id = i12(3,i)
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            do j = 1, maxndi
               if (kdi(j)(1:3) .eq. pa) then
                  pt(1) = pa//pb//pc//pd
                  pt(2) = pa//pb//pd//pc
                  pt(3) = pa//pc//pb//pd
                  pt(4) = pa//pc//pd//pb
                  pt(5) = pa//pd//pb//pc
                  pt(6) = pa//pd//pc//pb
                  do k = 1, 6
                     if (kdi(j) .eq. pt(k)) then
                        niprop = niprop + 1
                        iiprop(1,niprop) = ia
                        if (k.eq.1 .or. k.eq.2)  iiprop(2,niprop) = ib
                        if (k.eq.3 .or. k.eq.4)  iiprop(2,niprop) = ic
                        if (k.eq.5 .or. k.eq.6)  iiprop(2,niprop) = id
                        if (k.eq.3 .or. k.eq.5)  iiprop(3,niprop) = ib
                        if (k.eq.1 .or. k.eq.6)  iiprop(3,niprop) = ic
                        if (k.eq.2 .or. k.eq.4)  iiprop(3,niprop) = id
                        if (k.eq.4 .or. k.eq.6)  iiprop(4,niprop) = ib
                        if (k.eq.2 .or. k.eq.5)  iiprop(4,niprop) = ic
                        if (k.eq.1 .or. k.eq.3)  iiprop(4,niprop) = id
                        kprop(niprop) = dcon(j)
                        vprop(niprop) = tdi(j)
                        goto 60
                     end if
                  end do
               end if
            end do
   60       continue
         end if
      end do
c
c     if no improper dihedrals are used, turn off the potential
c
      if (niprop .eq. 0)  use_improp = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kimptor  --  improper torsional parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kimptor" assigns torsional parameters to each improper
c     torsion in the structure and processes any changed values
c
c
      subroutine kimptor
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'imptor.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kitors.i'
      include 'math.i'
      include 'potent.i'
      integer i,j,k,size,next
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer ft(6)
      real*8 angle,vt(6),st(6)
      character*3 pa,pb,pc,pd
      character*12 blank,pti,pt(6)
      character*20 keyword
      character*80 record,string
      logical header
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      data blank  / '            ' /
c
c
c     process keywords containing improper torsion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'IMPTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do j = 1, 6
               vt(j) = 0.0d0
               st(j) = 0.0d0
               ft(j) = 0
            end do
            string = record(next:80)
            read (string,*,err=10,end=10)  ia,ib,ic,id,
     &                                     (vt(j),st(j),ft(j),j=1,6)
   10       continue
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            pti = pa//pb//pc//pd
            call torphase (ft,vt,st)
            if (header) then
               header = .false.
               if (maswrk) write (iout,20)
   20          format (/,' Additional Improper Torsion Parameters :',
     &                 //,5x,'Atom Classes',12x,'1-Fold',12x,'2-Fold',
     &                    12x,'3-Fold',/)
            end if
            if (maswrk) write (iout,30) ia,ib,ic,id,(vt(j),st(j),j=1,3)
   30       format (1x,4i4,2x,3(f11.4,f7.1))
            do j = 1, maxnti
               if (kti(j).eq.blank .or. kti(j).eq.pti) then
                  kti(j) = pti
                  ti1(1,j) = vt(1)
                  ti1(2,j) = st(1)
                  ti2(1,j) = vt(2)
                  ti2(2,j) = st(2)
                  ti3(1,j) = vt(3)
                  ti3(2,j) = st(3)
                  goto 50
               end if
            end do
            if (maswrk) write (iout,40)
   40       format (/,' KIMPTOR  --  Too many Improper Torsion',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     assign improper torsional parameters for each improper
c     angle by putting parameter values into "itors" arrays
c
      nitors = 0
      do i = 1, n
         if (n12(i) .eq. 3) then
            ia = i12(1,i)
            ib = i12(2,i)
            ic = i
            id = i12(3,i)
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
            size = 3
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            do j = 1, maxnti
               if (kti(j)(7:9) .eq. pc) then
                  pt(1) = pa//pb//pc//pd
                  pt(2) = pb//pa//pc//pd
                  pt(3) = pa//pd//pc//pb
                  pt(4) = pd//pa//pc//pb
                  pt(5) = pb//pd//pc//pa
                  pt(6) = pd//pb//pc//pa
                  do k = 1, 6
                     if (kti(j) .eq. pt(k)) then
                        nitors = nitors + 1
                        if (k.eq.1 .or. k.eq.3)  iitors(1,nitors) = ia
                        if (k.eq.2 .or. k.eq.5)  iitors(1,nitors) = ib
                        if (k.eq.4 .or. k.eq.6)  iitors(1,nitors) = id
                        if (k.eq.2 .or. k.eq.4)  iitors(2,nitors) = ia
                        if (k.eq.1 .or. k.eq.6)  iitors(2,nitors) = ib
                        if (k.eq.3 .or. k.eq.5)  iitors(2,nitors) = id
                        iitors(3,nitors) = ic
                        if (k.eq.5 .or. k.eq.6)  iitors(4,nitors) = ia
                        if (k.eq.3 .or. k.eq.4)  iitors(4,nitors) = ib
                        if (k.eq.1 .or. k.eq.2)  iitors(4,nitors) = id
                        itors1(1,nitors) = ti1(1,j)
                        itors1(2,nitors) = ti1(2,j)
                        itors2(1,nitors) = ti2(1,j)
                        itors2(2,nitors) = ti2(2,j)
                        itors3(1,nitors) = ti3(1,j)
                        itors3(2,nitors) = ti3(2,j)
                        goto 60
                     end if
                  end do
               end if
            end do
   60       continue
         end if
      end do
c
c     find the cosine and sine of the phase angle for each torsion
c
      do i = 1, nitors
         angle = itors1(2,i) / radian
         itors1(3,i) = cos(angle)
         itors1(4,i) = sin(angle)
         angle = itors2(2,i) / radian
         itors2(3,i) = cos(angle)
         itors2(4,i) = sin(angle)
         angle = itors3(2,i) / radian
         itors3(3,i) = cos(angle)
         itors3(4,i) = sin(angle)
      end do
c
c     if no improper torsions are used, turn off the potential
c
      if (nitors .eq. 0)  use_imptor = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kmpole  --  multipole parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kmpole" assigns atomic multipole moments to the atoms of
c     the structure and processes any new or changed values
c
c
      subroutine kmpole
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kmulti.i'
      include 'mpole.i'
      include 'potent.i'
      include 'units.i'
      integer i,j,k,m,it,mt,imp,nmp
      integer size,next,number
      integer kz,kx,ztyp,xtyp
      integer mpt(maxnmp),mpkey(maxnmp)
      integer mpz(maxnmp),mpx(maxnmp)
      integer start(maxtyp),stop(maxtyp)
      real*8 mpl(13)
      character*3 pa,pb,pc
      character*8 axt
      character*9 pt
      character*20 keyword
      character*80 record,string
      logical header
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing atomic multipole parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            k = 0
            kz = 0
            kx = 0
            axt = 'Z-then-X'
            do j = 1, 13
               mpl(j) = 0.0d0
            end do
            string = record(next:80)
            read (string,*,err=50,end=50)  k,kz,kx,mpl(1)
            if (k .gt. 0) then
               if (kx .lt. 0)  axt = 'Bisector'
               kz = abs(kz)
               kx = abs(kx)
               record = keyline(i+1)
               read (record,*,err=50,end=50)  (mpl(j),j=2,4)
               record = keyline(i+2)
               read (record,*,err=50,end=50)  mpl(5)
               record = keyline(i+3)
               read (record,*,err=50,end=50)  (mpl(j),j=8,9)
               record = keyline(i+4)
               read (record,*,err=50,end=50)  (mpl(j),j=11,13)
               mpl(6) = mpl(8)
               mpl(7) = mpl(11)
               mpl(10) = mpl(12)
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,10)
   10             format (/,' Additional Atomic Multipole Parameters :',
     &                    //,6x,'Atom Type',5x,'Local Axes Definition',
     &                       10x,'Multipole Moments',/)
               end if
               if (maswrk) write (iout,20) k,kz,kx,axt,(mpl(j),j=1,5),
     &                             mpl(8),mpl(9),(mpl(j),j=11,13)
   20          format (6x,i6,5x,i6,1x,i6,3x,a8,4x,f9.5,/,45x,3f9.5,
     &                       /,45x,f9.5,/,45x,2f9.5,/,45x,3f9.5)
               size = 3
               call numeral (k,pa,size)
               call numeral (kz,pb,size)
               call numeral (kx,pc,size)
               pt = pa//pb//pc
               do j = 1, maxnmp
                  if (kmp(j).eq.'         ' .or. kmp(j).eq.pt) then
                     kmp(j) = pt
                     mpaxis(j) = axt
                     do m = 1, 13
                        multip(m,j) = mpl(m)
                     end do
                     goto 40
                  end if
               end do
               if (maswrk) write (iout,30)
   30          format (/,' KMPOLE  --  Too many Atomic Multipole',
     &                    ' Parameters')
               abort = .true.
   40          continue
            end if
   50       continue
         end if
      end do
c
c     zero out local axis definitions and multipoles for each atom
c
      do i = 1, n
         zaxis(i) = 0
         xaxis(i) = 0
         polaxe(i) = '        '
         do j = 1, 13
            pole(j,i) = 0.0d0
         end do
      end do
c
c     set indices into the list of atomic multipoles
c
      nmp = 0
      dowhile (kmp(nmp+1) .ne. '         ')
         nmp = nmp + 1
         mpt(nmp) = number(kmp(nmp)(1:3))
         mpz(nmp) = number(kmp(nmp)(4:6))
         mpx(nmp) = number(kmp(nmp)(7:9))
      end do
      call sort3 (nmp,mpt,mpkey)
      do i = 1, maxtyp
         start(i) = 0
         stop(i) = 0
      end do
      do i = 1, nmp
         k = mpt(i)
         if (start(k) .eq. 0)  start(k) = i
      end do
      do i = nmp, 1, -1
         k = mpt(i)
         if (stop(k) .eq. 0)  stop(k) = i
      end do
c
c     assign atomic multipole parameters to each atom
c
      do i = 1, n
         it = type(i)
         if (start(it) .ne. 0) then
            do j = start(it), stop(it)
               imp = mpkey(j)
               ztyp = mpz(imp)
               xtyp = mpx(imp)
               if (n12(i) .eq. 0) then
                  zaxis(i) = 1
                  xaxis(i) = 2
                  if (i .eq. 1)  zaxis(i) = 3
                  if (i .eq. 2)  xaxis(i) = 3
                  goto 60
               end if
               do k = 1, n12(i)
                  if (type(i12(k,i)) .eq. ztyp) then
                     if (n12(i) .ge. 2) then
                        do m = 1, n12(i)
                           mt = type(i12(m,i))
                           if (mt.eq.xtyp .and. m.ne.k) then
                              zaxis(i) = i12(k,i)
                              xaxis(i) = i12(m,i)
                              goto 60
                           end if
                        end do
                     else if (n12(i) .eq. 1) then
                        do m = 1, n13(i)
                           mt = type(i13(m,i))
                           if (mt .eq. xtyp) then
                              zaxis(i) = i12(k,i)
                              xaxis(i) = i13(m,i)
                              goto 60
                           end if
                        end do
                     end if
                  end if
               end do
            end do
   60       continue
            if (zaxis(i) .ne. 0) then
               polaxe(i) = mpaxis(imp)
               do j = 1, 13
                  pole(j,i) = multip(j,imp)
               end do
            end if
         end if
      end do
c
c     process keywords with multipole parameters for specific atoms
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            k = 0
            kz = 0
            kx = 0
            axt = 'Z-then-X'
            do j = 1, 13
               mpl(j) = 0.0d0
            end do
            string = record(next:80)
            read (string,*,err=90,end=90)  k,kz,kx,mpl(1)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               if (kx .lt. 0)  axt = 'Bisector'
               kz = abs(kz)
               kx = abs(kx)
               record = keyline(i+1)
               read (record,*,err=90,end=90)  (mpl(j),j=2,4)
               record = keyline(i+2)
               read (record,*,err=90,end=90)  mpl(5)
               record = keyline(i+3)
               read (record,*,err=90,end=90)  (mpl(j),j=8,9)
               record = keyline(i+4)
               read (record,*,err=90,end=90)  (mpl(j),j=11,13)
               mpl(6) = mpl(8)
               mpl(7) = mpl(11)
               mpl(10) = mpl(12)
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,70)
   70             format (/,' Additional Atomic Multipoles',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',8x,'Local Axes Definition',
     &                       10x,'Multipole Moments',/)
               end if
               if (maswrk) write (iout,80)  k,kz,kx,axt,(mpl(j),j=1,5),
     &                              mpl(8),mpl(9),(mpl(j),j=11,13)
   80          format (6x,i6,5x,i6,1x,i6,3x,a8,4x,f9.5,/,45x,3f9.5,
     &                       /,45x,f9.5,/,45x,2f9.5,/,45x,3f9.5)
               zaxis(k) = kz
               xaxis(k) = kx
               polaxe(k) = axt
               do j = 1, 13
                  pole(j,k) = mpl(j)
               end do
            end if
   90       continue
         end if
      end do
c
c     convert the dipole and quadrupole moments to Angstroms,
c     quadrupole divided by 3 for Applequist-Dykstra method
c
      do i = 1, n
         do k = 2, 4
            pole(k,i) = pole(k,i) * bohr
         end do
         do k = 5, 13
            pole(k,i) = pole(k,i) * bohr**2 / 3.0d0
         end do
      end do
c
c     get the order of the multipole expansion at each site
c
      do i = 1, n
         size = 0
         do k = 1, 13
            if (pole(k,i) .ne. 0.0d0)  size = max(k,size)
         end do
         if (size .gt. 4) then
            size = 13
         else if (size .gt. 1) then
            size = 4
         end if
         polsiz(i) = size
      end do
c
c     if polarization not used, remove zero and undefined multipoles
c
      if (.not. use_polar) then
         npole = 0
         do i = 1, n
            if (polsiz(i) .ne. 0) then
               if (zaxis(i).ne.0 .and. xaxis(i).ne.0) then
                  npole = npole + 1
                  ipole(npole) = i
                  zaxis(npole) = zaxis(i)
                  xaxis(npole) = xaxis(i)
                  polaxe(npole) = polaxe(i)
                  polsiz(npole) = polsiz(i)
                  do j = 1, 13
                     pole(j,npole) = pole(j,i)
                  end do
               end if
            end if
         end do
c
c     if no atomic multipoles are used, turn off the potential
c
         if (npole .eq. 0)  use_mpole = .false.
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kopbend  --  out-of-plane bending parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kopbend" assigns the force constants for out-of-plane bends;
c     also processes any new or changed parameter values
c
c
      subroutine kopbend
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kopbnd.i'
      include 'opbend.i'
      include 'potent.i'
      integer i,j,next,size,number
      integer ia,ib,id,it,itb,itd
      real*8 fopb
      character*3 pa,pb,pd
      character*6 ps
      character*20 keyword
      character*80 record,string
      logical header,jopb(maxtyp)
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing out-of-plane bend parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            fopb = 0.0d0
            string = record(next:80)
            read (string,*,err=10,end=10)  ia,ib,fopb
   10       continue
            if (header) then
               header = .false.
               if (maswrk) write (iout,20)
   20          format (/,' Additional Out-of-Plane Bend Parameters :',
     &                 //,5x,'Atom Classes',7x,'K(OPB)',/)
            end if
            if (maswrk) write (iout,30)  ia,ib,fopb
   30       format (6x,2i4,4x,f12.3)
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            ps = pa//pb
            do j = 1, maxnopb
               if (kaopb(j).eq.'      ' .or. kaopb(j).eq.ps) then
                  kaopb(j) = ps
                  copb(j) = fopb
                  goto 50
               end if
            end do
            if (maswrk) write (iout,40)
   40       format (/,' KOPBEND --  Too many Out-of-Plane',
     &                 ' Angle Bending Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     make list of atom types using out-of-plane bending
c
      do i = 1, maxtyp
         jopb(i) = .false.
      end do
      do i = 1, maxnopb
         if (kaopb(i) .eq. '      ')  goto 60
         it = number (kaopb(i)(1:3))
         jopb(it) = .true.
      end do
   60 continue
c
c     assign out-of-plane bending parameters for each angle
c
      header = .true.
      nopbend = 0
      do i = 1, nangle
         angin(i) = .false.
         ib = iang(2,i)
         itb = class(ib)
         if (jopb(itb) .and. n12(ib).eq.3) then
            id = iang(4,i)
            itd = class(id)
            size = 3
            call numeral (itb,pb,size)
            call numeral (itd,pd,size)
            ps = pb//pd
            do j = 1, maxnopb
               if (kaopb(j) .eq. ps) then
                  nopbend = nopbend + 1
                  iopb(nopbend) = i
                  kopb(nopbend) = copb(j)
                  angin(i) = .true.
                  goto 90
               end if
            end do
            abort = .true.
            if (header) then
               header = .false.
               if (maswrk) write (iout,70)
   70          format (/,' Undefined Out-of-Plane Bending Parameters :',
     &                 //,' Type',13x,'Atom Names',11x,'Atom Classes',/)
            end if
            if (maswrk) write (iout,80)  ib,name(ib),id,name(id),itb,itd
   80       format (' Angle-OP',3x,i6,'-',a3,i6,'-',a3,7x,2i5)
   90       continue
         else
            iang(4,i) = ib
         end if
      end do
c
c     if no out-of-plane bends are used, turn off the potential
c
      if (nopbend .eq. 0)  use_opbend = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine korbit  --  pisystem parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "korbit" assigns pi-orbital parameters to conjugated
c     systems and processes any new or changed parameters
c
c
      subroutine korbit
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'korbs.i'
      include 'picalc.i'
      include 'pistuf.i'
      include 'piterm.i'
      include 'tors.i'
      include 'units.i'
      integer i,j,k,size,next
      integer ia,ib,it,ita,itb
      real*8 elect,ioniz,repuls
      real*8 sslop,tslop
      character*3 pa,pb
      character*6 pt
      character*20 keyword
      character*80 record,string
      logical header,done
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing pisystem atom parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'PIATOM ') then
            if (header) then
               header = .false.
               if (maswrk) write (iout,10)
   10          format (/,' Additional Pisystem Atom Parameters :',
     &                 //,6x,'Atom Type',3x,'Electron',3x,'Ionization',
     &                    3x,'Repulsion',/)
            end if
            ia = 0
            elect = 0.0d0
            ioniz = 0.0d0
            repuls = 0.0d0
            string = record(next:80)
            read (string,*,err=20)  ia,elect,ioniz,repuls
   20       continue
            if (maswrk) write (iout,30)  ia,elect,ioniz,repuls
   30       format (8x,i4,3x,f10.3,3x,f10.3,2x,f10.3)
            if (ia.gt.0 .and. ia.le.maxclass) then
               electron(ia) = elect
               ionize(ia) = ioniz
               repulse(ia) = repuls
            else
   40          format (/,' KORBIT  --  Too many Atom Classes;',
     &                    ' Increase MAXCLASS')
               abort = .true.
            end if
         end if
      end do
c
c     process keywords containing pisystem bond parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'PIBOND ') then
            if (header) then
               header = .false.
               if (maswrk) write (iout,50)
   50          format (/,' Additional Pisystem Bond Parameters :',
     &                 //,6x,'Atom Types',7x,'d Force',4x,'d Length',/)
            end if
            ia = 0
            ib = 0
            sslop = 0.0d0
            tslop = 0.0d0
            string = record(next:80)
            read (string,*,err=60)  ia,ib,sslop,tslop
   60       continue
            if (maswrk) write (iout,70)  ia,ib,sslop,tslop
   70       format (6x,2i4,5x,2f11.3)
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do j = 1, maxnpi
               if (kpi(j).eq.'      ' .or. kpi(j).eq.pt) then
                  kpi(j) = pt
                  sslope(j) = sslop
                  tslope(j) = tslop
                  goto 90
               end if
            end do
            if (maswrk) write (iout,80)
   80       format (/,' KORBIT  --  Too many Pi-System Bond',
     &                 ' Type Parameters')
            abort = .true.
   90       continue
         end if
      end do
c
c     assign the values characteristic of the piatom types;
c     count the number of filled pi molecular orbitals
c
      nfill = 0
      do i = 1, norbit
         it = type(iorbit(i))
         q(i) = electron(it)
         w(i) = ionize(it) / evolt
         em(i) = repulse(it) / evolt
         nfill = nfill + nint(q(i))
      end do
      nfill = nfill / 2
c
c     assign parameters for all bonds between piatoms;
c     store the original bond lengths and force constants
c
      do i = 1, npibond
         j = pibond(1,i)
         ia = pibond(2,i)
         ib = pibond(3,i)
         ita = class(iorbit(ia))
         itb = class(iorbit(ib))
         size = 3
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt = pa//pb
         else
            pt = pb//pa
         end if
         done = .false.
         k = 0
         dowhile (.not.done .and. k.lt.maxnpi)
            k = k + 1
            if (kpi(k) .eq. pt) then
               done = .true.
               bkpi(i) = bk(j)
               blpi(i) = bl(j)
               kslope(i) = sslope(k)
               lslope(i) = tslope(k)
            end if
         end do
         if (.not.done) then
            if (maswrk) write (iout,100)  ita,itb
  100       format (/,' KORBIT  --  No Parameters for Pi-Bond',
     &                 ' between Atom Types',2i4)
         end if
      end do
c
c     store the original torsional constants across pibonds
c
      do i = 1, npitors
         j = pitors(1,i)
         torspi(i) = tors2(1,j)
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c
      subroutine kpolar
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kpolr.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      integer i,k,next
      real*8 pol,proot,sixth
      character*20 keyword
      character*80 record,string
      logical header
c
c
c     process keywords containing polarizability parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            string = record(next:80)
            read (string,*,err=40,end=40)  k,pol
            if (k .gt. 0) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :',
     &                    //,6x,'Atom Type',10x,'Alpha',/)
               end if
               if (k .le. maxtyp) then
                  polr(k) = pol
                  write (iout,20)  k,pol
   20             format (6x,i6,6x,f12.3)
               else
                  write (iout,30)
   30             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' Polarizability Parameters')
                  abort = .true.
               end if
            end if
   40       continue
         end if
      end do
c
c     find and store all the atomic dipole polarizabilities
c
      do i = 1, n
         polarize(i) = polr(type(i))
      end do
c
c     process keywords containing atom specific polarizabilities
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            string = record(next:80)
            read (string,*,err=70,end=70)  k,pol
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               if (header) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Dipole Polarizabilities',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',13x,'Alpha',/)
               end if
               write (iout,60)  k,pol
   60          format (6x,i6,6x,f12.3)
               polarize(k) = pol
            end if
   70       continue
         end if
      end do
c
c     remove zero and undefined polarizable multipoles from the list
c
      npole = 0
      do i = 1, n
         if (polsiz(i).ne.0 .or. polarize(i).ne.0.0d0) then
            if (zaxis(i).ne.0 .and. xaxis(i).ne.0) then
               npole = npole + 1
               ipole(npole) = i
               zaxis(npole) = zaxis(i)
               xaxis(npole) = xaxis(i)
               polaxe(npole) = polaxe(i)
               polsiz(npole) = polsiz(i)
               do k = 1, maxpole
                  pole(k,npole) = pole(k,i)
               end do
               polarize(npole) = polarize(i)
            end if
         end if
      end do
c
c     set the values used in the damping of the polarizability
c
      sixth = 1.0d0 / 6.0d0
      proot = 0.0d0
      if (pradius .gt. 0.0d0)  proot = sqrt(pradius)
      do i = 1, npole
         pdamp(i) = proot * polarize(i)**sixth
      end do
c
c     count the final number of polarizable multipole sites
c
      npolar = 0
      do i = 1, npole
         if (polarize(i) .ne. 0.0d0)  npolar = npolar + 1
      end do
c
c     if no polarizable multipoles are used, turn off the potentials
c
      if (npole .eq. 0)  use_mpole = .false.
      if (npolar .eq. 0)  use_polar = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kstrbnd  --  assign stretch-bend parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kstrbnd" assigns the parameters for the stretch-bend
c     interactions and processes new or changed parameter values
c
c
      subroutine kstrbnd
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kstbnd.i'
      include 'potent.i'
      include 'strbnd.i'
      integer i,j,k
      integer it,itb,next
      integer ia,ib,ic
      integer nh,nb1,nb2
      real*8 sbk(3)
      character*20 keyword
      character*80 record,string
      logical header
c
c
c     process keywords containing stretch-bend parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'STRBND ') then
            it = 0
            do j = 1, 3
               sbk(j) = 0.0d0
            end do
            string = record(next:80)
            read (string,*,err=10,end=10)  it,(sbk(j),j=1,3)
   10       continue
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Stretch-Bend Parameters :',
     &                 //,5x,'Atom Class',8x,'K(SB) 1',4x,'K(SB) 2',
     &                    4x,'K(SB) 3',/)
            end if
            write (iout,30)  it,(sbk(j),j=1,3)
   30       format (9x,i3,7x,3f11.3)
            do j = 1, 3
               stbn(j,it) = sbk(j)
            end do
         end if
      end do
c
c     assign stretch-bend parameters for each angle
c
      nstrbnd = 0
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         itb = class(ib)
         do k = 1, n12(ib)
            if (i12(k,ib) .eq. ia)  nb1 = bndlist(k,ib)
            if (i12(k,ib) .eq. ic)  nb2 = bndlist(k,ib)
         end do
         nh = 1
         if (atomic(ia) .le. 1)  nh = nh + 1
         if (atomic(ic) .le. 1)  nh = nh + 1
         if (stbn(nh,itb) .ne. 0.0d0) then
            nstrbnd = nstrbnd + 1
            isb(1,nstrbnd) = i
            isb(2,nstrbnd) = nb1
            isb(3,nstrbnd) = nb2
            ksb(nstrbnd) = stbn(nh,itb)
         end if
      end do
c
c     if no stretch-bends are used, turn off the potential
c
      if (nstrbnd .eq. 0)  use_strbnd = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kstrtor  --  find stretch-torsion parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kstrtor" assigns stretch-torsion parameters to torsions
c     needing them, and processes any new or changed values
c
c
      subroutine kstrtor
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'ksttor.i'
      include 'potent.i'
      include 'strtor.i'
      include 'tors.i'
      integer i,j,k,size,next
      integer ib,ic,itb,itc
      real*8 bt1,bt2,bt3
      character*3 pb,pc
      character*6 pt
      character*20 keyword
      character*80 record,string
      logical header
c
c
c     process keywords containing stretch-torsion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'STRTORS ') then
            ib = 0
            ic = 0
            bt1 = 0.0d0
            bt2 = 0.0d0
            bt3 = 0.0d0
            string = record(next:80)
            read (string,*,err=10,end=10)  ib,ic,bt1,bt2,bt3
   10       continue
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Stretch-Torsion Parameters :',
     &                 //,5x,'Atom Classes',7x,'1-Fold',6x,'2-Fold',
     &                    6x,'3-Fold',/)
            end if
            write (iout,30)  ib,ic,bt1,bt2,bt3
   30       format (6x,2i4,4x,3f12.3)
            size = 3
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ib .le. ic) then
               pt = pb//pc
            else
               pt = pc//pb
            end if
            do j = 1, maxnbt
               if (kbt(j).eq.'      ' .or. kbt(j).eq.pt) then
                  kbt(j) = pt
                  btcon(1,j) = bt1
                  btcon(2,j) = bt2
                  btcon(3,j) = bt3
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KSTRTOR  --  Too many Stretch-Torsion',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     assign stretch-torsion parameters as required
c
      nstrtor = 0
      do i = 1, ntors
         ib = itors(2,i)
         ic = itors(3,i)
         itb = class(ib)
         itc = class(ic)
         size = 3
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (itb .le. itc) then
            pt = pb//pc
         else
            pt = pc//pb
         end if
         do j = 1, maxnbt
            if (kbt(j) .eq. pt) then
               nstrtor = nstrtor + 1
               kst(1,nstrtor) = btcon(1,j)
               kst(2,nstrtor) = btcon(2,j)
               kst(3,nstrtor) = btcon(3,j)
               ist(1,nstrtor) = i
               do k = 1, n12(ib)
                  if (i12(k,ib) .eq. ic) then
                     ist(2,nstrtor) = bndlist(k,ib)
                     goto 60
                  end if
               end do
            end if
         end do
   60    continue
      end do
c
c     if no stretch-torsions are used, turn off the potential
c
      if (nstrtor .eq. 0)  use_strtor = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ktors  --  torsional parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ktors" assigns torsional parameters to each torsion in
c     the structure and processes any new or changed values
c
c
      subroutine ktors
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'ktorsn.i'
      include 'math.i'
      include 'potent.i'
      include 'ring.i'
      include 'tors.i'
      integer i,j,k,size,next
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer iring,nuse,iuse
      integer kindex(maxnt)
      integer ft(6)
      real*8 angle,vt(6),st(6)
      character*3 pa,pb,pc,pd
casa  -modi_dihed- begin
      character*3 zeros
c      character*12 blank,pt,kuse(maxnt)
      character*12 blank,pt,kuse(maxnt), pt0
casa  -modi_dihed- end
      character*20 keyword
      character*80 record,string
      logical header,dummy
      logical use_ring5,use_ring4
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      data blank  / '            ' /
c
c
c     process keywords containing torsional angle parameters
c
casa  -modi_dihed- 
      zeros = '000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:8) .eq. 'TORSION ')  iring = 0
         if (keyword(1:9) .eq. 'TORSION5 ')  iring = 5
         if (keyword(1:9) .eq. 'TORSION4 ')  iring = 4
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do j = 1, 6
               vt(j) = 0.0d0
               st(j) = 0.0d0
               ft(j) = 0
            end do
            string = record(next:80)
            read (string,*,err=10,end=10)  ia,ib,ic,id,
     &                                     (vt(j),st(j),ft(j),j=1,6)
   10       continue
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            if (ib .lt. ic) then
               pt = pa//pb//pc//pd
            else if (ic .lt. ib) then
               pt = pd//pc//pb//pa
            else if (ia .le. id) then
               pt = pa//pb//pc//pd
            else if (id .lt. ia) then
               pt = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            if (header) then
               header = .false.
               if (maswrk) write (iout,20)
   20          format (/,' Additional Torsional Parameters :',
     &                 //,5x,'Atom Classes',4x,'1-Fold',4x,'2-Fold',
     &                    4x,'3-Fold',4x,'4-Fold',4x,'5-Fold',
     &                    4x,'6-Fold',/)
            end if
            if (iring .eq. 0) then
               if (maswrk) write (iout,30)  ia,ib,ic,id,
     &                          (vt(j),nint(st(j)),j=1,6)
   30          format (1x,4i4,1x,6(f6.2,i4))
            else if (iring .eq. 5) then
               if (maswrk) write (iout,40)  ia,ib,ic,id,
     &                          (vt(j),nint(st(j)),j=1,6)
   40          format (1x,4i4,1x,6(f6.2,i4),3x,'5-Ring')
            else if (iring .eq. 4) then
               if (maswrk) write (iout,50)  ia,ib,ic,id,
     &                          (vt(j),nint(st(j)),j=1,6)
   50          format (1x,4i4,1x,6(f6.2,i4),3x,'4-Ring')
            end if
            if (iring .eq. 0) then
               do j = 1, maxnt
                  if (kt(j).eq.blank .or. kt(j).eq.pt) then
                     kt(j) = pt
                     t1(1,j) = vt(1)
                     t1(2,j) = st(1)
                     t2(1,j) = vt(2)
                     t2(2,j) = st(2)
                     t3(1,j) = vt(3)
                     t3(2,j) = st(3)
                     t4(1,j) = vt(4)
                     t4(2,j) = st(4)
                     t5(1,j) = vt(5)
                     t5(2,j) = st(5)
                     t6(1,j) = vt(6)
                     t6(2,j) = st(6)
                     goto 70
                  end if
               end do
               if (maswrk) write (iout,60)
   60          format (/,' KTORS   --  Too many Torsional Angle',
     &                    ' Parameters')
               abort = .true.
   70          continue
            else if (iring .eq. 5) then
               do j = 1, maxnt5
                  if (kt5(j).eq.blank .or. kt5(j).eq.pt) then
                     kt5(j) = pt
                     t15(1,j) = vt(1)
                     t15(2,j) = st(1)
                     t25(1,j) = vt(2)
                     t25(2,j) = st(2)
                     t35(1,j) = vt(3)
                     t35(2,j) = st(3)
                     t45(1,j) = vt(4)
                     t45(2,j) = st(4)
                     t55(1,j) = vt(5)
                     t55(2,j) = st(5)
                     t65(1,j) = vt(6)
                     t65(2,j) = st(6)
                     goto 90
                  end if
               end do
               if (maswrk) write (iout,80)
   80          format (/,' KTORS   --  Too many 5-Ring Torsional',
     &                    ' Parameters')
               abort = .true.
   90          continue
            else if (iring .eq. 4) then
               do j = 1, maxnt4
                  if (kt4(j).eq.blank .or. kt4(j).eq.pt) then
                     kt4(j) = pt
                     t14(1,j) = vt(1)
                     t14(2,j) = st(1)
                     t24(1,j) = vt(2)
                     t24(2,j) = st(2)
                     t34(1,j) = vt(3)
                     t34(2,j) = st(3)
                     t44(1,j) = vt(4)
                     t44(2,j) = st(4)
                     t54(1,j) = vt(5)
                     t54(2,j) = st(5)
                     t64(1,j) = vt(6)
                     t64(2,j) = st(6)
                     goto 110
                  end if
               end do
               if (maswrk) write (iout,100)
  100          format (/,' KTORS   --  Too many 4-Ring Torsional',
     &                    ' Parameters')
               abort = .true.
  110          continue
            end if
         end if
      end do
c
c     use small rings if present and parameters are available
c
      if (nring4.eq.0 .or. kt4(1).eq.blank) then
         use_ring4 = .false.
      else
         use_ring4 = .true.
      end if
      if (nring5.eq.0 .or. kt5(1).eq.blank) then
         use_ring5 = .false.
      else
         use_ring5 = .true.
      end if
c
c     assign torsional parameters for each torsional angle
c     by putting the parameter values into the "tors" arrays
c
      header = .true.
      nuse = 0
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         dummy = .false.
         if (min(atomic(ia),atomic(ib),
     &           atomic(ic),atomic(id)) .eq. 0) then
            dummy = .true.
         end if
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         size = 3
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         call numeral (itd,pd,size)
         if (itb .lt. itc) then
            pt = pa//pb//pc//pd
         else if (itc .lt. itb) then
            pt = pd//pc//pb//pa
         else if (ita .le. itd) then
            pt = pa//pb//pc//pd
         else if (itd .lt. ita) then
            pt = pd//pc//pb//pa
         end if
casa -modi_dihed- 
         pt0 = zeros//pt(4:9)//zeros
         tors1(1,i) = 0.0d0
         tors1(2,i) = 0.0d0
         tors2(1,i) = 0.0d0
         tors2(2,i) = 0.0d0
         tors3(1,i) = 0.0d0
         tors3(2,i) = 0.0d0
         tors4(1,i) = 0.0d0
         tors4(2,i) = 0.0d0
         tors5(1,i) = 0.0d0
         tors5(2,i) = 0.0d0
         tors6(1,i) = 0.0d0
         tors6(2,i) = 0.0d0
c
c     make a check for torsions inside small rings
c
         iring = 0
         if (n12(ia).eq.1 .or. n12(id).eq.1)  goto 120
         if (use_ring4) then
            do j = 1, n12(id)
               if (ia .eq. i12(j,id)) then
                  iring = 4
                  goto 120
               end if
            end do
         end if
         if (use_ring5) then
            do j = 1, n12(id)
               if (ic .ne. i12(j,id)) then
                  do k = 1, n12(ia)
                     if (i12(j,id) .eq. i12(k,ia)) then
                        iring = 5
                        goto 120
                     end if
                  end do
               end if
            end do
         end if
  120    continue
c
c     find parameters for this torsion; first check "kuse"
c     to save time for angle types already located
c
         if (iring .eq. 0) then
            do j = 1, nuse
               if (kuse(j) .eq. pt) then
                  iuse = kindex(j)
                  tors1(1,i) = tors1(1,iuse)
                  tors1(2,i) = tors1(2,iuse)
                  tors2(1,i) = tors2(1,iuse)
                  tors2(2,i) = tors2(2,iuse)
                  tors3(1,i) = tors3(1,iuse)
                  tors3(2,i) = tors3(2,iuse)
                  tors4(1,i) = tors4(1,iuse)
                  tors4(2,i) = tors4(2,iuse)
                  tors5(1,i) = tors5(1,iuse)
                  tors5(2,i) = tors5(2,iuse)
                  tors6(1,i) = tors6(1,iuse)
                  tors6(2,i) = tors6(2,iuse)
                  goto 190
               end if
            end do
            do j = 1, maxnt
               if (kt(j) .eq. pt) then
                  nuse = nuse + 1
                  kuse(nuse) = pt
                  kindex(nuse) = i
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  tors4(1,i) = t4(1,j)
                  tors4(2,i) = t4(2,j)
                  tors5(1,i) = t5(1,j)
                  tors5(2,i) = t5(2,j)
                  tors6(1,i) = t6(1,j)
                  tors6(2,i) = t6(2,j)
                  goto 190
               end if
            end do
casa -modi_dihed- begin
            do j = 1, maxnt
               if (kt(j) .eq. pt0) then
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  tors4(1,i) = t4(1,j)
                  tors4(2,i) = t4(2,j)
                  tors5(1,i) = t5(1,j)
                  tors5(2,i) = t5(2,j)
                  tors6(1,i) = t6(1,j)
                  tors6(2,i) = t6(2,j)
                  goto 190
               end if
            end do
casa -modi_dihed-  end
            if (.not. dummy) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,130)
  130             format (/,' Undefined Torsional Parameters :',
     &                    //,' Type',24x,'Atom Names',24x,
     &                       'Atom Classes',/)
               end if
               if (maswrk) write (iout,140) ia,name(ia),ib,name(ib),
     &                      ic,name(ic),id,name(id),ita,itb,itc,itd
  140          format (' Torsion',4x,4(i6,'-',a3),5x,4i5)
            end if
c
c     find the parameters for a 5-ring torsion
c
         else if (iring .eq. 5) then
            do j = 1, maxnt5
               if (kt5(j) .eq. pt) then
                  tors1(1,i) = t15(1,j)
                  tors1(2,i) = t15(2,j)
                  tors2(1,i) = t25(1,j)
                  tors2(2,i) = t25(2,j)
                  tors3(1,i) = t35(1,j)
                  tors3(2,i) = t35(2,j)
                  tors4(1,i) = t45(1,j)
                  tors4(2,i) = t45(2,j)
                  tors5(1,i) = t55(1,j)
                  tors5(2,i) = t55(2,j)
                  tors6(1,i) = t65(1,j)
                  tors6(2,i) = t65(2,j)
                  goto 190
               end if
            end do
            if (.not. dummy) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,150)
  150             format (/,' Undefined Torsional Parameters :',
     &                    //,' Type',24x,'Atom Names',24x,
     &                       'Atom Classes',/)
               end if
      if (maswrk) write (iout,160)  ia,name(ia),ib,name(ib),ic,name(ic),
     &                           id,name(id),ita,itb,itc,itd
  160          format (' 5-Ring',5x,4(i6,'-',a3),5x,4i5)
            end if
c
c     find the parameters for a 4-ring torsion
c
         else if (iring .eq. 4) then
            do j = 1, maxnt4
               if (kt4(j) .eq. pt) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  tors4(1,i) = t44(1,j)
                  tors4(2,i) = t44(2,j)
                  tors5(1,i) = t54(1,j)
                  tors5(2,i) = t54(2,j)
                  tors6(1,i) = t64(1,j)
                  tors6(2,i) = t64(2,j)
                  goto 190
               end if
            end do
            if (.not. dummy) then
               abort = .true.
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,170)
  170             format (/,' Undefined Torsional Parameters :',
     &                    //,' Type',24x,'Atom Names',24x,
     &                       'Atom Classes',/)
               end if
      if (maswrk) write (iout,180)  ia,name(ia),ib,name(ib),ic,name(ic),
     &                           id,name(id),ita,itb,itc,itd
  180          format (' 4-Ring',5x,4(i6,'-',a3),5x,4i5)
            end if
         end if
  190    continue
      end do
c
c     find the cosine and sine of the phase angle for each torsion
c
      do i = 1, ntors
         angle = tors1(2,i) / radian
         tors1(3,i) = cos(angle)
         tors1(4,i) = sin(angle)
         angle = tors2(2,i) / radian
         tors2(3,i) = cos(angle)
         tors2(4,i) = sin(angle)
         angle = tors3(2,i) / radian
         tors3(3,i) = cos(angle)
         tors3(4,i) = sin(angle)
         angle = tors4(2,i) / radian
         tors4(3,i) = cos(angle)
         tors4(4,i) = sin(angle)
         angle = tors5(2,i) / radian
         tors5(3,i) = cos(angle)
         tors5(4,i) = sin(angle)
         angle = tors6(2,i) / radian
         tors6(3,i) = cos(angle)
         tors6(4,i) = sin(angle)
      end do
c
c     if no torsion terms are used, turn off the potential
c
      if (ntors .eq. 0)  use_tors = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kurey  --  Urey-Bradley parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kurey" assigns the force constants and ideal distances
c     for the Urey-Bradley 1-3 interactions; also processes any
c     new or changed parameter values
c
c
      subroutine kurey
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kurybr.i'
      include 'potent.i'
      include 'urey.i'
      integer i,j,size,next
      integer ia,ib,ic
      integer ita,itb,itc
      real*8 bb,tt
      character*3 pa,pb,pc
      character*9 pt
      character*20 keyword
      character*80 record,string
      logical header
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing Urey-Bradley parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'UREYBRAD ') then
            ia = 0
            ib = 0
            ic = 0
            bb = 0.0d0
            tt = 0.0d0
            string = record(next:80)
            read (string,*,err=10,end=10)  ia,ib,ic,bb,tt
   10       continue
            if (header) then
               header = .false.
               if (maswrk) write (iout,20)
   20          format (/,' Additional Urey-Bradley Parameters :',
     &                 //,5x,'Atom Classes',8x,'K(UB)',6x,'Distance',/)
            end if
            if (maswrk) write (iout,30)  ia,ib,ic,bb,tt
   30       format (4x,3i4,2x,2f12.3)
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            do j = 1, maxnu
               if (ku(j).eq.'         ' .or. ku(j).eq.pt) then
                  ku(j) = pt
                  ucon(j) = bb
                  dst13(j) = tt
                  goto 50
               end if
            end do
            if (maswrk) write (iout,40)
   40       format (/,' KUREY  --  Too many Urey-Bradley',
     &                 ' Interaction Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     loop over all angles assigning the necessary constants
c
      nurey = 0
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 3
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pt = pa//pb//pc
         else
            pt = pc//pb//pa
         end if
c
c     search for Urey-Bradley parameters for current angle
c
         do j = 1, maxnu
            if (ku(j) .eq. pt) then
               nurey = nurey + 1
               iury(1,nurey) = ia
               iury(2,nurey) = ic
               uk(nurey) = ucon(j)
               ul(nurey) = dst13(j)
               goto 60
            end if
         end do
   60    continue
      end do
c
c     if no Urey-Bradley terms are used, turn off the potential
c
      if (nurey .eq. 0)  use_urey = .false.
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kvdw  --  van der Waals parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kvdw" assigns the parameters to be used in computing the
c     van der Waals interactions and processes any new or changed
c     values for these parameters
c
c
      subroutine kvdw
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'khbond.i'
      include 'kvdws.i'
      include 'kvdwpr.i'
      include 'math.i'
      include 'potent.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,k,ia,ib,next
      integer size,number
      real*8 rd,ep,rdn
      real*8 srad(maxclass)
      real*8 seps(maxclass)
      character*3 pa,pb
      character*6 pt
      character*20 keyword
      character*80 record,string
      logical header
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     process keywords containing van der Waals parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:4) .eq. 'VDW ') then
            k = 0
            rd = 0.0d0
            ep = 0.0d0
            rdn = 0.0d0
            string = record(next:80)
            read (string,*,err=10,end=10)  k,rd,ep,rdn
   10       continue
            if (k.ge.1 .and. k.le.maxtyp) then
               if (header) then
                  header = .false.
                  if (maswrk) write (iout,20)
   20             format (/,' Additional van der Waals Parameters :',
     &                    //,6x,'Class',6x,'Radius',5x,'Epsilon',
     &                       4x,'Reduction',/)
               end if
               if (ep .eq. 0.0d0)  ep = eps(k)
               if (rdn .eq. 0.0d0)  rdn = reduct(k)
               rad(k) = rd
               eps(k) = ep
               reduct(k) = rdn
               if (maswrk) write (iout,30)  k,rd,ep,rdn
   30          format (4x,i6,1x,2f12.4,f11.3)
            else if (k .gt. maxclass) then
               if (maswrk) write (iout,40)  maxclass
   40          format (/,' KVDW  --  Only Atom Classes through',i4,
     &                    ' are Allowed')
               abort = .true.
            end if
         end if
      end do
c
c     process keywords containing specific pair vdw parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:6) .eq. 'VDWPR ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:80)
            read (string,*,err=80,end=80)  ia,ib,rd,ep
            if (header) then
               header = .false.
               if (maswrk) write (iout,50)
   50          format (/,' Additional van der Waals Parameters for',
     &                    ' Specific Pairs :',
     &                 //,6x,'Classes',5x,'Radii Sum',3x,'Epsilon',/)
            end if
            if (maswrk) write (iout,60)  ia,ib,rd,ep
   60       format (3x,2i4,2x,2f12.4)
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do k = 1, maxnvp
               if (kvpr(k).eq.'      ' .or. kvpr(k).eq.pt) then
                  kvpr(k) = pt
                  radpr(k) = rd
                  epspr(k) = ep
                  goto 80
               end if
            end do
            if (maswrk) write (iout,70)
   70       format (/,' KVDW  --  Too many Special VDW Pair',
     &                 ' Parameters')
            abort = .true.
   80       continue
         end if
      end do
c
c     process keywords containing hydrogen bonding vdw parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:6) .eq. 'HBOND ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:80)
            read (string,*,err=120,end=120)  ia,ib,rd,ep
            if (header) then
               header = .false.
               if (maswrk) write (iout,90)
   90          format (/,' Additional van der Waals Hydrogen Bonding',
     &                    ' Parameters :',
     &                 //,6x,'Classes',5x,'Radii Sum',3x,'Epsilon',/)
            end if
            if (maswrk) write (iout,100)  ia,ib,rd,ep
  100       format (3x,2i4,2x,2f12.4)
            size = 3
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do k = 1, maxnvp
               if (khb(k).eq.'      ' .or. khb(k).eq.pt) then
                  khb(k) = pt
                  radhb(k) = rd
                  epshb(k) = ep
                  goto 120
               end if
            end do
            if (maswrk) write (iout,110)
  110       format (/,' KVDW  --  Too many Hydrogen Bonding Pair',
     &                 ' Parameters')
            abort = .true.
  120       continue
         end if
      end do
c
c     get the vdw radius and well depth for each atom class
c
      do i = 1, maxclass
         if (radtyp .eq. 'SIGMA')  rad(i) = twosix * rad(i)
         if (radsiz .eq. 'DIAMETER')  rad(i) = 0.5d0 * rad(i)
         srad(i) = sqrt(rad(i))
         eps(i) = abs(eps(i))
         seps(i) = sqrt(eps(i))
      end do
c
c     use combination rules to set pairwise vdw radii sums
c
      do i = 1, maxclass
         do k = i, maxclass
            if (rad(i).eq.0.0d0 .and. rad(k).eq.0.0d0) then
               rd = 0.0d0
            else if (radrule(1:10) .eq. 'ARITHMETIC') then
               rd = rad(i) + rad(k)
            else if (radrule(1:9) .eq. 'GEOMETRIC') then
               rd = 2.0d0 * (srad(i) * srad(k))
            else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
               rd = 2.0d0 * (rad(i)**3+rad(k)**3)/(rad(i)**2+rad(k)**2)
            else
               rd = rad(i) + rad(k)
            end if
            radmin(i,k) = rd
            radmin(k,i) = rd
         end do
      end do
c
c     use combination rules to set pairwise well depths
c
      do i = 1, maxclass
         do k = i, maxclass
            if (eps(i).eq.0.0d0 .and. eps(k).eq.0.0d0) then
               ep = 0.0d0
            else if (epsrule(1:10) .eq. 'ARITHMETIC') then
               ep = 0.5d0 * (eps(i) + eps(k))
            else if (epsrule(1:9) .eq. 'GEOMETRIC') then
               ep = seps(i) * seps(k)
            else if (epsrule(1:8) .eq. 'HARMONIC') then
               ep = 2.0d0 * (eps(i)*eps(k)) / (eps(i)+eps(k))
            else if (epsrule(1:3) .eq. 'HHG') then
               ep = 4.0d0 * (eps(i)*eps(k)) / (seps(i)+seps(k))**2
            else
               ep = seps(i) * seps(k)
            end if
            epsilon(i,k) = ep
            epsilon(k,i) = ep
         end do
      end do
c
c     vdw reduction factor information for each individual atom
c
      do i = 1, n
         kred(i) = reduct(class(i))
         if (n12(i).ne.1 .or. kred(i).eq.0.0d0) then
            ired(i) = i
         else
            ired(i) = i12(1,i)
         end if
      end do
c
c     radii and well depths for special atom class pairs
c
      do i = 1, maxnvp
         if (kvpr(i) .eq. '      ')  goto 130
         ia = number (kvpr(i)(1:3))
         ib = number (kvpr(i)(4:6))
         if (rad(ia) .eq. 0.0d0)  rad(ia) = 0.001d0
         if (rad(ib) .eq. 0.0d0)  rad(ib) = 0.001d0
         radmin(ia,ib) = radpr(i)
         radmin(ib,ia) = radpr(i)
         epsilon(ia,ib) = abs(epspr(i))
         epsilon(ib,ia) = abs(epspr(i))
      end do
  130 continue
c
c     radii and well depths for hydrogen bonding pairs
c
      if (vdwtyp .eq. 'MM3-HBOND') then
         do i = 1, maxclass
            do k = 1, maxclass
               radhbnd(k,i) = 0.0d0
               epshbnd(k,i) = 0.0d0
            end do
         end do
         do i = 1, maxnhb
            if (khb(i) .eq. '      ')  goto 140
            ia = number (khb(i)(1:3))
            ib = number (khb(i)(4:6))
            radhbnd(ia,ib) = radhb(i)
            radhbnd(ib,ia) = radhb(i)
            epshbnd(ia,ib) = abs(epshb(i))
            epshbnd(ib,ia) = abs(epshb(i))
         end do
  140    continue
      end if
c
c     set the coefficients for a 2- or 4-Gaussian vdw fit
c
      if (vdwtyp .eq. 'GAUSSIAN') then
         if (gausstyp .eq. 'LJ-4') then
            ngauss = 4
            igauss(1,1) = 846706.7d0
            igauss(2,1) = 15.464405d0 * twosix**2
            igauss(1,2) = 2713.651d0
            igauss(2,2) = 7.346875d0 * twosix**2
            igauss(1,3) = -9.699172d0
            igauss(2,3) = 1.8503725d0 * twosix**2
            igauss(1,4) = -0.7154420d0
            igauss(2,4) = 0.639621d0 * twosix**2
         else if (gausstyp .eq. 'LJ-2') then
            ngauss = 2
            igauss(1,1) = 14487.1d0
            igauss(2,1) = 9.05148d0 * twosix**2
            igauss(1,2) = -5.55338d0
            igauss(2,2) = 1.22536d0 * twosix**2
         else if (gausstyp .eq. 'MM3-2') then
            ngauss = 2
            igauss(1,1) = 2438.886d0
            igauss(2,1) = 9.342616d0
            igauss(1,2) = -6.197368d0
            igauss(2,2) = 1.564486d0
         else if (gausstyp .eq. 'MM2-2') then
            ngauss = 2
            igauss(1,1) = 3423.562d0
            igauss(2,1) = 9.692821d0
            igauss(1,2) = -6.503760d0
            igauss(2,2) = 1.585344d0
         else if (gausstyp .eq. 'IN-PLACE') then
            ngauss = 2
            igauss(1,1) = 300.0d0
            igauss(2,1) = 5.432428d0 * twosix**2
            igauss(1,2) = -5.55338d0
            igauss(2,2) = 1.551167d0 * twosix**2
         end if
      end if
c
c     remove zero-sized atoms from the list of vdw sites
c
      nvdw = 0
      do i = 1, n
         if (rad(class(i)) .ne. 0.0d0) then
            nvdw = nvdw + 1
            ivdw(nvdw) = i
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine lattice  --  setup periodic boundary conditions  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "lattice" stores the periodic box dimensions, finds angles
c     to be used in computing fractional coordinates, and decides
c     between images and replicates for the boundary conditions
c
c
      subroutine lattice
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'cutoff.i'
      include 'dipole.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      integer i,j,k
      integer nxcell,nycell,nzcell
      real*8 xlimit,ylimit,zlimit
      real*8 maximage,maxcut,alpha_cos
c
c
c     set the periodic boundary box axis lengths and angles
c
      if (xbox .eq. 0.0d0) then
         use_bounds = .false.
         use_image = .false.
         use_replica = .false.
         return
      else
         use_bounds = .true.
         use_image = .true.
         if (ybox .eq. 0.0d0)  ybox = xbox
         if (zbox .eq. 0.0d0)  zbox = xbox
         if (alpha .eq. 0.0d0)  alpha = 90.0d0
         if (beta .eq. 0.0d0)  beta = 90.0d0
         if (gamma .eq. 0.0d0)  gamma = 90.0d0
      end if
c
c     set default values for the potential energy cutoffs
c
      if (chgcut.eq.0.0d0 .or. chgcut.ge.1000000.0d0)  chgcut = 10.0d0
      if (dplcut.eq.0.0d0 .or. dplcut.ge.1000000.0d0)  dplcut = 10.0d0
      if (vdwcut.eq.0.0d0 .or. vdwcut.ge.1000000.0d0)  vdwcut = 10.0d0
c
c     determine the general lattice type
c
      if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
     &                    .and. gamma.eq.90.0d0) then
         orthogonal = .true.
         monoclinic = .false.
         triclinic = .false.
      else if (alpha.eq.90.0 .and. gamma.eq.90.0) then
         orthogonal = .false.
         monoclinic = .true.
         triclinic = .false.
      else
         orthogonal = .false.
         monoclinic = .false.
         triclinic = .true.
      end if
c
c     check for use of truncated octahedron
c
      if (octahedron) then
         if (xbox.eq.ybox .and. xbox.eq.zbox .and.
     &       orthogonal) then
            orthogonal = .false.
         else
            octahedron = .false.
            write (iout,10)
   10       format (/,' LATTICE  --  Truncated Octahedron',
     &                 ' Incompatible with Defined Cell')
         end if
      end if
c
c     set the half width values for the periodic box
c
      xbox2 = 0.5d0 * xbox
      ybox2 = 0.5d0 * ybox
      zbox2 = 0.5d0 * zbox
      if (octahedron)  box34 = 0.75d0 * xbox
c
c     set values used in computation of fractional coordinates
c
      if (orthogonal .or. octahedron) then
         alpha_cos = 0.0d0
         beta_sin = 1.0d0
         beta_cos = 0.0d0
         gamma_sin = 1.0d0
         gamma_cos = 0.0d0
         beta_term = 0.0d0
         gamma_term = 1.0d0
      else if (monoclinic) then
         alpha_cos = 0.0d0
         beta_sin = sin(beta/radian)
         beta_cos = cos(beta/radian)
         gamma_sin = 1.0d0
         gamma_cos = 0.0d0
         beta_term = 0.0d0
         gamma_term = beta_sin
      else if (triclinic) then
         alpha_cos = cos(alpha/radian)
         beta_sin = sin(beta/radian)
         beta_cos = cos(beta/radian)
         gamma_sin = sin(gamma/radian)
         gamma_cos = cos(gamma/radian)
         beta_term = (alpha_cos-beta_cos*gamma_cos) / gamma_sin
         gamma_term = sqrt(1.0d0-alpha_cos**2-beta_cos**2-gamma_cos**2
     &                +2.0d0*alpha_cos*beta_cos*gamma_cos) / gamma_sin
      end if
c
c     find the maximum sphere radius inscribed in periodic box
c
      if (orthogonal) then
         xlimit = xbox2
         ylimit = ybox2
         zlimit = zbox2
      else if (monoclinic) then
         xlimit = xbox2 * beta_sin
         ylimit = ybox2
         zlimit = zbox2 * beta_sin
      else if (triclinic) then
         xlimit = xbox2 * beta_sin * gamma_sin
         ylimit = ybox2 * gamma_sin
         zlimit = zbox2 * beta_sin
      else if (octahedron) then
         xlimit = (sqrt(3.0d0)/4.0d0) * xbox
         ylimit = xlimit
         zlimit = xlimit
      end if
      maximage = min(xlimit,ylimit,zlimit)
c
c     use replicate method to handle cutoffs too large for images
c
      maxcut = vdwcut
      if (nion .ne. 0)  maxcut = max(maxcut,chgcut)
      if (ndipole .ne. 0)  maxcut = max(maxcut,dplcut)
      if (maxcut .le. maximage) then
         use_replica = .false.
      else
         use_replica = .true.
      end if
c
c     truncated octahedron cannot use the replicates method
c
      if (octahedron .and. use_replica) then
         write (iout,20)
   20    format (/,' LATTICE  --  Truncated Octahedron',
     &              ' cannot be Replicated')
         call fatal
      end if
c
c     find the number of replicates for orthogonal lattice
c
      nxcell = int(maxcut/xlimit)
      nycell = int(maxcut/ylimit)
      nzcell = int(maxcut/zlimit)
      if (maxcut .gt. nxcell*xlimit)  nxcell = nxcell + 1
      if (maxcut .gt. nycell*ylimit)  nycell = nycell + 1
      if (maxcut .gt. nzcell*zlimit)  nzcell = nzcell + 1
      xcell = nxcell * xbox
      ycell = nycell * ybox
      zcell = nzcell * zbox
c
c     check the total number of replicated unit cells
c
      ncell = nxcell*nycell*nzcell - 1
      if (ncell .gt. maxcell) then
         write (iout,30)
   30    format (/,' LATTICE  --  Increase MAXCELL or Decrease',
     &              ' the Interaction Cutoffs')
         call fatal
      end if
c
c     assign indices to the required cell replicates
c
      ncell = 0
      do k = 0, nzcell-1
         do j = 0, nycell-1
            do i = 0, nxcell-1
               if (k.ne.0 .or. j.ne.0 .or. i.ne.0) then
                  ncell = ncell + 1
                  icell(1,ncell) = i
                  icell(2,ncell) = j
                  icell(3,ncell) = k
               end if
            end do
         end do
      end do
c
c     set the half width values for the full replicated cell
c
      xcell2 = 0.5d0 * xcell
      ycell2 = 0.5d0 * ycell
      zcell2 = 0.5d0 * zcell
c
c     print a message indicating the number of replicates used
c
      if (debug) then
         write (iout,40)  nxcell,nycell,nzcell
   40    format (/,' LATTICE  --  Periodicity using ',i2,' x',i2,
     &              ' x',i2,' Set of Cell Replicates')
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine lights  --  get neighbors via method of lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "lights" computes the set of nearest neighbor interactions
c     using the method of lights algorithm
c
c     literature reference:
c
c     F. Sullivan, R. D. Mountain and J. O'Connell, "Molecular
c     Dynamics on Vector Computers", Journal of Computational
c     Physics, 61, 138-153 (1985)
c
c
      subroutine lights (nsite,map,xsort,ysort,zsort)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'shunt.i'
      integer i,j,k,nsite,map(maxlight)
      real*8 box,xcut,ycut,zcut
      real*8 xmove,ymove,zmove
      real*8 xsort(maxlight),ysort(maxlight),zsort(maxlight)
c
c
c     check that maximum number of replicates is not exceeded
c
      if (use_replica) then
         if (xcell2.gt.xbox .or. ycell2.gt.ybox
     &           .or. zcell2.gt.zbox) then
            write (iout,10)
   10       format (/,' LIGHTS  --  Cutoff Distance is Too Large',
     &                 ' for Method of Lights')
            call fatal
         end if
      end if
c
c     truncated octahedron periodicity is not handled at present
c
      if (use_bounds) then
         if (octahedron) then
            write (iout,20)
   20       format (/,' LIGHTS  --  Method of Lights unsuitable',
     &                 ' for Truncated Octahedron')
            call fatal
         end if
      end if
c
c     when using images, move coordinates into periodic cell
c
      if (use_image) then
         do i = 1, nsite
            zsort(i) = zsort(i) / gamma_term
            ysort(i) = (ysort(i) - zsort(i)*beta_term) / gamma_sin
            xsort(i) = xsort(i) - ysort(i)*gamma_cos - zsort(i)*beta_cos
            dowhile (xsort(i) .gt. xcell2)
               xsort(i) = xsort(i) - xcell
            end do
            dowhile (xsort(i) .lt. -xcell2)
               xsort(i) = xsort(i) + xcell
            end do
            dowhile (ysort(i) .gt. ycell2)
               ysort(i) = ysort(i) - ycell
            end do
            dowhile (ysort(i) .lt. -ycell2)
               ysort(i) = ysort(i) + ycell
            end do
            dowhile (zsort(i) .gt. zcell2)
               zsort(i) = zsort(i) - zcell
            end do
            dowhile (zsort(i) .lt. -zcell2)
               zsort(i) = zsort(i) + zcell
            end do
         end do
      end if
c
c     generate replica coordinates for the sort arrays
c
      if (use_replica) then
         k = nsite
         do j = 1, ncell
            xmove = icell(1,j) * xbox
            ymove = icell(2,j) * ybox
            zmove = icell(3,j) * zbox
            do i = 1, nsite
               k = k + 1
               map(k) = i
               xsort(k) = xsort(i) + xmove
               ysort(k) = ysort(i) + ymove
               zsort(k) = zsort(i) + zmove
               dowhile (xsort(k) .gt. xcell2)
                  xsort(k) = xsort(k) - xcell
               end do
               dowhile (xsort(k) .lt. -xcell2)
                  xsort(k) = xsort(k) + xcell
               end do
               dowhile (ysort(k) .gt. ycell2)
                  ysort(k) = ysort(k) - ycell
               end do
               dowhile (ysort(k) .lt. -ycell2)
                  ysort(k) = ysort(k) + ycell
               end do
               dowhile (zsort(k) .gt. zcell2)
                  zsort(k) = zsort(k) - zcell
               end do
               dowhile (zsort(k) .lt. -zcell2)
                  zsort(k) = zsort(k) + zcell
               end do
            end do
         end do
      end if
c
c     map image and replicate sites onto their parent sites
c
      nlight = (ncell+1) * nsite
      do i = 0, ncell
         j = i * nsite
         do k = 1, nsite
            map(j+k) = k
         end do
      end do
c
c     sort the coordinate components into ascending order
c
      call sort2 (nlight,xsort,locx)
      call sort2 (nlight,ysort,locy)
      call sort2 (nlight,zsort,locz)
c
c     use of replicates requires secondary sorting along x-axis
c
      if (use_replica) then
         j = 1
         do i = 1, nlight-1
            if (xsort(i+1) .ne. xsort(i)) then
               call sort5 (i-j+1,locx(j),nsite)
               j = i + 1
            end if
         end do
         call sort5 (nlight-j+1,locx(j),nsite)
      end if
c
c     index the position of each atom in the sorted coordinates
c
      do i = 1, nlight
         rgx(locx(i)) = i
         rgy(locy(i)) = i
         rgz(locz(i)) = i
      end do
c
c     set the light width based on interaction cutoff
c
      xcut = off
      ycut = off
      zcut = off
      if (use_image) then
         if (monoclinic) then
            zcut = zcut / beta_sin
            xcut = xcut + zcut*abs(beta_cos)
         else if (triclinic) then
            zcut = zcut / gamma_term
            ycut = (ycut + zcut*abs(beta_term)) / gamma_sin
            xcut = xcut + ycut*abs(gamma_cos) + zcut*abs(beta_cos)
         end if
         xcut = min(xcut,xcell2)
         ycut = min(ycut,ycell2)
         zcut = min(zcut,zcell2)
      end if
c
c     find the negative x-coordinate boundary for each atom
c
      do i = nlight, 1, -1
         k = locx(i)
         if (k .le. nsite) then
            kbx(k) = i
         end if
      end do
c
c     find the positive x-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locx(i)
         if (k .le. nsite) then
            dowhile (xsort(j)-xsort(i)+box .lt. xcut)
               if (j .eq. nlight) then
                  if (use_image) then
                     j = 0
                     box = xcell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 30
            end do
   30       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            kex(k) = j
         end if
      end do
c
c     find the negative y-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locy(i)
         if (k .le. nsite) then
            dowhile (ysort(i)-ysort(j)+box .le. ycut)
               if (j .eq. 1) then
                  if (use_image) then
                     j = nlight + 1
                     box = ycell
                  end if
               end if
               j = j - 1
               if (j .lt. 1)  goto 40
            end do
   40       continue
            j = j + 1
            if (j .gt. nlight) then
               j = 1
               box = 0.0d0
            end if
            kby(k) = j
         end if
      end do
c
c     find the positive y-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locy(i)
         if (k .le. nsite) then
            dowhile (ysort(j)-ysort(i)+box .lt. ycut)
               if (j .eq. nlight) then
                  if (use_image) then
                     j = 0
                     box = ycell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 50
            end do
   50       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            key(k) = j
         end if
      end do
c
c     find the negative z-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locz(i)
         if (k .le. nsite) then
            dowhile (zsort(i)-zsort(j)+box .le. zcut)
               if (j .eq. 1) then
                  if (use_image) then
                     j = nlight + 1
                     box = zcell
                  end if
               end if
               j = j - 1
               if (j .lt. 1)  goto 60
            end do
   60       continue
            j = j + 1
            if (j .gt. nlight) then
               j = 1
               box = 0.0d0
            end if
            kbz(k) = j
         end if
      end do
c
c     find the positive z-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locz(i)
         if (k .le. nsite) then
            dowhile (zsort(j)-zsort(i)+box .lt. zcut)
               if (j .eq. nlight) then
                  if (use_image) then
                     j = 0
                     box = zcell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 70
            end do
   70       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            kez(k) = j
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine lmqn  --  conjugate gradient optimization  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "lmqn" is a gradient optimization routine which implements
c     steepest descent, Fletcher-Reeves CG, Polak-Ribiere CG,
c     Hestenes-Stiefel CG, Powell-Beale CG with restarts, and a
c     limited memory BFGS quasi-Newton method
c
c     literature references:
c
c     David G. Luenberger, "Linear and Nonlinear Programming,
c     2nd Edition", Addison-Wesley, Reading, MA, 1984; see in
c     particular section 9-7
c
c     Stephen G. Nash and Ariela Sofer, "Linear and Nonlinear
c     Programming", McGraw-Hill, New York, 1996; see chapter 12
c
c     parameters used in the main iteration:
c
c     nvar     number of parameters in the objective function
c     method   steepest descent, FRCG, PRCG, HSCG, Powell or LMQN
c     fctmin   normal exit if function gets less than
c     grdmin   normal exit if gradient norm gets less than
c     maxiter  error return if number of iterations exceeds
c     period   restart at least every period iterations
c     iprint   print iteration results every iprint iterations
c     iwrite   call user-supplied output every iwrite iterations
c     fast     steepest descent while function decrease exceeds
c     slow     restart if relative function decrease drops below
c     epsln    error if total move is less than
c     scale    factor by which actual function has been multiplied
c     rms      factor to convert grad norm and movement to rms
c     minimum  final value of the objective function
c     ncalls   number of function/gradient (fgvalue) calls
c     niter    number of conjugate gradient iterations performed
c     status   string containing informative termination message
c
c     parameters used in the line search:
c
c     cappa    reduction in projected gradient for termination
c     stpmin   minimum allowed line search step size
c     stpmax   maximum allowed line search step size
c     angmax   maximum angle between search and -grad directions
c     intmax   maximum number of interpolations in line search
c
c     vectors stored by the routine:
c
c     x        current parameter values
c     x_old    previous parameter values
c     g        current gradient vector
c     g_old    previous gradient vector
c     s        current minus previous parameter values
c     d        current minus previous gradient values
c     p        conjugate direction search vector
c
c     requried external routines:
c
c     fgvalue    function to evaluate function and gradient values
c     writeout   subroutine to write out info about current status
c
c
ckk -minimize-      subroutine lmqn (nvar,x,minimum,grdmin,fgvalue,writeout)
      subroutine lmqn (method,nvar,x,minimum,grdmin,fgvalue,
     *                 writeout)
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'keys.i'
      include 'linmin.i'
      include 'math.i'
      include 'minima.i'
      include 'output.i'
      include 'potent.i'
      include 'scales.i'
      include 'solute.i'
casa
      include 'restrn.i'
      integer   atm
casa
      integer i,ncalls,nerror,next,ipow
      integer nvar,niter,period,nstart
      real*8 grdmin,fast,slow,epsln
      real*8 f,f_old,f_new,f_move
      real*8 rms,beta,gamma,angle,ratio
      real*8 minimum,x_move,g_rms
      real*8 fgvalue,gg,gg_old
      real*8 sg,dg,dp,sd,dd,dtg,dtpt
      real*8 x(maxvar),g(maxvar)
      real*8 x_old(maxvar),g_old(maxvar)
      real*8 p(maxvar),s(maxvar),d(maxvar)
      real*8 pt(maxvar),dt(maxvar)
      character*3 method
      character*9 blank,status
      character*20 keyword
      character*80 record
      logical restart,done
      external fgvalue,writeout
c
c
c     initialize some values to be used below
c
      if (nvar .gt. maxvar) then
         write (iout,10)
   10    format (' LMQN  --  Too many Parameters,',
     &           ' Increase value of MAXVAR')
         return
      end if
      ncalls = 0
      rms = sqrt(dble(nvar))
      if (coordtype .eq. 'cartesian') then
         rms = rms / sqrt(3.0d0)
      else if (coordtype .eq. 'rigidbody') then
         rms = rms / sqrt(6.0d0)
      end if
      blank = '         '
      status = blank
ckk -minimize-      method = 'mqn'
      if(method.eq.'   ') method = 'mqn'
      restart = .true.
      done = .false.
c
c     set default values for variable scale factors
c
      if (.not. set_scale) then
         do i = 1, nvar
            if (scale(i) .eq. 0.0d0)  scale(i) = 1.0d0
         end do
      end if
c
c     set default parameters for the optimization
c
      if (fctmin .eq. 0.0d0)  fctmin = -1000000.0d0
      if (maxiter .eq. 0)  maxiter = 1000000
      if (nextiter .eq. 0)  nextiter = 1
      if (iprint .lt. 0)  iprint = 1
      if (iwrite .lt. 0)  iwrite = 1
      fast = 0.5d0
      slow = 0.0d0
      epsln = 1.0d-16
      period = max(200,nvar)
c
c     set default parameters for the line search
c
      if (cappa .eq. 0.0d0)  cappa = 0.1d0
      if (stpmin .eq. 0.0d0)  stpmin = 1.0d-16
      if (stpmax .eq. 0.0d0)  stpmax = 5.0d0
      if (angmax .eq. 0.0d0)  angmax = 88.0d0
      if (intmax .eq. 0)  intmax = 5
c
c     search the keywords for optimization parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:17) .eq. 'STEEPEST-DESCENT ') then
            method = 'sd'
         else if (keyword(1:16) .eq. 'FLETCHER-REEVES ') then
            method = 'fr'
         else if (keyword(1:14) .eq. 'POLAK-RIBIERE ') then
            method = 'pr'
         else if (keyword(1:17) .eq. 'HESTENES-STIEFEL ') then
            method = 'hs'
         else if (keyword(1:13) .eq. 'POWELL-BEALE ') then
            method = 'pow'
         else if (keyword(1:7) .eq. 'FCTMIN ') then
            read (record(next:80),*,err=20)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (record(next:80),*,err=20)  maxiter
         else if (keyword(1:9) .eq. 'NEXTITER ') then
            read (record(next:80),*,err=20)  nextiter
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (record(next:80),*,err=20)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (record(next:80),*,err=20)  iwrite
         else if (keyword(1:8) .eq. 'STEPMAX ') then
            read (record(next:80),*,err=20)  stpmax
         else if (keyword(1:8) .eq. 'STEPMIN ') then
            read (record(next:80),*,err=20)  stpmin
         else if (keyword(1:6) .eq. 'CAPPA ') then
            read (record(next:80),*,err=20)  cappa
         else if (keyword(1:7) .eq. 'ANGMAX ') then
            read (record(next:80),*,err=20)  angmax
         else if (keyword(1:6) .eq. 'EPSLN ') then
            read (record(next:80),*,err=20)  epsln
         else if (keyword(1:5) .eq. 'FAST ') then
            read (record(next:80),*,err=20)  fast
         else if (keyword(1:5) .eq. 'SLOW ') then
            read (record(next:80),*,err=20)  slow
         else if (keyword(1:7) .eq. 'PERIOD ') then
            read (record(next:80),*,err=20)  period
         else if (keyword(1:7) .eq. 'INTMAX ') then
            read (record(next:80),*,err=20)  intmax
         end if
   20    continue
      end do
c
c     get initial function and gradient values
c
      niter = nextiter - 1
      maxiter = niter + maxiter
      ncalls = ncalls + 1
      f = fgvalue (x,g)
      g_rms = 0.0d0
      f_move = 0.0d0
casa
c     write(6,*) 'wwwf2',npfix,(pfix(i),i=1,npfix),(ipfix(i),i=1,npfix)
      do i = 1, npfix
         if (pfix(i) .ge. RESTRAINV) then
            atm = (ipfix(i) - 1) * 3
c       write(iout,*) 'atm and g is',  atm, g(atm+1), g(atm+2), g(atm+3)
            g(atm+1) = 0.0d0
            g(atm+2) = 0.0d0
            g(atm+3) = 0.0d0
         end if
      end do
casa
      do i = 1, nvar
         x_old(i) = x(i)
         g_old(i) = g(i)
         g_rms = g_rms + (g(i)*scale(i))**2
         f_move = f_move + g(i)**2
      end do
      g_rms = sqrt(g_rms) / rms
      f_move = 0.5d0 * stpmax * sqrt(f_move)
c
c     print header information prior to iterations
c
      if (iprint .gt. 0) then
         if (method .eq. 'sd') then
            write (iout,30)
   30       format (/,' Steepest Descent Gradient Optimization :')
         else if (method .eq. 'fr') then
            write (iout,40)
   40       format (/,' Fletcher-Reeves Conjugate Gradient',
     &                 ' Optimization :')
         else if (method .eq. 'pr') then
            write (iout,50)
   50       format (/,' Polak-Ribiere Conjugate Gradient',
     &                 ' Optimization :')
         else if (method .eq. 'hs') then
            write (iout,60)
   60       format (/,' Hestenes-Stiefel Conjugate Gradient',
     &                 ' Optimization :')
         else if (method .eq. 'pow') then
            write (iout,70)
   70       format (/,' Powell-Beale Conjugate Gradient Optimization :')
         else if (method .eq. 'mqn') then
            write (iout,80)
   80       format (/,' Memoryless BFGS Quasi-Newton Optimization :')
         end if
         write (iout,90)
   90    format (/,' CG Iter    F Value      G RMS     F Move',
     &             '    X Move    Angle  FG Call  Comment',/)
         if (f.lt.1.0d7 .and. f.gt.-1.0d6 .and. g_rms.lt.1.0d5) then
            write (iout,100)  niter,f,g_rms,ncalls
  100       format (i6,f13.4,f11.4,30x,i7)
         else
            write (iout,110)  niter,f,g_rms,ncalls
  110       format (i6,d13.4,d11.4,30x,i7)
         end if
      end if
c
c     tests of the various termination criteria
c
      if (niter .ge. maxiter) then
         status = 'IterLimit'
         done = .true.
      end if
      if (f .le. fctmin) then
         status = 'SmallFct '
         done = .true.
      end if
      if (g_rms .le. grdmin) then
         status = 'SmallGrad'
         done = .true.
      end if
c
c     start of a new conjugate gradient iteration
c
      dowhile (.not. done)
         niter = niter + 1
         if (status .eq. blank)  nerror = 0
c
c     if using GB/SA solvation, update Born radii occasionally
c
         if (use_gbsa) then
            if (mod(niter,50) .eq. 0) then
               reborn = 1
               call born
               f = fgvalue (x,g)
            end if
            reborn = 0
         end if
c
c     compute the next search direction using Steepest Descent,
c     Fletcher-Reeves CG, Polak-Ribiere CG or Memoryless BFGS
c
casa
c     write(6,*) 'wwwf1',npfix,(pfix(i),i=1,npfix),(ipfix(i),i=1,npfix)
         do i = 1, npfix
            if (pfix(i) .ge. RESTRAINV) then
               atm = (ipfix(i) - 1) * 3
c       write(iout,*) 'atm and g is',  atm, g(atm+1), g(atm+2), g(atm+3)
               g(atm+1) = 0.0d0
               g(atm+2) = 0.0d0
               g(atm+3) = 0.0d0
            end if
         end do
casa
  120    continue
         status = blank
         if (method.eq.'sd' .or. restart) then
            do i = 1, nvar
               p(i) = -g(i)
            end do
            nstart = niter
            restart = .false.
         else if (method .eq. 'fr') then
            gg = 0.0d0
            gg_old = 0.0d0
            do i = 1, nvar
               gg = gg + g(i)*g(i)
               gg_old = gg_old + g_old(i)*g_old(i)
            end do
            beta = gg / gg_old
            do i = 1, nvar
               p(i) = -g(i) + beta*p(i)
            end do
         else if (method .eq. 'pr') then
            dg = 0.0d0
            gg_old = 0.0d0
            do i = 1, nvar
               dg = dg + d(i)*g(i)
               gg_old = gg_old + g_old(i)*g_old(i)
            end do
            beta = dg / gg_old
            do i = 1, nvar
               p(i) = -g(i) + beta*p(i)
            end do
         else if (method .eq. 'hs') then
            dg = 0.0d0
            dp = 0.0d0
            do i = 1, nvar
               dg = dg + d(i)*g(i)
               dp = dp + d(i)*p(i)
            end do
            beta = dg / dp
            do i = 1, nvar
               p(i) = -g(i) + beta*p(i)
            end do
         else if (method .eq. 'pow') then
            dg = 0.0d0
            dp = 0.0d0
            do i = 1, nvar
               dg = dg + d(i)*g(i)
               dp = dp + d(i)*p(i)
            end do
            beta = dg / dp
            if (niter .eq. ipow) then
               gamma = 0.0d0
            else
               dtg = 0.0d0
               dtpt = 0.0d0
               do i = 1, nvar
                  dtg = dtg + dt(i)*g(i)
                  dtpt = dtpt + dt(i)*pt(i)
               end do
               gamma = dtg / dtpt
            end if
            do i = 1, nvar
               p(i) = -g(i) + beta*p(i) + gamma*pt(i)
            end do
         else if (method .eq. 'mqn') then
            sg = 0.0d0
            dg = 0.0d0
            dd = 0.0d0
            sd = 0.0d0
            do i = 1, nvar
               sg = sg + s(i)*g(i)
               dg = dg + d(i)*g(i)
               dd = dd + d(i)*d(i)
               sd = sd + s(i)*d(i)
            end do
            do i = 1, nvar
               p(i) = -g(i) + (d(i)*sg+s(i)*dg)/sd
     &                   - (1.0d0+dd/sd)*(s(i)*sg/sd)
            end do
         end if
c
c     perform line search along the new conjugate direction
c
         f_old = f
         call search (nvar,f,g,x,p,f_move,angle,ncalls,fgvalue,status)
         f_new = f
c
c     if angle between the search direction and the negative
c     gradient vector was too large, use steepest descent
c
         if (status .eq. 'WideAngle') then
            restart = .true.
            goto 120
         end if
c
c     special test for reset of the Powell method
c
         if (method .eq. 'pow') then
            restart = .false.
            gg = 0.0d0
            gg_old = 0.0d0
            do i = 1, nvar
               gg = gg + g(i)*g(i)
               gg_old = gg_old + g(i)*g_old(i)
            end do
            ratio = gg_old / gg
            if (niter.eq.1 .or. ratio.ge.0.2d0) then
               ipow = niter + 1
               do i = 1, nvar
                  pt(i) = p(i)
                  dt(i) = g(i) - g_old(i)
               end do
               nstart = niter
               status = 'Resetting'
            end if
         end if
c
c     update variables based on results of this iteration
c
         f_move = f_old - f_new
         x_move = 0.0d0
         g_rms = 0.0d0
         do i = 1, nvar
            s(i) = x(i) - x_old(i)
            d(i) = g(i) - g_old(i)
            x_move = x_move + (s(i)/scale(i))**2
            g_rms = g_rms + (g(i)*scale(i))**2
            x_old(i) = x(i)
            g_old(i) = g(i)
         end do
         x_move = sqrt(x_move) / rms
         if (coordtype .eq. 'internal') then
            x_move = radian * x_move
         end if
         g_rms = sqrt(g_rms) / rms
c
c     test for error/restart due to line search problems
c
         if (status .eq. 'BadIntpln') then
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
         end if
         if (status .eq. 'IntplnErr') then
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
            do i = 1, nvar
               x(i) = x_old(i)
               g(i) = g_old(i)
            end do
         end if
c
c     test for error/restart due to lack of movement
c
         if (x_move .lt. epsln) then
            status = 'SmallMove'
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
         end if
c
c     test for error/restart due to function increase
c
         if (f_move .lt. 0.0d0) then
            status = 'Increase '
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
            do i = 1, nvar
               x(i) = x_old(i)
               g(i) = g_old(i)
            end do
         end if
c
c     test for error/restart due to slow progress
c
         ratio = f_move / abs(f_new-fctmin)
         if (abs(ratio) .lt. slow) then
            status = 'SlowDecr '
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
         end if
c
c     test for restart due to fast progress
c
         if (ratio .gt. fast) then
            status = 'FastDecr '
            nerror = 0
            restart = .true.
         end if
c
c     test for restart based on iterations since last restart
c
         if (niter-nstart .ge. period) then
            status = 'Periodic '
            restart = .true.
         end if
c
c     test for too many total iterations
c
         if (niter .ge. maxiter) then
            status = 'IterLimit'
            done = .true.
         end if
c
c     test the normal termination criteria
c
         if (f .le. fctmin) then
            status = 'SmallFct '
            done = .true.
         end if
         if (g_rms .le. grdmin) then
            status = 'SmallGrad'
            done = .true.
         end if
c
c     print intermediate results every few iterations
c
         if (iprint .gt. 0) then
            if (done .or. mod(niter,iprint).eq.0) then
               if (f.lt.1.0d7 .and. f.gt.-1.0d6 .and.
     &             g_rms.lt.1.0d5 .and. f_move.lt.1.0d5) then
                  write (iout,130)  niter,f,g_rms,f_move,
     &                              x_move,angle,ncalls,status
  130             format (i6,f13.4,f11.4,f11.4,f10.4,f9.2,i7,3x,a9)
               else
                  write (iout,140)  niter,f,g_rms,f_move,
     &                              x_move,angle,ncalls,status
  140             format (i6,d13.4,d11.4,d11.4,f10.4,f9.2,i7,3x,a9)
               end if
            end if
         end if
c
c     write intermediate results every few iterations
c
         if (iwrite .gt. 0) then
            if (done .or. mod(niter,iwrite).eq.0) then
               call writeout (x,niter)
            end if
         end if
      end do
c
c     set final value of the objective function
c
      minimum = f
      if (iprint .gt. 0) then
         if (status.eq.'SmallGrad' .or. status.eq.'SmallFct ') then
            write (iout,150)  status
  150       format (/,' LMQN  --  Normal Termination due to ',a9)
         else
            write (iout,160)  status
  160       format (/,' LMQN  --  Incomplete Convergence due to ',a9)
         end if
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine lmstep  --  finds Levenberg-Marquardt step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "lmstep" computes the Levenberg-Marquardt step during a
c     nonlinear least squares calculation; this version is based
c     upon ideas from the Minpack routine LMPAR together with
c     with the internal doubling strategy of Dennis and Schnabel
c
c     arguments and variables :
c
c     n        dimension of the problem
c     ga       real vector of length n containing the gradient
c                of the residual vector
c     a        real n by n array which on input contains in the full
c                upper triangle the upper triangle of the matrix r
c                resulting from the QR factorization of the Jacobian;
c                on output the full upper triangle is unaltered, and
c                the strict lower triangle contains the strict lower
c                triangle of the lower triangular matrix l which is
c                equal to the Cholesky factor of (j**t)*j + amu*xscale
c     lda      leading dimension of a exactly as specified in the
c                dimension statement of the calling program
c     ipvt     integer array of length n containing the pivoting
c                infomation from the QR factorization routine
c     xscale   real vector of length n containing the diagonal
c                scaling matrix for the variables
c     qtf      real vector containing the first n elements of
c                Q(transpose) * (the scaled residual vector)
c     amu      scalar containing an initial estimate of the
c                Levenberg-Marquardt parameter on input; on
c                output, amu contains the final estimate of
c                the Levenberg-Marquardt parameter
c     first    logical variable set true only if this is the first
c                call to this routine in this iteration
c     sa       real vector of length n containing the
c                Levenberg-Marquardt step
c     gnstep   real vector of length n containing the
c                Gauss-Newton step
c     gauss    logical variable which is true if the Gauss-Newton
c                step is acceptable, and is false otherwise
c     diag     vector of length n containing the diagonal elements
c                of the Cholesky factor of (j**t)*j + amu*xscale
c
c
      subroutine lmstep (n,ga,a,lda,ipvt,xscale,qtf,stpmax,
     &                       delta,amu,first,sa,gauss)
      implicit none
      integer maxlsq,maxrsd
      parameter (maxlsq=50)
      parameter (maxrsd=100)
      integer i,j,k,n,lda,ipvt(maxlsq),nsing
      real*8 stpmax,delta,amu
      real*8 ga(maxlsq),a(lda,maxlsq),xscale(maxlsq)
      real*8 qtf(maxlsq),sa(maxlsq),gnstep(maxlsq)
      real*8 diag(maxlsq),work1(maxlsq),work2(maxlsq)
      real*8 alow,alpha,amulow,amuhi,beta,small,sum
      real*8 deltap,gnleng,phi,phip,phipi,sgnorm
      real*8 high,tiny,stplen,temp,precise
      logical first,gauss,done
      save deltap,phi,phip,nsing
      external precise
c
c
c     set smallest floating point magnitude and spacing
c
      tiny = precise (1)
      small = precise (2)
c
c     if initial trust region is not provided by the user,
c     compute and use the length of the Cauchy step given
c     by beta = norm2(r*trans(p)*d**(-2)*g)**2
c
      if (delta .eq. 0.0d0) then
         amu = 0.0d0
         do i = 1, n
            work1(i) = ga(i) / xscale(i)
         end do
         alpha = 0.0d0
         do i = 1, n
            alpha = alpha + work1(i)**2
         end do
         beta = 0.0d0
         do i = 1, n
            temp = 0.0d0
            do j = i, n
               k = ipvt(j)
               temp = temp + a(i,j)*ga(k)/xscale(k)**2
            end do
            beta = beta + temp**2
         end do
         if (beta .le. tiny) then
            delta = alpha * sqrt(alpha)
         else
            delta = alpha * sqrt(alpha)/beta
         end if
         delta = min(delta,stpmax)
      end if
c
c     the following is done only on the first time through
c     this iteration:  (1) compute the Gauss-Newton step;
c     if the Jacobian is rank-deficient, obtain a least
c     squares solution; (2) compute the length of the scaled
c     Gauss-Newton step; (3) compute the norm of the scaled
c     gradient used in computing an upper bound for "amu"
c
      if (first) then
         nsing = n
         do j = 1, n
            if (a(j,j).eq.0.0d0 .and. nsing.eq.n)  nsing = j - 1
            if (nsing .lt. n)  work1(j) = 0.0d0
         end do
         work1(nsing) = qtf(nsing) / a(nsing,nsing)
         do j = nsing-1, 1, -1
            sum = 0.0d0
            do i = j+1, nsing
               sum = sum + a(j,i)*work1(i)
            end do
            work1(j) = (qtf(j)-sum) / a(j,j)
         end do
         do j = 1, n
            gnstep(ipvt(j)) = -work1(j)
         end do
c
c     find the length of scaled Gauss-Newton step
c
         do j = 1, n
            work1(j) = xscale(j) * gnstep(j)
         end do
         gnleng = 0.0d0
         do j = 1, n
            gnleng = gnleng + work1(j)**2
         end do
         gnleng = sqrt(gnleng)
c
c     find the length of the scaled gradient
c
         do j = 1, n
            work1(j) = ga(j) / xscale(j)
         end do
         sgnorm = 0.0d0
         do j = 1, n
            sgnorm = sgnorm + work1(j)**2
         end do
         sgnorm = sqrt(sgnorm)
      end if
c
c     set the bounds on the computed step
c
      high = 1.5d0
      alow = 0.75d0
c
c     check to see if the Gauss-Newton step is acceptable
c
      if (gnleng .le. high*delta) then
         gauss = .true.
         do j = 1, n
            sa(j) = gnstep(j)
         end do
         amu = 0.0d0
         delta = min(delta,gnleng)
c
c     the Gauss-Newton step is rejected, find a nontrivial step;
c     first compute a starting value of "amu" if previous step
c     was not a Gauss-Newton step
c
      else
         gauss = .false.
         if (amu .gt. 0.0d0)
     &      amu = amu - ((phi+deltap)/delta)*(((deltap-delta)+phi)/phip)
         phi = gnleng - delta
c
c     if the Jacobian is not rank deficient, the Newton step
c     provides a lower bound for "amu"; else set bound to zero
c
         if (nsing .eq. n) then
            if (first) then
               first = .false.
               do j = 1, n
                  k = ipvt(j)
                  work1(j) = gnstep(k) * xscale(k)**2
               end do
c
c     obtain trans(r**-1)*(trans(p)*s) by solving the
c     system of equations trans(r)*work1 = work1
c
               work1(n) = work1(n) / a(n,n)
               do j = n-1, 1, -1
                  sum = 0.0d0
                  do i = j+1, n
                     sum = sum + a(j,i)*work1(i)
                  end do
                  work1(j) = (work1(j)-sum) / a(j,j)
               end do
               phipi = 0.0d0
               do j = 1, n
                  phipi = phipi - work1(j)**2
               end do
               phipi = phipi / gnleng
            end if
            amulow = -phi / phipi
         else
            first = .false.
            amulow = 0.0d0
         end if
         amuhi = sgnorm / delta
c
c     iterate until a satisfactory "amu" is generated
c
         done = .false.
         dowhile (.not. done)
            if (amu.lt.amulow .or. amu.gt.amuhi) then
               amu = max(sqrt(amulow*amuhi),0.001d0*amuhi)
            end if
            temp = sqrt(amu)
            do j = 1, n
               work1(j) = temp * xscale(j)
            end do
c
c     solve the damped least squares system for the value of the
c     Levenberg-Marquardt step using More's Minpack technique
c
            call qrsolv (n,a,lda,ipvt,work1,qtf,sa,diag,work2)
            do j = 1, n
               sa(j) = -sa(j)
            end do
            do j = 1, n
               work2(j) = xscale(j) * sa(j)
            end do
            stplen = 0.0d0
            do j = 1, n
               stplen = stplen + work2(j)**2
            end do
            stplen = sqrt(stplen)
            phi = stplen - delta
            do j = 1, n
               k = ipvt(j)
               work1(j) = xscale(k) * work2(k)
            end do
            do j = 1, n
               if (abs(diag(j)) .ge. tiny) then
                  work1(j) = work1(j) / diag(j)
               end if
               if (j .lt. n) then
                  do i = j+1, n
                     work1(i) = work1(i) - work1(j)*a(i,j)
                  end do
               end if
            end do
            phip = 0.0d0
            do j = 1, n
               phip = phip - work1(j)**2
            end do
            phip = phip / stplen
c
c     check to see if the step is acceptable; if not,
c     update amulow, amuhi and amu for next iteration
c
            if ((stplen.ge.alow*delta.and.stplen.le.high*delta)
     &            .or. (amuhi-amulow).le.small) then
               done = .true.
            else
               amulow = max(amulow,amu-(phi/phip))
               if (phi .lt. 0.0d0)  amuhi = amu
               amu = amu - (stplen/delta)*(phi/phip)
            end if
         end do
      end if
      deltap = delta
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine lowcase  --  convert string to all lower case  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "lowcase" converts a text string to all lower case letters
c
c
      subroutine lowcase (string)
      implicit none
      integer i,length,code,ichar
      character*1 char
      character*(*) string
c
c
c     move through the string one character at a time,
c     converting upper case letters to lower case
c
      length = len(string)
      do i = 1, length
         code = ichar(string(i:i))
         if (code.ge.65 .and. code.le.90)
     &      string(i:i) = char(code+32)
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine makeint  --  convert Cartesian to internal  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "makeint" converts Cartesian to internal coordinates where
c     selection of internal coordinates is controlled by "mode"
c
c        mode = 0     automatic internal coordinates
c        mode = 1     manual selection of coordinates
c        mode = 2     use existing structure as a template
c        mode = 3     use dihedral angles in all cases
c
c
      subroutine makeint (mode)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'math.i'
      include 'zclose.i'
      include 'zcoord.i'
      integer i,j,i1,i2,i3,i4,i5
      integer adjacent,trial,mode,next
      integer iz0(0:maxatm),iz1(maxatm)
      real*8 bndangle,dihedral,sign
      character*1 answer,default
      character*8 phrase
      character*80 record
      logical more
c
c
c     zero out the values of the local defining atoms
c
      i1 = 0
      i2 = 0
      i3 = 0
      i4 = 0
      i5 = 0
      iz0(0) = 0
      do i = 1, n
         iz0(i) = 0
         iz1(i) = 0
      end do
c
c     zero out the final internal coordinate values
c
      if (mode .ne. 2) then
         do i = 1, n
            do j = 1, 4
               iz(j,i) = 0
            end do
         end do
      end if
      do i = 1, n
         zbond(i) = 0.0d0
         zang(i) = 0.0d0
         ztors(i) = 0.0d0
      end do
c
c     first, decide which of the atoms to define next
c
      do i = 1, n
         if (mode .eq. 1) then
            trial = i1 + 1
   10       continue
            write (iout,20)  trial
   20       format (/,' Atom Number to be Defined [',i5,'] :  ',$)
            read (input,30,err=10)  i1
   30       format (i10)
            if (i1 .eq. 0)  i1 = trial
            if (iz0(i1) .ne. 0) then
               write (iout,40)
   40          format (/,' Already Defined that Atom; Choose Another')
               if (i1 .eq. trial)  trial = trial + 1
               goto 10
            end if
         else
            i1 = i
         end if
c
c     define the bond length for the current atom
c
         if (i .ge. 2) then
            if (mode .ne. 2) then
               i2 = adjacent (i1,0,mode,more,iz0,iz1)
               if (i2 .eq. 0) then
                  write (iout,50)  i1
   50             format (/,' MAKEINT  --  Connectivity Error',
     &                       ' in defining Atom',i6)
                  call fatal
               end if
            else if (mode .eq. 2) then
               i2 = iz(1,i1)
            end if
            zbond(i1) = sqrt((x(i1)-x(i2))**2 + (y(i1)-y(i2))**2
     &                             + (z(i1)-z(i2))**2)
         end if
c
c     define the bond angle for the current atom
c
         if (i .ge. 3) then
            if (mode .ne. 2) then
               i3 = adjacent (i2,i1,mode,more,iz0,iz1)
               if (i3 .eq. 0) then
                  write (iout,60)  i1
   60             format (/,' MAKEINT  --  Connectivity Error',
     &                       ' in defining Atom',i6)
                  call fatal
               end if
            else if (mode .eq. 2) then
               i3 = iz(2,i1)
            end if
            zang(i1) = bndangle (i1,i2,i3)
         end if
c
c     decide whether to use a dihedral or second bond angle;
c     then find the value of the angle
c
         if (i .ge. 4) then
            if (mode .eq. 1) then
               if (more) then
                  phrase = 'D or [B]'
                  default = 'B'
               else
                  phrase = '[D] or B'
                  default = 'D'
               end if
               write (iout,70)  phrase
   70          format (/,' Specify with Dihedral Angle or Second',
     &                    ' Bond Angle (',a8,') :  ',$)
               read (input,80)  record
   80          format (a80)
               next = 1
               call gettext (record,answer,next)
               call upcase (answer)
               if (answer.ne.'B' .and. answer.ne.'D')  answer = default
            else if (mode .eq. 0) then
               if (more) then
                  answer = 'B'
               else
                  answer = 'D'
               end if
            else if (mode .eq. 2) then
               if (iz(4,i1) .eq. 0) then
                  answer = 'D'
               else
                  answer = 'B'
               end if
            else if (mode .eq. 3) then
               answer = 'D'
            end if
            if (answer .eq. 'B') then
               if (mode .ne. 2) then
                  i4 = adjacent (i2,i3,mode,more,iz0,iz1)
                  if (i4 .eq. 0) then
                     write (iout,90)  i1
   90                format (/,' MAKEINT  --  Connectivity Error',
     &                          ' in defining Atom',i6)
                     call fatal
                  end if
               else if (mode .eq. 2) then
                  i4 = iz(3,i1)
               end if
               ztors(i1) = bndangle (i1,i2,i4)
               i5 = 1
               sign = dihedral (i1,i2,i3,i4)
               if (sign .gt. 0.0d0)  i5 = -1
            else if (answer .eq. 'D') then
               if (mode .ne. 2) then
                  i4 = adjacent (i3,i2,mode,more,iz0,iz1)
                  if (i4 .eq. 0) then
                     write (iout,100)  i1
  100                format (/,' MAKEINT  --  Connectivity Error',
     &                          ' in defining Atom',i6)
                     call fatal
                  end if
               else if (mode .eq. 2) then
                  i4 = iz(3,i1)
               end if
               i5 = 0
               ztors(i1) = dihedral (i1,i2,i3,i4)
            end if
         end if
c
c     transfer defining atoms to permanent array;
c     mark the current atom as finished
c
         iz(1,i1) = iz0(i2)
         iz(2,i1) = iz0(i3)
         iz(3,i1) = iz0(i4)
         iz(4,i1) = i5
         iz0(i1) = i
         iz1(i1) = i2
      end do
c
c     add any bonds needed to make ring closures
c
      nadd = 0
      do i = 1, n
         do j = 1, n12(i)
            if (iz0(i) .lt. iz0(i12(j,i)) .and.
     &          iz1(i12(j,i)) .ne. i) then
               nadd = nadd + 1
               iadd(1,nadd) = iz0(i)
               iadd(2,nadd) = iz0(i12(j,i))
            end if
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine makepdb  --  convert Cartesian to PDB format  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "makexyz" converts a set of Cartesian coordinates to Protein
c     Data Bank format with special handling for systems consisting
c     of polypeptide chains, ligands and water molecules
c
c
      subroutine makepdb
      implicit none
      include 'sizes.i'
c     include 'aminos.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'molcul.i'
      include 'pdb.i'
      include 'sequen.i'
      integer i,j,k,m,kp
      integer iseq,kseq,start,stop
      integer pdbnum,justify,freeunit
      integer noxy,nhydro,cbi
      integer ni(maxatm),cai(maxatm)
      integer ci(maxatm),oi(maxatm)
      character*3 resname
      character*4 atmname
      character*7 moltyp
      character*60 seqfile
      logical exist,nbone,obone
      logical water(maxatm),hetmol(maxatm)
      include 'aminos.i'
c
c
c     get atomic number of each atom and count the molecules
c
c     assume molecule is a polypeptide if sequence file exists
c
      iseq = freeunit ()
      seqfile = filename(1:leng)//'.seq'
      call version (seqfile,'old')
      inquire (file=seqfile,exist=exist)
      if (exist) then
         moltyp = 'peptide'
         open (unit=iseq,file=seqfile,status='old')
         rewind (unit=iseq)
         call readseq (iseq)
         close (unit=iseq)
      else
         moltyp = 'generic'
      end if
c
c     check each molecule to see if it is a water molecule
c
      do i = 1, nmol
         water(i) = .false.
         if (imol(2,i)-imol(1,i) .eq. 2) then
            noxy = 0
            nhydro = 0
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               if (atomic(k) .eq. 8)  noxy = noxy + 1
               if (atomic(k) .eq. 1)  nhydro = nhydro + 1
            end do
            if (noxy.eq.1 .and. nhydro.eq.2)  water(i) = .true.
         end if
      end do
c
c     for general structures, transfer each atom to the PDB format
c
      if (moltyp .eq. 'generic') then
         do i = 1, nmol
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               atmname = ' '//name(k)
               if (water(i)) then
                  resname = 'HOH'
               else
                  justify = 0
                  call numeral (type(k),resname,justify)
               end if
               pdbnum = i
               call pdbatm (atmname,resname,pdbnum,k)
               pdbtyp(npdb) = 'HETATM'
            end do
         end do
         do i = 1, nmol
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               kp = pdblist(k)
               npdb12(kp) = n12(k)
               do m = 1, n12(k)
                  ipdb12(m,kp) = pdblist(i12(m,k))
               end do
            end do
         end do
      end if
c
c     find the amide nitrogens and all other backbone atoms
c
      if (moltyp .eq. 'peptide') then
         call attach
         kseq = 1
         do i = 1, n
            resname =  amino(seqtyp(kseq))
            if (resname .eq. 'FOR') then
               if (atomic(i) .eq. 6) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     else if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(kseq) = 0
                     cai(kseq) = i
                     ci(kseq) = i
                     oi(kseq) = i + 1
                     kseq = kseq + 1
                  end if
               end if
            else if (resname .eq. 'ACE') then
               if (atomic(i) .eq. 6) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     else if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(kseq) = 0
                     cai(kseq) = i
                     ci(kseq) = i + 1
                     oi(kseq) = i + 2
                     kseq = kseq + 1
                  end if
               end if
            else if (resname .eq. 'NH2') then
               if (atomic(i) .eq. 7) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(kseq) = i
                     cai(kseq) = 0
                     ci(kseq) = 0
                     oi(kseq) = 0
                     kseq = kseq + 1
                  end if
               end if
            else if (resname .eq. 'NME') then
               if (atomic(i) .eq. 7) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(kseq) = i
                     cai(kseq) = i + 1
                     ci(kseq) = 0
                     oi(kseq) = 0
                     kseq = kseq + 1
                  end if
               end if
            else
               if (atomic(i) .eq. 7) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     else if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(kseq) = i
                     cai(kseq) = i + 1
                     ci(kseq) = i + 2
                     oi(kseq) = i + 3
                     kseq = kseq + 1
                  end if
               end if
            end if
            if (kseq .gt. nseq)  goto 10
         end do
   10    continue
      end if
c
c     get all of the atoms for each amino acid residue in order
c
      if (moltyp .eq. 'peptide') then
         npdb = 0
         do m = 1, nchain
            start = ichain(1,m)
            stop = ichain(2,m)
            do i = start, stop
               resname = amino(seqtyp(i))
               if (resname .eq. 'FOR') then
                  call pdbatm (' C  ',resname,i,ci(i))
                  call pdbatm (' O  ',resname,i,oi(i))
               else if (resname .eq. 'ACE') then
                  call pdbatm (' CH3',resname,i,cai(i))
                  call pdbatm (' C  ',resname,i,ci(i))
                  call pdbatm (' O  ',resname,i,oi(i))
               else if (resname .eq. 'NH2') then
                  call pdbatm (' N  ',resname,i,ni(i))
               else if (resname .eq. 'NME') then
                  call pdbatm (' N  ',resname,i,ni(i))
                  call pdbatm (' CH3',resname,i,cai(i))
               else
                  call pdbatm (' N  ',resname,i,ni(i))
                  call pdbatm (' CA ',resname,i,cai(i))
                  call pdbatm (' C  ',resname,i,ci(i))
                  call pdbatm (' O  ',resname,i,oi(i))
               end if
               call getside (resname,i,ci(i),cai(i),cbi)
               if (resname .eq. 'CYS') then
                  do j = 1, n13(cbi)
                     if (atomic(i13(j,cbi)) .eq. 16)  resname = 'CYX'
                  end do
               end if
               if (i.eq.stop .and. ci(i).ne.0) then
                  do j = 1, n12(ci(i))
                     k = i12(j,ci(i))
                     if (atomic(k).eq.8 .and. k.ne.oi(i)) then
                        call pdbatm (' OXT',resname,i,k)
                        goto 20
                     end if
                  end do
   20             continue
               end if
               call gethydro (resname,i,m,ni(i),cai(i),cbi)
            end do
         end do
      end if
c
c     get any ligands and waters following polypeptide chains
c
      if (moltyp .eq. 'peptide') then
         do i = 1, nmol
            hetmol(i) = .true.
         end do
         do i = 1, npdb
            hetmol(molcule(i)) = .false.
         end do
         do i = 1, nmol
            if (hetmol(i)) then
               do j = imol(1,i), imol(2,i)
                  k = kmol(j)
                  atmname = ' '//name(k)
                  if (water(i)) then
                     resname = 'HOH'
                  else
                     justify = 0
                     call numeral (type(k),resname,justify)
                  end if
                  pdbnum = nseq + i - 1
                  call pdbatm (atmname,resname,pdbnum,k)
                  pdbtyp(npdb) = 'HETATM'
               end do
            end if
         end do
         do i = 1, nmol
            if (hetmol(i)) then
               do j = imol(1,i), imol(2,i)
                  k = kmol(j)
                  kp = pdblist(k)
                  npdb12(kp) = n12(k)
                  do m = 1, n12(k)
                     ipdb12(m,kp) = pdblist(i12(m,k))
                  end do
               end do
            end if
         end do
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine pdbatm  --  add atom to the PDB file  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "pdbatm" adds an atom to the Protein Data Bank file
c
c
      subroutine pdbatm (atmname,resname,ires,icoord)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'pdb.i'
      integer ires,icoord
      character*3 resname
      character*4 atmname
c
c
c     for each atom set the sequential number, record type, atom
c     name, residue name, residue number and atomic coordinates
c
      if (icoord .ne. 0) then
         npdb = npdb + 1
         pdbtyp(npdb) = 'ATOM  '
         atmnam(npdb) = atmname
         resnam(npdb) = resname
         resnum(npdb) = ires
         xpdb(npdb) = x(icoord)
         ypdb(npdb) = y(icoord)
         zpdb(npdb) = z(icoord)
         npdb12(npdb) = 0
         pdblist(icoord) = npdb
      end if
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine getside  --  extract side chain atoms  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "getside" finds the side chain heavy atoms for a single
c     amino acid residue and copies the names and coordinates
c     to the Protein Data Bank file
c
c
      subroutine getside (resname,ires,ci,cai,cbi)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      integer i,j,ires,ci,cai,cbi
      character*3 resname
c
c
c     if residue is glycine or a cap, there is no side chain
c
      cbi = 0
      if (resname .eq. 'GLY')  return
      if (resname .eq. 'UNK')  return
      if (resname .eq. 'FOR')  return
      if (resname .eq. 'ACE')  return
      if (resname .eq. 'NH2')  return
      if (resname .eq. 'NME')  return
c
c     find the beta carbon atom for the current residue
c
      do i = 1, n
         if (i.ne.ci .and. atomic(i).eq.6) then
            do j = 1, 4
               if (i12(j,i) .eq. cai) then
                  cbi = i
                  if (resname .ne. 'AIB') then
                     call pdbatm (' CB ',resname,ires,cbi)
                  else
                     call pdbatm (' CB1',resname,ires,cbi)
                  end if
                  goto 10
               end if
            end do
         end if
      end do
   10 continue
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         continue
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         continue
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         call pdbatm (' CG1',resname,ires,cbi+1)
         call pdbatm (' CG2',resname,ires,cbi+2)
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD1',resname,ires,cbi+2)
         call pdbatm (' CD2',resname,ires,cbi+3)
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         call pdbatm (' CG1',resname,ires,cbi+1)
         call pdbatm (' CG2',resname,ires,cbi+2)
         call pdbatm (' CD1',resname,ires,cbi+3)
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         call pdbatm (' OG ',resname,ires,cbi+1)
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         call pdbatm (' OG1',resname,ires,cbi+1)
         call pdbatm (' CG2',resname,ires,cbi+2)
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         call pdbatm (' SG ',resname,ires,cbi+1)
c
c     cysteine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         call pdbatm (' SG ',resname,ires,cbi+1)
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD ',resname,ires,cbi+2)
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD1',resname,ires,cbi+2)
         call pdbatm (' CD2',resname,ires,cbi+3)
         call pdbatm (' CE1',resname,ires,cbi+4)
         call pdbatm (' CE2',resname,ires,cbi+5)
         call pdbatm (' CZ ',resname,ires,cbi+6)
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD1',resname,ires,cbi+2)
         call pdbatm (' CD2',resname,ires,cbi+3)
         call pdbatm (' CE1',resname,ires,cbi+4)
         call pdbatm (' CE2',resname,ires,cbi+5)
         call pdbatm (' CZ ',resname,ires,cbi+6)
         call pdbatm (' OH ',resname,ires,cbi+7)
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD1',resname,ires,cbi+2)
         call pdbatm (' CD2',resname,ires,cbi+3)
         call pdbatm (' NE1',resname,ires,cbi+4)
         call pdbatm (' CE2',resname,ires,cbi+5)
         call pdbatm (' CE3',resname,ires,cbi+6)
         call pdbatm (' CZ2',resname,ires,cbi+7)
         call pdbatm (' CZ3',resname,ires,cbi+8)
         call pdbatm (' CH2',resname,ires,cbi+9)
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' ND1',resname,ires,cbi+2)
         call pdbatm (' CD2',resname,ires,cbi+3)
         call pdbatm (' CE1',resname,ires,cbi+4)
         call pdbatm (' NE2',resname,ires,cbi+5)
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' ND1',resname,ires,cbi+2)
         call pdbatm (' CD2',resname,ires,cbi+3)
         call pdbatm (' CE1',resname,ires,cbi+4)
         call pdbatm (' NE2',resname,ires,cbi+5)
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' ND1',resname,ires,cbi+2)
         call pdbatm (' CD2',resname,ires,cbi+3)
         call pdbatm (' CE1',resname,ires,cbi+4)
         call pdbatm (' NE2',resname,ires,cbi+5)
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' OD1',resname,ires,cbi+2)
         call pdbatm (' OD2',resname,ires,cbi+3)
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' OD1',resname,ires,cbi+2)
         call pdbatm (' ND2',resname,ires,cbi+3)
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD ',resname,ires,cbi+2)
         call pdbatm (' OE1',resname,ires,cbi+3)
         call pdbatm (' OE2',resname,ires,cbi+4)
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD ',resname,ires,cbi+2)
         call pdbatm (' OE1',resname,ires,cbi+3)
         call pdbatm (' NE2',resname,ires,cbi+4)
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' SD ',resname,ires,cbi+2)
         call pdbatm (' CE ',resname,ires,cbi+3)
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD ',resname,ires,cbi+2)
         call pdbatm (' CE ',resname,ires,cbi+3)
         call pdbatm (' NZ ',resname,ires,cbi+4)
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD ',resname,ires,cbi+2)
         call pdbatm (' NE ',resname,ires,cbi+3)
         call pdbatm (' CZ ',resname,ires,cbi+4)
         call pdbatm (' NH1',resname,ires,cbi+5)
         call pdbatm (' NH2',resname,ires,cbi+6)
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD ',resname,ires,cbi+2)
         call pdbatm (' NE ',resname,ires,cbi+3)
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         call pdbatm (' CB2',resname,ires,cbi+1)
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         call pdbatm (' CG ',resname,ires,cbi+1)
         call pdbatm (' CD ',resname,ires,cbi+2)
         call pdbatm (' OE ',resname,ires,cbi+3)
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         continue
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine gethydro  --  extract hydrogen atoms  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gethydro" finds the side chain hydrogen atoms for
c     a single amino acid residue and copies the names
c     and coordinates to the Protein Data Bank file
c
c
      subroutine gethydro (resname,ires,jchain,ni,cai,cbi)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'sequen.i'
      integer i,nh,hca
      integer ires,jchain
      integer ni,cai,cbi
      character*3 resname
      character*4 atmname
      logical allatom
c
c
c     get any amide hydrogen atoms for non-N-terminal residues
c
      if (ires .ne. ichain(1,jchain)) then
         if (resname .ne. 'PRO') then
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  if (resname .eq. 'NH2') then
                     call pdbatm ('1H  ',resname,ires,i)
                     call pdbatm ('2H  ',resname,ires,i+1)
                  else
                     call pdbatm (' H  ',resname,ires,i)
                  end if
                  goto 10
               end if
            end do
         end if
c
c     get any amide hydrogen atoms for N-terminal residues
c
      else
         if (resname .eq. 'PRO') then
            nh = 0
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  nh = nh + 1
                  if (nh .eq. 1) then
                     atmname = '1H  '
                  else if (nh .eq. 2) then
                     atmname = '2H  '
                  end if
                  call pdbatm (atmname,resname,ires,i)
                  if (nh .eq. 2)  goto 10
               end if
            end do
         else if (resname .eq. 'PCA') then
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  atmname = ' H  '
                  call pdbatm (atmname,resname,ires,i)
                  goto 10
               end if
            end do
         else
            nh = 0
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  nh = nh + 1
                  if (nh .eq. 1) then
                     atmname = '1H  '
                  else if (nh .eq. 2) then
                     atmname = '2H  '
                  else if (nh .eq. 3) then
                     atmname = '3H  '
                  end if
                  call pdbatm (atmname,resname,ires,i)
                  if (nh .eq. 3)  goto 10
               end if
            end do
         end if
      end if
   10 continue
c
c     get the alpha hydrogen atom for the current residue
c
      hca = 0
      do i = 1, n
         if (atomic(i).eq.1 .and. i12(1,i).eq.cai) then
            hca = i
            if (resname .eq. 'GLY') then
               atmname = '1HA '
            else if (resname .eq. 'FOR') then
               atmname = ' H  '
            else if (resname .eq. 'ACE') then
               atmname = '1H  '
            else if (resname .eq. 'NME') then
               atmname = '1H  '
            else
               atmname = ' HA '
            end if
            call pdbatm (atmname,resname,ires,i)
            goto 20
         end if
      end do
   20 continue
c
c     if no alpha hydrogen, then united atom force field
c
      if (hca .ne. 0) then
         allatom = .true.
      else if (resname .eq. 'AIB') then
         if (n12(cbi) .eq. 1) then
            allatom = .false.
         else
            allatom = .true.
         end if
      else
         allatom = .false.
      end if
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         if (allatom) then
            call pdbatm ('2HA ',resname,ires,hca+1)
         end if
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+2)
            call pdbatm ('2HB ',resname,ires,hca+3)
            call pdbatm ('3HB ',resname,ires,hca+4)
         end if
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         if (allatom) then
            call pdbatm (' HB ',resname,ires,hca+4)
            call pdbatm ('1HG1',resname,ires,hca+5)
            call pdbatm ('2HG1',resname,ires,hca+6)
            call pdbatm ('3HG1',resname,ires,hca+7)
            call pdbatm ('1HG2',resname,ires,hca+8)
            call pdbatm ('2HG2',resname,ires,hca+9)
            call pdbatm ('3HG2',resname,ires,hca+10)
         end if
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+5)
            call pdbatm ('2HB ',resname,ires,hca+6)
            call pdbatm (' HG ',resname,ires,hca+7)
            call pdbatm ('1HD1',resname,ires,hca+8)
            call pdbatm ('2HD1',resname,ires,hca+9)
            call pdbatm ('3HD1',resname,ires,hca+10)
            call pdbatm ('1HD2',resname,ires,hca+11)
            call pdbatm ('2HD2',resname,ires,hca+12)
            call pdbatm ('3HD2',resname,ires,hca+13)
         end if
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         if (allatom) then
            call pdbatm (' HB ',resname,ires,hca+5)
            call pdbatm ('1HG1',resname,ires,hca+6)
            call pdbatm ('2HG1',resname,ires,hca+7)
            call pdbatm ('1HG2',resname,ires,hca+8)
            call pdbatm ('2HG2',resname,ires,hca+9)
            call pdbatm ('3HG2',resname,ires,hca+10)
            call pdbatm ('1HD1',resname,ires,hca+11)
            call pdbatm ('2HD1',resname,ires,hca+12)
            call pdbatm ('3HD1',resname,ires,hca+13)
         end if
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+3)
            call pdbatm ('2HB ',resname,ires,hca+4)
            call pdbatm (' HG ',resname,ires,hca+5)
         else
            call pdbatm (' HG ',resname,ires,cbi+2)
         end if
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         if (allatom) then
            call pdbatm (' HB ',resname,ires,hca+4)
            call pdbatm ('1HG ',resname,ires,hca+5)
            call pdbatm ('1HG2',resname,ires,hca+6)
            call pdbatm ('2HG2',resname,ires,hca+7)
            call pdbatm ('3HG2',resname,ires,hca+8)
         else
            call pdbatm ('1HG ',resname,ires,cbi+3)
         end if
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+3)
            call pdbatm ('2HB ',resname,ires,hca+4)
            call pdbatm (' HG ',resname,ires,hca+5)
         else
            call pdbatm (' HG ',resname,ires,cbi+2)
         end if
c
c     cystine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+3)
            call pdbatm ('2HB ',resname,ires,hca+4)
         end if
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+4)
            call pdbatm ('2HB ',resname,ires,hca+5)
            call pdbatm ('1HG ',resname,ires,hca+6)
            call pdbatm ('2HG ',resname,ires,hca+7)
            call pdbatm ('1HD ',resname,ires,hca+8)
            call pdbatm ('2HD ',resname,ires,hca+9)
         end if
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+8)
            call pdbatm ('2HB ',resname,ires,hca+9)
            call pdbatm (' HD1',resname,ires,hca+10)
            call pdbatm (' HD2',resname,ires,hca+11)
            call pdbatm (' HE1',resname,ires,hca+12)
            call pdbatm (' HE2',resname,ires,hca+13)
            call pdbatm (' HZ ',resname,ires,hca+14)
         end if
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+9)
            call pdbatm ('2HB ',resname,ires,hca+10)
            call pdbatm (' HD1',resname,ires,hca+11)
            call pdbatm (' HD2',resname,ires,hca+12)
            call pdbatm (' HE1',resname,ires,hca+13)
            call pdbatm (' HE2',resname,ires,hca+14)
            call pdbatm (' HH ',resname,ires,hca+15)
         else
            call pdbatm (' HH ',resname,ires,cbi+8)
         end if
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+11)
            call pdbatm ('2HB ',resname,ires,hca+12)
            call pdbatm (' HD1',resname,ires,hca+13)
            call pdbatm (' HE1',resname,ires,hca+14)
            call pdbatm (' HE3',resname,ires,hca+15)
            call pdbatm (' HZ2',resname,ires,hca+16)
            call pdbatm (' HZ3',resname,ires,hca+17)
            call pdbatm (' HH2',resname,ires,hca+18)
         else
            call pdbatm (' HE1',resname,ires,cbi+10)
         end if
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+7)
            call pdbatm ('2HB ',resname,ires,hca+8)
            call pdbatm (' HD1',resname,ires,hca+9)
            call pdbatm (' HD2',resname,ires,hca+10)
            call pdbatm (' HE1',resname,ires,hca+11)
            call pdbatm (' HE2',resname,ires,hca+12)
         else
            call pdbatm (' HD1',resname,ires,cbi+6)
            call pdbatm (' HE2',resname,ires,cbi+7)
         end if
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+7)
            call pdbatm ('2HB ',resname,ires,hca+8)
            call pdbatm (' HD1',resname,ires,hca+9)
            call pdbatm (' HD2',resname,ires,hca+10)
            call pdbatm (' HE1',resname,ires,hca+11)
         else
            call pdbatm (' HD1',resname,ires,cbi+6)
         end if
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+7)
            call pdbatm ('2HB ',resname,ires,hca+8)
            call pdbatm (' HD2',resname,ires,hca+9)
            call pdbatm (' HE1',resname,ires,hca+10)
            call pdbatm (' HE2',resname,ires,hca+11)
         else
            call pdbatm (' HE2',resname,ires,cbi+6)
         end if
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+5)
            call pdbatm ('2HB ',resname,ires,hca+6)
         end if
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+5)
            call pdbatm ('2HB ',resname,ires,hca+6)
            call pdbatm ('1HD2',resname,ires,hca+7)
            call pdbatm ('2HD2',resname,ires,hca+8)
         else
            call pdbatm ('1HD2',resname,ires,cbi+4)
            call pdbatm ('2HD2',resname,ires,cbi+5)
         end if
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+6)
            call pdbatm ('2HB ',resname,ires,hca+7)
            call pdbatm ('1HG ',resname,ires,hca+8)
            call pdbatm ('2HG ',resname,ires,hca+9)
         end if
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+6)
            call pdbatm ('2HB ',resname,ires,hca+7)
            call pdbatm ('1HG ',resname,ires,hca+8)
            call pdbatm ('2HG ',resname,ires,hca+9)
            call pdbatm ('1HE2',resname,ires,hca+10)
            call pdbatm ('2HE2',resname,ires,hca+11)
         else
            call pdbatm ('1HE2',resname,ires,cbi+5)
            call pdbatm ('2HE2',resname,ires,cbi+6)
         end if
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+5)
            call pdbatm ('2HB ',resname,ires,hca+6)
            call pdbatm ('1HG ',resname,ires,hca+7)
            call pdbatm ('2HG ',resname,ires,hca+8)
            call pdbatm ('1HE ',resname,ires,hca+9)
            call pdbatm ('2HE ',resname,ires,hca+10)
            call pdbatm ('3HE ',resname,ires,hca+11)
         end if
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+6)
            call pdbatm ('2HB ',resname,ires,hca+7)
            call pdbatm ('1HG ',resname,ires,hca+8)
            call pdbatm ('2HG ',resname,ires,hca+9)
            call pdbatm ('1HD ',resname,ires,hca+10)
            call pdbatm ('2HD ',resname,ires,hca+11)
            call pdbatm ('1HE ',resname,ires,hca+12)
            call pdbatm ('2HE ',resname,ires,hca+13)
            call pdbatm ('1HZ ',resname,ires,hca+14)
            call pdbatm ('2HZ ',resname,ires,hca+15)
            call pdbatm ('3HZ ',resname,ires,hca+16)
         else
            call pdbatm ('1HZ ',resname,ires,cbi+5)
            call pdbatm ('2HZ ',resname,ires,cbi+6)
            call pdbatm ('3HZ ',resname,ires,cbi+7)
         end if
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+8)
            call pdbatm ('2HB ',resname,ires,hca+9)
            call pdbatm ('1HG ',resname,ires,hca+10)
            call pdbatm ('2HG ',resname,ires,hca+11)
            call pdbatm ('1HD ',resname,ires,hca+12)
            call pdbatm ('2HD ',resname,ires,hca+13)
            call pdbatm (' HE ',resname,ires,hca+14)
            call pdbatm ('1HH1',resname,ires,hca+15)
            call pdbatm ('2HH1',resname,ires,hca+16)
            call pdbatm ('1HH2',resname,ires,hca+17)
            call pdbatm ('2HH2',resname,ires,hca+18)
         else
            call pdbatm (' HE ',resname,ires,cbi+7)
            call pdbatm ('1HH1',resname,ires,cbi+8)
            call pdbatm ('2HH1',resname,ires,cbi+9)
            call pdbatm ('1HH2',resname,ires,cbi+10)
            call pdbatm ('2HH2',resname,ires,cbi+11)
         end if
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+5)
            call pdbatm ('2HB ',resname,ires,hca+6)
            call pdbatm ('1HG ',resname,ires,hca+7)
            call pdbatm ('2HG ',resname,ires,hca+8)
            call pdbatm ('1HD ',resname,ires,hca+9)
            call pdbatm ('2HD ',resname,ires,hca+10)
            call pdbatm ('1HE ',resname,ires,hca+11)
            call pdbatm ('2HE ',resname,ires,hca+12)
            call pdbatm ('1HE ',resname,ires,hca+13)
         else
            call pdbatm ('1HE ',resname,ires,cbi+4)
            call pdbatm ('2HE ',resname,ires,cbi+5)
            call pdbatm ('3HE ',resname,ires,cbi+6)
         end if
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         if (allatom) then
            call pdbatm ('1HB1',resname,ires,cbi+2)
            call pdbatm ('2HB1',resname,ires,cbi+3)
            call pdbatm ('3HB1',resname,ires,cbi+4)
            call pdbatm ('1HB2',resname,ires,cbi+5)
            call pdbatm ('2HB2',resname,ires,cbi+6)
            call pdbatm ('3HB2',resname,ires,cbi+7)
         end if
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         if (allatom) then
            call pdbatm ('1HB ',resname,ires,hca+5)
            call pdbatm ('2HB ',resname,ires,hca+6)
            call pdbatm ('1HG ',resname,ires,hca+7)
            call pdbatm ('2HG ',resname,ires,hca+8)
         end if
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         if (allatom) then
            call pdbatm ('2HA ',resname,ires,hca+1)
         end if
c
c     N-terminal acetyl residue  (ACE)
c
      else if (resname .eq. 'ACE') then
         if (allatom) then
            call pdbatm ('2H  ',resname,ires,hca+1)
            call pdbatm ('3H  ',resname,ires,hca+2)
         end if
c
c     N-terminal formyl residue  (FOR)
c
      else if (resname .eq. 'FOR') then
         continue
c
c     C-terminal N-methylamide residue  (NME)
c
      else if (resname .eq. 'NME') then
         if (allatom) then
            call pdbatm ('2H  ',resname,ires,hca+1)
            call pdbatm ('3H  ',resname,ires,hca+2)
         end if
c
c     C-terminal amide residue  (NH2)
c
      else if (resname .eq. 'NH2') then
         continue
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine makeref  --  copy structure to reference area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "makeref" copies the information contained in the "xyz" file
c     of the current structure into corresponding reference areas
c
c
      subroutine makeref
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'refer.i'
      include 'titles.i'
      integer i,j
c
c
c     copy the filename and title line for the structure
c
      reffile = filename
      refleng = leng
      reftitle = ttitle
      refltitle = ltitle
c
c     copy the coordinates, type and connectivity of each atom
c
      nref = n
      do i = 1, n
         refnam(i) = name(i)
         xref(i) = x(i)
         yref(i) = y(i)
         zref(i) = z(i)
         reftyp(i) = type(i)
         n12ref(i) = n12(i)
         do j = 1, n12(i)
            i12ref(j,i) = i12(j,i)
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine makexyz  --  convert internal to Cartesian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "makexyz" generates a complete set of Cartesian coordinates
c     for a full structure from the internal coordinate values
c
c
      subroutine makexyz
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'zcoord.i'
      integer i,ia,ib,ic,chiral
      real*8 bond,angle1,angle2
c
c
c     loop over each atom in turn, finding its coordinates
c
      do i = 1, n
         ia = iz(1,i)
         ib = iz(2,i)
         ic = iz(3,i)
         chiral = iz(4,i)
         bond = zbond(i)
         angle1 = zang(i)
         angle2 = ztors(i)
         call xyzatm (i,ia,bond,ib,angle1,ic,angle2,chiral)
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  function maxwell  --  Maxwell-Boltzmann distribution value  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "maxwell" returns a speed in Ang/sec randomly selected from a
c     3-D Maxwell-Boltzmann distribution for the specified particle
c     mass and system temperature
c
c     literature reference:
c
c     P. W. Atkins, "Physical Chemistry, 4th Edition", W. H. Freeman,
c     New York, 1990; see section 24.2 for general discussion
c
c
      function maxwell (mass,temper)
      implicit none
      include 'units.i'
      real*8 maxwell,mass,temper
      real*8 rho,beta,random,erfinv
      real*8 xspeed,yspeed,zspeed
c
c
c     set normalization factor for cumulative velocity distribution
c
      beta = sqrt(mass / (2.0d0*boltzmann*temper))
c
c     pick a randomly distributed velocity along each of three axes
c
      rho = random ()
      xspeed = erfinv (rho) / beta
      rho = random ()
      yspeed = erfinv (rho) / beta
      rho = random ()
      zspeed = erfinv (rho) / beta
c
c     set the final value of the particle speed in 3-dimensions
c
      maxwell = sqrt(xspeed**2 + yspeed**2 + zspeed**2)
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine mdinit  --  initialize an MD trajectory  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "mdinit" initializes the velocities and accelerations
c     for a molecular dynamics trajectory, including restarts
c
c
      subroutine mdinit
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'files.i'
      include 'iounit.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,idyn
      integer lext,freeunit
      real*8 e,maxwell,vel
      real*8 vec(3)
      real*8 derivs(3,maxatm)
      character*7 ext
      character*60 dynfile
      logical exist
c
c
c     try to restart using prior velocities and accelerations
c
      dynfile = filename(1:leng)//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      if (exist) then
         idyn = freeunit ()
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
         call readdyn (idyn)
         close (unit=idyn)
         call lattice
c
c     initialize velocities and accelerations for first step
c
      else
         call gradient (e,derivs)
         do i = 1, n
            if (use(i)) then
               vel = maxwell (mass(i),kelvin)
               call ranvec (vec)
               do j = 1, 3
                  v(j,i) = vel * vec(j)
                  a(j,i) = -convert * derivs(j,i) / mass(i)
                  a_old(j,i) = a(j,i)
               end do
            else
               do j = 1, 3
                  v(j,i) = 0.0d0
                  a(j,i) = 0.0d0
                  a_old(j,i) = 0.0d0
               end do
            end if
         end do
         if (nuse .eq. n)  call mdrest
      end if
c
c     check for any prior dynamics coordinate sets
c
      i = 0
      exist = .true.
      dowhile (exist)
         i = i + 1
         lext = 3
         call numeral (i,ext,lext)
         dynfile = filename(1:leng)//'.'//ext(1:lext)
         inquire (file=dynfile,exist=exist)
         if (.not.exist .and. i.lt.100) then
            lext = 2
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
         if (.not.exist .and. i.lt.10) then
            lext = 1
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
      end do
      nprior = i - 1
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mdrest  --  stop translation & rotation of CM  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mdrest" computes and then removes any translational
c     or rotational kinetic energy of the center of mass
c
c
      subroutine mdrest
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moldyn.i'
      include 'units.i'
      integer i,j
      real*8 etrans,erot
      real*8 weight,total
      real*8 xcm,ycm,zcm
      real*8 xx,xy,xz,yy,yz,zz
      real*8 xdiff,ydiff,zdiff
      real*8 mang(3),vang(3)
      real*8 vcm(3),tensor(3,3)
c
c
c     compute the linear velocity of the center of mass
c
      total = 0.0d0
      do j = 1, 3
         vcm(j) = 0.0d0
      end do
      do i = 1, n
         weight = mass(i)
         total = total + weight
         do j = 1, 3
            vcm(j) = vcm(j) + v(j,i)*weight
         end do
      end do
c
c     compute translational kinetic energy of center of mass
c
      etrans = 0.0d0
      do j = 1, 3
         vcm(j) = vcm(j) / total
         etrans = etrans + vcm(j)**2
      end do
      etrans = 0.5d0 * etrans * total / convert
c
c     eliminate any translation of the center of mass
c
      do i = 1, n
         do j = 1, 3
            v(j,i) = v(j,i) - vcm(j)
         end do
      end do
c
c     print the linear velocity of the center of mass
c
      if (verbose) then
         write (iout,10)  (vcm(i),i=1,3),etrans
   10    format (/,' Linear Velocity of CM : ',3d12.2,
     &           /,' Translational Kinetic Energy :',f20.4,
     &              ' Kcal/mole')
      end if
c
c     compute the angular momentum about the center of mass
c
      if (.not. use_bounds) then
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         do j = 1, 3
            mang(j) = 0.0d0
         end do
         do i = 1, n
            weight = mass(i)
            xcm = xcm + x(i)*weight
            ycm = ycm + y(i)*weight
            zcm = zcm + z(i)*weight
            mang(1) = mang(1) + (y(i)*v(3,i)-z(i)*v(2,i))*weight
            mang(2) = mang(2) + (z(i)*v(1,i)-x(i)*v(3,i))*weight
            mang(3) = mang(3) + (x(i)*v(2,i)-y(i)*v(1,i))*weight
         end do
         xcm = xcm / total
         ycm = ycm / total
         zcm = zcm / total
         mang(1) = mang(1) - (ycm*vcm(3)-zcm*vcm(2))*total
         mang(2) = mang(2) - (zcm*vcm(1)-xcm*vcm(3))*total
         mang(3) = mang(3) - (xcm*vcm(2)-ycm*vcm(1))*total
c
c     calculate and then invert the inertial tensor
c
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         do i = 1, n
            weight = mass(i)
            xdiff = x(i) - xcm
            ydiff = y(i) - ycm
            zdiff = z(i) - zcm
            xx = xx + xdiff*xdiff*weight
            xy = xy + xdiff*ydiff*weight
            xz = xz + xdiff*zdiff*weight
            yy = yy + ydiff*ydiff*weight
            yz = yz + ydiff*zdiff*weight
            zz = zz + zdiff*zdiff*weight
         end do
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
         call invert (3,3,tensor)
c
c     compute angular velocity and rotational kinetic energy
c
         erot = 0.0d0
         do i = 1, 3
            vang(i) = 0.0d0
            do j = 1, 3
               vang(i) = vang(i) + tensor(i,j)*mang(j)
            end do
            erot = erot + vang(i)*mang(i)
         end do
         erot = 0.5d0 * erot / convert
c
c     eliminate any rotation about the center of mass
c
         do i = 1, n
            xdiff = x(i) - xcm
            ydiff = y(i) - ycm
            zdiff = z(i) - zcm
            v(1,i) = v(1,i) - vang(2)*zdiff + vang(3)*ydiff
            v(2,i) = v(2,i) - vang(3)*xdiff + vang(1)*zdiff
            v(3,i) = v(3,i) - vang(1)*ydiff + vang(2)*xdiff
         end do
c
c     print the angular velocity about the center of mass
c
         if (verbose) then
            write (iout,20)  (vang(i),i=1,3),erot
   20       format (/,' Angular Velocity of CM :',3d12.2,
     &              /,' Rotational Kinetic Energy :',3x,f20.4,
     &                 ' Kcal/mole')
         end if
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdstat  --  compute averages for MD trajectory  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdstat" is called at each molecular dynamics time step to
c     form statistics on various average values and fluctuations,
c     and to periodically save the state of the trajectory
c
c
      subroutine mdstat (istep,dt,e_tot,e_pot,e_kin,temp,pres,vol,ndump)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'files.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'units.i'
      include 'usage.i'
      include 'warp.i'
      integer istep,ixyz,iend
      integer idump,ndump,lext,freeunit
      integer period,modstep,moddump
      real*8 dt,e_tot,e_pot,e_kin
      real*8 temp,pres,vol,pico,dens
      real*8 fluctuate,fluctuate2,intfluct,intfluct2
      real*8 potfluct,potfluct2,kinfluct,kinfluct2
      real*8 tfluct,tfluct2,pfluct,pfluct2,dfluct,dfluct2
      real*8 e_tot_sum,e_tot2_sum,e_int_sum,e_int2_sum
      real*8 e_tot_ave,e_tot2_ave,e_int_ave,e_int2_ave
      real*8 e_pot_sum,e_pot2_sum,e_kin_sum,e_kin2_sum
      real*8 e_pot_ave,e_pot2_ave,e_kin_ave,e_kin2_ave
      real*8 temp_sum,temp2_sum,temp_ave,temp2_ave
      real*8 pres_sum,pres2_sum,pres_ave,pres2_ave
      real*8 dens_sum,dens2_sum,dens_ave,dens2_ave
      character*7 ext
      character*60 coordfile,endfile
      logical exist
      save e_tot_sum,e_tot2_sum,e_int_sum,e_int2_sum
      save e_pot_sum,e_pot2_sum,e_kin_sum,e_kin2_sum
      save temp_sum,temp2_sum,pres_sum,pres2_sum
      save dens_sum,dens2_sum
c
c
c     write energy, temperature and pressure for current step
c
      if (verbose) then
         if (isobaric) then
            write (iout,10)  istep,e_tot,e_pot,e_kin,temp,pres
   10       format (i10,2f14.4,f12.4,2f11.2)
         else
            write (iout,20)  istep,e_tot,e_pot,e_kin,temp
   20       format (i10,2f14.4,f12.4,f11.2)
         end if
      end if
c
c     set number of steps for average and coordinate output
c
      period = 100
      modstep = mod(istep,period)
      moddump = mod(istep,ndump)
c
c     zero out summation variables for new averaging period
c
      if (modstep .eq. 1) then
         e_tot_sum = 0.0d0
         e_tot2_sum = 0.0d0
         e_pot_sum = 0.0d0
         e_pot2_sum = 0.0d0
         e_kin_sum = 0.0d0
         e_kin2_sum = 0.0d0
         e_int_sum = 0.0d0
         e_int2_sum = 0.0d0
         temp_sum = 0.0d0
         temp2_sum = 0.0d0
         pres_sum = 0.0d0
         pres2_sum = 0.0d0
         dens_sum = 0.0d0
         dens2_sum = 0.0d0
      end if
c
c     print header for the averages over a group of recent steps
c
      if (modstep .eq. 0) then
         pico = 1.0d+12 * dt * dble(istep)
         write (iout,30)  period,istep,pico
   30    format (/,' Values over the last ',i3,' out of',i8,
     &              ' Dynamics Steps',
     &           //,' Cumulative Time :',1x,f15.4,' Picosecond')
      end if
c
c     compute total energy and fluctuation for recent steps
c
      e_tot_sum = e_tot_sum + e_tot
      e_tot2_sum = e_tot2_sum + e_tot**2
      if (modstep .eq. 0) then
         e_tot_ave = e_tot_sum / dble(period)
         e_tot2_ave = e_tot2_sum / dble(period)
         fluctuate2 = e_tot2_ave - e_tot_ave**2
         if (fluctuate2 .gt. 0.0d0) then
            fluctuate = sqrt(fluctuate2)
         else
            fluctuate = 0.0d0
         end if
         write (iout,40)  e_tot_ave,fluctuate
   40    format (' Total Energy :',4x,f15.4,' Kcal/mole',3x,
     &              '(+/-',f9.4,')')
      end if
c
c     compute average potential energy and its fluctuation
c
      e_pot_sum = e_pot_sum + e_pot
      e_pot2_sum = e_pot2_sum + e_pot**2
      if (modstep .eq. 0) then
         e_pot_ave = e_pot_sum / dble(period)
         e_pot2_ave = e_pot2_sum / dble(period)
         potfluct2 = e_pot2_ave - e_pot_ave**2
         if (potfluct2 .gt. 0.0d0) then
            potfluct = sqrt(potfluct2)
         else
            potfluct = 0.0d0
         end if
         write (iout,50)  e_pot_ave,potfluct
   50    format (' Potential Energy :',f15.4,' Kcal/mole',3x,
     &              '(+/-',f9.4,')')
      end if
c
c     compute average kinetic energy and its fluctuation
c
      e_kin_sum = e_kin_sum + e_kin
      e_kin2_sum = e_kin2_sum + e_kin**2
      if (modstep .eq. 0) then
         e_kin_ave = e_kin_sum / dble(period)
         e_kin2_ave = e_kin2_sum / dble(period)
         kinfluct2 = e_kin2_ave - e_kin_ave**2
         if (kinfluct2 .gt. 0.0d0) then
            kinfluct = sqrt(kinfluct2)
         else
            kinfluct = 0.0d0
         end if
         write (iout,60)  e_kin_ave,kinfluct
   60    format (' Kinetic Energy :',2x,f15.4,' Kcal/mole',3x,
     &              '(+/-',f9.4,')')
      end if
c
c     compute average intermolecular energy and its fluctuation
c
      if (nmol.ne.1 .and. nmol.ne.n) then
         e_int_sum = e_int_sum + einter
         e_int2_sum = e_int2_sum + einter**2
         if (modstep .eq. 0) then
            e_int_ave = e_int_sum / dble(period)
            e_int2_ave = e_int2_sum / dble(period)
            intfluct2 = e_int2_ave - e_int_ave**2
            if (intfluct2 .gt. 0.0d0) then
               intfluct = sqrt(intfluct2)
            else
               intfluct = 0.0d0
            end if
            write (iout,70)  e_int_ave,intfluct
   70       format (' Intermolecular :',2x,f15.4,' Kcal/mole',3x,
     &                 '(+/-',f9.4,')')
         end if
      end if
c
c     compute the average temperature and its fluctuation
c
      temp_sum = temp_sum + temp
      temp2_sum = temp2_sum + temp**2
      if (modstep .eq. 0) then
         temp_ave = temp_sum / dble(period)
         temp2_ave = temp2_sum / dble(period)
         tfluct2 = temp2_ave - temp_ave**2
         if (tfluct2 .gt. 0.0d0) then
            tfluct = sqrt(tfluct2)
         else
            tfluct = 0.0d0
         end if
         write (iout,80)  temp_ave,tfluct
   80    format (' Temperature :',5x,f15.2,' Kelvin',6x,
     &              '(+/-',f9.2,')')
      end if
c
c     compute the average pressure and its fluctuation
c
      if (isobaric) then
         pres_sum = pres_sum + pres
         pres2_sum = pres2_sum + pres**2
         if (modstep .eq. 0) then
            pres_ave = pres_sum / dble(period)
            pres2_ave = pres2_sum / dble(period)
            pfluct2 = pres2_ave - pres_ave**2
            if (pfluct2 .gt. 0.0d0) then
               pfluct = sqrt(pfluct2)
            else
               pfluct = 0.0d0
            end if
            write (iout,90)  pres_ave,pfluct
   90       format (' Pressure :',8x,f15.2,' Atmosphere',2x,
     &                 '(+/-',f9.2,')')
         end if
c
c     compute the average density and its fluctuation
c
         dens = (1.0d24/avogadro) * (totmass/vol)
         dens_sum = dens_sum + dens
         dens2_sum = dens2_sum + dens**2
         if (modstep .eq. 0) then
            dens_ave = dens_sum / dble(period)
            dens2_ave = dens2_sum / dble(period)
            dfluct2 = dens2_ave - dens_ave**2
            if (dfluct2 .gt. 0.0d0) then
               dfluct = sqrt(dfluct2)
            else
               dfluct = 0.0d0
            end if
            write (iout,100)  dens_ave,dfluct
  100       format (' Density :',9x,f15.4,' Grams/cc',4x,
     &                 '(+/-',f9.4,')')
         end if
      end if
c
c     current deformation value for diffusional smoothing
c
      if (use_deform) then
         if (modstep .eq. 0) then
            write (iout,110)  deform
  110       format (' Deformation :',5x,f15.3,' Sqr Angs')
         end if
      end if
c
c     save a structure along the trajectory every so often
c
      if (moddump .eq. 0) then
         idump = nprior + istep/ndump
         lext = 3
         call numeral (idump,ext,lext)
         ixyz = freeunit ()
         coordfile = filename(1:leng)//'.'//ext(1:lext)
         call version (coordfile,'new')
         open (unit=ixyz,file=coordfile,status='new')
         call prtxyz (ixyz)
         close (unit=ixyz)
c
c     save information needed for a restart to a file
c
         call prtdyn
c
c     check for requested termination of the dynamics run
c
         iend = freeunit ()
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            write (iout,120)
  120       format (/,' Molecular Dynamics run Terminating due to',
     &                 ' User Request')
            call fatal
         end if
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic
      implicit double precision (a-h,o-z)
cnone
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      include 'vdwpot.i'
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER me,master,nproc,ibtyp,iptim,IR,IW
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
c
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c
c     set the connectivity, active atoms, and atom groups
c
      call attach
      call active
      call cluster
c
c     find and store each bond, angle and torsion
c
      call bonds
      call angles
      call torsions
c
c     find cutoff values and small ring systems
c
      call cutoffs
      call rings
c
c     get the force field parameters and assign atom types
c
      call field
      call katom
c
c     find any pisystem atoms, bonds and torsional angles
c
      call orbital
c
c     assign local geometry potential function parameters
c
      if (use_bond .or. use_strbnd .or. use_strtor
     &    .or. (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
      if (use_angle .or. use_strbnd .or. use_angang)  call kangle
      if (use_strbnd)  call kstrbnd
      if (use_urey)  call kurey
      if (use_angang)  call kangang
      if (use_angle .or. use_opbend)  call kopbend
      if (use_improp)  call kimprop
      if (use_imptor)  call kimptor
      if (use_tors .or. use_strtor .or. use_tortor)  call ktors
      if (use_strtor)  call kstrtor
c     if (use_tortor)  call ktortor
c
c     assign nonbonded potential function parameters
c
      if (use_vdw)  call kvdw
      if (use_charge .or. use_chgdpl .or. use_solv)  call kcharge
      if (use_dipole .or. use_chgdpl)  call kdipole
      if (use_mpole .or. use_polar .or. use_rxnfld)  call kmpole
      if (use_polar)  call kpolar
c
c     assign macroscopic solvation and pisystem parameters
c
      if (use_solv)  call solvate
      if (use_orbit)  call korbit
c
c     count the number of molcules in the system
c
      call molecule
c
c     find any crystal lattice or periodic box parameters
c
      call unitcell
      call lattice
c
c     set any geometric restraint terms to be applied
c
      call restrain
c
c     set any deformation parameters for potential smoothing
c
      call smooth
c
c     set hybrid parameters for free energy perturbation
c
      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         if (maswrk) write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine merge  --  merge reference & current systems  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "merge" combines the reference and current structures into
c     a single new "current" structure containing the reference
c     atoms followed by the atoms of the current structure
c
c
      subroutine merge
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'refer.i'
      integer i,j,k,ntotal
c
c
c     check for too many total atoms in the combined system
c
      ntotal = n + nref
      if (ntotal .gt. maxatm) then
         write (iout,10)  maxatm
   10    format (/,' MERGE  --  The Maximum of',i6,' Atoms',
     &              ' has been Exceeded')
         call fatal
      end if
c
c     move the current structure to higher atom numbers
c
      do i = n, 1, -1
         k = i + nref
         x(k) = x(i)
         y(k) = y(i)
         z(k) = z(i)
         type(k) = type(i)
         name(k) = name(i)
         n12(k) = n12(i)
         do j = 1, n12(i)
            i12(j,k) = i12(j,i) + nref
         end do
      end do
c
c     place reference structure in the current structure
c
      call getref
      n = ntotal
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine molecule  --  assign atoms to molecules  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "molecule" counts the molecules, assigns each atom to
c     its molecule and computes the mass of each molecule
c
c
      subroutine molecule
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'molcul.i'
      integer i,j,k,mi,mj,mk
      integer iattach,list(maxatm)
c
c
c     zero number of molecules and molecule membership list
c
      nmol = 0
      do i = 1, n
         molcule(i) = 0
      end do
c
c     assign each atom to its respective molecule
c
      do i = 1, n
         if (molcule(i) .eq. 0) then
            nmol = nmol + 1
            molcule(i) = nmol
         end if
         mi = molcule(i)
         do iattach = 1, n12(i)
            j = i12(iattach,i)
            mj = molcule(j)
            if (mj .eq. 0) then
               molcule(j) = mi
            else if (mi .lt. mj) then
               nmol = nmol - 1
               do k = 1, n
                  mk = molcule(k)
                  if (mk .eq. mj) then
                     molcule(k) = mi
                  else if (mk .gt. mj) then
                     molcule(k) = mk - 1
                  end if
               end do
            else if (mi .gt. mj) then
               nmol = nmol - 1
               do k = 1, n
                  mk = molcule(k)
                  if (mk .eq. mi) then
                     molcule(k) = mj
                  else if (mk .gt. mi) then
                     molcule(k) = mk - 1
                  end if
               end do
               mi = mj
            end if
         end do
      end do
c
c     pack atoms of each molecule into a contiguous indexed list
c
      do i = 1, n
         list(i) = molcule(i)
      end do
      call sort3 (n,list,kmol)
c
c     find the first and last atom in each molecule
c
      k = 1
      imol(1,1) = 1
      do i = 2, n
         j = list(i)
         if (j .ne. k) then
            imol(2,k) = i - 1
            k = j
            imol(1,k) = i
         end if
      end do
      imol(2,nmol) = n
c
c     sort the list of atoms in each molecule by atom number
c
      do i = 1, nmol
         k = imol(2,i) - imol(1,i) + 1
         call sort (k,kmol(imol(1,i)))
      end do
c
c     if all atomic masses are zero, set them all to unity
c
      do i = 1, n
         if (mass(i) .ne. 0.0d0)  goto 10
      end do
      do i = 1, n
         mass(i) = 1.0d0
      end do
   10 continue
c
c     compute the mass of each molecule and the total mass
c
      totmass = 0.0d0
      do i = 1, nmol
         molmass(i) = 0.0d0
         do k = imol(1,i), imol(2,i)
            molmass(i) = molmass(i) + mass(kmol(k))
         end do
         totmass = totmass + molmass(i)
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutate  --  set parameters for hybrid system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutate" constructs the hybrid hamiltonian for a specified
c     initial state, final state and mutation parameter "lambda"
c
c
      subroutine mutate
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'keys.i'
      include 'mutant.i'
      integer i,ihyb,it0,it1,next
      character*20 keyword
      character*80 record,string
c
c
c     zero number of hybrid atoms and hybrid atom list
c
      nhybrid = 0
      do i = 1, n
         alter(i) = .false.
      end do
c
c     search keywords for free energy perturbation options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'LAMBDA ') then
            nhybrid = 0
            lambda = 0.0d0
            string = record(next:80)
            read (string,*,err=10)  lambda
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:80)
            read (string,*,err=10)  ihyb,it0,it1
            nhybrid = nhybrid + 1
            ihybrid(nhybrid) = ihyb
            type0(nhybrid) = it0
            type1(nhybrid) = it1
            class0(nhybrid) = atmcls(it0)
            class1(nhybrid) = atmcls(it1)
            alter(ihyb) = .true.
         end if
   10    continue
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nextarg  --  find next command line argument  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nextarg" finds the next unused command line argument
c     and returns it in the input character string
c
c
      subroutine nextarg (string,exist)
      implicit none
      include 'argue.i'
      integer i,length
      character*(*) string
      logical exist
c
c
c     initialize the command argument as a blank string
c
      string = '          '
      exist = .false.
c
c     get the next command line argument and mark it as used
c
      if (narg .ne. 0) then
         length = min(len(string),len(arg(maxarg)))
         do i = 1, narg
            if (listarg(i)) then
               listarg(i) = .false.
               string = arg(i)(1:length)
               exist = .true.
               goto 10
            end if
         end do
   10    continue
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  function nexttext  --  find next non-blank character  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "nexttext" finds and returns the location of the first
c     non-blank character within an input text string; zero
c     is returned if no such character is found
c
c
      function nexttext (string)
      implicit none
      integer i,size,len,nexttext
      character*(*) string
c
c
c     move forward through the string, one character
c     at a time, looking for first non-blank character
c
      nexttext = 0
      size = len(string)
      do i = 1, size
         if (string(i:i) .gt. ' ') then
            nexttext = i
            goto 10
         end if
      end do
   10 continue
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  function number  --  convert text string to number  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "number" converts a text numeral into an integer value;
c     the input string must contain only numeric characters
c
c
      function number (string)
      implicit none
      include 'inform.i'
      include 'iounit.i'
      integer i,j,number
      integer first,last,trimtext
      integer digit,place(10)
      character*1 letter
      character*(*) string
      data place  / 1, 10, 100, 1000, 10000, 100000, 1000000,
     &              10000000, 100000000, 1000000000 /
c
c
c     initialize the integer value of number to zero
c
      number = 0
c
c     get the first and last nonblank characters
c
      last = trimtext (string)
      if (last .gt. 10) then
         write (iout,10)
   10    format (' NUMBER  --  Input Text String is Too Long')
         return
      end if
      first = 1
      do i = 1, last
         letter = string(i:i)
         if (letter .ne. ' ') then
            first = i
            goto 20
         end if
      end do
   20 continue
c
c     convert the text numeral into an integer number
c
      j = 0
      do i = last, first, -1
         j = j + 1
         letter = string(i:i)
         if (letter .eq. '0') then
            digit = 0
         else if (letter .eq. '1') then
            digit = 1
         else if (letter .eq. '2') then
            digit = 2
         else if (letter .eq. '3') then
            digit = 3
         else if (letter .eq. '4') then
            digit = 4
         else if (letter .eq. '5') then
            digit = 5
         else if (letter .eq. '6') then
            digit = 6
         else if (letter .eq. '7') then
            digit = 7
         else if (letter .eq. '8') then
            digit = 8
         else if (letter .eq. '9') then
            digit = 9
         else
            if (debug) then
               write (iout,30)
   30          format (/,' NUMBER  --  Non-Numeric Characters Found',
     &                    ' in Numeral String')
            end if
            number = 0
            goto 40
         end if
         number = number + digit * place(j)
      end do
   40 continue
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine numeral  --  convert number to text string  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "numeral" converts an input integer number into the
c     corresponding right- or left-justified text numeral
c
c     number  integer value of the number to be transformed
c     string  text string to be filled with corresponding numeral
c     size    on input, the minimal acceptable numeral length, if
c               zero then output will be right justified, if
c               nonzero then numeral is left-justified and padded
c               with leading zeros as necessary; upon output, the
c               number of non-blank characters in the numeral
c
c
      subroutine numeral (number,string,size)
      implicit none
      integer i,number,size,multi,pos
      integer length,minsize,len
      integer million,hunthou,tenthou
      integer thousand,hundred,tens,ones
      character*1 digit(0:9)
      character*(*) string
      logical right,negative
      data digit / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     set justification and size bounds for numeral string
c
      if (size .eq. 0) then
         right = .true.
         size = 1
      else
         right = .false.
      end if
      minsize = size
      length = len(string)
c
c     test the sign of the original number
c
      if (number .ge. 0) then
         negative = .false.
      else
         negative = .true.
         number = -number
      end if
c
c     use modulo arithmetic to find place-holding digits
c
      million = number / 1000000
      multi = 1000000 * million
      hunthou = (number-multi) / 100000
      multi = multi + 100000*hunthou
      tenthou = (number-multi) / 10000
      multi = multi + 10000*tenthou
      thousand = (number-multi) / 1000
      multi = multi + 1000*thousand
      hundred = (number-multi) / 100
      multi = multi + 100*hundred
      tens = (number-multi) / 10
      multi = multi + 10*tens
      ones = number - multi
c
c     find the correct length to be used for the numeral
c
      if (million .ne. 0) then
         size = 7
      else if (hunthou .ne. 0) then
         size = 6
      else if (tenthou .ne. 0) then
         size = 5
      else if (thousand .ne. 0) then
         size = 4
      else if (hundred .ne. 0) then
         size = 3
      else if (tens .ne. 0) then
         size = 2
      else
         size = 1
      end if
      size = min(size,length)
      size = max(size,minsize)
c
c     convert individual digits to a string of numerals
c
      if (size .eq. 7) then
         string(1:1) = digit(million)
         string(2:2) = digit(hunthou)
         string(3:3) = digit(tenthou)
         string(4:4) = digit(thousand)
         string(5:5) = digit(hundred)
         string(6:6) = digit(tens)
         string(7:7) = digit(ones)
      else if (size .eq. 6) then
         string(1:1) = digit(hunthou)
         string(2:2) = digit(tenthou)
         string(3:3) = digit(thousand)
         string(4:4) = digit(hundred)
         string(5:5) = digit(tens)
         string(6:6) = digit(ones)
      else if (size .eq. 5) then
         string(1:1) = digit(tenthou)
         string(2:2) = digit(thousand)
         string(3:3) = digit(hundred)
         string(4:4) = digit(tens)
         string(5:5) = digit(ones)
      else if (size .eq. 4) then
         string(1:1) = digit(thousand)
         string(2:2) = digit(hundred)
         string(3:3) = digit(tens)
         string(4:4) = digit(ones)
      else if (size .eq. 3) then
         string(1:1) = digit(hundred)
         string(2:2) = digit(tens)
         string(3:3) = digit(ones)
      else if (size .eq. 2) then
         string(1:1) = digit(tens)
         string(2:2) = digit(ones)
      else
         string(1:1) = digit(ones)
      end if
c
c     right-justify if desired, and pad with blanks
c
      if (right) then
         do i = size, 1, -1
            pos = length - size + i
            string(pos:pos) = string(i:i)
         end do
         do i = 1, length-size
            string(i:i) = ' '
         end do
      else
         do i = size+1, length
            string(i:i) = ' '
         end do
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine numgrad  --  numerical gradient of a function  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "numgrad" computes the gradient of the objective function
c     "fvalue" with respect to Cartesian coordinates of the atoms
c     via a two-sided numerical differentiation
c
c
      subroutine numgrad (fvalue,g,eps)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i
      real*8 fvalue,g(3,maxatm)
      real*8 eps,old,f,f0
      external fvalue
c
c
c     compute two-sided numerical gradient from function values
c
      do i = 1, n
         old = x(i)
         x(i) = x(i) - 0.5d0*eps
         f0 = fvalue ()
         x(i) = x(i) + eps
         f = fvalue ()
         x(i) = old
         g(1,i) = (f - f0) / eps
         old = y(i)
         y(i) = y(i) - 0.5d0*eps
         f0 = fvalue ()
         y(i) = y(i) + eps
         f = fvalue ()
         y(i) = old
         g(2,i) = (f - f0) / eps
         old = z(i)
         z(i) = z(i) - 0.5d0*eps
         f0 = fvalue ()
         z(i) = z(i) + eps
         f = fvalue ()
         z(i) = old
         g(3,i) = (f - f0) / eps
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine orbital  --  setup for pisystem calculation  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "orbital" finds and organizes lists of atoms in a pisystem,
c     bonds connecting pisystem atoms and torsions whose two
c     central atoms are both pisystem atoms
c
c
      subroutine orbital
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'iounit.i'
      include 'keys.i'
      include 'pistuf.i'
      include 'potent.i'
      include 'tors.i'
      integer i,j,k,m,ib,ic
      integer iorb,jorb,next
      integer piatoms(maxpi)
      character*20 keyword
      character*80 record,string
c
c
c     set the default values for the pisystem variables
c
      do i = 1, maxpi
         piatoms(i) = 0
      end do
c
c     check the keywords for any lists of pisystem atoms
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'PISYSTEM ') then
            string = record(next:80)
            read (string,*,err=10,end=10)  (piatoms(k),k=1,maxpi)
         end if
   10    continue
      end do
c
c     quit if no pisystem was found for consideration
c
      if (piatoms(1) .eq. 0) then
         use_orbit = .false.
         return
      else
         use_orbit = .true.
      end if
c
c     organize and make lists of the pisystem atoms
c
      do i = 1, n
         listpi(i) = .false.
      end do
      i = 1
      dowhile (piatoms(i) .ne. 0)
         if (piatoms(i) .gt. 0) then
            listpi(piatoms(i)) = .true.
            i = i + 1
         else
            do j = -piatoms(i), piatoms(i+1)
               listpi(j) = .true.
            end do
            i = i + 2
         end if
      end do
      norbit = 0
      do i = 1, n
         if (listpi(i)) then
            norbit = norbit + 1
            iorbit(norbit) = i
         end if
      end do
c
c     quit if the molecule contains too many piorbitals
c
      if (norbit .gt. maxpi) then
         write (iout,20)
   20    format (' ORBITAL  --  Too many Pi-Orbitals;',
     &           ' Increase MAXPI')
         call fatal
      end if
c
c     find three atoms which define a plane
c     perpendicular to each orbital
c
      call piplane
c
c     find and store the pisystem bonds
c
      npibond = 0
      do i = 1, norbit-1
         iorb = iorbit(i)
         do j = i, norbit
            jorb = iorbit(j)
            do k = 1, n12(iorb)
               if (i12(k,iorb) .eq. jorb) then
                  npibond = npibond + 1
                  do m = 1, nbond
                     if (iorb.eq.ibnd(1,m) .and.
     &                   jorb.eq.ibnd(2,m)) then
                        pibond(1,npibond) = m
                        pibond(2,npibond) = i
                        pibond(3,npibond) = j
                        goto 30
                     end if
                  end do
   30             continue
               end if
            end do
         end do
      end do
c
c     find and store the pisystem torsions
c
      npitors = 0
      do i = 1, ntors
         ib = itors(2,i)
         ic = itors(3,i)
         if (listpi(ib) .and. listpi(ic)) then
            do j = 1, npibond
               k = pibond(1,j)
               if (ib.eq.ibnd(1,k).and.ic.eq.ibnd(2,k) .or.
     &             ib.eq.ibnd(2,k).and.ic.eq.ibnd(1,k)) then
                  npitors = npitors + 1
                  pitors(1,npitors) = i
                  pitors(2,npitors) = j
                  goto 40
               end if
            end do
   40       continue
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine orient  --  rigid body reference coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "orient" computes a set of reference Cartesian coordinates
c     in standard orientation for each rigid body atom group
c
c
      subroutine orient
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      include 'rigid.i'
      integer i,j,k
      integer init,stop
      real*8 xcm,ycm,zcm
      real*8 phi,theta,psi
      real*8 xterm,yterm,zterm
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 a(3,3)
c
c
c     use current coordinates as default reference coordinates
c
      do i = 1, n
         xrb(i) = x(i)
         yrb(i) = y(i)
         zrb(i) = z(i)
      end do
c
c     compute the rigid body coordinates for each atom group
c
      call xyzrigid
c
c     get the center of mass and Euler angles for each group     
c
      do i = 1, ngrp
         xcm = rbc(1,i)
         ycm = rbc(2,i)
         zcm = rbc(3,i)
         phi = rbc(4,i)
         theta = rbc(5,i)
         psi = rbc(6,i)
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
         cpsi = cos(psi)
         spsi = sin(psi)
c
c     construct the rotation matrix from Euler angle values
c
         a(1,1) = ctheta * cphi
         a(2,1) = spsi*stheta*cphi - cpsi*sphi
         a(3,1) = cpsi*stheta*cphi + spsi*sphi
         a(1,2) = ctheta * sphi
         a(2,2) = spsi*stheta*sphi + cpsi*cphi
         a(3,2) = cpsi*stheta*sphi - spsi*cphi
         a(1,3) = -stheta
         a(2,3) = ctheta * spsi
         a(3,3) = ctheta * cpsi
c
c     translate and rotate each atom group into inertial frame
c
         init = igrp(1,i)
         stop = igrp(2,i)
         do j = init, stop
            k = kgrp(j)
            xterm = x(k) - xcm
            yterm = y(k) - ycm
            zterm = z(k) - zcm
            xrb(k) = a(1,1)*xterm + a(1,2)*yterm + a(1,3)*zterm
            yrb(k) = a(2,1)*xterm + a(2,2)*yterm + a(2,3)*zterm
            zrb(k) = a(3,1)*xterm + a(3,2)*yterm + a(3,3)*zterm
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine xyzrigid  --  determine rigid body coordinates  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "xyzrigid" computes the center of mass and Euler angle rigid
c     body coordinates for each atom group in the system
c
c     literature reference:
c
c     Herbert Goldstein, "Classical Mechanics, 2nd Edition",
c     Addison-Wesley, Reading, MA, 1980; see the Euler angle
c     xyz convention in Appendix B
c
c
      subroutine xyzrigid
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'group.i'
      include 'rigid.i'
      integer i,j,k,m
      integer init,stop
      real*8 xcm,ycm,zcm
      real*8 phi,theta,psi
      real*8 weigh,total,dot
      real*8 xx,xy,xz,yy,yz,zz
      real*8 xterm,yterm,zterm
      real*8 moment(3),vec(3,3)
      real*8 work1(3),work2(3)
      real*8 tensor(3,3),a(3,3)
c
c
c     get the first and last atom in the current group
c
      do i = 1, ngrp
         init = igrp(1,i)
         stop = igrp(2,i)
c
c     compute the position of the group center of mass
c
         total = 0.0d0
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         do j = init, stop
            k = kgrp(j)
            weigh = mass(k)
            total = total + weigh
            xcm = xcm + x(k)*weigh
            ycm = ycm + y(k)*weigh
            zcm = zcm + z(k)*weigh
         end do
         if (total .ne. 0.0d0) then
            xcm = xcm / total
            ycm = ycm / total
            zcm = zcm / total
         end if
c
c     compute and then diagonalize the inertial tensor
c
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         do j = init, stop
            k = kgrp(j)
            weigh = mass(k)
            xterm = x(k) - xcm
            yterm = y(k) - ycm
            zterm = z(k) - zcm
            xx = xx + xterm*xterm*weigh
            xy = xy + xterm*yterm*weigh
            xz = xz + xterm*zterm*weigh
            yy = yy + yterm*yterm*weigh
            yz = yz + yterm*zterm*weigh
            zz = zz + zterm*zterm*weigh
         end do
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
         call tnk_jacobi (3,3,tensor,moment,vec,work1,work2)
c
c     select the direction for each principle moment axis
c
         do m = 1, 2
            do j = init, stop
               k = kgrp(j)
               xterm = vec(1,m) * (x(k)-xcm)
               yterm = vec(2,m) * (y(k)-ycm)
               zterm = vec(3,m) * (z(k)-zcm)
               dot = xterm + yterm + zterm
               if (dot .lt. 0.0d0) then
                  vec(1,m) = -vec(1,m)
                  vec(2,m) = -vec(2,m)
                  vec(3,m) = -vec(3,m)
               end if
               if (dot .ne. 0.0d0)  goto 10
            end do
   10       continue
         end do
c
c     moment axes must give a right-handed coordinate system
c
         xterm = vec(1,1) * (vec(2,2)*vec(3,3)-vec(2,3)*vec(3,2))
         yterm = vec(2,1) * (vec(1,3)*vec(3,2)-vec(1,2)*vec(3,3))
         zterm = vec(3,1) * (vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2))
         dot = xterm + yterm + zterm
         if (dot .lt. 0.0d0) then
            do j = 1, 3
               vec(j,3) = -vec(j,3)
            end do
         end if
c
c     principal moment axes form rows of Euler rotation matrix
c
         do k = 1, 3
            do j = 1, 3
               a(k,j) = vec(j,k)
            end do
         end do
c
c     compute Euler angles consistent with the rotation matrix
c
         call roteuler (a,phi,theta,psi)
c
c     set the rigid body coordinates for each atom group
c
         rbc(1,i) = xcm
         rbc(2,i) = ycm
         rbc(3,i) = zcm
         rbc(4,i) = phi
         rbc(5,i) = theta
         rbc(6,i) = psi
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine roteuler  --  rotation matrix to Euler angles   ##
c     ##                                                             ##
c     #################################################################
c
c
c     "roteuler" computes a set of Euler angle values consistent
c     with an input rotation matrix
c
c
      subroutine roteuler (a,phi,theta,psi)
      implicit none
      include 'math.i'
      integer i
      real*8 phi,theta,psi,eps
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 a(3,3),b(3)
      logical err(3)
c
c
c     set the tolerance for Euler angles and rotation elements
c
      eps = 1.0d-8
c
c     get a trial value of theta from a single rotation element
c
      theta = asin(min(1.0d0,max(-1.0d0,-a(1,3))))
      ctheta = cos(theta)
      stheta = -a(1,3)
c
c     set the phi/psi difference when theta is either 90 or -90
c
      if (abs(ctheta) .le. eps) then
         phi = 0.0d0
         if (abs(a(3,1)) .lt. eps) then
            psi = asin(min(1.0d0,max(-1.0d0,-a(2,1)/a(1,3))))
         else if (abs(a(2,1)) .lt. eps) then
            psi = acos(min(1.0d0,max(-1.0d0,-a(3,1)/a(1,3))))
         else
            psi = atan(a(2,1)/a(3,1))
         end if
c
c     set the phi and psi values for all other theta values
c
      else
         if (abs(a(1,1)) .lt. eps) then
            phi = asin(min(1.0d0,max(-1.0d0,a(1,2)/ctheta)))
         else if (abs(a(1,2)) .lt. eps) then
            phi = acos(min(1.0d0,max(-1.0d0,a(1,1)/ctheta)))
         else
            phi = atan(a(1,2)/a(1,1))
         end if
         if (abs(a(3,3)) .lt. eps) then
            psi = asin(min(1.0d0,max(-1.0d0,a(2,3)/ctheta)))
         else if (abs(a(2,3)) .lt. eps) then
            psi = acos(min(1.0d0,max(-1.0d0,a(3,3)/ctheta)))
         else 
            psi = atan(a(2,3)/a(3,3))
         end if
      end if
c
c     find sine and cosine of the trial phi and psi values
c
      cphi = cos(phi)
      sphi = sin(phi)
      cpsi = cos(psi)
      spsi = sin(psi)
c
c     reconstruct the diagonal of the rotation matrix
c
      b(1) = ctheta * cphi
      b(2) = spsi*stheta*sphi + cpsi*cphi
      b(3) = ctheta * cpsi
c
c     compare the correct matrix diagonal to rebuilt diagonal
c
      do i = 1, 3
         err(i) = .false.
         if (abs(a(i,i)-b(i)) .gt. eps)  err(i) = .true.
      end do
c
c     alter Euler angles to get correct rotation matrix values
c
      if (err(1) .and. err(2))  phi = phi - sign(pi,phi)
      if (err(1) .and. err(3))  theta = -theta + sign(pi,theta)
      if (err(2) .and. err(3))  psi = psi - sign(pi,psi)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine rigidxyz  --  rigid body to Cartesian coords  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "rigidxyz" computes Cartesian coordinates for a rigid body
c     group via rotation and translation of reference coordinates  
c
c     literature reference:
c
c     Herbert Goldstein, "Classical Mechanics, 2nd Edition",
c     Addison-Wesley, Reading, MA, 1980; see the Euler angle
c     xyz convention in Appendix B
c
c
      subroutine rigidxyz
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      include 'rigid.i'
      integer i,j,k
      integer init,stop
      real*8 xcm,ycm,zcm
      real*8 phi,theta,psi
      real*8 xterm,yterm,zterm
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 a(3,3)
c
c
c     get the center of mass and Euler angles for each group     
c
      do i = 1, ngrp
         xcm = rbc(1,i)
         ycm = rbc(2,i)
         zcm = rbc(3,i)
         phi = rbc(4,i)
         theta = rbc(5,i)
         psi = rbc(6,i)
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
         cpsi = cos(psi)
         spsi = sin(psi)
c
c     construct the rotation matrix from Euler angle values
c
         a(1,1) = ctheta * cphi
         a(2,1) = spsi*stheta*cphi - cpsi*sphi
         a(3,1) = cpsi*stheta*cphi + spsi*sphi
         a(1,2) = ctheta * sphi
         a(2,2) = spsi*stheta*sphi + cpsi*cphi
         a(3,2) = cpsi*stheta*sphi - spsi*cphi
         a(1,3) = -stheta
         a(2,3) = ctheta * spsi
         a(3,3) = ctheta * cpsi
c
c     rotate and translate reference coordinates into global frame
c
         init = igrp(1,i)
         stop = igrp(2,i)
         do j = init, stop
            k = kgrp(j)
            xterm = xrb(k)
            yterm = yrb(k)
            zterm = zrb(k)
            x(k) = a(1,1)*xterm + a(2,1)*yterm + a(3,1)*zterm + xcm
            y(k) = a(1,2)*xterm + a(2,2)*yterm + a(3,2)*zterm + ycm
            z(k) = a(1,3)*xterm + a(2,3)*yterm + a(3,3)*zterm + zcm
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine orthog  --  Gram-Schmidt orthogonalization  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "orthog" performs an orthogonalization of an input matrix
c     via the modified Gram-Schmidt algorithm
c
c       m    first logical dimension of matrix to orthogonalize
c       mp   first physical dimension of matrix storage area
c       n    second logical dimension of matrix to orthogonalize
c       a    matrix to orthogonalize; contains result on exit
c
c
      subroutine orthog (m,mp,n,a)
      implicit none
      integer i,j,k,m,n,mp
      real*8 a(mp,*),rkk,rkj
c
c
c     compute the modified Gram-Schmidt orthogonalization
c
      do k = 1, n
         rkk = 0.0d0
         do i = 1, m
            rkk = rkk + a(i,k)**2
         end do
         rkk = sqrt(rkk)
         do i = 1, m
            a(i,k) = a(i,k) / rkk
         end do
         do j = k+1, n
            rkj = 0.0d0
            do i = 1, m
               rkj = rkj + a(i,k)*a(i,j)
            end do
            do i = 1, m
               a(i,j) = a(i,j) - a(i,k)*rkj
            end do
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine overlap  --  p-orbital overlap for pisystem  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "overlap" computes the overlap for two parallel p-orbitals
c     given the atomic numbers and distance of separation
c
c
      subroutine overlap (atmnum1,atmnum2,rang,ovlap)
      implicit none
      include 'units.i'
      integer atmnum1,atmnum2
      integer na,nb,la,lb
      real*8 ovlap,zeta(18)
      real*8 rbohr,rang,za,zb,s(3)
      save zeta
c
c     Slater orbital exponents for hydrogen through argon
c
      data zeta  / 1.000, 1.700, 0.650, 0.975, 1.300, 1.625,
     &             1.950, 2.275, 2.600, 2.925, 0.733, 0.950,
     &             1.167, 1.383, 1.600, 1.817, 2.033, 2.250 /
c
c
c     principal quantum number from atomic number
c
      na = 2
      nb = 2
      if (atmnum1 .gt. 10)  na = 3
      if (atmnum2 .gt. 10)  nb = 3
c
c     azimuthal quantum number for p-orbitals
c
      la = 1
      lb = 1
c
c     orbital exponent from stored ideal values
c
      za = zeta(atmnum1)
      zb = zeta(atmnum2)
c
c     convert interatomic distance to bohrs
c
      rbohr = rang / bohr
c
c     get pi-overlap via generic overlap integral routine
c
      call slatert (na,la,za,nb,lb,zb,rbohr,s)
      ovlap = s(2)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine slater  --  find overlap integrals for STO's  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "slater" is a general routine for computing the overlap
c     integrals between two Slater-type orbitals
c
c     literature reference:
c
c     D. B. Cook, "Structures and Approximations for Electrons in
c     Molecules", Ellis Horwood Limited, Sussex, England, 1978;
c     this version is adapted from the code provided in Chapter 7
c
c     variables and parameters:
c
c     na   principle quantum number for first orbital
c     la   azimuthal quantum number for first orbital
c     za   orbital exponent for the first orbital
c     nb   principle quantum number for second orbital
c     lb   azimuthal quantum number for second orbital
c     zb   orbital exponent for the second orbital
c     r    interatomic distance in atomic units
c     s    vector containing the sigma-sigma, pi-pi
c            and delta-delta overlaps upon output
c
c
cjrs 
      subroutine slatert (na,la,za,nb,lb,zb,r,s)
cjrs
      implicit none
      real*8 rmin,eps
      parameter (rmin=0.000001d0)
      parameter (eps=0.00000001d0)
      integer j,k,m,na,nb,la,lb,ja,jb,nn,max,maxx
      integer novi,ia(200),ib(200),idsga(5),idsgb(5)
      integer icosa(2),icosb(2),isina(4),isinb(4)
      real*8 s(3),r,za,zb,fact(15),cjkm
      real*8 a(20),b(20),c(200),cbase(20),theta(6)
      real*8 cosa(2),cosb(2),sinab(4),dsiga(5),dsigb(5)
      real*8 an,ana,anb,anr,rhalf,coef,p,pt
      logical done
      save icosa,icosb,cosa,cosb
      save idsga,idsgb,dsiga,dsigb
      save isina,isinb,sinab,theta,fact
      external cjkm
      data icosa / 0, 1 /
      data icosb / 0, 1 /
      data cosa  /  1.0d0, 1.0d0 /
      data cosb  / -1.0d0, 1.0d0 /
      data idsga / 0, 1, 2, 2, 0 /
      data idsgb / 0, 1, 2, 0, 2 /
      data dsiga / 3.0d0, 4.0d0, 3.0d0, -1.0d0, -1.0d0 /
      data dsigb / 3.0d0,-4.0d0, 3.0d0, -1.0d0, -1.0d0 /
      data isina / 0, 2, 0, 2 /
      data isinb / 0, 0, 2, 2 /
      data sinab / -1.0d0, 1.0d0, 1.0d0, -1.0d0 /
      data theta / 0.7071068d0, 1.2247450d0, 0.8660254d0,
     &             0.7905694d0, 1.9364916d0, 0.9682458d0 /
      data fact  / 1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, 120.0d0,
     &             720.0d0, 5040.0d0, 40320.0d0, 362880.0d0,
     &             3628800.0d0, 39916800.0d0, 479001600.0d0,
     &             6227020800.0d0, 87178291200.0d0 /
c
c
c     zero out the overlap integrals
c
      done = .false.
      s(1) = 0.0d0
      s(2) = 0.0d0
      s(3) = 0.0d0
      ana = (2.0d0*za)**(2*na+1) / fact(2*na+1)
      anb = (2.0d0*zb)**(2*nb+1) / fact(2*nb+1)
c
c     orbitals are on the same atomic center
c
      if (r .lt. rmin) then
         anr = 1.0d0
         j = na + nb + 1
         s(1) = fact(j) / ((za+zb)**j)
         an = sqrt(ana*anb)
         do novi = 1, 3
            s(novi) = s(novi) * an * anr
         end do
         return
      end if
c
c     compute overlap integrals for general case
c
      rhalf = 0.5d0 * r
      p = rhalf * (za+zb)
      pt = rhalf * (za-zb)
      nn = na + nb
      call aset (p,nn,a)
      call bset (pt,nn,b)
      k = na - la
      m = nb - lb
      max = k + m + 1
      do j = 1, max
         ia(j) = j - 1
         ib(j) = max - j
         cbase(j) = cjkm(j-1,k,m)
         c(j) = cbase(j)
      end do
      maxx = max
      if (la .eq. 1) then
         call polyp (c,ia,ib,maxx,cosa,icosa,icosb,2)
      else if (la .eq. 2) then
         call polyp (c,ia,ib,maxx,dsiga,idsga,idsgb,5)
      end if
      if (lb .eq. 1) then
         call polyp (c,ia,ib,maxx,cosb,icosa,icosb,2)
      else if (lb .eq. 2) then
         call polyp (c,ia,ib,maxx,dsigb,idsga,idsgb,5)
      end if
      novi = 1
      dowhile (.not. done)
         do j = 1, maxx
            ja = ia(j) + 1
            jb = ib(j) + 1
            coef = c(j)
            if (abs(coef) .ge. eps) then
               s(novi) = s(novi) + coef*a(ja)*b(jb)
            end if
         end do
         ja = la*(la+1)/2 + novi
         jb = lb*(lb+1)/2 + novi
         s(novi) = s(novi) * theta(ja) * theta(jb)
         if (novi.eq.1 .and. la.ne.0 .and. lb.ne.0) then
            maxx = max
            do j = 1, maxx
               c(j) = cbase(j)
            end do
            call polyp (c,ia,ib,maxx,sinab,isina,isinb,4)
            if (la .eq. 2) then
               call polyp (c,ia,ib,maxx,cosa,icosa,icosb,2)
            end if
            if (lb .eq. 2) then
               call polyp (c,ia,ib,maxx,cosb,icosa,icosb,2)
            end if
            novi = 2
         else if (novi.eq.2 .and. la.eq.2 .and. lb.eq.2) then
            maxx = max
            do j = 1, maxx
               c(j) = cbase(j)
            end do
            call polyp (c,ia,ib,maxx,sinab,isina,isinb,4)
            call polyp (c,ia,ib,maxx,sinab,isina,isinb,4)
            novi = 3
         else
            anr = rhalf**(na+nb+1)
            an = sqrt(ana*anb)
            do novi = 1, 3
               s(novi) = s(novi) * an * anr
            end do
            done = .true.
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polyp  --  polynomial product for STO overlap  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine polyp (c,ia,ib,max,d,iaa,ibb,n)
      implicit none
      integer i,j,k,m,max,n
      integer ia(200),ib(200),iaa(n),ibb(n)
      real*8 c(200),d(n)
c
c
      do j = 1, max
         do k = 1, n
            i = n - k + 1
            m = (i-1)*max + j
            c(m) = c(j) * d(i)
            ia(m) = ia(j) + iaa(i)
            ib(m) = ib(j) + ibb(i)
         end do
      end do
      max = n * max
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function cjkm  --  coefficients of spherical harmonics  ##
c     ##                                                          ##
c     ##############################################################
c
c
      function cjkm (j,k,m)
      implicit none
      integer i,j,k,m,min,max,id,idd,ip1
      real*8 cjkm,fact(15),b1,b2,sum
      save fact
      data fact  / 1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, 120.0d0,
     &             720.0d0, 5040.0d0, 40320.0d0, 362880.0d0,
     &             3628800.0d0, 39916800.0d0, 479001600.0d0,
     &             6227020800.0d0, 87178291200.0d0 /
c
c
      min = 1
      if (j .gt. m)  min = j - m + 1
      max = j + 1
      if (k .lt. j)  max = k + 1
      sum = 0.0d0
      do ip1 = min, max
         i = ip1 - 1
         id = k - i + 1
         b1 = fact(k+1) / (fact(i+1)*fact(id))
         if (j .lt. i) then
            b2 = 1.0d0
         else
            id = m - (j-i) + 1
            idd = j - i + 1
            b2 = fact(m+1) / (fact(idd)*fact(id))
         end if
         sum = sum + b1*b2*(-1.0d0)**i
      end do
      cjkm = sum * (-1.0d0)**(m-j)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine aset  --  get "a" functions by recursion  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine aset (alpha,n,a)
      implicit none
      integer i,n
      real*8 alpha,a(20),alp
c
c
      alp = 1.0d0 / alpha
      a(1) = exp(-alpha) * alp
      do i = 1, n
         a(i+1) = a(1) + dble(i)*a(i)*alp
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine bset  --  get "b" functions by recursion  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine bset (beta,n,b)
      implicit none
      real*8 eps
      parameter (eps=0.000001d0)
      integer i,j,n
      real*8 beta,b(20),bmax
      real*8 betam,d1,d2
      external bmax
c
c
      if (abs(beta) .lt. eps) then
         do i = 1, n+1
            b(i) = 2.0d0 / dble(i)
            if ((i/2)*2 .eq. i)  b(i) = 0.0d0
         end do
      else if (abs(beta) .gt. (dble(n)/2.3d0)) then
         d1 = exp(beta)
         d2 = 1.0d0 / d1
         betam = 1.0d0 / beta
         b(1) = (d1-d2) * betam
         do i = 1, n
            d1 = -d1
            b(i+1) = (d1-d2+dble(i)*b(i)) * betam
         end do
      else
         b(n+1) = bmax(beta,n)
         d1 = exp(beta)
         d2 = 1.0d0 / d1
         if ((n/2)*2 .ne. n)  d1 = -d1
         do i = 1, n
            j = n - i + 1
            d1 = -d1
            b(j) = (d1+d2+beta*b(j+1)) / dble(j)
         end do
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function bmax  --  find maximum order of "b" functions  ##
c     ##                                                          ##
c     ##############################################################
c
c
      function bmax (beta,n)
      implicit none
      real*8 eps
      parameter (eps=0.0000001d0)
      integer n
      real*8 bmax,beta,b,top,bot
      real*8 sum,fi,sign,term
      logical done
c
c
      done = .false.
      b = beta**2
      top = dble(n) + 1.0d0
      sum = 1.0d0 / top
      fi = 2.0d0
      sign = 2.0d0
      if ((n/2)*2 .ne. n) then
         top = top + 1.0d0
         sum = beta / top
         fi = fi + 1.0d0
         sign = -2.0d0
      end if
      term = sum
      dowhile (.not. done)
         bot = top + 2.0d0
         term = term * b * top / (fi*(fi-1.0d0)*bot)
         sum = sum + term
         if (abs(term) .le. eps) then
            done = .true.
         else
            fi = fi + 2.0d0
            top = bot
         end if
      end do
      bmax = sign * sum
      return
      end
