C 11 May 10 - DGF - parallelise echarge1 
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine eangang  --  angle-angle cross term energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "eangang" calculates the angle-angle potential energy
c
c
      subroutine eangang
      implicit none
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'usage.i'
      integer i,k,iangang
      integer ia,ib,ic,id,ie
      real*8 e,dt1,dt2,fgrp
      real*8 angle,dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xab,yab,zab,rab2
      real*8 xcb,ycb,zcb,rcb2
      real*8 xdb,ydb,zdb,rdb2
      real*8 xeb,yeb,zeb,reb2
      logical proceed
c
c
c     zero out the angle-angle cross term energy
c
      eaa = 0.0d0
c
c     calculate the angle-angle interaction energy term
c
      do iangang = 1, nangang
         i = iaa(1,iangang)
         k = iaa(2,iangang)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(1,k)
         ie = iang(3,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,5,ia,ib,ic,id,ie)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
c
c     compute the values of the two bond angles
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            rab2 = xab*xab + yab*yab + zab*zab
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            xeb = xie - xib
            yeb = yie - yib
            zeb = zie - zib
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb
            if (rab2*rcb2*rdb2*reb2 .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt1 = angle - anat(i)
               dot = xdb*xeb + ydb*yeb + zdb*zeb
               cosine = dot / sqrt(rdb2*reb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt2 = angle - anat(k)
c
c     get the angle-angle interaction energy
c
               e = aaunit * kaa(iangang) * dt1 * dt2
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total angle-angle energy
c
               eaa = eaa + e
            end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangang1  --  angle-angle energy & derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangang1" calculates the angle-angle potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine eangang1
      implicit none
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'bath.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'usage.i'
      include 'virial.i'
      integer i,k,iangang
      integer ia,ib,ic,id,ie
      real*8 e,angle,dot,cosine,fgrp
      real*8 dt1,dt2,deddt1,deddt2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xab,yab,zab,rab2
      real*8 xcb,ycb,zcb,rcb2
      real*8 xdb,ydb,zdb,rdb2
      real*8 xeb,yeb,zeb,reb2
      real*8 xp,yp,zp,rp
      real*8 xq,yq,zq,rq
      real*8 terma,termc,termd,terme
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxie,dedyie,dedzie
      logical proceed
c
c
c     zero out the angle-angle energy and first derivatives
c
      eaa = 0.0d0
      do i = 1, n
         deaa(1,i) = 0.0d0
         deaa(2,i) = 0.0d0
         deaa(3,i) = 0.0d0
      end do
c
c     find the energy of each angle-angle interaction
c
      do iangang = 1, nangang
         i = iaa(1,iangang)
         k = iaa(2,iangang)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(1,k)
         ie = iang(3,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,5,ia,ib,ic,id,ie)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
c
c     compute the values of the two bond angles
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            rab2 = xab*xab + yab*yab + zab*zab
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            xeb = xie - xib
            yeb = yie - yib
            zeb = zie - zib
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            xq = yeb*zdb - zeb*ydb
            yq = zeb*xdb - xeb*zdb
            zq = xeb*ydb - yeb*xdb
            rq = sqrt(xq*xq + yq*yq + zq*zq)
            if (rp*rq .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt1 = angle - anat(i)
               dot = xdb*xeb + ydb*yeb + zdb*zeb
               cosine = dot / sqrt(rdb2*reb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt2 = angle - anat(k)
c
c     get the energy and master chain rule terms for derivatives
c
               e = aaunit * kaa(iangang) * dt1 * dt2
               deddt1 = radian * e / dt1
               deddt2 = radian * e / dt2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  deddt1 = deddt1 * fgrp
                  deddt2 = deddt2 * fgrp
               end if
c
c     find chain rule terms for the first bond angle deviation
c
               terma = -deddt1 / (rab2*rp)
               termc = deddt1 / (rcb2*rp)
               dedxia = terma * (yab*zp-zab*yp)
               dedyia = terma * (zab*xp-xab*zp)
               dedzia = terma * (xab*yp-yab*xp)
               dedxic = termc * (ycb*zp-zcb*yp)
               dedyic = termc * (zcb*xp-xcb*zp)
               dedzic = termc * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the second bond angle deviation
c
               termd = -deddt2 / (rdb2*rq)
               terme = deddt2 / (reb2*rq)
               dedxid = termd * (ydb*zq-zdb*yq)
               dedyid = termd * (zdb*xq-xdb*zq)
               dedzid = termd * (xdb*yq-ydb*xq)
               dedxie = terme * (yeb*zq-zeb*yq)
               dedyie = terme * (zeb*xq-xeb*zq)
               dedzie = terme * (xeb*yq-yeb*xq)
c
c     get the central atom derivative terms by difference
c
               dedxib = -dedxia - dedxic - dedxid - dedxie
               dedyib = -dedyia - dedyic - dedyid - dedyie
               dedzib = -dedzia - dedzic - dedzid - dedzie
c
c     increment the total angle-angle energy and derivatives
c
               eaa = eaa + e
               deaa(1,ia) = deaa(1,ia) + dedxia
               deaa(2,ia) = deaa(2,ia) + dedyia
               deaa(3,ia) = deaa(3,ia) + dedzia
               deaa(1,ib) = deaa(1,ib) + dedxib
               deaa(2,ib) = deaa(2,ib) + dedyib
               deaa(3,ib) = deaa(3,ib) + dedzib
               deaa(1,ic) = deaa(1,ic) + dedxic
               deaa(2,ic) = deaa(2,ic) + dedyic
               deaa(3,ic) = deaa(3,ic) + dedzic
               deaa(1,id) = deaa(1,id) + dedxid
               deaa(2,id) = deaa(2,id) + dedyid
               deaa(3,id) = deaa(3,id) + dedzid
               deaa(1,ie) = deaa(1,ie) + dedxie
               deaa(2,ie) = deaa(2,ie) + dedyie
               deaa(3,ie) = deaa(3,ie) + dedzie
c
c     increment the virial for use in pressure computation
c
               if (isobaric) then
                  virx = virx + xab*dedxia + xcb*dedxic
     &                           + xdb*dedxid + xeb*dedxie
                  viry = viry + yab*dedyia + ycb*dedyic
     &                           + ydb*dedyid + yeb*dedyie
                  virz = virz + zab*dedzia + zcb*dedzic
     &                           + zdb*dedzid + zeb*dedzie
               end if
            end if
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eangang2  --  angle-angle Hessian; numerical  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eangang2" calculates the angle-angle potential energy
c     second derivatives with respect to Cartesian coordinates
c     using finite difference methods
c
c
      subroutine eangang2 (i)
      implicit none
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'group.i'
      include 'hessn.i'
      integer i,j,k,iangang
      integer ia,ib,ic,id,ie
      real*8 eps,fgrp,old,term
      real*8 d0(3,maxatm)
      logical proceed
c
c
c     set stepsize for derivatives and default group weight
c
      eps = 1.0d-7
      fgrp = 1.0d0
c
c     compute numerical angle-angle Hessian for current atom
c
      do iangang = 1, nangang
         j = iaa(1,iangang)
         k = iaa(2,iangang)
         ia = iang(1,j)
         ib = iang(2,j)
         ic = iang(3,j)
         id = iang(1,k)
         ie = iang(3,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,5,ia,ib,ic,id,ie)
         if (proceed)  proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic
     &                               .or. i.eq.id .or. i.eq.ie)
c
c     eliminate any duplicate atoms in the pair of angles
c
         if (proceed) then
            if (id.eq.ia .or. id.eq.ic)  then
               id = ie
               ie = 0
            else if (ie.eq.ia .or. ie.eq.ic) then
               ie = 0
            end if
c
c     find first derivatives for the base structure
c
            term = fgrp / eps
            call eangang2b (iangang)
            do j = 1, 3
               d0(j,ia) = deaa(j,ia)
               d0(j,ib) = deaa(j,ib)
               d0(j,ic) = deaa(j,ic)
               d0(j,id) = deaa(j,id)
               if (ie .ne. 0)  d0(j,ie) = deaa(j,ie)
            end do
c
c     find numerical x-components via perturbed structures
c
            old = x(i)
            x(i) = x(i) + eps
            call eangang2b (iangang)
            x(i) = old
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + term*(deaa(j,ia)-d0(j,ia))
               hessx(j,ib) = hessx(j,ib) + term*(deaa(j,ib)-d0(j,ib))
               hessx(j,ic) = hessx(j,ic) + term*(deaa(j,ic)-d0(j,ic))
               hessx(j,id) = hessx(j,id) + term*(deaa(j,id)-d0(j,id))
               if (ie .ne. 0)
     &            hessx(j,ie) = hessx(j,ie) + term*(deaa(j,ie)-d0(j,ie))
            end do
c
c     find numerical y-components via perturbed structures
c
            old = y(i)
            y(i) = y(i) + eps
            call eangang2b (iangang)
            y(i) = old
            do j = 1, 3
               hessy(j,ia) = hessy(j,ia) + term*(deaa(j,ia)-d0(j,ia))
               hessy(j,ib) = hessy(j,ib) + term*(deaa(j,ib)-d0(j,ib))
               hessy(j,ic) = hessy(j,ic) + term*(deaa(j,ic)-d0(j,ic))
               hessy(j,id) = hessy(j,id) + term*(deaa(j,id)-d0(j,id))
               if (ie .ne. 0)
     &            hessy(j,ie) = hessy(j,ie) + term*(deaa(j,ie)-d0(j,ie))
            end do
c
c     find numerical z-components via perturbed structures
c
            old = z(i)
            z(i) = z(i) + eps
            call eangang2b (iangang)
            z(i) = old
            do j = 1, 3
               hessz(j,ia) = hessz(j,ia) + term*(deaa(j,ia)-d0(j,ia))
               hessz(j,ib) = hessz(j,ib) + term*(deaa(j,ib)-d0(j,ib))
               hessz(j,ic) = hessz(j,ic) + term*(deaa(j,ic)-d0(j,ic))
               hessz(j,id) = hessz(j,id) + term*(deaa(j,id)-d0(j,id))
               if (ie .ne. 0)
     &            hessz(j,ie) = hessz(j,ie) + term*(deaa(j,ie)-d0(j,ie))
            end do
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eangang2b  --  angle-angle interaction derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eangang2b" calculates the angle-angle first derivatives for
c     a single interaction with respect to Cartesian coordinates;
c     used in computation of finite difference second derivatives
c
c
      subroutine eangang2b (i)
      implicit none
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'math.i'
      integer i,j,k,ia,ib,ic,id,ie
      real*8 angle,dot,cosine,term
      real*8 dt1,dt2,deddt1,deddt2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xab,yab,zab,rab2
      real*8 xcb,ycb,zcb,rcb2
      real*8 xdb,ydb,zdb,rdb2
      real*8 xeb,yeb,zeb,reb2
      real*8 xp,yp,zp,rp
      real*8 xq,yq,zq,rq
      real*8 terma,termc,termd,terme
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxie,dedyie,dedzie
c
c
c     set the coordinates of the involved atoms
c
      j = iaa(1,i)
      k = iaa(2,i)
      ia = iang(1,j)
      ib = iang(2,j)
      ic = iang(3,j)
      id = iang(1,k)
      ie = iang(3,k)
      xia = x(ia)
      yia = y(ia)
      zia = z(ia)
      xib = x(ib)
      yib = y(ib)
      zib = z(ib)
      xic = x(ic)
      yic = y(ic)
      zic = z(ic)
      xid = x(id)
      yid = y(id)
      zid = z(id)
      xie = x(ie)
      yie = y(ie)
      zie = z(ie)
c
c     zero out the first derivative components
c
      deaa(1,ia) = 0.0d0
      deaa(2,ia) = 0.0d0
      deaa(3,ia) = 0.0d0
      deaa(1,ib) = 0.0d0
      deaa(2,ib) = 0.0d0
      deaa(3,ib) = 0.0d0
      deaa(1,ic) = 0.0d0
      deaa(2,ic) = 0.0d0
      deaa(3,ic) = 0.0d0
      deaa(1,id) = 0.0d0
      deaa(2,id) = 0.0d0
      deaa(3,id) = 0.0d0
      deaa(1,ie) = 0.0d0
      deaa(2,ie) = 0.0d0
      deaa(3,ie) = 0.0d0
c
c     compute the values of the two bond angles
c
      xab = xia - xib
      yab = yia - yib
      zab = zia - zib
      rab2 = xab*xab + yab*yab + zab*zab
      xcb = xic - xib
      ycb = yic - yib
      zcb = zic - zib
      rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
      xdb = xid - xib
      ydb = yid - yib
      zdb = zid - zib
      rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
      xeb = xie - xib
      yeb = yie - yib
      zeb = zie - zib
      reb2 = xeb*xeb + yeb*yeb + zeb*zeb
      xp = ycb*zab - zcb*yab
      yp = zcb*xab - xcb*zab
      zp = xcb*yab - ycb*xab
      rp = sqrt(xp*xp + yp*yp + zp*zp)
      xq = yeb*zdb - zeb*ydb
      yq = zeb*xdb - xeb*zdb
      zq = xeb*ydb - yeb*xdb
      rq = sqrt(xq*xq + yq*yq + zq*zq)
      if (rp*rq .ne. 0.0d0) then
         dot = xab*xcb + yab*ycb + zab*zcb
         cosine = dot / sqrt(rab2*rcb2)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
         dt1 = angle - anat(j)
         dot = xdb*xeb + ydb*yeb + zdb*zeb
         cosine = dot / sqrt(rdb2*reb2)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
         dt2 = angle - anat(k)
c
c     get the energy and master chain rule terms for derivatives
c
         term = radian * aaunit * kaa(i)
         deddt1 = term * dt2
         deddt2 = term * dt1
c
c     find chain rule terms for the first bond angle deviation
c
         terma = -deddt1 / (rab2*rp)
         termc = deddt1 / (rcb2*rp)
         dedxia = terma * (yab*zp-zab*yp)
         dedyia = terma * (zab*xp-xab*zp)
         dedzia = terma * (xab*yp-yab*xp)
         dedxic = termc * (ycb*zp-zcb*yp)
         dedyic = termc * (zcb*xp-xcb*zp)
         dedzic = termc * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the second bond angle deviation
c
         termd = -deddt2 / (rdb2*rq)
         terme = deddt2 / (reb2*rq)
         dedxid = termd * (ydb*zq-zdb*yq)
         dedyid = termd * (zdb*xq-xdb*zq)
         dedzid = termd * (xdb*yq-ydb*xq)
         dedxie = terme * (yeb*zq-zeb*yq)
         dedyie = terme * (zeb*xq-xeb*zq)
         dedzie = terme * (xeb*yq-yeb*xq)
c
c     get the central atom derivative terms by difference
c
         dedxib = -dedxia - dedxic - dedxid - dedxie
         dedyib = -dedyia - dedyic - dedyid - dedyie
         dedzib = -dedzia - dedzic - dedzid - dedzie
c
c     set the angle-angle interaction first derivatives
c
         deaa(1,ia) = deaa(1,ia) + dedxia
         deaa(2,ia) = deaa(2,ia) + dedyia
         deaa(3,ia) = deaa(3,ia) + dedzia
         deaa(1,ib) = deaa(1,ib) + dedxib
         deaa(2,ib) = deaa(2,ib) + dedyib
         deaa(3,ib) = deaa(3,ib) + dedzib
         deaa(1,ic) = deaa(1,ic) + dedxic
         deaa(2,ic) = deaa(2,ic) + dedyic
         deaa(3,ic) = deaa(3,ic) + dedzic
         deaa(1,id) = deaa(1,id) + dedxid
         deaa(2,id) = deaa(2,id) + dedyid
         deaa(3,id) = deaa(3,id) + dedzid
         deaa(1,ie) = deaa(1,ie) + dedxie
         deaa(2,ie) = deaa(2,ie) + dedyie
         deaa(3,ie) = deaa(3,ie) + dedzie
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eangang3  --  angle-angle energy and analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eangang3" calculates the angle-angle potential energy;
c     also partitions the energy among the atoms
c
c
      subroutine eangang3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'angang.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'usage.i'
      integer i,k,iangang
      integer ia,ib,ic,id,ie
      real*8 e,dt1,dt2,fgrp
      real*8 angle,dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xab,yab,zab,rab2
      real*8 xcb,ycb,zcb,rcb2
      real*8 xdb,ydb,zdb,rdb2
      real*8 xeb,yeb,zeb,reb2
      logical header,huge,proceed
c
c
c     zero out the angle-angle cross term energy
c
      neaa = 0
      eaa = 0.0d0
      do i = 1, n
         aeaa(i) = 0.0d0
      end do
      header = .true.
c
c     find the energy of each angle-angle interaction
c
      do iangang = 1, nangang
         i = iaa(1,iangang)
         k = iaa(2,iangang)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(1,k)
         ie = iang(3,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,5,ia,ib,ic,id,ie)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
c
c     compute the values of the two bond angles
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            rab2 = xab*xab + yab*yab + zab*zab
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            xeb = xie - xib
            yeb = yie - yib
            zeb = zie - zib
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb
            if (rab2*rcb2*rdb2*reb2 .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt1 = angle - anat(i)
               dot = xdb*xeb + ydb*yeb + zdb*zeb
               cosine = dot / sqrt(rdb2*reb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt2 = angle - anat(k)
c
c     get the angle-angle interaction energy
c
               e = aaunit * kaa(iangang) * dt1 * dt2
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total angle-angle energy
c
               neaa = neaa + 1
               eaa = eaa + e
               aeaa(ib) = aeaa(ib) + e
c
c     print a warning if the energy of this interaction is large
c
               huge = (e .gt. 5.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Angle-Angle Cross',
     &                          ' Term Interactions :',
     &                       //,' Type',7x,'Center',5x,'Angle1'
     &                          ,5x,'Angle2',4x,'dAngle1'
     &                          ,3x,'dAngle2',6x,'Energy',/)
                  end if
                  write (iout,20)  ib,name(ib),ia,ic,id,ie,dt1,dt2,e
   20             format (' AngAng   ',i5,'-',a3,1x,2i5,1x,2i5,
     &                    2f10.4,f12.4)
               end if
            end if
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine eangle  --  angle bending potential energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "eangle" calculates the angle bending potential energy;
c     for trigonal centers, projected in-plane angles as per
c     the Wilson-Decius-Cross formalism are optionally used
c
c
      subroutine eangle
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'usage.i'
      integer i,ia,ib,ic,id
      real*8 e,ideal,force
      real*8 dot,cosine
      real*8 angle,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab,rab2
      real*8 xcb,ycb,zcb,rcb2
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xip,yip,zip
      real*8 xap,yap,zap,rap2
      real*8 xcp,ycp,zcp,rcp2
      real*8 xt,yt,zt,rt2,delta
      logical proceed
c
c
c     zero out the angle bending energy component
c
      ea = 0.0d0
c
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         ideal = anat(i)
         force = acon(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
c
c     compute the value of the bond angle
c
            if (.not.angin(i)) then
               xab = xia - xib
               yab = yia - yib
               zab = zia - zib
               rab2 = xab*xab + yab*yab + zab*zab
               xcb = xic - xib
               ycb = yic - yib
               zcb = zic - zib
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
                  dot = xab*xcb + yab*ycb + zab*zcb
                  cosine = dot / sqrt(rab2*rcb2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
c
c     find the bond angle bending energy
c
                  dt = angle - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the total bond angle bending energy
c
                  ea = ea + e
               end if
c
c     compute the value of the projected in-plane angle
c
            else
               xid = x(id)
               yid = y(id)
               zid = z(id)
               xad = xia - xid
               yad = yia - yid
               zad = zia - zid
               xbd = xib - xid
               ybd = yib - yid
               zbd = zib - zid
               xcd = xic - xid
               ycd = yic - yid
               zcd = zic - zid
               xt = yad*zcd - zad*ycd
               yt = zad*xcd - xad*zcd
               zt = xad*ycd - yad*xcd
               rt2 = xt*xt + yt*yt + zt*zt
               delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
               xip = xib + xt*delta
               yip = yib + yt*delta
               zip = zib + zt*delta
               xap = xia - xip
               yap = yia - yip
               zap = zia - zip
               rap2 = xap*xap + yap*yap + zap*zap
               xcp = xic - xip
               ycp = yic - yip
               zcp = zic - zip
               rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
               if (rap2.ne.0.0d0 .and. rcp2.ne.0.0d0) then
                  dot = xap*xcp + yap*ycp + zap*zcp
                  cosine = dot / sqrt(rap2*rcp2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
c
c     find the in-plane angle bending energy
c
                  dt = angle - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the total bond angle bending energy
c
                  ea = ea + e
               end if
            end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangle1  --  angle bend energy and derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangle1" calculates the angle bending potential energy and
c     the first derivatives with respect to Cartesian coordinates;
c     for trigonal centers, projected in-plane angles as per the
c     Wilson-Decius-Cross formalism are optionally used
c
c
      subroutine eangle1
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'bath.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ia,ib,ic,id
      real*8 e,ideal,force
      real*8 dot,cosine
      real*8 angle,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 deddt,term
      real*8 terma,termc
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab,rab2
      real*8 xcb,ycb,zcb,rcb2
      real*8 xdb,ydb,zdb
      real*8 xp,yp,zp,rp
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xip,yip,zip
      real*8 xap,yap,zap,rap2
      real*8 xcp,ycp,zcp,rcp2
      real*8 xt,yt,zt,rt2,ptrt2
      real*8 xm,ym,zm,rm,delta,delta2
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxip,dedyip,dedzip
      real*8 dpdxia,dpdyia,dpdzia
      real*8 dpdxic,dpdyic,dpdzic
      logical proceed
c
c
c     zero out energy and first derivative components
c
      ea = 0.0d0
      do i = 1, n
         dea(1,i) = 0.0d0
         dea(2,i) = 0.0d0
         dea(3,i) = 0.0d0
      end do
c
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         ideal = anat(i)
         force = acon(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
c
c     compute the value of the bond angle
c
            if (.not.angin(i)) then
               xab = xia - xib
               yab = yia - yib
               zab = zia - zib
               rab2 = xab*xab + yab*yab + zab*zab
               xcb = xic - xib
               ycb = yic - yib
               zcb = zic - zib
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               if (rp .ne. 0.0d0) then
                  dot = xab*xcb + yab*ycb + zab*zcb
                  cosine = dot / sqrt(rab2*rcb2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
c
c     get the energy and master chain rule terms for derivatives
c
                  dt = angle - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                  deddt = angunit * force * dt * radian
     &                      * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2
     &                           + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     deddt = deddt * fgrp
                  end if
c
c     compute derivative components for this interaction
c
                  terma = -deddt / (rab2*rp)
                  termc = deddt / (rcb2*rp)
                  dedxia = terma * (yab*zp-zab*yp)
                  dedyia = terma * (zab*xp-xab*zp)
                  dedzia = terma * (xab*yp-yab*xp)
                  dedxic = termc * (ycb*zp-zcb*yp)
                  dedyic = termc * (zcb*xp-xcb*zp)
                  dedzic = termc * (xcb*yp-ycb*xp)
                  dedxib = -dedxia - dedxic
                  dedyib = -dedyia - dedyic
                  dedzib = -dedzia - dedzic
c
c     increment the total bond angle energy and derivatives
c
                  ea = ea + e
                  dea(1,ia) = dea(1,ia) + dedxia
                  dea(2,ia) = dea(2,ia) + dedyia
                  dea(3,ia) = dea(3,ia) + dedzia
                  dea(1,ib) = dea(1,ib) + dedxib
                  dea(2,ib) = dea(2,ib) + dedyib
                  dea(3,ib) = dea(3,ib) + dedzib
                  dea(1,ic) = dea(1,ic) + dedxic
                  dea(2,ic) = dea(2,ic) + dedyic
                  dea(3,ic) = dea(3,ic) + dedzic
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx + xab*dedxia + xcb*dedxic
                     viry = viry + yab*dedyia + ycb*dedyic
                     virz = virz + zab*dedzia + zcb*dedzic
                  end if
               end if
c
c     compute the value of the projected in-plane angle
c
            else
               xid = x(id)
               yid = y(id)
               zid = z(id)
               xad = xia - xid
               yad = yia - yid
               zad = zia - zid
               xbd = xib - xid
               ybd = yib - yid
               zbd = zib - zid
               xcd = xic - xid
               ycd = yic - yid
               zcd = zic - zid
               xt = yad*zcd - zad*ycd
               yt = zad*xcd - xad*zcd
               zt = xad*ycd - yad*xcd
               rt2 = xt*xt + yt*yt + zt*zt
               delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
               xip = xib + xt*delta
               yip = yib + yt*delta
               zip = zib + zt*delta
               xap = xia - xip
               yap = yia - yip
               zap = zia - zip
               rap2 = xap*xap + yap*yap + zap*zap
               xcp = xic - xip
               ycp = yic - yip
               zcp = zic - zip
               rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
               xm = ycp*zap - zcp*yap
               ym = zcp*xap - xcp*zap
               zm = xcp*yap - ycp*xap
               rm = sqrt(xm*xm + ym*ym + zm*zm)
               if (rm .ne. 0.0d0) then
                  dot = xap*xcp + yap*ycp + zap*zcp
                  cosine = dot / sqrt(rap2*rcp2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
c
c     get the in-plane bend energy and master chain rule terms
c
                  dt = angle - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                  deddt = angunit * force * dt * radian
     &                      * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2
     &                           + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     deddt = deddt * fgrp
                  end if
c
c     chain rule terms for first derivative components
c
                  terma = -deddt / (rap2*rm)
                  termc = deddt / (rcp2*rm)
                  dedxia = terma * (yap*zm-zap*ym)
                  dedyia = terma * (zap*xm-xap*zm)
                  dedzia = terma * (xap*ym-yap*xm)
                  dedxic = termc * (ycp*zm-zcp*ym)
                  dedyic = termc * (zcp*xm-xcp*zm)
                  dedzic = termc * (xcp*ym-ycp*xm)
                  dedxip = -dedxia - dedxic
                  dedyip = -dedyia - dedyic
                  dedzip = -dedzia - dedzic
c
c     chain rule components for the projection of the central atom
c
                  delta2 = 2.0d0 * delta
                  ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2
                  term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd)
                  dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2
                  term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd)
                  dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2
                  term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd)
                  dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2
                  term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad)
                  dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2
                  term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad)
                  dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2
                  term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad)
                  dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2
c
c     compute derivative components for this interaction
c
                  dedxia = dedxia + dpdxia
                  dedyia = dedyia + dpdyia
                  dedzia = dedzia + dpdzia
                  dedxib = dedxip
                  dedyib = dedyip
                  dedzib = dedzip
                  dedxic = dedxic + dpdxic
                  dedyic = dedyic + dpdyic
                  dedzic = dedzic + dpdzic
                  dedxid = -dedxia - dedxib - dedxic
                  dedyid = -dedyia - dedyib - dedyic
                  dedzid = -dedzia - dedzib - dedzic
c
c     increment the in-plane angle bend energy and derivatives
c
                  ea = ea + e
                  dea(1,ia) = dea(1,ia) + dedxia
                  dea(2,ia) = dea(2,ia) + dedyia
                  dea(3,ia) = dea(3,ia) + dedzia
                  dea(1,ib) = dea(1,ib) + dedxib
                  dea(2,ib) = dea(2,ib) + dedyib
                  dea(3,ib) = dea(3,ib) + dedzib
                  dea(1,ic) = dea(1,ic) + dedxic
                  dea(2,ic) = dea(2,ic) + dedyic
                  dea(3,ic) = dea(3,ic) + dedzic
                  dea(1,id) = dea(1,id) + dedxid
                  dea(2,id) = dea(2,id) + dedyid
                  dea(3,id) = dea(3,id) + dedzid
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     xab = xia - xib
                     yab = yia - yib
                     zab = zia - zib
                     xcb = xic - xib
                     ycb = yic - yib
                     zcb = zic - zib
                     xdb = xid - xib
                     ydb = yid - yib
                     zdb = zid - zib
                     virx = virx + xab*dedxia + xcb*dedxic + xdb*dedxid
                     viry = viry + yab*dedyia + ycb*dedyic + ydb*dedyid
                     virz = virz + zab*dedzia + zcb*dedzic + zdb*dedzid
                  end if
               end if
            end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangle2  --  atom-wise bend & str-bend Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangle2" calculates second derivatives of the angle bending
c     energy for a single atom using a mixture of analytical and
c     finite difference methods; for trigonal centers, projected
c     in-plane angles as per the Wilson-Decius-Cross formalism are
c     optionally used
c
c
      subroutine eangle2 (i)
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'group.i'
      include 'hessn.i'
      integer i,j,k,ia,ib,ic,id
      real*8 eps,fgrp,old,term
      real*8 d0(3,maxatm)
      logical proceed
c
c
c     compute analytical angle bending Hessian elements
c
      call eangle2a (i)
c
c     set stepsize for derivatives and default group weight
c
      eps = 1.0d-7
      fgrp = 1.0d0
c
c     compute numerical in-plane bend Hessian for current atom
c
      do k = 1, nangle
         proceed = (angin(k))
         if (proceed) then
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            id = iang(4,k)
            proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic .or. i.eq.id)
         end if
         if (proceed .and. use_group)  call groups (proceed,fgrp,4,
     &                                              ia,ib,ic,id,0)
c
c     find first derivatives for the base structure
c
         if (proceed) then
            term = fgrp / eps
            call eangle2b (k)
            do j = 1, 3
               d0(j,ia) = dea(j,ia)
               d0(j,ib) = dea(j,ib)
               d0(j,ic) = dea(j,ic)
               d0(j,id) = dea(j,id)
            end do
c
c     find numerical x-components via perturbed structures
c
            old = x(i)
            x(i) = x(i) + eps
            call eangle2b (k)
            x(i) = old
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + term*(dea(j,ia)-d0(j,ia))
               hessx(j,ib) = hessx(j,ib) + term*(dea(j,ib)-d0(j,ib))
               hessx(j,ic) = hessx(j,ic) + term*(dea(j,ic)-d0(j,ic))
               hessx(j,id) = hessx(j,id) + term*(dea(j,id)-d0(j,id))
            end do
c
c     find numerical y-components via perturbed structures
c
            old = y(i)
            y(i) = y(i) + eps
            call eangle2b (k)
            y(i) = old
            do j = 1, 3
               hessy(j,ia) = hessy(j,ia) + term*(dea(j,ia)-d0(j,ia))
               hessy(j,ib) = hessy(j,ib) + term*(dea(j,ib)-d0(j,ib))
               hessy(j,ic) = hessy(j,ic) + term*(dea(j,ic)-d0(j,ic))
               hessy(j,id) = hessy(j,id) + term*(dea(j,id)-d0(j,id))
            end do
c
c     find numerical z-components via perturbed structures
c
            old = z(i)
            z(i) = z(i) + eps
            call eangle2b (k)
            z(i) = old
            do j = 1, 3
               hessz(j,ia) = hessz(j,ia) + term*(dea(j,ia)-d0(j,ia))
               hessz(j,ib) = hessz(j,ib) + term*(dea(j,ib)-d0(j,ib))
               hessz(j,ic) = hessz(j,ic) + term*(dea(j,ic)-d0(j,ic))
               hessz(j,id) = hessz(j,id) + term*(dea(j,id)-d0(j,id))
            end do
         end if
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eangle2a  --  angle bending Hessian; analytical  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eangle2a" calculates bond angle bending potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine eangle2a (iatom)
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'group.i'
      include 'hessn.i'
      include 'math.i'
      integer i,iatom,ia,ib,ic
      real*8 ideal,force
      real*8 dot,cosine
      real*8 angle,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab,rab
      real*8 xcb,ycb,zcb,rcb
      real*8 xp,yp,zp,rp,rp2
      real*8 xrab,yrab,zrab,rab2
      real*8 xrcb,yrcb,zrcb,rcb2
      real*8 xabp,yabp,zabp
      real*8 xcbp,ycbp,zcbp
      real*8 deddt,d2eddt2,terma,termc
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxib,ddtdyib,ddtdzib
      real*8 ddtdxic,ddtdyic,ddtdzic
      real*8 dxiaxia,dxiayia,dxiazia
      real*8 dxibxib,dxibyib,dxibzib
      real*8 dxicxic,dxicyic,dxiczic
      real*8 dyiayia,dyiazia,dziazia
      real*8 dyibyib,dyibzib,dzibzib
      real*8 dyicyic,dyiczic,dziczic
      real*8 dxibxia,dxibyia,dxibzia
      real*8 dyibxia,dyibyia,dyibzia
      real*8 dzibxia,dzibyia,dzibzia
      real*8 dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic
      real*8 dzibxic,dzibyic,dzibzic
      real*8 dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic
      real*8 dziaxic,dziayic,dziazic
      logical proceed
c
c
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ideal = anat(i)
         force = acon(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,3,ia,ib,ic,0,0)
         if (proceed)  proceed = (iatom.eq.ia .or. iatom.eq.ib
     &                                .or. iatom.eq.ic)
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
c
c     compute the value of the bond angle
c
            if (.not.angin(i)) then
               xab = xia - xib
               yab = yia - yib
               zab = zia - zib
               rab = sqrt(xab*xab + yab*yab + zab*zab)
               xcb = xic - xib
               ycb = yic - yib
               zcb = zic - zib
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               if (rp .ne. 0.0d0) then
                  dot = xab*xcb + yab*ycb + zab*zcb
                  cosine = dot / (rab*rcb)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
c
c     get the angle bending master chain rule terms
c
                  dt = angle - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt3 * dt
                  deddt = angunit * force * dt * radian
     &                      * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2
     &                          + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)
                  d2eddt2 = angunit * force * radian**2
     &                     * (2.0d0 + 6.0d0*cang*dt + 12.0d0*qang*dt2
     &                         + 20.0d0*pang*dt3 + 30.0d0*sang*dt4)
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     deddt = deddt * fgrp
                     d2eddt2 = d2eddt2 * fgrp
                  end if
c
c     first derivatives of bond angle with respect to coordinates
c
                  terma = -1.0d0 / (rab*rab*rp)
                  termc = 1.0d0 / (rcb*rcb*rp)
                  ddtdxia = terma * (yab*zp-zab*yp)
                  ddtdyia = terma * (zab*xp-xab*zp)
                  ddtdzia = terma * (xab*yp-yab*xp)
                  ddtdxic = termc * (ycb*zp-zcb*yp)
                  ddtdyic = termc * (zcb*xp-xcb*zp)
                  ddtdzic = termc * (xcb*yp-ycb*xp)
                  ddtdxib = -ddtdxia - ddtdxic
                  ddtdyib = -ddtdyia - ddtdyic
                  ddtdzib = -ddtdzia - ddtdzic
c
c     abbreviations used in defining chain rule terms
c
                  rab2 = 2.0d0 / (rab*rab)
                  xrab = xab * rab2
                  yrab = yab * rab2
                  zrab = zab * rab2
                  rcb2 = 2.0d0 / (rcb*rcb)
                  xrcb = xcb * rcb2
                  yrcb = ycb * rcb2
                  zrcb = zcb * rcb2
                  rp2 = 1.0d0 / (rp*rp)
                  xabp = (yab*zp-zab*yp) * rp2
                  yabp = (zab*xp-xab*zp) * rp2
                  zabp = (xab*yp-yab*xp) * rp2
                  xcbp = (ycb*zp-zcb*yp) * rp2
                  ycbp = (zcb*xp-xcb*zp) * rp2
                  zcbp = (xcb*yp-ycb*xp) * rp2
c
c     second derivatives of bond angle with respect to coordinates
c
                  dxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab)
                  dxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab)
                  dxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab)
                  dyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab)
                  dyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab)
                  dziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab)
                  dxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb)
                  dxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb)
                  dxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb)
                  dyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb)
                  dyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb)
                  dziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb)
                  dxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp
                  dxiayic = -terma*xab*yab - ddtdxia*yabp
                  dxiazic = -terma*xab*zab - ddtdxia*zabp
                  dyiaxic = -terma*xab*yab - ddtdyia*xabp
                  dyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp
                  dyiazic = -terma*yab*zab - ddtdyia*zabp
                  dziaxic = -terma*xab*zab - ddtdzia*xabp
                  dziayic = -terma*yab*zab - ddtdzia*yabp
                  dziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp
c
c     more angle deviation derivatives resulting from symmetry
c
                  dxibxia = -dxiaxia - dxiaxic
                  dxibyia = -dxiayia - dyiaxic
                  dxibzia = -dxiazia - dziaxic
                  dyibxia = -dxiayia - dxiayic
                  dyibyia = -dyiayia - dyiayic
                  dyibzia = -dyiazia - dziayic
                  dzibxia = -dxiazia - dxiazic
                  dzibyia = -dyiazia - dyiazic
                  dzibzia = -dziazia - dziazic
                  dxibxic = -dxicxic - dxiaxic
                  dxibyic = -dxicyic - dxiayic
                  dxibzic = -dxiczic - dxiazic
                  dyibxic = -dxicyic - dyiaxic
                  dyibyic = -dyicyic - dyiayic
                  dyibzic = -dyiczic - dyiazic
                  dzibxic = -dxiczic - dziaxic
                  dzibyic = -dyiczic - dziayic
                  dzibzic = -dziczic - dziazic
                  dxibxib = -dxibxia - dxibxic
                  dxibyib = -dxibyia - dxibyic
                  dxibzib = -dxibzia - dxibzic
                  dyibyib = -dyibyia - dyibyic
                  dyibzib = -dyibzia - dyibzic
                  dzibzib = -dzibzia - dzibzic
c
c     now, increment diagonal and off-diagonal Hessian elements
c
                  if (ia .eq. iatom) then
                     hessx(1,ia) = hessx(1,ia) + deddt*dxiaxia
     &                                  + d2eddt2*ddtdxia*ddtdxia
                     hessx(2,ia) = hessx(2,ia) + deddt*dxiayia
     &                                  + d2eddt2*ddtdxia*ddtdyia
                     hessx(3,ia) = hessx(3,ia) + deddt*dxiazia
     &                                  + d2eddt2*ddtdxia*ddtdzia
                     hessy(1,ia) = hessy(1,ia) + deddt*dxiayia
     &                                  + d2eddt2*ddtdyia*ddtdxia
                     hessy(2,ia) = hessy(2,ia) + deddt*dyiayia
     &                                  + d2eddt2*ddtdyia*ddtdyia
                     hessy(3,ia) = hessy(3,ia) + deddt*dyiazia
     &                                  + d2eddt2*ddtdyia*ddtdzia
                     hessz(1,ia) = hessz(1,ia) + deddt*dxiazia
     &                                  + d2eddt2*ddtdzia*ddtdxia
                     hessz(2,ia) = hessz(2,ia) + deddt*dyiazia
     &                                  + d2eddt2*ddtdzia*ddtdyia
                     hessz(3,ia) = hessz(3,ia) + deddt*dziazia
     &                                  + d2eddt2*ddtdzia*ddtdzia
                     hessx(1,ib) = hessx(1,ib) + deddt*dxibxia
     &                                  + d2eddt2*ddtdxia*ddtdxib
                     hessx(2,ib) = hessx(2,ib) + deddt*dyibxia
     &                                  + d2eddt2*ddtdxia*ddtdyib
                     hessx(3,ib) = hessx(3,ib) + deddt*dzibxia
     &                                  + d2eddt2*ddtdxia*ddtdzib
                     hessy(1,ib) = hessy(1,ib) + deddt*dxibyia
     &                                  + d2eddt2*ddtdyia*ddtdxib
                     hessy(2,ib) = hessy(2,ib) + deddt*dyibyia
     &                                  + d2eddt2*ddtdyia*ddtdyib
                     hessy(3,ib) = hessy(3,ib) + deddt*dzibyia
     &                                  + d2eddt2*ddtdyia*ddtdzib
                     hessz(1,ib) = hessz(1,ib) + deddt*dxibzia
     &                                  + d2eddt2*ddtdzia*ddtdxib
                     hessz(2,ib) = hessz(2,ib) + deddt*dyibzia
     &                                  + d2eddt2*ddtdzia*ddtdyib
                     hessz(3,ib) = hessz(3,ib) + deddt*dzibzia
     &                                  + d2eddt2*ddtdzia*ddtdzib
                     hessx(1,ic) = hessx(1,ic) + deddt*dxiaxic
     &                                  + d2eddt2*ddtdxia*ddtdxic
                     hessx(2,ic) = hessx(2,ic) + deddt*dxiayic
     &                                  + d2eddt2*ddtdxia*ddtdyic
                     hessx(3,ic) = hessx(3,ic) + deddt*dxiazic
     &                                  + d2eddt2*ddtdxia*ddtdzic
                     hessy(1,ic) = hessy(1,ic) + deddt*dyiaxic
     &                                  + d2eddt2*ddtdyia*ddtdxic
                     hessy(2,ic) = hessy(2,ic) + deddt*dyiayic
     &                                  + d2eddt2*ddtdyia*ddtdyic
                     hessy(3,ic) = hessy(3,ic) + deddt*dyiazic
     &                                  + d2eddt2*ddtdyia*ddtdzic
                     hessz(1,ic) = hessz(1,ic) + deddt*dziaxic
     &                                  + d2eddt2*ddtdzia*ddtdxic
                     hessz(2,ic) = hessz(2,ic) + deddt*dziayic
     &                                  + d2eddt2*ddtdzia*ddtdyic
                     hessz(3,ic) = hessz(3,ic) + deddt*dziazic
     &                                  + d2eddt2*ddtdzia*ddtdzic
                  else if (ib .eq. iatom) then
                     hessx(1,ib) = hessx(1,ib) + deddt*dxibxib
     &                                  + d2eddt2*ddtdxib*ddtdxib
                     hessx(2,ib) = hessx(2,ib) + deddt*dxibyib
     &                                  + d2eddt2*ddtdxib*ddtdyib
                     hessx(3,ib) = hessx(3,ib) + deddt*dxibzib
     &                                  + d2eddt2*ddtdxib*ddtdzib
                     hessy(1,ib) = hessy(1,ib) + deddt*dxibyib
     &                                  + d2eddt2*ddtdyib*ddtdxib
                     hessy(2,ib) = hessy(2,ib) + deddt*dyibyib
     &                                  + d2eddt2*ddtdyib*ddtdyib
                     hessy(3,ib) = hessy(3,ib) + deddt*dyibzib
     &                                  + d2eddt2*ddtdyib*ddtdzib
                     hessz(1,ib) = hessz(1,ib) + deddt*dxibzib
     &                                  + d2eddt2*ddtdzib*ddtdxib
                     hessz(2,ib) = hessz(2,ib) + deddt*dyibzib
     &                                  + d2eddt2*ddtdzib*ddtdyib
                     hessz(3,ib) = hessz(3,ib) + deddt*dzibzib
     &                                  + d2eddt2*ddtdzib*ddtdzib
                     hessx(1,ia) = hessx(1,ia) + deddt*dxibxia
     &                                  + d2eddt2*ddtdxib*ddtdxia
                     hessx(2,ia) = hessx(2,ia) + deddt*dxibyia
     &                                  + d2eddt2*ddtdxib*ddtdyia
                     hessx(3,ia) = hessx(3,ia) + deddt*dxibzia
     &                                  + d2eddt2*ddtdxib*ddtdzia
                     hessy(1,ia) = hessy(1,ia) + deddt*dyibxia
     &                                  + d2eddt2*ddtdyib*ddtdxia
                     hessy(2,ia) = hessy(2,ia) + deddt*dyibyia
     &                                  + d2eddt2*ddtdyib*ddtdyia
                     hessy(3,ia) = hessy(3,ia) + deddt*dyibzia
     &                                  + d2eddt2*ddtdyib*ddtdzia
                     hessz(1,ia) = hessz(1,ia) + deddt*dzibxia
     &                                  + d2eddt2*ddtdzib*ddtdxia
                     hessz(2,ia) = hessz(2,ia) + deddt*dzibyia
     &                                  + d2eddt2*ddtdzib*ddtdyia
                     hessz(3,ia) = hessz(3,ia) + deddt*dzibzia
     &                                  + d2eddt2*ddtdzib*ddtdzia
                     hessx(1,ic) = hessx(1,ic) + deddt*dxibxic
     &                                  + d2eddt2*ddtdxib*ddtdxic
                     hessx(2,ic) = hessx(2,ic) + deddt*dxibyic
     &                                  + d2eddt2*ddtdxib*ddtdyic
                     hessx(3,ic) = hessx(3,ic) + deddt*dxibzic
     &                                  + d2eddt2*ddtdxib*ddtdzic
                     hessy(1,ic) = hessy(1,ic) + deddt*dyibxic
     &                                  + d2eddt2*ddtdyib*ddtdxic
                     hessy(2,ic) = hessy(2,ic) + deddt*dyibyic
     &                                  + d2eddt2*ddtdyib*ddtdyic
                     hessy(3,ic) = hessy(3,ic) + deddt*dyibzic
     &                                  + d2eddt2*ddtdyib*ddtdzic
                     hessz(1,ic) = hessz(1,ic) + deddt*dzibxic
     &                                  + d2eddt2*ddtdzib*ddtdxic
                     hessz(2,ic) = hessz(2,ic) + deddt*dzibyic
     &                                  + d2eddt2*ddtdzib*ddtdyic
                     hessz(3,ic) = hessz(3,ic) + deddt*dzibzic
     &                                  + d2eddt2*ddtdzib*ddtdzic
                  else if (ic .eq. iatom) then
                     hessx(1,ic) = hessx(1,ic) + deddt*dxicxic
     &                                  + d2eddt2*ddtdxic*ddtdxic
                     hessx(2,ic) = hessx(2,ic) + deddt*dxicyic
     &                                  + d2eddt2*ddtdxic*ddtdyic
                     hessx(3,ic) = hessx(3,ic) + deddt*dxiczic
     &                                  + d2eddt2*ddtdxic*ddtdzic
                     hessy(1,ic) = hessy(1,ic) + deddt*dxicyic
     &                                  + d2eddt2*ddtdyic*ddtdxic
                     hessy(2,ic) = hessy(2,ic) + deddt*dyicyic
     &                                  + d2eddt2*ddtdyic*ddtdyic
                     hessy(3,ic) = hessy(3,ic) + deddt*dyiczic
     &                                  + d2eddt2*ddtdyic*ddtdzic
                     hessz(1,ic) = hessz(1,ic) + deddt*dxiczic
     &                                  + d2eddt2*ddtdzic*ddtdxic
                     hessz(2,ic) = hessz(2,ic) + deddt*dyiczic
     &                                  + d2eddt2*ddtdzic*ddtdyic
                     hessz(3,ic) = hessz(3,ic) + deddt*dziczic
     &                                  + d2eddt2*ddtdzic*ddtdzic
                     hessx(1,ib) = hessx(1,ib) + deddt*dxibxic
     &                                  + d2eddt2*ddtdxic*ddtdxib
                     hessx(2,ib) = hessx(2,ib) + deddt*dyibxic
     &                                  + d2eddt2*ddtdxic*ddtdyib
                     hessx(3,ib) = hessx(3,ib) + deddt*dzibxic
     &                                  + d2eddt2*ddtdxic*ddtdzib
                     hessy(1,ib) = hessy(1,ib) + deddt*dxibyic
     &                                  + d2eddt2*ddtdyic*ddtdxib
                     hessy(2,ib) = hessy(2,ib) + deddt*dyibyic
     &                                  + d2eddt2*ddtdyic*ddtdyib
                     hessy(3,ib) = hessy(3,ib) + deddt*dzibyic
     &                                  + d2eddt2*ddtdyic*ddtdzib
                     hessz(1,ib) = hessz(1,ib) + deddt*dxibzic
     &                                  + d2eddt2*ddtdzic*ddtdxib
                     hessz(2,ib) = hessz(2,ib) + deddt*dyibzic
     &                                  + d2eddt2*ddtdzic*ddtdyib
                     hessz(3,ib) = hessz(3,ib) + deddt*dzibzic
     &                                  + d2eddt2*ddtdzic*ddtdzib
                     hessx(1,ia) = hessx(1,ia) + deddt*dxiaxic
     &                                  + d2eddt2*ddtdxic*ddtdxia
                     hessx(2,ia) = hessx(2,ia) + deddt*dyiaxic
     &                                  + d2eddt2*ddtdxic*ddtdyia
                     hessx(3,ia) = hessx(3,ia) + deddt*dziaxic
     &                                  + d2eddt2*ddtdxic*ddtdzia
                     hessy(1,ia) = hessy(1,ia) + deddt*dxiayic
     &                                  + d2eddt2*ddtdyic*ddtdxia
                     hessy(2,ia) = hessy(2,ia) + deddt*dyiayic
     &                                  + d2eddt2*ddtdyic*ddtdyia
                     hessy(3,ia) = hessy(3,ia) + deddt*dziayic
     &                                  + d2eddt2*ddtdyic*ddtdzia
                     hessz(1,ia) = hessz(1,ia) + deddt*dxiazic
     &                                  + d2eddt2*ddtdzic*ddtdxia
                     hessz(2,ia) = hessz(2,ia) + deddt*dyiazic
     &                                  + d2eddt2*ddtdzic*ddtdyia
                     hessz(3,ia) = hessz(3,ia) + deddt*dziazic
     &                                  + d2eddt2*ddtdzic*ddtdzia
                  end if
               end if
            end if
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangle2b  --  in-plane bend Hessian; numerical  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangle2b" computes projected in-plane bending first derivatives
c     for a single angle with respect to Cartesian coordinates;
c     used in computation of finite difference second derivatives
c
c
      subroutine eangle2b (i)
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'math.i'
      integer i,ia,ib,ic,id
      real*8 ideal,force
      real*8 dot,cosine,angle
      real*8 dt,dt2,dt3,dt4
      real*8 deddt,terma,termc,term
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xip,yip,zip
      real*8 xap,yap,zap,rap2
      real*8 xcp,ycp,zcp,rcp2
      real*8 xt,yt,zt,rt2,ptrt2
      real*8 xm,ym,zm,rm,delta,delta2
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxip,dedyip,dedzip
      real*8 dpdxia,dpdyia,dpdzia
      real*8 dpdxic,dpdyic,dpdzic
c
c
c     set the atom numbers and parameters for this angle
c
      ia = iang(1,i)
      ib = iang(2,i)
      ic = iang(3,i)
      id = iang(4,i)
      ideal = anat(i)
      force = acon(i)
c
c     get the coordinates of the atoms in the angle
c
      xia = x(ia)
      yia = y(ia)
      zia = z(ia)
      xib = x(ib)
      yib = y(ib)
      zib = z(ib)
      xic = x(ic)
      yic = y(ic)
      zic = z(ic)
      xid = x(id)
      yid = y(id)
      zid = z(id)
c
c     zero out the first derivative components
c
      dea(1,ia) = 0.0d0
      dea(2,ia) = 0.0d0
      dea(3,ia) = 0.0d0
      dea(1,ib) = 0.0d0
      dea(2,ib) = 0.0d0
      dea(3,ib) = 0.0d0
      dea(1,ic) = 0.0d0
      dea(2,ic) = 0.0d0
      dea(3,ic) = 0.0d0
      dea(1,id) = 0.0d0
      dea(2,id) = 0.0d0
      dea(3,id) = 0.0d0
c
c     compute the value of the projected in-plane angle
c
      xad = xia - xid
      yad = yia - yid
      zad = zia - zid
      xbd = xib - xid
      ybd = yib - yid
      zbd = zib - zid
      xcd = xic - xid
      ycd = yic - yid
      zcd = zic - zid
      xt = yad*zcd - zad*ycd
      yt = zad*xcd - xad*zcd
      zt = xad*ycd - yad*xcd
      rt2 = xt*xt + yt*yt + zt*zt
      delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
      xip = xib + xt*delta
      yip = yib + yt*delta
      zip = zib + zt*delta
      xap = xia - xip
      yap = yia - yip
      zap = zia - zip
      rap2 = xap*xap + yap*yap + zap*zap
      xcp = xic - xip
      ycp = yic - yip
      zcp = zic - zip
      rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
      xm = ycp*zap - zcp*yap
      ym = zcp*xap - xcp*zap
      zm = xcp*yap - ycp*xap
      rm = sqrt(xm*xm + ym*ym + zm*zm)
      if (rm .ne. 0.0d0) then
         dot = xap*xcp + yap*ycp + zap*zcp
         cosine = dot / sqrt(rap2*rcp2)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
c
c     get the in-plane bending master chain rule terms
c
         dt = angle - ideal
         dt2 = dt * dt
         dt3 = dt2 * dt
         dt4 = dt2 * dt2
         deddt = angunit * force * dt * radian
     &             * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2
     &                  + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)
c
c     chain rule terms for first derivative components
c
         terma = -deddt / (rap2*rm)
         termc = deddt / (rcp2*rm)
         dedxia = terma * (yap*zm-zap*ym)
         dedyia = terma * (zap*xm-xap*zm)
         dedzia = terma * (xap*ym-yap*xm)
         dedxic = termc * (ycp*zm-zcp*ym)
         dedyic = termc * (zcp*xm-xcp*zm)
         dedzic = termc * (xcp*ym-ycp*xm)
         dedxip = -dedxia - dedxic
         dedyip = -dedyia - dedyic
         dedzip = -dedzia - dedzic
c
c     chain rule components for the projection of the central atom
c
         delta2 = 2.0d0 * delta
         ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2
         term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd)
         dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2
         term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd)
         dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2
         term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd)
         dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2
         term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad)
         dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2
         term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad)
         dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2
         term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad)
         dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2
c
c     compute derivative components for this interaction
c
         dedxia = dedxia + dpdxia
         dedyia = dedyia + dpdyia
         dedzia = dedzia + dpdzia
         dedxib = dedxip
         dedyib = dedyip
         dedzib = dedzip
         dedxic = dedxic + dpdxic
         dedyic = dedyic + dpdyic
         dedzic = dedzic + dpdzic
         dedxid = -dedxia - dedxib - dedxic
         dedyid = -dedyia - dedyib - dedyic
         dedzid = -dedzia - dedzib - dedzic
c
c     set the in-plane angle bending first derivatives
c
         dea(1,ia) = dedxia
         dea(2,ia) = dedyia
         dea(3,ia) = dedzia
         dea(1,ib) = dedxib
         dea(2,ib) = dedyib
         dea(3,ib) = dedzib
         dea(1,ic) = dedxic
         dea(2,ic) = dedyic
         dea(3,ic) = dedzic
         dea(1,id) = dedxid
         dea(2,id) = dedyid
         dea(3,id) = dedzid
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
c     ##  subroutine eangle3  --  angle bending energy & analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eangle3" calculates the angle bending potential energy,
c     also partitions the energy among the atoms; for trigonal
c     centers, projected in-plane angles as per the formalism
c     of Wilson-Decius-Cross are optionally used
c
c
      subroutine eangle3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'usage.i'
      integer i,ia,ib,ic,id
      real*8 e,ideal,force
      real*8 dot,cosine
      real*8 angle,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab,rab2
      real*8 xcb,ycb,zcb,rcb2
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xip,yip,zip
      real*8 xap,yap,zap,rap2
      real*8 xcp,ycp,zcp,rcp2
      real*8 xt,yt,zt,rt2,delta
      logical header,huge,proceed
c
c
c     zero out the angle bending energy and partitioning terms
c
      nea = 0
      ea = 0.0d0
      do i = 1, n
         aea(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         ideal = anat(i)
         force = acon(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
c
c     compute the value of the bond angle
c
            if (.not.angin(i)) then
               xab = xia - xib
               yab = yia - yib
               zab = zia - zib
               rab2 = xab*xab + yab*yab + zab*zab
               xcb = xic - xib
               ycb = yic - yib
               zcb = zic - zib
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
                  dot = xab*xcb + yab*ycb + zab*zcb
                  cosine = dot / sqrt(rab2*rcb2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
c
c     find the bond angle bending energy
c
                  dt = angle - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the total bond angle bending energy
c
                  nea = nea + 1
                  ea = ea + e
                  aea(ib) = aea(ib) + e
c
c     print a warning if the energy of this angle is large
c
                  huge = (e .gt. 5.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Angle Bending',
     &                             ' Interactions :',
     &                          //,' Type',15x,'Atom Names',16x,
     &                             'Ideal',4x,'Actual',6x,'Energy',/)
                     end if
                     write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                                name(ic),ideal,angle,e
   20                format (' Angle    ',i5,'-',a3,1x,i5,'-',a3,
     &                         1x,i5,'-',a3,2x,2f10.4,f12.4)
                  end if
               end if
c
c     compute the value of the projected in-plane angle
c
            else
               xid = x(id)
               yid = y(id)
               zid = z(id)
               xad = xia - xid
               yad = yia - yid
               zad = zia - zid
               xbd = xib - xid
               ybd = yib - yid
               zbd = zib - zid
               xcd = xic - xid
               ycd = yic - yid
               zcd = zic - zid
               xt = yad*zcd - zad*ycd
               yt = zad*xcd - xad*zcd
               zt = xad*ycd - yad*xcd
               rt2 = xt*xt + yt*yt + zt*zt
               delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
               xip = xib + xt*delta
               yip = yib + yt*delta
               zip = zib + zt*delta
               xap = xia - xip
               yap = yia - yip
               zap = zia - zip
               rap2 = xap*xap + yap*yap + zap*zap
               xcp = xic - xip
               ycp = yic - yip
               zcp = zic - zip
               rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
               if (rap2.ne.0.0d0 .and. rcp2.ne.0.0d0) then
                  dot = xap*xcp + yap*ycp + zap*zcp
                  cosine = dot / sqrt(rap2*rcp2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
c
c     find the in-plane angle bending energy
c
                  dt = angle - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the total bond angle bending energy
c
                  nea = nea + 1
                  ea = ea + e
                  aea(ib) = aea(ib) + e
c
c     print a warning if the energy of this angle is large
c
                  huge = (e .gt. 5.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,30)
   30                   format (/,' Individual Angle Bending',
     &                             ' Interactions :',
     &                          //,' Type',15x,'Atom Names',16x,
     &                             'Ideal',4x,'Actual',6x,'Energy',/)
                     end if
                     write (iout,40)  ia,name(ia),ib,name(ib),ic,
     &                                name(ic),ideal,angle,e
   40                format (' Angle-IP ',i5,'-',a3,1x,i5,'-',a3,
     &                          1x,i5,'-',a3,2x,2f10.4,f12.4)
                  end if
               end if
            end if
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
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine ebond  --  bond stretch potential energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "ebond" calculates the bond stretching energy
c
c
      subroutine ebond
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bndpot.i'
      include 'bond.i'
      include 'energi.i'
      include 'group.i'
      include 'usage.i'
      integer i,ia,ib
      real*8 e,ideal,force
      real*8 expterm,bde
      real*8 dt,dt2,fgrp
      real*8 xab,yab,zab,rab
      logical proceed
c
c
c     zero out the bond stretching energy
c
      eb = 0.0d0
c
c     calculate the bond stretching energy term
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,ia,ib,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
c
c     Morse potential; energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0d0*dt)
               bde = 0.25d0 * bndunit * force
               e = bde * (1.0d0-expterm)**2
c
c     Taylor series expansion of Morse potential through
c     the fourth power of the bond length deviation
c
            else
               dt2 = dt * dt
               e = bndunit * force * dt2 * (1.0d0+cbnd*dt+qbnd*dt2)
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total bond stretching energy
c
            eb = eb + e
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ebond1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ebond1" calculates the bond stretching energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine ebond1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bndpot.i'
      include 'bond.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ia,ib
      real*8 e,ideal,force
      real*8 expterm,bde,fgrp
      real*8 dt,dt2,deddt
      real*8 de,dedx,dedy,dedz
      real*8 xab,yab,zab,rab
      logical proceed
c
c
c     zero out the bond energy and first derivatives
c
      eb = 0.0d0
      do i = 1, n
         deb(1,i) = 0.0d0
         deb(2,i) = 0.0d0
         deb(3,i) = 0.0d0
      end do
c
c     calculate the bond stretch energy and first derivatives
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,ia,ib,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
c
c     Morse potential; energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0d0*dt)
               bde = 0.25d0 * bndunit * force
               e = bde * (1.0d0-expterm)**2
               deddt = 4.0d0 * bde * (1.0d0-expterm) * expterm
c
c     Taylor series expansion of Morse potential through
c     the fourth power of the bond length deviation
c
            else
               dt2 = dt * dt
               e = bndunit * force * dt2 * (1.0d0+cbnd*dt+qbnd*dt2)
               deddt = 2.0d0 * bndunit * force * dt
     &                    * (1.0d0+1.5d0*cbnd*dt+2.0d0*qbnd*dt2)
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               deddt = deddt * fgrp
            end if
c
c     compute chain rule terms needed for derivatives
c
            if (rab .eq. 0.0d0) then
               de = 0.0d0
            else
               de = deddt / rab
            end if
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total bond energy and first derivatives
c
            eb = eb + e
            deb(1,ia) = deb(1,ia) + dedx
            deb(2,ia) = deb(2,ia) + dedy
            deb(3,ia) = deb(3,ia) + dedz
            deb(1,ib) = deb(1,ib) - dedx
            deb(2,ib) = deb(2,ib) - dedy
            deb(3,ib) = deb(3,ib) - dedz
c
c     increment the virial for use in pressure computation
c
            if (isobaric) then
               virx = virx + xab*dedx
               viry = viry + yab*dedy
               virz = virz + zab*dedz
            end if
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ebond2  --  atom-by-atom bond stretch Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ebond2" calculates second derivatives of the bond
c     stretching energy for a single atom at a time
c
c
      subroutine ebond2 (i)
      implicit none
      include 'sizes.i'
      include 'atmlst.i'
      include 'atoms.i'
      include 'bndpot.i'
      include 'bond.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      integer i,j,k,ia,ib
      real*8 ideal,force,fgrp
      real*8 dt,dt2,expterm,bde
      real*8 xab,yab,zab,rab,rab2
      real*8 term,termx,termy,termz
      real*8 de,deddt,d2eddt2,d2e(3,3)
      logical proceed
c
c
c     compute the Hessian elements of the bond stretch energy
c
      ia = i
      do k = 1, n12(ia)
         j = bndlist(k,ia)
         if (ibnd(1,j) .eq. ia) then
            ib = ibnd(2,j)
         else
            ib = ibnd(1,j)
         end if
         ideal = bl(j)
         force = bk(j)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,ia,ib,0,0,0)
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            rab2 = xab*xab + yab*yab + zab*zab
            rab = sqrt(rab2)
            dt = rab - ideal
c
c     Morse potential; energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0d0*dt)
               bde = 0.25d0 * bndunit * force
               deddt = 4.0d0 * bde * (1.0d0-expterm) * expterm
               d2eddt2 = -8.0d0 * bde * (1.0d0-2.0d0*expterm) * expterm
c
c     Taylor series expansion of Morse potential through
c     the fourth power of the bond length deviation
c
            else
               dt2 = dt * dt
               deddt = 2.0d0 * bndunit * force * dt
     &                    * (1.0d0+1.5d0*cbnd*dt+2.0d0*qbnd*dt2)
               d2eddt2 = 2.0d0 * bndunit * force
     &                      * (1.0d0+3.0d0*cbnd*dt+6.0d0*qbnd*dt2)
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               deddt = deddt * fgrp
               d2eddt2 = d2eddt2 * fgrp
            end if
c
c     set the chain rule terms for the Hessian elements
c
            if (rab2 .eq. 0.0d0) then
               de = 0.0d0
               term = 0.0d0
            else
               de = deddt / rab
               term = (d2eddt2-de) / rab2
            end if
            termx = term * xab
            termy = term * yab
            termz = term * zab
            d2e(1,1) = termx*xab + de
            d2e(1,2) = termx*yab
            d2e(1,3) = termx*zab
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yab + de
            d2e(2,3) = termy*zab
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + d2e(1,j)
               hessy(j,ia) = hessy(j,ia) + d2e(2,j)
               hessz(j,ia) = hessz(j,ia) + d2e(3,j)
               hessx(j,ib) = hessx(j,ib) - d2e(1,j)
               hessy(j,ib) = hessy(j,ib) - d2e(2,j)
               hessz(j,ib) = hessz(j,ib) - d2e(3,j)
            end do
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ebond3  --  bond stretch energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ebond3" calculates the bond stretching energy; also
c     partitions the energy among the atoms
c
c
      subroutine ebond3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bndpot.i'
      include 'bond.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'usage.i'
      integer i,ia,ib
      real*8 e,ideal,force
      real*8 expterm,bde
      real*8 dt,dt2,fgrp
      real*8 xab,yab,zab,rab
      logical header,huge,proceed
c
c
c     zero out the bond energy and partitioning terms
c
      neb = 0
      eb = 0.0d0
      do i = 1, n
         aeb(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the bond stretching energy term
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,ia,ib,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
c
c     Morse potential; energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0d0*dt)
               bde = 0.25d0 * bndunit * force
               e = bde * (1.0d0-expterm)**2
c
c     Taylor series expansion of Morse potential through
c     the fourth power of the bond length deviation
c
            else
               dt2 = dt * dt
               e = bndunit * force * dt2 * (1.0d0+cbnd*dt+qbnd*dt2)
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total bond energy and partition between atoms
c
            neb = neb + 1
            eb = eb + e
            aeb(ia) = aeb(ia) + 0.5d0*e
            aeb(ib) = aeb(ib) + 0.5d0*e
c
c     print a warning if the energy of this bond is large
c
            huge = (e .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Bond Stretching',
     &                       ' Interactions :',
     &                    //,' Type',11x,'Atom Names',20x,'Ideal',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,20)  ia,name(ia),ib,name(ib),ideal,rab,e
   20          format (' Bond     ',i5,'-',a3,1x,i5,'-',a3,
     &                    12x,2f10.4,f12.4)
            end if
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ebuck  --  Buckingham van der Waals energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ebuck" calculates the van der Waals interaction energy
c     using the Buckingham exp-6 formula
c
c
      subroutine ebuck
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,ii,kk,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rv,eps,rdn
      real*8 p,p2,p6,p12,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 rik,rik2,rik3,rik4,rik5,taper
      real*8 expcut,expcut2,expterm,expmerge
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical proceed,iuse
c
c
c     zero out the van der Waals energy contributions
c
      ev = 0.0d0
      e14 = 0.0d0
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find the van der Waals energy via double loop search
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i)  eps = eps / vdwscale
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - cterm*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (skip(k) .eq. -i) then
                     e14 = e14 + e
                  else
                     ev = ev + e
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  kt = class(k)
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
                  if (rik2 .le. off2) then
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        expterm = aterm * exp(-bterm/p)
                        e = eps * (expterm - cterm*p6)
                     else
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik = sqrt(rik2)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy component;
c     interaction of an atom with its own image counts half
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                  end if
               end do
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ebuck1  --  Buckingham energy & derivatives  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ebuck1" calculates the van der Waals energy and its first
c     derivatives with respect to Cartesian coordinates using
c     the Buckingham exp-6 formula
c
c
      subroutine ebuck1
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,k,ii,kk,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rv,eps,rdn
      real*8 p,p2,p6,p12,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv,rvterm
      real*8 dedx,dedy,dedz,de
      real*8 rik,rik2,rik3,rik4,rik5,taper,dtaper
      real*8 expcut,expcut2,expterm,expmerge
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical proceed,iuse
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      e14 = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
         de14(1,i) = 0.0d0
         de14(2,i) = 0.0d0
         de14(3,i) = 0.0d0
      end do
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find van der Waals energy and derivatives via double loop
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute energy and derivatives for this interaction
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i)  eps = eps / vdwscale
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bterm / rv
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - cterm*p6)
                     de = eps * (rvterm*expterm+6.0d0*cterm*p6/rik)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  if (skip(k) .eq. -i) then
                     e14 = e14 + e
                     if (i .eq. iv) then
                        de14(1,i) = de14(1,i) + dedx
                        de14(2,i) = de14(2,i) + dedy
                        de14(3,i) = de14(3,i) + dedz
                     else
                        de14(1,i) = de14(1,i) + dedx*redi
                        de14(2,i) = de14(2,i) + dedy*redi
                        de14(3,i) = de14(3,i) + dedz*redi
                        de14(1,iv) = de14(1,iv) + dedx*rediv
                        de14(2,iv) = de14(2,iv) + dedy*rediv
                        de14(3,iv) = de14(3,iv) + dedz*rediv
                     end if
                     if (k .eq. kv) then
                        de14(1,k) = de14(1,k) - dedx
                        de14(2,k) = de14(2,k) - dedy
                        de14(3,k) = de14(3,k) - dedz
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        de14(1,k) = de14(1,k) - dedx*redk
                        de14(2,k) = de14(2,k) - dedy*redk
                        de14(3,k) = de14(3,k) - dedz*redk
                        de14(1,kv) = de14(1,kv) - dedx*redkv
                        de14(2,kv) = de14(2,kv) - dedy*redkv
                        de14(3,kv) = de14(3,kv) - dedz*redkv
                     end if
                  else
                     ev = ev + e
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                     end if
                     if (k .eq. kv) then
                        dev(1,k) = dev(1,k) - dedx
                        dev(2,k) = dev(2,k) - dedy
                        dev(3,k) = dev(3,k) - dedz
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        dev(1,k) = dev(1,k) - dedx*redk
                        dev(2,k) = dev(2,k) - dedy*redk
                        dev(3,k) = dev(3,k) - dedz*redk
                        dev(1,kv) = dev(1,kv) - dedx*redkv
                        dev(2,kv) = dev(2,kv) - dedy*redkv
                        dev(3,kv) = dev(3,kv) - dedz*redkv
                     end if
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx + xr*dedx
                     viry = viry + yr*dedy
                     virz = virz + zr*dedz
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  kt = class(k)
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     compute energy and derivatives for this interaction
c
                  if (rik2 .le. off2) then
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     rik = sqrt(rik2)
                     if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        rvterm = -bterm / rv
                        expterm = aterm * exp(-bterm/p)
                        e = eps * (expterm - cterm*p6)
                        de = eps * (rvterm*expterm+6.0d0*cterm*p6/rik)
                     else
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                        de = -12.0d0 * e / rik
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                              + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        de = e*dtaper + de*taper
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                     end if
c
c     find the chain rule terms for derivative components
c
                     de = de / rik
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                     end if
                     if (i .ne. k) then
                        if (k .eq. kv) then
                           dev(1,k) = dev(1,k) - dedx
                           dev(2,k) = dev(2,k) - dedy
                           dev(3,k) = dev(3,k) - dedz
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           dev(1,k) = dev(1,k) - dedx*redk
                           dev(2,k) = dev(2,k) - dedy*redk
                           dev(3,k) = dev(3,k) - dedz*redk
                           dev(1,kv) = dev(1,kv) - dedx*redkv
                           dev(2,kv) = dev(2,kv) - dedy*redkv
                           dev(3,kv) = dev(3,kv) - dedz*redkv
                        end if
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*dedx
                        viry = viry + yr*dedy
                        virz = virz + zr*dedz
                     end if
                  end if
               end do
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ebuck2  --  atom-by-atom Buckingham Hessian  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ebuck2" calculates the van der Waals second derivatives for
c     a single atom at a time using the Buckingham exp-6 formula
c
c
      subroutine ebuck2 (iatom,xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'shunt.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer iatom,i,j,k,ii,kk,iv,kv
      integer it,kt,skip(maxatm)
      integer nuse,use(5),jcell
      real*8 e,de,d2e,fgrp
      real*8 p,p2,p6,p12,eps,rv
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 redi2,rediv2,rediiv
      real*8 redik,redivk,redikv,redivkv
      real*8 rik,rik2,rik3,rik4,rik5
      real*8 taper,dtaper,d2taper
      real*8 d2edx,d2edy,d2edz,term(3,3)
      real*8 expcut,expcut2
      real*8 expterm,expmerge
      real*8 rvterm,rvterm2
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical proceed
c
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     check to see if the atom of interest is a vdw site
c
      nuse = 0
      do k = 1, nvdw
         if (ivdw(k) .eq. iatom) then
            nuse = nuse + 1
            use(nuse) = iatom
            goto 10
         end if
      end do
      return
   10 continue
c
c     determine the atoms involved via reduction factors
c
      nuse = 1
      use(nuse) = iatom
      do k = 1, n12(iatom)
         i = i12(k,iatom)
         if (ired(i) .eq. iatom) then
            nuse = nuse + 1
            use(nuse) = i
         end if
      end do
c
c     find van der Waals Hessian elements for involved atoms
c
      do ii = 1, nuse
         i = use(ii)
         iv = ired(i)
         redi = kred(i)
         if (i .ne. iv) then
            rediv = 1.0d0 - redi
            redi2 = redi * redi
            rediv2 = rediv * rediv
            rediiv = redi * rediv
         end if
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         skip(i) = i
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         do kk = 1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute Hessian elements for this interaction
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i)  eps = eps / vdwscale
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bterm / rv
                     rvterm2 = rvterm * rvterm
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - cterm*p6)
                     de = eps * (rvterm*expterm+6.0d0*cterm*p6/rik)
                     d2e = eps * (rvterm2*expterm-42.0d0*cterm*p6/rik2)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                     d2e = 156.0d0 * e / rik2
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     e = eps * (p12 - 2.0d0 * p6)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     d2taper = 20.0d0*c5*rik3 + 12.0d0*c4*rik2
     &                            + 6.0d0*c3*rik + 2.0d0*c2
                     d2e = e*d2taper + 2.0d0*de*dtaper + d2e*taper
                     de = e*dtaper + de*taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     de = de * fgrp
                     d2e = d2e * fgrp
                  end if
c
c     get chain rule terms for van der Waals Hessian elements
c
                  de = de / rik
                  d2e = (d2e-de) / rik2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  if (i .eq. iatom) then
                     if (i.eq.iv .and. k.eq.kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)
                           hessy(j,i) = hessy(j,i) + term(2,j)
                           hessz(j,i) = hessz(j,i) + term(3,j)
                           hessx(j,k) = hessx(j,k) - term(1,j)
                           hessy(j,k) = hessy(j,k) - term(2,j)
                           hessz(j,k) = hessz(j,k) - term(3,j)
                        end do
                     else if (k .eq. kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                           hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                           hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                           hessx(j,k) = hessx(j,k) - term(1,j)*redi
                           hessy(j,k) = hessy(j,k) - term(2,j)*redi
                           hessz(j,k) = hessz(j,k) - term(3,j)*redi
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                        end do
                     else if (i .eq. iv) then
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)
                           hessy(j,i) = hessy(j,i) + term(2,j)
                           hessz(j,i) = hessz(j,i) + term(3,j)
                           hessx(j,k) = hessx(j,k) - term(1,j)*redk
                           hessy(j,k) = hessy(j,k) - term(2,j)*redk
                           hessz(j,k) = hessz(j,k) - term(3,j)*redk
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redkv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redkv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redkv
                        end do
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        redik = redi * redk
                        redikv = redi * redkv
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                           hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                           hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                           hessx(j,k) = hessx(j,k) - term(1,j)*redik
                           hessy(j,k) = hessy(j,k) - term(2,j)*redik
                           hessz(j,k) = hessz(j,k) - term(3,j)*redik
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediiv
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediiv
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediiv
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redikv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redikv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redikv
                        end do
                     end if
                  else if (iv .eq. iatom) then
                     if (k .eq. kv) then
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                           hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                           hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                           hessx(j,k) = hessx(j,k) - term(1,j)*rediv
                           hessy(j,k) = hessy(j,k) - term(2,j)*rediv
                           hessz(j,k) = hessz(j,k) - term(3,j)*rediv
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                        end do
                     else
                        redk = kred(k)
                        redkv = 1.0d0 - redk
                        redivk = rediv * redk
                        redivkv = rediv * redkv
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                           hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                           hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                           hessx(j,k) = hessx(j,k) - term(1,j)*redivk
                           hessy(j,k) = hessy(j,k) - term(2,j)*redivk
                           hessz(j,k) = hessz(j,k) - term(3,j)*redivk
                           hessx(j,iv) = hessx(j,iv) + term(1,j)*rediv2
                           hessy(j,iv) = hessy(j,iv) + term(2,j)*rediv2
                           hessz(j,iv) = hessz(j,iv) + term(3,j)*rediv2
                           hessx(j,kv) = hessx(j,kv) - term(1,j)*redivkv
                           hessy(j,kv) = hessy(j,kv) - term(2,j)*redivkv
                           hessz(j,kv) = hessz(j,kv) - term(3,j)*redivkv
                        end do
                     end if
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nuse
         i = use(ii)
         iv = ired(i)
         redi = kred(i)
         if (i .ne. iv) then
            rediv = 1.0d0 - redi
            redi2 = redi * redi
            rediv2 = rediv * rediv
            rediiv = redi * rediv
         end if
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do kk = 1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do jcell = 1, ncell
                  kt = class(k)
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,jcell)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     compute Hessian elements for this interaction
c
                  if (rik2 .le. off2) then
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     rik = sqrt(rik2)
                     if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        rvterm = -bterm / rv
                        rvterm2 = rvterm * rvterm
                        expterm = aterm * exp(-bterm/p)
                        e = eps * (expterm - cterm*p6)
                        de = eps * (rvterm*expterm+6.0d0*cterm*p6/rik)
                        d2e = eps * (rvterm2*expterm
     &                                  -42.0d0*cterm*p6/rik2)
                     else
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                        de = -12.0d0 * e / rik
                        d2e = 156.0d0 * e / rik2
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        e = eps * (p12 - 2.0d0 * p6)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        d2taper = 20.0d0*c5*rik3 + 12.0d0*c4*rik2
     &                             + 6.0d0*c3*rik + 2.0d0*c2
                        d2e = e*d2taper + 2.0d0*de*dtaper + d2e*taper
                        de = e*dtaper + de*taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        de = de * fgrp
                        d2e = d2e * fgrp
                     end if
c
c     get chain rule terms for van der Waals Hessian elements
c
                     de = de / rik
                     d2e = (d2e-de) / rik2
                     d2edx = d2e * xr
                     d2edy = d2e * yr
                     d2edz = d2e * zr
                     term(1,1) = d2edx*xr + de
                     term(1,2) = d2edx*yr
                     term(1,3) = d2edx*zr
                     term(2,1) = term(1,2)
                     term(2,2) = d2edy*yr + de
                     term(2,3) = d2edy*zr
                     term(3,1) = term(1,3)
                     term(3,2) = term(2,3)
                     term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                     if (i .eq. iatom) then
                        if (i.eq.iv .and. k.eq.kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)
                              hessy(j,i) = hessy(j,i) + term(2,j)
                              hessz(j,i) = hessz(j,i) + term(3,j)
                              hessx(j,k) = hessx(j,k) - term(1,j)
                              hessy(j,k) = hessy(j,k) - term(2,j)
                              hessz(j,k) = hessz(j,k) - term(3,j)
                           end do
                        else if (k .eq. kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                              hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                              hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                              hessx(j,k) = hessx(j,k) - term(1,j)*redi
                              hessy(j,k) = hessy(j,k) - term(2,j)*redi
                              hessz(j,k) = hessz(j,k) - term(3,j)*redi
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediiv
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediiv
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediiv
                           end do
                        else if (i .eq. iv) then
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)
                              hessy(j,i) = hessy(j,i) + term(2,j)
                              hessz(j,i) = hessz(j,i) + term(3,j)
                              hessx(j,k) = hessx(j,k) - term(1,j)*redk
                              hessy(j,k) = hessy(j,k) - term(2,j)*redk
                              hessz(j,k) = hessz(j,k) - term(3,j)*redk
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redkv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redkv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redkv
                           end do
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           redik = redi * redk
                           redikv = redi * redkv
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*redi2
                              hessy(j,i) = hessy(j,i) + term(2,j)*redi2
                              hessz(j,i) = hessz(j,i) + term(3,j)*redi2
                              hessx(j,k) = hessx(j,k) - term(1,j)*redik
                              hessy(j,k) = hessy(j,k) - term(2,j)*redik
                              hessz(j,k) = hessz(j,k) - term(3,j)*redik
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediiv
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediiv
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediiv
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redikv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redikv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redikv
                           end do
                        end if
                     else if (iv .eq. iatom) then
                        if (k .eq. kv) then
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                              hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                              hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                              hessx(j,k) = hessx(j,k) - term(1,j)*rediv
                              hessy(j,k) = hessy(j,k) - term(2,j)*rediv
                              hessz(j,k) = hessz(j,k) - term(3,j)*rediv
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediv2
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediv2
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediv2
                           end do
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           redivk = rediv * redk
                           redivkv = rediv * redkv
                           do j = 1, 3
                              hessx(j,i) = hessx(j,i) + term(1,j)*rediiv
                              hessy(j,i) = hessy(j,i) + term(2,j)*rediiv
                              hessz(j,i) = hessz(j,i) + term(3,j)*rediiv
                              hessx(j,k) = hessx(j,k) - term(1,j)*redivk
                              hessy(j,k) = hessy(j,k) - term(2,j)*redivk
                              hessz(j,k) = hessz(j,k) - term(3,j)*redivk
                              hessx(j,iv) = hessx(j,iv)
     &                                         + term(1,j)*rediv2
                              hessy(j,iv) = hessy(j,iv)
     &                                         + term(2,j)*rediv2
                              hessz(j,iv) = hessz(j,iv)
     &                                         + term(3,j)*rediv2
                              hessx(j,kv) = hessx(j,kv)
     &                                         - term(1,j)*redivkv
                              hessy(j,kv) = hessy(j,kv)
     &                                         - term(2,j)*redivkv
                              hessz(j,kv) = hessz(j,kv)
     &                                         - term(3,j)*redivkv
                           end do
                        end if
                     end if
                  end if
               end do
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
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine ebuck3  --  Buckingham energy & analysis  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "ebuck3" calculates the van der Waals interaction energy
c     using the Buckingham exp-6 formula and also partitions
c     the energy among the atoms
c
c
      subroutine ebuck3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,ii,kk,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rv,eps,rdn
      real*8 p,p2,p6,p12,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 rik,rik2,rik3,rik4,rik5,taper
      real*8 expcut,expcut2,expterm,expmerge
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical header,huge,proceed,iuse
c
c
c     zero out the van der Waals energy and partitioning terms
c
      nev = 0
      ne14 = 0
      ev = 0.0d0
      e14 = 0.0d0
      do i = 1, n
         aev(i) = 0.0d0
         ae14(i) = 0.0d0
      end do
      header = .true.
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find the van der Waals energy via double loop search
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               if (use_image)  call image (xr,yr,zr,0)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i)  eps = eps / vdwscale
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - cterm*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (skip(k) .eq. -i) then
                     ne14 = ne14 + 1
                     e14 = e14 + e
                     ae14(i) = ae14(i) + 0.5d0*e
                     ae14(k) = ae14(k) + 0.5d0*e
                  else
                     nev = nev + 1
                     ev = ev + e
                     aev(i) = aev(i) + 0.5d0*e
                     aev(k) = aev(k) + 0.5d0*e
                  end if
c
c     print a warning if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual van der Waals',
     &                             ' Interactions :',
     &                          //,' Type',11x,'Atom Names',
     &                             18x,'Minimum',4x,'Actual',
     &                             6x,'Energy',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),
     &                                rv,sqrt(rik2),e
   20                format (' VDW-Buck ',i5,'-',a3,1x,i5,'-',a3,
     &                          12x,2f10.4,f12.4)
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kv))
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  kt = class(k)
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call image (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
                  if (rik2 .le. off2) then
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     p2 = (rv*rv) / rik2
                     p6 = p2 * p2 * p2
                     if (p2 .le. expcut2) then
                        p = sqrt(p2)
                        expterm = aterm * exp(-bterm/p)
                        e = eps * (expterm - cterm*p6)
                     else
                        p12 = p6 * p6
                        e = expmerge * eps * p12
                     end if
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik = sqrt(rik2)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                                + c2*rik2 + c1*rik + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy component
c
                     nev = nev + 1
                     if (i .eq. k) then
                        ev = ev + 0.5d0*e
                        aev(i) = aev(i) + 0.5d0*e
                     else
                        ev = ev + e
                        aev(i) = aev(i) + 0.5d0*e
                        aev(k) = aev(k) + 0.5d0*e
                     end if
c
c     print a warning if the energy of this interaction is large
c
                     huge = (e .gt. 10.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual van der Waals',
     &                                ' Interactions :',
     &                             //,' Type',11x,'Atom Names',
     &                                18x,'Minimum',4x,'Actual',
     &                                6x,'Energy',/)
                        end if
                        write (iout,40)  i,name(i),k,name(k),
     &                                   rv,sqrt(rik2),e
   40                   format (' VDW-Buck ',i5,'-',a3,1x,i5,'-',a3,
     &                             '   (XTAL)   ',2f10.4,f12.4)
                     end if
                  end if
               end do
            end if
         end do
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ebuck4  --  Buckingham van der Waals energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ebuck4" calculates the Buckingham van der Waals interaction
c     energy using the method of lights to locate neighboring atoms
c
c
      subroutine ebuck4
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'iounit.i'
      include 'light.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,ii,kk,iv,it,kt
      integer kgy,kgz,start,stop
      integer kskip,skip(maxatm)
      integer kmap,kvmap,map(maxlight)
      real*8 e,rv,eps,rdn
      real*8 p,p2,p6,p12,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 rik,rik2,rik3,rik4,rik5,taper
      real*8 expcut,expcut2,expterm,expmerge
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      real*8 xsort(maxlight),ysort(maxlight),zsort(maxlight)
      logical proceed,iuse,repeat
c
c
c     zero out the van der Waals energy contributions
c
      ev = 0.0d0
      e14 = 0.0d0
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do j = 1, nvdw
         i = ivdw(j)
         iv = ired(i)
         rdn = kred(i)
         xred(j) = rdn*(x(i)-x(iv)) + x(iv)
         yred(j) = rdn*(y(i)-y(iv)) + y(iv)
         zred(j) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nvdw,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         it = class(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            if (kk .le. nvdw) then
               kskip = skip(ivdw(kk))
               if (kskip .eq. i)  goto 20
            else
               kskip = 0
            end if
            kmap = ivdw(map(kk))
            kvmap = ired(kmap)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,kmap,0,0,0)
            if (proceed)  proceed = (iuse .or. use(kmap)
     &                                .or. use(kvmap))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(kmap)
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_image) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (kskip .eq. -i)  eps = eps / vdwscale
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - cterm*p6)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  if (kskip .eq. -i) then
                     e14 = e14 + e
                  else
                     ev = ev + e
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ebuck5  --  Buckingham energy & derivatives  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ebuck5" calculates the Buckingham van der Waals interaction
c     energy and its first derivatives using the method of lights
c     to locate neighboring atoms
c
c
      subroutine ebuck5
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'iounit.i'
      include 'light.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,ii,kk,iv,it,kt
      integer kgy,kgz,start,stop
      integer kskip,skip(maxatm)
      integer kmap,kvmap,map(maxlight)
      real*8 e,rv,eps,rdn
      real*8 p,p2,p6,p12,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv,rvterm
      real*8 dedx,dedy,dedz,de
      real*8 rik,rik2,rik3,rik4,rik5,taper,dtaper
      real*8 expcut,expcut2,expterm,expmerge
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      real*8 xsort(maxlight),ysort(maxlight),zsort(maxlight)
      logical proceed,iuse,repeat
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      e14 = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
         de14(1,i) = 0.0d0
         de14(2,i) = 0.0d0
         de14(3,i) = 0.0d0
      end do
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     switch from exponential to R^12 at very short range
c
      expcut = 2.0d0
      expcut2 = expcut * expcut
      expmerge = (aterm*exp(-bterm/expcut) - cterm*(expcut**6))
     &                               / (expcut**12)
c
c     calculate the "reduced" atomic coordinates
c
      do j = 1, nvdw
         i = ivdw(j)
         iv = ired(i)
         rdn = kred(i)
         xred(j) = rdn*(x(i)-x(iv)) + x(iv)
         yred(j) = rdn*(y(i)-y(iv)) + y(iv)
         zred(j) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nvdw,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         iuse = (use(i) .or. use(iv))
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = class(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * vdw12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * vdw13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * vdw14use
         end do
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            if (kk .le. nvdw) then
               kskip = skip(ivdw(kk))
               if (kskip .eq. i)  goto 20
            else
               kskip = 0
            end if
            kmap = ivdw(map(kk))
            kvmap = ired(kmap)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,kmap,0,0,0)
            if (proceed)  proceed = (iuse .or. use(kmap)
     &                                .or. use(kvmap))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = class(kmap)
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_image) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute energy and derivatives for this interaction
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (kskip .eq. -i)  eps = eps / vdwscale
                  p2 = (rv*rv) / rik2
                  p6 = p2 * p2 * p2
                  rik = sqrt(rik2)
                  if (p2 .le. expcut2) then
                     p = sqrt(p2)
                     rvterm = -bterm / rv
                     expterm = aterm * exp(-bterm/p)
                     e = eps * (expterm - cterm*p6)
                     de = eps * (rvterm*expterm+6.0d0*cterm*p6/rik)
                  else
                     p12 = p6 * p6
                     e = expmerge * eps * p12
                     de = -12.0d0 * e / rik
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  if (kskip .eq. -i) then
                     e14 = e14 + e
                     if (i .eq. iv) then
                        de14(1,i) = de14(1,i) + dedx
                        de14(2,i) = de14(2,i) + dedy
                        de14(3,i) = de14(3,i) + dedz
                     else
                        de14(1,i) = de14(1,i) + dedx*redi
                        de14(2,i) = de14(2,i) + dedy*redi
                        de14(3,i) = de14(3,i) + dedz*redi
                        de14(1,iv) = de14(1,iv) + dedx*rediv
                        de14(2,iv) = de14(2,iv) + dedy*rediv
                        de14(3,iv) = de14(3,iv) + dedz*rediv
                     end if
                     if (kmap .eq. kvmap) then
                        de14(1,kmap) = de14(1,kmap) - dedx
                        de14(2,kmap) = de14(2,kmap) - dedy
                        de14(3,kmap) = de14(3,kmap) - dedz
                     else
                        redk = kred(kmap)
                        redkv = 1.0d0 - redk
                        de14(1,kmap) = de14(1,kmap) - dedx*redk
                        de14(2,kmap) = de14(2,kmap) - dedy*redk
                        de14(3,kmap) = de14(3,kmap) - dedz*redk
                        de14(1,kvmap) = de14(1,kvmap) - dedx*redkv
                        de14(2,kvmap) = de14(2,kvmap) - dedy*redkv
                        de14(3,kvmap) = de14(3,kvmap) - dedz*redkv
                     end if
                  else
                     ev = ev + e
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                     end if
                     if (kmap .eq. kvmap) then
                        dev(1,kmap) = dev(1,kmap) - dedx
                        dev(2,kmap) = dev(2,kmap) - dedy
                        dev(3,kmap) = dev(3,kmap) - dedz
                     else
                        redk = kred(kmap)
                        redkv = 1.0d0 - redk
                        dev(1,kmap) = dev(1,kmap) - dedx*redk
                        dev(2,kmap) = dev(2,kmap) - dedy*redk
                        dev(3,kmap) = dev(3,kmap) - dedz*redk
                        dev(1,kvmap) = dev(1,kvmap) - dedx*redkv
                        dev(2,kvmap) = dev(2,kvmap) - dedy*redkv
                        dev(3,kvmap) = dev(3,kvmap) - dedz*redkv
                     end if
                  end if
c
c     increment the total intermolecular energy
c
                  if (kk .le. nvdw) then
                     if (molcule(i) .ne. molcule(kmap)) then
                        einter = einter + e
                     end if
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx + xr*dedx
                     viry = viry + yr*dedy
                     virz = virz + zr*dedz
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge  --  charge-charge potential energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge" calculates the charge-charge interaction energy
c
c
      subroutine echarge
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,in,kn,ic,kc
      real*8 e,f,fi,fik,r,r2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xic,yic,zic,xc,yc,zc
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4,rc5,rc6,rc7
      logical proceed,iuse
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     calculate charge interaction energy term
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         skip(in) = i
         do j = 1, n12(in)
            skip(i12(j,in)) = i * chg12use
         end do
         do j = 1, n13(in)
            skip(i13(j,in)) = i * chg13use
         end do
         do j = 1, n14(in)
            skip(i14(j,in)) = i * chg14use
         end do
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
            if (proceed)  proceed = (skip(kn) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_image)  call image (xc,yc,zc,0)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
                  if (skip(kn) .eq. -i)  fik = fik / chgscale
                  e = fik / r
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            kc = kion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  call image (xc,yc,zc,j)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     xr = xc + xi - x(k) + x(kc)
                     yr = yc + yi - y(k) + y(kc)
                     zr = zc + zi - z(k) + z(kc)
                     r2 = xr*xr + yr*yr + zr*zr
                     r = sqrt(r2)
                     fik = fi * pchg(kk)
                     e = fik / r
c
c     use shifted energy switching if near the cutoff distance
c
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        rc6 = rc3 * rc3
                        rc7 = rc3 * rc4
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                                  + f3*rc3 + f2*rc2 + f1*rc + f0)
                        e = e * taper + trans
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component;
c     interaction of an atom with its own image counts half
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                  end if
               end do
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1  --  charge-charge energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine echarge1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,in,kn,ic,kc
      real*8 e,de,dc,f,fi,fik,r,r2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xic,yic,zic,xc,yc,zc
      real*8 dedx,dedy,dedz
      real*8 dedxc,dedyc,dedzc
      real*8 shift,taper,dtaper,trans,dtrans
      real*8 rc,rc2,rc3,rc4,rc5,rc6,rc7
      logical proceed,iuse
c
      LOGICAL GOPARR,DSKWRK,MASWRK,gopart,mmonly,qmmm
      real*8 einter0,virx0,viry0,virz0 
      integer ME,MASTER,NPROC,IBTYP,IPTIM,mparti
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /tinopt/ mparti,mmonly,qmmm
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     Static parallelisation (D. G. Fedorov)
c     It is important to do point by point parallelisation rather than
c     in consequent blocks of ii below.
c     To avoid double counting, save the pristine values not initialised here.
c
      gopart=goparr.and.iand(mparti,1).eq.0
      einter0=einter
      virx0=virx
      viry0=viry
      virz0=virz
c
c     compute the charge interaction energy and first derivatives
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         skip(in) = i
         do j = 1, n12(in)
            skip(i12(j,in)) = i * chg12use
         end do
         do j = 1, n13(in)
            skip(i13(j,in)) = i * chg13use
         end do
         do j = 1, n14(in)
            skip(i14(j,in)) = i * chg14use
         end do
c        skip is accumulated and should be got on all nodes for safety
c        (the exact skpping logic being unclear to DGF). 
         if(.not.gopart.or.mod(ii,nproc).eq.me) then
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
            if (proceed)  proceed = (skip(kn) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_image)  call image (xc,yc,zc,0)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
                  if (skip(kn) .eq. -i)  fik = fik / chgscale
                  e = fik / r
                  de = -fik / r2
                  dc = 0.0d0
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                           + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                               + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                             + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                     dc = (e * dtaper + dtrans) / rc
                     de = de * taper
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     dc = dc * fgrp
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedxc = dc * xc
                  dedyc = dc * yc
                  dedzc = dc * zc
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,ic) = dec(1,ic) + dedxc
                  dec(2,ic) = dec(2,ic) + dedyc
                  dec(3,ic) = dec(3,ic) + dedzc
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
                  dec(1,kc) = dec(1,kc) - dedxc
                  dec(2,kc) = dec(2,kc) - dedyc
                  dec(3,kc) = dec(3,kc) - dedzc
c                 write(6,*) '  wwwpar',me,ec,e,dec(1,i)
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx + xr*dedx + xc*dedxc
                     viry = viry + yr*dedy + yc*dedyc
                     virz = virz + zr*dedz + zc*dedzc
                  end if
               end if
            end if
         end do
         end if
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica)  then
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion
         if(.not.gopart.or.mod(ii,nproc).eq.me) then
         i = iion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            kc = kion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  call image (xc,yc,zc,j)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     xr = xc + xi - x(k) + x(kc)
                     yr = yc + yi - y(k) + y(kc)
                     zr = zc + zi - z(k) + z(kc)
                     r2 = xr*xr + yr*yr + zr*zr
                     r = sqrt(r2)
                     fik = fi * pchg(kk)
                     e = fik / r
                     de = -fik / r2
                     dc = 0.0d0
c
c     use shifted energy switching if near the cutoff distance
c
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        rc6 = rc3 * rc3
                        rc7 = rc3 * rc4
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                              + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                        trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                                  + f3*rc3 + f2*rc2 + f1*rc + f0)
                        dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                                  + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                                + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                        dc = (e * dtaper + dtrans) / rc
                        de = de * taper
                        e = e * taper + trans
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                        dc = dc * fgrp
                     end if
c
c     form the chain rule terms for derivative expressions
c
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     dedxc = dc * xc
                     dedyc = dc * yc
                     dedzc = dc * zc
c
c     increment the energy and gradient values
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                     dec(1,i) = dec(1,i) + dedx
                     dec(2,i) = dec(2,i) + dedy
                     dec(3,i) = dec(3,i) + dedz
                     dec(1,ic) = dec(1,ic) + dedxc
                     dec(2,ic) = dec(2,ic) + dedyc
                     dec(3,ic) = dec(3,ic) + dedzc
                     if (i .ne. k) then
                        dec(1,k) = dec(1,k) - dedx
                        dec(2,k) = dec(2,k) - dedy
                        dec(3,k) = dec(3,k) - dedz
                        dec(1,kc) = dec(1,kc) - dedxc
                        dec(2,kc) = dec(2,kc) - dedyc
                        dec(3,kc) = dec(3,kc) - dedzc
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*dedx + xc*dedxc
                        viry = viry + yr*dedy + yc*dedyc
                        virz = virz + zr*dedz + zc*dedzc
                     end if
                  end if
               end do
            end if
         end do
         endif
      end do
      endif
      if(gopart) then
         call ddi_gsumf(4000,dec,3*n)
         call ddi_gsumf(4001,ec,1)
         einter=einter-einter0
         call ddi_gsumf(4002,einter,1)
         einter=einter+einter0
         if (isobaric) then
            virx=virx-virx0
            viry=viry-viry0
            virz=virz-virz0
            call ddi_gsumf(4003,virx,1)
            call ddi_gsumf(4004,viry,1)
            call ddi_gsumf(4005,virz,1)
            virx=virx+virx0
            viry=viry+viry0
            virz=virz+virz0
         end if
      end if
c
c     write(6,*) 'wwwpar',ec,einter,((dec(i,j),i=1,1),j=1,1)
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
c     ##  subroutine echarge2  --  atom-wise charge-charge Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge2" calculates second derivatives of the
c     charge-charge interaction energy for a single atom
c
c
      subroutine echarge2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,kk,in,kn
      integer jcell,skip(maxatm)
      real*8 e,de,fi,fik,fgrp
      real*8 d2e,d2edx,d2edy,d2edz
      real*8 xi,yi,zi,xr,yr,zr,term(3,3)
      real*8 shift,taper,dtaper,d2taper
      real*8 trans,dtrans,d2trans
      real*8 r,r2,r3,r4,r5,r6,r7
      logical proceed
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            in = jion(k)
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
c
c     set the list of atoms to be skipped
c
      do k = 1, n
         skip(k) = 0
      end do
      skip(in) = i
      do j = 1, n12(in)
         skip(i12(j,in)) = i * chg12use
      end do
      do j = 1, n13(in)
         skip(i13(j,in)) = i * chg13use
      end do
      do j = 1, n14(in)
         skip(i14(j,in)) = i * chg14use
      end do
c
c     set cutoff distances and switching function coefficients
c
      call switch ('CHARGE')
c
c     calculate the charge interaction energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
         if (proceed)  proceed = (skip(kn) .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            if (use_image)  call image (xr,yr,zr,0)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               fik = fi * pchg(kk)
               if (skip(kn) .eq. -i)  fik = fik / chgscale
c
c     compute chain rule terms for Hessian matrix elements
c
               de = -fik / r2
               d2e = -2.0d0 * de/r
c
c     use shifted energy switching if near the cutoff distance
c
               if (r2 .gt. cut2) then
                  e = fik / r
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  r3 = r2 * r
                  r4 = r2 * r2
                  r5 = r2 * r3
                  r6 = r3 * r3
                  r7 = r3 * r4
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                        + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                  d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                         + 6.0d0*c3*r + 2.0d0*c2
                  trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                            + f3*r3 + f2*r2 + f1*r + f0)
                  dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                            + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                            + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                  d2trans = fik * (42.0d0*f7*r5 + 30.0d0*f6*r4
     &                             + 20.0d0*f5*r3 + 12.0d0*f4*r2
     &                             + 6.0d0*f3*r + 2.0d0*f2)
                  d2e = e*d2taper + 2.0d0*de*dtaper
     &                     + d2e*taper + d2trans
                  de = e*dtaper + de*taper + dtrans
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  de = de * fgrp
                  d2e = d2e * fgrp
               end if
c
c     form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do kk = 1, nion
         k = iion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            do jcell = 1, ncell
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr,jcell)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
c
c     compute chain rule terms for Hessian matrix elements
c
                  de = -fik / r2
                  d2e = -2.0d0 * de/r
c
c     use shifted energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik / r
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                            + 6.0d0*c3*r + 2.0d0*c2
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                               + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                               + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                     d2trans = fik * (42.0d0*f7*r5 + 30.0d0*f6*r4
     &                                + 20.0d0*f5*r3 + 12.0d0*f4*r2
     &                                + 6.0d0*f3*r + 2.0d0*f2)
                     d2e = e*d2taper + 2.0d0*de*dtaper
     &                        + d2e*taper + d2trans
                     de = e*dtaper + de*taper + dtrans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     de = de * fgrp
                     d2e = d2e * fgrp
                  end if
c
c     form the individual Hessian element components
c
                  de = de / r
                  d2e = (d2e-de) / r2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  do j = 1, 3
                     hessx(j,i) = hessx(j,i) + term(1,j)
                     hessy(j,i) = hessy(j,i) + term(2,j)
                     hessz(j,i) = hessz(j,i) + term(3,j)
                     hessx(j,k) = hessx(j,k) - term(1,j)
                     hessy(j,k) = hessy(j,k) - term(2,j)
                     hessz(j,k) = hessz(j,k) - term(3,j)
                  end do
               end if
            end do
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3  --  charge-charge energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3" calculates the charge-charge interaction energy;
c     also partitions the energy among the atoms
c
c
      subroutine echarge3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moment.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,in,kn,ic,kc
      real*8 e,f,fi,fik,r,r2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xic,yic,zic,xc,yc,zc
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4,rc5,rc6,rc7
      real*8 weight,xcenter,ycenter,zcenter
      real*8 totchg,xsum,ysum,zsum
      logical header,huge,proceed,iuse
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      header = .true.
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     find the center of mass of the total system
c
      weight = 0.0d0
      xcenter = 0.0d0
      ycenter = 0.0d0
      zcenter = 0.0d0
      do i = 1, n
         weight = weight + mass(i)
         xcenter = xcenter + x(i)*mass(i)
         ycenter = ycenter + y(i)*mass(i)
         zcenter = zcenter + z(i)*mass(i)
      end do
      xcenter = xcenter / weight
      ycenter = ycenter / weight
      zcenter = zcenter / weight
c
c     get net charge and dipole components relative to center of mass
c
      totchg = 0.0d0
      xsum = 0.0d0
      ysum = 0.0d0
      zsum = 0.0d0
      do ii = 1, nion
         i = iion(ii)
         totchg = totchg + pchg(ii)
         xsum = xsum + (x(i)-xcenter)*pchg(ii)
         ysum = ysum + (y(i)-ycenter)*pchg(ii)
         zsum = zsum + (z(i)-zcenter)*pchg(ii)
      end do
      netchg = netchg + totchg
      xdipole = xdipole + debye*xsum
      ydipole = ydipole + debye*ysum
      zdipole = zdipole + debye*zsum
c
c     compute and partition the charge interaction energy
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         skip(in) = i
         do j = 1, n12(in)
            skip(i12(j,in)) = i * chg12use
         end do
         do j = 1, n13(in)
            skip(i13(j,in)) = i * chg13use
         end do
         do j = 1, n14(in)
            skip(i14(j,in)) = i * chg14use
         end do
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
            if (proceed)  proceed = (skip(kn) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_image)  call image (xc,yc,zc,0)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
                  if (skip(kn) .eq. -i)  fik = fik / chgscale
                  e = fik / r
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  nec = nec + 1
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
c
c     print a warning if the energy of this charge pair is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Charge-Charge',
     &                             ' Interactions :',
     &                          //,' Type',11x,'Atom Names',
     &                             16x,'Charges',5x,'Distance',
     &                             5x,'Energy',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),
     &                                pchg(ii),pchg(kk),r,e
   20                format (' Charge   ',i5,'-',a3,1x,i5,'-',
     &                         a3,8x,2f7.2,f10.4,f12.4)
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         do kk = ii, nion
            k = iion(kk)
            kc = kion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  call image (xc,yc,zc,j)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     xr = xc + xi - x(k) + x(kc)
                     yr = yc + yi - y(k) + y(kc)
                     zr = zc + zi - z(k) + z(kc)
                     r2 = xr*xr + yr*yr + zr*zr
                     r = sqrt(r2)
                     fik = fi * pchg(kk)
                     e = fik / r
c
c     use shifted energy switching if near the cutoff distance
c
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        rc6 = rc3 * rc3
                        rc7 = rc3 * rc4
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                                  + f3*rc3 + f2*rc2 + f1*rc + f0)
                        e = e * taper + trans
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                     nec = nec + 1
                     if (i .eq. k) then
                        ec = ec + 0.5d0*e
                        aec(i) = aec(i) + 0.5d0*e
                     else
                        ec = ec + e
                        aec(i) = aec(i) + 0.5d0*e
                        aec(k) = aec(k) + 0.5d0*e
                     end if
c
c     print a warning if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 100.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Charge-Charge',
     &                                ' Interactions :',
     &                             //,' Type',11x,'Atom Names',
     &                                16x,'Charges',5x,'Distance',
     &                                5x,'Energy',/)
                        end if
                        write (iout,40)  i,name(i),k,name(k),
     &                                   pchg(ii),pchg(kk),r,e
   40                   format (' Charge   ',i5,'-',a3,1x,i5,'-',
     &                            a3,' (XTAL) ',2f7.2,f10.4,f12.4)
                     end if
                  end if
               end do
            end if
         end do
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge4  --  charge-charge potential energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge4" calculates the charge-charge interaction energy
c     using the method of lights to locate neighboring atoms
c
c
      subroutine echarge4
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'iounit.i'
      include 'light.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,ii,kk,in,ic
      integer kgy,kgz,start,stop
      integer kskip,skip(maxatm)
      integer kmap,kcmap,map(maxlight)
      real*8 e,f,fi,fik,rik,rik2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xic,yic,zic,xc,yc,zc
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4,rc5,rc6,rc7
      real*8 xsort(maxlight),ysort(maxlight),zsort(maxlight)
      logical proceed,iuse,repeat
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nion
         k = kion(i)
         xsort(i) = x(k)
         ysort(i) = y(k)
         zsort(i) = z(k)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nion,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         skip(in) = i
         do j = 1, n12(in)
            skip(i12(j,in)) = i * chg12use
         end do
         do j = 1, n13(in)
            skip(i13(j,in)) = i * chg13use
         end do
         do j = 1, n14(in)
            skip(i14(j,in)) = i * chg14use
         end do
         xic = xsort(rgx(ii))
         yic = ysort(rgy(ii))
         zic = zsort(rgz(ii))
         xi = x(i) - x(ic)
         yi = y(i) - y(ic)
         zi = z(i) - z(ic)
         fi = f * pchg(ii)
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            if (kk .le. nion) then
               kskip = skip(jion(kk))
               if (kskip .eq. i)  goto 20
            else
               kskip = 0
            end if
            kmap = iion(map(kk))
            kcmap = kion(map(kk))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,kmap,0,0,0)
            if (proceed)  proceed = (iuse .or. use(kmap)
     &                                .or. use(kcmap))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - xsort(j)
               yc = yic - ysort(kgy)
               zc = zic - zsort(kgz)
               if (use_image) then
                  if (abs(xc) .gt. xcell2)  xc = xc - sign(xcell,xc)
                  if (abs(yc) .gt. ycell2)  yc = yc - sign(ycell,yc)
                  if (abs(zc) .gt. zcell2)  zc = zc - sign(zcell,zc)
                  if (monoclinic) then
                     xc = xc + zc*beta_cos
                     zc = zc * beta_sin
                  else if (triclinic) then
                     xc = xc + yc*gamma_cos + zc*beta_cos
                     yc = yc*gamma_sin + zc*beta_term
                     zc = zc * gamma_term
                  end if
               end if
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(kmap) + x(kcmap)
                  yr = yc + yi - y(kmap) + y(kcmap)
                  zr = zc + zi - z(kmap) + z(kcmap)
                  rik2 = xr*xr + yr*yr + zr*zr
                  rik = sqrt(rik2)
                  fik = fi * pchg(map(kk))
                  if (kskip .eq. -i)  fik = fik / chgscale
                  e = fik / rik
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge5  --  charge-charge energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge5" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using the method of lights to locate neighboring atoms
c
c
      subroutine echarge5
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'iounit.i'
      include 'light.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,ii,kk,in,ic
      integer kgy,kgz,start,stop
      integer kskip,skip(maxatm)
      integer kmap,kcmap,map(maxlight)
      real*8 e,de,dc,f,fi,fik,rik,rik2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xic,yic,zic,xc,yc,zc
      real*8 dedx,dedy,dedz,dedxc,dedyc,dedzc
      real*8 shift,taper,dtaper,trans,dtrans
      real*8 rc,rc2,rc3,rc4,rc5,rc6,rc7
      real*8 xsort(maxlight),ysort(maxlight),zsort(maxlight)
      logical proceed,iuse,repeat
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nion
         k = kion(i)
         xsort(i) = x(k)
         ysort(i) = y(k)
         zsort(i) = z(k)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (nion,map,xsort,ysort,zsort)
c
c     now, loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         iuse = (use(i) .or. use(ic))
         skip(in) = i
         do j = 1, n12(in)
            skip(i12(j,in)) = i * chg12use
         end do
         do j = 1, n13(in)
            skip(i13(j,in)) = i * chg13use
         end do
         do j = 1, n14(in)
            skip(i14(j,in)) = i * chg14use
         end do
         xic = xsort(rgx(ii))
         yic = ysort(rgy(ii))
         zic = zsort(rgz(ii))
         xi = x(i) - x(ic)
         yi = y(i) - y(ic)
         zi = z(i) - z(ic)
         fi = f * pchg(ii)
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            if (kk .le. nion) then
               kskip = skip(jion(kk))
               if (kskip .eq. i)  goto 20
            else
               kskip = 0
            end if
            kmap = iion(map(kk))
            kcmap = kion(map(kk))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,kmap,0,0,0)
            if (proceed)  proceed = (iuse .or. use(kmap)
     &                                .or. use(kcmap))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - xsort(j)
               yc = yic - ysort(kgy)
               zc = zic - zsort(kgz)
               if (use_image) then
                  if (abs(xc) .gt. xcell2)  xc = xc - sign(xcell,xc)
                  if (abs(yc) .gt. ycell2)  yc = yc - sign(ycell,yc)
                  if (abs(zc) .gt. zcell2)  zc = zc - sign(zcell,zc)
                  if (monoclinic) then
                     xc = xc + zc*beta_cos
                     zc = zc * beta_sin
                  else if (triclinic) then
                     xc = xc + yc*gamma_cos + zc*beta_cos
                     yc = yc*gamma_sin + zc*beta_term
                     zc = zc * gamma_term
                  end if
               end if
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(kmap) + x(kcmap)
                  yr = yc + yi - y(kmap) + y(kcmap)
                  zr = zc + zi - z(kmap) + z(kcmap)
                  rik2 = xr*xr + yr*yr + zr*zr
                  rik = sqrt(rik2)
                  fik = fi * pchg(map(kk))
                  if (kskip .eq. -i)  fik = fik / chgscale
                  e = fik / rik
                  de = -fik / rik2
                  dc = 0.0d0
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                           + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                               + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                             + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                     dc = (e * dtaper + dtrans) / rc
                     de = de * taper
                     e = e * taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     dc = dc * fgrp
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedxc = dc * xc
                  dedyc = dc * yc
                  dedzc = dc * zc
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,ic) = dec(1,ic) + dedxc
                  dec(2,ic) = dec(2,ic) + dedyc
                  dec(3,ic) = dec(3,ic) + dedzc
                  dec(1,kmap) = dec(1,kmap) - dedx
                  dec(2,kmap) = dec(2,kmap) - dedy
                  dec(3,kmap) = dec(3,kmap) - dedz
                  dec(1,kcmap) = dec(1,kcmap) - dedxc
                  dec(2,kcmap) = dec(2,kcmap) - dedyc
                  dec(3,kcmap) = dec(3,kcmap) - dedzc
c
c     increment the total intermolecular energy
c
                  if (kk .le. nion) then
                     if (molcule(i) .ne. molcule(kmap)) then
                        einter = einter + e
                     end if
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx + xr*dedx + xc*dedxc
                     viry = viry + yr*dedy + yc*dedyc
                     virz = virz + zr*dedz + zc*dedzc
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge6  --  smoothed charge-charge energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge6" calculates the smoothed charge-charge
c     interaction energy
c
c
      subroutine echarge6
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'warp.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,in,kn
      real*8 e,f,fi,fik,r,r2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 erf,defterm,time
      logical proceed,iuse
c
c
c     zero out the charge interaction energy
c
      ec = 0.0d0
      if (nion .eq. 0)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     scale deformation time by the diffusion coefficient
c
      if (deform .ne. 0.0d0) then
         time = diffc * deform
         defterm = 0.5d0 / sqrt(time)
      end if
c
c     calculate charge interaction energy term
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         iuse = (use(i))
         skip(in) = i
         do j = 1, n12(in)
            skip(i12(j,in)) = i * chg12use
         end do
         do j = 1, n13(in)
            skip(i13(j,in)) = i * chg13use
         end do
         do j = 1, n14(in)
            skip(i14(j,in)) = i * chg14use
         end do
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k))
            if (proceed)  proceed = (skip(kn) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
                  if (skip(kn) .eq. -i)  fik = fik / chgscale
                  e = fik / r
c
c     transform the potential via diffusional smoothing
c
                  if (deform .ne. 0.0d0)  e = e * erf(defterm*r)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-charge energy component
c
                  ec = ec + e
               end if
            end if
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge7  --  smoothed charge energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge7" calculates the smoothed charge-charge interaction
c     energy and first derivatives with respect to Cartesian coordinates
c    
c
c
      subroutine echarge7
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'math.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'warp.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,in,kn
      real*8 e,de,f,fi,fik,r,r2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 dedx,dedy,dedz
      real*8 erf,erfterm
      real*8 expcut,expterm
      real*8 defterm,time,time2
      logical proceed,iuse
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     scale deformation time by the diffusion coefficient
c
      if (deform .ne. 0.0d0) then
         expcut = -50.0d0
         time = diffc * deform
         time2 = sqrt(time)
         defterm = 0.5d0 / time2
      end if
c
c     compute charge interaction energy and first derivatives
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         iuse = (use(i))
         skip(in) = i
         do j = 1, n12(in)
            skip(i12(j,in)) = i * chg12use
         end do
         do j = 1, n13(in)
            skip(i13(j,in)) = i * chg13use
         end do
         do j = 1, n14(in)
            skip(i14(j,in)) = i * chg14use
         end do
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k))
            if (proceed)  proceed = (skip(kn) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
                  if (skip(k) .eq. -i)  fik = fik / chgscale
                  e = fik / r
                  de = -fik / r2
c
c     transform the potential via diffusional smoothing
c
                  if (deform .ne. 0.0d0) then
                     erfterm = erf(defterm*r)
                     expterm = -0.25d0 * r2 / time
                     if (expterm .gt. expcut) then
                        expterm = exp(expterm) / (sqrtpi*time2)
                     else
                        expterm = 0.0d0
                     end if
                     e = e * erfterm
                     de = de*erfterm + fik*expterm/r
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
               end if
            end if
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echarge8  --  atom-wise smoothed charge Hessian  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge8" calculates second derivatives of the smoothed
c     charge-charge interaction energy for a single atom
c
c
      subroutine echarge8 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'math.i'
      include 'shunt.i'
      include 'units.i'
      include 'warp.i'
      integer i,j,k,skip(maxatm)
      integer kk,in,kn
      real*8 fi,fik,de,d2e,fgrp
      real*8 d2edx,d2edy,d2edz
      real*8 r,r2,term(3,3)
      real*8 xi,yi,zi,xr,yr,zr
      real*8 erf,erfterm
      real*8 expcut,expterm
      real*8 defterm,time,time2
      logical proceed
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            in = jion(k)
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
c
c     set the list of atoms to be skipped
c
      do k = 1, n
         skip(k) = 0
      end do
      skip(in) = i
      do j = 1, n12(in)
         skip(i12(j,in)) = i * chg12use
      end do
      do j = 1, n13(in)
         skip(i13(j,in)) = i * chg13use
      end do
      do j = 1, n14(in)
         skip(i14(j,in)) = i * chg14use
      end do
c
c     set the coefficients of the switching function
c
      call switch ('CHARGE')
c
c     scale deformation time by the diffusion coefficient
c
      if (deform .ne. 0.0d0) then
         expcut = -50.0d0
         time = diffc * deform
         time2 = sqrt(time)
         defterm = 0.5d0 / time2
      end if
c
c     calculate the charge interaction energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
         if (proceed)  proceed = (skip(kn) .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               fik = fi * pchg(kk)
               if (skip(kn) .eq. -i)  fik = fik / chgscale
c
c     compute chain rule terms for Hessian matrix elements
c
               de = -fik / r2
               d2e = -2.0d0 * de/r
c
c     transform the potential via diffusional smoothing
c
               if (deform .ne. 0.0d0) then
                  erfterm = erf(defterm*r)
                  expterm = -0.25d0 * r2 / time
                  if (expterm .gt. expcut) then
                     expterm = exp(expterm) / (sqrtpi*time2)
                  else
                     expterm = 0.0d0
                  end if
                  de = de*erfterm + fik*expterm/r
                  d2e = -2.0d0*de/r - 0.5d0*fik*expterm/time
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  de = de * fgrp
                  d2e = d2e * fgrp
               end if
c
c     form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echarge9  --  smoothed charge energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge9" calculates the smoothed charge-charge interaction
c     energy; also partitions the energy among the atoms
c
c
      subroutine echarge9
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'warp.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,in,kn
      real*8 e,f,fi,fik,r,r2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 erf,defterm,time
      logical header,huge,proceed,iuse
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      header = .true.
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('CHARGE')
c
c     scale deformation time by the diffusion coefficient
c
      if (deform .ne. 0.0d0) then
         time = diffc * deform
         defterm = 0.5d0 / sqrt(time)
      end if
c
c     calculate charge interaction energy term
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         iuse = (use(i))
         skip(in) = i
         do j = 1, n12(in)
            skip(i12(j,in)) = i * chg12use
         end do
         do j = 1, n13(in)
            skip(i13(j,in)) = i * chg13use
         end do
         do j = 1, n14(in)
            skip(i14(j,in)) = i * chg14use
         end do
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. use(k))
            if (proceed)  proceed = (skip(kn) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  fik = fi * pchg(kk)
                  if (skip(k) .eq. -i)  fik = fik / chgscale
                  e = fik / r
c
c     transform the potential via diffusional smoothing
c
                  if (deform .ne. 0.0d0)  e = e * erf(defterm*r)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the charge-charge energy and partition between atoms
c
                  nec = nec + 1
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
c
c     print a warning if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Charge-Charge',
     &                             ' Interactions :',
     &                          //,' Type',11x,'Atom Names',
     &                             16x,'Charges',5x,'Distance',
     &                             5x,'Energy',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),pchg(ii),
     &                                pchg(kk),r,e
   20                format (' Charge   ',i5,'-',a3,1x,i5,'-',
     &                         a3,8x,2f7.2,f10.4,f12.4)
                  end if
               end if
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echgdpl  --  charge-dipole potential energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echgdpl" calculates the charge-dipole interaction energy
c
c
      subroutine echgdpl
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'dipole.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,i1,k1,k2
      integer skip(maxatm)
      real*8 e,rk2,rkr3,dotk
      real*8 f,fi,fik,fgrp
      real*8 r,r2,r3,r4,r5,taper
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xr,yr,zr
      logical proceed
c
c
c     zero out the overall charge-dipole interaction energy
c
      ecd = 0.0d0
      if (ndipole.eq.0 .or. nion.eq.0)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye * dielec)
      call switch ('CHGDPL')
c
c     get the total energy by looping over each charge-dipole pair
c
      do i = 1, nion
         i1 = iion(i)
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
            if (proceed)  proceed = (skip(k1).ne.i1 .and.
     &                                 skip(k2).ne.i1)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = xi - x(k1) - xk*sdpl(k)
               yr = yi - y(k1) - yk*sdpl(k)
               zr = zi - z(k1) - zk*sdpl(k)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  fik = fi * bdpl(k)
                  rk2 = xk*xk + yk*yk + zk*zk
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
                  e = fik * dotk / rkr3
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-dipole energy component
c
                  ecd = ecd + e
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, nion
         i1 = iion(i)
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  xr = xi - x(k1) - xk*sdpl(k)
                  yr = yi - y(k1) - yk*sdpl(k)
                  zr = zi - z(k1) - zk*sdpl(k)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     fik = fi * bdpl(k)
                     rk2 = xk*xk + yk*yk + zk*zk
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
                     e = fik * dotk / rkr3
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall charge-dipole energy component
c
                     ecd = ecd + e
                  end if
               end do
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echgdpl1  --  charge-dipole energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echgdpl1" calculates the charge-dipole interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine echgdpl1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'dipole.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,i1,k1,k2
      integer skip(maxatm)
      real*8 e,rk2,rkr3,dotk
      real*8 f,fi,fik,fgrp,sk1,sk2
      real*8 term,term2,term3
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xr,yr,zr
      real*8 termx,termy,termz
      real*8 termxk,termyk,termzk
      real*8 dedxi1,dedyi1,dedzi1
      real*8 dedxk1,dedyk1,dedzk1
      real*8 dedxk2,dedyk2,dedzk2
      real*8 r,r2,r3,r4,r5,taper,dtaper
      real*8 dtaperx,dtapery,dtaperz
      logical proceed
c
c
c     zero out the overall charge-dipole interaction energy
c     and set up the constants for the calculation
c
      ecd = 0.0d0
      do i = 1, n
         decd(1,i) = 0.0d0
         decd(2,i) = 0.0d0
         decd(3,i) = 0.0d0
      end do
      if (ndipole.eq.0 .or. nion.eq.0)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye * dielec)
      call switch ('CHGDPL')
c
c     get energy and derivs by looping over each charge-dipole pair
c
      do i = 1, nion
         i1 = iion(i)
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
            if (proceed)  proceed = (skip(k1).ne.i1 .and.
     &                                 skip(k2).ne.i1)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = xi - x(k1) - xk*sk2
               yr = yi - y(k1) - yk*sk2
               zr = zi - z(k1) - zk*sk2
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  fik = fi * bdpl(k)
                  rk2 = xk*xk + yk*yk + zk*zk
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
c
c     form the energy and master chain rule term for derivatives
c
                  e = fik * dotk / rkr3
                  term = fik / rkr3
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     term = term * fgrp
                  end if
c
c     secondary chain rule terms for derivative expressions
c
                  term2 = -3.0d0 * dotk / r2
                  term3 = dotk / rk2
                  termx = term * (xk+xr*term2)
                  termy = term * (yk+yr*term2)
                  termz = term * (zk+zr*term2)
                  termxk = term * (xr-xk*term3)
                  termyk = term * (yr-yk*term3)
                  termzk = term * (zr-zk*term3)
                  dedxi1 = termx
                  dedyi1 = termy
                  dedzi1 = termz
                  dedxk1 = -sk1*termx - termxk
                  dedyk1 = -sk1*termy - termyk
                  dedzk1 = -sk1*termz - termzk
                  dedxk2 = -sk2*termx + termxk
                  dedyk2 = -sk2*termy + termyk
                  dedzk2 = -sk2*termz + termzk
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     dtaper = dtaper * e/r
                     dtaperx = xr * dtaper
                     dtapery = yr * dtaper
                     dtaperz = zr * dtaper
                     e = e * taper
                     dedxi1 = dedxi1*taper + dtaperx
                     dedyi1 = dedyi1*taper + dtapery
                     dedzi1 = dedzi1*taper + dtaperz
                     dedxk1 = dedxk1*taper - sk1*dtaperx
                     dedyk1 = dedyk1*taper - sk1*dtapery
                     dedzk1 = dedzk1*taper - sk1*dtaperz
                     dedxk2 = dedxk2*taper - sk2*dtaperx
                     dedyk2 = dedyk2*taper - sk2*dtapery
                     dedzk2 = dedzk2*taper - sk2*dtaperz
                  end if
c
c     increment the overall energy and derivative expressions
c
                  ecd = ecd + e
                  decd(1,i1) = decd(1,i1) + dedxi1
                  decd(2,i1) = decd(2,i1) + dedyi1
                  decd(3,i1) = decd(3,i1) + dedzi1
                  decd(1,k1) = decd(1,k1) + dedxk1
                  decd(2,k1) = decd(2,k1) + dedyk1
                  decd(3,k1) = decd(3,k1) + dedzk1
                  decd(1,k2) = decd(1,k2) + dedxk2
                  decd(2,k2) = decd(2,k2) + dedyk2
                  decd(3,k2) = decd(3,k2) + dedzk2
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx - (xr+xk*sk2)*dedxk1
     &                              - (xr-xk*sk1)*dedxk2
                     viry = viry - (yr+yk*sk2)*dedyk1
     &                              - (yr-yk*sk1)*dedyk2
                     virz = virz - (zr+zk*sk2)*dedzk1
     &                              - (zr-zk*sk1)*dedzk2
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, nion
         i1 = iion(i)
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               do j = 1, ncell
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  xr = xi - x(k1) - xk*sk2
                  yr = yi - y(k1) - yk*sk2
                  zr = zi - z(k1) - zk*sk2
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     fik = fi * bdpl(k)
                     rk2 = xk*xk + yk*yk + zk*zk
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
c
c     form the energy and master chain rule term for derivatives
c
                     e = fik * dotk / rkr3
                     term = fik / rkr3
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        term = term * fgrp
                     end if
c
c     secondary chain rule terms for derivative expressions
c
                     term2 = -3.0d0 * dotk / r2
                     term3 = dotk / rk2
                     termx = term * (xk+xr*term2)
                     termy = term * (yk+yr*term2)
                     termz = term * (zk+zr*term2)
                     termxk = term * (xr-xk*term3)
                     termyk = term * (yr-yk*term3)
                     termzk = term * (zr-zk*term3)
                     dedxi1 = termx
                     dedyi1 = termy
                     dedzi1 = termz
                     dedxk1 = -sk1*termx - termxk
                     dedyk1 = -sk1*termy - termyk
                     dedzk1 = -sk1*termz - termzk
                     dedxk2 = -sk2*termx + termxk
                     dedyk2 = -sk2*termy + termyk
                     dedzk2 = -sk2*termz + termzk
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        dtaper = dtaper * e/r
                        dtaperx = xr * dtaper
                        dtapery = yr * dtaper
                        dtaperz = zr * dtaper
                        e = e * taper
                        dedxi1 = dedxi1*taper + dtaperx
                        dedyi1 = dedyi1*taper + dtapery
                        dedzi1 = dedzi1*taper + dtaperz
                        dedxk1 = dedxk1*taper - sk1*dtaperx
                        dedyk1 = dedyk1*taper - sk1*dtapery
                        dedzk1 = dedzk1*taper - sk1*dtaperz
                        dedxk2 = dedxk2*taper - sk2*dtaperx
                        dedyk2 = dedyk2*taper - sk2*dtapery
                        dedzk2 = dedzk2*taper - sk2*dtaperz
                     end if
c
c     increment the overall energy and derivative expressions
c
                     ecd = ecd + e
                     decd(1,i1) = decd(1,i1) + dedxi1
                     decd(2,i1) = decd(2,i1) + dedyi1
                     decd(3,i1) = decd(3,i1) + dedzi1
                     decd(1,k1) = decd(1,k1) + dedxk1
                     decd(2,k1) = decd(2,k1) + dedyk1
                     decd(3,k1) = decd(3,k1) + dedzk1
                     decd(1,k2) = decd(1,k2) + dedxk2
                     decd(2,k2) = decd(2,k2) + dedyk2
                     decd(3,k2) = decd(3,k2) + dedzk2
                  end if
c
c     increment the total intermolecular energy
c
                  einter = einter + e
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx - (xr+xk*sk2)*dedxk1
     &                              - (xr-xk*sk1)*dedxk2
                     viry = viry - (yr+yk*sk2)*dedyk1
     &                              - (yr-yk*sk1)*dedyk2
                     virz = virz - (zr+zk*sk2)*dedzk1
     &                              - (zr-zk*sk1)*dedzk2
                  end if
               end do
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echgdpl2  --  atom-wise charge-dipole Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echgdpl2" calculates second derivatives of the
c     charge-dipole interaction energy for a single atom
c
c
      subroutine echgdpl2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'dipole.i'
      include 'group.i'
      include 'hessn.i'
      include 'shunt.i'
      include 'units.i'
      integer i,k,ii,i1,k1,k2,jcell
      integer skip(maxatm),omit(maxatm)
      real*8 f,fi,fk,fik,fgrp,sk1,sk2
      real*8 xi,yi,zi,xk,yk,zk,xr,yr,zr,xq,yq,zq
      real*8 e,r2,rk2,rkr3,dotk,term,term2,part,part2
      real*8 termx,termy,termz,termxk,termyk,termzk
      real*8 xrr2,yrr2,zrr2,xkrk2,ykrk2,zkrk2
      real*8 dotk2,dotkr2,dotkrk2,factor,factork
      real*8 dedxi1,dedyi1,dedzi1,dedxk1,dedyk1
      real*8 dedzk1,dedxk2,dedyk2,dedzk2
      real*8 dtdxi1,dtdyi1,dtdzi1,dtdxk1,dtdyk1
      real*8 dtdzk1,dtdxk2,dtdyk2,dtdzk2
      real*8 dtxdxi1,dtxkdxi1,dtxdxk1,dtxkdxk1,dtxdxk2,dtxkdxk2
      real*8 dtydxi1,dtykdxi1,dtydxk1,dtykdxk1,dtydxk2,dtykdxk2
      real*8 dtzdxi1,dtzkdxi1,dtzdxk1,dtzkdxk1,dtzdxk2,dtzkdxk2
      real*8 dtxdyi1,dtxkdyi1,dtxdyk1,dtxkdyk1,dtxdyk2,dtxkdyk2
      real*8 dtydyi1,dtykdyi1,dtydyk1,dtykdyk1,dtydyk2,dtykdyk2
      real*8 dtzdyi1,dtzkdyi1,dtzdyk1,dtzkdyk1,dtzdyk2,dtzkdyk2
      real*8 dtxdzi1,dtxkdzi1,dtxdzk1,dtxkdzk1,dtxdzk2,dtxkdzk2
      real*8 dtydzi1,dtykdzi1,dtydzk1,dtykdzk1,dtydzk2,dtykdzk2
      real*8 dtzdzi1,dtzkdzi1,dtzdzk1,dtzkdzk1,dtzdzk2,dtzkdzk2
      real*8 r,r3,r4,r5,taper,dtaper,d2taper
      real*8 dtaperx,dtapery,dtaperz,d2taperxx,d2taperyy
      real*8 d2taperzz,d2taperxy,d2taperxz,d2taperyz
      logical proceed
c
c
c     zero out the lists of atoms to be skipped
c
      if (ndipole.eq.0 .or. nion.eq.0)  return
      do k = 1, n
         skip(k) = 0
         omit(k) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye * dielec)
      call switch ('CHGDPL')
c
c     first see if the atom of interest carries a charge
c
      do ii = 1, nion
         i1 = iion(ii)
         if (i1 .ne. i)  goto 10
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(ii)
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
            if (proceed)  proceed = (skip(k1).ne.i1 .and.
     &                                 skip(k2).ne.i1)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               xk = x(k1) - x(k2)
               yk = y(k1) - y(k2)
               zk = z(k1) - z(k2)
               xr = xi - x(k1) + xk*sk2
               yr = yi - y(k1) + yk*sk2
               zr = zi - z(k1) + zk*sk2
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = -fi * bdpl(k)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  fik = fik * fgrp
c
c     some abbreviations used in various chain rule terms
c
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  dotk2 = 2.0d0 * dotk
                  dotkr2 = dotk / r2
c
c     now, form the chain rule terms for first derivatives
c
                  term = fik / rkr3
                  term2 = -3.0d0 * dotk
                  termx = term * (xk+xrr2*term2)
                  termy = term * (yk+yrr2*term2)
                  termz = term * (zk+zrr2*term2)
                  termxk = term * (xr-dotk*xkrk2)
                  termyk = term * (yr-dotk*ykrk2)
                  termzk = term * (zr-dotk*zkrk2)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik * dotk / rkr3
                     dedxi1 = termx
                     dedyi1 = termy
                     dedzi1 = termz
                     dedxk1 = -sk1*termx + termxk
                     dedyk1 = -sk1*termy + termyk
                     dedzk1 = -sk1*termz + termzk
                     dedxk2 = -sk2*termx - termxk
                     dedyk2 = -sk2*termy - termyk
                     dedzk2 = -sk2*termz - termzk
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                            + 6.0d0*c3*r + 2.0d0*c2
                     dtaper = dtaper / r
                     dtaperx = xr * dtaper
                     dtapery = yr * dtaper
                     dtaperz = zr * dtaper
                     d2taper = e * (d2taper-dtaper)
                     dtaper = e * dtaper
                     d2taperxx = xr*xrr2*d2taper + dtaper
                     d2taperxy = xr*yrr2*d2taper
                     d2taperxz = xr*zrr2*d2taper
                     d2taperyy = yr*yrr2*d2taper + dtaper
                     d2taperyz = yr*zrr2*d2taper
                     d2taperzz = zr*zrr2*d2taper + dtaper
                     term = term * taper
                     termx = termx * taper
                     termy = termy * taper
                     termz = termz * taper
                     termxk = termxk * taper
                     termyk = termyk * taper
                     termzk = termzk * taper
                  end if
c
c     next, find the second derivative chain rule terms
c
                  dtdxi1 = -3.0d0 * xrr2
                  part = xk - dotk2*xrr2
                  factor = -3.0d0 * (dotkr2 + xrr2*part)
                  factork = 1.0d0 - xk*xkrk2
                  dtxdxi1 = dtdxi1*termx + term*factor
                  dtxkdxi1 = dtdxi1*termxk + term*factork
                  factor = -3.0d0 * yrr2 * part
                  factork = -yk * xkrk2
                  dtydxi1 = dtdxi1*termy + term*factor
                  dtykdxi1 = dtdxi1*termyk + term*factork
                  factor = -3.0d0 * zrr2 * part
                  factork = -zk * xkrk2
                  dtzdxi1 = dtdxi1*termz + term*factor
                  dtzkdxi1 = dtdxi1*termzk + term*factork
c
                  dtdyi1 = -3.0d0 * yrr2
                  part = yk - dotk2*yrr2
                  factor = -3.0d0 * xrr2 * part
                  factork = -xk * ykrk2
                  dtxdyi1 = dtdyi1*termx + term*factor
                  dtxkdyi1 = dtdyi1*termxk + term*factork
                  factor = -3.0d0 * (dotkr2 + yrr2*part)
                  factork = 1.0d0 - yk*ykrk2
                  dtydyi1 = dtdyi1*termy + term*factor
                  dtykdyi1 = dtdyi1*termyk + term*factork
                  factor = -3.0d0 * zrr2 * part
                  factork = -zk * ykrk2
                  dtzdyi1 = dtdyi1*termz + term*factor
                  dtzkdyi1 = dtdyi1*termzk + term*factork
c
                  dtdzi1 = -3.0d0 * zrr2
                  part = zk - dotk2*zrr2
                  factor = -3.0d0 * xrr2 * part
                  factork = -xk * zkrk2
                  dtxdzi1 = dtdzi1*termx + term*factor
                  dtxkdzi1 = dtdzi1*termxk + term*factork
                  factor = -3.0d0 * yrr2 * part
                  factork = -yk * zkrk2
                  dtydzi1 = dtdzi1*termy + term*factor
                  dtykdzi1 = dtdzi1*termyk + term*factork
                  factor = -3.0d0 * (dotkr2 + zrr2*part)
                  factork = 1.0d0 - zk*zkrk2
                  dtzdzi1 = dtdzi1*termz + term*factor
                  dtzkdzi1 = dtdzi1*termzk + term*factork
c
c     now, increment diagonal and off-diagonal Hessian elements
c
                  hessx(1,i1) = hessx(1,i1) + dtxdxi1
                  hessx(2,i1) = hessx(2,i1) + dtydxi1
                  hessx(3,i1) = hessx(3,i1) + dtzdxi1
                  hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi1 + dtxkdxi1
                  hessx(2,k1) = hessx(2,k1) - sk1*dtydxi1 + dtykdxi1
                  hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi1 + dtzkdxi1
                  hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi1 - dtxkdxi1
                  hessx(2,k2) = hessx(2,k2) - sk2*dtydxi1 - dtykdxi1
                  hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi1 - dtzkdxi1
                  hessy(1,i1) = hessy(1,i1) + dtxdyi1
                  hessy(2,i1) = hessy(2,i1) + dtydyi1
                  hessy(3,i1) = hessy(3,i1) + dtzdyi1
                  hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi1 + dtxkdyi1
                  hessy(2,k1) = hessy(2,k1) - sk1*dtydyi1 + dtykdyi1
                  hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi1 + dtzkdyi1
                  hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi1 - dtxkdyi1
                  hessy(2,k2) = hessy(2,k2) - sk2*dtydyi1 - dtykdyi1
                  hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi1 - dtzkdyi1
                  hessz(1,i1) = hessz(1,i1) + dtxdzi1
                  hessz(2,i1) = hessz(2,i1) + dtydzi1
                  hessz(3,i1) = hessz(3,i1) + dtzdzi1
                  hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi1 + dtxkdzi1
                  hessz(2,k1) = hessz(2,k1) - sk1*dtydzi1 + dtykdzi1
                  hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi1 + dtzkdzi1
                  hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi1 - dtxkdzi1
                  hessz(2,k2) = hessz(2,k2) - sk2*dtydzi1 - dtykdzi1
                  hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi1 - dtzkdzi1
c
c     more energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     hessx(1,i1) = hessx(1,i1) + dtaperx*dedxi1
     &                             + dtaperx*dedxi1 + d2taperxx
                     hessx(2,i1) = hessx(2,i1) + dtaperx*dedyi1
     &                             + dtapery*dedxi1 + d2taperxy
                     hessx(3,i1) = hessx(3,i1) + dtaperx*dedzi1
     &                             + dtaperz*dedxi1 + d2taperxz
                     hessx(1,k1) = hessx(1,k1) + dtaperx*dedxk1
     &                             - sk1*(dtaperx*dedxi1+d2taperxx)
                     hessx(2,k1) = hessx(2,k1) + dtaperx*dedyk1
     &                             - sk1*(dtapery*dedxi1+d2taperxy)
                     hessx(3,k1) = hessx(3,k1) + dtaperx*dedzk1
     &                             - sk1*(dtaperz*dedxi1+d2taperxz)
                     hessx(1,k2) = hessx(1,k2) + dtaperx*dedxk2
     &                             - sk2*(dtaperx*dedxi1+d2taperxx)
                     hessx(2,k2) = hessx(2,k2) + dtaperx*dedyk2
     &                             - sk2*(dtapery*dedxi1+d2taperxy)
                     hessx(3,k2) = hessx(3,k2) + dtaperx*dedzk2
     &                             - sk2*(dtaperz*dedxi1+d2taperxz)
                     hessy(1,i1) = hessy(1,i1) + dtapery*dedxi1
     &                             + dtaperx*dedyi1 + d2taperxy
                     hessy(2,i1) = hessy(2,i1) + dtapery*dedyi1
     &                             + dtapery*dedyi1 + d2taperyy
                     hessy(3,i1) = hessy(3,i1) + dtapery*dedzi1
     &                             + dtaperz*dedyi1 + d2taperyz
                     hessy(1,k1) = hessy(1,k1) + dtapery*dedxk1
     &                             - sk1*(dtaperx*dedyi1+d2taperxy)
                     hessy(2,k1) = hessy(2,k1) + dtapery*dedyk1
     &                             - sk1*(dtapery*dedyi1+d2taperyy)
                     hessy(3,k1) = hessy(3,k1) + dtapery*dedzk1
     &                             - sk1*(dtaperz*dedyi1+d2taperyz)
                     hessy(1,k2) = hessy(1,k2) + dtapery*dedxk2
     &                             - sk2*(dtaperx*dedyi1+d2taperxy)
                     hessy(2,k2) = hessy(2,k2) + dtapery*dedyk2
     &                             - sk2*(dtapery*dedyi1+d2taperyy)
                     hessy(3,k2) = hessy(3,k2) + dtapery*dedzk2
     &                             - sk2*(dtaperz*dedyi1+d2taperyz)
                     hessz(1,i1) = hessz(1,i1) + dtaperz*dedxi1
     &                             + dtaperx*dedzi1 + d2taperxz
                     hessz(2,i1) = hessz(2,i1) + dtaperz*dedyi1
     &                             + dtapery*dedzi1 + d2taperyz
                     hessz(3,i1) = hessz(3,i1) + dtaperz*dedzi1
     &                             + dtaperz*dedzi1 + d2taperzz
                     hessz(1,k1) = hessz(1,k1) + dtaperz*dedxk1
     &                             - sk1*(dtaperx*dedzi1+d2taperxz)
                     hessz(2,k1) = hessz(2,k1) + dtaperz*dedyk1
     &                             - sk1*(dtapery*dedzi1+d2taperyz)
                     hessz(3,k1) = hessz(3,k1) + dtaperz*dedzk1
     &                             - sk1*(dtaperz*dedzi1+d2taperzz)
                     hessz(1,k2) = hessz(1,k2) + dtaperz*dedxk2
     &                             - sk2*(dtaperx*dedzi1+d2taperxz)
                     hessz(2,k2) = hessz(2,k2) + dtaperz*dedyk2
     &                             - sk2*(dtapery*dedzi1+d2taperyz)
                     hessz(3,k2) = hessz(3,k2) + dtaperz*dedzk2
     &                             - sk2*(dtaperz*dedzi1+d2taperzz)
                  end if
               end if
            end if
         end do
   10    continue
      end do
c
c     now, see if the atom of interest is part of a dipole
c
      do k = 1, ndipole
         k1 = idpl(1,k)
         k2 = idpl(2,k)
         if (k1.ne.i .and. k2.ne.i)  goto 20
         do ii = 1, n12(k1)
            omit(i12(ii,k1)) = k
         end do
         do ii = 1, n12(k2)
            omit(i12(ii,k2)) = k
         end do
         sk1 = 1.0d0 - sdpl(k)
         sk2 = sdpl(k)
         xk = x(k1) - x(k2)
         yk = y(k1) - y(k2)
         zk = z(k1) - z(k2)
         rk2 = xk*xk + yk*yk + zk*zk
         xq = x(k1) - xk*sk2
         yq = y(k1) - yk*sk2
         zq = z(k1) - zk*sk2
         fk = -f * bdpl(k)
         do ii = 1, nion
            i1 = iion(ii)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
            if (proceed)  proceed = (omit(i1) .ne. k)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(i1) - xq
               yr = y(i1) - yq
               zr = z(i1) - zq
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fk * pchg(ii)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  fik = fik * fgrp
c
c     some abbreviations used in various chain rule terms
c
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  dotk2 = 2.0d0 * dotk
                  dotkr2 = dotk / r2
                  dotkrk2 = dotk / rk2
c
c     now, form the chain rule terms for first derivatives
c
                  term = fik / rkr3
                  term2 = -3.0d0 * dotk
                  termx = term * (xk+xrr2*term2)
                  termy = term * (yk+yrr2*term2)
                  termz = term * (zk+zrr2*term2)
                  termxk = term * (xr-dotk*xkrk2)
                  termyk = term * (yr-dotk*ykrk2)
                  termzk = term * (zr-dotk*zkrk2)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik * dotk / rkr3
                     dedxi1 = termx
                     dedyi1 = termy
                     dedzi1 = termz
                     dedxk1 = -sk1*termx + termxk
                     dedyk1 = -sk1*termy + termyk
                     dedzk1 = -sk1*termz + termzk
                     dedxk2 = -sk2*termx - termxk
                     dedyk2 = -sk2*termy - termyk
                     dedzk2 = -sk2*termz - termzk
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                            + 6.0d0*c3*r + 2.0d0*c2
                     dtaper = dtaper / r
                     dtaperx = xr * dtaper
                     dtapery = yr * dtaper
                     dtaperz = zr * dtaper
                     d2taper = e * (d2taper-dtaper)
                     dtaper = e * dtaper
                     d2taperxx = xr*xrr2*d2taper + dtaper
                     d2taperxy = xr*yrr2*d2taper
                     d2taperxz = xr*zrr2*d2taper
                     d2taperyy = yr*yrr2*d2taper + dtaper
                     d2taperyz = yr*zrr2*d2taper
                     d2taperzz = zr*zrr2*d2taper + dtaper
                     term = term * taper
                     termx = termx * taper
                     termy = termy * taper
                     termz = termz * taper
                     termxk = termxk * taper
                     termyk = termyk * taper
                     termzk = termzk * taper
                  end if
c
c     next, find the second derivative chain rule terms
c
                  if (k1 .eq. i) then
                     dtdxk1 = 3.0d0*sk1*xrr2 - xkrk2
                     part = sk1*xk - xr
                     part2 = sk1*dotk2*xrr2 - part
                     factor = 1.0d0 - 3.0d0*xrr2*part2
     &                           + 3.0d0*sk1*dotkr2
                     factork = -sk1 + dotk2*xkrk2*xkrk2
     &                            + xkrk2*part - dotkrk2
                     dtxdxk1 = dtdxk1*termx + term*factor
                     dtxkdxk1 = dtdxk1*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part2
                     factork = dotk2*ykrk2*xkrk2 + ykrk2*part
                     dtydxk1 = dtdxk1*termy + term*factor
                     dtykdxk1 = dtdxk1*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part2
                     factork = dotk2*zkrk2*xkrk2 + zkrk2*part
                     dtzdxk1 = dtdxk1*termz + term*factor
                     dtzkdxk1 = dtdxk1*termzk + term*factork
c
                     dtdyk1 = 3.0d0*sk1*yrr2 - ykrk2
                     part = sk1*yk - yr
                     part2 = sk1*dotk2*yrr2 - part
                     factor = -3.0d0 * xrr2 * part2
                     factork = dotk2*xkrk2*ykrk2 + xkrk2*part
                     dtxdyk1 = dtdyk1*termx + term*factor
                     dtxkdyk1 = dtdyk1*termxk + term*factork
                     factor = 1.0d0 - 3.0d0*yrr2*part2
     &                           + 3.0d0*sk1*dotkr2
                     factork = -sk1 + dotk2*ykrk2*ykrk2
     &                            + ykrk2*part - dotkrk2
                     dtydyk1 = dtdyk1*termy + term*factor
                     dtykdyk1 = dtdyk1*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part2
                     factork = dotk2*zkrk2*ykrk2 + zkrk2*part
                     dtzdyk1 = dtdyk1*termz + term*factor
                     dtzkdyk1 = dtdyk1*termzk + term*factork
c
                     dtdzk1 = 3.0d0*sk1*zrr2 - zkrk2
                     part = sk1*zk - zr
                     part2 = sk1*dotk2*zrr2 - part
                     factor = -3.0d0 * xrr2 * part2
                     factork = dotk2*xkrk2*zkrk2 + xkrk2*part
                     dtxdzk1 = dtdzk1*termx + term*factor
                     dtxkdzk1 = dtdzk1*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part2
                     factork = dotk2*ykrk2*zkrk2 + ykrk2*part
                     dtydzk1 = dtdzk1*termy + term*factor
                     dtykdzk1 = dtdzk1*termyk + term*factork
                     factor = 1.0d0 - 3.0d0*zrr2*part2
     &                           + 3.0d0*sk1*dotkr2
                     factork = -sk1 + dotk2*zkrk2*zkrk2
     &                            + zkrk2*part - dotkrk2
                     dtzdzk1 = dtdzk1*termz + term*factor
                     dtzkdzk1 = dtdzk1*termzk + term*factork
c
                  else if (k2 .eq. i) then
                     dtdxk2 = 3.0d0*sk2*xrr2 + xkrk2
                     part = sk2*xk + xr
                     part2 = sk2*dotk2*xrr2 - part
                     factor = -1.0d0 - 3.0d0*xrr2*part2
     &                           + 3.0d0*sk2*dotkr2
                     factork = -sk2 - dotk2*xkrk2*xkrk2
     &                            + xkrk2*part + dotkrk2
                     dtxdxk2 = dtdxk2*termx + term*factor
                     dtxkdxk2 = dtdxk2*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part2
                     factork = -dotk2*ykrk2*xkrk2 + ykrk2*part
                     dtydxk2 = dtdxk2*termy + term*factor
                     dtykdxk2 = dtdxk2*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part2
                     factork = -dotk2*zkrk2*xkrk2 + zkrk2*part
                     dtzdxk2 = dtdxk2*termz + term*factor
                     dtzkdxk2 = dtdxk2*termzk + term*factork
c
                     dtdyk2 = 3.0d0*sk2*yrr2 + ykrk2
                     part = sk2*yk + yr
                     part2 = sk2*dotk2*yrr2 - part
                     factor = -3.0d0 * xrr2 * part2
                     factork = -dotk2*xkrk2*ykrk2 + xkrk2*part
                     dtxdyk2 = dtdyk2*termx + term*factor
                     dtxkdyk2 = dtdyk2*termxk + term*factork
                     factor = -1.0d0 - 3.0d0*yrr2*part2
     &                           + 3.0d0*sk2*dotkr2
                     factork = -sk2 - dotk2*ykrk2*ykrk2
     &                            + ykrk2*part + dotkrk2
                     dtydyk2 = dtdyk2*termy + term*factor
                     dtykdyk2 = dtdyk2*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part2
                     factork = -dotk2*zkrk2*ykrk2 + zkrk2*part
                     dtzdyk2 = dtdyk2*termz + term*factor
                     dtzkdyk2 = dtdyk2*termzk + term*factork
c
                     dtdzk2 = 3.0d0*sk2*zrr2 + zkrk2
                     part = sk2*zk + zr
                     part2 = sk2*dotk2*zrr2 - part
                     factor = -3.0d0 * xrr2 * part2
                     factork = -dotk2*xkrk2*zkrk2 + xkrk2*part
                     dtxdzk2 = dtdzk2*termx + term*factor
                     dtxkdzk2 = dtdzk2*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part2
                     factork = -dotk2*ykrk2*zkrk2 + ykrk2*part
                     dtydzk2 = dtdzk2*termy + term*factor
                     dtykdzk2 = dtdzk2*termyk + term*factork
                     factor = -1.0d0 - 3.0d0*zrr2*part2
     &                           + 3.0d0*sk2*dotkr2
                     factork = -sk2 - dotk2*zkrk2*zkrk2
     &                            + zkrk2*part + dotkrk2
                     dtzdzk2 = dtdzk2*termz + term*factor
                     dtzkdzk2 = dtdzk2*termzk + term*factork
                  end if
c
c     now, increment diagonal and off-diagonal Hessian elements
c
                  if (i .eq. k1) then
                     hessx(1,i1) = hessx(1,i1) + dtxdxk1
                     hessx(2,i1) = hessx(2,i1) + dtydxk1
                     hessx(3,i1) = hessx(3,i1) + dtzdxk1
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxk1 + dtxkdxk1
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxk1 + dtykdxk1
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxk1 + dtzkdxk1
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxk1 - dtxkdxk1
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxk1 - dtykdxk1
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxk1 - dtzkdxk1
                     hessy(1,i1) = hessy(1,i1) + dtxdyk1
                     hessy(2,i1) = hessy(2,i1) + dtydyk1
                     hessy(3,i1) = hessy(3,i1) + dtzdyk1
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyk1 + dtxkdyk1
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyk1 + dtykdyk1
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyk1 + dtzkdyk1
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyk1 - dtxkdyk1
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyk1 - dtykdyk1
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyk1 - dtzkdyk1
                     hessz(1,i1) = hessz(1,i1) + dtxdzk1
                     hessz(2,i1) = hessz(2,i1) + dtydzk1
                     hessz(3,i1) = hessz(3,i1) + dtzdzk1
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzk1 + dtxkdzk1
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzk1 + dtykdzk1
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzk1 + dtzkdzk1
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzk1 - dtxkdzk1
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzk1 - dtykdzk1
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzk1 - dtzkdzk1
                  else if (i .eq. k2) then
                     hessx(1,i1) = hessx(1,i1) + dtxdxk2
                     hessx(2,i1) = hessx(2,i1) + dtydxk2
                     hessx(3,i1) = hessx(3,i1) + dtzdxk2
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxk2 + dtxkdxk2
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxk2 + dtykdxk2
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxk2 + dtzkdxk2
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxk2 - dtxkdxk2
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxk2 - dtykdxk2
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxk2 - dtzkdxk2
                     hessy(1,i1) = hessy(1,i1) + dtxdyk2
                     hessy(2,i1) = hessy(2,i1) + dtydyk2
                     hessy(3,i1) = hessy(3,i1) + dtzdyk2
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyk2 + dtxkdyk2
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyk2 + dtykdyk2
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyk2 + dtzkdyk2
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyk2 - dtxkdyk2
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyk2 - dtykdyk2
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyk2 - dtzkdyk2
                     hessz(1,i1) = hessz(1,i1) + dtxdzk2
                     hessz(2,i1) = hessz(2,i1) + dtydzk2
                     hessz(3,i1) = hessz(3,i1) + dtzdzk2
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzk2 + dtxkdzk2
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzk2 + dtykdzk2
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzk2 + dtzkdzk2
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzk2 - dtxkdzk2
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzk2 - dtykdzk2
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzk2 - dtzkdzk2
                  end if
c
c     more energy switching if near the cutoff distance
c
                  if (r2.gt.cut2 .and. i.eq.k1) then
                     hessx(1,i1) = hessx(1,i1) - sk1*dtaperx*dedxi1
     &                      + dtaperx*dedxk1 - sk1*d2taperxx
                     hessx(2,i1) = hessx(2,i1) - sk1*dtaperx*dedyi1
     &                      + dtapery*dedxk1 - sk1*d2taperxy
                     hessx(3,i1) = hessx(3,i1) - sk1*dtaperx*dedzi1
     &                      + dtaperz*dedxk1 - sk1*d2taperxz
                     hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxk1
     &                      - sk1*dtaperx*dedxk1 + sk1*sk1*d2taperxx
                     hessx(2,k1) = hessx(2,k1) - sk1*dtaperx*dedyk1
     &                      - sk1*dtapery*dedxk1 + sk1*sk1*d2taperxy
                     hessx(3,k1) = hessx(3,k1) - sk1*dtaperx*dedzk1
     &                      - sk1*dtaperz*dedxk1 + sk1*sk1*d2taperxz
                     hessx(1,k2) = hessx(1,k2) - sk1*dtaperx*dedxk2
     &                      - sk2*dtaperx*dedxk1 + sk1*sk2*d2taperxx
                     hessx(2,k2) = hessx(2,k2) - sk1*dtaperx*dedyk2
     &                      - sk2*dtapery*dedxk1 + sk1*sk2*d2taperxy
                     hessx(3,k2) = hessx(3,k2) - sk1*dtaperx*dedzk2
     &                      - sk2*dtaperz*dedxk1 + sk1*sk2*d2taperxz
                     hessy(1,i1) = hessy(1,i1) - sk1*dtapery*dedxi1
     &                      + dtaperx*dedyk1 - sk1*d2taperxy
                     hessy(2,i1) = hessy(2,i1) - sk1*dtapery*dedyi1
     &                      + dtapery*dedyk1 - sk1*d2taperyy
                     hessy(3,i1) = hessy(3,i1) - sk1*dtapery*dedzi1
     &                      + dtaperz*dedyk1 - sk1*d2taperyz
                     hessy(1,k1) = hessy(1,k1) - sk1*dtapery*dedxk1
     &                      - sk1*dtaperx*dedyk1 + sk1*sk1*d2taperxy
                     hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyk1
     &                      - sk1*dtapery*dedyk1 + sk1*sk1*d2taperyy
                     hessy(3,k1) = hessy(3,k1) - sk1*dtapery*dedzk1
     &                      - sk1*dtaperz*dedyk1 + sk1*sk1*d2taperyz
                     hessy(1,k2) = hessy(1,k2) - sk1*dtapery*dedxk2
     &                      - sk2*dtaperx*dedyk1 + sk1*sk2*d2taperxy
                     hessy(2,k2) = hessy(2,k2) - sk1*dtapery*dedyk2
     &                      - sk2*dtapery*dedyk1 + sk1*sk2*d2taperyy
                     hessy(3,k2) = hessy(3,k2) - sk1*dtapery*dedzk2
     &                      - sk2*dtaperz*dedyk1 + sk1*sk2*d2taperyz
                     hessz(1,i1) = hessz(1,i1) - sk1*dtaperz*dedxi1
     &                      + dtaperx*dedzk1 - sk1*d2taperxz
                     hessz(2,i1) = hessz(2,i1) - sk1*dtaperz*dedyi1
     &                      + dtapery*dedzk1 - sk1*d2taperyz
                     hessz(3,i1) = hessz(3,i1) - sk1*dtaperz*dedzi1
     &                      + dtaperz*dedzk1 - sk1*d2taperzz
                     hessz(1,k1) = hessz(1,k1) - sk1*dtaperz*dedxk1
     &                      - sk1*dtaperx*dedzk1 + sk1*sk1*d2taperxz
                     hessz(2,k1) = hessz(2,k1) - sk1*dtaperz*dedyk1
     &                      - sk1*dtapery*dedzk1 + sk1*sk1*d2taperyz
                     hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzk1
     &                      - sk1*dtaperz*dedzk1 + sk1*sk1*d2taperzz
                     hessz(1,k2) = hessz(1,k2) - sk1*dtaperz*dedxk2
     &                      - sk2*dtaperx*dedzk1 + sk1*sk2*d2taperxz
                     hessz(2,k2) = hessz(2,k2) - sk1*dtaperz*dedyk2
     &                      - sk2*dtapery*dedzk1 + sk1*sk2*d2taperyz
                     hessz(3,k2) = hessz(3,k2) - sk1*dtaperz*dedzk2
     &                      - sk2*dtaperz*dedzk1 + sk1*sk2*d2taperzz
                  else if (r2.gt.cut2 .and. i.eq.k2) then
                     hessx(1,i1) = hessx(1,i1) - sk2*dtaperx*dedxi1
     &                      + dtaperx*dedxk2 - sk2*d2taperxx
                     hessx(2,i1) = hessx(2,i1) - sk2*dtaperx*dedyi1
     &                      + dtapery*dedxk2 - sk2*d2taperxy
                     hessx(3,i1) = hessx(3,i1) - sk2*dtaperx*dedzi1
     &                      + dtaperz*dedxk2 - sk2*d2taperxz
                     hessx(1,k1) = hessx(1,k1) - sk2*dtaperx*dedxk1
     &                      - sk1*dtaperx*dedxk2 + sk1*sk2*d2taperxx
                     hessx(2,k1) = hessx(2,k1) - sk2*dtaperx*dedyk1
     &                      - sk1*dtapery*dedxk2 + sk1*sk2*d2taperxy
                     hessx(3,k1) = hessx(3,k1) - sk2*dtaperx*dedzk1
     &                      - sk1*dtaperz*dedxk2 + sk1*sk2*d2taperxz
                     hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxk2
     &                      - sk2*dtaperx*dedxk2 + sk2*sk2*d2taperxx
                     hessx(2,k2) = hessx(2,k2) - sk2*dtaperx*dedyk2
     &                      - sk2*dtapery*dedxk2 + sk2*sk2*d2taperxy
                     hessx(3,k2) = hessx(3,k2) - sk2*dtaperx*dedzk2
     &                      - sk2*dtaperz*dedxk2 + sk2*sk2*d2taperxz
                     hessy(1,i1) = hessy(1,i1) - sk2*dtapery*dedxi1
     &                      + dtaperx*dedyk2 - sk2*d2taperxy
                     hessy(2,i1) = hessy(2,i1) - sk2*dtapery*dedyi1
     &                      + dtapery*dedyk2 - sk2*d2taperyy
                     hessy(3,i1) = hessy(3,i1) - sk2*dtapery*dedzi1
     &                      + dtaperz*dedyk2 - sk2*d2taperyz
                     hessy(1,k1) = hessy(1,k1) - sk2*dtapery*dedxk1
     &                      - sk1*dtaperx*dedyk2 + sk1*sk2*d2taperxy
                     hessy(2,k1) = hessy(2,k1) - sk2*dtapery*dedyk1
     &                      - sk1*dtapery*dedyk2 + sk1*sk2*d2taperyy
                     hessy(3,k1) = hessy(3,k1) - sk2*dtapery*dedzk1
     &                      - sk1*dtaperz*dedyk2 + sk1*sk2*d2taperyz
                     hessy(1,k2) = hessy(1,k2) - sk2*dtapery*dedxk2
     &                      - sk2*dtaperx*dedyk2 + sk2*sk2*d2taperxy
                     hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyk2
     &                      - sk2*dtapery*dedyk2 + sk2*sk2*d2taperyy
                     hessy(3,k2) = hessy(3,k2) - sk2*dtapery*dedzk2
     &                      - sk2*dtaperz*dedyk2 + sk2*sk2*d2taperyz
                     hessz(1,i1) = hessz(1,i1) - sk2*dtaperz*dedxi1
     &                      + dtaperx*dedzk2 - sk2*d2taperxz
                     hessz(2,i1) = hessz(2,i1) - sk2*dtaperz*dedyi1
     &                      + dtapery*dedzk2 - sk2*d2taperyz
                     hessz(3,i1) = hessz(3,i1) - sk2*dtaperz*dedzi1
     &                      + dtaperz*dedzk2 - sk2*d2taperzz
                     hessz(1,k1) = hessz(1,k1) - sk2*dtaperz*dedxk1
     &                      - sk1*dtaperx*dedzk2 + sk1*sk2*d2taperxz
                     hessz(2,k1) = hessz(2,k1) - sk2*dtaperz*dedyk1
     &                      - sk1*dtapery*dedzk2 + sk1*sk2*d2taperyz
                     hessz(3,k1) = hessz(3,k1) - sk2*dtaperz*dedzk1
     &                      - sk1*dtaperz*dedzk2 + sk1*sk2*d2taperzz
                     hessz(1,k2) = hessz(1,k2) - sk2*dtaperz*dedxk2
     &                      - sk2*dtaperx*dedzk2 + sk2*sk2*d2taperxz
                     hessz(2,k2) = hessz(2,k2) - sk2*dtaperz*dedyk2
     &                      - sk2*dtapery*dedzk2 + sk2*sk2*d2taperyz
                     hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzk2
     &                      - sk2*dtaperz*dedzk2 + sk2*sk2*d2taperzz
                  end if
               end if
            end if
         end do
   20    continue
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion
         i1 = iion(ii)
         if (i1 .ne. i)  goto 30
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(ii)
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               do jcell = 1, ncell
                  xk = x(k1) - x(k2)
                  yk = y(k1) - y(k2)
                  zk = z(k1) - z(k2)
                  xr = xi - x(k1) + xk*sk2
                  yr = yi - y(k1) + yk*sk2
                  zr = zi - z(k1) + zk*sk2
                  call image (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rk2 = xk*xk + yk*yk + zk*zk
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = -fi * bdpl(k)
c
c     scale the interaction based on its group membership
c
                     if (use_group)  fik = fik * fgrp
c
c     some abbreviations used in various chain rule terms
c
                     xrr2 = xr / r2
                     yrr2 = yr / r2
                     zrr2 = zr / r2
                     xkrk2 = xk / rk2
                     ykrk2 = yk / rk2
                     zkrk2 = zk / rk2
                     dotk2 = 2.0d0 * dotk
                     dotkr2 = dotk / r2
c
c     now, form the chain rule terms for first derivatives
c
                     term = fik / rkr3
                     term2 = -3.0d0 * dotk
                     termx = term * (xk+xrr2*term2)
                     termy = term * (yk+yrr2*term2)
                     termz = term * (zk+zrr2*term2)
                     termxk = term * (xr-dotk*xkrk2)
                     termyk = term * (yr-dotk*ykrk2)
                     termzk = term * (zr-dotk*zkrk2)
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        e = fik * dotk / rkr3
                        dedxi1 = termx
                        dedyi1 = termy
                        dedzi1 = termz
                        dedxk1 = -sk1*termx + termxk
                        dedyk1 = -sk1*termy + termyk
                        dedzk1 = -sk1*termz + termzk
                        dedxk2 = -sk2*termx - termxk
                        dedyk2 = -sk2*termy - termyk
                        dedzk2 = -sk2*termz - termzk
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                               + 6.0d0*c3*r + 2.0d0*c2
                        dtaper = dtaper / r
                        dtaperx = xr * dtaper
                        dtapery = yr * dtaper
                        dtaperz = zr * dtaper
                        d2taper = e * (d2taper-dtaper)
                        dtaper = e * dtaper
                        d2taperxx = xr*xrr2*d2taper + dtaper
                        d2taperxy = xr*yrr2*d2taper
                        d2taperxz = xr*zrr2*d2taper
                        d2taperyy = yr*yrr2*d2taper + dtaper
                        d2taperyz = yr*zrr2*d2taper
                        d2taperzz = zr*zrr2*d2taper + dtaper
                        term = term * taper
                        termx = termx * taper
                        termy = termy * taper
                        termz = termz * taper
                        termxk = termxk * taper
                        termyk = termyk * taper
                        termzk = termzk * taper
                     end if
c
c     next, find the second derivative chain rule terms
c
                     dtdxi1 = -3.0d0 * xrr2
                     part = xk - dotk2*xrr2
                     factor = -3.0d0 * (dotkr2 + xrr2*part)
                     factork = 1.0d0 - xk*xkrk2
                     dtxdxi1 = dtdxi1*termx + term*factor
                     dtxkdxi1 = dtdxi1*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part
                     factork = -yk * xkrk2
                     dtydxi1 = dtdxi1*termy + term*factor
                     dtykdxi1 = dtdxi1*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part
                     factork = -zk * xkrk2
                     dtzdxi1 = dtdxi1*termz + term*factor
                     dtzkdxi1 = dtdxi1*termzk + term*factork
c
                     dtdyi1 = -3.0d0 * yrr2
                     part = yk - dotk2*yrr2
                     factor = -3.0d0 * xrr2 * part
                     factork = -xk * ykrk2
                     dtxdyi1 = dtdyi1*termx + term*factor
                     dtxkdyi1 = dtdyi1*termxk + term*factork
                     factor = -3.0d0 * (dotkr2 + yrr2*part)
                     factork = 1.0d0 - yk*ykrk2
                     dtydyi1 = dtdyi1*termy + term*factor
                     dtykdyi1 = dtdyi1*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part
                     factork = -zk * ykrk2
                     dtzdyi1 = dtdyi1*termz + term*factor
                     dtzkdyi1 = dtdyi1*termzk + term*factork
c
                     dtdzi1 = -3.0d0 * zrr2
                     part = zk - dotk2*zrr2
                     factor = -3.0d0 * xrr2 * part
                     factork = -xk * zkrk2
                     dtxdzi1 = dtdzi1*termx + term*factor
                     dtxkdzi1 = dtdzi1*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part
                     factork = -yk * zkrk2
                     dtydzi1 = dtdzi1*termy + term*factor
                     dtykdzi1 = dtdzi1*termyk + term*factork
                     factor = -3.0d0 * (dotkr2 + zrr2*part)
                     factork = 1.0d0 - zk*zkrk2
                     dtzdzi1 = dtdzi1*termz + term*factor
                     dtzkdzi1 = dtdzi1*termzk + term*factork
c
c     now, increment diagonal and off-diagonal Hessian elements
c
                     hessx(1,i1) = hessx(1,i1) + dtxdxi1
                     hessx(2,i1) = hessx(2,i1) + dtydxi1
                     hessx(3,i1) = hessx(3,i1) + dtzdxi1
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi1 + dtxkdxi1
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi1 + dtykdxi1
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi1 + dtzkdxi1
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi1 - dtxkdxi1
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi1 - dtykdxi1
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi1 - dtzkdxi1
                     hessy(1,i1) = hessy(1,i1) + dtxdyi1
                     hessy(2,i1) = hessy(2,i1) + dtydyi1
                     hessy(3,i1) = hessy(3,i1) + dtzdyi1
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi1 + dtxkdyi1
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi1 + dtykdyi1
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi1 + dtzkdyi1
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi1 - dtxkdyi1
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi1 - dtykdyi1
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi1 - dtzkdyi1
                     hessz(1,i1) = hessz(1,i1) + dtxdzi1
                     hessz(2,i1) = hessz(2,i1) + dtydzi1
                     hessz(3,i1) = hessz(3,i1) + dtzdzi1
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi1 + dtxkdzi1
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi1 + dtykdzi1
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi1 + dtzkdzi1
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi1 - dtxkdzi1
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi1 - dtykdzi1
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi1 - dtzkdzi1
c
c     more energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        hessx(1,i1) = hessx(1,i1) + dtaperx*dedxi1
     &                                + dtaperx*dedxi1 + d2taperxx
                        hessx(2,i1) = hessx(2,i1) + dtaperx*dedyi1
     &                                + dtapery*dedxi1 + d2taperxy
                        hessx(3,i1) = hessx(3,i1) + dtaperx*dedzi1
     &                                + dtaperz*dedxi1 + d2taperxz
                        hessx(1,k1) = hessx(1,k1) + dtaperx*dedxk1
     &                                - sk1*(dtaperx*dedxi1+d2taperxx)
                        hessx(2,k1) = hessx(2,k1) + dtaperx*dedyk1
     &                                - sk1*(dtapery*dedxi1+d2taperxy)
                        hessx(3,k1) = hessx(3,k1) + dtaperx*dedzk1
     &                                - sk1*(dtaperz*dedxi1+d2taperxz)
                        hessx(1,k2) = hessx(1,k2) + dtaperx*dedxk2
     &                                - sk2*(dtaperx*dedxi1+d2taperxx)
                        hessx(2,k2) = hessx(2,k2) + dtaperx*dedyk2
     &                                - sk2*(dtapery*dedxi1+d2taperxy)
                        hessx(3,k2) = hessx(3,k2) + dtaperx*dedzk2
     &                                - sk2*(dtaperz*dedxi1+d2taperxz)
                        hessy(1,i1) = hessy(1,i1) + dtapery*dedxi1
     &                                + dtaperx*dedyi1 + d2taperxy
                        hessy(2,i1) = hessy(2,i1) + dtapery*dedyi1
     &                                + dtapery*dedyi1 + d2taperyy
                        hessy(3,i1) = hessy(3,i1) + dtapery*dedzi1
     &                                + dtaperz*dedyi1 + d2taperyz
                        hessy(1,k1) = hessy(1,k1) + dtapery*dedxk1
     &                                - sk1*(dtaperx*dedyi1+d2taperxy)
                        hessy(2,k1) = hessy(2,k1) + dtapery*dedyk1
     &                                - sk1*(dtapery*dedyi1+d2taperyy)
                        hessy(3,k1) = hessy(3,k1) + dtapery*dedzk1
     &                                - sk1*(dtaperz*dedyi1+d2taperyz)
                        hessy(1,k2) = hessy(1,k2) + dtapery*dedxk2
     &                                - sk2*(dtaperx*dedyi1+d2taperxy)
                        hessy(2,k2) = hessy(2,k2) + dtapery*dedyk2
     &                                - sk2*(dtapery*dedyi1+d2taperyy)
                        hessy(3,k2) = hessy(3,k2) + dtapery*dedzk2
     &                                - sk2*(dtaperz*dedyi1+d2taperyz)
                        hessz(1,i1) = hessz(1,i1) + dtaperz*dedxi1
     &                                + dtaperx*dedzi1 + d2taperxz
                        hessz(2,i1) = hessz(2,i1) + dtaperz*dedyi1
     &                                + dtapery*dedzi1 + d2taperyz
                        hessz(3,i1) = hessz(3,i1) + dtaperz*dedzi1
     &                                + dtaperz*dedzi1 + d2taperzz
                        hessz(1,k1) = hessz(1,k1) + dtaperz*dedxk1
     &                                - sk1*(dtaperx*dedzi1+d2taperxz)
                        hessz(2,k1) = hessz(2,k1) + dtaperz*dedyk1
     &                                - sk1*(dtapery*dedzi1+d2taperyz)
                        hessz(3,k1) = hessz(3,k1) + dtaperz*dedzk1
     &                                - sk1*(dtaperz*dedzi1+d2taperzz)
                        hessz(1,k2) = hessz(1,k2) + dtaperz*dedxk2
     &                                - sk2*(dtaperx*dedzi1+d2taperxz)
                        hessz(2,k2) = hessz(2,k2) + dtaperz*dedyk2
     &                                - sk2*(dtapery*dedzi1+d2taperyz)
                        hessz(3,k2) = hessz(3,k2) + dtaperz*dedzk2
     &                                - sk2*(dtaperz*dedzi1+d2taperzz)
                     end if
                  end if
               end do
            end if
         end do
   30    continue
      end do
c
c     now, see if the atom of interest is part of a dipole
c
      do k = 1, ndipole
         k1 = idpl(1,k)
         k2 = idpl(2,k)
         if (k1.ne.i .and. k2.ne.i)  goto 40
         sk1 = 1.0d0 - sdpl(k)
         sk2 = sdpl(k)
         xk = x(k1) - x(k2)
         yk = y(k1) - y(k2)
         zk = z(k1) - z(k2)
         rk2 = xk*xk + yk*yk + zk*zk
         xq = x(k1) - xk*sk2
         yq = y(k1) - yk*sk2
         zq = z(k1) - zk*sk2
         fk = -f * bdpl(k)
         do ii = 1, nion
            i1 = iion(ii)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do jcell = 1, ncell
                  xr = x(i1) - xq
                  yr = y(i1) - yq
                  zr = z(i1) - zq
                  call image (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = fk * pchg(ii)
c
c     scale the interaction based on its group membership
c
                     if (use_group)  fik = fik * fgrp
c
c     some abbreviations used in various chain rule terms
c
                     xrr2 = xr / r2
                     yrr2 = yr / r2
                     zrr2 = zr / r2
                     xkrk2 = xk / rk2
                     ykrk2 = yk / rk2
                     zkrk2 = zk / rk2
                     dotk2 = 2.0d0 * dotk
                     dotkr2 = dotk / r2
                     dotkrk2 = dotk / rk2
c
c     now, form the chain rule terms for first derivatives
c
                     term = fik / rkr3
                     term2 = -3.0d0 * dotk
                     termx = term * (xk+xrr2*term2)
                     termy = term * (yk+yrr2*term2)
                     termz = term * (zk+zrr2*term2)
                     termxk = term * (xr-dotk*xkrk2)
                     termyk = term * (yr-dotk*ykrk2)
                     termzk = term * (zr-dotk*zkrk2)
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        e = fik * dotk / rkr3
                        dedxi1 = termx
                        dedyi1 = termy
                        dedzi1 = termz
                        dedxk1 = -sk1*termx + termxk
                        dedyk1 = -sk1*termy + termyk
                        dedzk1 = -sk1*termz + termzk
                        dedxk2 = -sk2*termx - termxk
                        dedyk2 = -sk2*termy - termyk
                        dedzk2 = -sk2*termz - termzk
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                               + 6.0d0*c3*r + 2.0d0*c2
                        dtaper = dtaper / r
                        dtaperx = xr * dtaper
                        dtapery = yr * dtaper
                        dtaperz = zr * dtaper
                        d2taper = e * (d2taper-dtaper)
                        dtaper = e * dtaper
                        d2taperxx = xr*xrr2*d2taper + dtaper
                        d2taperxy = xr*yrr2*d2taper
                        d2taperxz = xr*zrr2*d2taper
                        d2taperyy = yr*yrr2*d2taper + dtaper
                        d2taperyz = yr*zrr2*d2taper
                        d2taperzz = zr*zrr2*d2taper + dtaper
                        term = term * taper
                        termx = termx * taper
                        termy = termy * taper
                        termz = termz * taper
                        termxk = termxk * taper
                        termyk = termyk * taper
                        termzk = termzk * taper
                     end if
c
c     next, find the second derivative chain rule terms
c
                     if (k1 .eq. i) then
                        dtdxk1 = 3.0d0*sk1*xrr2 - xkrk2
                        part = sk1*xk - xr
                        part2 = sk1*dotk2*xrr2 - part
                        factor = 1.0d0 - 3.0d0*xrr2*part2
     &                              + 3.0d0*sk1*dotkr2
                        factork = -sk1 + dotk2*xkrk2*xkrk2
     &                               + xkrk2*part - dotkrk2
                        dtxdxk1 = dtdxk1*termx + term*factor
                        dtxkdxk1 = dtdxk1*termxk + term*factork
                        factor = -3.0d0 * yrr2 * part2
                        factork = dotk2*ykrk2*xkrk2 + ykrk2*part
                        dtydxk1 = dtdxk1*termy + term*factor
                        dtykdxk1 = dtdxk1*termyk + term*factork
                        factor = -3.0d0 * zrr2 * part2
                        factork = dotk2*zkrk2*xkrk2 + zkrk2*part
                        dtzdxk1 = dtdxk1*termz + term*factor
                        dtzkdxk1 = dtdxk1*termzk + term*factork
c
                        dtdyk1 = 3.0d0*sk1*yrr2 - ykrk2
                        part = sk1*yk - yr
                        part2 = sk1*dotk2*yrr2 - part
                        factor = -3.0d0 * xrr2 * part2
                        factork = dotk2*xkrk2*ykrk2 + xkrk2*part
                        dtxdyk1 = dtdyk1*termx + term*factor
                        dtxkdyk1 = dtdyk1*termxk + term*factork
                        factor = 1.0d0 - 3.0d0*yrr2*part2
     &                              + 3.0d0*sk1*dotkr2
                        factork = -sk1 + dotk2*ykrk2*ykrk2
     &                               + ykrk2*part - dotkrk2
                        dtydyk1 = dtdyk1*termy + term*factor
                        dtykdyk1 = dtdyk1*termyk + term*factork
                        factor = -3.0d0 * zrr2 * part2
                        factork = dotk2*zkrk2*ykrk2 + zkrk2*part
                        dtzdyk1 = dtdyk1*termz + term*factor
                        dtzkdyk1 = dtdyk1*termzk + term*factork
c
                        dtdzk1 = 3.0d0*sk1*zrr2 - zkrk2
                        part = sk1*zk - zr
                        part2 = sk1*dotk2*zrr2 - part
                        factor = -3.0d0 * xrr2 * part2
                        factork = dotk2*xkrk2*zkrk2 + xkrk2*part
                        dtxdzk1 = dtdzk1*termx + term*factor
                        dtxkdzk1 = dtdzk1*termxk + term*factork
                        factor = -3.0d0 * yrr2 * part2
                        factork = dotk2*ykrk2*zkrk2 + ykrk2*part
                        dtydzk1 = dtdzk1*termy + term*factor
                        dtykdzk1 = dtdzk1*termyk + term*factork
                        factor = 1.0d0 - 3.0d0*zrr2*part2
     &                              + 3.0d0*sk1*dotkr2
                        factork = -sk1 + dotk2*zkrk2*zkrk2
     &                               + zkrk2*part - dotkrk2
                        dtzdzk1 = dtdzk1*termz + term*factor
                        dtzkdzk1 = dtdzk1*termzk + term*factork
c
                     else if (k2 .eq. i) then
                        dtdxk2 = 3.0d0*sk2*xrr2 + xkrk2
                        part = sk2*xk + xr
                        part2 = sk2*dotk2*xrr2 - part
                        factor = -1.0d0 - 3.0d0*xrr2*part2
     &                              + 3.0d0*sk2*dotkr2
                        factork = -sk2 - dotk2*xkrk2*xkrk2
     &                               + xkrk2*part + dotkrk2
                        dtxdxk2 = dtdxk2*termx + term*factor
                        dtxkdxk2 = dtdxk2*termxk + term*factork
                        factor = -3.0d0 * yrr2 * part2
                        factork = -dotk2*ykrk2*xkrk2 + ykrk2*part
                        dtydxk2 = dtdxk2*termy + term*factor
                        dtykdxk2 = dtdxk2*termyk + term*factork
                        factor = -3.0d0 * zrr2 * part2
                        factork = -dotk2*zkrk2*xkrk2 + zkrk2*part
                        dtzdxk2 = dtdxk2*termz + term*factor
                        dtzkdxk2 = dtdxk2*termzk + term*factork
c
                        dtdyk2 = 3.0d0*sk2*yrr2 + ykrk2
                        part = sk2*yk + yr
                        part2 = sk2*dotk2*yrr2 - part
                        factor = -3.0d0 * xrr2 * part2
                        factork = -dotk2*xkrk2*ykrk2 + xkrk2*part
                        dtxdyk2 = dtdyk2*termx + term*factor
                        dtxkdyk2 = dtdyk2*termxk + term*factork
                        factor = -1.0d0 - 3.0d0*yrr2*part2
     &                              + 3.0d0*sk2*dotkr2
                        factork = -sk2 - dotk2*ykrk2*ykrk2
     &                               + ykrk2*part + dotkrk2
                        dtydyk2 = dtdyk2*termy + term*factor
                        dtykdyk2 = dtdyk2*termyk + term*factork
                        factor = -3.0d0 * zrr2 * part2
                        factork = -dotk2*zkrk2*ykrk2 + zkrk2*part
                        dtzdyk2 = dtdyk2*termz + term*factor
                        dtzkdyk2 = dtdyk2*termzk + term*factork
c
                        dtdzk2 = 3.0d0*sk2*zrr2 + zkrk2
                        part = sk2*zk + zr
                        part2 = sk2*dotk2*zrr2 - part
                        factor = -3.0d0 * xrr2 * part2
                        factork = -dotk2*xkrk2*zkrk2 + xkrk2*part
                        dtxdzk2 = dtdzk2*termx + term*factor
                        dtxkdzk2 = dtdzk2*termxk + term*factork
                        factor = -3.0d0 * yrr2 * part2
                        factork = -dotk2*ykrk2*zkrk2 + ykrk2*part
                        dtydzk2 = dtdzk2*termy + term*factor
                        dtykdzk2 = dtdzk2*termyk + term*factork
                        factor = -1.0d0 - 3.0d0*zrr2*part2
     &                              + 3.0d0*sk2*dotkr2
                        factork = -sk2 - dotk2*zkrk2*zkrk2
     &                               + zkrk2*part + dotkrk2
                        dtzdzk2 = dtdzk2*termz + term*factor
                        dtzkdzk2 = dtdzk2*termzk + term*factork
                     end if
c
c     now, increment diagonal and off-diagonal Hessian elements
c
                     if (i .eq. k1) then
                        hessx(1,i1) = hessx(1,i1) + dtxdxk1
                        hessx(2,i1) = hessx(2,i1) + dtydxk1
                        hessx(3,i1) = hessx(3,i1) + dtzdxk1
                        hessx(1,k1) = hessx(1,k1) - sk1*dtxdxk1
     &                                   + dtxkdxk1
                        hessx(2,k1) = hessx(2,k1) - sk1*dtydxk1
     &                                   + dtykdxk1
                        hessx(3,k1) = hessx(3,k1) - sk1*dtzdxk1
     &                                   + dtzkdxk1
                        hessx(1,k2) = hessx(1,k2) - sk2*dtxdxk1
     &                                   - dtxkdxk1
                        hessx(2,k2) = hessx(2,k2) - sk2*dtydxk1
     &                                   - dtykdxk1
                        hessx(3,k2) = hessx(3,k2) - sk2*dtzdxk1
     &                                   - dtzkdxk1
                        hessy(1,i1) = hessy(1,i1) + dtxdyk1
                        hessy(2,i1) = hessy(2,i1) + dtydyk1
                        hessy(3,i1) = hessy(3,i1) + dtzdyk1
                        hessy(1,k1) = hessy(1,k1) - sk1*dtxdyk1
     &                                   + dtxkdyk1
                        hessy(2,k1) = hessy(2,k1) - sk1*dtydyk1
     &                                   + dtykdyk1
                        hessy(3,k1) = hessy(3,k1) - sk1*dtzdyk1
     &                                   + dtzkdyk1
                        hessy(1,k2) = hessy(1,k2) - sk2*dtxdyk1
     &                                   - dtxkdyk1
                        hessy(2,k2) = hessy(2,k2) - sk2*dtydyk1
     &                                   - dtykdyk1
                        hessy(3,k2) = hessy(3,k2) - sk2*dtzdyk1
     &                                   - dtzkdyk1
                        hessz(1,i1) = hessz(1,i1) + dtxdzk1
                        hessz(2,i1) = hessz(2,i1) + dtydzk1
                        hessz(3,i1) = hessz(3,i1) + dtzdzk1
                        hessz(1,k1) = hessz(1,k1) - sk1*dtxdzk1
     &                                   + dtxkdzk1
                        hessz(2,k1) = hessz(2,k1) - sk1*dtydzk1
     &                                   + dtykdzk1
                        hessz(3,k1) = hessz(3,k1) - sk1*dtzdzk1
     &                                   + dtzkdzk1
                        hessz(1,k2) = hessz(1,k2) - sk2*dtxdzk1
     &                                   - dtxkdzk1
                        hessz(2,k2) = hessz(2,k2) - sk2*dtydzk1
     &                                   - dtykdzk1
                        hessz(3,k2) = hessz(3,k2) - sk2*dtzdzk1
     &                                   - dtzkdzk1
                     else if (i .eq. k2) then
                        hessx(1,i1) = hessx(1,i1) + dtxdxk2
                        hessx(2,i1) = hessx(2,i1) + dtydxk2
                        hessx(3,i1) = hessx(3,i1) + dtzdxk2
                        hessx(1,k1) = hessx(1,k1) - sk1*dtxdxk2
     &                                   + dtxkdxk2
                        hessx(2,k1) = hessx(2,k1) - sk1*dtydxk2
     &                                   + dtykdxk2
                        hessx(3,k1) = hessx(3,k1) - sk1*dtzdxk2
     &                                   + dtzkdxk2
                        hessx(1,k2) = hessx(1,k2) - sk2*dtxdxk2
     &                                   - dtxkdxk2
                        hessx(2,k2) = hessx(2,k2) - sk2*dtydxk2
     &                                   - dtykdxk2
                        hessx(3,k2) = hessx(3,k2) - sk2*dtzdxk2
     &                                   - dtzkdxk2
                        hessy(1,i1) = hessy(1,i1) + dtxdyk2
                        hessy(2,i1) = hessy(2,i1) + dtydyk2
                        hessy(3,i1) = hessy(3,i1) + dtzdyk2
                        hessy(1,k1) = hessy(1,k1) - sk1*dtxdyk2
     &                                   + dtxkdyk2
                        hessy(2,k1) = hessy(2,k1) - sk1*dtydyk2
     &                                   + dtykdyk2
                        hessy(3,k1) = hessy(3,k1) - sk1*dtzdyk2
     &                                   + dtzkdyk2
                        hessy(1,k2) = hessy(1,k2) - sk2*dtxdyk2
     &                                   - dtxkdyk2
                        hessy(2,k2) = hessy(2,k2) - sk2*dtydyk2
     &                                   - dtykdyk2
                        hessy(3,k2) = hessy(3,k2) - sk2*dtzdyk2
     &                                   - dtzkdyk2
                        hessz(1,i1) = hessz(1,i1) + dtxdzk2
                        hessz(2,i1) = hessz(2,i1) + dtydzk2
                        hessz(3,i1) = hessz(3,i1) + dtzdzk2
                        hessz(1,k1) = hessz(1,k1) - sk1*dtxdzk2
     &                                   + dtxkdzk2
                        hessz(2,k1) = hessz(2,k1) - sk1*dtydzk2
     &                                   + dtykdzk2
                        hessz(3,k1) = hessz(3,k1) - sk1*dtzdzk2
     &                                   + dtzkdzk2
                        hessz(1,k2) = hessz(1,k2) - sk2*dtxdzk2
     &                                   - dtxkdzk2
                        hessz(2,k2) = hessz(2,k2) - sk2*dtydzk2
     &                                   - dtykdzk2
                        hessz(3,k2) = hessz(3,k2) - sk2*dtzdzk2
     &                                   - dtzkdzk2
                     end if
c
c     more energy switching if near the cutoff distance
c
                     if (r2.gt.cut2 .and. i.eq.k1) then
                        hessx(1,i1) = hessx(1,i1) - sk1*dtaperx*dedxi1
     &                         + dtaperx*dedxk1 - sk1*d2taperxx
                        hessx(2,i1) = hessx(2,i1) - sk1*dtaperx*dedyi1
     &                         + dtapery*dedxk1 - sk1*d2taperxy
                        hessx(3,i1) = hessx(3,i1) - sk1*dtaperx*dedzi1
     &                         + dtaperz*dedxk1 - sk1*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxk1
     &                         - sk1*dtaperx*dedxk1 + sk1*sk1*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtaperx*dedyk1
     &                         - sk1*dtapery*dedxk1 + sk1*sk1*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperx*dedzk1
     &                         - sk1*dtaperz*dedxk1 + sk1*sk1*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk1*dtaperx*dedxk2
     &                         - sk2*dtaperx*dedxk1 + sk1*sk2*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk1*dtaperx*dedyk2
     &                         - sk2*dtapery*dedxk1 + sk1*sk2*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk1*dtaperx*dedzk2
     &                         - sk2*dtaperz*dedxk1 + sk1*sk2*d2taperxz
                        hessy(1,i1) = hessy(1,i1) - sk1*dtapery*dedxi1
     &                         + dtaperx*dedyk1 - sk1*d2taperxy
                        hessy(2,i1) = hessy(2,i1) - sk1*dtapery*dedyi1
     &                         + dtapery*dedyk1 - sk1*d2taperyy
                        hessy(3,i1) = hessy(3,i1) - sk1*dtapery*dedzi1
     &                         + dtaperz*dedyk1 - sk1*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtapery*dedxk1
     &                         - sk1*dtaperx*dedyk1 + sk1*sk1*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyk1
     &                         - sk1*dtapery*dedyk1 + sk1*sk1*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtapery*dedzk1
     &                         - sk1*dtaperz*dedyk1 + sk1*sk1*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk1*dtapery*dedxk2
     &                         - sk2*dtaperx*dedyk1 + sk1*sk2*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk1*dtapery*dedyk2
     &                         - sk2*dtapery*dedyk1 + sk1*sk2*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk1*dtapery*dedzk2
     &                         - sk2*dtaperz*dedyk1 + sk1*sk2*d2taperyz
                        hessz(1,i1) = hessz(1,i1) - sk1*dtaperz*dedxi1
     &                         + dtaperx*dedzk1 - sk1*d2taperxz
                        hessz(2,i1) = hessz(2,i1) - sk1*dtaperz*dedyi1
     &                         + dtapery*dedzk1 - sk1*d2taperyz
                        hessz(3,i1) = hessz(3,i1) - sk1*dtaperz*dedzi1
     &                         + dtaperz*dedzk1 - sk1*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperz*dedxk1
     &                         - sk1*dtaperx*dedzk1 + sk1*sk1*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtaperz*dedyk1
     &                         - sk1*dtapery*dedzk1 + sk1*sk1*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzk1
     &                         - sk1*dtaperz*dedzk1 + sk1*sk1*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk1*dtaperz*dedxk2
     &                         - sk2*dtaperx*dedzk1 + sk1*sk2*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk1*dtaperz*dedyk2
     &                         - sk2*dtapery*dedzk1 + sk1*sk2*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk1*dtaperz*dedzk2
     &                         - sk2*dtaperz*dedzk1 + sk1*sk2*d2taperzz
                     else if (r2.gt.cut2 .and. i.eq.k2) then
                        hessx(1,i1) = hessx(1,i1) - sk2*dtaperx*dedxi1
     &                         + dtaperx*dedxk2 - sk2*d2taperxx
                        hessx(2,i1) = hessx(2,i1) - sk2*dtaperx*dedyi1
     &                         + dtapery*dedxk2 - sk2*d2taperxy
                        hessx(3,i1) = hessx(3,i1) - sk2*dtaperx*dedzi1
     &                         + dtaperz*dedxk2 - sk2*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk2*dtaperx*dedxk1
     &                         - sk1*dtaperx*dedxk2 + sk1*sk2*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk2*dtaperx*dedyk1
     &                         - sk1*dtapery*dedxk2 + sk1*sk2*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk2*dtaperx*dedzk1
     &                         - sk1*dtaperz*dedxk2 + sk1*sk2*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxk2
     &                         - sk2*dtaperx*dedxk2 + sk2*sk2*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtaperx*dedyk2
     &                         - sk2*dtapery*dedxk2 + sk2*sk2*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperx*dedzk2
     &                         - sk2*dtaperz*dedxk2 + sk2*sk2*d2taperxz
                        hessy(1,i1) = hessy(1,i1) - sk2*dtapery*dedxi1
     &                         + dtaperx*dedyk2 - sk2*d2taperxy
                        hessy(2,i1) = hessy(2,i1) - sk2*dtapery*dedyi1
     &                         + dtapery*dedyk2 - sk2*d2taperyy
                        hessy(3,i1) = hessy(3,i1) - sk2*dtapery*dedzi1
     &                         + dtaperz*dedyk2 - sk2*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk2*dtapery*dedxk1
     &                         - sk1*dtaperx*dedyk2 + sk1*sk2*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk2*dtapery*dedyk1
     &                         - sk1*dtapery*dedyk2 + sk1*sk2*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk2*dtapery*dedzk1
     &                         - sk1*dtaperz*dedyk2 + sk1*sk2*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtapery*dedxk2
     &                         - sk2*dtaperx*dedyk2 + sk2*sk2*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyk2
     &                         - sk2*dtapery*dedyk2 + sk2*sk2*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtapery*dedzk2
     &                         - sk2*dtaperz*dedyk2 + sk2*sk2*d2taperyz
                        hessz(1,i1) = hessz(1,i1) - sk2*dtaperz*dedxi1
     &                         + dtaperx*dedzk2 - sk2*d2taperxz
                        hessz(2,i1) = hessz(2,i1) - sk2*dtaperz*dedyi1
     &                         + dtapery*dedzk2 - sk2*d2taperyz
                        hessz(3,i1) = hessz(3,i1) - sk2*dtaperz*dedzi1
     &                         + dtaperz*dedzk2 - sk2*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk2*dtaperz*dedxk1
     &                         - sk1*dtaperx*dedzk2 + sk1*sk2*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk2*dtaperz*dedyk1
     &                         - sk1*dtapery*dedzk2 + sk1*sk2*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk2*dtaperz*dedzk1
     &                         - sk1*dtaperz*dedzk2 + sk1*sk2*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperz*dedxk2
     &                         - sk2*dtaperx*dedzk2 + sk2*sk2*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtaperz*dedyk2
     &                         - sk2*dtapery*dedzk2 + sk2*sk2*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzk2
     &                         - sk2*dtaperz*dedzk2 + sk2*sk2*d2taperzz
                     end if
                  end if
               end do
            end if
         end do
   40    continue
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echgdpl3  --  charge-dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echgdpl3" calculates the charge-dipole interaction energy;
c     also partitions the energy among the atoms
c
c
      subroutine echgdpl3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'dipole.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,i1,k1,k2
      integer skip(maxatm)
      real*8 e,rk2,rkr3,dotk
      real*8 f,fi,fik,fgrp
      real*8 r,r2,r3,r4,r5,taper
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xr,yr,zr
      logical header,huge,proceed
c
c
c     zero out the overall charge-dipole interaction energy
c     and partitioning; set up constants for the calculation
c
      necd = 0
      ecd = 0.0d0
      do i = 1, n
         aecd(i) = 0.0d0
      end do
      if (ndipole.eq.0 .or. nion.eq.0)  return
      header = .true.
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye * dielec)
      call switch ('CHGDPL')
c
c     get the total energy by looping over each charge-dipole pair
c
      do i = 1, nion
         i1 = iion(i)
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
            if (proceed)  proceed = (skip(k1).ne.i1 .and.
     &                                 skip(k2).ne.i1)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = xi - x(k1) - xk*sdpl(k)
               yr = yi - y(k1) - yk*sdpl(k)
               zr = zi - z(k1) - zk*sdpl(k)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  fik = fi * bdpl(k)
                  rk2 = xk*xk + yk*yk + zk*zk
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
                  e = fik * dotk / rkr3
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall charge-dipole energy component
c
                  necd = necd + 1
                  ecd = ecd + e
                  aecd(i1) = aecd(i1) + 0.5d0*e
                  aecd(k1) = aecd(k1) + 0.25d0*e
                  aecd(k2) = aecd(k2) + 0.25d0*e
c
c     print a warning if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 25.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Charge-Dipole',
     &                             ' Interactions :',
     &                          //,' Type',8x,'Charge',11x,'Dipole',
     &                              17x,'Distance',6x,'Energy',/)
                     end if
                     write (iout,20)  i1,name(i1),k1,name(k1),
     &                                k2,name(k2),sqrt(r2),e
   20                format (' Chg-Dpl  ',i5,'-',a3,' / ',i5,'-',
     &                          a3,1x,i5,'-',a3,10x,f10.4,f12.4)
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, nion
         i1 = iion(i)
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,3,i1,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  xr = xi - x(k1) - xk*sdpl(k)
                  yr = yi - y(k1) - yk*sdpl(k)
                  zr = zi - z(k1) - zk*sdpl(k)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     fik = fi * bdpl(k)
                     rk2 = xk*xk + yk*yk + zk*zk
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
                     e = fik * dotk / rkr3
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall charge-dipole energy component
c
                     necd = necd + 1
                     ecd = ecd + e
                     aecd(i1) = aecd(i1) + 0.5d0*e
                     aecd(k1) = aecd(k1) + 0.25d0*e
                     aecd(k2) = aecd(k2) + 0.25d0*e
c
c     print a warning if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 25.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Charge-Dipole',
     &                                ' Interactions :',
     &                             //,' Type',8x,'Charge',11x,'Dipole',
     &                                17x,'Distance',6x,'Energy',/)
                        end if
                        write (iout,40)  i1,name(i1),k1,name(k1),
     &                                   k2,name(k2),sqrt(r2),e
   40                   format (' Chg-Dpl  ',i5,'-',a3,' / ',i5,'-',a3,
     &                          1x,i5,'-',a3,'  (XTAL)  ',f10.4,f12.4)
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
