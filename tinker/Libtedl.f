C 11 May 10 - DGF - parallelise elj1
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine edipole  --  dipole-dipole potential energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "edipole" calculates the dipole-dipole interaction energy
c
c
      subroutine edipole
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'dipole.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,i1,i2,k1,k2
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xq,yq,zq,xr,yr,zr
      real*8 f,fi,fik,fgrp
      real*8 e,ri2,rk2,rirkr3
      real*8 doti,dotk,dotp
      real*8 r,r2,r3,r4,r5,taper
      logical proceed
c
c
c     zero out the overall dipole interaction energy
c     and set up the constants for the calculation
c
      ed = 0.0d0
      if (ndipole .eq. 0)  return
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye**2 * dielec)
      call switch ('DIPOLE')
c
c     calculate the pairwise dipole interaction energy term
c
      do i = 1, ndipole-1
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*sdpl(i)
         yq = y(i1) + yi*sdpl(i)
         zq = z(i1) + zi*sdpl(i)
         fi = f * bdpl(i)
         do k = i+1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,4,i1,i2,k1,k2,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
            if (proceed)  proceed = (k1.ne.i1 .and. k1.ne.i2 .and.
     &                                 k2.ne.i1 .and. k2.ne.i2)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = xq - x(k1) - xk*sdpl(k)
               yr = yq - y(k1) - yk*sdpl(k)
               zr = zq - z(k1) - zk*sdpl(k)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr3 = sqrt(ri2*rk2*r2) * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(k)
                  e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
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
c     increment the overall dipole-dipole energy component
c
                  ed = ed + e
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
      do i = 1, ndipole
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*sdpl(i)
         yq = y(i1) + yi*sdpl(i)
         zq = z(i1) + zi*sdpl(i)
         fi = f * bdpl(i)
         do k = i, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,4,i1,i2,k1,k2,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  xr = xq - x(k1) - xk*sdpl(k)
                  yr = yq - y(k1) - yk*sdpl(k)
                  zr = zq - z(k1) - zk*sdpl(k)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rk2 = xk*xk + yk*yk + zk*zk
                     rirkr3 = sqrt(ri2*rk2*r2) * r2
                     dotp = xi*xk + yi*yk + zi*zk
                     doti = xi*xr + yi*yr + zi*zr
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = fi * bdpl(k)
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
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
c     increment the overall dipole-dipole energy component
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ed = ed + e
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
c     ##  subroutine edipole1  --  dipole-dipole energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "edipole1" calculates the dipole-dipole interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine edipole1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
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
      integer i,j,k,i1,i2,k1,k2
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xq,yq,zq,xr,yr,zr
      real*8 f,fi,fik,fgrp
      real*8 e,r2,ri2,rk2,rirkr3
      real*8 doti,dotk,dotp
      real*8 si1,si2,sk1,sk2
      real*8 de,dedr,dedrirk
      real*8 deddoti,deddotk,deddotp
      real*8 termx,termy,termz
      real*8 dedrirkri2,dedrirkrk2
      real*8 termxi,termyi,termzi
      real*8 termxk,termyk,termzk
      real*8 dedxi1,dedyi1,dedzi1
      real*8 dedxi2,dedyi2,dedzi2
      real*8 dedxk1,dedyk1,dedzk1
      real*8 dedxk2,dedyk2,dedzk2
      real*8 r,r3,r4,r5,taper,dtaper
      real*8 dtaperx,dtapery,dtaperz
      logical proceed
c
c
c     zero out the overall dipole interaction energy and derivs,
c     then set up the constants for the calculation
c
      ed = 0.0d0
      do i = 1, n
         ded(1,i) = 0.0d0
         ded(2,i) = 0.0d0
         ded(3,i) = 0.0d0
      end do
      if (ndipole .eq. 0)  return
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye**2 * dielec)
      call switch ('DIPOLE')
c
c     compute the dipole interaction energy and first derivatives
c
      do i = 1, ndipole-1
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         si1 = 1.0d0 - sdpl(i)
         si2 = sdpl(i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*si2
         yq = y(i1) + yi*si2
         zq = z(i1) + zi*si2
         fi = f * bdpl(i)
         do k = i+1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,4,i1,i2,k1,k2,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
            if (proceed)  proceed = (k1.ne.i1 .and. k1.ne.i2 .and.
     &                                 k2.ne.i1 .and. k2.ne.i2)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = xq - x(k1) - xk*sk2
               yr = yq - y(k1) - yk*sk2
               zr = zq - z(k1) - zk*sk2
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr3 = sqrt(ri2*rk2*r2) * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(k)
c
c     form the energy and master chain rule term for derivatives
c
                  e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
                  de = -fik / (rirkr3*r2)
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     secondary chain rule terms for derivative expressions
c
                  deddotp = -de * r2
                  deddoti = de * 3.0d0*dotk
                  deddotk = de * 3.0d0*doti
                  dedr = de * (3.0d0*dotp-15.0d0*doti*dotk/r2)
                  dedrirk = -e
                  dedrirkri2 = dedrirk / ri2
                  dedrirkrk2 = dedrirk / rk2
c
c     more chain rule terms for derivative expressions
c
                  termx = dedr*xr + deddoti*xi + deddotk*xk
                  termy = dedr*yr + deddoti*yi + deddotk*yk
                  termz = dedr*zr + deddoti*zi + deddotk*zk
                  termxi = dedrirkri2*xi + deddotp*xk + deddoti*xr
                  termyi = dedrirkri2*yi + deddotp*yk + deddoti*yr
                  termzi = dedrirkri2*zi + deddotp*zk + deddoti*zr
                  termxk = dedrirkrk2*xk + deddotp*xi + deddotk*xr
                  termyk = dedrirkrk2*yk + deddotp*yi + deddotk*yr
                  termzk = dedrirkrk2*zk + deddotp*zi + deddotk*zr
c
c     finally, the individual first derivative components
c
                  dedxi1 = si1*termx - termxi
                  dedyi1 = si1*termy - termyi
                  dedzi1 = si1*termz - termzi
                  dedxi2 = si2*termx + termxi
                  dedyi2 = si2*termy + termyi
                  dedzi2 = si2*termz + termzi
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
                     dedxi1 = dedxi1*taper + si1*dtaperx
                     dedyi1 = dedyi1*taper + si1*dtapery
                     dedzi1 = dedzi1*taper + si1*dtaperz
                     dedxi2 = dedxi2*taper + si2*dtaperx
                     dedyi2 = dedyi2*taper + si2*dtapery
                     dedzi2 = dedzi2*taper + si2*dtaperz
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
                  ed = ed + e
                  ded(1,i1) = ded(1,i1) + dedxi1
                  ded(2,i1) = ded(2,i1) + dedyi1
                  ded(3,i1) = ded(3,i1) + dedzi1
                  ded(1,i2) = ded(1,i2) + dedxi2
                  ded(2,i2) = ded(2,i2) + dedyi2
                  ded(3,i2) = ded(3,i2) + dedzi2
                  ded(1,k1) = ded(1,k1) + dedxk1
                  ded(2,k1) = ded(2,k1) + dedyk1
                  ded(3,k1) = ded(3,k1) + dedzk1
                  ded(1,k2) = ded(1,k2) + dedxk2
                  ded(2,k2) = ded(2,k2) + dedyk2
                  ded(3,k2) = ded(3,k2) + dedzk2
c
c     increment the total intermolecular energy
c
                  if (molcule(i1) .ne. molcule(k1)) then
                     einter = einter + e
                  end if
c
c     increment the virial for use in pressure computation
c
                  if (isobaric) then
                     virx = virx + xr*(dedxi1+dedxi2)
                     viry = viry + yr*(dedyi1+dedyi2)
                     virz = virz + zr*(dedzi1+dedzi2)
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
      do i = 1, ndipole
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         si1 = 1.0d0 - sdpl(i)
         si2 = sdpl(i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*si2
         yq = y(i1) + yi*si2
         zq = z(i1) + zi*si2
         fi = f * bdpl(i)
         do k = i, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,4,i1,i2,k1,k2,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
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
                  xr = xq - x(k1) - xk*sk2
                  yr = yq - y(k1) - yk*sk2
                  zr = zq - z(k1) - zk*sk2
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rk2 = xk*xk + yk*yk + zk*zk
                     rirkr3 = sqrt(ri2*rk2*r2) * r2
                     dotp = xi*xk + yi*yk + zi*zk
                     doti = xi*xr + yi*yr + zi*zr
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = fi * bdpl(k)
c
c     form the energy and master chain rule term for derivatives
c
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
                     de = -fik / (rirkr3*r2)
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                     end if
c
c     secondary chain rule terms for derivative expressions
c
                     deddotp = -de * r2
                     deddoti = de * 3.0d0*dotk
                     deddotk = de * 3.0d0*doti
                     dedr = de * (3.0d0*dotp-15.0d0*doti*dotk/r2)
                     dedrirk = -e
                     dedrirkri2 = dedrirk / ri2
                     dedrirkrk2 = dedrirk / rk2
c
c     more chain rule terms for derivative expressions
c
                     termx = dedr*xr + deddoti*xi + deddotk*xk
                     termy = dedr*yr + deddoti*yi + deddotk*yk
                     termz = dedr*zr + deddoti*zi + deddotk*zk
                     termxi = dedrirkri2*xi + deddotp*xk + deddoti*xr
                     termyi = dedrirkri2*yi + deddotp*yk + deddoti*yr
                     termzi = dedrirkri2*zi + deddotp*zk + deddoti*zr
                     termxk = dedrirkrk2*xk + deddotp*xi + deddotk*xr
                     termyk = dedrirkrk2*yk + deddotp*yi + deddotk*yr
                     termzk = dedrirkrk2*zk + deddotp*zi + deddotk*zr
c
c     finally, the individual first derivative components
c
                     dedxi1 = si1*termx - termxi
                     dedyi1 = si1*termy - termyi
                     dedzi1 = si1*termz - termzi
                     dedxi2 = si2*termx + termxi
                     dedyi2 = si2*termy + termyi
                     dedzi2 = si2*termz + termzi
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
                        dedxi1 = dedxi1*taper + si1*dtaperx
                        dedyi1 = dedyi1*taper + si1*dtapery
                        dedzi1 = dedzi1*taper + si1*dtaperz
                        dedxi2 = dedxi2*taper + si2*dtaperx
                        dedyi2 = dedyi2*taper + si2*dtapery
                        dedzi2 = dedzi2*taper + si2*dtaperz
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
                     if (i .eq. k)  e = 0.5d0 * e
                     ed = ed + e
                     ded(1,i1) = ded(1,i1) + dedxi1
                     ded(2,i1) = ded(2,i1) + dedyi1
                     ded(3,i1) = ded(3,i1) + dedzi1
                     ded(1,i2) = ded(1,i2) + dedxi2
                     ded(2,i2) = ded(2,i2) + dedyi2
                     ded(3,i2) = ded(3,i2) + dedzi2
                     if (i .ne. k) then
                        ded(1,k1) = ded(1,k1) + dedxk1
                        ded(2,k1) = ded(2,k1) + dedyk1
                        ded(3,k1) = ded(3,k1) + dedzk1
                        ded(1,k2) = ded(1,k2) + dedxk2
                        ded(2,k2) = ded(2,k2) + dedyk2
                        ded(3,k2) = ded(3,k2) + dedzk2
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*(dedxi1+dedxi2)
                        viry = viry + yr*(dedyi1+dedyi2)
                        virz = virz + zr*(dedzi1+dedzi2)
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edipole2  --  atom-wise dipole-dipole Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edipole2" calculates second derivatives of the
c     dipole-dipole interaction energy for a single atom
c
c
      subroutine edipole2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'dipole.i'
      include 'group.i'
      include 'hessn.i'
      include 'units.i'
      include 'shunt.i'
      integer i,i1,i2,k1,k2
      integer jcell,idipole,kdipole
      real*8 f,fi,fik,fgrp
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xq,yq,zq,xr,yr,zr
      real*8 e,r2,ri2,rk2,rirkr3
      real*8 doti,dotk,dotp
      real*8 si1,si2,sk1,sk2
      real*8 de,dedr,dedrirk
      real*8 deddoti,deddotk,deddotp
      real*8 termx,termy,termz
      real*8 termxi,termyi,termzi
      real*8 termxk,termyk,termzk
      real*8 enum,r2inv,ri2inv,dotik,xrr2,yrr2,zrr2
      real*8 xiri2,yiri2,ziri2
      real*8 xkrk2,ykrk2,zkrk2
      real*8 xixr,xiyr,xizr,yixr,yiyr,yizr,zixr,ziyr,zizr
      real*8 xkxr,xkyr,xkzr,ykxr,ykyr,ykzr,zkxr,zkyr,zkzr
      real*8 xixk,xiyk,xizk,yixk,yiyk,yizk,zixk,ziyk,zizk
      real*8 xrxr,xryr,xrzr,yryr,yrzr,zrzr
      real*8 xidotk,yidotk,zidotk,xkdoti,ykdoti,zkdoti
      real*8 factor,factori,factork,part,partik
      real*8 dedxi1,dedyi1,dedzi1,dedxi2,dedyi2,dedzi2
      real*8 dedxk1,dedyk1,dedzk1,dedxk2,dedyk2,dedzk2
      real*8 dtdxi1,dtdyi1,dtdzi1,dtdxi2,dtdyi2,dtdzi2
      real*8 dtxdxi1,dtxidxi1,dtxkdxi1,dtxdxi2,dtxidxi2,dtxkdxi2
      real*8 dtydxi1,dtyidxi1,dtykdxi1,dtydxi2,dtyidxi2,dtykdxi2
      real*8 dtzdxi1,dtzidxi1,dtzkdxi1,dtzdxi2,dtzidxi2,dtzkdxi2
      real*8 dtxdyi1,dtxidyi1,dtxkdyi1,dtxdyi2,dtxidyi2,dtxkdyi2
      real*8 dtydyi1,dtyidyi1,dtykdyi1,dtydyi2,dtyidyi2,dtykdyi2
      real*8 dtzdyi1,dtzidyi1,dtzkdyi1,dtzdyi2,dtzidyi2,dtzkdyi2
      real*8 dtxdzi1,dtxidzi1,dtxkdzi1,dtxdzi2,dtxidzi2,dtxkdzi2
      real*8 dtydzi1,dtyidzi1,dtykdzi1,dtydzi2,dtyidzi2,dtykdzi2
      real*8 dtzdzi1,dtzidzi1,dtzkdzi1,dtzdzi2,dtzidzi2,dtzkdzi2
      real*8 r,r3,r4,r5,taper,dtaper,d2taper
      real*8 dtaperx,dtapery,dtaperz
      real*8 d2taperxx,d2taperyy,d2taperzz
      real*8 d2taperxy,d2taperxz,d2taperyz
      logical proceed
c
c
c     set conversion factor and switching function coefficients
c
      if (ndipole .eq. 0)  return
      f = electric / (debye**2 * dielec)
      call switch ('DIPOLE')
c
c     calculate the dipole interaction energy Hessian elements
c
      do idipole = 1, ndipole
         i1 = idpl(1,idipole)
         i2 = idpl(2,idipole)
         si1 = 1.0d0 - sdpl(idipole)
         si2 = sdpl(idipole)
         if (i1.ne.i .and. i2.ne.i)  goto 10
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*si2
         yq = y(i1) + yi*si2
         zq = z(i1) + zi*si2
         fi = f * bdpl(idipole)
         do kdipole = 1, ndipole
            k1 = idpl(1,kdipole)
            k2 = idpl(2,kdipole)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,4,i1,i2,k1,k2,0)
            if (proceed)  proceed = (k1.ne.i1 .and. k1.ne.i2 .and.
     &                                 k2.ne.i1 .and. k2.ne.i2)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(kdipole)
               sk2 = sdpl(kdipole)
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = xq - x(k1) - xk*sk2
               yr = yq - y(k1) - yk*sk2
               zr = zq - z(k1) - zk*sk2
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr3 = sqrt(ri2*rk2*r2) * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(kdipole)
c
c     some abbreviations used in various chain rule terms
c
                  dotik = doti * dotk
                  enum = dotp*r2 - 3.0d0*dotik
                  r2inv = 15.0d0 / r2
                  ri2inv = 1.0d0 / ri2
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xiri2 = xi / ri2
                  yiri2 = yi / ri2
                  ziri2 = zi / ri2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  xixr = xi * xr
                  xiyr = xi * yr
                  xizr = xi * zr
                  yixr = yi * xr
                  yiyr = yi * yr
                  yizr = yi * zr
                  zixr = zi * xr
                  ziyr = zi * yr
                  zizr = zi * zr
                  xkxr = xk * xr
                  xkyr = xk * yr
                  xkzr = xk * zr
                  ykxr = yk * xr
                  ykyr = yk * yr
                  ykzr = yk * zr
                  zkxr = zk * xr
                  zkyr = zk * yr
                  zkzr = zk * zr
                  xixk = xi * xk
                  xiyk = xi * yk
                  xizk = xi * zk
                  yixk = yi * xk
                  yiyk = yi * yk
                  yizk = yi * zk
                  zixk = zi * xk
                  ziyk = zi * yk
                  zizk = zi * zk
                  xrxr = 3.0d0 * xr * xr
                  xryr = 3.0d0 * xr * yr
                  xrzr = 3.0d0 * xr * zr
                  yryr = 3.0d0 * yr * yr
                  yrzr = 3.0d0 * yr * zr
                  zrzr = 3.0d0 * zr * zr
                  xidotk = xi * dotk
                  yidotk = yi * dotk
                  zidotk = zi * dotk
                  xkdoti = xk * doti
                  ykdoti = yk * doti
                  zkdoti = zk * doti
c
c     scale the interaction based on its group membership
c
                  if (use_group)  fik = fik * fgrp
c
c     form the master chain rule term for derivatives
c
                  de = -fik / (rirkr3*r2)
c
c     form the chain rule terms for first derivatives
c
                  deddotp = -de * r2
                  deddoti = de * 3.0d0*dotk
                  deddotk = de * 3.0d0*doti
                  dedr = de * (3.0d0*dotp-15.0d0*dotik/r2)
                  dedrirk = de * enum
c
c     more first derivative chain rule expressions
c
                  termx = dedr*xr + deddoti*xi + deddotk*xk
                  termy = dedr*yr + deddoti*yi + deddotk*yk
                  termz = dedr*zr + deddoti*zi + deddotk*zk
                  termxi = dedrirk*xiri2 + deddotp*xk + deddoti*xr
                  termyi = dedrirk*yiri2 + deddotp*yk + deddoti*yr
                  termzi = dedrirk*ziri2 + deddotp*zk + deddoti*zr
                  termxk = dedrirk*xkrk2 + deddotp*xi + deddotk*xr
                  termyk = dedrirk*ykrk2 + deddotp*yi + deddotk*yr
                  termzk = dedrirk*zkrk2 + deddotp*zi + deddotk*zr
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
                     dedxi1 = si1*termx - termxi
                     dedyi1 = si1*termy - termyi
                     dedzi1 = si1*termz - termzi
                     dedxi2 = si2*termx + termxi
                     dedyi2 = si2*termy + termyi
                     dedzi2 = si2*termz + termzi
                     dedxk1 = -sk1*termx - termxk
                     dedyk1 = -sk1*termy - termyk
                     dedzk1 = -sk1*termz - termzk
                     dedxk2 = -sk2*termx + termxk
                     dedyk2 = -sk2*termy + termyk
                     dedzk2 = -sk2*termz + termzk
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
                     de = de * taper
                     termx = termx * taper
                     termy = termy * taper
                     termz = termz * taper
                     termxi = termxi * taper
                     termyi = termyi * taper
                     termzi = termzi * taper
                     termxk = termxk * taper
                     termyk = termyk * taper
                     termzk = termzk * taper
                  end if
c
c     next, find the second derivative chain rule terms
c
                  if (i .eq. i1) then
                     dtdxi1 = -5.0d0*si1*xrr2 + xiri2
                     part = si1*xkdoti - dotk*xr + si1*xidotk
     &                         - 2.0d0*si1*dotik*xrr2
                     partik = -xk*r2 + 2.0d0*si1*dotp*xr
     &                           - 3.0d0*si1*xkdoti + 3.0d0*xr*dotk
     &                           - 3.0d0*si1*xidotk
                     factor = 3.0d0*si1*dotp - 6.0d0*xkxr
     &                           + 6.0d0*si1*xixk - 3.0d0*dotk
     &                           - r2inv*(xr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*xkxr + xiri2*partik
     &                            - enum*(ri2inv-2.0d0*xiri2*xiri2)
                     factork = r2 + 3.0d0*si1*doti + si1*xixr
     &                            - xrxr + xkrk2*partik
                     dtxdxi1 = dtdxi1*termx + de*factor
                     dtxidxi1 = dtdxi1*termxi + de*factori
                     dtxkdxi1 = dtdxi1*termxk + de*factork
                     factor = -3.0d0*xkyr - 3.0d0*ykxr + 3.0d0*si1*xiyk
     &                           + 3.0d0*si1*yixk - r2inv*yr*part
                     factori = -2.0d0*si1*ykxr + 3.0d0*si1*xkyr
     &                           + yiri2*partik + 2.0d0*enum*yiri2*xiri2
                     factork = -2.0d0*si1*yixr - xryr + 3.0d0*si1*xiyr
     &                            + ykrk2*partik
                     dtydxi1 = dtdxi1*termy + de*factor
                     dtyidxi1 = dtdxi1*termyi + de*factori
                     dtykdxi1 = dtdxi1*termyk + de*factork
                     factor = -3.0d0*xkzr - 3.0d0*zkxr + 3.0d0*si1*xizk
     &                           + 3.0d0*si1*zixk - r2inv*zr*part
                     factori = -2.0d0*si1*zkxr + 3.0d0*si1*xkzr
     &                           + ziri2*partik + 2.0d0*enum*ziri2*xiri2
                     factork = -2.0d0*si1*zixr - xrzr + 3.0d0*si1*xizr
     &                            + zkrk2*partik
                     dtzdxi1 = dtdxi1*termz + de*factor
                     dtzidxi1 = dtdxi1*termzi + de*factori
                     dtzkdxi1 = dtdxi1*termzk + de*factork
c
                     dtdyi1 = -5.0d0*si1*yrr2 + yiri2
                     part = si1*ykdoti - dotk*yr + si1*yidotk
     &                         - 2.0d0*si1*dotik*yrr2
                     partik = -yk*r2 + 2.0d0*si1*dotp*yr
     &                           - 3.0d0*si1*ykdoti + 3.0d0*yr*dotk
     &                           - 3.0d0*si1*yidotk
                     factor = -3.0d0*ykxr - 3.0d0*xkyr + 3.0d0*si1*yixk
     &                           + 3.0d0*si1*xiyk - r2inv*xr*part
                     factori = -2.0d0*si1*xkyr + 3.0d0*si1*ykxr
     &                           + xiri2*partik + 2.0d0*enum*xiri2*yiri2
                     factork = -2.0d0*si1*xiyr - xryr + 3.0d0*si1*yixr
     &                            + xkrk2*partik
                     dtxdyi1 = dtdyi1*termx + de*factor
                     dtxidyi1 = dtdyi1*termxi + de*factori
                     dtxkdyi1 = dtdyi1*termxk + de*factork
                     factor = 3.0d0*si1*dotp - 6.0d0*ykyr
     &                           + 6.0d0*si1*yiyk - 3.0d0*dotk
     &                           - r2inv*(yr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*ykyr + yiri2*partik
     &                            - enum*(ri2inv-2.0d0*yiri2*yiri2)
                     factork = r2 + 3.0d0*si1*doti + si1*yiyr
     &                            - yryr + ykrk2*partik
                     dtydyi1 = dtdyi1*termy + de*factor
                     dtyidyi1 = dtdyi1*termyi + de*factori
                     dtykdyi1 = dtdyi1*termyk + de*factork
                     factor = -3.0d0*ykzr - 3.0d0*zkyr + 3.0d0*si1*yizk
     &                           + 3.0d0*si1*ziyk - r2inv*zr*part
                     factori = -2.0d0*si1*zkyr + 3.0d0*si1*ykzr
     &                           + ziri2*partik + 2.0d0*enum*ziri2*yiri2
                     factork = -2.0d0*si1*ziyr - yrzr + 3.0d0*si1*yizr
     &                            + zkrk2*partik
                     dtzdyi1 = dtdyi1*termz + de*factor
                     dtzidyi1 = dtdyi1*termzi + de*factori
                     dtzkdyi1 = dtdyi1*termzk + de*factork
c
                     dtdzi1 = -5.0d0*si1*zrr2 + ziri2
                     part = si1*zkdoti - dotk*zr + si1*zidotk
     &                         - 2.0d0*si1*dotik*zrr2
                     partik = -zk*r2 + 2.0d0*si1*dotp*zr
     &                           - 3.0d0*si1*zkdoti + 3.0d0*zr*dotk
     &                           - 3.0d0*si1*zidotk
                     factor = -3.0d0*zkxr - 3.0d0*xkzr + 3.0d0*si1*zixk
     &                           + 3.0d0*si1*xizk - r2inv*xr*part
                     factori = -2.0d0*si1*xkzr + 3.0d0*si1*zkxr
     &                           + xiri2*partik + 2.0d0*enum*xiri2*ziri2
                     factork = -2.0d0*si1*xizr - xrzr + 3.0d0*si1*zixr
     &                            + xkrk2*partik
                     dtxdzi1 = dtdzi1*termx + de*factor
                     dtxidzi1 = dtdzi1*termxi + de*factori
                     dtxkdzi1 = dtdzi1*termxk + de*factork
                     factor = -3.0d0*zkyr - 3.0d0*ykzr + 3.0d0*si1*ziyk
     &                           + 3.0d0*si1*yizk - r2inv*yr*part
                     factori = -2.0d0*si1*ykzr + 3.0d0*si1*zkyr
     &                           + yiri2*partik + 2.0d0*enum*yiri2*ziri2
                     factork = -2.0d0*si1*yizr - yrzr + 3.0d0*si1*ziyr
     &                            + ykrk2*partik
                     dtydzi1 = dtdzi1*termy + de*factor
                     dtyidzi1 = dtdzi1*termyi + de*factori
                     dtykdzi1 = dtdzi1*termyk + de*factork
                     factor = 3.0d0*si1*dotp - 6.0d0*zkzr
     &                           + 6.0d0*si1*zizk - 3.0d0*dotk
     &                           - r2inv*(zr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*zkzr + ziri2*partik
     &                            - enum*(ri2inv-2.0d0*ziri2*ziri2)
                     factork = r2 + 3.0d0*si1*doti + si1*zizr
     &                            - zrzr + zkrk2*partik
                     dtzdzi1 = dtdzi1*termz + de*factor
                     dtzidzi1 = dtdzi1*termzi + de*factori
                     dtzkdzi1 = dtdzi1*termzk + de*factork
c
                  else if (i .eq. i2) then
                     dtdxi2 = -5.0d0*si2*xrr2 - xiri2
                     part = si2*xkdoti + dotk*xr + si2*xidotk
     &                         - 2.0d0*si2*dotik*xrr2
                     partik = xk*r2 + 2.0d0*si2*dotp*xr
     &                           - 3.0d0*si2*xkdoti - 3.0d0*xr*dotk
     &                           - 3.0d0*si2*xidotk
                     factor = 3.0d0*si2*dotp + 6.0d0*xkxr
     &                           + 6.0d0*si2*xixk + 3.0d0*dotk
     &                           - r2inv*(xr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*xkxr + xiri2*partik
     &                            + enum*(ri2inv-2.0d0*xiri2*xiri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*xixr
     &                            + xrxr + xkrk2*partik
                     dtxdxi2 = dtdxi2*termx + de*factor
                     dtxidxi2 = dtdxi2*termxi + de*factori
                     dtxkdxi2 = dtdxi2*termxk + de*factork
                     factor = 3.0d0*xkyr + 3.0d0*ykxr + 3.0d0*si2*xiyk
     &                           + 3.0d0*si2*yixk - r2inv*yr*part
                     factori = -2.0d0*si2*ykxr + 3.0d0*si2*xkyr
     &                           + yiri2*partik - 2.0d0*enum*yiri2*xiri2
                     factork = -2.0d0*si2*yixr + xryr + 3.0d0*si2*xiyr
     &                            + ykrk2*partik
                     dtydxi2 = dtdxi2*termy + de*factor
                     dtyidxi2 = dtdxi2*termyi + de*factori
                     dtykdxi2 = dtdxi2*termyk + de*factork
                     factor = 3.0d0*xkzr + 3.0d0*zkxr + 3.0d0*si2*xizk
     &                           + 3.0d0*si2*zixk - r2inv*zr*part
                     factori = -2.0d0*si2*zkxr + 3.0d0*si2*xkzr
     &                           + ziri2*partik - 2.0d0*enum*ziri2*xiri2
                     factork = -2.0d0*si2*zixr + xrzr + 3.0d0*si2*xizr
     &                            + zkrk2*partik
                     dtzdxi2 = dtdxi2*termz + de*factor
                     dtzidxi2 = dtdxi2*termzi + de*factori
                     dtzkdxi2 = dtdxi2*termzk + de*factork
c
                     dtdyi2 = -5.0d0*si2*yrr2 - yiri2
                     part = si2*ykdoti + dotk*yr + si2*yidotk
     &                         - 2.0d0*si2*dotik*yrr2
                     partik = yk*r2 + 2.0d0*si2*dotp*yr
     &                           - 3.0d0*si2*ykdoti - 3.0d0*yr*dotk
     &                           - 3.0d0*si2*yidotk
                     factor = 3.0d0*ykxr + 3.0d0*xkyr + 3.0d0*si2*yixk
     &                           + 3.0d0*si2*xiyk - r2inv*xr*part
                     factori = -2.0d0*si2*xkyr + 3.0d0*si2*ykxr
     &                           + xiri2*partik - 2.0d0*enum*xiri2*yiri2
                     factork = -2.0d0*si2*xiyr + xryr + 3.0d0*si2*yixr
     &                            + xkrk2*partik
                     dtxdyi2 = dtdyi2*termx + de*factor
                     dtxidyi2 = dtdyi2*termxi + de*factori
                     dtxkdyi2 = dtdyi2*termxk + de*factork
                     factor = 3.0d0*si2*dotp + 6.0d0*ykyr
     &                           + 6.0d0*si2*yiyk + 3.0d0*dotk
     &                           - r2inv*(yr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*ykyr + yiri2*partik
     &                            + enum*(ri2inv-2.0d0*yiri2*yiri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*yiyr
     &                            + yryr + ykrk2*partik
                     dtydyi2 = dtdyi2*termy + de*factor
                     dtyidyi2 = dtdyi2*termyi + de*factori
                     dtykdyi2 = dtdyi2*termyk + de*factork
                     factor = 3.0d0*ykzr + 3.0d0*zkyr + 3.0d0*si2*yizk
     &                           + 3.0d0*si2*ziyk - r2inv*zr*part
                     factori = -2.0d0*si2*zkyr + 3.0d0*si2*ykzr
     &                           + ziri2*partik - 2.0d0*enum*ziri2*yiri2
                     factork = -2.0d0*si2*ziyr + yrzr + 3.0d0*si2*yizr
     &                            + zkrk2*partik
                     dtzdyi2 = dtdyi2*termz + de*factor
                     dtzidyi2 = dtdyi2*termzi + de*factori
                     dtzkdyi2 = dtdyi2*termzk + de*factork
c
                     dtdzi2 = -5.0d0*si2*zrr2 - ziri2
                     part = si2*zkdoti + dotk*zr + si2*zidotk
     &                         - 2.0d0*si2*dotik*zrr2
                     partik = zk*r2 + 2.0d0*si2*dotp*zr
     &                           - 3.0d0*si2*zkdoti - 3.0d0*zr*dotk
     &                           - 3.0d0*si2*zidotk
                     factor = 3.0d0*zkxr + 3.0d0*xkzr + 3.0d0*si2*zixk
     &                           + 3.0d0*si2*xizk - r2inv*xr*part
                     factori = -2.0d0*si2*xkzr + 3.0d0*si2*zkxr
     &                           + xiri2*partik - 2.0d0*enum*xiri2*ziri2
                     factork = -2.0d0*si2*xizr + xrzr + 3.0d0*si2*zixr
     &                            + xkrk2*partik
                     dtxdzi2 = dtdzi2*termx + de*factor
                     dtxidzi2 = dtdzi2*termxi + de*factori
                     dtxkdzi2 = dtdzi2*termxk + de*factork
                     factor = 3.0d0*zkyr + 3.0d0*ykzr + 3.0d0*si2*ziyk
     &                           + 3.0d0*si2*yizk - r2inv*yr*part
                     factori = -2.0d0*si2*ykzr + 3.0d0*si2*zkyr
     &                           + yiri2*partik - 2.0d0*enum*yiri2*ziri2
                     factork = -2.0d0*si2*yizr + yrzr + 3.0d0*si2*ziyr
     &                            + ykrk2*partik
                     dtydzi2 = dtdzi2*termy + de*factor
                     dtyidzi2 = dtdzi2*termyi + de*factori
                     dtykdzi2 = dtdzi2*termyk + de*factork
                     factor = 3.0d0*si2*dotp + 6.0d0*zkzr
     &                           + 6.0d0*si2*zizk + 3.0d0*dotk
     &                           - r2inv*(zr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*zkzr + ziri2*partik
     &                            + enum*(ri2inv-2.0d0*ziri2*ziri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*zizr
     &                            + zrzr + zkrk2*partik
                     dtzdzi2 = dtdzi2*termz + de*factor
                     dtzidzi2 = dtdzi2*termzi + de*factori
                     dtzkdzi2 = dtdzi2*termzk + de*factork
                  end if
c
c     now, increment diagonal and off-diagonal Hessian elements
c
                  if (i .eq. i1) then
                     hessx(1,i1) = hessx(1,i1) + si1*dtxdxi1 - dtxidxi1
                     hessx(2,i1) = hessx(2,i1) + si1*dtydxi1 - dtyidxi1
                     hessx(3,i1) = hessx(3,i1) + si1*dtzdxi1 - dtzidxi1
                     hessx(1,i2) = hessx(1,i2) + si2*dtxdxi1 + dtxidxi1
                     hessx(2,i2) = hessx(2,i2) + si2*dtydxi1 + dtyidxi1
                     hessx(3,i2) = hessx(3,i2) + si2*dtzdxi1 + dtzidxi1
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi1 - dtxkdxi1
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi1 - dtykdxi1
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi1 - dtzkdxi1
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi1 + dtxkdxi1
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi1 + dtykdxi1
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi1 + dtzkdxi1
                     hessy(1,i1) = hessy(1,i1) + si1*dtxdyi1 - dtxidyi1
                     hessy(2,i1) = hessy(2,i1) + si1*dtydyi1 - dtyidyi1
                     hessy(3,i1) = hessy(3,i1) + si1*dtzdyi1 - dtzidyi1
                     hessy(1,i2) = hessy(1,i2) + si2*dtxdyi1 + dtxidyi1
                     hessy(3,i2) = hessy(3,i2) + si2*dtzdyi1 + dtzidyi1
                     hessy(2,i2) = hessy(2,i2) + si2*dtydyi1 + dtyidyi1
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi1 - dtxkdyi1
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi1 - dtykdyi1
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi1 - dtzkdyi1
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi1 + dtxkdyi1
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi1 + dtykdyi1
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi1 + dtzkdyi1
                     hessz(1,i1) = hessz(1,i1) + si1*dtxdzi1 - dtxidzi1
                     hessz(2,i1) = hessz(2,i1) + si1*dtydzi1 - dtyidzi1
                     hessz(3,i1) = hessz(3,i1) + si1*dtzdzi1 - dtzidzi1
                     hessz(1,i2) = hessz(1,i2) + si2*dtxdzi1 + dtxidzi1
                     hessz(2,i2) = hessz(2,i2) + si2*dtydzi1 + dtyidzi1
                     hessz(3,i2) = hessz(3,i2) + si2*dtzdzi1 + dtzidzi1
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi1 - dtxkdzi1
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi1 - dtykdzi1
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi1 - dtzkdzi1
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi1 + dtxkdzi1
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi1 + dtykdzi1
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi1 + dtzkdzi1
                  else if (i .eq. i2) then
                     hessx(1,i1) = hessx(1,i1) + si1*dtxdxi2 - dtxidxi2
                     hessx(2,i1) = hessx(2,i1) + si1*dtydxi2 - dtyidxi2
                     hessx(3,i1) = hessx(3,i1) + si1*dtzdxi2 - dtzidxi2
                     hessx(1,i2) = hessx(1,i2) + si2*dtxdxi2 + dtxidxi2
                     hessx(2,i2) = hessx(2,i2) + si2*dtydxi2 + dtyidxi2
                     hessx(3,i2) = hessx(3,i2) + si2*dtzdxi2 + dtzidxi2
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi2 - dtxkdxi2
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi2 - dtykdxi2
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi2 - dtzkdxi2
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi2 + dtxkdxi2
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi2 + dtykdxi2
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi2 + dtzkdxi2
                     hessy(1,i1) = hessy(1,i1) + si1*dtxdyi2 - dtxidyi2
                     hessy(2,i1) = hessy(2,i1) + si1*dtydyi2 - dtyidyi2
                     hessy(3,i1) = hessy(3,i1) + si1*dtzdyi2 - dtzidyi2
                     hessy(1,i2) = hessy(1,i2) + si2*dtxdyi2 + dtxidyi2
                     hessy(2,i2) = hessy(2,i2) + si2*dtydyi2 + dtyidyi2
                     hessy(3,i2) = hessy(3,i2) + si2*dtzdyi2 + dtzidyi2
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi2 - dtxkdyi2
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi2 - dtykdyi2
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi2 - dtzkdyi2
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi2 + dtxkdyi2
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi2 + dtykdyi2
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi2 + dtzkdyi2
                     hessz(1,i1) = hessz(1,i1) + si1*dtxdzi2 - dtxidzi2
                     hessz(2,i1) = hessz(2,i1) + si1*dtydzi2 - dtyidzi2
                     hessz(3,i1) = hessz(3,i1) + si1*dtzdzi2 - dtzidzi2
                     hessz(1,i2) = hessz(1,i2) + si2*dtxdzi2 + dtxidzi2
                     hessz(2,i2) = hessz(2,i2) + si2*dtydzi2 + dtyidzi2
                     hessz(3,i2) = hessz(3,i2) + si2*dtzdzi2 + dtzidzi2
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi2 - dtxkdzi2
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi2 - dtykdzi2
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi2 - dtzkdzi2
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi2 + dtxkdzi2
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi2 + dtykdzi2
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi2 + dtzkdzi2
                  end if
c
c     more energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     if (i .eq. i1) then
                        hessx(1,i1) = hessx(1,i1) + si1*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxi1 + si1*si1*d2taperxx
                        hessx(2,i1) = hessx(2,i1) + si1*dtapery*dedxi1
     &                         + si1*dtaperx*dedyi1 + si1*si1*d2taperxy
                        hessx(3,i1) = hessx(3,i1) + si1*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzi1 + si1*si1*d2taperxz
                        hessx(1,i2) = hessx(1,i2) + si2*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxi2 + si2*si1*d2taperxx
                        hessx(2,i2) = hessx(2,i2) + si2*dtapery*dedxi1
     &                         + si1*dtaperx*dedyi2 + si2*si1*d2taperxy
                        hessx(3,i2) = hessx(3,i2) + si2*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzi2 + si2*si1*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxk1 - sk1*si1*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtapery*dedxi1
     &                         + si1*dtaperx*dedyk1 - sk1*si1*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzk1 - sk1*si1*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxk2 - sk2*si1*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtapery*dedxi1
     &                         + si1*dtaperx*dedyk2 - sk2*si1*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzk2 - sk2*si1*d2taperxz
                        hessy(1,i1) = hessy(1,i1) + si1*dtaperx*dedyi1
     &                         + si1*dtapery*dedxi1 + si1*si1*d2taperxy
                        hessy(2,i1) = hessy(2,i1) + si1*dtapery*dedyi1
     &                         + si1*dtapery*dedyi1 + si1*si1*d2taperyy
                        hessy(3,i1) = hessy(3,i1) + si1*dtaperz*dedyi1
     &                         + si1*dtapery*dedzi1 + si1*si1*d2taperyz
                        hessy(1,i2) = hessy(1,i2) + si2*dtaperx*dedyi1
     &                         + si1*dtapery*dedxi2 + si2*si1*d2taperxy
                        hessy(2,i2) = hessy(2,i2) + si2*dtapery*dedyi1
     &                         + si1*dtapery*dedyi2 + si2*si1*d2taperyy
                        hessy(3,i2) = hessy(3,i2) + si2*dtaperz*dedyi1
     &                         + si1*dtapery*dedzi2 + si2*si1*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtaperx*dedyi1
     &                         + si1*dtapery*dedxk1 - sk1*si1*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyi1
     &                         + si1*dtapery*dedyk1 - sk1*si1*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtaperz*dedyi1
     &                         + si1*dtapery*dedzk1 - sk1*si1*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtaperx*dedyi1
     &                         + si1*dtapery*dedxk2 - sk2*si1*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyi1
     &                         + si1*dtapery*dedyk2 - sk2*si1*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtaperz*dedyi1
     &                         + si1*dtapery*dedzk2 - sk2*si1*d2taperyz
                        hessz(1,i1) = hessz(1,i1) + si1*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxi1 + si1*si1*d2taperxz
                        hessz(2,i1) = hessz(2,i1) + si1*dtapery*dedzi1
     &                         + si1*dtaperz*dedyi1 + si1*si1*d2taperyz
                        hessz(3,i1) = hessz(3,i1) + si1*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzi1 + si1*si1*d2taperzz
                        hessz(1,i2) = hessz(1,i2) + si2*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxi2 + si2*si1*d2taperxz
                        hessz(2,i2) = hessz(2,i2) + si2*dtapery*dedzi1
     &                         + si1*dtaperz*dedyi2 + si2*si1*d2taperyz
                        hessz(3,i2) = hessz(3,i2) + si2*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzi2 + si2*si1*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxk1 - sk1*si1*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtapery*dedzi1
     &                         + si1*dtaperz*dedyk1 - sk1*si1*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzk1 - sk1*si1*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxk2 - sk2*si1*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtapery*dedzi1
     &                         + si1*dtaperz*dedyk2 - sk2*si1*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzk2 - sk2*si1*d2taperzz
                     else if (i .eq. i2) then
                        hessx(1,i1) = hessx(1,i1) + si1*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxi1 + si1*si2*d2taperxx
                        hessx(2,i1) = hessx(2,i1) + si1*dtapery*dedxi2
     &                         + si2*dtaperx*dedyi1 + si1*si2*d2taperxy
                        hessx(3,i1) = hessx(3,i1) + si1*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzi1 + si1*si2*d2taperxz
                        hessx(1,i2) = hessx(1,i2) + si2*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxi2 + si2*si2*d2taperxx
                        hessx(2,i2) = hessx(2,i2) + si2*dtapery*dedxi2
     &                         + si2*dtaperx*dedyi2 + si2*si2*d2taperxy
                        hessx(3,i2) = hessx(3,i2) + si2*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzi2 + si2*si2*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxk1 - sk1*si2*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtapery*dedxi2
     &                         + si2*dtaperx*dedyk1 - sk1*si2*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzk1 - sk1*si2*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxk2 - sk2*si2*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtapery*dedxi2
     &                         + si2*dtaperx*dedyk2 - sk2*si2*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzk2 - sk2*si2*d2taperxz
                        hessy(1,i1) = hessy(1,i1) + si1*dtaperx*dedyi2
     &                         + si2*dtapery*dedxi1 + si1*si2*d2taperxy
                        hessy(2,i1) = hessy(2,i1) + si1*dtapery*dedyi2
     &                         + si2*dtapery*dedyi1 + si1*si2*d2taperyy
                        hessy(3,i1) = hessy(3,i1) + si1*dtaperz*dedyi2
     &                         + si2*dtapery*dedzi1 + si1*si2*d2taperyz
                        hessy(1,i2) = hessy(1,i2) + si2*dtaperx*dedyi2
     &                         + si2*dtapery*dedxi2 + si2*si2*d2taperxy
                        hessy(2,i2) = hessy(2,i2) + si2*dtapery*dedyi2
     &                         + si2*dtapery*dedyi2 + si2*si2*d2taperyy
                        hessy(3,i2) = hessy(3,i2) + si2*dtaperz*dedyi2
     &                         + si2*dtapery*dedzi2 + si2*si2*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtaperx*dedyi2
     &                         + si2*dtapery*dedxk1 - sk1*si2*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyi2
     &                         + si2*dtapery*dedyk1 - sk1*si2*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtaperz*dedyi2
     &                         + si2*dtapery*dedzk1 - sk1*si2*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtaperx*dedyi2
     &                         + si2*dtapery*dedxk2 - sk2*si2*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyi2
     &                         + si2*dtapery*dedyk2 - sk2*si2*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtaperz*dedyi2
     &                         + si2*dtapery*dedzk2 - sk2*si2*d2taperyz
                        hessz(1,i1) = hessz(1,i1) + si1*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxi1 + si1*si2*d2taperxz
                        hessz(2,i1) = hessz(2,i1) + si1*dtapery*dedzi2
     &                         + si2*dtaperz*dedyi1 + si1*si2*d2taperyz
                        hessz(3,i1) = hessz(3,i1) + si1*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzi1 + si1*si2*d2taperzz
                        hessz(1,i2) = hessz(1,i2) + si2*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxi2 + si2*si2*d2taperxz
                        hessz(2,i2) = hessz(2,i2) + si2*dtapery*dedzi2
     &                         + si2*dtaperz*dedyi2 + si2*si2*d2taperyz
                        hessz(3,i2) = hessz(3,i2) + si2*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzi2 + si2*si2*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxk1 - sk1*si2*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtapery*dedzi2
     &                         + si2*dtaperz*dedyk1 - sk1*si2*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzk1 - sk1*si2*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxk2 - sk2*si2*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtapery*dedzi2
     &                         + si2*dtaperz*dedyk2 - sk2*si2*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzk2 - sk2*si2*d2taperzz
                     end if
                  end if
               end if
            end if
         end do
   10    continue
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do idipole = 1, ndipole
         i1 = idpl(1,idipole)
         i2 = idpl(2,idipole)
         si1 = 1.0d0 - sdpl(idipole)
         si2 = sdpl(idipole)
         if (i1.ne.i .and. i2.ne.i)  goto 30
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*si2
         yq = y(i1) + yi*si2
         zq = z(i1) + zi*si2
         fi = f * bdpl(idipole)
         do kdipole = 1, ndipole
            k1 = idpl(1,kdipole)
            k2 = idpl(2,kdipole)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,4,i1,i2,k1,k2,0)
c
c     compute the energy contribution for this interaction
c
            if (.not. proceed)  goto 20
            sk1 = 1.0d0 - sdpl(kdipole)
            sk2 = sdpl(kdipole)
            do jcell = 1, ncell
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = xq - x(k1) - xk*sk2
               yr = yq - y(k1) - yk*sk2
               zr = zq - z(k1) - zk*sk2
               call image (xr,yr,zr,jcell)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr3 = sqrt(ri2*rk2*r2) * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(kdipole)
c
c     some abbreviations used in various chain rule terms
c
                  dotik = doti * dotk
                  enum = dotp*r2 - 3.0d0*dotik
                  r2inv = 15.0d0 / r2
                  ri2inv = 1.0d0 / ri2
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xiri2 = xi / ri2
                  yiri2 = yi / ri2
                  ziri2 = zi / ri2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  xixr = xi * xr
                  xiyr = xi * yr
                  xizr = xi * zr
                  yixr = yi * xr
                  yiyr = yi * yr
                  yizr = yi * zr
                  zixr = zi * xr
                  ziyr = zi * yr
                  zizr = zi * zr
                  xkxr = xk * xr
                  xkyr = xk * yr
                  xkzr = xk * zr
                  ykxr = yk * xr
                  ykyr = yk * yr
                  ykzr = yk * zr
                  zkxr = zk * xr
                  zkyr = zk * yr
                  zkzr = zk * zr
                  xixk = xi * xk
                  xiyk = xi * yk
                  xizk = xi * zk
                  yixk = yi * xk
                  yiyk = yi * yk
                  yizk = yi * zk
                  zixk = zi * xk
                  ziyk = zi * yk
                  zizk = zi * zk
                  xrxr = 3.0d0 * xr * xr
                  xryr = 3.0d0 * xr * yr
                  xrzr = 3.0d0 * xr * zr
                  yryr = 3.0d0 * yr * yr
                  yrzr = 3.0d0 * yr * zr
                  zrzr = 3.0d0 * zr * zr
                  xidotk = xi * dotk
                  yidotk = yi * dotk
                  zidotk = zi * dotk
                  xkdoti = xk * doti
                  ykdoti = yk * doti
                  zkdoti = zk * doti
c
c     scale the interaction based on its group membership
c
                  if (use_group)  fik = fik * fgrp
c
c     form the master chain rule term for derivatives
c
                  de = -fik / (rirkr3*r2)
c
c     form the chain rule terms for first derivatives
c
                  deddotp = -de * r2
                  deddoti = de * 3.0d0*dotk
                  deddotk = de * 3.0d0*doti
                  dedr = de * (3.0d0*dotp-15.0d0*dotik/r2)
                  dedrirk = de * enum
c
c     more first derivative chain rule expressions
c
                  termx = dedr*xr + deddoti*xi + deddotk*xk
                  termy = dedr*yr + deddoti*yi + deddotk*yk
                  termz = dedr*zr + deddoti*zi + deddotk*zk
                  termxi = dedrirk*xiri2 + deddotp*xk + deddoti*xr
                  termyi = dedrirk*yiri2 + deddotp*yk + deddoti*yr
                  termzi = dedrirk*ziri2 + deddotp*zk + deddoti*zr
                  termxk = dedrirk*xkrk2 + deddotp*xi + deddotk*xr
                  termyk = dedrirk*ykrk2 + deddotp*yi + deddotk*yr
                  termzk = dedrirk*zkrk2 + deddotp*zi + deddotk*zr
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
                     dedxi1 = si1*termx - termxi
                     dedyi1 = si1*termy - termyi
                     dedzi1 = si1*termz - termzi
                     dedxi2 = si2*termx + termxi
                     dedyi2 = si2*termy + termyi
                     dedzi2 = si2*termz + termzi
                     dedxk1 = -sk1*termx - termxk
                     dedyk1 = -sk1*termy - termyk
                     dedzk1 = -sk1*termz - termzk
                     dedxk2 = -sk2*termx + termxk
                     dedyk2 = -sk2*termy + termyk
                     dedzk2 = -sk2*termz + termzk
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
                     de = de * taper
                     termx = termx * taper
                     termy = termy * taper
                     termz = termz * taper
                     termxi = termxi * taper
                     termyi = termyi * taper
                     termzi = termzi * taper
                     termxk = termxk * taper
                     termyk = termyk * taper
                     termzk = termzk * taper
                  end if
c
c     next, find the second derivative chain rule terms
c
                  if (i .eq. i1) then
                     dtdxi1 = -5.0d0*si1*xrr2 + xiri2
                     part = si1*xkdoti - dotk*xr + si1*xidotk
     &                         - 2.0d0*si1*dotik*xrr2
                     partik = -xk*r2 + 2.0d0*si1*dotp*xr
     &                           - 3.0d0*si1*xkdoti + 3.0d0*xr*dotk
     &                           - 3.0d0*si1*xidotk
                     factor = 3.0d0*si1*dotp - 6.0d0*xkxr
     &                           + 6.0d0*si1*xixk - 3.0d0*dotk
     &                           - r2inv*(xr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*xkxr + xiri2*partik
     &                            - enum*(ri2inv-2.0d0*xiri2*xiri2)
                     factork = r2 + 3.0d0*si1*doti + si1*xixr
     &                            - xrxr + xkrk2*partik
                     dtxdxi1 = dtdxi1*termx + de*factor
                     dtxidxi1 = dtdxi1*termxi + de*factori
                     dtxkdxi1 = dtdxi1*termxk + de*factork
                     factor = -3.0d0*xkyr - 3.0d0*ykxr + 3.0d0*si1*xiyk
     &                           + 3.0d0*si1*yixk - r2inv*yr*part
                     factori = -2.0d0*si1*ykxr + 3.0d0*si1*xkyr
     &                           + yiri2*partik + 2.0d0*enum*yiri2*xiri2
                     factork = -2.0d0*si1*yixr - xryr + 3.0d0*si1*xiyr
     &                            + ykrk2*partik
                     dtydxi1 = dtdxi1*termy + de*factor
                     dtyidxi1 = dtdxi1*termyi + de*factori
                     dtykdxi1 = dtdxi1*termyk + de*factork
                     factor = -3.0d0*xkzr - 3.0d0*zkxr + 3.0d0*si1*xizk
     &                           + 3.0d0*si1*zixk - r2inv*zr*part
                     factori = -2.0d0*si1*zkxr + 3.0d0*si1*xkzr
     &                           + ziri2*partik + 2.0d0*enum*ziri2*xiri2
                     factork = -2.0d0*si1*zixr - xrzr + 3.0d0*si1*xizr
     &                            + zkrk2*partik
                     dtzdxi1 = dtdxi1*termz + de*factor
                     dtzidxi1 = dtdxi1*termzi + de*factori
                     dtzkdxi1 = dtdxi1*termzk + de*factork
c
                     dtdyi1 = -5.0d0*si1*yrr2 + yiri2
                     part = si1*ykdoti - dotk*yr + si1*yidotk
     &                         - 2.0d0*si1*dotik*yrr2
                     partik = -yk*r2 + 2.0d0*si1*dotp*yr
     &                           - 3.0d0*si1*ykdoti + 3.0d0*yr*dotk
     &                           - 3.0d0*si1*yidotk
                     factor = -3.0d0*ykxr - 3.0d0*xkyr + 3.0d0*si1*yixk
     &                           + 3.0d0*si1*xiyk - r2inv*xr*part
                     factori = -2.0d0*si1*xkyr + 3.0d0*si1*ykxr
     &                           + xiri2*partik + 2.0d0*enum*xiri2*yiri2
                     factork = -2.0d0*si1*xiyr - xryr + 3.0d0*si1*yixr
     &                            + xkrk2*partik
                     dtxdyi1 = dtdyi1*termx + de*factor
                     dtxidyi1 = dtdyi1*termxi + de*factori
                     dtxkdyi1 = dtdyi1*termxk + de*factork
                     factor = 3.0d0*si1*dotp - 6.0d0*ykyr
     &                           + 6.0d0*si1*yiyk - 3.0d0*dotk
     &                           - r2inv*(yr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*ykyr + yiri2*partik
     &                            - enum*(ri2inv-2.0d0*yiri2*yiri2)
                     factork = r2 + 3.0d0*si1*doti + si1*yiyr
     &                            - yryr + ykrk2*partik
                     dtydyi1 = dtdyi1*termy + de*factor
                     dtyidyi1 = dtdyi1*termyi + de*factori
                     dtykdyi1 = dtdyi1*termyk + de*factork
                     factor = -3.0d0*ykzr - 3.0d0*zkyr + 3.0d0*si1*yizk
     &                           + 3.0d0*si1*ziyk - r2inv*zr*part
                     factori = -2.0d0*si1*zkyr + 3.0d0*si1*ykzr
     &                           + ziri2*partik + 2.0d0*enum*ziri2*yiri2
                     factork = -2.0d0*si1*ziyr - yrzr + 3.0d0*si1*yizr
     &                            + zkrk2*partik
                     dtzdyi1 = dtdyi1*termz + de*factor
                     dtzidyi1 = dtdyi1*termzi + de*factori
                     dtzkdyi1 = dtdyi1*termzk + de*factork
c
                     dtdzi1 = -5.0d0*si1*zrr2 + ziri2
                     part = si1*zkdoti - dotk*zr + si1*zidotk
     &                         - 2.0d0*si1*dotik*zrr2
                     partik = -zk*r2 + 2.0d0*si1*dotp*zr
     &                           - 3.0d0*si1*zkdoti + 3.0d0*zr*dotk
     &                           - 3.0d0*si1*zidotk
                     factor = -3.0d0*zkxr - 3.0d0*xkzr + 3.0d0*si1*zixk
     &                           + 3.0d0*si1*xizk - r2inv*xr*part
                     factori = -2.0d0*si1*xkzr + 3.0d0*si1*zkxr
     &                           + xiri2*partik + 2.0d0*enum*xiri2*ziri2
                     factork = -2.0d0*si1*xizr - xrzr + 3.0d0*si1*zixr
     &                            + xkrk2*partik
                     dtxdzi1 = dtdzi1*termx + de*factor
                     dtxidzi1 = dtdzi1*termxi + de*factori
                     dtxkdzi1 = dtdzi1*termxk + de*factork
                     factor = -3.0d0*zkyr - 3.0d0*ykzr + 3.0d0*si1*ziyk
     &                           + 3.0d0*si1*yizk - r2inv*yr*part
                     factori = -2.0d0*si1*ykzr + 3.0d0*si1*zkyr
     &                           + yiri2*partik + 2.0d0*enum*yiri2*ziri2
                     factork = -2.0d0*si1*yizr - yrzr + 3.0d0*si1*ziyr
     &                            + ykrk2*partik
                     dtydzi1 = dtdzi1*termy + de*factor
                     dtyidzi1 = dtdzi1*termyi + de*factori
                     dtykdzi1 = dtdzi1*termyk + de*factork
                     factor = 3.0d0*si1*dotp - 6.0d0*zkzr
     &                           + 6.0d0*si1*zizk - 3.0d0*dotk
     &                           - r2inv*(zr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*zkzr + ziri2*partik
     &                            - enum*(ri2inv-2.0d0*ziri2*ziri2)
                     factork = r2 + 3.0d0*si1*doti + si1*zizr
     &                            - zrzr + zkrk2*partik
                     dtzdzi1 = dtdzi1*termz + de*factor
                     dtzidzi1 = dtdzi1*termzi + de*factori
                     dtzkdzi1 = dtdzi1*termzk + de*factork
c
                  else if (i .eq. i2) then
                     dtdxi2 = -5.0d0*si2*xrr2 - xiri2
                     part = si2*xkdoti + dotk*xr + si2*xidotk
     &                         - 2.0d0*si2*dotik*xrr2
                     partik = xk*r2 + 2.0d0*si2*dotp*xr
     &                           - 3.0d0*si2*xkdoti - 3.0d0*xr*dotk
     &                           - 3.0d0*si2*xidotk
                     factor = 3.0d0*si2*dotp + 6.0d0*xkxr
     &                           + 6.0d0*si2*xixk + 3.0d0*dotk
     &                           - r2inv*(xr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*xkxr + xiri2*partik
     &                            + enum*(ri2inv-2.0d0*xiri2*xiri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*xixr
     &                            + xrxr + xkrk2*partik
                     dtxdxi2 = dtdxi2*termx + de*factor
                     dtxidxi2 = dtdxi2*termxi + de*factori
                     dtxkdxi2 = dtdxi2*termxk + de*factork
                     factor = 3.0d0*xkyr + 3.0d0*ykxr + 3.0d0*si2*xiyk
     &                           + 3.0d0*si2*yixk - r2inv*yr*part
                     factori = -2.0d0*si2*ykxr + 3.0d0*si2*xkyr
     &                           + yiri2*partik - 2.0d0*enum*yiri2*xiri2
                     factork = -2.0d0*si2*yixr + xryr + 3.0d0*si2*xiyr
     &                            + ykrk2*partik
                     dtydxi2 = dtdxi2*termy + de*factor
                     dtyidxi2 = dtdxi2*termyi + de*factori
                     dtykdxi2 = dtdxi2*termyk + de*factork
                     factor = 3.0d0*xkzr + 3.0d0*zkxr + 3.0d0*si2*xizk
     &                           + 3.0d0*si2*zixk - r2inv*zr*part
                     factori = -2.0d0*si2*zkxr + 3.0d0*si2*xkzr
     &                           + ziri2*partik - 2.0d0*enum*ziri2*xiri2
                     factork = -2.0d0*si2*zixr + xrzr + 3.0d0*si2*xizr
     &                            + zkrk2*partik
                     dtzdxi2 = dtdxi2*termz + de*factor
                     dtzidxi2 = dtdxi2*termzi + de*factori
                     dtzkdxi2 = dtdxi2*termzk + de*factork
c
                     dtdyi2 = -5.0d0*si2*yrr2 - yiri2
                     part = si2*ykdoti + dotk*yr + si2*yidotk
     &                         - 2.0d0*si2*dotik*yrr2
                     partik = yk*r2 + 2.0d0*si2*dotp*yr
     &                           - 3.0d0*si2*ykdoti - 3.0d0*yr*dotk
     &                           - 3.0d0*si2*yidotk
                     factor = 3.0d0*ykxr + 3.0d0*xkyr + 3.0d0*si2*yixk
     &                           + 3.0d0*si2*xiyk - r2inv*xr*part
                     factori = -2.0d0*si2*xkyr + 3.0d0*si2*ykxr
     &                           + xiri2*partik - 2.0d0*enum*xiri2*yiri2
                     factork = -2.0d0*si2*xiyr + xryr + 3.0d0*si2*yixr
     &                            + xkrk2*partik
                     dtxdyi2 = dtdyi2*termx + de*factor
                     dtxidyi2 = dtdyi2*termxi + de*factori
                     dtxkdyi2 = dtdyi2*termxk + de*factork
                     factor = 3.0d0*si2*dotp + 6.0d0*ykyr
     &                           + 6.0d0*si2*yiyk + 3.0d0*dotk
     &                           - r2inv*(yr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*ykyr + yiri2*partik
     &                            + enum*(ri2inv-2.0d0*yiri2*yiri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*yiyr
     &                            + yryr + ykrk2*partik
                     dtydyi2 = dtdyi2*termy + de*factor
                     dtyidyi2 = dtdyi2*termyi + de*factori
                     dtykdyi2 = dtdyi2*termyk + de*factork
                     factor = 3.0d0*ykzr + 3.0d0*zkyr + 3.0d0*si2*yizk
     &                           + 3.0d0*si2*ziyk - r2inv*zr*part
                     factori = -2.0d0*si2*zkyr + 3.0d0*si2*ykzr
     &                           + ziri2*partik - 2.0d0*enum*ziri2*yiri2
                     factork = -2.0d0*si2*ziyr + yrzr + 3.0d0*si2*yizr
     &                            + zkrk2*partik
                     dtzdyi2 = dtdyi2*termz + de*factor
                     dtzidyi2 = dtdyi2*termzi + de*factori
                     dtzkdyi2 = dtdyi2*termzk + de*factork
c
                     dtdzi2 = -5.0d0*si2*zrr2 - ziri2
                     part = si2*zkdoti + dotk*zr + si2*zidotk
     &                         - 2.0d0*si2*dotik*zrr2
                     partik = zk*r2 + 2.0d0*si2*dotp*zr
     &                           - 3.0d0*si2*zkdoti - 3.0d0*zr*dotk
     &                           - 3.0d0*si2*zidotk
                     factor = 3.0d0*zkxr + 3.0d0*xkzr + 3.0d0*si2*zixk
     &                           + 3.0d0*si2*xizk - r2inv*xr*part
                     factori = -2.0d0*si2*xkzr + 3.0d0*si2*zkxr
     &                           + xiri2*partik - 2.0d0*enum*xiri2*ziri2
                     factork = -2.0d0*si2*xizr + xrzr + 3.0d0*si2*zixr
     &                            + xkrk2*partik
                     dtxdzi2 = dtdzi2*termx + de*factor
                     dtxidzi2 = dtdzi2*termxi + de*factori
                     dtxkdzi2 = dtdzi2*termxk + de*factork
                     factor = 3.0d0*zkyr + 3.0d0*ykzr + 3.0d0*si2*ziyk
     &                           + 3.0d0*si2*yizk - r2inv*yr*part
                     factori = -2.0d0*si2*ykzr + 3.0d0*si2*zkyr
     &                           + yiri2*partik - 2.0d0*enum*yiri2*ziri2
                     factork = -2.0d0*si2*yizr + yrzr + 3.0d0*si2*ziyr
     &                            + ykrk2*partik
                     dtydzi2 = dtdzi2*termy + de*factor
                     dtyidzi2 = dtdzi2*termyi + de*factori
                     dtykdzi2 = dtdzi2*termyk + de*factork
                     factor = 3.0d0*si2*dotp + 6.0d0*zkzr
     &                           + 6.0d0*si2*zizk + 3.0d0*dotk
     &                           - r2inv*(zr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*zkzr + ziri2*partik
     &                            + enum*(ri2inv-2.0d0*ziri2*ziri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*zizr
     &                            + zrzr + zkrk2*partik
                     dtzdzi2 = dtdzi2*termz + de*factor
                     dtzidzi2 = dtdzi2*termzi + de*factori
                     dtzkdzi2 = dtdzi2*termzk + de*factork
                  end if
c
c     now, increment diagonal and off-diagonal Hessian elements
c
                  if (i .eq. i1) then
                     hessx(1,i1) = hessx(1,i1) + si1*dtxdxi1 - dtxidxi1
                     hessx(2,i1) = hessx(2,i1) + si1*dtydxi1 - dtyidxi1
                     hessx(3,i1) = hessx(3,i1) + si1*dtzdxi1 - dtzidxi1
                     hessx(1,i2) = hessx(1,i2) + si2*dtxdxi1 + dtxidxi1
                     hessx(2,i2) = hessx(2,i2) + si2*dtydxi1 + dtyidxi1
                     hessx(3,i2) = hessx(3,i2) + si2*dtzdxi1 + dtzidxi1
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi1 - dtxkdxi1
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi1 - dtykdxi1
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi1 - dtzkdxi1
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi1 + dtxkdxi1
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi1 + dtykdxi1
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi1 + dtzkdxi1
                     hessy(1,i1) = hessy(1,i1) + si1*dtxdyi1 - dtxidyi1
                     hessy(2,i1) = hessy(2,i1) + si1*dtydyi1 - dtyidyi1
                     hessy(3,i1) = hessy(3,i1) + si1*dtzdyi1 - dtzidyi1
                     hessy(1,i2) = hessy(1,i2) + si2*dtxdyi1 + dtxidyi1
                     hessy(3,i2) = hessy(3,i2) + si2*dtzdyi1 + dtzidyi1
                     hessy(2,i2) = hessy(2,i2) + si2*dtydyi1 + dtyidyi1
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi1 - dtxkdyi1
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi1 - dtykdyi1
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi1 - dtzkdyi1
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi1 + dtxkdyi1
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi1 + dtykdyi1
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi1 + dtzkdyi1
                     hessz(1,i1) = hessz(1,i1) + si1*dtxdzi1 - dtxidzi1
                     hessz(2,i1) = hessz(2,i1) + si1*dtydzi1 - dtyidzi1
                     hessz(3,i1) = hessz(3,i1) + si1*dtzdzi1 - dtzidzi1
                     hessz(1,i2) = hessz(1,i2) + si2*dtxdzi1 + dtxidzi1
                     hessz(2,i2) = hessz(2,i2) + si2*dtydzi1 + dtyidzi1
                     hessz(3,i2) = hessz(3,i2) + si2*dtzdzi1 + dtzidzi1
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi1 - dtxkdzi1
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi1 - dtykdzi1
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi1 - dtzkdzi1
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi1 + dtxkdzi1
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi1 + dtykdzi1
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi1 + dtzkdzi1
                  else if (i .eq. i2) then
                     hessx(1,i1) = hessx(1,i1) + si1*dtxdxi2 - dtxidxi2
                     hessx(2,i1) = hessx(2,i1) + si1*dtydxi2 - dtyidxi2
                     hessx(3,i1) = hessx(3,i1) + si1*dtzdxi2 - dtzidxi2
                     hessx(1,i2) = hessx(1,i2) + si2*dtxdxi2 + dtxidxi2
                     hessx(2,i2) = hessx(2,i2) + si2*dtydxi2 + dtyidxi2
                     hessx(3,i2) = hessx(3,i2) + si2*dtzdxi2 + dtzidxi2
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi2 - dtxkdxi2
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi2 - dtykdxi2
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi2 - dtzkdxi2
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi2 + dtxkdxi2
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi2 + dtykdxi2
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi2 + dtzkdxi2
                     hessy(1,i1) = hessy(1,i1) + si1*dtxdyi2 - dtxidyi2
                     hessy(2,i1) = hessy(2,i1) + si1*dtydyi2 - dtyidyi2
                     hessy(3,i1) = hessy(3,i1) + si1*dtzdyi2 - dtzidyi2
                     hessy(1,i2) = hessy(1,i2) + si2*dtxdyi2 + dtxidyi2
                     hessy(2,i2) = hessy(2,i2) + si2*dtydyi2 + dtyidyi2
                     hessy(3,i2) = hessy(3,i2) + si2*dtzdyi2 + dtzidyi2
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi2 - dtxkdyi2
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi2 - dtykdyi2
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi2 - dtzkdyi2
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi2 + dtxkdyi2
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi2 + dtykdyi2
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi2 + dtzkdyi2
                     hessz(1,i1) = hessz(1,i1) + si1*dtxdzi2 - dtxidzi2
                     hessz(2,i1) = hessz(2,i1) + si1*dtydzi2 - dtyidzi2
                     hessz(3,i1) = hessz(3,i1) + si1*dtzdzi2 - dtzidzi2
                     hessz(1,i2) = hessz(1,i2) + si2*dtxdzi2 + dtxidzi2
                     hessz(2,i2) = hessz(2,i2) + si2*dtydzi2 + dtyidzi2
                     hessz(3,i2) = hessz(3,i2) + si2*dtzdzi2 + dtzidzi2
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi2 - dtxkdzi2
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi2 - dtykdzi2
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi2 - dtzkdzi2
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi2 + dtxkdzi2
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi2 + dtykdzi2
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi2 + dtzkdzi2
                  end if
c
c     more energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     if (i .eq. i1) then
                        hessx(1,i1) = hessx(1,i1) + si1*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxi1 + si1*si1*d2taperxx
                        hessx(2,i1) = hessx(2,i1) + si1*dtapery*dedxi1
     &                         + si1*dtaperx*dedyi1 + si1*si1*d2taperxy
                        hessx(3,i1) = hessx(3,i1) + si1*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzi1 + si1*si1*d2taperxz
                        hessx(1,i2) = hessx(1,i2) + si2*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxi2 + si2*si1*d2taperxx
                        hessx(2,i2) = hessx(2,i2) + si2*dtapery*dedxi1
     &                         + si1*dtaperx*dedyi2 + si2*si1*d2taperxy
                        hessx(3,i2) = hessx(3,i2) + si2*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzi2 + si2*si1*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxk1 - sk1*si1*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtapery*dedxi1
     &                         + si1*dtaperx*dedyk1 - sk1*si1*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzk1 - sk1*si1*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxk2 - sk2*si1*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtapery*dedxi1
     &                         + si1*dtaperx*dedyk2 - sk2*si1*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzk2 - sk2*si1*d2taperxz
                        hessy(1,i1) = hessy(1,i1) + si1*dtaperx*dedyi1
     &                         + si1*dtapery*dedxi1 + si1*si1*d2taperxy
                        hessy(2,i1) = hessy(2,i1) + si1*dtapery*dedyi1
     &                         + si1*dtapery*dedyi1 + si1*si1*d2taperyy
                        hessy(3,i1) = hessy(3,i1) + si1*dtaperz*dedyi1
     &                         + si1*dtapery*dedzi1 + si1*si1*d2taperyz
                        hessy(1,i2) = hessy(1,i2) + si2*dtaperx*dedyi1
     &                         + si1*dtapery*dedxi2 + si2*si1*d2taperxy
                        hessy(2,i2) = hessy(2,i2) + si2*dtapery*dedyi1
     &                         + si1*dtapery*dedyi2 + si2*si1*d2taperyy
                        hessy(3,i2) = hessy(3,i2) + si2*dtaperz*dedyi1
     &                         + si1*dtapery*dedzi2 + si2*si1*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtaperx*dedyi1
     &                         + si1*dtapery*dedxk1 - sk1*si1*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyi1
     &                         + si1*dtapery*dedyk1 - sk1*si1*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtaperz*dedyi1
     &                         + si1*dtapery*dedzk1 - sk1*si1*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtaperx*dedyi1
     &                         + si1*dtapery*dedxk2 - sk2*si1*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyi1
     &                         + si1*dtapery*dedyk2 - sk2*si1*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtaperz*dedyi1
     &                         + si1*dtapery*dedzk2 - sk2*si1*d2taperyz
                        hessz(1,i1) = hessz(1,i1) + si1*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxi1 + si1*si1*d2taperxz
                        hessz(2,i1) = hessz(2,i1) + si1*dtapery*dedzi1
     &                         + si1*dtaperz*dedyi1 + si1*si1*d2taperyz
                        hessz(3,i1) = hessz(3,i1) + si1*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzi1 + si1*si1*d2taperzz
                        hessz(1,i2) = hessz(1,i2) + si2*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxi2 + si2*si1*d2taperxz
                        hessz(2,i2) = hessz(2,i2) + si2*dtapery*dedzi1
     &                         + si1*dtaperz*dedyi2 + si2*si1*d2taperyz
                        hessz(3,i2) = hessz(3,i2) + si2*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzi2 + si2*si1*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxk1 - sk1*si1*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtapery*dedzi1
     &                         + si1*dtaperz*dedyk1 - sk1*si1*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzk1 - sk1*si1*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxk2 - sk2*si1*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtapery*dedzi1
     &                         + si1*dtaperz*dedyk2 - sk2*si1*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzk2 - sk2*si1*d2taperzz
                     else if (i .eq. i2) then
                        hessx(1,i1) = hessx(1,i1) + si1*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxi1 + si1*si2*d2taperxx
                        hessx(2,i1) = hessx(2,i1) + si1*dtapery*dedxi2
     &                         + si2*dtaperx*dedyi1 + si1*si2*d2taperxy
                        hessx(3,i1) = hessx(3,i1) + si1*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzi1 + si1*si2*d2taperxz
                        hessx(1,i2) = hessx(1,i2) + si2*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxi2 + si2*si2*d2taperxx
                        hessx(2,i2) = hessx(2,i2) + si2*dtapery*dedxi2
     &                         + si2*dtaperx*dedyi2 + si2*si2*d2taperxy
                        hessx(3,i2) = hessx(3,i2) + si2*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzi2 + si2*si2*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxk1 - sk1*si2*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtapery*dedxi2
     &                         + si2*dtaperx*dedyk1 - sk1*si2*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzk1 - sk1*si2*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxk2 - sk2*si2*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtapery*dedxi2
     &                         + si2*dtaperx*dedyk2 - sk2*si2*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzk2 - sk2*si2*d2taperxz
                        hessy(1,i1) = hessy(1,i1) + si1*dtaperx*dedyi2
     &                         + si2*dtapery*dedxi1 + si1*si2*d2taperxy
                        hessy(2,i1) = hessy(2,i1) + si1*dtapery*dedyi2
     &                         + si2*dtapery*dedyi1 + si1*si2*d2taperyy
                        hessy(3,i1) = hessy(3,i1) + si1*dtaperz*dedyi2
     &                         + si2*dtapery*dedzi1 + si1*si2*d2taperyz
                        hessy(1,i2) = hessy(1,i2) + si2*dtaperx*dedyi2
     &                         + si2*dtapery*dedxi2 + si2*si2*d2taperxy
                        hessy(2,i2) = hessy(2,i2) + si2*dtapery*dedyi2
     &                         + si2*dtapery*dedyi2 + si2*si2*d2taperyy
                        hessy(3,i2) = hessy(3,i2) + si2*dtaperz*dedyi2
     &                         + si2*dtapery*dedzi2 + si2*si2*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtaperx*dedyi2
     &                         + si2*dtapery*dedxk1 - sk1*si2*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyi2
     &                         + si2*dtapery*dedyk1 - sk1*si2*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtaperz*dedyi2
     &                         + si2*dtapery*dedzk1 - sk1*si2*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtaperx*dedyi2
     &                         + si2*dtapery*dedxk2 - sk2*si2*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyi2
     &                         + si2*dtapery*dedyk2 - sk2*si2*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtaperz*dedyi2
     &                         + si2*dtapery*dedzk2 - sk2*si2*d2taperyz
                        hessz(1,i1) = hessz(1,i1) + si1*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxi1 + si1*si2*d2taperxz
                        hessz(2,i1) = hessz(2,i1) + si1*dtapery*dedzi2
     &                         + si2*dtaperz*dedyi1 + si1*si2*d2taperyz
                        hessz(3,i1) = hessz(3,i1) + si1*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzi1 + si1*si2*d2taperzz
                        hessz(1,i2) = hessz(1,i2) + si2*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxi2 + si2*si2*d2taperxz
                        hessz(2,i2) = hessz(2,i2) + si2*dtapery*dedzi2
     &                         + si2*dtaperz*dedyi2 + si2*si2*d2taperyz
                        hessz(3,i2) = hessz(3,i2) + si2*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzi2 + si2*si2*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxk1 - sk1*si2*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtapery*dedzi2
     &                         + si2*dtaperz*dedyk1 - sk1*si2*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzk1 - sk1*si2*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxk2 - sk2*si2*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtapery*dedzi2
     &                         + si2*dtaperz*dedyk2 - sk2*si2*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzk2 - sk2*si2*d2taperzz
                     end if
                  end if
               end if
            end do
   20       continue
         end do
   30    continue
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
c     ##  subroutine edipole3  --  dipole-dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edipole3" calculates the dipole-dipole interaction energy;
c     also partitions the energy among the atoms
c
c
      subroutine edipole3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'dipole.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moment.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,i1,i2,k1,k2
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xq,yq,zq,xr,yr,zr
      real*8 f,fi,fik,fgrp
      real*8 e,ri2,rk2,rirkr3
      real*8 doti,dotk,dotp
      real*8 r,r2,r3,r4,r5,taper
      real*8 xsum,ysum,zsum
      logical header,huge,proceed
c
c
c     zero out the overall dipole interaction energy contribution
c     and set up the constants for the calculation
c
      ned = 0
      ed = 0.0d0
      do i = 1, n
         aed(i) = 0.0d0
      end do
      if (ndipole .eq. 0)  return
      header = .true.
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye**2 * dielec)
      call switch ('DIPOLE')
c
c     compute the components of the dipole moment
c
      xsum = 0.0d0
      ysum = 0.0d0
      zsum = 0.0d0
      do i = 1, ndipole
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         xi = x(i1) - x(i2)
         yi = y(i1) - y(i2)
         zi = z(i1) - z(i2)
         fik = bdpl(i) / sqrt(xi*xi + yi*yi + zi*zi)
         xsum = xsum + fik*xi
         ysum = ysum + fik*yi
         zsum = zsum + fik*zi
      end do
      xdipole = xdipole + xsum
      ydipole = ydipole + ysum
      zdipole = zdipole + zsum
c
c     compute and partition the dipole interaction energy
c
      do i = 1, ndipole-1
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*sdpl(i)
         yq = y(i1) + yi*sdpl(i)
         zq = z(i1) + zi*sdpl(i)
         fi = f * bdpl(i)
         do k = i+1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,4,i1,i2,k1,k2,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
            if (proceed)  proceed = (k1.ne.i1 .and. k1.ne.i2 .and.
     &                                 k2.ne.i1 .and. k2.ne.i2)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = xq - x(k1) - xk*sdpl(k)
               yr = yq - y(k1) - yk*sdpl(k)
               zr = zq - z(k1) - zk*sdpl(k)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr3 = sqrt(ri2*rk2*r2) * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(k)
                  e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
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
c     increment the overall dipole-dipole energy component
c
                  ned = ned + 1
                  ed = ed + e
                  aed(i1) = aed(i1) + 0.25d0*e
                  aed(i2) = aed(i2) + 0.25d0*e
                  aed(k1) = aed(k1) + 0.25d0*e
                  aed(k2) = aed(k2) + 0.25d0*e
c
c     print a warning if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Dipole-Dipole',
     &                             ' Interactions :',
     &                          //,' Type',12x,'Dipole 1',14x,
     &                            'Dipole 2',6x,'Distance',
     &                             6x,'Energy',/)
                     end if
                     write (iout,20)  i1,name(i1),i2,name(i2),
     &                                k1,name(k1),k2,name(k2),
     &                                sqrt(r2),e
   20                format (' Dipole   ',i5,'-',a3,1x,i5,'-',a3,
     &                          ' / ',i5,'-',a3,1x,i5,'-',a3,
     &                          f10.4,f12.4)
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
      do i = 1, ndipole
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*sdpl(i)
         yq = y(i1) + yi*sdpl(i)
         zq = z(i1) + zi*sdpl(i)
         fi = f * bdpl(i)
         do k = i, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,4,i1,i2,k1,k2,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  xr = xq - x(k1) - xk*sdpl(k)
                  yr = yq - y(k1) - yk*sdpl(k)
                  zr = zq - z(k1) - zk*sdpl(k)
                  call image (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rk2 = xk*xk + yk*yk + zk*zk
                     rirkr3 = sqrt(ri2*rk2*r2) * r2
                     dotp = xi*xk + yi*yk + zi*zk
                     doti = xi*xr + yi*yr + zi*zr
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = fi * bdpl(k)
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
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
c     increment the overall dipole-dipole energy component
c
                     ned = ned + 1
                     if (i .eq. k) then
                        ed = ed + 0.5d0*e
                        aed(i1) = aed(i1) + 0.25d0*e
                        aed(i2) = aed(i2) + 0.25d0*e
                     else
                        ed = ed + e
                        aed(i1) = aed(i1) + 0.25d0*e
                        aed(i2) = aed(i2) + 0.25d0*e
                        aed(k1) = aed(k1) + 0.25d0*e
                        aed(k2) = aed(k2) + 0.25d0*e
                     end if
c
c     print a warning if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 10.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Dipole-Dipole',
     &                                ' Interactions :',
     &                             //,' Type',12x,'Dipole 1',14x,
     &                               'Dipole 2',6x,'Distance',
     &                                6x,'Energy',/)
                        end if
                        write (iout,40)  i1,name(i1),i2,name(i2),
     &                                   k1,name(k1),k2,name(k2),
     &                                   sqrt(r2),e
  40                   format (' Dipole   ',i5,'-',a3,1x,i5,'-',a3,
     &                            ' / ',i5,'-',a3,1x,i5,'-',a3,
     &                            f10.4,f12.4)
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
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine egauss  --  Gaussian van der Waals energy  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "egauss" calculates the van der Waals interaction energy
c     using a Gaussian expansion approximation
c
c
      subroutine egauss
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer i,j,k,ii,kk,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rik2,rdn,eps,rad2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 expcut,expterm
      real*8 width,wterm
      real*8 t1,t2,a(4),b(4)
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
c     set cutoff distances and switching function coefficients
c
      call switch ('VDW')
      expcut = -50.0d0
c
c     scale deformation time by the diffusion coefficient
c
      if (use_deform)  width = 4.0d0 * diffv * deform
      if (use_gda)  wterm = (2.0d0/3.0d0) * diffv
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
               rik2 = xr*xr + yr*yr + zr*zr
               if (rik2 .le. off2) then
                  eps = epsilon(kt,it)
                  rad2 = radmin(kt,it)**2
                  do j = 1, ngauss
                     a(j) = igauss(1,j) * eps
                     b(j) = igauss(2,j) / rad2
                     if (skip(k) .eq. -i)  a(j) = a(j) / vdwscale
                  end do
                  e = 0.0d0
c
c     potential smoothing via diffusion equation or density annealing
c
                  if (use_gda)  width = wterm * (m2(i)+m2(k))
                  do j = 1, ngauss
                     t1 = 1.0d0 + b(j)*width
                     t2 = sqrt(t1**3)
                     expterm = -b(j) * rik2 / t1
                     if (expterm .gt. expcut)
     &                  e = e + (a(j)/t2)*exp(expterm)
                  end do
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
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egauss1  --  Gaussian vdw energy & derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egauss1" calculates the van der Waals interaction energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a Gaussian expansion approximation
c
c
      subroutine egauss1
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer i,j,k,ii,kk,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rik,rik2,rdn,eps,rad2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 dedx,dedy,dedz,de
      real*8 expcut,expterm
      real*8 width,wterm
      real*8 t1,t2,a(4),b(4)
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
      expcut = -50.0d0
c
c     scale deformation time by the diffusion coefficient
c
      if (use_deform)  width = 4.0d0 * diffv * deform
      if (use_gda)  wterm = (2.0d0/3.0d0) * diffv
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
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute energy and derivatives for this interaction
c
               if (rik2 .le. off2) then
                  eps = epsilon(kt,it)
                  rad2 = radmin(kt,it)**2
                  do j = 1, ngauss
                     a(j) = igauss(1,j) * eps
                     b(j) = igauss(2,j) / rad2
                     if (skip(k) .eq. -i)  a(j) = a(j) / vdwscale
                  end do
                  e = 0.0d0
                  de = 0.0d0
                  rik = sqrt(rik2)
c
c     potential smoothing via diffusion equation or density annealing
c
                  if (use_gda)  width = wterm * (m2(i)+m2(k))
                  do j = 1, ngauss
                     t1 = 1.0d0 + b(j)*width
                     t2 = sqrt(t1**3)
                     expterm = -b(j) * rik2 / t1
                     if (expterm .gt. expcut) then
                        expterm = (a(j)/t2)*exp(expterm)
                        e = e + expterm
                        de = de - (2.0d0*b(j)*rik/t1)*expterm
                     end if
                  end do
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egauss2  --  atom-by-atom Gaussian vdw Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egauss2" calculates the van der Waals second derivatives
c     for a single atom at a time using a Gaussian approximation
c
c
      subroutine egauss2 (iatom,xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'group.i'
      include 'hessn.i'
      include 'shunt.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer iatom,i,j,k,ii,kk,iv,kv
      integer it,kt,skip(maxatm)
      integer nuse,use(5)
      real*8 de,d2e,rik,rik2,eps,rad2,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 redi2,rediv2,rediiv
      real*8 redik,redivk,redikv,redivkv
      real*8 d2edx,d2edy,d2edz,term(3,3)
      real*8 expcut,expterm
      real*8 width,wterm
      real*8 t1,t2,a(4),b(4)
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
      expcut = -50.0d0
c
c     scale deformation time by the diffusion coefficient
c
      if (use_deform)  width = 4.0d0 * diffv * deform
      if (use_gda)  wterm = (2.0d0/3.0d0) * diffv
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
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute Hessian elements for this interaction
c
               if (rik2 .le. off2) then
                  eps = epsilon(kt,it)
                  rad2 = radmin(kt,it)**2
                  do j = 1, ngauss
                     a(j) = igauss(1,j) * eps
                     b(j) = igauss(2,j) / rad2
                     if (skip(k) .eq. -i)  a(j) = a(j) / vdwscale
                  end do
                  de = 0.0d0
                  d2e = 0.0d0
                  rik = sqrt(rik2)
c
c     potential smoothing via diffusion equation or density annealing
c
                  if (use_gda)  width = wterm * (m2(i)+m2(k))
                  do j = 1, ngauss
                     t1 = 1.0d0 + b(j)*width
                     t2 = sqrt(t1**3)
                     expterm = -b(j) * rik2 / t1
                     if (expterm .gt. expcut) then
                        expterm = (a(j)*b(j)/(t2*t1))*exp(expterm)
                        de = de - 2.0d0*rik*expterm
                        d2e = d2e + (4.0d0*b(j)*rik2/t1-2.0d0)*expterm
                     end if
                  end do
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
c     ##  subroutine egauss3  --  Gaussian vdw energy & analysis  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "egauss3" calculates the van der Waals interaction energy
c     using a Gaussian approximation and also partitions the
c     energy among the atoms
c
c
      subroutine egauss3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'warp.i'
      integer i,j,k,ii,kk,iv,kv
      integer it,kt,skip(maxatm)
      real*8 e,rik2,rdn,eps,rad2,rv,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 expcut,expterm
      real*8 width,wterm
      real*8 t1,t2,a(4),b(4)
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
      expcut = -50.0d0
c
c     scale deformation time by the diffusion coefficient
c
      if (use_deform)  width = 4.0d0 * diffv * deform
      if (use_gda)  wterm = (2.0d0/3.0d0) * diffv
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
               rik2 = xr*xr + yr*yr + zr*zr
c
c     compute the energy contribution for this interaction
c
               if (rik2 .le. off2) then
                  eps = epsilon(kt,it)
                  rad2 = radmin(kt,it)**2
                  do j = 1, ngauss
                     a(j) = igauss(1,j) * eps
                     b(j) = igauss(2,j) / rad2
                     if (skip(k) .eq. -i)  a(j) = a(j) / vdwscale
                  end do
                  e = 0.0d0
c
c     potential smoothing via diffusion equation or density annealing
c
                  if (use_gda)  width = wterm * (m2(i)+m2(k))
                  do j = 1, ngauss
                     t1 = 1.0d0 + b(j)*width
                     t2 = sqrt(t1**3)
                     expterm = -b(j) * rik2 / t1
                     if (expterm .gt. expcut)
     &                  e = e + (a(j)/t2)*exp(expterm)
                  end do
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the total van der Waals energy components
c
                  if (skip(k) .eq. -i) then
                     ne14 = ne14 + 1
                     ae14(i) = ae14(i) + 0.5d0*e
                     ae14(k) = ae14(k) + 0.5d0*e
                     e14 = e14 + e
                  else
                     nev = nev + 1
                     aev(i) = aev(i) + 0.5d0*e
                     aev(k) = aev(k) + 0.5d0*e
                     ev = ev + e
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
                     rv = radmin(kt,it)
                     write (iout,20)  i,name(i),k,name(k),
     &                                rv,sqrt(rik2),e
   20                format (' VDW-Gaus ',i5,'-',a3,1x,i5,'-',a3,
     &                         12x,2f10.4,f12.4)
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
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine egeom  --  geometric restraint energy terms  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "egeom" calculates the energy due to restraints on atomic
c     positions, interatomic distances, dihedral angles, and
c     Gaussian weighted molecular compactness
c
c
      subroutine egeom
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'energi.i'
      include 'math.i'
      include 'molcul.i'
      include 'restrn.i'
      include 'usage.i'
      integer i,j,k,ia,ib,ic,id
      real*8 e,xr,yr,zr,rik,rik2,dt,dt2
      real*8 angle,target,force,cosine,sine
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 df1,df2,tf1,tf2,t1,t2
      real*8 xi,yi,zi,ri
      real*8 a,b,buffer
      real*8 r,r2,r6,r12
      logical intermol
c
c
c     zero out the geometric restraint energy terms
c
      eg = 0.0d0
c
c     compute the pseudoenergy for position restraints
c
      do j = 1, npfix
         i = ipfix(j)
         if (use(i)) then
            xr = x(i) - xpfix(j)
            yr = y(i) - ypfix(j)
            zr = z(i) - zpfix(j)
            dt2 = xr*xr + yr*yr + zr*zr
            force = pfix(j)
            e = force * dt2
            eg = eg + e
         end if
      end do
c
c     compute the pseudoenergy for distance restraints
c
      do j = 1, ndfix
         i = idfix(1,j)
         k = idfix(2,j)
         if (use(i) .or. use(k)) then
            xr = x(i) - x(k)
            yr = y(i) - y(k)
            zr = z(i) - z(k)
            intermol = (molcule(i) .ne. molcule(k))
            if (use_bounds .and. intermol)  call image (xr,yr,zr,0)
            rik2 = xr*xr + yr*yr + zr*zr
            rik = sqrt(rik2)
            df1 = dfix(1,j)
            df2 = dfix(2,j)
            target = rik
            if (rik .lt. df1)  target = df1
            if (rik .gt. df2)  target = df2
            force = dfix(3,j)
            dt = rik - target
            dt2 = dt * dt
            e = force * dt2
            eg = eg + e
         end if
      end do
c
c     compute the pseudoenergy for dihedral angle restraints
c
      do i = 1, ntfix
         ia = itfix(1,i)
         ib = itfix(2,i)
         ic = itfix(3,i)
         id = itfix(4,i)
         if (use(ia) .or. use(ib) .or. use(ic) .or. use(id)) then
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
               tf1 = tfix(1,i)
               tf2 = tfix(2,i)
               if (angle.gt.tf1 .and. angle.lt.tf2) then
                  target = angle
               else if (angle.gt.tf1 .and. tf1.gt.tf2) then
                  target = angle
               else if (angle.lt.tf2 .and. tf1.gt.tf2) then
                  target = angle
               else
                  t1 = angle - tf1
                  t2 = angle - tf2
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     target = tf1
                  else
                     target = tf2
                  end if
               end if
               force = tfix(3,i)
               dt = angle - target
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
               dt2 = dt * dt
               e = force * dt2
               eg = eg + e
            end if
         end if
      end do
c
c     compute the energy for a shallow Gaussian basin restraint
c
      if (use_basin) then
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do k = i+1, n
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               rik2 = xr*xr + yr*yr + zr*zr
               e = depth*exp(-width*rik2) - depth
               eg = eg + e
            end do
         end do
      end if
c
c     compute the energy for a spherical droplet boundary restraint
c
      if (use_wall) then
         buffer = 2.5d0
         a = 2048.0d0
         b = 64.0d0
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = sqrt(xi**2 + yi**2 + zi**2)
            r = rwall + buffer - ri
            r2 = r * r
            r6 = r2 * r2 * r2
            r12 = r6 * r6
            e = a/r12 - b/r6
            eg = eg + e
         end do
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine egeom1  --  restraint energy & derivatives  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "egeom1" calculates the potential energy and first derivatives
c     with respect to Cartesian coordinates for restraints on atom
c     positions, interatomic distances and dihedral angles
c
c
      subroutine egeom1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inter.i'
      include 'molcul.i'
      include 'math.i'
      include 'restrn.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,ia,ib,ic,id
      real*8 e,xr,yr,zr,rik,rik2
      real*8 de,dedx,dedy,dedz,dt,dt2,deddt
      real*8 angle,target,force,cosine,sine
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb,dedphi
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 df1,df2,tf1,tf2,t1,t2
      real*8 xi,yi,zi,ri
      real*8 a,b,buffer
      real*8 r,r2,r6,r12
      logical intermol
c
c
c     zero out the restraint energy term and first derivatives
c
      eg = 0.0d0
      do i = 1, n
         deg(1,i) = 0.0d0
         deg(2,i) = 0.0d0
         deg(3,i) = 0.0d0
      end do
c
c     compute the pseudoenergy for position restraints
c
      do j = 1, npfix
         i = ipfix(j)
         if (use(i)) then
            xr = x(i) - xpfix(j)
            yr = y(i) - ypfix(j)
            zr = z(i) - zpfix(j)
            dt2 = xr*xr + yr*yr + zr*zr
            force = pfix(j)
            e = force * dt2
            deddt = 2.0d0 * force
            dedx = deddt * xr
            dedy = deddt * yr
            dedz = deddt * zr
c
c     increment the total energy and first derivatives
c
            eg = eg + e
            deg(1,i) = deg(1,i) + dedx
            deg(2,i) = deg(2,i) + dedy
            deg(3,i) = deg(3,i) + dedz
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
c
c     compute the pseudoenergy for distance restraints
c
      do j = 1, ndfix
         i = idfix(1,j)
         k = idfix(2,j)
         if (use(i) .or. use(k)) then
            xr = x(i) - x(k)
            yr = y(i) - y(k)
            zr = z(i) - z(k)
            intermol = (molcule(i) .ne. molcule(k))
            if (use_bounds .and. intermol)  call image (xr,yr,zr,0)
            rik2 = xr*xr + yr*yr + zr*zr
            rik = sqrt(rik2)
            df1 = dfix(1,j)
            df2 = dfix(2,j)
            target = rik
            if (rik .lt. df1)  target = df1
            if (rik .gt. df2)  target = df2
            force = dfix(3,j)
            dt = rik - target
            dt2 = dt * dt
            e = force * dt2
            deddt = 2.0d0 * force * dt
c
c     compute chain rule terms needed for derivatives
c
            de = deddt / rik
            dedx = de * xr
            dedy = de * yr
            dedz = de * zr
c
c     increment the total energy and first derivatives
c
            eg = eg + e
            deg(1,i) = deg(1,i) + dedx
            deg(2,i) = deg(2,i) + dedy
            deg(3,i) = deg(3,i) + dedz
            deg(1,k) = deg(1,k) - dedx
            deg(2,k) = deg(2,k) - dedy
            deg(3,k) = deg(3,k) - dedz
c
c     increment the total intermolecular energy
c
            if (intermol) then
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
      end do
c
c     compute the pseudoenergy and derivs for dihedral restraints;
c     for each angle we first calculate the cosine and sine of the
c     dihedral between atoms ia-ib-ic-id
c
      do i = 1, ntfix
         ia = itfix(1,i)
         ib = itfix(2,i)
         ic = itfix(3,i)
         id = itfix(4,i)
         if (use(ia) .or. use(ib) .or. use(ic) .or. use(id)) then
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     calculate the pseudoenergy contribution for this angle
c
               tf1 = tfix(1,i)
               tf2 = tfix(2,i)
               if (angle.gt.tf1 .and. angle.lt.tf2) then
                  target = angle
               else if (angle.gt.tf1 .and. tf1.gt.tf2) then
                  target = angle
               else if (angle.lt.tf2 .and. tf1.gt.tf2) then
                  target = angle
               else
                  t1 = angle - tf1
                  t2 = angle - tf2
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     target = tf1
                  else
                     target = tf2
                  end if
               end if
               force = tfix(3,i)
               dt = angle - target
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
               dt2 = dt * dt
               e = force * dt2
               dedphi = 2.0d0 * radian * force * dt
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = yca*dphidzt - zca*dphidyt
     &                       + zdc*dphidyu - ydc*dphidzu
               dphidyib = zca*dphidxt - xca*dphidzt
     &                       + xdc*dphidzu - zdc*dphidxu
               dphidzib = xca*dphidyt - yca*dphidxt
     &                       + ydc*dphidxu - xdc*dphidyu
               dphidxic = zba*dphidyt - yba*dphidzt
     &                       + ydb*dphidzu - zdb*dphidyu
               dphidyic = xba*dphidzt - zba*dphidxt
     &                       + zdb*dphidxu - xdb*dphidzu
               dphidzic = yba*dphidxt - xba*dphidyt
     &                       + xdb*dphidyu - ydb*dphidxu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     compute first derivative components for this angle
c
               dedxia = dedphi * dphidxia
               dedyia = dedphi * dphidyia
               dedzia = dedphi * dphidzia
               dedxib = dedphi * dphidxib
               dedyib = dedphi * dphidyib
               dedzib = dedphi * dphidzib
               dedxic = dedphi * dphidxic
               dedyic = dedphi * dphidyic
               dedzic = dedphi * dphidzic
               dedxid = dedphi * dphidxid
               dedyid = dedphi * dphidyid
               dedzid = dedphi * dphidzid
c
c     increment the overall extra energy term and derivatives
c
               eg = eg + e
               deg(1,ia) = deg(1,ia) + dedxia
               deg(2,ia) = deg(2,ia) + dedyia
               deg(3,ia) = deg(3,ia) + dedzia
               deg(1,ib) = deg(1,ib) + dedxib
               deg(2,ib) = deg(2,ib) + dedyib
               deg(3,ib) = deg(3,ib) + dedzib
               deg(1,ic) = deg(1,ic) + dedxic
               deg(2,ic) = deg(2,ic) + dedyic
               deg(3,ic) = deg(3,ic) + dedzic
               deg(1,id) = deg(1,id) + dedxid
               deg(2,id) = deg(2,id) + dedyid
               deg(3,id) = deg(3,id) + dedzid
c
c     increment the total intermolecular energy
c
               if (molcule(ia).ne.molcule(ib) .or.
     &             molcule(ia).ne.molcule(ic) .or.
     &             molcule(ia).ne.molcule(id)) then
                  einter = einter + e
               end if
c
c     increment the virial for use in pressure computation
c
               if (isobaric) then
                  virx = virx - xba*dedxia + xdc*dedxid
     &                      + xcb*(dedxic-dedxid)
                  viry = viry - yba*dedyia + ydc*dedyid
     &                      + ycb*(dedyic-dedyid)
                  virz = virz - zba*dedzia + zdc*dedzid
     &                      + zcb*(dedzic-dedzid)
               end if
            end if
         end if
      end do
c
c     compute the energy for a shallow Gaussian basin restraint
c
      if (use_basin) then
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do k = i+1, n
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               rik2 = xr*xr + yr*yr + zr*zr
               e = depth * exp(-width*rik2)
               de = -2.0d0 * width * e
               e = e - depth
               dedx = de * xr
               dedy = de * yr
               dedz = de * zr
c
c     increment the total energy and first derivatives
c
               eg = eg + e
               deg(1,i) = deg(1,i) + dedx
               deg(2,i) = deg(2,i) + dedy
               deg(3,i) = deg(3,i) + dedz
               deg(1,k) = deg(1,k) - dedx
               deg(2,k) = deg(2,k) - dedy
               deg(3,k) = deg(3,k) - dedz
c
c     increment the virial for use in pressure computation
c
               if (isobaric) then
                  virx = virx + xr*dedx
                  viry = viry + yr*dedy
                  virz = virz + zr*dedz
               end if
            end do
         end do
      end if
c
c     compute the energy for a spherical droplet boundary restraint
c
      if (use_wall) then
         buffer = 2.5d0
         a = 2048.0d0
         b = 64.0d0
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = sqrt(xi**2 + yi**2 + zi**2)
            r = rwall + buffer - ri
            r2 = r * r
            r6 = r2 * r2 * r2
            r12 = r6 * r6
            e = a/r12 - b/r6
            if (ri .eq. 0.0d0)  ri = 1.0d0
            de = (12.0d0*a/r12 - 6.0d0*b/r6) / (r*ri)
            dedx = de * xi
            dedy = de * yi
            dedz = de * zi
c
c     increment the total energy and first derivatives
c
            eg = eg + e
            deg(1,i) = deg(1,i) + dedx
            deg(2,i) = deg(2,i) + dedy
            deg(3,i) = deg(3,i) + dedz
c
c     increment the virial for use in pressure computation
c
            if (isobaric) then
               xr = r * xi/ri
               yr = r * yi/ri
               zr = r * zi/ri
               virx = virx + xr*dedx
               viry = viry + yr*dedy
               virz = virz + zr*dedz
            end if
         end do
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine egeom2  --  atom-by-atom restraint Hessian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "egeom2" calculates second derivatives of any restraints on
c     atom positions, interatomic distances and dihedral angles
c
c
      subroutine egeom2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'deriv.i'
      include 'hessn.i'
      include 'math.i'
      include 'molcul.i'
      include 'restrn.i'
      integer i,j,k,ia,ib,ic,id
      integer ipos,idist,itors
      real*8 xr,yr,zr,target,force
      real*8 rik,rik2,dt,dt2,deddt,d2eddt2
      real*8 de,term,termx,termy,termz,d2e(3,3)
      real*8 dedphi,d2edphi2,cosine,sine,angle
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 df1,df2,tf1,tf2,t1,t2
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 dphidxibt,dphidyibt,dphidzibt
      real*8 dphidxibu,dphidyibu,dphidzibu
      real*8 dphidxict,dphidyict,dphidzict
      real*8 dphidxicu,dphidyicu,dphidzicu
      real*8 dxiaxia,dyiayia,dziazia,dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic,dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia,dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic,dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib,dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib,dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic,dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid,dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid,dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic,dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid,dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid,dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid,dzicxid,dzicyid,dziczid
      real*8 xi,yi,zi,dedr,d2edr2,expterm
      real*8 ri,ri2,r,r2,r6,r12,a,b,buffer
      logical intermol
c
c
c     compute the Hessian elements for position restraints
c
      do ipos = 1, npfix
         if (ipfix(ipos) .eq. i) then
            xr = x(i) - xpfix(ipos)
            yr = y(i) - ypfix(ipos)
            zr = z(i) - zpfix(ipos)
            dt2 = xr*xr + yr*yr + zr*zr
            force = pfix(ipos)
            deddt = 2.0d0 * force
            hessx(1,i) = hessx(1,i) + deddt
            hessy(2,i) = hessy(2,i) + deddt
            hessz(3,i) = hessz(3,i) + deddt
         end if
      end do
c
c     compute the Hessian elements for distance restraints;
c     note that the Hessian is discontinuous when an upper and
c     lower bound range is used instead of a single distance
c
      do idist = 1, ndfix
         ia = idfix(1,idist)
         ib = idfix(2,idist)
         if (i.eq.ia .or. i.eq.ib) then
            if (i .eq. ib) then
               ib = ia
               ia = i
            end if
            xr = x(ia) - x(ib)
            yr = y(ia) - y(ib)
            zr = z(ia) - z(ib)
            intermol = (molcule(ia) .ne. molcule(ib))
            if (use_bounds .and. intermol)  call image (xr,yr,zr,0)
            rik2 = xr*xr + yr*yr + zr*zr
            rik = sqrt(rik2)
            df1 = dfix(1,idist)
            df2 = dfix(2,idist)
            target = rik
            if (rik .lt. df1)  target = df1
            if (rik .gt. df2)  target = df2
            force = dfix(3,idist)
            dt = rik - target
            deddt = 2.0d0 * force * dt
            d2eddt2 = 2.0d0 * force
c
c     set the chain rule terms for the Hessian elements
c
            de = deddt / rik
            term = (d2eddt2-de) / rik2
            termx = term * xr
            termy = term * yr
            termz = term * zr
            d2e(1,1) = termx*xr + de
            d2e(1,2) = termx*yr
            d2e(1,3) = termx*zr
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yr + de
            d2e(2,3) = termy*zr
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zr + de
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
c
c     compute the Hessian elements for dihedral restraints;
c     for each angle we first calculate the cosine and sine
c     of the dihedral between atoms ia-ib-ic-id
c
      do itors = 1, ntfix
         ia = itfix(1,itors)
         ib = itfix(2,itors)
         ic = itfix(3,itors)
         id = itfix(4,itors)
         if (i.ne.ia .and. i.ne.ib .and.
     &       i.ne.ic .and. i.ne.id)  goto 10
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
         xba = xib - xia
         yba = yib - yia
         zba = zib - zia
         xcb = xic - xib
         ycb = yic - yib
         zcb = zic - zib
         xdc = xid - xic
         ydc = yid - yic
         zdc = zid - zic
         xt = yba*zcb - ycb*zba
         yt = zba*xcb - zcb*xba
         zt = xba*ycb - xcb*yba
         xu = ycb*zdc - ydc*zcb
         yu = zcb*xdc - zdc*xcb
         zu = xcb*ydc - xdc*ycb
         xtu = yt*zu - yu*zt
         ytu = zt*xu - zu*xt
         ztu = xt*yu - xu*yt
         rt2 = xt*xt + yt*yt + zt*zt
         ru2 = xu*xu + yu*yu + zu*zu
         rtru = sqrt(rt2 * ru2)
         if (rtru .eq. 0.0d0)  goto 10
         rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
         cosine = (xt*xu + yt*yu + zt*zu) / rtru
         sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
         if (sine .lt. 0.0d0)  angle = -angle
c
c     calculate the pseudoenergy master chain rule terms
c
         tf1 = tfix(1,itors)
         tf2 = tfix(2,itors)
         if (angle.gt.tf1 .and. angle.lt.tf2) then
            target = angle
         else if (angle.gt.tf1 .and. tf1.gt.tf2) then
            target = angle
         else if (angle.lt.tf2 .and. tf1.gt.tf2) then
            target = angle
         else
            t1 = angle - tf1
            t2 = angle - tf2
            if (t1 .gt. 180.0d0) then
               t1 = t1 - 360.0d0
            else if (t1 .lt. -180.0d0) then
               t1 = t1 + 360.0d0
            end if
            if (t2 .gt. 180.0d0) then
               t2 = t2 - 360.0d0
            else if (t2 .lt. -180.0d0) then
               t2 = t2 + 360.0d0
            end if
            if (abs(t1) .lt. abs(t2)) then
               target = tf1
            else
               target = tf2
            end if
         end if
         force = tfix(3,itors)
         dt = angle - target
         if (dt .gt. 180.0d0) then
            dt = dt - 360.0d0
         else if (dt .lt. -180.0d0) then
            dt = dt + 360.0d0
         end if
         dedphi = 2.0d0 * radian * force * dt
         d2edphi2 = 2.0d0 * radian**2 * force
c
c     abbreviations for first derivative chain rule terms
c
         xca = xic - xia
         yca = yic - yia
         zca = zic - zia
         xdb = xid - xib
         ydb = yid - yib
         zdb = zid - zib
         dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
         dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
         dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
         dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
         dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
         dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
         xycb2 = xcb*xcb + ycb*ycb
         xzcb2 = xcb*xcb + zcb*zcb
         yzcb2 = ycb*ycb + zcb*zcb
         rcbxt = -2.0d0 * rcb * dphidxt
         rcbyt = -2.0d0 * rcb * dphidyt
         rcbzt = -2.0d0 * rcb * dphidzt
         rcbt2 = rcb * rt2
         rcbxu = 2.0d0 * rcb * dphidxu
         rcbyu = 2.0d0 * rcb * dphidyu
         rcbzu = 2.0d0 * rcb * dphidzu
         rcbu2 = rcb * ru2
         dphidxibt = yca*dphidzt - zca*dphidyt
         dphidxibu = zdc*dphidyu - ydc*dphidzu
         dphidyibt = zca*dphidxt - xca*dphidzt
         dphidyibu = xdc*dphidzu - zdc*dphidxu
         dphidzibt = xca*dphidyt - yca*dphidxt
         dphidzibu = ydc*dphidxu - xdc*dphidyu
         dphidxict = zba*dphidyt - yba*dphidzt
         dphidxicu = ydb*dphidzu - zdb*dphidyu
         dphidyict = xba*dphidzt - zba*dphidxt
         dphidyicu = zdb*dphidxu - xdb*dphidzu
         dphidzict = yba*dphidxt - xba*dphidyt
         dphidzicu = xdb*dphidyu - ydb*dphidxu
c
c     chain rule terms for first derivative components
c
         dphidxia = zcb*dphidyt - ycb*dphidzt
         dphidyia = xcb*dphidzt - zcb*dphidxt
         dphidzia = ycb*dphidxt - xcb*dphidyt
         dphidxib = dphidxibt + dphidxibu
         dphidyib = dphidyibt + dphidyibu
         dphidzib = dphidzibt + dphidzibu
         dphidxic = dphidxict + dphidxicu
         dphidyic = dphidyict + dphidyicu
         dphidzic = dphidzict + dphidzicu
         dphidxid = zcb*dphidyu - ycb*dphidzu
         dphidyid = xcb*dphidzu - zcb*dphidxu
         dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     chain rule terms for second derivative components
c
         dxiaxia = rcbxt*dphidxia
         dxiayia = rcbxt*dphidyia - zcb*rcb/rt2
         dxiazia = rcbxt*dphidzia + ycb*rcb/rt2
         dxiaxib = rcbxt*dphidxibt + xcb*(zca*ycb-yca*zcb)/rcbt2
         dxiayib = rcbxt*dphidyibt + dphidzt
     &                + (xca*zcb*xcb+zca*yzcb2)/rcbt2
         dxiazib = rcbxt*dphidzibt - dphidyt
     &                - (xca*ycb*xcb+yca*yzcb2)/rcbt2
         dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2
         dxiayic = rcbxt*dphidyict - dphidzt
     &                - (xba*zcb*xcb+zba*yzcb2)/rcbt2
         dxiazic = rcbxt*dphidzict + dphidyt
     &                + (xba*ycb*xcb+yba*yzcb2)/rcbt2
         dxiaxid = 0.0d0
         dxiayid = 0.0d0
         dxiazid = 0.0d0
         dyiayia = rcbyt*dphidyia
         dyiazia = rcbyt*dphidzia - xcb*rcb/rt2
         dyiaxib = rcbyt*dphidxibt - dphidzt
     &                - (yca*zcb*ycb+zca*xzcb2)/rcbt2
         dyiayib = rcbyt*dphidyibt + ycb*(xca*zcb-zca*xcb)/rcbt2
         dyiazib = rcbyt*dphidzibt + dphidxt
     &                + (yca*xcb*ycb+xca*xzcb2)/rcbt2
         dyiaxic = rcbyt*dphidxict + dphidzt
     &                + (yba*zcb*ycb+zba*xzcb2)/rcbt2
         dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2
         dyiazic = rcbyt*dphidzict - dphidxt
     &                - (yba*xcb*ycb+xba*xzcb2)/rcbt2
         dyiaxid = 0.0d0
         dyiayid = 0.0d0
         dyiazid = 0.0d0
         dziazia = rcbzt*dphidzia
         dziaxib = rcbzt*dphidxibt + dphidyt
     &                + (zca*ycb*zcb+yca*xycb2)/rcbt2
         dziayib = rcbzt*dphidyibt - dphidxt
     &                - (zca*xcb*zcb+xca*xycb2)/rcbt2
         dziazib = rcbzt*dphidzibt + zcb*(yca*xcb-xca*ycb)/rcbt2
         dziaxic = rcbzt*dphidxict - dphidyt
     &                - (zba*ycb*zcb+yba*xycb2)/rcbt2
         dziayic = rcbzt*dphidyict + dphidxt
     &                + (zba*xcb*zcb+xba*xycb2)/rcbt2
         dziazic = rcbzt*dphidzict + zcb*zt/rcbt2
         dziaxid = 0.0d0
         dziayid = 0.0d0
         dziazid = 0.0d0
         dxibxic = -xcb*dphidxib/(rcb*rcb)
     &       - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &       - 2.0d0*(yt*zba-yba*zt)*dphidxibt/rt2
     &       - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &       + 2.0d0*(yu*zdb-ydb*zu)*dphidxibu/ru2
         dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
     &       - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &       - 2.0d0*(zt*xba-zba*xt)*dphidxibt/rt2
     &       + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &       + 2.0d0*(zu*xdb-zdb*xu)*dphidxibu/ru2
         dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2
         dxibyid = rcbyu*dphidxibu - dphidzu
     &                - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
         dxibzid = rcbzu*dphidxibu + dphidyu
     &                + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
         dyibzib = ycb*dphidzib/(rcb*rcb)
     &       - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &       - 2.0d0*(xt*zca-xca*zt)*dphidzibt/rt2
     &       + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &       + 2.0d0*(xu*zdc-xdc*zu)*dphidzibu/ru2
         dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu
     &       + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &       - 2.0d0*(yt*zba-yba*zt)*dphidyibt/rt2
     &       - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &       + 2.0d0*(yu*zdb-ydb*zu)*dphidyibu/ru2
         dyibyic = -ycb*dphidyib/(rcb*rcb)
     &       - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &       - 2.0d0*(zt*xba-zba*xt)*dphidyibt/rt2
     &       - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &       + 2.0d0*(zu*xdb-zdb*xu)*dphidyibu/ru2
         dyibxid = rcbxu*dphidyibu + dphidzu
     &                + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
         dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2
         dyibzid = rcbzu*dphidyibu - dphidxu
     &                - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
         dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu
     &       - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &       - 2.0d0*(yt*zba-yba*zt)*dphidzibt/rt2
     &       + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &       + 2.0d0*(yu*zdb-ydb*zu)*dphidzibu/ru2
         dzibzic = -zcb*dphidzib/(rcb*rcb)
     &       - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &       - 2.0d0*(xt*yba-xba*yt)*dphidzibt/rt2
     &       - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &       + 2.0d0*(xu*ydb-xdb*yu)*dphidzibu/ru2
         dzibxid = rcbxu*dphidzibu - dphidyu
     &                - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
         dzibyid = rcbyu*dphidzibu + dphidxu
     &                + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
         dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2
         dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
         dxicyid = rcbyu*dphidxicu + dphidzu
     &                + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
         dxiczid = rcbzu*dphidxicu - dphidyu
     &                - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
         dyicxid = rcbxu*dphidyicu - dphidzu
     &                - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
         dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
         dyiczid = rcbzu*dphidyicu + dphidxu
     &                + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
         dzicxid = rcbxu*dphidzicu + dphidyu
     &                + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
         dzicyid = rcbyu*dphidzicu - dphidxu
     &                - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
         dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
         dxidxid = rcbxu*dphidxid
         dxidyid = rcbxu*dphidyid + zcb*rcb/ru2
         dxidzid = rcbxu*dphidzid - ycb*rcb/ru2
         dyidyid = rcbyu*dphidyid
         dyidzid = rcbyu*dphidzid + xcb*rcb/ru2
         dzidzid = rcbzu*dphidzid
c
c     get some second derivative chain rule terms by difference
c
         dxibxib = -dxiaxib - dxibxic - dxibxid
         dxibyib = -dyiaxib - dxibyic - dxibyid
         dxibzib = -dxiazib - dzibxic - dzibxid
         dxibzic = -dziaxib - dxibzib - dxibzid
         dyibyib = -dyiayib - dyibyic - dyibyid
         dyibzic = -dziayib - dyibzib - dyibzid
         dzibzib = -dziazib - dzibzic - dzibzid
         dzibyic = -dyiazib - dyibzib - dzibyid
         dxicxic = -dxiaxic - dxibxic - dxicxid
         dxicyic = -dyiaxic - dyibxic - dxicyid
         dxiczic = -dziaxic - dzibxic - dxiczid
         dyicyic = -dyiayic - dyibyic - dyicyid
         dyiczic = -dziayic - dzibyic - dyiczid
         dziczic = -dziazic - dzibzic - dziczid
c
c     now, increment diagonal and off-diagonal Hessian elements
c
         if (i .eq. ia) then
            hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                        + d2edphi2*dphidxia*dphidxia
            hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                        + d2edphi2*dphidxia*dphidyia
            hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                        + d2edphi2*dphidxia*dphidzia
            hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                        + d2edphi2*dphidxia*dphidyia
            hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                        + d2edphi2*dphidyia*dphidyia
            hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                        + d2edphi2*dphidyia*dphidzia
            hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                        + d2edphi2*dphidxia*dphidzia
            hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                        + d2edphi2*dphidyia*dphidzia
            hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                        + d2edphi2*dphidzia*dphidzia
            hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                        + d2edphi2*dphidxia*dphidxib
            hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                        + d2edphi2*dphidyia*dphidxib
            hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                        + d2edphi2*dphidzia*dphidxib
            hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                        + d2edphi2*dphidxia*dphidyib
            hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                        + d2edphi2*dphidyia*dphidyib
            hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                        + d2edphi2*dphidzia*dphidyib
            hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                        + d2edphi2*dphidxia*dphidzib
            hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                        + d2edphi2*dphidyia*dphidzib
            hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                        + d2edphi2*dphidzia*dphidzib
            hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                        + d2edphi2*dphidxia*dphidxic
            hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                        + d2edphi2*dphidyia*dphidxic
            hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                        + d2edphi2*dphidzia*dphidxic
            hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                        + d2edphi2*dphidxia*dphidyic
            hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                        + d2edphi2*dphidyia*dphidyic
            hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                        + d2edphi2*dphidzia*dphidyic
            hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                        + d2edphi2*dphidxia*dphidzic
            hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                        + d2edphi2*dphidyia*dphidzic
            hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                        + d2edphi2*dphidzia*dphidzic
            hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                        + d2edphi2*dphidxia*dphidxid
            hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                        + d2edphi2*dphidyia*dphidxid
            hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                        + d2edphi2*dphidzia*dphidxid
            hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                        + d2edphi2*dphidxia*dphidyid
            hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                        + d2edphi2*dphidyia*dphidyid
            hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                        + d2edphi2*dphidzia*dphidyid
            hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                        + d2edphi2*dphidxia*dphidzid
            hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                        + d2edphi2*dphidyia*dphidzid
            hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                        + d2edphi2*dphidzia*dphidzid
         else if (i .eq. ib) then
            hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                        + d2edphi2*dphidxib*dphidxib
            hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                        + d2edphi2*dphidxib*dphidyib
            hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                        + d2edphi2*dphidxib*dphidzib
            hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                        + d2edphi2*dphidxib*dphidyib
            hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                        + d2edphi2*dphidyib*dphidyib
            hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                        + d2edphi2*dphidyib*dphidzib
            hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                        + d2edphi2*dphidxib*dphidzib
            hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                        + d2edphi2*dphidyib*dphidzib
            hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                        + d2edphi2*dphidzib*dphidzib
            hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                        + d2edphi2*dphidxib*dphidxia
            hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                        + d2edphi2*dphidyib*dphidxia
            hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                        + d2edphi2*dphidzib*dphidxia
            hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                        + d2edphi2*dphidxib*dphidyia
            hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                        + d2edphi2*dphidyib*dphidyia
            hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                        + d2edphi2*dphidzib*dphidyia
            hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                        + d2edphi2*dphidxib*dphidzia
            hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                        + d2edphi2*dphidyib*dphidzia
            hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                        + d2edphi2*dphidzib*dphidzia
            hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                        + d2edphi2*dphidxib*dphidxic
            hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                        + d2edphi2*dphidyib*dphidxic
            hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                        + d2edphi2*dphidzib*dphidxic
            hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                        + d2edphi2*dphidxib*dphidyic
            hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                        + d2edphi2*dphidyib*dphidyic
            hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                        + d2edphi2*dphidzib*dphidyic
            hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                        + d2edphi2*dphidxib*dphidzic
            hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                        + d2edphi2*dphidyib*dphidzic
            hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                        + d2edphi2*dphidzib*dphidzic
            hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                        + d2edphi2*dphidxib*dphidxid
            hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                        + d2edphi2*dphidyib*dphidxid
            hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                        + d2edphi2*dphidzib*dphidxid
            hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                        + d2edphi2*dphidxib*dphidyid
            hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                        + d2edphi2*dphidyib*dphidyid
            hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                        + d2edphi2*dphidzib*dphidyid
            hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                        + d2edphi2*dphidxib*dphidzid
            hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                        + d2edphi2*dphidyib*dphidzid
            hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                        + d2edphi2*dphidzib*dphidzid
         else if (i .eq. ic) then
            hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                        + d2edphi2*dphidxic*dphidxic
            hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                        + d2edphi2*dphidxic*dphidyic
            hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                        + d2edphi2*dphidxic*dphidzic
            hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                        + d2edphi2*dphidxic*dphidyic
            hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                        + d2edphi2*dphidyic*dphidyic
            hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                        + d2edphi2*dphidyic*dphidzic
            hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                        + d2edphi2*dphidxic*dphidzic
            hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                        + d2edphi2*dphidyic*dphidzic
            hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                        + d2edphi2*dphidzic*dphidzic
            hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                        + d2edphi2*dphidxic*dphidxia
            hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                        + d2edphi2*dphidyic*dphidxia
            hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                        + d2edphi2*dphidzic*dphidxia
            hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                        + d2edphi2*dphidxic*dphidyia
            hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                        + d2edphi2*dphidyic*dphidyia
            hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                        + d2edphi2*dphidzic*dphidyia
            hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                        + d2edphi2*dphidxic*dphidzia
            hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                        + d2edphi2*dphidyic*dphidzia
            hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                        + d2edphi2*dphidzic*dphidzia
            hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                        + d2edphi2*dphidxic*dphidxib
            hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                        + d2edphi2*dphidyic*dphidxib
            hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                        + d2edphi2*dphidzic*dphidxib
            hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                        + d2edphi2*dphidxic*dphidyib
            hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                        + d2edphi2*dphidyic*dphidyib
            hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                        + d2edphi2*dphidzic*dphidyib
            hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                        + d2edphi2*dphidxic*dphidzib
            hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                        + d2edphi2*dphidyic*dphidzib
            hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                        + d2edphi2*dphidzic*dphidzib
            hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                        + d2edphi2*dphidxic*dphidxid
            hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                        + d2edphi2*dphidyic*dphidxid
            hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                        + d2edphi2*dphidzic*dphidxid
            hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                        + d2edphi2*dphidxic*dphidyid
            hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                        + d2edphi2*dphidyic*dphidyid
            hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                        + d2edphi2*dphidzic*dphidyid
            hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                        + d2edphi2*dphidxic*dphidzid
            hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                        + d2edphi2*dphidyic*dphidzid
            hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                        + d2edphi2*dphidzic*dphidzid
         else if (i .eq. id) then
            hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                        + d2edphi2*dphidxid*dphidxid
            hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                        + d2edphi2*dphidxid*dphidyid
            hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                        + d2edphi2*dphidxid*dphidzid
            hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                        + d2edphi2*dphidxid*dphidyid
            hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                        + d2edphi2*dphidyid*dphidyid
            hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                        + d2edphi2*dphidyid*dphidzid
            hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                        + d2edphi2*dphidxid*dphidzid
            hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                        + d2edphi2*dphidyid*dphidzid
            hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                        + d2edphi2*dphidzid*dphidzid
            hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                        + d2edphi2*dphidxid*dphidxia
            hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                        + d2edphi2*dphidyid*dphidxia
            hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                        + d2edphi2*dphidzid*dphidxia
            hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                        + d2edphi2*dphidxid*dphidyia
            hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                        + d2edphi2*dphidyid*dphidyia
            hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                        + d2edphi2*dphidzid*dphidyia
            hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                        + d2edphi2*dphidxid*dphidzia
            hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                        + d2edphi2*dphidyid*dphidzia
            hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                        + d2edphi2*dphidzid*dphidzia
            hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                        + d2edphi2*dphidxid*dphidxib
            hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                        + d2edphi2*dphidyid*dphidxib
            hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                        + d2edphi2*dphidzid*dphidxib
            hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                        + d2edphi2*dphidxid*dphidyib
            hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                        + d2edphi2*dphidyid*dphidyib
            hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                        + d2edphi2*dphidzid*dphidyib
            hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                        + d2edphi2*dphidxid*dphidzib
            hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                        + d2edphi2*dphidyid*dphidzib
            hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                        + d2edphi2*dphidzid*dphidzib
            hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                        + d2edphi2*dphidxid*dphidxic
            hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                        + d2edphi2*dphidyid*dphidxic
            hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                        + d2edphi2*dphidzid*dphidxic
            hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                        + d2edphi2*dphidxid*dphidyic
            hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                        + d2edphi2*dphidyid*dphidyic
            hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                        + d2edphi2*dphidzid*dphidyic
            hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                        + d2edphi2*dphidxid*dphidzic
            hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                        + d2edphi2*dphidyid*dphidzic
            hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                        + d2edphi2*dphidzid*dphidzic
         end if
   10    continue
      end do
c
c     compute Hessian elements for shallow Gaussian basin restraint
c
      if (use_basin) then
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = 1, n
            if (k .ne. i) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               rik2 = xr*xr + yr*yr + zr*zr
               expterm = depth * width * exp(-width*rik2)
               dedr = -2.0d0 * expterm
               d2edr2 = (4.0d0*width*rik2-2.0d0) * expterm
c
c     set the chain rule terms for the Hessian elements
c
               if (rik2 .eq. 0.0d0) then
                  term = 0.0d0
               else
                  term = (d2edr2-dedr) / rik2
               end if
               termx = term * xr
               termy = term * yr
               termz = term * zr
               d2e(1,1) = termx*xr + dedr
               d2e(1,2) = termx*yr
               d2e(1,3) = termx*zr
               d2e(2,1) = d2e(1,2)
               d2e(2,2) = termy*yr + dedr
               d2e(2,3) = termy*zr
               d2e(3,1) = d2e(1,3)
               d2e(3,2) = d2e(2,3)
               d2e(3,3) = termz*zr + dedr
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + d2e(1,j)
                  hessy(j,i) = hessy(j,i) + d2e(2,j)
                  hessz(j,i) = hessz(j,i) + d2e(3,j)
                  hessx(j,k) = hessx(j,k) - d2e(1,j)
                  hessy(j,k) = hessy(j,k) - d2e(2,j)
                  hessz(j,k) = hessz(j,k) - d2e(3,j)
               end do
            end if
         end do
      end if
c
c     compute Hessian elements for a spherical boundary restraint
c
      if (use_wall) then
         buffer = 2.5d0
         a = 2048.0d0
         b = 64.0d0
         xi = x(i)
         yi = y(i)
         zi = z(i)
         ri2 = xi**2 + yi**2 + zi**2
         ri = sqrt(ri2)
         r = rwall + buffer - ri
         r2 = r * r
         r6 = r2 * r2 * r2
         r12 = r6 * r6
         if (ri .eq. 0.0d0) then
            ri = 1.0d0
            ri2 = 1.0d0
         end if
         dedr = (12.0d0*a/r12 - 6.0d0*b/r6) / (r*ri)
         d2edr2 = (156.0d0*a/r12 - 42.0d0*b/r6) / (r2*ri2)
c
c     set the chain rule terms for the Hessian elements
c
         d2edr2 = d2edr2 - dedr/ri2
         termx = d2edr2 * xi
         termy = d2edr2 * yi
         termz = d2edr2 * zi
         d2e(1,1) = termx*xi + dedr
         d2e(1,2) = termx*yi
         d2e(1,3) = termx*zi
         d2e(2,1) = d2e(1,2)
         d2e(2,2) = termy*yi + dedr
         d2e(2,3) = termy*zi
         d2e(3,1) = d2e(1,3)
         d2e(3,2) = d2e(2,3)
         d2e(3,3) = termz*zi + dedr
c
c     increment diagonal and non-diagonal Hessian elements
c
         do j = 1, 3
            hessx(j,i) = hessx(j,i) + d2e(1,j)
            hessy(j,i) = hessy(j,i) + d2e(2,j)
            hessz(j,i) = hessz(j,i) + d2e(3,j)
         end do
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
c     ##  subroutine egeom3  --  restraint energy terms & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egeom3" calculates the energy due to restraints on atomic
c     positions, interatomic distances, dihedral angles and Gaussian
c     weighted molecular size; also partitions energy among the atoms
c
c
      subroutine egeom3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'energi.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'restrn.i'
      include 'usage.i'
      integer i,j,k,ia,ib,ic,id
      real*8 e,xr,yr,zr,rik,rik2,dt,dt2
      real*8 angle,target,force,cosine,sine
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 df1,df2,tf1,tf2,t1,t2
      real*8 xi,yi,zi,ri
      real*8 a,b,buffer
      real*8 r,r2,r6,r12
      logical header,huge,intermol
c
c
c     zero out the restraint energy and partitioning terms
c
      neg = 0
      eg = 0.0d0
      do i = 1, n
         aeg(i) = 0.0d0
      end do
c
c     compute the pseudoenergy for position restraints
c
      header = .true.
      do j = 1, npfix
         i = ipfix(j)
         if (use(i)) then
            xr = x(i) - xpfix(j)
            yr = y(i) - ypfix(j)
            zr = z(i) - zpfix(j)
            dt2 = xr*xr + yr*yr + zr*zr
            force = pfix(j)
            e = force * dt2
            neg = neg + 1
            eg = eg + e
            aeg(i) = aeg(i) + e
            huge = (e .gt. 10.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Atomic Position Restraint',
     &                       ' Terms :',
     &                    //,' Type',6x,'Atom Name',11x,'Target',
     &                       ' Position',7x,'Distance',6x,'Energy',/)
               end if
               dt = sqrt(dt2)
               write (iout,20)  i,name(i),xpfix(j),ypfix(j),
     &                          zpfix(j),dt,e
   20          format (' Position ',i5,'-',a3,2x,4f10.4,f12.4)
            end if
         end if
      end do
c
c     compute the pseudoenergy for distance restraints
c
      header = .true.
      do j = 1, ndfix
         i = idfix(1,j)
         k = idfix(2,j)
         if (use(i) .or. use(k)) then
            xr = x(i) - x(k)
            yr = y(i) - y(k)
            zr = z(i) - z(k)
            intermol = (molcule(i) .ne. molcule(k))
            if (use_bounds .and. intermol)  call image (xr,yr,zr,0)
            rik2 = xr*xr + yr*yr + zr*zr
            rik = sqrt(rik2)
            df1 = dfix(1,j)
            df2 = dfix(2,j)
            target = rik
            if (rik .lt. df1)  target = df1
            if (rik .gt. df2)  target = df2
            force = dfix(3,j)
            dt = rik - target
            dt2 = dt * dt
            e = force * dt2
            neg = neg + 1
            eg = eg + e
            aeg(i) = aeg(i) + 0.5d0*e
            aeg(k) = aeg(k) + 0.5d0*e
            huge = (e .gt. 10.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Individual Interatomic Distance',
     &                       ' Restraint Terms :',
     &                    //,' Type',11x,'Atom Names',14x,'Ideal Range',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,40)  i,name(i),k,name(k),df1,df2,rik,e
   40          format (' Distance ',i5,'-',a3,1x,i5,'-',a3,
     &                    6x,2f8.2,f10.4,f12.4)
            end if
         end if
      end do
c
c     compute the pseudoenergy for dihedral angle restraints
c
      header = .true.
      do i = 1, ntfix
         ia = itfix(1,i)
         ib = itfix(2,i)
         ic = itfix(3,i)
         id = itfix(4,i)
         if (use(ia) .or. use(ib) .or. use(ic) .or. use(id)) then
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
               tf1 = tfix(1,i)
               tf2 = tfix(2,i)
               if (angle.gt.tf1 .and. angle.lt.tf2) then
                  target = angle
               else if (angle.gt.tf1 .and. tf1.gt.tf2) then
                  target = angle
               else if (angle.lt.tf2 .and. tf1.gt.tf2) then
                  target = angle
               else
                  t1 = angle - tf1
                  t2 = angle - tf2
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     target = tf1
                  else
                     target = tf2
                  end if
               end if
               force = tfix(3,i)
               dt = angle - target
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
               dt2 = dt * dt
               e = force * dt2
               neg = neg + 1
               eg = eg + e
               aeg(ib) = aeg(ib) + 0.5d0*e
               aeg(ic) = aeg(ic) + 0.5d0*e
               huge = (e .gt. 10.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,50)
   50                format (/,' Individual Dihedral Angle Restraint',
     &                          ' Terms :',
     &                       //,' Type',10x,'Atom Numbers',13x,'Ideal',
     &                          ' Range',4x,'Actual',6x,'Energy',/)
                  end if
                  write (iout,60)  ia,ib,ic,id,tf1,tf2,angle,e
   60             format (' Dihedral ',4i5,5x,2f8.2,f10.4,f12.4)
               end if
            end if
         end if
      end do
c
c     compute the energy for a shallow Gaussian basin restraint
c
      if (use_basin) then
         header = .true.
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do k = i+1, n
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               rik2 = xr*xr + yr*yr + zr*zr
               e = depth*exp(-width*rik2) - depth
               neg = neg + 1
               eg = eg + e
               aeg(i) = aeg(i) + 0.5d0*e
               aeg(k) = aeg(k) + 0.5d0*e
               huge = (e .gt. 10.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,70)
   70                format (/,' Individual Gaussian Basin',
     &                          ' Restraint Terms :',
     &                       //,' Type',11x,'Atom Names',20x,'Ideal',
     &                          4x,'Actual',6x,'Energy',/)
                  end if
                  rik = sqrt(rik2)
                  write (iout,80)  i,name(i),k,name(k),0.0d0,rik,e
   80             format (' Distance ',i5,'-',a3,1x,i5,'-',a3,
     &                       12x,2f10.4,f12.4)
               end if
            end do
         end do
      end if
c
c     compute the energy for a spherical droplet boundary restraint
c
      if (use_wall) then
         header = .true.
         buffer = 2.5d0
         a = 2048.0d0
         b = 64.0d0
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = sqrt(xi**2 + yi**2 + zi**2)
            r = rwall + buffer - ri
            r2 = r * r
            r6 = r2 * r2 * r2
            r12 = r6 * r6
            e = a/r12 - b/r6
            neg = neg + 1
            eg = eg + e
            aeg(i) = aeg(i) + e
            huge = (e .gt. 10.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,90)
   90             format (/,' Individual Spherical Boundary',
     &                       ' Restraint Terms :',
     &                    //,' Type',11x,'Atom Name',28x,'Distance',
     &                       6x,'Energy',/)
               end if
               write (iout,100)  i,name(i),ri,e
  100          format (' Wall',10x,i5,'-',a3,27x,f10.4,f12.4)
            end if
         end do
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
c     ##  subroutine ehal  --  Buffered 14-7 van der Waals energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ehal" calculates the van der Waals interaction energy using
c     Halgren's Buffered 14-7 formula
c
c
      subroutine ehal
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
      real*8 e,eps,rdn,fgrp
      real*8 rv,rv7,taper
      real*8 xi,yi,zi,xr,yr,zr
      real*8 sigma,kappa,kappa7,rho,tau
      real*8 rik,rik2,rik3,rik4,rik5,rik7
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
c     set the values of the two buffering coefficients
c
      sigma = 1.12d0
      kappa = 1.07d0
      kappa7 = kappa**7
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
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i)  eps = eps / vdwscale
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + (sigma-1.0d0)*rv7
                  tau = rik + (kappa-1.0d0)*rv
                  e = eps * kappa7 * (rv7/tau**7)
     &                     * (sigma*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
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
                     rik = sqrt(rik2)
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + (sigma-1.0d0)*rv7
                     tau = rik + (kappa-1.0d0)*rv
                     e = eps * kappa7 * (rv7/tau**7)
     &                        * (sigma*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
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
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal1  --  Buffered 14-7 energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal1" calculates the van der Waals interaction energy and
c     its first derivatives with respect to Cartesian coordinates
c     using Halgren's Buffered 14-7 formula
c
c
      subroutine ehal1
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
      real*8 e,eps,rdn,fgrp
      real*8 rv,rv7,rv14
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 de,dedx,dedy,dedz,taper,dtaper
      real*8 sigma,kappa,kappa7,rho,tau
      real*8 rik,rik2,rik3,rik4,rik5,rik6,rik7
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
c     set the values of the two buffering coefficients
c
      sigma = 1.12d0
      kappa = 1.07d0
      kappa7 = kappa**7
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
                  rik = sqrt(rik2)
                  rv7 = rv**7
                  rv14 = rv7 * rv7
                  rik6 = rik2**3
                  rik7 = rik6 * rik
                  rho = rik7 + (sigma-1.0d0)*rv7
                  tau = rik + (kappa-1.0d0)*rv
                  e = eps * kappa7 * (rv7/tau**7)
     &                     * (sigma*rv7/rho-2.0d0)
                  de = -7.0d0 * eps * kappa7 * (rv7/tau**8)
     &                      * (sigma*rv7/rho-2.0d0)
     &                 -7.0d0 * eps * kappa7 * sigma * rv14
     &                      * rik6 / (rho**2 * tau**7)
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
                     rik = sqrt(rik2)
                     rv7 = rv**7
                     rv14 = rv7 * rv7
                     rik6 = rik2**3
                     rik7 = rik6 * rik
                     rho = rik7 + (sigma-1.0d0)*rv7
                     tau = rik + (kappa-1.0d0)*rv
                     e = eps * kappa7 * (rv7/tau**7)
     &                        * (sigma*rv7/rho-2.0d0)
                     de = -7.0d0 * eps * kappa7 * (rv7/tau**8)
     &                         * (sigma*rv7/rho-2.0d0)
     &                    -7.0d0 * eps * kappa7 * sigma * rv14
     &                         * rik6 / (rho**2 * tau**7)
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
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal2  --  atom-by-atom Buffered 14-7 Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal2" calculates the van der Waals second derivatives for a
c     single atom at a time using Halgren's Buffered 14-7 formula
c
c
      subroutine ehal2 (iatom,xred,yred,zred)
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
      real*8 eps,rv,rv7,rv14
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 redi2,rediv2,rediiv
      real*8 redik,redivk,redikv,redivkv
      real*8 rik,rik2,rik3,rik4
      real*8 rik5,rik6,rik7,rik12
      real*8 taper,dtaper,d2taper
      real*8 d2edx,d2edy,d2edz,term(3,3)
      real*8 sigma,kappa,kappa7,rho,tau
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
c     set the values of the two buffering coefficients
c
      sigma = 1.12d0
      kappa = 1.07d0
      kappa7 = kappa**7
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
                  rik = sqrt(rik2)
                  rv7 = rv**7
                  rv14 = rv7 * rv7
                  rik6 = rik2**3
                  rik5 = rik6 / rik
                  rik7 = rik6 * rik
                  rik12 = rik6 * rik6
                  rho = rik7 + (sigma-1.0d0)*rv7
                  tau = rik + (kappa-1.0d0)*rv
                  de = -7.0d0 * eps * kappa7 * (rv7/tau**8)
     &                           * (sigma*rv7/rho-2.0d0)
     &                 -7.0d0 * eps * kappa7 * sigma * rv14
     &                           * rik6 / (rho**2 * tau**7)
                  d2e = 56.0d0 * eps * kappa7 * (rv7/tau**9)
     &                           * (sigma*rv7/rho-2.0d0)
     &                  + 98.0d0 * eps * kappa7 * sigma * rv14
     &                           * rik6 / (rho**2 * tau**8)
     &                  + 98.0d0 * eps * kappa7 * sigma * rv14
     &                           * rik12 / (rho**3 * tau**7)
     &                  - 42.0d0 * eps * kappa7 * sigma * rv14
     &                           * rik5 / (rho**2 * tau**7)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     e = eps * kappa7 * (rv7/tau**7)
     &                      * (sigma*rv7/rho-2.0d0)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
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
                     rik = sqrt(rik2)
                     rv7 = rv**7
                     rv14 = rv7 * rv7
                     rik6 = rik2**3
                     rik5 = rik6 / rik
                     rik7 = rik6 * rik
                     rik12 = rik6 * rik6
                     rho = rik7 + (sigma-1.0d0)*rv7
                     tau = rik + (kappa-1.0d0)*rv
                     de = -7.0d0 * eps * kappa7 * (rv7/tau**8)
     &                              * (sigma*rv7/rho-2.0d0)
     &                    -7.0d0 * eps * kappa7 * sigma * rv14
     &                              * rik6 / (rho**2 * tau**7)
                     d2e = 56.0d0 * eps * kappa7 * (rv7/tau**9)
     &                              * (sigma*rv7/rho-2.0d0)
     &                     + 98.0d0 * eps * kappa7 * sigma * rv14
     &                              * rik6 / (rho**2 * tau**8)
     &                     + 98.0d0 * eps * kappa7 * sigma * rv14
     &                              * rik12 / (rho**3 * tau**7)
     &                     - 42.0d0 * eps * kappa7 * sigma * rv14
     &                              * rik5 / (rho**2 * tau**7)
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        e = eps * kappa7 * (rv7/tau**7)
     &                         * ((sigma*rv7/rho)-2.0d0)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
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
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ehal3  --  Buffered 14-7 energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ehal3" calculates the van der Waals interaction energy using
c     Halgren's Buffered 14-7 formula and also partitions the energy
c     among the atoms
c
c
      subroutine ehal3
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
      real*8 e,eps,rdn,fgrp
      real*8 rv,rv7,taper
      real*8 xi,yi,zi,xr,yr,zr
      real*8 sigma,kappa,kappa7,rho,tau
      real*8 rik,rik2,rik3,rik4,rik5,rik7
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
c     set the values of the two buffering coefficients
c
      sigma = 1.12d0
      kappa = 1.07d0
      kappa7 = kappa**7
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
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i)  eps = eps / vdwscale
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + (sigma-1.0d0)*rv7
                  tau = rik + (kappa-1.0d0)*rv
                  e = eps * kappa7 * (rv7/tau**7)
     &                   * (sigma*rv7/rho-2.0d0)
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
   20                format (' VDW-MMFF ',i5,'-',a3,1x,i5,'-',a3,
     &                         12x,2f10.4,f12.4)
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
                     rik = sqrt(rik2)
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + (sigma-1.0d0)*rv7
                     tau = rik + (kappa-1.0d0)*rv
                     e = eps * kappa7 * (rv7/tau**7)
     &                      * (sigma*rv7/rho-2.0d0)
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
   40                   format (' VDW-MMFF ',i5,'-',a3,1x,i5,'-',a3,
     &                            '   (XTAL)   ',2f10.4,f12.4)
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal4  --  Buffered 14-7 van der Waals energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal4" calculates the van der Waals interaction energy using
c     Halgren's Buffered 14-7 formula and the method of lights to
c     locate neighboring atoms
c
c
      subroutine ehal4
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
      real*8 e,rv,rv7,eps,rdn,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 rik,rik2,rik3,rik4,rik5,rik7,taper
      real*8 sigma,kappa,kappa7,rho,tau
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
c     set the values of the two buffering coefficients
c
      sigma = 1.12d0
      kappa = 1.07d0
      kappa7 = kappa**7
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
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
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (kskip .eq. -i)  eps = eps / vdwscale
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + (sigma-1.0d0)*rv7
                  tau = rik + (kappa-1.0d0)*rv
                  e = eps * kappa7 * (rv7/tau**7)
     &                     * (sigma*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal5  --  Buffered 14-7 energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal5" calculates the van der Waals interaction energy and
c     its first derivatives using Halgren's Buffered 14-7 formula
c     and the method of lights to locate neighboring atoms
c
c
      subroutine ehal5
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
      real*8 e,rv,rv7,rv14,eps,rdn,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 dedx,dedy,dedz,de,taper,dtaper
      real*8 rik,rik2,rik3,rik4,rik5,rik6,rik7
      real*8 sigma,kappa,kappa7,rho,tau
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
c     set the values of the two buffering coefficients
c
      sigma = 1.12d0
      kappa = 1.07d0
      kappa7 = kappa**7
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
                  rik = sqrt(rik2)
                  rv7 = rv**7
                  rv14 = rv7 * rv7
                  rik6 = rik2**3
                  rik7 = rik6 * rik
                  rho = rik7 + (sigma-1.0d0)*rv7
                  tau = rik + (kappa-1.0d0)*rv
                  e = eps * kappa7 * (rv7/tau**7)
     &                     * (sigma*rv7/rho-2.0d0)
                  de = -7.0d0 * eps * kappa7 * (rv7/tau**8)
     &                      * (sigma*rv7/rho-2.0d0)
     &                 -7.0d0 * eps * kappa7 * sigma * rv14 * rik6
     &                      / (rho**2 * tau**7)
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
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine eimprop  --  improper dihedral energy  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "eimprop" calculates the improper dihedral potential energy
c
c
      subroutine eimprop
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'improp.i'
      include 'math.i'
      include 'usage.i'
      integer i,ia,ib,ic,id
      real*8 e,ideal,force
      real*8 dt,cosine,sine
      real*8 angle,fgrp
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb,rcb
      real*8 xdc,ydc,zdc
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      logical proceed
c
c
c     zero out improper dihedral energy
c
      eid = 0.0d0
c
c     calculate the improper dihedral angle energy term
c
      do i = 1, niprop
         ia = iiprop(1,i)
         ib = iiprop(2,i)
         ic = iiprop(3,i)
         id = iiprop(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the improper dihedral angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     set the improper dihedral parameters for this angle
c
               ideal = vprop(i)
               force = kprop(i)
               dt = angle - ideal
               dowhile (dt .gt. 180.0d0)
                  dt = dt - 360.0d0
               end do
               dowhile (dt .lt. -180.0d0)
                  dt = dt + 360.0d0
               end do
               dt = dt / radian
c
c     calculate the improper dihedral energy
c
               e = force * dt**2
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total improper dihedral energy
c
               eid = eid + e
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
c     ##  subroutine eimprop1  --  impr. dihedral energy & gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eimprop1" calculates improper dihedral energy and its
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine eimprop1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'improp.i'
      include 'math.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ia,ib,ic,id
      real*8 e,dedphi
      real*8 ideal,force
      real*8 dt,cosine,sine
      real*8 angle,fgrp
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb,rcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      logical proceed
c
c
c     zero out energy and first derivative components
c
      eid = 0.0d0
      do i = 1, n
         deid(1,i) = 0.0d0
         deid(2,i) = 0.0d0
         deid(3,i) = 0.0d0
      end do
c
c     calculate the improper dihedral angle energy term
c
      do i = 1, niprop
         ia = iiprop(1,i)
         ib = iiprop(2,i)
         ic = iiprop(3,i)
         id = iiprop(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the improper dihedral angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     set the improper dihedral parameters for this angle
c
               ideal = vprop(i)
               force = kprop(i)
               dt = angle - ideal
               dowhile (dt .gt. 180.0d0)
                  dt = dt - 360.0d0
               end do
               dowhile (dt .lt. -180.0d0)
                  dt = dt + 360.0d0
               end do
               dt = dt / radian
c
c     calculate improper energy and master chain rule term
c
               e = force * dt**2
               dedphi = 2.0d0 * force * dt
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedphi = dedphi * fgrp
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = yca*dphidzt - zca*dphidyt
     &                       + zdc*dphidyu - ydc*dphidzu
               dphidyib = zca*dphidxt - xca*dphidzt
     &                       + xdc*dphidzu - zdc*dphidxu
               dphidzib = xca*dphidyt - yca*dphidxt
     &                       + ydc*dphidxu - xdc*dphidyu
               dphidxic = zba*dphidyt - yba*dphidzt
     &                       + ydb*dphidzu - zdb*dphidyu
               dphidyic = xba*dphidzt - zba*dphidxt
     &                       + zdb*dphidxu - xdb*dphidzu
               dphidzic = yba*dphidxt - xba*dphidyt
     &                       + xdb*dphidyu - ydb*dphidxu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     compute first derivative components for this angle
c
               dedxia = dedphi * dphidxia
               dedyia = dedphi * dphidyia
               dedzia = dedphi * dphidzia
               dedxib = dedphi * dphidxib
               dedyib = dedphi * dphidyib
               dedzib = dedphi * dphidzib
               dedxic = dedphi * dphidxic
               dedyic = dedphi * dphidyic
               dedzic = dedphi * dphidzic
               dedxid = dedphi * dphidxid
               dedyid = dedphi * dphidyid
               dedzid = dedphi * dphidzid
c
c     calculate improper dihedral energy and derivatives
c
               eid = eid + e
               deid(1,ia) = deid(1,ia) + dedxia
               deid(2,ia) = deid(2,ia) + dedyia
               deid(3,ia) = deid(3,ia) + dedzia
               deid(1,ib) = deid(1,ib) + dedxib
               deid(2,ib) = deid(2,ib) + dedyib
               deid(3,ib) = deid(3,ib) + dedzib
               deid(1,ic) = deid(1,ic) + dedxic
               deid(2,ic) = deid(2,ic) + dedyic
               deid(3,ic) = deid(3,ic) + dedzic
               deid(1,id) = deid(1,id) + dedxid
               deid(2,id) = deid(2,id) + dedyid
               deid(3,id) = deid(3,id) + dedzid
c
c     increment the virial for use in pressure computation
c
               if (isobaric) then
                  virx = virx - xba*dedxia + xdc*dedxid
     &                      + xcb*(dedxic+dedxid)
                  viry = viry - yba*dedyia + ydc*dedyid
     &                      + ycb*(dedyic+dedyid)
                  virz = virz - zba*dedzia + zdc*dedzid
     &                      + zcb*(dedzic+dedzid)
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eimprop2  --  atom-wise imp. dihedral Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eimprop2" calculates second derivatives of the improper
c     dihedral angle energy for a single atom
c
c
      subroutine eimprop2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      include 'hessn.i'
      include 'improp.i'
      include 'math.i'
      integer i,ia,ib,ic,id,kiprop
      real*8 angle,ideal,force
      real*8 dedphi,d2edphi2
      real*8 dt,cosine,sine,fgrp
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 dphidxibt,dphidyibt,dphidzibt
      real*8 dphidxibu,dphidyibu,dphidzibu
      real*8 dphidxict,dphidyict,dphidzict
      real*8 dphidxicu,dphidyicu,dphidzicu
      real*8 dxiaxia,dyiayia,dziazia,dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic,dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia,dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic,dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib,dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib,dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic,dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid,dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid,dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic,dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid,dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid,dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid,dzicxid,dzicyid,dziczid
      logical proceed
c
c
c     calculate the improper dihedral angle energy term
c
      do kiprop = 1, niprop
         ia = iiprop(1,kiprop)
         ib = iiprop(2,kiprop)
         ic = iiprop(3,kiprop)
         id = iiprop(4,kiprop)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (i.eq.ia .or. i.eq.ib .or.
     &                              i.eq.ic .or. i.eq.id)
c
c     compute the value of the improper dihedral angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     set the improper dihedral parameters for this angle
c
               ideal = vprop(kiprop)
               force = kprop(kiprop)
               dt = angle - ideal
               dowhile (dt .gt. 180.0d0)
                  dt = dt - 360.0d0
               end do
               dowhile (dt .lt. -180.0d0)
                  dt = dt + 360.0d0
               end do
               dt = dt / radian
c
c     calculate the improper torsion master chain rule terms
c
               dedphi = 2.0d0 * force * dt
               d2edphi2 = 2.0d0 * force
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
               xycb2 = xcb*xcb + ycb*ycb
               xzcb2 = xcb*xcb + zcb*zcb
               yzcb2 = ycb*ycb + zcb*zcb
               rcbxt = -2.0d0 * rcb * dphidxt
               rcbyt = -2.0d0 * rcb * dphidyt
               rcbzt = -2.0d0 * rcb * dphidzt
               rcbt2 = rcb * rt2
               rcbxu = 2.0d0 * rcb * dphidxu
               rcbyu = 2.0d0 * rcb * dphidyu
               rcbzu = 2.0d0 * rcb * dphidzu
               rcbu2 = rcb * ru2
               dphidxibt = yca*dphidzt - zca*dphidyt
               dphidxibu = zdc*dphidyu - ydc*dphidzu
               dphidyibt = zca*dphidxt - xca*dphidzt
               dphidyibu = xdc*dphidzu - zdc*dphidxu
               dphidzibt = xca*dphidyt - yca*dphidxt
               dphidzibu = ydc*dphidxu - xdc*dphidyu
               dphidxict = zba*dphidyt - yba*dphidzt
               dphidxicu = ydb*dphidzu - zdb*dphidyu
               dphidyict = xba*dphidzt - zba*dphidxt
               dphidyicu = zdb*dphidxu - xdb*dphidzu
               dphidzict = yba*dphidxt - xba*dphidyt
               dphidzicu = xdb*dphidyu - ydb*dphidxu
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = dphidxibt + dphidxibu
               dphidyib = dphidyibt + dphidyibu
               dphidzib = dphidzibt + dphidzibu
               dphidxic = dphidxict + dphidxicu
               dphidyic = dphidyict + dphidyicu
               dphidzic = dphidzict + dphidzicu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     chain rule terms for second derivative components
c
               dxiaxia = rcbxt*dphidxia
               dxiayia = rcbxt*dphidyia - zcb*rcb/rt2
               dxiazia = rcbxt*dphidzia + ycb*rcb/rt2
               dxiaxib = rcbxt*dphidxibt + xcb*(zca*ycb-yca*zcb)/rcbt2
               dxiayib = rcbxt*dphidyibt + dphidzt
     &                      + (xca*zcb*xcb+zca*yzcb2)/rcbt2
               dxiazib = rcbxt*dphidzibt - dphidyt
     &                      - (xca*ycb*xcb+yca*yzcb2)/rcbt2
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2
               dxiayic = rcbxt*dphidyict - dphidzt
     &                      - (xba*zcb*xcb+zba*yzcb2)/rcbt2
               dxiazic = rcbxt*dphidzict + dphidyt
     &                      + (xba*ycb*xcb+yba*yzcb2)/rcbt2
               dxiaxid = 0.0d0
               dxiayid = 0.0d0
               dxiazid = 0.0d0
               dyiayia = rcbyt*dphidyia
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2
               dyiaxib = rcbyt*dphidxibt - dphidzt
     &                      - (yca*zcb*ycb+zca*xzcb2)/rcbt2
               dyiayib = rcbyt*dphidyibt + ycb*(xca*zcb-zca*xcb)/rcbt2
               dyiazib = rcbyt*dphidzibt + dphidxt
     &                      + (yca*xcb*ycb+xca*xzcb2)/rcbt2
               dyiaxic = rcbyt*dphidxict + dphidzt
     &                      + (yba*zcb*ycb+zba*xzcb2)/rcbt2
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2
               dyiazic = rcbyt*dphidzict - dphidxt
     &                      - (yba*xcb*ycb+xba*xzcb2)/rcbt2
               dyiaxid = 0.0d0
               dyiayid = 0.0d0
               dyiazid = 0.0d0
               dziazia = rcbzt*dphidzia
               dziaxib = rcbzt*dphidxibt + dphidyt
     &                      + (zca*ycb*zcb+yca*xycb2)/rcbt2
               dziayib = rcbzt*dphidyibt - dphidxt
     &                      - (zca*xcb*zcb+xca*xycb2)/rcbt2
               dziazib = rcbzt*dphidzibt + zcb*(yca*xcb-xca*ycb)/rcbt2
               dziaxic = rcbzt*dphidxict - dphidyt
     &                      - (zba*ycb*zcb+yba*xycb2)/rcbt2
               dziayic = rcbzt*dphidyict + dphidxt
     &                      + (zba*xcb*zcb+xba*xycb2)/rcbt2
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2
               dziaxid = 0.0d0
               dziayid = 0.0d0
               dziazid = 0.0d0
               dxibxic = -xcb*dphidxib/(rcb*rcb)
     &             - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidxibt/rt2
     &             - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidxibu/ru2
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
     &             - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidxibt/rt2
     &             + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidxibu/ru2
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2
               dxibyid = rcbyu*dphidxibu - dphidzu
     &                      - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
               dxibzid = rcbzu*dphidxibu + dphidyu
     &                      + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
               dyibzib = ycb*dphidzib/(rcb*rcb)
     &             - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &             - 2.0d0*(xt*zca-xca*zt)*dphidzibt/rt2
     &             + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &             + 2.0d0*(xu*zdc-xdc*zu)*dphidzibu/ru2
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu
     &             + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidyibt/rt2
     &             - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidyibu/ru2
               dyibyic = -ycb*dphidyib/(rcb*rcb)
     &             - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidyibt/rt2
     &             - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidyibu/ru2
               dyibxid = rcbxu*dphidyibu + dphidzu
     &                      + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2
               dyibzid = rcbzu*dphidyibu - dphidxu
     &                      - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu
     &             - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidzibt/rt2
     &             + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidzibu/ru2
               dzibzic = -zcb*dphidzib/(rcb*rcb)
     &             - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &             - 2.0d0*(xt*yba-xba*yt)*dphidzibt/rt2
     &             - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &             + 2.0d0*(xu*ydb-xdb*yu)*dphidzibu/ru2
               dzibxid = rcbxu*dphidzibu - dphidyu
     &                      - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
               dzibyid = rcbyu*dphidzibu + dphidxu
     &                      + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
               dxicyid = rcbyu*dphidxicu + dphidzu
     &                      + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
               dxiczid = rcbzu*dphidxicu - dphidyu
     &                      - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
               dyicxid = rcbxu*dphidyicu - dphidzu
     &                      - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
               dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
               dyiczid = rcbzu*dphidyicu + dphidxu
     &                      + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
               dzicxid = rcbxu*dphidzicu + dphidyu
     &                      + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
               dzicyid = rcbyu*dphidzicu - dphidxu
     &                      - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
               dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
               dxidxid = rcbxu*dphidxid
               dxidyid = rcbxu*dphidyid + zcb*rcb/ru2
               dxidzid = rcbxu*dphidzid - ycb*rcb/ru2
               dyidyid = rcbyu*dphidyid
               dyidzid = rcbyu*dphidzid + xcb*rcb/ru2
               dzidzid = rcbzu*dphidzid
c
c     get some second derivative chain rule terms by difference
c
               dxibxib = -dxiaxib - dxibxic - dxibxid
               dxibyib = -dyiaxib - dxibyic - dxibyid
               dxibzib = -dxiazib - dzibxic - dzibxid
               dxibzic = -dziaxib - dxibzib - dxibzid
               dyibyib = -dyiayib - dyibyic - dyibyid
               dyibzic = -dziayib - dyibzib - dyibzid
               dzibzib = -dziazib - dzibzic - dzibzid
               dzibyic = -dyiazib - dyibzib - dzibyid
               dxicxic = -dxiaxic - dxibxic - dxicxid
               dxicyic = -dyiaxic - dyibxic - dxicyid
               dxiczic = -dziaxic - dzibxic - dxiczid
               dyicyic = -dyiayic - dyibyic - dyicyid
               dyiczic = -dziayic - dzibyic - dyiczid
               dziczic = -dziazic - dzibzic - dziczid
c
c     now, increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                              + d2edphi2*dphidxia*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                              + d2edphi2*dphidyia*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                              + d2edphi2*dphidzia*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxia*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                              + d2edphi2*dphidyia*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                              + d2edphi2*dphidzia*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                              + d2edphi2*dphidxia*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                              + d2edphi2*dphidyia*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                              + d2edphi2*dphidzia*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                              + d2edphi2*dphidxia*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                              + d2edphi2*dphidyia*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                              + d2edphi2*dphidzia*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxia*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                              + d2edphi2*dphidyia*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                              + d2edphi2*dphidzia*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                              + d2edphi2*dphidxia*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                              + d2edphi2*dphidyia*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                              + d2edphi2*dphidzia*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                              + d2edphi2*dphidxia*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                              + d2edphi2*dphidyia*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                              + d2edphi2*dphidzia*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxia*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                              + d2edphi2*dphidyia*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                              + d2edphi2*dphidzia*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                              + d2edphi2*dphidxia*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                              + d2edphi2*dphidyia*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                              + d2edphi2*dphidzia*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                              + d2edphi2*dphidxia*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                              + d2edphi2*dphidyia*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                              + d2edphi2*dphidzia*dphidzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                              + d2edphi2*dphidxib*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                              + d2edphi2*dphidyib*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                              + d2edphi2*dphidzib*dphidzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxib*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                              + d2edphi2*dphidyib*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                              + d2edphi2*dphidzib*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                              + d2edphi2*dphidxib*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                              + d2edphi2*dphidyib*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                              + d2edphi2*dphidzib*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                              + d2edphi2*dphidxib*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                              + d2edphi2*dphidyib*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                              + d2edphi2*dphidzib*dphidzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                              + d2edphi2*dphidxib*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                              + d2edphi2*dphidyib*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                              + d2edphi2*dphidzib*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                              + d2edphi2*dphidxib*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                              + d2edphi2*dphidyib*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                              + d2edphi2*dphidzib*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                              + d2edphi2*dphidxib*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                              + d2edphi2*dphidyib*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                              + d2edphi2*dphidzib*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                              + d2edphi2*dphidxib*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                              + d2edphi2*dphidyib*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                              + d2edphi2*dphidzib*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                              + d2edphi2*dphidxib*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                              + d2edphi2*dphidyib*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                              + d2edphi2*dphidzib*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                              + d2edphi2*dphidxib*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                              + d2edphi2*dphidyib*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                              + d2edphi2*dphidzib*dphidzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                              + d2edphi2*dphidxic*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                              + d2edphi2*dphidyic*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                              + d2edphi2*dphidzic*dphidzic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxic*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                              + d2edphi2*dphidyic*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                              + d2edphi2*dphidzic*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                              + d2edphi2*dphidxic*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                              + d2edphi2*dphidyic*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                              + d2edphi2*dphidzic*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                              + d2edphi2*dphidxic*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                              + d2edphi2*dphidyic*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                              + d2edphi2*dphidzic*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                              + d2edphi2*dphidxic*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                              + d2edphi2*dphidyic*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                              + d2edphi2*dphidzic*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                              + d2edphi2*dphidxic*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                              + d2edphi2*dphidyic*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                              + d2edphi2*dphidzic*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                              + d2edphi2*dphidxic*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                              + d2edphi2*dphidyic*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                              + d2edphi2*dphidzic*dphidzib
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                              + d2edphi2*dphidxic*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                              + d2edphi2*dphidyic*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                              + d2edphi2*dphidzic*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                              + d2edphi2*dphidxic*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                              + d2edphi2*dphidyic*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                              + d2edphi2*dphidzic*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                              + d2edphi2*dphidxic*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                              + d2edphi2*dphidyic*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                              + d2edphi2*dphidzic*dphidzid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                              + d2edphi2*dphidxid*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                              + d2edphi2*dphidyid*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                              + d2edphi2*dphidzid*dphidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxid*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                              + d2edphi2*dphidyid*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                              + d2edphi2*dphidzid*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                              + d2edphi2*dphidxid*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                              + d2edphi2*dphidyid*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                              + d2edphi2*dphidzid*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                              + d2edphi2*dphidxid*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                              + d2edphi2*dphidyid*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                              + d2edphi2*dphidzid*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                              + d2edphi2*dphidxid*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                              + d2edphi2*dphidyid*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                              + d2edphi2*dphidzid*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                              + d2edphi2*dphidxid*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                              + d2edphi2*dphidyid*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                              + d2edphi2*dphidzid*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                              + d2edphi2*dphidxid*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                              + d2edphi2*dphidyid*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                              + d2edphi2*dphidzid*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                              + d2edphi2*dphidxid*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                              + d2edphi2*dphidyid*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                              + d2edphi2*dphidzid*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                              + d2edphi2*dphidxid*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                              + d2edphi2*dphidyid*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                              + d2edphi2*dphidzid*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                              + d2edphi2*dphidxid*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                              + d2edphi2*dphidyid*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                              + d2edphi2*dphidzid*dphidzic
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eimprop3  --  imp. dihedral energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eimprop3" calculates the improper dihedral potential
c     energy; also partitions the energy terms among the atoms
c
c
      subroutine eimprop3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'improp.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'usage.i'
      integer i,ia,ib,ic,id
      real*8 e,ideal,force
      real*8 dt,cosine,sine
      real*8 angle,fgrp
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb,rcb
      real*8 xdc,ydc,zdc
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      logical header,huge,proceed
c
c
c     zero out improper dihedral energy and partitioning terms
c
      neid = 0
      eid = 0.0d0
      do i = 1, n
         aeid(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the improper dihedral angle energy term
c
      do i = 1, niprop
         ia = iiprop(1,i)
         ib = iiprop(2,i)
         ic = iiprop(3,i)
         id = iiprop(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the improper dihedral angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     set the improper dihedral parameters for this angle
c
               ideal = vprop(i)
               force = kprop(i)
               dt = angle - ideal
               dowhile (dt .gt. 180.0d0)
                  dt = dt - 360.0d0
               end do
               dowhile (dt .lt. -180.0d0)
                  dt = dt + 360.0d0
               end do
               dt = dt / radian
c
c     calculate the improper dihedral energy
c
               e = force * dt**2
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total improper dihedral energy
c
               neid = neid + 1
               eid = eid + e
               aeid(ib) = aeid(ib) + 0.5d0*e
               aeid(ic) = aeid(ic) + 0.5d0*e
c
c     print a warning if the energy of this angle is large
c
               huge = (e .gt. 5.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Improper Dihedral',
     &                          ' Interactions :',
     &                       //,' Type',21x,'Atom Names',20x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                             name(ic),id,name(id),angle,e
   20             format (' Improper ',i5,'-',a3,3(1x,i5,'-',a3),
     &                       2x,f10.4,f12.4)
               end if
            end if
         end if
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
c     #########################################################
c     ##                                                     ##
c     ##  subroutine eimptor  --  improper torsional energy  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "eimptor" calculates the improper torsional potential energy
c
c
      subroutine eimptor
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'imptor.i'
      include 'torpot.i'
      include 'usage.i'
      integer i,ia,ib,ic,id
      real*8 e,fgrp
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 cosine,cosine2,cosine3
      real*8 sine,sine2,sine3
      real*8 phi1,phi2,phi3
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb,rcb
      real*8 xdc,ydc,zdc
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      logical proceed
c
c
c     zero out improper torsional energy
c
      eit = 0.0d0
c
c     calculate the improper torsional angle energy term
c
      do i = 1, nitors
         ia = iitors(1,i)
         ib = iitors(2,i)
         ic = iitors(3,i)
         id = iitors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the improper torsional parameters for this angle
c
               v1 = itors1(1,i)
               c1 = itors1(3,i)
               s1 = itors1(4,i)
               v2 = itors2(1,i)
               c2 = itors2(3,i)
               s2 = itors2(4,i)
               v3 = itors3(1,i)
               c3 = itors3(3,i)
               s3 = itors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
c
c     calculate the improper torsional energy for this angle
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total torsional angle energy
c
               eit = eit + e
            end if
         end if
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eimptor1  --  impr. torsion energy & gradient  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eimptor1" calculates improper torsional energy and its
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine eimptor1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'deriv.i'
      include 'energi.i'
      include 'group.i'
      include 'imptor.i'
      include 'torpot.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ia,ib,ic,id
      real*8 e,dedphi,fgrp
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 cosine,cosine2,cosine3
      real*8 sine,sine2,sine3
      real*8 phi1,phi2,phi3
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb,rcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      real*8 dphi1,dphi2,dphi3
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      logical proceed
c
c
c     zero out energy and first derivative components
c
      eit = 0.0d0
      do i = 1, n
         deit(1,i) = 0.0d0
         deit(2,i) = 0.0d0
         deit(3,i) = 0.0d0
      end do
c
c     calculate the improper torsional angle energy term
c
      do i = 1, nitors
         ia = iitors(1,i)
         ib = iitors(2,i)
         ic = iitors(3,i)
         id = iitors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the improper torsional parameters for this angle
c
               v1 = itors1(1,i)
               c1 = itors1(3,i)
               s1 = itors1(4,i)
               v2 = itors2(1,i)
               c2 = itors2(3,i)
               s2 = itors2(4,i)
               v3 = itors3(1,i)
               c3 = itors3(3,i)
               s3 = itors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
c
c     calculate improper torsion energy and master chain rule term
c
               e = torsunit * (v1*phi1+v2*phi2+v3*phi3)
               dedphi = torsunit * (v1*dphi1+v2*dphi2+v3*dphi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedphi = dedphi * fgrp
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = yca*dphidzt - zca*dphidyt
     &                       + zdc*dphidyu - ydc*dphidzu
               dphidyib = zca*dphidxt - xca*dphidzt
     &                       + xdc*dphidzu - zdc*dphidxu
               dphidzib = xca*dphidyt - yca*dphidxt
     &                       + ydc*dphidxu - xdc*dphidyu
               dphidxic = zba*dphidyt - yba*dphidzt
     &                       + ydb*dphidzu - zdb*dphidyu
               dphidyic = xba*dphidzt - zba*dphidxt
     &                       + zdb*dphidxu - xdb*dphidzu
               dphidzic = yba*dphidxt - xba*dphidyt
     &                       + xdb*dphidyu - ydb*dphidxu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     compute first derivative components for this angle
c
               dedxia = dedphi * dphidxia
               dedyia = dedphi * dphidyia
               dedzia = dedphi * dphidzia
               dedxib = dedphi * dphidxib
               dedyib = dedphi * dphidyib
               dedzib = dedphi * dphidzib
               dedxic = dedphi * dphidxic
               dedyic = dedphi * dphidyic
               dedzic = dedphi * dphidzic
               dedxid = dedphi * dphidxid
               dedyid = dedphi * dphidyid
               dedzid = dedphi * dphidzid
c
c     increment the improper torsion energy and gradient
c
               eit = eit + e
               deit(1,ia) = deit(1,ia) + dedxia
               deit(2,ia) = deit(2,ia) + dedyia
               deit(3,ia) = deit(3,ia) + dedzia
               deit(1,ib) = deit(1,ib) + dedxib
               deit(2,ib) = deit(2,ib) + dedyib
               deit(3,ib) = deit(3,ib) + dedzib
               deit(1,ic) = deit(1,ic) + dedxic
               deit(2,ic) = deit(2,ic) + dedyic
               deit(3,ic) = deit(3,ic) + dedzic
               deit(1,id) = deit(1,id) + dedxid
               deit(2,id) = deit(2,id) + dedyid
               deit(3,id) = deit(3,id) + dedzid
c
c     increment the components of the virial
c
               if (isobaric) then
                  virx = virx - xba*dedxia + xdc*dedxid
     &                      + xcb*(dedxic+dedxid)
                  viry = viry - yba*dedyia + ydc*dedyid
     &                      + ycb*(dedyic+dedyid)
                  virz = virz - zba*dedzia + zdc*dedzid
     &                      + zcb*(dedzic+dedzid)
               end if
            end if
         end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eimptor2  --  atom-wise imp. torsional Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eimptor2" calculates second derivatives of the improper
c     torsional energy for a single atom
c
c
      subroutine eimptor2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'group.i'
      include 'hessn.i'
      include 'imptor.i'
      include 'torpot.i'
      integer i,ia,ib,ic,id,kitors
      real*8 dedphi,d2edphi2,fgrp
      real*8 v1,v2,v3,c1,c2,c3,s1,s2,s3
      real*8 cosine,cosine2,cosine3
      real*8 sine,sine2,sine3
      real*8 xia,yia,zia,xib,yib,zib
      real*8 xic,yic,zic,xid,yid,zid
      real*8 xba,yba,zba,xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 xt,yt,zt,xu,yu,zu,xtu,ytu,ztu
      real*8 rt2,ru2,rtru,rcb
      real*8 dphi1,dphi2,dphi3
      real*8 d2phi1,d2phi2,d2phi3
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 dphidxibt,dphidyibt,dphidzibt
      real*8 dphidxibu,dphidyibu,dphidzibu
      real*8 dphidxict,dphidyict,dphidzict
      real*8 dphidxicu,dphidyicu,dphidzicu
      real*8 dxiaxia,dyiayia,dziazia,dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic,dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia,dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic,dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib,dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib,dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic,dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid,dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid,dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic,dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid,dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid,dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid,dzicxid,dzicyid,dziczid
      logical proceed
c
c
c     calculate the improper torsional angle energy term
c
      do kitors = 1, nitors
         ia = iitors(1,kitors)
         ib = iitors(2,kitors)
         ic = iitors(3,kitors)
         id = iitors(4,kitors)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (i.eq.ia .or. i.eq.ib .or.
     &                              i.eq.ic .or. i.eq.id)
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the improper torsional parameters for this angle
c
               v1 = itors1(1,kitors)
               c1 = itors1(3,kitors)
               s1 = itors1(4,kitors)
               v2 = itors2(1,kitors)
               c2 = itors2(3,kitors)
               s2 = itors2(4,kitors)
               v3 = itors3(1,kitors)
               c3 = itors3(3,kitors)
               s3 = itors3(4,kitors)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
               d2phi1 = -(cosine*c1 + sine*s1)
               d2phi2 = -4.0d0 * (cosine2*c2 + sine2*s2)
               d2phi3 = -9.0d0 * (cosine3*c3 + sine3*s3)
c
c     calculate the improper torsion master chain rule terms
c
               dedphi = torsunit * (v1*dphi1+v2*dphi2+v3*dphi3)
               d2edphi2 = torsunit * (v1*d2phi1+v2*d2phi2+v3*d2phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
               end if
c
c     abbreviations for first derivative chain rule terms
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
               xycb2 = xcb*xcb + ycb*ycb
               xzcb2 = xcb*xcb + zcb*zcb
               yzcb2 = ycb*ycb + zcb*zcb
               rcbxt = -2.0d0 * rcb * dphidxt
               rcbyt = -2.0d0 * rcb * dphidyt
               rcbzt = -2.0d0 * rcb * dphidzt
               rcbt2 = rcb * rt2
               rcbxu = 2.0d0 * rcb * dphidxu
               rcbyu = 2.0d0 * rcb * dphidyu
               rcbzu = 2.0d0 * rcb * dphidzu
               rcbu2 = rcb * ru2
               dphidxibt = yca*dphidzt - zca*dphidyt
               dphidxibu = zdc*dphidyu - ydc*dphidzu
               dphidyibt = zca*dphidxt - xca*dphidzt
               dphidyibu = xdc*dphidzu - zdc*dphidxu
               dphidzibt = xca*dphidyt - yca*dphidxt
               dphidzibu = ydc*dphidxu - xdc*dphidyu
               dphidxict = zba*dphidyt - yba*dphidzt
               dphidxicu = ydb*dphidzu - zdb*dphidyu
               dphidyict = xba*dphidzt - zba*dphidxt
               dphidyicu = zdb*dphidxu - xdb*dphidzu
               dphidzict = yba*dphidxt - xba*dphidyt
               dphidzicu = xdb*dphidyu - ydb*dphidxu
c
c     chain rule terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = dphidxibt + dphidxibu
               dphidyib = dphidyibt + dphidyibu
               dphidzib = dphidzibt + dphidzibu
               dphidxic = dphidxict + dphidxicu
               dphidyic = dphidyict + dphidyicu
               dphidzic = dphidzict + dphidzicu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     chain rule terms for second derivative components
c
               dxiaxia = rcbxt*dphidxia
               dxiayia = rcbxt*dphidyia - zcb*rcb/rt2
               dxiazia = rcbxt*dphidzia + ycb*rcb/rt2
               dxiaxib = rcbxt*dphidxibt + xcb*(zca*ycb-yca*zcb)/rcbt2
               dxiayib = rcbxt*dphidyibt + dphidzt
     &                      + (xca*zcb*xcb+zca*yzcb2)/rcbt2
               dxiazib = rcbxt*dphidzibt - dphidyt
     &                      - (xca*ycb*xcb+yca*yzcb2)/rcbt2
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2
               dxiayic = rcbxt*dphidyict - dphidzt
     &                      - (xba*zcb*xcb+zba*yzcb2)/rcbt2
               dxiazic = rcbxt*dphidzict + dphidyt
     &                      + (xba*ycb*xcb+yba*yzcb2)/rcbt2
               dxiaxid = 0.0d0
               dxiayid = 0.0d0
               dxiazid = 0.0d0
               dyiayia = rcbyt*dphidyia
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2
               dyiaxib = rcbyt*dphidxibt - dphidzt
     &                      - (yca*zcb*ycb+zca*xzcb2)/rcbt2
               dyiayib = rcbyt*dphidyibt + ycb*(xca*zcb-zca*xcb)/rcbt2
               dyiazib = rcbyt*dphidzibt + dphidxt
     &                      + (yca*xcb*ycb+xca*xzcb2)/rcbt2
               dyiaxic = rcbyt*dphidxict + dphidzt
     &                      + (yba*zcb*ycb+zba*xzcb2)/rcbt2
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2
               dyiazic = rcbyt*dphidzict - dphidxt
     &                      - (yba*xcb*ycb+xba*xzcb2)/rcbt2
               dyiaxid = 0.0d0
               dyiayid = 0.0d0
               dyiazid = 0.0d0
               dziazia = rcbzt*dphidzia
               dziaxib = rcbzt*dphidxibt + dphidyt
     &                      + (zca*ycb*zcb+yca*xycb2)/rcbt2
               dziayib = rcbzt*dphidyibt - dphidxt
     &                      - (zca*xcb*zcb+xca*xycb2)/rcbt2
               dziazib = rcbzt*dphidzibt + zcb*(yca*xcb-xca*ycb)/rcbt2
               dziaxic = rcbzt*dphidxict - dphidyt
     &                      - (zba*ycb*zcb+yba*xycb2)/rcbt2
               dziayic = rcbzt*dphidyict + dphidxt
     &                      + (zba*xcb*zcb+xba*xycb2)/rcbt2
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2
               dziaxid = 0.0d0
               dziayid = 0.0d0
               dziazid = 0.0d0
               dxibxic = -xcb*dphidxib/(rcb*rcb)
     &             - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidxibt/rt2
     &             - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidxibu/ru2
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
     &             - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidxibt/rt2
     &             + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidxibu/ru2
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2
               dxibyid = rcbyu*dphidxibu - dphidzu
     &                      - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
               dxibzid = rcbzu*dphidxibu + dphidyu
     &                      + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
               dyibzib = ycb*dphidzib/(rcb*rcb)
     &             - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &             - 2.0d0*(xt*zca-xca*zt)*dphidzibt/rt2
     &             + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &             + 2.0d0*(xu*zdc-xdc*zu)*dphidzibu/ru2
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu
     &             + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidyibt/rt2
     &             - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidyibu/ru2
               dyibyic = -ycb*dphidyib/(rcb*rcb)
     &             - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidyibt/rt2
     &             - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidyibu/ru2
               dyibxid = rcbxu*dphidyibu + dphidzu
     &                      + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2
               dyibzid = rcbzu*dphidyibu - dphidxu
     &                      - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu
     &             - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidzibt/rt2
     &             + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidzibu/ru2
               dzibzic = -zcb*dphidzib/(rcb*rcb)
     &             - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &             - 2.0d0*(xt*yba-xba*yt)*dphidzibt/rt2
     &             - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &             + 2.0d0*(xu*ydb-xdb*yu)*dphidzibu/ru2
               dzibxid = rcbxu*dphidzibu - dphidyu
     &                      - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
               dzibyid = rcbyu*dphidzibu + dphidxu
     &                      + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
               dxicyid = rcbyu*dphidxicu + dphidzu
     &                      + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
               dxiczid = rcbzu*dphidxicu - dphidyu
     &                      - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
               dyicxid = rcbxu*dphidyicu - dphidzu
     &                      - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
               dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
               dyiczid = rcbzu*dphidyicu + dphidxu
     &                      + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
               dzicxid = rcbxu*dphidzicu + dphidyu
     &                      + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
               dzicyid = rcbyu*dphidzicu - dphidxu
     &                      - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
               dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
               dxidxid = rcbxu*dphidxid
               dxidyid = rcbxu*dphidyid + zcb*rcb/ru2
               dxidzid = rcbxu*dphidzid - ycb*rcb/ru2
               dyidyid = rcbyu*dphidyid
               dyidzid = rcbyu*dphidzid + xcb*rcb/ru2
               dzidzid = rcbzu*dphidzid
c
c     get some second derivative chain rule terms by difference
c
               dxibxib = -dxiaxib - dxibxic - dxibxid
               dxibyib = -dyiaxib - dxibyic - dxibyid
               dxibzib = -dxiazib - dzibxic - dzibxid
               dxibzic = -dziaxib - dxibzib - dxibzid
               dyibyib = -dyiayib - dyibyic - dyibyid
               dyibzic = -dziayib - dyibzib - dyibzid
               dzibzib = -dziazib - dzibzic - dzibzid
               dzibyic = -dyiazib - dyibzib - dzibyid
               dxicxic = -dxiaxic - dxibxic - dxicxid
               dxicyic = -dyiaxic - dyibxic - dxicyid
               dxiczic = -dziaxic - dzibxic - dxiczid
               dyicyic = -dyiayic - dyibyic - dyicyid
               dyiczic = -dziayic - dzibyic - dyiczid
               dziczic = -dziazic - dzibzic - dziczid
c
c     now, increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                              + d2edphi2*dphidxia*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                              + d2edphi2*dphidxia*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                              + d2edphi2*dphidyia*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                              + d2edphi2*dphidxia*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                              + d2edphi2*dphidyia*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                              + d2edphi2*dphidzia*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxia*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                              + d2edphi2*dphidyia*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                              + d2edphi2*dphidzia*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                              + d2edphi2*dphidxia*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                              + d2edphi2*dphidyia*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                              + d2edphi2*dphidzia*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                              + d2edphi2*dphidxia*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                              + d2edphi2*dphidyia*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                              + d2edphi2*dphidzia*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxia*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                              + d2edphi2*dphidyia*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                              + d2edphi2*dphidzia*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                              + d2edphi2*dphidxia*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                              + d2edphi2*dphidyia*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                              + d2edphi2*dphidzia*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                              + d2edphi2*dphidxia*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                              + d2edphi2*dphidyia*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                              + d2edphi2*dphidzia*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxia*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                              + d2edphi2*dphidyia*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                              + d2edphi2*dphidzia*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                              + d2edphi2*dphidxia*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                              + d2edphi2*dphidyia*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                              + d2edphi2*dphidzia*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                              + d2edphi2*dphidxia*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                              + d2edphi2*dphidyia*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                              + d2edphi2*dphidzia*dphidzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                              + d2edphi2*dphidxib*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                              + d2edphi2*dphidxib*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                              + d2edphi2*dphidyib*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                              + d2edphi2*dphidxib*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                              + d2edphi2*dphidyib*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                              + d2edphi2*dphidzib*dphidzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                              + d2edphi2*dphidxib*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                              + d2edphi2*dphidyib*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                              + d2edphi2*dphidzib*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                              + d2edphi2*dphidxib*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                              + d2edphi2*dphidyib*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                              + d2edphi2*dphidzib*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                              + d2edphi2*dphidxib*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                              + d2edphi2*dphidyib*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                              + d2edphi2*dphidzib*dphidzia
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                              + d2edphi2*dphidxib*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                              + d2edphi2*dphidyib*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                              + d2edphi2*dphidzib*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                              + d2edphi2*dphidxib*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                              + d2edphi2*dphidyib*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                              + d2edphi2*dphidzib*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                              + d2edphi2*dphidxib*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                              + d2edphi2*dphidyib*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                              + d2edphi2*dphidzib*dphidzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                              + d2edphi2*dphidxib*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                              + d2edphi2*dphidyib*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                              + d2edphi2*dphidzib*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                              + d2edphi2*dphidxib*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                              + d2edphi2*dphidyib*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                              + d2edphi2*dphidzib*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                              + d2edphi2*dphidxib*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                              + d2edphi2*dphidyib*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                              + d2edphi2*dphidzib*dphidzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                              + d2edphi2*dphidxic*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                              + d2edphi2*dphidxic*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                              + d2edphi2*dphidyic*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                              + d2edphi2*dphidxic*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                              + d2edphi2*dphidyic*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                              + d2edphi2*dphidzic*dphidzic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                              + d2edphi2*dphidxic*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                              + d2edphi2*dphidyic*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                              + d2edphi2*dphidzic*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                              + d2edphi2*dphidxic*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                              + d2edphi2*dphidyic*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                              + d2edphi2*dphidzic*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                              + d2edphi2*dphidxic*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                              + d2edphi2*dphidyic*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                              + d2edphi2*dphidzic*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                              + d2edphi2*dphidxic*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                              + d2edphi2*dphidyic*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                              + d2edphi2*dphidzic*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                              + d2edphi2*dphidxic*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                              + d2edphi2*dphidyic*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                              + d2edphi2*dphidzic*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                              + d2edphi2*dphidxic*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                              + d2edphi2*dphidyic*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                              + d2edphi2*dphidzic*dphidzib
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                              + d2edphi2*dphidxic*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                              + d2edphi2*dphidyic*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                              + d2edphi2*dphidzic*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                              + d2edphi2*dphidxic*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                              + d2edphi2*dphidyic*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                              + d2edphi2*dphidzic*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                              + d2edphi2*dphidxic*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                              + d2edphi2*dphidyic*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                              + d2edphi2*dphidzic*dphidzid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                              + d2edphi2*dphidxid*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                              + d2edphi2*dphidxid*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                              + d2edphi2*dphidyid*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                              + d2edphi2*dphidxid*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                              + d2edphi2*dphidyid*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                              + d2edphi2*dphidzid*dphidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                              + d2edphi2*dphidxid*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                              + d2edphi2*dphidyid*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                              + d2edphi2*dphidzid*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                              + d2edphi2*dphidxid*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                              + d2edphi2*dphidyid*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                              + d2edphi2*dphidzid*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                              + d2edphi2*dphidxid*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                              + d2edphi2*dphidyid*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                              + d2edphi2*dphidzid*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                              + d2edphi2*dphidxid*dphidxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                              + d2edphi2*dphidyid*dphidxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                              + d2edphi2*dphidzid*dphidxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                              + d2edphi2*dphidxid*dphidyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                              + d2edphi2*dphidyid*dphidyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                              + d2edphi2*dphidzid*dphidyib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                              + d2edphi2*dphidxid*dphidzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                              + d2edphi2*dphidyid*dphidzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                              + d2edphi2*dphidzid*dphidzib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                              + d2edphi2*dphidxid*dphidxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                              + d2edphi2*dphidyid*dphidxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                              + d2edphi2*dphidzid*dphidxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                              + d2edphi2*dphidxid*dphidyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                              + d2edphi2*dphidyid*dphidyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                              + d2edphi2*dphidzid*dphidyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                              + d2edphi2*dphidxid*dphidzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                              + d2edphi2*dphidyid*dphidzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                              + d2edphi2*dphidzid*dphidzic
               end if
            end if
         end if
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eimptor3  --  impr. torsion energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eimptor3" calculates the improper torsional potential
c     energy; also partitions the energy terms among the atoms
c
c
      subroutine eimptor3
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'energi.i'
      include 'group.i'
      include 'imptor.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'torpot.i'
      include 'usage.i'
      integer i,ia,ib,ic,id
      real*8 e,angle,fgrp
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 cosine,cosine2,cosine3
      real*8 sine,sine2,sine3
      real*8 phi1,phi2,phi3
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb,rcb
      real*8 xdc,ydc,zdc
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      logical header,huge,proceed
c
c
c     zero out the torsional energy and partitioning terms
c
      neit = 0
      eit = 0.0d0
      do i = 1, n
         aeit(i) = 0.0d0
      end do
      header = .true.
c
c     calculate the improper torsional angle energy term
c
      do i = 1, nitors
         ia = iitors(1,i)
         ib = iitors(2,i)
         ic = iitors(3,i)
         id = iitors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,4,ia,ib,ic,id,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     set the improper torsional parameters for this angle
c
               v1 = itors1(1,i)
               c1 = itors1(3,i)
               s1 = itors1(4,i)
               v2 = itors2(1,i)
               c2 = itors2(3,i)
               s2 = itors2(4,i)
               v3 = itors3(1,i)
               c3 = itors3(3,i)
               s3 = itors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
c
c     calculate the improper torsional energy for this angle
c
               e = torsunit * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total torsional angle energy
c
               neit = neit + 1
               eit = eit + e
               aeit(ib) = aeit(ib) + 0.5d0*e
               aeit(ic) = aeit(ic) + 0.5d0*e
c
c     print a warning if the energy of this angle is large
c
               huge = (e .gt. 5.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,10)
   10                format (/,' Individual Improper Torsional',
     &                          ' Interactions :',
     &                       //,' Type',21x,'Atom Names',20x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,20)  ia,name(ia),ib,name(ib),ic,
     &                             name(ic),id,name(id),angle,e
   20             format (' Improper ',i5,'-',a3,3(1x,i5,'-',a3),
     &                       2x,f10.4,f12.4)
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
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine elj  --  Lennard-Jones van der Waals energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "elj" calculates the van der Waals interaction energy
c     using the Lennard-Jones 6-12 formalism
c
c
      subroutine elj
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
      real*8 e,p6,p12,rv,eps,rdn,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 rik,rik2,rik3,rik4,rik5,taper
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
c    compute the energy contribution for this interaction
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (skip(k) .eq. -i)  eps = eps / vdwscale
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0 * p6)
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
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0 * p6)
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj1  --  Lennard-Jones energy & derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj1" calculates the van der Waals energy and its first
c     derivatives with respect to Cartesian coordinates using
c     the Lennard-Jones 6-12 formalism
c
c
      subroutine elj1
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
      real*8 e,p6,p12,rv,eps,rdn,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 dedx,dedy,dedz,de
      real*8 rik,rik2,rik3,rik4,rik5,taper,dtaper
      real*8 xred(maxatm),yred(maxatm),zred(maxatm)
      logical proceed,iuse
c
      LOGICAL GOPARR,DSKWRK,MASWRK,gopart,mmonly,qmmm
      real*8 einter0,virx0,viry0,virz0
      integer ME,MASTER,NPROC,IBTYP,IPTIM,mparti
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /tinopt/ mparti,mmonly,qmmm
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
c        skip is accumulated and should be got on all nodes for safety
c        (the exact skipping logic being unclear to DGF).
         if(.not.gopart.or.mod(ii,nproc).eq.me) then
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
                  rik = sqrt(rik2)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0 * p6)
                  de = eps * (p12 - p6) * (-12.0d0/rik)
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
      do ii = 1, nvdw
         if(.not.gopart.or.mod(ii,nproc).eq.me) then
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
                     rik = sqrt(rik2)
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0 * p6)
                     de = eps * (p12 - p6) * (-12.0d0/rik)
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
         end if
      end do
      end if
      if(gopart) then
         call ddi_gsumf(4100,dev,3*n)
         call ddi_gsumf(4101,de14,3*n)
         call ddi_gsumf(4102,ev,1)
         call ddi_gsumf(4103,e14,1)
         einter=einter-einter0
         call ddi_gsumf(4104,einter,1)
         einter=einter+einter0
         if (isobaric) then
            virx=virx-virx0
            viry=viry-viry0
            virz=virz-virz0
            call ddi_gsumf(4105,virx,1)
            call ddi_gsumf(4106,viry,1)
            call ddi_gsumf(4107,virz,1)
            virx=virx+virx0
            viry=viry+viry0
            virz=virz+virz0
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
c     ##  subroutine elj2  --  atom-by-atom Lennard-Jones Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj2" calculates the van der Waals second derivatives for a
c     single atom at a time using the Lennard-Jones 6-12 formalism
c
c
      subroutine elj2 (iatom,xred,yred,zred)
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
      real*8 p6,p12,eps,rv
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 redi2,rediv2,rediiv
      real*8 redik,redivk,redikv,redivkv
      real*8 rik,rik2,rik3,rik4,rik5
      real*8 taper,dtaper,d2taper
      real*8 d2edx,d2edy,d2edz,term(3,3)
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
                  rik = sqrt(rik2)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  de = eps * (p12 - p6) * (-12.0d0/rik)
                  d2e = eps * (13.0d0*p12 - 7.0d0*p6) * (12.0d0/rik2)
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
                     rik = sqrt(rik2)
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     de = eps * (p12 - p6) * (-12.0d0/rik)
                     d2e = eps * (13.0d0*p12 - 7.0d0*p6) * (12.0d0/rik2)
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
c     ############################################################
c     ##                                                        ##
c     ##  subroutine elj3  --  Lennard-Jones energy & analysis  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "elj3" calculates the van der Waals interaction energy
c     using the Lennard-Jones 6-12 formalism and also partitions
c     the energy among the atoms
c
c
      subroutine elj3
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
      real*8 e,p6,p12,eps,rv,rdn,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 rik,rik2,rik3,rik4,rik5,taper
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
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0 * p6)
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
c     increment the total van der Waals energy components
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
   20                format (' VDW-LJ   ',i5,'-',a3,1x,i5,'-',a3,
     &                         12x,2f10.4,f12.4)
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
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0 * p6)
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
c     increment the total van der Waals energy components
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
   40                   format (' VDW-LJ   ',i5,'-',a3,1x,i5,'-',a3,
     &                            '   (XTAL)   ',2f10.4,f12.4)
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
c     ##  subroutine elj4  --  Lennard-Jones van der Waals energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj4" calculates the Lennard-Jones van der Waals interaction
c     energy using the method of lights to locate neighboring atoms
c
c
      subroutine elj4
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
      real*8 e,p6,p12,rv,eps,rdn,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 rik,rik2,rik3,rik4,rik5,taper
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
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0 * p6)
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj5  --  Lennard-Jones energy & derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj5" calculates the Lennard-Jones van der Waals interaction
c     energy and its first derivatives using the method of lights to
c     locate neighboring atoms
c
c
      subroutine elj5
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
      real*8 e,p6,p12,rv,eps,rdn,fgrp
      real*8 xi,yi,zi,xr,yr,zr
      real*8 redi,rediv,redk,redkv
      real*8 dedx,dedy,dedz,de
      real*8 rik,rik2,rik3,rik4,rik5,taper,dtaper
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
         map(i) = i
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
                  rik = sqrt(rik2)
                  p6 = rv**6 / rik2**3
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0 * p6)
                  de = eps * (p12 - p6) * (-12.0d0/rik)
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
