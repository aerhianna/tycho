      subroutine dynam(k,itest,ks)

c..   explicit second-order general relativistic hydrodynamics
c..   valid for up to mildly relativistic flows
c..   fails for large values of the relativistic kinematic factor gamma
c..   pseudoviscosity for shocks
c..   turbulent viscosity for explosions (wda 4-12-08)

      implicit none

      include 'dimenfile'
      include 'comod'
      include 'compu'
      include 'cphot'
      include 'cgtr.h'
      include 'cconst'
      include 'cnabla'

      real*8    ftest(4)

      integer*4 modek
ccccccccccccccccccc

      integer*4 j,k,n, nnn,ks,ki,itest,jj
      
      real*8    depth(kdm)
      real*8    hom(kdm),qt(kdm),vt2,rbar,ubar,delr
      real*8    ctime,dumi,tau,duphot,drphot,cth
      real*8    aom,fact
      real*8    t3p,t3m,t4p,t4m,akin 
      real*8    fako,pbar,aek,aet
      real*8    drdel,dtdel,dq,dp,fg1,fg2,fg3
      real*8    deltu

      real*8 fa,fb,fc,fd,fden,fak,ffak
      real*8 altmk,altpk,altmo,altpo, altm,altp
      real*8 ee(kdm),ff(kdm) ,akave(kdm)
c-------------------------------------------------------------------
c     gtr hydro, explicit mode
c..   compute new radii and specific volume  r(2,k),v(2,k)
c..   from velocities u(2,k)=u(3/2,k) which is centered between n=1,2

c..   peeling as an approximation for a recombination front
c..   matter outside the front is omitted from the grid
c..   assuming it is transparent and without radioactive heating
c..   other options are a two temperature solution in which radiation
c..   is treated as a bose-einstein gas with nonzero chemican potential
c..   or a full radiation transfer solution (both much more complex)
      if( modes .ne. 2 .and. it .eq. 1 )then
        depth(kk) = (r(1,kk)-r(1,kk-1))*ak(kk)/v(1,kk)
c       write(*,*)depth(kk)
        if( depth(kk) .lt. 0.20d0 )then
          do j = 1, kk-1
             k = kk - j +1
c.. optical depth
             depth(k) = (r(1,k)-r(1,k-1))*ak(k)/v(1,k)
             if( depth(k) .gt. 0.20d0 )then
                jj = k
                write(*,*)'MSG(DYNAM): peeling kk to ',jj,' from ',kk
                kk = jj
                go to 100
             endif
           enddo
         else
c..no change
         endif
      endif
100   continue
c..   end of peel algorithm

      if( t(1,kk) .lt. 2.01d3 )then
        write(*,'(i5,1p8e12.4)')kk,t(1,kk),depth(kk),ak(kk),
     1 tl(1,kk), xm(kk)/sol
        stop'a'
      endif
cccccccccccccccccc

      modek = 1

      if( it .eq. 1 )then
c..kinematic equation, Q and density
         do k = 2, kk
            td(k)  = t(2,k) - t(1,k)
            ctime  = dti(1) * efi(k)
            if( nouter .eq. 1 .and. k .eq. kk )then
               u(2,k) = 0
               dr(k)  = 0
            else
               u(2,k) = u(1,k) + du(k)*ctime
               dr(k)  =         u(2,k)*ctime
            endif
            r(2,k)  = r(1,k) + dr(k)
            a(k)    = pi4 * r(2,k)**2
            g(k)    = grav * xm(k) / r(2,k)**2
            v(2,k)  = pi43*(r(2,k)**3 - r(2,k-1)**3)/dmh(k)
     1           * 2.0d0/( gam(k) + gam(k-1) )
            dv(2,k) = v(2,k) - v(1,k)
            dumi    = dabs(dv(2,k)/v(1,k))

            if( dumi .gt. 5.0d-2  )then
               write(*,'(a12,2i5,1p11e11.3)')'dynam:vel ',
     2              k,it,r(1,k),dr(k),u(1,k),u(2,k),du(k),
     3              v(1,k),dv(2,k),ctime,dti(1)
            endif

            if( dumi .gt. 0.25d0 )then
               nnn = 3
               if( iter .gt. 1 )then
                  ftest(3) = dumi
                  go to 1000
               else
                  if( v(2,k) .le. 0.0d0 )then
                     write(*,*)'Dynam: negative V at ', k, l
                     goto 1000
                  endif
               endif
            endif
         enddo

c..   calculate pseudo-viscosity
c..   not updated in iteration
         do k = 3, kk
            if( dv(2,k) .ge. 0.0d0 )then
               q(2,k) = 0
            else
               deltu = u(2,k) - u(2,k-1)
c..   preserve homologous contraction
               if( u(2,k)/r(2,k) - u(2,k-1)/r(2,k-1) .lt. 0.0 )then
                  q(2,k) = qcons * deltu**2 / v(2,k)
               endif
            endif
c..   viscosity only in strong shocks
c            if( q(2,k) .lt. 1.0d-2*p(1,k)  ) q(2,k) = 0
         enddo

      endif

c..   turbulent pseudoviscosity = qt
c..   turbulent velocity = h
c..   homology = hom
      hom(1) = 0
      do k = 2, kk
         hom(k) = u(2,k)/r(2,k)
         qt(k) = 0
      enddo
      h(1) = 0
      h(kk) = 0
      do k = 2, kk-1
         rbar = 0.5d0*(r(2,k)+r(2,k-1))
         ubar = 0.5d0*(u(2,k)+u(2,k-1))
         delr = r(2,k) - r(2,k-1)
c         vt2 = - rbar*ubar*( rbar*( hom(k)-hom(k-1) )/delr 
c     1        + 0.5d0*( hom(k) + hom(k-1))  )
c         vt2 = -ubar*( u(2,k)-u(2,k-1) )
         vt2 = -ubar*( hom(k) - hom(k-1) )*rbar * 0.1d0
         if( vt2 .gt. 0.0d0 )then
            h(k) =sqrt( vt2 )
            if( dth(2) .gt. 0.5d0*delr/h(k) )then
               write(*,*)'dynam: dth too big ',k,dth(2),0.5d0*delr/h(k)
            endif

            qt(k) = vt2/v(2,k)
c            if( q(2,k) .eq. 0.0d0 )then
c               q(2,k) = qt(k)
c            endif

c            q(2,k) = q(2,k) + qt(k)
ccccccccccc
         endif

      enddo

c..   all values of it iteration

c..   evaluate ks = k at photosphere
      ks  = kk
      tau = 0
      do ki = 2, kk
         k = kk - ki + 2
         tau = ( r(2,k) - r(2,k-1) )*ak(k)/v(2,k) + tau
         if( tau .gt. 0.7d0 )then
            ks = k
            rphot  = r(2,ks) - v(2,ks)/ak(ks)
            rphot  = dmax1( rphot, r(2,ks-1) )
            duphot = u(2,ks) - u(2,ks-1)
            drphot = r(2,ks) - r(2,ks-1)
            uphot  = u(2,ks)
     1           + ( rphot - r(2,ks) ) * duphot / drphot
            goto 1001
         endif
      enddo
c..   total depth less than 0.7
      ks = 1
      rphot = r(2,1)
      uphot = u(2,1)
 1001 continue


c..   energy equation
c..   compute new temperature and luminosity implicitly 
c..   using Thomas tridiagonal method
      altp  = 0
      altm  = 0
      ee(1) = 0
      ff(1) = 0

      if( r(1,1) .gt. 0.0d0 )then
         a(1) = pi4*r(1,1)**2
         f(1) = tl(1,1)/a(1)
      else
         a(1) = 0.0d0
         f(1) = 0.0d0
      endif

      do k = 2, kk
c..save previous derivatives as we step through zone loop
         altpo = altp
         altmo = altm

         cth = dth(2)*efi(k)
         aom   = a(k)/dmi(k)
c..   optically thick/thin, flux-limiting and ionization
         t3p = t(2,k+1)**3
         t3m = t(2,k)**3
         t4p = t(2,k+1) * t3p
         t4m = t(2,k)   * t3m
c..keep opacity fixed during iteration
         if( it .le. 1 )then
            if( modek .eq. 0 )then
c..linear average
               akave(k) = 0.5d0*( ak(k) + ak(k+1) )
            elseif( modek .eq. 1 )then
c..logarithmic average
               akave(k) = sqrt( ak(k) * ak(k+1) )
            else
c..christy average: R.F.Christy, 1964, Rev.Mod.Phys. 36, 555
c..has correct limiting behavior for 
c..pseudo Kramers (kappa ~ T**-4 not T**-3.5)
c..recombination (kappa ~  T**12
c..and constant ak(k) = ak(k+1)
c..akave = 1/average( 1/kappa ) so this is inverse of RFC's expression
               akave(k) = (t4p + t4m)/( t4p/ak(k+1) + t4m/ak(k) )
            endif
         endif

c..flux limiting
c         akin =  1.0d0/( akave(k) + a(k)/(dmi(k)*3.0d0) )
         akin = 1.0d0/akave(k)
c         b(k)  = 0.0d0
c         h(k)  = 0.0d0
c         hp(k) = 0.0d0
c         ic(k) = 0
ccccccccccccccccccccccccccccccccccc

         ffak  = aom * cflux * (t4m - t4p)
         f(k)  = ffak * akin

         altpk =  aom * 4.0d0 * cflux * a(k) * t3p * akin
         altmk = -aom * 4.0d0 * cflux * a(k) * t3m * akin
         altp = altpk
         altm = altmk

         tl(2,k) = (f(k) + b(k))*a(k)

         dtl(k)  = tl(2,k) - tl(1,k)
         fako    = ( tl(1,k) - tl(1,k-1) )/dmh(k)
         fak     = ( tl(2,k) - tl(2,k-1) )/dmh(k)
         pbar    = 0.5d0*( p(1,k) + p(2,k) + q(1,k) + q(2,k) )

c..approximate energy generation by a power law in iterative loop
         fact = 
     1        ss(k) *(t(2,k)/t(1,k))**sa(k)  *(v(2,k)/v(1,k))**sb(k)
     2        + 
     3        snu(k)*(t(2,k)/t(1,k))**snua(k)*(v(2,k)/v(1,k))**snub(k)
c..   second order in time
c..   s(6,k) is actual source used in dynam and hstat
c..   s(3,k) is past value of s(4,k) + s(5,k) + s(7,k)
         s(6,k) = (1.0d0 - alphae)*s(3,k) + fact*alphae
         st(k) = ss(k)*sa(k)/t(1,k) + snu(k)*snua(k)/t(1,k)
         sv(k) = ss(k)*sb(k)/v(1,k) + snu(k)*snub(k)/v(1,k)

c..define tridiagonal coefficients ee and ff
c         if( q(1,k) .lt. 0.2d0*p(1,k) )then
            aek   = e(2,k) - e(1,k)
     1           + pbar*dv(2,k)
     2           + cth*fako*( 1.0d0 - alphae)
     3           + cth*fak * alphae    
     3           - cth*s(6,k)
            aet     = et(k) + 0.5d0*pt(k)*dv(2,k)
     1           - cth*st(k)*alphae
            fd    = aek
            fa    = -altp *cth/dmh(k)*alphae
            fc    =  altmo*cth/dmh(k)*alphae
            fb    = -aet + cth *(altm - altpo)/dmh(k)*alphae
c         else
c..   inside pseudoviscous shock - keep adiabatic
c..   no change in energy flow in pseudoviscous region
c            tl(2,k) = tl(1,k)
c            aek   = e(2,k) - e(1,k)
c     1           + pbar*dv(2,k)
c            aet     = et(k) + 0.5d0*pt(k)*dv(2,k)
c            fd    = aek
c            fa    = 0
c            fc    = 0
c            fb    = -aet 
c            write(*,'(i5,1p8e12.3)')k,aek,u(2,k),q(2,k),p(2,k)
ccccccccc
c         endif
         fden  = fb - fc*ee(k-1)
         ee(k) = fa/fden
         ff(k) = (fd + fc*ff(k-1))/fden
      enddo
c...........................end of k loop...............

c..   initialize arrays for saving worst test values
      do nnn = 1, 4
         wtest(nnn)  = 0
         nwtest(nnn) = 1
         ftest(nnn) = 0.0d0
      enddo

c..   reverse sweep      
      drdel = 0
      dtdel = 0
      do j = 1, kk-1
         k = kk+1 -j
         dtdel   = ee(k)*dtdel +ff(k)
         dt(2,k) = dt(2,k) + dtdel
         td(k)   = dtdel
         dtl(k)  = tl(2,k) - tl(1,k)

         if( abs( dt(2,k)/t(1,k) ) .gt. 0.1d0 )then
             write(*,'(a12,2i5,1p10e11.3)')'dynam: dT',
     1        k,it,t(1,k),dt(2,k),td(k),dt(2,k)/t(1,k),
     1           td(k)/t(1,k),tl(1,k),tl(2,k),akin,q(2,k)/p(2,k)
          endif
       enddo
c..............................end of loop in j (reverse k)

      do k = 2, kk
         t(2,k)  = t(1,k)  + dt(2,k)
         v(2,k)  = v(1,k)  + dv(2,k)
c..   use midstep estimate for pressure to get acceleration
         p(2,k)  = p(1,k) + pt(k)*dt(2,k) + pv(k)*dv(2,k)
         r(2,k)  = r(1,k)  + dr(k)
         tl(2,k) = tl(1,k) + dtl(k)
         s(1,k)  = (f(k-1)*a(k-1) - f(k)*a(k))/dmh(k)
         s(2,k)  = (b(k-1)*a(k-1) - b(k)*a(k))/dmh(k)
      enddo

      do k = 2, kk
c..   use new values to estimate new accelerations du(k) and
c..   new velocity u(2,k)
         dq      = q(2,k) - q(2,k+1)
         dp      = p(2,k) - p(2,k+1) + dq

c..   gtr hydro
         fg1   = grav  * gm(k) / r(2,k)**2
         fg2   =  dp *  a(k)/dmi(k) * gam(k) / grw(k)

         if( k .eq. 2 )then
            pbar = p(2, 2) + q(2, 2)
         elseif( k .eq. kk)then
            pbar = p(2,kk) + q(2,kk)
         else
            pbar = 0.5d0*(p(2,k+1) + p(2,k) + q(2,k+1) + q(2,k))
         endif

         fg3   = pi4 * grav * pbar * r(2,k) / crad2
c..   keep outer zone boundary fixed for nouter = 1
         if( nouter .eq. 1 .and. k .eq. kk )then
            du(k) = 0
         else
            du(k) = -fg1 + fg2 - fg3
         endif

c..   worst zones
         nnn = 1
         drdel    = u(2,k)*dth(2)*efi(k)
         ftest(1) = dabs(drdel/r(1,k))
         if( ftest(nnn) .gt. wtest(nnn) )then
            wtest(nnn)  = ftest(nnn)
            nwtest(nnn) = k
         endif
         if( ftest(1) .gt. 0.67d0 )then
            if( iter .gt. 1 )then
               go to 1000
            else
               if(      drdel+r(1,k) .gt. r(1,k+1)
     1              .or.  drdel+r(1,k) .le. r(1,k-1)
     2              .and. k .lt. kk )then
                  go to 1000
               endif
            endif
         endif
         if( ftest(nnn) .gt. 2.0d0*cdelt )then
            itest = nwtest(nnn)
         endif

         nnn = 3
         ftest(nnn) = dabs(td(k)/t(1,k))
         if( ftest(nnn) .gt. wtest(nnn) )then
            wtest(nnn)  = ftest(nnn)
            nwtest(nnn) = k
         endif
         if( ftest(nnn) .gt. resid )then
            itest = nwtest(nnn)
         endif
      enddo

c.........................................................

      if( it .eq. iter )then
         write(*,*)'DYNAM NONCONVERGENCE'
         write(*,'(2(a12,i5))')'kk',kk,'it',it
         write(*,'(a20,4i11))')'nwtest(n)',(nwtest(n),n=1,4)
         write(*,'(a20,1p4e11.3)')'wtest(n)',(wtest(n),n=1,4)

         k = nwtest(3)
         write(*,*)'K = nwtest(3)'
         write(*,'(a20,3i12)')'nwtest(3)',nwtest(3)-1,nwtest(3),
     1        nwtest(3)+1
         write(*,'(a20,1p3e12.4)')'q(1,k-1,k,k+1)',q(1,k-1),q(1,k),
     1        q(1,k+1)
         write(*,'(a20,1p3e12.4)')'q(2,k-1,k,k+1)',q(2,k-1),q(2,k),
     1        q(2,k+1)

         write(*,'(a20,1p3e12.4)')'q/p(1)',q(1,k-1)/p(1,k-1),
     1        q(1,k)/p(1,k),q(1,k+1)/p(1,k+1)
         write(*,'(a20,1p3e12.4)')'q/p(2)',q(2,k-1)/p(2,k-1),
     1        q(2,k)/p(2,k),q(2,k+1)/p(2,k+1)

         write(*,'(a20,1p3e12.4)')'ak(k)',ak(k-1),ak(k),ak(k+1)
         write(*,'(a20,1p8e12.4)')'mfp',v(2,k-1)/ak(k-1),v(2,k)/ak(k),
     1        v(2,k+1)/ak(k+1)
         write(*,'(a20,1p3e12.4)')'r(k+1)-r(k)',r(2,k-1)-r(2,k-2),
     1        r(2,k)-r(2,k-1),r(2,k+1)-r(2,k)
         write(*,'(a20,1p3e12.4)')'T(1)',t(1,k-1),t(1,k),t(1,k+1)
         write(*,'(a20,1p3e12.4)')'T(2)',t(2,k-1),t(2,k),t(2,k+1)

         write(*,'(a20,1p3e12.4)')'P(1)',p(1,k-1),p(1,k),p(1,k+1)
         write(*,'(a20,1p3e12.4)')'P(2)',p(2,k-1),p(2,k),p(2,k+1)

         write(*,'(a20,1p3e12.4)')'rho',1.0d0/v(2,k-1),1.0d0/v(2,k),
     1        1.0d0/v(2,k+1)        
         write(*,'(a20,6x,1p3e12.4)')'L(1)',tl(1,k-1),tl(1,k),tl(1,k+1)
         write(*,'(a20,6x,1p3e12.4)')'L(2)',tl(2,k-1),tl(2,k),tl(2,k+1)

         write(*,'(a20,6x,1p3e12.4)')'AF(2)',
     1        -a(k-1)*cflux*( t(2,k  )**4-t(2,k-1)**4 ),
     1        -a(k  )*cflux*( t(2,k+1)**4-t(2,k  )**4 ),
     1        -a(k+1)*cflux*( t(2,k+2)**4-t(2,k+1)**4 )



         write(*,'(a20,6x,1p3e12.4)')'R',r(2,k-1),r(2,k),r(2,k+1)
         write(*,'(a20,6x,1p3e12.4)')'u',u(2,k-1),u(2,k),u(2,k+1)
         write(*,'(a20,6x,1p3e12.4)')'du*dt',du(k-1)*dth(2),
     1        du(k)*dth(2),du(k+1)*dth(2)
         write(*,'(a20,6x,1p3e12.4)')'-Adp/dm',
     1        -a(k-1)*(p(2,k)-p(2,k-1))/dmi(k-1),
     2        -a(k)*(p(2,k+1)-p(2,k))/dmi(k),
     3        -a(k+1)*(p(2,k+2)-p(2,k+1))/dmi(k+1)
         write(*,'(a20,6x,1p3e12.4)')'du',du(k-1),du(k),du(k+1)
         write(*,'(a20,6x,1p3e12.4)')'g',g(k-1),g(k),g(k+1)
         write(*,'(a20,6x,1p3e12.4)')'AF(2)',
     1        -arad*( t(2,k  )**4-t(2,k-1)**4 )*a(k-1)/dmi(k-1)/3.0d0,
     2        -arad*( t(2,k+1)**4-t(2,k  )**4 )*a(k  )/dmi(k  )/3.0d0,
     3        -arad*( t(2,k+2)**4-t(2,k+1)**4 )*a(k+1)/dmi(k+1)/3.0d0


         stop'bbb-dynam'
      endif

      return

c------------------------------------panic exit, fact .gt. 2/3
 1000 continue
      itest = -1
      write(3,77) nnn, ftest(nnn),nwtest(nnn)
      write(3,70)j,k,nnn,it
      write(3,73)
      write(3,74)(ftest(n),n=1,4), resid

 70   format(4i5,' sweep index, zone, variable index,',
     1     ' iteration')
 72   format(' worst zones ', 4i5)
 73   format(' relative changes at bad zone, ',
     1     ' ftest(r,t,v,l), residual limit ' )
 74   format(1p10e12.4)
 77   format(' panic stop in iteration, variable =', i3,
     1     ',  value =', 1pe12.4, ', zone =', i4)

      return

      end

