      subroutine hydro(ks,no,mix,kxmx,dthc)

c--------------------------------------------------------
c..driver for hydrostatic or hydrodynamic evolution
c..controls envelope integrations
c..controls mass loss
c..mixing occurs prior to iteration
c--------------------------------------------------------
c..   revised 10-29-03
      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'csurface'
      include 'cphot'
      include 'cgtr.h'
      include 'cconst'
      include 'cnabla'
      include 'cburn'
      include 'cbug'
      include 'cenv'
      include 'ceoset'

      integer*4 mix,kxmx
      real*8 dthc

      integer*4 j,k,n, ks,no,nc,itest

      real*8    pbar,delv,delt
      real*8    etbar,evbar,dmasdt,dmass,dmfrak,aom
      real*8    tlum, atlmax
      integer*4 katlmax

      data cwtest/'r','L','T','V'/

c..   mode = 0 means explicit hydrodynamics
c..   mode = 1 means hydrostatics (heavily damped hydrodynamics)
c..   mode > 1 means implicit hydrodynamics of various sorts

c..   input:   l,ks,no
c..   output:    ks,no
c-----------------------------------------------------
c     write(*,*)'ENTERING HYDRO dth kk l ',dth(2),kk,l
c    1     ,x(nnuc,2),t(1,kk),t(1,kk+1)
c-----------------------------------------------------


c..   initialize, do not extrapolate
      do k = 1, kk+1
         r( 2,k) = r( 1,k)
         tl(2,k) = tl(1,k)
         dr( k ) = 0.0d0
         dtl(k ) = 0.0d0
         t( 2,k) = t( 1,k)
         v( 2,k) = v( 1,k)
         dt(2,k) = 0.0d0
         dv(2,k) = 0.0d0
      enddo
c      t(1,2) = 1.0d9
c      tl(1,2) = 1.0d40
c      t(2,2) = 1.0d9
c      tl(2,2) = 1.0d40
      dthc = 0.0d0

c..set past energy generation rate for initial timestep
      if( l .eq. 1 )then
         do k = 2,kk
            s(6,k) = ss(k) + snu(k)
         enddo
      endif

      do k = 2, kk
c..   keep old ssum in s(3,k) for centering in time
         if( no .ne. 0 )then
            s(3,k) = 0.0d0
         else
c..   s(3,k) = s(4,k) + s(5,k) + s(7,k)
            s(3,k) = s(6,k)
         endif
      enddo

c..   use thermal damping for first step, after nonconvergence,
c..   or with hydrostatic mode
      if( l .le. 1 .or. no .ne. 0 .or. mode .eq. 1 
     1     .or. mode .eq. 4 )then
         alphae = 1.0d0
      else
         alphae = 0.5d0
      endif

c..   reset error flag "no" after extrapolation tests
      no = 0

c..   reset value of it for tests in burn
      it   = 0

c..   set initial adiabatic gradient from past solution
      nc = 1
      do k = 2, kk
         pbar = 0.5d0 * (p(nc,k+1) + p(nc,k))
         delv = v(nc,k+1) - v(nc,k)
c..   y(k) = e(nc,k+1) - e(nc,k) + pbar*delv,    for Ledoux
c..   y(k)  = etbar*delt +  (pbar + evbar)*delv, for Schwarzschild
c..   z(k) must be defined consistently in hstat.f
         delt  =       t(nc,k+1) - t(nc,k)
         etbar = 0.5d0 * ( et(k+1) + et(k) )
         evbar = 0.5d0 * ( ev(k+1) + ev(k) )
         y(k)  =  e(nc,k+1) - e(nc,k) + pbar*delv
      enddo
      if( modes .eq. 2 )then
         y(kk)  =  0.0d0
      endif
      y(kk+1) = 0.0d0

c..   surface boundary condition
      if( modes .eq. 2 )then

c-------------------------------------------------------
c..   accretion or mass loss can be included (mloss .ne. 0)
c..   envelope lies between radii kk and kk+1
c..   put photospheric values outside kk+1 for T,P,rho
c..   integration in from photosphere

         if( mloss .ne. 0 )then

c..   update envelope mass for changes due to mass loss
            dmasdt = peryear * sol / secpy
            dmass  = dth(2)  * dmasdt
            dmfrak = 0.05d0  * dmh(kk+1)
            if( dmass .gt. dmfrak .or. dmass .lt. -dmfrak )then
               write(*,'(a30,1pe12.3)')
     1              'PANIC IN HYDRO: mass loss: dth =',dth(2)
               write(*,'(4(a15,1pe12.3))')
     1              'envelop mass',dmh(kk+1)/sol,
     2              'dmass/sol',dmass/sol,
     3              'peryear',peryear
               write(*,*)'reduce time step to keep d(Menv) small'
               write(*,'(4(a15,1pe12.3))')
     1              'dth est', -dmfrak/dmasdt,
     2              'new dmh(kk+)',(dmh(kk+1)-dmfrak)/sol
c..   update
               dth(2) = -dmfrak/dmasdt
c..   peryear is fixed, dth(2) reduced
               dmh(kk+1) = dmh(kk+1) - dmfrak
               write(*,'(a30,3(a12,1pe12.3))')'HYDRO: MAJOR MASS LOSS ',
     1              'new dth(2)',dth(2),'peryear',peryear,
     2              'new Menv',dmh(kk+1)/sol
c               stop'HYDRO: MASS LOSS'
ccccccccc

c               dth(2)    = 0.5d0 * dth(2)
c               dmh(kk+1) = dmh(kk+1) + 0.5d0 *  dmass
            else
               dmh(kk+1) = dmh(kk+1) + dmass
            endif

            xm(kk+1) = xm(kk) + dmh(kk+1)

         endif

c..   fits envelope to outer boundary kk
         
         call fitenv(1)

c..set new (n=2) values to initial values (iteration targets)
         r(2,kk+1)  = r(1,kk+1)
         tl(2,kk+1) = tl(1,kk+1)

c..   radiative flux
         f(kk) = cflux*a(kk)/dmi(kk)*(t(1,kk)**4 -t(1,kk+1)**4)
     1        * 2.0d0/( ak(kk) + ak(kk+1) )

c..   conserve angular momentum as envelope expands and contracts
         omeg(kk+1) = ajay(kk+1)/r(1,kk+1)**2
         rotshear(kk+1) = omeg(kk+1)*r(1,kk+1)

      elseif( modes .eq. 0 .or. modes .eq. 1 )then
c..   modes = 0, photosphere at kk
c..   modes = 1, photosphere moves through grid
c..   zero boundary
         dmh(kk+1) = 0.0d0
         dmi(kk)   = dmh(kk) * 0.5d0
         xm(kk+1)  = xm(kk) + dmh(kk+1)
         r(1,kk+1) = 0.0d0
         r(2,kk+1) = 0.0d0
         u(1,kk+1) = 0.0d0
         u(2,kk+1) = 0.0d0
         v(1,kk+1) = 0.0d0
         v(2,kk+1) = 0.0d0
         p(1,kk+1) = 0.0d0
         p(2,kk+1) = 0.0d0
         t(1,kk+1) = 0.0d0
         t(2,kk+1) = 0.0d0
         ak(kk+1)  = 0.0d0

         tlum = pi4*grav*xm(kk)/ak(kk)*cflux*t(1,kk)**4/p(1,kk)

         if( l .le. 1 )then
            tl(1,kk) = tlum
            tl(2,kk) = tlum
         endif
         tl(1,kk+1) = tl(1,kk)
         tl(2,kk+1) = tl(1,kk)

         etau   = 0.0d0
         epi    = 0.0d0
         dpidl  = 0.0d0
         dpidr  = 0.0d0
         dtaudl = 0.0d0
         dtaudr = 0.0d0
         enab   = 0.4d0
c..   average with envelope
         do n = 1, ndim
            x(n,kk+1) = x(n,kk)
            x(n,kk+2) = x(n,kk)
         enddo
         
      else
c..   other values
         write(*,*)'modes ',modes
         stop' define surface for this modes value'

      endif

c..   update as constant in iteration
      p(2,kk+1)  = p(1,kk+1)
      t(2,kk+1)  = t(1,kk+1)
      v(2,kk+1)  = v(1,kk+1)
      e(2,kk+1)  = e(1,kk+1)
      r(2,kk+1)  = r(1,kk+1)
      u(2,kk+1)  = u(1,kk+1)
      tl(2,kk+1) = tl(1,kk+1)

      if( modes .eq. 2 )then
         p(2,kk+2)  = p(1,kk+2)
         t(2,kk+2)  = t(1,kk+2)
         v(2,kk+2)  = v(1,kk+2)
         e(2,kk+2)  = e(1,kk+2)
         r(2,kk+2)  = r(1,kk+2)
c..   moves with envelope as a unit
         u(2,kk+2)  = u(1,kk+1)
         tl(2,kk+2) = tl(1,kk+2)
      endif

c..   mix the composition from converged (initial) solution
      if( mode .ne. 0 )call cmix(mix,kxmx,dthc,nc)

c..   set t,v values at end of time step interval, use cmix x
      nc   = 2
      call state(2,kk,nc)

c..   keep burn network solutions out of iterative loop
      call burn(nc,0,0,2,kk)

c..   only on first step so constant through iteration
      call gtrs(p,q,e,v,r,u,xm,dmi,dmh,kk,nc,newt,mode)

c..   update convective flags = ic(k) array and doux(k), dnab(k),
c..   etc., from n=2 values, just obtained from state and burn
      if( mode .ne. 0 )call cinit(2,kk,2)

      do it = 1, iter
c..   set t,v values at end of time step interval, use beginning x
c..   use updated t,v when available
         nc   = 2

         if( it .gt. 1 )call state(2,kk,nc)

c..   update convective flags = ic(k) array and doux(k), dnab(k),
c..   etc., from n=2 values, just obtained from state and burn
         call cinit(2,kk,nc)
c..   opacity from edd. boundary cond.
         ak(1)    = ak(2)
         a(kk)    = pi4 * r(2,kk)**2
c..   define surface opacity for dynamic modes
         if( modes .eq. 0 )then
c..   photosphere at kk
            aom      = a(kk)/dmi(kk)
            ak(kk+1) = aom / 1.5d0
         elseif( modes .eq. 1 )then
c..   photosphere moves through the grid
            ak(kk+1) = ak(kk)
         endif

c..   initialize flag to test iteration
         itest = 0

         if( mode .eq. 0) then

            call dynam(k,itest,ks)


         elseif( mode .eq. 1 .or. mode .eq. 2 )then

            call hstat2(k,itest,ks)


         else

            write(*,*)'error in hydro: mode=',mode
            stop'hydro: mode'

         endif

c..   debugging options........................................
         if( istep .ge. 2 .and. istep .le. kk )then
c..   writes zone istep, iteration it, convection flag ic(k), and
c..   selected zone variables R,L,T,V and their cumulatives changes
c..   for one zone istep, at each iteration
c..   this give lots of output, so use sparingly
            k = istep
            write(*,'(/a9,3i5,1p8e12.3)')'HYDRO',k,it,ic(k),
     1           r(2,k),dr(k),tl(2,k),dtl(k),
     2           t(2,k),dt(2,k),v(2,k),dv(2,k)
c..   add additional writes as needed here
            write(*,'(24x,1p8e12.3)')tl(2,k)/dmh(k),
     1           s(5,k),st(k)*dt(2,k),
     1           sv(k)*dv(2,k),a(k)*b(k)/dmh(k),x(nnuc-1,k),
     2           xd(nnuc-1,k)*dth(2)
         endif

c..   plots for debugging

         if( nbug .eq. 1 )then
c..realtime on iterations
            call dbplt(no)
         elseif( nbug .eq. 2 )then
c..pause after each iteration
            call dbplt(no)
            pause
         elseif( nbug .ne. 0 )then
            write(*,*)'HYDRO: nbug ',nbug
            stop'nbug error'
         endif


         if( itest .eq. 0 .or. iter .le. 1 )then
            atlmax = 0
            do k = 2, kk
               if( atlmax .lt.  abs( tl(2,k) ) )then
                  atlmax = abs( tl(2,k) )
                  katlmax = k
               endif
            enddo

            if( modes .eq. 2 )then
c..update target luminosity
               tl(2,kk+1) = tl(2,kk)
               tl(2,kk+2) = tl(2,kk)
c..   first order 2D interpolation to get photospheric radius r(2,kk+1)
c..   from a grid of join radii vr(nvmax(n),n) in a 5-pt stencil
c..   in luminosity and rjoin variables
c..   (correct a roughly 1 percent effect due to time staggering)
               r(2,kk+1) = vr(1,3)
     1              +(vr(1,2)-vr(1,5))/(vr(nvmax(2),2)-vr(nvmax(5),5))
     2              *(
     3              r(2,kk) - vr(nvmax(3),3) 
     4              -(vr(nvmax(1),1)-vr(nvmax(4),4))/(vl(1,1)-vl(1,4))
     5              *( tl(2,kk) - vl(1,3) )
     6              )
            endif

c-----------------------------------------------------
c      write(*,*)'LEAVING HYDRO dth kk L ',dth(2),kk,l
c     1     ,r(1,kk+2),r(2,kk+2)
c-----------------------------------------------------


            return
         endif


c..   failure for iterative solution
         if( itest .lt. 0 ) go to 1000

      enddo

c---------------------------------------------end of it loop

 1000 continue
c..   nonconvergence, panic on failure exit
      no = 1
      write(*,'( a25,4(a10,i6))')'HYDRO: nonconvergence at',
     1     'iteration',it,'step',l,'model',model
      write(3,'( a25,4(a10,i6))')'HYDRO: nonconvergence at',
     1     'iteration',it,
     1     'step',l,'model',model

      write(3,'(3a4,a3,8a11)')
     1     'j','k ','kk-k','ic','r','dr','t','dt','v','dv','tl','dtl'

      do j = 1, 4
         k = nwtest(j)
         if( k .gt. 0 .and. k .le. kk )then
            write(3,'(3i4,i3,1p8e11.3)') j, k,kk-k,ic(k), 
     1           r(1,k),dr(k),   t(1,k),dt(2,k),
     2           v(1,k),dv(2,k), tl(1,k),dtl(k)
            
         endif
      enddo

      write(*,*)
     1     'Leaving HYDRO in confusion! No update; will re-try.'

      return

      end


