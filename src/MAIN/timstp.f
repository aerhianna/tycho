c
c     
c     This subroutine calculates the new time step.
c     See notes below for limits considered.
c     
c     
      subroutine timstp(dthc,period,ktype,kkman,nc,kxmx)

c     
c     Choose new timestep
c     
c     REVISION HISTORY:
c     7-29-06
c     
c     
c     THIS ROUTINE NEEDS SERIOUS SPRING CLEANING!
c     (C.MEAKIN, 9-29-2008)
c     
c     
      implicit none
      
      include 'dimenfile'
      include 'comod'
      include 'ctmstp'
      include 'cburn'
      include 'compu'
      include 'cconst'
      include 'cenrchk'
      include 'cnabla'
      include 'cenv'
      
      real*8 period, ether, expnuc, taunuc, dthc, dlumin,alumin
      real*8 fact, fak, tman
      real*8 scr(kdm)
      real*8 zcno,dlnzcno
      
      integer*4 ktype, kkman, nc, i, kxmx
      integer*4 lo, knd, n, k, ix, ij
      integer*4 kchange
      
c..   restrict time step for many iterations
      integer*4 it1, it2
      
      real*8 tiny, huge
      parameter( tiny = 1.0d-40, huge = 1.0d99 )
      
c..   identifers for time step conditions
      data (dthflag(i),i=1,22)
     .     /'  ','CT','cx',' T',' V',' H','a3','He',' C',' N',' O',
     .     'Ne','Mg','Si','dm','dL','dt','it','mx','fx',
     .     'dR', 'zz' /
      data (dthflag(i),i=23,condim)
     .     /8*'  '/



c-------------------------------------------------------------------
c..   the limits are:
c..   CT is sound speed (Courant) condition
c..   cx is compression time implied from existing velocity
c..   T is fractional change in temperature
c..   V is fractional change in volume
c..   H,...,Ye are fractional changes in abundance
c..   tr is thermal runaway
c..   dm is mass change
c..   hk is Helmholz-Kelvin
c..   dt is the fractional change in timesep
c..   it is the iteration count
c..   mx is mixing from cmix.f and thoul.f
c..   fx is fixed dth
c..   zz is change in metallicity (CNO)
c--------------------------------------------------------------------
c..   nc = 2 is normal cycle, nc = 1 is first call to this subroutine
c..   use courant condition if others are 0
c..   mode=0  means hydrodynamic case
c..   mode=1  means hydrostatic case
c..   mode=2 means implicit hydrodynamics
c--------------------------------------------------------------------


c..   define lower bound on loop to omit hydro tests in hydrostatic case
      if( mode .ge. 1 )then
c..   compressional time enabled
         lo = 3
      else
         lo = 1
      endif
      

      
c..   tries courant condition for explicit hd
      do knd = 1, condim
         kkmin(knd)= 0
         tmin(knd) = 0.0d0
      enddo



c..   cond(1) is courant condition (sound travel time over zone)
c..   sound speed c(k) is defined by state.f
c..   "period" is sound travel time over radius segment in grid
      period = 0.0d0
      do k = 2, kk
         if( sound(k) .gt. 0.0d0 )then
            scr(k) = dabs(  ( r(nc,k) - r(nc,k-1) )*0.75d0/sound(k)  )
         else
            write(*,*)'timstp: zero sound speed in zone ',k,sound(k)
            stop'timstp'
         endif
         period = period + (r(nc,k) - r(nc,k-1))/sound(k)
      enddo
      tmin(1) = huge
      kkmin(1) = 0
      do k = 2, kk
         if( tmin(1) .gt. scr(k) )then
            tmin(1) = scr(k)
            kkmin(1) = k
         endif
      enddo
      

c..   cond(2) is compression time from existing velocities
c..   (avoid Lagrange radii crossing!)
c..	revised indexing wda 12-11-08 for zone centering of src::
      ncond = 2
      do k = 2, kk
         if( u(nc,k)-u(nc,k-1) .lt. 0.0d0 )then
c..   compression (target change is 0.02 of thickness)
            scr(k) = -0.02d0*(r(nc,k)-r(nc,k-1) )/(u(nc,k)-u(nc,k-1) )
         else
c..   expansion
            scr(k) = huge
         endif
c     write(*,'(i5,1p8e12.3)')k,scr(k),u(nc,k),r(nc,k),
c     1        (u(nc,k)-u(nc,k-1)),(r(nc,k)-r(nc,k-1))
      enddo
      scr(1) = huge

c..   find minimum estimate of time step
      tmin(ncond)  = scr(2)
      kkmin(ncond) = 2
      do k = 3, kk
         if( tmin(ncond) .gt. scr(k) .and. scr(k) .gt. 0.0d0)then
            tmin(ncond) = scr(k)
            kkmin(ncond) = k
         endif
      enddo
      
      
c..   tmin(3) is the temperature change condition
      do k = 2, kk
         if( dt(2,k) .ne. 0.0d0 )then
            scr(k) = dabs( dth(2)*(t(1,k) + temin)*cdelt/dt(2,k) )
         else
            scr(k) = 0.0d0
         endif
      enddo
c..   find minimum value, set initial value for test
      tmin(3) = scr(2)
      kkmin(3) = 2
      do k = 2, kk
         if( tmin(3) .gt. scr(k) )then
            tmin(3) = scr(k)
            kkmin(3) = k
         endif
      enddo
      
c..   tmin(4) limits change in zone volume
      do k = 2, kk
         if( dv(2,k) .ne. 0.0d0 )then
            scr(k) = dabs( v(1,k)*dth(2)*cdelv/dv(2,k) )
         else
            scr(k) = 0.0d0
         endif
      enddo
      tmin(4)  = huge
      kkmin(4) = 0
      do k = 2, kk
         if( tmin(4) .gt. scr(k) )then
            tmin(4) = scr(k)
            kkmin(4) = k
         endif
      enddo      
      

c..   composition changes: mass fractions
c..   tmin(n+4) is composition change for composition n
c..   nnxdt is number of nuclei tested, see build.f
      do n = 1, nnxdt
c..   nxdt() contains indices in network of chosen nuclei
c..   see build.f for definition
c..   choose nucleus ix to test
         ix = nxdt(n)
         ij = n+4
         do k = 2, kk
c..   fact is mass fraction, shifted down by tiny amount
            fact = ( x(ix,k) - tiny )
            if( fact .gt. 0.0d0 )then
               fak = xd(ix,k)
c..   this reduces sensitivity to small abundance nuclei
               fact = fact + 0.1d0
c..   use both positive and negative changes
               if(  fak .ne. 0.0d0 )then
                  scr(k) = cdeln * dabs( fact/fak )
               else
                  scr(k) = 0.0d0
               endif
            else
               scr(k) = 0.0d0
            endif
         enddo
         tmin(ij) = huge
         kkmin(ij) = 0
         do k = 2, kk
            if( tmin(ij) .gt. scr(k) .and. scr(k) .gt. 0.0d0)then
               tmin(ij) = scr(k)
               kkmin(ij) = k
            endif
         enddo
      enddo

c..   change in stellar mass
c..   tmin(5+nnxdt) limits change in envelope mass to 0.01d0
c..   per time step
      ncond = 5 + nnxdt
      if( mloss .lt. 0 .and. peryear .lt. 0.0d0 )then
c..0.02 implies 35 steps to half the mass in the envelope
c..0.03 implies 23
         tmin(ncond) = - 0.05d0*dmh(kk+1)/(peryear/secpy*sol)
         kkmin(ncond) = kk+1
      else
         tmin(ncond) = 0.0d0
      endif

c..fractional change in luminosity
c..only use in large time step limit (dth(2) > 1.0d10 seconds)
      ncond = 6 + nnxdt
      dlumin = abs(tl(2,kk) - tl(1,kk))
      alumin = dmax1( abs( tl(2,kk) ), abs( tl(1,kk) )  )
      if( dlumin .gt. 0.0d0 .and. nc .eq. 2 .and. 
     1     dth(2) .gt. 1.0d10 )then
         tmin(ncond) = 0.05d0 * alumin / dlumin * dth(2)
         kkmin(ncond) = kk
      else
         tmin(ncond) = huge
         kkmin(ncond) = 0
      endif

c..   restrict increase in time step to a factor of fak
      ncond = 7 + nnxdt
      if( nc .eq. 1 )then
c..no increase on first step (nc = 1 on first call)
         fak = 1.0d0
      elseif( mode .eq. 0 )then
         fak = 1.09d0
      else
         fak = 1.414d0
      endif
      tmin(ncond)  = fak*dth(2)
      kkmin(ncond) = 0

c..   reduce time step if iterations (it) approach iter
      ncond = 8 + nnxdt
      it1 = iter/3 + 1
      it2 = iter/2 + 1
      if( it .gt. it1 )then
c..   large
         if( it .gt. it2 )then
c..   adjust downward
            fact = 1.0d0*dth(2)
         else
c..   no change
            fact = 1.41d0*dth(2)
         endif
      else
c..   few iterations, allow increase
         fact = huge
c         fact = 1.41d0*dth(2)
      endif
      if( fact .gt. 0.0d0 )then
         tmin(ncond) = fact
      endif

c..diffusive mix; see cmix.f and dthoul.f
c..dth limit to keep diffusive scheme stable gives dthc
c..only exceed dthc slightly
      ncond = 9 + nnxdt
c..uses same fak as in timestep size change above
      if( nc .eq. 2 )then
c..increase timestep slowly
c..diffusion is only important for large timesteps
c         tmin(ncond)  = dthc * fak
c..exceeds previous time step by 10 percent
c..if diffusion is not restrictive, dtdiff >> dth(2) (see dthoul.f)
         tmin(ncond)  = dtdiff * 1.05d0
         kkmin(ncond) = kxmx
      else
c..no constraint
         tmin(ncond)  = 0.0d0
         kkmin(ncond) = 0.0d0
      endif
      tmin(ncond) = huge

c..   fractional change in outer radius
      if( modes .eq. 2 )then
         ncond = 11 + nnxdt
         if(  r(2,kk+1) .ne. r(1,kk+1) .and. r(2,kk+1) .gt. 0.0d0 )then
            tmin(ncond) = dth(2)*0.05d0/abs( r(2,kk+1)/r(1,kk+1)-1.0d0 )
            kkmin(ncond) = kk+1
         else
c..new radius not yet defined on first call
            tmin(ncond) = huge
            kkmin(ncond) = 0
         endif
      endif

c..change in metallicity
      ncond = 12 + nnxdt
      do k = 2, kk
         scr(k)  = 0.0d0
         zcno    = 0.0d0
         dlnzcno = 0.0d0
         do n = 1, nnuc
            if( lz(n) .eq. 6 .or. lz(n) .eq. 7 .or. lz(n) .eq. 8 )then
               zcno = zcno + x(n,k)*xa(n)
               dlnzcno = dlnzcno + xd(n,k)*xa(n)
            endif
         enddo
         if( dlnzcno .ne. 0.0d0 )then
            dlnzcno = abs( zcno / dlnzcno )
         else
            dlnzcno = 0.0d0
         endif
         if( nburn .eq. 0 )then
            scr(k) = dlnzcno * cdeln
         else
            scr(k) = dlnzcno 
         endif
      enddo
      tmin(ncond) = huge
      kkmin(ncond) = 0
c      do k = 2, kk
c         if( tmin(ncond) .gt. scr(k) .and. scr(k) .gt. 0.0d0)then
c            tmin(ncond) = scr(k)
c            kkmin(ncond) = k
c         endif
c      enddo    

c..safety check
      do i = 1, ncond
         if( tmin(i) .le. 0.0d0 )then
            tmin(i) = huge
         endif
      enddo

c..   knd is the condition type
c..   tmin(knd) are the minima for each condition type
c..   kkmin(knd) are zones where minimum conditions occurred
c..   tman is the minimum of tmin(knd)  knd=1,2,3, ...
c..   kkman is zone index where tman occurred
c..initialize < 0 to test if any estimates are sane
      kkman = -1

c..   choose most restrictive of the types - save as ktype
c..   arbitrary enormous value
      tman = huge
      do knd = lo, ncond
         if( tmin(knd) .gt. 0.0d0 )then
            if( tmin(knd) .lt. tman )then
               tman  = tmin(knd)
               kkman = kkmin(knd)
               ktype = knd
            endif
         endif
      enddo


      if( kkman .lt. 0 )then
         if( mode .eq. 0 )then
c..   hydrodynamic
            write(3,*)'error, all dt criteria = 0 '
            write(*,*)'error, all dt criteria = 0 '
            stop'timstp: hydro'
         else
c..   hydrostatic
            tman  = tmin(1)
            kkman = kkmin(1)
            ktype = ncond + 1
            write(*,*)' hydrostatic:  use courant cond, others zero'
            write(3,*)' hydrostatic:  use courant cond, others zero'
         endif
      endif

c..   idt .ne. 0 gives constant time step
      if(idt .eq. 0)then
c..   update to new timestep
         dth(1) = dth(2)
         dth(2) = tman
      else
c..   leaves dth(1) and dth(2) unchanged
         ncond =  nnxdt + 10
         ktype  = ncond
         kkmin(ktype) = 0
         tmin(ktype) = dth(2)
      endif


c     SUCCESS
      return
c     
      end


