      

      subroutine rezone(m)

c     
c     REVISION HISTORY
c     04-01-08   wda, no rezone for dtime < 1.0d7 seconds
c     07-30-06
c     
c     rezone adds/deletes an outer henyey zone,
c     forces uniform mass zoning,
c     oversees normal hydrostatic rezoning in mass,
c     does actual interpolation of (2,k) onto new grid (1,k)
      
c
c     NEEDS MAJOR SPRING CLEANING!
c     (C.MEAKIN, 9-29-2008)
      
      
      implicit none

      include 'dimenfile'
      include 'czone'
      include 'cgtr.h'
      include 'cburn'
      include 'cconst'
      include 'comod'
      include 'compu'
      include 'cenv'
      include 'csurface'
      include 'cnabla'
      include 'ceoset'
      include 'cgtintp'

      real*8    scr(kdm)
      real*8    scr0(kdm),scr1(kdm), scr2(kdm), scrx(ndim,kdm)
      real*8    xs(ndim),velm(kdm)
      real*8    scr5(kdm)
      
      real*8    dtiter
      real*8    erro,errsum,errsum0,errlim
      
      integer*4 m,j,jm,jmax,n,leat,nj,jmn,jlim
      integer*4 i,k
      
      real*8    ts,vs,ps,dmhs,dmis,tls,rs,xms,omegs,ajays,rotshears
      real*8    xmp,rb,tb,vb,pb,deltam
      real*8    fqqn,fepslum,fact,dum, xsum
      real*8    fak,faklo,fakhi,difp,difpo,difplo,difphi
      
      real*8    ytot(kdm)
      

c----------------------------------------------------------------
c     ktot .gt. 1 in zone gives normal logic for rezoning
c     ktot .eq. 0 in main gives skip of rezone call
c     ktot .lt. 0 in main following rezone call gives ktot = 0
c     ktot .lt. -1 here gives equal mass zoning (-ktot zones)
c--------------------------------------------------------
c      write(*,*)'entering REZONE'
c--------------------------------------------------------

c..   no rezoning cases:
c..   set njj as flag for no rezoning 
      njj = 0
c..no rezoning at all
c..not adapted for explicit hydro (mode = 0
c..silent return)
      if( mode .eq. 0 )then
         return
      endif
c..no rezoning if step is order of thermal time (roughly)
c      if( ktot .eq. 0 .or. dth(2) .lt. 1.0dr42 )then
      if( ktot .eq. 0 .or. dth(2) .lt. 1.0d-2 )then
c      if( ktot .eq. 0 )then
        if( ktot .gt. 0 )write(*,*)
     1'MSG(REZONE): no rezone due to small time step'
c  1.0d7 was used previously for 1 Msun, which was too large
c  for 50Msun with dM/dt
c  this works and corresponds to several hours
         return
      endif

c..   rezoning will occur
c..   initialize interpolation grids
      do j = 1, kdm
         qqo(j)  = 0
         qq2o(j) = 0
         qqn(j)  = 0
         qq2n(j) = 0
      enddo

      jm   = jj-1
      jmax = jj+1
c..   save boundary values
c..   zone center outside boundary kk(jj)
c..   use updated values (2,k)
      ts   =  t(2,jmax)
      vs   =  v(2,jmax)
      ps   =  p(2,jmax)
      dmhs =  dmh(jmax)
      do n = 1, ndim
         xs(n) = x(n,jmax)
      enddo

c..   zone boundary outside jj
c      tls  = tl(2,jmax)
      rs   =  r(2,jmax)
c..   with correct luminosity and rjoin, fitenv.f will construct the envelope
c..   and radius, so use converged luminosity tl(2,jj) at join
      tls  = tl(2,jj)

      xms  =   xm(jmax)
      omegs= omeg(jmax)
      ajays= ajay(jmax)
      rotshears = rotshear(jmax)
c..   photospheric values
      xmp  =   xm(jmax+1)
      rb   =   r(m,jmax+1)
      tb   =   t(m,jmax+1)
      vb   =   v(m,jmax+1)
      pb   =   p(m,jmax+1)

c..   use epslum to divide flame zones
      epslum = 0
      do j = 2, jj
         epslum = epslum + dmh(j)*s(5,j)
      enddo

c..   Will we move the join point in mass coordinate?
      if( modes .eq. 2 )then
c..   default is do nothing:         leat = 0
c..   move join inward:  leat = -1  njoin = -1
c..   move join outward: leat = +1  njoin = +1

c..   constraints on changes/match at join
         dum = abs( tl(2,kk+1) - tl(2,kk) )
     1        / (abs( tl(2,kk+1) )+abs( tl(2,kk) ) )
         fact =  abs( rarget - r(2,kk) )
     1        / (rarget + r(2,kk) )
            
         if( dum .gt. 0.02d0 .or. fact .gt. 0.10d0 
     1      .or. dth(2) .lt. 1.0d4)then
            write(*,*)'NO SURFACE REZONE: dL/dk ',dum,' dRin/dk ',fact
            leat = 0
         else
c..   total envelope mass constraint fmnenv
            if( dmh(jj+1)/(dmh(jj+1)+xm(jj)) .lt. fmnenv )then
               leat = -1
            elseif( dmh(jj+1)/(dmh(jj+1)+xm(jj)) .gt. fmxenv )then
               leat = 1
            else
               leat = 0
            endif
            
c..   override, temperature (burning) takes priority
c            if( t(1,jj) .gt. 2.2d6 .and. dmh(jj+1).gt.dmh(jj)
            if( t(1,jj) .gt. 0.5d6 .and. dmh(jj+1).gt.dmh(jj)
     1           .and. dmh(jj+1)/xm(jj) .gt. 1.0d-3 )then
c..   except for small envelope fractions
c..   too hot, move join out (lithium 7)
c..   solar Li7 is near 0.9e-8 by mass, so this is down by 1/100 or so
               if( x(lli7,jj) .gt. 0.9d-10 )then
c..   use lower value for Li7 because of its use as tracer
                  leat = 1
               endif
            endif
c            if( t(1,jj) .gt. 0.8d6 .and. dmh(jj+1).gt.dmh(jj)
            if( t(1,jj) .gt. 0.5d6 .and. dmh(jj+1).gt.dmh(jj)
     1           .and. dmh(jj+1)/xm(jj) .gt. 1.0d-3 )then
c..   too hot, envelop mass exceeds grid zone,
c..   envelope fraction is not tiny, so move join out (deuterium)
               if( x(ldeut,jj) .gt. 1.0d-8 )then
                  leat = 1
               endif
            endif

c..   override, avoid iterating on Hydrogen partial ionization regions 
c..   if join temperature is smaller the code jumps back and forth
c..   between solutions, reducing the timestep
c..   override, mass loss takes priority
c..   tmenv is minimum join temperature; set in params.d (gen.f)

            if( -peryear*sol*dth(2)/secpy .gt. 1.0d-2*dmh(jj+1)
     2           .or. t(1,jj) .lt. tmenv )then

c..   too much decrease in envelope mass, move boundary in,
c..   until join temperature is near network boundary tmenv
c..
c..   or, join temperature is too low, so that 
c..   "thermal pulses" due to partial ionization/opacity decrease
c..   can give bad convergence for hydrostatic models
               if( t(1,jj) .lt. 0.95d7 )then
c..   move in leat zones, up to a maximum of 4
                  do i = 1, 4
                     leat = -i
c..   conditions for exiting (stop increasing envelope mass)
c..   keep envelope to a fraction of its original mass
                     if( xm(jj+1)-xm(jj-i) .gt. fmxenv*xm(jj+1) .or.
c..   mass lost in this step is small compared to new envelope mass
     1                    -peryear*sol*dth(2)/secpy .lt. 
     2                    0.01d0*(xm(jj+1)-xm(jj-i)) .and.
c..   limit to cooler regions to avoid burning (Li, D above have already 
c..   been tested)
c..   temperature is above the partial ionization region
     4                    t(1,jj+1-i) .gt. tmenv )then
                        goto 190
                     endif
                  enddo
                  if( -peryear*dth(2)/secpy .gt. 1.0d-2* dmh(jj+1)/sol 
     1                 )then
                     write(*,'(a35,1pe12.3,a22)')
     1                    'REZONE warning: strong mass loss',
     2                    -peryear*dth(2)/secpy,
     3                    'M(sol) this timestep'
                     write(*,'(a35,1pe12.3,a8)')
     1                    'REZONE warning: envelope mass is',
     2                    dmh(jj+1)/sol,
     3                    'M(sol)'
                  endif
                  write(*,'(a35,i12,a29)')
     1                 'REZONE warning: wants more than ',
     1                 -leat,' zones; will continue anyway'
 190              continue

                  write(*,'(a5,i5, (3(a12,1pe12.3)) )')
     1                 'leat',leat,'new Menv',xm(jj+1)-xm(jj+leat),
     2                 'M lost',-peryear*sol*dth(2)/secpy,
     3                 'old Menv',dmh(jj+1)
               endif
            endif
         endif
      else
         leat = 0
      endif
c..logic for join motion is finished, leat defined

c..revise henyey arrays for new join position
      if( leat .le. -1 )then
         write(*,'(a30,2i5,4(a5,1pe13.5))')
     1        'REZONE: RESET JOIN FURTHER IN',
     2        jj+leat,-leat,'BY',(xm(jj)-xm(jj+leat))/sol,
     3        'TO',xm(jj+leat)/sol,'FROM',xm(jj)/sol,
     4        'ENV',dmh(jj+1)/sol
c         write(3,'(a30,2i5,4(a5,1pe13.5))')
c     1        'REZONE: RESET JOIN FURTHER IN',
c     1        jj+leat,-leat,'BY',(xm(jj)-xm(jj+leat))/sol,
c     2        'TO',xm(jj+leat)/sol,'FROM',xm(jj)/sol,
c     4        'ENV',dmh(jj+1)/sol

c---- new zoning
c.... njj........njj+1.njj+2
c.....|.....x.....|..x..|

c---- old zoning
c.....|..x..|..x..|..x..|
c.... jj-1..jj...jj+1..jj+2
c.......jj....jj+1..jj+2

c---- fitenv zoning
c.....|..|..|..|..|..
c..   jmaxz -1 -2 -3 -4

c..   limits for new zoning (eat in -leat zones)
         njj   = jj + leat

c..   revise boundary values
         nj    = njj+1
c..   these correspond to -leat zones inward from jj+1
         ts    =  t(m,nj)
         vs    =  v(m,nj)
         ps    =  p(m,nj)
c..   envelope mass is increased
         dmhs = dmh(jj+1)
c..   WARNING: the envelope moment of inertia should come from fitenv.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         dmis = dmh(jj+1)
         do i = 1, -leat
            j = jj + 1 - i
            dmhs = dmhs + dmh(j)
            dmis = dmis + dmi(j)
         enddo

c..   abundances preserved by rezoning: old zone abundances added
c..   to envelope 
         do n = 1,ndim
            xs(n) = (x(n,jj)*dmh(jj) + x(n,jj+1)*dmh(jj+1))/
     1           (dmh(jj) + dmh(jj+1))
         enddo
         do n = 1, ndim
            x(n,jj+1) = xs(n)
         enddo
c.. 
      rs    =  r(m,jj+1)
      tls   = tl(m,jj+1)
      xms   =   xm(jj+1)
      omegs = omeg(jj+1)

c..   angular momentum is conserved
c..   old envelope values
         do n = 1, ndim
            xs(n) = x(n,jj+1)*dmh(jj+1)
         enddo
         ajays = ajay(jj+1)*dmi(jj+1)

c..   mix newly eaten zones into envelope
         do i = 1, -leat
            do n = 1, ndim
               xs(n) = xs(n) + x(n,jj+1-i)*dmh(jj+1-i)
            enddo
c            ajays = ajays + ajay(jj+1-i)*dmi(jj+1-i)
            ajays = ajays + ajay(jj+1-i)*dmi(jj+1-i)*
     1              (r(m,jj+1-i)/r(m,jj+1))**2.0d0
         enddo

c..   rescale to per mass unit
         do n = 1, ndim
            xs(n) = xs(n) / dmhs
         enddo
         ajays = ajays / dmis
         omegs = ajays/r(m,njj)**2.0d0
c..surface values saved so adjust outer zone+1 = envelope
         xm(njj)     = xm(kk+leat)
         xm(njj+1)   = xms
         dmh(njj+1)  = xm(njj+1) - xm(njj)
         dmi(njj+1)  = dmh(njj)
         r(m,njj+1)  = rs
         tl(m,njj+1) = tls
         ajay(njj+1) = ajays
         omeg(njj+1) = omegs
         rotshear(njj+1) = omegs*rs

c..zero array as grid recede
c..zero array as grid recedes to smaller jj (kk+leat)
         dmh(njj+2) = dmh(jj+2)
         dmh(njj+3) = 0.0d0
         do j = njj+2,njj+4
            xm(j) = 0.0d0
         enddo

c..   define primary array 
         qqo(1) = xm(1)
         do j = 2, njj+1
            qqo(j) = qqo(j-1) + dmh(j)
         enddo

c..   join boundary moves inward
c..   innermost mass coordinates are unchanged (njj = jj-(-leat) )
         do j = 1, njj-1
            qqn(j) = qqo(j)
         enddo
         qqn(njj) = xm(jj+leat)

c..reset for gtintp
         jj = njj
         kk = jj

c........................
      elseif( leat .eq. 1 )then
c..........................
         write(*,'(a30,2i5,4(a5,1pe13.5))')
     1        'REZONE: RESET JOIN FURTHER OUT',
     1        jj,jj+1,'BY',dmh(jj)/sol,
     2        'TO',(xm(jj)+dmh(jj))/sol,'FROM',xm(jj)/sol,
     4        'ENV',(dmh(jj+1)-dmh(jj))/sol
c         write(3,'(a30,2i5,4(a5,1pe13.5))')
c     1        'REZONE: RESET JOIN FURTHER OUT',
c     1        jj,jj+1,'BY',dmh(jj)/sol,
c     2        'TO',(xm(jj)+dmh(jj))/sol,'FROM',xm(jj)/sol,
c     4        'ENV',(dmh(jj+1)-dmh(jj))/sol

c---- old zoning
c.....jj.........jj+1..jj+2
c.....|.....x.....|..x..|

c---- new zoning
c.....|..x..|..x..|..x..|
c...  njj-1..njj...njj+1..njj+2
c.......njj...njj+1..njj+2

c---- fitenv zoning
c.....|..|..|..|..|..
c..   jmaxz -1 -2 -3 -4

c..   limits for new zoning (add one zone by extrapolation)
         njj   = jj+1

c..   move outer zone values one zone outward, highest first(njj+2)
         r(m,njj+2) = r(m,jj+2)
         tl(m,njj+2)= tl(m,jj)
         xm(njj+2)  = xm(jj+2)
         dmi(njj+2) = dmi(jj+2)
c..   zone center quantities (njj+2)
         dmh(njj+2) = dmh(jj+2)
         t(m,njj+2) = t(m,jj+2)
         v(m,njj+2) = v(m,jj+2)
         p(m,njj+2) = p(m,jj+2)
         do n = 1, ndim
            x(n,njj+2) = x(n,jj+2)
         enddo

c..   move next zone values one zone outward(njj+1)
         r(m,njj+1) = r(m,jj+1)
         tl(m,njj+1)= tl(m,jj)
         xm(njj+1)  = xm(jj+1)

c..   use converged value for luminosity
         tl(m,njj)  = tl(m,jj)
c..   extrapolation uses equal mass increments, so
         dmh(njj)   = dmh(jj)
c..   move new surface mass coordinate out accordingly
         xm(njj)    = xm(jj) + dmh(njj)
c..   decrement envelope mass 
         dmh(njj+1) = xm(njj+1) - xm(njj)

c..   simply extrapolate one step (shares error)
c..   boundaries enclose equal mass zones

c..   zone centers may correspond to different mass zones
         fak = 2.0d0*dmh(jj)/( dmh(jj) + dmh(jj-1) )
         p(m,njj) = (1.0d0 + fak)*p(m,jj)-p(m,jj-1)*fak
         t(m,njj) = (1.0d0 + fak)*t(m,jj)-t(m,jj-1)*fak
c..same extrapolation as P and T
         v(m,njj) = 1.0d0/( (1.0d0+fak)/v(m,jj) -fak/v(m,jj-1))
c..mass conservation implies new radius
         r(m,njj) = ( r(m,njj-1)**3 + v(m,njj)/pi43*dmh(njj) 
     1        )**(1.0d0/3.0d0)

         if( modes .eq. 2 )then
c..   for consistency with fitenv and hstat
            dmi(njj+1) = dmh(njj)
            dmi(njj) = dmh(njj)
         else
            dmi(njj+1) = 0.5d0*( dmh(njj) + dmh(njj+1) )
         endif

c..   abundances preserved: new zone has envelope abundances
         do n = 1,ndim
            x(n,njj)   = xs(n)
            x(n,njj+1) = xs(n)
         enddo

c..   new zone has specific angular momentum of envelope 
         ajays = ajay(jj+1)
c         rotshears = rotshear(jj+1)
         omegs = ajays/(r(m,jj+1)**2.0d0)
         rotshears= omegs*r(m,jj+1)
c..   new zone outside njj
         ts     = t(m,njj+1)
         ps     = p(m,njj+1)
         vs     = v(m,njj+1)
         dmhs   = dmh(njj+1)
c..   zone boundary outside njj
         tls = tl(m,njj+1)
         rs  = r(m,njj+1)
         xms = xm(njj+1)
c         omegs = ajays*(r(m,jj+2)/rs)**2.0d0
c..   photospheric values
         xmp = xm(njj+2)
         rb  = r(m,njj+2)
         tb  = t(m,njj+2)
         vb  = v(m,njj+2)
         pb  = p(m,njj+2)
         
c..   define primary array
         qqo(1) = xm(1)
         do j = 2, njj+1
            qqo(j) = qqo(j-1) + dmh(j)
         enddo

c..   model adjusted, reset for gtintp
         jj = njj
         kk = njj
c.......................
      else
c.......................

         njj   = jj
c..   define primary array 
         qqo(1) = xm(1)
         do j = 2, jmax
            qqo(j) = qqo(j-1) + dmh(j)
         enddo

c..   same zone number as yet
c      jj = njj
c.......................
      endif
c.......................
c..the arrays are now adjusted out to boundary jj 
c..in old mass coordinates qqo

c----------------------------------------------------------------------
      if( ismoo .lt. 0 )then
c----------------------------------------------------------------------

c..   rezone with -kk uniform mass zones
         write(*,*)'REZONE: ismoo = ',ismoo
         if( -ismoo .lt. kk/2 )then
            write(*,*)'too few zones ', -ismoo,' old value ',kk
            stop'rezone ismoo'
         endif
         kk = -ismoo
         write(*,*)'REZONE: ',kk,' equal mass zones'
         njj    = kk - 1
         jmn    =  njj  - 1
         deltam = ( xm(jj) - xm(1) )/dble( njj - 1 )
         write(*,'(a20,1p8e12.3)')'zone mass',deltam

c..   redefine primary array
         qqn(1) = xm(1)
         do j = 2, njj
            qqn(j) = qqn(j-1) + deltam
         enddo
c..   consequently:
         do j = 1,jm
            qq2o(j+1) = 0.5d0 * (qqo(j) + qqo(j+1) )
         enddo

         do  j = 1,jmn
            qq2n(j+1) = 0.5d0 * (qqn(j) + qqn(j+1) )
         enddo

c..   reset flag for no automatic rezoning to allow relaxation
         ktot = 0

         write(3,11) njj+1, deltam
         write(6,11) njj+1, deltam
 11      format(i5,1x,'Equal mass rezoning: zonesize =',1pe13.5)

c----------------------------------------------------------------------
      elseif( ismoo .eq. 0 )then
c----------------------------------------------------------------------

         if( mode .ne. 0 )then
c..   gtintp defines new grid in mass, qqn(1:njj), for interpolation

            call gtintp
            
         endif

      else
c,,,,,,.see gtintp.f for smoothing as well???

c..     ismoo positive
c     smooths zone sizes by relaxation (critical diffusion)

       dqqn(1) = dmh(1)
       dqqn(2) = dmh(2)
c.. inner boundary
       dqqn(3) = 0.25d0*dmh(4) + 0.75d0*dmh(3)
       do j = 4, jj-2
          dqqn(j) = 0.25d0*dmh(j+1) + 0.5d0*dmh(j) 
     1      +0.25d0*dmh(j-1)
       enddo
       dqqn(jj-1) = 0.25d0*dmh(jj-2) + 0.75d0*dmh(jj-1)
       dqqn(jj) = dmh(jj)
       dqqn(jj+1) = dmh(jj+1)

       qqn(1) = qqo(1)
       do j = 2, jj+1
         qqn(j) = qqn(j-1) + dqqn(j)
       enddo
       
c..   consequently:
      do j = 1,jj-1
        qq2o(j+1) = 0.5d0 * (qqo(j) + qqo(j+1) )
      enddo

      do  j = 1,jj-1
        qq2n(j+1) = 0.5d0 * (qqn(j) + qqn(j+1) )
      enddo
      qq2o(1) = qq2o(2)
      qq2n(1) = qq2n(2)

c----------------------------------------------------------------------
      endif
c----------------------------------------------------------------------
c..   ready to interpolate and update variables

      do j = 1, jj+1
         ytot(j) = 0.0d0
         do n = 1, nnuc+1
            ytot(j) = ytot(j) + x(n,j)
         enddo
      enddo

c..   redefine masses from new grid
      do j = 1, njj
         xm(j) = qqn(j)
      enddo

c..   round-off ??
      do j = 2, njj
         dmh(j) = xm(j) - xm(j-1)
      enddo
      dmh(1) = dmh(2)

c..   uses new value at outer zone for ghostzone
      dmh(njj+1) = dmhs

      do j = 1, njj
         dmi(j) = 0.5d0*( dmh(j) + dmh(j+1) )
      enddo
c..   consistency with fitenv and hstat for modes=2
      if( modes .eq. 2 )then
         dmi(njj+1) = dmh(njj)
         dmi(njj)   = dmh(njj)
      else
         dmi(njj+1) = dmh(njj+1)
      endif
c..   abundances preserved: new zone has envelope abundances
         do n = 1,ndim
            x(n,njj)   = xs(n)
            x(n,njj+1) = xs(n)
         enddo
c..   reset iteration count
         it = 1

c..   interpolation for new values of secondary arrays
c..   use only linear interpolation (lintr subroutines, not aintr)
c..   to avoid Gibbs ringing
      call lintr2(qqo,qqn,tl,scr,jj,njj,m)

      call lintr2(qqo,qqn, u,scr,jj,njj,m)
      do j = 1, njj
         dr(j) =  u(m,j)*dth(m)
      enddo

c..   use linear interpolation for h.ge. 0
      call lintr1(qqo,qqn,  h,scr,jj,njj)

      do j = 1, njj
         if( ic(j) .eq. 0 .and. h(j) .gt. 1.0d0 )then
            ic(j) = 1
         endif
      enddo

      call lintr1(qqo,qqn,  y,scr,jj,njj)
      call lintr1(qqo,qqn, du,scr,jj,njj)
      call lintr1(qqo,qqn,gam,scr,jj,njj)
      call lintr1(qqo,qqn,dtl,scr,jj,njj)

      call lintr1(qqo,qqn,dnab,scr,jj,njj)
      call lintr1(qqo,qqn,dnad,scr,jj,njj)


c..   rotational variables (defined at boundary here!)
      call lintr1(qqo,qqn,omeg,scr,jj,njj)
      call lintr1(qqo,qqn,ajay,scr,jj,njj)
      call lintr1(qqo,qqn,rotshear,scr,jj,njj)
c..   zone center quantities
c..   use linear interpolation to avoid Gibbs wiggles at steep gradient

c..   interpolate abundances
c..   save initial values
      do j = 1, jj+1
         do n = 1, nnuc+1
            xold(n,j) = x(n,j)
         enddo
      enddo

      do j = 1, jj+1
         scr5(j) = doux(j)
      enddo

c..   lintrx modified to give a linear interpolation of abundances (2/2/06)
      call lintrx(qqo,qqn,x,scrx,jj,njj)

c..energy generation rates
      call lintr1(qq2o,qq2n,ss,scr,jj,njj)
      call lintr1(qq2o,qq2n,snu,scr,jj,njj)
      do j = 2, njj
         s(5,j) = ss(j)
         s(4,j) = snu(j)
         s(6,j) = ss(j) + snu(j)
         s(3,j) = s(6,j)
      enddo

c..   tiny errors from interpolation can mess up network solutions
c..   make x positive or zero
      do i = 2, njj+1
         do j = 1, ndim
            x(j,i) = dmax1( x(j,i), 0.0d0 )
         enddo
      enddo
c..   renormalize and reset Ye
      do i = 2, njj+1
         xsum = 0.0d0
         do j = 1, ndim-1
            xsum = xsum + x(j,i)
         enddo
         do j = 1, ndim-1
            x(j,i) = x(j,i)/xsum
         enddo
         x(ndim,i) = 0.0d0
         do j = 1, ndim-1
            x(ndim,i) = x(ndim,i) + x(j,i)/xa(j)*dble( lz(j) )
         enddo
      enddo

c..   mole number of particles
      do i = 2, njj+1
         ytot(i) = 0.0d0
         do j = 1, ndim
            ytot(i) = ytot(i) + x(j,i)/xa(j)
         enddo
      enddo

c..interpolate for specific volume v
      call lintr2(qq2o,qq2n,v,scr,jj,njj,m)

c..   get r from dmh and v : gtr version
c..   new inner radius is old r(m,1)
c..   this uses nucleon conservation
      do j = 2, njj
         fact   = dmh(j)*v(m,j)/pi43*( gam(j) + gam(j-1) )/2.0d0
         r(m,j) = (r(m,j-1)**3 + fact )**(1.0d0/3.0d0)
      enddo

      do k = 2, njj
         a(k) = pi4*r(m,k)**2
         g(k) = grav*xm(k)/r(m,k)**2
      enddo


c..   abundances may be changed by iterpolation above
c..   mass conservation has already been applied by V(r) calculation
c..   use pressure rather than temperature to get hydrostatic equilibrium

c..   interpolate for Temperature t to get first guess for iteration
c..   in t(m,k)
      call lintr2(qq2o,qq2n,t,scr,jj,njj,m)
      t(m,njj+1) = ts

c..   interpolate for pressure P; desired value in p(m,k) and src(k)
      call lintr2(qq2o,qq2n,p,scr,jj,njj,m)
      p(m,njj+1) = ps

c..   state needs kk 
c..   update index limit on henyey grid
      kk = njj

c..   fak is criterion for pressure convergence
      if( ismoo .lt. 0 )then
         fak = 1.0d-5
      else
         fak = 1.0d-8
      endif

c..jlim is j index for debug output for slow convergence
      jlim = 80
      do k = 2, kk
         do j = 1, 100
c..   state gets P(m,k) for this T,V pair
            call state(k,k,m)
c..   scr(k) is desired pressure for this V and composition
            dtiter = (scr(k)-p(m,k))/pt(k)
c..   iterated Pressure    to desired limit of accuracy
c..   iterated Temperature to desired limit of accuracy

            if( abs(dtiter) .lt. fak*t(m,k) .and.
     1           abs(scr(k)-p(m,k)) .lt. fak*p(m,k) )goto 300
            if( j .gt. jlim )then
               write(*,'(a12,2i5,1p3e10.2,1p4e15.7
     1              )')'rezone: P',k,j,
     1              dtiter,dtiter/t(m,k),scr(k)/p(m,k)-1.0d0,
     2          t(m,k),scr(k),p(m,k),pt(k)
            endif
c..slow convergence, so avoid oscillation
            if( j .ge. 80 )dtiter = 0.5d0*dtiter
            t(m,k) = t(m,k) + dtiter
         enddo

         write(*,*)'REZONE P(T) iteration error, j = ',j

               write(*,'(a12,2i5,1p3e10.2,1p4e15.7
     1              )')'rezone: P',k,j,
     1              dtiter,dtiter/t(m,k),scr(k)/p(m,k)-1.0d0,
     2          t(m,k),scr(k),scr(k)-p(m,k),pt(k)
         write(*,'(12x,2a5,3a10,4a15)')'k','j','dT','dlnT','dlnP',
     1        'T(m,k)', 'P(m,k)','scr-P','PT'
         write(*,*)'kk ',kk
         stop'rezone.f error P(T)'

 300     continue
         if( j .gt. jlim )then
               write(*,'(a12,2i5,1p3e10.2,1p4e15.7
     1              )')'rezone: P it',k,j,
     1              dtiter,dtiter/t(m,k),scr(k)/p(m,k)-1.0d0,
     2          t(m,k),scr(k),scr(k)-p(m,k),pt(k)
         endif
      enddo



      errlim = 1.0d-3
c..   evaluate HSE error
      do j = 3, njj-1
c..   scr(j) is error at boundary j in HSE
         scr(j) = ( p(m,j) - p(m,j+1)
     1        - grav*xm(j)*dmi(j)/(pi4*r(m,j)**4) )
     2        *2.0d0/(p(m,j)+p(m,j+1))
      enddo
      scr(2) = 0.0d0
      scr(njj) = 0.0d0

      errsum = 0.0d0
      do j = 2, njj-1
c..find worst zones
         erro = scr(j)
         errsum = errsum + erro**2
      enddo
      errsum = sqrt(errsum)/dble(njj-2)
      errsum0 = errsum
      if( errsum .gt. 1.0d-5 )
     1     write(*,'(a20,1pe12.3)')'rms error ini ',errsum



      do k = 2, njj
         dt(m,k) = 0.0d0
         dv(n,k) = 0.0d0
      enddo
c set dt dv to zero rather than interpolate
c      call lintr2(qq2o,qq2n,dt,scr,jj,njj,m)
c      call lintr2(qq2o,qq2n,dv,scr,jj,njj,m)

c..   set old boundary values at new index..................
c..   zone center outside njj, new value of jj(kk)
      t( m,njj+1)  = ts
      v( m,njj+1)  = vs
      p( m,njj+1)  = ps
      dmh(njj+1)   = dmhs
      do n = 1, ndim
         x(n,njj+1)   = xs(n)
      enddo
c..   zone boundary outside njj
      tl(m,njj+1)  = tls
      r( m,njj+1)  = rs
      xm(njj+1)    = xms
      omeg(njj+1)  = omegs
      ajay(njj+1)  = ajays
      rotshear(njj+1) = rotshears
c..   photospheric values
      xm(njj+2)  = xmp
      r(m,njj+2) = rb
      tl(m,njj+2)  = tls
      t(m,njj+2) = tb
      v(m,njj+2) = vb
      p(m,njj+2) = pb
c      omeg(njj+2) = omegs
c      ajay(njj+2) = omegs*r(m,njj+2)**2.0d0
      do n = 1, ndim
         x(n,njj+1+1) = xs(n)
      enddo
c..   reset changes
      dt(m,njj+1)  = 0
      dv(m,njj+1)  = 0
      dtl( njj+1)  = 0
      dr(  njj+1)  = 0
      du(  njj+1)  = 0

      n = 0
      do k = 2, kk
         fak = p(m,k+1)-p(m,k)
         if( fak .gt. 0.0d0 )then
            write(*,'(a25,3i5,1p8e12.3)')'REZONE: Finish delp',
     1           k,m,ic(k),fak,p(m,k),p(m,k+1),p(2,k+1)-p(2,k),p(2,k),
     2           p(2,k+1),xa(nnuc)*x(nnuc,k),xa(nnuc)*x(nnuc,k+1)
  
            n = n + 1
         endif
      enddo
      if( n .ne. 0 )then
         write(*,'(25x,3a5,8a12)')'k','m','ic','dp','p(m,k)','p(m,k+1)',
     1        'dp(2)','p(2,k)','p(2,k+1)','Xhe4(k)','Xhe4(k+1)'
         write(*,'(/a25,i5,1p8e12.3)')
     1        'REZONE: P inversions',n
c         stop'P inversion  in REZONE'
      endif
      
      
      
c     SUCCESS
      return
     
      end
      













