      subroutine burn(nc,mbflag,idebug,kin,ktop)
c..   1-29-08
c..   nuclear burning with path integration
c..   uses mass fractions

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'cgtr.h'
      include 'cconst'
      include 'cburn'
      include 'caeps'
      include 'cnabla'
      include 'ceoset'
      include 'crate'

c..   nburn is passed in common comod
c..   nburn < 0  gives no burning
c..   nburn = 0  gives path integration
c..   nburn = 1  gives explicit one step
c..   nburn = 2  gives operator split
c..
c..   nc = 1 is old model, nc = 2 is new model index in time
c..
c..   aeps(k) is the array for energy generation (erg/g/sec)-->s(5,k)
c..
c..   mbflag = 0  gives aeps(k)-->sa(k) and aex(n,k) --> x(n,k),
c..   an energy generation and an abundance update,
c..   and an estimate of sa(k)=dss(k)/dT and sb(k)=dss(k)/dV
c..
c..   mbflag = 1  gives only aeps(k)-->ss(k) and no finite difference
c..   estimate of sa(k)=dss(k)/dT and sb(k)=dss(k)/dV
c..
c..   idebug gives verbose diagnostic messages
c..   kin is innner zone number on burn loops
c..   kout is outer zone number 
c..
c..   xold(n,k) is array for abundances from the last model
c..   xold(n,k) is defined in tycho.f before cmix.f is called
c..   x(n,k)    is the array for mixed abundances from cmix.f (xold-->x)
c..   x(n,k)    is the array for replacing old abundances in tycho.f
c..   xd(n,k)   is the time derivative for purely nuclear changes
c..   x0(n,k)   is the initial abundance array (local) fed to burn.f
c..   aex(n,k)  is the array defining abundances for solven.f (x-->aex)
c..   solven.f updates these to its internal array y(n) (aex-->y) for each k
c..
c..   if leq = 0 then new values are saved y(n) --> aex(n,k)

      real*8 x0(ndim,kdm)

c..   itbrn(k) is the temperature index for the solven mode, 
c..   finite derivative estimate (for sanity at T boundaries)
c..   (0 = decay, 1 = low T, 2 = full network, 3=qse and nse)
c..
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..   ieburn(k) is the burning index for the solven code
c..   defines updates in iterative loop (it > 1)
c..   ieburn(k) = 0 => no update, explicit step; 
c..   ieburn(k) = 1 => aeps-->s(5,k), y--> aex, iterated step)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c..   nuck(n) is the zone k of most rapid depletion for nucleus n
c..   xdep(n) is value of that change

c..   working vectors for finite derivatives
      real*8    s4(kdm),s5(kdm)
      real*8    s4t(kdm),s4v(kdm),s5t(kdm),s5v(kdm),s51(kdm),s52(kdm)
      real*8    s5expt(kdm),s5expv(kdm)

      real*8    dxnucl,dxnuclmin
c..   sumex needed for inverse mapping of mass vs nucleon fractions
      real*8    sumex

      integer*4 nxnucl,kxnucl

      integer*4 idebug,nc,mbflag,kin,ktop, leq
      integer*4 k,n,m
      integer*4 inegs

      real*8 dtn, cfak, s1, s2, t9, rho, fact
      parameter( cfak = 1.01d0 )

      real*8 s5max
      integer*4 ks5max, kncyc, ncycmx

      real*8 dtc(kdm)
      real*8 sumx,sumdx

      integer*4 klast

      save
c----------------------------------------------------------------
c      write(*,*)'ENTERING BURN ',x(nnuc,2),xold(nnuc,2),
c     1     -x(nnuc,2)+xold(nnuc,2),nc,kk
ccccccccccccccccccccccccccccccccccccccc

c..   time this operation and accumulate it in runburn

      call seconds(runburn0)

c..   burn and sneut turned off
      if( nburn .eq. -1 .and. neutro .eq. -1 )return
      if( nburn .gt. 1 .or. nburn .lt. -1 )then
         write(*,*)'Only options -1, 0, 1 now implemented in burn.f'
         write(*,*)'nburn = ',nburn
         stop'Burn: nburn error'
      endif
      if( neutro .gt. 0 .or. neutro .lt. -1 )then
         write(*,*)'Only options -1, 0 now implemented in burn.f'
         write(*,*)'neutro = ',neutro
         stop'Burn: neutro error'
      endif

      if( it .le. 1 )then
c..   first iteration in hydro.f
c..   initialize
c..   hide some global arrays from solven.f by using include 'caeps' 
c..   initialize dtc(k) to default time step size
         do k = 1, kk
            dtc(k) = dth(2)
            ncyc(k) = 0
            xidk(k) = '     '
            bk(k) = 0.0d0
            yk(k) = 0.0d0
         enddo

c..   setup for path integration..................................
         if( nburn .eq. 0 )then
c..   radiative-convective boundaries from cinit.f
c..   define time for each path segment
            if( nrczones .ge. 1 )then
               do n = 1, nrczones
c..   grid zones only kk
                  klast = min(kend(n),kk)
                  if( kbeg(n) .lt. kk-1 )then
                     do k = kbeg(n),klast
                        dtc(k) = dth(2)*dmh(k)
     1                       /( xm(klast) - xm(kbeg(n)) )
                     enddo
                  endif
               enddo
            endif
         endif
c..   end of path integration setup...............................

c..   total number of subcycles per time step
         ncytot = 0

         do k = kin, ktop+1
c..   envelope zone too
            ss(k)  = 0.0d0
            snu(k) = 0.0d0
            s(4,k) = 0.0d0
            s(5,k) = 0.0d0
            s(7,k) = 0.0d0
            a(k)   = pi4*r(nc,k)**2
            st(k)  = 0.0d0
            sv(k)  = 0.0d0
            sa(k)  = 0.0d0
            sb(k)  = 0.0d0
            snua(k)= 0.0d0
            snub(k)= 0.0d0

            epsnuc(k) = 0.0d0

            do n = 1, ndim
               xd( n,k) = 0.0d0
               x0(n,k)  = x(n,k)
               axdcon(n,k) = 0.0d0
cccccccccccccccccccc
            enddo

            aeps(k)     = 0.0d0
            aepst(k)    = 0.0d0
            aepsv(k)    = 0.0d0
            aenucnu(k)  = 0.0d0

c..   choose network for finite difference evaluation
c..   uses anchor value tem(k) not finite difference value temp(k)
            t9 = tem(k)*1.0d-9
            if( t9 .lt. tburnd )then
c..   decay only
               itbrn(k) = 0
            elseif( t9 .lt. tburnlo .and. t9 .ge. tburnd )then
c..   loburn
               itbrn(k) = 1
            elseif( t9 .lt. tnse .and. t9 .ge. tburnlo )then
c..   full network
               itbrn(k) = 2
            elseif( t9 .ge. tnse )then
c..   qse and nse
               itbrn(k) = 3
            endif
         enddo

         if( nburn .ne. -1 )then
c..   network burn
c..   initial, obligatory network solution, energy generation
            s5max  = 0.0d0
            ks5max = 0
            do k = kin, ktop
c..   must account for mass excesses to get nucleon number*amu
c..   this was done in state.f and passed in include 'ceoset.h' 
c..   rhoz = nucleon number/ avagadro
c..   rhom = mass / unit volume
               temp(k) = tem(k) 
               den(k)  = rhoz(k)
c..   initial abundances for solven
               aex(ndim,k) = x0(ndim,k)
               do n = 1, ndim-1
c..moles from mass fractions
                  aex(n,k) = x0(n,k)/xa(n)
               enddo

c..   leq = 0 gives update of abundances y(n)-->aex(n,k)
               if( mset .eq. 0 )then
                  leq = 0
               else
c..   force no abundance update (for testing)
                  leq = 1
               endif

               call solven(dtc(k),dtn,k,
     1              idebug,ncytest,ncymax,inegs,it,kk,leq,jnb)

c..   save values for unperturbed T, V
               s5(k) = aeps(k)
c..   energy generation rate
               ss(k) = aeps(k)
c..   nuclear neutrino/antineutrino emissivity is -aenucnu(k) in erg/g/s
               epsnuc(k) = aenucnu(k)

c..   renormalize to new mass for mass fraction
               sumex = 0
               do n = 1, ndim-1
                  sumex = sumex + aex(n,k)*xa(n)
               enddo
               do n = 1, ndim-1
                  aex(n,k) = aex(n,k)/sumex
               enddo
               sumex = 0
               do n = 1, ndim-1
                  sumex = sumex + aex(n,k)*xa(n)
               enddo
               if( abs(sumex-1.0d0) .gt. 1.0d-8 )
     1              write(*,*)'burn ',sumex-1.0d0

               x(ndim,k) = 0.0d0

               do n = 1, ndim-1
                  x(n,k) = aex(n,k)*xa(n)
                  x(ndim,k) = x(ndim,k) + aex(n,k)*dble(lz(n))
               enddo

               if( abs(s5max)*dmh(ks5max) .lt. abs(ss(k))*dmh(k) )then
                  s5max  = ss(k)
                  ks5max = k
               endif
            enddo

c...........................................................
            if( nburn .eq. 0 )then
c               write(*,*)'PATH INTEGRAL NUCLEOSYNTHESIS'
c..   path integral nucleosynthesis
               if( nrczones .ge. 1 )then
c..   there are some convection zones
                  do m = 1, nrczones
c..   use grid zones only
                     klast = min(kend(m),kk)
                     if( kbeg(m) .lt. kk )then
                        sumex = 0.0d0
                        do n = 1, nnuc
                           sumdx = 0.0d0
                           sumx  = 0.0d0
                           do k = kbeg(m)+1,klast
c..   sum of contributions of each zone in CZ
                              sumdx = sumdx + (x(n,k)-x0(n,k))
c..   use nucleon fraction?
cccccccccccccccccccccccccc

c..   sums to unity over CZ
                              sumx  = sumx  + dtc(k)/dth(2)
c..   dtc(k)/dth(2) = dmh(k)/( xm(kend) - xm(kbeg) )
                           enddo
c                           write(*,*)'SUMX ',m,n,dtc(kbeg(n))/dth(2),
c     1                          sumx,sumx-1.0d0,x(n,kbeg(m))
ccccccccccc

c..   update X over CZ, keep positive
                           do k  = kbeg(m)+1,klast
                              if( x0(n,k) .gt. 0.0d0 )then
c..   positive x0
                                 if( x0(n,k) + sumdx .gt. 0.0d0 )then
c..   x(n,k) > 0.5 x0(n,k)
                                    x(n,k) = x0(n,k) + sumdx
                                 else
c..   avoid x(n,k) < 0
                                    x(n,k) = 0.5d0 * x0(n,k)
                                 endif
                              elseif( abs(x0(n,k)) .lt. 1.0d-10 )then
c..   negative x0 but small absolute value; reset to 0
                                 x(n,k) = 0.0d0
                              else
c..   error
                                 write(*,'(2i5,a5,1p8e12.3)')k,n,
     1                                cnuc(n),x0(n,k)
                                 stop'BURN: error negative x0'
                              endif
                           enddo
c..   check positive abundance
                           do k = kbeg(m)+1,klast
                              if( x(n,k) .lt. 0.0d0 )then
                                 write(*,'(2i5,a5,1p8e12.3)')k,n,
     1                                cnuc(n),x(n,k)
                                 stop'BURN: X negative'
                              endif
                           enddo
                           sumex = sumex + sumdx
cccccccccccccc
                        enddo
                     endif
                  enddo
               endif
            endif

c..   renormalize over zones
c..   nuclear energy/c**2 causes masses to wander slightly
c..   (modify to use baryon number)
            do k = 2, kk
               sumex = 0.0d0
               do n = 1, nnuc
                  sumex = sumex + x(n,k)
               enddo
c               if( abs( sumex-1.0d0 ) .gt. 1.0d-14 )
c     1              write(*,'(a5,i5,1p8e12.3)')'burn',k,sumex-1.0d0
               do n = 1, nnuc
                  x(n,k) = x(n,k)/sumex
               enddo
            enddo

c............................................................
            if( mbflag .eq. 0 )then
c..   only compute finite derivatives for most active zones 
c..   2Msol H burning converges with no derivatives!
               do k = kin, ktop

c..most vigorous burning and avoid no burning (preMSequence)
c..and avoid for hydrodynamics (si burning)
                  if( abs(s5(k)) .gt. 1.0d-2*abs(s5max) .and.
     2                abs(s5max)*dmh(k) .gt. 1.0d-3*abs(tl(1,kk)) 
     3                 )then

                     do m = 1, 2
                        if( m .eq. 1 )then
                           temp(k) = tem(k)  * cfak
                           den(k)  = rhoz(k)
                        else
                           temp(k) = tem(k) 
                           den(k)  = rhoz(k) * cfak
                        endif
c..   initial abundances for solven
c..   convert from mass fraction to moles
                        aex(ndim,k) = x0(ndim,k)
                        do n = 1, ndim-1
                           aex(n,k) = x0(n,k)/xa(n)
                        enddo
c..   leq .eq. 1 avoids overwrite of abundances y(n)-->aex(n,k)

                        call solven(dtc(k),dtn,k,
     1                       idebug,ncytest,ncymax,inegs,it,kk,1,jnb)

c..   construct derivitives (dT, dV) which include nuclear neutrino loss
                        if( m .eq. 1 )then
                           s51(k) = aeps(k) 
                        else
                           s52(k) = aeps(k)
                        endif
                     enddo

                     s5t(k) = (s51(k) - s5(k))/(cfak-1.0d0)/tem(k)
                     s5v(k) = (s5(k)  - s52(k))/(cfak-1.0d0)*rhoz(k)

c..   save derivatives as d ln eps / d ln T, and d ln eps / d ln V
                     if( s5(k) .ne. 0.0d0 )then
                        s5expt(k) = s5t(k)*tem(k)/s5(k)
                        s5expv(k) = s5v(k)/rhoz(k)/s5(k)
                     else
                        s5expt(k) = 0.0d0
                        s5expv(k) = 0.0d0
                     endif
                     sa(k) = s5expt(k)
                     sb(k) = s5expv(k)
                  else
                     sa(k) = 0.0d0
                     sb(k) = 0.0d0
                  endif
               enddo
            endif
         endif
c.................................................................
         if( neutro .eq. 0 )then
c..   neutrino cooling by e+e- processes
c..   finite differences in T and V for ST, SV derivatives
            do k = kin, ktop
               aye(k) = aex(ndim,k)
               if( tem(k) .gt. 1.0d7 )then
                  s1 = 0.0d0
                  s2 = 0.0d0
                  do m = 1,3
                     t9  = tem(k)*1.0d-9
                     rho = rhoz(k)
                     if( m.eq.1 ) t9 = t9 *cfak
                     if( m.eq.2 ) rho= rho*cfak
c..   replace with newer subroutine
                     call sneut(k,t9,rho,aye,aenu)
               
                     if( m .eq. 1 )then
                        s1 = aenu(k)
                     elseif( m .eq. 2 )then
                        s2 = aenu(k)
                     else
                        s4(k) = aenu(k)
                     endif
                  enddo
                  fact   =  cfak - 1.0d0
                  s4t(k) = (s1      - s4(k))/(fact*t(nc,k))
                  s4v(k) = (s4(k) - s2)*rho/fact
               else
c..   below threshold temperature
                  s4(k)  = 0.0d0
                  s4t(k) = 0.0d0
                  s4v(k) = 0.0d0
               endif
c..   neutrino energy change rate and logarithmic derivatives
               snu(k) = s4(k)
               if( s4(k) .ne. 0.0d0 )then
                  snua(k) = s4t(k)*tem(k)/s4(k)
                  snub(k) = s4v(k)*v(1,k)/s4(k)
               else
                  snua(k) = 0.0d0
                  snub(k) = 0.0d0
               endif
            enddo
         else
c..   no emission
            s4(k)   = 0.0d0
            s4t(k)  = 0.0d0
            s4v(k)  = 0.0d0
            snua(k) = 0.0d0
            snub(k) = 0.0d0
         endif
c.....end of non-nuclear neutrino emission......................
      endif
c---------------------------------------------------------------

c..define energy generation and neutrino loss rates, and their
c..combined derivatives for hydro.f
      do k = kin, ktop
         s(4,k) = s4(k)
         s(5,k) = s5(k)
         st(k)  = s4t(k) + s5t(k)
         sv(k)  = s4v(k) + s5v(k)
      enddo

c..   most active heating by nuclear burning
      s5max  = 0.0d0
      ks5max = 0
c..   find most rapid depletion zone for each nucleus
      kxnucl = 0
      dxnuclmin = 1.0d30
c..   greater than Hubble time in seconds
      do n = 1, ndim
c..   for each nucleus n, find the zone k with the most depletion
c..   nuck(n) is zone
c..   xdep(n) is mass fraction change
         xdep(n) = 0.0d0
         nuck(n) = 0
c..   for each nucleus n, find the zone k with the most production
c..   nucpk(n) is zone
c..   xprod(n) is mass fraction change
c..   uses mass fractions
         xprod(n) = 0.0d0
         nucpk(n) = 0
c..time derivative
         do k = kin, ktop
            xd(n,k) = 0
         enddo
      enddo
c..   estimate nuclear time step for each nucleus
      dxnucl = 1.0d40
      nxnucl = 0
c..   find most rapidly depeleted nucleus (nucdep) and the
c..   extent of change (xxdep)
      xxdep = 0.0d0
      nucdep = 0
c..   nucleus with maxinum cycles, corresponding zone k
      ncycmax = 0
      kcycmax = 0

      if( nburn .ne. -1 )then

c..   most active heating by nuclear burning
         do k = kin, ktop
            if( s5max .lt. s(5,k) )then
               s5max = s(5,k)
               ks5max = k
            endif
         enddo

c..   most subcycling in solven
         ncycmx = 0
         kncyc  = 0
         do k = kin, ktop
            if( ncyc(k) .gt. ncycmx )then
               ncycmx = ncyc(k)
               kncyc = k
            endif
         enddo
         if( ncycmx .gt. 10 )then
            write(*,'(a20,4(a12,i5),2(a12,1pe12.3))')
     1           'BURN: max leqs','kncyc',kncyc,'ncycmx',ncycmx,
     2           'ncytot',ncytot,'ks5max',ks5max,
     3           'T',tem(kncyc),'rho',den(kncyc)
         endif


c..   find most rapid depletion zone for each nucleus
c..   greater than Hubble time in seconds
         do n = 1, ndim
c..   for each nucleus n, find the zone k with the most depletion
c..   nuck(n) is zone
c..   xdep(n) is nucleon fraction change
c..   uses mass fractions/atomic mass
            xdep(n) = 0.0d0
            nuck(n) = 0
            do k = kin, ktop
               if( xdep(n) .gt. (x(n,k) - x0(n,k)) )then
                  xdep(n) = (x(n,k) - x0(n,k))
                  nuck(n) = k
               endif
            enddo

c..   for each nucleus n, find the zone k with the most production
c..   nucpk(n) is zone
c..   xprod(n) is nucleon fraction change
c..   uses mass fractions/A
            xprod(n) = 0.0d0
            nucpk(n) = 0
            do k = kin, ktop
               if( xprod(n) .lt. (x(n,k) - x0(n,k)) )then
                  xprod(n) = (x(n,k) - x0(n,k)) 
                  nucpk(n) = k
               endif
            enddo

c..   estimate nuclear time step for each nucleus
            dxnucl = 1.0d40
            nxnucl = 0
            if( xdep(n) .ne. 0.0d0 .and. nuck(n) .gt. 0 )then
ccc   dxnucl = dtc(k)* xa(n)*x(n,nuck(n))/abs( xdep(n) )
               dxnucl = dth(2)*x(n,nuck(n)) / abs( xdep(n) )
            endif
            if( dxnucl .lt. dxnuclmin .and. 
     1           x(n,nuck(n)) .gt. 1.0d-20 )then
               dxnuclmin = dxnucl
               nxnucl    = n
            endif
         enddo
c..   find most rapidly depeleted nucleus (nucdep) and the
c..   extent of change (xxdep)
         xxdep = 0.0d0
         nucdep = 0
         do n = 1, ndim
            if( xxdep .gt. xdep(n) )then
               nucdep = n
               xxdep = xdep(n)
            endif
         enddo
         kxnucl = nuck(nxnucl)
      
         if( nucdep .gt. 0 .and. nucdep .le. ndim 
     1        .and. kxnucl .gt. 1 .and. kxnucl .le. kdm )then
            write(*,'(a14,1pe12.3,a5,i5,a5,3(a10,i5))')'BURN: X change',
     1           xxdep,'k',kxnucl,cnuc(nucdep),'ic(k)',ic(kxnucl),
     2           'ic(k-1)',ic(kxnucl-1),'it',it
         endif

c..   update abundance time derivative
         do k = kin, ktop
            do n = 1,ndim
c..   xd() is used for time step control
c..   this would give only nuclear burning contribution to change
c..   uses mass fractions/A
               xd(n,k) = (x(n,k) - x0(n,k))/dth(2)
ccc   xd(n,k) = (x(n,k) - x0(n,k))/dtc(k)
c..   total change would use xold(n,k) instead of x0(n,k) here
            enddo
         enddo

c..   envelope is not forced to be homogeneous in composition 
c..   with outer zone
c..   envelope abundances (kk+1) are left alone here (see cmix.f)
ccccccccccccccccccccccccccccccccccccccccc

         ncycmax = 0
         do k = 2, kk
            if( ncyc(k) .gt. ncycmax )then
               kcycmax = k
               ncycmax = ncyc(k)
            endif
         enddo

      endif

c      write(*,*)'LEAVING BURN ',x(nnuc-1,2),x(nnuc,2),x(nnuc+1,2),
c     1     x(nnuc-1,535),x(nnuc,535)

c..timing burn, sum of elapsed time in runburn
         call seconds(runburn1)
         runburn = runburn + runburn1-runburn0


      return
      end

