
      subroutine cmix(mix,kxmx,dthc,nc)
c..   2-3-2008
c..   hydrodynamic convective mixing
c..   uses previous convection structure ( h(k) )
c..   "extra" (PAY) convective mixing is applied here
c..   "exact" or scaled Thoul diffusion (diffroutine.f with minimal 
c..   modification)
c..   crude rotational mixing
c..   (2,k) ---> (1,k) after mixing to avoid convective "flashing"
c..   uses mass fractions

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'ceoset'
      include 'ceos.h'
      include 'cconst'
      include 'cburn'
      include 'crate'
      include 'comdif'
      include 'comfixx'

c     input: kk, modec, modes, modex,
c     x, dmh, dmi, h, a, dth, v, xold, cdeln, tenvelop
c     xa(cburn), omeg,ajay (initial values), nsemi

c     output: mix, kxmx, x, dthc, dth(if necessary to change)
c     omeg,ajay (updated values)

      integer*4 i,j,k,n, mix,kxmx,nmix,kxmxo,negx,negk,nxmx
      integer*4 nc
      integer*4 kmax,kxloop
      integer*4 iraptor

      integer*4 nexcs,kexcs

      integer*4 ktry(kdm)
      real*8 delxtry

      real*8 delxmx,delxmxo,dmp,sum,sumz

      real*8 fact,dfact

      real*8 fak,dthc

      real*8 ttta(kdm),tttb(kdm),tttc(kdm),tttd(kdm),tttu(kdm)
      real*8 delx(kdm),tttbmax

      real*8 hydtest,hydmax,hydamax

c..use of cdeln causes some problems in semiconvective zone
c..a new variable cdelnc is used, locally defined for now
      real*8 cdelnc
      integer*4 khyd, kzmix

      real*8 schw(kdm), scho(kdm)
cccccccccccccccccccccccccccccccccccccccccccccccccc

      data cdelnc/1.0d-2/
c..   nmix is number of times cmix tries to evaluate changes in X
      data nmix/50/

      save
c---------------------------------------------------------------
c     write(*,*)'ENTERING CMIX'
c------------------------------
c..   paranoid check and
c..   initialize Ye for dthoul.f
      do k = 2, kk+1
         sum = 0.0d0
         sumz = 0.0d0
         do n = 1, nnuc
            sum = sum + x(n,k)
            sumz = sumz + x(n,k)*dble(lz(n))/xa(n)
         enddo
         x(nnuc+1,k) = sumz
         if( abs(sum-1.0d0) .gt. 1.0d-10 )then
            write(*,'(a50,2i5,1p8e12.3)')
     1        'CMIX: fixing normalization error on entry',k,ic(k),
     2           sum-1.0d0
         endif
c..   renormalize
            do n = 1, nnuc
               x(n,k) =  x(n,k)/sum
            enddo
c         endif
c..check renormalizations
            sum = 0.0d0
            sumz = 0.0d0
            do n = 1, nnuc
               sum = sum + x(n,k)
               sumz = sumz + x(n,k)*dble(lz(n))/xa(n)
            enddo
            x(nnuc+1,k) = sumz
            if( abs(sum-1.0d0) .gt. 1.0d-15 )then
               write(*,'(a30,2i5,1p8e12.3)')
     1        'renormalize on entry',k,ic(k),sum-1.0d0
         endif
      enddo

c..designed for hydrostatic mode
      if( mode .eq. 0 )return

      if( mixmode .lt. 0 .or. mixmode .gt. 3 )then
         write(*,*)'error: mixmode = ',mixmode
         stop'cmix: mixmode'
      endif

c..initialize for cinit
c..n=1 to use initial, converged values to set convection flags    
       if( l .le. 1 )then
c..   first time through, set cinit arrays
          call cinit(2,kk,1)
       endif

      do k = 2, kk
         schw(k) = dnrad(k) - dnad(k) - doux(k)
         scho(k) = schw(k)
      enddo

c..set boundary zone for recursion loop to avoid nan's
      if( modes .eq. 2 )then
c..envelope zone kk+1 included
         kmax = kk+1
      else
c..no mass at kk+1
         kmax = kk
      endif

c..   paranoid tests to avoid growing errors from alternating network
c..   and transport solutions
      do k = 2, kmax
         do j = 1,ndim
            if( x(j,k) .lt. 0.0d0 )then
               write(*,'(2a5,a6,8a12)')'k','n','nuc','A','X','T'
               write(*,'(2i5,a6,1p8e12.3)')k,j,cnuc(j),xa(j),
     1              x(j,k), t(1,k)
               stop' CMIX test on entry; x'
            endif
         enddo
      enddo
      do k = 2, kmax
         do j = 1,ndim
            if( xold(j,k) .lt. 0.0d0 )then
               write(*,'(2a5,a6,8a12)')'k','n','nuc','A','Y','T'
               write(*,'(2i5,a6,1p8e12.3)')k,j,cnuc(j),xa(j),
     1              xold(j,k), t(1,k)
               stop' CMIX test on entry; xold'
            endif
         enddo
      enddo
      
c..   save abundances before mixing
      do k = 1, kk+1
         do j = 1,ndim
            scrx(j,k) = x(j,k)
         enddo
      enddo
c..   pad boundary with envelope abundances
      do j = 1,ndim
         scrx(j,kk+2) = scrx(j,kk+1)
      enddo
c..   moles of free particles
      do k = 1, kk+1
         ytoto(k) = 0.0d0
         do j = 1,nnuc
            ytoto(k) = ytoto(k) + x(j,k)/xa(j)*(1.0d0 + dble(nz(j) ) )
         enddo
      enddo

c..   PAY implementation of Press(1982) wave mixing
      call pay(nc)

c..   rotational mixing effects
      call rotate

      do k = 2, kmax
         delx(k) = 0.0d0
      enddo
c..warning: dth(2) above is not adjusted for mix
c..this may give errors in rotation algorithm

      mix = 0

c..   initialize timestep for dthoul.f
      dthc = dth(2)

c..   need same coef.s for all loops on nuclei for nucleon conservation
      dmom(1) = 0.0d0
      dif(1)  = 0.0d0

c..   outer boundary condition: no flow through kk+1
c..   so dmom(kk+1) and dif(kk+1) = 0

      if( modes .ne. 2 )then
c..   revise boundary for no envelope
         dif(kk)  = 0.0d0
         dmom(kk) = 0.0d0
      else
c  boundary with fitenv solution
         dmom(kk+1) = 0.0d0
         dif(kk+1)  = 0.0d0
      endif

c..Thoul diffusion (gravitational settling)
      call dthoul(dthc)

      call xcheck(kk+1,x)


c-------------------------------------------------------------------
c..   use rho * velocity * area centered on boundary for conservation
c--------------------------------------------------------------------
c..   setup mass fluxes
      do k = 2, kk
c..h(k) is convective velocity or eddington-sweet current velocity
         
         dmom(k) =a(k)*h(k)*dthc*0.5d0*
     1       (1.0d0/v(1,k)+1.0d0/v(1,k+1))
         
      enddo

c..   the isolated convection zone is assumed to weakly mix because
c..   of its small depth
c..   isolated convective zone
      do k = 2, kk-1
c         if( dmom(k) .le. 0.0d0 .and. dmom(k-1) .gt. 0.0d0 
c     1        .and. dmom(k+1) .gt. 0.0d0 )then
c..   single zone mixes at half the critical (causal) rate
c            dmom(k)   = 0.125d0*dmi(k)
c            dmom(k)   = 0.0d0*dmi(k)
c            dmom(k) = 0.5d0*(dmom(k-1)+dmom(k+1))
c            ic(k-1) = 1
c            ic(k) = 1
c            ic(k+1) = 1
c         endif
         if( dmom(k) .gt. 0.0d0 .and. dmom(k-1) .le. 0.0d0 
     1        .and. dmom(k+1) .le. 0.0d0 )then
c..   single zone mixes at half the critical (causal) rate
            dmom(k)   = 0.125d0*dmi(k)
c            dmom(k)   = 0.0d0*dmi(k)
         endif
      enddo

c..update to initial convection guess; includes possible diffusion
      do n = 1, nnuc+1
         do k = 1, kk+1
            scrx(n,k) = x(n,k)
         enddo
      enddo

c..   advection and diffusion coefficients are now defined
c..   number of zones with excess abundance change due to mixing
      kexcs = 0
c..   return here for scaled mixing.....................................
 2000 continue

c..   exit for maximum try
      if( mix .gt. nmix ) goto 1010
c..   nmix=0 gives no reduction in time step
c..   count number of tries for mixing
      mix = mix + 1

c..more paranoia about nucleon conservation
      do k = 2, kmax
         sum = -1.0d0
         do j = 1, ndim
            if( x(j,k) .lt. -1.0d-15 )then
               write(*,'(2a5,a6,8a12)')'k','n','nuc','A','Y','T','oldY'
               write(*,'(2i5,a6,1p8e12.3)')k,j,cnuc(j),xa(j),
     1              x(j,k),t(1,k),scrx(j,k)
               stop'CMIX x before convection triad'
            endif
c            sum = sum + x(j,k)*xa(j)
            sum = sum + x(j,k)
         enddo

      enddo

      call xcheck(kk+1,x)

      nxmx = 0
      kxmx = 0
c-----------------------------------------------------loop on nuclei
c..   set up for test of maximum change
      kxmx     = 0
      kxmxo    = 0
      delxmx   = 0.0d0
      delxmxo  = 0.0d0
      if( modes .eq. 2 )then
         kxloop = kk+1
      else
         kxloop = kk
      endif
  
c..   outer boundary condition: no flow through kxloop
c..   so dmom(kxloop) and dif(kxloop) = 0
      dmom(kxloop) = 0.0d0
      dif(kxloop)  = 0.0d0

c..   nearest neighbor factors are same for all nuclei
c..   in one dimension geometry
      do k = 2, kxloop
c..use dthc here
c.......................advection...........semiconv
         ttta(k-1) = (- dmom(k)   /dmh(k) - dif(k)   )*dthc/dth(2)
         tttc(k-1) = (- dmom(k-1) /dmh(k) - dif(k-1) )*dthc/dth(2)
         tttb(k-1) = 1.0d0 - ttta(k-1) - tttc(k-1)

      enddo

c..   try to keep roundoff moderate
c..   choose a shorter "effective time" that still gives homogeneity
      tttbmax = 0.0d0
      do k = 2, kxloop
         tttbmax = dmax1( tttbmax, tttb(k-1) )
      enddo

      if( tttbmax .gt. 1.0d9 )then
c..keep diagonal dominance greater than 1.0d0/1.0d8 to minimize roundoff
         fak = 1.0d8/tttbmax
         write(*,'(4(a20,1pe12.3))')'CMIX: turnovers',
     1        tttbmax,'actual timestep',dthc,'correction',fak,
     2        'mix turnovers',fak*tttbmax -1.0d0
         do k = 2, kxloop
            ttta(k-1) = ttta(k-1) * fak
            tttc(k-1) = tttc(k-1) * fak
            tttb(k-1) = 1.0d0 - ttta(k-1) - tttc(k-1)
         enddo
      endif


      do n = 1, ndim
         do k = 2, kxloop
c..   use value updated for microscopic diffusion
            tttd(k-1) = scrx(n,k)
         enddo

c..   solve tridiagonal by Thomas method to get revised values

         call triad(ttta,tttb,tttc,tttd,tttu,kxloop-1)

c..   update x and find most restrictive zone
         if( modec .ne. -1 )then
            negx = 0
            negk = 0
            do k = 2, kxloop
               x(n,k)  = tttu(k-1)
c..   truncate roundoff to zero, else terminate on error
               if( x(n,k) .lt. -1.0d-16 )then
                  negx = n
                  negk = k
               elseif( x(n,k) .lt. 0.0d0 )then
                  x(n,k) = 0.0d0
               endif
c..   scrx are original abundances before mix
c..   delxmx is largest change in abundance
c..   kxmx is zone with that change
c..   nxmx is nucleus with that change
               delx(k) = dabs( x(n,k) - scrx(n,k) )
               if( delxmx .lt. delx(k) ) then
                  delxmx = delx(k)
                  kxmx = k
                  nxmx = n
               endif
            enddo
c..   make sure zone 1 is updated. Make = zone 2
               x(n,1)  = x(n,2)
c..   truncate roundoff to zero, else terminate on error
               if( x(n,1) .lt. -1.0d-16 )then
                  negx = n
                  negk = 1
               elseif( x(n,1) .lt. 0.0d0 )then
                  x(n,1) = 0.0d0
               endif

            if( negx .ne. 0 )goto 999

         else
c..   else skip update of abundances
         endif
      enddo

      call xcheck(kk+1,x)

c..   fix roundoff 
      do k = 2, kk+1
         sum = -1.0d0
         do n = 1, nnuc
            sum = sum + x(n,k)
         enddo
         if( abs(sum) .gt. 1.0d-15 )then
c..   notification if check error is not small
            if( abs(sum) .gt. 1.0d-8 )
     1           write(*,'(a30,2i5,1p8e12.3)')'CMIX: renormalize',
     2           k,ic(k),sum
c..   reset Ye too
            x(nnuc+1,k) = 0.0d0
            do n = 1, nnuc
               x(n,k) = x(n,k)/( 1.0d0 + sum )
               x(nnuc+1,k) = x(nnuc+1,k) + dble( lz(n) )*x(n,k)/xa(n)
            enddo
         endif
c..   paranoid check
         sum = -1.0d0
         do n = 1, nnuc
            sum = sum + x(n,k)
         enddo
         if( abs(sum) .gt. 1.0d-14 ) write(*,'(a30,2i5,1p8e12.3)')
     1        'err after triad',k,ic(k),sum
      enddo

c-------------------------------------------------end loop on nuclei
c..   moles of free particles
c..   used to define zones which need mixing at constant pressure
c..   and therefore temperature adjustment
      do k = 1, kk+1
         ytot(k) = 0.0d0
         do j = 1,nnuc
            ytot(k) = ytot(k) + x(j,k)/xa(j)*(1.0d0 + dble(nz(j) ) )
         enddo
      enddo
      do k = 1, kk+1
         if( abs( ytot(k)-ytoto(k)) .gt. 2.0d0*cdeln )then
            write(*,'(a20,2i5,1p12e12.3)')'CMIX: Ytot change',k,ic(k),
     2           ytot(k),ytoto(k),ytot(k)-ytoto(k)
     3           ,dmom(k)/dmi(k),dif(k),h(k),x(nnuc-1,k)-xold(nnuc-1,k),
     4           xd(nnuc-1,k)*dth(2),cdeln,entropy(k),doux(k)
         endif
      enddo

c..   count number of zones with excessive changes in mass fraction (kexcs)
      kexcs = 0
      do k = 2, kk
         nexcs = 0
c..sum over all nuclei
         do n = 1, nnuc
            if( abs(  x(n,k) - scrx(n,k) ) .ge. cdelnc )then
               nexcs = nexcs + 1
            endif
         enddo
         if( nexcs .gt. 0 )kexcs = kexcs + 1
      enddo

c..   monitor large changes in abundances
c..   scrx are original abundances before mix
c..   delxmx is largest change in abundance
c..   kxmx is zone with that change
c..   nxmx is nucleus with that change
      if( nxmx .gt. 0 .and. nxmx .le. ndim .and. kxmx .gt. 0 .and.
     1     kxmx .le. kdm .and. abs(delxmx) .gt. cdelnc )then
         write(*,'(a14,0pf10.4,a5,i5,a5,3(a8,i3),a6,1pe11.3,
     1        a6,i5,4(a6,1pe10.2))')
     1        'CMIX: X change',delxmx,
     2        'k',kxmx,cnuc(nxmx),'ic(k)',ic(kxmx),
     3        'ic(k-1)',ic(kxmx-1),'mix',mix,'dthc',dthc,'excs',kexcs,
     4        'dif',dif(kxmx),'dmom',dmom(kxmx)/dmh(kxmx),
     5        'dif-',dif(kxmx-1),'dmom-',dmom(kxmx-1)/dmh(kxmx-1)
      endif

c..   adjust the temperature for mixing at constant pressure
c..   for all zones (strong hydrostatic balance)
      do k = 2, kk
         nc   = 2
         do j = 1, 20
            
            call state(k,k,nc)

            fact   = p(2,k)-p(1,k)
c..   note that dV is zero as there is not net mass motion
            dfact  = -fact/pt(k)
            t(2,k) = t(2,k) + dfact

            if( abs(fact) .lt. 1.0d-7*p(1,k) )goto 100
            
         enddo
         write(*,*)'cmix pressure iteration error ',j
         write(*,*)'at zone ',k,' kk= ',kk
         write(*,'(4(a12,1pe12.3))')'T(2,k)',t(2,k),'p(1,k)',p(1,k),
     1           'p(2,k)',p(2,k),'P2-P1',fact
c         stop'cmix P'

 100     continue
      enddo


c..   n=2 to use new mixed values to set convection flags
c..   for all zones

      do k = 2, kk
         schw(k) = dnrad(k) - dnad(k) - doux(k)
      enddo

 
      if( mix .ge. nmix )then
c..   over the budget of cycles allowed
         write(*,*)'cmix: nmix exceeded ',mix,nmix
         stop'cmix nmix exceeded'
      endif

c..   monitor large change in abundance gradient (semiconvective flicker)
c..   reduce advection coefficient if change is too big
c..   loop over all zones
      do k = 2, kk-1
         delx(k) = 0.0d0
         ktry(k) = 0
c..   loop over all nuclei
         do n = 1, ndim
            delxtry = dabs( x(n,k+1) - scrx(n,k+1) 
     1           - x(n,k) + scrx(n,k) )
            if( delxtry .gt. delx(k) )then
               delx(k) = delxtry
            endif
         enddo

         if( delx(k) .gt. cdelnc*0.1d0 )ktry(k) = 1
      enddo
      
c..   test proton ingestion into helium burning zones but
c..   not carbon or later burning zones, which produce free protons
c..   at these levels
c..   set index to hydrogen (nnuc-1)
      n       = nnuc-1
      khyd    = 0
      kzmix   = 0
      hydmax  = 0.0d0
      hydamax = 0.0d0
      do j = 2, kmax-1
c..   avoid major hydrogen burning regions by using .lt. 1.0d-4
         if( scrx(n,j) .gt. 0.0d0 .and. scrx(n,j) .lt. 1.0d-4
     1        .and. x(n,j) .gt. 1.0d-15)then
c..change due to convection
            hydtest = ( x(n,j) - scrx(n,j) )
            if( abs(hydtest) .gt. cdelnc )then
c..test for large region of significant change
               kzmix = kzmix + 1
               ic(j) = 2
            endif
c..save most extreme case
            if( abs(hydtest) .gt. abs(hydmax) )then
               khyd = j
               hydmax = hydtest
               hydamax = x(n,j) - scrx(n,j)
            endif
         endif
      enddo
c..reduce convective timestep if excessive h mix into he
      if( khyd .gt. 0 )then
         if( abs(hydmax) .gt. cdelnc .or. kzmix .gt. 10 )then
             write(*,'(a20,1pe12.3,2(a6,i5),2(a5,1pe12.3),2(a5,i5))')
     1           'CMIX: H into He',hydmax,'k',khyd,'kzmix',kzmix,
     2           'H',x(nnuc-1,khyd),'T',t(1,khyd),
     3           'ic(k)',ic(khyd),'ic(k-1)',ic(khyd-1)
            do j = 1, kmax
               dmom(j) = 0.5d0 * dmom(j)
               dif(j)  = 0.5d0 * dif(j)
            enddo
            dthc    = 0.5d0 * dthc
            goto 2000
         endif
      endif

 1010 continue
c..   success

c..still more paranoia
      do k = 2, kk
         do j = 1, ndim
            if( x(j,k) .lt. 0.0d0 )then
               write(*,'(2a5,a6,8a12)')'k','n','nuc','A','X','T'
               write(*,'(2i5,a6,1p8e12.3)')k,j,cnuc(j),xa(j),
     1              x(j,k),t(1,k)
               stop'cmix x after triad'
            endif
         enddo
      enddo

      call xcheck(kk+1,x)

c..save convective+diffusive change in abundance
      do k = 1, kk
         do n = 1, ndim
            xdcon(n,k) = x(n,k) - xold(n,k)
         enddo
      enddo
      do n = 1,ndim
         xdcon(n,1) = xdcon(n,2)
      enddo
      iraptor = 0
      do k = 2, kk-1
c..avoid join boundary
         if( p(2,k+1) .ge. p(2,k) )then
           iraptor = iraptor + 1
            write(*,'(a25,a5,i5,4(a8,1pe11.3),4(a6,0pf7.3))')
     1           'CMIX: pressure inverted','k',k,
     2           'dP(2,k)',p(2,k+1)-p(2,k),'dP(1,k)',p(1,k+1)-p(1,k),
     3           'dT(2,k)',t(2,k+1)-t(2,k),'dT(1,k)',t(1,k+1)-t(1,k),
     4           'xold',xold(nnuc,k),'x',x(nnuc,k),
     5           'xold+',xold(nnuc,k+1),
     6           'x+',x(nnuc,k+1)
            write(*,'(2i5,4(a10,1pe12.3))')k,ic(k),'dif',dif(k),
     1           'dmom/dm',dmom(k)/dmi(k)
         endif
      enddo

c..   use convective limit on time step
c..   this is propagated through the rest of the solution
      dth(2) = dthc


c..put revised thermodynamic variables in arrays E(1,k), etc.
c      write(*,*)' CMIX: (2,k) ---> (1,k)'
      do k = 1, kk
         e(1,k) = e(2,k)
         p(1,k) = p(2,k)
         t(1,k) = t(2,k)
         v(1,k) = v(2,k)
      enddo

c--------------------------------
c     write(*,*)'LEAVING CMIX'
c--------------------------------

      return

 999  continue
      write(*,*)' CMIX terminates on X<0 error for ',cnuc(negx)
      write(*,'(5(a12,i5))')'bad zone',negk,'of kk=',kk,'nuc index',negx
      write(*,'(5(a10,1pe12.3))')'x(k-1)',x(negx,negk-1),
     1     'x(k)',x(negx,negk),'x(k+1)',x(negx,negk+1),
     1     'xd*dth', dth(2)*xd(negx,negk)
      write(*,'(5(a10,1pe12.3))')'xold(k-1)',xold(negx,negk-1),
     1     'xold(k)',xold(negx,negk),'xold(k+1)',xold(negx,negk+1)
      write(*,'(a5,a34,i5)')cnuc(negx),
     1     ' -- all abundances in bad zone ',negk
      do n = 1, ndim-4, 5
         write(*,'(5(i5,a5,1pe12.3))')(n+j-1,cnuc(n+j-1),x(n+j-1,negk),
     1        j=1,5)
      enddo
      write(*,*)cnuc(negx),' is negative in these zones:'
      write(*,'(2a5,10a11)')'zone','ic','x','xd*dth','xold','T(k)',
     1     'ttta','tttb','tttc','tttd','b-a-c','tttu'
      do k = 2, kmax
         if( x(negx,k) .lt. 0.0d0 )
     1       write(*,'(2i5,1p10e11.3)')k,ic(k),
     2        x(negx,k),xd(negx,k)*dth(2), xold(negx,k),t(2,k),
     3        ttta(k-1),tttb(k-1),tttc(k-1),tttd(k-1),
     4        tttb(k-1)-ttta(k-1)-tttc(k-1),tttu(k-1)
      enddo
      write(*,'(2a5,10a11)')'zone','ic','x','xd*dth','xold','T(k)',
     1     'ttta','tttb','tttc','tttd','b-a-c'
      stop'cmix test of -x'

      end


      subroutine fixx(k)

      implicit none

      include 'dimenfile'

      include 'cconst'
      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'comdif'
      include 'comfixx'

      integer*4 k,n,j,nc
      real*8 dfact,fact
c---------------------------------------------------
      fakmu = dmh(k)*dmh(k+1)/(dmh(k)+dmh(k+1))

      do n = 1, nnuc+1
c..   fakx is mix parameter between zones k and k+1
         x(n,k) = xold(n,k)*(1.0d0 -fakx*fakmu/dmh(k))
     1        +xold(n,k+1)*fakx*fakmu/dmh(k)
         x(n,k+1) = xold(n,k)*fakx*fakmu/dmh(k+1)
     1        +xold(n,k+1)*(1.0d0-fakx*fakmu/dmh(k+1))
      enddo

      nc   = 2
      do j = 1, 20
         call state(k,k,nc)
         
         fact = p(2,k)-p(1,k)
         dfact = -fact/pt(k)
         t(2,k) = t(2,k) + dfact
         if( abs(fact) .lt. 1.0d-9*p(1,k) )goto 110
      enddo
      write(*,*)'cmix pressure iteration error ',j,k,l
      stop'cmix: fixx 1'

 110  continue
      do j = 1, 20
         call state(k+1,k+1,nc)
         
         fact = p(2,k+1)-p(1,k+1)
         dfact = -fact/pt(k+1)
         t(2,k+1) = t(2,k+1) + dfact
         if( abs(fact) .lt. 1.0d-9*p(1,k+1) )goto 111
      enddo
      write(*,*)'cmix pressure iteration error ',j,k+1
      stop'cmix fixx 2'

 111  continue

      ytot(k)   = 0
      ytot(k+1) = 0
      do n = 1,nnuc+1
         ytot(k)   = ytot(k)   + x(n,k)
         ytot(k+1) = ytot(k+1) + x(n,k+1)
      enddo

c..   reset nablas and convection indices
      call cinit(k,k+1,2)

      return
      end
