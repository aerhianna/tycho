      subroutine advec(n,ks)

      implicit none

c..   mass loss uses kudritski formulation Teff>7500
c..   else it uses 
c..   mloss = -1 no RSG mass loss
c..   mloss = -2 reimers
c..   mloss = -3 bloecker
c..	  mloss = -4 schroeder (added this line 3/25)
c
c..	  amanda additions 3/26/15
c..   mloss = -5 vink, need to change temperature definitions/boundaries,
c..   since this is all for very high temperature cases
c..   mloss = -6 mokiem

      include 'dimenfile'
      include 'cburn'
      include 'comod'
      include 'compu'
      include 'cgtr.h'
      include 'cconst'
      include 'cadvec'
      include 'ceoset'
      include 'cenv'

      integer*4 n, ks, ncall, namix, k, kwr, j, m, kintenv
      parameter(kintenv = 8000)

      real*8 adtl,astl,flux,tlum,dmasdt,yhe,yh,dmasdy
      real*8 xlogl,tmass,te,yps,xmdcm,zmet,totmass

      real*8 area, eddflux, aog, dkappadr, dmwr, anu
      real*8 wrr(kintenv), wrt(kintenv), wrl(kintenv), wrdmi(kintenv),
     1       wrxm(kintenv), wrsound(kintenv), wrak(kintenv),
     1       wrx(nnuc,kintenv), wrye(kintenv)

      integer*4 nagnuc
      parameter( nagnuc = 286 )
      integer*4 i,  nagz(nagnuc),nagn(nagnuc),naga(nagnuc)
      character*2 cag(nagnuc)
      real*8    xag(nagnuc)
      common /cag/ xag,nagz,nagn,naga,cag
      real*8 xin(ndim),oldx(ndim),dmix,dacc,sum,xave
      save xin, ncall

      data ncall/0/
c--------------------------------------------------------------

      if( mode  .eq. 0 )then
c..no mass loss with explicit hydro
         peryear = 0.0d0
         return
      endif
	write(*,*)'mloss', mloss
      if( mloss .lt. 0 )then
c..   MASS LOST FROM STAR----------------------------
         adtl = abs( tl(1,ks) - tl(1,ks-1) )
         astl = abs( tl(1,ks) )  + abs( tl(1,ks-1) )
         if( adtl .gt. 1.0d-1*astl )then
c..   model not settled
            write(*,'(1p2e11.3,a40)')    adtl,astl,
     1           ' = adtl,astl, mass loss SUPPRESSED'
            peryear = 0.0d0
            return
         endif

c..envelope integration is assumed, so surface at r(kk+1)
         flux  = tl(1,kk+1)/(4.0d0*pi*r(1,kk+1)**2)
         te    = ( flux/sigma )**0.25d0
         tlum  = tl(1,kk+1)/sollum
         totmass = xm(kk+1)/sol
c..   mole fractions 3-6-2008
            yhe    = x(nnuc  ,kk)/xa(nnuc)
            yh     = x(nnuc-1,kk)/xa(nnuc-1)
            yps    = yhe/( yh + 1.0d-5 )
c..   relative to "solar"
            zmet   = ( 1.0d0 - yh - 4.0d0*yhe )/0.014d0
c..   reduction of coupling of pure H, He is about 1e-4 of solar ??
            zmet   = dmax1( zmet, 1.0d-4 )

         if( te .le. 7.500d3 )then
            if( mloss .eq. -1 )then
c..   no RSG mass loss
               peryear = 0.0d0
               vinf    = 0.0d0

            elseif( mloss .eq. -2 )then
               call reimers(te,tlum,r(1,kk+1),totmass,peryear,vinf)

	    elseif( mloss .eq. -3 )then
	       call bloecker(te,tlum,r(1,kk+1),totmass,peryear,vinf,
     1              zmet)

	    elseif( mloss .eq. -4 )then
	       call schroeder(te,tlum,r(1,kk+1),totmass,peryear,vinf)

	    elseif( mloss .eq. -5 )then
c.. Te < 7500 so Vink not used, just use Bloecker instead:
		   call bloecker(te,tlum,r(1,kk+1),totmass,peryear,vinf,
     1              zmet)
	    endif

c..   calculate velocity of Lagrangian grid
c..   d(mass)/dt = accretion rate (solar/year) * g/solar * year/sec
            dmasdt = peryear * sol /secpy

         else
c..BSG mass loss (Kudritski routines)
            xlogl  = log10( tlum )
c..   mole fractions 3-6-2008
            yhe    = x(nnuc  ,kk)/xa(nnuc)
            yh     = x(nnuc-1,kk)/xa(nnuc-1)
            yps    = yhe/( yh + 1.0d-5 )
c..   relative to "solar"
            zmet   = ( 1.0d0 - yh - 4.0d0*yhe )/0.014d0
c..   reduction of coupling of pure H, He is about 1e-4 of solar ??
            zmet   = dmax1( zmet, 1.0d-4 )
            tmass  = xm(kk)/sol

c            if(te .le. 4.0d4)then
               call radwind(xlogl,tmass,te,yps,xmdcm,vinf,zmet)
               dmasdy = -xmdcm
               peryear = dmasdy * secpy/sol
c            else
c               call extowind(te,tlum,xmdcm,vinf,zmet,r(1,kk+1),totmass)
c               peryear = -xmdcm
c               dmasdy = peryear*sol/secpy
c            endif

            if( te .ge. 12.5d3 )then
               if( mloss .eq. -5 )then
                 call vink(te,tlum,r(1,kk+1),totmass,peryear,vinf,zmet)

c              elseif( mloss .eq. -6 )then
c			     call mokiem(te,tlum,r(1,kk+1),totmass,peryear,vinf)

			   endif
            endif

         endif

c     ------------------
c..   AMANDA 10/03/16-5/18/17
c..   Blackman & Owen magnetic field mass loss
c     ------------------
         if( bomloss .eq. 1 )then
            if( peryear .gt. boperyear )then
			   peryear = boperyear
               boflag = 1
			   print *,'USING BLACKOWEN MASS LOSS'
            else
               boflag = 0
            endif
         endif

         write(*,'(2(a20,1pe12.3))')'dM/dt(sol/yr)',peryear,
     1        'v(inf)',vinf
c         if((altloss .eq. 2) .or. (altloss .eq. 3))then
c            call pulses2(n)
c            peryear = peryear + mdotpuls*secpy/sol
c         endif
c         write(*,*)'mdot2',peryear,mdotpuls*secpy/sol
cc         return

         wrflag = 0
         write(*,*)'wrflag=',wrflag

         kwr = kk + jmaxz 
         do m=2, kk
            wrr(m) = r(2,m)
            wrt(m) = t(1,m)
            wrdmi(m) = dmi(m)
            wrxm(m) = xm(m)
            wrak(m) = ak(m)
            wrl(m) = tl(1,m)
            wrsound(m) = sound(m)
            do j=1, nnuc
              wrx(j,m) = x(j,m)
            enddo
            wrye(m) = ye(m)
         enddo

         do m=1,jmaxz-1
            wrr(m+kk) = zr(jmaxz+1-m)
            wrt(m+kk) = ztem(jmaxz+1-m)
            wrdmi(m+kk) = zm(jmaxz-m)-zm(jmaxz+1-m)
            wrxm(m+kk) = wrxm(kk+m-1)+wrdmi(m+kk)
c            write(*,*)m+kk,wrxm(m+kk),wrdmi(kk+m),zm(jmaxz+1-m)
            wrak(m+kk) = zak(jmaxz+1-m)
            wrl(m+kk) = zl(jmaxz+1-m)
            wrsound(m+kk) = zsound(jmaxz+1-m)
            do j=1, nnuc
              wrx(j,m+kk) = x(j,kk)
            enddo
            wrye(m) = ye(kk)
         enddo
         
         if(wrflag .eq. 0)then
           ksonic = 0
           do k = 2, kwr-1
c..   Checks for conditions for Lamers & Nugis 02 optically
c..   thick radiation driven winds if condition not satisfied 
c..   in envel.f.

c..   Check for arad/g > 0.9
             area = pi4*wrr(k)**2
             eddflux = arad*crad/3.0d0 * (wrt(k)**4-wrt(k+1)**4)
     1           * area/wrdmi(k)
             aog = area*eddflux/(wrxm(k)*pi4*grav*crad) 
c             write(*,'(a20,i5,1p8e12.3)')'WR interior ',k,
c     1            aog,eddflux,xm(k),t(2,k),dmi(k),area         

c..   Check for dkappa/dr > 0
             dkappadr = (wrak(k+1) - wrak(k))/(wrr(k+1)-wrr(k))

             if(aog .gt. 1.0d0)then
c               if(dkappadr .gt. 0.0d0)then
c                 if(wrxm(k) .gt. 0.999d0*wrxm(kwr))then
                 if(wrt(k) .lt. 2.0d6)then
                   wrflag = 1
                   ksonic = k
                   rsonic = wrr(ksonic)
                   tsonic = wrt(ksonic)
                   csonic = wrsound(ksonic)
                   msonic = wrxm(ksonic)
                   yesonic = ye(ksonic)
                   kappsonic = wrak(ksonic)
c                   write(*,'(a25,2i5,1p8e12.3)')
c     1                  'WR mass loss in interior ', k,ks, tsonic, 
c     1                 rsonic/solrad, msonic/sol, aog, dkappadr,
c     1                 eddflux,t(2,k+1),dmi(k)
                 else
                   ksonic = k
                 endif
c               else
c                 ksonic = k
c               endif
c               else
c                 ksonic = k
c               endif
             else
               ksonic = k
             endif
           enddo
         endif

c..  Lamers & Nugis 02 mass loss rate for radiation driven winds
c..  in WR stars.
         anu = 0.0d0
         if( (altloss .eq. 1) .or. (altloss .eq. 3))then
           if(wrflag .eq. 1)then
             do i=1,nnuc
                anu = anu + wrx(i,ksonic)
             enddo
             dmwr = -1.267d-13*(8.314d7*(yesonic+(1.0d0/anu))
     1            )**0.5d0*tsonic**4.5d0* 
     1               (rsonic/solrad)**3.0d0*kappsonic/0.75
     1               / (grav*msonic/sol)/sol*secpy 
             
             write(*,'(4(a12,1pe12.3),a12,i5)')
     1            'dmwr',dmwr,'rsonic/solr',rsonic/solrad,'Ts',tsonic,
     2            'Ms/sol',msonic/sol,'ksonic',ksonic
             peryear = peryear + 1.0d0*dmwr
c             if(tsonic .gt. 2.0d5)then
c                peryear = max(peryear,-1.0d-2)
c             endif
             dmasdt = peryear * sol /secpy
             vinf = 2.0d0 * (2.0d0*grav*msonic/rsonic)**0.5d0
             write(*,'(a20,i5,1p12e12.3)')
     1            'WR mass loss on ',ksonic,peryear,vinf, 
     1               tsonic, rsonic/solrad, csonic, msonic/sol,
     2           peryear*sol/secpy*grav*wrxm(ksonic)/rsonic,wrl(ksonic),
     3            wrl(ksonic)*rsonic/(grav*wrxm(ksonic))*secpy/sol
c            stop'ADVEC'
cccccccccccccccccccc

           endif
         endif
         wrflag = 0

c         dmasdt = -1.0d-2 * sol/secpy
c         peryear = dmasdt*secpy/sol
         peryear = peryear*3.0d-1
         
         return

c..   mass changes done in envelope (see hydro.f)

      elseif( mloss .gt. 0 )then
c..   ACCRETION OF MASS ONTO STAR-------------------------------
         if( ncall .eq. 0 )then
c..   set up accreted abundances
            ncall = 1
            do i = 1, ndim
              xin(i) = 0.0d0
            enddo
            do i = 1, netsize
c..   mass fractions
               xin(i) = solarx(i)
            enddo
c..   burn deuterium
            xin(lhe3) = xin(lhe3) + xin(ldeut)
            xin(ldeut) = 0.0d0
            write(*,'(a20,2(a5,1pe12.3))')'D burned on infall',
     1           cnuc(ldeut),xin(ldeut),cnuc(lhe3),xin(lhe3)

            xin(ndim) = 0.0d0
            do i = 1, ndim-1
              xin(ndim) = xin(ndim) + dble( lz(i) )*xin(i)/xa(i)
            enddo
            write(*,*)'accretion of no-D solar abundances, Ye=',
     1           xin(ndim)
         endif
c..   executes this for all calls with accretion
c..   rate of mass addition in cgs units
         dmasdt  = peryear * sol /secpy
c..   amount of mass added this time step
         dacc    = dmasdt * dth(2)
         write(*,'(a12,4(a10,1pe12.3))')'advec ',
     1        'dmh(kk+1)',dmh(kk+1),'dacc', dacc,
     2        'H',xold(nnuc-1,kk+1),'He', xold(nnuc,kk+1)

c..   add mass to zone kk+1
         do n = 1, ndim
c..   reset
            x(n,kk+1) = ( dacc*xin(n) + xold(n,kk+1)*dmh(kk+1) )
     1           /(dacc + dmh(kk+1))
         enddo
         dmh(kk+1) = dmh(kk+1) + dacc
         sum = 0.0d0
         do n = 1, ndim-1
            sum = sum + x(n,kk+1)
         enddo
         if( abs(sum-1.0d0) .lt. 1.0d-5 )then
c..   renormalize
            do n = 1, ndim-1
               x(n,kk+1) = x(n,kk+1)/sum
            enddo
            sum = 0.0d0
            do n = 1, ndim-1
               sum = sum + x(n,kk+1)
            enddo
c..   Ye
            x(ndim,kk+1) = 0.0d0
            do n = 1, ndim-1
               x(ndim,kk+1) = x(ndim,kk+1) 
     1              + x(n,kk+1)*dble( lz(n) )/xa(n)
            enddo
         else
            write(*,*)'ADVEC: renormalization error ',sum-1.0d0
            stop'advec'
         endif

         write(*,'(4(a12,1pe12.4))')'sum-1 ',sum - 1.0d0,
     1        ' added ',dacc/dmh(kk+1), 'He(new)',x(nnuc,kk+1),
     2        'He(old)',xold(nnuc,kk+1)

      else

c..   NO ACCRETION OR MASS LOSS
         peryear = 0.0d0
         return

      endif

      end
