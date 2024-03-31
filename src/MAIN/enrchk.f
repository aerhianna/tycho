c     
c     
c     
      
      subroutine enrchk(no)
      
      implicit none

c     energy check : purely diagnostic
c     
c     uses GR (newt=0) or newtonian (newt=1) version of total energy 
c     called by online.f  tycho.f
c     
c     following quantities are only computed and printed 
c     to check code, they are not needed in the computational 
c     stream:
c     
c     ekin    = kinetic energy
c     pote    = potential energy
c     eint    = internal energy
c     esou    = net energy gained from internal sources.
c     elum    = net energy radiated from surface
c     echeck  = energy tally as shown in formula below
c     echk at t(n+1)
c     enucl   = nuclear energy relative to amu (12C)
c     
      
      
      include 'dimenfile'
      include 'cconst'
      include 'comod'
      include 'compu'
      include 'cgtr.h'
      include 'cburn'
      include 'cenrchk'

      real*8 dmiscl,ubar,uave,pboc2,gtrfak,dmhscl
      real*8 dpvm
      real*8 ddkin

      real*8 denucl, d2nucl, e2nucl
      real*8 tl1,tl2,delum,desou,dwork,fak,fact1,diss(kdm),work

      real*8 enewt
      
      integer*4 k,n,no
c     integer*4 k,l,n,no


c     

      if( l .eq. 0 )then
c..   first call only
         pdvol = 0.0d0
c..   l = 0 is initial call, prior to time step loop in tycho.f
         n = 1
c..   scale energies to ergs/gram; escl is mass (in grams) on grid
         escl = xm(kk) - xm(1)
c..   integrate over stellar structure
         eint = 0.0d0
         ekin = 0.0d0
         pote = 0.0d0
         pvm  = 0.0d0
         diske = 0.0d0
c..   zone boundary quantities
         do k = 1, kk
            dmiscl = dmi(k)/escl
            if( mode .eq. 0 )then
c     for u, n = 3/2 and 1/2
c     extrapolate to get n=5/2
c     average to get n=2
               ubar = u(n,k)
               uave = 0.5d0* u(1,k)**2
            else
c     for u, n = 2 and 1
               ubar = u(n,k)
               uave = 0.5d0*ubar**2
            endif
            ekin = ekin + dmiscl*uave
            pboc2 = 0.5d0*(p(n,k) + p(n,k+1) + q(n,k) + q(n,k+1))/crad2
            if( r(n,k) .gt. 0.0d0 )then
               if( newt .eq. 0 )then
c..   GR version
                  pote = pote + grav*( gm(k) * dmiscl / r(n,k) )
     1                 * (1.0d0 + pi4*pboc2*r(n,k)**3 /gm(k) )
               else
c..   newtonian version
                  pote = pote + grav*( xm(k) * dmiscl / r(n,k) )
               endif
            endif
         enddo
c..   zone center quantities
         do k = 2, kk
            gtrfak = ( gam(k)/grw(k) + gam(k-1)/grw(k-1) )/2.0d0
            dmhscl = dmh(k)/escl * gtrfak
            dpvm   = 3.0d0*( p(n,k) + q(n,k) )*v(n,k)*dmhscl
            pvm    = pvm + dpvm
            if( newt .eq. 0 )then
c..   GR version
               eint   = eint + dmhscl*e(n,k)
            else
c..   newtonian
               eint   = eint + dmh(k)/escl*e(n,k)
            endif
         enddo
         dmgc2  = ( xm(kk) - gm(kk) )/escl
         dmgc2  = dmgc2*crad2
         elum  = 0.0d0
         esou  = 0.0d0
         pdvol = 0.0d0
         work  = 0.0d0
         desou = 0.0d0
         dwork = 0.0d0
         enewt  = eint + ekin - pote
         enewtz = enewt
      else

c..   in time step loop, use updated value (eg., t(n,k))
         n = 2
c..   scale energies to ergs/gram; escl is mass (in grams) on grid
         escl = xm(kk) - xm(1)
c..   integrate over stellar structure
         eint = 0.0d0
         ekin = 0.0d0
         pote = 0.0d0
         pvm  = 0.0d0
         ddkin = 0.0d0
c..   zone boundary quantities
         do k = 1, kk
            dmiscl = dmi(k)/escl
            if( mode .eq. 0 )then
c..   for u, n = 3/2 and 1/2
c..   extrapolate to get n=5/2
c..   average to get n=2
               ubar = u(n,k) + du(k)*dti(1)*efi(k)
               uave = 0.25d0*(ubar**2 + u(1,k)**2)
            else
c..   for u, n = 2 and 1
               ubar = u(n,k)
               uave = 0.5d0*ubar**2
            endif

c..   dissipation required to damp to HSE
            if( mode .eq. 1 .or. mode .eq. 4 )then
               ddkin = dmiscl * u(n,k)*(r(2,k)-r(1,k))/dth(2)*efi(k) 
     1              + ddkin
            else
               ddkin = 0.0d0
            endif

            ekin = ekin + dmiscl*uave
            pboc2 = 0.5d0*(p(n,k) + p(n,k+1) + q(n,k) + q(n,k+1))/crad2
            if( r(n,k) .gt. 0.0d0 )then
               if( newt .eq. 0 )then
c..   GR version
                  pote = pote + grav*( gm(k) * dmiscl / r(n,k) )
     1                 * (1.0d0 + pi4*pboc2*r(n,k)**3 /gm(k) )
               else
c..   newtonian version
                  pote = pote + grav*( xm(k) * dmiscl / r(n,k) )
               endif
            endif
         enddo

c..   zone center quantities
         do k = 2, kk
            gtrfak = ( gam(k)/grw(k) + gam(k-1)/grw(k-1) )/2.0d0
            dmhscl = dmh(k)/escl * gtrfak
            dpvm   = 3.0d0*( p(n,k) + q(n,k) )*v(n,k)*dmhscl
            pvm    = pvm + dpvm
            if( newt .eq. 0 )then
c..   GR version
               eint   = eint + dmhscl*e(n,k)
            else
c..   newtonian
               eint   = eint + dmh(k)/escl*e(n,k)
            endif
         enddo
         dmgc2  = ( xm(kk) - gm(kk) )/escl
         dmgc2  = dmgc2*crad2
         enewt  = eint + ekin - pote
c..   integrate over structure and time, if step converged
         desou = 0.0d0
         dwork = 0.0d0
         if( no .eq. 0 )then
c..   surface work done over this run
c..   subsequent updates
            pdvol = p(2,kk+1)*a(kk)*u(2,kk)*dth(2)/escl + pdvol
c..centering at time n+1/2 is worse (envelope values?)
c..   echk up to t(n) = t(2) for successful iterations
c..   alphae: time centering of energy equation in dynam and hstat
            tl1   = alphae * tl(2, 1) + (1.0d0-alphae) * tl(1, 1) 
            tl2   = alphae * tl(2,kk) + (1.0d0-alphae) * tl(1,kk) 
            delum = (tl2*efi(kk) - tl1*efi(1))*dth(2)/escl
            elum  = elum + delum
            do k = 2, kk
               fak     = dmh(k) / escl
               fact1    = q(2,k) * dv(2,k) * fak / dth(2)
               diss(k) = diss(k) + fact1 * dth(2)
               dwork   = dwork + fact1
               desou   = desou 
     1              + (ss(k)+snu(k))*dth(2)*efi(k)*( dmh(k)/escl )
            enddo
c..   stress = stress + dstres
            esou = esou + desou
            work = work + dwork * dth(2)
            diske = diske + ddkin
         endif
      endif

c..   uses gtr version of total energy
      if( newt .eq. 0 )then
c..   GR version
         echeck = - dmgc2            - esou + elum + pdvol + diske
      else
c..   newtonian
         echeck = ekin + eint - pote - esou + elum + pdvol + diske
      endif

      if( l .eq. 0 )then
c..   initialize and save echkz, which is initial value of echeck
         echkz = echeck
         write(3,'(a12,1pe15.8)')"echeckz is ",echeck
         write(6,'(a12,1pe15.8)')"echeckz is ",echeck
      endif

c..   nuclear energy relative to amu (pure C12)
c..   uses nuclear mass excesses, converts to ergs/gm
      enucl = 0.0d0
      e2nucl= 0.0d0
      do k = 2, kk
         denucl = 0.0d0
         d2nucl = 0.0d0
         do n = 1, nnuc
            denucl = denucl + qex(n)*x(n,k)    * ergspermev * avagadro
            d2nucl = d2nucl + qex(n)*xold(n,k) * ergspermev * avagadro
         enddo
         enucl  = enucl  + denucl * dmh(k) /escl
         e2nucl = e2nucl + d2nucl * dmh(k) /escl
      enddo

      if( l .eq. 0 )then
c..   initialize and save initial value of enucl
         enuclz = enucl
      endif

c      if( esou .ne. 0.0d0 .and. echkz .ne. 0.0d0 )then
c..   monitor energy check and nuclear check
c         write(*,'(i5,5(a10,1pe11.2))')l,'eps-nuc',esou+enucl-enuclz,
c     1        'step eps',enucl-e2nucl+desou,
c     2        'nuc/eps-1', -(enucl-enuclz)/esou-1.0d0,
c     3        'Ech/z-1',echeck/echkz-1.0d0, 'nuc rnd',
c     4        enuclz*1.0d-16
c..   enucl-e2nucl is subject to round-off problems
c      endif



c     SUCCESS
      return
c     
      end








