      subroutine fitenv(nc)

      implicit none
c..   cleaned 12-24-08, last revised 10-29-03, rechecked join 12-18-05
c---------------------------------------------------------------------
c..   returns a stencil of 5 envelope models at
 
c                           (1) L*(1+epsill),R
c..  (5) L,R*(1-epsirr)     (3) L,R                    (2) L,R*(1+epsirr)
c                           (4) L*(1-epsill),R

c..   The central value L,R has an inner radius which matches
c..   the outer radius of the henyey grid, r(nc,kk)

c..   the size of the stencil span is scaled by the parameter epsilon

c..   L is the surface luminosity tl(1,kk+1) and is adjusted in the 
c..   iteration loop in hydro.f
c..   R is the envelope outer radius r(1,kk+1) and is adjusted here 
c..   to be consistent with the outer grid radius r(nc,kk)
c..   r(nc,kk) is adjusted in the iteration loop in hydro.f

c..   input:
c..   epsilon (stencil scale)
c..   x(i,kk+1) (abundances in envelope)
c..   xmass (mass enclosed by photosphere)
c..   xmenv (mass in envelope)
c..   the inner temperature of the envelope correspond to T(1,kk)
c..   tl(1,kk+1) (luminosity through envelope)
c..   r(1,kk+1) (initial guess at radius of photosphere)
c..   l = step number (for identification of diagnostic output)
c..   alphaml = ratio of mixing length to pressure scale height
c..   defined as parameter in cnabla

c..   output:
c..   envelope variables for each stencil: vtem, vrho, vp, vr,
c..   vm, vl, vnab, vnrad, vnad, vak, vyef, vel (temperature, density,
c..   pressure, radius, mass, luminosity, nabla, nabla(radiative),
c..   nabla(adiabatic), opacity, Ye(ionized), velocity)
c..   nvmax is number of array elements for each stencil

c..   epi is pressure    at kk (at join radius)
c..   tau is temperature at kk (at join radius)
c..   their derivatives for henyey iteration:
c..   dpidl is d(epi)/dL(join) with R(join) constant
c..   dpidl is d(epi)/dL(join) with R(join) constant
c..   dtaudr is d(tau)/dR(join) with L(join) constant
c..   dtaudr is d(tau)/dR(join) with L(join) constant

c..   xtlum is target luminosity used in hstat.f
c..   rarget is target radius at join, used in hstat.f

c---------------------------------------------------------------------
c..   integrates by numerical recipes variable stepsize Runge-Kutta 
c..   Allows consistency with second-order hydrodynamic
c..   finite differencing.
c..   uses (P,T,r,L), this requires an iteration on equation of state
c..   subroutine. Note: P < P(rad) = arad*T**$/3 is nonphysical (require
c..   positive densities)
c---------------------------------------------------------------------
c..   fitting at boundary: KK   = kk
c..   envelope surface at boundary:   KK+1
c..   envelope surface  T,V,P,ak stored in KK+1 = "kk+3/2"
c..   envelope surface  R,L stored in KK+1 = kk+1
c..   inner boundary of envelope for  R,L at KK = kk
c..   inner boundary of envelope for  T,P at KK = kk (etau,epi)
c..   mass of envelope is xm(kk+1) - xm(kk) = dmh(kk+1)

c..   variable:
c.......TVP.......LR.....TVP.....LR.....TVP.....LR
c........*........|.........*.....|.......*.......|
c......kk-1/2..........kk+1/2..........kk+3/2.....
c...............kk.............kk+1............kk+2
c..   code indices:
c........KK.............KK+1...........KK+2.......(centers)
c...............KK.............KK+1............KK+2(boundaries)
c     
c.................|<- envelope -->|<-envelope surface
c...............TLUM............TLUM.............
c...............rjoin...........RADIUS...........
c.......................................Te,Pe
c----------------------------------------------------------------------

      include 'dimenfile'
c..   store arrays from integration in cenv
      include 'cenv'
      include 'comod'
      include 'compu'
      include 'csurface'
      include 'cconst'
      include 'cnabla'
      include 'ceoset'
      include 'cadvec'

c..   loop index loo defined in cenv

      integer*4 modeg,i,itfinal
      integer*4 j, nc

      real*8    xmass,delr,xmenv,teff,accel,tlhayashi
      real*8    sm,sr,stl,sp,st0,sent,sv0,sd,stau ,db,tb
      real*8    radius0,xtlum0

      real*8    ppd,ppave,ttd,ttave,rrd,rrave
      real*8    chid, chiave, eentd, eentave
      real*8    d2lnrjdl,d2lnrjdr,d2lnsjdl,d2lnsjdr

      real*8    akb
      real*8    omb,xtledd,eps, testeps
      real*8    tlf(nenv),routf(nenv),rjoin(nenv),pf(nenv),tf(nenv)
      real*8    chif(nenv), dnabf(nenv), entf(nenv)

      real*8    aal(ndim,ndim),aab(ndim),aascale,pstar,tstar
      real*8    eentkk
      integer*4 lou

      save

      data modeg/0/

c----------------------------------------------------------------------     
c..   time this operation and accumulate it in runenvrl
      write(*,*)'ENTERING FITENV'
c----------------------------------------------------------------------

      call seconds(runenvel0)

c..   keep tolerance for envelope integration at 0.03 (or less) of that
c..   for model convergence
      eps = 3.0d-3 * resid

c..   epsilon is finite difference step size
c..   epsilon = 0.015 was original choice; read by gen
      epsirr = 0.5d0*epsilon
c..   scale L as R**2
      epsill = epsilon

c..   always update envelope models for beginning of run
c..   kk+2 has abundance from last envelope call
c..   kk+1 has abundance from this envelope call
      if( l .ge. 1 )then
         if( abs( tl(nc,kk)-xtlum ) .lt. epsill*abs( tl(nc,kk) )
     1        .and. abs( r(nc,kk+1)-radius ) .lt. 
     2        epsirr*abs( r(nc,kk+1) )             .and.
     3        abs(dmh(kk+1)-xmenv) .lt. 1.0d-2*abs(dmh(kk+1)) 
     4        .and. abs( x(nnuc,kk+1) + x(nnuc-1,kk+1) - x(nnuc,kk+2)
     5        - x(nnuc-1,kk+2) ) .lt. 1.0d-2 )then
c..	    no major changes; use old stencil of envelopes
ccccccccccccccccccccccccccccccccccccccccc
            write(*,*)'USING OLD STENCIL'
ccccccccccccccccccccccccccccccccccccccccc
            return
         else
            if( l .gt. 1 )then
             write(*,*)'FITENV --- OUT OF BOX, calculate new stencil ',l
c            write(*,'(10x,8a14)')'dL','dR','dS','drin','dmenv','dHkk'
c            write(*,'(a10,1p8e14.6)')'FITENV 5',
c    1          tl(nc,kk)/xtlum-1.0d0,
c    1          r(nc,kk+1)/radius-1.0d0,entropy(kk)/entf(3)-1.0d0,
c    2          r(nc,kk)/rjoin(3)-1.0d0,
c    3          dmh(kk+1)/xmenv-1.0d0,
c    4          x(nnuc,kk+1) + x(nnuc-1,kk+1) - x(nnuc,kk+2)
c    5          - x(nnuc-1,kk+2) 
             write(*,'(12a12)')'tl(nc,kk)','xtlum','r(nc,kk+1)',
     1          'radius','r(nc,kk)','rjoin','dmh(kk+1)','xmenv'
             write(*,'(1p12e12.4)')tl(nc,kk),xtlum, 
     1          r(nc,kk+1),radius, r(nc,kk),rjoin(3),
     2          dmh(kk+1),xmenv,
     3          epsirr, epsill,entropy(kk)

            elseif( l .eq. 1 )then
c..	initial model for this run
               write(*,*)'FITENV --- SET UP BOX ',l
               write(*,'(12a12)')'tl(nc,kk)','xtlum','r(nc,kk+1)',
     1          'radius','r(nc,kk)','rjoin','dmh(kk+1)','xmenv',
     2          'epsirr','epsill','entropy(kk)'
               write(*,'(1p12e12.4)')tl(nc,kk),xtlum, 
     1          r(nc,kk+1),radius, r(nc,kk),rjoin(3),
     2            dmh(kk+1),xmenv,
     3            epsirr, epsill,entropy(kk)
            endif
         endif
      endif

c..   set values in zone kk+1                        
c..   opacity required for mass dmh to be photosphere
c..   would be  ak(kk+2)  =   0.7d0 * a(kk)/dmh(kk)

c..   mixing with envelope done in cmix (wda 10-15-01)
c..   set composition in envelope zone and photosphere zone
      do i = 1, ndim
         x(i,kk+2) = x(i,kk+1)
      enddo

c..   mass interior to photosphere
      xmass = xm(kk) + dmh(kk+1)

c..   integrate in to the zone boundary kk
      xmenv = dmh(kk+1)

c..   adjust luminosity if it is wildly wrong
      omb    = arad*t(nc,kk)**4/( 3.0d0 * p(nc,kk) )
      xtledd = pi4*crad*grav*xm(kk)/ak(kk) * omb

c..   luminosity at midpoint of stencil for envelopes
      xtlum = tl(nc,kk)

      if( xtlum .lt. 0.1d0*xtledd )then
         write(*,*)'RESET xtlum ',xtlum,' -> ', 0.1d0*xtledd
         xtlum = 0.1d0*xtledd
      endif

c..   increase outer radius if it is less than rjoin boundary!
      if( r(nc,kk+1) .le. r(nc,kk) )then
         delr =  v(nc,kk) * dmh(kk+1)/( pi4 * r(nc,kk)**2 )
         r(nc,kk+1) = r(nc,kk) + delr
      endif
c..   auter radius at midpoint of stencil
      radius  = r(nc,kk+1)

c..   define join conditions as seen from grid
      g(kk) = grav*xm(kk)/r(nc,kk)**2
      pstar  = p(nc,kk) - g(kk)/a(kk)*dmh(kk)*0.5d0
c..   use values on grid not envel
      tstar = t(nc,kk) *( 1.0d0 -
     1     (t(nc,kk)-t(nc,kk-1))/(t(nc,kk)+t(nc,kk-1))/
     1    ( (p(nc,kk)-p(nc,kk-1))/(p(nc,kk)+p(nc,kk-1)))
     1     * ( g(kk)/a(kk)*dmh(kk)*0.5d0 /p(nc,kk) )
     1     )

c..assumes entropy(kk) = entropy(kk-1/2)
      eentkk = entropy(kk)

c..count loopbacks
      lou = 0
 1001 continue
      lou = lou + 1

c----------------------------------------------------------
c..   compute numerical derivatives and reference atmosphere
c----------------------------------------------------------
c..   nc=2: move stencil to RADIUS and XTLUM values from hydro.f
c..   matching central value of stencil (3) to radius and entropy
c..   at rjoin

c..   save converged values for reference
      radius0 = radius
      xtlum0  = xtlum


      do i = 1, nenv

c..   radius and luminosity of photosphere
         if( i .eq. 1 )then
            xtlum  = xtlum0 * (1.0d0 + epsill)
            radius = radius0
         elseif( i .eq. 2 )then
            xtlum  = xtlum0
            radius = radius0 * ( 1.0d0 + epsirr)
         elseif( i .eq. 3 )then
c..   next guess is previous value
            xtlum  = xtlum0
            radius = radius0
         elseif( i .eq. 4 )then
            xtlum = xtlum0 * (1.0d0 - epsill)
            radius = radius0
         elseif( i .eq. 5 )then
            xtlum = xtlum0
            radius = radius0 * (1.0d0 - epsirr)
         endif


         call envmod(sm,sr,stl,sp,st0,sent,sv0,sd,stau,teff,xmass,xmenv,
     1        accel,db,tb,eps,nc,itfinal,modeg)

         tlf(i)    = xtlum/sollum
         routf(i)  = radius
         rjoin(i)  = sr
         pf(i)     = sp
         tf(i)     = st0
         chif(i)   = grav*(xmass-xmenv)/(pi4*sr**4*sp)
         dnabf(i)  =  znab(jmaxz)
         entf(i)   = sent
         nvmax(i) = jmaxz


         do j = 1, nvmax(i)
            vtem(j,i) = ztem(j)
            vrho(j,i) = zrho(j)
            vp(j,i) = zp(j)
            vr(j,i) = zr(j)
            vm(j,i) = zm(j)
            vl(j,i) = zl(j)
            vnab(j,i) = znab(j)
            vnrad(j,i) = znrad(j)
            vnad(j,i) = znad(j)
            vak(j,i) = zak(j)
            vyef(j,i) = zyef(j)
            vvel(j,i) = zvel(j)
            vsound(j,i) = zsound(j)
            ve(j,i)   = ze(j)
c..compressibility
            vcmp(j,i) = -ztem(j)*zpt(j)/zpv(j)*zrho(j)
            ventr(j,i) = zentropy(j)
         enddo
      enddo

c..   join radius for reference envelope i=3
      rarget = rjoin(3)

c..   construct logarithimic derivatives with respect to L and 
c..   to photospheric Radius
c..   first construct factor: delta ln(variable)/epsilon
c..   pressure at inner boundary of envelope
      ppd   =  pf(1) - pf(4)
      ppave = (pf(1) + pf(4))
      dppdl = ppd/(epsill * ppave)

      ppd   =  pf(2) - pf(5)
      ppave = (pf(2) + pf(5))
      dppdr = ppd/(epsirr * ppave)

c..   temperature at inner boundary of envelope
      ttd   =  tf(1) - tf(4)
      ttave = (tf(1) + tf(4))
      dttdl = ttd/(epsill * ttave)

      ttd   =  tf(2) - tf(5)
      ttave = (tf(2) + tf(5))
      dttdr = ttd/(epsirr * ttave)

c..   radius at inner boundary of envelope (rjoin)
      rrd   =  rjoin(1) - rjoin(4)
      rrave = (rjoin(1) + rjoin(4))
      drrdl = rrd/(epsill * rrave)

      rrd   =  rjoin(2) - rjoin(5)
      rrave = (rjoin(2) + rjoin(5))
      drrdr = rrd/(epsirr * rrave)

c..   entropy at inner boundary
      eentd  = entf(1) - entf(4)
      eentave= entf(1) + entf(4)
      deentdl= eentd/(epsill*eentave)

      eentd  = entf(2) - entf(5)
      eentave= entf(2) + entf(5)
      deentdr= eentd/(epsirr*eentave)

c..   chi factor (Gm/4 pi r**4 P) at inner boundary
      chid   = chif(1) - chif(4)
      chiave = chif(1) + chif(4)
      dccdl  = chid/(epsill*chiave)

      chid   = chif(2) - chif(5)
      chiave = chif(2) + chif(5)
      dccdr  = chid/(epsirr*chiave)


      if( drrdr .eq. 0.0d0 )then
         write(*,*)'singularity in envelope transformation ',drrdr
         stop'fitenv error in drrdr'
      endif

c..   construct partials of epi,etau with respect to L(inner), R(inner)
      epi    = pf(3)
      etau   = tf(3)
      enchi  = chif(3)
      enab   = dnabf(3)
      eent   = entf(3)

      if( pf(3) .gt. 0.0d0 )then
         piln = dlog( pf(3) )
      else
         write(*,*)'FITENV: negative join pressure ',pf(3)
         stop'fitenv'
      endif

      if( tf(3) .gt. 0.0d0 )then
         tauln = dlog( tf(3) )
      else
         write(*,*)'FITENV: negative join temperature ',tf(3)
         stop'fitenv'
      endif

c..   d ln(pi)/dL at constant rjoin
      dlnpidl  = (dppdl - dppdr*drrdl/drrdr) /xtlum
c..   d(pi)/dL
      dpidl    = dlnpidl* epi
c..   d ln(pi)/drjoin at constant L
      dlnpidr  =  dppdr/drrdr /rarget
      dpidr    =  dlnpidr            * epi

c..   d ln(tau)/dL at constant rjoin
      dlntaudl = (dttdl - dttdr*drrdl/drrdr)/xtlum
      dtaudl   = dlntaudl* etau
c..   d ln(tau)/drjoin at constant L
      dlntaudr = dttdr/drrdr /rarget
      dtaudr   =  dlntaudr             * etau

c..   d (enchi)/dL at constant rjoin
      enchidl  = (dccdl - dccdr*drrdl/drrdr)*enchi/xtlum
c..   d (enchi)/drjoin at constant L
      enchidr  =  dccdr/drrdr * enchi / rarget

c..   d (eent)/dL at constant rjoin
      eentdl   =  (deentdl - deentdr*drrdl/drrdr)/xtlum
c..   d (eent)/drjoin at constant L
      eentdr   =  deentdr/drrdr/ rarget

c..   d (enab)/dL at constant rjoin
      enabdl   =  (denabdl - denabdr*drrdl/drrdr)*enab/xtlum
c..   d (enab)/drjoin at constant L
      enabdr   =  denabdr/drrdr * enab / rarget

      g(kk) = grav*xm(kk)/r(nc,kk)**2

c..   initialize arrays
      do i = 1,ndim
         do j = 1,ndim
            aal(i,j) = 0
         enddo
         aab(i) = 0
      enddo

c..   ln rjoin derivatives in Radius and L
      aal(1,1) = drrdr
      aal(1,2) = drrdl

c..   change to entropy formulation
c..   entropy derivatives in Radius and L(with rjoin constant)
      aal(2,2) = eentdl * xtlum

c      write(*,*)'Determinant ',aal(1,1)*aal(2,2)-aal(1,2)*aal(2,1)

c..   Right hand sides
c..   ln of rjoin
      aab(1) = dlog( r(nc,kk)/rarget )
c..   ln of entropy
      aab(2) = dlog( eentkk/eent )


      call leqs(aal,aab,2,ndim)

c..   aab now contains new Radius and L changes (ln)
      aascale = sqrt( aab(1)**2+aab(2)**2)

c..   new guess
      if( lou .lt. 30 )then
         if( sqrt(aab(1)**2+aab(2)**2) .lt. 0.02d0 )then
c..   linear correction in natural log should work
            radius = radius0*exp(aab(1))
            xtlum  =  xtlum0*exp(aab(2))
         else
c..   limit large change to avoid jump outside circle of convergence
            radius = radius0*exp(aab(1)*0.02d0/aascale)
            xtlum  =  xtlum0*exp(aab(2)*0.02d0/aascale)
         endif
      else
c..slow convergence, use bisection
         radius = radius0*exp(0.5d0*aab(1))
         xtlum  =  xtlum0*exp(0.5d0*aab(2))
      endif
      testeps = sqrt( ((r(nc,kk)-rarget)/r(nc,kk))**2.0d0
     1          +((eentkk-eent)/eentkk)**2.0d0 )

      if( lou .le. 1 )write(*,'(2a5,2a14,14a12)')' ','lou','xtlum',
     1  'radius','testeps','eent','rarget','aab(1)','aab(2)',
     2     'dr/r','ds/s','Det','aascale'
      write(*,'(a5,i5,1p2e14.6,1p14e12.4)')'lou a',lou,xtlum,radius,
     1  testeps,eent,rarget,aab(1),aab(2),(r(nc,kk)-rarget)/r(nc,kk),
     2  (eentkk-eent)/eentkk,aal(1,1)*aal(2,2)-aal(1,2)*aal(2,1),
     3  aascale

c       write(*,'(2(a8,1p5e12.3))')'rjoin',rjoin,'entf',entf
c       write(*,'(1p8e12.3)')
c    1  (rjoin(1)-2.0d0*rjoin(3)+rjoin(4))/rjoin(3),
c    2  (rjoin(1)-rjoin(4))/rjoin(3),
c    3  (entf(1)-2.0d0*entf(3)+entf(4))/entf(3),
c    4  (entf(1)-entf(4))/entf(3),aab(1)

c       write(*,'(1p8e12.3)')
c    1  (rjoin(2)-2.0d0*rjoin(3)+rjoin(5))/rjoin(3),
c    2  (rjoin(2)-rjoin(5))/rjoin(3),
c    3  (entf(2)-2.0d0*entf(3)+entf(5))/entf(3),
c    4  (entf(2)-entf(5))/entf(3),aab(2)



c..test for natural log linearity in interpolation on stencil
      d2lnrjdl =  rjoin(1)/rjoin(3)-rjoin(3)/rjoin(4)
      d2lnrjdr =  rjoin(2)/rjoin(3)-rjoin(3)/rjoin(5)
      d2lnsjdl =  entf(1)/entf(3)  -entf(3)/entf(4)
      d2lnsjdr =  entf(2)/entf(3)  -entf(3)/entf(5)

      if( abs(d2lnrjdl) .gt. 0.02d0 .or.
     1    abs(d2lnrjdr) .gt. 0.02d0 .or.
     2    abs(d2lnsjdl) .gt. 0.02d0 .or.
     3    abs(d2lnsjdr) .gt. 0.02d0 )then
        write(*,*)'MSG(fitenv.f): reduce epsilon ?',epsilon
        write(*,'(4(a12,1pe12.4))')'d2lnrjdl',d2lnrjdl,
     1  'd2lnrjdr',d2lnrjdr,'d2lnsjdl',d2lnsjdl,
     2  'dwlnsjdr',d2lnsjdr
      endif

       
      if( lou .ge. 60 )then
         write(*,*)lou,' steps in fitenv (nonconvergence?)'
         stop'fitenv'
      endif

c..   numerical accuracy of converged procedure is 1.0d-8
c..   keep error small compared to stencil stride epsilon
c      if( testeps .lt. resid )goto 1002
      if( testeps .lt. 1.0d-2 )goto 1002
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      go to 1001

 1002 continue



c..   monitor for trouble in finding a join solution
      if( lou .gt. 4 )
     1     write(*,'(a20,i5,1p8e12.4)')'FITENV: (steps) lou =',
     2     lou,testeps,tl(1,kk),tl(2,kk),vl(1,3)
c     write(*,'(1p5e12.4)')rjoin,entf
ccccccc

c..   envelope is converged; R and L have been adjusted to
c..   fit rjoin and entropy

c..   make envelope luminosity consistent with reference integration 
      tl(nc,kk+1) = vl(1,3)
      tl(nc,kk+2) = vl(1,3)
      xtlum       = vl(1,3)

c..   make envelope radius consistent with reference integration
      radius      = vr(1,3)

c..   radius for optical depth = 2/3
      r(nc,kk+1)  = radius
c..   "radius for zero temperature" by extrapolation
      akb = vak(1,3)
      r(nc,kk+2)  = radius + 4.0d0/(3.0d0*db*akb)

      dnab( kk) = vnab( nvmax(3),3)
      dnrad(kk) = vnrad(nvmax(3),3)
      dnad( kk) = vnad( nvmax(3),3)

c..   assumes no semiconvection at join of grid with envelope
c      if( dnrad(kk) .gt. dnad(kk) )then
c         ic(kk) = 1
c      else
c         ic(kk) = 0
c      endif
      if(ic(kk) .ne. 1)ic(kk) = 0

      enab  =  vnab(nvmax(3),3)
      enrad = vnrad(nvmax(3),3)
      enad  =  vnad(nvmax(3),3)

c..   set convective velocity to be consistent with envelope
c..   used in cmix.f
      hp(kk)     = vvel(nvmax(3),3)
c..center about k=kk
      p(nc,kk+1) = 2.0d0*vp(nvmax(3),3)-p(nc,kk)

      t(nc,kk+1) = 2.0d0*vtem(nvmax(3),3)-t(nc,kk)
      v(nc,kk+1) = 2.0d0/vrho(nvmax(3),3)-v(nc,kk)

c------------------------------------------------------------
       write(*,*)'LEAVING FITENV'
c------------------------------------------------------------
c..sum of elapsed time in fitenv.f
      call seconds(runenvel1)
      runenvel = runenvel + runenvel1-runenvel0

      return
      end



      subroutine envmod(sm,sr,stl,sp,st0,sent,sv0,sd,stau,teff,xmass,
     1        xmenv,accel,db,tb,eps,nc,itfinal,modeg)
c..   driver for getsurf.f and envel2.f
      implicit none

      include 'dimenfile'
      include 'cenv'
      include 'comod'
      include 'compu'
      include 'csurface'
      include 'cconst'
      include 'cnabla'

      integer*4 nc,itfinal,modeg

      real*8 sm,sr,stl,sp,st0,sent,sv0,sd,stau
      real*8 teff,xmass,xmenv,accel,db,tb,pb,eps

      real*8 tmass
      integer*4 kkp,nnc,mmodeg,jj
      common/cderivs/tmass,kkp,nnc,mmodeg,jj
c--------------------------------------------------------
c..   trick for data transfer to nr routines
      kkp = kk+1
      nnc = nc
      tmass = xmass
      mmodeg = modeg

c..   structure variables at outer boundary
c..   sm is mass variable inward from photosphere
      sm  = 0.0d0
      sr  = radius
      stl = xtlum

c..   dtry is initial guess at photospheric density (g/cc)
      if( loo .le. 1 )then
c..   first attempt for this time step
         v(nc,kk+2) = 1.0d+8
c..   all subsequent steps use previous converged value for guess
      endif
c..   set temperature to boundary (tau=0) value
      if( xtlum .gt. 0.0d0 )then
         teff = ( xtlum/(pi4*sigma*radius**2) )**0.25d0
      else
         write(*,*)'fitenv: envmod: xtlum ',xtlum
         stop'envmod 1'
      endif
      t(nc,kk+2) = teff
c..   g=gravitational acceleration (cm/s**2)
      g(kk+2) = grav*xmass/radius**2
      a(kk+2) = pi4*radius**2
      r(nc,kk+2) = radius

c..   calculates boundary values (tau=0) for P,kappa,...
      call getsurf(g,kk+2,nc,itfinal,accel)

c..   itfinal is number of iterations used in getsurf
c..   at outer boundary
      db       = 1.0d0/v(nc,kk+2)
      tb       = t(nc,kk+2)
      pb       = p(nc,kk+2)

c..   integrates from outer boundary in to xminner = xmass - xmenv
c..   stores final values in zone (1,kp)
c..   zone center index for envelope is photosphere index - 1

c..   uses numerical recipes runge-kutta adjustable stepsize and
c..   T,P variables

       call envel4(kk+1,nc,db,tb,pb,xmass,xmenv,
     1    stl,sp,sr,st0,sent,sm,sv0,sd,stau,modeg,eps)

       return
       end

