      subroutine envel2(kp,nc,db,tb,xmass,xmenv,
     1     stl,sp,sr,st,sm,sv,sd,stau,modeg,eps)

c..   written 6-1-06 wda
c..   Runge Kutta integration from photosphere inward
c..   terminates with envelope mass equal to xmenv
c..   numerical recipes adjustable step
c..   input 
c..   density at boundary db
c     temperature         tb
c     radius              sr
c     mass                sm
c     luminosity          stl
c     target mass         xminner
c..   total stellar mass  xmass
c..   envelope mass       xmenv
c..   output 
c..   zvariables, common cenv

      implicit none

      include 'dimenfile'
      include 'cconst'
      include 'cenv'
      include 'comod'
      include 'cnabla'
      include 'ceoset'

c..   index for integration
      integer*4 j,jj
c..   last is flag for starting integration
c..   iphotflg is flag for photosphere
      integer*4 last, iphotflg
c..   kp is tycho zone for photosphere/envelope
c..   nc is tycho time indes
c..   modeg is mode for nabla calculation
      integer*4 kp,nc,modeg
c..   index loo defined in cenv
      real*8    db,tb,xmass,xmenv,stl,sp,sr,st
      real*8    sm,sv,sd,stau

      INTEGER nvar
      parameter( nvar=4 )
      real*8 x1,x2,y(nvar),dydx(nvar),x10
      real*8 xfinal,yfinal
      real*8 yscal(nvar),yscal0(nvar)

      external derivs

      integer*4 jmax, i
      parameter( jmax = 5000 )
      real*8 xx(jmax),yy1(jmax),yy2(jmax)
      real*8 htry,eps,hdid,hnext,yold,xold
      real*8 test1,test2

      character*5 cvar(nvar)
      data cvar/'r','L','T','V'/
c      parameter( pi = 3.141592653589793d0 )
      real*8 tmass
c..   kp is tycho zone for photosphere/envelope
c..   nc is tycho time indes
c..   modeg is mode for nabla calculation
c..   this trick gets values into derivs.f without global change
      integer*4 kkp,nnc,mmodeg
      common/cderivs/tmass,kkp,nnc,mmodeg,jj

      integer*4 nmax,nstpmx,ndum,ndumstep
      real*8 hdum
      PARAMETER (NMAX=50,NSTPMX=200)
      DOUBLE PRECISION xdum(NSTPMX),ydum(NMAX,NSTPMX),sigma2
      COMMON /path/ xdum,ydum,sigma2
c----------------------------------------------------------------
      write(*,'(a30,1pe15.7,i15)')'entering envel2 ',dmh(kk),kk
ccccccccccccccccccccccccccccccccccccc
c..flag for smooth mass change to join
      iphotflg = 0
      
      kkp = kp
      nnc = nc
      mmodeg = modeg
      tmass = xmass
c...........................................
c..   y(1) = r
c..   y(2) = L/Lsun
c..   y(3) = T
c..   y(4) = V
c..start at photosphere and integrate inward
c...........................................
      y(1) = sr
      y(2) = stl/sollum
      y(3) = tb
      if( db .gt. 0.0d0 )then
         y(4) = 1.0d0/db
      else
         write(*,*)'ENVEL2.f: density error db=',db
         stop'envel2.f'
      endif

      x1 = 0.0d0
      call derivs(x1,y,dydx)

c..integrate in to radius xfinal
      xfinal = -xmenv
c..normalized to initial values
      do i = 1, nvar
         yscal(i) = y(i)
         yscal0(i) = y(i)
      enddo

      call getfirst(kp,nc,htry)

      j = 1
      ztem(j) = y(3)
      zrho(j) = 1.0d0/y(4)
      zl(j)   = y(2) * sollum
      zr(j)   = y(1)
      zm(j)   = x1

      do j = 2, jmax
c..save index for z arrays
         jj = j


         call derivs(x1,y,dydx)

c      write(*,'(a20,i5,1p12e11.3)')'after derivs ',j,y,dydx,htry,x1
cccccccccccc

c..avoid very large steps in smooth region
c         if( -htry .gt. 1.0d-2*xmenv )then
c            htry = -1.0d-2*xmenv
c         endif
c         if( -htry .gt. dmh(kk) )then
c            htry = -dmh(kk)
c         endif

         do i = 1, nvar
            yscal(i) = yscal0(i)
c            yscal(i) = y(i)
         enddo

         
         call rkqs(y,dydx,nvar,x1,htry,eps,yscal,hdid,hnext,derivs)

c..   save arrays
         ztem(j) = y(3)
         zrho(j) = 1.0d0/y(4)
         zl(j)   = y(2) * sollum
         zr(j)   = y(1)
         zm(j)   = x1
         zv(j)   = y(4)

c..stop at final mass coordinate
         if( x1 + 4.0d0*hdid -xfinal .gt. 0.0d0 
     1        .and. iphotflg .eq. 0 )then
c..not near boundary
c..use 4 instead of 5 for smoother increase than nr routine
            htry = 4.0d0*hdid
         else
c..   last ndumstep steps for smoother join to henyey mesh
c            write(*,*)j,x1-xfinal
cccccccccccc
c..   rkqs tries to take a step no more than 4*hdid
c..   this attempts to make a smooth join to equal mass zoning
            ndumstep = 4
            x2 = xfinal

            call rkdumb(y,nvar,x1,x2,ndumstep,derivs)

            do i = 1, ndumstep+1
               jj = j + i -1
               zr(jj)   = ydum(1,i)
               zl(jj)   = ydum(2,i)*sollum
               ztem(jj) = ydum(3,i)
               zv(jj)   = ydum(4,i)
               zm(jj)   = xdum(i)
               zrho(jj) = 1.0d0/zv(jj)
               t(nc,kp) = ztem(jj)
               v(nc,kp) = zv(jj)

               call state(kp,kp,nc)

               zp(jj) = p(nc,kp)
               ze(jj) = e(nc,kp)
               zak(jj) = ak(kp)
            enddo

c            do i = jj-ndumstep*2,jj
c               write(*,'(i5,1p8e12.4)')i,zm(i),zr(i),zp(i),
c     1              zm(i)-zm(i-1)
c            enddo

            goto 100
         endif

      enddo
 100  continue
      j = jj
      sr  = zr(j)
      sm  = zm(j)
      stl = zl(j)
      jmaxz = j
      st  = ztem(j)
      sv  = zv(j)
      t(nc,kp) = st
      v(nc,kp) = sv

c..   get consistent pressure at join from actual T and V
c..   used as join pressure

      call state(kp,kp,nc)

c..   save values at end of stepsize
      zv(j)   = v(nc,kp)
      ztem(j) = t(nc,kp)
      zrho(j) = 1.0d0/zv(j)
      zp(j)   = p(nc,kp)
      sp      = p(nc,kp)
      zr(j)   = sr
      zm(j)   = sm
      zl(j)   = stl
      znab(j) = dnab(kp)
      znad(j) = dnad(kp)
      znrad(j)= dnrad(kp)
      zvel(j) = hp(kp)
      zak(j)  = ak(kp)
      zyef(j) = yef(kp)
      zsound(j) = sound(kp)
      ze(j)   = e(nc,kp)

      zv(j)   = sv              
      zev(j)  = ev(kp)          
      zet(j)  = et(kp)          
      zpt(j)  = pt(kp)          
      zpv(j)  = pv(kp)   

      

c      write(*,'(i5,1p8e12.4)')j,sp,sr,sm,stl,sv
c      stop
cccccccccccccccccccccccccccccccccccccccccc

      return
      end




      subroutine getfirst(kp,nc,htry)
      implicit none
c..gets starting value (photospheric) for mass step in envel2.f
      include 'dimenfile'
      include 'comod'
      include 'ceoset'
      include 'cenv'
      include 'cconst'
      include 'cadvec'

      integer*4 kp,nc
      real*8 htry,sak,sr
c--------------------------------------------------------

      call state(kp,kp,nc)

      sak = ak(kp)
      sr  = r(nc,kp)
c..   sign of dmass determines direction of integration
c..   positive outward, negative inward
c..   use fraction of optical mean free path to begin
      htry = -pi4*sr**2/sak * 0.001d0
c..   htry is step size in mass coordinate

      return
      end



      subroutine derivs(xx,vv,dvv)
      implicit none

      include 'dimenfile'
      include 'comod'
      include 'ceoset'
      include 'cenv'
      include 'cconst'
      include 'cadvec'

      integer*4 nmax
      parameter(nmax=4)
      real*8 xx, dvv(NMAX), vv(NMAX)
      real*8 area, gravity, dpdm
      real*8 sp,sak,spt,spv,set,sev,sv,sd,dvdt,dnad,stl,sr
      real*8 dnab,dnrad,st,svel,dtdm,fakv,dvdm,dldm,drdm
      real*8 tmass

      integer*4 kp,nc,modeg,i,j
      common/cderivs/tmass,kp,nc,modeg,j
c--------------------------------------------------------
c..   y(1) = r
c..   y(2) = L
c..   y(3) = T
c..   y(4) = V
c
      sr      = vv(1)
      area    = pi4* sr**2
      gravity = grav*(xx + tmass)/sr**2
      sv      = vv(4)
      st      = vv(3) 
      stl     = vv(2) * sollum
c..set arrays for using state.f and fnab.f
      t(nc,kp)  = st
      v(nc,kp)  = sv
      r(nc,kp)  = sr
      tl(nc,kp) = stl

c..   get pressure, opacity, derivatives

      call state(kp,kp,nc)

c..define nablas at boundary using boundary variables
c..on mesh in TYCHO dnad is averaged across a boundary
      sp   = p(nc,kp)
      sak  = ak(kp)
      spt  = pt(kp)
      spv  = pv(kp)
      set  = et(kp)
      sev  = ev(kp)
      sd   = 1.0d0/sv
      dvdt = - et(kp)/( p(nc,kp) + ev(kp) )
      dnad = p(nc,kp)/( t(nc,kp)*( pt(kp) + pv(kp)*dvdt ))

      dpdm = - gravity/area
      st   = t(nc,kp)
      dtdm   = -stl/area**2*0.75d0*sak/(arad*crad*st**3)
      dnrad  = sp/dpdm*dtdm/st

      if( modeg .eq. 0 )then
c..   boehm-vitense

         call fnab(dnab,dnad,dnrad,
     1        st,sv,sp,gravity,spt,spv,set,sev,sak,svel,0.0d0,
     2        alphaml,modec,0)

         dtdm   = dnab * st/sp * dpdm
         fakv   = (1.0d0 - dnab*st/sp*spt)/spv*sp/sv
         dvdm   = fakv*sv/sp*dpdm


      elseif( modeg .eq. 1 )then
c..   adibatic, real eos
         dnab    = dmin1(dnad,dnrad)
         dtdm   = dnab * st/sp * dpdm
         fakv   = (1.0d0 - dnab*st/sp*spt)/spv*sp/sv
         dvdm   = fakv*sv/sp*dpdm
      else
c..   adiabatic, ideal eos
         dnab    = dmin1(dnad,dnrad)
         fakv   = - 0.75d0
         dvdm   = - 0.75d0*sv/sp*dpdm
         dtdm   =   0.25d0*st/sp*dpdm
      endif

      dldm   = 0
      drdm   = sv/area
c..   define derivatives for return
      dvv(1) = drdm
      dvv(2) = dldm
      dvv(3) = dtdm
      dvv(4) = dvdm

c..   save
      ztem(j) = st
      zrho(j) = sd
      zp(j)   = sp
      zr(j)   = sr
      zm(j)   = xx
      zl(j)   = stl
      znab(j) = dnab
      znad(j) = dnad
      znrad(j)= dnrad
      zvel(j) = svel
      zak(j)  = ak(kp)
      zyef(j) = yef(kp)
      zsound(j) = sound(kp)
      ze(j)   = e(nc,kp)
      zv(j)   = sv  
      zev(j)  = sev  
      zet(j)  = set 
      zpt(j)  = spt  
      zpv(j)  = spv 

      return
      end

