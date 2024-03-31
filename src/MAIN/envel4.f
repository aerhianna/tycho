      subroutine envel4(kp,nc,db,tb,pb,xmass,xmenv,
     1     stl,sp,sr,st,sent,sm,sv,sd,stau,modeg,eps)

c..   written 6-1-06 wda, revised 8-27-07
c..   Runge Kutta integration from photosphere inward
c..   terminates with envelope mass equal to xmenv
c..   numerical recipes adjustable step
c..   added constraints on mass step and temperature step

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
      real*8    db,tb,pb,xmass,xmenv,stl,sp,sr,st,sent
      real*8    sm,sv,sd,stau

      INTEGER nvar
      parameter( nvar=4 )
      real*8 x1,x2,y(nvar),dydx(nvar),x10
      real*8 xfinal,yfinal
      real*8 yscal(nvar),yscal0(nvar)

      external derivs4

      integer*4 jmax, i
      parameter( jmax = 4000 )
      real*8 xx(jmax),yy1(jmax),yy2(jmax)
      real*8 htry,eps,hdid,hnext,yold,xold
      real*8 test1,test2,htrym,htryt,htryv,htry0

      character*5 cvar(nvar)
      data cvar/'r','L','T','P'/

      real*8 tmass
      real*8 f,df
c..   kp is tycho zone for photosphere/envelope
c..   nc is tycho time index
c..   modeg is mode for nabla calculation
c..   this trick gets values into derivs.f without global change
      integer*4 kkp,nnc,mmodeg
      common/cderivs/tmass,kkp,nnc,mmodeg,jj

      integer*4 nmax,nstpmx,ndum,ndumstep
      real*8 hdum
c..   revise Numerical Recipes default dimension
      parameter( nmax=4, nstpmx=50000)
      DOUBLE PRECISION xdum(NSTPMX),ydum(NMAX,NSTPMX),sigma2
      COMMON /path/ xdum,ydum,sigma2
c----------------------------------------------------------------
c..flag for smooth mass change to join
      iphotflg = 0
      kkp = kk+1
      nnc = nc
      mmodeg = modeg
      tmass = xmass
c...........................................
c..   y(1) = r
c..   y(2) = L/Lsun
c..   y(3) = T
c..   y(4) = P
c..start at photosphere (T=Teff) and integrate inward
c...........................................
      y(1) = sr
      y(2) = stl/sollum
      y(3) = tb
      y(4) = pb
      
c..   v(nc,kkp) is used in derivs4 to start integration
c..   db is consistent with pb and tb
      v(nc,kkp) = 1.0d0/db
      jj = 1
      x1 = 0.0d0

      call derivs4(x1,y,dydx)

c..integrate in to radius xfinal
      xfinal = -xmenv
c..normalized to initial values
      do i = 1, nvar
         yscal(i) = y(i)
         yscal0(i) = y(i)
      enddo


      call getfirst4(kp,nc,htry)

      j = 1
      ztem(j) = y(3)
      zp(j)   = y(4)
      zl(j)   = y(2) * sollum
      zr(j)   = y(1)
      zm(j)   = x1
      zrho(j) = db


c.........................................................
      do j = 2, jmax
c..save index for z arrays
         jj = j
         
         call derivs4(x1,y,dydx)

         do i = 1, nvar
            yscal(i) = y(i)
         enddo
         htry0 = -htry
c..   keep changes in mass. T, V moderate
c..   Numerical Recipes stepper is too optimistic
c..   minimum magnitude of change
         htrym = abs(1.0d-2* y(1)/dydx(1))
         htryt = abs(0.1d0* y(3)/dydx(3) )
         htryv = abs(0.1d0* y(4)/dydx(4) )
c..set sign for inward integration
         htry = - dmin1(htrym,htryt,htryv,htry0)


         call rkqs(y,dydx,nvar,x1,htry,eps,yscal,hdid,hnext,derivs4)

c..   save arrays
         ztem(j) = y(3)
         zp(j)   = y(4)
         zl(j)   = y(2) * sollum
         zr(j)   = y(1)
         zm(j)   = x1
         
c..logic for exit of adaptive RK to make smooth join with grid
          if( abs(hdid) .gt. dmh(kk) 
c        if( abs(hdid) .gt. dmh(kk) + dmh(kk-1)
     1        .or.  x1 + hnext -xfinal -4.0d0*dmh(kk) .lt. 0.0d0 )then
c        if(  x1 + hnext -xfinal -50.0d0*dmh(kk) .lt. 0.0d0 )then
            iphotflg = 1

         endif
    

c..stop at final mass coordinate
c         if( x1 + hnext -xfinal -4.0d0*dmh(kk) .gt. 0.0d0 
c..   if stepsize is larger than join zone, use that stepsize
c..   as constant in rkdumb.f
         if( iphotflg .eq. 0 )then
c..not near boundary
            htry = hnext
c..returns to adjustable rk integration
         else

c..   last ndumstep steps for smoother join
            x2 = xfinal
            ndumstep = (x1-x2)/dmh(kk)
c..   avoid excessive ndumstep
            
            if( ndumstep .gt. 200 )ndumstep = 200
            

            call rkdumb(y,nvar,x1,x2,ndumstep,derivs4,ydum)

            do i = 1, ndumstep+1
               jj = j + i -1
               zr(jj)   = ydum(1,i)
               zl(jj)   = ydum(2,i)*sollum
               ztem(jj) = ydum(3,i)
               zp(jj)   = ydum(4,i)
               zm(jj)   = xdum(i)
               v(nc,kp) = zv(jj)
               zrho(jj) = 1.0d0/zv(jj)

            enddo           
            goto 100
         endif

      enddo
 100  continue
      
      jmaxz = jj
      
c..   get consistent specific volumes at join from actual T and P
      do j = 1, jmaxz
         t(nc,kp) = ztem(j)
         v(nc,kp) = 1.0d0/zrho(j)
         p(nc,kp) = zp(j)
c..   get density, opacity, derivatives by iteration
         do i = 1,30

            call state(kp,kp,nc)

            f = p(nc,kp) - zp(j)
            df = pv(kp)
            if( v(nc,kp)-f/df .gt. 0.5d0*v(nc,kp) )then
               v(nc,kp) = v(nc,kp) -f/df
            else
               v(nc,kp) = 0.5d0*v(nc,kp)
               write(*,'(a40,i5,1p8e12.3)')
     1 'ENVEL4.f V iter jump restricted to V/2',i,v(nc,kp)
            endif
            if( abs(f/df) .lt. 1.0d-8*v(nc,kp) )goto 1000
         enddo

         write(*,*)'iteration exceeded ',i,' in testenvel2.f'
         stop'testenvel2.f'
 1000     continue


c..   save values at each step using values from T and P
          zv(j)   = v(nc,kp)
          zrho(j) = 1.0d0/zv(j)
c          znab(j) = dnab(kp)
c          znad(j) = dnad(kp)
c          znrad(j)= dnrad(kp)
c          zvel(j) = hp(kp)
          zak(j)  = ak(kp)
          zyef(j) = yef(kp)
          zsound(j) = sound(kp)
          zentropy(j) = entropy(kp)
          ze(j)   = e(nc,kp)
          zev(j)  = ev(kp)          
          zet(j)  = et(kp)          
          zpt(j)  = pt(kp)          
          zpv(j)  = pv(kp)   

      enddo
c..join values
      sm  = zm(jmaxz)
      sr  = zr(jmaxz)
      stl = zl(jmaxz)
      sp  = zp(jmaxz)
      sv  = zv(jmaxz)
      st  = ztem(jmaxz)
      sent = zentropy(jmaxz)

      return
      end




      subroutine getfirst4(kp,nc,htry)
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



      subroutine derivs4(xx,vv,dvv)
      implicit none

      include 'dimenfile'
      include 'comod'
      include 'ceoset'
      include 'cenv'
      include 'cconst'
      include 'cadvec'
      include 'cnabla'

      integer*4 nmax
      parameter(nmax=4)
      real*8 xx, dvv(NMAX), vv(NMAX)
      real*8 area, gravity, dpdm
      real*8 sp,sak,spt,spv,set,sev,sv,sd,dvdt,dnad1,stl,sr
      real*8 dnab1,dnrad1,st,svel,dtdm,fakv,dvdm,dldm,drdm
      real*8 tmass
      real*8 f,df

      integer*4 kp,nc,modeg,i,j
      common/cderivs/tmass,kp,nc,modeg,j
c--------------------------------------------------------
c..   y(1) = r
c..   y(2) = L
c..   y(3) = T
c..   y(4) = P

c..   svel overwritten in fnab.f if convective
      svel    = 0.0d0
      sr      = vv(1)
      area    = pi4* sr**2
      gravity = grav*(xx + tmass)/sr**2
      sp      = vv(4)
      st      = vv(3) 
      stl     = vv(2) * sollum
c..set arrays for using state.f and fnab.f
      t(nc,kp)  = st
      r(nc,kp)  = sr
      tl(nc,kp) = stl
      
c..   guess last value
c..   set v(nc,kp) consistently with sp and st in envel4.f
      sv        = v(nc,kp)
c..   get density, opacity, derivatives by iteration
      do i = 1,50
     
         call state(kp,kp,nc)

         f = p(nc,kp) - sp
         df = pv(kp)

         if( sp .le. arad*t(nc,kp)**4/3.0d0 )then
c..pressure exceeds rad pressure, wants rho < 0 error
            write(*,'(/a40)')'envel4.f/derivs.f: error P<Prad'
            write(*,*)'This occurs when the stencil of envelopes'
            write(*,*)'intrudes into the unstable region where '
            write(*,*)'mass loss becomes catastrophic'
            write(*,*)'Try reducing epsilon in params.d'

            write(*,'(a12,i12)')'L',l
            write(*,'(a12,i12)')'iteration',i
            write(*,'(a12,i12)')'env zone',j
            write(*,'(a12,i12)')'kk+1',kp
            write(*,'(a12,1pe12.4)')'epsilon',epsilon
            write(*,'(a12,1pe12.4)')'pressure',sp
            write(*,'(a12,1pe12.4)')'rad P',arad*t(nc,kp)**4/3.0d0
            write(*,'(a12,1pe12.4)')'L/Ledd',
     1   tl(nc,kp) * ak(kp)*3.0d0*p(nc,kp)
     2   /(pi4*crad*grav*xm(kp)*arad*t(nc,kp)**4) 
            write(*,'(a12,1pe12.4)')'g',grav*xm(kp)/r(nc,kp)**2
            write(*,'(a12,1pe12.4)')'E',e(nc,kp)
            write(*,'(a12,1pe12.4)')'GM/R',grav*xm(kp)/r(nc,kp)
            write(*,'(a12,1pe12.4)')'R',r(nc,kp)
            write(*,'(a12,1pe12.4)')'rho',1.0d0/v(nc,kp)
            write(*,'(a12,1pe12.4)')'T',t(nc,kp)
            write(*,'(a12,1pe12.4)')'ak',ak(kp)
            write(*,'(a12,1pe12.4)')'xm',xx+tmass
            write(*,'(a12,1pe12.4)')'r',sr
            write(*,'(a12,1pe12.4)')'xx',xx

c            stop'envel4.f: arad err'
         endif

c..limit change in V for convergence with Newton-Raphson
         if( -f/df .gt. v(nc,kp) )then
            v(nc,kp) = v(nc,kp)*2.0d0
         elseif( -f/df .lt. -0.5d0*v(nc,kp) )then
            v(nc,kp) = v(nc,kp)*0.5d0
         else
            v(nc,kp) = v(nc,kp) -f/df
         endif

         if( abs(f/df) .lt. 1.0d-8*v(nc,kp) )goto 100
      enddo
      write(*,*)'iteration exceeded ',i,' in derivs4.f'
      write(*,'(a20,1pe12.4)')'f/df/V',f/df/v(nc,kp)
      write(*,'(a20,1pe12.4)')'f',f
      write(*,'(a20,1pe12.4)')'df',df
      write(*,'(a20,2i5)')'nc and kp',nc,kp
      write(*,'(a20,1pe12.4)')'V',v(nc,kp)
      write(*,'(a20,1pe12.4)')'T',t(nc,kp)
      write(*,'(a20,1pe12.4)')'P',p(nc,kp)
      write(*,'(a20,1pe12.4)')'Prad',arad*t(nc,kp)**4/3.0d0
      write(*,'(a20,1pe12.3)')'Ropal',1.0d0/
     1  ( v(nc,kp)*( t(nc,kp)*1.0d-6 )**3 )
      stop'derivs4.f'
 100  continue

c..define nablas at boundary using boundary variables
c..on mesh in TYCHO dnad is averaged across a boundary
      sp   = p(nc,kp)
      sv   = v(nc,kp)

      sak  = ak(kp)
      spt  = pt(kp)
      spv  = pv(kp)
      set  = et(kp)
      sev  = ev(kp)
      sd   = 1.0d0/sv
      dvdt = - et(kp)/( p(nc,kp) + ev(kp) )
      dnad1 = p(nc,kp)/( t(nc,kp)*( pt(kp) + pv(kp)*dvdt ))

      dpdm = - gravity/area
      st   = t(nc,kp)
      dtdm   = -stl/area**2*0.75d0*sak/(arad*crad*st**3)
      dnrad1  = sp/dpdm*dtdm/st

      if( modeg .eq. 0 )then
c..   boehm-vitense

         call fnab(dnab1,dnad1,dnrad1,
     1        st,sv,sp,gravity,spt,spv,set,sev,sak,svel,0.0d0,
     2        alphaml,uuml,modec,0)

         dtdm   = dnab1 * st/sp * dpdm
         fakv   = (1.0d0 - dnab1*st/sp*spt)/spv*sp/sv
         dvdm   = fakv*sv/sp*dpdm

      elseif( modeg .eq. 1 )then
c..   adibatic, real eos
         dnab1    = dmin1(dnad1,dnrad1)
         dtdm   = dnab1 * st/sp * dpdm
         fakv   = (1.0d0 - dnab1*st/sp*spt)/spv*sp/sv
         dvdm   = fakv*sv/sp*dpdm
      else
c..   adiabatic, ideal eos
         dnab1    = dmin1(dnad1,dnrad1)
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
      dvv(4) = dpdm

c..   save
      ztem(j) = st
      zrho(j) = sd
      zp(j)   = sp
      zr(j)   = sr
      zm(j)   = xx
      zl(j)   = stl
      znab(j) = dnab1
      znad(j) = dnad1
      znrad(j)= dnrad1
      dnab(kp) = dnab1
      dnad(kp) = dnad1
      dnrad(kp) = dnrad1
      zvel(j) = svel
      zak(j)  = ak(kp)
      zyef(j) = yef(kp)
      zsound(j) = sound(kp)
      zentropy(j) = entropy(kp)
      ze(j)   = e(nc,kp)
      zv(j)   = sv  
      zev(j)  = sev  
      zet(j)  = set 
      zpt(j)  = spt  
      zpv(j)  = spv 

      return
      end

