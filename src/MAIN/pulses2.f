      subroutine pulses2(nc)

c..   evaluates radial pulsation period
c..   Eddington, The Internal Constitution of The Stars, p. 186-188
c..   psi is (r+delta r)/r is r is unperturbed radius
c..   omega2 is circular frequence omega**2

c..   estimates initial values of the eigenfrequence from the
c..   sound travel time, and calculates up to nharm .le. nhdim
c..   wavefunctions, iterating to a zero derivative outer boundary
c..   wavefunctions normalized to 1 at r=0

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'czone'
      include 'cgen'
      include 'cconst'

      include 'compu'
      include 'ceoset'
      include 'ceos.h'
      include 'cburn'
      include 'cenv' !new
      include 'cadvec'
      real*8 xmsol
  
      real*8 gam(kdm),c2(kdm),pvbar(kdm),psi(kdm),delri(kdm)
      real*8 gamh(kdm),xmuor2(kdm),omega2,dtdv,factor,travel
      real*8 game(kdme),c2e(kdme),pvbare(kdme),psie(kdme),zetae(kdme)
      real*8 delrie(kdme),dtdve
      real*8 gamhe(kdme),xmuor2e(kdme)

c..   number of harmonics + 1
      integer*4 nhdim,nharm,itom
      parameter( nhdim = 10, itom=30 )
      real*8 omega(nhdim), om2guess(itom),zetaguess(itom)

c      real*4    xx(kdm),yy(kdm),xmin,xmax,ymin,ymax
c      real*4    xmin0,xmax0,ymin0,ymax0
c      real*4    xminf,xmaxf,yminf,ymaxf
c      real*4    tmass,sing

c      integer*4 izbu(ndim),inbu(ndim),ibu,idummy
c      character*5 cbu

      real*4   o2r3gm,must0,gam3,betta1,must0_1
      double precision deltaT_T,dT1dr,inertiaJ,dave,Knumerator
      double precision deltaT_T2,dLr2,dLr,depsdr,nucnu, nuclam
      parameter(nucnu = 1)
      integer*4 k, nc,  i, j, ii, nzeros
      character*2 marker
c      character*72 text
c      character*44 txt
c      character*2  txtxt
c      character*8  labl, labl1

c      character*12 cxlabel, cylabel
c      character*20 ctlabel
c      character*10  cnum

c      character*7  cmod
c      character*5 device

c      data linesty/1/,ipgflag/0/
c-----------------------------------------------------------------
c      write(6,*)"alphaml?,epsilon?"
c      read(5,*)alphaml, epsilon
     
c      open(2,file='pulses.in')

c..   input parameters for model modification from gentran.in
c..   dummy reads
c      read (2,*)text
c      read (2,*)text
c      read (2,*)text

c..   pgplot graphics output device (eg, /xwin, /cps)
c      read (2,*)txt,labl,device
c      labl1 = 'device'
c      if( labl .ne. labl1 )goto 2000

c..   model identity
c      read (2,*)txt,labl,cmod
c      labl1 = 'cmod'
c      if( labl .ne. labl1 )goto 2000

c..   lines or points?
c      read (2,*)txt,labl,linesty
c      labl1 = 'linesty'
c      if( labl .ne. labl1 )goto 2000

c..   force graphical limits?
c      read (2,*)txt,labl,xminf
c      labl1 = 'xminf'
c      if( labl .ne. labl1 )goto 2000

c      read (2,*)txt,labl,xmaxf
c      labl1 = 'xmaxf'
c      if( labl .ne. labl1 )goto 2000

c      if( xminf .ne. xmaxf )then
c         write(*,*)'      x override'
c      endif

c      read (2,*)txt,labl,yminf
c      labl1 = 'yminf'
c      if( labl .ne. labl1 )goto 2000

c      read (2,*)txt,labl,ymaxf
c      labl1 = 'ymaxf'
c      if( labl .ne. labl1 )goto 2000

c      if( yminf .ne. ymaxf )then
c         write(*,*)'      y override'
c      endif

c..   number of eigenvalues
c      read (2,*)txt,labl,nharm
c      labl1 = 'nharm'
c      if( labl .ne. labl1 )goto 2000
c      if( nharm .le. 0 .or. nharm .gt. nhdim )then
c         write(*,*)'pulses input: nhram, nhdim ',nharm,nhdim
c         stop' nharm out of range '
c      endif

      nharm = 3

c..   flag to show wavefunctions during iteration
c      read (2,*)txt,labl,iwave
c      labl1 = 'iwave'
c      if( labl .ne. labl1 )goto 2000
c      read (2,*)txt,labl,alphaml
c      labl1 = 'alphaml'
c      if( labl .ne. labl1 )goto 2000
c      read (2,*)txt,labl,epsilon
c      labl1 = 'epsilon'
c      if( labl .ne. labl1 )goto 2000
c      read (2,*)txt,labl,netsize
c      labl1 = 'netsize'
c      if( labl .ne. labl1 )goto 2000
c      read (2,*)text

      
c      close(2)

         open(51,file='period.dat')
         open(52,file='period1.dat')
         write(51,*)"mode  age soundtime  period  must0_1 kqa"
         write(52,*)"mode, period, age, kqa, logtem loglum"
c..   get the model, formatted read in tycho format
c      open(8,file=cmod)
c      rewind 8
c      nc = 1
c      call nreadf(note,nc)
c      close(8)

c      open(30,file='net.rc',status='old')
c      ibu = 1      
c..initializing qex and solarx
c 100  read(30,'(3i5,a5,0pf10.4,1pe12.4)',end=101)
c     1     idummy,izbu(ibu),inbu(ibu),cbu,qex(ibu),solarx(ibu)
c      lz(ibu)   = izbu(ibu)
c      ln(ibu)   = inbu(ibu)
c      cnuc(ibu) = cbu
c      write(*,'(3i5,a5)')ibu,izbu(ibu),inbu(ibu),cbu

c      if( izbu(ibu) .eq. 2 .and. inbu(ibu) .eq. 2 )goto 101
c      ibu = ibu + 1
c      goto 100
c 101  continue

c..   network size for this net.rc
c      netsize = ibu
c..   net.rc is consistent
c      read(30,'(10i5)')nucp

c      write(*,*)'net.rc network has netsize = ',netsize

c..   determine nonzero entries
c      mnucpg = 0
c      do j = 1, nucpg
c         if( nucp(j) .ne. 0 )then
c            mnucpg = mnucpg + 1
c         endif
c      enddo

c      write(*,*)'pulses1: net.rc exists and is left unchanged'
c      close(30)

c      call build(0)

c..   define mass coordinate
c      dmi(1)    = dmh(2)
c      dmi(kk+1) = dmh(kk+1)
c      do k = 2, kk
c         xm(k)  = xm(k-1) + dmh(k)
c         dmi(k) = 0.5*( dmh(k+1) + dmh(k) )
c      enddo
c      xm(kk+1) = xm(kk) + dmh(kk+1)
      xmsol = xm(kk+1)/sol

c      write(*,'(a8,a6,a10,0pf12.6,a10,i5,a8)')'Model  ',cmod,
c     1     'MASS', xmsol, 'nnuc', nnuc,'nuclei'

c      write(*,'( 3(a12,0pf12.6) )')
c     1     'envelop mass',dmh(kk+1)/sol,'r(env)/R',
c     1     r(1,kk)/r(1,kk+1),'R/Rsun',r(1,kk+1)/solrad


c..   renormalize abundances to full double precision
c      do k = 2, kk
c         sum = 0.0d0
c         do n = 1, nnuc
c            xa(n) = lz(n) + ln(n)
c            sum = sum + x(n,k)*xa(n)
c         enddo
c         if( abs(sum-1.0d0) .gt. 1.0d-2 )then
c            stop'nucleon check in pulses'
c         else
c            do n = 1, nnuc
c               x(n,k) = x(n,k)/sum
c            enddo
c         endif
c         ye(k)   = x(ndim,k)
c         write(*,*)k,ye(k)
c      enddo
      
c..   fill in thermodynamic and opacity variables
c      call state(2,kk,1)
      
c..   normalized radius for plotting
c      do k = 1, kk+1
c         rr(k) = r(1,k)/r(1,kk+1)
c      enddo
c      write(6,*)"line180",r(1,kk),rr(k),r(1,kk+1)

      dlr2 = 0.0d0
      deltaT_T2=0.0d0
cccccccccccccccccccccccccccc

c..   begin loop over harmonics
      mdotpuls = 0.0d0
      do i = 1, nharm
      
c..   initial guess for harmonic 
         ii = 1
c..   set up gamma at zone center
         do k = 2, kk
c            write(*,*)k,nc,p(nc,k),ev(k),et(k)
            dtdv    = - ( p(nc,k) + ev(k) )/et(k)
            gamh(k) = - v(nc,k)/p(nc,k)*( pv(k) + pt(k)*dtdv )
         enddo
         
c..   sound travel time
         travel = (r(nc,2)-r(nc,1))/dsqrt(gamh(2)*p(nc,2)*v(nc,2))
       
         do k = 2, kk-1
            travel    = travel +
     1           (r(nc,k+1)-r(nc,k))/dsqrt(gamh(k)*p(nc,k)*v(nc,k))
         enddo
cccccccccccccccc

        
c..   get envelope 
c..   l in comod now ccccccccccccccccccccccccccccccccccc
      call fitenv(nc)    
c    travel velocity calculation at the boundary between env and inner star
       
       travel    = travel +
     1    (zr(jmaxz-1)-r(nc,kk))/dsqrt(gamh(kk)*p(nc,kk)*v(nc,kk))
c       write(*,*)travel

c..   envelope
c..   map envelope to ascending index order with increasing radius
      do j = 1, jmaxz 
         k = jmaxz - j + 1
         dtdve    = - ( zp(j) + zev(j) )/zet(j)
         gamhe(k) = - zv(j)/zp(j)*( zpv(j) + zpt(j)*dtdve )
c         re(k) = zr(j)/r(1,kk+1)
      enddo         

      do k=2, jmaxz-1 ! index of z* , bigger j means inner radius
      j=jmaxz - k + 1   ! index of z* , bigger j means inner radius

         travel    = travel +
     1        (zr(j-1)-zr(j))/dsqrt(gamhe(k)*zp(j)*zv(j))
c         write(*,*)travel
      enddo
ccccccccccccccccccccccccccccc


 

c..   estimate eigenvalue frequency, increasing frequency for
c..   higher overtones (factor of i)
         om2guess(ii) = (dble(i)*pi/travel)**2

 1000    continue
c..   begin integration of wave equation
         omega2 = om2guess(ii)
         omega(i) = sqrt( omega2 )
        
c..   starting values (at center)
c..   psi = delta r /rzero
         psi(1)    = 1.0d0
         dtdv      = - ( p(nc,2) + ev(2) )/et(2)
         gam(1)    = - v(nc,2)/p(nc,2)*( pv(2) + pt(2)*dtdv )
         pvbar(1)  = p(nc,2)*v(nc,2)
         c2(1)     = gam(1)*pvbar(1)
         xmuor2(1) = grav*pi43/(p(nc,2)*v(nc,2)**2)
         factor    = 0.0d0
         delri(1)  = 0.25d0*( r(nc,2) - r(nc,1) )
c..   zeta = d(psi)/d(rzero)
         zeta(2)   = -0.1d0*((4.0d0/gam(1)-3.0d0)*xmuor2(1)
     1        + omega2/c2(1)) * 2.0d0 * delri(1)

ccccccccccccccccccccc
      dT1dr=0.0   !dT1/dr where T1=(deltaT)/To
      dLr=0.0	!delta Lr      
      inertiaJ=0.0  !integral(r**2.*psi**2.) over dM
 
      Knumerator=0.0  !numerator of Kappa, denominator is inertiaJ*omega2
  
      open(14,file='pulparameter.dat')
ccccccccccccccccccccccccccccccc
         do k = 2, kk-1
            psi(k)    = psi(k-1) + zeta(k)*( r(nc,k) - r(nc,k-1) )
            dtdv      = - ( p(nc,k+1) + ev(k+1) )/et(k+1)
            c2(k)     = 0.5d0*( gamh(k  )*p(nc,k  )*v(nc,k  )
     1           + gamh(k+1)*p(nc,k+1)*v(nc,k+1))
            pvbar(k)  = 0.5d0*( p(nc,k)*v(nc,k) + p(nc,k+1)*v(nc,k+1) )
            gam(k)    = c2(k)/pvbar(k)
            xmuor2(k) = grav*xm(k)/( r(nc,k)**3 * pvbar(k) )
            delri(k)  = 0.5d0*( r(nc,k+1) - r(nc,k-1) )
            factor    = delri(k)/r(nc,k)
     1           *(4.0d0-0.5d0*xmuor2(k)*r(nc,k)**2)
            zeta(k+1) = zeta(k)*(1.0d0-factor)/(1.0d0+factor)
     1           - psi(k)*( (4.0d0/gam(k)-3.0d0)*xmuor2(k)
     2           + omega2/c2(k) )*delri(k)/(1.0d0+factor)
cccccccccccccccccccccccccccccccc
        dT1dr=(1-gam(k-1))*(3.*zeta(k)+r(nc,k)*(zeta(k+1)
     &    -zeta(k-1))/delri(k)/2)  
ccc Add epsilon mechanism
        nuclam = 50.8d0 / (t(nc,k)/1.0d6)**0.333 - 2.0d0/3.0d0
        depsdr=nuclam*(3.*zeta(k)+r(nc,k)*(zeta(k+1)
     &    -zeta(k-1))/delri(k)/2) * nucnu*dT1dr
c        dT1dr=(1-gam(k-1))*(4.*zeta(k)+r(nc,k)*(zeta(k+1)*(1+factor)
c     &    -zeta(k))/delri(k))  
      deltaT_T=(1-gam(k-1))*(3.*psi(k)+r(nc,k)*zeta(k))  !T1
      dLr=(tl(nc,k)+tl(nc,k-1))/2.*(4.0*psi(k)-(akt(k-1)+akt(k))
     & /(ak(k-1)+ak(k))*deltaT_T*(t(nc,k-1)+t(nc,k))/2.
     & +4.0*deltaT_T+(r(nc,k)-r(nc,k-1))/(t(nc,k)-t(nc,k-1))*dT1dr)  !deltaLr
      if(ak(k-1).gt.0.) then
      Knumerator=Knumerator+(depsdr-dLr)*dT1dr*(r(nc,k)-r(nc,k-1))
      else  !to prevent 'NAN' 
c         write(6,*)"skipping ",k,"th shell "
      endif
      inertiaJ=inertiaJ+(r(nc,k)*psi(k))**2.*dmh(k) 

      if(k.eq.2) then  
c  the value at k=1 was NAN. We need this to put as a boundary condition in
c  the partial interation
         deltaT_T2=deltaT_T
         dLr2=dLr
      endif   

      write(14,*)deltaT_T,dT1dr,inertiaJ,dave,dLr,Knumerator
cccccccccccccccccccccccccccccccccccccc

         enddo
         psi(kk)      = psi(kk-1) + zeta(kk)*( r(nc,kk) - r(nc,kk-1) )
cccccccccccccccccccc
c NEW!!! 
c get the paramters at the boundary between the inner part and env
         c2(kk)     = 0.5d0*( gamh(kk  )*p(nc,kk  )*v(nc,kk  )
     1        + gamhe(2)*zp(jmaxz-1)*zv(jmaxz-1))

c jmaxz== corrsponds to kk in terms of radius
c jmaxz-1=> k=2

      pvbar(kk)  = 0.5d0*( p(nc,kk)*v(nc,kk) + zp(jmaxz-1)*zv(jmaxz-1) )
         gam(kk)    = c2(kk)/pvbar(kk)
         xmuor2(kk) = grav*xm(kk)/( r(nc,kk)**3 * pvbar(kk) )

      zetae(1)=zeta(kk)     
      psie(1)=psi(kk)
      c2e(1)=c2(kk)
      pvbare(1)=pvbar(kk)
      game(1)=gam(kk)
      xmuor2e(1)=xmuor2(kk)
      delrie(1)  = 0.5d0*( zr(jmaxz-1) - r(nc,kk-1) )
      factor= delrie(1)/zr(jmaxz)*(4.0d0-0.5d0*xmuor2e(1)*zr(jmaxz)**2)
      zetae(2)   =   zetae(1)*(1.0d0-factor)/(1.0d0+factor)
     1        - psie(1)*( (4.0d0/game(1)-3.0d0)*xmuor2e(1)
     2        + omega2/c2e(1) )*delrie(1)/(1.0d0+factor)
      
      dT1dr=(1-gam(kk-1))*(4.*zetae(1)+zr(jmaxz)*(zetae(2)
     &    -zeta(kk-1))/delrie(1)/2)
c               dT1dr=(1-gam(kk-1))*(4.*zetae(1)+zr(jmaxz)*(zetae(2)
c     &*(1+factor)-zeta(kk))/delrie(1))
      deltaT_T=(1-gam(kk-1))*(3.*psie(1)+zr(jmaxz)*zetae(1))    
      inertiaJ=inertiaJ+(zr(jmaxz)*psie(1))**2.*
     & (xm(kk+1)+zm(jmaxz)-xm(kk-1))
      dLr=(tl(nc,kk-1)+zl(jmaxz))/2.*(4.0*psi(kk)-(akt(kk-1)+akt(jmaxz))
     &/(ak(kk-1)+ak(jmaxz))*deltaT_T*(t(nc,kk-1)+ztem(jmaxz))/2.
     & +4.0*deltaT_T+(zr(jmaxz)-r(nc,kk-1))/(t(nc,kk)-t(nc,kk-1))*dT1dr)  !deltaLr
      Knumerator=Knumerator+dLr*dT1dr*(zr(jmaxz)-r(nc,kk-1))


      write(14,*)deltaT_T,dT1dr,inertiaJ,dave,dLr,Knumerator
c this omega is omega in Hansen & Kawaler p 366, eq 10.25
c      write(6,*)"omega",omega2*r(nc,kk+1)**3./grav/xm(kk+1)
      o2r3gm=omega2*r(nc,kk+1)**3./grav/xm(kk+1)

c now envelope calculation
      do k = 2, jmaxz-1
      j=jmaxz - k + 1   ! index of z* , bigger j means inner radius

         if(k.eq.2) then
         psie(k)    = psie(k-1) + zetae(k)*( zr(j) - r(nc,kk) )
         else
         psie(k)    = psie(k-1) + zetae(k)*( zr(j) - zr(j+1) )
         endif
         dtdve      = - ( zp(j-1) + zev(j-1) )/zet(j-1)
         c2e(k)     = 0.5d0*( gamhe(k  )*zp(j )*zv(j)
     1        + gamhe(k+1)*zp(j-1)*zv(j-1))
         pvbare(k)  = 0.5d0*( zp(j)*zv(j) + zp(j-1)*zv(j-1) )
         game(k)    = c2e(k)/pvbare(k)
         xmuor2e(k) = grav*(xm(kk+1)+zm(j))/( zr(j)**3 * pvbare(k) )

c zm(j) is negative value. xm(kk+1)+zm(j) equals the enclosed mass

         if(k.eq.2) then
         delrie(k)  = 0.5d0*( zr(j-1) - r(nc,kk) )
         else
         delrie(k)  = 0.5d0*( zr(j-1) - zr(j+1) )
         endif
       

         zetae(k+1) = zetae(k)*(1.0d0-factor)/(1.0d0+factor)
     1        - psie(k)*( (4.0d0/game(k)-3.0d0)*xmuor2e(k)
     2        + omega2/c2e(k) )*delrie(k)/(1.0d0+factor)

         dT1dr=(1-game(k-1))*(4.*zetae(k)+zr(j)*(zetae(k+1)
     &    -zetae(k-1))/delrie(k)/2)
c         dT1dr=(1-game(k-1))*(4.*zetae(k)+zr(j)*(zetae(k+1)*(1+factor)-
c     &    -zetae(k))/delrie(k))

      deltaT_T=(1-game(k-1))*(3.*psie(k)+zr(j)*zetae(k))    

      inertiaJ=inertiaJ+(zr(j)*psie(k))**2.*(zm(j)-zm(j+1))

      dLr=(zl(j+1)+zl(j))/2.*(4.0*psie(k)-(zakt(j+1)+zakt(j))/
     &(zak(j+1)+zak(j))*deltaT_T*(ztem(j+1)+ztem(j))/2.
     & +4.0*deltaT_T+(zr(j)-zr(j+1))/(ztem(j)-ztem(j+1))*dT1dr)  !deltaLr

      Knumerator=Knumerator+dLr*dT1dr*(zr(j)-zr(j+1))

      write(14,*)deltaT_T,dT1dr,inertiaJ,dave,dLr,Knumerator

      enddo

c  calculation for boundary values
      psie(jmaxz)= psie(jmaxz-1) + zetae(jmaxz)*( zr(1)-zr(2))
         factor    = delrie(k)/zr(j)*(4.0d0-0.5d0*xmuor2e(k)*zr(j)**2)
      dT1dr=(1-game(jmaxz-1))*(4.*zetae(jmaxz)+zr(1)*(zetae(jmaxz)
     &    -zetae(jmaxz-1))/(zr(1)-zr(2)))  ! approximately
c      dT1dr=(1-game(jmaxz-1))*(4.*zetae(jmaxz)+zr(1)*(zetae(jmaxz)*
c     &(1+factor)-zetae(jmaxz-1))/(zr(1)-zr(2)))  ! approximately
      deltaT_T=(1-game(jmaxz-1))*(3.*psie(jmaxz)+zr(1)*
     & zetae(jmaxz))

      inertiaJ=inertiaJ+zr(1)**2.*psie(jmaxz)**2.*(zm(1)-zm(2))

      dLr=(zl(2)+zl(1))/2.*(4.0*psie(jmaxz)-(zakt(2)+zakt(1))/
     & (zak(2)+zak(1))*deltaT_T*(ztem(2)+ztem(1))/2.
     & +4.0*deltaT_T+(zr(1)-zr(2))/(ztem(1)-ztem(2))*dT1dr)  !deltaLr
      Knumerator=Knumerator+dLr*dT1dr*(zr(1)-zr(2))

      write(14,*)deltaT_T,dT1dr,inertiaJ,dave,dLr,Knumerator   

      close(14)
    
      Knumerator=(deltaT_T*dLr-deltaT_T2*dLr2-Knumerator)/inertiaJ
c    now kappa before divided by omega2
c the other big changes are made from line 660
cccccccccccccccccccccccccccccccccccc

c..   save zeta at "surface" (scaled by R) for iteration on omega
c      zetaguess(ii) = zeta(kk)*r(nc,kk)
      zetaguess(ii) = zetae(jmaxz)*zr(1)

c         write(*,'(2(a15,i5),2(a15,1pe12.3)))')
c     1        'eigenvalue = ',i,' iteration = ',ii,'omega2',omega2,
c     2        'zeta*R',zeta(kk)*r(nc,kk)
c     2        'zeta*R',zetae(jmaxz)*zr(1)

c         write(*,'(a5,8a12)')'k','r/R','psi','zeta*R'
c..   check convergence
c         do k = kk-1, kk
c            write(*,'(i5,1p8e12.3)')k,rr(k),psi(k),zeta(k)*r(nc,kk)
c         enddo

c..   do we have the right overtone?
         nzeros = 0
         do k = 2, kk-1
            if( psi(k)*psi(k+1) .lt. 0 )then
c..   sign change
               nzeros = nzeros + 1
            endif
         enddo
cccccccc envelope contribution
         do k = 1, jmaxz-1
            if( psie(k)*psie(k+1) .lt. 0 )then
c..   sign change
               nzeros = nzeros + 1
            endif
         enddo
ccccccccccccccccccc


c         write(*,*)'nzeros ',nzeros,' i-1 ',i-1,
c     1        ' zeta*R ',zetaguess(ii)
            
         if( ii .eq. 1 )then
c..first attempt
            ii = 2
            if( nzeros .le. i-1 )then
c..increase wiggle
               om2guess(ii) = omega2 * 1.1d0
            else
c..decrease wiggle
               om2guess(ii) = omega2 * 0.9d0
            endif

         else
            ii = ii+1
c..usual cycle
            if( nzeros .lt. i-1 )then
c..increase wiggle
               om2guess(ii) = omega2 * 1.1d0
            elseif( nzeros .gt. i-1 )then
c..decrease wiggle
               om2guess(ii) = omega2 * 0.9d0
            else
c..try to converge
               om2guess(ii) = dabs(om2guess(ii-1)-zetaguess(ii-1) *
     1              (om2guess(ii-1)-om2guess(ii-2))/
     2              (zetaguess(ii-1)-zetaguess(ii-2)))
            endif
c            write(*,*)'OMEGA GUESS ',om2guess(ii),ii
         endif



c         if( ipgflag .eq. 0  )then
c            ipgflag = 1
c..   scratch file for character conversion
c            open(10,file='scratch.pg')

c            IF ( pgbeg(0,device,1,1) .NE. 1 ) STOP' pgbeg error'

c..   no query for device
c            call PGASK (.FALSE.)
c..   roman font
c            call PGSCF(2)
c..   font scaling (size)
c            call PGSCH(1.0)
c..   new page
c            call PGPAGE
c            do k=1,kk
c               rr(k)=r(nc,k)/r(nc,kk+1)
c            enddo
c            do k=1,jmaxz
c               j=jmaxz-k+1
c               re(k)=zr(j)/r(nc,kk+1)
c            enddo   
c            xmin = rr(1)
c            xmax = rr(kk+1)
      
c             xmax = re(jmaxz)  !new
c            ymin = 0.0d0
c            ymax = 0.0d0
c            do k = 2, kk

c               sing = psi(k)
c               ymin = amin1(ymin,sing)
c               ymax = amax1(ymax,sing)
c            enddo
ccccccccccccccccccc
c      do k = 2, jmaxz
c         write(21,*)re(k),psie(k)
c         sing = psie(k)
c         ymin = amin1(ymin,sing)
c         ymax = amax1(ymax,sing)
c      enddo
cccccccccccccccccccccc

c            cylabel = 'psi=del(r)/r(0)'

c            cxlabel = 'r/R'
c..   override x
c            if( xminf .ne. xmaxf )then
c               xmin = xminf
c               xmax = xmaxf
c            endif
c..   override y
c            if( yminf .ne. ymaxf )then
c               ymin = yminf
c               ymax = ymaxf
c            endif

c..   adjust for more aesthetic margins
c            xmin0 = xmin - (xmax-xmin)*0.05
c            xmax0 = xmax + (xmax-xmin)*0.05
c            ymin0 = ymin - (ymax-ymin)*0.05
c            ymax0 = ymax + (ymax-ymin)*0.05

c..   window
c            call PGSWIN(xmin0,xmax0,ymin0,ymax0)
c..   set color
c            call PGSCI(1)
c..   tickmarks on x(bottom=B,top=C) and y(left)
c            call PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
c..   Left label
c            call PGMTXT('L',2.0,0.5,0.5,cylabel)
c..   bottom label
c            call PGMTXT('B',2.0,0.5,0.5,cxlabel)

c..   title
c            call PGSCI(1)
c            tmass   = xm(kk+1)/sol
c            call ftoc(tmass,cnum)
c            ctlabel = 'M='//cnum
c            call PGTEXT(xmin0+0.05*(xmax-xmin),
c     1           ymax0+0.05*(ymax-ymin),ctlabel)
           
c            call PGTEXT(xmin0+0.95*(xmax-xmin),
c     1           ymax0+0.05*(ymax-ymin),cmod)



c..   zero line to emphasize zeros in wavefunction
c            call PGSCI(15)  
c            xx(1) = xmin0
c            xx(2) = xmax0
c            yy(1) = 0.
c            yy(2) = 0.
c            call PGLINE(2,xx,yy)

c..end of graphics initialization
c         endif

c         if( iwave .ne. 0 )then
c           call PGSCI(i+1)
c            do k = 1, kk
c               xx(k) = rr(k)
c               yy(k) = psi(k)
c            enddo
        
c            if( linesty .eq. 0 )then
c               call PGLINE(kk,xx,yy)
c            else
c               call PGPT(kk,xx,yy,20)
c            endif
c      do k = 1, jmaxz
c         xx(k) = re(k)
c         yy(k) = psie(k)
c         write(6,*)xx(k),yy(k)
c      enddo

c..   draw curve
c            if( linesty .eq. 0 )then
c               call PGLINE(jmaxz,xx,yy)
c            else
c               call PGPT(jmaxz,xx,yy,20)
c            endif
c..   label curve
c            sing = omega(i)
c            call ftoc(sing,cnum)
c            call PGTEXT(xx(jmaxz),yy(jmaxz),cnum)
c         endif

c         if( ii .lt. itom .and.
c     1        abs( zetaguess(ii)-zetaguess(ii-1)) .gt. 1.0d-4)then
c            goto 1000
c         else


c            do k = 2, kk

c               sing = psi(k)
c               ymin = 1000*amin1(ymin,sing)
c               ymax = 1000*amax1(ymax,sing)
c            enddo
ccccccccccccccccccc
c      do k = 2, jmaxz
c         write(21,*)re(k),psie(k)
c         sing = psie(k)
c         ymin = amin1(ymin,sing)
c         ymax = amax1(ymax,sing)
c      enddo
cccccccccccccccccccccc


c..we have convergence
cccccccccccccccccccccc
c New part

 
c      write(6,*)"####kga", Knumerator/2./omega2
c      write(6,*)"####1/kqa in days",2*pi/(Knumerator/2.0/omega2)/24/3600

      mdotpuls = mdotpuls+4.0d0*Knumerator*r(nc,kk+1)
     1           /(2.0d0*grav*xm(kk+1)) * inertiaJ

c      write(*,*)mdotpuls

c must0 is a boundary condition calculation from Hansen & Kawaler, p364, 
c     eq 10.14
      must0=(4.0-3.0*gamhe(jmaxz)+o2r3gm)*psie(jmaxz)-gamhe(jmaxz)*
     1 zetae(jmaxz)*r(nc,1+kk)
      dtdve=-( zp(1) + zev(1) )/zet(1)
         gam3=1.0-zv(1)/ztem(1)*dtdve
c must0_1 is a general boundary condition from J.P. Cox, p80,eq8.15
c  betta1 is Pg/(Pg+Prad)
        betta1=(zp(1)-ztem(1)**4.*arad/3.)/zp(1)
        must0_1=psie(jmaxz)*(o2r3gm*betta1-
     & 3.*gamhe(jmaxz)+4.*betta1+12.*
     & (gam3-1.)*(1.-betta1))/(gamhe(jmaxz)-4.*(gam3-1.)*(1.-betta1))
     & -r(nc,kk+1)*zetae(jmaxz)

c calculates how close to the instability strip
c      write(6,*)"close to instability? ",48.40-12.03*log10(ztem(1))-
c     & log10(tl(nc,kk+2)/sollum)

c to see that ztem(1) is diffrent from tl(nc,kk+2)
c      write(6,*)"TEMP, LUM",log10(ztem(1)),log10(zl(1)/sollum)
c      write(6,*)"TEMP, LUM",log10(t(nc,kk+2)),log10(tl(nc,kk+2)/sollum)
c      write(6,*) "general boundary condition, must0_1",must0_1   
c      write(*,*)'Boundary condition for d(dp)/dr=0',must0
ccccccccccccccccccccccccccccccc
c            write(*,'(a20,5x,2(a10,i5)/)')'converged values',
c     1           'harmonic',i,'iteration',ii
c         write(*,'(a20,2(1pe12.3,a5))')'2 pi/omega',
c     1        2.0d0*pi/dsqrt(omega2)/(24.0d0*3600.0d0), 'days',
c     2        2.0d0*pi/dsqrt(omega2),'sec'
c..   factor of 2 for round trip
c         write(*,'(a20,2(1pe12.3,a5))')'sound travel time',
c     1        2.0d0 * travel/(24.0d0*3600.0d0), 'days',
c     2        2.0d0 * travel,'sec'
c         write(*,'(a20,2(1pe12.3,a5))')'eddington period',
c     1        0.290d0*dsqrt ( v(nc,2)/(3.0d0*gam(1)-4.0d0 )), 'days',
c     2        0.290d0*dsqrt ( v(nc,2)/(3.0d0*gam(1)-4.0d0 ))
c     3        * 24.0d0 * 3600.0d0,'sec'
c         write(*,'(a20,1pe12.3)')'rho(c)/rho(ave)',
c     1        (pi43*r(nc,kk+1)**3)/(v(nc,2)*xm(kk+1))
c         write(*,'( 2(a15,1pe12.3) )') 'omega used', omega(i),
c     1        ' estimated ', dble(i)*pi/travel
c         write(*,'(a40/)')'Larger omega causes more wiggles '

cccccccccccccccccccccccc
         if(i.eq.1) then
            marker='e1'
            elseif(i.eq.2) then
               marker='e2'
              elseif(i.eq.3) then
                 marker='e3'
         endif   
c         write(52,*)marker, time,2.0d0 * travel/(24.0d0*3600.0d0),
c     1  2.0d0*pi/dsqrt(omega2)/(24.0d0*3600.0d0),must0_1,
c     1    Knumerator/2./omega2
c         write(51,*)marker, 2.0d0 * travel/(24.0d0*3600.0d0),
c     1  2.0d0*pi/dsqrt(omega2)/(24.0d0*3600.0d0),
c     1    Knumerator/2./omega2,time,must0_1
cccccccccccccccccccccccc        
c..   draw curve
c         call PGSCI(i+1)
c         open(58,file='zeta.dat')
c         do k = 1, kk
c            xx(k) = rr(k)
c            yy(k) = psi(k)
c            write(58,*)rr(k),psi(k)
c         enddo
          
c            if( linesty .eq. 0 )then
c               call PGLINE(kk,xx,yy)
c            else
c               call PGPT(kk,xx,yy,20)
c            endif
c      do k = 1, jmaxz
c         xx(k) = re(k)
c         yy(k) = psie(k)
c         write(6,*)xx(k),yy(k)
c            write(58,*)re(k),psie(k)
c      enddo
c      close(58)
c         if( linesty .eq. 0 )then
c            call PGLINE(kk,xx,yy)
c         call PGLINE(jmaxz,xx,yy)
c         else
c            call PGPT(kk,xx,yy,20)
c         call PGPT(jmaxz,xx,yy,20)
c         endif
c..   label curve
cccccccccccccc

c         sing = omega(i)
c         call ftoc(sing,cnum)
c         call PGTEXT(xx(jmaxz),yy(jmaxz),cnum)


c         endif

      enddo
c..   end of loop over eigenvalues


c      pause

      close(10)

c      call PGEND

c      stop'successful termination'

c 2000 continue
      close(51)
      close(52)
c      write(*,*)'pulses: error in  input file: pulses.in'
c      stop'pulses error'

      return
      end



