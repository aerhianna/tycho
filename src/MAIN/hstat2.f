      subroutine hstat2(k,itest,ks)

c------------------------------------------------------------------
c     implicit hydrostatics, damped and undampded hydrodynamics
c     assumes radiative diffusion (no flux limiting)
c..   4-1-08 wda
c     uses bisection if improvements alternate is sign
c
c     explicit snuc for low energy generation rates 
c     (needed for 200 Msol)
c
c..   last revised 5-16-05, rechecked join 12-18-05
c------------------------------------------------------------------
      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'csurface'
      include 'cphot'
      include 'cbug'
      include 'cgtr.h'
      include 'cconst'
      include 'cenv'
      include 'cburn'
      include 'ceoset'

      include 'caeps'

c..   henyey matrices
      integer*4 jsq,jyey2
      parameter( jsq = jyey*jyey, jyey2 = 2*jyey )

      real*8    ah(jsq), bh(jsq), ch(jsq), qlh(jsq), qli(jsq)
      real*8    qh(jyey), qkh(jyey), qki(jyey)
      real*8    ftest(jyey2)

c...  used in eos subroutine
      integer*4 ikd(kdm),ikt(kdm), intab(kdm)
      common/ineos/ikd,ikt,intab

      real*8    dqdv,dqdvm,deltr,delnv,fgtr,fact,fak,fako,pbar
      real*8    dp,dq,aom,part1,fg1,fg2,fg3,agm,dum,drm
      real*8    akitm,akivm,aki
      real*8    t3p,t3m
      real*8    amk,amrp,amlp,amt,amv,amlm,amrm
      real*8    aek,aerp,aelp,aet,aev,aelm,aerm
      real*8    ark,artp,arvp,arr,arl,artm,arvm
      real*8    alk,altp,alvp,alr,all,altm,alvm
      real*8    pplus,tplus,tbar,delp,ebar
c..   variable for short time step
      real*8    bold(kdm),bcons(kdm),cp(kdm)
      real*8    pdel,tdel,betta,bettb,vbar,fact1,fact2,fact3
      real*8    ombeta
      real*8    snucmin

      real*8 tconv,eff,omeff,tbaro,pbaro,delto,delpo,dnabo(kdm)
      real*8 dnma,dlmn,beta,ytot(kdm),dsc(kdm),lsc(kdm),alphasc
      real*8 radk,cpp

      integer*4 i, j, k, n, ks, nnn
      integer*4 itest, kdum

      real*8 akbar,tledding,delt,dvdt,dltdlp,dvdt2,dltdlp2
      real*8 dadvp,dadvm, dbdtp,dbdtm,dbdvp,dbdvm
      real*8 dcdtp,dcdtm,dcdvp,dcdvm, dddtp,dddtm,dddvp,dddvm
      real*8 dnabdtp,dnabdtm,dnabdvp,dnabdvm

      real*8 delv,etbar,evbar
      real*8    dndl,dndr,dndtp,dndvp,dndtm,dndvm
      real*8 lnuc(kdm),lgrav(kdm)
      real*8 dlnp,dlnt

      real*8 dmhfak
      real*8 fg4

      character*4 chvari(jyey2+1)
      data chvari/ '   R','   L','   T','   V','iter'/
      data alphasc/0.1d0/

      parameter( snucmin = 1.0d-5 )

      save
c-------------------------------------------------------------------
c     write(*,'(a25,2i5,1p8e12.3)')'ENTERING HSTAT2 ',
c    1     l,it,dnab(112),dnad(112),doux(112),dnabv(112)
c-----------------------------------------------

c..   arrays for henyey solution of staggered band matrix
      do i = 1, jyey
         qh(i) = 0.0d0
      enddo
      do i = 1,jsq
         ah(i) = 0.0d0
         bh(i) = 0.0d0
         ch(i) = 0.0d0
      enddo

c..   initialize  at n=2 level
      do k = 2, kk
         a(k) = pi4 * r(2,k)**2
         g(k) = grav * xm(k) /r(2,k)**2
      enddo

      do k = 2, kk
         ytot(k) = 0
         do n = 1, nnuc
            ytot(k) = ytot(k) + x(n,k)*(1.0d0+dble(lz(n)))/xa(n)
         enddo
      enddo

      if ( modes .eq. 2 )then
c..   force photosphere at kk+1 = ks (envelope integration)
         ks     = kk+1
         rphot  = r(2,ks)
         uphot  = u(2,ks)

      elseif (modes .eq. 0 )then
c..   force photosphere at kk = ks
         ks     = kk
         rphot  = r(2,ks)
         uphot  = u(2,ks)

      else
         write(*,'(/a25,i5,a15,i5)')'HSTAT(ERR): modes = ',
     1      modes,'for mode = ',mode
         write(*,*)
     1  'photosphere must be at edge of grid (modes = 2 or 0)'
         write(*,*)' for mode .ne. 0 (implicit dynamics)'
         stop'hstat 1'

      endif
c      t(1,2) = 5.0d8
c      r(1,2) = 1.1d7
c      tl(1,2) = 1.0d40

c..   keep ak(kk+1) fixed during iteration
      akt(kk+1) = 0.0d0
      akv(kk+1) = 0.0d0

      if( modes .eq. 2 )then
c..   assumed in join of envelope to interior
         dmi(kk) = dmh(kk)
      else
c..   dmh(kk+1) may have been updated, so update dmi(kk)
         dmi(kk) = 0.5d0*( dmh(kk) + dmh(kk+1) )
      endif

c..   central pseudoviscosity is 0.0d0, as are its derivatives
      dqdv = 0

c..   boundary condition at inner zone: Radius, Luminosity constant 
      do i = 1,jyey
         qki(i) = 0.0d0
         do j = i,jsq,jyey
            qli(j) = 0.0d0
         enddo
      enddo

      if( it .le. 1 )then
c..   initializations prior to iterative loop
         do k  = 1, kk+2
            do i = 1, 4
               do j = 1, 3
                  aa(i,j,k) = 0.0d0
               enddo
               dda(i,k)  = 0.0d0
               dda(i,k)  = 0.0d0
               aerr(i,k) = 0.0d0
            enddo
         enddo
c..   define pseudoviscosity at boundaries
         q(2,2)    = 0
         q(2,kk+1) = 0
      endif

c----------------------------------------------------------------------
c..   main computational loop on zone number k
c..   loop from inner boundary +1 to outer boundary -1
      do k = 2, kk-1

c..   pseudoviscosity for hydrodynamic modes
         dqdvm = dqdv
         if( dv(2,k+1) .ge. 0.0d0 
     1     .or. mode .eq. 1 .or. mode .eq. 2 )then
            q(2,k+1) = 0
            dqdv     = 0
         else
            deltr    = ( r(1,k+1) - r(1,k) )/dth(2)
            delnv    = ( v(2,k+1) - v(1,k+1) )/v(1,k+1)
            q(2,k+1) = qcons * ( delnv * deltr )**2 / v(1,k+1)
            dqdv     = 2.0d0 * qcons * delnv * ( deltr /v(1,k+1) )**2
         endif

c..   mass convervation equation

         fgtr = 2.0d0/( gam(k) + gam(k-1) )
         fact = fgtr*pi43*( r(2,k)**3 - r(2,k-1)**3 )/v(2,k)
         amk  = dmh(k) - fact
         amrp =-fgtr*a(k)/v(2,k)
         amlp = 0.0d0
         amt  = 0.0d0
         amv  = fact/v(2,k)
         amrm = fgtr*a(k-1)/v(2,k)
         amlm = 0.0d0

c..   first law (energy conservation equation)

         fak  = (tl(2,k) - tl(2,k-1))/dmh(k)
         fako = (tl(1,k) - tl(1,k-1))/dmh(k)
         pbar = 0.5d0*( p(1,k) + p(2,k) + q(1,k) + q(2,k) )

c..avoid implicit solution for low and negative nuclear heating
c..needed for 200 Msol in d burning; scales down for lower masses
         if( ss(k) .lt. snucmin*xm(kk)/sol )then
            sa(k) = 0.0d0
            sb(k) = 0.0d0
         endif
         if( snu(k) .lt. snucmin*xm(kk)/sol )then
            snua(k) = 0.0d0
            snub(k) = 0.0d0
         endif
c..approximate energy generation by a power law in iterative loop
         fact = 
     1        ss(k) *(t(2,k)/t(1,k))**sa(k)  *(v(2,k)/v(1,k))**sb(k)
     2        + 
     3        snu(k)*(t(2,k)/t(1,k))**snua(k)*(v(2,k)/v(1,k))**snub(k)
c..   second order in time
c..   s(6,k) is actual source used in dynam and hstat
c..   s(3,k) is past value of s(4,k) + s(5,k) + s(7,k)
         s(6,k) = (1.0d0 - alphae)*s(3,k) + fact*alphae
         st(k) = ss(k)*sa(k)/t(1,k) + snu(k)*snua(k)/t(1,k)
         sv(k) = ss(k)*sb(k)/v(1,k) + snu(k)*snub(k)/v(1,k)

c..   compositional change would be included here with e2-e1 form
c..   aek  = et(k)*dt(2,k) + ev(k)*dv(2,k)
ccccccccc

c        if( modec .eq. 1 )then
c..   ignore compositional effects for schwarzschild case
c           aek  = et(k)*dt(2,k) + ev(k)*dv(2,k)
c    1           + pbar   * dv(2,k)
c    2           + dth(2)*efi(k) * fak  * alphae
c    3           + dth(2)*efi(k) * fako * (1.0d0 - alphae)
c    4           - dth(2)*efi(k) * s(6,k)
c    5           - dth(2)*efi(k) * epsnuc(k)
c           aelp =   dth(2)*efi(k)/dmh(k)*alphae
c           aelm = - dth(2)*efi(k)/dmh(k)*alphae
c           aerp = 0
c           aerm = 0
c           aet  = et(k) + 0.5d0*pt(k)*dv(2,k)
c    1           - dth(2)*efi(k)*st(k)
c           aev  = ev(k) + 0.5d0*pv(k)*dv(2,k) + pbar
c    1           + 0.5d0*dqdv *dv(2,k)
c    2           - dth(2)*efi(k)*sv(k)

         if( mode .eq. 1 )then

            aek  = e(2,k) - e(1,k)
     1           + pbar   * dv(2,k)
     2           + dth(2)*efi(k) * fak  * alphae
     3           + dth(2)*efi(k) * fako * (1.0d0 - alphae)
     4           - dth(2)*efi(k) * s(6,k)
     5           - dth(2)*efi(k) * epsnuc(k)
            aelp =   dth(2)*efi(k)/dmh(k)*alphae
            aelm = - dth(2)*efi(k)/dmh(k)*alphae
            aerp = 0.0d0
            aerm = 0.0d0
            aet  = et(k) + 0.5d0*pt(k)*dv(2,k)
     1           - dth(2)*efi(k)*st(k)
            aev  = ev(k) + 0.5d0*pv(k)*dv(2,k) + pbar
     1           + 0.5d0*dqdv *dv(2,k)
     2           - dth(2)*efi(k)*sv(k)

        elseif( mode .eq. 2)then
c.. for low dth tests  -  only rad luminosity
c..

c.. short-time version
c            aek  = e(2,k) - e(1,k)
c..better for round-off ? 
            alphae = 1.0d-6

           aek  = et(k) * dt(2,k) + ev(k) * dv(2,k)
     1           + pbar   * dv(2,k)
     2           + dth(2)*efi(k) * fak  * alphae
     3           + dth(2)*efi(k) * fako * (1.0d0 - alphae)
     4           - dth(2)*efi(k) * s(6,k)
     5           - dth(2)*efi(k) * epsnuc(k)
            aelp =   dth(2)*efi(k)/dmh(k)*alphae
            aelm = - dth(2)*efi(k)/dmh(k)*alphae
            aerp = 0.0d0
            aerm = 0.0d0
            aet  = et(k) + 0.5d0*pt(k)*dv(2,k)
     1           - dth(2)*efi(k)*st(k)
            aev  = ev(k) + 0.5d0*pv(k)*dv(2,k) + pbar
     1           + 0.5d0*dqdv *dv(2,k)
     2           - dth(2)*efi(k)*sv(k)

c           aek  = e(2,k) - e(1,k)
c           aek  = et(k) * dt(2,k)
c     1           + pbar   * dv(2,k)
c     2           + dth(2)*efi(k) * fak 
c            aelp =   dth(2)*efi(k)/dmh(k)
c            aelm = - dth(2)*efi(k)/dmh(k)
c            aerp = 0
c            aerm = 0
c            aet  = et(k) + 0.5d0*pt(k)*dv(2,k)
c            aev  = ev(k) + 0.5d0*pv(k)*dv(2,k) + pbar
c     1           + 0.5d0*dqdv *dv(2,k)

c        write(*,'(i5,1p12e12.4)')k,aek,aelp*dtl(k+1),
c    1     aelm*dtl(k),aet*dt(2,k),aev*dv(2,k),
c    1     fak*dth(2)*efi(k),dtl(k)/tl(2,k),dr(k)/r(2,k),
c    3     dt(2,k)/t(2,k),dv(2,k)/v(2,k),tl(2,k),dtl(k)



         else

c.. long-time version
            aek  = e(2,k) - e(1,k)
     1           + pbar   * dv(2,k)
     2           + dth(2)*efi(k) * fak  * alphae
     3           + dth(2)*efi(k) * fako * (1.0d0 - alphae)
     4           - dth(2)*efi(k) * s(6,k)
     5           - dth(2)*efi(k) * epsnuc(k)
            aelp =   dth(2)*efi(k)/dmh(k)*alphae
            aelm = - dth(2)*efi(k)/dmh(k)*alphae
            aerp = 0.0d0
            aerm = 0.0d0
            aet  = et(k) + 0.5d0*pt(k)*dv(2,k)
     1           - dth(2)*efi(k)*st(k)
            aev  = ev(k) + 0.5d0*pv(k)*dv(2,k) + pbar
     1           + 0.5d0*dqdv *dv(2,k)
     2           - dth(2)*efi(k)*sv(k)

         endif
         amk  = amk*etar
         aek  = aek*etar

c===============================================================
c..   recurrsion coeficients at zone center: (R,L,T,V,R,L)
         ah(1) = amrp
         ah(2) = aerp
         ah(3) = amlp
         ah(4) = aelp

         bh(1) = amrm
         bh(2) = aerm
         bh(3) = amlm
         bh(4) = aelm

         ch(1) = amt
         ch(2) = aet
         ch(3) = amv
         ch(4) = aev

         qh(1) = amk
         qh(2) = aek

         aerr(1,k) = amk
         aerr(2,k) = aek

         call solve(ah,bh,ch,qh,qli,qlh,qki,qkh,jyey,jsq,k)

         call saver(jyey,jsq,k,1,4,3,kk+1,aa,qkh,qlh)

c..   momentum (Euler) equation: GTR hydro
         dq    = - q(2,k) + q(2,k+1)
         dp    = - p(2,k) + p(2,k+1) + dq
         aom   = a(k) / dmi(k)
         part1 = aom * dp
         fg1   = grav*gm(k)/r(2,k)**2
         fg2   = part1*gam(k)/grw(k)

         if( k .eq. 2 )then
            pbar = p(2,2 ) + q(2,2 )
         else
            pbar = 0.5d0*(p(2,k+1) + p(2,k) + q(2,k+1) + q(2,k))
         endif
         fg3   = pi4 * grav * pbar * r(2,k) / crad2
c..   centrifugal acceleration
         if( mrot .ne. 0 )then
            fg4 = 0.6666666667d0*omeg(k)**2*r(2,k)
         else
            fg4 = 0.0d0
         endif

         du(k) = -fg1 - fg2 - fg3 + fg4

         agm   = aom * gam(k) / grw(k)

         ark  = -du(k)
         artp =  agm*pt(k+1)
         arvp =  agm*( pv(k+1) + dqdv  )
         arr  =( agm*dp/r(2,k) - fg1/r(2,k) )*2.0d0
     1        -fg4/r(2,k)

         arl  =  0.0d0
         artm = -agm*pt(k)
         arvm = -agm*( pv(k)   + dqdvm )

c..   mode = 1, strict hydrostatic
c.....implicit hydrodynamics modes..........................
c     mode = 2, damps kinetic energy with dth e-folding time
c     mode = 3, first order backward differencing hydrodynamics
c     (mildly damped)
c     for no damping,  mode = 0 is more accurate
c............................................................
         dum = ( dth(2)*efi(k) )**2
         drm = r(2,k) - r(1,k)

         if( mode .eq. 1 )then
c..   implicit heavily damped hydro (e-fold = dth)
            ark = ark + drm/dum
            arr = arr + 1.0d0/dum

         elseif( mode .eq. 2 )then
c..   implicit mildy damped hydro (first order backwards difference)
            ark = ark + 2.0d0*( drm - u(1,k)*dth(2)*efi(k) ) / dum
            arr = arr + 2.0d0/dum

         elseif( mode .eq. 3 )then

            write(*,*)'hstat: error, mode = ',mode
            stop
         endif

         ark = ark*etar

c..   Luminosity equation: nabla formulation

c..   radiative nabla
        
         akbar    = 0.5d0*( ak(k+1) + ak(k) )
c..   better centering; avoids intermittent nonconvergence
         ombeta = 0.5d0*arad/3.0d0*( t(2,k+1)**4/p(2,k+1) 
     1        + t(2,k)**4/p(2,k) )
c         ombeta = 0.5d0*( t(2,k+1)**4 + t(2,k)**4 )*arad/3.0d0/
c     1            (0.5d0*(p(2,k+1) + p(2,k)))
         tledding = pi4*crad*grav*xm(k)*ombeta/akbar
c         if(dabs(0.25d0*tl(2,k)/tledding/dnrad(k)-1.0d0) .lt. 1.5d0)then
         dnrad(k) = 0.25d0 * tl(2,k) / tledding 
c         endif
         
c..   actual (structural) nabla
         tbar     = 0.5d0*( t(2,k+1) + t(2,k) )
         pbar     = 0.5d0*( p(2,k+1) + p(2,k) )
         delt     =         t(2,k+1) - t(2,k)
         delp     =         p(2,k+1) - p(2,k)
         if( delp .eq. 0.0d0 )then
            write(*,*)'hstat error: delp=0',
     1           ' at k,it,l = ',k,it,l
            stop'HSTAT: delp error'
c         elseif( delt .gt. 0.0d0 )then
c..temperature inversion
c            dnab(k)  = -dabs(delt*pbar/(delp*tbar))
c         else
          if(abs((delt*pbar/(delp*tbar)-dnab(k))/dnab(k)).le.1.0d-1)then
            dnab(k)  = delt*pbar/(delp*tbar)
          else
            dnab(k) = dnab(k)
          endif
         endif

c..   adiabatic nabla
         dvdt     = - et(k)/( p(2,k) + ev(k) )
         dltdlp   = p(2,k)/( t(2,k)*( pt(k) + pv(k)*dvdt ))
         dvdt2    = - et(k+1)/( p(2,k+1) + ev(k+1) )
         dltdlp2  = p(2,k+1)/( t(2,k+1)*( pt(k+1)+pv(k+1)*dvdt2 )) 
         dnad(k)  = 0.5d0*( dltdlp + dltdlp2)

c..   entropy gradient per zone
         delv     = v(2,k+1) - v(2,k)
         etbar    = 0.5d0*(  et(k+1) + et(k) )
         evbar    = 0.5d0*(  ev(k+1) + ev(k) )
c     z(k)     = etbar*delt +  (pbar + evbar)*delv
         z(k)     = e(2,k+1) - e(2,k) +  pbar*delv

c..   derivatives of dnab(k)....a=delt,b=tbar,c=pbar,d=delp
         dadvp =  0.0d0
         dadvm =  0.0d0
         dddtp =  pt(k+1)/delp
         dddtm = -pt(k  )/delp
         dddvp =  pv(k+1)/delp
         dddvm = -pv(k  )/delp
         dbdtp =  0.5d0/tbar
         dbdtm =  0.5d0/tbar
         dbdvp =  0.0d0
         dbdvm =  0.0d0
         dcdtp =  0.5d0*pt(k+1)/pbar
         dcdtm =  0.5d0*pt(k  )/pbar
         dcdvp =  0.5d0*pv(k+1)/pbar
         dcdvm =  0.5d0*pv(k  )/pbar
c..   rewritten to avoid division by zero at isothermal limit
c      dnabdtp = dnab(k)*( dadtp - dbdtp + dcdtp - dddtp )
c      dnabdtm = dnab(k)*( dadtm - dbdtm + dcdtm - dddtm )
         dnabdtp = dnab(k)*( - dbdtp + dcdtp - dddtp )
     1        + pbar/(tbar*delp)
         dnabdtm = dnab(k)*( - dbdtm + dcdtm - dddtm )
     1        - pbar/(tbar*delp)

         dnabdvp = dnab(k)*( dadvp - dbdvp + dcdvp - dddvp )
         dnabdvm = dnab(k)*( dadvm - dbdvm + dcdvm - dddvm )

c         dnabdtp = dnrad(k)*( - dbdtp + dcdtp - dddtp )
c     1        + pbar/(tbar*delp)
c         dnabdtm = dnrad(k)*( - dbdtm + dcdtm - dddtm )
c     1        - pbar/(tbar*delp)

c         dnabdvp = dnrad(k)*( dadvp - dbdvp + dcdvp - dddvp )
c         dnabdvm = dnrad(k)*( dadvm - dbdvm + dcdvm - dddvm )

         if( k .eq. 2 )then
c..   inner boundary condition 
            lnuc(1)  = 0.0d0
            lgrav(1) = 0.0d0
         endif
c..   gravitational luminosity
         lgrav(k) = lgrav(k-1) - ( e(2,k) - e(1,k)
     1        + 0.5d0*(p(2,k)+p(1,k)) * (v(2,k)-v(1,k)) )/dth(2)
     2        *dmh(k)
c..   nuclear luminosity
         lnuc(k)  = lnuc(k-1) + dmh(k)*( s(5,k)+ s(4,k) )

c..   subnab returns convective velocity hp(k)
         if( modec .ge. 0 )then
            call subnab(k,dndl,dndr,dndtp,dndvp,dndtm,dndvm,ytot)
         else
            hp(k) =  0
         endif

         if( modec .eq. 1 )then
           dnma = dnab(k) - dnad(k)
           dlmn = dnad(k) + doux(k) - dnab(k)
c           dnma = dnrad(k) - dnad(k)
c           dlmn = dnad(k) + doux(k) - dnrad(k)
           beta = rgas*ytot(k)*t(2,k)/v(2,k)/p(2,k)
           radk = 4.0d0*arad*crad*t(2,k)**3*v(2,k)/(3.0d0*ak(k))
           fact = beta*(8.0d0-3.0d0*beta)/
     1             (32.0d0 - 24.0d0*beta -beta**2)
           cpp = et(k) - (p(2,k)+ev(k))*pt(k)/pv(k)
           if( ic(k) .eq. 2 )then
c.. lsc is L(semiconvection)/L(rad) according to N. Langer,
c.. M. El Eid, K.J. Fricke, 1985, A&A, 145, 179-191 
             lsc(k) = alphasc*
     1             0.5d0*dnma/(dnab(k)*dlmn)*(dnma -fact*doux(k) )
             dsc(k) = alphasc*dnma/dlmn*radk*v(2,k)/(6.0d0*cpp)
cccccccccccccccccc
c          write(*,'(2i5,1p12e12.3)')k,ic(k),h(k),tl(2,k),
c    1   doux(k),dnma,dlmn,dnma/dlmn,beta,ytot(k),fact,lsc(k),
c    2   dsc(k)
           else
             lsc(k) = 0
             dsc(k) = 0
           endif
c          if( k .ge. 260 )then
c            write(*,'(2a5,12a12)')'k','ic','h','tl','doux',
c    1       'dnma','dlmn','dnma/dlmn','beta','ytot','fact',
c    2       'lsc/lr','dsc'
c            stop'hstat2'
c          endif
         endif


c..   ignores variation in dnad, doux over iteration

cccccccccccccccc
c..   time centered 4-22-06
c         eff = 0.5d0
c         omeff = 1.0d0 - eff
c         eff = 1.0d0
c         omeff = 0.0d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..time centered for large time steps, goes to zero for small
       eff = dth(2)/(dth(1) + 2.0d0*dth(2) )
c        eff = 0.5d0
c       eff = 1.0d0
       omeff = 1.0d0 - eff

c..   construct previous gradient
         tbaro = 0.5d0*( t(1,k+1) + t(1,k) )
         pbaro = 0.5d0*( p(1,k+1) + p(1,k) )
         delto =         t(1,k+1) - t(1,k)
         delpo =         p(1,k+1) - p(1,k)
         dnabo(k) = delto*pbaro/ ( delpo * tbaro )

         alk  = dnab(k) - dnabv(k)*eff - dnabo(k)*omeff

c         alk  = dnab(k) - dnabv(k)
c     1          - dnabv(k)* eff
c         write(*,*)dnab(k),dnabv(k),dnabo(k),dnrad(k)
         altp = dnabdtp - dndtp   *eff
         alvp = dnabdvp - dndvp   *eff
         alr  = 0.0d0   - dndr    *eff
         all  = 0.0d0   - dndl    *eff
         altm = dnabdtm - dndtm   *eff
         alvm = dnabdvm - dndvm   *eff
         alk  = alk*etar 
        
c..   recurrsion coeficients at zone boundary (T,V,R,L,T,V)
         ah(1) = artp
         ah(2) = altp
         ah(3) = arvp
         ah(4) = alvp

         bh(1) = artm
         bh(2) = altm
         bh(3) = arvm
         bh(4) = alvm

         ch(1) = arr
         ch(2) = alr
         ch(3) = arl
         ch(4) = all

         qh(1) = ark
         qh(2) = alk

         aerr(3,k) = ark
         aerr(4,k) = alk

         call solve(ah,bh,ch,qh,qlh,qli,qkh,qki,jyey,jsq,k)

         call saver(jyey,jsq,k,0,4,3,kk+1,aa,qki,qli)

      enddo

c------------------------------------------------end of main k loop


c..   now do "outer" boundary (join of envelope and interior)
      k = kk

c..   mass convervation equation (kk-1/2)

      fgtr = 2.0d0/( gam(k) + gam(k-1) )
      fact = fgtr*pi43*( r(2,k)**3 - r(2,k-1)**3 )/v(2,k)
      amk  = dmh(k) - fact
      amrp =-fgtr*a(k)/v(2,k)
      amlp = 0.0d0
      amt  = 0.0d0
      amv  = fact/v(2,k)
      amrm = fgtr*a(k-1)/v(2,k)
      amlm = 0.0d0

c..   first law = energy conservation equation (kk-1/2)
      fak  = (tl(2,k) - tl(2,k-1))/dmh(k)
      fako = (tl(1,k) - tl(1,k-1))/dmh(k)
      pbar = 0.5d0*( p(1,k) + p(2,k) + q(1,k) + q(2,k) )
      fact = s(4,k) + s(5,k) + s(7,k)
c..   s(6,k) is actual source used in dynam and hstat
c..   s(3,k) is past value of s(4,k) + s(5,k) + s(7,k)
      s(6,k) = alphae * fact  + (1.0d0 - alphae) * s(3,k)

c..   compositional change would be included here with e2-e1 form
c..   aek  = et(k)*dt(2,k) + ev(k)*dv(2,k)
      aek  = e(2,k) - e(1,k)
     1     + pbar   * dv(2,k)
     2     + dth(2)*efi(k) * fak  * alphae
     3     + dth(2)*efi(k) * fako * (1.0d0 - alphae)
     4     - dth(2)*efi(k) * s(6,k)
      aelp =   dth(2)*efi(k)/dmh(k)*alphae
      aelm = - dth(2)*efi(k)/dmh(k)*alphae
      aerp = 0
      aerm = 0
      aet  = et(k) + 0.5d0*pt(k)*dv(2,k)
     1     - dth(2)*efi(k)*st(k)*alphae
      aev  = ev(k) + 0.5d0*pv(k)*dv(2,k) + pbar
     1     - dth(2)*efi(k)*sv(k)*alphae
     4     + 0.5d0*dqdv *dv(2,k)
      amk  = amk*etar
      aek  = aek*etar

c=============================================================
c..   recurrsion coeficients at zone center: (R,L,T,V,R,L)
      ah(1) = amrp
      ah(2) = aerp
      ah(3) = amlp
      ah(4) = aelp

      bh(1) = amrm
      bh(2) = aerm
      bh(3) = amlm
      bh(4) = aelm

      ch(1) = amt
      ch(2) = aet
      ch(3) = amv
      ch(4) = aev

      qh(1) = amk
      qh(2) = aek

      aerr(1,k) = amk
      aerr(2,k) = aek

      call solve(ah,bh,ch,qh,qli,qlh,qki,qkh,jyey,jsq,k)

      call saver(jyey,jsq,k,1,4,3,kk+1,aa,qkh,qlh)

c..   Luminosity equations (kk)
      if( modes .eq. 2 )then
c..   envelope integration used to extrapolate T and P inward to kk-1/2
c..   convective state at kk is determined in envelope integration
c..   factor 1/2 is for extrapolation from kk-1/2 to kk, hence
c..   half of dmh(kk)
         dmhfak = 0.5d0

c..   long timestep
         all  = -dlntaudl -(enchidl*enab +enchi*enabdl)*dmh(k)*dmhfak
         alr  = -dlntaudr -(enchidr*enab +enchi*enabdr)*dmh(k)*dmhfak
         alk  = dlog( t(2,k) ) - tauln -enchi*enab*dmh(k)*dmhfak
     1        +all*(tl(2,k)-xtlum) +alr*(r(2,k)-rarget)
         altp = 0.0d0
         alvp = 0.0d0
         altm = 1.0d0/t(2,k)
         alvm = 0.0d0
         alk  = alk*etar
c..   for diagnostic purposes
         if(  enrad .gt. enab .and. enrad .gt. 0.0d0)then
            b(k) = (1.0d0-enab/enrad)*tl(2,k)/a(k)
         else
            b(k) = 0.0d0
         endif

c..   momentum (Euler) equation (kk)
c..   uses envelope join condition from fitenv integration
c..   for structure using homogeneous compsition,static ode's
c..   envelope join
         arl = -dlnpidl -enchidl* dmh(k)*dmhfak
         arr = -dlnpidr -enchidr* dmh(k)*dmhfak
         ark = dlog( p(2,k) ) - piln - enchi* dmh(k)*dmhfak
     1        +arl*(tl(2,k)-xtlum) + arr*(r(2,k)-rarget)

c..   set outer pressure, needed for diagnostics
         p(2,k+1) = p(2,k) * dexp( -enchi*dmh(k) )
         aom   = a(k)/( dmh(k) * dmhfak )
c..centered at kk with assumed symmetry kk+-1/2
c..this is the pressure gradient term
         fg1   = aom * epi *( dlog(p(2,k)) - piln
     1        -dlnpidl*(tl(2,k)-xtlum) -dlnpidr*(r(2,k)-rarget)  
     2        ) 
c..this is the usual gravity term
c         fg2   = grav * xm(k) / r(2,k)**2 
c..join version
         fg2 = (enchi
     1        + enchidl*(tl(2,k)-xtlum) + enchidr*(r(2,k)-rarget))
     2        * a(k) * epi
         du(k) = aom * epi * ark
c..   du(k) = fg1 - fg2
         artm = pt(k)/p(2,k)
         arvm = pv(k)/p(2,k)
         artp = 0.0d0
         arvp = 0.0d0

         ark = ark*etar

      else
c..   radiative zero boundary for modes .ne. 2
c..   use center value dmh(k) in aom at join (3-14-99)
         aom      = a(k)/dmh(k)
c..zero kk+1 for safety
         p(2,k+1) = 0.0d0
         t(2,k+1) = 0.0d0
         dpidl    = 0.0d0
         dpidr    = 0.0d0
         dtaudl   = 0.0d0
         dtaudr   = 0.0d0
         b(k)     = 0.0d0
         ic(k)    = 0
         aki      = ak(k)
         akitm    = akt(k)
         akivm    = akv(k)

c..   plane parallel gray atmosphere
         f(k) = arad*crad/3.0d0*t(2,k)**4/( aki/aom + 0.7d0 )

         alk  = tl(2,k)/a(k) - f(k)
c..   fact = dF/d(tau)
         fact = f(k)/(aki/aom + 0.7d0)
         altp = 0.0d0
         alvp = 0.0d0
         alr  = -2.0d0*tl(2,k)/a(k)/r(2,k) -fact*(-2.0d0*aki/aom/r(2,k))
         all  = 1.0d0/a(k)
         altm = -4.0d0*f(k)/t(2,k) -fact*(akt(k)/aom)
         alvm =                   -fact*(akv(k)/aom)        

         alk  = alk*etar

         dnad(k)  = dnad(k-1)
         dnabv(k) = dnabv(k-1)
         doux(k)  = 0.0d0
         h(k)     = 0.0d0
         hp(k)    = 0.0d0
c..   momentum (Euler) equation: GTR hydro
         dq    = - q(2,k) + q(2,k+1)
         dp    = - p(2,k) + p(2,k+1) + dq
c..uses dmh(k) = dmh(k+1), do effective dmi(k) = dmh(k) at k = kk
c..not   aom   = a(k) / dmi(k)
         aom   = a(k) / dmh(k)
         part1 = aom * dp
         fg1   = grav*gm(k)/r(2,k)**2
         fg2   = part1*gam(k)/grw(k)
         pbar = 0.5d0*(p(2,k+1) + p(2,k) + q(2,k+1) + q(2,k))
         fg3   = pi4 * grav * pbar * r(2,k) / crad2
         du(k) = -fg1 - fg2 - fg3
         agm   = aom * gam(k) / grw(k)
         ark  = -du(k)
         artp =  agm*pt(k+1)
         arvp =  agm*( pv(k+1) + dqdv  )
         arr  =( agm*dp/r(2,k) - fg1/r(2,k) )*2.0d0
         arl  =  0.0d0
         artm = -agm*pt(k)
         arvm = -agm*( pv(k)   + dqdvm )
c..   mode = 1, strict hydrostatic
c.....implicit hydrodynamics modes..........................
c     mode = 2, damps kinetic energy with dth e-folding time
c     mode = 3, first order backward differencing hydrodynamics
c     (mildly damped)
c     for no damping,  mode = 0 is more accurate
c............................................................
         dum = ( dth(2)*efi(k) )**2
         drm = r(2,k) - r(1,k)

         if( mode .eq. 2 )then
c..   implicit heavily damped hydro (e-fold = dth)
            ark = ark + drm/dum
            arr = arr + 1.0d0/dum
         elseif( mode .eq. 3 )then
c..   implicit mildy damped hydro (first order backwards difference)
            ark = ark + 2.0d0*( drm - u(1,k)*dth(2)*efi(k) ) / dum
            arr = arr + 2.0d0/dum
         endif

         ark = ark*etar

      endif

c=====================================
c..   recurrsion coeficients at zone boundary (T,V,R,L,T,V)
      ah(1) = artp
      ah(2) = altp
      ah(3) = arvp
      ah(4) = alvp

      bh(1) = artm
      bh(2) = altm
      bh(3) = arvm
      bh(4) = alvm

      ch(1) = arr
      ch(2) = alr
      ch(3) = arl
      ch(4) = all

      qh(1) = ark
      qh(2) = alk

      aerr(3,k) = ark
      aerr(4,k) = alk

      call solve(ah,bh,ch,qh,qlh,qli,qkh,qki,jyey,jsq,k)

      call saver(jyey,jsq,k,0,4,3,kk+1,aa,qki,qli)

c..   "outer" boundary done

c..   begin recursive sweep

      do k = 1, kk+2
         do j = 1, jyey2
c..   save previous change, zero for recursive sweep
            ddao(j,k) = dda(j,k)
            dda(j,k)  = 0.0d0
         enddo
      enddo


c..   initialize arrays for saving worst test values
      nnn = jyey2+1
      do j = 1, jyey2
         wtest(j)  = 0
         nwtest(j) = 1
      enddo

c..........................................back solution loop
      do j = 1, kk-1
         k = kk+1 - j
         dda(1,k) = aa(1,1,k)-aa(1,2,k)*dda(3,k+1)-aa(1,3,k)*dda(4,k+1)
         dda(2,k) = aa(2,1,k)-aa(2,2,k)*dda(3,k+1)-aa(2,3,k)*dda(4,k+1)
         dda(3,k) = aa(3,1,k)-aa(3,2,k)*dda(1,k  )-aa(3,3,k)*dda(2,k  )
         dda(4,k) = aa(4,1,k)-aa(4,2,k)*dda(1,k  )-aa(4,3,k)*dda(2,k  )
c..   staggered mesh
c..   drdel = aa(1,1,k) - aa(1,2,k)*dtdel - aa(1,3,k)*dvdel
c..   dldel = aa(2,1,k) - aa(2,2,k)*dtdel - aa(2,3,k)*dvdel
c..   dtdel = aa(3,1,k) - aa(3,2,k)*drdel - aa(3,3,k)*dldel
c..   dvdel = aa(4,1,k) - aa(4,2,k)*drdel - aa(4,3,k)*dldel
      enddo

c..........................................test extreme changes
      do k = 2, kk
c..   RADIUS
         ftest(1) = dda(1,k)/r(1,k)
         if( abs( ftest(1) ) .gt. abs( wtest(1) ) )then
            wtest(1)  = ftest(1)
            nwtest(1) = k
         endif
c..   LUMINOSITY
c..   small thermal energy content (dominates at small dth)
         ebar  = (et(k)*t(2,k)*dmh(k) + et(k+1)*t(2,k+1)*dmh(k+1))
     1        /dth(2) 
c..   small fractional change in luminosity from last converged value
         fak   = dabs(tl(1,k))
         fak   = dmax1(ebar,fak)
cccccccccccccccccccccccccccccccccc



         ftest(2) = dda(2,k)/fak

c        write(*,'(2i5,1p8e12.3)')k,it,tl(1,k),dda(2,k),fak,ebar,
c    1     ftest(2)
ccccccccccc


         if( abs( ftest(2) ) .gt. abs( wtest(2) ) )then
            wtest(2) = ftest(2)
            nwtest(2) = k
         endif
c..   TEMPERATURE
         ftest(3) = dda(3,k)/t(1,k)
         if( abs( ftest(3) ) .gt. abs( wtest(3) ) )then
            wtest(3) = ftest(3)
            nwtest(3) = k
         endif

c..   SPECIFIC VOLUME = 1/DENSITY
         ftest(4) = dda(4,k)/v(1,k)
         if( abs( ftest(4) ) .gt. abs( wtest(4) ) )then
            wtest(4) = ftest(4)
            nwtest(4) = k
            if( dda(4,k)+dv(2,k)+v(1,k) .lt. 0.0d0 )then
c..test for negative v(2,k) requires cumulative not just last change
               write(*,'(a16,a5,a10,i5,1pe12.3)')
     1              'HSTAT: negative ',chvari(4),
     1              ' zone ',nwtest(4),wtest(4)
               write(*,'(4(a15,1pe15.7))')'dda(4,k)',dda(4,k),
     1              'v(1,k)',v(1,k)
               write(*,'(4(a15,1pe15.7))')'dv(2,k)',dv(2,k),
     1              'v(2,k)',v(2,k),'sumv',dda(4,k)+dv(2,k)+v(1,k)
               write(*,'(4(a15,i15))')'ic ',ic(k)

               write(*,'(4(a15,1pe15.7))')'R(2,k)',r(2,k),
     1              'R(2,k+1)',r(2,k+1),'R(k+1)-R(k)',r(2,k+1)-r(2,k)
               write(*,'(4(a15,1pe15.7))')'U(2,k)',u(2,k),
     1              'U(2,k+1)',u(2,k+1),'U(k+1)-U(k)',u(2,k+1)-u(2,k)
               write(*,'(4(a15,1pe15.7))')'du(2,k)',du(k),
     1              'du(k+1)',du(k+1),'g(k)',g(k)
               write(*,'(4(a15,1pe15.7))')'P(2,k)',p(2,k),
     1              'P(2,k+1)',p(2,k+1),'P(k+1)-P(k)',p(2,k+1)-p(2,k)
               write(*,'(4(a15,1pe15.7))')'a(k)',a(k),
     1              'dmi(k)',dmi(k),'dmh(k)',dmh(k)
               write(*,'(4(a15,1pe15.7))')'T(2,k)',t(2,k),
     1              'tl(2,k)',tl(2,k),'tl(2,k+1)',tl(2,k+1)
               write(*,'(4(a15,1pe15.7))')'T(1,k)',t(1,k),
     1              'tl(1,k)',tl(1,k),'tl(1,k+1)',tl(1,k+1)
               write(*,'(4(a15,1pe15.7))')'s(3,k)',s(3,k),
     1              's(4,k)',s(4,k),'s(5,k)',s(5,k)
               n=nnuc-1
               write(*,'(a5,a10,1pe15.7,3(a15,1pe15.7))')
     1              cnuc(n), 'x(n,k)',x(n,k)*xa(n),
     2              'xold(n,k)',xold(n,k)*xa(n),
     3              'x(n,k+1)',x(n,k+1)*xa(n),
     4              'xold(n,k+1)',xold(n,k+1)*xa(n)

               write(*,'(a5,a10,1pe15.7,3(a15,1pe15.7))') cnuc(n),
     1              'xd(n,k)',xd(n,k)*xa(n),
     2              'xd(n,k+1)',xd(n,k+1)*xa(n)

               n=nnuc
               write(*,'(a5,a10,1pe15.7,3(a15,1pe15.7))') cnuc(n),
     1              'x(n,k)',x(n,k)*xa(n),'xold(n,k)',xold(n,k)*xa(n),
     2              'x(n,k+1)',x(n,k+1)*xa(n),
     3              'xold(n,k+1)',xold(n,k+1)*xa(n)

               nnn = 4

               goto 1000

            endif
         endif
      enddo

c..   find most restrictive zone for all conditions
      do n = 1, jyey2
         if( abs( wtest(n) ) .gt. resid ) then
            itest = nwtest(n)
         endif
      enddo

      if( it .gt. 9 .and. itest .ge. 1 .and. itest .le. kk )then
c..   if the most difficult zone has alternating signs on changes,
c..   use bisection
         if(  dda(1,itest)*ddao(1,itest) .lt. 0.0d0 .and.
     1        dda(2,itest)*ddao(2,itest) .lt. 0.0d0 .and.
     2        dda(3,itest)*ddao(3,itest) .lt. 0.0d0 .and.
     3        dda(4,itest)*ddao(4,itest) .lt. 0.0d0 )then
c            if(dabs(dda(1,itest)) .ge. dabs(ddao(1,itest)) .and.
c     1         dabs(dda(2,itest)) .ge. dabs(ddao(2,itest)) .and.
c     2         dabs(dda(3,itest)) .ge. dabs(ddao(3,itest)) .and.
c     3         dabs(dda(4,itest)) .ge. dabs(ddao(4,itest)) )then
               do k = 2,kk
                  do j = 1,4
                     dda(j,k) = dda(j,k)*0.5d0
                  enddo
               enddo
c            endif
            k = itest
            write(*,'(a20,a5,i5,a5,i5,1p8e14.6)')'using bisection,',
     1           'it =',it,'k =',itest,
     2           r(2,k),tl(2,k),t(2,k),1.0d0/v(2,k),
     3           (tl(2,k+1)-tl(2,k))/dmh(k),ak(k),
     4           dlog10(1.0d0/v(2,k)*(1.0d6/t(2,k))**3),h(k) 
         endif
      endif

c..   update to new guess at variables
      do k = 2, kk
         dr(k) = dr(k) + dda(1,k)
         if( mode .eq. 1 .or. mode .eq. 4 )then
            u(2,k)= dr(k)/( dth(2) * efi(k) )
         else
            u(2,k) = u(1,k) + du(k) * dth(2) * efi(k)
         endif
         dtl(k) = dtl(k)  + dda(2,k)
         dt(2,k)= dt(2,k) + dda(3,k)
         td(k)  = dda(3,k)
         dv(2,k)= dv(2,k) + dda(4,k)

         t(2,k)  = t(1,k) + dt(2,k)
         v(2,k)  = v(1,k) + dv(2,k)
         r(2,k)  = r(1,k) + dr(k)
         tl(2,k) = tl(1,k)+ dtl(k)

         if( k .gt. 2 )then
            s(1,k) = (f(k-1)*a(k-1) - f(k)*a(k))/dmh(k)
            s(2,k) = (b(k-1)*a(k-1) - b(k)*a(k))/dmh(k)
         else
            s(1,k) = (tl(2,1)   - f(2)*a(2))/dmh(2)
            s(2,k) = (b(1)*a(1) - b(2)*a(2))/dmh(2)
         endif

      enddo
c      t(2,2) = 5.0d8
c      r(2,2) = 1.1d7
c      tl(2,2) = 1.0d40
ccccccccccccccccccccccccccccccc
c      write(*,'(a12,i5,a12,1pe12.3)')'it = ',it,'dth(2)',dth(2)
c      write(*,'(2a5,13a12)')'k','ic','r','tl','t','v','h',
c     1  'dr','dtl','dt','dv','n-a-x','dnabv(k)','dnab(k)','doux(k)'
c      do k = 106,120
c      write(*,'(2i5,1p13e12.3)')k,ic(k),r(2,k),tl(2,k),t(2,k),v(2,k),
c     1   h(k),dr(k)/r(1,k),dtl(k)/tl(1,k),dt(2,k)/t(1,k),dv(2,k)/v(1,k),
c     2   dnab(k)-dnad(k)-doux(k),dnabv(k),dnab(k),doux(k)
c      enddo
c      if( it .ge. 4 ) stop'hstat2-aaa'
cccccccccccccccccc

      do n = 1, jyey2
c..   excessive variable change
         if( abs( wtest(n) ).gt. 0.8d0 )then
            write(3,'(a36,a5,a10,i5,1pe12.3)')
     1           'HSTAT: excessive variable change in',chvari(n),
     1           ' zone ',nwtest(n),wtest(n)
            write(*,'(a36,a5,a10,i5,1pe12.3)')
     1           'HSTAT: excessive variable change in',chvari(n),
     1           ' zone ',nwtest(n),wtest(n)
            nnn = n
            goto 1000
         endif
      enddo

c..........................................end of update loops

      if( it .eq. iter )then
         write(*,*)'hstat: too many iterations (',it,') no update '
         write(*,'(a12,i5,a12,a5,a12,1pe12.4)')
     1        'iteration:',it,'cause:',chvari(jyey2+1),'goal:',resid
c         write(*,'(5(a12,f12.3))')'mass',xm(kk+1)/sol,'Luminosity',
c     1        tl(2,kk)/sollum
         write(*,'(a20,4a12)')'variables',(chvari(n),n=1,jyey2)
         write(*,'(a20,4i12)')'k worst',(nwtest(n),n=1,jyey2)
         write(*,'(a20,1p4e12.4))')'worst values',(wtest(n),n=1,jyey2)
         write(*,'(a20,1p4e12.4))')'fraction',
     1        dr(nwtest(1))/r(1,nwtest(1)),
     2        dtl(nwtest(2))/tl(1,nwtest(2)),
     3        dt(2,nwtest(3))/t(1,nwtest(3)),
     4        dv(2,nwtest(4))/v(1,nwtest(4))
         write(*,'(a20,4i12)')'ic(kworst)',(ic(nwtest(n)),n=1,jyey2)
         write(*,'(a20,1p4e12.3)')
     1        'Ytot(kworst)',(x(ndim,nwtest(n)),n=1,jyey2)
         write(*,'(a20,1p4e12.3)')
     1        'Yold(kworst)',(xold(ndim,nwtest(n)),n=1,jyey2)
         write(*,'(a20,1p4e12.3)')
     1        'h(kworst)',(h(nwtest(n)),n=1,jyey2)
         write(*,'(a20,1p4e12.3)')
     1        'hp(kworst)',(hp(nwtest(n)),n=1,jyey2)

         write(*,'(a20,1p4e12.3)')
     1        'r,tl,t,v(1)',r(1,nwtest(1)),tl(1,nwtest(2)),
     2        t(1,nwtest(3)),v(1,nwtest(4))
         write(*,'(a20,1p4e12.3)')
     1        'r,tl,t,v(2)',r(2,nwtest(1)),tl(2,nwtest(2)),
     2        t(2,nwtest(3)),v(2,nwtest(4))
         write(*,'(a20,1p4e12.3)')
     1        'dr,dtl,dt,dv',dr(nwtest(1)),dtl(nwtest(2)),
     2        dt(2,nwtest(3)),dv(2,nwtest(4))

         write(3,'(2(a12,i5),a12,a5,a12,1pe12.4)')'zone',k,
     1        'iteration',it,'cause',chvari(jyey2+1),'goal',resid
         write(3,'(a20,4a12)')'variables',(chvari(n),n=1,jyey2)
         write(3,'(a20,4i12)')'worst zones',(nwtest(n),n=1,jyey2)
         write(3,'(a20,1p4e12.4))')'worst values',(wtest(n),n=1,jyey2)
         write(3,'(a20,1p4e12.4))')'fraction',
     1        dr(nwtest(1))/r(1,nwtest(1)),
     2        dtl(nwtest(2))/tl(1,nwtest(2)),
     3        dt(2,nwtest(3))/t(1,nwtest(3)),
     4        dv(2,nwtest(4))/v(1,nwtest(4))

      endif

      return

c------------------------------------panic exit
 1000 continue
      itest = -1

      write(*,'(20x,5a12)')cwtest,'resid'
      write(*,'(a20,4i12)')'worst zones',(nwtest(n),n=1,jyey2)
      write(*,'(a20,1p5e12.4))')'frac. change',(wtest(n),n=1,jyey2),
     1     resid

      k = nwtest(2)
      write(*,'(i5,1p8e12.3)')k,tl(1,k),tl(2,k),ak(k),
     1     dnab(k),dnad(k),dnrad(k)
ccccccccccccccccccc

      write(3,'(20x,5a12)')cwtest,'resid'
      write(3,'(a20,4i12)')'worst zones',(nwtest(n),n=1,jyey2)
      write(3,'(a20,1p5e12.4))')'worst values',(wtest(n),n=1,jyey2),
     1     resid

      if( nnn .eq. jyey2 + 2 )then
         kdum = kk-2
      else
         kdum = max0(nwtest(nnn)-2,2)
         kdum = min0(kdum,kk+5)
      endif


      write(3,'(a5,10a12)')'k',
     1     'R1','R2','T1','T2','V1','V2','TL1','TL2','P1','P2'
      do j = kdum,kdum+5
         write(3,'(i5,1p10e12.3)')j,r(1,j),r(2,j),t(1,j),t(2,j),
     1        v(1,j),v(2,j),tl(1,j),tl(2,j),p(1,j),p(2,j)
      enddo

      write(3,'(a5,8a12)')'k',cnuc(netsize-2),cnuc(netsize-1),
     1     cnuc(netsize),cnuc(lc12),cnuc(ln14),cnuc(lo16),'entropy'
      do j = kdum,kdum+5
         write(3,'(i5,1p8e12.3)')j,x(netsize-2,j),x(netsize-1,j),
     1        x(netsize,j),x(lc12,j),x(ln14,j),x(lo16,j),entropy(j)
      enddo

      write(3,'(a5,8a12)')'k','dnab','dnad','dnabv','doux','dnrad',
     1  'bv-a-x','n-a'
      do j = kdum,kdum+5
         write(3,'(i5,1p8e12.3)')j,dnab(j),dnad(j),dnabv(j),doux(j),
     1   dnrad(j),dnabv(j)-doux(j)-dnad(j),dnab(j)-dnad(j)
      enddo

      write(3,*)'OLD VALUES'
      write(3,'(a5,8a12)')'k',cnuc(netsize-2),cnuc(netsize-1),
     1     cnuc(netsize),cnuc(lc12),cnuc(ln14),cnuc(lo16)
      do j = kdum,kdum+5
         write(3,'(i5,1p8e12.3)')j,xold(netsize-2,j),xold(netsize-1,j),
     1        xold(netsize,j),xold(lc12,j),xold(ln14,j),xold(lo16,j)
      enddo

      write(3,'(2a5,8a12)')'k','ic',
     1     'h','hp','y','z','b','f','s4','s5'
      do j = kdum,kdum+5
         write(3,'(2i5,1p8e12.3)')j,ic(j),h(j),hp(j),y(j),z(j),
     1        b(j),f(j),s(4,j),s(5,j)
      enddo

      write(3,'(2a5,8a12)')'k','ic',
     1     's4','s5','st','sv','deut','t*st','v*sv'
      do j = kdum,kdum+5
         write(3,'(2i5,1p8e12.3)')j,ic(j),s(4,j),s(5,j),st(j),sv(j),
     1        x(1,j),t(1,j)*st(j),v(1,j)*sv(j)
      enddo


c      stop'hstat2 panic exit'

      return

      end

