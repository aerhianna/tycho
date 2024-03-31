      subroutine gtrs(p,q,e,v,r,u,xm,dmi,dmh,kk,n,newt,mode)

c..gtrr hydro structure(geometry) factors
c..modification (wda) of ken van riper (1979) ap. j. 232, 558
c..revised 5/15/90

      implicit none

      include 'dimenfile'

      real*8 pote,uave,ekin,eint,pqm,fpast,pq,f1,f2,pqp
      real*8 dpq,fact,dphi,bk,thek,gamh,dmhg

      integer*4 kk,n,newt,mode,k,kj

      real*8 p(2,kdm),q(2,kdm),e(2,kdm),r(2,kdm),u(2,kdm),v(2,kdm)
      real*8 xm(kdm),dmi(kdm),dmh(kdm)

      include 'cgtr.h'
      include 'cconst'

c..   input:   p,q,e,v,r,u,dmi,dmh,xm,kk,grav (dmi required 5/20/90)
c..   output:  grw,efi,gam,gm

c..   gtr hydro structure(geometry) factors
c..   modification (wda) of ken van riper (1979) ap. j. 232, 558
c---------------------------------------------------------------------- 
      pqp   = 0.0d0
      dpq   = 0.0d0

      if( newt .eq. 1 )then
c..   newtonian mechanics........................................
         if( r(n,1) .gt. 0.0d0 )then
            pote = dmi(1)*(grav * xm(1) / crad2) / r(n,1)
            uave = 0.25d0*(u(n,1)**2 + u(1,1)**2)
            ekin = uave/crad2*dmi(1)
         else
            pote = 0.0d0
            ekin = 0.0d0
         endif
         eint   = 0.0d0
         efi(1) = 1.0d0
         grw(1) = 1.0d0
         gam(1) = 1.0d0
         gm(1)  = xm(1)
         do k = 2, kk
            efi(k) = 1.0d0
            grw(k) = 1.0d0
            gam(k) = 1.0d0
            eint = eint + dmh(k)*e(n,k)/crad2
            if( mode .eq. 0 )then
               uave = 0.25d0*(u(n,k)**2 + u(1,k)**2)
            else
               uave = 0.5d0*u(n,k)**2
            endif
            ekin = ekin + uave/crad2*dmi(k)
            if( r(n,k) .gt. 0.0d0 )then
               pote = pote + dmi(k)*(grav * xm(k) / crad2) / r(n,k)
            endif
            gm(k) = xm(k) + eint + ekin - pote
         enddo

         return
      else
c..   general theory of relativity.................................
         pqm   = 0
         fpast = 0
         pq    = p(n,2) + q(n,2)
         do k = 2, kk
            f1 = e(n,k) + v(n,k)*( p(n,k)+q(n,k) )
            f2 = 1.0d0 + f1/crad2
            if(k .eq. 2 ) grw(1)   = f2
            if(k .eq. kk) grw(kk)  = f2
            if(k .gt. 2 ) grw(k-1) = 0.5d0*(f2 + fpast)
            fpast = f2
         enddo
c..   synchronize clocks with co-moving observer at edge (kk)
         efi(kk) = 1.0d0
         do kj = 2, kk
            k = kk - kj + 2
            if(k .lt. kk           ) pqp = p(n,k+1) + q(n,k+1)
            pq = p(n,k) + q(n,k)
            if(k .gt. 2            ) pqm = p(n,k-1) + q(n,k-1)
            if(k .eq. 2            ) dpq = pqp - pq
            if(k .eq. kk           ) dpq = pq - pqm
            if(k.ne.2 .and. k.ne.kk) dpq = 0.5d0*(pqp - pqm)
            fact     =   crad2*0.5d0*(grw(k) + grw(k-1))
            dphi     = - v(n,k)*dpq/fact
            efi(k-1) =   efi(k)*exp( - dphi )
         enddo
c..   solves for gamma and grav. mass together for better accuracy
         gam(1) = 1.0d0
         gm(1)  = xm(1)
         do k = 2, kk
c..   at k-1/2
            bk   = 1.0d0 + e(n,k)/crad2
            thek = grav*bk*dmh(k)/( crad2*r(n,k)*2.0d0 )
c..   this uses u at t(n), so for mode = 0
c..   it should be u(n,k) + dth(2)*du(k)
c..   to move from t(n-1/2) to t(n)
            if( mode .eq. 0 .or. mode .eq. 3 )then
               fact = 1.0d0 + (u(n,k)/crad)**2 + thek**2
            else
               fact = 1.0d0 + thek**2
            endif
            fact = fact - 2.0d0*grav*gm(k-1)/(crad2*r(n,k))
            fact = fact - 2.0d0*gam(k-1)*thek
            if(fact .le. 0.0d0 )then
               write(*,'(a30,i6,1p8e12.3)')'coordinate singularity at', 
     1              k,fact,gam(k-1),thek
               write(*,*)'newt,mode,n,u(n,k),crad,grav,bk,dmh,r(n,k)'
               write(*,*)newt,mode,n,u(n,k),crad,grav,bk,dmh(k),r(n,k)
               stop'gtrs.f'
            endif
            gam(k) = dsqrt( fact ) - thek
            gamh   = 0.5d0*(gam(k-1) + gam(k))
            dmhg   = dmh(k)*gamh*(1.0d0 + e(n,k)/crad2)
            gm(k)  = dmhg + gm(k-1)
         enddo
c..   gtr structure factors are now set up
         return
      endif

      end



