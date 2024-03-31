      subroutine rotate
c..7-30-06
      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'cconst'

      real*8 alp,dhy,ahy,ues(kdm),vofj,omegnew,omeg2
      real*8 xip,xil,omegbar,fact,domeg,ajayold,ajaynew,ajayedge
c..critical angular velocity
c      real*8 omegcrit

      real*8 srot(kdm),rotf(kdm),ajnew(kdm)

      integer*4 k
c--------------------------------------------------------------
c..   Eddington-Sweet currents (Tassoul,p. 197 approx.)
      do k = 2, kk
         g(k) = grav * xm(k) / r(1,k)**2
         alp  = omeg(k)**2 * r(1,k) / g(k)
         dhy  = abs( x(ndim-2,k+1) - x(ndim-2,k) )
         ahy  =      x(ndim-2,k+1) + x(ndim-2,k)
         ues(k)  = tl(1,k)/xm(k)/g(k)*alp
c..   add ES velocity only in homogeneous, radiative
c..   hydrogen regions as advection
c         if( ues .gt. 0.0d0 .and. dhy .lt. 1.0d-3
c     1        .and. ahy .gt. 1.0d-2
c     2        .and. ic(k) .eq. 0 )then
c            h(k) = max1( h(k), ues(k))
c         endif
      enddo

c..define angular velocity and momentum
      do k = 1, kk+1
c..initial omega
         srot(k)  = omeg(k)
c..ajnew has mass factor dmi(k)
         ajnew(k) = omeg(k) * r(1,k)**2 * dmi(k)
      enddo

c..two options for rotation evolution
      if( mrot .eq. 1 )then
c..   no change for uniform rotation
c..   loop defines velocity of angular momentum migration
         do k = 2, kk+1
            vofj    = omeg(k-1) - omeg(k)
            vofj    = max( vofj, 0.0d0 )
            rotf(k) = dth(2)*vofj
         enddo
         rotf(1) = 0

c..   loop modifies omeg from 1 to kk
         ajnew(1) = 0.0d0
         do k = 2, kk+1
c..   convective mixing tends to uniform rotation
c..   back reaction on radius ignored
c..   envelope is assumed coupled to zone kk
            srot(k) = (srot(k) + rotf(k)*srot(k-1))
     1           /( 1.0d0 + rotf(k))
         enddo

         do k = 2, kk+1
            omegnew = srot(k)
            omeg2 = (omeg(k) - omeg(k-1))**2 * dth(2)
            ajnew(k) = srot(k) * r(1,k)**2 * dmi(k)
            write(*,'(i5,1p8e12.3)')k,omeg(k),omegnew,omegnew*r(1,k),
     1           ajnew(k),ajay(k)*dmi(k),rotf(k),omeg2
         enddo
         write(*,'(a5,8a12)')"k","omega","omnew","vrot","ajnew",
     1        "ajay*dmi","rotf","omeg**2 dt"

      elseif( mrot .eq. 2 )then

c..   loop modifies omeg from 1 to kk
         do k = 2, kk+1
c..   convective mixing tends to uniform rotation
c..   at a rate determined by velocity h(k)
c..   back reaction on radius ignored
c..   convective velocity, including Eddington-Sweet currents
c..   carry angular momentum outward
            if( k .le. kk+1 )then
               rotf(k) = h(k)*dth(2)/( r(1,k) - r(1,k-1) )
            else
c..   relax to assumed rigid rotation for envelope (put in fitenv?)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               rotf(k) = 0.5d0
            endif
c..   moments of inertia
            xip       = dmi(  k)*r(1,k  )**2 /sol
            xil       = dmi(k-1)*r(1,k-1)**2 /sol
c            omegbar   = (xip*omeg(k) +  xil*omeg(k-1) )/(xip + xil)
            omegbar   = (omeg(k)+omeg(k-1))/2.0d0
            fact      = rotf(k)/( 1.0d0 + rotf(k) )
            domeg     = ( omegbar - srot(k) )*fact
            omeg(k)   = srot(k) + domeg
            if( xil .gt. 0.0d0 )then
c..   conserve angular momentum pairwise
               omeg(k-1) = omeg(k-1) - xip/xil*domeg
            else
               omeg(k-1) = omeg(k)
            endif
            if( k .gt. kk-3 )then
            write(*,'(i5,1p9e12.3)')k,rotf(k),srot(k),omeg(k),domeg,
     1           h(k),dth(2),r(1,k),ajay(k),omeg(k)*r(1,k)**2
            endif
         enddo
         write(*,'(a5,8a12)')'k','rotf','srot','omeg','domeg',
     1        'h','dth','r(1,k)','akay(k)'
         k = kk+1
         if( omeg(k) .ne. 0.0d0 )then
            write(*,'(a20,2i5,1p9e11.3)')"CMIX: KMAX ",k,ic(k),
     1           srot(k), omeg(k),r(1,k),r(1,k)*omeg(k),
     2           omeg(k)*r(1,k)**2*dmi(k) ,dmi(k),dmh(k)
         endif
      endif

      ajayold = 0
      do k = 1, kk+1
         ajayold = ajayold + dmi(k)*srot(k)*r(1,k)**2
      enddo

c..   update angular momentum from new rotation
      do k = 1, kk+1
         ajay(k) = omeg(k)* r(1,k)**2
         vrot(k) = omeg(k)* r(1,k)
c         write(*,*)"vrot ",vrot(k)
         rotshear(k) = max(0.0d0,sqrt(vrot(k)**2 + ues(k)**2))
      enddo

      ajaynew = 0
      do k = 1, kk+1
         ajaynew = ajaynew + dmi(k)*ajay(k)
      enddo
      ajayedge =  omeg(kk+1)*r(1,kk+1)**2*dmi(kk)


c      if( ajaynew .ne. 0.0d0 )then
cc         write(*,'(a10,1p9e11.3)')" J new ",ajayold,ajaynew,ajayedge,
c     1        ajaynew-ajayold,omeg(kk+1)*r(1,kk+1),omeg(kk+1) 
c         omegcrit =  sqrt( grav*xm(kk+1)/r(1,kk+1)**3 )
c         write(*,'(a15,i5,10(a8,1pe11.3))')'ROTATE: mrot',mrot,
c     1        'Jold',ajayold,'Jnew',ajaynew,
c     1        'delJ',ajaynew-ajayold,
c     2        'v(km/s)',omeg(kk+1)*r(1,kk+1)*1.0d-5,
c     2        'omsurf',omeg(kk+1),'omcrit',
c     4        omegcrit,'ratio',omeg(kk+1)/omegcrit,
c     5        'R',r(1,kk+1)
c         write(*,'(a15,1p9e12.3)')'ROTATE: omega ',
c     1        omeg(kk+1),omeg(kk),omeg(kk-1),
c     1        omeg(kk+1)*r(1,kk+1)**2,omeg(kk)*r(1,kk)**2,
c     3        omeg(kk-1)*r(1,kk-1)**2,
c     4        ajay(kk+1),ajay(kk),ajay(kk-1)
c      endif

      return
      end
