
      subroutine pay
c..7-30-06
c..pay routine for mixing according to W. Press

      implicit none

      include 'dimenfile'
      include 'cconst'
      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'ceoset'
      include 'cenv'

      real*8 v2up(kdm),v2down(kdm),bv2(kdm)
      real*8 dellr(kdm),hml(kdm),kh(kdm),omegc,sigmath(kdm),sigmam(kdm),
     1     sigmamup(kdm),sigmamdown(kdm),dumvar
      real*8     qfak(kdm),kv(kdm),epssqr(kdm)
      real*8 fact,fak,veq
      real*8 ybar
      integer*4 n

      integer*4 i,j,k,m
      integer*4 mm,reflec
c----------------------------------------------------------------
c..   test for no convection (so no wave sources)
      j = 0
      do k = 1, kk
         if( ic(k) .eq. 1 )j=j+1
      enddo

      if( j .eq. 0 )then
         do k = 2, kk
            vover(k)  = 0.0d0
            vunder(k) = 0.0d0
            dif(k)    = 0.0d0
         enddo

      else
c..   some convective sources
c..   dellr is distance between zone centers
c..   hml is pressure scale height
         do k = 2,kk
            if( k .lt. kk )then
               dellr(k) = 0.5d0*( r(1,k+1) - r(1,k-1))
            else
               dellr(k) =         r(1,k) - r(1,k-1)
            endif
            if( p(1,k) .gt. p(1,k+1) )then
               hml(k) = 0.5d0*dellr(k)*(p(1,k+1) + p(1,k))
     1              /(p(1,k) - p(1,k+1))
            else
c..   for safety if pressure inverted
               hml(k) = 1.0d1*dellr(k)
            endif
c..   over- and under-shooting speeds
            vover(k)  = h(k)
            vunder(k) = h(k)
         enddo
         hml(1)  = hml(2)
         hml(kk) = hml(kk-1)
ccccccccccccccccccccccccccc
         bv2(1) = 0.0d0
         do k = 2, kk
            bv2(k) = ( dnad(k) - dnab(k) )*g(k)**2/(
     4              0.5d0*(p(1,k)*v(1,k)+p(1,k+1)*v(1,k+1)) )
         enddo
c..downwhelling loop
c..get envelope convective velocity as boundary condition
c..use middle envelope in R,L space
         v2down(kk+1) = vvel(nvmax(3),3)**2
         if( bv2(kk) .lt. 0.0d0 )then
            v2down(kk+1) = -bv2(kk)/8.0d0 *hml(kk)**2
         else
            v2down(kk+1) = 0
         endif

         v2up(kk+1)   = v2down(kk+1)
         dellr(1) = dellr(2)
         dellr(kk+1) = dellr(kk)

         do j = 1, kk
            k = kk +1 -j
c            v2down(k) = v2down(k+1) -bv2(k)*dellr(k)*hml(k)

            v2down(k) = ( v2down(k+1) -bv2(k)*dellr(k)*hml(k)*2.0d0 )
     1           *hml(k)/( hml(k) + 16.0d0*dellr(k) )

            if( v2down(k) .lt. 0.0d0 )then
               v2down(k) = 0.0d0
            endif
         enddo
c..boundary condition at origin
         v2up(1) = v2down(1)
         v2up(2) = v2down(2)
         do k = 3, kk
c            v2up(k) = v2up(k-1) -bv2(k)*dellr(k)*hml(k)
            v2up(k) = ( v2up(k-1) -bv2(k)*dellr(k)*hml(k)*2.0d0)
     1           *hml(k)/(hml(k) + 16.0d0*dellr(k) )

            if( v2up(k) .lt. 0.0d0 )then
               v2up(k) = 0.0d0
            endif
         enddo


c..revise convective velocity h(k) for nonlocal integration
         if( nsweep .eq. 0 )then
c            write(*,'(10x,2a5,12a12)')'k','ic','h','veq','vdown',
c     1           'vup','nsqr','bv2','v2up','v2down','dnab','dnad',
c     2           'doux'
cccccccccccccccccccccccccccccccccccc
            do k = 1, kk
               fact = dsqrt( v2down(k) )
               fak  = dsqrt( v2up(k)   )
c               if( bv2(k) .lt. 0.0d0 )then
c                  veq = sqrt( -bv2(k)/8.0d0 )*hml(k)
c               else
c                  veq = 0
c               endif

c               write(*,'(a10,2i5,1p14e12.4)')'pay',k,ic(k),h(k),veq,
c     1              fact,fak,nsqr(k),bv2(k),v2up(k),v2down(k),dnab(k)-
c     2              dnad(k),doux(k), dnab(k)-dnad(k)-doux(k)
cccccccccccccccccccccccccc
               h(k) = 0.5d0*( fact + fak )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            enddo
         endif

c         write(*,'(5x,2a5,12a12)')'k','ic','h','veq',
c     1        'nsqr','bv2','v2up','v2down','nab-nad','doux','ab-ad-x',
c     2        'dely'
         do k = 2, kk
            if( bv2(k) .lt. 0.0d0 )then
               veq = sqrt( -bv2(k)/8.0d0 )*hml(k)
            else
               veq = 0
            endif
ccccccccccccccccccccccccccccccc4
            fak =  pi4*r(1,k)**4*0.5d0*( p(1,k)+p(1,k+1) )/
     1           ( grav*xm(k)*dmi(k) )
c            fact =  rgas*(x(nnuc-1,k)+x(nnuc,k)+x(nnuc+1,k))*t(1,k)/
c     6          (v(1,k)*p(1,k))
            ybar = 0.0d0
            do n = 1, nnuc+1
               ybar = ybar + x(n,k)
            enddo

            if( dellr(k)**2*nsqr(k) .gt. 0.0d0 )then
               fact = sqrt(  dellr(k)**2*nsqr(k) )
            else
               fact = 0.0d0
            endif

c            write(*,'(a5,2i5,1p14e12.4)')'pay',k,ic(k),h(k),veq,
c     1           nsqr(k),bv2(k),v2up(k),v2down(k),dnab(k)-
c     2           dnad(k),doux(k), dnab(k)-dnad(k)-doux(k),dely(k),
c     3           x(nnuc-1,k),fact
        enddo

c        write(*,'(5x,2a5,12a12)')'k','ic','h','veq',
c     1       'nsqr','bv2','v2up','v2down','nab-nad','doux','ab-ad-x',
c     2       'dely'
ccccccccccccccccccccccccccccccccccccc
c         write(*,'(1p8e12.3)') vvel(nvmax(3),3), vvel(nvmax(3),3)**2
c         write(*,'(1p8e12.3)') v2up(1),v2down(1),v2down(2)

c         stop'pay aaa'
ccccccccccccccccccccccc

         dellr(1) = dellr(2)
         hml(1)   = hml(2)

c================================================================
c..   PAY overshooting starts here..........................
         if( mixmode .eq. 0 .or. mixmode .eq. 1 )then
c..   force to value consistent with simulations
            m = 1
            reflec = 0
            do k = 2,  kk+2
               dif(k) = 0.0d0
            enddo
            do i=2,kk
               kh(i) = 0.0d0
            enddo
            omegc = 0.0d0
c..   up sweep
            do k = 2, kk
               astab(k) = astab(k)/hml(k)
 
c..   Detect edge of Kippenhan-Weigert convective region, using ic flag
               if( ic(k) .eq. 1 )then
c..   convective
                  m = k
                  reflec = k
                  sigmath(k) = 4.0d0*(gamma1(k)-1.0d0)/gamma1(k)*arad/
     1                 3.0d0*t(1,k)**4.0d0/p(1,k)*crad*v(1,k)/ak(k)
                  sigmamup(k) = sigmath(k)
                  kv(k) =1.0d0/( hml(k)*2.1d0)
                  omegc = ( h(k) * grav*xm(k)/r(1,k)**2 /
     1                 (0.5d0*( p(1,k)*v(1,k) + p(1,k+1)*v(1,k+1) )) )
c..   Calculate diffusion coefficient outside convective region based
c..   on Press(1981) 
               else
c..   radiative
                  dellr(k) = dabs(r(1,k) - r(1,m))
c..sqrt( 2*3 ) = sqrt(6)
c..sqrt( 7*8 ) = sqrt(56) = 7.48
                  kh(k) = 50.48d0/r(1,k)
c     kh(k) = (6.0d0/r(1,k)**2.0d0)**0.5d0
c..   ????????
c     if(x(3,nnuc-1) .lt. 1.0d-5)then
c..   order of indices is reversed; should be x(nucleus,zone)
c                  if(x(nnuc-1,3) .lt. 1.0d-5)then
c                     kh(k) = (12.0d0/r(1,k)**2.0d0)**0.5d0
c                  else
c                     kh(k) = (6.0d0/r(1,k)**2.0d0)**0.5d0
c                  endif
                  sigmath(k) = 4.0d0*(gamma1(k)-1.0d0)/gamma1(k)*arad/
     1                 3.0d0*t(1,k)**4.0d0/p(1,k)*crad*v(1,k)/ak(k)
                  qfak(k) = (2.0d0*nsqr(k)**0.5d0/(sigmath(k)*
     1                 kh(k)**2.0d0)*(omegc/nsqr(k)**0.5d0)**3.0d0*
     1                 (nsqr(k) - omegc**2.0d0)/nsqr(k))
                  if(nsqr(k) .le. 0.0d0)then
                     qfak(k) = 1.0d-14
                  endif
                  if((nsqr(k)/omegc**2.0d0 -1.0d0) .le. 0.0d0)then
                     qfak(k) = 1.0d0 
                  endif
                  kv(k) = kh(k)*(nsqr(k)/omegc**2.0d0 -1.0d0)**0.5d0
c                  write(*,'(a20,1p8e12.3)')'kv 1',kv(k),
c     1                 nsqr(k),kh(k),omegc
cccccccccccccccccccc

                  if(nsqr(k)-omegc**2.0d0 .gt. 0.0d0)then
                     reflec = k
                  else
                     kv(k)=kv(reflec-1)*dexp((r(1,k)-r(1,reflec-1))
     1                    *kv(reflec-1))
                  endif
c                  write(*,'(a20,1p8e12.3)')'kv 2',kv(k)
cccccccccccccccccccc
                  epssqr(k) = kv(k)**2.0d0*(omegc/kh(k))**2.0d0/
     1                 (nsqr(k)-omegc**2.0d0)
                  if((nsqr(k)/omegc**2.0d0 -1.0d0) .le. 0.0d0)then
                     epssqr(k) = dabs(epssqr(k))
                  endif
                  if(nsqr(k) .le. 0.0d0)then
                     epssqr(k) = 1.0d-14
                  endif
                  if(epssqr(k) .eq. 1.0d0)then
                     dmom(k) = 0.5d0
                  endif
                  hmlfak = dexp(-pi*(dellr(k)*kh(k)/qfak(k)))
                  sigmamup(k) = hmlfak*sigmath(k)
     1                 *epssqr(k)**2.0d0/qfak(k)
c                  write(*,'(a20,i5,1p8e12.3)')'sigmamup 1',k,sigmamup(k)
cccccccccccccccccccccccccc

                  if(nsqr(k)-omegc**2.0d0 .gt. 0.0d0)then
                     reflec = k
                  else
                     dumvar = kv(k)/ dabs(qfak(k))
                     sigmamup(k) = sigmamup(reflec-1)*dexp(-pi*
     1                    ((r(1,k)-r(1,reflec))*dumvar))
c                  write(*,'(a20,i5,1p8e12.3)')'sigmamup 2',k,sigmamup(k)
cccccccccccccccccccccccccc

                  endif
                  if(nsqr(k) .gt. 0.0d0)then
                     reflec = k
                  else
                     sigmamup(k) = sigmamup(reflec-1)*dexp(-pi*
     1                    ((r(1,k)-r(1,reflec))*kv(k)/
     1                    dabs(qfak(k))))
c                  write(*,'(a20,i5,1p8e12.3)')'sigmamup 3',k,sigmamup(k)
ccccccccccccccccccccccccccc

                  endif
                  if(omegc .eq. 0.0d0)then
                     sigmamup(k) = 0.0d0
c                  write(*,'(a20,i5,1p8e12.3)')'sigmamup 4',k,sigmamup(k)
cccccccccccccccccccccccccc

                  endif

c                  stop'pay'
cccccccccccc

               endif
            enddo

c..   down sweep
            omegc = 0.0d0
            mm = kk
            do j = 2, kk            
               k = kk +1 -j
               astab(k) = astab(k)/hml(k)
               
c..   Detect edge of Kippenhan-Weigert convective region, using
c..   ic flags
               if( ic(k) .eq. 1 )then
                  mm = k
                  reflec = k
                  sigmamdown(k) = sigmath(k)
                  omegc = ( h(k) * grav*xm(k)/r(1,k)**2 /
     1                 (0.5d0*( p(1,k)*v(1,k) + p(1,k+1)*v(1,k+1) )) )
c..   Calculate diffusion coefficient for internal wave mixing from 
c..   dissipation according to Press(1981)
               else
                  dellr(k) = dabs(r(1,k) - r(1,mm))
c..   pay 10/20/04
c                  kh(k) = (6.0d0/r(1,k)**2)**0.5d0
c..   sqrt(7*8) = 7.48
                  kh(k) = 50.48d0/r(1,k)

                  qfak(k) = (2.0d0*nsqr(k)**0.5d0/(sigmath(k)*
     1                 kh(k)**2.0d0)*(omegc/nsqr(k)**0.5d0)**3.0d0*
     1                 (nsqr(k) - omegc**2.0d0)/nsqr(k))
                  if(nsqr(k) .le. 0.0d0)then
                     qfak(k) = 1.0d-14
                  endif
                  if((nsqr(k)/omegc**2.0d0 -1.0d0) .le. 0.0d0)then
                     qfak(k) = 1.0d0 
                  endif
                  kv(k) = kh(k)*(nsqr(k)/omegc**2.0d0 -1.0d0)**0.5d0
                  if(nsqr(k)-omegc**2.0d0 .gt. 0.0d0)then
                     reflec = k
                  else
                     kv(k)=kv(reflec+1)*dexp(dabs(r(1,k)-r(1,reflec+1))
     1                    *kv(reflec+1))
                  endif
                  epssqr(k) = kv(k)**2.0d0*(omegc/kh(k))**2.0d0/
     1                 (nsqr(k)-omegc**2.0d0)
                  if((nsqr(k)/omegc**2.0d0 -1.0d0) .le. 0.0d0)then
                     epssqr(k) = dabs(epssqr(k))
                  endif
                  if(nsqr(k) .le. 0.0d0)then
                     epssqr(k) = 1.0d-14
                  endif
                  if(epssqr(k) .eq. 1.0d0)then
                     dmom(k) = 0.5d0
                  endif
                  hmlfak = dexp(-pi*(dellr(k)*kh(k)/qfak(k)))
                  sigmamdown(k) = hmlfak*sigmath(k)
     1                 *epssqr(k)**2.0d0/qfak(k)
                  if(nsqr(k)-omegc**2.0d0 .gt. 0.0d0)then
                     reflec = k
                  else
                     dumvar = kv(k)/ dabs(qfak(k))
                     sigmamdown(k) = sigmamdown(reflec+1)*dexp(-pi*
     1                    (dabs(r(1,k)-r(1,reflec))*dumvar))
                  endif
                  if(nsqr(k) .gt. 0.0d0)then
                     reflec = k
                  else
                     sigmamdown(k) = sigmamdown(reflec+1)*dexp(-pi*
     1                    (dabs(r(1,k)-r(1,reflec))*kv(k)/
     1                    dabs(qfak(k))))
                  endif
                  if(omegc .eq. 0.0d0)then
                     sigmamdown(k) = 0.0d0
                  endif
                  
               endif
            enddo

c..   Define diffusion coefficient using Press coefficient
            do k = 2, kk
               if(ic(k) .ne. 1)then
                  sigmam(k) = dmax1(sigmamup(k),sigmamdown(k))
                  dif(k) = sigmam(k)*dth(2)/(r(2,k) - r(2,k-1))**2.0d0
               else
c..   dif(kk)=0, protect against dellr(kk)=0
                  dif(k) = 0.0d0
               endif
            enddo
         else
            do k = 2, kk
               vover(k)  = 0.0d0
               vunder(k) = 0.0d0
               dif(k)    = 0.0d0
            enddo
         endif
c..   PAY overshoot ends here..................................
c=============================================================
         
      endif


ccccccccccccccccccccccccccccccccccccccccccccc
c      write(*,'(2a5,12a12)')'k','ic','h','bv2','nsqr','dif','Xhe4'
c      stop'pay'
cccccccccccccccccccccccccccccccccccccccccccccccc

      return
      end

