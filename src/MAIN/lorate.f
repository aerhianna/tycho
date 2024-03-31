      subroutine lorate(t9,rho,ye)

c..   fkt rate formula with T9, V derivatives (wda 5/9/99)

      implicit none

      include 'dimenfile'
      include 'cdeuter'
      include 'crate'

      real*8 t9, rho, t913, t953, t9ln, expon
      real*8 t92, t923, t943, factor, third
      real*8 ye

      parameter( third = 1.0d0/3.0d0 )

      integer*4 k, k1, k2

c----------------------------------------------------

      t913 = t9**third
      t953 = t9 * t913 *t913
      t9ln = dlog(t9)

      t92  = t9*t9
      t923 = t913*t913
      t943 = t913*t9

      do k = 1,lnreac
         expon = locoef(1,k) + locoef(2,k)/t9 + locoef(3,k)/t913
     1        + locoef(4,k)*t913 + locoef(5,k)*t9
     2        + locoef(6,k)*t953 + locoef(7,k)*t9ln

cccccccccccccccccccccccccccccccccccccccccccccccc
c     this is not general: casey 2/9/03
c         if( expon .gt. -300.0d0 .and. k .ne. 192 )then
c..he4 be9 n c12 is screwy (192)
c     fixed wda 2/22/03
ccccccccccccccccccccccccccccccccccccccccccccc

         if( expon .gt. -300.0d0 )then
            sig(k) = exp(expon)
            factor = -locoef(2,k)/t92 -third*locoef(3,k)/t943
     1           +third*locoef(4,k)/t923 + locoef(5,k)
     2           +5.0d0*third*locoef(6,k)*t923 + locoef(7,k)/t9

            sigt(k) = sig(k) * factor
         else
            sig(k)  = 0.0d0
            sigt(k) = 0.0d0
         endif
      enddo

c..sanity check at low temperature; some fkt rates blow up
      do k = 1, lnreac
         if( sig(k) .gt. 1.0d5 )then
            sig(k) = 1.0d5
            sigt(k) = 0.0d0
         endif
      enddo

c..   electron capture density*Ye factor
      k1 = l1deck(1)
      k2 = l2deck(2)
      do k = k1, k2
         if( lorlkh(k) .eq. ' ffn')then
c            logrho = dlog(rho)
c            call weakread(t9,logrho,locoef(1,k),expon)
c            write(*,*)rcoef(1,k),expon
            expon = locoef(1,k)
            if( expon .gt. -300.0d0 )then
               sig(k) = exp(expon)
               sigt(k) = 0.0d0
            endif
c            else
c               sig(k)  = 0.0d0
c               sigt(k) = 0.0d0
c            endif
c            if(t9 .le. 1.0d-2)then
c               sig(k)  = exp(locoef(2,k))
c               sigt(k) = 0.0d0
c            endif
c            if(logrho .le. 1.0d0)then
c               sig(k)  = exp(locoef(2,k))
c               sigt(k) = 0.0d0
c            endif
            if(locoef(1,k) .eq. 0.0d0)then
               sig(k)=0.0d0
               sigt(k)=0.0d0
            endif
c            write(*,*)rcoef(1,k),expon,sig(k)
c            sig(k)  = ye * rho * sig(k)
c            sigt(k) = ye * rho * sigt(k)
c            sigv(k) = -sig(k) * rho
            if( sig(k) .gt. 1.0d5 )then
               write(*,*)sig(k),ye,rho,locoef(1,k)
            endif
c            sig(k) = 0.0d0
c            sigt(k) = 0.0d0
         endif
      enddo
      
      do k = k1, k2
c         if( lorlkh(k) .eq. '  ec' .or. lorlkh(k) .eq. ' bec')then
         if(loec(k) .eq. 1)then
            sig(k)  = ye * rho * sig(k)
            sigt(k) = ye * rho * sig(k)
            sigv(k) = -sig(k) * rho
         endif
      enddo

c..   decks 4 to 7 have binary interaction in entrance channel
      k1 = l1deck(4)
      k2 = l2deck(7)
      do k = k1, k2
         if( lonrr(1,k) .eq. lonrr(2,k) )then
c..   pairs reduced for identical particles to avoid double counting
c..   sigv is constructed from sig
            sig(k)  = 0.5d0 * sig(k)
            sigt(k) = 0.5d0 * sigt(k)
         endif
c..   this is another power of rho
         sigv(k)  = -rho * rho * sig(k)
         sig(k)   =  rho * sig(k)
         sigt(k)  =  rho * sigt(k)
      enddo

c..   decks 8 and 9 have triple (tertiary collision) in entrance channel
      k1 = l1deck(8)
      k2 = l2deck(9)
      do k = k1, k2
         if( lonrr(1,k) .eq. lonrr(2,k) .and. 
     1       lonrr(2,k) .eq. lonrr(3,k))then
c..   waf, grc, baz (1967) ann rev a & a 5, 558 definition
c     has an extra factor of 1/2 for 3-alpha
c..   from 3 He4 destroyed per reaction, n/3! triples --> 1/2
c..   d/dV gives a compensating factor of 2
c..   factor of 1/3 moved from rhside and jacob to here 8/20/00
            sigv(k)  = -rho * rho * rho * sig(k) /3.0d0
            sig(k)   =        rho * rho * sig(k) /6.0d0
            sigt(k)  =        rho * rho * sigt(k)/6.0d0
         elseif(lonrr(1,k) .eq. lonrr(2,k) .or. 
     1          lonrr(1,k) .eq. lonrr(3,k)
     1          .or. lonrr(2,k) .eq. lonrr(3,k))then
            sigv(k)  = -rho * rho * rho * sig(k) 
            sig(k)   =        rho * rho * sig(k) /2.0d0
            sigt(k)  =        rho * rho * sigt(k)/2.0d0
         else
            sigv(k) = -rho * rho * rho * sig(k) * 2.0d0
            sig(k)  =        rho * rho * sig(k)
            sigt(k) =        rho * rho * sigt(k)
         endif
      enddo

c..   deck 10 has quadruple (quaternary collision) in entrance channel
      k1 = l1deck(10)
      k2 = l2deck(10)
      do k = k1, k2
         if( lonrr(1,k) .eq. lonrr(2,k) .and. 
     1       lonrr(2,k) .eq. lonrr(3,k)
     1       .and. lonrr(3,k) .eq. lonrr(4,k))then
c..   waf, grc, baz (1967) ann rev a & a 5, 558 definition
c     has an extra factor of 1/2 for 3-alpha
c..   from 3 He4 destroyed per reaction, n/3! triples --> 1/2
c..   d/dV gives a compensating factor of 2
c..   factor of 1/3 moved from rhside and jacob to here 8/20/00
            sigv(k)  = -rho * rho * rho * rho * sig(k) /8.0d0
            sig(k)   =        rho * rho * rho * sig(k) /24.0d0
            sigt(k)  =        rho * rho * rho * sigt(k)/24.0d0
         elseif(lonrr(1,k) .eq. lonrr(2,k) .and. 
     1          lonrr(2,k) .eq. lonrr(3,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /3.0d0
            sig(k)  =        rho * rho * rho * sig(k) /6.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/6.0d0
         elseif(lonrr(1,k) .eq. lonrr(2,k) .and. 
     1          lonrr(2,k) .eq. lonrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /3.0d0
            sig(k)  =        rho * rho * rho * sig(k) /6.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/6.0d0
         elseif(lonrr(2,k) .eq. lonrr(3,k) .and. 
     1          lonrr(3,k) .eq. lonrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /3.0d0
            sig(k)  =        rho * rho * rho * sig(k) /6.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/6.0d0
         elseif(lonrr(1,k) .eq. lonrr(3,k) .and. 
     1          lonrr(3,k) .eq. lonrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /3.0d0
            sig(k)  =        rho * rho * rho * sig(k) /6.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/6.0d0
         elseif(lonrr(1,k) .eq. lonrr(2,k) .and. 
     1          lonrr(3,k) .eq. lonrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /2.0d0
            sig(k)  =        rho * rho * rho * sig(k) /4.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/4.0d0
         elseif(lonrr(1,k) .eq. lonrr(2,k) .or. 
     1          lonrr(1,k) .eq. lonrr(3,k)
     1          .or. lonrr(1,k) .eq. lonrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) 
            sig(k)  =        rho * rho * rho * sig(k) /2.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/2.0d0
         elseif(lonrr(2,k) .eq. lonrr(3,k) .or. 
     1          lonrr(2,k) .eq. lonrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) 
            sig(k)  =        rho * rho * rho * sig(k) /2.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/2.0d0
         elseif(lonrr(3,k) .eq. lonrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) 
            sig(k)  =        rho * rho * rho * sig(k) /2.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/2.0d0
         else
            sigv(k) = -rho * rho * rho * rho * sig(k) * 2.0d0
            sig(k)  =        rho * rho * rho * sig(k)
            sigt(k) =        rho * rho * rho * sigt(k)
         endif
      enddo
      
c..   deck 11 has no density dependence
      k1 = l1deck(11)
      k2 = l2deck(11)
      do k = k1, k2
         sigv(k) = 0.0d0
      enddo

c..   deck 8 has triple in entrance channel
c      k1 = l1deck(8)
c      k2 = l2deck(8)
c      do k = k1, k2
c         if( lonrr(5,k) .eq. 0 )then
c..   waf, grc, baz (1967) ann rev a & a 5, 558 definition
c     has an extra factor of 1/2 for 3-alpha
c..   from 3 He4 destroyed per reaction, n/3! triples --> 1/2
c..   d/dV gives a compensating factor of 2
c..   factor of 1/3 moded from rhside and jacob to here 8/20/00
c            sigv(k)  = -rho * rho * rho * sig(k)/3.0d0
c            sig(k)   =  rho * rho * sig(k)  /6.0d0
c            sigt(k)  =  rho * rho * sigt(k) /6.0d0
c         else
c            sigv(k) = -rho * rho * rho * sig(k) * 2.0d0
c            sig(k)  =        rho * rho * sig(k)
c            sigt(k) =        rho * rho * sigt(k)
c         endif
c      enddo


      return
      end








