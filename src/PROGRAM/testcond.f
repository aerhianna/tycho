      program testcond

      implicit none

c..   driver for testing conduction opacity
c..   A. V. Sweigart, 1973, A and Ap 42, 459
c..   fit to Hubbard, W. B. and Lampe, M. 1969, ApJS 18, 297
c..   and Cannuto, V., 1970, ApJ 159, 641

c..   wda fit is similar to Sweigart, to 50 percent or better

      real*8 akch, xh, zmet, rho, tem
      real*8 logakch, logrho, logtem
      real*8 logakcr, akcr, xx, ww, lgd1, lgd2, akc
      real*8 dlogtem
      real*8 ye,xhe,x2,akc2

      integer*4 n

      data lgd1/5.0d0/, lgd2/7.0d0/
      data dlogtem/0.5d0/
c---------------------------------------------------------
      do n = 1,5

c         if( n .eq. 1 )then
            zmet   = 0.02d0
            xh     = 0.70d0
c            rho = 1.0d4
            rho = 10.0d0**(4.265d0+dble(n-1)*dlogtem)
c            tem = 1.0d1**(7.0d0+dble(n-1)*dlogtem)
            tem = 3.891d7
c         elseif( n .eq. 2 )then
c            zmet = 0.02d0
c            xh   = 0.0d0
c            rho = 1.0d6
c            tem = 2.0d8
c         else
c            zmet = 1.0d0
c            xh = 0.0d0
c            rho = 1.0d8
c            tem = 1.0d9
c         endif

         logrho = dlog10( rho )
         logtem = dlog10( tem )

c..   fit to Hubbard and Lampe
         logakch = dlog10( 1.0d0 + zmet - 0.6d0*xh) - 14.6196d0
     1        -( 3.5853d0 + 0.1386d0*logrho)* logrho
     2        +( 5.1324d0 - 0.3219d0*logtem)* logtem
     3        +  0.3901d0*logrho*logtem
         akch = 10.0d0**logakch

c..   fit to Cannuto
c..   does not apply to H-rich mixtures
         logakcr = dlog10( 1.0d0 + zmet ) - 14.04d0 - 1.02d0*logrho
     1        +2.275d0*logtem
         akcr = 10.0d0**logakcr

         xx = 0.0d0

         if( logrho .le. lgd1 )then
            ww = 0.0d0
         elseif( logrho .ge. lgd2 )then
            ww = 1.0d0
         else
            xx = ( logrho - lgd1)/(lgd2 - lgd1)
            ww = xx**3 * (3.0d0 - 2.0d0*xx)
         endif

         akc = (1.0d0-ww)*akch + ww*akcr
         xhe = 1.0d0 - zmet - xh

c..   wda approximation in state.f
         ye = 0.5d0*( 1.0d0 + xh )
         x2      = ( 1.015d-6 * rho * ye )**0.666667d0
         akc2     = 1.77d-6*( 1.0d0 + x2 )*( tem / rho )**2
     1        * ( 0.4d0*xhe + xh*0.1d0 + 0.8d0*zmet )
c..   last factor adjusted from hubbard-lampe (1969, apjs 18,297)
c..   tables for H, He, C

         write(*,'(1p8e12.3)')rho,tem,akch,akcr,xx,ww,akc,akc2

      enddo

      end
