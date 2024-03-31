      subroutine stabil(nc,gamma,del)
      implicit none

      include 'dimenfile'

c      integer*4 k, kk, nc
c      real*8    p(2,kdm),v(2,kdm),t(2,kdm),e(2,kdm)
c      real*8    dmh(kdm),pv(kdm),pt(kdm),ev(kdm),et(kdm)
c      real*8    xm(kdm)
      include 'comod'

      integer*4 k, nc
      real*8    scale,sum1,sum2,dm,c,gamma,del
      real*8    dtdv,gama
c---------------------------------------------------------
c..   calculates gamma1 = d ln P / d ln rho at constant entropy,
c..   weighted by P*V * dm, over whole star, for stability test
      sum1  = 0.0d0
      sum2  = 0.0d0
      scale = 1.0d0/( xm(kk) - xm(1) )
      do k = 2, kk
         dm   = dmh(k)*scale
         c    = p(nc,k)*v(nc,k)*dm
c..   construct gamma
         dtdv = - (p(nc,k) + ev(k))/et(k)
         gama  = - v(nc,k)/p(nc,k)*( pv(k) + pt(k)*dtdv )
         sum1 = sum1 +      c
         sum2 = sum2 + gama*c
      enddo
      gamma = sum2/sum1
      del   = gamma*0.75d0 - 1.0d0
      return
      end
