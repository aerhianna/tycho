      subroutine triad( a, b, c, d, u, j)
      implicit none

c..   tridiagonal solver

      include 'dimenfile'

      real*8 a(kdm),b(kdm),c(kdm),d(kdm),u(kdm),
     1     x(kdm),y(kdm)

      integer*4 j,k,m
      real*8 fact
c--------------------------------------------------
c     u(1)=t(2) in wda scheme
c     u(j)=t(kk)
c     aj*uj+1 + bj*uj + cj*uj-1 = dj

      x(j-1) = -c(j)/b(j)
      y(j-1) =  d(j)/b(j)

      do k = 2, j-1
         m    = j - k + 1
         fact = a(m)*x(m) + b(m)
         x(j-k) =  -c(m)              /fact
         y(j-k) = ( d(m) - a(m)*y(m) )/fact
      enddo

c      u(1) = ( d(1)-a(1)*y(1) )/( b(1) + a(1)*x(1) )
      u(1) = ( d(1)-a(1)*y(1) )/( b(1) + a(1)*x(1) + c(1) )

c..   general boundary condition
c..   u(1) = ( d(1)-a(1)*y(1) )/( b(1) + a(1)*x(1) + c(1) )
c..   but using spherical symmetry gives
c..   u(1) = -y(1)/( x(1) - 1.0d0 )

      do k = 1, j-1
         u(k+1) = x(k)*u(k) + y(k)
      enddo

      return
      end
