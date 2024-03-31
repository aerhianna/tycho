      subroutine minmax(yy,k1,kk,vmin1,vmax1)

c finds maximum and minimum of real*4 array,
c for index range: k1 .le. k .le. kk
      implicit none

      real*4 yy( 1 )
      real*4 vmin1, vmax1
      real*4 vmin, vmax, del, offset

      integer*4 k1, kk
      integer*4 k2, k
c-------------------------------------------
      vmin = yy(k1)
      vmax = yy(k1)
      k2 = k1 + 1

      do  k = k2, kk
        vmin = amin1( yy(k), vmin )
        vmax = amax1( yy(k), vmax )
      enddo

c..add offset
      del    = vmax - vmin
      offset = 0.05*abs(del)
      vmin1  = vmin - offset
      vmax1  = vmax + offset

      return
      end



