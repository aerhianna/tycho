      subroutine vlog(n1,n2,ymin,ymax,yy)

      implicit none

      include 'dimenfile'

      real*4 yy(kdm)
      real*4 ymin, ymax, ten, y10min, y10max

      integer*4 n1, n2, i
c------------------------------------------

      ten = 1.0e1
      y10min = ten**ymin
      y10max = ten**ymax

      do i = n1, n2

        if( yy(i) .gt. y10min .and. yy(i) .le. y10max )then
          yy(i) = alog10( yy(i) )
        elseif( yy(i) .le.  y10min )then
          yy(i) = ymin
        else
          yy(i) = ymax
        endif

      enddo

      return
      end


