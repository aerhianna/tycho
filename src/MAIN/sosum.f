      subroutine sosum(kk,s,dmh,xm,dth,e,u)

      implicit none

      include 'dimenfile'
      include 'cgtr.h'

c..   sums sources and sinks over mass, over time evolution
      real*8 s(7,kdm),dmh(kdm),ss(7),xm(kdm)
      real*8 dth(2),e(2,kdm),u(2,kdm),emscl,dele,delu

      integer*4 kk, n, k

c-------------------------------------------------------------------
      emscl = xm(kk) - xm(1)
      do n = 1, 7
         ss(n) = 0.0d0
         do k = 2, kk
            ss(n) = ss(n) + s(n,k)*(dmh(k)/emscl)*dth(1)*efi(k)
         enddo
      enddo

      dele = 0
      delu = 0
      do k = 2, kk
         dele = dele + (e(2,k) - e(1,k))*(dmh(k)/emscl)
         delu = delu
     1        + (u(2,k)**2 - u(1,k)**2)*(0.5d0*dmh(k)/emscl)
      enddo

      write(3,*)'Breakdown of net source terms (erg/g)'
      write(3,'(9a12)')'rad','conv',' ','neutrino','nuclear',
     1     'total',' ','del E','del Kine'
      write(3,'(1p9e12.4)')( ss(n), n=1,7 ),dele,delu

      return
      end
