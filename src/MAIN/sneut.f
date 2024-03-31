      subroutine sneut(k,t9,rho,aye,aenu)

c..beaudet, petrosian and salpeter neutrino approximations

      implicit none

      include 'dimenfile'

      integer k, n

c..external
      real*8 t9,rho
c..internal
      real*8 d,y,y2,y4,y6,y8,d93,e,e2,e3,
     1 ym1,ym2,ym3,ex,anum,aden,g,q

      real*8 aye(kdm),aenu(kdm)
      real*8 a(3,3),b(3,3),c(3),f(3)

      save a,b,c

      data a/
     1  6.002d+19,4.886d+10,2.320d-07,
     2  2.084d+20,7.580d+10,8.449d-08,
     3  1.872d+21,6.023d+10,1.787d-08/
      data b/
     1  9.383d-01,6.290d-03,2.581d-02,
     2 -4.141d-01,7.483d-03,1.734d-02,
     3  5.829d-02,3.061d-04,6.990d-04/
      data c/
     1  5.5924d+00,1.5654d+00,0.56457d+00/
c---------------------------------------------------------------------

      if( t9 .lt. 1.0d-02 )then
        aenu(k) = 0
      else
        d   = aye(k) * rho
        y   = t9 / 5.9302d0
        y2  = y*y
        y4  = y2*y2
        y6  = y4*y2
        y8  = y4*y4
        d93 = (d*1.0d-09)**0.33333333d0
        e   = d93/y
        e2  = e*e
        e3  = e2*e
        ym1 = 1.0d0/y
        ym2 = ym1**2
        ym3 = ym2*ym1
        do n = 1, 3
          ex   = dexp( -c(n)*e )
          anum = ( a(n,1) + a(n,2)*e + a(n,3)*e2 )*ex
          aden = e3 + b(n,1)*ym1 + b(n,2)*ym2 + b(n,3)*ym3
          f(n) = anum/aden
        enddo
        g = 1.0d0 - 13.04d0*y2 + 133.5d0*y4
     1            + 1534.0d0*y6 + 940.0d0*y8
        q = d**3*f(3) + d*y4*y*f(2) + g*dexp(-2.0d0*ym1)*f(1)
        aenu(k) = -q/rho
      endif
      return
      end







