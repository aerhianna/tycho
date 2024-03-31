      function fion(f,y)
      implicit none
      real*8 fion, f, y, x, fact
c..ionization function------------------------
c     solves z**2 + x*z - x = 0
c     where fion = z and x = f/y
c           input : f, y
c           output : fion = fraction ionized
c---------------------------------------------
      if( y .gt. 1.0d-12 )then
        x = f/y
        if( x .le. 1.0d3 )then
          fact = x*x + 4.0d0*x
          fion = 0.5d0*( dsqrt(fact) - x )
        else
c..avoid round-off for large x
          fion = 1.0d0 - 1.0d0/x
        endif
      else
c..protects against zero divide
        fion = 1.0d0
      endif
      return
      end

