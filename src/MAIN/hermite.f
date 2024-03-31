c     
c     Hermite Interpolation
c     
      subroutine hermite( f0her,f1her,d0her,d1her,dmher,fhher)
c..   hermite interpolation
      implicit none
      real*8 f0her,f1her,d0her,d1her,dmher,fhher
      real*8 a,b,c,d,m,f,h
c.....
      a = f0her
      b = d0her
      m = dmher
      d = ( d1her + d0her + (f1her-f0her)/m  )/m**2
      c = 3.0d0*( (f1her - f0her)/m**2 - d1her/m )
      h = m/2.0d0
      f = a + b*h + c*h**2 + d*h**3
      
      fhher = f
      
      return
      
      end
      
