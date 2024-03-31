      subroutine locate(xx,n,x,j)
      implicit none
c..   numerical recipes bisection algorithm
c..   locate index j for variable x lying between j, j+1
c..   in array xx of dimension n
c..   j  = desired index
c..   xx = array of tabulated values
c..   n  = dimension of xx
c..   x  = variable value to be interpolated

      integer*4 j,n
      real*8    x, xx(n)
      integer*4 jl,jm,ju
c---------------------------------------------------------------------
      jl = 0
      ju = n+1
 10   if( ju-jl .gt. 1 )then
         jm = (ju+jl)/2
         if( (xx(n) .gt. xx(1)) .eqv. (x .gt. xx(jm)) )then
            jl = jm
         else
            ju = jm
         endif
         goto 10
      endif

      j = jl

      return
      end
