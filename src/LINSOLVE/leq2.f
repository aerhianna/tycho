      subroutine leq2(aa,bb,n,ndim)
      implicit none

c..   external
      integer*4 ndim,n
      real*8 aa(ndim,ndim),bb(ndim)

c..   internal
      real*8 x(ndim), a,b,c,d,e,f,disc
c-----------------------------------------------------
c.....leq2 performs matrix inversion
c     2-dimension of a
      a = aa(1,1)
      b = aa(1,2)
      c = bb(1)
      d = aa(2,1)
      e = aa(2,2)
      f = bb(2)

      disc = a*e - b*d

      x(1) = (c*e - b*f)/disc
      x(2) = (a*f - c*d)/disc
      bb(1) = x(1)
      bb(2) = x(2)

      return
      end

