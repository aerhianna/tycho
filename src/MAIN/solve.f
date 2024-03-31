      subroutine solve( a, b, c, q,qlo,qln,qko,qkn,jj,jsq,k)
      implicit none

      integer*4 jj, jsq, k, i,iflg
      real*8 a(jsq),b(jsq),c(jsq),q(jj),qlo(jsq),qln(jsq)
      real*8 qko(jj),qkn(jj)
      real*8 cmbl(9),wrk(9),qpbk(3),d(3)
      integer*4 ipivot(3)
c..solves henyey matrix equation
c.............................................................
c     uses henyey conventions

c     q + a*w(+) + b*w(-) + c*w = 0
c     w = k - l*w(+)

c     jsq .le. jyey*jyey
c     jj  .le. jyey
c............................................................

c..set up c - b*qlo array................


      call mult(b,qlo,cmbl,jj,jj,jj)

      do i = 1,jsq
         cmbl(i) = c(i) - cmbl(i)
      enddo
c..get qln array.........................

      call divide(cmbl,a,qln,wrk,d,ipivot,jj,jsq,jsq,iflg,k)

c..get qpbk vector.......................

      call mult(b,qko,qpbk,jj,1,jj)

      do i = 1,jj
         qpbk(i) = - (qpbk(i) + q(i))
      enddo

c..get qkn vector........................

      call divide(cmbl,qpbk,qkn,wrk,d,ipivot,jj,jsq,jj,iflg,k)

      return
      end



