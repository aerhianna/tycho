      subroutine divide(a,b,c,w,d,ipivot,n,mn,mb,iflag,k)
      implicit none

      integer*4 n,m,mn,mb,iflag,ibeg,i,k

      integer*4 ipivot(n)

      real*8 a(mn),b(mb),c(mb),w(mn),d(n)

c  matrix division ac = b
c    or c = inverse of a * b
c  mn = dimension of square matrix a
c  mb = no. of columns in b, c

c     note that a(i,j) is
c     1,1     1,2
c     2,1     2,2
c     and is also a(k) where
c     1     3
c     2     4
c-----------------------------------------------------------
      m = mb/n

            call factor(a,w,ipivot,d,n,iflag)

      if( iflag .ne. 0 )then
         write(*,'(10a10)')'n','ipivot(n)','iflag','k'
         write(*,'(10i10)')n,ipivot(n),iflag,k
         stop' in divide.f'
      endif
      ibeg = 1

      do i = 1, m
            call subst(w,b(ibeg),c(ibeg),ipivot,n)
        ibeg = ibeg + n
      enddo

      return
      end



