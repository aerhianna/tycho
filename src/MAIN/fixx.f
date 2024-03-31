c...................................................................................
      
      subroutine fixx(k)

      implicit none

      include 'dimenfile'

      include 'cconst'
      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'comdif'
      include 'comfixx'

      integer*4 k,n,j,nc
      real*8 dfact,fact
      
c.....
      
      fakmu = dmh(k)*dmh(k+1)/(dmh(k)+dmh(k+1))
      
      print*,'MSG(cmix::fixx): fixx(k) CALLED FROM cmix: k,xm=',
     .     k,xm(k)/sol
      
      do n = 1, nnuc+1
c..   fakx is mix parameter between zones k and k+1
         x(n,k) = xold(n,k)*(1.0d0 -fakx*fakmu/dmh(k))
     1        +xold(n,k+1)*fakx*fakmu/dmh(k)
         x(n,k+1) = xold(n,k)*fakx*fakmu/dmh(k+1)
     1        +xold(n,k+1)*(1.0d0-fakx*fakmu/dmh(k+1))
      enddo
      
      nc   = 2
      do j = 1, 20
         call state(k,k,nc)
         
         fact = p(2,k)-p(1,k)
         dfact = -fact/pt(k)
         t(2,k) = t(2,k) + dfact
         if( abs(fact) .lt. 1.0d-9*p(1,k) )goto 110
      enddo
      write(*,*)'cmix pressure iteration error ',j,k,l
      stop'cmix: fixx 1'
      
 110  continue
      do j = 1, 20
         call state(k+1,k+1,nc)
         
         fact = p(2,k+1)-p(1,k+1)
         dfact = -fact/pt(k+1)
         t(2,k+1) = t(2,k+1) + dfact
         if( abs(fact) .lt. 1.0d-9*p(1,k+1) )goto 111
      enddo
      write(*,*)'cmix pressure iteration error ',j,k+1
      stop'cmix fixx 2'

 111  continue

      ytot(k)   = 0
      ytot(k+1) = 0
      do n = 1,nnuc+1
         ytot(k)   = ytot(k)   + x(n,k)
         ytot(k+1) = ytot(k+1) + x(n,k+1)
      enddo

c..   reset nablas and convection indices
      call cinit(k,k+1,2)
      
      return
      end
