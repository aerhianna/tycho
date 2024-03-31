

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..   example of polint use
c     i = k - 1
c     i = max0(i,1)
c     call polint(qq2o(i),scr1(i),3,qq2n(k),scr(k),scr2(k))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      

      SUBROUTINE polint(xa,ya,n,x,y,dy)
c..   numerical recipes subroutine for polynomial interpolation
      INTEGER*4 n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER*4 i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

      ns=1
      dif=dabs(x-xa(1))
      do i=1,n
         dift=dabs(x-xa(i))
         if (dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
         do i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0.D0)stop 'failure in polint'
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         enddo
         if (2*ns.lt.n-m)then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
      enddo


c     SUCCESS
      return
c     
      END
      
