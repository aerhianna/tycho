      subroutine lintr1(qqo,qqn,f,scr,jj,njj)
      implicit none

c..   for single arrays
c     interpolates for each  qqn(njj)   in old  qqo(jj) array
c     to find              new f(njj) from old    f(jj)
c     scr is scratch space used during interpolation
c     assumes that qqo and qqn are monotonic increasing functions of j
c     does simple linear interpolation
c     will extrapolate below qqo(1) or above qqo(jj) if need be

c     note that f is changed upon return

      include 'dimenfile'

      integer*4 jj,njj,ki,n,k
      real*8    qqo(kdm),qqn(kdm),f(kdm),scr(kdm)
      real*8    fak1,fakp,fakm
c--------------------------------------------------------------
c     ki = 2 forces extrapolation in inner zone

      ki = 2

      do n = ki, njj
c..   find zone neighbors
         do k = ki, jj-1
            if(  qqn(n) .le. qqo(k+1) )  go to 1000
         enddo
         k = jj-1
 1000    continue
c..   linear interpolation
         fak1   =   qqo(k+1) - qqo(k)
         fakp   = ( qqn(n)   - qqo(k) )/fak1
         fakm   = ( qqo(k+1) - qqn(n) )/fak1
         scr(n) = f(k+1)*fakp + f(k)*fakm
      enddo

      do n = ki, njj
         f(n) = scr(n)
      enddo

      return
      end


