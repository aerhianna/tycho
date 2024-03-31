      subroutine lintrx(qqo,qqn,x,scrx,jj,njj)
      implicit none

c..   wda 9/24/04
c..   for abundance arrays
c     interpolates for each  qqn(njj)   in old  qqo(jj) array
c     to find              new f(njj) from old    f(jj)
c     scr is scratch space used during interpolation
c     assumes that qqo and qqn are monotonic increasing functions of j

c..   intmod=0 does simple histogram assignment
c     will extrapolate below qqo(1) or above qqo(jj) if need be

c     interpolates in x
c     note that x is changed upon return

      include 'dimenfile'

      integer*4 intmod
      integer*4 jj,njj,ki,n,k,i
      real*8    qqo(kdm),qqn(kdm),x(ndim,kdm),scrx(ndim,kdm)
      real*8    fak1,fakp,fakm
      data intmod/1/
c--------------------------------------------------------------
c     ki = 2 forces extrapolation in inner zone
      ki = 2

      do i = ki, njj
c..   find zone neighbors
         do k = ki, jj-1
c..trigger on first time qqo(k+1) increases to .ge. qqn(i) 
            if(  qqo(k+1) .gt. qqn(i) )  go to 1000
         enddo
         k = jj-1
 1000    continue
         if( intmod .eq. 0 )then
c..   nearest neighbors found: qqo(k) .lt. qqn(i) .le. qqo(k+1)
c..   abundancees on zone centers x(n,k+1/2) = x(n,"k+1")
c..   use zone center to avoid logic error by roundoff
c..   histogram logic
            if( 0.5d0*(qqn(i) + qqn(i-1)) .lt. qqo(k) )then
               do n = 1, ndim
                  scrx(n,i) = x(n,k)
               enddo
            else
               do n = 1, ndim
                  scrx(n,i) = x(n,k+1)
               enddo
            endif

         else
c..   linear interpolation
            fak1   =   qqo(k+1) - qqo(k)
            fakp   = ( qqn(i)   - qqo(k) )/fak1
            fakm   = ( qqo(k+1) - qqn(i) )/fak1
c..   nearest neighbors found: qqo(k) .lt. qqn(i) .le. qqo(k+1)
c..   abundancees on zone centers x(n,k+1/2) = x(n,"k+1")
c..   use zone center to avoid logic error by roundoff
c..   histogram logic
            if( 0.5d0*(qqn(i) + qqn(i-1)) .lt. qqo(k) )then
               do n = 1, ndim
                  scrx(n,i) = x(n,k+1)*fakp + x(n,k)*fakm
               enddo
            else
               do n = 1, ndim
                  scrx(n,i) = x(n,k+1)*fakp + x(n,k)*fakm
               enddo
            endif
         endif
      enddo

c..interpolation complete, update original array x(n,i)
      do i = ki, njj
         do n = 1, ndim
            x(n,i) = scrx(n,i)
         enddo
      enddo

      return
      end

