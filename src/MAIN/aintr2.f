      subroutine aintr2(qqo,qqn,f,scr,jj,njj,m)
      implicit none

c..   for double arrays
c     interpolates for each  qqn(njj)   in old  qqo(jj) array
c     to find              new f(m,njj) from old    f(m,jj)
c     scr is scratch space used during interpolation
c     assumes that qqo and qqn are monotonic increasing functions of j
c     does simple quadratic interpolation
c     will extrapolate below qqo(1) or above qqo(jj) if need be

c     note that f is changed upon return

      include 'dimenfile'

      integer*4 jj,njj,m,n,k,ki
      real*8    qqo(kdm),qqn(kdm),f(2,kdm),scr(kdm)
      real*8    fact1,fact2,fact3
c--------------------------------------------------------------

      ki = 1

      do n = ki, njj
c..   find zone neighbors
         do k = ki, jj
            if(       qqo(k) .le. qqn(n)
     1           .and. qqn(n) .lt. qqo(k+1) )then
               go to 1000
            endif
         enddo
 1000    continue
         if( qqn(n) .gt. 0.5d0*(qqo(k) + qqo(k+1))   )then
            k = k+1
         endif
         if( k .gt. jj-2 )then
            k = jj-2
         endif
c..   quadratic interpolation
         fact1 = ( (qqn(n)   - qqo(k+1))*(qqn(n)   - qqo(k+2)) )/
     1        ( (qqo(k)   - qqo(k+1))*(qqo(k)   - qqo(k+2)) )
         fact2 = ( (qqn(n)   - qqo(k)  )*(qqn(n)   - qqo(k+2)) )/
     1        ( (qqo(k+1) - qqo(k)  )*(qqo(k+1) - qqo(k+2)) )
         fact3 = ( (qqn(n)   - qqo(k+1))*(qqn(n)   - qqo(k)  ) )/
     1        ( (qqo(k+2) - qqo(k+1))*(qqo(k+2) - qqo(k)  ) )
         scr(n) = fact1*f(m,k) + fact2*f(m,k+1) + fact3*f(m,k+2)
      enddo


      do n = ki, njj
         f(m,n) = scr(n)
      enddo

      return
      end



