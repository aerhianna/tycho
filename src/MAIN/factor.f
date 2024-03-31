      subroutine factor(a,w,ipivot,d,n,iflag)
      implicit none

c..gauss elimination, scaled partial pivoting, triangular factorization 
c         a  is input matrix
c         w is factored result
c         ipivot contains pivoting strategy
c         d is working space for keeping size of row elements
c     follow by subst to get solution of ax = b
c     iflag = 0 means ok
c             1 means error
c     algorithm from conte and de boor, p 132-133
c     *** they did not define w initially ***

      integer*4 n,iflag,i,j,nm1,k,kp1,ip,ipk
      integer*4 ipivot(n)
      real*8    a(n,n),w(n,n),d(n)
      real*8    rowmax, colmax, awikov, ratio
c----------------------------------------------------------

      iflag = 0
c..initialize w, ipivot, d
      do i = 1, n
        ipivot(i) = i
        rowmax    = 0.0d0

        do j = 1, n
          w(i,j) = a(i,j)
          rowmax = dmax1( rowmax, abs(w(i,j)) )
        enddo

        if( rowmax .eq. 0.0d0 )then
           write(*,*)'FACTOR.f MSG: rowmax = ',rowmax
           write(*,'(a5,1p8e12.4)')'a',a
           write(*,'(a5,1p8e12.4)')'w',w
           go to 999
        endif
        d(i) = rowmax
      enddo

c..gauss elim. with partial pivoting
      ip  = 1
      nm1 = n - 1
      if( nm1 .eq. 0 )return

      do k = 1, nm1
        j      = k
        kp1    = k + 1
        ip     = ipivot(k)
        colmax = dabs( w(ip,k) )/d(ip)

        do i = kp1, n
          ip     = ipivot(i)
          awikov = dabs( w(ip,k) )/d(ip)
          if( awikov .gt. colmax )then
            colmax = awikov
            j      = i
          endif
        enddo
        if ( colmax .eq. 0.0d0 )  go to 999
        ipk       = ipivot(j)
        ipivot(j) = ipivot(k)
        ipivot(k) = ipk
        do i = kp1, n
          ip      = ipivot(i)
          w(ip,k) = w(ip,k)/w(ipk,k)
          ratio   = - w(ip,k)
          do j = kp1, n
            w(ip,j) = ratio * w(ipk,j) + w(ip,j)
          enddo
        enddo
      enddo

      if( w(ip,n) .eq. 0.0d0 ) go to 999
      return

  999 continue
      iflag = 1
      return
      end

