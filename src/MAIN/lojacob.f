      subroutine lojacob(sig,ylo,a,
     1      nrr,k1deck,k2deck,loz,idebug)

c..   jacobian: d( dYi/dt vector )/dYj = (a)ij
c..   specify array elements of matrix a for ax = b

      implicit none

      include 'dimenfile'
c      include 'crate'
c      include 'cdeuter'

      real*8 sig(nreac), ylo(lodim), a(ndim,ndim)
      integer*4 nrr(6,nreac),k1deck(8),k2deck(8),loz(lodim)
      real*8 weaksum

      integer*4 i, j, k,idebug
      integer*4 k1,k2, i1,i2,i3,i4,i5,i6
      real*8    sigk,y2sig,y1sig
      real*8    y12sig,y23sig,y13sig
c---------------------------------------------------------------
      do i=1,ndim
         do j=1,ndim
            a(j,i) = 0.0d0
         enddo
      enddo

c..   beta decays, positron decays, electron captures
      weaksum = 0.0d0
      k1 = k1deck(1)
      k2 = k2deck(1)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            a(i1,i1) = a(i1,i1) + sig(k)
            a(i2,i1) = a(i2,i1) - sig(k)
            weaksum = weaksum + ( loz(i2)-loz(i1) )*ylo(i1)*sig(k)
c..   dYe/dt is weaksum
         enddo
      endif

c----------------------------------------------------
c..   deck 2: dissociations, i ---> j + k
      k1 = k1deck(2)
      k2 = k2deck(2)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            sigk =       sig(k)
            a(i1,i1) = a(i1,i1) + sigk
            a(i2,i1) = a(i2,i1) - sigk
            a(i3,i1) = a(i3,i1) - sigk
         enddo
      endif
c----------------------------------------------------
c..   deck 3: dissociations, i ---> j + k + l
      k1 = k1deck(3)
      k2 = k2deck(3)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            sigk =       sig(k)
            a(i1,i1) = a(i1,i1) + sigk
            a(i2,i1) = a(i2,i1) - sigk
            a(i3,i1) = a(i3,i1) - sigk
            a(i4,i1) = a(i4,i1) - sigk
         enddo
      endif
c----------------------------------------------------
c..   deck 4: captures, i + j ---> k
      k1 = k1deck(4)
      k2 = k2deck(4)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            y2sig =       ylo(i2)*sig(k)
            y1sig = ylo(i1)*      sig(k)
            a(i1,i1) = a(i1,i1) + y2sig
            a(i1,i2) = a(i1,i2) + y1sig
            a(i2,i2) = a(i2,i2) + y1sig
            a(i2,i1) = a(i2,i1) + y2sig
            a(i3,i1) = a(i3,i1) - y2sig
            a(i3,i2) = a(i3,i2) - y1sig
         enddo
      endif
c----------------------------------------------------
c..   deck 5: exchange, i + j ---> k + l
      k1 = k1deck(5)
      k2 = k2deck(5)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            y2sig =       ylo(i2)*sig(k)
            y1sig = ylo(i1)*      sig(k)
            a(i1,i1) = a(i1,i1) + y2sig
            a(i1,i2) = a(i1,i2) + y1sig
            a(i2,i2) = a(i2,i2) + y1sig
            a(i2,i1) = a(i2,i1) + y2sig
            a(i3,i1) = a(i3,i1) - y2sig
            a(i3,i2) = a(i3,i2) - y1sig
            a(i4,i1) = a(i4,i1) - y2sig
            a(i4,i2) = a(i4,i2) - y1sig
         enddo
      endif
c--------------------------------------------------------------
c..   deck 6: exchange, i + j ---> k + l + m
      k1 = k1deck(6)
      k2 = k2deck(6)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            y2sig =       ylo(i2)*sig(k)
            y1sig = ylo(i1)*      sig(k)
            a(i1,i1) = a(i1,i1) + y2sig
            a(i1,i2) = a(i1,i2) + y1sig
            a(i2,i2) = a(i2,i2) + y1sig
            a(i2,i1) = a(i2,i1) + y2sig
            a(i3,i1) = a(i3,i1) - y2sig
            a(i3,i2) = a(i3,i2) - y1sig
            a(i4,i1) = a(i4,i1) - y2sig
            a(i4,i2) = a(i4,i2) - y1sig
            a(i5,i1) = a(i5,i1) - y2sig
            a(i5,i2) = a(i5,i2) - y1sig
         enddo
      endif
c-------------------------------------------------------------------
c..   deck 7: exchange, i + j ---> k + l + m + n
      k1 = k1deck(7)
      k2 = k2deck(7)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            y2sig =       ylo(i2)*sig(k)
            y1sig = ylo(i1)*      sig(k)
            a(i1,i1) = a(i1,i1) + y2sig
            a(i1,i2) = a(i1,i2) + y1sig
            a(i2,i2) = a(i2,i2) + y1sig
            a(i2,i1) = a(i2,i1) + y2sig
            a(i3,i1) = a(i3,i1) - y2sig
            a(i3,i2) = a(i3,i2) - y1sig
            a(i4,i1) = a(i4,i1) - y2sig
            a(i4,i2) = a(i4,i2) - y1sig
            a(i5,i1) = a(i5,i1) - y2sig
            a(i5,i2) = a(i5,i2) - y1sig
            a(i6,i1) = a(i6,i1) - y2sig
            a(i6,i2) = a(i6,i2) - y1sig
         enddo
      endif
c----------------------------------------------------------
c..   deck 8: exchange, i + j + k ---> l
c..   and  i + j + k ---> l + m
      k1 = k1deck(8)
      k2 = k2deck(8)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            y12sig   = ylo(i1)*ylo(i2)      *sig(k)
            y13sig   = ylo(i1)      *ylo(i3)*sig(k)
            y23sig   =       ylo(i2)*ylo(i3)*sig(k)
            a(i1,i1) = a(i1,i1) + y23sig
            a(i1,i2) = a(i1,i2) + y13sig
            a(i1,i3) = a(i1,i3) + y12sig
            a(i2,i1) = a(i2,i1) + y23sig
            a(i2,i2) = a(i2,i2) + y13sig
            a(i2,i3) = a(i2,i3) + y12sig
            a(i3,i1) = a(i3,i1) + y23sig
            a(i3,i2) = a(i3,i2) + y13sig
            a(i3,i3) = a(i3,i3) + y12sig
            a(i4,i1) = a(i4,i1) - y23sig
            a(i4,i2) = a(i4,i2) - y13sig
            a(i4,i3) = a(i4,i3) - y12sig
            if( i5 .ne. 0 )then
               a(i5,i1) = a(i5,i1) - y23sig
               a(i5,i2) = a(i5,i2) - y13sig
               a(i5,i3) = a(i5,i3) - y12sig
            endif
         enddo
      endif

      return
      end





