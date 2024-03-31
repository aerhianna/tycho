c     
c     
c     
      subroutine xcheck(kk,x)
c     -----------------------------
c     REVISION HISTORY:
c     7-28-06
c     -----------------------------
c     checks conservation of masss 
c     in composition changes
c     -----------------------------
            
      implicit none

      include 'dimenfile'
      
      integer*4 k,kk,n,kworst

      real*8    x(ndim,kdm)
      real*8    dum1(kdm), dum2(kdm)
      real*8    xworst
      include 'cburn'
      
      
c     -------------------------------------------------
c     check convervation over k (nucleon and charge)
c     prints worst value to standard i/o if .gt. 1.0e-8
c     -------------------------------------------------
      
      do k = 2, kk
         dum1(k) = -1.0d0
         dum2(k) = -x(nnuc+1,k)
         do n = 1, nnuc
            dum1(k) = dum1(k) + x(n,k)
            dum2(k) = dum2(k) + dble( lz(n) )*x(n,k)/xa(n)
         enddo
      enddo

      xworst = 0.0d0
      kworst = 0
      do k = 2, kk
         if( abs( dum1(k) ) .gt. abs( xworst ) )then
            xworst = dum1(k)
            kworst = k
         endif
      enddo

      if( abs( xworst ) .gt. 1.0d-8  )then
c     this level is more that truncation in model 
c     format would imply
         write(*,'(a12,a8,i5,2(a12,1pe12.3))') 'xcheck: ',
     .        'zone',kworst,'err(mass)', dum1(kworst), 
     .        'd(Ye)', dum2(kworst)
      endif
      
      return      
      end

