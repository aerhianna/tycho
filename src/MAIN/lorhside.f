      subroutine lorhside(sig,sigt,sigv,b,bt,bv,ylo,
     1     nrr,k1deck,k2deck,loz,idebug,lnreac,xidlo,lorlkh)
c..   defines right hand side for instantaneous rate
c     and for solution
c..   specify b array elements for ax = b
c..   derivatives of rhs with respect to T9, V (wda 5/9/99
c..   revised 1/8/00

      implicit none

      include 'dimenfile'

c..xidlo is identifier of nuclei
      character*5 xidlo(lodim)
c..lorlkh is reaction type flag
      character*4	lorlkh(lrdim)
      

      real*8 sig(nreac),sigt(nreac),sigv(nreac)
      real*8 b(ndim),bt(ndim),bv(ndim),ylo(lodim)
      integer*4 nrr(6,nreac),k1deck(11),k2deck(11),loz(lodim)
      integer*4 lnreac

      integer*4 idebug, i, k
      integer*4 k1,k2, i1,i2,i3,i4,i5,i6
      real*8    yysig, ysig, y3sig, weaksum
      data i1/0/,i2/0/,i3/0/,i4/0/,i5/0/,i6/0/
c---------------------------------------------------------------
      do i=1,ndim
         b(i)  = 0.0d0
         bt(i) = 0.0d0
         bv(i) = 0.0d0
      enddo

c..   deck 1: i ---> j
c..   beta decays, positron decays,
c..   i + e ---> j for electron captures (see rate.f)
      weaksum = 0.0d0
      k1 = k1deck(1)
      k2 = k2deck(1)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            b(i1)    = b(i1)    - ylo(i1)*sig(k)
            b(i2)    = b(i2)    + ylo(i1)*sig(k)
            weaksum = weaksum + ( loz(i2)-loz(i1) )*ylo(i1)*sig(k)
c..   dYe/dt is weaksum
            bt(i1)    = bt(i1)    - ylo(i1)*sigt(k)
            bv(i1)    = bv(i1)    - ylo(i1)*sigv(k)
            bt(i2)    = bt(i2)    + ylo(i1)*sigt(k)
            bv(i2)    = bv(i2)    + ylo(i1)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a27,6x,6a12)')'weak: e- e+, i + e ---> j',
     1        'y1*sig','y1','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            write(*,'(i5,a5,2x,a5,2x,a5,5x,a4,1p6e12.3)')
     1           k,xidlo(i1)," --> ",xidlo(i2),lorlkh(k),
     2           ylo(i1)*sig(k),ylo(i1),sig(k)
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
            ysig = ylo(i1)*sig(k)
            b(i1)    = b(i1)    - ysig
            b(i2)    = b(i2)    + ysig
            b(i3)    = b(i3)    + ysig
            bt(i1)   = bt(i1)   - ylo(i1)*sigt(k)
            bv(i1)   = bv(i1)   - ylo(i1)*sigv(k)
            bt(i2)   = bt(i2)   + ylo(i1)*sigt(k)
            bv(i2)   = bv(i2)   + ylo(i1)*sigv(k)
            bt(i3)   = bt(i3)   + ylo(i1)*sigt(k)
            bv(i3)   = bv(i3)   + ylo(i1)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a28,2x,a4,7a10)')'dissociations i ---> j + k',
     1        'src','y*sig','y1','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            ysig = ylo(i1)*sig(k)
            write(*,'(i5,3(a5,1x),a5,2x,a4,1p7e10.2)') k,
     1           xidlo(i1), " --> ", xidlo(i2), xidlo(i3),lorlkh(k),
     1            ysig,ylo(i1),sig(k)
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
            ysig = ylo(i1)*sig(k)
            b(i1)    = b(i1)    - ysig
            b(i2)    = b(i2)    + ysig
            b(i3)    = b(i3)    + ysig
            b(i4)    = b(i4)    + ysig

            bt(i1)    = bt(i1)    - ylo(i1)*sigt(k)
            bv(i1)    = bv(i1)    - ylo(i1)*sigv(k)
            bt(i2)    = bt(i2)    + ylo(i1)*sigt(k)
            bv(i2)    = bv(i2)    + ylo(i1)*sigv(k)
            bt(i3)    = bt(i3)    + ylo(i1)*sigt(k)
            bv(i3)    = bv(i3)    + ylo(i1)*sigv(k)
            bt(i4)    = bt(i4)    + ylo(i1)*sigt(k)
            bv(i4)    = bv(i4)    + ylo(i1)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a33,6x,7a10)')'dissociations, i ---> j + k + l',
     1        'ysig','y1','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            ysig = ylo(i1)*sig(k)
            write(*,'(i4,4(a5,1x),a5,2x,a4,1p7e10.2)') k,
     1           xidlo(i1), " --> ", xidlo(i2), xidlo(i3), xidlo(i4),
     2           lorlkh(k),  ysig,ylo(i1),sig(k)
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
            yysig = ylo(i1)*ylo(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    + yysig

            bt(i1)    = bt(i1)    - ylo(i1)*ylo(i2)*sigt(k)
            bv(i1)    = bv(i1)    - ylo(i1)*ylo(i2)*sigv(k)
            bt(i2)    = bt(i2)    - ylo(i1)*ylo(i2)*sigt(k)
            bv(i2)    = bv(i2)    - ylo(i1)*ylo(i2)*sigv(k)
            bt(i3)    = bt(i3)    + ylo(i1)*ylo(i2)*sigt(k)
            bv(i3)    = bv(i3)    + ylo(i1)*ylo(i2)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a24,9x,7a10)')'captures, i + j ---> k',
     1        'y1*y2*sig','y1','y2','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)
            write(*,'(i4,3(a5,1x),a5,2x,a4,1p7e10.2)') k,
     1           xidlo(i1), xidlo(i2), " --> ", xidlo(i3),lorlkh(k),
     1           yysig,ylo(i1),ylo(i2),sig(k)
c            write(*,'(i4,3(a5,1x),a5,2x,a4,2a1,1p7e10.2)')
c     1           k,
c     1           xid(i1), xid(i2), " --> ", xid(i3),rlkh(k),rnr(k),
c     1           rvw(k),yysig,ylo(i1),ylo(i2),sig(k)
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
            yysig = ylo(i1)*ylo(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    + yysig
            b(i4)    = b(i4)    + yysig

            bt(i1)    = bt(i1)    - ylo(i1)*ylo(i2)*sigt(k)
            bv(i1)    = bv(i1)    - ylo(i1)*ylo(i2)*sigv(k)
            bt(i2)    = bt(i2)    - ylo(i1)*ylo(i2)*sigt(k)
            bv(i2)    = bv(i2)    - ylo(i1)*ylo(i2)*sigv(k)
            bt(i3)    = bt(i3)    + ylo(i1)*ylo(i2)*sigt(k)
            bv(i3)    = bv(i3)    + ylo(i1)*ylo(i2)*sigv(k)
            bt(i4)    = bt(i4)    + ylo(i1)*ylo(i2)*sigt(k)
            bv(i4)    = bv(i4)    + ylo(i1)*ylo(i2)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a33,6x,5a10)')'exchange, i + j ---> k + l',
     1        'yysig','y1','y2','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)
            write(*,'(i4,4(a5,1x),a5,2x,a4,1p5e10.2)') k,
     1           xidlo(i1),xidlo(i2), " --> ", xidlo(i3),xidlo(i4),
     1           lorlkh(k),
     1           yysig,ylo(i1),ylo(i2),sig(k)
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
            yysig = ylo(i1)*ylo(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    + yysig
            b(i4)    = b(i4)    + yysig
            b(i5)    = b(i5)    + yysig

            bt(i1)    = bt(i1)    - ylo(i1)*ylo(i2)*sigt(k)
            bt(i2)    = bt(i2)    - ylo(i1)*ylo(i2)*sigt(k)
            bt(i3)    = bt(i3)    + ylo(i1)*ylo(i2)*sigt(k)
            bt(i4)    = bt(i4)    + ylo(i1)*ylo(i2)*sigt(k)
            bt(i5)    = bt(i5)    + ylo(i1)*ylo(i2)*sigt(k)
            bv(i1)    = bv(i1)    - ylo(i1)*ylo(i2)*sigv(k)
            bv(i2)    = bv(i2)    - ylo(i1)*ylo(i2)*sigv(k)
            bv(i3)    = bv(i3)    + ylo(i1)*ylo(i2)*sigv(k)
            bv(i4)    = bv(i4)    + ylo(i1)*ylo(i2)*sigv(k)
            bv(i5)    = bv(i5)    + ylo(i1)*ylo(i2)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a39,6x,6a10)')'exchange, i + j ---> k + l + m',
     1        'y1*y2*sig','y1','y2','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)
            write(*,'(i4,5(a5,1x),a5,2x,a4,1p6e10.2)') k,
     1           xidlo(i1), xidlo(i2), " --> ", xidlo(i3),xidlo(i4),
     1           xidlo(i5), lorlkh(k),
     1           yysig,ylo(i1),ylo(i2),sig(k)
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
            yysig = ylo(i1)*ylo(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    + yysig
            b(i4)    = b(i4)    + yysig
            b(i5)    = b(i5)    + yysig
            b(i6)    = b(i6)    + yysig


            bt(i1)    = bt(i1)    - ylo(i1)*ylo(i2)*sigt(k)
            bt(i2)    = bt(i2)    - ylo(i1)*ylo(i2)*sigt(k)
            bt(i3)    = bt(i3)    + ylo(i1)*ylo(i2)*sigt(k)
            bt(i4)    = bt(i4)    + ylo(i1)*ylo(i2)*sigt(k)
            bt(i5)    = bt(i5)    + ylo(i1)*ylo(i2)*sigt(k)
            bt(i6)    = bt(i6)    + ylo(i1)*ylo(i2)*sigt(k)
            bv(i1)    = bv(i1)    - ylo(i1)*ylo(i2)*sigv(k)
            bv(i2)    = bv(i2)    - ylo(i1)*ylo(i2)*sigv(k)
            bv(i3)    = bv(i3)    + ylo(i1)*ylo(i2)*sigv(k)
            bv(i4)    = bv(i4)    + ylo(i1)*ylo(i2)*sigv(k)
            bv(i5)    = bv(i5)    + ylo(i1)*ylo(i2)*sigv(k)
            bv(i6)    = bv(i6)    + ylo(i1)*ylo(i2)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a45,6x,7a9)')'exchange, i + j ---> k + l + m + n',
     1        'y1*y2*sig','y1','y2','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)
            write(*,'(i4,6(a5,1x),a5,2x,a4,1p7e9.2)') k,
     1           xidlo(i1), xidlo(i2), " --> ", xidlo(i3),xidlo(i4),
     1           xidlo(i5),xidlo(i6),lorlkh(k),
     1           yysig,ylo(i1),ylo(i2),sig(k)
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

c..   triple alpha checked
            y3sig = ylo(i1)*ylo(i2)*ylo(i3)*sig(k)
            b(i1)    = b(i1)    - y3sig
            b(i2)    = b(i2)    - y3sig
            b(i3)    = b(i3)    - y3sig
            b(i4)    = b(i4)    + y3sig

            bt(i1)    = bt(i1)    - ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)
            bt(i2)    = bt(i2)    - ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)
            bt(i3)    = bt(i3)    - ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)
            bt(i4)    = bt(i4)    + ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)

            bv(i1)    = bv(i1)    - ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)
            bv(i2)    = bv(i2)    - ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)
            bv(i3)    = bv(i3)    - ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)
            bv(i4)    = bv(i4)    + ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)

         enddo
      endif
c----------------------------------------------------
      if( idebug .ne. 0 )then
         write(*,'(a39,6x,6a10)')'exchange, i + j + k ---> l, l+m',
     1        'y123*sig','y1','y2','y3'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            if( i5 .eq. 0 )then
               y3sig = ylo(i1)*ylo(i2)*ylo(i3)*sig(k)
               write(*,'(i4,4(a5,1x),a5,8x,a4,1p6e10.2)') k,
     1              xidlo(i1), xidlo(i2), xidlo(i3), " --> ",xidlo(i4),
     1              lorlkh(k),
     1              y3sig,ylo(i1),ylo(i2),ylo(i3),sig(k)
            else
               y3sig = ylo(i1)*ylo(i2)*ylo(i3)*sig(k)
               write(*,'(i4,5(a5,1x),a5,2x,a4,1p6e10.2)') k,
     1              xidlo(i1), xidlo(i2), xidlo(i3), " --> ",
     1              xidlo(i4),xidlo(i5), lorlkh(k),
     1              y3sig,ylo(i1),ylo(i2),ylo(i3),sig(k)
            endif
         enddo
c..   final idebug command
         stop'rhside'
      endif
c--------------------------------------------------------------
c..   deck 9: exchange, i + j + k ---> l + m
      k1 = k1deck(9)
      k2 = k2deck(9)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            yysig = ylo(i1)*ylo(i2)*ylo(i3)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    - yysig
            b(i4)    = b(i4)    + yysig
            b(i5)    = b(i5)    + yysig

            bt(i1)    = bt(i1)    - ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)
            bt(i2)    = bt(i2)    - ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)
            bt(i3)    = bt(i3)    - ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)
            bt(i4)    = bt(i4)    + ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)
            bt(i5)    = bt(i5)    + ylo(i1)*ylo(i2)*ylo(i3)*sigt(k)
            bv(i1)    = bv(i1)    - ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)
            bv(i2)    = bv(i2)    - ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)
            bv(i3)    = bv(i3)    - ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)
            bv(i4)    = bv(i4)    + ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)
            bv(i5)    = bv(i5)    + ylo(i1)*ylo(i2)*ylo(i3)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a39,6x,6a10)')'exchange, i + j ---> k + l + m',
     1        'y1*y2*sig','y1','y2','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            yysig = ylo(i1)*ylo(i2)*ylo(i3)*sig(k)
            write(*,'(i4,5(a5,1x),a5,2x,a4,1p6e10.2)') k,
     1           xidlo(i1), xidlo(i2), " --> ", xidlo(i3),xidlo(i4),
     1           xidlo(i5), lorlkh(k),
     1           yysig,ylo(i1),ylo(i2),sig(k)
         enddo
      endif

c-------------------------------------------------------------------
c..   deck 10: exchange, i + j + k + l ---> m + n
      k1 = k1deck(10)
      k2 = k2deck(10)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    - yysig
            b(i3)    = b(i3)    - yysig
            b(i4)    = b(i4)    - yysig
            b(i5)    = b(i5)    + yysig
            b(i6)    = b(i6)    + yysig


            bt(i1)    = bt(i1) - ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigt(k)
            bt(i2)    = bt(i2) - ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigt(k)
            bt(i3)    = bt(i3) - ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigt(k)
            bt(i4)    = bt(i4) - ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigt(k)
            bt(i5)    = bt(i5) + ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigt(k)
            bt(i6)    = bt(i6) + ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigt(k)
            bv(i1)    = bv(i1) - ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigv(k)
            bv(i2)    = bv(i2) - ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigv(k)
            bv(i3)    = bv(i3) - ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigv(k)
            bv(i4)    = bv(i4) - ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigv(k)
            bv(i5)    = bv(i5) + ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigv(k)
            bv(i6)    = bv(i6) + ylo(i1)*ylo(i2)*ylo(i3)*ylo(i4)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a45,6x,7a9)')'exchange, i + j ---> k + l + m + n',
     1        'y1*y2*sig','y1','y2','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            yysig = ylo(i1)*ylo(i2)**ylo(i3)*ylo(i4)*sig(k)
            write(*,'(i4,6(a5,1x),a5,2x,a4,1p7e9.2)') k,
     1           xidlo(i1), xidlo(i2), " --> ", xidlo(i3),xidlo(i4),
     1           xidlo(i5),xidlo(i6),lorlkh(k),
     1           yysig,ylo(i1),ylo(i2),sig(k)
         enddo
      endif

c--------------------------------------------------------------
c..   deck 11: decay, i ---> j + k + l + m
      k1 = k1deck(11)
      k2 = k2deck(11)
      if( k1 .ne. 0 .and. k2 .ne. 0 )then
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            yysig = ylo(i1)*sig(k)
            b(i1)    = b(i1)    - yysig
            b(i2)    = b(i2)    + yysig
            b(i3)    = b(i3)    + yysig
            b(i4)    = b(i4)    + yysig
            b(i5)    = b(i5)    + yysig

            bt(i1)    = bt(i1)    - ylo(i1)*sigt(k)
            bt(i2)    = bt(i2)    + ylo(i1)*sigt(k)
            bt(i3)    = bt(i3)    + ylo(i1)*sigt(k)
            bt(i4)    = bt(i4)    + ylo(i1)*sigt(k)
            bt(i5)    = bt(i5)    + ylo(i1)*sigt(k)
            bv(i1)    = bv(i1)    - ylo(i1)*sigv(k)
            bv(i2)    = bv(i2)    + ylo(i1)*sigv(k)
            bv(i3)    = bv(i3)    + ylo(i1)*sigv(k)
            bv(i4)    = bv(i4)    + ylo(i1)*sigv(k)
            bv(i5)    = bv(i5)    + ylo(i1)*sigv(k)
         enddo
      endif
      if( idebug .ne. 0 )then
         write(*,'(a39,6x,6a10)')'exchange, i + j ---> k + l + m',
     1        'y1*y2*sig','y1','y2','sig'
         do k = k1,k2
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            yysig = ylo(i1)*sig(k)
            write(*,'(i4,5(a5,1x),a5,2x,a4,1p6e10.2)') k,
     1           xidlo(i1), xidlo(i2), " --> ", xidlo(i3),xidlo(i4),
     1           xidlo(i5), lorlkh(k),
     1           yysig,ylo(i1),ylo(i2),sig(k)
         enddo
      endif

      return
      end


