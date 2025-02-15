      subroutine loburn(dtstar,k,it,kk,leq)

      implicit none

c..   smaller network for low temperatures
c..   automatically finds rates and links using data
c..   for larger network
c..   always does network solution

      include 'dimenfile'
      include 'cconst'
      include 'crate'
      include 'comcsolve'
      include 'caeps'
      include 'cdtnuc'
      include 'cdeuter'

      integer*4 k, j, i, i1, i2, i3, i4, i5, i6
      integer*4 idebug, leq
      integer*4 idk
      integer*4 it, kk

      real*8    t9, rho, sum, ye
      real*8 eb,dn
      real*8 dtstar
      real*8 fak,ylold(lodim)

      logical isit

      data idebug/0/

      save
c------------------------------------------------------

      t9  = 1.0d-9 * temp(k)
      rho = den(k)

c..   low temperature network solution

      do j = 1, ndim
         b(j) = 0.0d0
         y(j) = aex(j,k)
      enddo

c..   ye is electrons per nucleon
      ye = 0.0d0
      do j = 1, ndim-1
         ye = ye + y(j)*dble(nz(j))
      enddo
c..   initialize; b's assumed zero
      do j = 1, ndim
         b(j) = 0.0d0
         bt(j) = 0.0d0
         bv(j) = 0.0d0
      enddo


c..   idlo(1) is set to -1 initially, in solven.f
c..   if initial call, look for reduced set of nuclei 
c..   defined by loz, lon
c..   and define idlo(j) which gives index of nucleus in larger set
c..   for gather-scatter operations
      if( idlo(1) .lt. 0 )then

c         partfun(t9)

         do j = 1, ndim
            do i = 1, lodim-1
               if( nz(j) .eq. loz(i) .and. nn(j) .eq. lon(i) )then
                  idlo(i)  = j
                  xidlo(i) = xid(j)
                  lopf(i) = pf(j)
               endif

            enddo
         enddo
         idlo(lodim)  = ndim
         xidlo(lodim) = '  Ye'

         do i = 1, lodim
            if( idlo(i) .le. 0 )then
               write(*,*)'species ',i,' z n ',loz(i),lon(i),
     1              ' not in network'
               stop'loburn species test'
            endif
         enddo

c..   nuclear abundance: short list = ylo
         do i = 1, lodim
            ylo(i) = y(idlo(i))
         enddo

c..   write(*,*)'search for reactions among nuclei in short list'
c..   leave out reactions that should be small at low temperatures
c..   this is to avoid bad extrapolation of rates
         do j = 1,11
            l1deck(j) = 0
            l2deck(j) = 0
         enddo
         lnreac = 0
         do idk = 1,11
c..   write(*,*)'deck ',idk
            do j = k1deck(idk),k2deck(idk)
               isit = .false.
               i1 = 0
c..   is first nucleus in short list?
               do i = 1, lodim
                  if( idlo(i) .eq. nrr(1,j) )then
                     isit = .true.
                     i1 = i
                  endif
               enddo
c..   is second nucleus in short list?
               i2 = 0
               do i = 1, lodim
                  if( idlo(i) .eq. nrr(2,j) )then
                     i2 = i
                  endif
               enddo
c..   and so on...
               i3 = 0
               do i = 1, lodim
                  if( idlo(i) .eq. nrr(3,j) )then
                     i3 = i
                  endif
               enddo
               i4 = 0
               do i = 1, lodim
                  if( idlo(i) .eq. nrr(4,j) )then
                     i4 = i
                  endif
               enddo
               i5 = 0
               do i = 1, lodim
                  if( idlo(i) .eq. nrr(5,j) )then
                     i5 = i
                  endif
               enddo
               i6 = 0
               do i = 1, lodim
                  if( idlo(i) .eq. nrr(6,j) )then
                     i6 = i
                  endif
               enddo
               if( idk .eq. 1 .and. i1*i2 .ne. 0 )then
                  lnreac = lnreac + 1
                  lodeck(lnreac) = 1
                  lorr(lnreac) = j
                  lonrr(1,lnreac) = i1
                  lonrr(2,lnreac) = i2
c      write(*,'(2i5,7a5)')j,lnreac,
c     1                 rname(1,j),' --> ',rname(2,j)
               elseif( idk .eq. 2 .and. i1*i2*i3 .ne. 0 )then
                  lnreac = lnreac + 1
                  lodeck(lnreac) = 2
                  lorr(lnreac) = j
                  lonrr(1,lnreac) = i1
                  lonrr(2,lnreac) = i2
                  lonrr(3,lnreac) = i3

c      write(*,'(2i5,7a5)')j,lnreac,
c     1                 rname(1,j),' --> ',rname(2,j),rname(3,j)
               elseif( idk .eq. 3 .and. i1*i2*i3*i4 .ne. 0 )then
c..no photodissociation of Z .ge. 6 = carbon and above
                  if( loz(i1) .lt. 6 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 3
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
                  endif
c      write(*,'(2i5,7a5)')j,lnreac,
c     1                 rname(1,j),' --> ',rname(2,j),rname(3,j),
c     2                 rname(4,j)
               elseif( idk .eq. 4 .and. i1*i2*i3 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
                  if( loz(i1) .ne. 2 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 4
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
c      write(*,'(2i5,7a5)')j,lnreac,
c     1                 rname(1,j),rname(2,j),' --> ',rname(3,j)
c..no fusion of He and Z>2, with species heavier than Li Z=3
                  elseif( loz(i2) .lt. 4 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 4
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
c      write(*,'(2i5,7a5)')j,lnreac,
c     1                 rname(1,j),rname(2,j),' --> ',rname(3,j)
                  endif
               elseif( idk .eq. 5 .and. i1*i2*i3*i4 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
c..avoid C+C reactions

                  if( loz(i1) .ne. 2 .and. loz(i1) .lt. 6 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 5
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
c      write(*,'(2i5,7a5)')j,lnreac,
c     1                 rname(1,j),rname(2,j),' --> ',rname(3,j),
c     2                 rname(4,j)
                  elseif( loz(i2) .lt. 4 )then
c..gets He3, He4 reactions on light nuclei
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 5
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
c      write(*,'(2x,2i5,7a5)')j,lnreac,
c     1                 rname(1,j),rname(2,j),' --> ',rname(3,j),
c     2                 rname(4,j)
                  endif
               elseif( idk .eq. 6 .and. i1*i2*i3*i4*i5 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
                  if( loz(i1) .ne. 2 )then
                  lnreac = lnreac + 1
                  lodeck(lnreac) = 6
                  lorr(lnreac) = j
                  lonrr(1,lnreac) = i1
                  lonrr(2,lnreac) = i2
                  lonrr(3,lnreac) = i3
                  lonrr(4,lnreac) = i4
                  lonrr(5,lnreac) = i5
c                  write(*,'(2i5,7a5)')j,lnreac,
c     1                 rname(1,j),rname(2,j),' --> ',rname(3,j),
c     2                 rname(4,j),rname(5,j)
                  elseif( loz(i2) .lt. 4 )then
                  lnreac = lnreac + 1
                  lodeck(lnreac) = 6
                  lorr(lnreac) = j
                  lonrr(1,lnreac) = i1
                  lonrr(2,lnreac) = i2
                  lonrr(3,lnreac) = i3
                  lonrr(4,lnreac) = i4
                  lonrr(5,lnreac) = i5
c                  write(*,'(2i5,7a5)')j,lnreac,
c     1                 rname(1,j),rname(2,j),' --> ',rname(3,j),
c     2                 rname(4,j),rname(5,j)
                  endif
               elseif( idk .eq. 7 .and. i1*i2*i3*i4*i5*i6 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
                  if( loz(i1) .ne. 2 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 7
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
                     lonrr(5,lnreac) = i5
                     lonrr(6,lnreac) = i6
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),' --> ',rname(3,j),
c     2                    rname(4,j),rname(5,j),rname(6,j)
                  elseif( loz(i2) .lt. 4 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 7
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
                     lonrr(5,lnreac) = i5
                     lonrr(6,lnreac) = i6                    
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),' --> ',rname(3,j),
c     2                    rname(4,j),rname(5,j),rname(6,j)
                  endif
               elseif( idk .eq. 8 .and. i1*i2*i3*i4 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
                  if( loz(i1) .ne. 2 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 8
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),' --> ',rname(3,j),
c     2                    rname(4,j),rname(5,j),rname(6,j)
                  elseif( loz(i2) .lt. 4 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 8
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),' --> ',rname(3,j),
c     2                    rname(4,j),rname(5,j),rname(6,j)
                  endif
               elseif( idk .eq. 9 .and. i1*i2*i3*i4*i5 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
                  if( loz(i1) .ne. 2 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 9
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
                     lonrr(5,lnreac) = i5
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),rname(3,j),' --> ',
c     2                    rname(4,j),rname(5,j)
                  elseif( loz(i2) .lt. 4 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 9
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
                     lonrr(5,lnreac) = i5
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),rname(3,j),' --> ',
c     2                    rname(4,j),rname(5,j)
                  endif
               elseif( idk .eq. 10 .and. i1*i2*i3*i4*i5*i6 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
                  if( loz(i1) .ne. 2 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 10
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
                     lonrr(5,lnreac) = i5
                     lonrr(6,lnreac) = i6
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),rname(3,j),' --> ',
c     2                    rname(4,j),rname(5,j)
                  elseif( loz(i2) .lt. 4 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 10
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
                     lonrr(5,lnreac) = i5
                     lonrr(6,lnreac) = i6
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),rname(3,j),' --> ',
c     2                    rname(4,j),rname(5,j)
                  endif
               elseif( idk .eq. 11 .and. i1*i2*i3*i4*i5 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
                  if( loz(i1) .ne. 2 )then
                     lnreac = lnreac + 1
                     lodeck(lnreac) = 11
                     lorr(lnreac) = j
                     lonrr(1,lnreac) = i1
                     lonrr(2,lnreac) = i2
                     lonrr(3,lnreac) = i3
                     lonrr(4,lnreac) = i4
                     lonrr(5,lnreac) = i5
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),rname(3,j),' --> ',
c     2                    rname(4,j),rname(5,j)
c                  elseif( loz(i2) .lt. 4 )then
c                     lnreac = lnreac + 1
c                     lodeck(lnreac) = 11
c                     lorr(lnreac) = j
c                     lonrr(1,lnreac) = i1
c                     lonrr(2,lnreac) = i2
c                     lonrr(3,lnreac) = i3
c                     lonrr(4,lnreac) = i4
c                     lonrr(5,lnreac) = i5
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),rname(3,j),' --> ',
c     2                    rname(4,j),rname(5,j)
                  endif
               endif
c               elseif( idk .eq. 8 .and. i5 .eq. 0 .and.
c     1                 i1*i2*i3*i4 .ne. 0 )then
c..restrict alpha reactions to low Z targets only
c                  if( loz(i1) .ne. 2 )then
c                     lnreac = lnreac + 1
c                     lodeck(lnreac) = 8
c                     lorr(lnreac)    = j
c                     lonrr(1,lnreac) = i1
c                     lonrr(2,lnreac) = i2
c                     lonrr(3,lnreac) = i3
c                     lonrr(4,lnreac) = i4
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),rname(3,j),' --> ',
c     2                    rname(4,j)
c                  elseif( loz(i2) .lt. 4 )then
c                     lnreac = lnreac + 1
c                     lodeck(lnreac) = 8
c                     lorr(lnreac)    = j
c                     lonrr(1,lnreac) = i1
c                     lonrr(2,lnreac) = i2
c                     lonrr(3,lnreac) = i3
c                     lonrr(4,lnreac) = i4
c                     write(*,'(2i5,7a5)')j,lnreac,
c     1                    rname(1,j),rname(2,j),rname(3,j),' --> ',
c     2                    rname(4,j)
c                  endif
c               endif
            enddo
            if( lnreac .gt. 0 )then
               if( idk .eq. 1 )then
                  l1deck(idk) = idk
                  l2deck(idk) = lnreac
               else
                  l1deck(idk) = l2deck(idk-1)+1
                  l2deck(idk) = lnreac
               endif
            endif

c            if( idk .eq. 1 )then
c               write(*,'(a40)')'rates chosen for loburn network'
c            endif
c            write(*,'(a10,i5)')'ideck',idk

            do i = l1deck(idk),l2deck(idk)
c..   definitions and tests both

c      if( idk .eq. 1 )then
c      write(*,'(5i5,3a5,1p7e10.2)')i,lodeck(i),lorr(i),
c     1              lonrr(1,i),lonrr(2,i),rname(1,lorr(i)),
c     2              ' --> ',rname(2,lorr(i)),
c     3              rcoef(1,lorr(i))
c      elseif( idk .eq. 2 )then
c      write(*,'(5i5,4a5,1p7e10.2)')i,lodeck(i),lorr(i),
c     1              lonrr(1,i),lonrr(2,i),rname(1,lorr(i)),
c     2              ' --> ',rname(2,lorr(i)),rname(3,lorr(i)),
c     3              (rcoef(j,lorr(i)),j=1,7)
c      elseif( idk .eq. 3 )then
c      write(*,'(5i5,5a5,1p7e10.2)')i,lodeck(i),lorr(i),
c     1              lonrr(1,i),lonrr(2,i),rname(1,lorr(i)),
c     2              ' --> ',rname(2,lorr(i)),rname(3,lorr(i)),
c     3              rname(4,lorr(i)),rcoef(1,lorr(i))
c      elseif( idk .eq. 4 )then
c      write(*,'(5i5,4a5,1p7e10.2)')i,lodeck(i),lorr(i),
c     1              lonrr(1,i),lonrr(2,i),rname(1,lorr(i)),
c     2              rname(2,lorr(i)),' --> ',rname(3,lorr(i)),
c     3              rcoef(1,lorr(i))
c      elseif( idk .eq. 5 )then
c      write(*,'(5i5,5a5,1p7e10.2)')i,lodeck(i),lorr(i),
c     1              lonrr(1,i),lonrr(2,i),rname(1,lorr(i)),
c     2              rname(2,lorr(i)),' --> ',rname(3,lorr(i)),
c     3              rname(4,lorr(i)),rcoef(1,lorr(i))
c      elseif( idk .eq. 6 )then
c      write(*,'(5i5,6a5,1p7e10.2)')i,lodeck(i),lorr(i),
c     1              lonrr(1,i),lonrr(2,i),rname(1,lorr(i)),
c     2              rname(2,lorr(i)),' --> ',rname(3,lorr(i)),
c     3              rname(4,lorr(i)),rname(5,lorr(i)),rcoef(1,lorr(i))
c      elseif( idk .eq. 7 )then
c      write(*,'(5i5,7a5,1p7e10.2)')i,lodeck(i),lorr(i),
c     1              lonrr(1,i),lonrr(2,i),rname(1,lorr(i)),
c     2              rname(2,lorr(i)),' --> ',rname(3,lorr(i)),
c     3              rname(4,lorr(i)),rname(5,lorr(i)),rname(6,lorr(i)),
c     4              rcoef(1,lorr(i))
c      elseif( idk .eq. 8 )then
c      write(*,'(5i5,6a5,1p7e10.2)')i,lodeck(i),lorr(i),
c     1              lonrr(1,i),lonrr(2,i),rname(1,lorr(i)),
c     2              rname(2,lorr(i)),rname(3,lorr(i)),' --> ',
c     3              rname(4,lorr(i)),rname(5,lorr(i)),
c     4              rcoef(1,lorr(i))
c      endif

               do j = 1,7
                  locoef(j,i) = rcoef(j,lorr(i))
               enddo
               lorlkh(i) = rlkh(lorr(i))
               loec(i) = ec(lorr(i))
            enddo
         enddo

c..revise al26 to lab value at low temperatures
         write(*,*)'revising al26 to lab value at low temperatures'
         do i = l1deck(1), l2deck(1)
            if(  loz(lonrr(1,i)) .eq. 13 .and. lon(lonrr(1,i)) .eq.
     1           13 )then
               do j = 1,7
                  locoef(j,i) = 0.0d0
               enddo
c..ln( 1/meanlife )
               locoef(1,i) = -33.423d0
            endif
         enddo

c..   echo low burn nuclei
         write(*,*)'nuclei in low Temperature network'
         write(*,'(10(2x,i3,a5))')(i,xidlo(i),i=1,lodim)

      endif

c..   end of initialization...................................

c..   update short list array of abundances
      do i = 1, lodim
         ylo(i) = y(idlo(i))
c..   save initial value
         ylold(i) = ylo(i)
      enddo

      call lorate(t9,rho,ye)

      call lorhside(sig,sigt,sigv,b,bt,bv,ylo,
     1     lonrr,l1deck,l2deck,loz,idebug,lnreac,xidlo,lorlkh)

c..   solve reduced network
      call lojacob(sig,ylo,a,
     1     lonrr,l1deck,l2deck,loz,idebug)

      do j = 1,lonuc
         a(j,j) = a(j,j) + 1.0d0/dtstar
      enddo

c..   loleqs solves linear equations(leqa with lower dimension)
      call leqs(a,b,lonuc,lonuc)

c..   find fastest changing nucleus
c..   test for excessive changes
      sum = 0.0d0
      do j = 1, lonuc
         sum = sum + b(j)*dble( loz(j) + lon(j) )
      enddo

ccccccccccccccccccccccccccccccccccccccc
c      if( abs(sum) .gt. 1.0d-14 )then
ccccccccccccccccccccccccccccccccccccccc
      if( abs(sum) .gt. 1.0d-6 )then
         write(*,'(a35,a10,i5,2(a10,1pe10.2))')
     1        'loburn warning: (sum Xi .ne. 1) ',
     1        'zone',k,'(sum x)-1',sum,
     2        'T(k)',t9*1.0d9
      endif

c..   update
      do i = 1, lodim-1
         ylo(i) = ylo(i) + b(i)
         ylo(i) = dmax1(0.0d0,ylo(i))
c..   store updated value in full array
         j    = idlo(i)
         y(j) = ylo(i)
      enddo

      do j = 1, itot
         if( y(j) .lt. 1.0d-50 ) y(j) = 0.0d0
      enddo
c..   calculate energy generation rate using mass excesses
c..   DNE is 3/2 kT for new particles
      eb   = 0
      dn   = 0
      do j = 1, lonuc
         fak = ylo(j) - ylold(j)
         eb  = eb - qq(idlo(j))*fak
         dn  = dn +       fak
      enddo
      eb     = eb *   9.65d+17    /dtstar
      dn     = dn * 1.2476d+17 *t9/dtstar
c..   eb      = eb * ergspermev * avagadro    /dtstar
c..   dn      = dn * 0.12919d0 * ergspermev * avagadro *t9/dtstar

      aeps(k) = eb - dn

c..   save new values
      if( leq .eq. 0 )then
         do j = 1, itot
            aex(j,k) = y(j)
         enddo
         aex(ndim,k) = 0.0d0
         do j = 1, itot
            aex(ndim,k) = aex(ndim,k)
     1           + aex(j,k)*dble( nz(j) )
         enddo
      endif

      return
      end


