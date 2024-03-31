      subroutine select(xx,tt,vv,epsilon,jnbfkt,kk,nc,iabflag,iequil,
     1     jnb)

c..select faster flows and list them

      implicit none

      include 'dimenfile'
      include 'crate'
      include 'comcsolve'
      include 'comlink'
      include 'cdeuter'
      include 'cburn'
      include 'ceoset'
      include 'cconst'

      integer*4 jnb
      integer*4 j, iequil
      integer*4 i, k, kk, nc, maxk, i1,i2,i3,i4,i5,i6, iabflag
      real*8    xx(ndim,kdm),tt(2,kdm),dd(2,kdm),vv(2,kdm)
      real*8    rlmax,epsilon,ysig, yysig, y3sig
      real*8    erate
      real*8    sum,dif,fak,ratio,raty,equil,a1,a2,a3,a4,t9,rho
c..   fkt version of jnb rates
      real*8    jnbfkt(13)
      real*8    t953,t943,t923,t913,t932,t976,t965,cf88
      real*8    t9a,t9a13,t9a56

      character*5 blank
      data blank/'     '/
      save
c--------------------------------------------------------------------
c..   get abundances (mole fractions)
c..   wda march 5, 2008
      do k = 1, netsize
         y(k) = xx(k,kk)/xa(k)
         x(k) = xx(k,kk)
      enddo

      y(ndim) = 0.0d0
      do k = 1, netsize
         y(ndim) = y(ndim) + xx(k,kk)*dble(nz(k) )/xa(k)
      enddo

c..   rhoz is the effective nucleon density in grams/cc used by rate.f
c..   use rhoz from state.f; assumes a previous call to state.f
      dd(nc,kk) = rhoz(kk)

      if( tt(nc,kk) .gt. tburnlo*1.0d9 )then

c..   get flows for reaction links for full network

         do k = k1deck(1),k2deck(1)
            i1 = nrr(1,k)
            rlink(k) = y(i1)*sig(k)
         enddo
         do k = k1deck(2),k2deck(2)
            i1 = nrr(1,k)
            ysig = y(i1)*sig(k)
            rlink(k) = ysig
         enddo
         do k = k1deck(3),k2deck(3)
            i1 = nrr(1,k)
            ysig = y(i1)*sig(k)
            rlink(k) = ysig
         enddo
         do k = k1deck(4),k2deck(4)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            yysig = y(i1)*y(i2)*sig(k)  
            rlink(k) = yysig
         enddo
         do k = k1deck(5),k2deck(5)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            yysig = y(i1)*y(i2)*sig(k)  
            rlink(k) = yysig
         enddo
         do k = k1deck(6),k2deck(6)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            yysig = y(i1)*y(i2)*sig(k)   
            rlink(k) = yysig
         enddo
         do k = k1deck(7),k2deck(7)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            yysig = y(i1)*y(i2)*sig(k)   
            rlink(k) = yysig
         enddo
         do  k = k1deck(8),k2deck(8)
            i1 = nrr(1,k)
            y3sig = y(i1)**3*sig(k)/3.0d0
            rlink(k) = y3sig
         enddo

c..   fastest reaction link
         rlmax = 0.0d0
         do k = k1deck(1),k2deck(8)
            if( rlink(k) .gt. rlmax )then
               rlmax = rlink(k)
               maxk = k
            endif
         enddo
         write(*,'(a15,2a8,a12,15x,a25)')
     1        'fastest flow','index','deck','reaction',
     2        'energy release (erg/g/s)'

         erate = rlink(maxk)*qval(maxk)*cergs

         write(*,'(1pe15.3,2i8,6a5,1pe12.3)')rlmax,maxk,ideck(maxk),
     1        (rname(i,maxk),i=1,6),erate

         write(*,'(/2(a15,1pe11.3),a15,i5/)')'T(K)',tt(nc,kk),
     1        'rho(g/cc)',dd(nc,kk),'zone',kk


c..   analysis for equilibria
         if( iequil .eq. 1 .and. tt(nc,kk) .gt. 1.0d9 )then
            t9  = tt(nc,kk)*1.0d-9
            rho = dd(nc,kk)
c..   three body reactions and reverse reactions
            write(*,*)'decks 2 (1-->2) and 4 (2-->1)'

            write(*,'(2a6,a15,2x,a15,12a10)')'j','k',
     1           'forward','reverse',
     1           'f rate', 'r rate','sum','dif','r/f',
     2           'y2*y3/y1','sigr/sigf','equil','rat/eq',
     3           'f1*f2/f3'

            do j = k1deck(4), k2deck(4)
               k = irev(j)

c..   k is reverse, j forward so 1+2-->3
               if( k .gt. 0 )then
                  i1 = nrr(1,j)
                  i2 = nrr(2,j)
                  i3 = nrr(3,j)
                  a1 = anuc(i1)
                  a2 = anuc(i2)
                  a3 = anuc(i3)
                  if( y(i1) .ne. 0.0d0 )then
                     raty = y(i1)*y(i2)/y(i3)
                  else
                     raty = 0
                  endif
                  sum = rlink(j) + rlink(k)
                  dif = rlink(j) - rlink(k)
                  if( rlink(j) .ne. 0.0d0 )then
                     fak = rlink(k)/rlink(j)
                  else
                     fak = 0.0d0
                  endif
                  if( sig(j) .ne. 0.0d0 )then
                     ratio = sig(k)/sig(j)
                  else
                     ratio = 0.0d0
                  endif
                  equil = 9.8678d9
     1                 *pf(i1)*pf(i2)/pf(i3) 
     2                 *( a1*a2/a3 * t9 )**1.5d0/ rho
     3                 *exp( -11.605d0*qval(j)/t9 )
                  
                  if( rlink(j) .gt. epsilon*rlmax )then
                     write(*,'(2i6,3a5,2x,3a5,1p12e10.2)')
     1                    k,j,
     1                    (rname(i,j),i=1,3),
     2                    (rname(i,k),i=1,3),
     3                    rlink(j),rlink(k),sum,dif,fak,raty,
     4                    ratio,equil,ratio/equil,
     5                    pf(i1)*pf(i2)/pf(i3)
cccccccccccccccc
                  endif
               endif
            enddo

            write(*,'(2a6,a15,2x,a15,12a10)')'j','k',
     1           'forward','reverse',
     1           'f rate', 'r rate','sum','dif','r/f',
     2           'y2*y3/y1','sigr/sigf','equil','rat/eq',
     3           'f1*f2/f3'

c..   four body reactions and reverse reactions
c..   four body binary
c..   deck 5: 1+2-->3+4

c..   four body binary reactions and reverse reactions
            write(*,*)'deck 5 (2-->2) '

            write(*,'(2a6,a20,2x,a20,12a10)')'k','j','forward',
     1           'reverse',
     1           'f rate', 'r rate','sum','dif','r/f',
     2           'y12/y34','sigr/sigf','equil','rat/eq',
     3           'f12/f34'
            do k = k1deck(5), k2deck(5)
               j = irev(k)
               if( j .gt. 0 )then
c..   k is forward j is reverse,so 1+2-->3+4
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  i3 = nrr(3,k)
                  i4 = nrr(4,k)
                  a1 = dble( nz(nrr(1,k)) + nn(nrr(1,k) ) )
                  a2 = dble( nz(nrr(2,k)) + nn(nrr(2,k) ) )
                  a3 = dble( nz(nrr(3,k)) + nn(nrr(3,k) ) )
                  a4 = dble( nz(nrr(4,k)) + nn(nrr(4,k) ) )

                  if( y(i3) .ne. 0.0d0 .and. y(i4) .ne. 0.0d0)then
                     raty = y(i1)*y(i2)/y(i3)/y(i4)
                  else
                     raty = 12345.0d20
                  endif
                  sum = rlink(k) + rlink(j)
                  dif = rlink(k) - rlink(j)
                  if( rlink(j) .ne. 0.0d0 )then
                     fak = rlink(j)/rlink(k)
                  else
                     fak = 0.0d0
                  endif
                  if( sig(k) .ne. 0.0d0 )then
                     ratio = sig(j)/sig(k)
                  else
                     ratio = 0.0d0
                  endif
                  equil = pf(i1)*pf(i2)/pf(i3)/pf(i4)
     2                 *( a1*a2/a3/a4 )**1.5d0
     3                 *exp( -11.605d0*qval(k)/t9 )
               
                  if( rlink(j) .gt. epsilon*rlmax )then
                     write(*,'(2i6,4a5,2x,4a5,1p12e10.2)')k,j,
     1                    (rname(i,k),i=1,4),
     2                    (rname(i,j),i=1,4),
     3                    rlink(k),rlink(j),sum,dif,fak,raty,
     4                    ratio,equil,ratio/equil,
     5                    pf(i1)*pf(i2)/pf(i3)/pf(i4)
                  endif
               endif
            enddo
         
            write(*,'(2a6,a20,2x,a20,12a10)')'k','j','forward',
     1           'reverse',
     1           'f rate', 'r rate','sum','dif','r/f',
     2           'y12/y34','sigr/sigf','equil','rat/eq',
     3           'f12/f34'
            
c..   four body asymmetric number reactions
c..   deck 3: 1--> 2+3+4
c..   deck 8: 1+2+3-->4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         elseif( iequil .eq. 2 )then
            write(*,'(/a40/)')'steady state not implimented yet'
            stop'select 1'
         else
            if( iequil .ne. 0 )then
               write(*,'(/a70)')
     1              'temperature too low for iequil .ne. 0'
               write(*,'(a70/)')
     1              'no equilibrium analysis attempted'
            endif
             
         endif




c..   edited list of rates
         write(*,'(a4,a2,a10,12x,a8,22x,a6,
     1        5(a7,2x),2a10)') 'i ','dk','flow','reaction','source',
     2        'sig','y1','y2',' ',' ','erate','Q(Mev)'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c..   re-check conversion from dN/dt to epsilon for 3-alpha (Iben)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do k = 1, ireac

            in(k)   = 0
            ldot(k) = 1
            if( rlink(k) .gt. epsilon*rlmax )then
               if( ideck(k) .eq. 1 )then
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  rlink(k) = y(i1)*sig(k)
                  erate = rlink(k)*qval(k)*cergs

                  ldot(k) = 4
                  in(k)   = ideck(k)
                  x1(1,k) = nn(i1)
                  y1(1,k) = nz(i1)
                  x1(2,k) = nn(i2)
                  y1(2,k) = nz(i2)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)

                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p2e9.2,27x,1p2e10.3)') 
     1                 k,ideck(k),rlink(k),blank,
     1                 xid(i1), " --> ", xid(i2),blank,blank,blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),y(i1),erate,qval(k)


               elseif( ideck(k) .eq. 2 )then
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  i3 = nrr(3,k)
                  rlink(k) = y(i1)*sig(k)
                  erate = rlink(k)*qval(k)*cergs

                  in(k) = ideck(k)
                  x1(1,k) = nn(i1)
                  y1(1,k) = nz(i1)
                  x1(2,k) = nn(i3)
                  y1(2,k) = nz(i3)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)

                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p2e9.2,27x,1p2e10.3)') 
     1                 k,ideck(k),rlink(k),blank,
     1                 xid(i1), " --> ", xid(i2), xid(i3),blank,
     1                 blank,rlkh(k),rnr(k),rvw(k),
     1                 sig(k),y(i1),erate,qval(k)

               elseif( ideck(k) .eq. 3 )then
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  i3 = nrr(3,k)
                  i4 = nrr(4,k)
                  rlink(k) = y(i1)*sig(k)
                  erate = rlink(k)*qval(k)*cergs

                  in(k) = ideck(k)
                  x1(1,k) = nn(i1)
                  y1(1,k) = nz(i1)
                  x1(2,k) = nn(i4)
                  y1(2,k) = nz(i4)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)

                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p2e9.2,27x,1p2e10.3)') 
     1                 k,ideck(k),rlink(k),blank,
     1                 xid(i1), " --> ", xid(i2), xid(i3), xid(i4),
     1                 blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),y(i1),erate,qval(k)

               elseif( ideck(k) .eq. 4 )then
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  i3 = nrr(3,k)
                  rlink(k) = y(i1)*y(i2)*sig(k)
                  erate = rlink(k)*qval(k)*cergs

                  in(k) = ideck(k)
                  x1(1,k) = nn(i2)
                  y1(1,k) = nz(i2)
                  x1(2,k) = nn(i3)
                  y1(2,k) = nz(i3)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)
                  if( rlkh(k).eq.'  ec' .or. rlkh(k).eq.'bet+')then
                     ldot(k) = 4
                  endif
                  
                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p3e9.2,18x,1p2e10.3)') 
     1                 k,ideck(k),rlink(k),
     1                 xid(i1), xid(i2), " --> ", xid(i3),
     1                 blank,blank,blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),y(i1),y(i2),erate,qval(k)

               elseif( ideck(k) .eq. 5 )then
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  i3 = nrr(3,k)
                  i4 = nrr(4,k)
                  rlink(k) = y(i1)*y(i2)*sig(k)
                  erate = rlink(k)*qval(k)*cergs

                  in(k) = ideck(k)
                  x1(1,k) = nn(i2)
                  y1(1,k) = nz(i2)
                  x1(2,k) = nn(i4)
                  y1(2,k) = nz(i4)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)
                  
                  write(*,'(i4,i2,1pe10.2,7(a5,1x)a4,2a1,
     1                 1p3e9.2,18x,1p2e10.3)') 
     1                 k,ideck(k),rlink(k),
     1                 xid(i1), xid(i2), " --> ", xid(i3),xid(i4),
     1                 blank,blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),y(i1),y(i2),erate,qval(k)

               elseif( ideck(k) .eq. 6 )then
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  i3 = nrr(3,k)
                  i4 = nrr(4,k)
                  i5 = nrr(5,k)
                  rlink(k) = y(i1)*y(i2)*sig(k)
                  erate = rlink(k)*qval(k)*cergs

                  in(k) = ideck(k)
                  x1(1,k) = nn(i2)
                  y1(1,k) = nz(i2)
                  x1(2,k) = nn(i5)
                  y1(2,k) = nz(i5)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)
                  
                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p3e9.2,18x,1p2e10.3)') 
     1                 k,ideck(k),rlink(k),
     1                 xid(i1), xid(i2), " --> ", xid(i3),xid(i4),
     1                 xid(i5),blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),y(i1),y(i2),erate,qval(k)

               elseif( ideck(k) .eq. 7 )then
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  i3 = nrr(3,k)
                  i4 = nrr(4,k)
                  i5 = nrr(5,k)
                  i6 = nrr(6,k)
                  rlink(k) = y(i1)*y(i2)*sig(k)
                  erate = rlink(k)*qval(k)*cergs

                  in(k) = ideck(k)
                  x1(1,k) = nn(i2)
                  y1(1,k) = nz(i2)
                  x1(2,k) = nn(i6)
                  y1(2,k) = nz(i6)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)
                  
                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p3e9.2,18x,1p2e10.3)') 
     1                 k,ideck(k),rlink(k),
     1                 xid(i1), xid(i2), " --> ", xid(i3),xid(i4),
     1                 xid(i5),xid(i6),
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),y(i1),y(i2),erate,qval(k)
               elseif( ideck(k) .eq. 8 )then
                  i1 = nrr(1,k)
                  i2 = nrr(2,k)
                  i3 = nrr(3,k)
                  i4 = nrr(4,k)
                  i5 = nrr(5,k)
                  if( i5 .eq. 0 )then
                     rlink(k) = y(i1)**3*sig(k)/3.0d0
                     erate = rlink(k)*qval(k)*cergs

                     in(k) = ideck(k)
                     x1(1,k) = nn(i3)
                     y1(1,k) = nz(i3)
                     x1(2,k) = nn(i4)
                     y1(2,k) = nz(i4)
                     ltyp(k) = 10 + log10( rlink(k)/rlmax )
                     ltyp(k) = max0(ltyp(k),0)

                     write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                    1p4e9.2,9x,1p2e10.3)') 
     1                    k,ideck(k),rlink(k),
     1                    xid(i1), xid(i2), xid(i3), " --> ",xid(i4),
     1                    blank,blank,rlkh(k),rnr(k),rvw(k),
     1                    sig(k),y(i1),y(i2),y(i3),erate,qval(k)
                  else
                     rlink(k) = y(i1)**3*sig(k)/3.0d0
                     erate = rlink(k)*qval(k)*cergs

                     in(k) = ideck(k)
                     x1(1,k) = nn(i3)
                     y1(1,k) = nz(i3)
                     x1(2,k) = nn(i5)
                     y1(2,k) = nz(i5)
                     ltyp(k) = 10 + log10( rlink(k)/rlmax )
                     ltyp(k) = max0(ltyp(k),0)

                     write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                    1p4e9.2,9x,1p2e10.3)') 
     1                    k,ideck(k),rlink(k),
     1                    xid(i1), xid(i2), xid(i3), " --> ",xid(i4),
     1                    xid(i5),blank,
     1                    rlkh(k),rnr(k),rvw(k),
     1                    sig(k),y(i1),y(i2),y(i3),erate,qval(k)
                  endif
               else
                  stop'select'
               endif
            endif
         enddo

         if( iabflag .ne. 0 )then
            write(*,*)' nucleon fractions'
            do k = 1, ((netsize+1)/5)*5, 5
               write(*,'(5(a5,1pe11.3))')xid(k),x(k),xid(k+1),x(k+1),
     1              xid(k+2),x(k+2),xid(k+3),x(k+3),xid(k+4),x(k+4)
            enddo
            x(netsize+1) = y(netsize)
            write(*,'(5(a5,1pe11.3))')(xid(k),x(k),
     1           k= ((netsize+1)/5)*5+1,netsize+1 )
         endif

         if( jnb .ne. 0 )then
c..   look for jnb rates
            write(*,*)'Looking for JNB rates in select.f'
            write(*,'(7a5,3a12)')'k','i1','i2','i3','i4','i5','rlkh',
     1           'sig','signue'
c..   p+p ---> d
            do k = 1, 13
               jnbfkt(k) = 0.0d0
            enddo
            do k = k1deck(1),k2deck(1)
               i1 = nrr(1,k)
               i2 = nrr(2,k)
               if( lz(i1) .eq. 7 .and. ln(i1) .eq. 6 )then
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),' ',' ',rlkh(k),
     2                 sig(k),signue(k)
               elseif( lz(i1) .eq. 8 .and. ln(i1) .eq. 7 )then
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),' ',' ',rlkh(k),
     2                 sig(k),signue(k)
               elseif( lz(i1) .eq. 9 .and. ln(i1) .eq. 8 )then
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),' ',' ',rlkh(k),
     2                 sig(k),signue(k)
               elseif( lz(i1) .eq. 9 .and. ln(i1) .eq. 9 )then
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),' ',' ',rlkh(k),
     2                 sig(k),signue(k)

               endif
            enddo

            write(*,'(7a5,3a12)')'k','i1','i2','i3','i4','i5','rlkh',
     1           'sig','sum of sig','signue'
            do k = k1deck(4),k2deck(4)
               i1 = nrr(1,k)
               i2 = nrr(2,k)
               i3 = nrr(3,k)
               if( lz(i1) .eq. 1 .and. ln(i1) .eq. 0 .and.
     1              lz(i2) .eq. 1 .and. ln(i2) .eq. 0  .and.
     2              lz(i3) .eq. 1 .and. ln(i3) .eq. 1 )then
                  jnbfkt(1) = jnbfkt(1) + sig(k)
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),xid(i3),' ',rlkh(k),
     2                 sig(k),jnbfkt(1),signue(k)
               elseif( lz(i1) .eq. 2 .and. ln(i1) .eq. 1 .and.
     1                 lz(i2) .eq. 2 .and. ln(i2) .eq. 2  .and.
     2                 lz(i3) .eq. 4 .and. ln(i3) .eq. 3 )then
                  jnbfkt(3) = jnbfkt(3) + sig(k)
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),xid(i3),' ',rlkh(k),
     2                 sig(k),jnbfkt(3),signue(k)
               elseif( lz(i1) .eq. 1 .and. ln(i1) .eq. 0 .and.
     1                 lz(i2) .eq. 6 .and. ln(i2) .eq. 6  .and.
     2                 lz(i3) .eq. 7 .and. ln(i3) .eq. 6 )then
                  jnbfkt(4) = jnbfkt(4) + sig(k)
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),xid(i3),' ',rlkh(k),
     2                 sig(k),jnbfkt(4),signue(k)
               elseif( lz(i1) .eq. 1 .and. ln(i1) .eq. 0 .and.
     1                 lz(i2) .eq. 6 .and. ln(i2) .eq. 7  .and.
     2                 lz(i3) .eq. 7 .and. ln(i3) .eq. 7 )then
                  jnbfkt(5) = jnbfkt(5) + sig(k)
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),xid(i3),' ',rlkh(k),
     2                 sig(k),jnbfkt(5),signue(k)
               elseif( lz(i1) .eq. 1 .and. ln(i1) .eq. 0 .and.
     1                 lz(i2) .eq. 7 .and. ln(i2) .eq. 7  .and.
     2                 lz(i3) .eq. 8 .and. ln(i3) .eq. 7 )then
                  jnbfkt(6) = jnbfkt(6) + sig(k)
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),xid(i3),' ',rlkh(k),
     2                 sig(k),jnbfkt(6),signue(k)
               elseif( lz(i1) .eq. 1 .and. ln(i1) .eq. 0 .and.
     1                 lz(i2) .eq. 8 .and. ln(i2) .eq. 8  .and.
     2                 lz(i3) .eq. 9 .and. ln(i3) .eq. 8 )then
                  jnbfkt(7) = jnbfkt(7) + sig(k)
                  write(*,'(i5,4a5,6x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),xid(i3),' ',rlkh(k),
     2                 sig(k),jnbfkt(7),signue(k)
               endif
            enddo
            do k = k1deck(6),k2deck(6)
               i1 = nrr(1,k)
               i2 = nrr(2,k)
               i3 = nrr(3,k)
               i4 = nrr(4,k)
               i5 = nrr(5,k)
               if( lz(i1) .eq. 2 .and. ln(i1) .eq. 1 .and.
     1              lz(i2) .eq. 2 .and. ln(i2) .eq. 1   )then
                  jnbfkt(2) = jnbfkt(2) + sig(k)         
                  write(*,'(i5,5a5,1x,a4,1p3e12.3)')
     1                 k,xid(i1),xid(i2),xid(i3),xid(i4),
     2                 xid(i5),rlkh(k), sig(k),jnbfkt(2),signue(k)
               endif
            enddo


c..   hardwired cf88 rates
            write(*,'(/a30)')'CF88 rates, unscreened'
            t9   = tt(nc,kk)*1.0d-9
            t913 = t9**(1.0d0/3.0d0)
            t923 = t913**2
            t943 = t923**2
            t953 = t943*t913
            t976 = t9**(7.0d0/6.0d0)
            t965 = t9**(6.0d0/5.0d0)
            t932 = t9**1.5d0
            
            T9A = T9/(1.0+4.95e-02*T9)
            T9A13 = T9A**(1.0/3.0)
            T9A56 = T9A**(5.0/6.0)

            write(*,'(a20,8a16)')'reaction','NA<sig v>',
     1           '*rho','/1+delta01'
            cf88 = 4.01e-15/T923*exp(-3.380/T913)*
     1           (1.0+0.123*T913+1.09*T923+0.938*T9)
            write(*,'(a20,1p8e16.4)') 'H1(P,E+NU)H2',
     1           cf88,cf88*rhoz(kk),cf88*rhoz(kk)*0.5d0

            cf88 = 1.36e-20/T976*exp(-3.380/T913)
     |                 *(1.0-0.729*T913+9.82*T923)
            write(*,'(a20,1p8e16.4)')'H1(E-P,NU)H2' ,
     1           cf88,cf88*rhoz(kk),cf88*rhoz(kk)*0.5d0

            cf88 = 6.04e+10/T923*exp(-12.276/T913)*
     1           (1.0+0.034*T913-0.522*T923-0.124*T9+0.353*T943
     2           +0.213*T953)
            write(*,'(a20,1p8e16.4)') 'HE3(HE3,2P)HE4',
     1           cf88,cf88*rhoz(kk),cf88*rhoz(kk)*0.5d0


            cf88 = 5.61e+06*T9A56/T932*exp(-12.826/T9A13)
            write(*,'(a20,1p8e16.4)') 'HE4(HE3,G)BE7',
     1           cf88,cf88*rhoz(kk)


            cf88 = 2.04e+07/T923*exp(-13.690/T913-(T9/1.500)**2)*
     1           (1.00+0.030*T913+1.19*T923+0.254*T9+2.06*T943
     2           +1.12*T953)+1.08e+05/T932*exp(-4.925/T9)
     3           +2.15e+05/T932*exp(-18.179/T9)
            write(*,'(a20,1p8e16.4)')'C12(P,G)N13' ,
     1           cf88,cf88*rhoz(kk)

            cf88 = 8.01e+07/T923*exp(-13.717/T913-(T9/2.000)**2)*
     1           (1.0+0.030*T913+0.958*T923+0.204*T9+1.39*T943
     2           +0.753*T953)+1.21e+06/T965*exp(-5.701/T9)
            write(*,'(a20,1p8e16.4)')'C13(P,G)N14' ,
     1           cf88,cf88*rhoz(kk)


            cf88 = 4.90e+07/T923*exp(-15.228/T913-(T9/3.294)**2)*
     1           (1.0+0.027*T913-0.778*T923-0.149*T9+0.261*T943
     2           +0.127*T953)+2.37e+03/T932*exp(-3.011/T9)
     3           +2.19e+04*exp(-12.530/T9)
            write(*,'(a20,1p8e16.4)') 'N14(P,G)O15',
     1           cf88,cf88*rhoz(kk)


            cf88 =  1.50e+08/(T923*(1.0+2.13*(1.0-exp(-0.728*T923))))
     |      *exp(-16.692/T913)
            write(*,'(a20,1p8e16.4)')'O16(P,G)F17' ,
     1           cf88,cf88*rhoz(kk)


         endif
c         stop'3'
c--------------------------------------------------------------------
      elseif(  tt(nc,kk) .gt. tburnd*1.0d9 )then


c..lodeck................................................................
c..   get flows for reaction links
         do k = l1deck(1),l2deck(1)
            i1 = lonrr(1,k)
            rlink(k) = ylo(i1)*sig(k)
         enddo
         do k = l1deck(2),l2deck(2)
            i1 = lonrr(1,k)
            ysig = ylo(i1)*sig(k)
            rlink(k) = ysig
         enddo
         do k = l1deck(3),l2deck(3)
            i1 = lonrr(1,k)
            ysig = ylo(i1)*sig(k)
            rlink(k) = ysig
         enddo
         do k = l1deck(4),l2deck(4)
            i1 = lonrr(1,k)
            i2 = lonrr(2,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)  
            rlink(k) = yysig
         enddo
         do k = l1deck(5),l2deck(5)
            i1 = lonrr(1,k)
            i2 = lonrr(2,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)  
            rlink(k) = yysig
         enddo
         do k = l1deck(6),l2deck(6)
            i1 = lonrr(1,k)
            i2 = lonrr(2,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)   
            rlink(k) = yysig
         enddo
         do k = l1deck(7),l2deck(7)
            i1 = lonrr(1,k)
            i2 = lonrr(2,k)
            yysig = ylo(i1)*ylo(i2)*sig(k)   
            rlink(k) = yysig
         enddo
         do  k = l1deck(8),l2deck(8)
            i1 = lonrr(1,k)
            y3sig = ylo(i1)**3*sig(k)/3.0d0
            rlink(k) = y3sig
         enddo

c..   fastest reaction link
         rlmax = 0.0d0
         maxk = 0
         do k = 1, l2deck(8)
            if( rlink(k) .gt. rlmax )then
               rlmax = rlink(k)
               maxk = k
            endif
         enddo

         if( maxk .le. 0 .or. maxk .gt. l2deck(8) )then
            write(*,*)'select.f: maxk error ',maxk
            stop'select.f'
         endif
         write(*,'(a15,a8,a12)')'fastest rate','index','reaction'
         write(*,'(1pe15.3,2i8,6a5)')rlmax,maxk,lodeck(maxk),
     1        (rname(i,lorr(maxk)),i=1,6)

         write(*,'(/2(a15,1pe11.3),a15,i5/)')'T(K)',tt(nc,kk),
     1        'rho(g/cc)',dd(nc,kk),'zone',kk

c..   edited list of rates
         write(*,'(a4,a2,a10,12x,a8,22x,a6,
     1        6(a7,2x))') 'i ','dk','flow','reaction','source','sig',
     2        'y1','y2'

         do k = 1, l2deck(8)
            in(k)   = 0
            ldot(k) = 1
            if( rlink(k) .gt. epsilon*rlmax )then
               if( lodeck(k) .eq. 1 )then
                  i1 = lonrr(1,k)
                  i2 = lonrr(2,k)
                  rlink(k) = ylo(i1)*sig(k)

                  ldot(k) = 4
                  in(k)   = lodeck(k)
                  x1(1,k) = lon(i1)
                  y1(1,k) = loz(i1)
                  x1(2,k) = lon(i2)
                  y1(2,k) = loz(i2)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)

                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p6e9.2)') 
     1                 k,lodeck(k),rlink(k),blank,
     1                 xidlo(i1), " --> ", xidlo(i2),blank,blank,
     1                 blank,lorlkh(k),rnr(lorr(k)),rvw(lorr(k)),
     1                 sig(k),ylo(i1)

               elseif( lodeck(k) .eq. 2 )then
                  i1 = lonrr(1,k)
                  i2 = lonrr(2,k)
                  i3 = lonrr(3,k)
                  rlink(k) = ylo(i1)*sig(k)
                  in(k) = lodeck(k)
                  x1(1,k) = lon(i1)
                  y1(1,k) = loz(i1)
                  x1(2,k) = lon(i3)
                  y1(2,k) = loz(i3)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)

                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p7e9.2)') 
     1                 k,lodeck(k),rlink(k),blank,
     1                 xidlo(i1), " --> ", xidlo(i2), xidlo(i3),blank,
     1                 blank,lorlkh(k),rnr(lorr(k)),rvw(lorr(k)),
     1                 sig(k),ylo(i1)

               elseif( lodeck(k) .eq. 3 )then
                  i1 = lonrr(1,k)
                  i2 = lonrr(2,k)
                  i3 = lonrr(3,k)
                  i4 = lonrr(4,k)
                  rlink(k) = ylo(i1)*sig(k)

                  in(k) = lodeck(k)
                  x1(1,k) = lon(i1)
                  y1(1,k) = loz(i1)
                  x1(2,k) = lon(i4)
                  y1(2,k) = loz(i4)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)

                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p7e9.2)') 
     1                 k,lodeck(k),rlink(k),blank,
     1                 xidlo(i1), " --> ", xidlo(i2), xidlo(i3), 
     1                 xidlo(i4),blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),ylo(i1)

               elseif( lodeck(k) .eq. 4 )then
                  i1 = lonrr(1,k)
                  i2 = lonrr(2,k)
                  i3 = lonrr(3,k)

                  rlink(k) = ylo(i1)*ylo(i2)*sig(k)

                  in(k) = lodeck(k)
                  x1(1,k) = lon(i2)
                  y1(1,k) = loz(i2)
                  x1(2,k) = lon(i3)
                  y1(2,k) = loz(i3)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)
                  if( rlkh(k).eq.'  ec' .or. rlkh(k).eq.'bet+')then
                     ldot(k) = 4
                  endif

                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p7e9.2)') 
     1                 k,lodeck(k),rlink(k),
     1                 xidlo(i1), xidlo(i2), " --> ", xidlo(i3),
     1                 blank,blank,blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),ylo(i1),ylo(i2)
               elseif( lodeck(k) .eq. 5 )then
                  i1 = lonrr(1,k)
                  i2 = lonrr(2,k)
                  i3 = lonrr(3,k)
                  i4 = lonrr(4,k)
                  rlink(k) = ylo(i1)*ylo(i2)*sig(k)

                  in(k) = lodeck(k)
                  x1(1,k) = lon(i2)
                  y1(1,k) = loz(i2)
                  x1(2,k) = lon(i4)
                  y1(2,k) = loz(i4)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)
                  
                  write(*,'(i4,i2,1pe10.2,7(a5,1x)a4,2a1,
     1                 1p5e9.2)') 
     1                 k,lodeck(k),rlink(k),
     1                 xidlo(i1), xidlo(i2), " --> ", xidlo(i3),
     1                 xidlo(i4),blank,blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),ylo(i1),ylo(i2)
               elseif( lodeck(k) .eq. 6 )then
                  i1 = lonrr(1,k)
                  i2 = lonrr(2,k)
                  i3 = lonrr(3,k)
                  i4 = lonrr(4,k)
                  i5 = lonrr(5,k)
                  rlink(k) = ylo(i1)*ylo(i2)*sig(k)

                  in(k) = lodeck(k)
                  x1(1,k) = lon(i2)
                  y1(1,k) = loz(i2)
                  x1(2,k) = lon(i5)
                  y1(2,k) = loz(i5)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)
                  
                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p5e9.2)') 
     1                 k,lodeck(k),rlink(k),
     1                 xidlo(i1), xidlo(i2), " --> ", xidlo(i3),
     1                 xidlo(i4),xidlo(i5),blank,
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),ylo(i1),ylo(i2)
               elseif( lodeck(k) .eq. 7 )then
                  i1 = lonrr(1,k)
                  i2 = lonrr(2,k)
                  i3 = lonrr(3,k)
                  i4 = lonrr(4,k)
                  i5 = lonrr(5,k)
                  i6 = lonrr(6,k)
                  rlink(k) = ylo(i1)*ylo(i2)*sig(k)

                  in(k) = lodeck(k)
                  x1(1,k) = lon(i2)
                  y1(1,k) = loz(i2)
                  x1(2,k) = lon(i6)
                  y1(2,k) = loz(i6)
                  ltyp(k) = 10 + log10( rlink(k)/rlmax )
                  ltyp(k) = max0(ltyp(k),0)
                  
                  write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                 1p7e9.2)') 
     1                 k,lodeck(k),rlink(k),
     1                 xidlo(i1), xidlo(i2), " --> ", xidlo(i3),
     1                 xidlo(i4),xidlo(i5),xidlo(i6),
     1                 rlkh(k),rnr(k),rvw(k),
     1                 sig(k),ylo(i1),ylo(i2)
               elseif( lodeck(k) .eq. 8 )then
                  i1 = lonrr(1,k)
                  i2 = lonrr(2,k)
                  i3 = lonrr(3,k)
                  i4 = lonrr(4,k)
                  i5 = lonrr(5,k)
                  if( i5 .eq. 0 )then
                     rlink(k) = ylo(i1)**3*sig(k)/3.0d0

                     in(k) = lodeck(k)
                     x1(1,k) = lon(i3)
                     y1(1,k) = loz(i3)
                     x1(2,k) = lon(i4)
                     y1(2,k) = loz(i4)
                     ltyp(k) = 10 + log10( rlink(k)/rlmax )
                     ltyp(k) = max0(ltyp(k),0)

                     write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                    1p5e9.2)') 
     1                    k,lodeck(k),rlink(k),
     1                    xidlo(i1), xidlo(i2), xidlo(i3), " --> ",
     1                    xidlo(i4),
     1                    blank,blank,rlkh(k),rnr(k),rvw(k),
     1                    sig(k),ylo(i1),ylo(i2),ylo(i3)
                  else
                     rlink(k) = ylo(i1)**3*sig(k)/3.0d0

                     in(k) = lodeck(k)
                     x1(1,k) = lon(i3)
                     y1(1,k) = loz(i3)
                     x1(2,k) = lon(i5)
                     y1(2,k) = loz(i5)
                     ltyp(k) = 10 + log10( rlink(k)/rlmax )
                     ltyp(k) = max0(ltyp(k),0)

                     write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1                    1p6e9.2)') 
     1                    k,lodeck(k),rlink(k),
     1                    xidlo(i1), xidlo(i2), xidlo(i3), " --> ",
     1                    xidlo(i4),xidlo(i5),blank,
     1                    rlkh(k),rnr(k),rvw(k),
     1                    sig(k),ylo(i1),ylo(i2),ylo(i3)
                  endif
               else
                  stop'select'
               endif
            endif
         enddo

         write(*,*)'nuclei in low Temperature network'
c         write(*,'(10(2x,i3,a5))')(i,xidlo(i),i=1,lodim)
         if( iabflag .ne. 0 )then
            write(*,*)' nucleon fractions'
            do k = 1, ndim/5*5, 5
               write(*,'(5(a5,1pe11.3))')xid(k),x(k),
     1              xid(k+1),x(k+1),
     1              xid(k+2),x(k+2),xid(k+3),x(k+3),
     1              xid(k+4),x(k+4)
            enddo
            do k = ndim/5*5+1,ndim
               if( k .ne. ndim )then
                  write(*,'(5(a5,1pe11.3))')xid(k),x(k)
               else
                  write(*,'(5(a5,1pe11.3))')xid(k),y(k)
               endif
            enddo
         endif

      else

c..   get flows for decays only (lowest T)..............................
         do k = k1deck(1), k2deck(1)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            rlink(k) = y(i1) * sig(k)
         enddo

c..   fastest reaction link
         rlmax = 0.0d0
         do k = k1deck(1),k2deck(1)
            if( rlink(k) .gt. rlmax )then
               rlmax = rlink(k)
               maxk = k
            endif
         enddo
         write(*,*)'decays only for this low temperature'

         write(*,'(a15,2a8,a12)')'fastest rate','index','deck',
     1        'reaction'

         write(*,'(1pe15.3,2i8,6a5)')rlmax,maxk,ideck(maxk),
     1        (rname(i,maxk),i=1,6)

         write(*,'(/2(a15,1pe11.3),a15,i5/)')'T(K)',tt(nc,kk),
     1        'rho(g/cc)',dd(nc,kk),'zone',kk

c..   edited list of rates
         write(*,'(a4,a2,a10,12x,a8,22x,a6,
     1        6(a7,2x))') 'i ','dk','flow','reaction','source','sig',
     2        'y1','y2'

         do k = k1deck(1),k2deck(1)
            in(k)   = 0
            ldot(k) = 1
            if( rlink(k) .gt. epsilon*rlmax )then
               i1 = nrr(1,k)
               i2 = nrr(2,k)
               rlink(k) = y(i1)*sig(k)
c..solid lines for decays only zones
               ldot(k) = 1
               in(k)   = ideck(k)
               x1(1,k) = nn(i1)
               y1(1,k) = nz(i1)
               x1(2,k) = nn(i2)
               y1(2,k) = nz(i2)
               ltyp(k) = 10 + log10( rlink(k)/rlmax )
               ltyp(k) = max0(ltyp(k),0)
               write(*,'(i4,i2,1pe10.2,7(a5,1x),a4,2a1,
     1              1p6e9.2)') 
     1              k,ideck(k),rlink(k),blank,
     1              xid(i1), " --> ", xid(i2),blank,blank,blank,
     1              rlkh(k),rnr(k),rvw(k),
     1              sig(k),y(i1),y(i2)
            endif

         enddo

      endif

      end
