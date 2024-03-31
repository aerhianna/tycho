      subroutine rate(t90,rho,dtstar,kzone,jnb)

c..   fkt rate formula with T9, V derivatives (wda 5/9/99)
c..   weak rates added (pay 11/14/2004)
c..   weak screening (wda 11/14/04)
c..   statistical inverse rates (wda 11/14/04)
c..   jnb (Bahcall) correction factor for solar models 
c..   (177 nuclei network only)

      implicit none

      include 'dimenfile'
      include 'crate'
      include 'comcsolve'
      include 'ceoset'
      include 'cburn'
      include 'cconst'

      integer*4 tdex, rhodex
      real*8 t9array(tsize), rhoarray(rhosize)

      real*8 dtstar, linklim
      real*8 t9, rho, t913, t953, t9ln, expon, exponu
      real*8 t92, t923, t943, factor, third
      real*8 yee, t90, t9lim, logrho

      integer*4 i1, i2, i3, i4, i5, i6

c..   irev  is the index of the reverse rate
c..   iline is the ordering number for subrates of a given reaction
c..   (this includes nonresonant and resonant rates)
c..   isum  is a dummy variable for counting
c      integer*4 isum
      real*8 fak

c..t9lim is the same value as tburnlo parameter in cdeuter
      parameter( third = 1.0d0/3.0d0, t9lim = 0.01d0 )

      integer*4 k, k1, k2, j, kzone, jnb

      integer*4 i

      data t9array/0.01d0,0.10d0,0.20d0,0.40d0,0.70d0,
     1     1.0d0,1.5d0,2.0d0,3.0d0,5.0d0,10.0d0,30.0d0,100.0d0/
      data rhoarray/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,
     1     7.0d0,8.0d0,9.0d0,10.0d0,11.0d0/
      data linklim/1.0d14/
c----------------------------------------------------
c..   avoid bad fkt extrapolations to very low temperature by setting
c..   a basement value 
      if( t90 .gt. t9lim )then
         t9 = t90
      else
         t9 = t9lim
      endif
      
c..   electron mole number Ye = Ne/rho*avagadro
      yee = 0.0d0
      do j = 1, nnuc
         yee = yee + y(j) * dble( nz(j) )
      enddo

      t913 = t9**third
      t953 = t9 * t913 *t913
      t9ln = dlog(t9)

      t92  = t9*t9
      t923 = t913*t913
      t943 = t913*t9

c..   calculate nuclear partition functions
      logrho = dlog10(rho)
      if(logrho .le. 1.0d0) logrho = 1.0d0

      call locate(t9array,tsize,t9,tdex)
      call locate(rhoarray,rhosize,logrho,rhodex)

      logrho = dlog10(rho)

      call partfun(t9)
      
c..   define sig, sigt (sigv is calculated from sig below)
      do k = 1,ireac
         if( irev(k) .ge. 0 )then
c..   avoid reverse ffnu rates
            expon = rcoef(1,k) + rcoef(2,k)/t9 + rcoef(3,k)/t913
     1           + rcoef(4,k)*t913 + rcoef(5,k)*t9
     2           + rcoef(6,k)*t953 + rcoef(7,k)*t9ln

            if( expon .gt. -300.0d0 )then
               sig(k) = exp(expon)
               factor = -rcoef(2,k)/t92 -third*rcoef(3,k)/t943
     1              +third*rcoef(4,k)/t923 + rcoef(5,k)
     2              +5.0d0*third*rcoef(6,k)*t923 + rcoef(7,k)/t9
               sigt(k) = sig(k) * factor
            else
               sig(k)  = 0.0d0
               sigt(k) = 0.0d0
            endif
         endif
      enddo
      
c..   initialize
      do k = 1, ireac
         signue(k) = 0.0d0
      enddo
      
c..   calculate density factors, deck by deck
c..   deck 1: electron capture density*Ye factor
      k1 = k1deck(1)
      k2 = k2deck(1)
      do k = k1, k2
c         if( rlkh(k) .eq. ' ffn' )then
         if( iffn(k) .ne. 0 )then
            if(t9 .gt. t9lim .and. logrho .gt. 1.0d0)then

               call weakintp(t9,logrho,iffn(k),expon,exponu,
     1              tdex,rhodex,k)

               if( expon .gt. -300.0d0 )then
                  sig(k)  = 10.0d0**(expon)
                  sigt(k) = 0.0d0
               else
                  sig(k)  = 0.0d0
                  sigt(k) = 0.0d0
               endif
               if( exponu .gt. -300.0d0 )then
                  signue(k) = 10.0d0**(exponu)
               else
                  signue(k) = 0.0d0
               endif
            else
               if( irev(k) .lt. 0 )then
c..   ffnu reverse rates set to zero at low temperature
                  sigt(k)   = 0.0d0
                  sig(k)    = 0.0d0
                  signue(k) = 0.0d0
               else
c..   only lab forward rates
c..   pass through
               endif
            endif
         else
c..   no fuller,fowler,newman
            signue(k) = 0.0d0
         endif
      enddo
      
      do k = k1, k2
c         if( rlkh(k) .eq. '  ec' .or. rlkh(k) .eq. ' bec'
c     1       .or. rlkh(k) .eq. ' ecw')then
         if(ec(k) .eq. 1)then
            sig(k)  = yee * rho * sig(k)
            sigt(k) = yee * rho * sigt(k)
            sigv(k) = -sig(k) * rho
         endif
      enddo

c..   decks 2 and 3 have no density dependence
      k1 = k1deck(2)
      k2 = k2deck(3)
      do k = k1, k2
         sigv(k) = 0.0d0
      enddo
c..   decks 4 to 7 have binary interaction in entrance channel
      k1 = k1deck(4)
      k2 = k2deck(7)
      do k = k1, k2
         if( nrr(1,k) .eq. nrr(2,k) )then
c..   pairs reduced for identical particles to avoid double counting
c..   sigv is constructed from sig
            sig(k)  = 0.5d0 * sig(k)
            sigt(k) = 0.5d0 * sigt(k)
         endif
c         if(rname(1,k) .eq. '  c12' .and. rname(2,k) .eq. '  he4')then
c            sig(k) = 0.5d0 * sig(k)
c            sigt(k) = 0.5d0 * sigt(k)
c         endif
c         if(rname(2,k) .eq. '  c12' .and. rname(1,k) .eq. '  he4')then
c            sig(k) = 0.5d0 * sig(k)
c            sigt(k) = 0.5d0 * sigt(k)
c         endif
c..   this is another power of rho
         sigv(k)  = -rho * rho * sig(k)
         sig(k)   =  rho * sig(k)
         sigt(k)  =  rho * sigt(k)

      enddo
      
c..   decks 8 and 9 have triple (tertiary collision) in entrance channel
      k1 = k1deck(8)
      k2 = k2deck(9)
      do k = k1, k2
         if( nrr(1,k) .eq. nrr(2,k) .and. nrr(2,k) .eq. nrr(3,k))then
c..   waf, grc, baz (1967) ann rev a & a 5, 558 definition
c     has an extra factor of 1/2 for 3-alpha
c..   from 3 He4 destroyed per reaction, n/3! triples --> 1/2
c..   d/dV gives a compensating factor of 2
c..   factor of 1/3 moved from rhside and jacob to here 8/20/00
            sigv(k)  = -rho * rho * rho * sig(k) /3.0d0
            sig(k)   =        rho * rho * sig(k) /6.0d0
            sigt(k)  =        rho * rho * sigt(k)/6.0d0
         elseif(nrr(1,k) .eq. nrr(2,k) .or. nrr(1,k) .eq. nrr(3,k)
     1          .or. nrr(2,k) .eq. nrr(3,k))then
            sigv(k)  = -rho * rho * rho * sig(k) 
            sig(k)   =        rho * rho * sig(k) /2.0d0
            sigt(k)  =        rho * rho * sigt(k)/2.0d0
         else
            sigv(k) = -rho * rho * rho * sig(k) * 2.0d0
            sig(k)  =        rho * rho * sig(k)
            sigt(k) =        rho * rho * sigt(k)
         endif
      enddo

c..   deck 10 has quadruple (quaternary collision) in entrance channel
      k1 = k1deck(10)
      k2 = k2deck(10)
      do k = k1, k2
         if( nrr(1,k) .eq. nrr(2,k) .and. nrr(2,k) .eq. nrr(3,k)
     1       .and. nrr(3,k) .eq. nrr(4,k))then
c..   waf, grc, baz (1967) ann rev a & a 5, 558 definition
c     has an extra factor of 1/2 for 3-alpha
c..   from 3 He4 destroyed per reaction, n/3! triples --> 1/2
c..   d/dV gives a compensating factor of 2
c..   factor of 1/3 moved from rhside and jacob to here 8/20/00
            sigv(k)  = -rho * rho * rho * rho * sig(k) /8.0d0
            sig(k)   =        rho * rho * rho * sig(k) /24.0d0
            sigt(k)  =        rho * rho * rho * sigt(k)/24.0d0
         elseif(nrr(1,k) .eq. nrr(2,k) .and. nrr(2,k) .eq. nrr(3,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /2.0d0
            sig(k)  =        rho * rho * rho * sig(k) /6.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/6.0d0
         elseif(nrr(1,k) .eq. nrr(2,k) .and. nrr(2,k) .eq. nrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /2.0d0
            sig(k)  =        rho * rho * rho * sig(k) /6.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/6.0d0
         elseif(nrr(2,k) .eq. nrr(3,k) .and. nrr(3,k) .eq. nrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /2.0d0
            sig(k)  =        rho * rho * rho * sig(k) /6.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/6.0d0
         elseif(nrr(1,k) .eq. nrr(3,k) .and. nrr(3,k) .eq. nrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) /2.0d0
            sig(k)  =        rho * rho * rho * sig(k) /6.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/6.0d0
         elseif(nrr(1,k) .eq. nrr(2,k) .and. nrr(3,k) .eq. nrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) *3.0d0/4.0d0
            sig(k)  =        rho * rho * rho * sig(k) /4.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/4.0d0
         elseif(nrr(1,k) .eq. nrr(2,k) .or. nrr(1,k) .eq. nrr(3,k)
     1          .or. nrr(1,k) .eq.nrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) *1.5d0
            sig(k)  =        rho * rho * rho * sig(k) /2.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/2.0d0
         elseif(nrr(2,k) .eq. nrr(3,k) .or. nrr(2,k) .eq. nrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) *1.5d0
            sig(k)  =        rho * rho * rho * sig(k) /2.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/2.0d0
         elseif(nrr(3,k) .eq. nrr(4,k))then
            sigv(k) = -rho * rho * rho * rho * sig(k) *1.5d0
            sig(k)  =        rho * rho * rho * sig(k) /2.0d0
            sigt(k) =        rho * rho * rho * sigt(k)/2.0d0
         else
            sigv(k) = -rho * rho * rho * rho * sig(k) * 3.0d0
            sig(k)  =        rho * rho * rho * sig(k)
            sigt(k) =        rho * rho * rho * sigt(k)
         endif
      enddo
      
c..   deck 11 has no density dependence
      k1 = k1deck(11)
      k2 = k2deck(11)
      do k = k1, k2
         sigv(k) = 0.0d0
      enddo

c..   inverses not revised for screening in nse? ?????????????
c..
c      write(*,*)'check inverses'
c      do k= k1deck(6),k2deck(8)
c         if( irev(k) .gt. 1 )then
c            j = irev(k)
c            write(*,'(2i6,7a6,5x,a6,2a2,i6,2i6,7a6,5x,a6,2a2,i6)')
c     1           k,ideck(k),
c     1           (rname(i,k),i=1,7),
c     1           rlkh(k),rvw(k),rnr(k),iline(k),
c     2           j,ideck(j),
c     1           (rname(i,j),i=1,7),
c     1           rlkh(j),rvw(j),rnr(j),iline(j)
c         else
c            write(*,'(2i6,7a6,5x,a6,2a2,2i6)')k,ideck(k),
c     1        (rname(i,k),i=1,7),
c     1           rlkh(k),rvw(k),rnr(k),irev(k),iline(k)
c         endif
c      enddo

c..   check for forward rates, and add corresponding reverse rate
c..   deck 1:  i ---> j
c      write(*,*)'deck 1  i ---> j, no forward rates?'
c      do k= k1deck(1),k2deck(1)
c         if( irev(k) .eq. 0 )then
c            write(*,'(a30,2i5,2a6,5x,a6,2a2,2i5)')
c     1           'deck 1  i ---> j, lonely rates?',
c     1           k,ideck(k),
c     1           (rname(i,k),i=1,2),
c     1           rlkh(k),rvw(k),rnr(k),irev(k),iline(k)
c         endif
c      enddo
c deck 1: weak decays, i ---> j
c      do k= k1deck(1),k2deck(1)
c         if ( irev(k) .gt. 0 )then
c..forward rate, overwrite fkt value with detailed balance result
c            i1 = nrr(1,k)
c            i2 = nrr(2,k)
c            sig(irev(k)) = sig(irev(k))*(pf(i2)/pf(i1))
c            fak = pf(i2)/pf(i1)
c           fak =1.d0/( pf(i2)/(pf(i1))
c     1           *( anuc(i2)/anuc(i1))**1.5d0
c     2           *exp( -11.605d0*qval(k)/t9 ) )
c            sig(irev(k)) = fak * sig(irev(k))
c            sigt(irev(k)) = fak * sigt(irev(k))
c            sigv(irev(k)) = fak * sigv(irev(k))
c..   overwrite fkt value
c            sig(irev(k)) = fak * sig(k)
c            sigt(irev(k)) = fak * sigt(k)
c            sigv(irev(k)) = fak * sigv(k)
c         endif
c         write(*,*)k,irev(k),sig(k),sig(irev(k)),qval(k),fak
c      enddo
      
c.. deck 2: dissociations, i ---> j + k,
c      write(*,*)'deck 2  i ---> j + k, no forward rates?'
c      do k= k1deck(2),k2deck(2)
c            i1 = nrr(1,k)
c            i2 = nrr(2,k)
c            i3 = nrr(3,k)
c..   forward rates k
c..   irev = 0 for no reverse rate (eg. weak decays)
c..   irev < 0 for reverse rate k 
c..   (to be dealt with from forward rate elsewhere)
c            if( irev(k) .gt. 0 )then
c..   k is forward rate, irev(k) is reverse rate
c               fak =1.0d0 * (pf(i3)*pf(i2)/pf(i1) / rho * 9.8678d9
c     1              *( anuc(i2)*anuc(i3)/anuc(i1))**1.5d0*t9**1.5d0
c     2              *dexp( -11.605d0*qval(k)/t9 ))
c               fak =pf(i1)*pf(i2)/pf(i3) 
c     1              *( anuc(i1)*anuc(i2)/anuc(i3) )**1.5d0
c     2              *exp( -11.605d0*qval(k)/t9 )
c..   overwrite fkt value
c               sig(irev(k)) = fak * sig(k)
c               sigt(irev(k)) = fak * sigt(k)
c               sigv(irev(k)) = fak * sigv(k)
c               if( nrr(1,k) .eq. nrr(2,k) )then
c                  sig(irev(k)) = 2.0d0*sig(irev(k))
c                  sigt(irev(k)) = 2.0d0*sigt(irev(k))
c                  sigv(irev(k)) = 1.0d0*sigv(irev(k))
c               endif
c               if(irev(k) .eq. 810)then
c                  write(*,*)k,irev(k),qval(k),qval(irev(k)),
c     1                      dexp(-11.605d0*qval(k)/t9) 
c               endif
c            elseif( irev(k) .lt. 0 )then
c               write(*,*)k,irev(k),' irev ERROR',
c     1              ideck(k),ideck(-irev(k))
c               write(*,'(2i5,3a6,5x,a5,a1,1p3e12.3,3i5,1p6e11.3)')
c     1              k,ideck(k),
c     1              (rname(i,k),i=1,3),rlkh(k),rvw(k),sig(k),
c     2              sig(irev(k))*pf(i1)*pf(i2)/pf(i3) ,
c     2              fak,k,irev(k),iline(k)
c               stop'rate 1'
c            endif
c         if( irev(k) .gt. 0 )then
c            write(*,'(a20,2i5,3a6,5x,a6,2a2,2i5)')
c     1           'deck 2  i ---> j + k, forward rates?',
c     1           k,ideck(k),
c     1           (rname(i,k),i=1,3),
c     1           rlkh(k),rvw(k),rnr(k),irev(k),iline(k)
c         endif
c                write(*,*)k,sig(k),sig(irev(k)),k,irev(k),fak,qval(k)
c      enddo
     
ccc
c..   deck 3: i ---> j + k + l
c      write(*,*)'deck 3: i ---> j + k + l, no forward rates?'
c      do k= k1deck(3),k2deck(3)
c         if( irev(k) .gt. 0 )then
c      do k= k1deck(3),k2deck(3)
c            i1 = nrr(1,k)
c            i2 = nrr(2,k)
c            i3 = nrr(3,k)
c            i4 = nrr(4,k)
c..   forward rates k
c..   irev = 0 for no reverse rate (eg. weak decays)
c..   irev < 0 for reverse rate k 
c..   (to be dealt with from forward rate elsewhere)
c            if( irev(k) .gt. 0 )then
c..   k is forward rate, irev(k) is reverse rate
c               fak =pf(i4)*pf(i3)*pf(i2)/pf(i1) * (9.8678d9/rho)**2.0d0
c     1          *( anuc(i2)*anuc(i3)*anuc(i4)/anuc(i1))**1.5d0*t9**3.0d0
c     2              *dexp( -11.605d0*qval(k)/t9 )
c               fak =pf(i1)*pf(i2)/pf(i3) 
c     1              *( anuc(i1)*anuc(i2)/anuc(i3) )**1.5d0
c     2              *exp( -11.605d0*qval(k)/t9 )
c..   overwrite fkt value
c               sig(irev(k)) = fak * sig(k)
c               sigt(irev(k)) = fak * sigt(k)
c               sigv(irev(k)) = fak * sigv(k)
c               if( nrr(1,k) .eq. nrr(2,k) )then
c                  sig(irev(k)) = 2.0d0*sig(irev(k))
c                  sigt(irev(k)) = 2.0d0*sigt(irev(k))
c                  sigv(irev(k)) = 1.0d0*sigv(irev(k))
c               endif
c               if(irev(k) .eq. 810)then
c                  write(*,*)k,irev(k),qval(k),qval(irev(k)),
c     1                      dexp(-11.605d0*qval(k)/t9) 
c               endif
c            elseif( irev(k) .lt. 0 )then
c               write(*,*)k,irev(k),' irev ERROR',
c     1              ideck(k),ideck(-irev(k))
c               write(*,'(2i5,3a6,5x,a5,a1,1p3e12.3,3i5,1p6e11.3)')
c     1              k,ideck(k),
c     1              (rname(i,k),i=1,3),rlkh(k),rvw(k),sig(k),
c     2              sig(irev(k))*pf(i1)*pf(i2)/pf(i3) ,
c     2              fak,k,irev(k),iline(k)
c               stop'rate 1'
c            endif
c         if( irev(k) .gt. 0 )then
c            write(*,'(a20,2i5,3a6,5x,a6,2a2,2i5)')
c     1           'deck 2  i ---> j + k, forward rates?',
c     1           k,ideck(k),
c     1           (rname(i,k),i=1,3),
c     1           rlkh(k),rvw(k),rnr(k),irev(k),iline(k)
c         endif
c                write(*,*)k,sig(k),sig(irev(k)),k,irev(k),fak,qval(k)
c      enddo
c            write(*,'(a20,2i5,4a6,5x,a1)')
c     1           'deck 3: i ---> j + k + l, forward rates?',
c     1           k,ideck(k),(rname(i,k),i=1,4),rvw(k)
c         endif
c      enddo
c.. deck 4: captures, i + j ---> k
      
      do k= k1deck(4),k2deck(4)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
c..   forward rates k
c..   irev = 0 for no reverse rate (eg. weak decays)
c..   irev < 0 for reverse rate k 
c..   (to be dealt with from forward rate elsewhere)
            if( irev(k) .gt. 0 )then
c..   k is forward rate, irev(k) is reverse rate
               fak =9.8678d-9 * pf(i3)/pf(i1)/pf(i2) * rho
     1              *( anuc(i3)/anuc(i2)/anuc(i3))**1.5d0
     2              *dexp( -11.605d0*qval(k)/t9 )/t9**1.5
c               fak =pf(i1)*pf(i2)/pf(i3) 
c     1              *( anuc(i1)*anuc(i2)/anuc(i3) )**1.5d0
c     2              *exp( -11.605d0*qval(k)/t9 )
c..   overwrite fkt value
               sig(irev(k)) = fak * sig(k)
               sigt(irev(k)) = fak * sigt(k)
               sigv(irev(k)) = fak * sigv(k)
               if( nrr(1,k) .eq. nrr(2,k) )then
                  sig(irev(k)) = 2.0d0*sig(irev(k))
                  sigt(irev(k)) = 2.0d0*sigt(irev(k))
                  sigv(irev(k)) = 1.0d0*sigv(irev(k))
               endif
c               if(irev(k) .eq. 810)then
c                  write(*,*)k,irev(k),qval(k),qval(irev(k)),
c     1                      dexp(-11.605d0*qval(k)/t9) 
c               endif
c            elseif( irev(k) .lt. 0 )then
c               write(*,*)k,irev(k),' irev ERROR',
c     1              ideck(k),ideck(-irev(k))
c               write(*,'(2i5,3a6,5x,a5,a1,1p3e12.3,3i5,1p6e11.3)')
c     1              k,ideck(k),
c     1              (rname(i,k),i=1,3),rlkh(k),rvw(k),sig(k),
c     2              sig(irev(k))*pf(i1)*pf(i2)/pf(i3) ,
c     2              fak,k,irev(k),iline(k)
c               stop'rate 1'
            else
c..   for 176 nucleus network, these are weak interactions only
c               write(*,*)k,irev(k),' irev ERROR',
c     1              ideck(k),ideck(-irev(k))
c               write(*,'(2i5,3a6,5x,a5,a1,1p3e12.3,3i5,1p6e11.3)')
c     1              k,ideck(k),
c     1              (rname(i,k),i=1,3),rlkh(k),rvw(k),sig(k),
c     2              sig(irev(k))*pf(i1)*pf(i2)/pf(i3) ,
c     2              fak,k,irev(k),iline(k)
            endif
c                write(*,*)k,sig(k),sig(irev(k)),k,irev(k),fak,qval(k)    
    
      enddo
      
c deck 5: exchange, i + j ---> k + l
      do k= k1deck(5),k2deck(5)
         if ( irev(k) .lt. 0 )then
c..forward rate, overwrite fkt value with detailed balance result
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
c            sig(irev(k)) = sig(irev(k))*(pf(i3)*pf(i4)/pf(i1)/pf(i2))
            fak = pf(i3)*pf(i4)/(pf(i1)*pf(i2))
     1           *( anuc(i4)*anuc(i3)/(anuc(i1)*anuc(i2)))**1.5d0
     2           *exp( -11.605d0*qval(k)/t9 ) 

c..   overwrite fkt value
            sig(irev(k)) = fak * sig(k)
            sigt(irev(k)) = fak * sigt(k)
            sigv(irev(k)) = fak * sigv(k)
            if( nrr(1,k) .eq. nrr(2,k) )then
               sig(irev(k)) = 2.0d0*sig(irev(k))
               sigt(irev(k)) = 2.0d0*sigt(irev(k))
               sigv(irev(k)) = 1.0d0*sigv(irev(k))
            endif
c         elseif( irev(k) .eq. 0 )then
c..   heavy ion reactions have no inverse implemented so irev=0 for them
c            write(*,'(2i5,4a6,5x,a5,a1,1pe12.3,2i5,1pe12.3)')k,
c     1           ideck(k),
c     1           (rname(i,k),i=1,4),rlkh(k),rvw(k),qval(k),irev(k)
c..   for irev <0, the inverses are the forward rates; do nothing
         endif
      enddo

c deck 6: exchange, i + j ---> k + l + m
c..use fkt for now ccccccccccccccccccccccccccccccccccccccccccccc
      do k= k1deck(6),k2deck(6)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            if( irev(k) .gt. 0 )then
c..add partition function ratio for inverse reaction a la friedel
               fak=(pf(i3)*pf(i4)*pf(i5)/pf(i1)/pf(i2))/rho*9.8678d9
     1          *(anuc(i3)*anuc(i4)*(anuc(i5)/anuc(i1)*anuc(i2)))**1.5d0
     2          *exp( -11.605d0*qval(k)/t9 ) * t9**1.5d0
               sig(irev(k)) = fak * sig(k)
               sigt(irev(k)) = fak * sigt(k)
               sigv(irev(k)) = fak * sigv(k)
c            endif
            if( nrr(1,k) .eq. nrr(2,k) )then
               sig(irev(k)) = 2.0d0*sig(irev(k))
               sigt(irev(k)) = 2.0d0*sigt(irev(k))
               sigv(irev(k)) = 1.0d0*sigv(irev(k))
            endif
            endif
c            write(*,*)k,irev(k),sig(k),sig(irev(k)),fak,qval(k),
c     1            qval(irev(k)),exp( -11.605d0*qval(k)/t9 ) / t9**1.5d0
c            endif
c            fak = etanuc(i1,kzone) + etanuc(i2,kzone) 
c     1           - etanuc(i3,kzone)-etanuc(i4,kzone)-etanuc(i5,kzone)
c     2           +11.605d0*qval(k)/t9
c            if( irev(k) .gt. 0 )then
c               write(*,'(2i5,5a6,5x,a5,a1,1p3e11.3,2i5,1p8e11.3)')
c     1              k,ideck(k),
c     1              (rname(i,k),i=1,5),rlkh(k),rvw(k),
c     2              sig(k),sig(irev(k)),fak,irev(k),iline(k),
c     3              qval(k),
c     4              y(i1)*y(i2)*sig(k)/(y(i3)*y(i4)*y(i5)*sig(irev(k)))
c            else
c               write(*,'(2i5,5a6,5x,a5,a1,1p3e11.3,2i5,1p8e11.3)')
c     1              k,ideck(k),
c     1              (rname(i,k),i=1,5),rlkh(k),rvw(k),
c     2              sig(k),sig(-irev(k)),fak,irev(k),iline(k),
c     3              11.605d0*qval(k)/t9,
c     4              y(i1)*y(i2)*sig(k)/(y(i3)*y(i4)*y(i5)*sig(-irev(k)))
c            endif

      enddo
      
c deck 7: exchange, i + j ---> k + l + m + n
c..use fkt for now ccccccccccccccccccccccccccccccccccccccccccccccccccc
      do k= k1deck(7),k2deck(7)
c..no reverse rates in 176 nuc network
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            if( irev(k) .gt. 0 )then
              fak=(pf(i3)*pf(i4)*pf(i5)*pf(i6)*(9.8678d9/rho)**2.0d0
     1        /pf(i1)/pf(i2))*(anuc(i3)*anuc(i4)*
     2        anuc(i5)*anuc(i6)/(anuc(i1)*anuc(i2)))**1.5d0
     2        *exp( -11.605d0*qval(k)/t9 ) *t9**3.0d0
              sig(irev(k)) = fak * sig(k)
              sigt(irev(k)) = fak * sigt(k)
              sigv(irev(k)) = fak * sigv(k)
c            endif
            if( nrr(1,k) .eq. nrr(2,k) )then
               sig(irev(k)) = 2.0d0*sig(irev(k))
               sigt(irev(k)) = 2.0d0*sigt(irev(k))
               sigv(irev(k)) = 1.0d0*sigv(irev(k))
            endif
            endif
c            write(*,'(2i5,6a6,5x,a5,a1,1p3e12.3,2i5)')
c     1           k,ideck(k),(rname(i,k),i=1,6),rlkh(k),rvw(k),
c     2           sig(k),sig(irev(k)),qval(k),irev(k),iline(k)

      enddo

c deck 8:  3 ---> 1
c..use fkt for now, it looks ok for triple alpha
      do k= k1deck(8),k2deck(8)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
c            i5 = nrr(5,k)
            if( irev(k) .gt. 0 )then
c               if (i5 .eq. 0) then
                fak=(pf(i4)/pf(i1)/pf(i2)/pf(i3))*(rho/9.8678d9)**2.0d0
     1          *(anuc(i4)/(anuc(i1)*anuc(i2)*anuc(i3)))**1.5d0
     2          *exp( -11.605d0*qval(k)/t9 ) /t9**3.0d0
                sig(irev(k)) = fak * sig(k)
                sigt(irev(k)) = fak * sigt(k)
                sigv(irev(k)) = fak * sigv(k)
c               else
c                fak=(pf(i4)*pf(i5)/pf(i1)/pf(i2)/pf(i3))*9.8678d9/rho
c     1          *(anuc(i1)*anuc(i2)*anuc(i3)/(anuc(i4)*anuc(i5)))**1.5d0
c     2          *exp( -11.605d0*qval(k)/t9 )*t9**1.5d0
c                sig(irev(k)) = fak * sig(k)
c                sigt(irev(k)) = fak * sigt(k)
c                sigv(irev(k)) = fak * sigv(k)
c               endif
c           endif
           if( nrr(1,k) .eq. nrr(2,k) .and. 
     1          nrr(2,k) .eq. nrr(3,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*6.0d0
              sigt(irev(k)) = sigt(irev(k))*6.0d0
              sigv(irev(k)) = sigv(irev(k))*3.0d0
           elseif(nrr(1,k) .eq. nrr(2,k) .or. nrr(1,k) .eq. nrr(3,k)
     1            .or. nrr(2,k) .eq. nrr(3,k) .and. irev(k) .gt. 0)then
               sig(irev(k)) = 2.0d0*sig(irev(k))
               sigt(irev(k)) = 2.0d0*sigt(irev(k))
               sigv(irev(k)) = 1.0d0*sigv(irev(k))
           endif
           endif
      enddo
c deck 9:  3 ---> 2
      do k= k1deck(9),k2deck(9)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            if( irev(k) .gt. 0 )then
                fak=(pf(i4)*pf(i5)/pf(i1)/pf(i2)/pf(i3))/9.8678d9*rho
     1          *(anuc(i4)*anuc(i5)/(anuc(i1)*anuc(i2)*anuc(i3)))**1.5d0
     2          *exp( -11.605d0*qval(k)/t9 )/t9**1.5d0
                sig(irev(k)) = fak * sig(k)
                sigt(irev(k)) = fak * sigt(k)
                sigv(irev(k)) = fak * sigv(k)    
c           endif
           if( nrr(1,k) .eq. nrr(2,k) .and. 
     1          nrr(2,k) .eq. nrr(3,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*6.0d0
              sigt(irev(k)) = sigt(irev(k))*6.0d0
              sigv(irev(k)) = sigv(irev(k))*3.0d0
           elseif(nrr(1,k) .eq. nrr(2,k) .or. nrr(1,k) .eq. nrr(3,k)
     1            .or. nrr(2,k) .eq. nrr(3,k) .and. irev(k) .gt. 0)then
               sig(irev(k)) = 2.0d0*sig(irev(k))
               sigt(irev(k)) = 2.0d0*sigt(irev(k))
               sigv(irev(k)) = 2.0d0*sigv(irev(k))
           elseif(nrr(1,k) .eq. nrr(2,k) .and. nrr(1,k) .eq. nrr(3,k)
     1            .or. nrr(2,k) .eq. nrr(3,k) .and. irev(k) .gt. 0)then
               sig(irev(k)) = 4.0d0*sig(irev(k))
               sigt(irev(k)) = 4.0d0*sigt(irev(k))
               sigv(irev(k)) = 2.0d0*sigv(irev(k))
           endif
c             write(*,*)k,irev(k),sig(k),fak,qval(k),qval(irev(k)),
c     1                exp( -11.605d0*qval(k)/t9 ) / t9**1.5d0
           endif
           
      enddo

c deck 10:  4 ---> 2
      do k= k1deck(10),k2deck(10)
            i1 = nrr(1,k)
            i2 = nrr(2,k)
            i3 = nrr(3,k)
            i4 = nrr(4,k)
            i5 = nrr(5,k)
            i6 = nrr(6,k)
            if( irev(k) .gt. 0 )then
               fak=(pf(i5)*pf(i6)/pf(i1)/pf(i2)/pf(i3)/pf(i4))
     1              *(rho/9.8678d9)**2.0d0
     1        *(anuc(i5)*anuc(i6)/(anuc(i1)*anuc(i2)*anuc(i3)*anuc(i4)))
     2         **1.5d0
     2         *exp( -11.605d0*qval(k)/t9 )/t9**3.0d0
               sig(irev(k)) = fak * sig(k)
               sigt(irev(k)) = fak * sigt(k)
               sigv(irev(k)) = fak * sigv(k)
c           endif
           if( nrr(1,k) .eq. nrr(2,k) .and. nrr(2,k) .eq. nrr(3,k) .and.
     1          nrr(3,k) .eq. nrr(4,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*24.0d0
              sigt(irev(k)) = sigt(irev(k))*24.0d0
              sigv(irev(k)) = sigv(irev(k))*8.0d0
           elseif( nrr(1,k) .eq. nrr(2,k) .and. 
     1          nrr(2,k) .eq. nrr(3,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*6.0d0
              sigt(irev(k)) = sigt(irev(k))*6.0d0
              sigv(irev(k)) = sigv(irev(k))*3.0d0
           elseif( nrr(1,k) .eq. nrr(3,k) .and. 
     1          nrr(3,k) .eq. nrr(4,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*6.0d0
              sigt(irev(k)) = sigt(irev(k))*6.0d0
              sigv(irev(k)) = sigv(irev(k))*3.0d0
           elseif( nrr(1,k) .eq. nrr(2,k) .and. 
     1          nrr(2,k) .eq. nrr(4,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*6.0d0
              sigt(irev(k)) = sigt(irev(k))*6.0d0
              sigv(irev(k)) = sigv(irev(k))*3.0d0
           elseif( nrr(2,k) .eq. nrr(3,k) .and. 
     1          nrr(3,k) .eq. nrr(4,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*6.0d0
              sigt(irev(k)) = sigt(irev(k))*6.0d0
              sigv(irev(k)) = sigv(irev(k))*3.0d0
           elseif( nrr(1,k) .eq. nrr(2,k) .and. 
     1          nrr(3,k) .eq. nrr(4,k) .and. irev(k) .gt. 0 )then
              sig(irev(k)) = sig(irev(k))*4.0d0
              sigt(irev(k)) = sigt(irev(k))*4.0d0
              sigv(irev(k)) = sigv(irev(k))*2.0d0
           elseif(nrr(1,k) .eq. nrr(2,k) .or. nrr(1,k) .eq. nrr(3,k)
     1           .or. nrr(1,k) .eq. nrr(4,k) .or. nrr(2,k) .eq. nrr(3,k)
     1           .or. nrr(2,k) .eq. nrr(4,k) .or. nrr(3,k) .eq. nrr(4,k)
     1           .and. irev(k) .gt. 0)then
               sig(irev(k)) = 2.0d0*sig(irev(k))
               sigt(irev(k)) = 2.0d0*sigt(irev(k))
               sigv(irev(k)) = 1.0d0*sigv(irev(k))
           endif
           endif
      enddo

c      do k=k1deck(6),k2deck(6)
c         write(*,*)sig(k)
c      enddo
c      do k=k1deck(9),k2deck(9)
c         write(*,*)sig(k)
c      enddo
c..limit rate links to mitigate roundoff
      do k= k1deck(4),k2deck(10)
         if( irev(k) .eq. 0 )then
c..no reverse rate so only forward rate needs scaling
            fak = sig(k) * dtstar
            if( fak .gt. linklim )then
               sig(k) = linklim/dtstar
            endif
c            write(*,'(2i6,7a6,1p8e12.3)')k,irev(k),
c     1           (rname(i,k),i=1,7),sig(k),sig(k)*dtstar,
c     2           fak,linklim/fak,linklim/dtstar
         elseif( irev(k) .gt. 0 )then
c..reverse rate has been defined, so it will need scaling too
            if( sig(k)*dtstar .gt. linklim .or. 
     1           sig(irev(k))*dtstar .gt. linklim )then
c               write(*,'(2i6,7a6,1p8e12.3)')k,irev(k),
c     1              (rname(i,k),i=1,7),sig(k),sig(irev(k)),
c     2              linklim/fak,linklim/dtstar,
c     3              sig(k)/( sig(k) + sig(irev(k)))*linklim/dtstar,
c     3              (1.0d0-sig(k)/( sig(k) + sig(irev(k))))
c     4              *linklim/dtstar
               fak = sig(k)/( sig(k) + sig(irev(k)))
c..   scale both rates, preserving ratio
               sig(k)       = linklim/dtstar* fak
               sig(irev(k)) = linklim/dtstar*(1.0d0 - fak)
            endif
         endif
c         write(*,*)k,sig(k),rname(1,k),rname(2,k),rname(3,k),linklim
      enddo
c      stop
c      write(*,*)'large rates scaled back'
c      do  k= k1deck(1),k2deck(8)
c         fak = sig(k) * dtstar
c         if( fak .gt. linklim )then
c            write(*,'(3i6,7a6,1p8e12.3)')k,ideck(k),irev(k),
c     1           (rname(i,k),i=1,7),sig(k),sig(k)*dtstar,
c     2           fak,linklim/fak,linklim/dtstar
c         endif
c      enddo
c............................................................
      
      return
      
      end


      subroutine screen(t9,rho,kz)

c..   weak, intermediate and strong screening wda 12/22/04
c..   after Graboske, et al, ApJ 181, page 465 (1973), table 4
c..   see also DeWitt, Graboske, Cooper, 1973, Apj companion paper
c..   and John Bahcall's exportenergy.f subroutine

c..   kz is zone index

      implicit none

      include 'dimenfile'
      include 'crate'
      include 'cburn'
      include 'comcsolve'

      real*8 t9,rho,opbexp,z1,z2,z12,z3,z123,z4,amu,emu,xtr,zet,pfmc2
      real*8 efmkt,fprf,degd,t9m32
      real*8 weakscrn
      parameter( opbexp= 1.86d0 , weakscrn = 0.03d0 )
      real*8 xxl,xxl6,xxl8,zcurl,zbar,z58,z28,z33,tm1,uwk
      real*8 uint,ustr

      real*8 f0(nreac),zprd(nreac),z86(nreac)
      real*8 zprd3(nreac),zprd4(nreac),z863(nreac),z864(nreac)
      real*8 z53(nreac),z43(nreac),z23(nreac)
      real*8 z533(nreac),z433(nreac),z233(nreac)
      real*8 z534(nreac),z434(nreac),z234(nreac)
      real*8 utot(nreac),dscr(nreac),dsct(nreac)

      integer*4 kz,i,k,calls
      integer*4 j

      save

      data calls/0/
c---------------------------------------------------------------

      if( calls .eq. 0 )then
c..only for first time through
         calls = calls + 1
c         write(*,*) 'first time in screen ',calls,kz
c..   set up constants
         do k = 1, nreac
            zprd(k) = 0.0d0
            z86(k)  = 0.0d0
            z53(k)  = 0.0d0
            z43(k)  = 0.0d0
            z23(k)  = 0.0d0
         enddo
c..deck 4 a + b ---> c
c..deck 5 a + b ---> c + d
c..deck 6 a + b ---> c + d + e
c..deck 7 a + b ---> c + d + e + f
         do k = k1deck(4),k2deck(7)
            z1      = nz(nrr(1,k))
            z2      = nz(nrr(2,k)) 
            zprd(k) = z1*z2 
            z86(k)  = (z1+z2)**opbexp - z1**opbexp - z2**opbexp
            z53(k)  = (z1+z2)**(5.0d0/3.0d0) - z1**(5.0d0/3.0d0) 
     1           - z2**(5.0d0/3.0d0)
            z43(k)  = (z1+z2)**(4.0d0/3.0d0) - z1**(4.0d0/3.0d0) 
     1           - z2**(4.0d0/3.0d0)
            z23(k)  = (z1+z2)**(2.0d0/3.0d0) - z1**(2.0d0/3.0d0) 
     1           - z2**(2.0d0/3.0d0)
         enddo
c..deck 8 a + b + c = d 
c..deck 9 a + b + c = d + e
         do k = k1deck(8),k2deck(9)
            z1      = nz(nrr(1,k))
            z2      = nz(nrr(2,k))
            z3      = nz(nrr(3,k))
            z12     = nz(nrr(1,k))+nz(nrr(2,k))
            zprd(k) = z1*z2 
            zprd3(k)= z12*z3
            z86(k)  = (z1+z2)**opbexp - z1**opbexp - z2**opbexp
            z863(k) = (z12+z3)**opbexp-z12**opbexp - z3**opbexp
            z53(k)  = (z1+z2)**(5.0d0/3.0d0) - z1**(5.0d0/3.0d0) 
     1           - z2**(5.0d0/3.0d0)
            z533(k) = (z12+z3)**(5.0d0/3.0d0) - z12**(5.0d0/3.0d0) 
     1           - z3**(5.0d0/3.0d0)
            z43(k)  = (z1+z2)**(4.0d0/3.0d0) - z1**(4.0d0/3.0d0) 
     1           - z2**(4.0d0/3.0d0)
            z433(k) = (z12+z3)**(4.0d0/3.0d0) - z12**(4.0d0/3.0d0) 
     1           - z3**(4.0d0/3.0d0)
            z23(k)  = (z1+z2)**(2.0d0/3.0d0) - z1**(2.0d0/3.0d0) 
     1           - z2**(2.0d0/3.0d0)
            z233(k)  = (z12+z3)**(2.0d0/3.0d0) - z12**(2.0d0/3.0d0) 
     1           - z3**(2.0d0/3.0d0)
         enddo
c      endif

c..deck 10 a + b + c +d = e + f 
         do k = k1deck(10),k2deck(10)
            z1      = nz(nrr(1,k))
            z2      = nz(nrr(2,k))
            z3      = nz(nrr(3,k))
            z4      = nz(nrr(4,k))
            z12     = nz(nrr(1,k))+nz(nrr(2,k))
            z123    = nz(nrr(1,k))+nz(nrr(2,k))+nz(nrr(3,k))
            zprd(k) = z1*z2 
            zprd3(k)= z12*z3
            zprd4(k)= z123*z4
            z86(k)  = (z1+z2)**opbexp - z1**opbexp - z2**opbexp
            z863(k) = (z12+z3)**opbexp-z12**opbexp - z3**opbexp
            z864(k) = (z123+z4)**opbexp-z123**opbexp - z4**opbexp
            z53(k)  = (z1+z2)**(5.0d0/3.0d0) - z1**(5.0d0/3.0d0) 
     1           - z2**(5.0d0/3.0d0)
            z533(k) = (z12+z3)**(5.0d0/3.0d0) - z12**(5.0d0/3.0d0) 
     1           - z3**(5.0d0/3.0d0)
            z534(k) = (z123+z4)**(5.0d0/3.0d0) - z123**(5.0d0/3.0d0) 
     1           - z4**(5.0d0/3.0d0)
            z43(k)  = (z1+z2)**(4.0d0/3.0d0) - z1**(4.0d0/3.0d0) 
     1           - z2**(4.0d0/3.0d0)
            z433(k) = (z12+z3)**(4.0d0/3.0d0) - z12**(4.0d0/3.0d0) 
     1           - z3**(4.0d0/3.0d0)
            z434(k) = (z123+z4)**(4.0d0/3.0d0) - z123**(4.0d0/3.0d0) 
     1           - z4**(4.0d0/3.0d0)
            z23(k)  = (z1+z2)**(2.0d0/3.0d0) - z1**(2.0d0/3.0d0) 
     1           - z2**(2.0d0/3.0d0)
            z233(k)  = (z12+z3)**(2.0d0/3.0d0) - z12**(2.0d0/3.0d0) 
     1           - z3**(2.0d0/3.0d0)
            z234(k)  = (z123+z4)**(2.0d0/3.0d0) - z123**(2.0d0/3.0d0) 
     1           - z4**(2.0d0/3.0d0)
         enddo
      endif


      do i=1, nreac
         f0(i)   = 1.0d0
         utot(i) = 0.0d0
      enddo
c..   amu=Ytot,emu=Ye,xtr is an interpolation variable,zet=ave Z**2
      amu = 0.0d0
      emu = 0.0d0
      xtr = 0.0d0
      zet = 0.0d0
      do i = 1, inuc
         amu = amu + y(i)
         emu = emu + y(i)*dble( nz(i) )
         xtr = xtr + y(i)*dble( nz(i) )**1.58d0
         zet = zet + y(i)*dble( nz(i) )**2
      enddo

      pfmc2 = 1.017677d-4*( rho*emu )**0.666666667d0
      efmkt = 5.92986d0/t9*( sqrt( 1.0d0 + pfmc2 ) - 1.0d0 )
      if( efmkt .le. 1.0d-2 )then
         fprf = 1.0d0
      else
         degd = log10( efmkt )
         if( degd .ge. 1.5d0 )then
            fprf = 0.0d0
         else
            fprf =0.75793-0.54621*DEGD-0.30964*DEGD**2+0.12535*DEGD**3+
     1      0.1203*DEGD**4-0.012857*DEGD**5-0.014768*DEGD**6
         endif
      endif

      t9m32 = 1.0d0/t9**1.5d0
      XXL   = 5.9426d-6*T9M32*DSQRT(RHO*AMU)
      XXL6  = XXL**0.666667d0
      XXL8  = XXL**0.86d0
      ZCURL = DSQRT((ZET+FPRF*EMU)/AMU)
      ZBAR  = EMU/AMU
      Z58   = ZCURL**0.58d0
      Z28   = ZBAR**0.28d0
      Z33   = ZBAR**0.333333d0
      TM1   = XXL*ZCURL

c..for each reaction
c..deck 4 a + b ---> c
c..deck 5 a + b ---> c + d
c..deck 6 a + b ---> c + d + e
c..deck 7 a + b ---> c + d + e + f
      do i = k1deck(4),k2deck(7)
         uwk  = tm1*zprd(i)
         uint = 0
         ustr = 0
         if( uwk .le. weakscrn )then
            utot(i) = uwk
            dscr(i) =  0.5d0*uwk
            dsct(i) = -1.5d0*uwk
         else
            uint = 0.38d0*xxl8*xtr*z86(i)/(amu*z58*z28)
            if( uint .le. 2.0d0 )then
               utot(i) = uint
               dscr(i) =  0.43d0*uint
               dsct(i) = -1.29d0*uint
             else
                ustr = 0.624d0*z33*xxl6*(z53(i)+0.316d0*z33*z43(i)
     1               +0.737d0*z23(i)/(zbar*xxl6))
                if( ustr .lt. uint .or. uwk .ge. 5.0d0 )then
                   utot(i) = ustr
                   dscr(i) = 0.208d0*z33*(z53(i) + 0.316d0*z33*z43(i))
     1                  *xxl6
                   dsct(i) = -3.0d0*dscr(i)
                else
                   utot(i) = uint
                   dscr(i) =  0.43d0*uint
                   dsct(i) = -1.29d0*uint
                endif
             endif
         endif
      enddo

c..deck 8 a + b + c ---> d
c..deck 9 a + b + c ---> d + e
      do i = k1deck(8),k2deck(9)
         uwk  = tm1*zprd(i)+tm1*zprd3(i)
         uint = 0
         ustr = 0
         if( uwk .le. weakscrn )then
            utot(i) = uwk
            dscr(i) =  0.5d0*uwk
            dsct(i) = -1.5d0*uwk
         else
            uint = 0.38d0*xxl8*xtr*z86(i)/(amu*z58*z28)
     1             +0.38d0*xxl8*xtr*z863(i)/(amu*z58*z28)
            if( uint .le. 2.0d0 )then
               utot(i) = uint
               dscr(i) =  0.43d0*uint
               dsct(i) = -1.29d0*uint
             else
                ustr = 0.624d0*z33*xxl6*(z53(i)+0.316d0*z33*z43(i)
     1               +0.737d0*z23(i)/(zbar*xxl6))
     2               +0.624d0*z33*xxl6*(z533(i)+0.316d0*z33*z433(i)
     1               +0.737d0*z233(i)/(zbar*xxl6))    
                if( ustr .lt. uint .or. uwk .ge. 5.0d0 )then
                   utot(i) = ustr
                   dscr(i) = 0.208d0*z33*(z53(i) + 0.316d0*z33*z43(i))
     1                  *xxl6
     2                  +0.208d0*z33*(z533(i) + 0.316d0*z33*z433(i))
     1                  *xxl6
                   dsct(i) = -3.0d0*dscr(i)
                else
                   utot(i) = uint
                   dscr(i) =  0.43d0*uint
                   dsct(i) = -1.29d0*uint
                endif
             endif
         endif
      enddo

c..deck 10 a + b + c + d ---> e + f
      do i = k1deck(10),k2deck(10)
         uwk  = tm1*zprd(i)+tm1*zprd3(i)+tm1*zprd4(i)
         uint = 0
         ustr = 0
         if( uwk .le. weakscrn )then
            utot(i) = uwk
            dscr(i) =  0.5d0*uwk
            dsct(i) = -1.5d0*uwk
         else
            uint = 0.38d0*xxl8*xtr*z86(i)/(amu*z58*z28)
     1             +0.38d0*xxl8*xtr*z863(i)/(amu*z58*z28)
     2             +0.38d0*xxl8*xtr*z864(i)/(amu*z58*z28)
            if( uint .le. 2.0d0 )then
               utot(i) = uint
               dscr(i) =  0.43d0*uint
               dsct(i) = -1.29d0*uint
             else
                ustr = 0.624d0*z33*xxl6*(z53(i)+0.316d0*z33*z43(i)
     1               +0.737d0*z23(i)/(zbar*xxl6))
     2               +0.624d0*z33*xxl6*(z533(i)+0.316d0*z33*z433(i)
     1               +0.737d0*z233(i)/(zbar*xxl6))   
     2               +0.624d0*z33*xxl6*(z534(i)+0.316d0*z33*z434(i)
     1               +0.737d0*z234(i)/(zbar*xxl6))   
                if( ustr .lt. uint .or. uwk .ge. 5.0d0 )then
                   utot(i) = ustr
                   dscr(i) = 0.208d0*z33*(z53(i) + 0.316d0*z33*z43(i))
     1                  *xxl6
     2                  +0.208d0*z33*(z533(i) + 0.316d0*z33*z433(i))
     1                  *xxl6
     2                  +0.208d0*z33*(z534(i) + 0.316d0*z33*z434(i))
     1                  *xxl6
                   dsct(i) = -3.0d0*dscr(i)
                else
                   utot(i) = uint
                   dscr(i) =  0.43d0*uint
                   dsct(i) = -1.29d0*uint
                endif
             endif
         endif
      enddo

      do i = k1deck(4),k2deck(10)         
            f0(i) = exp (utot(i))
c         if( i .eq. 619 .or. i .eq. 620 .or. i .eq. 627 .or.
c     1        i .eq. 643 .or. i .eq. 644 .or. i .eq. 647 .or.
c     2        i .eq. 648 .or. i .eq. 652 .or. i .eq. 653 .or.
c     3        i .eq. 664 .and. kz .eq. 2 )then
c            write(*,'(i5, 3a5, 1p8e12.4)') i,
c     1           (xid(nrr(j,i)),j=1,3), zprd(i),f0(i),sig(i),
c     2           f0(i)*sig(i)
c         endif
      enddo


c..   adjust for screening
c..   can add T and V derivatives from f0 here
      do i = 1,nreac
           sig(i)  = sig(i)  * f0(i)
           sigt(i) = sigt(i) * f0(i)
           sigv(i) = sigv(i) * f0(i)
      enddo


c.....for solar model studies:..................................
c..   hardwired for 176 element network
c..   revise fkt rates to jnb consistency 
c..   used genrate to get factor for bare NA sig v ratio
c..   use genrate to identify reaction indices for new network
c      write(*,'(//a50)')
c     1     'screen.f: pp and cno reactions, full screening'
c      write(*,'(3a6,7a6,8a12)')'ireac','deck','irev',
c     1     '1','2','3','4','5','6',
c     2     'rlkh','NA sigma v','log life(s)'
c..   p + p ---> d by electron capture
c..   according to jn bahcall 1995
c      k   = 619
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + p ---> d by positron decay
c      k = 620
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + d ---> he3
c      k = 621
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   he3 + he3 ---> p + p + he4
c      k = 1840
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   he3 + he4 ---> be7
c      k = 627
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   be7 ---> li7
c      k = 3
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + li7 --> he4 + he4
c      k = 1086
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + c12
c      k = 643
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   n13 --> c13
c      k = 5
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + c13
c      k = 647
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + n14
c      k = 652
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + n15 --> he4 + c12
c      k = 1139
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + n15 --> o16
c      k = 658
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + o16
c      k = 664
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   f17 -->o17
c      k = 9
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c...  p+o17 -->He4 + n14
c      k = 1165
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   p + o17
c      k = 669
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c..   f18 -->o18
c      k = 10
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c...  p+o18 -->He4 + n15
c      k = 1171
c      write(*,'(3i6,7a6,1pe12.3,0p2f10.5)')k,ideck(k),irev(k),
c     1     (rname(i,k),i=1,6),rlkh(k),sig(k),-log10( sig(k)*100.0d0 )
c     2     ,f0(k)
c......................................................................

      return
      end



