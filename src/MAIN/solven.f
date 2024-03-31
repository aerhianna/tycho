      subroutine solven(dtstar,dth,k,
     1     idebug,ncytest,ncymax,inegs,it,kk,leq,jnb)

c..   evaluates rates, sets up matrix equations,
c..   and solves reaction network over time interval dth,
c..   updating mole fractions y and nucleon fractions x
c..   derivatives with respect to T, V (5/9/99)
c..   revised 7-7-00 for low T rates

      implicit none

      include 'dimenfile'

      include 'crate'
      include 'comcsolve'
      include 'caeps'
      include 'cdtnuc'
      include 'cdeuter'
      include 'cburn'
      include 'cconst'

      real*8  yold(ndim)
c..   instantaneous rate of energy generation
      real*8  ssi(kdm)

      real*8 tiny, t9, rho, dtstar, dth
      real*8 time, dtleft, dth0, rdt, sumd, suma, eb, dn, en
      real*8 fak, dne, yneu, ypro, xalf,  xnucleu
      real*8 dxnucl
      real*8 yestim,worsty
      real*8 enc,yoldsum,ynewsum,yesum, t9old(kdm)

      integer*4 inegs
      integer*4 idebug, ncycle, k, j, negflag, leq, jnb
      integer*4 it, kk, kreac

      integer*4 indx(ndim), i
      real*8 foo,alud(ndim,ndim)
cccccccccccccccccc


      integer*4 k1,k2,i1,i2
      integer*4 ncytest,ncymax
      integer*4 icall

c..   sparse solver variables
      integer*4 nelem,iflag
      parameter(nelem = nreac)
      real*8 altb(nnuc),berr
      

c..   loburn network
c..   idlo is initialized negative as a flag to loburn.f
c..   note strange fortran form for repetition of variable below (*-1)
      data idlo/lodim*-1/
      data loz/1,2,3,3,4,4,5,5,6,6, 7,7,7,8,8,8,8,9,9,9, 
     1     10,10,10,11,11,11,12,12,12,13, 13, 0,1,2,0/
      data lon/1,1,3,4,3,5,5,6,6,7, 6,7,8,7,8,9,10,8,9,10,
     1     10,11,12,10,11,12,12,13,14,13, 14, 1,0,2,0/
      data icall/0/
      save
c---------------------------------------------------------------
c..   floor on abundance variable used in network
      parameter( tiny=1.0d-100 )

c..   input: idebug,t9,rho,dtstar,dyc,modeb,k
c..   output: dth,yold,ncycle
c     aeps,y               (via commons caeps, comcsolve)
c     aepst,aepsv

c-------------------------------------------------------------

      dth0   = 0.0d0
      inegs  = 0

c---------------------------------------------------------------
c..   full network solution with abundance update (implicit), or
c..   energy generation rate from right hand side evaluation (explicit)

      ncycle = 0
c..   arrays passed through common caeps, set in burn.f
      t9     = temp(k)*1.0d-9
      rho    = den(k)
      if(icall .eq. 0)then
         t9old(k) = t9
      endif
      if( itbrn(k) .eq. 2 )then
c..   full network
         do j = 1, itot
            y(j)    = aex(j,k)
            b(j)    = 0.0d0
            altb(j) = 0.0d0
            bt(j)   = 0.0d0
            bv(j)   = 0.0d0
            yold(j) = y(j)
            ydcon(j) = axdcon(j,k)
         enddo
         
c..   evaluate rates and put in sig array

         call rate(t9,rho,dtstar,k,jnb)

c..   al26 rate: fix at low T
c..   reaclib fit is bad below 2.9d7 K
c         if( t9 .lt. 0.029d0 )then
c            sig(lkal26) = 3.051d-14 
c         endif
c         write(*,'(2i5,1p8e12.3)')k,lkal26,t9*1.0d9,rho,sig(lkal26),
c     1        y(lal26),dlog(2.0d0)/(7.2d5*secpy)
c         if( k .ge. 480 )stop'solven asd'
cccccccccccccc


c..   weak, intermediate, and strong screening
         call screen(t9,rho,k)
         
         time   = 0
         dtleft = dtstar
c..   most optimistic choice for time step
         dth    = dtleft

c..   ncycle: subcycle over shorter timesteps for network if needed
 100     continue

         ncycle = ncycle + 1

         if( ncycle .gt. ncymax )goto 101

c..   insure time interval does not overshoot
         if( dth .gt. dtleft )then
            dth = dtleft
         endif

         if( dth .le. 0.0d0 )then
            write(*,'(5a12)')
     1           'dth','dth0','dtleft','ncycle','k'
            write(*,'(1p3e12.3,2i12)')dth,dth0,dtleft,ncycle,k
            stop'solven dth'
         endif

         dth0   = dth

c..   setup matrices for linear equation solver
c..   need rates, abundances Y, partition functions (2J+1=w factors)
         call rhside(idebug)
c         do j=1,itot
c            altb(j) = b(j)
c         enddo
c         write(*,*)y
c         write(*,*)' '
c         write(*,*)b
         if( ncycle .le. 1 )then
c..   estimate first time step
c     call dtnuc(dth0,dth,0)
c..   first guess to be used for core collapse conditions
c     dth = dmax1(dth,0.5d0*dtstar)
            dth = dtstar
            if( dth .gt. dtleft )then
               dth = dtleft
            endif
            dth0   = dth


         endif
c..   calculate the instantaneous rate
         eb   = 0
         dn   = 0
         do j = 1, itot
            fak =  b(j)
            eb  = eb - qq(j)*fak
            dn  = dn +       fak
         enddo
         eb      = eb *   9.65d+17    
         dne     = dn * 1.2476d+17 
         ssi(k) = eb - dne

ccccccccccccccccccccccccccccccccccccccc
c..   add convective change
c     do j = 1, itot
c     b(j) = b(j) + ydcon(j)/dtstar
c     enddo
ccccccccccccccccccccccccccccccccccccccccccc


c..   continue with network solution
         call jacob

c..   update diagonal terms
         rdt = 1.0d0/dth
         do j = 1, itot
            a(j,j) = a(j,j) + rdt
         enddo
         
c..   SPARSE SOLVER SECTION
c..   put nonzero entries of A matrix into sparse_dfdy vector
         do j=1,nlinks
            sparse_dfdy(j) = a(iloc(j),jloc(j))
         enddo

c..   sparseu controls bias towards sparsity or numerical pivoting
c..   see sparse_ma28.f
         sparseu = 0.1d0

         call ma28bd(nnuc,nlinks,sparse_dfdy,nelem,iloc,jloc,jvect,
     1               ikeep,iw,sparsew,iflag)
         if (iflag .lt. 0) then
            write(*,*) 'error in ma28bd flag',iflag
            stop 'error in ma28bd in solven'
        endif

        call ma28cd(nnuc,sparse_dfdy,nelem,jvect,ikeep,b,
     1               sparsew,1)
         
c     leqs solves linear equations
c          call leqs(a,b,itot,nnuc)

c     call leqsnr(a,b,itot,nnuc)

c..   LU decomposition
c     call ludcmp(a,itot,itot,indx,foo)
c..   save LU decomposition in alud
c     do i = 1, nnuc
c     do j = 1, nnuc
c     alud(i,j) = a(i,j)
c     enddo
c     enddo
c..   back substitution
c     call lubksb(alub,nnuc,nnuc,indx,b)
c..   d is original rhs
c..   c is original a matrix
c..   alud is LU decomposition of original a matrix
c..   f is the solution vector from lubksb
c..   b is new solution for abundance x
c     call mprove(c,alud,nnuc,ndim,indx,d,f)

c..   evaluate next time step
         call dtnuc(dth0,dth,1)

c..   check result for excessive changes
         sumd = 0
         suma = 0
         sumx = 0
         sumz = 0
         do  j = 1, itot
            sumd = sumd +         b(j)  * dble( nn(j) + nz(j) )
            suma = suma +     abs(b(j)) * dble( nn(j) + nz(j) )
            sumx = sumx + (y(j) + b(j)) * dble( nz(j) + nn(j) )
         enddo
         
c..   catches worst negative abundance
         negflag = 0
         worsty  = 0.0d0
         do j = 1, itot
            yestim = y(j) + b(j)
            if( yestim .lt. -3.0d-16 )then
               negflag = j
               worsty = yestim
            endif
         enddo
         
c..   adjust time step
         if( negflag .gt. 0 )then
            if( b(negflag) .lt. 0.0d0 )then
c               if( ncycle .gt. 1 .and. b(negflag) .lt. -1.0d-15 )then
ccccccccccccccccccccccccccccccc
               if( ncycle .gt. 2 .and. b(negflag) .lt. -1.0d-14 )then

c..   only write if there is continued trouble and
c..   value is significantly different from roundoff
c                  write(*,'(a20,3i5,a5,1p8e12.3)')'negflag',k,ncycle,
c     1                 negflag,cnuc(negflag),dth,y(negflag),b(negflag)
c                  write(*,'(a20,a5,a3,i5,a7,i3,3(a10,1pe11.3))')
c     1                 'SOLVEN: negflag',cnuc(negflag),'k',k,
c     2                 'cycle',ncycle,'dth',dth,'y(neg)',
c     3                 y(negflag),'dy(neg)',b(negflag)

cccccccccccccccccccccccccccc
               endif

c..   would give 0.5 of initial value if linear decline
               if( y(negflag) .lt. 0.0d0 )then
c..   avoid zero abundance
                  dth = -0.5d0*dth0*y(negflag)/b(negflag)
               endif

            else
               write(*,*)'solven error: negflag',
     1              negflag,b(negflag)
               stop' solven negflag'
            endif
         endif

         if(  abs(sumd) .lt. 1.0d-7
     1        .and. negflag .eq. 0 
     2        .and. dth .ge. dth0*0.5d0 )then

c..   tests were successful
c..   so update particle densities, time elapsed
            sumx = 0
            sumz = 0
            eta  = 0
            do  j = 1, itot
               y(j) = y(j) + b(j)
c..   positive or zero
               y(j) = dmax1( y(j), 0.0d0 )
               x(j) = y(j) * dble( nz(j)+nn(j) )
               sumx = sumx + x(j)
               sumz = sumz + y(j)* dble( nz(j) )
               eta  = eta  + y(j)* dble(nn(j) - nz(j))

            enddo

c..   elapsed time
            time     = time + dth0
            dtleft   = dtstar - time
            checksum = sumx - 1.0d0

c..   output if excessive subcycles
            if( ncycle .ge. ncytest )then
               if( nucleu .gt. 0 .and. nucleu .le. itot )then

                  if( ncycle .eq. ncytest )
     1                 write(*,'(/2a4,a5,a11,10a9)')
     1                 'k','ncy','nuc','y','tau1','tau2','tau3',
     1                 'time','dth0','dth','dtleft'
     2                 ,'suma','sumd','check'
                  write(*,'(2i4,a5,1pe11.4,1p10e9.2)')
     1                 k,ncycle,xid(nucleu),y(nucleu),tau,time,
     1                 dth0,dth,dtleft
     2                 ,suma,sumd,checksum
                  write(3,'(2i4,a5,1pe11.4,1p10e9.2)')
     1                 k,ncycle,xid(nucleu),y(nucleu),tau,time,
     1                 dth0,dth,dtleft
     2                 ,suma,sumd,checksum
               else
                  write(*,'(2i6,12x,1p10e11.3)')
     1                 k,nucleu,tau,time,dth0,dth,dtleft
     2                 ,suma,sumd
                  write(3,'(2i6,12x,1p10e11.3)')
     1                 k,nucleu,tau,time,dth0,dth,dtleft
     2                 ,suma,sumd
               endif
            endif

c..   calculate energy generation rate using mass excesses
c..   DNE is 3/2 kT for new particles
            eb   = 0.0d0
            dn   = 0.0d0
            en   = 0.0d0
            do j = 1, itot
               fak =  y(j) - yold(j)
               eb  = eb - qq(j)*fak
               dn  = dn +       fak
               k1 = k1deck(1)
               k2 = k2deck(2)
               do kreac = k1, k2
                  if(nrr(1,kreac) .eq. j)then
                     en = en - avagadro*1.6021d-6*signue(kreac)*fak
                  endif
               enddo
            enddo
            eb      = eb *   9.65d+17    /(dtstar-dtleft)
            dne     = dn * 1.2476d+17 *t9/(dtstar-dtleft)
            en      = en / (dtstar-dtleft)
            aeps(k) = eb - dne + en

c..   calculate the neutrino emissivity
            aenucnu(k) = 0.0d0

            if( jnb .eq. 0 )then

c..   ffnu rates
               do j = k1deck(1), k2deck(1)
c..   sum the neutrino and antineutrino emission from ffnu rates
                  if( irev(j) .gt. 0 )then
                     aenucnu(k) = aenucnu(k) - 
     1                    ( y(nrr(1,j))*signue(j)
     2                    + y(nrr(2,j))*signue(irev(j)) )*cergs
                  endif
               enddo
            else

c..   Bahcall options
c..   check rhside.f for consistency in aenucnu component terms
c..   calculate the total nuclear neutrino energy emission
               aenucnu(k) = 0.0d0
c..   deck 1 has decays for CNO cycle
               do  j = k1deck(1), k2deck(1)
                  aenucnu(k) = aenucnu(k) - 
     1                 y(nrr(1,j))*signue(j)*cergs
c     write(*,'(i5,2a5,1p8e12.3)')j,xid(nrr(1,j)),
c     1                 xid(nrr(2,j)),
c     1                 signue(j), y(nrr(1,j)),sig(j),
c     2                 signue(j)*y(nrr(1,j))*cergs
               enddo
c..   deck 4 has pp--->d
               do j = k1deck(4),k2deck(4)
                  aenucnu(k) = aenucnu(k) 
     1                 - signue(j)*y(nrr(1,j))*y(nrr(2,j))*cergs
                  
c     write(*,'(i5,2a5,1p8e12.3)')j,xid(nrr(1,j)),
c     1                 xid(nrr(2,j)),
c     1                 signue(j), y(nrr(1,j)),y(nrr(2,j)),sig(j),
c     2                 signue(j)*y(nrr(1,j))*y(nrr(2,j))*cergs
               enddo

            endif

         else
c..   tests were unsuccessful
            dth = dmin1( dth, 0.5d0*dth0 )
            if( ncycle .ge. ncytest )then
               if( ncycle .eq. ncytest )
     1              write(*,'(a7,a10,a4,a10,
     1              3a5,7a10,a5)')' ','dth0',
     1              ' ','dth', 'k','ncy','xid','X','dX','dtlef',
     1              'yp','xa','suma','sumd','negf'
               yneu = y(itot-2) + b(itot-2)
               ypro = y(itot-1) + b(itot-1)
               xalf = 4.0d0*( y(itot) + b(itot) )
               if( nucleu .le. 0 .or. nucleu .gt. itot )then
                  write(*,'(2(a7,1pe9.1),3i5,1p7e10.2,i5)')
     1                 'reject',dth0,' try',dth, k,ncycle,nucleu, 
     2                 t9,dtleft,yneu,ypro,xalf,suma,sumd,negflag 
               else
                  xnucleu = y(nucleu)*
     1                 dble( nz(nucleu) + nn(nucleu))
                  dxnucl  = b(nucleu)*
     1                 dble( nz(nucleu) + nn(nucleu))
                  write(*,'(a7,1pe10.3,a4,1pe10.3,
     1                 2i5,a5,1p8e10.2,i5)')
     1                 'reject',dth0,' try',dth, k,ncycle,
     1                 xid(nucleu),xnucleu,dxnucl,dtleft,
     2                 yneu,ypro,xalf,suma,sumd,negflag
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  if( negflag .gt. 0 )
     1                 write(*,*)xid(negflag),y(negflag),
     1                 b(negflag),y(negflag)+b(negflag)

                  do j = 1, nreac
                     if( nrr(1,j) .eq. negflag .and.
     1                    sig(j) .gt. 0.0d0 )then
                        write(*,'(4i5,2a5,1p8e12.3)')negflag,nrr(1,j),
     1                       j,ideck(j),xid(nrr(1,j)),
     1                       xid(nrr(2,j)),
     1                       sig(j),sig(j)*dth,y(nrr(1,j))
                     endif
                  enddo
                  write(3,'(a7,1pe10.3,a4,1pe10.3,2i5,a5,
     1                 1p7e10.2,i5)')
     1                 'reject',dth0,' try',dth, k,ncycle,
     1                 xid(nucleu),xnucleu,dxnucl,dtleft,
     2                 ypro,xalf,suma,sumd,negflag
                  stop'solven: no success'
               endif
            endif
         endif
         
         if( dtleft .gt. 1.0d-14*dtstar )goto 100

         if( ncycle .gt. ncytest )then
            yneu = y(itot-2) + b(itot-2)
            ypro = y(itot-1) + b(itot-1)
            xalf = 4.0d0*( y(itot) + b(itot) )
            if( nucleu .le. 0 .or. nucleu .gt. itot )then
               write(*,'(2(a7,1pe9.1),3i5,1p7e10.2,i5)')
     1              'last',dth0,' try',dth, k,ncycle,nucleu,
     2              t9,dtleft,yneu,ypro,xalf,suma,sumd,negflag
            else
               xnucleu = y(nucleu)*dble( nz(nucleu) + nn(nucleu)) 
               dxnucl  = b(nucleu)*dble( nz(nucleu) + nn(nucleu)) 
               write(*,'(2(a7,1pe11.3),2i5,a5,1p7e10.2,i5)')
     1              'last',dth0,' try',dth, k,ncycle,xid(nucleu), 
     2              xnucleu,dxnucl,dtleft,ypro,xalf,suma,sumd,
     3              negflag
               write(3,'(a5,1pe10.2,a4,1pe10.2,2i5,a5,1p7e10.2,
     1              i5)')
     1              'last',dth0,' try',dth, k,ncycle,xid(nucleu), 
     2              xnucleu,dxnucl,dtleft,ypro,xalf,suma,sumd,
     3              negflag
            endif
         endif

c..   total matrix solutions
         ncytot = ncytot + ncycle


      elseif( itbrn(k) .eq. 3 )then
c..   itbrn(i) = 3 ...............................................
c..   nse solver plus weak interactions......................
         
         ncycle = 0
         rho = den(k)
         do j = 1, itot
            y(j)    = aex(j,k)
            b(j)    = 0.0d0
            yold(j) = y(j)
         enddo
c     write(*,*)t9,t9old(k)
c     if(t9 .lt. t9old(k))then
c     t9 = t9old(k)
c     endif
c..   evaluate rates and put in sig array
         call rate(t9,rho,dtstar,k,jnb)

         call screen(t9,rho,k)

         write(*,'(a20,3(a12,1pe12.3))')'qse/nse solver',
     1        't9',t9,'tnse',tnse
         write(*,*)'logic not complete'
         stop'solven: nse,qse'
cccccccccccccccccccccc

c     if(t9 .gt. tnse)then
         call nse(t9,rho,y,enc,k)
c     else
c     call sqse(t9,rho,y,enc,k)
c     endif

c     write(*,*)'ye',y(ndim)
c     call state(k,k+1,1)

c..   most optimistic choice for time step
         dth    = dtstar

c..   deck 1: i ---> j
c..   beta decays, positron decays,
c..   i + e ---> j for electron captures (see rate.f)
         weaksum = 0.0d0
         k1 = k1deck(1)
         k2 = k2deck(1)
         if( k1 .ne. 0 .and. k2 .ne. 0 )then
            do j = k1,k2
               i1 = nrr(1,j)
               i2 = nrr(2,j)
               b(i1)    = b(i1)    - y(i1)*sig(j)
               b(i2)    = b(i2)    + y(i1)*sig(j)
               weaksum = weaksum + ( nz(i2)-nz(i1) )*y(i1)*sig(j)
c..   dYe/dt is weaksum
               bt(i1)    = bt(i1)    - y(i1)*sigt(j)
               bv(i1)    = bv(i1)    - y(i1)*sigv(j)
               bt(i2)    = bt(i2)    + y(i1)*sigt(j)
               bv(i2)    = bv(i2)    + y(i1)*sigv(j)
            enddo
         endif

c..   setup matrices for linear equation solver
c..   need rates, abundances Y, partition functions (2J+1=w factors)
         call rhside(idebug)

c..   continue with network solution
         call jacob

c..   update diagonal terms
         rdt = 1.0d0/dth
         do j = 1, itot
            a(j,j) = a(j,j) + rdt
         enddo

c..   leqs solves linear equations
         call leqs(a,b,itot,nnuc)

c..   total matrix solutions done
         ncytot = ncytot + 1

c..   update particle densities, time elapsed
         sumx = 0
         sumz = 0
         eta  = 0
         do  j = 1, itot
            y(j) = y(j) + b(j)
c..   positive or zero
            y(j) = dmax1( y(j), 0.0d0 )
            x(j) = y(j) * dble( nz(j) + nn(j) )
            sumx = sumx + x(j)
            sumz = sumz + y(j)* dble( nz(j) )
            eta  = eta  + y(j)* dble(nn(j) - nz(j))
         enddo
         checksum = sumx - 1.0d0
         write(*,*)'checksum in qse ',checksum
c..   Brand new thesis paranoia
c..   renormalize abundances
         sumx = -1.0d0
         do j = 1, nnuc
            x(j) = y(j) * dble( nz(j) + nn(j) )
            sumx = sumx + x(j)
         enddo

         do j = 1, ndim-1
            x(j) = x(j)/( 1.0d0 + sumx )
         enddo
c..   reset Ye
         yesum = 0.0d0
         sumz   = -1.0d0
         do j = 1, ndim-1
            yesum = yesum + y(j) * dble( nz(j) )
            sumz = sumz + x(j)
         enddo
         x(ndim) = yesum
         y(ndim) = yesum

c..   calculate energy generation rate using mass excesses
c..   DNE is 3/2 kT for new particles
         eb   = 0
         dn   = 0
         yoldsum = 0.0d0
         ynewsum = 0.0d0
         do j = 1, itot
            yoldsum = yoldsum+yold(j)
            ynewsum = ynewsum+y(j)
            fak =  y(j) - yold(j)
            eb  = eb - qq(j)*fak
            dn  = dn +       fak
         enddo

         eb      = eb *   9.65d+17   /dtstar
         dne      = dn * 1.2476d+17 *t9/dtstar
         aeps(k) = (eb - dne)/rho
c     aeps(k) = enc/dtstar

c     write(*,*)'enc',enc/rho,aeps(k)
c     temp(k) = ynewsum/yoldsum*temp(k)
c     aeps(k) = enc
c     call state(k,k+1,1)

      elseif( itbrn(k) .eq. 0 )then
c..   itbrn(i) = 0 ...............................................
c..   network solution for decays only......................
         ncycle = 0
         rho = den(k)
         do j = 1, itot
            y(j)    = aex(j,k)
            b(j)    = 0.0d0
            yold(j) = y(j)
         enddo

c..   evaluate rates and put in sig array
         call rated(t9,rho)

c..   most optimistic choice for time step
         dth    = dtstar

c..   deck 1: i ---> j
c..   beta decays, positron decays,
c..   i + e ---> j for electron captures (see rate.f)
         weaksum = 0.0d0
         k1 = k1deck(1)
         k2 = k2deck(1)
         if( k1 .ne. 0 .and. k2 .ne. 0 )then
            do j = k1,k2
               i1 = nrr(1,j)
               i2 = nrr(2,j)
               b(i1)    = b(i1)    - y(i1)*sig(j)
               b(i2)    = b(i2)    + y(i1)*sig(j)
               weaksum = weaksum + ( nz(i2)-nz(i1) )*y(i1)*sig(j)
c..   dYe/dt is weaksum
               bt(i1)    = bt(i1)    - y(i1)*sigt(j)
               bv(i1)    = bv(i1)    - y(i1)*sigv(j)
               bt(i2)    = bt(i2)    + y(i1)*sigt(j)
               bv(i2)    = bv(i2)    + y(i1)*sigv(j)
            enddo
         endif

c..   setup matrices for linear equation solver
c..   need rates, abundances Y, partition functions (2J+1=w factors)
         call rhside(idebug)

c..   continue with network solution
         call jacob

c..   update diagonal terms
         rdt = 1.0d0/dth
         do j = 1, itot
            a(j,j) = a(j,j) + rdt
         enddo

c..   leqs solves linear equations
         call leqs(a,b,itot,nnuc)

c..   total matrix solutions done
         ncytot = ncytot + 1

c..   update particle densities, time elapsed
         sumx = 0
         sumz = 0
         eta  = 0
         do  j = 1, itot
            y(j) = y(j) + b(j)
c..   positive or zero
            y(j) = dmax1( y(j), 0.0d0 )
            x(j) = y(j) * dble( nz(j) + nn(j) )
            sumx = sumx + x(j)
            sumz = sumz + y(j)* dble( nz(j) )
            eta  = eta  + y(j)* dble(nn(j) - nz(j))
         enddo
         checksum = sumx - 1.0d0

c..   calculate energy generation rate using mass excesses
c..   DNE is 3/2 kT for new particles
         eb   = 0
         dn   = 0
         do j = 1, itot
            fak =  y(j) - yold(j)
            eb  = eb - qq(j)*fak
            dn  = dn +       fak
         enddo
         eb      = eb *   9.65d+17    /dtstar
         dne      = dn * 1.2476d+17 *t9/dtstar
         aeps(k) = eb - dne

      else
c..   itbrn(k) = 1 .....................................

         if( dlibeb .eq. 0 )then

c..   default burns dlibeb
            call loburn(dtstar,k,it,kk,leq)

         else

c..   move old value for update below
            do j = 1, ndim
               y(j) = aex(j,k)
            enddo

         endif
c..   loburn matrix solutions not added to ncytot

      endif


      if( leq .eq. 0 )then
c..   only update abundances y(j) to aex(j,k) for leq=0
c..   this is for calls which only require aeps(k)-->s(5,k) in burn.f
         do j = 1, itot
            aex(j,k) = y(j)
         enddo
         aex(ndim,k) = 0.0d0
         do j = 1, itot
            aex(ndim,k) = aex(ndim,k)
     1           + aex(j,k)*dble( nz(j) )
         enddo
      endif

c..   add base and finite difference states
      ncyc(k) = ncycle + ncyc(k)
      if( nucleu .gt. 0 .and. nucleu .le. itot )then
         xidk(k) = xid(nucleu)
         nxid(k) = nucleu
         bk(k) = b(nucleu)
         yk(k) = y(nucleu)
      else
         xidk(k) = '     '
         nxid(k) = 0
         bk(k) = 0.0d0
         yk(k) = 0.0d0
      endif


      icall = 1
      return

 101  continue
      write(*,*)'more than ',ncymax,' subcycles'
      write(*,*)ncycle,' subcycles in zone ',k
      write(*,'(a5,8a12)')'k','t9','rho','dtstar','dtleft','dth',
     1     'checksum','eta','sumd'
      write(*,'(i5,1p8e12.3)') k,t9,rho,dtstar,dtleft,dth,
     1     checksum,eta,sumd
      if( nucleu .gt. 0 .and. nucleu .le. itot )then
         write(*,*)'fastest nucleus is ',xid(nucleu)
      endif

      stop'solven too many subcycles'

      end


