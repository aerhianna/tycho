      subroutine sqse(t9,rho,yeq,enc,k)

c..   QUASISTATISTICAL EQUILIBRIUM (modified: wda 11/19/04)
c..   input
c     itot = total number of n,p,alpha,nuclei
c     inuc = itot-3 = number of "nuclei"
c     idebug
c     t9,rho,eta
c     uqse  = shift in chemical potential for heavy nuclei
c..   output
c     niter = number of iterations used
c     uaaa  = 2*( mu(n) + mu(p) ) 
c     uhat  =     mu(n) - mu(p)
c     ierr  = flag for nonconvergent return

      implicit none

      include 'dimenfile'
      include 'crate'
      include 'cburn'
      include 'comcsolve'
      include 'caeps'
      include 'ceoset'

c..   u = chemical potential, yeq = nse mole number,
c..   xeq = nucleon fraction, 
c..   qeq is used for alpha net (fakes mu(p) = mu(n) )

      real*8  u(nnuc), yeq(ndim), xeq(ndim), qeq(nnuc)
      real*8 yeq0(ndim)
c      character*5 plus, goes
c      character*5 blank

c..   nloop = number of loops in iteration
c..   delmx = maximum fractional step in interation variable
c..   tiny  = convergence needed in nucleon, charge conservation
      integer*4 nloop
      real*8    tiny
      parameter( nloop = 500, tiny = 1.0d-13 )

c..   niso = number of isotopes in FKT list
c      integer*4 niso
c      parameter(niso = 806) 
c      real*8 p(7)
c      character*4 lkh                                     
c      character*5 nam(6),inam(niso),rnam(6)                      
c      character*1 vw,nr

c..   cxid = character id for isotope for summary i/o
c..   xxeq = nucleon fraction in nse
      character*5 cxid(nnuc)
      real*8      xxeq(nnuc), xxeq0(nnuc)
      
c      data blank/'     '/
c      data plus /'  +  '/,goes/' <-> '/
      integer*4 icall,iachain
      data icall/0/,iachain/1/
      real*8 delmx
      data delmx/0.2d0/

      integer*4 ib,in,ip,ia,ni56,nc12,nsi28,nfe54,pivot
      integer*4 niter,ierr,nucn,nucpr,nuc0,iqse
      integer*4 i,k,iitot,numa,na,n
      real*8 bpern,t9,rho,fact,rho9,t932,theta,tk,xal,ualow
      real*8 aqse,rat,x28,term,uqse,uhlow,uhat,uaaa,uhat0,u0,aa,zz
      real*8 up,un,yheold,maxx,ww(ndim),aansi28
      real*8 ueff,fuma,fex,arat,omxa,uhhi,xeqm,dedmu,etaq
      real*8 dxdmh,dedmh,dudmu,dydmu,dudmh,det,errx,yeeq,errz
      real*8 anum,hnum,ddu,delua,ddh,deluh,dfracu,dfrach,dfrac,df0
      real*8 reduct,enc,pnc,gam4,xpp,xnn,xaa,ynn,ypp,xxl,dxdmu,dydmh
      real*8 deluq,adeluq

      character*7  filein
      logical tobe

c..   icall flags first call for initiation
c..   iachain flags lack of N .ne. Z nuclei
c..   ib is index of most bound nucleus
c..   bpern is binding energy of most bound nucleus
c..   ni56, nc12 are indices for special nuclei
      save icall,iachain,ib,bpern
      save in,ip,ia,ni56,nc12,nsi28,nfe54
      save delmx
c---------------------------------------------------------------------
      write(*,*)'ENTERING QSE',icall

      ierr = 0

      xal = 0.0d0
ccccccccccccccccccc

      if( icall .eq. 0 )then
c..   only on initial call
         icall = 1
c..   save initial values
         do i = 1, nnuc
            yeq0(i) = yeq(i)
         enddo
c..   determine if network is alpha-chain
         nucn = 0
         nucpr = 0
         nuc0 = 0
         pivot = 0
         maxx = 0.0d0
c..count neutron-rich, Z=N, and neutron-poor nuclei
         do i = 1, inuc
            if( nn(i) .gt. nz(i) )then
               nucn = nucn + 1
            elseif( nn(i) .lt. nz(i) )then
               nucpr = nucpr + 1
            else
               nuc0 = nuc0 + 1
            endif
            if(( yeq(i)*dble(nn(i)+nz(i)) ) .gt. maxx)then
c..   find most abundant nucleus
               maxx = (yeq(i)*dble(nn(i)+nz(i)))
               pivot = i
            endif
         enddo
         write(*,'(4(a15,i5))')"n-rich nuclei",nucn,
     1        "p-rich nuclei",nucpr,
     2        "Z=N nuclei",nuc0,
     3        "total nuclei",inuc

         if( inuc .eq. nuc0 )then
c..   true, it is alpha chain
            iachain = 0
            filein = 'data.a'
            inquire(file=filein,exist=tobe)
            if( tobe )then
               open(33,file=filein,form='formatted',status='old')
            else
               write(*,*)'sqse: no file ',filein,' in directory'
               stop'sqse: data.a'
            endif
c            open(33,file='data.a')
         else
c..   false, it is general network
            iachain = 1
         endif

c..   define "particles"
         ia   = itot
         ip   = itot - 1
         in   = itot - 2

c..   find most bound nucleus in network
c..   start with He4
         bpern = (-qq(ia) + qq(ip)*2.0 + qq(in)*2.0)/4.0d0
c..   find qse pivotal nucleus
         iqse = 0
         do i = 1, itot
            fact = (-qq(i) + qq(ip)*dble( nz(i) ) 
     1           + qq(in)*dble( nn(i) ) )
     2           /float( nz(i) + nn(i) ) 

c..   remember index of special nuclei
            if( nz(i) .eq. 28  .and. nn(i) .eq. 28  )ni56  = i
            if( nz(i) .eq.  6  .and. nn(i) .eq.  6  )nc12  = i
            if( nz(i) .eq. 26  .and. nn(i) .eq. 28  )nfe54 = i
            if( nz(i) .eq. 14  .and. nn(i) .eq. 14  )nsi28 = i
            if( fact .gt. bpern )then
               bpern = fact
               ib = i
            endif
         enddo

         write(*,'(i5,a5,f9.4,a60)')ib,xid(ib),bpern,
     1        " excess Mev/nucleon; most bound nucleus in network"

c         write(*,'(a8,a5,a25/)')
c     1        "Using ", xid(pivot), " as pivotal QSE nucleus"

      endif

c-----------------------------------------------------------
c..   executed once on each call
      
      rho9  = rho * 1.0d-9
      t932  = t9 * sqrt( t9 )
      theta = 0.1013d0 * rho9 /t932
      tk    = t9/11.605d0
c..   neutron excess eta: ye    = (1.0d0 - eta )*0.5d0
      eta = 1.0d0 - 2.0d0*ye(k)   
c..   nuclear partition functions; ground state times enhancement  
      do i=1,itot
         ww(i)   = w(i) * pf(i)
         aa      = dble( nn(i) + nz(i) )
         aansi28 = dble( nn(nsi28) + nz(nsi28) )
         if( yeq(i)*dble(nn(i)+nz(i)) .gt. maxx )then
               maxx  = yeq(i)*dble( nn(i) + nz(i) )
               pivot = i
         endif
      enddo

cccccccccccccccccccccccccccccccccc
      write(*,'(2i5,a5,1p8e12.3)')iachain,pivot,xid(pivot),maxx
c..overwriting here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      pivot = nsi28

      write(*,'(a8,a5,a25,a5,1pe12.3/)')
     1     "Using ", xid(pivot), " as pivotal QSE nucleus" ,
     2     'X=',yeq(pivot)*(nn(pivot)+nz(pivot))
c-------------------------------------------------
      if( iachain .ne. 0 )then
c-------------------------------------------------
         if( ni56 .eq. 0 )then
            write(*,*)'error in qse, no Ni56 in network'
            stop
         endif
         if( nfe54 .eq. 0 )then
            write(*,*)'error in qse, no Fe54 in network'
            stop
         endif
c..   QUASISTATISTICAL EQUILIBRIUM
         yheold = yeq(nnuc)

c      yeq(nnuc) = (yeq(ni56)/yeq(nsi28))**(1.0d0/7.0d0) * 8.0d0 * theta/
c     1    2.0d0**(3.0d0/14.0d0)*dexp(-1.0d0/7.0d0*49.382d0/tk)+yeq(nnuc)

         xal = yeq(nnuc) * dble(nz(nnuc) + nn(nnuc))

c..   try dissociated value for Yn/Yp for first uhat guess
c..   initial value for mu(alpha) = u(ia)
c         u(ia) = qq(itot) + tk*dlog( xal*theta/32.0d0 )
ccccccccccccccccccccccccccccccccc
         u(ia) = qq(nnuc) + tk*etanuc(nnuc,k)
         ualow = u(ia)

c..   first guess for chemical potential of basis nucleus
         aqse = nz(pivot) + nn(pivot)
c..   ratio of ni56/si28 by mass
c         rat     = qq(ni56)/qq(nsi28) * 2.0d0**2.5
c     1        * exp( ( qq(nsi28) + 7.0d0*u(ia) - qq(ni56) )/tk )
c         x28     = (1.0d0 - xal)/(1.0d0 + rat)
c         term    = x28 * theta / ( ww(nsi28) * aqse**2.5d0 )
c         u(nsi28)= qq(nsi28) + tk*dlog( term )
c         u(nsi28) = qq(nsi28) + tk*(etanuc(nsi28,k)+dlog(ww(nsi28)))

c         u(pivot) = qq(pivot) + tk*(etanuc(pivot,k)-dlog(ww(pivot)))
         u(pivot) = qq(pivot) + tk*etanuc(pivot,k)
         uqse     = u(pivot)

c..   initial guess for muhat = mu(n) - mu(p)
c..   high T value = free n have all neutron excess
         omxa  = 1.0d0 
         uhhi  = ( qq(in) - qq(ip) )
     1        + tk*dlog( (omxa + eta)/(omxa - eta) )
c..   low T value = Fe54 has all neutron excess
         uhlow = qq(nfe54) + tk*dlog( eta*theta*0.5d0/54.0d0**1.5 )
     1        - uqse     - 26.0d0/4.0d0*ualow
c..   compromise value
         xal = xal
         uhat  = uhlow*(1.0d0 - xal) + uhhi*xal
c..   use chemical potential values from state.f
         up   = qq(nnuc-1) + tk*etanuc(nnuc-1,k)
         un   = qq(nnuc-2) + tk*etanuc(nnuc-2,k)
         uhat = (un - up)/1.0d0
         uaaa = 2.0*(un + up)/1.0d0
c..   uses u(ia) = u + m here, so it is "mu" not "u"
c         uaaa = u(ia)
c..   save first guesses for convergence diagnostics
         uhat0 = uhat
         u0 = u(pivot)

         write(*,'(a30,1p8e12.3)')
     1  "quess: ua, uaaa, uhat, uqse ",u(ia),uaaa,uhat, uqse
ccccccccccccccccccccccccccc

c..   uses u(nsi28) = u + m here
c..   save guess for convergence diagnostics
c         u0 = u(nsi28)
c..   interate for nucleon number (sum Xi = 1)........
         do n = 1, nloop
c..   generate trial nse abundance values
            do i = 1, itot  
               if( i .eq. ia )then
                  aa     =  nz(i) + nn(i) 
                  zz     =  nz(i) 
                  ueff   = uaaa - qq(i)
                  arat   = ww(i) / theta
                  yeq(i) = xal/aa
                  xeq(i) = xal
c     yeq(i) = arat * exp( ueff / tk )
c     xeq(i) = yeq(i)*aa
                  fuma   = ( aa - aqse )*0.25d0
                  fex    = 0
                  u(i)   = qq(i) + ueff
               elseif( i .eq. ip )then
                  aa     =  nz(i) 
                  fuma   = 0
                  fex    = -1
                  arat   = ww(i) / theta
                  ueff   = 0.25d0 * uaaa
     1                 - 0.50d0 * uhat
     2                 - qq(i)
                  yeq(i) = arat * exp( ueff / tk )
                  xeq(i) = yeq(i)
                  u(i)   = qq(i) + ueff
               elseif( i .eq. in )then
                  aa     = nn(i) 
                  fuma   = 0
                  fex    = 1
                  arat   = ww(i) / theta
                  ueff   = 0.25d0 * uaaa
     1                 + 0.50d0 * uhat
     2                 - qq(i)
                  yeq(i) = arat * exp( ueff / tk )
                  xeq(i) = yeq(i)
                  u(i)   = qq(i) + ueff
c     aa     =  nz(i) + nn(i) 
c     fuma   = ( aa )*0.25d0
c     fex    = nn(i) - nz(i)
c     arat   = aa**1.5d0 * w(i) / theta
c     ueff   =  fuma * uaaa
c     1                    + fex *0.5d0* uhat
c     2                    - qq(i)
c     yeq(i) = arat * exp( ueff / tk )
c     xeq(i) = yeq(i)*aa
c     u(i)   = qq(i) + ueff
               else
                  aa     =  nz(i) + nn(i) 
                  fuma   = ( aa - aqse )*0.25d0
                  fex    = nn(i) - nz(i)
                  arat   = aa**1.5d0 * ww(i) / theta
                  ueff   = u(nsi28)     + fuma * uaaa
     1                 + fex *0.5d0* uhat
     2                 - qq(i)
                  yeq(i) = arat * exp( ueff / tk )
                  xeq(i) = yeq(i)*aa
                  u(i)   = qq(i) + ueff
               endif
            enddo
c..   generate Ye and sum Xi, and their derivatives for Newton-Raphson
            xeqm    = 0
            dxdmu   = 0
            dedmu   = 0
            etaq    = 0
            dxdmh   = 0
            dedmh   = 0
            do i = 1, itot  
               if( i .ne. pivot )then
c..   complete sum for all nucleons
                  if( nz(i) .gt. 7 )then
c..   z >= 6
                     aa    = nz(i) + nn(i) 
                     xeqm  = xeqm + aa * yeq(i)
                     dudmu = ( aa - aqse )*0.25d0 
                     dydmu = yeq(i) / tk * dudmu
                     dxdmu = dxdmu + aa*dydmu
                     dudmh = dble( nn(i) - nz(i) )*0.5d0
                     dydmh = yeq(i) / tk * dudmh
                     dxdmh = dxdmh + aa*dydmh
                     dedmu = dedmu + dble( nn(i) - nz(i) )*dydmu
                     dedmh = dedmh + dble( nn(i) - nz(i) )*dydmh
                     etaq  = etaq  + dble( nn(i) - nz(i) )*yeq(i)
                  elseif(i .eq. ip .or. i .eq. in .or. i .eq. ia)then
c..   z < 6
                     aa    = nz(i) + nn(i) 
                     xeqm  = xeqm + aa * yeq(i)
                     dudmu = aa*0.25d0
                     dydmu = yeq(i) / tk * dudmu
                     dxdmu = dxdmu + aa*dydmu
                     dudmh = dble( nn(i) - nz(i) )*0.5d0
                     dydmh = yeq(i) / tk * dudmh
                     dxdmh = dxdmh + aa*dydmh
                     dedmu = dedmu + dble( nn(i) - nz(i) )*dydmu
                     dedmh = dedmh + dble( nn(i) - nz(i) )*dydmh
                     etaq  = etaq  + dble( nn(i) - nz(i) )*yeq(i)
                  else
                     aa    = nz(i) + nn(i) 
                     xeqm  = xeqm  + aa * yeq(i)
                     etaq  = etaq  + dble( nn(i) - nz(i) )*yeq(i)
                  endif
               else
                  aa    = nz(i) + nn(i) 
                  xeqm  = xeqm  + aa * yeq(i)
                  etaq  = etaq  + dble( nn(i) - nz(i) )*yeq(i)
               endif
            enddo 
            det    = dxdmu*dedmh - dxdmh*dedmu
c..   safety check for redundant equations:
c..   attempts to avoid unrealistically large corrections
            if( n .eq. 1 )then
               if( abs(det) .le. 1.0d-1 )then
                  det = dsign( 1.0d-1, det ) 
               endif
            endif
            errx   = xeqm - 1.0d0
            yeeq   = (1.0d0 - etaq)*0.5d0
            errz   = etaq - eta
            anum   = -errx*dedmh + errz*dxdmh
            hnum   =  errx*dedmu - errz*dxdmu
            ddu    = anum/det
            delua  = -errx/dxdmu
            ddh    = hnum/det
            deluh  = -errz/dedmh
c..   safety check for large changes:
c..   reduces increments simultaneously
c..   uses n*kT as measure of reasonable step size
            dfracu = abs( ddu )
            dfrach = abs( ddh )
            dfrac  = dmax1( dfracu, dfrach )
            df0    = 3.0d0*tk
            if( dfrac .gt. df0 )then
               reduct = df0/dfrac
               ddu = ddu * reduct
               ddh = ddh * reduct
            else
               reduct = 1.0d0
            endif
            uaaa  = uaaa + ddu
            uhat  = uhat + ddh
c            u(nsi28) = uqse
            u(pivot) = uqse
            if( abs(xeqm - 1.0d0) .lt. tiny )goto 998
         enddo

         write(*,'(a33,i4,4(a8,1pe11.3))')
     1        "NONCONVERGENCE in QSE: iteration ",n-1,
     2        " d(X)", xeqm-1,  " d(eta)", errz,
     3        " mu(a)",   uaaa,  " mu(hat)",  uhat
         ierr  = 1
         niter = n
         return

 998     continue
c...................................end of iteration logic......
         write(*,'(a30,1p4e12.3,a12,i5)')
     1  "final: ua, uaaa, uhat, uqse ",u(ia),uaaa,uhat, uqse,
     2        'iterates',n
c..   QSE values are now determined
         niter = n
c..   uses mass excesses for energy calculation....
c..   energy relative to C12 nuclei
         enc  = 0
         pnc  = 0
         yeeq = 0
         do i = 1, itot
c            write(*,*)i,yeq(i)/yeq(nsi28),dexp((u(i)-ww(i))-
c     1            (u(nsi28)-ww(nsi28)))
            enc  = enc  + yeq(i)*(1.2476d+17*t9 +  qq(i)*9.65d+17 )
            pnc  = pnc  + yeq(i)* 1.2476d+17*t9 * 1.5d0
            yeeq = yeeq + dble( nz(i) )*yeq(i)
         enddo
         pnc  = pnc * rho
         gam4 = pnc/(enc * rho) + 1.0d0
         yeq(ndim) = yeeq
         
         write(*,'(a5,3a11,a15,7a11)')'n','t9','rho','eta','enc',
     1        'xeq(ia)','uaaa','uhat','uqse','det','yeeq','gam4'
         write(*,'(i5,1p3e11.2,1pe15.7,1p7e11.2)')n,t9,rho,eta,enc,
     1        xeq(ia),uaaa,uhat,uqse,det,yeeq,gam4

c..   sum n-rich,alpha,p-rich isotopes
         xpp = 0
         xnn = 0
         xaa = 0
         ynn = 0
         ypp = 0
         do i = 1, itot
            aa = nz(i) + nn(i)
            zz = nz(i)
            if( nz(i) .eq. nn(i) )then
               xaa = xaa + xeq(i)
            elseif( nz(i) .lt. nn(i) )then
               xnn = xnn + xeq(i)
               ynn = ynn + yeq(i)* dble( nn(i) - nz(i) )
            else
               xpp = xpp + xeq(i)
               ypp = ypp + yeq(i)* dble( nz(i) - nn(i) )
            endif
         enddo
c..   most abundant nuclei (above xxl)
         iitot = 0
         xxl   = 1.0d-3
         do i = 1, itot-3
            if( xeq(i) .gt. xxl )then
               iitot        = iitot + 1
               cxid(iitot)  = xid(i)
               xxeq(iitot)  = xeq(i)
               xxeq0(iitot) = yeq0(i)*dble( nz(i) + nn(i) )
            endif
         enddo
         if( iitot .gt. 0 )
     1        write(*,'(5(a5,1p2e10.2,2x))')(cxid(i),xxeq(i),
     1        xxeq0(i), i=1,iitot)

c..   n,p,alphas
         write(*,'(3(a5,1p2e10.2,2x))')(xid(i),xeq(i),
     1        yeq0(i)*dble( nz(i) + nn(i)), i=itot-2,itot)
         write(*,*)' '

c         write(33,'(2i5,1p13e11.3)')ierr,niter,t9,rho,eta,
c     1        uaaa,uaaa0,uhat,uhat0,xeq(itot-1),xeq(itot-2),xeq(itot)
c     2        ,tk,enc*1.0d-18,pnc/rho*1.0d-18

c         do i = 1, itot
c            write(34,'(3i5,1p16e11.3)')nz(i),nn(i),nz(i)+nn(i),xeq(i)
c         enddo
c         do i = 1, itot
c            aa = nz(i) + nn(i)
c            zz = nz(i)
c            zn = nn(i)
c            utst0 = uaaa/4.0d0*aa + (zn - zz)/2.0d0*uhat
c           utst1 = u(nsi28) + uaaa/4.0d0*(aa -28.0d0)
c     1           + (zn - zz)/2.0d0*uhat
c            write(35,'(3i5,1p8e12.3)')nz(i),nn(i),nz(i)+nn(i),xeq(i),
c     1           u(i),utst0,utst1,u(i)-utst0,u(i)-utst1
c         enddo

         return

c-------------------------------
      else
c-------------------------------
         do i = 1, itot
            qeq(i) = qq(i)
         enddo
c..   represent n,p as "average nucleons" so Yp = Yn
         qeq(in) = 0.5*( qq(in) + qq(ip) )
         qeq(ip) = qeq(in)


c..   alpha chain quasiequilibrium with X(He4) fixed

c..   fixed value for mu(alpha)  
         uaaa = qeq(ia) + tk*dlog( xal*theta/32.0d0 )
c..   compromise value
         u(ia) = uaaa

c..   first guess for chemical potential of basis nucleus
c..   take si28 as basis (MUST REWRITE FOR OTHER CHOICE)
         iqse    = nsi28
         aqse    = nz(iqse) + nn(iqse)
c..   ratio of ni56/si28 by mass
         rat     = qq(ni56)/qq(nsi28) * 2.0d0**2.5
     1        * exp( ( qq(nsi28) + 7.0d0*uaaa - qq(ni56) )/tk )
         x28     = (1.0d0 - xal)/(1.0d0 + rat)
         term    =  x28 * theta / ( ww(iqse) * aqse**2.5d0 )
         u(iqse) =  qeq(iqse) + tk*dlog( term ) 

         write(*,*)'ia,iqse,aqse,term,theta,u(iqse),u(ia)'

         write(*,'(2i5,1p8e10.2)')ia,iqse,aqse,term,theta,u(iqse),u(ia)


c..   uses u(nsi28) = u + m here

c..   save guess for convergence diagnostics
         u0 = u(iqse)

c..   interate for nucleon number (sum Xi = 1).....................

         do n = 1, nloop
c..   generate trial yeq values
            do i = 1, itot  
               if( i .eq. ip .or. i .eq. in )then
                  xeq(i) = 0
                  aa     = 1
                  yeq(i) = xeq(i)/aa
                  numa   = 0
                  na = aa
               elseif( i .eq. ia )then
                  aa     = nz(i) + nn(i)
                  numa   = ( aa - aqse )*0.25d0
                  xeq(i) = xal
                  yeq(i) = xeq(i)/aa
                  arat   = (aa)**1.5d0 * ww(i) / theta
                  ueff   = tk * dlog( yeq(i)/arat )
c     yeq(i) = arat * exp( ueff / tk )
                  na = aa
               else
                  aa     = nz(i) + nn(i)
                  numa   = ( aa - aqse )*0.25d0
                  arat   = (aa)**1.5d0 * ww(i) / theta
                  ueff   = u(iqse) + dble(numa)*uaaa - qeq(i)
                  yeq(i) = arat * exp( ueff / tk )
                  xeq(i) = yeq(i)*aa
                  na = aa
               endif
            enddo

            xeqm  = -1.0d0
            dxdmu = 0
            do i = 1, itot  
               if( i .ne. ia )then
                  aa    =  dble( nz(i) + nn(i) )
                  xeqm  = xeqm + aa * yeq(i)
                  dudmu = 1.0d0
                  dydmu = yeq(i) / tk * dudmu
                  dxdmu = dxdmu + aa*dydmu
               else
                  aa    =  dble( nz(i) + nn(i) )
                  xeqm  = xeqm + aa * yeq(i)
               endif
            enddo             
            delmx = 3.0d0*tk
            if( dxdmu .ne. 0.0d0 )then
               deluq  = - xeqm / dxdmu
c..   limit size of change in mu(alpha)
               adeluq  = abs( deluq )
               if( adeluq .gt. delmx )then
                  deluq = deluq/adeluq* delmx
               endif
               u(nsi28) = u(nsi28) + deluq
            else
               write(*,*)' error: dxdmu = ',dxdmu
               stop
            endif
            if( n .eq. 1 )write(*,*)'n,xeqm,dxdmu,del,',
     1           'adel,u(si28),u0,delmx,X(a),si28,ni56'
            write(*,'(i5,1p13e10.2)')n,xeqm,dxdmu,deluq,
     1           adeluq,u(nsi28),u0,delmx,xeq(ia),xeq(nsi28),xeq(ni56)

            if( abs(xeqm) .lt. tiny )goto 999
         enddo
         write(*,'(a25,i5,2(a10,1pe12.3))')
     1        "nonconvergence (equil), n =",n," mu(a)=",uaaa,
     1        "xeq - 1",xeqm
         ierr  = 1
         niter = n
         return

 999     continue
c...................................end of iteration logic......
c..   QSE values are now determined
         niter = n
         write(*,*)'QSE values determined in ',n,' iterates'
c..   uses mass excesses for energy calculation....
c..   set zero energy to cold, pure most bound nucleus in network
c     enc =  bpern * 9.65d+17
c..   relative to C12 nuclei
         enc = 0
         pnc = 0
         do i = 1, itot
            enc = enc + yeq(i)*( 1.2476d+17*t9 + qeq(i)*9.65d+17 )
            pnc = pnc + yeq(i)*  1.2476d+17*t9 *1.5d0 
         enddo
         pnc  = pnc * rho
         gam4 = pnc/(enc*rho) + 1.0d0

         write(*,'(i5,1p13e10.2)')n,t9,rho,enc*1.0d-18,xeqm,
     1        xeq(ia),xeq(ni56),xeq(ip),uaaa,pnc*1.0d-18/rho,u0

         write(33,'(i5,1p17e10.2)')n,t9,rho,uaaa
     1        ,(xeq(i),i=1,13),xeq(itot)

         write(*,'(6(a5,1pe10.2))') (xid(i),xeq(i),i=1,itot)
         write(*,*)' '

c------------------------
      endif
c------------------------

      return
      end
