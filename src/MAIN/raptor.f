
      subroutine raptor(m)

c..   finds newly rezoned regions with unphysical gradients
c..   adjusts pairwise the abundances about the boundary to
c..   reduce the gradient to a realistic value

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'czone'
      include 'cburn'
      include 'cconst'
      include 'cgtintp'

      real*8 scrx(ndim,kdm)
      real*8 x1(ndim),x2(ndim),fmix,fmx1,fmx2,dn1,dn2,dnm,fmax
      real*8 delf,fmid,fsing
      real*8 delt,delp,tbar,pbar,ybar
      real*8 delt1,delp1
      real*8 beta,doufak,dbar,vbar,dvdt,dltdlp
      real*8 akbar,tledd,tb4,ombeta,tle,tlum
      real*8 xx(2),yy(2),yi,dfmix,ofmix

      integer*4 krat, m, k, n, i, loop
c--------------------------------------------------------------
      krat = 0

c..save original abundances for caution
      do k = 2, kk
         do n = 1, nnuc
            scrx(n,k) = x(n,k)
         enddo
      enddo

c..find a bad zone
      do k = 2, kk
         delp = p(m,k+1)-p(m,k)

         if( dnab(k) .gt. 1.0d0 .or. dnab(k) .lt. -1.0d0 
     2        .and. iadd(k) .ne. 0 
     3        .and. doux(k) .gt. 1.0d-3 
     4        .or. delp .gt. 0.0d0 )then
c..   useless to mix if doux = 0 (zones already mixed)
c..   newly added zones have interpolated abundance which may be poor
c..   only pick zones with large nabla as needing improvement
c..
            krat = k
            write(*,'(a16,4i6,1p8e13.5)')'RAPTOR strikes: ',
     1           k,krat,iadd(k),ic(k),xm(k)/sol,dnab(k),doux(k),delp

c            write(*,'(9(a8,0pf10.3))')'oldnab',dnab(krat),
c     1           'ad',dnad(krat),
c     1           'nrad',dnrad(krat),'doux',doux(krat),
c     2           'r-a-x',dnrad(krat)-dnad(krat)-doux(krat),
c     3           'n-a-x',dnab(krat)-dnad(krat)-doux(krat),
c     4           'n-r',dnab(krat)-dnrad(krat),
c     5           'dx(he)',xa(nnuc)*(x(nnuc,krat+1)-x(nnuc,krat)),
c     6           'odx(he)',xa(nnuc)*(xold(nnuc,krat+1)-xold(nnuc,krat))

c..initial guess at limits of physical region (dP < 0) for HSE
c..   fmix = mixing parameter
c..   fmix = 0   gives old value
c..   fmix = fmax gives complete mixing (uniform x)
c..   complete mixing of unequal mass zones
c..   choose smaller value to avoid overshooting complete mix
            if( dmh(krat) .lt. dmh(krat+1) )then
               fmax = dmh(krat)/( dmh(krat+1) + dmh(krat) )
            else
               fmax =  dmh(krat+1)/( dmh(krat+1) + dmh(krat) )
            endif
c..this assumes only two zones in resevoir
            fmax = 0.5d0
cccccccccccccccccccccccccccccccccccccccccccccccccccc

c..   (1-fmax) is coef. of krat+1 value
            fmx1 = 0.0d0
            fmx2 = fmax
c..   a bad zone boundary found at krat
c..   interpolation arrays
            do n = 1, nnuc+1
               x1(n) = x(n,krat)
               x2(n) = x(n,krat+1)
            enddo

c..   look for bounding values
c..   no mixing
            fmix = fmx1
            do n = 1, nnuc+1
               x(n,krat)   = fmix*x2(n) + (1.0d0-fmix)*x1(n)
               x(n,krat+1) = fmix*x1(n) + (1.0d0-fmix)*x2(n)
            enddo

            call state(krat,krat+1,m)

            delt = t(m,krat+1)-t(m,krat)
            delp = p(m,krat+1)-p(m,krat)
            
            tbar = (t(m,krat+1) + t(m,krat))*0.5d0
            pbar = (p(m,krat+1) + p(m,krat))*0.5d0
            dn1  = delt/tbar*pbar/delp
c            write(*,'(5x,8a12)')'fmix','delt','delp','nabla',
c     1           'delp','delp+','delp-'
c            write(*,'(5x,1p8e12.3)')fmix,delt,delp,dn1,delp,
c     1           p(m,krat+2)-p(m,krat+1),p(m,krat)-p(m,krat-1)
c..   complete mixing
            fmix = fmx2
            do n = 1, nnuc+1
               x(n,krat)   = fmix*x2(n) + (1.0d0-fmix)*x1(n)
               x(n,krat+1) = fmix*x1(n) + (1.0d0-fmix)*x2(n)
            enddo

            call state(krat,krat+1,m)

            delt1 = t(m,krat+1)-t(m,krat)
            delp1 = p(m,krat+1)-p(m,krat)
            
            tbar = (t(m,krat+1) + t(m,krat))*0.5d0
            pbar = (p(m,krat+1) + p(m,krat))*0.5d0
            dn2  = delt1/tbar*pbar/delp1
c            write(*,'(5x,1p8e12.3)')fmix,delt1,delp1,dn2,delp1,
c     1           p(m,krat+2)-p(m,krat+1),p(m,krat)-p(m,krat-1)

            if( delp .ge. 0.0d0 .and. delp1 .lt. 0.0d0 )then
c..do bisection to find fmix value (fsing) for singularity in nabla
               delf = fmx2 - fmx1
               fmid = fmx1
               loop = 0
 300           continue
               loop = loop + 1

               if( loop .gt. 100 )stop'cc'

               delf = delf * 0.5d0
               fmix = fmid + delf
               
               do n = 1, nnuc+1
                  x(n,krat)   = fmix*x2(n) + (1.0d0-fmix)*x1(n)
                  x(n,krat+1) = fmix*x1(n) + (1.0d0-fmix)*x2(n)
               enddo

               call state(krat,krat+1,m)

               delt = t(m,krat+1)-t(m,krat)
               delp = p(m,krat+1)-p(m,krat)
            
               tbar = (t(m,krat+1) + t(m,krat))*0.5d0
               pbar = (p(m,krat+1) + p(m,krat))*0.5d0
               dnm  = delt/tbar*pbar/delp
c               write(*,'(i5,1p8e12.3)')loop,fmix,delt,delp,dnm,delf,fmid
c            write(*,'(i5,1p8e12.3)')loop,fmix,delt,delp,dnm,delp,
c     1           p(m,krat+2)-p(m,krat+1),p(m,krat)-p(m,krat-1)

               if( delp .gt. 0.0d0 )then
                  fmid = fmix
               endif

               if( delf .lt. 1.0d-4 )then
c..1.0e-4 give nabla of order 1.0e4 which is probably big enough
                  write(*,*)'singular point is fmix = ',fmix
c..   insure negative delp, so increase fmix slightly to be on the
c..   correct side of the singularity
                  fsing = fmix * 1.0002
                  fmix = fsing
                  goto 301
               endif

               goto 300
c...................................................................
            elseif( delp .lt. 0.0d0 .and. delp1 .lt. 0.0d0 )then
c..both are hydrostatic (HSE)
c..check if nabla is in range
c               write(*,*)'both are hydrostatic (HSE), no singular point'
c               write(*,'(4(a10,1pe11.3))')
c     1              'fmx', fmx1,'delp',delp,'nabla',dn1
c               write(*,'(4(a10,1pe11.3))')
c     1              'fmx2',fmx2,'delp',delp1,'nabla',dn2
c..choose plausible first interpolant
               fmix = 0.5d0*( fmx1 + fmx2 )
c               fmix = fmx1 + (fmx2-fmx1)/(dn2-dn1)*( dnad(krat)-dn1)

c               write(*,*)'trying fmix = ',fmix
c               write(*,'(1p8e12.3)')dnad(krat),dn1,dn2,fmx1,fmx2

c               stop'ppp'
c.................................................................
            elseif( delp .lt. 0.0d0 .and. delp1 .ge. 0.0d0 )then
c..mixed value is NOT hydrostatic but unmixed is
               write(*,*)'mixed value NOT hydrostatic but unmixed is'
               write(*,'(4(a10,1pe11.3))')
     1              'fmx', fmx1,'delp',delp,'nabla',dn1
               write(*,'(4(a10,1pe11.3))')
     1              'fmx2',fmx2,'delp',delp1,'nabla',dn2
               stop'-+'
c.................................................................
            else
               write(*,*)'nonphysical values at boundaries'
               write(*,'(4(a10,1pe11.3))')
     1              'fmx', fmx1,'delp',delp,'nabla',dn1
               write(*,'(4(a10,1pe11.3))')
     1              'fmx2',fmx2,'delp',delp1,'nabla',dn2
c..   both values of NOT in HSE
c..   choose closest value and punt
c               stop'++'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




            endif
            
 301        continue

c..   finite difference Newton method, 
c..   iterate on fmix, the degree of mixing

            xx(1) = fmix
            yy(1) = 0.0d0

            do i = 1, 40
c..   tentative mixed abundances about boundary krat
               do n = 1, nnuc+1
                  x(n,krat)   = fmix*x2(n) + (1.0d0-fmix)*x1(n)
                  x(n,krat+1) = fmix*x1(n) + (1.0d0-fmix)*x2(n)
               enddo

               call state(krat,krat+1,m)

               delt = t(m,krat+1)-t(m,krat)
               delp = p(m,krat+1)-p(m,krat)

               tbar = (t(m,krat+1) + t(m,krat))*0.5d0
               pbar = (p(m,krat+1) + p(m,krat))*0.5d0

c..adiabatic gradient
               dvdt    = - et(krat)/( p(m,krat) + ev(krat) )
               dltdlp  =
     1              p(m,krat)/( t(m,krat)*( pt(krat) + pv(krat)*dvdt ))
               dnad(krat) = 0.5d0*dltdlp
               dvdt    = - et(krat+1)/( p(m,krat+1) + ev(krat+1) )
               dltdlp  =
     1              p(m,krat+1)/( t(m,krat+1)*
     2              ( pt(krat+1) + pv(krat+1)*dvdt ))
               dnad(krat) = dnad(krat) + 0.5d0*dltdlp

c..   desired value of nabla: neutral stability
               yi =  dnad(krat) + 0.4

c..   hydrostatic balance requires a negative pressure gradient, so
c..   we may have to jump a singularity in nabla to get to the
c..   physical solution
               if( delp .gt. 0.0d0 )then
                  write(*,'(a30,i5,1pe12.3)')
     1                 'RAPTOR: POSITIVE p GRADIENT',krat,delp
               endif

c..   new nabla and mix parameter used to get it
               xx(2) = fmix
               yy(2) = delt/tbar*pbar/delp

               if( i .ge. 2 )then
c..   expand about recent point
                  dfmix = (yi-yy(2))/(yy(1)-yy(2))*(xx(1)-xx(2))

c                  write(*,'(i3,10(a6,1pe10.2))')i,
c     1                 'x2',xx(2),'y2',yy(2),
c     1                 'x1',xx(1),'y1',yy(1),
c     1                 'yi-y2', yi-yy(2), 'y2-y1', yy(2)-yy(1),
c     3                 'x2-x1',(xx(2)-xx(1)),
c     4                 'dfmix',dfmix,
c     5                 'dp',delp,'dT',delt

               else
c..   first time through
                  dfmix = 0.01d0

c                  write(*,'(i3,9(a6,1pe10.2))')i,
c     1                 'x2',xx(2),'y2',yy(2),
c     1                 'x1',xx(1),'y1',yy(1),
c     1                 'yi-y2', yi-yy(2), 'y2-y1', yy(2)-yy(1),
c     3                 'x2-x1',(xx(2)-xx(1)),
c     4                 'dfmix',dfmix,'yi',yi
               endif

c            write(*,'(i5,1p8e12.3)')i,fmix,delt,xx(2),yy(2),delp,
c     1           p(m,krat+2)-p(m,krat+1),p(m,krat)-p(m,krat-1)


               if( abs( dfmix ) .le. 1.0d-8 )goto 200

               ofmix = fmix
               fmix  = fmix + dfmix

c..adjust fmix near boundaries of allowed range
               if( fmix .lt. 0.0d0 .or. 
     1              fmix .gt. fmax )then
c..   unmixing (fmix <0 )not allowed 
c..   (2nd law of thermodynamics is violated)
c..   reverse mixing (fmix > 1/2 ) inverts the sense of the
c..   abundance gradient, and is unphysical

c..bisection to approach physical limiting cases without overshoot
                  if( fmix .lt. 0.0d0 )then
                     fmix = 0.5d0*ofmix
                  endif
                  if( fmix .gt. fmax )then
                     fmix = ofmix + 0.5d0*( fmax - ofmix )
                  endif
                  dfmix = fmix-ofmix
c..bail out; near edge
                  if( abs( dfmix ) .lt. 1.0d-4 )goto 200

c                  write(*,'(a30,f8.3,a10,f8.3,a10,f8.3,a10,i5,a15,i5,
c     1                 a10,f8.3,a10,1pe11.3)')
c     2                 'RAPTOR: limits: dfmix',dfmix,
c     3                 'old fmix',ofmix,'fmax',fmax,
c     4                 'krat',krat,'iterations',i,
c     4                 'fmix-->',fmix, 'nabla',yy(2)

c                  write(*,'(4(a10,1pe12.3))')'delt',delt,'tbar',tbar,
c     1                 'delp',delp,'pbar',pbar
c                  write(*,'(4(a10,1pe12.3))')'Xa',x(nnuc,krat)*xa(nnuc),
c     1                 'Xa+',x(nnuc,krat+1)*xa(nnuc),
c     1                 'old Xa ',scrx(nnuc,krat)*xa(nnuc),
c     1                 'old Xa+',scrx(nnuc,krat+1)*xa(nnuc)

               endif

c..   set latest value to old value
               xx(1) = xx(2)
               yy(1) = yy(2)
            enddo
            write(*,*)'RAPTOR: nonconvergence in ',i,
     1           ' steps, krat ',krat
            stop'RAPTOR nonconvergence'
c...................................................................
c..   best values determined, so update
 200        continue

c..reset nablas
            dnab(krat) = yy(2)
c..   radiative gradient
            akbar  = 0.5d0*( ak(krat) + ak(krat+1) )
            tledd  = pi4 * crad * grav * xm(krat) / akbar
            tb4    = 0.5d0*( t(m,krat+1)**4 + t(m,krat)**4 )
            ombeta = arad*tb4/(3.0d0*pbar)
            tle    = tledd * ombeta * 4.0d0 *dnab(i)
            if( tl(m,krat) .le. tle )then
               tlum = tle
            else
               tlum = tl(m,krat)
            endif
            dnrad(krat) = tlum / ( tledd * ombeta * 4.0d0 )

c..   evaluate new doux
            dely(krat) = 0
            ybar = 0
c..   loop includes free electron contribution
            do n = 1, nnuc+1
               dely(krat) = dely(krat) + ( x(n,krat+1) - x(n,krat) )
               ybar = ybar + (  x(n,krat+1) + x(n,krat) )* 0.5d0
            enddo
c     dbar = 0.5d0*( 1.0d0/v(nc,i+1) + 1.0d0/v(nc,i))
c..   avoids overflow if v(nc,i+1) is zero (modes=0)
            vbar = 0.5d0*(v(m,krat+1) + v(m,krat))
            dbar = 1.0d0/vbar
            if( dely(krat) .ge. 1.0d-14 .and. 
     1           r(m,krat) .gt. 2.0d0*r(m,2) .and.
     2           r(m,1) .le. 1.0d0 )then
c..   no semiconvection in sphere to shell transition
c..   1 cm is zero to avoid roundoff error
c..   semiconvection has lighter gas on top
               beta = rgas*ybar*tbar*dbar/pbar
               doufak = beta/(4.0d0 - 3.0d0*beta)
               doux(krat) = doufak * dely(krat)/ybar * 
     1              pi4*r(m,krat)**4*0.5d0*( p(m,krat)+p(m,krat+1) )/
     1              ( grav*xm(krat)*dmi(krat) )
            else
               doux(krat) = 0
            endif


c            write(*,'(9(a8,0pf10.3))')'newnab',dnab(krat),
c     1           'ad',dnad(krat),
c     1           'nrad',dnrad(krat),'doux',doux(krat),
c     2           'r-a-x',dnrad(krat)-dnad(krat)-doux(krat),
c     3           'n-a-x',dnab(krat)-dnad(krat)-doux(krat),
c     4           'n-r',dnab(krat)-dnrad(krat),
c     5           'dx(he)',xa(nnuc)*(x(nnuc,krat+1)-x(nnuc,krat)),
c     6           'odx(he)',xa(nnuc)*(xold(nnuc,krat+1)-xold(nnuc,krat))

         endif
      enddo

c..index of bad zone is krat
c      if( krat .le. 0 )then
c         write(*,*)'RAPTOR: no bad zones sighted'
c      endif

ccccccccccccccccccccccccccccccccccccccccccccc
c      do k = 2, kk
c         delp = p(m,k+1)-p(m,k)
c         if( delp .gt. 0.0d0 )then
c            write(*,'(a8,2i5,1p9e11.3)')'LV Rapt',
c     1           k,m,delp
c            write(*,'(a8,2i5,1p9e11.3)')'LV Pm',
c     1           k,m,
c     2           p(m,k-1),p(m,k),p(m,k+1),p(m,k+2),
c     3           p(m,k+3),p(m,k+4)
c            write(*,'(a8,2i5,1p9e11.3)')'LV dPm',
c     1           k,m,
c     2           p(m,k)-p(m,k-1),   p(m,k+1)-p(m,k),  p(m,k+2)-p(m,k+1),
c     3           p(m,k+3)-p(m,k+2), p(m,k+4)-p(m,k+3),p(m,k+5)-p(m,k+4)

c            write(*,'(a8,2i5,1p9e11.3)')'LV P2',
c     1           k,m,
c     2           p(2,k-1),p(2,k),p(2,k+1),p(2,k+2),
c     3           p(2,k+3),p(2,k+4)
c            write(*,'(a8,2i5,1p9e11.3)')'LV dP2',
c     1           k,m,
c     2           p(2,k)-p(2,k-1),   p(2,k+1)-p(2,k),  p(2,k+2)-p(2,k+1),
c     3           p(2,k+3)-p(2,k+2), p(2,k+4)-p(2,k+3),p(2,k+5)-p(2,k+4)
c            write(*,'(a8,2i5,1p9e11.3)')'LV Gm/R4',
c     1           k,m,
c     2           grav*xm(k-2)*dmi(k-2)/(pi4*r(2,k-2)**4),
c     2           grav*xm(k-1)*dmi(k-1)/(pi4*r(2,k-1)**4),
c     2           grav*xm(k  )*dmi(k  )/(pi4*r(2,k  )**4),
c     2           grav*xm(k+1)*dmi(k+1)/(pi4*r(2,k+1)**4),
c     2           grav*xm(k+2)*dmi(k+2)/(pi4*r(2,k+2)**4),
c     2           grav*xm(k+3)*dmi(k+3)/(pi4*r(2,k+3)**4)

c           write(*,'(a8,2i5,1p9e11.3)')'LV ox',
c     1           k,m,
c     2           xold(nnuc,k-1),xold(nnuc,k),xold(nnuc,k+1),
c     3           xold(nnuc,k+2),xold(nnuc,k+3),xold(nnuc,k+4)
c            write(*,'(a8,2i5,1p9e11.3)')'LV x',
c     1           k,m,
c     2           x(nnuc,k-1),x(nnuc,k),x(nnuc,k+1),
c     3           x(nnuc,k+2),x(nnuc,k+3),x(nnuc,k+4)

c            write(*,'(a8,2i5,1p9e11.3)')'LVscrx',
c     1           k,m,
c     2           scrx(nnuc,k-1),scrx(nnuc,k),scrx(nnuc,k+1),
c     3           scrx(nnuc,k+2),scrx(nnuc,k+3),
c     4           scrx(nnuc,k+4),scrx(nnuc,k+5)

c            write(*,'(a8,2i5,1p9e11.3)')'LV Tm',
c     1           k,m,
c     2           t(m,k-1),t(m,k),t(m,k+1),t(m,k+2),
c     3           t(m,k+3),t(m,k+4)

c            write(*,'(a8,2i5,1p9e11.3)')'LV T2',
c     1           k,m,
c     2           t(2,k-1),t(2,k),t(2,k+1),t(2,k+2),
c     3           t(2,k+3),t(2,k+4)
c         endif
c      enddo
cccccccccccccccccccccccccccccccccccccccccc

      return
      end
