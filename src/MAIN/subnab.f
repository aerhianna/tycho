
      subroutine subnab(k,dndl,dndr,dndtp,dndvp,dndtm,dndvm,ytot)
c..   revised 2-22-09
c     uses fnab.f
c..   determines boehm-vitense nabla and derivatives
c..   called in hstat.f

      implicit none

c..   implements inefficient convection by Bohm-Vitense mixing length
c..   theory as in Kippenhahn and Weigert, p. 51, section 7.2
c..   added compositional difference (dnad --> dnad+doux)
c..   for modec = 1 --> 2

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'cenv'
      include 'cfnab'
      include 'cconst'

      integer*4 i,k,j

c     real*8    uu,psi,coef
      real*8    fact1,fact2,fact3
      real*8    tbar,vbar,pbar,akbar,ombeta
      real*8    ptbar,pvbar,etbar,evbar,beta

      real*8    tiny,fun,dfun,dpsi,w,betta,bettb,xx,vc
      real*8    tledd, dpsidu,dpsidw,denom
      real*8    dwdl,dwdr,dwdtp,dwdvp,dwdtm,dwdvm
      real*8    dudl,dudr,dudtp,dudvp,dudtm,dudvm
      real*8    dndl,dndr,dndtp,dndvp,dndtm,dndvm
      real*8    ytot(kdm),abpsi,alphasc,zetasn

c.. standard value from Langer, El Eid, Fricke 1985 A&A,145,179 is 0.1
      data alphasc/0.1d0/
c---------------------------------------------------------------------

c..   dnab's are defined in hstat.f before this call
c     if( modec .eq. 1 )then
c..   schwarzschild
c        w    = dnrad(k) - dnad(k)
c     elseif( modec .eq. 2 .or. modec .eq. 3 )then
c..   ledoux plus mixing
c     else
c        write(*,*)'modec ',modec
c        stop'subnab.f modec '
c     endif
c     coef = 8.0d0/9.0d0
      j = k+1
c..   boundary values of thermodynamic quantities
      tbar  = 0.5d0*( t(2,k)+t(2,j))
      pbar  = 0.5d0*( p(2,k)+p(2,j))
      vbar  = 0.5d0*( v(2,k)+v(2,j))
      akbar = 0.5d0*( ak(k) + ak(j))
      ombeta = 0.5d0*arad/3.0d0*(t(2,k)**4/p(2,k)+t(2,j)**4/p(2,j))
      
      ptbar = 0.5d0*( pt(k) + pt(j))
      pvbar = 0.5d0*( pv(k) + pv(j))
      etbar = 0.5d0*( et(k) + et(j))
      evbar = 0.5d0*( ev(k) + ev(j))

      tledd = pi4*crad*grav*xm(k)*ombeta/akbar

      call fnab(dnabv(k),dnad(k),dnrad(k),tbar,vbar,pbar,g(k),
     1 ptbar,pvbar,etbar,evbar,akbar,hp(k),doux(k),
     2 alphaml,uuml,modec,k)

      w    = dnrad(k) - dnad(k) - doux(k)

c..   Ledoux
c      dnabv(k) = psi**2 - uu**2 + dnad(k) + doux(k)

c..   dnabv determined, calculate velocity and fluxes
c..   formulation valid for both Schwarzschild and Ledoux

c..   coefficient of radiative heat flux
c..   (* dnab)
c      cfrad =  4.0d0*arad*crad/3.0d0*tbar**4/pbar/akbar*g(k)
      cfrad = 4.0d0*crad*g(k)/akbar*ombeta
      f(k)  = cfrad * dnabv(k)

      if( w .gt. 0.0d0 .and. uu .ne. 0.0d0)then
c..   convective
c..   boundary value of beta = - dlnP/dlnT / dlnP/dlnV
         betta = t(2,k)*pt(k)/(-v(2,k)*pv(k))
         bettb = t(2,j)*pt(j)/(-v(2,j)*pv(j))
c..   setup boundary value of Cp * beta
         fact1 = 
     1        (t(2,k)*et(k)/(p(2,k)*v(2,k)) 
     2        + betta*(1.0d0+ev(k)/p(2,k)))
     3        * betta
         fact2 =
     4        (t(2,j)*et(j)/(p(2,j)*v(2,j)) 
     5        + bettb*(1.0d0+ev(j)/p(2,j)))
     6        * bettb 
         fact3 = 0.5d0*( fact1 + fact2 )
c..   coefficients of convective velocity
c..   vcon**2 = cvcon*(nab - nabe) 
c..   = cvcon*(coef*uu*(rad-nab))**(2/3)
c..   8U/9 = f0/b0
         cfcon = pbar * sqrt( pbar*vbar*0.5d0 )
     1        *(alphaml/2.0d0)**2 * fact3
c..   coefficient of convective heat flux 
         xx     = coef*uu*(dnrad(k) - dnabv(k))
         b(k)   = cfcon * xx

c..   derivatives of psi
         denom  = 3.0d0*(psi-uu)**2 + 2.0d0*coef*uu*psi
         dpsidu = -(-3.0d0*(psi-uu)**2 
     1        +coef*(psi**2 -3.0d0*uu**2 -w) )/denom
         dpsidw = coef*uu /denom
c..   derivatives of w
         dwdl  = 0.25d0/tledd
         dwdr  = 0.0d0
c..   opacity and (1-beta) variation included
         dwdtp = dnrad(k) * 0.5d0 * (  akt(j)/akbar
     1        + (pt(j)/p(2,j) - 4.0d0/t(2,j))
     2        * arad*t(2,j)**4/p(2,j)/3.0d0/ombeta   )
         dwdvp = dnrad(k) * 0.5d0 * ( akv(j)/akbar
     1        +  pv(j)/p(2,j)
     2        * arad*t(2,j)**4/p(2,j)/3.0d0/ombeta   )
         dwdtm = dnrad(k)* 0.5d0 * ( akt(k)/akbar
     1        + (pt(k)/p(2,k) - 4.0d0/t(2,k))
     2        * arad*t(2,k)**4/p(2,k)/3.0d0/ombeta   )
         dwdvm = dnrad(k)* 0.5d0 * ( akv(k)/akbar
     1        + pv(k)/p(2,k)
     2        * arad*t(2,k)**4/p(2,k)/3.0d0/ombeta   )
c..   derivatives of uu
         dudl  = 0.0d0
         dudr  = -2.0d0*uu/r(2,k)
         dudtp = uu*0.5d0*( 4.0d0/tbar-2.5d0*pt(j)/pbar -akt(j)/akbar )
         dudvp = uu*0.5d0*(-0.5d0/vbar-2.5d0*pv(j)/pbar -akv(j)/akbar )
         dudtm = uu*0.5d0*( 4.0d0/tbar-2.5d0*pt(k)/pbar -akt(k)/akbar )
         dudvm = uu*0.5d0*(-0.5d0/vbar-2.5d0*pv(k)/pbar -akv(k)/akbar )
c..   finally, these are the derivatives of dnabv(k)
         dndl  = 2.0d0*psi*(dpsidu * dudl  + dpsidw * dwdl)
     1        - 2.0d0*uu*dudl
         dndr  = 2.0d0*psi*(dpsidu * dudr  + dpsidw * dwdr)
     1        - 2.0d0*uu*dudr
         dndtp = 2.0d0*psi*(dpsidu * dudtp + dpsidw * dwdtp)
     1        - 2.0d0*uu*dudtp
c         write(*,*)uu,k
         dndvp = 2.0d0*psi*(dpsidu * dudvp + dpsidw * dwdvp)
     1        - 2.0d0*uu*dudvp
         dndtm = 2.0d0*psi*(dpsidu * dudtm + dpsidw * dwdtm)
     1        - 2.0d0*uu*dudtm
         dndvm = 2.0d0*psi*(dpsidu * dudvm + dpsidw * dwdvm)
     1        - 2.0d0*uu*dudvm


c        write(*,'(3i5,1p12e12.3)')
c    1  it,k,ic(k),f(k),b(k),dndl,dndr,dndtp,dndvp,dndtm,dndvm,
c    2  w,w+doux(k),doux(k)
c        stop'subnab conv'
cccccccccccccc
      elseif( w .gt. -doux(k) )then
c..semiconvective
c..   radiative case so dnabv(k) = dnrad(k)
c..   w = dnrad(k) - dnad(k), so dw = d(dnrad(k)) so same
c..   as d(dnad(k)) taken as zero above
         beta = 0.5d0*rgas*(ytot(k)*t(2,k)/(p(2,k)*v(2,k)) 
     1        + ytot(j)*t(2,j)/(p(2,j)*v(2,j)) )
         zetasn = alphasc*0.5d0*( dnab(k) - dnad(k))
     1           /(dnad(k)+doux(k)-dnad(k))
         zetasn = dmax1(0.0d0,zetasn)
         abpsi =( 32.0d0 - 36.0d0*beta + 9.0d0*beta**2)
     1         /( 32.0d0 - 24.0d0*beta - 3.0d0*beta**2)


         dnabv(k) = ( dnrad(k) + zetasn*( dnad(k) + abpsi*doux(k)) )
     1           /( 1.0d0 + zetasn )

         hp(k) = zetasn * 0.5d0 * arad * 
     1   ( t(2,k)**3*v(2,k)/(et(k)-pt(k)*(ev(k)+p(2,k))/pv(k) )
     2    +t(2,j)**3*v(2,j)/(et(j)-pt(j)*(ev(j)+p(2,j))/pv(j) ) )
     3   *crad * 4.0d0/9.0d0 * a(j)/(akbar * dmi(j))


c         write(*,'(1p8e12.3)')beta,dnad(k),zetasn,abpsi,dnabv(k)-dnrad(k),
c     1  hp(k)

         uu       = 0
         psi      = 0
         cvcon    = 0
         cfcon    = 0

c..   derivatives of w
         dwdl  = 0.25d0/tledd
         dwdr  = 0.0d0
c..   opacity and (1-beta) variation included
         dwdtp = dnrad(k) * 0.5d0 * (  akt(j)/akbar
     1        + (pt(j)/p(2,j) - 4.0d0/t(2,j))
     2        * arad*t(2,j)**4/p(2,j)/3.0d0/ombeta   )
         dwdvp = dnrad(k) * 0.5d0 * ( akv(j)/akbar
     1        +  pv(j)/p(2,j)
     2        * arad*t(2,j)**4/p(2,j)/3.0d0/ombeta   )
         dwdtm = dnrad(k)* 0.5d0 * ( akt(k)/akbar
     1        + (pt(k)/p(2,k) - 4.0d0/t(2,k))
     2        * arad*t(2,k)**4/p(2,k)/3.0d0/ombeta   )
         dwdvm = dnrad(k)* 0.5d0 * ( akv(k)/akbar
     1        + pv(k)/p(2,k)
     2        * arad*t(2,k)**4/p(2,k)/3.0d0/ombeta   )

c         hp(k) = 0.0d0
         b(k)  = 0.0d0
         dndl  = dwdl  
         dndr  = dwdr   
         dndtp = dwdtp  
         dndvp = dwdvp  
         dndtm = dwdtm  
         dndvm = dwdvm  


c        write(*,'(3i5,1p12e12.3)')
c     1  it,k,ic(k),f(k),b(k),dndl,dndr,dndtp,dndvp,dndtm,dndvm,
c     2  w, w+doux(k),doux(k)
c         stop"semi"

      else
c..   radiative case so dnabv(k) = dnrad(k)
c..   w = dnrad(k) - dnad(k), so dw = d(dnrad(k)) so same
c..   as d(dnad(k)) taken as zero above
         dnabv(k) = dnrad(k)
         uu       = 0
         psi      = 0
         cvcon    = 0
         cfcon    = 0

c..   derivatives of w
         dwdl  = 0.25d0/tledd
         dwdr  = 0.0d0
c..   opacity and (1-beta) variation included
         dwdtp = dnrad(k) * 0.5d0 * (  akt(j)/akbar
     1        + (pt(j)/p(2,j) - 4.0d0/t(2,j))
     2        * arad*t(2,j)**4/p(2,j)/3.0d0/ombeta   )
         dwdvp = dnrad(k) * 0.5d0 * ( akv(j)/akbar
     1        +  pv(j)/p(2,j)
     2        * arad*t(2,j)**4/p(2,j)/3.0d0/ombeta   )
         dwdtm = dnrad(k)* 0.5d0 * ( akt(k)/akbar
     1        + (pt(k)/p(2,k) - 4.0d0/t(2,k))
     2        * arad*t(2,k)**4/p(2,k)/3.0d0/ombeta   )
         dwdvm = dnrad(k)* 0.5d0 * ( akv(k)/akbar
     1        + pv(k)/p(2,k)
     2        * arad*t(2,k)**4/p(2,k)/3.0d0/ombeta   )

         hp(k) = 0.0d0
         b(k)  = 0.0d0
         dndl  = dwdl  
         dndr  = dwdr   
         dndtp = dwdtp  
         dndvp = dwdvp  
         dndtm = dwdtm  
         dndvm = dwdvm  

c        write(*,'(3i5,1p12e12.3)')
c    1 it,k,ic(k),f(k),b(k),dndl,dndr,
c    1  dndtp,dndvp,dndtm,dndvm,w,w+doux(k),doux(k)

c        stop'rad'
      endif

      return

      end


