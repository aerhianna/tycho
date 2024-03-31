      subroutine getvec(xx,xmin,xmax,cvarx,cxlabel)
c..   8-2-06
c..   get vectors of variable from model and graph them

      implicit none

      include 'dimenfile'

      include 'comod'
      include 'czone'
      include 'cgen'
      include 'cconst'
      include 'cburn'

      include 'cenv'
      include 'compu'
      include 'csurface'
      include 'cnabla'
      include 'ceoset'

      real*4      xx(kdm),xmin,xmax
      real*8      rhobar,pvbar
      real*8      area, flux, aog, sum
      real*8      fact,fak,fak2,fak3,fak4
      real*8      metal,hydrogen,helium
      real*8      scre(kdm),scrp(kdm),depth(kdm)
      real*8      rich

      integer*4 k, i, n, ndum, kmax
      integer*4 j

      character*5 cvarx
      character*12 cxlabel

      logical rightflush
c----------------------------------------------------------------------

      if( modes .eq. 2 )then
         kmax = kk+1
      else
         kmax = kk
      endif

      if( cvarx .eq. 'lgt' )then
         cxlabel = 'log T(K)'
         do k = 2, kmax+1
            if( t(1,k) .gt. 0.0d0 )then
               xx(k-1) =  dlog10( t(1,k) )
            else
               xx(k-1) = 0.0d0
            endif
         enddo
      elseif( cvarx .eq. 'lgd' )then
         do k = 2, kmax+1
            cxlabel = 'log rho'
            if( v(1,k) .gt. 0.0d0 )then
               xx(k-1) = -dlog10( v(1,k) )
            else
               xx(k-1) = 0.0d0
            endif
         enddo
      elseif( cvarx .eq. 'rho' )then
         do k = 2, kmax+1
            cxlabel = 'rho(g/cc)'
            if( v(1,k) .gt. 0.0d0 )then
               xx(k-1) = 1.0d0/ v(1,k) 
            else
               xx(k-1) = 0.0d0
            endif
         enddo
      elseif( cvarx .eq. 'tem' )then
         cxlabel = 'T(K)'
         do k = 2, kmax+1
            xx(k-1) = t(1,k) 
         enddo
      elseif( cvarx .eq. 'pres' )then
         cxlabel = 'P(erg/cc)'
         do k = 2, kmax+1
            xx(k-1) = p(1,k) 
         enddo

      elseif( cvarx .eq. 'entro' )then
         cxlabel = 'Entropy'
         do k = 2, kmax+1
            xx(k-1) = entropy(k)
         enddo
      elseif( cvarx .eq. 'nrg' )then
         cxlabel = 'Energy/g'
         do k = 2, kmax+1
            xx(k-1) = e(2,k)
         enddo


      elseif( cvarx .eq. 'lgp' )then
         cxlabel = 'log P'
         do k = 2, kmax+1
            if( p(1,k) .gt. 0.0d0 )then
               xx(k-1) = dlog10( p(1,k) )
            else
               xx(k-1) = 0.0d0
            endif
         enddo
      elseif( cvarx .eq. 'lnp' )then
         cxlabel = 'ln P'
         do k = 2, kmax+1
            if( p(1,k) .gt. 0.0d0 )then
               xx(k-1) = dlog( p(1,k) )
            else
               xx(k-1) = 0.0d0
            endif
         enddo
         
      elseif( cvarx .eq. 'Hp' .or. cvarx .eq. 'hp' )then
c..   pressure scale height
         cxlabel = 'Hp'
         xx(1) = 0.
         do k = 3, kmax-1
            xx(k-1) = 0.5d0*(p(1,k) +p(1,k+1) )/(p(1,k) -p(1,k+1) )
     1           *0.5d0*( r(1,k+1) -r(1,k-1) )
         enddo
         do k = kmax, kmax+1
            xx(k-1) = xx(k-2)
         enddo

      elseif( cvarx .eq. 'pscal' )then
         xx(1) = 0.
         do k = 3, kmax-1
            cxlabel = 'P scale/R '
               xx(k-1) = 0.5d0*(p(1,k) +p(1,k+1) )/(p(1,k) -p(1,k+1) )
     1           *0.5d0*( r(1,k+1) -r(1,k-1) )
     2           + xx(k-2)
         enddo
         do k = kmax, kmax+1
            xx(k-1) = xx(k-2)
         enddo
         do k = 2, kmax
            xx(k-1) = xx(k-1)/r(1,kmax-1)
         enddo

      elseif( cvarx .eq. 'u' )then
         do k = 2, kmax+1
            cxlabel = 'u(cm/s)'
            xx(k-1) = u(1,k)
         enddo

      elseif( cvarx .eq. 'gradu' )then
         do k = 2, kmax+1
            cxlabel = 'grad u'
            xx(k-1) = (u(1,k)-u(1,k-1))/(r(1,k)-r(1,k-1))
         enddo



      elseif( cvarx .eq. 'sound' )then
         do k = 2, kmax+1
            cxlabel = 'sound(cm/s)'
            xx(k-1) = sound(k)
         enddo
      elseif( cvarx .eq. 'machc' )then
         cxlabel = 'mach cnv'
         do k = 2, kmax+1
            if( sound(k) .gt. 0.0d0 )then
               xx(k-1) = h(k)/sound(k)
            else
               xx(k-1) = 0.0
            endif
         enddo

      elseif(cvarx .eq. 'dddt' )then
c..   constant pressure: d ln rho/d ln T
         cxlabel = 'dln\\gr/dlnT'
         do k = 2, kmax+1
            xx(k-1) = -v(1,k)*pv(k)/(t(1,k)*pt(k))
         enddo

      elseif( cvarx .eq. 'cour' )then
         do k = 2, kmax+1
            cxlabel = 'courant(s)'
            if( sound(k) .gt. 0.0d0 .and. k .gt. 1 
     1           .and. k .lt. kmax  )then
c               xx(k-1) = ( r(1,k) - r(1,k-1) )/sound(k)
 
               xx(k-1) = 0.5d0*(p(1,k) +p(1,k+1) )/(p(1,k) -p(1,k+1) )
     1              *0.5d0*( r(1,k+1) -r(1,k-1) )/sound(k)

            else
               xx(k-1) = 0.0
            endif
         enddo


      elseif( cvarx .eq. 'U' )then
         cxlabel = 'U'
         do k = 2, kmax-1
            rhobar = ( 1.0d0/v(1,k) + 1.0d0/v(1,k+1) )*0.5d0
            xx(k-1) = pi4*r(1,k)**3*rhobar/xm(k)
         enddo
c..   edge values: join and photosphere
         if( v(1,kmax+1) .gt. 0.0d0 )then
            xx(kmax) = pi4*r(1,kmax)**3/xm(kmax)/v(1,kmax+1)
         else
            xx(kmax) = 0.0
         endif
         if( v(1,kmax-1) .gt. 0.0d0 )then
            xx(kmax-1) = pi4*r(1,kmax-1)**3/xm(kmax-1)/v(1,kmax-1)
         else
            xx(kmax-1) = 0.0
         endif
ccccccccccccccccccccccccccccccccc
         write(*,'(i5,1p8e12.3)')kmax,xx(kmax),v(1,kmax),
     1        r(1,kmax),xm(kmax),v(1,kmax+1)
         write(*,'(i5,1p8e12.3)')kmax-1,xx(kmax-1),v(1,kmax-1),
     1        r(1,kmax-1),xm(kmax-1)
cccccccccccccccccccccccccccccccccccccc

      elseif( cvarx .eq. 'V' )then
         cxlabel = 'V'
         do k = 2, kmax
            pvbar = ( p(1,k)*v(1,k) + p(1,k+1)*v(1,k+1) )*0.5d0
            xx(k-1) = grav*xm(k)/( r(1,k) * pvbar )
         enddo
c..   edge values: join and photosphere
         pvbar = p(1,kmax)*v(1,kmax)
         if( pvbar .gt. 0.0d0 )then
            xx(kmax) = grav*xm(kmax)/r(1,kmax)/pvbar
         else
            xx(kmax) = 0.0d0
         endif

      elseif( cvarx .eq. 'h' )then
         do k = 2, kmax+1
            cxlabel = 'vconv(cm/s)'
            xx(k-1) = h(k)
         enddo
      elseif( cvarx .eq. 'xm' )then
         cxlabel = 'M(sol)'
         do k = 2, kmax
            xx(k-1) = 0.5d0*(xm(k-1) + xm(k))/sol
         enddo
         xx(kmax) = xm(kmax)/sol
      elseif( cvarx .eq. 'lgr' )then
         cxlabel = 'log R(cm)'
         do k = 2, kmax
            if( r(1,k) .gt. 0.0d0 )then
               xx(k-1) = dlog10( 0.5d0*(r(1,k-1) + r(1,k)) )
            else
               xx(k-1) = 0.0d0
            endif
            xx(kmax) = dlog10( r(1,kmax) )
         enddo

      elseif( cvarx .eq. 'dr' )then
         cxlabel = 'dR(cm)'
         do k = 2, kmax
            xx(k-1) = 0.5d0*(r(1,k+1) - r(1,k-1))
         enddo
         xx(kmax) = r(1,kmax)-r(1,kmax-1)
      elseif( cvarx .eq. 'lgdr' )then
         cxlabel = 'log dR(cm)'
         do k = 2, kmax
            xx(k-1) = dlog10( 0.5d0*(r(1,k+1) - r(1,k-1)) )
         enddo
         xx(kmax) = dlog10( r(1,kmax)-r(1,kmax-1) )

      elseif( cvarx .eq. 'r' .or. cvarx .eq. 'R' )then
c..centered at boundaries
         cxlabel = 'R(cm)'
         do k = 2, kmax
            xx(k-1) = 0.5d0*(r(1,k-1) + r(1,k))
         enddo
         xx(kmax) = r(1,kmax)

      elseif( cvarx .eq. 'rh' )then
c..centered on zones
         cxlabel = 'Rh(cm)'
         xx(1) = 0.5d0*(r(1,2)-r(1,1))
         do k = 3, kmax
            xx(k-1) = dsqrt(r(1,k-1) * r(1,k))
         enddo
         xx(kmax) = r(1,kmax)
         
      elseif( cvarx .eq. 'k'  )then
         cxlabel = 'zones'
         do k = 1, kmax 
            xx(k) = float(k)
         enddo
      elseif( cvarx .eq. 'khalf'  )then
         cxlabel = 'zones'
         do k = 1, kmax 
            xx(k) = float(k) - 0.5
         enddo

      elseif( cvarx .eq. 'xl' )then
         cxlabel = 'L/sol'
         do k = 2, kmax +1
            xx(k-1) = tl(1,k) / sollum
         enddo
      elseif( cvarx .eq. 'lgkap' )then
c..   log10 of opacity ak(k)
         cxlabel = 'lgkap'
         do k = 2, kmax +1
            if( ak(k) .gt. 0.0d0 )then
               xx(k-1) = dlog10( ak(k) )
            else
               xx(k-1) = -3.0
            endif
         enddo
      elseif( cvarx .eq. 'kappa' )then
c..   radiative opacity ak(k)
         cxlabel = 'kappa'
         do k = 2, kmax +1
            xx(k-1) = ak(k) 
         enddo
      elseif( cvarx .eq. 'aloak' )then
c..   effective opacity from structure and luminosity
c..   "luminosity opacity" / state.f opacity
c..   unity if model is consistent with state.f
         cxlabel = 'aloak'
         do k = 2, kmax -1
            a(k) = pi4*r(1,k)**2
            xx(k-1) = arad*crad/3.0d0 * (a(k)**2/tl(1,k))
     1           *( t(1,k)**4 - t(1,k+1)**4)/dmi(k)
     2           *2.0d0/( ak(k) + ak(k+1) )
         enddo
         xx(kmax) = 0.
         xx(kmax+1) = 0.

      elseif( cvarx .eq. 'arad' )then
         cxlabel = 'arad/g'
         do k = 2, kmax
            area = pi4*r(1,k)**2
c..kappa omitted in flux definition because it cancels in aog [a/g]
            if( k .lt. kmax-1 )then
            flux = arad*crad/3.0d0 * (t(1,k)**4-t(1,k+1)**4)
     1           * area/dmi(k)
            elseif( k .eq. kmax-1 )then
c..outer zone on grid, at join
            flux = arad*crad/3.0d0 * (t(1,k)**4-t(1,k+1)**4)
     1           * area/dmh(k)
            else
c..envelope: photospheric condition
               flux = tl(1,k)/area /1.5d0
            endif 
            aog = area*flux/xm(k)/(pi4*grav*crad)
            xx(k-1) = aog
         enddo
         xx(kmax) = xx(kmax-1)
         xx(kmax+1) = xx(kmax-1)

      elseif( cvarx .eq. 'lgrop' )then
         cxlabel = 'lg Ropal'
         do k = 2, kmax +1
            if( t(1,k) .gt. 0.0d0 .and. v(1,k) .gt. 0.0d0 )then
               xx(k-1) = -3.0d0*dlog10( t(1,k)*1.0d-6) 
     1              - dlog10(v(1,k) )
            else
               xx(k-1) = -20.0
            endif
         enddo
      elseif( cvarx .eq. 'lgg' )then
         cxlabel = 'log g'
         do k = 2, kmax
            if( r(1,k) .gt. 0.0d0 .and.xm(k) .gt. 0.0d0 )then
               xx(k-1) = dlog10( grav*xm(k)/r(1,k)**2 )
            else
               xx(k-1) = -20.0
            endif
         enddo
      elseif( cvarx .eq. 'g' )then
         cxlabel = 'g'
         do k = 2, kmax
            if( r(1,k) .gt. 0.0d0 .and.xm(k) .gt. 0.0d0 )then
               xx(k-1) = ( grav*xm(k)/r(1,k)**2 )
            else
               xx(k-1) = 0.0
            endif
         enddo
      elseif( cvarx .eq. 'dmh' )then
         cxlabel = 'zone mass'
         do k = 2, kmax
            xx(k-1) = dmh(k)/sol
         enddo

      elseif( cvarx .eq. 'dmi' )then
         cxlabel = 'bnd. mass'
         do k = 2, kmax
            xx(k-1) = dmi(k)/sol
         enddo

      elseif( cvarx .eq. 'lgdmh' )then
         cxlabel = 'lg zone mass'
         do k = 2, kmax
            if( dmh(k) .gt. 0.0d0 )then
               xx(k-1) = dlog10( dmh(k)/sol )
            else
               xx(k-1) = -20.0d0
            endif
         enddo

      elseif( cvarx .eq. 'dmdk' )then
         cxlabel = 'dm/dk'
         do k = 2, kmax
            if( k .lt. 1 .or. k .gt. kmax-1 )then
               xx(k-1) = 0.0d0
            else
               xx(k-1) = (
     1              ( dmh(k+1)-dmh(k) ) /
     2              (dmh(k+1)+dmh(k) )*2.0d0
     3              )

            endif
         enddo


      elseif( cvarx .eq. 'd2mdk' )then
         cxlabel = 'd2m/dk2'
         do k = 2, kmax
            if( k .lt. 1 .or. k .gt. kmax )then
               xx(k-1) = 0.0d0
            elseif( k .eq. 2 )then
               xx(k-1) = ( dmh(k+1)-dmh(k) )/(dmh(k+1)+dmh(k))*2.0d0
            elseif( k .eq. kmax )then
               xx(k-1) = ( dmh(k)-dmh(k-1) )/(dmh(k)+dmh(k-1))*2.0d0
            else
               xx(k-1) = (
     1              (dmh(k+1)-2.0d0*dmh(k)+dmh(k-1))/
     2              (dmh(k+1)+dmh(k)+dmh(k-1))*3.0d0
     3              )
     4              *(dmh(k+1) - dmh(k-1))
     5              /(dmh(k+1)+dmh(k)+dmh(k-1))*3.0d0
            endif
         enddo


      elseif( cvarx .eq. 'd3mdk' )then
         cxlabel = 'd3m/dk3'
         do k = 2, kmax
            if( k .le. 2 .or. k .ge. kmax )then
               xx(k-1) = 0.0d0
            elseif( k .eq. 3 )then
               xx(k-1) = 0.0d0
            elseif( k .eq. kmax-1 )then
               xx(k-1) = 0.0d0
            else
               xx(k-1) = abs(
     1              (dmh(k+1)-3.0d0*dmh(k)+3.0d0*dmh(k-1)-dmh(k-2))/
     2              (dmh(k+1)+dmh(k)+dmh(k-1)+dmh(k-2))*4.0d0
     3              )
            endif
         enddo

      elseif( cvarx .eq. 'dmerr' )then
         cxlabel = 'dm err'
         do k = 2, kmax
            if( k .le. 2 .or. k .ge. kmax )then
               xx(k-1) = 0.0d0
            else
               xx(k-1) = ( dmh(k+1)**2 - dmh(k)**2 ) * 0.25d0
     1              /dmi(k)**2
            endif
         enddo



      elseif( cvarx .eq. 'xlrad' )then
         cxlabel = 'Lrad/sol'
         do k = 2, kmax-1
            xx(k-1) = (pi4*r(1,k)**2)**2/dmi(k)*cflux*
     1           (t(1,k)**4-t(1,k+1)**4)/( 0.5d0*( ak(k) + ak(k+1)) )
     2           /sollum
         enddo

      elseif( cvarx .eq. 'xlcnv' )then
c..   implied convective luminosity = L - L(radiative)
         cxlabel = 'Lcnv/sol'
         xx(kmax) = 0
         xx(kmax-1) = 0
         do k = 2, kmax-1
            xx(k-1) = ( tl(1,k) 
     1           -(pi4*r(1,k)**2)**2/dmi(k)*cflux
     2           *(t(1,k)**4-t(1,k+1)**4)/( 0.5d0*( ak(k) + ak(k+1)) ) )
     3           /sollum
         enddo         
         
         elseif( cvarx .eq. 'lcke' )then
c..   equivalent luminosity of kinetic energy
            cxlabel = 'L(cke)'
            do k = 2, kmax-1
               xx(k-1) = pi4*r(1,k)**2*h(k)**3*2.0d0/(v(1,k)+v(1,k+1))
     1              /sollum
            enddo
ccccccccccccccccccccccccccccccccccccccccccccccccc
c            sum = 0.0d0
c            do k = 2, kmax-1
c               sum = sum + dmi(k)*h(k)**2 *0.5d0
cc               write(*,'(i5,1p8e12.3)')k,r(2,k),h(k),dmi(k),sum
c            enddo
c            write(*,'(a20,1pe12.3)')'convective KE',sum
c            sum = 0.0d0
c            do k = 2, kmax-1
c               if( h(k) .gt. 1.0d0 )then
c                  sum = sum + dmh(k)*e(1,k)
c               endif
cc               write(*,'(i5,1p10e12.3)')k,dmh(k),r(1,k),
cc     1              r(2,k),h(k),e(1,k),sum
ccccccccc
c            enddo
c            write(*,'(a20,1pe12.3)')'envelope E',sum
cccccccccccc

      elseif( cvarx .eq. 'du' )then
         cxlabel = 'du'
         do k = 2, kmax-1
            xx(k-1) = -grav*xm(k)/r(1,k)**2 
     1           -pi4*r(1,k)**2 * (p(1,k+1)-p(1,k))/dmi(k) 
         enddo 

       elseif( cvarx .eq. 'Adpdm' )then
         cxlabel = '-AdP/dm'
         do k = 2, kmax-1
            xx(k-1) = 
     1           -pi4*r(1,k)**2 * (p(1,k+1)-p(1,k))/dmi(k) 
         enddo 
         

      elseif( cvarx .eq. 'dnab' )then
c..actual gradient
         cxlabel = 'dnab'
         do k = 2, kmax-1
            xx(k-1) = dnab(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0
         
      elseif( cvarx .eq. 'dnad' )then
c..adiabatic gradient
         cxlabel = 'dnad'
         do k = 2, kmax-1
            xx(k-1) = dnad(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0
         
      elseif( cvarx .eq. 'dnrad' )then
c..radiative gradient
         cxlabel = 'dnrad'
         do k = 2, kmax-1
            xx(k-1) = dnrad(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0

      elseif( cvarx .eq. 'doux' )then
c..ledoux gradient
         cxlabel = 'doux'
         do k = 2, kmax-1
            xx(k-1) = doux(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0

      elseif( cvarx .eq. 'r-ad' )then
c..radiative gradient - adiabatic gradient
         cxlabel = 'r-ad'
         do k = 2, kmax-1
            xx(k-1) = dnrad(k)- dnad(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0

      elseif( cvarx .eq. 'r-a-x' )then
c..radiative gradient - adiabatic gradient
         cxlabel = 'r-a-x'
         do k = 2, kmax-1
            xx(k-1) = dnrad(k)- dnad(k) -doux(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0

      elseif( cvarx .eq. 'n-a' )then
c..radiative gradient - adiabatic gradient
         cxlabel = 'n-a'
         do k = 2, kmax-1
            xx(k-1) = dnab(k)- dnad(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0

      elseif( cvarx .eq. 'n-a-x' )then
c..radiative gradient - adiabatic gradient - ledoux
         cxlabel = 'n-a-x'
         do k = 2, kmax-1
            xx(k-1) = dnab(k)- dnad(k) -doux(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0


      elseif( cvarx .eq. 'gwfak' )then
c..gravity wave opacity enchancement
         cxlabel = 'a/n-1'
         do k = 2, kmax-1
            if( dnab(k) .lt. dnad(k) )then
               xx(k-1) = dnad(k)/dnab(k) - 1.0d0
            else
               xx(k-1) = 0.0d0
            endif
         enddo

      elseif( cvarx .eq. 'bvai2' )then
c..brunt vaisala frequency squared
         cxlabel = 'N(BV)**2'
         do k = 2, kmax-1
c            xx(k-1) = ( -dnab(k)+ dnad(k) +doux(k) )*
c     1           (grav*xm(k)/r(1,k)**2)**2/
c     1           (0.5d0*( p(1,k)*v(1,k) + p(1,k+1)*v(1,k+1) ))
            xx(k-1) = nsqr(k)
         enddo 
         xx(kmax) = 0
         xx(kmax-1) = 0
         
      elseif( cvarx .eq. 'Ri' )then

         cxlabel = 'Ri'

c         write(*,'(a5,11a11)')'k','r','xm','shear','nsqr','ventr',
c     1        'Froude','H','h','xx'
         do k = 2, kmax-1
c..velocity shear**2
            fact = 0.5d0*( ( h(k+1)-h(k))/(r(1,k+1)-r(1,k)) )**2
     1           + 0.5d0*( ( h(k)-h(k-1))/(r(1,k)-r(1,k-1)) )**2
c..buoyancy frequency**2
            fak = nsqr(k)
            if( fact .gt. 0.0d0 )then
c..Richardson number= 1 / Froude number
               rich = fak/fact
            else
c..   no mixing for Ri > 0.25
               rich = 2.0d0
            endif
            xx(k-1) = rich
c..entrainment velocity
c            if( rich .gt. 0.25d0 )then
c..   nonconvective
c               fak2 = (r(1,k+1)-r(1,k-1))/dth(2) * 0.05d0
c..keep sane for tiny dth(2)
c               fak2 = dmin1( fak2, 1.0d-2 * sound(k) )
c            else
c..   convective
c               fak2 = h(k)
c            endif

c            if( k .ge. 630 .and. k .le. 640 )then
c               write(*,'(i5,1p11e11.3)')k,r(1,k),xm(k),fact,nsqr(k),
c     1              fak2,fact/nsqr(k),x(nnuc-1,k),h(k),xx(k-1)
c            endif
         enddo 
c         write(*,'(a5,11a11)')'k','r','xm','shear','nsqr','ventr',
c     1        'Froude','H','h','xx'

c         write(*,'(a5,11a11)')'k',
c     1        'nab','nad', 'doux','dnabn','dnrad'

c         do k = 630, 640
c            write(*,'(i5,1p12e11.3)')k,
c     1           dnab(k),dnad(k),doux(k),dnabv(k),dnrad(k)
c         enddo
c         write(*,'(a5,12a11)')'k',
c     1        'nab','nad', 'doux','dnabv','dnrad'
c         write(*,*)kmax,kk
c         stop'getvec'
ccccccccccccccccccccccc







      elseif( cvarx .eq. 'vc/l2' )then
c..Richardson number = N**2/(du/dx)**2
         cxlabel = 'vc/l2'
         do k = 2, kmax-1
            xx(k-1) = ( h(k) * grav*xm(k)/r(1,k)**2 /
     1           (0.5d0*( p(1,k)*v(1,k) + p(1,k+1)*v(1,k+1) )) )**2
         enddo 

      elseif( cvarx .eq. 'l/r' )then
c..pressure scale height / radius
         cxlabel = 'ml/r'
         do k = 2, kmax-1
            xx(k-1) = 0.5d0*( p(1,k)*v(1,k) + p(1,k+1)*v(1,k+1) )
     1           * r(1,k)/( grav *xm(k) )
         enddo 

      elseif( cvarx .eq. 'ss' )then
c..nuclear energy generation rate (erg/g/sec)
         cxlabel = 'nuclear'
         do k = 2, kmax-1
            xx(k-1) = ss(k)
         enddo 
       elseif( cvarx .eq. 'sa' )then
c.. d ln eps /d ln T
         cxlabel = 'dlneps/dlnT'
         do k = 2, kmax-1
            xx(k-1) = sa(k)
         enddo 
            elseif( cvarx .eq. 'sb' )then
c.. d ln eps /d ln V
         cxlabel = 'dlneps/dlnV'
         do k = 2, kmax-1
            xx(k-1) = sb(k)
         enddo 
    
      elseif( cvarx .eq. 'ss' )then
c..nuclear energy generation rate (erg/g/sec)
         cxlabel = 'nuclear'
         do k = 2, kmax-1
            xx(k-1) = ss(k)
         enddo 
       elseif( cvarx .eq. 'sa' )then
c.. d ln eps /d ln T
         cxlabel = 'dlneps/dlnT'
         do k = 2, kmax-1
            xx(k-1) = sa(k)
         enddo 
            elseif( cvarx .eq. 'sb' )then
c.. d ln eps /d ln V
         cxlabel = 'dlneps/dlnV'
         do k = 2, kmax-1
            xx(k-1) = sb(k)
         enddo 
    
      elseif( cvarx .eq. 'ss' )then
c..nuclear energy generation rate (erg/g/sec)
         cxlabel = 'nuclear'
         do k = 2, kmax-1
            xx(k-1) = ss(k)
         enddo 
       elseif( cvarx .eq. 'sa' )then
c.. d ln eps /d ln T
         cxlabel = 'dlneps/dlnT'
         do k = 2, kmax-1
            xx(k-1) = sa(k)
         enddo 
            elseif( cvarx .eq. 'sb' )then
c.. d ln eps /d ln V
         cxlabel = 'dlneps/dlnV'
         do k = 2, kmax-1
            xx(k-1) = sb(k)
         enddo 
    
      elseif( cvarx .eq. 'snu' )then
c..nuclear energy generation rate (erg/g/sec)
         cxlabel = 'neutrino'
         do k = 2, kmax-1
            xx(k-1) = -snu(k)
         enddo 
       elseif( cvarx .eq. 'snua' )then
c.. d ln eps /d ln T
         cxlabel = 'dlnsnu/dlnT'
         do k = 2, kmax-1
            xx(k-1) = snua(k)
         enddo 
            elseif( cvarx .eq. 'snub' )then
c.. d ln eps /d ln V
         cxlabel = 'dlnsnu/dlnV'
         do k = 2, kmax-1
            xx(k-1) = snub(k)
         enddo 
           
      elseif( cvarx .eq. 'ss-nu' )then
c..nuclear energy generation rate (erg/g/sec)
         cxlabel = 'nuc-neut'
         do k = 2, kmax-1
            xx(k-1) = ss(k)+snu(k)
         enddo 

      elseif( cvarx .eq. 'ln-lu' )then
c..nuclear energy generation rate (erg/g/sec)
         cxlabel = 'Lnuc-Lneut'
         xx(1) = 0.0
         do k = 2, kmax-1
            xx(k) = xx(k-1) + (ss(k)+snu(k))*dmh(k)/sol
         enddo 
      elseif( cvarx .eq. 'ln' )then
c..nuclear energy generation rate (erg/g/sec)
         cxlabel = 'Lnuc'
         xx(1) = 0.0
         do k = 2, kmax-1
            xx(k) = xx(k-1) + ss(k)*dmh(k)/sol
         enddo 
      elseif( cvarx .eq. 'lu' )then
c..nuclear energy generation rate (erg/g/sec)
         cxlabel = '-Lneut'
         xx(1) = 0.0
         do k = 2, kmax-1
            xx(k) = xx(k-1) -snu(k)*dmh(k)/sol
         enddo 

      elseif( cvarx .eq. 'd2adm' )then
c..
         cxlabel = 'd2ad4'
         xx(1) = 0.0
c..   second derivative of He4 abundance, with respect to 
c..   mass coordinate in solar units
         xx(2) =  abs( x(nnuc,3)-x(nnuc,2) ) / dmi(2)*sol
     1        /dmh(2)*sol *4.0d0
         do k = 3, kmax-1
            xx(k) = abs( 
     1       ( x(nnuc,k+1) - x(nnuc,k)   ) /dmi(k)*sol
     2     - ( x(nnuc,k)   - x(nnuc,k-1) ) /dmi(k-1)*sol
     3           )/dmh(k)*sol *4.0d0
         enddo 
         xx(kmax) = 0.0
ccccccccccccccccccccccccccccccccccccccccc

      elseif( cvarx .eq. 'd2he4' )then
c..
         cxlabel = 'd2he4'
         xx(1) = 0.0
c..   second derivative of He4 abundance, with respect to 
c..   zone index
         xx(2) =  abs( x(nnuc,3)-x(nnuc,2) )*4.0d0
         do k = 3, kmax-1
            xx(k) = abs( 
     1       ( x(nnuc,k+1) - 2.0d0*x(nnuc,k)  +x(nnuc,k-1) )
     2           )*4.0d0
         enddo 
         xx(kmax) = 0.0
ccccccccccccccccccccccccccccccccccccccccc

      elseif( cvarx .eq. 'm/r4' )then
c..
         cxlabel = 'HSE test'
         do k = 2, kmax-1
            if( r(1,k) .gt. 0.0d0 )then
               xx(k) = -grav/pi4*xm(k)*dmi(k)/r(1,k)**4
     1              /(p(1,k+1) - p(1,k) ) -1.0d0
            else
               xx(k) = 0.0
            endif
         enddo 
         xx(1)    = 0.0
         xx(kmax) = 0.0
ccccccccccccccccccccccccccccccccccccccccc

      elseif( cvarx .eq. 'plin' )then
c..
         cxlabel = 'P linear'
         g(1) = 0.0d0
         a(1) = 0.0d0
         do k=2, kmax
            g(k) = grav*xm(k)/r(1,k)**2
            a(k) = pi4*r(1,k)**2
         enddo
         p(1,1) = p(1,2)
         do k = 2, kmax-1
            fact = ( 
     1           ( p(1,k+1)-p(1,k)   )/dmi(k)
     2           -(  p(1,k)-p(1,k-1) )/dmi(k-1)  )/dmh(k)
            fak = fact*2.0d0/( a(k)*g(k)-a(k-1)*g(k-1) )
     1           *( dmi(k)**2 - dmi(k-1)**2 )
            xx(k-1) = fak
            write(*,'(i5,1p12e12.3)')k,fact,fak,dmi(k),dmh(k),
     1           g(k),xx(k-1)

c     2           ( p(1,k+1)-p(1,k)   )/dmi(k),
c     3           ( p(1,k)  -p(1,k-1) )/dmi(k-1) 

c            if( k .gt. 5 )stop'getvec'
         enddo
         stop'getvec plin'
ccccccccccccccccccccccccccccccccccccccccc

      elseif( cvarx .eq. 'r3' )then
c..
         cxlabel = 'V lin.'
         xx(1) = 0.0
         do k = 2, kmax-2
            xx(k) = pi43*(  
     1             r(1,k+1)**3/(dmh(k+1)*v(1,k+1))
     2           - r(1,k  )**3/(dmh(k+1)*v(1,k+1)) 
     3           - r(1,k  )**3/(dmh(k  )*v(1,k  )) 
     4           + r(1,k-1)**3/(dmh(k  )*v(1,k  )) 
     5           )
         enddo 
         write(*,*)kmax,xx(1),xx(kmax-1),xx(kmax)
         xx(kmax-1) = 0.0
         xx(kmax) = 0.0

      elseif( cvarx .eq. 'kec' )then
c..convective kinetic energy (mixing length)
         cxlabel = 'KE(conv)'
         do k = 1, kmax-2
            xx(k) = dmi(k)*0.5d0*h(k)**2
         enddo
         xx(kmax) = 0.0d0

      elseif( cvarx .eq. 'lgh1' )then
c..
         cxlabel = 'log H1'
c..   log10 of proton abundance
         do k = 2, kmax
            if( x(nnuc-1,k) .gt. 0.0d0 )then
               xx(k) = dlog10( x(nnuc-1,k) )
            else
               xx(k) = -20.0
            endif
         enddo 
         
c..reflection about origin assumed
         xx(1) = xx(2)

      elseif( cvarx .eq. 'metal' )then
c..
         cxlabel = 'metals'
c..   "metallicity, helium, hydrogen"
         do k = 2, kmax
            xx(k) = 0
            metal = 0.0d0
            helium = 0.0d0
            hydrogen = 0.0d0
            do n = 1, nnuc
               if( lz(n) .gt. 2 )then
                  metal = metal +  x(n,k)
               elseif( lz(n) .eq. 2 )then
                  helium = helium + x(n,k)
               elseif( lz(n) .eq. 1 )then
                  hydrogen = hydrogen + x(n,k)
               endif
            enddo
            xx(k) = metal
c            write(*,'(i5,1p8e12.3)')k,xx(k),metal,helium,hydrogen,
c     1           metal+helium+hydrogen-1.0d0
         enddo 

c..reflection about origin assumed
         xx(1) = xx(2)
      elseif( cvarx .eq. 'heli' )then
c..
         cxlabel = 'Helium'
c..   "metallicity, helium, hydrogen"
         do k = 2, kmax
            xx(k) = 0
            metal = 0.0d0
            helium = 0.0d0
            hydrogen = 0.0d0
            do n = 1, nnuc
               if( lz(n) .gt. 2 )then
                  metal = metal +  x(n,k)
               elseif( lz(n) .eq. 2 )then
                  helium = helium + x(n,k)
               elseif( lz(n) .eq. 1 )then
                  hydrogen = hydrogen + x(n,k)
               endif
            enddo
            xx(k) = helium
c            write(*,'(i5,1p8e12.3)')k,xx(k),metal,helium,hydrogen,
c     1           metal+helium+hydrogen-1.0d0
         enddo 

c..reflection about origin assumed
         xx(1) = xx(2)
      elseif( cvarx .eq. 'hydr' )then
c..
         cxlabel = 'Hydrogen'
c..   "metallicity, helium, hydrogen"
         do k = 2, kmax
            xx(k) = 0
            metal = 0.0d0
            helium = 0.0d0
            hydrogen = 0.0d0
            do n = 1, nnuc
               if( lz(n) .gt. 2 )then
                  metal = metal +  x(n,k)
               elseif( lz(n) .eq. 2 )then
                  helium = helium + x(n,k)
               elseif( lz(n) .eq. 1 )then
                  hydrogen = hydrogen + x(n,k)
               endif
            enddo
            xx(k) = hydrogen
c            write(*,'(i5,1p8e12.3)')k,xx(k),metal,helium,hydrogen,
c     1           metal+helium+hydrogen-1.0d0
         enddo 

c..reflection about origin assumed
         xx(1) = xx(2)

      elseif( cvarx .eq. 'cn' )then
c..   CN cycle nuclei, mole fraction
         cxlabel = 'Y(CN)'

         do k = 2, kmax+1
            sum = 0
            do i = 1, nnuc
               if( lz(i) .eq. 6 .or. lz(i).eq. 7 )then
                  sum = sum + x(i,k)/xa(i)
               endif
            enddo
            xx(k-1) = sum
         enddo

      elseif( cvarx .eq. 'cno' )then
c..   CNO cycle nuclei, mole fraction
         cxlabel = 'Y(CNO)'

         do k = 2, kmax+1
            sum = 0
            do i = 1, nnuc
               if( lz(i) .eq. 6 .or. lz(i).eq. 7 
     1              .or. lz(i) .eq. 8 )then
                  sum = sum + x(i,k)/xa(i)
               endif
            enddo
            xx(k-1) = sum
         enddo

      elseif( cvarx .eq. 'cno6' )then
c..   CNO cycle nuclei, mole fraction, only O16
c..   (testing JNBahcall)
         cxlabel = 'Y(CNO6)'

         do k = 2, kmax+1
            sum = 0
            do i = 1, nnuc
               if( lz(i) .eq. 6 .or. lz(i).eq. 7 )then
                  sum = sum + x(i,k)/xa(i)
               endif
c..only O16
               if( lz(i) .eq. 8 .and. ln(i).eq. 8 )then
                  sum = sum + x(i,k)/xa(i)
               endif
            enddo
            xx(k-1) = sum
         enddo

      elseif( cvarx .eq. 'cn7' )then
c..   CNO cycle nuclei, mole fraction, only O17
c..   (testing JNBahcall)
         cxlabel = 'Y(CN7)'

         do k = 2, kmax+1
            sum = 0
            do i = 1, nnuc
               if( lz(i) .eq. 6 .or. lz(i).eq. 7 )then
                  sum = sum + x(i,k)/xa(i)
               endif
c..only O17
               if( lz(i) .eq. 8 .and. ln(i).eq. 9 )then
                  sum = sum + x(i,k)/xa(i)
               endif
            enddo
            xx(k-1) = sum
         enddo

      elseif( cvarx .eq. 'ysum' )then
c..   sum of free particle mode fractions
         cxlabel = 'Ye+Yion'
         do k = 2, kmax+1
            sum = 0
            do i = 1, nnuc
               sum = sum + x(i,k)/xa(i)*(1.0d0 + dble(lz(i)))
            enddo
            xx(k-1) = sum
         enddo

      elseif( cvarx .eq. 'dlydm' )then
c..d ln Y /dm
         cxlabel = 'dlnY/dm'
c..calculate Ytotal=Ye+Yion
         do k = 2, kmax+1
            sum = 0
            do i = 1, nnuc
               sum = sum + x(i,k)/xa(i)*(1.0d0 + dble(lz(i)))
            enddo
            xx(k-1) = sum
         enddo
c..calculate the gradient (dimensionless units)
         do k = 2, kmax-2
            xx(k-1) = (xx(k+1)-xx(k))/(xx(k+1)+xx(k))*2.0d0
     1           *(xm(kmax)/dmi(k))
         enddo
         xx(kmax) = 0.
         xx(kmax-1) = 0.

      elseif( cvarx .eq. 'yef' )then
c..   sum of free particle mode fractions
         cxlabel = 'Ye(ionized)'
         do k = 2, kmax
            xx(k-1) = yef(k)
         enddo
         xx(kmax) = 0

      elseif( cvarx .eq. 'gam1' )then
c..   adiabatic exponent
         cxlabel = 'gamma1'
         do k = 2, kmax+1
            xx(k-1) = gamma1(k)
         enddo

      elseif( cvarx .eq. 'mfp' )then
c..   photon mean free path
         cxlabel = 'mfp'
         do k = 2, kmax
            xx(k-1) = v(1,k)/ak(k)
         enddo
         
         

      elseif( cvarx .eq. 'opdep' )then
c..optical depth inward from photosphere
         cxlabel = 'op. depth'
         depth(kmax) = 0.7d0
         do j = 2, kmax
            k = kmax - j +2
            depth(k-1) = depth(k) + (r(1,k)-r(1,k-1))/v(1,k)*ak(k)
         enddo
         do k = 1,kmax
            xx(k) = depth(k)
         enddo
         
      elseif( cvarx .eq. 'lgodp' )then
c..log of optical depth inward from photosphere
         cxlabel = 'lg op. depth'
         depth(kmax) = 0.7d0
         do j = 2, kmax
            k = kmax - j +2
            depth(k-1) = depth(k) + (r(1,k)-r(1,k-1))/v(1,k)*ak(k)
         enddo
         do k = 1,kmax
            if( depth(k) .gt. 0.0d0)then
               xx(k) = dlog10( depth(k) )
            else  
               write(*,*)k,kmax,depth(k)
               stop'getvec: sign error in optical depth'
            endif
         enddo

      elseif( cvarx .eq. 'depth' )then
c..radial distance inward from photosphere
         cxlabel = 'depth'
         depth(kmax) = 0.0d0
c         write(*,*)kmax
         do j = 2, kmax
            k = kmax - j +2
            depth(k-1) = r(1,kmax)-r(1,k-1)
         enddo
         do k = 1, kmax
            xx(k) = depth(k)
         enddo

      elseif( cvarx .eq. 'dpthh' )then
c..radial distance inward from photosphere, at zone centers
         cxlabel = 'dpthh'
         depth(kmax) = 0.0d0
         do j = 2, kmax
            k = kmax - j +2
            depth(k-1) = r(1,kmax)- 0.5d0*(r(1,k)+r(1,k-1))
         enddo
         do k = 1, kmax
            xx(k) = depth(k)
         enddo


      elseif( cvarx .eq. 'bind' )then
c..   grav. binding energy
         cxlabel = 'Binding'
         scre(1) = 0.0d0
         scrp(1) = 0.0d0
         do k = 2, kmax
            scre(k) = scre(k-1) + (dmh(k)/sol)*e(1,k)
            scrp(k) = scrp(k-1) - (dmh(k)/sol)
     1           *3.0d0*p(1,k)*v(1,k) 
            xx(k-1) = scre(k) + scrp(k) 
c            write(*,'(i5,1p8e12.3)')k,scre(k),scrp(k),scre(k)+scrp(k)
         enddo
         write(*,*)'BINDING ENERGY SUMS'
         k = kmax
         write(*,'(20x,2a5,8a12)')'k1','k2','E','3PV','B'
         write(*,'(a20,a5,i5,1p8e12.3)')
     1  'total star','1',k,scre(k),scrp(k),scre(k)+scrp(k)
         k = nvmax(3)
ccccccc:
         write(*,'(a20,2i5,1p8e12.3)')
     1  'envelope',k,kmax,scre(k)-scre(kmax),
     2	scrp(k)-scrp(kmax),scre(k)+scrp(k)-scre(kmax)-scrp(kmax)
cccccccccccccccccccc
      elseif( cvarx .eq. 'pturb' )then
         cxlabel = 'Turb. P'
         do k = 1, kmax
           if( v(1,k) .gt. 0.0d0 )then
             xx(k) = 0.5d0*h(k)**2/v(1,k)
           else
             xx(k) = 0.0d0
           endif
        enddo

      else

c..search for cvarx in nuclei list
         do i = 1, nnuc
            if( rightflush( cvarx, cnuc(i) ) )then
               ndum = i
               cxlabel = 'X('//cnuc(ndum)//')'
               if( ndum .lt. 1 .or. ndum .gt. ndim )then
                  write(*,'(a11)')cvarx//' error'
                  stop'getvar error in '
               endif
c..   ignores inner ghost zone
               do k = 2, kmax + 1
                  xx(k-1) = x(ndum,k) 
               enddo
               go to 1000
            endif
         enddo


c..   unsuccessful search for variable
         write(*,'(/a15,5x,a5,5x,a18)') 'THE SYMBOL', cvarx,
     1        'IS NOT RECOGNIZED'

         write(*,*)'The recognized symbols are:'
         write(*,'(10a6)')'lgt','lgd','lgp','rho','tem','pres',
     1        'lnp','u','sound','cour','h','xm','lgr','dr','lgdr',
     2        'r','R','rh','k','xl','lgkap','kappa','aloak',
     3        'lgrop','arad','dmh','dmi',
     4        'lgdmh','dmdk','d2mdk','d3mdk','dmerr',
     5        'xlrad','xlcnv','du','dnab','dnad','dnrad',
     6        'doux','r-ad','r-a-x','n-a-x','n-a',
     7        'ss','sa','sb','snu','snua',
     8        'snub','ss-nu','ln-lu','ln','lu','d2he4',
     9        'm/r4','r3','bvai2','Ri','vc/l2','l/r',

     *        'U','V','lcke','machc','pscal','gwfak','kec',
     1        'lgh1','metal','heli','hydr','cn','cno','cno6','cn7',
     2        'ysum','dlydm','yef','gam1','plin','hp','Hp','mfp',
     3         'depth','lgodp','opdep','dephh'

         write(*,*)'and nuclei:'
         write(*,'(10a6)')cnuc
         write(*,*)'are defined at present (see getvec.f)'
         stop'x label'
      endif

c..succesful exit from abundance search
 1000 continue

      if( cvarx .eq. 'hyd' .or. cvarx .eq. 'he4' )then
         xmax =  1.05
         xmin = -0.05
      else
         call minmax(xx,1,kmax,xmin,xmax)
      endif

      return
      end

      function rightflush(cin5,cout5)

c..remove trailing blanks from cin5, and put them in front
c..nondestructive, modified value is local
      implicit none
      integer*4 i,j
      character*5 cin5,cout5,cmod,cdummy
      character*1 blank,shift
      logical rightflush
      data blank/' '/
c----------------------------------------------------------
c..fill dummy variable cmod
      cmod = cin5
c..look for blanks and shift as needed

c      if( cout5 .eq. '  o18' )write(*,'(a1,5(a1,a1))')'_',
c     1     (cmod(i:i),'_',i=1,5)

      do j = 1, 4
         if( cmod(5:5) .eq. blank )then
            shift = cmod(5:5)
            do i = 1, 4
            cdummy(i+1:i+1)=cmod(i:i)
            enddo
            cdummy(1:1) = shift
            cmod=cdummy
         endif

c      if( cout5 .eq. '  o18' )write(*,'(a1,5(a1,a1))')'_',
c     1     (cmod(i:i),'_',i=1,5)

      enddo

         if( cmod(4:4) .eq. blank )then
            shift = cmod(4:4)
            cdummy(5:5) = cmod(5:5)
            do i = 1, 3
            cdummy(i+1:i+1)=cmod(i:i)
            enddo
            cdummy(1:1) = shift
            cmod=cdummy
         endif


c      if( cout5 .eq. '  o18' )write(*,'(a1,5(a1,a1))')'_',
c     1     (cmod(i:i),'_',i=1,5)


         if( cmod(3:3) .eq. blank )then
            shift = cmod(3:3)
            cdummy(5:5) = cmod(5:5)
            cdummy(4:4) = cmod(4:4)
            do i = 1, 2
            cdummy(i+1:i+1)=cmod(i:i)
            enddo
            cdummy(1:1) = shift
            cmod=cdummy
         endif


         if( cmod(2:2) .eq. blank )then
            shift = cmod(2:2)
            cdummy(5:5) = cmod(5:5)
            cdummy(4:4) = cmod(4:4)
            cdummy(3:3) = cmod(3:3)
            cdummy(2:2) = cmod(1:1)
            cdummy(1:1) = shift
            cmod=cdummy
         endif

c..test if shifted string is the same as the second argument
      if( cmod .eq. cout5 )then
         rightflush =  .TRUE.
      else
         rightflush = .FALSE.
      endif

      return
      end
