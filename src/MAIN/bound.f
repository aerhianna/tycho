      subroutine bound
c..   7-30-06
      implicit none

      include 'dimenfile'
      include 'cbound'
      include 'cconst'
      include 'comod'
      include 'compu'

      real*8 t87a,power,usuper,work,sumni,sumke
      integer*4 j
c---------------------------------------------------------------------

      if( ibnd .eq. 0 )then
c     inner boundary conditions
c     simple time independent form
         if( r(1,1) .le. 0.0d0 )then
            g(1)    = 0.0d0
            a(1)    = 0.0d0
            tl(2,1) = 0.0d0
            r(1,1)  = 0.0d0
         else
            g(1)    = grav*xm(1)/r(1,1)**2
            a(1)    = pi4*r(1,1)**2
            tl(2,1) = tl(1,1)
         endif
      elseif( ibnd .eq. 1 )then
         if( nfall .eq. 1 )then
c     falls inward on SN1987A time scale
            t87a = 2.3d0
            u(2,1) = - r(1,1)/t87a
            r(2,1) =  r(1,1) + u(2,1)*dth(2)
            du(1)  = - u(2,1)/t87a
         else
            u(1,1) = 0.0d0
            du(1)  = 0.0d0
            r(2,1) = r(1,1)
         endif
         dr(1) = r(2,1) - r(1,1)

      elseif( ibnd .eq. 2 )then
c..   piston explosion
         if( esuper .gt. 1.0d0 )then
            power  = esuper/3.0d-2
            a(1)   = pi4*r(1,1)**2
            g(1)   = grav*xm(1)/r(1,1)**2
            usuper = ( power * v(1,1)/a(1) )**0.3333333d0
c..   work done
            work = power * dth(2)
            u(2,1) = usuper
            r(2,1) = r(1,1) + u(2,1)*dti(1)
            esuper = esuper - work
            esuper = dmax1( 0.0d0, esuper)
            sumni = 0
            sumke = 0
            do j = 2, kk
               if( u(1,j) .gt. 1.0d0 )then
                  sumni = sumni + dmh(j)*x(10,j)*56.0d0
                  sumke = sumke + dmh(j)*0.5d0*u(1,j)**2
               endif
            enddo
            sumni = sumni/1.987d33
            sumke = sumke*1.0d-51
            write(*,'(a9,1p10e11.3)')"ESUPER ",esuper,r(2,1),
     1           usuper,u(1,1),power,work,a(1),g(1),sumni,sumke
            write(3,'(a9,1p10e11.3)')"ESUPER ",esuper,r(2,1),
     1           usuper,u(1,1),power,work,a(1),g(1),sumni,sumke
         else
            u(1,1) = 0.0d0
            du(1)  = 0.0d0
            esuper = 0
            usuper = 0
            power  = 0
            work   = 0
            a(1)   = pi4*r(1,1)**2
            g(1)   = grav*xm(1)/r(1,1)**2
            r(2,1) = r(1,1)
            sumni  = 0
            sumke  = 0
            do j = 2, kk
               if( u(1,j) .gt. 1.0d8 )then
                  sumni = sumni + dmh(j)*x(10,j)*56.0d0
                  sumke = sumke + dmh(j)*0.5d0*u(1,j)**2
               endif
            enddo
            sumni = sumni/1.987d33
            sumke = sumke*1.0d-51
            write(*,'(a9,1p10e11.3)')"ESUPER ",esuper,r(2,1),
     1           usuper,u(1,1),power,work,a(1),g(1),sumni,sumke
            write(3,'(a9,1p10e11.3)')"ESUPER ",esuper,r(2,1),
     1           usuper,u(1,1),power,work,a(1),g(1),sumni,sumke
         endif
         dr(1) = r(2,1) - r(1,1)
      else
         write(*,*)'ibnd is ', ibnd, ' in bound'
         stop'bound'
      endif

      return
      end
