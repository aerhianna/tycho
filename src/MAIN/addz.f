c
c
c     
      subroutine addz(scr,k,m)
c     
c     
c     

      implicit none

      include 'dimenfile'
      include 'czone'
      include 'cgtr.h'
      include 'cburn'
      include 'cconst'
      include 'comod'
      include 'compu'
      include 'cenv'
      include 'csurface'
      include 'cnabla'
      include 'ceoset'
      include 'cgtintp'

      real*8    scr(kdm)
      real*8    scx(ndim,2)

      real*8    fak,dfak,del
      real*8    pbar,tbar,pave
      real*8    pp,pm,rk,vvm,vvp,deltt,tm,tp,dnabb

      integer*4 i,j,k,m,n
      
c----------------------------------------------------------------
c..   found problem zones via lack of hydrostatic equilibrium
c..   twiddle radius to improve HSE condition, keeping 
c..   abundance gradients the same
c     
      
      write(*,'(a25,i5,a5,i5,5(a8,1pe11.3))')
     1     'REZONE: ADDZ: bnd.',k,'of',njj,
     1     'P',p(m,k),'error',scr(k),
     2     'vcnv(k)',h(k), 
     4     'He(k)',  xa(nnuc)*x(nnuc,k), 
     3     'He(k+1)',xa(nnuc)*x(nnuc,k+1)

      if( k .lt. njj-1 )then
c..   extra io for bad zones
         if( abs(scr(k) ) .gt. 1.0d-2 )then
            write(*,'(a20,i5,1pe12.3)')'Excessive error ',k,scr(k)
            write(*,'(2a5,2a14,8a12)')'k','ic','M/sol','dM','err(HSE)',
     1           'P','X(He4)','r','t','v'
            if( k .le. 2 )then
               do i = 2,4
                  write(*,'(2i5,1p2e14.7,1p8e12.3)')i,ic(i),xm(i)/sol,
     1                 dmh(k)/sol,scr(i),p(m,i),
     2                 x(nnuc,i)*xa(nnuc),r(m,i),t(m,i),v(m,i)
               enddo
            elseif( k .ge. njj-2 )then
               do i = njj-4,njj
                  write(*,'(2i5,1p2e14.7,1p8e12.3)')i,ic(i),xm(i)/sol,
     1                 dmh(k)/sol,scr(i),p(m,i),
     2                 x(nnuc,i)*xa(nnuc),r(m,i),t(m,i),v(m,i)
               enddo
            else
               do i = k-2,k+2
                  write(*,'(2i5,1p2e14.7,1p8e12.3)')i,ic(i),xm(i)/sol,
     1                 dmh(k)/sol,scr(i),p(m,i),
     2                 x(nnuc,i)*xa(nnuc),r(m,i),t(m,i),v(m,i)
               enddo
            endif
            write(*,*)'k ',k,' kk ',kk,' njj ',njj,' m ',m
c            stop'addz in rezone'
         endif

c..save all abundances for zones k and k+1
         do i = 1,2
            do j = 1,nnuc+1
               scx(j,i) = x(j,k+i-1)
            enddo
         enddo

c..   estimate pressures from HSE and neighboring zones
         pm = p(m,k-1) -grav*xm(k-1)/(pi4*r(m,k-1)**4)*dmi(k-1)

         write(*,'(i5,4(a6,1pe12.4))')k,'pm',pm,'p(m,k)',p(m,k),
     1        'err',pm/p(m,k)-1.0d0

         if( k .lt. njj-1)then
            pp = p(m,k+2) +grav*xm(k+1)/(pi4*r(m,k+1)**4)*dmi(k+1)
         elseif( k .eq. njj-1)then
c..boundary value
            pp = p(m,k+2) +grav*xm(k+1)/(pi4*r(m,k+1)**4)*dmi(k+1)
            write(*,*)'boundary value'
         else
            write(*,*)k
            stop'k error in addz'
         endif

c         write(*,'(4(a12,1pe12.4))')'pp',pp,'p(m,k+1)',p(m,k+1),
c     1        'err',pp/p(m,k+1)-1.0d0


c..   compute radius from HSE with estimated pressures
         if( pm .gt. pp )then
            rk = ( grav* xm(k)*dmi(k)/(pi4*(pm-pp)))**0.25d0
         else
            rk = 0
            write(*,'(a25,i5,1p8e12.3)')'addz: rk: P inversion ',k,
     1           pm-pp,pm,pp
            write(*,'(4(a12,1pe12.4))')'pp',pp,'p(m,k+1)',p(m,k+1),
     1           'err',pp-p(m,k+1)
            write(*,*)'mass/sol ',xm(k)/sol,' dmh/sol ',dmh(k)/sol

            stop'r error in addz'
         endif

c         rk = sqrt( rk*r(m,k) )
cccccccccccccccc

c         write(*,'(4(a12,1pe12.4))')'rk',rk,'r(m,k)',r(m,k),
c     1        'err',rk/r(m,k)-1.0d0
c         write(*,'(4(a12,1pe12.4))')'r(m,k-1)',r(m,k-1),
c     1        'r(m,k+1)',r(m,k+1),'r+-r-',r(m,k+1)-r(m,k-1)

c            write(*,'(a20,1p8e12.3)')'dp',pm-pp,
c     1        grav*xm(k)*dmi(k)/(pi4*rk**4)
cccccccccccccccccccccccc

c..twiddle r(m,k)
         r(m,k) = rk


c..   compute specific volume (1/rho) from mass conservation
         vvm = pi43*( rk**3 - r(m,k-1)**3 )/dmh(k)
         vvp = pi43*( r(m,k+1)**3 - rk**3 )/dmh(k+1)

c         write(*,'(4(a12,1pe12.4))')'vvm',vvm,'v(m,k)',v(m,k),
c     1        'err',vvm/v(m,k)-1.0d0
c         write(*,'(4(a12,1pe12.4))')'vvp',vvp,'v(m,k+1)',v(m,k+1),
c     1        'err',vvp/v(m,k+1)-1.0d0


         if( vvm .lt. 0.0d0 .or. vvp .lt. 0.0d0 )stop'addz'
ccccccccccccccccccccccc

c..   revise specific volumes
         v(m,k)   = vvm
         v(m,k+1) = vvp

         pbar     = 0.5d0*( p(m,k+1) + p(m,k) )
         tbar     = 0.5d0*( t(m,k+1) + t(m,k) )
         deltt = (t(m,k+1) - t(m,k))/tbar*pbar/(p(m,k+1)-p(m,k))
c         write(*,'(4(a12,1pe12.4))')'dlnT/dlnP',deltt,
c     1        'dnab(k)',dnab(k),'err',deltt-dnab(k)
c         pave =  0.5d0*( p(m,k+1) + p(m,k-1)  - grav/pi4*(
c     2        xm(k-1)*dmi(k-1)/r(m,k-1)**4  
c     3        - xm(k)*dmi(k  )/r(m,k)**4  )
c     4        )
c         write(*,'(4(a12,1pe12.4))')'pave0',
c     1        pave, 'p(m,k)',p(m,k),'err',pave/p(m,k)-1.0d0


c..   smooth X (fak = 1 is unmixed, fak = 0.5 is uniform at mass average)
c..   check convective state
         if( h(k) .gt. 1.0d0 )then
c..convective, mass average mix
            fak = 0.5d0
            do j = 1, nnuc+1
               x(j,k) = (fak*dmh(k)*scx(j,1) 
     1              + (1.0d0-fak)*dmh(k+1)*scx(j,2))/
     2              (fak*dmh(k)+(1.0d0-fak)*dmh(k+1))
               x(j,k+1) = ((1.0d0-fak)*dmh(k)*scx(j,1) 
     1              + fak*dmh(k+1)*scx(j,2))/
     2              ((1.0d0-fak)*dmh(k)+fak*dmh(k+1))
            enddo
         else
c..radiative 
            fak = 1.0d0
            do j = 1, nnuc+1
               x(j,k)   = scx(j,1) 
               x(j,k+1) = scx(j,2)
            enddo
         endif


c         write(*,'(a20,i5,1p8e12.3)')'before state k',k,
c     1        v(m,k),p(m,k),r(m,k), v(m,k+1),p(m,k+1)
ccccccccccccccccccccc

c..   iterate on P for T, using new V and X, with estimate P's
c..   as the targets, for zones k and k+1
         do i = 1, 20

            call state(k,k,m)
      
            fak = p(m,k) - pm 
            dfak = pt(k)
            del = - fak/dfak
            t(m,k) = t(m,k) + del
            if( abs(fak) .lt. 1.0d-12*pm )goto 100
         enddo
         write(*,*)'addz: error in iteration, k=',k
         stop'addz'
 100     continue

c         write(*,'(a20,i5,1p8e12.3)')'before state k',k+1,
c     1        v(m,k),p(m,k),r(m,k), v(m,k+1),p(m,k+1)
ccccccccccccccccccccc

         do i = 1, 20

            call state(k+1,k+1,m)
      
            fak = p(m,k+1) - pp
            dfak = pt(k+1)
            del = - fak/dfak
            t(m,k+1) = t(m,k+1) + del
            if( abs(fak) .lt. 1.0d-12*pp )goto 101
         enddo
         write(*,*)'addz: error in iteration, k=',k
         stop'addz'
 101     continue

c         write(*,'(4(a12,1pe12.4))')'pm',pm,'p(m,k)',p(m,k),
c     1        'err',pm-p(m,k)
c         write(*,'(4(a12,1pe12.4))')'pp',pp,'p(m,k+1)',p(m,k+1),
c     1        'err',pp-p(m,k+1)

c         pave =  0.5d0*( p(m,k+1) + p(m,k-1)  - grav/pi4*(
c     2        xm(k-1)*dmi(k-1)/r(m,k-1)**4  
c     3        - xm(k)*dmi(k  )/r(m,k)**4  )
c     4        )
c         write(*,'(4(a12,1pe12.4))')'pave1',
c     1        pave, 'p(m,k)',p(m,k),'err',pave/p(m,k)-1.0d0

         pbar  = 0.5d0*( p(m,k+1) + p(m,k) )
         tbar  = 0.5d0*( t(m,k+1) + t(m,k) )
         dnabb = (t(m,k+1) - t(m,k))/tbar*pbar/(p(m,k+1)-p(m,k))

c         write(*,'(5(a12,1pe12.4))')
c     1        'dnab(k-1)',dnab(k-1),
c     1        'dnab(k)'  ,dnab(k),
c     2        'dnab(k+1)',dnab(k+1),
c     3        'dnabb'    ,dnabb,
c     1        'err',dnabb-dnab(k)

      elseif( k .ge. njj-1 )then
         do i = njj-3, njj+1
            write(*,'(i5,1p3e15.7,1p8e12.3)')i,r(m,i),qqo(i),qqn(i),
     1           qqo(i)-qqn(i),
     1           p(m,i),t(m,i),1.0d0/v(m,i),x(nnuc,i)
         enddo
         write(*,'(a5,3a15,8a12)')'k','r','qqo','qqn','qqo-qqn',
     1        'p','T','rho','Yhe4'
         write(*,'(i5,1p8e12.3)')nvmax(3),vp(nvmax(3),3),
     1        vtem(nvmax(3),3),vrho(nvmax(3),3)

      else

         write(*,*)k,njj
         stop'addz error ddd'

      endif

      return
      end

