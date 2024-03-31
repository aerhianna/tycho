      subroutine fnab(dnab,dnad,dnrad,st,sv,sp,gravity,spt,spv,
     1     set,sev,sak,svel,doux,alphaml,uuml,modec,k) 

c..   1-15-2008

c..real*8 function fnab for envelope integration

      implicit none

c..implements inefficient convection by Bohm-Vitense mixing length
c..theory as in Kippenhahn and Weigert, p. 51, section 7.2
c..called in envel.f and cinit.f
c..only uses scalars
c..rechecked versus KW, wda, 12-17-05
c..revised geometry factor uuml 8-24-07

      include 'cenv'
      include 'cfnab'
      include 'cconst'

      integer*4 j,modec,k
      
      integer*4 dumflag
      data dumflag/0/

      real*8    w,dnad,dnrad,st,sv,sp,gravity,spt,spv,set,sev,sak
      real*8    hp,hml,cp,delta,eratio,fact,fact2
      real*8    tiny,fun,dfun,dpsi,dnab,svel,cvcon
      real*8    doux,alphaml,uuml
      real*8    aa, bb, cc, dnabsc, dnled, krad, discr
      real*8    alphasc, svelsc
c---------------------------------------------------------------------
      
      dnled = dnad + doux
      
      if( dnrad .le. dnad ) then
c..   RADIATIVE GRADIENT
         dnab = dnrad
         svel = 0.0d0
         
      else if( dnrad .gt. dnad) then
c..   CONVECTIVE AND SEMI-CONVECTIVE GRADIENT
         
         w = dnrad - dnad      !ignore doux b/c of fast mixing
         
         
         hp     = sp*sv/gravity
         hml    = alphaml * hp
         cp     = set - spt*(sp + sev)/spv
         delta  = -st*spt/(sv*spv)
         eratio = arad*st**3*sv/cp
         fact   = 3.0d0*crad * sv /(hml**2 * sak)
         fact2  = dsqrt( 8.0d0 * hp/( gravity * delta) )
         uu     = eratio * fact * fact2 * uuml
c         write(*,*)eratio,fact,fact2,uu

c..   uuml is scaling for radiative damping in mlt


c..   iterate for solution of cubic equation in psi
         coef = 8.0d0/(9.0d0 * uuml)
         psi  = uu + (coef*uu*w)*0.33333d0
         tiny = 1.0d-9
         
         do j = 1, 300
            fun  = (psi-uu)**3+coef*uu*(psi**2-uu**2-w)
            dfun = 3.0d0*(psi-uu)**2 + 2.0d0*coef*uu*psi
            dpsi = - fun/dfun
            
            if( dpsi/psi .gt. 3.0d0 )then
               psi = 3.0d0 * psi
            elseif( dpsi/psi .lt. -0.75d0 )then
               psi = 0.25d0*psi               
            else
               psi = psi + dpsi
            endif
            
            if(j.gt.290 )write(*,'(2i5,1p8e12.4)')j,k,psi,fun,dfun,
     .           dpsi,dpsi/psi
            
            if(abs(dpsi).lt.tiny )goto 100
         enddo
         
         print*,'ERR(fnab): nonconvergence in fnab.'
         write(*,'(2a5,8a12)')'j','k','psi','fun','dfun','dpsi',
     .        'alphaml','uuml'
         write(*,'(2i5,1p8e12.3)')j,k,psi,fun,dfun,dpsi,alphaml,uuml
         write(*,'(5(a8,1pe12.3))')'dnab',dnab,'dnrad',dnrad,
     .        'dnad',dnad,'doux',doux,'r-a-x',w
         write(*,'(5(a8,1pe12.3))')'T',st,'V',sv,
     .        'P',sp,'g',gravity,'ak',sak
         write(*,'(5(a8,1pe12.3))')'sign',sign(1.0d0,dpsi),
     .        'uu',uu,'coef',coef
         stop' nonconvergence in function fnab.'
 100     continue
         
         dnab  = psi**2 - uu**2 + dnad
         cvcon = sqrt( delta * sp * sv * 0.125d0) * alphaml
         svel  = cvcon*( coef*uu*( dnrad - dnab ) )**0.33333333d0
       
         
         if(dnrad .lt. dnled) then !SEMICONVECTIVE GRADIENT
            
c     !!! USE SCHWARZCHILD VALUES FROM ABOVE !!!
            
         endif
         
      endif
      

ccccccccccccccccccccccccccccc
      if( dnab .gt. 10.0d0 )then
         write(*,*)'MSG(fnab.f): large dnab'
         write(*,'(2a5,12a12)')'k','j','dnab','svel','cvcon',
     1    'uu','dnrad','dnad','psi','doux','w'
         write(*,'(2i5,1p12e12.3)')k,j,dnab,svel,cvcon,uu,
     1     dnrad,dnad,psi,doux,w
c          write(*,'(8a12)')'sp','sv','sound'
c          write(*,'(1p8e12.3)')sp,sv,dsqrt(sp*sv)
c        stop'fnab'
      endif
     
      return

      end


