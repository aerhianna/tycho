c
c
c
      subroutine cinit(kin,ktop,nc)
c----------------------------------------------------
c..   defines nablas for convection
c..   sets flags for type of convection
c..   calculates convective velocities, buoyancy frequencies,
c..   boundary conditions

c..   requires that fitenv.f has already run (uses envelope variables)

      implicit none

c..   last revised 8-4-07

      include 'dimenfile'
      include 'comod'
      include 'compu'
      include 'cnabla'
      include 'cconst'
      include 'cburn'
      include 'cenv'
      include 'ceoset'

c     real*8 douxlim
c     parameter( douxlim = 1.0d-2 )

c..   kin and ktop are zone index k limits
c..   k1 and k2 are dummy indices for dnad calculation
c..   nc is time step index
c..   l is step number index for this run
c..   i is working index for loops
      integer*4 i,n,k1,k2,nc,kin,ktop,j,k
      integer*4 istep0,istep1,nrem,nex

c..   local working arrays
c..   ia(k) saved entry from nc=1 call (from cmix.f)
c..   ib(k) working value for this iteration (in hstat.f)
c..   ic(k) global convection flag
c..   ip(k) local value, saved from previous iteration on it
c..   ik(k) global flag for changed zones

c..   output: ic(i), doux(i),dnab(i),dnad(i),dnrad(i),dnabv(i),
c..   h(i),nsqr(i),nrczones

      integer*4 ib(kdm),ip(kdm)
c..   number of changed zones
      integer*4 kchange,ichange

      real*8 tbar,delt,vbar,dbar,pbar,delp,doufak,beta
      real*8 dvdt,dltdlp,akbar,tledd,tb4,tlum,tle,ombeta
      real*8 cpm,cpp,cpbar,tlc,ptbar,pvbar,etbar,evbar,gbar
      real*8 dnadbar,dnradbar,svel
      real*8 dnradold, dnabold, dnadold
      real*8 ybar(kdm),deldt(kdm),deldy(kdm),deld
      
      real*8 ykm,ykp

      real*8 fak,xmconv(kdm),rconv(kdm),sumq(kdm),sumke(kdm),sumtau
      real*8 rbot,rtop,xmbot,xmtop,pbot,ptop,pratio(kdm)
      real*8 entrczt(kdm),entrczb(kdm),nsqrb(kdm),nsqrt(kdm)
      real*8 sum1(kdm),sum2(kdm)

      real*8 dnabtiny, hsmall
      integer*4 ishear


c..   moved to cnabla
      integer*4 kcbeg,kcend,ibold
      
      save


c-----------------------------------------------------------
c     write(*,*)'MSG(cinit.f): ENTERING ',kin,ktop,nc
c--------------------------------------------
c..   no convection over-ride
      if( modec .lt. 0 )then
         return
      endif
      
      print*,'MSG(cinit): executing. l,model =',l,model
      

c..   SOME FLOW CONTROL FLAGS
      
      ishear = 1                !modify conv spd by shear/nsqr (1=yes)
      
      
      
c     mass exccess and mass fractions over zones
      do k = kin, ktop
c     use -1 for mass EXCESS
        zeta(k) = -1.0d0
c     Ye = charge per nucleon
        x(nnuc+1,k) = 0.0d0
        do j = 1, nnuc
          zeta(k) = zeta(k) + x(j,k)*xa(j)/nuca(j)
          x(nnuc+1,k) = x(nnuc+1,k) + x(j,k)*nucz(j)/nuca(j)
        enddo
      enddo      
c..   define current values of nablas, store old values
      do i = kin,ktop
c     
         dnrado(i) = dnrad(i)   !store old values
         dnado (i) = dnad (i)   !      "
         douxo (i) = doux (i)   !      "
c     
         tbar = 0.5d0*(t(nc,i+1) + t(nc,i))
         delt =        t(nc,i+1) - t(nc,i)
         vbar = 0.5d0*(v(nc,i+1) + v(nc,i))
         dbar = 1.0d0/vbar
         deld = ( 1.0d0/v(nc,i+1) - 1.0d0/v(nc,i))
         pbar = 0.5d0*(p(nc,i+1) + p(nc,i))
         delp =        p(nc,i+1) - p(nc,i)
         dely(i) = 0
         ybar(i) = 0
c..   includes free electron contribution
c     update Ye
         x(ndim,i) = 0
         do n = 1, nnuc
           x(ndim,i) = x(ndim,i) + x(n,i)*nucz(n)/nuca(n)
         enddo
         x(ndim,i+1) = 0
         do n = 1, nnuc
           x(ndim,i+1) = x(ndim,i+1) + x(n,i+1)*nucz(n)/nuca(n)
         enddo
c     total free particle number
         ykm = x(ndim,i)
         do n = 1, ndim-1
            ykm = ykm + x(n,i)/xa(n)
         enddo
         ykp = x(ndim,i+1)
         do n = 1, ndim-1
            ykp = ykp + x(n,i+1)/xa(n)
         enddo
         dely(i) = ykp - ykm
         ybar(i) = 0.5d0*( ykp + ykm )
c     
c..   call from cinit.f in rezone.f has zero P,V at kk+1
c     
         if( i .eq. ktop .and. i .ne. kin )then
c..   edge of grid  
            beta = rgas*
     1           t(1,i  )*ykm/(p(1,i  )*v(1,i  ) )
         else
c..   this form reduces roundoff
            beta = rgas*
     1           ( t(1,i  )*ykm/(p(1,i  )*v(1,i  ) )
     2           + t(1,i+1)*ykp/(p(1,i+1)*v(1,i+1) ) )*0.5d0
         endif

         
c..   TODO: THIS EXPRESSION FOR doufak IS ONLY CORRECT 
c..   FOR IDEAL GAS w/RAD PRESS, RIGHT?
         
         doufak = beta/(4.0d0 - 3.0d0*beta)
         
         if( dely(i) .ge. 1.0d-14 .and. 
     1        r(1,i) .gt. 2.0d0*r(1,2) .and.
     2        r(1,1) .le. 1.0d0 )then
c              doux(i) = doufak * dely(i) * 
c     1           pi4*r(1,i)**4*0.5d0*( p(1,i)/ykm + p(1,i+1)/ykp )/
c     1           ( grav*xm(i)*dmi(i) )
c     
            doux(i) = -doufak*(p(1,i+1)+p(1,i))/(ykp+ykm)*
     .           dely(i)/(p(1,i+1)-p(1,i))
            
         else
c..   expand about origin
c..   this form reduces roundoff
            doux(i) = doufak * dely(i)
     1           * 3.0d0 * r(1,i)/( grav*dmi(i) ) 
     2           * ( p(1,i)*v(1,i)/ykm +p(1,i+1)*v(1,i+1)/ykp)*0.5d0
         endif
         

c..   defines conventional nablas ( gradients = d ln t / d ln p )
c..   actual gradient
         dnab(i) = pbar*delt/( tbar*delp )
c..   dlnd/dlnT
         deldt(i) = tbar*deld / (dbar*delt)
c..   dlnd/dlny
         if( dely(i) .eq. 0.0d0 )then
            deldy(i) = 0.0d0
         else
            deldy(i) = ybar(i)*deld / (dbar*dely(i))
         endif
c..   adiabatic gradient
         dnad(i)  = 0.0d0
         if( i .lt. kk )then
            do k1 = 1, 2
               k2 = i -1 + k1
               dvdt    = - et(k2)/( p(nc,k2) + ev(k2) )
               dltdlp  =
     1              p(nc,k2)/( t(nc,k2)*( pt(k2) + pv(k2)*dvdt ))
               dnad(i) = dnad(i) + 0.5d0*dltdlp
            enddo
         else
            dvdt    = - et(kk)/( p(nc,kk) + ev(kk) )
            dltdlp  =
     1           p(nc,kk)/( t(nc,kk)*( pt(kk) + pv(kk)*dvdt ))
            dnad(i) = dltdlp
         endif

c..   kk+1 contains photospheric values, 
c..   so use kk as boundary value rather than average

c..   radiative gradient
         akbar  = 0.5d0*( ak(i) + ak(i+1) )
         tledd  = pi4 * crad * grav * xm(i) / akbar
         tb4    = (0.5d0*( t(1,i+1)**4.0d0/p(1,i+1) 
     1            + t(1,i)**4.0d0/p(1,i) ))
         ombeta = arad*tb4/(3.0d0)
         tle    = tledd * ombeta * 4.0d0 *dnab(i)

         f(i) = tle/a(i)

c     NOT SURE WHERE THIS LIMIT IS NEEDED
c     USING tl(n,ci) FOR LUMINOSITY IN THE MEANTIME
c     ...................................
c      tlum = min(tl(nc,i), 2.0d0*tl(1,i))
         tlum = tl(nc,i)

c..   b, the convective flux, is purely diagnostic
         b(i)     = (tlum - tle)/a(i)
         dnrad(i)  = tlum /(tledd * ombeta * 4.0d0)
      
         
c..   specific heat at constant pressure
         cpm = et(i) - ( p(nc,i) + ev(i) ) * pt(i) / pv(i)
         cpp = et(i+1) - ( p(nc,i+1) + ev(i+1) )*pt(i+1)/pv(i+1)
         cpbar = 0.5d0 * ( cpm + cpp )

c..   convective limiting luminosity
         tlc      = sound(i) * a(i) * tbar * cpbar / vbar
         edobv(i) = tle * 0.25d0 / tlc
cccccccccccccccccccccccccccccccccccccccc

c     if( tbar .lt. 1.0d6 )then
c..   defines boehm-vitense effective nabla
         ptbar = 0.5d0*( pt(i) + pt(i+1) )
         pvbar = 0.5d0*( pv(i) + pv(i+1) )
         etbar = 0.5d0*( et(i) + et(i+1) )
         evbar = 0.5d0*( ev(i) + ev(i+1) )
         gbar     = g(i)
         dnadbar  = dnad(i)
         dnradbar = dnrad(i)

c..   modec = 1 for schwarzschild, ddnab = dnab - dnad
c..   modec = 2 for ledoux, 3 for richardson, 
c..   ddnab = dnab - dnad - doux in fnab

         call fnab(dnabv(i), dnadbar,dnradbar,tbar,vbar,
     .        pbar,gbar,ptbar,pvbar,etbar,evbar,akbar,h(i),
     .        doux(i),alphaml,uuml,modec,i)

c         dnabv(i) = pbar*delt/( tbar*delp )
         if(dnrad(i).le.dnad(i)) then

         endif

         
         if( dnabv(i) .gt. 1.0d1 )then
            print*,'ERR(cinit): dnabv too big: xm = ',xm(i)/sol
            print*,'ERR(cinit): (new) dnrad,dnad,doux =',
     .           dnrad(i),dnad(i),doux(i)
            print*,'ERR(cinit): (old) dnrad,dnad,doux =',
     .           dnrado(i),dnado(i),douxo(i)
            print*,'ERR(cinit): ak(ii),ii=(i-2,i-1,i,i+1,i+2):',
     .           ak(i-2),ak(i-1),ak(i),ak(i+1),ak(i+2)
            print*,'ERR(cinit): tl(ii),ii=(i-2,i-1,i,i+1,i+2):',
     .           tl(nc,i-2),tl(nc,i-1),tl(nc,i),tl(nc,i+1),tl(nc,i+2)
            print*,'ERR(cinit): p(ii),ii=(i-2,i-1,i,i+1,i+2):',
     .           p(nc,i-2),p(nc,i-1),p(nc,i),p(nc,i+1),p(nc,i+2)
            print*,'ERR(cinit): t(ii),ii=(i-2,i-1,i,i+1,i+2):',
     .           t(nc,i-2),t(nc,i-1),t(nc,i),t(nc,i+1),t(nc,i+2)
            print*,'ERR(cinit): tlum, tle, tl(nc,i)=',tlum,tle,tl(nc,i)
            print*,'ERR(cinit): dnrad(tl) = ',
     .           tl(nc,i)/( tledd * ombeta * 4.0d0 )
c     print*,'ERR(cinit): 3.0/16*pi*G*a*c = ', 
c     .           3.0/(4.0*pi4*grav*arad*crad)
c            stop'cinit'           
         endif


c..   (Kippenhan & Weigert, p.39)
c..   calculate  Brunt-Vaisala frequency**2
         if( i .lt. kk )then
            nsqr(i) = ( -dnrad(i)+ dnad(i) +doux(i) )*
     1           (grav*xm(i)/r(nc,i)**2)**2/
     1           (0.5d0*( p(nc,i)*v(nc,i) + p(nc,i+1)*v(nc,i+1) ))
     1           *(v(nc,i)-v(nc,i+1))*(t(nc,i)+t(nc,i+1))*0.5d0/
     1           ((t(nc,i)-t(nc,i+1))*(v(nc,i)+v(nc,i+1))*0.5d0)
c            write(*,*)"N^2 ", nsqr(i)
c            nsqr(i) = 1.0d0/v(nc,i)*(v(nc,i)*(v(nc,i)-v(nc,i+1))/
c     1           (r(nc,i)-r(nc,i+1))-(p(nc,i)-p(nc,i+1))/
c     1           (r(nc,i)-r(nc,i+1))/(p(nc,i)*gamma1(i)))*
c     1           grav*xm(i)/r(nc,i)**2.0d0
            j = 0
c            write(*,*)i,nsqr(i)

         elseif( i .eq. kk )then
c..   at join
            j = nvmax(3)
            nsqr(i) = ( -vnab(j,3)+vnad(j,3))*
     .           (grav*(xm(kk+1)+vm(j,3))/vr(j,3)**2)**2
     .           *vrho(j,3)/vp(j,3)
         else
            stop'error in nsqr calculation in CINIT.f'
         endif
      enddo
      


c..   define edge values
      dnab(kk)  = dnab(kk-1)
      dnad(kk)  = dnad(kk-1)
      dnrad(kk) = dnrad(kk-1)
      doux(kk)  = 0.0d0
      dnabv(kk) = dnabv(kk-1)
      
c     astab(kk) = astab(kk-1)
      
c..   ------------------------------
c..   NABLAS HAVE NOW BEEN REDEFINED
c..   ------------------------------
      
      



      

c..   ----------------------------
c..   SETUP LOGIC FLAGS FOR MIXING
c..   ----------------------------
      
      if( modec .ne. 1 .and. modec .ne. 2 .and. modec .ne. 3 )then
         write(*,*)'ERR(cinit): modec = ',modec
         stop'CINIT: modec error'
      endif
      
      
c..   DEFINE WORKING ib ARRAY
c..   **NOW BASED ON MODEL GRADIENTS (dnrad,dnad,dndoux)
c..   **NOT RESULTS OF FNAB (i.e, NOT dnabv)
      
c     dnabtiny = 1.d-3
      do i = kin,ktop-1
         if(dnrad(i)-(dnad(i)+doux(i)).gt.1.0d-15)then
            ib(i) = 1           !CONVECTIVE (MLT) TEMP GRADIENT
         elseif( (dnrad(i) -dnad(i))  .gt.1.0d-15 )then
            ib(i) = 2           !SEMICONVECTIVE TEMP GRADIENT
         else
            ib(i) = 0           !RADIATIVE TEMP GRADIENT
            if(abs(h(i)).gt.0.d0) then
               if(nsqr(i) .le. 0.0d0)then
               print*,'ERR(cinit):i,xm,ib,h,dnabla=',
     .              i,xm(i)/sol,ib(i),h(i),(dnrad(i)-dnad(i))
               endif
            endif
         endif
      enddo


c..   BOUNDARY FLAGS: 
c..   1. in order to count convective core properly
c..   2. and zero out old values due to outer edge of mesh moving
      ib(1) = 0
      ib(kk+nvmax(3)) = 0






c..   CALCULATE A GRADIENT RICHARDSON NUMBER OF SORTS
c..   AT CONVECTIVE BOUNDARIES BASED ON CONVECTION SPEED
c..   
      hsmall = 1.d2             !units: [cm/s]
      do i = kin,ktop
         shear(i) = 0.0d0        !only at conv bndry, otherwise 0
c         write(*,*)"rotational velocity ", rotshear(i)
      

         if(ishear .eq. 1) then            
c            if( ib(i) .ne. 0 )then
c            if( ib(i) .ne. 1 .and. (h(i+1)-h(i-1)) .ge. hsmall)then
             if(nsqr(i) .le. 0.0d0 )then
c .and.(h(i+1)-h(i-1)) .ge. hsmall)then
               shear(i) = 0.5d0*(
     .              (h(i+1)-h(i))/(r(nc,i+1)-r(nc,i))**2.0d0 +
     .              (h(i)-h(i-1))/(r(nc,i)-r(nc,i-1))**2.0d0 )
c     .       +  (rotshear(i+1)-rotshear(i))/(r(nc,i+1)-r(nc,i))**2.0d0 +
c     .          (rotshear(i)-rotshear(i-1))/(r(nc,i)-r(nc,i-1))**2.0d0 )
c..   satisfy Schwarzschild convection but zero convective velocity
c..   Froude # = 1/ Richardson # = shear/nsqr avoids 1/zero
               if( shear(i) .ge. 4.0d0*dabs(nsqr(i)) )then
c     mix on convective time scale of adjacent convective zones
                  h(i) = 1.0d0*( h(i+1) + h(i-1) )*0.5d0
c              h(i) = h(i)+dabs((rotshear(i+1) - rotshear(i-1))*0.5d0)
               elseif( shear(i) .ge. nsqr(i) )then
c     mix on 5 times the convective time scale (fast mixing)
                  h(i) = ( h(i+1)+h(i-1))/2.0d0*0.2d0
c                  write(*,*)'before ',i,  h(i)
c              h(i) = h(i)+dabs(rotshear(i+1)-rotshear(i-1))
c     1             /(r(nc,i+1)-r(nc,i))
c/(0.006*(t(nc,i))**0.5d0))
c                  h(i) = 0.0d0
c                  write(*,*)'after ', i, h(i)
               else
c     erosion rate depends upon stratification energy and convective power
c     (differential Richardson number)
                  h(i)=0.5d0*shear(i)/nsqr(i)*(h(i+1)+h(i-1)) 
c                  h(i) = 0.0d0
               endif
            endif
c            endif
         endif
      enddo
      
      do i = ktop-1,kin
c         shear(i) = 0.0d0        !only at conv bndry, otherwise 0

         if(ishear .eq. 1) then            
c            if( ib(i) .ne. 0 )then
c            if( ib(i) .ne. 1 .and. (h(i+1)-h(i-1)) .ge. hsmall)then
             if(nsqr(i) .le. 0.0d0 )then
c                write(*,*)"rotational velocity ",rotshear(i)
c .and.(h(i+1)-h(i-1)) .ge. hsmall)then
c      shear(i) = 0.5d0*(
                shear(i) = 0.5d0 * (
     .          (h(i+1)-h(i))/(r(nc,i+1)-r(nc,i))**2.0d0 +
     .              (h(i)-h(i-1))/(r(nc,i)-r(nc,i-1))**2.0d0 )
c     .          (rotshear(i+1)-rotshear(i))/(r(nc,i+1)-r(nc,i))**2.0d0 +
c     .          (rotshear(i)-rotshear(i-1))/(r(nc,i)-r(nc,i-1))**2.0d0 )
c..   satisfy Schwarzschild convection but zero convective velocity
c..   Froude # = 1/ Richardson # = shear/nsqr avoids 1/zero
               if( shear(i) .ge. 4.0d0*dabs(nsqr(i)) )then
c     mix on convective time scale of adjacent convective zones
                  h(i) = max(h(i),1.0d0*( h(i+1) + h(i-1) )*0.5d0)
c                  h(i) = max(h(i),(rotshear(i+1) + rotshear(i-1))*0.5d0)
               elseif( shear(i) .ge. dabs(nsqr(i)) )then
c     mix on 5 times the convective time scale (fast mixing)
                  h(i) = max(h(i),( h(i+1)+h(i-1))/2.0d0*0.2d0)
c             h(i) = max(h(i),(rotshear(i+1)+rotshear(i-1))/2.0d0*0.2d0)
c                  h(i) = 0.0d0
               else
c     erosion rate depends upon stratification energy and convective power
c     (differential Richardson number)
                  h(i)=max(h(i),0.5d0*shear(i)/nsqr(i)*(h(i+1)+h(i-1))) 
c                  h(i) = max(h(i),0.0d0)
               endif
c            endif
            endif
         endif
      enddo
c..   
c..   CHECK JOIN BOUNDARY STATE. IF CONV USE RESID TO AVOID ERRORS
c..   AT CONVERGENCE LEVEL: envelope error 0.001
c..   
      if( (dnab(kk)-dnad(kk))*(vnab(nvmax(3),3)-vnad(nvmax(3),3)) .lt.
     .     -1.0d-3 )then
         write(*,*)
     .        'ERR(cinit.f) Inconsistent boundary condition at join',kk
         write(*,'(2a5,8a12)')'kk','ib','h','ddnab','dvnab','criterion',
     .        'resid'
         write(*,'(2i5,1p8e12.4)')kk, ib(kk),h(kk),dnab(kk)-dnad(kk),
     .        vnab(nvmax(3),3)-vnad(nvmax(3),3),
     .        (dnab(kk)-dnad(kk))*(vnab(nvmax(3),3)-vnad(nvmax(3),3)),
     .        resid
c     stop'cinit'
      endif
c..   loop over envelope
      do j = 1, nvmax(3)
c..   loop over envelope, invert order for i
         i = nvmax(3)-j+1
c..   effective k index for envelope zones (i,k) matching
         k = nvmax(3)-i+kk
         if( vnab(i,3)-vnad(i,3) .gt. 1.0d-15 )then
c..   convective
            h(k) = vvel(i,3)
            ib(k) = 1
         else
            ib(k) = 0
         endif
      enddo


c..   --------------------------------------
c..   LOCATE RADIATIVE-CONVECTIVE BOUNDARIES
c..   --------------------------------------
      nrczones  = 0             !convective zones
      nrsczones = 0             !semi-convective zones
c..   
      kcbeg     = 0
      kcend     = 0
      
c..   NOTE:
c..   kbeg+1 = first convectively unstable boundary
c..   kend   = last  convectively unstable boundary
c..   Therefore loop (kbeg, kend) runs over nrczones CZs
      
      do j = 2, kk+nvmax(3)
         
c..   BOTTOM OF CONVECTION ZONE
         if( ib(j-1) .ne. 1 .and. ib(j) .eq. 1 )then
            nrczones = nrczones + 1
            kcbeg           = j-1
            kbeg(nrczones) = kcbeg
         endif
c..   TOP OF CONVECTION ZONE
         if( (ib(j) .eq. 1 .and. ib(j+1) .ne. 1) )then 
            kcend          = j+1
            kend(nrczones) = kcend
         endif
         
c..   BOTTOM OF SEMICONVECTIVE ZONE
         if(ib(j-1) .ne. 2 .and. ib(j). eq. 2) then
            nrsczones = nrsczones +1
            kscbeg(nrsczones) = j-1
         endif
         
c..   TOP OF SEMICONVECTIVE ZONE
         if(ib(j) .eq. 2 .and. ib(j+1) .ne. 2) then
            kscend(nrsczones) = j+1
         endif
         
c..   END LOOP OVER ZONES
      enddo
      
      if( nrczones .gt. 0 )then
c..   terminate unfinished outer convective zone:
         if( kend(nrczones) .le. 0 )then
c..   set kend to equivalent k outer zone (kk counted twice)
            kend(nrczones) = kk+nvmax(3)-1
            ib(kend(nrczones)) = 1
            kbeg(nrczones+1) = 0
            ib( kbeg(nrczones+1) ) = 0
         endif
      endif
      
c..   **WARN: THERE SHOULD BE NO SEMI-CONVECTIVE ZONES 
c..   TO THE PHOTOSPHERE! STRANGER THINGS HAVE 
c..   PROLLY HAPPENED THO
c..   -------------------------------------------------
c..   END CONVECTION ZONES LIMITS SEARCH
c..   -------------------------------------------------





c..   ---------------------------------------
c..   INTEGRAL PROPERTIES OF CONVECTION ZONES
c..   ---------------------------------------

      do k=1, kk
         ld(k) = 0.0d0
      enddo
      if( nrczones .gt. 0 )then
         
         do j = 1, nrczones

c..   RADIUS, MASS, PRESSURE AT CONVECTIVE BASE
c..   (deal with inner, outer boundaries if needed)
            if( kbeg(j) .eq. 1 )then
               rbot        = r(nc,1)
               xmbot      = xm(1)
               pbot       = p(nc,2)
               entrczb(1) = entropy(2)
            elseif( kbeg(j) .le. kk )then
               rbot  = r(nc,kbeg(j))
               xmbot = xm(kbeg(j))
               pbot  = p(nc,kbeg(j)+1)
               entrczb(j) = entropy(kbeg(j)+1)
            else
               i = nvmax(3)-(kbeg(j) - kk)
               rbot  = vr(i,3)
               xmbot = xm(kk+1) + vm(i,3)
               pbot  = vp(i,3)
               entrczb(j) = ventr(i,3)
            endif 
            
            
c..   RADIUS, MASS, PRESSURE AT OUTER CONVECTIVE BNDRY
            if( kend(j) .le. kk )then
               rtop  = r(nc,kend(j))
               xmtop = xm(kend(j))
               ptop  = p(nc,kend(j))
               entrczt(j) = entropy(kend(j))
            else
               i = nvmax(3)-(kend(j) - kk)+1
               rtop  = vr(i,3)
               xmtop = xm(kk+1) + vm(i,3)
               ptop  = vp(i,3)
               entrczt(j) = ventr(i,3)
            endif 
            if( j .eq. nrczones )then
c..               rsnzy = rtop                  !amanda:commented out 2/26/18, don't need anymore
            endif
            
c..   CONVECTION ZONES EXTENTS (M,R) AND PRESSURE CHANGE
            xmconv(j) = xmtop - xmbot
            rconv (j) = rtop  - rbot
            
            if( ptop .gt. 0.0d0 .and. pbot .gt. 0.0d0 )then
               pratio(j) = pbot/ptop
            else
               write(*,'(/a30)')'pratio error in CINIT.f'
               write(*,'(2(a10,1pe12.4))')'Pbot',pbot,'Ptop',ptop
               write(*,'(4(a10,i5))')'j',j,'kk',kk
               write(*,'(4(a10,i5))')'kbeg',kbeg(j),'kend',kend(j),
     .              'ib(kbeg)',ib(kbeg(j)),'ib(kend)',ib(kend(j))
               do k = kbeg(j)-5,kk+nvmax(3)
                  i = nvmax(3) - k + kk
                  write(*,'(3i5,1p8e12.3)')k,i,ib(k),
     .                 vnab(i,3)-vnad(i,3),vvel(i,3),vr(i,3),vtem(i,3),
     .                 vrho(i,3),vp(i,3),p(1,k),p(2,k)
               enddo
               write(*,*)kk+nvmax(3)
               stop'error p in cinit.f'
            endif 
            
c..   INTEGRATE OVER CONVECTION ZONE j
            sumq(j)  = 0
            sumke(j) = 0
            sum1(j)  = 0
            sum2(j)  = 0
            
            do k = kbeg(j)+1,kend(j)
               if( k .lt. kk )then
                  sumke(j) = sumke(j)+0.5d0*h(k)**2*dmi(k)/xmconv(j)
                  fak = -0.5d0*g(k)*t(nc,k)*pt(k)/(v(nc,k)*pv(k))
     .                 *(dnab(k)-dnad(k)-doux(k))*dmi(k)/xmconv(j)
c..   avoid semiconvective regions
                  if( fak .gt. 0.0d0 )then
                     sumq(j)  = sumq(j) + fak
                  endif
                  ld(k) = rconv(j)
                  
               else
                  i = nvmax(3) -(k-kk)
c..   kinetic energy of convective motion per unit mass
                  sumke(j) = sumke(j) + 
     .                 0.25d0*(vvel(i,3)**2+vvel(i-1,3)**2)*
     .                 (vm(i-1,3)-vm(i,3))/xmconv(j)
                  
c..   sumq is twice the specific KE in balance
c..   dmi is not centered on i
                  sumq(j) = sumq(j) + vcmp(i,3)*
     .                 grav*(xm(kk+1)+vm(i,3))/vr(i,3)**2*
     .                 (vnab(i,3)-vnad(i,3))
     .                 *(vm(i-1,3)-vm(i,3))/xmconv(j)*rconv(j)
                  
                  sum1(j) = sum1(j) + dmh(i)*vvel(i,3)**2/xmconv(j)
                  
c..   buoyancy frequency
                  nsqr(k) = -(grav*(xm(kk+1)+vm(i,3))/vr(i,3)**2)**2
     .                 *vrho(i,3)/vp(i,3)*(vnab(i,3)-vnad(i,3))*
     .                 vcmp(i,3)
                  ld(k) = rconv(j)
               endif
               if( kbeg(j) .ge. kk )then
                  i = nvmax(3) -(kbeg(j)-kk)
                  nsqr(kbeg(j)) = 
     .                 -(grav*(xm(kk+1)+vm(i,3))/vr(i,3)**2)**2 *
     .                 vrho(i,3)/vp(i,3)
     .                 *(vnab(i,3)-vnad(i,3))*vcmp(i,3)
               endif
            enddo
            
         enddo                  !END LOOP OVER CONVECTION ZONES
         ld(kk+1) = ld(kk)
      endif                     !END CZ LOGIC
c..   ---------------------------------------------------------





      
c..   RESET INNER BOUNDARY ZONE CONDITION
c..   NOW THAT COUNTING IS DONE
      ib(1) = ib(2)
      
      
c..   
c..   !!!! DEBUG:  ELIMINATING THIS CRITERIA
c..   
c..   DO NOT ALLOW FOR ISOLATED CONVECTION ZONES
      do i = kin,ktop-1
        if( ib(i) .eq. 0 )then
          if( ib(i-1) .eq. 1 .and. ib(i+1) .eq. 1 )then
            ib(i) = 1
          endif
        endif
      enddo

      do i = kin,ktop-1
        if( ib(i) .eq. 1 )then
          if( ib(i-1) .eq. 0 .and. ib(i+1) .eq. 0 )then
            ib(i) = 0
          endif
        endif
      enddo
      
      
c..   SAVE ENTRIES

c..   FIRST RUN, FIRST STEP
      if(l.le.1.and.nc.eq.1)then
         do i = kin,ktop+1
            ia(i) = ib(i)
            ip(i) = ib(i)
         enddo
      endif


c..   COUNT NUMBER OF ZONES WHICH CHANGED CONVECTIVE STATE
      ichange = 0
      do i = kin,ktop
         ik(i) = 0
         if( ip(i) .ne. ib(i) )then
            ichange = ichange + 1
            ik(i) = 1
         endif
      enddo
      
      
c..   COPY REVISED CONVECTIVE STATE ARRAY TO ic(k)
c..   ib -> ic
      do i = kin,ktop
         ic(i) = ib(i)
      enddo
      
      
c..   COUNT ZONES CHANGING CONVECT STATE THIS STEP
c..   (IF CALLED FROM HSTAT.f)
      if( it .ne. 0 )then
         kchange = 0
         do i = kin,ktop
            if( ic(i) .ne. ia(i) )then
               kchange = kchange + 1
            endif
         enddo
      endif
      

      if( modec .lt. 0 )then
c..   MODE: NO CONVECTION, ONLY RAD DIFF
         do i = kin,ktop
            ic(i) = 0
            h(i)  = 0.0d0
         enddo
      endif
      


c     SUCCESS
      return
c     
      end


