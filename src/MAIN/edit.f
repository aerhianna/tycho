      subroutine edit(gamma,del,tim2,no)
c..   7-30-06
      implicit none

      include 'dimenfile'

      include 'comod'
      include 'compu'
      include 'ctmstp'
      include 'cburn'
      include 'cruntm.h'
      include 'cphot'
      include 'cconst'
      include 'conline'
      include 'ceoset'
      include 'cenrchk'
      include 'cnabla'

      real*8 d(kdm)
      real*8 tim2, xmbol,dtdv,gam,ropal
      real*8 dechk, ddechk, bind, etot, eergs, ekerg
      real*8 gamma, del, dmas, xsum
      real*8 dlepst, dlepsv
      integer*4 iwrite, no, knd, n, k, nii

c----------------------------------------------------------------

c..   write flag enforced here rather than in main
      if( mod(l,l2) .eq. 0 .or. l .ge. ll )then
c..   write option
         iwrite = 1
      else
c..   nowrite option
         iwrite = 0
      endif

      if( iwrite .eq. 1 )then
         write(3,14) model,mode,l,it,no,time,dth(2),etar
         write(3,42) modec, modes, modex, tenvelop
         write(3,16) runt,tim2
      endif

      xmbol = -2.5d0*xlol + 4.72d0
      if( iwrite .eq. 1 )then
         write(3,15)xlol,telog,xmbol,rphot,uphot
      endif

c...............................................................

      if( iwrite .eq. 1 )then
         write(3,*)'time step criterion, zone, value(sec):'
         write(3,'(4(2x,a5,i4,1pe9.2))')(dthflag(knd+1),
     1        kkmin(knd),tmin(knd), knd = 1, ncond)

         if( float(it) .gt. 0.2*float(iter) )then
            write(3,*)'last zone to converge: K=> ',kitmax,
     1           ' variable=> ',cnit(nitest)
            write(3,'(8a12)')'T','rhoz','rhoz*Ye','xm','R'
            write(3,'(1p8e12.3)')tem(kitmax),
     1           rhoz(kitmax),x(ndim,kitmax)/v(2,kitmax),xm(kitmax),
     2           r(2,kitmax)
            write(3,'(8a12)') 'P','PT','PV','E','ET','EV'
            write(3,'(1p8e12.3)')p(2,kitmax),pt(kitmax),
     1           pv(kitmax),e(2,kitmax),et(kitmax),ev(kitmax)
            write(3,'(8a12)')'dlnP/dlnT','dlnP/dlnV','dlnE/dlnT',
     1           'dlnE/dlnV'
            write(3,'(1p8e12.3)')
     1           pt(kitmax)*t(2,kitmax)/p(2,kitmax),
     2           pv(kitmax)*v(2,kitmax)/p(2,kitmax),
     3           et(kitmax)*t(2,kitmax)/e(2,kitmax),
     4           ev(kitmax)*v(2,kitmax)/e(2,kitmax)
            write(3,'(8a12)')'ak','akt','akv','dlnak/dlnT','dlnak/dlnV' 
            write(3,'(1p8e12.3)')ak(kitmax),
     1           akt(kitmax),akv(kitmax),
     2           akt(kitmax)*t(2,kitmax)/ak(kitmax),
     3           akv(kitmax)*v(2,kitmax)/ak(kitmax)
            dtdv = - (p(2,kitmax) + ev(kitmax))/et(kitmax)
            gam  =-v(2,kitmax)/p(2,kitmax)*(pv(kitmax)+pt(kitmax)*dtdv) 
            ropal = rhoz(kitmax)/(tem(kitmax)*1.0d-6)**3
            write(3,'(8a12)')'Ye','Yef','dlnYef/dlnT','dlnYef/dlnV',
     1           'dlnT/dlnV','Ye-Yef'
            write(3,'(1p8e12.3)')ye(kitmax),
     1           yef(kitmax),yeft(kitmax)*t(2,kitmax)/yef(kitmax),
     2           yefv(kitmax)*v(2,kitmax)/yef(kitmax),
     3           dtdv*v(2,kitmax)/t(2,kitmax),ye(kitmax)-yef(kitmax)

            write(3,'(8a12)')'gam','log ropal'
            write(3,'(1p8e12.3)')gam,dlog10( ropal )

            write(3,'(8a12)')
     1 'eps(nu)','eps(nuc)','deps/dT','deps/dV','dlneps/dT','dlneps/dV'
            if( s(5,kitmax)+s(4,kitmax) .ne. 0.0d0 )then
               dlepst = st(kitmax)*t(2,kitmax)/(s(5,kitmax)+s(4,kitmax))
               dlepsv = sv(kitmax)*v(2,kitmax)/(s(5,kitmax)+s(4,kitmax))
            else
               dlepst = 0.0d0
               dlepsv = 0.0d0
            endif
            write(3,'(1p8e12.3)')
     1 s(4,kitmax),s(5,kitmax),st(kitmax),sv(kitmax),dlepst,dlepsv
         endif
c         write(3,'(a5,10a12)')'K','xm','r','tl','T','rho','h','ss',
c     1   'n-a-x'
c         do k = 100, 120
c            write(3,'(i5,1p10e12.4)')k,xm(k),r(2,k),tl(2,k),t(2,k),
c     1	    1.0d0/v(2,k),h(k),ss(k),dnab(k)-dnad(k)-doux(k)
c         enddo
      endif

c     radiation (neutrino) stress not defined
c     echeck = echeck - stress

      bind   = -eint + pote
      etot   =  ekin + eint - pote
      dechk  = echeck - echkz
c..define fraction to include integrated energy release, loss
      if( echeck .ne. 0.0d0 )then
         ddechk =  dechk/(abs(etot) + abs(elum) + abs(esou) )
      else
         ddechk = 0
      endif



      if( iwrite .eq. 1 )then
         write(3,*)'Energy conservation checks (erg/gram):'
         write(3,'(6a13)')'frac. err','error','echeck','kinetic',
     1        'internal', 'potential'
         write(3,'(1p6e13.5)') ddechk,dechk,echeck,ekin,eint,
     1        pote
         write(3,'(6a13)')'Sources','Radiation','HD damping', 
     1        'GR M def.','Binding'
         write(3,'(1p6e13.5)')esou,elum,diske, dmgc2,bind

         eergs = etot*(xm(kk)-xm(1))
         ekerg = ekin*(xm(kk)-xm(1))
         write(3,'(a15,1pe23.15,a15,1pe23.15,a5)')'E(total) = ',eergs,
     1        ' ergs, E(kin) = ',ekerg,' ergs' 
      endif

c.....calculates gamma

      call stabil(2,gamma,del)

      if( iwrite .eq. 1 )then
         write(3,40)gamma,del
         if( del .le. 0.0 )  write(3,41)
      endif

c.....sums sources and sinks

      if( iwrite .eq. 1 )then
         call sosum(kk,s,dmh,xm,dth,e,u)
      endif

c.....mass integrated compositions

      dmas = xm(kk) - xm(1)
      do n = 1, ndim
         d(n) = 0.0d0
         do k = 2,kk
            d(n) = d(n) + dmh(k)*x(n,k) / dmas
         enddo
      enddo

      if( iwrite .eq. 1 )then
c..   abundances summed over zones; width is 6*(5+4+9)=108 columns
         write(3,*)' A and total mass fraction X'
         do nii = 1, netsize, 6
            if( nii+5 .lt. netsize )then
               write(3,'(6(a5,1pe9.2))')(cnuc(n),d(n),n=nii,nii+5)
            else
               write(3,'(6(a5,1pe9.2))')(cnuc(n),d(n),n=nii,netsize)
            endif
         enddo
      endif


c.....check of nucleon conservation
      xsum = 0
      do  n = 2, ndim
         xsum = xsum + d(n)
      enddo
      if( iwrite .eq. 1 )then
         write(3,10)xsum
      endif



 10   format(' nucleon check =',f14.10,//)
 11   format(1x,  1p8e13.5)
 12   format(1x,7(i3,i4,1pe9.1))
 14   format(1h //' EDIT: summary of time step..................'/,
     1     ' model=',i7,' mode=',i2,' l=',i7,' it=',i3,' no=',
     2     i2, 2x,'time', 1pe15.7, 3x,'dth(2)',
     3     1pe11.3,3x, 'etar',1pe11.3)
 15   format(/1x,'log(l/sol) =',f6.3,3x,' log(te(K)) =',f6.3,
     1     3x,' mbol =',f7.3, 3x, 'rphot =', 1pe11.3,
     2     3x,'uphot = ', 1pe11.3)
 16   format(' elapsed  cpu time is',f10.3, 1x,
     1     'seconds, elapsed star time is',1pe15.7,' seconds')
 18   format(1x,8a13)


 40   format(' gamma =',  f15.7,' del =',  f15.7)
 41   format(' !!!!! hydrodynamic instability !!!!!')
 42   format( ' modeC =', i2, '  modeS =', i2,
     1     '  modeX =', i2, '  Tenvelop =', 1pe11.3)


      return
      end







