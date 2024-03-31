      subroutine eos (kin,kk,tem,rho,ye,pe,pet,pev,ee,eet,eev)

c..   equation of state: fermi-dirac gas of electrons and positrons
c..   thermal part only
c..   rho = nucleons * amu

      implicit none

      include 'dimenfile'
      include 'ceos.h'
      include 'cconst'

      integer*4 kin, kk
      real*8    ff(4), f(3),df(3), g(4,3),dg(4,3)
      real*8    d(kdm),t(kdm)
      real*8    tem(kdm), rho(kdm),pe(kdm),pet(kdm),pev(kdm),
     1          ee(kdm),eet(kdm),eev(kdm),ye(kdm)

c..   local variables
      integer*4 k,ii,kd,kkt,maxd,mind,mint,maxt,n1,n2,n3,n4,jd
      real*8    del1,del2,del3,del4
      real*8    del12,del23,del24,del34,del41,del13
      real*8    b1,b2,b3,b4,c1,c2,c3,c4,alf1,alf2,alf3,alf4
      real*8    bet1,bet2,bet3,bet4,pee,ped,peet,erho,erhod,erhot,eta

      integer ikd(kdm),ikt(kdm), intab(kdm)
      common/ineos/ikd,ikt,intab
c---------------------------------------------------------------------- 
      if( kin .lt. 2 .or. kk .ge. kdm-1 )then
         write(*,*)'NEOS error: kin,kk ',kin,kk,kdm
         stop'eos'
      endif

      do k = kin, kk
         if( tem(k) .le. 0.0d0 )then
            write(*,*)'NEOS error: tem ',k,tem(k)
            stop'eos'
         endif
         if( rho(k) .le. 0.0d0 )then
            write(*,*)'NEOS error: rho ',k,rho(k)
            stop'eos'
         endif
         if( ye(k) .le. 0.0d0 )then
            write(*,*)'NEOS error: ye ',k,ye(k)
            stop'eos'
         endif
         d(k) = rho(k)*ye(k)
         t(k) = tem(k)
      enddo

c..   indices ikd,ikt and region flag intab are computed in state.f
      do ii = kin, kk
         if( intab(ii) .eq. 0 )then
c..   thermal part: from tables
c..   get indices
            kd  = ikd(ii)
            kkt = ikt(ii)
            if( kd .ge. ndd-1 )then
               maxd = ndd
               mind = maxd - 3
            elseif( kd .le. 2 )then
               mind = 1
               maxd = mind + 3
            else
               maxd = kd + 2
               mind = maxd - 3
            endif
c..   maxd = k + nd/2 and mind = maxd - nd for interp. order nd
c..   temperature
            if( kkt .ge. ntt-1 )then
               maxt = ntt
               mint = maxt - 3
            elseif( kkt .le. 2 )then
               mint = 1
               maxt = mint + 3
            else
               maxt = kkt + 2
               mint = maxt - 3
            endif
c..   maxt = k + nt/2 and mint = maxt - nt for interp. order nt
c..   interpolate in temperature
            n1 = mint
            n2 = mint + 1
            n3 = mint + 2
            n4 = maxt

            del1 = t(ii) - tt(n1)
            del2 = t(ii) - tt(n2)
            del3 = t(ii) - tt(n3)
            del4 = t(ii) - tt(n4)

            del12 = del1*del2
            del23 = del2*del3
            del24 = del2*del4
            del34 = del3*del4
            del41 = del4*del1
            del13 = del1*del3

            b1 = del2 * del34
            b2 = del1 * del34
            b3 = del12 * del4
            b4 = del12 * del3

            c1 = del23 + del24 + del34
            c2 = del13 + del41 + del34
            c3 = del12 + del41 + del24
            c4 = del12 + del13 + del23

            alf1 = b1/ttd(1,kkt)
            alf2 = b2/ttd(2,kkt)
            alf3 = b3/ttd(3,kkt)
            alf4 = b4/ttd(4,kkt)

            bet1 = c1/ttd(1,kkt)
            bet2 = c2/ttd(2,kkt)
            bet3 = c3/ttd(3,kkt)
            bet4 = c4/ttd(4,kkt)

c..   loop over densities
            do jd = mind, maxd
c..   set up variables (P=1,rho*E=2,mu/kt=3)
               ff(1) = array(jd,n1,1)
               ff(2) = array(jd,n2,1)
               ff(3) = array(jd,n3,1)
               ff(4) = array(jd,n4,1)
               f(1)  = alf1*ff(1) + alf2*ff(2)
     1               + alf3*ff(3) + alf4*ff(4)
               df(1) = bet1*ff(1) + bet2*ff(2)
     1               + bet3*ff(3) + bet4*ff(4)

               ff(1) = array(jd,n1,2)
               ff(2) = array(jd,n2,2)
               ff(3) = array(jd,n3,2)
               ff(4) = array(jd,n4,2)
               f(2)  = alf1*ff(1) + alf2*ff(2)
     1               + alf3*ff(3) + alf4*ff(4)
               df(2) = bet1*ff(1) + bet2*ff(2)
     1               + bet3*ff(3) + bet4*ff(4)

               ff(1) = array(jd,n1,3)
               ff(2) = array(jd,n2,3)
               ff(3) = array(jd,n3,3)
               ff(4) = array(jd,n4,3)
               f(3)  = alf1*ff(1) + alf2*ff(2)
     1               + alf3*ff(3) + alf4*ff(4)
c..   derivative of mue not needed?

c..   save values at density points for density interpolation
               g (jd-mind+1,1) = f (1)
               dg(jd-mind+1,1) = df(1)
               g (jd-mind+1,2) = f (2)
               dg(jd-mind+1,2) = df(2)
               g (jd-mind+1,3) = f (3)
c     dg(jd-mind+1,1) = df(1)
            enddo


c..   interpolate in density
            n1 = mind
            n2 = mind + 1
            n3 = mind + 2
            n4 = maxd

            del1 = d(ii) - dd(n1)
            del2 = d(ii) - dd(n2)
            del3 = d(ii) - dd(n3)
            del4 = d(ii) - dd(n4)

            del12 = del1*del2
            del23 = del2*del3
            del24 = del2*del4
            del34 = del3*del4
            del41 = del4*del1
            del13 = del1*del3

            b1 = del2 * del34
            b2 = del1 * del34
            b3 = del12 * del4
            b4 = del12 * del3

            c1 = del23 + del24 + del34
            c2 = del13 + del41 + del34
            c3 = del12 + del41 + del24
            c4 = del12 + del13 + del23

            alf1 = b1/ddd(1,kd)
            alf2 = b2/ddd(2,kd)
            alf3 = b3/ddd(3,kd)
            alf4 = b4/ddd(4,kd)

            bet1 = c1/ddd(1,kd)
            bet2 = c2/ddd(2,kd)
            bet3 = c3/ddd(3,kd)
            bet4 = c4/ddd(4,kd)

c..   set up variables (P=1,rho*E=2,mu/kt=3)
c..   Pressure
            ff(1) = g(1,1)
            ff(2) = g(2,1)
            ff(3) = g(3,1)
            ff(4) = g(4,1)
            f(1)  = alf1*ff(1) + alf2*ff(2) + alf3*ff(3) + alf4*ff(4)
            pee   = f(1)
c..   dP/d(rho)
            df(1) = bet1*ff(1) + bet2*ff(2) + bet3*ff(3) + bet4*ff(4)
            ped   = df(1)
c..   dP/dT
            ff(1) = dg(1,1)
            ff(2) = dg(2,1)
            ff(3) = dg(3,1)
            ff(4) = dg(4,1)
            f(1)  = alf1*ff(1) + alf2*ff(2) + alf3*ff(3) + alf4*ff(4)
            peet  = f(1)
c..   rho*E
            ff(1) = g(1,2)
            ff(2) = g(2,2)
            ff(3) = g(3,2)
            ff(4) = g(4,2)
            f(2)  = alf1*ff(1) + alf2*ff(2) + alf3*ff(3) + alf4*ff(4)
            erho  = f(2)
c..   d(rho*E)/d(rho)
            df(2) = bet1*ff(1) + bet2*ff(2) + bet3*ff(3) + bet4*ff(4)
            erhod = df(2)
c..   d(rho*E)/dT
            ff(1) = dg(1,2)
            ff(2) = dg(2,2)
            ff(3) = dg(3,2)
            ff(4) = dg(4,2)
            f(2)  = alf1*ff(1) + alf2*ff(2) + alf3*ff(3) + alf4*ff(4)
            erhot = f(2)
c..   mu/kT
            ff(1) = g(1,3)
            ff(2) = g(2,3)
            ff(3) = g(3,3)
            ff(4) = g(4,3)
            f(3)  = alf1*ff(1) + alf2*ff(2) + alf3*ff(3) + alf4*ff(4)
            eta   = f(3)
c..   derivative of mue not sought (df(3) not calculated)

            pe(ii)   = pee
            pet(ii)  = peet
            pev(ii)  = - rho(ii) * d(ii) * ped
c..   ee in erg/gm not erg/cc
            ee(ii)   = erho  / rho(ii)
            eet(ii)  = erhot / rho(ii)
            eev(ii)  = erho - d(ii) * erhod
         else
c..   out of table
            pe(ii)  = 0.0d0
            pet(ii) = 0.0d0
            pev(ii) = 0.0d0
c..   ete and eve are for E, elec here is rho*E as is e
            ee(ii)   = 0.0d0
            eet(ii)  = 0.0d0
            eev(ii)  = 0.0d0
         endif
      enddo

 1000 format(2i3,1p7e15.7)

      return
      end




