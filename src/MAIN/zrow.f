      subroutine zrow

c..initializes arrays
c-----------------------------------------------------------------------
      implicit none

      include 'dimenfile'

      include 'compu'
      include 'comod'
      include 'ctmstp'
      include 'cqdv.h'
      include 'cruntm.h'
      include 'cbug'
      include 'czone'
      include 'cnabla'
      include 'cphot'
      include 'cgtr.h'
      include 'cconst'

      integer*4 k, n
c-----------------------------------------------------------------

        do k = 1,kdm
          do n = 1, ndim
            x(n,k)    = 0.0d0
            xd(n,k)   = 0.0d0
            xold(n,k) = 0.0d0
            xdcon(n,k) = 0.0d0
          enddo

          do n = 1,7
            s(n,k) = 0.0d0
          enddo
c..nuclear burning energy generation rate and ln derivatives (T,V)
          ss(k)   = 0.0d0
          sa(k)   = 0.0d0
          sb(k)   = 0.0d0
c..neutrino energy generation rate (negative) and ln derivatives (T,V)
          snu(k)  = 0.0d0
          snua(k) = 0.0d0
          snub(k) = 0.0d0

          epsnuc(k) = 0.0d0
cccccccccccccccccccccccccccc
        enddo

      do k = 1,kdm
        do  n = 1,2
          tl(n,k) = 0.0d0
          r(n,k)  = 0.0d0
          u(n,k)  = 0.0d0
          e(n,k)  = 0.0d0
          p(n,k)  = 0.0d0
          t(n,k)  = 0.0d0
          v(n,k)  = 0.0d0
          q(n,k)  = 0.0d0
          dv(n,k) = 0.0d0
          dt(n,k) = 0.0d0
        enddo
      enddo

      do k = 1,kdm
      a(k)   = 0.0d0
      f(k)   = 0.0d0
      g(k)   = 0.0d0
      b(k)   = 0.0d0
      du(k)  = 0.0d0
      duold(k) = 0.0d0
      sound(k) = 0.0d0
      dtl(k) = 0.0d0
      h(k)   = 0.0d0
      y(k)   = 0.0d0
      z(k)   = 0.0d0
      xm(k)  = 0.0d0
      td(k)  = 0.0d0
      ak(k)  = 0.0d0
      et(k)  = 0.0d0
      ev(k)  = 0.0d0
      pt(k)  = 0.0d0
      pv(k)  = 0.0d0
      dr(k)  = 0.0d0
      dmi(k) = 0.0d0
      akt(k) = 0.0d0
      akv(k) = 0.0d0
      st(k)  = 0.0d0
      sv(k)  = 0.0d0
      dmh(k) = 0.0d0
      ssum(k) = 0.0d0
      diss(k) = 0.0d0
      dnab(k) = 0.0d0
      dnad(k) = 0.0d0
      dnrad(k) = 0.0d0
      ic(k)   = 0
      enddo

      do k = 1,kdm
        omeg(k)  = 0.0d0
        ajay(k)  = 0.0d0
      enddo

      do n = 1, 2
      dth(n) = 0.0d0
      dti(n) = 0.0d0
      enddo


      do n = 1, 4
        wtest(n)  = 0.0d0
        nwtest(n) = 0
      enddo


      return
      end

