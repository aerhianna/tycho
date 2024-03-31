      subroutine homol(fm,fr)

      implicit none

      include 'dimenfile'
      include 'cconst'
      include 'comod'

      real*8    fm,fr,fd,fv,fp,ft,fu,ftl,beta

      integer*4 k,kmax
c---------------------------------------------------------------
c..redefine model homologously
c  fm = desired mass/ mass of model
c  fr = radius desired/radius of model
c..some hints for scaling
c    It seldom if ever makes sense to scale an inhomogeneous model
c    The equation of state is assumed to be ideal gas plus radiation
c      so that degeneracy pressure is ignored in the scaling
c    Small changes, with fm/fr or sqrt(fm)/fr near unity (to keep
c      temperature changes small, work best
c---------------------------------------------------------------

c..scaling for hydrostatic equilibrium and mass conservation
c..thermal equilibrium must come from evolutionary settling
      fd   = fm / fr**3
      fv   = 1.0d0 / fd
      fp   = fd * fm / fr
      beta = 1.0d0 - arad*t(1,2)**4/(3.0d0*p(1,2))

c..add radiation pressure effect
      ft   = fp**( 1.0d0/(4.0d0-3.0d0*beta) ) / 
     1     fd**( beta / (4.0d0-3.0d0*beta)  )
      fu   = dsqrt(ft)
c..this is bad for low masses (electron degeneracy?)
      ftl     = ( ft * fr )**4 / fm

      write(*,'(8a12)')'fd','fp','ft','fu','ftl','fm**4','beta'
      write(*,'(1p8e12.3)')fd,fp,ft,fu,ftl,fm**4,beta

      kmax = kk+1
      do  k = 1, kmax
        t(2,k) = ft   * t(1,k)
        v(2,k) = fv   * v(1,k)
        p(2,k) = fp   * p(1,k)
        dmh(k) = fm   * dmh(k)
      enddo

      dmi(1)    = dmi(1)   *fm
      dmi(kmax) = dmi(kmax)*fm

      do k = 2, kk
        dmi(k)  =       0.5d0*( dmh(k) + dmh(k+1) )
        r(2,k)  = fr  * r(1,k)
        u(2,k)  = fu  * u(1,k)
        xm(k)   =       xm(k-1) + dmh(k)
        tl(2,k) = ftl * tl(1,k)
      enddo

      tl(2,kmax) = ftl * tl(1,kmax)
      xm(kmax)   = fm * xm(kmax)
      r(2,kmax)  = fr * r(1,kmax)
      r(2,kk+2)  = fr * r(1,kk+1)

c..   u(1,kmax) is unchanged
c..   scaled model has nc=2, eg, in t(2,k),...
c----------------------------------------------------------
c..end of redefine

      write(*,10)fm,fr
   10 format(1x,'scaled mass by',1pe14.6,', radii by',1pe14.6,
     1 ' in homol ')
      write(*,'(a20,1p8e12.3)')" scaled P, rho, T, L by ",
     1  fp, fd, ft, ftl

      return
      end
