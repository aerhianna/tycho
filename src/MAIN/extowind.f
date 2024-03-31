
      subroutine extowind(te,tlum,xmdot,vinfin,zmet,rph,tmass)
c..   Mass loss routine for evolved stars
c     Uses the mass loss routine of Dupree, A. K. & Reimers, 
c     D. (1987, in Exploring the Universe with the IUE Satellite
c     eds. Y Kondo et al. (Dordrecht: Reidel) p. 321) for RGB 
c     stars (Teff < 7500 K) Modified for AGB stars (Period >
c     100d) to use mass loss formulation of Bloeker, T. (1995,
c     A&A 297, 727)
c---------------------------------------------------------------

      implicit none
      
      include 'cconst'

      real*8 te, tlum, xmdot, vinfin, zmet, rph, tmass
      real*8 q0, qmin, zmin, q1
      real*8 a0(3), a1(3), a2(3), ltilde
      real*8 a0q0in(3),a0q0out(3),a0qmin(3),a0qmout(3),
     1       a0zmin(3),a0zmout(3)
      real*8 a1q0in(3),a1q0out(3),a1qmin(3),a1qmout(3),
     1       a1zmin(3),a1zmout(3)
      real*8 a2q0in(3),a2q0out(3),a2qmin(3),a2qmout(3),
     1       a2zmin(3),a2zmout(3)
      real*8 tegrid(3), dumd, zbracket
      real*8 alpha, delta, f1, f2, f3, vesc

      data tegrid(1)/6.0d4/,tegrid(2)/5.0d4/,tegrid(3)/4.0d4/
      data a0q0in(1)/-5.99d0/,a0q0in(2)/-4.85d0/,a0q0in(3)/-5.2d0/
      data a1q0in(1)/1.0d0/,a1q0in(2)/0.5d0/,a1q0in(3)/0.93d0/
      data a2q0in(1)/1.5d0/,a2q0in(2)/1.0d0/,a2q0in(3)/0.85d0/
      data a0qmin(1)/-8.0d0/a0qmin(2)/-10.35d0/,a0qmin(3)/-11.75d0/
      data a1qmin(1)/-1.2d0/,a1qmin(2)/3.25d0/,a1qmin(3)/3.65d0/
      data a2qmin(1)/2.15d0/,a2qmin(2)/0.0d0/,a2qmin(3)/0.0d0/
      data a0zmin(1)/-3.4d0/,a0zmin(2)/-3.85d0/,a0zmin(3)/-4.45d0/
      data a1zmin(1)/-0.4d0/,a1zmin(2)/-0.05d0/,a1zmin(3)/0.35d0/
      data a2zmin(1)/-0.65d0/,a2zmin(2)/-0.6d0/,a2zmin(3)/-0.8d0/
c-------------------------------------------------------
      call spline(tegrid,a0q0in,3,a0q0out)
      call spline(tegrid,a1q0in,3,a1q0out)
      call spline(tegrid,a2q0in,3,a2q0out)

      call spline(tegrid,a0qmin,3,a0qmout)
      call spline(tegrid,a1qmin,3,a1qmout)
      call spline(tegrid,a2qmin,3,a2qmout)

      call spline(tegrid,a0zmin,3,a0zmout)
      call spline(tegrid,a1zmin,3,a1zmout)
      call spline(tegrid,a2zmin,3,a2zmout)

      ltilde = dlog10(tlum) - 6.0d0

      if( zmet .gt. 0.0d0 )then
         zbracket = dlog10(zmet)
      else
         write(*,*)'extowind error: zmet ',zmet
      endif

      call splint(tegrid,a0q0in,3,a0q0out,te,a0(1),dumd)
      call splint(tegrid,a1q0in,3,a1q0out,te,a1(1),dumd)
      call splint(tegrid,a2q0in,3,a2q0out,te,a2(1),dumd)

      q0 = a0(1) + a1(1)*ltilde + a2(1)*ltilde**2.0d0

      call splint(tegrid,a0qmin,3,a0qmout,te,a0(2),dumd)
      call splint(tegrid,a1qmin,3,a1qmout,te,a1(2),dumd)
      call splint(tegrid,a2qmin,3,a2qmout,te,a2(2),dumd)

      qmin = a0(2) + a1(2)*ltilde + a2(2)*ltilde**2.0d0

      call splint(tegrid,a0zmin,3,a0zmout,te,a0(3),dumd)
      call splint(tegrid,a1zmin,3,a1zmout,te,a1(3),dumd)
      call splint(tegrid,a2zmin,3,a2zmout,te,a2(3),dumd)

      zmin = a0(3) + a1(3)*ltilde + a2(3)*ltilde**2.0d0

      q1 = (q0-qmin)*(-zmin)**(-0.5d0)
c      write(*,*)'extO ',q0, qmin, zmin, zbracket, ltilde, tlum

      xmdot = 0.44d0*10.0d0**(q1*(zbracket-zmin)**0.5d0+qmin)

      vesc = (2.0d0*grav*tmass*sol/rph)**0.5

      if(zbracket .ge. -2)then
         alpha = 0.5d0
      elseif(zbracket .ge. -4)then
         alpha = (zbracket + 4.0d0)/(-2.0d0) * 0.1d0+ 0.3d0
cccccccccccccccccccccccccccccccccccccccccccc
      else
         alpha = (zbracket + 10.0d0)/(-6.0d0) * 0.3d0+ 0.0d0
cccccccccccccccccccccccccccccccccccccccccccc
      endif

      delta = 0.075d0

      f1 = 8.0d0/5.0d0*(1.0d0-0.75d0*alpha)
      f2 = dexp(-2.0d0*delta)
      f3 = 1.0d0 - 0.3d0*dexp(-vesc/3.0d7)

      vinfin = 2.25d0*alpha/(1.0d0-alpha)*vesc*f1*f2*f3

      write(*,*)"extreme O mass loss ",xmdot,vinfin

      return

      end
