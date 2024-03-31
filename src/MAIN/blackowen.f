      subroutine blackowen(doprint)

      implicit none

      include 'dimenfile'
      include 'comod'
      include 'cconst'
      include 'cnabla'
      include 'cadvec'
      include 'ceoset'

c     Integers:
      integer*4 doprint,i                                     !input and counting variables
      integer*4 kcbeg,kcend,kmax,kbegnrc,kendnrc              !various indices

c     Constants:
      real*8    sunage,secpd,m_H                              !various constants
      real*8    shearsol,rossbysol                            !Sun shear & Rossby values
      real*8    tcxsun,Tnott_sun,Fx_sun,Lx_sun,L_sun          !Sun temp/flux/luminosity values
      real*8    Mdot_sun,Omega_sol,Br_sol                     !Sun mass loss, rotation, B-field values
      real*8    qalpha,thetas,q_inertia,om4pi,uA_sun          !various values used in code
      real*8    lambda                                        !lambda values

c     Star internal stuff:
      real*8    rkbeg,rkend,rkk,rkk1,rkk2,xmkk,xmkk1,xmkk2    !radius, mass
      real*8    bo_omega,bo_omega2                            !used to print omega
      real*8    vconv                                         !used to print conv. velocity (h)

c     Star computed stuff:
      real*8    psh(kdm)                                      !pressure scale height
      real*8    rsnorm,xshear,rossby,alphanot                 !various computed values
	  real*8    tauc,taug,taued                               !tau values
      real*8    gL,newgL                                      !gL
      real*8    br,newbr                                      !br
      real*8    lx,lratio                                     !lx & lratio
      real*8    fx                                            !fx
      real*8    fc,lc,fw,lw                                   !fc/lc & fw/lw
	  real*8    bigLx                                         !bigLx
      real*8    pnot,pnot17                                   !pnot
      real*8    Br_star,Br_star_check                         !Br_star
      real*8    omegadot,omegadot_check                       !omegadot
      real*8    uA,uA_check                                   !uA
      real*8    anot,Mdot_star                                !other various computed values
      real*8    LxLconv,lx25,Lalpha,L_conv(kdm),Lconv         !L stuff for B.O. eq25
      real*8    rhoconv,havg,hcount                           !more BO eq25 stuff

c     Tnot stuff:
      real*8    Tnot,Tnott                                    !Tnot & Tnott
      real*8    Tnot_a,Tnot_b,Tnot_c                          !Tnot_a/b/c
      real*8    Tnot_in                                       !Tnot_in
      real*8    Tnot_out                                      !Tnot_out

c     Reiners stuff:
      real*8    Omega_crit,ReinersC,dJdt,mdot,rdot
      real*8    deltaJ,delta_Omega,omegadot2

c     kk = outer index of interior zones, from comod
c     kk+1 = either outer index of exterior zones, or edge of star?
c     kk+2 = similar to kk+1? maybe the edge of star?
c     kdm = max zone size, in dimenfile
c     kbeg(nrczones) = index of beginning of convective zone?
c     kend(nrczones) = index of end of convective zone? (acts weird?)

c     gravitational constant G = grav from cconst [cm^3/g/s^2]
c     mass = xm(i) from comod = mass within radius i? [g]
c     radius = r(nc,i) from comod = radius at i at start (nc=1) or end (nc=2)? = cm?

c     sound speed = sound(i) from comod, state.f? = speed at radius i? [cm/s]
c     NOTE: sound = adiabatic sound speed.
c           want isothermal? see timmes.f for csound, but not sure what it is exactly?

c     shear = shear(i) from cnabla = shear at radius i? [unitless?]
c     angular velocity = omeg(i) from comod [rad/s]
c     rotational period = 2pi/angular velocity [s]
c     (also have rotational velocity vrot, in comod)

c     sol = solar mass in grams, from cconst
c     solrad = solar radius in cm, from cconst
c     secpy = seconds per year

c     rsnzy = radius of outer boundary of outer convective zone = rtop in cinit.f (commented out)
c     NOTE: tmloss from cadvec doesn't seem to be anywhere else in .f files, only .o binary files?
c     bomloss = yes/no flag set in params.d

c     start:
      if (doprint .eq. 1) then
         print *,' '
         print *,
     1      'START BLACKOWEN BLACKOWEN BLACKOWEN BLACKOWEN BLACKOWEN'
         print *,' '
      endif

      if (log(tl(1,kk+1)) .le. -0.15) then
         if (log(t(1,kk+1)) .ge. 3.7446) then
            if (zamstime .eq. 0.0) then
               zamstime = time
               zamsomega = omeg(kk)
            endif
         endif
      endif

      if (time/secpy .le. 5.d7) then
	     bosurfomeg = omeg(kk+1)
      endif

c     initialize variables:
      tcxsun = 1.5d6              !sun's temperature at base of the corona
      qalpha = 0.5                !(arbitrarily chosen)   0 < qalpha < 1
      thetas = 0.d0               !fiducial polar angle   ?does this apply if 1D?
      lambda = 0.16               !(0 < lambda < 1/3) -> 0.16  [Fig. 5 in B.O. they use 1/3]

c      print *,'lambda (norm/min/max)....= ',lambda,lambdamin,lambdamax

      sunage = 4.57d9 * secpy     !solar age in seconds, from bonanno et al 2002 arxiv:0204331
      shearsol = 8.3d0            !section 6.1 of the B+O paper, it discusses using shear = 8.3
      rossbysol = 2.d0            !rossby number for sun
      om4pi = 0.1                 !Big Theta = omega/4pi < 1, conduction down from corona is non neg.
                                  !along field lines perp. to surface (after Eq.22, defined as 1/10)
      q_inertia = 0.8             !inertial parameter q >= 1, num. factor, reduce moment of inertia
                                  !by an amount that depends on internal angular momentum transport
                                  !In B.O. Figure 5 they use q=1.2
      Br_sol = 2.                 !solar value from B+O paper: surface radial field, G(auss)
      Omega_sol = 2.97d-6         !solar value: seconds^-2
      Tnott_sun = tcxsun/3.d6     !Tnott_sun = 0.5, normalized temp (Tnot~_sun = Tcxsun/3d6)
      L_sun = 4.d33               !solar luminosity, erg/sec [from B.O.]
      Lx_sun = 6.d-7 * L_sun      !solar x-ray luminosity value [from B.O.]
      Fx_sun = Lx_sun/4./PI/solrad**2.
      m_H = 1.67d-24              !mass of proton (grams)
      Mdot_sun = 1.3d12           !Mass loss rate for sun (gram/sec)
      secpd = secpy/365.25        !seconds per day
      uA_sun = 2.64d7             ! cm/s

c     Tnot = temperature of corona [K]
c     Tnott = Tnot~ = T0_tilda = temperature of corona [3.d6 K] = Tnot/3.d6
c     Regime I: T0_tilda < T0sol = 0.5, solar type stars older than the sun
c     Regime II: T0_tilda > T0sol = 0.5, solar type stars younger than sun

c     print out some details about the star and indexes:
      kbegnrc = kbeg(nrczones)  !index of beginning of convective zone
      if(kbegnrc .lt. 3)kbegnrc = 3
      write(*,*)'kbegnrc = ',kbegnrc
      kendnrc = kend(nrczones)     !index of end of conv. zone (weird?)
      rkk = r(1,kk)/solrad         !radius at kk [outer index of interior zones]
      rkk1 = r(1,kk+1)/solrad      !radius at kk+1 [outer radius of star? or conv. zones?]
      rkk2 = r(1,kk+2)/solrad      !radius at kk+2 [not sure?]
c      print *,'r(kk)  [R_sun]...........= ',rkk
      print *,'r(kk1) [R_sun]...........=  ',rkk1
c      print *,'r(kk2) [R_sun]...........= ',rkk2
      rkbeg = r(1,kbegnrc)/solrad  !radius at beginning of conv. zone
      rkend = r(1,kendnrc)/solrad  !radius at end of conv. zone
c      print *,'r(kbegnrc) [R_sun].......= ',rkbeg
c      print *,'r(kendnrc) [R_sun].......= ',rkend
      xmkk = xm(kk)/sol            !mass within interior zones
      xmkk1 = xm(kk+1)/sol         !mass at kk+1 [mass of star?]
      xmkk2 = xm(kk+2)/sol         !mass at kk+2
c      print *,'xm(kk)  [M_sun]..........= ',xmkk
      print *,'xm(kk1) [M_sun]..........=  ',xmkk1
c      print *,'xm(kk2) [M_sun]..........= ',xmkk2
      print *,'time (age of * in Gyr)...= ',time/secpy/1.0d9
      print *,'Omega(kk)/Omega_sol .....= ',omeg(kk)/Omega_sol
      print *,'Omega(kk1)/Omega_sol ....= ',omeg(kk+1)/Omega_sol
      Lalpha = 2./5.*r(1,kbegnrc)  !L_alpha from below eq.25 in B.O.

c	  equation 4 from Blackman and Owen paper, to find the sonic radius:
c     (this value isn't used anywhere else)
      rsnorm = grav * xm(kk+1) / 2.d0 / sound(kk+1)**2./ r(1,kk+1)

c     equation for pressure scale height, H_p, taken from line 108 of getvec.f
      if( modes .eq. 2 )then
         kmax = kk+1
      else
         kmax = kk
      endif
      psh(1) = 0.
      do i = 3, kmax-1
         psh(i-1) = 0.5d0 * (p(1,i)+p(1,i+1)) / (p(1,i)-p(1,i+1))
     1           * 0.5d0 * (r(1,i+1)-r(1,i-1))
      enddo
      do i = kmax, kmax+1
         psh(i-1) = psh(i-2)
      enddo
c     print *,'H_p(kbegnrc) [cm?].......= ',psh(kbegnrc)

c     equation for convective turnover time, tau_c, taken from Landin et al. 2010
c                 (from the text before eq. 4)                  arxiv 1001.2754
c	  (tau_c = alpha * H_p / v)
c     alpha = mixing length parameter in pressure scale heights = alphaml in comod
      
c     (alphaml = 1.6?)   (Landin et al. use alpha=1.5 I think?)
c     v = convective velocity = h in comod   (note: hp = convective velocity at end of timestep)
c     (seems to work when psh taken at beginning of convective zone, but h taken at edge of star?)

      tauc = alphaml * psh(kbegnrc)/h(kk+1)
c     print *,'h(kk1) [cm/s?]...........= ',h(kk+1)
      print *,'tauc [day] ..............= ',tauc/secpd

c     equation for global convective turnover time, tau_g, taken from Landin eq. 4
c     (tau_g = integral(r_0 to r_star) 1/v(r) dr) (r_0 = radius at beggining of convective zone)
c     (v = convective velocity)                   (r_star = radius of star)
	  taug = 0.
	  kcbeg = kbegnrc
      kcend = kk+1            !kend, or kk?
	  do i=kcbeg+1,kcend
	     taug = taug + (r(1,i)-r(1,i-1))/h(i)
	  enddo
      print *,'taug [day]...............= ',taug/secpd

c     using tau_g to tau_c conversion factor computed from Landin et al.
c     print *,'tauc/taug ratio should be near 0.41 [Landin et al.]'
c     print *,'tauc/taug ...............= ',tauc/taug

c     ROSSBY NUMBER (eqn just below eq.1 in B.O., or fig.8 in Landin)
      rossby = 2.d0*pi/omeg(kk+1)/tauc
      rossby = 2.d0*pi/omeg(kk)/tauc
      rossby = 2.d0*pi/bosurfomeg/tauc
      print *,'rossby...................= ',rossby

c     equation for shear from just below eq.9 in B.O.
      xshear = omeg(kk+1)/ABS(omeg(kk+1)-omeg(kbegnrc))
      xshear = omeg(kk)/ABS(omeg(kk)-omeg(kbegnrc))
      xshear = 8.3
      print *,'B.O. use fixed shear of 8.3 in Section 6.1'
      print *,'xshear...................= ',xshear

c	  equation 9 from Blackman and Owen paper for Tau_ed (turbulent correlation time):
c     (this values is only used for alphanot equation below)
      taued = xshear * 2.d0*PI/omeg(kk+1) /
     1        (1.d0 + xshear * rossby)

c	  equation 8 from Blackman and Owen paper for kinetic helicity (magnetic fluctuations):
c     (this value isn't used anywhere else)
      alphanot = qalpha/6.d0 * taued**2. * omeg(kk+1) * h(kk+1)**2.
     1           / r(1,kbegnrc) * COS(thetas)

c	  gL, From below equation 26 in B.O., using different values for lambda:
      gL = (1.4 - (0.4 * time/sunage))**((lambda-1.)/4.)
      print *,'gL.......................=  ',gL, time, lambda, sunage

c	  equation 14 in B.O.: normalized surface radial magnetic field magnitude
	  br = gL * (xshear/shearsol)**(1./6.)
     1     * (SQRT(1. + (shearsol * rossbysol))
     2     / (SQRT(1. + (xshear * rossby))))
c     br = gL * (solrad/r(1,kk+1))**2.
      print *,'br [radial mag field]....= ',br

c	  equation 26 in B.O. gives the x-ray luminosity as function of magnetic field strength:
	  lx = br**(4./(1. - lambda))
c     print *,'lx ......................= ',lx

c     computing convective luminosity, taken from line 468 in getvec.f
c..   implied convective luminosity = L - L(radiative)
c     (seems to be in units of L_sun)
      L_conv(kmax) = 0
      L_conv(kmax-1) = 0
      do i = 2, kmax-1
         L_conv(i-1) = ( tl(1,i)
     1      -(pi4*r(1,i)**2)**2/dmi(i)*cflux
     2      *(t(1,i)**4-t(1,i+1)**4)/( 0.5d0*( ak(i) + ak(i+1)) ) )
     3      /sollum
      enddo
      Lconv = 0.
      do i=1,kmax
         if (L_conv(i) .gt. Lconv) Lconv = L_conv(i)
      enddo

c     equation 25 in B.O., L_x/L_conv:
      LxLconv = (Lalpha/r(1,kbegnrc))**(2./3.)
     1        * (xshear**(1./3.)/(1.+xshear*rossby))**2.
     2        * om4pi*4.*PI*(3.*PI/8.)**(1./3.)
     3        * (qalpha*COS(thetas)/6.)**(2./3.)
      lx25 = LxLconv * Lconv * L_sun / Lx_sun
c     print *,'lx25 [Tycho].............= ',lx25
      rhoconv = (xm(kk+1)-xm(kbegnrc))/
     1       ((4./3.*PI*r(1,kk+1)**3.)-(4./3.*PI*r(1,kbegnrc)**3.))
      havg = 0.
      hcount = 0.
      do i=kbegnrc,kk+1
         havg = havg + h(i)
         hcount = hcount + 1.
      enddo
      havg = havg/hcount
      Lconv = 4.*PI*r(1,kbegnrc)**2.*
     1        rhom(kbegnrc)*h(kbegnrc)**3.
      lx25 = LxLconv * Lconv / Lx_sun
c     print *,'lx25 [B.O.]..............= ',lx25

c     print out implied mass loss from initial lx values (see B.O. eq21):
c     print *,'Mdot [M_sun/yr]..........= ',
c    1      -lx * Mdot_sun/sol*secpy

c     opening files to output any test stuff
      open(54,file='test_output_P.txt')
      open(56,file='test_output_L.txt')
      open(57,file='test_output_PTdata.txt')

c     write output test file headers:
      write(56,*) 'star mass,radius in solar units ..= ',
     1       xm(kk+1)/sol,r(1,kk+1)/solrad
      write(56,*) 'star age in Gyr ..................= ',
     1       time/secpy/1.d9
      write(56,*)
     1       'i -- lx -- bigLx/Lx_sun -- lratio -- pnot -- Tnott'

c     while loop starting with lx from eq.26
c     computing T from eq.21, computing p from eq.20, computing fx from eq.18
c     then inputting resulting lx back into start of loop
c     until the fractional difference between lx's is < 0.001
      i = 0              ! keeps track of number of loops
      lratio = 100.      ! starting lratio, just needs to be high
c     print *,' '
c     print *,'Starting lx convergence loop!'
      do while (lratio .gt. 0.001)
c     print *,'-------------------------------------------'

c     Tnot~ calculated by inverting eq. 21 in B.O.
c     (NOTE: B.O. equation 22 approximation does NOT produce a correct value!)
      Tnot_a = LOG(lx)
      Tnot_b = (7.8 * xm(kk+1)/sol * solrad/r(1,kk+1)) / Tnott_sun
      Tnot_c = -(7.8 * xm(kk+1)/sol * solrad/r(1,kk+1))
      Tnot_in = -Tnot_c * exp(-(Tnot_a - Tnot_b))
      Tnot_in = 1.0
      call prodlog(Tnot_in,Tnot_out)         !used Mathematica to solve eq.21 for Tnott
      Tnott = -Tnot_c / Tnot_out
      Tnot = Tnott * 3.d6
c     print *,'Tnott ...................= ',Tnott

c	  equation 16 in B.O.: energy per unit time advected away by mass loss within
c     coronal scale height
c     Mdot_star = Mdot_sun * lx (because lx = Lx_star/Lx_sun = Mdot_star/Mdot_sun?)
c     (this goes into pnot17 below, but that's it)
      fw = ((Mdot_sun * lx) * Tnott*3.d6 * ((3 * boltz)/ m_H))
     1              / (4 * pi * (r(1,kk+1))**2.)

c     equation 17 in B.O. used to solve for pnot
c     (decided this equation was inaccurate, so it isn't used anywhere else)
      pnot17 = fw / (3.1d6 * Tnott**(0.5)
     1         * exp(3.9 * xm(kk+1)/sol * solrad/r(1,kk+1)
     2         * (1. - Tnott**(-1.))))
c     print *,'pnot_eq17................= ',pnot17

c     check pnot *pressure* value of Eq.20 in B.O. using our calculated Tnott (Tnot~) value:
      pnot = (xm(kk+1)/sol * (solrad/r(1,kk+1))**2.
     1      * 1.6 * om4pi*4.*PI * Tnott**(29./12.)) +
     2      (xm(kk+1)/sol * (solrad/r(1,kk+1))**2. * 0.75
     3      * Tnott**(13./6.) * exp(3.9 * xm(kk+1)/sol *
     4      solrad/r(1,kk+1) * (1.-(1./Tnott)))) +
     5      (xm(kk+1)/sol * (solrad/r(1,kk+1))**2. * 2.34
     6      * Tnott**(7./6.) * exp(3.9 * xm(kk+1)/sol *
     7      solrad/r(1,kk+1) * (1.-(1./Tnott))))

c     print *,'pnot [g/(cm*s^2)]........= ',pnot

c     equation 23: relationship between total integrated wind power and wind flux
c     (this value isn't used anywhere else [I don't think?])
      lw = 4. * pi * (r(1,kk+1))**2. * fw * 5.1 * xm(kk+1)/sol *
     1       1./Tnott * solrad/r(1,kk+1)

c     equation 18: x-ray radiation flux
c     Fx should increase with increasing Tnot (and pressure):

      fx = 1.24d6 * (pnot)**2.*(1./Tnott)**(5./3.)
     1         * (r(1,kk+1)/solrad)**2. * sol/xm(kk+1)
c     print *,'fx/fx_sun ...............= ',fx/Fx_sun

c     Lx == 4*pi*R^2 * Fx
      bigLx = 4. * PI * (r(1,kk+1))**2. * fx
c     lx == Lx/Lx_sun
c     print *,'bigLx/Lx_sun ............= ',bigLx/Lx_sun

      lratio = ABS(bigLx/Lx_sun - lx)/lx

      write(56,*) i,lx,bigLx/Lx_sun,lratio,pnot,Tnott

c     send new lx value back into start of convergence loop:
c      lx = bigLx/Lx_sun
      lratio = 0.0001

      i = i + 1
      enddo

      close(54)
      close(56)
      close(57)
      print *,'-------------------------------------------'
      print *,'Final values:'
c     print *,'i .......................= ',i
      print *,'Tnott > 0.5 --> younger than Sun'
      print *,'Tnott (0.5 for Sun)......= ',Tnott
      print *,'pnot [g/(cm*s^2)?].......= ',pnot
      print *,'Fx/Fx_sun ...............= ',fx/Fx_sun
      print *,'lx ......................= ',lx

c     print *,'Tnott from eq4 in J.G. ..= ',0.11 * fx**(0.26) / 3.

c     ****** ORGANIZED/CLEANED CODE UP TO HERE ******

c     equation 19/23: conductive loss flux/luminosity
      fc = 4.26d6 * pnot * (Tnott)**(3./4.) * om4pi
      print *,'Fc [erg/cm^2*s]..........= ',fc

      lc = 4. * pi * (r(1,kk+1))**2. * fc
      print *,'Lc [erg/s]...............= ',lc

c     calculating sound speed for base of corona (NOT TAKING FROM TYCHO):

      anot = ((2. * boltz * Tnot) / m_H)**0.5

c     compute new gL, br values using converged lx:

c      newbr = lx / exp(4./(1. - lambda))

c      newgL = newbr / (xshear/shearsol)**(1./6.)
c     1     / (sqrt(1. + (shearsol * rossbysol))
c     2     / (sqrt(1. + (xshear * rossby))))

c      print *,'newgL....................= ',newgL
c      print *,'newbr....................= ',newbr

c     From eq. 14, we find:

      Br_star = Br_sol * solrad**2. / r(1,kk+1)**2.  !eq 28
c     print *,'Br_star_28...............= ',Br_star
      Br_star = br * Br_sol / gL    !eq 14
      Br_star_check = newbr * Br_sol/ newgL

c     turns out the ratio br/gL and newbr/newgL are the same!
      print *,'Br_star..................= ',Br_star
c     print *,'Br_star_check............= ',Br_star_check

c     Eq. 30: The radial Alfven speed

      uA = Br_star / (4.*pi)**(0.5) /
     1        (pnot/(anot**2.))**(0.5)
c     print *,'uA.......................= ',uA

c     uA_check = Br_star_check / (4.*pi)**(0.5) /
c     1        (pnot/(anot**2.))**(0.5)

c     Eq. 40: Time evolution of angular velocity
c     omegadot = (-q_inertia * (omeg(kk+1)
c     1          * (Br_star**2.) * (r(1,kk+1)**2.))
c     2          / (0.059 * xm(kk+1) * uA))
      omegadot = -q_inertia * bosurfomeg
     1          * Br_star**2. * r(1,kk+1)**2.
     2          / (0.059 * xm(kk+1) * uA)

      print *,'omegadot [rad/s^2].......= ',omegadot

c     I am using a different equation for omegadot now, Reiners...2012, eq 6:

      Omega_crit = 8.56d-6  ! s^-1
      ReinersC = 2.66d3     ! g^5/3 cm^-10/3 s

      if (bosurfomeg .ge. Omega_crit) then
         dJdt = -ReinersC*bosurfomeg*r(1,kk+1)**(16./3.)
     1          *xm(kk+1)**(-2./3.)
      else
         dJdt = -ReinersC*bosurfomeg**5.*Omega_crit**(-4.)
     1          *r(1,kk+1)**(16./3.)*xm(kk+1)**(-2./3.)
      endif

      deltaJ = dJdt*dth(2)
      delta_Omega = deltaJ/dmi(kk+1)/r(1,kk+1)**2.
c     print *,'delta_Omega/dth(2).......= ',delta_Omega/dth(2)
c     omegadot = delta_Omega/dth(2)

      rdot = (r(2,kk+1)-r(1,kk+1))/dth(2)
      mdot = peryear*sol/secpy

      omegadot2 = (dJdt-mdot*r(1,kk+1)**2.*bosurfomeg
     1              -2.*dmi(kk+1)*r(1,kk+1)*rdot*bosurfomeg)
     2              /(dmi(kk+1)*r(1,kk+1)**2.)
      print *,'omegadot TWO.............= ',omegadot2

c      omegadot_check = (-q_inertia * (omeg(kk+1)/Omega_sol)
c     1          * (Br_star_check**2.) * ((r(1,kk+1)/solrad)**2.))
c     2          / (0.059 * xm(kk+1) * uA_check)

      bosurfomeg = bosurfomeg + omegadot*dth(2)
      print *,'bosurfomeg...............= ',bosurfomeg

      print *,'omegadot [rad/s^2].......= ',omegadot
      print *,'(angular loss) [deg/yr^2]= ',omegadot
     1           * (180./pi) * (3.15d7)**2.
c      print *,'omegadot_CHECK [rad/s^2] = ',omegadot_check
c      print *,'angloss_CHECK [deg/yr^2] = ',omegadot_check
c     1           * (180./pi) * (3.15d7)**2.

      Mdot_star = -lx * Mdot_sun
c     print *,'m_dot...[M*_star/M*_sun].= ',-lx
      print *,'M_dot_star [g/s].........= ',Mdot_star
      print *,'M_dot_star [sol/yr]......= ',Mdot_star/sol*secpy
      print *,'peryear .................= ',peryear
      print *,'peryear*timestep [sol]...= ',peryear*dth(2)/secpy
      print *,'timestep [10^14 sec].....= ',dth(2)/1.d14
      boperyear = Mdot_star/sol*secpy
      bodomeg = omegadot

      open(59,file='blackowen_vals.txt')
      bo_omega = omeg(kk+1)
      bo_omega2 = omeg(kk)
      vconv = h(kk+1)
      write(59,58) Tnot,lx,omegadot,Br_star,Mdot_star,rossby,
     1             pnot,bo_omega,bo_omega2,tauc,taug,vconv,
     2             bosurfomeg
 58   format(1p13e11.3)
 66   format(1p2e11.3)
      close(59)

      if (doprint .eq. 1) then
         print *,' '
         print *,
     1      'END BLACKOWEN BLACKOWEN BLACKOWEN BLACKOWEN BLACKOWEN'
         print *,' '
      endif

      return

      end
