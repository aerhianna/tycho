      subroutine vink(te,tlum,rph,totmass,vinkmdot,vinfin,zmetal)
c..   Mass loss routine for evolved stars
c     Uses the mass loss routine of Dupree, A. K. & Reimers,
c     D. (1987, in Exploring the Universe with the IUE Satellite
c     eds. Y Kondo et al. (Dordrecht: Reidel) p. 321) for RGB
c     stars (Teff < 7500 K) Modified for AGB stars (Period >
c     100d) to use mass loss formulation of Bloeker, T. (1995,
c     A&A 297, 727.
c
c	  3/26/15, Modified for O and B stars (with metallicity
c	  ranges between 0.01 and 10 Z/Zsol) from J.S. Vink,
c	  A. de Koter, and H.J.G.L.M. Lamers (2001, A&A 369, 574);
c	  and Modified for O and B stars from the LMC, from M.R. Mokiem
c	  et al. (2007, A&A 465, 1003).
c---------------------------------------------------------------
c Below copied from reimers.f

	  implicit none

	  include 'cconst'

	  real*8 te,tlum,rph,totmass,vinkmdot,vinfin,zmetal
	  real*8 vinkmdot_hot,vinkmdot_cool,Te_jump
	  real*8 vesc,M,R,g,c,a,b,vratio

	  parameter(c = 5.5d-13)
	  parameter(a = -1.9777d-1)
	  parameter(b = 1.1681d0)

c---------------------------------------------------------------
c.. section to calculate terminal velocity of wind, using linear
c.. fit to empirical data for 7 AGB stars in Dupree & Reimers

	  M = totmass*sol
	  
c.. escape velocity
	  vesc = (2.0d0*grav*M/rph)**0.5d0
c.. vinfinity/vesc
	  vratio = a*log10(tlum) + b
c.. terminal velocity
	  vinfin = vratio * vesc

c--------------------------------------------------------------

	  R = rph / solrad

c--------------------------------------------------------------

c.. zmetal is Z/Zsol?

c..for Gamma = 0.206

c	  vinkmdot1 = 10**((0.851)*log(zmetal) - (5.732))

c..for Gamma = 0.130

c	  vinkmdot2 = 10**((0.842)*log(zmetal) - (6.439))

c..for Gamma = 0.434

c	  vinkmdot3 = 10**((0.878)*log(zmetal) - (4.84))

c--------------------------------------------------------------

c.. In the critical temperature range between 22 500 ≤ Teff ≤ 27 500 K,
c	either Eq. (24) or Eq. (25) should be used depending on the position of
c	the bi-stability jump given by Eq. (15) -- below:

	  Te_jump = (61.2d0 + (2.59d0*(0.889d0*log10(zmetal)
     1              - 13.636d0)))*(1.0d3)

c.. Te_jump is given in units of kK, so I'm multiplying by 1000 to get K.

c--------------------------------------------------------------

c.. eq 24, For temperatures 27500 < Te <(eq) 50000:

      vinkmdot_hot = 10.0d0**(-6.697d0
     1              + (2.194d0*log10(tlum/(10.0d0**5.0d0)))
     1              - (1.313d0*log10(totmass/30.0d0))
     1              - (1.226d0*log10(2.6d0/2.0d0))
     1              + (0.933d0*log10(te/4.0d4))
     1              - (10.92d0*((log10(te/4.0d4))**2.0d0))
     1              + (0.85d0*log10(zmetal)))

c.. eq 25, For temperatures 12500 <(eq) Te <(eq) 22500, (here, ratio of vinfin/vesc = 1.3):

      vinkmdot_cool = 10.0d0**(-6.688d0
     1              + (2.210d0*log10(tlum/(10.0d0**5.0d0)))
     1              - (1.339d0*log10(totmass/30.0d0))
     1              - (1.601d0*log10(1.3d0/2.0d0))
     1              + (1.07d0*log10(te/2.0d4))
     1              + (0.85d0*log10(zmetal)))

c.. complete mass-loss recipe in- cluding the metallicity dependence,
c	(from paper) "Mdot is in M⊙ yr−1, L∗ and M∗ are in solar units and Teff is in K.
c	In this range the Galactic ratio of v∞/vesc = vratio = 2.6."

      if(te .ge. 1.25d4 .and. te .le. Te_jump)then
	     vinkmdot = -vinkmdot_cool

      endif

c      if(te .gt. Te_jump .and. te .le. 5.0d4)then
      if(te .gt. Te_jump)then
         vinkmdot = -vinkmdot_hot
         
c     Adjustments for testing hot star mass loss variations
      vinkmdot = 0.33*vinkmdot
         
      endif
c	write(*,*)vinkmdot,te,zmetal,Te_jump,tlum,M
      return
      end



