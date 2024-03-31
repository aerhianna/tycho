      subroutine reimers(te,tlum,rph,totmass,xmdot,vinfin)
c..Mass loss routine for evolved stars
c     Uses the mass loss routine of Dupree, A. K. & Reimers, 
c     D. (1987, in Exploring the Universe with the IUE Satellite
c     eds. Y Kondo et al. (Dordrecht: Reidel) p. 321) for RGB 
c     stars (Teff < 7500 K) Modified for AGB stars (Period >
c     100d) to use mass loss formulation of Bloeker, T. (1995,
c     A&A 297, 727
c---------------------------------------------------------------

      implicit none

      
      include 'cconst'
      real*8 te,tlum,rph,totmass,xmdot,vinfin
      real*8 P,Plog,reimersmdot
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
c.. log(period in days)
      Plog = -1.92d0 - 0.73d0 * log10(totmass) 
     1  + 1.86d0 * log10(R)

      P = 10.0d0**Plog
c.. surface gravity in solar units
      g = totmass / (R**2)
c.. Dupree & Reimers mass loss
      reimersmdot = - c * tlum / R / g
c.. Bloeker pulsational mass loss condition      
c      if(P .gt. 100)then
c        xmdot = 4.83d-9 * totmass**(-2.1) * tlum**2.7 
c     1    * reimersmdot
c      else
        xmdot = reimersmdot
c      endif     
cccccccccccccccccccccccccccc force reimers for massive star

      return
      end



