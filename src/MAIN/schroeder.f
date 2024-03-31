
      subroutine schroeder(te,tlum,rph,totmass,xmdot,vinfin)
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

      real*8 te,tlum,rph,totmass,xmdot,vinfin
      real*8 P,Plog,reimersmdot,Prat,Trat
      real*8 vesc,M,R,g,c,a,b,vratio, gsol
      
      parameter(c = 8.0d-14)
      parameter(a = -1.9777d-1)
      parameter(b = 1.1681d0)
      
      
c---------------------------------------------------------------
c.. section to calculate terminal velocity of wind, using linear
c.. fit to empirical data for 7 AGB stars in Dupree & Reimers
      gsol = (grav * sol / solrad**2.0d0)
      M = totmass*sol
c.. escape velocity
      vesc = (2.0d0*grav*M/rph)**0.5
c.. vinfinity/vesc
      vratio = a*log10(tlum) + b
c.. terminal velocity
      vinfin = vratio * vesc

c--------------------------------------------------------------

      R = rph / solrad 
c.. log(period in days)
      
c.. surface gravity in solar units
      g = totmass /  (R**2)
c.. Dupree & Reimers mass loss
      xmdot = - c * tlum / R / g * (te/4.0d3) * 
     1               (1+ gsol/(4.3d3*grav*M/rph**2.0d0))
c     write(*,'(a20,1p8e12.3)')"SC mass loss ",
c    1  xmdot,te,tlum,gsol,grav*M/rph**2.0d0,R

      return
      end

