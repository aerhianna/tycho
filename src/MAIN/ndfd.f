      subroutine ndfd(v1,d,t,p,pt,vpv,e,et,vev,s,si,k)
      implicit none

c..   fowler-hoyle for relativistic nondegenerate electron-positron gas 
c..   v1 = 1/rho
c..   d  = rho * Ye
c..   t  = temperature in Kelvin
c..   p, e, etc in cgs

      include 'cconst'

c..   use modern constants
c..pairc is me*c**2*avagadro in ergs
      real*8    pairc, tpairc, dcons
      parameter( pairc = avagadro * evrestm * ergspermev,
     1     tpairc = pairc/rgas, 
     2     dcons  = 2.0d0*pi4*(egrestm * crad/planck)**3/avagadro )

      integer*4 k
      real*8    v1,d,t,p,pt,vpv,e,et,vev,s,si
      real*8    pgas,ptgas,pvgas,egas,etgas,evgas,sgas
      real*8    b,ppair,ptpair,pvpair,epair,etpair,evpair,spair
      real*8    bkf,bk1,bk2,rhofak,bfact,cfact,afact,yp,ypt,gam
c---------------------------------------------------------------------- 
c..   si is ue/kT
c..   s  is S/avagadro*k

      pgas  =  rgas  * t * d
      ptgas =  rgas  * d
      pvgas = -pgas  / v1
      egas  =  pgas  * v1 * 1.5d0
      etgas =  ptgas * v1 * 1.5d0
      evgas =  0.0d0
c..   sgas is NOT ACCURATE
      sgas  =  0.0d0
c      sgas  =  (pgas/d+egas)/(rgas*t)
cccccccccccccccccccc

      b   = tpairc / t
      if( b .gt. 100.0d0 )then
         ppair  = 0.0d0
         ptpair = 0.0d0
         pvpair = 0.0d0
         epair  = 0.0d0
         etpair = 0.0d0
         evpair = 0.0d0
         spair  = 0.0d0
      else
c..   pair gas component (Chandrasekhar; Fowler and Hoyle)
c..   1.25331 is sqrt(pi/2)
         bkf    = 1.25331d0* sqrt(b) * exp(-b)
         bk2    = bkf * b * (1.0d0 + 1.875d0/b + 0.82031d0/b**2 )
         bk2    = dmin1( bk2, 2.0d0 )
         bk1    = bkf * (1.0d0 + 0.375d0/b )
         bk1    = dmin1( bk1, 1.0d0 )
         rhofak = dcons * bk2/b**3
         bfact  = d / rhofak
         cfact  = sqrt( bfact**2 + 4.0d0)
         afact  = 0.5d0*( bfact + cfact )
c..   afact is e(mu(e)/kt)
c..   yp is Y*rho for pairs * 2 = rho*( Y+ + Y- )
         yp = 2.0d0 * rhofak / afact
c..   J3(0)/[3 J2(0)] = 1.05051
c     yp = 1.05051d0*yp
c..   APPROXIMATION TO DERIVATIVE USED HERE
         ypt    = yp* ( b + 3.0d0 )/t
         ppair  = rgas * yp * t
         ptpair = (ypt*t + yp)*rgas
         pvpair = ppair / v1 * bfact/cfact
         gam    = 3.0d0
         epair  = yp*( gam*rgas*t + pairc ) * v1
         etpair = ypt*(gam*rgas*t + pairc ) + gam*rgas*yp
         etpair = etpair * v1
         evpair = yp*( gam*rgas*t + pairc )
         spair  = (epair/v1 + ppair)/(rgas*d*t)
      endif

      p    = pgas  + ppair
      pt   = ptgas +ptpair
      vpv  = pvgas +pvpair
c..   energy per unit mass
      e    = egas  + epair
      et   = etgas + etpair
      vev  = evgas + evpair
c..   entropy per unit Ye
      s    = sgas + spair

      return
      end









